"""The ring resonator's copy of the monitor-grid bug (dc_validation.md item 10).

WHAT WAS WRONG
--------------
`ring_resonator_hybrid.py` and `ring_resonator_hybrid_cascaded.py` build their
meep simulations by hand rather than through `FDTDSolver` -- the two ring ports
sit on y-normal planes and nothing in `beef/solvers/fdtd.py` draws those. So
each script carried its own copy of the defect item 1 fixed:

    fwidth = (1.0 / WL_MIN - 1.0 / WL_MAX) * 0.5
    sim.add_mode_monitor(fcen, fwidth, N_WL_FDTD, ...)

Character for character the same line, in raw meep, which is exactly why the
solver-side fix never reached it. meep spreads `nfreq` points evenly over
`[fcen - df/2, fcen + df/2]` **in frequency**, so with that `df` a request for
1.50-1.60 um actually sampled 1.5254-1.5754 um: half the requested span, and
ordered ascending in frequency, i.e. **descending in wavelength**. Indexing
that against an ascending-lambda array mislabels every point but the centre.

WHY IT SURVIVED ITEMS 1-3
-------------------------
The cache. These scripts fingerprint their own runs through a hand-rolled
`_sim_params()` dict, so nothing in items 1-3 changed their hash and all six
of their FDTD cache entries stayed *reachable* -- a re-run would have loaded
pre-fix data and reported it as a fresh result. They were deliberately left on
disk rather than deleted, since deleting only forces re-simulation of equally
wrong numbers.

`_sim_params()` now emits `resolved_monitor_freqs` / `resolved_source_fwidth`,
mirroring `FDTDSolver._resolved_spectral_config()`. Both are derived from
`wavelengths`, which the config already hashed, so they change nothing about
correctness today -- their job is to put the *derivation* in the fingerprint,
which is what strands the old entries. Verified against the real cache: all 6
FDTD entries now miss, and all 39 `solve_neff` / `solve_body_supermodes`
entries still hit (those are MPB mode solves, never touched by this bug, and
stranding them would force pointless re-simulation).

The tests below are source-level because importing either script *runs* it --
they are top-level scripts, not modules, and a bare import would kick off a
multi-hour meep job.
"""
from __future__ import annotations

import pathlib
import re

import numpy as np
import pytest

from beef.solvers.fdtd import monitor_freq_grid, source_fwidth

_SCRIPTS = pathlib.Path(__file__).parent.parent / "examples" / "ring_resonator" / "python"
# ring_resonator_hybrid_cascaded.py was guarded here too, until it moved to
# github.com/jakeh13/BEEF_archive on 2026-09-18. Frozen there, so the defect
# this class guards cannot be newly introduced into it; reaching across to a
# sibling checkout would break this suite on a clean clone of BEEF alone.
_RING_SCRIPTS = ("ring_resonator_hybrid.py",)

# What both scripts declare, and what run_*_fdtd() sweeps.
WL_MIN, WL_MAX = 1.50, 1.60
N_WL = 101

# The band the pre-fix triple actually sampled, from
# fcen = 1/1.55, df = (1/1.50 - 1/1.60) * 0.5.
PREFIX_BAND = (1.52537, 1.57544)


def _source(name: str) -> str:
    return (_SCRIPTS / name).read_text()


def _body(name: str) -> str:
    """Source with docstrings/comments stripped, so a prose mention of the old
    line (there is one, describing the defect) cannot pass for live code."""
    src = _source(name)
    src = re.sub(r'"""(?:.|\n)*?"""', "", src)
    return "\n".join(ln.split("#", 1)[0] for ln in src.splitlines())


class TestScriptsNoLongerCarryTheDefect:
    """Fast, no meep. These are the regression guards -- the physics tests
    below cost hours and nobody runs them casually."""

    @pytest.mark.parametrize("name", _RING_SCRIPTS)
    def test_no_fcen_df_nfreq_triple(self, name):
        """The monitor call must take an explicit frequency list."""
        body = _body(name)
        assert "fwidth, N_WL_FDTD" not in body, (
            f"{name} still passes meep's (fcen, df, nfreq) triple to "
            f"add_mode_monitor(). That triple cannot reproduce a wavelength "
            f"array linear in lambda and orders points descending in "
            f"wavelength -- use monitor_freq_grid()."
        )

    @pytest.mark.parametrize("name", _RING_SCRIPTS)
    def test_no_half_span_fwidth(self, name):
        body = _body(name)
        assert not re.search(r"1\.0\s*/\s*WL_MIN\s*-\s*1\.0\s*/\s*WL_MAX\s*\)\s*\*\s*0\.5", body), (
            f"{name} still computes fwidth as half the requested span. Even "
            f"once the monitors are fixed this starves the band edges -- use "
            f"source_fwidth(wavelengths, SOURCE_BANDWIDTH_FACTOR)."
        )

    @pytest.mark.parametrize("name", _RING_SCRIPTS)
    def test_uses_the_library_helpers(self, name):
        body = _body(name)
        for helper in ("monitor_freq_grid", "source_fwidth"):
            assert helper in body, f"{name} does not use {helper}()"

    @pytest.mark.parametrize("name", _RING_SCRIPTS)
    def test_results_are_reordered_onto_the_caller_axis(self, name):
        """monitor_freq_grid() is only half the fix. Its frequencies are
        ascending, so a result read off the monitor is still in descending-
        wavelength order until [wl_order] puts it back."""
        body = _body(name)
        assert "wl_order" in body, (
            f"{name} calls monitor_freq_grid() but never applies wl_order. "
            f"The monitor result is still in monitor order -- reversed "
            f"relative to the wavelength array the script reports against."
        )

    @pytest.mark.parametrize("name", _RING_SCRIPTS)
    def test_sim_params_fingerprints_the_resolved_grid(self, name):
        """Without this the six pre-fix cache entries stay reachable and a
        re-run silently serves them."""
        body = _body(name)
        assert "def _sim_params(wavelengths)" in body, (
            f"{name}'s _sim_params() takes no wavelengths, so it cannot "
            f"fingerprint the resolved monitor grid."
        )
        for key in ("resolved_monitor_freqs", "resolved_source_fwidth"):
            assert key in body, f"{name}'s _sim_params() omits {key}"


class TestTheGridItselfIsRight:
    """Numeric, no meep -- pins what the two grids actually are, so the prose
    above can be checked rather than believed."""

    def test_prefix_triple_really_did_sample_a_reversed_half_band(self):
        """Reconstructs the old behaviour from meep's documented rule, to show
        the band in PREFIX_BAND is derived rather than asserted."""
        fcen = 1.0 / 1.55
        df = (1.0 / WL_MIN - 1.0 / WL_MAX) * 0.5
        freqs = np.linspace(fcen - df / 2, fcen + df / 2, N_WL)   # meep's rule
        wl = 1.0 / freqs

        assert wl[0] > wl[-1], "meep's triple orders ascending in frequency"
        assert wl.max() == pytest.approx(PREFIX_BAND[1], abs=1e-4)
        assert wl.min() == pytest.approx(PREFIX_BAND[0], abs=1e-4)
        # Half the requested span, and up to 75 nm off at the edges.
        assert (wl.max() - wl.min()) == pytest.approx((WL_MAX - WL_MIN) / 2, abs=2e-3)
        assert abs(wl[0] - WL_MIN) > 0.02 and abs(wl[-1] - WL_MAX) > 0.02

    def test_new_grid_hits_the_requested_wavelengths_exactly(self):
        wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
        mon_freqs, wl_order = monitor_freq_grid(wavelengths)

        assert np.all(np.diff(mon_freqs) > 0), "meep requires ascending frequencies"
        np.testing.assert_allclose(1.0 / mon_freqs[wl_order], wavelengths, rtol=0, atol=1e-12)

    def test_wl_order_undoes_the_reversal(self):
        """The ascending-lambda case is exactly a flip. Stated explicitly so a
        refactor cannot quietly drop the reorder and still pass."""
        wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
        _, wl_order = monitor_freq_grid(wavelengths)
        np.testing.assert_array_equal(wl_order, np.arange(N_WL)[::-1])

        monitor_ordered = np.arange(N_WL, dtype=float)     # stand-in result
        np.testing.assert_array_equal(monitor_ordered[wl_order], monitor_ordered[::-1])

    def test_source_now_covers_the_full_band(self):
        """The pre-fix fwidth looked adequate only because the monitors shared
        it and sampled the same too-narrow band. Widening the monitors without
        widening the source would have starved the new edges."""
        wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
        full_span = 1.0 / WL_MIN - 1.0 / WL_MAX
        prefix_fwidth = full_span * 0.5

        assert source_fwidth(wavelengths, 2.0) == pytest.approx(2.0 * full_span)
        # Exactly 4x the old pulse: 2x for covering the full span rather than
        # half of it, and 2x again for SOURCE_BANDWIDTH_FACTOR putting the band
        # edges at 0.27 of peak instead of 0.005.
        assert source_fwidth(wavelengths, 2.0) == pytest.approx(4.0 * prefix_fwidth)
