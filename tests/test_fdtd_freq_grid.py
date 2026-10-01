"""Regression tests for FDTDSolver's DFT-monitor frequency grid and the
Gaussian source bandwidth that has to cover it.

The bug these lock down: monitors used to be built from meep's
``(fcen, df, nfreq)`` triple with ``df = (1/lam_min - 1/lam_max) * 0.5``.
Meep spreads nfreq points evenly over ``[fcen - df/2, fcen + df/2]`` *in
frequency*, ascending — so the result covered half the requested span, ran
backwards in wavelength, and was then indexed against the caller's
ascending-lambda array. A 1.50-1.60 um 11-point sweep actually sampled
1.5254-1.5754 um, up to 75 nm off, with only the centre point correct.
See dc_validation.md.
"""
from __future__ import annotations

import numpy as np
import pytest

from beef.solvers.fdtd import (
    FDTDSolver, check_source_coverage, monitor_freq_grid, source_fwidth,
)

WAVELENGTHS = np.linspace(1.50, 1.60, 11)   # the standard C-band example sweep
FCEN = 1.0 / float(np.mean(WAVELENGTHS))
SPAN = 1.0 / WAVELENGTHS.min() - 1.0 / WAVELENGTHS.max()


def test_monitor_freqs_are_ascending_and_exact():
    """meep requires ascending frequencies; they must be exactly 1/lambda —
    no (fcen, df, nfreq) triple can reproduce a grid linear in lambda."""
    mon_freqs, _ = monitor_freq_grid(WAVELENGTHS)

    assert np.all(np.diff(mon_freqs) > 0), "meep needs ascending frequencies"
    np.testing.assert_allclose(np.sort(1.0 / mon_freqs), np.sort(WAVELENGTHS))
    # The old triple's band, for contrast: it spanned only 1.5254-1.5754 um.
    assert 1.0 / mon_freqs.max() == pytest.approx(WAVELENGTHS.min())
    assert 1.0 / mon_freqs.min() == pytest.approx(WAVELENGTHS.max())


def test_wl_order_round_trips_monitor_results():
    """A monitor-ordered result reordered by wl_order must land back on the
    caller's wavelengths — this is the indexing the bug got wrong."""
    mon_freqs, wl_order = monitor_freq_grid(WAVELENGTHS)

    # Stand in for alpha[]: tag each monitor slot with the wavelength meep
    # actually computed it at, then reorder.
    monitor_ordered = 1.0 / mon_freqs
    np.testing.assert_allclose(monitor_ordered[wl_order], WAVELENGTHS)


def test_ascending_wavelengths_are_reversed_on_unpack():
    """For the usual ascending-lambda input the reorder is a plain flip —
    stated explicitly so a future refactor can't quietly drop it."""
    _, wl_order = monitor_freq_grid(WAVELENGTHS)
    np.testing.assert_array_equal(wl_order, np.arange(len(WAVELENGTHS))[::-1])


def test_unsorted_wavelengths_round_trip():
    """The caller's order is honoured whatever it is, sorted or not."""
    wl = np.array([1.55, 1.50, 1.60, 1.53])
    mon_freqs, wl_order = monitor_freq_grid(wl)

    assert np.all(np.diff(mon_freqs) > 0)
    np.testing.assert_allclose((1.0 / mon_freqs)[wl_order], wl)


def test_meep_returns_exactly_the_requested_wavelengths():
    """End-to-end against meep itself: the frequencies the monitor reports
    back must be the ones asked for. No FDTD run needed."""
    mp = pytest.importorskip("meep")

    mon_freqs, wl_order = monitor_freq_grid(WAVELENGTHS)
    fcen = 1.0 / float(np.mean(WAVELENGTHS))
    sim = mp.Simulation(
        cell_size=mp.Vector3(2, 2), resolution=10,
        # meep allocates field components from the sources; a monitor cannot
        # be added to a fieldless simulation.
        sources=[mp.Source(mp.GaussianSource(fcen, fwidth=0.1), mp.Ez,
                           center=mp.Vector3())],
    )
    mon = sim.add_mode_monitor(
        mon_freqs, mp.FluxRegion(center=mp.Vector3(0.5, 0), size=mp.Vector3(0, 1))
    )

    got = 1.0 / np.array(mp.get_flux_freqs(mon))
    np.testing.assert_allclose(got[wl_order], WAVELENGTHS, atol=1e-12)


# --------------------------------------------------------------------------
# Source bandwidth
# --------------------------------------------------------------------------

def test_fwidth_scales_the_full_span():
    """The factor multiplies the whole span, not half of it — the old
    `span * 0.5` is exactly the value that starved the band."""
    assert source_fwidth(WAVELENGTHS, 2.0) == pytest.approx(2.0 * SPAN)
    assert source_fwidth(WAVELENGTHS, 1.0) == pytest.approx(SPAN)


def test_single_wavelength_gets_a_real_bandwidth():
    """A one-wavelength request has zero span; fwidth must not collapse to 0
    (meep's Gaussian width is 1/fwidth, so 0 is not a usable pulse)."""
    fwidth = source_fwidth(np.array([1.55]), 2.0)
    assert fwidth > 0
    assert fwidth == pytest.approx(0.2 / 1.55)   # meep's usual 0.2 * fcen


def test_default_factor_clears_the_coverage_floor():
    """Every requested wavelength must carry >=10% of peak amplitude at the
    shipped default, with margin."""
    pytest.importorskip("meep")

    frac = check_source_coverage(
        FCEN, source_fwidth(WAVELENGTHS, FDTDSolver().source_bandwidth_factor),
        WAVELENGTHS,
    )
    assert frac.min() > 0.25
    assert frac.argmax() == len(WAVELENGTHS) // 2   # peak sits at fcen


def test_old_bandwidth_would_be_flagged_as_starved(caplog):
    """Lock in the regression: the pre-fix `span * 0.5` leaves the band edges
    ~5 orders down, and that must now warn rather than pass silently."""
    pytest.importorskip("meep")

    with caplog.at_level("WARNING"):
        frac = check_source_coverage(FCEN, 0.5 * SPAN, WAVELENGTHS)

    assert frac.min() < 1e-6
    assert "starved" in caplog.text


def test_solver_exposes_the_knob():
    """It is a real parameter, not a module constant: physics knobs must not
    live in defaults nothing can see."""
    assert FDTDSolver().source_bandwidth_factor == 2.0
    assert FDTDSolver(source_bandwidth_factor=3.0).source_bandwidth_factor == 3.0
