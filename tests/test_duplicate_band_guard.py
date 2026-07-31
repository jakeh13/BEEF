"""Fast, meep-free unit tests for FDTDSolver._iter_distinct_bands() -- the
duplicate-eigenmode guard (validation_markdown/modal_multiplexer_validation.md
gap 6).

meep is faked out entirely (a minimal stand-in module registered in
sys.modules, not the real package) so this exercises the guard's own retry/
threshold logic in isolation from any real eigensolve, the same "guard the
checker itself" pattern tests/test_mode_identity.py's
TestModeCorrespondenceLogic and tests/test_mode_identity_overlap.py already
use for the underlying overlap math. FDTDSolver() itself never imports meep
at construction time (see test_fdtd_portspec.py's own docstring, confirmed
there directly), so a default-constructed instance is enough to call the
method on.
"""
from __future__ import annotations

import sys
import types

import numpy as np
import pytest

from beef.solvers.fdtd import FDTDSolver


class _FakeVector3:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x, self.y, self.z = x, y, z


class _FakeEigenmode:
    def __init__(self, k_x, ey_fn):
        self.k = _FakeVector3(x=k_x)
        self._ey_fn = ey_fn

    def amplitude(self, point, component):
        if component == "Ey":
            return complex(self._ey_fn(point.y))
        if component == "Ez":
            return 0j
        raise AssertionError(f"unexpected component {component!r}")


class _FakeSim:
    """band_modes: {band_num: _FakeEigenmode}. Missing band_num -> the next
    higher band_num's mode is reused (simulates get_eigenmode always
    returning SOMETHING even past the last distinct mode, matching how a
    real over-request would behave in the exhaustion test)."""

    def __init__(self, band_modes):
        self.band_modes = band_modes
        self.requested_band_nums = []
        self.requested_seeds = []  # kpoint.x / f == the seed n_eff passed in

    def get_eigenmode(self, f, direction, where, band_num, kpoint, parity=None):
        self.requested_band_nums.append(band_num)
        self.requested_seeds.append(kpoint.x / f)
        return self.band_modes[band_num]


@pytest.fixture(autouse=True)
def fake_meep(monkeypatch):
    fake = types.ModuleType("meep")
    fake.Vector3 = _FakeVector3
    fake.X = "X"
    fake.Ey = "Ey"
    fake.Ez = "Ez"
    fake.NO_PARITY = "NO_PARITY"
    monkeypatch.setitem(sys.modules, "meep", fake)
    yield


@pytest.fixture
def solver():
    return FDTDSolver()


def _mode(k_x, sine_order):
    return _FakeEigenmode(k_x, lambda y, n=sine_order: np.sin((n + 1) * np.pi * y))


class TestNoDuplicates:
    def test_yields_n_modes_in_band_order_when_all_distinct(self, solver):
        sim = _FakeSim({1: _mode(1.0, 0), 2: _mode(0.9, 1), 3: _mode(0.8, 2)})
        where = None  # unused by the fake sim
        modes = list(solver._iter_distinct_bands(sim, where, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=3, n_clad=0.0))
        assert len(modes) == 3
        assert [m.k.x for m in modes] == [1.0, 0.9, 0.8]
        assert sim.requested_band_nums == [1, 2, 3]  # no retries needed


class TestOneDuplicate:
    def test_skips_a_duplicate_band_and_uses_the_next_one(self, solver):
        # band_num=2 is an EXACT copy of band_num=1's mode (the failure mode
        # this guard exists for); band_num=3 is genuinely the next distinct
        # mode and must become the 2nd yielded mode instead.
        sim = _FakeSim({
            1: _mode(1.0, 0),
            2: _mode(1.0, 0),   # duplicate of band 1
            3: _mode(0.9, 1),
            4: _mode(0.8, 2),
        })
        modes = list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=3, n_clad=0.0))
        assert [round(m.k.x, 6) for m in modes] == [1.0, 0.9, 0.8]
        assert sim.requested_band_nums == [1, 2, 3, 4]  # band 2 tried and rejected

    def test_logs_a_warning_on_the_skip(self, solver, caplog):
        import logging
        sim = _FakeSim({1: _mode(1.0, 0), 2: _mode(1.0, 0), 3: _mode(0.9, 1)})
        with caplog.at_level(logging.WARNING, logger="beef.solvers.fdtd"):
            list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=2, n_clad=0.0))
        assert any("duplicates an already-accepted eigenmode" in r.message for r in caplog.records)


class TestExhaustion:
    def test_raises_if_every_extra_attempt_is_also_a_duplicate(self, solver):
        # Every band beyond the first is a duplicate of band 1 -- the guard
        # must give up cleanly (not loop forever) once its extra-attempt
        # budget (_DUP_BAND_MAX_EXTRA) is exhausted.
        modes = {1: _mode(1.0, 0)}
        for b in range(2, 30):
            modes[b] = _mode(1.0, 0)  # all duplicates of band 1
        sim = _FakeSim(modes)
        with pytest.raises(RuntimeError, match="could not find .* distinct eigenmodes"):
            list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=2, n_clad=0.0))


class TestThreshold:
    def test_similar_but_genuinely_different_modes_are_not_treated_as_duplicates(self, solver):
        """Regression guard for the threshold itself: two DIFFERENT sine
        orders must never trip the duplicate check, however close their
        eigenvalue (k_x) happens to be -- only field SHAPE decides this,
        not neff proximity (that's the entire point of this check existing
        instead of an n_eff-based one)."""
        sim = _FakeSim({1: _mode(1.0, 0), 2: _mode(0.999999, 1)})  # near-identical k_x, different shape
        modes = list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=2, n_clad=0.0))
        assert len(modes) == 2
        assert sim.requested_band_nums == [1, 2]


class TestGuidedEdgeStop:
    """fdtd_vs_EME_mode_validation.md item 14: a band_num with no solution at
    the target frequency does not raise -- MPB's k-search limit-cycles
    forever instead of erroring (measured: 161 solves on one band before
    SLURM killed the job). The only defence is to never ask for that
    band_num, which means stopping as soon as a raw solve's n_eff crosses the
    cladding index -- nothing above it is guided, and asking further is what
    hangs."""

    def test_raises_once_a_raw_solve_crosses_the_cladding_index(self, solver):
        # Bands 1-2 guided (n_eff 1.0, 0.9); band 3 is unguided (0.5, below
        # the n_clad=0.8 threshold) -- requesting n_modes=3 must fail loudly
        # rather than silently accept the unguided band as the 3rd channel.
        sim = _FakeSim({1: _mode(1.0, 0), 2: _mode(0.9, 1), 3: _mode(0.5, 2)})
        with pytest.raises(RuntimeError, match="supplies only 2 guided mode"):
            list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=3, n_clad=0.8))
        # Only 3 attempts -- the walk stopped at the unguided band rather
        # than continuing toward max_band_num (where a real cross-section's
        # non-convergent k-search would live).
        assert sim.requested_band_nums == [1, 2, 3]

    def test_does_not_fire_when_every_requested_mode_is_guided(self, solver):
        # Same fake as TestNoDuplicates, but now with a real cladding value
        # below all three fake n_eff -- confirms the new check is inert when
        # nothing crosses it, not just when n_clad is an extreme override.
        sim = _FakeSim({1: _mode(1.0, 0), 2: _mode(0.9, 1), 3: _mode(0.8, 2)})
        modes = list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=3, n_clad=0.79))
        assert len(modes) == 3

    def test_cladding_check_runs_before_the_duplicate_overlap_sampling(self, solver):
        # band_num=2 has the SAME field shape as band 1 (would be flagged a
        # duplicate) AND a lower k_x that puts it at/below n_clad -- the
        # cladding check must fire first (it's the physically correct
        # reading: an unguided "duplicate" isn't a real second look at band
        # 1, it's just noise below cutoff) rather than being silently
        # skipped as a duplicate and letting the walk continue to band 3,
        # which doesn't exist in this fake at all.
        sim = _FakeSim({1: _mode(1.0, 0), 2: _mode(0.9, 0)})
        with pytest.raises(RuntimeError, match="supplies only 1 guided mode"):
            list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=2, n_clad=0.95))


class TestSeedLaddering:
    """fdtd_vs_EME_mode_validation.md item 12: seeding every attempt from the
    same fixed guess is the textbook shape of the duplicate this guard works
    around (two solves launched from an identical point can converge to the
    same eigenpair). Seeding attempt k from attempt (k-1)'s own raw converged
    n_eff instead was measured (diagnose_kpoint_ladder_v2.py, job 11588309)
    to recover a genuinely distinct mode the fixed seed silently drops, at no
    extra cost. These tests pin the MECHANISM (what seed gets passed) --  the
    physical claim (that laddering recovers more real modes) is what the real
    meep/MPB job already verified."""

    def test_first_attempt_seeds_from_the_fixed_default(self, solver):
        from beef.solvers.fdtd import _KPOINT_SEED_NEFF
        sim = _FakeSim({1: _mode(1.0, 0)})
        list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=1, n_clad=0.0))
        assert sim.requested_seeds == [_KPOINT_SEED_NEFF]

    def test_later_attempts_seed_from_the_previous_raw_neff_not_the_fixed_default(self, solver):
        from beef.solvers.fdtd import _KPOINT_SEED_NEFF
        sim = _FakeSim({1: _mode(1.0, 0), 2: _mode(0.9, 1), 3: _mode(0.8, 2)})
        list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=3, n_clad=0.0))
        assert sim.requested_seeds == [_KPOINT_SEED_NEFF, 1.0, 0.9]

    def test_a_skipped_bands_raw_neff_still_seeds_the_next_attempt(self, solver):
        # band 2 has the SAME field shape as band 1 (sine_order=0, so it gets
        # skipped as a duplicate) but a DIFFERENT k_x (0.95, not 1.0) -- this
        # is what actually distinguishes policy B (ladder off the raw
        # previous solve, what shipped) from policy D (ladder off the last
        # ACCEPTED band): if band 3 were seeded from band 1's accepted 1.0,
        # this test would fail. It must be seeded from band 2's own raw
        # 0.95, even though band 2 itself was rejected.
        from beef.solvers.fdtd import _KPOINT_SEED_NEFF
        sim = _FakeSim({1: _mode(1.0, 0), 2: _mode(0.95, 0), 3: _mode(0.9, 1)})
        list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=2, n_clad=0.0))
        assert sim.requested_seeds == [_KPOINT_SEED_NEFF, 1.0, 0.95]
