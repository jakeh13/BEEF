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

    def get_eigenmode(self, f, direction, where, band_num, kpoint, parity=None):
        self.requested_band_nums.append(band_num)
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
        modes = list(solver._iter_distinct_bands(sim, where, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=3))
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
        modes = list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=3))
        assert [round(m.k.x, 6) for m in modes] == [1.0, 0.9, 0.8]
        assert sim.requested_band_nums == [1, 2, 3, 4]  # band 2 tried and rejected

    def test_logs_a_warning_on_the_skip(self, solver, caplog):
        import logging
        sim = _FakeSim({1: _mode(1.0, 0), 2: _mode(1.0, 0), 3: _mode(0.9, 1)})
        with caplog.at_level(logging.WARNING, logger="beef.solvers.fdtd"):
            list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=2))
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
            list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=2))


class TestThreshold:
    def test_similar_but_genuinely_different_modes_are_not_treated_as_duplicates(self, solver):
        """Regression guard for the threshold itself: two DIFFERENT sine
        orders must never trip the duplicate check, however close their
        eigenvalue (k_x) happens to be -- only field SHAPE decides this,
        not neff proximity (that's the entire point of this check existing
        instead of an n_eff-based one)."""
        sim = _FakeSim({1: _mode(1.0, 0), 2: _mode(0.999999, 1)})  # near-identical k_x, different shape
        modes = list(solver._iter_distinct_bands(sim, None, f=1.0, mea_c=0.5, mea_w=1.0, n_modes=2))
        assert len(modes) == 2
        assert sim.requested_band_nums == [1, 2]
