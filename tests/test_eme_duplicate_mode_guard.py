"""Fast, meow-free unit tests for EMESolver.compute_modes()'s duplicate-mode
guard (validation_markdown/fdtd_vs_EME_mode_validation.md item 0).

meow is never imported: `_solve_modes_raw()` is monkeypatched to return
synthetic mode objects, so this exercises the guard's own detection/skip logic
in isolation from any real eigensolve -- the same "guard the checker" pattern
tests/test_duplicate_band_guard.py uses for the meep-side equivalent.

The bug being guarded is real and measured, not hypothetical: at the modal
multiplexer's interface 0, 16 modes at res 40, meow's own modes 10 and 11 came
back as the identical eigenvector (self-overlap 1.0000, row completeness
2.0003). It is the same failure class FDTDSolver._iter_distinct_bands() was
built for, on the other engine.
"""
from __future__ import annotations

import logging
from types import SimpleNamespace

import numpy as np
import pytest

from beef import EMESolver, SOI_160

NY, NZ = 64, 1


def _mode(order: int, neff: float = 2.0):
    """Synthetic meow-Mode-shaped object. stack_meow_mode_fields() reads
    Ex..Hz off it and remaps into the device frame; only the transverse E
    components are actually compared, so the rest are zeros of matching shape.

    Field is sin((order+1) pi u) on u in [0, 1) -- distinct orders are
    orthogonal, equal orders identical, which is exactly the distinction the
    guard must make.
    """
    u = (np.arange(NY) + 0.5) / NY
    prof = np.sin((order + 1) * np.pi * u)[:, None].astype(complex)
    zero = np.zeros((NY, NZ), dtype=complex)
    return SimpleNamespace(
        neff=complex(neff),
        Ex=prof, Ey=zero.copy(), Ez=zero.copy(),   # meow Ex -> device Ey (compared)
        Hx=zero.copy(), Hy=zero.copy(), Hz=zero.copy(),
    )


def _solver(**kw):
    return EMESolver(num_modes=[4], num_cells=[1], lengths=[1.0],
                     material_stack=SOI_160, **kw)


def _fake_raw(orders):
    """Returns a _solve_modes_raw stand-in yielding `orders` in sequence."""
    def raw(self, cross_section, num_modes):
        return [_mode(o, neff=2.0 - 0.01 * i) for i, o in enumerate(orders[:num_modes])]
    return raw


class TestNoDuplicates:
    def test_returns_the_first_n_when_all_distinct(self, monkeypatch):
        s = _solver()
        monkeypatch.setattr(EMESolver, "_solve_modes_raw", _fake_raw([0, 1, 2, 3, 4, 5]))
        modes = s.compute_modes(cross_section=None, num_modes=4)
        assert len(modes) == 4


class TestDuplicateIsSkipped:
    def test_skips_the_duplicate_and_reaches_further_up(self, monkeypatch):
        """The measured failure: index 2 repeats index 0. The guard must drop
        it and pull in the next distinct mode, so the basis still carries 4
        INDEPENDENT modes rather than 3 plus a copy."""
        s = _solver()
        monkeypatch.setattr(EMESolver, "_solve_modes_raw",
                            _fake_raw([0, 1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
        modes = s.compute_modes(cross_section=None, num_modes=4)
        assert len(modes) == 4
        # Pairwise-distinct check on what actually came back.
        from beef import field_overlap_matrix, stack_meow_mode_fields
        eta = field_overlap_matrix(stack_meow_mode_fields(modes),
                                   stack_meow_mode_fields(modes), dy=1.0, dz=1.0)
        off = eta - np.diag(np.diag(eta))
        assert float(np.max(np.abs(off))) < 0.9

    def test_logs_a_warning_naming_the_duplicated_index(self, monkeypatch, caplog):
        s = _solver()
        monkeypatch.setattr(EMESolver, "_solve_modes_raw",
                            _fake_raw([0, 1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
        with caplog.at_level(logging.WARNING, logger="beef.solvers.eme"):
            s.compute_modes(cross_section=None, num_modes=4)
        assert any("duplicates already-accepted mode" in r.message for r in caplog.records)


class TestNearDegenerateModesSurvive:
    def test_distinct_modes_with_equal_neff_are_kept(self, monkeypatch):
        """The critical negative case. A symmetric cross-section has genuinely
        near-degenerate PAIRS -- measured at 0.10-0.19 mutual overlap on the
        modal multiplexer -- and those are real, independent modes. Identical
        n_eff must not cause a drop: only field SHAPE decides, and the
        threshold sits well above that band."""
        s = _solver()

        def raw(self, cross_section, num_modes):
            return [_mode(o, neff=2.0) for o in range(num_modes)]  # all same neff

        monkeypatch.setattr(EMESolver, "_solve_modes_raw", raw)
        assert len(s.compute_modes(cross_section=None, num_modes=4)) == 4


class TestExhaustion:
    def test_raises_when_every_extra_mode_is_also_a_duplicate(self, monkeypatch):
        s = _solver()
        monkeypatch.setattr(EMESolver, "_solve_modes_raw", _fake_raw([0] * 40))
        with pytest.raises(RuntimeError, match="could not find 4 distinct modes"):
            s.compute_modes(cross_section=None, num_modes=4)


class TestGuardIsSwitchable:
    def test_disabling_it_restores_the_unguarded_behaviour(self, monkeypatch):
        """Off, the duplicate comes straight through -- which is what the
        before/after comparison on a real device needs, and what every result
        recorded before 2026-07-30 was computed with."""
        s = _solver(guard_duplicate_modes=False)
        monkeypatch.setattr(EMESolver, "_solve_modes_raw", _fake_raw([0, 1, 0, 2]))
        modes = s.compute_modes(cross_section=None, num_modes=4)
        from beef import field_overlap_matrix, stack_meow_mode_fields
        eta = field_overlap_matrix(stack_meow_mode_fields(modes),
                                   stack_meow_mode_fields(modes), dy=1.0, dz=1.0)
        assert eta[0, 2] > 0.99  # the duplicate is still there, unguarded

    def test_disabled_guard_does_not_over_request(self, monkeypatch):
        """Off, it must issue exactly the caller's num_modes -- no extra solve
        cost for anyone who opts out."""
        seen = {}

        def raw(self, cross_section, num_modes):
            seen["n"] = num_modes
            return [_mode(o) for o in range(num_modes)]

        monkeypatch.setattr(EMESolver, "_solve_modes_raw", raw)
        _solver(guard_duplicate_modes=False).compute_modes(cross_section=None, num_modes=4)
        assert seen["n"] == 4
