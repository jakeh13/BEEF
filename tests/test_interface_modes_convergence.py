"""Tests for beef.interface_modes_convergence.interface_modes_convergence()
-- fdtd_vs_EME_mode_validation.md item 6 (layer 2, completeness: does the
final S-matrix plateau as interface_modes grows).

Pure numpy/SMatrix: operates on already-built SMatrix objects, no solver or
geometry involved, so this needs neither meep nor meow. Structured as the
direct sibling of test_cut_position_invariance.py (item 7's own tests) --
same validation contract, same "sorted regardless of input order" guarantee
-- but the convergence quantity is a STEP toward a plateau, not a fitted
slope toward zero, since a mode-count sweep has no reason to move linearly.
"""
from __future__ import annotations

import numpy as np
import pytest

from beef import SMatrix
from beef.interface_modes_convergence import interface_modes_convergence

WAVELENGTHS = np.array([1.55])


def _smat(coeff: complex, port: str = "out", in_mode: int = 1, out_mode: int = 1) -> SMatrix:
    s = SMatrix(wavelengths=WAVELENGTHS)
    s.set(np.array([coeff]), port=port, in_mode=in_mode, out_mode=out_mode)
    return s


class TestBasicConvergence:
    def test_already_plateaued_has_zero_step(self):
        counts = [8, 10, 12, 14]
        # Same |coeff|^2 = 0.25 at every mode count -- adding modes changes
        # nothing, the basis was already complete.
        smats = [_smat(0.5 + 0j) for _ in counts]
        result = interface_modes_convergence(smats, counts)
        key = ("out", 1, 1)
        assert result["keys"][key]["last_step"] == pytest.approx(0.0, abs=1e-9)
        assert result["keys"][key]["max_step"] == pytest.approx(0.0, abs=1e-9)
        assert result["worst_key"] == key
        assert result["worst_last_step"] == pytest.approx(0.0, abs=1e-9)

    def test_still_converging_shows_shrinking_steps(self):
        counts = [4, 8, 12, 16]
        # transmission approaches 0.5 with shrinking increments -- the shape
        # a real basis-completeness sweep should have (big early jump as the
        # dominant mode is added, then diminishing returns).
        transmissions = np.array([0.30, 0.46, 0.49, 0.495])
        smats = [_smat(np.sqrt(t) + 0j) for t in transmissions]
        result = interface_modes_convergence(smats, counts)
        key = ("out", 1, 1)
        steps = result["keys"][key]["steps"]
        np.testing.assert_allclose(steps, np.diff(transmissions))
        # Steps shrink -- each is smaller in magnitude than the last.
        assert abs(steps[-1]) < abs(steps[0])
        assert result["keys"][key]["last_step"] == pytest.approx(steps[-1])
        assert result["keys"][key]["max_step"] == pytest.approx(abs(steps[0]))

    def test_values_and_mode_counts_are_returned_sorted_regardless_of_input_order(self):
        counts = [12, 8, 16]   # deliberately out of order
        transmissions = {8: 0.30, 12: 0.46, 16: 0.49}
        smats = [_smat(np.sqrt(transmissions[c]) + 0j) for c in counts]
        result = interface_modes_convergence(smats, counts)
        key = ("out", 1, 1)
        assert list(result["keys"][key]["mode_counts"]) == [8, 12, 16]
        assert result["keys"][key]["values"] == pytest.approx([0.30, 0.46, 0.49])


class TestMultipleKeys:
    def test_worst_key_picks_the_largest_magnitude_last_step(self):
        counts = [8, 16]
        smats = [SMatrix(wavelengths=WAVELENGTHS) for _ in counts]
        # "bar": 0.50 -> 0.50 (plateaued). "cross": 0.30 -> 0.45 (still moving
        # a lot between the only two counts tried).
        for i, c in enumerate(counts):
            smats[i].set(np.array([np.sqrt(0.5) + 0j]), port="bar", in_mode=1, out_mode=1)
            cross_t = 0.30 if c == 8 else 0.45
            smats[i].set(np.array([np.sqrt(cross_t) + 0j]), port="cross", in_mode=1, out_mode=1)
        result = interface_modes_convergence(smats, counts)
        assert result["worst_key"] == ("cross", 1, 1)
        assert result["worst_last_step"] == pytest.approx(0.15, abs=1e-9)
        assert result["keys"][("bar", 1, 1)]["last_step"] == pytest.approx(0.0, abs=1e-9)


class TestValidation:
    def test_raises_on_mismatched_lengths(self):
        with pytest.raises(ValueError, match="same length"):
            interface_modes_convergence([_smat(1.0), _smat(1.0)], [8])

    def test_raises_on_fewer_than_two_points(self):
        with pytest.raises(ValueError, match="at least 2"):
            interface_modes_convergence([_smat(1.0)], [8])

    def test_raises_on_mismatched_key_sets(self):
        a = _smat(1.0, port="out", in_mode=1, out_mode=1)
        b = SMatrix(wavelengths=WAVELENGTHS)
        b.set(np.array([1.0]), port="out", in_mode=1, out_mode=2)  # different key
        with pytest.raises(ValueError, match="same .* keys"):
            interface_modes_convergence([a, b], [8, 12])
