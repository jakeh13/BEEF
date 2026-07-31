"""Tests for beef.cut_position_invariance.cut_position_invariance() --
fdtd_vs_EME_mode_validation.md item 7 (generalizing
directional_coupler_hybrid_stub_sweep.py's one-off 50:50-split sweep into a
device-agnostic check).

Pure numpy/SMatrix: operates on already-built SMatrix objects, no solver or
geometry involved, so this needs neither meep nor meow.
"""
from __future__ import annotations

import numpy as np
import pytest

from beef import SMatrix
from beef.cut_position_invariance import cut_position_invariance

WAVELENGTHS = np.array([1.55])


def _smat(coeff: complex, port: str = "out", in_mode: int = 1, out_mode: int = 1) -> SMatrix:
    s = SMatrix(wavelengths=WAVELENGTHS)
    s.set(np.array([coeff]), port=port, in_mode=in_mode, out_mode=out_mode)
    return s


class TestBasicDrift:
    def test_constant_transmission_has_zero_slope_and_drift(self):
        stubs = [0.0, 0.25, 0.5, 0.75]
        # Same |coeff|^2 = 0.25 at every stub -- a perfectly converged sweep.
        smats = [_smat(0.5 + 0j) for _ in stubs]
        result = cut_position_invariance(smats, stubs)
        key = ("out", 1, 1)
        assert result["keys"][key]["slope"] == pytest.approx(0.0, abs=1e-9)
        assert result["keys"][key]["max_drift"] == pytest.approx(0.0, abs=1e-9)
        assert result["worst_key"] == key
        assert result["worst_slope"] == pytest.approx(0.0, abs=1e-9)

    def test_linearly_drifting_transmission_recovers_the_known_slope(self):
        stubs = np.array([0.0, 0.5, 1.0, 1.5])
        # transmission = 0.2 + 0.1*stub, exactly -- coeff = sqrt(that).
        transmissions = 0.2 + 0.1 * stubs
        smats = [_smat(np.sqrt(t) + 0j) for t in transmissions]
        result = cut_position_invariance(smats, stubs)
        key = ("out", 1, 1)
        assert result["keys"][key]["slope"] == pytest.approx(0.1, abs=1e-9)
        assert result["keys"][key]["intercept"] == pytest.approx(0.2, abs=1e-9)
        assert result["keys"][key]["max_drift"] == pytest.approx(0.1 * 1.5, abs=1e-9)

    def test_values_and_stubs_are_returned_sorted_regardless_of_input_order(self):
        stubs = [0.5, 0.0, 1.0]   # deliberately out of order
        transmissions = {0.0: 0.1, 0.5: 0.3, 1.0: 0.5}
        smats = [_smat(np.sqrt(transmissions[s]) + 0j) for s in stubs]
        result = cut_position_invariance(smats, stubs)
        key = ("out", 1, 1)
        assert list(result["keys"][key]["stubs"]) == [0.0, 0.5, 1.0]
        assert result["keys"][key]["values"] == pytest.approx([0.1, 0.3, 0.5])


class TestMultipleKeys:
    def test_worst_key_picks_the_largest_magnitude_slope(self):
        stubs = [0.0, 1.0]
        smats = [SMatrix(wavelengths=WAVELENGTHS) for _ in stubs]
        # "through": transmission 0.5 -> 0.5 (flat). "cross": 0.1 -> 0.4 (steep).
        for i, s in enumerate(stubs):
            smats[i].set(np.array([np.sqrt(0.5) + 0j]), port="through", in_mode=1, out_mode=1)
            cross_t = 0.1 + 0.3 * s
            smats[i].set(np.array([np.sqrt(cross_t) + 0j]), port="cross", in_mode=1, out_mode=1)
        result = cut_position_invariance(smats, stubs)
        assert result["worst_key"] == ("cross", 1, 1)
        assert result["worst_slope"] == pytest.approx(0.3, abs=1e-9)
        assert result["keys"][("through", 1, 1)]["slope"] == pytest.approx(0.0, abs=1e-9)


class TestValidation:
    def test_raises_on_mismatched_lengths(self):
        with pytest.raises(ValueError, match="same length"):
            cut_position_invariance([_smat(1.0), _smat(1.0)], [0.0])

    def test_raises_on_fewer_than_two_points(self):
        with pytest.raises(ValueError, match="at least 2"):
            cut_position_invariance([_smat(1.0)], [0.0])

    def test_raises_on_mismatched_key_sets(self):
        a = _smat(1.0, port="out", in_mode=1, out_mode=1)
        b = SMatrix(wavelengths=WAVELENGTHS)
        b.set(np.array([1.0]), port="out", in_mode=1, out_mode=2)  # different key
        with pytest.raises(ValueError, match="same .* keys"):
            cut_position_invariance([a, b], [0.0, 0.5])


class TestSMatrixKeys:
    def test_keys_returns_every_stored_entry(self):
        s = SMatrix(wavelengths=WAVELENGTHS)
        s.set(np.array([1.0]), port="out", in_mode=1, out_mode=1)
        s.set(np.array([0.5]), port="cross", in_mode=1, out_mode=1)
        assert set(s.keys()) == {("out", 1, 1), ("cross", 1, 1)}
