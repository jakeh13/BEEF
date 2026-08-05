"""Fast unit tests for beef.results.smatrix.SMatrix — no solver/meep needed."""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from beef.results.smatrix import SMatrix


def _wavelengths():
    return np.array([1.53, 1.55, 1.57])


class TestBasics:
    def test_set_and_coeff_roundtrip(self):
        s = SMatrix(wavelengths=_wavelengths())
        coeff = np.array([0.1 + 0.2j, 0.3 - 0.1j, 0.0 + 0.5j])
        s.set(coeff, port="through", in_mode=1, out_mode=1)
        np.testing.assert_allclose(s.coeff(port="through"), coeff)

    def test_transmission_is_magnitude_squared(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set(np.full(3, 0.5 + 0.0j), port="out")
        np.testing.assert_allclose(s.transmission(port="out"), np.full(3, 0.25))

    def test_insertion_loss_db(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set(np.full(3, np.sqrt(0.5) + 0.0j), port="out")
        np.testing.assert_allclose(s.insertion_loss_db(port="out"), np.full(3, -3.0103), atol=1e-3)

    def test_set_power_stores_sqrt_as_real_coeff(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.25), port="out")
        np.testing.assert_allclose(s.transmission(port="out"), np.full(3, 0.25))

    def test_splitting_ratio_fifty_fifty(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.4), port="a")
        s.set_power(np.full(3, 0.4), port="b")
        np.testing.assert_allclose(s.splitting_ratio("a", "b"), np.full(3, 0.5))

    def test_ports_lists_unique_sorted_labels(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set(np.zeros(3, dtype=complex), port="cross")
        s.set(np.zeros(3, dtype=complex), port="through")
        assert s.ports() == ["cross", "through"]


class TestPassivity:
    def test_check_passive_ok_no_warning(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.5), port="a", in_mode=1)
        s.set_power(np.full(3, 0.5), port="b", in_mode=1)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert s.check_passive() is True

    def test_check_passive_violation_warns_and_returns_false(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.9), port="a", in_mode=1)
        s.set_power(np.full(3, 0.9), port="b", in_mode=1)  # total 1.8 >> 1
        with pytest.warns(UserWarning, match="check_passive"):
            assert s.check_passive() is False

    def test_total_power_missing_in_mode_raises(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.5), port="a", in_mode=1)
        with pytest.raises(KeyError):
            s.total_power(in_mode=2)


class TestCoherentReceiverMetrics:
    """imbalance_db/excess_loss_db/cmrr_db/cross_excitation_phase_diff_deg --
    added for the 4x4 MMI 90-degree hybrid (mmi_4x4_90_hybrid_validation.md
    finding A8). Synthetic data only, hand-computable, no solver needed."""

    def test_imbalance_db_balanced_is_zero(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.25), port="a")
        s.set_power(np.full(3, 0.25), port="b")
        np.testing.assert_allclose(s.imbalance_db(["a", "b"]), np.zeros(3), atol=1e-10)

    def test_imbalance_db_unbalanced_matches_hand_calc(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.4), port="a")
        s.set_power(np.full(3, 0.1), port="b")
        # -10*log10(0.1/0.4) = -10*log10(0.25) = 6.0206 dB
        np.testing.assert_allclose(s.imbalance_db(["a", "b"]), np.full(3, 6.0206), atol=1e-3)

    def test_imbalance_db_generalizes_to_more_than_two_ports(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.4), port="a")
        s.set_power(np.full(3, 0.25), port="b")
        s.set_power(np.full(3, 0.1), port="c")
        np.testing.assert_allclose(s.imbalance_db(["a", "b", "c"]), np.full(3, 6.0206), atol=1e-3)

    def test_imbalance_db_needs_at_least_two_ports(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.4), port="a")
        with pytest.raises(ValueError):
            s.imbalance_db(["a"])

    def test_excess_loss_db_matches_hand_calc(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.4), port="a")
        s.set_power(np.full(3, 0.4), port="b")
        # -10*log10(0.8/1.0) = 0.9691 dB, positive = loss (paper's own sign
        # convention -- the OPPOSITE of insertion_loss_db()'s).
        np.testing.assert_allclose(s.excess_loss_db(["a", "b"]), np.full(3, 0.9691), atol=1e-3)

    def test_excess_loss_db_zero_when_lossless(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.6), port="a")
        s.set_power(np.full(3, 0.4), port="b")
        np.testing.assert_allclose(s.excess_loss_db(["a", "b"]), np.zeros(3), atol=1e-10)

    def test_cmrr_db_matches_hand_calc(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.6), port="a")
        s.set_power(np.full(3, 0.4), port="b")
        # -20*log10(|0.6-0.4|/(0.6+0.4)) = -20*log10(0.2) = 13.9794 dBe
        np.testing.assert_allclose(s.cmrr_db("a", "b"), np.full(3, 13.9794), atol=1e-3)

    def test_cmrr_db_near_perfect_rejection_is_large(self):
        s = SMatrix(wavelengths=_wavelengths())
        s.set_power(np.full(3, 0.5), port="a")
        s.set_power(np.full(3, 0.5), port="b")
        assert np.all(s.cmrr_db("a", "b") > 100.0)

    def test_cross_excitation_phase_diff_simple_offset(self):
        s_in1 = SMatrix(wavelengths=_wavelengths())
        s_in3 = SMatrix(wavelengths=_wavelengths())
        s_in1.set(np.exp(1j * np.radians(0.0)) * np.ones(3), port="out6")
        s_in3.set(np.exp(1j * np.radians(30.0)) * np.ones(3), port="out6")
        phi = s_in1.cross_excitation_phase_diff_deg(s_in3, port="out6")
        np.testing.assert_allclose(phi, np.full(3, 30.0), atol=1e-6)

    def test_cross_excitation_phase_diff_wraps_to_180(self):
        s_in1 = SMatrix(wavelengths=_wavelengths())
        s_in3 = SMatrix(wavelengths=_wavelengths())
        s_in1.set(np.exp(1j * np.radians(170.0)) * np.ones(3), port="out6")
        s_in3.set(np.exp(1j * np.radians(-170.0)) * np.ones(3), port="out6")
        # raw difference is -340 deg == +20 deg once wrapped
        phi = s_in1.cross_excitation_phase_diff_deg(s_in3, port="out6")
        np.testing.assert_allclose(phi, np.full(3, 20.0), atol=1e-6)

    def test_cross_excitation_phase_diff_matches_paper_eq15_targets_in_range(self):
        # phi~_5 = -45 deg, phi~_7 = 45 deg (both already inside (-180,180],
        # unlike phi~_6/phi~_8 -- see the method's own mod-360 caveat).
        s_in1 = SMatrix(wavelengths=_wavelengths())
        s_in3 = SMatrix(wavelengths=_wavelengths())
        s_in1.set(np.ones(3, dtype=complex), port="out5")
        s_in3.set(np.exp(1j * np.radians(-45.0)) * np.ones(3), port="out5")
        s_in1.set(np.ones(3, dtype=complex), port="out7")
        s_in3.set(np.exp(1j * np.radians(45.0)) * np.ones(3), port="out7")
        np.testing.assert_allclose(
            s_in1.cross_excitation_phase_diff_deg(s_in3, port="out5"), np.full(3, -45.0), atol=1e-6)
        np.testing.assert_allclose(
            s_in1.cross_excitation_phase_diff_deg(s_in3, port="out7"), np.full(3, 45.0), atol=1e-6)


class TestSerialization:
    def test_save_load_roundtrip(self, tmp_path):
        s = SMatrix(wavelengths=_wavelengths())
        s.set(np.array([0.1 + 0.2j, 0.3, 0.4 - 0.1j]), port="out", in_mode=1, out_mode=1)
        path = tmp_path / "result.pkl"
        s.save(str(path))
        loaded = SMatrix.load(str(path))
        np.testing.assert_allclose(loaded.coeff(port="out"), s.coeff(port="out"))
        np.testing.assert_allclose(loaded.wavelengths, s.wavelengths)

    def test_from_legacy_powers(self):
        s = SMatrix.from_legacy_powers(
            wavelengths=_wavelengths(),
            powers_ports=[[0.5, 0.5, 0.5], [0.25, 0.25, 0.25]],
            port_names=["through", "cross"],
        )
        np.testing.assert_allclose(s.transmission(port="through"), np.full(3, 0.5))
        np.testing.assert_allclose(s.transmission(port="cross"), np.full(3, 0.25))
