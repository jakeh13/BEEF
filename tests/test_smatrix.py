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
