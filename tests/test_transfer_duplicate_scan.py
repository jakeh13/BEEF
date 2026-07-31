"""Fast, meep/meow-free tests for beef.mode_identity.transfer_matrix_duplicate_scan()
/ excited_mode_axis() / warn_on_transfer_duplicate() -- the post-hoc
solve_section() duplicate-mode scan (fdtd_vs_EME_mode_validation.md item 8).

Pure numpy: these operate on plain solve_section()-shaped transfer matrices
([n_wl, n_out, n_in] complex arrays), not a real eigensolve, the same
"guard the checker itself" pattern tests/test_mode_identity_overlap.py and
tests/test_hybrid_mode_correspondence.py already use for the field-overlap
machinery.
"""
from __future__ import annotations

import logging

import numpy as np
import pytest

from beef.mode_identity import (
    excited_mode_axis, transfer_matrix_duplicate_scan, warn_on_transfer_duplicate,
)


class TestExcitedModeAxis:
    def test_fewer_input_modes_excites_the_input_axis(self):
        assert excited_mode_axis(n_modes_out=16, n_modes_in=3) == 2

    def test_fewer_output_modes_excites_the_output_axis(self):
        assert excited_mode_axis(n_modes_out=3, n_modes_in=16) == 1

    def test_equal_counts_default_to_the_input_axis(self):
        """solve_section()'s own tie-break (n_out <= n_in favors exciting
        input) -- matched here so the two never disagree about which side a
        square transfer matrix's duplicate scan should run over."""
        assert excited_mode_axis(n_modes_out=4, n_modes_in=4) == 2


def _random_vector(rng, n_wl, n_other):
    return rng.normal(size=(n_wl, n_other)) + 1j * rng.normal(size=(n_wl, n_other))


class TestTransferMatrixDuplicateScan:
    def test_distinct_random_columns_are_not_flagged(self):
        rng = np.random.default_rng(0)
        n_wl, n_out, n_in = 3, 8, 3
        T = rng.normal(size=(n_wl, n_out, n_in)) + 1j * rng.normal(size=(n_wl, n_out, n_in))
        scan = transfer_matrix_duplicate_scan(T)
        assert scan["excited_axis"] == 2
        assert scan["n_excited"] == n_in
        assert not scan["has_duplicate"]
        assert scan["max_similarity"] < 0.9

    def test_two_parallel_columns_are_flagged_as_duplicate(self):
        rng = np.random.default_rng(1)
        n_wl, n_out = 4, 10
        base = _random_vector(rng, n_wl, n_out)
        other = _random_vector(rng, n_wl, n_out)
        # T[:, :, 0] and T[:, :, 2] are the SAME vector up to a per-wavelength
        # complex scale -- exactly what a duplicate excited eigenmode would
        # produce (same physical mode excited twice under different indices).
        T = np.stack([base, other, 2.5j * base], axis=2)  # [n_wl, n_out, n_in=3]
        scan = transfer_matrix_duplicate_scan(T)
        assert scan["has_duplicate"]
        assert scan["worst_pair"] == (0, 2)
        assert scan["max_similarity"] == pytest.approx(1.0, abs=1e-9)

    def test_similarity_matrix_is_symmetric_with_unit_diagonal(self):
        rng = np.random.default_rng(2)
        T = rng.normal(size=(2, 5, 4)) + 1j * rng.normal(size=(2, 5, 4))
        scan = transfer_matrix_duplicate_scan(T)
        sim = scan["similarity"]
        assert sim.shape == (4, 4)
        assert np.allclose(sim, sim.T)
        assert np.allclose(np.diag(sim), 1.0)

    def test_a_zero_column_does_not_produce_a_spurious_duplicate(self):
        """A mode with zero response (e.g. genuinely uncoupled) must not
        divide-by-zero its way into a false 'duplicate' -- similarity against
        a zero vector is defined as 0, not NaN/inf."""
        rng = np.random.default_rng(3)
        base = _random_vector(rng, 3, 6)
        zero = np.zeros_like(base)
        T = np.stack([base, zero], axis=2)
        scan = transfer_matrix_duplicate_scan(T)
        assert not scan["has_duplicate"]
        assert np.isfinite(scan["similarity"]).all()

    def test_single_excited_mode_has_no_pair_to_report(self):
        T = np.ones((2, 3, 1), dtype=complex)
        scan = transfer_matrix_duplicate_scan(T)
        assert scan["n_excited"] == 1
        assert scan["worst_pair"] is None
        assert scan["max_similarity"] == 0.0
        assert not scan["has_duplicate"]

    def test_custom_threshold_is_honored(self):
        rng = np.random.default_rng(4)
        base = _random_vector(rng, 3, 6)
        # A mild correlation (~0.5), not a real duplicate, but enough to trip
        # a deliberately loose custom threshold.
        other = base + 3.0 * _random_vector(rng, 3, 6)
        T = np.stack([base, other], axis=2)
        scan_default = transfer_matrix_duplicate_scan(T)
        scan_loose = transfer_matrix_duplicate_scan(T, threshold=0.01)
        assert not scan_default["has_duplicate"]
        assert scan_loose["has_duplicate"]

    def test_rejects_non_3d_input(self):
        with pytest.raises(ValueError, match="3-D"):
            transfer_matrix_duplicate_scan(np.zeros((3, 3)))


class TestWarnOnTransferDuplicate:
    def _logger(self):
        return logging.getLogger("beef.test_transfer_duplicate_scan")

    def test_logs_a_warning_when_a_duplicate_is_found(self, caplog):
        rng = np.random.default_rng(5)
        base = _random_vector(rng, 3, 6)
        T = np.stack([base, 1.7 * base], axis=2)
        with caplog.at_level(logging.WARNING, logger="beef.test_transfer_duplicate_scan"):
            scan = warn_on_transfer_duplicate(T, self._logger(), x0=1.0, x1=2.0)
        assert scan["has_duplicate"]
        assert any("possible duplicate excited mode" in r.message for r in caplog.records)
        assert any("x=[1, 2]" in r.message for r in caplog.records)

    def test_does_not_log_when_no_duplicate_is_found(self, caplog):
        rng = np.random.default_rng(6)
        T = rng.normal(size=(3, 8, 3)) + 1j * rng.normal(size=(3, 8, 3))
        with caplog.at_level(logging.WARNING, logger="beef.test_transfer_duplicate_scan"):
            scan = warn_on_transfer_duplicate(T, self._logger(), x0=0.0, x1=1.0)
        assert not scan["has_duplicate"]
        assert not any("possible duplicate excited mode" in r.message for r in caplog.records)

    def test_returns_the_scan_dict_for_further_inspection(self):
        rng = np.random.default_rng(7)
        T = rng.normal(size=(2, 4, 2)) + 1j * rng.normal(size=(2, 4, 2))
        scan = warn_on_transfer_duplicate(T, self._logger(), x0=0.0, x1=1.0)
        assert set(scan) == {
            "excited_axis", "n_excited", "similarity", "worst_pair",
            "max_similarity", "has_duplicate",
        }
