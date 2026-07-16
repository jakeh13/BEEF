"""Fast unit tests for feces.solvers.hybrid.HybridSolver.

Uses fake solvers with a scripted solve_section() — no real FDTD/EME, no
meep/meow import — to check the matrix-chain cascade math, geometry_kwargs
fan-out, and the section cache independently of any real physics.
"""
from __future__ import annotations

import numpy as np
import pytest

from feces.solvers.hybrid import HybridSolver

WAVELENGTHS = np.array([1.53, 1.55, 1.57])


class FakeSolver:
    """Returns a fixed [M, n_out, n_in] matrix (broadcast over wavelength)
    from solve_section(), and counts how many times it's actually called."""

    def __init__(self, matrix):
        self.matrix = np.asarray(matrix, dtype=complex)
        self.call_count = 0

    def solve_section(self, geometry, wavelengths, x0, x1,
                       n_modes_in, n_modes_out, stub_left=0.0, stub_right=0.0,
                       **kwargs):
        self.call_count += 1
        M = len(wavelengths)
        return np.broadcast_to(self.matrix, (M,) + self.matrix.shape).copy()


class TestConstruction:
    def test_requires_at_least_two_solvers(self):
        with pytest.raises(ValueError, match="at least 2 solvers"):
            HybridSolver([FakeSolver([[1.0]])], boundaries=[0.0, 1.0], interface_modes=[])

    def test_boundaries_length_mismatch_raises(self):
        s1, s2 = FakeSolver([[1.0]]), FakeSolver([[1.0]])
        with pytest.raises(ValueError, match="boundaries must have"):
            HybridSolver([s1, s2], boundaries=[0.0, 1.0], interface_modes=[1])

    def test_interface_modes_length_mismatch_raises(self):
        s1, s2 = FakeSolver([[1.0]]), FakeSolver([[1.0]])
        with pytest.raises(ValueError, match="interface_modes must have"):
            HybridSolver([s1, s2], boundaries=[0.0, 1.0, 2.0], interface_modes=[1, 1])


class TestCascade:
    def test_two_section_scalar_cascade(self):
        # single-mode-in/out chain: full = mat2 @ mat1 (scalar case, order
        # doesn't matter numerically but exercises the einsum contraction).
        s1 = FakeSolver([[2.0]])
        s2 = FakeSolver([[3.0]])
        hybrid = HybridSolver([s1, s2], boundaries=[0.0, 5.0, 10.0], interface_modes=[1])
        result = hybrid.solve(geometry=None, wavelengths=WAVELENGTHS)
        coeff = result.coeff(port="out", in_mode=1, out_mode=1)
        np.testing.assert_allclose(coeff, np.full(3, 6.0 + 0j))

    def test_three_section_mode_count_change_at_each_interface(self):
        # interface_modes=[1, 2]: section 0 is 1-in/1-out, section 1 expands
        # 1 mode -> 2 modes, section 2 is a non-diagonal 2x2 mixing matrix.
        # Regression guard on the cascade's per-section shapes and the
        # einsum contraction order (mat_k @ ... @ mat_1).
        s1 = FakeSolver([[1.0]])
        m2 = np.array([[1.0], [0.5]])            # shape (n_out=2, n_in=1)
        m3 = np.array([[1.0, 2.0], [0.0, 1.0]])  # shape (n_out=2, n_in=2)
        s2 = FakeSolver(m2)
        s3 = FakeSolver(m3)
        hybrid = HybridSolver(
            [s1, s2, s3],
            boundaries=[0.0, 1.0, 2.0, 3.0],
            interface_modes=[1, 2],
        )
        result = hybrid.solve(geometry=None, wavelengths=WAVELENGTHS, in_modes=1, out_modes=2)

        expected = m3 @ m2 @ np.array([[1.0]])  # [2,1]
        for out_mode in (1, 2):
            coeff = result.coeff(port="out", in_mode=1, out_mode=out_mode)
            np.testing.assert_allclose(coeff, np.full(3, expected[out_mode - 1, 0]))

    def test_result_stored_under_requested_port_label(self):
        s1, s2 = FakeSolver([[1.0]]), FakeSolver([[1.0]])
        hybrid = HybridSolver([s1, s2], boundaries=[0.0, 1.0, 2.0], interface_modes=[1])
        result = hybrid.solve(geometry=None, wavelengths=WAVELENGTHS, port="cross")
        assert result.ports() == ["cross"]


class TestGeometryKwargsFanout:
    def test_single_dict_applied_to_every_section(self):
        seen = []

        class RecordingSolver(FakeSolver):
            def solve_section(self, geometry, wavelengths, x0, x1,
                               n_modes_in, n_modes_out, stub_left=0.0, stub_right=0.0,
                               **kwargs):
                seen.append(kwargs)
                return super().solve_section(geometry, wavelengths, x0, x1,
                                              n_modes_in, n_modes_out, **kwargs)

        s1, s2 = RecordingSolver([[1.0]]), RecordingSolver([[1.0]])
        hybrid = HybridSolver([s1, s2], boundaries=[0.0, 1.0, 2.0], interface_modes=[1])
        hybrid.solve(geometry=None, wavelengths=WAVELENGTHS, geometry_kwargs={"foo": "bar"})
        assert seen == [{"foo": "bar"}, {"foo": "bar"}]

    def test_list_gives_per_section_kwargs(self):
        seen = []

        class RecordingSolver(FakeSolver):
            def solve_section(self, geometry, wavelengths, x0, x1,
                               n_modes_in, n_modes_out, stub_left=0.0, stub_right=0.0,
                               **kwargs):
                seen.append(kwargs)
                return super().solve_section(geometry, wavelengths, x0, x1,
                                              n_modes_in, n_modes_out, **kwargs)

        s1, s2 = RecordingSolver([[1.0]]), RecordingSolver([[1.0]])
        hybrid = HybridSolver([s1, s2], boundaries=[0.0, 1.0, 2.0], interface_modes=[1])
        hybrid.solve(geometry=None, wavelengths=WAVELENGTHS,
                     geometry_kwargs=[{"output_idx": 0}, {"output_idx": 1}])
        assert seen == [{"output_idx": 0}, {"output_idx": 1}]

    def test_geometry_kwargs_list_wrong_length_raises(self):
        s1, s2 = FakeSolver([[1.0]]), FakeSolver([[1.0]])
        hybrid = HybridSolver([s1, s2], boundaries=[0.0, 1.0, 2.0], interface_modes=[1])
        with pytest.raises(ValueError, match="geometry_kwargs list must have"):
            hybrid.solve(geometry=None, wavelengths=WAVELENGTHS, geometry_kwargs=[{}])


class TestCache:
    def test_shared_cache_reuses_unchanged_sections_across_calls(self):
        s1, s2 = FakeSolver([[1.0]]), FakeSolver([[1.0]])
        hybrid = HybridSolver([s1, s2], boundaries=[0.0, 1.0, 2.0], interface_modes=[1])
        cache: dict = {}

        hybrid.solve(geometry=None, wavelengths=WAVELENGTHS,
                     geometry_kwargs=[{}, {"output_idx": 0}], cache=cache)
        hybrid.solve(geometry=None, wavelengths=WAVELENGTHS,
                     geometry_kwargs=[{}, {"output_idx": 1}], cache=cache)

        # s1's kwargs never changed between calls -> only ever simulated once.
        assert s1.call_count == 1
        # s2's kwargs differ each call -> simulated every time.
        assert s2.call_count == 2

    def test_without_cache_every_call_resimulates(self):
        s1, s2 = FakeSolver([[1.0]]), FakeSolver([[1.0]])
        hybrid = HybridSolver([s1, s2], boundaries=[0.0, 1.0, 2.0], interface_modes=[1])

        hybrid.solve(geometry=None, wavelengths=WAVELENGTHS)
        hybrid.solve(geometry=None, wavelengths=WAVELENGTHS)

        assert s1.call_count == 2
        assert s2.call_count == 2
