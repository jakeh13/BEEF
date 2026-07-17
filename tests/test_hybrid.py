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
LENGTHS = np.array([0.5, 1.0, 2.0])


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


class FakeSweepSolver(FakeSolver):
    """Adds a scripted solve_section_length_sweep: slice L is (1+L) * base
    matrix, so per-length expectations are trivially computable in tests."""

    def __init__(self, matrix):
        super().__init__(matrix)
        self.sweep_call_count = 0

    def solve_section_length_sweep(self, geometry, wavelengths, x0, x1, lengths,
                                   n_modes_in, n_modes_out,
                                   stub_left=0.0, stub_right=0.0, **kwargs):
        self.sweep_call_count += 1
        M = len(wavelengths)
        base = np.broadcast_to(self.matrix, (M,) + self.matrix.shape)
        return np.stack([(1.0 + L) * base for L in lengths]).astype(complex)


class FakeGeometry:
    """Straight waveguide of `width` centered at `center`, symmetric about
    x=0 by construction — for exercising mirror_pairs' geometric guardrail
    without meep/meow. `asymmetric_beyond` deliberately breaks symmetry past
    a given x, to test that the guardrail catches a bad pairing."""

    def __init__(self, width=0.45, center=0.0, asymmetric_beyond=None):
        self.width = width
        self.center = center
        self.asymmetric_beyond = asymmetric_beyond

    def width_at(self, x, **kwargs):
        if self.asymmetric_beyond is not None and x > self.asymmetric_beyond:
            return self.width * 2
        return self.width

    def center_at(self, x, **kwargs):
        return self.center


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


class TestMirrorPairsConstruction:
    def _solvers(self, n):
        return [FakeSolver([[1.0]]) for _ in range(n)]

    def test_out_of_range_index_raises(self):
        with pytest.raises(ValueError, match="valid section indices"):
            HybridSolver(self._solvers(3), boundaries=[0.0, 1.0, 2.0, 3.0],
                         interface_modes=[1, 1], mirror_pairs={5: 0})

    def test_self_pair_raises(self):
        with pytest.raises(ValueError, match="cannot mirror itself"):
            HybridSolver(self._solvers(3), boundaries=[0.0, 1.0, 2.0, 3.0],
                         interface_modes=[1, 1], mirror_pairs={1: 1})

    def test_forward_reference_raises(self):
        with pytest.raises(ValueError, match="must be < "):
            HybridSolver(self._solvers(3), boundaries=[0.0, 1.0, 2.0, 3.0],
                         interface_modes=[1, 1], mirror_pairs={0: 2})

    def test_chained_mirror_raises(self):
        with pytest.raises(ValueError, match="chained/transitive"):
            HybridSolver(self._solvers(4), boundaries=[0.0, 1.0, 2.0, 3.0, 4.0],
                         interface_modes=[1, 1, 1], mirror_pairs={1: 0, 2: 1})

    def test_mode_count_mismatch_raises_at_construction(self):
        with pytest.raises(ValueError, match="mode-count mismatch"):
            HybridSolver(self._solvers(4), boundaries=[0.0, 1.0, 2.0, 3.0, 4.0],
                         interface_modes=[2, 5, 3], mirror_pairs={2: 1})

    def test_mirror_pairs_none_is_noop(self):
        s1, s2 = FakeSolver([[2.0]]), FakeSolver([[3.0]])
        hybrid = HybridSolver([s1, s2], boundaries=[0.0, 5.0, 10.0], interface_modes=[1])
        assert hybrid.mirror_pairs == {}
        result = hybrid.solve(geometry=None, wavelengths=WAVELENGTHS)
        np.testing.assert_allclose(
            result.coeff(port="out", in_mode=1, out_mode=1), np.full(3, 6.0 + 0j)
        )


class TestMirrorPairs:
    def _make_hybrid(self, mirror_pairs=None):
        s0 = FakeSolver([[1.0], [0.5]])             # section 0: 1 -> 2 modes
        s1 = FakeSolver([[1.0, 0.0], [0.0, 1.0]])   # body: 2 -> 2 modes (identity)
        s2 = FakeSolver([[99.0, 99.0]])             # section 2: must never be called
        hybrid = HybridSolver(
            [s0, s1, s2], boundaries=[-2.0, -1.0, 1.0, 2.0],
            interface_modes=[2, 2], mirror_pairs=mirror_pairs,
        )
        return hybrid, s0, s1, s2

    def test_mirror_skips_solve_section_and_matches_transpose(self):
        geom = FakeGeometry(width=0.45)
        hybrid, s0, s1, s2 = self._make_hybrid(mirror_pairs={2: 0})
        result = hybrid.solve(geometry=geom, wavelengths=WAVELENGTHS)

        assert s2.call_count == 0
        # manual expectation: mats[2] == mats[0].T, cascade = mats[2] @ mats[1] @ mats[0]
        expected = s0.matrix.T @ s1.matrix @ s0.matrix  # [1,1]
        coeff = result.coeff(port="out", in_mode=1, out_mode=1)
        np.testing.assert_allclose(coeff, np.full(3, expected[0, 0]))

    def test_mirror_guardrail_raises_on_asymmetric_geometry(self):
        geom = FakeGeometry(width=0.45, asymmetric_beyond=0.5)
        hybrid, s0, s1, s2 = self._make_hybrid(mirror_pairs={2: 0})
        with pytest.raises(ValueError, match="width_at mismatch"):
            hybrid.solve(geometry=geom, wavelengths=WAVELENGTHS)

    def test_mirror_reused_section_populates_cache(self):
        geom = FakeGeometry(width=0.45)
        hybrid, s0, s1, s2 = self._make_hybrid(mirror_pairs={2: 0})
        cache: dict = {}

        hybrid.solve(geometry=geom, wavelengths=WAVELENGTHS, cache=cache)
        assert len(cache) == 3
        assert s0.call_count == 1 and s1.call_count == 1 and s2.call_count == 0

        hybrid.solve(geometry=geom, wavelengths=WAVELENGTHS, cache=cache)
        assert s0.call_count == 1 and s1.call_count == 1 and s2.call_count == 0


class TestLengthSweep:
    def _make_chain(self, mirror_pairs=None, s2_matrix=None):
        s0 = FakeSolver([[1.0], [0.5]])                     # 1 -> 2 modes
        s1 = FakeSweepSolver([[1.0, 2.0], [0.0, 1.0]])      # 2 -> 2 modes, swept
        s2 = FakeSolver(s2_matrix if s2_matrix is not None else [[3.0, 4.0]])  # 2 -> 1
        hybrid = HybridSolver(
            [s0, s1, s2], boundaries=[-2.0, -1.0, 1.0, 2.0],
            interface_modes=[2, 2], mirror_pairs=mirror_pairs,
        )
        return hybrid, s0, s1, s2

    def test_sweep_matches_manual_broadcast_cascade(self):
        geom = FakeGeometry()
        hybrid, s0, s1, s2 = self._make_chain()
        grid = hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS, sweep_section=1)

        assert grid.shape == (3, 3, 1, 1)
        for L_idx, L in enumerate(LENGTHS):
            expected = s2.matrix @ ((1.0 + L) * s1.matrix) @ s0.matrix
            np.testing.assert_allclose(grid[L_idx, :, 0, 0], np.full(3, expected[0, 0]))

    def test_lless_sections_cache_share_with_solve(self):
        geom = FakeGeometry()
        hybrid, s0, s1, s2 = self._make_chain()
        cache: dict = {}

        hybrid.solve(geometry=geom, wavelengths=WAVELENGTHS, in_modes=1, out_modes=1,
                     cache=cache)
        hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS, sweep_section=1,
                                  cache=cache)

        assert s0.call_count == 1   # cache hit in the sweep call
        assert s2.call_count == 1   # cache hit in the sweep call
        assert s1.call_count == 1        # solve()'s single-length entry
        assert s1.sweep_call_count == 1  # sweep's lengths-keyed entry, distinct

    def test_swept_cache_key_includes_lengths(self):
        geom = FakeGeometry()
        hybrid, s0, s1, s2 = self._make_chain()
        cache: dict = {}

        hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS, sweep_section=1, cache=cache)
        hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS, sweep_section=1, cache=cache)
        assert s1.sweep_call_count == 1  # identical lengths -> cache hit

        hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS * 2, sweep_section=1, cache=cache)
        assert s1.sweep_call_count == 2  # different lengths -> new entry

        hybrid.solve(geometry=geom, wavelengths=WAVELENGTHS, cache=cache)
        assert s1.call_count == 1  # solve() never sees the lengths-keyed 4-D entries

    def test_sweep_section_without_sweep_method_raises(self):
        geom = FakeGeometry()
        hybrid, s0, s1, s2 = self._make_chain()
        with pytest.raises(TypeError, match="solve_section_length_sweep"):
            hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS, sweep_section=0)

    def test_sweep_section_out_of_range_raises(self):
        geom = FakeGeometry()
        hybrid, s0, s1, s2 = self._make_chain()
        with pytest.raises(ValueError, match="valid section index"):
            hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS, sweep_section=5)

    def test_sweep_section_in_mirror_pairs_raises(self):
        geom = FakeGeometry()
        # all-sweep-capable solvers so the capability check can't fire first
        s0 = FakeSweepSolver([[1.0], [0.5]])
        s1 = FakeSweepSolver([[1.0, 2.0], [0.0, 1.0]])
        s2 = FakeSweepSolver([[3.0, 4.0]])
        hybrid = HybridSolver(
            [s0, s1, s2], boundaries=[-2.0, -1.0, 1.0, 2.0],
            interface_modes=[2, 2], mirror_pairs={2: 0},
        )
        with pytest.raises(ValueError, match="mirror_pairs"):
            hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS, sweep_section=2)
        with pytest.raises(ValueError, match="mirror_pairs"):
            hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS, sweep_section=0)

    def test_mirror_partner_stays_lless_in_sweep(self):
        geom = FakeGeometry()
        hybrid, s0, s1, s2 = self._make_chain(mirror_pairs={2: 0},
                                              s2_matrix=[[99.0, 99.0]])
        grid = hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS, sweep_section=1)

        assert s2.call_count == 0
        for L_idx, L in enumerate(LENGTHS):
            expected = s0.matrix.T @ ((1.0 + L) * s1.matrix) @ s0.matrix
            np.testing.assert_allclose(grid[L_idx, :, 0, 0], np.full(3, expected[0, 0]))

    def test_nonuniform_sweep_region_raises(self):
        geom = FakeGeometry(width=0.45, asymmetric_beyond=0.0)  # width doubles mid-span
        hybrid, s0, s1, s2 = self._make_chain()
        with pytest.raises(ValueError, match="not uniform"):
            hybrid.solve_length_sweep(geom, WAVELENGTHS, LENGTHS, sweep_section=1)
        assert s1.sweep_call_count == 0  # validation precedes solving
