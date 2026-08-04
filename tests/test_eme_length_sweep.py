"""Tests for EMESolver.solve_length_sweep() — the mode-basis-reuse length
sweep (one mode solve per wavelength serves every device length).

The helper tests are meep/meow-free (EMESolver constructs without importing
either). The equivalence test is marked slow: it runs a real meow mode solve
and checks that a length swept via solve_length_sweep() reproduces, to
floating-point tolerance, what a per-length solve() computes from scratch.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

from beef import EMESolver, SOI_160, Taper, mesh_cell_centers


def _solver(num_cells, lengths, **kw):
    num_modes = [1] * len(num_cells)
    return EMESolver(num_modes=num_modes, num_cells=num_cells, lengths=lengths,
                     material_stack=SOI_160, **kw)


class TestLengthSweepHelpers:
    def test_flat_cell_lengths_divides_regions_evenly(self):
        s = _solver(num_cells=[1, 4, 1], lengths=[0.5, 4.0, 0.5])
        assert s._flat_cell_lengths() == [0.5, 1.0, 1.0, 1.0, 1.0, 0.5]

    def test_flat_cell_lengths_accepts_override(self):
        s = _solver(num_cells=[1, 4, 1], lengths=[0.5, 4.0, 0.5])
        assert s._flat_cell_lengths([0.5, 8.0, 0.5]) == [0.5, 2.0, 2.0, 2.0, 2.0, 0.5]

    def test_region_cell_span(self):
        s = _solver(num_cells=[1, 4, 1], lengths=[0.5, 4.0, 0.5])
        assert s._region_cell_span(0) == (0, 1)
        assert s._region_cell_span(1) == (1, 5)
        assert s._region_cell_span(2) == (5, 6)

    def test_infer_sweep_region_unique(self):
        s = _solver(num_cells=[1, 8, 1], lengths=[0.5, 4.0, 0.5])
        assert s._infer_sweep_region() == 1

    def test_infer_sweep_region_ambiguous_raises(self):
        s = _solver(num_cells=[1, 4, 4, 1], lengths=[0.5, 2.0, 2.0, 0.5])
        with pytest.raises(ValueError):
            s._infer_sweep_region()

    def test_infer_sweep_region_none_raises(self):
        s = _solver(num_cells=[1, 1, 1], lengths=[0.5, 4.0, 0.5])
        with pytest.raises(ValueError):
            s._infer_sweep_region()


@pytest.mark.slow
class TestLengthSweepEquivalence:
    """solve_length_sweep() at length L must equal solve() on a solver built
    directly at length L — same cells, same modes, only cell_lengths differ."""

    def test_matches_per_length_solve_for_taper(self):
        wavelengths = np.linspace(1.30, 1.32, 2)
        stub, n_mid = 0.5, 3
        sweep_lengths = [2.0, 5.0]

        # Reference solver drawn at an arbitrary length (1.0) -- its swept
        # region length must not affect the result.
        ref_taper = Taper(length=1.0, width_start=0.35, width_end=0.8)
        ref_solver = EMESolver(
            num_modes=[1, 3, 1], num_cells=[1, n_mid, 1],
            lengths=[stub, 1.0, stub],
            oxide_width=4.0, oxide_height=1.5, res=20.0,
            material_stack=SOI_160, mode_solver="tidy3d",
        )
        swept = ref_solver.solve_length_sweep(
            ref_taper, wavelengths, sweep_lengths, ports=["through"],
        )

        for L, smat_swept in zip(sweep_lengths, swept):
            direct_taper = Taper(length=L, width_start=0.35, width_end=0.8)
            direct_solver = EMESolver(
                num_modes=[1, 3, 1], num_cells=[1, n_mid, 1],
                lengths=[stub, L, stub],
                oxide_width=4.0, oxide_height=1.5, res=20.0,
                material_stack=SOI_160, mode_solver="tidy3d",
            )
            smat_direct = direct_solver.solve(direct_taper, wavelengths, ports=["through"])
            np.testing.assert_allclose(
                smat_swept.coeff(port="through"),
                smat_direct.coeff(port="through"),
                rtol=1e-6, atol=1e-9,
                err_msg=f"length-sweep result at L={L} disagrees with a direct solve()",
            )


@pytest.mark.slow
class TestSolveFingerprintsMaterial:
    """Real, `.solve()`-executing confirmation of the fix in
    `tests/test_cache_config.py::TestEMESolveFingerprintsMaterial` (which is
    meep/meow-free by that file's own convention and only checks the fix is
    wired into the source). This is the end-to-end evidence: two real
    `EMESolver.solve(save_dir=...)` calls, differing ONLY in
    `material_stack.core_kind` at the SAME wg_height/sz -- the exact
    collision `_solver_config()`'s bare `{wg_height, sz}` would not have
    caught, and which `solve()` had no OTHER fingerprint contribution to catch
    with either, before 2026-08-03 (`mmi_1x4_sin300_validation.md` step 3).
    """

    def test_different_core_kind_same_dims_does_not_cache_collide(self, tmp_path):
        from beef.materials.stack import MaterialStack

        geometry = Taper(length=1.0, width_start=0.8, width_end=0.8)  # straight
        wavelengths = np.linspace(1.50, 1.51, 2)

        si = MaterialStack(wg_height=0.30, sz=2.0, core_kind="si")
        sin = MaterialStack(wg_height=0.30, sz=2.0, core_kind="si3n4")
        assert (si.wg_height, si.sz) == (sin.wg_height, sin.sz), (
            "the two stacks must agree on everything _solver_config() alone "
            "records, or this test isn't isolating the fixed gap"
        )

        save_dir = str(tmp_path)
        solver_si = EMESolver(num_modes=[1], num_cells=[1], lengths=[1.0],
                              oxide_width=3.0, oxide_height=2.0, res=20.0,
                              material_stack=si, mode_solver="tidy3d")
        solver_sin = EMESolver(num_modes=[1], num_cells=[1], lengths=[1.0],
                               oxide_width=3.0, oxide_height=2.0, res=20.0,
                               material_stack=sin, mode_solver="tidy3d")

        solver_si.solve(geometry, wavelengths, ports=["through"], save_dir=save_dir)
        solver_sin.solve(geometry, wavelengths, ports=["through"], save_dir=save_dir)

        run_dirs = [d for d in os.listdir(save_dir)
                    if os.path.isdir(os.path.join(save_dir, d))]
        assert len(run_dirs) == 2, (
            f"expected two distinct cache entries (si-core, si3n4-core), got "
            f"{len(run_dirs)}: {run_dirs} -- the material collision is back"
        )

    def test_cache_hit_on_the_second_identical_call(self, tmp_path):
        """Companion to the test above: confirms the fix didn't overshoot into
        making EVERY call a cache miss (e.g. by leaking something
        non-deterministic like a timestamp into the fingerprint)."""
        geometry = Taper(length=1.0, width_start=0.8, width_end=0.8)
        wavelengths = np.linspace(1.50, 1.51, 2)
        save_dir = str(tmp_path)

        def _solver():
            return EMESolver(num_modes=[1], num_cells=[1], lengths=[1.0],
                             oxide_width=3.0, oxide_height=2.0, res=20.0,
                             material_stack=SOI_160, mode_solver="tidy3d")

        smat_a = _solver().solve(geometry, wavelengths, ports=["through"], save_dir=save_dir)
        smat_b = _solver().solve(geometry, wavelengths, ports=["through"], save_dir=save_dir)
        np.testing.assert_array_equal(
            smat_a.coeff(port="through"), smat_b.coeff(port="through"))

        run_dirs = [d for d in os.listdir(save_dir)
                    if os.path.isdir(os.path.join(save_dir, d))]
        assert len(run_dirs) == 1, (
            f"identical solver/geometry/wavelengths produced {len(run_dirs)} "
            "cache entries instead of one hit -- the fingerprint gained "
            "nondeterminism"
        )


@pytest.mark.slow
class TestSectionModes:
    """EMESolver.section_modes() -- the field-sampling counterpart
    solve_section() computes internally but never surfaces (see
    HybridSolver's check_mode_correspondence, beef/solvers/hybrid.py)."""

    def test_returns_sane_modes_on_its_own_mesh(self):
        geometry = Taper(length=2.0, width_start=0.45, width_end=0.45)  # straight
        solver = EMESolver(
            num_modes=[2], num_cells=[1], lengths=[2.0],
            oxide_width=3.0, oxide_height=1.5, res=20.0,
            material_stack=SOI_160, mode_solver="tidy3d",
        )
        modes, mesh = solver.section_modes(geometry, wavelength=1.55, x0=-1.0, x1=1.0, n_modes=2)

        assert len(modes) == 2
        neff = sorted((float(np.real(m.neff)) for m in modes), reverse=True)
        assert neff[0] > 1.0  # a 450nm SOI strip guides at 1.55 um
        assert neff[0] > neff[1]

        y, z = mesh_cell_centers(mesh)
        assert np.asarray(modes[0].Ex).shape == (len(y), len(z))

    def test_accepts_a_caller_supplied_shared_mesh(self):
        """The mesh param exists so a caller (HybridSolver's correspondence
        check) can hand in ITS OWN diagnostic mesh instead of this solver's
        production self._build_mesh() -- confirm the returned mesh really is
        the one passed in, not silently rebuilt."""
        geometry = Taper(length=2.0, width_start=0.45, width_end=0.45)
        solver = EMESolver(
            num_modes=[2], num_cells=[1], lengths=[2.0],
            oxide_width=3.0, oxide_height=1.5, res=20.0,
            material_stack=SOI_160, mode_solver="tidy3d",
        )
        own_mesh = solver._build_mesh()
        modes, returned_mesh = solver.section_modes(
            geometry, wavelength=1.55, x0=-1.0, x1=1.0, n_modes=1, mesh=own_mesh)
        assert returned_mesh is own_mesh
        assert len(modes) == 1
