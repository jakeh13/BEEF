"""Fast, meep/meow-free tests for HybridSolver's check_mode_correspondence
guard (validation_markdown/modal_multiplexer_validation.md gaps 1/6).

meow is faked out entirely (a minimal sys.modules stand-in exposing only
Mesh2D, mirroring tests/test_duplicate_band_guard.py's fake-meep pattern) so
this exercises the guard's own dispatch/threshold logic, not a real
eigensolve. Mode fields are synthetic sine-order profiles (the same toy
pattern tests/test_mode_identity_overlap.py uses to pin
field_overlap_matrix()/best_correspondence()'s behavior) so correspondence
correctness/mismatch is exact and grid-independent, not a numerical
coincidence.
"""
from __future__ import annotations

import sys
import types

import numpy as np
import pytest

from beef.solvers.hybrid import HybridSolver

WAVELENGTHS = np.array([1.55])


class _FakeMesh2D:
    def __init__(self, x, y):
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)


@pytest.fixture(autouse=True)
def fake_meow(monkeypatch):
    fake = types.ModuleType("meow")
    fake.Mesh2D = _FakeMesh2D
    monkeypatch.setitem(sys.modules, "meow", fake)
    yield


class FakeGeometry:
    def width_at(self, x, **kwargs):
        return 0.45

    def center_at(self, x, **kwargs):
        return 0.0


class FakeSolver:
    """Minimal solve_section()-only stand-in, matching tests/test_hybrid.py's
    own FakeSolver -- duplicated here (not imported) since that file's fakes
    are private to it, same convention as every other tests/test_*.py file
    in this repo."""

    def __init__(self, matrix):
        self.matrix = np.asarray(matrix, dtype=complex)

    def solve_section(self, geometry, wavelengths, x0, x1,
                       n_modes_in, n_modes_out, stub_left=0.0, stub_right=0.0,
                       **kwargs):
        M = len(wavelengths)
        return np.broadcast_to(self.matrix, (M,) + self.matrix.shape).copy()


def _sine_fields(n_modes: int, y_grid: np.ndarray, z_grid: np.ndarray, orders) -> dict:
    """Synthetic transverse-E field set: returned mode i has shape
    Ey = sin((orders[i]+1) pi u), u in [0, 1] over y_grid's own span -- so
    two modes sharing a sine order are identical (overlap 1) and two modes
    with different orders are orthogonal (overlap ~0), independent of the
    caller's actual grid extent."""
    y_grid = np.asarray(y_grid, dtype=float)
    ny, nz = len(y_grid), len(z_grid)
    u = (y_grid - y_grid.min()) / (y_grid.max() - y_grid.min())
    Ey = np.zeros((n_modes, ny, nz), dtype=complex)
    for i, order in enumerate(orders):
        Ey[i] = np.sin((order + 1) * np.pi * u)[:, None]
    zeros = np.zeros((n_modes, ny, nz), dtype=complex)
    return dict(Ex=zeros.copy(), Ey=Ey, Ez=zeros.copy(),
                Hx=zeros.copy(), Hy=zeros.copy(), Hz=zeros.copy())


class FakeFieldSolver(FakeSolver):
    """FakeSolver + a scripted section_mode_fields() returning the sine-order
    toy fields above. `orders` (default identity) lets a test simulate a
    genuine correspondence mismatch by requesting a shuffled order list."""

    def __init__(self, matrix, sim_width=2.0, sim_height=1.0, orders=None):
        super().__init__(matrix)
        self.sim_width = sim_width
        self.sim_height = sim_height
        self.orders = orders
        # Default n_eff are all well above any plausible cladding index, so
        # the guided-mode filter is a no-op unless a test overrides these.
        self.neff_values = None

    def section_mode_fields(self, geometry, wavelength, x0, x1, n_modes,
                             y_grid, z_grid, **kwargs):
        orders = list(range(n_modes)) if self.orders is None else self.orders
        fields = _sine_fields(n_modes, y_grid, z_grid, orders)
        if self.neff_values is None:
            neff = np.array([2.0 - 0.01 * k for k in range(n_modes)])
        else:
            neff = np.asarray(self.neff_values, dtype=float)
        return dict(neff=neff, **fields)


def _chain(solver_a, solver_b, n_interface=3, **kwargs):
    return HybridSolver(
        [solver_a, solver_b], boundaries=[0.0, 1.0, 2.0], interface_modes=[n_interface],
        check_mode_correspondence=True, **kwargs,
    )


class TestOptInDefault:
    def test_check_is_off_by_default(self):
        """No sim_width/section_mode_fields on either fake -- if the check
        ran, it would have nothing to work with. Confirms the default
        (check_mode_correspondence=False) never even attempts the check."""
        a, b = FakeSolver(np.eye(3)), FakeSolver(np.eye(3))
        solver = HybridSolver([a, b], boundaries=[0.0, 1.0, 2.0], interface_modes=[3])
        result = solver.solve(FakeGeometry(), WAVELENGTHS, in_modes=3, out_modes=3)
        assert result is not None


class TestMatchingBases:
    def test_identical_bases_pass_silently(self):
        a = FakeFieldSolver(np.eye(3))
        b = FakeFieldSolver(np.eye(3))
        solver = _chain(a, b, n_interface=3)
        solver.solve(FakeGeometry(), WAVELENGTHS, in_modes=1, out_modes=1)  # must not raise


class TestMismatchedBases:
    def test_permuted_basis_raises_by_default(self):
        a = FakeFieldSolver(np.eye(3), orders=[0, 1, 2])
        b = FakeFieldSolver(np.eye(3), orders=[1, 0, 2])  # modes 0/1 swapped
        solver = _chain(a, b, n_interface=3)
        with pytest.raises(RuntimeError, match="mode-correspondence check failed"):
            solver.solve(FakeGeometry(), WAVELENGTHS, in_modes=1, out_modes=1)

    def test_permuted_basis_only_warns_when_action_is_warn(self, caplog):
        import logging
        a = FakeFieldSolver(np.eye(3), orders=[0, 1, 2])
        b = FakeFieldSolver(np.eye(3), orders=[1, 0, 2])
        solver = _chain(a, b, n_interface=3, mode_correspondence_action="warn")
        with caplog.at_level(logging.WARNING, logger="beef.solvers.hybrid"):
            solver.solve(FakeGeometry(), WAVELENGTHS, in_modes=1, out_modes=1)
        assert any("mode-correspondence check failed" in r.message for r in caplog.records)

    def test_low_overlap_below_threshold_raises(self):
        """Same order everywhere (a trivially valid identity permutation)
        but one mode's shape is deliberately corrupted below the overlap
        threshold -- order_preserved alone must not be enough to pass."""
        a = FakeFieldSolver(np.eye(3), orders=[0, 1, 2])
        b = FakeFieldSolver(np.eye(3), orders=[0, 1, 2])
        orig = FakeFieldSolver.section_mode_fields

        def corrupted(self, geometry, wavelength, x0, x1, n_modes, y_grid, z_grid, **kwargs):
            result = orig(self, geometry, wavelength, x0, x1, n_modes, y_grid, z_grid, **kwargs)
            result["Ey"][1] = np.sin(9.7 * np.asarray(y_grid))[:, None]  # not a matching sine order
            return result

        b.section_mode_fields = corrupted.__get__(b, FakeFieldSolver)
        solver = _chain(a, b, n_interface=3, mode_correspondence_min_overlap=0.9)
        with pytest.raises(RuntimeError, match="mode-correspondence check failed"):
            solver.solve(FakeGeometry(), WAVELENGTHS, in_modes=1, out_modes=1)


class TestUnsupportedSolverIsSkipped:
    def test_skips_with_log_when_one_side_lacks_the_hook(self, caplog):
        import logging
        a = FakeFieldSolver(np.eye(3))
        b = FakeSolver(np.eye(3))  # no section_mode_fields, no sim_width
        solver = _chain(a, b, n_interface=3)
        with caplog.at_level(logging.INFO, logger="beef.solvers.hybrid"):
            solver.solve(FakeGeometry(), WAVELENGTHS, in_modes=1, out_modes=1)  # must not raise
        assert any("skipping the interface-0 correspondence check" in r.message
                   for r in caplog.records)


class _FakeStack:
    def __init__(self, n_clad):
        self._n = n_clad

    def sio2_index(self, wl):
        return self._n


class TestGuidedModeFilter:
    """Regression guards for a real false-alarm source: radiation modes
    (n_eff below the cladding index) are not discrete physical modes -- their
    shape is set by the box and each engine's own discretization, so two
    engines' radiation modes overlap ~0 however well the real basis
    corresponds. Measured on the modal multiplexer: including them drove the
    reported worst overlap to 0.0000 purely from bands 13-15, while the 13
    guided bands were unaffected. validation_markdown/README.md's "Finding 1"
    records that EVERY interface basis in this repo carries at least one
    unguided mode, so without this filter the check false-alarms everywhere.
    """

    def _solvers(self, orders_a, orders_b, neff_a, neff_b, n_clad=1.45):
        a = FakeFieldSolver(np.eye(len(orders_a)), orders=orders_a)
        b = FakeFieldSolver(np.eye(len(orders_b)), orders=orders_b)
        a.neff_values, b.neff_values = neff_a, neff_b
        a.material_stack = b.material_stack = _FakeStack(n_clad)
        return a, b

    def test_unguided_modes_are_excluded_from_the_comparison(self):
        """Modes 0-1 correspond; mode 2 is unguided on both sides and
        deliberately mismatched. It must not fail the check."""
        a, b = self._solvers([0, 1, 2], [0, 1, 5], [2.0, 1.8, 1.2], [2.0, 1.8, 1.2])
        solver = _chain(a, b, n_interface=3)
        solver.solve(FakeGeometry(), WAVELENGTHS, in_modes=1, out_modes=1)  # must not raise

    def test_a_guided_mode_mismatch_is_still_caught(self):
        """The filter must not become a blanket excuse: a genuine
        disagreement among GUIDED modes still has to fail."""
        a, b = self._solvers([0, 1, 2], [1, 0, 2], [2.0, 1.8, 1.6], [2.0, 1.8, 1.6])
        solver = _chain(a, b, n_interface=3)
        with pytest.raises(RuntimeError, match="mode-correspondence check failed"):
            solver.solve(FakeGeometry(), WAVELENGTHS, in_modes=1, out_modes=1)

    def test_filter_is_skipped_when_no_material_stack_is_available(self):
        """Without a stack there is no cladding index to filter on; the
        check must still run rather than guess a threshold."""
        a = FakeFieldSolver(np.eye(2), orders=[0, 1])
        b = FakeFieldSolver(np.eye(2), orders=[0, 1])
        assert HybridSolver._cladding_index(a, b, 1.55) is None
        solver = _chain(a, b, n_interface=2)
        solver.solve(FakeGeometry(), WAVELENGTHS, in_modes=1, out_modes=1)  # must not raise

    def test_cladding_index_is_read_from_either_solver(self):
        a = FakeFieldSolver(np.eye(1))
        b = FakeFieldSolver(np.eye(1))
        b.material_stack = _FakeStack(1.4468)
        assert HybridSolver._cladding_index(a, b, 1.31) == pytest.approx(1.4468)


class TestDiagnosticMeshResolution:
    """Regression guards for a real bug: the diagnostic mesh was originally a
    fixed 24-pixel-wide grid, which works out to only ~4.8 px/um on a 5 um
    domain -- leaving SOI_160's 0.16 um core 0.77 cells tall. Since an EME
    side SOLVES its modes on this mesh, that made one side's modes physically
    wrong, not just coarsely sampled: measured 0.0004 worst overlap with
    modes 7-15 scrambled on the modal multiplexer, against 0.973-0.998 and
    clean ordering for the same interface at 40 px/um."""

    def test_mesh_pixels_are_exactly_square(self):
        """meow's compute_modes_meep() rejects a non-square mesh outright, so
        this must hold for any width/height aspect ratio."""
        mesh = HybridSolver._interface_diagnostic_mesh(5.0, 1.5, res=40.0)
        dx = mesh.x[1] - mesh.x[0]
        dy = mesh.y[1] - mesh.y[0]
        assert dx == pytest.approx(dy)
        assert dx == pytest.approx(1.0 / 40.0)

    def test_mesh_pixels_square_for_awkward_aspect_ratio(self):
        mesh = HybridSolver._interface_diagnostic_mesh(3.7, 1.15, res=37.0)
        assert (mesh.x[1] - mesh.x[0]) == pytest.approx(mesh.y[1] - mesh.y[0])

    def test_mesh_actually_resolves_a_thin_core(self):
        """The specific failure: a 0.16 um SOI_160 core must land on several
        cells, not a fraction of one."""
        mesh = HybridSolver._interface_diagnostic_mesh(5.0, 1.5, res=40.0)
        dz = mesh.y[1] - mesh.y[0]
        assert 0.16 / dz >= 4.0

    def test_res_defaults_to_the_coarser_of_the_two_solvers(self, solver_pair_res):
        a, b = solver_pair_res
        chain = _chain(a, b, n_interface=1)
        assert chain._interface_diagnostic_res(a, b) == pytest.approx(40.0)

    def test_explicit_res_overrides_the_derived_one(self, solver_pair_res):
        a, b = solver_pair_res
        chain = _chain(a, b, n_interface=1, mode_correspondence_res=25.0)
        assert chain._interface_diagnostic_res(a, b) == pytest.approx(25.0)

    def test_falls_back_when_neither_solver_exposes_res(self):
        a, b = FakeFieldSolver(np.eye(1)), FakeFieldSolver(np.eye(1))
        chain = _chain(a, b, n_interface=1)
        # FakeFieldSolver has no `res` attribute at all.
        assert chain._interface_diagnostic_res(a, b) == pytest.approx(40.0)


@pytest.fixture
def solver_pair_res():
    a = FakeFieldSolver(np.eye(1))
    b = FakeFieldSolver(np.eye(1))
    a.res, b.res = 80.0, 40.0  # coarser of the two must win
    return a, b


class TestModeCountMismatch:
    def test_raises_if_a_side_returns_the_wrong_mode_count(self):
        a = FakeFieldSolver(np.eye(3))
        b = FakeFieldSolver(np.eye(3))
        orig = FakeFieldSolver.section_mode_fields

        def short(self, geometry, wavelength, x0, x1, n_modes, y_grid, z_grid, **kwargs):
            result = orig(self, geometry, wavelength, x0, x1, n_modes, y_grid, z_grid, **kwargs)
            return {k: (v[:-1] if k != "neff" else v[:-1]) for k, v in result.items()}

        b.section_mode_fields = short.__get__(b, FakeFieldSolver)
        solver = _chain(a, b, n_interface=3)
        with pytest.raises(RuntimeError, match="needed 3 comparable modes"):
            solver.solve(FakeGeometry(), WAVELENGTHS, in_modes=1, out_modes=1)
