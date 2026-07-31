"""Guards beef.mode_identity's overlap+assignment math itself, with
synthetic fields -- fast, no meep/meow, no simulation.

Mirrors tests/test_mode_identity.py's TestModeCorrespondenceLogic pattern
(guard the checker before trusting what it says about real solves): if this
cannot detect a shuffle or distinguish orthogonal mode shapes, it would pass
forever and mean nothing.
"""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from beef.mode_identity import (
    overlap_diagnostics, self_overlap_diagnostics,
    best_correspondence, field_overlap_matrix, mesh_cell_centers, stack_meow_mode_fields,
)

N_Y = 200


def _toy_field_set(n_modes: int, y: np.ndarray) -> dict:
    """Toy transverse E fields on a 1-D-in-y grid (single z point): mode k
    is Ey = sin((k+1) pi y), Ez = 0. field_overlap_matrix() only reads
    Ey/Ez, so self-overlap is integral(sin^2) > 0 and cross-overlap between
    different k is the standard sine orthogonality integral (~0 on a fine
    enough grid) -- exactly the "different physical mode" case this module
    must distinguish, without any n_eff involved at all. Ex/Hx/Hy/Hz are
    present (matching section_mode_fields()'s real return shape) but
    unused zeros -- see field_overlap_matrix()'s docstring for why H isn't
    needed at all.
    """
    ny = len(y)
    zeros = np.zeros((n_modes, ny, 1), dtype=complex)
    Ey = np.stack([np.sin((k + 1) * np.pi * y) for k in range(n_modes)])[:, :, None].astype(complex)
    return dict(Ex=zeros.copy(), Ey=Ey, Ez=zeros.copy(), Hx=zeros.copy(), Hy=zeros.copy(), Hz=zeros.copy())


@pytest.fixture
def y_grid():
    # Cell-center-style midpoints on [0, 1], matching mesh_cell_centers()'s
    # own convention, so dy below is exactly right for the Riemann sum.
    edges = np.linspace(0.0, 1.0, N_Y + 1)
    return 0.5 * (edges[:-1] + edges[1:])


class TestFieldOverlapMatrix:
    def test_identical_fields_are_identity_and_unit_diagonal(self, y_grid):
        fields = _toy_field_set(4, y_grid)
        dy = 1.0 / N_Y
        eta = field_overlap_matrix(fields, fields, dy=dy, dz=1.0)
        np.testing.assert_allclose(np.diag(eta), 1.0, atol=1e-9)
        # Off-diagonal: different sine orders are orthogonal on this grid.
        off_diag = eta - np.diag(np.diag(eta))
        assert np.max(np.abs(off_diag)) < 1e-3

    def test_orthogonal_modes_have_near_zero_cross_overlap(self, y_grid):
        fields = _toy_field_set(3, y_grid)
        dy = 1.0 / N_Y
        eta = field_overlap_matrix(fields, fields, dy=dy, dz=1.0)
        assert eta[0, 1] < 1e-3
        assert eta[0, 2] < 1e-3
        assert eta[1, 2] < 1e-3

    def test_normalization_is_independent_of_each_engine_s_own_scale(self, y_grid):
        """eta must not depend on either side's absolute field amplitude --
        meep and meow are not required to agree on mode-power normalization,
        only on shape (see module docstring)."""
        fields_a = _toy_field_set(2, y_grid)
        fields_b = _toy_field_set(2, y_grid)
        for name in fields_b:
            fields_b[name] = fields_b[name] * 37.0  # arbitrary rescale
        dy = 1.0 / N_Y
        eta = field_overlap_matrix(fields_a, fields_b, dy=dy, dz=1.0)
        np.testing.assert_allclose(np.diag(eta), 1.0, atol=1e-9)

    def test_immune_to_mismatched_cross_engine_e_to_h_scaling(self, y_grid):
        """Regression test for a real bug: an earlier version of this
        module used the full E-x-H reciprocity overlap, mixing engine A's E
        with engine B's H. meep and meow were found to report E and H in
        mutually inconsistent absolute units (~375x apart, confirmed via a
        real solved mode) despite each being internally self-consistent --
        that drove the true-match diagonal down to ~0.05 across every
        device/resolution tried, silently, since the shape match was
        genuinely correct and only the cross-engine impedance normalization
        was wrong. Fields_b here has an Hz wildly rescaled relative to
        fields_a's Hz (mimicking that mismatch) while its Ey/Ez shape is
        identical to fields_a's -- eta must be unaffected, since Hz/Hy are
        not read by field_overlap_matrix() at all."""
        fields_a = _toy_field_set(3, y_grid)
        fields_b = _toy_field_set(3, y_grid)
        fields_b["Hz"] = fields_b["Hz"] + 999.0  # would have poisoned the old E-x-H* formula
        fields_b["Hy"] = fields_b["Hy"] - 12345.0
        dy = 1.0 / N_Y
        eta = field_overlap_matrix(fields_a, fields_b, dy=dy, dz=1.0)
        np.testing.assert_allclose(np.diag(eta), 1.0, atol=1e-9)


class TestBestCorrespondence:
    def test_recovers_identity_when_nothing_is_shuffled(self, y_grid):
        fields = _toy_field_set(5, y_grid)
        eta = field_overlap_matrix(fields, fields, dy=1.0 / N_Y, dz=1.0)
        correspondence, margin = best_correspondence(eta)
        assert correspondence == list(range(5))
        assert np.all(margin > 0.9)  # unambiguous: best ~1, next-best ~0

    def test_recovers_a_shuffle_that_nearest_neighbour_style_matching_would_need_to_survive(self, y_grid):
        """The whole point of this module: matching by best total overlap,
        not by nearest index/scalar. Shuffle engine B's mode order and
        confirm the assignment un-shuffles it exactly."""
        n = 5
        fields_a = _toy_field_set(n, y_grid)
        perm = [2, 0, 4, 1, 3]  # fields_b[j] holds fields_a's mode perm[j]
        fields_b = {name: arr[perm] for name, arr in fields_a.items()}

        eta = field_overlap_matrix(fields_a, fields_b, dy=1.0 / N_Y, dz=1.0)
        correspondence, margin = best_correspondence(eta)

        inv_perm = list(np.argsort(perm))
        assert correspondence == inv_perm
        assert np.all(margin > 0.9)

    def test_margin_is_small_for_a_genuine_degeneracy(self, y_grid):
        """Two IDENTICAL mode shapes at different indices (a true physical
        degeneracy, not a solver artifact) must show a small margin -- the
        assignment still has to pick one, but the low margin is the signal
        that this pairing is ambiguous rather than confidently resolved."""
        fields_a = _toy_field_set(3, y_grid)
        fields_b = _toy_field_set(3, y_grid)
        for name in fields_b:
            arr = fields_b[name].copy()
            arr[1] = fields_b[name][0]  # mode 1 is now a copy of mode 0
            fields_b[name] = arr

        eta = field_overlap_matrix(fields_a, fields_b, dy=1.0 / N_Y, dz=1.0)
        _, margin = best_correspondence(eta)
        assert margin[0] < 0.1

    def test_detects_a_missing_or_extra_mode_via_rectangular_eta(self, y_grid):
        fields_a = _toy_field_set(4, y_grid)
        fields_b = {name: arr[:3] for name, arr in fields_a.items()}
        eta = field_overlap_matrix(fields_a, fields_b, dy=1.0 / N_Y, dz=1.0)
        assert eta.shape == (4, 3)
        correspondence, _ = best_correspondence(eta)
        assert len(correspondence) == 4
        matched = [c for c in correspondence if c != -1]
        assert sorted(matched) == [0, 1, 2]  # every B column used exactly once
        assert correspondence.count(-1) == 1  # exactly one A row has no B counterpart


class TestOverlapDiagnostics:
    """The layer-1 upgrade from "test the diagonal" to "test the matrix"
    (validation_markdown/fdtd_vs_EME_mode_validation.md). If two bases span
    the same subspace and correspond, eta is close to a PERMUTATION matrix;
    the matched diagonal alone tests a strict subset of that and is blind to
    a mode present in one basis and missing from the other."""

    def test_identical_bases_are_a_clean_permutation(self, y_grid):
        fields = _toy_field_set(4, y_grid)
        eta = field_overlap_matrix(fields, fields, dy=1.0 / N_Y, dz=1.0)
        d = overlap_diagnostics(eta)
        assert d["order_preserved"] is True
        assert d["worst_diagonal"] > 0.999
        assert d["worst_row_completeness"] > 0.999
        assert d["worst_col_completeness"] > 0.999
        assert np.nanmax(d["off_diagonal_max"]) < 1e-3

    def test_a_pure_shuffle_keeps_completeness_at_one(self, y_grid):
        """A mis-ORDERING is not a basis deficiency: every mode is still
        present and representable, so completeness must stay ~1 even though
        the correspondence is not the identity. This is the distinction the
        diagonal alone cannot draw."""
        n = 4
        fields_a = _toy_field_set(n, y_grid)
        perm = [2, 0, 3, 1]
        fields_b = {name: arr[perm] for name, arr in fields_a.items()}
        d = overlap_diagnostics(field_overlap_matrix(fields_a, fields_b, dy=1.0 / N_Y, dz=1.0))
        assert d["order_preserved"] is False           # ordering is wrong ...
        assert d["worst_diagonal"] > 0.999             # ... but every mode found a partner
        assert d["worst_row_completeness"] > 0.999     # ... and nothing is missing
        assert d["worst_col_completeness"] > 0.999

    def test_a_missing_mode_collapses_row_completeness(self, y_grid):
        """The failure the diagonal-only check cannot see. Engine B's basis is
        missing engine A's mode 2 entirely (it carries a duplicate of mode 0
        in its place). A's row 2 then has no partner anywhere in B's span, so
        its completeness collapses toward 0 while OTHER rows still look
        perfect -- exactly the signature of truncation / basis mismatch, as
        opposed to a mis-ordering."""
        n = 4
        fields_a = _toy_field_set(n, y_grid)
        fields_b = {name: arr.copy() for name, arr in fields_a.items()}
        for name in fields_b:
            fields_b[name][2] = fields_a[name][0]  # mode 2 replaced by a copy of mode 0
        d = overlap_diagnostics(field_overlap_matrix(fields_a, fields_b, dy=1.0 / N_Y, dz=1.0))
        assert d["row_completeness"][2] < 0.05   # A's mode 2 is unrepresentable in B
        assert d["row_completeness"][1] > 0.99   # while its neighbours are fine
        assert d["row_completeness"][3] > 0.99
        assert d["worst_row_completeness"] < 0.05

    def test_rectangular_eta_flags_the_unmatched_row(self, y_grid):
        fields_a = _toy_field_set(4, y_grid)
        fields_b = {name: arr[:3] for name, arr in fields_a.items()}
        d = overlap_diagnostics(field_overlap_matrix(fields_a, fields_b, dy=1.0 / N_Y, dz=1.0))
        assert d["correspondence"].count(-1) == 1
        # The dropped mode has no column to be represented by at all.
        assert d["worst_row_completeness"] < 0.05


class TestSelfOverlapDiagnostics:
    """Layer 4 -- the same-engine control. One engine's modes against
    themselves MUST be the identity by construction; anything else is a fault
    in the measurement, not a finding about the device."""

    def test_a_well_resolved_set_is_the_identity(self, y_grid):
        fields = _toy_field_set(5, y_grid)
        d = self_overlap_diagnostics(fields, dy=1.0 / N_Y, dz=1.0)
        assert d["is_identity"] is True
        assert d["max_off_diagonal"] < 0.1

    def test_a_duplicated_mode_breaks_the_control(self, y_grid):
        """A duplicate eigenvector in one engine's own set -- the real
        FDTDSolver failure _iter_distinct_bands() guards -- must show up here
        as a non-identity self-overlap, before any cross-engine comparison
        gets a chance to blame the other engine."""
        fields = _toy_field_set(3, y_grid)
        for name in fields:
            fields[name][1] = fields[name][0]
        d = self_overlap_diagnostics(fields, dy=1.0 / N_Y, dz=1.0)
        assert d["is_identity"] is False
        assert d["max_off_diagonal"] > 0.9  # modes 0 and 1 are the same vector

    def test_a_non_orthogonal_mode_set_breaks_the_control(self, y_grid):
        """The other real failure mode, stated precisely. An under-resolved
        EME *solve* does not return coarsely-sampled eigenmodes -- it returns
        vectors that are not eigenmodes at all, and therefore not mutually
        orthogonal. Mode 1 here is a superposition rather than a pure order,
        mimicking that, and the control catches it.

        Note what this test deliberately does NOT claim: coarse *sampling* of
        an already-orthogonal set is not itself enough to fail. Verified while
        writing this -- three midpoint samples of sin(k*pi*y) preserve the
        orthogonality integral exactly (a quadrature coincidence), so
        `is_identity` still passes there. So the control's power against the
        fixed-pixel-count mesh bug comes from the solver's output degrading,
        not from the comparison grid being coarse per se.
        """
        fields = _toy_field_set(3, y_grid)
        for name in ("Ey",):
            fields[name][1] = fields[name][1] + 0.8 * fields[name][0]
        d = self_overlap_diagnostics(fields, dy=1.0 / N_Y, dz=1.0)
        assert d["is_identity"] is False
        assert d["max_off_diagonal"] > 0.1

    def test_coarse_sampling_of_orthogonal_modes_still_passes(self):
        """Pins the negative result above so the control's guarantee is not
        overstated in future: this grid is far too coarse to resolve anything,
        yet the control passes, because the *set* is still orthogonal."""
        edges = np.linspace(0.0, 1.0, 4)
        coarse = 0.5 * (edges[:-1] + edges[1:])
        d = self_overlap_diagnostics(_toy_field_set(3, coarse), dy=1.0 / 3, dz=1.0)
        assert d["is_identity"] is True


class TestMeshCellCenters:
    def test_returns_midpoints_not_nodes(self):
        mesh = SimpleNamespace(x=np.array([0.0, 1.0, 2.0, 3.0]), y=np.array([0.0, 2.0, 4.0]))
        y, z = mesh_cell_centers(mesh)
        np.testing.assert_allclose(y, [0.5, 1.5, 2.5])
        np.testing.assert_allclose(z, [1.0, 3.0])

    def test_output_length_is_one_less_than_input(self):
        mesh = SimpleNamespace(x=np.linspace(-1, 1, 31), y=np.linspace(-0.5, 0.5, 16))
        y, z = mesh_cell_centers(mesh)
        assert len(y) == 30
        assert len(z) == 15


class TestStackMeowModeFields:
    def test_stacks_each_component_across_modes(self):
        mode_a = SimpleNamespace(Ex=np.ones((2, 3)), Ey=np.zeros((2, 3)), Ez=np.zeros((2, 3)),
                                  Hx=np.zeros((2, 3)), Hy=np.zeros((2, 3)), Hz=np.zeros((2, 3)))
        mode_b = SimpleNamespace(Ex=np.full((2, 3), 5.0), Ey=np.zeros((2, 3)), Ez=np.zeros((2, 3)),
                                  Hx=np.zeros((2, 3)), Hy=np.zeros((2, 3)), Hz=np.zeros((2, 3)))
        stacked = stack_meow_mode_fields([mode_a, mode_b])
        assert set(stacked) == {"Ex", "Ey", "Ez", "Hx", "Hy", "Hz"}
        # mode_a/mode_b's nonzero component is meow's Ex (device-frame Ey
        # after the remap below) -- shape/stacking check only, remap
        # direction is covered separately.
        assert stacked["Ey"].shape == (2, 2, 3)
        np.testing.assert_allclose(stacked["Ey"][0], 1.0)
        np.testing.assert_allclose(stacked["Ey"][1], 5.0)

    def test_remaps_meow_s_local_axes_into_device_frame(self):
        """meow's z is the propagation axis (device x); x, y are transverse
        (device y, z) -- confirmed directly against a real solved mode via
        sum(Px)==sum(Py)==0, sum(Pz)!=0 (see mode_identity.py's docstring).
        A mock mode with a distinct value per component pins the exact
        remap so a future edit can't silently swap two axes back."""
        mode = SimpleNamespace(
            Ex=np.full((1, 1), 1.0), Ey=np.full((1, 1), 2.0), Ez=np.full((1, 1), 3.0),
            Hx=np.full((1, 1), 4.0), Hy=np.full((1, 1), 5.0), Hz=np.full((1, 1), 6.0),
        )
        stacked = stack_meow_mode_fields([mode])
        # device Ex (propagation) = meow Ez; device Ey = meow Ex; device Ez = meow Ey.
        assert stacked["Ex"][0, 0, 0] == 3.0
        assert stacked["Ey"][0, 0, 0] == 1.0
        assert stacked["Ez"][0, 0, 0] == 2.0
        # Same pattern for H.
        assert stacked["Hx"][0, 0, 0] == 6.0
        assert stacked["Hy"][0, 0, 0] == 4.0
        assert stacked["Hz"][0, 0, 0] == 5.0
