"""solve_section_length_sweep() persists the length-independent neff primitive,
never the derived transfer matrix.

The section this method builds is single-cell by construction
(`mw.create_cells(..., Ls=[1.0])`), so its propagation S-matrix is exactly
diagonal in the mode basis and `_propagate_diagonal()` rebuilds T at ANY length
from neff alone. Persisting T instead cost 342 MB for one 7601-length grid and,
because `lengths` is hashed verbatim into the cache key, made every other grid
a full miss.

These tests drive the `neff_override` path, which takes that analytic route
with no eigenmode solve at all — so they exercise the real caching behaviour
without importing meep or meow.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

from beef import EMESolver, SOI_160, Taper
from beef.solvers.eme import _propagate_diagonal

WL = np.array([1.55, 1.56])
# [n_wl, n_modes] -- two well-separated indices, no degeneracy to trip guards.
NEFF = np.array([[2.40 + 0j, 1.90 + 0j],
                 [2.39 + 0j, 1.89 + 0j]])


def _solver():
    return EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0],
                     oxide_width=3.0, oxide_height=2.0, res=20.0,
                     material_stack=SOI_160, mode_solver="tidy3d")


def _geometry():
    return Taper(length=1.0, width_start=0.8, width_end=0.8)   # straight


def _sweep(solver, lengths, save_dir):
    return solver.solve_section_length_sweep(
        # Taper is centred on x=0, so the section must span it whole:
        # a cut strictly inside a varying cross-section is refused.
        _geometry(), WL, -0.5, 0.5, np.asarray(lengths),
        n_modes_in=2, n_modes_out=2, neff_override=NEFF, save_dir=save_dir,
    )


def _entries(save_dir):
    return sorted(d for d in os.listdir(save_dir)
                  if os.path.isdir(os.path.join(save_dir, d))
                  and not d.startswith("_"))


class TestOnlyTheNeffPrimitiveIsPersisted:
    def test_a_fresh_sweep_writes_no_derived_entry(self, tmp_path):
        _sweep(_solver(), [5.0, 10.0, 20.0], str(tmp_path))
        entries = _entries(tmp_path)
        assert entries, "nothing was cached at all"
        assert not any("solve_section_length_sweep" in e for e in entries), (
            f"the derived transfer matrix was persisted: {entries}"
        )
        assert all("solve_section_neff" in e for e in entries), (
            f"expected only the length-independent neff primitive, got {entries}"
        )

    def test_a_different_length_grid_adds_no_entry(self, tmp_path):
        """The point of the change: `lengths` no longer keys anything."""
        solver = _solver()
        _sweep(solver, [5.0, 10.0], str(tmp_path))
        after_first = _entries(tmp_path)

        # A grid that shares no length with the first, and a different size.
        _sweep(solver, [7.5, 12.5, 33.0, 41.25], str(tmp_path))
        assert _entries(tmp_path) == after_first, (
            "a second length grid created a new cache entry -- lengths is still "
            "part of a cache key somewhere"
        )

    def test_reconstruction_matches_the_closed_form_it_claims(self, tmp_path):
        lengths = np.array([0.0, 3.0, 137.5])
        T = _sweep(_solver(), lengths, str(tmp_path))
        for li, L in enumerate(lengths):
            for wi, wl in enumerate(WL):
                for m in range(2):
                    expected = np.exp(2j * np.pi * NEFF[wi, m] / wl * L)
                    assert T[li, wi, m, m] == pytest.approx(expected, abs=0, rel=1e-12)
        off = T.copy()
        for m in range(2):
            off[:, :, m, m] = 0
        assert np.max(np.abs(off)) == 0.0, "a uniform section must be diagonal"

    def test_the_grids_agree_where_they_overlap(self, tmp_path):
        """Two calls, different grids, same physical lengths -> same numbers."""
        solver = _solver()
        a = _sweep(solver, [4.0, 9.0], str(tmp_path))
        b = _sweep(solver, [1.0, 4.0, 9.0, 16.0], str(tmp_path))
        np.testing.assert_array_equal(a[0], b[1])   # L = 4.0
        np.testing.assert_array_equal(a[1], b[2])   # L = 9.0


class TestPropagateDiagonalIsLengthIndependent:
    def test_any_length_grid_is_reconstructible_from_one_neff(self):
        fine = np.linspace(0.0, 400.0, 801)
        T_fine = _propagate_diagonal(NEFF, WL, fine, 2, 2)
        subset = fine[::7]
        T_sub = _propagate_diagonal(NEFF, WL, subset, 2, 2)
        np.testing.assert_array_equal(T_fine[::7], T_sub)


@pytest.mark.slow
class TestWholeDeviceSweepStoresTheFactorization:
    """solve_length_sweep() stores the cascade factorization, not the sweep.

    Real mode solves, so slow. The factorization is the interface S-matrices
    (a function of the modes alone) plus each cell's neff; it replays
    mw.compute_s_matrix() at any length grid, so a second, disjoint grid needs
    no new cache entry and no new mode solve.
    """

    def _solver(self):
        return EMESolver(num_modes=[2, 2, 2], num_cells=[1, 2, 1],
                         lengths=[0.5, 2.0, 0.5], oxide_width=3.0,
                         oxide_height=2.0, res=20.0, material_stack=SOI_160,
                         mode_solver="tidy3d")

    def _geom(self):
        return Taper(length=3.0, width_start=0.7, width_end=0.9)

    def test_one_entry_serves_two_disjoint_grids(self, tmp_path):
        solver, geom = self._solver(), self._geom()
        wl = np.array([1.55, 1.56])
        save_dir = str(tmp_path)

        a = solver.solve_length_sweep(geom, wl, np.array([1.0, 2.0, 3.0]),
                                      ports=["through"], save_dir=save_dir)
        entries = _entries(save_dir)
        assert entries == [e for e in entries if "factorization" in e], (
            f"expected only a factorization entry, got {entries}"
        )

        b = solver.solve_length_sweep(geom, wl, np.array([0.5, 1.0, 2.0, 3.0, 7.25]),
                                      ports=["through"], save_dir=save_dir)
        assert _entries(save_dir) == entries, (
            "a second length grid created a new entry -- the swept lengths are "
            "still part of a cache key"
        )

        # The lengths the two grids share must agree exactly, not approximately.
        for ia, ib in ((0, 1), (1, 2), (2, 3)):
            va = np.asarray(a[ia].transmission(port="through"))
            vb = np.asarray(b[ib].transmission(port="through"))
            np.testing.assert_array_equal(va, vb)

    def test_replay_reproduces_a_direct_meow_cascade(self):
        """The factorization's own claim, checked against meow directly."""
        import meow as mw
        from beef.cache.modes import extract_factorization, s_matrix_from_factorization

        solver, geom = self._solver(), self._geom()
        wl = 1.55
        cells = solver._build_cells(geom, side="through", wavelengths=np.array([wl]))
        env = mw.Environment(wl=wl, T=solver.temperature)
        css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
        modes = [solver.compute_modes(cs, nm)
                 for cs, nm in zip(css, solver._modes_per_cell())]
        fac = extract_factorization(modes)

        for Ls in ([0.5, 1.0, 1.0, 0.5], [0.5, 3.7, 3.7, 0.5], [0.5, 11.25, 11.25, 0.5]):
            S_ref, _ = mw.compute_s_matrix(modes, cells=cells, cell_lengths=Ls)
            S_fac, _ = s_matrix_from_factorization(fac, Ls)
            np.testing.assert_array_equal(np.asarray(S_ref), np.asarray(S_fac))
