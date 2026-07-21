"""Fast unit tests for the mode-count superset-reuse mechanism (Solver
recommendation #2): beef.cache.store.lookup_superset()/update_family_pointer()
plus the FDTDSolver/EMESolver family-config construction that feeds them.

No meep/meow import anywhere: FDTDSolver()/EMESolver() are constructed
(cheap) but solve_section()/solve_section_length_sweep() are never called
(both do `import meep`/`import meow` unconditionally, even down a
cache_only=True/superset-hit path — see test_cache_config.py's module
docstring for why that would break CI). Family-config construction is
instead exercised via section_fdtd_layout() + build_run_config() directly,
the same two calls solve_section() itself makes before touching meep.
"""
from __future__ import annotations

import json
import os

import numpy as np

from beef.cache.config import build_run_config
from beef.cache.store import (
    _family_pointer_path, load_result, lookup_superset, save_run, update_family_pointer,
)
from beef.geometry.taper import Taper
from beef.solvers.eme import EMESolver
from beef.solvers.fdtd import FDTDSolver, section_fdtd_layout


def _fdtd_family_config(geometry, solver, monitor_center, monitor_size, n_modes_in, n_modes_out):
    layout = section_fdtd_layout(
        geometry, -2.0, 2.0, 0.0, 0.0,
        n_modes_in, n_modes_out, solver.sim_width,
        solver.source_margin, solver.monitor_margin, solver.half_domain_aperture,
        {}, monitor_center=monitor_center, monitor_size=monitor_size,
    )
    return build_run_config(
        solver, "solve_section",
        dict(x0=-2.0, x1=2.0, wavelengths=[1.55],
             stub_left=0.0, stub_right=0.0,
             monitor_center=monitor_center, monitor_size=monitor_size,
             geometry_kwargs={}, reverse=layout.reverse),
        geometry,
    )


def _eme_family_config(geometry, solver):
    return build_run_config(
        solver, "solve_section",
        dict(x0=-2.0, x1=2.0, wavelengths=[1.55], stub_left=0.0, stub_right=0.0,
             geometry_kwargs={}),
        geometry,
    )


class TestStorePrimitives:
    """Pure beef.cache.store primitives — no solver classes needed."""

    def test_no_pointer_is_a_miss(self, tmp_path):
        assert lookup_superset(str(tmp_path), {"a": 1}, (6, 1)) is None

    def test_dominating_pointer_is_a_hit(self, tmp_path):
        save_dir = str(tmp_path)
        folder = save_run(save_dir, {"a": 1, "n": 12}, np.zeros((3, 12, 1)))
        update_family_pointer(save_dir, {"a": 1}, (12, 1), os.path.basename(folder))

        hit = lookup_superset(save_dir, {"a": 1}, (6, 1))
        assert hit is not None
        got_folder, stored_dims = hit
        assert got_folder == folder
        assert stored_dims == (12, 1)

    def test_equal_dims_is_a_hit(self, tmp_path):
        save_dir = str(tmp_path)
        folder = save_run(save_dir, {"a": 1}, np.zeros((3, 6, 1)))
        update_family_pointer(save_dir, {"a": 1}, (6, 1), os.path.basename(folder))
        assert lookup_superset(save_dir, {"a": 1}, (6, 1)) is not None

    def test_non_dominating_pointer_is_a_miss(self, tmp_path):
        save_dir = str(tmp_path)
        folder = save_run(save_dir, {"a": 1}, np.zeros((3, 6, 1)))
        update_family_pointer(save_dir, {"a": 1}, (6, 1), os.path.basename(folder))
        assert lookup_superset(save_dir, {"a": 1}, (12, 1)) is None  # exceeds stored "out"
        assert lookup_superset(save_dir, {"a": 1}, (6, 4)) is None   # exceeds stored "in"

    def test_missing_result_folder_is_a_miss(self, tmp_path):
        save_dir = str(tmp_path)
        update_family_pointer(save_dir, {"a": 1}, (12, 1), "nonexistent_hash")
        assert lookup_superset(save_dir, {"a": 1}, (6, 1)) is None

    def test_dims_length_mismatch_is_a_miss(self, tmp_path):
        save_dir = str(tmp_path)
        folder = save_run(save_dir, {"a": 1}, np.zeros((3, 6, 1)))
        update_family_pointer(save_dir, {"a": 1}, (6, 1), os.path.basename(folder))
        assert lookup_superset(save_dir, {"a": 1}, (6,)) is None

    def test_update_creates_pointer_when_none_exists(self, tmp_path):
        save_dir = str(tmp_path)
        update_family_pointer(save_dir, {"a": 1}, (2, 1), "hash_a")
        # Not yet backed by a real result folder, so lookup still misses on
        # the "does the folder exist" check — but the pointer file itself
        # must now exist and be valid JSON.
        path = _family_pointer_path(save_dir, {"a": 1})
        assert os.path.isfile(path)
        with open(path) as f:
            data = json.load(f)
        assert data == {"dims": [2, 1], "result_hash": "hash_a"}

    def test_update_overwrites_when_new_dominates(self, tmp_path):
        save_dir = str(tmp_path)
        f1 = save_run(save_dir, {"a": 1}, np.zeros((3, 2, 1)))
        f2 = save_run(save_dir, {"a": 1, "x": 2}, np.zeros((3, 6, 1)))
        update_family_pointer(save_dir, {"a": 1}, (2, 1), os.path.basename(f1))
        update_family_pointer(save_dir, {"a": 1}, (6, 1), os.path.basename(f2))

        hit = lookup_superset(save_dir, {"a": 1}, (4, 1))
        assert hit is not None
        assert hit[0] == f2
        assert hit[1] == (6, 1)

    def test_update_does_not_overwrite_when_new_does_not_dominate(self, tmp_path):
        save_dir = str(tmp_path)
        f_big_out = save_run(save_dir, {"a": 1}, np.zeros((3, 12, 1)))
        f_big_in = save_run(save_dir, {"a": 1, "x": 2}, np.zeros((3, 1, 6)))
        update_family_pointer(save_dir, {"a": 1}, (12, 1), os.path.basename(f_big_out))
        # (1, 6) doesn't dominate (12, 1) [1 < 12] and isn't dominated by it
        # either [6 > 1] -- incomparable, so the pointer must stay put.
        update_family_pointer(save_dir, {"a": 1}, (1, 6), os.path.basename(f_big_in))

        hit = lookup_superset(save_dir, {"a": 1}, (10, 1))
        assert hit is not None
        assert hit[0] == f_big_out
        assert hit[1] == (12, 1)


class TestSlicingRoundTrip:
    """The full save -> point -> lookup -> load -> slice path, with a real
    (synthetic) result array, no solver classes needed."""

    def test_sliced_result_matches_the_stored_array(self, tmp_path):
        save_dir = str(tmp_path)
        T = np.arange(3 * 12 * 4, dtype=complex).reshape(3, 12, 4)
        folder = save_run(save_dir, {"a": 1}, T)
        update_family_pointer(save_dir, {"a": 1}, (12, 4), os.path.basename(folder))

        hit = lookup_superset(save_dir, {"a": 1}, (6, 2))
        assert hit is not None
        got_folder, _ = hit
        sliced = load_result(got_folder)[:, :6, :2]
        np.testing.assert_array_equal(sliced, T[:, :6, :2])


class TestFDTDFamilyConfigConsistency:
    """FDTDSolver's family-config construction (mirrored from
    solve_section()'s own logic) must group requests that only differ in
    mode count, and must NOT group requests that cross the reverse
    (excited/measured side) boundary."""

    def test_same_family_when_only_n_out_grows_and_reverse_is_unchanged(self):
        geometry = Taper(length=4.0, width_start=0.45, width_end=0.45)
        solver = FDTDSolver(res=20.0, sim_width=4.0)

        cfg_a = _fdtd_family_config(geometry, solver, 0.0, 4.0, n_modes_in=1, n_modes_out=2)
        cfg_b = _fdtd_family_config(geometry, solver, 0.0, 4.0, n_modes_in=1, n_modes_out=12)
        assert cfg_a == cfg_b

    def test_different_family_when_reverse_flips(self):
        geometry = Taper(length=4.0, width_start=0.45, width_end=0.45)
        solver = FDTDSolver(res=20.0, sim_width=4.0)

        # n_out=1 -> reverse=True (n_out<=n_in); n_out=2 -> reverse=False.
        cfg_reverse_true = _fdtd_family_config(geometry, solver, 0.0, 4.0, n_modes_in=1, n_modes_out=1)
        cfg_reverse_false = _fdtd_family_config(geometry, solver, 0.0, 4.0, n_modes_in=1, n_modes_out=2)
        assert cfg_reverse_true != cfg_reverse_false

    def test_different_family_for_different_monitor_center(self):
        geometry = Taper(length=4.0, width_start=0.45, width_end=0.45)
        solver = FDTDSolver(res=20.0, sim_width=4.0)

        cfg_a = _fdtd_family_config(geometry, solver, 0.0, 4.0, n_modes_in=1, n_modes_out=2)
        cfg_b = _fdtd_family_config(geometry, solver, 1.0, 4.0, n_modes_in=1, n_modes_out=2)
        assert cfg_a != cfg_b


class TestEMEFamilyConfigConsistency:
    """EMESolver has no excited/measured-side asymmetry, so its family
    config should be identical across any mode-count split."""

    def test_same_family_across_any_mode_split(self):
        geometry = Taper(length=4.0, width_start=0.45, width_end=0.45)
        solver = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0])

        cfg_a = _eme_family_config(geometry, solver)
        cfg_b = _eme_family_config(geometry, solver)
        assert cfg_a == cfg_b

    def test_different_family_for_different_oxide_width(self):
        geometry = Taper(length=4.0, width_start=0.45, width_end=0.45)
        solver_a = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0], oxide_width=5.0)
        solver_b = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0], oxide_width=6.0)

        assert _eme_family_config(geometry, solver_a) != _eme_family_config(geometry, solver_b)
