"""Fast unit tests for beef.cache.store — no solver/meep needed."""
from __future__ import annotations

import json
import os

import numpy as np

from beef.cache.store import _hash_config, load_result, lookup, save_run


def _config(**call_kwargs_overrides):
    call_kwargs = dict(
        x0=-24.0, x1=-15.0, n_modes_in=1, n_modes_out=4,
        stub_left=0.5, stub_right=0.0,
        wavelengths=[1.50, 1.51, 1.52, 1.53, 1.54, 1.55],
        geometry_kwargs={"input_idx": 0},
    )
    call_kwargs.update(call_kwargs_overrides)
    return {
        "solver_type": "FDTDSolver",
        "solver_config": {"res": 40.0},
        "method": "solve_section",
        "call_kwargs": call_kwargs,
        "geometry_class": "DirectionalCoupler",
        "geometry_polygons": {"2,0": [[[0, 0], [1, 0], [1, 1]]]},
    }


class TestSaveRun:
    def test_round_trips_through_lookup(self, tmp_path):
        save_dir = str(tmp_path)
        config = _config()
        result = np.array([1.0 + 2.0j])
        save_run(save_dir, config, result)

        hit = lookup(save_dir, config)
        assert hit is not None
        np.testing.assert_allclose(load_result(hit), result)

    def test_writes_config_result_and_geometry_files(self, tmp_path):
        save_dir = str(tmp_path)
        config = _config()
        folder = save_run(save_dir, config, np.array([0.0]))

        assert (tmp_path / os.path.basename(folder) / "config.json").is_file()
        assert (tmp_path / os.path.basename(folder) / "result.pkl").is_file()
        assert _hash_config(config) in os.path.basename(folder)
        assert os.path.basename(folder).startswith("FDTDSolver_solve_section_DirectionalCoupler_")


class TestManifest:
    def test_appends_one_line_per_save(self, tmp_path):
        save_dir = str(tmp_path)
        save_run(save_dir, _config(n_modes_out=2), np.array([0.0]))
        save_run(save_dir, _config(n_modes_out=4), np.array([0.0]))

        lines = (tmp_path / "index.jsonl").read_text().splitlines()
        assert len(lines) == 2
        for line in lines:
            json.loads(line)  # each line is independently valid JSON

    def test_manifest_entry_summarizes_config(self, tmp_path):
        save_dir = str(tmp_path)
        config = _config(n_modes_out=6)
        save_run(save_dir, config, np.array([0.0]))

        entry = json.loads((tmp_path / "index.jsonl").read_text().splitlines()[0])
        assert entry["hash"] == _hash_config(config)
        assert entry["solver_type"] == "FDTDSolver"
        assert entry["method"] == "solve_section"
        assert entry["geometry_class"] == "DirectionalCoupler"
        assert entry["call_kwargs"]["n_modes_out"] == 6
        assert entry["call_kwargs"]["geometry_kwargs"] == {"input_idx": 0}

    def test_long_numeric_arrays_are_summarized_not_dumped(self, tmp_path):
        save_dir = str(tmp_path)
        config = _config(wavelengths=list(np.linspace(1.50, 1.60, 201)))
        save_run(save_dir, config, np.array([0.0]))

        entry = json.loads((tmp_path / "index.jsonl").read_text().splitlines()[0])
        wl_summary = entry["call_kwargs"]["wavelengths"]
        assert wl_summary == {"count": 201, "min": 1.50, "max": 1.60}

    def test_short_arrays_are_kept_as_is(self, tmp_path):
        save_dir = str(tmp_path)
        config = _config()  # wavelengths has 6 entries, at the >4 boundary
        save_run(save_dir, _config(wavelengths=[1.55, 1.56]), np.array([0.0]))

        entry = json.loads((tmp_path / "index.jsonl").read_text().splitlines()[0])
        assert entry["call_kwargs"]["wavelengths"] == [1.55, 1.56]

    def test_resaving_same_config_appends_duplicate_line(self, tmp_path):
        save_dir = str(tmp_path)
        config = _config()
        save_run(save_dir, config, np.array([0.0]))
        save_run(save_dir, config, np.array([1.0]))  # e.g. re-running the same script

        lines = (tmp_path / "index.jsonl").read_text().splitlines()
        assert len(lines) == 2
        assert json.loads(lines[0])["hash"] == json.loads(lines[1])["hash"]

    def test_manifest_failure_does_not_break_save_run(self, tmp_path, monkeypatch):
        import beef.cache.store as store

        def _boom(*args, **kwargs):
            raise RuntimeError("simulated manifest failure")

        monkeypatch.setattr(store, "_append_manifest", _boom)
        save_dir = str(tmp_path)
        config = _config()

        folder = save_run(save_dir, config, np.array([0.0]))  # must not raise

        assert (tmp_path / os.path.basename(folder) / "config.json").is_file()
        assert lookup(save_dir, config) == folder
        assert not (tmp_path / "index.jsonl").exists()


class TestSlugNaming:
    """Solver recommendation #3's folder-naming half: new writes get a
    human-readable slug prefix, and lookup() still finds cache entries
    written before this change (bare-hash folder name, no slug)."""

    def test_legacy_bare_hash_folder_is_still_found(self, tmp_path):
        # Simulate a cache entry written before slug-prefixed naming existed:
        # a bare-hash folder with no slug, holding config.json + result.pkl.
        save_dir = str(tmp_path)
        config = _config()
        legacy_folder = tmp_path / _hash_config(config)
        legacy_folder.mkdir()
        (legacy_folder / "config.json").write_text(json.dumps(config, sort_keys=True))
        import pickle
        with open(legacy_folder / "result.pkl", "wb") as f:
            pickle.dump(np.array([42.0]), f)

        hit = lookup(save_dir, config)
        assert hit == str(legacy_folder)
        np.testing.assert_allclose(load_result(hit), [42.0])

    def test_new_write_prefers_slug_folder_over_legacy(self, tmp_path):
        # If BOTH a legacy bare-hash folder and a new slug-prefixed folder
        # happen to exist for the same config, the slug-prefixed one (the
        # current naming scheme) wins.
        save_dir = str(tmp_path)
        config = _config()
        legacy_folder = tmp_path / _hash_config(config)
        legacy_folder.mkdir()
        (legacy_folder / "config.json").write_text(json.dumps(config, sort_keys=True))
        import pickle
        with open(legacy_folder / "result.pkl", "wb") as f:
            pickle.dump(np.array([1.0]), f)

        folder = save_run(save_dir, config, np.array([2.0]))
        assert folder != str(legacy_folder)

        hit = lookup(save_dir, config)
        assert hit == folder
        np.testing.assert_allclose(load_result(hit), [2.0])
