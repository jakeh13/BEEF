"""Fast unit tests for the multi-rank safety of beef.cache.store's writes
(dc_validation.md item 4).

Every rank of an `srun -n 48` job runs the whole solve redundantly, so before
this every rank also wrote save_dir. A pointer file ended up empty under that
concurrency and raised out of json.load() — stub-sweep tasks 4 and 5 died that
way and then hung for their full 5 h wall clock. These tests pin the two
independent defences: only one rank writes at all, and what a damaged pointer
does is miss rather than raise.

No meep import: the rank gate is exercised by monkeypatching _is_master(), and
_mpi_rank()'s environment fallback (the EME-only case, where meep is never
imported yet srun still launches 48 ranks) by setting the launcher variables
directly.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pytest

from beef import _utils
from beef.cache import store
from beef.cache.store import (
    _LAUNCHER_RANK_VARS, _atomic_write_json, _family_pointer_path, _is_master,
    _read_family_pointer, lookup, lookup_superset, save_run, update_family_pointer,
)

_LAUNCHER_VARS = _LAUNCHER_RANK_VARS


@pytest.fixture
def serial_env(monkeypatch):
    """Strip launcher variables so a test that doesn't set one looks serial —
    otherwise running the suite itself under srun would leak a rank in.

    Also hide meep and mpi4py from `sys.modules`. `_is_master()` consults
    three sources in order of authority and these tests exercise the third
    (the launcher environment), so the first two have to be absent for the
    assertion to mean anything. Without this the tests pass or fail on
    *import order*: any earlier test that pulls in meow silently loads meep
    too — `import meow` alone puts it in `sys.modules` — and then
    `am_really_master()` answers instead, returning True in a serial pytest
    process and demoting these to vacuous. Pinning it here keeps the fallback
    genuinely covered no matter what else the suite has imported.
    """
    for var in _LAUNCHER_VARS:
        monkeypatch.delenv(var, raising=False)
    for mod in ("meep", "mpi4py.MPI"):
        monkeypatch.delitem(sys.modules, mod, raising=False)


def _config(**kwargs):
    base = {"solver_type": "FDTDSolver", "method": "solve_section",
            "geometry_class": "Taper", "call_kwargs": {"n_modes_out": 2}}
    base.update(kwargs)
    return base


class TestRankDetection:
    def test_no_mpi_is_master(self, serial_env):
        assert _is_master()

    @pytest.mark.parametrize("var", _LAUNCHER_VARS)
    def test_launcher_env_demotes_nonzero_rank(self, serial_env, monkeypatch, var):
        """The EME-only path imports neither meep nor mpi4py, so without this
        every one of 48 ranks would call itself master — the original bug."""
        monkeypatch.setenv(var, "7")
        assert not _is_master()

    def test_launcher_env_rank_zero_is_master(self, serial_env, monkeypatch):
        monkeypatch.setenv("SLURM_PROCID", "0")
        assert _is_master()

    def test_unparseable_rank_falls_through(self, serial_env, monkeypatch):
        monkeypatch.setenv("SLURM_PROCID", "not-a-number")
        assert _is_master()  # must not raise

    def test_meep_wins_over_env(self, serial_env, monkeypatch):
        """meep owns MPI wherever it is loaded, so it outranks the environment.
        am_really_master(), not am_master(): the latter is per-group under
        divide_parallel_processes and would elect several writers."""
        class FakeMeep:
            @staticmethod
            def am_really_master():
                return False

        monkeypatch.setenv("SLURM_PROCID", "0")
        monkeypatch.setitem(store.sys.modules, "meep", FakeMeep)
        assert not _is_master()

    def test_am_master_is_not_used(self, serial_env, monkeypatch):
        """Guards the group-vs-global distinction: under
        divide_parallel_processes am_master() is true once per group, so using
        it would elect several concurrent writers."""
        class GroupMasterButNotGlobal:
            @staticmethod
            def am_really_master():
                return False

            @staticmethod
            def am_master():
                return True

        monkeypatch.setitem(store.sys.modules, "meep", GroupMasterButNotGlobal)
        assert not _is_master()


class TestOnlyMasterWrites:
    def test_save_run_writes_nothing_off_master(self, tmp_path, monkeypatch):
        monkeypatch.setattr(store, "_is_master", lambda: False)
        save_dir = str(tmp_path)
        config = _config()

        folder = save_run(save_dir, config, np.array([1.0]))

        assert not os.path.exists(folder)
        assert os.listdir(save_dir) == []  # no index.jsonl either

    def test_save_run_still_returns_the_folder_path_off_master(self, tmp_path, monkeypatch):
        """Callers feed this straight into update_family_pointer(), so the
        return value has to be right on every rank, not just the writer."""
        config = _config()
        monkeypatch.setattr(store, "_is_master", lambda: True)
        master_folder = save_run(str(tmp_path), config, np.array([1.0]))
        monkeypatch.setattr(store, "_is_master", lambda: False)
        worker_folder = save_run(str(tmp_path), config, np.array([1.0]))

        assert worker_folder == master_folder

    def test_master_write_is_readable_by_lookup(self, tmp_path, monkeypatch):
        monkeypatch.setattr(store, "_is_master", lambda: True)
        config = _config()
        save_run(str(tmp_path), config, np.array([1.0]))
        assert lookup(str(tmp_path), config) is not None

    def test_family_pointer_not_written_off_master(self, tmp_path, monkeypatch):
        monkeypatch.setattr(store, "_is_master", lambda: False)
        update_family_pointer(str(tmp_path), {"a": 1}, (6, 1), "somehash")
        assert not os.path.isfile(_family_pointer_path(str(tmp_path), {"a": 1}))

    def test_index_jsonl_gets_one_line_per_job_not_per_rank(self, tmp_path, monkeypatch):
        """48 ranks appending to one shared file is what interleaved it."""
        save_dir = str(tmp_path)
        for rank in range(48):
            monkeypatch.setattr(store, "_is_master", lambda r=rank: r == 0)
            save_run(save_dir, _config(), np.array([1.0]))

        with open(os.path.join(save_dir, "index.jsonl")) as f:
            lines = [line for line in f if line.strip()]
        assert len(lines) == 1
        json.loads(lines[0])  # intact, not interleaved


class TestTempFileNaming:
    def test_temp_name_is_host_qualified(self, tmp_path, monkeypatch):
        """PIDs are allocated per host, so a PID-only temp name is not unique
        on a filesystem shared between nodes.

        The naming now lives in `beef._utils.atomic_replace_write`, which
        `_atomic_write_json` delegates to — the example sweeps needed the same
        guarantee for their per-task outputs (dc_validation.md item 6), and one
        hardened implementation beats eight copies. This test still covers the
        cache's own call site; `tests/test_atomic_write.py` covers the helper
        directly."""
        tmps = []
        real_replace = os.replace

        def spy_replace(src, dst):
            tmps.append(str(src))
            return real_replace(src, dst)

        monkeypatch.setattr(_utils.os, "replace", spy_replace)
        target = str(tmp_path / "pointer.json")
        for host in ("node-a", "node-b"):
            monkeypatch.setattr(
                _utils, "socket",
                type("S", (), {"gethostname": staticmethod(lambda h=host: h)}),
            )
            _atomic_write_json(target, {"dims": [1], "result_hash": "h"})

        assert len(tmps) == 2
        assert tmps[0] != tmps[1], "same temp path on two hosts — the race is still open"
        assert "node-a" in tmps[0] and "node-b" in tmps[1]
        assert str(os.getpid()) in tmps[0]

    def test_write_is_still_atomic_and_leaves_no_temp(self, tmp_path):
        target = str(tmp_path / "pointer.json")
        _atomic_write_json(target, {"dims": [4, 1], "result_hash": "h"})
        with open(target) as f:
            assert json.load(f) == {"dims": [4, 1], "result_hash": "h"}
        assert os.listdir(tmp_path) == ["pointer.json"]


class TestCorruptPointerIsAMiss:
    """A pointer is a reuse hint, never data — the result it names is also
    reachable through the plain lookup() path. So damage left by the pre-fix
    race must degrade to a cache miss, never an exception."""

    def _corrupt(self, tmp_path, contents):
        path = _family_pointer_path(str(tmp_path), {"a": 1})
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write(contents)
        return path

    @pytest.mark.parametrize("contents", ["", '{"dims": [6, 1], "result_h', "   ", "null", "[]"])
    def test_unreadable_pointer_reads_as_absent(self, tmp_path, contents):
        self._corrupt(tmp_path, contents)
        assert _read_family_pointer(_family_pointer_path(str(tmp_path), {"a": 1})) is None

    def test_missing_keys_read_as_absent(self, tmp_path):
        self._corrupt(tmp_path, '{"dims": [6, 1]}')  # no result_hash
        assert _read_family_pointer(_family_pointer_path(str(tmp_path), {"a": 1})) is None

    def test_lookup_superset_misses_instead_of_raising(self, tmp_path):
        self._corrupt(tmp_path, "")
        assert lookup_superset(str(tmp_path), {"a": 1}, (2, 1)) is None

    def test_empty_pointer_is_healed_by_the_next_save(self, tmp_path, monkeypatch):
        """This is the crash site: update_family_pointer() raised here on an
        empty file and took the whole 48-rank job down with it."""
        monkeypatch.setattr(store, "_is_master", lambda: True)
        path = self._corrupt(tmp_path, "")

        folder = save_run(str(tmp_path), _config(), np.zeros((3, 6, 1)))
        update_family_pointer(str(tmp_path), {"a": 1}, (6, 1), os.path.basename(folder))

        with open(path) as f:
            assert json.load(f)["dims"] == [6, 1]
        assert lookup_superset(str(tmp_path), {"a": 1}, (2, 1)) is not None

    def test_valid_pointer_still_blocks_a_non_dominating_update(self, tmp_path, monkeypatch):
        """The healing path must not have turned every update into a clobber."""
        monkeypatch.setattr(store, "_is_master", lambda: True)
        save_dir = str(tmp_path)
        folder = save_run(save_dir, _config(), np.zeros((3, 12, 1)))
        update_family_pointer(save_dir, {"a": 1}, (12, 1), os.path.basename(folder))
        update_family_pointer(save_dir, {"a": 1}, (2, 1), "smaller_hash")

        with open(_family_pointer_path(save_dir, {"a": 1})) as f:
            assert json.load(f)["dims"] == [12, 1]
