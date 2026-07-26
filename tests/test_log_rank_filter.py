"""Fast unit tests for beef._log — the rank filter on the `beef` logger.

Under `srun -n 48` every rank runs the whole solve, so every INFO beef emits
appeared 48 times: in a real batch log (job 11480098) 48 of 107 lines were one
repeated "cache saved at ..." message. Those describe one shared cache, not
per-rank state, so one line is the honest count.

Warnings must survive on every rank — they are the records that can genuinely
be rank-specific — so the tests below pin both directions, not just the
suppression.
"""
from __future__ import annotations

import logging

import pytest

from beef import _log
from beef._log import _LAUNCHER_RANK_VARS, RankFilter, install, mpi_rank


@pytest.fixture
def serial_env(monkeypatch):
    """Strip launcher variables so a test that doesn't set one looks serial —
    otherwise running the suite under srun would leak a rank in."""
    for var in _LAUNCHER_RANK_VARS:
        monkeypatch.delenv(var, raising=False)
    monkeypatch.delitem(_log.sys.modules, "mpi4py.MPI", raising=False)
    monkeypatch.delitem(_log.sys.modules, "meep", raising=False)


def _record(level=logging.INFO, msg="cache saved at /x"):
    return logging.LogRecord("beef.solvers.fdtd", level, __file__, 1, msg, None, None)


class TestRankDetection:
    def test_serial_is_rank_zero(self, serial_env):
        assert mpi_rank() == 0

    @pytest.mark.parametrize("var", _LAUNCHER_RANK_VARS)
    def test_launcher_env_gives_the_rank(self, serial_env, monkeypatch, var):
        monkeypatch.setenv(var, "17")
        assert mpi_rank() == 17

    def test_unparseable_rank_falls_through(self, serial_env, monkeypatch):
        monkeypatch.setenv("SLURM_PROCID", "not-a-number")
        assert mpi_rank() == 0  # must not raise

    def test_meep_breaks_the_tie_when_env_is_absent(self, serial_env, monkeypatch):
        class FakeMeep:
            @staticmethod
            def am_really_master():
                return False

        monkeypatch.setitem(_log.sys.modules, "meep", FakeMeep)
        assert mpi_rank() != 0


class TestFiltering:
    def test_master_keeps_info(self, serial_env):
        assert RankFilter().filter(_record()) is True

    def test_worker_drops_info(self, serial_env, monkeypatch):
        monkeypatch.setenv("SLURM_PROCID", "7")
        assert RankFilter().filter(_record()) is False

    def test_worker_drops_debug(self, serial_env, monkeypatch):
        monkeypatch.setenv("SLURM_PROCID", "7")
        assert RankFilter().filter(_record(logging.DEBUG)) is False

    @pytest.mark.parametrize("level", [logging.WARNING, logging.ERROR, logging.CRITICAL])
    def test_worker_keeps_warnings_and_worse(self, serial_env, monkeypatch, level):
        """Suppressing these to de-duplicate INFO would trade noise for
        blindness — a mode solver failing on one rank's slice must still show."""
        monkeypatch.setenv("SLURM_PROCID", "7")
        assert RankFilter().filter(_record(level)) is True

    def test_surviving_worker_records_are_rank_tagged(self, serial_env, monkeypatch):
        monkeypatch.setenv("SLURM_PROCID", "7")
        rec = _record(logging.WARNING, "source starved at 1.60 um")
        RankFilter().filter(rec)
        assert "[rank 7]" in rec.msg

    def test_master_records_are_not_tagged(self, serial_env):
        rec = _record(logging.WARNING, "source starved")
        RankFilter().filter(rec)
        assert "rank" not in rec.msg

    def test_tagging_is_not_applied_twice(self, serial_env, monkeypatch):
        """Two filters on one record (or a re-emitted record) must not stack
        prefixes."""
        monkeypatch.setenv("SLURM_PROCID", "7")
        rec = _record(logging.WARNING, "boom")
        RankFilter().filter(rec)
        RankFilter().filter(rec)
        assert rec.msg.count("[rank 7]") == 1

    def test_rank_is_read_per_record_not_cached(self, serial_env, monkeypatch):
        """beef is imported before meep initializes MPI, so at construction
        time every process still looks like rank 0."""
        f = RankFilter()
        assert f.filter(_record()) is True
        monkeypatch.setenv("SLURM_PROCID", "3")
        assert f.filter(_record()) is False


class TestInstall:
    def test_install_is_idempotent(self):
        name = "beef.test.idempotent"
        install(name)
        install(name)
        n = sum(isinstance(f, RankFilter) for f in logging.getLogger(name).filters)
        assert n == 1, "stacked filters would re-tag the message once each"

    @pytest.mark.parametrize("module", [
        "beef.solvers.fdtd", "beef.solvers.eme", "beef.solvers.hybrid",
        "beef.cache.store", "beef.cache.config", "beef.sweep.parallel",
    ])
    def test_every_beef_module_logger_is_filtered(self, module):
        """Must be per-module. A single filter on the parent `beef` logger does
        nothing for these: a logger's filters are consulted only for records
        logged on that logger, never for records propagating up from a child.
        beef.solvers.fdtd's "cache saved at ..." was the 48x offender in job
        11480098's log, and a parent-only filter would have let all 48 through.
        """
        __import__(module)
        assert any(isinstance(f, RankFilter)
                   for f in logging.getLogger(module).filters), module

    def test_root_logger_is_untouched(self):
        """A library configuring global logging is its own kind of rude."""
        assert not any(isinstance(f, RankFilter) for f in logging.getLogger().filters)

    def test_solver_info_is_suppressed_off_master(self, serial_env, monkeypatch, caplog):
        """End to end through the real logger: the exact message that appeared
        48 times in job 11480098's log."""
        import beef.solvers.fdtd  # noqa: F401
        monkeypatch.setenv("SLURM_PROCID", "9")
        logger = logging.getLogger("beef.solvers.fdtd")
        with caplog.at_level(logging.INFO, logger="beef"):
            logger.info("solve: cache saved at /x")
        assert "cache saved" not in caplog.text

    def test_solver_info_survives_on_master(self, serial_env, caplog):
        import beef.solvers.fdtd  # noqa: F401
        logger = logging.getLogger("beef.solvers.fdtd")
        with caplog.at_level(logging.INFO, logger="beef"):
            logger.info("solve: cache saved at /x")
        assert "cache saved" in caplog.text

    def test_solver_warning_survives_off_master(self, serial_env, monkeypatch, caplog):
        """check_source_coverage()'s warning must reach the log from any rank."""
        import beef.solvers.fdtd  # noqa: F401
        monkeypatch.setenv("SLURM_PROCID", "9")
        logger = logging.getLogger("beef.solvers.fdtd")
        with caplog.at_level(logging.INFO, logger="beef"):
            logger.warning("source is starved at 1.60 um")
        assert "source is starved" in caplog.text
        assert "[rank 9]" in caplog.text
