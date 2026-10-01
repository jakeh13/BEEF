# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Georgia Tech Research Corporation
"""Rank-aware logging for the `beef` logger.

Every rank of an `srun -n 48` job runs the whole solve redundantly, so every
INFO beef emits comes out 48 times. In a real batch log that was 48 of 107
lines — one message, repeated once per rank, crowding out the results the log
exists to show.

The messages are duplicated because they describe *shared* state: "cache saved
at ...", "cache hit at ...". There is one cache and one write (writes are
rank-0 only — see beef.cache.store), so one line is the honest count.

Policy, applied to the `beef` logger only — never to the root logger or anyone
else's, since a library configuring global logging is its own kind of rude:

    master rank      everything, unchanged
    other ranks      WARNING and above only, tagged with the rank

Warnings and errors survive on every rank deliberately. They are the records
that *can* be rank-specific (a mode solver failing to converge on one rank's
slice, say), so suppressing them to de-duplicate INFO would trade noise for
blindness. They get a rank tag so a genuinely one-rank problem is visible as
one rank's problem.
"""
from __future__ import annotations

import logging
import os
import sys

_LAUNCHER_RANK_VARS = ("SLURM_PROCID", "OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK")


def mpi_rank() -> int:
    """This process's global MPI rank, or 0 when not running under MPI.

    Unlike beef.cache.store._is_master() — which only needs a bool and takes
    meep's am_really_master() as authoritative — this wants the actual number,
    to tag records with. meep cannot cheaply provide a *global* rank (my_rank()
    is per-group under divide_parallel_processes), so the launcher's
    environment is the primary source here and meep only breaks the tie between
    "rank 0" and "some other rank" when no launcher variable is set.

    Imports nothing that isn't already loaded: importing mpi4py into a serial
    process would call MPI_Init behind meep's back.
    """
    for var in _LAUNCHER_RANK_VARS:
        value = os.environ.get(var)
        if value is not None:
            try:
                return int(value)
            except ValueError:
                pass
    mpi4py_mpi = sys.modules.get("mpi4py.MPI")
    if mpi4py_mpi is not None:
        try:
            return int(mpi4py_mpi.COMM_WORLD.Get_rank())
        except Exception:
            pass
    mp = sys.modules.get("meep")
    if mp is not None:
        try:
            return 0 if mp.am_really_master() else 1
        except Exception:
            pass
    return 0


class RankFilter(logging.Filter):
    """Drop sub-WARNING records off the master rank; tag what survives.

    The rank is read per-record rather than cached at construction because
    beef is imported before meep initializes MPI — at import time every process
    still looks like rank 0.
    """

    def filter(self, record: logging.LogRecord) -> bool:
        rank = mpi_rank()
        if rank == 0:
            return True
        if record.levelno < logging.WARNING:
            return False
        if not getattr(record, "_beef_rank_tagged", False):
            record.msg = f"[rank {rank}] {record.msg}"
            record._beef_rank_tagged = True
        return True


def install(name: str = "beef") -> logging.Logger:
    """Attach the filter to one logger, idempotently, and return it.

    **Must be called per module, not once on the `beef` parent.** A logger's
    filters are consulted only for records logged directly on *that* logger;
    records propagating up from a child are seen by the ancestor's *handlers*,
    never its filters. So a single filter on "beef" would let every
    "beef.solvers.fdtd" record straight through — which is exactly the 48x
    duplication this exists to stop. Attaching a handler instead is not an
    option either: handlers belong to the application, not to a library.

    Idempotent because stacking filters would re-tag the message once each.
    """
    logger = logging.getLogger(name)
    if not any(isinstance(f, RankFilter) for f in logger.filters):
        logger.addFilter(RankFilter())
    return logger


def get_logger(name: str) -> logging.Logger:
    """`logging.getLogger(name)` with the rank filter attached. Every module
    under beef/ should use this in place of logging.getLogger(__name__)."""
    return install(name)
