"""On-disk lookup/save for the persistent solve_section() cache (see
beef/cache/__init__.py). Each matching config gets its own subfolder, named
{slug}_{hash(config)} — a human-readable descriptor prefix (solver type,
method, geometry class) plus a hash suffix for O(1) lookup/uniqueness (see
_slug()/_folder()) — holding config.json (human-readable, also re-checked by
full equality against a hash collision or schema drift), result.pkl, and
geometry.png. lookup() also checks the pre-slug bare-hash folder name as a
fallback, so a save_dir populated before folder names got a slug prefix
doesn't lose anything (Solver recommendation #3). save_dir also accumulates
a top-level index.jsonl (one line per save_run() call) so the whole cache
can be grep/jq-scanned without opening every subfolder's config.json — see
_manifest_entry()'s docstring. lookup_superset()/update_family_pointer() add
a second, orthogonal lookup path — see their own docstrings — that lets a
smaller mode-count request reuse a larger mode-count result already on disk.

Reads happen on every rank; **writes happen only on the master rank** (see
_is_master()/_barrier()). Under `srun -n 48` every rank runs the whole solve
redundantly and so reached every write path here, which raced on the shared
filesystem — dc_validation.md item 4.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import pickle
import re
import socket
import sys
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
from .._log import get_logger

logger = get_logger(__name__)

_FAMILY_POINTER_DIR = "_mode_family_pointers"


_LAUNCHER_RANK_VARS = ("SLURM_PROCID", "OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK")


def _is_master() -> bool:
    """True on exactly one process per job — the only one allowed to write to
    save_dir.

    Deliberately a bool and not a rank: nothing here needs to know *which*
    rank it is, and meep cannot cheaply say (see below), so returning a number
    would mean inventing one that later callers could misuse.

    Sources are consulted in order of authority, and **nothing is imported
    that isn't already loaded** — importing mpi4py into a serial process would
    call MPI_Init behind meep's back:

    1. meep, if the process has imported it. It owns MPI wherever it is
       loaded, which is every FDTD run. `am_really_master()`, not
       `am_master()`: the latter is per-group under
       `divide_parallel_processes` and would elect one writer per group.
    2. mpi4py, if something else already initialized it.
    3. The launcher's environment. EME-only runs import neither of the above
       (meow/MPB don't pull meep in) yet are still launched under `srun -n 48`,
       so without this every rank would call itself master — exactly the bug
       this exists to prevent. srun sets SLURM_PROCID; mpirun sets the
       OMPI/PMI variables.
    """
    mp = sys.modules.get("meep")
    if mp is not None:
        try:
            return bool(mp.am_really_master())
        except Exception:  # meep present but MPI not initialized yet
            pass
    mpi4py_mpi = sys.modules.get("mpi4py.MPI")
    if mpi4py_mpi is not None:
        try:
            return int(mpi4py_mpi.COMM_WORLD.Get_rank()) == 0
        except Exception:
            pass
    for var in _LAUNCHER_RANK_VARS:
        value = os.environ.get(var)
        if value is not None:
            try:
                return int(value) == 0
            except ValueError:
                pass
    return True


def _barrier() -> None:
    """Hold all ranks until the master's cache write is visible, so a later
    lookup() cannot return a hit on one rank and a miss on another. That
    divergence is the dangerous failure: ranks disagreeing about whether to
    simulate leaves meep's collective calls unmatched and the job hangs until
    it burns its wall clock, which is exactly what happened to stub-sweep
    tasks 4 and 5.

    Only meep can do this — a run that never imported it has no collective
    calls to desynchronize, so skipping the barrier there is safe (the worst
    case is one rank redundantly recomputing an EME result it could have
    loaded). Callers must therefore invoke this from *all* ranks or none: every
    save_run()/update_family_pointer() call site is guarded by conditions
    derived from the solver config, never from per-rank data, so they do.
    """
    mp = sys.modules.get("meep")
    if mp is None:
        return
    try:
        mp.all_wait()
    except Exception:
        logger.debug("cache: MPI barrier unavailable", exc_info=True)


class CacheMissError(RuntimeError):
    """Raised when cache_only=True but no cached run matches the requested
    config — the caller asked to never fall back to actually simulating."""


def _hash_config(config: Dict[str, Any]) -> str:
    blob = json.dumps(config, sort_keys=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:16]


def _slug(config: Dict[str, Any]) -> str:
    """Short filesystem-safe descriptor derived from config — e.g.
    "FDTDSolver_solve_section_DirectionalCoupler" — prefixed onto a cache
    folder's name so `ls save_dir` is browsable without opening
    index.jsonl/config.json first (Solver recommendation #3). Purely
    cosmetic: the hash suffix _folder() appends is still what's actually
    trusted for uniqueness/collision-safety, so a slug collision between two
    different configs is harmless (they still land in different folders)."""
    raw = f"{config.get('solver_type', 'Solver')}_{config.get('method', 'method')}_{config.get('geometry_class', 'Geometry')}"
    return re.sub(r"[^A-Za-z0-9_]+", "_", raw)


def _folder(save_dir: str, config: Dict[str, Any]) -> str:
    """Current (slug-prefixed) folder-naming scheme for new writes."""
    return os.path.join(save_dir, f"{_slug(config)}_{_hash_config(config)}")


def _legacy_folder(save_dir: str, config: Dict[str, Any]) -> str:
    """Pre-slug folder-naming scheme (bare hash only) — kept so cache
    entries written before folder names got a human-readable prefix stay
    reachable; save_run() only ever writes new entries under _folder()'s
    naming, never this one (see lookup()'s docstring)."""
    return os.path.join(save_dir, _hash_config(config))


def lookup(save_dir: str, config: Dict[str, Any]) -> Optional[str]:
    """Return the matching run's folder path, or None if no cached run
    matches config (hash miss, or a hash hit whose config.json doesn't
    actually equal config — a collision/schema-drift guard, not just a
    hash-trust lookup). Checks the current slug-prefixed folder name first,
    then falls back to the pre-slug bare-hash name — so renaming the
    on-disk naming scheme doesn't invalidate anything a save_dir already
    had cached; only new writes get the new name (Solver recommendation
    #3)."""
    for folder in (_folder(save_dir, config), _legacy_folder(save_dir, config)):
        config_path = os.path.join(folder, "config.json")
        if not os.path.isfile(config_path):
            continue
        with open(config_path) as f:
            cached_config = json.load(f)
        if cached_config != config:
            logger.warning(
                "cache: hash match at %s but config.json differs — treating as "
                "a miss (hash collision or schema drift).", folder,
            )
            continue
        return folder
    return None


def load_result(folder: str) -> np.ndarray:
    with open(os.path.join(folder, "result.pkl"), "rb") as f:
        return pickle.load(f)


def _save_geometry_plot(config: Dict[str, Any], path: str) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon

    fig, ax = plt.subplots(figsize=(8, 3))
    for layer_key, polys in config["geometry_polygons"].items():
        is_wg = layer_key == "2,0"
        for poly in polys:
            ax.add_patch(MplPolygon(
                np.asarray(poly), closed=True,
                facecolor="#4a4a4a" if is_wg else "#d8d8d8",
                edgecolor="none", alpha=1.0 if is_wg else 0.5,
                zorder=1 if is_wg else 0,
            ))
    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (µm, section-local)")
    ax.set_ylabel("y (µm)")
    ax.set_title(f"{config['solver_type']}.{config['method']} — {config['geometry_class']} (cached)")
    ax.grid(False)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _summarize_call_kwargs(call_kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Compact call_kwargs for one manifest line: long numeric arrays
    (wavelengths, lengths) are summarized as {count, min, max} instead of
    dumped in full, so each manifest line stays short and grep-friendly;
    everything else (x0/x1/n_modes_in/n_modes_out/stub_left/stub_right/
    geometry_kwargs/...) is kept as-is since it's already compact."""
    out: Dict[str, Any] = {}
    for key, value in call_kwargs.items():
        if isinstance(value, list) and len(value) > 4 and all(isinstance(v, (int, float)) for v in value):
            out[key] = {"count": len(value), "min": min(value), "max": max(value)}
        else:
            out[key] = value
    return out


def _manifest_entry(config: Dict[str, Any], folder: str) -> Dict[str, Any]:
    """One compact, human-readable summary of a cached config — everything
    needed to answer "what is this hash?" from a single manifest line,
    without repeating the full config.json (geometry_polygons/solver_config
    are deliberately omitted here; open that folder's own config.json for
    the complete fingerprint)."""
    return {
        "hash": _hash_config(config),
        "folder": os.path.basename(folder),
        "solver_type": config.get("solver_type"),
        "method": config.get("method"),
        "geometry_class": config.get("geometry_class"),
        "call_kwargs": _summarize_call_kwargs(config.get("call_kwargs", {})),
    }


def _append_manifest(save_dir: str, config: Dict[str, Any], folder: str) -> None:
    """Append one JSONL line to save_dir/index.jsonl summarizing this cache
    entry — purely additive bookkeeping alongside the per-folder config.json,
    so `ls`/individual folders are unaffected and existing cache entries
    remain reachable exactly as before. Not deduplicated: re-saving an
    already-cached config (e.g. re-running the same script) appends another
    line for the same hash — harmless for a browsing aid, and avoiding it
    would mean reading the whole file (and racing concurrent writers) on
    every save just to check."""
    entry = _manifest_entry(config, folder)
    with open(os.path.join(save_dir, "index.jsonl"), "a") as f:
        f.write(json.dumps(entry, sort_keys=True) + "\n")


def save_run(save_dir: str, config: Dict[str, Any], result: np.ndarray) -> str:
    """Write config.json + result.pkl + geometry.png to a new/overwritten
    subfolder of save_dir, named {slug}_{hash(config)} (Solver recommendation
    #3 — see _folder()); also append a summary line to save_dir/index.jsonl.
    Returns the folder path.

    **Only the master rank writes.** Every rank runs the whole solve
    redundantly under `srun`, so without this all 48 would race on the same
    config.json/result.pkl and interleave their appends to the single shared
    index.jsonl. Non-master ranks still get the folder path back, since callers
    pass it to update_family_pointer(). Collective: called from all ranks or
    none — see _barrier()."""
    folder = _folder(save_dir, config)
    if _is_master():
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, "config.json"), "w") as f:
            json.dump(config, f, indent=2, sort_keys=True)
        with open(os.path.join(folder, "result.pkl"), "wb") as f:
            pickle.dump(result, f)
        try:
            _save_geometry_plot(config, os.path.join(folder, "geometry.png"))
        except Exception:
            logger.warning("cache: failed to save geometry plot at %s", folder, exc_info=True)
        try:
            _append_manifest(save_dir, config, folder)
        except Exception:
            logger.warning("cache: failed to append manifest entry at %s", save_dir, exc_info=True)
    _barrier()
    return folder


def _atomic_write_json(path: str, data: Dict[str, Any]) -> None:
    """write-to-tmp-then-os.replace so a reader never sees a half-written
    file — same pattern the example sweep scripts already use for their own
    aggregate outputs, needed here because update_family_pointer() can be
    called by concurrent SLURM array tasks writing the same family's pointer
    file (see that function's docstring).

    The temp name carries host as well as PID because PIDs are allocated per
    host, so a PID-only name is not unique on a shared filesystem spanning
    nodes. That collision is one way stub-sweep tasks 4 and 5 could have found
    `path` empty; concurrent os.replace/read against the shared filesystem from
    48 ranks is an equally good explanation, and a 48-rank 2-node probe on
    2026-07-25 saw no PID reuse, so the specific mechanism is unproven. The
    rank gate on the callers removes every version of it; this is belt-and-
    braces for any future concurrent writer."""
    tmp = f"{path}.tmp.{socket.gethostname()}.{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(data, f)
    os.replace(tmp, path)


def _family_pointer_path(save_dir: str, family_config: Dict[str, Any]) -> str:
    name = f"{_slug(family_config)}_{_hash_config(family_config)}.json"
    return os.path.join(save_dir, _FAMILY_POINTER_DIR, name)


def _read_family_pointer(path: str) -> Optional[Dict[str, Any]]:
    """The pointer file's contents, or None if it is absent, empty, truncated,
    or missing the keys it should have.

    A pointer is a reuse hint, never data — the result it names is also
    reachable through the plain lookup() path — so an unreadable one is a
    cache miss, exactly like lookup()'s config-mismatch guard, and never an
    exception. Empty and truncated files are the specific damage the pre-fix
    multi-rank write race left behind, so any left on disk now heal on the
    next write instead of killing the job that reads them."""
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as f:
            pointer = json.load(f)
    except (json.JSONDecodeError, OSError):
        logger.warning(
            "cache: family pointer at %s is unreadable — treating as absent. "
            "This is the fingerprint of the pre-fix concurrent-write race; the "
            "file will be rewritten by the next save.", path,
        )
        return None
    if not isinstance(pointer, dict) or "dims" not in pointer or "result_hash" not in pointer:
        logger.warning(
            "cache: family pointer at %s is missing dims/result_hash — "
            "treating as absent (schema drift).", path,
        )
        return None
    return pointer


def lookup_superset(
    save_dir: str, family_config: Dict[str, Any], dims: Sequence[int],
) -> Optional[Tuple[str, Tuple[int, ...]]]:
    """Look up whether a cached run already computed at least `dims` modes
    for this family — a growing group of calls that are physically
    interchangeable except for how many modes they carried (same solver
    config, geometry, and — for FDTDSolver — same excited/measured side and
    a pinned monitor_size, so mode count is the *only* thing that varies the
    underlying simulation; see solve_section()'s docstring for why those two
    preconditions are the caller's responsibility, not this function's).

    `family_config` must already exclude the mode-count fields themselves
    (n_modes_in/n_modes_out or equivalent) — this is exactly what
    build_run_config()/build_whole_device_config() return when handed
    call_kwargs with those fields omitted, so no separate config builder is
    needed for this.

    `dims` is the *requested* mode counts, in the same order/axes as the
    result array's trailing dimensions (e.g. (n_modes_out, n_modes_in) to
    match a [wavelengths, n_modes_out, n_modes_in]-shaped result — the
    caller, not this function, knows and enforces that convention).

    Returns (folder, stored_dims) if a cached run's own dims dominate the
    request component-wise (stored >= requested in every position) — the
    caller loads that folder's result and slices its trailing axes down to
    `dims`. Returns None on any miss: no pointer yet, dims don't dominate, or
    the pointer refers to a folder that's gone missing.
    """
    path = _family_pointer_path(save_dir, family_config)
    pointer = _read_family_pointer(path)
    if pointer is None:
        return None
    stored_dims = tuple(pointer["dims"])
    if len(stored_dims) != len(dims):
        logger.warning(
            "cache: family pointer at %s has %d dims but %d were requested — "
            "treating as a miss (schema drift).", path, len(stored_dims), len(dims),
        )
        return None
    if not all(s >= r for s, r in zip(stored_dims, dims)):
        return None
    folder = os.path.join(save_dir, pointer["result_hash"])
    if not os.path.isfile(os.path.join(folder, "result.pkl")):
        logger.warning(
            "cache: family pointer at %s refers to missing folder %s — "
            "treating as a miss.", path, folder,
        )
        return None
    return folder, stored_dims


def update_family_pointer(
    save_dir: str, family_config: Dict[str, Any], dims: Sequence[int], result_hash: str,
) -> None:
    """Record that `result_hash` (a folder already written by save_run())
    was computed at `dims` modes for this family — but only if `dims`
    dominates whatever's already recorded (component-wise >=), or nothing is
    recorded yet. A result that doesn't dominate the existing pointer (e.g.
    the family's dims grew on one axis and shrank on another) is simply not
    recorded as the new family maximum; it's still saved under its own exact
    hash by save_run() and reachable via the plain lookup()/exact-config
    path, just not eligible to serve future superset lookups for this
    family. This is a read-then-atomic-write, not a fully atomic
    check-and-set.

    **Only the master rank writes** (see save_run()); collective, so call it
    from all ranks or none. Across *separate processes* — the case this cannot
    serialize, e.g. a SLURM array sweeping mode count, where each task is its
    own job — the worst case is a missed reuse opportunity (a smaller dims
    value briefly recorded after a larger one, later overwritten by the next
    call that dominates it), never an incorrect one: lookup_superset() only
    trusts a pointer whose target folder actually exists and whose recorded
    dims came from that exact result_hash, so a stale or lost update can only
    cause an unnecessary resimulation, not a wrong answer.

    That guarantee did *not* previously extend to ranks within a single job,
    which is what killed stub-sweep tasks 4 and 5: they raised here reading a
    pointer file that had ended up empty under 48 concurrent writers (see
    _atomic_write_json() for why the exact mechanism is unproven). An
    unreadable pointer is now treated as absent and overwritten rather than
    raised on, so a cache can never be the thing that kills a run.
    """
    if _is_master():
        path = _family_pointer_path(save_dir, family_config)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        existing = _read_family_pointer(path)
        if existing is None or all(
            new >= old for new, old in zip(dims, tuple(existing["dims"]))
        ):
            _atomic_write_json(path, {"dims": list(dims), "result_hash": result_hash})
    _barrier()
