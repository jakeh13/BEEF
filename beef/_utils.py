# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Georgia Tech Research Corporation
import os
import socket
from typing import Callable, Dict, List, Tuple

import numpy as np


def atomic_replace_write(path: str, write_fn: Callable[[str], None]) -> None:
    """Call ``write_fn(tmp_path)``, then atomically rename onto ``path``.

    Use this for **any** file a concurrent reader might pick up — most of all
    the per-task outputs a SLURM array's ``aggregate()`` step reads from its
    peers. A plain write is not atomic: a reader that opens the path between
    the writer's ``open()`` and its last flush sees a truncated file. For
    ``.npz`` that surfaces as ``zipfile.BadZipFile``, which is what killed 7 of
    9 tasks in job 11517245 (dc_validation.md item 6). ``os.replace`` is atomic
    on POSIX, so a reader sees either the old file or the new one, never a
    partial one.

    The temp name carries **hostname as well as PID** because PIDs are
    allocated per host, so a PID-only name is not unique on a shared
    filesystem spanning nodes — see dc_validation.md item 4, where that was one
    candidate explanation for stub-sweep tasks 4 and 5 finding a path empty.

    **``np.savez`` needs the open-file-object form here.** Given a *path*
    without a ``.npz`` suffix it appends one, so it would write to
    ``<tmp>.npz`` and leave ``<tmp>`` missing, making the ``os.replace`` fail
    after the work was already done. Pass a file object instead — and close it
    with ``with``, because ``zipfile`` does **not** close a file object it was
    handed, so a bare ``open()`` would leave the flush to refcounting and could
    let ``os.replace`` move a truncated file, which is the very failure this
    helper exists to prevent::

        def _write(tmp):
            with open(tmp, "wb") as fh:
                np.savez(fh, **payload)

        atomic_replace_write(path, _write)

    Args:
        path: final destination path.
        write_fn: callable receiving the temp path to write to. It must create
            exactly that path — see the ``np.savez`` note above.
    """
    tmp = f"{path}.tmp.{socket.gethostname()}.{os.getpid()}"
    write_fn(tmp)
    if not os.path.exists(tmp):
        raise RuntimeError(
            f"write_fn did not create the temp path it was given ({tmp!r}). "
            f"The usual cause is np.savez() called with a path rather than an "
            f"open file object: it appends '.npz' to a suffix-less name, so "
            f"the data landed at '{tmp}.npz' and {path!r} was never written."
        )
    os.replace(tmp, path)


def extrude_polygons(gds_polys: Dict[Tuple[int, int], List], extrusion_rules: dict) -> list:
    """Extrude a {layer: [polygons]} dict given meow GdsExtrusionRules.

    Reimplements meow.extrude_gds()'s own loop (see
    meow_MEEP/src/meow/gds_structures.py) directly over a plain polygon dict
    instead of a gdsfactory/gdstk Cell — BEEF geometries build their shapes
    as raw vertex lists already (see DeviceGeometry.polygons()), so there's
    no need to round-trip them through a GDS Component just to hand them to
    meow; each GdsExtrusionRule is directly callable on a bare polygon.

    Args:
        gds_polys: layer tuple -> list of polygons, each a list of (x, y)
            vertex tuples (µm) — the format DeviceGeometry.to_gds_eme()/
            section_gds() return.
        extrusion_rules: layer tuple -> list[meow.GdsExtrusionRule], e.g.
            MaterialStack.extrusion_rules.

    Returns:
        List of meow Structure3D objects.
    """
    structs = []
    for layer, polys in gds_polys.items():
        if layer not in extrusion_rules:
            continue
        for poly in polys:
            arr = np.asarray(poly, dtype=float)
            structs.extend(rule(arr) for rule in extrusion_rules[layer])
    return structs


def resolve_output_folder(base: str) -> str:
    """Return base if it does not exist, otherwise base_00, base_01, …"""
    if not os.path.exists(base):
        return base
    for i in range(100):
        candidate = f"{base}_{i:02d}"
        if not os.path.exists(candidate):
            return candidate
    raise RuntimeError(f"Could not find a free output folder after 100 attempts for base '{base}'")


def make_output_folders(output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/layout", exist_ok=True)
    os.makedirs(f"{output_dir}/layout/cross_sections", exist_ok=True)
    os.makedirs(f"{output_dir}/layout/monitors", exist_ok=True)
    os.makedirs(f"{output_dir}/data", exist_ok=True)
    os.makedirs(f"{output_dir}/logs", exist_ok=True)
    os.makedirs(f"{output_dir}/fields", exist_ok=True)
