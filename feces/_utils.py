import os
from typing import Dict, List, Tuple

import numpy as np


def extrude_polygons(gds_polys: Dict[Tuple[int, int], List], extrusion_rules: dict) -> list:
    """Extrude a {layer: [polygons]} dict given meow GdsExtrusionRules.

    Reimplements meow.extrude_gds()'s own loop (see
    meow_MEEP/src/meow/gds_structures.py) directly over a plain polygon dict
    instead of a gdsfactory/gdstk Cell — FECES geometries build their shapes
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
