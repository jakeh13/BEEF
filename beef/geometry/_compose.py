# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Georgia Tech Research Corporation
"""Offset-translation helpers for composing DeviceGeometry sub-instances
end-to-end along x (MMI's per-port Tapers, DirectionalCoupler's SBend arms,
CompoundGeometry's segment stitching all use these)."""
from __future__ import annotations

from typing import List, Tuple


def shift_polys(polys: List[List[Tuple[float, float]]], x_off: float, y_off: float = 0.0
                 ) -> List[List[Tuple[float, float]]]:
    """Translate a list of (x, y) polygons by (x_off, y_off)."""
    return [[(x + x_off, y + y_off) for x, y in poly] for poly in polys]


def shift_prisms(prisms, x_off: float, y_off: float = 0.0):
    """Return new mp.Prism objects with every vertex shifted by (x_off, y_off)."""
    import meep as mp

    return [
        mp.Prism(
            vertices=[mp.Vector3(v.x + x_off, v.y + y_off, v.z) for v in p.vertices],
            height=p.height,
            material=p.material,
        )
        for p in prisms
    ]
