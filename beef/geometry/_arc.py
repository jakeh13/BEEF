# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Georgia Tech Research Corporation
"""Quarter-annulus wedge polygon — shared primitive for any constant-radius
bend built from a fixed center/radius/angle-span rather than lateral
endpoint offsets (contrast beef.geometry.sbend, whose "circular"/"euler"
modes are parametrized by (dx, dy) transitions, not a fixed sweep angle).

Same construction (and signature) as the arc_polygon() duplicated in both
examples/ring_resonator/python/ring_resonator_hybrid*.py scripts and
beef/geometry/bend.py's private _arc_polygon — promoted here so a future
caller of any of the three has one implementation to match against.
"""
from __future__ import annotations

from typing import List, Tuple

import numpy as np


def arc_polygon(
    cx: float, cy: float, r_inner: float, r_outer: float,
    theta0_deg: float, theta1_deg: float, n_pts: int = 64,
) -> List[Tuple[float, float]]:
    """Quarter- (or any-angle-) annulus wedge vertex list, device frame.

    Samples `n_pts` points along the outer radius from theta0_deg to
    theta1_deg, then the inner radius back from theta1_deg to theta0_deg,
    giving one closed polygon (outer arc forward + inner arc reversed).

    Args:
        cx, cy: arc center (µm).
        r_inner, r_outer: inner/outer radius of the wedge (µm).
        theta0_deg, theta1_deg: sweep start/end angle (degrees, standard
            math convention — 0 = +x axis, increasing counterclockwise).
        n_pts: samples per arc (both inner and outer use the same count).
    """
    t0, t1 = np.radians(theta0_deg), np.radians(theta1_deg)
    thetas = np.linspace(t0, t1, n_pts)
    outer = [(cx + r_outer * np.cos(t), cy + r_outer * np.sin(t)) for t in thetas]
    inner = [(cx + r_inner * np.cos(t), cy + r_inner * np.sin(t)) for t in reversed(thetas)]
    return outer + inner
