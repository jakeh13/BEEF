# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Georgia Tech Research Corporation
"""Lateral S-bend waveguide — a standalone DeviceGeometry, usable on its own
or as the building block DirectionalCoupler composes for each of its 4 arms.

Centerline/extrusion math lives here (relocated from coupler.py, single
owner) — DirectionalCoupler no longer defines this math itself."""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from .base import DeviceGeometry, _oxide_box, polygons_to_prisms

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


def _circular_sbend_centerline(dx: float, dy: float, n_pts: int = 80):
    """Centerline of a matched-circular-arc S-bend: two constant-curvature
    arcs of equal radius R and equal sweep angle theta, curvature sign
    flipping at the midpoint, horizontal tangent at both ends.

    Closed form: dx = 2R sin(theta), dy = 2R(1 - cos(theta))
      => theta = 2*atan2(|dy|, dx),  R = dx / (2 sin(theta))
    """
    theta = 2.0 * np.arctan2(abs(dy), dx)
    R = dx / (2.0 * np.sin(theta))
    s = np.sign(dy)
    t = np.linspace(0.0, 1.0, n_pts)
    x = np.empty(n_pts)
    y = np.empty(n_pts)

    m1 = t <= 0.5
    a1 = (t[m1] / 0.5) * theta
    x[m1] = R * np.sin(a1)
    y[m1] = s * R * (1.0 - np.cos(a1))

    tp = 1.0 - t[~m1]
    a2 = (tp / 0.5) * theta
    x[~m1] = dx - R * np.sin(a2)
    y[~m1] = dy - s * R * (1.0 - np.cos(a2))

    return x, y


def _bezier_sbend_centerline(dx: float, dy: float, n_pts: int = 80):
    """Centerline of a cubic-Bezier S-bend with horizontal tangent at both
    ends — control points (0,0), (dx/2,0), (dx/2,dy), (dx,dy), the same
    layout a standard S-bend Bezier (e.g. gdsfactory's bend_s/path.bezier)
    uses. Expanding the cubic Bernstein basis at those control points gives
    the closed form below (y(t) is the familiar 3t^2-2t^3 smoothstep curve).
    """
    t = np.linspace(0.0, 1.0, n_pts)
    x = dx * (1.5 * t - 1.5 * t**2 + t**3)
    y = dy * (3.0 * t**2 - 2.0 * t**3)
    return x, y


_SBEND_CENTERLINES = {
    "euler": _bezier_sbend_centerline,
    "circular": _circular_sbend_centerline,
}


def _extrude_centerline(x: np.ndarray, y: np.ndarray, width: float) -> List[Tuple[float, float]]:
    """Constant-width ribbon polygon around a centerline, via finite-difference normals."""
    dxp = np.gradient(x)
    dyp = np.gradient(y)
    # Both S-bend centerlines (_bezier_sbend_centerline, _circular_sbend_centerline)
    # have an exactly horizontal analytic tangent at both endpoints, but
    # np.gradient's one-sided boundary estimate is only approximately
    # zero-slope (O(1/n_pts) error) — that tilts the endpoint normal
    # slightly off-vertical and leaves a small jog where the ribbon meets
    # the flush-vertical straight-section cap. Force the exact tangent at
    # the two endpoints; interior points keep the finite-difference estimate.
    dxp[0], dyp[0] = 1.0, 0.0
    dxp[-1], dyp[-1] = 1.0, 0.0
    norm = np.hypot(dxp, dyp)
    norm[norm == 0] = 1.0
    nx, ny = -dyp / norm, dxp / norm
    hw = width / 2.0
    left = list(zip(x + nx * hw, y + ny * hw))
    right = list(zip(x - nx * hw, y - ny * hw))
    return left + right[::-1]


def _sbend_polygon(
    x0: float, y0: float, x1: float, y1: float, width: float, bend_type: str,
    n_pts: int = 80,
) -> List[List[Tuple[float, float]]]:
    """Waveguide polygon(s) for a lateral S-bend from (x0,y0) to (x1,y1),
    horizontal (+x) tangent at both ends, constant width throughout.

    bend_type:
        "euler"    - cubic-Bezier centerline (see _bezier_sbend_centerline):
                     continuously-varying curvature, no curvature
                     discontinuity at the midpoint. Named to match
                     BendGeometry's "euler" (smooth curve) vs. "circular"
                     (constant-curvature arc) split; not a literal
                     Euler-spiral/clothoid.
        "circular" - two matched constant-curvature arcs (closed form above).
    """
    dx, dy = x1 - x0, y1 - y0
    if abs(dy) < 1e-9 or dx <= 0:
        hw = width / 2.0
        return [[(x0, y0 - hw), (x1, y0 - hw), (x1, y0 + hw), (x0, y0 + hw)]]

    centerline_fn = _SBEND_CENTERLINES.get(bend_type)
    if centerline_fn is None:
        raise ValueError(f"unknown bend_type {bend_type!r}, expected 'euler' or 'circular'")
    cx, cy = centerline_fn(dx, dy, n_pts=n_pts)
    poly = _extrude_centerline(cx, cy, width)
    return [[(x0 + px, y0 + py) for px, py in poly]]


class SBend(DeviceGeometry):
    """S-bend centered at x=0: (-length/2, -dy/2) -> (+length/2, +dy/2),
    horizontal tangent at both ends, constant width.

    The cross-section center varies everywhere along its length, so
    section_prisms()/section_gds() must include it whole or not at all.
    """

    def __init__(
        self,
        length: float,
        dy: float,
        width: float,
        bend_type: str = "euler",
        waveguide_height: float = 0.22,
        n_pts: int = 80,
    ):
        self.length = length
        self.dy = dy
        self.width = width
        self.bend_type = bend_type
        self.waveguide_height = waveguide_height
        self.n_pts = n_pts

    @property
    def total_length(self) -> float:
        return self.length

    def polygons(self, stub_length: float = 0.0) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        L, dy, w = self.length, self.dy, self.width
        wg = list(_sbend_polygon(-L/2, -dy/2, L/2, dy/2, w, self.bend_type, n_pts=self.n_pts))
        if stub_length > 0:
            sl = stub_length
            wg.append([(-L/2-sl, -dy/2-w/2), (-L/2, -dy/2-w/2), (-L/2, -dy/2+w/2), (-L/2-sl, -dy/2+w/2)])
            wg.append([(L/2, dy/2-w/2), (L/2+sl, dy/2-w/2), (L/2+sl, dy/2+w/2), (L/2, dy/2+w/2)])
        return {_WG_LAYER: wg}

    def to_gds_eme(
        self, straight_length: float, oxide_width: float,
        x_start: float = 0.0, **kwargs
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        L = self.length
        x_off = x_start + L / 2 + straight_length
        wg_polys = [
            [(x + x_off, y) for x, y in verts]
            for verts in self.polygons(stub_length=straight_length)[_WG_LAYER]
        ]
        ox_poly = _oxide_box(x_start, x_start + L + 2 * straight_length, oxide_width)
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    def _check_cut(self, x: float) -> None:
        if -self.length / 2 < x < self.length / 2:
            raise ValueError(
                f"Cut/stub at x={x:.4f} falls strictly inside an SBend "
                "(cross-section center varies everywhere along its length) — "
                "an SBend must be included whole or not at all in any section."
            )

    def _section_polygons(self, x0: float, x1: float,
                           stub_left: float = 0.0, stub_right: float = 0.0):
        L, dy, w = self.length, self.dy, self.width
        lo, hi = x0 - stub_left, x1 + stub_right
        self._check_cut(lo)
        self._check_cut(hi)

        polys: List[List[Tuple[float, float]]] = []

        def _rect(xa, xb, yc):
            if xb - xa > 1e-9:
                polys.append([(xa, yc-w/2), (xb, yc-w/2), (xb, yc+w/2), (xa, yc+w/2)])

        if lo < -L/2:
            _rect(lo, min(hi, -L/2), -dy/2)
        if lo <= -L/2 and hi >= L/2:
            polys.extend(_sbend_polygon(-L/2, -dy/2, L/2, dy/2, w, self.bend_type, n_pts=self.n_pts))
        if hi > L/2:
            _rect(max(lo, L/2), hi, dy/2)
        return polys

    def section_prisms(self, x0: float, x1: float, si_medium,
                        stub_left: float = 0.0, stub_right: float = 0.0,
                        wg_height: Optional[float] = None, **_kwargs):
        return polygons_to_prisms(
            self._section_polygons(x0, x1, stub_left, stub_right),
            si_medium, wg_height=wg_height,
        )

    def section_gds(self, x0: float, x1: float, oxide_width: float,
                     stub_left: float = 0.0, stub_right: float = 0.0, **_kwargs
                     ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        lo, hi = x0 - stub_left, x1 + stub_right
        x_off = -lo
        wg_polys = [
            [(x + x_off, y) for x, y in verts]
            for verts in self._section_polygons(x0, x1, stub_left, stub_right)
        ]
        ox_poly = _oxide_box(0.0, hi - lo, oxide_width)
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    def width_at(self, x: float, **_kwargs) -> float:
        return self.width

    def center_at(self, x: float, **_kwargs) -> float:
        L, dy = self.length, self.dy
        if x <= -L/2:
            return -dy/2
        if x >= L/2:
            return dy/2
        frac = (x + L/2) / L  # linear approx; only ever evaluated at/near a boundary
        return -dy/2 + frac * dy

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.length, straight_length]

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.length, self.width + abs(self.dy)
