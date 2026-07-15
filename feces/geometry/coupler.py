"""Directional coupler geometry: two waveguides brought into proximity over a
straight coupling length, connected to widely-separated access waveguides by
lateral S-bends on each side.

    y=+offset/2  ____                              ____
                     \\                            /
     "isolated"       \\                          /   "cross"
                        \\________________________/
                         (near pitch = width + gap)
                        /                          \\
     "in"              /                            \\
                  ____/                              \\____
    y=-offset/2                                            "through"

                |<-- bend_length -->|<- coupling_length ->|<-- bend_length -->|

Only the bottom-left ("in") port is ever excited; "isolated" is present in
the geometry (it is physically coupled to "in" through the same bends and
straight section) but carries no source. This mirrors MMI's n_inputs=1
convention for a symmetric N-port device.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from .base import DeviceGeometry

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


def _oxide_box(x_start: float, total_x: float, oxide_width: float) -> List[Tuple[float, float]]:
    return [
        (x_start, -oxide_width / 2),
        (total_x,  -oxide_width / 2),
        (total_x,   oxide_width / 2),
        (x_start,  oxide_width / 2),
    ]


def _arm_sign(idx: int) -> float:
    """0 -> bottom arm (-y), 1 -> top arm (+y). Matches MMI._port_ys ordering."""
    return -1.0 if idx == 0 else 1.0


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


class DirectionalCoupler(DeviceGeometry):
    """Symmetric 2-waveguide directional coupler with S-bend I/O transitions.

    Layout (device frame, x=0 at the coupling-region centre):
        access wg | S-bend | straight coupling region | S-bend | access wg
    Both arms run the full length; only the bottom-left arm ("in") is meant
    to be excited, giving 3 physical ports of interest: "in", "through"
    (bottom-right, same arm as "in"), and "cross" (top-right, the coupled
    arm). The top-left arm ("isolated") is present in the geometry but
    unexcited, matching a standard 4-port coupler where the 4th port carries
    negligible power by design.

    Args:
        waveguide_width: single-mode strip width (µm).
        gap: edge-to-edge separation between the two waveguides in the
            coupling region (µm).
        coupling_length: straight parallel interaction length (µm).
        bend_length: horizontal extent of each S-bend transition (µm).
        waveguide_offset: centre-to-centre pitch between the two arms at the
            access-waveguide ends (µm) - must be large enough that the two
            arms are effectively uncoupled there (several waveguide widths
            beyond the evanescent tail).
        bend_type: "euler" or "circular" (see _sbend_polygon).
        waveguide_height: Si layer thickness (µm).
    """

    def __init__(
        self,
        waveguide_width: float = 0.45,
        gap: float = 0.20,
        coupling_length: float = 8.0,
        bend_length: float = 12.0,
        waveguide_offset: float = 3.0,
        bend_type: str = "euler",
        waveguide_height: float = 0.22,
    ):
        self.waveguide_width = waveguide_width
        self.gap = gap
        self.coupling_length = coupling_length
        self.bend_length = bend_length
        self.waveguide_offset = waveguide_offset
        self.bend_type = bend_type
        self.waveguide_height = waveguide_height

    @property
    def coupling_pitch(self) -> float:
        """Centre-to-centre spacing between the two arms in the coupling region (µm)."""
        return self.waveguide_width + self.gap

    @property
    def total_length(self) -> float:
        return self.coupling_length + 2 * self.bend_length

    @property
    def input_ys(self) -> List[float]:
        """y-centres of the two input-side access waveguides, low to high:
        index 0 = bottom arm ("in", the one meant to be excited), index 1 =
        top arm ("isolated") — same ordering as _arm_sign()/MMI._port_ys.
        FDTDSolver.solve() reads this to place its source: port="bot" picks
        index 0 ("in"), port="top" index 1 ("isolated")."""
        return [-self.waveguide_offset / 2, +self.waveguide_offset / 2]

    # ------------------------------------------------------------------
    # Whole-device polygons (plain FDTD / geometry preview)
    # ------------------------------------------------------------------

    def polygons(
        self, stub_length: float = 0.0
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        hc = self.coupling_length / 2
        bl = self.bend_length
        ww = self.waveguide_width
        p_near = self.coupling_pitch / 2
        p_far = self.waveguide_offset / 2

        wg: List[List[Tuple[float, float]]] = []
        if hc > 1e-9:  # coupling_length=0 -> bends back-to-back, no straight section
            for s in (-1.0, 1.0):
                y = s * p_near
                wg.append([(-hc, y - ww/2), (hc, y - ww/2), (hc, y + ww/2), (-hc, y + ww/2)])
        for s in (-1.0, 1.0):
            wg.extend(_sbend_polygon(-hc - bl, s*p_far, -hc, s*p_near, ww, self.bend_type))
            wg.extend(_sbend_polygon(hc, s*p_near, hc + bl, s*p_far, ww, self.bend_type))

        if stub_length > 0:
            sl = stub_length
            ml2 = self.total_length / 2
            for s in (-1.0, 1.0):
                y = s * p_far
                wg.append([(-ml2 - sl, y - ww/2), (-ml2, y - ww/2),
                           (-ml2, y + ww/2), (-ml2 - sl, y + ww/2)])
                wg.append([(ml2, y - ww/2), (ml2 + sl, y - ww/2),
                           (ml2 + sl, y + ww/2), (ml2, y + ww/2)])
        return {_WG_LAYER: wg}

    def to_gds_eme(
        self, straight_length: float, oxide_width: float,
        x_start: float = 0.0, **kwargs,
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        ml = self.total_length
        x_off = x_start + ml / 2 + straight_length
        wg_polys = [
            [(x + x_off, y) for x, y in verts]
            for verts in self.polygons(stub_length=straight_length)[_WG_LAYER]
        ]
        ox_poly = _oxide_box(x_start, x_start + ml + 2 * straight_length, oxide_width)
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    # ------------------------------------------------------------------
    # Arbitrary-section interface (HybridSolver)
    # ------------------------------------------------------------------
    #
    # Device frame: x=0 at the coupling-region centre.
    #   x in [-hc-bl, -hc] : left S-bend  (near pitch -> far pitch, going left)
    #   x in [-hc,    +hc] : coupling region (uniform, both arms present)
    #   x in [+hc, +hc+bl] : right S-bend (near pitch -> far pitch, going right)
    # Anything beyond +-(hc+bl) is the (virtual, arbitrary-length) access
    # waveguide, uniform width ww, at y = +-waveguide_offset/2.
    #
    # input_idx/output_idx select which arm (0=bottom, 1=top) a cut/stub on
    # the left/right side resolves to once it reaches into the bend or
    # access-waveguide region - irrelevant while the cut stays inside the
    # (always-both-arms) coupling region. Only the coupling region needs
    # both arms simultaneously: past either bend, the two arms are far
    # enough apart (waveguide_offset) that each is modelled independently,
    # the same approximation MMI already makes for its output branches.

    def _bend_bounds(self):
        hc = self.coupling_length / 2
        bl = self.bend_length
        return (-hc - bl, -hc), (hc, hc + bl)

    def _check_cut(self, x: float) -> None:
        (in0, in1), (out0, out1) = self._bend_bounds()
        if in0 < x < in1 or out0 < x < out1:
            raise ValueError(
                f"Cut/stub at x={x:.4f} falls strictly inside an S-bend "
                "(cross-section centre varies there) - cuts must land on a "
                "bend endpoint or inside the uniform coupling region / "
                "access waveguide."
            )

    def _section_polygons(self, x0: float, x1: float,
                           stub_left: float = 0.0, stub_right: float = 0.0,
                           input_idx: int = 0, output_idx: int = 0):
        """Raw vertex lists (device frame) covering [x0-stub_left, x1+stub_right].

        Draws *both* arms everywhere, not just the input_idx/output_idx one —
        those two only select which arm width_at()/center_at() report for
        monitor sizing (i.e. which one HybridSolver is exciting/measuring);
        the small FDTD box a section occupies is a real physical cross-
        section of the whole coupler, and the other arm (e.g. "isolated" on
        the input side) really is sitting right there in it, unexcited.
        """
        (in0, in1), (out0, out1) = self._bend_bounds()
        ww = self.waveguide_width
        p_near = self.coupling_pitch / 2
        p_far = self.waveguide_offset / 2

        lo, hi = x0 - stub_left, x1 + stub_right
        self._check_cut(lo)
        self._check_cut(hi)

        polys: List[List[Tuple[float, float]]] = []

        def _rect(xa, xb, yc):
            if xb - xa > 1e-9:
                polys.append([(xa, yc-ww/2), (xb, yc-ww/2), (xb, yc+ww/2), (xa, yc+ww/2)])

        # left access waveguide (both arms, before the bend)
        if lo < in0:
            for s in (-1.0, 1.0):
                _rect(lo, min(hi, in0), s * p_far)
        # left S-bend (both arms, all-or-nothing - never cut mid-bend)
        if lo <= in0 and hi >= in1:
            for s in (-1.0, 1.0):
                polys.extend(_sbend_polygon(in0, s*p_far, in1, s*p_near, ww, self.bend_type))
        # coupling region (both arms, arbitrary sub-range)
        body_lo, body_hi = max(lo, in1), min(hi, out0)
        if body_hi - body_lo > 1e-9:
            _rect(body_lo, body_hi,  p_near)
            _rect(body_lo, body_hi, -p_near)
        # right S-bend (both arms)
        if lo <= out0 and hi >= out1:
            for s in (-1.0, 1.0):
                polys.extend(_sbend_polygon(out0, s*p_near, out1, s*p_far, ww, self.bend_type))
        # right access waveguide (both arms)
        if hi > out1:
            for s in (-1.0, 1.0):
                _rect(max(lo, out1), hi, s * p_far)

        return polys

    def section_prisms(self, x0: float, x1: float, si_medium,
                        stub_left: float = 0.0, stub_right: float = 0.0,
                        input_idx: int = 0, output_idx: int = 0,
                        wg_height: float | None = None):
        from .base import polygons_to_prisms
        return polygons_to_prisms(
            self._section_polygons(x0, x1, stub_left, stub_right, input_idx, output_idx),
            si_medium, wg_height=wg_height,
        )

    def section_gds(self, x0: float, x1: float, oxide_width: float,
                     stub_left: float = 0.0, stub_right: float = 0.0,
                     input_idx: int = 0, output_idx: int = 0
                     ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        lo, hi = x0 - stub_left, x1 + stub_right
        x_off = -lo
        wg_polys = [
            [(x + x_off, y) for x, y in verts]
            for verts in self._section_polygons(x0, x1, stub_left, stub_right, input_idx, output_idx)
        ]
        ox_poly = _oxide_box(0.0, hi - lo, oxide_width)
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    def width_at(self, x: float, input_idx: int = 0, output_idx: int = 0) -> float:
        (in0, in1), (out0, out1) = self._bend_bounds()
        if in1 <= x <= out0:
            return 2 * self.waveguide_width + self.gap  # full coupled-pair span
        return self.waveguide_width

    def center_at(self, x: float, input_idx: int = 0, output_idx: int = 0) -> float:
        (in0, in1), (out0, out1) = self._bend_bounds()
        if in1 <= x <= out0:
            return 0.0
        p_near = self.coupling_pitch / 2
        p_far = self.waveguide_offset / 2
        if x < in1:
            s = _arm_sign(input_idx)
            if x <= in0:
                return s * p_far
            frac = (x - in0) / (in1 - in0)  # linear approx; only ever evaluated eps from a boundary
            return s * (p_far + frac * (p_near - p_far))
        s = _arm_sign(output_idx)
        if x >= out1:
            return s * p_far
        frac = (x - out0) / (out1 - out0)
        return s * (p_near + frac * (p_far - p_near))

    def source_aperture_at(self, x: float, input_idx: int = 0,
                            output_idx: int = 0) -> float | None:
        """Widest eigenmode source / norm-monitor window (full y-span, µm),
        centered on center_at(x), that stays clear of the OTHER arm's core
        at this cross-section. Both arms are always physically present in
        any cross-section of the coupler, so a full-domain launch aperture
        centered on one arm would also reach the unexcited companion arm —
        corrupting the "single-mode" MPB solve into seeing the arm pair's
        near-degenerate even/odd supermodes instead of the intended isolated
        arm. Returns None inside the straight coupling region, where a
        launch/measurement legitimately must see both arms (the coupled-pair
        supermode basis); callers clamp the window to their own non-PML
        domain. Both arms mirror each other about y=0 everywhere, so the
        window is the same size for either input_idx/output_idx choice."""
        (in0, in1), (out0, out1) = self._bend_bounds()
        if in1 <= x <= out0:
            return None  # coupled-pair region — both arms belong in the basis
        idx = input_idx if x < in1 else output_idx
        c = self.center_at(x, input_idx=idx, output_idx=idx)
        # other arm's centre is the mirror image at -c: edge-to-window
        # half-span = |c - (-c)| - other arm's half width
        return 2.0 * (2.0 * abs(c) - self.waveguide_width / 2)

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.bend_length, self.coupling_length,
                self.bend_length, straight_length]

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.total_length, self.waveguide_offset
