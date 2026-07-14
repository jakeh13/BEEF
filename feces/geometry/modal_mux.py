"""Modal multiplexer geometry: an asymmetric directional coupler used for
mode-division multiplexing (MDM) — two waveguides of *different* widths
brought into proximity over a straight coupling length, connected to
widely-separated access waveguides by lateral S-bends, exactly like
DirectionalCoupler's shape.

    y=+offset/2  ____                              ____
                     \\                            /
     "access"         \\                          /   "access_out"
     (narrow,            \\________________________/   (TE0 only)
      TE0 only)          /                          \\
    y=-offset/2  ____   /                            \\____
    "bus"            \\/                                   "bus_out"
    (wide,                                                 (TE0 + TE1)
     TE0 + TE1)

                |<-- bend_length -->|<- coupling_length ->|<-- bend_length -->|

Unlike DirectionalCoupler (equal-width arms, used for power splitting), the
two arms here have deliberately *different* widths:

  - the "bus" arm is wide enough to guide both the fundamental (TE0) and
    first higher-order (TE1) mode,
  - the "access" arm is narrow enough to guide only TE0.

With access_width chosen so the access arm's TE0 effective index is close to
the bus arm's TE1 effective index (the phase-matching condition for
asymmetric directional couplers), coupling in the straight region is
mode-selective: access TE0 converts into bus TE1 while bus TE0 (phase
mismatched with the access arm) passes through largely undisturbed. That is
the multiplexing action — combining an independent signal from the access
arm onto a second mode of the shared bus waveguide.

This module reuses DirectionalCoupler's generic S-bend/oxide-box helpers
from .coupler (none of them are coupler-specific — they just extrude a
ribbon of a given width along a centerline, or build an oxide rectangle)
rather than duplicating that logic here.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from .base import DeviceGeometry
from .coupler import (
    _arm_sign,
    _oxide_box,
    _sbend_polygon,
    _OX_LAYER,
    _WG_LAYER,
)


class ModalMultiplexer(DeviceGeometry):
    """Asymmetric 2-waveguide directional coupler used as a mode-division
    multiplexer, with S-bend I/O transitions.

    Layout (device frame, x=0 at the coupling-region centre):
        access wg | S-bend | straight coupling region | S-bend | access wg
    Both arms run the full length. Arm 0 ("bus", bottom) is wide (carries
    TE0 + TE1); arm 1 ("access", top) is narrow (carries TE0 only). The
    bus arm is meant to be excited with TE0 (pass-through channel) and/or
    measured for TE0+TE1 (multiplexed output); the access arm is excited
    with TE0 (the channel being added onto the bus's TE1).

    Args:
        bus_width: wide-arm strip width (µm) — must support TE0 + TE1.
        access_width: narrow-arm strip width (µm) — single-mode (TE0 only),
            chosen so its TE0 effective index is close to the bus arm's TE1
            effective index (phase matching for mode-selective coupling).
        gap: edge-to-edge separation between the two waveguides in the
            coupling region (µm).
        coupling_length: straight parallel interaction length (µm).
        bend_length: horizontal extent of each S-bend transition (µm).
        waveguide_offset: centre-to-centre pitch between the two arms at the
            access-waveguide ends (µm) - must be large enough that the two
            arms are effectively uncoupled there.
        bend_type: "euler" or "circular" (see coupler._sbend_polygon).
        waveguide_height: Si layer thickness (µm).
    """

    def __init__(
        self,
        bus_width: float = 0.90,
        access_width: float = 0.44,
        gap: float = 0.20,
        coupling_length: float = 8.0,
        bend_length: float = 12.0,
        waveguide_offset: float = 3.0,
        bend_type: str = "euler",
        waveguide_height: float = 0.22,
    ):
        self.bus_width = bus_width
        self.access_width = access_width
        self.gap = gap
        self.coupling_length = coupling_length
        self.bend_length = bend_length
        self.waveguide_offset = waveguide_offset
        self.bend_type = bend_type
        self.waveguide_height = waveguide_height

    def _arm_width(self, idx: int) -> float:
        """0 -> bus (wide) arm, 1 -> access (narrow) arm."""
        return self.bus_width if idx == 0 else self.access_width

    def _p_near(self, idx: int) -> float:
        """Centre-line offset from y=0 for arm idx in the coupling region —
        chosen so the edge-to-edge gap between the two (different-width)
        arms is exactly self.gap."""
        return self.gap / 2 + self._arm_width(idx) / 2

    @property
    def total_length(self) -> float:
        return self.coupling_length + 2 * self.bend_length

    # ------------------------------------------------------------------
    # Whole-device polygons (plain FDTD / geometry preview)
    # ------------------------------------------------------------------

    def polygons(
        self, stub_length: float = 0.0
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        hc = self.coupling_length / 2
        bl = self.bend_length
        p_far = self.waveguide_offset / 2

        wg: List[List[Tuple[float, float]]] = []
        for idx in (0, 1):
            s = _arm_sign(idx)
            ww = self._arm_width(idx)
            y = s * self._p_near(idx)
            wg.append([(-hc, y - ww/2), (hc, y - ww/2), (hc, y + ww/2), (-hc, y + ww/2)])
            wg.extend(_sbend_polygon(-hc - bl, s*p_far, -hc, y, ww, self.bend_type))
            wg.extend(_sbend_polygon(hc, y, hc + bl, s*p_far, ww, self.bend_type))

        if stub_length > 0:
            sl = stub_length
            ml2 = self.total_length / 2
            for idx in (0, 1):
                s = _arm_sign(idx)
                ww = self._arm_width(idx)
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
    # Device frame: x=0 at the coupling-region centre. Same bend/coupling
    # layout as DirectionalCoupler — see coupler.py's comment block — the
    # only difference is that arm 0 ("bus") and arm 1 ("access") each carry
    # their own width and near-pitch offset instead of a single shared one.

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
        section of the whole device, and the other arm really is sitting
        right there in it. Same convention as DirectionalCoupler.
        """
        (in0, in1), (out0, out1) = self._bend_bounds()
        p_far = self.waveguide_offset / 2

        lo, hi = x0 - stub_left, x1 + stub_right
        self._check_cut(lo)
        self._check_cut(hi)

        polys: List[List[Tuple[float, float]]] = []

        def _rect(xa, xb, yc, ww):
            if xb - xa > 1e-9:
                polys.append([(xa, yc-ww/2), (xb, yc-ww/2), (xb, yc+ww/2), (xa, yc+ww/2)])

        # left access waveguide (both arms, before the bend)
        if lo < in0:
            for idx in (0, 1):
                s = _arm_sign(idx)
                _rect(lo, min(hi, in0), s * p_far, self._arm_width(idx))
        # left S-bend (both arms, all-or-nothing - never cut mid-bend)
        if lo <= in0 and hi >= in1:
            for idx in (0, 1):
                s = _arm_sign(idx)
                polys.extend(_sbend_polygon(
                    in0, s*p_far, in1, s*self._p_near(idx),
                    self._arm_width(idx), self.bend_type,
                ))
        # coupling region (both arms, arbitrary sub-range)
        body_lo, body_hi = max(lo, in1), min(hi, out0)
        if body_hi - body_lo > 1e-9:
            for idx in (0, 1):
                s = _arm_sign(idx)
                _rect(body_lo, body_hi, s * self._p_near(idx), self._arm_width(idx))
        # right S-bend (both arms)
        if lo <= out0 and hi >= out1:
            for idx in (0, 1):
                s = _arm_sign(idx)
                polys.extend(_sbend_polygon(
                    out0, s*self._p_near(idx), out1, s*p_far,
                    self._arm_width(idx), self.bend_type,
                ))
        # right access waveguide (both arms)
        if hi > out1:
            for idx in (0, 1):
                s = _arm_sign(idx)
                _rect(max(lo, out1), hi, s * p_far, self._arm_width(idx))

        return polys

    def section_prisms(self, x0: float, x1: float, si_medium,
                        stub_left: float = 0.0, stub_right: float = 0.0,
                        input_idx: int = 0, output_idx: int = 0):
        import meep as mp
        return [
            mp.Prism(
                vertices=[mp.Vector3(x, y, 0) for x, y in verts],
                height=mp.inf, material=si_medium,
            )
            for verts in self._section_polygons(x0, x1, stub_left, stub_right, input_idx, output_idx)
        ]

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
            return self.bus_width + self.access_width + self.gap  # full coupled-pair span
        idx = input_idx if x < in1 else output_idx
        return self._arm_width(idx)

    def center_at(self, x: float, input_idx: int = 0, output_idx: int = 0) -> float:
        (in0, in1), (out0, out1) = self._bend_bounds()
        if in1 <= x <= out0:
            return 0.0
        p_far = self.waveguide_offset / 2
        if x < in1:
            idx = input_idx
            s = _arm_sign(idx)
            p_near = self._p_near(idx)
            if x <= in0:
                return s * p_far
            frac = (x - in0) / (in1 - in0)  # linear approx; only ever evaluated eps from a boundary
            return s * (p_far + frac * (p_near - p_far))
        idx = output_idx
        s = _arm_sign(idx)
        p_near = self._p_near(idx)
        if x >= out1:
            return s * p_far
        frac = (x - out0) / (out1 - out0)
        return s * (p_near + frac * (p_far - p_near))

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.bend_length, self.coupling_length,
                self.bend_length, straight_length]

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.total_length, self.waveguide_offset
