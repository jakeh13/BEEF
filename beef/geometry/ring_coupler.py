# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Georgia Tech Research Corporation
"""RingResonatorCouplerSection: bus + one near-bend + ring stub, the
FDTD-relevant piece of a racetrack ring resonator's coupling region.

Reproduces examples/ring_resonator/python/ring_resonator_hybrid_cascaded.py's
build_left_polys() exactly (side="left"; side="right" is its analytically-
verified x-mirror — see that script's own module docstring for the
theta-angle derivation reused here). The body (straight bus+ring pair) and
the far bends are deliberately NOT covered by this class — see module
docstring below and beef.geometry.ring_bend for the standalone bend piece.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from ._arc import arc_polygon
from .base import DeviceGeometry, _oxide_box, polygons_to_prisms

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


@dataclass(frozen=True)
class RingPortGeometry:
    """Local geometry of the ring-stub port RingResonatorCouplerSection
    draws, solver-agnostic (plain floats, no meep/beef.solvers import — see
    beef.solvers.fdtd.PortSpec, which a caller builds from these fields
    directly: PortSpec(axis="y", sign=sign, center=x, width=width), with
    y0/y1 spanning up to y_boundary on whichever side the stub's own root
    (at the near bend) sits.

    Attributes:
        x: device-frame x of the stub column (transverse center for the
            y-normal port).
        width: stub aperture span (µm), matching PORT_APERTURE in the
            reference scripts.
        y_boundary: device-frame y up to which section_prisms() draws the
            stub — the y0 or y1 a caller passes to solve_section() on this
            side (before that solver's own dpml padding).
        sign: -1 for this device (a stub whose drawn material sits BELOW
            y_boundary, so its own port launches/measures toward -y, into
            the domain) — matches PortSpec.sign's meaning directly.
    """
    x: float
    width: float
    y_boundary: float
    sign: int


class RingResonatorCouplerSection(DeviceGeometry):
    """Bus (straight, y=0) + one near-bend (radius bend_radius) + ring stub,
    for one side ("left" or "right") of a racetrack ring's coupling region.

    Device frame, x=0 at the coupling-region centre (mirrors
    DirectionalCoupler's convention). x spans [outer_x, iface_x] (or the
    reverse for side="right") — outer_x is the true physical bus port
    (where an external caller's own stub_left/waveguide_length padding
    would be added, matching HybridSolver._stubs()' convention: this class
    does NOT bake extra access-waveguide length into its own geometry the
    way the monolithic reference script's INNER_X_HALF does), iface_x is
    the interior interface cut just past the near bend, matching
    IFACE_BUFFER past the bend's own endpoint.

    x0 must be < x1 in every section_prisms()/section_gds() call, in the
    usual ascending device-frame sense — for side="right" that means
    x0=iface_x (the smaller value) and x1=outer_x (the larger), the
    reverse pairing from side="left". Not defensively checked (matches
    DirectionalCoupler's own section_prisms(), which relies on the same
    implicit contract); a reversed pair silently draws nothing rather than
    raising, since every region guard below is an x0-to-x1 overlap test.

    section_prisms()/section_gds()/width_at()/center_at() describe ONLY
    this x-range along the shared x-axis; they say nothing about the ring
    stub once it leaves the x-axis via the near-bend — width_at()'s x-in/
    y-out contract cannot honestly describe a fixed-x slice through a
    curved wedge (not a function of x there, the same reason
    SBend._check_cut/DirectionalCoupler._check_cut forbid an interior bend
    cut at all). Use ring_port() for the stub's own geometry instead.

    Args:
        waveguide_width: single-mode strip width (µm). Default 0.45.
        gap: bus<->ring edge-to-edge coupling gap (µm). Default 0.25.
        bend_radius: near-bend centerline radius (µm). Default 5.0.
        near_length: coupling-region straight length, L_C (µm). Default 4.0.
        ring_stub_margin: ring-port reference-plane offset above the near
            bend's arc end, S_MON (µm). Default 0.5.
        iface_buffer: straight buffer past the near bend before the
            interface cut, IFACE_BUFFER (µm). Default 1.0.
        bus_clearance: extra bus length beyond the stub column, before its
            own outer (true port) end — matches the reference script's
            Y_BELOW, reused there for both this and the y-domain's own
            lower bound. Default 2.0.
        waveguide_height: Si layer thickness (µm). Default 0.22.
        side: "left" (default) or "right" — selects which near-bend/stub
            (and hence which mirror sign) this instance draws.
    """

    def __init__(
        self,
        waveguide_width: float = 0.45,
        gap: float = 0.25,
        bend_radius: float = 5.0,
        near_length: float = 4.0,
        ring_stub_margin: float = 0.5,
        iface_buffer: float = 1.0,
        bus_clearance: float = 2.0,
        waveguide_height: float = 0.22,
        side: str = "left",
    ):
        if side not in ("left", "right"):
            raise ValueError(f"side must be 'left' or 'right', got {side!r}.")
        if not (0.0 < iface_buffer < near_length / 2):
            raise ValueError(
                f"iface_buffer ({iface_buffer}) must be in (0, near_length/2={near_length/2}) "
                "so the interface cut lands strictly between the coupling centre and the bend."
            )
        self.waveguide_width = waveguide_width
        self.gap = gap
        self.bend_radius = bend_radius
        self.near_length = near_length
        self.ring_stub_margin = ring_stub_margin
        self.iface_buffer = iface_buffer
        self.bus_clearance = bus_clearance
        self.waveguide_height = waveguide_height
        self.side = side
        self._s = -1.0 if side == "left" else 1.0   # sign mapping canonical magnitudes -> device frame

    # ------------------------------------------------------------------
    # Geometry breakpoints (device frame)
    # ------------------------------------------------------------------

    @property
    def y_ring(self) -> float:
        """Ring near-straight centerline y (bus is at y=0)."""
        return self.waveguide_width + self.gap

    @property
    def hc(self) -> float:
        """x where the near-bend meets the coupling-region straight (the
        DirectionalCoupler-style coupling-region half-width)."""
        return self._s * self.near_length / 2

    @property
    def stub_x(self) -> float:
        """x of the stub column / bend's outer end."""
        return self._s * (self.near_length / 2 + self.bend_radius)

    @property
    def iface_x(self) -> float:
        """Interior interface cut x, iface_buffer past the near bend."""
        return self._s * (self.near_length / 2 - self.iface_buffer)

    @property
    def outer_x(self) -> float:
        """True physical bus port x (this class's own x0/x1 boundary)."""
        return self._s * (self.near_length / 2 + self.bend_radius + self.bus_clearance)

    @property
    def total_length(self) -> float:
        return abs(self.iface_x - self.outer_x)

    def ring_port(self) -> RingPortGeometry:
        y_port = self.y_ring + self.bend_radius + self.ring_stub_margin
        return RingPortGeometry(
            x=self.stub_x, width=3.0, y_boundary=y_port, sign=-1,
        )

    def _check_cut(self, x: float) -> None:
        bend_lo, bend_hi = sorted((self.stub_x, self.hc))
        if bend_lo < x < bend_hi:
            raise ValueError(
                f"Cut/stub at x={x:.4f} falls strictly inside the near-bend arc "
                "(cross-section is not a function of x there) - cuts must land "
                "on the arc endpoint (stub_x or hc) or in the uniform bus/"
                "coupled-pair region on either side."
            )

    # ------------------------------------------------------------------
    # Arbitrary-section interface (HybridSolver)
    # ------------------------------------------------------------------

    def _section_polygons(self, x0: float, x1: float,
                           stub_left: float = 0.0, stub_right: float = 0.0,
                           y_pml: float = 0.0) -> List[List[Tuple[float, float]]]:
        ww = self.waveguide_width
        y_ring = self.y_ring
        R = self.bend_radius

        lo, hi = x0 - stub_left, x1 + stub_right
        self._check_cut(lo)
        self._check_cut(hi)

        bend_lo, bend_hi = sorted((self.stub_x, self.hc))

        # The caller's stub_left/stub_right padding continues the waveguide
        # OUTWARD past this class's own [outer_x, iface_x] range, exactly as
        # DirectionalCoupler._section_polygons() draws its access arms out to
        # lo/hi and ring_bend.RingBendSection draws its input straight from lo.
        # Both ends are legitimately uniform: past outer_x is bare access bus,
        # and past iface_x is the mirrored section, where the same bus+sliver
        # pair continues straight. Clamping to outer_x/iface_x instead (as this
        # method did until 2026-08-02) left the padded domain filled with bare
        # cladding, so an eigenmode source or monitor placed on the padded
        # boundary -- which is where solve_section() puts them, at
        # lo + source_margin and hi - monitor_margin, NOT at x0/x1 -- sat on a
        # cross-section with no waveguide in it at all.
        #
        # The hc end of the sliver is the one boundary that must still clamp:
        # past it is the bend, not more straight guide.
        seam_out = max(self.iface_x, hi) if self.iface_x > self.hc else min(self.iface_x, lo)
        sliver_lo, sliver_hi = sorted((self.hc, seam_out))

        polys: List[List[Tuple[float, float]]] = []

        def _rect(xa, xb, yc, w):
            if xb - xa > 1e-9:
                polys.append([(xa, yc - w/2), (xb, yc - w/2), (xb, yc + w/2), (xa, yc + w/2)])

        # bus (continuous the whole way, y=0, padding included)
        _rect(lo, hi, 0.0, ww)

        # near-bend + stub column: all-or-nothing, same guard as _check_cut
        if lo <= bend_lo and hi >= bend_hi:
            theta0, theta1 = (180.0, 270.0) if self.side == "left" else (270.0, 360.0)
            polys.append(arc_polygon(self.hc, y_ring + R, R - ww/2, R + ww/2, theta0, theta1))
            port = self.ring_port()
            y_top = port.y_boundary + y_pml
            _rect_y = [(self.stub_x - ww/2, y_ring + R), (self.stub_x + ww/2, y_ring + R),
                       (self.stub_x + ww/2, y_top), (self.stub_x - ww/2, y_top)]
            polys.append(_rect_y)

        # near-straight ring sliver (between the bend's inner end and the interface cut)
        _rect(max(lo, sliver_lo), min(hi, sliver_hi), y_ring, ww)

        return polys

    def section_prisms(self, x0: float, x1: float, si_medium,
                        stub_left: float = 0.0, stub_right: float = 0.0,
                        y_pml: float = 0.0, wg_height: Optional[float] = None,
                        **_kwargs):
        return polygons_to_prisms(
            self._section_polygons(x0, x1, stub_left, stub_right, y_pml),
            si_medium, wg_height=wg_height,
        )

    def section_gds(self, x0: float, x1: float, oxide_width: float,
                     stub_left: float = 0.0, stub_right: float = 0.0,
                     y_pml: float = 0.0, **_kwargs
                     ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        lo, hi = x0 - stub_left, x1 + stub_right
        x_off = -lo
        wg_polys = [
            [(x + x_off, y) for x, y in verts]
            for verts in self._section_polygons(x0, x1, stub_left, stub_right, y_pml)
        ]
        ox_poly = _oxide_box(0.0, hi - lo, oxide_width)
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    def width_at(self, x: float, **_kwargs) -> float:
        self._check_cut(x)
        sliver_lo, sliver_hi = sorted((self.hc, self.iface_x))
        if sliver_lo <= x <= sliver_hi:
            return 2 * self.waveguide_width + self.gap  # bus + ring sliver, coupled-pair span
        return self.waveguide_width

    def center_at(self, x: float, **_kwargs) -> float:
        self._check_cut(x)
        sliver_lo, sliver_hi = sorted((self.hc, self.iface_x))
        return self.y_ring / 2 if sliver_lo <= x <= sliver_hi else 0.0

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.bend_radius, self.iface_buffer]

    def to_gds_eme(self, straight_length: float, oxide_width: float,
                    x_start: float = 0.0, **kwargs):
        raise NotImplementedError(
            "RingResonatorCouplerSection.to_gds_eme() is not implemented: this "
            "class covers only the mixed-orientation coupler section, never a "
            "whole-device EME cell layout — the straight bus+ring 'body' cross-"
            "section is built directly (raw polygons + mw.compute_modes_meep), "
            "not through a DeviceGeometry, matching "
            "ring_resonator_hybrid_cascaded.py's solve_body_supermodes()."
        )

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.total_length, self.ring_port().y_boundary
