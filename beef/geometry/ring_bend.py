# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Georgia Tech Research Corporation
"""RingResonatorBendSection: standalone 90-degree bend, x-normal in /
y-normal out, the loss/phase-reference sub-simulation reused by mirror
symmetry for both far corners of a racetrack ring.

Reproduces examples/ring_resonator/python/ring_resonator_hybrid_cascaded.py's
build_bend_polys() (identical in ring_resonator_hybrid.py), minus that
script's own baked-in access stub lengths — this class's own x-range is the
bare arc [0, bend_radius], matching every other DeviceGeometry's convention
that stub_left/stub_right (not the geometry class itself) supply access-
waveguide padding.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from ._arc import arc_polygon
from .base import DeviceGeometry, _oxide_box, polygons_to_prisms

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


@dataclass(frozen=True)
class BendPortGeometry:
    """Local geometry of the y-normal output port RingResonatorBendSection
    draws — see RingResonatorCouplerSection.RingPortGeometry for the same
    pattern and how a caller turns this into a
    beef.solvers.fdtd.PortSpec + y0/y1.

    Attributes:
        x: device-frame x of the output stub column.
        width: output aperture span (µm).
        y_boundary: device-frame y up to which section_prisms() draws the
            output stub (before a caller's own y_pml padding).
        sign: +1 — the drawn material sits ABOVE y_boundary, so this port's
            own launch/measurement direction is toward -y, into the domain
            from above (matches PortSpec.sign's meaning directly, in
            contrast with the coupler section's ring port).
    """
    x: float
    width: float
    y_boundary: float
    sign: int


class RingResonatorBendSection(DeviceGeometry):
    """x-normal input (straight run ending at x=0) -> 90-degree arc (radius
    bend_radius, center (0, bend_radius)) -> y-normal output stub at
    x=bend_radius. Both ports are explicit: input via the ABC's ordinary
    x0/x1 contract, output via out_port() — width_at()/center_at() cannot
    describe the output side (same limitation as
    RingResonatorCouplerSection.ring_port(), see that class's docstring).

    Args:
        waveguide_width: strip width (µm). Default 0.45.
        bend_radius: bend centerline radius (µm). Default 5.0.
        output_margin: reference-plane offset past the arc's own end where
            the output port sits, matching RingResonatorCouplerSection's
            ring_stub_margin. Default 0.5.
        waveguide_height: Si layer thickness (µm). Default 0.22.
    """

    def __init__(
        self,
        waveguide_width: float = 0.45,
        bend_radius: float = 5.0,
        output_margin: float = 0.5,
        waveguide_height: float = 0.22,
    ):
        self.waveguide_width = waveguide_width
        self.bend_radius = bend_radius
        self.output_margin = output_margin
        self.waveguide_height = waveguide_height

    @property
    def total_length(self) -> float:
        return self.bend_radius

    def out_port(self) -> BendPortGeometry:
        return BendPortGeometry(
            x=self.bend_radius, width=3.0,
            y_boundary=self.bend_radius + self.output_margin, sign=1,
        )

    def _check_cut(self, x: float) -> None:
        if 0.0 < x < self.bend_radius:
            raise ValueError(
                f"Cut/stub at x={x:.4f} falls strictly inside the bend arc "
                "(cross-section is not a function of x there) - cuts must "
                "land at x=0 (arc start) or x=bend_radius (arc end)."
            )

    def _section_polygons(self, x0: float, x1: float,
                           stub_left: float = 0.0, stub_right: float = 0.0,
                           y_pml: float = 0.0) -> List[List[Tuple[float, float]]]:
        ww = self.waveguide_width
        R = self.bend_radius
        lo, hi = x0 - stub_left, x1 + stub_right
        self._check_cut(lo)
        self._check_cut(hi)

        polys: List[List[Tuple[float, float]]] = []

        def _rect(xa, xb, yc, w):
            if xb - xa > 1e-9:
                polys.append([(xa, yc - w/2), (xb, yc - w/2), (xb, yc + w/2), (xa, yc + w/2)])

        # input straight, before the arc
        _rect(lo, min(hi, 0.0), 0.0, ww)
        # the arc itself + output stub: all-or-nothing
        if lo <= 0.0 and hi >= R:
            polys.append(arc_polygon(0.0, R, R - ww/2, R + ww/2, -90.0, 0.0))
            port = self.out_port()
            y_top = port.y_boundary + y_pml
            polys.append([(R - ww/2, R), (R + ww/2, R), (R + ww/2, y_top), (R - ww/2, y_top)])

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
        return self.waveguide_width

    def center_at(self, x: float, **_kwargs) -> float:
        return 0.0

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.bend_radius, straight_length]

    def to_gds_eme(self, straight_length: float, oxide_width: float,
                    x_start: float = 0.0, **kwargs):
        raise NotImplementedError(
            "RingResonatorBendSection.to_gds_eme() is not implemented: this "
            "class is an FDTD-only loss/phase reference sub-simulation, never "
            "part of a whole-device EME cell layout."
        )

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.total_length, self.out_port().y_boundary
