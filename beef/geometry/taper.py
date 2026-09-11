"""Linear-width waveguide taper — a standalone DeviceGeometry, usable on its
own or as the building block MMI composes for each port."""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from .base import DeviceGeometry, _oxide_box, polygons_to_prisms

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


class Taper(DeviceGeometry):
    """Straight-line linear width transition, centered at x=0:
    width_start at x=-length/2, width_end at x=+length/2, linear between.

    The centerline may ALSO be slanted: `y_start` at x=-length/2, `y_end`
    at x=+length/2, linear between (both default 0.0 — a straight,
    non-slanted taper, byte-identical to every use of this class before
    2026-09-09). Added for `MMIPort.y_far`
    (modal_multiplexer_validation.md gap 9 item 2): the fabricated seed
    mux's top branch taper measures MMI-side end at y=+1.0375 but
    waveguide-side end at y=+1.0480 — a real, if sub-pixel-at-res-40,
    asymmetry this repo's geometry model could not previously express at
    all (every taper was forced fully straight).

    The cross-section varies everywhere along its length, so
    section_prisms()/section_gds() must include it whole or not at all.
    """

    def __init__(
        self,
        length: float,
        width_start: float,
        width_end: float,
        waveguide_height: float = 0.16,
        y_start: float = 0.0,
        y_end: float = 0.0,
    ):
        self.length = length
        self.width_start = width_start
        self.width_end = width_end
        self.waveguide_height = waveguide_height
        self.y_start = y_start
        self.y_end = y_end

    @property
    def total_length(self) -> float:
        return self.length

    def polygons(self, stub_length: float = 0.0) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        L = self.length
        wa, wb = self.width_start, self.width_end
        ya, yb = self.y_start, self.y_end
        wg = [[(-L/2, ya-wa/2), (L/2, yb-wb/2), (L/2, yb+wb/2), (-L/2, ya+wa/2)]]
        if stub_length > 0:
            sl = stub_length
            wg.append([(-L/2-sl, ya-wa/2), (-L/2, ya-wa/2), (-L/2, ya+wa/2), (-L/2-sl, ya+wa/2)])
            wg.append([(L/2, yb-wb/2), (L/2+sl, yb-wb/2), (L/2+sl, yb+wb/2), (L/2, yb+wb/2)])
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
                f"Cut/stub at x={x:.4f} falls strictly inside a Taper "
                "(cross-section varies everywhere along its length) — a "
                "Taper must be included whole or not at all in any section."
            )

    def _section_polygons(self, x0: float, x1: float,
                           stub_left: float = 0.0, stub_right: float = 0.0):
        L = self.length
        wa, wb = self.width_start, self.width_end
        ya, yb = self.y_start, self.y_end
        lo, hi = x0 - stub_left, x1 + stub_right
        self._check_cut(lo)
        self._check_cut(hi)

        polys: List[List[Tuple[float, float]]] = []

        def _rect(xa, xb, w, y):
            if xb - xa > 1e-9:
                polys.append([(xa, y-w/2), (xb, y-w/2), (xb, y+w/2), (xa, y+w/2)])

        if lo < -L/2:
            _rect(lo, min(hi, -L/2), wa, ya)
        if lo <= -L/2 and hi >= L/2:
            polys.append([(-L/2, ya-wa/2), (L/2, yb-wb/2), (L/2, yb+wb/2), (-L/2, ya+wa/2)])
        if hi > L/2:
            _rect(max(lo, L/2), hi, wb, yb)
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
        L = self.length
        if x <= -L/2:
            return self.width_start
        if x >= L/2:
            return self.width_end
        frac = (x + L/2) / L
        return self.width_start + frac * (self.width_end - self.width_start)

    def center_at(self, x: float, **_kwargs) -> float:
        L = self.length
        if x <= -L/2:
            return self.y_start
        if x >= L/2:
            return self.y_end
        frac = (x + L/2) / L
        return self.y_start + frac * (self.y_end - self.y_start)

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.length, straight_length]

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.length, max(self.width_start, self.width_end)
