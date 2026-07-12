from __future__ import annotations

import itertools
from typing import Dict, List, Optional, Tuple

import gdsfactory as gf

from .base import DeviceGeometry

_ctr = itertools.count()

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


def _reset_gf():
    gf.gpdk.PDK.activate()


def _component(name: str) -> gf.Component:
    return gf.Component(f"{name}_{next(_ctr)}")


def _add_oxide(c: gf.Component, x_start: float, total_x: float, oxide_width: float) -> None:
    c.add_polygon([
        (x_start, -oxide_width / 2),
        (total_x,  -oxide_width / 2),
        (total_x,   oxide_width / 2),
        (x_start,  oxide_width / 2),
    ], layer=_OX_LAYER)


# ---------------------------------------------------------------------------
# 1×2 MMI
# ---------------------------------------------------------------------------

class MMI1x2(DeviceGeometry):
    """1×2 multimode interferometer splitter.

    waveguide_offset is the center-to-center distance between the two output
    waveguides.
    """

    def __init__(
        self,
        mm_length: float,
        mm_width: float,
        waveguide_width: float,
        taper_length: float,
        taper_width: float,
        waveguide_offset: float,
        waveguide_height: float = 0.16,
    ):
        self.mm_length = mm_length
        self.mm_width = mm_width
        self.waveguide_width = waveguide_width
        self.taper_length = taper_length
        self.taper_width = taper_width
        self.waveguide_offset = waveguide_offset
        self.waveguide_height = waveguide_height

    @property
    def total_length(self) -> float:
        return self.mm_length + 2 * self.taper_length

    def polygons(
        self, stub_length: float = 0.0
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """1×2 MMI polygons in device frame (center at x=0)."""
        L, W = self.mm_length, self.mm_width
        tl, tw = self.taper_length, self.taper_width
        ww = self.waveguide_width
        wo = self.waveguide_offset / 2
        ml = self.total_length  # L + 2*tl

        wg = [
            # MMI body
            [(-L/2, -W/2), (L/2, -W/2), (L/2, W/2), (-L/2, W/2)],
            # Input taper
            [(-L/2-tl, -ww/2), (-L/2, -tw/2), (-L/2, tw/2), (-L/2-tl, ww/2)],
            # Output top taper
            [(L/2, -tw/2+wo), (L/2+tl, -ww/2+wo), (L/2+tl, ww/2+wo), (L/2, tw/2+wo)],
            # Output bottom taper
            [(L/2, -tw/2-wo), (L/2+tl, -ww/2-wo), (L/2+tl, ww/2-wo), (L/2, tw/2-wo)],
        ]
        if stub_length > 0:
            sl = stub_length
            wg.append([(-ml/2-sl, -ww/2), (-ml/2, -ww/2), (-ml/2, ww/2), (-ml/2-sl, ww/2)])
            for sign in (+1, -1):
                wg.append([
                    (ml/2,    -ww/2 + sign*wo), (ml/2+sl, -ww/2 + sign*wo),
                    (ml/2+sl,  ww/2 + sign*wo), (ml/2,     ww/2 + sign*wo),
                ])
        return {_WG_LAYER: wg}

    def to_gds_eme(
        self, straight_length: float, oxide_width: float,
        x_start: float = 0.0, **kwargs
    ) -> gf.Component:
        _reset_gf()
        ml = self.total_length
        x_off = x_start + ml / 2 + straight_length
        c = _component(f"mmi_1x2_eme_{id(self)}")
        for verts in self.polygons(stub_length=straight_length)[_WG_LAYER]:
            c.add_polygon([(x + x_off, y) for x, y in verts], layer=_WG_LAYER)
        _add_oxide(c, x_start, x_start + ml + 2 * straight_length, oxide_width)
        return c

    # ------------------------------------------------------------------
    # Arbitrary-section interface (HybridSolver)
    # ------------------------------------------------------------------
    #
    # Device frame: x=0 at the MMI body centre (same as polygons()).
    #   x in [-half_body-tl, -half_body] : input taper (width ww -> tw)
    #   x in [-half_body, +half_body]    : MMI body (uniform, width W)
    #   x in [+half_body, +half_body+tl] : output taper (width tw -> ww,
    #                                      offset to y=+wo or y=-wo per branch)
    # Anything beyond +-(half_body+tl) is the (virtual, arbitrary-length)
    # access waveguide, uniform width ww.
    #
    # Cuts / stubs may land anywhere except strictly inside a taper (its
    # cross-section is not constant there). port_idx (0=top, 1=bot) selects
    # which output branch to include once a range reaches the output taper
    # or access waveguide; irrelevant (ignored) otherwise.

    def _taper_bounds(self):
        half_body = self.mm_length / 2
        tl = self.taper_length
        return (-half_body - tl, -half_body), (half_body, half_body + tl)

    def _check_cut(self, x: float) -> None:
        (in0, in1), (out0, out1) = self._taper_bounds()
        if in0 < x < in1 or out0 < x < out1:
            raise ValueError(
                f"Cut/stub at x={x:.4f} falls strictly inside a taper "
                "(cross-section varies there) — cuts must land on a taper "
                "endpoint or inside the uniform MMI body / access waveguide."
            )

    def _section_polygons(self, x0: float, x1: float,
                           stub_left: float = 0.0, stub_right: float = 0.0,
                           port_idx: Optional[int] = None):
        """Raw vertex lists (device frame) covering [x0-stub_left, x1+stub_right]."""
        (in0, in1), (out0, out1) = self._taper_bounds()
        ww, tw, W = self.waveguide_width, self.taper_width, self.mm_width
        wo = self.waveguide_offset / 2
        y_c = 0.0 if port_idx is None else (+wo if port_idx == 0 else -wo)

        lo, hi = x0 - stub_left, x1 + stub_right
        self._check_cut(lo)
        self._check_cut(hi)

        polys = []

        def _rect(xa, xb, half_w, yc=0.0):
            if xb - xa > 1e-9:
                polys.append([(xa, -half_w + yc), (xb, -half_w + yc),
                              (xb, half_w + yc), (xa, half_w + yc)])

        # left access waveguide (uniform, width ww)
        if lo < in0:
            _rect(lo, min(hi, in0), ww / 2)
        # input taper (all-or-nothing — never cut mid-taper)
        if lo <= in0 and hi >= in1:
            polys.append([(in0, -ww/2), (in1, -tw/2), (in1, tw/2), (in0, ww/2)])
        # MMI body (uniform, arbitrary sub-range)
        _rect(max(lo, in1), min(hi, out0), W / 2)
        # output taper (all-or-nothing, offset to the selected branch)
        if lo <= out0 and hi >= out1:
            polys.append([(out0, -tw/2 + y_c), (out1, -ww/2 + y_c),
                          (out1, ww/2 + y_c), (out0, tw/2 + y_c)])
        # right access waveguide (uniform, width ww, offset to the branch)
        if hi > out1:
            _rect(max(lo, out1), hi, ww / 2, yc=y_c)

        return polys

    def section_prisms(self, x0: float, x1: float, si_medium,
                        stub_left: float = 0.0, stub_right: float = 0.0,
                        port_idx: Optional[int] = None):
        import meep as mp
        return [
            mp.Prism(
                vertices=[mp.Vector3(x, y, 0) for x, y in verts],
                height=mp.inf, material=si_medium,
            )
            for verts in self._section_polygons(x0, x1, stub_left, stub_right, port_idx)
        ]

    def section_gds(self, x0: float, x1: float, oxide_width: float,
                     stub_left: float = 0.0, stub_right: float = 0.0,
                     port_idx: Optional[int] = None) -> gf.Component:
        _reset_gf()
        lo, hi = x0 - stub_left, x1 + stub_right
        x_off = -lo
        c = _component(f"mmi_1x2_section_{id(self)}")
        for verts in self._section_polygons(x0, x1, stub_left, stub_right, port_idx):
            c.add_polygon([(x + x_off, y) for x, y in verts], layer=_WG_LAYER)
        _add_oxide(c, 0.0, hi - lo, oxide_width)
        return c

    def width_at(self, x: float, port_idx: Optional[int] = None) -> float:
        (in0, in1), (out0, out1) = self._taper_bounds()
        ww, tw, W = self.waveguide_width, self.taper_width, self.mm_width
        if in1 <= x <= out0:
            return W
        if in0 <= x <= in1:
            frac = (x - in0) / self.taper_length
            return ww + frac * (tw - ww)
        if out0 <= x <= out1:
            frac = (x - out0) / self.taper_length
            return tw + frac * (ww - tw)
        return ww

    def center_at(self, x: float, port_idx: Optional[int] = None) -> float:
        half_body = self.mm_length / 2
        if x <= half_body or port_idx is None:
            return 0.0
        wo = self.waveguide_offset / 2
        return +wo if port_idx == 0 else -wo

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.taper_length, self.mm_length,
                self.taper_length, straight_length]

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.total_length, self.mm_width


# ---------------------------------------------------------------------------
# 2×2 MMI
# ---------------------------------------------------------------------------

class MMI2x2(DeviceGeometry):
    """2×2 multimode interferometer splitter.

    waveguide_offset is the full center-to-center distance between the two
    input (and two output) waveguides.
    """

    def __init__(
        self,
        mm_length: float,
        mm_width: float,
        waveguide_width: float,
        taper_length: float,
        taper_width: float,
        waveguide_offset: float,
        waveguide_height: float = 0.16,
    ):
        self.mm_length = mm_length
        self.mm_width = mm_width
        self.waveguide_width = waveguide_width
        self.taper_length = taper_length
        self.taper_width = taper_width
        self.waveguide_offset = waveguide_offset
        self.waveguide_height = waveguide_height

    @property
    def total_length(self) -> float:
        return self.mm_length + 2 * self.taper_length

    def polygons(
        self, stub_length: float = 0.0
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """2×2 MMI polygons in device frame (center at x=0)."""
        L, W = self.mm_length, self.mm_width
        tl, tw = self.taper_length, self.taper_width
        ww = self.waveguide_width
        wo = self.waveguide_offset / 2
        ml = self.total_length

        wg = [
            # MMI body
            [(-L/2, -W/2), (L/2, -W/2), (L/2, W/2), (-L/2, W/2)],
            # Input tapers (top + bottom)
            [(-L/2-tl, -ww/2+wo), (-L/2, -tw/2+wo), (-L/2, tw/2+wo), (-L/2-tl, ww/2+wo)],
            [(-L/2-tl, -ww/2-wo), (-L/2, -tw/2-wo), (-L/2, tw/2-wo), (-L/2-tl, ww/2-wo)],
            # Output tapers (top + bottom)
            [(L/2, -tw/2+wo), (L/2+tl, -ww/2+wo), (L/2+tl, ww/2+wo), (L/2, tw/2+wo)],
            [(L/2, -tw/2-wo), (L/2+tl, -ww/2-wo), (L/2+tl, ww/2-wo), (L/2, tw/2-wo)],
        ]
        if stub_length > 0:
            sl = stub_length
            for sign in (+1, -1):
                wg.append([
                    (-ml/2-sl, -ww/2 + sign*wo), (-ml/2, -ww/2 + sign*wo),
                    (-ml/2,     ww/2 + sign*wo), (-ml/2-sl, ww/2 + sign*wo),
                ])
                wg.append([
                    (ml/2,    -ww/2 + sign*wo), (ml/2+sl, -ww/2 + sign*wo),
                    (ml/2+sl,  ww/2 + sign*wo), (ml/2,     ww/2 + sign*wo),
                ])
        return {_WG_LAYER: wg}

    def to_gds_eme(
        self, straight_length: float, oxide_width: float,
        x_start: float = 0.0, **kwargs
    ) -> gf.Component:
        _reset_gf()
        ml = self.total_length
        x_off = x_start + ml / 2 + straight_length
        c = _component(f"mmi_2x2_eme_{id(self)}")
        for verts in self.polygons(stub_length=straight_length)[_WG_LAYER]:
            c.add_polygon([(x + x_off, y) for x, y in verts], layer=_WG_LAYER)
        _add_oxide(c, x_start, x_start + ml + 2 * straight_length, oxide_width)
        return c

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.taper_length, self.mm_length,
                self.taper_length, straight_length]

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.total_length, self.mm_width


# ---------------------------------------------------------------------------
# 1×4 MMI
# ---------------------------------------------------------------------------

class MMI1x4(DeviceGeometry):
    """1×4 multimode interferometer splitter.

    waveguide_offset is the center-to-center spacing between adjacent output
    waveguides; the four outputs are at ±wo/2 and ±3wo/2.
    """

    def __init__(
        self,
        mm_length: float,
        mm_width: float,
        waveguide_width: float,
        taper_length: float,
        taper_width: float,
        waveguide_offset: float,
        waveguide_height: float = 0.16,
    ):
        self.mm_length = mm_length
        self.mm_width = mm_width
        self.waveguide_width = waveguide_width
        self.taper_length = taper_length
        self.taper_width = taper_width
        self.waveguide_offset = waveguide_offset
        self.waveguide_height = waveguide_height

    @property
    def total_length(self) -> float:
        return self.mm_length + 2 * self.taper_length

    def polygons(
        self, stub_length: float = 0.0
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """1×4 MMI polygons in device frame (center at x=0)."""
        L, W = self.mm_length, self.mm_width
        tl, tw = self.taper_length, self.taper_width
        ww = self.waveguide_width
        wo = self.waveguide_offset  # adjacent spacing; outputs at ±wo/2 and ±3wo/2
        ml = self.total_length
        out_ys = [3*wo/2, wo/2, -wo/2, -3*wo/2]

        wg = [
            [(-L/2, -W/2), (L/2, -W/2), (L/2, W/2), (-L/2, W/2)],
            [(-L/2-tl, -ww/2), (-L/2, -tw/2), (-L/2, tw/2), (-L/2-tl, ww/2)],
        ]
        for yo in out_ys:
            wg.append([(L/2, -tw/2+yo), (L/2+tl, -ww/2+yo), (L/2+tl, ww/2+yo), (L/2, tw/2+yo)])

        if stub_length > 0:
            sl = stub_length
            wg.append([(-ml/2-sl, -ww/2), (-ml/2, -ww/2), (-ml/2, ww/2), (-ml/2-sl, ww/2)])
            for yo in out_ys:
                wg.append([
                    (ml/2, -ww/2+yo), (ml/2+sl, -ww/2+yo),
                    (ml/2+sl, ww/2+yo), (ml/2, ww/2+yo),
                ])
        return {_WG_LAYER: wg}

    def to_gds_eme(
        self, straight_length: float, oxide_width: float,
        x_start: float = 0.0, **kwargs
    ) -> gf.Component:
        _reset_gf()
        ml = self.total_length
        x_off = x_start + ml / 2 + straight_length
        c = _component(f"mmi_1x4_eme_{id(self)}")
        for verts in self.polygons(stub_length=straight_length)[_WG_LAYER]:
            c.add_polygon([(x + x_off, y) for x, y in verts], layer=_WG_LAYER)
        _add_oxide(c, x_start, x_start + ml + 2 * straight_length, oxide_width)
        return c

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.taper_length, self.mm_length,
                self.taper_length, straight_length]

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.total_length, self.mm_width
