from __future__ import annotations

import gdsfactory as gf

from .base import DeviceGeometry

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


def _reset_gf():
    # Activate PDK without clearing the cache -- clearing destroys previously
    # returned Component objects, making them unusable after the next build.
    gf.gpdk.PDK.activate()


def _xs(width: float):
    return gf.cross_section.cross_section(width=width, offset=0, layer=_WG_LAYER)


def _straight(length: float, width: float) -> gf.Component:
    """Straight waveguide with explicit layer and width (no width= kwarg clash)."""
    xs = _xs(width)
    return gf.components.straight(length=length, cross_section=xs)


# ---------------------------------------------------------------------------
# 1×2 MMI
# ---------------------------------------------------------------------------

class MMI1x2(DeviceGeometry):
    """1×2 multimode interferometer splitter.

    Parameters mirror the legacy GeoParamsMMI fields.  `waveguide_offset` is
    the center-to-center distance between the two output waveguides.
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

    def _draw_polygons(self, c: gf.Component, wo: float) -> None:
        """Add 1×2 MMI polygon set to c.  wo = half waveguide_offset."""
        L, W = self.mm_length, self.mm_width
        tl, tw = self.taper_length, self.taper_width
        ww = self.waveguide_width

        c.add_polygon([(-L/2, -W/2), (L/2, -W/2), (L/2, W/2), (-L/2, W/2)], layer=_WG_LAYER)
        c.add_polygon([(-L/2-tl, -ww/2), (-L/2, -tw/2), (-L/2, tw/2), (-L/2-tl, ww/2)], layer=_WG_LAYER)
        for sign in (+1, -1):
            c.add_polygon([
                (L/2,    -tw/2 + sign*wo), (L/2+tl, -ww/2 + sign*wo),
                (L/2+tl,  ww/2 + sign*wo), (L/2,     tw/2 + sign*wo),
            ], layer=_WG_LAYER)

    def to_gds_eme(self, straight_length: float, oxide_width: float,
                   x_start: float = 0.0, **kwargs) -> gf.Component:
        _reset_gf()
        wo = self.waveguide_offset / 2
        ml = self.total_length
        c = gf.Component(f"mmi_1x2_eme_{id(self)}")
        self._draw_polygons(c, wo)
        c.move((x_start + ml/2 + straight_length, 0))
        in_wg = c.add_ref(_straight(straight_length, self.waveguide_width))
        in_wg.move((x_start, 0))
        out_wg = c.add_ref(_straight(straight_length, self.waveguide_width))
        out_wg.move((x_start + ml + straight_length, wo))
        c.add_polygon([
            (x_start, -oxide_width/2), (x_start + ml + 2*straight_length, -oxide_width/2),
            (x_start + ml + 2*straight_length, oxide_width/2), (x_start, oxide_width/2),
        ], layer=_OX_LAYER)
        return c

    def to_gds_fdtd(self, domain_sx: float, dpml: float, waveguide_length: float,
                    side: str = "input", **kwargs) -> gf.Component:
        _reset_gf()
        ml = self.total_length
        wo = self.waveguide_offset / 2
        name = "mmi_1x2_fdtd_in" if side == "input" else "mmi_1x2_fdtd_out"
        c = gf.Component(name)
        self._draw_polygons(c, wo)
        if side == "input":
            c.move((ml/2 - domain_sx/2 + dpml + waveguide_length, 0))
        else:
            c.move((-ml/2 + domain_sx/2 - dpml - waveguide_length, 0))
        return c

    @property
    def bbox(self):
        return self.total_length, self.mm_width


# ---------------------------------------------------------------------------
# 2×2 MMI
# ---------------------------------------------------------------------------

class MMI2x2(DeviceGeometry):
    """2×2 multimode interferometer splitter.

    `waveguide_offset` is the full center-to-center distance between the two
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

    def _draw_polygons(self, c: gf.Component, wo: float, name: str = "mmi_2x2") -> None:
        """Add 2×2 MMI polygon set to c.  wo = half waveguide_offset."""
        L, W = self.mm_length, self.mm_width
        tl, tw = self.taper_length, self.taper_width
        ww = self.waveguide_width

        c.add_polygon([(-L/2, -W/2), (L/2, -W/2), (L/2, W/2), (-L/2, W/2)], layer=_WG_LAYER)
        for sign in (+1, -1):
            c.add_polygon([
                (-L/2-tl, -ww/2 + sign*wo), (-L/2, -tw/2 + sign*wo),
                (-L/2,     tw/2 + sign*wo), (-L/2-tl,  ww/2 + sign*wo),
            ], layer=_WG_LAYER)
            c.add_polygon([
                (L/2,    -tw/2 + sign*wo), (L/2+tl, -ww/2 + sign*wo),
                (L/2+tl,  ww/2 + sign*wo), (L/2,     tw/2 + sign*wo),
            ], layer=_WG_LAYER)

    def to_gds_eme(self, straight_length: float, oxide_width: float,
                   x_start: float = 0.0, side: str = "top", **kwargs) -> gf.Component:
        """Build GDS for EME simulation.

        Args:
            side: "top" (through port: input top → output top) or
                  "bot" (cross port: input top → output bottom).
        """
        _reset_gf()
        wo = self.waveguide_offset / 2
        ml = self.total_length
        out_yo = wo if side == "top" else -wo
        name = f"mmi_2x2_eme_{side}_{id(self)}"

        c = gf.Component(name)
        self._draw_polygons(c, wo)
        c.move((x_start + ml/2 + straight_length, 0))

        in_wg = c.add_ref(_straight(straight_length, self.waveguide_width))
        in_wg.move((x_start, wo))
        out_wg = c.add_ref(_straight(straight_length, self.waveguide_width))
        out_wg.move((x_start + ml + straight_length, out_yo))

        c.add_polygon([
            (x_start, -oxide_width/2), (x_start + ml + 2*straight_length, -oxide_width/2),
            (x_start + ml + 2*straight_length, oxide_width/2), (x_start, oxide_width/2),
        ], layer=_OX_LAYER)
        return c

    def to_gds_fdtd(self, domain_sx: float, dpml: float, waveguide_length: float,
                    side: str = "input", name: str = "mmi_2x2_fdtd", **kwargs) -> gf.Component:
        """Build GDS positioned for a MEEP simulation domain."""
        _reset_gf()
        ml = self.total_length
        wo = self.waveguide_offset / 2
        c = gf.Component(name)
        self._draw_polygons(c, wo)
        if side == "input":
            c.move((ml/2 - domain_sx/2 + dpml + waveguide_length, 0))
        else:
            c.move((-ml/2 + domain_sx/2 - dpml - waveguide_length, 0))
        return c

    @property
    def bbox(self):
        return self.total_length, self.mm_width


# ---------------------------------------------------------------------------
# 1×4 MMI
# ---------------------------------------------------------------------------

class MMI1x4(DeviceGeometry):
    """1×4 multimode interferometer splitter.

    `waveguide_offset` is the center-to-center spacing between adjacent output
    waveguides; outputs are at ±wo/2 and ±3wo/2.
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

    def _draw_polygons(self, c: gf.Component) -> None:
        L, W = self.mm_length, self.mm_width
        tl, tw = self.taper_length, self.taper_width
        ww = self.waveguide_width
        wo = self.waveguide_offset

        c.add_polygon([(-L/2, -W/2), (L/2, -W/2), (L/2, W/2), (-L/2, W/2)], layer=_WG_LAYER)
        c.add_polygon([(-L/2-tl, -ww/2), (-L/2, -tw/2), (-L/2, tw/2), (-L/2-tl, ww/2)], layer=_WG_LAYER)
        for yo in [3*wo/2, wo/2, -wo/2, -3*wo/2]:
            c.add_polygon([
                (L/2,    -tw/2 + yo), (L/2+tl, -ww/2 + yo),
                (L/2+tl,  ww/2 + yo), (L/2,     tw/2 + yo),
            ], layer=_WG_LAYER)

    def to_gds_eme(self, straight_length: float, oxide_width: float,
                   x_start: float = 0.0, output_idx: int = 0, **kwargs) -> gf.Component:
        _reset_gf()
        ml = self.total_length
        wo = self.waveguide_offset
        output_offsets = [3*wo/2, wo/2, -wo/2, -3*wo/2]
        out_yo = output_offsets[output_idx]
        name = f"mmi_1x4_eme_{output_idx}"

        c = gf.Component(name)
        self._draw_polygons(c)
        c.move((x_start + ml/2 + straight_length, 0))

        in_wg = c.add_ref(_straight(straight_length, self.waveguide_width))
        in_wg.move((x_start, 0))
        out_wg = c.add_ref(_straight(straight_length, self.waveguide_width))
        out_wg.move((x_start + ml + straight_length, out_yo))

        c.add_polygon([
            (x_start, -oxide_width/2), (x_start + ml + 2*straight_length, -oxide_width/2),
            (x_start + ml + 2*straight_length, oxide_width/2), (x_start, oxide_width/2),
        ], layer=_OX_LAYER)
        return c

    def to_gds_fdtd(self, domain_sx: float, dpml: float, waveguide_length: float,
                    side: str = "input", **kwargs) -> gf.Component:
        _reset_gf()
        ml = self.total_length
        name = "mmi_1x4_fdtd_in" if side == "input" else "mmi_1x4_fdtd_out"
        c = gf.Component(name)
        self._draw_polygons(c)
        if side == "input":
            c.move((ml/2 - domain_sx/2 + dpml + waveguide_length, 0))
        else:
            c.move((-ml/2 + domain_sx/2 - dpml - waveguide_length, 0))
        return c

    @property
    def bbox(self):
        return self.total_length, self.mm_width
