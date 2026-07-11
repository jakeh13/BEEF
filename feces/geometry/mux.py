"""MMI MUX geometry: cascaded 1×3 + 3×3 MMI multiplexer.

Polygon logic migrated from engine/helpers/meow.py and engine/helpers/meep.py.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import gdsfactory as gf

from .base import DeviceGeometry

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


@dataclass
class MuxGeoParams:
    """Geometry parameters for the cascaded MUX device."""
    mm_1_length: float = 10.0
    mm_2_length: float = 10.0
    taper_length: float = 1.0
    between_waveguide_length: float = 2.0
    mm_width: float = 3.0
    waveguide_offset: float = 1.0
    mm_waveguide_width: float = 1.0
    input_taper_width: float = 1.5
    top_bot_taper_width: float = 0.8
    middle_taper_width: float = 0.8
    top_waveguide_width: float = 0.5
    waveguide_width: float = 0.35


def _reset_gf():
    gf.clear_cache()
    gf.gpdk.PDK.activate()


def _xs(width: float):
    return gf.cross_section.cross_section(width=width, offset=0, layer=_WG_LAYER)


def _add_cascade_polygons(c: gf.Component, geo: MuxGeoParams, x_start: float = 0.0) -> None:
    """MMI1_body | 3×out_tapers | 3×conn_wgs | 3×in_tapers2 | MMI2_body (no input/output taper)."""
    L1, L2 = geo.mm_1_length, geo.mm_2_length
    Lt, Lc = geo.taper_length, geo.between_waveguide_length
    mw_val, wo = geo.mm_width, geo.waveguide_offset
    tbw, midw = geo.top_bot_taper_width, geo.middle_taper_width
    tww, ww = geo.top_waveguide_width, geo.waveguide_width

    x = x_start
    c.add_polygon([(x, -mw_val/2), (x+L1, -mw_val/2), (x+L1, mw_val/2), (x, mw_val/2)], layer=_WG_LAYER)
    x += L1
    for (yo, mmi_w, wg_w) in [(wo, tbw, tww), (0, midw, ww), (-wo, tbw, ww)]:
        c.add_polygon([(x, -mmi_w/2+yo), (x+Lt, -wg_w/2+yo), (x+Lt, wg_w/2+yo), (x, mmi_w/2+yo)], layer=_WG_LAYER)
    x += Lt
    for (yo, w) in [(wo, tww), (0, ww), (-wo, ww)]:
        c.add_polygon([(x, -w/2+yo), (x+Lc, -w/2+yo), (x+Lc, w/2+yo), (x, w/2+yo)], layer=_WG_LAYER)
    x += Lc
    for (yo, wg_w, mmi_w) in [(wo, tww, tbw), (0, ww, midw), (-wo, ww, tbw)]:
        c.add_polygon([(x, -wg_w/2+yo), (x+Lt, -mmi_w/2+yo), (x+Lt, mmi_w/2+yo), (x, wg_w/2+yo)], layer=_WG_LAYER)
    x += Lt
    c.add_polygon([(x, -mw_val/2), (x+L2, -mw_val/2), (x+L2, mw_val/2), (x, mw_val/2)], layer=_WG_LAYER)


def _add_mux_polygons(c: gf.Component, geo: MuxGeoParams, x_start: float = 0.0) -> None:
    """Full polygon sequence: input_taper | MMI1 | out_tapers | conn_wgs | in_tapers | MMI2 | out_tapers."""
    L1, L2 = geo.mm_1_length, geo.mm_2_length
    Lt, Lc = geo.taper_length, geo.between_waveguide_length
    mw_val, wo = geo.mm_width, geo.waveguide_offset
    iw, itw = geo.mm_waveguide_width, geo.input_taper_width
    tbw, midw = geo.top_bot_taper_width, geo.middle_taper_width
    tww, ww = geo.top_waveguide_width, geo.waveguide_width

    x = x_start
    c.add_polygon([(x, -iw/2), (x+Lt, -itw/2), (x+Lt, itw/2), (x, iw/2)], layer=_WG_LAYER)
    x += Lt
    c.add_polygon([(x, -mw_val/2), (x+L1, -mw_val/2), (x+L1, mw_val/2), (x, mw_val/2)], layer=_WG_LAYER)
    x += L1
    for (yo, mmi_w, wg_w) in [(wo, tbw, tww), (0, midw, ww), (-wo, tbw, ww)]:
        c.add_polygon([(x, -mmi_w/2+yo), (x+Lt, -wg_w/2+yo), (x+Lt, wg_w/2+yo), (x, mmi_w/2+yo)], layer=_WG_LAYER)
    x += Lt
    for (yo, w) in [(wo, tww), (0, ww), (-wo, ww)]:
        c.add_polygon([(x, -w/2+yo), (x+Lc, -w/2+yo), (x+Lc, w/2+yo), (x, w/2+yo)], layer=_WG_LAYER)
    x += Lc
    for (yo, wg_w, mmi_w) in [(wo, tww, tbw), (0, ww, midw), (-wo, ww, tbw)]:
        c.add_polygon([(x, -wg_w/2+yo), (x+Lt, -mmi_w/2+yo), (x+Lt, mmi_w/2+yo), (x, wg_w/2+yo)], layer=_WG_LAYER)
    x += Lt
    c.add_polygon([(x, -mw_val/2), (x+L2, -mw_val/2), (x+L2, mw_val/2), (x, mw_val/2)], layer=_WG_LAYER)
    x += L2
    for (yo, mmi_w, wg_w) in [(wo, tbw, tww), (0, midw, ww), (-wo, tbw, ww)]:
        c.add_polygon([(x, -mmi_w/2+yo), (x+Lt, -wg_w/2+yo), (x+Lt, wg_w/2+yo), (x, mmi_w/2+yo)], layer=_WG_LAYER)


class MMIMux(DeviceGeometry):
    """Cascaded 1×3 + 3×3 MMI multiplexer geometry."""

    def __init__(self, geo: MuxGeoParams):
        self.geo = geo

    @property
    def total_length(self) -> float:
        g = self.geo
        return 4*g.taper_length + g.mm_1_length + g.between_waveguide_length + g.mm_2_length

    def to_gds_eme(self, straight_length: float, oxide_width: float,
                   x_start: float = 0.0, port_idx: int | None = None, **kwargs) -> gf.Component:
        """Build full MUX GDS for EME.

        Args:
            port_idx: if None, include all 3 output arms; if 0/1/2 include only that arm.
        """
        _reset_gf()
        geo = self.geo
        sl = straight_length
        ow = oxide_width
        device_len = self.total_length
        total_len = sl + device_len + sl

        if port_idx is None:
            c = gf.Component("mmi_mux_full_eme")
            _add_mux_polygons(c, geo, x_start=sl)
            inp_xs = _xs(geo.mm_waveguide_width)
            inp = c.add_ref(gf.components.straight(length=sl, width=geo.mm_waveguide_width, cross_section=inp_xs))
            inp.move((x_start, 0))
            for (yo, w) in [(geo.waveguide_offset, geo.top_waveguide_width),
                            (0, geo.waveguide_width), (-geo.waveguide_offset, geo.waveguide_width)]:
                out = c.add_ref(gf.components.straight(length=sl, width=w, cross_section=_xs(w)))
                out.move((x_start + sl + device_len, yo))
        else:
            port_y = [geo.waveguide_offset, 0.0, -geo.waveguide_offset][port_idx]
            port_w = [geo.top_waveguide_width, geo.waveguide_width, geo.waveguide_width][port_idx]
            c = gf.Component(f"mmi_mux_port{port_idx}_eme")
            _add_mux_polygons(c, geo, x_start=sl)
            inp = c.add_ref(gf.components.straight(length=sl, width=geo.mm_waveguide_width, cross_section=_xs(geo.mm_waveguide_width)))
            inp.move((x_start, 0))
            out = c.add_ref(gf.components.straight(length=sl, width=port_w, cross_section=_xs(port_w)))
            out.move((x_start + sl + device_len, port_y))

        c.add_polygon([
            (x_start, -ow/2), (x_start + total_len, -ow/2),
            (x_start + total_len, ow/2), (x_start, ow/2),
        ], layer=_OX_LAYER)
        return c

    def to_gds_eme_cascade(self, oxide_width: float, mm_monitor_position: float, **kwargs) -> gf.Component:
        """Build cascade-only GDS for HYBRID EME (no access waveguides, aligned at monitor)."""
        _reset_gf()
        geo = self.geo
        x_start = -mm_monitor_position
        total_gds_len = geo.mm_1_length + 2*geo.taper_length + geo.between_waveguide_length + geo.mm_2_length
        c = gf.Component("mmi_mux_cascade_eme")
        _add_cascade_polygons(c, geo, x_start=x_start)
        c.add_polygon([
            (x_start, -oxide_width/2), (x_start + total_gds_len, -oxide_width/2),
            (x_start + total_gds_len, oxide_width/2), (x_start, oxide_width/2),
        ], layer=_OX_LAYER)
        return c

    def to_gds_fdtd(self, domain_sx: float, dpml: float, waveguide_length: float,
                    side: str = "input", port_idx: int = 0, mm_sim_length: float = 1.0,
                    **kwargs) -> gf.Component:
        """Build partial MUX GDS for FDTD (input or output side)."""
        _reset_gf()
        geo = self.geo
        Lt = geo.taper_length
        x_shift = (-domain_sx/2 + dpml + waveguide_length) if side == "input" else None

        if side == "input":
            c = gf.Component("mmi_mux_fdtd_in")
            iw, itw = geo.mm_waveguide_width, geo.input_taper_width
            c.add_polygon([(0, -iw/2), (Lt, -itw/2), (Lt, itw/2), (0, iw/2)], layer=_WG_LAYER)
            c.add_polygon([(Lt, -geo.mm_width/2), (Lt + mm_sim_length, -geo.mm_width/2),
                           (Lt + mm_sim_length, geo.mm_width/2), (Lt, geo.mm_width/2)], layer=_WG_LAYER)
            c.move((x_shift, 0))
        else:
            mw_val = geo.mm_width
            port_params = [
                (geo.waveguide_offset,  geo.top_bot_taper_width, geo.top_waveguide_width),
                (0,                     geo.middle_taper_width,  geo.waveguide_width),
                (-geo.waveguide_offset, geo.top_bot_taper_width, geo.waveguide_width),
            ]
            yo, mmi_w, wg_w = port_params[port_idx]
            c = gf.Component(f"mmi_mux_fdtd_out_{port_idx}")
            c.add_polygon([(-mm_sim_length, -mw_val/2), (0, -mw_val/2),
                           (0, mw_val/2), (-mm_sim_length, mw_val/2)], layer=_WG_LAYER)
            c.add_polygon([(0, -mmi_m/2+yo), (Lt, -wg_w/2+yo), (Lt, wg_w/2+yo), (0, mmi_m/2+yo)], layer=_WG_LAYER)
            x_shift_out = (domain_sx/2 - dpml - waveguide_length) - Lt
            c.move((x_shift_out, 0))
        return c

    @property
    def bbox(self):
        return self.total_length, self.geo.mm_width
