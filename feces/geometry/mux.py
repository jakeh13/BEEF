"""MMI MUX geometry: cascaded 1×3 + 3×3 MMI multiplexer.

Polygon logic migrated from engine/helpers/meow.py and engine/helpers/meep.py.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

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


def _straight_rect(x0: float, y_center: float, length: float, width: float) -> List[Tuple[float, float]]:
    return [(x0, y_center - width/2), (x0 + length, y_center - width/2),
            (x0 + length, y_center + width/2), (x0, y_center + width/2)]


def _add_cascade_polygons(geo: MuxGeoParams, x_start: float = 0.0) -> List[List[Tuple[float, float]]]:
    """MMI1_body | 3×out_tapers | 3×conn_wgs | 3×in_tapers2 | MMI2_body (no input/output taper)."""
    L1, L2 = geo.mm_1_length, geo.mm_2_length
    Lt, Lc = geo.taper_length, geo.between_waveguide_length
    mw_val, wo = geo.mm_width, geo.waveguide_offset
    tbw, midw = geo.top_bot_taper_width, geo.middle_taper_width
    tww, ww = geo.top_waveguide_width, geo.waveguide_width

    polys: List[List[Tuple[float, float]]] = []
    x = x_start
    polys.append([(x, -mw_val/2), (x+L1, -mw_val/2), (x+L1, mw_val/2), (x, mw_val/2)])
    x += L1
    for (yo, mmi_w, wg_w) in [(wo, tbw, tww), (0, midw, ww), (-wo, tbw, ww)]:
        polys.append([(x, -mmi_w/2+yo), (x+Lt, -wg_w/2+yo), (x+Lt, wg_w/2+yo), (x, mmi_w/2+yo)])
    x += Lt
    for (yo, w) in [(wo, tww), (0, ww), (-wo, ww)]:
        polys.append([(x, -w/2+yo), (x+Lc, -w/2+yo), (x+Lc, w/2+yo), (x, w/2+yo)])
    x += Lc
    for (yo, wg_w, mmi_w) in [(wo, tww, tbw), (0, ww, midw), (-wo, ww, tbw)]:
        polys.append([(x, -wg_w/2+yo), (x+Lt, -mmi_w/2+yo), (x+Lt, mmi_w/2+yo), (x, wg_w/2+yo)])
    x += Lt
    polys.append([(x, -mw_val/2), (x+L2, -mw_val/2), (x+L2, mw_val/2), (x, mw_val/2)])
    return polys


def _add_mux_polygons(geo: MuxGeoParams, x_start: float = 0.0) -> List[List[Tuple[float, float]]]:
    """Full polygon sequence: input_taper | MMI1 | out_tapers | conn_wgs | in_tapers | MMI2 | out_tapers."""
    L1, L2 = geo.mm_1_length, geo.mm_2_length
    Lt, Lc = geo.taper_length, geo.between_waveguide_length
    mw_val, wo = geo.mm_width, geo.waveguide_offset
    iw, itw = geo.mm_waveguide_width, geo.input_taper_width
    tbw, midw = geo.top_bot_taper_width, geo.middle_taper_width
    tww, ww = geo.top_waveguide_width, geo.waveguide_width

    polys: List[List[Tuple[float, float]]] = []
    x = x_start
    polys.append([(x, -iw/2), (x+Lt, -itw/2), (x+Lt, itw/2), (x, iw/2)])
    x += Lt
    polys.append([(x, -mw_val/2), (x+L1, -mw_val/2), (x+L1, mw_val/2), (x, mw_val/2)])
    x += L1
    for (yo, mmi_w, wg_w) in [(wo, tbw, tww), (0, midw, ww), (-wo, tbw, ww)]:
        polys.append([(x, -mmi_w/2+yo), (x+Lt, -wg_w/2+yo), (x+Lt, wg_w/2+yo), (x, mmi_w/2+yo)])
    x += Lt
    for (yo, w) in [(wo, tww), (0, ww), (-wo, ww)]:
        polys.append([(x, -w/2+yo), (x+Lc, -w/2+yo), (x+Lc, w/2+yo), (x, w/2+yo)])
    x += Lc
    for (yo, wg_w, mmi_w) in [(wo, tww, tbw), (0, ww, midw), (-wo, ww, tbw)]:
        polys.append([(x, -wg_w/2+yo), (x+Lt, -mmi_w/2+yo), (x+Lt, mmi_w/2+yo), (x, wg_w/2+yo)])
    x += Lt
    polys.append([(x, -mw_val/2), (x+L2, -mw_val/2), (x+L2, mw_val/2), (x, mw_val/2)])
    x += L2
    for (yo, mmi_w, wg_w) in [(wo, tbw, tww), (0, midw, ww), (-wo, tbw, ww)]:
        polys.append([(x, -mmi_w/2+yo), (x+Lt, -wg_w/2+yo), (x+Lt, wg_w/2+yo), (x, mmi_w/2+yo)])
    return polys


class MMIMux(DeviceGeometry):
    """Cascaded 1×3 + 3×3 MMI multiplexer geometry."""

    def __init__(self, geo: MuxGeoParams):
        self.geo = geo

    @property
    def total_length(self) -> float:
        g = self.geo
        return 4*g.taper_length + g.mm_1_length + g.between_waveguide_length + g.mm_2_length

    def to_gds_eme(self, straight_length: float, oxide_width: float,
                   x_start: float = 0.0, port_idx: int | None = None, **kwargs
                   ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """Build full MUX layout for EME.

        Args:
            port_idx: if None, include all 3 output arms; if 0/1/2 include only that arm.
        """
        geo = self.geo
        sl = straight_length
        ow = oxide_width
        device_len = self.total_length
        total_len = sl + device_len + sl

        wg_polys = _add_mux_polygons(geo, x_start=sl)
        wg_polys.append(_straight_rect(x_start, 0.0, sl, geo.mm_waveguide_width))
        if port_idx is None:
            for (yo, w) in [(geo.waveguide_offset, geo.top_waveguide_width),
                            (0, geo.waveguide_width), (-geo.waveguide_offset, geo.waveguide_width)]:
                wg_polys.append(_straight_rect(x_start + sl + device_len, yo, sl, w))
        else:
            port_y = [geo.waveguide_offset, 0.0, -geo.waveguide_offset][port_idx]
            port_w = [geo.top_waveguide_width, geo.waveguide_width, geo.waveguide_width][port_idx]
            wg_polys.append(_straight_rect(x_start + sl + device_len, port_y, sl, port_w))

        ox_poly = [
            (x_start, -ow/2), (x_start + total_len, -ow/2),
            (x_start + total_len, ow/2), (x_start, ow/2),
        ]
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    def to_gds_eme_cascade(self, oxide_width: float, mm_monitor_position: float, **kwargs
                            ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """Build cascade-only layout for hybrid EME (no access waveguides, aligned at monitor)."""
        geo = self.geo
        x_start = -mm_monitor_position
        total_gds_len = geo.mm_1_length + 2*geo.taper_length + geo.between_waveguide_length + geo.mm_2_length
        wg_polys = _add_cascade_polygons(geo, x_start=x_start)
        ox_poly = [
            (x_start, -oxide_width/2), (x_start + total_gds_len, -oxide_width/2),
            (x_start + total_gds_len, oxide_width/2), (x_start, oxide_width/2),
        ]
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    def to_gds_fdtd(self, domain_sx: float, dpml: float, waveguide_length: float,
                    side: str = "input", port_idx: int = 0, mm_sim_length: float = 1.0,
                    **kwargs) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """Build partial MUX layout for FDTD (input or output side)."""
        geo = self.geo
        Lt = geo.taper_length

        if side == "input":
            x_shift = -domain_sx/2 + dpml + waveguide_length
            iw, itw = geo.mm_waveguide_width, geo.input_taper_width
            wg_polys = [
                [(x_shift, -iw/2), (x_shift+Lt, -itw/2), (x_shift+Lt, itw/2), (x_shift, iw/2)],
                [(x_shift+Lt, -geo.mm_width/2), (x_shift+Lt+mm_sim_length, -geo.mm_width/2),
                 (x_shift+Lt+mm_sim_length, geo.mm_width/2), (x_shift+Lt, geo.mm_width/2)],
            ]
        else:
            mw_val = geo.mm_width
            port_params = [
                (geo.waveguide_offset,  geo.top_bot_taper_width, geo.top_waveguide_width),
                (0,                     geo.middle_taper_width,  geo.waveguide_width),
                (-geo.waveguide_offset, geo.top_bot_taper_width, geo.waveguide_width),
            ]
            yo, mmi_w, wg_w = port_params[port_idx]
            x_shift_out = (domain_sx/2 - dpml - waveguide_length) - Lt
            wg_polys = [
                [(x_shift_out-mm_sim_length, -mw_val/2), (x_shift_out, -mw_val/2),
                 (x_shift_out, mw_val/2), (x_shift_out-mm_sim_length, mw_val/2)],
                [(x_shift_out, -mmi_w/2+yo), (x_shift_out+Lt, -wg_w/2+yo),
                 (x_shift_out+Lt, wg_w/2+yo), (x_shift_out, mmi_w/2+yo)],
            ]

        return {_WG_LAYER: wg_polys}

    def eme_lengths(self, straight_length: float):
        g = self.geo
        return [
            straight_length,
            g.taper_length,
            g.mm_1_length,
            g.taper_length,
            g.between_waveguide_length,
            g.taper_length,
            g.mm_2_length,
            g.taper_length,
            straight_length,
        ]

    @property
    def bbox(self):
        return self.total_length, self.geo.mm_width
