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
    """Cascaded 1×3 + 3×3 MMI multiplexer geometry.

    Device frame (used by polygons() and the arbitrary-section interface):
    x=0 at the device centre, so the layout spans [-total_length/2,
    +total_length/2] with regions left→right:

        input taper | MMI 1 body | out tapers | conn wgs | in tapers |
        MMI 2 body | out tapers

    The single multi-mode input waveguide sits on-axis (y=0); the three
    inter-stage / output branches sit at y = +waveguide_offset ("top", index
    0), y = 0 ("mid", index 1) and y = -waveguide_offset ("bot", index 2) —
    the same port_idx convention as to_gds_fdtd() and the legacy master.
    """

    def __init__(self, geo: MuxGeoParams):
        self.geo = geo

    @property
    def total_length(self) -> float:
        g = self.geo
        return 4*g.taper_length + g.mm_1_length + g.between_waveguide_length + g.mm_2_length

    def _branches(self) -> List[Tuple[float, float, float]]:
        """(y_center, taper_width_at_MMI_face, waveguide_width) per branch,
        index 0=top, 1=mid, 2=bot."""
        g = self.geo
        return [
            (+g.waveguide_offset, g.top_bot_taper_width, g.top_waveguide_width),
            (0.0,                 g.middle_taper_width,  g.waveguide_width),
            (-g.waveguide_offset, g.top_bot_taper_width, g.waveguide_width),
        ]

    def polygons(
        self, stub_length: float = 0.0
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """Full MUX polygons in device frame (device centre at x=0)."""
        g = self.geo
        half = self.total_length / 2
        wg = _add_mux_polygons(g, x_start=-half)
        if stub_length > 0:
            sl = stub_length
            wg.append(_straight_rect(-half - sl, 0.0, sl, g.mm_waveguide_width))
            for (yo, _, w) in self._branches():
                wg.append(_straight_rect(half, yo, sl, w))
        return {_WG_LAYER: wg}

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

    # ------------------------------------------------------------------
    # Arbitrary-section interface (HybridSolver)
    # ------------------------------------------------------------------
    #
    # Device frame: x=0 at the device centre (see class docstring). Cuts /
    # stubs may land anywhere except strictly inside one of the four taper
    # bands (their cross-section varies with x); the uniform regions — the
    # access waveguides, both MMI bodies and the connecting waveguides —
    # accept arbitrary sub-ranges.

    def _region_bounds(self):
        """Region (x_start, x_end) edges left→right in the device frame:
        (t_in, body1, t_out1, conn, t_in2, body2, t_out2)."""
        g = self.geo
        x = -self.total_length / 2
        bounds = []
        for L in (g.taper_length, g.mm_1_length, g.taper_length,
                  g.between_waveguide_length, g.taper_length, g.mm_2_length,
                  g.taper_length):
            bounds.append((x, x + L))
            x += L
        return tuple(bounds)

    def _check_cut(self, x: float) -> None:
        t_in, _, t_out1, _, t_in2, _, t_out2 = self._region_bounds()
        for (a, b) in (t_in, t_out1, t_in2, t_out2):
            if a < x < b:
                raise ValueError(
                    f"Cut/stub at x={x:.4f} falls strictly inside a taper "
                    "(cross-section varies there) — cuts must land on a taper "
                    "endpoint or inside a uniform region (access waveguide, "
                    "MMI body, or connecting waveguides)."
                )

    def _section_polygons(self, x0: float, x1: float,
                           stub_left: float = 0.0, stub_right: float = 0.0,
                           port_idx: int | None = None):
        """Raw vertex lists (device frame) covering [x0-stub_left, x1+stub_right].

        Every branch is always drawn in full — the inter-stage structure
        (MMI 1 output tapers, connecting waveguides, MMI 2 input tapers)
        AND the output side (MMI 2 output tapers + output access
        waveguides) — so the simulated geometry never differs from the real
        device. port_idx is accepted (HybridSolver forwards one
        geometry_kwargs dict to every geometry hook) but never changes what
        is drawn: per-port isolation for the output FDTD sections happens
        in the eigenmode *launch window* instead — see source_aperture_at().
        (The legacy reference implementation, FECES_Legacy/engine/helpers/
        meep.py's create_temp_gds_for_mmi_mux_output, instead dropped the
        two unselected output branches from the drawing; shrinking the
        source aperture achieves the same clean single-branch mode solve
        without altering the structure.)
        """
        t_in, body1, t_out1, conn, t_in2, body2, t_out2 = self._region_bounds()
        g = self.geo
        branches = self._branches()

        lo, hi = x0 - stub_left, x1 + stub_right
        self._check_cut(lo)
        self._check_cut(hi)

        polys: List[List[Tuple[float, float]]] = []

        def _rect(xa, xb, yc, w):
            if xb - xa > 1e-9:
                polys.append([(xa, yc - w/2), (xb, yc - w/2),
                              (xb, yc + w/2), (xa, yc + w/2)])

        def _taper(xa, xb, yc, w_a, w_b):
            polys.append([(xa, yc - w_a/2), (xb, yc - w_b/2),
                          (xb, yc + w_b/2), (xa, yc + w_a/2)])

        # input MM access waveguide (uniform, on-axis)
        if lo < t_in[0]:
            _rect(lo, min(hi, t_in[0]), 0.0, g.mm_waveguide_width)
        # input taper (all-or-nothing — never cut mid-taper)
        if lo <= t_in[0] and hi >= t_in[1]:
            _taper(t_in[0], t_in[1], 0.0, g.mm_waveguide_width, g.input_taper_width)
        # MMI 1 body (uniform, arbitrary sub-range)
        _rect(max(lo, body1[0]), min(hi, body1[1]), 0.0, g.mm_width)
        # MMI 1 output tapers (all-or-nothing, all three branches)
        if lo <= t_out1[0] and hi >= t_out1[1]:
            for (yo, mmi_w, wg_w) in branches:
                _taper(t_out1[0], t_out1[1], yo, mmi_w, wg_w)
        # connecting waveguides (uniform, all three branches)
        for (yo, _, wg_w) in branches:
            _rect(max(lo, conn[0]), min(hi, conn[1]), yo, wg_w)
        # MMI 2 input tapers (all-or-nothing, all three branches)
        if lo <= t_in2[0] and hi >= t_in2[1]:
            for (yo, mmi_w, wg_w) in branches:
                _taper(t_in2[0], t_in2[1], yo, wg_w, mmi_w)
        # MMI 2 body (uniform, arbitrary sub-range)
        _rect(max(lo, body2[0]), min(hi, body2[1]), 0.0, g.mm_width)
        # MMI 2 output tapers (all-or-nothing, all three branches)
        if lo <= t_out2[0] and hi >= t_out2[1]:
            for (yo, mmi_w, wg_w) in branches:
                _taper(t_out2[0], t_out2[1], yo, mmi_w, wg_w)
        # output access waveguides (uniform, all three branches)
        if hi > t_out2[1]:
            for (yo, _, wg_w) in branches:
                _rect(max(lo, t_out2[1]), hi, yo, wg_w)

        return polys

    def section_prisms(self, x0: float, x1: float, si_medium,
                        stub_left: float = 0.0, stub_right: float = 0.0,
                        port_idx: int | None = None,
                        wg_height: float | None = None):
        from .base import polygons_to_prisms
        return polygons_to_prisms(
            self._section_polygons(x0, x1, stub_left, stub_right, port_idx),
            si_medium, wg_height=wg_height,
        )

    def section_gds(self, x0: float, x1: float, oxide_width: float,
                     stub_left: float = 0.0, stub_right: float = 0.0,
                     port_idx: int | None = None
                     ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        lo, hi = x0 - stub_left, x1 + stub_right
        x_off = -lo
        wg_polys = [
            [(x + x_off, y) for x, y in verts]
            for verts in self._section_polygons(x0, x1, stub_left, stub_right, port_idx)
        ]
        ox_poly = [
            (0.0, -oxide_width/2), (hi - lo, -oxide_width/2),
            (hi - lo, oxide_width/2), (0.0, oxide_width/2),
        ]
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    def width_at(self, x: float, port_idx: int | None = None) -> float:
        """Local cross-section width at device-frame x (µm). port_idx selects
        which branch (0=top, 1=mid, 2=bot) to report inside the multi-branch
        inter-stage/output regions; None defaults to the on-axis mid branch."""
        t_in, body1, t_out1, conn, t_in2, body2, t_out2 = self._region_bounds()
        g = self.geo
        _, mmi_w, wg_w = self._branches()[1 if port_idx is None else port_idx]

        if x <= t_in[0]:
            return g.mm_waveguide_width
        if x <= t_in[1]:
            frac = (x - t_in[0]) / g.taper_length
            return g.mm_waveguide_width + frac * (g.input_taper_width - g.mm_waveguide_width)
        if x <= body1[1]:
            return g.mm_width
        if x <= t_out1[1]:
            frac = (x - t_out1[0]) / g.taper_length
            return mmi_w + frac * (wg_w - mmi_w)
        if x <= conn[1]:
            return wg_w
        if x <= t_in2[1]:
            frac = (x - t_in2[0]) / g.taper_length
            return wg_w + frac * (mmi_w - wg_w)
        if x <= body2[1]:
            return g.mm_width
        if x <= t_out2[1]:
            frac = (x - t_out2[0]) / g.taper_length
            return mmi_w + frac * (wg_w - mmi_w)
        return wg_w

    def center_at(self, x: float, port_idx: int | None = None) -> float:
        """Local waveguide y-centre at device-frame x (µm) — see width_at()
        for the port_idx convention."""
        _, body1, _, _, _, body2, _ = self._region_bounds()
        yo = self._branches()[1 if port_idx is None else port_idx][0]
        if x <= body1[1]:
            return 0.0   # input wg / input taper / MMI 1 body — all on-axis
        if x < body2[0]:
            return yo    # inter-stage branches (tapers + connecting wgs)
        if x <= body2[1]:
            return 0.0   # MMI 2 body
        return yo        # output tapers / output access waveguides

    def source_aperture_at(self, x: float, port_idx: int | None = None) -> float | None:
        """Widest eigenmode source / norm-monitor window (full y-span, µm),
        centered on center_at(x, port_idx), that stays clear of every OTHER
        branch's core at this cross-section. With three branches only
        waveguide_offset (~1 µm) apart, a full-domain launch aperture
        centered on one port would also reach its neighbors, corrupting the
        "single-mode" MPB solve into seeing the near-degenerate supermodes
        of the branch trio (and half_domain_aperture can't isolate the mid
        port, which sits exactly on the y=0 axis) — this window keeps the
        drawn geometry complete and shrinks the launch instead. Returns
        None where the cross-section holds no neighboring branch to avoid
        (the on-axis input waveguide, the MMI bodies) or when port_idx is
        None; FDTDSolver.solve_section() additionally clamps the window to
        its own non-PML domain."""
        if port_idx is None:
            return None
        _, body1, _, _, _, body2, _ = self._region_bounds()
        if not (body1[1] < x < body2[0] or x > body2[1]):
            return None  # single on-axis structure here — no neighbor to avoid
        branches = self._branches()
        y_k = branches[port_idx][0]
        half = min(
            abs(yo - y_k) - self.width_at(x, port_idx=j) / 2
            for j, (yo, _, _) in enumerate(branches)
            if j != port_idx
        )
        return 2.0 * half

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
