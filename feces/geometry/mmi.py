"""Unified N-input / M-output multimode interferometer geometry."""
from __future__ import annotations

import itertools
from typing import Dict, List, Tuple

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


class MMI(DeviceGeometry):
    """N-input / M-output multimode interferometer.

    waveguide_offset is the center-to-center spacing between adjacent ports
    on either side; ports are placed symmetrically about y=0. n_inputs=1,
    n_outputs=2 gives a 1x2 splitter; n_inputs=2, n_outputs=2 a 2x2;
    n_inputs=1, n_outputs=4 a 1x4; etc.
    """

    def __init__(
        self,
        n_inputs: int,
        n_outputs: int,
        mm_length: float,
        mm_width: float,
        waveguide_width: float,
        taper_length: float,
        taper_width: float,
        waveguide_offset: float,
        waveguide_height: float = 0.16,
    ):
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs
        self.mm_length = mm_length
        self.mm_width = mm_width
        self.waveguide_width = waveguide_width
        self.taper_length = taper_length
        self.taper_width = taper_width
        self.waveguide_offset = waveguide_offset
        self.waveguide_height = waveguide_height

    @staticmethod
    def _port_ys(n: int, spacing: float) -> List[float]:
        """y-centers for n ports, spaced `spacing` apart, symmetric about y=0."""
        return [spacing * (i - (n - 1) / 2) for i in range(n)]

    @property
    def input_ys(self) -> List[float]:
        return self._port_ys(self.n_inputs, self.waveguide_offset)

    @property
    def output_ys(self) -> List[float]:
        return self._port_ys(self.n_outputs, self.waveguide_offset)

    @property
    def total_length(self) -> float:
        return self.mm_length + 2 * self.taper_length

    def polygons(
        self, stub_length: float = 0.0
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """MMI polygons in device frame (center at x=0)."""
        L, W = self.mm_length, self.mm_width
        tl, tw = self.taper_length, self.taper_width
        ww = self.waveguide_width
        ml = self.total_length

        wg = [[(-L/2, -W/2), (L/2, -W/2), (L/2, W/2), (-L/2, W/2)]]
        for yo in self.input_ys:
            wg.append([(-L/2-tl, -ww/2+yo), (-L/2, -tw/2+yo),
                       (-L/2, tw/2+yo), (-L/2-tl, ww/2+yo)])
        for yo in self.output_ys:
            wg.append([(L/2, -tw/2+yo), (L/2+tl, -ww/2+yo),
                       (L/2+tl, ww/2+yo), (L/2, tw/2+yo)])

        if stub_length > 0:
            sl = stub_length
            for yo in self.input_ys:
                wg.append([(-ml/2-sl, -ww/2+yo), (-ml/2, -ww/2+yo),
                           (-ml/2, ww/2+yo), (-ml/2-sl, ww/2+yo)])
            for yo in self.output_ys:
                wg.append([(ml/2, -ww/2+yo), (ml/2+sl, -ww/2+yo),
                           (ml/2+sl, ww/2+yo), (ml/2, ww/2+yo)])
        return {_WG_LAYER: wg}

    def to_gds_eme(
        self, straight_length: float, oxide_width: float,
        x_start: float = 0.0, **kwargs
    ) -> gf.Component:
        _reset_gf()
        ml = self.total_length
        x_off = x_start + ml / 2 + straight_length
        c = _component(f"mmi_{self.n_inputs}x{self.n_outputs}_eme_{id(self)}")
        for verts in self.polygons(stub_length=straight_length)[_WG_LAYER]:
            c.add_polygon([(x + x_off, y) for x, y in verts], layer=_WG_LAYER)
        _add_oxide(c, x_start, x_start + ml + 2 * straight_length, oxide_width)
        return c

    # ------------------------------------------------------------------
    # Arbitrary-section interface (HybridSolver)
    # ------------------------------------------------------------------
    #
    # Device frame: x=0 at the MMI body centre (same as polygons()).
    #   x in [-half_body-tl, -half_body] : input tapers (width ww -> tw)
    #   x in [-half_body, +half_body]    : MMI body (uniform, width W)
    #   x in [+half_body, +half_body+tl] : output tapers (width tw -> ww)
    # Anything beyond +-(half_body+tl) is the (virtual, arbitrary-length)
    # access waveguide region, uniform width ww.
    #
    # Cuts / stubs may land anywhere except strictly inside a taper (its
    # cross-section is not constant there). input_idx / output_idx select
    # which port to use once a range reaches into the input / output taper
    # or access-waveguide region; irrelevant (ignored) otherwise. Both
    # default to 0, which is exactly right for a single-port side (e.g. the
    # lone input of a 1xN device).

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
                           input_idx: int = 0, output_idx: int = 0):
        """Raw vertex lists (device frame) covering [x0-stub_left, x1+stub_right]."""
        (in0, in1), (out0, out1) = self._taper_bounds()
        ww, tw, W = self.waveguide_width, self.taper_width, self.mm_width
        y_in = self.input_ys[input_idx]
        y_out = self.output_ys[output_idx]

        lo, hi = x0 - stub_left, x1 + stub_right
        self._check_cut(lo)
        self._check_cut(hi)

        polys = []

        def _rect(xa, xb, half_w, yc=0.0):
            if xb - xa > 1e-9:
                polys.append([(xa, -half_w + yc), (xb, -half_w + yc),
                              (xb, half_w + yc), (xa, half_w + yc)])

        # left access waveguide (uniform, width ww, offset to the input port)
        if lo < in0:
            _rect(lo, min(hi, in0), ww / 2, yc=y_in)
        # input taper (all-or-nothing — never cut mid-taper)
        if lo <= in0 and hi >= in1:
            polys.append([(in0, -ww/2 + y_in), (in1, -tw/2 + y_in),
                          (in1, tw/2 + y_in), (in0, ww/2 + y_in)])
        # MMI body (uniform, arbitrary sub-range)
        _rect(max(lo, in1), min(hi, out0), W / 2)
        # output taper (all-or-nothing, offset to the selected branch)
        if lo <= out0 and hi >= out1:
            polys.append([(out0, -tw/2 + y_out), (out1, -ww/2 + y_out),
                          (out1, ww/2 + y_out), (out0, tw/2 + y_out)])
        # right access waveguide (uniform, width ww, offset to the output port)
        if hi > out1:
            _rect(max(lo, out1), hi, ww / 2, yc=y_out)

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
                     input_idx: int = 0, output_idx: int = 0) -> gf.Component:
        _reset_gf()
        lo, hi = x0 - stub_left, x1 + stub_right
        x_off = -lo
        c = _component(f"mmi_section_{id(self)}")
        for verts in self._section_polygons(x0, x1, stub_left, stub_right, input_idx, output_idx):
            c.add_polygon([(x + x_off, y) for x, y in verts], layer=_WG_LAYER)
        _add_oxide(c, 0.0, hi - lo, oxide_width)
        return c

    def width_at(self, x: float, input_idx: int = 0, output_idx: int = 0) -> float:
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

    def center_at(self, x: float, input_idx: int = 0, output_idx: int = 0) -> float:
        half_body = self.mm_length / 2
        if x < -half_body:
            return self.input_ys[input_idx]
        if x > half_body:
            return self.output_ys[output_idx]
        return 0.0

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.taper_length, self.mm_length,
                self.taper_length, straight_length]

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.total_length, self.mm_width
