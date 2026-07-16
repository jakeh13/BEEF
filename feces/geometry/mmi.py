"""Unified N-input / M-output multimode interferometer geometry."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from ._compose import shift_polys
from .base import DeviceGeometry, _oxide_box, polygons_to_prisms
from .taper import Taper

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


@dataclass
class MMIPort:
    """One input or output port: y-position plus its own taper geometry.

    The taper always spans from `waveguide_width` (at the access-waveguide
    end) to `taper_width` (at the MMI body face) over `taper_length`."""
    y: float
    waveguide_width: float
    taper_length: float
    taper_width: float


class MMI(DeviceGeometry):
    """N-input / M-output multimode interferometer.

    Ports default to a single shared taper_length/taper_width/waveguide_width,
    spaced `waveguide_offset` apart and symmetric about y=0 — n_inputs=1,
    n_outputs=2 gives a 1x2 splitter; n_inputs=2, n_outputs=2 a 2x2;
    n_inputs=1, n_outputs=4 a 1x4; etc. Pass `inputs`/`outputs` (lists of
    MMIPort) instead to give individual ports their own taper size and
    arbitrary y-placement.
    """

    def __init__(
        self,
        n_inputs: int,
        n_outputs: int,
        mm_length: float,
        mm_width: float,
        waveguide_width: float = 0.35,
        taper_length: float = 1.0,
        taper_width: float = 0.8,
        waveguide_offset: float = 1.0,
        inputs: Optional[List[MMIPort]] = None,
        outputs: Optional[List[MMIPort]] = None,
        waveguide_height: float = 0.16,
        use_neighbor_aperture: bool = False,
    ):
        self.mm_length = mm_length
        self.mm_width = mm_width
        self.waveguide_height = waveguide_height
        self.use_neighbor_aperture = use_neighbor_aperture

        def _default_ports(n: int) -> List[MMIPort]:
            return [
                MMIPort(y=y, waveguide_width=waveguide_width,
                        taper_length=taper_length, taper_width=taper_width)
                for y in self._port_ys(n, waveguide_offset)
            ]

        self.inputs = inputs if inputs is not None else _default_ports(n_inputs)
        self.outputs = outputs if outputs is not None else _default_ports(n_outputs)
        if len(self.inputs) != n_inputs:
            raise ValueError(f"len(inputs)={len(self.inputs)} != n_inputs={n_inputs}")
        if len(self.outputs) != n_outputs:
            raise ValueError(f"len(outputs)={len(self.outputs)} != n_outputs={n_outputs}")
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs

    @staticmethod
    def _port_ys(n: int, spacing: float) -> List[float]:
        """y-centers for n ports, spaced `spacing` apart, symmetric about y=0."""
        return [spacing * (i - (n - 1) / 2) for i in range(n)]

    @property
    def input_ys(self) -> List[float]:
        return [p.y for p in self.inputs]

    @property
    def output_ys(self) -> List[float]:
        return [p.y for p in self.outputs]

    @property
    def total_length(self) -> float:
        in_max = max((p.taper_length for p in self.inputs), default=0.0)
        out_max = max((p.taper_length for p in self.outputs), default=0.0)
        return self.mm_length + in_max + out_max

    def polygons(
        self, stub_length: float = 0.0, output_indices: List[int] = None
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """MMI polygons in device frame (center at x=0).

        output_indices: if given, only draw the output waveguide stub(s) at
        these indices instead of all n_outputs branches (e.g. to build an
        EME cross-section that monitors a single output port in isolation,
        rather than a coupled multi-waveguide supermode). Output tapers are
        always drawn for every branch — they're part of the MMI body's
        transition and stay regardless of which stub is monitored.
        """
        half = self.mm_length / 2
        W = self.mm_width
        in_max = max((p.taper_length for p in self.inputs), default=0.0)
        out_max = max((p.taper_length for p in self.outputs), default=0.0)
        stub_out_idxs = range(len(self.outputs)) if output_indices is None else output_indices

        wg: List[List[Tuple[float, float]]] = [
            [(-half, -W/2), (half, -W/2), (half, W/2), (-half, W/2)]
        ]

        for p in self.inputs:
            ww = p.waveguide_width
            x_taper_start = -half - p.taper_length
            t = Taper(length=p.taper_length, width_start=ww, width_end=p.taper_width)
            wg.extend(shift_polys(t.polygons()[_WG_LAYER], -half - p.taper_length/2, p.y))
            extra = in_max - p.taper_length
            if extra > 1e-9:
                wg.append([(x_taper_start-extra, -ww/2+p.y), (x_taper_start, -ww/2+p.y),
                           (x_taper_start, ww/2+p.y), (x_taper_start-extra, ww/2+p.y)])
            if stub_length > 0:
                sl = stub_length
                x_far = -half - in_max
                wg.append([(x_far-sl, -ww/2+p.y), (x_far, -ww/2+p.y),
                           (x_far, ww/2+p.y), (x_far-sl, ww/2+p.y)])

        for i, p in enumerate(self.outputs):
            ww = p.waveguide_width
            x_taper_end = half + p.taper_length
            t = Taper(length=p.taper_length, width_start=p.taper_width, width_end=ww)
            wg.extend(shift_polys(t.polygons()[_WG_LAYER], half + p.taper_length/2, p.y))
            extra = out_max - p.taper_length
            if extra > 1e-9:
                wg.append([(x_taper_end, -ww/2+p.y), (x_taper_end+extra, -ww/2+p.y),
                           (x_taper_end+extra, ww/2+p.y), (x_taper_end, ww/2+p.y)])
            if stub_length > 0 and i in stub_out_idxs:
                sl = stub_length
                x_far = half + out_max
                wg.append([(x_far, -ww/2+p.y), (x_far+sl, -ww/2+p.y),
                           (x_far+sl, ww/2+p.y), (x_far, ww/2+p.y)])
        return {_WG_LAYER: wg}

    def to_gds_eme(
        self, straight_length: float, oxide_width: float,
        x_start: float = 0.0, output_idx: int = None, **kwargs
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """Build the EME layout. output_idx restricts the drawn output branch(es)
        to a single output index (e.g. output_idx=1 for a top-only output
        stub) instead of all n_outputs — pass None to keep all of them."""
        in_max = max((p.taper_length for p in self.inputs), default=0.0)
        left_edge = -(self.mm_length / 2 + in_max)
        x_off = x_start - left_edge + straight_length
        output_indices = None if output_idx is None else [output_idx]
        wg_polys = [
            [(x + x_off, y) for x, y in verts]
            for verts in self.polygons(stub_length=straight_length, output_indices=output_indices)[_WG_LAYER]
        ]
        total_span = self.total_length + 2 * straight_length
        ox_poly = _oxide_box(x_start, x_start + total_span, oxide_width)
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    # ------------------------------------------------------------------
    # Arbitrary-section interface (HybridSolver)
    # ------------------------------------------------------------------
    #
    # Device frame: x=0 at the MMI body centre (same as polygons()).
    #   x in [-half_body-in_max, -half_body] : input tapers band (width
    #       varies for at least one port everywhere in this band — a port
    #       with a taper_length shorter than in_max=max(taper_length across
    #       inputs) presents plain straight waveguide before its own taper
    #       starts, but the band as a whole is still all-or-nothing)
    #   x in [-half_body, +half_body]        : MMI body (uniform, width W)
    #   x in [+half_body, +half_body+out_max]: output tapers band (same
    #       all-or-nothing convention, out_max=max(taper_length across
    #       outputs))
    # Anything beyond +-(half_body+in_max/out_max) is the (virtual,
    # arbitrary-length) access waveguide region, per-port uniform width.
    #
    # Cuts / stubs may land anywhere except strictly inside a taper band.
    # input_idx / output_idx select which port to use once a range reaches
    # into the input / output taper band or access-waveguide region;
    # irrelevant (ignored) otherwise. Both default to 0, which is exactly
    # right for a single-port side (e.g. the lone input of a 1xN device).

    def _taper_bounds(self):
        half_body = self.mm_length / 2
        in_max = max((p.taper_length for p in self.inputs), default=0.0)
        out_max = max((p.taper_length for p in self.outputs), default=0.0)
        return (-half_body - in_max, -half_body), (half_body, half_body + out_max)

    def _check_cut(self, x: float) -> None:
        (in0, in1), (out0, out1) = self._taper_bounds()
        if in0 < x < in1 or out0 < x < out1:
            raise ValueError(
                f"Cut/stub at x={x:.4f} falls strictly inside a taper band "
                "(cross-section varies there for at least one port) — cuts "
                "must land on a band endpoint or inside the uniform MMI "
                "body / access waveguide."
            )

    def _section_polygons(self, x0: float, x1: float,
                           stub_left: float = 0.0, stub_right: float = 0.0,
                           input_idx: int = 0, output_idx: int = 0):
        """Raw vertex lists (device frame) covering [x0-stub_left, x1+stub_right].

        Draws *every* input/output branch, not just input_idx/output_idx —
        those two only select which branch width_at()/center_at() report for
        monitor sizing (i.e. which one HybridSolver is exciting/measuring);
        the small FDTD box a section occupies is a real physical cross-section
        of the whole device, and every other input/output waveguide the
        device actually has is still sitting right there in it.
        """
        (in0, in1), (out0, out1) = self._taper_bounds()
        W = self.mm_width

        lo, hi = x0 - stub_left, x1 + stub_right
        self._check_cut(lo)
        self._check_cut(hi)

        polys: List[List[Tuple[float, float]]] = []

        def _rect(xa, xb, half_w, yc=0.0):
            if xb - xa > 1e-9:
                polys.append([(xa, -half_w + yc), (xb, -half_w + yc),
                              (xb, half_w + yc), (xa, half_w + yc)])

        # left access waveguides (uniform, one per input port)
        if lo < in0:
            for p in self.inputs:
                _rect(lo, min(hi, in0), p.waveguide_width / 2, yc=p.y)
        # input tapers band (all-or-nothing — never cut mid-band)
        if lo <= in0 and hi >= in1:
            for p in self.inputs:
                x_taper_start = in1 - p.taper_length
                extra = x_taper_start - in0
                if extra > 1e-9:
                    _rect(in0, x_taper_start, p.waveguide_width / 2, yc=p.y)
                t = Taper(length=p.taper_length, width_start=p.waveguide_width, width_end=p.taper_width)
                polys.extend(shift_polys(t.polygons()[_WG_LAYER], x_taper_start + p.taper_length/2, p.y))
        # MMI body (uniform, arbitrary sub-range)
        _rect(max(lo, in1), min(hi, out0), W / 2)
        # output tapers band (all-or-nothing, one per output port)
        if lo <= out0 and hi >= out1:
            for p in self.outputs:
                x_taper_end = out0 + p.taper_length
                t = Taper(length=p.taper_length, width_start=p.taper_width, width_end=p.waveguide_width)
                polys.extend(shift_polys(t.polygons()[_WG_LAYER], x_taper_end - p.taper_length/2, p.y))
                extra = out1 - x_taper_end
                if extra > 1e-9:
                    _rect(x_taper_end, out1, p.waveguide_width / 2, yc=p.y)
        # right access waveguides (uniform, one per output port)
        if hi > out1:
            for p in self.outputs:
                _rect(max(lo, out1), hi, p.waveguide_width / 2, yc=p.y)

        return polys

    def section_prisms(self, x0: float, x1: float, si_medium,
                        stub_left: float = 0.0, stub_right: float = 0.0,
                        input_idx: int = 0, output_idx: int = 0,
                        wg_height: Optional[float] = None, **_kwargs):
        return polygons_to_prisms(
            self._section_polygons(x0, x1, stub_left, stub_right, input_idx, output_idx),
            si_medium, wg_height=wg_height,
        )

    def section_gds(self, x0: float, x1: float, oxide_width: float,
                     stub_left: float = 0.0, stub_right: float = 0.0,
                     input_idx: int = 0, output_idx: int = 0, **_kwargs
                     ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        lo, hi = x0 - stub_left, x1 + stub_right
        x_off = -lo
        wg_polys = [
            [(x + x_off, y) for x, y in verts]
            for verts in self._section_polygons(x0, x1, stub_left, stub_right, input_idx, output_idx)
        ]
        ox_poly = _oxide_box(0.0, hi - lo, oxide_width)
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    def width_at(self, x: float, input_idx: int = 0, output_idx: int = 0, **_kwargs) -> float:
        half = self.mm_length / 2
        if -half <= x <= half:
            return self.mm_width
        if x < -half:
            p = self.inputs[input_idx]
            x0 = -half - p.taper_length
            if x <= x0:
                return p.waveguide_width
            frac = (x - x0) / p.taper_length
            return p.waveguide_width + frac * (p.taper_width - p.waveguide_width)
        p = self.outputs[output_idx]
        x1 = half + p.taper_length
        if x >= x1:
            return p.waveguide_width
        frac = (x - half) / p.taper_length
        return p.taper_width + frac * (p.waveguide_width - p.taper_width)

    def center_at(self, x: float, input_idx: int = 0, output_idx: int = 0, **_kwargs) -> float:
        half_body = self.mm_length / 2
        if x < -half_body:
            return self.inputs[input_idx].y
        if x > half_body:
            return self.outputs[output_idx].y
        return 0.0

    def source_aperture_at(self, x: float, input_idx: int = 0, output_idx: int = 0,
                            **_kwargs) -> Optional[float]:
        """Widest eigenmode source / norm-monitor window (full y-span, µm),
        centered on center_at(x), that stays clear of every OTHER port's
        core at this cross-section. Opt-in via use_neighbor_aperture — off
        by default so half_domain_aperture (a fixed half-domain rule) keeps
        driving existing single-taper-band devices unchanged; turn this on
        for devices with several closely-spaced ports (e.g. one sitting
        exactly on the y=0 axis) that the fixed rule can't isolate."""
        if not self.use_neighbor_aperture:
            return None
        (in0, in1), (out0, out1) = self._taper_bounds()
        if in1 <= x <= out0:
            return None  # MMI body — no single port to isolate
        if x < in1:
            idx, ports = input_idx, self.inputs
        else:
            idx, ports = output_idx, self.outputs
        if len(ports) < 2:
            return None
        y_k = ports[idx].y
        w = lambda j: self.width_at(x, input_idx=j) if x < in1 else self.width_at(x, output_idx=j)
        half = min(abs(ports[j].y - y_k) - w(j) / 2 for j in range(len(ports)) if j != idx)
        return 2.0 * half

    def eme_lengths(self, straight_length: float) -> List[float]:
        in_max = max((p.taper_length for p in self.inputs), default=0.0)
        out_max = max((p.taper_length for p in self.outputs), default=0.0)
        return [straight_length, in_max, self.mm_length, out_max, straight_length]

    @property
    def bbox(self) -> Tuple[float, float]:
        return self.total_length, self.mm_width
