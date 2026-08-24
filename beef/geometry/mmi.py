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

    **Tapered multimode body** (`body_taper_length` / `body_start_width` /
    `body_end_width`, added 2026-08-20, made independently-sided 2026-08-21 —
    `mmi_1x4_tapered_body_validation.md` step 2). The body may taper from
    `body_start_width` at the INPUT face, out to `mm_width` across the middle,
    and back to `body_end_width` at the OUTPUT face, instead of running at one
    constant width:

        |<-Lt->|<---- uniform mm_width ---->|<-Lt->|
        body_start_width                      body_end_width
        (input face)                          (output face)

    `mm_width` therefore always means *the width of the uniform middle*, and
    the two face widths (each defaulting to `None` → `mm_width`) the widths
    where the port tapers attach. Either direction is expressible on either
    side — a face width `< mm_width` means the body widens inward from that
    face, `>` means it narrows inward — and the direction is deliberately NOT
    baked in here; that is a device design choice, to be stated by the device
    script rather than assumed by the geometry.

    **The two faces are independent on purpose, and the 1x4 is why.** A real
    device's port fan-out constrains each face separately: the 1x4's four
    output tapers need 11.5 µm of its 12.0 µm output face, leaving 4.2% of
    room there, while its single on-axis input port needs 2.2 µm and leaves
    81.7%. A single symmetric face width is capped by the tighter side and
    throws away the headroom the other side actually has. Measured directly
    from the production constants — see that doc's own blocker section.

    **Defaults reproduce the straight body exactly** (`body_taper_length=0.0`,
    both face widths `None`): same polygons, same `width_at`, same 5-region
    `eme_lengths`, so every existing device is untouched.
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
        body_taper_length: float = 0.0,
        body_start_width: Optional[float] = None,
        body_end_width: Optional[float] = None,
    ):
        self.mm_length = mm_length
        self.mm_width = mm_width
        self.waveguide_height = waveguide_height
        self.use_neighbor_aperture = use_neighbor_aperture

        if body_taper_length < 0.0:
            raise ValueError(
                f"body_taper_length must be >= 0, got {body_taper_length}")
        if 2 * body_taper_length > mm_length + 1e-9:
            raise ValueError(
                f"body_taper_length={body_taper_length} µm x2 exceeds "
                f"mm_length={mm_length} µm — the two body taper bands would "
                f"overlap, leaving no uniform middle. Use at most "
                f"mm_length/2 = {mm_length / 2} µm."
            )
        for _name, _val in (("body_start_width", body_start_width),
                            ("body_end_width", body_end_width)):
            if _val is not None and _val <= 0.0:
                raise ValueError(f"{_name} must be > 0, got {_val}")
        self.body_taper_length = float(body_taper_length)
        # Resolved, never None, so every consumer can read them unconditionally;
        # `is_body_tapered` is the one place the straight/tapered distinction
        # is decided, rather than each call site re-deriving it.
        self.body_start_width = float(mm_width if body_start_width is None
                                      else body_start_width)
        self.body_end_width = float(mm_width if body_end_width is None
                                    else body_end_width)

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

    @property
    def frame_offset(self) -> float:
        """MMI's frame is centered on the MM body (x=0 at the body's own
        center), not on its bbox — its span is
        [-mm_length/2 - in_max, +mm_length/2 + out_max], which is only
        bbox-centered when in_max == out_max. See DeviceGeometry.frame_offset."""
        in_max = max((p.taper_length for p in self.inputs), default=0.0)
        out_max = max((p.taper_length for p in self.outputs), default=0.0)
        return (in_max - out_max) / 2

    def polygons(
        self, stub_length: float = 0.0, output_indices: List[int] = None,
        input_indices: List[int] = None,
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """MMI polygons in device frame (center at x=0).

        output_indices: if given, only draw the output waveguide stub(s) at
        these indices instead of all n_outputs branches (e.g. to build an
        EME cross-section that monitors a single output port in isolation,
        rather than a coupled multi-waveguide supermode). Output tapers are
        always drawn for every branch — they're part of the MMI body's
        transition and stay regardless of which stub is monitored.

        input_indices: the symmetric restriction on the input side — needed
        for a multi-input device (2x2, etc.) whose input access cross-section
        would otherwise contain several waveguides, so an EME solve there
        excites a two-guide supermode instead of a single input port. Input
        tapers, like output tapers, are always drawn for every branch.
        """
        half = self.mm_length / 2
        in_max = max((p.taper_length for p in self.inputs), default=0.0)
        out_max = max((p.taper_length for p in self.outputs), default=0.0)
        stub_out_idxs = range(len(self.outputs)) if output_indices is None else output_indices
        stub_in_idxs = range(len(self.inputs)) if input_indices is None else input_indices

        wg: List[List[Tuple[float, float]]] = list(self._body_polys(-half, half))

        for i, p in enumerate(self.inputs):
            ww = p.waveguide_width
            x_taper_start = -half - p.taper_length
            t = Taper(length=p.taper_length, width_start=ww, width_end=p.taper_width)
            wg.extend(shift_polys(t.polygons()[_WG_LAYER], -half - p.taper_length/2, p.y))
            extra = in_max - p.taper_length
            if extra > 1e-9:
                wg.append([(x_taper_start-extra, -ww/2+p.y), (x_taper_start, -ww/2+p.y),
                           (x_taper_start, ww/2+p.y), (x_taper_start-extra, ww/2+p.y)])
            if stub_length > 0 and i in stub_in_idxs:
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
        x_start: float = 0.0, output_idx: int = None, input_idx: int = None, **kwargs
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """Build the EME layout. output_idx restricts the drawn output branch(es)
        to a single output index (e.g. output_idx=1 for a top-only output
        stub) instead of all n_outputs — pass None to keep all of them.
        input_idx does the same on the input side (needed for a multi-input
        device so the EME input cross-section is one guide, not a supermode)."""
        in_max = max((p.taper_length for p in self.inputs), default=0.0)
        left_edge = -(self.mm_length / 2 + in_max)
        x_off = x_start - left_edge + straight_length
        output_indices = None if output_idx is None else [output_idx]
        input_indices = None if input_idx is None else [input_idx]
        wg_polys = [
            [(x + x_off, y) for x, y in verts]
            for verts in self.polygons(stub_length=straight_length,
                                       output_indices=output_indices,
                                       input_indices=input_indices)[_WG_LAYER]
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

    @property
    def is_body_tapered(self) -> bool:
        """True when the multimode body's own width varies along x.

        The single place the straight/tapered distinction is decided — a
        zero-length band, or BOTH face widths equal to the middle, is a
        *straight* body and must stay byte-identical to a device built
        without these parameters at all. Either face differing alone is
        enough to make it tapered (the 1x4's asymmetric profile is exactly
        that case: input face narrowed, output face left at mm_width)."""
        return (self.body_taper_length > 1e-12
                and (abs(self.body_start_width - self.mm_width) > 1e-12
                     or abs(self.body_end_width - self.mm_width) > 1e-12))

    def _body_taper_bounds(self):
        """The two body taper bands as ((in0, in1), (out0, out1)), or None
        for a straight body. Distinct from `_taper_bounds()`, which is the
        PORT tapers outside the body."""
        if not self.is_body_tapered:
            return None
        half = self.mm_length / 2
        Lt = self.body_taper_length
        return (-half, -half + Lt), (half - Lt, half)

    def body_width_at(self, x: float) -> float:
        """The multimode body's own width at device-frame x (µm).

        Only meaningful inside [-mm_length/2, +mm_length/2]; callers clamp
        first. Constant `mm_width` for a straight body."""
        if not self.is_body_tapered:
            return self.mm_width
        (bin0, bin1), (bout0, bout1) = self._body_taper_bounds()
        W, Lt = self.mm_width, self.body_taper_length
        Ws, We = self.body_start_width, self.body_end_width
        if x <= bin0:
            return Ws
        if x < bin1:
            return Ws + (x - bin0) / Lt * (W - Ws)
        if x <= bout0:
            return W
        if x < bout1:
            return W + (x - bout0) / Lt * (We - W)
        return We

    def _body_polys(self, lo: float, hi: float) -> List[List[Tuple[float, float]]]:
        """Body polygons over the drawn range [lo, hi], clipped to the body.

        Shared by `polygons()` (whole device) and `_section_polygons()` (a
        cascade section) so the two cannot disagree about what the body looks
        like — the same single-implementation reasoning behind
        `assert_uniform_cross_section()` and `section_fdtd_layout()`.

        Taper bands are drawn all-or-nothing, matching the port tapers'
        existing convention; `_check_cut()` is what makes that safe, since it
        refuses a cut landing strictly inside one.
        """
        half = self.mm_length / 2
        a, b = max(lo, -half), min(hi, half)
        if b - a <= 1e-9:
            return []

        def _quad(xa, xb, wa, wb):
            return [(xa, -wa / 2), (xb, -wb / 2), (xb, wb / 2), (xa, wa / 2)]

        if not self.is_body_tapered:
            # Byte-identical to the pre-2026-08-20 single rectangle.
            return [_quad(a, b, self.mm_width, self.mm_width)]

        (bin0, bin1), (bout0, bout1) = self._body_taper_bounds()
        W = self.mm_width
        Ws, We = self.body_start_width, self.body_end_width
        polys: List[List[Tuple[float, float]]] = []
        if lo <= bin0 and hi >= bin1:
            polys.append(_quad(bin0, bin1, Ws, W))
        mid_a, mid_b = max(a, bin1), min(b, bout0)
        if mid_b - mid_a > 1e-9:
            polys.append(_quad(mid_a, mid_b, W, W))
        if lo <= bout0 and hi >= bout1:
            polys.append(_quad(bout0, bout1, W, We))
        return polys

    def _check_cut(self, x: float) -> None:
        (in0, in1), (out0, out1) = self._taper_bounds()
        if in0 < x < in1 or out0 < x < out1:
            raise ValueError(
                f"Cut/stub at x={x:.4f} falls strictly inside a taper band "
                "(cross-section varies there for at least one port) — cuts "
                "must land on a band endpoint or inside the uniform MMI "
                "body / access waveguide."
            )
        body_bands = self._body_taper_bounds()
        if body_bands is not None:
            W = self.mm_width
            # Only a band whose width ACTUALLY varies is off-limits. With an
            # asymmetric profile one face commonly sits at mm_width (the 1x4's
            # output face, pinned by its port fan-out), making that band a
            # plain rectangle — refusing cuts there would forbid section
            # boundaries on a genuinely uniform cross-section, which is the
            # opposite of this guard's purpose.
            for (b0, b1), face_w in zip(body_bands,
                                        (self.body_start_width, self.body_end_width)):
                if abs(face_w - W) <= 1e-12:
                    continue
                if b0 < x < b1:
                    raise ValueError(
                        f"Cut/stub at x={x:.4f} falls strictly inside a BODY "
                        f"taper band ([{b0:.4f}, {b1:.4f}], width "
                        f"{face_w} -> {W} µm) — the multimode body's own width "
                        "varies there, so a section boundary cannot sit on a "
                        "single cross-section. Cuts must land on a band "
                        "endpoint or inside the uniform middle."
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
        # MMI body — uniform middle, optionally tapered at both faces
        # (`_body_polys` is a no-op outside [in1, out0] via its own clamping,
        # so passing the full [lo, hi] here is correct, not just convenient).
        polys.extend(self._body_polys(lo, hi))
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
            return self.body_width_at(x)
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
