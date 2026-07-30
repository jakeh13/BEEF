"""CompoundGeometry: concatenates DeviceGeometry devices (and plain waveguide
connectors) end-to-end along x into a single DeviceGeometry.

HybridSolver.solve()/solve_section() take exactly one geometry object for an
entire cascade (section boundaries and geometry_kwargs are the only
per-section variation — see beef.solvers.hybrid). Building a device out of
several independent sub-devices (e.g. two MMIs joined by a straight
waveguide) therefore requires one composite geometry that stitches them
together and forwards every DeviceGeometry hook to whichever sub-geometry
owns the queried x — that's what this module provides.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

from ._compose import shift_polys, shift_prisms
from .base import DeviceGeometry, _oxide_box, polygons_to_prisms

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


@dataclass
class WaveguideRun:
    """Uniform straight multi-port waveguide connector, safe to cut anywhere.

    ports: (y_center, width) per parallel waveguide carried through the run.
    """
    length: float
    ports: List[Tuple[float, float]]


Segment = Union[DeviceGeometry, WaveguideRun]


class CompoundGeometry(DeviceGeometry):
    """Concatenates `segments` (DeviceGeometry and/or WaveguideRun instances)
    end-to-end along x, left to right, into a single DeviceGeometry.

    Each segment keeps its own local frame (x=0 at its longitudinal center);
    CompoundGeometry places them back-to-back and translates every query
    (width_at/center_at/section_prisms/section_gds/...) into whichever
    segment's local frame owns that x.
    """

    def __init__(self, segments: List[Segment]):
        if len(segments) < 1:
            raise ValueError("CompoundGeometry needs at least 1 segment.")
        self.segments = list(segments)

    @staticmethod
    def _seg_length(seg: Segment) -> float:
        return seg.length if isinstance(seg, WaveguideRun) else seg.total_length

    @property
    def total_length(self) -> float:
        return sum(self._seg_length(s) for s in self.segments)

    def _placements(self):
        """Yield (segment, seg_lo, seg_hi) in device frame, left to right."""
        x = -self.total_length / 2
        for seg in self.segments:
            L = self._seg_length(seg)
            yield seg, x, x + L
            x += L

    def _locate(self, x: float, tol: float = 1e-6):
        for seg, lo, hi in self._placements():
            if lo - tol <= x <= hi + tol:
                return seg, lo, hi
        raise ValueError(f"x={x:.4f} outside CompoundGeometry (total_length={self.total_length:.4f})")

    def _overlaps(self, x0: float, x1: float, stub_left: float, stub_right: float):
        """(segment, seg_lo, seg_hi, overlap_lo, overlap_hi, stub_left, stub_right)
        for every segment overlapping [x0-stub_left, x1+stub_right], with the
        outer stub folded into only the first/last overlapping segment.

        The returned per-segment stubs are only the part of each stub that
        falls OUTSIDE that segment. This window is already stub-expanded
        (`lo`/`hi` below), and `ov_lo`/`ov_hi` are clipped against it — so
        wherever a stub does not actually reach past the segment's own edge,
        `ov_lo`/`ov_hi` already cover it and handing the raw stub down to
        the segment on top of that would apply it a second time. That was a
        real bug (drawn geometry 2x the stub too long, or a spurious
        "cut falls strictly inside a taper" from the doubled coordinate)
        whenever the padded window landed fully inside one segment with
        nothing to clip it back; see
        tests/test_geometry.py::TestCompoundGeometryStubIsAppliedExactlyOnce.

        The resulting invariant, which both callers' offset arithmetic
        relies on: ``ov_lo - stub_left == lo`` for the first overlapping
        segment and ``ov_hi + stub_right == hi`` for the last, in every
        case — clipped or not.
        """
        lo, hi = x0 - stub_left, x1 + stub_right
        raw = [
            (seg, seg_lo, seg_hi, max(lo, seg_lo), min(hi, seg_hi))
            for seg, seg_lo, seg_hi in self._placements()
        ]
        raw = [o for o in raw if o[4] - o[3] > 1e-9]
        result = []
        for i, (seg, seg_lo, seg_hi, ov_lo, ov_hi) in enumerate(raw):
            sl = max(0.0, seg_lo - lo) if i == 0 else 0.0
            sr = max(0.0, hi - seg_hi) if i == len(raw) - 1 else 0.0
            result.append((seg, seg_lo, seg_hi, ov_lo, ov_hi, sl, sr))
        return result, lo, hi

    # ------------------------------------------------------------------
    # Whole-device polygons (plain FDTD / geometry preview) — stub_length
    # not supported here (no generic way to know per-port stub width for an
    # arbitrary segment type); use to_gds_eme()/section_gds() instead, which
    # get their stub width from each end segment's own ports.
    # ------------------------------------------------------------------

    def polygons(self, stub_length: float = 0.0) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        if stub_length > 0:
            raise NotImplementedError(
                "CompoundGeometry.polygons(stub_length>0) is not supported — "
                "use to_gds_eme()/section_gds(), whose stub handling reads "
                "port width from the end segments themselves."
            )
        wg: List[List[Tuple[float, float]]] = []
        for seg, seg_lo, seg_hi in self._placements():
            if isinstance(seg, WaveguideRun):
                for (y, w) in seg.ports:
                    wg.append([(seg_lo, y-w/2), (seg_hi, y-w/2), (seg_hi, y+w/2), (seg_lo, y+w/2)])
            else:
                center = (seg_lo + seg_hi) / 2
                wg.extend(shift_polys(seg.polygons(stub_length=0.0)[_WG_LAYER], center))
        return {_WG_LAYER: wg}

    @staticmethod
    def _end_ports(seg: Segment, side: str, output_idx: Optional[int] = None
                    ) -> List[Tuple[float, float]]:
        """(y, width) pairs to extend as an access-waveguide stub at one end."""
        if isinstance(seg, WaveguideRun):
            return list(seg.ports)
        half = seg.total_length / 2
        eps = 1e-6
        if side == "input":
            ys = getattr(seg, "input_ys", None)
            if ys is None:
                return [(seg.center_at(-half + eps, input_idx=0), seg.width_at(-half + eps, input_idx=0))]
            return [(y, seg.width_at(-half + eps, input_idx=j)) for j, y in enumerate(ys)]
        ys = getattr(seg, "output_ys", None)
        if ys is None:
            ys = [seg.center_at(half - eps, output_idx=0)]
        idxs = range(len(ys)) if output_idx is None else [output_idx]
        return [(ys[j], seg.width_at(half - eps, output_idx=j)) for j in idxs]

    def to_gds_eme(
        self, straight_length: float, oxide_width: float,
        x_start: float = 0.0, output_idx: Optional[int] = None, **kwargs
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        placements = list(self._placements())
        last_seg = placements[-1][0]
        wg: List[List[Tuple[float, float]]] = []
        for seg, seg_lo, seg_hi in placements:
            if isinstance(seg, WaveguideRun):
                for (y, w) in seg.ports:
                    wg.append([(seg_lo, y-w/2), (seg_hi, y-w/2), (seg_hi, y+w/2), (seg_lo, y+w/2)])
            else:
                seg_kwargs = dict(kwargs)
                if seg is last_seg and output_idx is not None:
                    seg_kwargs["output_idx"] = output_idx
                sub = seg.to_gds_eme(0.0, oxide_width, x_start=0.0, **seg_kwargs)
                wg.extend(shift_polys(sub[_WG_LAYER], seg_lo))

        total = self.total_length
        x_off = x_start + total / 2 + straight_length
        wg_polys = shift_polys(wg, x_off)

        if straight_length > 0:
            sl = straight_length
            first_lo = placements[0][1]
            last_hi = placements[-1][2]
            xa = first_lo + x_off
            for y, w in self._end_ports(placements[0][0], "input"):
                wg_polys.append([(xa-sl, y-w/2), (xa, y-w/2), (xa, y+w/2), (xa-sl, y+w/2)])
            xb = last_hi + x_off
            for y, w in self._end_ports(last_seg, "output", output_idx=output_idx):
                wg_polys.append([(xb, y-w/2), (xb+sl, y-w/2), (xb+sl, y+w/2), (xb, y+w/2)])

        ox_poly = _oxide_box(x_start, x_start + total + 2 * straight_length, oxide_width)
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    # ------------------------------------------------------------------
    # Arbitrary-section interface (HybridSolver)
    # ------------------------------------------------------------------
    # A section spanning multiple segments (e.g. the tail of one MMI's
    # body, its output tapers, a connecting WaveguideRun, and the next
    # MMI's input tapers) is split across every overlapping segment and
    # each piece delegated to that segment's own section_prisms/section_gds
    # — which already raises ValueError on an interior cut inside a taper/
    # bend, so CompoundGeometry inherits that invariant for free.

    def section_prisms(self, x0: float, x1: float, si_medium,
                        stub_left: float = 0.0, stub_right: float = 0.0,
                        wg_height: Optional[float] = None, **kwargs):
        overlaps, _, _ = self._overlaps(x0, x1, stub_left, stub_right)
        prisms = []
        for seg, seg_lo, seg_hi, ov_lo, ov_hi, sl, sr in overlaps:
            if isinstance(seg, WaveguideRun):
                a, b = ov_lo - sl, ov_hi + sr
                for (y, w) in seg.ports:
                    prisms.extend(polygons_to_prisms(
                        [[(a, y-w/2), (b, y-w/2), (b, y+w/2), (a, y+w/2)]],
                        si_medium, wg_height=wg_height))
            else:
                center = (seg_lo + seg_hi) / 2
                sub = seg.section_prisms(ov_lo - center, ov_hi - center, si_medium,
                                          stub_left=sl, stub_right=sr, wg_height=wg_height, **kwargs)
                prisms.extend(shift_prisms(sub, center))
        return prisms

    def section_gds(self, x0: float, x1: float, oxide_width: float,
                     stub_left: float = 0.0, stub_right: float = 0.0, **kwargs
                     ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        overlaps, lo, hi = self._overlaps(x0, x1, stub_left, stub_right)
        wg_polys: List[List[Tuple[float, float]]] = []
        for seg, seg_lo, seg_hi, ov_lo, ov_hi, sl, sr in overlaps:
            if isinstance(seg, WaveguideRun):
                a, b = ov_lo - sl - lo, ov_hi + sr - lo
                for (y, w) in seg.ports:
                    wg_polys.append([(a, y-w/2), (b, y-w/2), (b, y+w/2), (a, y+w/2)])
            else:
                center = (seg_lo + seg_hi) / 2
                sub_dict = seg.section_gds(ov_lo - center, ov_hi - center, oxide_width,
                                            stub_left=sl, stub_right=sr, **kwargs)
                offset = (ov_lo - sl) - lo
                wg_polys.extend(shift_polys(sub_dict.get(_WG_LAYER, []), offset))
        ox_poly = _oxide_box(0.0, hi - lo, oxide_width)
        return {_WG_LAYER: wg_polys, _OX_LAYER: [ox_poly]}

    def width_at(self, x: float, **kwargs) -> float:
        seg, lo, hi = self._locate(x)
        if isinstance(seg, WaveguideRun):
            idx = kwargs.get("output_idx", kwargs.get("input_idx", 0))
            return seg.ports[idx][1]
        center = (lo + hi) / 2
        return seg.width_at(x - center, **kwargs)

    def center_at(self, x: float, **kwargs) -> float:
        seg, lo, hi = self._locate(x)
        if isinstance(seg, WaveguideRun):
            idx = kwargs.get("output_idx", kwargs.get("input_idx", 0))
            return seg.ports[idx][0]
        center = (lo + hi) / 2
        return seg.center_at(x - center, **kwargs)

    def source_aperture_at(self, x: float, **kwargs) -> Optional[float]:
        seg, lo, hi = self._locate(x)
        if isinstance(seg, WaveguideRun):
            return None
        center = (lo + hi) / 2
        return seg.source_aperture_at(x - center, **kwargs)

    def eme_lengths(self, straight_length: float) -> List[float]:
        lengths: List[float] = []
        for seg in self.segments:
            if isinstance(seg, WaveguideRun):
                lengths.append(seg.length)
            else:
                lengths.extend(seg.eme_lengths(0.0)[1:-1])
        return [straight_length] + lengths + [straight_length]

    @property
    def input_ys(self) -> Optional[List[float]]:
        return getattr(self.segments[0], "input_ys", None)

    @property
    def output_ys(self) -> Optional[List[float]]:
        return getattr(self.segments[-1], "output_ys", None)

    @property
    def bbox(self) -> Tuple[float, float]:
        widths = []
        for seg in self.segments:
            if isinstance(seg, WaveguideRun):
                widths.append(max((w for _, w in seg.ports), default=0.0))
            else:
                widths.append(seg.bbox[1])
        return self.total_length, max(widths, default=0.0)
