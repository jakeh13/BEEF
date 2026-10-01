# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Georgia Tech Research Corporation
"""UniformCrossSection: N parallel straight strips, constant along x.

Why this exists. The ring resonator's EME half builds its two cross-sections
from raw polygons and calls `mw.compute_modes_meep` directly — the bare single
strip in `solve_neff()` and the coupled bus+ring pair in
`solve_body_supermodes()` (see examples/ring_resonator/python/
ring_resonator_hybrid_cascaded.py). Neither goes through a DeviceGeometry, and
`RingResonatorCouplerSection.to_gds_eme()` deliberately raises: it covers the
mixed-orientation coupler section, never a cross-section.

The consequence: there is **no geometry object both engines can consume**, so
"engine agreement" and "mode identity" are not merely unrun, they are
unrunnable — you cannot compare an EMESolver against an FDTDSolver on a
geometry only one of them can see. This class is that shared object for the
uniform cross-sections, which is what both of those checks actually need.

It is also the simplest possible DeviceGeometry, and deliberately so: being
longitudinally uniform, **every cut is legal**. There is no `_check_cut` here,
unlike Taper (varies everywhere), SBend / DirectionalCoupler (no interior bend
cut) or RingResonatorCouplerSection (no cut inside the arc). A geometry whose
cross-section is a function of nothing is the one place a cut-position
invariance sweep has a known-exact answer.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from .base import DeviceGeometry, _oxide_box, polygons_to_prisms

_WG_LAYER = (2, 0)
_OX_LAYER = (1, 0)


class UniformCrossSection(DeviceGeometry):
    """`len(strips)` parallel straight strips, each constant in width and
    y-centre over the whole length.

    Args:
        strips: (y_center, width) per strip, in µm. Order is the strip index
            used by `width_at()`/`center_at()`'s `strip_idx` kwarg. Not sorted
            or de-duplicated — the caller's order is the port order.
        length: longitudinal extent (µm). For a pure mode solve the value is
            irrelevant (the cross-section is what is solved) but it still has
            to be a real number for the section/EME hooks, so it is explicit
            rather than defaulted to something arbitrary.
        waveguide_height: Si core thickness (µm), only used by the 3-D prism
            extrusion path.

    Examples matching the ring's own two cross-sections:

        # solve_neff()'s bare strip
        UniformCrossSection([(0.0, 0.45)], length=1.0)

        # solve_body_supermodes()'s coupled bus+ring pair (gap 0.25)
        UniformCrossSection([(0.0, 0.45), (0.70, 0.45)], length=1.0)
    """

    def __init__(
        self,
        strips: Sequence[Tuple[float, float]],
        length: float,
        waveguide_height: float = 0.22,
    ):
        if not strips:
            raise ValueError("UniformCrossSection needs at least one strip.")
        for i, (yc, w) in enumerate(strips):
            if w <= 0:
                raise ValueError(
                    f"UniformCrossSection strip {i} has width {w}; widths must be > 0."
                )
        if length <= 0:
            raise ValueError(f"UniformCrossSection length must be > 0, got {length}.")
        self.strips: List[Tuple[float, float]] = [(float(y), float(w)) for y, w in strips]
        self.length = float(length)
        self.waveguide_height = waveguide_height

    # ------------------------------------------------------------------
    @property
    def total_length(self) -> float:
        return self.length

    @property
    def bbox(self) -> Tuple[float, float]:
        lo = min(yc - w / 2 for yc, w in self.strips)
        hi = max(yc + w / 2 for yc, w in self.strips)
        return self.length, hi - lo

    def _rects(self, xa: float, xb: float) -> List[List[Tuple[float, float]]]:
        """One rectangle per strip spanning [xa, xb]. Uniform, so this is the
        only shape function this class needs."""
        if xb - xa <= 1e-12:
            return []
        return [
            [(xa, yc - w / 2), (xb, yc - w / 2), (xb, yc + w / 2), (xa, yc + w / 2)]
            for yc, w in self.strips
        ]

    def polygons(self, stub_length: float = 0.0
                 ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        half = self.length / 2
        return {_WG_LAYER: self._rects(-half - stub_length, half + stub_length)}

    def to_gds_eme(self, straight_length: float, oxide_width: float,
                   x_start: float = 0.0, **_kwargs
                   ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        total = self.length + 2 * straight_length
        return {
            _WG_LAYER: self._rects(x_start, x_start + total),
            _OX_LAYER: [_oxide_box(x_start, x_start + total, oxide_width)],
        }

    # ------------------------------------------------------------------
    # Arbitrary-section interface. No _check_cut: uniform in x, so any cut and
    # any stub is legal by construction -- that is the whole point of this
    # class (see module docstring).
    # ------------------------------------------------------------------
    def _section_polygons(self, x0: float, x1: float,
                          stub_left: float = 0.0, stub_right: float = 0.0):
        return self._rects(x0 - stub_left, x1 + stub_right)

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
        wg = [[(x - lo, y) for x, y in verts]
              for verts in self._section_polygons(x0, x1, stub_left, stub_right)]
        return {_WG_LAYER: wg, _OX_LAYER: [_oxide_box(0.0, hi - lo, oxide_width)]}

    # ------------------------------------------------------------------
    def _strip(self, strip_idx: int) -> Tuple[float, float]:
        try:
            return self.strips[strip_idx]
        except IndexError:
            raise ValueError(
                f"UniformCrossSection has {len(self.strips)} strip(s); "
                f"strip_idx={strip_idx} is out of range."
            ) from None

    def width_at(self, x: float, strip_idx: int = 0, **_kwargs) -> float:
        """Width of strip `strip_idx`, independent of x.

        Multi-strip geometries have no single width, so which one is reported
        is the caller's choice — same pattern as DirectionalCoupler's
        input_idx/output_idx selecting which arm width_at() describes.
        """
        return self._strip(strip_idx)[1]

    def center_at(self, x: float, strip_idx: int = 0, **_kwargs) -> float:
        return self._strip(strip_idx)[0]

    def source_aperture_at(self, x: float, strip_idx: int = 0, **_kwargs
                           ) -> Optional[float]:
        """Widest launch window on `strip_idx` that stays clear of its nearest
        neighbour, as a full y-span centred on `center_at()`.

        Halfway to the nearest neighbour's CENTRELINE. For equal widths that is
        the gap midpoint, so the window cannot contain any neighbouring core;
        for unequal widths it is conservative on the wide-neighbour side, which
        is the safe direction. Returns None for a single strip — nothing to
        avoid, so no geometry-imposed restriction (matching the base-class
        contract).

        This is what a coupled pair needs instead of
        FDTDSolver.half_domain_aperture: that rule splits the domain at y=0 and
        so assumes the excited guide sits on one side of the axis, which is not
        true of a pair whose one strip IS on the axis.
        """
        if len(self.strips) == 1:
            return None
        yc, _ = self._strip(strip_idx)
        nearest = min(abs(other_y - yc)
                      for j, (other_y, _w) in enumerate(self.strips) if j != strip_idx)
        return nearest          # 2 * (nearest / 2)

    def eme_lengths(self, straight_length: float) -> List[float]:
        return [straight_length, self.length, straight_length]
