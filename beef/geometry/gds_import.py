"""GDS file import: load raw polygons from an on-disk .gds layout and wrap
them as a DeviceGeometry for EME.

BEEF's other geometries (MMI, Taper, ...) build their polygons themselves,
parametrically, in Python. This module instead reads polygons already drawn
by an external tool (gdspy, gdsfactory, ...) and adapts them to the same
Dict[(layer, datatype), List[polygon]] convention every other geometry uses,
so they flow through extrude_polygons()/EMESolver unchanged.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from .base import (
    DeviceGeometry, _oxide_box, assert_uniform_cross_section, polygons_to_prisms,
)

_OX_LAYER = (1, 0)


def load_gds_polygons(
    gds_path: str,
    layers: Optional[Sequence[Tuple[int, int]]] = None,
    cell_name: Optional[str] = None,
) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
    """Read a GDS file and return its polygons as BEEF's standard
    Dict[(layer, datatype), List[polygon]] dict (polygon = list of (x, y)
    vertex tuples, µm — gdstk reports coordinates in µm regardless of the
    file's own unit/precision, so no manual unit conversion is needed).

    Args:
        gds_path: path to the .gds file.
        layers: if given, only these (layer, datatype) tuples are extracted
            (everything else in the file is ignored). Extracts every layer
            present otherwise.
        cell_name: name of the cell to read; defaults to the library's sole
            top-level cell (raises if there is more than one and no name is
            given).

    Returns:
        Dict[(layer, datatype), List[List[(x, y)]]], vertices in the GDS
        file's own coordinate frame (no shift/centering applied).
    """
    import gdstk

    lib = gdstk.read_gds(str(gds_path))

    if cell_name is not None:
        matches = [c for c in lib.cells if c.name == cell_name]
        if not matches:
            raise ValueError(f"cell {cell_name!r} not found in {gds_path}")
        cell = matches[0]
    else:
        top = lib.top_level()
        if len(top) != 1:
            names = [c.name for c in top]
            raise ValueError(
                f"{gds_path} has {len(top)} top-level cells {names} — "
                "pass cell_name= to disambiguate"
            )
        cell = top[0]

    polys = cell.get_polygons(apply_repetitions=True, include_paths=True, depth=None)

    out: Dict[Tuple[int, int], List[List[Tuple[float, float]]]] = {}
    for p in polys:
        layer = (p.layer, p.datatype)
        if layers is not None and layer not in layers:
            continue
        out.setdefault(layer, []).append([(float(x), float(y)) for x, y in p.points])
    return out


class GdsGeometry(DeviceGeometry):
    """DeviceGeometry backed by polygons imported directly from a GDS file.

    A straight-through (single path, no branching) device: the imported
    polygons already contain whatever access waveguides the GDS layout
    drew, so unlike Taper/MMI/SBend (which synthesize their own stub
    extensions from the `straight_length` argument on demand),
    to_gds_eme() places the imported polygons as-is and appends nothing.

    `region_lengths` is a caller-supplied breakdown of the imported
    geometry's own total span (e.g. [stub, taper, stub]), used only as the
    return value of eme_lengths() — a convenience for building the
    `lengths=` list passed to EMESolver. It must sum to the bounding-box
    span of `core_layers`' polygons (checked at construction).

    Since 2026-08-20 this class also implements the arbitrary-section hooks
    `HybridSolver` needs (`width_at`/`center_at`/`source_aperture_at`,
    `section_prisms`/`section_gds`), so an imported layout can occupy a
    section of a cascade rather than only a whole-device single-solver run.
    Every one of them is answered by intersecting the imported polygons with
    a vertical line — there is no parametric taper/body/branch description
    to consult — which is both the point and the source of this class's two
    stated limitations: the cut-validity check is sampled rather than exact
    (see `_check_cut`), and branch selection is positional (see
    `_branch_at`). See
    validation_markdown/resolved/gds_geometry_hybrid_cascade_validation.md.
    """

    def __init__(
        self,
        gds_path: str,
        core_layers: Sequence[Tuple[int, int]],
        region_lengths: Sequence[float],
        cell_name: Optional[str] = None,
        input_width: Optional[float] = None,
        input_y: float = 0.0,
        access_extension: bool = False,
        neighbor_aperture: bool = False,
        split_x: float = 0.0,
        cut_check_span: float = 0.05,
        cut_check_samples: int = 9,
        cut_check_tol: float = 1e-3,
        boolean_precision: float = 1e-7,
    ):
        """
        Args:
            input_width: if given, the width (µm) of a straight waveguide the
                GDS itself draws at the geometry's own LEFT edge — enables
                polygons(stub_length>0) to extend that one waveguide further
                left by stub_length (see polygons()'s docstring for why this
                is one-sided, not the base class's usual symmetric stub).
                None (default): stub_length stays a no-op, as before.
            input_y: y-centre (µm) of that same input waveguide.
            access_extension: synthesize straight access waveguides beyond
                either end facet, by extruding the cross-section measured
                there, whenever section_prisms()/section_gds()/to_gds_eme()
                are asked for geometry outside the imported polygons' own
                x-span. Off by default, which keeps this class's stated
                "places the imported polygons as-is and appends nothing"
                behaviour — including to_gds_eme()'s long-standing
                indifference to `straight_length`.

                Turn it ON for a GDS holding the bare device, whose ports
                end at their own facets: a HybridSolver end section ALWAYS
                asks for a stub past the physical port (`hybrid.py:369-370`
                takes it from `solver.waveguide_length`), and FDTDSolver.
                solve_section() adds another `dpml` on top so the waveguide
                runs through the absorbing layer rather than presenting it
                an abrupt core→cladding step. Leave it OFF for a GDS that
                draws its own access waveguides; a section reaching past
                the drawn geometry then raises instead of silently ending
                the core in mid-air. This generalizes to both ends — and to
                every branch at the facet, not one — the one-sided
                input_width/input_y mechanism above, which stays as-is for
                polygons()/solve_reflection()'s whole-device path.
            neighbor_aperture: opt in to source_aperture_at() reporting a
                measured neighbour-clearance launch window, mirroring MMI's
                own `use_neighbor_aperture`. Off by default so
                `half_domain_aperture` keeps driving existing runs.
            split_x: device-frame x separating the branches `input_idx`
                indexes from those `output_idx` indexes (see _branch_at()).
                Defaults to 0.0 — the imported bbox's own centre, since
                __init__ centres it there — which is right whenever the
                device's midpoint sits in a region with a single branch
                (an MMI body, a coupler's interaction region). Move it for
                a layout where it does not.
            cut_check_span, cut_check_samples, cut_check_tol: the sampled
                stand-in for a parametric cut-validity check — see
                _check_cut(), including what sampling cannot see.
            boolean_precision: µm grid gdstk's clipping boolean snaps to
                (section_prisms/section_gds). 1e-7 µm = 0.1 pm, four orders
                finer than gdstk's own 1 nm default, which would visibly
                move a slanted taper edge.
        """
        self.gds_path = str(gds_path)
        self.core_layers = list(core_layers)
        self.region_lengths = list(region_lengths)
        self.input_width = input_width
        self.input_y = input_y
        self.access_extension = bool(access_extension)
        self.neighbor_aperture = bool(neighbor_aperture)
        self.split_x = float(split_x)
        self.cut_check_span = float(cut_check_span)
        self.cut_check_samples = int(cut_check_samples)
        self.cut_check_tol = float(cut_check_tol)
        self._boolean_precision = float(boolean_precision)

        raw = load_gds_polygons(self.gds_path, layers=self.core_layers, cell_name=cell_name)

        xs = [x for layer in self.core_layers for poly in raw.get(layer, []) for x, _y in poly]
        if not xs:
            raise ValueError(
                f"no polygons found on layers {self.core_layers} in {gds_path}"
            )
        x_min, x_max = min(xs), max(xs)
        span = x_max - x_min
        lengths_sum = sum(self.region_lengths)
        if abs(span - lengths_sum) > 1e-6:
            raise ValueError(
                f"region_lengths sum to {lengths_sum:.4f} µm but the imported "
                f"geometry spans {span:.4f} µm on layers {self.core_layers} "
                f"({gds_path})"
            )

        self._span = span
        x_off = -(x_min + x_max) / 2.0  # center the bbox at x=0 (device frame)
        self._polys = {
            layer: [[(x + x_off, y) for x, y in poly] for poly in polys]
            for layer, polys in raw.items()
        }

        # Device-frame bbox, and the two end facets' cross-sections. Both
        # facets are measured once here rather than per query: they are what
        # width_at()/center_at() report beyond the drawn geometry and what
        # access_extension extrudes, so they must not drift between the two.
        self._x_min = min(x for polys in self._polys.values()
                          for poly in polys for x, _y in poly)
        self._x_max = max(x for polys in self._polys.values()
                          for poly in polys for x, _y in poly)
        all_ys = [y for polys in self._polys.values()
                  for poly in polys for _x, y in poly]
        self._y_min, self._y_max = min(all_ys), max(all_ys)

        self._facet_left = self._scan(self._x_min)
        self._facet_right = self._scan(self._x_max)
        self._facet_left_by_layer = {
            layer: self._scan_layer(layer, self._x_min) for layer in self._polys
        }
        self._facet_right_by_layer = {
            layer: self._scan_layer(layer, self._x_max) for layer in self._polys
        }

    @property
    def total_length(self) -> float:
        return self._span

    def eme_lengths(self, straight_length: float) -> List[float]:
        return list(self.region_lengths)

    def polygons(
        self, stub_length: float = 0.0
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """Waveguide polygons as imported from the GDS, centered at x=0 (see
        __init__), optionally with the LEFT (input) waveguide extended by
        stub_length.

        Deliberately ONE-SIDED, unlike every other geometry's polygons()
        (which extends symmetrically): this class's own device has a real,
        known-cross-section waveguide only on its input side (the GDS's own
        straight stub) — the far/right edge is an arbitrary cut through a
        freeform region with no well-defined width_at() there, so extending
        it the same way would draw a waveguide where none physically exists.
        See A3/step 8 in validation_markdown/high_contrast_grating_validation.md
        for why this matters: without this extension, the input waveguide
        ends in mid-air short of (or right at) the PML, and that clean-cut
        Si/oxide facet reflects part of the guided mode right there —
        contaminating exactly the reflection measurement
        FDTDSolver.solve_reflection() exists to make. Extending the
        waveguide THROUGH the PML (stub_length = edge_clearance + dpml, the
        distance from the geometry's own edge to the simulation cell's true
        boundary) lets the PML absorb it as an ordinary waveguide mode
        instead.

        stub_length > 0 requires input_width (set at construction) — without
        it there is no known cross-section to extrude, and this raises
        rather than silently drawing nothing (the pre-2026-08-04 behavior,
        which is still what stub_length=0 or input_width=None gives).
        """
        out = {layer: [list(poly) for poly in polys] for layer, polys in self._polys.items()}
        if stub_length > 0:
            if self.input_width is None:
                raise ValueError(
                    "polygons(stub_length>0) needs input_width set at "
                    "GdsGeometry construction — there is no way to infer a "
                    "freeform region's cross-section automatically"
                )
            x_left = -self._span / 2.0
            y0, y1 = self.input_y - self.input_width / 2.0, self.input_y + self.input_width / 2.0
            ext = [
                (x_left - stub_length, y0), (x_left, y0),
                (x_left, y1), (x_left - stub_length, y1),
            ]
            core_layer = self.core_layers[0]
            out.setdefault(core_layer, []).append(ext)
        return out

    def to_gds_eme(
        self, straight_length: float, oxide_width: float,
        x_start: float = 0.0, **kwargs
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """Whole-device EME layout, left edge at `x_start`.

        `straight_length` is ignored unless `access_extension` is on — the
        imported polygons normally already contain whatever access
        waveguides the layout drew, so appending more would draw waveguide
        where none exists. With access_extension=True the GDS instead holds
        the bare device, and this appends `straight_length` of each end
        facet's own cross-section on both sides, matching what every
        parametric geometry's to_gds_eme() does and what
        plot_hybrid_geometry() assumes when it places section boundaries at
        `total_length/2 + straight_length`.
        """
        ext = straight_length if (self.access_extension and straight_length > 0) else 0.0
        x_off = x_start + ext + self._span / 2.0
        out = {
            layer: [[(x + x_off, y) for x, y in poly] for poly in polys]
            for layer, polys in self._polys.items()
        }
        if ext > 0:
            for layer in self._polys:
                out[layer].extend(
                    self._rects(x_start, x_start + ext,
                                [(y0, y1) for y0, y1 in self._facet_left_by_layer[layer]])
                )
                out[layer].extend(
                    self._rects(x_start + ext + self._span,
                                x_start + 2 * ext + self._span,
                                [(y0, y1) for y0, y1 in self._facet_right_by_layer[layer]])
                )
        total_span = self._span + 2 * ext
        out[_OX_LAYER] = [_oxide_box(x_start, x_start + total_span, oxide_width)]
        return out

    @property
    def bbox(self) -> Tuple[float, float]:
        return self._span, self._y_max - self._y_min

    # ------------------------------------------------------------------
    # Cross-section queries, read off the polygons (HybridSolver hooks)
    # ------------------------------------------------------------------
    #
    # Device frame: bbox centred at x=0 (see __init__). Every query below
    # answers "what does a vertical line at x cut through?" from geometry
    # alone — there is no parametric taper/body/branch knowledge to consult,
    # which is the whole point: an imported layout has none.
    #
    #   x <  bbox left edge   : the left end facet's cross-section, held
    #   x >  bbox right edge  : the right end facet's cross-section, held
    #   otherwise             : scanline intersection at x
    #
    # Holding the facet cross-section outside the bbox mirrors what every
    # parametric geometry reports in its own (virtual, arbitrary-length)
    # access-waveguide region — MMI.width_at() returns waveguide_width for
    # any x beyond its taper band — and is what makes access_extension's
    # synthesized stubs and the cut check below well-defined out there.

    _X_TOL = 1e-9      # µm; vertex/scanline coincidence tolerance
    _MERGE_TOL = 1e-9  # µm; two y-intervals closer than this are one branch

    @staticmethod
    def _poly_crossings(poly: List[Tuple[float, float]], x: float) -> List[float]:
        """y values where the vertical line at `x` crosses `poly`'s boundary.

        Closed edge test (``x_lo - tol <= x <= x_hi + tol``) rather than the
        usual half-open one: the queries that matter most here land exactly
        ON a vertical facet — an end facet, or the seam where a taper's wide
        end abuts an MMI body — and a half-open rule silently drops the
        polygon on one side of it (reporting the taper's 1.3 µm instead of
        the body's 2.8 µm at that seam, say). The tolerance also absorbs the
        sub-femtometre vertex noise a GDS round-trip leaves behind.

        Limitation, stated rather than papered over: at a vertex where the
        boundary *touches* the line without crossing it (a local extremum in
        x), the closed test yields a coincident pair, i.e. a zero-width
        interval — dropped by _intervals_at(). Orthogonal/trapezoidal
        photonic layouts, which is what this class is fed, have no such
        vertices except at end facets, where the pair is the real
        cross-section.
        """
        tol = GdsGeometry._X_TOL
        ys: List[float] = []
        n = len(poly)
        for i in range(n):
            xa, ya = poly[i]
            xb, yb = poly[(i + 1) % n]
            dx = xb - xa
            if abs(dx) <= tol:
                continue  # vertical edge: parallel to the scanline
            lo, hi = (xa, xb) if xa < xb else (xb, xa)
            if lo - tol <= x <= hi + tol:
                t = min(max((x - xa) / dx, 0.0), 1.0)
                ys.append(ya + t * (yb - ya))
        ys.sort()
        return ys

    def _intervals_at(self, x: float) -> List[Tuple[float, float]]:
        """Core-material y-intervals at device-frame `x`, sorted by y and
        merged across layers/polygons — one entry per physically distinct
        branch at this cross-section.

        Merged, not per-polygon: an imported layout draws a body and the
        tapers that abut it as separate, touching polygons, and reporting
        those as separate branches would make width_at() answer 1.3 µm at a
        2.8 µm-wide MMI face. Every core layer contributes, since they are
        all core by this class's own definition.
        """
        if x < self._x_min:
            return self._facet_left
        if x > self._x_max:
            return self._facet_right
        return self._scan(x)

    def _scan(self, x: float) -> List[Tuple[float, float]]:
        """Scanline over every core layer at once — the guided cross-section
        as the mode solver sees it, regardless of which layer drew it."""
        return self._merge(
            [iv for layer in self._polys for iv in self._raw_intervals(layer, x)]
        )

    def _scan_layer(self, layer: Tuple[int, int], x: float) -> List[Tuple[float, float]]:
        """Scanline over ONE layer, for access_extension's per-layer stub
        rectangles — a two-material import must extend each material in its
        own layer, not extrude their union into whichever layer came first."""
        return self._merge(self._raw_intervals(layer, x))

    def _raw_intervals(self, layer: Tuple[int, int], x: float) -> List[Tuple[float, float]]:
        raw: List[Tuple[float, float]] = []
        for poly in self._polys[layer]:
            ys = self._poly_crossings(poly, x)
            for i in range(0, len(ys) - 1, 2):
                if ys[i + 1] - ys[i] > self._MERGE_TOL:
                    raw.append((ys[i], ys[i + 1]))
        return raw

    @classmethod
    def _merge(cls, raw: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
        if not raw:
            return []
        raw = sorted(raw)
        merged = [raw[0]]
        for lo, hi in raw[1:]:
            plo, phi = merged[-1]
            if lo <= phi + cls._MERGE_TOL:
                merged[-1] = (plo, max(phi, hi))
            else:
                merged.append((lo, hi))
        return merged

    def _branch_at(self, x: float, input_idx: int = 0, output_idx: int = 0,
                    branch_idx: Optional[int] = None, **_kwargs) -> Tuple[float, float]:
        """The one y-interval at `x` that width_at()/center_at() report.

        Branch selection follows MMI's convention so the same
        `geometry_kwargs` a parametric run already passes work unchanged:
        `input_idx` selects among the branches left of `split_x`,
        `output_idx` among those right of it; `branch_idx` overrides both.
        A cross-section with a single branch is unambiguous, so both are
        ignored there — again matching MMI, whose width_at() returns
        mm_width inside the body regardless of input_idx/output_idx.
        """
        ivs = self._intervals_at(x)
        if not ivs:
            raise ValueError(
                f"{type(self).__name__}: no core material at device-frame "
                f"x={x:.4g} µm (bbox spans [{self._x_min:.4g}, "
                f"{self._x_max:.4g}]) — nothing to report a cross-section for."
            )
        if len(ivs) == 1:
            return ivs[0]
        idx = branch_idx if branch_idx is not None else (
            input_idx if x < self.split_x else output_idx
        )
        if not 0 <= idx < len(ivs):
            side = "input_idx" if (branch_idx is None and x < self.split_x) else (
                "branch_idx" if branch_idx is not None else "output_idx")
            raise ValueError(
                f"{type(self).__name__}: {side}={idx} is out of range at "
                f"device-frame x={x:.4g} µm, which cuts {len(ivs)} branch(es) "
                f"at y-centres "
                f"{[round((a + b) / 2, 4) for a, b in ivs]} (ordered by y)."
            )
        return ivs[idx]

    def width_at(self, x: float, **kwargs) -> float:
        lo, hi = self._branch_at(x, **kwargs)
        return hi - lo

    def center_at(self, x: float, **kwargs) -> float:
        lo, hi = self._branch_at(x, **kwargs)
        return (lo + hi) / 2.0

    def source_aperture_at(self, x: float, **kwargs) -> Optional[float]:
        """Widest launch window centred on the selected branch that stays
        clear of its neighbours — measured directly off the neighbouring
        y-intervals the scanline already produced, rather than computed from
        known port spacing as MMI does.

        Opt-in via `neighbor_aperture`, exactly like MMI's own
        `use_neighbor_aperture`, and off by default for the same reason:
        `half_domain_aperture`'s fixed half-domain rule already drives the
        devices that work today, and silently changing which aperture wins
        would change their results.
        """
        if not self.neighbor_aperture:
            return None
        ivs = self._intervals_at(x)
        if len(ivs) < 2:
            return None
        lo, hi = self._branch_at(x, **kwargs)
        c = (lo + hi) / 2.0
        gaps = [min(abs(c - a), abs(c - b)) for a, b in ivs if (a, b) != (lo, hi)]
        return 2.0 * min(gaps)

    # ------------------------------------------------------------------
    # Section extraction (HybridSolver hooks)
    # ------------------------------------------------------------------

    def _check_cut(self, x: float, kwargs: Dict) -> None:
        """Raise if `x` lands strictly inside a region whose cross-section
        varies — the measured stand-in for the band bookkeeping every
        parametric geometry does (MMI._check_cut against _taper_bounds(),
        Taper._check_cut against its own length, ...). An imported layout
        has no bands to consult, so the test has to sample.

        Legal iff the cross-section holds over a `cut_check_span`
        neighbourhood on **at least one** side of x. That "one side" rule is
        what reproduces the parametric geometries' own contract rather than
        being stricter than it: they permit a cut exactly ON a band endpoint
        (MMI allows x = ±mm_length/2, where a taper varies on one side and
        the uniform body sits on the other) and reject only a cut strictly
        inside the band, where both sides vary.

        Sampling is `assert_uniform_cross_section` — deliberately the shared
        helper `HybridSolver._validate_uniform_sweep_section` uses, not a
        third private sampler (see this plan's step 4 and
        mmi_1x4_tapered_body_validation.md step 1).

        Caveat, and a real limitation of having no parametric knowledge
        rather than a bug: a feature narrower than the sample spacing
        (`cut_check_span / (cut_check_samples - 1)`, 6.25 nm at the
        defaults), or a taper gentle enough that it moves the width by less
        than `cut_check_tol` across the whole span, passes unnoticed. Widen
        `cut_check_span` or tighten `cut_check_tol` for a layout with
        features at that scale.
        """
        d = self.cut_check_span
        for x0, x1 in ((x - d, x), (x, x + d)):
            try:
                assert_uniform_cross_section(
                    self, x0, x1,
                    what="", why="",
                    geometry_kwargs=kwargs,
                    samples=self.cut_check_samples,
                    tol=self.cut_check_tol,
                )
                return  # uniform on this side — cut is legal
            except ValueError:
                continue
        raise ValueError(
            f"Cut/stub at x={x:.4f} falls strictly inside a varying "
            f"cross-section: sampling {self.cut_check_samples} points over "
            f"±{d:.4g} µm finds width_at/center_at changing by more than "
            f"{self.cut_check_tol:.2g} µm on BOTH sides of it. A cut must "
            f"land where the cross-section holds on at least one side (a "
            f"uniform run, or the boundary of a tapered region). Measured, "
            f"not looked up — see GdsGeometry._check_cut for what that "
            f"cannot see."
        )

    def _clip(self, polys: List[List[Tuple[float, float]]], lo: float, hi: float
               ) -> List[List[Tuple[float, float]]]:
        """`polys` intersected with the vertical band [lo, hi], one gdstk
        boolean per polygon so abutting features stay separate polygons
        (the same shape MMI._section_polygons() returns) instead of being
        unioned into one."""
        import gdstk

        pad = 1.0 + (self._y_max - self._y_min)
        window = gdstk.rectangle((lo, self._y_min - pad), (hi, self._y_max + pad))
        out: List[List[Tuple[float, float]]] = []
        for poly in polys:
            if max(x for x, _ in poly) <= lo + self._X_TOL:
                continue
            if min(x for x, _ in poly) >= hi - self._X_TOL:
                continue
            for piece in gdstk.boolean([gdstk.Polygon(poly)], [window], "and",
                                       precision=self._boolean_precision):
                verts = [(float(px), float(py)) for px, py in piece.points]
                if len(verts) >= 3:
                    out.append(verts)
        return out

    def _section_polygons(self, x0: float, x1: float,
                           stub_left: float = 0.0, stub_right: float = 0.0,
                           **kwargs) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """Imported polygons clipped to [x0-stub_left, x1+stub_right], per
        core layer, in the device frame — plus, when `access_extension` is
        on, the straight access waveguides needed to reach past either end
        facet.

        Layers are kept apart (unlike section_prisms(), which has one
        medium for all of them) because section_gds()'s consumer,
        extrude_polygons(), applies a different material rule per layer.
        """
        lo, hi = x0 - stub_left, x1 + stub_right
        if hi - lo <= self._X_TOL:
            raise ValueError(
                f"section spans nothing: [{lo:.4g}, {hi:.4g}] µm "
                f"(x0={x0:.4g}, x1={x1:.4g}, stubs {stub_left:.4g}/{stub_right:.4g})."
            )
        self._check_cut(lo, kwargs)
        self._check_cut(hi, kwargs)

        needs_left = lo < self._x_min - self._X_TOL
        needs_right = hi > self._x_max + self._X_TOL
        if (needs_left or needs_right) and not self.access_extension:
            raise ValueError(
                f"section [{lo:.4g}, {hi:.4g}] µm reaches outside the imported "
                f"geometry's own span [{self._x_min:.4g}, {self._x_max:.4g}] µm, "
                f"but access_extension=False, so there is nothing drawn there. "
                f"A HybridSolver end section always asks for a stub past the "
                f"device's physical port. Either draw the access waveguides in "
                f"the GDS itself, or pass access_extension=True to synthesize "
                f"them from the end facet's own measured cross-section."
            )

        out: Dict[Tuple[int, int], List[List[Tuple[float, float]]]] = {}
        for layer, polys in self._polys.items():
            clipped = self._clip(polys, lo, hi)
            if needs_left:
                clipped.extend(
                    self._rects(lo, min(hi, self._x_min), self._facet_left_by_layer[layer]))
            if needs_right:
                clipped.extend(
                    self._rects(max(lo, self._x_max), hi, self._facet_right_by_layer[layer]))
            if clipped:
                out[layer] = clipped
        return out

    @staticmethod
    def _rects(xa: float, xb: float, intervals: List[Tuple[float, float]]
                ) -> List[List[Tuple[float, float]]]:
        if xb - xa <= 1e-9:
            return []
        return [[(xa, y0), (xb, y0), (xb, y1), (xa, y1)] for y0, y1 in intervals]

    def section_prisms(self, x0: float, x1: float, si_medium,
                        stub_left: float = 0.0, stub_right: float = 0.0,
                        wg_height: Optional[float] = None, **kwargs):
        """Every core layer is extruded in the SAME `si_medium` — the base
        API hands section_prisms() one medium, and this class calls all of
        `core_layers` core. A genuinely multi-material import (Si + SiN,
        say) is therefore an EME-path device: section_gds() keeps the layers
        separate for extrude_polygons()' per-layer rules."""
        verts = [poly for polys in
                 self._section_polygons(x0, x1, stub_left, stub_right, **kwargs).values()
                 for poly in polys]
        return polygons_to_prisms(verts, si_medium, wg_height=wg_height)

    def section_gds(self, x0: float, x1: float, oxide_width: float,
                     stub_left: float = 0.0, stub_right: float = 0.0, **kwargs
                     ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        lo, hi = x0 - stub_left, x1 + stub_right
        sections = self._section_polygons(x0, x1, stub_left, stub_right, **kwargs)
        out = {
            layer: [[(x - lo, y) for x, y in poly] for poly in polys]
            for layer, polys in sections.items()
        }
        out[_OX_LAYER] = [_oxide_box(0.0, hi - lo, oxide_width)]
        return out
