# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Georgia Tech Research Corporation
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple


def _oxide_box(x_start: float, total_x: float, oxide_width: float) -> List[Tuple[float, float]]:
    return [
        (x_start, -oxide_width / 2),
        (total_x,  -oxide_width / 2),
        (total_x,   oxide_width / 2),
        (x_start,  oxide_width / 2),
    ]


def assert_uniform_cross_section(
    geometry,
    x0: float,
    x1: float,
    *,
    what: str,
    why: str,
    geometry_kwargs: Optional[Dict] = None,
    samples: int = 25,
    tol: float = 1e-3,
) -> None:
    """Raise if `geometry`'s cross-section varies anywhere across [x0, x1].

    THE one implementation of "is this span uniform?", deliberately shared
    rather than re-derived per call site — the same reasoning
    `section_fdtd_layout()` is shared with `beef/plot/geometry.py`'s own
    marker drawing: this repo has already been bitten by duplicated geometry
    logic drifting apart.

    Exists because several code paths are only physical over a constant
    cross-section and produce a *plausible wrong answer*, not an error, when
    handed a varying one:

      - `EMESolver.solve_section(neff_override=...)` writes a diagonal
        ``T[i,i] = exp(2j*pi*neff/wl*L)`` outright, skipping the mode solve.
        That diagonal form follows from left/right ports sharing ONE
        cross-section's mode basis (single cell, `Ls=[length]`) — true for a
        uniform span, false for a tapered one, where the transfer matrix
        depends on the mode *fields* at every internal interface and is not
        diagonal at all.
      - `EMESolver.solve_section_length_sweep()` / `HybridSolver.
        solve_length_sweep()` extrapolate an analytic propagation phase
        across swept lengths, which means the same thing.

    Sampling `width_at`/`center_at` is the same test
    `HybridSolver._validate_uniform_sweep_section` has always used; both
    are needed, since a section can hold its width while its centre walks
    (an S-bend) or vice versa.

    Args:
        geometry: any `DeviceGeometry` — only `width_at`/`center_at` are used.
        x0, x1: device-frame span to check, inclusive of both endpoints.
        what: who is asking, for the message's leading clause (e.g.
            ``"sweep_section=1"``, ``"section 2/3 (EMESolver)"``).
        why: what breaks on a non-uniform span, for the trailing clause —
            state the requirement, not just the refusal.
        geometry_kwargs: forwarded to `width_at`/`center_at` (e.g.
            ``output_idx``), since those are per-branch on a multi-port device.
        samples: how many points across the span to test.
        tol: µm tolerance; a span varying by less than this counts as uniform.

    Raises:
        ValueError: naming the first offending x, so the caller can see
            *where* the span stops being uniform rather than only that it does.
    """
    kwargs = geometry_kwargs or {}
    n = max(int(samples), 2)
    # Plain arithmetic rather than np.linspace: this module is stdlib-only
    # today and there is no reason for a geometry predicate to change that.
    span = x1 - x0
    xs = [x0 + span * i / (n - 1) for i in range(n)]
    w0 = geometry.width_at(xs[0], **kwargs)
    c0 = geometry.center_at(xs[0], **kwargs)
    for x in xs[1:]:
        if (abs(geometry.width_at(x, **kwargs) - w0) > tol
                or abs(geometry.center_at(x, **kwargs) - c0) > tol):
            raise ValueError(
                f"{what}: cross-section is not uniform over "
                f"[{x0:.4g}, {x1:.4g}] (width_at/center_at varies at "
                f"x={x:.4g}); {why}."
            )


def polygons_to_prisms(vertex_lists, si_medium, wg_height: Optional[float] = None):
    """meep Prisms from raw (x, y) vertex lists — the one place the 2-D/3-D
    extrusion convention lives.

    wg_height None (2-D): infinite extrusion at z=0, for a planar (x, y)
    simulation. wg_height h (3-D): a real silicon core of thickness h
    centered on z=0 (vertices placed at z=-h/2, extruded +h along z), with
    the oxide above/below coming from the simulation's default_material.
    """
    import meep as mp

    h = mp.inf if wg_height is None else wg_height
    z0 = 0.0 if wg_height is None else -wg_height / 2
    return [
        mp.Prism(
            vertices=[mp.Vector3(x, y, z0) for x, y in verts],
            height=h,
            material=si_medium,
        )
        for verts in vertex_lists
    ]


def polygons_to_material_grid(
    vertex_lists,
    oxide_medium,
    si_medium,
    cell_size_xy: Tuple[float, float],
    grid_resolution: float,
    wg_height: Optional[float] = None,
):
    """meep MaterialGrid (rasterized) alternative to polygons_to_prisms()'s
    analytic Prisms, for isolating geometry-construction method as a
    variable independent of everything else the solver does.

    Exists because reproducing the modal multiplexer manuscript's own
    reference FDTD script (mux_sim_test.py, a standalone, non-beef
    script built from the actual archived reference-generation code)
    found that its meep-native, gdstk.inside()-rasterized MaterialGrid
    construction reproduces the paper's numbers to ~0.02-0.03 dB, while
    this repo's own analytic-Prism construction (polygons_to_prisms())
    leaves a 0.77/1.85 dB shortfall on the same device, same resolution,
    same eps_averaging, same dispersion model -- see
    modal_multiplexer_validation.md gap 9 item 6. This function lets
    that same device be rebuilt through beef's own polygons() with only
    the construction method swapped, to test whether the correlation
    holds on beef's own geometry (not just the reference script's) and,
    eventually, isolate which part of the rasterization pipeline matters.

    grid_resolution is the MaterialGrid's own internal rasterization
    density (pixels/µm) -- independent of, and per the reference script
    several times finer than, the FDTD Yee-grid `resolution` passed to
    `mp.Simulation`. meep's own subpixel averaging (`eps_averaging=True`)
    then interpolates this fine mask down onto the coarser Yee grid,
    exactly as it does for an analytic Prism boundary.
    """
    import gdstk
    import numpy as np
    import meep as mp

    Sx, Sy = cell_size_xy
    Nx = max(int(round(Sx * grid_resolution)), 1)
    Ny = max(int(round(Sy * grid_resolution)), 1)
    x = np.linspace(-Sx / 2, Sx / 2, Nx)
    y = np.linspace(-Sy / 2, Sy / 2, Ny)
    X, Y = np.meshgrid(x, y)
    query_points = np.vstack((X.flatten(), Y.flatten())).T

    gdstk_polygons = [gdstk.Polygon(verts) for verts in vertex_lists]
    mask = np.array(gdstk.inside(query_points, gdstk_polygons)).reshape((Ny, Nx)).T

    h = mp.inf if wg_height is None else wg_height
    design_grid = mp.MaterialGrid(
        grid_size=mp.Vector3(Nx, Ny, 1),
        medium1=oxide_medium,
        medium2=si_medium,
        weights=mask.astype(float).flatten(),
    )
    return [mp.Block(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, h), material=design_grid)]


class DeviceGeometry(ABC):
    """Abstract base for all BEEF device geometries.

    Subclasses define their physical shape via polygons(), which is consumed by
    to_meep_prisms() (FDTD) and to_gds_eme() (EME). The GDS-intermediate path
    for FDTD is deprecated; new geometries should implement polygons() only.
    """

    # ------------------------------------------------------------------
    # Polygon interface (preferred for new geometries)
    # ------------------------------------------------------------------

    def polygons(
        self, stub_length: float = 0.0
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        """Waveguide polygon definitions in device frame (device center at x=0).

        Args:
            stub_length: access waveguide stub length to include on each side (µm).
                         Pass 0 for the bare MMI body + tapers only.

        Returns:
            Mapping from GDS layer tuple → list of polygons.
            Each polygon is a list of (x_µm, y_µm) vertex tuples (CCW winding).

        Note:
            The device frame places the longitudinal center of the MMI body at
            x=0.  For a standard FDTD domain
            ``Sx = total_length + 2*stub_length + 2*dpml``, the meep origin also
            coincides with x=0, so polygon coordinates map directly to meep
            coordinates without any shift.
        """
        raise NotImplementedError(
            f"{type(self).__name__} has not implemented polygons(); "
            "override it or use to_gds_fdtd() for the legacy GDS path."
        )

    def to_meep_prisms(self, si_medium, stub_length: float = 0.0,
                        wg_height: Optional[float] = None):
        """Convert waveguide polygons to a list of meep Prism objects.

        No GDS files are written.  Polygons from the (2, 0) layer are extruded
        with the given Si medium; the SiO₂ background is handled by meep's
        ``default_material`` setting.

        Args:
            si_medium: ``mp.Medium`` for silicon at the simulation wavelength.
            stub_length: access waveguide stub length passed to polygons() (µm).
            wg_height: silicon core thickness (µm) for a 3-D simulation —
                the core is extruded z=-wg_height/2..+wg_height/2; None
                (default) keeps the legacy 2-D infinite extrusion. See
                polygons_to_prisms().

        Returns:
            List of ``mp.Prism`` objects ready for ``mp.Simulation(geometry=...)``.
        """
        return polygons_to_prisms(
            self.polygons(stub_length=stub_length).get((2, 0), []),
            si_medium, wg_height=wg_height,
        )

    def to_meep_material_grid(
        self,
        oxide_medium,
        si_medium,
        cell_size_xy: Tuple[float, float],
        grid_resolution: float,
        stub_length: float = 0.0,
        wg_height: Optional[float] = None,
    ):
        """Convert waveguide polygons to a rasterized meep MaterialGrid,
        the same device geometry to_meep_prisms() builds analytically.
        See polygons_to_material_grid() for why this exists.
        """
        return polygons_to_material_grid(
            self.polygons(stub_length=stub_length).get((2, 0), []),
            oxide_medium, si_medium, cell_size_xy, grid_resolution,
            wg_height=wg_height,
        )

    # ------------------------------------------------------------------
    # GDS interface (still required for meow / EME)
    # ------------------------------------------------------------------

    @abstractmethod
    def to_gds_eme(
        self,
        straight_length: float,
        oxide_width: float,
        x_start: float = 0.0,
        **kwargs,
    ):
        """Build the waveguide/oxide polygon layout for EME (meow) simulation.

        Returns a Dict[layer_tuple, List[polygon]] (same shape as polygons()),
        with an extra oxide-cladding layer, covering
        [x_start, x_start + total_length + 2*straight_length] in the layout frame.

        Args:
            straight_length: length of access waveguides on each side (µm).
            oxide_width: y-span of the SiO₂ cladding rectangle (µm).
            x_start: left-edge x-coordinate of the layout frame (µm).
        """
        ...

    def to_gds_fdtd(
        self,
        domain_sx: float,
        dpml: float,
        waveguide_length: float,
        side: str = "input",
        **kwargs,
    ):
        """Build a waveguide-polygon layout (Dict[layer_tuple, List[polygon]])
        positioned for a MEEP simulation domain.

        Deprecated: prefer to_meep_prisms() for new-style geometries. Every
        in-tree geometry now implements polygons(); this hook is retained
        only for out-of-tree geometries still on the legacy GDS path.
        """
        raise NotImplementedError(
            f"{type(self).__name__} has not implemented to_gds_fdtd(); "
            "use to_meep_prisms() for new-style geometries."
        )

    # ------------------------------------------------------------------
    # Arbitrary-section interface (required for HybridSolver)
    # ------------------------------------------------------------------
    #
    # A device that wants to support HybridSolver's generic N-way cascade
    # must implement these three hooks. They let any solver ask for "the
    # geometry between x0 and x1" without either side knowing what kind of
    # feature (waveguide / taper / MMI body / ...) actually lives there.
    # Device frame: same coordinates as polygons() (device center at x=0).

    def section_prisms(
        self, x0: float, x1: float, si_medium,
        stub_left: float = 0.0, stub_right: float = 0.0, **kwargs
    ):
        """Meep Prisms covering the device geometry in [x0, x1] (device frame),
        extended by stub_left / stub_right on either side.

        Stubs and interior cut points (x0, x1) must land in a region of
        constant cross-section (e.g. a straight waveguide or the MMI body) —
        cutting inside a region whose cross-section varies with x (e.g. a
        taper) is not supported and should raise ValueError.
        """
        raise NotImplementedError(
            f"{type(self).__name__} has not implemented section_prisms(); "
            "required for use with HybridSolver."
        )

    def section_gds(
        self, x0: float, x1: float, oxide_width: float,
        stub_left: float = 0.0, stub_right: float = 0.0, **kwargs
    ):
        """Waveguide/oxide polygon layout covering [x0, x1] (+ stubs), for EME sections.

        Same semantics as section_prisms(), but returns a Dict[layer_tuple,
        List[polygon]] (same shape as polygons()) with the layout's left edge
        at x=0 (i.e. shifted so x0 - stub_left -> 0).
        """
        raise NotImplementedError(
            f"{type(self).__name__} has not implemented section_gds(); "
            "required for use with HybridSolver."
        )

    def width_at(self, x: float, **kwargs) -> float:
        """Local waveguide cross-section width at device-frame position x (µm)."""
        raise NotImplementedError(
            f"{type(self).__name__} has not implemented width_at(); "
            "required for use with HybridSolver."
        )

    def center_at(self, x: float, **kwargs) -> float:
        """Local waveguide y-center at device-frame position x (µm). Default 0."""
        return 0.0

    def source_aperture_at(self, x: float, **kwargs) -> Optional[float]:
        """Widest single-mode eigenmode-launch window (full y-span, µm),
        centered on center_at(x), that stays clear of neighboring structure
        at device-frame position x — or None for no geometry-imposed
        restriction (the default).

        FDTDSolver.solve_section() consults this whenever it excites exactly
        one mode: a non-None value shrinks the eigenmode source and its
        co-located normalization monitor to this window (further clamped to
        the solver's own non-PML domain) so the MPB mode solve neither sees
        a nearby unexcited waveguide nor reaches into the PML — while the
        drawn geometry stays complete. Override in geometries whose ports
        have close neighbors that half_domain_aperture's fixed half-domain
        rule can't isolate (e.g. an MMI with several closely-spaced output
        branches, one of which sits exactly on the y=0 axis — see MMI's
        use_neighbor_aperture)."""
        return None

    # ------------------------------------------------------------------
    # Shared properties / helpers
    # ------------------------------------------------------------------

    @property
    @abstractmethod
    def total_length(self) -> float:
        """Device length excluding access waveguides (µm)."""
        ...

    @abstractmethod
    def eme_lengths(self, straight_length: float) -> List[float]:
        """Region lengths for EME cell construction, left-to-right (µm).

        The list length must match the num_modes / num_cells lists passed to
        EMESolver.  For a standard MMI:
            [straight_length, taper_length, mm_length, taper_length, straight_length]
        """
        ...

    @property
    def bbox(self) -> Tuple[float, float]:
        """(width, height) bounding box of the device body only (µm)."""
        return self.total_length, 0.0

    @property
    def frame_offset(self) -> float:
        """Device-frame x-offset from this geometry's bbox midpoint to its
        own local-frame origin (x=0 in polygons()/width_at()/etc.).

        0.0 (the default) means the two coincide, i.e. this geometry's local
        frame spans exactly [-total_length/2, +total_length/2] — true for
        every geometry in this package except MMI. CompoundGeometry reads
        this to correctly place a segment placed at device-frame
        [seg_lo, seg_hi] whose own frame isn't centered on that bbox: the
        segment's local-frame origin sits at device x =
        (seg_lo + seg_hi) / 2 + frame_offset, not at the bbox midpoint
        itself.
        """
        return 0.0
