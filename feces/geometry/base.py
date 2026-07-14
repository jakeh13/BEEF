from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, List, Tuple


class DeviceGeometry(ABC):
    """Abstract base for all FECES device geometries.

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

    def to_meep_prisms(self, si_medium, stub_length: float = 0.0):
        """Convert waveguide polygons to a list of meep Prism objects.

        No GDS files are written.  Polygons from the (2, 0) layer are extruded
        with the given Si medium; the SiO₂ background is handled by meep's
        ``default_material`` setting.

        Args:
            si_medium: ``mp.Medium`` for silicon at the simulation wavelength.
            stub_length: access waveguide stub length passed to polygons() (µm).

        Returns:
            List of ``mp.Prism`` objects ready for ``mp.Simulation(geometry=...)``.
        """
        import meep as mp

        return [
            mp.Prism(
                vertices=[mp.Vector3(x, y, 0) for x, y in verts],
                height=mp.inf,
                material=si_medium,
            )
            for verts in self.polygons(stub_length=stub_length).get((2, 0), [])
        ]

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

        Deprecated: prefer to_meep_prisms() for new-style geometries.
        Retained for MMIMux and other complex geometries that have not yet
        migrated to the polygon-first interface.
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
