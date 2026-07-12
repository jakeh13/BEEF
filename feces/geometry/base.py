from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Tuple


class DeviceGeometry(ABC):
    """Abstract base for all FECES device geometries.

    Subclasses provide canonical polygon builders consumed by both the
    EME solver (via meow/GDS extrusion) and the FDTD solver (via MEEP GDS prisms).
    This eliminates the per-solver geometry duplication in the legacy codebase.
    """

    @abstractmethod
    def to_gds_eme(self, straight_length: float, oxide_width: float,
                   x_start: float = 0.0, **kwargs):
        """Build a gdsfactory Component for EME simulation.

        The component spans [x_start, x_start + total_length + 2*straight_length]
        in the layout frame (left edge at x_start).

        Args:
            straight_length: length of access waveguides on each side (µm).
            oxide_width: y-span of the SiO2 cladding rectangle (µm).
            x_start: left edge x-coordinate (µm).
        """
        ...

    @abstractmethod
    def to_gds_fdtd(self, domain_sx: float, dpml: float, waveguide_length: float,
                    side: str = "input", **kwargs):
        """Build a gdsfactory Component positioned for a MEEP simulation domain.

        The MEEP domain is centered at the origin: x ∈ [-domain_sx/2, +domain_sx/2].
        PML occupies the outermost `dpml` on each side; `waveguide_length` is the
        straight access waveguide between the PML and the device taper.

        Args:
            domain_sx: total x-span of the MEEP simulation domain (µm).
            dpml: PML thickness (µm).
            waveguide_length: access waveguide length inside the domain (µm).
            side: "input" (left half of device) or "output" (right half).
        """
        ...

    @property
    @abstractmethod
    def total_length(self) -> float:
        """Device length without access waveguides (µm)."""
        ...

    @abstractmethod
    def eme_lengths(self, straight_length: float) -> List[float]:
        """Region lengths for EME cell construction, left to right (µm).

        The length of this list must match the num_modes / num_cells lists
        passed to EMESolver.  For a standard MMI this is:
            [straight_length, taper_length, mm_length, taper_length, straight_length]
        """
        ...

    @property
    def bbox(self) -> Tuple[float, float]:
        """(width, height) bounding box of device body only (µm)."""
        return self.total_length, 0.0
