import numpy as np
from dataclasses import dataclass, field, asdict
from enum import Enum, auto
import json

@dataclass(eq=False, kw_only=True)
class GeoParamsMMI:
    waveguide_height: float = 0.16    # Height of waveguide layer
    waveguide_width: float = 0.35     # Width of single-mode waveguide
    mm_length: float = 4.0            # Length of multimode section
    mm_width: float = 2.0             # Width of multimode section
    taper_length: float = 1.0         # Length of input and output tapers
    taper_width: float = 0.8          # Width of taper at MMI interface
    waveguide_offset: float = 1.0     # Distance between centers of output waveguides
    gds_file: str = ""                # Path to pre-existing GDS file; used when sim_type="FDTD_GDS"

@dataclass(eq=False, kw_only=True)
class SimParams:
    oxide_height: float = 1.5         # Simulation height
    oxide_width: float = 5            # Simulation width
    res: float = 40                   # Simulation resolution in pixels per micron
    wavelength_min: float = 1.26      # Min wavelength of simulation in microns
    wavelength_max: float = 1.36      # Max wavelength of simulation in microns
    wavelength_pts: int = 10          # Number of wavelength points
    num_modes_MMI: int = 15           # Number of modes calculated in the multimode section
    dpml: float = 0.5                 # Thickness of PML layer in microns
    mm_sim_length: float = 1.0        # Length of the MMI section to simulate
    waveguide_length: float = 0.5     # Length of waveguide between PML and taper
    source_position: float = 0.0      # x_position from the PML where the source is placed
    sm_monitor_position: float = 0.25 # x_position from the PML where the single-mode waveguide mode monitor is
    mm_monitor_position: float = 0.25 # x_position from the PML where the multi-mode waveguide mode monitor is
    sim_type: str = "FDTD"            # Simulation type: "FDTD", "EME", or "HYBRID"
    num_cells_taper: int = 15         # Number of cells the taper regions are broken up into
    num_modes_taper: int = 7          # Number of modes computed in the taper regions

@dataclass(eq=False, kw_only=True)
class HyperParamsMMI:
    geo: GeoParamsMMI = field(default_factory=GeoParamsMMI)
    sim: SimParams = field(default_factory=SimParams)

    outputfolder: str = "./out/test"