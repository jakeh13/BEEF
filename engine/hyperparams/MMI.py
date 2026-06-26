from dataclasses import dataclass, field
from typing import List


def waveguide_params_from_gds(gds_file: str, edge: str = "right", port_zone: float = 0.05) -> tuple[float, float]:
    """Measure (waveguide_offset, waveguide_width) from the two single-mode
    ports at one edge of a GDS device (e.g. a 2x2 MMI splitter's output side).
    """
    import pya

    layout = pya.Layout()
    layout.read(gds_file)
    dbu = layout.dbu
    cell = layout.top_cell()

    region = pya.Region()
    for li in layout.layer_indices():
        region += pya.Region(cell.shapes(li))

    bb = region.bbox()
    x = bb.right if edge == "right" else bb.left
    half_zone = int(round(port_zone / dbu))
    pad = bb.height()
    slab = pya.Region(pya.Box(x - half_zone, bb.bottom - pad, x + half_zone, bb.top + pad))

    intervals = []
    for poly in (region & slab).merged().each():
        ys = [p.y for p in poly.each_point_hull()]
        intervals.append((min(ys) * dbu, max(ys) * dbu))
    intervals.sort()

    if len(intervals) != 2:
        raise ValueError(f"Expected 2 ports at the {edge} edge of {gds_file}, found {len(intervals)}")

    (y1_min, y1_max), (y2_min, y2_max) = intervals
    waveguide_width = ((y1_max - y1_min) + (y2_max - y2_min)) / 2
    waveguide_offset = (y2_min + y2_max) / 2 - (y1_min + y1_max) / 2
    return waveguide_offset, waveguide_width


@dataclass(eq=False, kw_only=True)
class GeoParamsMMI:
    waveguide_height: float = 0.16    # Height of waveguide layer
    waveguide_width: float = 0.35     # Width of single-mode waveguide
    mm_length: float = 4.0            # Length of multimode section
    mm_width: float = 2.0             # Width of multimode section
    taper_length: float = 1.0         # Length of input and output tapers
    taper_width: float = 0.8          # Width of taper at MMI interface
    waveguide_offset: float = 1.0     # Distance between centers of adjacent waveguides
    gds_file: str = ""                # Path to pre-existing GDS file; used when sim_type="FDTD_GDS"
    dilation: float = 0.0             # GDS body dilation in µm (positive=expand, negative=erode); 0 disables
    bend_radius: float = 5.0          # Bend radius in microns (FDTD_GDS_BEND only)
    bend_type: str = "euler"          # "euler" or "circular" (FDTD_GDS_BEND only)
    bend_locations: List[str] = field(default_factory=list)  # ["input"], ["output"], or ["input","output"]
    bend_wg_length: float = 1.0       # Straight waveguide length after/before bend for monitor placement
    bend_stub_length: float = 0.5     # Length of waveguide between MMI device edge and bend start

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
    sim_type: str = "FDTD"            # Simulation type: "FDTD", "EME", "HYBRID", or "FDTD_GDS"
    num_cells_taper: int = 15         # Number of cells the taper regions are broken up into
    num_modes_taper: int = 7          # Number of modes computed in the taper regions
    input_modes: List[int] = field(default_factory=lambda: [1])  # Waveguide modes to excite at the input (1-based)
    output_modes: List[int] = field(default_factory=lambda: [1]) # Waveguide modes to measure at each output arm
    input_ports: List[str] = field(default_factory=lambda: ["top"]) # Input ports to excite: "top", "bot", or both

@dataclass(eq=False, kw_only=True)
class HyperParamsMMI:
    geo: GeoParamsMMI = field(default_factory=GeoParamsMMI)
    sim: SimParams = field(default_factory=SimParams)
    outputfolder: str = "./out/test"
