"""BEEF — Bridged Engine for EME and FDTD.

Core public API:

    Geometry:
        MMI, MMIPort               — N-input/M-output MMI geometry (1x2, 2x2, 1x4, ...),
                                      with optional per-port taper size/placement
        BendGeometry               — 90° bend for non-collinear FDTD I/O
        DirectionalCoupler         — 2-waveguide directional coupler w/ S-bend I/O
        Taper, SBend               — standalone taper / S-bend building-block geometries
        CompoundGeometry, WaveguideRun — concatenate several geometries end-to-end
        GdsGeometry, load_gds_polygons — import polygons from an on-disk .gds file
        arc_polygon                — quarter-annulus wedge, for constant-radius bends
        RingResonatorCouplerSection, RingResonatorBendSection — ring resonator geometry
                                      primitives (not yet wired to HybridSolver — see
                                      validation_markdown/ring_resonator_validation.md)

    Materials:
        MaterialStack              — Si/SiO2 layer definitions
        SOI_160, SOI_220           — pre-built stacks for common platforms

    Solvers:
        EMESolver                  — meow-based eigenmode expansion
        FDTDSolver                 — meep-based FDTD
        PortSpec                   — explicit source/monitor orientation override for
                                      FDTDSolver.solve_section() (x- or y-normal ports)
        HybridSolver               — FDTD/EME cascaded (the core BEEF method)

    Results:
        SMatrix                    — unified S-parameter result type

    Network:
        SMatrixBlock                — dense multi-port S-parameter block
        connect, close_loop        — Redheffer star product / self-loop elimination
        Circuit                    — stateful multi-block network builder
        from_transfer_matrix, from_smatrix, to_smatrix — SMatrixBlock <-> SMatrix/HybridSolver conversions

    Mode identity (device- and resolution-agnostic — see beef.mode_identity):
        field_overlap_matrix, best_correspondence — mode-shape overlap + optimal
                                      assignment, replacing nearest-neighbor n_eff
                                      matching for checking whether HybridSolver's
                                      index-i-to-index-i contraction is valid
        mesh_cell_centers, stack_meow_mode_fields — grid/field-array helpers tying
                                      FDTDSolver.section_mode_fields() (meep) and
                                      EMESolver.compute_modes() (meow) together

    Cut-position invariance (device-agnostic — see beef.cut_position_invariance):
        cut_position_invariance    — drift/slope of a solve()'d SMatrix across a
                                      sequence of FDTD/EME cut positions; the
                                      authoritative correspondence check
                                      (fdtd_vs_EME_mode_validation.md layer 3)

    Utilities:
        extrude_polygons           — {layer: [polygons]} dict -> meow Structure3D list
        atomic_replace_write       — write-to-tmp-then-os.replace, for files a
                                     concurrent SLURM array task may be reading

    Sweep:
        ParallelSweep              — parallel parameter grid sweep

    Plot:
        plot_eme_fields            — EME field visualization
        plot_eme_geometry          — EME device + region geometry plot
        plot_fdtd_geometry         — FDTD device + source/monitor/PML geometry plot
        plot_hybrid_geometry       — Hybrid device + section geometry plot (+ sub-region plots)
        plot_performance           — transmission/insertion loss plots
        plot_splitting_ratio       — splitting ratio plot
"""

__version__ = "0.1.0"

# Under `srun -n 48` every rank runs the whole solve, so every INFO beef emits
# came out 48 times — those messages describe one shared cache, not per-rank
# state. Each module gets its logger from beef._log.get_logger(), which keeps
# WARNING and above only off the master rank (rank-tagged). Nothing is done to
# the root logger. It has to be per-module: a logger's filters are consulted
# only for records logged on that logger, never for records propagating up from
# a child. See beef/_log.py.

from ._utils import atomic_replace_write, extrude_polygons
from .cut_position_invariance import cut_position_invariance
from .geometry import (
    MMI, MMIPort, BendGeometry, DirectionalCoupler, Taper, SBend,
    CompoundGeometry, WaveguideRun, GdsGeometry, load_gds_polygons,
    arc_polygon, RingResonatorCouplerSection, RingResonatorBendSection,
)
from .materials import MaterialStack, SOI_160, SOI_220
from .mode_identity import (
    field_overlap_matrix, best_correspondence, mesh_cell_centers, stack_meow_mode_fields,
    overlap_diagnostics, self_overlap_diagnostics,
)
from .network import (
    SMatrixBlock, connect, close_loop, Circuit,
    from_transfer_matrix, from_smatrix, to_smatrix,
)
from .results import SMatrix
from .solvers import EMESolver, FDTDSolver, PortSpec, HybridSolver
from .sweep import ParallelSweep
from .plot import (
    plot_eme_fields, plot_eme_geometry, plot_fdtd_geometry, plot_hybrid_geometry,
    plot_performance, plot_splitting_ratio,
)

__all__ = [
    # geometry
    "MMI", "MMIPort",
    "BendGeometry",
    "DirectionalCoupler",
    "Taper", "SBend",
    "CompoundGeometry", "WaveguideRun",
    "GdsGeometry", "load_gds_polygons",
    "arc_polygon", "RingResonatorCouplerSection", "RingResonatorBendSection",
    # materials
    "MaterialStack", "SOI_160", "SOI_220",
    # network
    "SMatrixBlock", "connect", "close_loop", "Circuit",
    "from_transfer_matrix", "from_smatrix", "to_smatrix",
    # mode identity
    "field_overlap_matrix", "best_correspondence", "mesh_cell_centers", "stack_meow_mode_fields",
    "overlap_diagnostics", "self_overlap_diagnostics",
    # cut-position invariance
    "cut_position_invariance",
    # results
    "SMatrix",
    # solvers
    "EMESolver", "FDTDSolver", "PortSpec", "HybridSolver",
    # utilities
    "extrude_polygons", "atomic_replace_write",
    # sweep
    "ParallelSweep",
    # plot
    "plot_eme_fields", "plot_eme_geometry", "plot_fdtd_geometry", "plot_hybrid_geometry",
    "plot_performance", "plot_splitting_ratio",
]
