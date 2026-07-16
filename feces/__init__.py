"""FECES — FDTD/EME Cascaded Electromagnetic Solver.

Core public API:

    Geometry:
        MMI, MMIPort               — N-input/M-output MMI geometry (1x2, 2x2, 1x4, ...),
                                      with optional per-port taper size/placement
        BendGeometry               — 90° bend for non-collinear FDTD I/O
        DirectionalCoupler         — 2-waveguide directional coupler w/ S-bend I/O
        Taper, SBend               — standalone taper / S-bend building-block geometries
        CompoundGeometry, WaveguideRun — concatenate several geometries end-to-end

    Materials:
        MaterialStack              — Si/SiO2 layer definitions
        SOI_160, SOI_220           — pre-built stacks for common platforms

    Solvers:
        EMESolver                  — meow-based eigenmode expansion
        FDTDSolver                 — meep-based FDTD
        HybridSolver               — FDTD/EME cascaded (the core FECES method)

    Results:
        SMatrix                    — unified S-parameter result type

    Utilities:
        extrude_polygons           — {layer: [polygons]} dict -> meow Structure3D list

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

from ._utils import extrude_polygons
from .geometry import (
    MMI, MMIPort, BendGeometry, DirectionalCoupler, Taper, SBend,
    CompoundGeometry, WaveguideRun,
)
from .materials import MaterialStack, SOI_160, SOI_220
from .results import SMatrix
from .solvers import EMESolver, FDTDSolver, HybridSolver
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
    # materials
    "MaterialStack", "SOI_160", "SOI_220",
    # results
    "SMatrix",
    # solvers
    "EMESolver", "FDTDSolver", "HybridSolver",
    # utilities
    "extrude_polygons",
    # sweep
    "ParallelSweep",
    # plot
    "plot_eme_fields", "plot_eme_geometry", "plot_fdtd_geometry", "plot_hybrid_geometry",
    "plot_performance", "plot_splitting_ratio",
]
