"""FECES — FDTD/EME Combined Electromagnetic Solver.

Core public API:

    Geometry:
        MMI1x2, MMI2x2, MMI1x4   — MMI device geometry classes
        MMIMux, MuxGeoParams       — cascaded MUX geometry
        BendGeometry               — 90° bend for non-collinear FDTD I/O

    Materials:
        MaterialStack              — Si/SiO2 layer definitions
        SOI_160, SOI_220           — pre-built stacks for common platforms

    Solvers:
        EMESolver                  — meow-based eigenmode expansion
        FDTDSolver                 — meep-based FDTD
        HybridSolver               — FDTD/EME combined (the core FECES method)

    Results:
        SMatrix                    — unified S-parameter result type

    Sweep:
        ParallelSweep              — parallel parameter grid sweep

    Plot:
        plot_eme_fields            — EME field visualization
        plot_performance           — transmission/insertion loss plots
        plot_splitting_ratio       — splitting ratio plot
"""

__version__ = "0.1.0"

from .geometry import MMI1x2, MMI2x2, MMI1x4, MMIMux, MuxGeoParams, BendGeometry
from .materials import MaterialStack, SOI_160, SOI_220
from .results import SMatrix
from .solvers import EMESolver, FDTDSolver, HybridSolver
from .sweep import ParallelSweep
from .plot import plot_eme_fields, plot_performance, plot_splitting_ratio

__all__ = [
    # geometry
    "MMI1x2", "MMI2x2", "MMI1x4",
    "MMIMux", "MuxGeoParams",
    "BendGeometry",
    # materials
    "MaterialStack", "SOI_160", "SOI_220",
    # results
    "SMatrix",
    # solvers
    "EMESolver", "FDTDSolver", "HybridSolver",
    # sweep
    "ParallelSweep",
    # plot
    "plot_eme_fields", "plot_performance", "plot_splitting_ratio",
]
