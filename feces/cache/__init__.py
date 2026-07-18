"""Persistent save-folder cache for FDTDSolver/EMESolver section solves.

build_run_config() (config.py) fingerprints a solve_section() call (solver
params, call args, and the exact section geometry via section_gds()) into a
JSON-safe dict; lookup()/save_run() (store.py) hash that dict to a save-folder
subdirectory holding config.json + result.pkl + geometry.png. Wired into
FDTDSolver.solve_section(), EMESolver.solve_section()/
solve_section_length_sweep(), and threaded through HybridSolver via a new
save_dir parameter — empty/None disables the feature everywhere, matching the
HybridSolver.solve(cache=None) / ParallelSweep(output_dir=None) convention.
"""
from .config import build_run_config
from .store import lookup, load_result, save_run

__all__ = ["build_run_config", "lookup", "load_result", "save_run"]
