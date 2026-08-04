"""Persistent save-folder cache for FDTDSolver/EMESolver solves.

build_run_config() (config.py) fingerprints a solve_section()/
solve_section_length_sweep() call (solver params, call args, and the exact
section geometry via section_gds()) into a JSON-safe dict;
build_whole_device_config() does the same for the plain whole-device solve()
methods, using geometry.polygons()/to_gds_eme() instead since those calls
have no x0/x1 section to fingerprint via section_gds(). lookup()/save_run()
(store.py) hash either dict to a save-folder subdirectory holding
config.json + result.pkl + geometry.png. Wired into FDTDSolver.solve()/
solve_section(), EMESolver.solve()/solve_section()/
solve_section_length_sweep(), and threaded through HybridSolver via a new
save_dir parameter — empty/None disables the feature everywhere, matching the
HybridSolver.solve(cache=None) / ParallelSweep(output_dir=None) convention.

cache_only=True (same call sites) demands a cache hit instead of falling back
to simulating on a miss — raises CacheMissError. Notebooks (not scripts) set
this so re-running for documentation/rendering purposes can never silently
kick off a real (possibly multi-hour) simulation.

lookup_superset()/update_family_pointer() (store.py) add a second, orthogonal
lookup path: a request for fewer modes than a prior call already computed
(same solver/geometry/monitor "family", mode counts aside) is served by
slicing that larger result's mode axes instead of resimulating. See
lookup_superset()'s docstring for the physical preconditions this relies on
(same excited/measured side, pinned monitor geometry) — solve_section()/
solve_section_length_sweep() are responsible for only attempting it when
those hold.
"""
from .config import build_run_config, build_whole_device_config
from .store import (
    CacheMissError, is_master, lookup, lookup_superset, load_result, save_run,
    update_family_pointer,
)

__all__ = [
    "build_run_config", "build_whole_device_config", "CacheMissError", "is_master",
    "lookup", "lookup_superset", "load_result", "save_run", "update_family_pointer",
]
