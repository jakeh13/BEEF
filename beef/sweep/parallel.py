"""ParallelSweep: generic parallel parameter sweep over any BEEF solver.

Replaces the device-specific sweep scripts (mmi_2x2_sweep_parallel_efficient.py, etc.)
with a reusable parallel executor that works with any Solver + DeviceGeometry pair.
"""
from __future__ import annotations

import gc
import itertools
import logging
import os
import pickle
from multiprocessing import Pool
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from ..geometry.base import DeviceGeometry
from ..results.smatrix import SMatrix
from .._utils import make_output_folders

logger = logging.getLogger(__name__)


def _worker_fn(args):
    """Top-level (picklable) worker function for multiprocessing.Pool."""
    solver_fn, geometry_fn, params, wavelengths, out_dir, save = args
    try:
        geom = geometry_fn(**params)
        solver = solver_fn()
        smat = solver.solve(geom, wavelengths)
        if save and out_dir:
            tag = "_".join(f"{k}={v:.4g}" for k, v in sorted(params.items()))
            pt_dir = os.path.join(out_dir, tag)
            make_output_folders(pt_dir)
            smat.save(os.path.join(pt_dir, "data", "smatrix.pkl"))
        return (params, smat)
    except Exception as exc:
        logger.error("Worker failed for params=%s: %s", params, exc)
        raise


class ParallelSweep:
    """Run a parameter grid sweep in parallel using multiprocessing.Pool.

    Args:
        solver_fn: zero-arg callable returning a configured Solver instance.
                   Must be a factory (not the instance itself) so each worker
                   gets its own Solver without sharing state across processes.
        geometry_fn: callable(**params) -> DeviceGeometry for each grid point.
        param_grid: dict of {param_name: [values]}.  All combinations are swept.
        wavelengths: [M] wavelength array (µm).
        n_workers: number of parallel worker processes.
        output_dir: base directory for per-point output folders. None = no disk writes.
        save_per_point: if True and output_dir is set, write SMatrix pickle per point.
    """

    def __init__(
        self,
        solver_fn: Callable[[], Any],
        geometry_fn: Callable[..., DeviceGeometry],
        param_grid: Dict[str, List[Any]],
        wavelengths: np.ndarray,
        n_workers: int = 4,
        output_dir: Optional[str] = None,
        save_per_point: bool = True,
    ):
        self.solver_fn = solver_fn
        self.geometry_fn = geometry_fn
        self.param_grid = param_grid
        self.wavelengths = np.asarray(wavelengths)
        self.n_workers = n_workers
        self.output_dir = output_dir
        self.save_per_point = save_per_point

    def _param_combinations(self) -> List[Dict[str, Any]]:
        keys = list(self.param_grid.keys())
        values = list(self.param_grid.values())
        return [dict(zip(keys, combo)) for combo in itertools.product(*values)]

    def run(self) -> Dict[Tuple, SMatrix]:
        """Execute the full sweep and return {param_tuple: SMatrix} dict.

        param_tuple is a sorted tuple of (key, value) pairs for the grid point.
        """
        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)

        combos = self._param_combinations()
        tasks = [
            (self.solver_fn, self.geometry_fn, params,
             self.wavelengths, self.output_dir, self.save_per_point)
            for params in combos
        ]

        results: Dict[Tuple, SMatrix] = {}
        with Pool(processes=self.n_workers) as pool:
            for params, smat in pool.imap_unordered(_worker_fn, tasks):
                key = tuple(sorted(params.items()))
                results[key] = smat
                logger.info("Done: %s", params)
                gc.collect()

        if self.output_dir:
            summary_path = os.path.join(self.output_dir, "sweep_results.pkl")
            with open(summary_path, "wb") as f:
                pickle.dump(results, f)
            logger.info("Sweep complete. Results saved to %s", summary_path)

        return results
