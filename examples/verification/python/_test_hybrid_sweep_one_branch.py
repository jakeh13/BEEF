"""TEMPORARY -- one branch, full 281-length grid, real (cache_only=False)
solve. Isolates whether the full grid alone (not the 4-branch loop) is what
triggered job 11752875's fresh FDTD solve.
"""
import time

import numpy as np

from beef import MMI, EMESolver, FDTDSolver, HybridSolver, SIN_300

RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 2.4, 12.0
MM_LENGTH = 206.0
N_MODES = 12

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30, use_neighbor_aperture=True,
)

_fdtd_common = dict(
    res=RESOLUTION, sim_width=MM_WIDTH + 8.0, sim_height=SIN_300.sz, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SIN_300, decay_tol=1e-11,
)
input_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=MM_WIDTH + 8.0, oxide_height=SIN_300.sz,
    res=RESOLUTION, material_stack=SIN_300, mode_solver="tidy3d",
)
output_solver = FDTDSolver(**_fdtd_common)

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
body_stub = 0.25
boundaries = [-half_total, -half_body + body_stub, +half_body - body_stub, +half_total]
interface_modes = [N_MODES, N_MODES]
hybrid = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries, interface_modes=interface_modes,
)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
mm_lengths = np.linspace(20.0, 300.0, 281)
save_dir = "/storage/scratch1/8/jhiesener3/duty/example_sim_cache"

print(f"Real run (cache_only=False), FULL {len(mm_lengths)}-length grid, branch 0...", flush=True)
t0 = time.time()
coeffs = hybrid.solve_length_sweep(
    mmi, wavelengths, mm_lengths, sweep_section=1,
    geometry_kwargs={2: {"output_idx": 0}},
    cache=dict(), save_dir=save_dir, cache_only=False,
)
dt = time.time() - t0
print(f"DONE in {dt:.1f} s. shape: {coeffs.shape}", flush=True)
T = np.abs(coeffs[:, :, 0, 0]) ** 2
print("T band-mean at first/mid/last length:", T.mean(axis=1)[[0, 140, -1]], flush=True)
