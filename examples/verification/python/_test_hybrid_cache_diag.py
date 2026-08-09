"""TEMPORARY -- diagnoses why the full 4-branch/281-length hybrid sweep
(job 11752875) fell through to a fresh FDTD solve on section 0, when the
earlier 3-length/1-branch bounded test (job 11752541) hit cache in 341s with
an IDENTICAL section-0/section-2 configuration. cache_only=True throughout
so this can never accidentally trigger a real multi-hour solve -- any miss
raises CacheMissError immediately instead.
"""
import numpy as np

from beef import MMI, EMESolver, FDTDSolver, HybridSolver, SIN_300
from beef.cache.store import CacheMissError

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
save_dir = "/storage/scratch1/8/jhiesener3/duty/example_sim_cache"

print("=" * 60)
print("TEST A: 3-length grid (matches the SUCCEEDED bounded test), branch 0, cache_only=True")
try:
    coeffs = hybrid.solve_length_sweep(
        mmi, wavelengths, np.array([190.0, 206.0, 220.0]), sweep_section=1,
        geometry_kwargs={2: {"output_idx": 0}},
        cache=dict(), save_dir=save_dir, cache_only=True,
    )
    print("HIT -- shape:", coeffs.shape)
except CacheMissError as e:
    print("MISS:", str(e)[:300])

print("=" * 60)
print("TEST B: FULL 281-point grid, branch 0, cache_only=True")
try:
    coeffs = hybrid.solve_length_sweep(
        mmi, wavelengths, np.linspace(20.0, 300.0, 281), sweep_section=1,
        geometry_kwargs={2: {"output_idx": 0}},
        cache=dict(), save_dir=save_dir, cache_only=True,
    )
    print("HIT -- shape:", coeffs.shape)
except CacheMissError as e:
    print("MISS:", str(e)[:300])

print("=" * 60)
print("TEST C: 3-length grid, branch 1 (output_idx=1), cache_only=True")
try:
    coeffs = hybrid.solve_length_sweep(
        mmi, wavelengths, np.array([190.0, 206.0, 220.0]), sweep_section=1,
        geometry_kwargs={2: {"output_idx": 1}},
        cache=dict(), save_dir=save_dir, cache_only=True,
    )
    print("HIT -- shape:", coeffs.shape)
except CacheMissError as e:
    print("MISS:", str(e)[:300])
