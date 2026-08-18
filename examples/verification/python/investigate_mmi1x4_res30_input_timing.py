"""TIMING PROBE, not a full run -- times mmi_1x4_hybrid.py's own input-section
FDTD solve at RESOLUTION=30 for real, before committing the full 4-branch
prewarm blind. This IS the real production call (n_modes_in=1, n_modes_out=
N_MODES=12, same geometry/decay_tol as production) with save_dir set, so its
result is genuinely cached and reusable by the real run afterward -- not a
throwaway probe, unlike the calibration-band timing script run alongside it.

Purpose: directly test whether the (30/40)^4 ~= 3.16x FDTD cost-scaling
projection (grid cells x timesteps) actually holds for this device's own
weakly-confined nitride platform, by comparing this measured res=30 decay
time against the KNOWN res=40 baseline (~11-12 hrs historically for this
input/output section class, and this session's own job 11939350 measured
12,706 timesteps / ~7.3 hrs for a similarly-sized calibration-band domain at
res=40).

Run:
    python .../investigate_mmi1x4_res30_input_timing.py
"""
import os
import time
import numpy as np

from beef import MMI, FDTDSolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

RESOLUTION = 30.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 3.1, 12.0
MM_LENGTH = 39.90
N_MODES = 12
STRAIGHT_LENGTH = 0.5
body_stub = 0.25
DECAY_TOL = 1e-11

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30, use_neighbor_aperture=True,
)

_fdtd_common = dict(
    res=RESOLUTION, sim_width=MM_WIDTH + 8.0, sim_height=SIN_300.sz,
    dpml=0.5, source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SIN_300, decay_tol=DECAY_TOL,
)
input_solver = FDTDSolver(**_fdtd_common)

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
boundaries = [-half_total, -half_body + body_stub, +half_body - body_stub, +half_total]
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

print(f"INPUT SECTION TIMING PROBE, res={RESOLUTION}, decay_tol={DECAY_TOL}, "
      f"N_MODES={N_MODES}, N_WL={N_WL}", flush=True)
t0 = time.time()
M_in1 = input_solver.solve_section(
    mmi, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=input_solver.waveguide_length, stub_right=0.0,
    input_idx=0, save_dir=save_dir,
)
elapsed = time.time() - t0
print(f"\nINPUT SECTION DONE: {elapsed:.1f} s ({elapsed/3600:.3f} hrs)", flush=True)
print(f"Comparison: res=40 historical baseline for this section class was "
      f"~11-12 hrs (input/output taper+stub, pre-calibration era) and this "
      f"session's own job 11939350 measured 7.3 hrs for a similarly-sized "
      f"CALIBRATION domain at res=40. Projected (30/40)^4~=3.16x speedup "
      f"would put res=30 at ~2.3-3.8 hrs for a comparable domain.", flush=True)
