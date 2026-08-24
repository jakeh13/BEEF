"""TIMING PROBE -- tests decay_method="energy" on the real production
input-section solve, at the exact config currently running in job 12031838
(res=30, sim_width=16.0, sim_height=4.0). See
investigate_mmi1x4_energy_decay_calibration.py's module docstring for the
full 4x4-90-hybrid-bug context this is checking against.

This section is the strongest candidate for the same pathology: its own
decay-timestep count was IDENTICAL (23,388 steps, t=389.8) across two very
different domain sizes (sim_width/height 21x6 and 17x5) -- exactly the kind
of "insensitive to changes that should matter" symptom a single stuck DFT
component (not genuine bulk-energy physics) would produce.

Direct comparison point: job 12031838's own input-section result at this
identical config (decay_method="dft", the default) -- 23,388 timesteps,
t=389.8.

save_dir is NOT set here (unlike the earlier res=30 input timing probe) --
decay_method is part of this solver's own cache fingerprint
(beef/cache/config.py, added 2026-08-18 alongside the flag itself), so a
save_dir=None probe avoids writing a second, differently-keyed cache entry
before this method is validated for correctness, not just speed.

Run:
    python .../investigate_mmi1x4_energy_decay_input.py
"""
import os
import time
import numpy as np

from beef import MMI, FDTDSolver, SIN_300

RESOLUTION = 30.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 3.1, 12.0
MM_LENGTH = 39.90
N_MODES = 12
STRAIGHT_LENGTH = 0.5
body_stub = 0.25
DECAY_TOL = 1e-11

FDTD_SIM_WIDTH = MM_WIDTH + 4.0   # 16.0, current production value
FDTD_SIM_HEIGHT = 4.0

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30, use_neighbor_aperture=True,
)

_fdtd_common = dict(
    res=RESOLUTION, sim_width=FDTD_SIM_WIDTH, sim_height=FDTD_SIM_HEIGHT,
    dpml=0.5, source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SIN_300, decay_tol=DECAY_TOL,
    decay_method="energy", energy_decay_dt=20.0,
)
input_solver = FDTDSolver(**_fdtd_common)

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
boundaries = [-half_total, -half_body + body_stub, +half_body - body_stub, +half_total]
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

print(f"ENERGY-DECAY INPUT-SECTION PROBE, res={RESOLUTION}, decay_tol={DECAY_TOL}, "
      f"decay_method=energy, energy_decay_dt=20.0, N_MODES={N_MODES}, "
      f"sim_width={FDTD_SIM_WIDTH}, sim_height={FDTD_SIM_HEIGHT}", flush=True)
print(f"Comparison point (decay_method='dft', SAME config): job 12031838's own "
      f"input-section result -- 23,388 timesteps, t=389.8", flush=True)

t0 = time.time()
M_in1 = input_solver.solve_section(
    mmi, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=input_solver.waveguide_length, stub_right=0.0,
    input_idx=0, save_dir=None,
)
elapsed = time.time() - t0
print(f"\nINPUT SECTION (energy decay) DONE: {elapsed:.1f} s ({elapsed/3600:.3f} hrs)",
      flush=True)
print(f"Transfer matrix |T| range: [{np.abs(M_in1).min():.4e}, {np.abs(M_in1).max():.4e}]",
      flush=True)
