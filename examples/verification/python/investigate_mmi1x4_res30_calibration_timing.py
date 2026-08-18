"""TIMING PROBE, not a full run -- times ONE excited mode's real
section_mode_beta_details() cost (seed + merged-domain decay + diagonal
decompose) at RESOLUTION=30, before committing the full N_MODES=12
calibration blind.

n_modes=1 here, NOT production's own N_MODES=12 -- this makes the probe
fast (one excited-mode simulation instead of twelve) and, since seed/length
auto-selection depends on beta_max (dominated by the fundamental mode
either way), should closely approximate the FIRST mode's own real cost at
res=30. It does NOT capture mode-to-mode variation across the full 12-mode
spectrum (near-cutoff higher-order modes could plausibly decay differently)
-- flagged as a real, unresolved source of uncertainty in the resulting
projection, not swept under the rug. NOT cached (n_modes=1 doesn't match
production's own n_modes=12 cache key) -- purely a timing measurement,
thrown away after this prints its number.

Run:
    python .../investigate_mmi1x4_res30_calibration_timing.py
"""
import os
import time
import numpy as np

from beef import MMI, FDTDSolver, SIN_300

RESOLUTION = 30.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 3.1, 12.0
MM_LENGTH = 39.90
STRAIGHT_LENGTH = 0.5
body_stub = 0.25
DECAY_TOL = 1e-11
N_MODES_PROBE = 1   # see module docstring -- NOT production's N_MODES=12

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
input_solver = FDTDSolver(**_fdtd_common)   # the calibration SOURCE, same solver class as production

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
boundaries = [-half_total, -half_body + body_stub, +half_body - body_stub, +half_total]
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

x0, x1 = boundaries[1], boundaries[2]   # the body's own span, exactly what _calibrated_neff() passes

print(f"CALIBRATION BAND TIMING PROBE, res={RESOLUTION}, decay_tol={DECAY_TOL}, "
      f"n_modes(probe)={N_MODES_PROBE} (production uses N_MODES=12), N_WL={N_WL}",
      flush=True)
print(f"Body span: x0={x0:.4f}, x1={x1:.4f} um (span={x1-x0:.4f} um)", flush=True)

t0 = time.time()
calib = input_solver.section_mode_beta_details(
    mmi, wavelengths, x0, x1, N_MODES_PROBE,
    save_dir=None,   # NOT cached -- see module docstring
)
elapsed = time.time() - t0
print(f"\nCALIBRATION PROBE DONE: {elapsed:.1f} s ({elapsed/3600:.3f} hrs) "
      f"for {N_MODES_PROBE} mode(s)", flush=True)
print(f"lengths auto-selected: {calib.lengths}", flush=True)
print(f"seed_wavelengths: {calib.seed_wavelengths}", flush=True)
print(f"fitted beta (rad/um): {calib.beta}", flush=True)
print(f"\nProjection to production N_MODES=12 (LINEAR extrapolation, real "
      f"per-mode variance unmeasured): {elapsed*12/3600:.2f} hrs for the full "
      f"calibration's decay+decompose work alone (seed step scales with "
      f"n_modes too but is cheap -- MPB-only, no time-stepping).", flush=True)
