"""Step 1 of validation_markdown/section_mode_beta_multi_monitor_validation.md.

Times ONE real FDTDSolver.solve_section() call -- one excited band
(n_modes_in=1), one length, at the MMI 2x2 body's actual calibration window
(mid=0.0, L=2.2046 um -- the real middle length measured in
manuscript_figures/solver_figures/calibrate_body_neff_mmi_2x2_geometry.ipynb).
n_modes_out=N_MODES=14 kept full to match production's real output-side
decomposition cost -- only the excited (input) side is isolated to 1 band,
since that is the loop unit section_mode_beta()/solve_section() actually
iterates over internally.

Projects the old (3-lengths-independent) vs. new (merged multi-monitor,
FDTDSolver.solve_section_multi_length()) path costs from this single real
measurement. See the validation doc's step 1 for the recorded result
(3609.53 s on a local M3 MacBook, 2026-08-12) and step 4 for what still
needs a real run at this scale: the actual numerical equivalence test
between the two paths.

Run:
    python /path/to/examples/verification/python/time_section_mode_beta_single_sim.py
"""
import time

import numpy as np

from beef import MMI, FDTDSolver, SOI_220

RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.530, 1.565, 11
N_MODES = 14

mmi = MMI(
    n_inputs=2, n_outputs=2,
    mm_length=28.0, mm_width=2.8,
    waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
    waveguide_offset=1.5, waveguide_height=0.22,
)

half_body = mmi.mm_length / 2
body_stub = 0.25
boundaries = [None, -half_body + body_stub, +half_body - body_stub, None]
x0, x1 = boundaries[1], boundaries[2]
mid = 0.5 * (x0 + x1)

input_solver = FDTDSolver(
    res=RESOLUTION, sim_width=5.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SOI_220, decay_tol=1e-11, half_domain_aperture=True,
)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

# Real calibration length already measured
# (calibrate_body_neff_mmi_2x2_geometry.ipynb): beta_max=11.4001 rad/um ->
# lengths=[1.102, 2.205, 3.307] um. Using the middle one.
L = 2.2046
sec_x0, sec_x1 = mid - L / 2.0, mid + L / 2.0

print(f"Body span: [{x0:.3f}, {x1:.3f}] um, mid={mid:.3f} um")
print(f"Calibration length: L={L:.4f} um -> section [{sec_x0:.4f}, {sec_x1:.4f}] um")
print(f"Timing: n_modes_in=1 (one excited band), n_modes_out={N_MODES}, "
      f"{N_WL} wavelengths, save_dir=None (no cache)\n")

t0 = time.time()
T = input_solver.solve_section(
    mmi, wavelengths, sec_x0, sec_x1,
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=0.0, stub_right=0.0,
    monitor_center=0.0, monitor_size=input_solver.sim_width,
    save_dir=None,
)
elapsed = time.time() - t0

print(f"\nT shape: {T.shape}")
print(f"Elapsed (one excited band, one length): {elapsed:.2f} s")

n_bands_mmi = N_MODES  # MMI 2x2 body calibration excites all 14 bands
n_lengths_mmi = 3
old_total = n_bands_mmi * n_lengths_mmi * elapsed
new_total_projected = n_bands_mmi * elapsed  # merged: 1 domain-equivalent build per band, not 3
print(f"\nProjected OLD path total (3 lengths x {n_bands_mmi} bands): {old_total:.1f} s "
      f"({old_total/60:.1f} min)")
print(f"Projected NEW (merged) path total ({n_bands_mmi} bands, 1 domain each): "
      f"{new_total_projected:.1f} s ({new_total_projected/60:.1f} min)")
print(f"Projected speedup: {old_total/new_total_projected:.2f}x "
      f"(NOTE: assumes per-band cost at the longest merged domain equals this "
      f"single-length measurement -- not yet verified, see doc step 1's own caveat)")
