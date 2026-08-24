"""TIMING + CORRECTNESS PROBE -- tests decay_method="energy" on ONE
calibration mode, at the exact config currently running in the real
production prewarm (job 12031838: res=30, sim_width=16.0, sim_height=4.0),
to check whether this device's own unexpectedly-long/non-monotonic decay
behavior shares the SAME root cause as the 4x4 90-degree hybrid's own
diagnosed bug (validation_markdown/mmi_4x4_90_hybrid_validation.md): a case
that NEVER converged past 660,000+ timesteps under decay_method="dft"
(meep's stop_when_dft_decayed, which requires EVERY field component/
frequency on EVERY DFT monitor to individually decay, no ceiling) converged
completely normally (11,207 steps) under decay_method="energy"
(stop_when_energy_decayed, a single bulk-energy check) -- the bulk energy
had decayed normally the whole time; one specific DFT component was stuck.

This device's own symptom is milder (real convergence, just more steps than
expected and non-monotonic with resolution/domain size) but the same family
of DFT-criterion pathology is a plausible explanation, especially given
this device's own densely-packed near-cutoff/near-degenerate modes.

Direct comparison point: job 12031838's own mode-1 result at this identical
config (decay_method="dft", the default) was 11,467 timesteps, t=191.117.

energy_decay_dt=20.0 -- same starting guess the 4x4 device's own
investigation used, no more principled a choice here than there (its own
docstring: "No principled default yet").

Per decay_method's own docstring caveat: "Investigative, not yet validated
as a general replacement" -- this probe checks BOTH whether it terminates
faster AND whether the resulting fitted beta is consistent with the
decay_method="dft" result already in hand, not speed alone.

Run:
    python .../investigate_mmi1x4_energy_decay_calibration.py
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
N_MODES_PROBE = 1   # matches the earlier res=30 calibration probe, for a clean comparison

# CURRENT production domain (2026-08-19 finding, mmi_1x4_sin_body_domain_size_convergence.py)
FDTD_SIM_WIDTH = MM_WIDTH + 4.0   # 16.0
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

x0, x1 = boundaries[1], boundaries[2]

print(f"ENERGY-DECAY CALIBRATION PROBE, res={RESOLUTION}, decay_tol={DECAY_TOL}, "
      f"decay_method=energy, energy_decay_dt=20.0, n_modes(probe)={N_MODES_PROBE}, "
      f"sim_width={FDTD_SIM_WIDTH}, sim_height={FDTD_SIM_HEIGHT}", flush=True)
print(f"Comparison point (decay_method='dft', SAME config): job 12031838's own "
      f"mode-1 result -- 11,467 timesteps, t=191.117", flush=True)
print(f"Body span: x0={x0:.4f}, x1={x1:.4f} um (span={x1-x0:.4f} um)", flush=True)

t0 = time.time()
calib = input_solver.section_mode_beta_details(
    mmi, wavelengths, x0, x1, N_MODES_PROBE,
    save_dir=None,   # NOT cached -- probe only, and decay_method is fingerprinted separately anyway
)
elapsed = time.time() - t0
print(f"\nENERGY-DECAY PROBE DONE: {elapsed:.1f} s ({elapsed/3600:.3f} hrs) "
      f"for {N_MODES_PROBE} mode(s)", flush=True)
print(f"lengths auto-selected: {calib.lengths}", flush=True)
print(f"fitted beta (rad/um): {calib.beta}", flush=True)
print(f"guard_failures: {calib.guard_failures}", flush=True)
for wl_i, wl in enumerate(wavelengths):
    for m in range(N_MODES_PROBE):
        fit = calib.fits[wl_i][m]
        print(f"  wl={wl:.4f} mode={m}: beta={fit.beta:.6f} resid={fit.resid:.4e} "
              f"margin={fit.margin:.4f} seed_excess={fit.seed_excess:.6f}", flush=True)
print(f"\nProjection to production N_MODES=12 (linear): {elapsed*12/3600:.2f} hrs "
      f"for calibration's decay+decompose work alone, IF this speedup and "
      f"correctness both hold across all 12 modes (unverified, single-mode probe).",
      flush=True)
