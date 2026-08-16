"""Step 4 of validation_markdown/resolved/section_mode_beta_multi_monitor_validation.md.

Numerical equivalence test: old path (3 independent solve_section() calls)
vs. new path (one FDTDSolver.solve_section_multi_length() call), same MMI
2x2 body, same lengths, same seed_beta. REDUCED scale per the doc's step 4
budget (n_modes=2, 1 wavelength, real res=40/decay_tol=1e-11/real lengths --
NOT full N_MODES=14/11-wavelength production scale). Real cost estimate:
old ~3.5-4h, merged ~1.2-1.5h, ~5-6h combined -- cluster-gated, not a
login-node run. See batch/verification/run_equivalence_section_multi_length.batch.

Part A compares raw per-length phases and fitted beta between the two
paths. Part B checks the linearity/branch-margin guards still fire
correctly: B1 is a pure-math sanity check of beta_fit_from_phases() with no
FDTD involved; B2 injects a deliberate outlier into Part A's real phases
(both paths) and confirms the residual guard fires on both; B3 calls the
real section_mode_beta() entry point with deliberately too-widely-spaced
lengths and confirms its pre-simulation branch-margin guard (Guard 2) still
raises before any FDTD runs. B3 only exercises the old path's wired guard --
solve_section_multi_length() has no guard wrapper of its own yet (that
wiring is step 7, not done), which this script calls out explicitly rather
than pretending symmetric coverage exists.

Run:
    python /path/to/examples/verification/python/equivalence_test_solve_section_multi_length.py
"""
import time

import numpy as np

from beef import MMI, FDTDSolver, SOI_220
from beef.solvers.fdtd import (
    beta_fit_from_phases,
    _BETA_MIN_BRANCH_FRACTION,
    _BETA_FIT_RESIDUAL_MAX,
    _BETA_BRANCH_MARGIN_RAISE,
)

RESOLUTION = 40.0
WL_CENTER = 0.5 * (1.530 + 1.565)  # production band centre, single wavelength here
N_MODES = 2  # reduced from production's 14, per step 4's cost budget

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
span = x1 - x0

solver = FDTDSolver(
    res=RESOLUTION, sim_width=5.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SOI_220, decay_tol=1e-11, half_domain_aperture=True,
)

wavelengths = np.array([WL_CENTER])
iface = dict(monitor_center=0.0, monitor_size=solver.sim_width)

print(f"Body span: [{x0:.4f}, {x1:.4f}] um, mid={mid:.4f} um")
print(f"n_modes={N_MODES}, wavelengths={wavelengths.tolist()}, res={RESOLUTION}, "
      f"decay_tol={solver.decay_tol}\n")

# ---------------------------------------------------------------------------
# Pick lengths exactly the way section_mode_beta() does internally (its own
# `if lengths is None` block, fdtd.py:2026-2032) -- so Part A's lengths are
# the real ones this n_modes=2 call would auto-pick, not hand-picked values.
# ---------------------------------------------------------------------------
seed_span = min(2.0, span - 2.0 * solver.dpml)
seed = solver.section_mode_neff(
    mmi, wavelengths, mid - seed_span / 2.0, mid + seed_span / 2.0, N_MODES, **iface,
)
seed_beta = 2.0 * np.pi * seed / wavelengths[:, None]
beta_max = float(np.max(np.abs(seed_beta)))
max_safe_spacing = 2.0 * np.pi / (_BETA_MIN_BRANCH_FRACTION * beta_max)
usable_span = span - 2.0 * solver.dpml
spacing = min(usable_span / 3.0, max_safe_spacing / 2.0)
lengths = [spacing, 2.0 * spacing, 3.0 * spacing]

print(f"seed n_eff={seed.tolist()}, beta_max={beta_max:.4f} rad/um, "
      f"max_safe_spacing={max_safe_spacing:.4f} um")
print(f"lengths (auto-picked, n_modes={N_MODES}): {[f'{L:.4f}' for L in lengths]} um\n")

# ---------------------------------------------------------------------------
# Part A: old path (3 independent solve_section() calls)
# ---------------------------------------------------------------------------
print("=== Part A: old path (3 independent solve_section() calls) ===")
phases_old = np.zeros((len(lengths), len(wavelengths), N_MODES))
T_old = [None] * len(lengths)
times_old = []
for j, L in enumerate(lengths):
    t0 = time.time()
    T = solver.solve_section(
        mmi, wavelengths, mid - L / 2.0, mid + L / 2.0,
        n_modes_in=N_MODES, n_modes_out=N_MODES,
        save_dir=None, **iface,
    )
    elapsed = time.time() - t0
    times_old.append(elapsed)
    T_old[j] = T
    phases_old[j] = np.angle(np.diagonal(T, axis1=1, axis2=2))
    print(f"  length[{j}]={L:.4f} um: {elapsed:.1f} s")
old_total = sum(times_old)
print(f"Old path total: {old_total:.1f} s ({old_total/3600:.2f} h)\n")

# ---------------------------------------------------------------------------
# Part A: merged path (one solve_section_multi_length() call)
# ---------------------------------------------------------------------------
print("=== Part A: merged path (one solve_section_multi_length() call) ===")
t0 = time.time()
T_multi = solver.solve_section_multi_length(
    mmi, wavelengths, mid, lengths, N_MODES, **iface,
)
merged_total = time.time() - t0
print(f"Merged path total: {merged_total:.1f} s ({merged_total/3600:.2f} h)\n")

phases_merged = np.zeros((len(lengths), len(wavelengths), N_MODES))
for j in range(len(lengths)):
    phases_merged[j] = np.angle(np.diagonal(T_multi[j], axis1=1, axis2=2))

print(f"Measured speedup (this reduced scale): {old_total/merged_total:.2f}x "
      f"(doc step 1 projected 3.00x at full production scale)\n")

# ---------------------------------------------------------------------------
# Part A: compare raw phases and |T_ii| between the two paths
# ---------------------------------------------------------------------------
phase_diff = np.angle(np.exp(1j * (phases_old - phases_merged)))  # wrap to (-pi, pi]
max_phase_diff = float(np.max(np.abs(phase_diff)))

amp_old = np.stack([np.abs(np.diagonal(T, axis1=1, axis2=2)) for T in T_old])
amp_merged = np.stack([np.abs(np.diagonal(T_multi[j], axis1=1, axis2=2)) for j in range(len(lengths))])
max_amp_diff = float(np.max(np.abs(amp_old - amp_merged)))

print("=== Part A: raw phase/amplitude agreement ===")
print(f"max|phase_old - phase_merged| = {max_phase_diff:.6e} rad "
      f"(guard's own straight-line tolerance _BETA_FIT_RESIDUAL_MAX = {_BETA_FIT_RESIDUAL_MAX} rad, "
      f"for scale)")
print(f"max||T_ii|_old - |T_ii|_merged| = {max_amp_diff:.6e}\n")

# ---------------------------------------------------------------------------
# Part A: fit beta on both phase sets with the SAME seed, compare
# ---------------------------------------------------------------------------
print("=== Part A: fitted beta agreement (beta_fit_from_phases, shared seed) ===")
beta_old = np.zeros((len(wavelengths), N_MODES))
beta_merged = np.zeros((len(wavelengths), N_MODES))
resid_old = np.zeros((len(wavelengths), N_MODES))
resid_merged = np.zeros((len(wavelengths), N_MODES))
margin_old = np.zeros((len(wavelengths), N_MODES))
margin_merged = np.zeros((len(wavelengths), N_MODES))

for i in range(len(wavelengths)):
    for m in range(N_MODES):
        fit_o = beta_fit_from_phases(lengths, phases_old[:, i, m], seed_beta[i, m])
        fit_n = beta_fit_from_phases(lengths, phases_merged[:, i, m], seed_beta[i, m])
        beta_old[i, m], resid_old[i, m], margin_old[i, m] = fit_o.beta, fit_o.resid, fit_o.margin
        beta_merged[i, m], resid_merged[i, m], margin_merged[i, m] = fit_n.beta, fit_n.resid, fit_n.margin
        print(f"  wl[{i}] mode[{m}]: beta_old={fit_o.beta:.6f}  beta_merged={fit_n.beta:.6f}  "
              f"|diff|={abs(fit_o.beta - fit_n.beta):.3e} rad/um  "
              f"(resid_old={fit_o.resid:.2e}, resid_merged={fit_n.resid:.2e}, "
              f"margin_old={fit_o.margin:.3f}, margin_merged={fit_n.margin:.3f})")

max_beta_diff = float(np.max(np.abs(beta_old - beta_merged)))
max_beta_diff_pct = float(np.max(np.abs((beta_old - beta_merged) / beta_old)) * 100.0)
print(f"\nmax|beta_old - beta_merged| = {max_beta_diff:.6e} rad/um "
      f"({max_beta_diff_pct:.4f}%)")
print(f"worst residual either path: {max(float(resid_old.max()), float(resid_merged.max())):.3e} rad "
      f"(guard fires above {_BETA_FIT_RESIDUAL_MAX})")
print(f"worst margin either path: {max(float(margin_old.max()), float(margin_merged.max())):.3f} "
      f"(guard fires above {_BETA_BRANCH_MARGIN_RAISE})\n")

# ---------------------------------------------------------------------------
# Part B1: pure-math sanity check of beta_fit_from_phases() -- no FDTD.
# Synthetic phases from a KNOWN beta at the SAME real lengths; recovered
# beta should match to float precision.
# ---------------------------------------------------------------------------
print("=== Part B1: synthetic sanity check (no FDTD) ===")
true_beta_synth = 9.4  # rad/um, representative silicon-strip scale
L_arr = np.array(lengths)
phi_synth = np.angle(np.exp(1j * true_beta_synth * L_arr))
fit_synth = beta_fit_from_phases(lengths, phi_synth, true_beta_synth * 0.98)  # seed off by 2%
b1_pass = abs(fit_synth.beta - true_beta_synth) < 1e-9
print(f"true beta={true_beta_synth}, recovered beta={fit_synth.beta:.12f}, "
      f"|diff|={abs(fit_synth.beta - true_beta_synth):.3e}")
print(f"B1 {'PASS' if b1_pass else 'FAIL'}: beta_fit_from_phases() recovers a known slope "
      f"at these lengths to float precision.\n")

# ---------------------------------------------------------------------------
# Part B2: residual guard, reusing Part A's REAL phases (no extra FDTD).
# Inject a deliberate outlier at the middle length for wl=0, mode=0 in both
# phase sets, confirm the residual guard fires on data from EACH path.
# ---------------------------------------------------------------------------
print("=== Part B2: residual guard fires on real (outlier-injected) data ===")
phases_old_bad = phases_old.copy()
phases_merged_bad = phases_merged.copy()
phases_old_bad[1, 0, 0] = np.angle(np.exp(1j * (phases_old_bad[1, 0, 0] + 1.0)))
phases_merged_bad[1, 0, 0] = np.angle(np.exp(1j * (phases_merged_bad[1, 0, 0] + 1.0)))

fit_old_bad = beta_fit_from_phases(lengths, phases_old_bad[:, 0, 0], seed_beta[0, 0])
fit_merged_bad = beta_fit_from_phases(lengths, phases_merged_bad[:, 0, 0], seed_beta[0, 0])
b2_old_fires = fit_old_bad.resid > _BETA_FIT_RESIDUAL_MAX
b2_merged_fires = fit_merged_bad.resid > _BETA_FIT_RESIDUAL_MAX
print(f"old-path data + outlier:    resid={fit_old_bad.resid:.4f} rad -> "
      f"guard {'FIRES (correct)' if b2_old_fires else 'DID NOT FIRE (unexpected)'}")
print(f"merged-path data + outlier: resid={fit_merged_bad.resid:.4f} rad -> "
      f"guard {'FIRES (correct)' if b2_merged_fires else 'DID NOT FIRE (unexpected)'}")
print(f"B2 {'PASS' if (b2_old_fires and b2_merged_fires) else 'FAIL'}: residual guard "
      f"behaves consistently across both paths' real data.\n")

# ---------------------------------------------------------------------------
# Part B3: section_mode_beta()'s own pre-simulation branch-margin guard
# (Guard 2, fdtd.py ~2038-2047) -- cheap, only needs the seed solve, no
# time-stepping FDTD. Deliberately spaces lengths beyond max_safe_spacing.
# Only exercises the OLD path's wired entry point -- solve_section_multi_
# length() has no equivalent wrapper yet (step 7 wiring not done).
# ---------------------------------------------------------------------------
print("=== Part B3: section_mode_beta()'s static branch-margin guard (Guard 2) ===")
bad_lengths = [max_safe_spacing * 2.0, max_safe_spacing * 4.0]
print(f"Deliberately too-widely-spaced lengths: {[f'{L:.4f}' for L in bad_lengths]} um "
      f"(max_safe_spacing={max_safe_spacing:.4f} um)")
try:
    solver.section_mode_beta(
        mmi, wavelengths, mid - bad_lengths[-1], mid + bad_lengths[-1],
        n_modes=N_MODES, lengths=bad_lengths, save_dir=None,
    )
    b3_fires = False
    b3_msg = "(no exception raised)"
except ValueError as e:
    b3_fires = True
    b3_msg = str(e)
print(f"B3 {'PASS' if b3_fires else 'FAIL'}: {'guard raised ValueError as expected' if b3_fires else 'guard did NOT raise (unexpected)'}")
print(f"  message: {b3_msg}\n")
print("Note: B3 only exercises the OLD path's wired guard (section_mode_beta() -> "
      "solve_section()); solve_section_multi_length() has no guard wrapper of its "
      "own yet -- that wiring is step 7, not this step.\n")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("=== SUMMARY ===")
print(f"Old path: {old_total:.1f} s ({old_total/3600:.2f} h) | "
      f"Merged path: {merged_total:.1f} s ({merged_total/3600:.2f} h) | "
      f"Speedup: {old_total/merged_total:.2f}x")
print(f"max|phase_old - phase_merged| = {max_phase_diff:.6e} rad")
print(f"max||T_ii| diff| = {max_amp_diff:.6e}")
print(f"max|beta_old - beta_merged| = {max_beta_diff:.6e} rad/um ({max_beta_diff_pct:.4f}%)")
print(f"B1 (synthetic sanity): {'PASS' if b1_pass else 'FAIL'}")
print(f"B2 (residual guard, both paths): {'PASS' if (b2_old_fires and b2_merged_fires) else 'FAIL'}")
print(f"B3 (branch-margin guard, old path entry point): {'PASS' if b3_fires else 'FAIL'}")
