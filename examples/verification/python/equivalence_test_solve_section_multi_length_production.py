"""Production-scale equivalence test for the merged multi-monitor beta path.

resolved/section_mode_beta_multi_monitor_validation.md's step 4 established
old-vs-merged equivalence at REDUCED scale (n_modes=2, 1 wavelength):
`max|beta diff| = 3.50e-05 rad/um`, which that doc explicitly labels "strong
evidence, not a production-scale guarantee". Steps 5-6's convergence studies
carry the same reduced-scale caveat. This script closes that gap on the MMI
2x2 body at the real shipped scale: N_MODES=14, the 11-point C-band sweep,
res=40, decay_tol=1e-11 -- the exact configuration mmi_2x2_hybrid.py ships.

Both paths are reached through `section_mode_beta_details()`, not by hand:

    per-length  save_dir=CACHE, cache_only=True  -> merged_active False
    merged      save_dir=None                    -> merged_active True

so each gets its real Guard 2, its real fit, and its real post-fit guards,
and the comparison is between two things production could actually run
rather than two hand-rolled approximations of them.

**The per-length side costs zero FDTD.** Its three section solves were paid
for by job 11872191 (the step-8 cold run, 2026-08-15/16) and are still in
example_sim_cache/; `cache_only=True` demands them rather than allowing a
silent 30-hour re-simulation. Only the merged side actually time-steps.

`strict=False` on both calls: a guard failure is precisely the case worth
comparing rather than aborting on, and any failure lands in
`.guard_failures` and is reported below. Neither returned beta is used to
propagate anything here.

Run (48 ranks, matching job 11872191 so the timing comparison is
apples-to-apples):
    srun --mpi=pmi2 python3 examples/verification/python/equivalence_test_solve_section_multi_length_production.py
"""

import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, FDTDSolver, SOI_220, plot_beta_calibration
from beef.solvers.fdtd import _BETA_FIT_RESIDUAL_MAX, _BETA_BRANCH_MARGIN_RAISE

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
CACHE = os.path.join(_ROOT, "example_sim_cache")
output_dir = os.path.join(_ROOT, "example_output", "verification",
                          "mmi_2x2_beta_multilength_production")
os.makedirs(output_dir, exist_ok=True)

# --- production configuration, copied from mmi_2x2_hybrid.py ---------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.530, 1.565, 11
N_MODES = 14
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

mmi = MMI(
    n_inputs=2, n_outputs=2,
    mm_length=28.0, mm_width=2.8,
    waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
    waveguide_offset=1.5, waveguide_height=0.22,
)

half_body = mmi.mm_length / 2
body_stub = 0.25
x0, x1 = -half_body + body_stub, +half_body - body_stub   # the EME body span
body_length_nominal = mmi.mm_length - 2 * body_stub       # 27.5 um, what the phase rides

solver = FDTDSolver(
    res=RESOLUTION, sim_width=5.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SOI_220, decay_tol=1e-11, half_domain_aperture=True,
)

# Per-length wall-clock baseline for the speedup number, read off the cache
# folders job 11872191 wrote rather than off its log (which prints only one
# total elapsed): the three FDTDSolver_solve_section_MMI_* beta entries are
# stamped 2026-08-15 22:21, 2026-08-16 08:03, and 2026-08-16 18:15, against
# an 11:37 start for the first stage -- 10.7 h, 9.7 h, 10.2 h, i.e. 30.6 h
# for the calibration alone. Nearly flat across the three lengths, which is
# itself the reason merging should win ~3x: the decay wait dominates, not
# the propagation distance.
PER_LENGTH_BASELINE_H = 30.6

print("=" * 78)
print("Production-scale equivalence: per-length vs merged multi-monitor beta")
print("=" * 78)
print(f"Body span [{x0:.4f}, {x1:.4f}] um  (nominal body length {body_length_nominal:.2f} um)")
print(f"N_MODES={N_MODES}, {N_WL} wavelengths {WL_MIN}-{WL_MAX} um, "
      f"res={RESOLUTION}, decay_tol={solver.decay_tol}")
print(f"Cache: {CACHE}\n")

# ---------------------------------------------------------------------------
# Path 1 -- per-length, served from cache. Zero FDTD.
# ---------------------------------------------------------------------------
print("=== Path 1: per-length (3 independent domains), from cache ===")
t0 = time.time()
cal_old = solver.section_mode_beta_details(
    mmi, wavelengths, x0, x1, N_MODES,
    save_dir=CACHE, cache_only=True, strict=False,
)
t_old = time.time() - t0
print(f"  merged flag on returned calibration: {cal_old.merged}  (expected False)")
print(f"  lengths: {[f'{L:.6f}' for L in cal_old.lengths]} um")
print(f"  wall: {t_old:.1f} s  (cache hits + the MPB seed, which is never cached)")
print(f"  guard_failures: {list(cal_old.guard_failures) or 'none'}\n")

# ---------------------------------------------------------------------------
# Path 2 -- merged. This is the part that actually simulates.
# ---------------------------------------------------------------------------
print("=== Path 2: merged (one domain, one monitor per length) ===")
t0 = time.time()
cal_new = solver.section_mode_beta_details(
    mmi, wavelengths, x0, x1, N_MODES,
    save_dir=None, use_merged_domain=True, strict=False,
)
t_new = time.time() - t0
print(f"  merged flag on returned calibration: {cal_new.merged}  (expected True)")
print(f"  lengths: {[f'{L:.6f}' for L in cal_new.lengths]} um")
print(f"  wall: {t_new:.1f} s ({t_new/3600:.2f} h)")
print(f"  guard_failures: {list(cal_new.guard_failures) or 'none'}\n")

# Bank the ~10 h of time-stepping before anything downstream can throw. The
# analysis below is pure numpy/matplotlib, but a crash there would otherwise
# discard the only expensive thing this job produced.
np.savez(
    os.path.join(output_dir, "merged_raw.npz"),
    beta=cal_new.beta, phases=np.asarray(cal_new.phases),
    lengths=np.asarray(cal_new.lengths), seed_beta=cal_new.seed_beta,
    wavelengths=wavelengths, wall_s=t_new,
)
print(f"  raw merged arrays banked -> {output_dir}/merged_raw.npz\n")

same_lengths = np.allclose(cal_old.lengths, cal_new.lengths, rtol=0, atol=1e-12)
print(f"Both paths picked identical fit lengths: {same_lengths}  "
      f"(they share the same seed solve and the same selection code; a False "
      f"here would invalidate every comparison below)\n")

# ---------------------------------------------------------------------------
# Agreement
# ---------------------------------------------------------------------------
d_beta = cal_new.beta - cal_old.beta
max_d_beta = float(np.max(np.abs(d_beta)))
i_w, i_m = np.unravel_index(int(np.argmax(np.abs(d_beta))), d_beta.shape)
max_d_beta_pct = float(np.max(np.abs(d_beta / cal_old.beta)) * 100.0)

d_phase = np.angle(np.exp(1j * (np.asarray(cal_new.phases) - np.asarray(cal_old.phases))))
max_d_phase = float(np.max(np.abs(d_phase)))

# What the disagreement is worth where it is actually used: the body phase.
d_body_phase = max_d_beta * body_length_nominal

# ...and whether it would even change the EME body sweep's cache key, which
# stores neff rounded to 9 decimals (beef/cache/config.py _json_safe).
d_neff = np.abs(d_beta) * wavelengths[:, None] / (2 * np.pi)
max_d_neff = float(np.max(d_neff))

print("=== Agreement (merged - per-length) ===")
print(f"max|d_beta|            = {max_d_beta:.6e} rad/um  ({max_d_beta_pct:.6f}%)"
      f"  at lambda={wavelengths[i_w]:.4f} um, mode {i_m}")
print(f"max|d_phase| (raw)     = {max_d_phase:.6e} rad")
print(f"max|d_neff|            = {max_d_neff:.6e}"
      f"   -> cache key (9 dp) {'WOULD' if max_d_neff >= 5e-10 else 'would NOT'} change")
print(f"accumulated over the {body_length_nominal:.1f} um body: {d_body_phase:.6e} rad"
      f"  ({np.degrees(d_body_phase):.4f} deg)\n")

print("  for scale:")
print(f"    step 4's reduced-scale equivalence floor : 3.50e-05 rad/um")
print(f"    step 6's domain-length effect           : 5.78e-04 rad/um")
print(f"    this run                                : {max_d_beta:.2e} rad/um\n")

print("=== Per-mode worst disagreement across the band ===")
print(f"  {'mode':>4}  {'beta_per_length':>16}  {'beta_merged':>14}  {'max|d_beta|':>12}  {'%':>10}")
for m in range(N_MODES):
    col = np.abs(d_beta[:, m])
    k = int(np.argmax(col))
    print(f"  {m:>4}  {cal_old.beta[k, m]:>16.6f}  {cal_new.beta[k, m]:>14.6f}  "
          f"{col[k]:>12.3e}  {col[k]/cal_old.beta[k, m]*100:>10.6f}")
print()

# ---------------------------------------------------------------------------
# Guards, both paths
# ---------------------------------------------------------------------------
def _grids(cal):
    n_wl, n_modes = cal.beta.shape
    resid = np.array([[cal.fits[i][m].resid for m in range(n_modes)] for i in range(n_wl)])
    margin = np.array([[cal.fits[i][m].margin for m in range(n_modes)] for i in range(n_wl)])
    return resid, margin

resid_old, margin_old = _grids(cal_old)
resid_new, margin_new = _grids(cal_new)

print("=== Guards ===")
print(f"  worst residual  per-length {resid_old.max():.3e} | merged {resid_new.max():.3e} "
      f"rad   (guard fires above {_BETA_FIT_RESIDUAL_MAX})")
print(f"  worst margin    per-length {margin_old.max():.3f} | merged {margin_new.max():.3f} "
      f"       (guard fires above {_BETA_BRANCH_MARGIN_RAISE})")
resid_ratio = resid_new.max() / max(resid_old.max(), 1e-30)
print(f"  merged/per-length worst-residual ratio: {resid_ratio:.2f}x"
      f"  (step 5 saw residual grow with monitor count at reduced scale;"
      f" this is that effect at production scale)\n")

# ---------------------------------------------------------------------------
# Cost
# ---------------------------------------------------------------------------
print("=== Cost ===")
print(f"  per-length, if actually simulated : {PER_LENGTH_BASELINE_H:.1f} h "
      f"(job 11872191, measured from cache-folder timestamps)")
print(f"  merged, this run                  : {t_new/3600:.2f} h")
print(f"  speedup                           : {PER_LENGTH_BASELINE_H/(t_new/3600):.2f}x"
      f"  (doc step 1 projected 3.00x: 42.1 h -> 14.0 h)\n")

# ---------------------------------------------------------------------------
# Artifacts
# ---------------------------------------------------------------------------
np.savez(
    os.path.join(output_dir, "beta_comparison.npz"),
    wavelengths=wavelengths, lengths=np.asarray(cal_new.lengths),
    beta_per_length=cal_old.beta, beta_merged=cal_new.beta,
    phases_per_length=np.asarray(cal_old.phases), phases_merged=np.asarray(cal_new.phases),
    seed_beta=cal_old.seed_beta,
    resid_per_length=resid_old, resid_merged=resid_new,
    margin_per_length=margin_old, margin_merged=margin_new,
    wall_merged_s=t_new, per_length_baseline_h=PER_LENGTH_BASELINE_H,
)
print(f"Arrays saved -> {output_dir}/beta_comparison.npz")

plot_beta_calibration(
    cal_new, mmi, solver, os.path.join(output_dir, "calibration_merged"),
    monitor_center=0.0, monitor_size=solver.sim_width,
)
print(f"Merged-path calibration plots -> {output_dir}/calibration_merged/"
      f"  (compare against example_output/mmi_2x2/mmi_2x2_hybrid/calibration/)")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
colors = plt.cm.viridis(np.linspace(0, 1, N_MODES))

ax = axes[0]
for m in range(N_MODES):
    ax.plot(wavelengths, d_beta[:, m], color=colors[m], marker="o", ms=3, lw=1.0)
ax.axhline(0.0, color="gray", lw=0.8)
for y, lab in [(3.50e-05, "step 4 floor"), (5.78e-04, "step 6 domain effect")]:
    ax.axhline(y, color="crimson", ls=":", lw=0.9)
    ax.axhline(-y, color="crimson", ls=":", lw=0.9)
    ax.annotate(lab, (wavelengths[0], y), fontsize=7, color="crimson",
                va="bottom", ha="left")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel(r"$\beta_{merged} - \beta_{per-length}$  (rad/µm)")
ax.set_title("Merged vs per-length β, all 14 modes")
ax.grid(True, alpha=0.4)

ax = axes[1]
im = ax.imshow(np.abs(d_beta), aspect="auto", origin="lower", cmap="magma",
               extent=[-0.5, N_MODES - 0.5, wavelengths[0], wavelengths[-1]])
ax.set_xlabel("mode index")
ax.set_ylabel("Wavelength (µm)")
ax.set_title(r"$|\Delta\beta|$ (rad/µm)")
plt.colorbar(im, ax=ax)

plt.suptitle(
    f"MMI 2×2 body β — merged multi-monitor vs per-length, production scale "
    f"(N_MODES={N_MODES}, {N_WL} λ)  |  max|Δβ| = {max_d_beta:.2e} rad/µm",
    fontsize=10,
)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "beta_merged_vs_per_length.png"), dpi=300)
plt.close("all")
print(f"Comparison figure -> {output_dir}/beta_merged_vs_per_length.png\n")

# ---------------------------------------------------------------------------
print("=== SUMMARY ===")
print(f"identical fit lengths      : {same_lengths}")
print(f"max|d_beta|                : {max_d_beta:.6e} rad/um ({max_d_beta_pct:.6f}%)")
print(f"max|d_neff|                : {max_d_neff:.6e}")
print(f"body-phase consequence     : {d_body_phase:.6e} rad over {body_length_nominal:.1f} um")
print(f"worst residual (merged)    : {resid_new.max():.3e} rad  vs guard {_BETA_FIT_RESIDUAL_MAX}")
print(f"worst margin   (merged)    : {margin_new.max():.3f}      vs guard {_BETA_BRANCH_MARGIN_RAISE}")
print(f"guard failures             : per-length {list(cal_old.guard_failures) or 'none'} | "
      f"merged {list(cal_new.guard_failures) or 'none'}")
print(f"merged wall                : {t_new/3600:.2f} h  "
      f"(vs {PER_LENGTH_BASELINE_H:.1f} h per-length, "
      f"{PER_LENGTH_BASELINE_H/(t_new/3600):.2f}x)")
