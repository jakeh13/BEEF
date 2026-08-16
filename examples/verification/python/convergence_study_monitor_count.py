"""Step 5 of validation_markdown/resolved/section_mode_beta_multi_monitor_validation.md.

Convergence study #1 -- number of monitors/lengths. Gated on step 4's
numerical equivalence PASS (2026-08-13, job 11813378): with merging
validated, adding a monitor is cheap (no new domain build/decay wait), so
this sweeps N in {3, 4, 5, 6} monitors in ONE merged
FDTDSolver.solve_section_multi_length() run and re-fits beta on growing
prefixes of that single run's phase data -- checking whether the fitted
beta stabilizes as N grows, and how the linearity-residual/branch-margin
guards' sensitivity changes with N. This was prohibitively expensive to
check under the old N-separate-sims architecture (each extra length was
another full sim); merging is what makes it checkable at all.

Reduced scale, matching step 4's own scale (n_modes=2, 1 wavelength, real
res=40/decay_tol=1e-11, real auto-picked spacing) -- not full N_MODES=14/
11-wavelength production scale. The domain here is 2x longer than step 4's
(N_max=6 vs N=3 lengths, same spacing unit), so expect roughly double
step 4's merged-path cost (~2h) -- cluster-gated, not a login-node run.
See batch/verification/run_convergence_monitor_count.batch.

Run:
    python /path/to/examples/verification/python/convergence_study_monitor_count.py
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
N_MODES = 2  # reduced from production's 14, matching step 4's scale
N_MAX = 6  # sweep N in {3, 4, 5, 6} -- doc's own example set

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
      f"decay_tol={solver.decay_tol}, N_MAX={N_MAX}\n")

# ---------------------------------------------------------------------------
# Same seed/spacing formula as section_mode_beta() and step 4's equivalence
# test -- spacing is a single unit shared by every N in the sweep, so all
# N_MAX lengths sit in one shared domain rather than being recomputed per N.
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
lengths_all = [k * spacing for k in range(1, N_MAX + 1)]

print(f"seed n_eff={seed.tolist()}, beta_max={beta_max:.4f} rad/um, "
      f"max_safe_spacing={max_safe_spacing:.4f} um, spacing={spacing:.4f} um")
print(f"lengths_all (N={N_MAX}): {[f'{L:.4f}' for L in lengths_all]} um")
print(f"domain length: {lengths_all[-1]:.4f} um "
      f"(step 4's merged domain was {3 * spacing:.4f} um, N=3)\n")

# ---------------------------------------------------------------------------
# ONE merged run at the full N_MAX -- every smaller N below is a free
# re-fit of a prefix of this same phase data, no extra simulation.
# ---------------------------------------------------------------------------
print(f"=== One solve_section_multi_length() call, N={N_MAX} monitors ===")
t0 = time.time()
T_multi = solver.solve_section_multi_length(
    mmi, wavelengths, mid, lengths_all, N_MODES, **iface,
)
elapsed = time.time() - t0
print(f"Elapsed: {elapsed:.1f} s ({elapsed/3600:.2f} h)\n")

phases_all = np.zeros((N_MAX, len(wavelengths), N_MODES))
for k in range(N_MAX):
    phases_all[k] = np.angle(np.diagonal(T_multi[k], axis1=1, axis2=2))

# ---------------------------------------------------------------------------
# Re-fit beta for every N in {3, ..., N_MAX} from prefixes of phases_all --
# free (pure arithmetic), the whole point of doing this with one merged run.
# ---------------------------------------------------------------------------
print("=== Fitted beta vs. N (free re-fits of one run's phase data) ===")
results = {}  # N -> (beta[wl,mode], resid[wl,mode], margin[wl,mode])
for N in range(3, N_MAX + 1):
    beta_N = np.zeros((len(wavelengths), N_MODES))
    resid_N = np.zeros((len(wavelengths), N_MODES))
    margin_N = np.zeros((len(wavelengths), N_MODES))
    for i in range(len(wavelengths)):
        for m in range(N_MODES):
            fit = beta_fit_from_phases(lengths_all[:N], phases_all[:N, i, m], seed_beta[i, m])
            beta_N[i, m] = fit.beta
            resid_N[i, m] = fit.resid
            margin_N[i, m] = fit.margin
    results[N] = (beta_N, resid_N, margin_N)
    for i in range(len(wavelengths)):
        for m in range(N_MODES):
            print(f"  N={N}  wl[{i}] mode[{m}]: beta={beta_N[i, m]:.8f} rad/um  "
                  f"resid={resid_N[i, m]:.3e} rad  margin={margin_N[i, m]:.4f}")

print()
print("=== Convergence: |beta(N) - beta(N_MAX)| per cell ===")
beta_final, resid_final, margin_final = results[N_MAX]
for N in range(3, N_MAX + 1):
    beta_N, resid_N, margin_N = results[N]
    diff = np.abs(beta_N - beta_final)
    print(f"  N={N}: max|beta(N) - beta({N_MAX})| = {float(diff.max()):.6e} rad/um  "
          f"| worst resid={float(resid_N.max()):.3e} rad "
          f"(guard fires above {_BETA_FIT_RESIDUAL_MAX})  "
          f"| worst margin={float(margin_N.max()):.4f} "
          f"(guard fires above {_BETA_BRANCH_MARGIN_RAISE})")

print()
print("=== Convergence: |beta(N) - beta(N-1)| per cell (successive-N step size) ===")
for N in range(4, N_MAX + 1):
    beta_N, _, _ = results[N]
    beta_Nm1, _, _ = results[N - 1]
    diff = np.abs(beta_N - beta_Nm1)
    print(f"  N={N} vs N={N-1}: max|delta beta| = {float(diff.max()):.6e} rad/um")

print()
print("=== SUMMARY ===")
print(f"Single merged run (N={N_MAX}): {elapsed:.1f} s ({elapsed/3600:.2f} h)")
print(f"beta at N=3 (today's production default): "
      f"{results[3][0].flatten().tolist()}")
print(f"beta at N={N_MAX}: {beta_final.flatten().tolist()}")
print(f"max|beta(3) - beta({N_MAX})| = "
      f"{float(np.abs(results[3][0] - beta_final).max()):.6e} rad/um "
      f"-- today's 3-monitor default vs. the most-monitors fit measured here.")
