"""Step 6 of validation_markdown/resolved/section_mode_beta_multi_monitor_validation.md.

Convergence study #2 -- total simulation/domain length, INDEPENDENT of
monitor count. Step 5 (job 11838953) swept how many of ONE run's 6
already-measured monitors went into the beta fit, but that run's domain
was fixed at 6x spacing throughout -- it never actually varied domain
length. This script does: hold the 3 REAL fit lengths fixed (the same
[spacing, 2*spacing, 3*spacing] production/step4/step5 lengths), and vary
ONLY the total domain size by adding one extra, ignored "padding" length
that forces solve_section_multi_length()'s L_domain = max(lengths) bigger
without being used in the beta fit. Compares whether the 3 real lengths'
own fitted phase/beta shift as the shared domain around them grows --
exactly the near-field/launch-transient contamination question step 2
flagged and step 5's (corrected) writeup deferred to this step.

Reuses two domain sizes already measured for free, zero extra compute:
  - domain = 3*spacing (3.3069 um): step 4's actual merged-path run,
    job 11813378 -- the 3 real lengths with NO padding monitor at all.
  - domain = 6*spacing (6.6138 um): step 5's own N=3 row, job 11838953 --
    same 3 real lengths, but 3 OTHER monitors also present (unused in that
    row's fit), forcing a domain built for the full 6-length set.
This script adds two NEW domain sizes (pad = 4*spacing, 8*spacing) via two
genuinely separate solve_section_multi_length() calls (each a real domain
build + decay wait -- NOT free, unlike step 5's single-run reuse), and
prints a combined 4-point table (3x, 4x, 6x, 8x spacing) with all values in
one place for the doc writeup.

Reduced scale, matching steps 4-5 (n_modes=2, 1 wavelength, real
res=40/decay_tol=1e-11) -- not full N_MODES=14/11-wavelength production
scale. Cost projected from steps 4-5's own two data points (linear fit
time(L) ~= 3227 + 1243*L s): pad=4 ~2.4h, pad=8 ~3.9h, ~6.3h combined --
cluster-gated. See batch/verification/run_convergence_domain_length.batch.

Run:
    python /path/to/examples/verification/python/convergence_study_domain_length.py
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
N_MODES = 2  # matches steps 4-5's reduced scale
PAD_MULTIPLIERS = [4, 8]  # NEW domain sizes (in units of spacing); 3x and 6x are free/reused

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
# Same seed/spacing formula as steps 4-5 -- the 3 REAL fit lengths never
# change across this whole script; only the (ignored) padding length does.
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
lengths_real = [spacing, 2.0 * spacing, 3.0 * spacing]

print(f"seed n_eff={seed.tolist()}, beta_max={beta_max:.4f} rad/um, spacing={spacing:.4f} um")
print(f"lengths_real (always fit, never change): {[f'{L:.4f}' for L in lengths_real]} um\n")

# ---------------------------------------------------------------------------
# Free reference points -- already measured, zero new compute. Hardcoded
# from the recorded job logs (11813378, 11838953), not recomputed.
# ---------------------------------------------------------------------------
free_points = {
    3: dict(
        domain_L=3.0 * spacing, job="11813378 (step 4 merged path, no padding)",
        beta=[11.571281, 11.407538], resid=[8.00e-05, 7.66e-05], margin=[0.060, 0.060],
        elapsed_s=7339.3,
    ),
    6: dict(
        domain_L=6.0 * spacing, job="11838953 (step 5 N=3 row, 3 extra unused monitors present)",
        beta=[11.571207, 11.407468], resid=[1.069e-04, 1.024e-04], margin=[0.0601, 0.0603],
        elapsed_s=11451.1,
    ),
}

# ---------------------------------------------------------------------------
# New domain sizes: 3 real lengths + 1 ignored padding length, each its own
# genuine solve_section_multi_length() call (real domain build + decay wait).
# ---------------------------------------------------------------------------
new_points = {}
for pad_mult in PAD_MULTIPLIERS:
    pad_length = pad_mult * spacing
    lengths_call = lengths_real + [pad_length]
    print(f"=== New domain: pad={pad_mult}x spacing = {pad_length:.4f} um "
          f"(domain L_domain={pad_length:.4f} um) ===")
    t0 = time.time()
    T = solver.solve_section_multi_length(
        mmi, wavelengths, mid, lengths_call, N_MODES, **iface,
    )
    elapsed = time.time() - t0
    print(f"Elapsed: {elapsed:.1f} s ({elapsed/3600:.2f} h)")

    phases_real = np.zeros((3, len(wavelengths), N_MODES))
    for k in range(3):
        phases_real[k] = np.angle(np.diagonal(T[k], axis1=1, axis2=2))
    # T[3] (the padding monitor) is deliberately never read -- it exists
    # only to force L_domain bigger, not to contribute to the fit.

    beta_arr = np.zeros((len(wavelengths), N_MODES))
    resid_arr = np.zeros((len(wavelengths), N_MODES))
    margin_arr = np.zeros((len(wavelengths), N_MODES))
    for i in range(len(wavelengths)):
        for m in range(N_MODES):
            fit = beta_fit_from_phases(lengths_real, phases_real[:, i, m], seed_beta[i, m])
            beta_arr[i, m] = fit.beta
            resid_arr[i, m] = fit.resid
            margin_arr[i, m] = fit.margin
            print(f"  wl[{i}] mode[{m}]: beta={fit.beta:.6f} rad/um  "
                  f"resid={fit.resid:.3e} rad  margin={fit.margin:.4f}")
    new_points[pad_mult] = dict(
        domain_L=pad_length, beta=beta_arr.flatten().tolist(),
        resid=resid_arr.flatten().tolist(), margin=margin_arr.flatten().tolist(),
        elapsed_s=elapsed,
    )
    print()

# ---------------------------------------------------------------------------
# Combined table: 4 domain sizes (3x free, 4x new, 6x free, 8x new), same
# 3 real fit lengths throughout.
# ---------------------------------------------------------------------------
print("=== SUMMARY: same 3 real lengths, domain length swept independently ===")
all_domains = sorted(list(free_points.keys()) + list(new_points.keys()))
print(f"{'mult':>5} {'domain_L(um)':>13} {'beta_mode0':>12} {'beta_mode1':>12} "
      f"{'worst_resid':>12} {'worst_margin':>13} {'source':>10}")
for mult in all_domains:
    if mult in free_points:
        p = free_points[mult]
        source = "free/reused"
    else:
        p = new_points[mult]
        source = "new"
    beta0, beta1 = p["beta"][0], p["beta"][1]
    worst_resid = max(p["resid"]) if isinstance(p["resid"], list) else p["resid"]
    worst_margin = max(p["margin"]) if isinstance(p["margin"], list) else p["margin"]
    print(f"{mult:>5} {p['domain_L']:>13.4f} {beta0:>12.6f} {beta1:>12.6f} "
          f"{worst_resid:>12.3e} {worst_margin:>13.4f} {source:>10}")

beta0_all = [
    (free_points[mult]["beta"][0] if mult in free_points else new_points[mult]["beta"][0])
    for mult in all_domains
]
beta1_all = [
    (free_points[mult]["beta"][1] if mult in free_points else new_points[mult]["beta"][1])
    for mult in all_domains
]
print(f"\nbeta mode0 range across domain sweep: {max(beta0_all) - min(beta0_all):.6e} rad/um")
print(f"beta mode1 range across domain sweep: {max(beta1_all) - min(beta1_all):.6e} rad/um")
print(f"(guard fires above resid={_BETA_FIT_RESIDUAL_MAX} rad, margin={_BETA_BRANCH_MARGIN_RAISE})")
