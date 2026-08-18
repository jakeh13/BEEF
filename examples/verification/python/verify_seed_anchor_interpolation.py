"""Verify the anchored MPB seed against the real 11-wavelength MPB seed.

resolved/section_mode_beta_multi_monitor_validation.md step 9: the seed
(`FDTDSolver.section_mode_neff`) is now solved at `seed_anchors=3` wavelengths
and interpolated to the caller's grid, because its two jobs -- picking each
measured phase's 2*pi branch, and setting the fit-length spacing via beta_max --
are coarse decisions. The error budget for that default was sized against the
2x2 body's *calibrated* n_eff curve (a stand-in, since the MPB seed itself was
never cached). This measures the real thing.

Independent of job 11937894 by construction: it touches only
FDTDSolver_section_mode_neff_* cache entries and writes to its own output
folder, so it can run alongside the merged calibration.

Cost: 154 solves for the 11-point truth + 42 for the production 3-anchor set
= 196 eigensolves, ~40 min at 48 ranks. The 5- and 7-anchor convergence points
are on-grid subsets of the truth, so they cost nothing extra.

Run:
    srun --mpi=pmi2 python3 examples/verification/python/verify_seed_anchor_interpolation.py
"""

import json
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, FDTDSolver, SOI_220
from beef.solvers.fdtd import _BETA_MIN_BRANCH_FRACTION

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
CACHE = os.path.join(_ROOT, "example_sim_cache")
OUT = os.path.join(_ROOT, "example_output", "verification", "seed_anchor_interpolation")
os.makedirs(OUT, exist_ok=True)

N_MODES = 14
wl = np.linspace(1.530, 1.565, 11)

mmi = MMI(n_inputs=2, n_outputs=2, mm_length=28.0, mm_width=2.8,
          waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
          waveguide_offset=1.5, waveguide_height=0.22)
solver = FDTDSolver(res=40.0, sim_width=5.0, sim_height=1.5, dpml=0.5,
                    source_margin=0.0, monitor_margin=0.25, three_d=True,
                    material_stack=SOI_220, decay_tol=1e-11,
                    half_domain_aperture=True)

# exactly the window section_mode_beta_details() seeds on for this device:
# seed_span = min(2.0, 27.5 - 2*dpml) = 2.0 um centred on the body midpoint
X0, X1 = -1.0, 1.0
iface = dict(monitor_center=0.0, monitor_size=solver.sim_width)

print("=" * 76)
print("Anchored MPB seed vs the real 11-wavelength seed — MMI 2x2 body")
print("=" * 76)
print(f"window [{X0}, {X1}] um, {N_MODES} modes, {len(wl)} wavelengths "
      f"{wl[0]}-{wl[-1]} um, res={solver.res}\n")

t0 = time.time()
truth = solver.section_mode_neff(mmi, wl, X0, X1, N_MODES, save_dir=CACHE, **iface)
t_truth = time.time() - t0
print(f"11-point truth : {len(wl) * N_MODES} eigensolves, {t_truth/60:.1f} min")

anchors = np.linspace(float(wl.min()), float(wl.max()), 3)
t0 = time.time()
seed3 = solver.section_mode_neff(mmi, anchors, X0, X1, N_MODES, save_dir=CACHE, **iface)
t_anchor = time.time() - t0
print(f"3-anchor seed  : {len(anchors) * N_MODES} eigensolves, {t_anchor/60:.1f} min"
      f"   -> {t_truth/max(t_anchor, 1e-9):.1f}x cheaper\n")


def interp(wl_seed, neff_seed):
    """Exactly section_mode_beta_details()'s own fit."""
    x_mid = float(np.mean(wl_seed))
    out = np.empty((len(wl), neff_seed.shape[1]))
    for m in range(neff_seed.shape[1]):
        co = np.polyfit(wl_seed - x_mid, neff_seed[:, m], len(wl_seed) - 1)
        out[:, m] = np.polyval(co, wl - x_mid)
    return out


# sanity: the anchors must agree with the truth at the wavelengths they share
shared = [int(np.argmin(np.abs(wl - a))) for a in anchors]
agree = np.max(np.abs(seed3 - truth[shared]))
print(f"anchors vs truth at the 3 shared wavelengths: max|d| = {agree:.2e}"
      f"   (must be ~0 — same solver, same grid)\n")

fit3 = interp(anchors, seed3)
err3 = fit3 - truth

half_branch = np.pi / 1.085983652     # this device's own fit spacing
print("=== 3-anchor interpolation error vs the real seed, per mode ===")
print(f"  {'mode':>4}  {'n_eff@1.530':>11}  {'max|dn_eff|':>12}  {'d(beta) rad/um':>15}  "
      f"{'/branch half':>13}")
for m in range(N_MODES):
    e = float(np.abs(err3[:, m]).max())
    db = 2 * np.pi * e / 1.5475
    print(f"  {m:>4}  {truth[0, m]:>11.6f}  {e:>12.2e}  {db:>15.2e}  {db/half_branch:>13.2e}")
emax = float(np.abs(err3).max())
print(f"\n  worst overall: {emax:.2e} in n_eff = {2*np.pi*emax/1.5475:.2e} rad/um"
      f" = {2*np.pi*emax/1.5475/half_branch:.2e} of the branch half-spacing")
print(f"  production branch margins for comparison: 0.05 - 0.095\n")

# --- does beta_max move? it sets the fit lengths, and hence every downstream
#     cache key. The anchors include both band edges, so the maximum of
#     2*pi*n/lambda sits ON an anchor and should be reproduced exactly.
def lengths_from(seed_neff):
    sb = 2.0 * np.pi * seed_neff / wl[:, None]
    bmax = float(np.max(np.abs(sb)))
    max_safe = 2.0 * np.pi / (_BETA_MIN_BRANCH_FRACTION * bmax)
    spacing = min((27.5 - 2 * solver.dpml) / 3.0, max_safe / 2.0)
    return bmax, [spacing, 2 * spacing, 3 * spacing]

b_t, L_t = lengths_from(truth)
b_i, L_i = lengths_from(fit3)
print("=== the consequence that matters for caching: do the fit lengths move? ===")
print(f"  beta_max  truth {b_t:.9f}   interpolated {b_i:.9f}   d = {abs(b_t-b_i):.2e}")
print(f"  lengths   truth {[round(x, 9) for x in L_t]}")
print(f"            interp{[round(x, 9) for x in L_i]}")
same = all(round(a, 9) == round(b, 9) for a, b in zip(L_t, L_i))
print(f"  identical to the 9 dp the cache rounds at: {same}"
      f"   <- {'downstream cache entries stay valid' if same else 'DOWNSTREAM CACHE INVALIDATED'}\n")

# --- effect on seed_excess, the one reported quantity that shifts. beta_FDTD
#     from job 11872191's own calibrated neff_override, read off the cache.
ent = os.path.join(CACHE, "EMESolver_solve_section_length_sweep_MMI_4d8b56cdf51c5cd1",
                   "config.json")
if os.path.isfile(ent):
    n_fdtd = np.array(json.load(open(ent))["call_kwargs"]["neff_override"])[..., 0]
    beta_fdtd = 2 * np.pi * n_fdtd / wl[:, None]
    ex_t = beta_fdtd / (2 * np.pi * truth / wl[:, None]) - 1.0
    ex_i = beta_fdtd / (2 * np.pi * fit3 / wl[:, None]) - 1.0
    d_pp = np.abs(ex_t - ex_i) * 100.0
    print("=== effect on the reported seed_excess (beta/beta_seed - 1) ===")
    print(f"  excess range, true seed : {100*ex_t.min():.3f}% .. {100*ex_t.max():.3f}%")
    print(f"  worst shift from interpolating: {d_pp.max():.4f} percentage points"
          f" (mode {int(np.unravel_index(d_pp.argmax(), d_pp.shape)[1])})\n")
else:
    print("(calibrated neff entry absent — skipping the seed_excess comparison)\n")

# --- convergence in the number of anchors, free: on-grid subsets of the truth
print("=== convergence in seed_anchors (on-grid subsets, no extra solves) ===")
subsets = {3: [0, 5, 10], 5: [0, 3, 5, 7, 10], 7: [0, 2, 4, 5, 6, 8, 10]}
conv = {}
for n, idx in subsets.items():
    e = float(np.abs(interp(wl[idx], truth[idx]) - truth).max())
    conv[n] = e
    print(f"  {n} anchors (grid idx {idx}): worst |dn_eff| = {e:.2e}"
          f"  = {2*np.pi*e/1.5475/half_branch:.2e} of branch half")
print()

np.savez(os.path.join(OUT, "seed_anchor_verification.npz"),
         wavelengths=wl, anchors=anchors, truth=truth, seed3=seed3, fit3=fit3,
         err3=err3, beta_max_truth=b_t, beta_max_interp=b_i,
         lengths_truth=L_t, lengths_interp=L_i,
         t_truth_s=t_truth, t_anchor_s=t_anchor,
         conv_anchors=np.array(sorted(conv)), conv_err=np.array([conv[k] for k in sorted(conv)]))

fig, axes = plt.subplots(1, 3, figsize=(16, 4.4))
colors = plt.cm.viridis(np.linspace(0, 1, N_MODES))
ax = axes[0]
for m in range(N_MODES):
    ax.plot(wl, truth[:, m], color=colors[m], lw=1.0)
ax.plot(anchors, seed3[:, 0], "ko", ms=5, label="anchors solved")
for m in range(1, N_MODES):
    ax.plot(anchors, seed3[:, m], "ko", ms=3)
ax.set_xlabel("Wavelength (µm)"); ax.set_ylabel("$n_{eff}$ (MPB seed)")
ax.set_title("The seed: 14 bands, 3 solved wavelengths"); ax.legend(fontsize=8); ax.grid(alpha=.4)

ax = axes[1]
for m in range(N_MODES):
    ax.plot(wl, err3[:, m], color=colors[m], marker="o", ms=3, lw=1.0)
ax.axhline(0, color="gray", lw=.8)
ax.set_xlabel("Wavelength (µm)"); ax.set_ylabel(r"$n_{eff}^{interp} - n_{eff}^{true}$")
ax.set_title(f"Interpolation error (worst {emax:.1e})"); ax.grid(alpha=.4)

ax = axes[2]
ks = sorted(conv)
ax.semilogy(ks, [conv[k] for k in ks], "o-", color="seagreen")
ax.set_xlabel("seed_anchors"); ax.set_ylabel(r"worst $|\Delta n_{eff}|$")
ax.set_title("Convergence in anchor count"); ax.grid(alpha=.4, which="both")
ax.set_xticks(ks)

plt.suptitle("MPB seed — 3-anchor interpolation vs the real 11-wavelength solve "
             "(MMI 2×2 body, 14 modes)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(OUT, "seed_anchor_verification.png"), dpi=300)
plt.close("all")

print("=== SUMMARY ===")
print(f"anchors agree with truth where they overlap : {agree:.2e}")
print(f"worst interpolation error                   : {emax:.2e} in n_eff")
print(f"  as a fraction of the branch half-spacing  : {2*np.pi*emax/1.5475/half_branch:.2e}")
print(f"fit lengths unchanged at 9 dp               : {same}")
print(f"seed cost                                   : {t_truth/60:.1f} min -> {t_anchor/60:.1f} min")
print(f"artifacts -> {OUT}/")
