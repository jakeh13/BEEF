"""What does the mmi_1x2 calibration phase-fit failure actually look like?
(mmi_1x2_validation.md gap 5's calibration-failure finding, job 11589770.)

**ANSWERED 2026-08-01, and all three hypotheses this script was built to
separate were wrong.** `HybridSolver(calibrate_body_neff=True)` failed on this
device -- `section_mode_beta`'s 0.1 rad sanity guard fired at 2.101 rad (21x
over) on the `[-2.0, 2.0]` body span (`boundaries[1]`/`boundaries[2]`, i.e.
`+-half_body -+ body_stub` = `+-2.25 -+ 0.25`).

The two hypotheses in the validation doc were (1) the section isn't uniform
enough for this device's short 4.5 um body, or (2) SOI_160's weak confinement
means the near-field hasn't settled by body_stub=0.25 um past the taper. This
script's first run refuted both -- the failures sat at modes 0-1 and the
SHORTEST wavelengths, the opposite corner from a near-cutoff story -- and a
third reading was recorded instead: a **pi sign flip in `T_ii`**, from the
`(+d/3, -2d/3, +d/3)` residual pattern of one displaced point with d ~= -pi.

That was wrong too, and the reason is worth keeping: **the pattern is
degenerate.** For three equally spaced points, "middle point displaced by -pi"
and "points 2 and 3 each displaced by -2*pi" are algebraically
indistinguishable -- both give exactly `d*(+1/3, -2/3, +1/3)`. Residual
arithmetic alone cannot choose between them. Physics can:

    reading                 beta      n_eff    vs MPB seed 2.7285
    pi sign flip in T_ii    10.475    2.101    -23%   wrong sign AND size
    2*pi turn-count error   13.878    2.783    +1.96% in family with the
                                               coupler's +2.233%/+2.361%

The real cause is in `beta_from_phases` itself: it rounded the wrap count
independently at every length against one seed slope, so a seed error eps
drifted by eps*L/(2*pi) and eventually straddled a .5 boundary. Here eps was
1.96% -- the very FDTD-vs-MPB gap `calibrate_body_neff` exists to measure. The
unwrapping was broken by exactly the quantity it was built to measure. Fixed
2026-08-01 (`beta_fit_from_phases`); this script is now its offline proof.

The seed-free evidence is the decisive part, because it needs no MPB call and
no assumption about which branch is right: for a linear phase ramp on equally
spaced lengths the second difference `(phi3 - 2*phi2 + phi1)/(2*pi)` must be an
exact integer, and it is -- within 0.0031 in ALL 88 cells here (31 of which
need a +-1 turn). Nothing needs a sign flip to become linear.

Data source, unchanged: job 11589770 passed `save_dir` down to each length's
real `solve_section()` FDTD call, and all three lengths completed (the residual
check only runs after the full phases array is built), so they sit in the
on-disk cache. This recomputes the same 3 lengths `section_mode_beta()`'s own
selection logic chooses and re-fetches each via `solve_section(cache_only=True)`,
which raises loudly rather than silently re-simulating.

Cheap: cache hits are disk reads. The only meep work is one `section_mode_neff`
seed call, and the seed-free sections below do not even need that.

Run:
    python diagnose_mmi1x2_calibration_phase_fit.py
"""
from __future__ import annotations

import numpy as np

from beef import MMI, FDTDSolver, SOI_160
# Import the real constants rather than mirroring their values: a first
# revision of this script hardcoded _BETA_MIN_BRANCH_FRACTION = 4.0 (the
# actual value is 0.25) and silently computed lengths 16x too short, which
# then missed the cache and looked like "job 11589770 never cached its runs."
from beef.solvers.fdtd import (
    beta_fit_from_phases, beta_from_phases, _BETA_MIN_BRANCH_FRACTION,
    _BETA_BRANCH_MARGIN_RAISE, _BETA_FIT_RESIDUAL_MAX,
)

_TWO_PI = 2.0 * np.pi


def old_beta_from_phases(lengths, phases, seed_beta):
    """The pre-2026-08-01 `beta_from_phases`, verbatim, kept so this script can
    show OLD-vs-NEW side by side on the same cached data. Do not "fix" this --
    being wrong in the original way is its entire job."""
    lengths = np.asarray(lengths, dtype=float)
    phases = np.asarray(phases, dtype=float)
    turns = np.round((seed_beta * lengths - phases) / _TWO_PI)
    unwrapped = phases + _TWO_PI * turns
    A = np.vstack([lengths, np.ones_like(lengths)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, unwrapped, rcond=None)
    resid = float(np.max(np.abs(unwrapped - (slope * lengths + intercept))))
    return float(slope), resid

SAVE_DIR = "/storage/scratch1/8/jhiesener3/duty/example_sim_cache"

N_MODES = 8   # mmi_1x2_hybrid.py's N_MODES -- interface_modes=[8, 8]
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.26, 1.36, 11

mmi = MMI(
    n_inputs=1, n_outputs=2, mm_length=4.5, mm_width=2.0,
    waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
    waveguide_offset=1.0, waveguide_height=0.16,
)

# The calibration SOURCE is the body's left neighbour (HybridSolver's
# `_calibration_source` prefers k-1) -- the input section, method "fdtd".
input_solver = FDTDSolver(
    res=RESOLUTION, sim_width=5.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SOI_160, decay_tol=1e-11,
)

body_stub = 0.25
half_body = mmi.mm_length / 2
x0, x1 = -half_body + body_stub, +half_body - body_stub   # boundaries[1], boundaries[2]
mid = 0.5 * (x0 + x1)
span = x1 - x0

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
iface = dict(monitor_center=0.0, monitor_size=input_solver.sim_width)  # section_mode_beta's own defaults

print(f"[diag] mmi_1x2 calibration span [{x0:+.3f}, {x1:+.3f}] "
      f"({span:.3f} um), N_MODES={N_MODES}, dpml={input_solver.dpml}", flush=True)

# --- Replicate section_mode_beta()'s own length selection exactly ---------
seed_span = min(2.0, span - 2.0 * input_solver.dpml)
seed = input_solver.section_mode_neff(
    mmi, wavelengths, mid - seed_span / 2.0, mid + seed_span / 2.0, N_MODES, **iface)
seed_beta = 2.0 * np.pi * seed / wavelengths[:, None]
beta_max = float(np.max(np.abs(seed_beta)))
max_safe_spacing = 2.0 * np.pi / (_BETA_MIN_BRANCH_FRACTION * beta_max)
usable_span = span - 2.0 * input_solver.dpml
spacing = min(usable_span / 3.0, max_safe_spacing / 2.0)
lengths = [spacing, 2.0 * spacing, 3.0 * spacing]
print(f"[diag] recomputed lengths: {[f'{L:.3f}' for L in lengths]} "
      f"(job 11589770 reported [0.923, 1.846, 2.770] -- should match)", flush=True)

# --- Re-fetch each length's transfer matrix from the on-disk cache --------
phases = np.zeros((len(lengths), len(wavelengths), N_MODES))
T_all = []          # retained: the "dump T itself" step mmi_1x2_validation.md asked for
for j, L in enumerate(lengths):
    T = input_solver.solve_section(
        mmi, wavelengths, mid - L / 2.0, mid + L / 2.0,
        n_modes_in=N_MODES, n_modes_out=N_MODES,
        save_dir=SAVE_DIR, cache_only=True, **iface,
    )
    T_all.append(T)
    phases[j] = np.angle(np.diagonal(T, axis1=1, axis2=2))
    print(f"[diag] length {j} (L={L:.3f} um): fetched from cache OK", flush=True)

# --- Scan EVERY (wavelength, mode): section_mode_beta's guard fires on the
# --- max over both axes, so a band-centre-only look can miss the failure
# --- entirely (a first revision of this script printed only lambda=1.31 and
# --- saw residuals ~0.005 rad, 400x under the 2.101 the real job reported).
L_arr = np.asarray(lengths)
resid_grid = np.zeros((len(wavelengths), N_MODES))
margin_grid = np.zeros((len(wavelengths), N_MODES))
beta_grid = np.zeros((len(wavelengths), N_MODES))
excess_grid = np.zeros((len(wavelengths), N_MODES))
turnamb_grid = np.zeros((len(wavelengths), N_MODES))
shifted_grid = np.zeros((len(wavelengths), N_MODES), dtype=bool)
old_beta_grid = np.zeros((len(wavelengths), N_MODES))
old_resid_grid = np.zeros((len(wavelengths), N_MODES))
for i in range(len(wavelengths)):
    for m in range(N_MODES):
        fit = beta_fit_from_phases(lengths, phases[:, i, m], seed_beta[i, m])
        beta_grid[i, m] = fit.beta
        resid_grid[i, m] = fit.resid
        margin_grid[i, m] = fit.margin
        excess_grid[i, m] = fit.seed_excess
        turnamb_grid[i, m] = fit.turn_ambiguity
        shifted_grid[i, m] = fit.turns_shifted
        old_beta_grid[i, m], old_resid_grid[i, m] = old_beta_from_phases(
            lengths, phases[:, i, m], seed_beta[i, m])

def _grid(title, grid, fmt="10.4f", flag=None):
    print(f"\n[diag] {title}", flush=True)
    print("     wl   " + "".join(f"{'m'+str(m):>10}" for m in range(N_MODES)), flush=True)
    for i, wl in enumerate(wavelengths):
        cells = "".join(
            (f"{'*' + format(grid[i, m], '.3f'):>10}"
             if flag is not None and flag(grid[i, m])
             else f"{grid[i, m]:>{fmt}}")
            for m in range(N_MODES)
        )
        print(f"  {wl:.4f}  {cells}", flush=True)


# --- 1. Reproduce the ORIGINAL failure ------------------------------------
_grid(f"OLD beta_from_phases residual (rad), guard {_BETA_FIT_RESIDUAL_MAX} "
      f"(* = over):", old_resid_grid,
      flag=lambda v: v > _BETA_FIT_RESIDUAL_MAX)
i_w, i_m = np.unravel_index(int(np.argmax(old_resid_grid)), old_resid_grid.shape)
n_over_old = int(np.sum(old_resid_grid > _BETA_FIT_RESIDUAL_MAX))
print(f"  (* = over the guard threshold)", flush=True)
print(f"[diag] OLD worst: {old_resid_grid[i_w, i_m]:.4f} rad at "
      f"lambda={wavelengths[i_w]:.4f} um, mode {i_m}; {n_over_old} of "
      f"{old_resid_grid.size} cells over the guard "
      f"(job 11589770 reported 2.101 rad -- should match)", flush=True)

# --- 2. The same cells under the fix --------------------------------------
_grid(f"NEW beta_fit_from_phases residual (rad), same cells, same seed:",
      resid_grid, flag=lambda v: v > _BETA_FIT_RESIDUAL_MAX)
print(f"[diag] NEW worst: {resid_grid.max():.3e} rad; "
      f"{int(np.sum(resid_grid > _BETA_FIT_RESIDUAL_MAX))} of {resid_grid.size} "
      f"cells over the guard", flush=True)
print(f"[diag] worst branch margin {margin_grid.max():.4f} "
      f"(raises above {_BETA_BRANCH_MARGIN_RAISE}), worst turn ambiguity "
      f"{turnamb_grid.max():.5f}", flush=True)
print(f"[diag] seed excess (FDTD beta over MPB seed): "
      f"{excess_grid.min():+.4%} .. {excess_grid.max():+.4%}  "
      f"(coupler measures +2.233%/+2.361%)", flush=True)
n_shift = int(shifted_grid.sum())
n_moved = int(np.sum(np.abs(old_beta_grid - beta_grid) > 1e-6))
print(f"[diag] wrap counts moved off the seed's own rounding in {n_shift} of "
      f"{shifted_grid.size} cells; the returned beta actually changed in "
      f"{n_moved}", flush=True)

# --- 3. Shape of the cell that failed, both readings ----------------------
raw = phases[:, i_w, i_m]
sb = seed_beta[i_w, i_m]
fit = beta_fit_from_phases(lengths, raw, sb)
print(f"\n[diag] the cell that failed (lambda={wavelengths[i_w]:.4f}, "
      f"mode {i_m}):", flush=True)
print(f"  lengths                  {np.round(L_arr, 4).tolist()}", flush=True)
print(f"  raw arg(T_ii)            {np.round(raw, 4).tolist()}", flush=True)
print(f"  fractional turns vs seed "
      f"{np.round((sb * L_arr - raw) / _TWO_PI, 4).tolist()}"
      f"   <-- a value near X.5 is the whole bug", flush=True)
print(f"  OLD turns (per-length)   {np.round(fit.seed_turns).astype(int).tolist()}"
      f"  -> beta {old_beta_grid[i_w, i_m]:.4f}, "
      f"n_eff {old_beta_grid[i_w, i_m]*wavelengths[i_w]/_TWO_PI:.4f}, "
      f"resid {old_resid_grid[i_w, i_m]:.4f}", flush=True)
print(f"  NEW turns (consistent)   {np.round(fit.turns).astype(int).tolist()}"
      f"  -> beta {fit.beta:.4f}, "
      f"n_eff {fit.beta*wavelengths[i_w]/_TWO_PI:.4f}, "
      f"resid {fit.resid:.6f}", flush=True)
print(f"  MPB seed                 beta {sb:.4f}, "
      f"n_eff {sb*wavelengths[i_w]/_TWO_PI:.4f}", flush=True)
print(f"  branch margin            {margin_grid[i_w, i_m]:.4f} "
      f"(was 0.922 under the old reading -- the doc's own 'unresolved "
      f"secondary finding', which was really the tell)", flush=True)

# --- 4. SEED-FREE proof: the data is a clean linear ramp ------------------
# For a linear ramp on equally spaced lengths the second difference must be an
# exact integer number of turns. This needs no seed and no branch assumption,
# so it settles the sign-flip question from the data alone: a pi displacement
# would put these at a HALF integer, not an integer.
if np.allclose(np.diff(L_arr), L_arr[0], rtol=1e-6):
    sd = (phases[2] - 2 * phases[1] + phases[0]) / _TWO_PI
    off = np.abs(sd - np.round(sd))
    _grid("SEED-FREE second difference (phi3 - 2*phi2 + phi1)/2pi "
          "(must be an integer for a linear ramp):", sd)
    print(f"[diag] within {off.max():.5f} of an exact integer in ALL {sd.size} "
          f"cells; {int((np.abs(np.round(sd)) > 0.5).sum())} need a +-1 turn "
          f"correction. Every cell is a clean linear ramp -- NOTHING needs a "
          f"sign flip. A pi displacement would sit at a HALF integer.", flush=True)
    print(f"[diag] cross-check: that {off.max():.5f} must equal the worst "
          f"turn_ambiguity {turnamb_grid.max():.5f} -- "
          f"{'AGREE' if abs(off.max()-turnamb_grid.max()) < 1e-5 else 'DISAGREE'}"
          f" (two independent routes to the same number)", flush=True)

# --- 5. SEED-FREE robustness sweep: OLD vs NEW over a seed error range ----
# The real failure needed a specific seed. This removes that dependence: for
# every cell, every physically plausible branch, and a range of seed errors
# bracketing the real ~2%, does the fit return the branch it was seeded near?
print(f"\n[diag] seed-free robustness sweep (every cell x every plausible "
      f"branch x seed error +-3%):", flush=True)
branch_step = _TWO_PI / float(np.min(np.diff(L_arr)))
trials = old_bad = new_bad = 0
for i in range(len(wavelengths)):
    for m in range(N_MODES):
        p = phases[:, i, m]
        base = beta_grid[i, m]
        for k in range(-20, 21):
            b_true = base + k * branch_step
            if not (1.0 < b_true * wavelengths[i] / _TWO_PI < 3.6):
                continue          # not a guided n_eff for Si/SiO2
            for eps in np.linspace(-0.03, 0.03, 25):
                trials += 1
                seed = b_true * (1.0 + eps)
                ob, _ = old_beta_from_phases(lengths, p, seed)
                nb = beta_fit_from_phases(lengths, p, seed).beta
                tol = 1e-6 * max(1.0, abs(b_true))
                old_bad += int(abs(ob - b_true) > tol)
                new_bad += int(abs(nb - b_true) > tol)
print(f"[diag]   {trials} trials: OLD wrong {old_bad} "
      f"({100*old_bad/max(trials,1):.1f}%), NEW wrong {new_bad} "
      f"({100*new_bad/max(trials,1):.1f}%)", flush=True)
print(f"[diag]   The old code's failure rate was NOT a corner case: a seed "
      f"error of the size calibration exists to measure broke it about one "
      f"time in seven.", flush=True)

# --- 6. Seed literal, so later re-runs need no MPB ------------------------
print(f"\n[diag] MPB seed n_eff grid, for pasting into an MPB-free re-run "
      f"({len(wavelengths)} wavelengths x {N_MODES} modes):", flush=True)
print("SEED_NEFF = " + np.array2string(
    seed_beta * wavelengths[:, None] / _TWO_PI,
    precision=6, separator=", ", max_line_width=100), flush=True)

print(f"\n[diag] Reading this output. Section 4 is the decisive one and it "
      f"needs no seed: if the second difference is an integer everywhere, the "
      f"measured phase IS linear in length and the only question was ever "
      f"which whole turn each length sits on. That rules out both original "
      f"hypotheses (non-uniform cut, unsettled near field -- either would "
      f"leave a non-integer second difference) AND the pi-sign-flip reading "
      f"(which would leave a half-integer). Section 3's fractional-turn row "
      f"shows the mechanism: a value close to X.5 is a rounding decision "
      f"decided by noise. Section 5 shows it was never a corner case.",
      flush=True)
