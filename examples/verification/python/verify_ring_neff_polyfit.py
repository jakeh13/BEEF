"""Is the ring's degree-2 `n_eff` polynomial fit good enough?

`ring_resonator_validation.md` lists the interpolation strategy as load-bearing
and untested. `interp_complex` (the coefficient resampler) was closed on
2026-08-03 by `beef/spectral.py` + `tests/test_spectral.py`. This covers the
*other* half, which is a different function at different call sites:

    neff_fit = np.polyfit(wl_solve, neff_solve, deg=2)
    ... np.polyval(neff_fit, wl_dense)

The ring solves `n_eff` at a handful of wavelengths and fits a degree-2
polynomial through them, then evaluates that polynomial on the 101-point FDTD
grid and the 20001-point dense grid. A fit that is too low-order biases
`n_eff`, and `n_eff` enters the round-trip phase directly, so the error lands
on **resonance position** — the one quantity this device exists to report.

MEASURED IN THE CURRENCY THAT MATTERS
--------------------------------------
Reporting a residual in `n_eff` units invites the same mistake the resolution
sweep made before its lever arm was corrected: 1e-5 sounds small, but whether
it *is* small depends entirely on the path length multiplying it. So every
number here is converted the same way check 3 converts its own:

    dphi = 2*pi * dn_eff * L_OUTER / lambda        (L_OUTER = 7.0 um)

and reported as a fraction of 2*pi, i.e. of one FSR. That makes it directly
comparable to check 3's headline result — the shipped FDTD resolution moves the
comb by 12.5% of an FSR — instead of being a number in isolation.

NO NEW SIMULATION
-----------------
This reads `n_eff` values already in `example_sim_cache/` (26 bare-strip
wavelengths at res 100, 13 coupled-pair). The test is leave-one-out: drop an
interior point, fit the remaining ones, predict the dropped point, compare
against what the mode solver actually returned there. That is a genuine
held-out prediction, not a goodness-of-fit statistic on the training points —
a low-order fit can hug its own nodes and still be wrong between them, which is
exactly the failure mode in question.
"""
import glob
import json
import os
import pickle
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
CACHE = os.path.join(_ROOT, "example_sim_cache")
OUT = os.environ.get("RING_PF_OUT") or os.path.join(
    _ROOT, "example_output", "verification", "ring_neff_polyfit")

# == ring_resonator_hybrid.py:805 / _cascaded.py:864. The lever arm that turns
# an n_eff error into a phase error; see check 3's own lever-arm correction for
# what happens when this is taken on trust.
L_OUTER = 7.0
LAMBDA0 = 1.55
MODE_RES = 100.0                 # the ring's shipped MODE_RESOLUTION
SHIPPED_DEGREE = 2

# The bar: check 3 measured the shipped FDTD resolution as moving the comb by
# 12.5% of an FSR. An interpolation choice contributing a comparable amount
# would be a real defect; one contributing orders less is not worth changing.
CHECK3_FRACTION = 0.125


def _load(pattern, resolution, value_fn):
    rows = []
    for d in glob.glob(os.path.join(CACHE, pattern)):
        cfg_p = os.path.join(d, "config.json")
        res_p = os.path.join(d, "result.pkl")
        if not (os.path.exists(cfg_p) and os.path.exists(res_p)):
            continue
        cfg = json.load(open(cfg_p))
        if cfg.get("resolution") != resolution:
            continue
        rows.append((cfg["wavelength"], value_fn(pickle.load(open(res_p, "rb")))))
    rows.sort()
    if not rows:
        raise SystemExit(
            f"[pf] no cached entries matching {pattern!r} at resolution "
            f"{resolution} in {CACHE}. Run the ring once to populate it.")
    return np.array([w for w, _ in rows]), np.array([v for _, v in rows])


def dphi_fraction(dneff):
    """n_eff error -> round-trip phase error, as a fraction of one FSR."""
    return (2 * np.pi * dneff * L_OUTER / LAMBDA0) / (2 * np.pi)


def leave_one_out(wl, n, deg):
    """Max/mean held-out prediction error over the interior points.

    Endpoints are excluded deliberately: dropping one turns the remaining fit
    into an extrapolation, which is a different (and harsher) question than the
    interpolation the ring actually performs.
    """
    if len(wl) < deg + 3:
        return None
    errs = []
    for i in range(1, len(wl) - 1):
        m = np.ones(len(wl), bool)
        m[i] = False
        f = np.polyfit(wl[m], n[m], deg=deg)
        errs.append(abs(np.polyval(f, wl[i]) - n[i]))
    return np.array(errs)


def report(label, wl, n):
    print(f"\n[pf] ===== {label}: {len(wl)} cached wavelengths, "
          f"{wl.min():.3f}-{wl.max():.3f} um =====", flush=True)
    print(f"[pf] {'deg':>4} {'LOO max |dn|':>14} {'LOO mean |dn|':>14} "
          f"{'max dphi (rad)':>15} {'% of one FSR':>14}", flush=True)
    out = {}
    for deg in (1, 2, 3, 4):
        errs = leave_one_out(wl, n, deg)
        if errs is None:
            print(f"[pf] {deg:>4}   (too few points for a held-out fit)", flush=True)
            continue
        frac = dphi_fraction(errs.max())
        out[deg] = (errs.max(), errs.mean(), frac)
        mark = "  <- shipped" if deg == SHIPPED_DEGREE else ""
        print(f"[pf] {deg:>4} {errs.max():14.3e} {errs.mean():14.3e} "
              f"{2 * np.pi * errs.max() * L_OUTER / LAMBDA0:15.5f} "
              f"{frac * 100:13.4f}%{mark}", flush=True)
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"[pf] lever arm L_OUTER = {L_OUTER} um at lambda = {LAMBDA0} um",
          flush=True)
    print(f"[pf] bar for comparison: check 3 found the shipped FDTD resolution "
          f"moves the comb by {CHECK3_FRACTION * 100:.1f}% of an FSR", flush=True)

    wl_s, n_s = _load("*solve_neff*", MODE_RES, lambda r: float(r))
    strip = report("bare strip (solve_neff)", wl_s, n_s)

    # Supermodes come back as (neff_even, neff_odd); the even branch carries the
    # body phase that P_body's diagonal uses, so fit that one.
    try:
        wl_b, n_b = _load("*solve_body_supermodes*", MODE_RES,
                          lambda r: float(r[0]))
        body = report("coupled pair, even supermode", wl_b, n_b)
    except SystemExit as e:
        print(f"[pf] {e}", flush=True)
        wl_b, n_b, body = np.array([]), np.array([]), {}

    # --- the production configuration, end to end ------------------------
    # The ring fits 13 points and evaluates everywhere. Reproduce that: fit
    # degree 2 on a 13-point subset, then measure the residual against ALL
    # cached points, including the ones the fit never saw.
    wl13 = np.linspace(wl_s.min(), wl_s.max(), 13)
    idx = sorted({int(np.argmin(abs(wl_s - w))) for w in wl13})
    f2 = np.polyfit(wl_s[idx], n_s[idx], deg=SHIPPED_DEGREE)
    resid = n_s - np.polyval(f2, wl_s)
    prod_frac = dphi_fraction(abs(resid).max())
    print(f"\n[pf] PRODUCTION-LIKE: degree {SHIPPED_DEGREE} fit on "
          f"{len(idx)} points, residual against all {len(wl_s)}:", flush=True)
    print(f"[pf]   max |residual| = {abs(resid).max():.3e} in n_eff", flush=True)
    print(f"[pf]   -> {prod_frac * 100:.4f}% of one FSR", flush=True)

    # --- verdict ---------------------------------------------------------
    shipped = strip.get(SHIPPED_DEGREE)
    print("\n[pf] =========================== VERDICT ===========================",
          flush=True)
    if shipped is None:
        print("[pf] INCONCLUSIVE: not enough cached points to hold one out.",
              flush=True)
        return 1
    frac = max(shipped[2], prod_frac)
    ratio = CHECK3_FRACTION / frac if frac > 0 else float("inf")
    if frac < 0.01 * CHECK3_FRACTION:
        print(f"[pf] ✓ The degree-{SHIPPED_DEGREE} n_eff fit is NOT a limiting "
              f"error source.", flush=True)
        print(f"[pf]   Worst held-out phase error {frac * 100:.4f}% of an FSR, "
              f"about {ratio:.0f}x smaller than", flush=True)
        print(f"[pf]   the {CHECK3_FRACTION * 100:.1f}% the shipped FDTD "
              f"resolution contributes (check 3). Raising the", flush=True)
        print("[pf]   fit degree would change no resonance position anyone "
              "could measure.", flush=True)
        print(f"[pf]   Even degree 1 comes in at "
              f"{strip[1][2] * 100:.4f}% of an FSR, so degree "
              f"{SHIPPED_DEGREE} has real headroom.", flush=True)
    else:
        print(f"[pf] ✗ The degree-{SHIPPED_DEGREE} fit contributes "
              f"{frac * 100:.4f}% of an FSR — not negligible against check 3's "
              f"{CHECK3_FRACTION * 100:.1f}%.", flush=True)
        print("[pf]   Raise the degree or solve n_eff on a denser grid before "
              "trusting resonance positions.", flush=True)
    print("[pf] SCOPE: this validates the FIT, not the n_eff values it is fitted "
          "to. Those are", flush=True)
    print(f"[pf]   check 3's business, and at res {MODE_RES:g} they are "
          "converged.", flush=True)

    np.savez(os.path.join(OUT, "neff_polyfit.npz"),
             wl_strip=wl_s, neff_strip=n_s, wl_body=wl_b, neff_body=n_b,
             degrees=np.array(sorted(strip)),
             loo_max=np.array([strip[d][0] for d in sorted(strip)]),
             loo_frac_fsr=np.array([strip[d][2] for d in sorted(strip)]),
             production_residual=resid, production_frac_fsr=prod_frac,
             l_outer=L_OUTER)
    print(f"\n[pf] saved -> {os.path.join(OUT, 'neff_polyfit.npz')}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
