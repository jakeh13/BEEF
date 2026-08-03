"""Check 4, the comparison: full-device FDTD vs both hybrid routes.

`ring_resonator_fdtd.py` (job 11655143) solves the whole device in one meep run
— no decomposition, no EME, no network algebra, no reference planes. The two
hybrid scripts assemble the same device from parts and **disagree with each
other by ~0.05 FSR**. Nothing internal to those two can say which is right. This
script is where the direct solve adjudicates.

WHAT IT IS REALLY TESTING, GIVEN THE RES-60 RESULT
---------------------------------------------------
Job 11653509 swept `FDTD_RESOLUTION` 20/30/40/60 on the cascaded route with the
mode side pinned (`neff_solve[0]` bit-identical at all four, so everything that
moved was the FDTD half). Two quantities converged very differently:

    FSR                 fitted order p ~ 1.4-2.1   res 60 within 0.16% of an FSR
    absolute comb pos.  fitted order p ~ 0.70      res 60 still 29% of an FSR out

A *difference* of positions converging at ~h^2 while the *absolute* position
crawls at h^0.7 has one natural explanation: a **wavelength-flat phase offset**,
which shifts every dip equally and cancels exactly in the FSR. That is the same
signature as the unexplained flat 0.976*pi in the ring-excited coefficients, so
those are likely one defect and not two.

The full-device run is the discriminator, because it shares none of the hybrid's
phase bookkeeping:

  * If the full device's comb ALSO drifts at p ~ 0.7, the flat phase is physical
    /discretization — meep resolving the staircased arcs — and the hybrid is
    faithfully reproducing it. Resolution is then the whole story and no amount
    of bookkeeping work helps.
  * If the full device's comb converges CLEANLY while the hybrid's does not, the
    flat phase is in the **decomposition's reference-plane arithmetic**. That
    localizes 0.976*pi to code rather than physics, after two mechanism
    hypotheses were refuted by direct test.

Either outcome is informative, which is the point of running it.

WHY THE PARTNER RUNS USE RING_THREE_D=0 AND NOT RING_TEST_MODE=1
-----------------------------------------------------------------
The full device is only affordable in 2-D (3-D is ~167 h at res 40), so the
hybrid side must be 2-D too. `RING_TEST_MODE=1` would have been the obvious way
and is the wrong one: it moves FDTD resolution, `MODE_RESOLUTION`, `N_WL_FDTD`,
`N_WL_SOLVE` and `N_WL_DENSE` all at once, and `N_WL_DENSE=2001` is already
known to be grid-limited for Q. A disagreement measured against it could be a
decomposition error or any of four confounds. `RING_THREE_D=0` (split out of
`TEST_MODE` on 2026-08-03) moves exactly one knob.

Reads only saved `.npz`; runs no simulation.

Env:
    RING_C4_RES   comma-separated resolutions to compare (default 20,30,40)
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
_EO = os.path.join(_ROOT, "example_output", "ring_resonator")

RESOLUTIONS = [float(v) for v in os.environ.get("RING_C4_RES", "20,30,40").split(",")]
OUT = os.environ.get("RING_C4_OUT") or os.path.join(
    _ROOT, "example_output", "verification", "ring_check4")

FULL = os.path.join(_EO, "ring_resonator_fdtd", "2d", "res_{res:g}", "full_fdtd.npz")
HYB = {
    "monolithic": os.path.join(_EO, "ring_resonator_hybrid",
                               "res_sweep", "res_{res:g}", "two_d", "ring_data.npz"),
    "cascaded": os.path.join(_EO, "ring_resonator_hybrid_cascaded",
                             "res_sweep", "res_{res:g}", "two_d", "ring_data.npz"),
}


def _dips(path):
    """Comb wavelengths (um) from either script's npz.

    The hybrid scripts store (wavelength, T_min) pairs; the full-device script
    stores bare wavelengths. Normalise rather than assume.
    """
    d = np.load(path, allow_pickle=True)
    raw = d["dips"]
    return np.array([x[0] if np.ndim(x) else x for x in raw], dtype=float)


def _full_comb(path):
    """Full-device comb from GROUP DELAY, recomputed from the stored complex
    amplitudes rather than read from `dips`.

    Recomputed on purpose: `dips` in older full_fdtd.npz files holds AMPLITUDE
    minima, which for an all-pass ring measure loss rather than resonance (see
    ring_resonator_fdtd.py section 4). Deriving the comb here from
    `a_device`/`a_busonly` means this comparison gives the same answer against
    old and new npz files alike, instead of silently comparing a group-delay
    comb on one side against an amplitude comb on the other.
    """
    from scipy.signal import find_peaks
    d = np.load(path, allow_pickle=True)
    wl = d["wavelengths"]
    ratio = np.asarray(d["a_device"]) / np.asarray(d["a_busonly"])
    c = 299.792458                                   # um/ps
    tau_g = (wl ** 2 / (2 * np.pi * c)) * np.gradient(np.unwrap(np.angle(ratio)), wl)
    x = np.abs(tau_g - np.median(tau_g))
    if np.std(x) == 0:
        return np.array([])
    # One peak per resonance: without a minimum separation the shoulders of
    # each resonance are picked up too (measured: 22 "peaks" for 7 resonances
    # at res 40, in triplets 1.0-1.2 nm apart).
    fsr_guess = (1.55 ** 2) / (4.3 * 43.4159)
    dist = max(1, int(0.5 * fsr_guess / (wl[1] - wl[0])))
    pk, _ = find_peaks(x, prominence=np.std(x), distance=dist)
    return wl[pk]


def align(a, b):
    """Match two combs that may have a different number of lines.

    A comb near a window edge gains or loses a line under a small shift -- that
    is exactly what res 60 did in job 11653509, where the first dip slid below
    1.500 um and the naive mean FSR jumped by 0.13 nm for purely bookkeeping
    reasons. So never compare `mean(diff(a))` against `mean(diff(b))` directly.

    Returns (a_sub, b_sub) of equal length, choosing the offset that minimises
    total squared distance.
    """
    if len(a) == 0 or len(b) == 0:
        return a[:0], b[:0]
    n = min(len(a), len(b))
    best, best_cost = None, np.inf
    for ia in range(len(a) - n + 1):
        for ib in range(len(b) - n + 1):
            cost = float(np.sum((a[ia:ia + n] - b[ib:ib + n]) ** 2))
            if cost < best_cost:
                best_cost, best = cost, (a[ia:ia + n], b[ib:ib + n])
    return best


def fsr_of(c):
    return (c[-1] - c[0]) / (len(c) - 1) * 1000 if len(c) >= 2 else float("nan")


def main():
    os.makedirs(OUT, exist_ok=True)
    print("[c4] full-device FDTD vs hybrid routes, 2-D, matched physics", flush=True)
    print("[c4] (RING_THREE_D=0 on the hybrid side: one knob, not five)", flush=True)

    rows, missing = [], []
    for res in RESOLUTIONS:
        fp = FULL.format(res=res)
        if not os.path.exists(fp):
            missing.append(fp)
            continue
        full = _full_comb(fp)
        print(f"\n[c4] ================= res {res:g} =================", flush=True)
        print(f"[c4] full device : {len(full)} dips, FSR {fsr_of(full):.3f} nm",
              flush=True)
        for name, tmpl in HYB.items():
            hp = tmpl.format(res=res)
            if not os.path.exists(hp):
                missing.append(hp)
                print(f"[c4] {name:11s}: MISSING ({hp})", flush=True)
                continue
            hyb = _dips(hp)
            fa, fb = align(full, hyb)
            if len(fa) < 2:
                print(f"[c4] {name:11s}: too few dips to compare", flush=True)
                continue
            shift = float(np.mean(fb - fa)) * 1000            # nm
            # An offset approaching half an FSR is ambiguous: the aligner would
            # rather pair each line with its NEIGHBOUR than admit a big shift.
            # Say so instead of reporting a confidently wrong small number.
            if abs(shift) > 0.4 * fsr_of(fa):
                print(f"[c4] {name:11s}: ⚠ offset is {abs(shift) / fsr_of(fa) * 100:.0f}% "
                      f"of an FSR — approaching the half-FSR aliasing limit, so "
                      f"the line pairing itself may be off by one.", flush=True)
            fsr_full, fsr_hyb = fsr_of(fa), fsr_of(fb)
            frac = abs(shift) / fsr_full * 100 if fsr_full else float("nan")
            print(f"[c4] {name:11s}: {len(hyb)} dips, FSR {fsr_of(hyb):.3f} nm  "
                  f"| aligned on {len(fa)} lines", flush=True)
            print(f"[c4]              comb offset vs full device: {shift:+7.3f} nm "
                  f"= {frac:5.2f}% of an FSR", flush=True)
            print(f"[c4]              FSR difference            : "
                  f"{fsr_hyb - fsr_full:+7.3f} nm", flush=True)
            rows.append(dict(res=res, route=name, shift_nm=shift,
                             frac_fsr=frac, fsr_full=fsr_full, fsr_hyb=fsr_hyb,
                             n_aligned=len(fa)))

    if missing:
        print(f"\n[c4] {len(missing)} input(s) not yet present:", flush=True)
        for m in missing[:8]:
            print(f"[c4]   {m}", flush=True)
        print("[c4] (jobs 11655143 = full device, 11655349 = 2-D hybrid partner)",
              flush=True)
    if not rows:
        print("\n[c4] nothing to compare yet.", flush=True)
        return 2

    # --- the discriminator ------------------------------------------------
    print("\n[c4] ===================== VERDICT =====================", flush=True)
    for name in HYB:
        sub = sorted([r for r in rows if r["route"] == name], key=lambda r: r["res"])
        if len(sub) < 2:
            continue
        offs = [r["shift_nm"] for r in sub]
        print(f"[c4] {name}: offset vs full device by resolution: "
              + ", ".join(f"res {r['res']:g} {r['shift_nm']:+.3f} nm"
                          for r in sub), flush=True)
        shrinking = abs(offs[-1]) < 0.5 * abs(offs[0])
        if shrinking:
            print(f"[c4]   -> the offset SHRINKS with resolution. The hybrid and "
                  f"the direct solve are", flush=True)
            print(f"[c4]      converging to the same device; the residual is "
                  f"discretization, not bookkeeping.", flush=True)
        else:
            print(f"[c4]   -> the offset does NOT shrink with resolution "
                  f"({offs[0]:+.3f} -> {offs[-1]:+.3f} nm).", flush=True)
            print(f"[c4]      That is a RESOLUTION-INDEPENDENT disagreement "
                  f"between the decomposition and", flush=True)
            print(f"[c4]      the direct solve — i.e. the flat phase lives in "
                  f"the hybrid's reference-plane", flush=True)
            print(f"[c4]      arithmetic, not in meep's treatment of the arcs. "
                  f"That localizes the 0.976pi.", flush=True)

    mono = [r for r in rows if r["route"] == "monolithic"]
    casc = [r for r in rows if r["route"] == "cascaded"]
    if mono and casc:
        print("\n[c4] WHICH ROUTE IS RIGHT (the question neither route can answer "
              "about itself):", flush=True)
        for res in sorted({r["res"] for r in rows}):
            m = next((r for r in mono if r["res"] == res), None)
            c = next((r for r in casc if r["res"] == res), None)
            if m and c:
                winner = "monolithic" if abs(m["shift_nm"]) < abs(c["shift_nm"]) else "cascaded"
                print(f"[c4]   res {res:g}: monolithic {m['shift_nm']:+.3f} nm, "
                      f"cascaded {c['shift_nm']:+.3f} nm  -> closer: {winner}",
                      flush=True)

    np.savez(os.path.join(OUT, "check4_comparison.npz"),
             resolutions=np.asarray([r["res"] for r in rows]),
             routes=np.asarray([r["route"] for r in rows]),
             shift_nm=np.asarray([r["shift_nm"] for r in rows]),
             frac_fsr=np.asarray([r["frac_fsr"] for r in rows]),
             fsr_full=np.asarray([r["fsr_full"] for r in rows]),
             fsr_hyb=np.asarray([r["fsr_hyb"] for r in rows]))
    print(f"\n[c4] saved -> {os.path.join(OUT, 'check4_comparison.npz')}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
