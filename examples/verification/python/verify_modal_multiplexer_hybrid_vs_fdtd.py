"""Modal multiplexer: hybrid cascade vs full-device FDTD vs the manuscript seed.

The comparison `modal_multiplexer_validation.md` has wanted since gap 3 closed
and that has only ever been done by hand, at 2 wavelengths, in a session
transcript (gap 7's 2026-08-09 entry). Both engines now produce 11-wavelength
S-matrices on geometry verified against the fabricated GDS, so this makes it a
real, repeatable check.

THREE COMPARISONS, deliberately separated because they answer different things:

  1. hybrid vs FDTD, per wavelength and channel — an ENGINE-agreement check,
     the modal-multiplexer analogue of the sibling MMIs' gap 4. Both sides are
     this repo's own code on identical geometry, so a disagreement is a bug or
     a basis-truncation effect, not physics.
  2. each engine vs the MANUSCRIPT's seed numbers — a PHYSICS check against an
     external reference.
  3. the size of (1) against the size of (2) — the discriminator for gap 9's
     open hypothesis that the manuscript's simulated numbers, being Lumerical
     EME, systematically over-predict transmission on the radiating channel.
     If the hybrid (EME bodies) sits between the FDTD and the paper, and does
     so most on TE01, that hypothesis survives. If hybrid and FDTD agree and
     both miss the paper, it does not.

⚠ THE MANUSCRIPT'S NUMBERS ARE TWO PASSES. The measured test structure is a
back-to-back pair of identical mirrored multiplexers (Fig. 3d), so S13 and S24
each traverse the device TWICE. A single device is half the dB. This is not a
footnote: comparing a single-pass result to the two-pass figure produced an
apparent 0.07 dB "match" during this device's own history that did not exist
(gap 9). The halving is applied here once, explicitly, so no caller repeats it.

⚠ PORT MAPPING, from Fig. 3(a): the three single-mode ports are labelled
NC / A / B top to bottom, so this repo's `out_0` = top = NC (drawn but
unused), `out_1` = mid = A, `out_2` = bot = B. TE00 is routed to A and TE01
to B; `top` carries no channel and its power is lost, not an output.

⚠ BAND. The manuscript measures 1285-1325 nm (line 166), NOT the full
1260-1360 O-band both engines sweep. Means are reported over both, but the
measurement band is the one comparable to the paper.

Run:
    python verification/python/verify_modal_multiplexer_hybrid_vs_fdtd.py
"""
import argparse
import os
import sys

import numpy as np


# --- manuscript seed, Hiesener et al., lines 268-272 -----------------------
# "compared with -2.79 (-3.30) dB, -5.99 (-5.51) dB, and -32.02 (-19.49) dB,
#  respectively, for the seed" -- measured (simulated), TWO passes.
PAPER_2PASS = {
    "TE00_measured": -2.79, "TE00_simulated": -3.30,
    "TE01_measured": -5.99, "TE01_simulated": -5.51,
    "XT_measured": -32.02,  "XT_simulated": -19.49,
}
PAPER_1PASS = {k: v / 2.0 for k, v in PAPER_2PASS.items()}

BAND_MEAS = (1.285, 1.325)      # manuscript measurement band (µm)
BAND_FULL = (1.26, 1.36)        # what both engines sweep

CHANNELS = [                     # (label, input mode index, output port index)
    ("TE00 -> A (mid)", 1, 1),
    ("TE01 -> B (bot)", 2, 2),
]
CROSSTALK = [
    ("XT TE01 -> A", 2, 1),
    ("XT TE00 -> B", 1, 2),
]
NC_PORT = 0


def db(x):
    return 10.0 * np.log10(np.maximum(np.asarray(x, dtype=float), 1e-12))


# THE TWO ENGINES NAME THEIR PORTS DIFFERENTLY, and the mapping is verified
# from source rather than assumed -- a silent swap here would invert the whole
# comparison, which is exactly the false alarm this device already produced
# once (gap 7, 2026-08-09, an apparent top/bot "transposition" that turned out
# to be the FDTD reference being wrong).
#
#   modal_multiplexer_fdtd.py   -> "out_0" / "out_1" / "out_2", and its own
#       section-5 comment says these are "ordered by output_monitor_centers"
#       and "relabelled top/mid/bot ... matching modal_multiplexer_hybrid.py's
#       PORTS convention".
#   modal_multiplexer_hybrid.py -> PORTS = [("top", 0), ("mid", 1), ("bot", 2)]
#
# Both derive from BRANCH_SPECS ordered top(+WO) / mid(0) / bot(-WO), so
# index k is the same physical branch in both. Ordered top -> mid -> bot.
PORT_ALIASES = (("out_0", "top"), ("out_1", "mid"), ("out_2", "bot"))


def load(path):
    from beef import SMatrix
    s = SMatrix.load(path)
    wl = np.asarray(s.wavelengths)
    keys = set(s.keys())

    def port_name(k):
        for cand in PORT_ALIASES[k]:
            if (cand, 1, 1) in keys:
                return cand
        raise KeyError(
            f"{path}: no recognised name for output branch {k}; "
            f"keys are {sorted(keys)}. Expected one of {PORT_ALIASES[k]}.")

    names = [port_name(k) for k in range(3)]
    P = {m: np.array([[s.transmission(names[k], m, 1)[j] for k in range(3)]
                      for j in range(len(wl))]) for m in (1, 2)}
    return wl, P, names


def band_mean(wl, vals, band):
    m = (wl >= band[0] - 1e-9) & (wl <= band[1] + 1e-9)
    if not m.any():
        return float("nan"), 0
    return float(db(vals[m]).mean()), int(m.sum())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    root = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                        "..", "..", "..", ".."))
    d = os.path.join(root, "example_output", "modal_multiplexer")
    ap.add_argument("--fdtd", default=os.path.join(
        d, "modal_multiplexer_fdtd", "smatrix.pkl"))
    ap.add_argument("--hybrid", default=os.path.join(
        d, "modal_multiplexer_hybrid", "smatrix_hybrid.pkl"))
    args = ap.parse_args()

    missing = [p for p in (args.fdtd, args.hybrid) if not os.path.exists(p)]
    if missing:
        print("missing S-matrix:\n" + "\n".join("  " + m for m in missing))
        return 1

    wl_f, P_f, names_f = load(args.fdtd)
    wl_h, P_h, names_h = load(args.hybrid)
    print(f"FDTD   : {len(wl_f)} wavelengths {wl_f.min():.3f}-{wl_f.max():.3f} µm"
          f"   ports {names_f} (top, mid, bot)")
    print(f"hybrid : {len(wl_h)} wavelengths {wl_h.min():.3f}-{wl_h.max():.3f} µm"
          f"   ports {names_h} (top, mid, bot)")

    shared = np.array(sorted(set(np.round(wl_f, 6)) & set(np.round(wl_h, 6))))
    print(f"shared : {len(shared)} wavelengths\n")
    if len(shared) == 0:
        print("no shared wavelengths — nothing to compare.")
        return 1

    fi = {round(w, 6): i for i, w in enumerate(wl_f)}
    hi = {round(w, 6): i for i, w in enumerate(wl_h)}

    # --- 1. per-wavelength engine agreement -------------------------------
    print("=" * 74)
    print("1. ENGINE AGREEMENT  (both this repo, identical geometry)")
    print("=" * 74)
    worst = (0.0, None)
    for label, mi, po in CHANNELS + CROSSTALK:
        print(f"\n  {label}")
        print("     wl      FDTD(dB)  hybrid(dB)   delta      FDTD(T)  hybrid(T)")
        for w in shared:
            a = P_f[mi][fi[round(w, 6)], po]
            b = P_h[mi][hi[round(w, 6)], po]
            dd = db(b) - db(a)
            print(f"   {w:.3f}   {db(a):+8.2f}  {db(b):+8.2f}  {dd:+8.2f}   "
                  f"{a:8.4f}  {b:8.4f}")
            if abs(dd) > abs(worst[0]):
                worst = (float(dd), f"{label} at {w:.3f} µm")
    print(f"\n  worst engine disagreement: {worst[0]:+.2f} dB ({worst[1]})")

    # --- 2. each engine vs the manuscript ---------------------------------
    print("\n" + "=" * 74)
    print("2. VS MANUSCRIPT SEED  (paper figures HALVED — they are two passes)")
    print("=" * 74)
    for band, bl in ((BAND_MEAS, "measurement band 1285-1325 nm"),
                     (BAND_FULL, "full O-band 1260-1360 nm")):
        print(f"\n  -- {bl} --")
        print(f"  {'channel':18s} {'FDTD':>8s} {'hybrid':>8s} "
              f"{'paper sim':>10s} {'paper meas':>11s}   {'F-gap':>7s} {'H-gap':>7s}")
        for label, mi, po in CHANNELS:
            key = "TE00" if mi == 1 else "TE01"
            f_m, nf = band_mean(wl_f, P_f[mi][:, po], band)
            h_m, nh = band_mean(wl_h, P_h[mi][:, po], band)
            psim = PAPER_1PASS[f"{key}_simulated"]
            pmea = PAPER_1PASS[f"{key}_measured"]
            print(f"  {label:18s} {f_m:+8.2f} {h_m:+8.2f} {psim:+10.2f} "
                  f"{pmea:+11.2f}   {f_m-psim:+7.2f} {h_m-psim:+7.2f}")
        for label, mi, po in CROSSTALK:
            f_m, _ = band_mean(wl_f, P_f[mi][:, po], band)
            h_m, _ = band_mean(wl_h, P_h[mi][:, po], band)
            print(f"  {label:18s} {f_m:+8.2f} {h_m:+8.2f} "
                  f"{PAPER_1PASS['XT_simulated']:+10.2f} "
                  f"{PAPER_1PASS['XT_measured']:+11.2f}")
        print(f"  (F-gap / H-gap = engine minus paper SIMULATED, "
              f"positive = engine transmits MORE)")

    # --- 3. the discriminator ---------------------------------------------
    print("\n" + "=" * 74)
    print("3. GAP 9 DISCRIMINATOR  — is the residual a solver-class effect?")
    print("=" * 74)
    print("  The manuscript's simulated numbers come from Lumerical EME, which")
    print("  on a truncated basis cannot radiate freely and should therefore")
    print("  OVER-predict transmission -- most on the channel that radiates most.")
    print("  The hybrid propagates both MMI bodies in EME; the FDTD does not.\n")
    verdict = []
    for label, mi, po in CHANNELS:
        key = "TE00" if mi == 1 else "TE01"
        f_m, _ = band_mean(wl_f, P_f[mi][:, po], BAND_MEAS)
        h_m, _ = band_mean(wl_h, P_h[mi][:, po], BAND_MEAS)
        psim = PAPER_1PASS[f"{key}_simulated"]
        # A nan hybrid mean (no samples inside the band) must NOT read as a
        # failed test -- it is no test at all. Distinguish the two.
        if not np.isfinite(h_m) or not np.isfinite(f_m):
            verdict.append(None)
            print(f"  {label:18s} FDTD {f_m:+7.2f}  hybrid {h_m:+7.2f}  "
                  f"paper {psim:+7.2f}   hybrid between: NO DATA in band")
            continue
        between = (f_m < h_m <= psim + 1e-9) or (psim <= h_m < f_m)
        verdict.append(between)
        print(f"  {label:18s} FDTD {f_m:+7.2f}  hybrid {h_m:+7.2f}  "
              f"paper {psim:+7.2f}   hybrid between: {'YES' if between else 'no'}")
    # radiated power, the mechanism's own observable
    for eng, wl_, P_ in (("FDTD", wl_f, P_f), ("hybrid", wl_h, P_h)):
        m = (wl_ >= BAND_MEAS[0] - 1e-9) & (wl_ <= BAND_MEAS[1] + 1e-9)
        if m.any():
            print(f"  {eng:6s} total recovered power: "
                  f"TE0 {P_[1][m].sum(1).mean():.3f}  TE1 {P_[2][m].sum(1).mean():.3f}"
                  f"   (NC port takes TE0 {P_[1][m, NC_PORT].mean():.3f} / "
                  f"TE1 {P_[2][m, NC_PORT].mean():.3f})")
    print()
    tested = [v for v in verdict if v is not None]
    if not tested:
        print("  => VERDICT WITHHELD: the hybrid has no wavelength samples inside")
        print("     the manuscript's measurement band, so this test did not run.")
        return 0
    verdict = tested
    if all(verdict):
        print("  => hybrid sits between FDTD and the paper on BOTH channels:")
        print("     consistent with the EME-over-prediction hypothesis.")
    elif any(verdict):
        print("  => hybrid sits between on ONE channel only: partial support,")
        print("     look at which channel and whether it is the radiating one.")
    else:
        print("  => hybrid does NOT sit between on either channel: the")
        print("     solver-class hypothesis is NOT supported by this data.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
