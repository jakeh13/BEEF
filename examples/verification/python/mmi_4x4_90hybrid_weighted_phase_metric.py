"""4x4 MMI 90-degree hybrid -- POWER-WEIGHTED hybrid-vs-FDTD phase agreement,
and what the resolution probes say about where the residual comes from.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-09-11.

WHY WEIGHTED. The unweighted mean |d phase| over the four output arms is a
misleading figure of merit for this device: at the reduced probe lengths two
of the four arms sit 3-4 decades below the dominant arm (|T|^2 ~ 5e-5 vs
~4.5e-1), and the "phase" of a field that small is not a physical quantity --
it is the argument of numerical noise. An unweighted mean lets those two
noise arms dominate the score (they carry the 33-60 deg residuals) and hides
the fact that every arm actually carrying light agrees to ~2 deg. Weighting
each arm by its share of the output power is the honest aggregate: it scores
the phase error the device's actual signal would experience.

    weighted_err = sum_j ( |T_j|^2 / sum_k |T_k|^2 ) * |d phase_j|

Reported both per-wavelength and pooled, for calibrated and uncalibrated, at
each probe length. Also expressed as a PERCENTAGE, with the denominator
stated explicitly -- "1.5%" of a phase is meaningless until you say of what.
Two readings are given because both are defensible and they differ by 4x:
  * of 360 deg (a full cycle -- the natural period of a phase), and
  * of 90 deg (this device's quadrature spacing -- the scale on which a
    hybrid's phase errors actually matter, and the scale the <= 5 deg spec
    is implicitly set against).

Run:
    python mmi_4x4_90hybrid_weighted_phase_metric.py
"""
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
IN_NPZ = os.path.join(_ROOT, "example_output", "verification",
                      "mmi_4x4_90hybrid_calib2_phase_check", "calib2_phase_check.npz")
PROBE_DIR = os.path.join(_ROOT, "example_output", "mmi_4x4_90hybrid", "calib_probe")

LENGTHS = [5.00, 10.00, 15.00]
ARMS = ["out5", "out6", "out7", "out8"]
WL = np.array([1.29, 1.31, 1.33])

# Full-device FDTD |T|^2 per arm, top port (paper port 4), per wavelength.
# Measured: jobs 13069782 (L=5), 13069783 (L=10), 13069784 (L=15).
FDTD_POWER = {
    5.00: dict(out5=[3.01573918e-05, 8.30517757e-05, 2.48077488e-05],
               out6=[5.02001100e-06, 7.10677569e-05, 1.84884508e-04],
               out7=[0.00526921, 0.00707336, 0.00513688],
               out8=[0.44720259, 0.46024323, 0.44437211]),
    10.00: dict(out5=[3.79802359e-06, 4.84054721e-05, 1.89428767e-04],
                out6=[0.00048889, 0.0008474, 0.00137776],
                out7=[0.05853651, 0.05838946, 0.06177298],
                out8=[0.55209692, 0.5560155, 0.55220692]),
    15.00: dict(out5=[6.39388625e-05, 3.37006280e-04, 6.07003288e-04],
                out6=[0.01311959, 0.01882186, 0.02060494],
                out7=[0.03559528, 0.04183967, 0.03510486],
                out8=[0.37710189, 0.38971857, 0.3202675]),
}


def weighted(L, key):
    """Per-wavelength power-weighted |d phase| (deg) at one length."""
    d = np.load(IN_NPZ)
    err = np.array([np.abs(d[f"L{L:.2f}_{a}_{key}"]) for a in ARMS])   # (arm, wl)
    pw = np.array([FDTD_POWER[L][a] for a in ARMS])                    # (arm, wl)
    w = pw / pw.sum(axis=0, keepdims=True)
    return (w * err).sum(axis=0), (w * err).sum() / len(WL), err, w


def main():
    print("=" * 78)
    print("POWER-WEIGHTED hybrid-vs-FDTD phase error, top port (paper port 4)")
    print("=" * 78)
    print(f"\n{'L (um)':>7} {'wl':>6} {'uncal (deg)':>12} {'CAL (deg)':>11}"
          f" {'cal % of 360':>13} {'cal % of 90':>12}")
    pooled = {}
    for L in LENGTHS:
        wu, pu, _, _ = weighted(L, "d_uncal")
        wc, pc, err_c, w = weighted(L, "d_cal")
        for i, wl in enumerate(WL):
            print(f"{L:7.0f} {wl:6.2f} {wu[i]:12.2f} {wc[i]:11.2f}"
                  f" {100*wc[i]/360:13.3f} {100*wc[i]/90:12.2f}")
        print(f"{L:7.0f} {'POOL':>6} {pu:12.2f} {pc:11.2f}"
              f" {100*pc/360:13.3f} {100*pc/90:12.2f}")
        pooled[L] = (pu, pc)
        # Which arm supplies the residual, after weighting.
        contrib = (w * err_c).mean(axis=1)
        share = contrib / contrib.sum()
        print("        residual budget (calibrated): "
              + "  ".join(f"{a}={100*s:4.1f}%" for a, s in zip(ARMS, share)))

    allc = np.mean([pooled[L][1] for L in LENGTHS])
    allu = np.mean([pooled[L][0] for L in LENGTHS])
    print(f"\nACROSS ALL THREE LENGTHS: uncal {allu:.2f} deg -> CAL {allc:.2f} deg")
    print(f"   calibrated = {100*allc/360:.3f} % of 360 deg  |  "
          f"{100*allc/90:.2f} % of 90 deg")
    print(f"   target 1.5% of 360 deg = {0.015*360:.2f} deg  -> "
          f"{'ALREADY MET' if allc < 0.015*360 else 'NOT met'}")
    print(f"   target 1.5% of  90 deg = {0.015*90:.2f} deg  -> "
          f"{'ALREADY MET' if allc < 0.015*90 else 'NOT met'}")

    # ---- resolution scaling of the calibration offset --------------------
    print("\n" + "=" * 78)
    print("RESOLUTION SCALING of the n_eff offset (mode 1)")
    print("=" * 78)
    rows = []
    for fn in sorted(os.listdir(PROBE_DIR)):
        if not fn.endswith("_mode01.npz"):
            continue
        d = np.load(os.path.join(PROBE_DIR, fn))
        if abs(float(d["wl_min"]) - 1.29) > 1e-9 or abs(float(d["wl_max"]) - 1.33) > 1e-9:
            continue          # band probes are a different axis; skip here
        if abs(float(d["bw_factor"]) - 2.0) > 1e-9:
            continue          # bandwidth probes likewise
        i = int(np.argmin(np.abs(d["wavelengths"] - 1.31)))
        rows.append((float(d["res"]), float(d["neff_modesolver"][i]),
                     float(d["neff_calibrated"][i]), float(d["offset"][i]),
                     float(d["resid"][i]), str(d["tag"])))
    rows.sort()
    print(f"\n{'res':>5} {'n_eff modesolver':>17} {'n_eff calibrated':>17}"
          f" {'offset':>9} {'fit resid':>10}")
    for r, ms, cal, off, res_, tag in rows:
        print(f"{r:5.0f} {ms:17.6f} {cal:17.6f} {off:+9.6f} {res_:10.2e}")
    if len(rows) >= 2:
        print("\npairwise scaling exponent p, assuming offset ~ res^-p:")
        for (r0, _, _, o0, _, _), (r1, _, _, o1, _, _) in zip(rows, rows[1:]):
            p = np.log(o0 / o1) / np.log(r1 / r0)
            print(f"   res {r0:.0f} -> {r1:.0f}:  offset {o0:+.6f} -> {o1:+.6f}"
                  f"   ratio {o0/o1:.2f}   p = {p:.2f}"
                  f"   (p=1 first-order, p=2 second-order/Yee)")
        print("\nBOTH engines move with resolution -- this is NOT a case of one")
        print("being right and the other wrong:")
        for r, ms, cal, off, _, _ in rows:
            print(f"   res {r:.0f}:  modesolver {ms:.6f}   FDTD {cal:.6f}"
                  f"   gap {off:+.6f}")


if __name__ == "__main__":
    main()
