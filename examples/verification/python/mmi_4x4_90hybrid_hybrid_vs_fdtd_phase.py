"""4x4 90 hybrid -- hybrid vs full-device FDTD, per-arm amplitude and phase,
at every body length where both engines have an 11-wavelength result.

WHAT THIS IS FOR. Engine agreement, not device performance. At these short
bodies the device is nowhere near its 4-way imaging point, so CMRR and
quadrature are meaningless here; what matters is whether the FDTD/EME cascade
reproduces what a full-device FDTD computes for the same geometry.

TWO PHASE READINGS, and the second is the one to quote.

  ABSOLUTE  d_j = arg(hybrid_j) - arg(fdtd_j) per arm. This carries a COMMON
            MODE shared by all four arms -- an overall reference-plane /
            propagation-phase offset between the two engines. It is real, but
            it is one number repeated four times, and it says nothing about
            whether the arms are correctly related to each other.

  RELATIVE  d_j minus the power-weighted circular mean of d over the arms.
            A 90-degree hybrid's output is a statement about the phases of the
            arms RELATIVE to one another, so this is the quantity the device
            is actually specified on. Removing the common mode previously cut
            a res-60 score by 90x -- see the hybrid-phase-metric note; scoring
            the absolute number instead mostly scores the reference plane.

WEIGHTING. Both readings are aggregated as sum_j w_j |d_j| with w_j the arm's
FDTD power share. An unweighted mean is actively misleading here: arms sitting
3-4 decades down carry the argument of numerical noise, not a phase, and they
dominate an unweighted average. The per-arm table is printed unweighted so the
weighting can never hide a disagreement -- read both.

MIRROR. The hybrid length sweep drives paper port 1; these FDTD references
drive port 4. Port 4's arm pattern is port 1's mirror (out5<->out8,
out6<->out7). This script sidesteps the issue entirely by reading the sweep's
saved `coeff[input_idx=3]`, i.e. the hybrid's OWN port-4 excitation -- it never
mirrors anything by hand. If `coeff` is missing, re-run the hybrid example.

Run:  python mmi_4x4_90hybrid_hybrid_vs_fdtd_phase.py
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
OUT = os.path.join(_ROOT, "example_output", "mmi_4x4_90hybrid")
SWEEP = os.path.join(OUT, "mmi_4x4_90hybrid_hybrid", "hybrid_length_sweep.npz")
PORT4_IN = 3                      # paper port 4 -> input_idx 3
ARMS = ["out5", "out6", "out7", "out8"]
POWER_FLOOR = 0.05                # "carries real light" bar, as in the sweep


def wrap(a):
    return ((np.asarray(a) + 180.0) % 360.0) - 180.0


def weighted_circular_mean(d_deg, w):
    """Power-weighted mean of angles -- circular, not arithmetic, so it stays
    correct when the arms straddle the +-180 branch cut."""
    z = np.sum(w * np.exp(1j * np.deg2rad(d_deg)))
    return float(np.rad2deg(np.angle(z)))


if not os.path.exists(SWEEP):
    sys.exit(f"missing {SWEEP}")
z = np.load(SWEEP, allow_pickle=True)
if "coeff" not in z.files:
    sys.exit("hybrid_length_sweep.npz has no 'coeff' -- re-run "
             "mmi_4x4_90hybrid_hybrid.py to save complex per-arm coefficients")
L = z["lengths"]
wl_h = z["wavelengths"]
coeff = z["coeff"]                # [in, out, length, wavelength] complex

rows = []
print(f"hybrid sweep: {len(L)} lengths {L[0]:.2f}-{L[-1]:.2f} um, "
      f"{len(wl_h)} wavelengths {wl_h[0]:.3f}-{wl_h[-1]:.3f} um\n")

for Lr in (5, 10, 15, 20, 25, 30, 35, 40):
    d = os.path.join(OUT, f"mmi_4x4_90hybrid_fdtd_L{Lr:.2f}_p4", "fdtd_results.npz")
    if not os.path.exists(d):
        print(f"L={Lr:<3} no FDTD reference yet"); continue
    f = np.load(d)
    wl_f = f["wavelengths"]
    if len(wl_f) != len(wl_h) or not np.allclose(wl_f, wl_h):
        print(f"L={Lr:<3} SKIP -- FDTD has {len(wl_f)} wavelengths, hybrid has "
              f"{len(wl_h)}; not comparable"); continue
    if Lr < L[0] - 1e-9:
        print(f"L={Lr:<3} SKIP -- below the sweep's {L[0]:.2f} um floor "
              f"(shared grid starts there); FDTD exists but hybrid does not")
        continue
    li = int(np.argmin(np.abs(L - Lr)))

    Ph = np.abs(coeff[PORT4_IN, :, li, :]) ** 2           # [4, nwl]
    phih = np.angle(coeff[PORT4_IN, :, li, :], deg=True)
    Pf = np.stack([f[f"T_port4_{a}"] for a in ARMS])
    phif = np.stack([f[f"phase_port4_{a}"] for a in ARMS])

    dabs = wrap(phih - phif)                               # [4, nwl]
    w = Pf / np.clip(Pf.sum(axis=0, keepdims=True), 1e-30, None)
    cm = np.array([weighted_circular_mean(dabs[:, k], w[:, k])
                   for k in range(dabs.shape[1])])
    drel = wrap(dabs - cm[None, :])

    wabs = float(np.mean(np.sum(w * np.abs(dabs), axis=0)))
    wrel = float(np.mean(np.sum(w * np.abs(drel), axis=0)))
    wrel_max = float(np.max(np.sum(w * np.abs(drel), axis=0)))
    nstrong = int((Pf.mean(axis=1) >= POWER_FLOOR).sum())

    print(f"=== L = {L[li]:.2f} um   ({nstrong}/4 arms >= {POWER_FLOOR:.0%} FDTD power) ===")
    print(f"{'arm':>5} {'|T|^2 hyb':>10} {'|T|^2 fdtd':>11} {'dP/P':>8} "
          f"{'|dphi| abs':>11} {'|dphi| rel':>11}   weight")
    for a in range(4):
        dp = (Ph[a].mean() - Pf[a].mean()) / max(Pf[a].mean(), 1e-30)
        print(f"{ARMS[a]:>5} {Ph[a].mean():10.5f} {Pf[a].mean():11.5f} "
              f"{dp:+7.1%} {np.mean(np.abs(dabs[a])):11.2f} "
              f"{np.mean(np.abs(drel[a])):11.2f}   {w[a].mean():.3f}"
              + ("" if Pf[a].mean() >= POWER_FLOOR else "   <- weak, phase is noise"))
    print(f"  common mode (band mean): {weighted_circular_mean(cm, np.ones_like(cm)/len(cm)):+.2f} deg")
    print(f"  POWER-WEIGHTED  absolute {wabs:6.2f} deg   RELATIVE {wrel:6.2f} deg "
          f"(worst wl {wrel_max:.2f})\n")
    rows.append((L[li], nstrong, wabs, wrel, wrel_max))

if rows:
    print("=" * 74)
    print(f"{'L (um)':>8} {'arms>=5%':>9} {'abs (deg)':>11} {'REL (deg)':>11} {'rel worst':>11}")
    for r in rows:
        print(f"{r[0]:8.2f} {r[1]:9d} {r[2]:11.2f} {r[3]:11.2f} {r[4]:11.2f}")
    print("\nRELATIVE is the column to quote -- see this module's docstring.")
