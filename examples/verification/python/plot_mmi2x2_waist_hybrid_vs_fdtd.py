"""Transmission and splitting-ratio spectra: waisted 2x2 hybrid vs full FDTD.

Both legs are already measured at the SAME geometry (the 30.50 um design point,
Ww=2.24 Lw=6 Lt=6 L=21.5) and the same band:

  full-device FDTD : job 12797278, 51 wavelengths -> smatrix_waisted_*.pkl
  hybrid cascade   : job 12892459, 11 wavelengths (one shard per wavelength),
                     calibrate_body_neff=False (REFUSED for a waisted body)

The hybrid's 11 wavelengths are an exact subset of the FDTD's 51, so the FDTD
is drawn as a continuous reference curve and the hybrid as markers on top --
no interpolation anywhere in the comparison.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import SMatrix

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                    "example_output"))
# Geometry-parameterised: two design points now have hybrid legs (the 30.50 um
# incumbent and the 31.00 um confirmed point), and they must not be confused.
_WW = os.environ.get("MMI2X2_PLOT_WW", "2.24")
_LW = os.environ.get("MMI2X2_PLOT_LW", "6")
_LT = os.environ.get("MMI2X2_PLOT_LT", "6")
_L  = os.environ.get("MMI2X2_PLOT_L",  "21.5")
_NM = os.environ.get("MMI2X2_PLOT_NM", "10")
_GEO = f"Ww{_WW}_Lw{_LW}_Lt{_LT}_L{_L}"
HYB = os.path.join(_OUT, "mmi_2x2_tapered",
                   f"hybrid_vs_fdtd_{_GEO}_res40_ct8_nm{_NM}.json")
FDTD = os.path.join(_OUT, "mmi_2x2_tapered", f"smatrix_waisted_{_GEO}.pkl")
outdir = os.path.join(_OUT, "mmi_2x2_tapered")

d = json.load(open(HYB))
wl_h = np.array(d["wavelengths"])
Tb_h = np.array(d["hybrid"]["T_bar"])
Tc_h = np.array(d["hybrid"]["T_cross"])

sm = SMatrix.load(FDTD)
wl_f = np.linspace(1.530, 1.565, 51)
cb = sm.coeff(port="out_0", in_mode=1, out_mode=1)
cc = sm.coeff(port="out_1", in_mode=1, out_mode=1)
Tb_f, Tc_f = np.abs(cb) ** 2, np.abs(cc) ** 2

print(f"hybrid {len(wl_h)} wl, fdtd {len(wl_f)} wl")

# ---- 1. transmission spectra ---------------------------------------------
fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                       gridspec_kw={"height_ratios": [2, 1]})
ax[0].plot(wl_f, Tb_f, "-", color="C0", lw=1.8, label="bar — full FDTD")
ax[0].plot(wl_f, Tc_f, "-", color="C1", lw=1.8, label="cross — full FDTD")
ax[0].plot(wl_f, Tb_f + Tc_f, "-", color="C2", lw=1.4, label="total — full FDTD")
ax[0].plot(wl_h, Tb_h, "o", color="C0", ms=6, mfc="white", mew=1.6,
           label="bar — hybrid")
ax[0].plot(wl_h, Tc_h, "s", color="C1", ms=6, mfc="white", mew=1.6,
           label="cross — hybrid")
ax[0].plot(wl_h, Tb_h + Tc_h, "^", color="C2", ms=6, mfc="white", mew=1.6,
           label="total — hybrid")
ax[0].axhline(0.5, color="0.6", ls=":", lw=1)
ax[0].set_ylabel("transmission")
ax[0].set_title(f"Waisted 2×2 MMI — Ww{_WW} Lw{_LW} Lt{_LT} L{_L} "
                f"({float(_L)+9.0:.2f} µm total)\n"
                "hybrid FDTD/EME cascade vs full-device FDTD  —  "
                "no body-β calibration (refused for a waisted body)",
                fontsize=10)
ax[0].legend(fontsize=8, ncol=2, loc="center left")
ax[0].grid(alpha=0.3)

# residuals at the shared wavelengths only
idx = [int(np.argmin(np.abs(wl_f - w))) for w in wl_h]
ax[1].plot(wl_h, Tb_h - Tb_f[idx], "o-", color="C0", ms=5, label="bar")
ax[1].plot(wl_h, Tc_h - Tc_f[idx], "s-", color="C1", ms=5, label="cross")
ax[1].axhline(0, color="0.4", lw=1)
ax[1].set_xlabel("wavelength (µm)")
ax[1].set_ylabel("hybrid − FDTD")
ax[1].legend(fontsize=8)
ax[1].grid(alpha=0.3)
plt.tight_layout()
p1 = os.path.join(outdir, f"spectra_transmission_hybrid_vs_fdtd_{_GEO}.png")
plt.savefig(p1, dpi=160)
plt.close()
print(f"written -> {p1}")

# ---- 2. splitting ratio ---------------------------------------------------
sr_f = Tb_f / (Tb_f + Tc_f)
sr_h = Tb_h / (Tb_h + Tc_h)
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                       gridspec_kw={"height_ratios": [2, 1]})
ax[0].plot(wl_f, 100 * sr_f, "-", color="C3", lw=1.8, label="full FDTD")
ax[0].plot(wl_h, 100 * sr_h, "o", color="C3", ms=6, mfc="white", mew=1.6,
           label="hybrid")
ax[0].axhline(50, color="0.6", ls=":", lw=1, label="ideal 50/50")
ax[0].set_ylabel("bar fraction  T_bar/(T_bar+T_cross)  [%]")
ax[0].set_title(f"Splitting ratio — waisted 2×2 Lw{_LW} Lt{_LT} L{_L} "
                f"({float(_L)+9.0:.2f} µm), hybrid vs full-device FDTD",
                fontsize=10)
ax[0].legend(fontsize=8)
ax[0].grid(alpha=0.3)
ax[1].plot(wl_h, 100 * (sr_h - sr_f[idx]), "o-", color="C3", ms=5)
ax[1].axhline(0, color="0.4", lw=1)
ax[1].set_xlabel("wavelength (µm)")
ax[1].set_ylabel("Δ split [%]")
ax[1].grid(alpha=0.3)
plt.tight_layout()
p2 = os.path.join(outdir, f"spectra_splitting_ratio_hybrid_vs_fdtd_{_GEO}.png")
plt.savefig(p2, dpi=160)
plt.close()
print(f"written -> {p2}")

print(f"\nsplit at band centre: FDTD {100*sr_f[25]:.2f}%  "
      f"hybrid {100*sr_h[5]:.2f}%")
print(f"split range over band: FDTD {100*sr_f.min():.2f}-{100*sr_f.max():.2f}%  "
      f"hybrid {100*sr_h.min():.2f}-{100*sr_h.max():.2f}%")
