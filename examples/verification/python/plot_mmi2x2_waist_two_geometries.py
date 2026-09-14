"""Hybrid vs full-device FDTD on BOTH waisted 2x2 design points, side by side.

Why a two-geometry figure rather than two separate ones. Every hybrid-vs-FDTD
conclusion in this study came from one body (Lw6 Lt6 L21.5) until 2026-09-13.
Putting the confirmed 31.00 um point beside it shows the thing a single-device
plot cannot: the CROSS deficit reproduces on both, while the BAR error is
geometry-dependent and reverses sign. That distinction is the whole reason the
second geometry was worth running.

Reads only merged comparison JSONs and the stored full-device SMatrices -- no
solving.
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
                                    "example_output", "mmi_2x2_tapered"))

CASES = [
    ("incumbent  Lw6 Lt6 L21.5", "Ww2.24_Lw6_Lt6_L21.5", 30.50),
    ("CONFIRMED  Lw6 Lt5 L22.0", "Ww2.24_Lw6_Lt5_L22", 31.00),
]

fig, axes = plt.subplots(3, 2, figsize=(13, 11), sharex=True,
                         gridspec_kw={"height_ratios": [2, 1, 1.3]})

for col, (label, geo, total) in enumerate(CASES):
    d = json.load(open(os.path.join(
        _OUT, f"hybrid_vs_fdtd_{geo}_res40_ct8_nm10.json")))
    wl_h = np.array(d["wavelengths"])
    Tb_h = np.array(d["hybrid"]["T_bar"])
    Tc_h = np.array(d["hybrid"]["T_cross"])

    sm = SMatrix.load(os.path.join(_OUT, f"smatrix_waisted_{geo}.pkl"))
    wl_f = np.linspace(1.530, 1.565, 51)
    Tb_f = np.abs(sm.coeff(port="out_0", in_mode=1, out_mode=1)) ** 2
    Tc_f = np.abs(sm.coeff(port="out_1", in_mode=1, out_mode=1)) ** 2
    idx = [int(np.argmin(np.abs(wl_f - w))) for w in wl_h]

    # --- transmission
    a = axes[0, col]
    a.plot(wl_f, Tb_f, "-", color="C0", lw=1.8, label="bar — full FDTD")
    a.plot(wl_f, Tc_f, "-", color="C1", lw=1.8, label="cross — full FDTD")
    a.plot(wl_f, Tb_f + Tc_f, "-", color="C2", lw=1.4, label="total — full FDTD")
    a.plot(wl_h, Tb_h, "o", color="C0", ms=6, mfc="white", mew=1.6, label="bar — hybrid")
    a.plot(wl_h, Tc_h, "s", color="C1", ms=6, mfc="white", mew=1.6, label="cross — hybrid")
    a.plot(wl_h, Tb_h + Tc_h, "^", color="C2", ms=6, mfc="white", mew=1.6,
           label="total — hybrid")
    a.axhline(0.5, color="0.6", ls=":", lw=1)
    a.set_ylim(0.43, 1.0)
    a.set_title(f"{label}   ({total:.2f} µm)", fontsize=11)
    if col == 0:
        a.set_ylabel("transmission")
        a.legend(fontsize=7.5, ncol=2, loc="center left")
    a.grid(alpha=0.3)

    # --- residuals, SHARED y-scale so the two columns are comparable by eye
    a = axes[1, col]
    a.plot(wl_h, Tb_h - Tb_f[idx], "o-", color="C0", ms=5, label="bar")
    a.plot(wl_h, Tc_h - Tc_f[idx], "s-", color="C1", ms=5, label="cross")
    a.axhline(0, color="0.4", lw=1)
    a.set_ylim(-0.025, 0.025)
    if col == 0:
        a.set_ylabel("hybrid − FDTD")
        a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # --- splitting ratio
    sr_f = 100 * Tb_f / (Tb_f + Tc_f)
    sr_h = 100 * Tb_h / (Tb_h + Tc_h)
    a = axes[2, col]
    a.plot(wl_f, sr_f, "-", color="C3", lw=1.8, label="full FDTD")
    a.plot(wl_h, sr_h, "o", color="C3", ms=6, mfc="white", mew=1.6, label="hybrid")
    a.axhline(50, color="0.6", ls=":", lw=1, label="ideal 50/50")
    a.set_ylim(49.4, 52.8)
    a.set_xlabel("wavelength (µm)")
    if col == 0:
        a.set_ylabel("bar fraction  [%]")
        a.legend(fontsize=8)
    a.grid(alpha=0.3)
    mb = 100 * np.mean(np.abs(Tb_h - Tb_f[idx])) / np.mean(Tb_f)
    mc = 100 * np.mean(np.abs(Tc_h - Tc_f[idx])) / np.mean(Tc_f)
    a.text(0.03, 0.06, f"mean |ΔT|:  bar {mb:.2f}%   cross {mc:.2f}%",
           transform=a.transAxes, fontsize=9,
           bbox=dict(fc="white", ec="0.7", alpha=0.9))

fig.suptitle(
    "Waisted 2×2 MMI — hybrid FDTD/EME cascade vs full-device FDTD, two geometries\n"
    "identical config (N_MODES=10, C_TRANS=8, res 40, 11 wl); neither body can be "
    "β-calibrated.  The CROSS deficit reproduces on both; the BAR error does not.",
    fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.955])
out = os.path.join(_OUT, "spectra_hybrid_vs_fdtd_two_geometries.png")
fig.savefig(out, dpi=160)
print(f"written -> {out}")
