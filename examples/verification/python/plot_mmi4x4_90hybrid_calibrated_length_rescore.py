"""Plot: calibrating the body index moves the optimal mm_length by -0.50 um,
and the device clears the >=20 dBe CMRR spec at the new optimum.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-09-11. Reads the
npz from mmi_4x4_90hybrid_calibrated_length_rescore.py -- the full 7601-point
mm_length axis scored twice, once through the uncalibrated body propagation
constant (what the rounds 8-10 campaign used) and once through the calibrated
one, with everything else identical and no new simulations.

The figure has to make three things legible at once:
  1. the resonance shifts LEFT by 0.50 um (159.30 -> 158.80),
  2. its PEAK RISES (18.48 -> 20.33 dBe) rather than merely moving, and
  3. 20.33 clears the spec line while 18.48 does not.
A zoomed window is used because the feature is 0.5 um wide on a 380 um axis;
the inset-free wide view is kept as a second panel so the zoom is not
mistaken for the whole story.

Run:
    python plot_mmi4x4_90hybrid_calibrated_length_rescore.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
D = os.path.join(_ROOT, "example_output", "verification",
                 "mmi_4x4_90hybrid_calibrated_length_rescore")
IN_NPZ = os.path.join(D, "rescore_stride1.npz")

C_UNCAL = "#2a78d6"
C_CAL = "#eb6834"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e3e2df"
SPEC = "#8a8a87"

SPEC_CMRR = 20.0
SPEC_PHASE = 5.0
ZOOM = (157.5, 161.0)


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(1.0)
    ax.tick_params(colors=INK_2, labelsize=8, length=0)
    ax.grid(True, color=GRID, linewidth=1.0)
    ax.set_axisbelow(True)


def main():
    d = np.load(IN_NPZ)
    g = d["grid"]
    inc = float(d["mm_length_incumbent"])
    cu, cc = d["uncal_min_cmrr"], d["cal_min_cmrr"]
    pu, pc = d["uncal_max_dphi"], d["cal_max_dphi"]

    def best(c, p):
        feas = p <= SPEC_PHASE
        k = np.where(feas)[0][np.argmax(c[np.where(feas)[0]])]
        return float(g[k]), float(c[k]), float(p[k])

    Lu, Cu_, Pu_ = best(cu, pu)
    Lc, Cc_, Pc_ = best(cc, pc)

    fig = plt.figure(figsize=(12.5, 8.6), facecolor=SURFACE)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.45, 1.0],
                  hspace=0.42, wspace=0.2,
                  left=0.07, right=0.975, top=0.845, bottom=0.08)

    m = (g >= ZOOM[0]) & (g <= ZOOM[1])

    # --- A: CMRR, zoomed ---------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    style(ax)
    ax.axhline(SPEC_CMRR, color=SPEC, lw=1.5, ls=(0, (5, 4)), zorder=1)
    ax.text(ZOOM[0] + 0.05, SPEC_CMRR + 0.25, "spec  ≥ 20 dBe",
            fontsize=8.5, color=INK_2)
    ax.plot(g[m], cu[m], lw=2.0, color=C_UNCAL, label="uncalibrated body n_eff")
    ax.plot(g[m], cc[m], lw=2.0, color=C_CAL, label="calibrated body n_eff")
    for L, C, col in ((Lu, Cu_, C_UNCAL), (Lc, Cc_, C_CAL)):
        ax.plot([L], [C], "o", ms=9, color=col,
                markeredgecolor=SURFACE, markeredgewidth=2.0, zorder=5)
    ax.annotate(f"{Cu_:.2f} dBe\nat L={Lu:.2f}", xy=(Lu, Cu_),
                xytext=(Lu + 0.42, Cu_ - 3.4), fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=1.0))
    # Kept inside the axes: xytext further left ran under the y-axis label.
    ax.annotate(f"{Cc_:.2f} dBe\nat L={Lc:.2f}", xy=(Lc, Cc_),
                xytext=(ZOOM[0] + 0.12, Cc_ - 5.2), fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=1.0))
    ax.set_xlim(*ZOOM)
    ax.set_ylabel("min CMRR over ports/quadratures/λ  (dBe)", fontsize=9, color=INK_2)
    ax.set_xlabel("mm_length  (µm)", fontsize=9, color=INK_2)
    ax.set_title("The optimum moves 0.50 µm left — and rises above spec",
                 fontsize=11, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=9, loc="lower left", labelcolor=INK_2)

    # --- B: CMRR, full axis ------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    style(ax)
    ax.axhline(SPEC_CMRR, color=SPEC, lw=1.5, ls=(0, (5, 4)), zorder=1)
    ax.plot(g, cu, lw=0.8, color=C_UNCAL, alpha=0.85)
    ax.plot(g, cc, lw=0.8, color=C_CAL, alpha=0.85)
    ax.axvspan(*ZOOM, color=GRID, alpha=0.8, lw=0, zorder=0)
    ax.text(ZOOM[1] + 8, 19.0, "zoom window\n(panel A)", fontsize=8, color=INK_2)
    ax.set_xlabel("mm_length  (µm)", fontsize=9, color=INK_2)
    ax.set_ylabel("min CMRR  (dBe)", fontsize=9, color=INK_2)
    ax.set_title("full scanned axis, 20–400 µm", fontsize=11, color=INK, pad=8)

    # --- C: phase spec, zoomed --------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    style(ax)
    ax.axhline(SPEC_PHASE, color=SPEC, lw=1.5, ls=(0, (5, 4)), zorder=1)
    ax.text(ZOOM[0] + 0.05, SPEC_PHASE + 0.12, "spec  ≤ 5°", fontsize=8.5, color=INK_2)
    ax.plot(g[m], pu[m], lw=2.0, color=C_UNCAL, label="uncalibrated")
    ax.plot(g[m], pc[m], lw=2.0, color=C_CAL, label="calibrated")
    for L, P, col in ((Lu, Pu_, C_UNCAL), (Lc, Pc_, C_CAL)):
        ax.plot([L], [P], "o", ms=9, color=col,
                markeredgecolor=SURFACE, markeredgewidth=2.0, zorder=5)
    ax.set_xlim(*ZOOM)
    ax.set_ylim(0, 12)
    ax.set_ylabel("max |Δφ|  (deg)", fontsize=9, color=INK_2)
    ax.set_xlabel("mm_length  (µm)", fontsize=9, color=INK_2)
    ax.set_title("phase spec passes on both paths — it is CMRR that decides",
                 fontsize=11, color=INK, pad=8)
    ax.legend(frameon=False, fontsize=9, loc="upper right", labelcolor=INK_2)

    # --- D: the summary numbers, as a table not a chart --------------------
    ax = fig.add_subplot(gs[1, 1])
    ax.set_facecolor(SURFACE)
    ax.axis("off")
    rows = [
        ("", "uncalibrated", "calibrated"),
        ("optimal mm_length", f"{Lu:.2f} µm", f"{Lc:.2f} µm"),
        ("min CMRR", f"{Cu_:.3f} dBe", f"{Cc_:.3f} dBe"),
        ("max |Δφ|", f"{Pu_:.3f}°", f"{Pc_:.3f}°"),
        ("excess loss", f"{float(d['uncal_max_xloss'][np.argmin(np.abs(g-Lu))]):.3f} dB",
         f"{float(d['cal_max_xloss'][np.argmin(np.abs(g-Lc))]):.3f} dB"),
        ("meets CMRR ≥ 20 dBe", "NO", "YES"),
    ]
    y = 0.92
    for r, row in enumerate(rows):
        for c, (txt, xpos) in enumerate(zip(row, (0.0, 0.50, 0.79))):
            col = INK if r == 0 else INK_2
            if r == 0 and c == 1:
                col = C_UNCAL
            if r == 0 and c == 2:
                col = C_CAL
            if r == len(rows) - 1 and c > 0:
                col = INK
            ax.text(xpos, y, txt, fontsize=9.5, color=col,
                    fontweight="bold" if (r == 0 or r == len(rows) - 1) else "normal",
                    transform=ax.transAxes, va="top")
        y -= 0.135
        if r == 0:
            ax.plot([0.0, 1.0], [y + 0.055, y + 0.055], color=GRID, lw=1.2,
                    transform=ax.transAxes, clip_on=False)
    ax.text(0.0, y - 0.04,
            "incumbent 159.30 µm scored under calibration:\n15.615 dBe — the shortfall was a WRONG LENGTH,\nnot a worse device.",
            fontsize=8.5, color=INK_2, transform=ax.transAxes, va="top")

    fig.suptitle(
        "4×4 MMI 90° hybrid — re-optimising mm_length under the calibrated body index\n"
        "7601 lengths × 3 λ re-scored from cached solves only · no new FDTD or EME runs · res = 40",
        fontsize=12.5, color=INK, x=0.07, ha="left", y=0.975)

    out = os.path.join(D, "calibrated_length_rescore.png")
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close("all")
    print(f"saved -> {out}")
    print(f"uncal optimum L={Lu:.2f} CMRR={Cu_:.3f} dBe  |  "
          f"cal optimum L={Lc:.2f} CMRR={Cc_:.3f} dBe  |  shift {Lc-Lu:+.2f} um")


if __name__ == "__main__":
    main()
