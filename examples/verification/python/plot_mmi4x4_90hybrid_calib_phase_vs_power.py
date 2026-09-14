"""Plot: n_eff calibration collapses the hybrid-vs-FDTD phase error, and the
residual lives only in the arms that carry no power.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-09-11. Reads the
phase comparison written by mmi_4x4_90hybrid_calib2_phase_check.py (job
13102478) and pairs it with the full-device FDTD arm powers measured by
jobs 13069782/783/784, so the two halves of the claim are on one figure:

  - calibration drops mean |d phase| on the dominant arm from ~175 deg to
    ~2 deg at L=15 (and comparably at 5/10), and
  - the arms where calibration does NOT visibly help (out5, out6 at the
    short lengths) are exactly the arms sitting at |T|^2 ~ 1e-8..1e-3,
    where the phase of a near-zero field is not a physical quantity.

Deliberately NOT a dual-axis chart: power is the scatter's x-axis in the
lower panel and a printed annotation in the upper facets, never a second
y-scale opposite phase error.

Run:
    python plot_mmi4x4_90hybrid_calib_phase_vs_power.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
IN_NPZ = os.path.join(_ROOT, "example_output", "verification",
                      "mmi_4x4_90hybrid_calib2_phase_check", "calib2_phase_check.npz")
OUT_DIR = os.path.join(_ROOT, "example_output", "verification",
                       "mmi_4x4_90hybrid_calib2_phase_check")
os.makedirs(OUT_DIR, exist_ok=True)

# Validated categorical slots 1-2 (see the dataviz reference palette); this
# pair passes all-pairs CVD (dE 24.7 protan) and normal-vision (33.6) on the
# light surface, so the two series are safe without relying on position alone.
C_UNCAL = "#2a78d6"     # slot 1, blue
C_CAL = "#eb6834"       # slot 2, orange
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e3e2df"

LENGTHS = [5.00, 10.00, 15.00]
ARMS = ["out5", "out6", "out7", "out8"]

# Full-device FDTD |T|^2 per arm, top port (paper port 4), at 1.29/1.31/1.33 um.
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

NOISE_FLOOR = 1e-3   # |T|^2 below which "phase" is a near-zero field's phase


def load_phase_errors():
    d = np.load(IN_NPZ)
    out = {}
    for L in LENGTHS:
        for a in ARMS:
            du = np.abs(d[f"L{L:.2f}_{a}_d_uncal"])
            dc = np.abs(d[f"L{L:.2f}_{a}_d_cal"])
            out[(L, a)] = (float(np.mean(du)), float(np.mean(dc)))
    return out


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(1.0)
    ax.tick_params(colors=INK_2, labelsize=8, length=0)
    ax.grid(True, axis="y", color=GRID, linewidth=1.0, linestyle="-")
    ax.set_axisbelow(True)


def main():
    err = load_phase_errors()

    fig = plt.figure(figsize=(12.5, 8.4), facecolor=SURFACE)
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.0, 1.15],
                  hspace=0.52, wspace=0.22,
                  left=0.07, right=0.975, top=0.855, bottom=0.085)

    # ---- Row 1: per-length grouped bars, phase error only -----------------
    ymax = max(max(v) for v in err.values()) * 1.18
    for col, L in enumerate(LENGTHS):
        ax = fig.add_subplot(gs[0, col])
        style(ax)
        x = np.arange(len(ARMS))
        w = 0.36
        gap = 0.02          # 2px-equivalent surface gap between adjacent bars
        u = [err[(L, a)][0] for a in ARMS]
        c = [err[(L, a)][1] for a in ARMS]
        ax.bar(x - w / 2 - gap / 2, u, w, color=C_UNCAL, label="uncalibrated")
        ax.bar(x + w / 2 + gap / 2, c, w, color=C_CAL, label="calibrated")

        # Selective direct labels: the calibrated value only, since "how small
        # did it get" is the number being read. Uncalibrated is read off the axis.
        for xi, cv in zip(x, c):
            ax.text(xi + w / 2 + gap / 2, cv + ymax * 0.025, f"{cv:.0f}°",
                    ha="center", va="bottom", fontsize=8, color=INK_2)

        ax.set_xticks(x)
        ax.set_xticklabels(
            [f"{a}\n|T|²={np.mean(FDTD_POWER[L][a]):.1e}" for a in ARMS],
            fontsize=8, color=INK_2)
        ax.set_ylim(0, ymax)
        ax.set_title(f"L = {L:.0f} µm", fontsize=11, color=INK, pad=8)
        if col == 0:
            ax.set_ylabel("mean |Δφ| vs full-device FDTD  (deg)",
                          fontsize=9, color=INK_2)
        else:
            ax.set_yticklabels([])
        if col == 0:
            ax.legend(frameon=False, fontsize=9, loc="upper left",
                      labelcolor=INK_2)

    # ---- Row 2: the relationship -- phase error vs arm power --------------
    ax = fig.add_subplot(gs[1, :])
    style(ax)
    ax.grid(True, axis="x", color=GRID, linewidth=1.0)

    pw, eu, ec = [], [], []
    for L in LENGTHS:
        for a in ARMS:
            pw.append(float(np.mean(FDTD_POWER[L][a])))
            eu.append(err[(L, a)][0])
            ec.append(err[(L, a)][1])
    pw, eu, ec = np.array(pw), np.array(eu), np.array(ec)

    ax.axvspan(1e-6, NOISE_FLOOR, color=GRID, alpha=0.45, lw=0, zorder=0)
    ax.text(3e-6, 172, "arms carrying no power\n(|T|² < 10⁻³) — phase of a\nnear-zero field",
            fontsize=8.5, color=INK_2, va="top", ha="left")

    ax.scatter(pw, eu, s=78, color=C_UNCAL, label="uncalibrated",
               zorder=3, edgecolors=SURFACE, linewidths=2.0)
    ax.scatter(pw, ec, s=78, color=C_CAL, label="calibrated",
               zorder=4, edgecolors=SURFACE, linewidths=2.0)
    for p, a_, b_ in zip(pw, eu, ec):
        ax.plot([p, p], [a_, b_], color=INK_2, lw=1.0, alpha=0.28, zorder=2)

    # Selective direct labels: the single most extreme improvement only.
    imax = int(np.argmax(eu - ec))
    ax.annotate(f"{eu[imax]:.0f}° → {ec[imax]:.0f}°",
                xy=(pw[imax], ec[imax]), xytext=(pw[imax] * 0.42, ec[imax] + 34),
                fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=1.0))

    ax.set_xscale("log")
    ax.set_xlim(1.5e-6, 1.2)
    ax.set_ylim(-6, 190)
    ax.set_xlabel("full-device FDTD power in that arm,  |T|²   (log scale)",
                  fontsize=9, color=INK_2)
    ax.set_ylabel("mean |Δφ| vs full-device FDTD  (deg)", fontsize=9, color=INK_2)
    # Upper-right overlaps the four high-error uncalibrated points; the
    # left of the noise band is the only genuinely empty region.
    ax.legend(frameon=False, fontsize=9, loc="center left", labelcolor=INK_2)
    ax.set_title("Every arm that carries real power agrees to ~2° once calibrated",
                 fontsize=11, color=INK, pad=8)

    fig.suptitle(
        "4×4 MMI 90° hybrid — n_eff calibration removes the hybrid-vs-FDTD phase drift\n"
        "top port (paper port 4) driven · 12 arm×length points · Δn_eff = 0.0458 "
        "(mode solver → meep time-stepper)",
        fontsize=12.5, color=INK, x=0.07, ha="left", y=0.975)

    out = os.path.join(OUT_DIR, "calib_phase_vs_power.png")
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close("all")
    print(f"saved -> {out}")

    # Table view -- identity never color-alone, and the numbers stay readable
    # without the figure.
    print(f"\n| L (µm) | arm | FDTD |T|² | mean |Δφ| uncal | mean |Δφ| cal |")
    print("|---|---|---|---|---|")
    for L in LENGTHS:
        for a in ARMS:
            u, c = err[(L, a)]
            print(f"| {L:.0f} | {a} | {np.mean(FDTD_POWER[L][a]):.2e} | {u:.1f}° | {c:.1f}° |")


if __name__ == "__main__":
    main()
