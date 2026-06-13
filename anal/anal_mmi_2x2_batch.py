import argparse
import pickle
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

# ── styling config ────────────────────────────────────────────────────────────
wl_min = 1.26;  wl_max = 1.36      # x-axis wavelength limits (µm)
T_min  = 0.9;  T_max  = 1.0      # transmission y-axis limits
sr_min = 0.485; sr_max = 0.535   # splitting ratio y-axis limits

label_font     = "Arial"
label_fontsize = 14
tick_font      = "Arial"
tick_fontsize  = 12

color_pre_bend  = "#000000"
color_post_bend = "#000000"

fig_width      = 4.0
fig_height     = 3.0
legend_fontsize = 12
# ─────────────────────────────────────────────────────────────────────────────


def _label(key):
    port, in_m, out_m = key
    return f"Port {port} In {in_m} → Out {out_m}"


def _add_series(ax, top, bot, tag, linestyle, color, wavelengths, db=False):
    for key in sorted(top.keys()):
        total = np.array(top[key]) + np.array(bot[key])
        y = 10 * np.log10(total) if db else total
        ax.plot(wavelengths, y, linestyle=linestyle, color=color, label=f"{_label(key)} {tag}")


def _add_ratio(ax, top, bot, tag, linestyle, color, wavelengths):
    for key in sorted(top.keys()):
        t = np.array(bot[key])
        total = t + np.array(top[key])
        ax.plot(wavelengths, t / total, linestyle=linestyle, color=color, label=f"{_label(key)} {tag}")


def _style_ax(ax):
    ax.tick_params(direction="in", labelsize=tick_fontsize)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontname(tick_font)
    ax.set_xlabel(ax.get_xlabel(), fontname=label_font, fontsize=label_fontsize)
    ax.set_ylabel(ax.get_ylabel(), fontname=label_font, fontsize=label_fontsize)
    ax.grid(True)
    leg = ax.get_legend()
    if leg is not None:
        plt.setp(leg.get_texts(), fontsize=legend_fontsize)


def process(sim_dir: Path) -> None:
    with open(sim_dir / "data" / "sim_results.pkl", "rb") as f:
        data = pickle.load(f)

    wavelengths = np.array(data["wavelengths"])

    if "powers_taper_top" in data:
        taper_top = data["powers_taper_top"]
        taper_bot = data["powers_taper_bot"]
        bend_top  = data.get("powers_bend_top", {})
        bend_bot  = data.get("powers_bend_bot", {})
    else:
        taper_top = data["powers_top_by_mode"]
        taper_bot = data["powers_bot_by_mode"]
        bend_top  = {}
        bend_bot  = {}

    has_bend = bool(bend_top)

    for db in [False, True]:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        _add_series(ax, taper_top, taper_bot, "(pre-bend)", "-", color_pre_bend, wavelengths, db=db)
        if has_bend:
            _add_series(ax, bend_top, bend_bot, "(post-bend)", "--", color_post_bend, wavelengths, db=db)
        ax.set_xlabel("Wavelength (µm)")
        ax.set_ylabel("Total Transmission (dB)" if db else "Total Transmission")
        ax.set_xlim(wl_min, wl_max)
        if db:
            ax.set_ylim(10 * np.log10(T_min), 10 * np.log10(T_max))
        else:
            ax.set_ylim(T_min, T_max)
        handles = [Line2D([0], [0], color=color_pre_bend,  ls="-",  label="Pre-bend"),
                   Line2D([0], [0], color=color_post_bend, ls="--", label="Post-bend")]
        ax.legend(handles=handles if has_bend else handles[:1])
        _style_ax(ax)
        fig.tight_layout()
        fname = "total_transmission_dB.png" if db else "total_transmission.png"
        fig.savefig(sim_dir / fname, dpi=300)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    _add_ratio(ax, taper_top, taper_bot, "(pre-bend)", "-", color_pre_bend, wavelengths)
    if has_bend:
        _add_ratio(ax, bend_top, bend_bot, "(post-bend)", "--", color_post_bend, wavelengths)
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Splitting Ratio (Top / Total)")
    ax.set_xlim(wl_min, wl_max)
    ax.set_ylim(sr_min, sr_max)
    handles = [Line2D([0], [0], color=color_pre_bend,  ls="-",  label="Pre-bend"),
               Line2D([0], [0], color=color_post_bend, ls="--", label="Post-bend")]
    ax.legend(handles=handles if has_bend else handles[:1])
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(sim_dir / "splitting_ratio.png", dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Batch-generate MMI 2x2 plots for every simulation folder in a directory."
    )
    parser.add_argument(
        "parent_dir",
        type=Path,
        help="Parent directory containing simulation subfolders (each must have data/sim_results.pkl).",
    )
    args = parser.parse_args()

    parent = args.parent_dir
    if not parent.is_dir():
        parser.error(f"Not a directory: {parent}")

    subdirs = sorted(d for d in parent.iterdir() if d.is_dir())
    if not subdirs:
        print(f"No subdirectories found in {parent}")
        return

    for subdir in subdirs:
        pkl = subdir / "data" / "sim_results.pkl"
        if not pkl.exists():
            print(f"Skipping {subdir.name}: no data/sim_results.pkl")
            continue
        print(f"Processing {subdir.name} ...")
        process(subdir)
        print(f"  → saved plots to {subdir}")


if __name__ == "__main__":
    main()
