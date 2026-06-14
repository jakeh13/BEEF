import pickle
import re
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

SWEEP_DIR = Path("/storage/scratch1/8/jhiesener3/shit/out/mmi_2x2_jmh_dilation_sweep_60")

# ── styling ───────────────────────────────────────────────────────────────────
label_font      = "Arial"
label_fontsize  = 14
tick_font       = "Arial"
tick_fontsize   = 12
legend_fontsize = 12
fig_width       = 4.0
fig_height      = 3.0
color_avg       = "#000000"
color_band      = "#888888"
# ─────────────────────────────────────────────────────────────────────────────


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


def _load_powers(pkl_path: Path):
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)
    wavelengths = np.array(data["wavelengths"])
    if "powers_taper_top" in data:
        top = data["powers_taper_top"]
        bot = data["powers_taper_bot"]
    else:
        top = data["powers_top_by_mode"]
        bot = data["powers_bot_by_mode"]
    return wavelengths, top, bot


def main():
    if not SWEEP_DIR.is_dir():
        print(f"Sweep directory not found: {SWEEP_DIR}")
        return

    # Discover dilation_Nnm subdirectories, sorted by dilation value
    pattern = re.compile(r"^dilation_(\d+)nm$")
    runs = []
    for d in SWEEP_DIR.iterdir():
        m = pattern.match(d.name)
        if m and d.is_dir():
            pkl = d / "data" / "sim_results.pkl"
            if pkl.exists():
                runs.append((int(m.group(1)), d, pkl))
            else:
                print(f"Skipping {d.name}: no data/sim_results.pkl")
    runs.sort(key=lambda x: x[0])

    if not runs:
        print(f"No completed dilation runs found in {SWEEP_DIR}")
        return

    dilations_nm = []
    trans_min, trans_max, trans_avg = [], [], []
    ratio_min, ratio_max, ratio_avg = [], [], []

    for dil_nm, subdir, pkl in runs:
        wavelengths, top, bot = _load_powers(pkl)

        # Accumulate total transmission and splitting ratio over all keys
        all_total = []
        all_ratio = []
        for key in top:
            t = np.array(top[key])
            b = np.array(bot[key])
            total = t + b
            all_total.append(total)
            all_ratio.append(b / total)

        total_arr = np.concatenate(all_total)  # all wavelengths × all keys
        ratio_arr = np.concatenate(all_ratio)

        dilations_nm.append(dil_nm)
        trans_min.append(np.min(total_arr))
        trans_max.append(np.max(total_arr))
        trans_avg.append(np.mean(total_arr))
        ratio_min.append(np.min(ratio_arr))
        ratio_max.append(np.max(ratio_arr))
        ratio_avg.append(np.mean(ratio_arr))
        print(f"  dilation={dil_nm} nm  T_avg={trans_avg[-1]:.4f}  SR_avg={ratio_avg[-1]:.4f}")

    x = np.array(dilations_nm)
    trans_min = np.array(trans_min)
    trans_max = np.array(trans_max)
    trans_avg = np.array(trans_avg)
    ratio_min = np.array(ratio_min)
    ratio_max = np.array(ratio_max)
    ratio_avg = np.array(ratio_avg)

    legend_handles = [
        Line2D([0], [0], color=color_avg, ls="-",  label="Average"),
        Line2D([0], [0], color=color_band, ls="--", label="Min / Max"),
    ]

    # ── Total transmission (linear) ───────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.fill_between(x, trans_min, trans_max, color=color_band, alpha=0.3)
    ax.plot(x, trans_min, "--", color=color_band)
    ax.plot(x, trans_max, "--", color=color_band)
    ax.plot(x, trans_avg, "-",  color=color_avg, label="Average")
    ax.set_xlabel("Dilation (nm)")
    ax.set_ylabel("Total Transmission")
    ax.legend(handles=legend_handles)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(SWEEP_DIR / "transmission_vs_dilation.png", dpi=300)
    plt.close(fig)

    # ── Total transmission (dB) ───────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.fill_between(x, 10*np.log10(trans_min), 10*np.log10(trans_max), color=color_band, alpha=0.3)
    ax.plot(x, 10*np.log10(trans_min), "--", color=color_band)
    ax.plot(x, 10*np.log10(trans_max), "--", color=color_band)
    ax.plot(x, 10*np.log10(trans_avg), "-",  color=color_avg)
    ax.set_xlabel("Dilation (nm)")
    ax.set_ylabel("Total Transmission (dB)")
    ax.legend(handles=legend_handles)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(SWEEP_DIR / "transmission_vs_dilation_dB.png", dpi=300)
    plt.close(fig)

    # ── Splitting ratio ───────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.fill_between(x, ratio_min, ratio_max, color=color_band, alpha=0.3)
    ax.plot(x, ratio_min, "--", color=color_band)
    ax.plot(x, ratio_max, "--", color=color_band)
    ax.plot(x, ratio_avg, "-",  color=color_avg)
    ax.set_xlabel("Dilation (nm)")
    ax.set_ylabel("Splitting Ratio (Bot / Total)")
    ax.legend(handles=legend_handles)
    _style_ax(ax)
    fig.tight_layout()
    fig.savefig(SWEEP_DIR / "splitting_ratio_vs_dilation.png", dpi=300)
    plt.close(fig)

    print(f"\nSaved plots to {SWEEP_DIR}")


if __name__ == "__main__":
    main()
