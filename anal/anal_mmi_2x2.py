import numpy as np
import pickle
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

output_subfolder = "mmi_2x2_archana_bend"
# mmi_2x2_jmh_bend
# FDTD_GDS_BEND_res_40_gds_mmi2x2 FDTD_GDS_BEND_res_40_gds_splitter_2x2_f53

# mmi_2x2_archana_bend
# FDTD_GDS_BEND_res_40_gds_PSO_MMI FDTD_GDS_BEND_res_40_gds_Splitter_60
simulation_folder = "FDTD_GDS_BEND_res_40_gds_Splitter_60"
OUTPUT_FOLDER = f"/storage/scratch1/8/jhiesener3/shit/out/{output_subfolder}/{simulation_folder}"

wl_min = 1.5; wl_max = 1.6  # x-axis wavelength limits (µm)
T_min  = 0.9; T_max  = 1      # transmission y-axis limits
sr_min = 0.485; sr_max = 0.535  # splitting ratio y-axis limits

label_font     = "Arial"  # font for axis labels
label_fontsize = 14       # font size for axis labels
tick_font      = "Arial"  # font for tick numbers
tick_fontsize  = 12       # font size for tick numbers

color_pre_bend  = "#000000"  # line color for pre-bend series
color_post_bend = "#000000"  # line color for post-bend series

fig_width  = 4.0  # figure width (inches)
fig_height = 3.0  # figure height (inches)
legend_fontsize = 12  # font size for legend text

with open(f"{OUTPUT_FOLDER}/data/sim_results.pkl", "rb") as f:
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

keys = sorted(taper_top.keys())
has_bend = bool(bend_top)


def _label(key):
    port, in_m, out_m = key
    return f"Port {port} In {in_m} → Out {out_m}"


def _add_series(ax, top, bot, tag, linestyle, color, db=False):
    for key in sorted(top.keys()):
        total = np.array(top[key]) + np.array(bot[key])
        y = 10 * np.log10(total) if db else total
        ax.plot(wavelengths, y, linestyle=linestyle, color=color, label=f"{_label(key)} {tag}")


def _add_ratio(ax, top, bot, tag, linestyle, color):
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


for db in [False, True]:
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    _add_series(ax, taper_top, taper_bot, "(pre-bend)", "-", color_pre_bend, db=db)
    if has_bend:
        _add_series(ax, bend_top, bend_bot, "(post-bend)", "--", color_post_bend, db=db)
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Total Transmission (dB)" if db else "Total Transmission")
    ax.set_xlim(wl_min, wl_max)
    ax.set_ylim(10*np.log10(T_min), 10*np.log10(T_max)) if db else ax.set_ylim(T_min, T_max)
    handles = [Line2D([0], [0], color=color_pre_bend,  ls="-",  label="Pre-bend"),
               Line2D([0], [0], color=color_post_bend, ls="--", label="Post-bend")]
    ax.legend(handles=handles if has_bend else handles[:1])
    _style_ax(ax)
    fig.tight_layout()
    fname = "total_transmission_dB.png" if db else "total_transmission.png"
    fig.savefig(f"{OUTPUT_FOLDER}/{fname}", dpi=300)
    plt.close(fig)

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
_add_ratio(ax, taper_top, taper_bot, "(pre-bend)", "-", color_pre_bend)
if has_bend:
    _add_ratio(ax, bend_top, bend_bot, "(post-bend)", "--", color_post_bend)
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Splitting Ratio (Top / Total)")
ax.set_xlim(wl_min, wl_max)
ax.set_ylim(sr_min, sr_max)
handles = [Line2D([0], [0], color=color_pre_bend,  ls="-",  label="Pre-bend"),
           Line2D([0], [0], color=color_post_bend, ls="--", label="Post-bend")]
ax.legend(handles=handles if has_bend else handles[:1])
_style_ax(ax)
fig.tight_layout()
fig.savefig(f"{OUTPUT_FOLDER}/splitting_ratio.png", dpi=300)
plt.close(fig)
