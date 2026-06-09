"""
Efficient mm_length sweep for MMI 2x2 (hybrid simulation).

FDTD (top_simulation, bot_simulation) runs once — both are independent of
mm_length. Modes are computed once per wavelength from cross-sections that are
independent of mm_length (constant-width MMI rectangle). Each sweep value only
recomputes the S-matrix (cheap) with new cell lengths.
"""
import argparse
import dataconf
import json
import logging
import os
import pickle
from dataclasses import asdict

import meep as mp
import meow as mw
import numpy as np
from matplotlib import pyplot as plt

from engine.helpers.general import make_output_folders
from engine.helpers.meow import load_2x2_MMI_top, stack_extrusions_GF
from engine.hyperparams.MMI import HyperParamsMMI as HyperParams
from mmi_2x2_master import MMI_2x2

MM_LENGTH_SWEEP = np.arange(20.0, 36.5, 0.5)


def run_sweep(args):
    hp: HyperParams = dataconf.load(args.config, HyperParams)
    sweep_dir = os.path.join(hp.outputfolder, "sweep")
    os.makedirs(sweep_dir, exist_ok=True)

    if mp.am_master():
        handler = logging.FileHandler(os.path.join(sweep_dir, "sweep.log"))
        handler.setLevel("INFO")
        handler.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)-4s | %(filename)s | %(funcName)s:%(lineno)04d | %(message)s"
        ))
        logger = logging.getLogger("mmi_2x2_sweep")
        logger.addHandler(handler)
        logger.setLevel("INFO")
    else:
        logger = logging.getLogger("mmi_2x2_sweep")

    wavelengths = np.linspace(hp.sim.wavelength_min, hp.sim.wavelength_max, hp.sim.wavelength_pts)
    sweep_vals = [float(v) for v in MM_LENGTH_SWEEP]

    # --- Phase 1: Run FDTD simulations once (mm_length-independent) ---
    if mp.am_master(): logger.info("Running top and bottom FDTD simulations (once for all sweep values)")
    sim = MMI_2x2(args)
    top_coeffs = sim.top_simulation(wavelengths)  # (n_wl, n_modes)
    bot_coeffs = sim.bot_simulation(wavelengths)  # (n_wl, n_modes)

    # --- Phase 2: Build reference structs (cross-sections do not vary with mm_length) ---
    if mp.am_master(): logger.info("Building reference EME cells")
    mmi_gds = load_2x2_MMI_top(
        hp.geo.mm_length, 0.0, hp.geo.waveguide_width,
        hp.geo.waveguide_offset, hp.sim.oxide_width,
        -hp.geo.taper_length - hp.sim.mm_monitor_position,
        hp.geo.taper_length, hp.geo.mm_width, hp.geo.taper_width,
    )
    extrusion_rules = stack_extrusions_GF(hp.geo.waveguide_height, hp.sim.oxide_height)
    structs = mw.extrude_gds(mmi_gds, extrusion_rules)
    mesh = mw.Mesh2D(
        x=np.linspace(-hp.sim.oxide_width / 2, hp.sim.oxide_width / 2,
                      int(hp.sim.res * hp.sim.oxide_width) + 1),
        y=np.linspace(-hp.sim.oxide_height / 2, hp.sim.oxide_height / 2,
                      int(hp.sim.res * hp.sim.oxide_height) + 1),
    )
    ref_Ls = [hp.geo.mm_length - 2 * hp.sim.mm_monitor_position]
    cells_ref = mw.create_cells(structures=structs, mesh=mesh, Ls=ref_Ls)

    for mm_val in sweep_vals:
        make_output_folders(os.path.join(sweep_dir, f"mmi_{mm_val:.2f}"))

    powers_top    = {v: [] for v in sweep_vals}
    powers_bot    = {v: [] for v in sweep_vals}
    angles_top    = {v: [] for v in sweep_vals}
    angles_bot    = {v: [] for v in sweep_vals}
    mm_coeffs_all = {v: [] for v in sweep_vals}

    # --- Phase 3: Per-wavelength mode solve (once), then cheap S-matrix per sweep value ---
    for wl_idx, wavelength in enumerate(wavelengths):
        if mp.am_master(): logger.info(
            f"Computing modes at wavelength {wavelength:.4f} μm (reused for all mm_length values)"
        )
        env = mw.Environment(wl=float(wavelength), T=25.0)
        css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells_ref]
        modes = [mw.compute_modes_meep(cs, num_modes=hp.sim.num_modes_MMI) for cs in css]

        top_coeff = top_coeffs[wl_idx]
        bot_coeff = bot_coeffs[wl_idx]

        for mm_val in sweep_vals:
            cells_new = mw.create_cells(
                structures=structs, mesh=mesh,
                Ls=[mm_val - 2 * hp.sim.mm_monitor_position],
            )
            S, _ = mw.compute_s_matrix(modes, cells=cells_new)

            mm_coeffs_1_wl = [S[m, hp.sim.num_modes_MMI + m] for m in range(hp.sim.num_modes_MMI)]

            propagated = [x * y for x, y in zip(top_coeff, mm_coeffs_1_wl)]
            out_top = sum(x * y for x, y in zip(propagated, top_coeff))
            out_bot = sum(x * y for x, y in zip(propagated, bot_coeff))

            powers_top[mm_val].append(float(np.abs(out_top) ** 2))
            powers_bot[mm_val].append(float(np.abs(out_bot) ** 2))
            angles_top[mm_val].append(float(np.angle(out_top, deg=True)))
            angles_bot[mm_val].append(float(np.angle(out_bot, deg=True)))
            mm_coeffs_all[mm_val].append(mm_coeffs_1_wl)

            if mp.am_master():
                logger.info(
                    f"  mm_length={mm_val:.2f} wl={wavelength:.4f}: "
                    f"P_top={powers_top[mm_val][-1]:.4f} P_bot={powers_bot[mm_val][-1]:.4f}"
                )

    # --- Phase 4: Save per-point results and summary plot ---
    if mp.am_master():
        hp_dict = asdict(hp)
        for mm_val in sweep_vals:
            point_dir = os.path.join(sweep_dir, f"mmi_{mm_val:.2f}")

            hp_dict["geo"]["mm_length"] = mm_val
            hp_dict["outputfolder"] = point_dir
            with open(os.path.join(point_dir, "config.json"), "w") as f:
                print(json.dumps(hp_dict, indent=4, default=lambda x: x.name), file=f)

            with open(os.path.join(point_dir, "data", "sim_results.pkl"), "wb") as f:
                pickle.dump({
                    "top_coeffs": top_coeffs,
                    "mm_coeffs": mm_coeffs_all[mm_val],
                    "bot_coeffs": bot_coeffs,
                    "powers_top": powers_top[mm_val],
                    "powers_bot": powers_bot[mm_val],
                    "angles_top": angles_top[mm_val],
                    "angles_bot": angles_bot[mm_val],
                    "wavelengths": wavelengths,
                }, f)

            plt.figure()
            plt.plot(wavelengths, powers_top[mm_val], label="Top")
            plt.plot(wavelengths, powers_bot[mm_val], label="Bottom")
            plt.xlabel("Wavelength (μm)")
            plt.ylabel("Transmission")
            plt.title(f"MMI 2x2 — mm_length = {mm_val:.2f} μm")
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(point_dir, "performance.png"), dpi=300)
            plt.close("all")

        _plot_summary(sweep_dir, sweep_vals, wavelengths, powers_top, powers_bot)


def _plot_summary(sweep_dir, sweep_vals, wavelengths, powers_top, powers_bot):
    multi_wl = len(wavelengths) > 1
    fig, ax = plt.subplots()

    for label, powers in [("Top", powers_top), ("Bottom", powers_bot)]:
        vals = np.array([powers[v] for v in sweep_vals])  # (n_sweep, n_wl)
        if multi_wl:
            mean = vals.mean(axis=1)
            lo   = vals.min(axis=1)
            hi   = vals.max(axis=1)
            (line,) = ax.plot(sweep_vals, mean, marker="o", label=label)
            ax.fill_between(sweep_vals, lo, hi, alpha=0.25, color=line.get_color())
        else:
            ax.plot(sweep_vals, vals[:, 0], marker="o", label=label)

    ax.set_xlabel("mm_length (μm)")
    if multi_wl:
        ax.set_ylabel(
            f"Transmission (λ = {wavelengths[0]:.3f}–{wavelengths[-1]:.3f} μm)\n"
            "solid = mean, band = [min, max]"
        )
    else:
        ax.set_ylabel(f"Transmission at λ={wavelengths[0]:.3f} μm")
    ax.set_title("MMI 2x2 — mm_length sweep")
    ax.legend()
    fig.tight_layout()
    out_path = os.path.join(sweep_dir, "sweep_summary.png")
    fig.savefig(out_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Summary plot saved to {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file (sim_type should be HYBRID)", required=True)
    args = parser.parse_args()
    run_sweep(args)
