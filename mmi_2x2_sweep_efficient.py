"""
Efficient mm_length sweep for MMI 2x2 (EME only).

Modes are computed once per wavelength from cross-sections that are
independent of mm_length. Each sweep value only recomputes the S-matrix
(cheap) with new cell lengths, not the eigenmodes (expensive).

Through port: load_2x2_MMI_top  (input top → output top)
Cross  port: load_2x2_MMI_bot  (input top → output bottom)
"""
import argparse
import dataconf
import json
import logging
import os
import pickle
from dataclasses import asdict, replace as dc_replace

import meep as mp
import meow as mw
import numpy as np
from matplotlib import pyplot as plt
from meow import compute_modes_tidy3d as compute_modes

from engine.helpers.general import make_output_folders
from engine.helpers.meow import load_2x2_MMI_bot, load_2x2_MMI_top, stack_extrusions_GF
from engine.hyperparams.MMI import HyperParamsMMI as HyperParams
from engine.plotting.meow import plot_eme_fields

MM_LENGTH_SWEEP = np.arange(20.0, 36.5, 0.5)

PORT_NAMES = ["through", "cross"]


def run_efficient_sweep(args):
    hp: HyperParams = dataconf.load(args.config, HyperParams)
    save_plots = args.save_plots
    base_out = hp.outputfolder
    sweep_dir = os.path.join(base_out, "sweep")
    os.makedirs(sweep_dir, exist_ok=True)

    if mp.am_master():
        handler = logging.FileHandler(os.path.join(sweep_dir, "sweep.log"))
        handler.setLevel("INFO")
        handler.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)-4s | %(filename)s | %(funcName)s:%(lineno)04d | %(message)s"
        ))
        logger = logging.getLogger("mmi_2x2_sweep_efficient")
        logger.addHandler(handler)
        logger.setLevel("INFO")
    else:
        logger = logging.getLogger("mmi_2x2_sweep_efficient")

    wavelengths = np.linspace(hp.sim.wavelength_min, hp.sim.wavelength_max, hp.sim.wavelength_pts)
    extrusion_rules = stack_extrusions_GF(hp.geo.waveguide_height, hp.sim.oxide_height)
    taper_Ls = np.array([hp.geo.taper_length / hp.sim.num_cells_taper] * hp.sim.num_cells_taper)
    num_modes_list = (
        [1]
        + [hp.sim.num_modes_taper] * hp.sim.num_cells_taper
        + [hp.sim.num_modes_MMI]
        + [hp.sim.num_modes_taper] * hp.sim.num_cells_taper
        + [1]
    )
    mesh = mw.Mesh2D(
        x=np.linspace(-hp.sim.oxide_width / 2, hp.sim.oxide_width / 2,
                      int(hp.sim.res * hp.sim.oxide_width) + 1),
        y=np.linspace(-hp.sim.oxide_height / 2, hp.sim.oxide_height / 2,
                      int(hp.sim.res * hp.sim.oxide_height) + 1),
    )

    # Reference cells — cross-sections are independent of mm_length.
    # Through and cross geometries share identical input cross-sections; only the
    # output waveguide position differs, so input_cells_ref is taken from one call.
    ref_mm_length = float(MM_LENGTH_SWEEP[0])
    Ls_ref = np.concatenate((
        [hp.sim.waveguide_length], taper_Ls, [ref_mm_length], taper_Ls, [hp.sim.waveguide_length]
    ))

    load_fns = [load_2x2_MMI_top, load_2x2_MMI_bot]
    input_cells_ref = None
    output_cells_ref = []

    for port_idx, load_fn in enumerate(load_fns):
        mmi_gds = load_fn(
            ref_mm_length, hp.sim.waveguide_length, hp.geo.waveguide_width,
            hp.geo.waveguide_offset, hp.sim.oxide_width, 0.0,
            hp.geo.taper_length, hp.geo.mm_width, hp.geo.taper_width,
        )
        structs = mw.extrude_gds(mmi_gds, extrusion_rules)
        cells = mw.create_cells(structures=structs, mesh=mesh, Ls=Ls_ref)
        if port_idx == 0:
            input_cells_ref = cells[0:-1]
        output_cells_ref.append(cells[-1])

    # Pre-pass: create output folders, save GDS layout plots, cache cells.
    sweep_vals = [float(v) for v in MM_LENGTH_SWEEP]
    cells_per_mm_val = {}

    for mm_val in sweep_vals:
        point_dir = os.path.join(sweep_dir, f"mmi_{mm_val:.2f}")
        if save_plots:
            make_output_folders(point_dir)

        Ls_new = np.concatenate((
            [hp.sim.waveguide_length], taper_Ls, [mm_val], taper_Ls, [hp.sim.waveguide_length]
        ))

        port_cells_new = []
        for port_idx, load_fn in enumerate(load_fns):
            mmi_gds_new = load_fn(
                mm_val, hp.sim.waveguide_length, hp.geo.waveguide_width,
                hp.geo.waveguide_offset, hp.sim.oxide_width, 0.0,
                hp.geo.taper_length, hp.geo.mm_width, hp.geo.taper_width,
            )
            if save_plots and mp.am_master():
                plt.figure()
                mmi_gds_new.plot()
                plt.savefig(
                    os.path.join(point_dir, "layout", f"eme_gds_{PORT_NAMES[port_idx]}.png"),
                    dpi=300,
                )
                plt.close("all")
            structs_new = mw.extrude_gds(mmi_gds_new, extrusion_rules)
            cells_new = mw.create_cells(structures=structs_new, mesh=mesh, Ls=Ls_new)
            if port_idx == 0:
                input_cells_new = cells_new[0:-1]
            port_cells_new.append(cells_new[-1])

        cells_per_mm_val[mm_val] = (input_cells_new, port_cells_new)

    # Results storage: sweep_val → port → list-over-wavelengths
    powers  = {v: [[] for _ in range(2)] for v in sweep_vals}
    angles  = {v: [[] for _ in range(2)] for v in sweep_vals}
    coeffs  = {v: [[] for _ in range(2)] for v in sweep_vals}

    center_idx = int(np.floor(len(wavelengths) / 2))

    for wl_idx, wavelength in enumerate(wavelengths):
        if mp.am_master():
            logger.info(f"Computing modes at wavelength {wavelength:.4f} μm (reused for all mm_length values)")

        env = mw.Environment(wl=float(wavelength), T=25.0)

        # Compute modes ONCE per wavelength from reference cross-sections.
        css_input = [mw.CrossSection.from_cell(cell=cell, env=env) for cell in input_cells_ref]
        input_modes = []
        for cell_num, (cs, n) in enumerate(zip(css_input, num_modes_list)):
            if mp.am_master():
                logger.info(f"  Input cell {cell_num}: computing {n} mode(s)")
            input_modes.append(compute_modes(cs, num_modes=n))

        output_modes = []
        for port_idx, out_cell in enumerate(output_cells_ref):
            cs = mw.CrossSection.from_cell(cell=out_cell, env=env)
            if mp.am_master():
                logger.info(f"  Output port {PORT_NAMES[port_idx]}: computing 1 mode")
            output_modes.append(compute_modes(cs, num_modes=1))

        if wl_idx == center_idx and save_plots and mp.am_master():
            cs_out_per_port = [
                mw.CrossSection.from_cell(cell=out_cell, env=env)
                for out_cell in output_cells_ref
            ]

        # Sweep mm_length: reuse modes, only recompute S-matrix.
        for mm_val in sweep_vals:
            input_cells_new, port_cells_new = cells_per_mm_val[mm_val]
            point_dir = os.path.join(sweep_dir, f"mmi_{mm_val:.2f}")

            if wl_idx == center_idx and save_plots and mp.am_master():
                for cs_num, cs in enumerate(css_input):
                    plt.clf()
                    plt.figure()
                    mw.visualize(cs, show=False)
                    plt.savefig(
                        os.path.join(point_dir, "layout", "cross_sections", f"cs_input_{cs_num}.png"),
                        dpi=300,
                    )
                    plt.close("all")
                for port_idx, cs in enumerate(cs_out_per_port):
                    plt.clf()
                    plt.figure()
                    mw.visualize(cs, show=False)
                    plt.savefig(
                        os.path.join(point_dir, "layout", "cross_sections", f"cs_out_{PORT_NAMES[port_idx]}.png"),
                        dpi=300,
                    )
                    plt.close("all")

            for port_idx in range(2):
                full_modes = input_modes + [output_modes[port_idx]]
                full_cells = input_cells_new + [port_cells_new[port_idx]]
                S, _ = mw.compute_s_matrix(full_modes, cells=full_cells)
                c = S[0, 1]
                coeffs[mm_val][port_idx].append(c)
                powers[mm_val][port_idx].append(float(np.abs(c) ** 2))
                angles[mm_val][port_idx].append(float(np.angle(c, deg=True)))
                if mp.am_master():
                    logger.info(
                        f"  mm_length={mm_val:.2f} port={PORT_NAMES[port_idx]} "
                        f"wl={wavelength:.4f}: P={np.abs(c)**2:.4f}"
                    )

                if wl_idx == center_idx and save_plots and mp.am_master():
                    plot_eme_fields(dc_replace(hp, outputfolder=point_dir), full_cells, full_modes, PORT_NAMES[port_idx])

    if mp.am_master():
        if save_plots:
            _save_results(sweep_dir, sweep_vals, wavelengths, hp, powers, angles, coeffs)
        _plot_summary(sweep_dir, sweep_vals, wavelengths, powers)


def _save_results(sweep_dir, sweep_vals, wavelengths, hp, powers, angles, coeffs):
    hp_dict = asdict(hp)
    for mm_val in sweep_vals:
        point_dir = os.path.join(sweep_dir, f"mmi_{mm_val:.2f}")

        hp_dict["geo"]["mm_length"] = mm_val
        hp_dict["outputfolder"] = point_dir
        with open(os.path.join(point_dir, "config.json"), "w") as f:
            print(json.dumps(hp_dict, indent=4, default=lambda x: x.name), file=f)

        with open(os.path.join(point_dir, "data", "sim_results.pkl"), "wb") as f:
            pickle.dump({
                "coeffs_ports": coeffs[mm_val],
                "powers_ports": powers[mm_val],
                "angles_ports": angles[mm_val],
                "wavelengths": wavelengths,
            }, f)

        plt.figure()
        for port_idx, name in enumerate(PORT_NAMES):
            plt.plot(wavelengths, powers[mm_val][port_idx], label=name.capitalize())
        plt.xlabel("Wavelength (μm)")
        plt.ylabel("Transmission")
        plt.title(f"MMI 2x2 — mm_length = {mm_val:.2f} μm")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(point_dir, "performance.png"), dpi=300)
        plt.close("all")


def _plot_summary(sweep_dir, sweep_vals, wavelengths, powers):
    multi_wl = len(wavelengths) > 1
    fig, ax = plt.subplots()

    for port_idx, name in enumerate(PORT_NAMES):
        vals = np.array([powers[v][port_idx] for v in sweep_vals])  # (n_sweep, n_wl)
        if multi_wl:
            mean = vals.mean(axis=1)
            lo   = vals.min(axis=1)
            hi   = vals.max(axis=1)
            (line,) = ax.plot(sweep_vals, mean, marker="o", label=name.capitalize())
            ax.fill_between(sweep_vals, lo, hi, alpha=0.25, color=line.get_color())
        else:
            ax.plot(sweep_vals, vals[:, 0], marker="o", label=name.capitalize())

    ax.set_xlabel("mm_length (μm)")
    if multi_wl:
        ax.set_ylabel(
            f"Transmission (λ = {wavelengths[0]:.3f}–{wavelengths[-1]:.3f} μm)\n"
            "solid = mean, band = [min, max]"
        )
    else:
        ax.set_ylabel(f"Transmission at λ={wavelengths[0]:.3f} μm")
    ax.set_title("MMI 2x2 — mm_length sweep (EME)")
    ax.legend()
    fig.tight_layout()
    out_path = os.path.join(sweep_dir, "sweep_summary.png")
    fig.savefig(out_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Summary plot saved to {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file (sim_type must be EME)", required=True)
    parser.add_argument("--save-plots", action="store_true",
                        help="save layout/GDS plots, cross-sections, EME fields, and per-point performance plots")
    args = parser.parse_args()
    run_efficient_sweep(args)
