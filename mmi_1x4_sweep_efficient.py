"""
Efficient mm_length sweep for MMI 1x4 (EME only).

Modes are computed once per wavelength from cross-sections that are
independent of mm_length. Each sweep value only recomputes the S-matrix
(cheap) with new cell lengths, not the eigenmodes (expensive).
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
from engine.helpers.meow import load_1x4_MMI, stack_extrusions_GF
from engine.hyperparams.MMI_1x4 import HyperParamsMMI as HyperParams
from engine.plotting.meow import plot_eme_fields

MM_LENGTH_SWEEP = np.arange(5.0, 12.0, 0.1)  # 5.0–15.0 μm, 21 points


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
        logger = logging.getLogger("mmi_1x4_sweep_efficient")
        logger.addHandler(handler)
        logger.setLevel("INFO")
    else:
        logger = logging.getLogger("mmi_1x4_sweep_efficient")

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

    # Reference cells — cross-sections are identical for all mm_length values since
    # the cross-sectional geometry (WG width/height, taper profile, MMI width) is
    # independent of mm_length. Modes computed from these cells are reused for every
    # sweep point.
    ref_mm_length = float(MM_LENGTH_SWEEP[0])
    Ls_ref = np.concatenate((
        [hp.sim.waveguide_length], taper_Ls, [ref_mm_length], taper_Ls, [hp.sim.waveguide_length]
    ))

    input_cells_ref = None
    output_cells_ref = []
    for output_idx in range(4):
        mmi_gds = load_1x4_MMI(
            ref_mm_length, hp.sim.waveguide_length, hp.geo.waveguide_width,
            hp.geo.waveguide_offset, hp.sim.oxide_width, 0.0,
            hp.geo.taper_length, hp.geo.mm_width, hp.geo.taper_width, output_idx
        )
        structs = mw.extrude_gds(mmi_gds, extrusion_rules)
        cells = mw.create_cells(structures=structs, mesh=mesh, Ls=Ls_ref)
        if output_idx == 0:
            input_cells_ref = cells[0:-1]
        output_cells_ref.append(cells[-1])

    # Pre-pass: create output folders, save GDS layout plots, and build cells for
    # each sweep value. Separating this from the wavelength loop avoids redundant
    # GDS I/O and makes center-wavelength layout/field saving straightforward.
    sweep_vals = [float(v) for v in MM_LENGTH_SWEEP]
    cells_per_mm_val = {}

    for mm_val in sweep_vals:
        point_dir = os.path.join(sweep_dir, f"mmi_{mm_val:.2f}")
        if save_plots:
            make_output_folders(point_dir)

        Ls_new = np.concatenate((
            [hp.sim.waveguide_length], taper_Ls, [mm_val], taper_Ls, [hp.sim.waveguide_length]
        ))

        input_cells_new = None
        port_range = range(4) if save_plots else range(1)
        for output_idx in port_range:
            mmi_gds_new = load_1x4_MMI(
                mm_val, hp.sim.waveguide_length, hp.geo.waveguide_width,
                hp.geo.waveguide_offset, hp.sim.oxide_width, 0.0,
                hp.geo.taper_length, hp.geo.mm_width, hp.geo.taper_width, output_idx
            )
            if save_plots and mp.am_master():
                plt.figure()
                mmi_gds_new.plot()
                plt.savefig(os.path.join(point_dir, "layout", f"eme_gds_port{output_idx}.png"), dpi=300)
                plt.close("all")
            if output_idx == 0:
                structs_new = mw.extrude_gds(mmi_gds_new, extrusion_rules)
                cells_new = mw.create_cells(structures=structs_new, mesh=mesh, Ls=Ls_new)
                input_cells_new = cells_new[0:-1]
                output_cell_new = cells_new[-1]

        cells_per_mm_val[mm_val] = (input_cells_new, output_cell_new)

    # Results storage: sweep_val → port → list-over-wavelengths
    powers = {v: [[] for _ in range(4)] for v in sweep_vals}
    angles = {v: [[] for _ in range(4)] for v in sweep_vals}
    coeffs = {v: [[] for _ in range(4)] for v in sweep_vals}

    center_idx = int(np.floor(len(wavelengths) / 2))

    for wl_idx, wavelength in enumerate(wavelengths):
        if mp.am_master(): logger.info(f"Computing modes at wavelength {wavelength:.4f} μm (reused for all mm_length values)")

        env = mw.Environment(wl=float(wavelength), T=25.0)

        # --- Compute modes ONCE per wavelength ---
        css_input = [mw.CrossSection.from_cell(cell=cell, env=env) for cell in input_cells_ref]
        input_modes = []
        for cell_num, (cs, n) in enumerate(zip(css_input, num_modes_list)):
            if mp.am_master(): logger.info(f"  Input cell {cell_num}: computing {n} mode(s)")
            input_modes.append(compute_modes(cs, num_modes=n))

        output_modes = []
        for output_idx, out_cell in enumerate(output_cells_ref):
            cs = mw.CrossSection.from_cell(cell=out_cell, env=env)
            if mp.am_master(): logger.info(f"  Output port {output_idx}: computing 1 mode")
            output_modes.append(compute_modes(cs, num_modes=1))

        if wl_idx == center_idx and save_plots and mp.am_master():
            cs_out_per_port = [
                mw.CrossSection.from_cell(cell=out_cell, env=env)
                for out_cell in output_cells_ref
            ]

        # --- Sweep mm_length: use pre-built cells, reuse modes ---
        for mm_val in sweep_vals:
            input_cells_new, output_cell_new = cells_per_mm_val[mm_val]
            point_dir = os.path.join(sweep_dir, f"mmi_{mm_val:.2f}")

            if wl_idx == center_idx and save_plots and mp.am_master():
                for cs_num, cs in enumerate(css_input):
                    plt.clf()
                    plt.figure()
                    mw.visualize(cs, show=False)
                    plt.savefig(os.path.join(point_dir, "layout", "cross_sections", f"cs_input_{cs_num}.png"), dpi=300)
                    plt.close("all")
                for out_idx, cs in enumerate(cs_out_per_port):
                    plt.clf()
                    plt.figure()
                    mw.visualize(cs, show=False)
                    plt.savefig(os.path.join(point_dir, "layout", "cross_sections", f"cs_out_{out_idx}_0.png"), dpi=300)
                    plt.close("all")

            for output_idx in range(4):
                full_modes = input_modes + [output_modes[output_idx]]
                full_cells = input_cells_new + [output_cell_new]
                S, _ = mw.compute_s_matrix(full_modes, cells=full_cells)
                c = S[0, 1]
                coeffs[mm_val][output_idx].append(c)
                powers[mm_val][output_idx].append(float(np.abs(c) ** 2))
                angles[mm_val][output_idx].append(float(np.angle(c, deg=True)))
                if mp.am_master():
                    logger.info(f"  mm_length={mm_val:.2f} port={output_idx} wl={wavelength:.4f}: P={np.abs(c)**2:.4f}")

                if wl_idx == center_idx and save_plots and mp.am_master():
                    plot_eme_fields(dc_replace(hp, outputfolder=point_dir), full_cells, full_modes, output_idx)

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
        for port_idx, port_powers in enumerate(powers[mm_val]):
            plt.plot(wavelengths, port_powers, label=f"Port {port_idx}")
        plt.xlabel("Wavelength (μm)")
        plt.ylabel("Transmission")
        plt.title(f"MMI 1x4 — mm_length = {mm_val:.2f} μm")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(point_dir, "performance.png"), dpi=300)
        plt.close("all")


def _plot_summary(sweep_dir, sweep_vals, wavelengths, powers):
    center_idx = int(np.floor(len(wavelengths) / 2))
    plt.figure()
    for port_idx in range(4):
        port_vals = [powers[v][port_idx][center_idx] for v in sweep_vals]
        plt.plot(sweep_vals, port_vals, marker="o", label=f"Port {port_idx}")
    plt.xlabel("mm_length (μm)")
    plt.ylabel(f"Transmission at λ={wavelengths[center_idx]:.3f} μm")
    plt.title("MMI 1x4 — mm_length sweep")
    plt.legend()
    plt.tight_layout()
    out_path = os.path.join(sweep_dir, "sweep_summary.png")
    plt.savefig(out_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Summary plot saved to {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file (sim_type must be EME)", required=True)
    parser.add_argument("--save-plots", action="store_true",
                        help="save layout/GDS plots, cross-sections, EME fields, and per-point performance plots")
    args = parser.parse_args()
    run_efficient_sweep(args)
