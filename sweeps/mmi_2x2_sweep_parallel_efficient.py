"""
Efficient parallel mm_length x mm_width sweep for MMI 2x2 (EME only).

Cell stack: [input_wg] + [input_taper]*N + [MMI] + [output_taper]*N + [output_wg].

Cross-sections independent of BOTH sweep params:
  input_wg, input_tapers  -> shared left cells (main process, once per wavelength)
  output_tapers, output_wg (through & cross) -> shared right cells (main process, once per wavelength)

Cross-section dependent on mm_width only (not mm_length):
  MMI body -> each worker computes these once per (wavelength, mm_width)

S-matrix recomputed per (mm_length, mm_width) -- no mode solving, just propagation.

Two ports monitored: through (top input -> top output) and cross (top input -> bottom output).
"""
import sys as _sys, os as _os
_FECES_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _FECES_ROOT)
_sys.path.insert(0, _os.path.join(_FECES_ROOT, "masters"))
_os.environ.setdefault("OMP_NUM_THREADS", "1")
del _sys, _os

import argparse
import dataconf
import gc
import json
import logging
import os
import pickle
import tempfile
from dataclasses import asdict, replace as dc_replace
from multiprocessing import Pool

import meep as mp
import meow as mw
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from meow import compute_modes_tidy3d as compute_modes

from engine.helpers.general import make_output_folders
from engine.helpers.meow import load_2x2_MMI_bot, load_2x2_MMI_top, stack_extrusions_GF
from engine.hyperparams.MMI import HyperParamsMMI as HyperParams
from engine.plotting.meow import plot_eme_fields

N_PROCESSES = 4
WL_BATCH_SIZE = 1  # max wavelengths whose MMI modes are held in memory at once per worker
MIN_LW = 0.2

MM_LENGTH_SWEEP = np.linspace(20.0, 60.0, 41)   # 20:1:36 um
MM_WIDTH_SWEEP  = np.linspace(5.0, 7.0, 21)      # 5.0:0.2:7.0 um

PORT_NAMES = ["through", "cross"]
LOAD_FNS   = [load_2x2_MMI_top, load_2x2_MMI_bot]


def _build_cells(hp, mesh, extrusion_rules, taper_Ls, length, width, load_fn):
    taper_width = (width - MIN_LW)/2
    waveguide_offset = taper_width + MIN_LW
    mmi_gds = load_fn(
        length, hp.sim.waveguide_length, hp.geo.waveguide_width,
        waveguide_offset, hp.sim.oxide_width, 0.0,
        hp.geo.taper_length, width, taper_width,
    )
    structs = mw.extrude_gds(mmi_gds, extrusion_rules)
    Ls = np.concatenate((
        [hp.sim.waveguide_length], taper_Ls, [length], taper_Ls, [hp.sim.waveguide_length]
    ))
    return mw.create_cells(structures=structs, mesh=mesh, Ls=Ls)


def _init_worker(hp, mesh, extrusion_rules, taper_Ls,
                  mmi_num_modes, right_num_modes,
                  split_left, split_right,
                  shared_modes_path,
                  wavelengths, sweep_dir, save_plots, full_output):
    global _hp, _mesh, _extrusion_rules, _taper_Ls
    global _mmi_num_modes, _right_num_modes
    global _split_left, _split_right
    global _shared_left_modes, _shared_right_modes
    global _wavelengths, _sweep_dir, _save_plots, _full_output
    _hp = hp
    _mesh = mesh
    _extrusion_rules = extrusion_rules
    _taper_Ls = taper_Ls
    _mmi_num_modes = mmi_num_modes
    _right_num_modes = right_num_modes
    _split_left = split_left
    _split_right = split_right
    # Load shared modes from disk once per worker -- avoids pickling large arrays
    # N_PROCESSES times through initargs.
    with open(shared_modes_path, "rb") as f:
        _d = pickle.load(f)
    _shared_left_modes  = _d["left"]
    _shared_right_modes = _d["right"]
    del _d
    _wavelengths = wavelengths
    _sweep_dir = sweep_dir
    _save_plots = save_plots
    _full_output = full_output


def _run_width(mm_width):
    print(f"[sweep] worker starting mm_width = {mm_width:.3f} um")
    center_idx = int(np.floor(len(_wavelengths) / 2))
    center_length = float(MM_LENGTH_SWEEP[len(MM_LENGTH_SWEEP) // 2])
    ref_length = float(MM_LENGTH_SWEEP[0])

    ref_cells_top = _build_cells(_hp, _mesh, _extrusion_rules, _taper_Ls,
                                  ref_length, mm_width, load_2x2_MMI_top)
    mmi_cell_ref = ref_cells_top[_split_left]
    del ref_cells_top

    # Split wavelength indices into batches so only WL_BATCH_SIZE sets of MMI modes
    # are held in memory at once instead of all n_wavelengths simultaneously.
    all_wl = list(range(len(_wavelengths)))
    batches = [all_wl[i:i + WL_BATCH_SIZE] for i in range(0, len(all_wl), WL_BATCH_SIZE)]

    # Accumulators: wavelength results arrive per-batch and are appended in order,
    # so acc_powers[length][port][k] == power at global wavelength index k.
    acc_powers = {float(l): [[] for _ in range(2)] for l in MM_LENGTH_SWEEP}
    acc_angles = {float(l): [[] for _ in range(2)] for l in MM_LENGTH_SWEEP}

    for batch in batches:
        # Solve MMI modes for this batch only -- freed at end of each batch.
        mmi_modes_batch = {}
        for wl_idx in batch:
            env = mw.Environment(wl=float(_wavelengths[wl_idx]), T=25.0)
            cs = mw.CrossSection.from_cell(cell=mmi_cell_ref, env=env)
            mmi_modes_batch[wl_idx] = compute_modes(cs, num_modes=_mmi_num_modes)

        for length in MM_LENGTH_SWEEP:
            length = float(length)
            cells_per_port = [
                _build_cells(_hp, _mesh, _extrusion_rules, _taper_Ls, length, mm_width, fn)
                for fn in LOAD_FNS
            ]

            for wl_idx in batch:
                full_left = _shared_left_modes[wl_idx]
                full_mmi  = [mmi_modes_batch[wl_idx]]
                for port_idx in range(2):
                    full_right = _shared_right_modes[port_idx][wl_idx]
                    full_modes = full_left + full_mmi + full_right
                    S, _ = mw.compute_s_matrix(full_modes, cells=cells_per_port[port_idx])
                    c = S[0, 1]
                    acc_powers[length][port_idx].append(float(np.abs(c) ** 2))
                    acc_angles[length][port_idx].append(float(np.angle(c, deg=True)))
                print(
                    f"  mm_length={length:.2f} mm_width={mm_width:.3f} "
                    f"wl={_wavelengths[wl_idx]:.4f}: "
                    f"T={acc_powers[length][0][-1]:.4f} X={acc_powers[length][1][-1]:.4f}"
                )

            del cells_per_port
            gc.collect()

        del mmi_modes_batch
        gc.collect()

    # All batches done -- full wavelength arrays now in accumulators.
    # Solve center-wavelength MMI mode once per mm_width (not per length) for layout/field plots.
    mmi_mode_center = None
    env_c = None
    if _full_output:
        env_c = mw.Environment(wl=float(_wavelengths[center_idx]), T=25.0)
        cs_c = mw.CrossSection.from_cell(cell=mmi_cell_ref, env=env_c)
        mmi_mode_center = compute_modes(cs_c, num_modes=_mmi_num_modes)

    # Geometry derived from mm_width -- consistent with _build_cells.
    taper_width_w     = (mm_width - MIN_LW) / 2
    waveguide_offset_w = taper_width_w + MIN_LW

    result = {}
    for length in MM_LENGTH_SWEEP:
        length = float(length)
        powers_ports = acc_powers[length]
        angles_ports = acc_angles[length]

        result[length] = {
            "center_powers": [powers_ports[p][center_idx] for p in range(2)],
        }

        if _save_plots or _full_output:
            point_dir = os.path.join(_sweep_dir, f"len_{length:.2f}_wid_{mm_width:.3f}")
            make_output_folders(point_dir)

            hp_dict = asdict(_hp)
            hp_dict["geo"]["mm_length"] = length
            hp_dict["geo"]["mm_width"] = mm_width
            hp_dict["outputfolder"] = point_dir
            with open(os.path.join(point_dir, "config.json"), "w") as f:
                print(json.dumps(hp_dict, indent=4, default=lambda x: x.name), file=f)

            with open(os.path.join(point_dir, "data", "sim_results.pkl"), "wb") as f:
                pickle.dump({
                    "powers_ports": powers_ports,
                    "angles_ports": angles_ports,
                    "wavelengths": _wavelengths,
                }, f)

            plt.figure()
            for port_idx, name in enumerate(PORT_NAMES):
                plt.plot(_wavelengths, powers_ports[port_idx], label=name.capitalize())
            plt.xlabel("Wavelength (um)")
            plt.ylabel("Transmission")
            plt.title(f"MMI 2x2 -- mm_length={length:.2f} um, mm_width={mm_width:.3f} um")
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(point_dir, "performance.png"), dpi=300)
            plt.close("all")

            if _full_output:
                for fn, port_name in zip(LOAD_FNS, PORT_NAMES):
                    mmi_gds = fn(
                        length, _hp.sim.waveguide_length, _hp.geo.waveguide_width,
                        waveguide_offset_w, _hp.sim.oxide_width, 0.0,
                        _hp.geo.taper_length, mm_width, taper_width_w,
                    )
                    plt.figure()
                    mmi_gds.plot()
                    plt.savefig(
                        os.path.join(point_dir, "layout", f"eme_gds_{port_name}.png"), dpi=300
                    )
                    plt.close("all")

                cells_top = _build_cells(_hp, _mesh, _extrusion_rules, _taper_Ls,
                                         length, mm_width, LOAD_FNS[0])
                css_full = [mw.CrossSection.from_cell(cell=c, env=env_c) for c in cells_top]
                for cs_num, cs in enumerate(css_full):
                    plt.clf()
                    plt.figure()
                    mw.visualize(cs, show=False)
                    plt.savefig(
                        os.path.join(point_dir, "layout", "cross_sections", f"cs_{cs_num}.png"),
                        dpi=300,
                    )
                    plt.close("all")
                del css_full

                cells_per_port = [
                    _build_cells(_hp, _mesh, _extrusion_rules, _taper_Ls, length, mm_width, fn)
                    for fn in LOAD_FNS
                ]
                full_left_c = _shared_left_modes[center_idx]
                full_mmi_c  = [mmi_mode_center]
                for port_idx in range(2):
                    full_right_c = _shared_right_modes[port_idx][center_idx]
                    full_modes_c = full_left_c + full_mmi_c + full_right_c
                    plot_eme_fields(
                        dc_replace(_hp, outputfolder=point_dir),
                        cells_per_port[port_idx], full_modes_c, PORT_NAMES[port_idx],
                    )
                del cells_top, cells_per_port, full_left_c, full_mmi_c, full_modes_c
                gc.collect()

    if mmi_mode_center is not None:
        del mmi_mode_center
        gc.collect()

    return (float(mm_width), result)


def run_efficient_sweep(args):
    hp: HyperParams = dataconf.load(args.config, HyperParams)
    sweep_dir = os.path.join(hp.outputfolder, "sweep")
    os.makedirs(sweep_dir, exist_ok=True)

    if mp.am_master():
        handler = logging.FileHandler(os.path.join(sweep_dir, "sweep.log"))
        handler.setLevel("INFO")
        handler.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)-4s | %(filename)s | %(funcName)s:%(lineno)04d | %(message)s"
        ))
        logger = logging.getLogger("mmi_2x2_sweep_parallel_efficient")
        logger.addHandler(handler)
        logger.setLevel("INFO")
    else:
        logger = logging.getLogger("mmi_2x2_sweep_parallel_efficient")

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

    # split_left: index of first non-left cell (= input_wg + input_tapers)
    # split_right: index of first right cell (= after MMI body cell)
    split_left  = 1 + hp.sim.num_cells_taper
    split_right = split_left + 1
    mmi_num_modes   = hp.sim.num_modes_MMI
    right_num_modes = num_modes_list[split_right:]

    assert len(num_modes_list) == 2 * hp.sim.num_cells_taper + 3

    # Build reference cells for main-process shared-mode extraction.
    ref_length = float(MM_LENGTH_SWEEP[0])
    ref_width  = float(MM_WIDTH_SWEEP[0])

    # Left cells are identical across all (length, width) pairs -- extract from top GDS.
    ref_cells_top = _build_cells(hp, mesh, extrusion_rules, taper_Ls,
                                  ref_length, ref_width, load_2x2_MMI_top)
    left_cells_ref = ref_cells_top[:split_left]
    left_num_modes = num_modes_list[:split_left]

    # Right cells differ between ports but not across (length, width).
    right_cells_ref = []
    for fn in LOAD_FNS:
        ref_cells = _build_cells(hp, mesh, extrusion_rules, taper_Ls,
                                  ref_length, ref_width, fn)
        right_cells_ref.append(ref_cells[split_right:])

    if mp.am_master():
        logger.info("Computing shared left and right cell modes once per wavelength")

    shared_left_modes  = {}
    shared_right_modes = [{} for _ in range(2)]

    for wl_idx, wavelength in enumerate(wavelengths):
        env = mw.Environment(wl=float(wavelength), T=25.0)

        css_left = [mw.CrossSection.from_cell(cell=c, env=env) for c in left_cells_ref]
        shared_left_modes[wl_idx] = [
            compute_modes(cs, num_modes=n) for cs, n in zip(css_left, left_num_modes)
        ]

        for port_idx in range(2):
            css_right = [mw.CrossSection.from_cell(cell=c, env=env)
                         for c in right_cells_ref[port_idx]]
            shared_right_modes[port_idx][wl_idx] = [
                compute_modes(cs, num_modes=n)
                for cs, n in zip(css_right, right_num_modes)
            ]

    # Write shared modes to a temp file so workers load from disk instead of
    # receiving N_PROCESSES pickled copies through initargs.
    tmp = tempfile.NamedTemporaryFile(suffix=".pkl", delete=False)
    try:
        pickle.dump({"left": shared_left_modes, "right": shared_right_modes}, tmp)
        tmp.close()
        shared_modes_path = tmp.name

        if mp.am_master():
            logger.info(f"Dispatching {len(MM_WIDTH_SWEEP)} widths across {N_PROCESSES} processes")

        Z = {name: np.full((len(MM_LENGTH_SWEEP), len(MM_WIDTH_SWEEP)), np.nan)
             for name in PORT_NAMES}
        all_results = {}

        with Pool(
            processes=N_PROCESSES,
            maxtasksperchild=1,
            initializer=_init_worker,
            initargs=(hp, mesh, extrusion_rules, taper_Ls,
                      mmi_num_modes, right_num_modes,
                      split_left, split_right,
                      shared_modes_path,
                      wavelengths, sweep_dir, args.save_plots, args.full_output),
        ) as pool:
            for mm_width, length_dict in pool.imap_unordered(
                _run_width, [float(w) for w in MM_WIDTH_SWEEP]
            ):
                j = int(np.argmin(np.abs(MM_WIDTH_SWEEP - mm_width)))
                all_results[mm_width] = length_dict
                for length, d in length_dict.items():
                    i = int(np.argmin(np.abs(MM_LENGTH_SWEEP - length)))
                    for port_idx, name in enumerate(PORT_NAMES):
                        Z[name][i, j] = d["center_powers"][port_idx]

    finally:
        os.unlink(tmp.name)

    if mp.am_master():
        _plot_summary(sweep_dir, wavelengths, Z, all_results)


def _plot_summary(sweep_dir, wavelengths, Z, all_results):
    center_idx = int(np.floor(len(wavelengths) / 2))
    wl_c = wavelengths[center_idx]

    with open(os.path.join(sweep_dir, "sweep_results.pkl"), "wb") as f:
        pickle.dump({
            "mm_length":    MM_LENGTH_SWEEP,
            "mm_width":     MM_WIDTH_SWEEP,
            "transmission": Z,
            "wavelength_center": wl_c,
            "all_results":  all_results,
        }, f)

    X, Y = np.meshgrid(MM_WIDTH_SWEEP, MM_LENGTH_SWEEP)

    for name in PORT_NAMES:
        # Heatmap
        plt.figure()
        mesh_plot = plt.pcolormesh(X, Y, Z[name], shading="auto", cmap="viridis")
        plt.colorbar(mesh_plot, label=f"{name.capitalize()} transmission @ {wl_c:.3f} um")
        plt.xlabel("mm_width (um)")
        plt.ylabel("mm_length (um)")
        plt.title(f"MMI 2x2 -- {name} port -- length / width sweep (parallel, efficient)")
        plt.tight_layout()
        heatmap_path = os.path.join(sweep_dir, f"sweep_heatmap_{name}.png")
        plt.savefig(heatmap_path, dpi=300)
        plt.close("all")
        print(f"[sweep] Heatmap ({name}) saved to {heatmap_path}")

        # Surface plot
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        ax.plot_surface(X, Y, Z[name], cmap="viridis")
        ax.set_xlabel("mm_width (um)")
        ax.set_ylabel("mm_length (um)")
        ax.set_zlabel(f"{name.capitalize()} @ {wl_c:.3f} um")
        ax.set_title(f"MMI 2x2 -- {name} port -- length / width sweep")
        surface_path = os.path.join(sweep_dir, f"sweep_surface_{name}.png")
        plt.savefig(surface_path, dpi=300)
        plt.close("all")
        print(f"[sweep] Surface ({name}) saved to {surface_path}")

    # Total transmission: through + cross
    Z_total = Z["through"] + Z["cross"]

    plt.figure()
    mesh_plot = plt.pcolormesh(X, Y, Z_total, shading="auto", cmap="viridis",
                               vmin=0.0, vmax=1.0)
    plt.colorbar(mesh_plot, label=f"Through + Cross @ {wl_c:.3f} um")
    plt.xlabel("mm_width (um)")
    plt.ylabel("mm_length (um)")
    plt.title("MMI 2x2 -- total transmission (through + cross)")
    plt.tight_layout()
    total_heatmap_path = os.path.join(sweep_dir, "sweep_heatmap_total.png")
    plt.savefig(total_heatmap_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Total-transmission heatmap saved to {total_heatmap_path}")

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(X, Y, Z_total, cmap="viridis")
    ax.set_xlabel("mm_width (um)")
    ax.set_ylabel("mm_length (um)")
    ax.set_zlabel(f"Total transmission @ {wl_c:.3f} um")
    ax.set_title("MMI 2x2 -- total transmission (through + cross)")
    total_surface_path = os.path.join(sweep_dir, "sweep_surface_total.png")
    plt.savefig(total_surface_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Total-transmission surface saved to {total_surface_path}")

    # Splitting ratio: through / (through + cross), 0.5 = perfect 50:50 split.
    # Diverging colormap centered at 0.5 -- blue = cross-dominant, red = through-dominant.
    with np.errstate(invalid="ignore", divide="ignore"):
        Z_split = np.where(Z_total > 0, Z["through"] / Z_total, np.nan)

    plt.figure()
    mesh_plot = plt.pcolormesh(X, Y, Z_split, shading="auto", cmap="RdBu_r",
                               vmin=0.0, vmax=1.0)
    plt.colorbar(mesh_plot, label=f"Through / (Through + Cross) @ {wl_c:.3f} um")
    plt.xlabel("mm_width (um)")
    plt.ylabel("mm_length (um)")
    plt.title("MMI 2x2 -- splitting ratio (0.5 = ideal)")
    plt.tight_layout()
    split_heatmap_path = os.path.join(sweep_dir, "sweep_heatmap_splitting_ratio.png")
    plt.savefig(split_heatmap_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Splitting-ratio heatmap saved to {split_heatmap_path}")

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(X, Y, Z_split, cmap="RdBu_r", vmin=0.0, vmax=1.0)
    ax.set_xlabel("mm_width (um)")
    ax.set_ylabel("mm_length (um)")
    ax.set_zlabel(f"Splitting ratio @ {wl_c:.3f} um")
    ax.set_title("MMI 2x2 -- splitting ratio (0.5 = ideal)")
    split_surface_path = os.path.join(sweep_dir, "sweep_surface_splitting_ratio.png")
    plt.savefig(split_surface_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Splitting-ratio surface saved to {split_surface_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file (sim_type must be EME)", required=True)
    parser.add_argument("--save-plots", action="store_true",
                        help="save per-point config/data/performance.png in addition to the summary")
    parser.add_argument("--full-output", action="store_true",
                        help="save GDS layouts, cross-sections, and EME field plots at the "
                             "center (length, width) point. Implies --save-plots; much heavier I/O.")
    args = parser.parse_args()
    run_efficient_sweep(args)
