"""
Efficient parallel mm_length x waveguide_offset sweep for MMI 1x2 (EME only).

Cell stack: [input_wg] + [input_taper]*N + [MMI] + [output_taper]*N + [output_wg].
The leading N+2 cells (input_wg, input tapers, MMI) never depend on mm_length
(constant cross-section regardless of cell length) or waveguide_offset (only the
output side is shifted). So their modes are computed once per wavelength, globally,
in the main process, then handed to every worker.

The trailing N+1 cells (output tapers, output_wg) depend on waveguide_offset but not
mm_length, so each worker computes them once per (wavelength, offset) and reuses them
across all mm_length values it's assigned, only recomputing the cheap S-matrix per
length. No mode is solved more than once across the whole sweep.
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
from dataclasses import asdict, replace as dc_replace
from multiprocessing import Pool

import meep as mp
import meow as mw
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)
from meow import compute_modes_tidy3d as compute_modes

from engine.helpers.general import make_output_folders
from engine.helpers.meow import load_1x2_MMI, stack_extrusions_GF
from engine.hyperparams.MMI import HyperParamsMMI as HyperParams
from engine.plotting.meow import plot_eme_fields

N_PROCESSES = 8  # number of offset values to process concurrently, one core each

MM_LENGTH_SWEEP = np.linspace(3.5, 4.5, 21)          # 3.5:0.05:4.5 um
WAVEGUIDE_OFFSET_SWEEP = np.linspace(0.9, 1.1, 21)   # 0.9:0.01:1.1 um


def _build_cells(hp, mesh, extrusion_rules, taper_Ls, length, offset):
    mmi_gds = load_1x2_MMI(
        length, hp.sim.waveguide_length, hp.geo.waveguide_width,
        offset, hp.sim.oxide_width, 0.0,
        hp.geo.taper_length, hp.geo.mm_width, hp.geo.taper_width,
    )
    structs = mw.extrude_gds(mmi_gds, extrusion_rules)
    Ls = np.concatenate((
        [hp.sim.waveguide_length], taper_Ls, [length], taper_Ls, [hp.sim.waveguide_length]
    ))
    return mw.create_cells(structures=structs, mesh=mesh, Ls=Ls)


def _init_worker(hp, mesh, extrusion_rules, taper_Ls, num_modes_list, split,
                  shared_modes, shared_css, wavelengths, sweep_dir, save_plots, full_output):
    global _hp, _mesh, _extrusion_rules, _taper_Ls, _num_modes_list, _split
    global _shared_modes, _shared_css, _wavelengths, _sweep_dir, _save_plots, _full_output
    _hp = hp
    _mesh = mesh
    _extrusion_rules = extrusion_rules
    _taper_Ls = taper_Ls
    _num_modes_list = num_modes_list
    _split = split
    _shared_modes = shared_modes
    _shared_css = shared_css
    _wavelengths = wavelengths
    _sweep_dir = sweep_dir
    _save_plots = save_plots
    _full_output = full_output


def _run_offset(offset):
    print(f"[sweep] worker starting waveguide_offset = {offset:.3f} um")
    center_idx = int(np.floor(len(_wavelengths) / 2))
    # GDS/cross-section/field plots are for visual inspection only; generating them for
    # all 21 mm_length values per offset (441 total) repeatedly inside this same
    # long-lived worker process is what drives memory usage up without bound, since
    # there's no process recycling between mm_length values to hand memory back to the
    # OS (unlike mmi_1x2_sweep_parallel.py's maxtasksperchild=1). One representative
    # point per offset is enough to inspect the geometry/fields.
    center_length = float(MM_LENGTH_SWEEP[len(MM_LENGTH_SWEEP) // 2])

    ref_cells = _build_cells(_hp, _mesh, _extrusion_rules, _taper_Ls, float(MM_LENGTH_SWEEP[0]), offset)
    output_cells_ref = ref_cells[_split:]
    output_num_modes = _num_modes_list[_split:]

    output_modes = {}
    output_css = {}
    for wl_idx, wavelength in enumerate(_wavelengths):
        env = mw.Environment(wl=float(wavelength), T=25.0)
        css = [mw.CrossSection.from_cell(cell=c, env=env) for c in output_cells_ref]
        output_modes[wl_idx] = [
            compute_modes(cs, num_modes=n) for cs, n in zip(css, output_num_modes)
        ]
        output_css[wl_idx] = css

    result = {}
    for length in MM_LENGTH_SWEEP:
        length = float(length)
        cells_new = _build_cells(_hp, _mesh, _extrusion_rules, _taper_Ls, length, offset)

        powers, angles = [], []
        for wl_idx in range(len(_wavelengths)):
            full_modes = _shared_modes[wl_idx] + output_modes[wl_idx]
            S, _port_map = mw.compute_s_matrix(full_modes, cells=cells_new)
            c = S[0, 1]
            powers.append(float(np.abs(c) ** 2))
            angles.append(float(np.angle(c, deg=True)))
            print(
                f"  mm_length={length:.2f} waveguide_offset={offset:.3f} "
                f"wl={_wavelengths[wl_idx]:.4f}: P={powers[-1]:.4f}"
            )

        result[length] = {"powers": powers, "angles": angles}

        if _save_plots or _full_output:
            point_dir = os.path.join(_sweep_dir, f"len_{length:.2f}_off_{offset:.3f}")
            make_output_folders(point_dir)

            hp_dict = asdict(_hp)
            hp_dict["geo"]["mm_length"] = length
            hp_dict["geo"]["waveguide_offset"] = offset
            hp_dict["outputfolder"] = point_dir
            with open(os.path.join(point_dir, "config.json"), "w") as f:
                print(json.dumps(hp_dict, indent=4, default=lambda x: x.name), file=f)

            with open(os.path.join(point_dir, "data", "sim_results.pkl"), "wb") as f:
                pickle.dump({
                    "powers": powers,
                    "angles": angles,
                    "wavelengths": _wavelengths,
                }, f)

            plt.figure()
            plt.plot(_wavelengths, powers)
            plt.xlabel("Wavelength (μm)")
            plt.ylabel("Transmission")
            plt.title(f"MMI 1x2 — mm_length={length:.2f} μm, waveguide_offset={offset:.3f} μm")
            plt.tight_layout()
            plt.savefig(os.path.join(point_dir, "performance.png"), dpi=300)
            plt.close("all")

        if _full_output and length == center_length:
            mmi_gds = load_1x2_MMI(
                length, _hp.sim.waveguide_length, _hp.geo.waveguide_width,
                offset, _hp.sim.oxide_width, 0.0,
                _hp.geo.taper_length, _hp.geo.mm_width, _hp.geo.taper_width,
            )
            plt.figure()
            mmi_gds.plot()
            plt.savefig(os.path.join(point_dir, "layout", "eme_gds.png"), dpi=300)
            plt.close("all")

            css_full = _shared_css[center_idx] + output_css[center_idx]
            for cs_num, cs in enumerate(css_full):
                plt.clf()
                plt.figure()
                mw.visualize(cs, show=False)
                plt.savefig(
                    os.path.join(point_dir, "layout", "cross_sections", f"cs_{cs_num}.png"), dpi=300
                )
                plt.close("all")

            full_modes_center = _shared_modes[center_idx] + output_modes[center_idx]
            plot_eme_fields(dc_replace(_hp, outputfolder=point_dir), cells_new, full_modes_center)

            del mmi_gds, css_full, full_modes_center
            gc.collect()

        del cells_new
        gc.collect()

    return (float(offset), result)


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
        logger = logging.getLogger("mmi_1x2_sweep_parallel_efficient")
        logger.addHandler(handler)
        logger.setLevel("INFO")
    else:
        logger = logging.getLogger("mmi_1x2_sweep_parallel_efficient")

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
    split = hp.sim.num_cells_taper + 2
    mesh = mw.Mesh2D(
        x=np.linspace(-hp.sim.oxide_width / 2, hp.sim.oxide_width / 2,
                      int(hp.sim.res * hp.sim.oxide_width) + 1),
        y=np.linspace(-hp.sim.oxide_height / 2, hp.sim.oxide_height / 2,
                      int(hp.sim.res * hp.sim.oxide_height) + 1),
    )

    # Reference cells — input_wg/taper/MMI cross-sections are independent of both
    # mm_length and waveguide_offset, so any reference (length, offset) works here.
    ref_cells = _build_cells(hp, mesh, extrusion_rules, taper_Ls,
                              float(MM_LENGTH_SWEEP[0]), float(WAVEGUIDE_OFFSET_SWEEP[0]))
    assert len(ref_cells) == len(num_modes_list) == 2 * hp.sim.num_cells_taper + 3
    shared_cells_ref = ref_cells[:split]
    shared_num_modes = num_modes_list[:split]

    if mp.am_master(): logger.info("Computing shared (input/MMI) modes once per wavelength")
    shared_modes = {}
    shared_css = {}
    for wl_idx, wavelength in enumerate(wavelengths):
        env = mw.Environment(wl=float(wavelength), T=25.0)
        css = [mw.CrossSection.from_cell(cell=c, env=env) for c in shared_cells_ref]
        shared_modes[wl_idx] = [
            compute_modes(cs, num_modes=n) for cs, n in zip(css, shared_num_modes)
        ]
        shared_css[wl_idx] = css

    if mp.am_master():
        logger.info(f"Dispatching {len(WAVEGUIDE_OFFSET_SWEEP)} offsets across {N_PROCESSES} processes")

    with Pool(
        processes=N_PROCESSES,
        maxtasksperchild=1,
        initializer=_init_worker,
        initargs=(hp, mesh, extrusion_rules, taper_Ls, num_modes_list, split,
                  shared_modes, shared_css, wavelengths, sweep_dir,
                  args.save_plots, args.full_output),
    ) as pool:
        results = pool.map(_run_offset, [float(o) for o in WAVEGUIDE_OFFSET_SWEEP])

    if mp.am_master():
        _plot_summary(sweep_dir, wavelengths, results)


def _plot_summary(sweep_dir, wavelengths, results):
    center_idx = int(np.floor(len(wavelengths) / 2))
    Z = np.full((len(MM_LENGTH_SWEEP), len(WAVEGUIDE_OFFSET_SWEEP)), np.nan)
    all_results = {}

    for offset, length_dict in results:
        j = int(np.argmin(np.abs(WAVEGUIDE_OFFSET_SWEEP - offset)))
        all_results[offset] = length_dict
        for length, d in length_dict.items():
            i = int(np.argmin(np.abs(MM_LENGTH_SWEEP - length)))
            # 1x2 MMI is mirror-symmetric, so only one output arm is simulated;
            # total transmission is twice the simulated single-output-port power.
            Z[i, j] = 2 * d["powers"][center_idx]

    with open(os.path.join(sweep_dir, "sweep_results.pkl"), "wb") as f:
        pickle.dump({
            "mm_length": MM_LENGTH_SWEEP,
            "waveguide_offset": WAVEGUIDE_OFFSET_SWEEP,
            "total_transmission": Z,
            "wavelength_center": wavelengths[center_idx],
            "all_results": all_results,
        }, f)

    X, Y = np.meshgrid(WAVEGUIDE_OFFSET_SWEEP, MM_LENGTH_SWEEP)

    plt.figure()
    mesh_plot = plt.pcolormesh(X, Y, Z, shading="auto", cmap="viridis")
    plt.colorbar(mesh_plot, label=f"Total transmission @ λ={wavelengths[center_idx]:.3f} μm")
    plt.xlabel("waveguide_offset (μm)")
    plt.ylabel("mm_length (μm)")
    plt.title("MMI 1x2 — length / offset sweep (parallel, efficient)")
    plt.tight_layout()
    heatmap_path = os.path.join(sweep_dir, "sweep_heatmap.png")
    plt.savefig(heatmap_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Heatmap saved to {heatmap_path}")

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(X, Y, Z, cmap="viridis")
    ax.set_xlabel("waveguide_offset (μm)")
    ax.set_ylabel("mm_length (μm)")
    ax.set_zlabel(f"Total transmission @ λ={wavelengths[center_idx]:.3f} μm")
    ax.set_title("MMI 1x2 — length / offset sweep (parallel, efficient)")
    surface_path = os.path.join(sweep_dir, "sweep_surface.png")
    plt.savefig(surface_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Surface plot saved to {surface_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file (sim_type should be EME)", required=True)
    parser.add_argument("--save-plots", action="store_true",
                        help="save per-point config/data/performance.png in addition to the summary")
    parser.add_argument("--full-output", action="store_true",
                        help="match mmi_1x2_sweep_parallel.py's per-point output exactly: "
                             "config/data/performance.png plus layout/eme_gds.png, "
                             "layout/cross_sections/cs_*.png, and fields/*.png (at the center "
                             "wavelength). Implies --save-plots; much heavier I/O.")
    args = parser.parse_args()
    run_efficient_sweep(args)
