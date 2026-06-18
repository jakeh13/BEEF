import sys as _sys, os as _os
_FECES_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _FECES_ROOT)
_sys.path.insert(0, _os.path.join(_FECES_ROOT, "masters"))
del _sys, _os

import argparse
import json
import os
import pickle
import tempfile
from types import SimpleNamespace

import meep as mp
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

from mmi_1x2_master import MMI_1x2

MM_LENGTH_SWEEP = np.linspace(3.5, 4.5, 21)          # 3.5:0.05:4.5 um
WAVEGUIDE_OFFSET_SWEEP = np.linspace(0.9, 1.1, 21)   # 0.9:0.01:1.1 um


def run_sweep(args):
    with open(args.config) as f:
        base_cfg = json.load(f)

    base_out = base_cfg["outputfolder"]
    sweep_dir = os.path.join(base_out, "sweep")
    os.makedirs(sweep_dir, exist_ok=True)

    completed = []  # list of (mm_length, waveguide_offset)

    for length in MM_LENGTH_SWEEP:
        for offset in WAVEGUIDE_OFFSET_SWEEP:
            if mp.am_master():
                print(f"[sweep] mm_length = {length:.2f} um, waveguide_offset = {offset:.3f} um")

            cfg = json.loads(json.dumps(base_cfg))
            cfg["geo"]["mm_length"] = float(length)
            cfg["geo"]["waveguide_offset"] = float(offset)
            cfg["outputfolder"] = os.path.join(sweep_dir, f"len_{length:.2f}_off_{offset:.3f}")

            with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
                json.dump(cfg, f)
                tmp_path = f.name

            try:
                mmi = MMI_1x2(SimpleNamespace(config=tmp_path))
                mmi.run()
                completed.append((float(length), float(offset)))
            finally:
                os.unlink(tmp_path)

    if mp.am_master():
        _plot_summary(sweep_dir, base_cfg, completed)


def _plot_summary(sweep_dir, base_cfg, completed):
    sim_cfg = base_cfg.get("sim", {})
    wl_min = sim_cfg.get("wavelength_min", 1.26)
    wl_max = sim_cfg.get("wavelength_max", 1.36)
    wl_pts = sim_cfg.get("wavelength_pts", 10)
    wavelengths = np.linspace(wl_min, wl_max, wl_pts)
    center_idx = int(np.floor(len(wavelengths) / 2))

    Z = np.full((len(MM_LENGTH_SWEEP), len(WAVEGUIDE_OFFSET_SWEEP)), np.nan)

    for length, offset in completed:
        i = int(np.argmin(np.abs(MM_LENGTH_SWEEP - length)))
        j = int(np.argmin(np.abs(WAVEGUIDE_OFFSET_SWEEP - offset)))
        pkl_path = os.path.join(sweep_dir, f"len_{length:.2f}_off_{offset:.3f}", "data", "sim_results.pkl")
        if not os.path.exists(pkl_path):
            print(f"[sweep] WARNING: missing results for length={length:.2f}, offset={offset:.3f}, skipping")
            continue
        with open(pkl_path, "rb") as f:
            data = pickle.load(f)
        # 1x2 MMI is mirror-symmetric, so only one output arm is simulated;
        # total transmission is twice the simulated single-output-port power.
        Z[i, j] = 2 * data["powers"][center_idx]

    with open(os.path.join(sweep_dir, "sweep_results.pkl"), "wb") as f:
        pickle.dump({
            "mm_length": MM_LENGTH_SWEEP,
            "waveguide_offset": WAVEGUIDE_OFFSET_SWEEP,
            "total_transmission": Z,
            "wavelength_center": wavelengths[center_idx],
        }, f)

    X, Y = np.meshgrid(WAVEGUIDE_OFFSET_SWEEP, MM_LENGTH_SWEEP)

    plt.figure()
    mesh = plt.pcolormesh(X, Y, Z, shading="auto", cmap="viridis")
    plt.colorbar(mesh, label=f"Total transmission @ λ={wavelengths[center_idx]:.3f} μm")
    plt.xlabel("waveguide_offset (μm)")
    plt.ylabel("mm_length (μm)")
    plt.title("MMI 1x2 — length / offset sweep")
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
    ax.set_title("MMI 1x2 — length / offset sweep")
    surface_path = os.path.join(sweep_dir, "sweep_surface.png")
    plt.savefig(surface_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Surface plot saved to {surface_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="base config file (sim_type should be EME)", required=True)
    args = parser.parse_args()
    run_sweep(args)
