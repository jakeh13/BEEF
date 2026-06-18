import sys as _sys, os as _os
_FECES_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _FECES_ROOT)
_sys.path.insert(0, _os.path.join(_FECES_ROOT, "masters"))
# Each worker process runs a single-rank meep sim; pin it to one core so
# N_PROCESSES workers don't oversubscribe the machine via OpenMP threading.
_os.environ.setdefault("OMP_NUM_THREADS", "1")
del _sys, _os

import argparse
import json
import os
import pickle
import tempfile
from multiprocessing import Pool
from types import SimpleNamespace

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

from mmi_1x2_master import MMI_1x2

N_PROCESSES = 8  # number of sims to run concurrently, one processor each

MM_LENGTH_SWEEP = np.linspace(3.5, 4.5, 21)          # 3.5:0.05:4.5 um
WAVEGUIDE_OFFSET_SWEEP = np.linspace(0.9, 1.1, 21)   # 0.9:0.01:1.1 um


def _run_point(point):
    base_cfg, sweep_dir, length, offset = point

    cfg = json.loads(json.dumps(base_cfg))
    cfg["geo"]["mm_length"] = float(length)
    cfg["geo"]["waveguide_offset"] = float(offset)
    cfg["outputfolder"] = os.path.join(sweep_dir, f"len_{length:.2f}_off_{offset:.3f}")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(cfg, f)
        tmp_path = f.name

    try:
        print(f"[sweep] mm_length = {length:.2f} um, waveguide_offset = {offset:.3f} um")
        mmi = MMI_1x2(SimpleNamespace(config=tmp_path))
        mmi.run()
        return (float(length), float(offset), True)
    except Exception as e:
        print(f"[sweep] FAILED mm_length={length:.2f} waveguide_offset={offset:.3f}: {e}")
        return (float(length), float(offset), False)
    finally:
        os.unlink(tmp_path)


def run_sweep(args):
    with open(args.config) as f:
        base_cfg = json.load(f)

    sweep_dir = os.path.join(base_cfg["outputfolder"], "sweep")
    os.makedirs(sweep_dir, exist_ok=True)

    points = [
        (base_cfg, sweep_dir, length, offset)
        for length in MM_LENGTH_SWEEP
        for offset in WAVEGUIDE_OFFSET_SWEEP
    ]

    print(f"[sweep] Running {len(points)} points on {N_PROCESSES} processes (1 core each)")
    # maxtasksperchild=1 gives every sim its own fresh process, avoiding any
    # state (e.g. logger handlers) leaking across sims sharing a worker.
    with Pool(processes=N_PROCESSES, maxtasksperchild=1) as pool:
        results = pool.map(_run_point, points)

    completed = [(length, offset) for (length, offset, ok) in results if ok]
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
    plt.title("MMI 1x2 — length / offset sweep (parallel)")
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
    ax.set_title("MMI 1x2 — length / offset sweep (parallel)")
    surface_path = os.path.join(sweep_dir, "sweep_surface.png")
    plt.savefig(surface_path, dpi=300)
    plt.close("all")
    print(f"[sweep] Surface plot saved to {surface_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="base config file (sim_type should be EME)", required=True)
    args = parser.parse_args()
    run_sweep(args)
