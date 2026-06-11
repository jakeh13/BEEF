import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
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

from mmi_1x4_master import MMI_1x4

MM_LENGTH_SWEEP = np.arange(5.0, 15.5, 0.5)  # 21 points, 5.0–15.0 μm


def run_sweep(args):
    with open(args.config) as f:
        base_cfg = json.load(f)

    base_out = base_cfg["outputfolder"]
    sweep_dir = os.path.join(base_out, "sweep")
    os.makedirs(sweep_dir, exist_ok=True)

    completed_vals = []

    for val in MM_LENGTH_SWEEP:
        if mp.am_master():
            print(f"[sweep] mm_length = {val:.2f} μm")

        cfg = json.loads(json.dumps(base_cfg))
        cfg["geo"]["mm_length"] = float(val)
        cfg["outputfolder"] = os.path.join(sweep_dir, f"mml_{val:.2f}")

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(cfg, f)
            tmp_path = f.name

        try:
            mmi = MMI_1x4(SimpleNamespace(config=tmp_path))
            mmi.run()
            completed_vals.append(val)
        finally:
            os.unlink(tmp_path)

    if mp.am_master():
        _plot_summary(sweep_dir, base_cfg, completed_vals)


def _plot_summary(sweep_dir, base_cfg, completed_vals):
    sim_cfg = base_cfg.get("sim", {})
    wl_min = sim_cfg.get("wavelength_min", 1.26)
    wl_max = sim_cfg.get("wavelength_max", 1.36)
    wl_pts = sim_cfg.get("wavelength_pts", 10)
    wavelengths = np.linspace(wl_min, wl_max, wl_pts)
    center_idx = int(np.floor(len(wavelengths) / 2))

    powers_by_port = [[] for _ in range(4)]

    for val in completed_vals:
        pkl_path = os.path.join(sweep_dir, f"mml_{val:.2f}", "data", "sim_results.pkl")
        if not os.path.exists(pkl_path):
            print(f"[sweep] WARNING: missing results for mm_length={val:.2f}, skipping")
            continue
        with open(pkl_path, "rb") as f:
            data = pickle.load(f)
        for port_idx in range(4):
            powers_by_port[port_idx].append(data["powers_ports"][port_idx][center_idx])

    plt.figure()
    for port_idx, powers in enumerate(powers_by_port):
        plt.plot(completed_vals[: len(powers)], powers, marker="o", label=f"Port {port_idx}")
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
    parser.add_argument("--config", "-c", help="base config file", required=True)
    args = parser.parse_args()
    run_sweep(args)
