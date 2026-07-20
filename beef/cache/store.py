"""On-disk lookup/save for the persistent solve_section() cache (see
beef/cache/__init__.py). Each matching config gets its own subfolder, named
by a hash of the config for O(1) lookup, holding config.json (human-readable,
also re-checked by full equality against a hash collision or schema drift),
result.pkl, and geometry.png.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import pickle
from typing import Any, Dict, Optional

import numpy as np

logger = logging.getLogger(__name__)


def _hash_config(config: Dict[str, Any]) -> str:
    blob = json.dumps(config, sort_keys=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:16]


def _folder(save_dir: str, config: Dict[str, Any]) -> str:
    return os.path.join(save_dir, _hash_config(config))


def lookup(save_dir: str, config: Dict[str, Any]) -> Optional[str]:
    """Return the matching run's folder path, or None if no cached run
    matches config (hash miss, or a hash hit whose config.json doesn't
    actually equal config — a collision/schema-drift guard, not just a
    hash-trust lookup)."""
    folder = _folder(save_dir, config)
    config_path = os.path.join(folder, "config.json")
    if not os.path.isfile(config_path):
        return None
    with open(config_path) as f:
        cached_config = json.load(f)
    if cached_config != config:
        logger.warning(
            "cache: hash match at %s but config.json differs — treating as "
            "a miss (hash collision or schema drift).", folder,
        )
        return None
    return folder


def load_result(folder: str) -> np.ndarray:
    with open(os.path.join(folder, "result.pkl"), "rb") as f:
        return pickle.load(f)


def _save_geometry_plot(config: Dict[str, Any], path: str) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon

    fig, ax = plt.subplots(figsize=(8, 3))
    for layer_key, polys in config["geometry_polygons"].items():
        is_wg = layer_key == "2,0"
        for poly in polys:
            ax.add_patch(MplPolygon(
                np.asarray(poly), closed=True,
                facecolor="#4a4a4a" if is_wg else "#d8d8d8",
                edgecolor="none", alpha=1.0 if is_wg else 0.5,
                zorder=1 if is_wg else 0,
            ))
    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (µm, section-local)")
    ax.set_ylabel("y (µm)")
    ax.set_title(f"{config['solver_type']}.{config['method']} — {config['geometry_class']} (cached)")
    ax.grid(False)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_run(save_dir: str, config: Dict[str, Any], result: np.ndarray) -> str:
    """Write config.json + result.pkl + geometry.png to a new/overwritten
    subfolder of save_dir, named by hash(config). Returns the folder path."""
    folder = _folder(save_dir, config)
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "config.json"), "w") as f:
        json.dump(config, f, indent=2, sort_keys=True)
    with open(os.path.join(folder, "result.pkl"), "wb") as f:
        pickle.dump(result, f)
    try:
        _save_geometry_plot(config, os.path.join(folder, "geometry.png"))
    except Exception:
        logger.warning("cache: failed to save geometry plot at %s", folder, exc_info=True)
    return folder
