# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Georgia Tech Research Corporation
"""EME field visualization.

Migrated from engine/plotting/meow.py with updated signature (output_dir replaces hp).
"""
from __future__ import annotations

import os

import numpy as np
from matplotlib import pyplot as plt


def plot_eme_fields(
    cells,
    modes,
    output_dir: str,
    tag: str = "1",
    input_mode: int = 1,
) -> None:
    """Propagate modes through cells and save 6-component field plots.

    Args:
        cells: meow cell list (from mw.create_cells).
        modes: list of mode lists per cell cross-section.
        output_dir: directory where plots are saved (fields/ subfolder created).
        tag: filename tag to distinguish multiple calls.
        input_mode: 1-based mode index to excite at the left boundary.
    """
    import meow as mw

    fields_dir = os.path.join(output_dir, "fields")
    os.makedirs(fields_dir, exist_ok=True)

    z = np.linspace(cells[0].z_min, cells[-1].z_max, 500)
    exc = input_mode - 1  # 0-based

    for component, label in [
        ("Ex", "|Ex|"), ("Ey", "|Ey|"), ("Ez", "|Ez|"),
        ("Hx", "|Hx|"), ("Hy", "|Hy|"), ("Hz", "|Hz|"),
    ]:
        field, x = mw.propagate_modes(
            modes, cells, y=0.0, z=z, excite_mode_l=exc, component=component
        )
        fig, ax = plt.subplots(figsize=(10, 4))
        im = ax.pcolormesh(z, x, np.abs(field).T, shading="auto", cmap="jet")
        ax.set_xlabel("z [µm]")
        ax.set_ylabel("x [µm]")
        ax.set_title(f"{label}(x, z) — input mode {input_mode}")
        plt.colorbar(im, ax=ax, label=label)
        plt.tight_layout()
        plt.savefig(os.path.join(fields_dir, f"eme_{component}_{tag}.png"), dpi=300)
        plt.close("all")
