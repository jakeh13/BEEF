"""Top-down device geometry plot."""
from __future__ import annotations

import os
from typing import List, Optional, Tuple

import numpy as np


def plot_geometry(
    geometry,
    solver=None,
    show_eme_regions: bool = False,
    region_labels: Optional[List[str]] = None,
    show_monitors: bool = False,
    show_pml: bool = False,
    ax=None,
    title: Optional[str] = None,
    output_path: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 4),
    wg_layer: Tuple[int, int] = (2, 0),
    ox_layer: Tuple[int, int] = (1, 0),
):
    """Top-down device geometry plot.

    Works with any DeviceGeometry (EME, FDTD, or hybrid).

    Args:
        geometry: DeviceGeometry instance.
        solver: EMESolver or FDTDSolver (drives show_eme_regions / show_monitors).
        show_eme_regions: overlay EME region shading from solver.lengths.
        region_labels: custom per-region labels (len must equal len(solver.lengths)).
        show_monitors: draw source and DFT monitor lines (FDTDSolver only).
        show_pml: shade PML absorbing boundary regions (FDTDSolver only).
        ax: existing Axes; creates new figure if None.
        title: plot title (defaults to class name of geometry).
        output_path: full path for output PNG; directory created if needed.
        figsize: figure size (width, height) in inches.
        wg_layer: GDS layer tuple for waveguide polygons.
        ox_layer: GDS layer tuple for oxide cladding polygon.

    Returns:
        (fig, ax)
    """
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.patches import Polygon as MplPolygon

    # Geometry-building params from solver if available
    if solver is not None and hasattr(solver, "lengths"):
        straight_length = solver.lengths[0]
        oxide_width = solver.oxide_width
    else:
        straight_length = 0.5
        oxide_width = 5.0

    gds_comp = geometry.to_gds_eme(
        straight_length=straight_length,
        oxide_width=oxide_width,
    )

    # gdsfactory returns KLayout polygon objects; coords are in dbu (1000 dbu = 1 µm)
    poly_dict = gds_comp.get_polygons(by="tuple")
    dbu = gds_comp.kcl.dbu  # µm per dbu (typically 0.001)

    def _to_xy(kpoly):
        return np.array([[pt.x * dbu, pt.y * dbu] for pt in kpoly.each_point_hull()])

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Oxide background
    for kp in poly_dict.get(ox_layer, []):
        patch = MplPolygon(
            _to_xy(kp), closed=True,
            facecolor="#d8d8d8", edgecolor="none", alpha=0.5, zorder=0,
        )
        ax.add_patch(patch)

    # Waveguide core
    for kp in poly_dict.get(wg_layer, []):
        patch = MplPolygon(
            _to_xy(kp), closed=True,
            facecolor="royalblue", edgecolor="steelblue", lw=0.4, alpha=0.9, zorder=1,
        )
        ax.add_patch(patch)

    # EME region overlays
    if show_eme_regions and solver is not None and hasattr(solver, "lengths"):
        lens = solver.lengths
        boundaries = np.concatenate([[0.0], np.cumsum(lens)])
        palette = plt.cm.Set3.colors  # 12 pastel colors

        for i, (x0, x1) in enumerate(zip(boundaries[:-1], boundaries[1:])):
            color = palette[i % len(palette)]
            ax.axvspan(x0, x1, alpha=0.22, color=color, zorder=2)
            ax.axvline(x0, color="#444444", lw=0.8, linestyle="--", alpha=0.65, zorder=3)
            label = (
                region_labels[i]
                if region_labels and i < len(region_labels)
                else f"R{i + 1}"
            )
            ax.text(
                (x0 + x1) / 2, oxide_width / 2 * 0.88,
                label, ha="center", va="top", fontsize=7.5, color="#222222", zorder=4,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.65, edgecolor="none"),
            )
        ax.axvline(boundaries[-1], color="#444444", lw=0.8, linestyle="--", alpha=0.65, zorder=3)

    # FDTD source and monitor overlays
    # GDS x-coords: source_position from left, SM at sm_monitor_position,
    # MM at (total_gds_length - mm_monitor_position) from left.
    _fdtd_attrs = ("source_position", "sm_monitor_position", "mm_monitor_position",
                   "waveguide_length", "dpml")
    _has_fdtd = solver is not None and all(hasattr(solver, a) for a in _fdtd_attrs)

    is_fdtd = show_monitors and _has_fdtd
    if is_fdtd:
        sl = solver.waveguide_length
        tl = geometry.total_length
        total_gds = tl + 2 * sl

        monitors = [
            (solver.source_position,              "source",  "crimson",   "-"),
            (solver.sm_monitor_position,           "SM mon",  "darkorange","--"),
            (total_gds - solver.mm_monitor_position, "MM mon", "seagreen", "--"),
        ]
        for x_gds, label, color, ls in monitors:
            line = ax.axvline(x_gds, color=color, lw=1.5, linestyle=ls, alpha=0.9, zorder=5)
            line.set_path_effects([pe.withStroke(linewidth=3.5, foreground="white")])
            ax.text(
                x_gds, oxide_width / 2 * 0.62, label,
                ha="center", va="top", fontsize=7, color=color, zorder=6,
                rotation=90,
                bbox=dict(boxstyle="round,pad=0.1", facecolor="white", alpha=0.8, edgecolor=color, lw=0.5),
            )

    # PML absorbing boundary regions
    if show_pml and _has_fdtd:
        sl = solver.waveguide_length
        tl = geometry.total_length
        total_gds = tl + 2 * sl
        dpml = solver.dpml
        pml_kw = dict(alpha=0.18, color="#8B0000", zorder=2)
        # Left / right PML (vertical bands outside device x-range)
        ax.axvspan(-dpml, 0.0,                    **pml_kw)
        ax.axvspan(total_gds, total_gds + dpml,   **pml_kw)
        for x_pml, label_pml in [(-dpml / 2, "PML"), (total_gds + dpml / 2, "PML")]:
            ax.text(
                x_pml, oxide_width / 2 * 0.88, label_pml,
                ha="center", va="top", fontsize=7, color="#8B0000", zorder=7,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.65, edgecolor="none"),
            )
        # Top / bottom PML (horizontal bands at domain y-edges)
        ax.axhspan(-oxide_width / 2, -oxide_width / 2 + dpml, **pml_kw)
        ax.axhspan( oxide_width / 2 - dpml,  oxide_width / 2, **pml_kw)
        # Edge lines at every PML boundary
        edge_kw = dict(color="#5C0000", lw=0.8, alpha=0.75, zorder=8)
        for x_edge in (-dpml, 0.0, total_gds, total_gds + dpml):
            ax.axvline(x_edge, **edge_kw)
        for y_edge in (-oxide_width / 2, -oxide_width / 2 + dpml,
                        oxide_width / 2 - dpml,  oxide_width / 2):
            ax.axhline(y_edge, **edge_kw)

    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_title(title or type(geometry).__name__)
    ax.grid(False)

    if output_path:
        parent = os.path.dirname(output_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")

    return fig, ax
