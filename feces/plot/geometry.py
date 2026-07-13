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
    show_hybrid_sections: bool = False,
    ax=None,
    title: Optional[str] = None,
    output_path: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 4),
    wg_layer: Tuple[int, int] = (2, 0),
    ox_layer: Tuple[int, int] = (1, 0),
    **geometry_kwargs,
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
        show_hybrid_sections: shade the three hybrid solver sections (input FDTD /
            EME body / output FDTD) with boundary lines.  Requires the solver to
            be a HybridSolver and the geometry to expose taper_length + mm_length.
        ax: existing Axes; creates new figure if None.
        title: plot title (defaults to class name of geometry).
        output_path: full path for output PNG; directory created if needed.
        figsize: figure size (width, height) in inches.
        wg_layer: GDS layer tuple for waveguide polygons.
        ox_layer: GDS layer tuple for oxide cladding polygon.
        **geometry_kwargs: forwarded to geometry.to_gds_eme() (e.g. output_idx
            to draw only a single output branch, matching what the solver
            actually monitors).

    Returns:
        (fig, ax)
    """
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.patches import Polygon as MplPolygon, Rectangle

    # Geometry-building params from solver if available
    if solver is not None and hasattr(solver, "lengths"):
        # EME solver
        straight_length = solver.lengths[0]
        oxide_width = solver.oxide_width
    elif solver is not None and hasattr(solver, "solvers") and hasattr(solver, "boundaries"):
        # HybridSolver — delegate domain sizing to the first sub-solver in the chain
        _ref = solver.solvers[0]
        _wl = getattr(_ref, "waveguide_length", 0.5)
        _dpml = getattr(_ref, "dpml", 0.0)
        straight_length = _wl + _dpml
        oxide_width = getattr(_ref, "sim_width", 5.0) + 2 * _dpml
    elif solver is not None and hasattr(solver, "sim_width"):
        # FDTD solver — extend stubs through PML so waveguides reach domain edge
        _wl = getattr(solver, "waveguide_length", 0.5)
        _dpml = getattr(solver, "dpml", 0.0)
        straight_length = _wl + _dpml
        oxide_width = solver.sim_width + 2 * _dpml   # full domain y-span (inner + 2 PML)
    else:
        straight_length = 0.5
        oxide_width = 5.0

    gds_comp = geometry.to_gds_eme(
        straight_length=straight_length,
        oxide_width=oxide_width,
        **geometry_kwargs,
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
    _fdtd_attrs = ("source_margin", "monitor_margin", "waveguide_length", "dpml")
    _fdtd_ref = (solver.solvers[0]
                 if solver is not None and hasattr(solver, "solvers") and hasattr(solver, "boundaries")
                 else solver)
    _has_fdtd = _fdtd_ref is not None and all(hasattr(_fdtd_ref, a) for a in _fdtd_attrs)

    is_fdtd = show_monitors and _has_fdtd
    if is_fdtd:
        tl = geometry.total_length
        _dpml = _fdtd_ref.dpml
        total_gds = tl + 2 * straight_length  # = tl + 2*(wl+dpml), full domain width

        monitors = [
            (_dpml + _fdtd_ref.source_margin,               "source",  "crimson",    "-"),
            (_dpml + _fdtd_ref.monitor_margin,              "in mon",  "darkorange", "--"),
            (total_gds - _dpml - _fdtd_ref.monitor_margin,  "out mon", "seagreen",   "--"),
        ]
        y0, y1 = -oxide_width / 2, oxide_width / 2          # domain edges
        y_in0, y_in1 = y0 + _dpml, y1 - _dpml               # inner sim bounds (no PML)
        for x_gds, label, color, ls in monitors:
            (line,) = ax.plot([x_gds, x_gds], [y_in0, y_in1], color=color, lw=1.5, linestyle=ls, alpha=0.9, zorder=5)
            line.set_path_effects([pe.withStroke(linewidth=3.5, foreground="white")])
            ax.text(
                x_gds, oxide_width / 2 * 0.62, label,
                ha="center", va="top", fontsize=7, color=color, zorder=6,
                rotation=90,
                bbox=dict(boxstyle="round,pad=0.1", facecolor="white", alpha=0.8, edgecolor=color, lw=0.5),
            )

    # PML absorbing boundary regions
    if show_pml and _has_fdtd:
        tl = geometry.total_length
        dpml = _fdtd_ref.dpml
        total_gds = tl + 2 * straight_length  # full domain width (waveguides extend to edges)
        y0, y1 = -oxide_width / 2, oxide_width / 2           # domain edges
        y_in0, y_in1 = y0 + dpml, y1 - dpml                  # inner sim bounds
        pml_kw = dict(alpha=0.18, color="#8B0000", zorder=2)
        # Left / right PML strips — full domain height, waveguides run through them
        ax.add_patch(Rectangle((0.0,              y0), dpml, oxide_width, **pml_kw))
        ax.add_patch(Rectangle((total_gds - dpml, y0), dpml, oxide_width, **pml_kw))
        # Top / bottom PML strips — full domain x-span
        ax.add_patch(Rectangle((0.0, y0),    total_gds, dpml, **pml_kw))
        ax.add_patch(Rectangle((0.0, y_in1), total_gds, dpml, **pml_kw))
        for x_pml, label_pml in [(dpml / 2, "PML"), (total_gds - dpml / 2, "PML")]:
            ax.text(
                x_pml, y_in1 * 0.88, label_pml,
                ha="center", va="top", fontsize=7, color="#8B0000", zorder=7,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.65, edgecolor="none"),
            )
        # Edge lines at PML boundaries — clipped to domain
        edge_kw = dict(color="#5C0000", lw=0.8, alpha=0.75, zorder=8)
        for x_edge in (0.0, dpml, total_gds - dpml, total_gds):
            ax.plot([x_edge, x_edge], [y0, y1], **edge_kw)
        for y_edge in (y0, y_in0, y_in1, y1):
            ax.plot([0.0, total_gds], [y_edge, y_edge], **edge_kw)

    # Hybrid section overlays — generic over any number/type of solvers
    _is_hybrid = (
        show_hybrid_sections
        and solver is not None
        and hasattr(solver, "solvers")
        and hasattr(solver, "boundaries")
        and hasattr(geometry, "total_length")
    )
    if _is_hybrid:
        # device frame (x=0 at geometry centre) -> GDS plot frame
        x_shift = geometry.total_length / 2 + straight_length
        bounds_gds = [b + x_shift for b in solver.boundaries]
        y0_h, y1_h = -(oxide_width / 2), +(oxide_width / 2)
        h_h = oxide_width
        palette = plt.cm.Set2.colors

        for k, sub_solver in enumerate(solver.solvers):
            xa, xb = bounds_gds[k], bounds_gds[k + 1]
            color = palette[k % len(palette)]
            ax.add_patch(Rectangle((xa, y0_h), xb - xa, h_h, alpha=0.15, color=color, zorder=3))
            ax.text((xa + xb) / 2, y1_h * 0.55, type(sub_solver).__name__,
                    ha="center", va="top", fontsize=7, color=color, zorder=10)

        for k in range(1, len(solver.solvers)):
            x_b = bounds_gds[k]
            label_b = f"{type(solver.solvers[k-1]).__name__}|{type(solver.solvers[k]).__name__}"
            ax.plot([x_b, x_b], [y0_h, y1_h], color="#333333", lw=1.2, linestyle="--",
                    alpha=0.8, zorder=9)
            ax.text(x_b, y1_h * 0.80, label_b, ha="center", va="top",
                    fontsize=7, color="#333333", zorder=10,
                    bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                              alpha=0.75, edgecolor="none"))

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
