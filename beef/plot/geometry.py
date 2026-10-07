# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Georgia Tech Research Corporation
"""Top-down device geometry plots — one function per solver type.

plot_eme_geometry() / plot_fdtd_geometry() / plot_hybrid_geometry() each draw
the device's waveguide/oxide polygons plus an overlay of whatever that solver
actually measures: EME's per-region mode-solve cells, FDTD's source/monitor/
PML layout, or (for a HybridSolver) the chain of FDTD/EME sections plus a
zoomed-in sub-plot per section — reusing the same drawing primitives
(_draw_polygons/_draw_eme_regions/_draw_marker/_draw_pml) so a hybrid
section's sub-plot reads the same way the matching standalone plot would.
"""
from __future__ import annotations

import os
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from ..solvers.fdtd import section_fdtd_layout

__all__ = ["plot_eme_geometry", "plot_fdtd_geometry", "plot_section_fdtd_geometry", "plot_hybrid_geometry"]


# ---------------------------------------------------------------------------
# Shared low-level primitives
# ---------------------------------------------------------------------------

def _new_ax(ax, figsize):
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    return fig, ax


def _draw_polygons(ax, gds_polys, wg_layer, ox_layer, x_shift: float = 0.0):
    """Oxide + waveguide-core patches from a {layer: [polygons]} dict (the
    format DeviceGeometry.to_gds_eme()/section_gds() return).

    x_shift: added to every x-coordinate. Used when gds_polys was built with
    a *different* (larger) stub than the local frame the rest of a plot
    (markers, PML rects) is drawn in — e.g. a section rendered with its
    waveguide extended dpml further into the PML (see
    _plot_hybrid_sub_regions) shifts gds_polys' own x_local=0 left by dpml
    relative to that shared frame, so x_shift=dpml re-aligns it.
    """
    from matplotlib.patches import Polygon as MplPolygon

    def _to_xy(poly):
        return np.array([[x + x_shift, y] for x, y in poly])

    for poly in gds_polys.get(ox_layer, []):
        ax.add_patch(MplPolygon(
            _to_xy(poly), closed=True,
            facecolor="#d8d8d8", edgecolor="none", alpha=0.5, zorder=0,
        ))
    for poly in gds_polys.get(wg_layer, []):
        ax.add_patch(MplPolygon(
            _to_xy(poly), closed=True,
            facecolor="royalblue", edgecolor="steelblue", lw=0.4, alpha=0.9, zorder=1,
        ))


def _draw_eme_regions(ax, boundaries, oxide_width, labels=None, palette=None):
    """Shade/label each EME region [boundaries[i], boundaries[i+1])."""
    import matplotlib.pyplot as plt

    palette = palette if palette is not None else plt.cm.Set3.colors  # 12 pastel colors
    for i, (x0, x1) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        color = palette[i % len(palette)]
        ax.axvspan(x0, x1, alpha=0.22, color=color, zorder=2)
        ax.axvline(x0, color="#444444", lw=0.8, linestyle="--", alpha=0.65, zorder=3)
        label = labels[i] if labels and i < len(labels) else f"R{i + 1}"
        ax.text(
            (x0 + x1) / 2, oxide_width / 2 * 0.88,
            label, ha="center", va="top", fontsize=7.5, color="#222222", zorder=4,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.65, edgecolor="none"),
        )
    ax.axvline(boundaries[-1], color="#444444", lw=0.8, linestyle="--", alpha=0.65, zorder=3)


def _draw_marker(ax, x_gds, y_lo, y_hi, label, color, ls):
    """A single vertical source/monitor line spanning [y_lo, y_hi]."""
    import matplotlib.patheffects as pe

    (line,) = ax.plot([x_gds, x_gds], [y_lo, y_hi], color=color, lw=1.5, linestyle=ls, alpha=0.9, zorder=5)
    line.set_path_effects([pe.withStroke(linewidth=3.5, foreground="white")])
    ax.text(
        x_gds, y_hi - 0.08 * (y_hi - y_lo), label,
        ha="center", va="top", fontsize=7, color=color, zorder=6,
        rotation=90,
        bbox=dict(boxstyle="round,pad=0.1", facecolor="white", alpha=0.8, edgecolor=color, lw=0.5),
    )


def _draw_fdtd_monitors(ax, dpml, total_gds, oxide_width, source_margin, monitor_margin,
                         output_monitor_centers=None, output_monitor_widths=None,
                         sim_width=None, src_yo=0.0, src_width=None):
    """source/in-mon/out-mon markers matching FDTDSolver.solve()'s layout:
    single source + co-located input-normalization monitor (both centered on
    src_yo — the actual input waveguide position, see plot_fdtd_geometry's
    docstring — spanning src_width, or the full sim_width when src_width is
    None: PML edge to PML edge, exclusive of PML, exactly mirroring solve()'s
    own default; a non-None src_width is the geometry-provided launch window
    solve() shrinks to, see geometry.source_aperture_at), plus one or more
    output monitors on the right (one per output_monitor_centers/widths
    entry, each spanning only its own y-range, or a single full-height line
    if unset)."""
    y0, y1 = -oxide_width / 2, oxide_width / 2
    y_in0, y_in1 = y0 + dpml, y1 - dpml

    src_x = dpml + source_margin
    in_x = dpml + monitor_margin
    out_x = total_gds - dpml - monitor_margin

    w_src = sim_width if src_width is None else src_width
    y_lo = max(src_yo - w_src / 2, y_in0)
    y_hi = min(src_yo + w_src / 2, y_in1)
    _draw_marker(ax, src_x, y_lo, y_hi, "source", "crimson", "-")
    _draw_marker(ax, in_x, y_lo, y_hi, "in mon", "darkorange", "--")

    if output_monitor_centers:
        out_widths = (
            output_monitor_widths
            if output_monitor_widths is not None
            else [sim_width] * len(output_monitor_centers)
        )
        for i, (yc, yw) in enumerate(zip(output_monitor_centers, out_widths)):
            label = "out mon" if len(output_monitor_centers) == 1 else f"out mon {i}"
            _draw_marker(ax, out_x, yc - yw / 2, yc + yw / 2, label, "seagreen", "--")
    else:
        _draw_marker(ax, out_x, y_in0, y_in1, "out mon", "seagreen", "--")


def _draw_pml(ax, dpml, total_gds, oxide_width):
    """PML absorbing-boundary strips on all four sides of an FDTD domain."""
    from matplotlib.patches import Rectangle

    y0, y1 = -oxide_width / 2, oxide_width / 2
    y_in1 = y1 - dpml
    pml_kw = dict(alpha=0.18, color="#8B0000", zorder=2)
    ax.add_patch(Rectangle((0.0, y0), dpml, oxide_width, **pml_kw))
    ax.add_patch(Rectangle((total_gds - dpml, y0), dpml, oxide_width, **pml_kw))
    ax.add_patch(Rectangle((0.0, y0), total_gds, dpml, **pml_kw))
    ax.add_patch(Rectangle((0.0, y_in1), total_gds, dpml, **pml_kw))
    for x_pml, label_pml in [(dpml / 2, "PML"), (total_gds - dpml / 2, "PML")]:
        ax.text(
            x_pml, y_in1 * 0.88, label_pml,
            ha="center", va="top", fontsize=7, color="#8B0000", zorder=7,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.65, edgecolor="none"),
        )


def _finalize(fig, ax, title, geometry, output_path):
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
        fig.savefig(output_path, dpi=600, bbox_inches="tight")


# ---------------------------------------------------------------------------
# EME
# ---------------------------------------------------------------------------

def plot_eme_geometry(
    geometry,
    solver,
    region_labels: Optional[List[str]] = None,
    ax=None,
    title: Optional[str] = None,
    output_path: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 4),
    wg_layer: Tuple[int, int] = (2, 0),
    ox_layer: Tuple[int, int] = (1, 0),
    **geometry_kwargs,
):
    """Top-down geometry plot for an EMESolver.

    Draws the device polygons plus one shaded/labeled region per entry in
    solver.lengths — exactly the cells EMESolver.solve()/solve_sweep() mode-
    solves (e.g. [wg_in, taper_in, MMI body, taper_out, wg_out]).

    Args:
        geometry: DeviceGeometry instance.
        solver: EMESolver — its .lengths/.oxide_width drive the domain size
            and region shading.
        region_labels: custom per-region labels (defaults to "R1", "R2", ...).
        ax: existing Axes; creates a new figure if None.
        title: plot title (defaults to the geometry class name).
        output_path: full path for output PNG; directory created if needed.
        figsize: figure size (width, height) in inches.
        wg_layer/ox_layer: GDS layer tuples for waveguide/oxide polygons.
        **geometry_kwargs: forwarded to geometry.to_gds_eme() (e.g. output_idx
            to draw only a single output branch).

    Returns:
        (fig, ax)
    """
    straight_length = solver.lengths[0]
    oxide_width = solver.oxide_width

    gds_polys = geometry.to_gds_eme(
        straight_length=straight_length, oxide_width=oxide_width, **geometry_kwargs
    )

    fig, ax = _new_ax(ax, figsize)
    _draw_polygons(ax, gds_polys, wg_layer, ox_layer)

    boundaries = np.concatenate([[0.0], np.cumsum(solver.lengths)])
    _draw_eme_regions(ax, boundaries, oxide_width, labels=region_labels)

    _finalize(fig, ax, title, geometry, output_path)
    return fig, ax


# ---------------------------------------------------------------------------
# FDTD
# ---------------------------------------------------------------------------

def plot_fdtd_geometry(
    geometry,
    solver,
    port: str = "top",
    show_pml: bool = True,
    ax=None,
    title: Optional[str] = None,
    output_path: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 4),
    wg_layer: Tuple[int, int] = (2, 0),
    ox_layer: Tuple[int, int] = (1, 0),
    **geometry_kwargs,
):
    """Top-down geometry plot for an FDTDSolver.

    Draws the device polygons plus the source/monitor markers and PML
    shading matching what FDTDSolver.solve() actually sets up: one source +
    co-located input-normalization monitor near the left PML — centered on
    geometry.input_ys' last (port="top") or first (port="bot") entry (the
    same entry either way for a single-input device), spanning the full
    sim_width (PML edge to PML edge, exclusive of PML) unless the geometry
    itself reports a launch window at the source cross-section
    (source_aperture_at() non-None, e.g. DirectionalCoupler with its
    physically-present unexcited second arm), in which case they span that
    window clamped to the non-PML domain, exactly like solve() — and one
    output monitor per entry in solver.output_monitor_centers/widths near
    the right PML (each spanning only its own y-range).

    Args:
        geometry: DeviceGeometry instance.
        solver: FDTDSolver — its .sim_width/.dpml/.source_margin/
            .monitor_margin/.output_monitor_centers drive domain size and
            marker placement.
        port: which input port solve() would be called with — "top" selects
            geometry.input_ys' last (highest-y) entry, "bot" its first
            (lowest-y) entry (falls back to a single port at y=0.0 if the
            geometry has no input_ys attribute).
        show_pml: shade the PML absorbing-boundary regions.
        ax: existing Axes; creates a new figure if None.
        title: plot title (defaults to the geometry class name).
        output_path: full path for output PNG; directory created if needed.
        figsize: figure size (width, height) in inches.
        wg_layer/ox_layer: GDS layer tuples for waveguide/oxide polygons.
        **geometry_kwargs: forwarded to geometry.to_gds_eme() and (when the
            geometry defines it) geometry.source_aperture_at().

    Returns:
        (fig, ax)
    """
    _wl = getattr(solver, "waveguide_length", 0.5)
    dpml = getattr(solver, "dpml", 0.0)
    straight_length = _wl + dpml
    oxide_width = solver.sim_width + 2 * dpml

    gds_polys = geometry.to_gds_eme(
        straight_length=straight_length, oxide_width=oxide_width, **geometry_kwargs
    )

    fig, ax = _new_ax(ax, figsize)
    _draw_polygons(ax, gds_polys, wg_layer, ox_layer)

    total_gds = geometry.total_length + 2 * straight_length
    input_ys = getattr(geometry, "input_ys", [0.0])
    src_yo = input_ys[-1] if port == "top" else input_ys[0]
    # Same source-width rule as FDTDSolver.solve(): a geometry-provided
    # launch window (source_aperture_at() non-None at the source's device-
    # frame x, clamped to the non-PML domain) shrinks the source + norm
    # monitor; otherwise they span the full sim_width.
    src_width = None
    _aperture_at = getattr(geometry, "source_aperture_at", None)
    if _aperture_at is not None:
        x_src_dev = -(geometry.total_length / 2 + _wl) + solver.source_margin
        src_aperture = _aperture_at(x_src_dev, **geometry_kwargs)
        if src_aperture is not None:
            src_width = min(src_aperture, solver.sim_width - 2.0 * abs(src_yo))
    _draw_fdtd_monitors(
        ax, dpml, total_gds, oxide_width,
        source_margin=solver.source_margin, monitor_margin=solver.monitor_margin,
        output_monitor_centers=getattr(solver, "output_monitor_centers", None),
        output_monitor_widths=getattr(solver, "output_monitor_widths", None),
        sim_width=solver.sim_width,
        src_yo=src_yo,
        src_width=src_width,
    )
    if show_pml:
        _draw_pml(ax, dpml, total_gds, oxide_width)

    _finalize(fig, ax, title, geometry, output_path)
    return fig, ax


def plot_section_fdtd_geometry(
    geometry,
    solver,
    x0: float,
    x1: float,
    n_in: int = 1,
    n_out: int = 1,
    stub_left: float = 0.0,
    stub_right: float = 0.0,
    monitor_center: Optional[float] = None,
    monitor_size: Optional[float] = None,
    interface_monitors: Optional[List[Tuple[float, str]]] = None,
    show_pml: bool = True,
    ax=None,
    title: Optional[str] = None,
    output_path: Optional[str] = None,
    figsize: Tuple[float, float] = (6, 3.5),
    wg_layer: Tuple[int, int] = (2, 0),
    ox_layer: Tuple[int, int] = (1, 0),
    **geometry_kwargs,
):
    """Top-down geometry plot for a single FDTDSolver.solve_section() call,
    standalone — no HybridSolver/section-index required.

    Same rendering `_plot_hybrid_sub_regions()`'s FDTD branch uses internally
    for a HybridSolver's own sub-region plots, factored out so a caller that
    invokes `solve_section()` directly — e.g. `section_mode_beta()`'s
    per-length calibration windows, which aren't tied to any HybridSolver
    boundary index — can visualize exactly what it simulated without faking
    up a single-section HybridSolver just to get a plot.

    Args:
        geometry: DeviceGeometry instance.
        solver: FDTDSolver whose .sim_width/.dpml/.source_margin/
            .monitor_margin/.half_domain_aperture drive the layout.
        x0, x1: device-frame section bounds — same meaning as
            FDTDSolver.solve_section()'s own x0/x1.
        n_in, n_out: same as solve_section()'s n_modes_in/n_modes_out — only
            affects source-vs-monitor-side layout via section_fdtd_layout(),
            not a drawn mode count.
        stub_left, stub_right, monitor_center, monitor_size: same as
            solve_section()'s own args.
        interface_monitors: explicit ``[(x_local, label), ...]`` interface
            monitors, replacing the single one this draws by default —
            what a MERGED multi-length domain
            (`FDTDSolver.solve_section_multi_length()`) actually contains:
            one shared source, one monitor per length. Positions belong in
            the section-local frame and should come from
            `merged_mea_positions()`, not be re-derived by the caller.
        show_pml: shade the PML absorbing-boundary rectangles.
        ax: existing Axes; creates a new figure if None.
        title: plot title (defaults to a generic section description).
        output_path: full path for output PNG; directory created if needed.
        figsize, wg_layer, ox_layer: as plot_fdtd_geometry().
        **geometry_kwargs: forwarded to geometry.section_gds() and
            section_fdtd_layout()'s aperture lookups.

    Returns:
        (fig, ax)
    """
    dpml = getattr(solver, "dpml", 0.0) or 0.0
    oxide_width = solver.sim_width + 2 * dpml
    L = (x1 - x0) + stub_left + stub_right

    # solve_section() draws stub_left/right + dpml (not just stub_left/
    # right) so the waveguide continues through the full PML thickness
    # instead of stopping at its inner edge — match that here so the plot
    # shows what's actually simulated (see _plot_hybrid_sub_regions).
    draw_stub_left = stub_left + dpml
    draw_stub_right = stub_right + dpml
    gds_polys = geometry.section_gds(
        x0, x1, oxide_width, stub_left=draw_stub_left, stub_right=draw_stub_right,
        **geometry_kwargs
    )

    fig, ax = _new_ax(ax, figsize)
    _draw_polygons(ax, gds_polys, wg_layer, ox_layer, x_shift=-dpml)

    _draw_section_fdtd_markers(
        ax, geometry, solver, x0, x1, stub_left, stub_right,
        oxide_width, n_in, n_out, geometry_kwargs,
        monitor_center=monitor_center, monitor_size=monitor_size,
        interface_monitors=interface_monitors,
    )
    margin = max(0.5, 0.08 * L)
    if show_pml:
        margin = max(margin, dpml + 0.1)
        y0, y1 = -oxide_width / 2, oxide_width / 2
        _draw_pml_rect(ax, -dpml, 0.0, y0, y1, label="PML")
        _draw_pml_rect(ax, L, L + dpml, y0, y1, label="PML")
        _draw_pml_rect(ax, -dpml, L + dpml, y0, y0 + dpml)
        _draw_pml_rect(ax, -dpml, L + dpml, y1 - dpml, y1)

    ax.set_xlim(-margin, L + margin)
    ax.set_ylim(-oxide_width / 2, oxide_width / 2)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (µm, section-local)")
    ax.set_ylabel("y (µm)")
    ax.set_title(title or f"{type(solver).__name__} section [{x0:.2f}, {x1:.2f}] µm")
    ax.grid(False)

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        fig.savefig(output_path, dpi=600, bbox_inches="tight")
    return fig, ax


# ---------------------------------------------------------------------------
# Hybrid
# ---------------------------------------------------------------------------

def plot_hybrid_geometry(
    geometry,
    solver,
    in_modes: int = 1,
    out_modes: int = 1,
    ax=None,
    title: Optional[str] = None,
    output_path: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 4),
    wg_layer: Tuple[int, int] = (2, 0),
    ox_layer: Tuple[int, int] = (1, 0),
    sub_region_dir: Optional[str] = None,
    sub_region_pad: Optional[float] = None,
    sub_region_show_pml: bool = False,
    sub_region_show_eme_detail: bool = False,
    sub_region_configs: Optional[List[Tuple[str, Union[Dict[str, Any], List[Dict[str, Any]]]]]] = None,
    sub_region_sections: Optional[List[int]] = None,
    **geometry_kwargs,
):
    """Top-down geometry plot for a HybridSolver.

    Draws the device polygons plus one shaded/labeled block per section
    (matching solver.solvers/solver.boundaries) and the section boundary
    lines. PML is intentionally *not* drawn on this overview — each section
    is its own independent little simulation with its own PML (including at
    interior cuts, see sub_region_show_pml below), so there's no single
    "domain edge" PML for the assembled device as a whole; use
    sub_region_show_pml on the per-section plots instead.

    Args:
        geometry: DeviceGeometry instance.
        solver: HybridSolver.
        in_modes/out_modes: mode counts at the true left/right physical
            ports — used both to label sub-region plots and (for the sub-
            region source/monitor placement) to pick the correct excitation
            side; pass the same values you'll use in solver.solve().
        ax: existing Axes; creates a new figure if None.
        title: plot title (defaults to the geometry class name).
        output_path: full path for the overview PNG; directory created if needed.
        figsize: overview figure size (width, height) in inches.
        wg_layer/ox_layer: GDS layer tuples for waveguide/oxide polygons.
        sub_region_dir: if given, also render zoomed-in PNGs per section
            into this directory (created if needed) — an FDTD section gets
            plot_fdtd_geometry-style source/monitor markers (using the
            reciprocity rule FDTDSolver.solve_section() itself uses: excite
            from whichever side has fewer modes, so an output-like section
            correctly shows its source on the *output* side, labeled
            "source (reverse)"), sized/centered via geometry.width_at()/
            center_at() exactly like solve_section() itself; an EME section
            gets a single shaded region labeled with its mode count (or,
            with sub_region_show_eme_detail, its num_modes/num_cells
            breakdown).
        sub_region_pad: x-margin (µm) shown around each section in the
            sub-region plots; defaults to max(0.5, 8% of the section width).
        sub_region_show_pml: also shade each section's own PML (all four
            sides — left/right at its stub_left/stub_right + dpml, top/
            bottom at +-oxide_width/2) in the sub-region plots. solve_section()
            wraps *every* section (interior cuts included) in its own
            dpml-thick PML, so an interior section's left/right strip is
            drawn immediately at its cut — this can visually overlap a
            neighboring section's own crop by one PML width, which is an
            accurate reflection of the two sections being solved as fully
            independent small simulations, not a rendering bug.
        sub_region_show_eme_detail: for EME sections, also break the single
            shaded region down into its sub-solver's num_modes/num_cells/
            lengths — one shaded sub-region per lengths[i] (rescaled to
            fill the section's real span, since a HybridSolver EME sub-
            solver's own lengths are frequently placeholders — see the
            hybrid tutorials — so only their *relative* proportions are
            meaningful) labeled with num_modes[i], further subdivided by
            thin dotted lines into num_cells[i] slices (the individual
            meow cells EMESolver._build_cells() would mesh that region
            into in the standalone, non-hybrid solve/solve_sweep() path).
        sub_region_configs: list of (name, geometry_kwargs) or
            (name, geometry_kwargs, mirror_pairs) entries, where
            geometry_kwargs is either a single dict (applied to every
            section), a sparse {section_idx: kwargs} dict mapping only the
            sections that need kwargs (unnamed sections get {}; a dict with
            int keys is always read as this sparse form), or a list of
            len(solver.solvers) dicts (one per section) — all matching
            HybridSolver.solve()'s own geometry_kwargs convention. E.g.
            [("through", {2: {"output_idx": 0}}),
            ("cross", {2: {"output_idx": 1}})] for a device whose
            physical output port varies per hybrid.solve() call. A section
            is only re-rendered once per *distinct* kwargs value it actually
            takes across configs — e.g. an input section whose kwargs never
            change across configs is rendered once, unsuffixed, not once
            per config. The optional third element declares that config's
            mirror_pairs ({section: partner}, matching HybridSolver's own
            parameter) when it differs per config (e.g. one output port
            mirror-derived, the other genuinely solved — two HybridSolver
            instances); it defaults to solver.mirror_pairs. A section that
            a config mirror-derives is NOT rendered for that config — the
            sub-region plots stay 1:1 with the simulations actually run —
            and is instead hatched and labeled in the overview plot.
            Defaults to a single unnamed entry using **geometry_kwargs
            below, matching the plain (non-branching) case.
        sub_region_sections: indices into solver.solvers to render sub-region
            plots for; None (default) renders every section not excluded by
            a mirror declaration (see sub_region_configs above — mirror-
            derived sections are skipped automatically, so this is only
            needed for sections skipped for some *other* reason, e.g. a
            transfer matrix reused by a symmetry argument hand-rolled
            outside HybridSolver).
        **geometry_kwargs: forwarded to geometry.to_gds_eme() for the
            overview plot, and used as the default sub_region_configs entry.

    Returns:
        (fig, ax)
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    _ref = solver.solvers[0]
    _wl = getattr(_ref, "waveguide_length", 0.5)
    dpml = getattr(_ref, "dpml", 0.0)
    straight_length = _wl + dpml
    oxide_width = getattr(_ref, "sim_width", 5.0) + 2 * dpml

    gds_polys = geometry.to_gds_eme(
        straight_length=straight_length, oxide_width=oxide_width, **geometry_kwargs
    )

    fig, ax = _new_ax(ax, figsize)
    _draw_polygons(ax, gds_polys, wg_layer, ox_layer)

    # device frame (x=0 at geometry centre) -> GDS plot frame
    x_shift = geometry.total_length / 2 + straight_length
    bounds_gds = [b + x_shift for b in solver.boundaries]
    y0_h, y1_h = -(oxide_width / 2), +(oxide_width / 2)
    palette = plt.cm.Set2.colors

    # Normalize configs up front (not just for the sub-region plots): each
    # entry is (name, geometry_kwargs) or (name, geometry_kwargs,
    # mirror_pairs), and the mirror declarations also drive the overview's
    # mirrored-region annotations below. mirror_pairs defaults to the
    # solver's own, so passing a mirror_pairs-carrying HybridSolver alone
    # (no explicit configs) annotates and skips correctly too.
    K = len(solver.solvers)
    default_mirror = dict(getattr(solver, "mirror_pairs", {}) or {})
    configs = sub_region_configs if sub_region_configs is not None else [("", geometry_kwargs)]
    normalized = []
    for entry in configs:
        name, cfg = entry[0], entry[1]
        mirror = dict(entry[2]) if len(entry) > 2 else default_mirror
        if isinstance(cfg, dict) and any(isinstance(key, int) for key in cfg):
            # Sparse per-section mapping {section_idx: kwargs}, same
            # convention HybridSolver._kwargs_per_section() accepts —
            # unnamed sections get {}.
            if not all(isinstance(key, int) for key in cfg):
                raise ValueError(
                    f"sub_region_configs entry {name!r} mixes int section "
                    f"indices with other keys; use either a sparse "
                    f"{{section_idx: kwargs}} mapping or a single plain "
                    f"kwargs dict, not both."
                )
            bad = [key for key in cfg if not (0 <= key < K)]
            if bad:
                raise ValueError(
                    f"sub_region_configs entry {name!r}: section indices "
                    f"must be in [0, {K - 1}], got {sorted(bad)}."
                )
            per_section = [dict(cfg.get(i, {})) for i in range(K)]
        elif isinstance(cfg, dict):
            per_section = [cfg] * K
        else:
            per_section = list(cfg)
        if len(per_section) != K:
            raise ValueError(
                f"sub_region_configs entry {name!r} has {len(per_section)} "
                f"geometry_kwargs entries, expected a single dict, a sparse "
                f"{{section_idx: kwargs}} mapping, or a list of "
                f"len(solver.solvers) = {K}."
            )
        normalized.append((name, per_section, mirror))
    n_configs = len(normalized)

    # section index -> {partner index: [config names that mirror it]}
    mirror_ann = {}
    for name, _per_section, mirror in normalized:
        for mk, mj in mirror.items():
            mirror_ann.setdefault(mk, {}).setdefault(mj, []).append(name)

    for k, sub_solver in enumerate(solver.solvers):
        xa, xb = bounds_gds[k], bounds_gds[k + 1]
        color = palette[k % len(palette)]
        ax.add_patch(Rectangle((xa, y0_h), xb - xa, y1_h - y0_h, alpha=0.15, color=color, zorder=3))
        ax.text((xa + xb) / 2, y1_h * 0.55, type(sub_solver).__name__,
                ha="center", va="top", fontsize=7, color=color, zorder=10)
        if k in mirror_ann:
            ax.add_patch(Rectangle((xa, y0_h), xb - xa, y1_h - y0_h, fill=False,
                                   hatch="///", edgecolor=color, lw=0, alpha=0.5, zorder=4))
            parts = []
            for j, names in sorted(mirror_ann[k].items()):
                scope = "" if len(names) == n_configs else f" ({', '.join(n for n in names if n)})"
                parts.append(f"mirror of section {j + 1}{scope}")
            ax.text((xa + xb) / 2, y1_h * 0.30, "; ".join(parts),
                    ha="center", va="top", fontsize=6.5, color=color, zorder=10,
                    style="italic")

    for k in range(1, len(solver.solvers)):
        x_b = bounds_gds[k]
        label_b = f"{type(solver.solvers[k-1]).__name__}|{type(solver.solvers[k]).__name__}"
        ax.plot([x_b, x_b], [y0_h, y1_h], color="#333333", lw=1.2, linestyle="--",
                alpha=0.8, zorder=9)
        ax.text(x_b, y1_h * 0.80, label_b, ha="center", va="top",
                fontsize=7, color="#333333", zorder=10,
                bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                          alpha=0.75, edgecolor="none"))

    _finalize(fig, ax, title, geometry, output_path)

    if sub_region_dir:
        _plot_hybrid_sub_regions(
            geometry, solver, oxide_width, in_modes, out_modes,
            sub_region_dir, wg_layer, ox_layer, sub_region_pad,
            sub_region_show_pml, sub_region_show_eme_detail, normalized,
            sub_region_sections,
        )

    return fig, ax


def _draw_section_fdtd_markers(ax, geometry, sub_solver, x0, x1, stub_left, stub_right,
                                oxide_width, n_in, n_out, geometry_kwargs,
                                monitor_center=None, monitor_size=None,
                                interface_monitors=None):
    """Source/monitor markers for one HybridSolver FDTD section, in the
    section's own local frame — left edge of its section_gds()/
    section_prisms() domain (i.e. x0-stub_left) at x_local=0, matching the
    frame FDTDSolver.solve_section() itself builds its meep cell in.

    All position/size math (which side excites vs. measures, source/monitor
    apertures, the source_aperture_at()/half_domain_aperture priority order)
    comes from section_fdtd_layout() — the same function
    FDTDSolver.solve_section() calls — so this plot can't independently
    drift out of sync with what actually gets simulated; see that
    function's docstring in beef/solvers/fdtd.py for the full rules.
   
    `interface_monitors`, when given, REPLACES the single interface monitor
    with an explicit ``[(x_local, label), ...]`` set -- for a merged
    multi-length domain, which carries one interior monitor per length off a
    single shared source (`FDTDSolver.solve_section_multi_length()`). The
    caller supplies the whole set, boundary monitor included, so nothing is
    drawn twice; positions must come from `merged_mea_positions()` rather
    than be re-derived here, for the drift reason above.
    """
    dpml = getattr(sub_solver, "dpml", 0.0) or 0.0
    sim_width = getattr(sub_solver, "sim_width", oxide_width - 2 * dpml)

    layout = section_fdtd_layout(
        geometry, x0, x1, stub_left, stub_right, n_in, n_out, sim_width,
        getattr(sub_solver, "source_margin", 0.0),
        getattr(sub_solver, "monitor_margin", 0.25),
        getattr(sub_solver, "half_domain_aperture", False),
        geometry_kwargs, monitor_center=monitor_center, monitor_size=monitor_size,
    )

    y_top, y_bot = oxide_width / 2 - dpml, -(oxide_width / 2 - dpml)

    def _clamped_span(c, width):
        return max(c - width / 2, y_bot), min(c + width / 2, y_top)

    src_label = "source (reverse)" if layout.reverse else "source"
    y_lo, y_hi = _clamped_span(layout.src_c, layout.src_width)
    _draw_marker(ax, layout.src_x, y_lo, y_hi, src_label, "crimson", "-")
    _draw_marker(ax, layout.norm_x, y_lo, y_hi, "norm mon", "darkorange", "--")

    mea_lo, mea_hi = _clamped_span(layout.mea_c, layout.mea_w)
    if interface_monitors is None:
        _draw_marker(ax, layout.mea_x, mea_lo, mea_hi, "interface mon", "seagreen", "--")
    else:
        for _x, _label in interface_monitors:
            _draw_marker(ax, _x, mea_lo, mea_hi, _label, "seagreen", "--")


def _draw_pml_rect(ax, xa, xb, ya, yb, label=None):
    from matplotlib.patches import Rectangle

    ax.add_patch(Rectangle((xa, ya), xb - xa, yb - ya, alpha=0.18, color="#8B0000", zorder=2))
    if label:
        ax.text(
            (xa + xb) / 2, yb - 0.12 * (yb - ya), label,
            ha="center", va="top", fontsize=6.5, color="#8B0000", zorder=7,
            bbox=dict(boxstyle="round,pad=0.1", facecolor="white", alpha=0.65, edgecolor="none"),
        )


def _draw_eme_region_detail(ax, x0_g, x1_g, oxide_width, num_modes, num_cells, lengths):
    """Proportionally distribute a sub-solver's num_modes/num_cells/lengths
    across [x0_g, x1_g] — one shaded region per lengths[i] (rescaled to fill
    the section's real span; only the *relative* proportions of a
    HybridSolver EME sub-solver's lengths are meaningful, see the
    sub_region_show_eme_detail docstring above), labeled with num_modes[i]
    and subdivided into num_cells[i] thin dotted "slices" — the individual
    cells EMESolver._build_cells() would mesh that region into."""
    total_len = sum(lengths)
    scale = (x1_g - x0_g) / total_len if total_len > 0 else 0.0
    boundaries = [x0_g]
    for length in lengths:
        boundaries.append(boundaries[-1] + length * scale)

    labels = [f"{nm} modes" for nm in num_modes]
    _draw_eme_regions(ax, boundaries, oxide_width, labels=labels)

    for i, (xa, xb) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        n_cells = num_cells[i] if i < len(num_cells) else 1
        if n_cells > 1:
            cell_len = (xb - xa) / n_cells
            for j in range(1, n_cells):
                ax.axvline(xa + j * cell_len, color="#666666", lw=0.5, linestyle=":",
                           alpha=0.6, zorder=3.5)


def _section_kwargs_variants(normalized, k):
    """Group a section's geometry_kwargs across all configs by equality,
    returning [(names, kwargs), ...] — one entry per *distinct* kwargs value,
    each tagged with every config name that uses it. A section whose kwargs
    never change across configs comes back as a single entry naming every
    config (the caller renders that as one unsuffixed, config-agnostic plot);
    a section that genuinely branches (e.g. output_idx varying with the
    physical port) comes back as one entry per distinct value.

    A config that mirror-derives section k (k in its mirror_pairs) never
    simulates it, so it contributes no variant — the sub-region plots stay
    1:1 with the simulations actually run. If every config mirrors k, no
    variants come back and the section gets no sub-region plot at all (it's
    hatched/labeled in the overview instead). Because the caller compares
    len(names) against the total config count for the unsuffixed-name rule,
    a partially-mirrored section's surviving variant keeps its config
    suffix (e.g. section_2_..._cross.png), which is the informative name."""
    variants = []  # [[names], kwargs]
    for name, per_section, mirror in normalized:
        if k in mirror:
            continue
        kwargs_k = per_section[k]
        match = next((v for v in variants if v[1] == kwargs_k), None)
        if match is None:
            variants.append([[name], kwargs_k])
        else:
            match[0].append(name)
    return variants


def _plot_hybrid_sub_regions(geometry, solver, oxide_width, in_modes, out_modes,
                              output_dir, wg_layer, ox_layer, pad,
                              show_pml, show_eme_detail, normalized,
                              sections=None):
    """Render each section in its *own* local frame via geometry.section_gds()
    (left edge of that section's own stubbed domain at x_local=0) — the same
    geometry-construction call FDTDSolver.solve_section()/EMESolver.
    solve_section() themselves make, rather than cropping a view out of the
    whole-device to_gds_eme() drawing. This is what makes the "both arms
    present" fix visible here for free (section_prisms()/section_gds() now
    always draw every arm) and keeps the PML/marker placement pixel-exact
    to what solve_section() actually simulates."""
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)
    n_sections = len(solver.solvers)
    n_configs = len(normalized)

    # Same n_in/n_out bookkeeping as HybridSolver.solve() itself.
    n_in_list = [in_modes] + list(solver.interface_modes)
    n_out_list = list(solver.interface_modes) + [out_modes]

    for k, sub_solver in enumerate(solver.solvers):
        if sections is not None and k not in sections:
            continue
        x0, x1 = solver.boundaries[k], solver.boundaries[k + 1]
        n_in, n_out = n_in_list[k], n_out_list[k]
        is_eme = hasattr(sub_solver, "lengths")
        is_fdtd = hasattr(sub_solver, "sim_width")

        # Same stub_left/stub_right rule HybridSolver.solve() itself uses:
        # only the true first/last section gets the access-waveguide stub
        # (EMESolver has no waveguide_length, so an EME section's stub is
        # always 0 regardless of position).
        wl = getattr(sub_solver, "waveguide_length", 0.0)
        stub_left = wl if k == 0 else 0.0
        stub_right = wl if k == n_sections - 1 else 0.0
        L = (x1 - x0) + stub_left + stub_right
        dpml = getattr(sub_solver, "dpml", 0.0) or 0.0

        for names, kwargs_k in _section_kwargs_variants(normalized, k):
            # This kwargs value applies to every config -> config-agnostic,
            # render once, unsuffixed (e.g. an input/body section that
            # doesn't vary with which output port is being solved for).
            config_name = "" if len(names) == n_configs else "_".join(names)

            # A section's kwargs dict may mix true geometry kwargs with
            # FDTDSolver.solve_section()'s own monitor_center/monitor_size
            # parameters (they ride in the same per-section dict, exactly as
            # HybridSolver.solve() forwards them — there they bind to
            # solve_section's named parameters rather than reaching the
            # geometry hooks). Split them off here the same way, before
            # section_gds()/width_at()/center_at() see the dict.
            kwargs_k = dict(kwargs_k)
            monitor_center = kwargs_k.pop("monitor_center", None)
            monitor_size = kwargs_k.pop("monitor_size", None)

            # FDTDSolver.solve_section() draws stub_left/right + dpml (not
            # just stub_left/right) so the waveguide continues through the
            # full PML thickness instead of stopping at its inner edge (see
            # its docstring in beef/solvers/fdtd.py) — match that here so
            # the plot shows what's actually simulated. L (used below for
            # marker/axis placement) still keys off the true, unextended
            # stub_left/stub_right, matching the solver's own bookkeeping.
            draw_stub_left = stub_left + dpml if is_fdtd else stub_left
            draw_stub_right = stub_right + dpml if is_fdtd else stub_right
            gds_polys = geometry.section_gds(
                x0, x1, oxide_width, stub_left=draw_stub_left, stub_right=draw_stub_right, **kwargs_k
            )

            margin = pad if pad is not None else max(0.5, 0.08 * L)
            if show_pml and is_fdtd:
                margin = max(margin, dpml + 0.1)

            fig, ax = plt.subplots(figsize=(6, 3.5))
            # gds_polys was built with draw_stub_left/right (extended by dpml
            # for an FDTD section), which shifts its own x_local=0 left by
            # dpml relative to the true stub_left/right frame markers/PML
            # rects use below — x_shift=-dpml re-aligns it (see
            # _draw_polygons' docstring).
            _draw_polygons(ax, gds_polys, wg_layer, ox_layer, x_shift=-dpml if is_fdtd else 0.0)

            if is_eme:
                if show_eme_detail:
                    _draw_eme_region_detail(
                        ax, 0.0, L, oxide_width,
                        sub_solver.num_modes, sub_solver.num_cells, sub_solver.lengths,
                    )
                else:
                    label = f"{n_in} modes" if n_in == n_out else f"{n_in}→{n_out} modes"
                    _draw_eme_regions(ax, [0.0, L], oxide_width, labels=[label])
                kind = "EME"
            elif is_fdtd:
                _draw_section_fdtd_markers(ax, geometry, sub_solver, x0, x1, stub_left, stub_right,
                                            oxide_width, n_in, n_out, kwargs_k,
                                            monitor_center=monitor_center,
                                            monitor_size=monitor_size)
                if show_pml:
                    y0, y1 = -oxide_width / 2, oxide_width / 2
                    _draw_pml_rect(ax, -dpml, 0.0, y0, y1, label="PML")           # left
                    _draw_pml_rect(ax, L, L + dpml, y0, y1, label="PML")         # right
                    _draw_pml_rect(ax, -dpml, L + dpml, y0, y0 + dpml)           # bottom
                    _draw_pml_rect(ax, -dpml, L + dpml, y1 - dpml, y1)           # top
                kind = "FDTD"
            else:
                kind = type(sub_solver).__name__

            # adjustable="box" (not "datalim"): keep these explicit crop
            # limits fixed and let the axes box itself shrink/grow to match
            # the aspect ratio, rather than "equal" silently re-expanding
            # the data limits back out to the full (uncropped) polygon extent.
            ax.set_xlim(-margin, L + margin)
            ax.set_ylim(-oxide_width / 2, oxide_width / 2)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("x (µm, section-local)")
            ax.set_ylabel("y (µm)")
            title_suffix = f" [{config_name}]" if config_name else ""
            ax.set_title(f"Section {k}: {type(sub_solver).__name__} ({kind}){title_suffix}")
            ax.grid(False)

            fname = f"section_{k}_{type(sub_solver).__name__}"
            if config_name:
                fname += f"_{config_name}"
            fig.savefig(os.path.join(output_dir, f"{fname}.png"), dpi=600, bbox_inches="tight")
            plt.close(fig)
