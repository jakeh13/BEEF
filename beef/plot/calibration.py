"""Plots for `FDTDSolver.section_mode_beta_details()`'s β calibration.

Two things worth looking at after a body calibration, neither of which the
returned β alone can show:

- **What was actually simulated** — the per-length FDTD windows
  `section_mode_beta()` drew inside the body span. Rendered by
  `plot_beta_calibration_geometry()` via the existing
  `plot_section_fdtd_geometry()`, whose own docstring already names these
  calibration windows as its motivating case.
- **Whether the phase really is linear in length** — the fit behind every
  (wavelength, mode) cell. The residual guard reduces that to one number per
  cell and raises only past 0.1 rad; a fit sitting at 0.09 rad passes
  silently and looks identical to one at 1e-5 in the returned β. These plots
  are how that difference becomes visible.

All of these take an already-computed `BetaCalibration`, so they cost
nothing but disk. Getting one cheaply for a device already simulated: call
`section_mode_beta_details(..., save_dir=..., cache_only=True)`, which reads
the on-disk cache and raises loudly rather than silently re-simulating.

Colour: mode index is an *ordered* quantity, not a set of unrelated
categories, so series are stepped along one perceptually-uniform sequential
ramp (`viridis`) rather than given cycled categorical hues — at the
`N_MODES=14` these calibrations actually run at, cycled hues would repeat
and stop identifying anything.
"""
from __future__ import annotations

import os
from typing import List, Optional, Sequence, Tuple

import numpy as np
from matplotlib import pyplot as plt

from ..solvers.fdtd import _BETA_FIT_RESIDUAL_MAX, _fit_line
from .geometry import plot_section_fdtd_geometry


def _mode_colors(n_modes: int):
    """One colour per mode, stepped along a sequential ramp (see module doc)."""
    return plt.get_cmap("viridis")(np.linspace(0.0, 0.88, max(n_modes, 1)))


def plot_beta_calibration_fits(
    cal,
    output_dir: str,
    modes: Optional[Sequence[int]] = None,
    filename_prefix: str = "fit_mode",
    ncols: int = 4,
    panel_size: Tuple[float, float] = (3.0, 2.4),
) -> List[str]:
    """One PNG per mode: the phase-vs-length fit at every wavelength.

    Each panel is one (wavelength, mode) cell — the measured points are
    `arg(T_ii)` *unwrapped by the turn counts the fit chose*
    (``phase + 2π·turns``), not the raw wrapped values, because the raw
    values jump by 2π between lengths and would show no line at all. The
    fitted line is the β that calibration returns for that cell; the
    annotation carries β, the max residual, and the branch margin, with the
    residual flagged when it exceeds `_BETA_FIT_RESIDUAL_MAX`.

    Args:
        cal: `BetaCalibration` from `FDTDSolver.section_mode_beta_details()`.
        output_dir: directory for the PNGs; created if needed.
        modes: mode indices to plot; None (default) plots all of them.
        filename_prefix: output files are ``<prefix>_<mode:02d>.png``.
        ncols: panels per row in each mode's wavelength grid.
        panel_size: (width, height) inches per panel.

    Returns:
        The list of paths written, in mode order.
    """
    os.makedirs(output_dir, exist_ok=True)
    n_wl, n_modes = cal.beta.shape
    mode_list = list(range(n_modes)) if modes is None else list(modes)
    L = np.asarray(cal.lengths, dtype=float)
    L_line = np.linspace(0.0, L.max() * 1.05, 100)
    written: List[str] = []

    for m in mode_list:
        ncol = min(ncols, n_wl)
        nrow = int(np.ceil(n_wl / ncol))
        fig, axes = plt.subplots(
            nrow, ncol, figsize=(panel_size[0] * ncol, panel_size[1] * nrow),
            squeeze=False,
        )
        for i in range(n_wl):
            ax = axes[i // ncol][i % ncol]
            fit = cal.fits[i][m]
            unwrapped = cal.phases[:, i, m] + 2.0 * np.pi * np.asarray(fit.turns)
            slope, intercept = _fit_line(L, unwrapped)

            ax.plot(L_line, slope * L_line + intercept, "-", lw=1.2,
                    color="0.55", zorder=1,
                    label=f"fit  β={fit.beta:.4f}")
            ax.plot(L, unwrapped, "o", ms=5, color="#1f6f8b", zorder=2,
                    label="measured")

            over = fit.resid > _BETA_FIT_RESIDUAL_MAX

            # The line spans tens of radians and the deviations from it are
            # ~1e-4, so on the main axes every fit looks perfectly straight
            # whether it is or not — and "is the phase linear in L" is the
            # only question the residual guard actually asks. The inset is
            # where that becomes visible; scaled to the guard threshold when
            # the fit is anywhere near it, else to the data.
            resid_pts = unwrapped - (slope * L + intercept)
            inset = ax.inset_axes([0.56, 0.09, 0.40, 0.26])
            inset.axhline(0.0, color="0.6", lw=0.7)
            inset.plot(L, resid_pts, "o-", ms=3, lw=0.8,
                       color="#a11" if over else "#1f6f8b")
            span = max(float(np.max(np.abs(resid_pts))), 1e-12)
            inset.set_ylim(-1.6 * span, 1.6 * span)
            inset.set_title("residual (rad)", fontsize=6, pad=2)
            inset.tick_params(labelsize=5, length=2, pad=1)
            inset.set_xticks([])
            inset.patch.set_alpha(0.85)
            # matplotlib parks the shared "1e-5" exponent above the axes at
            # full size, where it collides with the inset title.
            inset.yaxis.get_offset_text().set_fontsize(5)
            ax.set_title(
                f"λ={cal.wavelengths[i]:.4f} µm\n"
                f"β={fit.beta:.4f}  resid={fit.resid:.2e}"
                f"{' ⚠' if over else ''}  margin={fit.margin:.3f}",
                fontsize=8,
            )
            ax.grid(True, alpha=0.25, lw=0.6)
            ax.tick_params(labelsize=7)
            if i // ncol == nrow - 1:
                ax.set_xlabel("section length L (µm)", fontsize=8)
            if i % ncol == 0:
                ax.set_ylabel("unwrapped arg($T_{ii}$) (rad)", fontsize=8)

        for k in range(n_wl, nrow * ncol):        # blank any unused cells
            axes[k // ncol][k % ncol].axis("off")

        # One legend for the figure: every panel draws the same two series,
        # so repeating it per panel would be n_wl copies of the same key.
        handles, labels = axes[0][0].get_legend_handles_labels()
        fig.legend([handles[1], handles[0]], [labels[1], "fit"],
                   loc="upper right", fontsize=8, framealpha=0.9)
        failed = bool(getattr(cal, "guard_failures", ()))
        fig.suptitle(
            f"β fit vs section length — mode {m}   "
            f"(span [{cal.x0:.3f}, {cal.x1:.3f}] µm, "
            f"{len(cal.lengths)} lengths, "
            f"{'merged domain' if cal.merged else 'per-length domains'})"
            + ("   ⚠ GUARD FAILED — inspection only" if failed else ""),
            fontsize=10,
            color="#a11" if failed else "black",
        )
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        path = os.path.join(output_dir, f"{filename_prefix}_{m:02d}.png")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        written.append(path)

    return written


def plot_beta_calibration_summary(
    cal,
    output_dir: str,
    filename: str = "calibration_summary.png",
) -> str:
    """Whole-calibration overview: β and seed-excess vs wavelength, plus the
    residual and branch-margin grids.

    Four panels, because they answer four different questions: is β smooth in
    wavelength (top left), is the FDTD-vs-MPB excess in the ~2% family this
    calibration exists to measure (top right), and do any (wavelength, mode)
    cells sit near either guard threshold (bottom row, where the colour scale
    is what makes an outlier cell visible at a glance).

    Args:
        cal: `BetaCalibration` from `FDTDSolver.section_mode_beta_details()`.
        output_dir: directory for the PNG; created if needed.
        filename: output filename.

    Returns:
        The path written.
    """
    os.makedirs(output_dir, exist_ok=True)
    n_wl, n_modes = cal.beta.shape
    wl = np.asarray(cal.wavelengths, dtype=float)
    colors = _mode_colors(n_modes)

    resid = np.array([[cal.fits[i][m].resid for m in range(n_modes)]
                      for i in range(n_wl)])
    margin = np.array([[cal.fits[i][m].margin for m in range(n_modes)]
                       for i in range(n_wl)])
    excess = np.array([[cal.fits[i][m].seed_excess for m in range(n_modes)]
                       for i in range(n_wl)])

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax_beta, ax_excess, ax_resid, ax_margin = axes.ravel()

    for m in range(n_modes):
        ax_beta.plot(wl, cal.beta[:, m], "-o", ms=3, lw=1.2, color=colors[m])
        ax_excess.plot(wl, 100.0 * excess[:, m], "-o", ms=3, lw=1.2, color=colors[m])
    ax_beta.set_xlabel("wavelength (µm)")
    ax_beta.set_ylabel("fitted β (rad/µm)")
    ax_beta.set_title("FDTD's own accumulated β")
    ax_beta.grid(True, alpha=0.25, lw=0.6)

    ax_excess.axhline(0.0, color="0.6", lw=0.8)
    ax_excess.set_xlabel("wavelength (µm)")
    ax_excess.set_ylabel("100·(β/β_seed − 1)  (%)")
    ax_excess.set_title("excess over the MPB seed — the discretization gap")
    ax_excess.grid(True, alpha=0.25, lw=0.6)

    # Mode identity via a colourbar rather than a 14-entry legend box: the
    # index is ordered, so the ramp itself is the key (see module docstring).
    sm = plt.cm.ScalarMappable(
        cmap="viridis", norm=plt.Normalize(vmin=0, vmax=max(n_modes - 1, 1)))
    # Integer ticks only — "mode 1.75" does not exist. Thinned so a 14-mode
    # calibration doesn't stack unreadable labels.
    mode_ticks = np.arange(0, n_modes, max(1, int(np.ceil(n_modes / 8))))
    for ax in (ax_beta, ax_excess):
        fig.colorbar(sm, ax=ax, label="mode index", pad=0.02, ticks=mode_ticks)

    for ax, grid, label, cmap in (
        (ax_resid, resid, f"fit residual (rad) — guard at {_BETA_FIT_RESIDUAL_MAX}", "magma_r"),
        (ax_margin, margin, "branch margin — guard at 0.9", "magma_r"),
    ):
        mesh = ax.pcolormesh(np.arange(n_modes + 1) - 0.5,
                             np.arange(n_wl + 1) - 0.5,
                             grid, cmap=cmap, shading="flat")
        fig.colorbar(mesh, ax=ax, label=label, pad=0.02)
        ax.set_xlabel("mode index")
        ax.set_ylabel("wavelength (µm)")     # the ticks below are µm, not indices
        ax.set_xticks(mode_ticks)            # integers only, as above
        ax.set_yticks(range(n_wl))
        ax.set_yticklabels([f"{v:.4f}" for v in wl], fontsize=7)
        ax.set_title(label.split(" — ")[0])

    # A calibration returned with failed guards (strict=False) must never
    # render as if it were usable — say so on the figure itself, not only in
    # the log the plot outlives.
    failed = bool(getattr(cal, "guard_failures", ()))
    fig.suptitle(
        f"β calibration summary — span [{cal.x0:.3f}, {cal.x1:.3f}] µm, "
        f"lengths {[f'{v:.3f}' for v in cal.lengths]} µm, "
        f"{'merged domain' if cal.merged else 'per-length domains'}"
        + ("\n⚠ GUARD FAILED — inspection only, this β is not usable for "
           "propagation" if failed else ""),
        fontsize=11,
        color="#a11" if failed else "black",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = os.path.join(output_dir, filename)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_beta_calibration_geometry(
    cal,
    geometry,
    solver,
    output_dir: str,
    n_modes: Optional[int] = None,
    monitor_center: Optional[float] = None,
    monitor_size: Optional[float] = None,
    filename_prefix: str = "length",
    **geometry_kwargs,
) -> List[str]:
    """One top-down FDTD geometry plot per calibration length.

    Draws exactly the windows `section_mode_beta()` simulated — each length
    centred on the body midpoint, i.e. ``[mid − L/2, mid + L/2]`` — using the
    same `plot_section_fdtd_geometry()` the hybrid examples already use for
    their own section plots.

    Note this renders the *per-length* windows, which is what the default
    (and every `save_dir`-passing production caller's) path actually
    simulates. A merged-domain calibration instead builds ONE domain sized to
    the longest length with interior monitors at each shorter one; when
    `cal.merged` is True these plots therefore show the equivalent geometry
    per length rather than the single shared domain — flagged in each title
    rather than silently implied.

    Args:
        cal: `BetaCalibration` from `FDTDSolver.section_mode_beta_details()`.
        geometry: the DeviceGeometry that was calibrated.
        solver: the FDTDSolver that ran it (drives the drawn layout).
        output_dir: directory for the PNGs; created if needed.
        n_modes: mode count for the source-vs-monitor side layout; defaults
            to the calibration's own (`cal.beta.shape[1]`), which is what
            `section_mode_beta()` passed as both n_in and n_out.
        monitor_center, monitor_size: as `section_mode_beta()`; pass the same
            values that call used so the drawn monitor matches.
        filename_prefix: output files are ``<prefix>_<index>_L<length>.png``.
        **geometry_kwargs: forwarded to `plot_section_fdtd_geometry()`.

    Returns:
        The list of paths written, in calibration-length order.
    """
    os.makedirs(output_dir, exist_ok=True)
    n = cal.beta.shape[1] if n_modes is None else n_modes
    written: List[str] = []

    for j, L in enumerate(cal.lengths):
        x0, x1 = cal.mid - L / 2.0, cal.mid + L / 2.0
        path = os.path.join(output_dir, f"{filename_prefix}_{j}_L{L:.4f}.png")
        plot_section_fdtd_geometry(
            geometry, solver, x0, x1,
            n_in=n, n_out=n,
            monitor_center=0.0 if monitor_center is None else monitor_center,
            monitor_size=solver.sim_width if monitor_size is None else monitor_size,
            title=(f"β calibration length {j}: L={L:.4f} µm  "
                   f"[{x0:.3f}, {x1:.3f}]"
                   + ("  (merged run — per-length equivalent shown)"
                      if cal.merged else "")),
            output_path=path,
            **geometry_kwargs,
        )
        written.append(path)

    return written


def plot_beta_calibration(
    cal,
    geometry,
    solver,
    output_dir: str,
    monitor_center: Optional[float] = None,
    monitor_size: Optional[float] = None,
    **geometry_kwargs,
) -> dict:
    """Everything above, into ``<output_dir>/geometry/`` and ``.../fits/``.

    The one call an example script wants: the simulated windows, the
    per-mode wavelength fits, and the whole-calibration summary, each in its
    own subfolder so a 14-mode calibration's ~15 PNGs don't land loose beside
    the device's own result plots.

    Args:
        cal: `BetaCalibration` from `FDTDSolver.section_mode_beta_details()`.
        geometry, solver: as `plot_beta_calibration_geometry()`.
        output_dir: parent directory; both subfolders are created under it.
        monitor_center, monitor_size: as `section_mode_beta()` — pass the
            same values that call used.
        **geometry_kwargs: forwarded to the geometry plots.

    Returns:
        ``{"geometry": [...], "fits": [...], "summary": path}`` — every file
        written, so a caller can log or assert on them.
    """
    geom_dir = os.path.join(output_dir, "geometry")
    fit_dir = os.path.join(output_dir, "fits")
    return {
        "geometry": plot_beta_calibration_geometry(
            cal, geometry, solver, geom_dir,
            monitor_center=monitor_center, monitor_size=monitor_size,
            **geometry_kwargs,
        ),
        "fits": plot_beta_calibration_fits(cal, fit_dir),
        "summary": plot_beta_calibration_summary(cal, fit_dir),
    }
