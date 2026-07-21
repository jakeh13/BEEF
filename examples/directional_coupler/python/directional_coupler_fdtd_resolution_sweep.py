"""
Tutorial: directional coupler — single-domain FDTD mesh-resolution convergence sweep.

The same 220 nm SOI, C-band device and single-domain method as
directional_coupler_fdtd.py (one broadband full-device meep simulation, fixed
coupling_length = 13.5 µm — the ~3 dB length for this cross-section), but
instead of running once at a fixed mesh resolution, this sweeps

    RESOLUTION = 20, 30, 40, ..., 120 px/µm   (11 values, 20:10:120)

to check mesh convergence: does the transmission spectrum actually stop
changing as the mesh gets finer, and if so, by what resolution? Everything
else (geometry, ports, wavelength sweep) is identical to
directional_coupler_fdtd.py — resolution is the only swept knob.

Outputs (example_output/directional_coupler/directional_coupler_fdtd_resolution_sweep/):

    spectrums/spectrum_R{res}.png   per-resolution transmission/IL spectrum
                                     (same layout as the single-resolution
                                     example's fdtd_results.png)
    data/R{res}.npz                 per-resolution raw sweep data
    resolution_convergence.png      aggregate two-panel convergence plot
    resolution_sweep.npz            aggregate grid over completed resolutions

Run:
    # one resolution per SLURM array task (see
    # run_directional_coupler_fdtd_resolution_sweep.batch):
    python /path/to/examples/directional_coupler/python/directional_coupler_fdtd_resolution_sweep.py <index 0..10>

    # or the whole sweep serially (slow — batch queue strongly preferred):
    python /path/to/examples/directional_coupler/python/directional_coupler_fdtd_resolution_sweep.py

Every invocation ends by (re)aggregating whatever resolutions have finished so
far, writing the aggregate plot/data atomically — so concurrent array tasks
can't corrupt them, the plot is usable (annotated "n/11 resolutions") while
the array is still draining, and the last task to finish leaves the complete
sweep behind. Re-running an index simply overwrites that resolution's files.

Note on cost: 3-D FDTD cost scales roughly as resolution^4 (cell count grows
~resolution^3, and the timestep shrinks ~1/resolution), so the res=120 task is
dramatically more expensive than the res=20/40 tasks — see the batch script's
own comment on walltime for this reason.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import DirectionalCoupler, FDTDSolver, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "directional_coupler", "directional_coupler_fdtd_resolution_sweep"))
spectrum_dir = os.path.join(output_dir, "spectrums")
data_dir = os.path.join(output_dir, "data")
for _d in (output_dir, spectrum_dir, data_dir):
    os.makedirs(_d, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters — identical to directional_coupler_fdtd.py,
# except RESOLUTION is swept instead of pinned.
# ---------------------------------------------------------------------------
WL_MIN = 1.50       # wavelength sweep start (µm)
WL_MAX = 1.60       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

COUPLING_LENGTH = 13.5   # µm, hard-coded — same ~3 dB length as directional_coupler_fdtd.py

RESOLUTION_SWEEP = np.arange(20.0, 120.0 + 1e-9, 10.0)   # px/µm — 20:10:120, 11 values

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)


def _tag(R: float) -> str:
    """Stable per-resolution file tag, e.g. R020 / R120."""
    return f"R{int(R):03d}"


def _atomic_replace_write(path: str, write_fn) -> None:
    """write_fn(tmp_path) then atomically rename onto path — concurrent array
    tasks all regenerate the aggregate files, and a plain overwrite could be
    read (or collide) half-written; os.replace is atomic on POSIX."""
    tmp = f"{path}.tmp.{os.getpid()}"
    write_fn(tmp)
    os.replace(tmp, path)


def make_coupler() -> DirectionalCoupler:
    return DirectionalCoupler(
        waveguide_width=0.45,             # single-mode strip width (µm)
        gap=0.25,                         # coupling-region edge-to-edge gap (µm)
        coupling_length=COUPLING_LENGTH,  # straight parallel interaction length (µm)
        bend_length=9.0,                  # horizontal extent of each S-bend (µm)
        waveguide_offset=3.0,             # far-end centre-to-centre port pitch (µm)
        bend_type="euler",
        waveguide_height=0.22,            # Si layer thickness (µm) — SOI_220
    )


def make_solver(R: float) -> FDTDSolver:
    return FDTDSolver(
        res=R,
        sim_width=6.0,
        sim_height=1.5,
        dpml=0.5,
        source_margin=0.0,
        monitor_margin=0.25,
        three_d=True,   # full 3-D FDTD: real 220 nm Si core in oxide, PML on all six faces
        num_input_modes=1,
        num_output_modes=1,
        # out_0 = bottom arm ("through", same arm as "in"), out_1 = top arm
        # ("cross") — half-domain windows, same as directional_coupler_fdtd.py.
        output_monitor_centers=[-1.5, +1.5],
        output_monitor_widths=[3.0, 3.0],
        material_stack=SOI_220,
        decay_tol=1e-11,
    )


# ---------------------------------------------------------------------------
# Per-resolution simulation: one broadband full-device FDTD run -> spectrum
# plot in spectrums/ + raw data in data/ for the aggregation step.
# ---------------------------------------------------------------------------
def simulate_resolution(R: float) -> None:
    dc = make_coupler()
    solver = make_solver(R)
    print(f"[sweep] resolution = {R:.0f} px/µm "
          f"(coupling_length={COUPLING_LENGTH:.1f} µm, single FDTD run, "
          f"{N_WL} wavelengths {WL_MIN:.2f}–{WL_MAX:.2f} µm)")

    smat = solver.solve(dc, wavelengths, port="bot", input_mode=1, save_dir=save_dir)
    smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength
    T_through = smat.transmission(port="out_0", in_mode=1, out_mode=1)
    T_cross = smat.transmission(port="out_1", in_mode=1, out_mode=1)
    T_total = T_through + T_cross

    np.savez(
        os.path.join(data_dir, f"{_tag(R)}.npz"),
        resolution=R, wavelengths=wavelengths,
        T_through=T_through, T_cross=T_cross,
    )

    idx = len(wavelengths) // 2
    print(f"[sweep] res={R:.0f} px/µm @ λ={wavelengths[idx]:.3f} µm: "
          f"through={T_through[idx]:.4f} cross={T_cross[idx]:.4f} "
          f"total={T_total[idx]:.4f} "
          f"ratio={T_cross[idx]/max(T_total[idx], 1e-20):.3f}")

    # -- Per-resolution spectrum plot (same layout as fdtd_results.png) ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax = axes[0]
    ax.plot(wavelengths, T_through, label="Through", color="steelblue")
    ax.plot(wavelengths, T_cross,   label="Cross",   color="darkorange", linestyle="--")
    ax.plot(wavelengths, T_total,   label="Total",   color="k", lw=1.5)
    ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Transmission")
    ax.set_title(f"Transmission vs Wavelength (resolution={R:.0f} px/µm)")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)

    ax = axes[1]
    for T, label, color, ls in [(T_through, "Through", "steelblue", "-"),
                                (T_cross,   "Cross",   "darkorange", "--"),
                                (T_total,   "Total",   "k", "-")]:
        ax.plot(wavelengths, 10 * np.log10(np.clip(T, 1e-20, None)),
                label=label, color=color, linestyle=ls, lw=1.5 if label == "Total" else None)
    ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Insertion Loss (dB)")
    ax.set_title(f"Insertion Loss vs Wavelength (resolution={R:.0f} px/µm)")
    ax.legend(fontsize=8)
    ax.grid(True)

    plt.suptitle("Directional Coupler — Single-Domain FDTD (resolution-sweep point)", fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(spectrum_dir, f"spectrum_{_tag(R)}.png"), dpi=300)
    plt.close("all")
    print(f"[sweep] Spectrum plot saved → {spectrum_dir}/spectrum_{_tag(R)}.png")


# ---------------------------------------------------------------------------
# Aggregation: collect every completed resolution (only tags from THIS
# sweep's RESOLUTION_SWEEP list — stale files from an older sweep definition
# are ignored) into the resolution-convergence plot + grid, mirroring the
# fdtd_sweep/hybrid_modes_sweep examples' aggregation pattern.
# ---------------------------------------------------------------------------
def aggregate() -> None:
    done, T_through_rows, T_cross_rows = [], [], []
    for R in RESOLUTION_SWEEP:
        path = os.path.join(data_dir, f"{_tag(R)}.npz")
        if not os.path.exists(path):
            continue
        d = np.load(path)
        done.append(R)
        T_through_rows.append(d["T_through"])
        T_cross_rows.append(d["T_cross"])

    n_done, n_all = len(done), len(RESOLUTION_SWEEP)
    if n_done < 2:
        print(f"[sweep] Aggregation skipped — only {n_done}/{n_all} resolutions done.")
        return

    done = np.asarray(done)
    T_through_grid = np.vstack(T_through_rows)   # [n_done, n_wl]
    T_cross_grid = np.vstack(T_cross_rows)       # [n_done, n_wl]
    T_total_grid = T_through_grid + T_cross_grid

    def _write_npz(tmp):
        with open(tmp, "wb") as f:
            np.savez(f, resolutions=done, wavelengths=wavelengths,
                     T_through=T_through_grid, T_cross=T_cross_grid)

    _atomic_replace_write(os.path.join(output_dir, "resolution_sweep.npz"), _write_npz)
    print(f"[sweep] Resolution-sweep grid saved ({n_done}/{n_all} resolutions) → {output_dir}/resolution_sweep.npz")

    # constrained_layout (not tight_layout) is what correctly reserves space
    # for the colorbar and suptitle together — see the hybrid modes-sweep
    # script's own note on this same layout bug.
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)

    # -- Panel 1: spectrum shape vs resolution -----------------------------
    ax = axes[0]
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=RESOLUTION_SWEEP.min(), vmax=RESOLUTION_SWEEP.max())
    for i, R in enumerate(done):
        ax.plot(wavelengths, T_cross_grid[i], color=cmap(norm(R)), lw=1.2)
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Cross transmission")
    ax.set_title(f"Spectrum shape vs resolution (coupling_length={COUPLING_LENGTH:.1f} µm)")
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Resolution (px/µm)")

    # -- Panel 2: convergence vs resolution --------------------------------
    ax = axes[1]
    for grid, label, color in [
        (T_through_grid, "Through", "steelblue"),
        (T_cross_grid,   "Cross",   "darkorange"),
        (T_total_grid,   "Total",   "k"),
    ]:
        mean, lo, hi = grid.mean(axis=1), grid.min(axis=1), grid.max(axis=1)
        ax.plot(done, mean, label=label, color=color, marker="o", markersize=4,
                lw=1.5 if label == "Total" else 1.2,
                linestyle="--" if label == "Cross" else "-")
        ax.fill_between(done, lo, hi, color=color, alpha=0.15, lw=0)
    ax.set_xlabel("Resolution (px/µm)")
    ax.set_ylabel("Transmission (mean over wavelength, shaded = min–max)")
    ax.set_title("Convergence vs mesh resolution")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(RESOLUTION_SWEEP)

    if n_done < n_all:
        fig.text(0.02, 0.97, f"{n_done}/{n_all} resolutions done", ha="left", va="top",
                  fontsize=8, color="#555555",
                  bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"))

    fig.suptitle("Directional Coupler — Single-Domain FDTD Resolution Convergence Study", fontsize=11)
    _atomic_replace_write(
        os.path.join(output_dir, "resolution_convergence.png"),
        lambda tmp: fig.savefig(tmp, dpi=300, format="png"),
    )
    plt.close("all")
    print(f"[sweep] Convergence plot saved ({n_done}/{n_all} resolutions) → {output_dir}/resolution_convergence.png")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        i = int(sys.argv[1])
        if not 0 <= i < len(RESOLUTION_SWEEP):
            raise SystemExit(
                f"index {i} out of range — sweep has {len(RESOLUTION_SWEEP)} "
                f"resolutions (0..{len(RESOLUTION_SWEEP) - 1})")
        simulate_resolution(RESOLUTION_SWEEP[i])
    else:
        print(f"[sweep] No index given — running all {len(RESOLUTION_SWEEP)} "
              "resolutions serially (the SLURM array batch script is much faster).")
        for R in RESOLUTION_SWEEP:
            simulate_resolution(R)
    aggregate()
