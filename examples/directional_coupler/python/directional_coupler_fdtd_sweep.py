"""
Tutorial: directional coupler — brute-force FDTD coupling-length sweep.

The same 220 nm SOI, C-band device as directional_coupler_hybrid.py /
directional_coupler_fdtd.py, characterized the expensive way: one
single-domain, full-device FDTD simulation (see directional_coupler_fdtd.py
for that method's anatomy) **per coupling length**, sweeping

    coupling_length = 0, 2, 4, ..., 60 µm   (31 lengths)

Contrast with the hybrid tutorial, which sweeps 201 lengths for the cost of
2 FDTD + 1 EME solve (lengths enter only as analytic EME phase factors, see
its module docstring): here every length is its own complete first-principles
simulation, with nothing decomposed and no truncated mode basis. That makes
this sweep the brute-force validation reference for the hybrid's
coupling_length_sweep.png curve — at roughly 31/3 ≈ 10× the solver cost of
the entire hybrid sweep, for 15% of its length resolution. The aggregate
plot overlays the hybrid's curves (when that example's length_sweep.npz
exists) so the two methods can be compared point by point.

Outputs (example_output/directional_coupler/directional_coupler_fdtd_sweep/):

    spectrums/spectrum_L{length}um.png   per-length transmission/IL spectrum
                                         (same layout as the single-length
                                         example's fdtd_results.png)
    data/L{length}um.npz                 per-length raw sweep data
    coupling_length_sweep.png            aggregate power-transfer-vs-length
                                         plot (same layout as the hybrid's)
    length_sweep.npz                     aggregate grid over completed lengths

Run:
    # one length per SLURM array task (see run_directional_coupler_fdtd_sweep.batch):
    python /path/to/examples/directional_coupler/directional_coupler_fdtd_sweep.py <index 0..30>

    # or the whole sweep serially (slow — batch queue strongly preferred):
    python /path/to/examples/directional_coupler/directional_coupler_fdtd_sweep.py

Every invocation ends by (re)aggregating whatever lengths have finished so
far, writing the aggregate plot/data atomically — so concurrent array tasks
can't corrupt them, the plot is usable (annotated "n/31 lengths") while the
array is still draining, and the last task to finish leaves the complete
sweep behind. Re-running an index simply overwrites that length's files.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from feces import DirectionalCoupler, FDTDSolver, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "directional_coupler", "directional_coupler_fdtd_sweep"))
spectrum_dir = os.path.join(output_dir, "spectrums")
data_dir = os.path.join(output_dir, "data")
for _d in (output_dir, spectrum_dir, data_dir):
    os.makedirs(_d, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "example_sim_cache"))

# Hybrid example's saved length sweep — overlaid on the aggregate plot when present.
_HYBRID_SWEEP_NPZ = os.path.abspath(os.path.join(
    output_dir, "..", "directional_coupler_hybrid", "length_sweep.npz"))

# ---------------------------------------------------------------------------
# Global simulation parameters — identical to directional_coupler_fdtd.py,
# except coupling_length is swept instead of pinned.
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm)
WL_MIN = 1.50       # wavelength sweep start (µm)
WL_MAX = 1.60       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

COUPLING_LENGTHS = np.arange(0.0, 60.0 + 1e-9, 2.0)   # µm — 0:2:60, 31 lengths

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)


def _tag(L: float) -> str:
    """Stable per-length file tag, e.g. L04.0um / L60.0um."""
    return f"L{L:04.1f}um"


def _atomic_replace_write(path: str, write_fn) -> None:
    """write_fn(tmp_path) then atomically rename onto path — concurrent array
    tasks all regenerate the aggregate files, and a plain overwrite could be
    read (or collide) half-written; os.replace is atomic on POSIX."""
    tmp = f"{path}.tmp.{os.getpid()}"
    write_fn(tmp)
    os.replace(tmp, path)


def make_coupler(L: float) -> DirectionalCoupler:
    return DirectionalCoupler(
        waveguide_width=0.45,       # single-mode strip width (µm)
        gap=0.25,                   # coupling-region edge-to-edge gap (µm)
        coupling_length=float(L),   # straight parallel interaction length (µm) — the swept knob
        bend_length=9.0,            # horizontal extent of each S-bend (µm)
        waveguide_offset=3.0,       # far-end centre-to-centre port pitch (µm)
        bend_type="euler",
        waveguide_height=0.22,      # Si layer thickness (µm) — SOI_220
    )


def make_solver() -> FDTDSolver:
    return FDTDSolver(
        res=RESOLUTION,
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
# Per-length simulation: one broadband full-device FDTD run -> spectrum plot
# in spectrums/ + raw data in data/ for the aggregation step.
# ---------------------------------------------------------------------------
def simulate_length(L: float) -> None:
    dc = make_coupler(L)
    solver = make_solver()
    print(f"[sweep] coupling_length = {L:.1f} µm "
          f"(device total {dc.total_length:.1f} µm, single FDTD run, "
          f"{N_WL} wavelengths {WL_MIN:.2f}–{WL_MAX:.2f} µm)")

    smat = solver.solve(dc, wavelengths, port="bot", input_mode=1, save_dir=save_dir)
    smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength
    T_through = smat.transmission(port="out_0", in_mode=1, out_mode=1)
    T_cross = smat.transmission(port="out_1", in_mode=1, out_mode=1)
    T_total = T_through + T_cross

    np.savez(
        os.path.join(data_dir, f"{_tag(L)}.npz"),
        coupling_length=L, wavelengths=wavelengths,
        T_through=T_through, T_cross=T_cross,
    )

    idx = len(wavelengths) // 2
    print(f"[sweep] L={L:.1f} µm @ λ={wavelengths[idx]:.3f} µm: "
          f"through={T_through[idx]:.4f} cross={T_cross[idx]:.4f} "
          f"total={T_total[idx]:.4f} "
          f"ratio={T_cross[idx]/max(T_total[idx], 1e-20):.3f}")

    # -- Per-length spectrum plot (same layout as fdtd_results.png) --------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax = axes[0]
    ax.plot(wavelengths, T_through, label="Through", color="steelblue")
    ax.plot(wavelengths, T_cross,   label="Cross",   color="darkorange", linestyle="--")
    ax.plot(wavelengths, T_total,   label="Total",   color="k", lw=1.5)
    ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Transmission")
    ax.set_title(f"Transmission vs Wavelength (coupling_length={L:.1f} µm)")
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
    ax.set_title(f"Insertion Loss vs Wavelength (coupling_length={L:.1f} µm)")
    ax.legend(fontsize=8)
    ax.grid(True)

    plt.suptitle("Directional Coupler — Single-Domain FDTD (length-sweep point)", fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(spectrum_dir, f"spectrum_{_tag(L)}.png"), dpi=300)
    plt.close("all")
    print(f"[sweep] Spectrum plot saved → {spectrum_dir}/spectrum_{_tag(L)}.png")


# ---------------------------------------------------------------------------
# Aggregation: collect every completed length (only tags from THIS sweep's
# COUPLING_LENGTHS list — stale files from an older sweep definition are
# ignored) into the coupling-length sweep plot + grid, mirroring the hybrid
# example's coupling_length_sweep.png / length_sweep.npz.
# ---------------------------------------------------------------------------
def aggregate() -> None:
    done, T_through_rows, T_cross_rows = [], [], []
    for L in COUPLING_LENGTHS:
        path = os.path.join(data_dir, f"{_tag(L)}.npz")
        if not os.path.exists(path):
            continue
        d = np.load(path)
        done.append(L)
        T_through_rows.append(d["T_through"])
        T_cross_rows.append(d["T_cross"])

    n_done, n_all = len(done), len(COUPLING_LENGTHS)
    if n_done < 2:
        print(f"[sweep] Aggregation skipped — only {n_done}/{n_all} lengths done.")
        return

    done = np.asarray(done)
    T_through_grid = np.vstack(T_through_rows)   # [n_done, n_wl]
    T_cross_grid = np.vstack(T_cross_rows)       # [n_done, n_wl]
    T_total_grid = T_through_grid + T_cross_grid

    def _write_npz(tmp):
        # open file object, not a path: np.savez would append ".npz" to the
        # tmp path and break the atomic-rename bookkeeping
        with open(tmp, "wb") as f:
            np.savez(f, coupling_lengths=done, wavelengths=wavelengths,
                     T_through=T_through_grid, T_cross=T_cross_grid)

    _atomic_replace_write(os.path.join(output_dir, "length_sweep.npz"), _write_npz)
    print(f"[sweep] Length-sweep grid saved ({n_done}/{n_all} lengths) → {output_dir}/length_sweep.npz")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for grid, label, color in [
        (T_through_grid, "Through", "steelblue"),
        (T_cross_grid,   "Cross",   "darkorange"),
        (T_total_grid,   "Total",   "k"),
    ]:
        mean, lo, hi = grid.mean(axis=1), grid.min(axis=1), grid.max(axis=1)
        ax.plot(done, mean, label=f"{label} (FDTD)", color=color,
                marker="o", markersize=3,
                lw=1.5 if label == "Total" else 1.2,
                linestyle="--" if label == "Cross" else "-")
        ax.fill_between(done, lo, hi, color=color, alpha=0.15, lw=0)

    # Overlay the hybrid sweep's mean curves for direct method comparison.
    if os.path.exists(_HYBRID_SWEEP_NPZ):
        h = np.load(_HYBRID_SWEEP_NPZ)
        for grid, color in [(h["T_through"], "steelblue"), (h["T_cross"], "darkorange")]:
            ax.plot(h["coupling_lengths"], grid.mean(axis=1), color=color,
                    lw=0.9, linestyle=":", alpha=0.8)
        ax.plot([], [], color="gray", lw=0.9, linestyle=":", label="Hybrid (mean)")

    if n_done < n_all:
        ax.text(0.02, 0.97, f"{n_done}/{n_all} lengths done", transform=ax.transAxes,
                ha="left", va="top", fontsize=8, color="#555555",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"))
    ax.set_xlabel("Coupling length (µm)")
    ax.set_ylabel("Transmission (mean over wavelength, shaded = min–max)")
    ax.set_title("Directional Coupler — Power Transfer vs Coupling Length (FDTD)")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    _atomic_replace_write(
        os.path.join(output_dir, "coupling_length_sweep.png"),
        lambda tmp: fig.savefig(tmp, dpi=300, format="png"),
    )
    plt.close("all")
    print(f"[sweep] Sweep plot saved ({n_done}/{n_all} lengths) → {output_dir}/coupling_length_sweep.png")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        i = int(sys.argv[1])
        if not 0 <= i < len(COUPLING_LENGTHS):
            raise SystemExit(
                f"index {i} out of range — sweep has {len(COUPLING_LENGTHS)} "
                f"lengths (0..{len(COUPLING_LENGTHS) - 1})")
        simulate_length(COUPLING_LENGTHS[i])
    else:
        print(f"[sweep] No index given — running all {len(COUPLING_LENGTHS)} "
              "lengths serially (the SLURM array batch script is much faster).")
        for L in COUPLING_LENGTHS:
            simulate_length(L)
    aggregate()
