"""
Tutorial: directional coupler — hybrid FDTD/EME mesh-resolution convergence sweep.

The same 220 nm SOI, C-band device as directional_coupler_hybrid.py, at its
default N_MODES=6, but instead of a single fixed mesh resolution this sweeps

    RESOLUTION = 20, 30, 40, ..., 120 px/µm   (11 values, 20:10:120)

— the FDTD/EME mesh density used by both sub-solvers — while still sweeping
the full 201-point coupling_length range at each resolution (still "free" per
resolution via EMESolver.solve_section_length_sweep(), since mode profiles
don't depend on section length — see its docstring in beef/solvers/eme.py).
This checks mesh convergence for the hybrid method: does the result actually
stop changing as the mesh gets finer, and does it match
directional_coupler_hybrid.py's own fixed-resolution=40 result?

Unlike the N_MODES sweep (directional_coupler_hybrid_modes_sweep.py), where
only the EME body solver needed rebuilding per sweep value (the shared FDTD
solver's resolution never changed), here BOTH sub-solvers depend on the swept
value, so both FDTDSolver and EMESolver are rebuilt fresh per resolution.

Outputs (example_output/directional_coupler/directional_coupler_hybrid_resolution_sweep/):

    plots/coupling_length_sweep_R{res}.png   per-resolution power-transfer-vs-length
                                              plot (same layout as hybrid.py's
                                              coupling_length_sweep.png)
    data/R{res}.npz                          per-resolution raw length-sweep grid
                                              (same schema as hybrid.py's own
                                              length_sweep.npz)
    resolution_convergence_sweep.png         aggregate two-panel convergence plot
    resolution_sweep.npz                     aggregate grid over completed resolutions

Run:
    # one resolution per SLURM array task (see
    # run_directional_coupler_hybrid_resolution_sweep.batch):
    python /path/to/examples/directional_coupler/python/directional_coupler_hybrid_resolution_sweep.py <index 0..10>

    # or the whole sweep serially (slow — batch queue strongly preferred):
    python /path/to/examples/directional_coupler/python/directional_coupler_hybrid_resolution_sweep.py

Every invocation ends by (re)aggregating whatever resolutions have finished so
far, writing the aggregate plot/data atomically — so concurrent array tasks
can't corrupt them, the plot is usable (annotated "n/11 resolutions done")
while the array is still draining, and the last task to finish leaves the
complete sweep behind. Re-running an index simply overwrites that
resolution's files.

Note on cost: 3-D FDTD cost scales roughly as resolution^4 (cell count grows
~resolution^3, and the timestep shrinks ~1/resolution), so the res=120 task
is dramatically more expensive than the res=20/40 tasks — see the batch
script's own comment on walltime for this reason.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import DirectionalCoupler, EMESolver, FDTDSolver, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "directional_coupler", "directional_coupler_hybrid_resolution_sweep"))
plot_dir = os.path.join(output_dir, "plots")
data_dir = os.path.join(output_dir, "data")
for _d in (output_dir, plot_dir, data_dir):
    os.makedirs(_d, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_sim_cache"))

# hybrid.py's own fixed-resolution=40 length sweep — overlaid on panel 1 of
# the aggregate plot when present, as a reference against the repo's existing
# tutorial example.
_HYBRID_SWEEP_NPZ = os.path.abspath(os.path.join(
    output_dir, "..", "directional_coupler_hybrid", "length_sweep.npz"))

# ---------------------------------------------------------------------------
# Global simulation parameters — identical to directional_coupler_hybrid.py,
# except RESOLUTION is swept instead of pinned.
# ---------------------------------------------------------------------------
WL_MIN = 1.50       # wavelength sweep start (µm)
WL_MAX = 1.60       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

N_MODES = 6   # even + odd supermodes carried at the two body interfaces — hybrid.py's default, pinned here

coupling_lengths = np.linspace(0.0, 60.0, 201)   # µm — same 201-pt range as hybrid.py
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

RESOLUTION_SWEEP = np.arange(20.0, 120.0 + 1e-9, 10.0)   # px/µm — 20:10:120, 11 values — the swept knob


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


# ---------------------------------------------------------------------------
# 1. Geometry — shared across every resolution value (device dimensions don't
#    depend on mesh resolution, only device layout).
# ---------------------------------------------------------------------------
L_mid = coupling_lengths[len(coupling_lengths) // 2]

dc = DirectionalCoupler(
    waveguide_width=0.45,    # single-mode strip width (µm)
    gap=0.25,                # coupling-region edge-to-edge gap (µm)
    coupling_length=L_mid,   # straight parallel interaction length (µm) — sweep reference value
    bend_length=9.0,         # horizontal extent of each S-bend (µm)
    waveguide_offset=3.0,    # far-end centre-to-centre port pitch (µm)
    bend_type="euler",       # "euler" (cubic-Bezier S-bend) or "circular" (matched arcs)
    waveguide_height=0.22,   # Si layer thickness (µm) — SOI_220
)

half_body = dc.coupling_length / 2
half_total = dc.total_length / 2
boundaries = [-half_total, -half_body, +half_body, +half_total]

# Interface-monitor override for the input-section FDTD runs — see
# directional_coupler_hybrid.py's own comment for the full reasoning: spans
# the full non-PML domain centered on the coupled pair's shared y=0 axis, so
# it matches the EME body's own y=0-centered supermode basis. Fixed
# monitor_size (not resolution-dependent) so every resolution's monitor
# geometry stays comparable.
_SIM_WIDTH = 6.0
_iface_mon = dict(monitor_center=0.0, monitor_size=_SIM_WIDTH)

print(f"Directional coupler total length: {dc.total_length:.2f} µm "
      f"(coupling pitch {dc.coupling_pitch:.3f} µm, far pitch {dc.waveguide_offset:.2f} µm)")
print(f"Resolution sweep: {RESOLUTION_SWEEP[0]:.0f}..{RESOLUTION_SWEEP[-1]:.0f} px/µm "
      f"({len(RESOLUTION_SWEEP)} values), each over {len(coupling_lengths)} coupling "
      f"lengths x {N_WL} wavelengths, N_MODES={N_MODES}")


# ---------------------------------------------------------------------------
# Per-resolution simulation: 2 FDTD section solves + 1 EME length-sweep solve,
# combined via the same matrix-chain contraction as hybrid.py's step 5.
# Unlike the N_MODES sweep, BOTH sub-solvers depend on the swept value, so
# both are rebuilt fresh here (not shared at module scope).
# ---------------------------------------------------------------------------
def simulate_resolution(R: float) -> None:
    R = float(R)
    print(f"[sweep] resolution = {R:.0f} px/µm ({N_WL} wavelengths, "
          f"{len(coupling_lengths)} coupling lengths)")

    fdtd_solver = FDTDSolver(
        res=R,
        sim_width=_SIM_WIDTH,
        sim_height=1.5,
        dpml=0.5,
        source_margin=0.0,
        monitor_margin=0.25,
        three_d=True,   # full 3-D FDTD: real 220 nm Si core in oxide, PML on all six faces
        half_domain_aperture=True,
        material_stack=SOI_220,
        decay_tol=1e-11,
    )
    body_solver = EMESolver(
        num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholder length; solve_section_length_sweep() sweeps its own
        oxide_width=6.0, oxide_height=1.5, res=R, material_stack=SOI_220,
        mode_solver="meep",
    )

    M_in0 = fdtd_solver.solve_section(
        dc, wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
        input_idx=0, save_dir=save_dir, **_iface_mon,
    )  # [n_wl, N_MODES, 1] — "in" arm -> body modes
    M_in1 = fdtd_solver.solve_section(
        dc, wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
        input_idx=1, save_dir=save_dir, **_iface_mon,
    )  # [n_wl, N_MODES, 1] — "isolated" arm -> body modes

    M_out_through = M_in0.transpose(0, 2, 1)  # mirror-symmetry reuse, not a new simulation
    M_out_cross = M_in1.transpose(0, 2, 1)

    M_body = body_solver.solve_section_length_sweep(
        dc, wavelengths, boundaries[1], boundaries[2], coupling_lengths,
        n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
    )  # [n_L, n_wl, N_MODES, N_MODES]

    tmp = np.einsum("Lwij,wjk->Lwik", M_body, M_in0)               # [n_L, n_wl, N_MODES, 1]
    full_through = np.einsum("wij,Lwjk->Lwik", M_out_through, tmp) # [n_L, n_wl, 1, 1]
    full_cross = np.einsum("wij,Lwjk->Lwik", M_out_cross, tmp)     # [n_L, n_wl, 1, 1]

    T_through_grid = np.abs(full_through[:, :, 0, 0]) ** 2   # [n_L, n_wl]
    T_cross_grid = np.abs(full_cross[:, :, 0, 0]) ** 2        # [n_L, n_wl]
    T_total_grid = T_through_grid + T_cross_grid

    np.savez(
        os.path.join(data_dir, f"{_tag(R)}.npz"),
        resolution=R, coupling_lengths=coupling_lengths, wavelengths=wavelengths,
        T_through=T_through_grid, T_cross=T_cross_grid,
    )

    ratio_mean = T_cross_grid.mean(axis=1) / np.clip(T_total_grid.mean(axis=1), 1e-20, None)
    idx_ref = int(np.argmin(np.abs(ratio_mean - 0.5)))
    print(f"[sweep] resolution={R:.0f} px/µm: 50:50 split near L={coupling_lengths[idx_ref]:.1f} µm "
          f"(ratio={ratio_mean[idx_ref]:.3f})")

    # -- Per-resolution sweep plot (same layout as hybrid.py's coupling_length_sweep.png) --
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for grid, label, color in [
        (T_through_grid, "Through", "steelblue"),
        (T_cross_grid,   "Cross",   "darkorange"),
        (T_total_grid,   "Total",   "k"),
    ]:
        mean, lo, hi = grid.mean(axis=1), grid.min(axis=1), grid.max(axis=1)
        ax.plot(coupling_lengths, mean, label=label, color=color,
                lw=1.5 if label == "Total" else 1.2,
                linestyle="--" if label == "Cross" else "-")
        ax.fill_between(coupling_lengths, lo, hi, color=color, alpha=0.15, lw=0)
    ax.set_xlabel("Coupling length (µm)")
    ax.set_ylabel("Transmission (mean over wavelength, shaded = min–max)")
    ax.set_title(f"Directional Coupler — Power Transfer vs Coupling Length (resolution={R:.0f} px/µm)")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(plot_dir, f"coupling_length_sweep_{_tag(R)}.png"), dpi=300)
    plt.close("all")
    print(f"[sweep] Plot saved → {plot_dir}/coupling_length_sweep_{_tag(R)}.png")


# ---------------------------------------------------------------------------
# Aggregation: collect every completed resolution into a two-panel
# convergence plot + stacked grid, mirroring the N_MODES sweep's aggregation
# pattern.
# ---------------------------------------------------------------------------
def aggregate() -> None:
    done, T_through_rows, T_cross_rows = [], [], []
    for R in RESOLUTION_SWEEP:
        path = os.path.join(data_dir, f"{_tag(R)}.npz")
        if not os.path.exists(path):
            continue
        d = np.load(path)
        done.append(float(R))
        T_through_rows.append(d["T_through"])
        T_cross_rows.append(d["T_cross"])

    n_done, n_all = len(done), len(RESOLUTION_SWEEP)
    if n_done < 2:
        print(f"[sweep] Aggregation skipped — only {n_done}/{n_all} resolutions done.")
        return

    done = np.asarray(done)
    T_through_stack = np.stack(T_through_rows)   # [n_R, n_L, n_wl]
    T_cross_stack = np.stack(T_cross_rows)       # [n_R, n_L, n_wl]
    T_total_stack = T_through_stack + T_cross_stack

    def _write_npz(tmp):
        with open(tmp, "wb") as f:
            np.savez(f, resolutions=done, coupling_lengths=coupling_lengths, wavelengths=wavelengths,
                     T_through=T_through_stack, T_cross=T_cross_stack)

    _atomic_replace_write(os.path.join(output_dir, "resolution_sweep.npz"), _write_npz)
    print(f"[sweep] Resolution-sweep grid saved ({n_done}/{n_all} resolutions) → {output_dir}/resolution_sweep.npz")

    # Reference coupling length for the convergence panel: the 50:50 cross/
    # through split length, computed from whichever completed resolution is
    # LARGEST so far — the closest available proxy to a converged result,
    # self-correcting as more resolutions land, robust to out-of-order SLURM
    # array completion.
    i_ref_R = int(np.argmax(done))
    ratio_ref = T_cross_stack[i_ref_R].mean(axis=1) / np.clip(T_total_stack[i_ref_R].mean(axis=1), 1e-20, None)
    idx_L_ref = int(np.argmin(np.abs(ratio_ref - 0.5)))
    L_ref = coupling_lengths[idx_L_ref]
    R_ref = done[i_ref_R]

    # constrained_layout (not tight_layout/subplots_adjust) is what correctly
    # reserves space for the colorbar and suptitle together — tight_layout
    # doesn't know about the colorbar's stolen axis space and mis-clips
    # panel 2's long y-label (see the N_MODES sweep's own note on this bug).
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)

    # -- Panel 1: S-curve shape vs resolution ------------------------------
    ax = axes[0]
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=RESOLUTION_SWEEP.min(), vmax=RESOLUTION_SWEEP.max())
    for i, R in enumerate(done):
        mean_cross = T_cross_stack[i].mean(axis=1)
        ax.plot(coupling_lengths, mean_cross, color=cmap(norm(R)), lw=1.2)

    if os.path.exists(_HYBRID_SWEEP_NPZ):
        h = np.load(_HYBRID_SWEEP_NPZ)
        ax.plot(h["coupling_lengths"], h["T_cross"].mean(axis=1), color="gray",
                lw=1.2, linestyle=":", label="Hybrid tutorial (resolution=40, mean)")

    ax.axvline(L_ref, color="gray", linestyle="--", lw=0.8,
               label=f"Reference L ({L_ref:.1f} µm, from R={R_ref:.0f})")
    ax.set_xlabel("Coupling length (µm)")
    ax.set_ylabel("Cross transmission (mean over wavelength)")
    ax.set_title("S-curve shape vs resolution")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Resolution (px/µm)")

    # -- Panel 2: convergence vs resolution at the fixed reference length --
    ax = axes[1]
    T_through_ref = T_through_stack[:, idx_L_ref, :]   # [n_R, n_wl]
    T_cross_ref = T_cross_stack[:, idx_L_ref, :]
    T_total_ref = T_total_stack[:, idx_L_ref, :]
    for grid, label, color in [
        (T_through_ref, "Through", "steelblue"),
        (T_cross_ref,   "Cross",   "darkorange"),
        (T_total_ref,   "Total",   "k"),
    ]:
        mean, lo, hi = grid.mean(axis=1), grid.min(axis=1), grid.max(axis=1)
        ax.plot(done, mean, label=label, color=color, marker="o", markersize=4,
                lw=1.5 if label == "Total" else 1.2,
                linestyle="--" if label == "Cross" else "-")
        ax.fill_between(done, lo, hi, color=color, alpha=0.15, lw=0)
    ax.set_xlabel("Resolution (px/µm)")
    ax.set_ylabel("Transmission (mean over wavelength, shaded = min–max)")
    ax.set_title(f"Convergence vs resolution (at coupling_length={L_ref:.1f} µm)")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(RESOLUTION_SWEEP)

    if n_done < n_all:
        fig.text(0.02, 0.97, f"{n_done}/{n_all} resolutions done", ha="left", va="top",
                  fontsize=8, color="#555555",
                  bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"))

    fig.suptitle("Directional Coupler — Hybrid FDTD/EME Resolution Convergence Study", fontsize=11)
    _atomic_replace_write(
        os.path.join(output_dir, "resolution_convergence_sweep.png"),
        lambda tmp: fig.savefig(tmp, dpi=300, format="png"),
    )
    plt.close("all")
    print(f"[sweep] Convergence plot saved ({n_done}/{n_all} resolutions) → {output_dir}/resolution_convergence_sweep.png")


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
