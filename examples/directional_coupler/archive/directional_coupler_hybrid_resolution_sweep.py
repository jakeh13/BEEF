"""
Tutorial: directional coupler — hybrid FDTD/EME mesh-resolution convergence sweep.

The same 220 nm SOI, C-band device and method as directional_coupler_hybrid.py
(2 FDTD S-bend sections + 1 EME body, N_MODES=5), at a single FIXED
coupling_length = 15.13 µm — the same 3 dB length used by
directional_coupler_fdtd_resolution_sweep.py (and directional_coupler_fdtd.py)
— so the two resolution-convergence studies are directly, apples-to-apples
comparable (same physical device length). Instead of a single fixed mesh
resolution, this sweeps

    RESOLUTION = 20, 30, 40, ..., 120 px/µm   (11 values, 20:10:120)

— the FDTD/EME mesh density used by both sub-solvers — to check mesh
convergence for the hybrid method: does the result actually stop changing as
the mesh gets finer, and does it match the single-domain ground truth at the
same resolution/length?

Unlike the N_MODES sweep (directional_coupler_hybrid_modes_sweep.py), where
only the EME body solver needed rebuilding per sweep value (the shared FDTD
solver's resolution never changed), here BOTH sub-solvers depend on the swept
value, so both FDTDSolver and EMESolver are rebuilt fresh per resolution.

Outputs (example_output/directional_coupler/directional_coupler_hybrid_resolution_sweep/):

    spectrums/spectrum_R{res}.png   per-resolution transmission/IL spectrum
                                     (same layout as directional_coupler_
                                     fdtd_resolution_sweep.py's spectrum plot)
    data/R{res}.npz                 per-resolution raw sweep data
    resolution_convergence.png      aggregate two-panel convergence plot
    resolution_sweep.npz            aggregate grid over completed resolutions

Run:
    # one resolution per SLURM array task (see
    # batch/sweep/run_directional_coupler_hybrid_resolution_sweep.batch):
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

from beef import DirectionalCoupler, EMESolver, FDTDSolver, SOI_220, atomic_replace_write

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "directional_coupler", "directional_coupler_hybrid_resolution_sweep"))
spectrum_dir = os.path.join(output_dir, "spectrums")
data_dir = os.path.join(output_dir, "data")
for _d in (output_dir, spectrum_dir, data_dir):
    os.makedirs(_d, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters — identical to directional_coupler_hybrid.py,
# except RESOLUTION is swept instead of pinned, and coupling_length is fixed
# to the same value directional_coupler_fdtd_resolution_sweep.py uses (rather
# than hybrid.py's own 200-point length sweep) so the two resolution studies
# are directly comparable.
# ---------------------------------------------------------------------------
WL_MIN = 1.50       # wavelength sweep start (µm)
WL_MAX = 1.60       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

N_MODES = 5   # even + odd supermodes carried at the two body interfaces — hybrid.py's default, pinned here.
              # Was 6: band 5 is a radiation mode, dropped at no cost (resolved/directional_coupler_validation.md check 2).

# Same 3 dB length as directional_coupler_fdtd_resolution_sweep.py — see
# directional_coupler_fdtd.py for the provenance (full-device FDTD ground truth
# crosses 0.5 at 15.129 µm). Was 13.5 µm until 2026-07-28, a stale pre-axis-fix
# value that is really a 43:57 split; directional_coupler_validation.md item 14.
#
# As with the FDTD sweep, this did not invalidate the convergence conclusions —
# one device across meshes — only the claim that 13.5 was a 3 dB length.
COUPLING_LENGTH = 15.13   # µm

# µm of the straight coupling region folded into each FDTD section, past the
# S-bend/coupling-region discontinuity — same value/reasoning as hybrid.py's
# own BODY_STUB (see its step-3 comments). Must satisfy
# 2*BODY_STUB <= COUPLING_LENGTH.
BODY_STUB = 0.25   # µm

RESOLUTION_SWEEP = np.arange(20.0, 120.0 + 1e-9, 10.0)   # px/µm — 20:10:120, 11 values — the swept knob

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)


def _tag(R: float) -> str:
    """Stable per-resolution file tag, e.g. R020 / R120."""
    return f"R{int(R):03d}"


# ---------------------------------------------------------------------------
# 1. Geometry — shared across every resolution value (device dimensions don't
#    depend on mesh resolution, only device layout).
# ---------------------------------------------------------------------------
dc = DirectionalCoupler(
    waveguide_width=0.45,             # single-mode strip width (µm)
    gap=0.25,                         # coupling-region edge-to-edge gap (µm)
    coupling_length=COUPLING_LENGTH,  # straight parallel interaction length (µm)
    bend_length=9.0,                  # horizontal extent of each S-bend (µm)
    waveguide_offset=3.0,             # far-end centre-to-centre port pitch (µm)
    bend_type="euler",                # "euler" (cubic-Bezier S-bend) or "circular" (matched arcs)
    waveguide_height=0.22,             # Si layer thickness (µm) — SOI_220
)

half_body = dc.coupling_length / 2
half_total = dc.total_length / 2
boundaries = [-half_total, -half_body + BODY_STUB, +half_body - BODY_STUB, +half_total]

# EME body length: each FDTD section already folds BODY_STUB µm of the
# straight coupling region into itself, so the three sections' straight
# coverage adds back up to exactly COUPLING_LENGTH.
BODY_LENGTH = COUPLING_LENGTH - 2 * BODY_STUB
assert BODY_LENGTH >= 0, (
    f"COUPLING_LENGTH ({COUPLING_LENGTH:.2f} µm) must be at least 2*BODY_STUB "
    f"({2 * BODY_STUB:.2f} µm)"
)

# Interface-monitor pin for the input-section FDTD runs — see
# directional_coupler_hybrid.py's own comment for the full reasoning: spans
# the full non-PML domain centered on the coupled pair's shared y=0 axis, so
# it matches the EME body's own y=0-centered supermode basis. Fixed
# monitor_size (not resolution-dependent) so every resolution's monitor
# geometry stays comparable.
_SIM_WIDTH = 6.0
_iface_mon = dict(monitor_center=0.0, monitor_size=_SIM_WIDTH)

print(f"Directional coupler total length: {dc.total_length:.2f} µm "
      f"(coupling length {dc.coupling_length:.2f} µm, "
      f"coupling pitch {dc.coupling_pitch:.3f} µm, far pitch {dc.waveguide_offset:.2f} µm)")
print(f"Coupling-length budget: 2 x {BODY_STUB:.2f} µm inside the FDTD sections "
      f"+ EME body {BODY_LENGTH:.2f} µm = {BODY_LENGTH + 2 * BODY_STUB:.2f} µm total "
      f"(nominal {COUPLING_LENGTH:.2f} µm; interface monitors at x = ±{abs(boundaries[1]):.2f} µm)")
print(f"Resolution sweep: {RESOLUTION_SWEEP[0]:.0f}..{RESOLUTION_SWEEP[-1]:.0f} px/µm "
      f"({len(RESOLUTION_SWEEP)} values), N_MODES={N_MODES}, {N_WL} wavelengths")


# ---------------------------------------------------------------------------
# Per-resolution simulation: 2 FDTD section solves + 1 EME body solve at the
# fixed COUPLING_LENGTH, combined via the same matrix-chain contraction as
# hybrid.py's step 5. Unlike the N_MODES sweep, BOTH sub-solvers depend on
# the swept value, so both are rebuilt fresh here (not shared at module scope).
# ---------------------------------------------------------------------------
def simulate_resolution(R: float) -> None:
    R = float(R)
    print(f"[sweep] resolution = {R:.0f} px/µm "
          f"(coupling_length={COUPLING_LENGTH:.1f} µm, {N_WL} wavelengths)")

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
    # Both engines are meshed at the SAME R here, unlike the other hybrid
    # scripts, which now pin FDTD 40 / EME 80 separately (dc_validation.md item
    # 11). That is intentional for a resolution sweep — but it means this sweep
    # measures the two convergences *convolved* and cannot separate them. It is
    # why a 7% EME error went unnoticed: the joint curve looked converged
    # because FDTD dominates the observable (T_cross) while EME dominates Δβ.
    # Read a flat curve here as "the pair is jointly converged", never as "each
    # engine is converged". For the per-engine answer see
    # examples/verification/python/verify_item_6_engine_agreement.py.
    body_solver = EMESolver(
        num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholder length; solve_section_length_sweep() below sets its own
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
        dc, wavelengths, boundaries[1], boundaries[2], np.array([BODY_LENGTH]),
        n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
    )[0]  # [n_wl, N_MODES, N_MODES] — single-element length axis dropped

    tmp = np.einsum("wij,wjk->wik", M_body, M_in0)             # [n_wl, N_MODES, 1]
    full_through = np.einsum("wij,wjk->wik", M_out_through, tmp)  # [n_wl, 1, 1]
    full_cross = np.einsum("wij,wjk->wik", M_out_cross, tmp)      # [n_wl, 1, 1]

    T_through = np.abs(full_through[:, 0, 0]) ** 2   # [n_wl]
    T_cross = np.abs(full_cross[:, 0, 0]) ** 2        # [n_wl]
    T_total = T_through + T_cross

    # Atomic — aggregate() reads every task's file, including one a concurrent
    # array task is mid-write. See directional_coupler_hybrid_stub_sweep.py,
    # where a plain np.savez here killed 7 of 9 tasks (dc_validation.md item 6).
    def _write_task_npz(tmp: str) -> None:
        # `with`, not a bare open(): zipfile does NOT close a file object
        # it was handed, so relying on refcounting to flush would let
        # os.replace move a truncated file — the very bug being fixed.
        with open(tmp, "wb") as fh:
            np.savez(
                fh,
                resolution=R, wavelengths=wavelengths,
                T_through=T_through, T_cross=T_cross,
            )

    atomic_replace_write(
        os.path.join(data_dir, f"{_tag(R)}.npz"), _write_task_npz)

    idx = len(wavelengths) // 2
    print(f"[sweep] res={R:.0f} px/µm @ λ={wavelengths[idx]:.3f} µm: "
          f"through={T_through[idx]:.4f} cross={T_cross[idx]:.4f} "
          f"total={T_total[idx]:.4f} "
          f"ratio={T_cross[idx]/max(T_total[idx], 1e-20):.3f}")

    # -- Per-resolution spectrum plot (same layout as directional_coupler_
    #    fdtd_resolution_sweep.py's spectrum_{tag}.png) --------------------
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

    plt.suptitle(f"Directional Coupler — Hybrid FDTD/EME (coupling_length={COUPLING_LENGTH:.1f} µm, resolution-sweep point)", fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(spectrum_dir, f"spectrum_{_tag(R)}.png"), dpi=300)
    plt.close("all")
    print(f"[sweep] Spectrum plot saved → {spectrum_dir}/spectrum_{_tag(R)}.png")


# ---------------------------------------------------------------------------
# Aggregation: collect every completed resolution (only tags from THIS
# sweep's RESOLUTION_SWEEP list — stale files from an older sweep definition
# are ignored) into the resolution-convergence plot + grid, mirroring
# directional_coupler_fdtd_resolution_sweep.py's aggregation exactly.
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

    atomic_replace_write(os.path.join(output_dir, "resolution_sweep.npz"), _write_npz)
    print(f"[sweep] Resolution-sweep grid saved ({n_done}/{n_all} resolutions) → {output_dir}/resolution_sweep.npz")

    # constrained_layout (not tight_layout) is what correctly reserves space
    # for the colorbar and suptitle together — see the N_MODES sweep script's
    # own note on this same layout bug.
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

    fig.suptitle("Directional Coupler — Hybrid FDTD/EME Resolution Convergence Study", fontsize=11)
    atomic_replace_write(
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
