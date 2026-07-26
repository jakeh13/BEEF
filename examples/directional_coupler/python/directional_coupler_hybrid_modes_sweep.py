"""
Tutorial: directional coupler — hybrid FDTD/EME N_MODES convergence sweep.

The same 220 nm SOI, C-band device as directional_coupler_hybrid.py, but
instead of sweeping coupling_length at a fixed N_MODES (=6, hybrid.py's
default), this sweeps

    N_MODES = 2, 3, 4, ..., 12   (11 values)

— the number of EME supermodes carried across the two coupling-region
interfaces — while still sweeping the full 200-point coupling_length range
at each N_MODES value. hybrid.py's own module docstring claims each arm is
single-mode by design, so 2 modes (even + odd) should fully capture the
coupling physics, and its N_MODES=6 default is set generously above that
minimum ("bump N_MODES up only if the arm width supports higher-order
per-arm modes too"). This sweep tests that claim directly: does the hybrid
method's result actually change as more modes are carried, or does it
converge by N_MODES=2 as expected?

Unlike the coupling_length axis (which solve_section_length_sweep() sweeps
for free — see its docstring in beef/solvers/eme.py — because mode profiles
don't depend on section length), N_MODES directly changes the eigenmode
solve itself, so each N_MODES value requires its own fresh 2 FDTD section
solves + 1 EME mode-solve-per-wavelength (the same per-task cost as one full
run of hybrid.py's step 5, just at a different N_MODES). The FDTD run count
doesn't grow with N_MODES (still 1 excited-mode run per section, by
reciprocity — see FDTDSolver.solve_section()'s docstring); only the
eigenmode-decomposition analysis and the EME mode count scale with N.

Outputs (example_output/directional_coupler/directional_coupler_hybrid_modes_sweep/):

    plots/coupling_length_sweep_N{NN}.png   per-N_MODES power-transfer-vs-length
                                             plot (same layout as hybrid.py's
                                             coupling_length_sweep.png)
    data/N{NN}.npz                          per-N_MODES raw length-sweep grid
                                             (same schema as hybrid.py's own
                                             length_sweep.npz)
    modes_convergence_sweep.png             aggregate two-panel convergence plot
    modes_sweep.npz                         aggregate grid over completed N_MODES

Run:
    # one N_MODES value per SLURM array task (see
    # batch/sweep/run_directional_coupler_hybrid_modes_sweep.batch):
    python /path/to/examples/directional_coupler/python/directional_coupler_hybrid_modes_sweep.py <index 0..10>

    # or the whole sweep serially (slow — batch queue strongly preferred):
    python /path/to/examples/directional_coupler/python/directional_coupler_hybrid_modes_sweep.py

Every invocation ends by (re)aggregating whatever N_MODES values have
finished so far, writing the aggregate plot/data atomically — so concurrent
array tasks can't corrupt them, the plot is usable (annotated "n/11 N_MODES
done") while the array is still draining, and the last task to finish leaves
the complete sweep behind. Re-running an index simply overwrites that
N_MODES value's files.
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
    _HERE, "..", "..", "..", "..", "example_output", "directional_coupler", "directional_coupler_hybrid_modes_sweep"))
plot_dir = os.path.join(output_dir, "plots")
data_dir = os.path.join(output_dir, "data")
for _d in (output_dir, plot_dir, data_dir):
    os.makedirs(_d, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_sim_cache"))

# hybrid.py's own fixed-N_MODES=6 length sweep — overlaid on panel 1 of the
# aggregate plot when present, as a reference against the repo's existing
# tutorial example.
_HYBRID_SWEEP_NPZ = os.path.abspath(os.path.join(
    output_dir, "..", "directional_coupler_hybrid", "length_sweep.npz"))

# ---------------------------------------------------------------------------
# Global simulation parameters — identical to directional_coupler_hybrid.py.
# ---------------------------------------------------------------------------
# The two engines get SEPARATE resolutions, and this is deliberate — a single
# shared RESOLUTION was how a 7% error hid for months (dc_validation.md item 11).
#
# meep's FDTD is converged at 40: its implied Δβ moves only 2.74% across res
# 30–90, non-monotonically (noise, not a trend). meow/MPB on its own Mesh2D is
# NOT — Δβ_EME at 1.55 µm runs 0.1125 (res 20) → 0.0957 (40) → 0.0926 (60) →
# 0.0894 (80) → 0.0895 (100), so 40 overstates the even/odd splitting by 7% and
# 80 is the first converged point. meep's subpixel smoothing is the likely
# reason it converges so much earlier on the same nominal grid.
#
# Measured directly in dc_validation.md Done item 6 (job 11489846) and pinned by
# tests/test_eme_convergence.py. Do not collapse these back into one knob.
FDTD_RESOLUTION = 40.0   # pixels/µm — converged; see above
EME_RESOLUTION = 80.0    # pixels/µm — 40 is NOT converged; see above
WL_MIN = 1.50       # wavelength sweep start (µm)
WL_MAX = 1.60       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

# µm of the straight coupling region folded into each FDTD section, past the
# S-bend/coupling-region discontinuity — same value/reasoning as hybrid.py's
# own BODY_STUB (see its step-3 comments). The sweep starts at
# 2*BODY_STUB because a device shorter than both stubs put together can't be
# split into these three sections at all.
BODY_STUB = 0.25   # µm

coupling_lengths = np.linspace(2 * BODY_STUB, 60.0, 200)   # µm — same 200-pt range as hybrid.py
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

NUM_MODES_SWEEP = np.arange(2, 13)   # 2..12 inclusive, 11 values — the swept knob


def _tag(N: int) -> str:
    """Stable per-N_MODES file tag, e.g. N02 / N12."""
    return f"N{int(N):02d}"


def _atomic_replace_write(path: str, write_fn) -> None:
    """write_fn(tmp_path) then atomically rename onto path — concurrent array
    tasks all regenerate the aggregate files, and a plain overwrite could be
    read (or collide) half-written; os.replace is atomic on POSIX."""
    tmp = f"{path}.tmp.{os.getpid()}"
    write_fn(tmp)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# 1. Geometry and FDTD solver — shared across every N_MODES value (neither
#    depends on the coupling-region mode count, only on device geometry).
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

_fdtd_common = dict(
    res=FDTD_RESOLUTION,
    sim_width=6.0,
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    half_domain_aperture=True,
    material_stack=SOI_220,
    decay_tol=1e-11,
)
fdtd_solver = FDTDSolver(**_fdtd_common)

half_body = dc.coupling_length / 2
half_total = dc.total_length / 2
boundaries = [-half_total, -half_body + BODY_STUB, +half_body - BODY_STUB, +half_total]

# Interface-monitor pin for the input-section FDTD runs — see
# directional_coupler_hybrid.py's own comment for the full reasoning: spans
# the full non-PML domain centered on the coupled pair's shared y=0 axis, so
# it matches the EME body's own y=0-centered supermode basis.
_iface_mon = dict(monitor_center=0.0, monitor_size=fdtd_solver.sim_width)

# EME body length per swept coupling_length: each FDTD section already folds
# BODY_STUB µm of the straight coupling region into itself, so the three
# sections' straight coverage adds back up to exactly coupling_length.
body_lengths = coupling_lengths - 2 * BODY_STUB
# >= 0, not > 0: the shortest swept length is exactly 2*BODY_STUB — the
# degenerate device whose whole coupling region lives in the two FDTD sections
# and whose EME body is zero-length (an identity propagation matrix).
assert body_lengths[0] >= 0, (
    f"coupling_lengths must start at >= 2*BODY_STUB ({2 * BODY_STUB:.2f} µm); "
    f"got {coupling_lengths[0]:.2f} µm"
)
assert np.allclose(body_lengths + 2 * BODY_STUB, coupling_lengths)

print(f"Directional coupler total length: {dc.total_length:.2f} µm "
      f"(coupling pitch {dc.coupling_pitch:.3f} µm, far pitch {dc.waveguide_offset:.2f} µm)")
print(f"Coupling-length budget: 2 x {BODY_STUB:.2f} µm inside the FDTD sections "
      f"+ EME body {body_lengths[0]:.2f}–{body_lengths[-1]:.2f} µm "
      f"= {coupling_lengths[0]:.2f}–{coupling_lengths[-1]:.2f} µm total "
      f"(interface monitors at x = ±{abs(boundaries[1]):.2f} µm)")
print(f"N_MODES sweep: {NUM_MODES_SWEEP[0]}..{NUM_MODES_SWEEP[-1]} "
      f"({len(NUM_MODES_SWEEP)} values), each over {len(coupling_lengths)} coupling "
      f"lengths x {N_WL} wavelengths")


# ---------------------------------------------------------------------------
# Per-N_MODES simulation: 2 FDTD section solves + 1 EME length-sweep solve,
# combined via the same matrix-chain contraction as hybrid.py's step 5.
# ---------------------------------------------------------------------------
def simulate_n_modes(N: int) -> None:
    N = int(N)
    print(f"[sweep] N_MODES = {N} ({N_WL} wavelengths, {len(coupling_lengths)} coupling lengths)")

    body_solver = EMESolver(
        num_modes=[N], num_cells=[1], lengths=[1.0],  # placeholder length; solve_section_length_sweep() sweeps its own
        oxide_width=6.0, oxide_height=1.5, res=EME_RESOLUTION, material_stack=SOI_220,
        mode_solver="meep",
    )

    M_in0 = fdtd_solver.solve_section(
        dc, wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N,
        stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
        input_idx=0, save_dir=save_dir, **_iface_mon,
    )  # [n_wl, N, 1] — "in" arm -> body modes
    M_in1 = fdtd_solver.solve_section(
        dc, wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N,
        stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
        input_idx=1, save_dir=save_dir, **_iface_mon,
    )  # [n_wl, N, 1] — "isolated" arm -> body modes

    M_out_through = M_in0.transpose(0, 2, 1)  # mirror-symmetry reuse, not a new simulation
    M_out_cross = M_in1.transpose(0, 2, 1)

    M_body = body_solver.solve_section_length_sweep(
        dc, wavelengths, boundaries[1], boundaries[2], body_lengths,
        n_modes_in=N, n_modes_out=N, save_dir=save_dir,
    )  # [n_L, n_wl, N, N] — plotted/saved against coupling_lengths (the physical total)

    tmp = np.einsum("Lwij,wjk->Lwik", M_body, M_in0)               # [n_L, n_wl, N, 1]
    full_through = np.einsum("wij,Lwjk->Lwik", M_out_through, tmp) # [n_L, n_wl, 1, 1]
    full_cross = np.einsum("wij,Lwjk->Lwik", M_out_cross, tmp)     # [n_L, n_wl, 1, 1]

    T_through_grid = np.abs(full_through[:, :, 0, 0]) ** 2   # [n_L, n_wl]
    T_cross_grid = np.abs(full_cross[:, :, 0, 0]) ** 2        # [n_L, n_wl]
    T_total_grid = T_through_grid + T_cross_grid

    np.savez(
        os.path.join(data_dir, f"{_tag(N)}.npz"),
        N_modes=N, coupling_lengths=coupling_lengths, wavelengths=wavelengths,
        T_through=T_through_grid, T_cross=T_cross_grid,
    )

    ratio_mean = T_cross_grid.mean(axis=1) / np.clip(T_total_grid.mean(axis=1), 1e-20, None)
    idx_ref = int(np.argmin(np.abs(ratio_mean - 0.5)))
    print(f"[sweep] N_MODES={N}: 50:50 split near L={coupling_lengths[idx_ref]:.1f} µm "
          f"(ratio={ratio_mean[idx_ref]:.3f})")

    # -- Per-N_MODES sweep plot (same layout as hybrid.py's coupling_length_sweep.png) --
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
    ax.set_title(f"Directional Coupler — Power Transfer vs Coupling Length (N_MODES={N})")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(plot_dir, f"coupling_length_sweep_{_tag(N)}.png"), dpi=300)
    plt.close("all")
    print(f"[sweep] Plot saved → {plot_dir}/coupling_length_sweep_{_tag(N)}.png")


# ---------------------------------------------------------------------------
# Aggregation: collect every completed N_MODES value into a two-panel
# convergence plot + stacked grid, mirroring the fdtd_sweep example's
# atomic-write aggregation pattern.
# ---------------------------------------------------------------------------
def aggregate() -> None:
    done, T_through_rows, T_cross_rows = [], [], []
    for N in NUM_MODES_SWEEP:
        path = os.path.join(data_dir, f"{_tag(N)}.npz")
        if not os.path.exists(path):
            continue
        d = np.load(path)
        done.append(int(N))
        T_through_rows.append(d["T_through"])
        T_cross_rows.append(d["T_cross"])

    n_done, n_all = len(done), len(NUM_MODES_SWEEP)
    if n_done < 2:
        print(f"[sweep] Aggregation skipped — only {n_done}/{n_all} N_MODES done.")
        return

    done = np.asarray(done)
    T_through_stack = np.stack(T_through_rows)   # [n_N, n_L, n_wl]
    T_cross_stack = np.stack(T_cross_rows)       # [n_N, n_L, n_wl]
    T_total_stack = T_through_stack + T_cross_stack

    def _write_npz(tmp):
        with open(tmp, "wb") as f:
            np.savez(f, N_modes=done, coupling_lengths=coupling_lengths, wavelengths=wavelengths,
                     T_through=T_through_stack, T_cross=T_cross_stack)

    _atomic_replace_write(os.path.join(output_dir, "modes_sweep.npz"), _write_npz)
    print(f"[sweep] Modes-sweep grid saved ({n_done}/{n_all} N_MODES) → {output_dir}/modes_sweep.npz")

    # Reference coupling length for the convergence panel: the 50:50 cross/
    # through split length, computed from whichever completed N_MODES is
    # LARGEST so far — the closest available proxy to a converged result,
    # self-correcting as more N_MODES values land, robust to out-of-order
    # SLURM array completion.
    i_ref_N = int(np.argmax(done))
    ratio_ref = T_cross_stack[i_ref_N].mean(axis=1) / np.clip(T_total_stack[i_ref_N].mean(axis=1), 1e-20, None)
    idx_L_ref = int(np.argmin(np.abs(ratio_ref - 0.5)))
    L_ref = coupling_lengths[idx_L_ref]
    N_ref = done[i_ref_N]

    # constrained_layout (not tight_layout/subplots_adjust) is what correctly
    # reserves space for the colorbar and suptitle together — tight_layout
    # doesn't know about the colorbar's stolen axis space and mis-clips
    # panel 2's long y-label.
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)

    # -- Panel 1: S-curve shape vs N_MODES ---------------------------------
    ax = axes[0]
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=NUM_MODES_SWEEP.min(), vmax=NUM_MODES_SWEEP.max())
    for i, N in enumerate(done):
        mean_cross = T_cross_stack[i].mean(axis=1)
        ax.plot(coupling_lengths, mean_cross, color=cmap(norm(N)), lw=1.2)

    if os.path.exists(_HYBRID_SWEEP_NPZ):
        h = np.load(_HYBRID_SWEEP_NPZ)
        ax.plot(h["coupling_lengths"], h["T_cross"].mean(axis=1), color="gray",
                lw=1.2, linestyle=":", label="Hybrid tutorial (N_MODES=6, mean)")

    ax.axvline(L_ref, color="gray", linestyle="--", lw=0.8,
               label=f"Reference L ({L_ref:.1f} µm, from N={N_ref})")
    ax.set_xlabel("Coupling length (µm)")
    ax.set_ylabel("Cross transmission (mean over wavelength)")
    ax.set_title("S-curve shape vs N_MODES")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="N_MODES")

    # -- Panel 2: convergence vs N_MODES at the fixed reference length -----
    ax = axes[1]
    T_through_ref = T_through_stack[:, idx_L_ref, :]   # [n_N, n_wl]
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
    ax.set_xlabel("N_MODES")
    ax.set_ylabel("Transmission (mean over wavelength, shaded = min–max)")
    ax.set_title(f"Convergence vs N_MODES (at coupling_length={L_ref:.1f} µm)")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(NUM_MODES_SWEEP)

    if n_done < n_all:
        fig.text(0.02, 0.97, f"{n_done}/{n_all} N_MODES done", ha="left", va="top",
                  fontsize=8, color="#555555",
                  bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"))

    fig.suptitle("Directional Coupler — Hybrid N_MODES Convergence Study", fontsize=11)
    _atomic_replace_write(
        os.path.join(output_dir, "modes_convergence_sweep.png"),
        lambda tmp: fig.savefig(tmp, dpi=300, format="png"),
    )
    plt.close("all")
    print(f"[sweep] Convergence plot saved ({n_done}/{n_all} N_MODES) → {output_dir}/modes_convergence_sweep.png")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        i = int(sys.argv[1])
        if not 0 <= i < len(NUM_MODES_SWEEP):
            raise SystemExit(
                f"index {i} out of range — sweep has {len(NUM_MODES_SWEEP)} "
                f"N_MODES values (0..{len(NUM_MODES_SWEEP) - 1})")
        simulate_n_modes(NUM_MODES_SWEEP[i])
    else:
        print(f"[sweep] No index given — running all {len(NUM_MODES_SWEEP)} "
              "N_MODES values serially (the SLURM array batch script is much faster).")
        for N in NUM_MODES_SWEEP:
            simulate_n_modes(N)
    aggregate()
