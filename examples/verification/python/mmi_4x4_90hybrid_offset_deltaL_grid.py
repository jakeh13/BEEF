"""4x4 MMI 90-degree hybrid -- joint (waveguide_offset x ΔL) grid, from
ALREADY-CACHED data, zero new FDTD/EME compute.

validation_markdown/mmi_4x4_90_hybrid_validation.md, "Open decisions" item
5's own follow-up. The single-ΔL offset screening sweep (job 11817429)
found offset alone (at the fixed ΔL=-2.75 tuned specifically for
offset=1.5) is dramatically worse at every other offset tested -- but that
sweep only ever evaluated the ONE ΔL slice each input matrix was combined
with. Two things are already sitting on disk, needing no new solves:

- The 4 real input-section matrices per offset (job 11817429's own 2 real
  + 2 sign-derived per new offset; job 11778128's 4 real for the
  offset=1.5 baseline).
- The FULL 81-point ΔL body length-sweep grid (job 11768847), which is
  offset-independent (validation doc: "waveguide_offset body-cache
  independence", confirmed via cache_only=True) -- the SAME grid applies
  to every offset, not just offset=1.5.

The original scripts only ever sliced that grid at ONE ΔL per run. This
script re-combines the SAME cached input matrices against EVERY point on
that grid instead of just one -- pure numpy post-processing, no meep/meow
call, no cluster time -- to answer whether ANY (offset, ΔL) combination in
the already-explored ranges clears the paper's cited system requirements
(CMRR >= 20 dBe, |Δφ_j| <= 5 deg), not just whether offset alone does at
the one ΔL already tried.

Run (login node, seconds):
    python mmi_4x4_90hybrid_offset_deltaL_grid.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import SMatrix

import mmi_4x4_90hybrid_offset_sweep as base

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "offset_deltaL_grid"))
os.makedirs(output_dir, exist_ok=True)

OFFSETS = [1.2, 1.35, 1.5, 1.65, 1.8]
DRIVEN = [0, 2]
OUTPUTS = [0, 1, 2, 3]
PORT_NAMES = ["out5", "out6", "out7", "out8"]
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}

_DELTA_L_GRID = np.linspace(-10.0, 10.0, 81)


def get_M_in(input_solver, mmi, boundaries, offset, sign):
    """All 4 input matrices for `offset`, entirely from cache -- real
    solves via cache_only=True (job 11778128 for the baseline, job
    11817429's own input_idx=0,1 for every other offset), ports 2/3
    derived via the (already-validated) sign vector for anything not the
    baseline, where all 4 are real."""
    if abs(offset - base.BASELINE_OFFSET) < 1e-9:
        return [base.solve_input(input_solver, mmi, boundaries, k, cache_only=True) for k in range(4)]
    M0 = base.solve_input(input_solver, mmi, boundaries, 0, cache_only=True)
    M1 = base.solve_input(input_solver, mmi, boundaries, 1, cache_only=True)
    M3 = M0 * sign[None, :, None]
    M2 = M1 * sign[None, :, None]
    return [M0, M1, M2, M3]


def metrics_at(M_in, M_body_one, wavelengths):
    """Same metric set as offset_sweep.py's solve_offset_point(), for one
    (offset, ΔL) point -- M_body_one shape (1, n_wl, N, N)."""
    M_out = [Mk.transpose(0, 2, 1) for Mk in M_in]
    coeffs = {}
    for i in DRIVEN:
        tmp = np.einsum("Lwab,wbk->Lwak", M_body_one, M_in[i])
        for j in OUTPUTS:
            c = np.einsum("wab,Lwbk->Lwak", M_out[j], tmp)
            coeffs[(i, j)] = c[0, :, 0, 0]

    smat = {i: SMatrix(wavelengths=wavelengths) for i in DRIVEN}
    for i in DRIVEN:
        for j in OUTPUTS:
            smat[i].set(coeffs[(i, j)], port=PORT_NAMES[j])

    imbalance = {i: smat[i].imbalance_db(PORT_NAMES) for i in DRIVEN}
    excess_loss = {i: smat[i].excess_loss_db(PORT_NAMES) for i in DRIVEN}
    cmrr_I = {i: smat[i].cmrr_db("out5", "out8") for i in DRIVEN}
    cmrr_Q = {i: smat[i].cmrr_db("out6", "out7") for i in DRIVEN}

    phi = {j: smat[2].cross_excitation_phase_diff_deg(smat[0], port=PORT_NAMES[j]) for j in OUTPUTS}
    dphi = {}
    for j in [1, 2, 3]:
        raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
        dphi[j] = ((raw + 180.0) % 360.0) - 180.0

    max_dphi = float(np.max(np.abs(np.concatenate([dphi[1], dphi[2], dphi[3]]))))
    min_cmrr = float(np.min(np.concatenate([cmrr_I[0], cmrr_Q[0], cmrr_I[2], cmrr_Q[2]])))
    max_excess_loss = float(np.max(np.concatenate([excess_loss[0], excess_loss[2]])))
    max_imbalance = float(np.max(np.concatenate([imbalance[0], imbalance[2]])))
    return max_dphi, min_cmrr, max_excess_loss, max_imbalance


def main():
    input_solver, body_solver = base.make_solvers()

    max_dphi_grid = np.full((len(OFFSETS), 81), np.nan)
    min_cmrr_grid = np.full((len(OFFSETS), 81), np.nan)
    max_xl_grid = np.full((len(OFFSETS), 81), np.nan)
    max_imb_grid = np.full((len(OFFSETS), 81), np.nan)

    sign, _ = base.derive_sign_vector(input_solver, body_solver)

    for oi, offset in enumerate(OFFSETS):
        mmi = base.make_mmi(offset)
        boundaries, body_stub = base.boundaries_for(mmi)
        M_in = get_M_in(input_solver, mmi, boundaries, offset, sign)

        _body_lengths_grid = (0.75 * base.L_PI + _DELTA_L_GRID) - 2 * body_stub
        M_body_grid = body_solver.solve_section_length_sweep(
            mmi, base.wavelengths, boundaries[1], boundaries[2], _body_lengths_grid,
            n_modes_in=base.N_MODES, n_modes_out=base.N_MODES, save_dir=base.save_dir, cache_only=True,
        )
        print(f"offset={offset}: body grid fetched from cache, shape={M_body_grid.shape}")

        for li in range(81):
            max_dphi, min_cmrr, max_xl, max_imb = metrics_at(M_in, M_body_grid[li:li + 1], base.wavelengths)
            max_dphi_grid[oi, li] = max_dphi
            min_cmrr_grid[oi, li] = min_cmrr
            max_xl_grid[oi, li] = max_xl
            max_imb_grid[oi, li] = max_imb
        best_li = int(np.argmin(max_dphi_grid[oi]))
        print(f"  best ΔL={_DELTA_L_GRID[best_li]:+.2f} (of the 81-pt grid): "
              f"max|Δφ|={max_dphi_grid[oi, best_li]:.3f} deg  min_CMRR={min_cmrr_grid[oi, best_li]:.3f} dBe")

    np.savez(os.path.join(output_dir, "grid_results.npz"),
             offsets=np.array(OFFSETS), delta_L_grid=_DELTA_L_GRID,
             max_dphi=max_dphi_grid, min_cmrr=min_cmrr_grid,
             max_excess_loss=max_xl_grid, max_imbalance=max_imb_grid)

    # Global best point (by max|Δφ| -- the metric currently failing spec
    # hardest / most load-bearing for the "coherent receiver" use case).
    flat_idx = np.nanargmin(max_dphi_grid)
    oi_best, li_best = np.unravel_index(flat_idx, max_dphi_grid.shape)
    print(f"\n=== Global best (by max|Δφ|) over the explored grid ===")
    print(f"  offset={OFFSETS[oi_best]}  ΔL={_DELTA_L_GRID[li_best]:+.2f}  "
          f"max|Δφ|={max_dphi_grid[oi_best, li_best]:.3f} deg  "
          f"min_CMRR={min_cmrr_grid[oi_best, li_best]:.3f} dBe  "
          f"max_XL={max_xl_grid[oi_best, li_best]:.3f} dB  "
          f"max_imb={max_imb_grid[oi_best, li_best]:.3f} dB")
    print(f"  (for reference, baseline offset=1.5 @ ΔL=-2.75: max|Δφ|=12.320, min_CMRR=1.531)")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    titles = ["max|Δφ| (deg, lower better, spec ≤5)", "min CMRR (dBe, higher better, spec ≥20)",
              "max excess loss (dB, lower better)", "max imbalance (dB, lower better)"]
    grids = [max_dphi_grid, min_cmrr_grid, max_xl_grid, max_imb_grid]
    for ax, grid, title in zip(axes.flat, grids, titles):
        for oi, offset in enumerate(OFFSETS):
            ax.plot(_DELTA_L_GRID, grid[oi], label=f"offset={offset}")
        ax.set_xlabel("ΔL (µm)")
        ax.set_title(title)
        ax.legend(fontsize=7)
        ax.grid(True)
    plt.suptitle("4×4 90° hybrid — joint (offset, ΔL) grid, from cached data only", fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "offset_deltaL_grid.png"), dpi=200)
    print(f"\nplot saved -> {output_dir}/offset_deltaL_grid.png")


if __name__ == "__main__":
    main()
