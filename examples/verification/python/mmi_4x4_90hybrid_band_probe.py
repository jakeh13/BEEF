"""Per-wavelength metric breakdown for an already-simulated geometry.

Answers two questions the campaign's aggregated metrics cannot, at ZERO new
simulation (MMI4X4_CACHE_ONLY is forced on, so a cache miss raises):

  1. HOW MUCH OF THE WORST CASE IS BAND EDGE. wg_eval reports max|dphi|, min
     CMRR, max excess loss and max imbalance aggregated over the band, so a
     number set by 1.29/1.33 um is indistinguishable from one that is flat.
     That distinction decides whether narrowing the optimisation band to
     1310 +- 10 nm buys real margin or just re-labels the same device.

  2. WHETHER EXCESS LOSS IS BASIS-CONVERGED. The N_MODES=20 confirmation
     (job 12752148) predates excess-loss recording, so phase and CMRR are
     known converged at W=8.0 and loss is not. Run this at N_MODES 16 and 20
     on the same geometry to find out.

Usage (env, same knobs as wg_eval):
    MMI4X4_W=8.5 MMI4X4_YGAP=0.825 MMI4X4_TAPER_LENGTH=4.0 \
        python mmi_4x4_90hybrid_band_probe.py

Reads the geometry's own best mm_length out of its eval npz so the probe scores
exactly the length wg_eval selected, rather than re-running the length search.
"""
import os

os.environ["MMI4X4_CACHE_ONLY"] = "1"

import numpy as np

import mmi_4x4_90hybrid_wg_eval as W
from beef import SMatrix


def npz_name():
    fn = f"eval_W{W.MM_WIDTH:.2f}_y{W.Y_GAP:.2f}"
    if abs(W.TAPER_LENGTH - 1.0) > 1e-9:
        fn += f"_tl{W.TAPER_LENGTH:.2f}"
    if W.N_MODES != 16:
        fn += f"_n{W.N_MODES}"
    return fn + ".npz"


def main():
    path = os.path.join(W.out_dir, npz_name())
    with np.load(path) as d:
        best_len = float(d["best_mm_length"])
        ref = {k: float(d[k]) for k in
               ("max_dphi", "min_cmrr", "objective") if k in d}
    print(f"PROBE: {npz_name()}  W={W.MM_WIDTH} y_gap={W.Y_GAP} "
          f"taper_width={W.TAPER_WIDTH:.3f} taper_length={W.TAPER_LENGTH} "
          f"N_MODES={W.N_MODES}", flush=True)
    print(f"PROBE: scoring at that npz's own best mm_length={best_len:.2f} um; "
          f"aggregated reference from the npz: {ref}", flush=True)

    input_solver, body_solver = W.make_solvers()
    mmi = W.make_mmi(W.L_CENTRE)
    boundaries, body_stub = W.boundaries_for(mmi)

    M_in = [input_solver.solve_section(
        mmi, W.wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=W.N_MODES,
        stub_left=input_solver.waveguide_length, stub_right=0.0,
        input_idx=k, save_dir=W.save_dir, cache_only=True) for k in range(4)]

    li = int(np.argmin(np.abs(W.MM_LENGTH_GRID - best_len)))
    body_lengths = W.MM_LENGTH_GRID - 2 * body_stub
    M_body = body_solver.solve_section_length_sweep(
        mmi, W.wavelengths, boundaries[1], boundaries[2], body_lengths,
        n_modes_in=W.N_MODES, n_modes_out=W.N_MODES, save_dir=W.save_dir,
        cache_only=True)[li:li + 1]

    # Same cascade as wg_eval.metrics_at, but keeping the wavelength axis.
    M_out = [Mk.transpose(0, 2, 1) for Mk in M_in]
    smat = {i: SMatrix(wavelengths=W.wavelengths) for i in W.DRIVEN}
    for i in W.DRIVEN:
        tmp = np.einsum("Lwab,wbk->Lwak", M_body, M_in[i])
        for j in W.OUTPUTS:
            c = np.einsum("wab,Lwbk->Lwak", M_out[j], tmp)
            smat[i].set(c[0, :, 0, 0], port=W.PORT_NAMES[j])

    cm = np.array([smat[i].cmrr_db(a, b) for i in W.DRIVEN
                   for a, b in (("out5", "out8"), ("out6", "out7"))])
    xl = np.array([smat[i].excess_loss_db(W.PORT_NAMES) for i in W.DRIVEN])
    ib = np.array([smat[i].imbalance_db(W.PORT_NAMES) for i in W.DRIVEN])
    phi = {j: np.asarray(smat[2].cross_excitation_phase_diff_deg(
        smat[0], port=W.PORT_NAMES[j])) for j in W.OUTPUTS}
    dphi = []
    for j in (1, 2, 3):
        raw = phi[j] - phi[0] - (W.TARGET_DEG[j] - W.TARGET_DEG[0])
        dphi.append(((raw + 180.0) % 360.0) - 180.0)
    dphi = np.abs(np.array(dphi))

    print(f"\n{'wl (um)':>8} {'max|dphi|':>10} {'min CMRR':>10} "
          f"{'max XL':>8} {'max imb':>9}")
    print("-" * 50)
    for w in range(len(W.wavelengths)):
        print(f"{W.wavelengths[w]:8.3f} {dphi[:, w].max():10.3f} "
              f"{cm[:, w].min():10.3f} {xl[:, w].max():8.3f} "
              f"{ib[:, w].max():9.3f}")
    print("-" * 50)
    print(f"{'BAND':>8} {dphi.max():10.3f} {cm.min():10.3f} "
          f"{xl.max():8.3f} {ib.max():9.3f}   <- what wg_eval reports")
    c = len(W.wavelengths) // 2
    print(f"{'CENTRE':>8} {dphi[:, c].max():10.3f} {cm[:, c].min():10.3f} "
          f"{xl[:, c].max():8.3f} {ib[:, c].max():9.3f}   <- "
          f"{W.wavelengths[c]:.3f} um alone (band-edge penalty = BAND - CENTRE)")


if __name__ == "__main__":
    main()
