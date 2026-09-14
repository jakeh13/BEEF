"""4x4 MMI 90-degree hybrid -- re-score the ENTIRE mm_length axis under the
CALIBRATED body propagation constant, using only simulations already on disk.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-09-11.

WHY THIS IS FREE, AND WHY IT IS THE FIRST THING TO DO. The n_eff calibration
is a property of the BODY CROSS-SECTION, not of the body LENGTH -- it is a
per-(wavelength, mode) propagation constant, and the length sweep just
propagates those betas further. So every quantity needed to re-score any
mm_length under calibration is already cached:

  * M_in for all four input ports at the production geometry (the FDTD
    input-section solves; length-INdependent -- the taper and stub geometry
    does not change when the body gets longer), and
  * the calibrated body length-sweep grid, 5.00-400.00 um @ 0.05 um, solved
    once by job 13098454 with neff_override=neff_calibrated.

That makes re-optimising mm_length pure cascade arithmetic over 7601 grid
points -- seconds of numpy, zero new FDTD or EME solves. Every OTHER knob
(mm_width, taper_width, taper_length, waveguide_width) changes the body
CROSS-SECTION and therefore invalidates the calibration, so each of those
needs a fresh 16-mode calibration run before it can be scored honestly.
Length is the one axis that comes for free, so it is the one to exhaust first.

WHY IT MATTERS. The rounds 8-10 campaign picked mm_length = 159.30 um by
scoring ~45 points through the UNCALIBRATED path, i.e. with a body index
0.0458 too low. A wrong beta does not merely shift the metrics at a fixed
length -- it shifts WHERE ALONG THE LENGTH AXIS the self-imaging condition
lands, because the accumulated phase is beta*L. To first order the optimum
moves to keep beta*L invariant, L_opt_cal ~ L_opt_uncal * n_uncal/n_cal
= 159.30 * 2.9632/3.0094 ~ 156.85 um -- a ~2.4 um shift, roughly 50 grid
steps. This script measures the real shift rather than trusting that
estimate, and re-scores both paths on the same grid so the comparison is
like-for-like.

Scoring replicates mmi_4x4_90hybrid_hybrid.py exactly (same SMatrix methods,
same A11 phase convention, same Eq. 15 targets) so the numbers are directly
comparable to every CMRR figure this campaign has quoted.

Run:
    MMI4X4_RESCORE_STRIDE=200 python mmi_4x4_90hybrid_calibrated_length_rescore.py
"""
import os
import time

import numpy as np

from beef import MMI, EMESolver, FDTDSolver, SMatrix, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
save_dir = os.path.join(_ROOT, "example_sim_cache")
OUT_DIR = os.path.join(_ROOT, "example_output", "verification",
                       "mmi_4x4_90hybrid_calibrated_length_rescore")
os.makedirs(OUT_DIR, exist_ok=True)
CAL_PATH = os.path.join(_ROOT, "example_output", "mmi_4x4_90hybrid",
                        "calib2_per_mode", "aggregated_beta_n16.npz")

# --- design point: pinned to mmi_4x4_90hybrid_hybrid.py -------------------
WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 8.5
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 1.30
TAPER_LENGTH = 4.0
MM_LENGTH = 159.30          # the campaign's UNCALIBRATED optimum -- the incumbent
N_MODES = 16
N_MODES_BODY = 16
N_MODES_TAPER = 12
BODY_STUB = 0.25
WL_MIN, WL_MAX, N_WL = 1.29, 1.33, 3
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

DRIVEN = [0, 2]                                  # paper's I = {1, 3}
OUTPUTS = [0, 1, 2, 3]
PORT_NAMES = ["out5", "out6", "out7", "out8"]
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}     # Eq. 15
SPEC_PHASE_DEG = 5.0
SPEC_CMRR_DBE = 20.0

STRIDE = int(os.environ.get("MMI4X4_RESCORE_STRIDE", 1))
# Common grid: the uncalibrated sweep starts at 20.0, the calibrated at 5.00.
# Scoring both on the SAME points is the whole point of the comparison.
GRID = np.round(np.arange(20.0, 400.0 + 1e-9, 0.05), 2)


def make_mmi(mm_length):
    return MMI(
        n_inputs=4, n_outputs=4,
        mm_length=mm_length, mm_width=MM_WIDTH, waveguide_width=WAVEGUIDE_WIDTH,
        taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
        waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
        use_neighbor_aperture=True,
    )


def main():
    cal = np.load(CAL_PATH)
    neff_calibrated = cal["neff_calibrated"]
    print(f"[rescore] calibration: neff_cal[1.31, mode1] = {neff_calibrated[1, 0]:.6f}")

    fdtd_common = dict(
        res=40.0, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_220, decay_tol=1e-4, decay_method="energy",
        energy_decay_dt=20.0,
    )
    input_solver = FDTDSolver(**fdtd_common)
    body_solver = EMESolver(
        num_modes=[N_MODES_BODY], num_cells=[1], lengths=[1.0],
        oxide_width=fdtd_common["sim_width"], oxide_height=fdtd_common["sim_height"],
        res=40.0, material_stack=SOI_220, mode_solver="meep",
    )

    mmi = make_mmi(MM_LENGTH)
    half_body = mmi.mm_length / 2
    half_total = mmi.total_length / 2
    boundaries = [-half_total, -half_body + BODY_STUB,
                  +half_body - BODY_STUB, +half_total]

    # Input sections -- length-independent, all four ports, cache-only.
    M_in = [
        input_solver.solve_section(
            mmi, wavelengths, boundaries[0], boundaries[1],
            n_modes_in=1, n_modes_out=N_MODES,
            stub_left=input_solver.waveguide_length, stub_right=0.0,
            input_idx=k, save_dir=save_dir, cache_only=True,
        )
        for k in range(4)
    ]
    M_out = [Mk.transpose(0, 2, 1) for Mk in M_in]
    print(f"[rescore] M_in cached OK: {M_in[0].shape}")

    # Body grids. The length-sweep cache is keyed on the REQUESTED lengths
    # array, so a strided request is a cache MISS even though every point in
    # it was solved -- always ask for the exact full grid the cached runs used
    # and stride the returned array instead. (The calibrated grid was solved
    # from 5.00; its own script's request array is reproduced verbatim below,
    # then sliced onto the common 20.0-400.0 range.)
    t0 = time.time()
    M_body_uncal = body_solver.solve_section_length_sweep(
        mmi, wavelengths, boundaries[1], boundaries[2],
        GRID - 2 * BODY_STUB, n_modes_in=N_MODES_BODY, n_modes_out=N_MODES_BODY,
        save_dir=save_dir, cache_only=True,
    )
    CAL_GRID = np.round(np.arange(5.00, 400.0 + 1e-9, 0.05), 2)
    M_body_cal_full = body_solver.solve_section_length_sweep(
        mmi, wavelengths, boundaries[1], boundaries[2],
        CAL_GRID - 2 * BODY_STUB, n_modes_in=N_MODES_BODY, n_modes_out=N_MODES_BODY,
        neff_override=neff_calibrated, save_dir=save_dir, cache_only=True,
    )
    i0 = int(np.argmin(np.abs(CAL_GRID - GRID[0])))
    assert abs(CAL_GRID[i0] - GRID[0]) < 1e-9, "calibrated grid does not contain 20.00"
    M_body_cal = M_body_cal_full[i0:i0 + len(GRID)]
    assert len(M_body_cal) == len(GRID), "calibrated grid shorter than common grid"
    print(f"[rescore] body grids loaded in {time.time()-t0:.1f}s: "
          f"uncal {M_body_uncal.shape}, cal {M_body_cal.shape}")

    grid = GRID[::STRIDE]
    M_body_uncal = M_body_uncal[::STRIDE]
    M_body_cal = M_body_cal[::STRIDE]
    print(f"[rescore] scoring {len(grid)} lengths (stride {STRIDE}) "
          f"from {grid[0]:.2f} to {grid[-1]:.2f} um")

    results = {}
    for tag, M_body in (("uncal", M_body_uncal), ("cal", M_body_cal)):
        t0 = time.time()
        # coeffs[(i, j)] -> (nL, nwl); one vectorised cascade over the length axis.
        coeffs = {}
        for i in DRIVEN:
            tmp = np.einsum("Lwab,wbk->Lwak", M_body, M_in[i])
            for j in OUTPUTS:
                c = np.einsum("wab,Lwbk->Lwak", M_out[j], tmp)
                coeffs[(i, j)] = c[:, :, 0, 0]

        nL = len(grid)
        min_cmrr = np.empty(nL)
        max_dphi = np.empty(nL)
        max_xloss = np.empty(nL)
        max_imb = np.empty(nL)
        for li in range(nL):
            smat = {}
            for i in DRIVEN:
                s = SMatrix(wavelengths=wavelengths)
                for j in OUTPUTS:
                    s.set(coeffs[(i, j)][li], port=PORT_NAMES[j])
                smat[i] = s
            cm = []
            xl = []
            imb = []
            for i in DRIVEN:
                cm.append(smat[i].cmrr_db("out5", "out8"))
                cm.append(smat[i].cmrr_db("out6", "out7"))
                xl.append(smat[i].excess_loss_db(PORT_NAMES))
                imb.append(smat[i].imbalance_db(PORT_NAMES))
            min_cmrr[li] = float(np.min(np.concatenate(cm)))
            max_xloss[li] = float(np.max(np.concatenate(xl)))
            max_imb[li] = float(np.max(np.abs(np.concatenate(imb))))
            # A11 convention: smat[2] relative to smat[0], exactly as production.
            phi = {j: smat[2].cross_excitation_phase_diff_deg(smat[0], port=PORT_NAMES[j])
                   for j in OUTPUTS}
            dmax = 0.0
            for j in [1, 2, 3]:
                raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
                d = ((raw + 180.0) % 360.0) - 180.0
                dmax = max(dmax, float(np.max(np.abs(d))))
            max_dphi[li] = dmax
        results[tag] = dict(min_cmrr=min_cmrr, max_dphi=max_dphi,
                            max_xloss=max_xloss, max_imb=max_imb)
        print(f"[rescore] {tag}: scored {nL} lengths in {time.time()-t0:.1f}s")

    # ---- report ----------------------------------------------------------
    def report(tag):
        r = results[tag]
        feas = r["max_dphi"] <= SPEC_PHASE_DEG
        print(f"\n=== {tag.upper()} ===")
        print(f"  lengths meeting |dphi| <= {SPEC_PHASE_DEG}deg: "
              f"{int(feas.sum())} / {len(grid)}")
        if feas.any():
            cand = np.where(feas)[0]
            best = cand[np.argmax(r["min_cmrr"][cand])]
            print(f"  BEST feasible: L = {grid[best]:.2f} um   "
                  f"min CMRR = {r['min_cmrr'][best]:.3f} dBe   "
                  f"max|dphi| = {r['max_dphi'][best]:.3f} deg   "
                  f"excess loss = {r['max_xloss'][best]:.3f} dB")
            order = cand[np.argsort(-r["min_cmrr"][cand])][:10]
            print("  top 10 feasible lengths:")
            for k in order:
                print(f"    L={grid[k]:8.2f}  CMRR={r['min_cmrr'][k]:7.3f} dBe"
                      f"  max|dphi|={r['max_dphi'][k]:6.3f}  xloss={r['max_xloss'][k]:6.3f} dB")
        gbest = int(np.argmax(r["min_cmrr"]))
        print(f"  unconstrained CMRR max: L = {grid[gbest]:.2f} um  "
              f"CMRR = {r['min_cmrr'][gbest]:.3f} dBe  "
              f"max|dphi| = {r['max_dphi'][gbest]:.3f} deg")
        # The incumbent, for reference.
        inc = int(np.argmin(np.abs(grid - MM_LENGTH)))
        print(f"  incumbent L={grid[inc]:.2f}: CMRR={r['min_cmrr'][inc]:.3f} dBe  "
              f"max|dphi|={r['max_dphi'][inc]:.3f} deg  "
              f"xloss={r['max_xloss'][inc]:.3f} dB")

    for tag in ("uncal", "cal"):
        report(tag)

    np.savez(os.path.join(OUT_DIR, f"rescore_stride{STRIDE}.npz"),
             grid=grid, wavelengths=wavelengths, mm_length_incumbent=MM_LENGTH,
             **{f"{tag}_{k}": v for tag, r in results.items() for k, v in r.items()})
    print(f"\n[rescore] saved -> {OUT_DIR}/rescore_stride{STRIDE}.npz")


if __name__ == "__main__":
    main()
