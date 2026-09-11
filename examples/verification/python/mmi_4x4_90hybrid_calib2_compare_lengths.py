"""4x4 MMI 90-degree hybrid -- with-vs-without n_eff calibration comparison,
at the CURRENT (optimised) design point, across the lengths this session
has been probing (5.00 / 10.00 / 15.00 / 159.30 um production).

validation_markdown/mmi_4x4_90_hybrid_validation.md follow-up, 2026-09-11.
Requires mmi_4x4_90hybrid_calib2_mode_worker.py (x16, SLURM array) +
mmi_4x4_90hybrid_calib2_aggregate.py to have already produced
example_output/mmi_4x4_90hybrid/calib2_per_mode/aggregated_beta_n16.npz.

Calibration is a pure per-cross-section, length-independent correction
(EMESolver.solve_section_length_sweep's own neff_override enters only the
propagation phase exp(2j*pi*neff/wl*L), never the mode basis) -- so ONE
aggregated neff_calibrated array applies at every length below without
re-measurement. This script:

  1. Runs ONE fresh (uncached -- neff_override changes the config
     fingerprint) calibrated body solve, on a length grid spanning ALL
     FOUR target lengths at once (5.00-400.00 um), so the expensive part
     (meow's mode-basis solve) happens exactly once and every target
     length is a free slice of the same result afterward.
  2. Reuses the UNCALIBRATED body solve's own EXISTING cache, one
     cache_only=True call per length, matching exactly the grid each
     length's own prior run (this session's short-body probes, or the
     production script at 159.30) already computed -- zero new compute.
  3. For each length, combines both M_body variants with the SAME cached
     M_in[0]/M_in[2] (ports 1 and 3 -- the paper's own metrics, also
     already cached at every one of these lengths from prior work this
     session) via the identical cascade einsum the production hybrid
     script uses, and prints CMRR/Delta-phi/excess-loss/imbalance
     side by side.

Run:
    python mmi_4x4_90hybrid_calib2_compare_lengths.py
"""
import os

import numpy as np

from beef import MMI, EMESolver, FDTDSolver, SMatrix, SOI_220
from beef.cache.store import CacheMissError

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "mmi_4x4_90hybrid_calib2_compare")
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.join(_ROOT, "example_sim_cache")

CAL_PATH = os.path.join(_ROOT, "example_output", "mmi_4x4_90hybrid", "calib2_per_mode", "aggregated_beta_n16.npz")

# Pinned to mmi_4x4_90hybrid_hybrid.py's own settled CURRENT design point.
WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 8.5
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 1.30
TAPER_LENGTH = 4.0
N_MODES = 16
BODY_STUB = 0.25
WL_MIN, WL_MAX, N_WL = 1.29, 1.33, 3
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

TARGET_LENGTHS = [5.00, 10.00, 15.00, 159.30]
PORT_NAMES = ["out5", "out6", "out7", "out8"]
DRIVEN = [0, 2]   # paper's ports 1 and 3
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}


def make_mmi(mm_length):
    return MMI(
        n_inputs=4, n_outputs=4,
        mm_length=mm_length, mm_width=MM_WIDTH, waveguide_width=WAVEGUIDE_WIDTH,
        taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
        waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
        use_neighbor_aperture=True,
    )


def make_solvers():
    fdtd_common = dict(
        res=40.0, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_220, decay_tol=1e-4, decay_method="energy", energy_decay_dt=20.0,
    )
    input_solver = FDTDSolver(**fdtd_common)
    body_solver = EMESolver(
        num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
        oxide_width=fdtd_common["sim_width"], oxide_height=fdtd_common["sim_height"],
        res=40.0, material_stack=SOI_220, mode_solver="meep",
    )
    return input_solver, body_solver


def metrics_from_cascade(M_in, M_out, M_body_slice, wavelengths):
    """M_body_slice: [1, n_wl, N_MODES, N_MODES] (one length). Returns dict
    of the paper's own metrics, identical math to the production script."""
    coeffs = {}
    for i in DRIVEN:
        tmp = np.einsum("Lwab,wbk->Lwak", M_body_slice, M_in[i])
        for j in range(4):
            c = np.einsum("wab,Lwbk->Lwak", M_out[j], tmp)
            coeffs[(i, j)] = c[0, :, 0, 0]

    smat = {i: SMatrix(wavelengths=wavelengths) for i in DRIVEN}
    for i in DRIVEN:
        for j in range(4):
            smat[i].set(coeffs[(i, j)], port=PORT_NAMES[j])
        smat[i].check_passive()

    imbalance = {i: smat[i].imbalance_db(PORT_NAMES) for i in DRIVEN}
    excess_loss = {i: smat[i].excess_loss_db(PORT_NAMES) for i in DRIVEN}
    cmrr_I = {i: smat[i].cmrr_db("out5", "out8") for i in DRIVEN}
    cmrr_Q = {i: smat[i].cmrr_db("out6", "out7") for i in DRIVEN}
    phi = {j: smat[2].cross_excitation_phase_diff_deg(smat[0], port=PORT_NAMES[j]) for j in range(4)}
    dphi = {}
    for j in [1, 2, 3]:
        raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
        dphi[j] = ((raw + 180.0) % 360.0) - 180.0

    max_dphi = float(np.max(np.abs([dphi[1], dphi[2], dphi[3]])))
    min_cmrr = float(np.min([cmrr_I[0], cmrr_I[2], cmrr_Q[0], cmrr_Q[2]]))
    max_xl = float(np.max([excess_loss[0], excess_loss[2]]))
    max_imb = float(np.max([imbalance[0], imbalance[2]]))
    return dict(max_dphi=max_dphi, min_cmrr=min_cmrr, max_xl=max_xl, max_imb=max_imb)


def main():
    if not os.path.exists(CAL_PATH):
        raise SystemExit(f"missing {CAL_PATH} -- run the calib2 worker array + aggregate first.")
    cal = np.load(CAL_PATH)
    neff_calibrated = cal["neff_calibrated"]   # [n_wl, N_MODES]
    assert neff_calibrated.shape == (N_WL, N_MODES), neff_calibrated.shape
    print(f"loaded calibration: neff_calibrated[0]={neff_calibrated[0]}")

    input_solver, body_solver = make_solvers()

    # Use the longest target length's geometry only to fix the (x0, x1)
    # calibration cross-section span (boundaries[1]/[2] scale with
    # mm_length, but the body cross-section itself -- and hence the mode
    # basis / calibration -- is identical regardless of which mm_length
    # the geometry object is built with; only the ABSOLUTE boundary
    # positions differ, and solve_section_length_sweep is called fresh
    # per target mm_length below anyway).
    mmi_cal = make_mmi(159.30)
    half_body_cal = mmi_cal.mm_length / 2
    x0_cal, x1_cal = -half_body_cal + BODY_STUB, +half_body_cal - BODY_STUB

    cal_grid = np.round(np.arange(5.00, 400.0 + 1e-9, 0.05), 2)
    cal_body_lengths = cal_grid - 2 * BODY_STUB
    print(f"solving CALIBRATED body over unified grid {cal_grid[0]}-{cal_grid[-1]} um "
          f"({len(cal_grid)} points) -- fresh, neff_override changes the cache fingerprint...")
    M_body_cal_grid = body_solver.solve_section_length_sweep(
        mmi_cal, wavelengths, x0_cal, x1_cal, cal_body_lengths,
        n_modes_in=N_MODES, n_modes_out=N_MODES, neff_override=neff_calibrated,
        save_dir=save_dir,
    )
    print(f"  calibrated body solved: {M_body_cal_grid.shape}")

    results = {}
    for L in TARGET_LENGTHS:
        mmi = make_mmi(L)
        half_body = mmi.mm_length / 2
        half_total = mmi.total_length / 2
        boundaries = [-half_total, -half_body + BODY_STUB, +half_body - BODY_STUB, +half_total]

        M_in = []
        for k in range(4):
            Mk = input_solver.solve_section(
                mmi, wavelengths, boundaries[0], boundaries[1],
                n_modes_in=1, n_modes_out=N_MODES,
                stub_left=input_solver.waveguide_length, stub_right=0.0,
                input_idx=k, save_dir=save_dir, cache_only=True,
            )
            M_in.append(Mk)
        M_out = [Mk.transpose(0, 2, 1) for Mk in M_in]

        # Calibrated: slice from the unified grid solved above.
        Li_cal = int(np.argmin(np.abs(cal_grid - L)))
        M_body_cal = M_body_cal_grid[Li_cal:Li_cal + 1]
        m_cal = metrics_from_cascade(M_in, M_out, M_body_cal, wavelengths)

        # Uncalibrated: reuse each length's own existing cache, same grid
        # convention every prior script at this length used.
        grid_lo = 20.0 if L >= 20.0 else L
        uncal_grid = np.round(np.arange(grid_lo, 400.0 + 1e-9, 0.05), 2)
        uncal_body_lengths = uncal_grid - 2 * BODY_STUB
        M_body_uncal_grid = body_solver.solve_section_length_sweep(
            mmi, wavelengths, boundaries[1], boundaries[2], uncal_body_lengths,
            n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir, cache_only=True,
        )
        Li_uncal = int(np.argmin(np.abs(uncal_grid - L)))
        M_body_uncal = M_body_uncal_grid[Li_uncal:Li_uncal + 1]
        m_uncal = metrics_from_cascade(M_in, M_out, M_body_uncal, wavelengths)

        results[L] = dict(uncal=m_uncal, cal=m_cal)
        print(f"\nL={L:.2f} um:")
        print(f"  UNCALIBRATED: max|dphi|={m_uncal['max_dphi']:.3f} deg  min_CMRR={m_uncal['min_cmrr']:.3f} dBe  "
              f"max_XL={m_uncal['max_xl']:.3f} dB  max_imb={m_uncal['max_imb']:.3f} dB")
        print(f"  CALIBRATED:   max|dphi|={m_cal['max_dphi']:.3f} deg  min_CMRR={m_cal['min_cmrr']:.3f} dBe  "
              f"max_XL={m_cal['max_xl']:.3f} dB  max_imb={m_cal['max_imb']:.3f} dB")

    print("\n\n| L (um) | max|dphi| uncal | max|dphi| cal | min CMRR uncal | min CMRR cal | max XL uncal | max XL cal | max imb uncal | max imb cal |")
    print("|---|---|---|---|---|---|---|---|---|")
    for L in TARGET_LENGTHS:
        u, c = results[L]["uncal"], results[L]["cal"]
        print(f"| {L:.2f} | {u['max_dphi']:.3f} | {c['max_dphi']:.3f} | {u['min_cmrr']:.3f} | {c['min_cmrr']:.3f} | "
              f"{u['max_xl']:.3f} | {c['max_xl']:.3f} | {u['max_imb']:.3f} | {c['max_imb']:.3f} |")

    np.savez(os.path.join(output_dir, "calib2_compare_lengths.npz"),
             lengths=np.array(TARGET_LENGTHS),
             **{f"L{L:.2f}_{variant}_{k}": v
                for L in TARGET_LENGTHS for variant in ("uncal", "cal") for k, v in results[L][variant].items()})
    print(f"\nsaved -> {output_dir}/calib2_compare_lengths.npz")


if __name__ == "__main__":
    main()
