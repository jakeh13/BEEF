"""4x4 MMI 90-degree hybrid -- re-check the CURRENT baseline design point
(offset=1.5, ΔL=-2.75, N_MODES_BODY=16, job 11778128) WITH calibrate_
body_neff-style correction, which this device has NEVER used.

validation_markdown/mmi_4x4_90_hybrid_validation.md, "Open decisions" item
5's own follow-up. `HybridSolver(calibrate_body_neff=True)` (beef/solvers/
hybrid.py:134) propagates an EME body on the beta its NEIGHBOURING FDTD
section actually accumulates (measured via `FDTDSolver.section_mode_beta()`)
instead of on MPB's own mode-solved n_eff -- the two differ by ~2.3%/mode
at res 40 (dc_validation.md item 6), and because coupling depends on Δβ
(near 9.4 rad/µm on the coupler, beta/Δβ ~ 99), that 2.3% lands as a ~10%
error on the coupling itself. This device's own hand-rolled cascade
(mmi_4x4_90hybrid_hybrid.py, mmi_4x4_90hybrid_offset_sweep.py, ...) never
goes through `HybridSolver._solve_sections()` at all -- HybridSolver is
only ever instantiated for the geometry plot -- so this calibration path
has structurally never been reached here, at any point (confirmed by grep:
zero hits for "calibrate_body_neff"/"CALIBRATE" across every 4x4-hybrid
script).

Still true as of 2026-08-24, though the REASON narrowed. `beef.MirrorSpec`
(landed that day) lets a mirror-reused matrix be expressed through the
supported API, so "HybridSolver can't express this" is no longer right as a
blanket statement. Checked against this device specifically: the ±1
restriction is not what blocks it -- MirrorSpec requires strictly ±1
mode_signs, and this device's mirror relations measured R≈1.0000 with ≤0.15
deg deviation from an exact 0/180 sign (2026-08-13, circular statistics,
16/16 modes, both pairs), because the 90-degree relationship is between
output ports for a given input, not between mirror-partner input
excitations. What blocks it is the N-PORT shape: `HybridSolver.solve()`
returns one SMatrix under one `port` label, while this device needs eight
port pairs (2 driven inputs x 4 outputs). So the hand-rolled cascade stays,
and with it this script's hand-rolled calibration. Given this device's headline failure IS phase error (Δφ_j), the
exact quantity most directly driven by Δβ over the body length, this is a
real, previously-unexamined candidate explanation -- possibly a modeling
artifact, not (only) a real device limitation.

Mechanism replicated by hand here, following HybridSolver._calibrated_neff()
(beef/solvers/hybrid.py:1250) exactly:
  1. beta = input_solver.section_mode_beta(mmi, wavelengths, boundaries[1],
     boundaries[2], N_MODES, lengths=None, save_dir=save_dir) -- measures
     the FDTD-consistent beta for the BODY's own cross-section (the input
     solver's settings -- res, decay_tol, material_stack -- apply, since
     it is asking "what would meep's time-stepper see propagating through
     THIS cross-section", not anything about the input taper). lengths=None
     lets section_mode_beta() pick its own safety-capped default (3 short,
     evenly-spaced lengths -- cheap relative to a full production input
     solve, since length is short even though the transverse domain is the
     full body width).
  2. neff = beta * wavelengths[:, None] / (2*pi) -- exact conversion
     HybridSolver itself uses.
  3. Re-run the body length-sweep WITH neff_override=neff -- mode SHAPES
     (hence every interface overlap/coupling coefficient) are untouched
     (EMESolver.solve_section_length_sweep()'s own docstring), only the
     per-length propagation phase changes. This IS a genuinely new EME
     solve (neff_override is part of every cache tier's key -- a
     calibrated and uncalibrated run never share a cache entry), so
     expect comparable cost to the original body solve (job 11768847,
     ~77-90 min), not free.

Input matrices are NOT re-solved -- neff_override only touches the body;
the offset=1.5 baseline's 4 real input-section solves (job 11778128) are
reused from cache unchanged, exactly like every other offset/ΔL
investigation this session.

Run (SLURM, real new FDTD/EME cost -- see the batch script):
    python mmi_4x4_90hybrid_calibrated_baseline.py
Env override for a cheap smoke test:
    MMI4X4_CAL_RES, MMI4X4_CAL_DECAY_TOL, MMI4X4_CAL_NWL
"""
import os
import numpy as np

from beef import SMatrix
from beef.cache.store import CacheMissError

import mmi_4x4_90hybrid_offset_sweep as base

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "calibrated_baseline"))
os.makedirs(output_dir, exist_ok=True)

RESOLUTION = float(os.environ.get("MMI4X4_CAL_RES", 40.0))
DECAY_TOL = float(os.environ.get("MMI4X4_CAL_DECAY_TOL", 1e-11))
N_WL = int(os.environ.get("MMI4X4_CAL_NWL", 3))
wavelengths = np.linspace(base.WL_MIN, base.WL_MAX, N_WL)

DRIVEN = [0, 2]
OUTPUTS = [0, 1, 2, 3]
PORT_NAMES = ["out5", "out6", "out7", "out8"]
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}


def make_solvers():
    fdtd_common = dict(
        res=RESOLUTION, sim_width=base.MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=base.SOI_220, decay_tol=DECAY_TOL,
    )
    from beef import FDTDSolver, EMESolver
    input_solver = FDTDSolver(**fdtd_common)
    # SEPARATE solver, decay_method="energy" -- validation doc 2026-08-18:
    # stop_when_dft_decayed() (input_solver's default) stalls indefinitely
    # on this body cross-section's mode 1 at production res=40 (past
    # 660,000+ timesteps, never converging -- confirmed not a mode-identity
    # issue, see the doc's own mode-1 investigation). stop_when_energy_
    # decayed() converges normally (t=140.1) at the identical settings.
    # MUST be a separate instance from input_solver: decay_method is now
    # part of the cache fingerprint (beef/cache/config.py, added the same
    # day) specifically so an "energy" run never silently reuses/pollutes
    # a "dft" cache entry -- setting it on input_solver directly would also
    # change the fingerprint for the already-cached (job 11778128) baseline
    # input-matrix fetches below, turning free cache hits into wasted
    # re-solves, the exact bug job 11938100/11981048 already hit once.
    calibration_solver = FDTDSolver(**{**fdtd_common, "decay_method": "energy", "energy_decay_dt": 20.0})
    body_solver = EMESolver(
        num_modes=[base.N_MODES], num_cells=[1], lengths=[1.0],
        oxide_width=fdtd_common["sim_width"], oxide_height=fdtd_common["sim_height"],
        res=RESOLUTION, material_stack=base.SOI_220, mode_solver="meep",
    )
    return input_solver, calibration_solver, body_solver


def metrics_for(M_in, M_body_one, wavelengths):
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
        smat[i].check_passive()

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
    return dict(max_dphi=max_dphi, min_cmrr=min_cmrr,
                max_excess_loss=max_excess_loss, max_imbalance=max_imbalance,
                imbalance=imbalance, excess_loss=excess_loss, cmrr_I=cmrr_I, cmrr_Q=cmrr_Q,
                dphi_6=dphi[1], dphi_7=dphi[2], dphi_8=dphi[3])


def main():
    input_solver, calibration_solver, body_solver = make_solvers()
    mmi = base.make_mmi(base.BASELINE_OFFSET)
    boundaries, body_stub = base.boundaries_for(mmi)

    print(f"=== Baseline (offset={base.BASELINE_OFFSET}, ΔL={base.DELTA_L}) "
          f"WITH FDTD-calibrated body neff ===")

    # 1. Real 4 input matrices -- reused from cache unchanged (job 11778128
    # at production settings; smoke-test settings will genuinely miss and
    # need their own cheap fresh solves).
    try:
        M_in = [base.solve_input(input_solver, mmi, boundaries, k, cache_only=True) for k in range(4)]
        print("  input matrices: cache HIT (all 4, no new FDTD)")
    except CacheMissError:
        print("  input matrices: cache MISS (expected at non-production settings) -- solving fresh")
        M_in = [base.solve_input(input_solver, mmi, boundaries, k, cache_only=False) for k in range(4)]

    # 2. Calibration: FDTD-measured beta for the BODY's own cross-section,
    # via the DEDICATED calibration_solver (decay_method="energy" -- see
    # make_solvers()'s own comment for why this must be a separate
    # instance from input_solver), following HybridSolver._calibrated_neff()
    # in every other respect.
    x0, x1 = boundaries[1], boundaries[2]
    print(f"  measuring FDTD-consistent beta for the body cross-section "
          f"x=[{x0:.4f},{x1:.4f}] (N_MODES={base.N_MODES}, "
          f"decay_method={calibration_solver.decay_method!r}) ...")
    beta = calibration_solver.section_mode_beta(
        mmi, wavelengths, x0, x1, base.N_MODES,
        lengths=None, save_dir=base.save_dir,
    )
    neff_calibrated = beta * wavelengths[:, None] / (2.0 * np.pi)
    print(f"  calibrated neff[0] (first mode, per wl): {neff_calibrated[:, 0]}")

    # 3. Uncalibrated body -- fetched from cache if available (production
    # baseline), for a direct side-by-side comparison.
    _DELTA_L_GRID = np.linspace(-10.0, 10.0, 81)
    _body_lengths_grid = (0.75 * base.L_PI + _DELTA_L_GRID) - 2 * body_stub
    _Li = int(np.argmin(np.abs(_DELTA_L_GRID - base.DELTA_L)))
    try:
        M_body_uncal_grid = body_solver.solve_section_length_sweep(
            mmi, wavelengths, x0, x1, _body_lengths_grid,
            n_modes_in=base.N_MODES, n_modes_out=base.N_MODES,
            save_dir=base.save_dir, cache_only=True,
        )
        uncal = metrics_for(M_in, M_body_uncal_grid[_Li:_Li + 1], wavelengths)
        print(f"  UNCALIBRATED (cache): max|Δφ|={uncal['max_dphi']:.3f} deg  "
              f"min_CMRR={uncal['min_cmrr']:.3f} dBe")
    except CacheMissError:
        uncal = None
        print("  UNCALIBRATED: cache MISS at these settings, skipping direct comparison")

    # 4. Calibrated body -- genuinely new EME solve (neff_override changes
    # the cache key; cannot reuse job 11768847's uncalibrated grid).
    print(f"  solving CALIBRATED body length sweep ({len(_body_lengths_grid)} lengths, "
          f"one mode solve/wavelength, neff_override applied) ...")
    M_body_cal_grid = body_solver.solve_section_length_sweep(
        mmi, wavelengths, x0, x1, _body_lengths_grid,
        n_modes_in=base.N_MODES, n_modes_out=base.N_MODES,
        neff_override=neff_calibrated, save_dir=base.save_dir,
    )
    cal = metrics_for(M_in, M_body_cal_grid[_Li:_Li + 1], wavelengths)
    print(f"  CALIBRATED @ ΔL={base.DELTA_L}: max|Δφ|={cal['max_dphi']:.3f} deg  "
          f"min_CMRR={cal['min_cmrr']:.3f} dBe  max_XL={cal['max_excess_loss']:.3f} dB  "
          f"max_imb={cal['max_imbalance']:.3f} dB")

    # 5. Also find calibrated-body's OWN best ΔL over the full grid (cheap,
    # already have the full grid) -- calibration could plausibly shift the
    # optimum length too, not just the value at the old optimum.
    best_dphi = np.inf
    best_li = _Li
    for li in range(len(_DELTA_L_GRID)):
        m = metrics_for(M_in, M_body_cal_grid[li:li + 1], wavelengths)
        if m["max_dphi"] < best_dphi:
            best_dphi = m["max_dphi"]
            best_li = li
    best = metrics_for(M_in, M_body_cal_grid[best_li:best_li + 1], wavelengths)
    print(f"  CALIBRATED best ΔL over full grid: ΔL={_DELTA_L_GRID[best_li]:+.2f}  "
          f"max|Δφ|={best['max_dphi']:.3f} deg  min_CMRR={best['min_cmrr']:.3f} dBe")

    np.savez(os.path.join(output_dir, "calibrated_baseline_results.npz"),
             delta_L_grid=_DELTA_L_GRID, neff_calibrated=neff_calibrated,
             uncal_max_dphi=(uncal["max_dphi"] if uncal else np.nan),
             uncal_min_cmrr=(uncal["min_cmrr"] if uncal else np.nan),
             cal_at_old_deltaL_max_dphi=cal["max_dphi"], cal_at_old_deltaL_min_cmrr=cal["min_cmrr"],
             cal_best_deltaL=_DELTA_L_GRID[best_li], cal_best_max_dphi=best["max_dphi"],
             cal_best_min_cmrr=best["min_cmrr"])
    print(f"\nsaved -> {output_dir}/calibrated_baseline_results.npz")


if __name__ == "__main__":
    main()
