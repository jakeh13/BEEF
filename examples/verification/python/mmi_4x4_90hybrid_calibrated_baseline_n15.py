"""4x4 MMI 90-degree hybrid -- calibrated-vs-uncalibrated body at N_MODES=15,
dropping mode 16 (whose own FDTD run does not converge at all).

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-20. The
16-mode calibration cannot be completed: mode 16's own isolated run does not
merely decay slowly, it DIVERGES (energy ratio bottoms at 1.95e-11 near
t~260, just above decay_tol=1e-11, then grows monotonically to 1.7e53 by
t~2081) -- an undamped-Lorentzian (`gamma=0`, beef/materials/stack.py:284)
late-time instability, not a stopping-criterion problem. Modes 1-15 all
converged normally. This script takes the user's decision to proceed on the
15 modes that exist rather than chase mode 16 further.

Two comparisons, deliberately BOTH, because "N=15 calibrated" differs from
the published baseline in TWO ways at once and conflating them would be a
wrong reading:

  A. calibrated-15 vs uncalibrated-15  -- isolates the CALIBRATION effect
     alone (identical 15-mode basis on both sides). This is the number that
     answers "does FDTD-consistent beta change this device's answer?"
  B. calibrated-15 vs uncalibrated-16  -- vs the published production
     baseline (job 11778128). Confounds calibration with the dropped mode,
     so it is reported as context, NOT as the calibration effect.

The 4 input FDTD matrices are NOT re-solved at 15 modes. They are fetched
from job 11778128's cache at 16 modes and TRUNCATED to the first 15
(`M[:, :15, :]`). This relies on the same per-band independence the repo
already relies on for `measure_diagonal_only` and for the per-mode
calibration split: meep decomposes each band separately, so bands 1-15 do
not depend on whether band 16 was also requested. Re-solving at
n_modes_out=15 would be a genuine cache miss costing ~10h of FDTD for a
result expected to be bit-comparable -- not spent without reason.

Run:
    python mmi_4x4_90hybrid_calibrated_baseline_n15.py
"""
import os
import sys
import numpy as np

from beef import EMESolver
from beef.cache.store import CacheMissError

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmi_4x4_90hybrid_offset_sweep as base
import mmi_4x4_90hybrid_calibrated_baseline as cal_base

N_CAL_MODES = int(os.environ.get("MMI4X4_CAL_MODES", 15))

_HERE = os.path.dirname(os.path.abspath(__file__))
BETA_PATH = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid",
    "calibration_per_mode", f"aggregated_beta_n{N_CAL_MODES:02d}.npz"))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid",
    "calibrated_baseline"))
os.makedirs(output_dir, exist_ok=True)


def main():
    if not os.path.exists(BETA_PATH):
        raise SystemExit(f"missing {BETA_PATH} -- run mmi_4x4_90hybrid_calibration_"
                         f"aggregate.py with MMI4X4_CAL_MODES={N_CAL_MODES} first.")
    d = np.load(BETA_PATH)
    neff_calibrated = d["neff_calibrated"]
    wavelengths = d["wavelengths"]
    assert neff_calibrated.shape == (len(wavelengths), N_CAL_MODES), \
        f"aggregated beta has shape {neff_calibrated.shape}, expected {(len(wavelengths), N_CAL_MODES)}"

    input_solver, _calibration_solver, _body16 = cal_base.make_solvers()
    body_solver = EMESolver(
        num_modes=[N_CAL_MODES], num_cells=[1], lengths=[1.0],
        oxide_width=base.MM_WIDTH + 3.0, oxide_height=2.0,
        res=cal_base.RESOLUTION, material_stack=base.SOI_220, mode_solver="meep",
    )

    mmi = base.make_mmi(base.BASELINE_OFFSET)
    boundaries, body_stub = base.boundaries_for(mmi)
    x0, x1 = boundaries[1], boundaries[2]

    print(f"=== Baseline (offset={base.BASELINE_OFFSET}, dL={base.DELTA_L}), "
          f"N_MODES={N_CAL_MODES} (mode 16 dropped) ===")
    print(f"  calibrated neff (first wl): {neff_calibrated[0]}")

    # 4 input matrices at 16 modes from cache, truncated to N_CAL_MODES.
    M_in16 = [base.solve_input(input_solver, mmi, boundaries, k, cache_only=True)
              for k in range(4)]
    print(f"  input matrices: cache HIT (all 4 at {base.N_MODES} modes, no new FDTD)")
    M_in = [M[:, :N_CAL_MODES, :] for M in M_in16]

    _DELTA_L_GRID = np.linspace(-10.0, 10.0, 81)
    _body_lengths_grid = (0.75 * base.L_PI + _DELTA_L_GRID) - 2 * body_stub
    _Li = int(np.argmin(np.abs(_DELTA_L_GRID - base.DELTA_L)))

    # Comparison B's reference: the published 16-mode uncalibrated grid.
    try:
        g16 = _body16.solve_section_length_sweep(
            mmi, wavelengths, x0, x1, _body_lengths_grid,
            n_modes_in=base.N_MODES, n_modes_out=base.N_MODES,
            save_dir=base.save_dir, cache_only=True,
        )
        uncal16 = cal_base.metrics_for(M_in16, g16[_Li:_Li + 1], wavelengths)
        print(f"  [B ref] UNCALIBRATED n={base.N_MODES} (cache): "
              f"max|dphi|={uncal16['max_dphi']:.3f} deg  min_CMRR={uncal16['min_cmrr']:.3f} dBe")
    except CacheMissError:
        uncal16 = None
        print(f"  [B ref] UNCALIBRATED n={base.N_MODES}: cache MISS, comparison B unavailable")

    # Comparison A, uncalibrated half -- genuinely new (15-mode basis).
    print(f"  solving UNCALIBRATED body sweep at n={N_CAL_MODES} "
          f"({len(_body_lengths_grid)} lengths) ...")
    g_uncal = body_solver.solve_section_length_sweep(
        mmi, wavelengths, x0, x1, _body_lengths_grid,
        n_modes_in=N_CAL_MODES, n_modes_out=N_CAL_MODES,
        save_dir=base.save_dir,
    )
    uncal = cal_base.metrics_for(M_in, g_uncal[_Li:_Li + 1], wavelengths)
    print(f"  UNCALIBRATED n={N_CAL_MODES}: max|dphi|={uncal['max_dphi']:.3f} deg  "
          f"min_CMRR={uncal['min_cmrr']:.3f} dBe  max_XL={uncal['max_excess_loss']:.3f} dB  "
          f"max_imb={uncal['max_imbalance']:.3f} dB")

    # Comparison A, calibrated half.
    print(f"  solving CALIBRATED body sweep at n={N_CAL_MODES} (neff_override) ...")
    g_cal = body_solver.solve_section_length_sweep(
        mmi, wavelengths, x0, x1, _body_lengths_grid,
        n_modes_in=N_CAL_MODES, n_modes_out=N_CAL_MODES,
        neff_override=neff_calibrated, save_dir=base.save_dir,
    )
    cal = cal_base.metrics_for(M_in, g_cal[_Li:_Li + 1], wavelengths)
    print(f"  CALIBRATED   n={N_CAL_MODES}: max|dphi|={cal['max_dphi']:.3f} deg  "
          f"min_CMRR={cal['min_cmrr']:.3f} dBe  max_XL={cal['max_excess_loss']:.3f} dB  "
          f"max_imb={cal['max_imbalance']:.3f} dB")

    # Calibration can shift the optimum length, not just the value at it.
    best_li, best_dphi = _Li, np.inf
    for li in range(len(_DELTA_L_GRID)):
        m = cal_base.metrics_for(M_in, g_cal[li:li + 1], wavelengths)
        if m["max_dphi"] < best_dphi:
            best_dphi, best_li = m["max_dphi"], li
    best = cal_base.metrics_for(M_in, g_cal[best_li:best_li + 1], wavelengths)
    print(f"  CALIBRATED best dL over full grid: dL={_DELTA_L_GRID[best_li]:+.2f}  "
          f"max|dphi|={best['max_dphi']:.3f} deg  min_CMRR={best['min_cmrr']:.3f} dBe")

    print(f"\n  [A] calibration effect at n={N_CAL_MODES}, same basis both sides:")
    print(f"      max|dphi|  {uncal['max_dphi']:.3f} -> {cal['max_dphi']:.3f} deg "
          f"({cal['max_dphi'] - uncal['max_dphi']:+.3f})")
    print(f"      min_CMRR   {uncal['min_cmrr']:.3f} -> {cal['min_cmrr']:.3f} dBe "
          f"({cal['min_cmrr'] - uncal['min_cmrr']:+.3f})")
    if uncal16 is not None:
        print(f"  [B] vs published n={base.N_MODES} uncalibrated (confounds calibration "
              f"with the dropped mode -- context only):")
        print(f"      max|dphi|  {uncal16['max_dphi']:.3f} -> {cal['max_dphi']:.3f} deg")
        print(f"      min_CMRR   {uncal16['min_cmrr']:.3f} -> {cal['min_cmrr']:.3f} dBe")
        print(f"  [B'] mode-count effect alone (uncal 16 -> uncal 15): "
              f"max|dphi| {uncal16['max_dphi']:.3f} -> {uncal['max_dphi']:.3f} deg")

    np.savez(os.path.join(output_dir, f"calibrated_baseline_n{N_CAL_MODES:02d}_results.npz"),
             n_modes=N_CAL_MODES, delta_L_grid=_DELTA_L_GRID,
             neff_calibrated=neff_calibrated,
             uncal15_max_dphi=uncal["max_dphi"], uncal15_min_cmrr=uncal["min_cmrr"],
             uncal15_max_XL=uncal["max_excess_loss"], uncal15_max_imb=uncal["max_imbalance"],
             cal15_max_dphi=cal["max_dphi"], cal15_min_cmrr=cal["min_cmrr"],
             cal15_max_XL=cal["max_excess_loss"], cal15_max_imb=cal["max_imbalance"],
             uncal16_max_dphi=(uncal16["max_dphi"] if uncal16 else np.nan),
             uncal16_min_cmrr=(uncal16["min_cmrr"] if uncal16 else np.nan),
             cal_best_deltaL=_DELTA_L_GRID[best_li], cal_best_max_dphi=best["max_dphi"],
             cal_best_min_cmrr=best["min_cmrr"])
    print(f"\nsaved -> {output_dir}/calibrated_baseline_n{N_CAL_MODES:02d}_results.npz")


if __name__ == "__main__":
    main()
