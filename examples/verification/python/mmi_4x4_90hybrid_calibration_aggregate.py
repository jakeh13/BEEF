"""4x4 MMI 90-degree hybrid -- combine the 16 per-mode calibration workers
(mmi_4x4_90hybrid_calibration_mode_worker.py) into the same BetaCalibration
result section_mode_beta_details() would have produced from one process.

validation_markdown/mmi_4x4_90_hybrid_validation.md follow-up, 2026-08-20.
See the worker script's own docstring for why per-mode isolation is exact,
not an approximation: `beta_fit_from_phases()` (the same, unmodified,
imported library function used here) is called once per (wavelength,
mode) cell using only that mode's own phase-vs-length data and seed_beta
-- reproducing this per mode in a separate process changes nothing about
what gets computed, only how many processes compute it.

Applies the SAME two post-fit guards section_mode_beta_details() applies
(residual, branch margin) using the SAME imported thresholds, so a
mode that would have failed those guards in a single-process run fails
them here too, not silently passed through.

Run (after all 16 mode_XX.npz files exist):
    python mmi_4x4_90hybrid_calibration_aggregate.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmi_4x4_90hybrid_offset_sweep as base
from beef.solvers.fdtd import (
    beta_fit_from_phases, _BETA_FIT_RESIDUAL_MAX, _BETA_BRANCH_MARGIN_RAISE,
    _BETA_SEED_EXCESS_WARN,
)

IN_DIR = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..",
    "example_output", "mmi_4x4_90hybrid", "calibration_per_mode"))

# Modes 1..N_CAL_MODES only. Defaults to the full basis; set to 15 to drop
# mode 16, whose own FDTD run does not converge at all (2026-08-20: an
# undamped-Lorentzian late-time instability, energy GROWING to 1.7e53 by
# t~2081 rather than decaying -- validation doc's own correction section).
N_CAL_MODES = int(os.environ.get("MMI4X4_CAL_MODES", base.N_MODES))


def main():
    n_modes = N_CAL_MODES
    wavelengths = None
    lengths_ref = None
    phases = None       # [n_lengths, n_wl, n_modes]
    seed_beta = None    # [n_wl, n_modes]

    missing = []
    for k in range(1, n_modes + 1):
        path = os.path.join(IN_DIR, f"mode_{k:02d}.npz")
        if not os.path.exists(path):
            missing.append(k)
            continue
        d = np.load(path)
        if wavelengths is None:
            wavelengths = d["wavelengths"]
            lengths_ref = d["lengths"]
            n_lengths = len(lengths_ref)
            phases = np.zeros((n_lengths, len(wavelengths), n_modes))
            seed_beta = np.zeros((len(wavelengths), n_modes))
        else:
            assert np.allclose(d["lengths"], lengths_ref), f"mode {k}: lengths mismatch"
            assert np.allclose(d["wavelengths"], wavelengths), f"mode {k}: wavelengths mismatch"
        phases[:, :, k - 1] = d["phases"]
        seed_beta[:, k - 1] = d["seed_beta"]

    if missing:
        raise SystemExit(f"missing per-mode files for modes {missing} -- run the worker for each first.")

    print(f"loaded all {n_modes} modes' calibration data, lengths={lengths_ref}")

    beta = np.zeros((len(wavelengths), n_modes))
    grid_resid = np.zeros((len(wavelengths), n_modes))
    grid_margin = np.zeros((len(wavelengths), n_modes))
    grid_excess = np.zeros((len(wavelengths), n_modes))
    grid_turn_amb = np.zeros((len(wavelengths), n_modes))

    for i in range(len(wavelengths)):
        for m in range(n_modes):
            fit = beta_fit_from_phases(list(lengths_ref), phases[:, i, m], seed_beta[i, m])
            beta[i, m] = fit.beta
            grid_resid[i, m] = fit.resid
            grid_margin[i, m] = fit.margin
            grid_excess[i, m] = fit.seed_excess
            grid_turn_amb[i, m] = fit.turn_ambiguity

    def _cell(grid):
        i, m = np.unravel_index(int(np.argmax(np.abs(grid))), grid.shape)
        return float(grid[i, m]), f"λ={wavelengths[i]:.4f} µm mode {m}"

    worst_resid, at_resid = _cell(grid_resid)
    worst_margin, at_margin = _cell(grid_margin)
    worst_excess, at_excess = _cell(grid_excess)

    print(f"worst residual {worst_resid:.3e} rad at {at_resid} (guard: {_BETA_FIT_RESIDUAL_MAX})")
    print(f"worst branch margin {worst_margin:.3f} at {at_margin} (guard: {_BETA_BRANCH_MARGIN_RAISE})")
    print(f"worst seed excess {100*worst_excess:+.2f}% at {at_excess} (warn: {100*_BETA_SEED_EXCESS_WARN:.1f}%)")

    guard_failures = []
    if worst_resid > _BETA_FIT_RESIDUAL_MAX:
        guard_failures.append(f"residual {worst_resid:.3f} exceeds {_BETA_FIT_RESIDUAL_MAX} at {at_resid}")
    if worst_margin > _BETA_BRANCH_MARGIN_RAISE:
        guard_failures.append(f"branch margin {worst_margin:.3f} exceeds {_BETA_BRANCH_MARGIN_RAISE} at {at_margin}")
    if guard_failures:
        for g in guard_failures:
            print(f"GUARD FAILED: {g}")
        raise SystemExit("aggregation guards failed -- beta below is NOT usable for propagation.")

    neff_calibrated = beta * wavelengths[:, None] / (2.0 * np.pi)
    out_path = os.path.join(IN_DIR, f"aggregated_beta_n{n_modes:02d}.npz")
    np.savez(out_path, beta=beta, neff_calibrated=neff_calibrated, wavelengths=wavelengths,
             grid_resid=grid_resid, grid_margin=grid_margin, grid_excess=grid_excess)
    print(f"beta[0] (first wavelength, all modes): {beta[0]}")
    print(f"neff_calibrated[0]: {neff_calibrated[0]}")
    print(f"saved -> {out_path}")


if __name__ == "__main__":
    main()
