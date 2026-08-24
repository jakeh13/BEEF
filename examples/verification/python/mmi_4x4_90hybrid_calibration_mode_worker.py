"""4x4 MMI 90-degree hybrid -- calibrate ONE excited mode, in its own
process (worker half of a one-process-per-mode calibration).

validation_markdown/mmi_4x4_90_hybrid_validation.md follow-up, 2026-08-20.
The full 16-mode calibration (job 12052385) got 15/16 modes decaying
cleanly under decay_method="energy" but stalled on whichever mode landed
16th in the sequential loop -- confirmed (mmi_4x4_90hybrid_mode16_
isolated.py, job 12066668) that mode 16 has NO intrinsic problem: it
decays perfectly normally (t=220.1) when it is the ONLY simulation in a
fresh process. The bug is specific to running 16 mp.Simulation builds/
resets sequentially in one long-lived process, not to any one mode's
physics. Fix: give every excited mode its OWN process (this script, run
once per mode via a SLURM array), instead of one process looping over
all 16.

Correctness note, checked against the library's own code, not assumed:
`section_mode_beta_details()`'s post-fit step
(`phases = np.angle(np.diagonal(T_multi, axis1=2, axis2=3))`, then one
`beta_fit_from_phases(lengths, phases[:, i, m], seed_beta[i, m])` call PER
MODE) never mixes information across modes -- each mode's own beta fit
depends only on that mode's own diagonal phase-vs-length data and that
mode's own seed_beta. Isolating modes into separate processes and merging
the results afterward (mmi_4x4_90hybrid_calibration_aggregate.py) is
therefore mathematically identical to the single-process run, not an
approximation.

This worker: computes the SAME seed/lengths section_mode_beta_details()
would (unpatched -- section_mode_neff() is cheap and already cache-safe,
not the part that was stalling), then calls a runtime-patched copy of
solve_section_multi_length() that excites ONLY the requested band
(same technique as mmi_4x4_90hybrid_mode16_isolated.py, generalized to
any mode index), and saves that mode's own phase-vs-length data + seed
beta to a per-mode .npz file for the aggregator to combine.

Run: python mmi_4x4_90hybrid_calibration_mode_worker.py <mode_index 1..16>
"""
import inspect
import textwrap
import types
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmi_4x4_90hybrid_offset_sweep as base
from beef import FDTDSolver, SOI_220
import beef.solvers.fdtd as fdtd_mod

MODE = int(sys.argv[1])   # 1-indexed, matches meep's own eig_band convention
assert 1 <= MODE <= base.N_MODES, f"mode index {MODE} out of range [1,{base.N_MODES}]"

OUT_DIR = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..",
    "example_output", "mmi_4x4_90hybrid", "calibration_per_mode"))
os.makedirs(OUT_DIR, exist_ok=True)


def make_patched_single_mode_method(k: int):
    """Runtime copy of FDTDSolver.solve_section_multi_length with the
    excite loop restricted to band k only -- see module docstring."""
    src = inspect.getsource(FDTDSolver.solve_section_multi_length)
    src = textwrap.dedent(src)
    assert "for k in range(1, n_excite + 1):" in src
    src = src.replace("for k in range(1, n_excite + 1):", f"for k in [{k}]:")
    src = src.replace(
        "def solve_section_multi_length(",
        f"def solve_section_multi_length_mode{k}_only(", 1,
    )
    ns = dict(fdtd_mod.__dict__)
    exec(compile(src, f"<mode{k}_isolated>", "exec"), ns)
    return ns[f"solve_section_multi_length_mode{k}_only"]


def main():
    mmi = base.make_mmi(base.BASELINE_OFFSET)
    boundaries, body_stub = base.boundaries_for(mmi)
    x0, x1 = boundaries[1], boundaries[2]
    mid = 0.5 * (x0 + x1)
    wavelengths = base.wavelengths
    n_modes = base.N_MODES

    solver = FDTDSolver(
        res=40.0, sim_width=base.MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_220, decay_tol=1e-11,
        decay_method="energy", energy_decay_dt=20.0,
    )

    # Seed + lengths, computed EXACTLY as section_mode_beta_details() does
    # (unpatched -- section_mode_neff() is cheap/cached, not the stalling
    # part). Kept in lockstep with the library's own logic deliberately,
    # not re-derived independently, so a future change there does not
    # silently diverge from what this worker assumes.
    span = float(x1 - x0)
    seed_span = min(2.0, span - 2.0 * solver.dpml)
    seed_anchors = 3
    wl_seed = wavelengths
    if (seed_anchors is not None and seed_anchors >= 2
            and len(wavelengths) > seed_anchors
            and float(np.ptp(wavelengths)) > 0.0):
        wl_seed = np.linspace(float(np.min(wavelengths)), float(np.max(wavelengths)), seed_anchors)

    seed_at = solver.section_mode_neff(
        mmi, wl_seed, mid - seed_span / 2.0, mid + seed_span / 2.0, n_modes,
        monitor_center=0.0, monitor_size=solver.sim_width,
        save_dir=base.save_dir,
    )
    if wl_seed is wavelengths:
        seed = seed_at
    else:
        x_mid = float(np.mean(wl_seed))
        seed = np.empty((len(wavelengths), n_modes))
        for m in range(n_modes):
            coeffs = np.polyfit(wl_seed - x_mid, seed_at[:, m], len(wl_seed) - 1)
            seed[:, m] = np.polyval(coeffs, wavelengths - x_mid)

    seed_beta = 2.0 * np.pi * seed / wavelengths[:, None]
    beta_max = float(np.max(np.abs(seed_beta)))
    from beef.solvers.fdtd import _BETA_MIN_BRANCH_FRACTION
    max_safe_spacing = 2.0 * np.pi / (_BETA_MIN_BRANCH_FRACTION * beta_max)
    usable_span = span - 2.0 * solver.dpml
    spacing = min(usable_span / 3.0, max_safe_spacing / 2.0)
    lengths = [spacing, 2.0 * spacing, 3.0 * spacing]

    print(f"mode {MODE}: seed_beta[:,mode-1]={seed_beta[:, MODE-1]}, lengths={lengths}")

    patched = make_patched_single_mode_method(MODE)
    solver.solve_section_multi_length_single = types.MethodType(patched, solver)

    T_multi = solver.solve_section_multi_length_single(
        mmi, wavelengths, mid, lengths, n_modes,
        measure_diagonal_only=True,
        monitor_center=0.0, monitor_size=solver.sim_width,
    )   # [n_lengths, n_wl, n_modes, n_modes], only row/col MODE-1 populated

    phases_mode = np.angle(T_multi[:, :, MODE - 1, MODE - 1])   # [n_lengths, n_wl]

    out_path = os.path.join(OUT_DIR, f"mode_{MODE:02d}.npz")
    np.savez(
        out_path, mode=MODE, lengths=np.array(lengths),
        phases=phases_mode, seed_beta=seed_beta[:, MODE - 1],
        wavelengths=wavelengths,
    )
    print(f"mode {MODE}: DONE, saved -> {out_path}")


if __name__ == "__main__":
    main()
