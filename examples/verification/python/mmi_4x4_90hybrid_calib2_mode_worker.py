"""4x4 MMI 90-degree hybrid -- n_eff calibration at the CURRENT (optimised)
design point, one excited body mode per process.

validation_markdown/mmi_4x4_90_hybrid_validation.md follow-up, 2026-09-11.
The device's only prior calibration attempt (job 12078959, 2026-08-20) was
measured at the SUPERSEDED geometry (mm_width=6.0, taper_width=0.8,
taper_length=1.0, mm_length=78.9396) -- everything about the design has
moved since (mm_width=8.5, taper_width=1.30, taper_length=4.0,
mm_length=159.30, the campaign documented 2026-08-28..09-05), and the old
result ("calibration does not rescue this device": phase +0.150 deg worse,
CMRR +2.547 dBe better but still 16 dBe short) has never been re-checked at
the point actually being validated.

Reuses the EXACT proven architecture from mmi_4x4_90hybrid_calibration_
mode_worker.py (2026-08-20): running all 16 modes' FDTD phase-vs-length
solves sequentially in ONE process stalled on whichever mode landed 16th
(measured: mode 16 decayed cleanly to t=220.1 in total isolation, so the
bug is specific to many sequential mp.Simulation builds in one long-lived
process, not to any one mode's physics -- see that script's own docstring
for the full diagnostic trail, including the undamped-Lorentzian blow-up
this was mistaken for at first, 1.95e-11 -> 1.7e53 by t~2081 on job
12069773_16). One mode per process (this script, run via a 16-task SLURM
array) sidesteps it entirely. `beta_fit_from_phases()` fits each
(wavelength, mode) cell from only that mode's own phase-vs-length data and
seed_beta (checked directly against the library's own code, not assumed --
see the original worker's docstring), so per-mode isolation is exact, not
an approximation.

ONE CHANGE from the 2026-08-20 worker, deliberate: decay_tol here is 1e-4
(decay_method="energy"), not that script's 1e-11. 1e-11 was already found
MEASURED UNREACHABLE for this device at production settings (validation
doc, mmi_1x2_decay_tol_measurement.md) well before this calibration
re-attempt, and 1e-4 is what every campaign point (the optimisation that
produced the current design) was actually evaluated at -- using 1e-11 here
would calibrate against a decay criterion the device's own reported
numbers never use.

Correctness note (identical to the 2026-08-20 worker, restated because this
script does not import it): `section_mode_beta_details()`'s post-fit step
never mixes information across modes -- each mode's own beta fit depends
only on that mode's own diagonal phase-vs-length data and that mode's own
seed_beta.

This worker: computes the SAME seed/lengths section_mode_beta_details()
would (unpatched -- section_mode_neff() is cheap and already cache-safe),
then calls a runtime-patched copy of solve_section_multi_length() that
excites ONLY the requested band, and saves that mode's own phase-vs-length
data + seed beta to a per-mode .npz file for the aggregator to combine.

Run: python mmi_4x4_90hybrid_calib2_mode_worker.py <mode_index 1..16>
"""
import inspect
import textwrap
import types
import sys
import os
import numpy as np

from beef import MMI, FDTDSolver, SOI_220
import beef.solvers.fdtd as fdtd_mod

# Pinned to mmi_4x4_90hybrid_hybrid.py's own settled CURRENT design point
# exactly (campaign rounds 8-10, 2026-08-28..09-05) -- NOT the superseded
# mm_width=6.0 point mmi_4x4_90hybrid_offset_sweep.py's own `base` module
# still carries (deliberately left there for that script's own historical
# reproducibility; this script is self-contained instead of importing it).
WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 8.5
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 1.30
TAPER_LENGTH = 4.0
MM_LENGTH = 159.30   # the production body length -- calibration is a pure
                      # per-cross-section, length-independent correction
                      # (checked directly: solve_section_length_sweep's
                      # neff_override enters only the propagation phase
                      # exp(2j*pi*neff/wl*L), never the mode basis itself),
                      # so this choice only affects where inside the body
                      # the 3 calibration lengths are sampled from, not
                      # what gets measured. Using the real production
                      # length keeps that sampling well inside a region
                      # this device's own body solves already run at.
N_MODES = 16
BODY_STUB = 0.25
WL_MIN, WL_MAX, N_WL = 1.29, 1.33, 3
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

MODE = int(sys.argv[1])   # 1-indexed, matches meep's own eig_band convention
assert 1 <= MODE <= N_MODES, f"mode index {MODE} out of range [1,{N_MODES}]"

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))
OUT_DIR = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "calib2_per_mode"))
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
    mmi = MMI(
        n_inputs=4, n_outputs=4,
        mm_length=MM_LENGTH, mm_width=MM_WIDTH, waveguide_width=WAVEGUIDE_WIDTH,
        taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
        waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
        use_neighbor_aperture=True,
    )
    half_body = mmi.mm_length / 2
    x0, x1 = -half_body + BODY_STUB, +half_body - BODY_STUB
    mid = 0.5 * (x0 + x1)
    n_modes = N_MODES

    solver = FDTDSolver(
        res=40.0, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_220, decay_tol=1e-4,
        decay_method="energy", energy_decay_dt=20.0,
    )

    # Seed + lengths, computed EXACTLY as section_mode_beta_details() does
    # (unpatched -- section_mode_neff() is cheap/cached, not the stalling
    # part). Kept in lockstep with the library's own logic deliberately.
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
        save_dir=save_dir,
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
