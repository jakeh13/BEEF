"""4x4 MMI 90-degree hybrid -- THE decisive test: does n_eff calibration fix
the hybrid-vs-FDTD phase drift at the reduced body lengths?

validation_markdown/mmi_4x4_90_hybrid_validation.md follow-up, 2026-09-11.

THE FINDING THIS TESTS. Comparing per-arm phase(S) between the hybrid
cascade and the full-device FDTD for the TOP port (paper port 4) at
mm_length = 5/10/15 um showed a systematic, wavelength-independent offset
that GROWS LINEARLY WITH LENGTH -- ~-55 deg at L=5, ~-118 deg at L=10,
~180 deg at L=15, and coherently across ALL FOUR output arms at once, not
just the dominant one. A linear-in-length phase offset is the signature of
a constant propagation-constant mismatch, Delta_phi = Delta_beta * L, and
back-solving the ~-13.6 deg/um slope gives Delta_n_eff ~ 0.05.

That prediction was then confirmed independently, from data already on
disk: the calibration (job 13098452/453) measures n_eff +0.0458 HIGHER
than the MPB/meow mode solver, uniformly across all 16 modes (spread only
0.00046 -- mode-INdependent to within 1%). Sign and magnitude both match
the phase-drift estimate.

So the hypothesis is: the uncalibrated hybrid propagates its body modes at
the mode solver's n_eff (2.9634 at 1.31 um, mode 1) while the full-device
FDTD propagates at what meep's own time-stepper actually realizes on its
Yee grid (3.0094) -- and the difference accumulates as phase over the body
length. Calibration exists precisely to measure the latter. If that is the
whole story, running the SAME cascade with neff_override=neff_calibrated
should collapse the phase offset toward zero at every length.

This script runs exactly that comparison, cheaply -- every input it needs
is already cached (M_in per port per length from the top-port probes; the
calibrated body grid from job 13098454; the uncalibrated body grid from
the top-port hybrid runs), so it is cascade arithmetic only, no new solves.

Run:
    python mmi_4x4_90hybrid_calib2_phase_check.py
"""
import os

import numpy as np

from beef import MMI, EMESolver, FDTDSolver, SMatrix, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "mmi_4x4_90hybrid_calib2_phase_check")
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.join(_ROOT, "example_sim_cache")
_CAL_DEFAULT = os.path.join(_ROOT, "example_output", "mmi_4x4_90hybrid",
                            "calib2_per_mode", "aggregated_beta_n16.npz")
# 2026-09-12: MMI4X4_CALPATH points this check at a DIFFERENT calibration (e.g.
# calib_probe/aggregated_res60cal_n16.npz, the res=60 family from job 13110658)
# to ask whether a better-converged calibration actually moves the measured
# hybrid-vs-FDTD phase agreement. Only the body's propagation constants change;
# M_in/M_out and the FDTD reference stay fixed, which is what isolates the
# calibration's own contribution.
CAL_PATH = os.environ.get("MMI4X4_CALPATH", _CAL_DEFAULT)
_CAL_IS_DEFAULT = os.path.abspath(CAL_PATH) == os.path.abspath(_CAL_DEFAULT)

WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 8.5
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 1.30
TAPER_LENGTH = 4.0
# 2026-09-13: MMI4X4_NMODES raises the interface basis for knob 3. Requires a
# matching calibration (MMI4X4_CALPATH) and pre-solved input sections at the
# same mode count -- mmi_4x4_90hybrid_input_presolve.py. The FDTD reference
# phases below are basis-independent (they are full-device measurements), so
# they stay the comparison target unchanged.
N_MODES = int(os.environ.get("MMI4X4_NMODES", 16))
BODY_STUB = 0.25

# True only for the exact configuration every cached grid on disk was built
# for. Gates cache_only, the full-length grids, and the untagged output name --
# so a probe at a different basis or calibration can neither fail on a cache
# miss it was always going to have, nor overwrite the production result.
_PRODUCTION = _CAL_IS_DEFAULT and N_MODES == 16
WL_MIN, WL_MAX, N_WL = 1.29, 1.33, 3
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

TOP_PORT_IDX = 3            # paper port 4 -- the port the FDTD reference drove
PORT_NAMES = ["out5", "out6", "out7", "out8"]
TARGET_LENGTHS = [5.00, 10.00, 15.00]
# MMI4X4_LENGTHS restricts the comparison to a subset. Needed for any non-16
# basis: the input sections must be pre-solved at that basis for EVERY length
# used (mmi_4x4_90hybrid_input_presolve.py), and M_in is fetched cache_only, so
# asking for a length that was not pre-solved fails rather than silently
# re-solving a multi-hour FDTD job on a login node.
if os.environ.get("MMI4X4_LENGTHS"):
    TARGET_LENGTHS = [float(x) for x in os.environ["MMI4X4_LENGTHS"].split(",")]

# Full-device FDTD reference phases (deg), port 4 driven, at 1.29/1.31/1.33.
# Measured: jobs 13069782 (L=5), 13069783 (L=10), 13069784 (L=15).
FDTD_PHASE = {
    5.00: dict(out5=[34.58649655, -131.17397169, 49.91945626],
               out6=[-58.70111646, 81.92947257, -45.01416398],
               out7=[49.93290295, -159.75841268, -49.65853321],
               out8=[160.71083898, -67.43651265, 62.62377765]),
    10.00: dict(out5=[-157.44224035, -53.76781916, 6.03966959],
                out6=[-39.0030333, -4.65247673, 39.09812191],
                out7=[-150.66049006, -109.12313724, -62.13633284],
                out8=[43.75948475, 89.39623992, 142.07973725]),
    15.00: dict(out5=[-62.17986379, 177.64162818, 122.16303988],
                out6=[52.48602574, -2.65407576, -42.48023759],
                out7=[67.77041603, 31.71407055, 6.97094881],
                out8=[-78.42020657, -120.39693215, -146.83156673]),
}


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
        material_stack=SOI_220, decay_tol=1e-4, decay_method="energy",
        energy_decay_dt=20.0,
    )
    input_solver = FDTDSolver(**fdtd_common)
    body_solver = EMESolver(
        num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
        oxide_width=fdtd_common["sim_width"], oxide_height=fdtd_common["sim_height"],
        res=40.0, material_stack=SOI_220, mode_solver="meep",
    )
    return input_solver, body_solver


def cascade_phases(M_in, M_out, M_body_slice):
    """Per-arm phase(S) in degrees for the TOP port, one body length."""
    out = {}
    tmp = np.einsum("Lwab,wbk->Lwak", M_body_slice, M_in[TOP_PORT_IDX])
    for j, name in enumerate(PORT_NAMES):
        c = np.einsum("wab,Lwbk->Lwak", M_out[j], tmp)
        out[name] = np.angle(c[0, :, 0, 0], deg=True)
    return out


def main():
    cal = np.load(CAL_PATH)
    neff_calibrated = cal["neff_calibrated"]
    print(f"calibration loaded: neff_calibrated[1.31] mode1={neff_calibrated[1,0]:.6f}")

    input_solver, body_solver = make_solvers()

    # Calibrated body, unified grid spanning every target length -- already
    # solved and cached by job 13098454, so this is a cache hit.
    mmi_cal = make_mmi(159.30)
    half_cal = mmi_cal.mm_length / 2
    # The production calibration's 5.00-400.00 grid is already solved and
    # cached (job 13098454), so reuse it verbatim -- the cache is keyed on the
    # requested lengths array, so asking for a subset would miss. A NON-default
    # calibration has no such grid, and only the three probe lengths are needed
    # here, so request exactly those rather than re-solving 7901 points.
    # _PRODUCTION is the 16-mode + default-calibration configuration that every
    # cached grid on disk was built for. Any departure from it (a different
    # calibration, or a different interface basis) has no cached grid at all,
    # so those runs request only the three probe lengths and are allowed to
    # solve them.
    cal_grid = (np.round(np.arange(5.00, 400.0 + 1e-9, 0.05), 2) if _PRODUCTION
                else np.array(TARGET_LENGTHS))
    M_body_cal_grid = body_solver.solve_section_length_sweep(
        mmi_cal, wavelengths, -half_cal + BODY_STUB, +half_cal - BODY_STUB,
        cal_grid - 2 * BODY_STUB, n_modes_in=N_MODES, n_modes_out=N_MODES,
        neff_override=neff_calibrated, save_dir=save_dir,
    )
    print(f"calibrated body grid: {M_body_cal_grid.shape}")

    rows = []
    for L in TARGET_LENGTHS:
        mmi = make_mmi(L)
        half_body = mmi.mm_length / 2
        half_total = mmi.total_length / 2
        boundaries = [-half_total, -half_body + BODY_STUB,
                      +half_body - BODY_STUB, +half_total]

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

        # Uncalibrated body: each length's own already-cached grid (production
        # basis only -- see _PRODUCTION above).
        grid_lo = 20.0 if L >= 20.0 else L
        uncal_grid = (np.round(np.arange(grid_lo, 400.0 + 1e-9, 0.05), 2)
                      if _PRODUCTION else np.array([L]))
        M_body_uncal_grid = body_solver.solve_section_length_sweep(
            mmi, wavelengths, boundaries[1], boundaries[2],
            uncal_grid - 2 * BODY_STUB, n_modes_in=N_MODES, n_modes_out=N_MODES,
            save_dir=save_dir, cache_only=_PRODUCTION,
        )
        iu = int(np.argmin(np.abs(uncal_grid - L)))
        ic = int(np.argmin(np.abs(cal_grid - L)))

        ph_uncal = cascade_phases(M_in, M_out, M_body_uncal_grid[iu:iu + 1])
        ph_cal = cascade_phases(M_in, M_out, M_body_cal_grid[ic:ic + 1])

        print(f"\n=== L={L:.2f} um -- phase(S) deg, top port (4) ===")
        print(f"{'arm':6s} {'FDTD':>28s} {'HYB uncal':>28s} {'HYB cal':>28s}"
              f" {'d(uncal)':>22s} {'d(cal)':>22s}")
        for name in PORT_NAMES:
            f = np.array(FDTD_PHASE[L][name])
            u = ph_uncal[name]
            c = ph_cal[name]
            du = ((u - f + 180) % 360) - 180
            dc = ((c - f + 180) % 360) - 180
            print(f"{name:6s} {np.round(f,1)!s:>28s} {np.round(u,1)!s:>28s}"
                  f" {np.round(c,1)!s:>28s} {np.round(du,1)!s:>22s} {np.round(dc,1)!s:>22s}")
            rows.append((L, name, f, u, c, du, dc))

        all_du = np.concatenate([((ph_uncal[n] - np.array(FDTD_PHASE[L][n]) + 180) % 360) - 180
                                 for n in PORT_NAMES])
        all_dc = np.concatenate([((ph_cal[n] - np.array(FDTD_PHASE[L][n]) + 180) % 360) - 180
                                 for n in PORT_NAMES])
        print(f"  SUMMARY L={L:.2f}: mean|d| uncal = {np.mean(np.abs(all_du)):7.2f} deg"
              f"   -> cal = {np.mean(np.abs(all_dc)):7.2f} deg")

    print("\n\n| L (um) | arm | mean |d phase| uncal | mean |d phase| cal |")
    print("|---|---|---|---|")
    for L in TARGET_LENGTHS:
        for name in PORT_NAMES:
            r = [x for x in rows if x[0] == L and x[1] == name][0]
            print(f"| {L:.2f} | {name} | {np.mean(np.abs(r[5])):.1f} | {np.mean(np.abs(r[6])):.1f} |")

    # Tagged when a non-default calibration is used, so a probe can never
    # overwrite the res=40 result the published figures are built from.
    _out_name = ("calib2_phase_check.npz" if _PRODUCTION else
                 "calib2_phase_check_n%d_" % N_MODES +
                 os.path.splitext(os.path.basename(CAL_PATH))[0] + ".npz")
    np.savez(os.path.join(output_dir, _out_name),
             lengths=np.array(TARGET_LENGTHS), wavelengths=wavelengths,
             **{f"L{L:.2f}_{n}_{k}": v for L, n, f, u, c, du, dc in rows
                for k, v in (("fdtd", f), ("uncal", u), ("cal", c),
                             ("d_uncal", du), ("d_cal", dc))})
    print(f"\nsaved -> {output_dir}/{_out_name}")


if __name__ == "__main__":
    main()
