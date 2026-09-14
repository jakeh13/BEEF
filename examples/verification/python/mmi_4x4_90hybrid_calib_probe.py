"""4x4 MMI 90-degree hybrid -- parameterised single-mode calibration probe,
for root-causing the +0.0458 n_eff offset between the mode solver and what
meep's time-stepper actually propagates.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-09-11 ("ROOT
CAUSE FOUND" section). That offset fully explains the hybrid-vs-FDTD phase
drift at the reduced body lengths (predicted vs observed agree to 2.3% at
L=5/10/15). This script probes WHY the offset is there and what it depends
on, along three independent axes -- steps 3/4/5 of the debug plan:

  STEP 3 (MMI4X4_CALIB_RES) -- resolution. If the offset is Yee-grid
      numerical dispersion, it must scale ~1/res^2: res=20 should show ~4x
      the res=40 offset, res=60 should show ~(40/60)^2 = 0.44x of it. If it
      does NOT scale that way, grid dispersion is the wrong mechanism and
      the whole explanation needs revisiting. This is the falsification
      test for the leading hypothesis, so it is worth more than the other
      two combined.

  STEP 4 (MMI4X4_CALIB_WLMIN/WLMAX) -- source band and centre. The measured
      n_eff varies by 0.027 (mode 1) to 0.046 (mode 16) across 1.29-1.33 --
      the SAME ORDER as the offset being chased. Under dispersive=False
      this repo's _media() freezes the material index at the BAND CENTRE,
      so a band whose centre is not the monitored wavelength evaluates the
      material at the wrong place. Probing beta at a FIXED monitored
      wavelength (1.31) under bands that put 1.31 at the centre vs at
      either edge measures how badly that bites -- and therefore how
      tightly the full-device FDTD and the hybrid's section solves have to
      agree on their band, which is currently an unguarded assumption
      (they happen to match today at 1.29-1.33; nothing enforces it).

  STEP 5 (MMI4X4_CALIB_BWFACTOR) -- pulse bandwidth. Phase at a DFT monitor
      is set by the PHASE index, not the group index (n_g = 3.90 vs
      n_eff = 3.01 here -- a 0.89 gap, 20x the offset, so a simple
      phase/group mix-up is already ruled out on magnitude alone). The
      residual way n_g could contaminate is through the beta fit itself:
      a shorter (broader-band) pulse carries more envelope structure, and
      if that biases the extracted phase the fitted beta would move with
      pulse bandwidth. If beta at 1.31 is invariant across bandwidth
      factors, envelope/group effects are ruled out definitively.

Everything else is pinned to the CURRENT design point and to the proven
one-mode-per-process architecture (see mmi_4x4_90hybrid_calib2_mode_worker
.py's docstring for why per-mode isolation is exact and why it is needed).
Default mode is 1: the offset was measured mode-INDEPENDENT to under 1%
(+0.04558 to +0.04604 across all 16), so one mode answers all three
questions at 1/16 the cost.

Run:
    MMI4X4_CALIB_RES=20 MMI4X4_CALIB_TAG=res20 \
        python mmi_4x4_90hybrid_calib_probe.py [mode]
"""
import inspect
import textwrap
import types
import sys
import os
import numpy as np

from beef import MMI, FDTDSolver, SOI_220
import beef.solvers.fdtd as fdtd_mod

# --- design point: pinned to mmi_4x4_90hybrid_hybrid.py, unchanged --------
WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 8.5
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 1.30
TAPER_LENGTH = 4.0
MM_LENGTH = 159.30
# 2026-09-13: knob 3 of the weighted-phase plan needs a calibration at a LARGER
# interface basis than production's 16. The 16-mode basis captures only
# 99.2-99.3% of the top port's input power (measured from the cached M_in), and
# the missing ~0.7% is the right order to explain the ~2 deg residual. Raising
# this re-seeds the whole band solve, so a 24-mode calibration is NOT the
# 16-mode one plus 8 more -- every mode must be re-measured.
N_MODES = int(os.environ.get("MMI4X4_CALIB_NMODES", 16))
BODY_STUB = 0.25

# --- the three probe axes -------------------------------------------------
RES = float(os.environ.get("MMI4X4_CALIB_RES", 40.0))
WL_MIN = float(os.environ.get("MMI4X4_CALIB_WLMIN", 1.29))
WL_MAX = float(os.environ.get("MMI4X4_CALIB_WLMAX", 1.33))
N_WL = int(os.environ.get("MMI4X4_CALIB_NWL", 3))
BW_FACTOR = float(os.environ.get("MMI4X4_CALIB_BWFACTOR", 2.0))
TAG = os.environ.get("MMI4X4_CALIB_TAG", "baseline")

MODE = int(sys.argv[1]) if len(sys.argv) > 1 else 1
assert 1 <= MODE <= N_MODES, f"mode index {MODE} out of range [1,{N_MODES}]"

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))
OUT_DIR = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "calib_probe"))
os.makedirs(OUT_DIR, exist_ok=True)


def make_patched_single_mode_method(k: int):
    """Runtime copy of solve_section_multi_length restricted to band k --
    same technique, and same reason, as mmi_4x4_90hybrid_calib2_mode_worker
    .py: 16 sequential mp.Simulation builds in one process stalls, one
    mode per process does not."""
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
    print(f"[probe] tag={TAG} mode={MODE} res={RES} band={WL_MIN}-{WL_MAX} "
          f"({N_WL} pts) bw_factor={BW_FACTOR}")
    print(f"[probe] wavelengths={wavelengths}")

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

    solver = FDTDSolver(
        res=RES, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_220, decay_tol=1e-4,
        decay_method="energy", energy_decay_dt=20.0,
        source_bandwidth_factor=BW_FACTOR,
    )

    # Seed + lengths, computed exactly as section_mode_beta_details() does.
    span = float(x1 - x0)
    seed_span = min(2.0, span - 2.0 * solver.dpml)
    seed_anchors = 3
    wl_seed = wavelengths
    if (seed_anchors is not None and seed_anchors >= 2
            and len(wavelengths) > seed_anchors
            and float(np.ptp(wavelengths)) > 0.0):
        wl_seed = np.linspace(float(np.min(wavelengths)), float(np.max(wavelengths)), seed_anchors)

    seed_at = solver.section_mode_neff(
        mmi, wl_seed, mid - seed_span / 2.0, mid + seed_span / 2.0, N_MODES,
        monitor_center=0.0, monitor_size=solver.sim_width,
        save_dir=save_dir,
    )
    if wl_seed is wavelengths:
        seed = seed_at
    else:
        x_mid = float(np.mean(wl_seed))
        seed = np.empty((len(wavelengths), N_MODES))
        for m in range(N_MODES):
            coeffs = np.polyfit(wl_seed - x_mid, seed_at[:, m], len(wl_seed) - 1)
            seed[:, m] = np.polyval(coeffs, wavelengths - x_mid)

    seed_beta = 2.0 * np.pi * seed / wavelengths[:, None]
    beta_max = float(np.max(np.abs(seed_beta)))
    from beef.solvers.fdtd import _BETA_MIN_BRANCH_FRACTION
    max_safe_spacing = 2.0 * np.pi / (_BETA_MIN_BRANCH_FRACTION * beta_max)
    usable_span = span - 2.0 * solver.dpml
    spacing = min(usable_span / 3.0, max_safe_spacing / 2.0)
    lengths = [spacing, 2.0 * spacing, 3.0 * spacing]

    # The MODE-SOLVER n_eff at these same wavelengths, for the direct
    # comparison this whole probe exists to make.
    seed_mode = seed[:, MODE - 1]
    print(f"[probe] mode-solver n_eff (seed) = {seed_mode}")
    print(f"[probe] lengths = {lengths}")

    patched = make_patched_single_mode_method(MODE)
    solver.solve_section_multi_length_single = types.MethodType(patched, solver)

    T_multi = solver.solve_section_multi_length_single(
        mmi, wavelengths, mid, lengths, N_MODES,
        measure_diagonal_only=True,
        monitor_center=0.0, monitor_size=solver.sim_width,
    )
    phases_mode = np.angle(T_multi[:, :, MODE - 1, MODE - 1])

    # Fit beta here too, so each probe reports its own answer standalone
    # rather than needing the aggregator (single mode, no cross-mode data).
    from beef.solvers.fdtd import beta_fit_from_phases
    beta = np.zeros(len(wavelengths))
    resid = np.zeros(len(wavelengths))
    for i in range(len(wavelengths)):
        fit = beta_fit_from_phases(list(lengths), phases_mode[:, i], seed_beta[i, MODE - 1])
        beta[i] = fit.beta
        resid[i] = fit.resid
    neff_cal = beta * wavelengths / (2.0 * np.pi)
    offset = neff_cal - seed_mode

    print(f"\n[probe] RESULT tag={TAG} mode={MODE} res={RES} band={WL_MIN}-{WL_MAX} bw={BW_FACTOR}")
    for i, w in enumerate(wavelengths):
        print(f"[probe]   wl={w:.4f}  n_eff(modesolver)={seed_mode[i]:.6f}  "
              f"n_eff(calibrated)={neff_cal[i]:.6f}  offset={offset[i]:+.6f}  resid={resid[i]:.3e}")

    out_path = os.path.join(OUT_DIR, f"probe_{TAG}_mode{MODE:02d}.npz")
    np.savez(out_path, tag=TAG, mode=MODE, res=RES, n_modes=N_MODES,
             wl_min=WL_MIN, wl_max=WL_MAX,
             bw_factor=BW_FACTOR, wavelengths=wavelengths, lengths=np.array(lengths),
             phases=phases_mode, seed_beta=seed_beta[:, MODE - 1],
             neff_modesolver=seed_mode, neff_calibrated=neff_cal,
             offset=offset, resid=resid, beta=beta)
    print(f"[probe] saved -> {out_path}")


if __name__ == "__main__":
    main()
