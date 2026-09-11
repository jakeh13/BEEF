"""TEMPORARY -- prewarms mmi_1x4_hybrid.py's HybridSolver cache, one output
branch per process (env var MMI1X4_OUTPUT_IDX), mirroring mmi_1x4_hybrid.py's
own solver/geometry construction exactly (same _fdtd_common, N_MODES=12,
CALIBRATE_BODY=True, boundaries/body_stub, and the SAME mm_lengths/
body_lengths grid) so the resulting cache entries under save_dir are hit by
the real production run afterward, rather than silently missing on some
unrelated config difference.

Re-created 2026-08-14 for the NEW production geometry (MM_LENGTH=39.90,
WAVEGUIDE_OFFSET=3.1, validation doc step 9c) -- the original version of
this file (same name) prewarmed the OLD 206.0/2.4 geometry, succeeded (jobs
11696233/11734486, ~1.7 days/branch), and was deleted 2026-08-07 once that
run was confirmed. This is genuinely new territory beyond just a new
geometry, though: CALIBRATE_BODY=True (step 9a) has NEVER been run for real
at ANY geometry on this device -- it adds len(calibrate_lengths) (default
2) extra FDTD section solves on top of the input/output sections' own
~12-24 hr-class decay_tol=1e-11 cost (dc_validation.md item 6's mechanism,
same class of solve). So this prewarm is genuinely more expensive than the
original, not just a geometry swap -- real cost is UNMEASURED for this
combination, budget for it accordingly (see the batch script's own comment).

Redundancy, same as the original run: each of the 4 parallel tasks
redundantly repeats the SHARED input-section + calibration-section FDTD
solves (identical across branches, since nothing enforces an ordering
between the 4 racing array tasks) -- not fixed here, reproducing the same
already-proven structure rather than engineering a leaner 2-stage pipeline
under time pressure. Delete once mmi_1x4_hybrid.py's real run against the
warm cache is confirmed.

Run:
    MMI1X4_OUTPUT_IDX=0 python .../_prewarm_mmi_1x4_hybrid.py
"""
import os
import numpy as np

from beef import MMI, EMESolver, FDTDSolver, HybridSolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

RESOLUTION = float(os.environ.get("MMI1X4_HPREWARM_RES", 30.0))
# DEFAULTS UNCHANGED (1.50/1.60/11 -- the original full-band grid). Env-
# overridable as of 2026-09-11 for a C-band sub-sweep aimed at pinning the
# TE10 cutoff (the kink observed between the 16- and 24-mode results near
# 1.555-1.56 um -- see mmi_1x4_sin300_validation.md item 9d). Changing ANY of
# these changes `wavelengths`, which is hashed directly into every
# solve_section() fingerprint (_resolved_spectral_config, fdtd.py:479) -- a
# new grid is a guaranteed cache MISS on the input/output sections, not an
# incremental add. There is no way around that; only the target grid should
# be chosen deliberately, not spread across many ad hoc re-runs.
WL_MIN = float(os.environ.get("MMI1X4_HPREWARM_WL_MIN", 1.50))
WL_MAX = float(os.environ.get("MMI1X4_HPREWARM_WL_MAX", 1.60))
N_WL = int(os.environ.get("MMI1X4_HPREWARM_NWL", 11))
DECAY_TOL = float(os.environ.get("MMI1X4_HPREWARM_DECAY_TOL", 1e-11))
# DEFAULT UNCHANGED (True). Env-overridable as of 2026-09-06, matching
# mmi_1x4_hybrid.py's own MMI1X4_CALIBRATE_BODY knob, for the interface_modes
# convergence sweep (item 9d).
#
# WHY THE SWEEP TURNS IT OFF. The calibration is this prewarm's entire cost.
# It is measure_diagonal_only, which time-steps once per excited band, so it is
# 16 runs at N_MODES=16 against the input section's 1 -- measured on the
# cancelled array 12853535: ~4.7 of its projected ~5.3 days. Everything else
# (input section + output section) is under a day.
#
# Dropping it does not compromise what the sweep measures. A mode-count
# convergence test asks whether the S-matrix plateaus as interface_modes grows;
# the calibration replaces the body's per-mode beta and was measured to move
# total T by ~0.4% (0.9454 calibrated vs 0.9495 uncalibrated, jobs 12227591 /
# 12227682) against the ~3.45% of interface power that sits outside the 12-mode
# basis. It is also mildly HARMFUL against the whole-device FDTD reference at 8
# of 11 wavelengths. So it is not the effect under study, and it is an order of
# magnitude smaller than the one that is.
#
# ⚠ The resulting convergence conclusion is established for the UNCALIBRATED
# chain. If the plateau lands near the production 12, re-confirm it calibrated
# -- much cheaper once the plateau location is known than sweeping blind.
CALIBRATE_BODY = os.environ.get("MMI1X4_CALIBRATE_BODY", "1") not in ("0", "false", "False")

WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 3.1, 12.0
MM_LENGTH = 39.90
N_MODES = int(os.environ.get("MMI1X4_HPREWARM_NMODES", 12))   # matches mmi_1x4_hybrid.py's own N_MODES

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30, use_neighbor_aperture=True,
)

# FDTD_SIM_WIDTH/HEIGHT, CHANGED 2026-08-19 from MM_WIDTH+8.0/SIN_300.sz --
# see mmi_1x4_hybrid.py's own _fdtd_common comment for the full measurement
# (mmi_1x4_sin_body_domain_size_convergence.py, job 12029224). Must match
# mmi_1x4_hybrid.py's own values exactly -- this prewarms ITS cache.
FDTD_SIM_WIDTH = MM_WIDTH + 4.0
FDTD_SIM_HEIGHT = 4.0

_fdtd_common = dict(
    res=RESOLUTION, sim_width=FDTD_SIM_WIDTH, sim_height=FDTD_SIM_HEIGHT,
    dpml=0.5, source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SIN_300, decay_tol=DECAY_TOL,
)
input_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=FDTD_SIM_WIDTH, oxide_height=FDTD_SIM_HEIGHT,
    res=RESOLUTION, material_stack=SIN_300, mode_solver="tidy3d",
)
output_solver = FDTDSolver(**_fdtd_common)

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
body_stub = 0.25

boundaries = [-half_total, -half_body + body_stub, +half_body - body_stub, +half_total]
interface_modes = [N_MODES, N_MODES]

# Mirrors mmi_1x4_hybrid.py's PIN_MONITOR block (2026-09-05) -- read that
# comment for why both keys are needed and why the BODY section must not get
# them. This file's whole purpose is to warm THAT script's keys, so the two
# must agree exactly: prewarming unpinned entries while the real run asks for
# pinned ones warms a cache nothing will ever read.
PIN_MONITOR = os.environ.get("MMI1X4_PIN_MONITOR", "1") not in ("0", "false", "False")
_MONITOR_PIN = (dict(monitor_center=0.0, monitor_size=FDTD_SIM_WIDTH)
                if PIN_MONITOR else {})

hybrid = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries, interface_modes=interface_modes,
    calibrate_body_neff=CALIBRATE_BODY,
)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
mm_lengths = np.linspace(20.0, 300.0, 2801)
body_lengths = mm_lengths - 2 * body_stub
BODY_SECTION = 1

k = int(os.environ["MMI1X4_OUTPUT_IDX"])
if k not in (0, 1, 2, 3):
    raise ValueError(f"MMI1X4_OUTPUT_IDX must be 0-3, got {k}")

print(f"Prewarming branch {k}: HybridSolver.solve_length_sweep(), "
      f"{len(body_lengths)} lengths, {len(wavelengths)} wavelengths, "
      f"N_MODES={N_MODES}, CALIBRATE_BODY={CALIBRATE_BODY}, "
      f"PIN_MONITOR={PIN_MONITOR}, decay_tol={DECAY_TOL}", flush=True)
coeffs = hybrid.solve_length_sweep(
    mmi, wavelengths, body_lengths, sweep_section=BODY_SECTION,
    geometry_kwargs={0: dict(_MONITOR_PIN),
                     2: dict(_MONITOR_PIN, output_idx=k)},
    save_dir=save_dir,
)
print(f"Branch {k} done -- cache warmed at {save_dir}", flush=True)
