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
WL_MIN, WL_MAX = 1.50, 1.60
N_WL = int(os.environ.get("MMI1X4_HPREWARM_NWL", 11))
DECAY_TOL = float(os.environ.get("MMI1X4_HPREWARM_DECAY_TOL", 1e-11))
CALIBRATE_BODY = True

WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 3.1, 12.0
MM_LENGTH = 39.90
N_MODES = int(os.environ.get("MMI1X4_HPREWARM_NMODES", 12))   # matches mmi_1x4_hybrid.py's own N_MODES

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30, use_neighbor_aperture=True,
)

_fdtd_common = dict(
    res=RESOLUTION, sim_width=MM_WIDTH + 8.0, sim_height=SIN_300.sz,
    dpml=0.5, source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SIN_300, decay_tol=DECAY_TOL,
)
input_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=MM_WIDTH + 8.0, oxide_height=SIN_300.sz,
    res=RESOLUTION, material_stack=SIN_300, mode_solver="tidy3d",
)
output_solver = FDTDSolver(**_fdtd_common)

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
body_stub = 0.25

boundaries = [-half_total, -half_body + body_stub, +half_body - body_stub, +half_total]
interface_modes = [N_MODES, N_MODES]

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
      f"decay_tol={DECAY_TOL}", flush=True)
coeffs = hybrid.solve_length_sweep(
    mmi, wavelengths, body_lengths, sweep_section=BODY_SECTION,
    geometry_kwargs={2: {"output_idx": k}}, save_dir=save_dir,
)
print(f"Branch {k} done -- cache warmed at {save_dir}", flush=True)
