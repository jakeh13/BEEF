"""Joint taper-offset x mm_length design sweep for the 1x4 MMI on 300 nm
Si3N4, all-EME (mode_solver="tidy3d") -- validation_markdown/
mmi_1x4_sin300_validation.md step 9c.

Combines two sweeps that were previously run separately and never crossed:
  - mmi_1x4_sin_body_length_sweep.py swept mm_length at ONE fixed
    waveguide_offset (2.4 -- today's production value).
  - mmi_1x4_sin_taper_offset_sweep.py swept waveguide_offset at ONE fixed
    mm_length (44.0 -- the corrected hybrid sweep's own spread minimum),
    finding a promising but geometry-constrained direction (best reachable
    point 3.163 um, total T~0.762, short of the ~3.34 um the field's own
    peak positions implied).

Neither answered the real question: for EACH offset, what is the BEST
mm_length? This script does, by running the mm_length sweep (free --
solve_length_sweep()'s one-mode-solve-per-wavelength trick, same as
mmi_1x4_sin_body_length_sweep.py) at every one of a grid of offsets (NOT
free -- a fresh EMESolver.solve_length_sweep() call per offset, since offset
changes the cross-section shape at every one of the 19 cells the same way
mmi_1x4_sin_taper_offset_sweep.py's own module docstring already explains).

mm_width is held FIXED at today's production value (12.0 um) across the
whole sweep -- user's explicit choice, 2026-08-13, over widening it to fit
larger offsets. This keeps N_BODY's mode-count basis (20 of 23 guided modes
at mm_width=12.0, job 11660226) valid at every point without needing to
re-derive it for a wider body. Consequence: the offset grid is capped by
this device's own geometric feasibility ceiling at mm_width=12.0,
offset < (mm_width - taper_width)/3 = 3.2667 um (see
mmi_1x4_sin_taper_offset_sweep.py's own derivation).

OFFSETS -- user-specified, 2026-08-13, after two rounds of feasibility
checks against this ceiling and the ADJACENT-taper non-overlap floor
(offset > taper_width = 2.2 um, or neighboring tapers overlap at the body
face):
    [2.2, 2.5, 2.8, 3.1, 3.2167] um
2.2 sits exactly ON the taper_width floor (0.0 um edge-to-edge clearance --
touching, not overlapping; allowed, unlike a value below it). The 5th point
is 3.2167, NOT the literal 3.2667 ceiling the user typed in chat -- adjusted
down by the same 0.05 um safety margin mmi_1x4_sin_taper_offset_sweep.py's
own _FEASIBLE_MAX already uses, because 3.2667 (a 4-decimal rounding of the
true 3.2666...7 repeating fraction) rounds UP across the true boundary and
would fail this script's own outermost_half_extent < mm_width/2 assertion
by ~5e-5 um. Both endpoints asserted explicitly at import time below, not
just documented, per this project's own "pinned value must land exactly on
its own grid" lesson (the coupler's 0.0208 um discrepancy, dc_validation.md).

COST. Each (offset, branch) pair is a fresh, real EMESolver.solve_length_sweep()
call over the FULL mm_length grid -- same per-call cost class as
mmi_1x4_sin_body_length_sweep.py's own real run (~5:49-6:29 hrs/branch at
N_WL=5, job array 11667114). 5 offsets x 4 branches = 20 independent tasks,
no caching reuse ACROSS offsets (genuinely different cross-sections; within
one (offset, branch) task, save_dir still caches the per-cell mode solves
the way every other EME script here does). Budget ~6-7 hrs/task with margin.

N_WL=5 (not the production 11), N_BODY=20, N_TAPER=10, C_TAPER=8 -- same
design-tool config as mmi_1x4_sin_body_length_sweep.py's own real run, for
the same reason (this is an exploratory sweep to find a candidate design,
not the final validated result -- step 6's own real device scripts still own
that at N_WL=11).

Run: one (offset_idx, branch) pair per invocation, selected by env vars
MMI1X4_TOL_OFFSET_IDX (0-4) and MMI1X4_OUTPUT_IDX (0-3) -- SLURM array task
ID decoded into both by the batch script (task_id = 4*offset_idx + branch).
    MMI1X4_TOL_OFFSET_IDX=0 MMI1X4_OUTPUT_IDX=0 \\
        python .../mmi_1x4_sin_taper_offset_length_sweep.py

mmi_1x4_sin_taper_offset_length_sweep_merge.py combines all 20 saved
per-(offset,branch) files into the final (offset, mm_length) pick + plots,
once every task has completed.
"""
import os

import numpy as np

from beef import MMI, EMESolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification",
    "mmi_1x4_sin_taper_offset_length_sweep"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Global simulation parameters -- pinned to mmi_1x4_sin_body_length_sweep.py's
# own real-run config (step 5), not the production mmi_1x4_eme.py config
# (N_WL=11). Env-var overridable for a smoke pass.
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI1X4_TOL_RES", 40.0))
WL_MIN, WL_MAX = 1.50, 1.60
N_WL = int(os.environ.get("MMI1X4_TOL_NWL", 5))
N_BODY = int(os.environ.get("MMI1X4_TOL_NBODY", 20))
N_TAPER = int(os.environ.get("MMI1X4_TOL_NTAPER", 10))
C_TAPER = int(os.environ.get("MMI1X4_TOL_CTAPER", 8))

# ---------------------------------------------------------------------------
# Geometry -- fixed except waveguide_offset (swept). mm_width held at
# today's production value; see module docstring for why.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
MM_WIDTH = 12.0
STRAIGHT_LENGTH = 0.5

# User-specified offset grid, 2026-08-13 -- see module docstring for the
# feasibility derivation of both endpoints.
OFFSETS = np.array([2.2, 2.5, 2.8, 3.1, 3.2167])

_ADJACENT_FLOOR = TAPER_WIDTH   # offset < this overlaps neighboring tapers
_CONTAINMENT_CEIL = (MM_WIDTH - TAPER_WIDTH) / 3.0   # exact geometric ceiling, 3.2666...7

assert np.all(OFFSETS >= _ADJACENT_FLOOR - 1e-9), (
    f"OFFSETS contains a value below the adjacent-taper floor {_ADJACENT_FLOOR} um: "
    f"{OFFSETS[OFFSETS < _ADJACENT_FLOOR - 1e-9]}"
)
assert np.all(OFFSETS < _CONTAINMENT_CEIL - 1e-6), (
    f"OFFSETS contains a value at/past the containment ceiling {_CONTAINMENT_CEIL:.6f} um: "
    f"{OFFSETS[OFFSETS >= _CONTAINMENT_CEIL - 1e-6]}"
)

# mm_length sweep grid -- same range/density as mmi_1x4_sin_body_length_sweep.py's
# own (step 5), so results are directly comparable across both scripts.
mm_lengths = np.linspace(20.0, 300.0, 281)

BODY_REGION = 2   # index into the 5-region chain: wg_in | taper_in | body | taper_out | wg_out

# ---------------------------------------------------------------------------
# Task selection -- one (offset, branch) pair per process invocation.
# ---------------------------------------------------------------------------
oi = int(os.environ["MMI1X4_TOL_OFFSET_IDX"])
if not (0 <= oi < len(OFFSETS)):
    raise ValueError(f"MMI1X4_TOL_OFFSET_IDX must be 0-{len(OFFSETS) - 1}, got {oi}")
k = int(os.environ["MMI1X4_OUTPUT_IDX"])
if k not in (0, 1, 2, 3):
    raise ValueError(f"MMI1X4_OUTPUT_IDX must be 0-3, got {k}")

offset = float(OFFSETS[oi])
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

mmi = MMI(
    n_inputs=1, n_outputs=4,
    mm_length=float(mm_lengths[len(mm_lengths) // 2]),   # reference length for cross-section sampling only
    mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH,
    taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH,
    waveguide_offset=offset,
    waveguide_height=0.30,
)

outermost_half_extent = max(abs(y) for y in mmi.output_ys) + TAPER_WIDTH / 2
assert outermost_half_extent < MM_WIDTH / 2, (
    f"offset={offset:.4f}: outermost taper face (half-extent "
    f"{outermost_half_extent:.4f} um) does not fit inside mm_width/2="
    f"{MM_WIDTH / 2:.4f} um"
)
adjacent_gap = offset - TAPER_WIDTH
assert adjacent_gap >= -1e-9, (
    f"offset={offset:.4f}: adjacent tapers overlap by {-adjacent_gap:.4f} um "
    f"(taper_width={TAPER_WIDTH})"
)

print(f"offset_idx={oi} (offset={offset:.4f} um, adjacent gap={adjacent_gap:.4f} um, "
      f"outermost clearance={MM_WIDTH / 2 - outermost_half_extent:.4f} um), branch={k}")
print(f"Config: res={RESOLUTION}, N_WL={N_WL}, N_BODY={N_BODY}, N_TAPER={N_TAPER}, "
      f"C_TAPER={C_TAPER}, mm_length grid [{mm_lengths[0]:.1f}, {mm_lengths[-1]:.1f}] "
      f"({len(mm_lengths)} pts)")

solver = EMESolver(
    num_modes=[1, N_TAPER, N_BODY, N_TAPER, 1],
    num_cells=[1, C_TAPER, 1, C_TAPER, 1],
    lengths=[STRAIGHT_LENGTH, TAPER_LENGTH, float(mm_lengths[len(mm_lengths) // 2]),
             TAPER_LENGTH, STRAIGHT_LENGTH],
    oxide_width=MM_WIDTH + 8.0, oxide_height=SIN_300.sz,
    res=RESOLUTION, material_stack=SIN_300, mode_solver="tidy3d",
)

results_k = solver.solve_length_sweep(
    mmi, wavelengths, mm_lengths,
    ports=[f"out{k}"], sweep_region=BODY_REGION, output_idx=k,
    save_dir=save_dir,
)   # list[SMatrix], aligned with mm_lengths
T_k = np.array([smat.transmission(f"out{k}") for smat in results_k])   # [n_L, n_wl]
print(f"  offset_idx={oi} branch={k} (y={mmi.output_ys[k]:+.4f} um) solved -- "
      f"band-mean T range [{T_k.mean(axis=1).min():.4f}, {T_k.mean(axis=1).max():.4f}]")

out_path = os.path.join(output_dir, f"sweep_offset{oi}_branch{k}.npz")
np.savez(
    out_path,
    mm_lengths=mm_lengths, wavelengths=wavelengths, T_k=T_k,
    offset_idx=oi, offset=offset, branch=k, output_y=mmi.output_ys[k],
    waveguide_width=WAVEGUIDE_WIDTH, taper_width=TAPER_WIDTH,
    taper_length=TAPER_LENGTH, mm_width=MM_WIDTH,
    n_body=N_BODY, n_taper=N_TAPER, c_taper=C_TAPER, offsets=OFFSETS,
)
print(f"Saved -> {out_path}")
print("Run mmi_1x4_sin_taper_offset_length_sweep_merge.py once all "
      f"{len(OFFSETS)}x4={len(OFFSETS) * 4} (offset, branch) files exist.")
