"""Local refinement sweep, following up on the coarse
mmi_1x4_sin_taper_offset_length_sweep.py result (validation_markdown/
mmi_1x4_sin300_validation.md step 9c, job 11817364, 2026-08-14): that 1 um/
0.3-ish um coarse grid found a sharp candidate at offset=3.1, mm_length=40 --
total T=0.9598, 0.77% relative imbalance, dominating every other design
point measured on this device (production: T=0.60/4.3%; the earlier 44 um
exploration: T=0.75/~15%). Neighbors on the coarse grid (39, 41 um) were
smooth and consistent, so this is a real peak, not a single-pixel artifact
-- but the grid itself is too coarse to trust 40.0/3.1 as the exact optimum.
This script narrows in.

WHAT CHANGED FROM THE COARSE SWEEP
-----------------------------------
  - offset grid: 7 points, [2.90, 2.95, ..., 3.20] um (0.05 um steps) --
    brackets the discovered peak (3.1) with margin on both sides, stops
    short of the coarse sweep's other (weaker) peak at 3.2167 -- that
    region is already reasonably characterized at 1 um length resolution
    and is the lower-throughput candidate, not the priority here.
  - mm_length grid: 201 points, linspace(30, 50) -- 0.1 um steps, vs the
    coarse sweep's 1 um. Free per solve_length_sweep()'s own cost model (one
    mode solve per wavelength, every length an analytic phase), so there is
    no reason to be stingy once a narrow range is picked.
  - N_WL=1 (center wavelength only), DOWN from the coarse sweep's N_WL=5.
    This is a LOCALIZATION pass, not a validation pass -- the goal is to
    find the true (offset, mm_length) peak cheaply, matching
    mmi_1x4_sin_taper_offset_sweep.py's own precedent for offset-focused
    exploratory sweeps (also N_WL=1 by default). A real N_WL=11
    confirmation at whatever point this pins down is the next step after
    this one, not folded in here.
  - N_BODY/N_TAPER/C_TAPER/res unchanged (20/10/8/40) -- these set the
    mode-solve basis fidelity, not what this pass is trying to save cost
    on; changing them would confound the comparison against the coarse
    sweep's own numbers.
  - Separate output directory (mmi_1x4_sin_taper_offset_length_refine, not
    the coarse sweep's own directory) -- avoids any collision with that
    already-complete, good 20-file dataset.

COST. N_WL 5->1 should cut per-task mode-solve cost roughly 5x versus the
coarse sweep's own measured 5:46-6:33 hrs/task (SLURM sacct, job 11817364) --
an ESTIMATE, not measured for this config yet, hence the smoke test before
the real 28-task (7 offsets x 4 branches) array.

Run: one (offset_idx, branch) pair per invocation, same convention as the
coarse sweep script.
    MMI1X4_TOL_OFFSET_IDX=0 MMI1X4_OUTPUT_IDX=0 \\
        python .../mmi_1x4_sin_taper_offset_length_refine.py

mmi_1x4_sin_taper_offset_length_refine_merge.py combines all 28 saved files
once every task has completed.
"""
import os

import numpy as np

from beef import MMI, EMESolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification",
    "mmi_1x4_sin_taper_offset_length_refine"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Global simulation parameters -- N_WL=1 (localization pass); everything
# else pinned to the coarse sweep's own config. Env-var overridable for a
# smoke pass.
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI1X4_TOL_RES", 40.0))
WL_MIN, WL_MAX = 1.50, 1.60
WL_CENTER = 0.5 * (WL_MIN + WL_MAX)
N_WL = int(os.environ.get("MMI1X4_TOL_NWL", 1))
N_BODY = int(os.environ.get("MMI1X4_TOL_NBODY", 20))
N_TAPER = int(os.environ.get("MMI1X4_TOL_NTAPER", 10))
C_TAPER = int(os.environ.get("MMI1X4_TOL_CTAPER", 8))

# ---------------------------------------------------------------------------
# Geometry -- fixed except waveguide_offset (swept). mm_width held at the
# same production value the coarse sweep used; see that script's own
# module docstring for why.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
MM_WIDTH = 12.0
STRAIGHT_LENGTH = 0.5

# Refined offset grid -- brackets the coarse sweep's discovered peak
# (offset=3.1) with margin. See module docstring for why 3.20, not 3.2167
# (the coarse sweep's own ceiling-adjacent second peak, lower priority here).
OFFSETS = np.arange(2.90, 3.2001, 0.05)

_ADJACENT_FLOOR = TAPER_WIDTH
_CONTAINMENT_CEIL = (MM_WIDTH - TAPER_WIDTH) / 3.0

assert np.all(OFFSETS >= _ADJACENT_FLOOR - 1e-9), (
    f"OFFSETS contains a value below the adjacent-taper floor {_ADJACENT_FLOOR} um"
)
assert np.all(OFFSETS < _CONTAINMENT_CEIL - 1e-6), (
    f"OFFSETS contains a value at/past the containment ceiling {_CONTAINMENT_CEIL:.6f} um"
)

# Refined mm_length grid -- brackets the coarse sweep's discovered peak
# (mm_length=40.0) at 0.1 um resolution, free per solve_length_sweep().
mm_lengths = np.linspace(30.0, 50.0, 201)

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
wavelengths = np.array([WL_CENTER]) if N_WL == 1 else np.linspace(WL_MIN, WL_MAX, N_WL)

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

out_path = os.path.join(output_dir, f"refine_offset{oi}_branch{k}.npz")
np.savez(
    out_path,
    mm_lengths=mm_lengths, wavelengths=wavelengths, T_k=T_k,
    offset_idx=oi, offset=offset, branch=k, output_y=mmi.output_ys[k],
    waveguide_width=WAVEGUIDE_WIDTH, taper_width=TAPER_WIDTH,
    taper_length=TAPER_LENGTH, mm_width=MM_WIDTH,
    n_body=N_BODY, n_taper=N_TAPER, c_taper=C_TAPER, offsets=OFFSETS,
)
print(f"Saved -> {out_path}")
print("Run mmi_1x4_sin_taper_offset_length_refine_merge.py once all "
      f"{len(OFFSETS)}x4={len(OFFSETS) * 4} (offset, branch) files exist.")
