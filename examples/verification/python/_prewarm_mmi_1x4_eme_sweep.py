"""TEMPORARY -- prewarms mmi_1x4_eme.py's solve_length_sweep()-based cache
entry (a genuinely different cache namespace from solve()'s own, since the
method name is part of the fingerprint -- see mmi_1x4_eme.py's own Section 4
comment). One branch per process, selected by MMI1X4_OUTPUT_IDX, mirroring
the original single-length prewarm's own per-branch-array pattern (job
11672738). Delete once mmi_1x4_eme.py's real run against the warm cache has
been confirmed.

UPDATED 2026-08-14 to the new production geometry (MM_LENGTH=39.90,
WAVEGUIDE_OFFSET=3.1, validation doc step 9c) -- supersedes this file's own
original purpose (prewarming the OLD MM_LENGTH=206.0/WAVEGUIDE_OFFSET=2.4
geometry's own N_WL=11 length-sweep plot, job array 11783204, which timed
out at 8h with zero progress markers and was left for "the next session to
pick up" without ever landing). That old geometry is no longer production;
prewarming it is no longer useful, so this file now targets the geometry
that actually matters.

mm_lengths NARROWED 2026-08-14 to [35, 45] um (101 points), NOT widened to
the full [20,300] range -- a real bug found the hard way, not a
cost-saving choice. A first version of this file DID widen to 2801 points
(matching mmi_1x4_eme.py's then-current full-range grid) on the assumption
that solve_length_sweep()'s "free length sweep" trick makes point count
irrelevant to cost. It is NOT irrelevant at this chain size: the real
prewarm attempt at 2801 points (job 11848462, N_WL=11) TIMED OUT at 10
hours with zero progress markers, because meow's compute_s_matrix() cascade
is called once per length across all 19 cells in this chain (mode SOLVES
are shared across lengths, the cascade APPLICATION is not) --
2801 lengths x 11 wavelengths = 30,811 calls, vs. the CONFIRMED design
point's own successful run at 41 lengths x 11 wavelengths = 451 calls
(mmi_1x4_sin_taper_offset_confirm.py, job 11840722, 4:17-4:41 hrs). See
mmi_1x4_eme.py's own module docstring for the full writeup -- this file's
mm_lengths must stay in lockstep with that script's own grid, since the
whole point is to warm its exact cache key.

Run:
    MMI1X4_OUTPUT_IDX=0 python .../_prewarm_mmi_1x4_eme_sweep.py
"""
import os
import numpy as np

from beef import MMI, EMESolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

RESOLUTION = float(os.environ.get("MMI1X4_PREWARM_RES", 40.0))
WL_MIN, WL_MAX = 1.50, 1.60
N_WL = int(os.environ.get("MMI1X4_PREWARM_NWL", 11))
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 3.1, 12.0
STRAIGHT_LENGTH = 0.5
MM_LENGTH = 39.90
N_BODY = int(os.environ.get("MMI1X4_PREWARM_NBODY", 20))
N_TAPER = int(os.environ.get("MMI1X4_PREWARM_NTAPER", 10))
C_TAPER = int(os.environ.get("MMI1X4_PREWARM_CTAPER", 8))
BODY_REGION = 2

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30,
)
solver = EMESolver(
    num_modes=[1, N_TAPER, N_BODY, N_TAPER, 1],
    num_cells=[1, C_TAPER, 1, C_TAPER, 1],
    lengths=[STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, STRAIGHT_LENGTH],
    oxide_width=MM_WIDTH + 8.0, oxide_height=SIN_300.sz,
    res=RESOLUTION, material_stack=SIN_300, mode_solver="tidy3d",
)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
mm_lengths = np.arange(35.0, 45.001, 0.1)   # must match mmi_1x4_eme.py's own grid exactly (same cache key)

k = int(os.environ["MMI1X4_OUTPUT_IDX"])
if k not in (0, 1, 2, 3):
    raise ValueError(f"MMI1X4_OUTPUT_IDX must be 0-3, got {k}")

print(f"Prewarming branch {k}: solve_length_sweep(), {len(mm_lengths)} lengths, "
      f"{len(wavelengths)} wavelengths, N_BODY={N_BODY}", flush=True)
results_k = solver.solve_length_sweep(
    mmi, wavelengths, mm_lengths, ports=[f"out{k}"],
    sweep_region=BODY_REGION, output_idx=k, save_dir=save_dir,
)
print(f"Branch {k} done -- cache warmed at {save_dir}", flush=True)
