"""TEMPORARY -- prewarms mmi_1x4_eme.py's new solve_length_sweep()-based
cache entry (a genuinely different cache namespace from solve()'s own, since
the method name is part of the fingerprint -- see mmi_1x4_eme.py's own
Section 4 comment). One branch per process, selected by MMI1X4_OUTPUT_IDX,
mirroring the original single-length prewarm's own per-branch-array pattern
(job 11672738). Delete once mmi_1x4_eme.py's real run against the warm
cache has been confirmed.

Run:
    MMI1X4_OUTPUT_IDX=0 python .../_prewarm_mmi_1x4_eme_sweep.py
"""
import os
import numpy as np

from beef import MMI, EMESolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 2.4, 12.0
STRAIGHT_LENGTH = 0.5
MM_LENGTH = 206.0
N_BODY, N_TAPER, C_TAPER = 20, 10, 8
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
mm_lengths = np.linspace(20.0, 300.0, 281)

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
