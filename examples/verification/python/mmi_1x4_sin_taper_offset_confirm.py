"""Full-fidelity confirmation of the design candidate found by the coarse
(mmi_1x4_sin_taper_offset_length_sweep.py, job 11817364) and refined
(mmi_1x4_sin_taper_offset_length_refine.py, job 11839158) sweeps --
validation_markdown/mmi_1x4_sin300_validation.md step 9c, 2026-08-14.

Both prior passes used reduced N_WL (5, then 1) as a design-tool/localization
shortcut. They converged cleanly: offset=3.1 um did not move between a 0.3
um-coarse and 0.05 um-fine offset grid, and mm_length converged to a BROAD
plateau (39.5-40.3 um all giving T=0.968-0.971, <0.7% imbalance) rather than
a fragile single point. This script is the promised next step: the same
device, at production wavelength sampling (N_WL=11, the full C-band every
other real device script in this repo actually reports), pinned at the
confirmed offset and a narrow length window around the confirmed plateau.

offset is FIXED at 3.1 um (not swept) -- this is a confirmation, not another
design search, so there is only one geometry per branch: 4 tasks total, not
4 x N_offsets.

mm_lengths: 38.0-42.0 um at 0.1 um steps (41 points) -- brackets the
confirmed plateau with margin, still free per solve_length_sweep()'s own
cost model.

COST. N_WL 1->11 is an 11x increase in per-wavelength mode-solve count vs.
the refine pass's own measured 55-58 min/task; the coarse sweep's own
N_WL=1->5 step scaled to ~350-393 min (a ~1.35x overhead above pure linear
scaling, not exactly proportional) -- so this confirmation is estimated at
roughly 10-14 hrs/task, not measured for N_WL=11 specifically yet. Smoke
test first, same two-step discipline as every other real run this session.

Run: one branch per invocation.
    MMI1X4_OUTPUT_IDX=0 python .../mmi_1x4_sin_taper_offset_confirm.py

mmi_1x4_sin_taper_offset_confirm_merge.py combines all 4 branch files once
they've all completed.
"""
import os

import numpy as np

from beef import MMI, EMESolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification",
    "mmi_1x4_sin_taper_offset_confirm"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Global simulation parameters -- N_WL=11, the full production C-band
# (matches mmi_1x4_eme.py/mmi_1x4_hybrid.py). Env-var overridable for a
# smoke pass.
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI1X4_TOC_RES", 40.0))
WL_MIN, WL_MAX = 1.50, 1.60
N_WL = int(os.environ.get("MMI1X4_TOC_NWL", 11))
N_BODY = int(os.environ.get("MMI1X4_TOC_NBODY", 20))
N_TAPER = int(os.environ.get("MMI1X4_TOC_NTAPER", 10))
C_TAPER = int(os.environ.get("MMI1X4_TOC_CTAPER", 8))

# ---------------------------------------------------------------------------
# Geometry -- offset FIXED at the confirmed candidate (not swept -- see
# module docstring). mm_width unchanged from both prior sweeps.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
MM_WIDTH = 12.0
STRAIGHT_LENGTH = 0.5

OFFSET = 3.1000   # um -- confirmed candidate, unchanged across both prior sweeps

mm_lengths = np.arange(38.0, 42.001, 0.1)   # brackets the confirmed 39.5-40.3 plateau

BODY_REGION = 2   # index into the 5-region chain: wg_in | taper_in | body | taper_out | wg_out

k = int(os.environ["MMI1X4_OUTPUT_IDX"])
if k not in (0, 1, 2, 3):
    raise ValueError(f"MMI1X4_OUTPUT_IDX must be 0-3, got {k}")

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

mmi = MMI(
    n_inputs=1, n_outputs=4,
    mm_length=float(mm_lengths[len(mm_lengths) // 2]),   # reference length for cross-section sampling only
    mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH,
    taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH,
    waveguide_offset=OFFSET,
    waveguide_height=0.30,
)

outermost_half_extent = max(abs(y) for y in mmi.output_ys) + TAPER_WIDTH / 2
assert outermost_half_extent < MM_WIDTH / 2, (
    f"offset={OFFSET:.4f}: outermost taper face (half-extent "
    f"{outermost_half_extent:.4f} um) does not fit inside mm_width/2="
    f"{MM_WIDTH / 2:.4f} um"
)

print(f"CONFIRMATION run: offset={OFFSET:.4f} um (fixed), branch={k}, "
      f"y={mmi.output_ys[k]:+.4f} um")
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
print(f"  branch={k} (y={mmi.output_ys[k]:+.4f} um) solved -- "
      f"band-mean T range [{T_k.mean(axis=1).min():.4f}, {T_k.mean(axis=1).max():.4f}]")

out_path = os.path.join(output_dir, f"confirm_branch{k}.npz")
np.savez(
    out_path,
    mm_lengths=mm_lengths, wavelengths=wavelengths, T_k=T_k,
    offset=OFFSET, branch=k, output_y=mmi.output_ys[k],
    waveguide_width=WAVEGUIDE_WIDTH, taper_width=TAPER_WIDTH,
    taper_length=TAPER_LENGTH, mm_width=MM_WIDTH,
    n_body=N_BODY, n_taper=N_TAPER, c_taper=C_TAPER,
)
print(f"Saved -> {out_path}")
print("Run mmi_1x4_sin_taper_offset_confirm_merge.py once all 4 branch files exist.")
