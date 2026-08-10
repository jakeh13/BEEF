"""Output taper offset (center-to-center port spacing) design sweep for the
1x4 MMI on 300 nm Si3N4, at a SHORT mm_length=44.0 um -- the corrected
hybrid body-length sweep's own global spread-minimum candidate
(examples/verification/python/mmi_1x4_hybrid_body_length_sweep.py, job
11782082, AFTER the body_stub fix -- see that script's own module docstring
for the bug this corrects), essentially exactly where visual inspection of
the (pre-fix) plot flagged an early output-power crossing point.

WHY THIS EXISTS
----------------
mm_length=206 was tuned (step 5) against the CURRENT waveguide_offset=2.4 um
port spacing. The body-length sweep shows several OTHER candidate lengths
with comparable or better split balance, including one far shorter (44 um
vs 206 um -- a >4x footprint reduction on the multimode section alone). But
that sweep only ever tested the CURRENT port spacing at every length -- it
never asked whether 44 um's own self-imaging pattern (a physically DIFFERENT
image spacing than 206 um's, since these are different points on the beat
pattern) might actually want a DIFFERENT waveguide_offset to align cleanly
with where the output tapers sit. This sweep asks that question directly:
fix mm_length=44.0, vary waveguide_offset, see whether some offset gives a
band-mean split at least as good as the CURRENT 206 um/2.4 um combination
(step 5's own chosen point) -- MEASURED, not assumed either way.

COST -- genuinely different from the length sweep
----------------------------------------------------
The mm_length sweep is analytically free (EMESolver.solve_length_sweep()) --
one mode solve, every length applied as an propagation phase. Offset is NOT
free: it changes the taper/body CROSS-SECTION SHAPE at every one of the
5-region chain's 19 cells (the whole reason the taper is subdivided into
C_TAPER cells at all is to track that shape changing), so every swept
offset needs a genuinely fresh EMESolver.solve() call, same cost class as
mmi_1x4_eme.py's own real per-branch prewarm (~2:45 hrs/branch at
mm_length=206, job 11672738) -- COST DOES NOT OBVIOUSLY SHRINK just because
mm_length is shorter here, since that per-branch cost is dominated by
mode-solve count/complexity across the 19 cross-sections, not by body
length. RESOLUTION/N_WL/N_BODY/N_TAPER/C_TAPER are env-var overridable for a
smoke pass -- run one BEFORE the real sweep, same two-step discipline this
whole session's other genuinely-expensive scripts have used.

GEOMETRIC FEASIBILITY
----------------------
waveguide_offset must stay > TAPER_WIDTH (2.2 um, unchanged) for adjacent
output tapers not to overlap at the body face, and small enough that the
outermost port's taper face still fits inside mm_width=12 (unchanged) --
MMI's own uniform-pitch convention (single waveguide_offset parameter, not
custom per-port MMIPort placement) keeps center-to-center spacing identical
across all 4 ports by construction, matching the requirement that adjacent
taper spacing stay uniform. Outermost port y = offset*(4-1)/2 = 1.5*offset;
non-overlap+containment requires 1.5*offset + TAPER_WIDTH/2 < mm_width/2,
i.e. offset < (mm_width - TAPER_WIDTH) / 3 = 3.2667 um. Feasible range:
(2.2, 3.2667) um -- swept from just above the lower bound to just below the
upper one.

Run (smoke pass first, recommended given the cost warning above):
    MMI1X4_OFFSET_RES=20 MMI1X4_OFFSET_NWL=1 MMI1X4_OFFSET_NBODY=6 \\
    MMI1X4_OFFSET_NTAPER=4 MMI1X4_OFFSET_CTAPER=2 MMI1X4_OFFSET_NPOINTS=2 \\
        python .../mmi_1x4_sin_taper_offset_sweep.py <branch 0..3>

Run (real config, once smoke timing is known):
    python .../mmi_1x4_sin_taper_offset_sweep.py <branch 0..3>
"""
import os
import sys

import numpy as np

from beef import MMI, EMESolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification", "mmi_1x4_sin_taper_offset_sweep"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Pinned to mmi_1x4/python/mmi_1x4_eme.py's own config, except MM_LENGTH
# (fixed at the corrected hybrid sweep's short-length candidate, not the
# production 206.0) and WAVEGUIDE_OFFSET (swept, not fixed at 2.4).
# Env-var overrides for a smoke pass -- see module docstring.
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI1X4_OFFSET_RES", 40.0))
WL_MIN, WL_MAX = 1.50, 1.60
N_WL = int(os.environ.get("MMI1X4_OFFSET_NWL", 11))
N_BODY = int(os.environ.get("MMI1X4_OFFSET_NBODY", 20))
N_TAPER = int(os.environ.get("MMI1X4_OFFSET_NTAPER", 10))
C_TAPER = int(os.environ.get("MMI1X4_OFFSET_CTAPER", 8))
N_POINTS = int(os.environ.get("MMI1X4_OFFSET_NPOINTS", 21))

WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
MM_WIDTH = 12.0
STRAIGHT_LENGTH = 0.5

MM_LENGTH = 44.0   # um -- corrected hybrid sweep's own spread minimum (job 11782082);
                   # NOT the production 206.0 -- this script asks whether a matching
                   # port spacing exists for THIS short length, not that one.

# Feasible waveguide_offset range (see module docstring's geometric derivation).
_OFFSET_MIN = TAPER_WIDTH + 0.05   # 0.05 um minimum edge-to-edge margin, not zero
_OFFSET_MAX = (MM_WIDTH - TAPER_WIDTH) / 3.0 - 0.05
offsets = np.linspace(_OFFSET_MIN, _OFFSET_MAX, N_POINTS)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

if len(sys.argv) != 2 or sys.argv[1] not in ("0", "1", "2", "3"):
    raise SystemExit(f"usage: {sys.argv[0]} <output branch 0-3>")
BRANCH = int(sys.argv[1])

print(f"MM_LENGTH={MM_LENGTH} um (fixed), sweeping waveguide_offset over "
      f"[{offsets[0]:.4f}, {offsets[-1]:.4f}] um ({N_POINTS} pts), branch {BRANCH}")
print(f"Config: res={RESOLUTION}, N_WL={N_WL}, N_BODY={N_BODY}, N_TAPER={N_TAPER}, "
      f"C_TAPER={C_TAPER}")

T_k = np.zeros((N_POINTS, N_WL))
output_ys_k = np.zeros(N_POINTS)

for i, offset in enumerate(offsets):
    mmi = MMI(
        n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
        waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
        taper_width=TAPER_WIDTH, waveguide_offset=float(offset),
        waveguide_height=0.30,
    )
    outermost_half_extent = max(abs(y) for y in mmi.output_ys) + TAPER_WIDTH / 2
    assert outermost_half_extent < MM_WIDTH / 2, (
        f"offset={offset:.4f}: outermost taper face (half-extent "
        f"{outermost_half_extent:.3f} um) does not fit inside mm_width/2="
        f"{MM_WIDTH/2:.3f} um"
    )

    solver = EMESolver(
        num_modes=[1, N_TAPER, N_BODY, N_TAPER, 1],
        num_cells=[1, C_TAPER, 1, C_TAPER, 1],
        lengths=[STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, STRAIGHT_LENGTH],
        oxide_width=MM_WIDTH + 8.0, oxide_height=SIN_300.sz,
        res=RESOLUTION, material_stack=SIN_300, mode_solver="tidy3d",
    )
    smat = solver.solve(mmi, wavelengths, ports=[f"out{BRANCH}"], output_idx=BRANCH,
                         save_dir=save_dir)
    T_k[i] = smat.transmission(f"out{BRANCH}")
    output_ys_k[i] = mmi.output_ys[BRANCH]
    print(f"  offset={offset:.4f} um (output_y={mmi.output_ys[BRANCH]:+.4f} um): "
          f"T band-mean={T_k[i].mean():.4f} solved")

np.savez(
    os.path.join(output_dir, f"taper_offset_sweep_branch{BRANCH}.npz"),
    offsets=offsets, wavelengths=wavelengths, T_k=T_k, output_ys=output_ys_k,
    branch=BRANCH, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_width=TAPER_WIDTH,
    taper_length=TAPER_LENGTH, n_body=N_BODY, n_taper=N_TAPER, c_taper=C_TAPER,
)
print(f"Branch {BRANCH} data saved -> {output_dir}/taper_offset_sweep_branch{BRANCH}.npz")
