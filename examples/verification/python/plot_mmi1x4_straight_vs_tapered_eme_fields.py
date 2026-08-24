"""Side-by-side EME field plots: the straight 1x4 body vs the tapered one.

`mmi_1x4_tapered_body_validation.md` -- the visual counterpart to step 3's
numeric evidence that a tapered body genuinely propagates differently.

WHY A NEW SCRIPT RATHER THAN RUNNING THE EXISTING ONE.
`plot_mmi1x4_eme_fields.py` (2026-08-07) is pinned to the SUPERSEDED
geometry -- MM_LENGTH=206.0, WAVEGUIDE_OFFSET=2.4, res=40 -- all three
replaced 2026-08-14 by the step-9c sweep (MM_LENGTH=39.90,
WAVEGUIDE_OFFSET=3.1) and 2026-08-18 (res 40 -> 30). It also has no saved
output anywhere in example_output/, so there is no straight-body field plot
of the CURRENT device to compare against. This script produces both, from
one config, so the only difference between the two runs is the body taper
itself -- which is the whole point of the comparison and is exactly what
running two separately-pinned scripts could not guarantee.

BOTH GEOMETRIES SHARE EVERY OTHER CONSTANT, pinned to
mmi_1x4/python/mmi_1x4_hybrid.py's production config (see that script for
each constant's provenance). The tapered variant differs ONLY by
body_taper_length/body_start_width.

THE TAPER PROFILE IS A PLACEHOLDER, NOT A DESIGN. body_start_width=6.0 /
body_taper_length=8.0 has had no sweep run against it -- step 6 picks the
real values. It is asymmetric (input face narrowed, output face left at
mm_width) because the 1x4's own output fan-out needs 11.5 of its 12.0 um
output face and leaves only 4.2% of room there, while the single on-axis
input port leaves 81.7% -- see that doc's blocker section.

THE BODY CELL COUNT IS WHAT MAKES THE TAPER VISIBLE. A tapered body drawn
as ONE cell would be sampled at a single cross-section and would look
straight; BODY_CELLS>1 is what resolves the width variation. Its value here
(24) is chosen for a legible picture, NOT from a convergence study -- step 5
is what establishes the converged count, and this script deliberately does
not pre-empt it. The straight body is drawn with the SAME cell count so the
two runs differ by geometry alone and not by discretization.

COST / SMOKE PASS. Every cell needs its own mode solve, so cost scales with
(cells x modes). Defaults here are already modest (MM_LENGTH=39.90 is 5.2x
shorter than the 206.0 the old script used), but run the smoke pass first --
same two-step discipline the 4x4's own field script documents:

    MMI1X4_FIELD_RES=15 MMI1X4_FIELD_NBODY=4 MMI1X4_FIELD_NTAPER=3 \\
        MMI1X4_FIELD_BODYCELLS=6 MMI1X4_FIELD_CTAPER=3 \\
        python .../plot_mmi1x4_straight_vs_tapered_eme_fields.py

    python .../plot_mmi1x4_straight_vs_tapered_eme_fields.py   # real config
"""
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")

import meow as mw

from beef import MMI, EMESolver, SIN_300, plot_eme_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_1x4",
    "mmi_1x4_straight_vs_tapered_eme"))
os.makedirs(output_dir, exist_ok=True)


def _env_f(name, default):
    return float(os.environ.get(name, default))


def _env_i(name, default):
    return int(os.environ.get(name, default))


# --- pinned to mmi_1x4_hybrid.py's production config ------------------------
RESOLUTION = _env_f("MMI1X4_FIELD_RES", 30.0)      # 40 -> 30, 2026-08-18
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 1.1, 2.2, 5.0
WAVEGUIDE_OFFSET, MM_WIDTH = 3.1, 12.0             # offset 2.4 -> 3.1, step 9c
MM_LENGTH = 39.90                                  # 206.0 -> 39.90, step 9c
STRAIGHT_LENGTH = 0.5

# --- what THIS script adds --------------------------------------------------
N_BODY = _env_i("MMI1X4_FIELD_NBODY", 12)          # = mmi_1x4_hybrid.py N_MODES
N_TAPER = _env_i("MMI1X4_FIELD_NTAPER", 10)
C_TAPER = _env_i("MMI1X4_FIELD_CTAPER", 8)
BODY_CELLS = _env_i("MMI1X4_FIELD_BODYCELLS", 24)  # legibility, not convergence
WL_PLOT = 1.55
OUTPUT_IDX = 0

# The taper profile under test -- placeholder, see module docstring.
BODY_TAPER_LENGTH = _env_f("MMI1X4_FIELD_TAPER_L", 8.0)
BODY_START_WIDTH = _env_f("MMI1X4_FIELD_START_W", 6.0)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

_geom_common = dict(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30,
)

CASES = [
    ("straight", dict()),
    ("tapered", dict(body_taper_length=BODY_TAPER_LENGTH,
                     body_start_width=BODY_START_WIDTH)),
]

print(f"res={RESOLUTION}  N_BODY={N_BODY}  N_TAPER={N_TAPER}  "
      f"C_TAPER={C_TAPER}  BODY_CELLS={BODY_CELLS}")
print(f"taper profile: start_width={BODY_START_WIDTH} -> mm_width={MM_WIDTH} "
      f"over Lt={BODY_TAPER_LENGTH} um (output face left at mm_width)")
print()

for name, extra in CASES:
    t0 = time.time()
    mmi = MMI(**_geom_common, **extra)

    # Sanity: the two cases must differ ONLY in the body profile.
    assert mmi.is_body_tapered == bool(extra), (
        f"{name}: is_body_tapered={mmi.is_body_tapered} contradicts kwargs {extra}"
    )
    half = mmi.mm_length / 2
    print(f"[{name}] body width: in-face={mmi.width_at(-half):.3f}  "
          f"mid={mmi.width_at(0.0):.3f}  out-face={mmi.width_at(half):.3f} um")

    solver = EMESolver(
        num_modes=[1, N_TAPER, N_BODY, N_TAPER, 1],
        # Body gets BODY_CELLS in BOTH cases -- see module docstring: equal
        # discretization is what makes this a geometry comparison.
        num_cells=[1, C_TAPER, BODY_CELLS, C_TAPER, 1],
        lengths=[STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH,
                 TAPER_LENGTH, STRAIGHT_LENGTH],
        oxide_width=MM_WIDTH + 4.0,     # 16.0 -- matches hybrid's FDTD_SIM_WIDTH
        oxide_height=4.0,               # matches hybrid's FDTD_SIM_HEIGHT
        res=RESOLUTION, material_stack=SIN_300, mode_solver="tidy3d",
    )

    n_cells = sum(solver.num_cells)
    print(f"[{name}] building {n_cells} cells (output branch {OUTPUT_IDX})...")
    cells = solver._build_cells(mmi, side=f"out{OUTPUT_IDX}",
                                wavelengths=wavelengths, output_idx=OUTPUT_IDX)
    modes_per_cell = solver._modes_per_cell()

    print(f"[{name}] solving modes at {WL_PLOT} um for {len(cells)} cells "
          f"({sum(modes_per_cell)} mode-solves total)...")
    env = mw.Environment(wl=float(WL_PLOT), T=solver.temperature)
    css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
    modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]

    tag = f"{name}_out{OUTPUT_IDX}_{WL_PLOT:.3f}um"
    plot_eme_fields(cells, modes, output_dir, tag=tag, input_mode=1)
    print(f"[{name}] done in {time.time() - t0:.1f}s -> "
          f"{output_dir}/fields/eme_*_{tag}.png")
    print()

print(f"Both cases written under {output_dir}/fields/")
# Component names are MEOW's, not the device frame's: meow's Ex is the
# device's transverse-horizontal field, i.e. the dominant TE component here.
# MEASURED on the smoke pass (job 12120999): |Ex| peaks at ~34, |Ey| at 0.63
# -- 54x weaker. eme_Ey_*.png is the sidewall-discontinuity component, which
# traces the body OUTLINE (a useful check that the taper is drawn where you
# think it is) but is not the propagating field.
print("Compare eme_Ex_straight_*.png against eme_Ex_tapered_*.png "
      "-- meow's Ex is the dominant TE component for this device "
      "(|Ex| ~ 34 vs |Ey| ~ 0.6, measured).")
print("eme_Ey_*.png traces the body outline instead -- handy for confirming "
      "the taper band sits where intended.")
