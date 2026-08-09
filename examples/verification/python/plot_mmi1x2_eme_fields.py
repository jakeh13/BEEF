"""Field visualization for the 1x2 MMI EME solve, one wavelength/branch.

Companion to plot_mmi1x4_eme_fields.py (see that script's own docstring for
the full rationale) -- beef.plot_eme_fields() exists but had zero call sites
anywhere in examples/tests before that script. Same pattern here, pinned
EXACTLY to mmi_1x2/python/mmi_1x2_eme.py's own construction (geometry,
EMESolver config, wavelength grid) -- see that script for the full
provenance of every constant below.

Port convention (mmi_1x2_eme.py's own solve() calls):
    ports=["out_top"], output_idx=1
    ports=["out_bot"], output_idx=0
OUTPUT_PORT below selects one of these two (side string AND output_idx must
match together, they are not a numeric pair the way mmi_1x4's out{k}/k are).

Cost: 19 cross-sections (wg_in | 8x taper_in | body | 8x taper_out | wg_out,
num_cells=[1,8,1,8,1]), tidy3d eigensolves only, no time-stepping -- same
"minutes" cost class as mmi_1x4's own field-plot job (11736246, ~18 min).

Run:
    python .../plot_mmi1x2_eme_fields.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")

import meow as mw

from beef import MMI, EMESolver, SOI_160, plot_eme_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_1x2", "mmi_1x2_eme"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Pinned to mmi_1x2/python/mmi_1x2_eme.py -- see that script for the full
# provenance of every constant here.
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.26, 1.36, 11

# What THIS script adds on top of mmi_1x2_eme.py's own config.
WL_PLOT = 1.31             # O-band centre -- the field snapshot wavelength
OUTPUT_PORT = "out_top"    # "out_top" (output_idx=1) or "out_bot" (output_idx=0)
_OUTPUT_IDX = {"out_top": 1, "out_bot": 0}[OUTPUT_PORT]

mmi = MMI(
    n_inputs=1, n_outputs=2,
    mm_length=4.5, mm_width=2.0, waveguide_width=0.35,
    taper_length=1.0, taper_width=0.8, waveguide_offset=1.0,
    waveguide_height=0.16,
)

solver = EMESolver(
    num_modes=[1, 5, 12, 5, 1],
    num_cells=[1, 8, 1, 8, 1],
    lengths=[0.5, 1.0, 4.5, 1.0, 0.5],
    oxide_width=5.0, oxide_height=1.5,
    res=RESOLUTION, material_stack=SOI_160, mode_solver="tidy3d",
)

# The FULL production wavelength array, not just WL_PLOT -- _build_cells()
# hands this to _extrusion_rules(), and matching mmi_1x2_eme.py's own call
# exactly keeps this a faithful reproduction of what solve() actually builds
# internally, not a reasoned-to-be-equivalent shortcut.
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

print(f"Building {sum(solver.num_cells)} cells for port {OUTPUT_PORT!r} "
      f"(output_idx={_OUTPUT_IDX}, mmi_1x2_eme.py's own geometry/solver config)...")
cells = solver._build_cells(mmi, side=OUTPUT_PORT, wavelengths=wavelengths,
                             output_idx=_OUTPUT_IDX)
modes_per_cell = solver._modes_per_cell()

print(f"Solving modes at lambda={WL_PLOT} um for all {len(cells)} cells "
      f"(modes/cell: {modes_per_cell})...")
env = mw.Environment(wl=float(WL_PLOT), T=solver.temperature)
css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]

tag = f"{OUTPUT_PORT}_{WL_PLOT:.3f}um"
plot_eme_fields(cells, modes, output_dir, tag=tag, input_mode=1)
print(f"Field plots saved -> {output_dir}/fields/eme_*_{tag}.png")
