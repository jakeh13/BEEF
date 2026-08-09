"""Field visualization for the 2x2 MMI EME solve, one wavelength/port pair.

Companion to plot_mmi1x4_eme_fields.py (see that script's own docstring for
the full rationale) -- beef.plot_eme_fields() exists but had zero call sites
anywhere in examples/tests before that script. Same pattern here, pinned
EXACTLY to mmi_2x2/python/mmi_2x2_eme.py's own construction (geometry,
EMESolver config, wavelength grid) -- see that script for the full
provenance of every constant below.

Port convention (mmi_2x2_eme.py's own solve() calls, always input_idx=1,
the top input -- "bar"/"cross" name which output_idx that top input reaches):
    ports=["bar"],   input_idx=1, output_idx=1   (top -> top)
    ports=["cross"], input_idx=1, output_idx=0   (top -> bottom)
PORT below selects one of these two.

Cost: 33 cross-sections (wg_in | 15x taper_in | body | 15x taper_out | wg_out,
num_cells=[1,15,1,15,1]), tidy3d eigensolves only, no time-stepping -- same
"minutes" cost class as mmi_1x4's own field-plot job (11736246, ~18 min for
19 cross-sections), scaled up somewhat for the extra taper cells.

Run:
    python .../plot_mmi2x2_eme_fields.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")

import meow as mw

from beef import MMI, EMESolver, SOI_220, plot_eme_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_2x2", "mmi_2x2_eme"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Pinned to mmi_2x2/python/mmi_2x2_eme.py -- see that script for the full
# provenance of every constant here.
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.530, 1.565, 11

# What THIS script adds on top of mmi_2x2_eme.py's own config.
WL_PLOT = 1.55       # C-band, close to the sweep centre -- the field snapshot wavelength
PORT = "bar"         # "bar" (input_idx=1, output_idx=1) or "cross" (input_idx=1, output_idx=0)
_IDX = {"bar": (1, 1), "cross": (1, 0)}
_INPUT_IDX, _OUTPUT_IDX = _IDX[PORT]

mmi = MMI(
    n_inputs=2, n_outputs=2,
    mm_length=28.0, mm_width=2.8, waveguide_width=0.45,
    taper_length=4.0, taper_width=1.3, waveguide_offset=1.5,
    waveguide_height=0.22,
)

solver = EMESolver(
    num_modes=[1, 7, 15, 7, 1],
    num_cells=[1, 15, 1, 15, 1],
    lengths=[0.5, 4.0, 28.0, 4.0, 0.5],
    oxide_width=5.0, oxide_height=1.5,
    res=RESOLUTION, material_stack=SOI_220, mode_solver="tidy3d",
)

# The FULL production wavelength array, not just WL_PLOT -- _build_cells()
# hands this to _extrusion_rules(), and matching mmi_2x2_eme.py's own call
# exactly keeps this a faithful reproduction of what solve() actually builds
# internally, not a reasoned-to-be-equivalent shortcut.
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

print(f"Building {sum(solver.num_cells)} cells for port {PORT!r} "
      f"(input_idx={_INPUT_IDX}, output_idx={_OUTPUT_IDX}, "
      f"mmi_2x2_eme.py's own geometry/solver config)...")
cells = solver._build_cells(mmi, side=PORT, wavelengths=wavelengths,
                             input_idx=_INPUT_IDX, output_idx=_OUTPUT_IDX)
modes_per_cell = solver._modes_per_cell()

print(f"Solving modes at lambda={WL_PLOT} um for all {len(cells)} cells "
      f"(modes/cell: {modes_per_cell})...")
env = mw.Environment(wl=float(WL_PLOT), T=solver.temperature)
css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]

tag = f"{PORT}_{WL_PLOT:.3f}um"
plot_eme_fields(cells, modes, output_dir, tag=tag, input_mode=1)
print(f"Field plots saved -> {output_dir}/fields/eme_*_{tag}.png")
