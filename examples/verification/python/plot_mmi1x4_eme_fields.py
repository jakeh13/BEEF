"""Field visualization for the 1x4 MMI EME solve, one wavelength/branch.

validation_markdown/mmi_1x4_sin300_validation.md -- ad hoc verification, not
one of step 8's four checks. Answers "what does the field actually look like
propagating through this device" via beef.plot_eme_fields(), which already
exists in the codebase (beef/plot/fields.py, exported from beef's top level)
but -- checked, `grep -rn plot_eme_fields` across examples/ and tests/ --
has ZERO call sites anywhere as of 2026-08-07. mmi_1x4_eme.py's own
solver.solve() builds exactly the `cells`/`modes` this needs internally
(beef/solvers/eme.py:765-771) but discards them once it has the S-matrix, so
plotting fields means rebuilding that same pair for one wavelength, not
flipping a flag on solve() -- there isn't one.

Pinned EXACTLY to mmi_1x4/python/mmi_1x4_eme.py's own construction (geometry,
EMESolver config, wavelength grid) -- see that script for the full
provenance of every constant below; nothing here is re-derived. Only
difference: solves at ONE wavelength (WL_PLOT, C-band centre by default) and
ONE output branch (OUTPUT_IDX) instead of the full 11-point x 4-branch
sweep, since a field snapshot is inherently single-wavelength/single-branch.

Cost: 19 distinct cross-sections (wg_in | 8x taper_in | body | 8x taper_out |
wg_out -- the same cell structure mmi_1x4_eme.py's solver builds via
num_cells=[1,8,1,8,1]), one tidy3d eigensolve each at res=40
(mode_solver="tidy3d", step 4's fix, unchanged here) -- no time-stepping
anywhere, same "minutes, not hours" cost class as
audit_interface_modes.py's own per-device audits. Does NOT touch
save_dir/example_sim_cache -- solve()'s persistent cache only ever stores
S-matrices, not per-cell mode fields, so this is always a fresh solve
regardless of whether mmi_1x4_eme.py's own cache is warm.

Run:
    python .../plot_mmi1x4_eme_fields.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")

import meow as mw

from beef import MMI, EMESolver, SIN_300, plot_eme_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_1x4", "mmi_1x4_eme"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Pinned to mmi_1x4/python/mmi_1x4_eme.py -- see that script for the full
# provenance of every constant here.
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 2.4, 12.0
STRAIGHT_LENGTH = 0.5
MM_LENGTH = 206.0
N_BODY, N_TAPER, C_TAPER = 20, 10, 8

# What THIS script adds on top of mmi_1x4_eme.py's own config.
WL_PLOT = 1.55       # C-band centre -- the field snapshot wavelength
OUTPUT_IDX = 0       # which of the 4 output branches to plot (0..3)

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

# The FULL production wavelength array, not just WL_PLOT -- _build_cells()
# hands this to _extrusion_rules(), and matching mmi_1x4_eme.py's own call
# exactly (rather than e.g. [WL_PLOT] alone) keeps this a faithful
# reproduction of what solve() actually builds internally, not a
# reasoned-to-be-equivalent shortcut.
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

print(f"Building {sum(solver.num_cells)} cells for output branch {OUTPUT_IDX} "
      f"(mmi_1x4_eme.py's own geometry/solver config)...")
cells = solver._build_cells(mmi, side=f"out{OUTPUT_IDX}", wavelengths=wavelengths,
                             output_idx=OUTPUT_IDX)
modes_per_cell = solver._modes_per_cell()

print(f"Solving modes at lambda={WL_PLOT} um for all {len(cells)} cells "
      f"(modes/cell: {modes_per_cell})...")
env = mw.Environment(wl=float(WL_PLOT), T=solver.temperature)
css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]

tag = f"out{OUTPUT_IDX}_{WL_PLOT:.3f}um"
plot_eme_fields(cells, modes, output_dir, tag=tag, input_mode=1)
print(f"Field plots saved -> {output_dir}/fields/eme_*_{tag}.png")
