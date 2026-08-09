"""Field visualization for the 4x4 MMI 90-degree hybrid EME solve, one
wavelength/port pair.

Companion to plot_mmi1x4_eme_fields.py (see that script's own docstring for
the full rationale) -- beef.plot_eme_fields() exists but had zero call sites
anywhere in examples/tests before that script. Same pattern here, pinned
EXACTLY to mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_eme.py's own construction
(geometry, EMESolver config, wavelength grid) -- see that script for the full
provenance of every constant below, including the still-open A9/A11 first-
pass-design caveats it documents (N_MODES_BODY=8, taper params not
independently optimized, etc.) -- unchanged here, this script only adds a
field snapshot on top of an existing config, it does not re-derive it.

**Cost warning, inherited directly from mmi_4x4_90hybrid_eme.py's own
docstring**: that script measured ~416 s/mode for this device's OWN body
cross-section (Step 3a) and flags its 8 solve() calls as expensive enough to
need a smoke pass first via RESOLUTION/N_MODES_TAPER/N_CELLS_TAPER env-var
overrides -- unlike mmi_1x2/mmi_2x2/mmi_1x4, whose own field-plot companions
are "minutes" jobs. This script needs ALL 63 cross-sections of the 5-region
chain (wg_in | 30x taper_in | body | 30x taper_out | wg_out, num_cells=
[1,30,1,30,1]) solved once -- comparable cost to ONE of that script's own 8
solve() calls, not a new or larger cost than what that pipeline already
pays routinely, but still large in absolute terms. Reuses the SAME env-var
override names that script defines (MMI4X4_EME_RES / _NWL / _NMODES_TAPER /
_NCELLS_TAPER) for a smoke pass -- run one before the real-resolution field
plot, same two-step discipline.

Port convention (mmi_4x4_90hybrid_eme.py's own solve() calls):
    ports=[f"out{5+j}_in{i}"], input_idx=i, output_idx=j
    for i in DRIVEN=[0, 2] (paper's ports 1/3), j in OUTPUTS=[0,1,2,3]
INPUT_IDX/OUTPUT_IDX below select one such (driven input, measured output)
pair; default is the paper's port 1 -> port 5 (input_idx=0, output_idx=0).

Run (smoke pass first, recommended given the cost warning above):
    MMI4X4_EME_RES=20 MMI4X4_EME_NWL=1 MMI4X4_EME_NMODES_TAPER=2 \\
        MMI4X4_EME_NCELLS_TAPER=5 python .../plot_mmi4x4_90hybrid_eme_fields.py

Run (real config, once the smoke pass's timing is known):
    python .../plot_mmi4x4_90hybrid_eme_fields.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")

import meow as mw

from beef import MMI, EMESolver, SOI_220, plot_eme_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "mmi_4x4_90hybrid_eme"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Pinned to mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_eme.py -- see that
# script for the full provenance of every constant here, including the
# same env-var overrides for a smoke pass (see module docstring above).
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI4X4_EME_RES", 40.0))
WL_MIN, WL_MAX = 1.26, 1.36
N_WL = int(os.environ.get("MMI4X4_EME_NWL", 3))
N_MODES_TAPER = int(os.environ.get("MMI4X4_EME_NMODES_TAPER", 4))
N_CELLS_TAPER = int(os.environ.get("MMI4X4_EME_NCELLS_TAPER", 30))

WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 6.0
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 0.8
TAPER_LENGTH = 1.0
L_PI = 108.9195
DELTA_L = -1.50
MM_LENGTH = 0.75 * L_PI + DELTA_L   # = 80.1896
N_MODES_BODY = 8

# What THIS script adds on top of mmi_4x4_90hybrid_eme.py's own config.
WL_PLOT = 1.31       # O-band centre -- the field snapshot wavelength
INPUT_IDX = 0        # paper's port 1 (DRIVEN[0] in the real script)
OUTPUT_IDX = 0       # paper's port 5 / "out5" (OUTPUTS[0] in the real script)

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
    use_neighbor_aperture=True,
)

solver = EMESolver(
    num_modes=[1, N_MODES_TAPER, N_MODES_BODY, N_MODES_TAPER, 1],
    num_cells=[1, N_CELLS_TAPER, 1, N_CELLS_TAPER, 1],
    lengths=[0.5, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, 0.5],
    oxide_width=MM_WIDTH + 3.0, oxide_height=2.0,
    res=RESOLUTION, material_stack=SOI_220, mode_solver="tidy3d",
)

# The FULL production wavelength array, not just WL_PLOT -- _build_cells()
# hands this to _extrusion_rules(), matching mmi_4x4_90hybrid_eme.py's own
# call exactly keeps this a faithful reproduction of what solve() actually
# builds internally, not a reasoned-to-be-equivalent shortcut.
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

port_label = f"out{5 + OUTPUT_IDX}_in{INPUT_IDX}"
print(f"Building {sum(solver.num_cells)} cells for port {port_label!r} "
      f"(input_idx={INPUT_IDX}, output_idx={OUTPUT_IDX}, res={RESOLUTION}, "
      f"N_MODES_TAPER={N_MODES_TAPER}, N_CELLS_TAPER={N_CELLS_TAPER})...")
cells = solver._build_cells(mmi, side=port_label, wavelengths=wavelengths,
                             input_idx=INPUT_IDX, output_idx=OUTPUT_IDX)
modes_per_cell = solver._modes_per_cell()

print(f"Solving modes at lambda={WL_PLOT} um for all {len(cells)} cells "
      f"(modes/cell: {modes_per_cell})...")
env = mw.Environment(wl=float(WL_PLOT), T=solver.temperature)
css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]

tag = f"{port_label}_{WL_PLOT:.3f}um"
plot_eme_fields(cells, modes, output_dir, tag=tag, input_mode=1)
print(f"Field plots saved -> {output_dir}/fields/eme_*_{tag}.png")
