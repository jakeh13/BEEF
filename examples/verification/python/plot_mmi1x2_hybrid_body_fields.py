"""Field visualization for the 1x2 MMI's HYBRID BODY SECTION (not pure EME).

Step 1 of validation_markdown/hybrid_field_extraction_validation.md.

WHY THIS EXISTS ALONGSIDE plot_mmi1x2_eme_fields.py. That script (2026-08-08)
already plots this device's fields, but it is pinned to mmi_1x2_eme.py's
PURE-EME configuration -- the whole five-region chain, tidy3d mode solver, 12
body modes, full 4.5 um body. The hybrid cascade propagates something
different, and the plan's own preamble disqualifies the pure-EME construction
as a stand-in for it ("would visualize the DISQUALIFIED method, not the
trusted hybrid-path the BO objective is actually built on" -- the 4x4
campaign measured pure-EME missing port crosstalk by 44-62 deg).

This script is pinned to mmi_1x2_hybrid.py's BODY SECTION instead:

    |                     | pure-EME script | THIS script (hybrid body) |
    | mode solver         | tidy3d          | meep                      |
    | body modes          | 12              | 8  (N_MODES)              |
    | body span           | 4.5 um (full)   | 4.0 um (+/- 2.0)          |
    | regions             | all 5           | body only                 |

mode_solver="meep" is not a preference -- it is what mmi_1x2_hybrid.py passes,
and compute_modes()'s own docstring gives the reason: it "matches the same
underlying eigenmode solver/conventions as FDTDSolver's own
get_eigenmode_coefficients() calls in HybridSolver chains, avoiding basis
mismatches between EME and FDTD sections."

WHAT THIS CAN AND CANNOT SHOW. The body is 62% of the device (4.0 of 6.50 um);
the input and output sections are FDTD and contribute no field data at all
today (zero add_dft_fields/get_array calls anywhere in fdtd.py -- that is
step 2 of the plan). So this is the body panel of an eventual three-panel
map, not the whole device.

It also cannot be read off the cascade's cache, and the reason is stronger
than "the cache stores arrays": with calibrate_body_neff=True (the shipped
default since 2026-08-17) EMESolver.solve_section() takes its neff_override
fast path, which SKIPS the eigenmode solve outright and builds the body
transfer matrix analytically as a diagonal exp(2j*pi*neff/wl*L). The shipped
hybrid run therefore never computes these mode fields at all. This script
pays for that one solve.

The flip side, and the reason step 3's body half is tractable: a diagonal
body transfer matrix means NO mode mixing across the body -- each interface
mode propagates independently, accumulating only phase.

Cost: one cross-section, 8 meep/MPB modes, one wavelength. The cheapest field
solve of any device in this repo.

Run:
    python .../plot_mmi1x2_hybrid_body_fields.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")

import meow as mw

from beef import MMI, EMESolver, SOI_160, plot_eme_fields
from beef._utils import extrude_polygons

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_1x2", "mmi_1x2_hybrid_body"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Pinned to mmi_1x2/python/mmi_1x2_hybrid.py -- every constant below is that
# script's own value. See it for the full provenance (N_MODES=8 comes from
# resolved/mmi_1x2_validation.md gap 3; RESOLUTION=40 from gap 2).
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.26, 1.36, 11
N_MODES = 8
BODY_STUB = 0.25

WL_PLOT = 1.31        # O-band centre -- the field snapshot wavelength
INPUT_MODE = 1        # 1-based; which body mode to excite at the left interface

mmi = MMI(
    n_inputs=1, n_outputs=2,
    mm_length=4.5, mm_width=2.0, waveguide_width=0.35,
    taper_length=1.0, taper_width=0.8, waveguide_offset=1.0,
    waveguide_height=0.16,
)

# The hybrid's body section: boundaries[1] .. boundaries[2].
half_body = mmi.mm_length / 2
x0 = -half_body + BODY_STUB
x1 = +half_body - BODY_STUB
# HybridSolver._stubs(k) gives a nonzero stub only to the FIRST (k=0) and LAST
# (k=K-1) sections; the body is k=1 of 3, so both are zero here.
stub_left = stub_right = 0.0
length = (x1 - x0) + stub_left + stub_right

# mmi_1x2_hybrid.py's body_solver, verbatim. num_modes/num_cells/lengths are
# placeholders on the solve_section() path (that method reads
# max(n_modes_in, n_modes_out) instead, and asserts len(num_cells)==1).
solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=5.0, oxide_height=1.5,
    res=RESOLUTION, material_stack=SOI_160,
    mode_solver="meep",
)

# The FULL production wavelength array, not just WL_PLOT -- _extrusion_rules()
# consumes it, and matching mmi_1x2_hybrid.py's own sweep keeps this a
# faithful reproduction of what the cascade builds rather than a
# reasoned-to-be-equivalent shortcut. (Same discipline as
# plot_mmi1x2_eme_fields.py's own note on this.)
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

print(f"1x2 hybrid BODY section: x in [{x0:+.2f}, {x1:+.2f}] um "
      f"(length {length:.2f} um, {N_MODES} modes, mode_solver='meep')")
print(f"device total_length = {mmi.total_length:.2f} um -- this section is "
      f"{100.0 * length / mmi.total_length:.0f}% of it")

# Replicates EMESolver.solve_section()'s own cell construction (eme.py's
# non-override branch) for this one section, rather than _build_cells(), which
# would build all five whole-device regions.
gds_polys = mmi.section_gds(x0, x1, solver.oxide_width,
                            stub_left=stub_left, stub_right=stub_right)
structs = extrude_polygons(gds_polys, solver._extrusion_rules(wavelengths))
mesh = solver._build_mesh()
cell_lengths = solver._flat_cell_lengths(region_lengths=[length])
modes_per_cell = [N_MODES] * len(cell_lengths)
cells = mw.create_cells(structures=structs, mesh=mesh, Ls=cell_lengths)

print(f"built {len(cells)} cell(s), lengths {cell_lengths}")
print(f"solving {N_MODES} modes at lambda={WL_PLOT} um...")

env = mw.Environment(wl=float(WL_PLOT), T=solver.temperature)
css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]

neff = [complex(m.neff).real for m in modes[0]]
print(f"body n_eff: {', '.join(f'{n:.6f}' for n in neff)}")
gaps = [neff[i] - neff[i + 1] for i in range(len(neff) - 1)]
print(f"smallest consecutive gap: {min(gaps):.6f}")

tag = f"hybrid_body_mode{INPUT_MODE}_{WL_PLOT:.3f}um"
plot_eme_fields(cells, modes, output_dir, tag=tag, input_mode=INPUT_MODE)
print(f"Field plots saved -> {output_dir}/fields/eme_*_{tag}.png")
