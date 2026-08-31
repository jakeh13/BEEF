"""MMI body length (mm_length) design sweep for the 1x4 MMI on 300 nm Si3N4 --
validation_markdown/mmi_1x4_sin300_validation.md step 5.

SECOND VERSION -- the first version of this script (all-EME via
EMESolver.solve_section() over the taper+access bands) was WRONG, not just
imprecise: it returned T=1.0 for every output, every mm_length, every
wavelength (job 11660322 smoke test). Root cause, confirmed empirically
(job 11661894) and by reading meow's own source (create_cells/CrossSection.
from_cell/compute_s_matrix):

    solve_section() with a single Ls=[length] cell samples the geometry at
    ONE z location (Cell.structures_2d: "z = 0.5*(z_min+z_max)") and computes
    ONE set of modes for BOTH the left and right "ports" of that cell.
    mw.compute_s_matrix() with only one cell has no adjacent cross-section to
    build a mode-overlap interface against, so the returned block is the
    IDENTITY (mode i -> mode i, unit magnitude, phase only) -- BY
    CONSTRUCTION, regardless of how physically complex the sampled cross
    section is. Verified directly: solve_section() on this device's own
    output taper+access band, n_modes=6, returned exactly np.eye(6).

This means solve_section() can only ever be valid on a genuinely UNIFORM
region (no shape change end to end) -- exactly what
solve_section_length_sweep()'s docstring already requires for the MMI body,
and exactly why the existing hybrid examples (mmi_1x2_hybrid.py,
mmi_2x2_hybrid.py) use FDTD -- not EME solve_section() -- for every taper
band. There is no "all-EME, cheaper than hybrid" shortcut for a device with
real tapers; the established all-EME device examples (mmi_1x2_eme.py) don't
take one either: they call EMESolver.solve()/solve_length_sweep() (the FULL
multi-cell device solve, num_cells=[1,8,1,8,1] staircasing each taper into 8
real sub-cells with genuine mode-mixing interfaces between them), never
solve_section() directly on a taper span.

This version follows mmi_1x2_eme.py's actual pattern: one solver, 5 regions
(wg_in | taper_in | MMI body | taper_out | wg_out), tapers subdivided into
real cells, output_idx restricting the drawn branch per call (four calls,
1x4 has no mirror-symmetry shortcut to reuse -- see the original module
docstring's reasoning, still correct). The body length sweep rides
solve_length_sweep()'s multi-cell generalization of the single-cell trick
(sweep_region=2, the body): one mode solve per wavelength, every length
applied as an analytic phase -- so the wide mm_length sweep is still free.

mode_solver="tidy3d" throughout (not "meep", despite mmi_1x2_hybrid.py/
mmi_2x2_hybrid.py using "meep" to match FDTD's own basis) -- step 4 found
compute_modes_meep silently wrong for this platform's weakly-confined modes,
and the body's near-cutoff bands (job 11660226: margins down to 0.004) are
exactly the regime that broke it. All-EME here, so there is no FDTD basis to
match against anyway.

Geometry, from steps 4-5 (unchanged from v1):
    waveguide_width = 1.1 um   -- step 4's revised (2-mode) access width
    sz              = 5.0 um  -- step 4's converged MaterialStack.sz
    taper_width     = 2.2 um  -- 2.0x waveguide_width, same ratio order as
                                 mmi_1x2 (2.29x) / mmi_2x2 (2.89x)
    waveguide_offset= 2.4 um  -- taper_width + 0.2 um gap, the SAME
                                 taper-face clearance convention mmi_1x2
                                 (1.0-0.8=0.2) and mmi_2x2 (1.5-1.3=0.2) use
    mm_width        = 12.0 um -- geometric minimum for 4 ports at this
                                 offset/taper_width is (N-1)*offset+taper_width
                                 = 9.40 um; this adds ~28% clearance margin.
                                 Verified non-overlapping by direct plot
                                 (outermost port half-extent 4.70 um against
                                 mm_width/2=6.0 um, 1.30 um clearance).

N_BODY=20: job 11660226 found 23 guided modes at mm_width=12.0, margins above
n_clad from +0.228 (mode 0) down to +0.004 (mode 22, barely guided). The
bottom 3 (margins <0.03) are excluded -- close enough to cutoff to be poorly
converged, same reasoning mmi_2x2_hybrid.py's own N_MODES comment uses
("was 15: band 14 sits below the cladding index -- a radiation mode, dropped
at no cost").

N_TAPER=10, C_TAPER=8 (cells per taper): a PRELIMINARY choice for this
design-tool script, not yet convergence-checked -- the taper's widest
cross-section (2.2 um, at the body face) supports more modes than the 1.1 um
access guide's 2 but fewer than the body's 20; 10 is a working guess sized
between those two anchors. Step 6's own instruction already flags this
("Decide N_MODES and RESOLUTION jointly, not by copying -- neither device's
numbers transfer") -- convergence-check N_TAPER/C_TAPER there before this
script's numbers are treated as final.

Run:
    python /path/to/examples/verification/python/mmi_1x4_sin_body_length_sweep.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, EMESolver, plot_eme_geometry, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification", "mmi_1x4_sin_body_length_sweep"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_sim_cache"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # pts/um -- this codebase's standing mode-solve default
WL_MIN = 1.50       # um, C-band (step 4c's settled choice)
WL_MAX = 1.60
N_WL = 5            # reduced from 11 -- see cost note in module docstring

# ---------------------------------------------------------------------------
# 1. Geometry -- built at a REFERENCE mm_length; the body region's own
#    length is swept separately below via solve_length_sweep(sweep_region=2).
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
# ⚠ Constants refreshed 2026-08-31 to the post-step-9c geometry
# (mm_length 206.0 -> 39.90, waveguide_offset 2.4 -> 3.1, set
# 2026-08-14). Any output in example_output/ predating that date was
# produced with the OLD values and does not correspond to this file.
WAVEGUIDE_OFFSET = 3.1
MM_WIDTH = 12.0
STRAIGHT_LENGTH = 0.5   # um -- short wg_in/wg_out stub, matches mmi_1x2_eme.py

MM_LENGTH_REF = 150.0   # reference value for the geometry plot only

mmi = MMI(
    n_inputs=1, n_outputs=4,
    mm_length=MM_LENGTH_REF,
    mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH,
    taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30,
)
print(f"MMI total length (at reference mm_length={MM_LENGTH_REF}): {mmi.total_length:.2f} um")
print(f"Output ports at y = {mmi.output_ys} um")

# Geometric non-overlap check, asserted (not just eyeballed off a plot) --
# same check v1 used, unaffected by this rewrite.
_outermost_half_extent = max(abs(y) for y in mmi.output_ys) + TAPER_WIDTH / 2
assert _outermost_half_extent < MM_WIDTH / 2, (
    f"outermost output port's taper face (half-extent {_outermost_half_extent:.3f} um) "
    f"does not fit inside the MMI body (mm_width/2={MM_WIDTH/2:.3f} um) -- "
    "ports would overlap the body edge"
)
print(f"Port clearance check: outermost half-extent {_outermost_half_extent:.3f} um "
      f"< mm_width/2 {MM_WIDTH/2:.3f} um -- OK "
      f"(clearance {MM_WIDTH/2 - _outermost_half_extent:.3f} um)")

# ---------------------------------------------------------------------------
# 2. EME solver -- 5 regions: wg_in | taper_in | MMI body | taper_out | wg_out.
#    Region index 2 (body) is what solve_length_sweep(sweep_region=2) varies;
#    every other region's length is fixed at its self.lengths value below.
# ---------------------------------------------------------------------------
N_BODY = 20     # body modes -- see module docstring
N_TAPER = 10    # taper-cell modes -- see module docstring (preliminary)
C_TAPER = 8     # cells per taper (cross-section staircase resolution)

solver = EMESolver(
    num_modes=[1,       N_TAPER,  N_BODY,      N_TAPER,  1],
    num_cells=[1,       C_TAPER,  1,           C_TAPER,  1],
    lengths=  [STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH_REF, TAPER_LENGTH, STRAIGHT_LENGTH],
    oxide_width=MM_WIDTH + 8.0,   # generous lateral clearance past the body edges
    oxide_height=SIN_300.sz,      # 5.0 um, step 4's converged value
    res=RESOLUTION,
    material_stack=SIN_300,
    mode_solver="tidy3d",         # step 4's fix -- compute_modes_meep silently
                                  # returns wrong answers for this weakly-
                                  # confined material; see that step's writeup
)
BODY_REGION = 2

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
plot_eme_geometry(
    mmi, solver=solver,
    region_labels=["wg in", "taper in", "MMI body", "taper out", "wg out"],
    title="1x4 MMI on 300 nm Si3N4 -- Body Length Design Sweep",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved -> {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Body length sweep grid -- free (one mode solve per wavelength) per
#    solve_length_sweep()'s multi-cell generalization of the single-cell
#    trick, regardless of how many points this grid has.
# ---------------------------------------------------------------------------
mm_lengths = np.linspace(20.0, 300.0, 281)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nBody length sweep: {mm_lengths[0]:.1f}-{mm_lengths[-1]:.1f} um "
      f"({len(mm_lengths)} pts), wavelengths {wavelengths[0]:.3f}-{wavelengths[-1]:.3f} um "
      f"({N_WL} pts), {N_BODY} body modes, {N_TAPER} taper modes x {C_TAPER} cells\n")

# ---------------------------------------------------------------------------
# 5. Solve ONE output branch per process invocation, selected by the
#    MMI1X4_OUTPUT_IDX env var. This script used to loop over all four
#    branches in one process; the smoke test (job 11662670, reduced scope)
#    took ~15 min/branch even at N_BODY=6/N_TAPER=4/C_TAPER=2/N_WL=3 -- the
#    real config (N_BODY=20/N_TAPER=10/C_TAPER=8) scales to roughly a full
#    day if run sequentially in one job. Splitting into four independent
#    SLURM jobs (one per output_idx, run in parallel) turns that into
#    wall-clock hours instead of a day, at the same total CPU cost. Each
#    invocation restricts the EME cross-sections to its one output branch
#    only (mmi_1x2_eme.py's own per-branch convention, "so the mode is the
#    isolated waveguide mode rather than a supermode") -- no mirror-symmetry
#    shortcut available (this device's one input side is structurally
#    different from its four-branch output side).
#    mmi_1x4_sin_body_length_sweep_merge.py combines the four saved
#    per-branch files into the final mm_length pick + plot.
# ---------------------------------------------------------------------------
k = int(os.environ["MMI1X4_OUTPUT_IDX"])
if k not in (0, 1, 2, 3):
    raise ValueError(f"MMI1X4_OUTPUT_IDX must be 0-3, got {k}")

results_k = solver.solve_length_sweep(
    mmi, wavelengths, mm_lengths,
    ports=[f"out{k}"], sweep_region=BODY_REGION, output_idx=k,
    save_dir=save_dir,
)  # list[SMatrix], aligned with mm_lengths
T_k = np.array([smat.transmission(f"out{k}") for smat in results_k])  # [n_L, n_wl]
print(f"  output branch {k} (y={mmi.output_ys[k]:+.2f} um) solved")

np.savez(
    os.path.join(output_dir, f"mm_length_sweep_branch{k}.npz"),
    mm_lengths=mm_lengths, wavelengths=wavelengths, T_k=T_k, output_idx=k,
    output_y=mmi.output_ys[k],
    waveguide_width=WAVEGUIDE_WIDTH, taper_width=TAPER_WIDTH,
    taper_length=TAPER_LENGTH, waveguide_offset=WAVEGUIDE_OFFSET,
    mm_width=MM_WIDTH, n_body=N_BODY, n_taper=N_TAPER, c_taper=C_TAPER,
)
print(f"Branch {k} data saved -> {output_dir}/mm_length_sweep_branch{k}.npz")
print("Run mmi_1x4_sin_body_length_sweep_merge.py once all four branches complete.")
