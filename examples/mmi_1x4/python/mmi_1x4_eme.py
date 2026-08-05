"""Tutorial: 1x4 MMI splitter on 300 nm Si3N4 -- EME simulation (BEEF method).

validation_markdown/mmi_1x4_sin300_validation.md step 6b.

Geometry: single-mode-ish 1.1 um input waveguide -> MMI body -> four 1.1 um
output waveguides. 300 nm Si3N4 core on SiO2 (SIN_300).

Method: Eigenmode Expansion (EME) via meow, mmi_1x2_eme.py's pattern exactly
(solver.solve() once per output branch, output_idx restricting the EME
cross-sections to that branch only -- "so the mode is the isolated waveguide
mode rather than a supermode", per that script's own comment). 1x4 means
FOUR solves, not two -- no mirror-symmetry shortcut (this device's single
input side is structurally different from its four-branch output side).

mode_solver="tidy3d", not "meep" -- step 4 found compute_modes_meep silently
wrong for this platform's weakly-confined modes (see
mmi_1x4_sin_body_length_sweep.py's module docstring for the full finding).

N_BODY=20, N_TAPER=10, C_TAPER=8: from step 5's design work (job 11660226:
23 guided body modes at mm_width=12.0, top 20 kept; job 11660092: only 5
guided modes at the taper's widest cross-section 2.2 um, comfortably under
the N_TAPER=10 budget). See mmi_1x4_sin_body_length_sweep.py for the full
provenance and the C_TAPER=8 discretization reasoning.

MM_LENGTH = 206.0 um -- MEASURED 2026-08-04 (job array 11667114 +
mmi_1x4_sin_body_length_sweep_merge.py), user-selected over the pure
argmin(spread) winner (280.0 um) for a better footprint/throughput/balance
tradeoff -- see the constant's own comment below and validation doc step 5
for the full three-way comparison.

Run:
    python /path/to/examples/mmi_1x4/python/mmi_1x4_eme.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")

from beef import MMI, EMESolver, SIN_300, plot_eme_geometry, plot_performance

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_1x4", "mmi_1x4_eme"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # pts/um -- this codebase's standing mode-solve default
WL_MIN = 1.50       # um, C-band (step 4c's settled choice)
WL_MAX = 1.60
N_WL = 11

# ---------------------------------------------------------------------------
# 1. Geometry -- steps 4-5's design.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
WAVEGUIDE_OFFSET = 2.4
MM_WIDTH = 12.0
STRAIGHT_LENGTH = 0.5   # um -- short wg_in/wg_out stub, matches mmi_1x2_eme.py

# MEASURED 2026-08-04, job array 11667114 + mmi_1x4_sin_body_length_sweep_merge.py.
# mm_length=206.0 is on that sweep's integer grid exactly (linspace(20,300,281)
# -> spacing 1.0), asserted below per the coupler-bug lesson (a "pinned" value
# argmin'd onto a grid that didn't contain it exactly, differing by 0.0208 um
# between two examples for weeks). NOT the pure argmin(spread) winner
# (spread=0.0013 at 280.0 um) -- deliberately chosen instead for a better
# balance of footprint, throughput, and split precision (spread=0.0052, 3.7%
# relative imbalance, IL=-2.44 dB, vs. 280 um's 1.0% imbalance/-2.83 dB/6.2x
# the length, and 45 um's 8.2% imbalance/-2.20 dB/shortest option). See
# validation doc step 5 for the full three-way comparison this was picked
# from -- a genuine design tradeoff, not a default.
MM_LENGTH = 206.0
_MM_LENGTH_GRID = np.linspace(20.0, 300.0, 281)   # mmi_1x4_sin_body_length_sweep.py's own grid
_nearest = _MM_LENGTH_GRID[np.argmin(np.abs(_MM_LENGTH_GRID - MM_LENGTH))]
assert _nearest == MM_LENGTH, (
    f"MM_LENGTH={MM_LENGTH} is not an exact point on the measured sweep grid "
    f"(nearest grid value: {_nearest}) -- do not silently snap, this is exactly "
    f"the coupler bug (0.0208 um discrepancy shipped for weeks)."
)

mmi = MMI(
    n_inputs=1, n_outputs=4,
    mm_length=MM_LENGTH,
    mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH,
    taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30,
)
print(f"MMI total length: {mmi.total_length:.2f} um")
print(f"Output ports at y = {mmi.output_ys} um")

# ---------------------------------------------------------------------------
# 2. Configure EME solver -- 5 regions: wg_in | taper_in | MMI body |
#    taper_out | wg_out, tapers subdivided into real cells (C_TAPER each).
# ---------------------------------------------------------------------------
N_BODY = 20     # body modes -- job 11660226
N_TAPER = 10    # taper-cell modes -- job 11660092 (5 guided at width=2.2)
C_TAPER = 8     # cells per taper (cross-section staircase resolution)

solver = EMESolver(
    num_modes=[1,       N_TAPER,  N_BODY,   N_TAPER,  1],
    num_cells=[1,       C_TAPER,  1,        C_TAPER,  1],
    lengths=  [STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, STRAIGHT_LENGTH],
    oxide_width=MM_WIDTH + 8.0,
    oxide_height=SIN_300.sz,
    res=RESOLUTION,
    material_stack=SIN_300,
    mode_solver="tidy3d",   # step 4's fix -- see module docstring
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
plot_eme_geometry(
    mmi, solver=solver,
    region_labels=["wg in", "taper in", "MMI body", "taper out", "wg out"],
    title="1x4 MMI on 300 nm Si3N4 -- EME",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved -> {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep -- one solve() call per output branch (four solves,
#    no mirror-symmetry shortcut -- see module docstring).
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nSweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um\n")

smat = None
for k in range(4):
    smat_k = solver.solve(mmi, wavelengths, ports=[f"out{k}"], output_idx=k, save_dir=save_dir)
    if smat is None:
        smat = smat_k
    else:
        smat._data.update(smat_k._data)
    print(f"  output {k} (y={mmi.output_ys[k]:+.2f} um) solved")
smat.check_passive()   # sanity: total output power must not exceed 1 at any wavelength

# ---------------------------------------------------------------------------
# 5. Print summary
# ---------------------------------------------------------------------------
T = [smat.transmission(f"out{k}") for k in range(4)]
T_total = sum(T)
wl_center = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_center))

print(f"\nResults at lambda = {wl_center:.3f} um:")
for k in range(4):
    print(f"  Output {k} (y={mmi.output_ys[k]:+.2f} um)  "
          f"T = {T[k][idx]:.4f}  IL = {10*np.log10(max(T[k][idx],1e-20)):.2f} dB")
print(f"  Total  T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  (ideal 1x4 split: T_each=0.25, IL_total=0.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved -> {output_dir}/smatrix.pkl")

plot_performance(smat, output_dir=output_dir, title="1x4 MMI Splitter on 300 nm Si3N4 (EME)")
print(f"Performance plot saved -> {output_dir}/performance.png")
