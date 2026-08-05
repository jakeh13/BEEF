"""Tutorial: 1x4 MMI splitter on 300 nm Si3N4 -- FDTD/EME hybrid simulation
(BEEF method).

validation_markdown/mmi_1x4_sin300_validation.md step 6d.

HybridSolver chains three sections -- [input] -> [body] -> [output] -- each
independently FDTD or EME (mmi_1x2_hybrid.py's/mmi_2x2_hybrid.py's pattern:
INPUT_METHOD = OUTPUT_METHOD = "fdtd", BODY_METHOD = "eme"). The MMI body is
a uniform mm_width slab, so EMESolver.solve_section_length_sweep() (used
inside the hybrid's own body section) would still be free across a length
sweep if this script grew one -- not needed here, a single mm_length run.

Unlike mmi_1x2_hybrid.py/mmi_2x2_hybrid.py, this device has NO left-right
mirror symmetry to reuse: 2x2's identical port arrangement on both body
faces lets it solve the input section once per input arm and reuse those
same matrices (transposed) as the output section, never separately
simulating it. A 1x4 has one input port and four structurally different
output ports -- input section solved once (the only input), output section
solved four times (output_idx=0..3), sharing the input+body simulations via
one `cache` dict across all four hybrid.solve() calls (mmi_1x2_hybrid.py's
own pattern for its two output ports, extended to four here).

Aperture: use_neighbor_aperture=True on the MMI geometry (not
half_domain_aperture on the FDTD solver) -- mmi_1x2_hybrid.py's/
mmi_2x2_hybrid.py's half_domain_aperture only works because a shared y=0
axis cleanly separates their two ports into "near"/"far" halves. This
device's four ports (y = [-3.6, -1.2, +1.2, +3.6] um) don't split that way:
a y=0 cut still groups ports 0&1 (both negative-y) and 2&3 (both positive-y)
together. use_neighbor_aperture instead sizes each port's aperture to the
distance to its OWN nearest neighbor specifically -- see
beef/geometry/mmi.py's source_aperture_at() docstring. Crosstalk between the
four output waveguides was checked separately and found negligible
(T_cross=0.0002 at the actual gap over the taper+stub span -- validation doc
step 6b/6c/6d, job 11667138), so this is a numerical-isolation choice, not
compensation for a real physical coupling problem. The single input port has
no neighbor (n_inputs=1) and keeps the default full-width aperture, same as
mmi_1x2_hybrid.py's own input section.

mode_solver="tidy3d" for the EME body sub-solver -- NOT "meep", despite
mmi_1x2_hybrid.py/mmi_2x2_hybrid.py using "meep" to match FDTD's own
eigenmode basis (their own comment: "matches FDTDSolver's own eigenmode
basis in this hybrid chain"). Step 4 found compute_modes_meep silently wrong
for this platform's weakly-confined modes, and the body's near-cutoff bands
(job 11660226: margins down to 0.004) are exactly the regime that broke it
-- using it here would risk a wrong body mode basis, not just the ~2.3%/mode
phase mismatch calibrate_body_neff exists to fix. Any resulting basis
mismatch between the body's tidy3d modes and the FDTD sections' own meep
eigenbasis should be checked against calibrate_body_neff=True (available,
unused by default here -- same reasoning mmi_1x2_hybrid.py's own
CALIBRATE_BODY flag documents) if step 8's engine-agreement check finds it
matters.

MM_LENGTH = 206.0 um -- MEASURED 2026-08-04, see mmi_1x4_eme.py's module
docstring and the constant's own comment below for the full provenance.

Run:
    python /path/to/examples/mmi_1x4/python/mmi_1x4_hybrid.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, EMESolver, FDTDSolver, HybridSolver, SIN_300, plot_hybrid_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_1x4", "mmi_1x4_hybrid"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN = 1.50
WL_MAX = 1.60
N_WL = 11

# ---------------------------------------------------------------------------
# 0. Simulation method per section
# ---------------------------------------------------------------------------
INPUT_METHOD = "fdtd"
BODY_METHOD = "eme"
OUTPUT_METHOD = "fdtd"

# ---------------------------------------------------------------------------
# 1. Geometry -- steps 4-5's design. use_neighbor_aperture=True -- see
#    module docstring.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
WAVEGUIDE_OFFSET = 2.4
MM_WIDTH = 12.0

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
    use_neighbor_aperture=True,
)
print(f"MMI total length: {mmi.total_length:.2f} um")
print(f"Output ports at y = {mmi.output_ys} um")

# ---------------------------------------------------------------------------
# 2. Configure sub-solvers -- one per section, per the METHOD switches above.
# ---------------------------------------------------------------------------
_fdtd_common = dict(
    res=RESOLUTION,
    sim_width=MM_WIDTH + 8.0,
    sim_height=SIN_300.sz,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    material_stack=SIN_300,
    decay_tol=1e-11,
)

N_MODES = 20   # interface modes at input/output <-> body boundaries -- job
               # 11660226 (23 guided body modes at mm_width=12.0, top 20 kept)


def make_section_solver(method: str):
    """Build a solver for one hybrid section. use_neighbor_aperture lives on
    the MMI geometry (set above), not on the FDTD solver -- unlike
    mmi_1x2_hybrid.py/mmi_2x2_hybrid.py's half_domain_aperture, no per-call
    solver flag is needed here."""
    if method == "fdtd":
        return FDTDSolver(**_fdtd_common)
    if method == "eme":
        return EMESolver(
            num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
            oxide_width=MM_WIDTH + 8.0, oxide_height=SIN_300.sz,
            res=RESOLUTION, material_stack=SIN_300,
            mode_solver="tidy3d",   # step 4's fix -- see module docstring
        )
    raise ValueError(f"unknown method {method!r}, expected 'fdtd' or 'eme'")


input_solver = make_section_solver(INPUT_METHOD)
body_solver = make_section_solver(BODY_METHOD)
output_solver = make_section_solver(OUTPUT_METHOD)

# ---------------------------------------------------------------------------
# 3. Build HybridSolver
# ---------------------------------------------------------------------------
half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
body_stub = 0.25   # um of MMI body folded into each FDTD section, past the
                   # taper/body cross-section discontinuity -- inherited from
                   # mmi_1x2_hybrid.py/mmi_2x2_hybrid.py's convention, not
                   # re-derived for this platform (flagged, step 6d note)

boundaries = [
    -half_total,
    -half_body + body_stub,
    +half_body - body_stub,
    +half_total,
]
interface_modes = [N_MODES, N_MODES]

hybrid = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
)

# ---------------------------------------------------------------------------
# 4. Geometry plot
# ---------------------------------------------------------------------------
plot_hybrid_geometry(
    mmi, solver=hybrid,
    in_modes=1, out_modes=1,
    title="1x4 MMI on 300 nm Si3N4 -- Hybrid Solver Sections",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True,
    sub_region_show_eme_detail=True,
    sub_region_configs=[
        (f"out{k}", {2: {"output_idx": k}}) for k in range(4)
    ],
)
print(f"Geometry plot saved -> {output_dir}/geometry_hybrid.png")
print(f"Sub-region plots saved -> {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run hybrid simulation -- once per output port. output_idx is passed
#    ONLY to the output section's kwargs (input/body sections don't branch,
#    so their kwargs are identical across all four calls) and a shared cache
#    is passed to every hybrid.solve() call, so the input/body sections are
#    only ever simulated once, not once per output port.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nHybrid sweep: {wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

section_cache = {}
smats = []
for k in range(4):
    smat_k = hybrid.solve(
        mmi, wavelengths, port=f"out{k}", cache=section_cache, save_dir=save_dir,
        geometry_kwargs={2: {"output_idx": k}},
    )
    smats.append(smat_k)
    print(f"  output {k} (y={mmi.output_ys[k]:+.2f} um) solved")

smat = smats[0]
for smat_k in smats[1:]:
    smat._data.update(smat_k._data)
smat.check_passive()   # sanity: total output power must not exceed 1 at any wavelength

T = [smat.transmission(port=f"out{k}") for k in range(4)]
T_total = sum(T)

wl_c = wavelengths[len(wavelengths) // 2]
idx = int(np.argmin(np.abs(wavelengths - wl_c)))

print(f"\nResults at lambda = {wl_c:.3f} um:")
for k in range(4):
    print(f"  Output {k} (y={mmi.output_ys[k]:+.2f} um)  "
          f"T = {T[k][idx]:.4f}  IL = {10*np.log10(max(T[k][idx],1e-20)):.2f} dB")
print(f"  Total  T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  (ideal 1x4 split: T_each=0.25, IL_total=0.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix saved -> {output_dir}/smatrix_hybrid.pkl")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
colors = ["steelblue", "darkorange", "seagreen", "firebrick"]

ax = axes[0]
for k in range(4):
    ax.plot(wavelengths, T[k], label=f"Output {k} (y={mmi.output_ys[k]:+.2f} um)", color=colors[k])
ax.plot(wavelengths, T_total, label="Total", color="k", lw=1.5)
ax.axhline(0.25, color="gray", linestyle=":", lw=0.8, label="Ideal per-port")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Transmission")
ax.set_title("Hybrid -- Transmission vs Wavelength")
ax.legend(fontsize=7)
ax.grid(True)
ax.set_ylim(0, 1.05)

ax = axes[1]
for k in range(4):
    IL = 10 * np.log10(np.clip(T[k], 1e-20, None))
    ax.plot(wavelengths, IL, label=f"Output {k}", color=colors[k])
IL_tot = 10 * np.log10(np.clip(T_total, 1e-20, None))
ax.plot(wavelengths, IL_tot, label="Total", color="k", lw=1.5)
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title("Hybrid -- Insertion Loss vs Wavelength")
ax.legend(fontsize=7)
ax.grid(True)

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- FDTD/EME Hybrid Simulation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved -> {output_dir}/hybrid_results.png")
