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

from beef import MMI, EMESolver, FDTDSolver, HybridSolver, SIN_300, SMatrix, plot_hybrid_geometry

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

# N_MODES=12, NOT the body's own 20-of-23-guided count (job 11660226) --
# a preliminary, unvalidated reduction, not a re-derived optimum. First
# attempt at N_MODES=20 (job array 11672739, 4 tasks, -t 24:00:00) hit the
# wall-time limit having only reached MPB band 14/20 on the FIRST of two
# interfaces (input<->body), with the second interface and all FDTD
# time-stepping still ahead -- true per-branch cost extrapolates to 60-80+
# hours, not the ~24h estimate that walltime was sized against. In
# HybridSolver, interface_modes IS the body's resolution for this cascade
# (drives solve_section_length_sweep()'s n_modes_in/out directly --
# EMESolver(num_modes=...) below is an unused placeholder for this call
# path, same as mmi_1x2_hybrid.py's/mmi_2x2_hybrid.py's own precedent
# comment already says), so there is no way to keep body accuracy at 20
# while only shrinking the FDTD interface count -- this genuinely trades
# some self-imaging accuracy for tractability, the same tradeoff
# mmi_1x2_hybrid.py (N_MODES=8) and mmi_2x2_hybrid.py (N_MODES=14) already
# accept relative to their own fully-guided counts. MPB's own log from the
# timed-out run showed near-cutoff modes taking far more iterations to
# converge than well-confined ones (up to 191 vs 8-10) -- dropping the
# marginal top modes should cut cost by more than the naive 12/20 ratio.
# Revisit at step 8's resolution/mode-count convergence check, alongside
# the still-open res=40 question that check already needs to answer.
N_MODES = 12


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
# 5. Run hybrid simulation across BOTH wavelength and mm_length in one pass.
#    The body is a uniform mm_width slab, so HybridSolver.solve_length_sweep()
#    solves it once per wavelength and applies every swept length as an
#    analytic propagation phase (beef/solvers/hybrid.py's own docstring) --
#    the whole 281-point length sweep costs the same as solving the single
#    MM_LENGTH=206.0 result alone would. The two FDTD sections dedupe against
#    any prior solve()/solve_length_sweep() call sharing this cache dict (same
#    config = same cache key regardless of which method built it), so this
#    costs no more than the single-length version did. output_idx is passed
#    ONLY to the output section's kwargs (input/body sections don't branch,
#    so their kwargs are identical across all four calls) and a shared cache
#    is passed to every call, so the input/body sections are only ever
#    simulated once, not once per output port -- same reuse the single-length
#    version already had, now extended across the length grid too.
#
#    MEASURED 2026-08-09 (verification/mmi_1x4_hybrid_body_length_sweep.py,
#    job 11770532): 25.3 s for all 4 branches x 281 lengths against the warm
#    cache -- confirmed cheap before folding into this script.
#
#    ⚠ body_lengths, NOT mm_lengths, is what gets swept -- MEASURED
#    2026-08-09/10: HybridSolver._stubs(k) returns (0.0, 0.0) for any
#    section that is neither the very first nor the very last (the body is
#    neither), so solve_section()'s own `length = (x1-x0) + stub_left +
#    stub_right` reduces to just `x1-x0` for the body -- which is
#    `mm_length - 2*body_stub`, NOT mm_length, since body_stub of MMI body is
#    folded into each adjacent FDTD section's own domain (module docstring's
#    "3. Build HybridSolver" comment) and boundaries[1]/boundaries[2] are
#    built accordingly. solve_section_length_sweep() applies whatever
#    `lengths` it is given directly as the propagation length with NO stub
#    adjustment (unlike solve_section()), so the caller must pre-subtract
#    2*body_stub to mean the same physical length solve_section() would at
#    that mm_length -- exactly what mmi_2x2_hybrid.py's own hand-rolled
#    cascade already does (`body_lengths = mm_lengths - 2 * body_stub`).
#    Missing this the first time through produced a real, measured ~6%
#    relative error on individual port T (0.1635 vs the correct 0.1541 at
#    output 0, λ=1.55, MM_LENGTH=206) -- large enough to matter because T(L)
#    oscillates rapidly with length here (see the sweep plot), so even a
#    0.5 µm propagation-length error lands on a materially different point
#    of the curve. Verified after the fix: solve(), solve_length_sweep()
#    with a single-element body_lengths, and this full 281-point sweep's own
#    slice at MM_LENGTH all agree exactly.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
mm_lengths = np.linspace(20.0, 300.0, 281)   # um; MM_LENGTH=206.0 lands exactly on this grid
body_lengths = mm_lengths - 2 * body_stub    # EME body section's own physical length
BODY_SECTION = 1   # index into solvers=[input(0), body(1), output(2)]
print(f"\nHybrid sweep: {wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um "
      f"({len(wavelengths)} pts, {N_MODES} interface modes), "
      f"mm_length {mm_lengths[0]:.1f} - {mm_lengths[-1]:.1f} um ({len(mm_lengths)} pts)\n")

section_cache = {}
coeffs_grids = np.zeros((4, len(mm_lengths), len(wavelengths)), dtype=complex)
for k in range(4):
    coeffs = hybrid.solve_length_sweep(
        mmi, wavelengths, body_lengths, sweep_section=BODY_SECTION,
        geometry_kwargs={2: {"output_idx": k}},
        cache=section_cache, save_dir=save_dir,
    )   # [n_L, n_wl, 1, 1] complex
    coeffs_grids[k] = coeffs[:, :, 0, 0]
    print(f"  output {k} (y={mmi.output_ys[k]:+.2f} um) solved")

T_grids = np.abs(coeffs_grids) ** 2   # [4, n_L, n_wl]

L_idx = int(np.argmin(np.abs(mm_lengths - MM_LENGTH)))
assert mm_lengths[L_idx] == MM_LENGTH, (
    f"MM_LENGTH={MM_LENGTH} is not an exact point on this sweep's own grid "
    f"(nearest: {mm_lengths[L_idx]}) -- do not silently snap, this is exactly "
    f"the coupler bug (0.0208 um discrepancy shipped for weeks)."
)

smat = SMatrix(wavelengths=wavelengths)
for k in range(4):
    smat.set(coeffs_grids[k, L_idx], port=f"out{k}")
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

# ---------------------------------------------------------------------------
# 7. Multimode-length sweep plot -- band-mean transmission per output +
#    max-min spread across the 4 outputs vs mm_length, same metric/format as
#    validation_markdown/mmi_1x4_sin300_validation.md step 5's own sweep
#    (mmi_1x4_sin_body_length_sweep.py) and this session's
#    mmi_1x4_hybrid_body_length_sweep.py verification script, now folded
#    into the real device script rather than living only in verification/.
#    A 1x4's natural sweep diagnostic is per-output balance, not a single
#    bar/cross ratio the way mmi_2x2_hybrid.py's own sweep plot uses -- that
#    device has exactly two outputs, this one has four.
# ---------------------------------------------------------------------------
T_band_mean = T_grids.mean(axis=2)   # [4, n_L]
T_total_grid = T_grids.sum(axis=0)   # [n_L, n_wl]
spread = T_band_mean.max(axis=0) - T_band_mean.min(axis=0)   # [n_L]
sweep_best_idx = int(np.argmin(spread))
mm_length_sweep_best = float(mm_lengths[sweep_best_idx])

print(f"\nLength-sweep optimum (min max-min spread across 4 outputs): "
      f"mm_length = {mm_length_sweep_best:.1f} um (spread={spread[sweep_best_idx]:.4f}); "
      f"chosen production value MM_LENGTH={MM_LENGTH} has spread={spread[L_idx]:.4f}")

np.savez(os.path.join(output_dir, "mm_length_sweep.npz"),
         mm_lengths=mm_lengths, wavelengths=wavelengths,
         T_grids=T_grids, T_total_grid=T_total_grid, spread=spread,
         mm_length_best=mm_length_sweep_best, mm_length_chosen=MM_LENGTH,
         output_ys=mmi.output_ys)
print(f"Length-sweep data saved -> {output_dir}/mm_length_sweep.npz")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

ax = axes[0]
for k in range(4):
    ax.plot(mm_lengths, T_band_mean[k],
             label=f"Output {k} (y={mmi.output_ys[k]:+.2f} um)", color=colors[k])
ax.plot(mm_lengths, T_total_grid.mean(axis=1), "k--", lw=1, label="Total")
ax.axvline(MM_LENGTH, color="purple", linestyle="-.", lw=1.2,
           label=f"chosen (production)={MM_LENGTH:.1f} um")
ax.axvline(mm_length_sweep_best, color="gray", linestyle=":", lw=1.2,
           label=f"sweep optimum={mm_length_sweep_best:.1f} um")
ax.axhline(0.25, color="gray", linestyle=":", lw=0.6)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("Band-mean transmission")
ax.set_title("Hybrid -- output power vs body length")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[1]
ax.semilogy(mm_lengths, spread, color="k")
ax.axvline(MM_LENGTH, color="purple", linestyle="-.", lw=1.2)
ax.axvline(mm_length_sweep_best, color="gray", linestyle=":", lw=1.2)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("max-min spread across 4 outputs")
ax.set_title("Hybrid -- equal-split quality vs body length")
ax.grid(True, which="both")

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- Hybrid body length sweep (analytic, no extra FDTD)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "mm_length_sweep.png"), dpi=300)
plt.close("all")
print(f"Length-sweep plot saved -> {output_dir}/mm_length_sweep.png")
