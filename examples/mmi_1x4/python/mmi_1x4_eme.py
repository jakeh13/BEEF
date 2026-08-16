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

MM_LENGTH = 39.90 um, WAVEGUIDE_OFFSET = 3.1 um -- CHANGED 2026-08-14,
validation doc step 9c. Supersedes the original MM_LENGTH=206.0/
WAVEGUIDE_OFFSET=2.4 pair (step 5, 2026-08-04): a user-requested
taper-offset x mm_length joint sweep (coarse grid, job 11817364; refined
grid, job 11839158) found this candidate, then CONFIRMED it at full
production N_WL=11 (job 11840722, this script's own real config): total
T=0.9620 (IL=-0.17 dB), 0.75% relative imbalance -- vs. the old point's
T=0.6014/4.26%. ~60% more throughput, ~5.2x shorter multimode section,
better balanced. mm_width=12.0 held fixed throughout that investigation
(user's explicit call, step 9c) -- not re-optimized, just not blocking this
change. See the constant's own comment below and validation doc step 9c for
the full candidate history (including the two weaker candidates measured
along the way: production's own 206.0/2.4, and an earlier 44.0 um/offset
~3.0-3.2 exploration at T~0.755/~15% imbalance).

mm_lengths window NARROWED 2026-08-14 to [35, 45] um (101 points) from the
original full-range [20, 300] (2801 points after the same-day grid-density
fix) -- NOT a cost-saving shortcut on the confirmed result, a real bug found
the hard way. solve_length_sweep() (the multi-cell, full-19-cell-chain
method this script uses) calls meow's compute_s_matrix() ONCE PER SWEPT
LENGTH inside the per-wavelength loop (beef/solvers/eme.py) -- the mode
SOLVE is genuinely shared across lengths as the docstring promises, but the
CASCADE (interface stitching + propagation across all 19 cells) is not, and
its own per-call cost is not negligible at this chain size. A real prewarm
attempt at the full 2801-point grid (job 11848462, N_WL=11, matching this
script's own config) TIMED OUT at 10 hours with zero progress markers
(tidy3d gives none) -- 30,811 total compute_s_matrix() calls (2801 lengths x
11 wavelengths) vs. the CONFIRMED candidate's own successful run
(mmi_1x4_sin_taper_offset_confirm.py, job 11840722, 41 lengths x 11
wavelengths = 451 calls, 4:17-4:41 hrs) is a ~68x jump in call count for the
same per-wavelength mode-solve cost. This is also consistent with, and now
explains, this script's OWN historical N_WL=11 length-sweep-plot attempt
(job 11783204, 281 lengths) timing out at 8h with the same symptom, which
was attributed only to N_WL scaling at the time -- the real driver was
ALSO this per-length cascade cost, compounding on top.

HybridSolver's own body-length sweep (mmi_1x4_hybrid.py) is NOT affected by
this: it calls solve_section_length_sweep() (the SINGLE-cell variant) for
the body region, where each compute_s_matrix() call is on one cell only --
genuinely cheap regardless of length count (MEASURED: 25.3s for all 4
branches x 281 lengths, job 11770532). The cost blowup here is specific to
the multi-cell chain this script's own EMESolver.solve_length_sweep() walks.

[35, 45] brackets the confirmed MM_LENGTH=39.90 with real margin (the
refine sweep's own plateau was 39.5-40.3) at a call count (101 x 11 = 1111)
close to the CONFIRMED-safe scale above, rather than re-deriving the full
design-space view this script used to show -- that context already lives in
mmi_1x4_sin_taper_offset_length_sweep.py's/..._refine.py's own plots
(step 9c), not duplicated here at prohibitive cost.

Run:
    python /path/to/examples/mmi_1x4/python/mmi_1x4_eme.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
WAVEGUIDE_OFFSET = 3.1   # um -- CHANGED 2026-08-14 from 2.4, see module docstring
MM_WIDTH = 12.0
STRAIGHT_LENGTH = 0.5   # um -- short wg_in/wg_out stub, matches mmi_1x2_eme.py

# CHANGED 2026-08-14 from 206.0 -- CONFIRMED at this script's own real
# config (N_WL=11, job 11840722), validation doc step 9c. mm_length=39.90 is
# on this file's own mm_lengths grid exactly (linspace(20,300,2801) ->
# spacing 0.1), asserted below per the coupler-bug lesson (a "pinned" value
# argmin'd onto a grid that didn't contain it exactly, differing by 0.0208 um
# between two examples for weeks -- the same failure mode this session's own
# taper-offset x mm_length sweeps re-guarded against explicitly). NOT the
# pure per-offset argmin(spread) winner at THIS offset either (the confirm
# run's own merge picks the highest-throughput point under a 2% imbalance
# cap, not the single tightest-balance point) -- see validation doc step 9c
# for the full candidate history this was picked from.
MM_LENGTH = 39.90
_MM_LENGTH_GRID = np.arange(35.0, 45.001, 0.1)   # 0.1 um spacing, narrow window -- see mm_lengths below
_nearest = _MM_LENGTH_GRID[np.argmin(np.abs(_MM_LENGTH_GRID - MM_LENGTH))]
# np.isclose, not == : linspace's own float accumulation puts _nearest ~7e-15 off
# MM_LENGTH at this grid density (unlike the old integer-spaced 206.0 grid, where
# == held exactly) -- same atol=1e-9/rtol=0 convention directional_coupler_hybrid.py
# already uses for this identical "pinned value must land exactly on its own grid"
# check, tight enough to still catch a real mis-snap (e.g. the coupler's own 0.0208
# um discrepancy) while tolerating genuine floating-point noise.
assert np.isclose(_nearest, MM_LENGTH, rtol=0, atol=1e-9), (
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
# 4. Wavelength + mm_length sweep -- one solve_length_sweep() call per output
#    branch (four solves, no mirror-symmetry shortcut -- see module
#    docstring). The body (region BODY_REGION=2 of the 5-region chain) is a
#    uniform mm_width slab, so EMESolver.solve_length_sweep() solves its
#    MODES once per wavelength and applies every swept length via meow's
#    compute_s_matrix() cascade -- BUT that cascade call itself is NOT free
#    per length at this chain size (19 cells: wg_in|8-cell taper_in|body|
#    8-cell taper_out|wg_out), unlike the single-cell body-only sweep
#    mmi_1x4_hybrid.py uses. See module docstring's 2026-08-14 finding for
#    the full story -- mm_lengths is deliberately a NARROW window now, not
#    the historical full-range grid.
#
#    NOTE: solve_length_sweep() caches under its own method name, a genuinely
#    different cache entry from solve()'s -- the first run after this change
#    needs a fresh prewarm even though save_dir is unchanged (see the
#    validation doc's own note on this, mmi_1x4_sin300_validation.md step 6).
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
mm_lengths = np.arange(35.0, 45.001, 0.1)   # um, 0.1 um spacing, narrow window; MM_LENGTH=39.90 lands exactly on this grid -- see module docstring
BODY_REGION = 2   # index into the 5-region chain: wg_in(0) taper_in(1) body(2) taper_out(3) wg_out(4)
print(f"\nSweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um, "
      f"mm_length {mm_lengths[0]:.1f} - {mm_lengths[-1]:.1f} um ({len(mm_lengths)} pts)\n")

L_idx = int(np.argmin(np.abs(mm_lengths - MM_LENGTH)))
assert np.isclose(mm_lengths[L_idx], MM_LENGTH, rtol=0, atol=1e-9), (   # see the MM_LENGTH grid-check comment above for why isclose, not ==
    f"MM_LENGTH={MM_LENGTH} is not an exact point on this sweep's own grid "
    f"(nearest: {mm_lengths[L_idx]}) -- do not silently snap, this is exactly "
    f"the coupler bug (0.0208 um discrepancy shipped for weeks)."
)

smat = None
T_grids = np.zeros((4, len(mm_lengths), len(wavelengths)))
for k in range(4):
    results_k = solver.solve_length_sweep(
        mmi, wavelengths, mm_lengths, ports=[f"out{k}"],
        sweep_region=BODY_REGION, output_idx=k, save_dir=save_dir,
    )   # list[SMatrix], aligned with mm_lengths
    T_grids[k] = np.array([sm.transmission(f"out{k}") for sm in results_k])
    print(f"  output {k} (y={mmi.output_ys[k]:+.2f} um) solved")

    smat_k = results_k[L_idx]
    if smat is None:
        smat = smat_k
    else:
        smat._data.update(smat_k._data)
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

# ---------------------------------------------------------------------------
# 7. Multimode-length sweep plot -- band-mean transmission per output +
#    max-min spread across the 4 outputs vs mm_length. Same metric/format as
#    validation_markdown/mmi_1x4_sin300_validation.md step 5's own sweep
#    (mmi_1x4_sin_body_length_sweep.py) and mmi_1x4_hybrid.py's own matching
#    sweep plot, folded into this script rather than living only in
#    verification/, since it comes free from step 4's own solve_length_sweep()
#    call above (already paid for, not an extra solve).
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
colors = ["steelblue", "darkorange", "seagreen", "firebrick"]

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
ax.set_title("EME -- output power vs body length")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[1]
ax.semilogy(mm_lengths, spread, color="k")
ax.axvline(MM_LENGTH, color="purple", linestyle="-.", lw=1.2)
ax.axvline(mm_length_sweep_best, color="gray", linestyle=":", lw=1.2)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("max-min spread across 4 outputs")
ax.set_title("EME -- equal-split quality vs body length")
ax.grid(True, which="both")

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- EME body length sweep (analytic, no extra mode solves)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "mm_length_sweep.png"), dpi=300)
plt.close("all")
print(f"Length-sweep plot saved -> {output_dir}/mm_length_sweep.png")
