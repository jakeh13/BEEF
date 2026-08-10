"""MMI body length (mm_length) design sweep for the 1x4 MMI on 300 nm Si3N4 --
HYBRID variant, companion to mmi_1x4_sin_body_length_sweep.py (the EME
variant, step 5's original design tool).

WHY THIS IS CHEAP, UNLIKE THE EME TOOL'S OWN SWEEP
---------------------------------------------------
The EME design tool's own sweep (job array 11667114) took ~5:49-6:29 hrs per
output branch because a pure-EME chain has no FDTD section to lean a cache
on -- every one of its 19 cross-sections (wg_in | 8x taper_in | body | 8x
taper_out | wg_out) is solved fresh.

HybridSolver.solve_length_sweep()'s own docstring says the opposite is true
here: "every other section is solved exactly as in solve(), with
byte-identical cache keys, so sections dedupe freely between solve() and
solve_length_sweep() calls sharing one cache dict". mmi_1x4_hybrid.py's real
run (job 11696233's prewarm + job 11734486's real run, 2026-08-05/07) already
populated example_sim_cache/ for BOTH FDTD sections (input, shared across all
four branches; output, once per branch) at N_MODES=12 -- and those sections'
geometry does NOT depend on the swept body length at all (the taper+stub
shape is fixed; only the uniform body's OWN internal propagation length
varies, which solve_length_sweep()'s analytic trick handles independently of
the drawn geometry, same trick the EME tool's own solve_length_sweep()
already uses on its body region). So this script should hit that EXISTING
cache for both FDTD sections and only pay for ONE fresh body cross-section
mode solve (tidy3d, N_MODES=12) per wavelength -- the same "single solve,
free length sweep" cost class as the real mmi_1x4_hybrid.py run itself (job
11734486: 23.5 s), not the EME tool's multi-hour-per-branch cost.

**For that cache reuse to actually hit**, this script is pinned BIT-FOR-BIT
to mmi_1x4_hybrid.py's real config -- same mm_length=206.0 (NOT a made-up
reference value the way the EME tool used 150.0, since here the FDTD
sections' absolute geometry coordinates depend on mm_length via `boundaries`,
and must match what was actually cached), same wavelengths array
(np.linspace(1.50, 1.60, 11) exactly -- resolved_material_config() folds the
full array into the cache fingerprint, not just its range), same N_MODES=12,
same body_stub=0.25. Any of these differing would cache-miss and re-trigger
the ~12 hr/section FDTD solves this script exists to avoid.

Unlike the EME tool (split into 4 SLURM array tasks + a merge script because
each branch alone took hours), this is cheap enough to do all 4 branches
sequentially in ONE process, sharing one cache dict across ports the same
way mmi_1x4_hybrid.py's own real run does -- no array/merge split needed.

Produces the same two-panel plot format as
mmi_1x4_sin_body_length_sweep_merge.py's own mm_length_sweep.png (band-mean
transmission per output vs length; log-scale spread vs length), for direct
visual comparison against the EME variant.

Run:
    python /path/to/examples/verification/python/mmi_1x4_hybrid_body_length_sweep.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, EMESolver, FDTDSolver, HybridSolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "verification", "mmi_1x4_hybrid_body_length_sweep"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Pinned BIT-FOR-BIT to mmi_1x4/python/mmi_1x4_hybrid.py -- see module
# docstring above for why this must match exactly, not just approximately,
# for the cache-reuse this script's whole cost argument depends on.
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 2.4, 12.0
MM_LENGTH = 206.0   # the REAL, cached device length -- NOT a sweep placeholder
N_MODES = 12

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30, use_neighbor_aperture=True,
)
print(f"MMI total length: {mmi.total_length:.2f} um")
print(f"Output ports at y = {mmi.output_ys} um")

_fdtd_common = dict(
    res=RESOLUTION, sim_width=MM_WIDTH + 8.0, sim_height=SIN_300.sz, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SIN_300, decay_tol=1e-11,
)
input_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=MM_WIDTH + 8.0, oxide_height=SIN_300.sz,
    res=RESOLUTION, material_stack=SIN_300, mode_solver="tidy3d",
)
output_solver = FDTDSolver(**_fdtd_common)

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
body_stub = 0.25
boundaries = [-half_total, -half_body + body_stub, +half_body - body_stub, +half_total]
interface_modes = [N_MODES, N_MODES]

hybrid = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries, interface_modes=interface_modes,
)
BODY_SECTION = 1   # index into solvers=[input(0), body(1), output(2)]

# ---------------------------------------------------------------------------
# Same 281-point grid as the EME variant, for direct comparability.
# solve_length_sweep()'s analytic trick makes the grid size free regardless.
#
# ⚠ body_lengths, NOT mm_lengths, is what gets swept -- MEASURED 2026-08-09/10:
# HybridSolver._stubs(k) returns (0.0, 0.0) for the body section (neither the
# first nor the last section), so solve_section()'s own `length = (x1-x0) +
# stub_left + stub_right` reduces to just `x1-x0` = mm_length - 2*body_stub
# for the body, NOT mm_length -- body_stub of MMI body is folded into each
# adjacent FDTD section's own domain instead. solve_section_length_sweep()
# applies whatever `lengths` it's given directly as the propagation length
# with NO such adjustment, so the caller must pre-subtract 2*body_stub --
# exactly what mmi_2x2_hybrid.py's own hand-rolled cascade already does
# (`body_lengths = mm_lengths - 2 * body_stub`). Missing this the first time
# through produced a real, measured ~6% relative error on individual port T
# (0.1635 vs the correct 0.1541 at output 0, λ=1.55, mm_length=206) -- large
# because T(L) oscillates rapidly with length here, so even a 0.5 µm
# propagation-length error lands on a materially different point of the
# curve. Verified after the fix: solve(), solve_length_sweep() with a
# single-element body_lengths, and this full 281-point sweep's own slice at
# mm_length=206 all agree exactly.
# ---------------------------------------------------------------------------
mm_lengths = np.linspace(20.0, 300.0, 281)
body_lengths = mm_lengths - 2 * body_stub   # EME body section's own physical length
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nBody length sweep: {mm_lengths[0]:.1f}-{mm_lengths[-1]:.1f} um "
      f"({len(mm_lengths)} pts), wavelengths {wavelengths[0]:.3f}-{wavelengths[-1]:.3f} um "
      f"({N_WL} pts), {N_MODES} interface modes\n")

# ---------------------------------------------------------------------------
# One shared cache dict across all 4 branches -- input(0) and body(1) are
# identical for every branch and should only ever be solved once between
# them; only output(2)'s geometry_kwargs (output_idx) differs per call.
# ---------------------------------------------------------------------------
section_cache = {}
T_grids = np.zeros((4, len(mm_lengths), len(wavelengths)))
for k in range(4):
    coeffs = hybrid.solve_length_sweep(
        mmi, wavelengths, body_lengths, sweep_section=BODY_SECTION,
        geometry_kwargs={2: {"output_idx": k}},
        cache=section_cache, save_dir=save_dir,
    )   # [n_L, n_wl, 1, 1] complex
    T_grids[k] = np.abs(coeffs[:, :, 0, 0]) ** 2
    print(f"  output {k} (y={mmi.output_ys[k]:+.2f} um) solved")

T_total_grid = T_grids.sum(axis=0)   # [n_L, n_wl]
T_band_mean = T_grids.mean(axis=2)   # [4, n_L]
spread = T_band_mean.max(axis=0) - T_band_mean.min(axis=0)   # [n_L]
best_idx = int(np.argmin(spread))
mm_length_best = float(mm_lengths[best_idx])

print(f"\nBest equal-split length (min max-min spread across 4 outputs): "
      f"mm_length = {mm_length_best:.3f} um")
print(f"  spread at this length: {spread[best_idx]:.4f}")
print(f"  per-output band-mean T: {[f'{t:.4f}' for t in T_band_mean[:, best_idx]]}")
print(f"  total T: {T_total_grid[best_idx].mean():.4f}")

# The actual chosen production point (206.0), for direct comparison against
# the sweep's own optimum -- exact grid hit (linspace(20,300,281) -> integer
# spacing), same assertion discipline mmi_1x4_hybrid.py's own MM_LENGTH uses.
_chosen_idx = int(np.argmin(np.abs(mm_lengths - MM_LENGTH)))
assert mm_lengths[_chosen_idx] == MM_LENGTH, (
    f"MM_LENGTH={MM_LENGTH} is not an exact point on this sweep's own grid "
    f"(nearest: {mm_lengths[_chosen_idx]}) -- do not silently snap."
)
print(f"\nAt the CHOSEN production mm_length={MM_LENGTH}: "
      f"spread={spread[_chosen_idx]:.4f}, "
      f"per-output band-mean T={[f'{t:.4f}' for t in T_band_mean[:, _chosen_idx]]}, "
      f"total T={T_total_grid[_chosen_idx].mean():.4f}")

np.savez(
    os.path.join(output_dir, "mm_length_sweep_hybrid.npz"),
    mm_lengths=mm_lengths, wavelengths=wavelengths,
    T_grids=T_grids, T_total_grid=T_total_grid, spread=spread,
    mm_length_best=mm_length_best, mm_length_chosen=MM_LENGTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_width=TAPER_WIDTH,
    taper_length=TAPER_LENGTH, waveguide_offset=WAVEGUIDE_OFFSET,
    mm_width=MM_WIDTH, n_modes=N_MODES, output_ys=mmi.output_ys,
)
print(f"\nData saved -> {output_dir}/mm_length_sweep_hybrid.npz")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

ax = axes[0]
colors = ["steelblue", "darkorange", "seagreen", "firebrick"]
for k in range(4):
    ax.plot(mm_lengths, T_band_mean[k], label=f"output {k} (y={mmi.output_ys[k]:+.2f} um)",
             color=colors[k])
ax.plot(mm_lengths, T_total_grid.mean(axis=1), "k--", lw=1, label="total")
ax.axvline(mm_length_best, color="gray", linestyle=":", lw=1.2,
           label=f"sweep optimum={mm_length_best:.1f} um")
ax.axvline(MM_LENGTH, color="purple", linestyle="-.", lw=1.2,
           label=f"chosen (production)={MM_LENGTH:.1f} um")
ax.axhline(0.25, color="gray", linestyle=":", lw=0.6)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("Band-mean transmission")
ax.set_title("1x4 MMI (hybrid) -- output power vs body length")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[1]
ax.semilogy(mm_lengths, spread, color="k")
ax.axvline(mm_length_best, color="gray", linestyle=":", lw=1.2)
ax.axvline(MM_LENGTH, color="purple", linestyle="-.", lw=1.2)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("max-min spread across 4 outputs")
ax.set_title("Equal-split quality vs body length (hybrid)")
ax.grid(True, which="both")

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- body length design sweep (HYBRID)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "mm_length_sweep_hybrid.png"), dpi=300)
plt.close("all")
print(f"Plot saved -> {output_dir}/mm_length_sweep_hybrid.png")
