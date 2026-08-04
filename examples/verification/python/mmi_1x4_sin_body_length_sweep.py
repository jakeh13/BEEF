"""MMI body length (mm_length) design sweep for the 1x4 MMI on 300 nm Si3N4 --
validation_markdown/mmi_1x4_sin300_validation.md step 5.

Same trick directional_coupler_hybrid.py uses for L_50 and mmi_2x2_hybrid.py
uses for its own mm_length: EMESolver.solve_section_length_sweep() solves the
MMI body's modes ONCE per wavelength (it is a uniform mm_width slab -- the
one thing that method requires) and evaluates every swept length
analytically, so a 200-point mm_length sweep costs one mode solve, not 200.

Unlike mmi_2x2_hybrid.py, this device has NO left-right mirror-symmetry
shortcut to reuse: 2x2 has the same port arrangement on both sides of the
body, so its output section is literally the x-mirror of its input section
(by reciprocity) and is never separately simulated. A 1x4 has one input port
and four output ports -- structurally different sides -- so this script
solves the input section once (input_idx=0, the only input) and each output
branch separately (output_idx=0..3), four solves. All-EME (no FDTD): this is
a design/measurement tool for finding mm_length, not the final hybrid
example (that is step 6's mmi_1x4_hybrid.py, built once mm_length is known).

Geometry, from steps 4-5:
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
                                 mm_width/2=6.0 um, 1.30 um clearance) before
                                 spending any solve on it.

Only the fundamental input mode is launched (n_modes_in=1) and only each
output's fundamental mode is measured (n_modes_out=1) -- the standard,
established way to find a self-imaging length, and consistent with every
existing MMI example. The access waveguide's SECOND guided mode (step 4's
follow-up measurement) is deliberately out of scope here; it is flagged in
the validation doc as a step 6 concern (mode identity/correspondence at a
genuinely 2-mode access waveguide), not resolved by this script.

N_MODES=20 for the body: a real mode-solve on the mm_width=12.0 body (job
11660226) found 23 guided modes, margins above n_clad from +0.228 (mode 0)
down to +0.004 (mode 22, barely guided). The bottom 3 (margins <0.03) are
excluded -- close enough to cutoff that they are poorly converged and
contribute negligible overlap with a comfortably-confined 1.1 um input mode
anyway (this is exactly why the 1x2/2x2 examples drop marginal bands too --
mmi_2x2_hybrid.py's own N_MODES comment: "was 15: band 14 sits below the
cladding index -- a radiation mode, dropped at no cost").

Run:
    python /path/to/examples/verification/python/mmi_1x4_sin_body_length_sweep.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, EMESolver, plot_eme_geometry
from beef.materials.stack import SIN_300, sellmeier_index

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification", "mmi_1x4_sin_body_length_sweep"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # pts/um -- this codebase's standing mode-solve default
WL_MIN = 1.50       # um, C-band (step 4c's settled choice)
WL_MAX = 1.60
N_WL = 11
N_MODES = 20        # body modes carried -- see module docstring

# ---------------------------------------------------------------------------
# 1. Geometry -- built at a REFERENCE mm_length for the input/output taper
#    polygons (those don't depend on mm_length); the body's own length is
#    swept separately below and applied analytically.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
WAVEGUIDE_OFFSET = 2.4
MM_WIDTH = 12.0

# Starting bracket only (plan's own instruction: do not use a formula alone).
# L_pi = 4*n_r*W_e^2/(3*lambda), n_r ~ 1.67 (this platform's confined-mode
# index, from step 4/5's own measurements), W_e ~ mm_width for a body this
# much wider than the confinement scale. Paired-interference 1xN imaging
# needs asymmetric consideration of L_pi -- this is a BRACKET, not a
# prediction, and the swept curve (not the formula) picks the real length.
_n_r_estimate = 1.674
_L_pi_estimate = 4 * _n_r_estimate * MM_WIDTH**2 / (3 * 1.55)
print(f"Starting-bracket estimate: L_pi ~ {_L_pi_estimate:.1f} um "
      f"(n_r~{_n_r_estimate}, W_e~{MM_WIDTH} um) -- a bracket seed, not a prediction")

MM_LENGTH_REF = 150.0   # reference value for geometry/taper polygons only

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
# the inequality the plan's step 5 explicitly asks to write down and check.
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
# 2. EME solver (shared config for input/body/output -- all sections use the
#    same mesh/material/resolution; only n_modes_in/out and x0/x1 differ per
#    call). All-EME: this is a design tool, not the final hybrid device.
# ---------------------------------------------------------------------------
solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholder -- solve_section(_length_sweep) sweeps its own
    oxide_width=MM_WIDTH + 8.0,   # generous lateral clearance past the body edges
    oxide_height=SIN_300.sz,      # 5.0 um, step 4's converged value
    res=RESOLUTION,
    material_stack=SIN_300,
    mode_solver="tidy3d",         # step 4's fix -- compute_modes_meep silently
                                  # returns wrong answers for this weakly-
                                  # confined material; see that step's writeup
)

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
# 4. Boundaries. body_stub=0 -- unlike the FDTD-adjacent hybrid examples,
#    every section here is EME, so there is no monitor-clearance reason to
#    fold a stub of body into the neighbouring section; the body/taper
#    interface (a legal cut -- a band endpoint, not mid-taper) is used
#    directly.
# ---------------------------------------------------------------------------
half_body_ref = MM_LENGTH_REF / 2
half_total_ref = mmi.total_length / 2
boundaries = [-half_total_ref, -half_body_ref, +half_body_ref, +half_total_ref]

# mm_length sweep grid. See _L_pi_estimate above for why this range -- wide
# on purpose, since the analytic cost is the same for 20 points or 200.
mm_lengths = np.linspace(20.0, 300.0, 281)
body_lengths = mm_lengths   # body_stub=0, so no subtraction needed

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nBody length sweep: {mm_lengths[0]:.1f}-{mm_lengths[-1]:.1f} um "
      f"({len(mm_lengths)} pts), wavelengths {wavelengths[0]:.3f}-{wavelengths[-1]:.3f} um "
      f"({N_WL} pts), {N_MODES} body modes\n")

# ---------------------------------------------------------------------------
# 5. Solve: 1 input section + 4 output sections (one per output_idx) +
#    1 body length sweep. No mirror-symmetry shortcut available (see module
#    docstring) -- 6 solves total, all cheap (single-cell EME, no time-
#    stepping).
# ---------------------------------------------------------------------------
M_in = solver.solve_section(
    mmi, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES, input_idx=0,
)  # [n_wl, N_MODES, 1]

M_out = []
for k in range(4):
    M_out_k = solver.solve_section(
        mmi, wavelengths, boundaries[2], boundaries[3],
        n_modes_in=N_MODES, n_modes_out=1, output_idx=k,
    )  # [n_wl, 1, N_MODES]
    M_out.append(M_out_k)
    print(f"  output section {k} (y={mmi.output_ys[k]:+.2f} um) solved")

M_body = solver.solve_section_length_sweep(
    mmi, wavelengths, boundaries[1], boundaries[2], body_lengths,
    n_modes_in=N_MODES, n_modes_out=N_MODES,
)  # [n_L, n_wl, N_MODES, N_MODES]

tmp = np.einsum("Lwij,wjk->Lwik", M_body, M_in)   # [n_L, n_wl, N_MODES, 1]
T_grids = []
for k in range(4):
    coeffs_k = np.einsum("wij,Lwjk->Lwik", M_out[k], tmp)   # [n_L, n_wl, 1, 1]
    T_grids.append(np.abs(coeffs_k[:, :, 0, 0]) ** 2)        # [n_L, n_wl]
T_grids = np.stack(T_grids, axis=0)   # [4, n_L, n_wl]
T_total_grid = T_grids.sum(axis=0)    # [n_L, n_wl]

# ---------------------------------------------------------------------------
# 6. Pick the best equal-split length: minimize the spread across the 4
#    outputs' band-mean transmission (the 1x4 analogue of the coupler's
#    50:50 crossing / the 2x2's argmin-to-nominal, generalized to N ports).
# ---------------------------------------------------------------------------
T_band_mean = T_grids.mean(axis=2)   # [4, n_L] -- band-mean per output per length
spread = T_band_mean.max(axis=0) - T_band_mean.min(axis=0)   # [n_L]
best_idx = int(np.argmin(spread))
mm_length_best = float(mm_lengths[best_idx])

print(f"\nBest equal-split length (min max-min spread across 4 outputs): "
      f"mm_length = {mm_length_best:.3f} um")
print(f"  spread at this length: {spread[best_idx]:.4f}")
print(f"  per-output band-mean T: {[f'{t:.4f}' for t in T_band_mean[:, best_idx]]}")
print(f"  total T: {T_total_grid[best_idx].mean():.4f}")

# Local minima worth reporting too -- paired-interference imaging can have
# more than one good length in range, and the plan's own step 5 note #3
# ("iterate mm_width if the equal-split point is unreachable or excess loss
# is poor") means a single global argmin could hide a better, shorter option.
_is_local_min = (
    (spread[1:-1] < spread[:-2]) & (spread[1:-1] < spread[2:])
)
_local_min_idxs = np.flatnonzero(_is_local_min) + 1
print(f"\nAll local minima in the swept range (length, spread):")
for i in _local_min_idxs:
    print(f"  {mm_lengths[i]:7.2f} um   spread={spread[i]:.4f}   "
          f"total T={T_total_grid[i].mean():.4f}")

# ---------------------------------------------------------------------------
# 7. Save + plot
# ---------------------------------------------------------------------------
np.savez(
    os.path.join(output_dir, "mm_length_sweep.npz"),
    mm_lengths=mm_lengths, wavelengths=wavelengths,
    T_grids=T_grids, T_total_grid=T_total_grid, spread=spread,
    mm_length_best=mm_length_best,
    waveguide_width=WAVEGUIDE_WIDTH, taper_width=TAPER_WIDTH,
    taper_length=TAPER_LENGTH, waveguide_offset=WAVEGUIDE_OFFSET,
    mm_width=MM_WIDTH, n_modes=N_MODES,
)
print(f"\nData saved -> {output_dir}/mm_length_sweep.npz")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

ax = axes[0]
colors = ["steelblue", "darkorange", "seagreen", "firebrick"]
for k in range(4):
    ax.plot(mm_lengths, T_band_mean[k], label=f"output {k} (y={mmi.output_ys[k]:+.2f} um)",
           color=colors[k])
ax.plot(mm_lengths, T_total_grid.mean(axis=1), "k--", lw=1, label="total")
ax.axvline(mm_length_best, color="gray", linestyle=":", lw=1.2,
          label=f"chosen mm_length={mm_length_best:.1f} um")
ax.axhline(0.25, color="gray", linestyle=":", lw=0.6)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("Band-mean transmission")
ax.set_title("1x4 MMI -- output power vs body length")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[1]
ax.semilogy(mm_lengths, spread, color="k")
ax.axvline(mm_length_best, color="gray", linestyle=":", lw=1.2)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("max-min spread across 4 outputs")
ax.set_title("Equal-split quality vs body length")
ax.grid(True, which="both")

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- body length design sweep", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "mm_length_sweep.png"), dpi=300)
plt.close("all")
print(f"Plot saved -> {output_dir}/mm_length_sweep.png")
