"""
Tutorial: 2×2 MMI — FDTD/EME hybrid simulation (FECES method).

HybridSolver is a generic N-section cascade: any solver type, in any order,
any count >= 2. Here we chain three sections —

    [input section]  ->  [body section]  ->  [output section]

— each independently configurable as FDTD or EME via the METHOD switches
below. HybridSolver itself has no idea any of these are FDTD/EME or that the
device is an MMI; it only knows each solver implements solve_section() and
combines the resulting per-section modal transfer matrices by matrix-chain
multiplication.  The MMI-specific geometry (tapers, body, branches) lives
entirely in MMI's section_prisms/section_gds/width_at/center_at methods.

Unlike the 1×2 splitter, a 2×2 has two ports on *both* sides sitting close
together (waveguide_offset = 1.5 µm), so both the input and the output FDTD
sections use half_domain_aperture to isolate the excited/measured port from
its neighbor. The top input is excited (input_idx=1); the top output is the
"bar" port (output_idx=1) and the bottom output the "cross" port (output_idx=0).

Sim count: 2 FDTD + 1 EME mode-solve (mirror reuse on the bar-output
section via HybridSolver's mirror_pairs). The MMI is left-right
mirror-symmetric, so the output section with output_idx=1 ("bar") is
exactly the x-mirror of the input section with input_idx=1 (the fixed
excitation) — its transfer matrix is just the input section's transposed,
not resimulated. The cross output (output_idx=0) does NOT mirror the fixed
input_idx=1 excitation (it would only mirror an input_idx=0 excitation,
which is never run here), so it's still solved directly; getting cross "for
free" too would need the directional_coupler_hybrid.py-style manual bypass
(solving both input branches directly and deriving every output via
transpose) — not done here since it saves nothing further for this
example's bar/cross metric.

On top of the single-length result, the multimode-section length is swept
via HybridSolver.solve_length_sweep() — the EME body is a uniform
mm_width slab, so its modes are solved once per wavelength and each swept
length is applied as an analytic propagation phase (same trick as
directional_coupler_hybrid.py's coupling-length sweep, but through the
generic HybridSolver API instead of a manual bypass). The whole sweep adds
zero FDTD simulations and zero extra mode solves over the single-length
run: still 2 FDTD + 1 EME mode-solve total.

Run:
    python /path/to/examples/mmi_2x2/mmi_2x2_hybrid.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from feces import MMI, EMESolver, FDTDSolver, HybridSolver, SMatrix, SOI_220, plot_hybrid_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_output", "mmi_2x2", "mmi_2x2_hybrid"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by both the EME and FDTD sub-solvers
WL_MIN = 1.530      # wavelength sweep start (µm)
WL_MAX = 1.565      # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

# ---------------------------------------------------------------------------
# 0. Simulation method per section — HybridSolver is solver-agnostic (any
#    mix of FDTDSolver/EMESolver, in any order), so each of the three
#    sections can independently be "fdtd" or "eme".
# ---------------------------------------------------------------------------
INPUT_METHOD = "fdtd"    # "fdtd" or "eme"
BODY_METHOD = "eme"      # "fdtd" or "eme"
OUTPUT_METHOD = "fdtd"   # "fdtd" or "eme"

# ---------------------------------------------------------------------------
# 1. Geometry
# ---------------------------------------------------------------------------
mmi = MMI(
    n_inputs=2,             # 2x2 splitter
    n_outputs=2,
    mm_length=28.0,
    mm_width=2.8,
    waveguide_width=0.45,
    taper_length=4.0,
    taper_width=1.3,
    waveguide_offset=1.5,
    waveguide_height=0.22,  # SOI_220
)
print(f"MMI total length: {mmi.total_length:.2f} µm")
print(f"Input ports at y = {mmi.input_ys} µm; outputs at y = {mmi.output_ys} µm")

# ---------------------------------------------------------------------------
# 2. Configure sub-solvers — one per section, per the METHOD switches above.
#
#    input section  : input access waveguide + input taper + a bit of the
#                      MMI body (enough for a clean interface monitor past
#                      the taper/body cross-section discontinuity)
#    body section   : propagates N modes through the remaining rectangular
#                      MMI body
#    output section : the mirror image of the input section, one instance
#                      per output port
# ---------------------------------------------------------------------------

_fdtd_common = dict(
    res=RESOLUTION,
    sim_width=5.0,
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,       # full 3-D FDTD: real 220 nm Si core in oxide, PML on all six faces
    material_stack=SOI_220,
    decay_tol=1e-11,
)

N_MODES = 15   # modes carried at the two input/output <-> body interfaces

# Both the input and output ports have a neighbor waveguide_offset = 1.5 µm
# away that is drawn in the same section (all branches are always physically
# present). half_domain_aperture confines each section's source/norm aperture
# to the near half of the domain (from the shared y=0 symmetry axis to the
# near PML edge), structurally excluding the opposite branch.


def make_section_solver(method: str, half_domain_aperture: bool = False):
    """Build a solver for one hybrid section. Both FDTDSolver and EMESolver
    implement the same solve_section() interface, so either can fill any
    slot in the chain — num_modes/num_cells/lengths on EMESolver only
    matter for its full-device solve()/solve_sweep(), not solve_section(),
    so placeholder values are fine here."""
    if method == "fdtd":
        kwargs = dict(_fdtd_common)
        kwargs["half_domain_aperture"] = half_domain_aperture
        return FDTDSolver(**kwargs)
    if method == "eme":
        return EMESolver(
            num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
            oxide_width=5.0, oxide_height=1.5, res=RESOLUTION, material_stack=SOI_220,
            mode_solver="meep",  # matches FDTDSolver's own eigenmode basis in this hybrid chain
        )
    raise ValueError(f"unknown method {method!r}, expected 'fdtd' or 'eme'")


# Both sides have two closely-spaced ports, so both the input and output
# sections need the narrower half-domain aperture (unlike the 1×2, whose lone
# input keeps the full sim_width).
input_solver = make_section_solver(INPUT_METHOD, half_domain_aperture=True)
body_solver = make_section_solver(BODY_METHOD)
output_solver = make_section_solver(OUTPUT_METHOD, half_domain_aperture=True)

# ---------------------------------------------------------------------------
# 3. Build HybridSolver — a generic 3-section chain here, but HybridSolver
#    itself works with any solver list of length >= 2, in any order.
# ---------------------------------------------------------------------------
half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
body_stub = 0.25   # µm of MMI body folded into each FDTD section, past the
                   # taper/body cross-section discontinuity

boundaries = [
    -half_total,               # physical input port (input taper tip)
    -half_body + body_stub,    # input section <-> body section interface
    +half_body - body_stub,    # body section <-> output section interface
    +half_total,               # physical output port (output taper tip)
]
interface_modes = [N_MODES, N_MODES]

# Multimode-length sweep grid. The geometry stays built at the nominal
# mm_length — only the body cross-section (a length-independent mm_width
# slab) feeds the one EME mode solve, and each swept length is applied
# analytically. Each FDTD section folds body_stub µm of MMI body into
# itself, so the EME body section's physical length is mm_length minus
# both stubs.
mm_lengths = np.linspace(16.0, 40.0, 201)   # µm; nominal 28.0 lands exactly on index 100
body_lengths = mm_lengths - 2 * body_stub   # EME body section's physical length

hybrid = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
)

# Second instance sharing the same solvers/boundaries/interface_modes, with
# the output section (idx 2) declared as the x-mirror of the input section
# (idx 0) — valid only when output_idx == input_idx (see module docstring).
# Used for the "bar" call only; "cross" still goes through the plain
# `hybrid` above since output_idx=0 != the fixed input_idx=1.
hybrid_mirror_bar = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
    mirror_pairs={2: 0},
)

# The top input is excited for every solve below (input_idx=1); only the
# output section's output_idx varies (1 = top = bar, 0 = bottom = cross).
IN_KW = {"input_idx": 1}

# ---------------------------------------------------------------------------
# 4. Geometry plot — shows which solver covers which section
# ---------------------------------------------------------------------------
plot_hybrid_geometry(
    mmi, solver=hybrid,
    in_modes=1, out_modes=1,
    title="2×2 MMI — Hybrid Solver Sections",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True,
    sub_region_show_eme_detail=True,
    # One set of sub-region plots per output port — matches the two
    # hybrid.solve() calls below. The input section (fixed input_idx=1) and
    # body section don't vary between ports, so they're rendered once.
    sub_region_configs=[
        ("bar",   [IN_KW, {}, {"output_idx": 1}]),
        ("cross", [IN_KW, {}, {"output_idx": 0}]),
    ],
)
print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")
print(f"Sub-region plots saved → {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run hybrid length sweep — once per output port. input_idx=1 (top input)
#    is fixed across both calls, so with a shared cache the input and body
#    sections are only ever simulated once. The bar call's output section
#    (output_idx=1, matching the fixed input_idx=1) is mirror-derived via
#    hybrid_mirror_bar instead of resimulated; the cross call's output
#    section (output_idx=0) is genuinely different from the excited port
#    and is solved directly through the plain (unmirrored) hybrid. The body
#    section is length-swept analytically from a single mode solve. Total:
#    2 FDTD (input, cross-output) + 1 EME mode-solve (body) — identical to
#    a single-length run, now yielding all 201 lengths.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nHybrid sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, {N_MODES} interface modes), "
      f"mm_length {mm_lengths[0]:.1f} – {mm_lengths[-1]:.1f} µm "
      f"({len(mm_lengths)} pts)\n")

section_cache = {}
coeffs_bar = hybrid_mirror_bar.solve_length_sweep(
    mmi, wavelengths, body_lengths, sweep_section=1, cache=section_cache,
    geometry_kwargs=[IN_KW, {}, {"output_idx": 1}])   # [n_L, n_wl, 1, 1]
coeffs_cross = hybrid.solve_length_sweep(
    mmi, wavelengths, body_lengths, sweep_section=1, cache=section_cache,
    geometry_kwargs=[IN_KW, {}, {"output_idx": 0}])

T_bar_grid   = np.abs(coeffs_bar[:, :, 0, 0]) ** 2    # [n_L, n_wl]
T_cross_grid = np.abs(coeffs_cross[:, :, 0, 0]) ** 2
T_total_grid = T_bar_grid + T_cross_grid
ratio_grid   = T_bar_grid / np.clip(T_total_grid, 1e-20, None)  # bar fraction

# Slice at the nominal mm_length for the wavelength-domain view — keep the
# complex coefficients for this slice so the saved SMatrix preserves phase.
L_plot_idx = int(np.argmin(np.abs(mm_lengths - mmi.mm_length)))  # exact: 28.0 is on-grid
L_plot = mm_lengths[L_plot_idx]

T_bar   = T_bar_grid[L_plot_idx]
T_cross = T_cross_grid[L_plot_idx]
T_total = T_total_grid[L_plot_idx]

wl_c = wavelengths[len(wavelengths) // 2]
idx  = np.argmin(np.abs(wavelengths - wl_c))

print(f"\nResults at λ = {wl_c:.3f} µm, mm_length = {L_plot:.1f} µm:")
print(f"  Bar output   T = {T_bar[idx]:.4f}  IL = {10*np.log10(max(T_bar[idx],1e-20)):.2f} dB")
print(f"  Cross output T = {T_cross[idx]:.4f}  IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total        T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  (ideal 50/50 split: T_each=0.5, IL_total=0.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat = SMatrix(wavelengths=wavelengths)
smat.set(coeffs_bar[L_plot_idx, :, 0, 0], port="bar")
smat.set(coeffs_cross[L_plot_idx, :, 0, 0], port="cross")
smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix (mm_length={L_plot:.1f} µm slice) saved → {output_dir}/smatrix_hybrid.pkl")

np.savez(os.path.join(output_dir, "mm_length_sweep.npz"),
         mm_lengths=mm_lengths, wavelengths=wavelengths,
         T_bar=T_bar_grid, T_cross=T_cross_grid)
print(f"Length-sweep grids saved → {output_dir}/mm_length_sweep.npz")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_bar,   label="Bar output",   color="steelblue")
ax.plot(wavelengths, T_cross, label="Cross output", color="darkorange", linestyle="--")
ax.plot(wavelengths, T_total, label="Total",         color="k", lw=1.5)
ax.axhline(0.5, color="gray", linestyle=":", lw=0.8, label="Ideal per-port")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title(f"Hybrid — Transmission vs Wavelength (mm_length={L_plot:.1f} µm)")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)

ax = axes[1]
IL_bar   = 10 * np.log10(np.clip(T_bar,   1e-20, None))
IL_cross = 10 * np.log10(np.clip(T_cross, 1e-20, None))
IL_tot   = 10 * np.log10(np.clip(T_total, 1e-20, None))
ax.plot(wavelengths, IL_bar,   label="Bar output",   color="steelblue")
ax.plot(wavelengths, IL_cross, label="Cross output", color="darkorange", linestyle="--")
ax.plot(wavelengths, IL_tot,   label="Total",         color="k", lw=1.5)
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title(f"Hybrid — Insertion Loss vs Wavelength (mm_length={L_plot:.1f} µm)")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle("2×2 MMI — FDTD/EME Hybrid Simulation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/hybrid_results.png")

# -- Multimode-length sweep plot: mean over wavelength, shaded min–max -----
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

ax = axes[0]
for grid, label, color in [
    (T_bar_grid,   "Bar",   "steelblue"),
    (T_cross_grid, "Cross", "darkorange"),
    (T_total_grid, "Total", "k"),
]:
    mean, lo, hi = grid.mean(axis=1), grid.min(axis=1), grid.max(axis=1)
    ax.plot(mm_lengths, mean, label=label, color=color,
            lw=1.5 if label == "Total" else 1.2,
            linestyle="--" if label == "Cross" else "-")
    ax.fill_between(mm_lengths, lo, hi, color=color, alpha=0.15, lw=0)
ax.axvline(L_plot, color="gray", linestyle=":", lw=0.8,
           label=f"Nominal ({L_plot:.1f} µm)")
ax.axhline(0.5, color="gray", linestyle=":", lw=0.8)
ax.set_xlabel("Multimode section length (µm)")
ax.set_ylabel("Transmission (mean over wavelength, shaded = min–max)")
ax.set_title("Transmission vs MM Length")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)

ax = axes[1]
r_mean, r_lo, r_hi = ratio_grid.mean(axis=1), ratio_grid.min(axis=1), ratio_grid.max(axis=1)
ax.plot(mm_lengths, r_mean, color="seagreen", lw=1.2, label="Bar fraction")
ax.fill_between(mm_lengths, r_lo, r_hi, color="seagreen", alpha=0.15, lw=0)
ax.axhline(0.5, color="gray", linestyle=":", lw=0.8, label="50:50")
ax.axvline(L_plot, color="gray", linestyle=":", lw=0.8,
           label=f"Nominal ({L_plot:.1f} µm)")
ax.set_xlabel("Multimode section length (µm)")
ax.set_ylabel("Bar fraction  T_bar / (T_bar + T_cross)")
ax.set_title("Splitting Ratio vs MM Length")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1)

plt.suptitle("2×2 MMI — Multimode-Length Sweep (analytic, no extra simulations)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "mm_length_sweep.png"), dpi=300)
plt.close("all")
print(f"Length-sweep plot saved → {output_dir}/mm_length_sweep.png")
