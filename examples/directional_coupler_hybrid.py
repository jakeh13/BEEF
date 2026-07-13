"""
Tutorial: standard directional coupler — FDTD/EME hybrid simulation (FECES method).

Geometry (220 nm SOI, C-band):

    "in"  ___                                      ____  "through"
             \\                                    /
              \\__________________________________/        (coupling gap)
              /                                    \\
    ___      /                                      \\___  "cross"
   "isolated"

    |<- bend_length ->|<--- coupling_length --->|<- bend_length ->|

Three HybridSolver sections, one per physically distinct feature:

    [input S-bend]  ->  [straight coupled pair]  ->  [output S-bend]
       FDTD                    EME                       FDTD

  - The S-bends that carry the excited waveguide from a well-isolated access
    pitch down into the tight coupling gap (and back out again) are curved,
    non-adiabatic features — exactly what FDTD is for.
  - The straight coupled-pair region is where the actual directional
    coupling happens: an eigenmode expansion over its 2 supermodes (even +
    odd) is both exact and far cheaper than FDTD here. Each individual arm
    is single-mode by design, so 2 modes fully capture the excited coupling
    physics (bump N_MODES up only if the arm width supports higher-order
    per-arm modes too).
  - The output section is solved once per output port ("through" and
    "cross"), exactly like the MMI hybrid tutorial. Since each output port
    is single-mode (n_modes_out=1) while the interface to the EME body
    carries N_MODES=2, FDTDSolver.solve_section() automatically excites
    from the *output* side and propagates backward, relying on reciprocity
    to get the same transfer matrix a forward source would — the "reverse
    simulated" FDTD run — for free, with no extra code (see its docstring
    in feces/solvers/fdtd.py).

Coupling-length sweep — the classic directional-coupler "beat length"
characterization — without any additional FDTD or EME simulations:

  - The input/output FDTD sections' transfer matrices don't depend on
    coupling_length at all (the S-bend/access-waveguide geometry is defined
    relative to the coupling-region edges, and FDTDSolver.solve_section()
    re-centers on its own padded domain, so it's translation-invariant).
    They're solved exactly once each, using one reference DirectionalCoupler
    instance built at the *middle* of the length sweep (also what the
    geometry/section preview plots use).
  - The EME body's mode profiles and effective indices depend only on the
    coupling cross-section (width/gap), not on how long the straight section
    is — so EMESolver.solve_section_length_sweep() computes the two
    supermodes once per wavelength and reuses them for every length in the
    sweep, each just a cheap diagonal phase-matrix evaluation
    (S_ii(L) = exp(2j*pi*neff_i/wavelength*L); see its docstring in
    feces/solvers/eme.py). Total solver cost is exactly what a single fixed
    coupling_length run already costs: 3 FDTD solves + 1 EME mode-solve per
    wavelength — regardless of sweeping 101 lengths instead of 1.

HybridSolver itself has no idea this is a directional coupler; it only
knows each solver implements solve_section() and combines the resulting
per-section modal transfer matrices by matrix-chain multiplication. The
coupler-specific geometry (S-bends, coupling region) lives entirely in
DirectionalCoupler's section_prisms/section_gds/width_at/center_at methods.
(HybridSolver.solve() itself has no notion of a length-swept section, so the
matrix-chain contraction below is done directly, the same way it does it
internally, broadcasting over the extra length axis.)

Run:
    python /path/to/examples/directional_coupler_hybrid.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from feces import (
    DirectionalCoupler, EMESolver, FDTDSolver, HybridSolver, SMatrix, SOI_220, plot_hybrid_geometry,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_output", "directional_coupler_hybrid"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Geometry — all device dimensions as variables (220 nm SOI).
#
#    coupling_lengths is the actual sweep; dc itself is built at the sweep's
#    middle value (bends/access-waveguide geometry don't depend on
#    coupling_length, so this single instance is equally valid as the
#    reference geometry for the FDTD sections at every length in the sweep,
#    and it's also exactly what the geometry/section preview plots should
#    show per the "plots at the middle length" convention).
# ---------------------------------------------------------------------------
coupling_lengths = np.linspace(0.0, 20.0, 101)   # µm, the swept coupling length
L_mid = coupling_lengths[len(coupling_lengths) // 2]
L_mid_idx = len(coupling_lengths) // 2

dc = DirectionalCoupler(
    waveguide_width=0.45,    # single-mode strip width (µm)
    gap=0.20,                # coupling-region edge-to-edge gap (µm)
    coupling_length=L_mid,   # straight parallel interaction length (µm) — sweep reference value
    bend_length=9.0,         # horizontal extent of each S-bend (µm)
    waveguide_offset=3.0,    # far-end centre-to-centre port pitch (µm)
    bend_type="euler",       # "euler" (gdsfactory bend_s) or "circular" (matched arcs)
    waveguide_height=0.22,   # Si layer thickness (µm) — SOI_220
)
print(f"Directional coupler total length: {dc.total_length:.2f} µm "
      f"(coupling pitch {dc.coupling_pitch:.3f} µm, far pitch {dc.waveguide_offset:.2f} µm)")
print(f"Coupling length sweep: {coupling_lengths[0]:.1f} – {coupling_lengths[-1]:.1f} µm "
      f"({len(coupling_lengths)} pts, reference/plot length = {L_mid:.2f} µm)")

# ---------------------------------------------------------------------------
# 2. Configure sub-solvers — one per section.
#
#    input section  : input access waveguide + input S-bend + a sliver of
#                      the coupling region (enough for a clean interface
#                      monitor away from the bend/coupling-region seam)
#    body section   : propagates N_MODES supermodes through the remaining
#                      straight coupled-pair region
#    output section : the mirror image of the input section, one instance
#                      per output port ("through" / "cross")
# ---------------------------------------------------------------------------
N_MODES = 2   # even + odd supermodes carried at the two body interfaces

_fdtd_common = dict(
    res=20.0,
    sim_width=6.0,        # >= waveguide_offset + margin, so both far-pitch arms sit inside PML clearance
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    material_stack=SOI_220,
    decay_tol=1e-11,
)

input_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholder length; solve_section_length_sweep() sweeps its own
    oxide_width=6.0, oxide_height=1.5, res=40.0, material_stack=SOI_220,
)
output_solver = FDTDSolver(**_fdtd_common)

# ---------------------------------------------------------------------------
# 3. Build HybridSolver — 3-section chain: FDTD bend | EME coupling | FDTD bend.
#    (Used below only for the geometry/section preview plots; the actual
#    numeric sweep bypasses HybridSolver.solve() — see step 5.)
# ---------------------------------------------------------------------------
half_body = dc.coupling_length / 2
half_total = dc.total_length / 2
body_stub = 0.5   # µm of straight coupling region folded into each FDTD
                  # section, past the S-bend/coupling-region seam — the
                  # cross-section is continuous there, but this keeps each
                  # interface monitor comfortably inside the always-both-arm
                  # region rather than exactly on it (same trick used in
                  # the MMI hybrid tutorial's taper/body seam).

boundaries = [
    -half_total,              # physical input port (access waveguide end)
    -half_body + body_stub,   # input section <-> body section interface
    +half_body - body_stub,   # body section <-> output section interface
    +half_total,              # physical output port (access waveguide end)
]
interface_modes = [N_MODES, N_MODES]

hybrid = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
)

# ---------------------------------------------------------------------------
# 4. Geometry plot — shows which solver covers which section, at dc's
#    (middle-of-sweep) coupling_length.
# ---------------------------------------------------------------------------
plot_hybrid_geometry(
    dc, solver=hybrid,
    in_modes=1, out_modes=1,
    title=f"Directional Coupler — Hybrid Solver Sections (coupling_length={L_mid:.1f} µm)",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True,
    sub_region_show_eme_detail=True,
    # One set of sub-region plots per output port — matches the two output
    # configurations solved in step 5 below; only the output section's
    # kwargs actually vary, so the input/body sections are each rendered
    # once (not once per output port).
    sub_region_configs=[
        ("through", [{}, {}, {"output_idx": 0}]),
        ("cross", [{}, {}, {"output_idx": 1}]),
    ],
)
print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")
print(f"Sub-region plots saved → {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run the sweep — 3 FDTD solves (input, output/through, output/cross) plus
#    one EME mode-solve per wavelength, all independent of coupling_length,
#    combined with the coupling_length-swept body matrix via the same
#    matrix-chain contraction HybridSolver.solve() uses internally (just
#    broadcast over the extra length axis).
# ---------------------------------------------------------------------------
wavelengths = np.linspace(1.50, 1.60, 11)
print(f"\nWavelength sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

M_in = input_solver.solve_section(
    dc, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=input_solver.waveguide_length, stub_right=0.0,
)  # [n_wl, N_MODES, 1]
M_out_through = output_solver.solve_section(
    dc, wavelengths, boundaries[2], boundaries[3],
    n_modes_in=N_MODES, n_modes_out=1,
    stub_left=0.0, stub_right=output_solver.waveguide_length,
    output_idx=0,
)  # [n_wl, 1, N_MODES]
M_out_cross = output_solver.solve_section(
    dc, wavelengths, boundaries[2], boundaries[3],
    n_modes_in=N_MODES, n_modes_out=1,
    stub_left=0.0, stub_right=output_solver.waveguide_length,
    output_idx=1,
)  # [n_wl, 1, N_MODES]
M_body = body_solver.solve_section_length_sweep(
    dc, wavelengths, boundaries[1], boundaries[2], coupling_lengths,
    n_modes_in=N_MODES, n_modes_out=N_MODES,
)  # [n_L, n_wl, N_MODES, N_MODES]

tmp = np.einsum("Lwij,wjk->Lwik", M_body, M_in)                          # [n_L, n_wl, N_MODES, 1]
full_through = np.einsum("wij,Lwjk->Lwik", M_out_through, tmp)          # [n_L, n_wl, 1, 1]
full_cross = np.einsum("wij,Lwjk->Lwik", M_out_cross, tmp)              # [n_L, n_wl, 1, 1]

T_through_grid = np.abs(full_through[:, :, 0, 0]) ** 2   # [n_L, n_wl]
T_cross_grid = np.abs(full_cross[:, :, 0, 0]) ** 2        # [n_L, n_wl]
T_total_grid = T_through_grid + T_cross_grid

# Slice at the reference (middle-of-sweep) coupling_length for the
# wavelength-domain view — keep the complex coefficients (not just power)
# for this slice so the saved SMatrix preserves phase, same as before.
coeff_through = full_through[L_mid_idx, :, 0, 0]
coeff_cross = full_cross[L_mid_idx, :, 0, 0]
T_through = T_through_grid[L_mid_idx]
T_cross = T_cross_grid[L_mid_idx]
T_total = T_total_grid[L_mid_idx]

wl_c = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_c))

print(f"Results at coupling_length = {L_mid:.2f} µm, λ = {wl_c:.3f} µm:")
print(f"  Through output T = {T_through[idx]:.4f}  IL = {10*np.log10(max(T_through[idx],1e-20)):.2f} dB")
print(f"  Cross output   T = {T_cross[idx]:.4f}  IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total          T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  Coupling ratio   = {T_cross[idx]/max(T_total[idx],1e-20):.3f}  "
      f"(0 = fully through, 1 = fully cross)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat = SMatrix(wavelengths=wavelengths)
smat.set(coeff_through, port="through")
smat.set(coeff_cross, port="cross")
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix saved (coupling_length={L_mid:.2f} µm slice) → {output_dir}/smatrix_hybrid.pkl")

np.savez(
    os.path.join(output_dir, "length_sweep.npz"),
    coupling_lengths=coupling_lengths, wavelengths=wavelengths,
    T_through=T_through_grid, T_cross=T_cross_grid,
)
print(f"Length-sweep grid saved → {output_dir}/length_sweep.npz")

# -- Wavelength-domain plot, at the reference coupling_length --------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_through, label="Through", color="steelblue")
ax.plot(wavelengths, T_cross,   label="Cross",   color="darkorange", linestyle="--")
ax.plot(wavelengths, T_total,   label="Total",   color="k", lw=1.5)
ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title(f"Transmission vs Wavelength (coupling_length={L_mid:.1f} µm)")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)

ax = axes[1]
IL_through = 10 * np.log10(np.clip(T_through, 1e-20, None))
IL_cross   = 10 * np.log10(np.clip(T_cross,   1e-20, None))
IL_total   = 10 * np.log10(np.clip(T_total,   1e-20, None))
ax.plot(wavelengths, IL_through, label="Through", color="steelblue")
ax.plot(wavelengths, IL_cross,   label="Cross",   color="darkorange", linestyle="--")
ax.plot(wavelengths, IL_total,   label="Total",   color="k", lw=1.5)
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title(f"Insertion Loss vs Wavelength (coupling_length={L_mid:.1f} µm)")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle("Directional Coupler — FDTD/EME Hybrid Simulation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Wavelength-domain results plot saved → {output_dir}/hybrid_results.png")

# -- Coupling-length sweep plot, mean/min/max across the wavelength sweep --
fig, ax = plt.subplots(figsize=(7, 4.5))

for grid, label, color in [
    (T_through_grid, "Through", "steelblue"),
    (T_cross_grid,   "Cross",   "darkorange"),
    (T_total_grid,   "Total",   "k"),
]:
    mean, lo, hi = grid.mean(axis=1), grid.min(axis=1), grid.max(axis=1)
    ax.plot(coupling_lengths, mean, label=label, color=color,
            lw=1.5 if label == "Total" else 1.2,
            linestyle="--" if label == "Cross" else "-")
    ax.fill_between(coupling_lengths, lo, hi, color=color, alpha=0.15, lw=0)

ax.axvline(L_mid, color="gray", linestyle=":", lw=0.8, label=f"Plotted length ({L_mid:.1f} µm)")
ax.set_xlabel("Coupling length (µm)")
ax.set_ylabel("Transmission (mean over wavelength, shaded = min–max)")
ax.set_title("Directional Coupler — Power Transfer vs Coupling Length")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "coupling_length_sweep.png"), dpi=300)
plt.close("all")
print(f"Coupling-length sweep plot saved → {output_dir}/coupling_length_sweep.png")
