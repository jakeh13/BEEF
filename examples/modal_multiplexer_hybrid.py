"""
Tutorial: modal (mode-division) multiplexer — FDTD/EME hybrid simulation (FECES method).

Geometry (220 nm SOI, C-band) — an *asymmetric* directional coupler:

    "access"  ___                                      ____  "access_out"
    (narrow,     \\                                    /      (TE0 only,
     TE0 only)     \\__________________________________/        leftover)
                   /                                    \\
    "bus"    ___  /                                      \\___  "bus_out"
    (wide,                                                      (TE0 + TE1,
     TE0+TE1)                                                    multiplexed)

    |<- bend_length ->|<--- coupling_length --->|<- bend_length ->|

A DirectionalCoupler (see directional_coupler_hybrid.py) has two *equal*-width
arms, so it only ever splits/combines power *within* a single mode — it's a
power multiplexer. This device instead gives the two arms two *different*
widths:

  - the "bus" arm is wide enough to guide both the fundamental mode (TE0)
    and the first higher-order mode (TE1),
  - the "access" arm is narrow enough to guide only TE0.

With the access arm's width chosen so its TE0 effective index sits close to
the bus arm's TE1 effective index (the phase-matching condition for
asymmetric directional couplers), coupling in the straight region becomes
*mode-selective*: power launched into the access arm's TE0 converts into the
bus arm's TE1, while power already on the bus arm's TE0 is phase-mismatched
with the access arm and mostly passes straight through undisturbed. That's
the multiplexing action this script characterizes — an independent signal
riding on the access waveguide gets added onto a *second mode* of the shared
bus waveguide, without disturbing whatever is already on its first mode.

BUS_WIDTH/ACCESS_WIDTH below are an approximate, literature-typical
phase-matched pair for 220 nm SOI at 1.55 µm — picked by hand, not derived
from a dedicated dispersion sweep. See directional_coupler_mode.py for the
style of a standalone meow mode-solve script; if the conversion-length sweep
plot below never shows TE1 conversion approaching ~1, that phase-matching
guess is the first thing to revisit.

As with the plain directional coupler, only **2 FDTD simulations + 1 EME
simulation total** are needed, thanks to the same two tricks:

    [input S-bend]  ->  [straight coupled pair]  ->  [output S-bend]
       FDTD (x2)               EME (x1)              reused, not resimulated

  - The S-bends (curved, non-adiabatic) are solved with FDTD; the straight
    coupled-pair region (uniform cross-section) is solved once via EME as an
    eigenmode expansion over N_MODES supermodes of the *asymmetric* coupled
    pair (need enough modes to resolve the hybridized bus(2-mode)+access
    (1-mode) system — see N_MODES below).
  - ModalMultiplexer has no left/right asymmetry (same bend_length, arm
    widths, waveguide_offset, bend_type on both sides), so — exactly as in
    directional_coupler_hybrid.py — the output section is never separately
    simulated: by reciprocity, the input section's transfer matrix doubles
    as the output section's transfer matrix once transposed.
  - Unlike the plain coupler, the "bus" port is *multi-mode* (TE0 + TE1) on
    both ends, so its input-side FDTD run launches both modes separately
    (n_modes_in=2) instead of just one — everything else about the
    section-solving machinery is unchanged.

Run:
    python /path/to/examples/modal_multiplexer_hybrid.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from feces import (
    ModalMultiplexer, EMESolver, FDTDSolver, HybridSolver, SMatrix, SOI_220, plot_hybrid_geometry,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_output", "modal_multiplexer_hybrid"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by both the EME and FDTD sub-solvers
WL_MIN = 1.50       # wavelength sweep start (µm)
WL_MAX = 1.60       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

# Coupling length used for the wavelength-domain ("hybrid_results.png") plot
# and the SMatrix that gets saved. Set to a specific value (µm) to pin it;
# leave as None to auto-pick the length (from coupling_lengths, defined
# below) that maximizes access->bus TE1 conversion efficiency.
PLOT_COUPLING_LENGTH = None   # µm, or None for auto (peak TE1 conversion)

# ---------------------------------------------------------------------------
# 1. Geometry — all device dimensions as variables (220 nm SOI).
#
#    coupling_lengths is the actual sweep; mux itself is built at the
#    sweep's middle value (bends/access-waveguide geometry don't depend on
#    coupling_length, so this single instance is equally valid as the
#    reference geometry for the FDTD sections at every length in the sweep).
# ---------------------------------------------------------------------------
coupling_lengths = np.linspace(0.0, 60.0, 201)   # µm, the swept coupling length
L_mid = coupling_lengths[len(coupling_lengths) // 2]

mux = ModalMultiplexer(
    bus_width=0.90,          # wide arm — guides TE0 + TE1 (µm)
    access_width=0.44,       # narrow arm — guides TE0 only, ~phase-matched to bus TE1 (µm)
    gap=0.20,                # coupling-region edge-to-edge gap (µm)
    coupling_length=L_mid,   # straight parallel interaction length (µm) — sweep reference value
    bend_length=12.0,        # horizontal extent of each S-bend (µm)
    waveguide_offset=3.0,    # far-end centre-to-centre port pitch (µm)
    bend_type="euler",       # "euler" (cubic-Bezier S-bend) or "circular" (matched arcs)
    waveguide_height=0.22,   # Si layer thickness (µm) — SOI_220
)
print(f"Modal multiplexer total length: {mux.total_length:.2f} µm "
      f"(bus width {mux.bus_width:.3f} µm, access width {mux.access_width:.3f} µm, "
      f"gap {mux.gap:.3f} µm, far pitch {mux.waveguide_offset:.2f} µm)")
print(f"Coupling length sweep: {coupling_lengths[0]:.1f} – {coupling_lengths[-1]:.1f} µm "
      f"({len(coupling_lengths)} pts, build reference length = {L_mid:.2f} µm)")

# ---------------------------------------------------------------------------
# 2. Configure sub-solvers — one per section.
#
#    input section  : input access waveguide + input S-bend, ending exactly
#                      at the S-bend/coupling-region seam
#    body section   : propagates N_MODES supermodes through the full
#                      straight coupled-pair region (length = coupling_length)
#    output section : the mirror image of the input section — never
#                      separately simulated (see module docstring)
# ---------------------------------------------------------------------------
N_MODES = 4   # supermodes carried at the two body interfaces — bus(2) + access(1)
              # hybridize into >= 3 supermodes; 4 gives margin (cf. DC's N_MODES=6 for a 2-mode system)

_fdtd_common = dict(
    res=RESOLUTION,
    sim_width=6.0,        # >= waveguide_offset + margin, so both far-pitch arms sit inside PML clearance
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    # A full-sim_width source/norm aperture at a true I/O port would reach
    # past the shared y=0 axis to the *other* arm (only waveguide_offset
    # away), since both arms are always physically present in a section -
    # corrupting the isolated single-arm excitation. half_domain_aperture
    # confines the aperture to the near half of the domain (from the
    # symmetry axis to the near PML edge) instead. Only applies when the
    # excited side is genuinely single-mode (see FDTDSolver docstring) — the
    # bus side's 2-mode excitation always keeps the full aperture.
    half_domain_aperture=True,
    material_stack=SOI_220,
    decay_tol=1e-11,
)

# One FDTDSolver instance covers both the input and output sections — they're
# mirror images of the same physical structure (see module docstring), so
# there's no need for two separately-configured objects.
fdtd_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholder length; solve_section_length_sweep() sweeps its own
    oxide_width=6.0, oxide_height=1.5, res=RESOLUTION, material_stack=SOI_220,
    mode_solver="meep",  # matches FDTDSolver's own eigenmode basis in this hybrid chain
)

# ---------------------------------------------------------------------------
# 3. Build HybridSolver — 3-section chain: FDTD bend | EME coupling | FDTD bend.
#    (Used below only for the geometry/section preview plots; the actual
#    numeric sweep bypasses HybridSolver.solve() — see step 5.)
# ---------------------------------------------------------------------------
half_body = mux.coupling_length / 2
half_total = mux.total_length / 2

boundaries = [
    -half_total,   # physical input port (access waveguide end)
    -half_body,    # input section <-> body section interface (S-bend/coupling-region seam)
    +half_body,    # body section <-> output section interface (S-bend/coupling-region seam)
    +half_total,   # physical output port (access waveguide end)
]
interface_modes = [N_MODES, N_MODES]

hybrid = HybridSolver(
    solvers=[fdtd_solver, body_solver, fdtd_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
)

# ---------------------------------------------------------------------------
# 4. Geometry plot — shows which solver covers which section, at mux's
#    (middle-of-sweep) coupling_length.
# ---------------------------------------------------------------------------
plot_hybrid_geometry(
    mux, solver=hybrid,
    in_modes=1, out_modes=N_MODES,
    title=f"Modal Multiplexer — Hybrid Solver Sections (coupling_length={L_mid:.1f} µm)",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True,
    sub_region_show_eme_detail=True,
    # One sub-region plot per simulation actually run in step 5: the input
    # section is solved once per input arm (input_idx=0 "bus", input_idx=1
    # "access"). sub_region_sections=[0, 1] skips the output section
    # (index 2), which is never separately simulated (see module docstring)
    # but is still visible in the full-device overview plot above.
    sub_region_configs=[
        ("bus_in", [{"input_idx": 0}, {}, {}]),
        ("access_in", [{"input_idx": 1}, {}, {}]),
    ],
    sub_region_sections=[0, 1],
)
print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")
print(f"Sub-region plots saved → {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run the sweep — 2 FDTD solves (input section, once per input arm) plus
#    one EME mode-solve per wavelength, combined with the coupling_length-
#    swept body matrix via the same matrix-chain contraction
#    HybridSolver.solve() uses internally (broadcast over the extra length
#    axis). The output section is never separately simulated: M_bus_in/
#    M_access_in double as the output-side transfer matrices once transposed
#    (see module docstring). Two excitation scenarios are evaluated at the
#    (multi-mode) bus_out port:
#      bus_in    -> bus_out   : pass-through check (want TE0 high, TE1 low)
#      access_in -> bus_out   : the multiplexing channel (want TE1 high, TE0 low)
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nWavelength sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

M_bus_in = fdtd_solver.solve_section(
    mux, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=2, n_modes_out=N_MODES,
    stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
    input_idx=0,
)  # [n_wl, N_MODES, 2] — bus arm (TE0, TE1) -> body modes
M_access_in = fdtd_solver.solve_section(
    mux, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
    input_idx=1,
)  # [n_wl, N_MODES, 1] — access arm (TE0) -> body modes

M_out_bus = M_bus_in.transpose(0, 2, 1)        # [n_wl, 2, N_MODES] — mirror-symmetry reuse, not a new simulation
M_out_access = M_access_in.transpose(0, 2, 1)  # [n_wl, 1, N_MODES]

# The FDTD input/output sections end exactly at the S-bend/coupling-region
# seam (boundaries[1]/[2] above capture none of the straight interior), so
# the EME body models the full nominal coupling_length directly.
body_lengths = coupling_lengths
M_body = body_solver.solve_section_length_sweep(
    mux, wavelengths, boundaries[1], boundaries[2], body_lengths,
    n_modes_in=N_MODES, n_modes_out=N_MODES,
)  # [n_L, n_wl, N_MODES, N_MODES]

# -- Scenario 1: bus_in TE0 -> bus_out (TE0 pass-through, TE0->TE1 crosstalk) --
tmp_bus = np.einsum("Lwij,wjk->Lwik", M_body, M_bus_in[:, :, 0:1])        # [n_L, n_wl, N_MODES, 1]
full_bus_from_bus = np.einsum("wij,Lwjk->Lwik", M_out_bus, tmp_bus)      # [n_L, n_wl, 2, 1]

# -- Scenario 2: access_in TE0 -> bus_out (the multiplexed channel) --------
tmp_access = np.einsum("Lwij,wjk->Lwik", M_body, M_access_in)             # [n_L, n_wl, N_MODES, 1]
full_bus_from_access = np.einsum("wij,Lwjk->Lwik", M_out_bus, tmp_access)  # [n_L, n_wl, 2, 1]

T_bus_pass_grid = np.abs(full_bus_from_bus[:, :, 0, 0]) ** 2       # [n_L, n_wl] — bus TE0 -> bus TE0
T_bus_crosstalk_grid = np.abs(full_bus_from_bus[:, :, 1, 0]) ** 2  # bus TE0 -> bus TE1 (undesired)
T_access_leak_grid = np.abs(full_bus_from_access[:, :, 0, 0]) ** 2  # access TE0 -> bus TE0 (undesired)
T_access_conv_grid = np.abs(full_bus_from_access[:, :, 1, 0]) ** 2  # access TE0 -> bus TE1 (the mux, headline metric)

# ---------------------------------------------------------------------------
# 5b. Choose the coupling length used for the wavelength-domain plot/SMatrix.
#     PLOT_COUPLING_LENGTH pins it explicitly (snapped to the nearest sampled
#     length); left as None, auto-pick the first coupling length at which
#     mean-over-wavelength access->bus TE1 conversion peaks (the "conversion
#     length", analogous to a plain coupler's 3 dB length).
# ---------------------------------------------------------------------------
if PLOT_COUPLING_LENGTH is not None:
    L_plot_idx = int(np.argmin(np.abs(coupling_lengths - PLOT_COUPLING_LENGTH)))
    print(f"Plot length: user-set {PLOT_COUPLING_LENGTH:.2f} µm -> "
          f"nearest sampled length {coupling_lengths[L_plot_idx]:.2f} µm")
else:
    conv_mean = T_access_conv_grid.mean(axis=1)
    d = np.diff(conv_mean)
    local_maxima = np.flatnonzero((d[:-1] > 0) & (d[1:] <= 0)) + 1
    L_plot_idx = int(local_maxima[0]) if local_maxima.size else int(np.argmax(conv_mean))
    print(f"Plot length: auto-picked {coupling_lengths[L_plot_idx]:.2f} µm "
          f"(peak TE1 conversion, efficiency = {conv_mean[L_plot_idx]:.3f})")
L_plot = coupling_lengths[L_plot_idx]

# Slice at the chosen coupling_length for the wavelength-domain view — keep
# the complex coefficients (not just power) for this slice so the saved
# SMatrix preserves phase, same as directional_coupler_hybrid.py.
coeff_bus_pass = full_bus_from_bus[L_plot_idx, :, 0, 0]
coeff_bus_crosstalk = full_bus_from_bus[L_plot_idx, :, 1, 0]
coeff_access_leak = full_bus_from_access[L_plot_idx, :, 0, 0]
coeff_access_conv = full_bus_from_access[L_plot_idx, :, 1, 0]
T_bus_pass = T_bus_pass_grid[L_plot_idx]
T_bus_crosstalk = T_bus_crosstalk_grid[L_plot_idx]
T_access_leak = T_access_leak_grid[L_plot_idx]
T_access_conv = T_access_conv_grid[L_plot_idx]

wl_c = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_c))

print(f"Results at coupling_length = {L_plot:.2f} µm, λ = {wl_c:.3f} µm:")
print(f"  bus_in  TE0 -> bus_out TE0 (pass-through) T = {T_bus_pass[idx]:.4f}  "
      f"IL = {10*np.log10(max(T_bus_pass[idx],1e-20)):.2f} dB")
print(f"  bus_in  TE0 -> bus_out TE1 (crosstalk)    T = {T_bus_crosstalk[idx]:.4f}  "
      f"IL = {10*np.log10(max(T_bus_crosstalk[idx],1e-20)):.2f} dB")
print(f"  access_in TE0 -> bus_out TE0 (leakage)    T = {T_access_leak[idx]:.4f}  "
      f"IL = {10*np.log10(max(T_access_leak[idx],1e-20)):.2f} dB")
print(f"  access_in TE0 -> bus_out TE1 (conversion) T = {T_access_conv[idx]:.4f}  "
      f"IL = {10*np.log10(max(T_access_conv[idx],1e-20)):.2f} dB")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat = SMatrix(wavelengths=wavelengths)
smat.set(coeff_bus_pass, port="bus_out_TE0__from_bus_in")
smat.set(coeff_bus_crosstalk, port="bus_out_TE1__from_bus_in")
smat.set(coeff_access_leak, port="bus_out_TE0__from_access_in")
smat.set(coeff_access_conv, port="bus_out_TE1__from_access_in")
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix saved (coupling_length={L_plot:.2f} µm slice) → {output_dir}/smatrix_hybrid.pkl")

np.savez(
    os.path.join(output_dir, "length_sweep.npz"),
    coupling_lengths=coupling_lengths, wavelengths=wavelengths,
    T_bus_pass=T_bus_pass_grid, T_bus_crosstalk=T_bus_crosstalk_grid,
    T_access_leak=T_access_leak_grid, T_access_conv=T_access_conv_grid,
)
print(f"Length-sweep grid saved → {output_dir}/length_sweep.npz")

# -- Wavelength-domain plot, at the reference coupling_length --------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_bus_pass, label="bus TE0 → bus TE0 (through)", color="steelblue")
ax.plot(wavelengths, T_bus_crosstalk, label="bus TE0 → bus TE1 (crosstalk)", color="steelblue", linestyle=":")
ax.plot(wavelengths, T_access_conv, label="access TE0 → bus TE1 (mux)", color="darkorange")
ax.plot(wavelengths, T_access_leak, label="access TE0 → bus TE0 (leakage)", color="darkorange", linestyle=":")
ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title(f"Transmission vs Wavelength (coupling_length={L_plot:.1f} µm)")
ax.legend(fontsize=7)
ax.grid(True)
ax.set_ylim(0, 1.05)

ax = axes[1]
IL_bus_pass = 10 * np.log10(np.clip(T_bus_pass, 1e-20, None))
IL_bus_crosstalk = 10 * np.log10(np.clip(T_bus_crosstalk, 1e-20, None))
IL_access_conv = 10 * np.log10(np.clip(T_access_conv, 1e-20, None))
IL_access_leak = 10 * np.log10(np.clip(T_access_leak, 1e-20, None))
ax.plot(wavelengths, IL_bus_pass, label="bus TE0 → bus TE0 (through)", color="steelblue")
ax.plot(wavelengths, IL_bus_crosstalk, label="bus TE0 → bus TE1 (crosstalk)", color="steelblue", linestyle=":")
ax.plot(wavelengths, IL_access_conv, label="access TE0 → bus TE1 (mux)", color="darkorange")
ax.plot(wavelengths, IL_access_leak, label="access TE0 → bus TE0 (leakage)", color="darkorange", linestyle=":")
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title(f"Insertion Loss vs Wavelength (coupling_length={L_plot:.1f} µm)")
ax.legend(fontsize=7)
ax.grid(True)

plt.suptitle("Modal Multiplexer — FDTD/EME Hybrid Simulation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Wavelength-domain results plot saved → {output_dir}/hybrid_results.png")

# -- Coupling-length sweep plot, mean/min/max across the wavelength sweep --
fig, ax = plt.subplots(figsize=(7, 4.5))

for grid, label, color, ls in [
    (T_access_conv_grid, "access TE0 → bus TE1 (mux)", "darkorange", "-"),
    (T_access_leak_grid, "access TE0 → bus TE0 (leakage)", "darkorange", ":"),
    (T_bus_pass_grid, "bus TE0 → bus TE0 (through)", "steelblue", "-"),
    (T_bus_crosstalk_grid, "bus TE0 → bus TE1 (crosstalk)", "steelblue", ":"),
]:
    mean, lo, hi = grid.mean(axis=1), grid.min(axis=1), grid.max(axis=1)
    ax.plot(coupling_lengths, mean, label=label, color=color, lw=1.3, linestyle=ls)
    ax.fill_between(coupling_lengths, lo, hi, color=color, alpha=0.12, lw=0)

_plot_len_note = "peak conversion, auto" if PLOT_COUPLING_LENGTH is None else "user-set"
ax.axvline(L_plot, color="gray", linestyle=":", lw=0.8,
           label=f"Plotted length ({L_plot:.1f} µm, {_plot_len_note})")
ax.set_xlabel("Coupling length (µm)")
ax.set_ylabel("Transmission (mean over wavelength, shaded = min–max)")
ax.set_title("Modal Multiplexer — Mode Conversion vs Coupling Length")
ax.legend(fontsize=7)
ax.grid(True)
ax.set_ylim(0, 1.05)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "coupling_length_sweep.png"), dpi=300)
plt.close("all")
print(f"Coupling-length sweep plot saved → {output_dir}/coupling_length_sweep.png")
