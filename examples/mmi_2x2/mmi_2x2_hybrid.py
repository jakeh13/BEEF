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

Run:
    python /path/to/examples/mmi_2x2/mmi_2x2_hybrid.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from feces import MMI, EMESolver, FDTDSolver, HybridSolver, SOI_220, plot_hybrid_geometry

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

hybrid = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
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
# 5. Run hybrid simulation — once per output port. input_idx=1 (top input)
#    is fixed across both calls, so with a shared cache the input and body
#    sections are only ever simulated once; only the output section is
#    re-run for the other output_idx.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nHybrid sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

section_cache = {}
smat_bar = hybrid.solve(mmi, wavelengths, port="bar", cache=section_cache,
                        geometry_kwargs=[IN_KW, {}, {"output_idx": 1}])
smat_cross = hybrid.solve(mmi, wavelengths, port="cross", cache=section_cache,
                          geometry_kwargs=[IN_KW, {}, {"output_idx": 0}])

T_bar   = smat_bar.transmission(port="bar")
T_cross = smat_cross.transmission(port="cross")
T_total = T_bar + T_cross

wl_c = wavelengths[len(wavelengths) // 2]
idx  = np.argmin(np.abs(wavelengths - wl_c))

print(f"\nResults at λ = {wl_c:.3f} µm:")
print(f"  Bar output   T = {T_bar[idx]:.4f}  IL = {10*np.log10(max(T_bar[idx],1e-20)):.2f} dB")
print(f"  Cross output T = {T_cross[idx]:.4f}  IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total        T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  (ideal 50/50 split: T_each=0.5, IL_total=0.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat = smat_bar
smat._data.update(smat_cross._data)
smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix_hybrid.pkl")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_bar,   label="Bar output",   color="steelblue")
ax.plot(wavelengths, T_cross, label="Cross output", color="darkorange", linestyle="--")
ax.plot(wavelengths, T_total, label="Total",         color="k", lw=1.5)
ax.axhline(0.5, color="gray", linestyle=":", lw=0.8, label="Ideal per-port")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title("Hybrid — Transmission vs Wavelength")
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
ax.set_title("Hybrid — Insertion Loss vs Wavelength")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle("2×2 MMI — FDTD/EME Hybrid Simulation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/hybrid_results.png")
