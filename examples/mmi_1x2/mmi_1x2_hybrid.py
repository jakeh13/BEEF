"""
Tutorial: 1×2 MMI — FDTD/EME hybrid simulation (FECES method).

HybridSolver is a generic N-section cascade: any solver type, in any order,
any count >= 2. Here we chain three sections —

    [input section]  ->  [body section]  ->  [output section]

— each independently configurable as FDTD or EME via the METHOD switches
below. HybridSolver itself has no idea any of these are FDTD/EME or that the
device is an MMI; it only knows each solver implements solve_section() and
combines the resulting per-section modal transfer matrices by matrix-chain
multiplication (the direct generalization of

    alpha_full[m] = sum_n alpha_in[m,n] * alpha_mm[m,n] * alpha_out[m,n]

to non-diagonal per-section matrices and an arbitrary number of sections).
The MMI-specific geometry (tapers, body, output branches) lives entirely in
MMI's section_prisms/section_gds/width_at/center_at methods.

Run:
    python /path/to/examples/mmi_1x2/mmi_1x2_hybrid.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from feces import MMI, EMESolver, FDTDSolver, HybridSolver, SOI_160, plot_hybrid_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_output", "mmi_1x2", "mmi_1x2_hybrid"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by both the EME and FDTD sub-solvers
WL_MIN = 1.26       # wavelength sweep start (µm)
WL_MAX = 1.36       # wavelength sweep end (µm)
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
    n_inputs=1,             # 1x2 splitter
    n_outputs=2,
    mm_length=4.5,
    mm_width=2.0,
    waveguide_width=0.35,
    taper_length=1.0,
    taper_width=0.8,
    waveguide_offset=1.0,
    waveguide_height=0.16,
)
print(f"MMI total length: {mmi.total_length:.2f} µm")

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
    three_d=True,       # full 3-D FDTD: real 160 nm Si core in oxide, PML on all six faces
    material_stack=SOI_160,
    decay_tol=1e-11,
)

N_MODES = 12   # modes carried at the two input/output <-> body interfaces

# The output section's single-mode port (the physical taper tip) sits only
# waveguide_offset=1.0 µm from the *other* output branch, which is drawn in
# the same section too (both branches are always physically present).
# half_domain_aperture confines the source/norm aperture to the near half of
# the domain (from the shared y=0 symmetry axis out to the near PML edge),
# so it structurally excludes the opposite branch regardless of how close
# the branches are. The real input port (single waveguide, no neighbor) is
# left untouched — its aperture keeps spanning the full sim_width.


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
            oxide_width=5.0, oxide_height=1.5, res=RESOLUTION, material_stack=SOI_160,
            mode_solver="meep",  # matches FDTDSolver's own eigenmode basis in this hybrid chain
        )
    raise ValueError(f"unknown method {method!r}, expected 'fdtd' or 'eme'")


# Input section keeps the default (full sim_width) aperture — a single
# waveguide port has no neighbor to avoid. Only the output section, whose
# two branches sit close together, needs the narrower aperture.
input_solver = make_section_solver(INPUT_METHOD)
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

# ---------------------------------------------------------------------------
# 4. Geometry plot — shows which solver covers which section
# ---------------------------------------------------------------------------
plot_hybrid_geometry(
    mmi, solver=hybrid,
    in_modes=1, out_modes=1,
    title="1×2 MMI — Hybrid Solver Sections",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True,
    sub_region_show_eme_detail=True,
    # One set of sub-region plots per output port — matches the two
    # hybrid.solve() calls below (geometry_kwargs={2: {"output_idx": ..}},
    # the sparse per-section form: only section 2's kwargs are named, the
    # input/body sections implicitly get {}): only the output section's
    # kwargs actually vary, so the input/body sections are each rendered
    # once (not once per output port).
    sub_region_configs=[
        ("out_top", {2: {"output_idx": 1}}),
        ("out_bot", {2: {"output_idx": 0}}),
    ],
)
print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")
print(f"Sub-region plots saved → {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run hybrid simulation — once per output port (output_idx selects the
#    output branch of the MMI geometry: output_ys = [-wo, +wo], so index 0
#    is the bottom port and index 1 is the top port). output_idx is passed
#    ONLY to the output section's kwargs (the input/body sections don't
#    branch, so their kwargs are identical across both calls) and a shared
#    cache is passed to both solve() calls, so the input/body sections are
#    only ever simulated once — not once per output port.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nHybrid sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

section_cache = {}
smat_top = hybrid.solve(mmi, wavelengths, port="out_top", cache=section_cache,
                         geometry_kwargs={2: {"output_idx": 1}})
smat_bot = hybrid.solve(mmi, wavelengths, port="out_bot", cache=section_cache,
                         geometry_kwargs={2: {"output_idx": 0}})

T_top   = smat_top.transmission(port="out_top")
T_bot   = smat_bot.transmission(port="out_bot")
T_total = T_top + T_bot

wl_c = wavelengths[len(wavelengths) // 2]
idx  = np.argmin(np.abs(wavelengths - wl_c))

print(f"\nResults at λ = {wl_c:.3f} µm:")
print(f"  Top output    T = {T_top[idx]:.4f}  IL = {10*np.log10(max(T_top[idx],1e-20)):.2f} dB")
print(f"  Bottom output T = {T_bot[idx]:.4f}  IL = {10*np.log10(max(T_bot[idx],1e-20)):.2f} dB")
print(f"  Total         T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  (ideal 50/50 split: T_each=0.5, IL_total=0.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat = smat_top
smat._data.update(smat_bot._data)
smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix_hybrid.pkl")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_top,   label="Top output",  color="steelblue")
ax.plot(wavelengths, T_bot,   label="Bot output",  color="darkorange", linestyle="--")
ax.plot(wavelengths, T_total, label="Total",        color="k", lw=1.5)
ax.axhline(0.5, color="gray", linestyle=":", lw=0.8, label="Ideal per-port")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title("Hybrid — Transmission vs Wavelength")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)

ax = axes[1]
IL_top = 10 * np.log10(np.clip(T_top,   1e-20, None))
IL_bot = 10 * np.log10(np.clip(T_bot,   1e-20, None))
IL_tot = 10 * np.log10(np.clip(T_total, 1e-20, None))
ax.plot(wavelengths, IL_top, label="Top output",  color="steelblue")
ax.plot(wavelengths, IL_bot, label="Bot output",  color="darkorange", linestyle="--")
ax.plot(wavelengths, IL_tot, label="Total",        color="k", lw=1.5)
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title("Hybrid — Insertion Loss vs Wavelength")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle("1×2 MMI — FDTD/EME Hybrid Simulation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/hybrid_results.png")
