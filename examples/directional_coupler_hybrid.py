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

HybridSolver itself has no idea this is a directional coupler; it only
knows each solver implements solve_section() and combines the resulting
per-section modal transfer matrices by matrix-chain multiplication. The
coupler-specific geometry (S-bends, coupling region) lives entirely in
DirectionalCoupler's section_prisms/section_gds/width_at/center_at methods.

Run:
    python /path/to/examples/directional_coupler_hybrid.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from feces import DirectionalCoupler, EMESolver, FDTDSolver, HybridSolver, SOI_220, plot_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_output", "directional_coupler_hybrid"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Geometry — all device dimensions as variables (220 nm SOI)
# ---------------------------------------------------------------------------
dc = DirectionalCoupler(
    waveguide_width=0.45,    # single-mode strip width (µm)
    gap=0.20,                # coupling-region edge-to-edge gap (µm)
    coupling_length=8.0,     # straight parallel interaction length (µm)
    bend_length=12.0,        # horizontal extent of each S-bend (µm)
    waveguide_offset=3.0,    # far-end centre-to-centre port pitch (µm)
    bend_type="euler",       # "euler" (gdsfactory bend_s) or "circular" (matched arcs)
    waveguide_height=0.22,   # Si layer thickness (µm) — SOI_220
)
print(f"Directional coupler total length: {dc.total_length:.2f} µm "
      f"(coupling pitch {dc.coupling_pitch:.3f} µm, far pitch {dc.waveguide_offset:.2f} µm)")

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
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholder length; solve_section() derives its own
    oxide_width=6.0, oxide_height=1.5, res=40.0, material_stack=SOI_220,
)
output_solver = FDTDSolver(**_fdtd_common)

# ---------------------------------------------------------------------------
# 3. Build HybridSolver — 3-section chain: FDTD bend | EME coupling | FDTD bend.
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
# 4. Geometry plot — shows which solver covers which section
# ---------------------------------------------------------------------------
plot_geometry(
    dc, solver=hybrid,
    show_hybrid_sections=True,
    show_pml=True,
    title="Directional Coupler — Hybrid Solver Sections",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
)
print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")

# ---------------------------------------------------------------------------
# 5. Run hybrid simulation — once per output port. output_idx selects which
#    arm the *output* section's S-bend follows (0 = bottom = "through",
#    1 = top = "cross"); input_idx stays at its default (0 = the excited
#    "in" arm) for both calls, so only the output section's kwargs differ.
#    A shared cache means the input/body sections are only ever simulated
#    once, not once per output port.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(1.50, 1.60, 11)
print(f"\nHybrid sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

section_cache = {}
smat_through = hybrid.solve(dc, wavelengths, port="through", cache=section_cache,
                            geometry_kwargs=[{}, {}, {"output_idx": 0}])
smat_cross = hybrid.solve(dc, wavelengths, port="cross", cache=section_cache,
                           geometry_kwargs=[{}, {}, {"output_idx": 1}])

T_through = smat_through.transmission(port="through")
T_cross   = smat_cross.transmission(port="cross")
T_total   = T_through + T_cross

wl_c = wavelengths[len(wavelengths) // 2]
idx  = np.argmin(np.abs(wavelengths - wl_c))

print(f"\nResults at λ = {wl_c:.3f} µm:")
print(f"  Through output T = {T_through[idx]:.4f}  IL = {10*np.log10(max(T_through[idx],1e-20)):.2f} dB")
print(f"  Cross output   T = {T_cross[idx]:.4f}  IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total          T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  Coupling ratio   = {T_cross[idx]/max(T_total[idx],1e-20):.3f}  "
      f"(0 = fully through, 1 = fully cross; tune coupling_length/gap to hit 0.5)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat = smat_through
smat._data.update(smat_cross._data)
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix_hybrid.pkl")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_through, label="Through", color="steelblue")
ax.plot(wavelengths, T_cross,   label="Cross",   color="darkorange", linestyle="--")
ax.plot(wavelengths, T_total,   label="Total",   color="k", lw=1.5)
ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title("Hybrid — Transmission vs Wavelength")
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
ax.set_title("Hybrid — Insertion Loss vs Wavelength")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle("Directional Coupler — FDTD/EME Hybrid Simulation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/hybrid_results.png")
