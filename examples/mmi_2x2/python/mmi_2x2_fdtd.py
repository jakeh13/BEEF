"""
Tutorial: simulate a 2x2 MMI splitter with the BEEF library using FDTD.

Geometry:
    - Two single-mode input waveguides → MMI body → two output waveguides
    - SOI platform: 220 nm Si on SiO2
    - Dimensions from a typical fabricated C-band 2x2 device (mm_length = 28 µm)

Method: FDTD via meep with GDS prism geometry and EigenModeSource.
One input port (the top one) is excited; multimode cross-section monitors on
each output half capture coupling into all guided modes. A 2x2 MMI at its
first self-imaging length is a 50/50 splitter: exciting one input sends ~50%
to the output on the same side ("bar") and ~50% to the opposite side ("cross"),
with a 90° relative phase.

Both inputs sit close together (waveguide_offset = 1.5 µm), so the input
source is confined to a single waveguide via MMI's neighbor-aware aperture
(use_neighbor_aperture=True) — otherwise a full-width source would launch into
both inputs at once.

For the fastest result, use EME (mmi_2x2_eme.py). FDTD is the reference for
validation and for devices where EME is inaccurate.

Run from any directory:
    python /path/to/examples/mmi_2x2/python/mmi_2x2_fdtd.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, FDTDSolver, SOI_220, plot_fdtd_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_2x2", "mmi_2x2_fdtd"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by the FDTD solver
WL_MIN = 1.530      # wavelength sweep start (µm)
WL_MAX = 1.565      # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

# ---------------------------------------------------------------------------
# 1. Define geometry
# ---------------------------------------------------------------------------
mmi = MMI(
    n_inputs=2,             # 2x2 splitter
    n_outputs=2,
    mm_length=28.0,         # MMI body length (µm) — first 2x2 self-imaging length
    mm_width=2.8,           # MMI body width (µm)
    waveguide_width=0.45,   # single-mode access waveguide width (µm)
    taper_length=4.0,       # input/output taper length (µm)
    taper_width=1.3,        # taper width at MMI junction (µm)
    waveguide_offset=1.5,   # center-to-center port spacing (µm)
    waveguide_height=0.22,  # Si layer thickness (µm) — SOI_220
    use_neighbor_aperture=True,  # confine the input source to one of the two inputs
)

print(f"MMI total length: {mmi.total_length:.2f} µm")
print(f"Input ports at y = {mmi.input_ys} µm; outputs at y = {mmi.output_ys} µm")

# ---------------------------------------------------------------------------
# 2. Configure FDTD solver
# ---------------------------------------------------------------------------
solver = FDTDSolver(
    res=RESOLUTION,            # mesh resolution (pixels/µm)
    sim_width=5.0,             # simulation y-span (µm) — clears mm_width=2.8 plus margin
    sim_height=1.5,            # slab z-height reference (µm)
    dpml=0.5,                  # PML thickness, all 4 sides (µm)
    source_margin=0.0,         # source x-offset from left PML inner edge (µm)
    monitor_margin=0.25,       # monitor offset from PML inner edges, both sides (µm)
    three_d=True,       # full 3-D FDTD: real 220 nm Si core in oxide, PML on all six faces
    num_input_modes=1,         # input modes for normalization
    num_output_modes=1,        # modes extracted per output waveguide
    # Each output monitor spans one half of the domain (shared y=0 axis to the
    # near PML edge), so the two output waveguides (±0.75 µm) are captured in
    # isolation rather than as a two-guide supermode.
    output_monitor_centers=[2.5/2, -2.5/2],   # out_0 = top half (bar), out_1 = bottom half (cross)
    output_monitor_widths=[2.5, 2.5],
    material_stack=SOI_220,    # 220 nm SOI stack
    decay_tol=1e-11,           # DFT convergence: lower = more accurate, slower
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
plot_fdtd_geometry(
    mmi, solver=solver,
    title="2×2 MMI Splitter — Top-down Geometry (FDTD)",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved → {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)  # µm, C-band

print(f"Sweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm")
print("(single broadband pulse; DFT monitors extract all wavelengths in one sim)")

# Excite the top input port; measure coupling to the fundamental mode at each
# output half.
smat = solver.solve(mmi, wavelengths, port="top", input_mode=1, save_dir=save_dir)
smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength

# ---------------------------------------------------------------------------
# 5. Print summary
# ---------------------------------------------------------------------------
wl_center = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_center))

# "out_0" / "out_1" match the order of output_monitor_centers: [+1.25, -1.25].
# Excited input is the top port, so out_0 (top) is the bar output and out_1
# (bottom) is the cross output.
T_bar   = smat.transmission(port="out_0", in_mode=1, out_mode=1)  # [M]
T_cross = smat.transmission(port="out_1", in_mode=1, out_mode=1)  # [M]
T_total = T_bar + T_cross

print(f"\nResults at λ = {wl_center:.3f} µm:")
print(f"  Bar output (top)     : {T_bar[idx]:.4f}")
print(f"  Cross output (bottom): {T_cross[idx]:.4f}")
print(f"  Total (sum)          : {T_total[idx]:.4f}")
print(f"  Insertion loss       : {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  (ideal 50/50 split → T_each=0.5, IL_total=0.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix.pkl")

# Transmission vs wavelength: both outputs + total
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(wavelengths, T_bar,   label="Bar output (top, mode 1)")
ax.plot(wavelengths, T_cross, label="Cross output (bottom, mode 1)", linestyle="--")
ax.plot(wavelengths, T_total, label="Total (bar + cross)", color="k", lw=1.5)
ax.axhline(0.5, color="gray", linestyle=":", lw=0.8, label="Ideal per-port (0.5)")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title("2×2 MMI — FDTD Transmission vs Wavelength")
ax.legend()
ax.grid(True)
ax.set_ylim(0, 1.05)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "transmission.png"), dpi=300)
plt.close("all")
print(f"Transmission plot saved → {output_dir}/transmission.png")
