"""
Tutorial: simulate a 1x2 MMI splitter with the BEEF library using FDTD.

Geometry:
    - Single-mode input waveguide → MMI body → two output waveguides
    - SOI platform: 160 nm Si on SiO2
    - Dimensions from typical fabricated device (mm_length = 4.5 µm)

Method: FDTD via meep with GDS prism geometry and EigenModeSource.
The full device is simulated in a single domain. A multimode cross-section
monitor is placed near the output to capture coupling into all guided modes.
Transmission to the dominant output mode (mode 1) gives the port IL.

For the fastest result, use EME (mmi_1x2_eme.py). FDTD is the reference
for validation and for devices where EME is inaccurate (high-contrast bends,
scattering geometries).

Run from any directory:
    python /path/to/examples/mmi_1x2/python/mmi_1x2_fdtd.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, FDTDSolver, SOI_160, plot_fdtd_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_1x2", "mmi_1x2_fdtd"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by the FDTD solver
WL_MIN = 1.26       # wavelength sweep start (µm)
WL_MAX = 1.36       # wavelength sweep end (µm)
N_WL = 11            # number of wavelength points

# ---------------------------------------------------------------------------
# 1. Define geometry
# ---------------------------------------------------------------------------
mmi = MMI(
    n_inputs=1,             # 1x2 splitter
    n_outputs=2,
    mm_length=4.5,          # MMI body length (µm)
    mm_width=2.0,           # MMI body width (µm)
    waveguide_width=0.35,   # single-mode access waveguide width (µm)
    taper_length=1.0,       # input/output taper length (µm)
    taper_width=0.8,        # taper width at MMI junction (µm)
    waveguide_offset=1.0,   # center-to-center output spacing (µm)
    waveguide_height=0.16,  # Si layer thickness (µm) — SOI_160
)

print(f"MMI total length: {mmi.total_length:.2f} µm")

# ---------------------------------------------------------------------------
# 2. Configure FDTD solver
# ---------------------------------------------------------------------------
solver = FDTDSolver(
    res=RESOLUTION,            # mesh resolution (pixels/µm)
    sim_width=5.0,             # simulation y-span (µm)
    sim_height=1.5,            # slab z-height reference (µm)
    dpml=0.5,                  # PML thickness, all 4 sides (µm)
    source_margin=0.0,         # source x-offset from left PML inner edge (µm)
    monitor_margin=0.25,       # monitor offset from PML inner edges, both sides (µm)
    three_d=True,       # full 3-D FDTD: real 160 nm Si core in oxide, PML on all six faces
    num_input_modes=1,         # input modes for normalization
    num_output_modes=1,        # modes extracted per output waveguide
    # Output monitor centers at ±sim_width/4 = ±1.25 µm (half-domain aperture,
    # each spanning from the shared y=0 axis to the near PML edge)
    output_monitor_centers=[2.5/2, -2.5/2],
    output_monitor_widths=[2.5, 2.5],
    material_stack=SOI_160,    # 160 nm SOI stack
    decay_tol=1e-11,           # DFT convergence: lower = more accurate, slower
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
plot_fdtd_geometry(
    mmi, solver=solver,
    title="1×2 MMI Splitter — Top-down Geometry (FDTD)",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved → {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)  # µm, O-band

print(f"Sweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm")
print("(single broadband pulse; DFT monitors extract all wavelengths in one sim)")

# Source from top input port; measure coupling to each mode at output monitor
smat = solver.solve(mmi, wavelengths, port="top", input_mode=1, save_dir=save_dir)
smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength

# ---------------------------------------------------------------------------
# 5. Print summary
# ---------------------------------------------------------------------------
wl_center = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_center))

# Power in fundamental mode at each output waveguide
# "out_0" / "out_1" match the order of output_monitor_centers: [+0.5, -0.5]
T_top = smat.transmission(port="out_0", in_mode=1, out_mode=1)  # [M]
T_bot = smat.transmission(port="out_1", in_mode=1, out_mode=1)  # [M]
T_total = T_top + T_bot  # insertion loss = sum of both outputs

print(f"\nResults at λ = {wl_center:.3f} µm:")
print(f"  Top output power    : {T_top[idx]:.4f}")
print(f"  Bottom output power : {T_bot[idx]:.4f}")
print(f"  Total (sum)         : {T_total[idx]:.4f}")
print(f"  Insertion loss      : {10*np.log10(T_total[idx]):.2f} dB")
print(f"  (ideal 50/50 split → T_each=0.5, IL_total=0.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix.pkl")

# Transmission vs wavelength: both outputs + total
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(wavelengths, T_top,   label="Top output (mode 1)")
ax.plot(wavelengths, T_bot,   label="Bottom output (mode 1)", linestyle="--")
ax.plot(wavelengths, T_total, label="Total (top + bot)", color="k", lw=1.5)
ax.axhline(0.5, color="gray", linestyle=":", lw=0.8, label="Ideal per-port (0.5)")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title("1×2 MMI — FDTD Transmission vs Wavelength")
ax.legend()
ax.grid(True)
ax.set_ylim(0, 1.05)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "transmission.png"), dpi=300)
plt.close("all")
print(f"Transmission plot saved → {output_dir}/transmission.png")
