"""
Tutorial: simulate a 1x2 MMI splitter with the FECES library using EME.

Geometry:
    - Single-mode input waveguide → MMI body → two output waveguides
    - SOI platform: 160 nm Si on SiO2
    - Dimensions from typical fabricated device (mm_length = 4.5 µm)

Method: Eigenmode Expansion (EME) via meow.  Fast, no FDTD required.
The EME models the path from the input to the top output waveguide.
For a balanced 1x2 splitter each output carries ~50% of the input power
(ideal insertion loss: −3 dB per port).

Run from any directory:
    python /path/to/examples/mmi_1x2/mmi_1x2_eme.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")           # headless — swap for "Qt5Agg" for interactive

from feces import MMI, EMESolver, SOI_160, plot_eme_geometry, plot_performance

# Output to a top-level example_output/ folder outside the repo
_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_output", "mmi_1x2", "mmi_1x2_eme"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by the EME solver
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
# 2. Configure EME solver
# ---------------------------------------------------------------------------
solver = EMESolver(
    # 5 regions: wg_in | taper_in | MMI body | taper_out | wg_out
    num_modes=[1,   5,   12,  5,   1  ],  # modes per region
    num_cells=[1,   8,  1,   8,  1  ],  # cells per region (tapers subdivided)
    lengths=  [0.5, 1.0, 4.5, 1.0, 0.5], # region lengths (µm)
    oxide_width=5.0,                       # cross-section mesh y-span (µm)
    oxide_height=1.5,                      # cross-section mesh z-span (µm)
    res=RESOLUTION,                        # mesh pts/µm
    material_stack=SOI_160,                # 160 nm SOI stack
    mode_solver="tidy3d",                  # tidy3d-backed mode solver (no FDTD section to stay consistent with)
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
os.makedirs(output_dir, exist_ok=True)
plot_eme_geometry(
    mmi, solver=solver,
    region_labels=["wg in", "taper in", "MMI body", "taper out", "wg out"],
    title="1×2 MMI Splitter — Top-down Geometry (top output only)",
    output_path=os.path.join(output_dir, "geometry.png"),
    output_idx=1,  # top output branch only — matches what solver.solve() below monitors
)
print(f"Geometry plot saved → {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)  # µm, O-band

print(f"Sweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm")

# "out" = top output port. output_idx=1 restricts the EME cross-sections to
# the top output branch only (output_ys[1] = +waveguide_offset/2), so the
# solved mode is the isolated waveguide mode rather than a supermode of both
# outputs. S[0,1] coeff gives coupling from input mode 1 into that waveguide.
smat = solver.solve(mmi, wavelengths, ports=["out"], output_idx=1)

# ---------------------------------------------------------------------------
# 5. Print summary
# ---------------------------------------------------------------------------
T = smat.transmission("out")
IL = smat.insertion_loss_db("out")
wl_center = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_center))

print(f"\nResults at λ = {wl_center:.3f} µm:")
print(f"  Transmission to top output : {T[idx]:.4f}")
print(f"  Insertion loss             : {IL[idx]:.2f} dB")
print(f"  (ideal 50/50 split → T=0.5, IL=−3.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------------------------

smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix.pkl")

# Plot insertion loss and transmission vs wavelength
plot_performance(smat, output_dir=output_dir, title="1×2 MMI Splitter (EME)")
print(f"Performance plot saved → {output_dir}/performance.png")
