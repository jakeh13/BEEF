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
    python /path/to/examples/mmi_1x2_eme.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")           # headless — swap for "Qt5Agg" for interactive
import matplotlib.pyplot as plt

from feces import MMI1x2, EMESolver, SOI_160, plot_geometry, plot_performance

# Output to a top-level example_output/ folder outside the repo
_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_output", "mmi_1x2_eme"))

# ---------------------------------------------------------------------------
# 1. Define geometry
# ---------------------------------------------------------------------------
mmi = MMI1x2(
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
    num_modes=[1,   7,   15,  7,   1  ],  # modes per region
    num_cells=[1,   12,  1,   12,  1  ],  # cells per region (tapers subdivided)
    lengths=  [0.5, 1.0, 4.5, 1.0, 0.5], # region lengths (µm)
    oxide_width=5.0,                       # cross-section mesh y-span (µm)
    oxide_height=1.5,                      # cross-section mesh z-span (µm)
    res=20.0,                              # mesh pts/µm
    material_stack=SOI_160,                # 160 nm SOI stack
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
os.makedirs(output_dir, exist_ok=True)
plot_geometry(
    mmi, solver=solver,
    show_eme_regions=True,
    region_labels=["wg in", "taper in", "MMI body", "taper out", "wg out"],
    title="1×2 MMI Splitter — Top-down Geometry",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved → {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep
# ---------------------------------------------------------------------------
wavelengths = np.linspace(1.26, 1.36, 5)  # µm, O-band

print(f"Sweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm")

# "out" = top output port. GDS draws one output waveguide;
# S[0,1] coeff gives coupling from input mode 1 into that waveguide.
smat = solver.solve(mmi, wavelengths, ports=["out"])

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

# Custom IL plot vs ideal −3 dB
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(wavelengths, IL, label="Top output")
ax.axhline(-3.0, color="k", linestyle="--", lw=0.8, label="Ideal −3 dB")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title("1×2 MMI — Insertion Loss vs Wavelength")
ax.legend()
ax.grid(True)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "insertion_loss.png"), dpi=300)
plt.close("all")
print(f"IL plot saved           → {output_dir}/insertion_loss.png")
