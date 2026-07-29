"""
Tutorial: simulate a 2x2 MMI splitter with the BEEF library using EME.

Geometry:
    - Two single-mode input waveguides → MMI body → two output waveguides
    - SOI platform: 220 nm Si on SiO2
    - Dimensions from a typical fabricated C-band 2x2 device (mm_length = 28 µm)

Method: Eigenmode Expansion (EME) via meow.  Fast, no FDTD required.
The top input is excited and each output is measured in turn.  Because both
inputs (and both outputs) sit only waveguide_offset = 1.5 µm apart, the EME
cross-sections are restricted to a single input branch (input_idx) and a
single output branch (output_idx) so the mode solver sees one isolated
waveguide at each end rather than a two-guide supermode.

A 2x2 MMI at its first self-imaging length is a balanced splitter: exciting
one input sends ~50% to the output on the same side ("bar") and ~50% to the
opposite side ("cross"), so each port carries an ideal −3 dB.

Run from any directory:
    python /path/to/examples/mmi_2x2/python/mmi_2x2_eme.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")           # headless — swap for "Qt5Agg" for interactive

from beef import MMI, EMESolver, SOI_220, plot_eme_geometry, plot_performance

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_2x2", "mmi_2x2_eme"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by the EME solver
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
)

print(f"MMI total length: {mmi.total_length:.2f} µm")
print(f"Input ports at y = {mmi.input_ys} µm; outputs at y = {mmi.output_ys} µm")

# ---------------------------------------------------------------------------
# 2. Configure EME solver
# ---------------------------------------------------------------------------
solver = EMESolver(
    # 5 regions: wg_in | taper_in | MMI body | taper_out | wg_out
    num_modes=[1,   7,   15,  7,   1  ],  # modes per region
    num_cells=[1,   15,  1,   15,  1  ],  # cells per region (tapers subdivided)
    lengths=  [0.5, 4.0, 28.0, 4.0, 0.5], # region lengths (µm)
    oxide_width=5.0,                       # cross-section mesh y-span (µm)
    oxide_height=1.5,                      # cross-section mesh z-span (µm)
    res=RESOLUTION,                        # mesh pts/µm
    material_stack=SOI_220,                # 220 nm SOI stack
    mode_solver="tidy3d",                  # tidy3d-backed mode solver (no FDTD section to stay consistent with)
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
os.makedirs(output_dir, exist_ok=True)
plot_eme_geometry(
    mmi, solver=solver,
    region_labels=["wg in", "taper in", "MMI body", "taper out", "wg out"],
    title="2×2 MMI Splitter — Top-down Geometry",
    output_path=os.path.join(output_dir, "geometry.png"),
    # output_idx omitted here — draws both output branches for the reference
    # plot (each solve() below restricts to one branch at a time).
)
print(f"Geometry plot saved → {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)  # µm, C-band

print(f"Sweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm")

# Excite the top input (input_idx=1) and measure each output branch in turn.
# input_idx/output_idx restrict the EME cross-sections to one input and one
# output waveguide so the mode solver sees an isolated single-mode guide at
# each end (the two branches are waveguide_offset - waveguide_width = 1.05 µm
# apart edge-to-edge past the taper, comfortably decoupled). output_idx=1 is
# the top output (bar), output_idx=0 the bottom output (cross).
smat_bar   = solver.solve(mmi, wavelengths, ports=["bar"],   input_idx=1, output_idx=1, save_dir=save_dir)
smat_cross = solver.solve(mmi, wavelengths, ports=["cross"], input_idx=1, output_idx=0, save_dir=save_dir)
smat = smat_bar
smat._data.update(smat_cross._data)
smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength

# ---------------------------------------------------------------------------
# 5. Print summary
# ---------------------------------------------------------------------------
T_bar, T_cross = smat.transmission("bar"), smat.transmission("cross")
T_total = T_bar + T_cross
wl_center = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_center))

print(f"\nResults at λ = {wl_center:.3f} µm:")
print(f"  Bar output   T = {T_bar[idx]:.4f}  IL = {10*np.log10(max(T_bar[idx],1e-20)):.2f} dB")
print(f"  Cross output T = {T_cross[idx]:.4f}  IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total        T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  (ideal 50/50 split → T_each=0.5, IL_total=0.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix.pkl")

# Plot insertion loss and transmission vs wavelength — plot_performance()
# defaults to every port in smat, so both outputs show up automatically.
plot_performance(smat, output_dir=output_dir, title="2×2 MMI Splitter (EME)")
print(f"Performance plot saved → {output_dir}/performance.png")
