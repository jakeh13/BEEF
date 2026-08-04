"""Tutorial: 1x4 MMI splitter on 300 nm Si3N4 -- FDTD simulation (BEEF method).

validation_markdown/mmi_1x4_sin300_validation.md step 6c.

Geometry: single-mode-ish 1.1 um input waveguide -> MMI body -> four 1.1 um
output waveguides. 300 nm Si3N4 core on SiO2 (SIN_300). Full device in a
single FDTD domain (mmi_1x2_fdtd.py's pattern), not the section-by-section
hybrid (that's mmi_1x4_hybrid.py).

Unlike mmi_1x2_fdtd.py's two half-domain monitor windows (each spanning from
the shared y=0 axis to the near PML edge -- a trick that only works for
TWO ports), this device has FOUR output ports at y = [-3.6, -1.2, +1.2,
+3.6] um (waveguide_offset=2.4 um spacing) -- "four monitors cannot each
take half the domain" (validation doc step 6c). Monitor widths are instead
sized directly off the port pitch: 2.0 um each (< 2.4 um spacing, leaving
0.2 um clearance per side to the neighboring monitor), centered on each
mmi.output_ys value.

Crosstalk between the four output waveguides was checked separately and
found negligible (T_cross=0.0002 at the actual 1.3 um edge-to-edge gap over
the taper+stub span -- validation doc step 6b/6c/6d, job 11667138), so this
monitor sizing is a numerical-isolation choice (don't let one monitor's
window catch a neighbor's core), not compensation for a real physical
coupling problem.

*** MM_LENGTH IS A PLACEHOLDER, NOT YET THE MEASURED VALUE -- see
mmi_1x4_eme.py's module docstring for the same note. ***

Run:
    python /path/to/examples/mmi_1x4/python/mmi_1x4_fdtd.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, FDTDSolver, SIN_300, plot_fdtd_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_1x4", "mmi_1x4_fdtd"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/um)
WL_MIN = 1.50       # C-band (step 4c's settled choice)
WL_MAX = 1.60
N_WL = 11

# ---------------------------------------------------------------------------
# 1. Define geometry -- steps 4-5's design.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
WAVEGUIDE_OFFSET = 2.4
MM_WIDTH = 12.0

MM_LENGTH = None   # TODO: fill in from mmi_1x4_sin_body_length_sweep_merge.py
                   # (job array 11667114) -- splice into that script's
                   # mm_lengths grid and assert the hit is exact (coupler-bug
                   # lesson) before trusting this value.
if MM_LENGTH is None:
    raise NotImplementedError(
        "MM_LENGTH is still a placeholder -- fill in the measured equal-split "
        "length from mmi_1x4_sin_body_length_sweep_merge.py (job array "
        "11667114) before running this example."
    )

mmi = MMI(
    n_inputs=1, n_outputs=4,
    mm_length=MM_LENGTH,
    mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH,
    taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30,
)
print(f"MMI total length: {mmi.total_length:.2f} um")
print(f"Output ports at y = {mmi.output_ys} um")

# ---------------------------------------------------------------------------
# 2. Configure FDTD solver
# ---------------------------------------------------------------------------
MONITOR_WIDTH = 2.0   # um -- < waveguide_offset=2.4, 0.2 um clearance per
                      # side to the neighboring monitor (see module docstring)

solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=MM_WIDTH + 8.0,   # generous lateral clearance past the body edges
    sim_height=SIN_300.sz,      # 5.0 um, step 4's converged value
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,               # full 3-D FDTD: real 300 nm Si3N4 core in oxide
    num_input_modes=1,
    num_output_modes=1,
    output_monitor_centers=list(mmi.output_ys),
    output_monitor_widths=[MONITOR_WIDTH] * 4,
    material_stack=SIN_300,
    decay_tol=1e-11,
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
plot_fdtd_geometry(
    mmi, solver=solver,
    title="1x4 MMI on 300 nm Si3N4 -- Top-down Geometry (FDTD)",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved -> {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep -- single broadband pulse, DFT monitors extract every
#    wavelength (and every output port) in one simulation.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nSweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um")
print("(single broadband pulse; DFT monitors extract all wavelengths in one sim)")

smat = solver.solve(mmi, wavelengths, port="top", input_mode=1, save_dir=save_dir)  # single input -- "top"/"bot" collapse to the same position
smat.check_passive()   # sanity: total output power must not exceed 1 at any wavelength

# ---------------------------------------------------------------------------
# 5. Print summary
# ---------------------------------------------------------------------------
wl_center = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_center))

# "out_0".."out_3" match the order of output_monitor_centers, i.e. mmi.output_ys
T = [smat.transmission(port=f"out_{k}", in_mode=1, out_mode=1) for k in range(4)]
T_total = sum(T)

print(f"\nResults at lambda = {wl_center:.3f} um:")
for k in range(4):
    print(f"  Output {k} (y={mmi.output_ys[k]:+.2f} um) power : {T[k][idx]:.4f}")
print(f"  Total (sum)              : {T_total[idx]:.4f}")
print(f"  Insertion loss           : {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  (ideal 1x4 split -> T_each=0.25, IL_total=0.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved -> {output_dir}/smatrix.pkl")

fig, ax = plt.subplots(figsize=(8, 4.5))
colors = ["steelblue", "darkorange", "seagreen", "firebrick"]
for k in range(4):
    ax.plot(wavelengths, T[k], label=f"Output {k} (y={mmi.output_ys[k]:+.2f} um)", color=colors[k])
ax.plot(wavelengths, T_total, label="Total", color="k", lw=1.5)
ax.axhline(0.25, color="gray", linestyle=":", lw=0.8, label="Ideal per-port (0.25)")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Transmission")
ax.set_title("1x4 MMI on 300 nm Si3N4 -- FDTD Transmission vs Wavelength")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "transmission.png"), dpi=300)
plt.close("all")
print(f"Transmission plot saved -> {output_dir}/transmission.png")
