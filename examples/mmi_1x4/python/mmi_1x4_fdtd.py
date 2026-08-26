"""Tutorial: 1x4 MMI splitter on 300 nm Si3N4 -- whole-device FDTD (BEEF
method). validation_markdown/mmi_1x4_sin300_validation.md step 8 check 4.

RE-OPENED 2026-08-19. This script was originally drafted (job history:
example_output/mmi_1x4/mmi_1x4_fdtd/geometry.png, a structural-validation
leftover) and then DROPPED, step 6c, 2026-08-04: at the then-production
MM_LENGTH=206.0 um, device total length was ~216 um, and "3-D time-domain
FDTD scales with volume... a different cost class than the ~5 um devices
mmi_1x2_fdtd.py/mmi_2x2_fdtd.py simulate whole." Step 9c's design-space
sweep (2026-08-14/15) since moved production to MM_LENGTH=39.90/
WAVEGUIDE_OFFSET=3.1 -- device total length is now 49.90 um (MEASURED,
mmi_1x4_eme.py's own real run), a 4.3x reduction -- worth reconsidering,
not previously possible at the old geometry.

Structural advantage over the hybrid/EME per-branch pattern: this device
has exactly ONE input port, so a whole-device FDTD run needs only ONE
excitation (not one per output branch) -- solver.solve() with 4 output
monitors reading simultaneously off the same broadband pulse, matching
mmi_1x2_fdtd.py's own 1-input/N-output pattern exactly, just with 4 output
monitors instead of 2.

Domain size, MEASURED via direct calculation (not assumed): whole-device
computational cell is 50.90 x 21.00 x 6.00 um -- 173.2M grid cells at
res=30, vs. the input-section timing probe's own 6.75 x 21 x 6 um
(21.7M cells) -- a 7.54x larger domain (linear in x only; y,z unchanged).
Per-step cost should scale roughly with cell count (~8x slower per step
than the input probe's own measured 0.995 s/step). The genuinely UNKNOWN
factor is how many total timesteps this longer domain needs to reach
decay_tol=1e-11 -- if similar to the input probe's own 23,388 steps, total
cost projects to ~52 hrs; if steps scale with the ~7.5x longer propagation
path, ~390 hrs. That is too wide a spread to commit to blind -- run
investigate_mmi1x4_res30_fdtd_wholedevice_timing.py (loose decay_tol probe
at the REAL x-extent) first, same two-step discipline as every other
expensive job this session, before trusting any number below.

num_input_modes=num_output_modes=1 -- same single-mode-launch convention
every other script for this device uses (WAVEGUIDE_WIDTH=1.1 um is a
genuinely 2-mode design point, per step 4's own finding; higher-order
content is treated as unwanted loss, a stated simplification, not assumed
away).

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
# Global simulation parameters -- pinned to mmi_1x4_eme.py's/
# mmi_1x4_hybrid.py's own production geometry (step 9c) and RESOLUTION
# choice (2026-08-18 measurement, see those scripts' own RESOLUTION
# comments -- res=30 gives real speedup on straight-waveguide-class domains,
# genuinely UNMEASURED for a whole-device run with bends/tapers present).
# ---------------------------------------------------------------------------
RESOLUTION = 30.0
WL_MIN = 1.50
WL_MAX = 1.60
N_WL = 11

# ---------------------------------------------------------------------------
# 1. Geometry -- identical to mmi_1x4_eme.py/mmi_1x4_hybrid.py, step 9c.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
WAVEGUIDE_OFFSET = 3.1
MM_WIDTH = 12.0
MM_LENGTH = 39.90

mmi = MMI(
    n_inputs=1, n_outputs=4,
    mm_length=MM_LENGTH,
    mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH,
    taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30,
    use_neighbor_aperture=True,
)
print(f"MMI total length: {mmi.total_length:.2f} um")
print(f"Output ports at y = {mmi.output_ys} um")

# ---------------------------------------------------------------------------
# 2. Configure FDTD solver -- one excitation (single input port), 4 output
#    monitors reading simultaneously off the same broadband pulse. Monitor
#    centers match mmi.output_ys exactly (ascending y order) so "out_k" here
#    is the SAME physical port as output_idx=k in mmi_1x4_eme.py/
#    mmi_1x4_hybrid.py, enabling a direct per-port comparison.
# ---------------------------------------------------------------------------
# Monitor width: 2.8 um per port -- less than the 3.1 um port pitch
# (WAVEGUIDE_OFFSET), leaving ~0.15 um clearance to each neighbor on both
# sides, comfortably containing the 1.1 um guided mode without overlapping
# the adjacent port's own monitor.
_MONITOR_WIDTH = 2.8

# FDTD_SIM_WIDTH/HEIGHT, CHANGED 2026-08-19 from MM_WIDTH+8.0/SIN_300.sz --
# see mmi_1x4_hybrid.py's own _fdtd_common comment for the measurement
# (mmi_1x4_sin_body_domain_size_convergence.py, job 12029224). Also checked
# here specifically for THIS script's own outermost-port PML clearance
# (a whole-device-only concern, since output_ys reach the domain edge
# directly rather than through a per-branch isolated section): at
# sim_width=16.0, PML inner edge sits at +-7.5 um; the outermost port
# (y=+-4.65) plus its own monitor half-width (1.4) reaches 6.05 um, leaving
# 1.45 um of clearance to PML -- comfortable for a 1.1 um guided mode.
FDTD_SIM_WIDTH = MM_WIDTH + 4.0
FDTD_SIM_HEIGHT = 4.0

solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=FDTD_SIM_WIDTH,
    sim_height=FDTD_SIM_HEIGHT,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    num_input_modes=1,
    num_output_modes=1,
    output_monitor_centers=list(mmi.output_ys),
    output_monitor_widths=[_MONITOR_WIDTH] * 4,
    material_stack=SIN_300,
    decay_tol=1e-11,
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
plot_fdtd_geometry(
    mmi, solver=solver,
    title="1x4 MMI on 300 nm Si3N4 -- Whole-Device FDTD",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved -> {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep -- single broadband pulse, one excitation, all 4
#    outputs read simultaneously.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nSweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um")
print("(single broadband pulse; DFT monitors extract all wavelengths in one sim)")

smat = solver.solve(mmi, wavelengths, port="top", input_mode=1, save_dir=save_dir)
smat.check_passive()   # sanity: total output power must not exceed 1 at any wavelength

# ---------------------------------------------------------------------------
# 5. Print summary
# ---------------------------------------------------------------------------
T = [smat.transmission(port=f"out_{k}", in_mode=1, out_mode=1) for k in range(4)]
T_total = sum(T)

wl_center = wavelengths[len(wavelengths) // 2]
idx = int(np.argmin(np.abs(wavelengths - wl_center)))

print(f"\nResults at lambda = {wl_center:.3f} um:")
for k in range(4):
    print(f"  Output {k} (y={mmi.output_ys[k]:+.2f} um)  "
          f"T = {T[k][idx]:.4f}  IL = {10*np.log10(max(T[k][idx], 1e-20)):.2f} dB")
print(f"  Total  T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx], 1e-20)):.2f} dB")
print(f"  (ideal 1x4 split: T_each=0.25, IL_total=0.0 dB)")
# Comparison references, CORRECTED 2026-08-26. This previously cited
# "EME (job 11840722, confirmed N_WL=11): T=[0.2414, 0.2396, 0.2396, 0.2414],
# total=0.9620" -- but 11840722 is `beef_mmi1x4_toc`, a step-9c DESIGN-SWEEP
# array that finished 2026-08-14, the same day the geometry changed. It is not
# a production EME run and its numbers predate the settled geometry. The
# authoritative EME production result is job 11887395 (and the res 40->60
# convergence grid, jobs 12323220/12331012, confirms it is converged).
print(f"\nComparison at band centre:")
print(f"  EME production   (job 11887395, res=40, 20x5, N_WL=11): "
      f"T=[0.2431, 0.2423, 0.2423, 0.2431], total=0.9707  (favours OUTER)")
print(f"  whole-device FDTD (job 12029541, res=30, 16x4, N_WL=1):  "
      f"T=[0.2400, 0.2424, 0.2424, 0.2400], total=0.9648  (favours INNER)")
print(f"  ^ THIS run is the N_WL=11 counterpart of 12029541 -- same method and")
print(f"    parameters, so its own band-centre row above should reproduce that")
print(f"    line. A mismatch there is a red flag about this run, not a physics")
print(f"    result. The EME-vs-FDTD gap (~3e-03/port) is check 1's real, and")
print(f"    resolution-converged, engine disagreement.")

# ---------------------------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved -> {output_dir}/smatrix.pkl")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
colors = ["steelblue", "darkorange", "seagreen", "firebrick"]

ax = axes[0]
for k in range(4):
    ax.plot(wavelengths, T[k], label=f"Output {k} (y={mmi.output_ys[k]:+.2f} um)", color=colors[k])
ax.plot(wavelengths, T_total, label="Total", color="k", lw=1.5)
ax.axhline(0.25, color="gray", linestyle=":", lw=0.8, label="Ideal per-port")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Transmission")
ax.set_title("Whole-device FDTD -- Transmission vs Wavelength")
ax.legend(fontsize=7)
ax.grid(True)
ax.set_ylim(0, 1.05)

ax = axes[1]
for k in range(4):
    IL = 10 * np.log10(np.clip(T[k], 1e-20, None))
    ax.plot(wavelengths, IL, label=f"Output {k}", color=colors[k])
IL_tot = 10 * np.log10(np.clip(T_total, 1e-20, None))
ax.plot(wavelengths, IL_tot, label="Total", color="k", lw=1.5)
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title("Whole-device FDTD -- Insertion Loss vs Wavelength")
ax.legend(fontsize=7)
ax.grid(True)

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- Whole-Device FDTD (check 4)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "transmission.png"), dpi=300)
plt.close("all")
print(f"Transmission plot saved -> {output_dir}/transmission.png")
