"""
Tutorial: standard directional coupler — single-domain full-device FDTD.

Geometry (220 nm SOI, C-band) — the same device as the hybrid tutorial
(directional_coupler_hybrid.py), at ONE hard-coded coupling length:

    "in"  ___                                      ____  "through"
             \\                                    /
              \\__________________________________/        (coupling gap)
              /                                    \\
    ___      /                                      \\___  "cross"
   "isolated"

    |<- bend_length ->|<--- coupling_length --->|<- bend_length ->|

Method: brute-force reference FDTD via meep. The WHOLE coupler — input
S-bends, straight coupled-pair region, output S-bends — sits in a single
simulation domain. One broadband GaussianSource pulse excites the
fundamental mode of the "in" arm; DFT mode monitors at the two output arms
extract the complex transmission at every wavelength in one run, so the
entire wavelength sweep costs exactly **1 FDTD simulation**.

Contrast with the hybrid tutorial, which decomposes the same device into
[FDTD bend | EME coupling region | FDTD bend]:

  - Here nothing is decomposed or approximated by a truncated mode basis —
    every scattering/radiation path in the real structure is captured in
    full. That makes this script the validation reference for the hybrid
    result at this coupling length.
  - The price is that the domain contains the entire coupling_length of
    essentially-uniform waveguide — the exact region EME propagates for
    free — and the run cost grows with coupling_length. Sweeping the
    coupling length this way means one full FDTD run *per length*, whereas
    the hybrid example sweeps 201 lengths for the cost of one (see its
    module docstring). Use FDTD to spot-check; use the hybrid to explore.

The coupling length is hard-coded to 15.15 µm — the measured 50:50 (3 dB)
splitting length for this cross-section, read off the hybrid example's
coupling_length_sweep.png (its own PLOT_COUPLING_LENGTH) and verified here.

Two details make the single-domain launch clean (both handled by the
library, shown here for the tutorial value):

  - The source is placed on the "in" arm via DirectionalCoupler.input_ys
    (port="bot" → index 0, the bottom arm at y = -waveguide_offset/2).
  - Both arms are physically present at the input cross-section, only
    waveguide_offset apart — a full-domain source aperture centered on the
    "in" arm would reach the unexcited "isolated" arm (corrupting the mode
    solve into the arm pair's even/odd supermodes) and the PML.
    FDTDSolver.solve() therefore shrinks the source + normalization monitor
    to the geometry's own launch window (DirectionalCoupler.
    source_aperture_at, clamped to the non-PML domain), exactly like the
    hybrid path does for its sections.

Run:
    python /path/to/examples/directional_coupler/python/directional_coupler_fdtd.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import DirectionalCoupler, FDTDSolver, SOI_220, plot_fdtd_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "directional_coupler", "directional_coupler_fdtd"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm)
WL_MIN = 1.50       # wavelength sweep start (µm)
WL_MAX = 1.60       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

# The 3 dB length for this cross-section — the β-calibrated hybrid's own
# auto-picked 50:50 split point (directional_coupler_hybrid.py,
# PLOT_COUPLING_LENGTH, job 11544656): 15.150754 µm, its nearest sampled
# length on a 200-point sweep grid. Hard-coded here to the same value so this
# script's plot title and the hybrid's agree to the decimal instead of two
# independent numbers that happen to be close. (The full-device FDTD ground
# truth — example_output/verification/ground_truth_sweep/length_sweep.npz,
# directional_coupler_validation.md item 5 — crosses 0.5 at 15.129 µm
# independently, 0.14% away; that was this constant's source from
# 2026-07-28 until this change moved it to the hybrid's own number instead.)
#
# WAS 13.5 µm until 2026-07-28, and that was wrong: the reference gives
# cross/(cross+through) = 0.4318 there, i.e. a 43:57 splitter, not 3 dB. It came
# from the pre-fix hybrid sweep, which reported L₅₀ = 13.65 µm *with the
# wavelength axis reversed* (item 1), and it survived both that fix and item 5's
# ground-truth rebuild because nothing rechecked it.
COUPLING_LENGTH = 15.15   # µm, hard-coded — the 3 dB length; see above

# ---------------------------------------------------------------------------
# 1. Geometry — identical to directional_coupler_hybrid.py except that
#    coupling_length is pinned instead of swept.
# ---------------------------------------------------------------------------
dc = DirectionalCoupler(
    waveguide_width=0.45,             # single-mode strip width (µm)
    gap=0.25,                         # coupling-region edge-to-edge gap (µm)
    coupling_length=COUPLING_LENGTH,  # straight parallel interaction length (µm)
    bend_length=9.0,                  # horizontal extent of each S-bend (µm)
    waveguide_offset=3.0,             # far-end centre-to-centre port pitch (µm)
    bend_type="euler",                # "euler" (cubic-Bezier S-bend) or "circular" (matched arcs)
    waveguide_height=0.22,            # Si layer thickness (µm) — SOI_220
)
print(f"Directional coupler total length: {dc.total_length:.2f} µm "
      f"(coupling length {dc.coupling_length:.2f} µm, "
      f"coupling pitch {dc.coupling_pitch:.3f} µm, far pitch {dc.waveguide_offset:.2f} µm)")

# ---------------------------------------------------------------------------
# 2. Configure the FDTD solver — one domain around the whole device.
#
#    The output arms sit at y = ±waveguide_offset/2 = ±1.5 µm; each output
#    monitor is a half-domain window (3 µm wide, from the shared y=0 axis to
#    the near PML edge) centered on its own arm, so each measures only its
#    arm's fundamental mode without reaching the other arm or the PML —
#    the same window the source/norm monitor get automatically from
#    DirectionalCoupler.source_aperture_at on the input side.
# ---------------------------------------------------------------------------
solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=6.0,        # >= waveguide_offset + margin, so both arms sit inside PML clearance
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,       # full 3-D FDTD: real 220 nm Si core in oxide, PML on all six faces
    num_input_modes=1,
    num_output_modes=1,   # fundamental mode per output arm
    # Order fixes the SMatrix port names: out_0 = bottom arm ("through",
    # same arm as "in"), out_1 = top arm ("cross", the coupled arm).
    output_monitor_centers=[-1.5, +1.5],
    output_monitor_widths=[3.0, 3.0],
    material_stack=SOI_220,
    decay_tol=1e-11,
)

# ---------------------------------------------------------------------------
# 3. Geometry plot — the full single simulation domain: source + norm
#    monitor on the "in" arm (shrunk to the geometry's launch window), one
#    output monitor per arm, PML on all four sides.
# ---------------------------------------------------------------------------
plot_fdtd_geometry(
    dc, solver=solver,
    port="bot",   # bottom input arm = "in" (see DirectionalCoupler.input_ys)
    title=f"Directional Coupler — Full-Device FDTD (coupling_length={COUPLING_LENGTH:.2f} µm)",
    output_path=os.path.join(output_dir, "geometry_fdtd.png"),
)
print(f"Geometry plot saved → {output_dir}/geometry_fdtd.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep — a single broadband run. port="bot" excites the
#    fundamental TE mode of the "in" arm; the DFT monitors resolve all
#    N_WL wavelengths from that one pulse.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nWavelength sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, single FDTD run)\n")

smat = solver.solve(dc, wavelengths, port="bot", input_mode=1, save_dir=save_dir)
smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength

# ---------------------------------------------------------------------------
# 5. Print summary at the centre wavelength
# ---------------------------------------------------------------------------
T_through = smat.transmission(port="out_0", in_mode=1, out_mode=1)  # bottom arm, same as "in"
T_cross = smat.transmission(port="out_1", in_mode=1, out_mode=1)    # top arm, coupled
T_total = T_through + T_cross

wl_c = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_c))

print(f"Results at coupling_length = {COUPLING_LENGTH:.2f} µm, λ = {wl_c:.3f} µm:")
print(f"  Through output T = {T_through[idx]:.4f}  IL = {10*np.log10(max(T_through[idx],1e-20)):.2f} dB")
print(f"  Cross output   T = {T_cross[idx]:.4f}  IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total          T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  Coupling ratio   = {T_cross[idx]/max(T_total[idx],1e-20):.3f}  "
      f"(0 = fully through, 1 = fully cross; ~0.5 expected at the 3 dB length)")

# ---------------------------------------------------------------------------
# 6. Save and plot — same two-panel layout as the hybrid example's
#    hybrid_results.png, for side-by-side comparison.
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix_fdtd.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix_fdtd.pkl")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_through, label="Through", color="steelblue")
ax.plot(wavelengths, T_cross,   label="Cross",   color="darkorange", linestyle="--")
ax.plot(wavelengths, T_total,   label="Total",   color="k", lw=1.5)
ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title(f"Transmission vs Wavelength (coupling_length={COUPLING_LENGTH:.2f} µm)")
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
ax.set_title(f"Insertion Loss vs Wavelength (coupling_length={COUPLING_LENGTH:.2f} µm)")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle("Directional Coupler — Single-Domain FDTD Simulation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "fdtd_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/fdtd_results.png")
