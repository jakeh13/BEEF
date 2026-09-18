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

The coupling length is hard-coded to 15.1 µm — a round tenth of a micron from
the 50:50 (3 dB) splitting length this cross-section shows under the dispersive
material model, and the same number directional_coupler_hybrid.py pins its
PLOT_COUPLING_LENGTH to, so the two engines are compared at one identical
length. See COUPLING_LENGTH below for the measurements behind it (and for why
the value is a rounded convenience, not itself a measured crossing).

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
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm)
WL_MIN = 1.50       # wavelength sweep start (µm)
WL_MAX = 1.60       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

# The 3 dB length for this cross-section, taken from THIS engine's own
# reference: the full-device FDTD ground-truth length sweep
# (verify_item_5_ground_truth.py, array 11480725 / job 11482062, frozen
# material model, 40 lengths × 11 λ). Its band-mean cross/(cross+through)
# crosses 0.5 at 15.129 µm, bracketed by L=15.0 (0.494617) and L=15.5
# (0.515476) — hence 15.13. That is a *band-mean*; at λ=1.55 alone the same
# sweep crosses at 15.313 µm.
#
# Deliberately NOT pinned to the hybrid's 15.15. The hybrid's auto-picked
# 50:50 point (directional_coupler_hybrid.py, PLOT_COUPLING_LENGTH, job
# 11544656) is 15.150754 µm — its nearest sample on a 200-point, 0.299 µm
# sweep grid. These are two engines' independent crossings of the same
# physical quantity, 0.14% apart, and that difference is a *result* this
# device exists to measure, not a drift to be reconciled. Collapsing them to
# one number hides it. Commit e85b350 (a ring-resonator CACHE cleanup) did
# exactly that with no mention; reverted 2026-07-31. The prose here claimed
# the pinning for two days after the value went back — corrected 2026-08-02.
#
# NB the live sweep file is gone: clear.sh cleared example_output/verification/,
# so the npz survives only at archive/example_output_20260728/verification/
# ground_truth_sweep/length_sweep.npz. Regenerate via
# batch/verification/run_verify_item_5_ground_truth.batch.
#
# WAS 13.5 µm until 2026-07-28, and that was wrong: the reference gives
# cross/(cross+through) = 0.4318 there, i.e. a 43:57 splitter, not 3 dB. It came
# from the pre-fix hybrid sweep, which reported L₅₀ = 13.65 µm *with the
# wavelength axis reversed* (item 1), and it survived both that fix and item 5's
# ground-truth rebuild because nothing rechecked it.
COUPLING_LENGTH = 15.1    # µm, hard-coded — the tenth-of-a-micron round number
                          # nearest the DISPERSIVE 50:50 point, set 2026-08-03.
                          # Everything above this line is the FROZEN-era record
                          # and is why the value used to read 15.13; dispersion
                          # defaulted 2026-08-01 and moved L₅₀ to 15.089958 µm
                          # (hybrid, job 11609658) / ~15.074 µm inferred for this
                          # engine (job 11609659), and 15.1 is the nearest tenth
                          # to both. The two engines' crossings are still ~0.1%
                          # apart and that gap is a *result* this device measures
                          # — pinning both examples to one length is for the
                          # side-by-side comparison, NOT a claim they agree.
                          # As of 2026-08-03 the hybrid simulates a true 15.1
                          # too: it used to argmin its pinned length onto a
                          # 0.299 µm sweep grid (which contains neither 15.13
                          # nor 15.1) and silently land on 15.150754; it now
                          # splices the pinned length into the grid, so the two
                          # examples are at the identical length. Documented in
                          # archive/dc_validation_20260727.md item 14. Commit
                          # e85b350 (a ring-resonator CACHE cleanup) moved
                          # this to 15.15 with no mention of it; restored
                          # 2026-07-31.

# DC_COUPLING_LENGTH overrides it without editing this file. Added 2026-08-01
# so a hybrid run can hand its own measured 50:50 length straight to this
# script: directional_coupler_hybrid.py writes coupling_length_50.txt, and
# run_dc_fdtd_at_hybrid_l50.batch reads it into this variable.
#
# Why that chaining exists rather than a second hardcoded constant: every
# number in the comment above (15.13, 15.129, 15.150754) was measured on the
# FROZEN material model. Material dispersion defaulted on 2026-08-01, and L50
# is exactly the kind of quantity it can move — the split ratio integrates
# Δβ over the band. Re-using a frozen constant for a dispersive run would
# assume the answer. Leave the constant alone as the documented frozen-era
# reference; override it per-run.
#
# It moved, and the override earned itself. Measured 2026-08-02, jobs
# 11609658 (hybrid, DC_PLOT_COUPLING_LENGTH=auto) → 11609659 (this script at
# that length): the hybrid's band-mean L50 goes 15.150754 → **15.089958** µm
# under dispersion. Running this script there gives band-mean ratio
# **0.500675** across 1.50–1.60 µm — 0.13% off a true 50:50, i.e. the
# dispersive FDTD crossing is ~15.074 µm (0.000675 / 0.0418 per µm, the
# hybrid's local slope). So dispersion moves BOTH engines down by ~0.05 µm
# and they still agree to ~0.1%: 15.129/15.151 frozen, 15.074/15.090
# dispersive. The engine gap survives the material model — it is not a
# dispersion artifact.
#
# The file default moved 15.13 -> 15.1 on 2026-08-03 (user's call), to sit a
# tenth of a micron from the dispersive 50:50 point instead of on the frozen-era
# one. KEEP THE CAVEAT IN VIEW: ~15.074 is not a measured crossing for THIS
# engine — it rests on one length plus the *hybrid's* slope, and there is still
# no dispersive full-device FDTD length sweep. 15.1 is defensible because it is
# within a tenth of this engine's own numbers too (15.129 frozen, ~15.074
# inferred), so it is not the e85b350 mistake of pinning this engine to the
# other's answer; but it is a rounded convenience value, not a measurement.
# Regenerate the item-5 sweep dispersively before treating 15.1 as this engine's
# measured 3 dB length.
_env_CL = os.environ.get("DC_COUPLING_LENGTH")
if _env_CL is not None:
    COUPLING_LENGTH = float(_env_CL)
    print(f"DC_COUPLING_LENGTH -> COUPLING_LENGTH = {COUPLING_LENGTH:.6f} µm "
          f"(file default 15.1 µm)")

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
