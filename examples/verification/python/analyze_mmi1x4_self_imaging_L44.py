"""Self-imaging field snapshot at mm_length=44.0 um -- the corrected hybrid
body-length sweep's own global spread-minimum candidate (job 11782082,
after the body_stub fix -- see mmi_1x4_hybrid_body_length_sweep.py's own
module docstring). Companion to analyze_mmi1x4_self_imaging.py (the
mm_length=206 version); same method, different fixed length.

Purpose: read off the ACTUAL transverse lobe positions of the self-imaging
pattern at the body's own end (z = body_start + 44.0), to estimate what
waveguide_offset would align the output tapers with where the light
actually is at this length -- the current production offset (2.4 um) was
tuned against mm_length=206's own pattern, a physically different point on
the beat curve, so there is no reason it should also fit 44 um's pattern.
This script's own peak-detection feeds mmi_1x4_sin_taper_offset_sweep.py's
choice of where to center its (coarse, 5-point) sweep, rather than that
script scanning the whole geometrically-feasible offset range blind.

Field pattern here is genuinely output_idx-independent -- the body and input
sections don't depend on output_idx at all (only the final taper_out/wg_out
stub does, which sits PAST the z-cut this script reads), so OUTPUT_IDX below
is an arbitrary, inconsequential choice, not a per-branch result.

Run:
    python .../analyze_mmi1x4_self_imaging_L44.py
"""
import os

import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import meow as mw

from beef import MMI, EMESolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification", "mmi_1x4_sin_taper_offset_sweep"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Pinned to mmi_1x4/python/mmi_1x4_eme.py -- see that script for the full
# provenance of every constant here, EXCEPT MM_LENGTH (fixed at the
# corrected hybrid sweep's short-length candidate, not the production 206.0).
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 2.4, 12.0
STRAIGHT_LENGTH = 0.5
MM_LENGTH = 44.0
N_BODY, N_TAPER, C_TAPER = 20, 10, 8

WL_PLOT = 1.55
OUTPUT_IDX = 0   # arbitrary -- see module docstring, the body's own field doesn't depend on this

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30,
)

solver = EMESolver(
    num_modes=[1, N_TAPER, N_BODY, N_TAPER, 1],
    num_cells=[1, C_TAPER, 1, C_TAPER, 1],
    lengths=[STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, STRAIGHT_LENGTH],
    oxide_width=MM_WIDTH + 8.0, oxide_height=SIN_300.sz,
    res=RESOLUTION, material_stack=SIN_300, mode_solver="tidy3d",
)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

print(f"Building {sum(solver.num_cells)} cells at MM_LENGTH={MM_LENGTH}...")
cells = solver._build_cells(mmi, side=f"out{OUTPUT_IDX}", wavelengths=wavelengths,
                             output_idx=OUTPUT_IDX)
modes_per_cell = solver._modes_per_cell()

print(f"Solving modes at lambda={WL_PLOT} um for all {len(cells)} cells...")
env = mw.Environment(wl=float(WL_PLOT), T=solver.temperature)
css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]

z_wg_in_end = STRAIGHT_LENGTH
z_body_start = z_wg_in_end + TAPER_LENGTH
z_body_end = z_body_start + MM_LENGTH
z_taper_out_end = z_body_end + TAPER_LENGTH
print(f"\nRegion boundaries (um): body [{z_body_start}, {z_body_end}], "
      f"taper_out end {z_taper_out_end}")
print(f"True output port y-positions (current offset={WAVEGUIDE_OFFSET}): {mmi.output_ys} um")

z = np.linspace(cells[0].z_min, cells[-1].z_max, 800)
field, x = mw.propagate_modes(modes, cells, y=0.0, z=z, excite_mode_l=0, component="Ex")
intensity = np.abs(field) ** 2

idx_mouth = int(np.argmin(np.abs(z - z_body_end)))
cut = intensity[idx_mouth, :]

# Peak detection on the transverse cut at the body's own end -- the actual
# self-imaging lobe positions at THIS length, not assumed to match the
# current (206-um-tuned) port positions at all.
peak_idx, _ = find_peaks(cut, height=cut.max() * 0.15, distance=int(0.3 * RESOLUTION))
peak_ys = x[peak_idx]
peak_heights = cut[peak_idx]
print(f"\nDetected {len(peak_ys)} peaks at z={z[idx_mouth]:.2f} um: "
      f"y = {[f'{y:.3f}' for y in peak_ys]} um (heights {[f'{h:.1f}' for h in peak_heights]})")

if len(peak_ys) >= 2:
    spacings = np.diff(np.sort(peak_ys))
    print(f"Adjacent peak spacings: {[f'{s:.3f}' for s in spacings]} um")
    print(f"Mean adjacent spacing (candidate waveguide_offset): {spacings.mean():.3f} um")

fig, axes = plt.subplots(1, 2, figsize=(16, 5), gridspec_kw={"width_ratios": [2.4, 1]})

ax = axes[0]
im = ax.pcolormesh(z, x, np.abs(field).T, shading="auto", cmap="jet")
for zb, label in [(z_body_start, "body start"), (z_body_end, "body end / taper mouth"),
                   (z_taper_out_end, "taper end")]:
    ax.axvline(zb, color="white", lw=1.2, ls="--")
    ax.text(zb, x.max() * 0.92, label, rotation=90, color="white", fontsize=7, ha="right", va="top")
for y_port in mmi.output_ys:
    ax.axhline(y_port, color="magenta", lw=0.8, ls=":", alpha=0.8)
plt.colorbar(im, ax=ax, label="|Ex|")
ax.set_xlabel("z [µm]")
ax.set_ylabel("x [µm]")
ax.set_title(f"|Ex|(x, z) -- mm_length={MM_LENGTH} um, annotated\n"
             f"(dashed white = region seams, dotted magenta = CURRENT 206-um-tuned y-positions)")

ax2 = axes[1]
ax2.plot(cut, x, color="firebrick", lw=1.2)
for y_port in mmi.output_ys:
    ax2.axhline(y_port, color="0.6", lw=0.8, ls=":")
for py, ph in zip(peak_ys, peak_heights):
    ax2.plot(ph, py, "kx", ms=8)
    ax2.text(ph * 1.02, py, f"{py:+.2f}", fontsize=7, va="center")
ax2.set_xlabel("|Ex|^2 (a.u.)")
ax2.set_ylabel("x [µm]")
ax2.set_title(f"Transverse cut at z={z[idx_mouth]:.1f} um\n"
              f"(x = detected peaks; dotted gray = current 206-um port y's)")
ax2.set_ylim(x.min(), x.max())

plt.tight_layout()
out_path = os.path.join(output_dir, f"self_imaging_annotated_L{MM_LENGTH:.0f}_{WL_PLOT:.3f}um.png")
fig.savefig(out_path, dpi=200)
plt.close("all")
print(f"\nAnnotated plot saved -> {out_path}")
