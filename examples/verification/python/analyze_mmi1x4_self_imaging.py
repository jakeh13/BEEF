"""Annotated self-imaging diagnostic for the 1x4 MMI EME field
(plot_mmi1x4_eme_fields.py's companion) -- answers "does the field plot
suggest a more optimal output-taper position and/or a shorter mm_length?"

WHY THIS EXISTS
---------------
plot_eme_fields() (beef/plot/fields.py) draws |field|(x, z) with plain axis
numbers -- no markers for where the body/taper boundaries or the four real
output ports actually sit. A bright interior self-imaging spot is easy to
misread as "the body could be cut shorter here", but for a 1x4 splitter the
four ports sit at y = [-3.6, -1.2, +1.2, +3.6] um -- NOT at y=0 -- so a
single bright lobe AT y=0 partway through the body is almost certainly the
classical "restored single image" every MMI shows at a fraction of its beat
length (a generic MMI phenomenon, independent of the N-way split intended),
not evidence of a better cut point for the actual 4-way split.

The rigorous answer to "would a different mm_length be better" already
exists on record: mmi_1x4_sin_body_length_sweep.py's 281-point sweep (step
5), which measured T_total/spread across ALL FOUR outputs simultaneously, at
every 1 um of mm_length from 20 to 300. That is a direct multi-port power
measurement, not a visual read of one field snapshot -- this script exists
to check the field plot AGAINST that data, with real markers, rather than
letting the field plot's texture stand in for it.

WHAT THIS ADDS on top of plot_mmi1x4_eme_fields.py
---------------------------------------------------
1. The same full (x, z) |Ex| map, but with vertical lines at the real region
   boundaries (taper_in/body seam, body/taper_out seam, taper_out/wg_out
   seam) and horizontal dashed lines at the four TRUE output y-positions --
   makes misreading an interior self-image spot as "near a port" or "near
   the taper mouth" much harder.
2. A transverse |Ex(y)|^2 cross-cut taken AT the body/taper_out seam (the
   taper mouth) -- the CORRECT place to check whether the self-imaging
   pattern lines up with the real output positions, not an arbitrary
   interior z. Ports marked with vertical lines.

Also confirms directly (not assumed) that output_idx restricts only the
FINAL straight stub, not the tapers: beef/geometry/mmi.py's own
MMI.polygons() docstring says "Output tapers are always drawn for every
branch"; output_indices only gates the stub_out_idxs check. So the
taper_out region of ANY output_idx's solve shows all four tapers at once,
and only the wg_out region isolates the one selected port.

Reuses plot_mmi1x4_eme_fields.py's exact cells/modes construction -- no new
solve, if that script's own field-plot job has already run and this is
being re-run standalone it re-solves (same "minutes" cost, not re-derived
here).

Run:
    python .../analyze_mmi1x4_self_imaging.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import meow as mw

from beef import MMI, EMESolver, SIN_300

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_1x4", "mmi_1x4_eme"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Pinned to mmi_1x4/python/mmi_1x4_eme.py -- identical to
# plot_mmi1x4_eme_fields.py's own config.
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 2.4, 12.0
STRAIGHT_LENGTH = 0.5
MM_LENGTH = 206.0
N_BODY, N_TAPER, C_TAPER = 20, 10, 8

WL_PLOT = 1.55
OUTPUT_IDX = 0

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

print(f"Building {sum(solver.num_cells)} cells for output branch {OUTPUT_IDX}...")
cells = solver._build_cells(mmi, side=f"out{OUTPUT_IDX}", wavelengths=wavelengths,
                             output_idx=OUTPUT_IDX)
modes_per_cell = solver._modes_per_cell()

print(f"Solving modes at lambda={WL_PLOT} um for all {len(cells)} cells...")
env = mw.Environment(wl=float(WL_PLOT), T=solver.temperature)
css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]

# ---------------------------------------------------------------------------
# Real region boundaries, computed the same way _build_cells() accumulates
# them (STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, STRAIGHT_LENGTH),
# not eyeballed off the plot.
# ---------------------------------------------------------------------------
z_wg_in_end = STRAIGHT_LENGTH
z_body_start = z_wg_in_end + TAPER_LENGTH
z_body_end = z_body_start + MM_LENGTH
z_taper_out_end = z_body_end + TAPER_LENGTH
z_device_end = z_taper_out_end + STRAIGHT_LENGTH
print(f"\nRegion boundaries (um): wg_in [0, {z_wg_in_end}], "
      f"taper_in [{z_wg_in_end}, {z_body_start}], "
      f"body [{z_body_start}, {z_body_end}], "
      f"taper_out [{z_body_end}, {z_taper_out_end}], "
      f"wg_out [{z_taper_out_end}, {z_device_end}]")
print(f"True output port y-positions: {mmi.output_ys} um "
      f"(output_idx={OUTPUT_IDX} -> y={mmi.output_ys[OUTPUT_IDX]:+.2f} um)")

z = np.linspace(cells[0].z_min, cells[-1].z_max, 800)
field, x = mw.propagate_modes(modes, cells, y=0.0, z=z, excite_mode_l=0, component="Ex")
intensity = np.abs(field) ** 2

fig, axes = plt.subplots(1, 2, figsize=(16, 5), gridspec_kw={"width_ratios": [2.4, 1]})

ax = axes[0]
im = ax.pcolormesh(z, x, np.abs(field).T, shading="auto", cmap="jet")
for zb, label in [(z_body_start, "body start"), (z_body_end, "body end / taper mouth"),
                   (z_taper_out_end, "taper end")]:
    ax.axvline(zb, color="white", lw=1.2, ls="--")
    ax.text(zb, x.max() * 0.92, label, rotation=90, color="white", fontsize=7,
            ha="right", va="top")
for y_port in mmi.output_ys:
    ax.axhline(y_port, color="magenta", lw=0.8, ls=":", alpha=0.8)
plt.colorbar(im, ax=ax, label="|Ex|")
ax.set_xlabel("z [µm]")
ax.set_ylabel("x [µm]")
ax.set_title(f"|Ex|(x, z) -- output branch {OUTPUT_IDX}, annotated\n"
             f"(dashed white = region seams, dotted magenta = TRUE output y-positions)")

# The correct place to judge "is the self-image where the ports are": a
# transverse cut AT the taper mouth (body end), not an arbitrary interior z.
idx_mouth = int(np.argmin(np.abs(z - z_body_end)))
ax2 = axes[1]
ax2.plot(intensity[idx_mouth, :], x, color="firebrick", lw=1.2)
for y_port in mmi.output_ys:
    ax2.axhline(y_port, color="0.4", lw=0.8, ls=":")
    ax2.text(intensity[idx_mouth, :].max() * 1.02, y_port, f"{y_port:+.2f}",
              fontsize=7, va="center")
ax2.set_xlabel("|Ex|^2 (a.u.)")
ax2.set_ylabel("x [µm]")
ax2.set_title(f"Transverse cut at z={z[idx_mouth]:.1f} um\n(taper mouth -- true port y's marked)")
ax2.set_ylim(x.min(), x.max())

plt.tight_layout()
out_path = os.path.join(output_dir, "fields", f"self_imaging_annotated_out{OUTPUT_IDX}_{WL_PLOT:.3f}um.png")
os.makedirs(os.path.dirname(out_path), exist_ok=True)
fig.savefig(out_path, dpi=200)
plt.close("all")
print(f"\nAnnotated plot saved -> {out_path}")
