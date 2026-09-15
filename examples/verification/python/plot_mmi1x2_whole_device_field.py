"""1x2 MMI — the WHOLE-DEVICE hybrid field map, all three sections.

Completes validation_markdown/hybrid_field_extraction_validation.md Step 3.
Assembles, on one device-frame axis and one absolute amplitude scale:

  FDTD input  [-3.25, -2.00]  captured (job 12736849)
  EME  body   [-2.00, +2.00]  driven by the captured interface field
  FDTD output [+2.00, +3.25]  8 FORWARD excites superposed (job 13189964)

⚠ WHAT IS AND IS NOT VALIDATED, because the distinction matters:

  VALIDATED
    * the input/body seam -- the 8-mode projection reproduces the FDTD
      interface profile (0.9818 of its power; the two curves overlie exactly)
    * modal content at that seam vs the transfer matrix, two independent
      routes agreeing to 2.1% max (job 12736849's Step 4 cross-check)
    * parity end to end -- the output branches come out 1.000 : 1.000, against
      1899 : 1 for the reverse-excited capture (job 13189964 CHECK 1)

  NOT VALIDATED
    * the PER-MODE normalisation of the output superposition. c_k come from an
      area-integral projection onto EME modes; E^(k) is meep's response to a
      unit EigenModeSource. The 50/50 branch check does NOT test this: the
      device is symmetric and only EVEN body modes are excited, so ANY
      even-mode superposition is symmetric. A wrong relative weighting between
      modes 0/2/4 would still give ratio 1.000.
    * therefore the output panel's INTERNAL STRUCTURE (how the field
      redistributes between x=+2.00 and the ports) is indicative, not
      quantitative. The check that would settle it is the per-branch
      transmission against the hybrid's own validated 0.3749 at 1.31 um,
      which needs a power normalisation not yet carried through.

Run:
    python .../plot_mmi1x2_whole_device_field.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import (MMI, EMESolver, SOI_160, solve_body_modes,
                  load_captured_section, project_interface_onto_body)

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
RES, WL, N_MODES, STUB = 40.0, 1.31, 8, 0.25
base = os.path.join(_ROOT, "example_output", "mmi_1x2",
                    f"mmi_1x2_hybrid_stitched_res{int(RES)}")

mmi = MMI(n_inputs=1, n_outputs=2, mm_length=4.5, mm_width=2.0,
          waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
          waveguide_offset=1.0, waveguide_height=0.16)
hb, ht = mmi.mm_length / 2, mmi.total_length / 2
x_b0, x_b1 = -hb + STUB, hb - STUB

# --- input (captured) ------------------------------------------------------
cap = load_captured_section(os.path.join(base, "captured_input_section",
                                         "section_fields_band1.npz"))
# --- body (projected + propagated) ----------------------------------------
eme = EMESolver(num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
                oxide_width=5.0, oxide_height=1.5, res=RES,
                material_stack=SOI_160, mode_solver="meep")
solved = solve_body_modes(mmi, eme, x_b0, x_b1, WL, N_MODES,
                          wavelengths=np.linspace(1.26, 1.36, 11))
proj = project_interface_onto_body(cap, solved, verbose=False)
comp = proj["component"]
F = solved["fields"][comp]
iz = int(np.argmax(np.abs(F).sum(axis=(0, 1))))
phi = (F / np.sqrt((np.abs(F) ** 2).sum(axis=(1, 2))
                   * solved["dy"] * solved["dz"])[:, None, None])[:, :, iz]
xs = np.linspace(0.0, x_b1 - x_b0, 700)
E_b = (np.exp(2j * np.pi * np.outer(xs, solved["neff"]) / WL)
       * proj["c"][None, :]) @ phi
xb, yb, Ib = x_b0 + xs, solved["y"], np.abs(E_b) ** 2
# --- output (forward superposition) ---------------------------------------
sup = np.load(os.path.join(base, "captured_output_forward_idx1",
                           "forward_superposition.npz"))
xo, yo, Io = sup["x"], sup["y"], sup["I"]

xi, yi = np.asarray(cap["x"]), np.asarray(cap["y"])
Ii = np.abs(np.asarray(cap[comp])) ** 2

# Shared scale across all three, since that is the whole point of projecting
# at the seam rather than normalising each panel.
vmax = max(Ii.max(), Ib.max(), Io.max())
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(15, 8),
                              gridspec_kw={"height_ratios": [2.1, 1]})
for X, Y, I in ((xi, yi, Ii), (xb, yb, Ib), (xo, yo, Io)):
    im = ax.pcolormesh(X, Y, I.T, shading="auto", cmap="inferno",
                       vmin=0, vmax=vmax)
for xs_ in (x_b0, x_b1):
    ax.axvline(xs_, color="deepskyblue", lw=1.0, ls=":", alpha=.9)
for p in mmi.output_ys:
    ax.axhline(p, color="white", lw=1.0, ls=":", alpha=.85)
ax.axhline(1.0, color="lime", lw=.8, alpha=.45); ax.axhline(-1.0, color="lime", lw=.8, alpha=.45)
for xc, lab in ((-2.62, "FDTD input\n(captured)"), (0.0, "EME body\n(real-field driven)"),
                (2.62, "FDTD output\n(8 forward excites)")):
    ax.text(xc, 1.32, lab, color="w", ha="center", fontsize=8.5)
ax.set_xlim(-ht, ht); ax.set_ylim(-1.6, 1.6)
ax.set_xlabel("device x — propagation (µm)"); ax.set_ylabel("device y (µm)")
ax.set_title(f"1×2 MMI — WHOLE-DEVICE hybrid field map, |{comp}|² at λ = {WL} µm  "
             f"(all {mmi.total_length:.2f} µm)\n"
             "one absolute amplitude scale across all three sections; "
             "cyan dotted = section seams, white dotted = output ports",
             fontsize=10)
plt.colorbar(im, ax=ax, pad=.01, label="|E|² (arb.)")

for X, I, lab, c in ((xi, Ii, "input", "steelblue"), (xb, Ib, "body", "crimson"),
                     (xo, Io, "output", "seagreen")):
    ax2.plot(X, I.sum(axis=1) / I.sum(axis=1).max(), color=c, label=lab, lw=1.4)
for xs_ in (x_b0, x_b1):
    ax2.axvline(xs_, color="deepskyblue", lw=1.0, ls=":")
ax2.set_xlim(-ht, ht); ax2.set_xlabel("device x (µm)")
ax2.set_ylabel("∫|E|²dy (per-section norm.)")
ax2.set_title("longitudinal power profile — each section normalised to its own max "
              "(the seam step is basis truncation, see docstring)", fontsize=9)
ax2.legend(fontsize=8); ax2.grid(alpha=.3)
fig.tight_layout()
out = os.path.join(base, "fields", "whole_device_res40_1.310um.png")
fig.savefig(out, dpi=160, bbox_inches="tight")
print(f"saved -> {out}")
print(f"coverage: x {-ht:+.2f} to {+ht:+.2f} um = 100% of the device")
pe = Io[-1]
print(f"output port-end lobes at |y| = "
      f"{abs(yo[int(np.argmax(np.where(yo>0, pe, 0)))]):.3f} um "
      f"(ports at {max(mmi.output_ys):.3f})")
