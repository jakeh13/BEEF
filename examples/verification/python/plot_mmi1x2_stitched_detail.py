"""Detailed per-section panels for the 1x2 stitched hybrid field map.

Renders the FDTD input section and the EME body separately and at proper size,
plus transverse cuts, from the res-40 Step 3 run's saved capture. The combined
map compresses 5.25 um into one wide strip and hides structure; this is the
readable version.

Needs the capture from stitch_mmi1x2_hybrid_fields.py. Re-solves the body
modes (~55 s at res 40) because the driver does not persist them.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from beef import (MMI, EMESolver, SOI_160, solve_body_modes,
                  load_captured_section, project_interface_onto_body)

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
RES = float(os.environ.get("MMI1X2_STITCH_RES", 40.0))
base = os.path.join(_ROOT, "example_output", "mmi_1x2",
                    f"mmi_1x2_hybrid_stitched_res{int(RES)}")
cap = load_captured_section(os.path.join(base, "captured_input_section",
                                         "section_fields_band1.npz"))

WL, N_MODES, STUB = 1.31, 8, 0.25
mmi = MMI(n_inputs=1, n_outputs=2, mm_length=4.5, mm_width=2.0,
          waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
          waveguide_offset=1.0, waveguide_height=0.16)
hb, ht = mmi.mm_length / 2, mmi.total_length / 2
x_b0, x_b1 = -hb + STUB, hb - STUB

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
E_body = (np.exp(2j * np.pi * np.outer(xs, solved["neff"]) / WL)
          * proj["c"][None, :]) @ phi
xb, yb, Ib = x_b0 + xs, solved["y"], np.abs(E_body) ** 2

xi, yi = np.asarray(cap["x"]), np.asarray(cap["y"])
Ii = np.abs(np.asarray(cap[comp])) ** 2

fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.9], hspace=0.42, wspace=0.22)

# (a) input section, its own scale
ax = fig.add_subplot(gs[0, :])
im = ax.pcolormesh(xi, yi, Ii.T, shading="auto", cmap="inferno")
ax.axhline(0.175, color="cyan", lw=.7, ls="--", alpha=.7)
ax.axhline(-0.175, color="cyan", lw=.7, ls="--", alpha=.7)
ax.axhline(1.0, color="lime", lw=.8, ls="-", alpha=.5)
ax.axhline(-1.0, color="lime", lw=.8, ls="-", alpha=.5)
ax.set_ylim(-1.6, 1.6); ax.set_ylabel("device y (µm)")
ax.set_title(f"(a) INPUT section — FDTD, captured  |{comp}|²   "
             f"x ∈ [{xi.min():.2f}, {xi.max():.2f}] µm\n"
             "cyan = 0.35 µm access waveguide edges, green = 2.0 µm MMI body "
             "wall (taper opens 0.35 → 0.8 µm)", fontsize=9)
plt.colorbar(im, ax=ax, pad=.01)

# (b) body, its own scale
ax = fig.add_subplot(gs[1, :])
im = ax.pcolormesh(xb, yb, Ib.T, shading="auto", cmap="inferno")
for p in mmi.output_ys:
    ax.axhline(p, color="white", lw=1.1, ls=":", alpha=.95)
ax.axhline(1.0, color="lime", lw=.8, alpha=.5); ax.axhline(-1.0, color="lime", lw=.8, alpha=.5)
ax.set_ylim(-1.6, 1.6); ax.set_ylabel("device y (µm)")
ax.set_title(f"(b) BODY — EME, driven by the captured interface field  |{comp}|²   "
             f"x ∈ [{xb.min():.2f}, {xb.max():.2f}] µm\n"
             "white dotted = output port centres ±0.5 µm, green = body walls "
             "±1.0 µm", fontsize=9)
plt.colorbar(im, ax=ax, pad=.01)

# (c) transverse cuts through the chain
ax = fig.add_subplot(gs[2, 0])
for xq, lab, c in ((xi[0], f"input entry x={xi[0]:+.2f}", "steelblue"),
                   (xi[-1], f"seam (FDTD side) x={xi[-1]:+.2f}", "darkorange")):
    p = Ii[int(np.argmin(np.abs(xi - xq)))]
    ax.plot(yi, p / p.max(), color=c, label=lab)
p0 = Ib[0]; ax.plot(yb, p0 / p0.max(), color="crimson", ls="--",
                    label=f"seam (body side) x={xb[0]:+.2f}")
pe = Ib[-1]; ax.plot(yb, pe / pe.max(), color="seagreen",
                     label=f"body exit x={xb[-1]:+.2f}")
for p in mmi.output_ys:
    ax.axvline(p, color="gray", lw=.8, ls=":")
ax.set_xlim(-1.4, 1.4); ax.set_xlabel("device y (µm)")
ax.set_ylabel("|E|² (each normalised)")
ax.set_title("(c) transverse profiles along the chain\n"
             "grey dotted = output ports", fontsize=9)
ax.legend(fontsize=7); ax.grid(alpha=.3)

# (d) modal content: real vs the Gaussian Step 1 assumed
ax = fig.add_subplot(gs[2, 1])
w = proj["weights"] / proj["weights"].sum()
gauss = np.array([0.9022, 0, 0.0935, 0, 0.0013, 0, 0, 0])
gauss = gauss / gauss.sum()
k = np.arange(N_MODES); bw = 0.38
ax.bar(k - bw/2, gauss, bw, label="Gaussian launch (Step 1)", color="lightgray",
       edgecolor="gray")
ax.bar(k + bw/2, w, bw, label="REAL captured field", color="seagreen")
for i in (0, 2, 4):
    if gauss[i] > 0:
        ax.text(i + bw/2, w[i] + .02, f"{w[i]/gauss[i]:.1f}×", ha="center",
                fontsize=8, color="darkgreen")
ax.set_xlabel("body mode index"); ax.set_ylabel("normalised |c$_k$|²")
ax.set_title("(d) what the body is actually driven with\n"
             "labels = real/Gaussian ratio", fontsize=9)
ax.set_xticks(k); ax.legend(fontsize=8); ax.grid(alpha=.3, axis="y")

fig.suptitle(f"1×2 MMI hybrid — input + body fields at res {RES:.0f}, "
             f"λ = {WL} µm  (output section NOT included: 81% of the device)",
             fontsize=12)
out = os.path.join(base, "fields", f"detail_input_body_res{int(RES)}.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=160, bbox_inches="tight")
print(f"saved -> {out}")
print(f"body exit peak at |y| = "
      f"{abs(yb[int(np.argmax(np.where(yb>0, Ib[-1], 0)))]):.3f} um "
      f"(ports at {max(mmi.output_ys):.3f})")
