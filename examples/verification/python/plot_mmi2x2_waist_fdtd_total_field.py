"""Waisted 2x2 -- FDTD total-device field map, for the same device-frame
comparison stitch_mmi2x2_waist_hybrid_fields.py builds on the hybrid side.

Reads the whole-device capture already produced by capture_mmi2x2_waist_fields.py
MMI2X2_CAP_WHAT=full (job already run -- see the validation markdown's
captured_full_device_res40_wl1.5475 entry). No FDTD re-run.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI
from beef.plot.hybrid_stitch import (dominant_component, draw_geometry_outline,
                                     load_captured_section)

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
tapered_dir = os.path.join(_ROOT, "example_output", "mmi_2x2_tapered")
cap_path = os.path.join(tapered_dir, "captured_full_device_res40_wl1.5475",
                        "section_fields_band1.npz")

# Pinned to the design point this capture was run at (see
# mmi_2x2_waist_hybrid_vs_fdtd.py's own defaults).
mmi = MMI(n_inputs=2, n_outputs=2, mm_length=21.5, mm_width=2.8,
          waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
          waveguide_offset=1.5, waveguide_height=0.22,
          body_waist_width=2.24, body_waist_length=6.0,
          body_transition_length=6.0)

cap = load_captured_section(cap_path)
comp, mags = dominant_component(cap)
print(f"component magnitudes: {', '.join(f'{k}={v:.3e}' for k, v in mags.items())}"
      f"  -> using {comp}")

x, y = np.asarray(cap["x"]), np.asarray(cap["y"])
I = np.abs(np.asarray(cap[comp])) ** 2

fig, ax = plt.subplots(figsize=(14, 4.0))
im = ax.pcolormesh(x, y, I.T, shading="auto", cmap="inferno")
draw_geometry_outline(ax, mmi, x, color="lime")
ax.set_xlabel("device x — propagation (µm)")
ax.set_ylabel("device y — in-plane (µm)")
ax.set_ylim(-2.0, 2.0)
ax.set_title(
    f"waisted 2×2 MMI — FULL-DEVICE FDTD, |{comp}|²   λ = {float(cap['wavelength']):.4f} µm\n"
    "single-domain brute-force run: no hybrid seam, no EME approximation "
    "(solid lime = core outline)",
    fontsize=10)
plt.colorbar(im, ax=ax, pad=.01, label=f"|{comp}|²")
out = os.path.join(tapered_dir, f"fdtd_total_field_{comp}_res40.png")
fig.savefig(out, dpi=200, bbox_inches="tight")
plt.close("all")
print(f"saved -> {out}")
