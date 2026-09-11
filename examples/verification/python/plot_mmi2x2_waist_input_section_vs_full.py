"""Input-section field: hybrid's FDTD section vs the SAME span of the full device.

The hybrid-vs-FDTD comparison (jobs 12892459 / 12797278) agrees to ~3% in
|dT| but carries a FLAT ~1.5% splitting-ratio offset -- flat across the whole
band, not growing with wavelength. A dispersion-free bias does not look like an
uncorrected body beta (which accumulates with propagation); it looks like
something at the interfaces. This is the direct test of the first interface.

Both maps come from FDTD at res 40, wl 1.5475 um, identical geometry:
  captured_input_section : solve_section over [-14.750, -10.500]  (the hybrid's
                           own section 0, half_domain_aperture=True)
  captured_full_device   : solve_section over [-14.750, +14.750]  (whole device,
                           no interior cut)
Cropping the full-device map to the section's x-range makes the comparison
exact -- same physics, same span, one cut at the seam and one not.

⚠ The two runs differ in TWO ways beyond the cut, and both are stated in the
figure: the section run uses half_domain_aperture (it must, to isolate one of
two closely-spaced inputs) and carries n_modes_out=14 monitors at its right
boundary, where the full-device run measures a single output mode. So a
difference here is "sectioning as configured" vs "not sectioned", not a clean
one-variable experiment.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                    "example_output", "mmi_2x2_tapered"))
WL = 1.5475
sec = np.load(os.path.join(BASE, f"captured_input_section_res40_wl{WL:.4f}",
                           "section_fields_band1.npz"), allow_pickle=True)
full = np.load(os.path.join(BASE, f"captured_full_device_res40_wl{WL:.4f}",
                            "section_fields_band1.npz"), allow_pickle=True)

COMP = os.environ.get("MMI2X2_CMP_COMP", "Ey")   # dominant TE for this device
xs, ys, Fs = sec["x"], sec["y"], np.abs(sec[COMP])
xf, yf, Ff = full["x"], full["y"], np.abs(full[COMP])

# crop the full-device map to the section's span
m = (xf >= xs.min() - 1e-9) & (xf <= xs.max() + 1e-9)
xf_c, Ff_c = xf[m], Ff[m]
print(f"section span  {xs.min():.3f}..{xs.max():.3f} um, {len(xs)} x-samples")
print(f"full  cropped {xf_c.min():.3f}..{xf_c.max():.3f} um, {len(xf_c)} x-samples")

# normalise each to its own max -- absolute scales differ (different source
# normalisation and monitor sets), so the comparison is of SHAPE
Fs_n, Ff_n = Fs / Fs.max(), Ff_c / Ff_c.max()

fig, ax = plt.subplots(3, 1, figsize=(9, 9))
vm = 1.0
for a, (F, x, ttl) in zip(ax[:2], [
        (Fs_n, xs, "hybrid input SECTION (cut at the seam, half-domain aperture)"),
        (Ff_n, xf_c, "FULL DEVICE, same x-range (no interior cut)")]):
    im = a.pcolormesh(x, ys, F.T, shading="auto", cmap="magma", vmin=0, vmax=vm)
    a.set_title(f"|{COMP}| — {ttl}", fontsize=9)
    a.set_ylabel("y (µm)")
    plt.colorbar(im, ax=a, label=f"|{COMP}| (norm.)")
    for xb in (-10.5,):
        a.axvline(xb, color="cyan", ls="--", lw=1)

# on-axis-of-each-port line cuts, so a shape difference is quantitative
iy_top = int(np.argmin(np.abs(ys - 0.75)))
iy_bot = int(np.argmin(np.abs(ys + 0.75)))
ax[2].plot(xs, Fs_n[:, iy_top], "-", color="C0", label="section, y=+0.75 (excited)")
ax[2].plot(xf_c, Ff_n[:, iy_top], "--", color="C0",
           label="full device, y=+0.75")
ax[2].plot(xs, Fs_n[:, iy_bot], "-", color="C3", label="section, y=−0.75 (neighbour)")
ax[2].plot(xf_c, Ff_n[:, iy_bot], "--", color="C3",
           label="full device, y=−0.75")
ax[2].axvline(-10.5, color="cyan", ls="--", lw=1, label="seam (body face + stub)")
ax[2].set_xlabel("x (µm)")
ax[2].set_ylabel(f"|{COMP}| (norm.)")
ax[2].legend(fontsize=7, ncol=2)
ax[2].grid(alpha=0.3)
fig.suptitle("Waisted 2×2 — input-section field, sectioned vs full device "
             f"(FDTD, res 40, λ={WL} µm)\n"
             "cyan dashed = the hybrid's section-0/body seam", fontsize=10)
plt.tight_layout(rect=[0, 0, 1, 0.96])
p = os.path.join(BASE, f"input_section_vs_full_{COMP}_res40.png")
plt.savefig(p, dpi=160)
plt.close()
print(f"written -> {p}")

# quantitative: how different are the two on the shared grid?
if len(xs) == len(xf_c):
    d = np.abs(Fs_n - Ff_n)
    print(f"\nshape difference on the shared span: max {d.max():.4f}, "
          f"mean {d.mean():.4f} (of a max-normalised field)")
else:
    print(f"\ngrids differ in length ({len(xs)} vs {len(xf_c)}); "
          "line cuts above are the comparison")
