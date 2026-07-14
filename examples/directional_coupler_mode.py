"""
Tutorial: directional-coupler coupling length via even/odd supermode analysis
(meep's own eigenmode/MPB solver — no time-domain FDTD, no FECES solvers).

Method, per Ansys/Lumerical's "Evanescent waveguide couplers":
https://optics.ansys.com/hc/en-us/articles/360042304694-Evanescent-waveguide-couplers

    Two identical, individually single-mode waveguides placed close enough
    to interact evanescently no longer have independent modes — the coupled
    cross-section instead supports two "supermodes":

        even (symmetric)      — effective index n_eff,even  (higher)
        odd  (antisymmetric)  — effective index n_eff,odd   (lower)

    Launching light into just one waveguide excites an equal superposition
    of both supermodes. Because they propagate at slightly different phase
    velocities, power beats back and forth between the two physical
    waveguides as they co-propagate. The distance for a full 100% transfer
    from one waveguide to the other — the coupling length — is

        L_c = lambda / (2 * |n_eff,even - n_eff,odd|)

    and the fraction of power transferred to the other waveguide after
    propagating a length L is

        T_cross(L) = sin^2(pi * L / (2 * L_c))

This script builds the SAME cross-section as the coupling region of
directional_coupler_hybrid.py (220 nm SOI, 0.45 um wide strips, 0.25 um
gap), solves for the two lowest-order eigenmodes of that cross-section with
meep's own mode solver (meow.compute_modes_meep, which wraps
mp.Simulation.get_eigenmode — i.e. MPB) at each wavelength, and reports
L_c(wavelength). It also sweeps the gap at a fixed wavelength to show the
classic exponential falloff of coupling strength with separation.

This is dramatically cheaper than a full propagation sweep (no time-domain
FDTD, no length dependence to resimulate at all — coupling length is a
per-wavelength/per-gap CROSS-SECTION property, solved once each) which is
exactly why the "solve two supermodes and take the beat length" method is
the standard way to characterize a directional coupler.

Only reuses feces.materials.SOI_220 for the Si/SiO2 dispersion + extrusion
rules (so the material model matches directional_coupler_hybrid.py exactly)
— everything else (GDS layout, meow calls, mode solving) is done directly,
without going through any FECES solver class.

Run:
    python /path/to/examples/directional_coupler_mode.py
"""

import itertools
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gdsfactory as gf
import meow as mw

from feces.materials import SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_output", "directional_coupler_mode"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Geometry / material constants — matches directional_coupler_hybrid.py's
#    coupling-region cross-section (220 nm SOI).
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 0.45     # um, single-mode strip width
GAP_NOMINAL = 0.25         # um, edge-to-edge coupling gap used for the wavelength sweep
MATERIAL_STACK = SOI_220   # 220 nm Si on SiO2 — same stack as the hybrid example

MESH_WIDTH = 4.0           # um, lateral (y) extent of the mode-solve mesh
MESH_HEIGHT = 2.0          # um, vertical (z) extent — matches MaterialStack.sz
RES = 40.0                 # pts/um, isotropic (required — meep's resolution is a single scalar)

WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 21   # matches directional_coupler_hybrid.py's C-band sweep
GAPS = [0.15, 0.20, 0.25, 0.30, 0.40, 0.50]   # um, secondary sweep at the center wavelength

_ctr = itertools.count()


# ---------------------------------------------------------------------------
# 2. Cross-section builder — two parallel silicon strips, `gap` apart.
# ---------------------------------------------------------------------------
def build_cross_section(gap: float, wl: float) -> mw.CrossSection:
    """gdsfactory GDS -> meow CrossSection for two coupled parallel strips.

    The GDS "length" here is arbitrary (1 um) — a CrossSection is a single
    transverse slice, so nothing about it depends on propagation length.
    """
    pitch = WAVEGUIDE_WIDTH + gap
    yc0, yc1 = -pitch / 2, +pitch / 2
    length = 1.0

    gf.gpdk.PDK.activate()
    c = gf.Component(f"coupler_xsec_{next(_ctr)}")
    for yc in (yc0, yc1):
        c.add_polygon([
            (0, yc - WAVEGUIDE_WIDTH / 2), (length, yc - WAVEGUIDE_WIDTH / 2),
            (length, yc + WAVEGUIDE_WIDTH / 2), (0, yc + WAVEGUIDE_WIDTH / 2),
        ], layer=(2, 0))
    c.add_polygon([
        (0, -MESH_WIDTH / 2), (length, -MESH_WIDTH / 2),
        (length, MESH_WIDTH / 2), (0, MESH_WIDTH / 2),
    ], layer=(1, 0))

    structs = mw.extrude_gds(c, MATERIAL_STACK.extrusion_rules)
    mesh = mw.Mesh2D(
        x=np.linspace(-MESH_WIDTH / 2, MESH_WIDTH / 2, int(RES * MESH_WIDTH) + 1),
        y=np.linspace(-MESH_HEIGHT / 2, MESH_HEIGHT / 2, int(RES * MESH_HEIGHT) + 1),
    )
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=float(wl), T=25.0)
    return mw.CrossSection.from_cell(cell=cells[0], env=env)


def supermode_indices(gap: float, wl: float):
    """Even/odd supermode effective indices via meep's own mode solver (MPB).

    compute_modes_meep already returns modes sorted by decreasing
    real(neff); for two coupled single-mode waveguides the top two bands
    are exactly the even (higher index, symmetric) and odd (lower index,
    antisymmetric) supermodes.
    """
    from meow import compute_modes_meep

    cs = build_cross_section(gap, wl)
    modes = compute_modes_meep(cs, num_modes=2)
    n_even, n_odd = sorted((float(np.real(m.neff)) for m in modes), reverse=True)
    return n_even, n_odd


def coupling_length(n_even: float, n_odd: float, wl: float) -> float:
    return wl / (2.0 * abs(n_even - n_odd))


# ---------------------------------------------------------------------------
# 3. Wavelength sweep at the nominal gap.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"Directional coupler supermode analysis — {WAVEGUIDE_WIDTH} um strips, "
      f"{GAP_NOMINAL} um gap, 220 nm SOI")
print(f"Wavelength sweep: {wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um ({len(wavelengths)} pts)\n")

n_even_wl = np.zeros(N_WL)
n_odd_wl = np.zeros(N_WL)
for i, wl in enumerate(wavelengths):
    n_even_wl[i], n_odd_wl[i] = supermode_indices(GAP_NOMINAL, wl)

Lc_wl = wavelengths / (2.0 * np.abs(n_even_wl - n_odd_wl))

wl_c = wavelengths[len(wavelengths) // 2]
idx_c = int(np.argmin(np.abs(wavelengths - wl_c)))
print(f"At lambda = {wl_c:.3f} um: n_eff,even = {n_even_wl[idx_c]:.5f}, "
      f"n_eff,odd = {n_odd_wl[idx_c]:.5f}, delta_n = {n_even_wl[idx_c]-n_odd_wl[idx_c]:.5f}")
print(f"  Coupling length L_c = {Lc_wl[idx_c]:.3f} um  "
      f"(100% cross-over at L = L_c; back to 'through' at L = 2*L_c)\n")

# ---------------------------------------------------------------------------
# 4. Gap sweep at the center wavelength.
# ---------------------------------------------------------------------------
print(f"Gap sweep at lambda = {wl_c:.3f} um: {GAPS} um")
n_even_gap = np.zeros(len(GAPS))
n_odd_gap = np.zeros(len(GAPS))
for i, gap in enumerate(GAPS):
    n_even_gap[i], n_odd_gap[i] = supermode_indices(gap, wl_c)

Lc_gap = wl_c / (2.0 * np.abs(n_even_gap - n_odd_gap))
for gap, Lc in zip(GAPS, Lc_gap):
    print(f"  gap = {gap:.2f} um  ->  L_c = {Lc:.3f} um")

# ---------------------------------------------------------------------------
# 5. Save data
# ---------------------------------------------------------------------------
np.savez(
    os.path.join(output_dir, "supermode_sweep.npz"),
    wavelengths=wavelengths, n_even_wl=n_even_wl, n_odd_wl=n_odd_wl, Lc_wl=Lc_wl,
    gaps=np.asarray(GAPS), n_even_gap=n_even_gap, n_odd_gap=n_odd_gap, Lc_gap=Lc_gap,
    waveguide_width=WAVEGUIDE_WIDTH, gap_nominal=GAP_NOMINAL,
)
print(f"\nData saved -> {output_dir}/supermode_sweep.npz")

# ---------------------------------------------------------------------------
# 6. Plots
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

ax = axes[0]
ax.plot(wavelengths, n_even_wl, label="n_eff, even (symmetric)", color="steelblue")
ax.plot(wavelengths, n_odd_wl, label="n_eff, odd (antisymmetric)", color="darkorange", linestyle="--")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Effective index")
ax.set_title(f"Supermode indices (gap={GAP_NOMINAL} um)")
ax.legend(fontsize=8)
ax.grid(True)

ax = axes[1]
ax.plot(wavelengths, Lc_wl, color="k", lw=1.5)
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Coupling length L_c (um)")
ax.set_title(f"Coupling length vs wavelength (gap={GAP_NOMINAL} um)")
ax.grid(True)

ax = axes[2]
ax.semilogy(GAPS, Lc_gap, "o-", color="firebrick")
ax.set_xlabel("Gap (um)")
ax.set_ylabel("Coupling length L_c (um)")
ax.set_title(f"Coupling length vs gap (lambda={wl_c:.2f} um)")
ax.grid(True, which="both")

plt.suptitle("Directional Coupler — Even/Odd Supermode Analysis (meep MPB)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "supermode_analysis.png"), dpi=300)
plt.close("all")
print(f"Supermode analysis plot saved -> {output_dir}/supermode_analysis.png")

# -- Predicted power transfer vs length, at the nominal gap/center wavelength,
#    over the same 0-10 um range directional_coupler_hybrid.py sweeps — lets
#    you sanity-check that script's full FDTD/EME result against this much
#    cheaper cross-section-only prediction.
L_range = np.linspace(0.0, 60.0, 201)
T_cross_pred = np.sin(np.pi * L_range / (2.0 * Lc_wl[idx_c])) ** 2
T_through_pred = 1.0 - T_cross_pred

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(L_range, T_through_pred, label="Through (predicted)", color="steelblue")
ax.plot(L_range, T_cross_pred, label="Cross (predicted)", color="darkorange", linestyle="--")
ax.axvline(Lc_wl[idx_c], color="gray", linestyle=":", lw=0.8,
           label=f"L_c = {Lc_wl[idx_c]:.2f} um")
ax.set_xlabel("Coupling length (um)")
ax.set_ylabel("Predicted power transfer")
ax.set_title(f"Predicted Through/Cross Power vs Length (lambda={wl_c:.2f} um, gap={GAP_NOMINAL} um)\n"
             f"T_cross(L) = sin^2(pi*L / (2*L_c)) — compare against directional_coupler_hybrid.py's "
             f"coupling_length_sweep.png")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "predicted_power_transfer.png"), dpi=300)
plt.close("all")
print(f"Predicted power-transfer plot saved -> {output_dir}/predicted_power_transfer.png")
