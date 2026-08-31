"""Output-waveguide crosstalk check for the 1x4 MMI on 300 nm Si3N4 --
validation_markdown/mmi_1x4_sin300_validation.md step 6b/6c/6d's own flagged
concern, checked in parallel with step 5's mm_length sweep (job 11667114)
since it does not depend on mm_length.

The plan's step 6b note: "the 'comfortably decoupled' justification in the
1x2's comment (waveguide_offset - waveguide_width = 0.65 um edge-to-edge)
must be RECOMPUTED for nitride, where a weakly-confined mode couples across
a gap that would be safe in silicon."

Method: directional_coupler_mode.py's even/odd supermode analysis (Ansys/
Lumerical "Evanescent waveguide couplers") applied to the ACTUAL 1x4 output
geometry -- two SIN_300 strips, waveguide_width=1.1 um, edge-to-edge gap =
waveguide_offset - waveguide_width = 2.4 - 1.1 = 1.3 um. Compute the
coupling length L_c = lambda / (2*|n_eff,even - n_eff,odd|), compare against
the actual propagation length the four output waveguides sit this close
together in the modeled device: taper_length (5.0 um) + the short access
stub each solver adds (0.5 um EME / body_stub-scale FDTD) ~= 5.5 um. If
L_c >> 5.5 um, crosstalk over that span is negligible (T_cross =
sin^2(pi*L/(2*L_c)) stays small) and the existing waveguide_offset=2.4 design
does not need to change for this reason.

mode_solver: compute_modes_tidy3d, NOT compute_modes_meep despite the
original directional_coupler_mode.py using meep -- step 4 found
compute_modes_meep silently wrong for this platform's weakly-confined modes.
A coupled well-guided single-mode pair (1.1 um wide, 2 guided modes each in
isolation) is a less extreme case than the near-cutoff window that broke it,
but tidy3d is used throughout this device's design work for consistency and
safety, not re-litigated per script.

Run:
    python /path/to/examples/verification/python/mmi_1x4_sin_output_crosstalk.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import meow as mw
from meow import compute_modes_tidy3d

from beef import extrude_polygons
from beef.materials.stack import MaterialStack

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification", "mmi_1x4_sin_output_crosstalk"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Geometry / material constants -- matches the 1x4 MMI's actual output
# access-waveguide design (steps 4-5).
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 1.1
# ⚠ Constants refreshed 2026-08-31 to the post-step-9c geometry
# (mm_length 206.0 -> 39.90, waveguide_offset 2.4 -> 3.1, set
# 2026-08-14). Any output in example_output/ predating that date was
# produced with the OLD values and does not correspond to this file.
#
# NOTE the direction this moves the conclusion: the gap widens from
# 2.4-1.1 = 1.3 um to 3.1-1.1 = 2.0 um, so the recorded pre-9c crosstalk
# result is CONSERVATIVE for the current device, not optimistic.
WAVEGUIDE_OFFSET = 3.1
GAP_NOMINAL = WAVEGUIDE_OFFSET - WAVEGUIDE_WIDTH   # 1.3 um, edge-to-edge
SZ = 5.0   # step 4's converged MaterialStack.sz
STACK = MaterialStack(wg_height=0.30, sz=SZ, core_kind="si3n4")

MESH_WIDTH = 8.0    # um, lateral mesh extent -- generous past two 1.1um
                    # strips at a 3.5um pitch plus cladding decay length
RES = 40.0          # pts/um, this codebase's standing mode-solve default

WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11   # matches this device's C-band sweep
GAPS = [0.5, 0.8, 1.0, 1.3, 1.6, 2.0]    # um, secondary sweep at center wavelength;
                                          # 1.3 is this device's actual gap

# The propagation length these 4 waveguides actually sit at close (<
# waveguide_offset) pitch in the modeled device: the taper band, plus the
# short straight stub each script's sections add on the far side of it.
RELEVANT_LENGTH = 5.0 + 0.5   # taper_length + STRAIGHT_LENGTH, mmi_1x4_sin_body_length_sweep.py


def _axis_grid(half_span: float, res: float) -> np.ndarray:
    n_half = int(round(half_span * res))
    step = 1.0 / res
    return step * np.arange(-n_half, n_half + 1)


def build_cross_section(gap: float, wl: float) -> mw.CrossSection:
    """Two parallel SIN_300 strips, `gap` apart (edge-to-edge) -- the
    coupled-waveguide analogue of build_strip_cross_section() in
    mmi_1x4_sin_cross_section.py (which builds one isolated strip)."""
    pitch = WAVEGUIDE_WIDTH + gap
    yc0, yc1 = -pitch / 2, +pitch / 2
    length = 1.0

    wg_polys = [
        [(0, yc - WAVEGUIDE_WIDTH / 2), (length, yc - WAVEGUIDE_WIDTH / 2),
         (length, yc + WAVEGUIDE_WIDTH / 2), (0, yc + WAVEGUIDE_WIDTH / 2)]
        for yc in (yc0, yc1)
    ]
    ox_poly = [
        (0, -MESH_WIDTH / 2), (length, -MESH_WIDTH / 2),
        (length, MESH_WIDTH / 2), (0, MESH_WIDTH / 2),
    ]
    gds_polys = {(2, 0): wg_polys, (1, 0): [ox_poly]}

    structs = extrude_polygons(gds_polys, STACK.extrusion_rules)
    mesh = mw.Mesh2D(x=_axis_grid(MESH_WIDTH / 2, RES), y=_axis_grid(SZ / 2, RES))
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=float(wl), T=25.0)
    return mw.CrossSection.from_cell(cell=cells[0], env=env)


def supermode_indices(gap: float, wl: float):
    """Even/odd supermode effective indices via compute_modes_tidy3d (step
    4's fix -- compute_modes_meep silently wrong for this platform)."""
    cs = build_cross_section(gap, wl)
    modes = compute_modes_tidy3d(cs, num_modes=2)
    n_even, n_odd = sorted((float(np.real(m.neff)) for m in modes), reverse=True)
    return n_even, n_odd


def coupling_length(n_even: float, n_odd: float, wl: float) -> float:
    return wl / (2.0 * abs(n_even - n_odd))


# ---------------------------------------------------------------------------
# Wavelength sweep at the device's actual gap.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"1x4 MMI output-waveguide crosstalk check -- {WAVEGUIDE_WIDTH} um strips, "
      f"{GAP_NOMINAL:.2f} um edge-to-edge gap (waveguide_offset={WAVEGUIDE_OFFSET}), "
      f"300 nm Si3N4 (sz={SZ})")
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
Lc_c = Lc_wl[idx_c]
T_cross_relevant = np.sin(np.pi * RELEVANT_LENGTH / (2.0 * Lc_c)) ** 2
print(f"  Coupling length L_c = {Lc_c:.3f} um")
print(f"  Relevant propagation length (taper + stub) = {RELEVANT_LENGTH:.2f} um "
      f"= {RELEVANT_LENGTH/Lc_c:.4f} x L_c")
print(f"  T_cross at that length = sin^2(pi*L/(2*L_c)) = {T_cross_relevant:.6f}")
if T_cross_relevant < 1e-3:
    print("  -> NEGLIGIBLE crosstalk over the modeled taper+stub span; "
          "waveguide_offset=2.4 um does not need revision for this reason.")
else:
    print("  -> NOT obviously negligible -- reconsider waveguide_offset or "
          "flag as a design risk before step 6.")

# ---------------------------------------------------------------------------
# Gap sweep at the center wavelength -- shows the falloff, brackets the
# device's actual gap=1.3 um against both tighter and looser spacings.
# ---------------------------------------------------------------------------
print(f"\nGap sweep at lambda = {wl_c:.3f} um: {GAPS} um")
n_even_gap = np.zeros(len(GAPS))
n_odd_gap = np.zeros(len(GAPS))
for i, gap in enumerate(GAPS):
    n_even_gap[i], n_odd_gap[i] = supermode_indices(gap, wl_c)

Lc_gap = wl_c / (2.0 * np.abs(n_even_gap - n_odd_gap))
for gap, Lc in zip(GAPS, Lc_gap):
    T_cross = np.sin(np.pi * RELEVANT_LENGTH / (2.0 * Lc)) ** 2
    marker = "  <-- this device's actual gap" if abs(gap - GAP_NOMINAL) < 1e-9 else ""
    print(f"  gap = {gap:.2f} um  ->  L_c = {Lc:9.2f} um   "
          f"T_cross({RELEVANT_LENGTH:.1f} um) = {T_cross:.6f}{marker}")

# ---------------------------------------------------------------------------
# Save data + plot
# ---------------------------------------------------------------------------
np.savez(
    os.path.join(output_dir, "crosstalk_sweep.npz"),
    wavelengths=wavelengths, n_even_wl=n_even_wl, n_odd_wl=n_odd_wl, Lc_wl=Lc_wl,
    gaps=np.asarray(GAPS), n_even_gap=n_even_gap, n_odd_gap=n_odd_gap, Lc_gap=Lc_gap,
    waveguide_width=WAVEGUIDE_WIDTH, gap_nominal=GAP_NOMINAL,
    relevant_length=RELEVANT_LENGTH,
)
print(f"\nData saved -> {output_dir}/crosstalk_sweep.npz")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

ax = axes[0]
ax.plot(wavelengths, Lc_wl, color="steelblue", marker="o", ms=3)
ax.axhline(RELEVANT_LENGTH, color="gray", linestyle=":", lw=1.2,
           label=f"taper+stub length ({RELEVANT_LENGTH:.1f} um)")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Coupling length L_c (um)")
ax.set_title(f"L_c vs wavelength, gap={GAP_NOMINAL:.2f} um")
ax.legend(fontsize=8)
ax.grid(True)

ax = axes[1]
ax.semilogy(GAPS, Lc_gap, color="darkorange", marker="o", ms=4)
ax.axvline(GAP_NOMINAL, color="gray", linestyle=":", lw=1.2,
           label=f"actual gap ({GAP_NOMINAL:.2f} um)")
ax.axhline(RELEVANT_LENGTH, color="gray", linestyle="--", lw=1.0,
           label=f"taper+stub length ({RELEVANT_LENGTH:.1f} um)")
ax.set_xlabel("Edge-to-edge gap (um)")
ax.set_ylabel("Coupling length L_c (um)")
ax.set_title(f"L_c vs gap, lambda={wl_c:.3f} um")
ax.legend(fontsize=8)
ax.grid(True, which="both")

plt.suptitle("1x4 MMI output-waveguide crosstalk -- 300 nm Si3N4", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "crosstalk_sweep.png"), dpi=300)
plt.close("all")
print(f"Plot saved -> {output_dir}/crosstalk_sweep.png")
