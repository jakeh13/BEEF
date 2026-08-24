"""Tutorial: directional_coupler_fdtd.py's single-domain full-device FDTD
method (no EME decomposition), re-derived on SOI_160 at duty/dc.gds's real
dimensions, for the same 4 "DC_test" 4C10-die devices directional_coupler_
soi160_hybrid_dc_test.py characterizes via the hybrid method (see that
script's own module docstring for the full id-to-coupling_length derivation:
1901/1902 -> 11.634 um, 1911/1912 -> 3.0 um, 1921/1922 -> 0.5 um,
1931/1932 -> 0.0 um).

WHY A SEPARATE FULL-FDTD SCRIPT, AND WHY IT EXISTS AT ALL (2026-08-10):

directional_coupler_soi160_hybrid_dc_test.py's own docstring flags
coupling_length=0 (device 1931/1932) as unreachable through the 3-section
HybridSolver decomposition: with no straight coupling region at all, there
is nothing for the EME "body" section to propagate through, and that
script's own BODY_STUB math makes the body length negative and asserts
rather than produce a nonsense near-zero-length EME solve. Its docstring
names the fix directly: "a single end-to-end FDTD solve of the whole
bend-to-bend geometry, skipping the EME body section entirely" -- that is
exactly this script.

Confirmed safe at coupling_length=0 by reading the actual geometry/solver
code (not assumed): DirectionalCoupler.polygons() guards the straight
coupling rectangle with `if hc > 1e-9` (beef/geometry/coupler.py) and
total_length = coupling_length + 2*bend_length stays strictly positive
(2*bend_length = 20 um here) even at coupling_length=0 -- the two S-bends
just sit back-to-back. FDTDSolver.solve() (unlike solve_section()/
HybridSolver) never touches width_at()/center_at()/section_prisms() at
all, so none of the section-boundary math the hybrid script's assertion
guards against is even in this code path.

Same device, solver params, and wavelength window as directional_coupler_
soi160_hybrid_dc_test.py (SOI_160, duty/dc.gds's own waveguide_width=0.5/
gap=0.25/bend_length=10.0/waveguide_offset=4.0/bend_type="euler",
FDTD_RESOLUTION=40 ASSUMED not re-measured on this platform, same as that
script) -- the only method difference is the whole device sits in ONE FDTD
domain (directional_coupler_fdtd.py's method), not decomposed into
[FDTD bend | EME body | FDTD bend].

CACHE: checked directly against example_sim_cache/ before writing this
docstring (2026-08-10) -- the only existing whole-device FDTDSolver.solve()
cache entries for DirectionalCoupler are the SOI_220 directional_coupler_
fdtd.py device (waveguide_width=0.45, bend_length=9.0); none match this
script's SOI_160 dims at any coupling_length. Unlike the hybrid method
(which shares its 2 FDTD bend-section solves across every coupling_length,
see the hybrid script's own docstring), full single-domain FDTD's domain
covers the ENTIRE device including the coupling region, so a different
coupling_length is a different domain length end to end -- no cache reuse
is possible between this script's 4 lengths, or with the hybrid script's
cache. Every run below is a genuinely new FDTD solve.

Run:
    python /path/to/directional_coupler_soi160_fdtd_dc_test.py <coupling_length>
    python /path/to/directional_coupler_soi160_fdtd_dc_test.py 11.634
    python /path/to/directional_coupler_soi160_fdtd_dc_test.py 0.0
"""

import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import DirectionalCoupler, FDTDSolver, SOI_160, plot_fdtd_geometry

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("coupling_length", type=float,
                     help="Straight coupled-pair length in um (e.g. 11.634, 3.0, 0.5, or 0.0 -- "
                          "see module docstring's id-to-length table). Unlike the hybrid dc_test "
                          "script, 0.0 is fully supported here.")
args = parser.parse_args()
COUPLING_LENGTH = args.coupling_length

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "directional_coupler",
    "directional_coupler_soi160_fdtd_dc_test", f"L{COUPLING_LENGTH:.3f}"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters -- identical to directional_coupler_soi160_
# hybrid_dc_test.py's production values (itself ASSUMED from the SOI_220
# convergence study, not re-measured on this platform; see that script's docstring).
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN = 1.50
WL_MAX = 1.60
N_WL = 11

# DC160_WL_MIN/DC160_WL_MAX/DC160_N_WL override the full-band window above --
# same env-var family, same all-or-nothing guard, and the same
# "the guard is the env var, not the value" rule as
# directional_coupler_soi160_hybrid.py's own block, so the hybrid and
# full-FDTD methods can be pointed at one band with one pair of exports and
# compared apples-to-apples. Unset (the default) reproduces job 11785063's
# original full-band behavior exactly, output_dir unchanged; any override
# redirects to L<len>/windows/<tag>/ so the existing per-length
# coupling_result.npz/smatrix_fdtd.pkl cannot be overwritten.
#
# There is no cache reuse to be had either way here: the wavelength array is
# hashed into the FDTD fingerprint, AND (per the module docstring) a
# whole-device solve() at a different coupling_length is a different domain
# end to end. Every window x length combination is a genuinely new solve.
_win_env = (os.environ.get("DC160_WL_MIN"), os.environ.get("DC160_WL_MAX"), os.environ.get("DC160_N_WL"))
if any(v is not None for v in _win_env):
    if not all(v is not None for v in _win_env):
        raise SystemExit(
            "DC160_WL_MIN/DC160_WL_MAX/DC160_N_WL must be set together (or "
            "not at all) -- got a partial override, which would silently "
            "mix a custom window with the default WL_MIN/WL_MAX/N_WL above."
        )
    WL_MIN, WL_MAX, N_WL = float(_win_env[0]), float(_win_env[1]), int(_win_env[2])
    output_dir = os.path.join(output_dir, "windows", f"{WL_MIN:.4f}-{WL_MAX:.4f}_n{N_WL}")
    print(f"DC160 wavelength window override: WL_MIN={WL_MIN} WL_MAX={WL_MAX} "
          f"N_WL={N_WL} -> output_dir={output_dir}")
    os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Geometry -- same real duty/dc.gds dimensions as the hybrid dc_test script.
# ---------------------------------------------------------------------------
dc = DirectionalCoupler(
    waveguide_width=0.5,
    gap=0.25,
    coupling_length=COUPLING_LENGTH,
    bend_length=10.0,
    waveguide_offset=4.0,
    bend_type="euler",
    waveguide_height=0.16,
)
print(f"Directional coupler (SOI_160, duty/dc.gds dims, full-device FDTD, "
      f"coupling_length={COUPLING_LENGTH:.3f}) total length: {dc.total_length:.2f} um "
      f"(coupling pitch {dc.coupling_pitch:.3f} um, far pitch {dc.waveguide_offset:.2f} um)")

# ---------------------------------------------------------------------------
# 2. Configure the FDTD solver -- one domain around the whole device.
#
#    Output arms sit at y = +/- waveguide_offset/2 = +/- 2.0 um (vs SOI_220's
#    +/- 1.5); each output monitor is a half-domain window (sim_width/2 = 3.5
#    um wide, from the shared y=0 axis to the near PML edge) centered on its
#    own arm -- same convention as directional_coupler_fdtd.py, scaled to
#    this platform's sim_width=7.0 (bumped from SOI_220's 6.0 for the same
#    reason directional_coupler_soi160_hybrid.py's section solves are: each
#    arm needs waveguide_offset/2 + PML clearance).
# ---------------------------------------------------------------------------
solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=7.0,
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    num_input_modes=1,
    num_output_modes=1,
    # Order fixes the SMatrix port names: out_0 = bottom arm ("through", same
    # arm as "in", port="bot"), out_1 = top arm ("cross", the coupled arm).
    output_monitor_centers=[-2.0, +2.0],
    output_monitor_widths=[3.5, 3.5],
    material_stack=SOI_160,
    decay_tol=1e-11,
)

# ---------------------------------------------------------------------------
# 3. Geometry plot.
# ---------------------------------------------------------------------------
plot_fdtd_geometry(
    dc, solver=solver,
    port="bot",
    title=f"Directional Coupler (SOI_160, coupling_length={COUPLING_LENGTH:.3f} um) — Full-Device FDTD",
    output_path=os.path.join(output_dir, "geometry_fdtd.png"),
)
print(f"Geometry plot saved -> {output_dir}/geometry_fdtd.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep -- a single broadband run, entire wavelength window
#    resolved from one pulse's DFT.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nWavelength sweep: {wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um "
      f"({len(wavelengths)} pts, single FDTD run)\n")

smat = solver.solve(dc, wavelengths, port="bot", input_mode=1, save_dir=save_dir)
smat.check_passive()

# ---------------------------------------------------------------------------
# 5. Print summary at the centre wavelength.
# ---------------------------------------------------------------------------
T_through = smat.transmission(port="out_0", in_mode=1, out_mode=1)
T_cross = smat.transmission(port="out_1", in_mode=1, out_mode=1)
T_total = T_through + T_cross

wl_c = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_c))

print(f"Results at coupling_length = {COUPLING_LENGTH:.3f} um, lambda = {wl_c:.3f} um:")
print(f"  Through output T = {T_through[idx]:.4f}  IL = {10*np.log10(max(T_through[idx],1e-20)):.2f} dB")
print(f"  Cross output   T = {T_cross[idx]:.4f}  IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total          T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  Coupling ratio   = {T_cross[idx]/max(T_total[idx],1e-20):.3f}  "
      f"(0 = fully through, 1 = fully cross)")

# ---------------------------------------------------------------------------
# 6. Save and plot.
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix_fdtd.pkl"))
print(f"\nSMatrix saved -> {output_dir}/smatrix_fdtd.pkl")

np.savez(
    os.path.join(output_dir, "coupling_result.npz"),
    wavelengths=wavelengths, T_through=T_through, T_cross=T_cross, T_total=T_total,
    coupling_length=COUPLING_LENGTH, fdtd_resolution=RESOLUTION,
)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_through, label="Through", color="steelblue")
ax.plot(wavelengths, T_cross,   label="Cross",   color="darkorange", linestyle="--")
ax.plot(wavelengths, T_total,   label="Total",   color="k", lw=1.5)
ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Transmission")
ax.set_title(f"Transmission vs Wavelength (coupling_length={COUPLING_LENGTH:.3f} um)")
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
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title(f"Insertion Loss vs Wavelength (coupling_length={COUPLING_LENGTH:.3f} um)")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle(f"Directional Coupler (SOI_160, coupling_length={COUPLING_LENGTH:.3f} um) — "
             f"Single-Domain FDTD", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "fdtd_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved -> {output_dir}/fdtd_results.png")
