"""Tutorial: directional_coupler_soi160_hybrid.py at an arbitrary coupling_length,
for the OTHER real "DC_test" devices on the 4C10 die (GDS labels 1901/1902,
1921/1922, 1931/1932 -- 1911/1912 IS the device directional_coupler_soi160_
hybrid.py already characterizes at coupling_length=3.0, see below).

Same device, solver, wavelength window, and platform (SOI_160, duty/dc.gds's
own waveguide_width=0.5/gap=0.25/bend_length=10.0/waveguide_offset=4.0/
bend_type="euler") as directional_coupler_soi160_hybrid.py -- the only
difference is coupling_length is a command-line argument instead of the
hardcoded 3.0 um.

WHY THESE 4 DEVICES, AND HOW THEIR coupling_length VALUES WERE DERIVED
(2026-08-10, manuscript_figures/fp_resonator_figures session):

The 4C10 die (examples/gds_files-style layout, but this die-level GDS lives
at .../4C10_measurements/gds/4C10_right.gds on the shared OneDrive, not in
this repo) carries 4 "DC_test" device pairs at GDS labels 1901/1902,
1911/1912, 1921/1922, 1931/1932 (see 4C10_right.json in that same folder).
Each was measured twice, under two DIFFERENT LabExT id schemes:
  - fab-time IL_sweep ids (1901/1902/1911/1912/1921/1922/1931/1932),
    folders right_chip_0/, right_chip_1/, C_right_chip_1/.
  - a later re-measurement under LUNA-sweep ids 3913/3914/3923/3924 for the
    HCG Fabry-Perot resonator characterization specifically -- these are
    the ones already normalized/plotted in manuscript_figures/
    fp_resonator_figures/HCG_plotting_LUNA.ipynb section 2, and the ones
    directional_coupler_soi160_hybrid.py already simulates at
    coupling_length=3.0.

Auto-measuring each device's real coupling_length from the die GDS directly
(gdstk, same method Step 0 used on the standalone dc.gds cutout: find the
two waveguide-width-scale rectangles sharing an x- or y-span with a small
edge-to-edge gap) worked for 3 of the 4 -- confirmed against known/expected
values -- but assigning which of several nearby candidate values belongs to
WHICH of the 4 labels was not fully reliable from geometry alone (the 4
devices sit close enough together that spatial nearest-label matching
overlapped between them, and heater/electrode comb structures on the same
GDS layer produced false candidates that needed filtering by rectangle
width).

The id-to-length MAPPING actually used below was instead confirmed by
comparing MEASURED spectral shape: each device's IL_sweep AIN0/AIN1 dB
split (loopback-normalized) was plotted, and the SIZE of that split at a
FIXED gap=0.25 monotonically encodes where coupling_length sits on the
sin^2(kappa*L) coupling oscillation --near a peak (~50/50) means the
smallest AIN0/AIN1 split, near a trough (weak coupling) means the largest.
1911/1912's own measured split (~10 dB) matched the ALREADY-KNOWN
coupling_length=3.0 device's own measured split (9.3-9.4 dB, from the LUNA
data) almost exactly -- strong evidence 1911/1912 simply IS that same
device under its original fab-time label. The other 3, ordered by
INCREASING split (== decreasing coupling_length, since larger split means
further from the near-peak point 1901/1902 sits at):

    device id     measured AIN0/AIN1 split    coupling_length
    1901/1902     ~1.3 dB   (near 50/50)       11.634 um
    1911/1912     ~10.0 dB  (matches 3913/3914/3923/3924)  3.0 um  <- already done
    1921/1922     ~14.6 dB                     0.5 um
    1931/1932     ~17.7 dB  (most asymmetric)   0 um (or very near it)

This mapping is a MEASUREMENT-BASED inference, not a re-derived GDS
dimension read the way Step 0's original 3.0 um value was -- treat the
11.634/0.5/0 values here as good working hypotheses, not verified GDS
numbers, unless/until someone re-derives them directly from the die GDS
with a location-disambiguated read.

BODY_STUB is NOT the fixed 0.25 um directional_coupler_soi160_hybrid.py
uses -- that value only leaves a positive EME body length for
coupling_length=3.0 (2*BODY_STUB=0.5 out of 3.0). For the two short target
lengths (0.5, ~0) that would make body_length negative. This script instead
scales BODY_STUB down for short coupling_length (see BODY_STUB_FRACTION
below) and asserts positivity exactly like the original script does.

coupling_length=0 (device 1931/1932) is a genuine edge case: with no
straight coupling region at all, there is nothing for the EME "body"
section to propagate through -- the 3-section HybridSolver architecture
(FDTD bend | EME body | FDTD bend) this script and
directional_coupler_soi160_hybrid.py both use doesn't have a length-0
degenerate form. Run coupling_length=0 through THIS script and it will
raise the "body_length must be positive" assertion deliberately, rather
than silently produce a nonsense near-zero-length EME solve. A true
coupling_length=0 result likely needs a different construction entirely
(e.g. a single end-to-end FDTD solve of the whole bend-to-bend geometry,
skipping the EME body section) -- not attempted here; flag this to
whoever runs the id-1931/1932 case before spending compute on it.

coupling_length=0.5 (device 1921/1922) is a SECOND, separate edge case,
found the hard way (job 11785062_1, 2026-08-10): it clears the body_length>0
guard above (body_length=0.35 um with the *0.15 BODY_STUB formula), but each
FDTD section's own domain is padded by dpml=0.5 um past its body-facing
boundary regardless of BODY_STUB (solve_section() draws stub+dpml, not just
stub), and that padding alone (0.5 um) is already bigger than the whole
0.35 um body -- it spills into the OPPOSITE S-bend. Without a guard this
surfaces as a ValueError from deep inside beef/geometry/coupler.py's own
_check_cut() (a real geometry-construction failure, not this script's own
message) after the 24-core allocation is already running. The
_body_stub_cap assertion below catches it up front instead, same pattern as
the coupling_length=0 case: coupling_length must exceed roughly 2*dpml for
this 3-section decomposition to have a valid interior interface AT ALL, not
just a formula tweak away -- directional_coupler_soi160_fdtd_dc_test.py
(full single-domain FDTD, no such constraint) already has a real result for
this device, job 11785063.

Each invocation writes to its own length-tagged subfolder, mirroring
mmi_2x2_fdtd_alt_length.py's convention:

    example_output/directional_coupler/directional_coupler_soi160_hybrid_dc_test/
        L<coupling_length>/
            geometry_hybrid.png, sections/, smatrix_hybrid.pkl,
            coupling_result.npz, coupling_pieces.npz, hybrid_results.png

CACHE NOTE (checked 2026-08-10, via solve_section(..., cache_only=True) --
zero real computation, confirmed safe): the FDTD S-bend section cache keys
on ABSOLUTE x0/x1 position, and DirectionalCoupler's geometry is built
symmetric about the coupling-region CENTRE, so a different coupling_length
shifts total_length/2 and therefore x0/x1 even though the S-bend's own
LOCAL shape is unchanged. Confirmed empirically: coupling_length=11.634 is
a cache MISS against the existing example_sim_cache/. Every run of this
script for a new length is therefore a REAL new FDTD solve, same
2-FDTD-solve + 1-EME-mode-solve-per-wavelength cost as the original
coupling_length=3.0 characterization (which ran as a real SLURM batch job,
`-t 720`, per validation_markdown/running_jobs.md) -- run this via the
cluster batch queue (see the `directional_coupler_soi160_hybrid_dc_test`
Makefile target and a `COUPLING_LENGTH=<value> make directional_coupler_
soi160_hybrid_dc_test`-style batch array), not on a login node.

Run:
    python /path/to/directional_coupler_soi160_hybrid_dc_test.py <coupling_length>
    python /path/to/directional_coupler_soi160_hybrid_dc_test.py 11.634
    python /path/to/directional_coupler_soi160_hybrid_dc_test.py 0.5
"""

import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import (
    DirectionalCoupler, EMESolver, FDTDSolver, HybridSolver, SMatrix, SOI_160, plot_hybrid_geometry,
)

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("coupling_length", type=float,
                     help="Straight coupled-pair length in um (e.g. 11.634 or 0.5 -- see module docstring's id-to-length table)")
args = parser.parse_args()
COUPLING_LENGTH = args.coupling_length

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "directional_coupler",
    "directional_coupler_soi160_hybrid_dc_test", f"L{COUPLING_LENGTH:.3f}"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters -- identical to directional_coupler_soi160_
# hybrid.py's production values (ASSUMED from the SOI_220 convergence
# study, not re-measured on this platform; see that script's own docstring).
# ---------------------------------------------------------------------------
FDTD_RESOLUTION = 40.0
EME_RESOLUTION = 80.0
WL_MIN = 1.50
WL_MAX = 1.60
N_WL = 11

# BODY_STUB scales down for short coupling_length instead of the original
# script's fixed 0.25 um (which only leaves a positive EME body length for
# coupling_length=3.0 -- see module docstring). Kept comfortably below
# coupling_length/2; the assertion below is the same hard guard the
# original script uses, not a new check.
#
# ALSO capped by _body_stub_cap (added 2026-08-10 after job 11785062_1's
# real crash at coupling_length=0.5): each FDTD section's own domain is
# padded by _DPML past its body-facing boundary regardless of BODY_STUB --
# solve_section() draws stub+dpml, not just stub, so the waveguide/PML
# continues through the full PML thickness (see plot_hybrid_geometry's
# draw_stub_left/right, beef/plot/geometry.py, which mirrors this exactly
# for the section preview plots) -- so that padding alone reaches
# BODY_STUB - coupling_length/2 + _DPML past the coupling-region CENTRE on
# each side. At coupling_length=0.5, the *0.15 formula alone gives
# BODY_STUB=0.075, and the "in" section's padded right edge lands at
# x=+0.325 -- inside the (+0.25, +10.25) right S-bend, not the uniform body
# it was supposed to land in. beef/geometry/coupler.py's own _check_cut()
# caught it (a real ValueError, not a silently wrong answer) but from deep
# inside geometry-construction code, not this script's own guard -- the
# _body_stub_cap assertion below turns that into a clear, expected failure
# at the top of the script instead, the same pattern coupling_length=0
# already gets.
_DPML = 0.5   # must match _fdtd_common's dpml below -- single source of truth
_STUB_MARGIN = 0.05   # um -- keep the padding safely clear of the opposite
                      # bend rather than exactly on the knife-edge boundary
_body_stub_cap = COUPLING_LENGTH - _DPML - _STUB_MARGIN
assert _body_stub_cap > 0, (
    f"coupling_length ({COUPLING_LENGTH:.3f} um) leaves no room for the FDTD "
    f"sections' own {_DPML:.2f} um PML padding on the body-facing side of "
    f"EACH section without the two paddings overlapping through the middle "
    f"-- needs coupling_length > 2*dpml (+ margin) = "
    f"{2 * _DPML + _STUB_MARGIN:.2f} um for this hybrid decomposition to "
    f"have a valid interior interface at all, regardless of how small "
    f"BODY_STUB is made. coupling_length=0.5 (device 1921/1922) is EXPECTED "
    f"to hit this -- use directional_coupler_soi160_fdtd_dc_test.py instead "
    f"(full single-domain FDTD has no such constraint; already run "
    f"successfully at this length, job 11785063)."
)
BODY_STUB = min(0.25, COUPLING_LENGTH * 0.15, _body_stub_cap)

N_MODES = 2   # fundamental-mode even/odd supermode pair -- same measured
              # simplification as directional_coupler_soi160_hybrid.py
              # (its own docstring point 2); unchanged by coupling_length.

CALIBRATE_BODY = True

# ---------------------------------------------------------------------------
# 1. Geometry -- same real duty/dc.gds dimensions, coupling_length varies.
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
print(f"Directional coupler (SOI_160, duty/dc.gds dims, coupling_length={COUPLING_LENGTH:.3f}) "
      f"total length: {dc.total_length:.2f} um (coupling pitch {dc.coupling_pitch:.3f} um, "
      f"far pitch {dc.waveguide_offset:.2f} um)")

# ---------------------------------------------------------------------------
# 2. Configure sub-solvers -- one per section.
# ---------------------------------------------------------------------------
_fdtd_common = dict(
    res=FDTD_RESOLUTION,
    sim_width=7.0,
    sim_height=1.5,
    dpml=_DPML,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    half_domain_aperture=True,
    material_stack=SOI_160,
    decay_tol=1e-11,
)

fdtd_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=7.0, oxide_height=1.5, res=EME_RESOLUTION, material_stack=SOI_160,
    mode_solver="meep",
)

# ---------------------------------------------------------------------------
# 3. Build the section boundaries -- FDTD bend | EME coupled body | FDTD bend.
# ---------------------------------------------------------------------------
half_body = dc.coupling_length / 2
half_total = dc.total_length / 2

boundaries = [
    -half_total,
    -half_body + BODY_STUB,
    +half_body - BODY_STUB,
    +half_total,
]
interface_modes = [N_MODES, N_MODES]

body_length = boundaries[2] - boundaries[1]
assert body_length > 0, (
    f"BODY_STUB ({BODY_STUB:.3f} um) too large for COUPLING_LENGTH "
    f"({COUPLING_LENGTH:.3f} um) -- the EME body length would be negative. "
    f"coupling_length=0 (device 1931/1932) is EXPECTED to hit this -- see "
    f"module docstring, that case needs a different construction entirely."
)
print(f"Coupling-length budget: 2 x {BODY_STUB:.3f} um inside the FDTD "
      f"sections + EME body {body_length:.3f} um = {COUPLING_LENGTH:.3f} um total")

_iface_mon = dict(monitor_center=0.0, monitor_size=fdtd_solver.sim_width)

# ---------------------------------------------------------------------------
# 4. Geometry plot.
# ---------------------------------------------------------------------------
hybrid = HybridSolver(
    solvers=[fdtd_solver, body_solver, fdtd_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
)
plot_hybrid_geometry(
    dc, solver=hybrid,
    in_modes=1, out_modes=1,
    title=f"Directional Coupler (SOI_160, coupling_length={COUPLING_LENGTH:.3f} um) — Hybrid Solver Sections",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True,
    sub_region_show_eme_detail=True,
    sub_region_configs=[
        ("in", {0: {"input_idx": 0, **_iface_mon}}, {2: 0}),
        ("isolated", {0: {"input_idx": 1, **_iface_mon}}, {2: 0}),
    ],
)
print(f"Geometry plot saved -> {output_dir}/geometry_hybrid.png")
print(f"Sub-region plots saved -> {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run -- 2 real FDTD solves + 1 EME mode-solve per wavelength. NOT a
#    cache hit against the coupling_length=3.0 run's cache entries -- see
#    module docstring's CACHE NOTE.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nWavelength sweep: {wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

M_in0 = fdtd_solver.solve_section(
    dc, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
    input_idx=0, save_dir=save_dir, **_iface_mon,
)
M_in1 = fdtd_solver.solve_section(
    dc, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
    input_idx=1, save_dir=save_dir, **_iface_mon,
)

M_out_through = M_in0.transpose(0, 2, 1)
M_out_cross = M_in1.transpose(0, 2, 1)

_body_neff = None
if CALIBRATE_BODY:
    _beta = fdtd_solver.section_mode_beta(
        dc, wavelengths, boundaries[1], boundaries[2], N_MODES,
        save_dir=save_dir, **_iface_mon,
    )
    _body_neff = _beta * wavelengths[:, None] / (2 * np.pi)
    _i = int(np.argmin(np.abs(wavelengths - 1.55)))
    _seed_half = min(1.0, half_body - 2.0 * fdtd_solver.dpml)
    _mpb = fdtd_solver.section_mode_neff(
        dc, wavelengths[_i:_i + 1], -_seed_half, +_seed_half, N_MODES, **_iface_mon)[0]
    print(f"Body beta calibration at lambda = {wavelengths[_i]:.3f} um "
          f"(even/odd supermodes):")
    print(f"    MPB mode index  n_eff {_mpb[0]:.6f} {_mpb[1]:.6f}   "
          f"delta_n_eff {_mpb[0] - _mpb[1]:.6f}")
    print(f"    FDTD own phase  n_eff {_body_neff[_i, 0]:.6f} {_body_neff[_i, 1]:.6f}   "
          f"delta_n_eff {_body_neff[_i, 0] - _body_neff[_i, 1]:.6f}   <- body propagates on this")

M_body = body_solver.solve_section_length_sweep(
    dc, wavelengths, boundaries[1], boundaries[2], np.array([body_length]),
    n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
    neff_override=_body_neff,
)

tmp = np.einsum("wij,wjk->wik", M_body[0], M_in0)
full_through = np.einsum("wij,wjk->wik", M_out_through, tmp)
full_cross = np.einsum("wij,wjk->wik", M_out_cross, tmp)

coeff_through = full_through[:, 0, 0]
coeff_cross = full_cross[:, 0, 0]
T_through = np.abs(coeff_through) ** 2
T_cross = np.abs(coeff_cross) ** 2
T_total = T_through + T_cross

wl_c = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_c))

print(f"\nResults at coupling_length = {COUPLING_LENGTH:.3f} um, lambda = {wl_c:.3f} um:")
print(f"  Through output T = {T_through[idx]:.4f}  IL = {10*np.log10(max(T_through[idx],1e-20)):.2f} dB")
print(f"  Cross output   T = {T_cross[idx]:.4f}  IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total          T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  Coupling ratio   = {T_cross[idx]/max(T_total[idx],1e-20):.3f}  "
      f"(0 = fully through, 1 = fully cross)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat = SMatrix(wavelengths=wavelengths)
smat.set(coeff_through, port="through")
smat.set(coeff_cross, port="cross")
smat.check_passive()
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix saved -> {output_dir}/smatrix_hybrid.pkl")

np.savez(
    os.path.join(output_dir, "coupling_result.npz"),
    wavelengths=wavelengths, T_through=T_through, T_cross=T_cross, T_total=T_total,
    coupling_length=COUPLING_LENGTH,
    fdtd_resolution=FDTD_RESOLUTION, eme_resolution=EME_RESOLUTION,
)
np.savez(
    os.path.join(output_dir, "coupling_pieces.npz"),
    wavelengths=wavelengths, M_in0=M_in0, M_in1=M_in1, beta=_beta,
    body_length=body_length, coupling_length=COUPLING_LENGTH,
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
ax.set_title("Insertion Loss vs Wavelength")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle(f"Directional Coupler (SOI_160, coupling_length={COUPLING_LENGTH:.3f} um) — FDTD/EME Hybrid", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved -> {output_dir}/hybrid_results.png")
