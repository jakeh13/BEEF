"""Tutorial: the SAME directional coupler method as directional_coupler_hybrid.py
(FDTD S-bends + EME coupled-pair body), re-derived on 160 nm SOI at a REAL
fabricated coupler's own dimensions --
validation_markdown/high_contrast_grating_fabry_perot_validation.md, step 2.

This is the second physical half (after the HCG mirror itself) of the
HCG Fabry-Perot resonator that doc plans. The HCG mirror is `SOI_160`;
directional_coupler_hybrid.py is `SOI_220` -- the two could not be cascaded
into one cavity as shipped (see that doc's step 0/A1). This script redesigns
the coupler on `SOI_160` instead of the HCG's platform, at the real
dimensions of `duty/dc.gds` (read once with `gdstk`, NOT imported as
`GdsGeometry` -- see the doc's step 0 table for exactly how each number
below was read and cross-checked against DirectionalCoupler's own
constructor semantics):

    waveguide_width=0.5, gap=0.25, coupling_length=3.0,
    bend_length=10.0, waveguide_offset=4.0, bend_type="euler" (assumed,
    not re-derived from the GDS's own curve family -- a modeling choice).

Two differences from directional_coupler_hybrid.py, both load-bearing:

1. coupling_length is FIXED at the GDS's own 3.0 um, not swept for a 50:50
   split -- this is a real fabricated device's own dimension, not a free
   design parameter (step 0's decision: use the real GDS numbers directly).
   So there is no coupling_length sweep, no L50 search, no
   PLOT_COUPLING_LENGTH pin -- just ONE length, measured once.

2. N_MODES=2, not 5. The real 0.5 um / 160 nm strip carries 2 guided modes
   (fundamental + a genuinely guided 2nd-order mode, both with real margin
   above the SiO2 cladding index -- step 0's own MPB check), so the coupled
   cross-section itself supports 4 supermodes, not 2. But a follow-up
   MEASUREMENT (a real FDTDSolver.solve_section() run of this exact S-bend,
   res=25, decay_tol=1e-10, mode 1 excited from the far access port,
   decomposed onto all 4 supermodes at the coupled-region interface) found
   less than 0.01% of the launched power reaches the 2nd-order-mode pair
   (power fractions: 0.4965/0.4953 in the fundamental pair, 5.5e-5/2.9e-5 in
   the 2nd-order pair). The large index gap between the two modes
   (n_eff ~2.15 vs ~1.50 at 1.55 um) suppresses bend-induced conversion at
   this bend_length, and that is now measured, not argued. N_MODES=2 here
   is therefore a MEASURED simplification of a genuinely multimode
   cross-section, not a redesigned single-mode waveguide and not an
   unverified assumption either.

FDTD_RESOLUTION=40 / EME_RESOLUTION=80 are ASSUMED from
directional_coupler_hybrid.py's own SOI_220 convergence study
(dc_validation.md item 11), NOT re-measured on this platform -- the user's
explicit instruction for this first pass. This is a real, stated
approximation: nothing in dc_validation.md's own convergence numbers says a
different core thickness (160 vs 220 nm) or different device dimensions
(this coupler's own width/gap/bend_length/offset all differ from the
SOI_220 script's) preserves the same converged resolution. If this
script's own result looks physically implausible (T_through+T_cross far
from what check_passive() expects, or wildly different from a bare
supermode-beat estimate), a real convergence sweep on THIS platform is the
first thing to distrust, not the coupling physics.

CALIBRATE_BODY=True for the same reason directional_coupler_hybrid.py uses
it and calls it non-optional (dc_validation.md item 6): meep's own
time-stepped phase and MPB's mode index disagree by a few percent on Δβ,
and Δβ is what sets the coupling ratio.

WAVELENGTH WINDOW (added 2026-08-07, high_contrast_grating_fabry_perot_
validation.md step 5's follow-up): the first run of this script (job
11716260, N_WL=11 across the full 1.50-1.60 um C-band) fed into the
Fabry-Perot closed-form check and turned up a real resonance near
1.570 um -- but 11 points at 10 nm spacing cannot resolve it. A finesse
estimate from that SAME real run's own numbers (|rho|^2=0.9714,
|t|^2=0.8989, round-trip power |rho^2 t^2|=0.8732) gives
F = pi*sqrt(0.8732)/(1-0.8732) ~ 23, and with Step 1's ~9.5 nm FSR
estimate that puts FWHM ~ 0.41 nm -- some 24x narrower than the original
10 nm grid, not just "a bit coarse." DC160_WL_MIN/DC160_WL_MAX/DC160_N_WL
override WL_MIN/WL_MAX/N_WL below (unset keeps the original full-band
production behavior exactly) so a narrow, dense window can bracket the
found resonance without re-paying the cost of densifying the whole
100 nm band: the same single broadband FDTD pulse captures every
requested wavelength's DFT coefficient "for free" regardless of how many
are requested, so this is not N_WL separate solves, but a real,
independent SLURM run all the same (a different wavelength array is a
different cache fingerprint, so nothing from job 11716260's cache
carries over except the geometry- and length-only sub-solves that don't
depend on wavelengths, e.g. section_mode_beta()'s calibration-length
FDTD runs, still valid unchanged). A non-default window redirects
output_dir to its own subfolder (see below) so it can never overwrite the
original full-band production result job 11716260 already produced.

Run:
    python /path/to/examples/directional_coupler/python/directional_coupler_soi160_hybrid.py
    DC160_WL_MIN=1.555 DC160_WL_MAX=1.585 DC160_N_WL=601 python .../directional_coupler_soi160_hybrid.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import (
    DirectionalCoupler, EMESolver, FDTDSolver, HybridSolver, SMatrix, SOI_160, plot_hybrid_geometry,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "directional_coupler", "directional_coupler_soi160_hybrid"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
# ASSUMED from directional_coupler_hybrid.py's own SOI_220 convergence study
# (dc_validation.md item 11) -- NOT independently measured on SOI_160, on
# this device's own dimensions, per the user's explicit instruction for this
# first pass. See module docstring.
FDTD_RESOLUTION = 40.0   # pixels/um -- ASSUMED, not re-measured here
EME_RESOLUTION = 80.0    # pixels/um -- ASSUMED, not re-measured here
WL_MIN = 1.50       # C-band, matches the HCG mirror's own band
WL_MAX = 1.60
N_WL = 11

# DC160_RES_SCALE (check 3's resolution-convergence sweep) scales BOTH
# engines together and redirects output -- see the DC160_BODY_STUB block
# below (same env-var override family) for the mutual-exclusivity guard and
# the single-wavelength cost reasoning.
_res_scale_env_early = os.environ.get("DC160_RES_SCALE")
if _res_scale_env_early is not None:
    _res_scale = float(_res_scale_env_early)
    FDTD_RESOLUTION *= _res_scale
    EME_RESOLUTION *= _res_scale
    output_dir = os.path.join(output_dir, "resolution_sweep", f"scale{_res_scale:g}")
    print(f"DC160 resolution sweep override: scale={_res_scale} -> "
          f"FDTD_RESOLUTION={FDTD_RESOLUTION} EME_RESOLUTION={EME_RESOLUTION} "
          f"-> output_dir={output_dir}")
    os.makedirs(output_dir, exist_ok=True)

# DC160_WL_MIN/DC160_WL_MAX/DC160_N_WL override the full-band production
# window above -- see module docstring's "WAVELENGTH WINDOW" section. Unset
# (the default) reproduces job 11716260's original full-band behavior
# exactly, output_dir unchanged. Any override redirects output_dir to its
# own subfolder -- "the guard is the env var, not the value" (same pattern
# ring_resonator_hybrid_cascaded.py's RING_FDTD_RES/RING_THREE_D use), so a
# dense zoom run can never silently overwrite the full-band production
# smatrix_hybrid.pkl this script already produced once.
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

# DC160_RES_SCALE (check 3, Step 2 item 2's resolution-convergence gap) and
# DC160_BODY_STUB (check 1's cut-point/stub-invariance test, validation
# doc's step 6) are two SEPARATE, mutually exclusive sweep axes -- each
# changes what the FDTD cache actually re-solves, so combining them would
# make a single run's cost/cache-key story ambiguous. Both force N_WL=1 at
# lambda=1.55 um (unless the wavelength window above was ALSO explicitly
# set) because these sweeps test a NUMERICAL PROPERTY (does resolution or
# cut-point choice change the answer), not the spectrum -- full N_WL=11
# decomposition cost would be paid for nothing new.
_res_scale_env = os.environ.get("DC160_RES_SCALE")
_stub_env = os.environ.get("DC160_BODY_STUB")
if _res_scale_env is not None and _stub_env is not None:
    raise SystemExit(
        "DC160_RES_SCALE and DC160_BODY_STUB are mutually exclusive -- each "
        "is its own sweep axis with its own cache/output-path story; set "
        "only one per run."
    )
if (_res_scale_env is not None or _stub_env is not None) and not any(v is not None for v in _win_env):
    WL_MIN = WL_MAX = 1.55
    N_WL = 1
    print("DC160 resolution/stub sweep: forcing single-wavelength (1.55 um) "
          "mode -- this tests a numerical property, not the spectrum.")

os.makedirs(output_dir, exist_ok=True)

# The real fabricated coupler's own straight coupled-pair length -- FIXED,
# not swept for a 50:50 split. See module docstring point 1.
COUPLING_LENGTH = 3.0   # um -- duty/dc.gds, step 0

# um of the straight coupling region folded into each FDTD section, past the
# S-bend/coupling-region discontinuity -- same convention (and same value)
# as directional_coupler_hybrid.py's BODY_STUB. 2*BODY_STUB=0.5 um leaves a
# body_length of 2.5 um out of the 3.0 um total -- comfortably positive.
BODY_STUB = 0.25   # um

# DC160_BODY_STUB overrides BODY_STUB for check 1's cut-point sweep --
# validation_markdown/high_contrast_grating_fabry_perot_validation.md
# step 6. Must stay < COUPLING_LENGTH/2 = 1.5 (checked below, same guard
# the default value already passes).
if _stub_env is not None:
    BODY_STUB = float(_stub_env)
    output_dir = os.path.join(output_dir, "stub_sweep", f"stub{BODY_STUB:.3f}")
    print(f"DC160 body-stub sweep override: BODY_STUB={BODY_STUB} -> output_dir={output_dir}")
    os.makedirs(output_dir, exist_ok=True)

N_MODES = 2   # fundamental-mode even/odd supermode pair ONLY -- see module
              # docstring point 2 for the measurement justifying this.

CALIBRATE_BODY = True   # see module docstring -- not optional, same
                        # reasoning as directional_coupler_hybrid.py's own flag.

# ---------------------------------------------------------------------------
# 1. Geometry -- the real duty/dc.gds dimensions, on SOI_160.
# ---------------------------------------------------------------------------
dc = DirectionalCoupler(
    waveguide_width=0.5,     # duty/dc.gds, step 0 -- single-mode strip width
    gap=0.25,                # duty/dc.gds, step 0 -- coupling-region edge-to-edge gap
    coupling_length=COUPLING_LENGTH,   # duty/dc.gds, step 0 -- FIXED, not swept
    bend_length=10.0,        # duty/dc.gds, step 0
    waveguide_offset=4.0,    # duty/dc.gds, step 0 -- far-end centre-to-centre port pitch
    bend_type="euler",       # ASSUMED -- not re-derived from the GDS's own curve family
    waveguide_height=0.16,   # SOI_160
)
print(f"Directional coupler (SOI_160, duty/dc.gds dimensions) total length: "
      f"{dc.total_length:.2f} um (coupling pitch {dc.coupling_pitch:.3f} um, "
      f"far pitch {dc.waveguide_offset:.2f} um)")

# ---------------------------------------------------------------------------
# 2. Configure sub-solvers -- one per section.
# ---------------------------------------------------------------------------
# waveguide_offset=4.0 here vs directional_coupler_hybrid.py's 3.0 --
# sim_width bumped from 6.0 to 7.0 to preserve that script's own PML
# clearance ratio (half_domain_aperture, so each arm needs waveguide_offset/2
# + clearance inside the non-PML half-domain).
_fdtd_common = dict(
    res=FDTD_RESOLUTION,
    sim_width=7.0,
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    half_domain_aperture=True,
    material_stack=SOI_160,
    decay_tol=1e-11,
)

fdtd_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholder length; solve_section_length_sweep() below sets the real one
    oxide_width=7.0, oxide_height=1.5, res=EME_RESOLUTION, material_stack=SOI_160,
    mode_solver="meep",
)

# ---------------------------------------------------------------------------
# 3. Build the section boundaries -- FDTD bend | EME coupled body | FDTD bend.
# ---------------------------------------------------------------------------
half_body = dc.coupling_length / 2
half_total = dc.total_length / 2

boundaries = [
    -half_total,               # physical input port (access waveguide end)
    -half_body + BODY_STUB,    # input section <-> body section interface
    +half_body - BODY_STUB,    # body section <-> output section interface
    +half_total,               # physical output port (access waveguide end)
]
interface_modes = [N_MODES, N_MODES]

body_length = boundaries[2] - boundaries[1]
assert body_length > 0, (
    f"BODY_STUB ({BODY_STUB:.2f} um) too large for COUPLING_LENGTH "
    f"({COUPLING_LENGTH:.2f} um) -- the EME body length would be negative"
)
print(f"Coupling-length budget: 2 x {BODY_STUB:.2f} um inside the FDTD "
      f"sections + EME body {body_length:.2f} um "
      f"= {COUPLING_LENGTH:.2f} um total")

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
    title=f"Directional Coupler (SOI_160, duty/dc.gds dims) — Hybrid Solver Sections",
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
# 5. Run -- 2 FDTD solves (input section, once per input arm) + 1 EME
#    mode-solve per wavelength, at the ONE real coupling_length. The output
#    section is never separately simulated: M_in0/M_in1 double as the
#    output-side transfer matrices once transposed (mirror-symmetry reuse,
#    identical reasoning to directional_coupler_hybrid.py).
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nWavelength sweep: {wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

M_in0 = fdtd_solver.solve_section(
    dc, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
    input_idx=0, save_dir=save_dir, **_iface_mon,
)  # [n_wl, N_MODES, 1] -- "in" arm -> body modes
M_in1 = fdtd_solver.solve_section(
    dc, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
    input_idx=1, save_dir=save_dir, **_iface_mon,
)  # [n_wl, N_MODES, 1] -- "isolated" arm -> body modes

M_out_through = M_in0.transpose(0, 2, 1)  # [n_wl, 1, N_MODES] -- mirror-symmetry reuse
M_out_cross = M_in1.transpose(0, 2, 1)    # [n_wl, 1, N_MODES]

_body_neff = None
if CALIBRATE_BODY:
    _beta = fdtd_solver.section_mode_beta(
        dc, wavelengths, boundaries[1], boundaries[2], N_MODES,
        save_dir=save_dir, **_iface_mon,
    )  # [n_wl, N_MODES] -- rad/um, fitted from the FDTD's own phase vs length
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
)  # [1, n_wl, N_MODES, N_MODES]

tmp = np.einsum("wij,wjk->wik", M_body[0], M_in0)          # [n_wl, N_MODES, 1]
full_through = np.einsum("wij,wjk->wik", M_out_through, tmp)   # [n_wl, 1, 1]
full_cross = np.einsum("wij,wjk->wik", M_out_cross, tmp)        # [n_wl, 1, 1]

coeff_through = full_through[:, 0, 0]
coeff_cross = full_cross[:, 0, 0]
T_through = np.abs(coeff_through) ** 2
T_cross = np.abs(coeff_cross) ** 2
T_total = T_through + T_cross

wl_c = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_c))

print(f"\nResults at coupling_length = {COUPLING_LENGTH:.2f} um (fixed, real GDS "
      f"dimension), lambda = {wl_c:.3f} um:")
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

# The intermediate transfer-matrix pieces (NOT just the final coeff_through/
# coeff_cross) -- needed by directional_coupler_soi160_hybrid_dense.py to
# reconstruct a dense spectrum via interpolation rather than a brute-force
# dense FDTD re-solve. M_in0/M_in1 (arm-excitation -> body-mode overlaps) and
# _beta (the FDTD-calibrated body propagation constant) each vary slowly and
# safely with wavelength on their own; only the FINAL combined coefficient's
# phase varies rapidly (the S-bend path is traversed twice by the
# mirror-symmetry reuse above, doubling its phase rate) -- see that script's
# own module docstring for the real, measured justification.
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
ax.set_title(f"Transmission vs Wavelength (coupling_length={COUPLING_LENGTH:.2f} um, fixed)")
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

plt.suptitle("Directional Coupler (SOI_160, duty/dc.gds dims) — FDTD/EME Hybrid", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved -> {output_dir}/hybrid_results.png")
