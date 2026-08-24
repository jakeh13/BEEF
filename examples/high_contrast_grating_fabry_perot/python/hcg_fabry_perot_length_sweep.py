"""Tutorial: 2-D design sweep of the HCG Fabry-Perot resonator over the
directional coupler's COUPLING LENGTH and the connecting waveguide's ARM
LENGTH, reporting extinction ratio, FSR and FWHM as three heatmaps.

This is hcg_fabry_perot_hybrid.py's device, swept. That script solves ONE
operating point (COUPLING_LENGTH=3.0 um, L_ARM=53.08 um -- the real
fabricated structure). This one asks the design question instead: how do the
cavity's three headline numbers move across a (coupling_length, arm_length)
plane?

WHY THIS COSTS (ALMOST) NO SIMULATION -- the whole point of the script.

Read hcg_fabry_perot_hybrid.py's closed form (its section 6):

    rho   = r * phase**2
    D     = 1 - rho**2 * t**2
    T     = |t * (1 + rho**2 * (k**2 - t**2)) / D|**2

and then note where the two swept variables actually enter:

  * L_ARM enters ONLY through `phase = exp(2j*pi*n_eff_arm/wl * L_ARM)`.
    Nothing simulated depends on it. n_eff_arm is measured ONCE on a bare
    6 um strip and then acts analytically over any arm length.

  * coupling_length enters ONLY through `t` and `k`, and within those ONLY
    through the EME body's diagonal phase,
        M_body[m, m] = exp(1j * beta_body[m] * body_length),
        body_length  = coupling_length - 2*BODY_STUB.
    The two FDTD S-bend excitations (M_in0/M_in1) and the body's beta
    calibration are properties of the BEND and the CROSS-SECTION, not of how
    much straight coupler sits between them.

So the entire 2-D plane is analytic post-processing over FOUR cached
simulations. This is the same construction
directional_coupler_soi160_hybrid_length_sweep.py uses for the coupler alone
(reference geometry at L_REF, body swept analytically), extended one axis.

THE FOUR PIECES, AND WHERE THEIR CACHE ENTRIES COME FROM
(all four are `save_dir=example_sim_cache/` lookups; set HCGFP_CACHE_ONLY=1
to confirm coverage with ZERO computation before spending a real allocation):

  1. HCG mirror r(lambda)   -- FDTDSolver.solve_reflection(), res RESOLUTION.
     Byte-identical to the call high_contrast_grating_fdtd.py makes, so it
     hits that script's entry at the matching res/band.
  2/3. Coupler S-bends M_in0/M_in1 and the body beta calibration --
     FDTDSolver.solve_section() x2 + section_mode_beta_details(), at the
     L_REF geometry. Byte-identical to directional_coupler_soi160_hybrid.py
     at the matching band.
  4. Arm n_eff        -- section_mode_beta_details() on a bare WG_WIDTH strip
     at ARM_RESOLUTION. Shared with hcg_fabry_perot_hybrid.py.

DEFAULTS ARE PINNED TO THE 2026-08-21 C-BAND MEASUREMENT RE-RUNS so that
pieces 1-3 are cache HITS out of the box: band 1.530-1.565 um / 36 pts (jobs
12093321 and 12093330), mirror res 80 (job 12093330), coupler res 40
(job 12093321, unchanged -- still that script's own ASSUMED value).

RESOLUTION IS NOT UNIFORM ACROSS THIS DEVICE, and that is inherited, not
introduced here. The mirror runs at 80, the arm at 80, the coupler's S-bends
at 40 and its EME body at 80. hcg_fabry_perot_hybrid.py's own docstring
already records why there is no single number to inherit and what the cost
is; raising the coupler to 80 as well is a real, unrun convergence job, not
a constant edit -- see "KNOWN OPEN GAPS" below.

WHAT THE THREE METRICS MEAN HERE

Each (coupling_length, arm_length) cell evaluates the closed form on a dense
wavelength grid, finds the resonance dips, and reduces them to one scalar:

  extinction ratio -- 10*log10(baseline / T_min) per dip, MEDIAN over dips.
      "How deep is a resonance", in dB. Baseline is the band median of T,
      the same reference the dip finder itself uses.
  FSR              -- mean spacing between adjacent dip minima, in nm.
      Needs >= 2 dips in the band; see the band-coverage limit below.
  FWHM             -- full width at half depth per dip, MEDIAN over dips, nm.
      Half-power crossings are interpolated between dense-grid samples
      (`_cross`), so this does not staircase on the grid step.

ER RIPPLE -- a real property of the metric, measured 2026-08-22, not noise
to be smoothed away silently. ER is a MEDIAN over whichever resonances fall
inside the band, and as coupling_length shifts the comb, dips enter and leave
that fixed window: over a 400-point Lc scan the in-band resonance count moved
between 16.2 and 23.6 and changed 392 times. The resulting jitter on ER is
std 0.223 dB, peak-to-peak 2.561 dB about a smoothed envelope.

CONSEQUENCE FOR READING THE ER MAP: individual cell-to-cell maxima are NOT
resolvable design points. On that same scan the apparent "global optimum"
(32.30 dB at Lc = 26.60) sat only +1.39 dB above the envelope -- inside the
ripple amplitude. Read the ENVELOPE, not the argmax. The physical result there
was a broad plateau peaking ~31.3 dB near Lc ~ 21.7 um and staying within
1 dB over Lc in [17.8, 36.6] um.

DENSE-GRID SENSITIVITY, measured on this device rather than assumed: going
2001 -> 8001 -> 20001 dense points moved ER 8.35 -> 10.90 -> 10.88 dB while
FSR held at 2.2455 nm throughout and FWHM converged by 8001 (0.1198 ->
0.1010 -> 0.1010 nm). ER is the sensitive metric because it samples the dip
BOTTOM; a coarse grid under-reports depth silently. The default 20001 is
converged for this cavity and the script prints a samples-per-resonance
adequacy check and warns below 20.

MEDIAN, not mean, for ER and FWHM: near a coupling length where the cavity
crosses critical coupling, one dip in the band can go anomalously deep/narrow
and a mean would let that single cell-internal outlier speak for the cell.

BAND-COVERAGE LIMIT, the one real artifact to know about. FSR shrinks as
1/L_ARM, so SHORT arms put fewer resonances in a fixed band. At the default
35 nm band, two dips need FSR <~ 17 nm, i.e. L_ARM >~ 7 um -- below that the
FSR cell is NaN by arithmetic, not by physics. The dip COUNT per cell is
saved alongside the metrics (`n_dips` in the .npz, and its own diagnostic
panel) precisely so a NaN is explainable rather than mysterious. Widen the
band, or raise LARM_MIN, rather than reading a hole as a result.

KNOWN OPEN GAPS -- inherited from hcg_fabry_perot_hybrid.py, not fixed here:
  - The coupler's FDTD_RESOLUTION=40 is ASSUMED, never converged for this
    device. It is the one piece NOT raised to 80 in the defaults.
  - The symmetric-cavity assumption: the real structure's two arms differ by
    3.0 um (54.583 / 51.583). This sweep varies ONE shared arm length, so it
    cannot see that asymmetry -- the real device is a point on this plane
    only under the same averaging Step 1 already applies.
  - Short coupling lengths reuse the L_REF S-bend solves. RESOLVED
    2026-08-21, and not in the direction expected: this was flagged as likely
    to degrade as coupling_length -> 2*BODY_STUB (the two bends end up ~0.5 um
    apart), but the direct check against independent single-domain FDTD
    (jobs 12093324 vs 12093328_0) found the 0.5 um agreement as good as the
    3.0 um agreement -- coupling ratio matching to 0.0000, cross to 1e-4. The
    residual is a small T_total deficit that is LENGTH-INDEPENDENT (+0.00170
    at 0.5 um vs +0.00145 at 3.0 um, with the bends 6x closer), i.e. a fixed
    N_MODES=2 interface-truncation cost, not a proximity failure. The left
    edge of the coupling-length axis is therefore usable, not provisional.

Run:
    python .../hcg_fabry_perot_length_sweep.py
    HCGFP_CACHE_ONLY=1 python .../hcg_fabry_perot_length_sweep.py   # coverage probe, no compute
    HCGFP_LC_N=60 HCGFP_LARM_N=60 python .../hcg_fabry_perot_length_sweep.py
"""

import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import (
    DirectionalCoupler, EMESolver, FDTDSolver, GdsGeometry, SOI_160,
    load_gds_polygons,
)
from beef.cache import is_master
from beef.geometry import UniformCrossSection
from beef.spectral import interp_complex, unwrap_is_safe

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
_GDS_PATH = os.path.abspath(os.path.join(_HERE, "..", "..", "gds_files", "HCG_seeded_TO.gds"))
save_dir = os.path.join(_ROOT, "example_sim_cache")

MASTER = is_master()

# CACHE_ONLY makes every solver call a pure lookup that RAISES rather than
# simulates on a miss. Run this first, on a login node, before submitting: it
# answers "is this sweep free?" in seconds and costs nothing if the answer is
# no. Same knob, same spelling, as hcg_fabry_perot_hybrid.py's.
CACHE_ONLY = bool(int(os.environ.get("HCGFP_CACHE_ONLY", "0")))

# ---------------------------------------------------------------------------
# Band -- defaults pinned to the 2026-08-21 C-band measurement re-runs so
# pieces 1-3 are cache hits. Overridable as a set (all or none), same
# all-or-nothing guard as the DC160_*/HCG_* families elsewhere in this repo.
# ---------------------------------------------------------------------------
WL_MIN, WL_MAX, N_WL = 1.530, 1.565, 36
_win_env = (os.environ.get("HCGFP_WL_MIN"), os.environ.get("HCGFP_WL_MAX"),
            os.environ.get("HCGFP_N_WL"))
if any(v is not None for v in _win_env):
    if not all(v is not None for v in _win_env):
        raise SystemExit(
            "HCGFP_WL_MIN/HCGFP_WL_MAX/HCGFP_N_WL must be set together (or "
            "not at all) -- a partial override would silently mix a custom "
            "window with the defaults above."
        )
    WL_MIN, WL_MAX, N_WL = float(_win_env[0]), float(_win_env[1]), int(_win_env[2])

# The mirror's own grid. hcg_fabry_perot_hybrid.py runs it finer (21 vs 11)
# and requires it to be a SUPERSET of the assembly grid. On the 2026-08-21
# re-run band the mirror and the coupler were run on the SAME 36-point grid,
# so this is the identity by default -- the superset assertion below still
# runs, because the day someone overrides one and not the other is exactly
# when a silent regrid would otherwise happen.
N_WL_HCG = int(os.environ.get("HCGFP_N_WL_HCG", N_WL))

RESOLUTION = float(os.environ.get("HCGFP_RES", 80.0))          # HCG mirror FDTD  -- job 12093330
FDTD_RESOLUTION = float(os.environ.get("HCGFP_DC_RES", 40.0))  # coupler S-bends  -- job 12093321
EME_RESOLUTION = 80.0                                          # coupler EME body
ARM_RESOLUTION = float(os.environ.get("HCGFP_ARM_RES", 80.0))  # connecting waveguide

N_MODES = 2          # coupler body: fundamental even/odd supermode pair only
BODY_STUB = 0.25     # um of the coupling length folded into each FDTD section
L_REF = 3.0          # um -- reference geometry the S-bend solves are done at.
                     # Deliberately the die's own real length, NOT the middle
                     # of the sweep: it makes those solves cache HITS against
                     # directional_coupler_soi160_hybrid.py. Same choice, same
                     # reasoning, as directional_coupler_soi160_hybrid_length_sweep.py.
WG_WIDTH = 0.5
ARM_CAL_SPAN = 6.0
MESH_WIDTH = 4.0
SIM_HEIGHT = 1.5
DPML = 0.5
CALIBRATE_BODY = True
CALIBRATE_ARM = bool(int(os.environ.get("HCGFP_CALIBRATE_ARM", "1")))

# ---------------------------------------------------------------------------
# Sweep grids.
#
# LC_MIN floor is 2*BODY_STUB: below it the EME body length goes negative and
# the three-section decomposition does not exist. Identical floor, identical
# reason, as directional_coupler_soi160_hybrid_length_sweep.py's SWEEP_MIN.
# ---------------------------------------------------------------------------
LC_MIN = float(os.environ.get("HCGFP_LC_MIN", 2 * BODY_STUB))
LC_MAX = float(os.environ.get("HCGFP_LC_MAX", 30.0))
LC_N = int(os.environ.get("HCGFP_LC_N", 40))
LARM_MIN = float(os.environ.get("HCGFP_LARM_MIN", 10.0))
LARM_MAX = float(os.environ.get("HCGFP_LARM_MAX", 100.0))
LARM_N = int(os.environ.get("HCGFP_LARM_N", 40))

if LC_MIN < 2 * BODY_STUB - 1e-12:
    raise SystemExit(
        f"HCGFP_LC_MIN={LC_MIN} is below 2*BODY_STUB={2*BODY_STUB} -- the EME "
        f"body length would be negative and the [FDTD bend | EME body | FDTD "
        f"bend] decomposition would not exist at that length. A true "
        f"coupling_length -> 0 device needs single-domain FDTD instead (see "
        f"directional_coupler_soi160_fdtd_dc_test.py)."
    )

# Dense grid for the resonance metrics. Sized from the NARROWEST resonance
# the sweep actually produces, not from the fabricated device's:
# measured 2026-08-21 over the default 40x40 plane, the fabricated cell's
# FWHM is 0.101 nm but the high-finesse corner reaches 0.035 nm. At the
# default 35 nm band, 40001 points is 0.000875 nm -- 40 samples across even
# that narrowest cell, clearing the adequacy guard below with headroom.
# (20001 would give exactly 20, sitting on the threshold.) This costs no
# simulation, only the analytic loop: a 40x40 plane is ~9 s.
N_WL_DENSE = int(os.environ.get("HCGFP_N_DENSE", 40001))

# Dip-finder threshold: a local minimum counts as a resonance only if it
# falls below this fraction of the band baseline.
DIP_THRESHOLD = 0.7

REAL_LC, REAL_LARM = 3.0, 53.08   # the fabricated structure, marked on each map

output_dir = os.path.join(_ROOT, "example_output", "high_contrast_grating_fabry_perot",
                          "hcg_fabry_perot_length_sweep")
if not CALIBRATE_ARM:
    output_dir += "_uncalibrated_arm"
output_dir = os.path.join(
    output_dir, f"{WL_MIN:.4f}-{WL_MAX:.4f}_n{N_WL}",
    f"res{RESOLUTION:g}_dc{FDTD_RESOLUTION:g}_arm{ARM_RESOLUTION:g}")

# SWEEP-GRID REDIRECT (added 2026-08-22, after a real overwrite). The path
# above encodes the band and the three resolutions -- i.e. everything that
# changes what was SIMULATED -- but NOT the sweep axes, which change only
# what was evaluated. That was wrong: a fine 1-D-style rescan
# (LC_N=400/LARM_N=5, to resolve the ER envelope) silently overwrote the
# 40x40 production heatmaps at this same path, because both are "the same
# simulations".
#
# Same rule as the DC160_*/HCG_* window overrides: THE GUARD IS THE ENV VAR,
# NOT THE VALUE. Unset reproduces the documented production path exactly;
# setting any sweep-axis knob redirects to its own grid-tagged folder, so a
# diagnostic rescan can never clobber the production plane again.
_grid_env = [os.environ.get(v) for v in
             ("HCGFP_LC_MIN", "HCGFP_LC_MAX", "HCGFP_LC_N",
              "HCGFP_LARM_MIN", "HCGFP_LARM_MAX", "HCGFP_LARM_N")]
if any(v is not None for v in _grid_env):
    output_dir = os.path.join(
        output_dir, "grids",
        f"Lc{LC_MIN:g}-{LC_MAX:g}_n{LC_N}__La{LARM_MIN:g}-{LARM_MAX:g}_n{LARM_N}")
    print(f"Sweep-grid override -> output_dir={output_dir}")
if MASTER:
    os.makedirs(output_dir, exist_ok=True)

wl = np.linspace(WL_MIN, WL_MAX, N_WL)
wl_hcg = np.linspace(WL_MIN, WL_MAX, N_WL_HCG)
M = len(wl)

print("HCG Fabry-Perot resonator -- 2-D (coupling_length, arm_length) design sweep")
print(f"  band         : {WL_MIN:.4f}-{WL_MAX:.4f} um, {N_WL} pts "
      f"(mirror grid {N_WL_HCG} pts)")
print(f"  resolutions  : mirror {RESOLUTION:g}, coupler S-bend {FDTD_RESOLUTION:g}, "
      f"coupler body {EME_RESOLUTION:g}, arm {ARM_RESOLUTION:g}")
print(f"  coupling_len : {LC_N} pts, {LC_MIN:.3f}-{LC_MAX:.3f} um")
print(f"  arm_length   : {LARM_N} pts, {LARM_MIN:.3f}-{LARM_MAX:.3f} um")
print(f"  dense grid   : {N_WL_DENSE} pts ({(WL_MAX-WL_MIN)/(N_WL_DENSE-1)*1000:.5f} nm)")
print(f"  cache_only   : {CACHE_ONLY}")
print(f"  output       : {output_dir}")

# ===========================================================================
# 1. HCG mirror -- piece 1 of 4.
# ===========================================================================
print(f"\n--- 1. HCG mirror (FDTD, res {RESOLUTION:g}, 3-D, {N_WL_HCG} pts) ---")
_t0 = time.time()

CORE_LAYERS = [(2, 0)]
EDGE_CLEARANCE = 0.3
# region_lengths derived from the file itself, NOT hardcoded -- byte-identical
# to high_contrast_grating_fdtd.py's own construction, which is what makes the
# solve_reflection() call below a cache hit rather than a 3-hour re-solve.
_raw = load_gds_polygons(_GDS_PATH, layers=CORE_LAYERS)
_xs = [x for poly in _raw[CORE_LAYERS[0]] for x, _y in poly]
hcg_geom = GdsGeometry(
    gds_path=_GDS_PATH,
    core_layers=CORE_LAYERS,
    region_lengths=[max(_xs) - min(_xs)],
    input_width=WG_WIDTH,
    input_y=0.0,
)
hcg_solver = FDTDSolver(
    res=RESOLUTION, sim_width=7.5, sim_height=SIM_HEIGHT, dpml=DPML,
    source_margin=0.45, monitor_margin=0.8, three_d=True,
    num_input_modes=1, material_stack=SOI_160, decay_tol=1e-11,
)
smat_hcg = hcg_solver.solve_reflection(
    hcg_geom, wl_hcg, port="top", input_mode=1,
    edge_clearance=EDGE_CLEARANCE, extend_input_waveguide=True,
    save_dir=save_dir, cache_only=CACHE_ONLY,
)
r_hcg = smat_hcg.coeff(port="in", in_mode=1, out_mode=1)
R_hcg = np.abs(r_hcg) ** 2
print(f"HCG mirror: R max {R_hcg.max():.4f} at {wl_hcg[np.argmax(R_hcg)]:.4f} um, "
      f"min {R_hcg.min():.4f}  ({time.time()-_t0:.1f} s)")

# The mirror grid must CONTAIN the assembly grid -- never silently regridded.
_idx = [int(np.argmin(np.abs(wl_hcg - w))) for w in wl]
if not np.allclose(wl_hcg[_idx], wl, atol=1e-9):
    raise SystemExit(
        f"The mirror's {N_WL_HCG}-point grid is not a superset of the "
        f"{N_WL}-point assembly grid -- r(lambda) would have to be resampled "
        f"before the per-factor dense reconstruction below, which is exactly "
        f"what that scheme exists to avoid. Set HCGFP_N_WL_HCG to match."
    )
r = r_hcg[_idx]

# ===========================================================================
# 2/3. Directional coupler at the REFERENCE length -- pieces 2 and 3 of 4.
#
#      Solved once, at L_REF. The sweep varies only the EME body's analytic
#      diagonal phase downstream of here.
# ===========================================================================
print(f"\n--- 2. Coupler S-bends + body beta at L_REF = {L_REF:.2f} um "
      f"({N_WL} pts) ---")
_t0 = time.time()

dc = DirectionalCoupler(
    waveguide_width=WG_WIDTH, gap=0.25, coupling_length=L_REF,
    bend_length=10.0, waveguide_offset=4.0, bend_type="euler",
    waveguide_height=SOI_160.wg_height,
)
dc_fdtd = FDTDSolver(
    res=FDTD_RESOLUTION, sim_width=7.0, sim_height=SIM_HEIGHT, dpml=DPML,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    half_domain_aperture=True, material_stack=SOI_160, decay_tol=1e-11,
)
dc_eme = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=7.0, oxide_height=SIM_HEIGHT, res=EME_RESOLUTION,
    material_stack=SOI_160, mode_solver="meep",
)

half_body = dc.coupling_length / 2
half_total = dc.total_length / 2
boundaries = [-half_total, -half_body + BODY_STUB, +half_body - BODY_STUB, +half_total]
_iface_mon = dict(monitor_center=0.0, monitor_size=dc_fdtd.sim_width)

M_in0 = dc_fdtd.solve_section(
    dc, wl, boundaries[0], boundaries[1], n_modes_in=1, n_modes_out=N_MODES,
    stub_left=dc_fdtd.waveguide_length, stub_right=0.0, input_idx=0,
    save_dir=save_dir, cache_only=CACHE_ONLY, **_iface_mon,
)
M_in1 = dc_fdtd.solve_section(
    dc, wl, boundaries[0], boundaries[1], n_modes_in=1, n_modes_out=N_MODES,
    stub_left=dc_fdtd.waveguide_length, stub_right=0.0, input_idx=1,
    save_dir=save_dir, cache_only=CACHE_ONLY, **_iface_mon,
)
cal_body = dc_fdtd.section_mode_beta_details(
    dc, wl, boundaries[1], boundaries[2], N_MODES,
    save_dir=save_dir, cache_only=CACHE_ONLY, **_iface_mon,
)
beta_body = cal_body.beta                       # [n_wl, N_MODES] rad/um
_i = int(np.argmin(np.abs(wl - 0.5 * (WL_MIN + WL_MAX))))
_neff_body = beta_body * wl[:, None] / (2 * np.pi)
_mpb_body = cal_body.seed_beta * wl[:, None] / (2 * np.pi)
print(f"Body beta calibration at {wl[_i]:.4f} um (even/odd supermodes):")
print(f"    MPB mode index  n_eff {_mpb_body[_i,0]:.6f} {_mpb_body[_i,1]:.6f}   "
      f"delta {_mpb_body[_i,0]-_mpb_body[_i,1]:.6f}")
print(f"    FDTD own phase  n_eff {_neff_body[_i,0]:.6f} {_neff_body[_i,1]:.6f}   "
      f"delta {_neff_body[_i,0]-_neff_body[_i,1]:.6f}   <- body propagates on this")
print(f"Coupler pieces: {time.time()-_t0:.1f} s")

# ===========================================================================
# 4. Connecting waveguide -- piece 4 of 4.
# ===========================================================================
print(f"\n--- 3. Connecting waveguide ({WG_WIDTH:.2f} um strip, "
      f"res {ARM_RESOLUTION:g}) ---")
_t0 = time.time()

arm_geom = UniformCrossSection(
    [(0.0, WG_WIDTH)], length=ARM_CAL_SPAN, waveguide_height=SOI_160.wg_height,
)
arm_solver = FDTDSolver(
    res=ARM_RESOLUTION, sim_width=MESH_WIDTH, sim_height=SIM_HEIGHT, dpml=DPML,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    num_input_modes=1, num_output_modes=1, output_monitor_centers=[0.0],
    material_stack=SOI_160, decay_tol=1e-11,
)
if CALIBRATE_ARM:
    cal_arm = arm_solver.section_mode_beta_details(
        arm_geom, wl, -ARM_CAL_SPAN / 2.0, +ARM_CAL_SPAN / 2.0, 1,
        save_dir=save_dir, cache_only=CACHE_ONLY,
    )
    neff_arm = cal_arm.beta[:, 0] * wl / (2 * np.pi)
    neff_arm_mpb = cal_arm.seed_beta[:, 0] * wl / (2 * np.pi)
    print(f"Arm n_eff at {wl[_i]:.4f} um: MPB {neff_arm_mpb[_i]:.6f} -> "
          f"FDTD {neff_arm[_i]:.6f}  ({neff_arm[_i]/neff_arm_mpb[_i]-1.0:+.4%})")
    print("    <- the arms propagate on the FDTD value; every other block in "
          "this cavity is already referenced to meep's own beta.")
else:
    # Narrow window, save_dir, cache_only -- byte-identical to the parent
    # script's own CALIBRATE_ARM=0 branch, so the two share a cache entry.
    neff_arm = arm_solver.section_mode_neff(
        arm_geom, wl, -1.0, +1.0, 1, save_dir=save_dir, cache_only=CACHE_ONLY,
    )[:, 0]
    neff_arm_mpb = neff_arm
    print(f"HCGFP_CALIBRATE_ARM=0: arms on the MPB index "
          f"({neff_arm[0]:.6f} to {neff_arm[-1]:.6f}) -- diagnostic, own folder.")
print(f"Arm calibration: {time.time()-_t0:.1f} s")

# ===========================================================================
# 5. Per-factor dense reconstruction.
#
#    Identical scheme, identical justification, as hcg_fabry_perot_hybrid.py
#    section 8: each FACTOR is resampled where its own phase winding is slow
#    enough to be safe, and the fast-winding COMBINATION (t, k, and the
#    cavity response) is then rebuilt analytically on the dense grid rather
#    than ever being interpolated. Every resampled factor is checked with
#    unwrap_is_safe rather than assumed safe.
#
#    The 2026-08-21 re-run band helps here: the assembly grid is 1 nm, not
#    the 10 nm this scheme was originally validated against, so every
#    per-factor unwrap has ~10x the margin it was accepted with.
# ===========================================================================
print(f"\n--- 4. Dense reconstruction ({N_WL_DENSE} pts) ---")
wl_dense = np.linspace(WL_MIN, WL_MAX, N_WL_DENSE)


def _check_safe(c, label, margin=0.5):
    safe, worst = unwrap_is_safe(np.asarray(c), margin=margin)
    flag = "ok" if safe else "UNSAFE"
    print(f"    unwrap {label:<34} worst {worst:6.3f} rad/step  [{flag}]")
    if not safe:
        print(f"      ^ {label} winds faster than the {margin} rad margin on the "
              f"{N_WL}-point grid. Its dense values are NOT trustworthy; "
              f"raise N_WL rather than reading the result.")
    return safe


def _fit_smooth(x, y, label, deg=2):
    """Degree-2 fit for quantities that are smooth in lambda (n_eff, beta).

    Fitted rather than interpolated for the same reason the parent script
    fits them: these are physically near-quadratic across a device band, so a
    fit both denoises and extrapolates safely, where interpolation would
    inherit every wiggle of the sparse grid.
    """
    p = np.polyfit(x, y, deg)
    resid = float(np.max(np.abs(np.polyval(p, x) - y)))
    print(f"    fit    {label:<34} max resid {resid:.3e}")
    return p


_all_safe = True
_all_safe &= _check_safe(r, "r (HCG mirror)")
for _m in range(N_MODES):
    _all_safe &= _check_safe(M_in0[:, _m, 0], f"M_in0 mode {_m}")
    _all_safe &= _check_safe(M_in1[:, _m, 0], f"M_in1 mode {_m}")

r_dense = interp_complex(wl, r, wl_dense, check_unwrap=0.5)
M_in0_dense = np.stack(
    [interp_complex(wl, M_in0[:, m, 0], wl_dense, check_unwrap=0.5) for m in range(N_MODES)],
    axis=1)[:, :, None]                                    # [n_dense, N_MODES, 1]
M_in1_dense = np.stack(
    [interp_complex(wl, M_in1[:, m, 0], wl_dense, check_unwrap=0.5) for m in range(N_MODES)],
    axis=1)[:, :, None]

beta_body_dense = np.stack(
    [np.polyval(_fit_smooth(wl, beta_body[:, m], f"beta_body mode {m}"), wl_dense)
     for m in range(N_MODES)], axis=1)                     # [n_dense, N_MODES]
neff_arm_dense = np.polyval(_fit_smooth(wl, neff_arm, "n_eff (arm)"), wl_dense)

if not _all_safe:
    print("  !! at least one factor failed its unwrap check -- see above.")

# Pre-transpose the output-side matrices once (mirror-symmetry reuse:
# M_out_through = M_in0^T, exactly as hcg_fabry_perot_hybrid.py derives).
M_in0_d_T = M_in0_dense.transpose(0, 2, 1)                 # [n_dense, 1, N_MODES]
M_in1_d_T = M_in1_dense.transpose(0, 2, 1)

# The arm phase factors out of L_ARM as a pure power: phase(L) = base**L,
# with base = exp(2j*pi*n_eff/wl). Precomputing the exponent lets the arm
# axis be evaluated without re-forming the exponential per cell.
_arm_exponent = 2j * np.pi * neff_arm_dense / wl_dense     # [n_dense]

# ===========================================================================
# 6. The sweep itself -- pure analytic evaluation, no simulation.
# ===========================================================================
coupling_lengths = np.linspace(LC_MIN, LC_MAX, LC_N)
arm_lengths = np.linspace(LARM_MIN, LARM_MAX, LARM_N)


def _cross(wl_grid, T, i_out, i_in, half):
    """Wavelength where T crosses `half`, linearly between samples i_out
    (on/above half) and i_in (below half). Falls back to the sample position
    when the two are numerically equal."""
    t_out, t_in = T[i_out], T[i_in]
    if t_out == t_in:
        return wl_grid[i_out]
    f = (half - t_out) / (t_in - t_out)
    f = min(max(f, 0.0), 1.0)
    return wl_grid[i_out] + f * (wl_grid[i_in] - wl_grid[i_out])


def _find_dips(wl_grid, T, threshold=DIP_THRESHOLD):
    """Local minima + FWHM, same convention as ring_resonator_hybrid.py's
    find_dips() and hcg_fabry_perot_hybrid.py's _find_dips().

    NOTE (2026-08-21): this is now the FIFTH copy of this routine in the
    examples tree (both ring scripts, ring_resonator_fdtd, the FP script,
    here). It is duplicated rather than shared because promoting it to
    beef.spectral is a shared-library change and there were jobs in the queue
    when this file was written. It belongs next to interp_complex/
    unwrap_is_safe -- promote it and collapse all five when the queue is
    clear.

    Returns (wl0, T_min, fwhm) per dip.
    """
    dips = []
    baseline = np.median(T)
    if not np.isfinite(baseline) or baseline <= 0:
        return dips, baseline
    # Vectorised local-minimum detection. The scalar scan this replaces was
    # O(N_dense) PYTHON iterations per cell -- 32M for a 40x40 map at the
    # default density, which dominated the whole sweep. The mask is the same
    # predicate, evaluated in numpy; only the handful of real dips then pay
    # for the (short, local) half-power walks below.
    core = T[1:-1]
    idx = np.flatnonzero((core < T[:-2]) & (core < T[2:]) &
                         (core < threshold * baseline)) + 1
    for i in idx:
        i = int(i)
        half = (T[i] + baseline) / 2.0
        lo = i
        while lo > 0 and T[lo] < half:
            lo -= 1
        hi = i
        while hi < len(T) - 1 and T[hi] < half:
            hi += 1
        if hi > lo:
            # Sub-grid half-power crossings. Without this the FWHM is
            # quantised to the dense-grid step, which is a 2-3 cell staircase
            # at coarse HCGFP_N_DENSE and still a bias at the default -- and
            # FWHM is one of the three reported metrics, so it should not
            # carry an avoidable grid artifact. `lo` sits at the last sample
            # AT OR ABOVE half on the left (crossing in [lo, lo+1]); `hi` at
            # the first AT OR ABOVE half on the right (crossing in [hi-1, hi]).
            wl_lo = _cross(wl_grid, T, lo, lo + 1, half)
            wl_hi = _cross(wl_grid, T, hi, hi - 1, half)
            dips.append((wl_grid[i], T[i], wl_hi - wl_lo))
    return dips, baseline


ER = np.full((LARM_N, LC_N), np.nan)      # dB      -- rows = arm length
FSR = np.full((LARM_N, LC_N), np.nan)     # nm
FWHM = np.full((LARM_N, LC_N), np.nan)    # nm
NDIPS = np.zeros((LARM_N, LC_N), dtype=int)

print(f"\n--- 5. Sweeping {LC_N} x {LARM_N} = {LC_N*LARM_N} cells "
      f"(analytic; no simulation) ---")
_t0 = time.time()

def _coupler_tk(Lc):
    """(t, k) on the dense grid at one coupling length.

    coupling_length enters the whole device ONLY here, and within here only
    through the EME body's diagonal phase.
    """
    body_length = Lc - 2 * BODY_STUB
    M_body_d = np.zeros((N_WL_DENSE, N_MODES, N_MODES), dtype=complex)
    for m in range(N_MODES):
        M_body_d[:, m, m] = np.exp(1j * beta_body_dense[:, m] * body_length)
    _tmp = np.einsum("wij,wjk->wik", M_body_d, M_in0_dense)
    t_d = np.einsum("wij,wjk->wik", M_in0_d_T, _tmp)[:, 0, 0]
    k_d = np.einsum("wij,wjk->wik", M_in1_d_T, _tmp)[:, 0, 0]
    return t_d, k_d


def _cavity_T(t_d, k_d, La):
    """Cavity output |T|^2 on the dense grid at one arm length.

    hcg_fabry_perot_hybrid.py's closed form (its section 6), which that
    script asserts against a real beef.network.Circuit assembly to
    floating-point precision. Arm length enters ONLY through `phase`.
    """
    phase = np.exp(_arm_exponent * La)
    rho2 = (r_dense * phase ** 2) ** 2
    t2, k2 = t_d ** 2, k_d ** 2
    return np.abs(t_d * (1 + rho2 * (k2 - t2)) / (1 - rho2 * t2)) ** 2


def _metrics(T_d):
    """(ER dB, FSR nm, FWHM nm, n_dips) from one dense spectrum."""
    dips, baseline = _find_dips(wl_dense, T_d)
    er = fsr = fwhm = np.nan
    if dips:
        tmin = np.array([d[1] for d in dips])
        er = float(np.median(10 * np.log10(baseline / np.maximum(tmin, 1e-300))))
        fwhm = float(np.median([d[2] for d in dips])) * 1000.0
    if len(dips) >= 2:
        fsr = float(np.mean(np.diff([d[0] for d in dips]))) * 1000.0
    return er, fsr, fwhm, len(dips)


for ic, Lc in enumerate(coupling_lengths):
    t_d, k_d = _coupler_tk(Lc)
    for ia, La in enumerate(arm_lengths):
        ER[ia, ic], FSR[ia, ic], FWHM[ia, ic], NDIPS[ia, ic] = \
            _metrics(_cavity_T(t_d, k_d, La))

    if MASTER and (ic + 1) % max(1, LC_N // 10) == 0:
        print(f"    {ic+1:3d}/{LC_N} coupling lengths  "
              f"({time.time()-_t0:.1f} s elapsed)")

print(f"Sweep complete: {time.time()-_t0:.1f} s")
print(f"  cells with >=1 dip : {int((NDIPS >= 1).sum())}/{NDIPS.size}")
print(f"  cells with >=2 dips: {int((NDIPS >= 2).sum())}/{NDIPS.size}  (FSR needs 2)")
# DENSE-GRID ADEQUACY GUARD. The metrics are read off `wl_dense`, so the
# grid has to resolve a resonance, not merely contain one. FWHM converges
# first (the sub-grid crossings in _cross() carry it); the EXTINCTION RATIO is
# the sensitive one, because it samples the dip BOTTOM rather than its
# flanks -- measured on this device, ER moved 8.35 -> 10.90 -> 10.88 dB going
# 2001 -> 8001 -> 20001 points while FWHM was already converged at 8001.
# An under-resolved grid therefore UNDER-reports depth, silently.
_step_nm = (WL_MAX - WL_MIN) / (N_WL_DENSE - 1) * 1000.0
if np.isfinite(FWHM).any():
    _spf = float(np.nanmin(FWHM)) / _step_nm
    print(f"\nDense-grid adequacy: {_step_nm:.5f} nm/step, "
          f"{_spf:.1f} samples across the NARROWEST resonance found.")
    if _spf < 20:
        print(f"  !! WARNING: fewer than 20 samples across the narrowest "
              f"resonance. Extinction ratio is UNDER-reported at this "
              f"density (it samples the dip minimum) and FWHM may be biased. "
              f"Raise HCGFP_N_DENSE (currently {N_WL_DENSE}) to at least "
              f"{int(np.ceil(20 * (WL_MAX-WL_MIN) * 1000 / float(np.nanmin(FWHM)))) + 1} "
              f"before reading the ER map.")

for name, arr, unit in [("extinction ratio", ER, "dB"), ("FSR", FSR, "nm"),
                        ("FWHM", FWHM, "nm")]:
    if np.isfinite(arr).any():
        print(f"  {name:<17}: {np.nanmin(arr):9.4f} to {np.nanmax(arr):9.4f} {unit}")
    else:
        print(f"  {name:<17}: NO FINITE CELLS -- see the dip-count panel")

# Anchor the sweep to the fabricated device, evaluated EXACTLY at
# (REAL_LC, REAL_LARM) rather than read off the nearest grid cell -- otherwise
# this "check" would silently drift with LC_N/LARM_N and could report a cell
# tens of microns away, which is not a check at all.
_er0, _fsr0, _fwhm0, _nd0 = _metrics(_cavity_T(*_coupler_tk(REAL_LC), REAL_LARM))
print(f"\nFabricated device, evaluated exactly "
      f"(Lc={REAL_LC:.3f} um, L_arm={REAL_LARM:.3f} um):")
print(f"  ER {_er0:.3f} dB   FSR {_fsr0:.4f} nm   FWHM {_fwhm0:.4f} nm   "
      f"({_nd0} dips)")
print("  (hcg_fabry_perot_hybrid.py's own single-point run at Lc=3.0/L_arm=53.08 "
      "on the 1.50-1.60 um band reported mean FSR 2.245 nm -- a DIFFERENT band "
      "and mirror resolution, so agreement in magnitude is the check here, not "
      "equality.)")

# ===========================================================================
# 7. Save + the three heatmaps.
# ===========================================================================
if MASTER:
    np.savez(
        os.path.join(output_dir, "fp_length_sweep.npz"),
        coupling_lengths=coupling_lengths, arm_lengths=arm_lengths,
        extinction_ratio_db=ER, fsr_nm=FSR, fwhm_nm=FWHM, n_dips=NDIPS,
        wavelengths=wl, wl_min=WL_MIN, wl_max=WL_MAX, n_wl=N_WL,
        n_wl_dense=N_WL_DENSE, resolution_hcg=RESOLUTION,
        resolution_dc=FDTD_RESOLUTION, resolution_eme=EME_RESOLUTION,
        resolution_arm=ARM_RESOLUTION, body_stub=BODY_STUB, l_ref=L_REF,
        calibrate_arm=CALIBRATE_ARM, r_hcg=r_hcg, neff_arm=neff_arm,
        real_lc=REAL_LC, real_larm=REAL_LARM,
        real_er_db=_er0, real_fsr_nm=_fsr0, real_fwhm_nm=_fwhm0,
        beta_body=beta_body,
    )
    print(f"\nSweep grid saved -> {output_dir}/fp_length_sweep.npz")

    _extent = [coupling_lengths[0], coupling_lengths[-1],
               arm_lengths[0], arm_lengths[-1]]
    _panels = [
        ("extinction_ratio", ER, "Extinction ratio (dB)", "magma",
         "median over resonances in band"),
        ("fsr", FSR, "FSR (nm)", "viridis",
         "mean adjacent-dip spacing; NaN where <2 dips fall in band"),
        ("fwhm", FWHM, "FWHM (nm)", "cividis",
         "median over resonances in band"),
    ]
    for fname, arr, label, cmap, sub in _panels:
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        im = ax.imshow(arr, origin="lower", aspect="auto", extent=_extent,
                       cmap=cmap, interpolation="nearest")
        cb = fig.colorbar(im, ax=ax)
        cb.set_label(label)
        # Mark the real fabricated device.
        ax.plot([REAL_LC], [REAL_LARM], marker="*", ms=14, mfc="white",
                mec="k", mew=1.0, ls="none", label="fabricated device")
        ax.set_xlabel("Directional coupler coupling length (um)")
        ax.set_ylabel("Connecting waveguide arm length (um)")
        ax.set_title(f"HCG Fabry-Perot — {label}\n{sub}", fontsize=10)
        ax.legend(fontsize=7, loc="upper right", framealpha=0.85)
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, f"heatmap_{fname}.png"), dpi=300)
        plt.close(fig)
        print(f"Heatmap saved -> {output_dir}/heatmap_{fname}.png")

    # Diagnostic panel: the dip count is what explains every NaN above.
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    im = ax.imshow(NDIPS, origin="lower", aspect="auto", extent=_extent,
                   cmap="Greys", interpolation="nearest")
    fig.colorbar(im, ax=ax).set_label("resonances found in band")
    ax.plot([REAL_LC], [REAL_LARM], marker="*", ms=14, mfc="white", mec="k",
            mew=1.0, ls="none")
    ax.set_xlabel("Directional coupler coupling length (um)")
    ax.set_ylabel("Connecting waveguide arm length (um)")
    ax.set_title("Resonances per cell — the reason any cell above is NaN\n"
                 f"(FSR needs >=2; band is {(WL_MAX-WL_MIN)*1000:.0f} nm wide)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "heatmap_dip_count.png"), dpi=300)
    plt.close(fig)
    print(f"Diagnostic saved -> {output_dir}/heatmap_dip_count.png")

print("\nDone.")
