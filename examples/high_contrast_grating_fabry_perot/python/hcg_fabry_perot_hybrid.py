"""Tutorial: the HCG Fabry-Perot resonator -- two HCG mirrors terminating a
cavity, tapped in the middle by a directional coupler.

This is the first script to actually cascade three INDEPENDENTLY solved
devices via beef.network.Circuit rather than composing sections within one
HybridSolver chain (that doc's A2) -- and the first to build a genuinely
1-port self-reflecting block (A3), since beef.network.convert.from_smatrix()
cannot represent one by its own documented design (self-terms are always
zero-filled there).

SELF-CONTAINED, 2026-08-18 (rewritten). Every simulation this cavity is made
of now runs HERE, through the library, against example_sim_cache/ -- nothing
is read out of another script's example_output/*.pkl any more. The three
physical pieces are rebuilt from their own geometry and solver settings:

  1. HCG mirror -- FDTDSolver.solve_reflection() on
     examples/gds_files/HCG_seeded_TO.gds, res 60, 3-D, 21 points over
     1.50-1.60 um. Byte-for-byte the call
     examples/high_contrast_grating/python/high_contrast_grating_fdtd.py
     makes, so it hits that script's own cache entry.
  2. Directional coupler -- two FDTDSolver.solve_section() S-bend
     excitations + a beta-CALIBRATED EMESolver body, on SOI_160 at
     duty/dc.gds's real dimensions. Same calls, same constants as
     examples/directional_coupler/python/directional_coupler_soi160_hybrid.py,
     so those hit cache too.
  3. Connecting waveguide -- NEW. A real FDTD beta calibration of the bare
     0.5 um strip (see CALIBRATE_ARM below), replacing the MPB mode index
     the arms used to propagate on.

Why this rewrite. Until now this script consumed two saved pickles, which
made the cavity's answer depend on which revision of two other scripts had
last written them -- and it hid the fact that the two 53.08 um arms, the
single longest optical path in the device, were the ONE piece still
propagating on MPB's mode index while everything they cascade with had been
moved onto meep's own beta. Rebuilding the pieces in place fixes both: the
provenance and the calibration are now visible in one file. As a regression
guard the rebuilt HCG r(lambda) and coupler t/k(lambda) are compared against
those saved pickles when they happen to still be on disk (best-effort, not
required) -- they must agree to machine precision.

CALIBRATE_ARM (default True) -- the point of the rewrite.

`section_mode_beta_details()` measures the beta meep's time-stepper actually
accumulates through a cross-section and finds it ~2.3% above MPB's mode
index at res 40 (see that method's own docstring, and dc_validation.md item
6). Every other block in this cavity is already referenced to that meep
beta: the HCG's r(lambda) is a meep mode-monitor measurement, the coupler's
S-bends are meep, and its EME body propagates on its own FDTD-calibrated
beta (CALIBRATE_BODY below). The arms were not, and they are 2 x 53.08 um
against the coupler's 23 um -- so ~106 um of the ~130 um one-way path was
carried by a different engine's index than the rest of the loop.

That is not a small correction here. The resonant denominator is
1 - rho^2 t^2 with rho = r * phase^2, so n_eff acts over a FULL 4 * L_ARM =
212.3 um round trip; a ~2% index change is worth many free spectral ranges
of comb position. The script prints the measured excess and exactly what it
is worth in FSRs rather than leaving it as an argument.

HONEST LIMIT, stated rather than buried: calibrating does NOT make the
absolute resonance positions converged. It makes the loop self-consistent --
one engine's discretization carrying every segment -- which is the thing a
closed round trip actually requires. Both MPB's index and meep's beta are
discretizations of the same continuum problem and neither is the truth; the
residual is a real, unmeasured error on absolute comb position. FSR, Q and
depth are all ratios/differences and move only by the same ~2%.
HCGFP_CALIBRATE_ARM=0 puts the arms back on an MPB index, into its own output
folder, so the two can be read side by side. It is NOT a bit-identical replay
of the pre-2026-08-18 script: that one solved the index with meow's own
Mesh2D at res 100, this one asks meep's eigensolver on the FDTD grid the
calibration also runs on. That is deliberate -- holding the grid fixed makes
the difference between the two curves a statement about the two ENGINES
rather than about two meshes. The archived original is at
archive/example_output_20260818/.

HCGFP_CALIBRATE_BODY=0 does the same for the coupler BODY, and the two
compose: with both off the round trip carries no FDTD-measured beta at all,
which is the pre-calibration baseline this device's whole rewrite was
measured against. Note the asymmetry in what that costs. Turning the ARM
calibration off makes the run CHEAPER (one MPB sweep replaces an FDTD fit).
Turning the BODY calibration off makes it DEARER: `neff_override=None` gives
EMESolver nothing to propagate on, so the diagonal fast path is gone and
the real eigenmode solve runs.

ARM_RESOLUTION defaults to 80, matching the HCG MIRROR's res rather than the
coupler's res 40 -- the user's explicit call. There is no single device resolution to
inherit here (the mirror runs at 60, the coupler at 40), and a bare strip is
cheap enough at 60 that the arm's beta -- the dominant lever on the comb --
is measured on the finer of the two meshes this device already uses. The
stated cost: the arm's discretization is then NOT identical to the coupler
body's, so a milder version of the inconsistency above survives, over 23 um
of coupler instead of 106 um of arm. HCGFP_RES / HCGFP_DC_RES /
HCGFP_ARM_RES override the three and redirect output to a res-stamped folder
("res80_dc40_arm80", named the way hcg_fabry_perot_length_sweep.py names its
own), so a real convergence sweep can be run without touching this.

Symmetric-coupler assumption, VERIFIED not just assumed (2026-08-07): a
real reciprocal, physically symmetric 4-port coupler has t'==t, k'==k
(the "isolated"-excited response equals the "in"-excited response) --
confirmed directly from this device's own M_in0/M_in1/M_body pieces
(k' matched k to machine precision, t' matched t to ~1e-5). So the "in"-
excited t/k alone are enough to build the full 4-port block; a separate
reconstruction of the "isolated" branch is not needed.

Two independent computations of the same spectrum, asserted to agree:
  1. The real beef.network.Circuit assembly (5 blocks, 4 connect() calls) --
     this device's actual topology.
  2. The hand-derived closed-form Fabry-Perot-with-DC-tap formula (the
     ring resonator's own closed-form-check convention, generalized to a
     4-port tap instead of one closed ring port; independently verified
     against a brute-force 4-port linear solve in
     tests/test_fabry_perot_network_crosscheck.py).
Matching to floating-point precision is the check itself.

The connecting waveguide block is built as a genuinely BIDIRECTIONAL
SMatrixBlock -- NOT via beef.network.convert.from_transfer_matrix(), which
is documented one-directional by design and was found (2026-08-05,
tests/test_fabry_perot_network_crosscheck.py) to silently zero out a
mirror's reflection when used for this.

KNOWN OPEN GAPS, not fixed by this script (see the validation doc):
  - FDTD_RESOLUTION=40/EME_RESOLUTION=80 (inherited from
    directional_coupler_soi160_hybrid.py) are ASSUMED for the coupler, not
    independently converged for this device.
  - L_ARM=53.08 um is the AVERAGE of the real structure's own two arms
    (54.583 and 51.583 um -- a real 3.0 um asymmetry this symmetric-cavity
    model does not capture, per Step 1's own framing).
  - The arm beta is measured at ONE resolution per run; nothing here
    converges it. HCGFP_ARM_RES is the hook for that sweep.
  - The resolution knobs do NOT move the wavelength grid. Pinning res to the
    2026-08-21 C-band re-runs (80/40/80) does not by itself make those runs'
    cache entries hit: they were measured on 1.530-1.565 um / 36 pts and this
    script asks for 1.50-1.60 um / 11 (21 for the mirror). Probe with
    HCGFP_CACHE_ONLY=1 before spending an allocation.

DENSE RESONANCE CHARACTERIZATION (section 8)

The 11-point/10 nm coupler grid is far coarser than the ~0.41 nm FWHM
estimated from this cavity's own finesse -- the sparse spectrum only shows
that a resonance exists, not its real depth/width/Q/position. Brute force
does not work: job 11726827 tried a 601-point re-solve of the coupler and
TIMED OUT at 12 h after 3 of ~7 sub-solves, because each solve's end-of-run
eigenmode decomposition scales with N_WL (11 -> 601 is ~55x on that piece,
not "the same run with more points").

Fixed with interpolation instead, and now done IN THIS FILE rather than by
loading a separate script's dense pickle. Each FACTOR is resampled where it
is safe to do so and recombined analytically -- the same scheme
directional_coupler_soi160_hybrid_dense.py derived and validated
per-factor against real 601-point ground truth (<=1.3% relative error):

    M_in0/M_in1 (arm -> body-mode overlap):  ~65-69 deg/10nm  -> interpolate
    beta_body (coupler body):                ~15-16 deg/10nm  -> degree-2 fit
    n_eff_arm (connecting waveguide):        smooth           -> degree-2 fit
    r (HCG mirror, on its own 21-pt grid):   ~7.5 deg/5nm     -> interpolate
    coeff_through/cross (their COMBINATION): ~153 deg/10nm    -> NEVER resampled

The last line is the whole reason this is done per-factor: the combined
coupler coefficient fails beef.spectral.unwrap_is_safe outright, because
this device's mirror-symmetry reuse (M_out_through = M_in0^T) counts the
input arm's phase twice. Every factor that IS resampled is checked with
unwrap_is_safe/check_unwrap rather than assumed safe.

Run:
    python /path/to/examples/high_contrast_grating_fabry_perot/python/hcg_fabry_perot_hybrid.py
    HCGFP_CALIBRATE_ARM=0 python .../hcg_fabry_perot_hybrid.py   # MPB-index arms, own folder
    HCGFP_ARM_RES=40 python .../hcg_fabry_perot_hybrid.py        # arm-beta res sweep, own folder
    HCGFP_CALIBRATE_BODY=0 python .../hcg_fabry_perot_hybrid.py  # MPB-index body, own folder
    HCGFP_RES=80 HCGFP_DC_RES=40 HCGFP_ARM_RES=80 \
        python .../hcg_fabry_perot_hybrid.py                     # -> .../res80_dc40_arm80/
"""

import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection

from beef import (
    DirectionalCoupler, EMESolver, FDTDSolver, GdsGeometry, HybridSolver,
    SMatrix, SOI_160, load_gds_polygons,
    plot_beta_calibration, plot_hybrid_geometry,
)
from beef.cache import is_master
from beef.geometry import UniformCrossSection
from beef.network import Circuit
from beef.network.block import SMatrixBlock
from beef.spectral import interp_complex, unwrap_is_safe

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
_GDS_PATH = os.path.abspath(os.path.join(_HERE, "..", "..", "gds_files", "HCG_seeded_TO.gds"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_sim_cache"))

# Every write below goes through this. From 2026-08-18 this script does real
# FDTD (it used to be pure post-processing of two pickles) and so runs under
# `srun -n 24` -- every rank executes the whole file, and unguarded savefig/
# save calls would have 24 processes writing the same paths on a shared
# filesystem. beef.cache elects exactly one writer; reuse that election
# rather than inventing a second one.
MASTER = is_master()

# cache_only: demand a cache hit instead of falling back to simulating.
# The companion notebook sets this so a
# documentation re-render can never kick off a multi-hour FDTD run; the
# script defaults it off so it can populate the cache in the first place.
CACHE_ONLY = bool(int(os.environ.get("HCGFP_CACHE_ONLY", "0")))

# ---------------------------------------------------------------------------
# Global parameters
# ---------------------------------------------------------------------------
# BAND, 2026-08-24: the real C-band operating window, 1.530-1.565 um / 36 pts,
# replacing the earlier 1.500-1.600 um / 11+21. The HCG mirror is a resonant
# reflector seeded for the C band (Hiesener/Pesch/Ralph, IEEE PTL 36(9), 2024) --
# NOT the O-band the rest of this repo's SOI_160 examples use.
#
# This is the SAME grid the 2026-08-21 measurement re-runs used, so the mirror and
# the calibrated arm are cache HITS out of the box rather than 4 h of re-simulation.
#
# STATED COST of the narrower window, so it is not discovered later: the mirror's
# own reflectivity PEAK sits at 1.580 um (|r|^2 = 0.9714), OUTSIDE this band. On
# 1.530-1.565 um the mirror is on the shoulder of its resonance, not the top of it.
# That is the correct choice for a device specified over the C band, but it means
# the cavity finesse reported here is NOT the best this mirror can do.
WL_MIN = 1.530
WL_MAX = 1.565
N_WL = 36           # ONE grid for every section now -- mirror, coupler and arm alike
N_WL_HCG = N_WL     # the mirror no longer runs a finer grid of its own
N_WL_DENSE = 20001  # dense post-processing grid (section 8) -- costs no simulation

# Every resolution is an env knob, same shape as hcg_fabry_perot_length_sweep.py's
# (HCGFP_RES / HCGFP_DC_RES / HCGFP_ARM_RES) so the two scripts can be pinned to the
# same grid and share cache entries.
#
# Defaults raised 60 -> 80 on the mirror and the arm (2026-08-24), so that an un-set
# environment now means "the current production resolutions" rather than "the
# 2026-08-18 ones". This makes the DEFAULTS match hcg_fabry_perot_length_sweep.py's
# own defaults exactly. The superseded res-60 outputs are archived under
# archive/example_output_20260824/.
RESOLUTION = float(os.environ.get("HCGFP_RES", 80.0))          # HCG mirror FDTD
FDTD_RESOLUTION = float(os.environ.get("HCGFP_DC_RES", 80.0))  # coupler S-bend FDTD
EME_RESOLUTION = float(os.environ.get("HCGFP_EME_RES", 80.0))  # coupler EME body
ARM_RESOLUTION = float(os.environ.get("HCGFP_ARM_RES", 80.0))  # arm beta calibration

N_MODES = 2   # coupler body: fundamental even/odd supermode pair ONLY. MEASURED
              # simplification (<0.01% of launched power reaches the 2nd-order
              # pair) -- see directional_coupler_soi160_hybrid.py's docstring.

# coupler EME body propagates on the FDTD's own beta. HCGFP_CALIBRATE_BODY=0 puts it
# back on the EME's own eigenmode solve -- and that is NOT free: with no neff_override
# the diagonal fast path is gone and EMESolver runs the real mode solve.
CALIBRATE_BODY = bool(int(os.environ.get("HCGFP_CALIBRATE_BODY", "1")))
CALIBRATE_ARM = bool(int(os.environ.get("HCGFP_CALIBRATE_ARM", "1")))

# Derive the "isolated"-arm excitation from the "in"-arm one instead of running
# a second FDTD solve, exactly as directional_coupler_hybrid.py does (that
# script's SYMMETRIC_MODE_REUSE; it is NOT an MMI feature). Launching "in" vs
# "isolated" is itself a mirror operation on a cross-section that is mirror-
# symmetric about y=0, so each body mode's coefficient can only keep or flip
# sign: M_in1[m] = MODE_SIGNS[m] * M_in0[m].
#
# MEASURED ON THIS DEVICE, 2026-08-25, against the real two-arm FDTD pair
# already cached at res 80 on 1.530-1.565 um / 36 pts -- not inherited from the
# SOI_220 script and not assumed:
#
#     mode 0 (even)  ratio +0.999963 - 0.000008j   reuse error 5.4e-05
#     mode 1 (odd)   ratio -1.000000 + 0.000001j   reuse error 1.1e-06
#
# with the ratio's argument spread ~1e-05 rad, i.e. a clean pure sign.
#
# WHY IT IS SAFER HERE THAN ON THE SOI_220 COUPLER: parity constrains the SIGN,
# never the amplitude, so a near-cutoff mode of ~1e-6 magnitude can disagree
# wildly while contributing nothing -- which is why that script's three
# higher-order modes are ~1% noisy. This device carries N_MODES=2 and BOTH are
# dominant (|M| ~ 0.707 each, the 50/50 split a symmetric coupler must give),
# so that failure mode is not in play at all.
#
# WHAT IT SAVES: one of the coupler's TWO S-bend solves. It does NOT touch the
# body beta calibration's three runs (a slope fit over three different
# propagation lengths, which no cross-section symmetry collapses), so the
# coupler goes 5 FDTD runs -> 4. In wall clock it is worth more than that ratio
# suggests: at res 80 each S-bend took 335 min against 228 min for all three
# body runs together, so this removes ~37% of the coupler's total.
#
# HCGFP_SYMMETRIC_MODE_REUSE=0 runs the second arm for real and re-derives the
# signs, which is how the numbers above were produced and how to re-derive them
# if the geometry ever changes.
SYMMETRIC_MODE_REUSE = bool(int(os.environ.get("HCGFP_SYMMETRIC_MODE_REUSE", "1")))
MODE_SIGNS = np.array([1, -1])

COUPLING_LENGTH = 3.0   # um -- duty/dc.gds's real straight coupled-pair length
BODY_STUB = 0.25        # um of that length folded into each FDTD section

L_ARM = 53.08   # um -- Step 1's REAL pin, 2026-08-08, traced from duty/Mux_lpbk_9$4.gds's own
                # real fabricated mirror-DC-mirror loopback structure (a real test structure with
                # this device's own topology), NOT the earlier 20 um placeholder. Path traced
                # instance-by-instance (gdstk reference transforms: origin/rotation/x_reflection),
                # left arm = straight(25.375) + bend(R=10.0, quarter-circle arc len=15.708,
                # user's instruction: treat bends as straight waveguide) + straight(13.5) =
                # 54.583 um; right arm = straight(25.375) + bend(15.708) + straight(10.5) =
                # 51.583 um. R=10.0 confirmed two independent ways: from the bend cell's own core
                # polygon extent, and from each bend's own port-to-port chord length (chord=R*sqrt(2)
                # for a 90 deg arc) -- both give exactly 10.0. The two real arms differ by exactly
                # 3.0 um (a real layout asymmetry, not modeled here -- this device's own symmetric-
                # cavity assumption, Step 1, uses one shared length for both arms); L_ARM here is
                # their average, 53.083 -> 53.08.
WG_WIDTH = 0.5   # um -- the same real strip width as the HCG's own access waveguide and the DC's arms

# Longitudinal budget the arm calibration cuts its fit lengths from. Same
# reasoning (and same value) as ring_resonator_hybrid.py's CALIBRATE_SPAN:
# section_mode_beta() picks [s, 2s, 3s] with s = min((span-2*dpml)/3,
# max_safe_spacing/2), and at beta ~ 9 rad/um the physics cap is ~1.4 um,
# which needs span >= 5.0 not to be the binding constraint. The strip is
# uniform, so this is a budget, not a claim about the real arm's geometry.
ARM_CAL_SPAN = 6.0
MESH_WIDTH = 4.0                     # transverse window for the arm's cross-section
SIM_HEIGHT = 1.5                     # slab z-height reference, shared by all three pieces
DPML = 0.5

output_dir = os.path.join(_ROOT, "example_output", "high_contrast_grating_fabry_perot",
                          "hcg_fabry_perot_hybrid")
# Escape hatches redirect output, same shape as mmi_1x2_hybrid.py's
# MMI1X2_CALIBRATE_BODY and ring_resonator_hybrid.py's RING_FDTD_RES: "the
# guard is the env var, not the value", so a diagnostic run can never
# overwrite the production spectrum.
if not CALIBRATE_ARM:
    output_dir += "_uncalibrated_arm"
if not CALIBRATE_BODY:
    output_dir += "_uncalibrated_body"
# Any resolution knob set at all -> its own res-stamped subfolder, named exactly the
# way hcg_fabry_perot_length_sweep.py names its own ("res80_dc40_arm80"). The guard is
# the env var, not the value, so a run pinned to the defaults still lands somewhere it
# cannot overwrite the production spectrum.
if any(os.environ.get(_v) is not None
       for _v in ("HCGFP_RES", "HCGFP_DC_RES", "HCGFP_EME_RES", "HCGFP_ARM_RES")):
    output_dir = os.path.join(
        output_dir, f"res{RESOLUTION:g}_dc{FDTD_RESOLUTION:g}_arm{ARM_RESOLUTION:g}")
if MASTER:
    os.makedirs(output_dir, exist_ok=True)

wl_hcg = np.linspace(WL_MIN, WL_MAX, N_WL_HCG)
wl = np.linspace(WL_MIN, WL_MAX, N_WL)          # the assembly grid
M = len(wl)

print(f"HCG Fabry-Perot resonator -- self-contained rebuild")
print(f"  CALIBRATE_BODY = {CALIBRATE_BODY} (coupler EME body, res {EME_RESOLUTION:g})")
print(f"  resolutions    : mirror {RESOLUTION:g}, coupler S-bend {FDTD_RESOLUTION:g}, "
      f"coupler body {EME_RESOLUTION:g}, arm {ARM_RESOLUTION:g}")
print(f"  CALIBRATE_ARM  = {CALIBRATE_ARM} (connecting waveguide, res {ARM_RESOLUTION:g})")
print(f"  SYMMETRIC_MODE_REUSE = {SYMMETRIC_MODE_REUSE} "
      f"({'1' if SYMMETRIC_MODE_REUSE else '2'} coupler S-bend FDTD solve"
      f"{'' if SYMMETRIC_MODE_REUSE else 's'})")
print(f"  CACHE_ONLY     = {CACHE_ONLY}")
print(f"  output -> {output_dir}")

# ---------------------------------------------------------------------------
# 1. The real HCG mirror -- FDTDSolver.solve_reflection(), run here.
#
#    Reproduces high_contrast_grating_fdtd.py's production call exactly
#    (geometry from the same GDS, res 60, 3-D, same PML/margins/clearance,
#    extend_input_waveguide=True), so it hits that script's own cache entry
#    rather than re-simulating. extend_input_waveguide runs the access
#    waveguide through the PML: without it the guide ends in a clean-cut
#    Si/oxide facet that reflects on its own, contaminating the very
#    quantity being measured.
# ---------------------------------------------------------------------------
CORE_LAYERS = [(2, 0)]
EDGE_CLEARANCE = 0.3    # plain-oxide gap between the GDS bbox and the PML, both sides
SOURCE_MARGIN = 0.45    # measured inward from the PML inner edge -- both land inside
MONITOR_MARGIN = 0.8    # the GDS's own 1.082 um input stub
HCG_SIM_WIDTH = 7.5     # y-span; the grating's arcs reach y = +-2.99 um

_raw = load_gds_polygons(_GDS_PATH, layers=CORE_LAYERS)
_xs = [x for poly in _raw[CORE_LAYERS[0]] for x, _y in poly]
hcg_geom = GdsGeometry(
    gds_path=_GDS_PATH,
    core_layers=CORE_LAYERS,
    region_lengths=[max(_xs) - min(_xs)],   # derived from the file, not hardcoded
    input_width=WG_WIDTH,
    input_y=0.0,
)
HCG_LENGTH = hcg_geom.total_length

hcg_solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=HCG_SIM_WIDTH,
    sim_height=SIM_HEIGHT,
    dpml=DPML,
    source_margin=SOURCE_MARGIN,
    monitor_margin=MONITOR_MARGIN,
    three_d=True,
    num_input_modes=1,
    material_stack=SOI_160,
    decay_tol=1e-11,
)

print(f"\n--- 1. HCG mirror (FDTD, res {RESOLUTION:g}, 3-D, {N_WL_HCG} pts) ---")
print(f"HCG total length (x-span of imported polygons): {HCG_LENGTH:.4f} um")
_t0 = time.time()
smat_hcg = hcg_solver.solve_reflection(
    hcg_geom, wl_hcg, port="top", input_mode=1,
    edge_clearance=EDGE_CLEARANCE, extend_input_waveguide=True,
    save_dir=save_dir, cache_only=CACHE_ONLY,
)
smat_hcg.check_passive()
r_hcg = smat_hcg.coeff(port="in", in_mode=1, out_mode=1)
R_hcg = np.abs(r_hcg) ** 2
print(f"HCG mirror: {N_WL_HCG} pts, {wl_hcg[0]:.3f}-{wl_hcg[-1]:.3f} um in "
      f"{(time.time()-_t0)/60:.1f} min -- peak |r|^2={R_hcg.max():.4f} at "
      f"{wl_hcg[np.argmax(R_hcg)]:.3f} um, min {R_hcg.min():.4f}")

# ---------------------------------------------------------------------------
# 2. The real directional coupler -- 2 FDTD S-bend excitations + a
#    beta-calibrated EME body, run here.
#
#    Same geometry and same solver constants as
#    directional_coupler_soi160_hybrid.py, so all three expensive pieces hit
#    that script's cache. The output section is never separately simulated:
#    M_in0/M_in1 double as the output-side transfer matrices once transposed
#    (mirror-symmetry reuse, identical reasoning to that script).
# ---------------------------------------------------------------------------
print(f"\n--- 2. Directional coupler (FDTD S-bends + "
      f"{'calibrated ' if CALIBRATE_BODY else ''}EME body, {N_WL} pts) ---")

dc = DirectionalCoupler(
    waveguide_width=WG_WIDTH,          # duty/dc.gds, step 0 -- single-mode strip width
    gap=0.25,                          # duty/dc.gds -- coupling-region edge-to-edge gap
    coupling_length=COUPLING_LENGTH,   # duty/dc.gds -- FIXED, not swept
    bend_length=10.0,                  # duty/dc.gds
    waveguide_offset=4.0,              # duty/dc.gds -- far-end centre-to-centre port pitch
    bend_type="euler",                 # ASSUMED -- not re-derived from the GDS's curve family
    waveguide_height=SOI_160.wg_height,
)
DC_LENGTH = dc.total_length
print(f"Coupler (SOI_160, duty/dc.gds dims) total length: {DC_LENGTH:.2f} um "
      f"(coupling pitch {dc.coupling_pitch:.3f} um, far pitch {dc.waveguide_offset:.2f} um)")

# waveguide_offset=4.0 needs sim_width 7.0 to keep directional_coupler_hybrid.py's
# own PML clearance ratio (half_domain_aperture -> each arm needs
# waveguide_offset/2 + clearance inside the non-PML half-domain).
dc_fdtd = FDTDSolver(
    res=FDTD_RESOLUTION,
    sim_width=7.0,
    sim_height=SIM_HEIGHT,
    dpml=DPML,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    half_domain_aperture=True,
    material_stack=SOI_160,
    decay_tol=1e-11,
)
dc_eme = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],   # length set per call below
    oxide_width=7.0, oxide_height=SIM_HEIGHT, res=EME_RESOLUTION,
    material_stack=SOI_160, mode_solver="meep",
)

half_body = dc.coupling_length / 2
half_total = dc.total_length / 2
boundaries = [
    -half_total,               # physical input port (access waveguide end)
    -half_body + BODY_STUB,    # input section <-> body section interface
    +half_body - BODY_STUB,    # body section <-> output section interface
    +half_total,               # physical output port
]
body_length = boundaries[2] - boundaries[1]
assert body_length > 0, (
    f"BODY_STUB ({BODY_STUB:.2f} um) too large for COUPLING_LENGTH "
    f"({COUPLING_LENGTH:.2f} um) -- the EME body length would be negative"
)
print(f"Coupling-length budget: 2 x {BODY_STUB:.2f} um inside the FDTD sections "
      f"+ EME body {body_length:.2f} um = {COUPLING_LENGTH:.2f} um total")

_iface_mon = dict(monitor_center=0.0, monitor_size=dc_fdtd.sim_width)

_t0 = time.time()
M_in0 = dc_fdtd.solve_section(
    dc, wl, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=dc_fdtd.waveguide_length, stub_right=0.0,
    input_idx=0, save_dir=save_dir, cache_only=CACHE_ONLY, **_iface_mon,
)  # [n_wl, N_MODES, 1] -- "in" arm -> body modes
if SYMMETRIC_MODE_REUSE:
    # The sign vector is a property of THIS mode count on THIS mirror-symmetric
    # cross-section. Raising N_MODES without re-deriving it would silently apply
    # a two-entry parity to modes nobody measured, so refuse rather than guess.
    if MODE_SIGNS.shape != (N_MODES,):
        raise ValueError(
            f"MODE_SIGNS has {MODE_SIGNS.shape[0]} entries but N_MODES={N_MODES} "
            f"-- re-derive it with HCGFP_SYMMETRIC_MODE_REUSE=0 before reusing."
        )
    M_in1 = M_in0 * MODE_SIGNS[None, :, None]
    print(f"Coupler S-bend excitations: {(time.time()-_t0)/60:.1f} min "
          f"(1 FDTD solve; the isolated arm reused via MODE_SIGNS="
          f"{MODE_SIGNS.tolist()})")
else:
    M_in1 = dc_fdtd.solve_section(
        dc, wl, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=dc_fdtd.waveguide_length, stub_right=0.0,
        input_idx=1, save_dir=save_dir, cache_only=CACHE_ONLY, **_iface_mon,
    )  # [n_wl, N_MODES, 1] -- "isolated" arm -> body modes
    print(f"Coupler S-bend excitations: {(time.time()-_t0)/60:.1f} min "
          f"(2 FDTD solves, both arms simulated)")
    # Re-derivation, printed so the constant above can be checked (or updated)
    # against a real two-arm run rather than trusted.
    _ratio = (M_in1[:, :, 0] / M_in0[:, :, 0])
    _signs = np.sign(_ratio.real.mean(axis=0)).astype(int)
    _err = np.abs(M_in0[:, :, 0] * _signs - M_in1[:, :, 0]) / np.abs(M_in1[:, :, 0])
    print(f"MODE_SIGNS re-derivation (HCGFP_SYMMETRIC_MODE_REUSE=0):")
    for _m in range(N_MODES):
        print(f"    mode {_m}: ratio {_ratio[:, _m].mean():+.6f}  "
              f"reuse error max {_err[:, _m].max():.2e}  "
              f"|M| {np.abs(M_in0[M//2, _m, 0]):.4f}")
    print(f"    implied MODE_SIGNS = {_signs.tolist()}  "
          f"(file has {MODE_SIGNS.tolist()})")

# The body's own beta. _details() rather than section_mode_beta(): same call,
# same cost, same guards -- it just also hands back the measurement behind the
# answer, so the calibration can be plotted below.
body_neff = None
cal_body = None
if CALIBRATE_BODY:
    _t0 = time.time()
    cal_body = dc_fdtd.section_mode_beta_details(
        dc, wl, boundaries[1], boundaries[2], N_MODES,
        save_dir=save_dir, cache_only=CACHE_ONLY, **_iface_mon,
    )
    body_neff = cal_body.beta * wl[:, None] / (2 * np.pi)   # [n_wl, N_MODES]
    # MPB's own answer comes back on the calibration object (`seed_beta` -- the
    # branch prior, never the result), so the comparison below costs NOTHING.
    # A separate section_mode_neff() sweep would be a second MPB solve of the
    # identical cross-section, and on this device those are not cheap: the
    # validation doc records a 4 h budget consumed inside two of them.
    _mpb_body = cal_body.seed_beta * wl[:, None] / (2 * np.pi)
    _i = int(np.argmin(np.abs(wl - 1.55)))
    print(f"Body beta calibration at lambda = {wl[_i]:.3f} um "
          f"({(time.time()-_t0)/60:.1f} min), even/odd supermodes:")
    print(f"    MPB mode index  n_eff {_mpb_body[_i, 0]:.6f} {_mpb_body[_i, 1]:.6f}   "
          f"delta_n_eff {_mpb_body[_i, 0] - _mpb_body[_i, 1]:.6f}")
    print(f"    FDTD own phase  n_eff {body_neff[_i, 0]:.6f} {body_neff[_i, 1]:.6f}   "
          f"delta_n_eff {body_neff[_i, 0] - body_neff[_i, 1]:.6f}   <- body propagates on this")

# With neff_override the transfer matrix is diagonal and determined entirely
# by that array, so EMESolver skips the eigenmode solve outright.
M_body = dc_eme.solve_section(
    dc, wl, boundaries[1], boundaries[2],
    n_modes_in=N_MODES, n_modes_out=N_MODES,
    neff_override=body_neff, save_dir=save_dir, cache_only=CACHE_ONLY,
)  # [n_wl, N_MODES, N_MODES]

M_out_through = M_in0.transpose(0, 2, 1)  # [n_wl, 1, N_MODES] -- mirror-symmetry reuse
M_out_cross = M_in1.transpose(0, 2, 1)

_tmp = np.einsum("wij,wjk->wik", M_body, M_in0)              # [n_wl, N_MODES, 1]
t = np.einsum("wij,wjk->wik", M_out_through, _tmp)[:, 0, 0]  # complex "through"
k = np.einsum("wij,wjk->wik", M_out_cross, _tmp)[:, 0, 0]    # complex "cross"
print(f"Coupler at band centre: |t|^2={np.abs(t[M//2])**2:.4f} "
      f"|k|^2={np.abs(k[M//2])**2:.4f} total={np.abs(t[M//2])**2+np.abs(k[M//2])**2:.4f}")

# ---------------------------------------------------------------------------
# 3. The connecting waveguide -- FDTD beta calibration of the bare strip.
#
#    WHY UniformCrossSection + FDTDSolver AND NOT A HAND-ROLLED MEEP RUN
#    (the same reasoning ring_resonator_hybrid.py records for its own
#    straights): beef/geometry/uniform.py exists precisely so a bare strip
#    is a DeviceGeometry both engines can consume. Going through
#    section_mode_beta_details() inherits, for free, the merged
#    multi-monitor domain, the diagonal-only decode, the cached 3-anchor MPB
#    seed, the 2*pi branch guard, the phase-linearity residual guard, and
#    plot_beta_calibration()'s figures. Re-deriving any of that here would
#    duplicate tested library code and produce a beta with none of the
#    guards attached.
#
#    Fundamental mode only: Step 0's own S-bend measurement found <0.01% of
#    excited power reaches the 2nd-order mode pair, so the arm only ever
#    carries mode 1.
# ---------------------------------------------------------------------------
print(f"\n--- 3. Connecting waveguide ({WG_WIDTH:.2f} um strip, "
      f"{L_ARM:.2f} um per arm) ---")

arm_geom = UniformCrossSection(
    [(0.0, WG_WIDTH)], length=ARM_CAL_SPAN, waveguide_height=SOI_160.wg_height,
)
arm_solver = FDTDSolver(
    res=ARM_RESOLUTION,
    # MESH_WIDTH, so the calibrated n_eff and the MPB n_eff below describe the
    # same transverse window and the excess between them is a statement about
    # the two ENGINES rather than about two domain sizes.
    sim_width=MESH_WIDTH,
    sim_height=SIM_HEIGHT,
    dpml=DPML,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    num_input_modes=1,
    num_output_modes=1,
    output_monitor_centers=[0.0],
    material_stack=SOI_160,
    decay_tol=1e-11,
)

cal_arm = None
if CALIBRATE_ARM:
    _t0 = time.time()
    cal_arm = arm_solver.section_mode_beta_details(
        arm_geom, wl, -ARM_CAL_SPAN / 2.0, +ARM_CAL_SPAN / 2.0, 1,
        save_dir=save_dir, cache_only=CACHE_ONLY,
    )
    neff_arm = cal_arm.beta[:, 0] * wl / (2 * np.pi)
    # MPB's index rides along on the calibration object (its branch prior),
    # so the comparison below costs no extra mode solve. It is solved at
    # seed_anchors=3 wavelengths and interpolated across the band -- accurate
    # to ~1e-4 in n_eff, which is a diagnostic-grade reference, not a second
    # production curve. The uncalibrated branch below solves it properly.
    neff_arm_mpb = cal_arm.seed_beta[:, 0] * wl / (2 * np.pi)
    _excess = neff_arm / neff_arm_mpb - 1.0
    _i = M // 2
    print(f"Arm beta calibration ({ARM_CAL_SPAN:.1f} um span, res "
          f"{ARM_RESOLUTION:g}, {(time.time()-_t0)/60:.1f} min):")
    print(f"    MPB mode index  n_eff {neff_arm_mpb[_i]:.6f}")
    print(f"    FDTD own phase  n_eff {neff_arm[_i]:.6f}   <- the arms propagate on this")
    print(f"    excess {100*_excess[_i]:+.3f}%  (band range "
          f"{100*_excess.min():+.3f}% .. {100*_excess.max():+.3f}%)")
    # What that excess is worth on the quantity this device reports. The
    # LEVER ARM IS 4 * L_ARM, not L_ARM: the resonant denominator is
    # 1 - rho^2 t^2 with rho = r * phase^2, so the arm index enters the round
    # trip through phase^4. Getting this wrong is easy and was measured
    # explicitly on the ring (job 11983058) after an L_OUTER-vs-net-arm slip.
    _dphi = 2 * np.pi * (neff_arm[_i] - neff_arm_mpb[_i]) * (4 * L_ARM) / wl[_i]
    print(f"    -> round-trip phase shift {_dphi:+.4f} rad "
          f"= {_dphi/(2*np.pi):+.2f} FSR (lever arm 4 x L_ARM = {4*L_ARM:.2f} um)")
    print(f"    NOTE: this makes the loop self-consistent, NOT convergent -- "
          f"absolute comb position stays uncertain at this order (see docstring).")
else:
    # The escape hatch pays for its own MPB sweep, at every wavelength rather
    # than the calibration's 3 anchors. Solved on a SHORT window at the
    # midpoint rather than the full span: section_mode_neff() extends its cut
    # by dpml on each side, and the cross-section is identical anywhere in a
    # uniform region, so a narrow window is the same modes, cheaper.
    _t0 = time.time()
    neff_arm_mpb = arm_solver.section_mode_neff(
        arm_geom, wl, -1.0, +1.0, 1, save_dir=save_dir, cache_only=CACHE_ONLY,
    )[:, 0]
    neff_arm = neff_arm_mpb
    print(f"CALIBRATE_ARM=0 -- the arms stay on MPB's mode index, and the "
          f"round trip mixes engines. {N_WL} pts in {(time.time()-_t0)/60:.1f} min, "
          f"{neff_arm_mpb[0]:.4f} to {neff_arm_mpb[-1]:.4f}")

# ---------------------------------------------------------------------------
# 4. Build the network blocks.
#
#    The HCG block is a genuine 1-port self-reflecting block, which
#    from_smatrix() cannot represent: it zero-fills a port's own self-term by
#    design, and this IS a port reflecting into itself (A3).
#
#    The waveguide block is genuinely BIDIRECTIONAL, built directly rather
#    than via from_transfer_matrix() -- see the module docstring for the real
#    bug that caused.
# ---------------------------------------------------------------------------
# The mirror's grid must be an exact superset of the assembly grid, or slicing
# would silently compare the wrong points. Since 2026-08-24 the two ARE the same
# 36-point grid, so this reduces to an identity check -- kept because it is the
# guard that makes re-splitting them later safe rather than silently wrong.
idx = [int(np.argmin(np.abs(wl_hcg - w))) for w in wl]
if not np.allclose(wl_hcg[idx], wl, atol=1e-9):
    raise ValueError(
        "The HCG's wavelength grid is not an exact superset of the coupler's -- "
        "slicing would silently compare the wrong points."
    )
r = r_hcg[idx]

hcg_block = SMatrixBlock(
    wavelengths=wl, ports=("p",), modes_per_port={"p": 1}, data=r[:, None, None])

phase = np.exp(2j * np.pi * neff_arm / wl * L_ARM)
wg_data = np.zeros((M, 2, 2), dtype=complex)
wg_data[:, 1, 0] = phase
wg_data[:, 0, 1] = phase
wg_block_left = SMatrixBlock(wavelengths=wl, ports=("in", "out"),
                             modes_per_port={"in": 1, "out": 1}, data=wg_data)
wg_block_right = SMatrixBlock(wavelengths=wl, ports=("in", "out"),
                              modes_per_port={"in": 1, "out": 1}, data=wg_data.copy())

ports = ("through_left", "through_right", "input_reflection", "output")
L_, R_, I_, O_ = 0, 1, 2, 3
dc_data = np.zeros((M, 4, 4), dtype=complex)
dc_data[:, L_, R_] = t; dc_data[:, R_, L_] = t
dc_data[:, I_, O_] = t; dc_data[:, O_, I_] = t
dc_data[:, L_, O_] = k; dc_data[:, O_, L_] = k
dc_data[:, R_, I_] = k; dc_data[:, I_, R_] = k
dc_block = SMatrixBlock(wavelengths=wl, ports=ports,
                        modes_per_port={p: 1 for p in ports}, data=dc_data)

# ---------------------------------------------------------------------------
# 5. Assemble -- the real Circuit wiring, step 4.
# ---------------------------------------------------------------------------
circuit = Circuit()
circuit.add_block("hcg_left", hcg_block)
circuit.add_block("wg_left", wg_block_left)
circuit.add_block("dc", dc_block)
circuit.add_block("wg_right", wg_block_right)
circuit.add_block("hcg_right", hcg_block)

circuit.connect("hcg_left", "p", "wg_left", "in")
circuit.connect("wg_left", "out", "dc", "through_left")
circuit.connect("dc", "through_right", "wg_right", "in")
circuit.connect("wg_right", "out", "hcg_right", "p")

cavity = circuit.reduce(external_ports={
    "dc.input_reflection": "input_reflection",
    "dc.output": "output",
})
R_circuit = np.abs(cavity.submatrix("input_reflection", "input_reflection")[:, 0, 0]) ** 2
T_circuit = np.abs(cavity.submatrix("output", "input_reflection")[:, 0, 0]) ** 2

# ---------------------------------------------------------------------------
# 6. Independent cross-check -- the hand-derived closed form (step 5),
#    verified separately against a brute-force 4-port solve in
#    tests/test_fabry_perot_network_crosscheck.py. Must agree with the
#    Circuit result to floating-point precision, or something in the
#    wiring/port-labeling above is wrong.
# ---------------------------------------------------------------------------
rho = r * phase ** 2
D = 1 - rho ** 2 * t ** 2
self_reflection = k ** 2 * rho / D
transmission = t * (1 + rho ** 2 * (k ** 2 - t ** 2)) / D
R_closed = np.abs(self_reflection) ** 2
T_closed = np.abs(transmission) ** 2

max_R_diff = float(np.max(np.abs(R_circuit - R_closed)))
max_T_diff = float(np.max(np.abs(T_circuit - T_closed)))
print(f"\nCircuit vs closed-form cross-check: max|dR|={max_R_diff:.2e}  max|dT|={max_T_diff:.2e}")
if max_R_diff > 1e-8 or max_T_diff > 1e-8:
    raise AssertionError(
        f"Circuit assembly and the closed form disagree by more than "
        f"floating-point precision (max|dR|={max_R_diff:.2e}, "
        f"max|dT|={max_T_diff:.2e}) -- a wiring or port-labeling bug, not a "
        f"physical result."
    )
print("Circuit and closed-form agree to floating-point precision.")

# ---------------------------------------------------------------------------
# 6b. Regression guard against the two pieces' own published results.
#
#     Best-effort, not required: if the constituent scripts' pickles are
#     still on disk, the rebuild above must reproduce them to machine
#     precision. That is what makes "rebuilt in place" a refactor rather
#     than a new device. A mismatch means this file's constants have drifted
#     from the script it mirrors -- which is exactly the drift that made
#     loading the pickles look safe in the first place.
# ---------------------------------------------------------------------------
def _compare_published(path, label, grid, pairs):
    if not os.path.exists(path):
        print(f"  {label}: no published result at {path} -- skipped (not required).")
        return
    ref = SMatrix.load(path)
    if len(ref.wavelengths) != len(grid) or not np.allclose(ref.wavelengths, grid, atol=1e-9):
        print(f"  {label}: published grid differs from this run's -- skipped.")
        return
    for port, in_mode, mine in pairs:
        theirs = ref.coeff(port=port, in_mode=in_mode, out_mode=1)
        d = float(np.max(np.abs(mine - theirs)))
        flag = "OK" if d < 1e-12 else "MISMATCH"
        print(f"  {label} '{port}': max|rebuilt - published| = {d:.3e}  [{flag}]")
        if d >= 1e-12:
            raise AssertionError(
                f"{label} port '{port}' rebuilt here differs from the published "
                f"{path} by {d:.3e}. This script's constants have drifted from "
                f"the script that produced that file."
            )


print("\nRegression guard vs the constituent scripts' published results:")
_compare_published(
    os.path.join(_ROOT, "example_output", "high_contrast_grating",
                 "high_contrast_grating_fdtd", "smatrix.pkl"),
    "HCG mirror", wl_hcg, [("in", 1, r_hcg)],
)
_compare_published(
    os.path.join(_ROOT, "example_output", "directional_coupler",
                 "directional_coupler_soi160_hybrid", "smatrix_hybrid.pkl"),
    "coupler", wl, [("through", 1, t), ("cross", 1, k)],
)

# ---------------------------------------------------------------------------
# 7. Report, save, plot -- the sparse grid.
# ---------------------------------------------------------------------------
print(f"\n{'lambda':>8} {'R (input_reflection)':>21} {'T (output)':>11}")
for i, w in enumerate(wl):
    print(f"{w:8.3f} {R_circuit[i]:21.6f} {T_circuit[i]:11.6f}")
i_peak = int(np.argmax(T_circuit))
print(f"\nPeak output transmission: {T_circuit[i_peak]:.4f} at {wl[i_peak]:.3f} um")
print("Passivity check (R+T <= 1):", float(np.max(R_circuit + T_circuit)))
print("\nThis grid (10 nm spacing) is far coarser than this cavity's own "
      "estimated ~0.41 nm resonance FWHM -- see the dense reconstruction "
      "below (section 8) for the actual resolved lineshape.")

def _body_beta_from_diagonal(diag, wl_arr, length, neff_anchor, min_margin=5.0):
    """The body's beta recovered from the EME transfer matrix's own diagonal.

    WHY THIS EXISTS: with CALIBRATE_BODY=0 there is no `cal_body`, so until
    2026-08-25 this script saved `beta_body` as ZEROS -- the same key name, the
    same shape and dtype as the calibrated run's real measurement, so a reader
    could not tell the placeholder from data. The beta was never lost (the run
    computed and propagated it), it was just never written down. This recovers
    it from what WAS propagated.

    M_body is diagonal and single-cell, so M[m, m] = exp(1j*beta_m*length) and
    beta follows from the phase -- except that beta*length is ~22 rad here, so
    the phase has wrapped ~3.5 times and angle() alone is 2*pi*k short of the
    answer. Unwrapping across wavelength fixes the SHAPE (measured: adjacent
    steps are <= 0.03 rad, far inside pi, so the unwrap is unambiguous) but not
    the offset: the guided window 1.4441 < n_eff < 3.4779 admits four candidate
    k, and all four have normal dispersion, so neither bound nor slope picks one.

    `neff_anchor` breaks the tie on physics this device already measured for
    free: the coupled pair's even/odd supermodes STRADDLE the isolated strip's
    own index, within ~0.05 of it, so the arm's MPB n_eff -- already solved in
    the CALIBRATE_ARM=0 branch -- identifies k with a wide margin (measured
    2026-08-25: the correct k sits 0.013 away, the runners-up 0.61 and 0.63, a
    factor of ~45). The margin is CHECKED, not assumed: if the best k is not
    `min_margin` times closer than the next, this returns NaN rather than a
    confidently wrong beta.

    Returns (beta [n_wl, n_modes], k per mode, source string).
    """
    n_wl, n_modes = diag.shape
    beta = np.full((n_wl, n_modes), np.nan)
    ks = []
    anchor = float(np.mean(neff_anchor))
    for m in range(n_modes):
        unwrapped = np.unwrap(np.angle(diag[:, m]))
        scored = []
        for k in range(0, 12):
            n_eff = (unwrapped + 2 * np.pi * k) * wl_arr / (2 * np.pi * length)
            if n_eff.min() <= 1.4441 or n_eff.max() >= 3.4779:
                continue          # not a guided mode of this Si/SiO2 stack
            scored.append((abs(float(np.mean(n_eff)) - anchor), k, n_eff))
        scored.sort()
        if not scored:
            print(f"    WARNING: mode {m} -- no candidate branch lands in the "
                  f"guided window; beta_body left as NaN")
            ks.append(-1)
            continue
        best = scored[0]
        if len(scored) > 1 and best[0] * min_margin > scored[1][0]:
            print(f"    WARNING: mode {m} -- branch ambiguous (best k={best[1]} at "
                  f"d={best[0]:.4f}, next k={scored[1][1]} at d={scored[1][0]:.4f}, "
                  f"margin < {min_margin:g}x); beta_body left as NaN")
            ks.append(-1)
            continue
        beta[:, m] = 2 * np.pi * best[2] / wl_arr
        ks.append(best[1])
    return beta, ks, "EME diagonal phase, 2pi branch anchored on the arm MPB index"


# beta_body is now written for BOTH branches -- measured directly when the body
# is calibrated, recovered from the propagated diagonal when it is not.
_M_body_diag = np.diagonal(M_body, axis1=1, axis2=2).copy()   # [n_wl, N_MODES]
if cal_body is not None:
    beta_body_saved = cal_body.beta
    beta_body_source = "FDTDSolver.section_mode_beta_details (calibrated)"
    beta_body_k = []
else:
    print("Recovering the uncalibrated body's beta from the EME diagonal "
          "(CALIBRATE_BODY=0, so there is no calibration object to read):")
    beta_body_saved, beta_body_k, beta_body_source = _body_beta_from_diagonal(
        _M_body_diag, wl, body_length, neff_arm_mpb)
    if np.all(np.isfinite(beta_body_saved)):
        _nb = beta_body_saved * wl[:, None] / (2 * np.pi)
        print(f"    branch k per mode: {beta_body_k}")
        print(f"    n_eff at {wl[M//2]:.4f} um: {_nb[M//2, 0]:.6f} {_nb[M//2, 1]:.6f}"
              f"   delta_n_eff {_nb[M//2, 0] - _nb[M//2, 1]:.6f}")

if MASTER:
    smat = SMatrix(wavelengths=wl)
    smat.set(self_reflection, port="input_reflection")
    smat.set(transmission, port="output")
    smat.save(os.path.join(output_dir, "smatrix_fabry_perot.pkl"))

    # The intermediate pieces, so a downstream analysis never has to re-derive
    # them (and so a future reader can see exactly which n_eff the arms used).
    np.savez(
        os.path.join(output_dir, "cavity_pieces.npz"),
        wavelengths=wl, wavelengths_hcg=wl_hcg,
        r_hcg=r_hcg, M_in0=M_in0, M_in1=M_in1,
        beta_body=beta_body_saved,
        # The exact complex diagonal that was PROPAGATED, saved unconditionally
        # and losslessly: beta_body above is derived from it in the uncalibrated
        # branch, and this is the primitive that derivation cannot corrupt.
        M_body_diag=_M_body_diag,
        body_length=body_length, coupling_length=COUPLING_LENGTH,
        neff_arm=neff_arm, neff_arm_mpb=neff_arm_mpb, L_arm=L_ARM,
        calibrate_arm=int(CALIBRATE_ARM), arm_resolution=ARM_RESOLUTION,
        # Self-describing provenance: without these, three folders hold a
        # `beta_body` of identical shape and dtype meaning different things, and
        # only the folder NAME distinguishes them -- which is not in the data.
        calibrate_body=int(CALIBRATE_BODY),
        beta_body_source=str(beta_body_source),
        eme_resolution=EME_RESOLUTION,
    )

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(wl, R_circuit, "o-", label="Reflection (input_reflection)", color="darkorange")
    ax.plot(wl, T_circuit, "o-", label="Transmission (output)", color="steelblue")
    ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Power")
    ax.set_title("HCG Fabry-Perot resonator -- every piece solved in this script\n"
                 "(10 nm grid -- known under-sampled near resonance)")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "spectrum.png"), dpi=300)
    plt.close("all")
    print(f"\nSMatrix saved -> {output_dir}/smatrix_fabry_perot.pkl")
    print(f"Intermediate pieces saved -> {output_dir}/cavity_pieces.npz")
    print(f"Spectrum plot saved -> {output_dir}/spectrum.png")

# ---------------------------------------------------------------------------
# 8. Dense resonance characterization -- per-factor resampling, no new FDTD.
#    See the module docstring for the full derivation and which factor gets
#    which treatment. Every resampled factor is CHECKED (unwrap_is_safe /
#    check_unwrap), never assumed safe.
# ---------------------------------------------------------------------------
print(f"\n--- 8. Dense resonance characterization ({N_WL_DENSE} pts) ---")
wl_dense = np.linspace(WL_MIN, WL_MAX, N_WL_DENSE)
print(f"Dense grid: {N_WL_DENSE} pts, {wl_dense[0]:.4f}-{wl_dense[-1]:.4f} um "
      f"({(wl_dense[-1]-wl_dense[0])/(N_WL_DENSE-1)*1000:.4f} nm spacing)")

print("Per-factor unwrap safety (margin=0.5, i.e. worst step < 90 deg):")
for name, arr in [("M_in0", M_in0[:, :, 0]), ("M_in1", M_in1[:, :, 0])]:
    for m in range(N_MODES):
        safe, worst = unwrap_is_safe(arr[:, m], margin=0.5)
        print(f"  {name} mode {m}: safe={safe}  worst_step={np.degrees(worst):.1f} deg")
        if not safe:
            raise ValueError(
                f"{name} mode {m} fails unwrap_is_safe at margin=0.5 "
                f"(worst_step={np.degrees(worst):.1f} deg) -- the per-factor "
                f"interpolation this section relies on is NOT safe for this "
                f"device/grid. Raise N_WL (a real re-solve) rather than "
                f"trusting this interpolation."
            )
safe_r, worst_r = unwrap_is_safe(r_hcg, margin=0.5)
print(f"  r (HCG, 21-pt grid): safe={safe_r}  worst_step={np.degrees(worst_r):.1f} deg")
if not safe_r:
    raise ValueError(
        f"The HCG's r(lambda) fails unwrap_is_safe at margin=0.5 on its own "
        f"{N_WL_HCG}-point grid (worst_step={np.degrees(worst_r):.1f} deg) -- "
        f"raise N_WL_HCG rather than interpolating it."
    )


def _fit_smooth(x, y, label, deg=2, tol=0.01):
    """degree-`deg` polynomial fit + a residual gate, in the units the caller
    cares about: the residual as a fraction of the quantity's OWN band-wide
    range. A fit that eats 1% of the range is not good enough to evaluate a
    phase from over a 212 um lever arm."""
    p = np.polyfit(x, y, deg)
    resid = y - np.polyval(p, x)
    max_resid = float(np.max(np.abs(resid)))
    rel = max_resid / float(np.ptp(y))
    print(f"  {label}: degree-{deg} fit max residual {max_resid:.2e} "
          f"({rel:.2%} of its own band-wide range)")
    if rel > tol:
        raise ValueError(
            f"{label}'s degree-{deg} fit residual is {rel:.2%} of its own "
            f"range -- too large to trust the analytic phase evaluation. Use "
            f"more sparse points or a higher fit degree before trusting this."
        )
    return p


print("Smooth factors, fitted analytically (never interpolated):")
neff_arm_fit = _fit_smooth(wl, neff_arm, "n_eff (connecting waveguide)")
beta_body_fit = None
if CALIBRATE_BODY:
    beta_body_fit = [
        _fit_smooth(wl, cal_body.beta[:, m], f"beta_body mode {m}")
        for m in range(N_MODES)
    ]

# --- resample the safe factors -------------------------------------------
r_dense = interp_complex(wl_hcg, r_hcg, wl_dense, check_unwrap=0.5)
M_in0_dense = np.stack(
    [interp_complex(wl, M_in0[:, m, 0], wl_dense, check_unwrap=0.5) for m in range(N_MODES)],
    axis=1)[:, :, None]
M_in1_dense = np.stack(
    [interp_complex(wl, M_in1[:, m, 0], wl_dense, check_unwrap=0.5) for m in range(N_MODES)],
    axis=1)[:, :, None]

# --- evaluate the smooth factors exactly ---------------------------------
M_body_dense = np.zeros((len(wl_dense), N_MODES, N_MODES), dtype=complex)
if CALIBRATE_BODY:
    for m in range(N_MODES):
        M_body_dense[:, m, m] = np.exp(1j * np.polyval(beta_body_fit[m], wl_dense) * body_length)
else:
    # Uncalibrated body: the EME's own diagonal, resampled the same way the
    # arm-excitation factors are (it is a phase factor of the same character).
    for m in range(N_MODES):
        M_body_dense[:, m, m] = interp_complex(
            wl, np.diagonal(M_body, axis1=1, axis2=2)[:, m], wl_dense, check_unwrap=0.5)

neff_arm_dense = np.polyval(neff_arm_fit, wl_dense)
phase_dense = np.exp(2j * np.pi * neff_arm_dense / wl_dense * L_ARM)

# --- recombine, with the SAME contraction section 2 used ------------------
_tmp_d = np.einsum("wij,wjk->wik", M_body_dense, M_in0_dense)
t_dense = np.einsum("wij,wjk->wik", M_in0_dense.transpose(0, 2, 1), _tmp_d)[:, 0, 0]
k_dense = np.einsum("wij,wjk->wik", M_in1_dense.transpose(0, 2, 1), _tmp_d)[:, 0, 0]
print(f"Coupler, dense: max(|t|^2+|k|^2) = {np.max(np.abs(t_dense)**2 + np.abs(k_dense)**2):.4f}")

# Sanity: the dense reconstruction must land back on the sparse points it was
# built from. This is the one place the whole per-factor scheme can be checked
# against data it did not fit -- the recombination is not part of any fit.
_at = [int(np.argmin(np.abs(wl_dense - w))) for w in wl]
_dt = float(np.max(np.abs(t_dense[_at] - t)))
_dk = float(np.max(np.abs(k_dense[_at] - k)))
print(f"Dense-vs-sparse coupler agreement at the {N_WL} measured points: "
      f"max|dt|={_dt:.2e} max|dk|={_dk:.2e}")

rho_dense = r_dense * phase_dense ** 2
D_dense = 1 - rho_dense ** 2 * t_dense ** 2
self_reflection_dense = k_dense ** 2 * rho_dense / D_dense
transmission_dense = t_dense * (1 + rho_dense ** 2 * (k_dense ** 2 - t_dense ** 2)) / D_dense
R_dense = np.abs(self_reflection_dense) ** 2
T_dense = np.abs(transmission_dense) ** 2
print(f"Passivity check, dense grid (max R+T, should be <=1): {np.max(R_dense + T_dense):.4f}")

# This dense evaluation is the closed form ONLY -- NOT a re-run of Circuit at
# thousands of points. The two already agree to floating-point precision on
# the sparse grid (section 6); re-verifying that same equivalence at every
# dense point would be redundant, not new information.


def _find_dips(wl_grid, T):
    """Numpy-only local-minima + FWHM/Q extraction -- same convention as
    ring_resonator_hybrid.py's own find_dips()."""
    dips = []
    baseline = np.median(T)
    for i in range(1, len(T) - 1):
        if T[i] < T[i - 1] and T[i] < T[i + 1] and T[i] < 0.7 * baseline:
            half = (T[i] + baseline) / 2.0
            lo = i
            while lo > 0 and T[lo] < half:
                lo -= 1
            hi = i
            while hi < len(T) - 1 and T[hi] < half:
                hi += 1
            if hi > lo:
                fwhm = wl_grid[hi] - wl_grid[lo]
                q = wl_grid[i] / fwhm if fwhm > 0 else np.nan
                dips.append((wl_grid[i], T[i], fwhm, q))
    return dips


dips = _find_dips(wl_dense, T_dense)
print(f"\nDense-grid assembly ({N_WL_DENSE} points): {len(dips)} resonance dip(s) found.")
if dips:
    print(f"  {'lambda (um)':>12}  {'T_min':>8}  {'FWHM (nm)':>10}  {'Q':>10}")
    for wl0, Tmin, fwhm, q in dips:
        print(f"  {wl0:12.4f}  {Tmin:8.4f}  {fwhm*1000:10.4f}  {q:10.1f}")
if len(dips) >= 2:
    fsr_meas = float(np.mean(np.diff([d[0] for d in dips])))
    print(f"  Measured mean FSR: {fsr_meas*1000:.3f} nm")

if MASTER:
    smat_dense = SMatrix(wavelengths=wl_dense)
    smat_dense.set(self_reflection_dense, port="input_reflection")
    smat_dense.set(transmission_dense, port="output")
    smat_dense.save(os.path.join(output_dir, "smatrix_fabry_perot_dense.pkl"))

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(wl_dense, R_dense, label="Reflection (input_reflection)", color="darkorange", lw=0.8)
    ax.plot(wl_dense, T_dense, label="Transmission (output)", color="steelblue", lw=0.8)
    ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
    for wl0, Tmin, fwhm, q in dips:
        ax.axvline(wl0, color="k", linestyle="--", lw=0.5, alpha=0.5)
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Power")
    ax.set_title(f"HCG Fabry-Perot resonator -- {len(dips)} resonances resolved\n"
                 f"{N_WL_DENSE}-pt per-factor reconstruction; arms on "
                 f"{'FDTD-calibrated beta' if CALIBRATE_ARM else 'MPB mode index'}")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "spectrum_dense.png"), dpi=200)
    plt.close("all")
    print(f"\nSMatrix saved -> {output_dir}/smatrix_fabry_perot_dense.pkl")
    print(f"Spectrum plot saved -> {output_dir}/spectrum_dense.png")

# ---------------------------------------------------------------------------
# 9. Every constituent piece's own geometry/diagnostic plots -- generated
#    here, from the same geometry/solver objects that were simulated, rather
#    than copied out of another script's output folder. The picture cannot
#    drift from what ran.
# ---------------------------------------------------------------------------
if MASTER:
    print(f"\n--- 9. Geometry and diagnostic plots ---")

    # --- 9a. HCG mirror domain -------------------------------------------
    # plot_fdtd_geometry() assumes a geometry that extends its own silicon by
    # stub_length to reach the PML -- true for every parametric geometry,
    # false for GdsGeometry -- so it would draw the source/monitor markers in
    # the wrong place. This mirrors solve_reflection()'s own domain formula.
    def plot_hcg_domain(geometry, solver, edge_clearance, stub_length, output_path):
        Sx = geometry.total_length + 2 * (edge_clearance + solver.dpml)
        Sy = solver.sim_width + 2 * solver.dpml
        polys = geometry.polygons(stub_length=stub_length)
        fig, ax = plt.subplots(figsize=(9, 8))
        pc = PatchCollection(
            [MplPolygon(p, closed=True) for p in polys.get((2, 0), [])],
            facecolor="steelblue", edgecolor="k", linewidths=0.3, alpha=0.85)
        ax.add_collection(pc)
        src_x = -Sx / 2 + solver.dpml + solver.source_margin
        mon_x = -Sx / 2 + solver.dpml + solver.monitor_margin
        ax.axvline(src_x, color="tab:green", lw=1.5, label="source")
        ax.axvline(mon_x, color="tab:cyan", lw=1.5, label="reflection monitor")
        ax.axvline(-geometry.total_length / 2, color="tab:red", ls="--", lw=1, label="GDS bbox edge")
        ax.axvline(geometry.total_length / 2, color="tab:red", ls="--", lw=1)
        ax.axvspan(-Sx / 2, -Sx / 2 + solver.dpml, color="orange", alpha=0.15, label="PML (x)")
        ax.axvspan(Sx / 2 - solver.dpml, Sx / 2, color="orange", alpha=0.15)
        ax.axhspan(-Sy / 2, -Sy / 2 + solver.dpml, color="purple", alpha=0.15, label="PML (y)")
        ax.axhspan(Sy / 2 - solver.dpml, Sy / 2, color="purple", alpha=0.15)
        ax.set_xlim(-Sx / 2 - 0.2, Sx / 2 + 0.2)
        ax.set_ylim(-Sy / 2 - 0.2, Sy / 2 + 0.2)
        ax.set_aspect("equal")
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
        ax.set_title("High-contrast grating — FDTD domain (solve_reflection)\n"
                     "input waveguide extended through PML")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.savefig(output_path, dpi=150)
        plt.close(fig)

    hcg_dir = os.path.join(output_dir, "hcg_mirror")
    os.makedirs(hcg_dir, exist_ok=True)
    plot_hcg_domain(hcg_geom, hcg_solver, EDGE_CLEARANCE,
                    EDGE_CLEARANCE + DPML,   # == solve_reflection(extend_input_waveguide=True)
                    os.path.join(hcg_dir, "geometry.png"))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(wl_hcg, R_hcg, marker="o", ms=3, color="k", label="R (reflection into input mode)")
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Reflection")
    ax.set_title("High-contrast grating mirror — FDTD reflection vs wavelength")
    ax.legend(); ax.grid(True); ax.set_ylim(0, 1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(hcg_dir, "reflection.png"), dpi=300)
    plt.close("all")
    print(f"HCG mirror geometry/reflection -> {hcg_dir}/")

    # --- 9b. Coupler geometry + its own through/cross ---------------------
    dc_dir = os.path.join(output_dir, "directional_coupler")
    os.makedirs(dc_dir, exist_ok=True)
    plot_hybrid_geometry(
        dc,
        solver=HybridSolver(solvers=[dc_fdtd, dc_eme, dc_fdtd],
                            boundaries=boundaries,
                            interface_modes=[N_MODES, N_MODES]),
        in_modes=1, out_modes=1,
        title="Directional Coupler (SOI_160, duty/dc.gds dims) — Hybrid Solver Sections",
        output_path=os.path.join(dc_dir, "geometry_hybrid.png"),
        sub_region_dir=os.path.join(dc_dir, "sections"),
        sub_region_show_pml=True,
        sub_region_show_eme_detail=True,
        sub_region_configs=[
            ("in", {0: {"input_idx": 0, **_iface_mon}}, {2: 0}),
            ("isolated", {0: {"input_idx": 1, **_iface_mon}}, {2: 0}),
        ],
    )
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(wl, np.abs(t) ** 2, "o-", label="through |t|^2")
    ax.plot(wl, np.abs(k) ** 2, "o-", label="cross |k|^2")
    ax.plot(wl, np.abs(t) ** 2 + np.abs(k) ** 2, "k:", label="total")
    ax.set_xlabel("Wavelength (um)"); ax.set_ylabel("Power")
    ax.set_title(f"Directional coupler (SOI_160, L_c={COUPLING_LENGTH:.1f} um)"
                 f"{' — beta-calibrated body' if CALIBRATE_BODY else ''}")
    ax.legend(fontsize=8); ax.grid(True); ax.set_ylim(0, 1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(dc_dir, "hybrid_results.png"), dpi=200)
    plt.close("all")
    print(f"Coupler geometry/sections/results -> {dc_dir}/")

    # --- 9c. The two calibrations, plotted -------------------------------
    # Free: no extra simulation, just the data section_mode_beta_details()
    # already produced. The residual guard reduces a whole calibration to one
    # number per cell and only raises past 0.1 rad, so a marginal fit passes
    # silently and would otherwise be invisible in the n_eff that was used.
    if cal_body is not None:
        plot_beta_calibration(cal_body, dc, dc_fdtd,
                              os.path.join(dc_dir, "calibration"), **_iface_mon)
        print(f"Coupler body calibration plots -> {dc_dir}/calibration/")
    if cal_arm is not None:
        arm_dir = os.path.join(output_dir, "connecting_waveguide")
        plot_beta_calibration(cal_arm, arm_geom, arm_solver,
                              os.path.join(arm_dir, "calibration"))
        fig, ax = plt.subplots(figsize=(6, 4))
        wl_fine = np.linspace(WL_MIN, WL_MAX, 200)
        ax.plot(wl, neff_arm_mpb, "s", color="0.55", mfc="none",
                label="MPB mode index (not propagated on)")
        ax.plot(wl, neff_arm, "o", label="FDTD-calibrated (propagated on)")
        ax.plot(wl_fine, np.polyval(neff_arm_fit, wl_fine), "-", label="degree-2 fit")
        ax.set_xlabel("Wavelength (um)"); ax.set_ylabel("n_eff")
        ax.set_title(f"Connecting waveguide n_eff\nFDTD-calibrated, res {ARM_RESOLUTION:g}")
        ax.legend(fontsize=8); ax.grid(True)
        plt.tight_layout()
        fig.savefig(os.path.join(arm_dir, "neff.png"), dpi=200)
        plt.close("all")
        print(f"Arm calibration plots -> {arm_dir}/")

    # --- 9d. Overall assembled-cavity schematic ---------------------------
    # NOT a real shared FDTD domain -- each piece above is independently
    # solved in its own coordinate frame (A2). A labeled block diagram to
    # x-scale, built from each piece's own measured real length.
    def _draw_block(ax, x0, width, label, color, height=1.0, y0=-0.5):
        ax.add_patch(plt.Rectangle((x0, y0), width, height, facecolor=color, edgecolor="k", lw=1.0))
        ax.annotate(label, (x0 + width / 2, y0 + height + 0.15), ha="center", fontsize=9)
        ax.annotate(f"{width:.2f} um", (x0 + width / 2, y0 - 0.3), ha="center", fontsize=8, color="gray")
        return x0 + width

    _arm_label = ("waveguide\n(analytic, FDTD-calibrated)" if CALIBRATE_ARM
                  else "waveguide\n(analytic, MPB index)")
    fig, ax = plt.subplots(figsize=(14, 3.5))
    x1 = _draw_block(ax, 0.0, HCG_LENGTH, "HCG mirror\n(FDTD)", "#4C72B0")
    x2 = _draw_block(ax, x1, L_ARM, _arm_label, "#DD8452", height=0.4, y0=-0.2)
    x3 = _draw_block(ax, x2, DC_LENGTH,
                     "directional coupler\n(FDTD + calibrated EME)" if CALIBRATE_BODY
                     else "directional coupler\n(FDTD + EME)", "#55A868")
    x4 = _draw_block(ax, x3, L_ARM, _arm_label, "#DD8452", height=0.4, y0=-0.2)
    x5 = _draw_block(ax, x4, HCG_LENGTH, "HCG mirror\n(FDTD)", "#4C72B0")

    ax.annotate("", xy=(x2, -1.4), xytext=(x2, -0.5), arrowprops=dict(arrowstyle="-|>", color="k"))
    ax.annotate("input_reflection", (x2, -1.6), ha="center", fontsize=8)
    ax.annotate("", xy=(x3, -1.4), xytext=(x3, -0.5), arrowprops=dict(arrowstyle="-|>", color="k"))
    ax.annotate("output", (x3, -1.6), ha="center", fontsize=8)

    ax.set_xlim(-2, x5 + 2)
    ax.set_ylim(-2.0, 1.8)
    ax.set_xlabel("Position along cavity (um)")
    ax.set_yticks([])
    ax.set_title(f"HCG Fabry-Perot resonator -- overall assembled geometry "
                 f"(schematic, to x-scale, NOT one shared FDTD domain)\n"
                 f"total cavity length {x5:.2f} um")
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "overall_geometry.png"), dpi=200)
    plt.close("all")
    print(f"Overall assembled-cavity schematic -> {output_dir}/overall_geometry.png")
