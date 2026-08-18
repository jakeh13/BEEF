"""
Tutorial: racetrack ring resonator — FDTD/EME hybrid simulation (BEEF method),
built with raw meep/meow instead of BEEF's solver classes.

Geometry (220 nm SOI, C-band):

                          L_FAR (not drawn — analytic neff phase only)
                    .------------------------.
                   /                          \\
              far bend                     far bend      <- standalone
                  |                              |           bend FDTD,
              L_SIDE (analytic neff)          L_SIDE          reused x2
                  |                              |
             near-left bend                near-right bend  <- THIS script's
                  |                              |             coupler FDTD
    "bus in" ----------------- L_C -------------------- "bus out"
                     (near straight, coupled to bus)

Why this device needs raw meep, not FDTDSolver/HybridSolver:
    FDTDSolver only ever builds x-normal (y-z plane) sources/monitors
    (`mp.Vector3(0, w, h)` everywhere in beef/solvers/fdtd.py). The two ring
    ports on top of the near-left/near-right bends are y-normal (x-z) planes
    — nothing in BEEF draws or measures those today, and there's no
    rotation helper anywhere in beef/geometry/_compose.py (translation
    only). So this example hand-builds the two FDTD simulations directly
    with meep, the same way examples/directional_coupler/directional_coupler_mode.py
    hand-builds a meow cross-section without going through EMESolver. See
    notes.md's Geometry/Solver recommendations for what a first-class
    RingResonator DeviceGeometry + y-normal FDTDSolver support would need.

Why this is a feedback-loop sum, not a HybridSolver matrix chain:
    Every other *_hybrid.py example threads light through a straight
    cascade of sections once. A ring is a closed loop: light re-enters the
    coupling region every round trip. That's the closed-form all-pass
    resonator sum below, not a product of transfer matrices.

Simulation count — 4 FDTD runs + 13 mode solves total:
    (a) Coupler-section FDTD, bus-excited  -> t (bus->bus), k1 (bus->ring)
    (b) Coupler-section FDTD, ring-excited -> t_ring (ring->ring), k2 (ring->bus)
    (c) Standalone 90 degree bend FDTD     -> b(lambda), reused for BOTH far
        corners (they're geometrically identical to each other and far
        enough from the bus to not need their own coupled-section run)
    Mode solve on the bare single-mode strip (meow direct, no BEEF solver
    class — same pattern as directional_coupler_mode.py) at 13 wavelengths
    -> neff(lambda), used for the 3 straight-propagation phases (2 side
    straights + 1 far straight) AND to reference-plane-correct b(lambda).
    (d) Bare-strip beta calibration, CALIBRATE_BODY (default on, 2026-08-18)
        -> the n_eff those straight phases are ACTUALLY propagated on. The
        meow/MPB curve above is still solved, but only as the reference the
        calibrated one is reported against; meep's own accumulated phase is
        what the assembly uses, because every other term in the round trip
        (b, t, k1, k2, t_ring) is a meep measurement. See the CALIBRATE_BODY
        block in section 1. One merged multi-monitor FDTD run, 1 mode.

Wavelength strategy:
    FDTD DFT monitors sample N_WL_FDTD points (moderate — dense DFT across
    5 monitors x hundreds of frequencies is expensive); the mode solve runs
    at N_WL_SOLVE=13 points only (each is a real MPB solve). Both are
    interpolated onto a dense N_WL_DENSE-point grid before assembling the
    resonance: FDTD coefficients via magnitude + unwrapped-phase
    interpolation (linear real/imag interpolation would badly distort the
    fast-winding propagation phase), neff via a degree-2 polynomial fit
    (piecewise-linear neff would put an artificial kink in dphi/dlambda at
    every solve node, corrupting resonance positions). See notes.md's new
    Solver recommendation for making this a proper library-level API.

Assembly (closed-form all-pass resonator, NOT a HybridSolver cascade):
    phi = 2*pi*neff(lambda)*L_outer/lambda        # L_outer = 2 side straights (short
                                                    # by the ring-port reference offset)
                                                    # + 1 far straight
    b_corrected = b(lambda) * (per-bend reference-plane phase correction)
    P = b_corrected**2 * exp(1j*phi)               # one full round trip minus the
                                                    # coupler section (2 far bends + 3 straights)
    E_thru = t + k1*k2*P / (1 - t_ring*P)
    T = |E_thru|**2

Run:
    python /path/to/examples/ring_resonator/python/ring_resonator_hybrid.py
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
import matplotlib.patches as mpatches
import matplotlib.path as mpath

import meep as mp
import meow as mw

from beef import extrude_polygons, plot_beta_calibration
from beef.cache import CacheMissError, load_result, lookup, save_run
from beef.geometry import UniformCrossSection
from beef.geometry.base import polygons_to_prisms
from beef.materials import SOI_220
from beef.materials.stack import resolved_material_config
from beef.results import SMatrix
from beef.spectral import interp_complex
from beef.solvers.fdtd import (
    FDTDSolver, check_source_coverage, monitor_freq_grid, source_fwidth,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "ring_resonator", "ring_resonator_hybrid"))
# A resolution sweep (RING_FDTD_RES set, see FDTD_RESOLUTION below) writes into
# its own subfolder. Without this every task overwrites the same ring_data.npz.
#
# ADDED 2026-08-03, and it was already costing data. Only the CASCADED script
# had this guard -- and only the cascaded script read RING_FDTD_RES at all.
# Job 11655349 swept both scripts over res 20/30/40; this one ignored the
# variable, ran three times at res 40, and overwrote its own output twice. Two
# of the six results in that job never existed. Nothing production was harmed
# (the two_d redirect below kept it out of the way) but the sweep was silently
# a third as informative as it looked. Keep the two scripts symmetric.
_res_sweep = os.environ.get("RING_FDTD_RES")
if _res_sweep:
    output_dir = os.path.join(output_dir, "res_sweep", f"res_{float(_res_sweep):g}")
# RING_THREE_D writes into its own subfolder: a 2-D run's ring_data.npz is not
# the production result and must never land on top of it. The guard is the env
# var, not its value, so an explicit RING_THREE_D=1 also stays out of the way.
_three_d_env = os.environ.get("RING_THREE_D")
if _three_d_env is not None:
    output_dir = os.path.join(output_dir,
                              "three_d" if bool(int(_three_d_env)) else "two_d")
# Propagate the analytic straights on the beta meep itself accumulates rather
# than on MPB's mode index -- see the CALIBRATE_SPAN block in section 1 for
# why this device needs it and what it moves. Defined here, ahead of
# os.makedirs, because turning it off has to redirect the output.
CALIBRATE_BODY = bool(int(os.environ.get("RING_CALIBRATE_BODY", "1")))
if not CALIBRATE_BODY:
    output_dir = output_dir + "_uncalibrated"
os.makedirs(output_dir, exist_ok=True)
# Per-simulation geometry previews live here, NOT at the top level: the output
# root is for whole-device results. Mirrors ring_resonator_hybrid_cascaded.py.
sections_dir = os.path.join(output_dir, "sections")
os.makedirs(sections_dir, exist_ok=True)

# Persistent on-disk simulation cache (see beef/cache/) — same shared
# example_sim_cache/ directory every other example's FDTDSolver/
# EMESolver save_dir points at. This script hand-rolls raw meep/meow sims
# instead of going through those solver classes (see module docstring), so
# it fingerprints+stores/loads each run itself via beef.cache's lookup()/
# save_run()/load_result() directly rather than the solver-specific
# build_run_config()/build_whole_device_config() helpers (those expect a
# DeviceGeometry + solver instance this script doesn't have).
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# 0. Test-mode knob — TEST_MODE=True runs in minutes (2-D, low res, few
#    wavelengths) to smoke-test the whole pipeline (geometry, monitor
#    orientation/sign, assembly) before committing to the multi-hour full-res
#    3-D production run. Set to False for the real result.
# ---------------------------------------------------------------------------
TEST_MODE = bool(int(os.environ.get("RING_TEST_MODE", "0")))

# ---------------------------------------------------------------------------
# 1. Global simulation / geometry parameters (220 nm SOI, C-band)
# ---------------------------------------------------------------------------
MATERIAL_STACK = SOI_220
WAVEGUIDE_WIDTH = 0.45     # um, single-mode strip width
GAP = 0.25                 # um, bus<->ring edge-to-edge coupling gap
R = 5.0                    # um, bend radius (centerline)
L_C = 4.0                  # um, near (coupling) straight length — the coupling-strength knob
L_SIDE = 2.0                # um, each vertical side straight (analytic neff only)
L_FAR = L_C                 # um, far straight (analytic neff only)

# ring_resonator_validation.md gap 3 / the coupler's P3 knob-split precedent.
# ONE shared RESOLUTION used to drive both the meep FDTD and the meow/MPB mode
# solve. The coupler split its own for a reason -- "shared RESOLUTION was how a
# 7% error hid for months (dc_validation.md item 11)" -- and this device had the
# same defect, measured 2026-08-01 (job 11616610):
#
#   res    n_eff@1.55   d(round-trip phase)
#    40     2.311706      --   (shipped)
#   100     2.339523    +1.58 rad = 25.1% OF A FULL FSR
#
# The mode solve is what sets resonance position (phi = 2*pi*n_eff*L_outer/lambda),
# it converges at 100, and it is CHEAP -- the whole six-resolution sweep ran in
# 4:34. So it goes to 100. The FDTD side has never been measured at any
# resolution and 3-D FDTD cost scales roughly as res^4, so it stays at 40 rather
# than paying ~39x for an unmeasured claim. Splitting is what makes those two
# independent statements expressible at all.
# RING_FDTD_RES override added 2026-08-03 to match _cascaded.py, which has had
# it since the gap-3 sweep. Without it this script silently ignored the sweep
# variable -- see the output_dir note above for what that cost.
FDTD_RESOLUTION = float(os.environ.get(
    "RING_FDTD_RES", 20.0 if TEST_MODE else 40.0))   # pts/um -- meep; NOT converged
MODE_RESOLUTION = 20.0 if TEST_MODE else 100.0   # pts/um -- meow/MPB; gap 3 converged

# ---------------------------------------------------------------------------
# n_eff CALIBRATION -- propagate the analytic straights on the beta meep
# itself accumulates, not on MPB's mode index.
#
# This is the same correctness fix every other hybrid example in this repo
# now takes (resolved/directional_coupler_validation.md item 6; the current
# mechanism and cost model are step 9 of resolved/section_mode_beta_multi_
# monitor_validation.md). It reaches the ring LAST, on 2026-08-18, because
# this device does not go through HybridSolver -- it is raw meep/meow, so
# `calibrate_body_neff=True` had nothing to attach to and the device was
# simply never audited for it.
#
# WHY IT APPLIES HERE, and in fact applies harder than on a coupler body.
# Two consumers of neff, both of them mixing an MPB index into a quantity
# meep measured:
#
#   1. phi_bend_correction. run_bend_fdtd()'s raw coefficient is measured
#      between monitor planes that sit BEND_D_IN/BEND_D_OUT of straight
#      waveguide away from the arc ends, and the assembly de-embeds that
#      straight by dividing out exp(1j*2*pi*neff*(D_IN+D_OUT)/lambda). The
#      phase being removed was accumulated by meep's time-stepper. Removing
#      it with MPB's index subtracts the wrong number from a measured
#      quantity -- there is no self-consistency argument on the other side
#      of this one.
#   2. phi = 2*pi*neff*L_outer/lambda, the analytic straight propagation.
#      Same argument as an EME body spliced between FDTD sections: the
#      round-trip sum is a closed loop through b (meep), t/k1/k2/t_ring
#      (meep) and phi (MPB), and a loop that changes answer depending on
#      which engine carried which segment is not measuring the device.
#
# SCALE, from this device's own numbers rather than the coupler's. The
# res-40 MPB-vs-meep gap measured on the coupler is ~2.2-2.4% in beta. On
# L_OUTER = 7.00 um at 1.55 um that is ~1.5 rad of round-trip phase, i.e.
# roughly a quarter of an FSR -- the same order as gap 3's mode-resolution
# finding, which the validation doc calls "not a caveat, it is a wrong
# answer for the one quantity this device exists to report". So this is not
# a refinement; it moves the comb.
#
# The switch itself is up with output_dir (it inverts the folder suffix, the
# same escape-hatch shape as mmi_1x2_hybrid.py's MMI1X2_CALIBRATE_BODY).
#
# Longitudinal span the calibration is allowed to cut its fit lengths from.
# section_mode_beta() picks [s, 2s, 3s] with s = min(usable_span/3,
# max_safe_spacing/2) and usable_span = span - 2*dpml, so the span only has
# to be large enough not to be the binding constraint: at beta ~ 9.5 rad/um
# the physics cap is max_safe_spacing/2 = 1.32 um, which needs span >= 5.0.
# 6.0 leaves room for the band edge without making the merged domain longer
# than it has to be. The strip is uniform, so the span is a budget, not a
# structure -- nothing about the ring's real geometry is being asserted here.
CALIBRATE_SPAN = 6.0
# 2-D/3-D is its OWN knob, split out of TEST_MODE 2026-08-03 -- see the same
# comment in ring_resonator_hybrid_cascaded.py. Short version: validation check
# 4 can only afford the full device in 2-D (3-D is 167 h at res 40), so it needs
# a 2-D run of this script to compare against, and RING_TEST_MODE=1 cannot be
# that partner because it moves five knobs at once. RING_THREE_D=0 moves exactly
# one. Unset keeps the previous behaviour.
THREE_D = (not TEST_MODE) if _three_d_env is None else bool(int(_three_d_env))
WL_MIN, WL_MAX = 1.50, 1.60                    # um, C-band
N_WL_FDTD = 21 if TEST_MODE else 101           # FDTD DFT points (moderate — 5 monitors x freqs is expensive)
N_WL_SOLVE = 5 if TEST_MODE else 13            # mode-solve (neff) points — each a real MPB solve
N_WL_DENSE = 2001 if TEST_MODE else 20001     # dense interpolated grid, to resolve resonance dips.
# Raised from 501/2001 on 2026-08-03 because Q was GRID-LIMITED, not physical: at 2001
# points the 1.50-1.60 um span steps 0.05 nm, and every reported FWHM came out an exact
# multiple of it (0.15/0.20/0.25/0.30 = 3/4/5/6 steps), so Q was set by this knob rather
# than by either resolution. 20001 points = 0.005 nm, ~30 samples across the narrowest
# FWHM seen so far. Pure post-processing on interpolated coefficients -- costs no FDTD and
# is in no cache fingerprint (the FDTD grid is `wavelengths`, N_WL_FDTD, which is separate).
# CONVERGENCE VERIFIED, not assumed: rebuilding the dense assembly from the same cached
# coefficients at 2001/5001/10001/20001/40001/80001 points gives first-dip FWHM
# 0.1500/0.1400/0.1100/0.1050/0.1050/0.1050 nm -- flat from 20001 on. The old 2001-point
# value under-reported Q by ~43% (10006 vs 14294). Do not lower this to "save time": it is
# an interpolation, and the whole reported Q spectrum moves with it.
DPML = 0.5
SIM_HEIGHT = 1.5            # um, inner z-span when THREE_D
DECAY_TOL = 1e-11
SOURCE_BANDWIDTH_FACTOR = 2.0   # multiplies the full span 1/WL_MIN - 1/WL_MAX;
                                # FDTDSolver's own default. Must exceed 1 — see
                                # beef.solvers.fdtd.source_fwidth().
SOURCE_MARGIN = 0.0
MONITOR_MARGIN = 0.25
S_MON = 0.5                 # um, ring-port reference plane offset above each near bend's arc end
STUB_LEN = 1.5               # um, inner (non-PML) vertical/horizontal stub length at each port
PORT_APERTURE = 3.0           # um, mode-decomposition window width at every port (generous — see
                               # docstring: no port sits near a neighboring waveguide at its
                               # own monitor plane, so cross-contamination isn't a concern here)

_cell_sz = (SIM_HEIGHT + 2 * DPML) if THREE_D else 0.0
_mon_sz = SIM_HEIGHT if THREE_D else 0.0
_wg_height = MATERIAL_STACK.wg_height if THREE_D else None

# Perimeter = 2 horizontal straights (L_C near + L_FAR far) + 2 vertical
# straights (2*L_SIDE) + 4 quarter bends (= one full circle, 2*pi*R).
PERIMETER = L_C + L_FAR + 2 * L_SIDE + 2 * np.pi * R
NG_GUESS = 4.3   # typical group index for a 450x220 Si strip, C-band — just for the printed FSR estimate
FSR_GUESS = (1.55 ** 2) / (NG_GUESS * PERIMETER)
print(f"Racetrack perimeter: {PERIMETER:.2f} um  ->  predicted FSR ~ {FSR_GUESS*1000:.2f} nm "
      f"(~{(WL_MAX-WL_MIN)/FSR_GUESS:.1f} resonances across {WL_MIN}-{WL_MAX} um)")
if TEST_MODE:
    print("RING_TEST_MODE=1: reduced res/2-D/few wavelengths — smoke test only, not a physical result.")


# ---------------------------------------------------------------------------
# 2. Polygon helpers
# ---------------------------------------------------------------------------
def arc_polygon(cx, cy, r_inner, r_outer, theta0_deg, theta1_deg, n_pts=64):
    """Quarter-annulus wedge vertex list — same construction as
    beef/geometry/bend.py's _arc_polygon (cribbed here rather than
    importing a private function from that legacy, non-DeviceGeometry
    module)."""
    t0, t1 = np.radians(theta0_deg), np.radians(theta1_deg)
    thetas = np.linspace(t0, t1, n_pts)
    outer = [(cx + r_outer * np.cos(t), cy + r_outer * np.sin(t)) for t in thetas]
    inner = [(cx + r_inner * np.cos(t), cy + r_inner * np.sin(t)) for t in reversed(thetas)]
    return outer + inner


def rect_poly(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


# ---------------------------------------------------------------------------
# 3. Coupler-section geometry: bus straight + ring near straight + the two
#    near bends that curve the near straight's ends up into vertical stubs.
#    Bus at y=0; near straight at y_ring = WAVEGUIDE_WIDTH + GAP (center-to-
#    center pitch). Both near bends share radius R, centered directly above
#    the near straight's ends, so a wave entering a stub moving -y exits the
#    near straight moving +x (co-directional coupling with the bus) — see
#    the docstring's tangent-direction derivation.
# ---------------------------------------------------------------------------
W = WAVEGUIDE_WIDTH
Y_RING = W + GAP
STUB_X = L_C / 2 + R                    # x of both vertical stub columns
Y_BELOW = 2.0                            # inner clearance below the bus (also reused below as the
                                          # ring-port monitors' PORT_APERTURE/2 + margin clearance)
INNER_X_HALF = STUB_X + Y_BELOW          # bus extends this far past the stub columns (inner, no PML) —
                                          # sized so ring_L/ring_R's transverse PORT_APERTURE/2 window
                                          # clears the PML, same margin as Y_BELOW gives the bus monitors
Y_TOP = Y_RING + R + STUB_LEN            # inner y of the ring-port stub tops

COUPLER_SX = 2 * (INNER_X_HALF + DPML)
COUPLER_SY = (Y_TOP + DPML) - (-Y_BELOW - DPML)
COUPLER_Y_CENTER = ((Y_TOP + DPML) + (-Y_BELOW - DPML)) / 2.0


def build_coupler_polys():
    polys = [
        rect_poly(-INNER_X_HALF - DPML, -W / 2, INNER_X_HALF + DPML, W / 2),          # bus, through x-PML
        rect_poly(-L_C / 2, Y_RING - W / 2, L_C / 2, Y_RING + W / 2),                  # near straight
        arc_polygon(-L_C / 2, Y_RING + R, R - W / 2, R + W / 2, 180, 270),             # left near bend
        arc_polygon(L_C / 2, Y_RING + R, R - W / 2, R + W / 2, 270, 360),              # right near bend
        rect_poly(-STUB_X - W / 2, Y_RING + R, -STUB_X + W / 2, Y_TOP + DPML),         # left stub, through y-PML
        rect_poly(STUB_X - W / 2, Y_RING + R, STUB_X + W / 2, Y_TOP + DPML),           # right stub, through y-PML
    ]
    return polys


# ---------------------------------------------------------------------------
# 4. Standalone 90-degree bend geometry: x-in straight -> arc -> y-out
#    straight. Same radius R as the coupler section's near bends, so its
#    transmission b(lambda) is representative of the racetrack's other two
#    (far) corners, which are geometrically identical and isolated (far
#    enough from the bus/near straight that they don't need their own
#    coupled-section simulation).
# ---------------------------------------------------------------------------
# Built in a natural frame: input stub x in [-STUB_LEN, 0], arc center (0, R), output stub y in [R, R+STUB_LEN].
# X_MAX/Y_MIN widened past the bare bend radius by Y_BELOW (the same
# PORT_APERTURE/2 + margin clearance the coupler section uses below its bus
# straight, section 3) — without it, the arc's own outer radius (R+W/2) and
# the out_mon/norm_mon transverse mode-decomposition apertures
# (PORT_APERTURE/2 each) clip into the PML.
BEND_X_MIN, BEND_X_MAX = -STUB_LEN - DPML, R + Y_BELOW + DPML
BEND_Y_MIN, BEND_Y_MAX = -Y_BELOW - DPML, R + STUB_LEN + DPML
BEND_SX = BEND_X_MAX - BEND_X_MIN
BEND_SY = BEND_Y_MAX - BEND_Y_MIN
BEND_X_CENTER = (BEND_X_MAX + BEND_X_MIN) / 2.0
BEND_Y_CENTER = (BEND_Y_MAX + BEND_Y_MIN) / 2.0


def build_bend_polys():
    return [
        rect_poly(-STUB_LEN - DPML, -W / 2, 0.0, W / 2),               # input stub, through x-PML
        arc_polygon(0.0, R, R - W / 2, R + W / 2, -90, 0),              # the bend itself
        rect_poly(R - W / 2, R, R + W / 2, R + STUB_LEN + DPML),        # output stub, through y-PML
    ]


# Reference-plane straight lengths the bend FDTD's raw coefficient must be
# corrected by, to isolate the arc's own contribution (see run_bend_fdtd()).
BEND_D_IN = STUB_LEN - MONITOR_MARGIN
BEND_D_OUT = STUB_LEN - MONITOR_MARGIN


# ---------------------------------------------------------------------------
# The straight guide the standalone-bend FDTD carries INSIDE its own measured
# coefficient: BEND_D_IN before the arc and BEND_D_OUT after it, per bend -- and
# b is squared, because there are two far bends. phi_bend_correction strips that
# PHASE out of b so b_corrected refers to the arc endpoints, and L_OUTER then
# supplies the full physical straights. Net straight phase = L_OUTER.
#
# Equivalent formulation, verified 2026-08-03 to 1.2e-13 in T and to identical
# dip wavelengths: never touch b, and use L_OUTER - BEND_STRAIGHT_TOTAL as the
# analytic length instead. That is the DirectionalCoupler's approach
# (`body_lengths = coupling_lengths - 2*BODY_STUB`). The two agree exactly because
# both strip and re-add with the same neff. Keeping the strip because it is
# MANDATORY whenever the difference below would go negative, and because
# b_corrected ("the arc alone") stays meaningful across ring geometries.
#
# THE GUARD. The DC has one (`assert body_lengths[0] >= 0`); this script did not
# until 2026-08-03. The equivalent analytic length must not go negative: if the
# bend's own stubs cover more straight than the device physically has between its
# arcs, the model attributes straight to the bends that does not exist -- and it
# would do so SILENTLY. Headroom here is only 2.0 um, and in terms of the knob that
# sets it, only 0.5 um: BEND_STRAIGHT_TOTAL = 4*(STUB_LEN - MONITOR_MARGIN), so
# STUB_LEN = 2.0 is exactly the boundary and anything above it trips this assert
# (verified by exercising this block against broken configs, 2026-08-03). A smaller
# ring trips it without touching the bend at all -- L_SIDE=1.0, L_FAR=2.0 at the
# shipped STUB_LEN already fails.
#
# Note this uses the SHIPPED L_outer expression. Whether that expression is itself
# plane-consistent is a separate question -- see the round-trip ledger, which
# reports it.
# ---------------------------------------------------------------------------
L_OUTER = 2 * (L_SIDE - S_MON) + L_FAR
BEND_STRAIGHT_TOTAL = 2 * (BEND_D_IN + BEND_D_OUT)          # 2 far bends
BEND_STRAIGHT_HEADROOM = L_OUTER - BEND_STRAIGHT_TOTAL
assert BEND_STRAIGHT_HEADROOM >= 0, (
    f"the far bends' own straight stubs ({BEND_STRAIGHT_TOTAL:.2f} um total) exceed "
    f"the physical straight path between the arcs (L_OUTER = {L_OUTER:.2f} um), so the "
    f"DC-style equivalent analytic length would be negative. Reduce STUB_LEN or "
    f"MONITOR_MARGIN, or lengthen L_SIDE / L_FAR."
)


# ---------------------------------------------------------------------------
# 5. Raw-meep simulation helpers
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# M1 FIX, 2026-08-03 — the FDTD half now tracks the band.
#
# Until today both media were built with `si_medium(wl_center)` / a SCALAR, so
# the whole 1.50-1.60 um sweep ran on one index frozen at band centre while the
# meow/EME half re-evaluated its SampledMaterials at every wavelength. The two
# engines were simulating different silicon: n_Si off by +0.004338 at 1.50 um
# and -0.003913 at 1.60 um, which is +1.42% / -1.20% of an FSR in round-trip
# phase.
#
# The library fixed this on 2026-08-01 inside `FDTDSolver._media()`, so every
# example that goes through the solver inherited it with no edit. These ring
# scripts are the ONLY non-library callers of si_medium/sio2_medium in the repo
# -- they hand-roll mp.Simulation and never construct an FDTDSolver, so they
# could not inherit a solver fix. This function is `FDTDSolver._media()`
# replicated, deliberately line-for-line, rather than a second opinion about
# what "dispersive" means.
#
# NOTE the reparameterisation matters and is not optional: a naive dispersive
# meep Medium LOSES 0.0496 in n_eff at res 40, because meep subpixel-averages
# only the instantaneous eps and meep's Sellmeier form puts eps_inf = 1, giving
# no Si/SiO2 contrast for smoothing to act on. `si_medium(dispersive=True)`
# routes through `meep_medium_reparam()`, which puts the contrast back in
# eps_inf. See MaterialStack.si_medium()'s docstring for the measurement.
#
# RING_FDTD_DISPERSIVE overrides for an A/B; unset inherits the stack (True).
_disp_env = os.environ.get("RING_FDTD_DISPERSIVE")
FDTD_DISPERSIVE = None if _disp_env is None else bool(int(_disp_env))


def _media(wavelengths):
    """(si_medium, sio2_medium) for this run's band -- FDTDSolver._media()."""
    dispersive = MATERIAL_STACK.resolve_dispersive(FDTD_DISPERSIVE)
    wl = (np.asarray(wavelengths, dtype=float) if dispersive
          else float(np.mean(wavelengths)))
    return (MATERIAL_STACK.si_medium(wl, dispersive=dispersive),
            MATERIAL_STACK.sio2_medium(wl, dispersive=dispersive))


def make_sim(polys, wavelengths, cell_size, geometry_center):
    """`wavelengths` is the run's whole BAND, not a centre wavelength -- see
    _media() above for why that distinction is the entire M1 fix."""
    si, sio2 = _media(wavelengths)
    prisms = polygons_to_prisms(polys, si, wg_height=_wg_height)
    return mp.Simulation(
        resolution=FDTD_RESOLUTION,
        cell_size=mp.Vector3(cell_size[0], cell_size[1], _cell_sz),
        geometry_center=mp.Vector3(geometry_center[0], geometry_center[1], 0),
        boundary_layers=[mp.PML(DPML)],
        geometry=prisms,
        default_material=sio2,
    )


def add_x_monitor(sim, mon_freqs, x, yc, width=PORT_APERTURE):
    """`mon_freqs` is an explicit ascending frequency list from
    monitor_freq_grid(), never meep's (fcen, df, nfreq) triple — that triple
    spreads its points evenly in *frequency* over [fcen - df/2, fcen + df/2],
    so it cannot reproduce a wavelength array linear in λ, and it orders them
    ascending in frequency, i.e. descending in wavelength. This script used to
    pass `fwidth = (1/WL_MIN - 1/WL_MAX) * 0.5` as df, so a request for
    1.50–1.60 µm actually sampled 1.5254–1.5754 µm: half the span, reversed.
    Same defect as dc_validation.md item 1, in raw meep rather than through
    FDTDSolver, which is why the solver fix never reached it (item 10)."""
    return sim.add_mode_monitor(
        mon_freqs,
        mp.FluxRegion(center=mp.Vector3(x, yc), size=mp.Vector3(0, width, _mon_sz)),
    )


def add_y_monitor(sim, mon_freqs, y, xc, width=PORT_APERTURE):
    """See add_x_monitor() — same contract, y-normal plane."""
    return sim.add_mode_monitor(
        mon_freqs,
        mp.FluxRegion(center=mp.Vector3(xc, y), size=mp.Vector3(width, 0, _mon_sz)),
    )


def get_coeffs(sim, mon, band=1, forward=True):
    """alpha[0,:,0] = mode traveling in the monitor's implied positive
    direction (+x for an x-normal/zero-x-extent monitor, +y for a
    y-normal/zero-y-extent monitor — meep infers this from the FluxRegion's
    own zero-extent axis, same as every x-normal monitor in
    beef/solvers/fdtd.py, none of which pass direction= explicitly either).
    alpha[0,:,1] is the opposite-traveling mode. forward=False is needed
    exactly once in this script: run B's source launches -y (down into the
    bend), so its own norm monitor (below the source, catching the leading
    edge of that downward pulse) must read the "backward" (-y) index."""
    coeff = sim.get_eigenmode_coefficients(mon, [band], eig_parity=mp.NO_PARITY)
    return coeff.alpha[0, :, 0 if forward else 1]


def _json_polys(polys):
    """Polygon vertex lists -> plain JSON-safe floats, for a cache config's
    geometry_polygons (mirrors beef.cache.config's own float-rounding, so
    a cache hit is exact-equality-safe against JSON round-tripping)."""
    return [[[round(float(x), 9), round(float(y), 9)] for x, y in poly] for poly in polys]


def _sim_params(wavelengths):
    """Every module-level constant that affects the physical result of
    run_bend_fdtd() and every excite of run_coupler_fdtd() — the raw-meep
    equivalent of FDTDSolver's vars(solver) (see beef/cache/config.py's
    _solver_config). Changing any of these must change the cache key, or a
    stale run could be served for a different simulation.

    L_C is NOT here: run_bend_fdtd()'s build_bend_polys() has no dependence
    on L_C at all, so putting it in this shared helper re-invalidated that
    cache on every L_C edit for no physical reason. run_coupler_fdtd()
    really does depend on L_C (build_coupler_polys() widens with it, and
    geometry_center stays fixed at x=0.0, so there's no translation-
    invariance trick available the way there is in
    ring_resonator_hybrid_cascaded.py's run_left_fdtd()) — it adds "L_C" to
    its own config dict directly instead.

    The two `resolved_*` entries mirror FDTDSolver's
    `_resolved_spectral_config()`. They are derived from `wavelengths`, which
    the config already hashes, so they add nothing to correctness *today* —
    their job is to make the derivation itself part of the fingerprint. That
    is what strands the six pre-fix ring entries: they were written when the
    same `wavelengths` resolved to a reversed half-span monitor grid, and
    without this they would keep matching and be served silently for a run
    that now simulates something different. Nothing in dc_validation.md items
    1–3 touched this script's hash, which is exactly why those entries
    survived (item 10). Rounded to plain floats because _hash_config() runs
    json.dumps on this dict and numpy arrays are not JSON-serializable."""
    mon_freqs, _ = monitor_freq_grid(wavelengths)
    return {
        "resolved_monitor_freqs": [round(float(f), 9) for f in mon_freqs],
        "resolved_source_fwidth": round(
            float(source_fwidth(wavelengths, SOURCE_BANDWIDTH_FACTOR)), 9),
        "waveguide_width": WAVEGUIDE_WIDTH,
        "gap": GAP,
        "R": R,
        "resolution": FDTD_RESOLUTION,
        "three_d": THREE_D,
        "sim_height": SIM_HEIGHT,
        "dpml": DPML,
        "decay_tol": DECAY_TOL,
        "source_margin": SOURCE_MARGIN,
        "monitor_margin": MONITOR_MARGIN,
        "port_aperture": PORT_APERTURE,
        "y_below": Y_BELOW,
        "stub_len": STUB_LEN,
        "s_mon": S_MON,
        "material_stack": {"wg_height": MATERIAL_STACK.wg_height, "sz": MATERIAL_STACK.sz},
        # M3 (ring_resonator_validation.md): geometry alone is NOT a material
        # fingerprint. Without the resolved indices, neither the 2026-07-28
        # Sellmeier change nor the 2026-08-01 dispersion-default flip moves this
        # key, so a stale entry keeps matching and is served with the old
        # physics -- the same failure the resolved_* spectral fields above exist
        # to stop, on the material axis. The third argument is the RESOLVED
        # value (M4), not the stack's flag.
        #
        # M1 FIXED 2026-08-03: it used to be a hardcoded `False`, because the
        # FDTD half really did simulate frozen media -- si_medium(wl_center)
        # with a SCALAR. It now passes the band (see _media()), so this
        # resolves to True and the key turns over correctly, exactly as the M4
        # note predicted it would need to.
        **resolved_material_config(MATERIAL_STACK, np.asarray(wavelengths),
                                  MATERIAL_STACK.resolve_dispersive(FDTD_DISPERSIVE)),
    }


def get_norm_calibrated(sim, mon, forward=True):
    """Normalization amplitude for a source-adjacent monitor, magnitude
    calibrated against raw Poynting flux (mp.get_fluxes), phase taken from
    the modal coefficient.

    Empirically (via a standalone reproduction, see the debugging session
    behind this fix), meep's get_eigenmode_coefficients() undercounts the
    true launched power specifically for the "backward" (index-1) branch
    next to a y-normal EigenModeSource (needed here for run B's -y launch
    into the ring) — raw Poynting flux at the same monitor plane, cross-
    checked against a fully independent power-conservation accounting
    (flux out every port summed vs. flux at this norm plane, ~96%, physical
    given real bend loss), showed the coefficient-based reading was too
    small by the same factor that inflated |t_ring| above 1. The forward
    (index-0, x-normal-source) case (run A, the bend run) does not show
    this effect — conservation there is already correct to four decimal
    places using the coefficient alone — so this calibration is applied
    only where it's needed (run B's norm monitor), not everywhere."""
    coeff = get_coeffs(sim, mon, forward=forward)
    flux = np.abs(np.asarray(mp.get_fluxes(mon)))
    return np.sqrt(flux) * np.exp(1j * np.angle(coeff))


# ---------------------------------------------------------------------------
# 6. Coupler-section FDTD — two excitations, four ports each run.
# ---------------------------------------------------------------------------
def run_coupler_fdtd(excite, wavelengths, save_dir=None, cache_only=False):
    """excite: "bus" (run A) or "ring" (run B, source on the left stub,
    launched -y into the bend). Returns a dict of complex coefficient
    arrays [len(wavelengths)] keyed by physical quantity name.

    save_dir: persistent on-disk cache directory (see beef/cache/).
        None/"" (default) disables the cache — every call re-simulates.
    cache_only: if True, never run the simulation — a cache hit returns as
        usual, a miss raises CacheMissError instead of falling through to
        meep. Intended for notebooks that must reuse the script's
        already-populated cache rather than silently kicking off a fresh
        (possibly multi-hour) run."""
    config = None
    if save_dir:
        config = {
            "solver_type": "ring_resonator_hybrid.run_coupler_fdtd",
            "method": excite,
            "geometry_class": "RingResonatorCoupler",
            "geometry_polygons": {"2,0": _json_polys(build_coupler_polys())},
            "wavelengths": [round(float(w), 9) for w in wavelengths],
            "L_C": L_C,
            **_sim_params(wavelengths),
        }
        hit = lookup(save_dir, config)
        if hit is not None:
            print(f"  [cache] run_coupler_fdtd({excite!r}) hit -> {hit}")
            return load_result(hit)
    if cache_only:
        raise CacheMissError(
            f"run_coupler_fdtd({excite!r}): cache_only=True but no cached run "
            f"found for this config (save_dir={save_dir!r}). Run once with "
            f"cache_only=False to populate the cache."
        )

    wl_center = float(np.mean(wavelengths))
    fcen = 1.0 / wl_center
    fwidth = source_fwidth(wavelengths, SOURCE_BANDWIDTH_FACTOR)
    mon_freqs, wl_order = monitor_freq_grid(wavelengths)
    check_source_coverage(fcen, fwidth, wavelengths)

    polys = build_coupler_polys()
    sim = make_sim(polys, wavelengths, (COUPLER_SX, COUPLER_SY), (0.0, COUPLER_Y_CENTER))

    bus_l_x = -INNER_X_HALF + MONITOR_MARGIN
    bus_r_x = INNER_X_HALF - MONITOR_MARGIN
    ring_y = Y_RING + R + S_MON

    if excite == "bus":
        src_x = -INNER_X_HALF + SOURCE_MARGIN
        sim.sources = [mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=fwidth),
            center=mp.Vector3(src_x, 0.0), size=mp.Vector3(0, PORT_APERTURE, _mon_sz),
            eig_band=1, eig_kpoint=mp.Vector3(1, 0, 0), eig_match_freq=True,
        )]
        norm_mon = add_x_monitor(sim, mon_freqs, bus_l_x, 0.0)
    elif excite == "ring":
        src_y = Y_TOP - SOURCE_MARGIN
        sim.sources = [mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=fwidth),
            center=mp.Vector3(-STUB_X, src_y), size=mp.Vector3(PORT_APERTURE, 0, _mon_sz),
            eig_band=1, eig_kpoint=mp.Vector3(0, -1, 0), eig_match_freq=True,
        )]
        norm_mon = add_y_monitor(sim, mon_freqs, Y_TOP - MONITOR_MARGIN, -STUB_X)
    else:
        raise ValueError(f"unknown excite={excite!r}")

    bus_r_mon = add_x_monitor(sim, mon_freqs, bus_r_x, 0.0)
    ring_r_mon = add_y_monitor(sim, mon_freqs, ring_y, STUB_X)

    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=DECAY_TOL))

    # Every get_coeffs()/get_norm_calibrated() array below is in *monitor*
    # order (ascending frequency). Take each ratio there — numerator and
    # denominator are aligned — then reorder once onto `wavelengths` with
    # [wl_order], exactly as FDTDSolver.solve_section() does.
    if excite == "bus":
        a_norm = get_coeffs(sim, norm_mon, forward=True)
        t = (get_coeffs(sim, bus_r_mon, forward=True) / a_norm)[wl_order]
        k1 = (get_coeffs(sim, ring_r_mon, forward=True) / a_norm)[wl_order]
        sim.reset_meep()
        result = {"t": t, "k1": k1}
    else:
        # launched -y — flux-calibrated normalization, see get_norm_calibrated()
        a_norm = get_norm_calibrated(sim, norm_mon, forward=False)
        t_ring = (get_coeffs(sim, ring_r_mon, forward=True) / a_norm)[wl_order]
        k2 = (get_coeffs(sim, bus_r_mon, forward=True) / a_norm)[wl_order]
        sim.reset_meep()
        result = {"t_ring": t_ring, "k2": k2}

    if save_dir:
        folder = save_run(save_dir, config, result)
        print(f"  [cache] run_coupler_fdtd({excite!r}) saved -> {folder}")
    return result


# ---------------------------------------------------------------------------
# 7. Standalone 90-degree bend FDTD
# ---------------------------------------------------------------------------
def run_bend_fdtd(wavelengths, save_dir=None, cache_only=False):
    """save_dir: persistent on-disk cache directory (see beef/cache/).
    None/"" (default) disables the cache — every call re-simulates.
    cache_only: if True, never run the simulation — see
    run_coupler_fdtd()'s cache_only docstring for the full semantics."""
    config = None
    if save_dir:
        config = {
            "solver_type": "ring_resonator_hybrid.run_bend_fdtd",
            "method": "run_bend_fdtd",
            "geometry_class": "RingResonatorBend",
            "geometry_polygons": {"2,0": _json_polys(build_bend_polys())},
            "wavelengths": [round(float(w), 9) for w in wavelengths],
            **_sim_params(wavelengths),
        }
        hit = lookup(save_dir, config)
        if hit is not None:
            print(f"  [cache] run_bend_fdtd hit -> {hit}")
            return load_result(hit)
    if cache_only:
        raise CacheMissError(
            f"run_bend_fdtd(): cache_only=True but no cached run found for "
            f"this config (save_dir={save_dir!r}). Run once with "
            f"cache_only=False to populate the cache."
        )

    wl_center = float(np.mean(wavelengths))
    fcen = 1.0 / wl_center
    fwidth = source_fwidth(wavelengths, SOURCE_BANDWIDTH_FACTOR)
    mon_freqs, wl_order = monitor_freq_grid(wavelengths)
    check_source_coverage(fcen, fwidth, wavelengths)

    polys = build_bend_polys()
    sim = make_sim(polys, wavelengths, (BEND_SX, BEND_SY), (BEND_X_CENTER, BEND_Y_CENTER))

    src_x = -STUB_LEN + SOURCE_MARGIN
    sim.sources = [mp.EigenModeSource(
        mp.GaussianSource(fcen, fwidth=fwidth),
        center=mp.Vector3(src_x, 0.0), size=mp.Vector3(0, PORT_APERTURE, _mon_sz),
        eig_band=1, eig_kpoint=mp.Vector3(1, 0, 0), eig_match_freq=True,
    )]
    norm_mon = add_x_monitor(sim, mon_freqs, -STUB_LEN + MONITOR_MARGIN, 0.0)
    out_mon = add_y_monitor(sim, mon_freqs, R + STUB_LEN - MONITOR_MARGIN, R)

    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=DECAY_TOL))

    # Ratio in monitor order, reordered once onto `wavelengths` — see
    # run_coupler_fdtd().
    a_norm = get_coeffs(sim, norm_mon, forward=True)
    b_raw = (get_coeffs(sim, out_mon, forward=True) / a_norm)[wl_order]
    sim.reset_meep()

    if save_dir:
        folder = save_run(save_dir, config, b_raw)
        print(f"  [cache] run_bend_fdtd saved -> {folder}")
    return b_raw


# ---------------------------------------------------------------------------
# 8. Mode solve — bare single-mode strip neff(lambda), meow direct (no
#    BEEF solver class), same pattern as directional_coupler_mode.py.
# ---------------------------------------------------------------------------
MESH_WIDTH = 4.0
MESH_HEIGHT = 2.0


def solve_neff(wl, save_dir=None, cache_only=False):
    """save_dir: persistent on-disk cache directory (see beef/cache/).
    None/"" (default) disables the cache — every call re-solves.
    cache_only: if True, never run the solve — see run_coupler_fdtd()'s
    cache_only docstring for the full semantics."""
    config = None
    if save_dir:
        config = {
            "solver_type": "ring_resonator_hybrid.solve_neff",
            "method": "solve_neff",
            "geometry_class": "BareStrip",
            "geometry_polygons": {
                "2,0": _json_polys([rect_poly(0.0, -W / 2, 1.0, W / 2)]),
                "1,0": _json_polys([rect_poly(0.0, -MESH_WIDTH / 2, 1.0, MESH_WIDTH / 2)]),
            },
            "waveguide_width": WAVEGUIDE_WIDTH,
            "resolution": MODE_RESOLUTION,
            "mesh_width": MESH_WIDTH,
            "mesh_height": MESH_HEIGHT,
            "material_stack": {"wg_height": MATERIAL_STACK.wg_height, "sz": MATERIAL_STACK.sz},
            # M3/M4: same reasoning as _sim_params(). dispersive=True here and
            # NOT False as on the FDTD side, because this is the meow half --
            # it builds media from MATERIAL_STACK.extrusion_rules (SampledMaterial,
            # re-evaluated per wavelength), so it genuinely IS exact-per-lambda.
            # That asymmetry is M2, and recording it truthfully in the key is
            # the point: the two halves currently simulate DIFFERENT media, and
            # a fingerprint that hid that would be the bug, not the fix.
            **resolved_material_config(MATERIAL_STACK, np.asarray([float(wl)]), True),
            "wavelength": round(float(wl), 9),
        }
        hit = lookup(save_dir, config)
        if hit is not None:
            return load_result(hit)
    if cache_only:
        raise CacheMissError(
            f"solve_neff(wl={wl!r}): cache_only=True but no cached run found "
            f"for this config (save_dir={save_dir!r}). Run once with "
            f"cache_only=False to populate the cache."
        )

    length = 1.0
    wg_poly = rect_poly(0.0, -W / 2, length, W / 2)
    ox_poly = rect_poly(0.0, -MESH_WIDTH / 2, length, MESH_WIDTH / 2)
    gds_polys = {(2, 0): [wg_poly], (1, 0): [ox_poly]}

    structs = extrude_polygons(gds_polys, MATERIAL_STACK.extrusion_rules)
    mesh = mw.Mesh2D(
        x=np.linspace(-MESH_WIDTH / 2, MESH_WIDTH / 2, int(MODE_RESOLUTION * MESH_WIDTH) + 1),
        y=np.linspace(-MESH_HEIGHT / 2, MESH_HEIGHT / 2, int(MODE_RESOLUTION * MESH_HEIGHT) + 1),
    )
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=float(wl), T=25.0)
    cs = mw.CrossSection.from_cell(cell=cells[0], env=env)

    from meow import compute_modes_meep
    modes = compute_modes_meep(cs, num_modes=1)
    neff = float(np.real(modes[0].neff))

    if save_dir:
        save_run(save_dir, config, neff)
    return neff


# ---------------------------------------------------------------------------
# 8b. n_eff calibration -- the SAME beta measurement every other hybrid example
#     in this repo uses, reached through the library rather than re-derived
#     here. See the CALIBRATE_BODY block in section 1 for why this device
#     needs it.
#
#     WHY UniformCrossSection + FDTDSolver AND NOT ANOTHER HAND-ROLLED MEEP RUN.
#     Everything else in this file is raw meep, for the reason the module
#     docstring gives: FDTDSolver only builds x-normal sources/monitors and the
#     ring's two ring-side ports are y-normal. That constraint does not apply
#     to a straight strip. The bare cross-section is x-normal end to end, which
#     is exactly the case beef/geometry/uniform.py was written for -- its own
#     module docstring names solve_neff()'s bare strip as motivating example #1
#     and gives `UniformCrossSection([(0.0, 0.45)], length=1.0)` verbatim.
#
#     So this goes through section_mode_beta_details() and inherits, for free,
#     everything step 9 of resolved/section_mode_beta_multi_monitor_
#     validation.md built: the merged multi-monitor domain, the diagonal-only
#     decode, the cached 3-anchor MPB seed, the 2*pi branch guard, the
#     phase-linearity residual guard, and plot_beta_calibration()'s figures.
#     Re-deriving any of that against this script's own conventions would be
#     duplicating tested library code and would produce a beta with none of
#     the guards attached.
#
#     The solver is configured from THIS script's constants, not FDTDSolver's
#     defaults -- resolution, PML, margins, decay tolerance, source bandwidth,
#     3-D/2-D and the dispersion resolution all come from the module-level
#     values the ring's own runs use. That is what makes the measured beta the
#     one this script's bend and coupler simulations actually accumulate,
#     rather than some other meep run's.
# ---------------------------------------------------------------------------
def build_calibration_geometry():
    """The bare strip as a DeviceGeometry, i.e. solve_neff()'s cross-section
    in the one form FDTDSolver can also consume.

    Single source for the calibration call and its geometry plot, same rule as
    build_bend_polys() -- the picture cannot drift from what was simulated.
    """
    return UniformCrossSection(
        [(0.0, W)], length=CALIBRATE_SPAN,
        waveguide_height=MATERIAL_STACK.wg_height,
    )


def build_calibration_solver():
    """FDTDSolver configured from THIS script's constants.

    Every value here is the module-level one the ring's own bend and coupler
    runs use. That is the whole point: the beta being measured has to be the
    beta those simulations accumulate, which it only is if the mesh, the PML,
    the monitor placement, the decay criterion, the source bandwidth, the
    2-D/3-D choice and the material resolution all match.
    """
    return FDTDSolver(
        res=FDTD_RESOLUTION,
        # MESH_WIDTH, so the calibrated n_eff and solve_neff()'s meow n_eff
        # describe the same transverse window and the excess between them is a
        # statement about the two ENGINES rather than about two domain sizes.
        sim_width=MESH_WIDTH,
        sim_height=SIM_HEIGHT,
        three_d=THREE_D,
        dpml=DPML,
        source_margin=SOURCE_MARGIN,
        monitor_margin=MONITOR_MARGIN,
        num_input_modes=1,
        num_output_modes=1,
        output_monitor_centers=[0.0],
        material_stack=MATERIAL_STACK,
        # This script's own A/B override (M4); None correctly inherits the
        # stack's flag rather than hardcoding either answer.
        dispersive=FDTD_DISPERSIVE,
        decay_tol=DECAY_TOL,
        source_bandwidth_factor=SOURCE_BANDWIDTH_FACTOR,
    )


def calibrate_strip_beta(wavelengths, save_dir=None, cache_only=False):
    """`BetaCalibration` for the bare single-mode strip, from meep's own phase.

    Returns the full details object (not just `.beta`) so the caller can plot
    the fit behind the answer: the residual guard reduces a whole calibration
    to one number per cell and only raises past 0.1 rad, so a marginal fit
    passes silently and is invisible in the n_eff the straights end up on.

    `save_dir`/`cache_only` have the same meaning as everywhere else in this
    script; they are handed straight to the solver, which fingerprints the
    call itself (beef/cache/) rather than reusing this file's _sim_params().
    """
    return build_calibration_solver().section_mode_beta_details(
        build_calibration_geometry(), np.asarray(wavelengths, dtype=float),
        -CALIBRATE_SPAN / 2.0, CALIBRATE_SPAN / 2.0, 1,
        save_dir=save_dir, cache_only=cache_only,
    )


# ---------------------------------------------------------------------------
# 9. Interpolation helpers
# ---------------------------------------------------------------------------
# interp_complex now comes from beef.spectral (imported above) rather than
# being defined here. It was duplicated verbatim in this script and
# ring_resonator_hybrid_cascaded.py, had no test, and is load-bearing: see
# tests/test_spectral.py, which measures why real/imag interpolation is the
# wrong tool (it damps the coefficients and made Q look non-convergent) and
# pins this strategy's own failure mode.


def find_dips(wl_dense, T):
    """Numpy-only local-minima + FWHM/Q/finesse extraction (no scipy, matching
    every other example in this repo)."""
    dips = []
    baseline = np.median(T)
    for i in range(1, len(T) - 1):
        if T[i] < T[i - 1] and T[i] < T[i + 1] and T[i] < 0.7 * baseline:
            half = (T[i] + baseline) / 2.0
            # walk outward to the half-depth crossings for FWHM
            lo = i
            while lo > 0 and T[lo] < half:
                lo -= 1
            hi = i
            while hi < len(T) - 1 and T[hi] < half:
                hi += 1
            if hi > lo:
                fwhm = wl_dense[hi] - wl_dense[lo]
                q = wl_dense[i] / fwhm if fwhm > 0 else np.nan
                dips.append((wl_dense[i], T[i], fwhm, q))
    return dips


# ---------------------------------------------------------------------------
# 10. Geometry preview plots
# ---------------------------------------------------------------------------
def plot_cell(polys, cell_size, geometry_center, monitors, path, title):
    fig, ax = plt.subplots(figsize=(7, 5))
    cx, cy = geometry_center
    sx, sy = cell_size

    outer = mpatches.Rectangle((cx - sx / 2, cy - sy / 2), sx, sy)
    inner = mpatches.Rectangle((cx - sx / 2 + DPML, cy - sy / 2 + DPML), sx - 2 * DPML, sy - 2 * DPML)
    pml_path = mpath.Path.make_compound_path(
        mpath.Path(outer.get_verts()[::-1]), mpath.Path(inner.get_verts())
    )
    ax.add_patch(mpatches.PathPatch(pml_path, facecolor="crimson", alpha=0.15, edgecolor="none",
                                     zorder=0, label="PML"))

    for poly in polys:
        ax.add_patch(mpatches.Polygon(poly, closed=True, facecolor="steelblue", edgecolor="k", lw=0.3))
    ax.add_patch(mpatches.Rectangle((cx - sx / 2, cy - sy / 2), sx, sy, fill=False, edgecolor="gray", lw=0.8))
    ax.add_patch(mpatches.Rectangle((cx - sx / 2 + DPML, cy - sy / 2 + DPML), sx - 2 * DPML, sy - 2 * DPML,
                                     fill=False, edgecolor="gray", linestyle=":", lw=0.8))
    for (kind, x0, y0, x1, y1, label, color) in monitors:
        ax.plot([x0, x1], [y0, y1], color=color, lw=2.0, solid_capstyle="butt")
        ax.annotate(label, (0.5 * (x0 + x1), 0.5 * (y0 + y1)), fontsize=7, color=color,
                    ha="center", va="bottom")
    ax.set_xlim(cx - sx / 2 - 0.5, cx + sx / 2 + 0.5)
    ax.set_ylim(cy - sy / 2 - 0.5, cy + sy / 2 + 0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(title)
    ax.legend(handles=[mpatches.Patch(facecolor="crimson", alpha=0.15, label=f"PML (d={DPML} um)")],
              loc="upper right", fontsize=7)
    plt.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close("all")


def _mon_seg(width, x=None, y=None, xc=None, yc=None):
    if x is not None:
        return (x, yc - width / 2, x, yc + width / 2)
    return (xc - width / 2, y, xc + width / 2, y)


plot_cell(
    build_coupler_polys(), (COUPLER_SX, COUPLER_SY), (0.0, COUPLER_Y_CENTER),
    [
        ("src", *_mon_seg(PORT_APERTURE, x=-INNER_X_HALF, yc=0.0), "bus src", "green"),
        ("mon", *_mon_seg(PORT_APERTURE, x=-INNER_X_HALF + MONITOR_MARGIN, yc=0.0), "bus_L/norm", "orange"),
        ("mon", *_mon_seg(PORT_APERTURE, x=INNER_X_HALF - MONITOR_MARGIN, yc=0.0), "bus_R", "orange"),
        ("mon", *_mon_seg(PORT_APERTURE, xc=-STUB_X, y=Y_RING + R + S_MON), "ring_L", "purple"),
        ("mon", *_mon_seg(PORT_APERTURE, xc=STUB_X, y=Y_RING + R + S_MON), "ring_R", "purple"),
    ],
    os.path.join(sections_dir, "fdtd_coupler_section.png"),
    "Ring Resonator — Coupler Section (bus + near straight + 2 near bends)",
)
plot_cell(
    build_bend_polys(), (BEND_SX, BEND_SY), (BEND_X_CENTER, BEND_Y_CENTER),
    [
        ("src", *_mon_seg(PORT_APERTURE, x=-STUB_LEN + SOURCE_MARGIN, yc=0.0), "src", "green"),
        ("mon", *_mon_seg(PORT_APERTURE, x=-STUB_LEN + MONITOR_MARGIN, yc=0.0), "norm", "orange"),
        ("mon", *_mon_seg(PORT_APERTURE, xc=R, y=R + STUB_LEN - MONITOR_MARGIN), "out", "purple"),
    ],
    os.path.join(sections_dir, "fdtd_standalone_bend.png"),
    "Ring Resonator — Standalone 90 deg Bend (loss/phase reference)",
)
print(f"Section previews saved -> {sections_dir}/ (fdtd_coupler_section, fdtd_standalone_bend)")


# ---------------------------------------------------------------------------
# 10b. Reference-plane table + round-trip ledger.
#
#      Counterpart to ring_resonator_hybrid_cascaded.py's
#      device_assembly_ledger(): for every contribution to the ring's optical
#      path, the length this script's algebra assigns vs the physical span
#      between the reference planes that actually bound it. A non-zero mismatch
#      is unaccounted round-trip phase, and resonance position is round-trip
#      phase, so it is a wrong answer for the one quantity this device reports.
#
#      Added 2026-08-03 after the cascaded ledger found a 0.5 um shortfall
#      there. This device has the opposite problem and a second one besides:
#      its ring port is ASYMMETRIC (in and out planes differ), which no single
#      lever-arm constant can express. Built as a table rather than patched
#      term-by-term deliberately -- correcting L_outer alone was measured to
#      improve FSR agreement while making comb-position agreement worse, which
#      is what happens when several offsets interact.
# ---------------------------------------------------------------------------
ARC_LEN = np.pi * R / 2.0                    # quarter-circle centerline

# Ring-stub reference planes, as heights ABOVE the near bend's arc end.
RING_IN_PLANE = STUB_LEN - MONITOR_MARGIN    # ring-excited norm monitor, LEFT stub  (Y_TOP - MONITOR_MARGIN)
RING_OUT_PLANE = S_MON                       # ring_r_mon,                RIGHT stub (Y_RING + R + S_MON)
BUS_L_PLANE = -INNER_X_HALF + MONITOR_MARGIN
BUS_R_PLANE = INNER_X_HALF - MONITOR_MARGIN


def coefficient_reference_planes():
    """Which plane pair each of the four coefficients is referenced to.

    The point of tabulating it: `k1` and `k2` are NOT referenced to the same
    ring plane, because this script launches the ring-excited run from the stub
    TOP and reads the ring output at S_MON. Reciprocity (`|k1| == |k2|`) is
    only expected for coefficients at MATCHING planes, so the reciprocity print
    below compares magnitudes across a plane mismatch -- and the phases differ
    by 0.712*pi, measured 2026-08-03.
    """
    return [
        ("t",      f"bus x={BUS_L_PLANE:+.2f}", f"bus x={BUS_R_PLANE:+.2f}", "symmetric"),
        ("k1",     f"bus x={BUS_L_PLANE:+.2f}", f"ring +{RING_OUT_PLANE:.2f} above arc", "—"),
        ("k2",     f"ring +{RING_IN_PLANE:.2f} above arc", f"bus x={BUS_R_PLANE:+.2f}", "—"),
        ("t_ring", f"ring +{RING_IN_PLANE:.2f} above arc", f"ring +{RING_OUT_PLANE:.2f} above arc",
         f"ASYMMETRIC by {RING_IN_PLANE - RING_OUT_PLANE:.2f} um"),
    ]


def plot_reference_planes(path):
    """The four coefficients' reference planes, drawn.

    ADDED 2026-08-03 because the plane mismatch this script has always had is
    the thing the k1/k2 phase residual is measured against, and it existed only
    as a table (`coefficient_reference_planes()`). A table makes the numbers
    checkable; it does not make the ASYMMETRY visible, and the asymmetry is the
    point.

    What the figure has to get right, and what a naive version would get wrong:
    the two ring planes are on OPPOSITE stubs. `RING_IN_PLANE` (+1.25, the
    ring-excited run's normalization) is on the LEFT stub; `RING_OUT_PLANE`
    (+0.50, `ring_r_mon`) is on the RIGHT. Drawing them both on one stub would
    suggest 0.75 um of guide physically separates them, which is false -- they
    are separated by the whole coupler section. The 0.75 um is a difference in
    how far each sits ABOVE its own arc end, and the -beta*d correction is valid
    because the two stubs are identical isolated straight guide, mirrored.

    That "identical isolated straight guide" claim is now measured rather than
    assumed: the stub cross-section solves to n_eff = 2.402279 at 1.50 um and
    2.339523 at 1.55 um, bit-identical to the bare strip (nearest neighbour is
    ~5 um away). So the transport n_eff is NOT wrong for cross-section reasons.
    """
    fig, ax = plt.subplots(figsize=(9, 7))
    for poly in build_coupler_polys():
        ax.add_patch(MplPolygon(poly, closed=True, facecolor="none",
                                edgecolor="0.35", lw=1.0, zorder=2))

    arc_end = Y_RING + R
    y_in, y_out = arc_end + RING_IN_PLANE, arc_end + RING_OUT_PLANE

    def seg(x0, y0, x1, y1, color, ls, label):
        ax.plot([x0, x1], [y0, y1], color=color, lw=1.8, ls=ls, zorder=6,
                label=label)

    seg(*_mon_seg(PORT_APERTURE, x=BUS_L_PLANE, yc=0.0), "orange", "-",
        f"bus planes (x = {BUS_L_PLANE:+.2f}, {BUS_R_PLANE:+.2f}) — symmetric")
    seg(*_mon_seg(PORT_APERTURE, x=BUS_R_PLANE, yc=0.0), "orange", "-", None)
    seg(*_mon_seg(PORT_APERTURE, xc=-STUB_X, y=y_in), "crimson", "-",
        f"ring IN plane, LEFT stub  (+{RING_IN_PLANE:.2f} above arc) — k2, t_ring")
    seg(*_mon_seg(PORT_APERTURE, xc=STUB_X, y=y_out), "purple", "--",
        f"ring OUT plane, RIGHT stub (+{RING_OUT_PLANE:.2f} above arc) — k1, t_ring")

    # Mark each arc end. Without this the "+1.25 above arc" / "+0.50 above arc"
    # labels have no visible datum, which was the first version's flaw.
    lbl_a = "arc end (the datum both offsets are measured from)"
    for x in (-STUB_X, STUB_X):
        ax.plot([x - PORT_APERTURE / 2, x + PORT_APERTURE / 2], [arc_end, arc_end],
                color="0.55", lw=1.1, ls=":", zorder=5, label=lbl_a)
        lbl_a = None
    # Offsets annotated on the INSIDE of each stub, so neither runs into the frame.
    ax.annotate("", xy=(-STUB_X + 1.0, y_in), xytext=(-STUB_X + 1.0, arc_end),
                arrowprops=dict(arrowstyle="<->", color="crimson", lw=1.3))
    ax.text(-STUB_X + 1.15, (arc_end + y_in) / 2, f"+{RING_IN_PLANE:.2f}",
            color="crimson", fontsize=9.5, va="center", fontweight="bold")
    ax.annotate("", xy=(STUB_X - 1.0, y_out), xytext=(STUB_X - 1.0, arc_end),
                arrowprops=dict(arrowstyle="<->", color="purple", lw=1.3))
    ax.text(STUB_X - 1.15, (arc_end + y_out) / 2, f"+{RING_OUT_PLANE:.2f}",
            color="purple", fontsize=9.5, va="center", ha="right", fontweight="bold")

    d = RING_IN_PLANE - RING_OUT_PLANE
    ax.text(0.0, arc_end + STUB_LEN + 0.6,
            f"ring planes differ by {d:.2f} um — k1 and k2 are NOT referenced alike.\n"
            f"Correcting by exp(-i*beta*{d:.2f}) leaves a residual that DRIFTS "
            f"0.0068pi across the band\n(101-pt fit: a free propagation term beats "
            f"a pure constant 23x), so it is not a sign\nand not a fixed length. "
            f"Stub cross-section = bare strip exactly, so it is not the n_eff either.",
            ha="center", va="bottom", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.45", fc="#fff4e6", ec="crimson", lw=1.0))

    ax.set_xlim(-INNER_X_HALF - 1.0, INNER_X_HALF + 1.0)
    ax.set_ylim(-2.0, arc_end + STUB_LEN + 3.6)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("Monolithic ring — coefficient reference planes "
                 "(the k1/k2 mismatch, drawn)")
    ax.legend(loc="lower center", fontsize=8, framealpha=0.95)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def device_assembly_ledger():
    """Round-trip ledger. Rows with assigned=None are measured by FDTD between
    two planes and so carry no length assumption."""
    L_outer_shipped = 2 * (L_SIDE - S_MON) + L_FAR
    rows = [
        dict(span="coupler section: ring_in -> ring_out", modelled_by="t/k1/k2/t_ring (one FDTD)",
             assigned=None,
             physical=RING_IN_PLANE + ARC_LEN + L_C + ARC_LEN + RING_OUT_PLANE),
        dict(span="right side straight", modelled_by="L_outer (analytic neff)",
             assigned=L_SIDE - S_MON, physical=L_SIDE - RING_OUT_PLANE),
        dict(span="far-right bend arc", modelled_by="b_corrected (planes stripped)",
             assigned=ARC_LEN, physical=ARC_LEN),
        dict(span="far straight", modelled_by="L_outer (analytic neff)",
             assigned=L_FAR, physical=L_FAR),
        dict(span="far-left bend arc", modelled_by="b_corrected (planes stripped)",
             assigned=ARC_LEN, physical=ARC_LEN),
        # The defect: L_outer subtracts S_MON on BOTH sides, but the ring INPUT
        # plane is RING_IN_PLANE (1.25), not S_MON (0.5).
        dict(span="left side straight", modelled_by="L_outer (analytic neff)",
             assigned=L_SIDE - S_MON, physical=L_SIDE - RING_IN_PLANE),
    ]
    for r in rows:
        r["mismatch"] = None if r["assigned"] is None else r["assigned"] - r["physical"]
    return rows, L_outer_shipped


def print_assembly_ledger():
    rows, L_shipped = device_assembly_ledger()
    print("\nCoefficient reference planes (why |k1|==|k2| is not a like-for-like test)")
    print(f"  {'coeff':8s}{'normalized at':34s}{'measured at':34s}{'note'}")
    print("  " + "-" * 100)
    for name, src, dst, note in coefficient_reference_planes():
        print(f"  {name:8s}{src:34s}{dst:34s}{note}")

    print("\nRound-trip ledger — assigned length vs physical span between reference planes")
    print(f"  {'span':40s}{'modelled by':32s}{'assigned':>10s}{'physical':>10s}{'mismatch':>10s}")
    print("  " + "-" * 102)
    total_mismatch = 0.0
    for r in rows:
        a = "measured" if r["assigned"] is None else f"{r['assigned']:.4f}"
        m = "     —" if r["mismatch"] is None else f"{r['mismatch']:+.4f}"
        flag = ""
        if r["mismatch"] is not None:
            total_mismatch += r["mismatch"]
            if abs(r["mismatch"]) > 1e-9:
                flag = "  <-- MISMATCH"
        print(f"  {r['span']:40s}{r['modelled_by']:32s}{a:>10s}{r['physical']:10.4f}{m:>10s}{flag}")
    print(f"\n  Bend straight-stub budget: b^2 carries {BEND_STRAIGHT_TOTAL:.4f} um of straight "
          f"(phase stripped, magnitude NOT)")
    print(f"    vs L_OUTER = {L_OUTER:.4f} um of physical straight  ->  DC-style equivalent "
          f"analytic length {BEND_STRAIGHT_HEADROOM:.4f} um")
    print(f"    HEADROOM {BEND_STRAIGHT_HEADROOM:.4f} um before the strip becomes MANDATORY "
          f"(guarded at import; the DC's own assert is `body_lengths[0] >= 0`)")
    physical_total = sum(r["physical"] for r in rows)
    print("  " + "-" * 102)
    print(f"  {'physical round trip (sum of spans)':40s}{'':32s}{'':>10s}{physical_total:10.4f}")
    print(f"  {'PERIMETER constant':40s}{'':32s}{'':>10s}{PERIMETER:10.4f}"
          f"   (delta {physical_total - PERIMETER:+.2e})")
    L_true = (L_SIDE - RING_OUT_PLANE) + (L_SIDE - RING_IN_PLANE) + L_FAR
    print(f"\n  L_outer as shipped = {L_shipped:.4f} um   "
          f"plane-consistent value = {L_true:.4f} um   (over-count {L_shipped - L_true:+.4f})")
    if abs(total_mismatch) > 1e-9:
        print(f"  ** TOTAL UNACCOUNTED PATH: {total_mismatch:+.4f} um ** — the modelled round trip is "
              f"LONGER than the device.")
        print("     NOT patched here: correcting L_outer alone was measured to improve route FSR "
              "agreement (2.9% -> 0.24%)\n     while making comb-position agreement WORSE "
              "(0.31 -> 0.45 FSR), i.e. at least one further offset\n     interacts with it — most "
              "likely the k1/k2 ring-plane mismatch in the table above.")
    else:
        print("  Every assigned length matches its physical span. Bookkeeping closes.")
    return total_mismatch


print_assembly_ledger()
plot_reference_planes(os.path.join(sections_dir, "reference_planes.png"))
print(f"Reference-plane diagram saved -> {sections_dir}/reference_planes.png")


# ---------------------------------------------------------------------------
# 11. Run the 3 FDTD simulations
# ---------------------------------------------------------------------------
wl_fdtd = np.linspace(WL_MIN, WL_MAX, N_WL_FDTD)

print("\n[1/3] Coupler-section FDTD, bus-excited...")
res_a = run_coupler_fdtd("bus", wl_fdtd, save_dir=save_dir)
t, k1 = res_a["t"], res_a["k1"]

print("[2/3] Coupler-section FDTD, ring-excited...")
res_b = run_coupler_fdtd("ring", wl_fdtd, save_dir=save_dir)
t_ring, k2 = res_b["t_ring"], res_b["k2"]

print("[3/3] Standalone bend FDTD...")
b_raw = run_bend_fdtd(wl_fdtd, save_dir=save_dir)

mid = N_WL_FDTD // 2
conserve_a = np.abs(t[mid]) ** 2 + np.abs(k1[mid]) ** 2
conserve_b = np.abs(t_ring[mid]) ** 2 + np.abs(k2[mid]) ** 2
print(f"\nPower-conservation check @ mid-lambda (should be <~1, both runs are lossy/leaky by design "
      f"since each only captures 2 of several outgoing ports):")
print(f"  run A (bus-excited): |t|^2+|k1|^2 = {conserve_a:.4f}")
print(f"  run B (ring-excited): |t_ring|^2+|k2|^2 = {conserve_b:.4f}")
print(f"  reciprocity check: |k1|={np.abs(k1[mid]):.4f} vs |k2|={np.abs(k2[mid]):.4f} (should be close)")

# ---------------------------------------------------------------------------
# 12. Mode solve — sparse neff(lambda), degree-2 polynomial fit
# ---------------------------------------------------------------------------
print(f"\nMode solve at {N_WL_SOLVE} wavelengths...")
wl_solve = np.linspace(WL_MIN, WL_MAX, N_WL_SOLVE)
neff_mpb = np.array([solve_neff(wl, save_dir=save_dir) for wl in wl_solve])

# The MPB curve is solved unconditionally even when the calibration is on --
# it costs nothing new (13 cached mode solves) and it is the only thing the
# calibrated curve can be read against. `seed_excess` on the calibration says
# the same thing per mode, but against the FDTD solver's OWN seed at
# MODE_RESOLUTION-independent settings, not against the res-100 curve this
# script would otherwise have propagated on.
neff_solve = neff_mpb
calibration = None
_calibrate_now = CALIBRATE_BODY and THREE_D
if CALIBRATE_BODY and not THREE_D:
    # A 2-D run has no (y, z) cross-section, so there is nothing here that
    # corresponds to the 3-D strip the calibration would measure. Skip loudly:
    # silently calibrating against a different structure is the failure mode
    # this device's own validation history is full of.
    print("  RING_CALIBRATE_BODY=1 but this is a 2-D run (RING_THREE_D=0) — "
          "calibration SKIPPED, straights stay on the meow/MPB n_eff.")
if _calibrate_now:
    print(f"  Calibrating the straights on meep's own beta "
          f"({CALIBRATE_SPAN:.1f} um span, res {FDTD_RESOLUTION:g})...")
    calibration = calibrate_strip_beta(wl_solve, save_dir=save_dir)
    neff_solve = calibration.beta[:, 0] * wl_solve / (2.0 * np.pi)
    _excess = neff_solve / neff_mpb - 1.0
    _i = len(wl_solve) // 2
    print(f"  Bare-strip n_eff at lambda = {wl_solve[_i]:.3f} um:")
    print(f"    MPB mode index  {neff_mpb[_i]:.6f}")
    print(f"    FDTD own phase  {neff_solve[_i]:.6f}   "
          f"<- the straights propagate on this")
    print(f"    excess {100*_excess[_i]:+.3f}%  "
          f"(band range {100*_excess.min():+.3f}% .. {100*_excess.max():+.3f}%)")
    # What that excess is worth on the quantity this device reports. L_OUTER is
    # the analytic lever arm; the bend correction rides the same n_eff but with
    # the opposite sign, so this is the straight-phase term alone.
    _dphi = 2 * np.pi * (neff_solve[_i] - neff_mpb[_i]) * L_OUTER / wl_solve[_i]
    print(f"    -> round-trip phase shift {_dphi:+.4f} rad "
          f"= {100*_dphi/(2*np.pi):+.1f}% of a full FSR")

    # The measurement behind the answer. Free -- no extra simulation, just the
    # data section_mode_beta_details() already produced. Same builders the
    # calibration itself used, so the drawn domain is the simulated one.
    plot_beta_calibration(
        calibration, build_calibration_geometry(), build_calibration_solver(),
        os.path.join(output_dir, "calibration"),
    )
    print(f"  Calibration plots saved -> {output_dir}/calibration/ "
          f"(geometry/ + fits/)")

neff_fit = np.polyfit(wl_solve, neff_solve, deg=2)
neff_resid = neff_solve - np.polyval(neff_fit, wl_solve)
print(f"  neff fit residual (max abs): {np.max(np.abs(neff_resid)):.2e}")

fig, ax = plt.subplots(figsize=(6, 4))
wl_fine = np.linspace(WL_MIN, WL_MAX, 200)
if _calibrate_now:
    ax.plot(wl_solve, neff_mpb, "s", color="0.55", mfc="none",
            label="meow/MPB mode index (not propagated on)")
ax.plot(wl_solve, neff_solve, "o",
        label="FDTD-calibrated (propagated on)" if _calibrate_now else "mode-solve points")
ax.plot(wl_fine, np.polyval(neff_fit, wl_fine), "-", label="degree-2 fit")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("neff")
ax.set_title("Bare strip effective index"
             + (" — FDTD-calibrated" if _calibrate_now else ""))
ax.legend()
ax.grid(True)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "neff_lambda.png"), dpi=200)
plt.close("all")

# ---------------------------------------------------------------------------
# 13. Dense-grid assembly
# ---------------------------------------------------------------------------
wl_dense = np.linspace(WL_MIN, WL_MAX, N_WL_DENSE)
neff_dense = np.polyval(neff_fit, wl_dense)

t_d = interp_complex(wl_fdtd, t, wl_dense)
k1_d = interp_complex(wl_fdtd, k1, wl_dense)
k2_d = interp_complex(wl_fdtd, k2, wl_dense)
t_ring_d = interp_complex(wl_fdtd, t_ring, wl_dense)
b_d = interp_complex(wl_fdtd, b_raw, wl_dense)

# Reference-plane correction: the bend FDTD's raw coefficient is measured
# between its norm/output monitor planes, not the true arc endpoints —
# strip out the straight-waveguide propagation phase over BEND_D_IN +
# BEND_D_OUT so b_corrected represents the arc alone.
phi_bend_correction = 2 * np.pi * neff_dense * (BEND_D_IN + BEND_D_OUT) / wl_dense
b_corrected = b_d * np.exp(-1j * phi_bend_correction)

L_outer = L_OUTER                              # module-level, guarded above
phi = 2 * np.pi * neff_dense * L_outer / wl_dense
P = (b_corrected ** 2) * np.exp(1j * phi)

E_thru = t_d + k1_d * k2_d * P / (1.0 - t_ring_d * P)
T = np.abs(E_thru) ** 2

dips = find_dips(wl_dense, T)
if len(dips) >= 2:
    fsr_meas = np.mean(np.diff([d[0] for d in dips]))
else:
    fsr_meas = np.nan
print(f"\nDense-grid assembly ({N_WL_DENSE} points): {len(dips)} resonance dip(s) found.")
print(f"  Measured mean FSR: {fsr_meas*1000:.3f} nm (predicted: {FSR_GUESS*1000:.2f} nm)")
if dips:
    print(f"  {'lambda (um)':>12}  {'T_min':>8}  {'FWHM (nm)':>10}  {'Q':>10}")
    for wl0, Tmin, fwhm, q in dips[:10]:
        print(f"  {wl0:12.4f}  {Tmin:8.4f}  {fwhm*1000:10.4f}  {q:10.1f}")

wl_c = wl_dense[len(wl_dense) // 2]
idx = np.argmin(np.abs(wl_dense - wl_c))
print(f"\nResults at lambda = {wl_c:.4f} um:")
print(f"  |t|={np.abs(t_d[idx]):.4f}  |k1|={np.abs(k1_d[idx]):.4f}  |k2|={np.abs(k2_d[idx]):.4f}  "
      f"|t_ring|={np.abs(t_ring_d[idx]):.4f}  |b|={np.abs(b_corrected[idx]):.4f}")
print(f"  T_through = {T[idx]:.4f}  IL = {10*np.log10(max(T[idx],1e-20)):.2f} dB")

# ---------------------------------------------------------------------------
# 14. Independent bend / coupler-region transmission plots
# ---------------------------------------------------------------------------
T_bend = np.abs(b_d) ** 2  # note: the reference-plane phase correction (-> b_corrected)
                            # only rotates the phase, so it doesn't change |b|^2 — this
                            # magnitude curve is the same whether read off b_d or b_corrected

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(wl_dense, T_bend, color="steelblue")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Transmission")
ax.set_title("Standalone 90 deg Bend — Transmission")
ax.grid(True)
ax.set_ylim(0, 1.05)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "bend_transmission.png"), dpi=300)
plt.close("all")

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(wl_dense, np.abs(t_d) ** 2, label="bus->bus |t|^2")
ax.plot(wl_dense, np.abs(k1_d) ** 2, linestyle="--", label="bus->ring |k1|^2")
ax.plot(wl_dense, np.abs(k2_d) ** 2, linestyle=":", label="ring->bus |k2|^2")
ax.plot(wl_dense, np.abs(t_ring_d) ** 2, linestyle="-.", label="ring->ring |t_ring|^2")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Transmission / coupling")
ax.set_title("Coupler Section — Transmission & Coupling Coefficients")
ax.legend()
ax.grid(True)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "coupler_transmission.png"), dpi=300)
plt.close("all")
print(f"Bend/coupler transmission plots saved -> {output_dir}/bend_transmission.png, coupler_transmission.png")

# ---------------------------------------------------------------------------
# 15. Save and plot
# ---------------------------------------------------------------------------
smat = SMatrix(wavelengths=wl_dense)
smat.set(E_thru, port="through")
smat.check_passive()
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix saved -> {output_dir}/smatrix_hybrid.pkl")

np.savez(
    os.path.join(output_dir, "ring_data.npz"),
    wl_fdtd=wl_fdtd, wl_dense=wl_dense, wl_solve=wl_solve,
    t=t, k1=k1, k2=k2, t_ring=t_ring, b_raw=b_raw,
    neff_solve=neff_solve, neff_dense=neff_dense,
    # The uncalibrated curve alongside the propagated one. Without it the npz
    # cannot say WHICH engine's n_eff produced the comb it holds, and the
    # calibrated/uncalibrated comparison would need a second full run to
    # reconstruct. Equal to neff_solve when the calibration is off or skipped,
    # which is itself the record that it was.
    calibrated=_calibrate_now, neff_mpb=neff_mpb,
    T=T, dips=np.array([(d[0], d[1], d[2], d[3]) for d in dips]) if dips else np.empty((0, 4)),
    perimeter=PERIMETER, fsr_guess=FSR_GUESS,
)
print(f"Data saved -> {output_dir}/ring_data.npz")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wl_dense, T, color="steelblue", lw=1.0)
for wl0, Tmin, fwhm, q in dips:
    ax.plot(wl0, Tmin, "rv", ms=4)
if len(dips) >= 2:
    ax.annotate("", xy=(dips[1][0], 1.0), xytext=(dips[0][0], 1.0),
                arrowprops=dict(arrowstyle="<->", color="gray"))
    ax.text(0.5 * (dips[0][0] + dips[1][0]), 1.02, f"FSR~{fsr_meas*1000:.1f}nm",
            ha="center", fontsize=8, color="gray")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Transmission")
ax.set_title("Ring Resonator — Through-port Transmission")
ax.grid(True)
ax.set_ylim(0, 1.1)

ax = axes[1]
IL = 10 * np.log10(np.clip(T, 1e-20, None))
ax.plot(wl_dense, IL, color="steelblue", lw=1.0)
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title("Ring Resonator — Insertion Loss")
ax.grid(True)

plt.suptitle("Racetrack Ring Resonator — Hybrid Simulation (raw-meep coupler+bend FDTD, analytic straights)",
             fontsize=10)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved -> {output_dir}/hybrid_results.png")
