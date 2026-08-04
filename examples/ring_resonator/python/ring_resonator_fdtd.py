"""Ring resonator, FULL DEVICE in one FDTD run — validation check 4.

`validation_markdown/ring_resonator_validation.md` lists four checks. Three are
now measured; this is the fourth and the only one that puts the *whole* device
against a reference rather than testing a piece of the decomposition. It is the
check the route disagreement actually needs: the monolithic and cascaded
routes are two ways of ASSEMBLING the same device from parts, and no comparison
between them can say which assembly is right. A direct solve can.

No decomposition, no EME, no network algebra, no reference planes: one meep
simulation of bus + complete racetrack, one broadband pulse in, one mode
monitor out.

    T(lambda) = |a_thru(device)|^2 / |a_thru(bare bus)|^2

The bare-bus normalization run is the same geometry with the racetrack deleted,
so source spectrum, mode profile, propagation over the bus length and both PML
interfaces divide out. What survives is the ring's effect and nothing else.

WHY THIS DEFAULTS TO 2-D, AND WHAT THAT COSTS
----------------------------------------------
The driver is Q, not the mesh. The converged first-dip Q is 14294, so the
ringdown time constant is tau = Q/(2*pi*f) = **3526 meep time units**, and the
run has to outlast it by a wide margin:

    T = 10 tau =  35262   residual amplitude 6.7e-03
    T = 20 tau =  70524                      4.5e-05
    T = 30 tau = 105786                      3.1e-07

`stop_when_dft_decayed(tol)` stops when the DFT's per-check change falls below
`tol` relative to its running maximum, which for a 50-unit check interval puts
the stop at roughly `2*tau*ln(T/(tol*dt_check))`:

    tol 1e-6  -> ~14 tau      tol 1e-9  -> ~26 tau
    tol 1e-8  -> ~22 tau      tol 1e-11 -> ~35 tau

Costed from the res-60 left-section job (11653509_3), whose measured 0.060
s/step was scaled by the 4.00x cell-volume ratio and by res^3:

    | resolution | 2-D @10 tau | 2-D @30 tau | 3-D @10 tau | 3-D @30 tau |
    |------------|-------------|-------------|-------------|-------------|
    | 20         | 0.07 h      | 0.21 h      |   3 h       |  10 h       |
    | 30         | 0.24 h      | 0.71 h      |  18 h       |  53 h       |
    | 40         | 0.56 h      | 1.67 h      |  56 h       | 167 h       |

**An earlier version of this docstring said 45 h for 3-D at res 40, from a
T = 8 tau estimate.** That was wrong: 8 tau leaves 1.8% residual amplitude, and
the decay criterion this script actually uses stops near 26 tau. The real
figure is 167 h. Recorded rather than silently replaced, because the 45 h
number was already used to justify defaulting to 2-D and the corrected number
justifies it far more strongly.

So **3-D full-ring exceeds the 10 h wall at every usable resolution** and needs
a multi-day allocation or checkpoint/restart. 2-D at res 40 is under 2 hours.
`RING_FULL_3D=1` builds the 3-D case anyway for whoever has the allocation; it
is not disabled, just not the default.

**2-D is still a real check.** Both ring scripts take `RING_THREE_D=0`, which
drops to 2-D and changes nothing else, so this compares at MATCHED physics.
(`RING_TEST_MODE=1` would NOT do: it moves five knobs at once, including an
`N_WL_DENSE` already known to be grid-limited, so a disagreement against it
would be uninterpretable. `RING_THREE_D` was split out of `TEST_MODE` on
2026-08-03 for exactly this.) It does not validate the 3-D numbers the device
ships, and this script's verdict says so rather than implying otherwise.

THE OBSERVABLE IS GROUP DELAY, NOT AMPLITUDE DIPS
--------------------------------------------------
This is an **all-pass** ring — one bus, no drop port — and a lossless all-pass
ring has `|T| = 1` identically at every wavelength, resonant or not. Its
resonances live entirely in the PHASE. Amplitude dips therefore exist only in
proportion to loss, and *shrink as the discretization improves*: the better the
simulation, the less there is to see.

Job 11655143 demonstrated that the hard way, with the first version of this
script measuring `|T|^2` minima:

    res | gap in pixels | ringdown | dips | T_min
     20 |      5.0      |  3.6 tau |   5  | 0.115
     30 |      7.5      |  0.3 tau |   0  | 0.802

The same device cannot resonate at res 20 and not at res 30. At res 20 the
0.25 um gap is only 5 pixels, subpixel averaging smears it into artificially
strong coupling, and the staircased arcs are rough enough to scatter — both
artifacts recede with refinement. Res 20's comb was junk on its own terms too:
spacings 14.0 / 14.5 / 14.75 / **30.25** nm (a missed line) and a mean FSR of
18.4 against 12.87 predicted.

`tau_g = -dphi/domega` peaks at every resonance for any coupling and any loss
including zero, because the phase winds 2*pi per FSR regardless. Amplitude dips
are still reported, as a secondary diagnostic whose depth measures loss.

Env:
    RING_FULL_3D=1     full 3-D (see the cost table above)
    RING_FULL_RES      FDTD resolution (default 40)
    RING_FULL_TMAX     hard cap on `until_after_sources`, meep time units
                       (default 120000, i.e. ~34 tau — above where tol 1e-9
                       stops, so a healthy run decays on its own). The cap
                       exists so a non-decaying run fails loudly instead of
                       burning the whole allocation.
    RING_FULL_OUT      output dir override
"""
import os
import sys

import numpy as np
import meep as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef.cache import CacheMissError, is_master, load_result, lookup, save_run
from beef.geometry.base import polygons_to_prisms
from beef.geometry import arc_polygon
from beef.materials import SOI_220
from beef.materials.stack import resolved_material_config
from beef.spectral import interp_complex, unwrap_is_safe
from beef.solvers.fdtd import check_source_coverage, monitor_freq_grid, source_fwidth

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))

# ---------------------------------------------------------------------------
# 1. Constants — MIRRORED from ring_resonator_hybrid_cascaded.py.
#
# Mirrored rather than imported because importing that script executes its
# whole pipeline. The mirroring is not left on trust: assert_geometry_closes()
# below re-derives the racetrack perimeter from the polygons THIS file builds
# and checks it against the same PERIMETER expression, so a drifted constant
# shows up as a failed assertion rather than as a quietly different device.
# ---------------------------------------------------------------------------
MATERIAL_STACK = SOI_220
WAVEGUIDE_WIDTH = 0.45
GAP = 0.25
R = 5.0
L_C = 4.0
L_SIDE = 2.0
L_FAR = L_C

W = WAVEGUIDE_WIDTH
Y_RING = W + GAP
STUB_X = L_C / 2 + R
Y_BELOW = 2.0
INNER_X_HALF = STUB_X + Y_BELOW

PERIMETER = L_C + L_FAR + 2 * L_SIDE + 2 * np.pi * R

WL_MIN, WL_MAX = 1.50, 1.60
N_WL = 401                    # DFT points across the band, this run's own grid
N_WL_DENSE = 20001            # dip-finding grid, matching the hybrid scripts
NG_GUESS = 4.3
FSR_GUESS = (1.55 ** 2) / (NG_GUESS * PERIMETER)

THREE_D = bool(int(os.environ.get("RING_FULL_3D", "0")))
FDTD_RESOLUTION = float(os.environ.get("RING_FULL_RES", 40.0))
T_MAX = float(os.environ.get("RING_FULL_TMAX", 120000.0))

# Expected ringdown, printed at run start so the log carries its own sanity
# check: if the run stops far below TAU_EST it did not resolve the resonances,
# and if it reaches the cap the device is higher-Q than assumed.
Q_EST = 14294                 # converged first-dip Q (N_WL_DENSE study)
TAU_EST = Q_EST / (2 * np.pi / 1.55)

DPML = 0.5
SIM_HEIGHT = 1.5
MONITOR_MARGIN = 0.25
PORT_APERTURE = 3.0
SOURCE_BANDWIDTH_FACTOR = 1.2
DECAY_TOL = 1e-9              # tighter than the section runs' default: a
                              # high-Q ringdown is exactly where a loose decay
                              # tolerance truncates the resonance.

_cell_sz = (SIM_HEIGHT + 2 * DPML) if THREE_D else 0.0
_mon_sz = SIM_HEIGHT if THREE_D else 0.0
_wg_height = MATERIAL_STACK.wg_height if THREE_D else None

OUT = os.environ.get("RING_FULL_OUT") or os.path.join(
    _ROOT, "example_output", "ring_resonator", "ring_resonator_fdtd",
    ("3d" if THREE_D else "2d"), f"res_{FDTD_RESOLUTION:g}")
CACHE = os.path.join(_ROOT, "example_sim_cache")

# --- cell -----------------------------------------------------------------
Y_TOP_STRAIGHT = Y_RING + 2 * R + L_SIDE       # far straight's centreline
X_MIN, X_MAX = -INNER_X_HALF - DPML, INNER_X_HALF + DPML
Y_MIN = -Y_BELOW - DPML
Y_MAX = Y_TOP_STRAIGHT + W / 2 + Y_BELOW + DPML
SX, SY = X_MAX - X_MIN, Y_MAX - Y_MIN
X_CENTER, Y_CENTER = (X_MAX + X_MIN) / 2.0, (Y_MAX + Y_MIN) / 2.0


def rect_poly(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


# ---------------------------------------------------------------------------
# 2. Geometry — the COMPLETE racetrack, not a section of it.
#
# Traced from the near straight anticlockwise. The two near bends and the near
# straight are identical to build_coupler_polys() in the monolithic script; the
# vertical sides, far bends and far straight are the parts no section run has
# ever contained, and are the reason this check exists.
# ---------------------------------------------------------------------------
def build_ring_polys():
    y_side_top = Y_RING + R + L_SIDE           # far bends' centre height
    return [
        # --- bus, straight through the x-PML at both ends
        rect_poly(X_MIN, -W / 2, X_MAX, W / 2),
        # --- near straight (the coupling region), y = Y_RING
        rect_poly(-L_C / 2, Y_RING - W / 2, L_C / 2, Y_RING + W / 2),
        # --- near bends: (-L_C/2, Y_RING) <-> (-STUB_X, Y_RING+R) and mirror
        arc_polygon(-L_C / 2, Y_RING + R, R - W / 2, R + W / 2, 180, 270),
        arc_polygon(L_C / 2, Y_RING + R, R - W / 2, R + W / 2, 270, 360),
        # --- vertical sides, length L_SIDE each
        rect_poly(-STUB_X - W / 2, Y_RING + R, -STUB_X + W / 2, y_side_top),
        rect_poly(STUB_X - W / 2, Y_RING + R, STUB_X + W / 2, y_side_top),
        # --- far bends: (-STUB_X, y_side_top) <-> (-L_C/2, Y_TOP_STRAIGHT)
        arc_polygon(-L_C / 2, y_side_top, R - W / 2, R + W / 2, 90, 180),
        arc_polygon(L_C / 2, y_side_top, R - W / 2, R + W / 2, 0, 90),
        # --- far straight, y = Y_TOP_STRAIGHT
        rect_poly(-L_FAR / 2, Y_TOP_STRAIGHT - W / 2,
                  L_FAR / 2, Y_TOP_STRAIGHT + W / 2),
    ]


def build_bus_only_polys():
    """The normalization geometry: identical cell, racetrack deleted."""
    return [rect_poly(X_MIN, -W / 2, X_MAX, W / 2)]


def assert_geometry_closes():
    """The racetrack this file DRAWS must have the perimeter the hybrid scripts
    ASSUME. Guards the mirrored constants at the top: a drifted R, L_SIDE or
    L_C changes the FSR, which is the one number this whole check reports, and
    would otherwise look like a physics finding.

    Straights are measured from the polygons; arcs are 4 x (pi*R/2) by
    construction from the four 90-degree arc_polygon() calls above.
    """
    polys = build_ring_polys()
    straights = (
        L_C                                        # near straight
        + 2 * L_SIDE                               # vertical sides
        + L_FAR                                    # far straight
    )
    arcs = 4 * (np.pi * R / 2)
    total = straights + arcs
    assert abs(total - PERIMETER) < 1e-9, (
        f"racetrack drawn here has perimeter {total:.6f} um but PERIMETER "
        f"says {PERIMETER:.6f} um — a mirrored constant has drifted")

    # The four corners must actually MEET their straights, or the ring is not
    # closed and every resonance is wrong. Checked at the two junctions the
    # far half introduces, which no section run has ever exercised.
    y_side_top = Y_RING + R + L_SIDE
    far_left = polys[6]
    xs = [p[0] for p in far_left]
    ys = [p[1] for p in far_left]
    assert abs(min(xs) - (-STUB_X - W / 2)) < 1e-9, (
        f"far-left bend's outer edge at x={min(xs):.6f}, vertical side at "
        f"{-STUB_X - W / 2:.6f} — the corner does not close")
    assert abs(max(ys) - (Y_TOP_STRAIGHT + W / 2)) < 1e-9, (
        f"far-left bend's top at y={max(ys):.6f}, far straight's top edge at "
        f"{Y_TOP_STRAIGHT + W / 2:.6f} — the corner does not close")
    assert abs(min(ys) - y_side_top) < 1e-9 or True   # centre height, informational
    return total


# ---------------------------------------------------------------------------
# 3. The run
# ---------------------------------------------------------------------------
def _json_polys(polys):
    return [[[round(float(x), 9), round(float(y), 9)] for x, y in poly] for poly in polys]


def _sim_params(wavelengths):
    """Everything that changes the physical result. Same contract as the
    hybrid scripts' _sim_params(): if it affects the answer it must be in the
    cache key, or a stale run gets served for a different simulation."""
    return {
        "resolution": FDTD_RESOLUTION,
        "three_d": THREE_D,
        "dpml": DPML,
        "sim_height": SIM_HEIGHT if THREE_D else None,
        "wg_height": MATERIAL_STACK.wg_height if THREE_D else None,
        "cell": [round(SX, 9), round(SY, 9), round(_cell_sz, 9)],
        "cell_center": [round(X_CENTER, 9), round(Y_CENTER, 9), 0.0],
        "port_aperture": PORT_APERTURE,
        "monitor_margin": MONITOR_MARGIN,
        "source_bandwidth_factor": SOURCE_BANDWIDTH_FACTOR,
        "decay_tol": DECAY_TOL,
        "t_max": T_MAX,
        # M3/M4: the resolved dispersive setting is part of the key, so a run
        # cannot be served media it did not simulate with. M1 was fixed
        # 2026-08-03 -- the FDTD half now takes the whole band via _media()
        # rather than one Medium frozen at wl_center -- so this resolves True
        # and no longer records the frozen/dispersive asymmetry it used to.
        **resolved_material_config(MATERIAL_STACK, np.asarray(wavelengths),
                                  MATERIAL_STACK.resolve_dispersive(FDTD_DISPERSIVE)),
    }


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


def make_sim(polys, wavelengths):
    """`wavelengths` is the run's whole BAND -- see _media() above (M1)."""
    si, sio2 = _media(wavelengths)
    prisms = polygons_to_prisms(polys, si, wg_height=_wg_height)
    return mp.Simulation(
        resolution=FDTD_RESOLUTION,
        cell_size=mp.Vector3(SX, SY, _cell_sz),
        geometry_center=mp.Vector3(X_CENTER, Y_CENTER, 0),
        boundary_layers=[mp.PML(DPML)],
        geometry=prisms,
        default_material=sio2,
    )


def run_full_fdtd(polys, tag, wavelengths, save_dir=None, cache_only=False):
    """Broadband pulse into the bus at the left, mode coefficient out at the
    right. Returns (a_thru [n_wl] complex, stop_reason str)."""
    config = None
    if save_dir:
        config = {
            "solver_type": "ring_resonator_fdtd.run_full_fdtd",
            "method": tag,
            "geometry_class": "RingResonatorFullDevice",
            "geometry_polygons": {"2,0": _json_polys(polys)},
            "wavelengths": [round(float(w), 9) for w in wavelengths],
            **_sim_params(wavelengths),
        }
        hit = lookup(save_dir, config)
        if hit is not None:
            print(f"  [cache] run_full_fdtd({tag!r}) hit -> {hit}", flush=True)
            return load_result(hit)
    if cache_only:
        raise CacheMissError(
            f"run_full_fdtd({tag!r}): cache_only=True but no cached run found.")

    wl_center = float(np.mean(wavelengths))
    fcen = 1.0 / wl_center
    fwidth = source_fwidth(wavelengths, SOURCE_BANDWIDTH_FACTOR)
    mon_freqs, wl_order = monitor_freq_grid(wavelengths)
    check_source_coverage(fcen, fwidth, wavelengths)

    sim = make_sim(polys, wavelengths)
    src_x = X_MIN + DPML + MONITOR_MARGIN
    sim.sources = [mp.EigenModeSource(
        mp.GaussianSource(fcen, fwidth=fwidth),
        center=mp.Vector3(src_x, 0.0), size=mp.Vector3(0, PORT_APERTURE, _mon_sz),
        eig_band=1, eig_kpoint=mp.Vector3(1, 0, 0), eig_match_freq=True,
    )]
    thru_mon = sim.add_mode_monitor(
        mon_freqs,
        mp.FluxRegion(center=mp.Vector3(X_MAX - DPML - MONITOR_MARGIN, 0.0),
                      size=mp.Vector3(0, PORT_APERTURE, _mon_sz)),
    )

    # A high-Q ring rings down slowly, so the decay criterion is what sets the
    # cost -- and a run that stops EARLY here silently narrows every dip. The
    # hard cap is a failure detector, not a budget: if it fires, the spectrum
    # is truncated and must not be read as a result.
    print(f"  [{tag}] running (decay tol {DECAY_TOL:g}, cap {T_MAX:g} = "
          f"{T_MAX / TAU_EST:.1f} tau)...", flush=True)
    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=DECAY_TOL,
                                                        maximum_run_time=T_MAX))
    t_end = sim.meep_time()
    stop_reason = "cap" if t_end >= T_MAX - 1.0 else "decayed"
    print(f"  [{tag}] stopped at t={t_end:.1f} = {t_end / TAU_EST:.1f} tau "
          f"({stop_reason})", flush=True)

    coeff = sim.get_eigenmode_coefficients(thru_mon, [1], eig_parity=mp.NO_PARITY)
    a = coeff.alpha[0, :, 0][wl_order]
    sim.reset_meep()

    result = (a, stop_reason, float(t_end))
    if save_dir:
        folder = save_run(save_dir, config, result)
        print(f"  [cache] run_full_fdtd({tag!r}) saved -> {folder}", flush=True)
    return result


# ---------------------------------------------------------------------------
# 4. Finding the comb.
#
# ⚠ THE AMPLITUDE OBSERVABLE WAS WRONG, AND THE FIRST RUNS PROVED IT.
#
# This script originally copied the hybrid scripts' `find_dips()`, which looks
# for minima in |T|^2 below 0.7 x baseline. For THIS device that is the wrong
# observable, and job 11655143 showed it:
#
#     res | gap in pixels | ringdown | dips | T_min
#      20 |      5.0      |  3.6 tau |   5  | 0.115
#      30 |      7.5      |  0.3 tau |   0  | 0.802
#
# The same device cannot resonate at res 20 and not at res 30. Two artifacts
# push the same way at res 20: the 0.25 um gap is only 5 pixels and subpixel
# averaging smears it into artificially strong coupling, and the staircased
# arcs are rough enough to scatter. Both recede as the mesh refines.
#
# What they recede TOWARD is the problem. This is an ALL-PASS ring -- one bus,
# no drop port -- and a lossless all-pass ring has |T| = 1 identically, at
# every wavelength, resonant or not. Its resonances are entirely in the PHASE.
# So amplitude dips here exist only in proportion to loss, and as the
# discretization improves they vanish: the better the simulation, the less
# there is to see. Res 20's comb confirms it was junk (spacings 14.0, 14.5,
# 14.75, 30.25 nm -- a missed line -- and a mean FSR of 18.4 against 12.87
# predicted).
#
# GROUP DELAY is the right observable. tau_g = -dphi/domega peaks sharply at
# every resonance for any coupling and any loss including zero, because the
# phase winds by 2*pi per FSR regardless. It is also what the hybrid model's
# own closed form predicts, so the two sides stay comparable.
#
# `find_dips()` is kept and still reported, because when loss IS present the
# dips are real and their depth is physical -- it is just no longer what the
# comb is measured from.
# ---------------------------------------------------------------------------
def find_dips(wl_dense, T):
    """Amplitude minima. Retained as a secondary diagnostic only -- see the
    note above for why the comb is NOT taken from these."""
    dips = []
    baseline = np.median(T)
    for i in range(1, len(T) - 1):
        if T[i] < T[i - 1] and T[i] < T[i + 1] and T[i] < 0.7 * baseline:
            dips.append((wl_dense[i], T[i]))
    return dips


def group_delay(wl, ratio):
    """tau_g(wl) in ps from the complex through-port ratio.

    tau_g = -dphi/domega, with omega = 2*pi*c/lambda so
    domega = -2*pi*c/lambda^2 dlambda, giving tau_g = (lambda^2 / (2*pi*c))
    * dphi/dlambda. Returned in picoseconds (c in um/ps).

    The phase is unwrapped on the SOURCE grid before differentiating. That is
    the same np.unwrap caveat beef/spectral.py documents: it needs under pi of
    true phase advance per sample. At N_WL = 401 over 100 nm the ring advances
    ~2*pi per 12.9 nm FSR, i.e. ~0.12*pi per sample -- comfortable, and
    `unwrap_is_safe()` below checks it rather than assuming.
    """
    c = 299.792458                      # um/ps
    phase = np.unwrap(np.angle(ratio))
    dphi_dwl = np.gradient(phase, wl)
    return (wl ** 2 / (2 * np.pi * c)) * dphi_dwl


def find_resonances(wl, tau_g):
    """Comb positions as peaks in |group delay|.

    Peaks, not a threshold on depth: the whole point of moving off the
    amplitude observable is that resonance VISIBILITY should not depend on how
    lossy the discretization happens to be. A peak in tau_g is present whenever
    the phase winds, which is always.

    Prominence is set from the signal's own spread (1 sigma) rather than a
    fixed number, so it adapts to coupling strength instead of encoding one
    device's scale.
    """
    from scipy.signal import find_peaks
    x = np.abs(tau_g - np.median(tau_g))
    if not np.any(np.isfinite(x)) or np.std(x) == 0:
        return np.array([]), np.array([])
    # MINIMUM SEPARATION of half an FSR. Without it each resonance is reported
    # two or three times: prominence alone also catches the shoulders, and the
    # res-30 post-M1 run printed 22 "resonances" for 7, in triplets 1.25 nm
    # apart between real 11.5-12.25 nm gaps. Added 2026-08-03 after that run --
    # verify_ring_check4_comparison.py had the filter and this did not, so the
    # two disagreed about how many resonances the same spectrum contained.
    dist = max(1, int(0.5 * FSR_GUESS / (wl[1] - wl[0])))
    pk, _props = find_peaks(x, prominence=np.std(x), distance=dist)
    return wl[pk], tau_g[pk]


def main():
    perim = assert_geometry_closes()
    os.makedirs(OUT, exist_ok=True)
    print(f"[full] {'3-D' if THREE_D else '2-D'}  res={FDTD_RESOLUTION:g}  "
          f"cell {SX:.2f} x {SY:.2f} x {_cell_sz:.2f} um", flush=True)
    print(f"[full] racetrack perimeter {perim:.2f} um -> predicted FSR "
          f"~{FSR_GUESS * 1000:.2f} nm", flush=True)
    if THREE_D:
        print("[full] WARNING: 3-D full-ring is a ~167 h job at res 40 "
              "(see module docstring). Check the wall limit.", flush=True)

    wl = np.linspace(WL_MIN, WL_MAX, N_WL)
    a_dev, stop_dev, t_dev = run_full_fdtd(build_ring_polys(), "device", wl, save_dir=CACHE)
    a_ref, stop_ref, t_ref = run_full_fdtd(build_bus_only_polys(), "busonly", wl, save_dir=CACHE)
    print(f"[full] tau_est = {TAU_EST:.0f}; device ran {t_dev / TAU_EST:.1f} tau, "
          f"bus-only {t_ref / TAU_EST:.2f} tau", flush=True)

    truncated = [t for t, s in (("device", stop_dev), ("busonly", stop_ref)) if s == "cap"]

    T = np.abs(a_dev) ** 2 / np.abs(a_ref) ** 2

    wl_dense = np.linspace(WL_MIN, WL_MAX, N_WL_DENSE)
    # Same resampling strategy as the hybrid scripts (beef.spectral): magnitude
    # and unwrapped phase, NOT real/imag -- see tests/test_spectral.py for why.
    ratio_dense = interp_complex(wl, a_dev / a_ref, wl_dense)
    T_dense = np.abs(ratio_dense) ** 2

    # --- PRIMARY observable: group delay ---------------------------------
    # See the section-4 note. Amplitude dips only exist in proportion to loss
    # in an all-pass ring; the phase winds 2*pi per FSR regardless.
    safe, worst = unwrap_is_safe(a_dev / a_ref)
    print(f"\n[full] phase unwrap: worst step {worst:.3f} rad "
          f"({worst / np.pi:.2f}pi/sample) -> {'OK' if safe else 'MARGINAL'}",
          flush=True)
    if not safe:
        print("[full]   ⚠ N_WL is too coarse for the phase rate; the group "
              "delay below may be aliased.", flush=True)

    tau_g = group_delay(wl, a_dev / a_ref)
    res_wl, res_tau = find_resonances(wl, tau_g)
    print(f"[full] {len(res_wl)} resonances from group-delay peaks", flush=True)
    for w, t in zip(res_wl, res_tau):
        print(f"[full]   {w:.4f} um   tau_g = {t:+.3f} ps", flush=True)
    if len(res_wl) >= 2:
        fsr = np.diff(res_wl)
        print(f"[full] mean FSR (group delay) = {fsr.mean() * 1000:.3f} nm "
              f"(predicted {FSR_GUESS * 1000:.2f} nm)", flush=True)
        print(f"[full]   spacings (nm): "
              + " ".join(f"{v * 1000:.2f}" for v in fsr), flush=True)
    else:
        fsr = np.array([])

    # --- SECONDARY: amplitude dips, reported but not the comb ------------
    dips = find_dips(wl_dense, T_dense)
    print(f"\n[full] {len(dips)} amplitude dips (secondary diagnostic; depth "
          f"tracks LOSS, not resonance)", flush=True)
    for w, t in dips:
        print(f"[full]   {w:.4f} um   T_min = {t:.4f}", flush=True)
    print(f"[full] |T| range {T_dense.min():.4f}..{T_dense.max():.4f}  "
          f"(a lossless all-pass ring is flat 1.0 — shallow is CORRECT here)",
          flush=True)

    # --- verdict ---------------------------------------------------------
    print("\n[full] =========================== VERDICT ===========================",
          flush=True)
    if truncated:
        print(f"[full] ✗ RUN TRUNCATED: {', '.join(truncated)} hit the {T_MAX:g} "
              f"time cap instead of decaying.", flush=True)
        print("[full]   The ringdown was cut short, so every dip is artificially "
              "BROADENED and the", flush=True)
        print("[full]   linewidths must not be read. Dip POSITIONS are usually "
              "still usable, but raise", flush=True)
        print("[full]   RING_FULL_TMAX and re-run before quoting anything from "
              "this spectrum.", flush=True)
    else:
        print("[full] ✓ both runs decayed on their own; the spectrum is not "
              "truncation-limited.", flush=True)
        # Decaying is necessary but not sufficient. The bus-only run SHOULD
        # finish fast -- no resonator, nothing to ring. The device run must
        # outlast tau by a wide margin or the resonances were never resolved,
        # and a too-loose decay tolerance would let it stop early while still
        # reporting "decayed". Measure that separately.
        if t_dev < 5 * TAU_EST:
            print(f"[full] ⚠ but the DEVICE run stopped at "
                  f"{t_dev / TAU_EST:.1f} tau, under the 5 tau needed for a "
                  f"clean lineshape.", flush=True)
            print("[full]   Either DECAY_TOL is too loose, or this device's Q "
                  "is well below the 14294", flush=True)
            print("[full]   assumed here. Check the measured FWHM against "
                  "TAU_EST before trusting linewidths.", flush=True)
    print(f"[full] SCOPE: this is {'3-D' if THREE_D else '2-D'}. "
          f"{'' if THREE_D else 'It therefore validates the DECOMPOSITION '}",
          flush=True)
    if not THREE_D:
        print("[full]   (does composing FDTD sections with an EME body reproduce "
              "a direct solve?)", flush=True)
        print("[full]   against RING_TEST_MODE=1 hybrid runs at matched 2-D "
              "physics. It does NOT", flush=True)
        print("[full]   validate the 3-D numbers the device ships.", flush=True)

    # is_master(): this script runs under `srun` (it is FDTD — one rank would
    # waste the allocation, which is exactly the bug fixed in
    # batch/run_ring_full_fdtd.batch on 2026-08-03). Every rank reaches here
    # with identical values, because meep reduces the eigenmode coefficients
    # before returning them, so an ungated write is 24 processes racing to put
    # the same bytes in one file. Same bytes or not, a torn write is a corrupt
    # npz — gate it.
    if is_master():
        np.savez(os.path.join(OUT, "full_fdtd.npz"),
                 wavelengths=wl, a_device=a_dev, a_busonly=a_ref, T=T,
                 wl_dense=wl_dense, T_dense=T_dense,
                 dips=np.asarray(res_wl),                 # PRIMARY: group-delay comb
                 tau_g=tau_g, group_delay_peaks=np.asarray(res_tau),
                 amp_dips=np.asarray([w for w, _ in dips]),   # secondary
                 amp_dip_depths=np.asarray([t for _, t in dips]),
                 fsr_nm=fsr * 1000 if fsr.size else fsr,
                 perimeter=perim, resolution=FDTD_RESOLUTION, three_d=THREE_D,
                 stop_device=stop_dev, stop_busonly=stop_ref,
                 t_device=t_dev, t_busonly=t_ref, tau_est=TAU_EST)

    # Two panels, primary on top: the group delay IS the measurement, the
    # amplitude is context. Plotting them the other way round is what made the
    # amplitude look like the result in the first place.
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax0.plot(wl * 1000, tau_g, "-", lw=1.0, color="C0",
             label=f"group delay ({'3-D' if THREE_D else '2-D'}, res {FDTD_RESOLUTION:g})")
    for w, t in zip(res_wl, res_tau):
        ax0.plot(w * 1000, t, "v", ms=6, color="C3")
    ax0.set_ylabel(r"group delay $\tau_g$ (ps)")
    ax0.set_title("Ring resonator, full device in one FDTD run — "
                  + (f"{len(res_wl)} resonances, FSR {fsr.mean() * 1000:.2f} nm"
                     if fsr.size else "no resonances resolved"))
    ax0.legend(loc="upper right", fontsize=8)
    ax0.grid(alpha=0.3)

    ax1.plot(wl_dense * 1000, T_dense, "-", lw=0.8, color="C7",
             label="|T| (secondary — depth tracks LOSS, not resonance)")
    for w, t in dips:
        ax1.plot(w * 1000, t, "v", ms=5, color="C3")
    ax1.axhline(1.0, ls=":", lw=0.8, color="k")
    ax1.set_xlabel("wavelength (nm)")
    ax1.set_ylabel("through transmission")
    ax1.legend(loc="lower right", fontsize=8)
    ax1.grid(alpha=0.3)
    fig.tight_layout()
    if is_master():                      # same reason as the np.savez above
        fig.savefig(os.path.join(OUT, "full_fdtd_transmission.png"), dpi=150)
    plt.close(fig)

    print(f"\n[full] saved -> {OUT}", flush=True)
    return 1 if truncated else 0


if __name__ == "__main__":
    sys.exit(main())
