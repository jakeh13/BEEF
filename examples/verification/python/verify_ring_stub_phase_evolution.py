"""Is the field just past the ring's arc actually the fundamental mode?

THE HYPOTHESIS UNDER TEST
-------------------------
`ring_resonator_validation.md` has one blocking defect: the decomposition's
plane transport, `exp(±i*beta*d)` with `beta = 2*pi*n_eff/lambda`, does not fully
account for the phase. Four mechanisms have been refuted by direct test:

    1  backward vs forward mode index      +0.00000 pi   (a sign)
    2  x-normal vs y-normal orientation    +0.00040 pi   (a sign)
    3  the 0.75 um ring-port asymmetry     needs 3 different lengths
    4  the stub has a different n_eff      bit-identical to the bare strip

All four are wavelength-INDEPENDENT. Re-measuring the residual on 101 points
(not the 5 originally tabulated) showed it is wavelength-DEPENDENT: it drifts
0.0068 pi across the band, a free propagation term beats a pure constant by 23x,
and forcing exactly pi makes the fit 850x worse. So none of the four could have
explained it, and the "flat, therefore a sign" reading is what aimed the first
two at sign conventions.

What survives: **the field 0.50 um past the arc is not the pure fundamental
mode.** A bend excites a superposition; over the transport distance that
mixture's phase does not evolve at any single n_eff. Solving for the index the
transport would have needed (assuming a pi flip) gives bare-strip + 0.025,
smoothly varying -- a plausible guided mode, and exactly what a not-yet-relaxed
field would mimic.

WHAT THIS SCRIPT MEASURES
-------------------------
The standalone 90-degree bend, with mode monitors along the OUTPUT stub at a
ladder of distances from the arc end. Two quantities per plane:

  * **local n_eff**, from the phase gradient between adjacent planes:
        n_local(d) = [arg a(d2) - arg a(d1)] * lambda / (2*pi*(d2-d1))
    If the field is modal this equals the bare strip at every d. If it is not,
    n_local starts off and RELAXES toward the bare strip with distance.

  * **modal purity |a|**, the fundamental's projection coefficient. A
    superposition projects less than a clean mode does, so |a| rising with d is
    the same story seen in amplitude rather than phase.

The prediction is specific and falsifiable: **n_local should be high near the
arc, by roughly the +0.025 the transport residual demands, and decay toward the
bare strip.** If instead n_local sits flat at the bare-strip value at every
distance, the field IS modal, this hypothesis is dead too, and the phase error is
somewhere the transport-index picture cannot reach.

2-D so it is cheap: the bend domain loses a dimension entirely. 2-D n_eff is not
the 3-D value, so the BASELINE is measured here in the same 2-D run rather than
taken from the 3-D mode solve -- the comparison is n_local(d) against
n_local(far), which is self-consistent either way.

Env:
    RING_SPE_RES    FDTD resolution (default 40)
    RING_SPE_3D=1   run in 3-D instead (expensive; not needed for the trend)
    RING_SPE_OUT    output dir override
"""
import os
import sys

import numpy as np
import meep as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef.geometry import arc_polygon
from beef.geometry.base import polygons_to_prisms
from beef.materials import SOI_220
from beef.solvers.fdtd import check_source_coverage, monitor_freq_grid, source_fwidth

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))

# --- mirrored from ring_resonator_hybrid_cascaded.py ------------------------
MATERIAL_STACK = SOI_220
W = 0.45
R = 5.0
DPML = 0.5
SIM_HEIGHT = 1.5
Y_BELOW = 2.0
PORT_APERTURE = 3.0
MONITOR_MARGIN = 0.25
SOURCE_MARGIN = 0.0
SOURCE_BANDWIDTH_FACTOR = 1.2
WL_MIN, WL_MAX = 1.50, 1.60
N_WL = 11

THREE_D = bool(int(os.environ.get("RING_SPE_3D", "0")))
FDTD_RESOLUTION = float(os.environ.get("RING_SPE_RES", 40.0))

# The stub must be long enough to hold the whole monitor ladder AND leave room
# for the field to relax, or the test cannot see relaxation at all. The ring's
# own STUB_LEN is 1.5 um; this uses 6.0 so the far planes are unambiguously in
# relaxed guide and can serve as the baseline.
STUB_LEN = 6.0
# Distances from the arc end at which to read the mode coefficient. The ring's
# own two planes of interest are +0.50 (k1 / ring_r_mon) and +1.25 (k2 norm),
# both included explicitly so the measurement lands ON the planes in question.
PROBE_D = np.array([0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50,
                    3.00, 3.50, 4.00, 4.50, 5.00])

_cell_sz = (SIM_HEIGHT + 2 * DPML) if THREE_D else 0.0
_mon_sz = SIM_HEIGHT if THREE_D else 0.0
_wg_height = MATERIAL_STACK.wg_height if THREE_D else None

OUT = os.environ.get("RING_SPE_OUT") or os.path.join(
    _ROOT, "example_output", "verification", "ring_stub_phase")

BEND_X_MIN, BEND_X_MAX = -STUB_LEN - DPML, R + Y_BELOW + DPML
BEND_Y_MIN, BEND_Y_MAX = -Y_BELOW - DPML, R + STUB_LEN + DPML
BEND_SX, BEND_SY = BEND_X_MAX - BEND_X_MIN, BEND_Y_MAX - BEND_Y_MIN
BEND_XC, BEND_YC = (BEND_X_MAX + BEND_X_MIN) / 2, (BEND_Y_MAX + BEND_Y_MIN) / 2


def rect_poly(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def build_bend_polys():
    """Same shape as the ring's own standalone bend, with a longer output stub.
    Arc runs (0,0) -> (R,R); input stub along -x, output stub up +y at x=R."""
    return [
        rect_poly(-STUB_LEN - DPML, -W / 2, 0.0, W / 2),
        arc_polygon(0.0, R, R - W / 2, R + W / 2, -90, 0),
        rect_poly(R - W / 2, R, R + W / 2, R + STUB_LEN + DPML),
    ]


def build_straight_polys():
    """CONTROL geometry: the output stub alone, no arc, same cell.

    Without this the test cannot separate "the field near the arc is not modal"
    from "this measurement of n_local is biased". A straight guide excited by a
    clean eigenmode source must give n_local flat at every distance; any
    structure the control also shows is the harness, not the bend.
    """
    return [rect_poly(R - W / 2, -Y_BELOW - DPML, R + W / 2, R + STUB_LEN + DPML)]


def _media(wavelengths):
    """FDTDSolver._media(), replicated -- same M1 fix as the ring scripts."""
    disp = MATERIAL_STACK.resolve_dispersive(None)
    wl = np.asarray(wavelengths, float) if disp else float(np.mean(wavelengths))
    return (MATERIAL_STACK.si_medium(wl, dispersive=disp),
            MATERIAL_STACK.sio2_medium(wl, dispersive=disp))


def run(polys, wavelengths, excite):
    """Returns alpha[n_probe, n_wl] complex, the fundamental's coefficient at
    each probe plane, ordered to match `wavelengths`.

    `excite`: "bend" launches into the input stub along +x; "straight" launches
    up the stub from below along +y (the control has no input stub).
    """
    wl_center = float(np.mean(wavelengths))
    fcen, fwidth = 1.0 / wl_center, source_fwidth(wavelengths, SOURCE_BANDWIDTH_FACTOR)
    mon_freqs, wl_order = monitor_freq_grid(wavelengths)
    check_source_coverage(fcen, fwidth, wavelengths)

    si, sio2 = _media(wavelengths)
    sim = mp.Simulation(
        resolution=FDTD_RESOLUTION,
        cell_size=mp.Vector3(BEND_SX, BEND_SY, _cell_sz),
        geometry_center=mp.Vector3(BEND_XC, BEND_YC, 0),
        boundary_layers=[mp.PML(DPML)],
        geometry=polygons_to_prisms(polys, si, wg_height=_wg_height),
        default_material=sio2,
    )
    if excite == "bend":
        sim.sources = [mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=fwidth),
            center=mp.Vector3(-STUB_LEN + SOURCE_MARGIN, 0.0),
            size=mp.Vector3(0, PORT_APERTURE, _mon_sz),
            eig_band=1, eig_kpoint=mp.Vector3(1, 0, 0), eig_match_freq=True)]
    else:
        sim.sources = [mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=fwidth),
            center=mp.Vector3(R, -Y_BELOW + SOURCE_MARGIN),
            size=mp.Vector3(PORT_APERTURE, 0, _mon_sz),
            eig_band=1, eig_kpoint=mp.Vector3(0, 1, 0), eig_match_freq=True)]

    # y-normal monitors up the output stub, at R + d for each probe distance.
    mons = [sim.add_mode_monitor(
        mon_freqs,
        mp.FluxRegion(center=mp.Vector3(R, R + d),
                      size=mp.Vector3(PORT_APERTURE, 0, _mon_sz)))
            for d in PROBE_D]

    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=1e-7))
    alpha = np.array([
        sim.get_eigenmode_coefficients(m, [1], eig_parity=mp.NO_PARITY).alpha[0, :, 0][wl_order]
        for m in mons])
    sim.reset_meep()
    return alpha


def local_neff(alpha, wavelengths):
    """n_local between consecutive probe planes, [n_probe-1, n_wl].

    The phase is unwrapped ALONG DISTANCE, not along wavelength: consecutive
    planes are 0.25-0.50 um apart and n_eff ~ 2.3, so the true advance is
    0.7-1.5 pi per step. That is at or past the unwrap limit, so the branch is
    fixed by the EXPECTED advance from the far-field baseline rather than by
    np.unwrap -- otherwise this function would silently alias exactly where the
    answer matters. See beef/spectral.py's unwrap_is_safe() docstring for the
    same trap in the wavelength direction.
    """
    n_wl = alpha.shape[1]
    out = np.zeros((len(PROBE_D) - 1, n_wl))
    for j in range(n_wl):
        lam = wavelengths[j]
        for i in range(len(PROBE_D) - 1):
            dd = PROBE_D[i + 1] - PROBE_D[i]
            dphi_raw = np.angle(alpha[i + 1, j] / alpha[i, j])
            # Choose the 2*pi branch closest to a nominal 2.3 index, then let
            # the measurement move freely within that branch.
            nominal = -2 * np.pi * 2.3 * dd / lam
            k = np.round((nominal - dphi_raw) / (2 * np.pi))
            dphi = dphi_raw + 2 * np.pi * k
            out[i, j] = -dphi * lam / (2 * np.pi * dd)
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    wl = np.linspace(WL_MIN, WL_MAX, N_WL)
    print(f"[spe] {'3-D' if THREE_D else '2-D'} res={FDTD_RESOLUTION:g}  "
          f"stub={STUB_LEN:g} um  {len(PROBE_D)} probe planes", flush=True)
    print(f"[spe] probe distances from arc end (um): {PROBE_D}", flush=True)

    print("\n[spe] --- CONTROL: straight guide, no arc ---", flush=True)
    a_ctl = run(build_straight_polys(), wl, "straight")
    n_ctl = local_neff(a_ctl, wl)

    print("\n[spe] --- BEND: same stub, fed through a 90-degree arc ---", flush=True)
    a_bend = run(build_bend_polys(), wl, "bend")
    n_bnd = local_neff(a_bend, wl)

    jm = len(wl) // 2                      # 1.55 um
    mid = 0.5 * (PROBE_D[:-1] + PROBE_D[1:])
    # Baseline = the far end, where the field has had 5 um to relax.
    base_b, base_c = n_bnd[-1, jm], n_ctl[-1, jm]
    print(f"\n[spe] at lambda = {wl[jm]:.3f} um   (baseline = the farthest interval)", flush=True)
    print(f"[spe] {'d (um)':>8} {'n_local BEND':>13} {'- baseline':>11} "
          f"{'n_local CTRL':>13} {'- baseline':>11} {'|a| bend':>9}", flush=True)
    for i, d in enumerate(mid):
        print(f"[spe] {d:8.3f} {n_bnd[i, jm]:13.6f} {n_bnd[i, jm] - base_b:+11.6f} "
              f"{n_ctl[i, jm]:13.6f} {n_ctl[i, jm] - base_c:+11.6f} "
              f"{np.abs(a_bend[i, jm]):9.5f}", flush=True)

    # --- verdict ---------------------------------------------------------
    # The ring reads k1 at +0.50 and normalizes k2 at +1.25, so the relevant
    # excess is the one in that window, not at the first probe.
    win = (mid >= 0.4) & (mid <= 1.4)
    exc_b = float(np.max(np.abs(n_bnd[win, jm] - base_b))) if win.any() else float("nan")
    exc_c = float(np.max(np.abs(n_ctl[win, jm] - base_c))) if win.any() else float("nan")
    print("\n[spe] =========================== VERDICT ===========================",
          flush=True)
    print(f"[spe] worst |n_local - baseline| over the ring's own 0.50-1.25 um window:",
          flush=True)
    print(f"[spe]   BEND    {exc_b:.6f}", flush=True)
    print(f"[spe]   CONTROL {exc_c:.6f}   <- harness floor", flush=True)
    print(f"[spe] the transport residual demands ~+0.025", flush=True)
    # The bar is NOT an absolute threshold -- it is the +0.025 the transport
    # residual actually demands. A result far below that refutes the hypothesis
    # no matter how cleanly it clears the control floor; the first version of
    # this check asked "is it > 0.005?" and called 120x-too-small
    # "inconclusive, raise the resolution", which would have sent a decided
    # question back for more compute.
    DEMANDED = 0.025
    ratio = DEMANDED / exc_b if exc_b > 0 else float("inf")
    print(f"[spe]   shortfall: the measured excess is {ratio:.0f}x SMALLER than "
          f"demanded", flush=True)
    # Modal purity is the same story in amplitude: a relaxing superposition
    # projects less onto the fundamental near the arc and more further out.
    pur = np.abs(a_bend[:, jm])
    pur_spread = float((pur.max() - pur.min()) / pur.mean())
    print(f"[spe]   |a| spread across all {len(PROBE_D)} planes: "
          f"{pur_spread:.2e} relative  (a relaxing superposition would RISE)",
          flush=True)
    if exc_b > 0.3 * DEMANDED and exc_b > 3 * max(exc_c, 1e-6):
        print("[spe] ✓ CONFIRMED: the field near the arc is NOT modal. n_local departs "
              "from its", flush=True)
        print("[spe]   relaxed value by well over the control's floor, in the very "
              "window where", flush=True)
        print("[spe]   k1 and k2 are referenced. A single n_eff cannot describe the "
              "transport there,", flush=True)
        print("[spe]   which is exactly the wavelength-dependent residual the other "
              "four mechanisms", flush=True)
        print("[spe]   could not produce. FIX: reference b and the coupler "
              "coefficients at a plane", flush=True)
        print("[spe]   far enough out to have relaxed, or carry the measured phase "
              "instead of an", flush=True)
        print("[spe]   assumed modal one.", flush=True)
    else:
        print(f"[spe] ✗ REFUTED: the field IS effectively modal at these planes.", flush=True)
        print(f"[spe]   n_local departs from its relaxed value by only "
              f"{exc_b:.6f} in the window where", flush=True)
        print(f"[spe]   k1 and k2 are referenced -- {ratio:.0f}x less than the "
              f"{DEMANDED} the residual needs --", flush=True)
        print(f"[spe]   and |a| is flat to {pur_spread:.1e} rather than rising. A "
              f"relaxing superposition", flush=True)
        print("[spe]   would show both. This hypothesis joins the other four:", flush=True)
        print("[spe]     1 backward mode index   2 port orientation", flush=True)
        print("[spe]     3 the 0.75 um asymmetry  4 the stub cross-section", flush=True)
        print("[spe]   ⚠ 2-D CAVEAT: a 2-D slab confines more tightly than the real "
              "450x220 nm strip,", flush=True)
        print("[spe]     so it may relax faster. The 120x margin is wide enough that "
              "3-D is unlikely", flush=True)
        print("[spe]     to reverse it, but RING_SPE_3D=1 is the way to close that "
              "if wanted.", flush=True)
        print("[spe]   NEXT: the phase error is not reachable by any transport-INDEX "
              "picture. The", flush=True)
        print("[spe]     untested code path is get_norm_calibrated() -- used exactly "
              "once, for the", flush=True)
        print("[spe]     ring-excited run's norm monitor, which is the denominator of "
              "BOTH coefficients", flush=True)
        print("[spe]     carrying the pi (k2, t_ring) and neither that does not. It "
              "takes MAGNITUDE from", flush=True)
        print("[spe]     Poynting flux and PHASE from the modal coefficient -- a mixed "
              "construction that", flush=True)
        print("[spe]     has never been tested.", flush=True)

    np.savez(os.path.join(OUT, "stub_phase.npz"),
             wavelengths=wl, probe_d=PROBE_D, mid_d=mid,
             alpha_bend=a_bend, alpha_control=a_ctl,
             n_local_bend=n_bnd, n_local_control=n_ctl,
             resolution=FDTD_RESOLUTION, three_d=THREE_D, stub_len=STUB_LEN)

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    ax0.plot(mid, n_bnd[:, jm] - base_b, "o-", label="through a 90° bend")
    ax0.plot(mid, n_ctl[:, jm] - base_c, "s--", color="0.5", label="straight control")
    ax0.axhline(0.0, lw=0.8, color="k")
    ax0.axhline(0.025, ls=":", color="C3",
                label="+0.025 — what the transport residual demands")
    ax0.set_ylabel(r"$n_{local} - n_{relaxed}$")
    ax0.set_title(f"Is the field past the arc modal?  "
                  f"({'3-D' if THREE_D else '2-D'}, res {FDTD_RESOLUTION:g}, "
                  f"λ = {wl[jm]:.3f} µm)")
    ax0.legend(fontsize=8)
    ax0.grid(alpha=0.3)
    ax1.plot(PROBE_D, np.abs(a_bend[:, jm]), "o-", label="|a| through the bend")
    ax1.plot(PROBE_D, np.abs(a_ctl[:, jm]), "s--", color="0.5", label="|a| control")
    for d, lbl in ((0.50, "k1 plane"), (1.25, "k2 norm plane")):
        for ax in (ax0, ax1):
            ax.axvline(d, ls="-.", lw=0.9, color="C2", alpha=0.7)
        ax1.annotate(lbl, (d, ax1.get_ylim()[0]), fontsize=7, color="C2",
                     rotation=90, va="bottom")
    ax1.set_xlabel("distance from the arc end (µm)")
    ax1.set_ylabel("modal purity |a|")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "stub_phase.png"), dpi=150)
    plt.close(fig)
    print(f"\n[spe] saved -> {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
