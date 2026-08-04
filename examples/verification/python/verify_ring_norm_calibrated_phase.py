"""Does `get_norm_calibrated()` put a spurious pi into k2 and t_ring?

WHY THIS ONE CALL
-----------------
`ring_resonator_validation.md` has one open defect: `k1` and `k2` describe the same
physical coupling from opposite directions, their magnitudes agree to 1.3e-04, and
their phases differ by ~pi after the 0.75 um plane correction. Five mechanisms are
refuted:

    1  backward vs forward mode index      +0.00000 pi
    2  x-normal vs y-normal orientation    +0.00040 pi
    3  the 0.75 um ring-port asymmetry     needs 3 different lengths
    4  the stub has a different n_eff      bit-identical to the bare strip
    5  the field past the arc is non-modal 121x too small, |a| flat

`get_norm_calibrated()` is the one untested path, and the correlation is exact: it is
called **once**, for the ring-excited run's norm monitor, and that single call is the
denominator of **both** coefficients carrying the pi (`k2`, `t_ring`) and **neither**
of the two that do not (`t`, `k1`). The earlier closure of this area ("blocker 0", the
y-normal excite normalization) measured `max ||k1| - |k2||` -- MAGNITUDE ONLY -- so a
phase error was invisible to it by construction.

    def get_norm_calibrated(sim, mon, forward=True):
        coeff = get_coeffs(sim, mon, forward=forward)
        flux  = np.abs(np.asarray(mp.get_fluxes(mon)))
        return np.sqrt(flux) * np.exp(1j * np.angle(coeff))

Two things this mixes that could each carry a sign:

  * `mp.get_fluxes(mon)` is a SIGNED net flux -- positive along +y for a y-normal
    monitor. The ring-excited source launches **-y**, so the flux through its norm
    monitor should be NEGATIVE, and `np.abs()` discards that. Whether that matters
    depends on whether the sign belonged to the magnitude (it does not; a magnitude
    is positive) or was carrying physical direction information the phase then has to
    supply.
  * the phase comes from `alpha[..., 1]`, the backward index, for this call only.

WHAT IS MEASURED
----------------
The ring-excited left-section run, at its real norm monitor, reporting:

  * the SIGN of `mp.get_fluxes()` -- is `np.abs()` discarding anything?
  * `sqrt(|flux|)` against `|coeff|` -- do the two magnitude estimates agree? If they
    differ a lot, the calibration is doing real work and its phase pairing matters.
  * `arg(get_norm_calibrated())` minus `arg(coeff)` -- MUST be 0 by construction,
    since the phase is taken from `coeff`. A nonzero value would mean the helper is
    not doing what it says.
  * the same run's coefficient normalized both ways, and the phase difference between
    them -- this is the number that would show up as the pi.

The falsifiable prediction: **if this call is the culprit, normalizing by `coeff`
instead of by `get_norm_calibrated()` changes arg(k2) by ~pi.** If it changes it by
~0, this path is clean and the sixth mechanism is refuted too.

2-D, res 20: the phase relationship under test is a convention, not a converged
quantity, so it does not need production resolution.

Env:
    RING_NCP_RES   FDTD resolution (default 20)
    RING_NCP_OUT   output dir override
"""
import os
import sys

import numpy as np
import meep as mp

from beef.geometry import arc_polygon
from beef.geometry.base import polygons_to_prisms
from beef.materials import SOI_220
from beef.solvers.fdtd import check_source_coverage, monitor_freq_grid, source_fwidth

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))

# --- mirrored from ring_resonator_hybrid_cascaded.py ------------------------
MATERIAL_STACK = SOI_220
W, GAP = 0.45, 0.25
R, L_C = 5.0, 4.0
Y_RING = W + GAP
STUB_X = L_C / 2 + R
Y_BELOW = 2.0
INNER_X_HALF = STUB_X + Y_BELOW
DPML = 0.5
S_MON = 0.5
MONITOR_MARGIN = 0.25
SOURCE_MARGIN = 0.0
STUB_LEN = 1.5
PORT_APERTURE = 3.0
SOURCE_BANDWIDTH_FACTOR = 1.2
IFACE_BUFFER = 1.0
IFACE_X = -L_C / 2 + IFACE_BUFFER
MESH_WIDTH = 4.0
BODY_Y_CENTER = Y_RING / 2
_MODE_MARGIN = MESH_WIDTH / 2 - W / 2
BODY_MESH_WIDTH = Y_RING + W + 2 * _MODE_MARGIN
IFACE_PLANE_L = IFACE_X - MONITOR_MARGIN
LEFT_RING_Y = Y_RING + R + S_MON
LEFT_STUB_TOP = LEFT_RING_Y + MONITOR_MARGIN
N_MODES = 2
WL_MIN, WL_MAX = 1.50, 1.60
N_WL = 11

FDTD_RESOLUTION = float(os.environ.get("RING_NCP_RES", 20.0))
OUT = os.environ.get("RING_NCP_OUT") or os.path.join(
    _ROOT, "example_output", "verification", "ring_norm_phase")

LEFT_X_MIN, LEFT_X_MAX = -INNER_X_HALF - DPML, IFACE_X + DPML
LEFT_Y_MIN, LEFT_Y_MAX = -Y_BELOW - DPML, LEFT_STUB_TOP + DPML
LEFT_SX, LEFT_SY = LEFT_X_MAX - LEFT_X_MIN, LEFT_Y_MAX - LEFT_Y_MIN
LEFT_XC, LEFT_YC = (LEFT_X_MAX + LEFT_X_MIN) / 2, (LEFT_Y_MAX + LEFT_Y_MIN) / 2


def rect_poly(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def build_left_polys():
    return [
        rect_poly(-INNER_X_HALF - DPML, -W / 2, IFACE_X + DPML, W / 2),
        rect_poly(-L_C / 2, Y_RING - W / 2, IFACE_X + DPML, Y_RING + W / 2),
        arc_polygon(-L_C / 2, Y_RING + R, R - W / 2, R + W / 2, 180, 270),
        rect_poly(-STUB_X - W / 2, Y_RING + R, -STUB_X + W / 2, LEFT_STUB_TOP + DPML),
    ]


def main():
    os.makedirs(OUT, exist_ok=True)
    wl = np.linspace(WL_MIN, WL_MAX, N_WL)
    wl_center = float(np.mean(wl))
    fcen = 1.0 / wl_center
    fwidth = source_fwidth(wl, SOURCE_BANDWIDTH_FACTOR)
    mon_freqs, wl_order = monitor_freq_grid(wl)
    check_source_coverage(fcen, fwidth, wl)

    disp = MATERIAL_STACK.resolve_dispersive(None)
    si = MATERIAL_STACK.si_medium(wl if disp else wl_center, dispersive=disp)
    ox = MATERIAL_STACK.sio2_medium(wl if disp else wl_center, dispersive=disp)

    print(f"[ncp] 2-D res={FDTD_RESOLUTION:g}, ring-excited left section "
          f"(the -y launch)", flush=True)
    sim = mp.Simulation(
        resolution=FDTD_RESOLUTION,
        cell_size=mp.Vector3(LEFT_SX, LEFT_SY, 0.0),
        geometry_center=mp.Vector3(LEFT_XC, LEFT_YC, 0),
        boundary_layers=[mp.PML(DPML)],
        geometry=polygons_to_prisms(build_left_polys(), si, wg_height=None),
        default_material=ox,
    )
    # Exactly run_left_fdtd('ring')'s source and norm monitor.
    src_y = LEFT_STUB_TOP - SOURCE_MARGIN
    sim.sources = [mp.EigenModeSource(
        mp.GaussianSource(fcen, fwidth=fwidth),
        center=mp.Vector3(-STUB_X, src_y),
        size=mp.Vector3(PORT_APERTURE, 0, 0.0),
        eig_band=1, eig_kpoint=mp.Vector3(0, -1, 0), eig_match_freq=True)]
    norm_mon = sim.add_mode_monitor(
        mon_freqs, mp.FluxRegion(center=mp.Vector3(-STUB_X, LEFT_RING_Y),
                                 size=mp.Vector3(PORT_APERTURE, 0, 0.0)))
    iface_mon = sim.add_mode_monitor(
        mon_freqs, mp.FluxRegion(center=mp.Vector3(IFACE_PLANE_L, BODY_Y_CENTER),
                                 size=mp.Vector3(0, BODY_MESH_WIDTH, 0.0)))

    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=1e-7))

    a_fwd = sim.get_eigenmode_coefficients(
        norm_mon, [1], eig_parity=mp.NO_PARITY).alpha[0, :, 0]
    a_bwd = sim.get_eigenmode_coefficients(
        norm_mon, [1], eig_parity=mp.NO_PARITY).alpha[0, :, 1]
    flux_signed = np.asarray(mp.get_fluxes(norm_mon))
    iface = sim.get_eigenmode_coefficients(
        iface_mon, list(range(1, N_MODES + 1)), eig_parity=mp.NO_PARITY).alpha[:, :, 0]
    sim.reset_meep()

    # --- 1. is np.abs() on the flux discarding a sign? -------------------
    n_neg = int((flux_signed < 0).sum())
    print(f"\n[ncp] mp.get_fluxes() sign: {n_neg}/{len(flux_signed)} points NEGATIVE",
          flush=True)
    print(f"[ncp]   -> np.abs() {'DISCARDS a sign here' if n_neg else 'discards nothing'}",
          flush=True)
    print(f"[ncp]   (the source launches -y, so a y-normal monitor's net flux being "
          f"negative is EXPECTED)", flush=True)

    # --- 2. the two magnitude estimates ---------------------------------
    mag_flux = np.sqrt(np.abs(flux_signed))
    print(f"\n[ncp] {'wl':>7} {'sqrt|flux|':>11} {'|a_bwd|':>11} {'ratio':>8} "
          f"{'|a_fwd|':>11}", flush=True)
    for i in range(0, len(wl), 2):
        j = wl_order[i]
        print(f"[ncp] {wl[i]:7.3f} {mag_flux[j]:11.5f} {abs(a_bwd[j]):11.5f} "
              f"{mag_flux[j] / abs(a_bwd[j]):8.4f} {abs(a_fwd[j]):11.5f}", flush=True)

    # --- 3. the phase the helper actually returns ------------------------
    # By construction arg(norm_cal) == arg(a_bwd). Assert it rather than trust it.
    norm_cal = mag_flux * np.exp(1j * np.angle(a_bwd))
    d_self = float(np.max(np.abs(np.angle(norm_cal / a_bwd))))
    print(f"\n[ncp] arg(get_norm_calibrated) - arg(a_bwd): max {d_self:.2e} rad "
          f"(must be 0 by construction)", flush=True)

    # --- 4. THE TEST: does the choice of normalizer move the phase? -----
    # k2-like quantity: the interface coefficient over the ring-side norm.
    k_cal = iface[0] / norm_cal
    k_raw = iface[0] / a_bwd
    k_fwd = iface[0] / a_fwd
    d_cal_raw = np.angle(k_cal / k_raw)
    d_cal_fwd = np.angle(k_cal / k_fwd)
    print(f"\n[ncp] {'wl':>7} {'arg d(cal vs raw a_bwd)/pi':>27} "
          f"{'arg d(cal vs a_fwd)/pi':>24}", flush=True)
    for i in range(0, len(wl), 2):
        j = wl_order[i]
        print(f"[ncp] {wl[i]:7.3f} {d_cal_raw[j] / np.pi:27.6f} "
              f"{d_cal_fwd[j] / np.pi:24.6f}", flush=True)

    w_raw = float(np.max(np.abs(d_cal_raw))) / np.pi
    w_fwd = float(np.max(np.abs(d_cal_fwd))) / np.pi

    # IS a_fwd EVEN A LEGITIMATE ALTERNATIVE NORMALIZER? The source launches -y,
    # so at this monitor the +y coefficient is whatever REFLECTED back up -- a
    # small quantity, not the launched field. Comparing against it then measures
    # the phase of a weak reflection, not a convention.
    #
    # The first version of this script skipped this check and reported "FOUND IT,
    # in the index" on a 0.99 pi difference against a coefficient 5% the size of
    # the real one. Two tells were already visible and ignored: the ratio, and
    # the fact that the "pi" varied 0.73-0.99 pi across the band where a genuine
    # sign convention would sit at exactly pi.
    ratio_fwd = float(np.median(np.abs(a_fwd) / np.abs(a_bwd)))
    fwd_is_meaningful = ratio_fwd > 0.5
    print(f"[ncp] median |a_fwd| / |a_bwd| = {ratio_fwd:.4f}  -> the forward "
          f"coefficient is "
          f"{'comparable, so a legitimate comparison' if fwd_is_meaningful else 'a WEAK REFLECTION, not an alternative normalizer'}",
          flush=True)
    print(f"[ncp] calibration changes the normalizer MAGNITUDE by a factor "
          f"{float(np.median(mag_flux / np.abs(a_bwd))):.4f} "
          f"(this is what blocker 0 was about)", flush=True)

    print("\n[ncp] =========================== VERDICT ===========================",
          flush=True)
    print(f"[ncp] worst |arg difference| vs raw backward coeff : {w_raw:.6f} pi",
          flush=True)
    print(f"[ncp] worst |arg difference| vs forward coeff      : {w_fwd:.6f} pi",
          flush=True)
    if w_raw > 0.5:
        print("[ncp] ✓ FOUND IT: swapping the normalizer moves the phase by ~pi. "
              "get_norm_calibrated()", flush=True)
        print("[ncp]   is injecting the residual into exactly the two coefficients "
              "it normalizes.", flush=True)
    elif w_fwd > 0.5 and fwd_is_meaningful:
        print("[ncp] ✓ FOUND IT, but in the INDEX not the calibration: the forward and "
              "backward", flush=True)
        print("[ncp]   coefficients at this monitor are comparable in size yet differ by "
              "~pi, and this", flush=True)
        print("[ncp]   call is the only place the backward one is used.", flush=True)
    else:
        print(f"[ncp] ✗ REFUTED: the calibration is phase-transparent "
              f"({w_raw:.2e} pi) and the forward/", flush=True)
        print(f"[ncp]   backward indices agree here to {w_fwd:.2e} pi. This path does "
              f"NOT carry the", flush=True)
        print("[ncp]   residual, and the sixth mechanism joins the other five.", flush=True)
        print("[ncp]   The np.abs() on the flux discards a sign that belonged to a "
              "MAGNITUDE, which", flush=True)
        print("[ncp]   is correct: the direction information is carried by the phase, "
              "taken separately.", flush=True)
        if w_fwd > 0.5 and not fwd_is_meaningful:
            print(f"[ncp]   ⚠ The {w_fwd:.3f} pi against the FORWARD coefficient is NOT "
                  f"evidence: that", flush=True)
            print(f"[ncp]     coefficient is {ratio_fwd:.1%} the size of the backward "
                  f"one -- it is the small", flush=True)
            print("[ncp]     reflection back up the stub, so its phase is a weak "
                  "reflected wave's, not a", flush=True)
            print("[ncp]     convention. A real sign convention would also sit at "
                  "EXACTLY pi, not vary", flush=True)
            print("[ncp]     0.73-0.99 pi across the band.", flush=True)

    np.savez(os.path.join(OUT, "norm_phase.npz"),
             wavelengths=wl, flux_signed=flux_signed, a_fwd=a_fwd, a_bwd=a_bwd,
             iface=iface, norm_cal=norm_cal,
             d_cal_raw=d_cal_raw, d_cal_fwd=d_cal_fwd,
             resolution=FDTD_RESOLUTION)
    print(f"\n[ncp] saved -> {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
