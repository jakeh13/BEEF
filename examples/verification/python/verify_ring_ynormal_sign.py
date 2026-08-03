"""Is meep's backward eigenmode index sign-consistent with the forward one, on a
y-normal port?

WHY THIS EXISTS. ring_resonator_hybrid.py's ring-excited run normalizes with
`get_coeffs(..., forward=False)` -> `alpha[0, :, 1]` (its docstring: needed
"exactly once in this script", because that source launches -y), while every
other coefficient in the script is normalized and measured with `alpha[0, :, 0]`.
Measured 2026-08-03: after correcting for the ring port's own plane asymmetry,
`k1` and `k2` differ by a WAVELENGTH-FLAT 0.976*pi +- 0.003 across 100 nm, with
magnitudes equal to four decimals -- i.e. a convention constant, not a length.
`k2` and `t_ring` are exactly the two coefficients normalized through that one
`forward=False` call.

Applying that pi (with the plane fix) collapses the monolithic-vs-cascaded comb
disagreement from ~0.55 FSR to ~0.05 FSR and makes the dip counts match. But the
underlying claim -- that `alpha[..., 1]` carries -1 relative to `alpha[..., 0]`
-- is a statement about meep's convention and had never been tested. This script
tests it, because ring_resonator_validation.md already records one confidently
wrong "pi sign flip in T_ii" diagnosis, and a bare `-1` inserted on the strength
of a curve fit would be that mistake again.

Blocker 0 (the y-normal excite normalization, closed 2026-08-01 at fitted
h^2 order p=1.988) does NOT cover this: that test compared `max ||k1| - |k2||`,
magnitudes only. A phase pi is invisible to it.

THE TEST. Two mirror-image runs on the same uniform y-oriented strip:

    run "+y":  source at y = -S launching +y,  norm monitor at y = -S + M,
               coefficient read as alpha[..., 0]   -> a_plus
    run "-y":  source at y = +S launching -y,  norm monitor at y = +S - M,
               coefficient read as alpha[..., 1]   -> a_minus

The geometry is uniform and the cell is symmetric about y=0, so the two runs are
exact mirror images: each monitor sits the same distance M downstream of its own
source in its own propagation direction. Physics therefore demands
|a_minus| == |a_plus| and arg(a_minus) == arg(a_plus).

    a_minus / a_plus ~= +1  -> conventions agree; the ring's pi is NOT here, and
                              nothing in get_coeffs() should change.
    a_minus / a_plus ~= -1  -> the backward index carries the sign. Then EVERY
                              `forward=False` caller inherits it, so the fix
                              belongs in get_coeffs(), not in the ring script.

Deliberately NOT a ring: no bend, no coupler, no feedback. A straight guide is
the only geometry where the expected answer is known exactly, so any departure
is the convention and nothing else. 2-D and single-band, so it runs on a login
node in seconds.

Run:
    python /path/to/examples/verification/python/verify_ring_ynormal_sign.py
"""

import os

import numpy as np
import meep as mp

# Matched to ring_resonator_hybrid.py's strip and band.
WG_WIDTH = 0.45
WL_MIN, WL_MAX = 1.50, 1.60
N_FREQ = 5
RES = 30.0
DPML = 1.0
SIM_WIDTH = 4.0          # transverse (x) span, non-PML
S = 3.0                  # source offset from y=0, each run on its own side
M = 0.25                 # monitor offset downstream of the source (MONITOR_MARGIN)
PORT_APERTURE = 3.0
N_SI, N_SIO2 = 3.4777, 1.444      # frozen at band centre; a constant index is all
                                  # this test needs -- it asks about a sign, not dispersion

OUT = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "..",
    "example_output", "verification", "ring_ynormal_sign"))


def run(sign):
    """One run. sign=+1 launches +y from y=-S; sign=-1 launches -y from y=+S.

    Returns the norm-monitor coefficient, read through the index matching the
    launch direction: alpha[..., 0] for +y, alpha[..., 1] for -y -- exactly the
    pairing ring_resonator_hybrid.py's get_coeffs(forward=...) selects.
    """
    sy = 2 * (S + M + DPML) + 2.0
    sx = SIM_WIDTH + 2 * DPML
    fcen = 0.5 * (1.0 / WL_MIN + 1.0 / WL_MAX)
    df = 1.0 / WL_MIN - 1.0 / WL_MAX
    freqs = np.linspace(1.0 / WL_MAX, 1.0 / WL_MIN, N_FREQ)

    sim = mp.Simulation(
        cell_size=mp.Vector3(sx, sy),
        resolution=RES,
        boundary_layers=[mp.PML(DPML)],
        default_material=mp.Medium(index=N_SIO2),
        geometry=[mp.Block(                       # uniform strip along y, through both PMLs
            size=mp.Vector3(WG_WIDTH, mp.inf),
            center=mp.Vector3(0.0, 0.0),
            material=mp.Medium(index=N_SI),
        )],
        sources=[mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(0.0, -sign * S),
            size=mp.Vector3(PORT_APERTURE, 0),
            eig_band=1,
            eig_kpoint=mp.Vector3(0, float(sign), 0),
            eig_match_freq=True,
        )],
        force_complex_fields=True,
    )
    # M downstream of the source, in that run's own propagation direction.
    mon = sim.add_mode_monitor(freqs, mp.FluxRegion(
        center=mp.Vector3(0.0, -sign * S + sign * M),
        size=mp.Vector3(PORT_APERTURE, 0),
    ))
    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=1e-10))
    a = sim.get_eigenmode_coefficients(mon, [1], eig_parity=mp.NO_PARITY)
    idx = 0 if sign > 0 else 1       # forward index for +y, backward index for -y
    other = 1 - idx
    return (np.asarray(a.alpha[0, :, idx]),
            np.asarray(a.alpha[0, :, other]),
            np.asarray(freqs))


def main():
    print("[ysign] mirror-image y-normal launch test — is alpha[...,1] sign-consistent "
          "with alpha[...,0]?", flush=True)
    print(f"[ysign] res={RES:g} 2-D, S={S} um, monitor {M} um downstream, "
          f"{N_FREQ} freqs over {WL_MIN}-{WL_MAX} um", flush=True)

    print("[ysign] run +y (source at y=%+.2f, read alpha[...,0])" % (-S), flush=True)
    a_plus, off_plus, freqs = run(+1)
    print("[ysign] run -y (source at y=%+.2f, read alpha[...,1])" % (+S), flush=True)
    a_minus, off_minus, _ = run(-1)

    wl = 1.0 / freqs
    ratio = a_minus / a_plus

    print("\n[ysign] the wrong-index amplitudes should be ~0 (nothing travels the other way):")
    print(f"[ysign]   max |wrong index| / |right index|: "
          f"+y {np.max(np.abs(off_plus) / np.abs(a_plus)):.2e}, "
          f"-y {np.max(np.abs(off_minus) / np.abs(a_minus)):.2e}")

    print(f"\n[ysign] {'wl (um)':>9}{'|a_+y|':>11}{'|a_-y|':>11}"
          f"{'|ratio|':>10}{'arg(ratio)/pi':>15}")
    for i in range(len(wl)):
        print(f"[ysign] {wl[i]:>9.4f}{abs(a_plus[i]):>11.5f}{abs(a_minus[i]):>11.5f}"
              f"{abs(ratio[i]):>10.5f}{np.angle(ratio[i])/np.pi:>15.5f}")

    mag_spread = float(np.max(np.abs(np.abs(ratio) - 1.0)))
    arg_over_pi = np.angle(ratio) / np.pi
    arg_mean = float(np.mean(arg_over_pi))
    arg_spread = float(np.max(arg_over_pi) - np.min(arg_over_pi))

    print(f"\n[ysign] |ratio| deviation from 1 : {mag_spread:.3e}   "
          f"(mirror symmetry demands 0)")
    print(f"[ysign] arg(ratio)/pi mean       : {arg_mean:+.5f}")
    print(f"[ysign] arg(ratio)/pi spread     : {arg_spread:.3e}   "
          f"(flat => a convention constant, drifting => a path-length difference)")

    near_zero = abs(arg_mean) < 0.02
    near_pi = abs(abs(arg_mean) - 1.0) < 0.02
    flat = arg_spread < 0.02
    print()
    if not flat:
        print("[ysign] VERDICT: INCONCLUSIVE — arg(ratio) is not flat, so the two runs "
              "are not the\n         mirror images this test assumes. Check S/M/cell "
              "symmetry before reading anything\n         into the mean.")
    elif near_zero:
        print("[ysign] VERDICT: CONVENTIONS AGREE (ratio ~ +1). alpha[...,1] is NOT "
              "sign-flipped relative\n         to alpha[...,0]. The ring's flat pi is "
              "therefore NOT this call, and get_coeffs()\n         must not be changed. "
              "Look instead at the bend rotation between the ring's\n         "
              "normalization and measurement planes.")
    elif near_pi:
        print("[ysign] VERDICT: BACKWARD INDEX CARRIES -1 (ratio ~ -1). This reproduces "
              "the ring's flat pi\n         on a geometry with no bend and no coupler, so "
              "the sign is the convention and not\n         the device. EVERY "
              "forward=False caller inherits it -> fix belongs in get_coeffs(),\n         "
              "not in the ring script.")
    else:
        print(f"[ysign] VERDICT: UNEXPECTED — flat at {arg_mean:+.4f} pi, which is "
              "neither 0 nor +-1.\n         Do not patch anything on this; find out what "
              "sets that value first.")

    os.makedirs(OUT, exist_ok=True)
    np.savez(os.path.join(OUT, "ynormal_sign.npz"),
             wavelengths=wl, a_plus=a_plus, a_minus=a_minus,
             off_plus=off_plus, off_minus=off_minus,
             res=RES, S=S, M=M)
    print(f"\n[ysign] saved -> {OUT}/ynormal_sign.npz")

    orientation_test()


# ---------------------------------------------------------------------------
# TEST 2 — cross-orientation reciprocity on a 90 degree bend.
#
# Test 1 above compares y-normal against y-normal and finds the conventions
# agree. It therefore CANNOT explain the ring's pi, because the ring's pi is
# between k1 (bus-excited, x-normal source) and k2 (ring-excited, y-normal
# source) -- a comparison across port ORIENTATIONS that test 1 never makes.
#
# This test makes exactly that comparison, on the ring's own standalone bend
# geometry, which is the simplest structure that has one port of each kind:
#
#   run A:  excite the x stub (+x),  measure at the y stub  -> c_A
#   run B:  excite the y stub (-y),  measure at the x stub  -> c_B
#
# The bend is reciprocal and both runs use the same two planes at the same two
# distances, so reciprocity demands c_A == c_B as COMPLEX numbers. This is
# precisely the ring's k1/k2 pairing (k1 = ring<-bus, k2 = bus<-ring) reduced to
# a geometry with no coupler, no feedback and no second bend -- so a pi here is
# the port convention and nothing else.
#
# Note what makes this stronger than the ring's own reciprocity print: that one
# compares |k1| vs |k2| (magnitudes), which is exactly the check that let a
# phase pi through. This compares the full complex ratio.
# ---------------------------------------------------------------------------
BEND_R = 5.0             # matches ring_resonator_hybrid.py's R
BEND_STUB = 1.5          # matches STUB_LEN
BEND_MON = 0.25          # matches MONITOR_MARGIN


def _bend_polys():
    """The ring's standalone-bend outline: x stub -> quarter arc -> y stub.
    Mirrors ring_resonator_hybrid.py's build_bend_polys() shape."""
    th = np.linspace(-np.pi / 2, 0.0, 96)
    ri, ro = BEND_R - WG_WIDTH / 2, BEND_R + WG_WIDTH / 2
    outer = [(ro * np.cos(t), BEND_R + ro * np.sin(t)) for t in th]
    inner = [(ri * np.cos(t), BEND_R + ri * np.sin(t)) for t in th[::-1]]
    return [
        [(-BEND_STUB - DPML, -WG_WIDTH / 2), (0.0, -WG_WIDTH / 2),
         (0.0, WG_WIDTH / 2), (-BEND_STUB - DPML, WG_WIDTH / 2)],
        outer + inner,
        [(BEND_R - WG_WIDTH / 2, BEND_R), (BEND_R + WG_WIDTH / 2, BEND_R),
         (BEND_R + WG_WIDTH / 2, BEND_R + BEND_STUB + DPML),
         (BEND_R - WG_WIDTH / 2, BEND_R + BEND_STUB + DPML)],
    ]


def run_bend(excite):
    """excite='x': launch +x on the x stub, measure on the y stub.
       excite='y': launch -y on the y stub, measure on the x stub.
    Returns (measured/norm) per frequency — the bend's transfer coefficient,
    normalized at its own input plane, both planes BEND_MON inside the stub."""
    x0, x1 = -BEND_STUB - DPML, BEND_R + 1.0 + DPML
    y0, y1 = -1.0 - DPML, BEND_R + BEND_STUB + DPML
    sx, sy = x1 - x0, y1 - y0
    cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    fcen = 0.5 * (1.0 / WL_MIN + 1.0 / WL_MAX)
    df = 1.0 / WL_MIN - 1.0 / WL_MAX
    freqs = np.linspace(1.0 / WL_MAX, 1.0 / WL_MIN, N_FREQ)

    prisms = [mp.Prism([mp.Vector3(px, py) for (px, py) in poly],
                       height=mp.inf, material=mp.Medium(index=N_SI))
              for poly in _bend_polys()]
    x_plane = -BEND_STUB + BEND_MON                 # on the x stub
    y_plane = BEND_R + BEND_STUB - BEND_MON         # on the y stub

    if excite == "x":
        src = mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(-BEND_STUB, 0.0), size=mp.Vector3(0, PORT_APERTURE),
            eig_band=1, eig_kpoint=mp.Vector3(1, 0, 0), eig_match_freq=True)
        norm_reg = mp.FluxRegion(center=mp.Vector3(x_plane, 0.0),
                                 size=mp.Vector3(0, PORT_APERTURE))
        out_reg = mp.FluxRegion(center=mp.Vector3(BEND_R, y_plane),
                                size=mp.Vector3(PORT_APERTURE, 0))
        norm_fwd, out_fwd = True, True               # +x in, +y out
    else:
        src = mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(BEND_R, BEND_R + BEND_STUB), size=mp.Vector3(PORT_APERTURE, 0),
            eig_band=1, eig_kpoint=mp.Vector3(0, -1, 0), eig_match_freq=True)
        norm_reg = mp.FluxRegion(center=mp.Vector3(BEND_R, y_plane),
                                 size=mp.Vector3(PORT_APERTURE, 0))
        out_reg = mp.FluxRegion(center=mp.Vector3(x_plane, 0.0),
                                size=mp.Vector3(0, PORT_APERTURE))
        norm_fwd, out_fwd = False, False             # -y in, -x out

    sim = mp.Simulation(
        cell_size=mp.Vector3(sx, sy), resolution=RES,
        geometry_center=mp.Vector3(cx, cy),
        boundary_layers=[mp.PML(DPML)],
        default_material=mp.Medium(index=N_SIO2),
        geometry=prisms, sources=[src], force_complex_fields=True)
    nm = sim.add_mode_monitor(freqs, norm_reg)
    om = sim.add_mode_monitor(freqs, out_reg)
    sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=1e-9))
    g = lambda m, f: np.asarray(sim.get_eigenmode_coefficients(
        m, [1], eig_parity=mp.NO_PARITY).alpha[0, :, 0 if f else 1])
    return g(om, out_fwd) / g(nm, norm_fwd), np.asarray(freqs)


def orientation_test():
    print("\n[ysign] === TEST 2: cross-orientation reciprocity on a 90 deg bend ===",
          flush=True)
    print("[ysign] run A: excite x stub (+x), measure y stub", flush=True)
    c_x, freqs = run_bend("x")
    print("[ysign] run B: excite y stub (-y), measure x stub", flush=True)
    c_y, _ = run_bend("y")

    wl = 1.0 / freqs
    ratio = c_y / c_x
    print(f"\n[ysign] {'wl (um)':>9}{'|c_A| (x->y)':>14}{'|c_B| (y->x)':>14}"
          f"{'|ratio|':>10}{'arg(ratio)/pi':>15}")
    for i in range(len(wl)):
        print(f"[ysign] {wl[i]:>9.4f}{abs(c_x[i]):>14.5f}{abs(c_y[i]):>14.5f}"
              f"{abs(ratio[i]):>10.5f}{np.angle(ratio[i])/np.pi:>15.5f}")

    arg_over_pi = np.angle(ratio) / np.pi
    mean, spread = float(np.mean(arg_over_pi)), float(np.ptp(arg_over_pi))
    mag = float(np.max(np.abs(np.abs(ratio) - 1.0)))
    print(f"\n[ysign] |ratio| max deviation from 1 : {mag:.3e}   "
          f"(reciprocity demands 0)")
    print(f"[ysign] arg(ratio)/pi mean          : {mean:+.5f}")
    print(f"[ysign] arg(ratio)/pi spread        : {spread:.3e}")
    print()
    if mag > 0.05:
        print("[ysign] VERDICT: INCONCLUSIVE — magnitudes disagree, so the two runs are not "
              "seeing the\n         same structure/planes. Fix that before reading the phase.")
    elif spread > 0.05:
        print("[ysign] VERDICT: INCONCLUSIVE — arg(ratio) drifts, so this is a path-length "
              "difference\n         between the two runs, not a convention constant.")
    elif abs(mean) < 0.05:
        print("[ysign] VERDICT: ORIENTATIONS AGREE (ratio ~ +1). An x-normal port and a "
              "y-normal port\n         carry the SAME phase convention. Combined with test 1, "
              "the ring's flat pi is NOT\n         a meep port-convention artifact at all — "
              "look in the ring script's own algebra.")
    elif abs(abs(mean) - 1.0) < 0.05:
        print("[ysign] VERDICT: CROSS-ORIENTATION pi (ratio ~ -1). Reciprocity fails by a SIGN "
              "between an\n         x-normal and a y-normal port, on a reciprocal bend with "
              "matched planes. That is the\n         ring's k1/k2 pi, reproduced on a geometry "
              "with no coupler and no feedback.")
    else:
        print(f"[ysign] VERDICT: UNEXPECTED — flat at {mean:+.4f} pi. Do not patch anything; "
              "find out what sets it.")
    return c_x, c_y, wl


if __name__ == "__main__":
    main()
