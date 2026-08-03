"""Ring wiring items 2-4 -- the 2-port aggregate block, cross-checked.

ring_resonator_validation.md "Wiring to-do", decided 2026-08-01: build the
legacy-shaped 2-port [bus, ring] block FIRST (close_loop() Branch A, which is
implemented), use it purely as the numerical cross-check, and only then
implement Branch B for the physical 4-port form.

Why this ordering has diagnostic value: close_loop()'s ALGEBRA is already
proven against the legacy closed form symbolically
(tests/test_network_ring_crosscheck.py). So if the two routes below disagree,
the fault is necessarily in the new INPUTS -- item 2's hand-assembled
SMatrixBlock or item 3's loop propagator P -- not in the network layer. One
unknown at a time.

Two routes, identical inputs:

  A (legacy closed form)   E_thru = t + k1*k2*P / (1 - t_ring*P)
  B (beef.network)         SMatrixBlock[bus, ring] -> close_loop(ring, ring, P)

The coefficients come from RingResonatorCouplerSection through
FDTDSolver.solve_section() -- the library path -- rather than from the
reference script's hand-rolled meep. That substitution is the actual point of
item 2, and it is where an index/port-convention error would live.

Blocker 0 is CLOSED (2026-08-01, jobs 11616920/11616921): the y-normal excite
path's apparent reciprocity violation was second-order discretization, fitted
order p=1.988. So the ring-excited run below can be trusted.
"""
import os
import numpy as np

from beef.geometry.ring_coupler import RingResonatorCouplerSection
from beef.materials import SOI_220
from beef.network import close_loop, connect
from beef.network.block import SMatrixBlock
from beef.solvers.fdtd import FDTDSolver, PortSpec

# --- geometry / sweep, matching ring_resonator_hybrid.py's constants --------
WAVEGUIDE_WIDTH, GAP, R, L_C = 0.45, 0.25, 5.0, 4.0
L_SIDE, L_FAR = 2.0, L_C
S_MON, STUB_LEN, MONITOR_MARGIN = 0.5, 1.5, 0.25
PORT_APERTURE = 3.0
BEND_D_IN = BEND_D_OUT = STUB_LEN - MONITOR_MARGIN

FDTD_RESOLUTION = float(os.environ.get("RING_XC_RES", "20"))
THREE_D = os.environ.get("RING_XC_2D") != "1"
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, int(os.environ.get("RING_XC_NWL", "5"))

OUT = os.path.abspath(os.path.join(os.path.dirname(__file__),
      "..", "..", "..", "..", "example_output", "verification", "ring_network_crosscheck"))


N_IFACE = int(os.environ.get("RING_XC_NIFACE", "4"))   # supermode basis at the seam

# Straight-bus control span (2026-08-02). Job 11620125's straight-guide
# control repeated with THIS device's actual half-domain aperture instead of a
# centred one: same excite window, same measure window, over a stretch of bus
# with no bend and no coupling in it.
#
# DERIVED, never hardcoded -- two jobs died learning why. 11627585 (-5.0) and
# 11628207 (-8.0) both crashed in geometry construction because the cut that
# RingResonatorCouplerSection._check_cut() validates is NOT x1 + stub_right:
# solve_section() pads both stubs by its own dpml first (fdtd.py:1978), so the
# tested cut is x1 + stub_right + dpml. A cut strictly inside the near-bend arc
# has no cross-section that is a function of x; legal cuts land on an arc
# endpoint (stub_x or hc) or outside the arc entirely.
#
# The bound therefore depends on near_length, bend_radius, stub_right AND the
# solver's dpml -- and dpml is a SOLVER attribute, invisible from the geometry
# constants, which is exactly what the second attempt's arithmetic missed.
# Hardcoding a value that satisfies today's four re-breaks the moment any one
# of them moves.
def _straight_control_cut(g, stub_pad, dpml):
    """Largest legal control cut: lands the padded stub end exactly on stub_x.

    Returns the x to pass as the control section's inner bound, on whichever
    side `g` draws. Strict-inequality guard, so landing ON stub_x is legal.
    """
    s = -1.0 if g.side == "left" else 1.0
    return g.stub_x + s * (stub_pad + dpml)


# Env override for a deliberately shorter/longer span; None means derive.
_XS_ENV = os.environ.get("RING_XC_XSTRAIGHT")
X_STRAIGHT = float(_XS_ENV) if _XS_ENV else None

# NOTE on span length, still an open design call. With bus_clearance=2.0 the
# derived cut leaves outer_x -> cut = 0.5 um of straight bus. The control's
# premise is that a clipped launch sheds over a propagation DISTANCE, so 0.5 um
# may not answer it. The clean knob for more is the geometry's own
# bus_clearance (span = bus_clearance - stub_right - dpml, and the legality
# bound above is independent of it) -- NOT stub_right=0.0 for this one call,
# which buys length by changing port padding away from every other measurement
# in this script. See ring_resonator_validation.md item 2.


def _side_block(side, wavelengths, solver, n_iface):
    """One per-side section as a 3-port SMatrixBlock: [bus, ring, iface].

    Ports, and why each aperture is shaped the way it is:

      bus   -- x-normal at outer_x, the TRUE physical bus port. Half-domain
               aperture, split at the bus<->ring midpoint, so it sees the
               isolated bus mode and not a two-guide supermode.
      ring  -- y-normal at the stub top (ring_port()). The bus is at y=0 and
               this plane is at y=6.2, so it is naturally isolated.
      iface -- x-normal at iface_x, the interior SEAM against the mirrored
               section. This is NOT a port in the physical sense: the geometry
               deliberately draws a ring sliver parallel to the bus 0.25 um
               away from hc to iface_x, so the local eigenmodes here ARE
               supermodes. It therefore gets the FULL aperture and n_iface
               modes -- the opposite of the bus rule, and correct for the same
               reason. Reading this plane as a 1-mode isolated bus port is what
               lost ~70% of the power on 2026-08-01.

    Only forward transmission is filled; reflection and reverse transmission
    are zero, matching what solve_section() actually models (the same
    approximation HybridSolver's own cascade makes, and the case where
    connect() reduces to ordinary transfer-matrix multiplication).
    """
    g = RingResonatorCouplerSection(
        waveguide_width=WAVEGUIDE_WIDTH, gap=GAP, bend_radius=R,
        near_length=L_C, ring_stub_margin=S_MON, side=side,
    )
    rp = g.ring_port()
    x0, x1 = sorted((g.outer_x, g.iface_x))
    y0 = -2.0
    common = dict(stub_left=1.0, stub_right=1.0, y0=y0, y1=rp.y_boundary,
                  y_pml=solver.dpml)

    # half-domain bus aperture (excludes the ring centreline at y_ring)
    split = g.y_ring / 2.0
    lo = y0 + 0.25
    bus_c, bus_w = (lo + split) / 2.0, split - lo
    # Seam aperture: cover BOTH guides and little else.
    #
    # This was originally the full domain height (y=[-1.75, 5.95], 7.7 um).
    # That was wrong in a way the power audit exposed: at the seam plane the
    # only structure is the bus (y=0) and the ring sliver (y=0.70), together
    # spanning ~1 um; the remaining ~7 um is empty cladding. Such a window
    # supports exactly TWO guided supermodes (even/odd) -- every band past
    # that is radiation, and asking for them returned |T|^2 up to 9.2, i.e.
    # 9x energy GAIN, because a radiation band's normalization is meaningless.
    # Requesting more seam modes therefore makes the basis worse, not better.
    guide_lo = 0.0 - WAVEGUIDE_WIDTH / 2 - 0.5
    guide_hi = g.y_ring + WAVEGUIDE_WIDTH / 2 + 0.5
    span_lo, span_hi = max(guide_lo, y0 + 0.25), min(guide_hi, rp.y_boundary - 0.25)
    if_c, if_w = (span_lo + span_hi) / 2.0, span_hi - span_lo

    # For side="left" the bus port is at LOW x (outer_x < iface_x) and the seam
    # at HIGH x; for side="right" it is mirrored, so the roles of the two x
    # edges swap. sign picks the edge, and excite/measure read it oppositely.
    left = (side == "left")
    bus_excite = PortSpec(axis="x", sign=1 if left else -1, center=bus_c, width=bus_w)
    bus_measure = PortSpec(axis="x", sign=-1 if left else 1, center=bus_c, width=bus_w)
    if_excite = PortSpec(axis="x", sign=-1 if left else 1, center=if_c, width=if_w)
    if_measure = PortSpec(axis="x", sign=1 if left else -1, center=if_c, width=if_w)
    ring_excite = PortSpec(axis="y", sign=rp.sign, center=rp.x, width=rp.width)
    ring_measure = PortSpec(axis="y", sign=-rp.sign, center=rp.x, width=rp.width)

    print(f"[xc] side={side}: bus aperture y=[{lo:.3f},{split:.3f}] "
          f"(ring centreline {g.y_ring:.3f} excluded); seam aperture "
          f"y=[{span_lo:.3f},{span_hi:.3f}] with {n_iface} modes", flush=True)

    # Resolve and validate the launch control's cut UP FRONT, before any solve.
    # Cut legality is a sub-second geometry check, but jobs 11627585/11628207
    # each burned ~3:20 on 24 cores discovering it only after the two expensive
    # bus solves had already run. Fail in the first second instead.
    x_straight = None
    if left:
        x_straight = (X_STRAIGHT if X_STRAIGHT is not None
                      else _straight_control_cut(g, common["stub_right"], solver.dpml))
        s = -1.0 if g.side == "left" else 1.0
        g._check_cut(x_straight - s * (common["stub_right"] + solver.dpml))
        print(f"[xc] side={side}: launch-control span [{x0:.2f},{x_straight:.2f}] "
              f"= {abs(x_straight - x0):.2f} um of straight bus "
              f"(cut validated against stub_x={g.stub_x:.2f})", flush=True)

    if left:
        # left: physical ports feed the seam
        T_if_bus = solver.solve_section(g, wavelengths, x0, x1,
            n_modes_in=1, n_modes_out=n_iface,
            excite_port=bus_excite, measure_port=if_measure, **common)
        T_if_ring = solver.solve_section(g, wavelengths, x0, x1,
            n_modes_in=1, n_modes_out=n_iface,
            excite_port=ring_excite, measure_port=if_measure, **common)
        # THE THIRD EXIT, added 2026-08-02. This section is a genuine 3-port
        # (bus, ring, iface), but the audit below only ever measured the two
        # paths INTO the seam. Power coupled from the bus up the ring stub
        # leaves through ring_measure, and was being counted as loss.
        #
        # That matters because every hypothesis formed on 2026-08-01 -- seven
        # of them, all wrong -- rested on the implicit premise that
        # capture(iface<-bus) alone had to approach 1. It never did. The
        # quantity that must approach 1 is the SUM printed below.
        #
        # Diagnostic only: deliberately NOT folded into `data`. Adding a
        # forward bus->ring term would change the very composition this
        # script exists to cross-check, and one unknown at a time is the
        # whole point of doing the 2-port aggregate form first.
        T_ring_bus = solver.solve_section(g, wavelengths, x0, x1,
            n_modes_in=1, n_modes_out=1,
            excite_port=bus_excite, measure_port=ring_measure, **common)
        # Per-section power audit. sum_n |T[n]|^2 over the seam's own modes is
        # the power actually captured by the n_iface-mode basis. For a lossless
        # section that must be ~1; anything well below it means the section
        # itself is losing power BEFORE any composition, which composing can
        # only multiply. This is the diagnostic that separates "the section is
        # lossy" from "the extraction/composition is wrong".
        for nm, T in (("iface<-bus", T_if_bus), ("iface<-ring", T_if_ring),
                      ("ring<-bus", T_ring_bus)):
            cap = (np.abs(T[:, :, 0]) ** 2).sum(axis=1)
            print(f"[xc]   {side} {nm}: captured power per wl = "
                  + " ".join(f"{c:.4f}" for c in cap), flush=True)
        # The check this run exists for: both exits from a bus excitation.
        cap_if = (np.abs(T_if_bus[:, :, 0]) ** 2).sum(axis=1)
        cap_ring = (np.abs(T_ring_bus[:, :, 0]) ** 2).sum(axis=1)
        total = cap_if + cap_ring
        print(f"[xc]   {side} TOTAL out of bus (iface + ring) = "
              + " ".join(f"{c:.4f}" for c in total), flush=True)
        print(f"[xc]   {side} deficit vs 1.0 = "
              + " ".join(f"{1.0 - c:+.4f}" for c in total), flush=True)

        # --- Where the bus power actually goes. Added 2026-08-02 after
        # ring<-bus came back 0.0000 and killed the third-exit hypothesis.
        # Both diagnostic; neither is folded into `data`.
        #
        # These two only discriminate TOGETHER. Forward-alone cannot separate
        # "reflected" from "radiated", and those two answers point at
        # different fixes -- the same trap as measuring iface<-bus without
        # ring<-bus.
        #
        # A -- REFLECTION at the bus port. Note bus_measure carries sign=-1,
        # and for a MEASURE port that is the LOW edge reading the backward
        # band (fdtd.py's excite/measure sign conventions are mirrored), so
        # on this left section it monitors light coming back out of outer_x.
        # _side_block fills forward transmission ONLY, so until now a
        # reflecting section and a radiating one were indistinguishable here.
        R_bus = solver.solve_section(g, wavelengths, x0, x1,
            n_modes_in=1, n_modes_out=1,
            excite_port=bus_excite, measure_port=bus_measure, **common)
        # B -- FORWARD power at the seam plane, in the bus's OWN basis. Same
        # window as bus_excite, but sign=+1 on a measure port puts it at the
        # HIGH edge (iface_x) reading forward. This asks whether the power
        # reaches the seam at all, independently of whether the seam's
        # n_iface supermodes can represent it once it arrives.
        bus_at_seam = PortSpec(axis="x", sign=1, center=bus_c, width=bus_w)
        T_bus_bus = solver.solve_section(g, wavelengths, x0, x1,
            n_modes_in=1, n_modes_out=1,
            excite_port=bus_excite, measure_port=bus_at_seam, **common)
        for nm, T in (("T bus<-bus (seam plane, bus basis)", T_bus_bus),
                      ("R bus<-bus (reflection at outer_x)", R_bus)):
            cap = (np.abs(T[:, :, 0]) ** 2).sum(axis=1)
            print(f"[xc]   {side} {nm}: power per wl = "
                  + " ".join(f"{c:.4f}" for c in cap), flush=True)
        budget = ((np.abs(T_bus_bus[:, :, 0]) ** 2).sum(axis=1)
                  + (np.abs(R_bus[:, :, 0]) ** 2).sum(axis=1)
                  + cap_ring)
        print(f"[xc]   {side} BUDGET fwd(bus basis) + refl + ring = "
              + " ".join(f"{c:.4f}" for c in budget), flush=True)
        print(f"[xc]   {side} unaccounted (-> radiation) = "
              + " ".join(f"{1.0 - c:+.4f}" for c in budget), flush=True)

        # C -- THE LAUNCH CONTROL. Added 2026-08-02 after A and B came back
        # 0.001-0.02 (reflection) and 0.05-0.33 (forward), leaving ~80%
        # unaccounted. Two readings survive that and they need different
        # fixes: the bend genuinely radiates, OR the source launches a
        # clipped non-eigenmode field that sheds into radiation immediately.
        # The bus half-domain window is y=[-1.750, 0.350] -- 1.75 um below a
        # centreline at y=0 but only 0.35 um above it, to clear the ring at
        # y=0.70. A 0.45 um strip mode's evanescent tail runs past 0.35 um.
        #
        # Same excite window, same measure window, same solver settings, over
        # a span with NO bend and no coupling. Job 11620125 could NOT settle
        # this: its offset-aperture cases used the SAME aperture at both ends,
        # where a clipped mode still transmits perfectly into itself. Here the
        # launch is tested against propagation instead.
        T_straight = solver.solve_section(g, wavelengths, x0, x_straight,
            n_modes_in=1, n_modes_out=1,
            excite_port=bus_excite, measure_port=bus_at_seam, **common)
        cap_st = (np.abs(T_straight[:, :, 0]) ** 2).sum(axis=1)
        print(f"[xc]   {side} CONTROL bus<-bus straight [{x0:.2f},{x_straight:.2f}]"
              f" ({abs(x_straight - x0):.2f} um, no bend): power per wl = "
              + " ".join(f"{c:.4f}" for c in cap_st), flush=True)
        print(f"[xc]   {side}   ~1 -> launch sound, loss is AT THE BEND; "
              f"~0.2 -> the half-domain aperture clips the bus mode and every "
              f"bus-excited number since defect 3 inherits it", flush=True)
        M = len(wavelengths); N = 2 + n_iface
        data = np.zeros((M, N, N), dtype=complex)
        data[:, 2:, 0:1] = T_if_bus      # iface <- bus
        data[:, 2:, 1:2] = T_if_ring     # iface <- ring
        return SMatrixBlock(wavelengths=np.asarray(wavelengths),
                            ports=("bus_L", "ring_L", "iface"),
                            modes_per_port={"bus_L": 1, "ring_L": 1, "iface": n_iface},
                            data=data)
    # right: the seam feeds the physical ports
    T_bus_if = solver.solve_section(g, wavelengths, x0, x1,
        n_modes_in=n_iface, n_modes_out=1,
        excite_port=if_excite, measure_port=bus_measure, **common)
    T_ring_if = solver.solve_section(g, wavelengths, x0, x1,
        n_modes_in=n_iface, n_modes_out=1,
        excite_port=if_excite, measure_port=ring_measure, **common)
    for nm, T in (("bus<-iface", T_bus_if), ("ring<-iface", T_ring_if)):
        cap = (np.abs(T[:, 0, :]) ** 2).sum(axis=1)
        print(f"[xc]   {side} {nm}: summed |T|^2 over seam modes = "
              + " ".join(f"{c:.4f}" for c in cap), flush=True)
    M = len(wavelengths); N = n_iface + 2
    data = np.zeros((M, N, N), dtype=complex)
    data[:, n_iface:n_iface+1, :n_iface] = T_bus_if     # bus  <- iface
    data[:, n_iface+1:n_iface+2, :n_iface] = T_ring_if  # ring <- iface
    return SMatrixBlock(wavelengths=np.asarray(wavelengths),
                        ports=("iface", "bus_R", "ring_R"),
                        modes_per_port={"iface": n_iface, "bus_R": 1, "ring_R": 1},
                        data=data)


def coupler_coefficients(wavelengths):
    """t, k1, k2, t_ring for the FULL coupler, by composing the two per-side
    sections across their shared seam.

    This is the correction made 2026-08-01: `iface_x` is a SEAM, not a port,
    and the legacy `t` is bus_left -> bus_right, which spans BOTH sides. It
    cannot be read off a single per-side section at all. connect() joins them;
    with zero self-terms it reduces to ordinary transfer-matrix
    multiplication through the n_iface-mode interface basis.
    """
    solver = FDTDSolver(
        res=FDTD_RESOLUTION, sim_width=PORT_APERTURE, sim_height=1.5, dpml=0.5,
        source_margin=0.0, monitor_margin=0.75,
        three_d=THREE_D, material_stack=SOI_220, decay_tol=1e-9,
    )
    print("[xc] left section  (bus_L, ring_L -> seam)", flush=True)
    left = _side_block("left", wavelengths, solver, N_IFACE)
    print("[xc] right section (seam -> bus_R, ring_R)", flush=True)
    right = _side_block("right", wavelengths, solver, N_IFACE)

    print("[xc] connect() across the seam", flush=True)
    full = connect(left, "iface", right, "iface")
    print(f"[xc] composed ports: {full.ports}", flush=True)

    def S(out_port, in_port):
        r = full.port_slice(out_port); c = full.port_slice(in_port)
        return full.data[:, r, c][:, 0, 0]

    return (S("bus_R", "bus_L"), S("ring_R", "bus_L"),
            S("ring_R", "ring_L"), S("bus_R", "ring_L"))


def loop_propagator(wavelengths, neff, b_raw):
    """Wiring item 3. The reference-plane correction is LOAD-BEARING for
    resonance position -- reproduce it exactly or every resonance shifts."""
    phi_bend = 2 * np.pi * neff * (BEND_D_IN + BEND_D_OUT) / wavelengths
    b_corr = b_raw * np.exp(-1j * phi_bend)
    L_outer = 2 * (L_SIDE - S_MON) + L_FAR
    phi = 2 * np.pi * neff * L_outer / wavelengths
    return (b_corr ** 2) * np.exp(1j * phi), b_corr, L_outer


def build_block(wavelengths, t, k1, k2, t_ring):
    """Wiring item 2. data[m][i, j] is S[out=i, in=j] over ports [bus, ring],
    matching close_loop()'s documented substitution: S_bus,bus=t,
    S_bus,ring=k2, S_ring,bus=k1, S_ring,ring=t_ring."""
    M = len(wavelengths)
    data = np.zeros((M, 2, 2), dtype=complex)
    data[:, 0, 0] = t
    data[:, 0, 1] = k2
    data[:, 1, 0] = k1
    data[:, 1, 1] = t_ring
    return SMatrixBlock(ports=["bus", "ring"], modes_per_port={"bus": 1, "ring": 1},
                        wavelengths=np.asarray(wavelengths), data=data)


def main():
    wl = np.linspace(WL_MIN, WL_MAX, N_WL)
    print(f"[xc] res={FDTD_RESOLUTION} three_d={THREE_D} n_wl={N_WL}", flush=True)

    t, k1, t_ring, k2 = coupler_coefficients(wl)
    # Bend coefficient: reuse the coupler's own ring-side self-term magnitude
    # as a stand-in ONLY if no bend run is supplied. A real assembly uses
    # run_bend_fdtd()'s b; this script's job is the ROUTE comparison, and both
    # routes consume whatever b is handed in, so the cross-check is valid
    # either way.
    b_raw = np.full(len(wl), 0.99 + 0.0j)
    neff = np.full(len(wl), 2.3395)   # gap 3 converged value at res 100

    P, b_corr, L_outer = loop_propagator(wl, neff, b_raw)

    # Route A -- legacy closed form
    E_A = t + k1 * k2 * P / (1.0 - t_ring * P)

    # Route B -- beef.network
    block = build_block(wl, t, k1, k2, t_ring)
    closed = close_loop(block, port_out="ring", port_in="ring", loop=P)
    E_B = closed.data[:, 0, 0]

    print(f"\n{'wl':>7} {'|E| legacy':>12} {'|E| network':>12} {'abs diff':>12}")
    worst = 0.0
    for i, w in enumerate(wl):
        d = abs(E_A[i] - E_B[i]); worst = max(worst, d)
        print(f"{w:7.3f} {abs(E_A[i]):12.8f} {abs(E_B[i]):12.8f} {d:12.3e}")

    # --- PHYSICS checks on the coefficients themselves --------------------
    # The route comparison above is nearly tautological: both consume the same
    # inputs, and close_loop()'s algebra is already proven, so it agrees to
    # machine precision even when the coefficients are nonsense. It did
    # exactly that on 2026-08-01 while |k1|=0.03 against |k2|=0.38. These
    # checks are what actually validate item 2's extraction.
    print(f"\n{'wl':>7} {'|t|':>9} {'|k1|':>9} {'|k2|':>9} {'|t_ring|':>9}"
          f" {'|t|^2+|k1|^2':>13}")
    for i, w in enumerate(wl):
        print(f"{w:7.3f} {abs(t[i]):9.5f} {abs(k1[i]):9.5f} {abs(k2[i]):9.5f} "
              f"{abs(t_ring[i]):9.5f} {abs(t[i])**2 + abs(k1[i])**2:13.5f}")

    unitarity = np.abs(t)**2 + np.abs(k1)**2
    recip = np.abs(np.abs(k1) - np.abs(k2)).max()
    print(f"\n[xc] bus-side power sum |t|^2+|k1|^2: min {unitarity.min():.4f} "
          f"max {unitarity.max():.4f}   (lossless coupler -> ~1)")
    print(f"[xc] reciprocity max ||k1|-|k2||: {recip:.4f}   (-> 0)")
    bad = []
    if unitarity.min() < 0.80:
        bad.append(f"power sum {unitarity.min():.4f} far below 1 -- the excited "
                   "power is leaving through a port this 2-port model is not "
                   "measuring (check each PortSpec's edge/aperture)")
    if recip > 0.10:
        bad.append(f"||k1|-|k2|| = {recip:.4f} -- a reciprocal passive coupler "
                   "must give |k1| == |k2|; a large gap means one of the two "
                   "runs is reading the wrong plane")
    if bad:
        print("\n[xc] PHYSICS CHECK FAILED:")
        for b in bad:
            print(f"[xc]   - {b}")
    else:
        print("[xc] physics checks PASS (power sum ~1, reciprocity holds)")

    print(f"\n[xc] worst |E_legacy - E_network| = {worst:.3e}")
    print(f"[xc] L_outer = {L_outer:.4f} um")
    if worst < 1e-12:
        print("[xc] VERDICT: routes AGREE to machine precision -- item 2's block "
              "assembly and item 3's P are consistent with the proven closed form.")
    else:
        print("[xc] VERDICT: routes DISAGREE. close_loop()'s algebra is already "
              "proven symbolically, so the fault is in the INPUTS: the block "
              "index convention (item 2) or P (item 3).")

    os.makedirs(OUT, exist_ok=True)
    np.savez(os.path.join(OUT, "crosscheck.npz"),
             wavelengths=wl, t=t, k1=k1, k2=k2, t_ring=t_ring, P=P,
             E_legacy=E_A, E_network=E_B, worst_abs_diff=worst,
             res=FDTD_RESOLUTION, three_d=THREE_D, L_outer=L_outer)
    print(f"[xc] saved -> {OUT}/crosscheck.npz")


if __name__ == "__main__":
    main()
