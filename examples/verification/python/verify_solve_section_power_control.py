"""CONTROL: does FDTDSolver.solve_section() conserve power on a trivial guide?

Written 2026-08-01 after six successive wrong hypotheses about why the ring
coupler's extracted coefficients lose ~70-99% of their power
(ring_resonator_validation.md, wiring item 2). Every one of those was inferred
from the composed result on a complicated geometry. This is the check that
should have come first.

A uniform straight waveguide, both ports x-normal, one mode, no ring, no
bend, no seam, no coupling. |T|^2 MUST be ~1 -- there is nowhere for the
power to go. Run at the SAME solver settings the ring cross-check uses, so a
low answer here indicts the configuration (resolution, decay_tol, aperture
convention, 2-D-ness) rather than the ring geometry, and a ~1 answer clears
the solver and localizes the fault to the ring's own ports/geometry.

Sweeps the aperture width too, because that is the one knob the ring work
kept changing without ever checking what it does on a case with a known
answer.
"""
import os
import numpy as np

from beef.geometry.taper import Taper
from beef.materials import SOI_220
from beef.solvers.fdtd import FDTDSolver, PortSpec

RES = float(os.environ.get("CTL_RES", "20"))
THREE_D = os.environ.get("CTL_2D") != "1"
WL = np.array([1.50, 1.55, 1.60])
WIDTH = 0.45
LENGTH = 8.0
SIM_WIDTH = 3.0


def run(aperture, use_portspec, y0=None, y1=None, center=0.0):
    g = Taper(length=LENGTH, width_start=WIDTH, width_end=WIDTH,
              waveguide_height=0.22)
    s = FDTDSolver(res=RES, sim_width=SIM_WIDTH, sim_height=1.5, dpml=0.5,
                   source_margin=0.0, monitor_margin=0.75, three_d=THREE_D,
                   material_stack=SOI_220, decay_tol=1e-9)
    kw = dict(n_modes_in=1, n_modes_out=1, stub_left=1.0, stub_right=1.0)
    if use_portspec:
        kw["excite_port"] = PortSpec(axis="x", sign=1, center=center, width=aperture)
        kw["measure_port"] = PortSpec(axis="x", sign=1, center=center, width=aperture)
    if y0 is not None:
        kw.update(y0=y0, y1=y1, y_pml=s.dpml)
    T = s.solve_section(g, WL, -LENGTH / 2, LENGTH / 2, **kw)
    return np.abs(T[:, 0, 0]) ** 2


def main():
    print(f"[ctl] res={RES} three_d={THREE_D}  uniform straight "
          f"{WIDTH} um x {LENGTH} um -- |T|^2 must be ~1\n", flush=True)
    rows = []

    print("[ctl] A: layout path (no PortSpec at all) -- the baseline every "
          "other device uses", flush=True)
    rows.append(("layout path, default aperture", run(None, False)))

    for ap in (1.0, 2.0, 3.0):
        print(f"[ctl] B: PortSpec, aperture {ap} um, default y-domain", flush=True)
        rows.append((f"PortSpec aperture {ap}", run(ap, True)))

    # tall y-domain, as the ring sections use (y0/y1 forced, guide near the floor)
    print("[ctl] C: PortSpec + TALL forced y-domain (as the ring uses)", flush=True)
    rows.append(("PortSpec ap 2.1, y=[-2.0, 6.2]", run(2.1, True, y0=-2.0, y1=6.2)))

    # D -- the one variable the original control never swept: an aperture
    # OFFSET from the guide it is meant to measure. Every case above centres
    # on the guide; the ring's half-domain bus aperture is centred 0.7 um off
    # it (window y=[-1.75,0.35] against a guide at y=0). If the offset is what
    # costs the power, this is where it shows, on a geometry whose answer is
    # known to be 1.
    # Run these in the RING's own y-domain (y=[-2.0, 6.2]); the control's
    # default 4.0 um cell is too short for an offset window and the aperture
    # guard rightly refuses it. Case C above already showed that this tall
    # domain costs nothing when the aperture is centred, so any drop here is
    # attributable to the OFFSET alone.
    for c, w in ((0.0, 2.1), (-0.35, 2.1), (-0.7, 2.1)):
        tag = "centred" if c == 0.0 else "OFFSET"
        print(f"[ctl] D: {tag} aperture centre {c} width {w} (tall domain)", flush=True)
        rows.append((f"{tag} centre {c}, width {w}, tall",
                     run(w, True, y0=-2.0, y1=6.2, center=c)))

    print(f"\n{'case':>36} " + " ".join(f"{w:>9.2f}" for w in WL))
    for label, v in rows:
        print(f"{label:>36} " + " ".join(f"{x:9.5f}" for x in v))

    print("\n[ctl] READING:")
    print("[ctl]   A ~1  -> solver+config are sound; the ring's ports/geometry are at fault")
    print("[ctl]   A low -> the configuration itself loses power; everything")
    print("[ctl]            built on it (ring item 2) is downstream of THIS")
    print("[ctl]   B vs A -> whether PortSpec itself costs power vs the layout path")
    print("[ctl]   C vs B -> whether a tall forced y-domain is what costs it")


if __name__ == "__main__":
    main()
