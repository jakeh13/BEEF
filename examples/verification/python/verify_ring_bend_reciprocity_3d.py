"""Ring blocker 0 -- y-normal excite reciprocity, in 3-D.

ring_resonator_validation.md's boxed blocker recorded |T| x->y vs y->x
differing by ~0.73 and localized it to the flux-calibrated normalization.
The 2026-08-01 investigation refuted that (see the doc): the conventions
agree once the norm monitor leaves the source near field, and the "validated"
x->y direction was itself 27% low because the y-normal port's aperture was
buried in PML.

Every one of those runs was 2-D (`three_d=False`). Production ring runs 3-D,
and in 2-D meep the TE/TM polarizations decouple -- an x-propagating and a
y-propagating mode do not carry the same field components, which reciprocity
would report as exactly this kind of constant ratio. So the open question is
whether the residual ~0.75 is a 2-D artifact that never applied to the ring.

Both fixes from that investigation are applied here:
  stub_right      -- x headroom so the y-normal port is not inside PML
  monitor_margin  -- norm plane out of the y-source's near field

Pass  => reciprocity holds in 3-D; blocker 0 is a 2-D artifact and the wiring
         is unblocked.
Fail  => a real orientation asymmetry survives into 3-D; chase it there.
"""
import os
import numpy as np

from beef.geometry.ring_bend import RingResonatorBendSection
from beef.materials import SOI_220
from beef.solvers.fdtd import FDTDSolver, PortSpec

RES = float(os.environ.get("RING_RECIP_RES", "15"))
THREE_D = os.environ.get("RING_RECIP_2D") != "1"
STUB_RIGHT = 2.0        # x headroom -- keeps the y port clear of PML
MONITOR_MARGIN = 0.75   # out of the y-source near field
APERTURE = 3.0
WAVELENGTHS = np.array([1.53, 1.55, 1.57])
OUT = os.path.abspath(os.path.join(os.path.dirname(__file__),
      "..", "..", "..", "..", "example_output", "verification", "ring_bend_reciprocity"))


def main():
    g = RingResonatorBendSection(waveguide_width=0.45, bend_radius=5.0, output_margin=0.5)
    s = FDTDSolver(res=RES, sim_width=3.0, sim_height=1.5, dpml=0.5,
                   source_margin=0.0, monitor_margin=MONITOR_MARGIN,
                   three_d=THREE_D, material_stack=SOI_220, decay_tol=1e-9)
    port = g.out_port()
    # y0=-2.0, NOT -1.0: at -1.0 the X-normal port's own aperture overruns the
    # PML in y (spans y=[-3.75,-0.75] against an inner edge at +/-3.25). In the
    # y->x direction that port is the MEASURE side, so the clip removes
    # collected power and reads as an orientation asymmetry. solve_section()'s
    # aperture guard now raises on it; the first 3-D run (11616611) predates
    # the guard and was contaminated exactly this way.
    common = dict(n_modes_in=1, n_modes_out=1, stub_left=1.0, stub_right=STUB_RIGHT,
                  y0=-2.0, y1=port.y_boundary, y_pml=s.dpml)

    print(f"[recip] res={RES} three_d={THREE_D} stub_right={STUB_RIGHT} "
          f"monitor_margin={MONITOR_MARGIN}", flush=True)

    T_xy = s.solve_section(g, WAVELENGTHS, 0.0, g.bend_radius,
        excite_port=PortSpec(axis="x", sign=1, center=0.0, width=APERTURE),
        measure_port=PortSpec(axis="y", sign=1, center=g.bend_radius, width=APERTURE),
        **common)
    print("[recip] x->y done", flush=True)

    T_yx = s.solve_section(g, WAVELENGTHS, 0.0, g.bend_radius,
        excite_port=PortSpec(axis="y", sign=-1, center=g.bend_radius, width=APERTURE),
        measure_port=PortSpec(axis="x", sign=-1, center=0.0, width=APERTURE),
        **common)
    print("[recip] y->x done", flush=True)

    t_xy, t_yx = np.abs(T_xy[:, 0, 0]), np.abs(T_yx[:, 0, 0])
    ratio = t_yx / t_xy
    print(f"\n{'wl':>6} {'|T| x->y':>10} {'|T| y->x':>10} {'ratio':>9}")
    for i, w in enumerate(WAVELENGTHS):
        print(f"{w:6.2f} {t_xy[i]:10.6f} {t_yx[i]:10.6f} {ratio[i]:9.6f}")

    print(f"\n[recip] mean ratio {ratio.mean():.6f}  (reciprocity requires 1.0)")
    print("[recip] x->y should also be ~0.98 -- a 5 um-radius 90-deg bend in "
          "450x220 SOI is near-lossless. A low value there means the port "
          "aperture is clipped again, not that the bend is lossy.")
    if np.allclose(ratio, 1.0, rtol=0.15):
        print("[recip] VERDICT: reciprocity HOLDS -> the 2-D ~0.75 was an "
              "artifact; blocker 0 does not apply to the production 3-D path.")
    else:
        print("[recip] VERDICT: reciprocity VIOLATED in 3-D -> real asymmetry, "
              "blocker 0 stands. Do not build the wiring on this path.")

    os.makedirs(OUT, exist_ok=True)
    np.savez(os.path.join(OUT, "reciprocity.npz"),
             wavelengths=WAVELENGTHS, T_xy=T_xy, T_yx=T_yx, ratio=ratio,
             res=RES, three_d=THREE_D, stub_right=STUB_RIGHT,
             monitor_margin=MONITOR_MARGIN)
    print(f"[recip] saved -> {OUT}/reciprocity.npz")


if __name__ == "__main__":
    main()
