"""Ring gap 3 -- bare-strip neff mesh convergence.

ring_resonator_validation.md gap 3: the ring's bare single-mode strip neff is
fed into THREE straight-propagation phases and the reference-plane correction
on b(lambda), so an neff error lands directly in resonance position. It has
never been measured at any resolution, for either engine.

This sweeps RESOLUTION through solve_neff()'s exact recipe (same meow
Mesh2D/CrossSection/compute_modes_meep path, same MATERIAL_STACK, same
MESH_WIDTH/HEIGHT) and reports both neff and the quantity that actually
matters -- the round-trip phase phi = 2*pi*neff*L_outer/lambda, in radians,
since that is what sets resonance position.

Mode solve only: no FDTD, no time stepping.
"""
import os
import numpy as np
import meow as mw
from meow import compute_modes_meep

from beef.materials import SOI_220
from beef import extrude_polygons

WG_WIDTH = 0.45
MESH_WIDTH, MESH_HEIGHT = 4.0, 2.0
MATERIAL_STACK = SOI_220
# L_outer is the analytic-phase path length, used only to convert a d(neff) into
# a d(phase) -- but its value sets the lever arm linearly, so it has to be the
# real one.
#
# FIXED 2026-08-03. This was `2 * 2.0 + 10.0` = 14.0 um, on a comment asserting
# L_C=10.0 was "the shipped scale". It is not: ring_resonator_hybrid.py:125 ships
# L_C = 4.0, and its round-trip phase (`:805`) uses a reference-plane pullback the
# old expression dropped entirely. Job 11634859 printed `L_outer = 7.0000 um`
# directly, confirming the value below. The old 14.0 was wrong twice over -- wrong
# L_C, and no -2*S_MON term -- so every Delta-phi this script reported was
# exactly 2x too large, and the headline "shipped res 40 is a QUARTER-FSR off
# res 100" is really an EIGHTH (12.6%).
#
# Mirrored from ring_resonator_hybrid.py rather than re-derived; keep in sync.
L_SIDE, S_MON, L_FAR = 2.0, 0.5, 4.0          # == ring_resonator_hybrid.py:152,207,153
L_OUTER = 2 * (L_SIDE - S_MON) + L_FAR        # == ring_resonator_hybrid.py:349  -> 7.0

# THE LEVER ARM IS NOT THE SAME FOR THE TWO ROUTES, as of 2026-08-03.
#
# What this script measures is a mode-solver mesh error in neff. That error only
# reaches the resonance position through the straight length the ANALYTIC phase
# covers. Whatever straight length is instead carried inside an FDTD-measured
# bend coefficient is immune to it -- it has its own, separate discretization
# error, which is what the res_sweep study measures.
#
# The two routes now divide that length differently:
#
#   monolithic (ring_resonator_hybrid.py:1148, phase-strip form)
#       b_corrected strips BEND_D_IN+BEND_D_OUT back to the bare arc, so the
#       analytic phase spans the full L_OUTER = 7.00 um.
#
#   cascaded (ring_resonator_hybrid_cascaded.py:1657, DC-style form, default
#   since 2026-08-03)
#       b_d KEEPS its 1.25 um on each side of each arc, so the analytic phase
#       spans only L_OUTER - BEND_STRAIGHT_TOTAL = 7.00 - 5.00 = 2.00 um.
#
# The two P assemblies are algebraically identical (max |P1-P2| = 7.2e-15), so
# this is NOT a physics difference -- but it IS a 3.5x difference in how much
# an neff mesh error matters, so a single dphi column would misreport one route.
#
# THIS IS THE THIRD TIME THIS LEVER ARM HAS BEEN WRONG. First it was 14.0 (wrong
# L_C, and missing the -2*S_MON pullback -> 2x too large). Fixed to 7.0. Then the
# bend accounting moved 5.00 um of straight out of the analytic phase and into b,
# which left 7.0 3.5x too large for the cascaded route without anything in this
# file changing. Report both, and derive them from the constants above rather
# than hardcoding, so the next accounting change shows up as a number instead of
# a stale comment.
BEND_D_IN = BEND_D_OUT = 1.5 - 0.25           # STUB_LEN - MONITOR_MARGIN -> 1.25
BEND_STRAIGHT_TOTAL = 2 * (BEND_D_IN + BEND_D_OUT)            # 2 far bends -> 5.0
LEVER_ARMS = {                                # route label -> analytic straight length
    "monolithic (phase-strip, L_OUTER)": L_OUTER,
    "cascaded (DC-style, L_OUTER-BEND_STRAIGHT_TOTAL)": L_OUTER - BEND_STRAIGHT_TOTAL,
}
RESOLUTIONS = [20.0, 30.0, 40.0, 60.0, 80.0, 100.0]
WAVELENGTHS = [1.50, 1.55, 1.60]
OUT = os.path.abspath(os.path.join(os.path.dirname(__file__),
      "..", "..", "..", "..", "example_output", "verification", "ring_gap3_neff"))


def rect_poly(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def solve_neff(wl, res):
    length = 1.0
    gds = {(2, 0): [rect_poly(0.0, -WG_WIDTH / 2, length, WG_WIDTH / 2)],
           (1, 0): [rect_poly(0.0, -MESH_WIDTH / 2, length, MESH_WIDTH / 2)]}
    structs = extrude_polygons(gds, MATERIAL_STACK.extrusion_rules)
    mesh = mw.Mesh2D(
        x=np.linspace(-MESH_WIDTH / 2, MESH_WIDTH / 2, int(res * MESH_WIDTH) + 1),
        y=np.linspace(-MESH_HEIGHT / 2, MESH_HEIGHT / 2, int(res * MESH_HEIGHT) + 1),
    )
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    cs = mw.CrossSection.from_cell(cell=cells[0], env=mw.Environment(wl=float(wl), T=25.0))
    return float(np.real(compute_modes_meep(cs, num_modes=1)[0].neff))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    neff = np.zeros((len(RESOLUTIONS), len(WAVELENGTHS)))
    for i, res in enumerate(RESOLUTIONS):
        for j, wl in enumerate(WAVELENGTHS):
            neff[i, j] = solve_neff(wl, res)
            print(f"[gap3] res={res:5.1f} wl={wl:.2f} neff={neff[i, j]:.6f}", flush=True)

    print("\n res |" + "".join(f"  neff@{w:.2f}" for w in WAVELENGTHS))
    for i, res in enumerate(RESOLUTIONS):
        print(f"{res:5.1f}|" + "".join(f"  {neff[i, j]:9.6f}"
                                      for j in range(len(WAVELENGTHS))))

    # One dphi table per route, because the two routes hand a different amount of
    # straight length to the analytic phase (see LEVER_ARMS above).
    for label, arm in LEVER_ARMS.items():
        print(f"\n[gap3] dphi, {label}: analytic straight = {arm:.2f} um")
        print(" res |" + "".join(f"   dphi@{w:.2f}(rad)" for w in WAVELENGTHS)
              + "   |  worst as frac of an FSR")
        for i, res in enumerate(RESOLUTIONS):
            if not i:
                continue
            d = [2 * np.pi * (neff[i, j] - neff[i - 1, j]) * arm / WAVELENGTHS[j]
                 for j in range(len(WAVELENGTHS))]
            worst = max(abs(x) for x in d) / (2 * np.pi)
            print(f"{res:5.1f}|" + "".join(f"   {x:+12.4f}" for x in d)
                  + f"   |  {worst:8.4f}")

    print("\n[gap3] dphi is the CHANGE in round-trip phase from the previous "
          "resolution row. A full 2*pi means the resonance comb has slid by an "
          "entire FSR between those two meshes.")
    print("[gap3] The cascaded row is the one that describes the shipped default. "
          "Its smaller lever arm is NOT a free improvement -- the 5.00 um it "
          "dropped now rides inside the FDTD-measured bend coefficient, whose own "
          "resolution error is what the res_sweep study measures separately.")
    np.savez(os.path.join(OUT, "neff_convergence.npz"),
             resolutions=np.array(RESOLUTIONS), wavelengths=np.array(WAVELENGTHS),
             neff=neff, L_outer=L_OUTER,
             lever_arm_labels=np.array(list(LEVER_ARMS.keys())),
             lever_arms=np.array(list(LEVER_ARMS.values())))
    print(f"[gap3] saved -> {OUT}/neff_convergence.npz")
