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
# ring_resonator_hybrid.py: L_SIDE=2.0, S_MON, L_FAR=L_C. L_outer is the
# analytic-phase path length; 2*(L_SIDE) + L_FAR with L_C=10.0 is the shipped
# scale. Only used to convert a d(neff) into a d(phase), so the exact value
# just sets the lever arm.
L_OUTER = 2 * 2.0 + 10.0
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

    print("\n res |" + "".join(f"  neff@{w:.2f}" for w in WAVELENGTHS)
          + "".join(f"   dphi@{w:.2f}(rad)" for w in WAVELENGTHS))
    for i, res in enumerate(RESOLUTIONS):
        row = f"{res:5.1f}|" + "".join(f"  {neff[i, j]:9.6f}" for j in range(len(WAVELENGTHS)))
        if i:
            d = [2 * np.pi * (neff[i, j] - neff[i - 1, j]) * L_OUTER / WAVELENGTHS[j]
                 for j in range(len(WAVELENGTHS))]
            row += "".join(f"   {x:+12.4f}" for x in d)
        print(row)

    print("\n[gap3] dphi is the CHANGE in round-trip phase from the previous "
          "resolution row. A full 2*pi means the resonance comb has slid by an "
          "entire FSR between those two meshes.")
    np.savez(os.path.join(OUT, "neff_convergence.npz"),
             resolutions=np.array(RESOLUTIONS), wavelengths=np.array(WAVELENGTHS),
             neff=neff, L_outer=L_OUTER)
    print(f"[gap3] saved -> {OUT}/neff_convergence.npz")
