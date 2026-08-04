"""Targeted follow-up to mmi_1x4_sin_cross_section.py's step 4 sweep --
validation_markdown/mmi_1x4_sin300_validation.md.

That sweep found the single-mode width bracket for 300 nm Si3N4 is only
0.2-0.35 um and weakly confined there (margin above n_clad ~0.003 or less,
sz not converged to 1e-4 even at 10 um). User's call, given that finding:
reconsider the design width at 1.1 um instead of chasing the narrow
single-mode window. This measures what is ACTUALLY there at 1.1 um -- guided
mode count, each mode's n_eff/confinement margin, and sz convergence at this
width -- rather than assuming.

Same cross-section builder, wavelength (1.55 um C-band), resolution (40
pts/um) and compute_modes_tidy3d backend as the main sweep script (switched
there after compute_modes_meep was found to silently return spurious
box-mode answers for this material's weaker index contrast -- see that
script's docstring for the full chain of evidence).
"""
import numpy as np
import meow as mw
from meow import compute_modes_tidy3d

from beef import extrude_polygons
from beef.materials.stack import MaterialStack, sellmeier_index

WG_HEIGHT = 0.30
WL = 1.55
RES = 40.0
MESH_WIDTH = 8.0
N_CLAD = float(sellmeier_index("sio2", WL))
WIDTH = 1.1
SZ_CANDIDATES = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
MAX_MODES_REQUESTED = 8


def axis_grid(half_span, res):
    n_half = int(round(half_span * res))
    step = 1.0 / res
    return step * np.arange(-n_half, n_half + 1)


def cross_section(width, sz, wl=WL):
    length = 1.0
    wg = [(0, -width/2), (length, -width/2), (length, width/2), (0, width/2)]
    ox = [(0, -MESH_WIDTH/2), (length, -MESH_WIDTH/2),
          (length, MESH_WIDTH/2), (0, MESH_WIDTH/2)]
    gds = {(2, 0): [wg], (1, 0): [ox]}
    stack = MaterialStack(wg_height=WG_HEIGHT, sz=sz, core_kind="si3n4")
    structs = extrude_polygons(gds, stack.extrusion_rules)
    mesh = mw.Mesh2D(x=axis_grid(MESH_WIDTH/2, RES), y=axis_grid(sz/2, RES))
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=float(wl), T=25.0)
    return mw.CrossSection.from_cell(cell=cells[0], env=env)


print("=" * 78)
print(f"Width={WIDTH} um follow-up -- 300 nm Si3N4, lambda={WL} um, n_clad={N_CLAD:.6f}")
print("=" * 78)

n_eff_top = np.zeros(len(SZ_CANDIDATES))
n_guided_list = []
for i, sz in enumerate(SZ_CANDIDATES):
    cs = cross_section(WIDTH, sz)
    modes = compute_modes_tidy3d(cs, num_modes=MAX_MODES_REQUESTED)
    neffs = sorted((float(np.real(m.neff)) for m in modes), reverse=True)
    guided = [n for n in neffs if n > N_CLAD]
    n_eff_top[i] = neffs[0]
    n_guided_list.append(guided)
    margins = [g - N_CLAD for g in guided]
    print(f"  sz={sz:5.2f} um  ->  {len(guided)} guided: "
          f"{[f'{g:.6f} (margin {m:+.5f})' for g, m in zip(guided, margins)]}")

d = np.abs(np.diff(n_eff_top))
print(f"\n  Successive |dn_eff| (top mode): {[f'{x:.2e}' for x in d]}")
for tol in (1e-4, 5e-4, 1e-3):
    idxs = np.flatnonzero(d < tol)
    if idxs.size:
        i = idxs[0]
        print(f"  clears {tol:.0e} first going {SZ_CANDIDATES[i]:.1f}->{SZ_CANDIDATES[i+1]:.1f} um")
    else:
        print(f"  never clears {tol:.0e} within the swept range")

print("\n" + "=" * 78)
print("END")
print("=" * 78)
