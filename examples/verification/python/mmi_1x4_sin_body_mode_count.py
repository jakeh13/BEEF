import numpy as np
import meow as mw
from meow import compute_modes_tidy3d
from beef import extrude_polygons
from beef.materials.stack import MaterialStack, sellmeier_index

WG_HEIGHT = 0.30
SZ = 5.0
WL = 1.55
RES = 40.0
N_CLAD = float(sellmeier_index("sio2", WL))

def axis_grid(half_span, res):
    n_half = int(round(half_span * res))
    step = 1.0 / res
    return step * np.arange(-n_half, n_half + 1)

def cross_section(width, mesh_width):
    length = 1.0
    wg = [(0, -width/2), (length, -width/2), (length, width/2), (0, width/2)]
    ox = [(0, -mesh_width/2), (length, -mesh_width/2), (length, mesh_width/2), (0, mesh_width/2)]
    gds = {(2, 0): [wg], (1, 0): [ox]}
    stack = MaterialStack(wg_height=WG_HEIGHT, sz=SZ, core_kind="si3n4")
    structs = extrude_polygons(gds, stack.extrusion_rules)
    mesh = mw.Mesh2D(x=axis_grid(mesh_width/2, RES), y=axis_grid(SZ/2, RES))
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=WL, T=25.0)
    return mw.CrossSection.from_cell(cell=cells[0], env=env)

for width in (12.0,):
    mesh_width = width + 8.0  # generous lateral clearance past the body edges
    cs = cross_section(width, mesh_width)
    modes = compute_modes_tidy3d(cs, num_modes=24)
    neffs = sorted((float(np.real(m.neff)) for m in modes), reverse=True)
    guided = [n for n in neffs if n > N_CLAD]
    print(f"width={width} um: {len(guided)} guided modes of 24 requested")
    for i, n in enumerate(guided):
        print(f"  mode {i}: n_eff={n:.6f}  margin={n-N_CLAD:+.6f}")
