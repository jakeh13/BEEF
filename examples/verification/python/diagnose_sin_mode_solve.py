"""Diagnostic: why compute_modes_meep() isn't returning a guided mode for the
300 nm Si3N4 strip at any (width, sz) tested so far in
mmi_1x4_sin_cross_section.py (all candidates returned n_eff < n_clad).

Sanity check first (a Si control at the identical box, which is known-good
territory for this exact solver everywhere else in this repo), then several
angles on the Si3N4 case: more requested modes (is the bound mode present but
mis-ranked?), a much wider/stronger structure, and higher resolution -- to
separate "solver search strategy problem" from "resolution/subpixel artifact"
from "genuinely no bound mode there".

Analytic cross-check (symmetric-slab TE0, infinite width, done separately):
n_eff = 1.673911, cladding decay length 0.291 um -- a real, tightly-confined
bound mode should exist for the vertical structure alone. Not reproduced by
compute_modes_meep at any (width<=2.4um, sz<=5um) tested manually so far.
"""
import numpy as np
import meow as mw
from meow import compute_modes_meep

from beef import extrude_polygons
from beef.materials.stack import MaterialStack, sellmeier_index

WL = 1.55
N_CLAD = float(sellmeier_index("sio2", WL))
N_SI = float(sellmeier_index("si", WL))
N_SIN = float(sellmeier_index("si3n4", WL))


def axis_grid(half_span, res):
    n_half = int(round(half_span * res))
    step = 1.0 / res
    return step * np.arange(-n_half, n_half + 1)


def cross_section(core_kind, width, sz, mesh_width, res, wl=WL):
    length = 1.0
    wg = [(0, -width/2), (length, -width/2), (length, width/2), (0, width/2)]
    ox = [(0, -mesh_width/2), (length, -mesh_width/2),
          (length, mesh_width/2), (0, mesh_width/2)]
    gds = {(2, 0): [wg], (1, 0): [ox]}
    stack = MaterialStack(wg_height=0.30, sz=sz, core_kind=core_kind)
    structs = extrude_polygons(gds, stack.extrusion_rules)
    mesh = mw.Mesh2D(x=axis_grid(mesh_width/2, res), y=axis_grid(sz/2, res))
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=float(wl), T=25.0)
    return mw.CrossSection.from_cell(cell=cells[0], env=env)


def report(label, core_kind, width, sz, mesh_width, res, num_modes, n_clad=N_CLAD):
    cs = cross_section(core_kind, width, sz, mesh_width, res)
    modes = compute_modes_meep(cs, num_modes=num_modes)
    neffs = sorted((float(np.real(m.neff)) for m in modes), reverse=True)
    guided = [n for n in neffs if n > n_clad]
    print(f"\n{label}")
    print(f"  core={core_kind} width={width} sz={sz} mesh_width={mesh_width} "
          f"res={res} num_modes={num_modes} n_clad={n_clad:.6f}")
    print(f"  neffs: {[f'{n:.5f}' for n in neffs]}")
    print(f"  -> {len(guided)} guided: {[f'{n:.5f}' for n in guided]}")
    return neffs, guided


print("=" * 78)
print("DIAGNOSTIC: compute_modes_meep() on the 300 nm Si3N4 strip")
print(f"n_clad(SiO2)={N_CLAD:.6f}  n_core(Si3N4)={N_SIN:.6f}  n_core(Si)={N_SI:.6f}")
print("Analytic infinite-width symmetric-slab TE0 for 300 nm Si3N4/SiO2: "
      "n_eff=1.673911, cladding decay length 0.291 um")
print("=" * 78)

# 1. Control: SAME box, Si core instead of Si3N4. Known-good territory --
#    directional_coupler_mode.py runs this exact solver on Si structures
#    routinely. If this ALSO fails to guide, the bug is generic to this
#    box/solver combination, not specific to Si3N4's index range.
report("1. CONTROL -- Si core, otherwise identical box", "si",
       width=1.0, sz=2.0, mesh_width=8.0, res=40.0, num_modes=3)

# 2. Si3N4, same as the failing case, but request MANY more modes -- is the
#    bound mode present in the spectrum but simply not among the top 6?
report("2. Si3N4, num_modes=25 -- is the bound mode just mis-ranked?", "si3n4",
       width=1.0, sz=2.0, mesh_width=8.0, res=40.0, num_modes=25)

# 3. Si3N4, much wider and taller -- deep into unambiguously-multimode
#    territory, where the guided-mode manifold is dense and the true
#    ground-state mode is index-separated from radiation modes by more.
report("3. Si3N4, wide/strong (width=4.0, sz=6.0) -- easier target for the search?",
      "si3n4", width=4.0, sz=6.0, mesh_width=10.0, res=40.0, num_modes=6)

# 4. Si3N4, higher resolution -- rules out a subpixel-averaging/resolution
#    artifact on the thin 300 nm (12-pixel at res=40) core.
report("4. Si3N4, res=80 (24 pixels across the 300 nm core)", "si3n4",
      width=1.0, sz=2.0, mesh_width=8.0, res=80.0, num_modes=6)

print("\n" + "=" * 78)
print("END DIAGNOSTIC")
print("=" * 78)
