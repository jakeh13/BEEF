"""Probe: is eval_from_basis()'s 69 s/point dominated by interface S-matrices?

Step 6's analytic grid (mmi_2x2_waist_body_validation.md) reuses ONE mode set
across the whole (s, Lt, Lw) grid, yet each point still costs ~69 s -- 384
points/Ww = 7.4 h, which is what timed out job 12630600.

meow's compute_s_matrix_sax() takes `interfaces_fn` as an argument and calls
it as interfaces_fn(modes). The interface S-matrices are a pure function of the
modes, and the modes are IDENTICAL for every grid point (that is the premise
the whole free-sweep rests on). Only compute_propagation_s_matrices() consumes
cell_lengths. So a memoised interfaces_fn is EXACT, not an approximation.

This probe measures the split and checks exactness. It does NOT change the
sweep script.
"""
import os
import sys
import time

import numpy as np
import meow as mw
from meow.eme.interface import compute_interface_s_matrices
from meow.eme.propagation import compute_propagation_s_matrices

from beef import MMI, EMESolver, SOI_220

WL = 1.5475
MM_WIDTH = 2.8
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 0.45, 1.3, 4.0
WAVEGUIDE_OFFSET, WAVEGUIDE_HEIGHT = 1.5, 0.22
STRAIGHT_LENGTH = 0.5
TAPER_CELLS = 15
C_TRANS = 8
RES = 40.0
OXW, OXH = 5.0, 1.5
N_PORT_MODES, N_TAPER_MODES, N_BODY_MODES = 1, 7, 15
WW, LW_REF, LT_REF, S_REF = 2.24, 8.0, 4.0, 1.0


def regions(lw, lt, s):
    return (
        [STRAIGHT_LENGTH, TAPER_LENGTH, s, lt, lw, lt, s, TAPER_LENGTH,
         STRAIGHT_LENGTH],
        [1, TAPER_CELLS, 1, C_TRANS, 1, C_TRANS, 1, TAPER_CELLS, 1],
        [N_PORT_MODES, N_TAPER_MODES, N_BODY_MODES, N_BODY_MODES,
         N_BODY_MODES, N_BODY_MODES, N_BODY_MODES, N_TAPER_MODES,
         N_PORT_MODES],
    )


def cell_lengths_for(cells, s, lt, lw):
    """Per-cell lengths for the 9-region layout (same rule as the sweep)."""
    lens, ncells, _ = regions(lw, lt, s)
    out = []
    for L, n in zip(lens, ncells):
        out.extend([L / n] * n)
    assert len(out) == len(cells), (len(out), len(cells))
    return out


lens, ncells, modes_per = regions(LW_REF, LT_REF, S_REF)
solver = EMESolver(num_modes=modes_per, num_cells=ncells, lengths=lens,
                   oxide_width=OXW, oxide_height=OXH, res=RES,
                   material_stack=SOI_220, mode_solver="tidy3d")
mmi = MMI(n_inputs=2, n_outputs=2, mm_length=2 * S_REF + 2 * LT_REF + LW_REF,
          mm_width=MM_WIDTH, waveguide_width=WAVEGUIDE_WIDTH,
          taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
          waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=WAVEGUIDE_HEIGHT,
          body_waist_width=WW, body_waist_length=LW_REF,
          body_transition_length=LT_REF)

wavelengths = np.array([WL])
mpc = solver._modes_per_cell()
cells = solver._build_cells(mmi, side="bar", wavelengths=wavelengths,
                            input_idx=1, output_idx=1)
print(f"{len(cells)} cells, modes-per-cell {mpc}", flush=True)
env = mw.Environment(wl=WL, T=solver.temperature)

t0 = time.time()
css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, mpc)]
print(f"mode solve (1 wl, 1 port): {time.time()-t0:.1f}s", flush=True)

cl = cell_lengths_for(cells, S_REF, LT_REF, LW_REF)

t0 = time.time()
S_ref, pm_ref = mw.compute_s_matrix(modes, cells=cells, cell_lengths=cl)
t_full = time.time() - t0
print(f"full compute_s_matrix          : {t_full:.2f}s", flush=True)

t0 = time.time()
ifaces = compute_interface_s_matrices(modes)
t_if = time.time() - t0
print(f"  interfaces alone             : {t_if:.2f}s", flush=True)

t0 = time.time()
compute_propagation_s_matrices(modes, cells, cell_lengths=cl)
t_pr = time.time() - t0
print(f"  propagations alone           : {t_pr:.2f}s", flush=True)
print(f"  cascade/sax remainder        : {t_full - t_if - t_pr:.2f}s", flush=True)


def cached_ifaces_fn(_modes):
    return ifaces


# exactness + speed of the cached path, at a DIFFERENT (s, lt, lw)
worst = 0.0
times = []
for (sv, ltv, lwv) in [(S_REF, LT_REF, LW_REF), (2.75, 4.0, 8.0),
                       (1.0, 6.0, 10.0), (3.5, 3.0, 6.0)]:
    clv = cell_lengths_for(cells, sv, ltv, lwv)
    t0 = time.time()
    Sc, pmc = mw.compute_s_matrix(modes, cells=cells, cell_lengths=clv,
                                  interfaces_fn=cached_ifaces_fn)
    tc = time.time() - t0
    t0 = time.time()
    Sd, pmd = mw.compute_s_matrix(modes, cells=cells, cell_lengths=clv)
    td = time.time() - t0
    assert pmc == pmd
    d = float(np.max(np.abs(np.asarray(Sc) - np.asarray(Sd))))
    worst = max(worst, d)
    times.append((tc, td))
    print(f"  s={sv:.2f} Lt={ltv:.1f} Lw={lwv:.1f}: cached {tc:.2f}s vs "
          f"uncached {td:.2f}s  |dS| = {d:.3e}", flush=True)

print(f"\nworst |dS| cached vs uncached: {worst:.3e}")
sp = np.mean([t[1] for t in times]) / max(np.mean([t[0] for t in times]), 1e-9)
print(f"mean speedup per S-matrix: {sp:.1f}x")
print(f"VERDICT: {'EXACT' if worst < 1e-13 else 'NOT EXACT -- do not use'}")
