"""mmi_1x2_validation.md gap 2 — is absolute n_eff the wrong thing to converge?

Same question, same recipe as verify_mmi_2x2_gap2_dbeta.py -- read that
script's docstring for the full argument (item 11's coupler lesson: sweep
the quantity the physics depends on, not an arbitrary one).

This device is flagged as the worst-converging of the three audited
(mmi_1x2_validation.md gap 2: still +0.0083 at res 100, non-monotone,
3x the 2x2's step over the same interval) and nothing has tested whether
that matters -- its self-imaging length (mm_length=4.5) is set by the beat
length Δβ_01 between its two lowest guided modes, same as any MMI, not by
either mode's own n_eff.

Cost: same as the 2x2 script -- EMESolver.compute_modes() solves the lowest
N bands in one call, so N_MODES=2 is cheap regardless of the device's real
N_MODES=8. Batch job, not login-node (2x2's smoke test measured minutes per
point on this class of cross-section).

Run:
    OMP_NUM_THREADS=1 python verify_mmi_1x2_gap2_dbeta.py
"""
from __future__ import annotations

import time

import numpy as np
import meow as mw

from beef import MMI, EMESolver, SOI_160, extrude_polygons

WL = 1.31
STACK = SOI_160
SECTION = (-1.0, 1.0)        # 1.25 um clear of each taper seam (mm_length=4.5)
SIM_W, SIM_H = 5.0, 1.5
N_MODES = 2                  # only the pair the beat length depends on
RESOLUTIONS = (40.0, 60.0, 80.0, 100.0, 120.0)


def _mmi():
    return MMI(
        n_inputs=1, n_outputs=2,
        mm_length=4.5, mm_width=2.0, waveguide_width=0.35,
        taper_length=1.0, taper_width=0.8, waveguide_offset=1.0,
        waveguide_height=0.16,
    )


def neff_pair(mmi, res):
    eme = EMESolver(
        num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
        oxide_width=SIM_W, oxide_height=SIM_H,
        res=float(res), material_stack=STACK, mode_solver="meep",
    )
    x0, x1 = SECTION
    gds = mmi.section_gds(x0, x1, eme.oxide_width)
    cells = mw.create_cells(
        structures=extrude_polygons(gds, eme.material_stack.extrusion_rules),
        mesh=eme._build_mesh(), Ls=[x1 - x0],
    )
    cs = mw.CrossSection.from_cell(
        cell=cells[0], env=mw.Environment(wl=WL, T=eme.temperature))
    modes = eme.compute_modes(cs, N_MODES)
    return np.array([float(np.real(m.neff)) for m in modes])


def main():
    mmi = _mmi()
    print(f"[gap2] mmi_1x2 body, section={SECTION}, wl={WL}, N_MODES={N_MODES}")
    print(f"{'res':>5}  {'n_eff[0]':>10}  {'n_eff[1]':>10}  {'dbeta (rad/um)':>14}  "
          f"{'L_pi (um)':>10}  {'step':>10}")

    prev_dbeta = None
    for res in RESOLUTIONS:
        t0 = time.time()
        neff = neff_pair(mmi, res)
        dbeta = 2 * np.pi / WL * (neff[0] - neff[1])
        l_pi = np.pi / dbeta
        step = "" if prev_dbeta is None else f"{dbeta - prev_dbeta:+.5f}"
        print(f"{res:5.0f}  {neff[0]:10.6f}  {neff[1]:10.6f}  {dbeta:14.5f}  "
              f"{l_pi:10.4f}  {step:>10}   ({(time.time()-t0):.1f}s)")
        prev_dbeta = dbeta


if __name__ == "__main__":
    main()
