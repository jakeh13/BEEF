"""mmi_2x2_validation.md gap 2 — is absolute n_eff the wrong thing to converge?

The question
------------
Gap 2's existing measurement (job 11500536, `audit_interface_modes.py
--convergence`) shows the MMI body's fundamental-mode n_eff still moving at
res 100 (steps +0.0121 / +0.0068 / +0.0027 at res 40/60/80/100) and explicitly
declines to conclude anything from that alone:

    "Whether that matters depends on which quantity this device's result
    actually depends on ... the coupler's lesson was that sweeping a quantity
    the physics is not sensitive to (T_cross) hid a 7% error in one it is
    (Δβ). Do not simply raise the resolution and assume it is fixed."

`directional_coupler_validation.md` item 11 (archived: `archive/
dc_validation_20260727.md`) is the precedent this is following: absolute
n_eff on that device was converged-looking noise; the coupling physics
actually depends on Δβ = 2π/λ·(n_eff_even − n_eff_odd), which was still 7%
off at the shipped res 40 and needed res 80.

This device is not a coupled pair, but the same logic applies to a different
pair of modes. A 2×2 MMI's self-imaging length (why `mm_length=28.0` was
chosen at all) is set by mode-beat lengths in the multimode body, and the
dominant one is L_π = π/Δβ_01 between the two lowest-order guided modes
(Soldano & Pennings-style MMI theory — the general-interference length is a
rational multiple of L_π, not a function of either mode's own n_eff). So the
analogous "does the physically load-bearing quantity converge" test here is
Δβ_01 = 2π/λ·(n_eff[0] − n_eff[1]), not n_eff[0] alone.

Cost
----
`EMESolver.compute_modes()` solves the lowest N bands in ONE simultaneous
MPB call, so asking for only the 2 modes this needs (not the device's real
N_MODES=14) is cheap by construction — no time-stepping, no FDTD, a single
2-D cross-section mesh per resolution. Login-node smoke run.

Run:
    OMP_NUM_THREADS=1 python verify_mmi_2x2_gap2_dbeta.py
"""
from __future__ import annotations

import time

import numpy as np
import meow as mw

from beef import MMI, EMESolver, SOI_220, extrude_polygons

WL = 1.55
STACK = SOI_220
SECTION = (-1.0, 1.0)        # 13 um clear of either taper seam (mm_length=28)
SIM_W, SIM_H = 5.0, 1.5
N_MODES = 2                  # only the pair the beat length depends on
RESOLUTIONS = (40.0, 60.0, 80.0, 100.0, 120.0)


def _mmi():
    return MMI(
        n_inputs=2, n_outputs=2,
        mm_length=28.0, mm_width=2.8, waveguide_width=0.45,
        taper_length=4.0, taper_width=1.3, waveguide_offset=1.5,
        waveguide_height=0.22,
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
    print(f"[gap2] mmi_2x2 body, section={SECTION}, wl={WL}, N_MODES={N_MODES}")
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
