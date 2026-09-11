"""4x4 MMI 90-degree hybrid -- Step 8 Check 2 diagnostic, part 2: does the
production duplicate-guard's own coarse diagnostic grid cause a false
positive on this device's wide slab?

validation_markdown/mmi_4x4_90_hybrid_validation.md Step 8, Check 2.
mmi_4x4_90hybrid_mode_field_dump.py (job 13062141) found ZERO duplicates
among the first 30 raw bands at this cross-section, sampled on a fine
231x41-point (0.05 um) grid -- yet the production audit
(audit_interface_modes.py, job 13060190) FAILED with
``FDTDSolver._iter_distinct_bands()``'s duplicate guard, which uses its own
FIXED 16x8-point diagnostic grid (``_DUP_BAND_DIAG_NY``/``_DUP_BAND_DIAG_NZ``
in beef/solvers/fdtd.py) regardless of caller resolution. 8 points across
the whole sim_height=2.0 um gives ~0.25 um z-spacing against this device's
0.22 um-thick core -- under one cell resolving the guided layer. Leading
hypothesis: that under-resolution aliases two genuinely-different,
similarly-shaped high-order modes into an apparent >0.9 overlap.

This script re-walks the SAME cross-section, SAME k-seed-laddering
algorithm as both the guard and the prior field-dump script, but samples
EACH raw band's field at BOTH grids in the same pass -- the guard's exact
coarse 16x8 grid (same construction as
``FDTDSolver._iter_distinct_bands()``: cell centers of ``_DUP_BAND_DIAG_NY``
edges over the monitor width, ``_DUP_BAND_DIAG_NZ`` edges over
``self._mon_sz``) and a fine 231x41 grid matching the prior dump script --
so the two overlap matrices are computed from the literal same live
eigenmode objects, not compared across two separate runs' worth of
(possibly not bit-identical) MPB convergence.

Reports, per band pair, whether the COARSE grid's overlap crosses the 0.9
duplicate threshold while the FINE grid's does not -- direct confirmation
(or refutation) of the aliasing hypothesis, not an inference from where an
exception fired.

Run:
    python examples/verification/python/mmi_4x4_90hybrid_mode_dup_grid_compare.py
"""
from __future__ import annotations

import os

import numpy as np

from beef import MMI, FDTDSolver, SOI_220, field_overlap_matrix

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "mmi_4x4_90hybrid_mode_dup_grid_compare")
os.makedirs(output_dir, exist_ok=True)

# Identical geometry/section/solver config to mmi_4x4_90hybrid_mode_field_dump.py
# and audit_interface_modes.py's _mmi_4x4_90hybrid() -- same cross-section as
# both the failing audit and the clean fine-grid dump.
MM_WIDTH = 8.5
WAVEGUIDE_WIDTH = 0.35
TAPER_LENGTH = 4.0
TAPER_WIDTH = 1.30
MM_LENGTH = 159.30
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
WL = 1.31
STACK = SOI_220
SECTION = (-1.0, 1.0)
SIM_WIDTH = MM_WIDTH + 3.0
SIM_HEIGHT = 2.0
RES = 40.0
DPML = 0.5

MAX_BAND_TRY = 27           # covers 1-24 (the guard's own attempt budget)
                            # plus a few past it, where mode_field_dump.py
                            # found spacing starting to collapse
_KPOINT_SEED_NEFF = 2.5

# The guard's OWN diagnostic grid -- must match beef/solvers/fdtd.py's
# _DUP_BAND_DIAG_NY / _DUP_BAND_DIAG_NZ exactly, not a value chosen here.
_DUP_BAND_DIAG_NY = 16
_DUP_BAND_DIAG_NZ = 8

# The fine grid mmi_4x4_90hybrid_mode_field_dump.py already used, for a
# same-run cross-check against that job's own (separately-solved) result.
FINE_Y = np.linspace(-SIM_WIDTH / 2, SIM_WIDTH / 2, int(SIM_WIDTH / 0.05) + 1)
FINE_Z = np.linspace(-SIM_HEIGHT / 2, SIM_HEIGHT / 2, int(SIM_HEIGHT / 0.05) + 1)
FINE_DY = float(FINE_Y[1] - FINE_Y[0])
FINE_DZ = float(FINE_Z[1] - FINE_Z[0])


def main():
    import meep as mp

    mmi = MMI(
        n_inputs=4, n_outputs=4,
        mm_length=MM_LENGTH, mm_width=MM_WIDTH, waveguide_width=WAVEGUIDE_WIDTH,
        taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
        waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
        use_neighbor_aperture=True,
    )
    solver = FDTDSolver(
        res=RES, sim_width=SIM_WIDTH, sim_height=SIM_HEIGHT, dpml=DPML,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=STACK, decay_tol=1e-11,
    )
    n_clad = float(solver.material_stack.sio2_index(WL))

    x0, x1 = SECTION
    sim, where, mea_c, mea_w = solver._build_section_mode_sim(
        mmi, WL, x0, x1, monitor_center=0.0, monitor_size=SIM_WIDTH)

    # Exactly FDTDSolver._iter_distinct_bands()'s own diag_y/diag_z
    # construction (beef/solvers/fdtd.py) -- cell centers, not edges.
    diag_y_edges = np.linspace(mea_c - mea_w / 2, mea_c + mea_w / 2, _DUP_BAND_DIAG_NY + 1)
    coarse_y = 0.5 * (diag_y_edges[:-1] + diag_y_edges[1:])
    diag_z_edges = np.linspace(-solver._mon_sz / 2, solver._mon_sz / 2, _DUP_BAND_DIAG_NZ + 1)
    coarse_z = 0.5 * (diag_z_edges[:-1] + diag_z_edges[1:])
    print(f"[compare] coarse grid: {_DUP_BAND_DIAG_NY} x {_DUP_BAND_DIAG_NZ}, "
          f"y-spacing={coarse_y[1]-coarse_y[0]:.4f} um, "
          f"z-spacing={coarse_z[1]-coarse_z[0]:.4f} um (core height 0.22 um)")
    print(f"[compare] fine grid: {len(FINE_Y)} x {len(FINE_Z)}, "
          f"y-spacing={FINE_DY:.4f} um, z-spacing={FINE_DZ:.4f} um")

    f = 1.0 / WL
    neffs = []
    coarse_Ey, coarse_Ez = [], []
    fine_Ey, fine_Ez = [], []
    last_raw_neff = None

    for band_num in range(1, MAX_BAND_TRY + 1):
        seed_neff = last_raw_neff if last_raw_neff is not None else _KPOINT_SEED_NEFF
        em = sim.get_eigenmode(
            f, mp.X, where, band_num, mp.Vector3(f * seed_neff, 0, 0),
            parity=mp.NO_PARITY,
        )
        n = float(np.real(em.k.x)) / f
        last_raw_neff = n
        guided = n > n_clad
        print(f"[compare] band_num={band_num:3d}  n_eff={n:.6f}  "
              f"{'guided' if guided else 'UNGUIDED -- stopping walk'}")
        if not guided:
            break
        neffs.append(n)

        cy = np.array([[em.amplitude(mp.Vector3(0, float(y), float(z)), mp.Ey)
                        for z in coarse_z] for y in coarse_y])
        cz = np.array([[em.amplitude(mp.Vector3(0, float(y), float(z)), mp.Ez)
                        for z in coarse_z] for y in coarse_y])
        coarse_Ey.append(cy)
        coarse_Ez.append(cz)

        fy = np.array([[em.amplitude(mp.Vector3(0, float(y), float(z)), mp.Ey)
                        for z in FINE_Z] for y in FINE_Y])
        fz = np.array([[em.amplitude(mp.Vector3(0, float(y), float(z)), mp.Ez)
                        for z in FINE_Z] for y in FINE_Y])
        fine_Ey.append(fy)
        fine_Ez.append(fz)

    sim.reset_meep()
    n_found = len(neffs)
    print(f"\n[compare] {n_found} guided raw bands found")

    coarse_fields = dict(Ey=np.stack(coarse_Ey), Ez=np.stack(coarse_Ez))
    fine_fields = dict(Ey=np.stack(fine_Ey), Ez=np.stack(fine_Ez))
    coarse_dy = float(coarse_y[1] - coarse_y[0])
    coarse_dz = float(coarse_z[1] - coarse_z[0])
    eta_coarse = field_overlap_matrix(coarse_fields, coarse_fields, dy=coarse_dy, dz=coarse_dz)
    eta_fine = field_overlap_matrix(fine_fields, fine_fields, dy=FINE_DY, dz=FINE_DZ)

    lines = []
    lines.append("| band_num | n_eff | coarse best off-diag | vs band | fine best off-diag | vs band | coarse flags dup, fine does not? |")
    lines.append("|---|---|---|---|---|---|---|")
    any_confirmed = False
    for i in range(n_found):
        rc = eta_coarse[i].copy(); rc[i] = -1.0
        jc = int(np.argmax(rc)); bc = rc[jc]
        rf = eta_fine[i].copy(); rf[i] = -1.0
        jf = int(np.argmax(rf)); bf = rf[jf]
        flag = bc > 0.9 and bf <= 0.9
        any_confirmed = any_confirmed or flag
        lines.append(f"| {i+1} | {neffs[i]:.6f} | {bc:.6f} | {jc+1} | {bf:.6f} | {jf+1} | {'YES' if flag else ''} |")
    table = "\n".join(lines)
    print("\n" + table)
    print(f"\n[compare] HYPOTHESIS {'CONFIRMED' if any_confirmed else 'NOT CONFIRMED'}: "
          f"{'at least one' if any_confirmed else 'no'} band pair flags a coarse-grid-only duplicate")

    with open(os.path.join(output_dir, "grid_compare_table.md"), "w") as fh:
        fh.write(f"# {os.path.basename(__file__)} -- coarse vs fine grid overlap comparison\n\n")
        fh.write(f"coarse grid: {_DUP_BAND_DIAG_NY}x{_DUP_BAND_DIAG_NZ} "
                 f"(dy={coarse_dy:.4f}, dz={coarse_dz:.4f} um); "
                 f"fine grid: {len(FINE_Y)}x{len(FINE_Z)} (dy={FINE_DY:.4f}, dz={FINE_DZ:.4f} um)\n\n")
        fh.write(table + "\n\n")
        fh.write(f"HYPOTHESIS {'CONFIRMED' if any_confirmed else 'NOT CONFIRMED'}\n")

    np.savez(os.path.join(output_dir, "grid_compare.npz"),
             neff=np.array(neffs), eta_coarse=eta_coarse, eta_fine=eta_fine,
             coarse_y=coarse_y, coarse_z=coarse_z, fine_y=FINE_Y, fine_z=FINE_Z)
    print(f"\n[compare] table + npz -> {output_dir}/")


if __name__ == "__main__":
    main()
