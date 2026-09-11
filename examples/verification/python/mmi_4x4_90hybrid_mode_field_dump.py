"""4x4 MMI 90-degree hybrid -- Step 8 Check 2 diagnostic: raw per-band meep
eigenmode field dump at the body cross-section, WITHOUT the duplicate-mode
guard.

validation_markdown/mmi_4x4_90_hybrid_validation.md Step 8, Check 2. Job
13060190 (audit_interface_modes.py mmi_4x4_90hybrid, N_MODES=16) FAILED:

    RuntimeError: FDTDSolver: could not find 16 distinct eigenmodes within
    24 band_num attempts at this cross-section -- meep/MPB returned a
    duplicate more times than this guard tolerates.

That guard (FDTDSolver._iter_distinct_bands(), beef/solvers/fdtd.py) SKIPS
any band_num whose field overlaps an already-accepted one above 0.9 and
raises once it runs out of attempts -- useful for keeping a real interface
basis rank-full, but it also throws away exactly the evidence needed to see
WHY duplicates are showing up here: the skipped fields are never returned to
the caller. This script re-implements the same raw per-band walk (same
_build_section_mode_sim() preamble, same k-seed laddering, same
guided-edge stop) but keeps EVERY raw band_num's result, plots each one's
transverse |Ey|^2 intensity in the y-z plane, and computes the full N x N
self-overlap matrix so actual duplicate pairs (not just "a guard tripped")
are visible directly.

No FDTDSolver private method is bypassed via subclassing (see the repo's own
"never subclass solvers to instrument" lesson) -- this calls
FDTDSolver._build_section_mode_sim() once (the same non-guarded preamble
section_mode_neff()/section_mode_fields() both already share) and then talks
to meep's `sim.get_eigenmode()` directly, exactly as
FDTDSolver._iter_distinct_bands() does internally, just without its
skip/raise logic.

Run:
    python examples/verification/python/mmi_4x4_90hybrid_mode_field_dump.py
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from beef import MMI, FDTDSolver, SOI_220, field_overlap_matrix

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "mmi_4x4_90hybrid_mode_field_dump")
os.makedirs(output_dir, exist_ok=True)

# Pinned to examples/mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_hybrid.py and to
# audit_interface_modes.py's _mmi_4x4_90hybrid() config exactly -- the same
# cross-section job 13060190 failed on.
MM_WIDTH = 8.5
WAVEGUIDE_WIDTH = 0.35
TAPER_LENGTH = 4.0
TAPER_WIDTH = 1.30
MM_LENGTH = 159.30
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
WL = 1.31
STACK = SOI_220
SECTION = (-1.0, 1.0)       # body-interior cut, same convention as every
                            # other device's own Check 2 audit config
SIM_WIDTH = MM_WIDTH + 3.0
SIM_HEIGHT = 2.0
RES = 40.0
DPML = 0.5

N_MODES_REQUESTED = 16      # what the real hybrid cascade asks for
MAX_BAND_TRY = 30           # walk this far past N_MODES_REQUESTED looking
                            # for the guided edge / enough distinct modes to
                            # see the failure clearly; audit_interface_modes
                            # itself only ever tried 24 (16 + its own +8 slack)

_KPOINT_SEED_NEFF = 2.5     # same starting k-seed FDTDSolver._iter_distinct_bands() uses

# Field sample grid for both the plots and the overlap matrix -- 0.05 um
# spacing, matching this device's own res=40 (1/40=0.025, so this is 2x
# coarser than the structure's own grid, plenty for an aggregate overlap
# integral and for visually resolving mode lobes on an 8.5 um body).
Y_GRID = np.linspace(-SIM_WIDTH / 2, SIM_WIDTH / 2, int(SIM_WIDTH / 0.05) + 1)
Z_GRID = np.linspace(-SIM_HEIGHT / 2, SIM_HEIGHT / 2, int(SIM_HEIGHT / 0.05) + 1)
DY = float(Y_GRID[1] - Y_GRID[0])
DZ = float(Z_GRID[1] - Z_GRID[0])


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
    print(f"[dump] mm_width={MM_WIDTH}, section={SECTION}, wl={WL} um, "
          f"n_cladding={n_clad:.6f}, requesting up to {N_MODES_REQUESTED} "
          f"modes (walking up to {MAX_BAND_TRY} raw band_num attempts)")

    x0, x1 = SECTION
    sim, where, mea_c, mea_w = solver._build_section_mode_sim(
        mmi, WL, x0, x1, monitor_center=0.0, monitor_size=SIM_WIDTH)

    f = 1.0 / WL
    neffs = []
    fields_Ey = []
    fields_Ez = []
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
        print(f"[dump] band_num={band_num:3d}  n_eff={n:.6f}  "
              f"{'guided' if guided else 'UNGUIDED -- stopping walk'}")
        if not guided:
            break

        ey = np.array([[em.amplitude(mp.Vector3(0, float(y), float(z)), mp.Ey)
                        for z in Z_GRID] for y in Y_GRID])
        ez = np.array([[em.amplitude(mp.Vector3(0, float(y), float(z)), mp.Ez)
                        for z in Z_GRID] for y in Y_GRID])
        neffs.append(n)
        fields_Ey.append(ey)
        fields_Ez.append(ez)

        intensity = np.abs(ey) ** 2 + np.abs(ez) ** 2
        fig, ax = plt.subplots(figsize=(6, 3))
        pc = ax.pcolormesh(Y_GRID, Z_GRID, intensity.T, shading="auto", cmap="inferno")
        fig.colorbar(pc, ax=ax, label="|Ey|^2 + |Ez|^2")
        ax.set_xlabel("y (um, device transverse)")
        ax.set_ylabel("z (um, vertical)")
        ax.set_title(f"raw band_num={band_num}, n_eff={n:.6f}"
                     f"{' (guided)' if guided else ' (UNGUIDED)'}")
        ax.set_aspect("equal")
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, f"mode_{band_num:02d}_neff{n:.4f}.png"), dpi=150)
        plt.close(fig)

    sim.reset_meep()

    n_found = len(neffs)
    print(f"\n[dump] {n_found} guided raw bands found before the walk stopped "
          f"(requested {N_MODES_REQUESTED})")

    fields = dict(Ey=np.stack(fields_Ey), Ez=np.stack(fields_Ez))
    eta = field_overlap_matrix(fields, fields, dy=DY, dz=DZ)

    # Table: band, n_eff, best off-diagonal overlap + which band it is with
    # (the direct duplicate signal -- >0.9 is FDTDSolver's own threshold).
    lines = []
    lines.append("| band_num | n_eff | best off-diag overlap | vs band | duplicate? |")
    lines.append("|---|---|---|---|---|")
    for i in range(n_found):
        row = eta[i].copy()
        row[i] = -1.0
        j = int(np.argmax(row))
        best = row[j]
        dup = "YES" if best > 0.9 else ("maybe" if best > 0.5 else "")
        lines.append(f"| {i+1} | {neffs[i]:.6f} | {best:.6f} | {j+1} | {dup} |")
    table = "\n".join(lines)
    print("\n" + table)

    with open(os.path.join(output_dir, "mode_overlap_table.md"), "w") as fh:
        fh.write(f"# {os.path.basename(__file__)} -- self-overlap table\n\n")
        fh.write(f"mm_width={MM_WIDTH}, section={SECTION}, wl={WL}, "
                 f"n_cladding={n_clad:.6f}, {n_found} guided raw bands found\n\n")
        fh.write(table + "\n")

    np.savez(os.path.join(output_dir, "mode_dump.npz"),
             neff=np.array(neffs), eta=eta, y_grid=Y_GRID, z_grid=Z_GRID)

    print(f"\n[dump] plots + table + npz -> {output_dir}/")


if __name__ == "__main__":
    main()
