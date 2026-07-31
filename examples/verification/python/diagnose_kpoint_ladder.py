"""SUPERSEDED by `diagnose_kpoint_ladder_v2.py` -- kept as-is because its
result (3 exact ties / 17 distinct, fixed seed; 1 tie / 19 distinct, naive
ladder) is quoted verbatim in notes.md Solver #18 and CHANGELOG.md
2026-07-30. This script is a flat 20-solve sweep with no accept/skip/retry
and only two policies (fixed, naive ladder); v2 runs a real guarded search
under four policies, including the two follow-ups this script's own findings
called for. Run v2 for anything beyond reproducing the historical numbers
above.

Original docstring, unchanged below --

Does laddering the k seed remove `_iter_distinct_bands()`'s duplicate
eigenmodes? (validation_markdown/fdtd_vs_EME_mode_validation.md item 12.)

The hypothesis
--------------
`_iter_distinct_bands()` (`beef/solvers/fdtd.py`) seeds EVERY `get_eigenmode()`
call from the same fixed guess, `_KPOINT_SEED_NEFF = 2.5` -- band 1 and band 18
alike. `get_eigenmode()` iterates that seed to the true k for the requested
`band_num`; if two nearby band_nums both converge from the same remote guess,
landing on the same eigenvector is the textbook failure, and it is confirmed
real (job 11575516: exact n_eff ties at band_num 9/10 and 13/14 on this exact
cross-section at res 40, before any guard).

The fix under test: instead of always seeding from 2.5, seed band k from band
(k-1)'s own CONVERGED n_eff -- so each search starts near the local mode
density instead of a fixed remote point. If duplicates disappear, both
existing duplicate-band guards (`_iter_distinct_bands()` on the meep side,
`EMESolver.compute_modes()`'s over-request on the meow side) become a backstop
rather than load-bearing, and item 11's grid-resolution concern stops
mattering as much for the meep side.

What it does
------------
Solves band_num 1..N_BANDS at the modal multiplexer's interface 0, res 40 --
the SAME cross-section and resolution where job 11575516 confirmed the raw
(pre-guard) duplicates -- under two seeding policies:

    fixed:   seed_neff = _KPOINT_SEED_NEFF for every band_num (today's code).
    laddered: seed_neff = neff[band_num - 1] (this run's own previous
              CONVERGED value; band 1 still seeds from _KPOINT_SEED_NEFF,
              since there is no "previous" band). Naive on purpose -- chains
              off whatever band_num-1 converged to, duplicate or not, since
              that is the actual sequential call pattern _iter_distinct_bands()
              would make if this were wired in for real.

For each policy: the raw n_eff sequence, exact ties (|dn_eff| < 1e-9 --
threshold-free, does not depend on any field-overlap grid), and the worst
off-diagonal field-overlap among all pairs on a FINE (128x16) sampling grid --
deliberately not the shipped 16x8 grid, so this measurement doesn't inherit
item 11's separate open question about grid resolution.

A risk worth watching for, not just the hoped-for fix: seeding band k directly
at band (k-1)'s converged value could instead make it MORE likely to re-find
band (k-1)'s own eigenvector, if the Newton iteration doesn't move far from its
start. The exact-tie and overlap numbers below distinguish "ladder helped" from
"ladder made it worse" -- both are answers, and the max-iterations-per-solve
printed alongside also flags a seed that produces a harder, not easier, search.
"""
from __future__ import annotations

import numpy as np

from beef import MMI, MMIPort, WaveguideRun, CompoundGeometry, FDTDSolver, SOI_160
from beef.mode_identity import field_overlap_matrix
from beef.solvers.fdtd import _KPOINT_SEED_NEFF

WL = 1.31
N_MODES = 16
N_BANDS = 20          # same as diagnose_dup_band_grid.py -- covers both known
                       # duplicate band_nums (11, 15) at res 40 with headroom
RES = 40.0             # production resolution -- matches the confirmed
                       # pre-guard reproduction (band_num 9/10, 13/14 exact ties)
BODY_STUB = 0.75
SIM_W, SIM_H = 5.0, 1.5
DPML = 0.5
HALF_SPAN = 0.05
FINE_GRID = (128, 16)  # ground-truth overlap grid, deliberately independent
                       # of item 11's shipped-grid question

TAPER_L, WO, BETWEEN, MM_WIDTH = 2.0, 1.0375, 3.0, 2.8
BRANCH_SPECS = [(+WO, 0.33, 0.725), (0.0, 0.35, 1.15), (-WO, 0.35, 0.725)]


def _branch_ports():
    return [MMIPort(y=y, waveguide_width=ww, taper_length=TAPER_L, taper_width=tw)
            for y, ww, tw in BRANCH_SPECS]


def _build():
    mmi1 = MMI(n_inputs=1, n_outputs=3, mm_length=16.92, mm_width=MM_WIDTH,
               inputs=[MMIPort(y=0.0, waveguide_width=0.7, taper_length=TAPER_L,
                               taper_width=1.4)],
               outputs=_branch_ports(), waveguide_height=0.16)
    mmi2 = MMI(n_inputs=3, n_outputs=3, mm_length=34.24, mm_width=MM_WIDTH,
               inputs=_branch_ports(), outputs=_branch_ports(),
               waveguide_height=0.16, use_neighbor_aperture=True)
    connector = WaveguideRun(length=BETWEEN, ports=[(y, ww) for y, ww, _ in BRANCH_SPECS])
    return mmi1, CompoundGeometry([mmi1, connector, mmi2])


def _sample(em, mp, mea_c, mea_w, mon_sz, ny, nz):
    y_edges = np.linspace(mea_c - mea_w / 2, mea_c + mea_w / 2, ny + 1)
    y = 0.5 * (y_edges[:-1] + y_edges[1:])
    z_edges = np.linspace(-mon_sz / 2, mon_sz / 2, nz + 1)
    z = 0.5 * (z_edges[:-1] + z_edges[1:])
    ey = np.array([[em.amplitude(mp.Vector3(0, float(yy), float(zz)), mp.Ey)
                    for zz in z] for yy in y])
    ez = np.array([[em.amplitude(mp.Vector3(0, float(yy), float(zz)), mp.Ez)
                    for zz in z] for yy in y])
    return ey, ez


def _solve_sequence(sim, mp, where, f, mea_c, mea_w, mon_sz, laddered: bool):
    """Returns (neff list, ey list, ez list) for band_num 1..N_BANDS, under
    the fixed or laddered seeding policy."""
    neff, eys, ezs = [], [], []
    seed_neff = _KPOINT_SEED_NEFF
    for band_num in range(1, N_BANDS + 1):
        em = sim.get_eigenmode(
            f, mp.X, where, band_num, mp.Vector3(f * seed_neff, 0, 0),
            parity=mp.NO_PARITY)
        n = float(np.real(em.k.x)) / f
        neff.append(n)
        ey, ez = _sample(em, mp, mea_c, mea_w, mon_sz, *FINE_GRID)
        eys.append(ey)
        ezs.append(ez)
        if laddered:
            seed_neff = n     # next band seeds from THIS band's converged value
    return neff, eys, ezs


def _report(name, neff, eys, ezs):
    print(f"\n===== {name} =====", flush=True)
    print(f"  n_eff: {[round(v, 6) for v in neff]}", flush=True)

    ties = [(i + 1, j + 1) for i in range(len(neff)) for j in range(i + 1, len(neff))
            if abs(neff[i] - neff[j]) < 1e-9]
    print(f"  exact ties (|dn_eff| < 1e-9): {ties if ties else 'none'}", flush=True)

    fields = dict(Ey=np.stack(eys), Ez=np.stack(ezs))
    eta = field_overlap_matrix(fields, fields, dy=1.0, dz=1.0)
    np.fill_diagonal(eta, 0.0)
    i, j = np.unravel_index(np.argmax(eta), eta.shape)
    print(f"  worst off-diagonal overlap (fine {FINE_GRID[0]}x{FINE_GRID[1]} grid): "
          f"band {i + 1} vs band {j + 1} = {eta[i, j]:.6f}", flush=True)
    return ties, float(eta[i, j])


def main():
    import meep as mp

    mmi1, _ = _build()
    x_loc = -(mmi1.mm_length / 2 - BODY_STUB)
    n_clad = SOI_160.sio2_index(WL)
    print(f"cross-section: mmi1 local x = {x_loc:.4f} um (interface 0), "
          f"wl = {WL} um, res {RES}, {N_BANDS} bands, cladding {n_clad:.6f}",
          flush=True)

    solver = FDTDSolver(res=RES, sim_width=SIM_W, sim_height=SIM_H, dpml=DPML,
                        source_margin=0.0, monitor_margin=0.25, three_d=True,
                        material_stack=SOI_160, decay_tol=1e-11)

    results = {}
    for name, laddered in (("fixed seed (today's code)", False),
                            ("laddered seed (band k <- band k-1's neff)", True)):
        sim, where, mea_c, mea_w = solver._build_section_mode_sim(
            mmi1, WL, x_loc - HALF_SPAN, x_loc + HALF_SPAN,
            monitor_center=0.0, monitor_size=SIM_W)
        mon_sz = solver._mon_sz
        f = 1.0 / WL
        neff, eys, ezs = _solve_sequence(sim, mp, where, f, mea_c, mea_w, mon_sz, laddered)
        sim.reset_meep()
        results[name] = _report(name, neff, eys, ezs)

    print("\n===== verdict =====", flush=True)
    print("  fixed-seed ties/overlap are the reference: known duplicates at "
          "band_num 9/10, 13/14\n  (job 11575516, res 40). If the laddered "
          "run has FEWER ties and a lower worst\n  overlap, laddering helps. "
          "If it has ties the fixed run didn't, laddering made\n  it worse -- "
          "re-finding the previous band rather than avoiding it.", flush=True)


if __name__ == "__main__":
    main()
