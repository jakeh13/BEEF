"""Does `_iter_distinct_bands()`'s duplicate check survive its own grid?
(validation_markdown/fdtd_vs_EME_mode_validation.md item 11.)

The question
------------
`FDTDSolver._iter_distinct_bands()` decides "same eigenvector or not" by
sampling each candidate band's transverse E on a FIXED `_DUP_BAND_DIAG_NY` x
`_DUP_BAND_DIAG_NZ` = 16 x 8 grid spanning the whole MONITOR, then thresholding
the normalized overlap at 0.9. The constant's justifying comment -- "only needs
to distinguish same eigenvector (0.9999...) from different (<= 0.1257)" -- was
measured on the SIX-mode directional coupler.

At the modal multiplexer's interface it is asked to judge SIXTEEN modes on a
2.8 um core inside a 5.0 um monitor: 0.3125 um transverse pitch, only ~8 of the
16 points inside the core at all, against a ~0.19 um lobe pitch for the top of
the basis. Idealized slab profiles sin(m*pi*y/W) sampled on that exact grid
alias mode m onto mode 18-m at eta ~ 0.9988 -- past the 0.9 threshold, i.e. a
GENUINE distinct mode would be discarded as a duplicate and silently skipped.

That is an argument about idealized profiles. This script settles it on the
real fields.

What it does
------------
Solves band_num 1..N_BANDS at the mux's interface 0 with the same
`get_eigenmode()` call and the same k seed the guard itself uses, samples each
band on SEVERAL diagnostic grids (the shipped 16x8 first, then progressively
finer), and for each grid replays the guard's own greedy accept/skip logic.

Three things to read off:

    which band_nums each grid SKIPS -- if the shipped grid skips a band the
        fine grids keep, that band is a real mode being thrown away, and item
        11 is a live bug rather than a latent one.
    the resulting n_eff sequence -- compare against meow's spectrum from job
        11575516 (recorded under item 1), which finds one guided mode meep
        does not.
    the MARGIN: worst off-diagonal overlap among ACCEPTED bands. The guard's
        0.9 threshold is only meaningful if genuinely-distinct modes sit well
        below it. On the coupler that number was 0.1257. If it approaches 0.9
        here, the threshold has no discriminating power left at this mode
        count regardless of which bands happen to be flagged.

Cost
----
Deliberately res 20 -- the resolution where the shipped guard was observed
firing twice (bands 11 and 16, job 11575516). Eigensolves are the cheap part
(`get_eigenmode()` does not time-step); the sampling dominates, so the grid
list is kept short. Minutes, not hours: a login-node smoke run, not a batch job.
"""
from __future__ import annotations

import numpy as np

from beef import (
    MMI, MMIPort, WaveguideRun, CompoundGeometry, FDTDSolver, SOI_160,
)
from beef.mode_identity import field_overlap_matrix
from beef.solvers.fdtd import _DUP_BAND_OVERLAP_THRESHOLD, _KPOINT_SEED_NEFF

WL = 1.31
N_MODES = 16          # what the device requests
N_BANDS = 20          # solve past it, so a skip has somewhere to go
RES = 20.0
BODY_STUB = 0.75
SIM_W, SIM_H = 5.0, 1.5
DPML = 0.5
HALF_SPAN = 0.05

# (ny, nz). First entry is what beef/solvers/fdtd.py ships today.
GRIDS = [(16, 8), (32, 8), (64, 16), (128, 16)]

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
    """Same midpoint-of-edges construction _iter_distinct_bands() uses."""
    y_edges = np.linspace(mea_c - mea_w / 2, mea_c + mea_w / 2, ny + 1)
    y = 0.5 * (y_edges[:-1] + y_edges[1:])
    z_edges = np.linspace(-mon_sz / 2, mon_sz / 2, nz + 1)
    z = 0.5 * (z_edges[:-1] + z_edges[1:])
    ey = np.array([[em.amplitude(mp.Vector3(0, float(yy), float(zz)), mp.Ey)
                    for zz in z] for yy in y])
    ez = np.array([[em.amplitude(mp.Vector3(0, float(yy), float(zz)), mp.Ez)
                    for zz in z] for yy in y])
    return ey, ez


def _replay_guard(eys, ezs, n_modes):
    """The accept/skip logic of _iter_distinct_bands(), on already-sampled
    fields. Returns (accepted band_nums, skipped band_nums, worst off-diagonal
    overlap among accepted)."""
    accepted, skipped = [], []
    acc_ey, acc_ez = [], []
    worst_accepted = 0.0
    for band_num in range(1, len(eys) + 1):
        if len(accepted) == n_modes:
            break
        ey, ez = eys[band_num - 1], ezs[band_num - 1]
        if accepted:
            existing = dict(Ey=np.stack(acc_ey), Ez=np.stack(acc_ez))
            candidate = dict(Ey=ey[None], Ez=ez[None])
            eta = field_overlap_matrix(existing, candidate, dy=1.0, dz=1.0)
            worst = float(np.max(eta))
            if worst > _DUP_BAND_OVERLAP_THRESHOLD:
                skipped.append((band_num, worst))
                continue
            worst_accepted = max(worst_accepted, worst)
        accepted.append(band_num)
        acc_ey.append(ey)
        acc_ez.append(ez)
    return accepted, skipped, worst_accepted


def main():
    import meep as mp

    mmi1, _ = _build()
    x_loc = -(mmi1.mm_length / 2 - BODY_STUB)
    n_clad = SOI_160.sio2_index(WL)
    print(f"cross-section: mmi1 local x = {x_loc:.4f} um (interface 0), "
          f"wl = {WL} um, res {RES}, solving {N_BANDS} bands, "
          f"guard would keep {N_MODES}")
    print(f"cladding index = {n_clad:.6f}", flush=True)

    solver = FDTDSolver(res=RES, sim_width=SIM_W, sim_height=SIM_H, dpml=DPML,
                        source_margin=0.0, monitor_margin=0.25, three_d=True,
                        material_stack=SOI_160, decay_tol=1e-11)
    sim, where, mea_c, mea_w = solver._build_section_mode_sim(
        mmi1, WL, x_loc - HALF_SPAN, x_loc + HALF_SPAN,
        monitor_center=0.0, monitor_size=SIM_W)
    mon_sz = solver._mon_sz
    f = 1.0 / WL

    print(f"\nmonitor: centre {mea_c:.4f}, width {mea_w:.4f} um; "
          f"z span {mon_sz:.4f} um; core is {MM_WIDTH} um wide", flush=True)

    # ---- one eigensolve per band, reused by every grid --------------------
    print(f"\n===== solving {N_BANDS} bands (unguarded, no skipping) =====",
          flush=True)
    ems, neff_raw = [], []
    for band_num in range(1, N_BANDS + 1):
        em = sim.get_eigenmode(
            f, mp.X, where, band_num, mp.Vector3(f * _KPOINT_SEED_NEFF, 0, 0),
            parity=mp.NO_PARITY)
        ems.append(em)
        neff_raw.append(float(np.real(em.k.x)) / f)
        print(f"  band {band_num:2d}: n_eff {neff_raw[-1]:.6f}"
              f"{'' if neff_raw[-1] > n_clad else '   (unguided)'}", flush=True)

    exact_ties = [(i + 1, j + 1) for i in range(N_BANDS) for j in range(i + 1, N_BANDS)
                  if abs(neff_raw[i] - neff_raw[j]) < 1e-9]
    print(f"\n  exact n_eff ties (|dn| < 1e-9): "
          f"{exact_ties if exact_ties else 'none'}", flush=True)
    print("  ^ an exact tie is INDEPENDENT evidence of a real duplicate solve, "
          "and needs\n    no field sampling at all -- the reference the grids "
          "below are judged against.", flush=True)

    # ---- replay the guard on each grid ------------------------------------
    for ny, nz in GRIDS:
        shipped = " (SHIPPED)" if (ny, nz) == GRIDS[0] else ""
        pitch_y = mea_w / ny
        in_core = int(np.sum(np.abs(
            0.5 * (np.linspace(mea_c - mea_w / 2, mea_c + mea_w / 2, ny + 1)[:-1]
                   + np.linspace(mea_c - mea_w / 2, mea_c + mea_w / 2, ny + 1)[1:])
        ) <= MM_WIDTH / 2))
        print(f"\n===== diagnostic grid {ny} x {nz}{shipped} =====", flush=True)
        print(f"  transverse pitch {pitch_y:.4f} um, {in_core}/{ny} points "
              f"inside the {MM_WIDTH} um core", flush=True)

        eys, ezs = zip(*[_sample(em, mp, mea_c, mea_w, mon_sz, ny, nz)
                         for em in ems])
        accepted, skipped, worst_acc = _replay_guard(list(eys), list(ezs), N_MODES)

        print(f"  skipped band_nums: "
              f"{[(b, round(w, 6)) for b, w in skipped] if skipped else 'none'}",
              flush=True)
        print(f"  accepted band_nums: {accepted}", flush=True)
        print(f"  worst off-diagonal overlap among ACCEPTED: {worst_acc:.4f}"
              f"   (coupler reference: 0.1257; threshold "
              f"{_DUP_BAND_OVERLAP_THRESHOLD})", flush=True)
        kept_neff = [round(neff_raw[b - 1], 5) for b in accepted]
        print(f"  resulting n_eff: {kept_neff}", flush=True)
        print(f"  guided: {sum(1 for v in kept_neff if v > n_clad)}/{len(kept_neff)}",
              flush=True)

    print("\n===== how to read this =====", flush=True)
    print("  Shipped grid skips a band the fine grids KEEP, and that band is "
          "not in the\n  exact-tie list  ->  a genuine mode is being discarded. "
          "Item 11 is live.", flush=True)
    print("  All grids skip the same bands, and those bands ARE exact ties  ->  "
          "the\n  shipped grid gets the right answer here; the margin line is "
          "then what says\n  whether it does so with any room to spare.",
          flush=True)

    sim.reset_meep()


if __name__ == "__main__":
    main()
