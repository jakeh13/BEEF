"""Do either of the two follow-up k-seed ladder variants close item 12?
(validation_markdown/fdtd_vs_EME_mode_validation.md item 12.)

Where this picks up
--------------------
`diagnose_kpoint_ladder.py` established two things by seeding each band from a
FIXED guess vs. from the RAW previous band's converged n_eff (job on
2026-07-30, res 40, mux interface 0):

    fixed (today's code): 3 exact ties in 20 raw solves, 17/20 distinct.
    naive ladder:          1 exact tie,                   19/20 distinct.

Naive laddering recovers real basis rank (2 modes the fixed seed silently
drops), but also reproduces the predicted risk: seeding band k AT band
(k-1)'s own answer can re-find that same eigenvector. It moved the duplicate,
it didn't remove it.

That script only ran 20 RAW solves with no accept/skip/retry -- a fair test of
the seeding hypothesis in isolation, but not what production would actually
do, since `_iter_distinct_bands()` retries past a detected duplicate. This
script closes that gap: it runs a REAL guarded search (accept/skip/retry,
target N distinct modes, same overlap threshold and same fine 128x16 grid the
item-11 measurement showed trustworthy) under four seeding policies, so
"attempts needed" and "does laddering ever get stuck re-finding the same
duplicate" are measured directly rather than inferred.

The four policies
------------------
    A. fixed        seed = _KPOINT_SEED_NEFF always (today's code).
    B. naive ladder  seed = whatever band_num-1 converged to, RAW (accepted or
                     skipped) -- same policy diagnose_kpoint_ladder.py tested,
                     now under a real accept/skip/retry loop instead of a flat
                     20-solve sweep, to see whether a skip traps it.
    C. step ladder   seed = (last RAW converged n_eff) - STEP. The basis here
                     is monotonically descending in n_eff, so subtracting a
                     small step aims the next search PAST the band just found
                     instead of landing exactly on it -- a direct fix for B's
                     observed self-reconvergence failure.
    D. accepted ladder  seed = the last ACCEPTED band's n_eff, not the raw
                     previous solve. A skip does not move the seed, so a
                     rejected duplicate can't poison the next attempt's
                     starting point the way it can under B.

STEP (policy C) is set well below the tightest genuine gap seen in this
device's basis so far (~0.02, between the two near-cutoff bands at 1.433726
and 1.413238) -- large enough to plausibly escape exact reconvergence, small
enough that it should not leap over a real neighbouring mode. Documented as a
first guess, not tuned.

What "resolved" means here
---------------------------
Not zero duplicates ever -- `_iter_distinct_bands()` already tolerates and
retries past those (that's what the guard is FOR). What matters for closing
item 12:

    attempts needed to reach N_TARGET distinct modes (fewer is better --
        every skip is a wasted eigensolve);
    whether any policy gets stuck retrying the same duplicate repeatedly
        (the real risk B's raw test couldn't see, since it never retried);
    whether all four policies converge on the SAME final set of physical
        modes -- if the search order changes which eigenvalue gets "found
        first" near closely-spaced bands, that would be a new, more subtle
        problem than any of items 0/6/11/12 have raised so far.

Reshaped 2026-07-31 after job 11586046 HUNG
--------------------------------------------
The first submission never printed a single policy block. It was not slow --
it was stuck. `N_TARGET = 20` is unattainable at this cross-section, and asking
for it walks `sim.get_eigenmode()` off the end of the spectrum:

    bands  1-15   guided (n_eff 2.7064 down to 1.4848, cladding 1.446804)
    bands 16-22   BELOW the cladding index -- unguided, but still convergent
    band  23      no k reaches the target frequency at all

Three duplicate pairs sit inside that range (10/11, 14/15, 19/20), so 20
DISTINCT modes needs band_num 23 -- and MPB's k-search for band 23 does not
converge, it limit-cycles between k = 0.2346 (f = 0.4775) and k = 1.7480
(f = 1.0283) either side of the target f = 0.7634. It burned 161 MPB solves,
59% of the whole job's eigensolve work, and would have spun until the wall
clock killed it. Bands 1-22 needed 3-8 solves each, so the original ~25 s/solve
budget was roughly right; the wall time was never the problem.

Two changes follow, and the first is the load-bearing one:

    N_TARGET is now 16 -- what the device actually requests (`modal_
        multiplexer_hybrid.py`'s `N_MODES`) -- not 20 chosen for comparability
        with the v1 raw sweep. v1's bands 16-20 were unguided anyway, so the
        comparability that was lost was with numbers that had no physical
        standing.
    The search STOPS when a raw solve returns n_eff <= the cladding index.
        Past that point the mode is not guided, there is nothing above it to
        find, and the next band_num is where the k-search stops converging.
        This is a PHYSICAL bound, not a timeout: MPB runs in C, so a Python
        signal/alarm cannot interrupt a non-converging solve once it has been
        entered. The only reliable defence is to not ask for the band.

Note the same exposure exists in production -- `_iter_distinct_bands()` allows
band_num up to `n_modes + _DUP_BAND_MAX_EXTRA` (24 for this device) and its
`RuntimeError` can only fire AFTER a solve returns, so a non-converging solve
never reaches it. Recorded as its own item; not fixed here.
"""
from __future__ import annotations

import numpy as np

from beef import MMI, MMIPort, WaveguideRun, CompoundGeometry, FDTDSolver, SOI_160
from beef.mode_identity import field_overlap_matrix
from beef.solvers.fdtd import _DUP_BAND_OVERLAP_THRESHOLD, _KPOINT_SEED_NEFF

WL = 1.31
N_TARGET = 16          # what modal_multiplexer_hybrid.py actually requests
                       # (N_MODES). Was 20 for comparability with v1's raw
                       # sweep; that walked off the end of the spectrum and hung
                       # job 11586046 -- see the module docstring. Expect the
                       # cladding stop below to fire first (15 guided at res 40),
                       # which is the honest answer to "how many distinct modes
                       # does this cross-section have", not a truncation.
MAX_ATTEMPTS = 24      # target + 8, mirroring _DUP_BAND_MAX_EXTRA's own headroom
                       # convention. A backstop only -- the cladding stop is
                       # what is expected to end each search, and it is the one
                       # that keeps band_num away from the non-convergent region.
                       # Wall time: bands 1-22 of job 11586046 took ~58 min for
                       # 22 attempts (~2.6 min each, deeper bands costing more),
                       # so a ~18-attempt search is ~45 min and four policies
                       # ~3 h.
RES = 40.0              # production resolution, matches the confirmed reproduction
BODY_STUB = 0.75
SIM_W, SIM_H = 5.0, 1.5
DPML = 0.5
HALF_SPAN = 0.05
FINE_GRID = (128, 16)   # item 11 showed this agrees with the shipped 16x8 grid
STEP = 0.005            # policy C's nudge -- see module docstring
TIE_TOL = 1e-9

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


# Each policy: (accepted_neff, last_raw_neff) -> seed n_eff for the NEXT attempt.
POLICIES = {
    "A. fixed (today's code)":
        lambda accepted, last_raw: _KPOINT_SEED_NEFF,
    "B. naive ladder (raw previous solve)":
        lambda accepted, last_raw: last_raw if last_raw is not None else _KPOINT_SEED_NEFF,
    "C. step ladder (raw previous - STEP)":
        lambda accepted, last_raw: (last_raw - STEP) if last_raw is not None else _KPOINT_SEED_NEFF,
    "D. accepted ladder (last ACCEPTED band)":
        lambda accepted, last_raw: accepted[-1] if accepted else _KPOINT_SEED_NEFF,
}


def _guarded_search(sim, mp, where, f, mea_c, mea_w, mon_sz, policy, n_clad):
    """Real accept/skip/retry search for up to N_TARGET distinct modes,
    mirroring `_iter_distinct_bands()`'s own logic, parameterized only by seed
    policy.

    Terminates on whichever comes first:

        "target"    N_TARGET distinct modes accepted;
        "unguided"  a raw solve came back at or below the cladding index, so
                    the guided spectrum is exhausted and the next band_num is
                    where MPB's k-search stops converging (job 11586046);
        "attempts"  MAX_ATTEMPTS backstop.

    Returns the stop reason alongside the results, because "found 15 of 16
    because there are only 15" and "found 15 of 16 because it ran out of
    attempts" are different findings and must not read the same.
    """
    accepted_neff, accepted_ey, accepted_ez = [], [], []
    skipped = []
    last_raw_neff = None
    band_num = 0
    stop = "attempts"
    unguided = None
    while True:
        if len(accepted_neff) >= N_TARGET:
            stop = "target"
            break
        if band_num >= MAX_ATTEMPTS:
            stop = "attempts"
            break
        band_num += 1
        seed_neff = policy(accepted_neff, last_raw_neff)
        em = sim.get_eigenmode(
            f, mp.X, where, band_num, mp.Vector3(f * seed_neff, 0, 0),
            parity=mp.NO_PARITY)
        n = float(np.real(em.k.x)) / f
        last_raw_neff = n

        # Guided-spectrum stop, BEFORE the overlap sampling: an unguided solve
        # is not a basis candidate, and asking for the band above it is what
        # hung job 11586046. Deliberately checked on the RAW n_eff -- whether
        # this band would have been accepted or skipped as a duplicate is
        # irrelevant once it is below cutoff.
        if n <= n_clad:
            unguided = (band_num, round(n, 6))
            stop = "unguided"
            break

        ey, ez = _sample(em, mp, mea_c, mea_w, mon_sz, *FINE_GRID)

        if accepted_ey:
            existing = dict(Ey=np.stack(accepted_ey), Ez=np.stack(accepted_ez))
            candidate = dict(Ey=ey[None], Ez=ez[None])
            eta = field_overlap_matrix(existing, candidate, dy=1.0, dz=1.0)
            worst = float(np.max(eta))
            if worst > _DUP_BAND_OVERLAP_THRESHOLD:
                skipped.append((band_num, round(worst, 6), round(n, 6)))
                continue

        accepted_neff.append(n)
        accepted_ey.append(ey)
        accepted_ez.append(ez)

    return accepted_neff, skipped, band_num, stop, unguided


def main():
    import meep as mp

    mmi1, _ = _build()
    x_loc = -(mmi1.mm_length / 2 - BODY_STUB)
    n_clad = SOI_160.sio2_index(WL)
    print(f"cross-section: mmi1 local x = {x_loc:.4f} um (interface 0), "
          f"wl = {WL} um, res {RES}, target {N_TARGET} distinct modes, "
          f"max {MAX_ATTEMPTS} attempts, cladding {n_clad:.6f}", flush=True)
    print(f"policy C step = {STEP}", flush=True)
    print(f"stopping when a raw solve returns n_eff <= {n_clad:.6f} "
          f"(guided spectrum exhausted -- see module docstring, job 11586046)",
          flush=True)

    solver = FDTDSolver(res=RES, sim_width=SIM_W, sim_height=SIM_H, dpml=DPML,
                        source_margin=0.0, monitor_margin=0.25, three_d=True,
                        material_stack=SOI_160, decay_tol=1e-11)

    results = {}
    for name, policy in POLICIES.items():
        sim, where, mea_c, mea_w = solver._build_section_mode_sim(
            mmi1, WL, x_loc - HALF_SPAN, x_loc + HALF_SPAN,
            monitor_center=0.0, monitor_size=SIM_W)
        mon_sz = solver._mon_sz
        f = 1.0 / WL
        neff, skipped, attempts, stop, unguided = _guarded_search(
            sim, mp, where, f, mea_c, mea_w, mon_sz, policy, n_clad)
        sim.reset_meep()

        print(f"\n===== {name} =====", flush=True)
        print(f"  attempts used: {attempts} (skips: {len(skipped)})", flush=True)
        print(f"  skipped (band_num, overlap, raw n_eff): "
              f"{skipped if skipped else 'none'}", flush=True)
        print(f"  distinct guided modes found: {len(neff)} of {N_TARGET} requested",
              flush=True)
        print(f"  stopped because: {stop}", flush=True)
        if stop == "unguided":
            print(f"    band_num {unguided[0]} came back at n_eff "
                  f"{unguided[1]} <= cladding {n_clad:.6f} -- guided spectrum "
                  f"exhausted, not a truncation", flush=True)
        elif stop == "attempts":
            print(f"    hit the MAX_ATTEMPTS={MAX_ATTEMPTS} backstop with the "
                  f"guided spectrum NOT exhausted -- this policy wasted more "
                  f"attempts on duplicates than the headroom allows", flush=True)
        print(f"  n_eff: {[round(v, 6) for v in neff]}", flush=True)
        results[name] = dict(neff=neff, attempts=attempts, skips=len(skipped),
                             stop=stop)

    print("\n===== cross-policy comparison =====", flush=True)
    ref_name = "A. fixed (today's code)"
    ref_set = np.array(sorted(results[ref_name]["neff"], reverse=True))
    print(f"  {'policy':38s} {'found':>6s} {'attempts':>9s} {'skips':>6s} "
          f"{'stop':>9s}  max |dn_eff| vs A", flush=True)
    for name, r in results.items():
        this_set = np.array(sorted(r["neff"], reverse=True))
        m = min(len(this_set), len(ref_set))
        # Compare over the common prefix: the counts are themselves a result
        # now (a policy that recovers more distinct guided modes is the whole
        # point), so a length mismatch must not collapse the comparison to nan.
        if m:
            max_dev = float(np.max(np.abs(this_set[:m] - ref_set[:m])))
            agrees = max_dev < 1e-4 and len(this_set) == len(ref_set)
        else:
            max_dev, agrees = float("nan"), False
        note = "same set" if agrees else (
            f"same first {m}" if max_dev < 1e-4 else "DIFFERENT modes")
        print(f"  {name:38s} {len(this_set):6d} {r['attempts']:9d} "
              f"{r['skips']:6d} {r['stop']:>9s}  {max_dev:.6f}  ({note})",
              flush=True)

    print("\n===== how to read this =====", flush=True)
    print("  MORE distinct guided modes found  ->  that policy recovers basis", flush=True)
    print("  rank the fixed seed silently loses; this is the primary number,", flush=True)
    print("  since a missing mode is a dropped HybridSolver channel.", flush=True)
    print("  Same count, FEWER attempts  ->  strict improvement: cheaper and", flush=True)
    print("  equally correct. Every skip is a wasted eigensolve.", flush=True)
    print("  'DIFFERENT modes' means seed order changed WHICH eigenvalue near a", flush=True)
    print("  crowded region gets found first -- a new finding, not covered by", flush=True)
    print("  diagnose_kpoint_ladder.py's raw-solve comparison.", flush=True)
    print("  stop='attempts' on any policy is a red flag: it means duplicates", flush=True)
    print("  ate the whole +8 headroom that production's guard also relies on.", flush=True)


if __name__ == "__main__":
    main()
