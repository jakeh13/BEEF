"""TEMPORARY -- one-off diagnostic, delete once
mode_crossing_correspondence_repair_validation.md's Stage 0 item (a) is
closed: "measure the duplicate-threshold's own self-overlap values
directly across domain size, to see exactly which pair of modes is
fighting the 0.90 threshold and whether that overlap value itself trends
with domain size or is genuinely noisy."

`FDTDSolver._iter_distinct_bands()` only ever LOGS a `worst` overlap value
when it crosses `_DUP_BAND_OVERLAP_THRESHOLD` (0.90, an outright skip). It
throws away every near-miss -- an attempt accepted at worst=0.85 leaves no
trace distinguishable from one accepted at worst=0.05. That is exactly the
number the 2026-08-10 walked-back conclusion needs: the domain sweep in
`diag_domain_convergence_neff.py` found meep's accepted 14-mode basis
reshuffling non-monotonically between sim_width 5.0-9.0 (a ~2.270 mode
drops out at 6.0, reappears at 9.0) and outright failing at 11.0 (6
duplicate-band hits). If that reshuffling is `_iter_distinct_bands()`'s own
0.90 threshold sitting on a knife-edge for a genuinely close-packed
cluster, the worst-overlap value at each attempt should sit close to 0.90
and wander across domain size. If it is something else, the value should
sit nowhere near 0.90.

This script re-implements `_iter_distinct_bands()`'s exact loop (same
k-seed laddering since 2026-07-31, same n_clad stop condition, same
shipped 16x8 diagnostic grid) verbatim -- not calling the guarded method
itself, because the method only yields accepted modes and only logs a
skip's worst value, not every attempt's full overlap row against every
already-accepted mode. Re-implementing it here is the same pattern
`diagnose_dup_band_grid.py` already used to inspect this guard's internals
(see its `_replay_guard()`), extended to run live eigensolves per domain
size (fields are NOT reusable across sim_width, unlike that script's
single-domain grid comparison) and to keep every attempt's full eta row,
not just the max.

Domain sweep matches `diag_domain_convergence_neff.py` exactly (sim_width
5.0-11.0 um, 1.0 um steps, fixed aspect ratio, same epsilon fix for the
float-truncation mesh bug), shipped N_MODES=14, wl=1.26 um (the originally
failing point), interface 0 (FDTD-in <-> EME-MMI1) -- reusing the real
production script's own construction via the same truncated-exec trick, so
this measures the actual device, not a standalone approximation of it.

Cost: `get_eigenmode()` does not time-step (seconds each); the 16x8 field
sample per attempt is the same cheap grid the guard itself uses in
production, not the fine grids `diagnose_dup_band_grid.py` explored. 7
domain points x ~14-25 attempts each (more than 14 whenever a skip fires)
x 1 eigensolve+sample each -- same order of cost as
`diag_domain_convergence_neff.py`'s meep-only half, comfortably a batch
job, not a login-node run.
"""
import logging
import os
import sys

import numpy as np

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "..", "modal_multiplexer",
                       "python", "modal_multiplexer_hybrid.py")
SCRIPT = os.path.abspath(SCRIPT)
MARKER = "hybrid = HybridSolver("

with open(SCRIPT) as f:
    src = f.read()

end_marker = "\n)\n"
start = src.index(MARKER)
end = src.index(end_marker, start) + len(end_marker)
truncated = src[:end]

print(f"truncated source: {len(truncated)} of {len(src)} chars", flush=True)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s", stream=sys.stdout)

ns = {"__name__": "__main__", "__file__": SCRIPT}
code = compile(truncated, SCRIPT, "exec")
exec(code, ns)

hybrid = ns["hybrid"]
mux = ns["mux"]
fdtd_solver = ns["fdtd_solver"]

from beef.mode_identity import field_overlap_matrix
from beef.solvers.fdtd import _DUP_BAND_OVERLAP_THRESHOLD, _DUP_BAND_DIAG_NY, _DUP_BAND_DIAG_NZ, \
    _DUP_BAND_MAX_EXTRA, _KPOINT_SEED_NEFF

N_MODES = 14  # matches diag_domain_convergence_neff.py's shipped-target value.
WL = 1.26
X_IFACE = hybrid.boundaries[1]  # interface 0: FDTD-in <-> EME-MMI1
HALF_SPAN = 0.05  # matches HybridSolver's own _CORR_HALF_SPAN
X0, X1 = X_IFACE - HALF_SPAN, X_IFACE + HALF_SPAN

print(f"interface 0 at x={X_IFACE:.4f}, query window [{X0:.4f}, {X1:.4f}], "
      f"N_MODES={N_MODES}, wl={WL}, threshold={_DUP_BAND_OVERLAP_THRESHOLD}", flush=True)

SIM_WIDTHS = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
ASPECT = 1.5 / 5.0
EPS = 1e-6


def safe_height(width: float) -> float:
    # Same float-truncation fix diag_domain_convergence_neff.py needed --
    # sim_width * ASPECT is not always exactly representable in binary
    # float, so meow's int()-truncating mesh construction can land one
    # pixel short of intent. Meep's grid here does not hit that specific
    # bug (no meow call in this script), but applying the identical nudge
    # keeps sim_height byte-identical to the already-computed
    # diag_domain_convergence_neff.py points, for a fair overlay.
    return round(width * ASPECT, 6) + EPS


def replay_with_full_overlap(sim, mp, where, f, mea_c, mea_w, mon_sz, n_modes, n_clad):
    """Verbatim re-implementation of FDTDSolver._iter_distinct_bands(), except
    every attempt's FULL eta row (overlap against every already-accepted
    mode, not just the max) is recorded, whether or not it crosses the
    real skip threshold."""
    diag_y_edges = np.linspace(mea_c - mea_w / 2, mea_c + mea_w / 2, _DUP_BAND_DIAG_NY + 1)
    diag_y = 0.5 * (diag_y_edges[:-1] + diag_y_edges[1:])
    diag_z_edges = np.linspace(-mon_sz / 2, mon_sz / 2, _DUP_BAND_DIAG_NZ + 1)
    diag_z = 0.5 * (diag_z_edges[:-1] + diag_z_edges[1:])

    accepted_ey, accepted_ez, accepted_band, accepted_raw_neff = [], [], [], []
    attempts = []  # dicts: band_num, raw_neff, eta_row (vs each accepted idx), worst, worst_vs, accepted(bool)
    band_num = 0
    n_accepted = 0
    last_raw_neff = None
    max_band_num = n_modes + _DUP_BAND_MAX_EXTRA
    while n_accepted < n_modes:
        band_num += 1
        if band_num > max_band_num:
            print(f"    !! max_band_num ({max_band_num}) exceeded, matching "
                  f"production's own RuntimeError -- stopping this domain point here", flush=True)
            break
        seed_neff = last_raw_neff if last_raw_neff is not None else _KPOINT_SEED_NEFF
        em = sim.get_eigenmode(
            f, mp.X, where, band_num, mp.Vector3(f * seed_neff, 0, 0),
            parity=mp.NO_PARITY,
        )
        n = float(np.real(em.k.x)) / f
        last_raw_neff = n
        if n <= n_clad:
            print(f"    band_num={band_num}: n_eff={n:.6f} <= n_clad={n_clad:.6f} "
                  f"-- guided edge reached, stopping (matches production's stop rule)", flush=True)
            break
        ey = np.array([[em.amplitude(mp.Vector3(0, float(y), float(z)), mp.Ey)
                        for z in diag_z] for y in diag_y])
        ez = np.array([[em.amplitude(mp.Vector3(0, float(y), float(z)), mp.Ez)
                        for z in diag_z] for y in diag_y])
        rec = dict(band_num=band_num, raw_neff=n, eta_row=None, worst=0.0, worst_vs=-1, accepted=True)
        if n_accepted > 0:
            existing = dict(Ey=np.stack(accepted_ey), Ez=np.stack(accepted_ez))
            candidate = dict(Ey=ey[None], Ez=ez[None])
            eta = field_overlap_matrix(existing, candidate, dy=1.0, dz=1.0)
            eta_row = eta[:, 0]  # overlap of this candidate against each already-accepted mode
            worst_idx = int(np.argmax(eta_row))
            worst = float(eta_row[worst_idx])
            rec.update(eta_row=eta_row.copy(), worst=worst,
                       worst_vs=accepted_band[worst_idx])
            if worst > _DUP_BAND_OVERLAP_THRESHOLD:
                rec["accepted"] = False
                attempts.append(rec)
                print(f"    band_num={band_num}: n_eff={n:.6f}  SKIPPED (dup of "
                      f"band {accepted_band[worst_idx]}, overlap {worst:.6f} > "
                      f"{_DUP_BAND_OVERLAP_THRESHOLD})", flush=True)
                continue
        attempts.append(rec)
        tag = f"worst={rec['worst']:.6f} vs band {rec['worst_vs']}" if rec['worst_vs'] >= 0 else "first mode"
        print(f"    band_num={band_num}: n_eff={n:.6f}  accepted (logical mode "
              f"{n_accepted}, {tag})", flush=True)
        accepted_ey.append(ey)
        accepted_ez.append(ez)
        accepted_band.append(band_num)
        accepted_raw_neff.append(n)
        n_accepted += 1
    return attempts, accepted_band, accepted_raw_neff


all_attempts = {}  # sim_width -> attempts list
all_accepted = {}  # sim_width -> (accepted_band, accepted_raw_neff)

for sim_width in SIM_WIDTHS:
    sim_height = safe_height(sim_width)
    print(f"\n{'='*20} sim_width={sim_width}, sim_height={sim_height:.6f} {'='*20}", flush=True)

    fdtd_solver.sim_width = sim_width
    fdtd_solver.sim_height = sim_height

    import meep as mp
    n_clad = fdtd_solver.material_stack.sio2_index(WL)
    sim, where, mea_c, mea_w = fdtd_solver._build_section_mode_sim(mux, WL, X0, X1, None, None)
    mon_sz = fdtd_solver._mon_sz
    f = 1.0 / WL

    attempts, acc_band, acc_neff = replay_with_full_overlap(
        sim, mp, where, f, mea_c, mea_w, mon_sz, N_MODES, n_clad)
    all_attempts[sim_width] = attempts
    all_accepted[sim_width] = (acc_band, acc_neff)

    sim.reset_meep()

    n_skipped = sum(1 for a in attempts if not a["accepted"])
    worst_among_accepted = max((a["worst"] for a in attempts if a["accepted"] and a["worst_vs"] >= 0),
                                default=0.0)
    print(f"  -> {len(acc_band)} accepted, {n_skipped} skipped, worst overlap among "
          f"ACCEPTED attempts: {worst_among_accepted:.6f}", flush=True)

print(f"\n{'='*20} SUMMARY: worst-overlap value at every attempt, by domain size {'='*20}", flush=True)
print(f"threshold = {_DUP_BAND_OVERLAP_THRESHOLD}\n", flush=True)
header = f"{'sim_width':>10s}  " + "  ".join(f"band{i:>2d}" for i in range(1, 26))
print(header, flush=True)
for sim_width in SIM_WIDTHS:
    attempts = all_attempts[sim_width]
    row = []
    for a in attempts:
        mark = "*" if not a["accepted"] else (" " if a["worst_vs"] >= 0 else "-")
        row.append(f"{a['worst']:.3f}{mark}" if a["worst_vs"] >= 0 else "  -- ")
    print(f"{sim_width:>10.1f}  " + "  ".join(f"{v:>6s}" for v in row), flush=True)
print("\n(* = SKIPPED as duplicate, blank = accepted despite a nonzero worst overlap, "
      "-- = first mode / no prior accepted mode to compare against)", flush=True)

print(f"\n{'='*20} which accepted-mode index each near-miss/skip fought {'='*20}", flush=True)
for sim_width in SIM_WIDTHS:
    attempts = all_attempts[sim_width]
    close_calls = [a for a in attempts if a["worst_vs"] >= 0 and a["worst"] > 0.5]
    print(f"\n-- sim_width={sim_width} -- close calls (worst overlap > 0.5):", flush=True)
    if not close_calls:
        print("   none", flush=True)
    for a in close_calls:
        status = "SKIPPED" if not a["accepted"] else "accepted"
        print(f"   band_num={a['band_num']:2d} (n_eff={a['raw_neff']:.6f}) vs "
              f"already-accepted band_num={a['worst_vs']:2d}: overlap={a['worst']:.6f}  [{status}]",
              flush=True)

print(f"\n{'='*20} accepted n_eff sequence by domain size (for cross-reference against "
      f"diag_domain_convergence_neff.py) {'='*20}", flush=True)
for sim_width in SIM_WIDTHS:
    acc_band, acc_neff = all_accepted[sim_width]
    print(f"  sim_width={sim_width}: {np.array2string(np.array(acc_neff), precision=6)}", flush=True)

print("\nDIAGNOSTIC DONE", flush=True)
