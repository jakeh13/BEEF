"""Batch Bayesian optimization driver for the (mm_width, y_gap) search,
replacing the 106-point grid (mmi_4x4_90hybrid_ygap_grid.py) with a
sequential-round, constant-liar-batch BO loop over the SAME per-point
objective (mmi_4x4_90hybrid_wg_eval.py's hybrid-path evaluator).

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-28.

WHY BO INSTEAD OF THE GRID. The grid's real per-point cost (job 12383245,
measured) is 4:43:45, not the ~45-75 min extrapolated earlier -- 106 points
at 8-way concurrency is 45+ wall-clock hours even with nodes free on demand,
and this session has seen 4-node requests queue 3.5-22+ hours. BO trades full
surface coverage for far fewer evaluations to find the optimum (typically a
few dozen vs 106 on a smooth 2D response), at the same per-point cost.

WHY mm_width IS CATEGORICAL, NOT CONTINUOUS. wg_eval's own cost breakdown:
body meep solve ~1-2h, FDTD input sections ~3-3.7h. The body depends ONLY on
mm_width (verified by fingerprint, see wg_eval.py) and Phase A already paid
for it at exactly 5 widths (5.0, 5.5, 6.0, 6.5, 7.0 -- job 12387804 + the
12383245 timing job backfilling 6.0). Restricting mm_width to that lattice
means every BO candidate's body solve is a cache hit; treating mm_width as
continuous would pay a fresh ~1-2h body solve at essentially every candidate.
y_gap stays continuous (its own axis costs nothing extra to search finely --
it only changes which FDTD input sections get solved, which BO would be
paying for at grid resolution anyway).

SEARCH SPACE, REPARAMETERIZED 2026-09-01. dimensions =
[Categorical(MM_WIDTHS), Real(TAPER_MIN, TAPER_MAX) taper_width], with
y_gap = mm_width/4 - taper_width, clipped into each width's own valid range.
It was [Categorical(MM_WIDTHS), Real(0, 1) y_gap frac] until 24 points showed
every width's optimum sitting at taper_width 0.85-1.12 -- a diagonal ridge in
(width, y_gap) space, axis-aligned in (width, taper_width) space. The change
is exact and free: y_gap and taper_width determine each other given the width,
so all existing results convert with no new simulation. EPS keeps rounding to
2 decimals (wg_eval's own npz-naming precision) from landing exactly on the
inverse-taper boundary.

CANDIDATE DISTINCTNESS. MIN_YGAP_SEP rejects a candidate within 0.03 um of an
evaluated point at the same width, ASK_OVERSAMPLE refills the round from the
GP's own ranking, and explore_fill() falls back to farthest-point exploration.
Without these, round 6 (2026-09-01) proposed W=6.00 at y_gap 0.50/0.51/0.52
against points already measured at 0.50/0.51 and shrank to a single job.

OBJECTIVE. wg_eval.py already returns a single scalar, lower-is-better,
0 = both specs met (see its own objective()) and already internally optimizes
body length per point (the 81-point DELTA_L_GRID) -- BO here only searches
(mm_width, y_gap), length is not a BO dimension.

ROUND STRUCTURE, synchronous by design. Given the measured node-contention
severity this session (single 4-node requests queueing hours), this driver
does NOT try to keep N candidates in flight continuously. Each invocation:
  1. Absorbs any eval_W*_y*.npz on disk not yet in the state file (free
     training data -- the nominal point from job 12383245 seeds round 1).
  2. Harvests any PENDING candidates whose SLURM job has finished.
  3. If anything is still RUNNING/PENDING in the queue, reports status and
     stops -- no piling on more jobs while a round is outstanding.
  4. Otherwise fits skopt's GP on all results so far and proposes the next
     BATCH_SIZE candidates via Optimizer.ask(n_points=..., strategy="cl_min")
     (constant-liar batch acquisition -- avoids the batch collapsing onto
     one point, no BoTorch dependency needed), submits one sbatch job per
     candidate, and exits. Re-invoke later (or under the loop skill) to
     harvest + propose the next round.

State lives in wg_ygap_opt/bo_state.json (results) and bo_pending.json (in-
flight jobs) -- both plain JSON so a session crash mid-round loses nothing;
re-running just re-checks squeue/sacct.

Usage:
    python mmi_4x4_90hybrid_bo.py status              # harvest + report, no submit
    python mmi_4x4_90hybrid_bo.py propose [--batch N] # harvest, then submit next round
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

from mmi_4x4_90hybrid_ygap_grid import (
    MM_WIDTHS, WAVEGUIDE_WIDTH, Y_GAP_MIN, LENGTH_MIN, LENGTH_MAX, y_gaps_for)

EPS = 0.01
BATCH_SIZE_DEFAULT = 4
# Minimum y_gap separation, at the SAME mm_width, between a proposed candidate
# and any point already evaluated or already picked this round. Measured reason,
# not a guess: rounds 5-6 proposed W=6.00 at y_gap 0.50/0.51/0.52 while
# W=6.00 y_gap=0.50 was already in hand at objective 1.0047 -- three of four
# 4-node/4 h jobs would have re-measured a known point to 0.01 um, inside this
# device's own objective noise, and the round shrank to a single job. 0.03 um is
# 3x the npz naming precision (2 decimals) and ~1/35 of the widest y_gap range.
MIN_YGAP_SEP = 0.03
# opt.ask() is oversampled so MIN_YGAP_SEP rejections can be replaced from the
# GP's own ranking rather than shrinking the round.
ASK_OVERSAMPLE = 4
NOMINAL = (6.0, 0.70)

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
OUT_DIR = os.path.join(_ROOT, "example_output", "mmi_4x4_90hybrid", "wg_ygap_opt")
BATCH_SCRIPT = os.path.join(_ROOT, "batch", "verification",
                             "run_mmi_4x4_90hybrid_bo_point.batch")
RESULTS_PATH = os.path.join(OUT_DIR, "bo_state.json")
PENDING_PATH = os.path.join(OUT_DIR, "bo_pending.json")
os.makedirs(OUT_DIR, exist_ok=True)


def y_gap_bounds(mm_width):
    hi = mm_width / 4.0 - WAVEGUIDE_WIDTH - EPS
    return Y_GAP_MIN, hi


# The GP searches (mm_width, taper_width), NOT (mm_width, y_gap) -- measured
# reason, 2026-09-01: every width's optimum in the 24-point table sits at
# taper_width 0.85-1.12, i.e. the good region is a line of constant TAPER
# width, which is a diagonal ridge in (width, y_gap) space and axis-aligned in
# (width, taper_width) space. An isotropic GP kernel can learn the second and
# not the first. y_gap = mm_width/4 - taper_width exactly, so every existing
# result reparameterizes for free.
TAPER_MIN = WAVEGUIDE_WIDTH + EPS
TAPER_MAX = max(MM_WIDTHS) / 4.0 - Y_GAP_MIN


def taper_bounds(mm_width):
    """Valid taper_width at this width: the mirror of y_gap_bounds()."""
    return TAPER_MIN, mm_width / 4.0 - Y_GAP_MIN


# Slack allowed between a GP-proposed taper_width and this width's valid range
# before the candidate is rejected outright. Clipping instead of rejecting is
# what makes many candidates pile onto the y_gap floor: the GP's Real dimension
# spans the WIDEST width's taper range (up to 1.65), so at W=5.0-6.0 a large
# fraction of it is out of range and clips to the same y_gap=0.10 point.
TAPER_CLIP_SLACK = 0.02


def taper_to_ygap(mm_width, taper_width):
    """y_gap for this (width, taper_width), or None if out of range.

    A taper_width within TAPER_CLIP_SLACK of the range is clipped in (rounding
    to 2 decimals would land there anyway); anything further out is rejected by
    the caller rather than silently collapsed onto the boundary.
    """
    lo, hi = taper_bounds(mm_width)
    t = float(taper_width)
    if t < lo - TAPER_CLIP_SLACK or t > hi + TAPER_CLIP_SLACK:
        return None
    t = min(max(t, lo), hi)
    return round(mm_width / 4.0 - t, 2)


def ygap_to_taper(mm_width, y_gap):
    return mm_width / 4.0 - y_gap


def frac_to_ygap(mm_width, frac):
    lo, hi = y_gap_bounds(mm_width)
    return round(lo + frac * (hi - lo), 2)


def ygap_to_frac(mm_width, y_gap):
    lo, hi = y_gap_bounds(mm_width)
    return float(np.clip((y_gap - lo) / (hi - lo), 0.0, 1.0))


def load_json(path, default):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return default


def save_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=2)
    os.replace(tmp, path)


def npz_path(w, y):
    return os.path.join(OUT_DIR, f"eval_W{w:.2f}_y{y:.2f}.npz")


def read_eval(w, y):
    p = npz_path(w, y)
    if not os.path.exists(p):
        return None
    with np.load(p) as d:
        return {
            "mm_width": float(w), "y_gap": float(y),
            "taper_width": float(d["taper_width"]),
            "objective": float(d["objective"]),
            "max_dphi": float(d["max_dphi"]),
            "min_cmrr": float(d["min_cmrr"]),
            "best_delta_L": float(d["best_delta_L"]),
            "best_mm_length": float(d["best_mm_length"]),
        }


def absorb_disk(results):
    known = {(r["mm_width"], r["y_gap"]) for r in results}
    added = 0
    for fn in os.listdir(OUT_DIR):
        if not (fn.startswith("eval_W") and fn.endswith(".npz")):
            continue
        try:
            w = float(fn[6:10])
            y = float(fn[12:16])
        except ValueError:
            continue
        # Reject anything with a suffix past the y_gap field -- e.g. round 8's
        # eval_W6.50_y0.56_tl4.00.npz. Those are a DIFFERENT taper_length at
        # the same (mm_width, y_gap); absorbing one would overwrite the
        # default-taper point's objective with another geometry's.
        if fn != f"eval_W{w:.2f}_y{y:.2f}.npz":
            continue
        if (w, y) in known:
            continue
        rec = read_eval(w, y)
        if rec is not None:
            rec["source"] = "disk"
            results.append(rec)
            known.add((w, y))
            added += 1
    return added


def dedupe_results(results):
    """Belt-and-suspenders against the absorb/harvest double-count class of
    bug -- keeps the LAST entry per (mm_width, y_gap) (a "bo"-sourced record
    carries a jobid, strictly more informative than a bare "disk" one)."""
    by_key = {}
    for r in results:
        by_key[(r["mm_width"], r["y_gap"])] = r
    return list(by_key.values())


def sacct_state(jobid):
    try:
        out = subprocess.run(
            ["sacct", "-j", str(jobid), "-o", "State", "-n", "-X", "-P"],
            capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception as exc:
        return f"UNKNOWN({exc})"
    return out.splitlines()[0].strip() if out else "UNKNOWN"


def harvest(results, pending):
    # known BEFORE this call, e.g. via absorb_disk() -- a job whose npz was
    # already picked up that way is fully accounted for; appending it again
    # here would double-count it (found the hard way: cmd_status's own
    # absorb_disk-then-harvest order silently duplicated 3/4 of round 1's
    # points, 1.75x-weighting them in the GP fit before this fix).
    known = {(r["mm_width"], r["y_gap"]) for r in results}
    still_pending = []
    finished_ok, finished_bad = 0, 0
    for p in pending:
        st = sacct_state(p["jobid"])
        if st in ("RUNNING", "PENDING", "REQUEUED", "SUSPENDED", "UNKNOWN", ""):
            p["last_state"] = st
            still_pending.append(p)
            continue
        if st == "COMPLETED":
            if (p["mm_width"], p["y_gap"]) in known:
                finished_ok += 1
                continue
            rec = read_eval(p["mm_width"], p["y_gap"])
            if rec is not None:
                rec["source"] = "bo"
                rec["jobid"] = p["jobid"]
                results.append(rec)
                known.add((p["mm_width"], p["y_gap"]))
                finished_ok += 1
                continue
            print(f"WARNING: job {p['jobid']} COMPLETED but no npz at "
                  f"W={p['mm_width']:.2f} y={p['y_gap']:.2f} -- dropping, "
                  f"not resubmitting automatically.", file=sys.stderr)
            finished_bad += 1
            continue
        # FAILED / TIMEOUT / CANCELLED / OUT_OF_MEMORY etc.
        print(f"job {p['jobid']} (W={p['mm_width']:.2f} y={p['y_gap']:.2f}) "
              f"ended {st} -- per standing policy, NOT auto-resubmitted. "
              f"Report and wait for direction.", file=sys.stderr)
        finished_bad += 1
    return still_pending, finished_ok, finished_bad


def build_optimizer(results):
    from skopt import Optimizer
    from skopt.space import Categorical, Real
    dims = [Categorical(MM_WIDTHS, name="mm_width"),
            Real(TAPER_MIN, TAPER_MAX, name="taper_width")]
    # random_state=None (not a fixed seed): the Optimizer is rebuilt fresh
    # every invocation (state isn't persisted, only replayed via tell()
    # below), and while under n_initial_points its random-design draws do
    # NOT depend on told data at all -- a fixed seed here reproduced round
    # 1's exact 4 candidates again for round 2 (caught before submitting,
    # all 4 skipped as duplicates, zero cluster cost). Reproducibility of
    # the random-exploration phase isn't a requirement this campaign needs;
    # getting new candidates every round is.
    opt = Optimizer(dims, base_estimator="GP", acq_func="EI",
                     n_initial_points=8, random_state=None)
    if results:
        X = [[r["mm_width"], ygap_to_taper(r["mm_width"], r["y_gap"])]
             for r in results]
        y = [r["objective"] for r in results]
        opt.tell(X, y)
    return opt


def too_close(mm_width, y_gap, taken):
    """True if (mm_width, y_gap) is within MIN_YGAP_SEP of an already-taken
    point at the same width. Widths are categorical, so only same-width
    neighbours collide."""
    return any(w == mm_width and abs(y - y_gap) < MIN_YGAP_SEP
               for (w, y) in taken)


def explore_fill(taken, n):
    """Farthest-point fill: pick n lattice candidates that maximise the minimum
    distance to everything already taken, in (width-index, y_gap-frac) space.

    Used only when the GP's oversampled ask() cannot supply MIN_YGAP_SEP-distinct
    candidates -- i.e. when the acquisition has collapsed onto an already-measured
    neighbourhood. Pure exploration is the right fallback there: a round that
    submits fewer jobs than asked costs a whole extra queue cycle, and the GP has
    nothing new to learn from a 0.01 um re-measurement anyway.
    """
    lattice = [(w, y) for w in MM_WIDTHS for y in y_gaps_for(w)]
    picked = []
    for _ in range(n):
        best, best_d = None, -1.0
        for (w, y) in lattice:
            if too_close(w, y, taken) or too_close(w, y, picked):
                continue
            wi = MM_WIDTHS.index(w) / max(1, len(MM_WIDTHS) - 1)
            t = (ygap_to_taper(w, y) - TAPER_MIN) / (TAPER_MAX - TAPER_MIN)
            d = min((((MM_WIDTHS.index(tw) / max(1, len(MM_WIDTHS) - 1)) - wi) ** 2
                     + ((ygap_to_taper(tw, ty) - TAPER_MIN)
                        / (TAPER_MAX - TAPER_MIN) - t) ** 2) ** 0.5
                    for (tw, ty) in list(taken) + picked) if (taken or picked) else 1.0
            if d > best_d:
                best, best_d = (w, y), d
        if best is None:
            break
        picked.append(best)
    return picked


def propose_candidates(opt, results, q, seeds=()):
    """Return up to q (mm_width, y_gap) candidates, each MIN_YGAP_SEP-distinct
    from every evaluated point and from each other.

    Order: hand-picked `seeds` first (physics-motivated points the acquisition
    cannot reach), then GP-ranked candidates, then farthest-point exploration
    only to fill a short round.
    """
    taken = [(r["mm_width"], r["y_gap"]) for r in results]
    picked = []
    for (w, y) in seeds:
        w, y = float(w), round(float(y), 2)
        if too_close(w, y, taken) or too_close(w, y, picked):
            print(f"BO: skip seed W={w:.2f} y_gap={y:.2f} -- within "
                  f"{MIN_YGAP_SEP} um of an evaluated/picked point.")
            continue
        lo, hi = y_gap_bounds(w)
        if not (lo - 1e-9 <= y <= hi + 1e-9):
            print(f"BO: skip seed W={w:.2f} y_gap={y:.2f} -- outside the "
                  f"valid range ({lo:.2f}, {hi:.2f}).")
            continue
        picked.append((w, y))
    if len(picked) >= q:
        return picked[:q]
    asked = opt.ask(n_points=q * ASK_OVERSAMPLE, strategy="cl_min")
    for mm_width, taper_width in asked:
        if len(picked) >= q:
            break
        mm_width = float(mm_width)
        y_gap = taper_to_ygap(mm_width, float(taper_width))
        if y_gap is None:
            print(f"BO: skip GP candidate W={mm_width:.2f} "
                  f"taper={float(taper_width):.3f} -- outside this width's "
                  f"valid taper range {taper_bounds(mm_width)[0]:.2f}-"
                  f"{taper_bounds(mm_width)[1]:.2f}.")
            continue
        if too_close(mm_width, y_gap, taken):
            print(f"BO: skip GP candidate W={mm_width:.2f} y_gap={y_gap:.2f} "
                  f"(taper={ygap_to_taper(mm_width, y_gap):.3f}) -- within "
                  f"{MIN_YGAP_SEP} um of an evaluated point.")
            continue
        if too_close(mm_width, y_gap, picked):
            continue
        picked.append((mm_width, y_gap))
    if len(picked) < q:
        short = q - len(picked)
        fill = explore_fill(taken + picked, short)
        if fill:
            print(f"BO: GP ask exhausted after {len(picked)}/{q} distinct "
                  f"candidate(s) -- filling {len(fill)} by farthest-point "
                  f"exploration.")
        picked.extend(fill)
    return picked


def submit_point(w, y):
    env = os.environ.copy()
    name = f"beef_bo_W{w:.2f}_y{y:.2f}"
    cmd = ["sbatch", f"--export=ALL,MMI4X4_W={w:.2f},MMI4X4_YGAP={y:.2f}",
           f"--job-name={name}", BATCH_SCRIPT]
    out = subprocess.run(cmd, cwd=_ROOT, capture_output=True, text=True, env=env)
    if out.returncode != 0:
        raise RuntimeError(f"sbatch failed for W={w} y={y}: {out.stderr}")
    # "Submitted batch job 12345678"
    jobid = out.stdout.strip().split()[-1]
    return jobid


def report_best(results):
    if not results:
        print("BO: no results yet.")
        return
    best = min(results, key=lambda r: r["objective"])
    print(f"BO: {len(results)} evaluated points so far.")
    print(f"BO: best  mm_width={best['mm_width']:.2f} y_gap={best['y_gap']:.2f} "
          f"taper_width={best['taper_width']:.3f}")
    print(f"BO:       mm_length={best['best_mm_length']:.4f} "
          f"(dL={best['best_delta_L']:+.2f})")
    print(f"BO:       max|dphi|={best['max_dphi']:.3f} deg  "
          f"min CMRR={best['min_cmrr']:.3f} dBe  objective={best['objective']:.4f}")
    print("BO:       spec: |dphi| <= 5 deg, CMRR >= 20 dBe")
    if (abs(best["best_mm_length"] - LENGTH_MIN) < 1e-6
            or abs(best["best_mm_length"] - LENGTH_MAX) < 1e-6):
        print(f"BO:       WARNING: best point's length sweep hit the "
              f"{LENGTH_MIN}-{LENGTH_MAX}um edge -- not a trustworthy optimum "
              f"at this width, widen MM_LENGTH_GRID.")
    if abs(best["y_gap"] - Y_GAP_MIN) < 1e-9:
        print("BO:       WARNING: best point sits on the y_gap floor (0.10), a "
              "modelling choice not a physical limit -- consider lowering it.")


def cmd_status(args):
    results = load_json(RESULTS_PATH, [])
    pending = load_json(PENDING_PATH, [])
    added = absorb_disk(results)
    pending, ok, bad = harvest(results, pending)
    results = dedupe_results(results)
    save_json(RESULTS_PATH, results)
    save_json(PENDING_PATH, pending)
    print(f"BO: absorbed {added} pre-existing result(s) from disk; "
          f"harvested {ok} completed, {bad} failed/dropped this check.")
    if pending:
        print(f"BO: {len(pending)} job(s) still outstanding:")
        for p in pending:
            print(f"    job {p['jobid']}  W={p['mm_width']:.2f} "
                  f"y_gap={p['y_gap']:.2f}  state={p.get('last_state','?')}")
    else:
        print("BO: no jobs outstanding -- ready to propose the next round.")
    report_best(results)


def cmd_propose(args):
    results = load_json(RESULTS_PATH, [])
    pending = load_json(PENDING_PATH, [])
    absorb_disk(results)
    pending, ok, bad = harvest(results, pending)
    results = dedupe_results(results)
    save_json(RESULTS_PATH, results)
    save_json(PENDING_PATH, pending)
    if pending:
        print(f"BO: {len(pending)} job(s) still outstanding from the current "
              f"round -- not proposing more until they clear. Run `status` "
              f"or re-run `propose` later.")
        for p in pending:
            print(f"    job {p['jobid']}  W={p['mm_width']:.2f} "
                  f"y_gap={p['y_gap']:.2f}  state={p.get('last_state','?')}")
        return

    q = args.batch
    seeds = []
    for spec in (args.seed or []):
        try:
            w_s, y_s = spec.split(",")
            seeds.append((float(w_s), float(y_s)))
        except ValueError:
            raise SystemExit(f"--seed expects W,y_gap (e.g. 7.0,0.65), got {spec!r}")
    opt = build_optimizer(results)
    candidates = propose_candidates(opt, results, q, seeds=seeds)
    if len(candidates) < q:
        print(f"BO: WARNING: only {len(candidates)}/{q} distinct candidate(s) "
              f"available -- the (mm_width, y_gap) lattice is nearly exhausted "
              f"at MIN_YGAP_SEP={MIN_YGAP_SEP}.")
    new_pending = []
    for mm_width, y_gap in candidates:
        if args.dry_run:
            print(f"BO: DRY RUN would submit  W={mm_width:.2f}  "
                  f"y_gap={y_gap:.2f}  "
                  f"(taper_width={ygap_to_taper(mm_width, y_gap):.3f})")
            continue
        jobid = submit_point(mm_width, y_gap)
        print(f"BO: submitted job {jobid}  W={mm_width:.2f}  "
              f"y_gap={y_gap:.2f}  "
              f"(taper_width={ygap_to_taper(mm_width, y_gap):.3f})")
        new_pending.append({"mm_width": mm_width, "y_gap": y_gap,
                             "jobid": jobid, "submitted_at": time.time(),
                             "last_state": "PENDING"})
    if args.dry_run:
        print("BO: DRY RUN -- nothing submitted, pending state untouched.")
        report_best(results)
        return
    save_json(PENDING_PATH, new_pending)
    print(f"BO: round submitted, {len(new_pending)} new job(s). "
          f"Re-run `status`/`propose` once they complete.")
    report_best(results)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sp = sub.add_parser("status")
    sp.set_defaults(func=cmd_status)
    pp = sub.add_parser("propose")
    pp.add_argument("--batch", type=int, default=BATCH_SIZE_DEFAULT)
    pp.add_argument("--seed", action="append", metavar="W,YGAP",
                    help="hand-picked point to include in this round, ahead of "
                         "the GP's own candidates; repeatable "
                         "(e.g. --seed 7.0,0.65)")
    pp.add_argument("--dry-run", action="store_true",
                    help="show the candidates this round would submit, "
                         "submit nothing, leave pending state untouched")
    pp.set_defaults(func=cmd_propose)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
