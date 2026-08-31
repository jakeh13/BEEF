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

SEARCH SPACE. dimensions = [Categorical(MM_WIDTHS), Real(0, 1) frac], with
y_gap = Y_GAP_MIN + frac * (mm_width/4 - WAVEGUIDE_WIDTH - Y_GAP_MIN - EPS)
-- the same decoupling trick the grid module's docstring explains, just
continuous in frac instead of lattice steps. EPS keeps rounding to 2 decimals
(wg_eval's own npz-naming precision) from landing exactly on the inverse-taper
boundary.

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
    MM_WIDTHS, WAVEGUIDE_WIDTH, Y_GAP_MIN, LENGTH_MIN, LENGTH_MAX)

EPS = 0.01
BATCH_SIZE_DEFAULT = 4
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
            Real(0.0, 1.0, name="y_gap_frac")]
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
        X = [[r["mm_width"], ygap_to_frac(r["mm_width"], r["y_gap"])]
             for r in results]
        y = [r["objective"] for r in results]
        opt.tell(X, y)
    return opt


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
    opt = build_optimizer(results)
    known = {(r["mm_width"], r["y_gap"]) for r in results}
    asked = opt.ask(n_points=q, strategy="cl_min")
    new_pending = []
    for mm_width, frac in asked:
        mm_width = float(mm_width)
        y_gap = frac_to_ygap(mm_width, float(frac))
        if (mm_width, y_gap) in known:
            print(f"BO: candidate W={mm_width:.2f} y_gap={y_gap:.2f} already "
                  f"evaluated (rounding collision) -- skipping resubmission.")
            continue
        jobid = submit_point(mm_width, y_gap)
        print(f"BO: submitted job {jobid}  W={mm_width:.2f}  y_gap={y_gap:.2f} "
              f"(frac={frac:.3f})")
        new_pending.append({"mm_width": mm_width, "y_gap": y_gap,
                             "jobid": jobid, "submitted_at": time.time(),
                             "last_state": "PENDING"})
        known.add((mm_width, y_gap))
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
    pp.set_defaults(func=cmd_propose)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
