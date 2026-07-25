"""TEMPORARY -- multi-rank check of the beef cache write gate
(dc_validation.md item 4). Delete with the rest of examples/verification/.

The unit tests in tests/test_cache_multirank.py pin the logic, but they cannot
reproduce the setting the bug lived in: 48 real ranks spread over 2 nodes,
writing one cache on shared storage. This does exactly that, importing meep so
am_really_master() and the barrier are the real ones, while running no
simulation at all -- seconds, not hours.

Must be run on >=2 nodes with save_dir on SHARED storage. A node-local path
(/tmp) silently defeats the whole point: each node then gets its own directory
and cannot race with the other.

    srun -A <account> -q inferno -N 2 --ntasks-per-node=24 -t 10 \\
         --mem-per-cpu=2G --mpi=pmi2 python \\
         verification/python/verify_cache_multirank.py <shared_dir>

Passes when index.jsonl has one line per save rather than one per rank, no
temp files survive, and the family pointer is readable. Exits non-zero
otherwise. Also reports whether any PID was reused across hosts -- the
mechanism originally blamed for the tasks-4/5 crash, unconfirmed as of
2026-07-25 (a 48-rank sample showed none), which is why the fix does not rely
on it.
"""
import glob
import json
import os
import socket
import sys
from collections import defaultdict

import numpy as np
import meep as mp

from beef.cache.store import lookup, lookup_superset, save_run, update_family_pointer

if len(sys.argv) < 2:
    raise SystemExit(f"usage: {os.path.basename(sys.argv[0])} <shared_save_dir>")
save_dir = sys.argv[1]
os.makedirs(save_dir, exist_ok=True)
rank_dir = os.path.join(save_dir, "_ranks")
N_ITER = 5
fam = {"solver_type": "FDTDSolver", "method": "solve_section", "geometry_class": "Taper"}

os.makedirs(rank_dir, exist_ok=True)
rank = int(mp.my_rank())
with open(os.path.join(rank_dir, f"{rank:03d}.json"), "w") as f:
    json.dump({"rank": rank, "host": socket.gethostname(), "pid": os.getpid()}, f)

for i in range(N_ITER):
    cfg = {"solver_type": "FDTDSolver", "method": "solve_section",
           "geometry_class": "Taper", "geometry_polygons": {},
           "call_kwargs": {"n_modes_out": 6, "iter": i}}
    folder = save_run(save_dir, cfg, np.zeros((3, 6, 1)))
    update_family_pointer(save_dir, fam, (6, 1), os.path.basename(folder))
    # Every rank must agree the entry exists, or ranks diverge on whether to
    # simulate next time and meep's collective calls go unmatched -> hang.
    assert lookup(save_dir, cfg) is not None, f"rank {rank} cannot see the save"
    assert lookup_superset(save_dir, fam, (2, 1)) is not None, f"rank {rank} lost the pointer"

mp.all_wait()
if mp.am_really_master():
    lines = [l for l in open(os.path.join(save_dir, "index.jsonl")) if l.strip()]
    for l in lines:
        json.loads(l)                      # raises if two ranks interleaved
    strays = glob.glob(os.path.join(save_dir, "**", "*.tmp.*"), recursive=True)
    pdir = os.path.join(save_dir, "_mode_family_pointers")
    with open(os.path.join(pdir, sorted(os.listdir(pdir))[0])) as f:
        pointer = json.load(f)             # raises if truncated -- the crash

    recs = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(rank_dir, "*.json")))]
    hosts = sorted({x["host"] for x in recs})
    bypid = defaultdict(set)
    for x in recs:
        bypid[x["pid"]].add(x["host"])
    dup = {p: h for p, h in bypid.items() if len(h) > 1}
    masters = sum(1 for _ in [1])          # only this rank got here

    print(f"\n  ranks           {mp.count_processors()} on {len(hosts)} hosts")
    print(f"  index.jsonl     {len(lines)} lines (expected {N_ITER}; pre-fix would be "
          f"{N_ITER * mp.count_processors()})")
    print(f"  stray tmp files {len(strays)}")
    print(f"  pointer dims    {pointer['dims']}")
    print(f"  PIDs shared across hosts: {len(dup)}"
          + (f"  e.g. {list(dup.items())[0]}" if dup else ""))
    ok = len(lines) == N_ITER and not strays and len(recs) == mp.count_processors()
    print("\nOK" if ok else "\nFAIL")
    sys.exit(0 if ok else 1)
