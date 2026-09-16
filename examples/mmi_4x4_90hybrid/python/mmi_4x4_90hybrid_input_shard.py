"""4x4 90 hybrid -- ONE FDTD input-section solve, so the four can run at once.

mmi_4x4_90hybrid_hybrid.py solves its four input sections in a plain `for k in
range(4)` loop. Each is an independent 3-D FDTD run costing ~1.75 h on 24
ranks (job 13024689: 6:58:47 total, the four solves evenly spaced through it
with the body a cache hit), so the loop serialises ~7 h of work that has no
data dependence at all -- M_in[k] never reads M_in[k-1].

This script runs exactly ONE of them and writes it to the shared production
cache, so four jobs can populate all four entries concurrently and the example
then cache-hits every one. Unlike the body's neff sharding
(mmi_4x4_90hybrid_body_neff_shard.py) there is nothing to merge: input_idx is
already part of solve_section()'s own cache key, so each shard writes the very
entry the example will look up.

Time-stepping scales ~98% with ranks (48->96 gave 1.95x, job 12217184), so
unlike the body shards these DO want many ranks each.

  --verify    Ask for the 3-wavelength solves with cache_only=True. They are
              known to be cached, so a hit proves this script's geometry and
              solver reconstruction are byte-identical to the example's --
              using the library's own key construction rather than a
              reimplementation of it. Costs nothing. RUN THIS FIRST: a miss
              means the shards would populate keys the example never reads,
              and the 20 h would simply be spent again.
  --shard K   Solve input section K at the 11-wavelength default.
"""
import argparse
import os
import sys

import numpy as np

from beef import MMI, FDTDSolver, SOI_220
from beef.cache import CacheMissError

_HERE = os.path.dirname(os.path.abspath(__file__))
SAVE_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# --- pinned to mmi_4x4_90hybrid_hybrid.py; --verify is what proves it ------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.29, 1.33, 11
N_MODES = 16
WAVEGUIDE_WIDTH, MM_WIDTH = 0.35, 8.5
TAPER_WIDTH, TAPER_LENGTH = 1.30, 4.0
MM_LENGTH = float(os.environ.get("MMI4X4_MM_LENGTH", 159.30))
BODY_STUB = 0.25

WL11 = np.linspace(WL_MIN, WL_MAX, N_WL)
WL3 = np.linspace(WL_MIN, WL_MAX, 3)

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=MM_WIDTH / 4.0,
    waveguide_height=0.22, use_neighbor_aperture=True,
)
input_solver = FDTDSolver(
    res=RESOLUTION, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SOI_220, decay_tol=1e-4, decay_method="energy",
    energy_decay_dt=20.0,
)
_half_body, _half_total = mmi.mm_length / 2, mmi.total_length / 2
B0, B1 = -_half_total, -_half_body + BODY_STUB


def solve(k, wavelengths, cache_only):
    return input_solver.solve_section(
        mmi, wavelengths, B0, B1,
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=input_solver.waveguide_length, stub_right=0.0,
        input_idx=k, save_dir=SAVE_DIR, cache_only=cache_only,
    )


def do_verify():
    print(f"input section x in [{B0:+.3f}, {B1:+.3f}] um, {N_MODES} out-modes")
    ok = True
    for k in range(4):
        try:
            M = solve(k, WL3, True)
            print(f"  input_idx={k}: 3-wl cache HIT {M.shape}")
        except CacheMissError:
            print(f"  input_idx={k}: 3-wl cache MISS")
            ok = False
    if ok:
        print("\ngeometry/solver reconstruction is byte-identical to the example's")
    else:
        print("\nMISS -- do NOT run the shards; they would populate unread keys")
    for k in range(4):
        try:
            solve(k, WL11, True)
            print(f"  input_idx={k}: 11-wl ALREADY CACHED")
        except CacheMissError:
            print(f"  input_idx={k}: 11-wl not yet solved (expected)")
    return 0 if ok else 1


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--verify", action="store_true")
    p.add_argument("--shard", type=int)
    a = p.parse_args()
    if a.verify:
        sys.exit(do_verify())
    if a.shard is None:
        p.error("pick --verify or --shard K")
    print(f"input shard {a.shard}: {N_WL} wavelengths {WL_MIN}-{WL_MAX} um, "
          f"res {RESOLUTION}", flush=True)
    M = solve(a.shard, WL11, False)
    print(f"input_idx={a.shard} solved: {M.shape}")
