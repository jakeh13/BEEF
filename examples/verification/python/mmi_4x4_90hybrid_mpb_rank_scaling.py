"""Does the MPB phase of calibration actually use the MPI ranks it is given?

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-22.

WHY. A phase breakdown of this device's own logs showed the full-device
FDTD run (job 12079104) spends only 4.2% of wall time outside time-stepping
-- which is why it scaled at 98% efficiency to 96 ranks -- while the
per-mode CALIBRATION workers spend 47-78% outside time-stepping, in MPB
mode solves and eigenmode decomposition. That independently reproduces
step 9 of section_mode_beta_multi_monitor_validation.md ("MPB eigenmode
decomposition is ~89% of calibration cost").

If meep runs those MPB eigenproblems redundantly on every rank -- each
process solving the identical problem -- then during that majority of the
runtime N-1 of N cores are duplicating work rather than contributing, and
calibration jobs should be asking for far fewer ranks than they do (the
16-task array used 24 ranks/task = 384 cores).

That redundancy is INFERRED, not yet measured. This measures it.

METHOD. Time `section_mode_neff()` -- pure MPB, no FDTD time-stepping at
all -- on this device's own body cross-section, at the same n_modes
production calibration uses. Run the same script at two rank counts and
compare:
    wall(8 ranks) ~= wall(24 ranks)  -> MPB is redundant, ranks do nothing
    wall(8 ranks) ~= 3x wall(24)     -> MPB scales, ranks are being used

**No `save_dir` is passed, deliberately.** Every other call site caches
this, and a cache hit would return instantly and measure nothing. This
recomputes from scratch every time, which is the whole point. It therefore
writes NOTHING to the cache either, so it cannot pollute or collide with
the production entries the calibration array already wrote.

Run:
    python mmi_4x4_90hybrid_mpb_rank_scaling.py
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmi_4x4_90hybrid_offset_sweep as base
from beef import FDTDSolver


def main():
    n_modes = base.N_MODES
    mmi = base.make_mmi(base.BASELINE_OFFSET)
    boundaries, _stub = base.boundaries_for(mmi)
    x0, x1 = boundaries[1], boundaries[2]
    mid = 0.5 * (x0 + x1)
    seed_span = min(0.5, (x1 - x0) * 0.25)

    solver = FDTDSolver(
        res=40.0, sim_width=base.MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=base.SOI_220, decay_tol=1e-11,
    )

    # 3 anchors, matching the production seed path's own seed_anchors=3.
    wl_seed = np.linspace(base.WL_MIN, base.WL_MAX, 3)

    try:
        import meep as mp
        nproc = mp.count_processors()
    except Exception:
        nproc = -1
    print(f"MPBSCALE: ranks={nproc} n_modes={n_modes} wavelengths={len(wl_seed)}",
          flush=True)

    t0 = time.time()
    neff = solver.section_mode_neff(
        mmi, wl_seed, mid - seed_span / 2.0, mid + seed_span / 2.0, n_modes,
        monitor_center=0.0, monitor_size=solver.sim_width,
        # NO save_dir -- see module docstring. Must not hit cache.
    )
    dt = time.time() - t0

    print(f"MPBSCALE: section_mode_neff took {dt:.1f} s ({dt/60:.2f} min) "
          f"on {nproc} ranks", flush=True)
    print(f"MPBSCALE: neff[0] = {np.asarray(neff)[0]}", flush=True)


if __name__ == "__main__":
    main()
