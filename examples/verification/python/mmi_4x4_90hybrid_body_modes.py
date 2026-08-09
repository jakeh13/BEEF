"""4x4 MMI 90-degree hybrid -- Step 3a: body mode count at the settled
mm_width=6.0 (Step 2), O-band, SOI_220.

validation_markdown/mmi_4x4_90_hybrid_validation.md Step 3: before the real
FDTD/EME sweep can be built, N_MODES for the interface transfer matrices
must be measured, not guessed -- every other device in this repo (2x2: 14
of 15, 1x4 SiN: 20 of 23, modal mux: 14 of 16) pins this from a real audit,
not a default. Step 1c only ever solved 2 modes at this cross-section (just
enough for L_pi); this asks for more and reports the real guided count.

Reuses mmi_4x4_90hybrid_cross_section.py's solve_modes() (same
UniformCrossSection + EMESolver.compute_modes() pipeline, no FDTD).

Run:
    python examples/verification/python/mmi_4x4_90hybrid_body_modes.py
"""
from __future__ import annotations

import time

import numpy as np

from mmi_4x4_90hybrid_cross_section import solve_modes, N_CLAD, guided, WL, RES

MM_WIDTH = 6.0   # Step 2's settled candidate
# Calibration point, not the final N_MODES: job 11706986 requested 16 and
# TIMED OUT at the full 60 min with ZERO output -- the 2-mode point (job
# 11693043, Step 1c) took 874.6s, so "16 modes" was not simply a few times
# more expensive as guessed. 8 is a data point to establish the real
# per-mode scaling before committing more compute to a final number.
#
# 2026-08-08: job 11707935/11707835 confirmed all 8 requested modes come
# back GUIDED with a healthy, still-widening n_eff margin (2.962->2.833,
# cladding 1.4468) -- meaning the true guided-mode ceiling for this
# cross-section was never actually found, only a lower bound established.
# The taper-mode-count investigation found the SAME open question at the
# taper's own 0.8um-wide end mattered a lot (real guided count 6, not the
# 2 an earlier buggy check wrongly suggested) -- extending here to find
# the real ceiling for the body too, at the ~416 s/mode rate measured
# above (8 modes, 3324.7s): 20 modes budgeted at ~2.3h, via batch (not the
# 60-min timeout that killed job 11706986).
NUM_MODES_REQUESTED = 20
SIM_W = MM_WIDTH + 14.0   # generous PML/margin clearance for a wide mode request
SIM_H = 2.0               # Step 1b's converged sz


def main():
    print(f"[3a] body mode count -- mm_width={MM_WIDTH} um, wl={WL} um, "
          f"res={RES}, requesting {NUM_MODES_REQUESTED} modes, "
          f"n_clad(SiO2)={N_CLAD:.6f}")
    t0 = time.time()
    neffs = solve_modes([(0.0, MM_WIDTH)], num_modes=NUM_MODES_REQUESTED,
                         sim_w=SIM_W, sim_h=SIM_H)
    g = guided(neffs)
    print(f"requested={NUM_MODES_REQUESTED}  guided={len(g)}  "
          f"({time.time()-t0:.1f}s)")
    for i, n in enumerate(neffs):
        flag = "guided" if n > N_CLAD else "UNGUIDED"
        print(f"  band {i:2d}  n_eff={n:.6f}  {flag}")


if __name__ == "__main__":
    main()
