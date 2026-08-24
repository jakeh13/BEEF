"""Measure the real per-solve() cost of mmi_4x4_90hybrid_eme.py at the
CORRECTED N_MODES_BODY=16 / DELTA_L=-2.75 point, so the re-run's batch
budget is measured rather than extrapolated.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-21. The EME
script and its field-plot companion were missed by the 2026-08-09 body-mode
redo and sat at the superseded N_MODES_BODY=8 / DELTA_L=-1.50 point; both
are now corrected, and every EME artifact on disk has to be regenerated.

WHY A PROBE AND NOT AN EXTRAPOLATION. The last real EME run (job 11736983)
took 1:13:27 at N_MODES_BODY=8. Doubling the body basis does not scale by
any factor this repo has measured, and this device has already burned
budget twice on guessed EME costs -- once by extrapolating a ~0.3-0.6 dB
num_cells effect from a coarse proxy that measured 0.08-0.10 dB for real,
and once by reading an ISOLATED solve() call's 8x slowdown as a real cost
when it was mostly one-time JAX JIT compilation amortized across the real
script's 8 same-shape calls.

That JIT caveat is why this probe times TWO calls, not one:
    call 1 = JIT compile + solve      (one-time)
    call 2 = solve at steady state    (the number that actually scales)
    estimated real 8-call total ~= call1 + 7 * call2

TWO DELIBERATE SAFETY PROPERTIES:
  * Writes NOTHING to the script's own output_dir. `mmi_4x4_90hybrid_eme.py`
    hardcodes that path and its env vars change only solver params, never
    the path -- on 2026-08-07 two coarse smoke invocations overwrote job
    11726819's official results exactly this way. Any output here goes to a
    scratch path.
  * Shares the real `save_dir`, at the REAL config, so the two solves it
    completes are genuine cache entries the full re-run then hits instead of
    recomputing. The probe is not throwaway work -- it is the first 2 of the
    8 calls the real run needs.

Run:
    python mmi_4x4_90hybrid_eme_timing_probe.py
"""
import os
import time

import numpy as np

from beef import MMI, EMESolver, SOI_220

# Pinned to mmi_4x4_90hybrid_eme.py's CORRECTED constants. Duplicated here
# deliberately (this is a throwaway measurement, not a third place for the
# design point to drift) -- if that script moves, this probe is stale and
# should be re-derived, not edited.
RESOLUTION = 40.0
N_WL = 3
WL_MIN, WL_MAX = 1.26, 1.36
WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 6.0
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 0.8
TAPER_LENGTH = 1.0
L_PI = 108.9195
DELTA_L = -2.75
MM_LENGTH = 0.75 * L_PI + DELTA_L    # 78.9396
N_MODES_BODY = 16
N_MODES_TAPER = 12
N_CELLS_TAPER = 30

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.22, use_neighbor_aperture=True,
)

solver = EMESolver(
    num_modes=[1, N_MODES_TAPER, N_MODES_BODY, N_MODES_TAPER, 1],
    num_cells=[1, N_CELLS_TAPER, 1, N_CELLS_TAPER, 1],
    lengths=[0.5, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, 0.5],
    oxide_width=MM_WIDTH + 3.0,
    oxide_height=2.0,
    res=RESOLUTION,
    material_stack=SOI_220,
    mode_solver="meep",
)

print(f"PROBE: N_MODES_BODY={N_MODES_BODY} DELTA_L={DELTA_L} "
      f"MM_LENGTH={MM_LENGTH:.4f} num_cells={N_CELLS_TAPER} res={RESOLUTION}",
      flush=True)
print(f"PROBE: total_length={mmi.total_length:.4f} um, {N_WL} wavelengths", flush=True)

times = []
for n, (i, j) in enumerate([(0, 0), (0, 1)], start=1):
    t0 = time.time()
    solver.solve(
        mmi, wavelengths, ports=[f"out{5+j}_in{i}"],
        input_idx=i, output_idx=j, save_dir=save_dir,
    )
    dt = time.time() - t0
    times.append(dt)
    print(f"PROBE: call {n} (input_idx={i}, output_idx={j}) took {dt:.1f} s "
          f"({dt/60:.2f} min)", flush=True)

if len(times) == 2:
    jit = times[0] - times[1]
    est = times[0] + 7 * times[1]
    print(f"PROBE: one-time JIT/first-call overhead ~= {jit:.1f} s", flush=True)
    print(f"PROBE: steady-state per call ~= {times[1]:.1f} s "
          f"({times[1]/60:.2f} min)", flush=True)
    print(f"PROBE: ESTIMATED full 8-call run = call1 + 7*call2 = "
          f"{est:.1f} s = {est/3600:.2f} h", flush=True)
    print(f"PROBE: (2 of those 8 calls are now cached by this probe, so the "
          f"real run should cost ~{6*times[1]/3600:.2f} h of fresh work)",
          flush=True)
