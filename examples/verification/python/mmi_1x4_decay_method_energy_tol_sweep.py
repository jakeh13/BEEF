"""Whole-device FDTD shutoff study: decay_method="energy" and a decay_tol
sweep, validated against the trusted decay_method="dft" reference.

validation_markdown/mmi_1x4_sin300_validation.md step 8 check 4 / check 3.

WHY. Every FDTD run on this device uses decay_method="dft" (the FDTDSolver
default -- no 1x4 script sets it) with decay_tol=1e-11.
mp.stop_when_dft_decayed() requires EVERY field component at EVERY frequency
on EVERY DFT monitor to individually fall below tol, and no maximum_run_time
ceiling is ever set. That is the exact configuration that stalled the 4x4
90-degree hybrid past 660,000+ timesteps with no sign of decaying, which is
why decay_method="energy" was added 2026-08-18 -- still marked
"investigative, not yet validated as a general replacement" in its own
docstring.

The risk is live right now: the band-centre whole-device reference (job
12029541) carried ~4 DFT components and finished in 40,420 steps; the
full-band check-4 run (job 12260631) carries ~44 (11 wavelengths x 4 output
monitors), all of which must decay independently.

WHAT MAKES THIS TESTABLE RATHER THAN A GUESS. We have a trusted reference at
exactly this geometry/resolution: job 12029541, decay_method="dft",
decay_tol=1e-11, res=30, N_WL=1, whole device --

    T per output = [0.23997513751844815, 0.2424475379509313,
                    0.24244750636842585, 0.2399752894983026]
    Total T      = 0.9648

so every run below is scored against a real number, not against itself.

THE decay_tol SUBTLETY, which is easy to get wrong. `decay_tol` does NOT
mean the same thing to the two methods:
  * "dft"    -> mp.stop_when_dft_decayed(tol=decay_tol): a RELATIVE
                threshold applied per DFT component.
  * "energy" -> mp.stop_when_energy_decayed(decay_by=decay_tol): the factor
                the BULK field energy in the whole cell must fall by.
1e-11 is a reasonable per-component DFT threshold but an absurdly strict
bulk-energy target; 1e-5 is the sensible scale for the energy criterion, and
is this study's default (DECAY_TOL_DEFAULT). Carrying 1e-11 across to the
energy method unchanged would be comparing two different demands and calling
it a method comparison.

**Library defaults are NOT changed by this script.** FDTDSolver's own
decay_tol=1e-11 / decay_method="dft" defaults stay exactly as they are --
changing them would silently move every other device in the repo and
invalidate their caches. 1e-5 is the default *for this study*.

WHAT IS SWEPT. decay_method="energy" at decay_tol in TOL_LIST (loosest
first, so the cheapest runs report before any expensive one). Each run
reports steps, wall time, T per output, total T, and the max per-port
deviation from the dft reference above -- which is the number that decides
whether the cheaper shutoff is trustworthy, not the speedup.

Run (see the batch file; the loosest tolerance alone is a good smoke pass):
    MMI1X4_DECAY_TOLS=1e-3 python .../mmi_1x4_decay_method_energy_tol_sweep.py
    python .../mmi_1x4_decay_method_energy_tol_sweep.py
"""
import os
import time

import numpy as np

from beef import MMI, FDTDSolver, SIN_300

# --- the trusted dft reference (job 12029541) -------------------------------
REF_T = np.array([0.23997513751844815, 0.2424475379509313,
                  0.24244750636842585, 0.2399752894983026])
REF_TOTAL = float(REF_T.sum())
REF_STEPS = 40420
REF_HOURS = 41.032

# --- pinned to the reference run's own config, unchanged --------------------
RESOLUTION = 30.0
WL_CENTER = 1.55
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 1.1, 2.2, 5.0
WAVEGUIDE_OFFSET, MM_WIDTH = 3.1, 12.0
MM_LENGTH = 39.90
FDTD_SIM_WIDTH, FDTD_SIM_HEIGHT = MM_WIDTH + 4.0, 4.0
_MONITOR_WIDTH = 2.8

# --- what this study varies -------------------------------------------------
DECAY_TOL_DEFAULT = 1e-5
_tols = os.environ.get("MMI1X4_DECAY_TOLS")
TOL_LIST = ([float(t) for t in _tols.split(",")] if _tols
            else [1e-3, DECAY_TOL_DEFAULT, 1e-7])
ENERGY_DT = float(os.environ.get("MMI1X4_ENERGY_DT", 20.0))

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30, use_neighbor_aperture=True,
)
wavelengths = np.array([WL_CENTER])

print("=" * 72)
print("DECAY-METHOD / DECAY-TOL STUDY -- whole-device 1x4, res=30, N_WL=1")
print(f"  reference (dft, decay_tol=1e-11, job 12029541):")
print(f"    T = {list(np.round(REF_T, 6))}")
print(f"    total T = {REF_TOTAL:.4f}   ({REF_STEPS} steps, {REF_HOURS:.2f} hrs)")
print(f"  sweeping decay_method='energy', energy_decay_dt={ENERGY_DT}")
print(f"  decay_tol values: {TOL_LIST}  (study default {DECAY_TOL_DEFAULT})")
print("  NOTE: decay_tol is a per-DFT-component threshold for 'dft' but a")
print("        BULK-ENERGY decay FACTOR for 'energy' -- not the same demand.")
print("=" * 72, flush=True)

rows = []
for tol in TOL_LIST:
    print(f"\n--- decay_method='energy', decay_tol={tol:g} ---", flush=True)
    solver = FDTDSolver(
        res=RESOLUTION, sim_width=FDTD_SIM_WIDTH, sim_height=FDTD_SIM_HEIGHT,
        dpml=0.5, source_margin=0.0, monitor_margin=0.25, three_d=True,
        num_input_modes=1, num_output_modes=1,
        output_monitor_centers=list(mmi.output_ys),
        output_monitor_widths=[_MONITOR_WIDTH] * 4,
        material_stack=SIN_300,
        decay_tol=tol,
        decay_method="energy",
        energy_decay_dt=ENERGY_DT,
    )
    t0 = time.time()
    try:
        smat = solver.solve(mmi, wavelengths, port="top", input_mode=1, save_dir=None)
    except Exception as exc:
        print(f"  FAILED: {type(exc).__name__}: {exc}", flush=True)
        rows.append((tol, None, None, None, None))
        continue
    elapsed = time.time() - t0
    T = np.array([float(np.abs(smat.transmission(port=f"out_{k}", in_mode=1,
                                                 out_mode=1))[0]) for k in range(4)])
    dev = float(np.max(np.abs(T - REF_T)))
    print(f"  {elapsed:.1f} s ({elapsed/3600:.3f} hrs)", flush=True)
    print(f"  T       = {list(np.round(T, 6))}")
    print(f"  total T = {T.sum():.4f}   (reference {REF_TOTAL:.4f}, "
          f"delta {T.sum()-REF_TOTAL:+.4f})")
    print(f"  MAX PER-PORT DEVIATION FROM REFERENCE: {dev:.2e}", flush=True)
    rows.append((tol, elapsed, T, T.sum(), dev))

print("\n" + "=" * 72)
print("SUMMARY -- is the cheaper shutoff trustworthy, and how much cheaper?")
print(f"{'decay_tol':>12} {'hrs':>8} {'speedup':>9} {'total T':>9} {'max dev':>10}")
print(f"{'dft 1e-11':>12} {REF_HOURS:8.2f} {'1.0x':>9} {REF_TOTAL:9.4f} "
      f"{'(reference)':>10}")
for tol, elapsed, T, tot, dev in rows:
    if elapsed is None:
        print(f"{tol:12g} {'FAILED':>8}")
        continue
    print(f"{tol:12g} {elapsed/3600:8.2f} {REF_HOURS/(elapsed/3600):8.1f}x "
          f"{tot:9.4f} {dev:10.2e}")
print("=" * 72)
print("A deviation at or below ~1e-4 per port would make 'energy' a safe")
print("substitute at this device; a large speedup with a large deviation")
print("means the run stopped before the fields had actually settled.")
