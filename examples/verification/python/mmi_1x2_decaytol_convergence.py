"""1x2 MMI — where does transmission ACTUALLY converge vs decay time?

validation_markdown/hybrid_field_extraction_validation.md, 2026-08-31. Direct
port of mmi_4x4_90hybrid_decaytol_convergence.py's method to this device.

WHY. `decay_tol=1e-11` is inherited across this repo with **no recorded
justification anywhere**, and it is the dominant cost knob: the fully-cold
Step 0 re-baseline (job 12627227) burned 5 h and hit its wall without
finishing, on a device whose previous run took 1:28. Before picking a looser
default by analogy with the modal multiplexer's 1e-9, measure what THIS device
actually needs.

METHOD (costs ONE run, not one run per tolerance). `FDTDSolver.solve(
snapshot_every=, snapshot_path=)` decomposes the monitors every
`SNAP_EVERY` meep time units mid-run and appends {t, energy, energy_ratio, T}
to a CSV. Because DFT monitors ACCUMULATE, the sample at time t equals exactly
what a run stopped at t would have returned -- so one run sweeps the whole
transmission-vs-decay_tol curve. `energy_ratio` IS the decay_tol that would
have stopped the run at that sample, which makes the readout reference-free:
find where per-step |dT| drops below the precision the conclusions need, and
that is the required tolerance, measured instead of inherited.

Snapshotting is observation only -- it does not touch the fields, the stopping
condition, or the returned SMatrix, and is deliberately NOT in the cache
fingerprint, so this run correctly shares a cache entry with a plain one.

decay_tol is left TIGHT (1e-11) on purpose: the point is to run PAST the
convergence knee so the curve covers every looser tolerance. The CSV -- not
the SMatrix -- is this script's deliverable, and it is flushed per sample, so
a wall-clock kill still leaves a usable curve.

Read it with: plot_mmi1x2_tolcurve.py (safe to run while this is still going).

Config pinned to mmi_1x2/python/mmi_1x2_fdtd.py -- see that script for the
provenance of every constant.

Run:
    srun --mpi=pmi2 python3 mmi_1x2_decaytol_convergence.py
"""
import os

import numpy as np

from beef import MMI, FDTDSolver, SOI_160

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                   "example_output", "mmi_1x2"))
os.makedirs(OUT, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                        "example_sim_cache"))

# --- pinned to mmi_1x2_fdtd.py ---------------------------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.26, 1.36, 11

# Sample interval. The 4x4 used 100.0 on an ~84 um device; this one is 6.5 um,
# an order of magnitude shorter, so a proportionally finer interval keeps the
# knee resolved rather than stepped over. Overridable.
SNAP_EVERY = float(os.environ.get("MMI1X2_SNAP_EVERY", "10.0"))

mmi = MMI(
    n_inputs=1, n_outputs=2,
    mm_length=4.5, mm_width=2.0, waveguide_width=0.35,
    taper_length=1.0, taper_width=0.8, waveguide_offset=1.0,
    waveguide_height=0.16,
)

# mmi_1x2_fdtd.py's solver, VERBATIM. The first attempt at this script
# (job 12636314) omitted the monitor geometry and died in 12 s on
# "output_monitor_centers is empty" -- these are not optional defaults on this
# device. resolved/mmi_1x2_validation.md notes the 2026-08-09 output-monitor-
# width bug explicitly EXCLUDED the 1x2 precisely because it passes its own
# explicit widths rather than falling back on the buggy default, so copying
# them across is required for this curve to describe the shipped device.
solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=5.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25,
    three_d=True,
    num_input_modes=1, num_output_modes=1,
    # +/- sim_width/4 = +/-1.25 um, each spanning the shared y=0 axis out to
    # the near PML edge.
    output_monitor_centers=[2.5 / 2, -2.5 / 2],
    output_monitor_widths=[2.5, 2.5],
    material_stack=SOI_160,
    decay_tol=1e-11,     # TIGHT deliberately -- see module docstring
)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
snap_path = os.path.join(OUT, "tolcurve_res40_top.csv")

print(f"1x2 full-device FDTD, res={RESOLUTION:.0f}, {N_WL} wavelengths "
      f"{WL_MIN}-{WL_MAX} um")
print(f"device length {mmi.total_length:.2f} um, snapshot_every={SNAP_EVERY}")
print(f"snapshots -> {snap_path}")
print("decay_tol=1e-11 deliberately tight: run past the knee so the curve "
      "covers every looser tolerance.\n")

smat = solver.solve(mmi, wavelengths, port="top", input_mode=1,
                    save_dir=save_dir,
                    snapshot_every=SNAP_EVERY, snapshot_path=snap_path)

# FDTDSolver.solve() labels its OUTPUT ports "out_0", "out_1", ... by
# output_monitor_centers order -- NOT "out_top"/"out_bot", which is
# mmi_1x2_hybrid.py's (HybridSolver) convention. The `port=` argument above is
# the INPUT selector, unrelated to these keys. Mixing the two conventions is
# what crashed job 12636328 with KeyError('out_top', 1, 1) at this line, AFTER
# the full simulation had completed -- the CSV was already written, so only
# this print was lost. output_ys = [-0.5, +0.5], so index 0 is the bottom
# branch and index 1 the top. (resolved/mmi_1x2_validation.md gap 5 documents
# this same out_0/out_1 vs out_top/out_bot reconciliation.)
T_bot = smat.transmission(port="out_0")
T_top = smat.transmission(port="out_1")
i_c = len(wavelengths) // 2
print(f"\nfinal at {wavelengths[i_c]:.3f} um: T_top={T_top[i_c]:.6f} "
      f"T_bot={T_bot[i_c]:.6f} total={T_top[i_c] + T_bot[i_c]:.6f}")
print(f"CSV: {snap_path}")
