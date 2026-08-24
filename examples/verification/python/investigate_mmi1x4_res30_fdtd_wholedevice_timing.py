"""TIMING PROBE, not a full run -- times the whole-device FDTD domain's own
decay to decay_tol=1e-11 at the REAL geometry/resolution, before committing
mmi_1x4_fdtd.py's real production run blind.

N_WL=1 (not production's 11) -- this does NOT reduce decay cost (decay is
governed by the structure/PML, not by how many wavelengths' worth of DFT
accumulation ride along for free), but it does skip the comparatively cheap
post-decay decomposition step, so this probe's own wall-clock time is (to
good approximation) purely the decay time -- exactly the number this probe
exists to measure. There is no way to get this number more cheaply: decay
time at the real tolerance is the unknown itself, not a quantity a loose-
tolerance extrapolation can be trusted for on this device (this platform's
own decay behavior has already surprised twice -- input-section
res=30 needing ~1.4x MORE steps than expected, not fewer).

Real cost is genuinely unknown going in: could be as low as ~52 hrs
(if decay step count is similar to the input-section-alone probe's own
23,388 steps) or as high as ~390 hrs (if steps scale with the ~7.5x longer
propagation path) -- see mmi_1x4_fdtd.py's own module docstring for the
domain-size arithmetic behind that range. Bounded by this probe's own
walltime; if it hasn't decayed by then, that IS the answer (check 4 may
simply not be tractable at this device's real cost, not a failure of this
script).

Run:
    python .../investigate_mmi1x4_res30_fdtd_wholedevice_timing.py
"""
import os
import time
import numpy as np

from beef import MMI, FDTDSolver, SIN_300

RESOLUTION = 30.0
WL_CENTER = 1.55
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH, WAVEGUIDE_OFFSET, MM_WIDTH = 1.1, 2.2, 5.0, 3.1, 12.0
MM_LENGTH = 39.90
DECAY_TOL = 1e-11
_MONITOR_WIDTH = 2.8

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30, use_neighbor_aperture=True,
)
print(f"MMI total length: {mmi.total_length:.2f} um (whole-device domain, "
      f"x={mmi.total_length + 1.0:.2f} um incl. 2x dpml)", flush=True)

FDTD_SIM_WIDTH = MM_WIDTH + 4.0   # CHANGED 2026-08-19 from +8.0 -- see mmi_1x4_fdtd.py's own comment
FDTD_SIM_HEIGHT = 4.0             # CHANGED 2026-08-19 from SIN_300.sz (5.0)

solver = FDTDSolver(
    res=RESOLUTION, sim_width=FDTD_SIM_WIDTH, sim_height=FDTD_SIM_HEIGHT,
    dpml=0.5, source_margin=0.0, monitor_margin=0.25, three_d=True,
    num_input_modes=1, num_output_modes=1,
    output_monitor_centers=list(mmi.output_ys),
    output_monitor_widths=[_MONITOR_WIDTH] * 4,
    material_stack=SIN_300, decay_tol=DECAY_TOL,
)

wavelengths = np.array([WL_CENTER])

print(f"WHOLE-DEVICE FDTD TIMING PROBE, res={RESOLUTION}, decay_tol={DECAY_TOL}, "
      f"N_WL=1 (center only)", flush=True)
t0 = time.time()
smat = solver.solve(mmi, wavelengths, port="top", input_mode=1, save_dir=None)
elapsed = time.time() - t0
print(f"\nWHOLE-DEVICE PROBE DONE: {elapsed:.1f} s ({elapsed/3600:.3f} hrs)", flush=True)
T = [float(np.abs(smat.transmission(port=f"out_{k}", in_mode=1, out_mode=1))[0]) for k in range(4)]
print(f"T per output at center wavelength: {T}", flush=True)
print(f"Total T: {sum(T):.4f}", flush=True)
