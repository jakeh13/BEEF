"""How expensive is `section_mode_beta()`'s seed step, really?
(running_jobs.md's 11592462 [0] entry / resolved/directional_coupler_validation.md.)

Job 11592462_0 (`beef_dc_examples` idx 0, `directional_coupler_hybrid.py`)
TIMED OUT at its full 5h. Reading its log line by line (all 512 lines) showed
neither of the script's own "2 FDTD solves" ever left a trace -- both
`Initializing structure...` blocks instead match `section_mode_neff()`'s own
documented shape exactly ("one structure, every wavelength"), not
`solve_section()`'s. `CALIBRATE_BODY = True` (hardcoded on) runs
`section_mode_beta()` first, whose *seed* step is one such call (N_MODES=5,
all 11 wavelengths); a second, separate single-wavelength
`section_mode_neff()` diagnostic call follows it later in the script. Both
blocks in the log **completed** -- but the 15 real, documented-as-expensive
`solve_section()` calibration runs that should run between them left zero
trace, while these two calls -- documented as "cheap... seconds, not
minutes" -- appear to have consumed nearly the whole 5h between them.

That is backwards from the docstring's own cost model. This script measures
`section_mode_neff()` directly, at the coupler's real production
configuration (N_MODES=5, the same cross-section, the same 3 x 7 x 2.5
computational cell), to find out whether it is really this expensive here --
before anyone resubmits the whole hybrid script a third time.

Cheap by design, same affordability argument as this session's other
`diagnose_*` scripts: no time-stepping, MPB eigensolves only. Reports
wall-clock per wavelength so a single login-node smoke run (1-2 wavelengths,
env `DC_DIAG_NWL=1`) can extrapolate the full 11-wavelength cost without
needing to run all 11 to find out.

Run (login node smoke check, per this workspace's convention):
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DC_DIAG_NWL=1 \\
        python diagnose_dc_calibration_seed_cost.py

Run (more wavelengths, still cheap -- eigensolve only):
    DC_DIAG_NWL=3 python diagnose_dc_calibration_seed_cost.py
"""
from __future__ import annotations

import os
import time

import numpy as np

from beef import DirectionalCoupler, FDTDSolver, SOI_220

N_WL_TOTAL = 11
N_WL = int(os.environ.get("DC_DIAG_NWL", "1"))
assert 1 <= N_WL <= N_WL_TOTAL, f"DC_DIAG_NWL must be in [1, {N_WL_TOTAL}]"

N_MODES = 5   # matches directional_coupler_hybrid.py's N_MODES exactly
BODY_STUB = 0.25
FDTD_RESOLUTION = 40.0
WL_MIN, WL_MAX = 1.50, 1.60

# Same reference geometry directional_coupler_hybrid.py itself uses: built at
# the middle of its 200-point coupling_length sweep, since the body cross-
# section (what section_mode_neff solves) doesn't depend on coupling_length
# at all -- only on gap/waveguide_width, both fixed below.
coupling_lengths = np.linspace(2 * BODY_STUB, 60.0, 200)
L_mid = float(coupling_lengths[len(coupling_lengths) // 2])

dc = DirectionalCoupler(
    waveguide_width=0.45, gap=0.25, coupling_length=L_mid, bend_length=9.0,
    waveguide_offset=3.0, bend_type="euler", waveguide_height=0.22,
)

fdtd_solver = FDTDSolver(
    res=FDTD_RESOLUTION, sim_width=6.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    half_domain_aperture=True, material_stack=SOI_220, decay_tol=1e-11,
)

half_body = dc.coupling_length / 2
_x0_full, _x1_full = -half_body + BODY_STUB, +half_body - BODY_STUB
# section_mode_beta() does NOT hand section_mode_neff() the full
# [_x0_full, _x1_full] span for its seed -- it narrows to a short window at
# the midpoint (seed_span, capped at 2.0 um) and section_mode_neff() then
# pads that by its own `pad` (default self.dpml) on each side when building
# the interior sim. Replicate that exactly, not the raw boundaries, or the
# geometry check rejects the cut (as an S-bend interior point) instead of
# measuring anything.
_mid = 0.5 * (_x0_full + _x1_full)
_span = _x1_full - _x0_full
_seed_span = min(2.0, _span - 2.0 * fdtd_solver.dpml)
x0, x1 = _mid - _seed_span / 2.0, _mid + _seed_span / 2.0
iface = dict(monitor_center=0.0, monitor_size=fdtd_solver.sim_width)

all_wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL_TOTAL)
# Same picks section_mode_beta()'s seed call would use, just truncated to
# N_WL points -- middle-out, so a 1-wavelength run tests the same band-centre
# point the script's own comparison-only diagnostic call uses (wavelengths[5]
# = 1.55 um), not an arbitrary edge point.
mid = N_WL_TOTAL // 2
order = sorted(range(N_WL_TOTAL), key=lambda i: abs(i - mid))[:N_WL]
wavelengths = all_wavelengths[sorted(order)]

print(f"[diag] N_MODES={N_MODES}, cross-section x in [{x0:+.3f}, {x1:+.3f}] "
      f"(coupling_length={L_mid:.2f} um mid-sweep reference)", flush=True)
print(f"[diag] timing {N_WL} of {N_WL_TOTAL} wavelengths: "
      f"{wavelengths.tolist()}", flush=True)

t0 = time.time()
per_wl_times = []
for i, wl in enumerate(wavelengths):
    t_wl0 = time.time()
    neff = fdtd_solver.section_mode_neff(
        dc, np.array([wl]), x0, x1, N_MODES, **iface)
    dt = time.time() - t_wl0
    per_wl_times.append(dt)
    print(f"[diag] wavelength {i+1}/{N_WL} (lambda={wl:.4f} um): "
          f"{dt:.1f} s, n_eff={neff[0].tolist()}", flush=True)

total = time.time() - t0
mean_per_wl = total / N_WL
print(f"\n[diag] {N_WL} wavelength(s), {total:.1f} s total, "
      f"{mean_per_wl:.1f} s/wavelength mean", flush=True)
print(f"[diag] extrapolated cost for all {N_WL_TOTAL} wavelengths "
      f"(section_mode_beta's seed step): {mean_per_wl * N_WL_TOTAL / 60:.1f} min "
      f"(one call) -- section_mode_beta also runs this seed step INSIDE "
      f"itself once, plus the separate single-wavelength comparison call "
      f"later in directional_coupler_hybrid.py costs ~{mean_per_wl:.1f} s more.",
      flush=True)
print(f"[diag] for comparison, the 15 real solve_section() calibration runs "
      f"section_mode_beta() also does (3 lengths x {N_MODES} excited modes) "
      f"are the ones documented as 'genuinely expensive' -- if the number "
      f"above is itself already minutes-to-hours, the docstring's 'seconds, "
      f"not minutes' claim for section_mode_neff does not hold at this "
      f"device's production N_MODES/domain size.", flush=True)
