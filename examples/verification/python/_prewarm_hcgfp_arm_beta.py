"""Prewarm ONLY the HCG Fabry-Perot's connecting-waveguide beta calibration.

Why this exists as its own file rather than "just run the sweep script":
hcg_fabry_perot_length_sweep.py needs FOUR cached pieces, and as of
2026-08-21 three of them are being produced by jobs already in the queue
(12093330 = HCG mirror at res 80; 12093321 = coupler S-bends + body beta).
Running that script with HCGFP_CACHE_ONLY=0 before those land would not wait
for them -- it would SIMULATE them, duplicating two multi-hour jobs. This
file solves the one genuinely missing piece and nothing else.

The arm is a bare WG_WIDTH strip. It shares no cache key with the mirror or
the coupler, so it is safe to run in parallel with both.

The call below must stay BYTE-IDENTICAL to hcg_fabry_perot_length_sweep.py's
own arm block (and to hcg_fabry_perot_hybrid.py's, which differs only in
res/band) -- that identity is the entire point. If you change one, change
all of them, or this prewarm silently warms a key nothing looks up.

Run (defaults are the 2026-08-21 C-band re-run pins):
    python .../_prewarm_hcgfp_arm_beta.py
    HCGFP_ARM_RES=60 python .../_prewarm_hcgfp_arm_beta.py
    HCGFP_CACHE_ONLY=1 python .../_prewarm_hcgfp_arm_beta.py   # probe only
"""

import os
import time

import numpy as np

from beef import FDTDSolver, SOI_160
from beef.geometry import UniformCrossSection

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
save_dir = os.path.join(_ROOT, "example_sim_cache")

CACHE_ONLY = bool(int(os.environ.get("HCGFP_CACHE_ONLY", "0")))

WL_MIN = float(os.environ.get("HCGFP_WL_MIN", 1.530))
WL_MAX = float(os.environ.get("HCGFP_WL_MAX", 1.565))
N_WL = int(os.environ.get("HCGFP_N_WL", 36))
ARM_RESOLUTION = float(os.environ.get("HCGFP_ARM_RES", 80.0))

WG_WIDTH = 0.5
ARM_CAL_SPAN = 6.0
MESH_WIDTH = 4.0
SIM_HEIGHT = 1.5
DPML = 0.5

wl = np.linspace(WL_MIN, WL_MAX, N_WL)

print(f"Prewarm: HCG-FP connecting-waveguide beta calibration")
print(f"  band  : {WL_MIN:.4f}-{WL_MAX:.4f} um, {N_WL} pts")
print(f"  res   : {ARM_RESOLUTION:g}")
print(f"  strip : {WG_WIDTH:.2f} um wide, {ARM_CAL_SPAN:.1f} um calibration span")
print(f"  cache_only: {CACHE_ONLY}")

arm_geom = UniformCrossSection(
    [(0.0, WG_WIDTH)], length=ARM_CAL_SPAN, waveguide_height=SOI_160.wg_height,
)
arm_solver = FDTDSolver(
    res=ARM_RESOLUTION, sim_width=MESH_WIDTH, sim_height=SIM_HEIGHT, dpml=DPML,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    num_input_modes=1, num_output_modes=1, output_monitor_centers=[0.0],
    material_stack=SOI_160, decay_tol=1e-11,
)

_t0 = time.time()
cal = arm_solver.section_mode_beta_details(
    arm_geom, wl, -ARM_CAL_SPAN / 2.0, +ARM_CAL_SPAN / 2.0, 1,
    save_dir=save_dir, cache_only=CACHE_ONLY,
)
neff = cal.beta[:, 0] * wl / (2 * np.pi)
neff_mpb = cal.seed_beta[:, 0] * wl / (2 * np.pi)
_i = int(np.argmin(np.abs(wl - 0.5 * (WL_MIN + WL_MAX))))

print(f"\nDone in {(time.time()-_t0)/60:.1f} min.")
print(f"  n_eff at {wl[_i]:.4f} um: MPB {neff_mpb[_i]:.6f} -> FDTD {neff[_i]:.6f} "
      f"({neff[_i]/neff_mpb[_i]-1.0:+.4%})")
print(f"  n_eff across band: {neff[0]:.6f} -> {neff[-1]:.6f}")
print("\nCache is now warm for hcg_fabry_perot_length_sweep.py at this "
      "band/resolution. Verify with:")
print("  HCGFP_CACHE_ONLY=1 python .../hcg_fabry_perot_length_sweep.py")
