"""Prewarm ONE piece of the HCG Fabry-Perot's DIRECTIONAL COUPLER at a time.

Companion to _prewarm_hcgfp_arm_beta.py, and it exists for the same reason:
hcg_fabry_perot_hybrid.py needs several cached pieces, and running the script
itself to get them would solve them one after another in a single job. The
coupler's two S-bend excitations and its body beta calibration share no cache
key with each other, so they can run as CONCURRENT SLURM array tasks and
collapse the critical path from the sum of the three to the longest of them.

Selected with HCGFP_PREWARM_PIECE:

    sbend0   solve_section(input_idx=0)          -- "in" arm -> body modes
    sbend1   solve_section(input_idx=1)          -- "isolated" arm -> body modes
    body     section_mode_beta_details(N_MODES)  -- the body's own beta
    arm      section_mode_neff on the bare strip -- the UNCALIBRATED arm branch
    emebody  EMESolver.solve_section, no neff_override -- the real eigenmode
             solve that CALIBRATE_BODY=0 forces (with the override it is a free
             diagonal and needs no prewarm at all)

WHY res 80 IS EXPENSIVE HERE, stated so it is not a surprise: the coupler
solves are three-dimensional. At res 40 the S-bend cell is 11.75 x 8 x 2.5 um
= 15.0 Mcells and the pair plus the body calibration took 2 h 01 m on 24 cores
(job 12093321). Doubling the resolution is 8x the cells and 2x the timesteps,
so ~16x the work: order 780 core-hours for the same set at res 80.

THE CALLS BELOW MUST STAY BYTE-IDENTICAL to hcg_fabry_perot_hybrid.py's own
coupler block. That identity is the entire point -- a prewarm that differs by
one kwarg warms a key nothing ever looks up, and the expensive job is wasted
silently. It is CHECKABLE rather than merely asserted: run this file at the
resolution the hybrid already has cached and demand a hit,

    HCGFP_CACHE_ONLY=1 HCGFP_DC_RES=40 HCGFP_PREWARM_PIECE=sbend0 python <this>

which passes only if this file builds the identical call. Do that before
spending an allocation on res 80.

Run:
    HCGFP_PREWARM_PIECE=sbend0 python .../_prewarm_hcgfp_coupler.py
    HCGFP_DC_RES=40 HCGFP_PREWARM_PIECE=body python .../_prewarm_hcgfp_coupler.py
"""

import os
import time

import numpy as np

from beef import EMESolver, FDTDSolver, SOI_160
from beef.geometry import DirectionalCoupler, UniformCrossSection

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
save_dir = os.path.join(_ROOT, "example_sim_cache")

CACHE_ONLY = bool(int(os.environ.get("HCGFP_CACHE_ONLY", "0")))
PIECE = os.environ.get("HCGFP_PREWARM_PIECE", "").strip().lower()

# Band and resolutions default to hcg_fabry_perot_hybrid.py's own defaults.
WL_MIN = float(os.environ.get("HCGFP_WL_MIN", 1.530))
WL_MAX = float(os.environ.get("HCGFP_WL_MAX", 1.565))
N_WL = int(os.environ.get("HCGFP_N_WL", 36))
FDTD_RESOLUTION = float(os.environ.get("HCGFP_DC_RES", 80.0))
EME_RESOLUTION = float(os.environ.get("HCGFP_EME_RES", 80.0))
ARM_RESOLUTION = float(os.environ.get("HCGFP_ARM_RES", 80.0))

N_MODES = 2
COUPLING_LENGTH = 3.0
BODY_STUB = 0.25
WG_WIDTH = 0.5
ARM_CAL_SPAN = 6.0
MESH_WIDTH = 4.0
SIM_HEIGHT = 1.5
DPML = 0.5

wl = np.linspace(WL_MIN, WL_MAX, N_WL)

dc = DirectionalCoupler(
    waveguide_width=WG_WIDTH,
    gap=0.25,
    coupling_length=COUPLING_LENGTH,
    bend_length=10.0,
    waveguide_offset=4.0,
    bend_type="euler",
    waveguide_height=SOI_160.wg_height,
)
dc_fdtd = FDTDSolver(
    res=FDTD_RESOLUTION,
    sim_width=7.0,
    sim_height=SIM_HEIGHT,
    dpml=DPML,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    half_domain_aperture=True,
    material_stack=SOI_160,
    decay_tol=1e-11,
)
dc_eme = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=7.0, oxide_height=SIM_HEIGHT, res=EME_RESOLUTION,
    material_stack=SOI_160, mode_solver="meep",
)

half_body = dc.coupling_length / 2
half_total = dc.total_length / 2
boundaries = [-half_total, -half_body + BODY_STUB,
              +half_body - BODY_STUB, +half_total]

_iface_mon = dict(monitor_center=0.0, monitor_size=dc_fdtd.sim_width)

print("Prewarm: HCG-FP directional coupler")
print(f"  piece : {PIECE or '(none -- set HCGFP_PREWARM_PIECE)'}")
print(f"  band  : {WL_MIN:.4f}-{WL_MAX:.4f} um, {N_WL} pts")
print(f"  res   : S-bend/body FDTD {FDTD_RESOLUTION:g}, EME {EME_RESOLUTION:g}, "
      f"arm {ARM_RESOLUTION:g}")
print(f"  cache_only: {CACHE_ONLY}")

_t0 = time.time()

if PIECE in ("sbend0", "sbend1"):
    idx = 0 if PIECE == "sbend0" else 1
    out = dc_fdtd.solve_section(
        dc, wl, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=dc_fdtd.waveguide_length, stub_right=0.0,
        input_idx=idx, save_dir=save_dir, cache_only=CACHE_ONLY, **_iface_mon,
    )
    print(f"\nDone in {(time.time()-_t0)/60:.1f} min. |M|^2 at band centre: "
          f"{np.abs(out[N_WL // 2, :, 0])**2}")

elif PIECE == "body":
    cal = dc_fdtd.section_mode_beta_details(
        dc, wl, boundaries[1], boundaries[2], N_MODES,
        save_dir=save_dir, cache_only=CACHE_ONLY, **_iface_mon,
    )
    neff = cal.beta * wl[:, None] / (2 * np.pi)
    mpb = cal.seed_beta * wl[:, None] / (2 * np.pi)
    i = N_WL // 2
    print(f"\nDone in {(time.time()-_t0)/60:.1f} min. Even/odd supermodes at "
          f"{wl[i]:.4f} um:")
    print(f"  MPB  {mpb[i, 0]:.6f} {mpb[i, 1]:.6f}  d={mpb[i, 0]-mpb[i, 1]:.6f}")
    print(f"  FDTD {neff[i, 0]:.6f} {neff[i, 1]:.6f}  d={neff[i, 0]-neff[i, 1]:.6f}")

elif PIECE == "arm":
    arm_geom = UniformCrossSection(
        [(0.0, WG_WIDTH)], length=ARM_CAL_SPAN, waveguide_height=SOI_160.wg_height,
    )
    arm_solver = FDTDSolver(
        res=ARM_RESOLUTION, sim_width=MESH_WIDTH, sim_height=SIM_HEIGHT, dpml=DPML,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        num_input_modes=1, num_output_modes=1, output_monitor_centers=[0.0],
        material_stack=SOI_160, decay_tol=1e-11,
    )
    neff = arm_solver.section_mode_neff(
        arm_geom, wl, -1.0, +1.0, 1, save_dir=save_dir, cache_only=CACHE_ONLY,
    )[:, 0]
    print(f"\nDone in {(time.time()-_t0)/60:.1f} min. "
          f"MPB n_eff {neff[0]:.6f} -> {neff[-1]:.6f}")

elif PIECE == "emebody":
    out = dc_eme.solve_section(
        dc, wl, boundaries[1], boundaries[2],
        n_modes_in=N_MODES, n_modes_out=N_MODES,
        neff_override=None, save_dir=save_dir, cache_only=CACHE_ONLY,
    )
    print(f"\nDone in {(time.time()-_t0)/60:.1f} min. "
          f"|diag|^2 at band centre: {np.abs(np.diagonal(out[N_WL // 2]))**2}")

else:
    raise SystemExit(
        "set HCGFP_PREWARM_PIECE to one of: sbend0, sbend1, body, arm, emebody"
    )
