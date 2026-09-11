"""Hybrid-path evaluation of the 4x4 hybrid at one (mm_width, y_gap) point.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-25.

WHY THE HYBRID PATH AND NOT PURE EME. The y_gap bias study (jobs 12352460-62
vs 12352475-78) measured the pure-EME error against full-device FDTD at two
y_gap values and it is NOT a constant offset:

    y_gap    EME max|dphi|   FDTD max|dphi|   EME error
    0.40        66.29 deg       22.44 deg      +43.86
    0.70        82.48 deg       20.46 deg      +62.02

The error MOVES by +18.2 deg over 0.3 um of y_gap while the true signal moves
only 2.0 deg -- the bias drifts ~9x faster than the quantity being optimized,
and the two methods RANK the two candidates in OPPOSITE order. Cause: EME
draws one input and one output branch per call and discards port crosstalk,
and y_gap IS the separation controlling that crosstalk. So a pure-EME
optimizer would be steered by the approximation. This path instead uses real
FDTD input sections (all 4 ports present simultaneously) cascaded with an EME
body -- the same construction mmi_4x4_90hybrid_body_sweep.py uses, which
check 4 showed tracks the full-device FDTD far more closely.

COST STRUCTURE, and the asymmetry this exploits. The body cross-section
depends ONLY on mm_width -- y_gap does not touch it. So the expensive meep
body mode solve is shared across every y_gap at a given width:

    * per mm_width : 1 body meep solve  (~77 min measured, job 11767748)
    * per y_gap    : FDTD input sections only (~45-75 min at 96 ranks,
                     extrapolated from the offset sweep's 2:10-3:02 at 24)

Driving the campaign width-major therefore costs several-fold less than
scattering points across mm_width. Everything goes through save_dir, so a
re-run of any point is a cache hit.

BODY BACKEND IS meep, NOT tidy3d -- deliberately. body_sweep.py's own comment:
"matches FDTDSolver's own eigenmode basis in this hybrid chain". The hybrid
cascade contracts an FDTD interface basis against the EME body basis, so the
two must be solved on the same engine; tidy3d here would reintroduce a
basis-mismatch error of exactly the kind this path exists to avoid.

Run:
    MMI4X4_W=6.0 MMI4X4_YGAP=0.70 python mmi_4x4_90hybrid_wg_eval.py
"""
import os
import time

import numpy as np

from beef import MMI, EMESolver, FDTDSolver, SMatrix, SOI_220
from beef.cache.store import CacheMissError
from mmi_4x4_90hybrid_ygap_grid import MM_LENGTH_GRID, LENGTH_MIN, LENGTH_MAX

# MMI4X4_POINT selects a grid point by index instead of naming W and y_gap
# separately -- the SLURM array passes its task id straight through, so the
# array's point numbering and the collector's come from the SAME enumeration
# (mmi_4x4_90hybrid_ygap_grid.py) rather than two hand-kept copies.
_POINT = os.environ.get("MMI4X4_POINT")
if _POINT not in (None, ""):
    from mmi_4x4_90hybrid_ygap_grid import POINTS
    MM_WIDTH, Y_GAP = POINTS[int(_POINT)]
    MM_WIDTH, Y_GAP = float(MM_WIDTH), float(Y_GAP)
else:
    MM_WIDTH = float(os.environ.get("MMI4X4_W", 6.0))
    Y_GAP = float(os.environ.get("MMI4X4_YGAP", 0.70))
# Phase A of the campaign: solve ONLY the body (which depends on mm_width and
# not on y_gap -- verified by direct fingerprint comparison, the tier-3
# neff_config is byte-identical across y_gap) so that every y_gap point at
# this width is later a cache hit. Without this, N points launched in
# parallel would each pay the ~77 min body solve independently.
BODY_ONLY = os.environ.get("MMI4X4_BODY_ONLY", "") not in ("", "0")
# Re-score a geometry that has already been simulated, without running anything:
# every FDTD section and the body length sweep must already be in the cache, and
# a miss raises rather than silently starting a 5 h solve. Added 2026-09-03 for
# the excess-loss audit -- the campaign's 33 points recorded only phase and CMRR,
# so loss had to be recovered from cache for geometries already paid for.
CACHE_ONLY = os.environ.get("MMI4X4_CACHE_ONLY", "") not in ("", "0")
RESOLUTION = float(os.environ.get("MMI4X4_EVAL_RES", 40.0))
DECAY_TOL = float(os.environ.get("MMI4X4_EVAL_TOL", 1e-4))
N_MODES = int(os.environ.get("MMI4X4_EVAL_NMODES", 16))
N_WL = int(os.environ.get("MMI4X4_EVAL_NWL", 3))
# 1.31 +- 0.02 um, per the user's spec for the optimization.
WL_MIN = float(os.environ.get("MMI4X4_WLMIN", 1.29))
WL_MAX = float(os.environ.get("MMI4X4_WLMAX", 1.33))

WAVEGUIDE_WIDTH = 0.35
# taper_length became a real search axis on 2026-09-01 (round 8). It was
# hardcoded at 1.0 um -- the MMI class default, Step 2's placeholder -- and the
# 2026-08-07/08 taper_sweep.py study closed it as null. That study was PURE EME,
# the path later disqualified for this device (no port crosstalk), and its
# "flat" answer was identical to 4 displayed digits across a 3x length change,
# which reads as solver insensitivity rather than physics. Under this
# evaluator's FDTD input sections (all 4 ports present) it is untested.
# Geometry: 0.35 -> taper_width over TAPER_LENGTH, so at the settled
# taper_width=1.065 the default is a 0.7 um width change in 1.0 um of length.
TAPER_LENGTH = float(os.environ.get("MMI4X4_TAPER_LENGTH", 1.0))
L_PI_REF = 108.9195          # measured at mm_width=6.0; L_pi/W^2 is ~constant
                             # (3.02-3.03 um^-1 over 3.5-16 um, Step 1c), so
                             # scale rather than re-measure per width.
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = WAVEGUIDE_OFFSET - Y_GAP
if TAPER_WIDTH <= WAVEGUIDE_WIDTH:
    raise SystemExit(
        f"INVERSE TAPER: y_gap={Y_GAP} at mm_width={MM_WIDTH} gives "
        f"taper_width={TAPER_WIDTH:.3f} <= waveguide_width={WAVEGUIDE_WIDTH}. "
        f"Valid y_gap < {WAVEGUIDE_OFFSET - WAVEGUIDE_WIDTH:.3f}.")

# L_pi scales as W^2; L_CENTRE (3/4 L_pi, the self-imaging condition the
# device was originally designed around) is used ONLY to build the geometry
# template below -- it fixes the taper/access-waveguide boundaries that the
# FDTD input-section solves are cached against (mm_length is part of that
# cache key, see the validation doc's mm_length/cache-fingerprint gotcha), so
# changing it would force 4 fresh ~1h FDTD re-solves for no benefit. The
# actual length SWEPT is the shared absolute MM_LENGTH_GRID (20-400 um @
# 50 nm, mmi_4x4_90hybrid_ygap_grid.py) -- independent of L_CENTRE, and free
# to widen because solve_section_length_sweep's cost never depends on how
# many lengths are requested (one mode solve total, see its own docstring).
L_PI = L_PI_REF * (MM_WIDTH / 6.0) ** 2
L_CENTRE = 0.75 * L_PI

DRIVEN = [0, 2]
OUTPUTS = [0, 1, 2, 3]
PORT_NAMES = ["out5", "out6", "out7", "out8"]
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                        "example_sim_cache"))
out_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                       "example_output", "mmi_4x4_90hybrid",
                                       "wg_ygap_opt"))
os.makedirs(out_dir, exist_ok=True)
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)


def make_mmi(mm_length):
    return MMI(n_inputs=4, n_outputs=4, mm_length=mm_length, mm_width=MM_WIDTH,
               waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
               taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
               waveguide_height=0.22, use_neighbor_aperture=True)


def make_solvers():
    fdtd_common = dict(
        res=RESOLUTION, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_220, decay_tol=DECAY_TOL,
        decay_method="energy", energy_decay_dt=20.0,
    )
    input_solver = FDTDSolver(**fdtd_common)
    body_solver = EMESolver(
        num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
        oxide_width=fdtd_common["sim_width"], oxide_height=fdtd_common["sim_height"],
        res=RESOLUTION, material_stack=SOI_220, mode_solver="meep",
    )
    return input_solver, body_solver


def boundaries_for(mmi):
    half_body = mmi.mm_length / 2
    half_total = mmi.total_length / 2
    body_stub = 0.25
    return [-half_total, -half_body + body_stub,
            +half_body - body_stub, +half_total], body_stub


def metrics_at(M_in, M_body_one):
    """Paper metrics for one body length. Same cascade the hybrid script uses."""
    M_out = [Mk.transpose(0, 2, 1) for Mk in M_in]
    smat = {i: SMatrix(wavelengths=wavelengths) for i in DRIVEN}
    for i in DRIVEN:
        tmp = np.einsum("Lwab,wbk->Lwak", M_body_one, M_in[i])
        for j in OUTPUTS:
            c = np.einsum("wab,Lwbk->Lwak", M_out[j], tmp)
            smat[i].set(c[0, :, 0, 0], port=PORT_NAMES[j])
    cm = [smat[i].cmrr_db(a, b) for i in DRIVEN for a, b in
          (("out5", "out8"), ("out6", "out7"))]
    # A11's corrected self/other order -- port3 as self, port1 as other.
    phi = {j: np.asarray(smat[2].cross_excitation_phase_diff_deg(
        smat[0], port=PORT_NAMES[j])) for j in OUTPUTS}
    dphi = {}
    for j in (1, 2, 3):
        raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
        dphi[j] = ((raw + 180.0) % 360.0) - 180.0
    max_dphi = float(max(np.abs(dphi[j]).max() for j in (1, 2, 3)))
    min_cmrr = float(min(np.asarray(c).min() for c in cm))
    # Excess loss and imbalance are RECORDED, not optimised -- the objective is
    # deliberately unchanged so every ranking in the validation doc stays
    # comparable. They are here because the campaign drove mm_width (and with
    # it mm_length, 91 -> 141 um) upward for 33 points while measuring only
    # phase and CMRR, which cannot see a device getting lossier as it grows.
    max_xl = float(max(np.asarray(smat[i].excess_loss_db(PORT_NAMES)).max()
                       for i in DRIVEN))
    max_imb = float(max(np.asarray(smat[i].imbalance_db(PORT_NAMES)).max()
                        for i in DRIVEN))
    return max_dphi, min_cmrr, max_xl, max_imb


def objective(max_dphi, min_cmrr, w_phase=0.5, w_cmrr=0.5):
    """Weighted, normalised, lower-is-better; 0 = both specs met.

    Each term scaled by its own spec so they are commensurate. CMRR is
    one-sided: exceeding 20 dBe earns no credit, it just stops penalising.
    """
    return (w_phase * (max_dphi / 5.0)
            + w_cmrr * max(0.0, (20.0 - min_cmrr) / 20.0))


def main():
    print(f"WGEVAL: mm_width={MM_WIDTH:.2f} y_gap={Y_GAP:.2f} -> "
          f"offset={WAVEGUIDE_OFFSET:.3f} taper_width={TAPER_WIDTH:.3f} "
          f"taper_length={TAPER_LENGTH:.2f}", flush=True)
    print(f"WGEVAL: L_pi={L_PI:.3f} centre={L_CENTRE:.3f} "
          f"length sweep {LENGTH_MIN}-{LENGTH_MAX} um, "
          f"{len(MM_LENGTH_GRID)} lengths, {N_WL} wl {WL_MIN}-{WL_MAX}", flush=True)
    print(f"WGEVAL: res={RESOLUTION} N_MODES={N_MODES} decay_tol={DECAY_TOL:.0e} "
          f"energy-decay", flush=True)

    input_solver, body_solver = make_solvers()
    mmi = make_mmi(L_CENTRE)
    boundaries, body_stub = boundaries_for(mmi)

    if BODY_ONLY:
        t1 = time.time()
        body_lengths = MM_LENGTH_GRID - 2 * body_stub
        body_solver.solve_section_length_sweep(
            mmi, wavelengths, boundaries[1], boundaries[2], body_lengths,
            n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir)
        print(f"WGEVAL: BODY_ONLY body solve + {len(body_lengths)}-length "
              f"sweep in {(time.time()-t1)/60:.1f} min", flush=True)
        return

    t0 = time.time()
    M_in = []
    for k in range(4):
        tk = time.time()
        M_in.append(input_solver.solve_section(
            mmi, wavelengths, boundaries[0], boundaries[1],
            n_modes_in=1, n_modes_out=N_MODES,
            stub_left=input_solver.waveguide_length, stub_right=0.0,
            input_idx=k, save_dir=save_dir, cache_only=CACHE_ONLY))
        print(f"WGEVAL:   input_idx={k} {time.time()-tk:.0f}s", flush=True)
    t_in = time.time() - t0
    print(f"WGEVAL: 4 FDTD input sections in {t_in/60:.1f} min", flush=True)

    t1 = time.time()
    body_lengths = MM_LENGTH_GRID - 2 * body_stub
    M_body_grid = body_solver.solve_section_length_sweep(
        mmi, wavelengths, boundaries[1], boundaries[2], body_lengths,
        n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
        cache_only=CACHE_ONLY)
    print(f"WGEVAL: body solve + {len(body_lengths)}-length sweep in "
          f"{(time.time()-t1)/60:.1f} min", flush=True)

    best = None
    for li in range(len(MM_LENGTH_GRID)):
        md, mc, xl, imb = metrics_at(M_in, M_body_grid[li:li + 1])
        obj = objective(md, mc)
        if best is None or obj < best[0]:
            best = (obj, float(MM_LENGTH_GRID[li]), md, mc, xl, imb)
    obj, best_mm_length, md, mc, xl, imb = best
    dl = best_mm_length - L_CENTRE  # informational only, vs the old design centre
    print(f"\nWGEVAL: BEST mm_length={best_mm_length:.2f} um "
          f"(dL vs 3/4-Lpi centre={dl:+.2f})", flush=True)
    print(f"WGEVAL:   max|dphi|={md:.3f} deg  min CMRR={mc:.3f} dBe", flush=True)
    print(f"WGEVAL:   max excess loss={xl:.3f} dB  max imbalance={imb:.3f} dB",
          flush=True)
    print(f"WGEVAL:   objective={obj:.4f}", flush=True)
    print(f"WGEVAL: total {(time.time()-t0)/60:.1f} min", flush=True)

    # The historical name is kept EXACTLY when taper_length is at its default,
    # so the 28 existing points and mmi_4x4_90hybrid_bo.py's own filename
    # parser see no change; any other taper_length gets a _tl tag, which that
    # parser now rejects rather than mis-reading as a default-taper point.
    fn = f"eval_W{MM_WIDTH:.2f}_y{Y_GAP:.2f}"
    if abs(TAPER_LENGTH - 1.0) > 1e-9:
        fn += f"_tl{TAPER_LENGTH:.2f}"
    # N_MODES belongs in the name for the same reason taper_length does: a
    # basis-convergence re-run is a DIFFERENT measurement of the same geometry,
    # not a replacement for it. Without this, round 9's N_MODES=20 confirmation
    # at (8.0, 0.80, tl=4.0) would have silently overwritten the N_MODES=16
    # result it exists to check against. Caught before the job's save step;
    # the 16-mode npz was also copied to *_nmodes16.npz.bak at that moment.
    if N_MODES != 16:
        fn += f"_n{N_MODES}"
    # Third instance of the same hazard (after _tl and _n): the wavelengths are
    # part of every FDTD cache key, so a different band is a DIFFERENT
    # measurement of the same geometry -- and without this tag the 1310+-10 nm
    # evaluation would overwrite the +-20 nm result it exists to be compared
    # against. Standing rule for this evaluator: every axis added here must also
    # appear in the output name.
    if abs(WL_MIN - 1.29) > 1e-9 or abs(WL_MAX - 1.33) > 1e-9:
        fn += f"_wl{WL_MIN:.3f}-{WL_MAX:.3f}"
    fn += ".npz"
    np.savez(os.path.join(out_dir, fn),
             mm_width=MM_WIDTH, y_gap=Y_GAP, taper_width=TAPER_WIDTH,
             taper_length=TAPER_LENGTH, n_modes=N_MODES,
             wl_min=WL_MIN, wl_max=WL_MAX, n_wl=N_WL,
             best_delta_L=dl, best_mm_length=best_mm_length,
             max_dphi=md, min_cmrr=mc, objective=obj,
             max_excess_loss_db=xl, max_imbalance_db=imb,
             wavelengths=wavelengths)
    print(f"WGEVAL: saved -> {fn}", flush=True)


if __name__ == "__main__":
    main()
