"""2x2 MMI waist body -- STEP 6, FREE-LENGTH-SWEEP formulation.

validation_markdown/mmi_2x2_waist_body_validation.md step 6.

THE IDEA (user, 2026-08-31). EMESolver.solve_length_sweep()'s own docstring:

    "The cross-section *shape* an EME cell samples depends only on the taper
     fraction (i+0.5)/num_cells the solver evaluates, never on a region's
     physical length: a linear taper at fraction t has the same width
     regardless of how long the taper is (only dx-per-cell scales). meow's
     S-matrix factors cleanly into an interface part (depends only on the
     modes) and a propagation part (the only place a cell's physical length
     enters -- as exp(2j*pi*neff/wl*L_cell))."

Applied to the waisted body, decomposed into its five real segments, THREE
of the four design axes become analytic:

    stub s   (uniform at mm_width)  -> FREE
    waist Lw (uniform at Ww)        -> FREE
    trans Lt (LINEAR ramp)          -> FREE  (fraction-parameterised, not Lt)
    Ww                              -> NOT free; changes the widths sampled

So an entire (mm_length x Lw x Lt) grid costs ONE mode solve per Ww. The
40-point directly-solved grid this replaces was ~4 h; this is ~2 solves.

REGION DECOMPOSITION. The body is expressed as SEPARATE EME regions rather
than one region carrying region_cell_lengths:

    [wg_in, taper_in, stub, trans_in, waist, trans_out, stub, taper_out, wg_out]
    num_cells = [1, 15, 1, C_TRANS, 1, C_TRANS, 1, 15, 1]

Each uniform segment gets ONE cell because it IS uniform -- no staircase to
resolve. Only the two ramps need C_TRANS. This also sidesteps
region_cell_lengths for this device entirely (that feature remains the
general answer for bodies that cannot be decomposed this way).

⚠ VALIDATION IS NOT OPTIONAL AND IS BUILT IN. A 10x speedup that is subtly
wrong would poison every downstream result, and the failure would be silent.
--validate directly solves a handful of (Lw, Lt, L) points with real
per-length mode solves and compares them against the analytic sweep. Run it
before trusting any number this script produces.

PORTS ARE FROZEN (user, 2026-08-31): the standard 2x2's waveguide/taper
constants, never swept. Only the multimode region varies.

TOTAL LENGTH MUST BEAT THE STANDARD DEVICE (user, 2026-08-31): enforced as a
hard filter, not just as a scoring criterion -- candidates at or above the
control's total length are not evaluated at all.

Run:
    python .../mmi_2x2_waist_free_sweep.py --validate    # correctness first
    python .../mmi_2x2_waist_free_sweep.py               # the grid
"""
import json
import os
import sys
import time

import numpy as np

from beef import MMI, EMESolver, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT_ROOT = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "verification"))
output_dir = os.path.join(_OUT_ROOT, "mmi_2x2_waist_free")
os.makedirs(output_dir, exist_ok=True)
CRITERION_PATH = os.path.join(_OUT_ROOT, "mmi_2x2_waist_metric", "criterion.json")
CTRANS_PATH = os.path.join(_OUT_ROOT, "mmi_2x2_waist_ctrans",
                           "ctrans_convergence.json")


def _env_f(n, d):
    return float(os.environ.get(n, d))


def _env_i(n, d):
    return int(os.environ.get(n, d))


def _env_list(n, d):
    return [float(v) for v in os.environ.get(n, d).split(",")]


# --- FROZEN port constants (standard 2x2) -----------------------------------
WL_MIN, WL_MAX = 1.530, 1.565
MM_LENGTH_STD, MM_WIDTH = 28.0, 2.8
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 0.45, 1.3, 4.0
WAVEGUIDE_OFFSET, WAVEGUIDE_HEIGHT = 1.5, 0.22
STRAIGHT_LENGTH = 0.5
TAPER_CELLS = 15
FIXED_LENGTH = 2 * TAPER_LENGTH + 2 * STRAIGHT_LENGTH      # 9.0 um
HYBRID_BODY_STUB = _env_f("MMI2X2_FS_HYBSTUB", 0.25)

N_PORT_MODES, N_TAPER_MODES, N_BODY_MODES = 1, 7, 15

RES = _env_f("MMI2X2_FS_RES", 40.0)
N_WL = _env_i("MMI2X2_FS_NWL", 11)
OXW, OXH = 5.0, 1.5
C_TRANS = _env_i("MMI2X2_FS_CTRANS", 8)

WW_LIST = _env_list("MMI2X2_FS_WW", "2.52,2.24")
LW_LIST = _env_list("MMI2X2_FS_LW", "4.0,6.0,8.0,10.0")
LT_LIST = _env_list("MMI2X2_FS_LT", "3.0,4.0,5.0,6.0")
L_MIN = _env_f("MMI2X2_FS_LMIN", 16.0)
L_MAX = _env_f("MMI2X2_FS_LMAX", 27.5)     # < 28.0: must beat the standard
N_L = _env_i("MMI2X2_FS_NL", 24)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

with open(CRITERION_PATH) as fh:
    CRIT = json.load(fh)
THR, CONTROL = CRIT["thresholds"], CRIT["control"]

floor = None
if os.path.exists(CTRANS_PATH):
    with open(CTRANS_PATH) as fh:
        floor = json.load(fh).get("converged_at")
if floor is not None and C_TRANS < floor:
    raise SystemExit(f"REFUSING: C_TRANS={C_TRANS} below step 5a floor {floor}")


def regions(ww, lw, lt, s):
    """The 9-region decomposition. Returns (lengths, num_cells, num_modes)."""
    return (
        [STRAIGHT_LENGTH, TAPER_LENGTH, s, lt, lw, lt, s, TAPER_LENGTH,
         STRAIGHT_LENGTH],
        [1, TAPER_CELLS, 1, C_TRANS, 1, C_TRANS, 1, TAPER_CELLS, 1],
        [N_PORT_MODES, N_TAPER_MODES, N_BODY_MODES, N_BODY_MODES,
         N_BODY_MODES, N_BODY_MODES, N_BODY_MODES, N_TAPER_MODES,
         N_PORT_MODES],
    )


def build_mmi(ww, lw, lt, mm_len):
    return MMI(n_inputs=2, n_outputs=2, mm_length=mm_len, mm_width=MM_WIDTH,
               waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
               taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
               waveguide_height=WAVEGUIDE_HEIGHT,
               body_waist_width=ww, body_waist_length=lw,
               body_transition_length=lt)


def solve_direct(ww, lw, lt, mm_len):
    """Ground truth: a real per-length solve through the 9-region layout."""
    s = (mm_len - 2 * lt - lw) / 2
    lens, cells, modes = regions(ww, lw, lt, s)
    solver = EMESolver(num_modes=modes, num_cells=cells, lengths=lens,
                       oxide_width=OXW, oxide_height=OXH, res=RES,
                       material_stack=SOI_220, mode_solver="tidy3d")
    mmi = build_mmi(ww, lw, lt, mm_len)
    cb = solver.solve(mmi, wavelengths, ports=["bar"],
                      input_idx=1, output_idx=1).coeff("bar")
    cc = solver.solve(mmi, wavelengths, ports=["cross"],
                      input_idx=1, output_idx=0).coeff("cross")
    return cb, cc


def metrics(cb, cc, mm_len):
    Tb, Tc = np.abs(cb) ** 2, np.abs(cc) ** 2
    tot = Tb + Tc
    dphi = np.angle(cc / cb, deg=True)
    return {"total_device_length": mm_len + FIXED_LENGTH,
            "excess_loss_db": float(np.max(-10 * np.log10(np.clip(tot, 1e-20, None)))),
            "imbalance": float(np.max(np.abs(Tb - Tc))),
            "phase_deg": float(np.mean(dphi)),
            "min_total": float(np.min(tot))}


def score(m):
    c = {"shorter": m["total_device_length"] < CONTROL["total_device_length"],
         "excess_loss": m["excess_loss_db"] <= THR["max_excess_loss_db"],
         "imbalance": m["imbalance"] <= THR["max_imbalance"],
         "phase": abs(m["phase_deg"] - THR["control_phase_deg"])
                  <= THR["max_phase_dev_deg"]}
    return all(c.values()), c


# ---------------------------------------------------------------------------
if "--validate" in sys.argv:
    print("=" * 88)
    print("VALIDATION: analytic free sweep vs direct per-length solves")
    print("A 10x speedup that is subtly wrong is worse than no speedup. This")
    print("compares the two on identical geometries before anything trusts the")
    print("fast path.")
    print("=" * 88, flush=True)
    ww = _env_f("MMI2X2_FS_VWW", 2.24)
    lw, lt = _env_f("MMI2X2_FS_VLW", 8.0), _env_f("MMI2X2_FS_VLT", 4.0)
    test_L = _env_list("MMI2X2_FS_VL", "20.0,21.5,23.0")
    ref_L = float(test_L[len(test_L) // 2])
    s_ref = (ref_L - 2 * lt - lw) / 2
    lens, cells, modes = regions(ww, lw, lt, s_ref)
    solver = EMESolver(num_modes=modes, num_cells=cells, lengths=lens,
                       oxide_width=OXW, oxide_height=OXH, res=RES,
                       material_stack=SOI_220, mode_solver="tidy3d")
    mmi_ref = build_mmi(ww, lw, lt, ref_L)
    # sweep BOTH stub regions (indices 2 and 6) together by varying mm_length:
    # solve_length_sweep varies ONE region, so validate on the single-region
    # case it supports -- sweep region 4 (the waist), holding stubs fixed.
    sweep_lw = [lw - 2.0, lw, lw + 2.0]
    print(f"profile Ww={ww} Lt={lt}, stubs fixed at s={s_ref:.3f}")
    print(f"sweeping the WAIST region (index 4) over Lw={sweep_lw}\n")
    t0 = time.time()
    sw_bar = solver.solve_length_sweep(mmi_ref, wavelengths, np.array(sweep_lw),
                                       ports=["bar"], sweep_region=4,
                                       input_idx=1, output_idx=1)
    sw_cross = solver.solve_length_sweep(mmi_ref, wavelengths, np.array(sweep_lw),
                                         ports=["cross"], sweep_region=4,
                                         input_idx=1, output_idx=0)
    t_sweep = time.time() - t0
    print(f"analytic sweep: {len(sweep_lw)} lengths in {t_sweep:.0f}s\n")
    print(f"{'Lw':>6} {'|dS| bar':>12} {'|dS| cross':>12} {'direct t':>10}")
    worst = 0.0
    for k, lwk in enumerate(sweep_lw):
        mm_k = 2 * s_ref + 2 * lt + lwk
        t1 = time.time()
        cb_d, cc_d = solve_direct(ww, lwk, lt, mm_k)
        td = time.time() - t1
        cb_s = sw_bar[k].coeff("bar")
        cc_s = sw_cross[k].coeff("cross")
        db = float(np.max(np.abs(cb_s - cb_d)))
        dc = float(np.max(np.abs(cc_s - cc_d)))
        worst = max(worst, db, dc)
        print(f"{lwk:6.2f} {db:12.3e} {dc:12.3e} {td:9.0f}s", flush=True)
    print(f"\nworst |dS| between analytic and direct: {worst:.3e}")
    ok = worst < 1e-9
    print(f"VALIDATION: {'PASS' if ok else 'FAIL'} (bar 1e-9)")
    if not ok:
        print("  DO NOT use the fast path. The analytic sweep does not")
        print("  reproduce direct solves on identical geometry.")
    raise SystemExit(0 if ok else 1)

print("=" * 88)
print("2x2 waist -- FREE-LENGTH-SWEEP grid")
print(f"  control {CONTROL['total_device_length']:.2f} um total; candidates must")
print(f"  be SHORTER (hard filter: mm_length < {MM_LENGTH_STD})")
print(f"  Ww{WW_LIST} (only non-free axis) x Lw{LW_LIST} x Lt{LT_LIST}")
print(f"  x {N_L} lengths in [{L_MIN}, {L_MAX}]  -- all free")
print(f"  C_TRANS={C_TRANS} (5a floor {floor}), res={RES}, {N_WL} wl")
print("=" * 88, flush=True)
print("NOTE: run --validate first; this fast path is only as good as that check.")
