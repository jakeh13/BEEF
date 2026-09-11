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
import gc
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
    """Parse a list-valued env var. Accepts ',' OR ':' as the separator --
    `sbatch --export` splits on commas (it uses them as its own VAR=VAL
    separator) but not on colons, so a colon-separated list survives being
    passed that way. See the batch file's note."""
    raw = os.environ.get(n, d).replace(":", ",")
    return [float(v) for v in raw.split(",") if v.strip()]


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
# Env-overridable so the step-6 RE-SCORE can run the grid on a different mode
# basis than the original 7/15 (step 5b: the design point passes the phase
# criterion only at N_TAPER=7, so the grid has to be re-scored on a settled
# basis). These feed SHARD_FP below, so a re-score lands on its own shards
# instead of being mistaken for already-done rows.
N_TAPER_MODES = _env_i("MMI2X2_FS_NTAPER", N_TAPER_MODES)
N_BODY_MODES = _env_i("MMI2X2_FS_NBODY", N_BODY_MODES)
TAPER_CELLS = _env_i("MMI2X2_FS_TAPERCELLS", TAPER_CELLS)

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

# --- shardable, resumable row storage --------------------------------------
# The grid is ~384 points per Ww at ~69 s each uncached, i.e. 7.4 h -- one Ww
# does not fit an 8 h walltime with its basis solve. Two mechanisms fix that:
# rows stream to disk as they are produced, and a run skips any (Ww, Lw, Lt, L)
# already present in ANY shard file. So the grid can be split across jobs
# (MMI2X2_FS_WW / MMI2X2_FS_LW) and re-submitted after a timeout with no
# repeated work.
SHARD_TAG = os.environ.get("MMI2X2_FS_TAG", "main")
# The fingerprint is in the FILENAME so a cheap low-res smoke run can never be
# picked up as already-done grid rows: shards are only ever read back for the
# exact basis they were produced at.
#
# ⚠ THE MODE BASIS IS PART OF THE FINGERPRINT, added 2026-09-10, and leaving it
# out would have quietly corrupted the step-6 re-score. Step 5b showed the
# scoring depends strongly on N_TAPER (0.1217 dB over 5..13, and the design
# point passes the phase criterion ONLY at the N_TAPER=7 the grid used). The
# re-score therefore runs at a DIFFERENT taper basis than the original grid --
# but with a fingerprint of only (res, n_wl, C_TRANS) both would share one
# shard filename, so:
#   * `_shard_files()` globs that pattern, so `load_all_rows()` would MIX
#     N_TAPER=7 rows with re-scored rows indistinguishably, and
#   * the resume path skips any (Ww, Lw, Lt, L) already present in ANY shard,
#     so the re-score would skip all 668 points as "already done" and hand
#     back the ORIGINAL N_TAPER=7 numbers labelled as re-scored.
# A silently wrong answer, not a crash. Same defect class as the six
# filename-collision incidents this study has already recorded.
SHARD_FP = (f"res{RES:g}_nwl{N_WL}_ct{C_TRANS}"
            f"_nt{N_TAPER_MODES}_tc{TAPER_CELLS}_nb{N_BODY_MODES}")
ROWS_PATH = os.path.join(output_dir,
                         f"free_sweep_rows_{SHARD_FP}_{SHARD_TAG}.jsonl")
IFACE_CACHE = _env_i("MMI2X2_FS_IFCACHE", 1) == 1
SELFCHECK_PTS = [(1.00, 6.0, 10.0), (2.75, 4.0, 8.0)]


def key_of(ww, lw, lt, L):
    return (round(float(ww), 4), round(float(lw), 4), round(float(lt), 4),
            round(float(L), 4))


def _shard_files():
    import glob
    return sorted(glob.glob(os.path.join(
        output_dir, f"free_sweep_rows_{SHARD_FP}_*.jsonl")))


def load_all_rows():
    """Every row from every shard, de-duplicated on (Ww, Lw, Lt, L)."""
    out = {}
    for fp in _shard_files():
        with open(fp) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue     # a torn final line from a killed job
                out[key_of(r["Ww"], r["Lw"], r["Lt"], r["mm_length"])] = r
    return list(out.values())


def load_done():
    return {key_of(r["Ww"], r["Lw"], r["Lt"], r["mm_length"])
            for r in load_all_rows()}


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
    # ⚠ criterion.json pre-registers the phase check as
    #   "max |bar/cross phase - control phase| <= 2.0 deg"
    # but `phase_deg` here is the band MEAN, and score() compares mean-to-mean.
    # Excess loss and imbalance are both max-over-band, so phase is the one
    # component scored leniently: equal-and-opposite excursions at the two band
    # edges average away. These two extra fields report the max reading of the
    # SAME pre-registered text without changing what score() does -- the
    # scoring stays as-pre-registered and the discrepancy is reported, rather
    # than the criterion being quietly redefined after the sweep.
    ph_ctrl = THR["control_phase_deg"]
    return {"total_device_length": mm_len + FIXED_LENGTH,
            "excess_loss_db": float(np.max(-10 * np.log10(np.clip(tot, 1e-20, None)))),
            "imbalance": float(np.max(np.abs(Tb - Tc))),
            "phase_deg": float(np.mean(dphi)),
            "phase_max_dev_deg": float(np.max(np.abs(dphi - ph_ctrl))),
            "phase_spread_deg": float(np.ptp(dphi)),
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



# ---------------------------------------------------------------------------
# MULTI-AXIS free sweep. solve_length_sweep() varies ONE region, but the
# mechanism underneath is more general: its inner loop (beef/solvers/eme.py)
# reuses one mode set and passes `cell_lengths=` to mw.compute_s_matrix(). So
# an arbitrary cell-length list works with the same modes -- which means
# stub s, waist Lw AND transition Lt can all vary from a single solve.
#
# WHY THAT IS EXACT, per region:
#   stub  (uniform mm_width) -- one cell, any length
#   waist (uniform Ww)       -- one cell, any length
#   ramp  (linear W->Ww)     -- cell j samples the width at fraction
#                               (j+0.5)/C_TRANS along the ramp, which does not
#                               depend on Lt; only dx-per-cell scales
# Ww is the one axis that changes the widths sampled, so it alone needs a
# fresh solve. Validated against direct solves at 1e-13 (job 12624884).
# ---------------------------------------------------------------------------
REGION_STUB_L, REGION_RAMP_IN, REGION_WAIST, REGION_RAMP_OUT, REGION_STUB_R = 2, 3, 4, 5, 6


def _cell_lengths_for(s, lt, lw):
    """Flat per-cell length list for the 9-region layout at (s, Lt, Lw)."""
    return ([STRAIGHT_LENGTH]
            + [TAPER_LENGTH / TAPER_CELLS] * TAPER_CELLS
            + [s]
            + [lt / C_TRANS] * C_TRANS
            + [lw]
            + [lt / C_TRANS] * C_TRANS
            + [s]
            + [TAPER_LENGTH / TAPER_CELLS] * TAPER_CELLS
            + [STRAIGHT_LENGTH])


def build_basis(ww, lw_ref, lt_ref, s_ref):
    """One mode solve per port for this Ww. Returns {port: (cells, modes_per_wl)}."""
    import meow as mw
    lens, cells_n, modes_n = regions(ww, lw_ref, lt_ref, s_ref)
    solver = EMESolver(num_modes=modes_n, num_cells=cells_n, lengths=lens,
                       oxide_width=OXW, oxide_height=OXH, res=RES,
                       material_stack=SOI_220, mode_solver="tidy3d")
    mmi_ref = build_mmi(ww, lw_ref, lt_ref, 2 * s_ref + 2 * lt_ref + lw_ref)
    mpc = solver._modes_per_cell()
    out = {}
    for port, (side, iidx, oidx) in {"bar": ("bar", 1, 1),
                                     "cross": ("cross", 1, 0)}.items():
        cells = solver._build_cells(mmi_ref, side=side, wavelengths=wavelengths,
                                    input_idx=iidx, output_idx=oidx)
        per_wl = []
        for wl in wavelengths:
            env = mw.Environment(wl=float(wl), T=solver.temperature)
            css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
            per_wl.append([solver.compute_modes(cs, nm)
                           for cs, nm in zip(css, mpc)])
        out[port] = (cells, per_wl)
    return out


# --- EXACT interface-S-matrix memoisation ----------------------------------
# meow's compute_s_matrix_sax(modes, cells, cell_lengths, interfaces_fn=...)
# splits into compute_interface_s_matrices(modes) -- a pure function of the
# MODES -- and compute_propagation_s_matrices(..., cell_lengths=...), the only
# consumer of length. Across this grid the modes never change (that is the
# premise the whole free sweep rests on), so the interface matrices are
# recomputed identically at every one of the ~384 points per Ww. Memoising
# them is EXACT, not an approximation -- and it is self-checked at startup
# below rather than trusted.
_IFACE_CACHE = {}


def _make_interfaces_fn(key):
    """Memoised interfaces_fn for one (port, wl) modal basis."""
    # the exact import path job 12674313's probe validated
    from meow.eme.interface import compute_interface_s_matrices

    def fn(modes):
        hit = _IFACE_CACHE.get(key)
        if hit is None:
            hit = compute_interface_s_matrices(modes)
            _IFACE_CACHE[key] = hit
        return hit

    return fn


def eval_from_basis(basis, s, lt, lw, iface_cache=True):
    """S-coefficients at (s, Lt, Lw) from an already-solved basis -- no solve."""
    import meow as mw
    cl = _cell_lengths_for(s, lt, lw)
    coeffs = {}
    for port, (cells, per_wl) in basis.items():
        amp = np.zeros(len(wavelengths), dtype=complex)
        for i, modes in enumerate(per_wl):
            kw = ({"interfaces_fn": _make_interfaces_fn((port, i))}
                  if iface_cache else {})
            S, _ = mw.compute_s_matrix(modes, cells=cells, cell_lengths=cl,
                                       **kw)
            amp[i] = S[0, 1]
        coeffs[port] = amp
    return coeffs["bar"], coeffs["cross"]


def selfcheck_iface_cache(basis, pts):
    """Cached vs uncached on identical geometry. Aborts the job if it differs.

    Cheap (a handful of points) against a whole grid produced on a fast path,
    and it runs in the SAME process as the grid, so it certifies this job's
    own cache rather than a previous job's.
    """
    worst = 0.0
    for (sv, ltv, lwv) in pts:
        _IFACE_CACHE.clear()
        t0 = time.time()
        cb_c, cc_c = eval_from_basis(basis, sv, ltv, lwv, iface_cache=True)
        t_c = time.time() - t0
        t0 = time.time()
        cb_u, cc_u = eval_from_basis(basis, sv, ltv, lwv, iface_cache=False)
        t_u = time.time() - t0
        d = max(float(np.max(np.abs(cb_c - cb_u))),
                float(np.max(np.abs(cc_c - cc_u))))
        worst = max(worst, d)
        # the pass above was COLD -- it POPULATED the cache, so it pays the
        # full interface cost exactly like the uncached one. The speedup only
        # exists warm, which is the regime the grid loop runs in.
        t0 = time.time()
        eval_from_basis(basis, sv, ltv, lwv, iface_cache=True)
        t_w = time.time() - t0
        print(f"   selfcheck s={sv:.2f} Lt={ltv:.1f} Lw={lwv:.1f}: "
              f"uncached {t_u:.1f}s | cache-cold {t_c:.1f}s | "
              f"cache-warm {t_w:.1f}s ({t_u/max(t_w,1e-9):.1f}x)  "
              f"|dS|={d:.3e}", flush=True)
    print(f"   worst |dS| cached vs uncached: {worst:.3e}", flush=True)
    if worst > 1e-13:
        raise SystemExit("REFUSING: interface cache is not exact "
                         f"(|dS|={worst:.3e} > 1e-13). Set "
                         "MMI2X2_FS_IFCACHE=0 and re-run.")
    return worst



# ---------------------------------------------------------------------------
if "--validate-multi" in sys.argv:
    print("=" * 88)
    print("VALIDATION (MULTI-AXIS): eval_from_basis() vs direct solves")
    print("The single-region check (--validate) proves solve_length_sweep's")
    print("mechanism. This checks MY reimplementation, which varies THREE")
    print("regions at once through hand-built cells -- a different code path.")
    print("Trusting the first check for the second would be unjustified.")
    print("=" * 88, flush=True)
    ww = _env_f("MMI2X2_FS_VWW", 2.24)
    lt_ref, lw_ref, s_ref = 4.0, 8.0, 2.75
    t0 = time.time()
    basis = build_basis(ww, lw_ref, lt_ref, s_ref)
    print(f"basis (both ports) solved in {time.time()-t0:.0f}s\n", flush=True)
    # deliberately vary ALL THREE axes away from the reference, including Lt
    trials = [(s_ref, lt_ref, lw_ref),          # the reference itself
              (1.50, lt_ref, lw_ref),           # stub only
              (s_ref, lt_ref, 6.00),            # waist only
              (s_ref, 6.00, lw_ref),            # RAMP length only -- the
                                                # least obvious of the three
              (1.00, 5.00, 6.00)]               # all three at once
    print(f"{'s':>6} {'Lt':>6} {'Lw':>6} {'|dS| bar':>12} {'|dS| cross':>12} {'t':>7}")
    worst = 0.0
    for (sv, ltv, lwv) in trials:
        cb_a, cc_a = eval_from_basis(basis, sv, ltv, lwv)
        mm = 2 * sv + 2 * ltv + lwv
        t1 = time.time()
        cb_d, cc_d = solve_direct(ww, lwv, ltv, mm)
        db = float(np.max(np.abs(cb_a - cb_d)))
        dc = float(np.max(np.abs(cc_a - cc_d)))
        worst = max(worst, db, dc)
        print(f"{sv:6.2f} {ltv:6.2f} {lwv:6.2f} {db:12.3e} {dc:12.3e} "
              f"{time.time()-t1:6.0f}s", flush=True)
    print(f"\nworst |dS|: {worst:.3e}")
    ok = worst < 1e-9
    print(f"MULTI-AXIS VALIDATION: {'PASS' if ok else 'FAIL'} (bar 1e-9)")
    if not ok:
        print("  DO NOT use the multi-axis fast path. Fall back to direct solves.")
    raise SystemExit(0 if ok else 1)

# ---------------------------------------------------------------------------
# --taper-convergence: step 5b PART 2 -- push N_TAPER past 13 and sweep the
# taper CELL COUNT with it.
#
# WHY PAIRS AND NOT TWO INDEPENDENT SWEEPS. Everywhere else in this study the
# rule is "sweep knobs separately or a joint curve hides which one is
# unconverged" (check 3). Here the hypothesis IS an interaction: 5b part 1 found
# N_TAPER still moving at 13 while TAPER_CELLS sat fixed at 15, so pushing more
# modes through the same 15 cells may itself be the limiter. Independent sweeps
# cannot test that. This mode therefore takes EXPLICIT (N_TAPER, TAPER_CELLS)
# pairs so the two edges and the diagonal are all measured:
#   edge B'  : N_TAPER 15..21 at TAPER_CELLS=15   (does the mode axis plateau?)
#   edge C   : TAPER_CELLS 10..30 at N_TAPER=13   (does the cell axis matter?)
#   diagonal : both raised together                (is it the interaction?)
#
# 5b part 1 established what this must beat: N_BODY is converged at 15 (span
# 0.0053 dB over 12..27), while N_TAPER spans 0.1217 dB over 5..13 and the
# design point passes the phase criterion ONLY at the N_TAPER=7 that step 6
# happened to use. Until this axis plateaus there is no defensible basis to
# re-score step 6 on.
#
# Pairs are given as "ntaper:ncells,ntaper:ncells,..." -- note ':' inside a
# pair and ',' between pairs, and see _env_list's note about sbatch --export
# eating commas.
# ---------------------------------------------------------------------------
if "--taper-convergence" in sys.argv:
    print("=" * 88)
    print("STEP 5b PART 2 -- N_TAPER past 13, and TAPER_CELLS with it")
    print("=" * 88, flush=True)
    ww = _env_f("MMI2X2_TC_WW", 2.24)
    lw = _env_f("MMI2X2_TC_LW", 8.0)
    lt = _env_f("MMI2X2_TC_LT", 6.0)
    L = _env_f("MMI2X2_TC_L", 20.5)
    spec = os.environ.get(
        "MMI2X2_TC_PAIRS",
        # edge B' (cells fixed at 15)      edge C (modes fixed at 13)
        "15:15,17:15,19:15,21:15,"
        "13:10,13:20,13:25,13:30,"
        # diagonal: both raised together
        "19:25,21:30")
    pairs = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        a, b = chunk.split(":")
        pairs.append((int(a), int(b)))
    s_stub = (L - 2 * lt - lw) / 2
    FDTD = {"excess_loss_db": 0.2036, "imbalance": 0.0976,
            "phase_deg": -91.8486}
    # measured in 5b part 1 (job 12780420), TAPER_CELLS=15 throughout
    PART1 = {(5, 15): 0.4676, (7, 15): 0.3910, (9, 15): 0.3658,
             (11, 15): 0.3459, (13, 15): 0.3536}
    print(f"design point Ww={ww} Lw={lw} Lt={lt} L={L} (stub {s_stub:.3f})")
    print(f"{len(pairs)} (N_TAPER, TAPER_CELLS) pairs: {pairs}")
    print(f"N_BODY held at {N_BODY_MODES} (5b part 1: converged at 15)")
    print(f"FDTD reference (check 4): IL {FDTD['excess_loss_db']:.4f} dB | "
          f"phase {FDTD['phase_deg']:.3f} deg")
    print("5b part 1, TAPER_CELLS=15: " +
          ", ".join(f"{k[0]}->{v:.4f}" for k, v in sorted(PART1.items())))
    print(flush=True)

    # Filename carries the GEOMETRY and the PAIRS, not just res. This mode is
    # meant to be run as a job array (one pair per task, so the heaviest pair
    # still fits the -ece 5h cap), and a res-only name would have every task
    # writing the same file -- the clobbering failure this study has already
    # hit six times, most recently with a geometry-less check4 name. One task
    # per pair therefore lands on its own shard; assemble by globbing.
    _pairs_tag = "_".join(f"{a}x{b}" for a, b in pairs)
    tpath = os.path.join(
        output_dir,
        f"taper_convergence_res{RES:g}_Ww{ww:g}_Lw{lw:g}_Lt{lt:g}_L{L:g}"
        f"_nb{N_BODY_MODES}_nwl{N_WL}_pairs{_pairs_tag}.json")
    out = []
    print(f"{'N_TAPER':>7} {'CELLS':>6} {'IL dB':>9} {'dIL vs FDTD':>11} "
          f"{'imb':>8} {'phase':>10} {'mean-dev':>9} {'basis':>7} {'crit':>5}")
    for (ntaper, ncells) in pairs:
        global_prev = (N_TAPER_MODES, TAPER_CELLS)
        globals()["N_TAPER_MODES"] = ntaper
        globals()["TAPER_CELLS"] = ncells
        _IFACE_CACHE.clear()
        gc.collect()
        t0 = time.time()
        b = build_basis(ww, lw, lt, max(HYBRID_BODY_STUB, 1.0))
        dt = time.time() - t0
        cb, cc = eval_from_basis(b, s_stub, lt, lw, iface_cache=IFACE_CACHE)
        m = metrics(cb, cc, float(L))
        passed, _ = score(m)
        mean_dev = abs(m["phase_deg"] - THR["control_phase_deg"])
        row = {"n_taper": ntaper, "taper_cells": ncells, "basis_seconds": dt,
               "passed": passed, "mean_dev_deg": mean_dev,
               "d_IL_vs_fdtd": m["excess_loss_db"] - FDTD["excess_loss_db"],
               **m}
        out.append(row)
        with open(tpath, "w") as fh:
            json.dump({"design": {"Ww": ww, "Lw": lw, "Lt": lt, "L": L},
                       "fdtd_reference": FDTD, "part1_taper_cells_15": {
                           f"{k[0]}": v for k, v in PART1.items()},
                       "res": RES, "n_wl": N_WL, "n_body": N_BODY_MODES,
                       "complete": False, "rows": out}, fh, indent=2)
        del b
        globals()["N_TAPER_MODES"], globals()["TAPER_CELLS"] = global_prev
        gc.collect()
        print(f"{ntaper:7d} {ncells:6d} {m['excess_loss_db']:9.4f} "
              f"{row['d_IL_vs_fdtd']:+11.4f} {m['imbalance']:8.4f} "
              f"{m['phase_deg']:10.3f} {mean_dev:9.3f} {dt:6.0f}s "
              f"{'PASS' if passed else 'FAIL':>5}", flush=True)

    with open(tpath, "w") as fh:
        json.dump({"design": {"Ww": ww, "Lw": lw, "Lt": lt, "L": L},
                   "fdtd_reference": FDTD,
                   "part1_taper_cells_15": {f"{k[0]}": v
                                            for k, v in PART1.items()},
                   "res": RES, "n_wl": N_WL, "n_body": N_BODY_MODES,
                   "complete": True, "rows": out}, fh, indent=2)
    print(f"\nwritten -> {tpath}")

    def show(label, rows_x, key):
        if len(rows_x) < 3:
            print(f"\n--- {label} --- only {len(rows_x)} points: NO VERDICT")
            return
        print(f"\n--- {label} ---")
        for r in rows_x:
            line = (f"  {key}={r[key]:>3}  IL {r['excess_loss_db']:.4f}  "
                    f"(vs FDTD {r['d_IL_vs_fdtd']:+.4f})")
            if r.get("from_part1"):
                # part 1 recorded IL only. Printing "phase nan ... FAIL" here
                # mislabelled N_TAPER=7 -- the ONE value that passed -- as a
                # failure, in the very table the conclusion is read off.
                line += "   [5b part 1: IL only, see that table for phase]"
            else:
                line += (f"  phase {r['phase_deg']:.3f}  "
                         f"mean-dev {r['mean_dev_deg']:.3f}"
                         f"  {'PASS' if r['passed'] else 'FAIL'}")
            print(line)
        span = (max(r["excess_loss_db"] for r in rows_x)
                - min(r["excess_loss_db"] for r in rows_x))
        print(f"  IL span: {span:.4f} dB "
              f"({'FLAT' if span < 0.01 else 'STILL MOVING'})")

    edgeB = sorted([r for r in out if r["taper_cells"] == 15],
                   key=lambda r: r["n_taper"])
    # fold in part 1's own TAPER_CELLS=15 points so the mode edge is whole
    for nt, il in sorted(PART1.items()):
        if not any(r["n_taper"] == nt[0] for r in edgeB):
            edgeB.append({"n_taper": nt[0], "taper_cells": 15,
                          "excess_loss_db": il, "from_part1": True,
                          "d_IL_vs_fdtd": il - FDTD["excess_loss_db"]})
    edgeB.sort(key=lambda r: r["n_taper"])
    show("edge B': N_TAPER at TAPER_CELLS=15 (incl. 5b part 1)",
         edgeB, "n_taper")
    show("edge C: TAPER_CELLS at N_TAPER=13",
         sorted([r for r in out if r["n_taper"] == 13],
                key=lambda r: r["taper_cells"]), "taper_cells")
    diag = [r for r in out if r["n_taper"] != 13 and r["taper_cells"] != 15]
    if diag:
        print("\n--- diagonal (both raised) ---")
        for r in diag:
            print(f"  N_TAPER={r['n_taper']} CELLS={r['taper_cells']}  "
                  f"IL {r['excess_loss_db']:.4f} (vs FDTD "
                  f"{r['d_IL_vs_fdtd']:+.4f})  mean-dev {r['mean_dev_deg']:.3f}"
                  f"  {'PASS' if r['passed'] else 'FAIL'}")
    print("\n⚠ A plateau at the largest value TRIED is not a ceiling.")
    print("⚠ Until an axis is FLAT, step 6 has no defensible basis to "
          "re-score on.")
    raise SystemExit(0)


# ---------------------------------------------------------------------------
# --modes-convergence: STEP 5b, the num_modes plateau for the waisted body.
#
# WHY THIS IS NOT A FORMALITY, and why it is worth running AFTER check 4.
# The entire step-6 grid ran at N_BODY_MODES=15 / N_TAPER_MODES=7 with no
# convergence evidence for this geometry (the straight device deliberately uses
# 14). Check 4 then measured EME over-estimating excess loss by **0.19-0.24 dB
# on BOTH geometries** vs a full-device FDTD reference -- a consistent,
# geometry-independent offset, which is exactly the signature of a basis too
# small to hold the power a brute-force solve keeps.
#
# TESTABLE PREDICTION: if truncation is the cause, raising the mode count
# lowers the EME excess loss TOWARD the FDTD value and then plateaus. If the
# loss is flat in mode count, truncation is refuted and the offset is something
# else. Both outcomes are informative, which is the point.
#
# The straight device's N_MODES=14 was chosen because band 14 sits below the
# cladding index -- i.e. the bands above are RADIATIVE, and radiation modes are
# precisely what an EME basis needs to represent power leaving the guide. So
# the body sweep deliberately runs well past 15.
#
# TWO INDEPENDENT 1-D SWEEPS, not a grid -- the same reason check 3 sweeps
# engines separately: a joint curve cannot tell you WHICH knob is unconverged.
#   A: N_TAPER fixed, N_BODY varied   (the 2.8 um faces and the 2.24 um waist)
#   B: N_BODY fixed, N_TAPER varied   (the 1.3->2.8 um port tapers)
#
# ⚠ The 4x4's warning, quoted in the 1x4 doc: a probe that hits its own request
# cap only ever measured a FLOOR, never found a ceiling. Push until the metric
# plateaus; do not stop at the largest value tried.
# ---------------------------------------------------------------------------
if "--modes-convergence" in sys.argv:
    print("=" * 88)
    print("STEP 5b -- num_modes convergence at the waisted design point")
    print("=" * 88, flush=True)
    ww = _env_f("MMI2X2_MC_WW", 2.24)
    lw = _env_f("MMI2X2_MC_LW", 8.0)
    lt = _env_f("MMI2X2_MC_LT", 6.0)
    L = _env_f("MMI2X2_MC_L", 20.5)
    body_list = [int(v) for v in _env_list("MMI2X2_MC_BODY", "12,15,18,21,24,27")]
    taper_list = [int(v) for v in _env_list("MMI2X2_MC_TAPER", "5,7,9,11,13")]
    body_ref = _env_i("MMI2X2_MC_BODYREF", 15)
    taper_ref = _env_i("MMI2X2_MC_TAPERREF", 7)
    s_stub = (L - 2 * lt - lw) / 2
    # check-4 full-device FDTD reference, job 12743476 (res 40, 51 wl -> the
    # criterion's own 11-point subset, decay_tol 1e-5, passive=True)
    FDTD = {"excess_loss_db": 0.2036, "imbalance": 0.0976,
            "phase_deg": -91.8486}
    print(f"design point Ww={ww} Lw={lw} Lt={lt} L={L} (stub {s_stub:.3f})")
    print(f"sweep A: N_TAPER={taper_ref} fixed, N_BODY {body_list}")
    print(f"sweep B: N_BODY={body_ref} fixed, N_TAPER {taper_list}")
    print(f"FDTD reference (check 4): IL {FDTD['excess_loss_db']:.4f} dB | "
          f"imb {FDTD['imbalance']:.4f} | phase {FDTD['phase_deg']:.3f} deg")
    print("PREDICTION UNDER TEST: if the EME-vs-FDTD loss gap is basis")
    print("truncation, IL falls toward the FDTD value and plateaus.\n",
          flush=True)

    mpath = os.path.join(output_dir, f"modes_convergence_res{RES:g}.json")
    out = []

    def run_point(nbody, ntaper, axis):
        # regions() reads these globals at call time, so rebinding them here
        # is what changes the basis the solve is built on.
        global N_BODY_MODES, N_TAPER_MODES
        prev_b, prev_t = N_BODY_MODES, N_TAPER_MODES
        N_BODY_MODES, N_TAPER_MODES = nbody, ntaper
        _IFACE_CACHE.clear()
        gc.collect()
        t0 = time.time()
        b = build_basis(ww, lw, lt, max(HYBRID_BODY_STUB, 1.0))
        dt = time.time() - t0
        cb, cc = eval_from_basis(b, s_stub, lt, lw, iface_cache=IFACE_CACHE)
        m = metrics(cb, cc, float(L))
        passed, _ = score(m)
        row = {"axis": axis, "n_body": nbody, "n_taper": ntaper,
               "basis_seconds": dt, "passed": passed,
               "d_IL_vs_fdtd": m["excess_loss_db"] - FDTD["excess_loss_db"],
               **m}
        out.append(row)
        # every point hits disk as it lands
        with open(mpath, "w") as fh:
            json.dump({"design": {"Ww": ww, "Lw": lw, "Lt": lt, "L": L},
                       "fdtd_reference": FDTD, "res": RES, "n_wl": N_WL,
                       "complete": False, "rows": out}, fh, indent=2)
        del b
        N_BODY_MODES, N_TAPER_MODES = prev_b, prev_t
        gc.collect()
        print(f"{axis:>1} {nbody:6d} {ntaper:7d} {m['excess_loss_db']:9.4f} "
              f"{row['d_IL_vs_fdtd']:+9.4f} {m['imbalance']:9.4f} "
              f"{m['phase_deg']:10.3f} {dt:7.0f}s "
              f"{'PASS' if passed else 'FAIL'}", flush=True)

    print(f"{'ax':>1} {'N_BODY':>6} {'N_TAPER':>7} {'IL dB':>9} "
          f"{'dIL vs FDTD':>9} {'imb':>9} {'phase':>10} {'basis':>8}")
    for nb in body_list:
        run_point(nb, taper_ref, "A")
    for nt in taper_list:
        if nt == taper_ref:
            continue          # already measured on axis A at body_ref
        run_point(body_ref, nt, "B")

    with open(mpath, "w") as fh:
        json.dump({"design": {"Ww": ww, "Lw": lw, "Lt": lt, "L": L},
                   "fdtd_reference": FDTD, "res": RES, "n_wl": N_WL,
                   "complete": True, "rows": out}, fh, indent=2)
    print(f"\nwritten -> {mpath}")

    for axis, label, lst in (("A", "N_BODY", body_list),
                             ("B", "N_TAPER", taper_list)):
        rows_a = [r for r in out if r["axis"] == axis]
        if axis == "B":
            # the (body_ref, taper_ref) point lives on axis A but IS an axis-B
            # point; without it the taper curve is missing its anchor and the
            # smoke printed a one-row axis B.
            shared = [r for r in out if r["axis"] == "A"
                      and r["n_body"] == body_ref and r["n_taper"] == taper_ref]
            rows_a = sorted(rows_a + shared, key=lambda r: r["n_taper"])
        if len(rows_a) < 2:
            continue
        print(f"\n--- axis {axis} ({label}) ---")
        ref = rows_a[-1]
        for r in rows_a:
            print(f"  {label}={r['n_body' if axis=='A' else 'n_taper']:>3} "
                  f"IL {r['excess_loss_db']:.4f} "
                  f"(d vs largest {r['excess_loss_db']-ref['excess_loss_db']:+.4f}, "
                  f"d vs FDTD {r['d_IL_vs_fdtd']:+.4f})")
        span = max(r["excess_loss_db"] for r in rows_a) - \
               min(r["excess_loss_db"] for r in rows_a)
        print(f"  IL span across this axis: {span:.4f} dB")
        # A verdict needs enough points to see a trend AND a plateau. The
        # res-15 smoke printed "truncation REFUTED" off TWO points -- the same
        # premature-conclusion bug the pre-registered stop rule had before it
        # was gated on grid completeness. Three points is the minimum that can
        # distinguish "flat" from "still moving".
        if len(rows_a) < 3:
            print(f"  -> only {len(rows_a)} points on this axis: NO VERDICT. "
                  "A trend needs at least 3.")
        elif span < 0.01:
            print(f"  -> FLAT in {label} (span < 0.01 dB): truncation on this "
                  "axis is REFUTED as the source of the EME-vs-FDTD offset.")
        else:
            print(f"  -> IL MOVES with {label}. Check whether it plateaus and "
                  "whether it moves TOWARD the FDTD value.")
    print("\n⚠ A plateau at the largest value TRIED is not a ceiling. If IL is "
          "still moving at the top of a list, extend it.")
    raise SystemExit(0)


# ---------------------------------------------------------------------------
# --res-convergence: check 3 (resolution) for the WAISTED body, EME side.
#
# WHY THIS IS NOT A FORMALITY. The whole step-6 grid ran at res 40, but the
# straight 2x2's own check 3 measured Delta-beta_01 needing **res 100** for
# 0.05% (job 11590785, resolved/mmi_2x2_validation.md). The design point's
# LENGTH is set by phase accumulation -- L_pi = pi/Delta-beta -- so an
# unconverged Delta-beta moves the optimum length, and a ~1% error is ~0.2 um,
# comparable to the grid's own 0.5 um step.
#
# The criterion's control was also measured at res 40, so the pass/fail
# COMPARISON is at matched resolution and partly self-cancelling. What is not
# protected is the absolute length. This mode measures both: the waist
# cross-section's n_eff_0/n_eff_1/Delta-beta_01/L_pi per resolution, and
# whether the best length in a window around the design point MOVES.
#
# One basis per resolution; lengths inside a resolution stay free.
# ---------------------------------------------------------------------------
if "--res-convergence" in sys.argv:
    print("=" * 88)
    print("CHECK 3 (resolution), EME side, at the waisted design point")
    print("=" * 88, flush=True)
    ww = _env_f("MMI2X2_RC_WW", 2.24)
    lw = _env_f("MMI2X2_RC_LW", 8.0)
    lt = _env_f("MMI2X2_RC_LT", 6.0)
    l_design = _env_f("MMI2X2_RC_L", 20.5)
    res_list = _env_list("MMI2X2_RC_RES", "40,60,80,100")
    win = _env_f("MMI2X2_RC_WIN", 1.0)
    step = _env_f("MMI2X2_RC_STEP", 0.25)
    wl_mid = float(wavelengths[len(wavelengths) // 2])
    Ls = [l_design + k * step
          for k in range(-int(win / step), int(win / step) + 1)]
    Ls = [L for L in Ls
          if (L - 2 * lt - lw) / 2 >= HYBRID_BODY_STUB - 1e-9
          and L < MM_LENGTH_STD]
    print(f"design point Ww={ww} Lw={lw} Lt={lt} L={l_design}")
    print(f"resolutions {res_list}; length window {Ls}")
    print(f"waist cross-section n_eff read at wl={wl_mid:.4f} um\n", flush=True)

    # cell index of the single waist cell in the 9-region layout
    _, cells_n, _ = regions(ww, lw, lt, 1.0)
    WAIST_CELL = sum(cells_n[:4])
    out = []
    basis = None
    for res in res_list:
        t0 = time.time()
        # build_basis() reads the module-level RES when it constructs its
        # EMESolver, so swapping it here is what changes the mode-solve grid.
        prev_res, RES = RES, res
        _IFACE_CACHE.clear()
        # Drop the PREVIOUS resolution's basis before building the next one.
        # Without this the old basis stays referenced through the assignment
        # below, so peak memory spans two resolutions -- which is what
        # OOM-killed job 12696267 at 96 GB during its res-100 basis (mode-field
        # arrays scale as res^2, and res 100 is 6.25x res 40 on its own).
        basis = None
        gc.collect()
        basis = build_basis(ww, lw, lt, max(HYBRID_BODY_STUB, 1.0))
        dt = time.time() - t0
        # Delta-beta_01 of the UNIFORM WAIST cross-section -- step 8's own note
        # says to state WHICH cross-section is being asked about, since a
        # waisted body has no single one. This is the waist.
        wl_i = len(wavelengths) // 2
        wmodes = basis["bar"][1][wl_i][WAIST_CELL]
        n0, n1 = float(np.real(wmodes[0].neff)), float(np.real(wmodes[1].neff))
        dbeta = 2 * np.pi / wl_mid * (n0 - n1)
        l_pi = np.pi / dbeta
        rows_r = []
        for L in Ls:
            s_stub = (L - 2 * lt - lw) / 2
            cb, cc = eval_from_basis(basis, s_stub, lt, lw,
                                     iface_cache=IFACE_CACHE)
            m = metrics(cb, cc, float(L))
            passed, _ = score(m)
            rows_r.append({"mm_length": float(L), "passed": passed, **m})
        best = min(rows_r, key=lambda r: r["excess_loss_db"])
        at_design = min(rows_r,
                        key=lambda r: abs(r["mm_length"] - l_design))
        print(f"-- res {res:g}  (basis {dt:.0f}s)")
        print(f"   waist n_eff_0={n0:.6f} n_eff_1={n1:.6f} "
              f"dbeta_01={dbeta:.6f} /um  L_pi={l_pi:.4f} um")
        print(f"   best-IL length in window: {best['mm_length']:.3f} um "
              f"(IL {best['excess_loss_db']:.4f})")
        print(f"   at design L={at_design['mm_length']:.3f}: "
              f"IL {at_design['excess_loss_db']:.4f} "
              f"imb {at_design['imbalance']:.4f} "
              f"phase {at_design['phase_deg']:.3f} "
              f"maxdev {at_design['phase_max_dev_deg']:.3f} "
              f"{'PASS' if at_design['passed'] else 'FAIL'}", flush=True)
        out.append({"res": res, "basis_seconds": dt, "n_eff_0": n0,
                    "n_eff_1": n1, "dbeta_01": dbeta, "L_pi": l_pi,
                    "best_length": best["mm_length"], "rows": rows_r,
                    "at_design": at_design})
        # Write after EVERY resolution, not after the loop. A res-100 basis is
        # ~3.5 h; overrunning the walltime with 40/60/80 finished and unwritten
        # is the same defect that cost job 12630600 its whole Ww=2.52 half.
        rtag_partial = "-".join(f"{r['res']:g}" for r in out)
        with open(os.path.join(output_dir,
                               f"res_convergence_res{rtag_partial}.json"),
                  "w") as fh:
            json.dump({"design": {"Ww": ww, "Lw": lw, "Lt": lt,
                                  "L": l_design},
                       "wl_mid": wl_mid, "complete": False, "rows": out},
                      fh, indent=2)
        RES = prev_res
    rtag = "-".join(f"{r:g}" for r in res_list)
    rpath = os.path.join(output_dir, f"res_convergence_res{rtag}.json")
    with open(rpath, "w") as fh:
        json.dump({"design": {"Ww": ww, "Lw": lw, "Lt": lt, "L": l_design},
                   "wl_mid": wl_mid, "complete": True, "rows": out}, fh,
                  indent=2)
    print(f"\nwritten -> {rpath}")
    print(f"\n{'res':>6} {'dbeta_01':>11} {'L_pi':>9} {'d(L_pi)%':>9} "
          f"{'best L':>8} {'IL@design':>10}")
    ref = out[-1]
    for r in out:
        dp = 100 * (r["L_pi"] - ref["L_pi"]) / ref["L_pi"]
        print(f"{r['res']:6.0f} {r['dbeta_01']:11.6f} {r['L_pi']:9.4f} "
              f"{dp:+9.3f} {r['best_length']:8.3f} "
              f"{r['at_design']['excess_loss_db']:10.4f}")
    print(f"\n(d(L_pi)% is relative to the highest resolution run, "
          f"res {ref['res']:g}.)")
    moved = {r["best_length"] for r in out}
    if len(moved) > 1:
        pairs = ", ".join(f"res {r['res']:g} -> {r['best_length']:.3f}"
                          for r in out)
        print(f"⚠ the best length in the window MOVES with resolution: {pairs}"
              " -- the design point is resolution-dependent and none of the "
              "resolutions tried fixes it.")
    else:
        print(f"best length is stable at {moved.pop():.3f} um across all "
              f"resolutions tried ({', '.join(f'{r:g}' for r in res_list)}).")
    raise SystemExit(0)


# ---------------------------------------------------------------------------
# --rescore-passers: re-evaluate every PASSING row and report the phase check
# under the MAX-over-band reading of the pre-registered text, alongside the
# mean reading that score() actually applied. Costs one basis per Ww present
# among the passers (~19 min each) plus ~7.5 s per passer -- cheap enough that
# there is no excuse for leaving the discrepancy unquantified.
#
# This does NOT re-rank anything on its own. It answers one question: does the
# lenient phase reading change WHICH points pass?
# ---------------------------------------------------------------------------
if "--rescore-passers" in sys.argv:
    print("=" * 88)
    print("RESCORE: passing rows under mean-phase (as scored) vs max-phase")
    print("=" * 88, flush=True)
    all_rows = load_all_rows()
    pass_rows = [r for r in all_rows if r["passed"]]
    if not pass_rows:
        raise SystemExit(f"no passers among {len(all_rows)} rows on disk")
    ph_ctrl, ph_tol = THR["control_phase_deg"], THR["max_phase_dev_deg"]
    print(f"{len(pass_rows)} passers of {len(all_rows)} rows; "
          f"control phase {ph_ctrl:.3f} deg, tol {ph_tol} deg\n", flush=True)
    out, flipped = [], []
    for ww in sorted({r["Ww"] for r in pass_rows}):
        grp = sorted([r for r in pass_rows if r["Ww"] == ww],
                     key=lambda r: r["total_device_length"])
        t0 = time.time()
        lt_ref, lw_ref = LT_LIST[len(LT_LIST) // 2], LW_LIST[len(LW_LIST) // 2]
        basis = build_basis(ww, lw_ref, lt_ref, max(HYBRID_BODY_STUB, 1.0))
        print(f"-- Ww={ww}: basis in {time.time()-t0:.0f}s, {len(grp)} passers",
              flush=True)
        if IFACE_CACHE:
            selfcheck_iface_cache(basis, SELFCHECK_PTS)
        print(f"{'Lw':>5} {'Lt':>4} {'L':>7} {'total':>6} {'mean-dev':>9} "
              f"{'MAX-dev':>8} {'spread':>7} {'max ok':>7}")
        for r in grp:
            cb, cc = eval_from_basis(basis, r["s"], r["Lt"], r["Lw"],
                                     iface_cache=IFACE_CACHE)
            m = metrics(cb, cc, r["mm_length"])
            mean_dev = abs(m["phase_deg"] - ph_ctrl)
            mx = m["phase_max_dev_deg"]
            ok = mx <= ph_tol
            if not ok:
                flipped.append({**r, "phase_max_dev_deg": mx,
                                "phase_spread_deg": m["phase_spread_deg"]})
            print(f"{r['Lw']:5.1f} {r['Lt']:4.1f} {r['mm_length']:7.3f} "
                  f"{r['total_device_length']:6.2f} {mean_dev:9.3f} "
                  f"{mx:8.3f} {m['phase_spread_deg']:7.3f} "
                  f"{('ok' if ok else 'FAILS'):>7}", flush=True)
            out.append({**r, "phase_max_dev_deg": mx,
                        "phase_spread_deg": m["phase_spread_deg"],
                        "mean_dev_deg": mean_dev, "passes_max_reading": ok})
    rpath = os.path.join(output_dir, f"rescore_maxphase_{SHARD_FP}.json")
    with open(rpath, "w") as fh:
        json.dump({"control_phase_deg": ph_ctrl, "tol_deg": ph_tol,
                   "n_passers": len(pass_rows), "n_fail_max_reading":
                   len(flipped), "rows": out}, fh, indent=2)
    print(f"\nwritten -> {rpath}")
    print(f"{len(flipped)} of {len(pass_rows)} passers FAIL the max-over-band "
          "reading of the same pre-registered phase text.")
    surv = [r for r in out if r["passes_max_reading"]]
    if surv:
        b = min(surv, key=lambda r: (r["total_device_length"],
                                     r["excess_loss_db"]))
        print(f"shortest passer surviving BOTH readings: Ww={b['Ww']} "
              f"Lw={b['Lw']} Lt={b['Lt']} L={b['mm_length']:.3f} -> "
              f"{b['total_device_length']:.2f} um "
              f"(IL {b['excess_loss_db']:.4f}, imb {b['imbalance']:.4f}, "
              f"max-phase-dev {b['phase_max_dev_deg']:.3f} deg)")
    else:
        print("NO passer survives the max reading -- the phase criterion, read "
              "as pre-registered in words, rejects the entire passer set.")
    raise SystemExit(0)


# ---------------------------------------------------------------------------
# --confirm: re-solve the winner DIRECTLY, plus its neighbours in L.
#
# The analytic path is validated to 1e-13 (12624884) and its interface cache to
# 0.000e+00 (12674313), so this is NOT a check on the fast path -- it is the
# check the script's own closing warning demands: that the design point being
# quoted is a real solve, and that it is not perched on a grid edge. The
# shortest passer recovered from 12630600 sat 1e-4 dB inside the loss
# threshold, which is exactly the situation a one-grid-step neighbour test is
# for.
#
# Points come from MMI2X2_FS_CONFIRM="Ww,Lw,Lt,L;..." if set; otherwise the
# shortest passer on disk and its two L-neighbours on the grid.
# ---------------------------------------------------------------------------
if "--confirm" in sys.argv:
    print("=" * 88)
    print("CONFIRM: direct per-length solves at the design point")
    print("=" * 88, flush=True)
    spec = os.environ.get("MMI2X2_FS_CONFIRM", "").strip()
    pts = []
    if spec:
        # Fields may be separated by "," or ":" -- `sbatch --export` treats a
        # comma as ITS OWN pair separator, so a comma-formatted spec passed
        # that way arrives shredded (it silently became MMI2X2_FS_CONFIRM=2.24
        # plus junk vars). Prefer ":" when going through --export, or set the
        # variable in the submitting shell with --export=ALL.
        for chunk in spec.replace(":", ",").split(";"):
            if chunk.strip():
                vals = [float(v) for v in chunk.split(",") if v.strip()]
                if len(vals) != 4:
                    raise SystemExit(
                        f"MMI2X2_FS_CONFIRM chunk {chunk!r} has {len(vals)} "
                        "fields, expected 4 (Ww,Lw,Lt,L). If this came through "
                        "`sbatch --export`, the commas were eaten -- use ':' "
                        "between fields or export the var in the shell.")
                pts.append(tuple(vals))
        print(f"points from MMI2X2_FS_CONFIRM: {pts}")
    else:
        all_rows = load_all_rows()
        passers_d = [r for r in all_rows if r["passed"]]
        if not passers_d:
            raise SystemExit(f"no passer among {len(all_rows)} rows on disk "
                             "-- nothing to confirm")
        best_d = min(passers_d, key=lambda r: (r["total_device_length"],
                                               r["excess_loss_db"]))
        dL = (L_MAX - L_MIN) / (N_L - 1)
        base = (best_d["Ww"], best_d["Lw"], best_d["Lt"])
        pts = [(*base, best_d["mm_length"] + k * dL) for k in (-1, 0, 1)]
        pts = [q for q in pts if L_MIN - 1e-9 <= q[3] < MM_LENGTH_STD
               and (q[3] - 2 * q[2] - q[1]) / 2 >= HYBRID_BODY_STUB - 1e-9]
        print(f"shortest passer of {len(passers_d)} on disk: Ww={base[0]} "
              f"Lw={base[1]} Lt={base[2]} L={best_d['mm_length']:.3f}")
        print(f"confirming it and its L-neighbours (dL={dL:.3f}): "
              f"{[f'{q[3]:.3f}' for q in pts]}")
    analytic = {key_of(*q[:3], q[3]): r
                for q in pts
                for r in load_all_rows()
                if key_of(r["Ww"], r["Lw"], r["Lt"], r["mm_length"])
                == key_of(*q[:3], q[3])}
    print(f"\n{'Ww':>6} {'Lw':>5} {'Lt':>5} {'L':>7} {'total':>7} {'IL':>8} "
          f"{'imb':>8} {'phase':>9} {'PASS':>5} {'|dS| vs analytic':>18}")
    out_rows = []
    for (ww, lw, lt, L) in pts:
        t0 = time.time()
        cb, cc = solve_direct(ww, lw, lt, L)
        m = metrics(cb, cc, float(L))
        passed, checks = score(m)
        k = key_of(ww, lw, lt, L)
        ar = analytic.get(k)
        dtxt = "-- (no analytic row)"
        if ar is not None:
            # re-derive the analytic coefficients is not possible without the
            # basis, so compare the METRICS the decision was made on
            dtxt = (f"dIL={m['excess_loss_db'] - ar['excess_loss_db']:+.2e} "
                    f"dphi={m['phase_deg'] - ar['phase_deg']:+.2e}")
        mx = m["phase_max_dev_deg"]
        print(f"{ww:6.2f} {lw:5.1f} {lt:5.1f} {L:7.3f} "
              f"{m['total_device_length']:7.2f} {m['excess_loss_db']:8.4f} "
              f"{m['imbalance']:8.4f} {m['phase_deg']:9.3f} {str(passed):>5} "
              f"{dtxt:>18}  ({time.time()-t0:.0f}s)", flush=True)
        print(f"       phase: mean-dev {abs(m['phase_deg'] - THR['control_phase_deg']):.3f} deg "
              f"(SCORED, <=2.0) | MAX-dev over band {mx:.3f} deg "
              f"({'ok' if mx <= THR['max_phase_dev_deg'] else 'EXCEEDS 2.0'}) "
              f"| band spread {m['phase_spread_deg']:.3f} deg", flush=True)
        out_rows.append({"Ww": ww, "Lw": lw, "Lt": lt, "mm_length": float(L),
                         "passed": passed, "checks": checks,
                         "analytic_row": ar, "direct": m})
    # tag-keyed: two confirm jobs in the same session wrote this file
    # untagged on 2026-09-01 and the second silently clobbered the first's
    # rows (both survived only in their SLURM logs).
    cpath = os.path.join(output_dir, f"confirm_{SHARD_FP}_{SHARD_TAG}.json")
    with open(cpath, "w") as fh:
        json.dump({"criterion_used": CRIT, "rows": out_rows}, fh, indent=2)
    print(f"\nwritten -> {cpath}")
    ok_pts = [r for r in out_rows if r["passed"]]
    print(f"{len(ok_pts)}/{len(out_rows)} confirm points pass on a DIRECT solve.")
    if len(ok_pts) < len(out_rows):
        print("  ⚠ the design point has a passing NEIGHBOUR that fails -- it is")
        print("  on a grid edge. Quote a point with margin, not this one.")
    raise SystemExit(0)


print("=" * 88)
print("2x2 waist -- FREE-LENGTH-SWEEP grid")
print(f"  control {CONTROL['total_device_length']:.2f} um total; candidates must")
print(f"  be SHORTER (hard filter: mm_length < {MM_LENGTH_STD})")
print(f"  Ww{WW_LIST} (only non-free axis) x Lw{LW_LIST} x Lt{LT_LIST}")
print(f"  x {N_L} lengths in [{L_MIN}, {L_MAX}]  -- all free")
print(f"  C_TRANS={C_TRANS} (5a floor {floor}), res={RES}, {N_WL} wl")
print("=" * 88, flush=True)
print("NOTE: --validate confirmed this path to 1e-13 (job 12624884).")
print(f"  shard tag {SHARD_TAG!r}; rows appended to {os.path.basename(ROWS_PATH)}")
print(f"  interface memoisation: {'ON (self-checked below)' if IFACE_CACHE else 'OFF'}")

# --summarize: aggregate every shard on disk and print the ranking. No solves,
# so the winner can be read off between jobs (and after a timeout).
if "--summarize" in sys.argv:
    WW_LIST = []
    print("  --summarize: no solving, aggregating shards on disk only")

done = load_done()
if done:
    print(f"  RESUME: {len(done)} points already on disk across all shards; "
          f"they will be skipped", flush=True)

rows = []
t_solve = 0.0
fh_rows = open(ROWS_PATH, "a", buffering=1)
print(f"\n{'Ww':>6} {'Lw':>5} {'Lt':>5} {'s':>6} {'L':>7} {'total':>7} "
      f"{'IL':>8} {'imb':>8} {'phase':>9} {'PASS':>5} {'fails':>20}")
for ww in WW_LIST:
    todo = [(lw, lt, float(L))
            for lw in LW_LIST for lt in LT_LIST
            for L in np.linspace(L_MIN, L_MAX, N_L)
            if (L - 2 * lt - lw) / 2 >= HYBRID_BODY_STUB - 1e-9
            and L < MM_LENGTH_STD
            and key_of(ww, lw, lt, float(L)) not in done]
    if not todo:
        print(f"\n-- Ww={ww}: all points already on disk, no basis solve",
              flush=True)
        continue
    t0 = time.time()
    # reference geometry for the basis: mid-grid, guaranteed legal
    lt_ref, lw_ref = LT_LIST[len(LT_LIST) // 2], LW_LIST[len(LW_LIST) // 2]
    s_ref = max(HYBRID_BODY_STUB, 1.0)
    basis = build_basis(ww, lw_ref, lt_ref, s_ref)
    dt = time.time() - t0
    t_solve += dt
    print(f"\n-- Ww={ww}: basis solved in {dt:.0f}s "
          f"(ONE solve per port; every length below is free); "
          f"{len(todo)} points to do", flush=True)
    if IFACE_CACHE:
        selfcheck_iface_cache(basis, SELFCHECK_PTS)
    for (lw, lt, L) in todo:
        s_stub = (L - 2 * lt - lw) / 2
        cb, cc = eval_from_basis(basis, s_stub, lt, lw,
                                 iface_cache=IFACE_CACHE)
        m = metrics(cb, cc, float(L))
        passed, checks = score(m)
        fails = ",".join(k for k, v in checks.items() if not v) or "-"
        row = {"Ww": ww, "Lw": lw, "Lt": lt, "s": float(s_stub),
               "mm_length": float(L), "passed": passed,
               "checks": checks, **m}
        rows.append(row)
        # EVERY row hits disk as it is produced. Job 12630600 timed out at 8 h
        # with the whole Ww=2.52 half of the grid finished and NOTHING written,
        # because the JSON was dumped only after the loop.
        fh_rows.write(json.dumps(row) + "\n")
        done.add(key_of(ww, lw, lt, L))
        if passed or len(rows) % 40 == 0:
            print(f"{ww:6.2f} {lw:5.1f} {lt:5.1f} {s_stub:6.2f} {L:7.3f} "
                  f"{m['total_device_length']:7.2f} {m['excess_loss_db']:8.4f} "
                  f"{m['imbalance']:8.4f} {m['phase_deg']:9.3f} "
                  f"{str(passed):>5} {fails:>20}", flush=True)
fh_rows.close()
rows = load_all_rows()

# How many legal points does the CURRENT grid definition contain? The stop
# rule is a pre-registered decision about the whole grid, so it must not be
# evaluated on a partial aggregate -- a --summarize over 3 smoke rows printed
# "STOP RULE TRIGGERED" before this guard existed, which is precisely how a
# negative gets published off an unfinished sweep.
# LITERAL, not the env-narrowed lists: a shard job runs one Ww, and reading
# WW_LIST here would let it declare the whole grid COMPLETE with the other Ww
# never simulated. This is the pre-registered grid of the validation doc.
PREREG_WW = [2.52, 2.24]
PREREG_LW = [4.0, 6.0, 8.0, 10.0]
PREREG_LT = [3.0, 4.0, 5.0, 6.0]
PREREG_L = np.linspace(16.0, 27.5, 24)
expected = [(ww, lw, lt, float(L))
            for ww in PREREG_WW for lw in PREREG_LW for lt in PREREG_LT
            for L in PREREG_L
            if (L - 2 * lt - lw) / 2 >= HYBRID_BODY_STUB - 1e-9
            and L < MM_LENGTH_STD]
have = {key_of(r["Ww"], r["Lw"], r["Lt"], r["mm_length"]) for r in rows}
missing = [q for q in expected if key_of(*q) not in have]
COMPLETE = not missing

print("\n" + "=" * 88)
print(f"evaluated {len(rows)} candidates from {t_solve:.0f}s of mode solves "
      f"({len(WW_LIST)} bases x 2 ports)")
print(f"grid coverage: {len(expected) - len(missing)}/{len(expected)} legal "
      f"points on disk -- {'COMPLETE' if COMPLETE else 'PARTIAL'}")
if not COMPLETE:
    print(f"  {len(missing)} points still owed; re-submit the same job to "
          "resume (finished points are skipped)")
passers = [r for r in rows if r["passed"]]
if passers:
    best = min(passers, key=lambda r: r["total_device_length"])
    print(f"{len(passers)} PASS all four criteria.")
    print(f"BEST (shortest passer): Ww={best['Ww']} Lw={best['Lw']} "
          f"Lt={best['Lt']} L={best['mm_length']:.3f} (stub {best['s']:.2f})")
    save = CONTROL["total_device_length"] - best["total_device_length"]
    print(f"  total {best['total_device_length']:.2f} um vs control "
          f"{CONTROL['total_device_length']:.2f} -> saves {save:.2f} um "
          f"({100*save/CONTROL['total_device_length']:.1f}%)")
    print(f"  IL {best['excess_loss_db']:.4f} dB | imb {best['imbalance']:.4f} "
          f"| phase {best['phase_deg']:.3f} deg")
    print("\n  ⚠ CONFIRM before quoting: re-solve this point DIRECTLY (the")
    print("  analytic path is validated to 1e-13, but confirm the winner), and")
    print("  check it is not sitting on a grid edge.")
else:
    print(f"NO candidate of {len(rows)} passes all four.")
    near = [r for r in rows
            if r["excess_loss_db"] <= 1.05 * CONTROL["max_excess_loss_db"]]
    if rows:
        closest = min(rows, key=lambda r: r["excess_loss_db"])
        print(f"  closest excess loss {closest['excess_loss_db']:.4f} dB "
              f"(control {CONTROL['max_excess_loss_db']:.4f})")
        print(f"  within 5% of control's loss: {len(near)}")
    if not near and COMPLETE:
        print("\n  ** PRE-REGISTERED STOP RULE TRIGGERED **")
        print("  A clean negative here is a MECHANISM RESULT: it refutes")
        print("  transition radiation as the 1x4's explanation, since this")
        print("  device removes both conditions that hypothesis depends on.")
    elif not near:
        print("\n  stop rule NOT evaluated: the grid is PARTIAL. A negative")
        print("  is only a mechanism result over the whole pre-registered")
        print("  grid; finish the missing points first.")

SUMMARY_PATH = os.path.join(output_dir,
                            f"free_sweep_{SHARD_FP}_{SHARD_TAG}.json")
with open(SUMMARY_PATH, "w") as fh:
    json.dump({"criterion_used": CRIT, "n_pass": len(passers),
               "mode_solve_seconds": t_solve,
               "shards": [os.path.basename(f) for f in _shard_files()],
               "iface_cache": IFACE_CACHE,
               "grid_complete": COMPLETE,
               "n_expected": len(expected), "n_missing": len(missing),
               "rows": rows}, fh, indent=2)
print(f"\nwritten -> {SUMMARY_PATH}")

