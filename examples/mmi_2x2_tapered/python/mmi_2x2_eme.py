"""
Tutorial: simulate the WAISTED 2x2 MMI splitter with the BEEF library using EME.

The pure-EME counterpart of this folder's other two scripts. No FDTD anywhere,
so it is the fast way to look at this device — seconds to minutes, against
hours for mmi_2x2_fdtd.py and ~27 core-hours for mmi_2x2_hybrid_alt.py.

Geometry: the validated design point — a 2.8 µm body that narrows to a
2.24 µm interior waist and ramps back out, reaching the same 50/50 split with
90 degree quadrature in 31.00 µm total against the straight control's 37.00 µm
(16.2 % shorter). See validation_markdown/resolved/mmi_2x2_waist_body_validation.md.

NINE REGIONS, NOT THE STRAIGHT DEVICE'S FIVE. `examples/mmi_2x2`'s EME script
uses wg | taper | body | taper | wg, because there the body is one uniform
slab. A waisted body is three different things, so the body expands into five
regions and the chain becomes:

    idx  region     length   cells  modes   what
     0   wg in        0.50     1      1     single-mode access waveguide
     1   taper in     4.00    15      7     port taper
     2   stub         3.00     1     15     uniform 2.80 µm
     3   ramp         5.00     8     15     2.80 -> 2.24 µm  (needs the staircase)
     4   waist        6.00     1     15     uniform 2.24 µm
     5   ramp         5.00     8     15     2.24 -> 2.80 µm
     6   stub         3.00     1     15     uniform 2.80 µm
     7   taper out    4.00    15      7     port taper
     8   wg out       0.50     1      1     single-mode access waveguide
                    -------
                     31.00 = the criterion's total_device_length

This is the SAME decomposition the step-6 design sweep scored the whole grid
on (`examples/verification/python/mmi_2x2_waist_free_sweep.py`), reused here
deliberately rather than re-invented, so this script's numbers are directly
comparable to the grid that chose the design point.

ONE CELL FOR EVERY UNIFORM REGION. Propagation through a uniform region is
analytic, so one cell is exact there; extra cells would only insert identity
interfaces between identical cross-sections, and each such interface sheds a
little |E|^2. Only the two ramps get the C_TRANS staircase.

⚠ N_TAPER = 7 IS NOT A CONVERGED VALUE, and this is the single most important
caveat on every number below. Step 5b measured the scoring to move 0.1217 dB
over N_TAPER = 5..13, and the design point passes the phase criterion ONLY at
7 — the value the grid happened to use. 7 is therefore the basis this device
was SELECTED on, not a basis it has been shown to be insensitive to. Keep it
to stay comparable with the grid; do not read it as converged. It is exposed
as MMI2X2T_EME_NTAPER so the sensitivity can be re-measured.

⚠ C_TRANS = 8 is the step-5a floor for resolving the ramps, and it is an
OPTIMUM rather than a floor for insertion loss: the hybrid's IL overshoots
full-device FDTD and climbs with cell count past ct16. Raising it does not
make the answer better.

⚠ N_BODY = 15 is over-requested in the WAIST specifically. The 2.24 µm waist
supports exactly 10 distinct modes (N=11 is an exact duplicate, overlap
1.000000); 15 is right for the 2.80 µm faces. It is kept at 15 across all five
body regions because that is what the design grid used, and changing it here
would silently score this device on a different basis than the one that chose
it. The duplicates are harmless, just wasted solve time.

mode_solver = "tidy3d": there is no FDTD section anywhere in this script for
the basis to stay consistent with, which is the reason the hybrid scripts pin
"meep" instead.

Run from any directory:
    python /path/to/examples/mmi_2x2_tapered/python/mmi_2x2_eme.py

    MMI2X2T_EME_NWL=3 python .../mmi_2x2_eme.py        # quick look
    MMI2X2T_EME_CTRANS=16 python .../mmi_2x2_eme.py    # staircase sensitivity
"""

import functools
import hashlib
import json
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")           # headless — swap for "Qt5Agg" for interactive

from beef import (MMI, EMESolver, SMatrix, SOI_220, plot_eme_geometry,
                  plot_performance)

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEVICE_OUT = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_2x2_tapered"))
output_dir = os.path.join(_DEVICE_OUT, "mmi_2x2_eme")
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_sim_cache"))
# ⚠ OVERRIDABLE SO THE MODE-MEMO A/B CAN BE MEASURED AT ALL. solve() caches
# the WHOLE call under save_dir, so a second run of this script with the same
# config is served entirely from disk and never reaches compute_modes(). A
# memo-off/memo-on A/B sharing one save_dir therefore measures the save_dir
# cache and not the memo, and its |dS| comparison is vacuous because both arms
# return the same stored array. Measured 2026-09-15: arm 1 591 s, arm 2 7 s,
# |dS| = 0 -- and the memo counters never fired at all, which is what exposed
# it. Point each arm at its own empty directory to measure the memo.
save_dir = os.environ.get("MMI2X2T_EME_SAVEDIR", save_dir)
# The whole-device FDTD reference written by mmi_2x2_fdtd.py in this folder.
FDTD_REF_DIR = os.path.join(_DEVICE_OUT, "mmi_2x2_fdtd")


def _env_f(name, default):
    return float(os.environ.get(name, default))


def _env_i(name, default):
    return int(os.environ.get(name, default))


# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by the EME solver
WL_MIN = 1.530      # wavelength sweep start (µm)
WL_MAX = 1.565      # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points — the criterion's own grid

N_PORT_MODES = 1    # single-mode access waveguide
N_TAPER_MODES = 7   # ⚠ NOT converged — see the module docstring
N_BODY_MODES = 15   # the 2.80 µm faces; over-requested in the 2.24 µm waist
TAPER_CELLS = 15    # a port taper genuinely needs ~15 cells to resolve
C_TRANS = 8         # cells per width ramp — the step-5a floor, and an OPTIMUM

RESOLUTION = _env_f("MMI2X2T_EME_RES", RESOLUTION)
N_WL = _env_i("MMI2X2T_EME_NWL", N_WL)
N_TAPER_MODES = _env_i("MMI2X2T_EME_NTAPER", N_TAPER_MODES)
N_BODY_MODES = _env_i("MMI2X2T_EME_NBODY", N_BODY_MODES)
TAPER_CELLS = _env_i("MMI2X2T_EME_TAPERCELLS", TAPER_CELLS)
C_TRANS = _env_i("MMI2X2T_EME_CTRANS", C_TRANS)

# --- the design point (identical to the two siblings in this folder) -------
WAIST_WIDTH = 2.24        # Ww — interior body width (µm)
WAIST_LENGTH = 6.0        # Lw — uniform waist run (µm)
TRANSITION_LENGTH = 5.0   # Lt — each width ramp (µm)
MM_LENGTH = 22.0          # L  — total body length (µm)
STRAIGHT_LENGTH = 0.5     # access waveguide drawn at each end (µm)
TAPER_LENGTH = 4.0

# ---------------------------------------------------------------------------
# 0. Content-addressed mode memo — the reuse this 9-region layout invites
#
# This device re-presents the SAME cross-section many times over, and nothing
# in the library caches that for us: EMESolver.compute_modes() is documented
# as "exposed for caching", i.e. it provides the hook and leaves the policy to
# the caller. Both hybrid scripts in this folder already install this memo;
# without it the pure-EME script was the only one paying full price.
#
# Where the duplicates come from, for a body that is mirror-symmetric about
# x = 0 and uniform through the waist:
#   * the two STUB regions are the same uniform 2.80 um slab;
#   * the two RAMP regions are mirror images, so cell i of ramp-in presents
#     exactly the cross-section that cell (C_TRANS-1-i) of ramp-out does;
#   * the two TAPER regions likewise, and the two access waveguides;
#   * and the whole BODY is identical between the bar and cross solves, which
#     differ only in which output branch they measure.
#
# ⚠ KEYED ON CONTENT, NOT ON WIDTH, and that distinction is load-bearing.
# Width-keyed reuse was measured BREAKING PHYSICS on the 1x4 (totalT 0.9703 ->
# 0.1507): equal width is NOT equal cross-section once a cell sits off-axis,
# which is exactly the case here — input_idx/output_idx put the taper and
# access-waveguide cells at y = +/-0.75 um, not on the axis. Hashing the
# actual refractive-index arrays makes the reuse exact by construction: the
# mirrored taper cells of the CROSS solve hash differently from the bar
# solve's and are correctly re-solved, while the on-axis body cells match and
# are not. The failure mode is "no speedup", never "wrong answer".
#
# functools.wraps, NOT a subclass: subclassing a solver in this repo breaks
# cache fingerprints and silently drops cache_only/save_dir forwarding.
#
# Set MMI2X2T_EME_MEMO=0 to disable (used by the A/B below).
# ---------------------------------------------------------------------------
MEMO_ON = _env_i("MMI2X2T_EME_MEMO", 1) == 1
_MODE_MEMO, _MEMO_STATS = {}, {"hits": 0, "misses": 0, "bypassed": 0}


def _install_mode_memo():
    orig = EMESolver.compute_modes

    def _key(cs, n_modes):
        h = hashlib.blake2b(digest_size=16)
        for comp in ("nx", "ny", "nz"):
            arr = np.ascontiguousarray(getattr(cs, comp))
            h.update(np.asarray(arr.shape, dtype=np.int64).tobytes())
            h.update(arr.tobytes())
        return (h.hexdigest(), float(cs.env.wl), float(cs.env.T), int(n_modes))

    @functools.wraps(orig)
    def compute_modes(self, cs, n_modes, *args, **kwargs):
        if not MEMO_ON:
            _MEMO_STATS["bypassed"] += 1
            return orig(self, cs, n_modes, *args, **kwargs)
        try:
            k = _key(cs, n_modes)
        except Exception:          # unknown cross-section shape -> just solve
            return orig(self, cs, n_modes, *args, **kwargs)
        if k in _MODE_MEMO:
            _MEMO_STATS["hits"] += 1
            return _MODE_MEMO[k]
        _MEMO_STATS["misses"] += 1
        res = orig(self, cs, n_modes, *args, **kwargs)
        _MODE_MEMO[k] = res
        return res

    EMESolver.compute_modes = compute_modes


_install_mode_memo()
print(f"mode memo: {'ON' if MEMO_ON else 'OFF'} "
      "(content-addressed on nx/ny/nz + wl + T + n_modes)")


# ---------------------------------------------------------------------------
# 1. Define geometry
# ---------------------------------------------------------------------------
mmi = MMI(
    n_inputs=2,             # 2x2 splitter
    n_outputs=2,
    mm_length=MM_LENGTH,    # MMI body length (µm) — shortened by the waist
    mm_width=2.8,           # MMI body width at the two faces (µm)
    body_waist_width=WAIST_WIDTH,
    body_waist_length=WAIST_LENGTH,
    body_transition_length=TRANSITION_LENGTH,
    waveguide_width=0.45,   # single-mode access waveguide width (µm)
    taper_length=TAPER_LENGTH,
    taper_width=1.3,        # taper width at MMI junction (µm)
    waveguide_offset=1.5,   # center-to-center port spacing (µm)
    waveguide_height=0.22,  # Si layer thickness (µm) — SOI_220
)

# The uniform stub between the body face and the start of each ramp.
stub = (MM_LENGTH - 2 * TRANSITION_LENGTH - WAIST_LENGTH) / 2      # 3.0 µm
if stub < 0:
    raise SystemExit(
        f"REFUSING: 2*Lt + Lw = {2*TRANSITION_LENGTH + WAIST_LENGTH} µm "
        f"exceeds mm_length = {MM_LENGTH} µm; there is no room for the body.")

# ⚠ LENGTH CONVENTION. `MMI.total_length` counts mm_length + both tapers and
# stops there (30.00 µm). The criterion's `total_device_length` also includes
# the 0.5 µm access waveguide drawn at each end (31.00 µm, against the
# straight control's 37.00), which is what the 16.2 %-shorter claim is built
# on. Scoring mmi.total_length against the control would report 18.9 %.
FIXED_LENGTH = 2 * TAPER_LENGTH + 2 * STRAIGHT_LENGTH      # 9.0 µm
total_device_length = mmi.mm_length + FIXED_LENGTH

half_body = mmi.mm_length / 2
print(f"MMI body {mmi.mm_length:.2f} µm; device {total_device_length:.2f} µm "
      f"(straight control 37.00 — "
      f"{100*(37.0 - total_device_length)/37.0:.1f}% shorter)")
print(f"Body width: face {mmi.width_at(-half_body):.3f} -> "
      f"waist {mmi.width_at(0.0):.3f} -> face {mmi.width_at(half_body):.3f} µm")
print(f"Input ports at y = {mmi.input_ys} µm; outputs at y = {mmi.output_ys} µm")

# ---------------------------------------------------------------------------
# 2. Configure EME solver — the 9-region decomposition
#
# The SAME layout the step-6 design sweep scored the whole grid on, so these
# numbers are directly comparable to the grid that chose this design point.
# ---------------------------------------------------------------------------
REGION_LABELS = ["wg in", "taper in", "stub", "ramp in", "waist",
                 "ramp out", "stub", "taper out", "wg out"]
lengths = [STRAIGHT_LENGTH, TAPER_LENGTH, stub, TRANSITION_LENGTH,
           WAIST_LENGTH, TRANSITION_LENGTH, stub, TAPER_LENGTH,
           STRAIGHT_LENGTH]
num_cells = [1, TAPER_CELLS, 1, C_TRANS, 1, C_TRANS, 1, TAPER_CELLS, 1]
num_modes = [N_PORT_MODES, N_TAPER_MODES, N_BODY_MODES, N_BODY_MODES,
             N_BODY_MODES, N_BODY_MODES, N_BODY_MODES, N_TAPER_MODES,
             N_PORT_MODES]

assert abs(sum(lengths) - total_device_length) < 1e-9, (
    f"regions sum to {sum(lengths)} µm but the device is "
    f"{total_device_length} µm — the decomposition does not tile it")

print(f"\n{'idx':>3} {'region':>9} {'length':>7} {'cells':>6} {'modes':>6}")
for i, lab in enumerate(REGION_LABELS):
    print(f"{i:>3} {lab:>9} {lengths[i]:7.2f} {num_cells[i]:6d} {num_modes[i]:6d}")
print(f"{'':>3} {'total':>9} {sum(lengths):7.2f} {sum(num_cells):6d} cells")

solver = EMESolver(
    num_modes=num_modes,
    num_cells=num_cells,
    lengths=lengths,
    oxide_width=5.0,                       # cross-section mesh y-span (µm)
    oxide_height=1.5,                      # cross-section mesh z-span (µm)
    res=RESOLUTION,                        # mesh pts/µm
    material_stack=SOI_220,                # 220 nm SOI stack
    mode_solver="tidy3d",                  # no FDTD section to stay consistent with
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
plot_eme_geometry(
    mmi, solver=solver,
    region_labels=REGION_LABELS,
    title="Waisted 2x2 MMI Splitter — Top-down Geometry (EME)",
    output_path=os.path.join(output_dir, "geometry.png"),
    # output_idx omitted — draws both output branches for the reference plot
    # (each solve() below restricts to one branch at a time).
)
print(f"\nGeometry plot saved -> {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)  # µm, C-band
print(f"Sweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} - {wavelengths[-1]:.3f} µm", flush=True)

# Excite the top input (input_idx=1) and measure each output branch in turn.
# input_idx/output_idx restrict the EME cross-sections to one input and one
# output waveguide, so the mode solver sees an isolated single-mode guide at
# each end rather than a two-guide supermode. output_idx=1 is the top output
# (bar), output_idx=0 the bottom output (cross).
t0 = time.time()
smat_bar = solver.solve(mmi, wavelengths, ports=["bar"],
                        input_idx=1, output_idx=1, save_dir=save_dir)
smat_cross = solver.solve(mmi, wavelengths, ports=["cross"],
                          input_idx=1, output_idx=0, save_dir=save_dir)
elapsed = time.time() - t0
smat = smat_bar
smat._data.update(smat_cross._data)
smat.check_passive()  # total output power must not exceed 1 at any wavelength
print(f"Solved in {elapsed:.0f}s")
# THREE distinct outcomes, and conflating them is how the first A/B attempt
# fooled itself. Zero hits+misses means the memo never SAW a call, which is
# either because it is switched off or because solve() served the whole call
# from save_dir and never reached the mode solver at all — and only the
# second of those makes the run's timing a cache-read rather than a solve.
_h, _m, _b = (_MEMO_STATS["hits"], _MEMO_STATS["misses"],
              _MEMO_STATS["bypassed"])
if _h + _m:
    print(f"mode memo: {_m} cross-sections solved, {_h} reused "
          f"({100*_h/(_h+_m):.0f}% avoided) across both the bar and cross "
          f"solves and all {N_WL} wavelengths")
elif _b:
    print(f"mode memo: DISABLED — {_b} cross-sections solved with no reuse. "
          f"This is the honest baseline cost; the timing above is a real "
          f"solve.")
else:
    print("mode memo: NOT EXERCISED — compute_modes() was never called, so "
          "this run did no mode solving at all and was served entirely from "
          "save_dir's whole-call cache. Timings here are cache-read times, "
          "NOT solve times, and an A/B built on them is vacuous.")

# ---------------------------------------------------------------------------
# 5. Scores — the four criteria this device is judged on
#
# Keep the COMPLEX coefficients: the bar/cross relative phase is this device's
# whole function, and |T| discards it.
# ---------------------------------------------------------------------------
c_bar, c_cross = smat.coeff("bar"), smat.coeff("cross")
T_bar, T_cross = np.abs(c_bar) ** 2, np.abs(c_cross) ** 2
T_total = T_bar + T_cross

# Copied from example_output/verification/mmi_2x2_waist_metric/criterion.json,
# which was written BEFORE any design sweep ran and is scored automatically
# rather than eyeballed. Two of the three thresholds are defined RELATIVE to
# the straight control's own measured performance, which is why they are not
# round numbers.
CONTROL_LENGTH_UM = 37.00
CONTROL_PHASE_DEG = -89.74957705211592
MAX_EXCESS_LOSS_DB = 0.48090621635824105   # control 0.3809 + 0.1
MAX_IMBALANCE = 0.1371280008043103         # 2x control 0.06856
MAX_PHASE_DEV_DEG = 2.0

dphi = np.angle(c_cross / c_bar, deg=True)
metrics = {
    "total_device_length": float(total_device_length),
    "excess_loss_db": float(np.max(-10 * np.log10(np.clip(T_total, 1e-20, None)))),
    "imbalance": float(np.max(np.abs(T_bar - T_cross))),
    "phase_deg": float(np.mean(dphi)),
    "phase_mean_dev_deg": float(abs(np.mean(dphi) - CONTROL_PHASE_DEG)),
    "phase_max_dev_deg": float(np.max(np.abs(dphi - CONTROL_PHASE_DEG))),
    "phase_spread_deg": float(np.ptp(dphi)),
    "min_total": float(np.min(T_total)),
}
# ⚠ The pre-registered text is AMBIGUOUS between a band-MEAN and a band-MAX
# reading of the phase deviation, and it is not academic: it killed the
# original 29.50 µm winner, and the 30.50 µm point passes one and fails the
# other. Both are always reported so the verdict never depends on the reader's
# assumption.
checks = {
    "shorter": metrics["total_device_length"] < CONTROL_LENGTH_UM,
    "excess_loss": metrics["excess_loss_db"] <= MAX_EXCESS_LOSS_DB,
    "imbalance": metrics["imbalance"] <= MAX_IMBALANCE,
    "phase(mean)": metrics["phase_mean_dev_deg"] <= MAX_PHASE_DEV_DEG,
    "phase(max)": metrics["phase_max_dev_deg"] <= MAX_PHASE_DEV_DEG,
}

wl_c = wavelengths[len(wavelengths) // 2]
idx = int(np.argmin(np.abs(wavelengths - wl_c)))
print(f"\nResults at lambda = {wl_c:.3f} µm:")
print(f"  Bar output   T = {T_bar[idx]:.4f}  "
      f"IL = {10*np.log10(max(T_bar[idx],1e-20)):.2f} dB")
print(f"  Cross output T = {T_cross[idx]:.4f}  "
      f"IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total        T = {T_total[idx]:.4f}  "
      f"IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  Relative phase = {dphi[idx]:+.3f} deg  (ideal +/-90)")
print(f"  (ideal 50/50 split -> T_each=0.5, IL_total=0.0 dB)")

print(f"\nBand scores ({N_WL} wl over {WL_MIN}-{WL_MAX} µm):")
print(f"  total length     {metrics['total_device_length']:8.2f} µm "
      f"(control {CONTROL_LENGTH_UM:.2f} -> "
      f"{100*(CONTROL_LENGTH_UM - metrics['total_device_length'])/CONTROL_LENGTH_UM:.1f}% shorter)")
print(f"  excess loss      {metrics['excess_loss_db']:8.4f} dB  (<= {MAX_EXCESS_LOSS_DB:.4f})")
print(f"  imbalance        {metrics['imbalance']:8.4f}     (<= {MAX_IMBALANCE:.4f})")
print(f"  phase            {metrics['phase_deg']:8.3f}     (control {CONTROL_PHASE_DEG:.3f})")
print(f"  phase mean-dev   {metrics['phase_mean_dev_deg']:8.3f}     (<= {MAX_PHASE_DEV_DEG})")
print(f"  phase MAX-dev    {metrics['phase_max_dev_deg']:8.3f}     (<= {MAX_PHASE_DEV_DEG})")
print(f"  min total T      {metrics['min_total']:8.4f}")
fails = ", ".join(k for k, v in checks.items() if not v) or "-"
print(f"  VERDICT: {'PASS' if all(checks.values()) else 'FAIL'}  (failing: {fails})")

# ---------------------------------------------------------------------------
# 6. Compare against the whole-device FDTD reference
#
# EME is an approximation of the same device the FDTD script brute-forces, and
# on THIS body it cannot be beta-calibrated, so the gap is the measurement.
# ---------------------------------------------------------------------------
_geo = (f"Ww{WAIST_WIDTH:g}_Lw{WAIST_LENGTH:g}"
        f"_Lt{TRANSITION_LENGTH:g}_L{MM_LENGTH:g}")
ref_path = os.path.join(FDTD_REF_DIR, f"smatrix_waisted_{_geo}.pkl")
fdtd = None
if os.path.exists(ref_path):
    sm_ref = SMatrix.load(ref_path)
    fb = sm_ref.coeff(port="out_0", in_mode=1, out_mode=1)
    fc = sm_ref.coeff(port="out_1", in_mode=1, out_mode=1)
    # Match each EME wavelength to its nearest reference wavelength and REFUSE
    # if any pair is not essentially identical — comparing two engines at
    # wavelengths that merely nearly agree would fold a wavelength offset into
    # the engine disagreement being measured.
    ref_wl = np.linspace(WL_MIN, WL_MAX, len(fb))
    sel = [int(np.argmin(np.abs(ref_wl - w))) for w in wavelengths]
    mismatch = float(np.max(np.abs(ref_wl[sel] - wavelengths)))
    if mismatch > 1e-9:
        raise SystemExit(
            f"wavelength grids do not align: worst |dwl| = {mismatch:.3e} µm "
            f"(EME {len(c_bar)} pts vs reference {len(fb)}). Use an N_WL whose "
            f"grid is a subset of the reference's, e.g. N_WL=11 against 51.")
    fb, fc = fb[sel], fc[sel]
    FTb, FTc = np.abs(fb) ** 2, np.abs(fc) ** 2
    fphi = np.angle(fc / fb, deg=True)
    dTb, dTc = np.abs(T_bar - FTb), np.abs(T_cross - FTc)
    fdtd = {
        "excess_loss_db": float(np.max(-10*np.log10(np.clip(FTb+FTc, 1e-20, None)))),
        "imbalance": float(np.max(np.abs(FTb - FTc))),
        "phase_deg": float(np.mean(fphi)),
        "phase_max_dev_deg": float(np.max(np.abs(fphi - CONTROL_PHASE_DEG))),
        "max_abs_dT": {"bar": float(dTb.max()), "cross": float(dTc.max())},
        "mean_pct": {"bar": float(100*dTb.mean()/FTb.mean()),
                     "cross": float(100*dTc.mean()/FTc.mean())},
    }
    print("\n" + "=" * 74)
    print("PURE EME vs FULL-DEVICE FDTD — identical geometry and wavelengths")
    print("=" * 74)
    print(f"{'port':>6} {'max|dT|':>10} {'mean|dT|':>10} {'% of mean T_fdtd':>17}")
    print(f"{'bar':>6} {dTb.max():10.4f} {dTb.mean():10.4f} "
          f"{fdtd['mean_pct']['bar']:17.2f}")
    print(f"{'cross':>6} {dTc.max():10.4f} {dTc.mean():10.4f} "
          f"{fdtd['mean_pct']['cross']:17.2f}")
    print(f"\n{'quantity':>18} {'EME':>10} {'full FDTD':>10} {'delta':>10}")
    for q in ("excess_loss_db", "imbalance", "phase_deg", "phase_max_dev_deg"):
        print(f"{q:>18} {metrics[q]:10.4f} {fdtd[q]:10.4f} "
              f"{metrics[q]-fdtd[q]:+10.4f}")
    print("\nFor scale, the two hybrid cascades in this folder against the same")
    print("reference: EME-body 0.76 % bar / 3.17 % cross, FDTD-ramp 5.76 % /")
    print("1.20 %. Pure EME has no FDTD anywhere, so it is the cheapest point")
    print("on that accuracy/cost curve, not a competitor to them.")
else:
    print(f"\n(no whole-device FDTD reference at {ref_path} to compare "
          "against; run mmi_2x2_fdtd.py in this folder for the head-to-head)")

# ---------------------------------------------------------------------------
# 7. Save results
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved -> {output_dir}/smatrix.pkl")

# The filename carries every parameter the result depends on, GEOMETRY FIRST —
# this study lost results to filename collisions eight separate times, and
# N_TAPER in particular changes the scoring by 0.1217 dB.
tag = (f"{_geo}_res{RESOLUTION:g}_nwl{N_WL}_ct{C_TRANS}"
       f"_nt{N_TAPER_MODES}_tc{TAPER_CELLS}_nb{N_BODY_MODES}")
with open(os.path.join(output_dir, f"eme_{tag}.json"), "w") as fh:
    json.dump({"config": {"res": RESOLUTION, "n_wl": N_WL, "c_trans": C_TRANS,
                          "n_taper_modes": N_TAPER_MODES,
                          "n_body_modes": N_BODY_MODES,
                          "taper_cells": TAPER_CELLS,
                          "mode_solver": "tidy3d",
                          "mode_memo": MEMO_ON,
                          "mode_solves": _MEMO_STATS["misses"],
                          "mode_reuses": _MEMO_STATS["hits"],
                          "mode_bypassed": _MEMO_STATS["bypassed"],
                          "regions": {"labels": REGION_LABELS,
                                      "lengths": lengths,
                                      "num_cells": num_cells,
                                      "num_modes": num_modes}},
               "design": {"Ww": WAIST_WIDTH, "Lw": WAIST_LENGTH,
                          "Lt": TRANSITION_LENGTH, "L": MM_LENGTH},
               "seconds": elapsed,
               "wavelengths": wavelengths.tolist(),
               "T_bar": T_bar.tolist(), "T_cross": T_cross.tolist(),
               "phase_deg_per_wl": dphi.tolist(),
               "metrics": metrics, "checks": checks,
               "passed": all(checks.values()),
               "vs_full_fdtd": fdtd}, fh, indent=2)
print(f"Scores saved -> eme_{tag}.json")

# Plot insertion loss and transmission vs wavelength — plot_performance()
# defaults to every port in smat, so both outputs show up automatically.
plot_performance(smat, output_dir=output_dir,
                 title="Waisted 2x2 MMI Splitter (EME)")
print(f"Performance plot saved -> {output_dir}/performance.png")
