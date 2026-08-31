"""2x2 MMI waist body -- STEP 6: the design sweep, scored against step 4's
PRE-REGISTERED criterion.

validation_markdown/mmi_2x2_waist_body_validation.md step 6.

GATED ON STEP 5a, which converged at C_TRANS=8 (job 12555775). Every point
here is measured THROUGH the transition, and the mechanism under test is
radiative loss AT the transition -- an under-resolved transition is a
staircase, and a staircase radiates for numerical reasons. A negative
produced on an unconverged C_TRANS could not distinguish "the profile fails"
from "the staircase radiates", so this script refuses to run below the
converged floor.

SCORING IS NOT EYEBALLED. The criterion is loaded from criterion.json,
written by step 4 BEFORE any of these results existed:
    shorter than the control AND excess loss <= bar AND imbalance <= bar
    AND bar/cross phase within bar.
A candidate must pass ALL FOUR. Ranking happens only among passers, shortest
wins. Ranking on a single component is exactly the 1x4's mistake -- its sweep
ranked on balance alone and surfaced two catastrophic-totalT points as
"BEST".

LENGTHS ARE CENTRED ON THE EQUAL-PHASE PREDICTION, not on the straight
device's own mm_length. L_pi ~ W^2, so a narrowed interior images SHORTER;
solving a waisted body at the control's 28.0 um puts it past its imaging
point (measured: 23.2% phase overshoot, 2.3-2.5 dB apparent loss, phase 38
deg off quadrature -- job 12555720). The prediction is ARGUED (the 1x4's
analogous formula was refuted outright), so it CENTRES the window and the
window is what measures.

STAGED, deliberately. At C_TRANS=8 a point costs ~21 min, so a full
cartesian grid is ~14 h. Stage 1 isolates (Ww, Lt) -- the pair that sets the
sidewall angle, i.e. the mechanism under test -- at fixed Lw, over a 3-length
window. Stage 2 (different env vars, same script) refines around whatever
survives.

Run:
    python .../mmi_2x2_waist_design_sweep.py                     # stage 1
    MMI2X2_DS_WW=2.24 MMI2X2_DS_LT=4.0 MMI2X2_DS_NLEN=5 ...      # stage 2
"""
import json
import os
import time

import numpy as np

from beef import MMI, EMESolver, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT_ROOT = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "verification"))
output_dir = os.path.join(_OUT_ROOT, "mmi_2x2_waist_design")
os.makedirs(output_dir, exist_ok=True)
CRITERION_PATH = os.path.join(_OUT_ROOT, "mmi_2x2_waist_metric", "criterion.json")
CTRANS_PATH = os.path.join(_OUT_ROOT, "mmi_2x2_waist_ctrans",
                           "ctrans_convergence.json")


def _env_f(name, default):
    return float(os.environ.get(name, default))


def _env_i(name, default):
    return int(os.environ.get(name, default))


def _env_list(name, default):
    return [float(v) for v in os.environ.get(name, default).split(",")]


# --- production constants (the resolved straight 2x2) -----------------------
WL_MIN, WL_MAX = 1.530, 1.565
MM_LENGTH, MM_WIDTH = 28.0, 2.8
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 0.45, 1.3, 4.0
WAVEGUIDE_OFFSET, WAVEGUIDE_HEIGHT = 1.5, 0.22
STRAIGHT_LENGTH = 0.5
NUM_MODES = [1, 7, 15, 7, 1]
TAPER_CELLS = 15
FIXED_LENGTH = 2 * TAPER_LENGTH + 2 * STRAIGHT_LENGTH   # 9.0 um of ports

RES = _env_f("MMI2X2_DS_RES", 40.0)
N_WL = _env_i("MMI2X2_DS_NWL", 11)
OXW, OXH = 5.0, 1.5

WW_LIST = _env_list("MMI2X2_DS_WW", "2.52,2.24")
LW_LIST = _env_list("MMI2X2_DS_LW", "8.0")
LT_LIST = _env_list("MMI2X2_DS_LT", "4.0,6.0")
N_LEN = _env_i("MMI2X2_DS_NLEN", 3)
LEN_HALFSPAN = _env_f("MMI2X2_DS_LSPAN", 2.0)      # +/- um about L*

C_STUB = _env_i("MMI2X2_DS_CSTUB", 1)
C_WAIST = _env_i("MMI2X2_DS_CWAIST", 4)
C_TRANS = _env_i("MMI2X2_DS_CTRANS", 8)            # step 5a's converged floor

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

# --- gates ------------------------------------------------------------------
with open(CRITERION_PATH) as fh:
    CRIT = json.load(fh)
THR = CRIT["thresholds"]
CONTROL = CRIT["control"]

converged_floor = None
if os.path.exists(CTRANS_PATH):
    with open(CTRANS_PATH) as fh:
        converged_floor = json.load(fh).get("converged_at")
if converged_floor is not None and C_TRANS < converged_floor:
    raise SystemExit(
        f"REFUSING TO RUN: C_TRANS={C_TRANS} is below step 5a's converged "
        f"floor of {converged_floor}. Results on an unconverged transition "
        "cannot distinguish a failing profile from a radiating staircase -- "
        "which is the mechanism this whole device was chosen to test.")


def equal_phase_length(ww, lw, lt, w=MM_WIDTH, l_straight=MM_LENGTH):
    """Body length matching the straight device's accumulated imaging phase.
    ARGUED -- centres the window, does not replace measuring it."""
    target = l_straight / w ** 2
    ramp = 2.0 * (lt / (w - ww)) * (1.0 / ww - 1.0 / w) if w != ww else 0.0
    stub_total = (target - ramp - lw / ww ** 2) * w ** 2
    if stub_total < 0:
        return float("nan")
    return stub_total + 2 * lt + lw


def build(ww, lw, lt, mm_len):
    mmi = MMI(n_inputs=2, n_outputs=2, mm_length=mm_len, mm_width=MM_WIDTH,
              waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
              taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
              waveguide_height=WAVEGUIDE_HEIGHT,
              body_waist_width=ww, body_waist_length=lw,
              body_transition_length=lt)
    half = mm_len / 2
    (din0, din1), (dout0, dout1) = mmi._body_waist_bounds()
    segs = [(-half, din0, C_STUB), (din0, din1, C_TRANS),
            (din1, dout0, C_WAIST), (dout0, dout1, C_TRANS),
            (dout1, half, C_STUB)]
    body = []
    for x0, x1, n in segs:
        if n <= 0 or x1 - x0 <= 1e-12:
            continue
        body.extend([(x1 - x0) / n] * n)
    return mmi, body


def evaluate(ww, lw, lt, mm_len):
    mmi, body = build(ww, lw, lt, mm_len)
    lengths = [STRAIGHT_LENGTH, TAPER_LENGTH, mm_len, TAPER_LENGTH,
               STRAIGHT_LENGTH]
    solver = EMESolver(
        num_modes=NUM_MODES,
        num_cells=[1, TAPER_CELLS, len(body), TAPER_CELLS, 1],
        lengths=lengths, region_cell_lengths={2: body},
        oxide_width=OXW, oxide_height=OXH, res=RES,
        material_stack=SOI_220, mode_solver="tidy3d")
    c_bar = solver.solve(mmi, wavelengths, ports=["bar"],
                         input_idx=1, output_idx=1).coeff("bar")
    c_cross = solver.solve(mmi, wavelengths, ports=["cross"],
                           input_idx=1, output_idx=0).coeff("cross")
    T_bar, T_cross = np.abs(c_bar) ** 2, np.abs(c_cross) ** 2
    total = T_bar + T_cross
    dphi = np.angle(c_cross / c_bar, deg=True)
    return {
        "total_device_length": mm_len + FIXED_LENGTH,
        "excess_loss_db": float(np.max(-10 * np.log10(np.clip(total, 1e-20, None)))),
        "imbalance": float(np.max(np.abs(T_bar - T_cross))),
        "phase_deg": float(np.mean(dphi)),
        "min_total": float(np.min(total)),
    }


def score(m):
    """All four must hold. Returns (passed, per-check dict)."""
    checks = {
        "shorter": m["total_device_length"] < CONTROL["total_device_length"],
        "excess_loss": m["excess_loss_db"] <= THR["max_excess_loss_db"],
        "imbalance": m["imbalance"] <= THR["max_imbalance"],
        "phase": abs(m["phase_deg"] - THR["control_phase_deg"]) <= THR["max_phase_dev_deg"],
    }
    return all(checks.values()), checks


print("=" * 100)
print("2x2 MMI waist -- STEP 6: design sweep (scored, not eyeballed)")
print(f"  control: {CONTROL['total_device_length']:.2f} um total, "
      f"IL {CONTROL['max_excess_loss_db']:.4f} dB, "
      f"imb {CONTROL['max_imbalance']:.4f}, phase {THR['control_phase_deg']:.3f} deg")
print(f"  PASS bars: len < {CONTROL['total_device_length']:.2f} | "
      f"IL <= {THR['max_excess_loss_db']:.4f} | "
      f"imb <= {THR['max_imbalance']:.4f} | "
      f"phase within {THR['max_phase_dev_deg']:.1f} deg")
print(f"  C_TRANS={C_TRANS} (step 5a floor {converged_floor}), res={RES}, {N_WL} wl")
print(f"  grid Ww{WW_LIST} x Lw{LW_LIST} x Lt{LT_LIST}, "
      f"{N_LEN} lengths +/-{LEN_HALFSPAN} um about L*")
print("=" * 100, flush=True)

rows = []
print(f"{'Ww':>6} {'Lw':>5} {'Lt':>5} {'deg':>5} {'L*':>7} {'L':>7} {'total':>7} "
      f"{'IL':>8} {'imb':>8} {'phase':>9} {'PASS':>5} {'fails':>22} {'t':>6}")
for ww in WW_LIST:
    for lw in LW_LIST:
        for lt in LT_LIST:
            lstar = equal_phase_length(ww, lw, lt)
            if not np.isfinite(lstar):
                print(f"{ww:6.2f} {lw:5.1f} {lt:5.1f}  -- no equal-phase solution")
                continue
            deg = np.degrees(np.arctan(((MM_WIDTH - ww) / 2) / lt))
            lens = ([lstar] if N_LEN == 1 else
                    np.linspace(lstar - LEN_HALFSPAN, lstar + LEN_HALFSPAN, N_LEN))
            for L in lens:
                if 2 * lt + lw > L - 1e-9:
                    print(f"{ww:6.2f} {lw:5.1f} {lt:5.1f} {deg:5.2f} {lstar:7.3f} "
                          f"{L:7.3f}  skipped (2*Lt+Lw exceeds body)")
                    continue
                t0 = time.time()
                try:
                    m = evaluate(ww, lw, lt, float(L))
                except Exception as exc:
                    print(f"{ww:6.2f} {lw:5.1f} {lt:5.1f} {deg:5.2f} {lstar:7.3f} "
                          f"{L:7.3f}  FAILED {type(exc).__name__}: {exc}", flush=True)
                    continue
                passed, checks = score(m)
                fails = ",".join(k for k, v in checks.items() if not v) or "-"
                rows.append({"Ww": ww, "Lw": lw, "Lt": lt, "sidewall_deg": float(deg),
                             "L_star": float(lstar), "mm_length": float(L),
                             "passed": passed, "checks": checks, **m})
                print(f"{ww:6.2f} {lw:5.1f} {lt:5.1f} {deg:5.2f} {lstar:7.3f} "
                      f"{L:7.3f} {m['total_device_length']:7.2f} "
                      f"{m['excess_loss_db']:8.4f} {m['imbalance']:8.4f} "
                      f"{m['phase_deg']:9.3f} {str(passed):>5} {fails:>22} "
                      f"{time.time()-t0:5.0f}s", flush=True)

# ---------------------------------------------------------------------------
# Verdict + pre-registered stop rule
# ---------------------------------------------------------------------------
print("\n" + "=" * 100)
passers = [r for r in rows if r["passed"]]
if passers:
    best = min(passers, key=lambda r: r["total_device_length"])
    print(f"{len(passers)} of {len(rows)} points PASS all four criteria.")
    print(f"BEST (shortest passer): Ww={best['Ww']} Lw={best['Lw']} Lt={best['Lt']} "
          f"L={best['mm_length']:.3f}")
    print(f"  total device {best['total_device_length']:.2f} um vs control "
          f"{CONTROL['total_device_length']:.2f} um  "
          f"-> saves {CONTROL['total_device_length']-best['total_device_length']:.2f} um "
          f"({100*(CONTROL['total_device_length']-best['total_device_length'])/CONTROL['total_device_length']:.1f}%)")
    print(f"  IL {best['excess_loss_db']:.4f} dB, imbalance {best['imbalance']:.4f}, "
          f"phase {best['phase_deg']:.3f} deg")
    print("\n  NEXT: re-run this point at production settings before quoting it,")
    print("  and confirm it is not sitting on the edge of the length window.")
else:
    print(f"NO point of {len(rows)} passes all four criteria.")
    # stop rule component 1: within 5% of control's excess loss
    if rows:
        near = [r for r in rows
                if r["excess_loss_db"] <= 1.05 * CONTROL["max_excess_loss_db"]]
        closest = min(rows, key=lambda r: r["excess_loss_db"])
        print(f"  closest on excess loss: {closest['excess_loss_db']:.4f} dB "
              f"(control {CONTROL['max_excess_loss_db']:.4f}, "
              f"bar {THR['max_excess_loss_db']:.4f})")
        print(f"  points within 5% of control's excess loss: {len(near)}")
        if not near:
            print("\n  ** PRE-REGISTERED STOP RULE TRIGGERED (clause 1). **")
            print("  Two Ww values x their length windows produced no point within")
            print("  5% of the control's excess loss. STOP and write the negative.")
            print("  A clean negative here is a MECHANISM RESULT, not a failure:")
            print("  it refutes transition radiation as the 1x4's explanation,")
            print("  since this device removes both conditions that hypothesis")
            print("  depends on. Do not keep sweeping around it.")
        else:
            print("\n  Stop rule NOT triggered: some points are close on loss.")
            print("  A refinement pass is justified -- but bound it, and do not")
            print("  move a threshold to manufacture a pass.")

out = {"criterion_used": CRIT, "settings": {
    "res": RES, "n_wl": N_WL, "C_TRANS": C_TRANS, "C_STUB": C_STUB,
    "C_WAIST": C_WAIST, "num_modes": NUM_MODES, "taper_cells": TAPER_CELLS,
    "ctrans_converged_floor": converged_floor},
    "rows": rows,
    "n_pass": len(passers)}
path = os.path.join(output_dir, "design_sweep.json")
with open(path, "w") as fh:
    json.dump(out, fh, indent=2)
print(f"\nwritten -> {path}")
