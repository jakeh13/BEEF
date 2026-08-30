"""2x2 MMI waist body -- STEP 5a: transition slice count (C_TRANS)
convergence. GATING: step 6 must not run until this passes.

validation_markdown/mmi_2x2_waist_body_validation.md step 5a.

WHY THIS GATES EVERYTHING. The mechanism under test on this device is
RADIATIVE LOSS AT THE WIDTH TRANSITIONS (the leading explanation for the
1x4's waist failure -- see that doc's re-closing note). An under-resolved
transition is a STAIRCASE, and a staircase radiates for numerical reasons
that have nothing to do with the real sidewall. Measuring transition loss on
an unconverged transition would therefore CONFIRM THE HYPOTHESIS FOR THE
WRONG REASON -- the exact trap mmi_1x4_sin300_validation.md records for the
EME taper staircase in its own check 1.

The 1x4 already measured this axis mattering: doubling C_TRANS 6->12 visibly
changed the field character (job 12446756), turning a single spurious refocus
into a progressive multi-lobe pattern. That is a large qualitative change from
a pure discretisation parameter. Note what it does and does not establish: it
shows NON-convergence at 6, not convergence at 12. This script settles the
question properly for the 2x2 rather than inheriting a number.

MEASURED AT THE STEEPEST INTENDED PROFILE, ON PURPOSE. Convergence is
measured at the most aggressive (Ww) and shortest (Lt) profile step 6 intends
to explore, i.e. the steepest sidewall in the grid. A count that resolves the
steepest ramp also resolves every gentler one, so the answer is a
CONSERVATIVE floor for the whole grid rather than a value tied to one gentle
point. If step 6 ever explores steeper than this, re-run 5a at that edge.

CONVERGENCE RULE, PRE-REGISTERED. Converged = successive change below
tolerance for TWO CONSECUTIVE refinements (one plateau can be coincidence).
Tolerances are 1/10 of step 4's own design margins, so discretisation error
is an order of magnitude below the differences the design sweep must resolve:

    excess loss   0.010 dB   (design margin 0.10 dB)
    imbalance     0.0069     (design margin 0.0686, the control's own)
    phase         0.20 deg   (design margin 2.0 deg)

Scored on step 4's FULL metric, not totalT -- a staircase can conserve power
while rotating phase, and this device's observable is a transfer matrix.

Run:
    MMI2X2_CT_NWL=3 MMI2X2_CT_LIST=2,4,8 python .../mmi_2x2_waist_ctrans_convergence.py
    python .../mmi_2x2_waist_ctrans_convergence.py      # real sweep
"""
import json
import os
import time

import numpy as np

from beef import MMI, EMESolver, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT_ROOT = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "verification"))
output_dir = os.path.join(_OUT_ROOT, "mmi_2x2_waist_ctrans")
os.makedirs(output_dir, exist_ok=True)
CRITERION_PATH = os.path.join(_OUT_ROOT, "mmi_2x2_waist_metric", "criterion.json")


def _env_f(name, default):
    return float(os.environ.get(name, default))


def _env_i(name, default):
    return int(os.environ.get(name, default))


# --- production constants (the resolved straight 2x2) -----------------------
WL_MIN, WL_MAX = 1.530, 1.565
MM_LENGTH, MM_WIDTH = 28.0, 2.8
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 0.45, 1.3, 4.0
WAVEGUIDE_OFFSET, WAVEGUIDE_HEIGHT = 1.5, 0.22
STRAIGHT_LENGTH = 0.5
NUM_MODES = [1, 7, 15, 7, 1]        # per REGION -- not the hybrid's N_MODES=14
TAPER_CELLS = 15                    # port tapers; 1 costs ~2 dB of staircase

RES = _env_f("MMI2X2_CT_RES", 40.0)
N_WL = _env_i("MMI2X2_CT_NWL", 11)
OXW, OXH = 5.0, 1.5

# Steepest profile step 6 intends to explore -- see module docstring.
WW = _env_f("MMI2X2_CT_WW", 2.24)          # most aggressive intended (80%)
LW = _env_f("MMI2X2_CT_LW", 8.0)
LT = _env_f("MMI2X2_CT_LT", 4.0)           # shortest intended transition
C_STUB = _env_i("MMI2X2_CT_CSTUB", 1)
C_WAIST = _env_i("MMI2X2_CT_CWAIST", 4)
C_TRANS_LIST = [int(v) for v in
                os.environ.get("MMI2X2_CT_LIST", "2,4,6,8,12,16,24").split(",")]

# Convergence tolerances: 1/10 of step 4's design margins.
TOL = {"excess_loss_db": 0.010, "imbalance": 0.0069, "phase_deg": 0.20}
N_CONSECUTIVE = 2

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
sidewall_deg = np.degrees(np.arctan(((MM_WIDTH - WW) / 2) / LT))


def equal_phase_length(ww, lw, lt, l_straight=MM_LENGTH, w=MM_WIDTH):
    """Body length whose accumulated imaging phase matches the straight
    device's, for this waist profile. ARGUED, see below.

    ⚠ WHY THIS MATTERS FOR A *CONVERGENCE* STUDY, not just for step 6.
    L_pi ~ W^2, so narrowing the interior SHORTENS the imaging length. Run at
    the straight device's own mm_length, a waisted body is well past its
    self-imaging point -- measured for this profile: the phase integral
    overshoots by 23.2% at L=28.0, and the first smoke run (job 12555720)
    duly reported 2.3-2.5 dB excess loss and a phase 38 deg off quadrature.
    Those numbers say nothing about the waist's merit; they say the device
    was being solved 6.5 um past where it images. Converging the
    discretisation at an operating point step 6 will never use is not
    representative, so this auto-centres the study instead.

    Equal accumulated phase, integrating 1/W(x)^2 over the 5 segments:
        L_straight / W^2 == integral dx / W(x)^2
    The two transition bands integrate analytically for a linear ramp.

    ⚠ ARGUED, NOT MEASURED. The 1x4's directly analogous argued formula was
    REFUTED outright (at its predicted length that device lost 44% of its
    power). This is used ONLY to place a representative operating point and
    to centre step 6's window -- never as a substitute for measuring.
    """
    target = l_straight / w ** 2
    ramp = 2.0 * (lt / (w - ww)) * (1.0 / ww - 1.0 / w) if w != ww else 0.0
    fixed = ramp + lw / ww ** 2
    # remaining phase must come from the two stubs, both at full width w
    stub_total = (target - fixed) * w ** 2
    if stub_total < 0:
        return float("nan")
    return stub_total + 2 * lt + lw


L_PRED = equal_phase_length(WW, LW, LT)
MM_LENGTH_USED = _env_f("MMI2X2_CT_L", L_PRED if np.isfinite(L_PRED) else MM_LENGTH)
lengths = [STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH_USED, TAPER_LENGTH,
           STRAIGHT_LENGTH]

print("=" * 84)
print("2x2 MMI waist -- STEP 5a: C_TRANS convergence (GATING for step 6)")
print(f"  profile Ww={WW} Lw={LW} Lt={LT}  -> sidewall {sidewall_deg:.2f} deg "
      f"(steepest in the intended grid)")
print(f"  mm_length {MM_LENGTH_USED:.3f} um  (equal-phase prediction "
      f"{L_PRED:.3f}; straight control is {MM_LENGTH})")
print(f"  -- ARGUED centring only, so convergence is measured near where step 6")
print(f"     will actually operate rather than 6.5 um past the imaging point")
print(f"  res={RES}  {N_WL} wl  C_STUB={C_STUB} C_WAIST={C_WAIST} (held fixed)")
print(f"  sweeping C_TRANS {C_TRANS_LIST}")
print(f"  converged = {N_CONSECUTIVE} consecutive steps under "
      f"{TOL['excess_loss_db']} dB / {TOL['imbalance']:.4f} / {TOL['phase_deg']} deg")
print("=" * 84, flush=True)

if os.path.exists(CRITERION_PATH):
    with open(CRITERION_PATH) as fh:
        criterion = json.load(fh)
    print(f"loaded step 4 criterion: control excess loss "
          f"{criterion['control']['max_excess_loss_db']:.4f} dB, "
          f"imbalance {criterion['control']['max_imbalance']:.4f}\n", flush=True)
else:
    criterion = None
    print(f"WARNING: no criterion at {CRITERION_PATH} -- run step 4 first.\n")


def build(c_trans):
    """Waisted MMI plus the segment-aligned body subdivision for this C_TRANS.

    Uses region_cell_lengths (added 2026-08-29) so cell boundaries land on the
    waist's internal edges. An even division cannot do that, and a cell that
    straddles a segment boundary sits on no single cross-section -- which is
    exactly the artefact this convergence study would otherwise be measuring.
    """
    mmi = MMI(n_inputs=2, n_outputs=2, mm_length=MM_LENGTH_USED, mm_width=MM_WIDTH,
              waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
              taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
              waveguide_height=WAVEGUIDE_HEIGHT,
              body_waist_width=WW, body_waist_length=LW,
              body_transition_length=LT)
    half = MM_LENGTH_USED / 2
    (din0, din1), (dout0, dout1) = mmi._body_waist_bounds()
    segs = [(-half, din0, C_STUB), (din0, din1, c_trans),
            (din1, dout0, C_WAIST), (dout0, dout1, c_trans),
            (dout1, half, C_STUB)]
    body = []
    for x0, x1, n in segs:
        if n <= 0 or x1 - x0 <= 1e-12:
            continue
        body.extend([(x1 - x0) / n] * n)
    return mmi, body


def solve_metrics(mmi, body_cells):
    solver = EMESolver(
        num_modes=NUM_MODES,
        num_cells=[1, TAPER_CELLS, len(body_cells), TAPER_CELLS, 1],
        lengths=lengths, region_cell_lengths={2: body_cells},
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
        "excess_loss_db": float(np.max(-10 * np.log10(np.clip(total, 1e-20, None)))),
        "imbalance": float(np.max(np.abs(T_bar - T_cross))),
        "phase_deg": float(np.mean(dphi)),
        "min_total": float(np.min(total)),
    }


rows = []
print(f"{'C_TRANS':>8} {'cells':>6} {'IL(dB)':>9} {'imbal':>9} {'phase':>9} "
      f"{'minT':>8} {'dIL':>9} {'dImb':>9} {'dPh':>8} {'time':>7}")
for c in C_TRANS_LIST:
    mmi, body = build(c)
    t0 = time.time()
    try:
        m = solve_metrics(mmi, body)
    except Exception as exc:
        print(f"{c:8d}  FAILED {type(exc).__name__}: {exc}", flush=True)
        continue
    if rows:
        p = rows[-1][1]
        d_il = abs(m["excess_loss_db"] - p["excess_loss_db"])
        d_im = abs(m["imbalance"] - p["imbalance"])
        d_ph = abs(m["phase_deg"] - p["phase_deg"])
    else:
        d_il = d_im = d_ph = float("nan")
    rows.append((c, m, (d_il, d_im, d_ph)))
    print(f"{c:8d} {len(body):6d} {m['excess_loss_db']:9.5f} "
          f"{m['imbalance']:9.5f} {m['phase_deg']:9.4f} {m['min_total']:8.5f} "
          f"{d_il:9.2e} {d_im:9.2e} {d_ph:8.3f} {time.time()-t0:6.0f}s",
          flush=True)

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------
print("\n" + "=" * 84)
converged_at = None
streak = 0
for c, m, (d_il, d_im, d_ph) in rows[1:]:
    ok = (d_il < TOL["excess_loss_db"] and d_im < TOL["imbalance"]
          and d_ph < TOL["phase_deg"])
    streak = streak + 1 if ok else 0
    if streak >= N_CONSECUTIVE and converged_at is None:
        converged_at = c

if converged_at is not None:
    print(f"CONVERGED at C_TRANS = {converged_at} "
          f"({N_CONSECUTIVE} consecutive steps inside tolerance).")
    print(f"  Use C_TRANS >= {converged_at} for step 6.")
    print(f"  ⚠ This is a FLOOR tied to the profile measured here "
          f"(Ww={WW}, Lt={LT}, sidewall {sidewall_deg:.2f} deg).")
    print("    It is conservative for GENTLER profiles. If step 6 explores a "
          "STEEPER\n    ramp than this, re-run 5a at that edge.")
else:
    print("NOT CONVERGED within the swept range.")
    print("  Step 6 is BLOCKED. Extend C_TRANS upward and re-run -- do not")
    print("  proceed on an unconverged transition: the mechanism under test")
    print("  (transition radiation) is exactly what a staircase fakes.")

out = {
    "profile": {"Ww": WW, "Lw": LW, "Lt": LT, "sidewall_deg": float(sidewall_deg),
                "C_STUB": C_STUB, "C_WAIST": C_WAIST,
                "mm_length_used": MM_LENGTH_USED,
                "equal_phase_prediction": float(L_PRED),
                "straight_control_mm_length": MM_LENGTH},
    "settings": {"res": RES, "n_wl": N_WL, "num_modes": NUM_MODES,
                 "taper_cells": TAPER_CELLS},
    "tolerances": TOL, "n_consecutive": N_CONSECUTIVE,
    "rows": [{"C_TRANS": c, **m} for c, m, _ in rows],
    "converged_at": converged_at,
}
path = os.path.join(output_dir, "ctrans_convergence.json")
with open(path, "w") as fh:
    json.dump(out, fh, indent=2)
print(f"\nwritten -> {path}")
