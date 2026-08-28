"""Tapered-body 1x4: find a real design point, and MEASURE the length payoff.

validation_markdown/mmi_1x4_tapered_body_validation.md step 6 (and the
"Length payoff" section, which is ARGUED and has never been simulated).

WHAT THIS REPLACES. The tapered device currently has NO design: the
body_taper_length=8.0 / body_start_width=6.0 pair used for the field plots is
a placeholder with no sweep behind it, and the payoff estimate is a 1-D
self-imaging scaling argument:

    L_tapered = L_straight - Lt * (W/Ws - 1)          [ARGUED, not simulated]

from equating the accumulated 1/W^2 imaging budget, exactly integrable for a
linear taper. It predicts ~20-30% shorter at Ws=5-6 um, Lt=8-12 um. It ignores
taper RADIATION LOSS and imaging QUALITY entirely -- the two terms that punish
an aggressive taper -- so it says the premise is worth testing, not that it
works.

GEOMETRY CONSTRAINT (measured, not assumed). The four output tapers need 11.5
of the 12.0 um output face; the single on-axis input port needs 2.2. So the
body taper is ASYMMETRIC: body_start_width is swept, body_end_width stays at
mm_width. All of the saving comes from the input side.

WHAT IS SWEPT, AND WHY mm_length IS SWEPT WITH IT. A tapered body changes the
self-imaging length, so comparing profiles at a FIXED mm_length would measure
"how detuned is this profile from the straight device's length", not the
profile's own merit. For each (Ws, Lt) the script therefore sweeps mm_length
around that profile's own PREDICTED optimum and reports the best it finds.
If the prediction is wrong the sweep still finds the real optimum, provided the
window is wide enough -- LEN_FRAC controls that, and the script WARNS if a
profile's best lands on a window edge (i.e. the optimum may be outside).

WHY NOT solve_length_sweep(). Its analytic phase trick is uniform-section-only
(step 3's own note), and on the whole-device path sweeping the body region's
length would rescale ALL its cells -- including the fixed-length taper bands.
Explicit per-length solve() calls keep the taper bands fixed at Lt, which is
the geometry actually intended.

THE STRAIGHT CONTROL IS RUN AT THE SAME DISCRETISATION. A straight body is
included with the SAME BODY_CELLS as every tapered case, so the comparison is
geometry-only. Drawing the straight baseline with num_cells=1 and the tapered
ones with num_cells=N would confound taper physics with cell-count -- the same
trap the straight-vs-tapered field-plot script already avoids.

METRICS. Per (profile, length): total T (captures taper radiation loss) and
BALANCE = max-min across the four outputs (captures imaging quality). A short
device that splits badly is not a win, so both are reported and the "best"
length per profile is chosen on balance, with total T shown alongside.

COST. Deliberately run at a DESIGN resolution (RES, default 30) and the 16x4
domain, not the EME production 40/20x5: this is a relative comparison between
profiles, and the confound grid (jobs 12323220/12331012) measured 760 s per
whole-device solve at 30/16x4 versus 3098 s at 40/20x5. The winning profile
should be re-run at production settings before being trusted as a number.
Smoke first -- the batch file does one profile.

Run:
    MMI1X4_TP_WS=6.0 MMI1X4_TP_LT=8.0 MMI1X4_TP_NLEN=3 python .../mmi_1x4_taper_profile_sweep.py
    python .../mmi_1x4_taper_profile_sweep.py
"""
import os
import time

import numpy as np

from beef import MMI, EMESolver, SIN_300

# --- production constants (step 9c) ----------------------------------------
WL = 1.55
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 1.1, 2.2, 5.0
WAVEGUIDE_OFFSET, MM_WIDTH = 3.1, 12.0
L_STRAIGHT = 39.90          # the straight device's settled mm_length
STRAIGHT_LENGTH = 0.5
N_BODY, N_TAPER, C_TAPER = 20, 10, 8

# --- sweep config -----------------------------------------------------------
RES = float(os.environ.get("MMI1X4_TP_RES", 30.0))
OXW = float(os.environ.get("MMI1X4_TP_OXW", MM_WIDTH + 4.0))
OXH = float(os.environ.get("MMI1X4_TP_OXH", 4.0))
BODY_CELLS = int(os.environ.get("MMI1X4_TP_BODYCELLS", 16))
_ws = os.environ.get("MMI1X4_TP_WS")
_lt = os.environ.get("MMI1X4_TP_LT")
WS_LIST = [float(v) for v in _ws.split(",")] if _ws else [8.0, 6.0, 5.0]
LT_LIST = [float(v) for v in _lt.split(",")] if _lt else [8.0, 12.0]
N_LEN = int(os.environ.get("MMI1X4_TP_NLEN", 5))     # lengths per profile
LEN_FRAC = float(os.environ.get("MMI1X4_TP_LENFRAC", 0.20))  # +/- window


# Sweep-centre override. Added 2026-08-27 after smoke job 12383125 REFUTED the
# ARGUED formula below: at its predicted L=31.90 (Ws=6, Lt=8) the device gave
# totalT=0.5544 / balance=4.98e-02 against the straight control's 0.9703 /
# 2.38e-03 -- a 44% power loss, i.e. imaging FAILURE, not taper radiation. The
# best point in that window also sat on the UPPER edge (38.28), so the true
# optimum is at or beyond the straight device's own 39.90 um. Centring the
# window on the formula would keep sweeping the wrong region.
L_CENTRE = os.environ.get("MMI1X4_TP_LCENTRE")


def predicted_length(ws, lt):
    """Sweep centre. MMI1X4_TP_LCENTRE overrides; otherwise the ARGUED
    L = L_straight - Lt*(W/Ws - 1), which smoke 12383125 showed is wrong for
    this device -- kept only so the refuted prediction stays reproducible
    alongside whatever the measurement says."""
    if L_CENTRE is not None:
        return float(L_CENTRE)
    return L_STRAIGHT - lt * (MM_WIDTH / ws - 1.0)


def build(mm_length, ws=None, lt=0.0):
    kw = {}
    if lt > 0 and ws is not None:
        kw = dict(body_taper_length=lt, body_start_width=ws)
    return MMI(
        n_inputs=1, n_outputs=4, mm_length=mm_length, mm_width=MM_WIDTH,
        waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
        taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
        waveguide_height=0.30, use_neighbor_aperture=True, **kw
    )


def solve(mmi):
    """Total T and balance for one device. Straight control uses the SAME
    BODY_CELLS as the tapered cases -- see module docstring."""
    solver = EMESolver(
        num_modes=[1, N_TAPER, N_BODY, N_TAPER, 1],
        num_cells=[1, C_TAPER, BODY_CELLS, C_TAPER, 1],
        lengths=[STRAIGHT_LENGTH, TAPER_LENGTH, mmi.mm_length,
                 TAPER_LENGTH, STRAIGHT_LENGTH],
        oxide_width=OXW, oxide_height=OXH, res=RES,
        material_stack=SIN_300, mode_solver="tidy3d",
    )
    T = np.array([
        float(solver.solve(mmi, np.array([WL]), ports=[f"out{k}"],
                           output_idx=k).transmission(port=f"out{k}")[0])
        for k in range(4)
    ])
    return T


print("=" * 84)
print("TAPERED 1x4 -- taper-profile sweep (design point + MEASURED length payoff)")
print(f"  res={RES}  domain {OXW:.0f}x{OXH:.0f}  BODY_CELLS={BODY_CELLS}  wl={WL}")
print(f"  Ws list {WS_LIST}   Lt list {LT_LIST}   {N_LEN} lengths/profile "
      f"(+/-{100*LEN_FRAC:.0f}% about prediction)")
print(f"  output face pinned at mm_width={MM_WIDTH} (fan-out needs 11.5)")
print("=" * 84, flush=True)

# --- straight control, same discretisation ---------------------------------
t0 = time.time()
Ts = solve(build(L_STRAIGHT))
bal_s = Ts.max() - Ts.min()
print(f"\nSTRAIGHT control  L={L_STRAIGHT:.2f}  totalT={Ts.sum():.4f}  "
      f"balance={bal_s:.2e}   ({time.time()-t0:.0f}s)")
print(f"  T = {list(np.round(Ts, 5))}", flush=True)

rows = []
for ws in WS_LIST:
    for lt in LT_LIST:
        Lp = predicted_length(ws, lt)
        if Lp <= 2 * lt + 1.0:
            print(f"\n-- SKIP Ws={ws} Lt={lt}: predicted L={Lp:.2f} leaves no "
                  f"uniform middle (needs > 2*Lt={2*lt:.1f})", flush=True)
            continue
        lens = np.linspace(Lp * (1 - LEN_FRAC), Lp * (1 + LEN_FRAC), N_LEN)
        print(f"\n-- Ws={ws:.1f} Lt={lt:.1f}  predicted L={Lp:.2f}  "
              f"sweeping {lens[0]:.2f}..{lens[-1]:.2f}", flush=True)
        best = None
        for L in lens:
            if L <= 2 * lt + 1.0:
                print(f"     L={L:6.2f}  skipped (no uniform middle)")
                continue
            t0 = time.time()
            try:
                T = solve(build(L, ws=ws, lt=lt))
            except Exception as exc:
                print(f"     L={L:6.2f}  FAILED {type(exc).__name__}: {exc}", flush=True)
                continue
            bal = T.max() - T.min()
            print(f"     L={L:6.2f}  totalT={T.sum():.4f}  balance={bal:.2e}  "
                  f"({time.time()-t0:.0f}s)", flush=True)
            if best is None or bal < best[2]:
                best = (L, T.sum(), bal, T)
        if best is None:
            continue
        L, tot, bal, T = best
        edge = abs(L - lens[0]) < 1e-9 or abs(L - lens[-1]) < 1e-9
        rows.append((ws, lt, Lp, L, tot, bal, edge))
        print(f"   BEST L={L:.2f}  totalT={tot:.4f}  balance={bal:.2e}"
              + ("   <-- ON WINDOW EDGE, true optimum may lie outside" if edge else ""),
              flush=True)

print("\n" + "=" * 84)
print(f"{'Ws':>5} {'Lt':>5} {'L_pred':>7} {'L_best':>7} {'totalT':>8} {'balance':>9} "
      f"{'saving':>8} {'vs straight':>12}")
print(f"{'--':>5} {'--':>5} {'--':>7} {L_STRAIGHT:7.2f} {Ts.sum():8.4f} {bal_s:9.2e} "
      f"{'--':>8} {'(control)':>12}")
for ws, lt, Lp, L, tot, bal, edge in rows:
    sav = L_STRAIGHT - L
    print(f"{ws:5.1f} {lt:5.1f} {Lp:7.2f} {L:7.2f} {tot:8.4f} {bal:9.2e} "
          f"{sav:7.2f}um {100*sav/L_STRAIGHT:6.1f}%" + ("  *edge" if edge else ""))
print("=" * 84)
print("A profile only wins if it is SHORTER *and* holds totalT/balance against")
print("the straight control. Re-run any winner at production settings")
print("(res=40, 20x5) before quoting its numbers -- this grid is a DESIGN sweep.")
