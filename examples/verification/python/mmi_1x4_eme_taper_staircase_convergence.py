"""Residual 2 of validation_markdown/resolved/mmi_1x4_sin300_validation.md:
does the EME taper STAIRCASE discretization (C_TAPER) explain the
EME-vs-HYBRID split-shape/total-T disagreement (Residual 1)?

C_TAPER=8 is the only path EME has for approximating the port tapers --
`num_cells=[1, C_TAPER, 1, C_TAPER, 1]` in mmi_1x4_eme.py subdivides each
taper into C_TAPER *uniform-width* cells, i.e. a staircase. FDTD meshes the
taper's continuous linear width ramp exactly, with no such discretization.
Residual 2 flagged this as the LEADING suspect for Residual 1 and said the
convergence sweep that would test it was never run, at the user's
instruction. This is that sweep, run now.

WHAT THIS TESTS. Everything else pinned to mmi_1x4_eme.py's own production
config (res=40, 20x5 domain -- CASE A of the split-shape confound script,
which already ruled OUT resolution/domain as the cause of the OUTER-vs-INNER
sign flip). Only C_TAPER moves. Band centre only (N_WL=1), matching the
confound script's own cost argument.

READ IT LIKE THIS.
  - If total T trends from EME's published 0.9707 toward hybrid's 0.9454
    (and/or the split flips OUTER -> INNER) as C_TAPER grows, the staircase
    was a real contributor to Residual 1 -- restate check 1 accordingly.
  - If T and split-favour are flat across C_TAPER (within noise), the
    staircase is NOT the explanation and Residual 2's leading suspect is
    excluded, same status as mode truncation already got.
  - Either way, report where (if anywhere) T/split PLATEAU -- that is the
    number relevant to Residual 3 (EME resolution) too, since both residuals
    are about EME-side discretization.

Cost. Each taper cell is a distinct cross-section (linear width ramp), so
unlike the body's N_BODY=20 (all identical cross-sections, heavily
mode-reused) taper cells do NOT cache-reuse across each other -- cost scales
close to linearly in C_TAPER. C_TAPER=8 (published) -> 2*8*N_TAPER=160
taper mode-solves per branch; C_TAPER=64 -> 1280. Expect the top of the grid
to cost roughly 8x the published run per branch, still far cheaper than any
FDTD comparator.

Run:
    MMI1X4_CTAPER_VALUES=8,16,32,64 python .../mmi_1x4_eme_taper_staircase_convergence.py
    python .../mmi_1x4_eme_taper_staircase_convergence.py
"""
import os
import time

import numpy as np

from beef import MMI, EMESolver, SIN_300

# --- published/reference numbers this is scored against --------------------
EME_PUB = np.array([0.2431, 0.2423, 0.2423, 0.2431])          # job 11887395, C_TAPER=8
EME_PUB_TOTAL = 0.9707
HYBRID_TOTAL = 0.9454                                          # RESOLVED header, Residual 1
FDTD_REF = np.array([0.23997513751844815, 0.2424475379509313,
                     0.24244750636842585, 0.2399752894983026])  # job 12029541
FDTD_TOTAL = 0.9648

# --- device geometry, pinned to production (step 9c) -----------------------
WL = 1.55
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 1.1, 2.2, 5.0
WAVEGUIDE_OFFSET, MM_WIDTH = 3.1, 12.0
MM_LENGTH = 39.90
STRAIGHT_LENGTH = 0.5
N_BODY, N_TAPER = 20, 10
RES, OXW, OXH = 40.0, MM_WIDTH + 8.0, 5.0   # mmi_1x4_eme.py's own production config (case A)

_want = os.environ.get("MMI1X4_CTAPER_VALUES")
C_TAPER_VALUES = [int(v) for v in _want.split(",")] if _want else [8, 16, 32, 64]

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30, use_neighbor_aperture=True,
)
wavelengths = np.array([WL])


def _shape(T):
    """Which pair is favoured, and by how much. Outer = ports 0,3."""
    outer = 0.5 * (T[0] + T[3])
    inner = 0.5 * (T[1] + T[2])
    return outer, inner, ("OUTER" if outer > inner else "INNER"), abs(outer - inner)


print("=" * 78)
print("EME TAPER STAIRCASE (C_TAPER) CONVERGENCE -- band centre only")
print(f"pinned config: res={RES:.0f}, domain {OXW:.0f}x{OXH:.0f} (case A, published)")
for tag, T, tot in (("EME published (C_TAPER=8)", EME_PUB, EME_PUB_TOTAL),
                    ("FDTD reference", FDTD_REF, FDTD_TOTAL)):
    o, i, who, d = _shape(T)
    print(f"  {tag:28s} total={tot:.4f}  favours {who}  |o-i|={d:.2e}")
print(f"  {'HYBRID (Residual 1)':28s} total={HYBRID_TOTAL:.4f}  favours INNER (per RESOLVED header)")
print("=" * 78, flush=True)

rows = []
for c_taper in C_TAPER_VALUES:
    print(f"\n--- C_TAPER={c_taper} ---", flush=True)
    solver = EMESolver(
        num_modes=[1, N_TAPER, N_BODY, N_TAPER, 1],
        num_cells=[1, c_taper, 1, c_taper, 1],
        lengths=[STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, STRAIGHT_LENGTH],
        oxide_width=OXW, oxide_height=OXH, res=RES,
        material_stack=SIN_300, mode_solver="tidy3d",
    )
    t0 = time.time()
    try:
        T = np.array([
            float(np.abs(solver.solve(mmi, wavelengths, ports=[f"out{k}"],
                                      output_idx=k).transmission(port=f"out{k}"))[0])
            for k in range(4)
        ])
    except Exception as exc:
        print(f"  FAILED: {type(exc).__name__}: {exc}", flush=True)
        rows.append((c_taper, None, None, None, None))
        continue
    el = time.time() - t0
    o, i, who, d = _shape(T)
    print(f"  {el:.1f}s   T = {list(np.round(T, 6))}")
    print(f"  total={T.sum():.4f}   favours {who}   |outer-inner|={d:.2e}")
    print(f"  vs EME pub : max per-port dev = {np.max(np.abs(T - EME_PUB)):.2e}")
    print(f"  vs FDTD ref: max per-port dev = {np.max(np.abs(T - FDTD_REF)):.2e}", flush=True)
    rows.append((c_taper, T, who, d, el))

print("\n" + "=" * 78)
print(f"{'C_TAPER':>8} {'total T':>9} {'favours':>8} {'|o-i|':>9} "
      f"{'d(hybrid)':>10} {'d(FDTD)':>9} {'min':>6}")
for c_taper, T, who, d, el in rows:
    if T is None:
        print(f"{c_taper:>8} {'FAILED':>9}")
        continue
    d_hybrid = abs(float(T.sum()) - HYBRID_TOTAL)
    d_fdtd = abs(float(T.sum()) - FDTD_TOTAL)
    print(f"{c_taper:8d} {T.sum():9.4f} {who:>8} {d:9.2e} "
          f"{d_hybrid:10.4f} {d_fdtd:9.4f} {el/60:6.1f}")
print("=" * 78)
print("If total T trends 0.9707 -> 0.9454 (toward hybrid) and/or favour flips")
print("OUTER->INNER as C_TAPER grows, the staircase is a real contributor to")
print("Residual 1. If T/favour are flat, the staircase is excluded as a cause,")
print("same status mode truncation already got.")
