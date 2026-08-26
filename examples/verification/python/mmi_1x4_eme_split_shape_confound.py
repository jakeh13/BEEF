"""Why does the whole-device EME disagree with the whole-device FDTD on the
SPLIT SHAPE? Decompose the confound: engine vs resolution vs domain size.

validation_markdown/mmi_1x4_sin300_validation.md step 8 check 1 / check 3.

THE PROBLEM WITH THE COMPARISON AS IT STANDS. Check 1 was reported as
"EME favours the OUTER port pair, the whole-device FDTD reference favours the
INNER pair, so EME has the sign wrong." Those two runs do share the same
DEVICE geometry -- but NOT the same simulation parameters:

    quantity          EME prod (11887395)   FDTD ref (12029541)
    resolution                        40                     30
    domain (y x z)                20 x 5                 16 x 4
    body modes                        20      none (whole-device)

So three things differ at once, and attributing the split-shape sign to the
ENGINE is not supported by those two runs alone. (The domain difference is
not incidental: 16 x 4 was adopted 2026-08-19 on a convergence check whose
own conclusion was "modes 19/22 never converge but are IRRELEVANT because
production keeps N_MODES=12" -- true for the hybrid's 12-mode interfaces, but
mmi_1x4_eme.py carries N_BODY=20, so it is the one caller that check does not
cover. EME staying at the older, larger 20 x 5 is arguably the more
appropriate choice for 20 modes, not an oversight.)

Check 3 (job 12275175) then found the EME/tidy3d DIFFERENTIAL (Delta n_eff,
the quantity MMI self-imaging actually depends on) does NOT converge with
resolution the way meep's does -- 2.17e-03 drift 30->40, no smaller than its
absolute drift, against meep's 1.99e-04. Split shape is a
Delta-beta-across-many-modes quantity, so an under-converged EME differential
is a live candidate cause. Hypothesis, not established.

WHAT THIS SCRIPT DOES. Whole-device EME at BAND CENTRE ONLY (N_WL=1, matching
the FDTD reference's own single wavelength), over a small grid of (res,
domain) so each variable moves alone:

    A  res=40, 20x5   reproduce mmi_1x4_eme.py's published numbers
    B  res=30, 16x4   match the FDTD reference's OWN parameters
    C  res=40, 16x4   isolate DOMAIN (vs A)
    D  res=50, 20x5   isolate RESOLUTION upward (vs A)
    E  res=60, 20x5   is the split shape still moving at all?

Scored against both published results:
    EME published   T = [0.2431, 0.2423, 0.2423, 0.2431]  total 0.9707  (OUTER)
    FDTD reference  T = [0.23998, 0.24245, 0.24245, 0.23998] total 0.9648 (INNER)

READ IT LIKE THIS. If A reproduces the published EME split and B/C/D/E move
it toward INNER, the sign was a resolution/domain artefact, not an engine
disagreement -- and check 1's conclusion needs restating. If the split stays
OUTER at res=60, the disagreement is real and the engine is the difference.
Either outcome is worth having; the current state (three variables moving at
once) supports neither.

Cost: band centre only, so each config is ~182 mode-solves per output branch
x 4 branches. Cheap relative to anything FDTD. Also benefits from the
2026-08-25 lazy duplicate-guard change in EMESolver.compute_modes() (solves
num_modes instead of num_modes+8 unless a duplicate actually appears), so
these numbers come from AFTER that change -- it is a pure cost optimisation
with measured max |dn_eff| = 0, so it does not affect comparability with the
pre-change published runs above.

Run:
    MMI1X4_SPLIT_CASES=A,B python .../mmi_1x4_eme_split_shape_confound.py
    python .../mmi_1x4_eme_split_shape_confound.py
"""
import os
import time

import numpy as np

from beef import MMI, EMESolver, SIN_300

# --- published results this is scored against ------------------------------
EME_PUB = np.array([0.2431, 0.2423, 0.2423, 0.2431])          # job 11887395, res40/20x5
FDTD_REF = np.array([0.23997513751844815, 0.2424475379509313,
                     0.24244750636842585, 0.2399752894983026])  # job 12029541, res30/16x4

# --- device geometry, pinned to production (step 9c) -----------------------
WL = 1.55
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 1.1, 2.2, 5.0
WAVEGUIDE_OFFSET, MM_WIDTH = 3.1, 12.0
MM_LENGTH = 39.90
STRAIGHT_LENGTH = 0.5
N_BODY, N_TAPER, C_TAPER = 20, 10, 8

# --- the (res, oxide_width, oxide_height) grid -----------------------------
CASES = {
    "A": (40.0, MM_WIDTH + 8.0, 5.0),   # mmi_1x4_eme.py as published
    "B": (30.0, MM_WIDTH + 4.0, 4.0),   # the FDTD reference's own parameters
    "C": (40.0, MM_WIDTH + 4.0, 4.0),   # domain isolated (vs A)
    "D": (50.0, MM_WIDTH + 8.0, 5.0),   # resolution isolated (vs A)
    "E": (60.0, MM_WIDTH + 8.0, 5.0),   # still moving?
}
_want = os.environ.get("MMI1X4_SPLIT_CASES")
ORDER = [c.strip() for c in _want.split(",")] if _want else ["A", "B", "C", "D", "E"]

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
print("EME SPLIT-SHAPE CONFOUND DECOMPOSITION -- band centre only")
for tag, ref in (("EME published (res40, 20x5)", EME_PUB),
                 ("FDTD reference (res30, 16x4)", FDTD_REF)):
    o, i, who, d = _shape(ref)
    print(f"  {tag:30s} total={ref.sum():.4f}  favours {who}  |o-i|={d:.2e}")
print("=" * 78, flush=True)

rows = []
for tag in ORDER:
    res, oxw, oxh = CASES[tag]
    print(f"\n--- case {tag}: res={res:.0f}, domain {oxw:.0f} x {oxh:.0f} ---", flush=True)
    solver = EMESolver(
        num_modes=[1, N_TAPER, N_BODY, N_TAPER, 1],
        num_cells=[1, C_TAPER, 1, C_TAPER, 1],
        lengths=[STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, STRAIGHT_LENGTH],
        oxide_width=oxw, oxide_height=oxh, res=res,
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
        rows.append((tag, res, oxw, oxh, None, None, None, None))
        continue
    el = time.time() - t0
    o, i, who, d = _shape(T)
    print(f"  {el:.1f}s   T = {list(np.round(T, 6))}")
    print(f"  total={T.sum():.4f}   favours {who}   |outer-inner|={d:.2e}")
    print(f"  vs FDTD ref: max per-port dev = {np.max(np.abs(T - FDTD_REF)):.2e}")
    print(f"  vs EME pub : max per-port dev = {np.max(np.abs(T - EME_PUB)):.2e}", flush=True)
    rows.append((tag, res, oxw, oxh, T, who, d, el))

print("\n" + "=" * 78)
print(f"{'case':>5} {'res':>5} {'domain':>9} {'total T':>9} {'favours':>8} "
      f"{'|o-i|':>9} {'devFDTD':>9} {'hrs':>6}")
for tag, res, oxw, oxh, T, who, d, el in rows:
    if T is None:
        print(f"{tag:>5} {res:5.0f} {f'{oxw:.0f}x{oxh:.0f}':>9} {'FAILED':>9}")
        continue
    print(f"{tag:>5} {res:5.0f} {f'{oxw:.0f}x{oxh:.0f}':>9} {T.sum():9.4f} {who:>8} "
          f"{d:9.2e} {np.max(np.abs(T-FDTD_REF)):9.2e} {el/3600:6.2f}")
print("=" * 78)
print("If A reproduces OUTER and B/C/D/E trend to INNER, the split-shape sign")
print("was a resolution/domain artefact and check 1 needs restating. If it")
print("stays OUTER at res=60, the engine disagreement is real.")
