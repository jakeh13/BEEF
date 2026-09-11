"""A/B validation of the content-addressed mode memo for the waisted 2x2 body.

The memo (see mmi_2x2_waist_hybrid_vs_fdtd.py) reuses already-solved modes when
a cell's cross-section is byte-identical to one seen before -- which the waisted
body offers in bulk, because it is mirror-symmetric about x=0 and uniform
through the waist (40 cells -> ~13 distinct).

⚠ THIS MUST BE VALIDATED, NOT ASSUMED. The obvious version of this idea --
keying mode reuse on WIDTH -- was measured BREAKING PHYSICS on the 1x4
(totalT 0.9703 -> 0.1507, reverted the same day) because equal width does not
mean equal cross-section for an off-axis cell. The memo here keys on the actual
nx/ny/nz arrays instead, which should make reuse exact; this script is what
turns "should" into a measurement.

Solves the SAME body section twice in one process -- memo ON, then memo OFF --
and compares the resulting transfer matrices element by element. Anything above
~1e-12 means the memo is not exact and must not be used.

Deliberately small (low res, few cells, one wavelength): the question is
exactness, not accuracy, and an exactness bug shows up at any resolution.
"""
import os
import time

import numpy as np

from beef import MMI, EMESolver, SOI_220

# ⚠ meep's mode solver requires a SQUARE mesh (its resolution is isotropic),
# so the oxide box must be an integer number of pixels in BOTH directions:
# oxide_width*res AND oxide_height*res must be integers. With
# oxide_height=1.5 that means res must be a multiple of 2 -- res=15 gives
# 1.5*15 = 22.5 px, which rounds to 22 and makes dy != dx:
#   ValueError: compute_modes_meep requires a square mesh ...
#              got dx=0.0667, dy=0.0682
# res 20: 5.0*20=100, 1.5*20=30 -> dx = dy = 0.05. res 40 (production) also OK.
RES = float(os.environ.get("MEMO_RES", 20.0))
CELLS = int(os.environ.get("MEMO_CELLS", 20))
WL = float(os.environ.get("MEMO_WL", 1.5475))
N_MODES = int(os.environ.get("MEMO_NMODES", 8))
WW, LW, LT, L = 2.24, 6.0, 6.0, 21.5
BODY_STUB = 0.25

mmi = MMI(n_inputs=2, n_outputs=2, mm_length=L, mm_width=2.8,
          waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
          waveguide_offset=1.5, waveguide_height=0.22,
          body_waist_width=WW, body_waist_length=LW,
          body_transition_length=LT)
half_body = mmi.mm_length / 2
x0, x1 = -half_body + BODY_STUB, half_body - BODY_STUB

print("=" * 78)
print("MODE MEMO A/B -- waisted 2x2 body section")
print(f"  body [{x0:.3f}, {x1:.3f}] um in {CELLS} cells, res {RES}, "
      f"wl {WL} um, {N_MODES} modes")
print(f"  widths: start {mmi.width_at(x0):.3f}  mid {mmi.width_at(0.0):.3f}  "
      f"end {mmi.width_at(x1):.3f} um")
print("=" * 78, flush=True)

_orig = EMESolver.compute_modes
STATS = {"hits": 0, "misses": 0}
MEMO = {}
MEMO_ON = [False]


def _install():
    import functools
    import hashlib

    def _key(cs, n):
        h = hashlib.blake2b(digest_size=16)
        for comp in ("nx", "ny", "nz"):
            a = np.ascontiguousarray(getattr(cs, comp))
            h.update(np.asarray(a.shape, dtype=np.int64).tobytes())
            h.update(a.tobytes())
        return (h.hexdigest(), float(cs.env.wl), float(cs.env.T), int(n))

    @functools.wraps(_orig)
    def compute_modes(self, cs, n_modes, *a, **kw):
        if not MEMO_ON[0]:
            return _orig(self, cs, n_modes, *a, **kw)
        k = _key(cs, n_modes)
        if k in MEMO:
            STATS["hits"] += 1
            return MEMO[k]
        STATS["misses"] += 1
        r = _orig(self, cs, n_modes, *a, **kw)
        MEMO[k] = r
        return r

    EMESolver.compute_modes = compute_modes


_install()


def run(memo_on):
    MEMO.clear()
    STATS["hits"] = STATS["misses"] = 0
    MEMO_ON[0] = memo_on
    solver = EMESolver(num_modes=[N_MODES], num_cells=[CELLS],
                       lengths=[x1 - x0], oxide_width=5.0, oxide_height=1.5,
                       res=RES, material_stack=SOI_220, mode_solver="meep")
    t0 = time.time()
    mat = solver.solve_section(mmi, np.array([WL]), x0, x1,
                               n_modes_in=N_MODES, n_modes_out=N_MODES)
    dt = time.time() - t0
    return np.asarray(mat), dt, dict(STATS)


print("[A] memo OFF ...", flush=True)
m_off, t_off, s_off = run(False)
print(f"[A] {t_off:.0f}s, {s_off['misses']} solves\n", flush=True)

print("[B] memo ON ...", flush=True)
m_on, t_on, s_on = run(True)
tot = s_on["hits"] + s_on["misses"]
print(f"[B] {t_on:.0f}s, {s_on['misses']} solved + {s_on['hits']} reused "
      f"of {tot} ({100*s_on['hits']/max(tot,1):.0f}% avoided)\n", flush=True)

d = np.abs(m_on - m_off)
rel = d / np.maximum(np.abs(m_off), 1e-30)
print(f"{'max |dS|':>14} {d.max():.3e}")
print(f"{'mean |dS|':>14} {d.mean():.3e}")
print(f"{'max rel':>14} {rel.max():.3e}")
print(f"{'speedup':>14} {t_off/max(t_on,1e-9):.2f}x")
ok = d.max() < 1e-12
print(f"\nVERDICT: {'EXACT -- memo is safe to use' if ok else 'NOT EXACT -- DO NOT USE THE MEMO'}")
if not ok:
    print("  The reuse changed the answer. Do not enable MMI2X2_HY_MEMO.")
raise SystemExit(0 if ok else 1)
