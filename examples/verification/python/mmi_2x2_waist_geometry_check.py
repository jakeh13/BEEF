"""2x2 MMI with a WAIST body -- steps 1 and 2 of
validation_markdown/mmi_2x2_waist_body_validation.md.

STEP 1 (geometry sanity). Renders the waisted 2x2 next to the straight
control from ONE config, so the only difference between the two pictures is
the body profile, and prints the width profile at the five segment
boundaries. Also re-derives the port-clearance constraint from output_ys
rather than hardcoding it.

STEP 2 (untapered limit -- the cheap, diagnostic one). A waisted MMI whose
waist parameters are degenerate must reproduce the STRAIGHT device's answer
exactly. This is checked two ways:

  a) GEOMETRICALLY, free: body_transition_length=0 and
     body_waist_width=mm_width must each give is_body_waisted=False and a
     uniform width profile.
  b) NUMERICALLY, one short EME solve: a body_waist_width == mm_width
     "waist" must return the same S-matrix as the straight device to
     numerical precision. This catches a whole class of bug the geometric
     check cannot -- e.g. cell boundaries being placed differently, or the
     waisted code path taking a different branch that happens to agree on
     width but not on discretisation.

(b) is what makes this worth cluster time at all; (a) alone is already
covered by tests/test_geometry.py::TestMMIBodyWaist, added 2026-08-29.

WHY THIS RUNS BEFORE ANY SWEEP. The 1x4's own equivalent step is named the
"cheapest and most diagnostic check" in its plan, and this device inherits a
control that is genuinely settled (resolved/mmi_2x2_validation.md -- all four
checks passed), so any deviation here is attributable to the waist rather
than to an unvalidated baseline.

Run:
    python .../mmi_2x2_waist_geometry_check.py            # steps 1 + 2a (free)
    MMI2X2_WG_SOLVE=1 python .../mmi_2x2_waist_geometry_check.py   # + step 2b
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")

from beef import MMI, EMESolver, SOI_220, plot_eme_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "verification",
    "mmi_2x2_waist_geometry"))
os.makedirs(output_dir, exist_ok=True)


def _env_f(name, default):
    return float(os.environ.get(name, default))


def _env_i(name, default):
    return int(os.environ.get(name, default))


# --- production constants, from the RESOLVED straight 2x2 -------------------
#     (examples/mmi_2x2/python/mmi_2x2_hybrid.py:98-107)
WL = 1.55
MM_LENGTH, MM_WIDTH = 28.0, 2.8
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 0.45, 1.3, 4.0
WAVEGUIDE_OFFSET = 1.5
WAVEGUIDE_HEIGHT = 0.22          # SOI_220
STRAIGHT_LENGTH = 0.5
N_MODES = 14                     # the straight device's settled count

# --- placeholder waist profile: NOT a design, step 6 picks the real one -----
WW = _env_f("MMI2X2_WG_WW", 2.24)         # 80% of mm_width
LW = _env_f("MMI2X2_WG_LW", 8.0)
LT = _env_f("MMI2X2_WG_LT", 6.0)

RES = _env_f("MMI2X2_WG_RES", 30.0)       # cheap; step 5 sets the real one
OXW = _env_f("MMI2X2_WG_OXW", 5.0)        # straight 2x2's own oxide_width
OXH = _env_f("MMI2X2_WG_OXH", 1.5)
DO_SOLVE = _env_i("MMI2X2_WG_SOLVE", 0)

_common = dict(
    n_inputs=2, n_outputs=2, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=WAVEGUIDE_HEIGHT,
)

print("=" * 78)
print("2x2 MMI, WAIST body -- steps 1 & 2")
print(f"  straight control: L={MM_LENGTH} W={MM_WIDTH} (resolved device)")
print(f"  placeholder waist: Ww={WW} Lw={LW} Lt={LT}  (NOT a design point)")
print("=" * 78)

straight = MMI(**_common)
waisted = MMI(**_common, body_waist_width=WW, body_waist_length=LW,
              body_transition_length=LT)

# ---------------------------------------------------------------------------
# STEP 1 -- geometry sanity
# ---------------------------------------------------------------------------
print("\n--- STEP 1: geometry ---")
print(f"straight  is_body_waisted={straight.is_body_waisted}")
print(f"waisted   is_body_waisted={waisted.is_body_waisted}")

half = MM_LENGTH / 2
(din0, din1), (dout0, dout1) = waisted._body_waist_bounds()
print(f"\nwaist bounds: transitions [{din0:.2f},{din1:.2f}] and "
      f"[{dout0:.2f},{dout1:.2f}]; waist [{din1:.2f},{dout0:.2f}]")
print(f"stub length either side: {din0 - (-half):.2f} um")

print("\nwidth profile at the five segment boundaries:")
for label, x in [("in face   ", -half), ("stub->tran", din0),
                 ("tran->wais", din1), ("waist mid ", 0.0),
                 ("wais->tran", dout0), ("tran->stub", dout1),
                 ("out face  ", half)]:
    print(f"  {label} x={x:+7.2f}  width={waisted.body_width_at(x):6.3f} um")

# Port clearance, re-derived (never hardcoded) -- see the validation doc.
needed = 2 * (max(abs(y) for y in waisted.output_ys) + TAPER_WIDTH / 2)
headroom = MM_WIDTH - needed
print(f"\nport clearance: ports+tapers span {needed:.3f} um of the "
      f"{MM_WIDTH:.3f} um face")
print(f"  -> headroom {headroom:.3f} um ({100*headroom/MM_WIDTH:.1f}%)")
if headroom <= 1e-9:
    print("  -> ZERO headroom: the FACE-TAPER profile is geometrically")
    print("     IMPOSSIBLE on this device. The interior waist is the only")
    print("     body profile it can host. (This is the premise of the whole")
    print("     validation doc -- it is re-derived here, not assumed.)")

# Render both from one config.
solver_kw = dict(oxide_width=OXW, oxide_height=OXH, res=RES,
                 material_stack=SOI_220, mode_solver="tidy3d")
lengths = [STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH,
           STRAIGHT_LENGTH]
labels = ["wg_in", "taper_in", "MMI body", "taper_out", "wg_out"]

for name, geom in [("straight", straight), ("waisted", waisted)]:
    solver = EMESolver(num_modes=[1, N_MODES, N_MODES, N_MODES, 1],
                       num_cells=[1, 1, 1, 1, 1], lengths=lengths,
                       **solver_kw)
    path = os.path.join(output_dir, f"geometry_{name}.png")
    plot_eme_geometry(geom, solver, region_labels=labels, output_path=path,
                      title=(f"2x2 MMI -- {name} body"
                             + ("" if name == "straight" else
                                f" (Ww={WW}, Lw={LW}, Lt={LT})")),
                      side="out0", output_idx=0)
    print(f"rendered -> {path}")

# ---------------------------------------------------------------------------
# STEP 2a -- untapered limit, geometric (free)
# ---------------------------------------------------------------------------
print("\n--- STEP 2a: untapered limit, geometric ---")
degenerate = {
    "body_transition_length=0": MMI(**_common, body_waist_width=WW,
                                    body_waist_length=LW,
                                    body_transition_length=0.0),
    "body_waist_width=mm_width": MMI(**_common, body_waist_width=MM_WIDTH,
                                     body_waist_length=LW,
                                     body_transition_length=LT),
}
ok = True
for label, m in degenerate.items():
    uniform = all(
        abs(m.body_width_at(x) - MM_WIDTH) < 1e-12
        for x in np.linspace(-half, half, 101)
    )
    flag = (not m.is_body_waisted) and uniform
    ok = ok and flag
    print(f"  {label:28s} is_body_waisted={m.is_body_waisted!s:5s} "
          f"uniform={uniform!s:5s}  {'OK' if flag else 'FAIL'}")
print(f"  step 2a: {'PASS' if ok else 'FAIL'}")

# ---------------------------------------------------------------------------
# STEP 2b -- untapered limit, numeric (one short EME solve)
# ---------------------------------------------------------------------------
if not DO_SOLVE:
    print("\n--- STEP 2b: SKIPPED (set MMI2X2_WG_SOLVE=1) ---")
    print("Steps 1 and 2a are free; 2b costs a real EME solve and belongs in")
    print("a batch job. Run it before trusting any waisted result.")
else:
    print("\n--- STEP 2b: untapered limit, numeric ---")
    print("A body_waist_width == mm_width 'waist' must reproduce the straight")
    print("device's S-matrix to numerical precision.")
    wavelengths = np.array([WL])
    results = {}
    for name, geom in [
        ("straight", straight),
        ("degenerate-waist", degenerate["body_waist_width=mm_width"]),
    ]:
        solver = EMESolver(num_modes=[1, N_MODES, N_MODES, N_MODES, 1],
                           num_cells=[1, 1, 1, 1, 1], lengths=lengths,
                           **solver_kw)
        S = solver.solve(geom, wavelengths=wavelengths)
        # SMatrix stores dict[(port, in_mode, out_mode)] -> complex[M];
        # compare EVERY stored coefficient, not a hand-picked one.
        results[name] = {k: np.asarray(v) for k, v in S._data.items()}
        print(f"  {name}: solved, {len(results[name])} S entries")

    a, b = results["straight"], results["degenerate-waist"]
    if set(a) != set(b):
        print(f"  FAIL: key sets differ")
        print(f"    only in straight:  {sorted(set(a) - set(b))[:5]}")
        print(f"    only in waisted:   {sorted(set(b) - set(a))[:5]}")
    else:
        dmax = max(float(np.max(np.abs(a[k] - b[k]))) for k in a)
        worst = max(a, key=lambda k: float(np.max(np.abs(a[k] - b[k]))))
        print(f"  compared {len(a)} coefficients")
        print(f"  max |S_straight - S_degenerate| = {dmax:.3e}  (worst key {worst})")
        print(f"  step 2b: {'PASS' if dmax < 1e-9 else 'FAIL'} (bar: 1e-9)")

print("\n" + "=" * 78)
print("Steps 1/2 are gates, not results. A design point is step 6's job, and")
print("step 4 fixes the pass/fail metric BEFORE any sweep runs -- see the doc.")
