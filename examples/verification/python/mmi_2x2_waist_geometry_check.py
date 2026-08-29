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


def waist_region_cell_lengths(mmi, c_stub, c_trans, c_waist):
    """The body's own 5-segment cell-length list, aligned exactly to
    stub/transition/waist/transition/stub.

    STEP 3. This is what `EMESolver(region_cell_lengths=...)` was added for
    (2026-08-29): the default "divide each region evenly by its num_cells"
    cannot place a boundary at din0/din1/dout0/dout1, so with an even division
    a single cell straddles a width discontinuity and the body's profile is
    misrepresented no matter how many cells are used.

    The 1x4's own sweep script solved this by BYPASSING
    EMESolver.solve()/_build_cells() entirely and hand-calling
    mw.create_cells(). That workaround is now retired: the library places the
    boundaries, so this device keeps caching, calibration and every other
    solve() feature instead of reimplementing them.
    """
    half = mmi.mm_length / 2
    if not mmi.is_body_waisted:
        return None
    (din0, din1), (dout0, dout1) = mmi._body_waist_bounds()
    segs = [(-half, din0, c_stub), (din0, din1, c_trans),
            (din1, dout0, c_waist), (dout0, dout1, c_trans),
            (dout1, half, c_stub)]
    lens = []
    for x0, x1, n in segs:
        if n <= 0 or x1 - x0 <= 1e-12:
            continue
        seg = (x1 - x0) / n
        lens.extend([seg] * n)
    return lens


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
    # NEAR-degenerate, not exactly degenerate. body_waist_width == mm_width
    # sets is_body_waisted=False, so that geometry takes the STRAIGHT code
    # path and trivially returns a bit-identical answer (measured: exactly
    # 0.000e+00, job 12482595) -- which proves the fallback works but never
    # exercises the waisted cell-construction path this gate exists to test.
    # A waist that is narrower by EPS instead is genuinely waisted
    # (is_body_waisted=True, real 5-segment cell construction) and must
    # CONVERGE to the straight answer as EPS -> 0.
    # A single EPS with a fixed bar cannot tell "physics" from "bug": an MMI's
    # self-imaging length goes as W^2, so a fractional width change of
    # 2*dW/W already shifts the beat phase by ~2*pi*2*dW/W over one beat
    # length -- at EPS=1e-3 that is ~4.5e-3 rad, the same order as a
    # measured |dS| of 1.1e-2 (job 12482618, which FAILED a 5e-3 bar).
    # So SWEEP EPS instead and test the SCALING:
    #   |dS| proportional to EPS  -> the waisted path converges correctly,
    #                                the discretisation is sound, and any
    #                                fixed bar was simply mis-set.
    #   |dS| flattening to a floor -> a real bug in the waisted cell
    #                                construction that does NOT vanish in
    #                                the limit.
    EPS_LIST = [float(v) for v in
                os.environ.get("MMI2X2_WG_EPS", "1e-2,1e-3,1e-4,1e-5").split(",")]

    wavelengths = np.array([WL])
    cases = [
        ("straight", straight),
        ("degenerate-waist", degenerate["body_waist_width=mm_width"]),
    ]
    for eps in EPS_LIST:
        near = MMI(**_common, body_waist_width=MM_WIDTH - eps,
                   body_waist_length=LW, body_transition_length=LT)
        assert near.is_body_waisted, f"EPS={eps} must take the waisted path"
        cases.append((f"near-eps{eps:g}", near))

    results = {}
    for name, geom in cases:
        solver = EMESolver(num_modes=[1, N_MODES, N_MODES, N_MODES, 1],
                           num_cells=[1, 1, 1, 1, 1], lengths=lengths,
                           **solver_kw)
        S = solver.solve(geom, wavelengths=wavelengths)
        # SMatrix stores dict[(port, in_mode, out_mode)] -> complex[M];
        # compare EVERY stored coefficient, not a hand-picked one.
        results[name] = {k: np.asarray(v) for k, v in S._data.items()}
        print(f"  {name}: solved, {len(results[name])} S entries")

    def _dmax(name):
        b = results[name]
        if set(a) != set(b):
            return None
        return max(float(np.max(np.abs(a[k] - b[k]))) for k in a)

    a = results["straight"]

    # -- exact degenerate: fallback check only --
    d_exact = _dmax("degenerate-waist")
    exact_ok = d_exact is not None and d_exact < 1e-9
    print(f"\n  exact degenerate (is_body_waisted=False -> STRAIGHT path):")
    print(f"    max |dS| = {d_exact:.3e}  -> {'PASS' if exact_ok else 'FAIL'}")
    print("    (proves only that the fallback works, not the waisted path)")

    # -- the real test: does |dS| scale with EPS? --
    print(f"\n  waisted-path convergence (the load-bearing check):")
    print(f"    {'EPS':>10} {'max|dS|':>12} {'ratio':>8}  {'|dS|/EPS':>10}")
    prev = None
    rows = []
    for eps in EPS_LIST:
        d = _dmax(f"near-eps{eps:g}")
        if d is None:
            print(f"    {eps:10.1e}  KEY-SET MISMATCH")
            continue
        ratio = (prev / d) if (prev is not None and d > 0) else float("nan")
        rows.append((eps, d))
        print(f"    {eps:10.1e} {d:12.3e} {ratio:8.2f}  {d/eps:10.3e}")
        prev = d

    # Linear convergence means |dS|/EPS is CONSTANT across the sweep.
    verdict = "INCONCLUSIVE"
    if len(rows) >= 3:
        slopes = [d / eps for eps, d in rows]
        spread = max(slopes) / min(slopes) if min(slopes) > 0 else float("inf")
        smallest_eps, smallest_d = rows[-1]
        if spread < 3.0:
            verdict = ("PASS -- |dS| scales linearly with EPS "
                       f"(|dS|/EPS constant to {spread:.2f}x). The waisted "
                       "cell construction CONVERGES to the straight body; "
                       "the earlier fixed-bar FAIL was a mis-set bar, not a bug.")
        elif smallest_d > 0.1 * rows[0][1]:
            verdict = ("FAIL -- |dS| FLOORS instead of scaling: shrinking EPS "
                       f"{rows[0][0]:g} -> {smallest_eps:g} ({rows[0][0]/smallest_eps:.0f}x) "
                       f"only moved |dS| {rows[0][1]:.2e} -> {smallest_d:.2e}. "
                       "A discretisation error in the waisted path that does "
                       "NOT vanish in the limit. DO NOT trust waisted results "
                       "-- and re-examine the 1x4's, which share this path.")
        else:
            verdict = (f"UNCLEAR -- |dS|/EPS varies {spread:.1f}x across the "
                       "sweep. Neither clean linear convergence nor a flat "
                       "floor; extend the sweep before drawing a conclusion.")
    print(f"\n  step 2b verdict: {verdict}")

# ---------------------------------------------------------------------------
# STEP 3 -- segment-aligned cell construction via region_cell_lengths
# ---------------------------------------------------------------------------
print("\n--- STEP 3: segment-aligned cell construction ---")
C_STUB = _env_i("MMI2X2_WG_CSTUB", 1)
C_TRANS = _env_i("MMI2X2_WG_CTRANS", 4)
C_WAIST = _env_i("MMI2X2_WG_CWAIST", 2)
body_cells = waist_region_cell_lengths(waisted, C_STUB, C_TRANS, C_WAIST)
n_body = len(body_cells)
print(f"body subdivided into {n_body} cells "
      f"(stub={C_STUB}, trans={C_TRANS}, waist={C_WAIST}, per side)")
print(f"  sums to {sum(body_cells):.6f} um vs mm_length {MM_LENGTH}")

aligned = EMESolver(
    num_modes=[1, N_MODES, N_MODES, N_MODES, 1],
    num_cells=[1, 1, n_body, 1, 1], lengths=lengths,
    region_cell_lengths={2: body_cells},   # region 2 is the MMI body
    **solver_kw)
flat = aligned._flat_cell_lengths()
print(f"  solver accepted it: {len(flat)} cells total, "
      f"sum {sum(flat):.6f} um vs device {sum(lengths):.6f} um")

# Every cell must sit on ONE cross-section: check no cell straddles a
# width discontinuity, by comparing the width at each cell's two edges.
x = -MM_LENGTH / 2
straddlers = []
for L in body_cells:
    w0, w1 = waisted.body_width_at(x), waisted.body_width_at(x + L)
    mid = waisted.body_width_at(x + L / 2)
    # inside a transition the width varies linearly, so edges differ; that is
    # expected and fine. What must NOT happen is a cell spanning a segment
    # BOUNDARY, i.e. part-uniform part-ramp. Detect via second difference.
    if abs((w0 + w1) / 2 - mid) > 1e-9:
        straddlers.append((x, x + L, w0, mid, w1))
    x += L
print(f"  cells straddling a segment boundary: {len(straddlers)} "
      f"({'OK' if not straddlers else 'PROBLEM'})")

# Contrast: what an EVEN division would have done with the same cell count.
even = [MM_LENGTH / n_body] * n_body
x, even_straddlers = -MM_LENGTH / 2, 0
for L in even:
    w0, w1 = waisted.body_width_at(x), waisted.body_width_at(x + L)
    mid = waisted.body_width_at(x + L / 2)
    if abs((w0 + w1) / 2 - mid) > 1e-9:
        even_straddlers += 1
    x += L
print(f"  (an EVEN division of the same {n_body} cells would straddle "
      f"{even_straddlers} -- this is what region_cell_lengths fixes)")
print(f"  step 3: {'PASS' if not straddlers else 'FAIL'}")

print("\n" + "=" * 78)
print("Steps 1/2/3 are gates, not results. A design point is step 6's job, and")
print("step 4 fixes the pass/fail metric BEFORE any sweep runs -- see the doc.")
