# 2×2 MMI with a WAIST body — build & validation plan

Opened **2026-08-29**, at the user's request: apply the interior-waist body
profile (`body_waist_width`/`body_waist_length`/`body_transition_length`,
added 2026-08-27/28 for the 1×4) to the **2×2 MMI on SOI_220**.

> **This doc is a PLAN, not a gap list.** Nothing below has been simulated for
> a waisted 2×2. Every number quoted for the *straight* 2×2 is **MEASURED** and
> cited to `resolved/mmi_2x2_validation.md`; everything about what a *waisted*
> 2×2 does is **ARGUED** and is a prediction to be checked. It converts to a
> normal gap list against the [four checks](README.md#the-four-checks) once a
> design point exists and the scripts run.
>
> **Read `mmi_1x4_tapered_body_validation.md` first.** The same profile was
> tried on the 1×4 and **failed outright** — twice closed, on three Ww values,
> doubled taper resolution, and a length search past 2× the device's own
> length. This plan exists because the 2×2 differs in the specific ways that
> would discriminate *why* it failed, not because the 1×4 result is being
> ignored.

---

## Why this device, and why it is a real test rather than a retry

The 1×4's waist failure is **not** fully root-caused. The best evidence
(`mmi_1x4_tapered_body_validation.md`, job 12449008) is that at two lengths the
device reached **better balance than its straight control** (1.22e-03/1.60e-03
vs 2.38e-03) while `totalT` stayed catastrophic (0.47/0.18). Good phase with
bad power argues for **real radiative loss at the width transitions** rather
than self-imaging beat-length desynchronization. That hypothesis makes a
falsifiable prediction, and the 2×2 is where it can be tested:

| | 1×4 (failed) | 2×2 (this plan) |
|---|---|---|
| platform | SIN_300, Δn² = **1.90** | SOI_220, Δn² = **10.01** |
| `mm_width` | 12.0 µm | **2.8 µm** |
| `mm_length` | 39.90 µm | **28.0 µm** |
| body aspect W/L | 0.301 | **0.100** |
| absolute narrowing at 25% | 3.0 µm ≈ **1.9 λ** | 0.7 µm ≈ **0.45 λ** |
| observable | 4-way power split + balance | **bar/cross transfer matrix — phase, not just power** |

**The discriminating argument (ARGUED, this is the hypothesis under test).**
Radiative loss at a taper scales with how non-adiabatic the width change is
and with how weakly the mode is confined. The 2×2 is favourable on *both*
axes at once: its index contrast is **5.3× stronger** in Δn² (so a given
sidewall angle radiates far less), and the absolute width excursion needed for
the same *fractional* narrowing is **4.3× smaller** in µm. For a fixed
transition length the local sidewall angle goes as Δ(width)/2·Lt, so matching
the 1×4's own (failed) 25%-narrowing geometry at Lt = 8 µm would put the 2×2 at
roughly **1/4 the sidewall angle in a 5× higher-contrast waveguide** — i.e.
several times more adiabatic on both counts simultaneously.

**So the two outcomes are both informative, which is the point of running it:**
- **If the waisted 2×2 works** → the 1×4's failure was transition radiation,
  the profile is sound, and it is a low-contrast/wide-body device that cannot
  host it. The mechanism question closes.
- **If the waisted 2×2 fails the same way** (totalT collapse at every length,
  balance sometimes fine) → transition radiation is **refuted** as the
  explanation, because this device removes both of the conditions that were
  supposed to be causing it. The failure would then be intrinsic to perturbing
  a self-imaging body's width at all, which is a much stronger and more
  general statement than the 1×4 alone can support.

⚠ **Do not let step 6 quietly become an unbounded parameter hunt.** The 1×4
consumed three Ww values, a resolution doubling and a 9-point length search
before re-closing. Step 6 below carries a **pre-registered stop rule** for
exactly this reason.

---

## What is already true — MEASURED, and what it buys

### The straight 2×2 is the strongest control in the repo

`resolved/mmi_2x2_validation.md` — **all four checks passed**, device resolved.
Directly usable here:

- **Geometry** (`examples/mmi_2x2/python/mmi_2x2_hybrid.py:98-107`):
  `n_inputs=2, n_outputs=2, mm_length=28.0, mm_width=2.8,
  waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
  waveguide_offset=1.5, waveguide_height=0.22`, `SOI_220`, band 1.530–1.565.
- **Check 2 (mode identity): passed, "cleanest device audited"** — and
  reconfirmed after both duplicate-mode guards existed (job 11588569).
- **Check 3 (resolution): converged at res 100 on Δβ_01** (job 11590785).
- **Check 4 (whole-device FDTD): measured**, 0.82%/0.62% max |ΔT|, with the
  gap **sign-flipping across the band** (job 11585844/11585849).
- **`N_MODES = 14`**, deliberately not 15 — band 14 sits below the cladding
  index (a radiation mode, dropped at no cost).
- **A β-calibration already exists and is cheap here**: 30.6 h → **4.14 h
  measured** (job 11937894), and ~3.4 h with the seed cached.

**This means the waisted 2×2 inherits a control that is genuinely settled**,
which the 1×4's tapered variant never had to the same standard. Any deviation
is attributable to the waist rather than to an unvalidated baseline.

### The geometry layer already supports this — MEASURED (read from source)

`beef/geometry/mmi.py` — the waist parameters are **port-count agnostic**. They
touch only the body: `_body_waist_bounds()` (line 351) derives the transition
bands from `body_transition_length`/`body_waist_length` about x=0, and
`body_width_at()` (line 361) is a pure function of x. Nothing in that path
reads `n_inputs`/`n_outputs`. The constructor already enforces the two guards
that matter (`mmi.py:128-145`): waist and face-taper parameters are **mutually
exclusive**, and `2·body_transition_length + body_waist_length ≤ mm_length`.

**No geometry change is needed for a 2×2 waist** — and this is now **MEASURED,
not just read from source**: a waisted 2×2 constructs correctly with the
expected bounds and widths (see step 1's measured block). A render and the
boundary assertions are still owed.

### ⚠ Two things that do NOT carry over from the 1×4 — read before reusing code

1. **`mmi_1x4_waist_profile_sweep.py` is not reusable as-is.** It hardcodes
   `n_inputs=1, n_outputs=4`, the 1×4's own constants, and — critically — a
   `solve_reuse()` whose **width-keyed cache is only valid on-axis**. Extending
   that cache across the on-axis/off-axis boundary was measured **breaking
   physics** on 2026-08-28 (job 12453222: totalT 0.9703 → 0.1507) and reverted;
   see that doc's "Efficiency attempt REVERTED" section and
   `width-only-mode-reuse-breaks-off-axis-cells` in memory. **The 2×2 has
   *both* ports off-axis** (`waveguide_offset=1.5` → inputs at y = ±0.75), so
   there is no on-axis chain at all here — the 1×4's body-only cache rule must
   be re-derived for this device, not copied.
2. **The 2×2 is a transfer matrix, not a power splitter.** Its figures of merit
   are **T_bar/T_cross and the phase between them**, not a 4-way balance. A
   waist that preserved total power but rotated the bar/cross phase would be a
   *failure* here and would not show up in the 1×4's own metric at all. Step 4
   defines the metric before any sweep runs, deliberately.

---

## The ordered steps

Each step names its **deliverable** and its **evidence** (what makes it done).
Steps 1–3 are cheap and must precede any expensive sweep.

### Step 1 — Geometry sanity for a waisted 2×2 (minutes, no cluster)

**Deliverable.** A rendered waisted-2×2 layout at the straight device's own
constants plus a placeholder waist, confirming: both faces sit at `mm_width`,
the interior narrows symmetrically, the two input and two output tapers still
land on `input_ys`/`output_ys` without overlapping the narrowed body, and
`is_body_waisted` is True.

**Evidence.** A geometry plot plus assertions on `body_width_at()` at the five
segment boundaries. Re-derive the port-clearance constraint from `output_ys`
rather than hardcoding it, the way the 1×4's own
`test_the_real_1x4_asymmetric_profile` does.

#### ✅ Partially MEASURED already, 2026-08-29 (login node, seconds)

Two of this step's claims were checked on the real `MMI` class rather than left
as arithmetic:

- **The 2×2's output face is exactly 100% consumed by its ports.** Ports at
  `output_ys = [-0.75, +0.75]`, each taper `taper_width=1.3` wide, so the pair
  spans **2.800 µm against `mm_width` = 2.800 µm — headroom 0.000 µm (0.0%)**.
  This is *tighter* than the 1×4's own blocker (which still had 4.2% at its
  output face) and it means **the face-taper profile is not merely unattractive
  here, it is geometrically impossible.** The interior waist is the only body
  profile this device can host at these constants — a stronger motivation than
  the 1×4 ever had.
- **The waist constructs correctly at 2×2 parameters.** `MMI(..., n_inputs=2,
  n_outputs=2, body_waist_width=2.24, body_waist_length=8.0,
  body_transition_length=6.0)` gives `is_body_waisted=True`, bounds
  `((-10.0, -4.0), (4.0, 10.0))`, and widths **2.800 / 2.240 / 2.800** µm at
  in-face / mid / out-face — confirming the port-count-agnostic claim above on
  the real object, not just from reading the source.

#### ✅ STEP 1 DONE 2026-08-29 — render + tests landed

`examples/verification/python/mmi_2x2_waist_geometry_check.py` renders the
straight control and the waisted device **from one config** (so the only
difference is the body profile) and prints the width profile at every segment
boundary. Output under
`example_output/verification/mmi_2x2_waist_geometry/`.

Measured at the placeholder profile (Ww=2.24, Lw=8.0, Lt=6.0 — **not** a
design point, step 6's job):

```
waist bounds: transitions [-10.00,-4.00] and [4.00,10.00]; waist [-4.00,4.00]
stub length either side: 4.00 um
  in face    x= -14.00  width= 2.800 um
  stub->tran x= -10.00  width= 2.800 um
  tran->wais x=  -4.00  width= 2.240 um
  waist mid  x=  +0.00  width= 2.240 um
  tran->stub x= +10.00  width= 2.800 um
  out face   x= +14.00  width= 2.800 um
port clearance: 2.800 um of the 2.800 um face -> headroom 0.000 um (0.0%)
```

Render visually confirms the intended symmetric hourglass with both faces at
`mm_width` and the port tapers landing exactly on them.

#### ⚠ FINDING: the waist profile had ZERO test coverage until today

Writing this step's assertions turned up something the 1×4 work did not:
**`grep -r waist tests/` returned nothing.** The profile shipped 2026-08-27/28
and has since produced several jobs' worth of real physics results (every 1×4
waist sweep), with no regression test anywhere. The claim recorded at the time
— *"780 passed, 0 failed"* — only established the feature is **inert to the
existing suite**, which it is by construction: every waist parameter defaults
to off. It never exercised the profile itself. The "verified standalone"
checking was real but ad-hoc and was never committed.

**Closed:** `tests/test_geometry.py::TestMMIBodyWaist`, **17 tests, all
passing** (123 pass across the geometry suite, no regressions). Covers the
untapered limit (three separate degenerate paths), bounds symmetry, the
`mm_length`-independence the sweep scripts rely on, five-boundary widths,
ramp linearity and mirror symmetry, all four constructor guards, the cut
validation contract (transition bands refuse; uniform waist and stubs allow),
and **both real devices** — the 1×4's production waist (retroactively) and
this doc's 2×2, with the zero-headroom port constraint re-derived from
`output_ys` rather than hardcoded, so a future port-layout change fails the
test instead of silently invalidating this doc's premise.

**Cost.** Minutes.

### Step 2 — Untapered limit: the waist must reproduce the straight 2×2 (cheap)

**Deliverable.** With `body_transition_length → 0` (or
`body_waist_width = mm_width`), the waisted code path must return the straight
device's answer to numerical precision.

**Evidence.** `is_body_waisted` False in both degenerate cases (already
asserted by the constructor's own logic, `mmi.py:337-349` — confirm it holds
at 2×2 parameters), and a single-wavelength EME solve matching the straight
2×2's cached result. **This is the cheapest and most diagnostic check in the
plan** — the 1×4 doc says the same thing about its own step 5 — and it must
pass before any sweep is submitted.

Split into two halves, because they cost very different amounts:

- **2a, geometric — free.** Both degenerate paths must give
  `is_body_waisted=False` and a uniform width profile.
- **2b, numeric — one short EME solve.** A `body_waist_width == mm_width`
  "waist" must return the same S-matrix as the straight device. This catches
  a class of bug 2a cannot: cell boundaries placed differently, or the
  waisted code path branching somewhere that agrees on *width* but not on
  *discretisation*. Compares **every** stored coefficient in the `SMatrix`,
  not a hand-picked one, against a 1e-9 bar.

#### ✅ STEP 2a DONE 2026-08-29 — PASS

```
  body_transition_length=0     is_body_waisted=False uniform=True   OK
  body_waist_width=mm_width    is_body_waisted=False uniform=True   OK
  step 2a: PASS
```

(uniformity checked at 101 points across the full body, not just the
boundaries.)

#### ⏳ STEP 2b IN FLIGHT — job 12482595

`batch/verification/run_mmi_2x2_waist_geometry_check.batch`
(`MMI2X2_WG_SOLVE=1`). **No waisted 2×2 result should be trusted until this
passes** — that is the entire point of the gate.

**Cost.** One short job.

### Step 3 — Cell-construction path for a 5-segment 2×2 body

**Deliverable.** A 2×2 solve that places EME cell boundaries exactly at the
waist's internal edges (stub/transition/waist/transition/stub).

**Evidence.** Cell count and per-cell widths printed and checked against
`body_width_at()` at each cell centre. **Note the known blocker, carried over
and still true:** `EMESolver._build_cells()` divides each region *uniformly* by
its cell count and cannot place boundaries at the waist edges, which is why the
1×4 script hand-builds the body's cell-length list and calls
`mw.create_cells()` directly. **Decide deliberately** whether to (a) repeat
that bypass for the 2×2, or (b) fix `_build_cells()` to accept explicit
segment boundaries — (b) is the better library outcome and would retire a
documented workaround, but it is a **shared-code change** and per this repo's
standard an "inert" refactor is an argument until a run says otherwise (the
directional coupler is the cheapest confirming re-run; see
`mmi_1x4_sin300_validation.md` step 2's Evidence #3 for the pattern).

⚠ **Do not carry over the 1×4's width-keyed reuse cache without re-deriving
it.** See the warning above — the 2×2 has no on-axis chain.

**Cost.** Hours of implementation; one short confirming job.

### Step 4 — Define the figure of merit BEFORE sweeping (no compute)

**Deliverable.** A written, pre-registered pass/fail criterion for a waisted
2×2, in terms of **T_bar, T_cross, excess loss and the bar/cross phase**,
scored against the straight control's own measured values across the
1.530–1.565 band.

**Evidence.** The criterion written into the sweep script itself, so results
are scored rather than eyeballed. **This step exists because of a specific
1×4 near-miss:** its sweep ranked candidates by `balance` alone, which is why
two points with *better-than-control balance* and catastrophic `totalT`
(0.47/0.18) surfaced as "BEST" in the summary table. A 2×2 has more ways to be
subtly wrong (phase), so the metric must be fixed in advance.

**Cost.** None.

### Step 5 — Convergence for the waisted body (cheap-to-moderate)

Two axes, neither inherited:

- **Transition cell count.** Sweep it and show the S-matrix plateaus. **The
  1×4 measured this mattering**: doubling `C_TRANS` 6→12 visibly changed the
  field character (job 12446756), though it did not rescue the device. Do not
  pick a count for legibility and quote results from it.
- **`num_modes` through the narrowed section.** The body *narrows*, so the
  guided-mode count **falls** in the waist — the opposite of the 1×4-tapered
  doc's widening case. The binding constraint is the **widest** cell (the
  faces, at `mm_width`), where the straight device's `N_MODES=14` already
  applies. ⚠ Heed the 4×4's warning quoted in the 1×4 doc: a probe that hits
  its own request cap only ever measured a floor, never found a ceiling —
  raise the cap until a real ceiling appears.

**Evidence.** Plateau tables for both axes, with the tolerance stated.

**Cost.** Moderate; cheaper than the 1×4's equivalent (smaller domain, fewer
body modes).

### Step 6 — Find a design point, with a PRE-REGISTERED STOP RULE

**Deliverable.** Either a (Ww, Lw, Lt, mm_length) point that beats the straight
2×2 on step 4's criterion, or a decisive negative.

**Search plan.** Sweep `mm_length` per profile, as the 1×4 did. Start
**gentle** — the 1×4's own hard-won lesson is that aggressive narrowing
(50%) fails catastrophically and is not informative. Suggested first grid:
`Ww/mm_width ∈ {0.90, 0.80}` (i.e. Ww ≈ 2.52, 2.24 µm), `Lt` chosen to match
or better the 1×4's sidewall angle, `Lw ≈ 6-10 µm`, ~5 lengths around 28.0.

**⚠ STOP RULE, pre-registered here so it cannot be moved later.** Stop and
write the negative result if **either**:
1. Two Ww values × their length windows produce **no point** within 5% of the
   straight control's excess loss, **or**
2. The best point is not **both** shorter than 28.0 µm **and** within
   step 4's criterion.

The 1×4 spent three Ww values, a resolution doubling and a 9-point length
extension past the stop-worthy point. **A clean negative here is a publishable
mechanism result** (it refutes the transition-radiation hypothesis, see the top
of this doc) — it is not a failure to be swept around.

**Cost.** ~1000-1300 s per (profile, length) point was the 1×4's rate at
design settings; expect **less** here (2.8 µm body vs 12.0 µm, smaller domain).
Budget a first grid at a few hours, and **run one smoke point first** — the
same two-step discipline every expensive job in this repo follows.

### Step 7 — Field plots at the decisive point (cheap, do this either way)

**Deliverable.** Straight-vs-waisted 2×2 EME field plots at whatever step 6
lands on, from **one config** so the only difference is the waist.

**Evidence.** The rendered comparison. `plot_mmi1x4_waist_profile_eme_fields.py`
is the working template (note its own docstring: it *copies* the sweep script's
body-segment logic rather than importing it, because that script has no
`__main__` guard and would run its whole grid as an import side effect).

**Why it is worth doing on a negative too:** on the 1×4 this is exactly what
surfaced the spurious mid-body refocus and motivated the mechanism question
now being tested.

### Step 8 — The four checks (only if step 6 succeeds)

Standard, per [README](README.md#the-four-checks) — and **not** to be started
before a design point exists. Two device-specific notes:

- **Check 1** is more interesting than usual: there is no single body
  cross-section, so state which one the check is asked about — probably the
  **uniform waist**, with the transitions covered by step 5's convergence
  instead. Same framing the 1×4-tapered doc uses.
- **`calibrate_body_neff` is not free here and its meaning changes.** The
  straight 2×2's calibration is 3.4-4.1 h *for a uniform body*; a waisted body
  has no single β to calibrate. Decide explicitly what is being calibrated
  before spending it — the 1×4-tapered doc records `calibrate_body_neff`
  being **refused** for exactly this reason.

### Step 9 — Wiring, only once the physics is settled

`examples/mmi_2x2_waist/python/mmi_2x2_waist_eme.py` and `_hybrid.py`, plus a
**`Makefile` target** — note that `mmi_4x4_90hybrid` has none, which is how
that device stayed unrun.

---

## Cross-references

- `mmi_1x4_tapered_body_validation.md` — the same profile's failure on the
  1×4: face-taper dead end, waist closed → reopened → re-closed, and the
  good-balance/bad-totalT finding this plan is built to test.
- `resolved/mmi_2x2_validation.md` — the straight control: all four checks,
  `N_MODES=14`, res-100 convergence, calibration cost model.
- `README.md` — the four checks.
- Memory: `width-only-mode-reuse-breaks-off-axis-cells` (why the 1×4's cache
  cannot be copied here), `mmi-1x4-tapered-body-both-profiles-closed`.
