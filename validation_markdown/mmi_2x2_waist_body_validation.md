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

#### ✅ STEP 2b DONE 2026-08-29 — PASS, but only after the test itself was fixed twice

**This step is worth reading as a method lesson, not just a result.** The
first two versions of it were both too weak to gate anything, and the second
one produced a false alarm.

**Attempt 1 (job 12482595) — vacuous.** Compared the straight device against
`body_waist_width == mm_width` and got `max |ΔS| = 0.000e+00`, exactly zero,
in 56 s. That looked like a perfect pass and was nearly meaningless: an exact
degenerate waist sets `is_body_waisted=False`, so it takes the **identical
straight-body code path**. Bit-identical zero is what taking the same branch
produces. It confirms the *fallback*, and never touches the waisted
cell-construction path the gate exists to test.

**Attempt 2 (job 12482618) — FAILED, and the failure was spurious.** Switched
to a genuinely waisted **near**-degenerate case (`Ww = mm_width − 1e-3`, which
does set `is_body_waisted=True`). Measured `|ΔS| = 1.127e-02` against a 5e-3
bar → FAIL. **The bar was the problem, not the code.** An MMI's self-imaging
length goes as W², so a fractional width change shifts the beat phase by
~2π·2δW/W over a beat length — at δW=1e-3 that is ~4.5e-3 rad, the same order
as the measured deviation. A single EPS with a fixed threshold **cannot
distinguish real physics from a discretisation bug.**

**Attempt 3 (job 12482636) — the right test: SCALING, not a threshold.**

| EPS (µm) | max abs ΔS | ratio | abs ΔS / EPS |
|---|---|---|---|
| 1.0e-02 | 9.338e-02 | — | 9.34 |
| 1.0e-03 | 1.127e-02 | 8.29 | 11.27 |
| 1.0e-04 | 1.148e-03 | 9.82 | 11.48 |
| 1.0e-05 | 1.150e-04 | 9.98 | 11.50 |

**Clean first-order convergence.** Each decade of EPS reduces |ΔS| by a
factor approaching exactly 10 (8.29 → 9.82 → **9.98**), and `|ΔS|/EPS` is
constant to **1.23×** across the sweep. The waisted 5-segment cell
construction genuinely converges to the straight body in the limit.

**Verdict: step 2 PASSES.** The waisted code path is sound.

**⚠ This also clears the 1×4.** That device's waist results share this exact
code path. Had |ΔS| floored instead of scaling, every 1×4 waist conclusion —
including the twice-closed "the waist profile fails" verdict in
`mmi_1x4_tapered_body_validation.md` — would have been resting on a
contaminated path. It scales, so those conclusions stand.

#### 📌 Number to carry into steps 4 and 6: dS/dW ≈ 11.5 µm⁻¹

The converged slope means **this device's S-matrix moves ~11.5 per µm of
waist width**. Two consequences worth stating before any sweep:

1. **The 2×2 is very width-sensitive**, and a candidate design's narrowing is
   ~0.3-0.6 µm — orders of magnitude outside this linear regime, so do not
   extrapolate with this slope. It bounds *numerical* tolerance, not physics.
2. **It sets the meaningful precision floor for step 4's metric.** Distinguishing
   two candidate profiles whose widths differ by <1e-3 µm is not physically
   meaningful at this resolution — the sweep grid should not be finer than that.

**Cost.** Three jobs, all under 1 minute each.

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

#### ✅ STEP 3 DONE 2026-08-29 — option (b) taken, the library was fixed

**Decision: (b).** `EMESolver(region_cell_lengths={region_idx: [lengths]})`,
an explicit per-region subdivision. `_build_cells()` and
`_flat_cell_lengths()` are now unified onto one code path instead of
duplicating the division logic, and **the 1×4's documented bypass is
retired** — a waisted body can go through `solve()` and keep caching,
calibration and everything else, rather than reimplementing them.

**Why (b) rather than repeating the bypass:** the bypass was already costing
real capability (the 1×4's sweep script gave up `solve()` entirely) and had
already produced one genuine physics bug when its hand-rolled mode cache was
extended (job 12453222, see the warning above). Fixing the library removes
the reason to hand-roll at all.

**The problem, quantified rather than asserted** — measured on this device's
placeholder profile: with the same 12 body cells, an **even division
straddles 4 segment boundaries** (a single cell spanning part-uniform,
part-ramp, so it sits on no single cross-section); the segment-aligned
subdivision straddles **0**. Even division does not converge out of this —
more cells still misplace boundaries.

**Guarded hard, because both misuse modes are silent:** a wrong cell *count*
desyncs `_modes_per_cell()` from the cells; a wrong *sum* changes the
device's physical length while every other number still looks plausible.
Both are refused at construction. Solving an overridden region at a
different length is also refused rather than silently rescaled — those
boundaries are tied to specific internal features and have no correct
automatic rescaling.

**Cache handling (the part that had to be right):** listed in
`_METHOD_ATTRS` for `solve`/`solve_length_sweep`, so a non-default
subdivision can never be served by an evenly-divided entry; **and** in
`_DEFAULT_OMITTED_ATTRS`, so at its `None` default the key does not appear
in the fingerprint at all and **every existing EME cache entry and golden
stays valid**. Same established pattern as `minimum_run_time`/
`decay_method`/`energy_decay_dt`. Pinned by
`test_default_is_omitted_so_existing_entries_stay_valid`.

**✅ MEASURED INERT — job 12484205.** Per this repo's standard, an "inert"
refactor is an argument until a run says otherwise. Step 2b's convergence
sweep was re-run through the rewired `_build_cells()` and reproduces job
12482636 **to every printed digit** on all four EPS rows
(9.338e-02 / 1.127e-02 / 1.148e-03 / 1.150e-04, `|ΔS|/EPS` constant to
1.23×). This is a real EME solve through exactly the changed code path, run
before and after, with identical output.

**Tests.** 13 new in `tests/test_cache_config.py`
(`TestRegionCellLengthsFingerprint`, `TestRegionCellLengthsValidation`).
Full suite **810 passed, 21 skipped**, up from 780 — +17 waist geometry
(step 1), +13 here, no existing test's outcome changed.

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

#### ✅ STEP 4 DONE 2026-08-29 — baseline measured, criterion written to disk

`examples/verification/python/mmi_2x2_waist_metric.py`, job **12487988**
(res=40, 11 wavelengths, the resolved device's own settings).

**Straight control baseline:**

| | value |
|---|---|
| T_bar | 0.4284 → 0.4454 across band |
| T_cross | 0.4876 → 0.5092 |
| total | 0.9160 → 0.9546 (rising monotonically) |
| max excess loss | **0.3809 dB** |
| max imbalance | **0.0686** |
| bar/cross phase | −89.86° → −90.00°, **max deviation from quadrature 0.41°** |
| Δβ₀₁ (band centre) | 0.161948 rad/µm |
| **L_π** | **19.3988 µm** |
| mm_length / L_π | **1.4434** |

**The phase is an independent validation, not just a metric component.**
Quadrature (±90°) is the theoretical signature of a 2×2 MMI, and this device
sits within **0.41°** of it across the whole band without anything having
been tuned to produce that. It is strong evidence the port conventions,
excitation and solver config are all right — which matters because two of
them were wrong an hour earlier (see the traps below).

**The pre-registered criterion**, written to
`example_output/verification/mmi_2x2_waist_metric/criterion.json`. A
candidate passes only if **all** hold:

- **shorter** — total device length < control's 37.0 µm
- **excess loss** ≤ 0.3809 + 0.10 = **0.4809 dB**
- **imbalance** ≤ 2× control = **0.1371**
- **phase** within **2.0°** of the control's

**Ranking: only among candidates that pass all four; shortest wins.** Never
rank on one component — that is precisely the 1×4's mistake.

⚠ **Step 6 must LOAD this JSON and score against it.** Moving a threshold
after seeing results is the thing pre-registration exists to prevent.

#### ⚠ Two config traps, caught only because a NUMBER looked wrong

The first baseline attempt (job 12487486) reported **~2 dB excess loss and
0.18 imbalance on the straight, resolved, well-behaved control**. The script
ran cleanly and produced a plausible-looking table; only the *implausibility
of the physics* flagged it. Both causes were mine:

1. **Port tapers must be subdivided.** The resolved example uses
   `num_cells=[1, 15, 1, 15, 1]` — **15 cells per 4 µm port taper**. I used
   1, making each taper a single-step staircase worth ~2 dB of purely
   numerical loss. Note the shape: the **body** is the region that gets 1
   cell (it is uniform); the **tapers** need many. For a *waisted* body that
   flips — the body then needs step 3's segment-aligned subdivision.
2. **`N_MODES = 14` is not the EME mode count.** 14 is
   `mmi_2x2_hybrid.py`'s **interface** count for the hybrid cascade. The EME
   example's per-region `num_modes` is `[1, 7, 15, 7, 1]`. Same device,
   similar-looking name, different quantity — and carrying one across to the
   other is silent.

Corrected, total went 0.61 → 0.92-0.95 and the phase snapped to quadrature.

#### 📌 L_π = 19.40 µm gives step 6 a predicted starting length — ARGUED

The device sits at **mm_length / L_π = 1.4434**. Since L_π ∝ W², a narrower
interior has a *shorter* L_π, so preserving the same accumulated imaging
phase over a **physically shorter** body is exactly the premise the waist is
built on. Equating accumulated phase gives a predicted design length rather
than a blind sweep:

```
L_straight / W²  =  ∫ dx / W(x)²      (over the 5-segment profile)
```

⚠ **Label this ARGUED and check it, do not trust it.** The 1×4's directly
analogous argued formula (`L = L_straight − Lt·(W/Ws − 1)`) was **REFUTED
outright** — at its predicted length that device lost 44% of its power. Use
this to *centre* step 6's length window, never to replace measuring it.

**Cost.** Two short jobs (~13 min total).

### Step 5 — Convergence for the waisted body (cheap-to-moderate)

Two axes, neither inherited. **5a is a hard gate on every later number in
this doc and must run before step 6, not alongside it.**

#### 5a — Transition slice count (`C_TRANS`) convergence ⚠ GATING

**This is the single most load-bearing convergence axis in the plan**, for a
measured reason: on the 1×4, doubling `C_TRANS` 6→12 **visibly changed the
field character** — job 12446756 turned a single spurious refocus into a
progressive multi-lobe interference pattern. That was a large qualitative
change from a pure discretisation parameter, on the same physical device.
It did not rescue that device, but it proves the transition region is
**under-resolved at counts people pick by eye**, and every waisted result is
downstream of it.

The transition is also exactly where the physics under test lives: the
leading hypothesis for the 1×4's failure is **radiative loss at the width
transitions** (see the top of this doc). A transition that is
under-discretised is a *staircase*, and a staircase radiates for numerical
reasons that have nothing to do with the real sidewall. **Measuring
transition loss on an unconverged transition would confirm the hypothesis
for the wrong reason** — the exact trap `mmi_1x4_sin300_validation.md`
records for the EME taper staircase in its own check 1.

**Deliverable.** A converged `C_TRANS` for this device, and the S-matrix
plateau table that establishes it.

**Procedure, pre-registered.**
- Sweep `C_TRANS ∈ {2, 4, 6, 8, 12, 16, 24}` at a **fixed** geometry — use a
  mid-aggressiveness profile from step 6's intended grid, not the gentlest
  (a gentle waist under-stresses the very region being resolved).
- Hold `C_STUB`/`C_WAIST` fixed while sweeping `C_TRANS`, so the axis is
  isolated. The stub and waist are *uniform*, so they are not the suspect;
  spot-check them once at the end rather than sweeping all three at once.
- Score on **step 4's full metric** (T_bar, T_cross, excess loss, bar/cross
  phase) — not on totalT alone. A staircase can conserve power while
  rotating phase, and this device's observable is a transfer matrix.
- **Converged** = successive |ΔS| below the tolerance step 4 sets, for **two
  consecutive** refinements, not one. One step can plateau by coincidence.

**Evidence.** The plateau table with the tolerance stated, plus the chosen
count. ⚠ **Report the converged count as a floor tied to the profile it was
measured on.** A more aggressive waist (larger Δwidth over the same `Lt`) is
a steeper ramp and may need more slices — if step 6 explores past the
aggressiveness used here, 5a must be re-run at the edge of that range, not
assumed to carry.

⚠ **Do not reuse the 1×4's answer.** Its `C_TRANS=12` was never shown
converged either — it was one doubling that changed the picture, which
demonstrates non-convergence at 6 rather than convergence at 12.

#### 5b — `num_modes` through the narrowed section

The body *narrows*, so the guided-mode count **falls** in the waist — the
opposite of the 1×4-tapered doc's widening case. The binding constraint is
the **widest** cell (the faces, at `mm_width`), where the straight device's
`N_MODES=14` already applies. ⚠ Heed the 4×4's warning quoted in the 1×4
doc: a probe that hits its own request cap only ever measured a floor, never
found a ceiling — raise the cap until a real ceiling appears.

**Evidence.** Plateau tables for both axes, with the tolerance stated.

**Cost.** Moderate; cheaper than the 1×4's equivalent (smaller domain, fewer
body modes). 5a is ~7 solves at one geometry — cheap next to step 6's grid,
and it is what makes that grid mean anything.

#### ✅ STEP 5a DONE 2026-08-30 — CONVERGED at C_TRANS = 8. Step 6 unblocked.

`examples/verification/python/mmi_2x2_waist_ctrans_convergence.py`, job
**12555775** (res=40, 11 wavelengths, 2h43m). Measured at the steepest
profile in the intended grid — Ww=2.24, Lw=8.0, Lt=4.0, **4.00° sidewall** —
at its own equal-phase length **21.500 µm** (see the centring note below).

| C_TRANS | cells | IL (dB) | imbalance | phase | ΔIL | Δimb | Δphase |
|---|---|---|---|---|---|---|---|
| 2 | 10 | 0.48805 | 0.04914 | −91.300 | — | — | — |
| 4 | 14 | 0.40951 | 0.09024 | −90.675 | 7.85e-02 | 4.11e-02 | 0.625 |
| 6 | 18 | 0.41773 | 0.09441 | −90.643 | 8.22e-03 | 4.17e-03 | 0.032 |
| **8** | 22 | 0.41851 | 0.09570 | −90.633 | **7.72e-04** | **1.29e-03** | **0.010** |
| 12 | 30 | 0.41732 | 0.09674 | −90.624 | 1.18e-03 | 1.04e-03 | 0.008 |
| 16 | 38 | 0.41860 | 0.09720 | −90.613 | 1.28e-03 | 4.61e-04 | 0.011 |
| 24 | 54 | 0.41857 | 0.09732 | −90.613 | 3.17e-05 | 1.15e-04 | 0.000 |

(tolerances: 0.010 dB / 0.0069 / 0.20° — one tenth of step 4's design margins)

**Converged at C_TRANS = 8** by the pre-registered two-consecutive rule.
Worth noting the sweep was carried out to 24 rather than stopped at first
satisfaction: **five** consecutive refinements stay inside tolerance and
nothing re-diverges, with the final step at 3.17e-05. A rule that stops at
the first plateau cannot see a later drift; this one was checked for.

⚠ **C_TRANS = 8 is a FLOOR tied to this profile** (4.00° sidewall). It is
conservative for gentler profiles — which is why 5a was deliberately run at
the *steepest* point in the intended grid rather than a mid or gentle one. If
step 6 ever explores a steeper ramp than 4.00°, re-run 5a at that edge.

**The 1×4's C_TRANS=12 is now contextualised, not inherited.** That number
came from a single doubling that changed the field character — which
demonstrates non-convergence at 6, not convergence at 12. This device
converges at 8 on its own measurement.

#### ⚠ The centring fix — a misleading result that pointed the *expected* way

The first 5a smoke (job 12555720) ran at the straight control's own
`mm_length = 28.0` and reported **2.3–2.5 dB excess loss with the phase 38°
off quadrature**. Read naively that says *"the waist fails on the 2×2 too"* —
a tidy confirmation of the 1×4 result, and therefore the easiest kind of
wrong answer to accept.

It was an artefact of the operating point. **L_π ∝ W², so a narrowed interior
images SHORTER**: the accumulated-phase integral overshoots by **23.2%** at
L=28.0, i.e. the device was being solved **6.5 µm past where it images**.
Re-run at the equal-phase length of 21.500 µm, the identical profile gives
**0.43–0.50 dB and phase −90.2° to −90.8°** — back at quadrature.

The script now auto-centres on the equal-phase prediction, so convergence is
measured near where step 6 will actually operate rather than at a point it
will never use. Predicted lengths across the intended grid (all clearing the
`2·Lt + Lw ≤ mm_length` guard):

| Ww/W | Lw | Lt | L* | saving | total device |
|---|---|---|---|---|---|
| 0.90 | 6.0 | 4.0 | 25.70 | 2.30 | 34.70 |
| 0.90 | 8.0 | 6.0 | 24.79 | 3.21 | 33.79 |
| 0.80 | 6.0 | 4.0 | 22.63 | 5.38 | 31.62 |
| 0.80 | 8.0 | 4.0 | 21.50 | 6.50 | 30.50 |
| 0.80 | 8.0 | 6.0 | 20.50 | 7.50 | 29.50 |

⚠ **Still ARGUED.** The 1×4's directly analogous formula was refuted outright.
This centres step 6's window; it does not replace measuring it.

#### 📌 Early signal, NOT a result

At the converged C_TRANS=8 this profile sits at IL 0.4185 dB, imbalance
0.0957, phase −90.63°, total device 30.50 µm — which would clear **all four**
of step 4's criteria (≤0.4809 dB, ≤0.1371, within 2.0°, <37.0 µm). That is
one point at one profile, measured for convergence rather than as a design
result, and step 6 is what actually tests it. But it is the first concrete
indication the 2×2 may behave differently from the 1×4 — which is the whole
reason this device was chosen.

### Step 6 — Find a design point, with a PRE-REGISTERED STOP RULE

**Deliverable.** Either a (Ww, Lw, Lt, mm_length) point that beats the straight
2×2 on step 4's criterion, or a decisive negative.

⚠ **BLOCKED ON STEP 5a.** Do not start this grid until the transition slice
count is converged. Every point here is measured through the transition, and
the mechanism under test (radiative loss at the transitions) is precisely
what an under-resolved transition fakes. A negative produced on an
unconverged `C_TRANS` would be uninterpretable — it could not distinguish
"the profile fails" from "the staircase radiates."

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

#### 📌 Standing design constraint (user, 2026-08-31): ONLY the multimode region is optimized

**The ports are frozen at the standard 2×2's values and are never swept** —
`waveguide_width=0.45`, `taper_width=1.3`, `taper_length=4.0`,
`waveguide_offset=1.5`, straight `0.5`, `num_modes=[1,7,15,7,1]`,
15 cells per port taper. The hybrid variant's FDTD sections must likewise be
identical to the standard device's. Only `body_waist_width` /
`body_waist_length` / `body_transition_length` / `mm_length` vary.

**Why it matters beyond tidiness:** it is what makes the comparison a
controlled one. Any performance difference is then attributable to the body
profile rather than to a re-tuned port, which is exactly the isolation the
mechanism test needs.

**✅ MEASURED 2026-08-31 — the port hardware really is invariant.**
`geometry.section_gds()` re-origins each section to `x_start=0`, so the input
FDTD section's polygons come out **byte-identical (same hash)** across
`mm_length` 28.0 vs 21.5 *and* across straight vs waisted. Verified directly
rather than assumed.

**⚠ The geometric constraint this imposes on the design grid.** The hybrid
folds `body_stub = 0.25 µm` of MMI body into each FDTD section
(`mmi_2x2_hybrid.py:195`). For those sections to stay identical, the waist
transition must not reach into that stub, so the uniform run before the
transition begins,

```
s = (mm_length − 2·Lt − Lw) / 2   must be  ≥ 0.25 µm
```

This is **strictly tighter** than the constructor's own `2·Lt + Lw ≤
mm_length` guard (which only demands `s ≥ 0`). Now enforced in
`mmi_2x2_waist_design_sweep.py`, which skips and reports any violating point.
Checked against the stage-1 grid: 2 of 40 points violate it, both with
*negative* stub — so they were already being skipped, but by luck rather than
by design, and one point (Ww=2.24, Lw=8.0, Lt=6.0 at L*=20.50) sits exactly
on the 0.25 boundary.

#### ⚠ Cost finding for the hybrid: identical geometry does NOT by itself give cache reuse

Worth knowing before step 9 is budgeted. The section fingerprint
(`build_run_config`, `beef/cache/config.py:287`) includes **`call_kwargs`,
and `x0`/`x1` are among them**. Those are absolute device-frame coordinates:
the input section runs `[−half_total, −half_body + 0.25]`, which is
`x0 = −18.000` at `mm_length=28.0` but `−14.750` at `21.5`. So two runs whose
section polygons hash identically still land on **different cache entries**
whenever `mm_length` differs.

**Consequence, and the fix.** Do *not* rebuild the geometry per swept length
for the hybrid. Follow the standard 2×2's own trick
(`mmi_2x2_hybrid.py:206-213`): keep the geometry at a **fixed reference
`mm_length`**, solve the FDTD sections **once**, and vary only the EME body
section. The waisted body cannot use that device's *analytic* length
substitution (its body is not a length-independent uniform slab, so each
profile/length needs its own body EME solve) — but the body solve is the
cheap half. The expensive FDTD sections stay at one fingerprint.

**And a real payoff if the reference is chosen as `mm_length = 28.0`:** the
sections' `x0`/`x1` then match the *standard 2×2's own completed runs*
exactly, making them genuine cache hits rather than fresh solves. Given that
device's calibration cost (~3.4–4.1 h, and FDTD sections are the expensive
part — `resolved/mmi_2x2_validation.md`), that is the difference between a
cheap hybrid variant and an expensive one. **Verify the hit before relying on
it** (compare fingerprints, or watch for time-stepping in the log) — a
"should be a cache hit" that silently misses is exactly the failure mode
`cache-only-probe-still-writes-outputs` and the 1×4's superset-cache probe
both record.

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
