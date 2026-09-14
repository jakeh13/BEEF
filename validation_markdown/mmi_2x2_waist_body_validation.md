# 2×2 MMI with a WAIST body — build & validation plan

Opened **2026-08-29**, at the user's request: apply the interior-waist body
profile (`body_waist_width`/`body_waist_length`/`body_transition_length`,
added 2026-08-27/28 for the 1×4) to the **2×2 MMI on SOI_220**.

> ## ✅ STATUS 2026-09-14 — the device is VALIDATED. Read this before anything below.
>
> **Design point: `Ww = 2.24, Lw = 6.0, Lt = 5.0, L = 22.0` — 31.00 µm total,
> 16.2 % shorter than the 37.00 µm straight control.**
> Full-device FDTD (51 wl, res 40): IL **0.2010 dB** · imbalance **0.0124** ·
> phase mean-dev **0.267°** · MAX-dev **0.447°** · min total T 0.9548.
> It passes all four criteria **under both the mean AND the max phase reading**,
> so it does not depend on the mean-vs-max ambiguity in the pre-registered text
> — the ambiguity that decided, and killed, the original winner.
>
> **The original 29.50 µm / 20.3 % headline is RETRACTED** (fails on both
> settled taper bases and on full-device FDTD; also never hybrid-buildable at
> s = 0.25 µm). The 30.50 µm interim point passes the mean reading but FAILS
> the max one.
>
> Grid re-scored at a settled basis and verified stable across two independent
> bases (660/668 verdicts agree, 70 passers). The 1×4's waist failure does
> **not** reproduce here: this device works.
>
> **Method status (separate from the design claim).** The hybrid FDTD/EME
> cascade reproduces full-device FDTD to **0.76 % on the bar port** at the
> design point — the calibrated straight device's quality, on a body that
> cannot be β-calibrated. A **cross-port deficit of ~1–2 %** persists on both
> geometries tested; five hypotheses for it were measured and four falsified
> (see "FIVE HYPOTHESES TESTED"). **None of this bears on the design claim** —
> check 4 is single-domain brute force and uses no EME.
>
> ---
>
> *(original framing, kept for provenance)* **This doc opened as a PLAN, not a
> gap list.** At the time nothing had been simulated for a waisted 2×2: every
> number for the *straight* 2×2 was **MEASURED** and cited to
> `resolved/mmi_2x2_validation.md`, while everything about a *waisted* 2×2 was
> **ARGUED**, a prediction to be checked. It has since converted to a normal
> gap list against the [four checks](README.md#the-four-checks).
>
> **Read `mmi_1x4_tapered_body_validation.md` first.** The same profile was
> tried on the 1×4 and **failed outright** — twice closed, on three Ww values,
> doubled taper resolution, and a length search past 2× the device's own
> length. This plan exists because the 2×2 differs in the specific ways that
> would discriminate *why* it failed, not because the 1×4 result is being
> ignored.


---

## 🤝 HANDOFF — state as of 2026-09-14

**NOTHING IS RUNNING.** No queued or running SLURM jobs for this device, no
background processes. Anything in the queue under this account belongs to the
1×4 C-band and 4×4 calibration studies, not to this one. A new session can
start clean.

### Where to pick up — open items, highest value first

1. **Non-uniform ramp cells — the one experiment that would close the method
   question.** The staircase attribution (each EME cell boundary sheds ~2 % of
   |E|², so MORE cells means MORE spurious loss) is supported by correlation,
   not causation. `EMESolver` already accepts `region_cell_lengths=`, measured
   INERT when it was added, so the test is cheap: at a FIXED total cell count,
   redistribute cells toward where the width changes fastest and see whether IL
   moves. If it does, placement matters and there is a real fix; if it does
   not, the error is purely per-boundary count. Either answer is publishable
   and it needs no new FDTD.

2. **The CROSS-port deficit is the remaining method error, and it is narrow.**
   ~1–2 % absolute on BOTH geometries tested, while the bar port reaches 0.76 %
   on the design point. It is NOT the parity shortcut (removing it made things
   worse), NOT the mode count, NOT C_TRANS convergence, and only marginally β
   (calibration moved it 0.0171 → 0.0153). Anything new should target the cross
   path specifically rather than re-testing the four falsified hypotheses —
   see "FIVE HYPOTHESES TESTED" before spending compute.

3. **Optional: extend the calibrated comparison to all 11 wavelengths.** Only
   band centre was run with per-cell β calibration. Mostly confirmatory now
   that the mechanism is understood (uniform offset cancels in a bar−cross
   difference), so low priority.

### What a new session must NOT redo

- ❌ Do not re-run the C_TRANS ladder hoping for convergence. Measured to ct32
  on two mode bases; IL diverges by construction. **C_TRANS is an OPTIMUM, not
  a floor.**
- ❌ Do not raise `N_MODES` above 10 on this body. The 2.24 µm waist supports
  exactly 10 distinct modes; N=11 is an EXACT duplicate (overlap 1.000000).
  14 was inherited from the 2.8 µm face and cost a 15× runtime penalty.
- ❌ Do not quote the 29.50 µm / 20.3 % headline. Retracted — fails on two
  settled bases and on full-device FDTD, and was never hybrid-buildable.
- ❌ Do not trust `MMI2X2_HY_REUSE=0` output produced before 2026-09-13. That
  branch returned a copy of the BAR result.

### Traps this study paid for — apply them, do not rediscover them

- **Filename collisions, eight instances.** An output name must carry every
  parameter the result depends on, **geometry first**. The most recent: the
  hybrid comparison was keyed on `res/ct/nm` only, so merging a second design
  point silently overwrote the first.
- **Read the reference's own `config.json`; never assume a script default
  reproduces it.** Cost 4h15m of 48-rank time — the check-4 references are at
  `decay_tol=1e-05`, the script default is `1e-11`, which is below the
  reachable floor and never converges.
- **Bit-identical output from a changed parameter is a BUG SIGNAL**, not a
  small effect. That is how the dead seam-modes knob was found.
- **One geometry cannot characterise a method.** The "systematic bar-high /
  cross-low bias" held across N_MODES, C_TRANS, parity and calibration
  variations — all on one body — and reversed sign on the first new body tried.
- **Divide per-rank log counts by the rank count** before treating them as a
  measurement.

### Key artifacts

| what | where |
|---|---|
| all outputs for this device | `example_output/mmi_2x2_tapered/` |
| two-geometry comparison figure | `spectra_hybrid_vs_fdtd_two_geometries.png` |
| per-geometry spectra | `spectra_*_hybrid_vs_fdtd_Ww2.24_Lw6_Lt{5,6}_L{22,21.5}.png` |
| stitched hybrid field + its replay cache | `hybrid_stitched_res40/` (`recon_cache_*.npz` re-renders in seconds, no re-solve) |
| β(width) calibration curve | `beta_width_curve_res40_nm10_nwl11_nw9.json` |
| full-device FDTD references | `check4_*.json`, `smatrix_waisted_*.pkl` |
| re-scored grid rows | `example_output/verification/mmi_2x2_waist_free/free_sweep_rows_*_nt{7,13,21}_tc*_nb15_*.jsonl` |

### Library changes landed this session (shared — other studies depend on them)

`EMESolver.solve_section(cell_neff_override=...)` and
`HybridSolver(cell_neff_fn=...)` add **per-cell β calibration for TAPERED
bodies**, which previously could not be calibrated at all. Backward
compatibility is *tested*: the new cache-fingerprint key is added only when the
override is non-None, with a test asserting an uncalibrated config does not
contain it — otherwise every cached `solve_section` entry in the repo would
have been invalidated, including the FDTD sections the modal-mux, 1×4 and 4×4
studies depend on. Full suite **814 passed**; 13 new tests in
`tests/test_eme_cell_neff_override.py`.

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

#### ❌ STEP 5b DONE 2026-09-02 — job 12780420 (3:38:23) — AND IT INVALIDATES STEP 6's SCORING

**Read this before trusting any step-6 number.** Two independent 1-D sweeps at
the design point (`Ww=2.24, Lw=8.0, Lt=6.0, L=20.5`), res 40, 11 wl, each point
streamed to disk.

##### Axis A — `N_BODY`: CONVERGED at 15, and my truncation hypothesis is REFUTED

| N_BODY | 12 | 15 | 18 | 21 | 24 | 27 |
|---|---|---|---|---|---|---|
| IL dB | 0.3881 | 0.3910 | 0.3934 | 0.3931 | 0.3928 | 0.3930 |

Span **0.0053 dB**; plateau from 18, last four within **0.0006 dB**. A real
plateau, not a cap artifact (the 4×4's warning) — it flattens well below the 27
tried. **`N_BODY=15` was adequate**, so step 6's body basis is fine.

I had predicted that if the check-4 EME-vs-FDTD loss offset were basis
truncation, raising the body count would drive IL toward FDTD. It does not —
it moves *away* by 0.005 dB. **Prediction falsified on measurement.**

##### ❌ Axis B — `N_TAPER`: NOT CONVERGED, and 7 is a lucky point

| N_TAPER | IL dB | gap to FDTD | phase° | mean-dev° | verdict |
|---|---|---|---|---|---|
| 5 | 0.4676 | +0.264 | −91.862 | 2.112 | FAIL |
| **7 (what step 6 used)** | 0.3910 | +0.187 | −91.446 | **1.696** | **PASS** |
| 9 | 0.3658 | +0.162 | −91.803 | 2.053 | FAIL |
| 11 | 0.3459 | +0.142 | −91.939 | 2.189 | FAIL |
| 13 | 0.3536 | +0.150 | −92.025 | 2.275 | FAIL |

Span **0.1217 dB — 23× the body axis.** Loss falls steeply then scatters around
~0.35 dB without a clean plateau, and never reaches FDTD's 0.2036 dB.

> ### ❌ THE DESIGN POINT PASSES THE PHASE CRITERION ONLY AT `N_TAPER = 7`
> It **fails at 5, 9, 11 and 13.** `7` sits between failures, and its phase
> (−91.446°) breaks a sequence that otherwise brackets the check-4 FDTD
> reference (**−91.849°**) from both sides. That is the profile of a lucky
> point on an unconverged axis, not of a converged one.

**Consequences, stated plainly:**

1. **Step 6's entire 56-passer set was scored at `N_TAPER=7`.** The ranking,
   the design point and the 20.3% headline all rest on a basis measured here as
   unconverged. **The 29.50 µm point must not be treated as validated**, and
   step 6 needs re-scoring at a converged taper basis before it can be closed.
2. **This corroborates check 4 independently.** Full-device FDTD gave a phase
   mean-dev of **2.099°**; a better-resolved EME gives **2.05–2.28°**. Two
   different routes now agree the phase margin is not there, and `N_TAPER=7`
   was the single configuration that hid it.
3. **Neither mode axis explains the loss offset.** At the best taper basis EME
   still sits **~0.15 dB above FDTD**, so ~80% of the check-4 gap has another
   cause — candidates: `C_TRANS` interaction, the 15-cell taper
   discretisation, or something structural in the cascade. **Open.**

⚠ **`N_TAPER` is not converged even at 13** (11→13 moves +0.0077 dB, having
moved −0.0199 before that). The axis needs extending past 13, and the taper
cell count (`TAPER_CELLS=15`) should be swept with it — more modes through the
same 15 cells may itself be the limiter. **Do not read 11 or 13 as "the
converged answer".**

#### ▶ IN FLIGHT 2026-09-10 — STEP 5b PART 2, as the prerequisite for re-scoring step 6

Asked to re-score step 6. **It cannot be validly re-scored yet, and the
reason is in 5b part 1's own table above**: `N_TAPER` is unconverged even at
13 (11→13 moves +0.0077 dB having moved −0.0199 before it), so picking any
particular value to re-score at — 13 included — would repeat step 6's
original mistake at a different number. 5b part 1 also left `TAPER_CELLS`
pinned at 15 throughout, so *more modes through the same 15 cells* may itself
be the limiter; that is an interaction, and independent sweeps cannot test it.

Job **13058925**, array 0-9, one `(N_TAPER, TAPER_CELLS)` pair per task via
the script's existing `--taper-convergence` mode:

- **edge B'** (tasks 0-3): N_TAPER 15, 17, 19, 21 at TAPER_CELLS=15 — does the
  mode axis plateau once pushed past 13?
- **edge C** (tasks 4-7): TAPER_CELLS 10, 20, 25, 30 at N_TAPER=13 — does the
  cell axis matter at all?
- **diagonal** (tasks 8-9): 19×25 and 21×30 — is it the interaction?

Held at 5b part 1's own design point (`Ww=2.24 Lw=8.0 Lt=6.0 L=20.5`) so these
points EXTEND that table instead of starting a second one. The port taper
(0.45→1.3 over 4.0 µm) is frozen and identical at every step-6 candidate, so a
taper basis that plateaus here plateaus for the grid.

⚠ **One-pair-per-task was required, not a convenience.** The mode wrote to
`taper_convergence_res{RES}.json` — res only, no geometry, no pair — so an
array would have had all ten tasks clobbering one file, the same defect this
study has already hit six times (most recently the geometry-less `check4`
name). Filename now carries geometry + pair + `nb`/`nwl`; assemble by
globbing. The mode already rewrites its JSON after every pair, so a task that
times out still leaves finished work on disk.

**Then, and only then**, re-score the 668-point grid at the plateaued basis.

##### ✅ EDGE B' IS ALREADY DONE — recovered, not re-run

The earlier killed 5b-part-2 attempt had in fact **completed all four edge-B'
points** before it was stopped; they were sitting in the res-only
`taper_convergence_res40.json` (the very filename defect noted above, which is
why they were easy to miss). Recovered 2026-09-10 and folded in with 5b part
1's own points, at the same res 40 / 11 wl / N_BODY=15 / design point:

| N_TAPER (cells=15) | IL dB | ΔIL vs prev | mean-dev° | verdict |
|---|---|---|---|---|
| 5 | 0.4676 | — | 2.112 | FAIL |
| **7 (what step 6 used)** | 0.3910 | −0.0766 | **1.696** | **PASS** |
| 9 | 0.3658 | −0.0252 | 2.053 | FAIL |
| 11 | 0.3459 | −0.0199 | 2.189 | FAIL |
| 13 | 0.3536 | +0.0077 | 2.275 | FAIL |
| 15 | 0.3542 | +0.0006 | 2.253 | FAIL |
| 17 | 0.3612 | +0.0070 | 2.304 | FAIL |
| 19 | 0.3693 | +0.0081 | 2.261 | FAIL |
| 21 | 0.3706 | +0.0013 | 2.239 | FAIL |

**The mode axis has essentially settled past 13.** Span 13→21 is
**0.0170 dB**, against **0.1217 dB** over 5→13 — a 7× reduction. It still
creeps upward rather than being flat, so "converged" is too strong; but it is
no longer the dominant uncertainty, and it moves AWAY from FDTD (0.2036 dB),
not toward it.

> ### ❌ THE DESIGN POINT FAILS THE PHASE CRITERION AT EVERY N_TAPER FROM 9 TO 21
> mean-dev sits in a tight **2.05°–2.30°** band across nine mode counts, against
> a **≤2.0°** criterion. `N_TAPER=7`'s 1.696° is an outlier low by ~0.5°, off
> a cluster that never comes back under threshold. This is now the THIRD
> independent route to the same verdict: converged EME says 2.24–2.30°,
> full-device FDTD (check 4) says **2.099°**, and only the single `N_TAPER=7`
> configuration step 6 happened to use ever said otherwise.

**Consequence: the 29.50 µm design point is not recoverable by refining the
taper basis, and step 6's 20.3%-shorter headline does not survive.** The
re-score is therefore not "re-rank the same 56 passers" — it has to answer
whether ANY grid point clears the phase criterion on a settled basis.

##### ⚠ The CELL axis matters more than the mode axis past 13

First edge-C point (job 13058925 task 4, the one shard that survived the OOM):

| N_TAPER=13 | TAPER_CELLS | IL dB | mean-dev° |
|---|---|---|---|
| | 10 | 0.3900 | 2.231 |
| | 15 | 0.3536 | 2.275 |

Dropping 15→10 cells moves IL **+0.0364 dB** — more than twice the entire
13→21 mode-axis span (0.0170 dB). So `TAPER_CELLS` is the live axis now, and
whether 15 is itself converged is unanswered. Job **13059401** (array 0-4)
covers the remaining `13:20`, `13:25`, `13:30` and the diagonals `19:25`,
`21:30`.

##### ✅ EDGE C DONE 2026-09-10 (job 13059401) — and mean-dev is BASIS-INVARIANT

| N_TAPER | CELLS | IL dB | ΔIL vs prev | imb | phase° | mean-dev° |
|---|---|---|---|---|---|---|
| 13 | 10 | 0.3900 | — | 0.0796 | −91.981 | 2.231 |
| 13 | 15 | 0.3536 | −0.0364 | — | — | 2.275 |
| 13 | 20 | 0.3327 | −0.0209 | 0.0807 | −91.987 | 2.238 |
| 13 | 25 | 0.3302 | −0.0025 | 0.0812 | −91.987 | 2.237 |
| 13 | 30 | 0.3220 | −0.0082 | 0.0812 | −91.977 | 2.227 |
| 19 | 25 | 0.3511 | — | 0.0777 | −91.929 | 2.180 |

Plus the two diagonals, which turn out to be the best points measured:

| N_TAPER × CELLS | IL dB | mean-dev° |
|---|---|---|
| 19 × 25 | 0.3511 | 2.180 |
| 21 × 30 | 0.3436 | **2.160** |

**The decisive observation is the last column.** Over all 13 settled points
(`N_TAPER` ≥ 9, cells 10→30) mean-dev spans only **2.053°–2.304°** — a 0.251°
spread — and **every single one FAILS the 2.0° criterion**, while IL swings
0.068 dB across the same space. `N_TAPER=7`'s 1.696° is the sole outlier in
the entire 15-point map.

⚠ **"Basis-invariant" is slightly too strong, and the diagonal is why.**
Raising BOTH axes together does help — the diagonal reaches 2.160/2.180,
below every edge point (edges give 2.23–2.30) — so 5b part 1's interaction
hypothesis is partly borne out, and the diagonal has **not flattened**
(19×25 → 21×30 still moves −0.020°, with IL still falling −0.008 dB). So the
honest claim is: **the design point fails phase on every basis actually
measured, by a margin (0.16°) far larger than the drift left in the axis**,
and closing the remaining 0.16° by refinement alone would need an
extrapolation this study's own record warns against (convergence ratios
stall — see the HCG lesson). Not "proven physical", but not a basis artifact
either.
IL, by contrast, is still creeping (cells 25→30 moves −0.008 dB, modes 13→21
+0.017 dB) and never approaches FDTD's 0.2036 dB — best here is 0.3220 at
13×30, still **+0.118 dB** high, so 5b part 1's unexplained loss offset
stands.

#### 📋 THE HYBRID-vs-FDTD GAP: FIVE HYPOTHESES TESTED, 2026-09-08 .. 09-12

Consolidated so a reader does not have to reconstruct the sequence. Each row
was a real measurement, not an argument, and four of five were falsified.

| # | hypothesis | test | verdict |
|---|---|---|---|
| 1 | rank-deficient interface basis (`N_MODES=14` on a 2.24 µm waist) | re-ran the whole ct4..ct32 ladder at the rank-clean `N_MODES=10` | ❌ **identical ladder** — IL 0.2415 at ct32 on both bases, agreeing to ~0.004 dB. Real defect (15× runtime penalty) but not the cause. |
| 2 | `C_TRANS` not converged | ladder to ct32 on two bases | ❌ **it is the STAIRCASE** — IL diverges because each added cell is another artificial step. C_TRANS has an OPTIMUM, not a floor. |
| 3 | `MODE_SIGNS` parity shortcut (cross port derived, not simulated) | `MMI2X2_HY_REUSE=0`, both arms actually solved | ❌ **made it WORSE** — imbalance 0.0339 → 0.0460 against FDTD's 0.0168. The shortcut was masking error, not causing it. |
| 4 | interface mode count starved at the FDTD↔EME seams | `interface_modes = [14,10,10,10,10,14]` | ❌ **bit-identical output** — and the reason is structural: section 1 has `n_out = interface_modes[1] = 10`, so the extra seam modes are discarded at the very next interface. The corrected experiment (14 everywhere) is just the nm14 run, which gives imbalance 0.0333 — still 2× FDTD. |
| 5 | uncalibrated body β | built per-cell calibration (below) and ran it | ⚠️ **PARTIAL** — improves all three axes but closes only ~11–27% of the gap. |

##### ⚠ Two traps inside those tests, both worth more than the tests themselves

**`MMI2X2_HY_REUSE=0` had NEVER produced a valid cross result on this device.**
The branch set `cross_kwargs = {0: {"input_idx": 0}, OUT: {"output_idx": 0}}`
— excite the BOTTOM arm, read the BOTTOM output — which on a y-mirror-symmetric
device is the BAR path reflected. It burned 1896 s of real FDTD and returned an
exact copy of bar (`T_bar == T_cross` to 4 dp, imbalance 0.0000, phase −0.001°).
Correcting only the `input_idx` then hits a hard refusal from the library:

```
ValueError: mirror_pairs[6]=0: center_at mismatch at x=10.85 vs
mirrored x=-10.85: -0.75 vs 0.75 um
```

because a plain x-mirror maps output onto input WITHOUT flipping y, and is only
legal when both ports share a y. So `input_idx=0` was not a typo — it was the
only way to keep that mirror legal, and it silently changed WHICH QUANTITY was
computed. The honest configuration is no mirror at all (`mirror_pairs={}`),
solving the output section directly. The script's own docstring offers
`SYMMETRIC_MODE_REUSE=0` as the way to "remove the assumption"; anyone who took
it up before this fix got a silently wrong answer.

**A knob that looks wired up but is discarded downstream.** The seam-modes test
produced bit-identical coefficients, which was the tell — genuine physics
rarely reproduces to the last digit. Checking `_mode_lists()` showed the extra
modes are truncated at the adjacent interface. **Bit-identical output from a
changed parameter means the parameter did not take effect**; treat it as a bug
signal, never as "the effect is small".

#### ✅ PER-CELL β CALIBRATION FOR A TAPERED BODY — new capability, 2026-09-11/12

`calibrate_body_neff` is what closed the engine-agreement gap on every
UNIFORM-body device here, and it is hard-refused for a tapered one for two
independent reasons:

1. **Mechanically** it takes `solve_section()`'s `neff_override` fast path,
   which builds a DIAGONAL `exp(2jπ·n_eff·L/λ)` and *skips the mode solve and
   `compute_s_matrix()` outright*. Diagonal means NO mode mixing — exact for a
   slab, and it deletes the entire physics of a taper.
2. **Practically** an intermediate ramp width has no neighbouring FDTD section
   to take β from, and `section_mode_beta()` measures β by fitting dφ/dL — a
   0.75 µm ramp slice is far too short to fit a slope on.

**Both are now solved.** `EMESolver.solve_section(cell_neff_override=...)` and
`HybridSolver(cell_neff_fn=...)` keep the multi-cell solve and every interface
S-matrix and substitute n_eff ONLY in the propagation within each cell. The
split is exact: interface overlaps integrate the mode FIELDS, propagation uses
n_eff, so calibrating one cannot disturb the other. β comes from a measured
β(width) curve (`measure_beta_width_curve.py`, 9 widths 2.24–2.80 µm on uniform
`UniformCrossSection` strips at the hybrid's own FDTD settings), interpolated
per cell by `beta_width_curve.py` — legitimate because each EME cell IS a
uniform strip.

**The offset being corrected is large, and independently confirms the 4×4's
finding.** FDTD minus mode-solver n_eff at res 40:

| width | modes 0–4 |
|---|---|
| 2.24 µm | +0.0417 … +0.0481 |
| 2.80 µm | +0.0417 … +0.0445 |

against the 4×4 study's **+0.0462** — a different device, wavelength (1.5475 vs
1.31) and material. **Grid dispersion is a property of the METHOD at a given
resolution, not of one device.**

**But the benefit here is modest, and the reason is quantitative.** At band
centre, calibrated vs uncalibrated vs FDTD:

| | uncalibrated | calibrated | full FDTD |
|---|---|---|---|
| excess loss dB | 0.1875 | 0.1877 | 0.2178 |
| imbalance | 0.0339 | **0.0321** | 0.0168 |
| phase ° | −90.953 | **−91.144** | −91.661 |

Every axis improves (phase gap 0.708 → 0.517°, imbalance 0.0171 → 0.0153) but
none closes. The offset is nearly UNIFORM across modes — spread only ~0.0006
among the dominant modes 0–2 — and these metrics are DIFFERENCES between bar
and cross, where a uniform offset cancels. Combined with a short 21.5 µm body,
a large correction yields a small net change. Expect the opposite where the 4×4
saw it: long bodies, absolute phase.

##### Two implementation traps, both regression-tested

- **NaN in a cache fingerprint makes the cache miss forever.** `json` writes a
  bare `NaN` and reads it back as `float('nan')`, and `NaN != NaN`, so
  `store.lookup()`'s `cached_config != config` collision guard can never match.
  Measured: the first calibrated run logged "hash match … but config.json
  differs" on all 8 ranks. NaN is now a string sentinel.
- **Mode IDENTITY is not stable across a β(width) curve.** The solver sorts by
  n_eff, so crossing modes swap indices — m6/m7 exchange between w=2.66 and
  2.73 µm here, and m8/m9 are near-degenerate at 2.45. Interpolating such a
  mode blends two physical modes. `well_ordered_modes()` detects this and
  leaves those modes uncalibrated via NaN rather than forcing all-or-nothing.
  It selects **modes 0–4 (99.70 % of the power)**; the 0.47 % remainder keeps
  its solved n_eff.

⚠ Backward compatibility is *tested*, not asserted: the new fingerprint key is
added ONLY when the override is non-None, with a test proving an uncalibrated
config does not contain it. Otherwise every cached `solve_section` entry in the
repo — the expensive FDTD sections the modal-mux, 1×4 and 4×4 studies depend on
— would have been invalidated at once. Full suite: **814 passed**.

#### ✅ HYBRID MEASURED ON A SECOND GEOMETRY 2026-09-13 — job 13153984, and it CORRECTS an earlier over-generalisation

Every hybrid-vs-FDTD conclusion in this file up to now came from ONE body
(`Lw6 Lt6 L21.5`). The confirmed 31.00 µm design point gives a second,
independent geometry, which is what separates "a property of the METHOD" from
"a property of that body". Same config both times (nm10, ct8, res 40, 11 wl):

| | NEW `Lw6 Lt5 L22.0` | incumbent `Lw6 Lt6 L21.5` |
|---|---|---|
| mean \|ΔT\| bar | **0.76 %** | 2.91 % |
| mean \|ΔT\| cross | 3.17 % | 3.32 % |
| bar bias range | −0.0087 … **+0.0041** | **+0.0090 … +0.0186** |
| cross bias range | −0.0201 … −0.0078 | −0.0214 … −0.0083 |
| Δ split | +0.36 … +0.99 % | +1.02 … +1.74 % |
| Δ IL dB | +0.0436 | −0.0109 |

❌ **The "systematic bar-high / cross-low bias" recorded earlier does NOT
generalise.** On the new geometry the bar residual straddles zero; only the
CROSS deficit reproduces, at nearly identical magnitude on both bodies. The
correct statement is therefore narrower: **the hybrid under-delivers to the
cross port by ~1–2 % absolute on both geometries, while the bar port's error
is geometry-dependent.** That the bar excess reversed sign between two bodies
is exactly why one geometry was not enough to characterise the method.

✅ **The new geometry's bar port is essentially at the calibrated straight
device's quality** (0.76 % against that device's 0.62 %/0.59 %), on a body
that cannot be β-calibrated at all.

The cross deficit is also the component that survived the other tests: it did
not improve when the parity shortcut was removed (simulating the cross port
properly made it WORSE — 0.0460 vs 0.0339 imbalance), and per-cell β
calibration moved it only from +0.0171 to +0.0153. So the cross path is where
the remaining method error lives, and it is not the parity vector, not the
mode count, and only marginally β.

⚠ Consistent with — but not proof of — the staircase attribution: the new body
has a SHORTER ramp (Lt=5.0 vs 6.0) at the same C_TRANS=8, hence shorter cells
and smaller per-boundary steps. Establishing that causally needs the
non-uniform-cell test (`EMESolver(region_cell_lengths=...)`) at fixed cell
count, which has not been run.

⚠ **Eighth filename collision, same class as the rest.** The comparison output
was named `hybrid_vs_fdtd_res{RES}_ct{C_TRANS}_nm{N_MODES}.json` — no
geometry — so merging the new device silently OVERWROTE the incumbent's
result. Both recovered by re-merging their (geometry-tagged) shards; the name
now carries `Ww/Lw/Lt/L`, and so do the plot filenames and titles.

#### ✅✅ CHECK 4 CONFIRMS A NEW DESIGN POINT 2026-09-12 — job 13124560, and it passes BOTH phase readings

Full-device FDTD (51 wl, res 40, decay_tol 1e-05 — identical config to the two
references already on disk) on step 6's best-margin re-scored candidate:

| candidate | µm | IL dB | imb | mean-dev° | MAX-dev° | saving | verdict |
|---|---|---|---|---|---|---|---|
| **NEW `Lw6 Lt5 L22.0`** | 31.00 | **0.2010** | **0.0124** | **0.2674** | **0.4474** | 16.2% | ✅ **PASS on BOTH readings** |
| incumbent `Lw6 Lt6 L21.5` | 30.50 | 0.2178 | 0.0168 | 1.9110 | 2.5644 | 17.6% | PASS mean / ❌ FAIL max |
| step 6's original `Lw8 Lt6 L20.5` | 29.50 | 0.2036 | 0.0976 | 2.0994 | 2.3230 | 20.3% | ❌ FAIL |

**The new point is better on EVERY axis simultaneously** — lower loss
(−0.0168 dB), better imbalance (−0.0043), higher minimum total transmission
(0.9548 vs 0.9511), and **5.7× the phase margin** — for 0.5 µm more length.

📌 **It also dissolves the criterion ambiguity that has run through this whole
study.** The pre-registered text says MAX-over-band phase deviation while
`score()` applies the band MEAN, and 4 of 56 passers flipped between the two
readings (see the phase-criterion note). The incumbent sits exactly in that
crack: 1.911° mean (pass) against 2.564° max (fail). **The new point passes
both with room to spare — 0.267° mean, 0.447° max, against a 2.0° budget.**
A design that does not depend on which reading governs is worth 1.4
percentage points of length saving.

**EME and FDTD agree on the margin.** Step 6 predicted 0.022° (mean, at
nt21/tc30); FDTD measures 0.267°. Both are an order of magnitude inside the
budget, so the two engines agree this candidate is well-centred — unlike the
original 29.50 µm winner, where EME's `N_TAPER=7` reading (1.696°) and FDTD
(2.099°) disagreed across the threshold and the EME value turned out to be
basis luck.

⚠ **Method note, recorded because it cost 4h15m of 48-rank time.** The first
submission (job 13119819) used the script's `decay_tol` DEFAULT of 1e-11,
having assumed it matched the references. It does not — both references were
measured at **1e-05**, which is in each one's own `config.json`. 1e-11 is
below this study's measured hard floor, so that run was still time-stepping at
t=739 with no decay in sight when it was cancelled. **Read the reference's
stored config; do not assume a script default reproduces it.**

#### ✅ STEP 6 RE-SCORED 2026-09-10 — job 13061906, 8/8 tasks, 668/668 points, no faults

Full grid re-scored at the settled basis **N_TAPER=13 / TAPER_CELLS=30 /
N_BODY=15** (8 tasks by (Ww, Lw), ~1h20m each; basis 2144 s per Ww then
~20 s/point warm, interface memo re-verified |dS| = 0.000e+00 at the new
basis). Result is **POSITIVE, and the opposite of what I predicted.**

**PASSERS: 56 (nt7/tc15) → 70 (nt13/tc30).** Churn: 46 kept, 10 lost, 24
gained. Refining the taper basis ADDS passers rather than removing them.

##### ⚠ TWO DIFFERENT "DESIGN POINTS" — I had been conflating them, and they diverge

| | Ww | Lw | Lt | L | s | total | nt7 verdict | nt13/tc30 verdict |
|---|---|---|---|---|---|---|---|---|
| step 6's original winner, and **5b's subject** | 2.24 | 8.0 | 6.0 | 20.5 | **0.25** | 29.50 | PASS (mdev 1.696) | ❌ **FAIL (2.227)** |
| **the point every hybrid run used** | 2.24 | 6.0 | 6.0 | 21.5 | 1.75 | 30.50 | PASS (mdev 1.505) | ✅ **PASS (1.846)** |

- **Step 6's 29.50 µm headline is dead**, now by direct re-score rather than
  extrapolation: mdev 2.227 > 2.0. That is the point 5b swept across the whole
  2-D basis map and found failing at every settled basis, and the re-score
  agrees. It was *already* excluded from hybrid work anyway — s=0.25 is below
  the s ≥ body_stub+dpml = 0.75 buildability floor.
- **The hybrid work's design point SURVIVES**, and improves on two axes at the
  better basis: IL 0.3880 → **0.3017 dB**, imbalance 0.0127 → **0.0088**, with
  phase mdev 1.505 → 1.846 (still inside 2.0). **So the hybrid-vs-FDTD
  comparison, the C_TRANS sweep and the stitched field maps were all measured
  on a design point that is still valid.**

#### ✅ STAGE 2 CONFIRMS IT 2026-09-10 — job 13065614, 668/668 at N_TAPER=21 / TAPER_CELLS=30

Re-ran the whole grid on the most phase-favourable basis in the 2-D map
(mdev 2.160° at the old point, vs 2.227° at 13×30). Basis cost 2737 s vs
2144 s at 13 modes — only 1.28×, so the full grid was affordable rather than
just a shortlist. Zero rank-deficiency escalations even at 21 modes (the port
taper is wide enough to fill the basis, unlike the waist cells).

| basis | passers |
|---|---|
| nt7 / tc15 (original) | 56 |
| nt13 / tc30 (stage 1) | **70** |
| nt21 / tc30 (stage 2) | **70** |

**The two settled bases agree on 660 of 668 verdicts (98.8%)** — 66 passers in
common, 4 unique to each. So the re-scored passer set is basis-STABLE, not an
artifact of picking 13×30. That is the check stage 2 existed to make, and it
passes.

##### Verdicts across all three bases

| point | s | total | nt7/tc15 | nt13/tc30 | nt21/tc30 |
|---|---|---|---|---|---|
| old 29.50 winner (5b's subject) `Lw8 Lt6 L20.5` | 0.25 | 29.50 | PASS 1.696 | ❌ FAIL 2.227 | ❌ FAIL 2.160 |
| **hybrid design point** `Lw6 Lt6 L21.5` | 1.75 | 30.50 | PASS 1.505 | ✅ PASS 1.846 | ✅ **PASS 1.716** |
| cand C `Lw8 Lt4 L21.5` | 2.75 | 30.50 | PASS 0.883 | PASS 1.206 | ✅ **PASS 1.138** |
| cand A `Lw6 Lt5 L22.0` | 3.00 | 31.00 | PASS 0.209 | PASS 0.103 | ✅ **PASS 0.022** |
| cand B `Lw8 Lt4 L22.0` | 3.00 | 31.00 | PASS 0.680 | PASS 0.187 | ✅ **PASS 0.115** |

- **The 29.50 µm headline is dead on BOTH settled bases** — not a one-basis
  fluke.
- **The hybrid design point passes on all three.** Everything built on it —
  the hybrid-vs-FDTD comparison, the C_TRANS sweep, the stitched field maps —
  stands.
- **Candidate A's phase margin improves monotonically as the basis is refined**
  (0.209 → 0.103 → 0.022°), which is what a genuinely well-centred design
  looks like, as opposed to the old winner whose 1.696° was basis luck.

##### Final ranking on the most phase-favourable basis (buildable, s ≥ 0.75)

| total µm | Lw | Lt | L | s | IL dB | imb | mdev° | note |
|---|---|---|---|---|---|---|---|---|
| 30.50 | 8.0 | 4.0 | 21.5 | 2.75 | 0.3307 | 0.0883 | **1.138** | shortest with real margin |
| 30.50 | 10.0 | 5.0 | 21.5 | 0.75 | 0.4310 | **0.1369** | 1.697 | imb at the 0.1371 limit |
| 30.50 | 6.0 | 6.0 | 21.5 | 1.75 | 0.3278 | **0.0115** | 1.716 | incumbent; best imbalance |
| 30.50 | 8.0 | 6.0 | 21.5 | 0.75 | 0.3691 | 0.1083 | 1.929 | |
| 31.00 | 6.0 | 5.0 | 22.0 | 3.00 | 0.3655 | 0.0126 | **0.022** | best overall balance |
| 31.00 | 8.0 | 4.0 | 22.0 | 3.00 | **0.2660** | 0.0991 | 0.115 | best IL in the grid |

68 of 70 passers are buildable. **Shortest buildable stays 30.50 µm = 17.6%**
shorter than the 37.00 µm control.

📌 **RECOMMENDATION — the design choice is now a real trade, not a search.**
Spending **0.5 µm more length (31.00 µm, 16.2% instead of 17.6%)** buys
`Lw6 Lt5 L22.0`: phase deviation **0.022° against the incumbent's 1.716°** —
a ~78× margin improvement — at comparable imbalance (0.0126 vs 0.0115) and
0.038 dB more loss. Given phase margin is what killed the old winner and what
check 4's full-device FDTD independently flagged, that trade looks strongly
worth taking. `Lw8 Lt4 L22.0` is the alternative if loss dominates: **IL
0.2660 dB, the best in the grid and within 0.062 dB of full-device FDTD**,
with mdev 0.115° but 8× worse imbalance.
⚠ Any switch needs a full-device FDTD confirmation at the new geometry
(check 4) before it replaces the incumbent — the hybrid sections are frozen
around s, and every candidate above has a different s.

##### Stage-1 passer shortlist (13×30, superseded by the table above)

| total µm | Ww | Lw | Lt | L | s | IL dB | imb | mdev° | buildable |
|---|---|---|---|---|---|---|---|---|---|
| 29.50 | 2.24 | 10.0 | 5.0 | 20.5 | 0.25 | 0.3536 | 0.1156 | 1.714 | **no** |
| 30.00 | 2.24 | 10.0 | 5.0 | 21.0 | 0.50 | 0.3496 | 0.1280 | 1.515 | **no** |
| 30.50 | 2.24 | 8.0 | 4.0 | 21.5 | 2.75 | 0.3077 | 0.0893 | 1.206 | yes |
| 30.50 | 2.24 | 6.0 | 5.0 | 21.5 | 2.75 | 0.4649 | 0.0121 | 1.388 | yes |
| 30.50 | 2.24 | 6.0 | 6.0 | 21.5 | 1.75 | 0.3017 | **0.0088** | 1.846 | yes (incumbent) |
| 31.00 | 2.24 | 6.0 | 5.0 | 22.0 | 3.00 | 0.3387 | 0.0124 | **0.103** | yes |
| 31.00 | 2.24 | 8.0 | 4.0 | 22.0 | 3.00 | **0.2379** | 0.0993 | 0.187 | yes |

68 of the 70 passers are hybrid-buildable. **Shortest buildable is 30.50 µm =
17.6% shorter than the 37.00 µm control** (not the retracted 20.3%); 29.50 µm
survives only if buildability is ignored.

**Two candidates stand out well beyond the incumbent**, both at 31.00 µm
(16.2%) for 0.5 µm more length: `Lw6.0 Lt5.0 L22.0` has mdev **0.103°** with
imbalance 0.0124, and `Lw8.0 Lt4.0 L22.0` has **IL 0.2379 dB — the best in the
grid, and within 0.034 dB of the full-device FDTD reference (0.2036)**, which
is notable given 5b's unexplained ~0.15 dB EME-vs-FDTD offset at the old
point. Either would buy far more phase margin than the incumbent's 1.846°.

##### ❌ MY PREDICTION OF A NEGATIVE RESULT WAS WRONG — the design space has real margin

Before re-scoring I predicted step 6 would likely close negative, reasoning
that a ~+0.5° phase correction would wipe out the 56 passers. **Checked
against the 668 rows already on disk, that is false.** The passer set's
mean-dev at `N_TAPER=7` is widely distributed (min **0.022°**, 13 points under
0.25°), not clustered near the 2.0° limit. Applying the measured
**−0.54°** shift uniformly:

- **53 of 56** passers still clear phase; **50** of those are hybrid-buildable
  (s ≥ 0.75)
- 6 previously phase-FAILING points would be promoted UNDER 2.0°

❗ **But the design point every piece of hybrid work used dies, narrowly.**
`Ww2.24 Lw6.0 Lt6.0 L21.5 s1.75` (30.50 µm) goes **md 1.505 → 2.045** — over
by 0.045°, i.e. exactly on the bubble. Its replacements at the same 30.50 µm
have far more phase margin but markedly worse imbalance:

| total µm | Lw | Lt | s | IL dB | imb | md@nt7 | md shifted |
|---|---|---|---|---|---|---|---|
| 30.50 | 6.0 | 6.0 | 1.75 | 0.3880 | **0.0127** | 1.505 | **2.045 FAIL** |
| 30.50 | 8.0 | 6.0 | 0.75 | 0.4226 | 0.0992 | 0.622 | 1.162 |
| 30.50 | 10.0 | 5.0 | 0.75 | 0.4701 | 0.1275 | 0.692 | 1.232 |
| 31.00 | 6.0 | 5.0 | 3.00 | 0.4442 | 0.0153 | 0.209 | **0.331** |
| 31.00 | 8.0 | 4.0 | 3.00 | 0.3517 | 0.0990 | 0.680 | **0.140** |

Shortest surviving buildable candidate is **30.50 µm (17.6% saving)**, not the
20.3% headline; 29.50 µm survives only if buildability is ignored (its best
candidate has s=0.25 < 0.75).

❌ **RETRACTED 2026-09-10 by the actual re-score — discard this screen
entirely, not just its caveats.** The measured `nt7 → nt13/tc30` phase shift
across all 668 points has median **+0.163°** and ranges **−94.9° to +3.9°** —
nothing resembling the uniform −0.54° assumed here. Concretely, the screen
predicted the 30.50 µm incumbent would fail at 2.045°; it actually passes at
**1.846°**, and the screen's "53 of 56 survive" bears no relation to the real
churn (46 kept, 10 lost, 24 gained, 70 total). The large shifts are at short
lengths (L 16.5-18.0) where phase sits near ±180°; `|phase − control|` does
not wrap, but control is 90° from the wrap boundary so no point can be both
near control and near wrapping, and those rows fail by 80°+ regardless.
**Kept as a record of a method that did not work: a single-geometry basis
shift does not extrapolate across this grid.**

⚠ (original caveat, superseded) The −0.54° shift was
measured at ONE geometry, and the taper basis interacts with the body's mode
content, so the shift certainly varies with `Lw`/`Lt`. Extrapolating it across
the grid is exactly the kind of assumption this study has been burned by. Its
real value is twofold: it refutes the negative-result expectation, and it shows
the current design point lands 0.045° over the line — close enough that its
fate depends on the per-geometry shift, which is precisely why the grid must
be RECOMPUTED rather than extrapolated. Note also that several strong
candidates sit at IL 0.42–0.47 dB against a 0.4809 dB threshold, so the
still-drifting IL axis can flip THEIR verdicts even though phase is settled —
argues for re-scoring at a high `TAPER_CELLS` (IL falls with cells: 0.3536 at
15 → 0.3220 at 30).

⚠ **First submission (13058925) was OOM-killed** on rank 4 at 8G/cpu: it asked
for `--ntasks-per-node=8` + `srun --mpi=pmi2`, but this sweep is mode-solve
dominated and MPB is rank-redundant, so 8 MPI ranks multiplied the memory
without buying speed. `run_mmi_2x2_waist_free.batch` — the known-good runner
for this same script, which already dispatches `--taper-convergence` — uses
`--ntasks=1 --cpus-per-task=8 --mem-per-cpu=12G` and plain `python`. Resource
pattern copied from it rather than invented.

#### (original plan text) 5b — `num_modes` through the narrowed section

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

#### ✅ FREE-SWEEP FORMULATION VALIDATED 2026-08-31 — 21× speedup, step 6 running

At the user's suggestion, step 6 uses the **analytic length sweep** rather
than re-solving every point. `solve_length_sweep()`'s docstring gives the
mechanism: a cell's cross-section depends on the taper **fraction**
`(i+0.5)/num_cells`, never on the region's physical length, and meow's
S-matrix factors into an interface part (modes only) and a propagation part
(the only place length enters). Its inner loop reuses one mode set and passes
`cell_lengths=` to `mw.compute_s_matrix()` — so an **arbitrary** cell-length
list works with the same modes.

Applied to the waisted body as **9 separate EME regions**
(`[wg_in, taper_in, stub, ramp, waist, ramp, stub, taper_out, wg_out]`, each
uniform segment getting ONE cell because it *is* uniform), **three of the four
axes become analytic**:

| axis | free? | why |
|---|---|---|
| stub `s` | ✅ | uniform at `mm_width` |
| waist `Lw` | ✅ | uniform at `Ww` |
| **ramp `Lt`** | ✅ | linear ramp — cells sample by *fraction*, not by `Lt` |
| `Ww` | ❌ | changes the widths sampled |

**One mode solve per (Ww, port) covers the entire (s, Lt, Lw) space** — 4
solves for a ~768-point grid, against ~21 min/point directly.

**Validated twice, deliberately.** A 21× speedup that is subtly wrong would
look like a *result*, not an error.

- **Single-region** (job 12624884): the library's own `solve_length_sweep`
  against direct solves — worst |ΔS| **1.168e-13**.
- **Multi-axis** (job 12627188): *my own* reimplementation, which varies three
  regions at once through hand-built cells — a different code path, so the
  first check does not cover it. Five cases, worst |ΔS| **1.109e-13**:

| s | Lt | Lw | \|ΔS\| bar | \|ΔS\| cross | |
|---|---|---|---|---|---|
| 2.75 | 4.00 | 8.00 | 1.109e-13 | 1.035e-13 | reference |
| 1.50 | 4.00 | 8.00 | 9.916e-14 | 9.324e-14 | stub only |
| 2.75 | 4.00 | 6.00 | 8.934e-14 | 1.011e-13 | waist only |
| 2.75 | 6.00 | 8.00 | 1.076e-13 | 8.261e-14 | **`Lt` alone** |
| 1.00 | 5.00 | 6.00 | 8.280e-14 | 9.650e-14 | all three |

The **`Lt`-alone** row is the load-bearing one: stub and waist being
length-free is intuitive, but the ramp claim is the one that would be easy to
get wrong and hard to notice. It holds at machine precision.

**Both of the user's hard constraints are enforced in the grid loop**, not
merely scored: candidates whose stub falls below the hybrid's folded 0.25 µm
are skipped (the FDTD sections must stay identical to the standard device's),
and candidates at or above the standard device's own `mm_length` are never
evaluated.

**Step 6 grid — job 12630600 TIMED OUT at 8 h. Re-running as 12674569
(2026-09-01).**

#### ⚠ 12630600: the Ww=2.52 half finished and was very nearly lost

The job died at the walltime **after** completing all of `Ww=2.52` and one row
of `Ww=2.24`. `free_sweep.json` was written only after the loop, so **nothing
reached disk**. The results survive only as filtered stdout in
`reports/12630600.out` (passers, plus every 40th row).

**✅ MEASURED — `Ww = 2.52` produces 34 passing candidates, 334/334 points, now
ON DISK** (job 12674569, 2026-09-01, re-run of the lost half). This is a
positive trend, not the pre-registered negative. Ranked by total length
(control 37.00 µm; thresholds IL ≤ 0.48091 dB, imb ≤ 0.1371,
|phase + 89.750| ≤ 2.0°):

| Ww | Lw | Lt | s | mm_length | total | IL dB | imb | phase° | IL margin |
|---|---|---|---|---|---|---|---|---|---|
| 2.52 | 10.0 | 6.0 | 1.00 | 24.000 | **33.00** | 0.4808 | 0.0065 | −90.749 | **9.27e-05** |
| 2.52 | 10.0 | 6.0 | 1.25 | 24.500 | 33.50 | **0.3197** | 0.0138 | −90.502 | 1.61e-01 |
| 2.52 | 10.0 | 5.0 | 2.25 | 24.500 | 33.50 | 0.3787 | 0.0176 | −90.434 | 1.02e-01 |
| 2.52 | 8.0 | 6.0 | 2.25 | 24.500 | 33.50 | 0.4600 | 0.0250 | −90.360 | 2.09e-02 |

**The re-run reproduces the timed-out job exactly**, which retires the question
of whether the stdout-recovered numbers were trustworthy: on the two rows that
12630600 printed, |ΔIL| ≤ 1.4e-05 dB, |Δimb| ≤ 3.0e-05, |Δphase| ≤ 4.3e-04° —
i.e. differences only in the 4th decimal place the old log printed to. The
whole Ww=2.52 passer sequence also appeared in the same order at the same row
counts.

⚠ **The shortest passer is grid-edge marginal on loss and must not be quoted as
the design point.** Its margin to the pre-registered threshold is
**9.27e-05 dB**. The `L = 24.5, Lw = 10, Lt = 6` point one grid step away costs
0.5 µm and clears the same threshold by **1.61e-01 dB — 1700× the margin**.
A criterion is a cliff, and a point sitting 9e-05 dB inside it is not a
measurement of "passes", it is a measurement of "the grid step landed there".
The `--confirm` path below re-solves both directly, with neighbours.

#### ✅ Two fixes landed 2026-09-01 before the re-run

1. **Rows now stream to disk.** Every candidate is appended to
   `free_sweep_rows_<res,nwl,ctrans>_<tag>.jsonl` as it is produced, and a run
   skips any `(Ww, Lw, Lt, L)` already present in any shard. A timeout now
   costs only the unfinished points, and the grid can be split across jobs
   (`MMI2X2_FS_WW` / `MMI2X2_FS_LW` / `MMI2X2_FS_TAG`). `--summarize`
   aggregates the shards with no solving. The settings fingerprint is in the
   **filename** so a low-res smoke run cannot be read back as done grid rows.

2. **The "free" evaluation was not actually free — now it is. ✅ MEASURED
   EXACT, job 12674313.** 69 s/point was the whole reason 384 points did not
   fit 8 h. `meow.compute_s_matrix` accepts `interfaces_fn`, and
   `compute_interface_s_matrices(modes)` is a **pure function of the modes** —
   which are identical at every point of this grid, by the same argument the
   free sweep itself rests on. It was being recomputed 384 times.

   | piece of one `compute_s_matrix` | time |
   |---|---|
   | interfaces (modes only) | 3.07 s |
   | **propagations (the only length-dependent part)** | **0.01 s** |
   | cascade/sax remainder | 4.65 s |

   Memoising `interfaces_fn` per `(port, wavelength)`: **0.34 s vs 3.44 s,
   9.9×, worst |ΔS| = 0.000e+00 — exactly zero**, at three off-reference
   `(s, Lt, Lw)` points. Not an approximation; the propagation factor is 0.1%
   of the cost and is still recomputed in full.

   **The speedup self-checks inside the grid job.** Per this repo's standard a
   10× that is subtly wrong would read as a result, so
   `selfcheck_iface_cache()` runs cached-vs-uncached at two points **in the
   same process as the grid** and `SystemExit`s if |ΔS| > 1e-13 — it certifies
   that job's cache, not a previous job's probe. `MMI2X2_FS_IFCACHE=0` falls
   back.

   Expected re-run cost: **~2 h for the full grid** (both bases 36 min + ~740
   points × 7.5 s), against 15 h for the old path.

#### ✅ Smoke PASSED 2026-09-01 (job 12674549) — the machinery, not the physics

res 15 / 3 wl / 3 lengths at `Lw=8, Lt=4`, its own shard fingerprint. All three
mechanisms confirmed live: the shard file
`free_sweep_rows_res15_nwl3_ct8_smoke.jsonl` was written **row by row**, the
resume skip-set loaded, and the in-process cache self-check reported
**|ΔS| = 0.000e+00**.

Its "NO candidate passes" line is **not** a result and does **not** trigger the
stop rule: 3 points at res 15 on an `(Lw=8, Lt=4)` line that already failed at
res 40. Read it as "the script ran", nothing more.

⚠ **Two script-hygiene bugs the smoke exposed, both fixed after 12674569 had
already started** (so they take effect from the next run, not this one):

1. The self-check's "cached" timing was meaningless — that pass is what
   *populates* the cache, so it pays the full interface cost exactly like the
   uncached one (smoke printed `cached 10.3s vs uncached 7.0s`, which reads as
   a slowdown). The speedup only exists warm, which is the regime the grid loop
   actually runs in. A warm third pass is now timed and reported.
2. The end-of-run summary was a fixed `free_sweep.json`, so the **smoke
   overwrote the grid's own summary file**. It is now
   `free_sweep_<fingerprint>_<tag>.json`. The smoke-written copy was deleted;
   12674569 still writes the old fixed name.

#### ✅ STEP 6 GRID COMPLETE 2026-09-01 — job 12674569, `COMPLETED` in 02:14:01

668/668 pre-registered points, **56 passers**, 0 faults, exit 0. Both bases
self-checked at `|ΔS| = 0.000e+00` **independently** (1124 s and 1143 s) — the
second Ww's cache is certified against its own mode set, not inherited.
Measured cost 2267 s of mode solves for the whole grid; the pre-memoisation
path would have needed ~15 h.

**`Ww = 2.24` — never simulated before this job — is the better profile**, and
it beats `Ww = 2.52` on length at every comparable point. Shortest passers of
the whole grid, with each component as a fraction of its pre-registered budget:

| Ww | Lw | Lt | s | L | total | IL% | imb% | Δphase(mean)° | ph% | worst |
|---|---|---|---|---|---|---|---|---|---|---|
| 2.24 | 8.0 | 6.0 | **0.25** | 20.500 | **29.50** | 81% | 62% | 1.696 | 85% | 85% |
| 2.24 | 8.0 | 6.0 | 0.50 | 21.000 | 30.00 | 78% | 69% | 1.084 | 54% | 78% |
| 2.24 | 6.0 | 6.0 | 1.75 | 21.500 | 30.50 | 81% | **9%** | 1.505 | 75% | 81% |
| 2.24 | 8.0 | 5.0 | 1.75 | 21.500 | 30.50 | 74% | 73% | 1.574 | 79% | 79% |

**29.50 µm vs the 37.00 µm control is −7.50 µm, 20.3%.** The stop rule does
**not** trigger; the deliverable is a design point, not a negative.

> ❌ **RETRACTED 2026-09-12.** This 20.3 % figure belongs to a point that fails
> on both settled taper bases and on full-device FDTD. The delivered saving is
> **16.2 % (31.00 µm)**. See the STATUS block at the top of this file.

⚠ **No point in this grid has broad headroom** — every one of the top six uses
78–88% of *some* budget. And the two most attractive points are both
**constraint-boundary, not interior**: the 29.50 µm winner has `s = 0.25`
*exactly* on the folded-hybrid stub floor (it cannot shorten along that axis at
all) and spends 85% of the phase budget. The `Lw=6/Lt=6/L=21.5` point at
30.50 µm costs 1 µm and is interior on the stub with **9% of the imbalance
budget against 62%**. That trade is the real decision, and it is the user's.

#### ⚠ FINDING: the phase criterion is scored MORE LENIENTLY than pre-registered

`criterion.json`'s own words: `"phase": "max |bar/cross phase - control phase|
<= 2.0 deg"`. What the code does (`metrics()`/`score()`): `phase_deg =
np.mean(dphi)` over the band, compared to `control_phase_deg`, itself a
`np.mean`. **Excess loss and imbalance are both max-over-band — phase alone is
a mean.** Equal-and-opposite excursions at the two band edges cancel, so a
candidate can swing several degrees across the band and still score as
sub-2°.

This is not academic here: the top candidates spend 75–85% of the phase budget
**on the lenient reading**.

**The criterion was NOT redefined after the sweep.** `score()` still applies
exactly what it applied before, and the 56-passer set above stands as scored.
Two reported-only fields were added (`phase_max_dev_deg`, `phase_spread_deg`).

##### ✅ MEASURED, job 12687020 — 4 of 56 passers fail the strict reading, and the LENGTH RESULT SURVIVES

| Ww | Lw | Lt | L | total | mean-dev° (scored) | **MAX-dev°** | spread° |
|---|---|---|---|---|---|---|---|
| 2.24 | 8.0 | 6.0 | 20.500 | 29.50 | 1.696 | **2.244** | 1.419 |
| 2.24 | 8.0 | 4.0 | 21.500 | 30.50 | **0.883** | **2.076** | 2.198 |
| 2.24 | 8.0 | 4.0 | 22.500 | 31.50 | 1.954 | **2.239** | 0.889 |
| 2.52 | 10.0 | 3.0 | 26.000 | 35.00 | 1.988 | **2.014** | 0.089 |

The `Lw=8, Lt=4, 30.50 µm` row is the case that justifies having checked at
all: its **mean** deviation is 0.883°, a comfortable 44% of budget, while the
band actually spreads 2.198° and peaks 2.076° — past the limit. Averaging
converted a failing point into a comfortable-looking one.

**The 20.3% length saving is unaffected by which reading governs.**

> ❌ **RETRACTED 2026-09-12 — and this claim is the one that aged worst.** The
> reading did not merely pick a different 29.50 µm point; on a SETTLED taper
> basis both 29.50 µm candidates fail outright, so the saving was never
> reading-independent — it was basis-dependent, and the basis was wrong. The
> delivered design point (`Lw6 Lt5 L22.0`, 31.00 µm, 16.2 %) is genuinely
> reading-independent: 0.267° mean and 0.447° max, both far inside 2.0°.

> **Shortest passer surviving BOTH readings: `Ww = 2.24, Lw = 10.0, Lt = 5.0,
> L = 20.500` → 29.50 µm total** — IL 0.4255 dB, imb 0.1095, max-phase-dev
> **1.516°**. Same 29.50 µm as the mean-reading winner: the strict reading
> changes *which* point wins, not what it achieves.

It is also the better point on its own terms — mean-dev 0.561° against the
displaced winner's 1.696°, so it is not surviving by an averaging accident. It
costs more of the loss and imbalance budgets (88%/80% vs 81%/62%) and, like the
displaced winner, sits on the stub floor at `s = 0.25`.

**✅ DECIDED by the user, 2026-09-01: the MEAN reading governs.** That is what
`score()` has always applied, so the 56-passer set and its ranking stand
unchanged, and the design point is the shortest passer under it:

> ### ❌ RETRACTED 2026-09-12 — THIS POINT FAILS. The design point is now `Lw = 6.0, Lt = 5.0, L = 22.0` (31.00 µm).
> Retained only as the record of what step 6 originally concluded. Do not quote
> the 20.3 % headline below; it does not survive.
>
> **Why it fails, by three independent routes:**
> * **Re-scored grid** (jobs 13061906 / 13065614) at BOTH settled taper bases:
>   mean-dev **2.227°** at nt13/tc30 and **2.160°** at nt21/tc30, against a
>   2.0° budget. The two bases agree on 660/668 verdicts, so this is not a
>   basis artifact.
> * **Full-device FDTD** (check 4): mean-dev **2.099°**, MAX-dev 2.323°.
> * **Step 5b's 2-D basis map**: FAILS at every one of 13 settled
>   `(N_TAPER, TAPER_CELLS)` points, 2.05–2.30°. Its original 1.696° at
>   `N_TAPER=7` is the sole outlier in the whole map — basis luck, not margin.
>
> It was also never hybrid-buildable: `s = 0.25 µm` is below the
> `s ≥ body_stub + dpml = 0.75 µm` floor.
>
> **The replacement is measured, not inferred** — see "CHECK 4 CONFIRMS A NEW
> DESIGN POINT": `Ww2.24 Lw6.0 Lt5.0 L22.0`, 31.00 µm (16.2 % shorter),
> IL 0.2010 dB · imbalance 0.0124 · mean-dev **0.267°** · MAX-dev **0.447°**,
> i.e. it passes on BOTH phase readings and so does not depend on the
> mean-vs-max ambiguity that decided this one.
>
> ### 📌 DESIGN POINT (as scored at N_TAPER=7) — `Ww = 2.24, Lw = 8.0, Lt = 6.0, L = 20.5` (stub 0.25)
> **29.50 µm total vs the 37.00 µm straight control — 7.50 µm shorter, 20.3%.**
> IL 0.3910 dB · imbalance 0.0852 · mean phase dev 1.696° · direct-solve
> confirmed (job 12687092, `dIL = −1.05e-14` against the analytic grid).

Recorded plainly rather than buried: **this point's max-over-band phase
deviation is 2.244°**, which would fail the 2.0° tolerance under the max
reading of the same pre-registered text (job 12687020). The user chose the mean
reading with that number in hand. The strict-reading survivor
(`Lw = 10.0, Lt = 5.0, L = 20.5`) reaches the **same 29.50 µm** and remains
available at no length cost if the reading is ever revisited — see the
comparison table below, which stands as measured either way.

##### ✅ STEP 6 CLOSED 2026-09-01 — 6/6 direct confirms pass (jobs 12687092, 12692011)

Every candidate re-solved with real per-length mode solves. **All six reproduce
the analytic grid at 1e-14 and all six pass on a direct solve.**

| Ww | Lw | Lt | s | L | total | IL dB | imb | MAX-dev° | spread° | strict | \|dS\| vs analytic |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2.24 | 8.0 | 6.0 | 0.25 | 20.5 | 29.50 | 0.3910 | 0.0852 | **2.244** | 1.419 | **FAIL** | dIL −1.05e-14 |
| 2.24 | 8.0 | 6.0 | 0.50 | 21.0 | 30.00 | 0.3728 | 0.0944 | 1.869 | 1.669 | ok | dIL −1.89e-14 |
| 2.24 | **10.0** | **5.0** | 0.25 | 20.5 | **29.50** | 0.4255 | 0.1095 | 1.516 | 1.732 | **ok** | dIL +1.28e-14 |
| 2.24 | 10.0 | 5.0 | 0.50 | 21.0 | 30.00 | 0.4486 | 0.1194 | 0.814 | 0.834 | ok | dIL −1.23e-14 |
| 2.24 | **6.0** | **6.0** | 1.75 | 21.5 | **30.50** | 0.3880 | **0.0127** | 1.640 | **0.237** | **ok** | dIL +2.43e-14 |
| 2.24 | 6.0 | 6.0 | 2.00 | 22.0 | 31.00 | 0.3591 | 0.0348 | 1.412 | 0.264 | ok | dIL −1.30e-14 |

**The grid-edge test comes back clean for both surviving candidates.** Each
passes at its own length *and* at its L-neighbour, so both are plateaux rather
than spikes. The one strict failure is the only point whose neighbour behaves
differently — further evidence it is an isolated band excursion, not a region.

**Two defensible design points, and the choice is a real engineering trade:**

| | **29.50 µm** `Lw=10, Lt=5, L=20.5` | **30.50 µm** `Lw=6, Lt=6, L=21.5` |
|---|---|---|
| saving vs 37.00 µm control | **20.3%** | 17.6% |
| excess loss | 0.4255 dB (88% of budget) | 0.3880 dB (81%) |
| imbalance | 0.1095 (**80%**) | 0.0127 (**9%**) |
| max phase dev | 1.516° | 1.640° |
| band spread | 1.732° | **0.237°** |
| stub | **0.25 µm — ON the hybrid floor** | 1.75 µm — interior |
| passes both phase readings | ✅ | ✅ |

The extra 1 µm buys an order of magnitude of imbalance headroom, an essentially
flat band, and clearance from the fabrication-relevant stub constraint. **The
shortest device is not automatically the right one here**, and per this doc's
own ranking rule (shortest among passers) the 29.50 µm point wins on the
letter — which is exactly why the trade is being written down rather than
silently resolved.

⚠ **Filename bug, found and fixed.** Both confirm jobs wrote the same untagged
`confirm_<fingerprint>.json`, so 12692011 clobbered 12687092's rows. The path
is now tag-keyed, and the lost rows were recovered from the SLURM log into
`confirm_res40_nwl11_ct8_job12687092.json`. Same failure mode as the timeout
that started this session: results with nowhere durable to land.

##### ✅ The analytic fast path is vindicated on a real solve

The first direct confirm (`Ww=2.24, Lw=8, Lt=6, L=20.5`, job 12687092) returned
**`dIL = −1.05e-14, dphi = +4.26e-14`** against the analytic row — machine
precision on a genuine per-length solve, not a self-comparison. Every number in
this section rests on a path that is now checked three independent ways:
`solve_length_sweep` vs direct (1e-13, 12624884), the multi-axis
reimplementation vs direct (1e-13, 12627188), the interface cache vs uncached
(0.000e+00, 12674313 + in-process each run), and now the finished grid vs
direct (1e-14).

#### ⚠ FINDING: the stop rule could fire off a PARTIAL grid — now guarded

Testing `--summarize` against the 3-row smoke shard printed
**`** PRE-REGISTERED STOP RULE TRIGGERED **`**. Three low-res smoke points on
one `(Lw, Lt)` line were enough to produce the sentence that this document
treats as a publishable mechanism result. Nothing had gone wrong physically —
the summary block simply never asked whether the grid it was summarising was
finished.

Fixed: the summary now reconstructs the **literal pre-registered grid** (Ww
{2.52, 2.24} × Lw {4,6,8,10} × Lt {3,4,5,6} × 24 lengths in [16.0, 27.5], minus
the stub and length filters) = **668 legal points**, 334 per Ww, and reports
coverage. Below full coverage it prints "stop rule NOT evaluated: the grid is
PARTIAL" instead. `grid_complete` / `n_missing` are recorded in the summary
JSON.

The grid definition used for this check is **literal, not read from
`MMI2X2_FS_WW`/`_LW`/`_LT`** — a single-Ww shard job would otherwise declare
the whole grid COMPLETE with the other Ww never simulated. That is the same
class of bug as the stop rule itself: a narrower run quietly certifying a
broader claim.

**Revised cost estimate from the real point count:** 334 points/Ww × ~7.5 s
warm + ~18 min basis + ~5 min in-process self-check ≈ **1 h per Ww**, ~2 h for
the grid.

#### The direct confirm step 6 owes is now wired, not just noted

`--confirm` (batch: `VALIDATE=confirm`) re-solves the design point **and its
two L-neighbours** with real per-length mode solves. This is **not** a check on
the analytic path — that is settled at 1e-13 (12624884) and its interface cache
at 0.000e+00 (12674313). It is the check the script's own closing warning
demands: that the quoted point is a real solve and is **not perched on a grid
edge**, which is exactly the risk the 1e-4 dB margin above flags. If a
neighbour of a passing point fails, the script says so and refuses to endorse
it. `--summarize` (`VALIDATE=summarize`) aggregates the shards with no solving,
so the standing best can be read between jobs.

**Job in flight (2026-09-01):** `12674569` — full grid, both Ww, res 40, 11 wl,
`-t 08:00:00`, `gts-sr108-gedc_main`, started 09:15. **Not landed. Every
Ww=2.52 number in the table above is recovered stdout from the timed-out
12630600 and is superseded by this job's own shard file once it lands.**

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

#### ✅ STEP 7 DONE 2026-09-01 — jobs 12696051 (bar) + 12696052 (cross), res 40

24 PNGs under `example_output/verification/mmi_2x2_waist_fields/fields/`
(straight/waisted × bar/cross × 6 components), band centre 1.5475 µm, both
cases from one config so the waist is the only difference. Geometry verified in
the log: waisted body `in-face 2.800 / mid 2.240 / out-face 2.800`, straight
uniform at 2.800.

**What the comparison shows** (`eme_Ex_*`, meow's Ex = the device's in-plane
transverse field, the dominant TE component):

1. **Nothing is shed at the ramps.** The cladding beyond the body edge stays
   dark through both transition bands (z ≈ 10.7 and 18.8 µm, exactly where the
   `stub 0.25 / Lt 6 / Lw 8` geometry puts them). Transition radiation — the
   hypothesised 1×4 failure mechanism — is **visibly absent here**, which is
   the qualitative counterpart of the measured 0.3910 dB.
2. **The waist is visible as a boundary pinch**, ±1.4 µm → ±1.12 µm, over
   exactly the predicted span.
3. **Self-imaging survives the two ramps intact.** The waisted device completes
   the same cross-over and reforms two clean lobes at ±0.75 µm at its output,
   7.5 µm earlier than the control.

Point 3 is the mechanism in a picture, and it is not truncation: **narrowing
the body raises Δβ₀₁, which contracts L_π, so the self-image forms sooner.**
First quantitative support the same day — the waist cross-section's
`L_π = 12.57 µm` against the full-width body's **19.40 µm** (a 35% contraction;
res-15 smoke, job 12696148, refined by the check-3 run below).

#### (superseded) IN PROGRESS note — `plot_mmi2x2_waist_eme_fields.py`, smoke 12695975

Straight control vs the design point above, one config, only the waist
differing. Design point overridable via `MMI2X2_WF_WW/LW/LT/L` (e.g. to render
the strict-reading survivor instead).

**This does NOT copy the 1×4's template, and the reason is structural.**
`plot_mmi1x4_waist_profile_eme_fields.py` bypasses `_build_cells()` and
hand-rolls the body segmentation, because the 1×4 sweep declares the body as
**one** region that `_build_cells()` would divide uniformly, unable to place
boundaries at the stub/transition/waist edges. Its own docstring warns the two
files' segment logic must be kept in sync by hand.

The 2×2 free sweep already avoided that problem by declaring the body as
**nine** regions, each uniform or a single linear ramp — for which an even
division per region is exactly right. So this script goes through the normal
solver path and **imports the sweep's own layout shape rather than copying
it**: the plotted geometry is the swept geometry by construction and cannot
drift out of sync. That removes the 1×4 script's standing maintenance hazard
rather than inheriting it.

**One-cell uniform regions render correctly** — `plot_eme_fields()` builds a
500-point `z` grid over the whole device and `mw.propagate_modes` evaluates
analytically *within* a cell, so a uniform region needs one cell for the
S-matrix and still renders at full spatial resolution. No cell-count inflation
is needed for the plot.

#### ⚠ `plot_hybrid_body_fields()` does NOT apply to this device — use the multi-cell path

A shared helper landed 2026-08-31 (`beef.plot_hybrid_body_fields`,
`beef/plot/hybrid_body.py`) that reconstructs a hybrid cascade's **EME body
section** field from one mode solve. **It cannot be pointed at a waisted body,
and the reason is physics rather than a missing feature.**

It assumes a **uniform** body cross-section, and two load-bearing consequences
follow: one cross-section's basis serves the whole section, and with
`calibrate_body_neff=True` the body transfer matrix is a **pure diagonal**
`exp(2j·π·n_eff/λ·L)`, so the field is just `Σ_k c_k φ_k(y) e^{iβ_k x}`. That
sum *is* the reconstruction.

A waisted body breaks both. It goes through `solve_section()`'s real
**multi-cell** branch — one cross-section per cell, real interface overlaps,
**non-diagonal** transfer matrix. `beef` already enforces this:
`assert_uniform_cross_section()` is called from both
`EMESolver.solve_section()`'s override path and `HybridSolver._calibrated_neff()`
(`hybrid.py:1534`), so the result is a **hard assertion, not a subtly wrong
plot** — the right failure mode.

**What to use instead** — and the library already has it:
`plot_eme_fields()` + `mw.propagate_modes()` handle multi-cell propagation
natively. `plot_mmi1x4_waist_profile_eme_fields.py` (written for the 1×4 waist)
is the working template: build cells for the body span with the **real
`num_cells`**, solve modes per cell, hand both lists to `plot_eme_fields()`.

**Two traps worth stealing from the helper's development**, both real bugs:

- **`plot_eme_fields()` labels in meow's frame, not the device's** —
  `device Ex = meow Ez`, `device Ey = meow Ex`, `device Ez = meow Ey`. The
  panel titled `|Ex|` is the device's **in-plane transverse** field and its
  `x [µm]` axis is device **y**. This affects how the 1×4 waist field plots
  already in this repo should be read (that script's own comment does record
  it — it is consistent, not wrong).
- **Any modal projection must be a 2-D area integral, never a 1-D cut through
  the core.** Measured on the 1×2: a 1-D slice gave captured power **1.5857**
  (impossible) and put **0.5876** on a parity-forbidden mode; the 2-D `dy·dz`
  integral gave 0.9970 and 1e-11. **The parity zeros survived both**, so a
  parity check does not catch this — assert `captured_power <= 1`.

**The standard (uniform-body) 2×2 is a different matter and the helper does
apply there** — now wired into `examples/mmi_2x2/python/mmi_2x2_hybrid.py`
and its notebook, with the off-axis launch (`input_ys = [-0.75, +0.75]`, not
0) and a `port_phases()` readout, since this device's own criterion includes
bar/cross quadrature.

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

#### ✅ Check 3 (resolution, EME side) — CLOSED, jobs 12696267 + 12708666

**Why this was gating rather than routine.** The whole step-6 grid ran at
**res 40**, but this device's straight sibling measured **Δβ₀₁ needing res 100**
for 0.05% (job 11590785, `resolved/mmi_2x2_validation.md` gap 2). The design
length is set by phase accumulation — `L_π = π/Δβ` — so an unconverged Δβ moves
the optimum, and ~1% is ~0.2 µm against a 0.5 µm grid step. The pass/fail
*comparison* is protected (the control was also measured at res 40); the
absolute 20.5 µm is not.

`--res-convergence` measures both at once, and states **which cross-section**
it is asking about — the **uniform waist**, per step 8's own note that a
waisted body has no single one. Transitions stay covered by step 5a's C_TRANS
convergence.

**MEASURED so far (2 of 4 blocks):**

| res | n_eff_0 | n_eff_1 | Δβ₀₁ /µm | L_π µm | best-IL L | IL@design | maxdev° |
|---|---|---|---|---|---|---|---|
| 40 | 2.832691 | 2.770181 | 0.253805 | 12.3780 | 20.750 | 0.3910 | 2.244 |
| 60 | 2.831598 | 2.769116 | 0.253691 | 12.3836 | 20.750 | 0.3745 | 2.340 |

**L_π moves 0.045% from res 40 → 60** — the straight body needed res 100 to
reach 0.05% on the same quantity. The waist's two lowest modes are better
separated (narrower guide), so **the quantity the design length depends on is
markedly less resolution-sensitive here than the gating concern assumed.** The
best-IL length is **stable at 20.750 µm** across both.

The design point stays PASS at both. Loss shifts 0.3910 → 0.3745 dB — 0.017 dB,
and in the *favourable* direction; imbalance and phase are unchanged to 3-4
decimals. Max-dev creeps 2.244 → 2.340°, so it would fail the strict phase
reading at **both** resolutions — consistent, not a new defect, and moot under
the mean reading the user adopted.

**Res 80 landed and Δβ₀₁ is CONVERGED — it oscillates, it does not trend:**
12.3780 → 12.3836 → 12.3730, a ±0.09% band. That is the signature of a
quantity sitting on solver noise rather than one still moving. **The best-IL
length is 20.750 µm at every resolution tried**, and the design point stays
PASS throughout with its loss drifting *favourably* (0.3910 → 0.3745 → 0.3601
dB) — so the res-40 grid was mildly pessimistic about the design point, never
optimistic.

⚠ **OOM, not a physics problem — job 12696267 was killed at res 100.**
`OUT_OF_MEMORY` at 4:41 with MaxRSS **100.7 GB** against a 96 GB request. The
walltime was never the binding constraint (I spent two hours tracking the wrong
resource). Cause is structural: the loop held the **previous** resolution's
basis alive through the assignment, so peak memory spanned two resolutions at
once, and mode-field arrays scale as res² (res 100 is 6.25× res 40 on its own).
Fixed — the old basis is dropped and `gc.collect()`ed before the next is built.

⚠ **And the JSON wrote only after the loop**, so res 40/60/80 produced no file
— the same defect that cost job 12630600 its Ww=2.52 half this morning, in a
code path I had added *after* fixing it in the grid path. The numbers survived
only because each block was transcribed here as it landed. Both fixed: the mode
now writes a partial JSON with a `complete` flag after **every** resolution, and
the three blocks were recovered from the SLURM log into
`res_convergence_res40-60-80.json` with provenance.

**✅ CHECK 3 CLOSED 2026-09-01 — CONVERGED (res 100 via job 12708666,
`COMPLETED`, basis 10952 s, MaxRSS 86.0 GB of 128 GB).**

| res | Δβ₀₁ /µm | L_π µm | ΔL_π vs res 100 | best-IL L | IL@design | |
|---|---|---|---|---|---|---|
| 40 | 0.253805 | 12.3780 | +0.024% | 20.750 | 0.3910 | PASS |
| 60 | 0.253691 | 12.3836 | +0.069% | 20.750 | 0.3745 | PASS |
| 80 | 0.253907 | 12.3730 | −0.016% | 20.750 | 0.3601 | PASS |
| 100 | 0.253867 | 12.3750 | — | 20.750 | 0.3627 | PASS |

**L_π spans 0.086% across a 2.5× resolution range and oscillates rather than
trends** — converged, and converged already at res 40. **The res-40 grid that
produced the design point was sound.** Against this device's straight sibling
needing **res 100** for 0.05% on the same quantity, the waist converges
dramatically faster: the narrower guide separates the two lowest modes better,
so the quantity the design length depends on is far less mesh-sensitive.

**The best-IL length is 20.750 µm at every resolution tried** — the optimum is
not resolution-dependent, which was the specific risk that made this check
gating. The design point holds PASS throughout, with loss landing at 0.3627 dB
at res 100 against 0.3910 dB at res 40: the grid was mildly **pessimistic**
about it, never optimistic.

Sizing notes for whoever runs this next: basis wall time scales as **res^2.44**
(1235 / 3103 / 6219 / 10952 s at res 40/60/80/100, both ports, 11 wl, 51
cells), and a single res-100 basis peaks at **86 GB** — 128 GB is the right
request, 96 GB is not, and the first attempt's 256 GB only bought queue delay.

📌 **The design point is pinned to a constraint, not to an interior optimum.**
At `Lw=8, Lt=6` the stub floor `s ≥ 0.25` forbids `L < 20.5` outright, so the
window is one-sided and 20.5 µm is the shortest length the geometry permits —
not a resolved minimum. Note also that **20.750 µm is the loss optimum**, beating
the design point by 0.042 dB; 20.5 wins only on the ranking rule (shortest among
passers). If a future run wants to go shorter, a different `(Lw, Lt)` is
required — this one cannot express it.

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

#### 🏃 Check 4 (full-device FDTD reference) — job 12736252, res 40, 51 wl

**Why this became the load-bearing check rather than a closing formality.**
`calibrate_body_neff=True` is **hard-refused for a waisted body** —
`assert_uniform_cross_section()` fires from `HybridSolver._calibrated_neff()`
(`hybrid.py:1532`): *"a calibrated body propagates on one per-mode beta, which
only describes a uniform cross-section"*. That flag is precisely what closed
check 1 elsewhere in this repo (MPB's mode index vs meep's time-stepped phase,
~2.3%/mode; on the coupler it moved L₅₀ from 5.37% off to **0.11%**). So a
waisted hybrid necessarily carries an **uncorrected** phase error in the
quantity this device's criterion scores most tightly. A full-device FDTD run
depends on neither the cascade nor the calibration, so it is the only available
way to confirm the design point against something that is not EME.

Config and rank layout are inherited **verbatim** from the straight device's
own check-4 reference (`mmi_2x2_fdtd.py`, 2 nodes × 24 = 48 ranks, job
11585844_1) so this is a comparable measurement, and the 51-wavelength grid
**contains the criterion's 11-point grid exactly** (every 5th point), so the
comparison is on identical wavelengths rather than interpolated ones. Both
cases are scored with the same metric functions the EME sweep used.

##### ✅ The STRAIGHT CONTROL lands, and the engines agree better than feared

| quantity | FDTD res 40 | EME (criterion control) | delta |
|---|---|---|---|
| excess loss | **0.1398 dB** | 0.3809 dB | **−0.241** |
| imbalance | 0.0655 | 0.0686 | −0.003 |
| phase (mean) | −90.461° | −89.750° | **−0.711°** |
| phase MAX-dev | 1.103° | — | — |
| band spread | 0.590° | — | — |

**Engine agreement on phase is 0.711°** — inside the 2.0° budget. The
`calibrate_body_neff` refusal is a real structural limitation but a small
numerical one on this device, which is the reassuring answer and could not have
been known without measuring.

⚠ **The one real disagreement is LOSS, and its direction matters: EME
over-estimates by 0.24 dB.** 0.1398 dB (FDTD) vs 0.3809 dB (EME) on identical
geometry. That is the expected sign for a truncated 15-mode EME basis shedding
power a brute-force solve retains. **If the offset carries to the waisted case,
the design point's true loss is ≈0.13 dB rather than 0.39 dB** — i.e. it has
far more loss margin than the criterion credits, and the 9.27e-05 dB
"margin" that made the *other* 29.50 µm candidate look knife-edge was an EME
artifact. Do not assume it carries; the waisted row will say.

⚠ **A res-15 smoke of this same script reported the straight phase gap as
−3.70°, and I flagged it as budget-breaking.** It was a coarse-mesh artifact;
at res 40 it is −0.711°. Recorded because the smoke's numbers are in the log
and should not be mistaken for a result.

##### 📌 Measured for free: the step-9 cache-hit prediction is CONFIRMED

`[straight] solved in 4s` — a **cache hit**, not a solve. Holding
`mm_length = 28.0` made the run land on the standard 2×2's own completed
reference entry, exactly as step 9's cost note predicted but explicitly left
unverified ("verify the hit before relying on it"). It is now verified, and the
straight control cost nothing.

⚠ Minor: the script's per-case prints appear **48 times**, once per MPI rank.
Cosmetic, but prints should be rank-0 guarded.

##### 📌 MEASURED: `decay_tol=1e-11` is ~7 decades stricter than this device needs

The repo inherits `decay_tol=1e-11` everywhere with **no recorded
justification** (that is a standing note in `decay-tol-has-a-hard-floor`).
Job 12736252's snapshots settle it for the 2×2, at no extra cost, because DFT
monitors accumulate: sample *t* equals what a run stopped at *t* would have
returned.

Band-centre transmission against the energy-ratio axis (waisted, res 40):

| t | energy ratio | T_bar | T_cross | ΔT_bar per 20 t |
|---|---|---|---|---|
| 360 | 1.97e-03 | 0.506004 | 0.428740 | 5.04e-02 |
| 380 | 4.23e-04 | 0.519244 | 0.436904 | 1.32e-02 |
| 400 | 1.50e-04 | 0.521083 | 0.439559 | 1.84e-03 |
| 420 | 4.75e-05 | 0.521367 | 0.439845 | **2.84e-04** |

**Transmission stops moving by ~1e-4 at an energy ratio of ~1e-4**, while
`decay_tol=1e-11` (dft) targets an energy ratio near **1e-13** — about
**seven decades and several hours of stepping past the point where the answer
is fixed.** The decay also eased monotonically after t≈360 (0.89 → 0.39
decades/20 t) with **no growth phase**, i.e. the 1×2's honest-slow-tail
pattern, not the 4×4's floor-then-diverge instability.

**Consequence, acted on.** 12736252 was cancelled at t≈510 and re-run as
**`12743476` at `decay_tol=1e-5`** (~1e-7 on the energy axis, two decades past
settling and already demonstrably reached on this device). The snapshots carry
only |T|, not complex amplitudes, so a walltime kill would have cost the
**phase** measurement — the one criterion component under scrutiny — which is
why waiting for a run projected to overrun its wall was the worse option.

⚠ **The two cases therefore carry different tolerances**: straight control at
1e-11 (a cache hit on the standard device's own reference), waisted at 1e-5.
Normally a mismatch to avoid; here the *control* is the stricter one and both
are converged in transmission far earlier, so the direction is safe. Recorded
rather than glossed.

---

#### 📌 TODO (user, 2026-09-06): per-section β calibration, AFTER the uncalibrated run lands

**Today the waisted hybrid runs with `calibrate_body_neff=False`, and that is
forced, not chosen** — `assert_uniform_cross_section()` refuses a calibrated
body whose cross-section varies (`hybrid.py:1532`). So this cascade carries an
**uncorrected body β**, i.e. the ~2.3%/mode MPB-vs-meep phase gap that
`calibrate_body_neff` exists to remove, in the one quantity this device's
criterion scores most tightly.

**Splitting the body changes what is possible.** With the body decomposed as

```
stub (uniform) | ramp | waist (uniform) | ramp | stub (uniform)
```

the **stub and waist sections are uniform**, so they are exactly the kind of
region calibration is valid for — and at the design point the waist alone is
6.0 µm of the 21.0 µm body, with the stubs adding 3.0 µm more. Only the two
ramps genuinely cannot be calibrated.

⚠ **The flag cannot express that today.** `calibrate_body_neff` is a single
global `bool` (`hybrid.py:400`, `443`) applied to **every** EME section that
supports `neff_override`, so turning it on would attempt the ramps too and trip
the assertion. Making it per-section is a small, additive library change —
accept a bool *or* a collection of section indices, defaulting to today's
behaviour.

**Order matters and is deliberate**: run the **uncalibrated** cascade first, so
there is a measured baseline, then calibrate the uniform regions. Doing it the
other way round would leave no way to attribute the difference.

📌 **The figure of merit is the gap to the FULL-DEVICE FDTD, not internal
self-consistency** (user, 2026-09-06). The reference already exists — job
12797278 at this exact geometry:

| quantity | full-device FDTD |
|---|---|
| excess loss | **0.2178 dB** |
| imbalance | **0.0168** |
| phase (mean) | **−91.661°** |

So the test is: does **\|hybrid − FDTD\|** shrink, per quantity, when the
uniform stub and waist sections are calibrated? Calibration exists precisely to
remove the MPB-vs-meep β discrepancy between the EME body and the FDTD sections,
so if it is doing its job the cascade should move **toward** the brute-force
answer. If the gap does not shrink, calibration is not buying what it is
assumed to buy on this device — which would itself be worth knowing, and is the
kind of claim the straight 2×2 could never test because its body is uniform
end-to-end and always calibratable.

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

#### ❌ HARD CONSTRAINT, MEASURED 2026-09-03: a hybrid-buildable waist needs `s ≥ 0.75 µm`

Found by trying to build the missing hybrid leg (the straight 2×2 has
hybrid-vs-FDTD as gap 5, 0.76%/0.59% max|ΔT|; the waisted one never did). The
library refused the design point outright:

```
ValueError: Cut/stub at x=-9.5000 falls strictly inside a BODY waist transition
band ([-10.0000, -4.0000], 2.8 <-> 2.24 µm) — the multimode body's own width
varies there, so a section boundary cannot sit on a single cross-section.
```

**Mechanism** (`beef/solvers/fdtd.py`, `solve_section`): every FDTD section is
its own little simulation with PML on both sides, so it draws its waveguide
`dpml` **further than the section functionally needs, at both ends** — deliberate,
so the PML never meets an abrupt core→cladding step that would back-reflect. At
the *interior* seam that padding reaches `dpml` into the body. With
`body_stub = 0.25` and `dpml = 0.5` the drawn geometry needs uniform
cross-section out to **0.75 µm** from the body face:

> **`s = (mm_length − 2·Lt − Lw)/2  ≥  body_stub + dpml = 0.75 µm`**

**This is 3× tighter than step 9's `s ≥ 0.25`**, and it is a *different* claim:
step 9's bound keeps the FDTD sections **identical** to the standard device's;
this one is what makes them **constructible at all**. Step 6's grid enforced
only the loose one.

**Effect on the passer set — measured against the step-6 rows:**

| | count | |
|---|---|---|
| step-6 passers | 56 | |
| clear `s ≥ 0.75` (hybrid-buildable) | **52 (93%)** | |
| **excluded** | 4 | **including BOTH 29.50 µm candidates** (`s = 0.25`) |

The constraint is not restrictive in general — it bites precisely the points
step 6 ranked first, because "shortest among passers" drives `s` to its floor.
The 29.50 µm candidates were already failing phase at a converged taper basis
(step 5b); they now also cannot be cascaded.

**Shortest hybrid-buildable passers** (all 30.50 µm):

| Ww | Lw | Lt | s | L | IL dB | imb | phase° |
|---|---|---|---|---|---|---|---|
| 2.24 | 8.0 | 5.0 | 1.75 | 21.5 | 0.3557 | 0.1002 | −91.323 |
| 2.24 | **6.0** | **6.0** | **1.75** | **21.5** | 0.3880 | **0.0127** | −91.255 |
| 2.24 | 8.0 | 4.0 | 2.75 | 21.5 | 0.4185 | 0.0957 | −90.632 |

📌 **`Lw=6, Lt=6, L=21.5` (30.50 µm) is now the candidate on three independent
grounds**: 9% of the imbalance budget against the 29.50 µm points' 62–80%, a
band spread of **0.237°** against 1.732°, and constructibility under the
cascade. It is the geometry the hybrid-vs-FDTD comparison is being run at
(jobs `12797278` full FDTD → `12797279` hybrid).

Both scripts now carry the bound as an **up-front guard naming the reason**
rather than failing ~200 frames deep in meep, and the reference S-matrices are
geometry-keyed so a second design point cannot overwrite the first.

#### ✅ FULL-DEVICE FDTD AT 30.50 µm — job 12797278, and it PASSES where 29.50 µm FAILED

Run at `Ww=2.24, Lw=6.0, Lt=6.0, L=21.5` (`s=1.75`, the shortest
hybrid-buildable passer), res 40, 51 wl, `decay_tol=1e-5` satisfied at
t=379.8, `passive=True`, 8078 s on 48 ranks.

| quantity | 29.50 µm (`Lw=8,Lt=6`) | **30.50 µm (`Lw=6,Lt=6`)** | threshold | budget used |
|---|---|---|---|---|
| excess loss | 0.2036 dB | **0.2178 dB** | 0.4809 | **45%** |
| imbalance | 0.0976 | **0.0168** | 0.1371 | **12%** |
| phase mean-dev | **2.099° FAIL** | **1.911° PASS** | 2.0 | 96% |
| phase MAX-dev | 2.323° | 2.564° | 2.0 | fails strict |
| **verdict** | **FAIL (phase)** | **✅ PASS (all four)** | | |

**This is the first candidate to pass the pre-registered criterion on a
brute-force reference rather than on EME alone.** The extra 1 µm over the
29.50 µm point buys the phase pass outright, and cuts imbalance use from 71%
to 12%.

**The engine offsets reproduce, which is what makes them a systematic rather
than noise** — EME vs FDTD on this geometry: IL 0.3880 → 0.2178 (**−0.17 dB**,
against −0.19/−0.24 on the other two geometries) and phase −91.255 → −91.661
(**−0.41°**, against −0.40/−0.71). Same sign, same magnitude, three
geometries.

⚠ **Phase uses 96% of its budget** (1.911° of 2.0°), and the engine offset
alone is ±0.4–0.7°. It passes, but with no real margin — and its MAX-dev
(2.564°) still fails the strict reading of the pre-registered text, as every
candidate has.

⚠ **Sixth filename collision, same class.** `check4_res40_nwl51.json` carried
no geometry, so this run overwrote the L=20.5 summary. Recovered from
`reports/12743476.out` into `check4_Ww2.24_Lw8_Lt6_L20.5_res40_nwl51.json`
with provenance, and the path is now geometry-keyed — as are the S-matrix and
decay-snapshot names, which had the same defect.

#### ⚠ FINDING 2026-09-04: `N_MODES = 14` is rank-deficient at the waist — but only for ~6% of cells

**Read the correction below before quoting any count from this section.**

Cells in the waisted body occasionally cannot supply 14 distinct modes:

```
EMESolver.compute_modes: meow mode 7 duplicates already-accepted mode 6
  (self-overlap 1.000000 > 0.90) -- skipping it and reaching further up the
  spectrum. A rank-deficient basis would silently drop one channel of the 14
  requested.
EMESolver.compute_modes: only 13 of 14 requested modes were distinct among 14
  solved -- re-solving with 22 to reach further up the spectrum.
```

**Extent, measured correctly: ~7 events out of 120 cell-solves (~6%)** in job
12826121 (40 cells × 3 wavelengths), split as ~2 cells short by two modes and
~5 short by one.

❌ **CORRECTION — I got this badly wrong three times, and the cause is worth
recording.** These messages are emitted by **every MPI rank**, so on this
48-rank job each event appears **48 times** in the log. I counted raw log lines
and reported 240 / 288 / 336 "re-solves", built a growth rate from them, and
concluded first that *every* cell was rank-deficient, then corrected to *a
minority*, then un-corrected. All three readings were artifacts of a 48×
inflated proxy. **The true number is the line count divided by the rank count.**

⚠ **A decision was made on that bad metric.** Job 12797279 (the 11-wavelength
hybrid) was CANCELLED at 3:00 because the inflated count projected an overrun
of the 5 h `-ece` wall by ~15 min. At the true rate (5 events, not 240) that
projection was meaningless and the run would very likely have completed,
delivering the full 11-wavelength comparison. The re-run at 3 wavelengths
(12826121) was an unnecessary cost and a coarser measurement.

**The physics point survives the correction, at reduced scope.** `N_MODES=14`
is inherited from the STRAIGHT 2×2, where band 14 already sits below the
cladding index — i.e. 14 is AT the radiation boundary for a **2.8 µm** body.
The waist is **2.24 µm**, so the narrowest cells genuinely cannot fill a
14-channel basis and must reach up to 22. That is a real statement about
sizing the interface basis from the waist rather than the face, and it bears on
check 2's mode-identity question — but it affects a small minority of cells,
not the whole body, and it is NOT a significant cost driver.

✅ The guard behaved correctly throughout: it detected each rank deficiency and
re-solved rather than silently dropping a channel — the
`fdtd_vs_EME_mode_validation.md` item-0 protection working as designed.

📌 **Lesson for this repo's logs:** any per-rank log line is multiplied by the
rank count (48 on the 2-node MMI jobs). The same duplication was noticed and
reported for `[straight] total 36.00 um` a day earlier and then not applied
here. **Divide by the rank count before treating a log-line count as a
measurement.**

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

#### ⚠ FINDING 2026-09-08: the hybrid's own `C_TRANS` is NOT converged by 32 — IL diverges from FDTD past 16

Step 5a converged `C_TRANS=8` for an **EME-only** body (line 610). That is a
different quantity from `C_TRANS` inside the **hybrid** chain (FDTD stub |
EME ramp×C_TRANS | EME waist | EME ramp×C_TRANS | FDTD stub), which cascades
the ramp against real FDTD port sections instead of EME stubs on both ends.
Single wavelength (index 5, λ=1.5475 µm, band centre), fixed design point
(Ww=2.24, Lw=6.0, Lt=6.0, L=21.5). Job 12904661 (ct4/6/8/12/16) +
12962996 (ct24/32, forced onto `gedc_main` — neither fit the `-ece` 5 h cap).

| C_TRANS | IL (dB) | imb | phase (°) | ΔIL | Δimb | Δphase |
|---|---|---|---|---|---|---|
| 4 | 0.1699 | 0.0528 | −90.314 | — | — | — |
| 6 | 0.1688 | 0.0423 | −90.780 | −0.0011 | −0.0106 | −0.466 |
| 8 | 0.1845 | 0.0333 | −90.987 | +0.0157 | −0.0090 | −0.207 |
| 12 | 0.2076 | 0.0298 | −91.104 | +0.0231 | −0.0035 | −0.117 |
| 16 | 0.2120 | 0.0294 | −91.254 | +0.0044 | −0.0004 | −0.150 |
| 24 | 0.2274 | 0.0274 | −91.372 | +0.0154 | −0.0020 | −0.118 |
| 32 | 0.2415 | 0.0251 | −91.596 | +0.0141 | −0.0023 | −0.224 |

FDTD reference (full-device, job 12797278): IL 0.2178 dB, imb 0.0168,
phase −91.661°.

**Phase gap to FDTD shrinks a lot (0.407° → 0.289° → 0.065° from 16→24→32)
but the step sizes do not monotonically shrink** — Δphase went
−0.150 (12→16), −0.118 (16→24), −0.224 (24→32), i.e. it got *bigger* on the
last step rather than continuing to halve. **imbalance is still converging
monotonically toward FDTD** (gap 0.0126 → 0.0106 → 0.0083 at 16/24/32) but has
not flattened by 32.

❌ **IL does not converge toward FDTD past C_TRANS=16 — it overshoots and
keeps climbing away from it.** IL was approaching 0.2178 dB through ct12
(0.2076, under by 0.010) and ct16 (0.2120, under by 0.0058), then crossed over
and diverged: ct24 = 0.2274 (over by 0.0096), ct32 = 0.2415 (over by 0.0237).
Two full extra ramp refinements made the excess-loss match to FDTD *worse*,
monotonically, in the wrong direction.

**Honest read: `C_TRANS` in the hybrid is not converged at 32, and the IL
trend rules out simple discretization refinement as the whole story.** A
quantity that is supposed to be approaching a fixed reference from below and
instead sails past it and keeps climbing is not consistent with "not enough
ramp slices yet" — that failure mode saturates, it doesn't overshoot and
diverge. Candidates worth checking before spending more C_TRANS: mode-order
drift between adjacent `C_TRANS` values changing which mode the interface
picks up at the ramp/waist boundary cell-by-cell (this repo's known
mode-crossing failure mode — see `mode-crossing-correspondence-repair-plan`
memory), or the `N_MODES=14` rank-deficiency at the waist (line "FINDING
2026-09-04" above) interacting differently with more, thinner ramp cells.
**Not diagnosed — flagging so the C_TRANS=8 hybrid-vs-FDTD numbers already
reported (2.90%/3.22% bar/cross) are read as measured at an unconverged
setting, and so more C_TRANS is not the next thing to try.**

#### ✅ RESOLVED 2026-09-10: the IL divergence is the STAIRCASE, not the basis — and C_TRANS has an OPTIMUM, not a limit

Two hypotheses for the C_TRANS divergence were tested and **both falsified**,
which leaves a third that is directly supported by measurement.

**(1) N_MODES=10 does not fix it.** Re-ran the whole ct4..ct32 ladder at the
rank-clean basis (job 13075410, 7 tasks; 12-17 s per rung because the FDTD
sections were already cached from the nm10 wavelength run and C_TRANS only
re-slices the EME body):

| C_TRANS | IL nm10 | IL nm14 | imb nm10 | imb nm14 | phase nm10 | phase nm14 |
|---|---|---|---|---|---|---|
| 4 | 0.1750 | 0.1699 | 0.0524 | 0.0528 | −90.259 | −90.314 |
| 6 | 0.1719 | 0.1688 | 0.0427 | 0.0423 | −90.741 | −90.780 |
| 8 | 0.1875 | 0.1845 | 0.0339 | 0.0333 | −90.953 | −90.987 |
| 12 | 0.2035 | 0.2076 | 0.0299 | 0.0298 | −91.157 | −91.104 |
| 16 | 0.2078 | 0.2120 | 0.0295 | 0.0294 | −91.308 | −91.254 |
| 24 | 0.2234 | 0.2274 | 0.0275 | 0.0274 | −91.424 | −91.372 |
| 32 | **0.2415** | **0.2415** | 0.0256 | 0.0251 | −91.609 | −91.596 |

The ladders agree to ~0.004 dB throughout and land on **exactly 0.2415 at
ct32**. The overshoot past FDTD's 0.2178 is real and reproduces on a
rank-clean basis. **The ~50-58% rank-deficiency was a genuine defect worth
fixing — it cost a 15x runtime penalty — but it explains neither the C_TRANS
divergence nor the hybrid-vs-FDTD gap.**

**(2) The three quantities diverge from each other, which is the actual clue:**

- **phase CONVERGES** — −91.609 against FDTD's −91.661, a gap of **0.052°**,
  still shrinking
- **imbalance converges but plateaus short** — 0.0256 vs 0.0168, still moving
  −0.0019/rung
- **IL DIVERGES** — crosses FDTD around ct12-16 (0.2035, 0.2078) and then
  climbs away to 0.2415

**(3) The staircase explains exactly this pattern.** The stitched-field
measurement above found a **~2% |E|² discontinuity at EVERY ramp cell
boundary**, flat past 8 modes, i.e. the artificial step itself rather than
basis truncation. ct4 has 8 such boundaries; ct32 has 64. **Every added cell
is another artificial discontinuity that sheds power, so refining the ramp
ADDS spurious loss instead of removing it** — which is why IL climbs without
saturating, why the ct12-16 agreement with FDTD is a coincidental crossing
rather than convergence, and why phase (which does not accumulate that way)
converges cleanly.

📌 **Consequence for practice: `C_TRANS` has an OPTIMUM, not a floor.** Step
5a's "converged at 8" and the instinct to push it higher are both wrong for
the hybrid — past ~ct16 more cells make the loss worse. Any future waisted-body
cascade should treat C_TRANS as a bias/resolution trade to be tuned, not a
convergence parameter to be maximised. Non-uniform cell placement (the
`region_cell_lengths` hook already in `EMESolver`) is the principled fix:
concentrate cells where the width actually changes fastest instead of adding
boundaries uniformly.

#### ⚠ FOLLOW-UP 2026-09-09: the rank-deficiency rate in the hybrid's ramp is ~50-58%, not ~6%

Chased the mode-order-drift / rank-deficiency candidates above by counting
`EMESolver.compute_modes: only N of 14 requested modes were distinct ...
re-solving with 22` events straight from rank 0's own (unprefixed) log lines
in the ct16/24/32 jobs — no rank-multiplication division needed, exact counts:

| C_TRANS | distinct cell-solves | deficient (re-solved to 22 modes) | fraction |
|---|---|---|---|
| 16 | 18 | 9 | 50.0% |
| 24 | 26 | 15 | 57.7% |
| 32 | 34 | 19 | 55.9% |

**This is an order of magnitude higher than the "~6%" figure recorded above
(FINDING 2026-09-04) — and the two are not in conflict, they measure
different things.** That 6% was the EME-only full body (job 12826121, 40
uniform cells × 3 wavelengths, no FDTD ports). This count is the hybrid
chain's ramp/waist cells specifically, at one wavelength. Roughly half of
every distinct cross-section solved in the ramp cannot fill a 14-mode basis
and silently re-solves to 22 before the interface stitch — and **the fraction
holds essentially constant as C_TRANS grows** (50→58→56%), so refining the
ramp does not shrink the problem, it just multiplies the deficient-cell count
in step with the total (9→15→19, tracking 18→26→34 almost exactly).

**This is now the leading candidate for the IL divergence**, ahead of
mode-order drift (which was speculative; this is measured). A basis that is
14-wide at some interfaces and re-solved to 22-then-truncated at others is an
inconsistent interface stitch by construction, cell to cell, and there is no
reason its residual error should shrink monotonically as more such
inconsistent cells are added — which matches the observed overshoot-and-climb
IL trend better than a plain discretization-error story does. **Not yet
proven causal** (would need e.g. a run with `N_MODES` sized to the 2.24 µm
waist rather than inherited from the 2.8 µm face, per the same FINDING
2026-09-04 already flagged this as worth trying) — but this is the concrete
mechanism to test next, not another C_TRANS sweep.

#### ✅ FOLLOW-UP 2026-09-09: the waist's own rank-clean mode count is 10, not 14

Direct diagnostic at the 2.24 µm waist cross-section (res 40, λ=1.5475,
`guard_duplicate_modes=False` so no escalation masks the result — see
`self_overlap_diagnostics`-style raw check, `_solve_modes_raw` via
`compute_modes(..., guard_duplicate_modes=False)`), sweeping requested mode
count N and reporting the worst off-diagonal self-overlap among the raw N
modes:

| N | max off-diagonal overlap | clean? |
|---|---|---|
| 8 | 0.0794 | yes |
| 9 | 0.0794 | yes |
| 10 | 0.0794 | yes |
| 11 | 1.0000 | **no** |
| 12 | 1.0000 | no |
| 13 | 1.0000 | no |
| 14 | 1.0000 | no |
| 15 | 1.0000 | no |
| 16 | 1.0000 | no |

**The waist supports exactly 10 distinct guided/quasi-guided modes at this
wavelength — mode 11 onward is an EXACT duplicate (overlap 1.000000), not a
gradual degradation.** `N_MODES=14`, inherited from the 2.8 µm face where
band 14 sits right at the cladding index, asks the 2.24 µm waist for 4 modes
past what it can supply — confirming FINDING 2026-09-04's face-vs-waist
sizing point directly rather than by inference.

**This is the concrete next test the IL-divergence finding called for.**
Rerunning the hybrid C_TRANS/wavelength sweep at `N_MODES=10` (not 14) is
the way to check whether the interface-basis inconsistency (14-wide at some
cells, 22-then-truncated at others, per the follow-up above) was the actual
driver of the IL overshoot — at N_MODES=10 every ramp/waist cell should
request a basis it can natively fill, with no escalation and no
inconsistency to test against. Not yet run at 11 wavelengths; the user
chose to hold the hybrid rerun for this fix rather than spend more compute
sweeping an unconverged N_MODES=14 further.

#### ✅ DELIVERABLE 2026-09-09: stitched hybrid field, real multi-cell EME body

Built `stitch_mmi2x2_waist_hybrid_fields.py` -- unlike
`stitch_mmi2x2_hybrid_fields.py` (the straight device, which propagates the
body as a single diagonal `exp(2j*pi*neff*x/wl)` per mode, exact only for a
uniform cross-section), this drives the REAL multi-cell body (stub|ramp
×C_TRANS|waist|ramp×C_TRANS|stub) with a genuine bidirectional EME cascade:
the captured FDTD interface field is projected onto the first body cell's
basis, then propagated cell-to-cell via meow's own
`compute_propagation_s_matrices`/`compute_interface_s_matrices`/`propagate`
(the same primitives `mw.propagate_modes()` uses internally), reconstructed
in beef's device frame via `stack_meow_mode_fields()` per cell. Every
individual piece is machinery this repo had already validated elsewhere; new
here is driving it from a real captured excitation across 5 non-uniform
regions instead of a canonical single-mode one.

Run at production res=40, `C_TRANS=8`, and **`N_MODES=10`** (the rank-clean
count measured above, not the broken 14) -- job 13026733, ~2.3h on 1 rank:

```
projected onto 10 body modes; basis accounts for 1.0005 of the interface
  field's power
SEAM CONTINUITY at x = -10.500 um
  integrated |E|^2: FDTD side 3.8136e+03  body side 3.7350e+03  ratio 0.9794
```

Only ONE rank-deficiency event fired across all 19 distinct cell mode-solves
(down from the smoke test's 11 out of 11 at N_MODES=8/res=20) -- direct
confirmation that N_MODES=10 is close to the waist's natural mode count
everywhere along the ramp, not just at the narrowest cross-section.

**Field plot**: `example_output/verification/mmi_2x2_tapered/
hybrid_stitched_res40/fields/hybrid_stitched_waisted_res40_ct8_nm10_1.5475um.png`.
Qualitatively matches the full-device FDTD capture
(`fdtd_total_field_Ey_res40.png`, from the already-captured
`captured_full_device_res40_wl1.5475` data, no re-run): power launched on the
top arm (y=+0.75) crosses to the bottom arm through the body, landing near
the bottom output by the body's far end in both pictures.

#### ⚠ MEASURED 2026-09-09: ~2% field discontinuity at EVERY ramp cell boundary, and it does NOT shrink with more modes

The user spotted per-cell magnitude variation in the ramp on the stitched
field map and asked whether it was physical. Chased it properly rather than
by eye: the reconstruction inputs (per-cell mode slabs, cascade amplitudes,
n_eff, cell lengths) are now cached to `recon_cache_*.npz` by the stitch
script, so rendering questions are answerable in seconds instead of re-paying
the mode solves. Metric = relative jump in the transverse |E|² profile across
each interior cell boundary,
`sum|I_right - I_left| / (0.5*sum(I_right+I_left))`.

❌ **Two rendering explanations tested and FALSIFIED.** (1) A
normalization-convention mismatch (meow normalizes to unit POWER FLUX at
`meow/fde/meep.py:169`, the script was renormalizing to unit ∫|E|²dA) is real
and now fixed -- but measured, the per-cell factor is flat to 0.1%
(0.847→0.846), far too small to be visible. (2) Rendering the OLD and NEW
variants from IDENTICAL cached physics puts no visible discontinuity at the
cell boundaries in either; the only genuine rendering defect was texture in
the long cells, where a fixed 40-samples-per-cell rule gave 0.15 µm quads
across the 6 µm waist against 0.019 µm in a 0.75 µm ramp slice (now
60 samples/µm). A per-cell `argmax` slice rule was also wandering (z index 15
instead of 14 in 4 of 19 cells, all in ramps) -- fixed by pinning, but small.

✅ **What survives is physical**, at res20/N_MODES=8/C_TRANS=8:

| modes kept | max jump | mean jump |
|---|---|---|
| 2 | 0.0664 | 0.0455 |
| 3 | 0.0534 | 0.0396 |
| 4 | 0.0383 | 0.0270 |
| 5 | 0.0316 | 0.0214 |
| 6 | 0.0547 | 0.0204 |
| 8 | 0.0547 | 0.0201 |

**The jump falls until ~5 modes and then FLATTENS at ~2% mean -- modes 6, 7,
8 change it by nothing (0.0204, 0.0204, 0.0201).** So the residual is not
basis truncation; it is the staircase itself, i.e. the real cross-section
discontinuity between adjacent uniform EME cells. Per-boundary, the ramp
steps sit at 1-2.5% each, with one outlier: **the boundary between the two
NARROWEST ramp cells, immediately before the waist, jumps 5.5% -- ~3x
typical** -- the same place the waist rank-deficiency bites (modes near
cutoff, natural mode count 10).

**Why this matters beyond the picture.** The production cascade uses the SAME
cell decomposition and the SAME interface S-matrices, so this ~2%-per-boundary
staircase mismatch is present in the S-matrix path too -- the field map is a
picture of the error mechanism, not a separate artifact. It also reframes the
C_TRANS puzzle: raising C_TRANS shrinks each step but ADDS more steps, and if
per-step error does not fall faster than the step count grows, the total need
not converge -- which is exactly the observed IL overshoot-and-climb.
✅ **Confirmed at production settings** (res40, N_MODES=10, C_TRANS=8; job
13043451, 2h09m, re-measured from that run's own cache in seconds):

| modes kept | max jump | mean jump |
|---|---|---|
| 2 | 0.0613 | 0.0460 |
| 4 | 0.0389 | 0.0267 |
| 6 | 0.0520 | 0.0206 |
| 8 | 0.0510 | 0.0192 |
| 10 | 0.0510 | 0.0192 |

Mean 0.0192 against 0.0201 at res20/nm8, flat past 8 modes (0.0192 at both 8
and 10), and the SAME worst boundary -- cell 7→8, the two narrowest ramp
cells before the waist, at 5.1% against 5.5%. Ramp-internal mean 0.0201.
**Resolution- and mode-count-independent, so the staircase reading holds and
this is not a truncation or meshing artifact.**

**Caveat carried forward**: this is `C_TRANS=8`, which the FOLLOW-UP above
this one shows is not converged (IL still diverging through 32) -- and that
sweep was itself run at the broken N_MODES=14. The seam-continuity
improvement (0.9387 at N_MODES=8/res=20 -> 0.9794 at N_MODES=10/res=40) is
suggestive but confounds N_MODES and resolution in the same comparison, so it
is not yet a clean before/after. The next real test remains: rerun the
C_TRANS/wavelength sweep at N_MODES=10 and see whether the IL divergence
goes away.

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
