# Multi-excitation cascades — build plan

Opened **2026-08-20**. `mmi_2x2_hybrid.py` and `directional_coupler_hybrid.py`
both bypass `HybridSolver.solve()` and hand-roll an `einsum` cascade, because
their cross path pairs one arm's transposed matrix with the *other* arm's
excitation — which `solve()`'s one-call contract cannot express. A side effect is
that **neither can use the `calibrate_body_neff` flag**, since the flag only fires
inside `_solve_sections()`.

This is the library-integration half of `notes.md` Solver recommendation #2.

> **This doc is a PLAN, not a gap list.** It converts to a regression suite
> rather than a four-checks list.
>
> **The headline finding is MEASURED** (read directly off the source, with line
> references): routing the 2×2 through the flag would **not change a single
> number**. The cost accounting in step 3 is **ARGUED** until measured.

---

## ✅ CLOSED 2026-08-21 — all four steps done, and the result beat the bar

**Neither the 2×2 nor the directional coupler hand-rolls its cascade any
more.** Both excitations of each device now go through
`HybridSolver.solve_length_sweep()`, and both migrated examples are
**bit-identical** to what they replaced, each at **the same 1 FDTD simulation**
the hand-rolled path cost:

| device | vs. | `T` grids | complex coeffs |
|---|---|---|---|
| 2×2 MMI | shipped run + step 1's golden | 0.000e+00 | **0.000e+00** |
| directional coupler | shipped `length_sweep.npz` | 0.000e+00 | — |

Two library additions made that possible: **`MirrorSpec`** (step 3a — declares
a cross-excitation mirror with a parity vector) and **`calibrate_kwargs`**
(step 4 — lets the calibration take kwargs the section solver must not see,
which is what kept the DC's body cache valid).

Per device, before vs. after:

| | hand-rolled (before) | migrated (after) |
|---|---|---|
| `FDTDSolver.solve_section()` calls, both excitations | 1 | **1** |
| `section_mode_beta_details()` calls | 1 | **1** |
| `calibrate_body_neff` usable | ❌ | ✅ |
| symmetry arguments geometrically re-validated per call | ❌ | ✅ |
| agreement with the golden | — | **bit-identical** |

**Two plan assumptions turned out to be wrong, both in the good direction** —
see the steps below for the measurements:

1. **`Circuit` was not needed, and would have been worse.** Step 2's whole
   premise was to compose the two excitations with `beef.network.Circuit`.
   That is unnecessary: once step 3's `MirrorSpec` exists, each excitation is
   an ordinary single-excitation `solve_length_sweep()` call and there is
   nothing left to compose. And `Circuit` would have *cost* accuracy —
   `connect()` contracts with BLAS `@` while `solve_length_sweep()` uses the
   same `np.einsum` the hand-rolled path did, which is exactly why the result
   is bit-identical instead of the 4-eps step 1 measured for the `Circuit`
   route.
2. **The bar did not need relaxing after all.** Step 1 relaxed
   "bit-identical" to "closely aligns (<1e-13)" because BLAS-vs-einsum
   rounding made bit-identity unreachable *through `Circuit`*. Through
   `solve_length_sweep()` it is reachable, and was reached. The relaxed gate
   is kept as the tests' enforcement bar (it is the right bar for the
   `Circuit` equivalence tests, which still exist), but the migration itself
   is held to the stricter standard it actually met.

---

## Summary

| | status |
|---|---|
| `calibrate_body_neff` flag on `HybridSolver` | ✅ implemented, tested, default `False` |
| The 2×2 propagating its body on FDTD-measured β | ✅ **already happening** — hand-rolled, not via the flag |
| `mirror_pairs` reuse within one excitation | ✅ implemented (`hybrid.py:263`, validated per-call at `:1380`) |
| Composing two excitations' matrices in one `solve()` call | ❌ **not supported** |
| A primitive that can compose them after the fact | ✅ **`beef.network.Circuit` — exists since 2026-07-29, never used for this** |
| `Circuit` composition == the hand-rolled einsum | ✅ **MEASURED 2026-08-20** — 4.0 eps on the real device matrices, 0.8–11 eps synthetic; algebraically exact under exact arithmetic on bar *and* cross. `tests/test_multi_excitation_cascade.py`, 12 tests. Bar is **"closely aligns" (<1e-13 rel)**, relaxed from bit-identical |
| The real device's captured cascade matrices | ✅ job 12070199, COMPLETED in 38 s (all cache hits) |
| Step 2's run count (was ARGUED) | ✅ **MEASURED** — 2 solves vs the hand-rolled 1; the shared cache already recovers the input-side duplicate, so step 3 must recover exactly **one output-side solve**, not the 2+ this plan first estimated |
| Step 3's "fix the cache/validate ordering" trap | ✅ **already fixed in commit 29519d0** (this plan's own commit) and tested; plan text was stale |
| Declaring one section as another's transpose **across** excitations | ✅ **`MirrorSpec`, implemented 2026-08-21** (`solvers/hybrid.py`, exported as `beef.MirrorSpec`) — preserves the single-FDTD-run count |
| Composing two excitations without a hand-rolled einsum | ✅ two `solve_length_sweep()` calls + one shared `cache=`; no `Circuit` needed |
| The 2×2 using `calibrate_body_neff` | ✅ **now actually via the flag**, not the hand-rolled equivalent |
| Migrated result vs the golden | ✅ **bit-identical** (0.000e+00) incl. complex phase, 201 × 11 grid |
| `directional_coupler_hybrid.py` migrated too | ✅ **bit-identical** vs its shipped `length_sweep.npz`, 1 FDTD solve for both excitations |
| Calibration kwargs separable from section kwargs | ✅ **`calibrate_kwargs`**, added 2026-08-21 — the DC's blocker |
| Regression tests | ✅ `tests/test_multi_excitation_cascade.py` 30 tests (12 from step 1 + **18 new**), full suite **766** green |

Ordered work: **establish numerical equivalence (step 1) → migrate onto
`Circuit` (step 2) → generalise the mirror relation to recover the run count
(step 3) → retire the hand-rolled cascade (step 4).** **All four are now
COMPLETE**, with step 2's `Circuit` premise measured to be unnecessary — see
the CLOSED block above.

---

## The finding: switching would change nothing numerically

`HybridSolver`'s calibration path and the 2×2's hand-rolled one compute the same
quantity from the same source.

**`HybridSolver`** — `_calibration_source(k=1)` (`hybrid.py:1248-1270`) skips any
solver that has `section_mode_beta` (so the EME body qualifies for calibration),
then walks `(k-1, k+1)` and returns the first neighbour that *has* it, preferring
the left. For the 2×2's `[FDTD, EME, FDTD]` chain that is **the input
`FDTDSolver`**. `_calibrated_neff` (`hybrid.py:1272-1321`) then calls
`source.section_mode_beta_details(geometry, wavelengths, x0, x1, n_modes, lengths=self.calibrate_lengths, ...)`
over `x0, x1 = self.boundaries[1], self.boundaries[2]` and returns

```python
beta * np.asarray(wavelengths, dtype=float)[:, None] / (2.0 * np.pi)
```

**The example** (`mmi_2x2_hybrid.py:279-283`):

```python
calibration = input_solver.section_mode_beta_details(
    mmi, wavelengths, boundaries[1], boundaries[2], N_MODES, save_dir=save_dir,
)
beta = calibration.beta
body_neff = beta * wavelengths[:, None] / (2 * np.pi)
```

Same object, same boundaries, same `n_modes`, same conversion. `calibrate_lengths`
defaults to `None` (`hybrid.py:268`), which reaches `section_mode_beta_details` as
`lengths=None` — identical to the example passing nothing. The only structural
difference is that `HybridSolver` forwards a `passthrough` of
`monitor_center`/`monitor_size` when present in the section kwargs
(`hybrid.py:1291-1293`); the 2×2 passes neither to the body.

**So the flag is not a physics question — it is an API question.** Do not expect
the device's numbers to move, and treat any movement as a bug.

*source: `beef/solvers/hybrid.py:1248-1321`; `examples/mmi_2x2/python/mmi_2x2_hybrid.py:279-302`*

---

## Why `solve()` cannot express it today

The 2×2's cross path is `M_in0ᵀ · M_body · M_in1` — excited at the **top** arm,
measured at the **bottom**. Three distinct savings are stacked
(`mmi_2x2_hybrid.py:20-43`):

1. **Output section never simulated.** By mirror symmetry plus reciprocity,
   `output_idx=k` is the x-mirror of `input_idx=k`, so output matrices are input
   matrices transposed (`:266-267`, `.transpose(0, 2, 1)`).
2. **Bottom arm never simulated.** Launching the bottom arm is itself a mirror of
   launching the top, so each body mode's coefficient matches or flips sign by
   parity — the measured `MODE_SIGNS` vector (`:140`), applied at `:257`
   (`M_in0 = M_in1 * MODE_SIGNS[None, :, None]`). This is what drops input-side
   FDTD runs **from two to one**.
3. **Analytic length sweep** on the uniform body (`:298-302`).

`mirror_pairs` covers only (1), and only *within* one excitation.
`_validate_mirror_pair` (`hybrid.py:373-431`) samples
`geometry.width_at(x_k, **kwargs_k)` against `width_at(x_j, **kwargs_j)` about
`mirror_axis`. For the cross path, section 0 carries `input_idx=1` (top) and
section 2 carries `output_idx=0` (bottom) — diagonal, not mirror-symmetric — so
the check would **correctly reject** the pair. Nothing expresses (2) at all:
there is no way to declare "section 0 under `input_idx=0` equals section 0 under
`input_idx=1`, times a parity vector."

Consequence: naively porting the example to `solve()` today forfeits (1) and (2)
for the cross path, pushing input-side FDTD runs from 1 to 2+ and adding a real
output-section simulation. At ~1.5–2.5 h per 3-D run for this device
(`running_jobs.md:116-121`), that is hours, not noise — and it would invalidate
the manuscript's claim that the whole response comes from a single FDTD
simulation.

---

## Step 1 — Equivalence harness

Before changing anything, pin the current answer. The existing cached section
results make this nearly free.

- Capture the current hand-rolled `M_in1`, `M_in0`, `M_body`, and the composed
  bar/cross matrices for the shipped configuration.
- These become the regression target for steps 2–4. **Bar relaxed 2026-08-20
  from "bit-identical" to "closely aligns"** — see the measured block below for
  why bit-identity is both unreachable and useless here.
- Follows the precedent of
  [`resolved/calibrate_body_neff_mode_solve_skip_validation.md`](calibrate_body_neff_mode_solve_skip_validation.md),
  which verified the mode-solve skip bit-exact against a reconstructed old path.
- Fake-solver level is enough for the composition algebra, so most of this
  belongs beside `tests/test_hybrid.py` (35 tests) and `tests/test_network.py`
  (29) and needs neither meep nor meow.

### MEASURED 2026-08-20 — DONE, and it relaxed this plan's bar

`tests/test_multi_excitation_cascade.py` (12 tests, meep/meow-free, all
passing; full fast suite green on the same tree).

**★ The "bit-identical, not within tolerance" bar above is NOT reachable as
literally written, and step 2's "done when" inherited the same error.** The
reason is not that the migration is wrong — it is that `connect()` contracts
with `@` (BLAS zgemm) while the example uses `np.einsum`, which sum in
different orders. `np.einsum` vs `@` on the bare 14-element contraction is
already not bit-identical (max |Δ| **3.79e-15**), so the discrepancy is
inherited from numpy, not introduced by `beef.network`.

What **is** true, and is the stronger claim:

- **`Circuit` is algebraically exact.** Fed matrices whose every intermediate
  product and sum is exactly representable in float64 (small integers, all
  intermediates ≤ 1.27e4 against float64's exact-integer range of 2⁵³), the
  `Circuit` → `reduce()` → `to_smatrix()` route reproduces the example's einsum
  chain **bit for bit** — on the bar path *and* the cross path, and across the
  201-length sweep grid. `network/ops.py:148-152`'s claim that the star product
  "reduces exactly to ordinary transfer-matrix multiplication" is now MEASURED,
  not argued.
- **On real float data the residual is a few ulp**: max relative difference
  **5.70e-16 (2.6 eps)** at the shipped `n_wl=11`/`N_MODES=14`, and 0.8–11 eps
  across `N_MODES` 4–64. It does not grow meaningfully with mode count.
- **It is not the star-product machinery.** `from_transfer_matrix()` zero-fills
  every self-term, so `connect()`'s systems are exactly `I - 0 @ 0 = I`, and
  `np.linalg.solve(I, X)` returns `X` bit-for-bit (pinned by a test, so nobody
  later goes hunting for a bug in `ops.py`).

**★ THE BAR, as relaxed 2026-08-20 at the user's direction: "closely aligns",
enforced as max relative difference < 1e-13** (`FLOAT_CHAIN_RTOL` in the test
module). Two measured reasons bit-identity was dropped rather than merely
missed:

1. **It is unreachable** — the BLAS-vs-einsum rounding above, inherited from
   numpy.
2. **It buys no detection power.** Mutation-tested against this harness, every
   realistic migration bug lands at *O(1) relative error*: using the cross
   composition where the bar belongs → 2.02, conjugating the input → 4.49,
   dropping one of 14 modes → 2.58. The gate at 1e-13 catches all three with
   ~13 orders to spare. The only thing bit-identity additionally catches is a
   1-ulp nudge, which is not a defect — while making the suite hostage to the
   BLAS build, thread count and numpy version.

Headroom: real device matrices agree to **4.0 eps** (job 12070199, full 201×11
grid), synthetic to 0.8–11 eps — about two orders inside the gate.

The exact-arithmetic tests are **retained but demoted from the bar to a failure
localizer** (`TestCircuitIsAlgebraicallyExact`). Under exact arithmetic both
routes are exact, so a difference there is necessarily algebraic rather than
numerical — if the gate fails, that class says which kind of failure it is.
Nothing should be blocked on it alone.

### A blind spot this harness closes

`tests/test_network.py::TestConnect::test_mmi_2x2_chain_matches_manual_einsum`
already existed and covers the **bar path only**. Measured: **the bar path
cannot detect a transposed body matrix.** It is `M_in1ᵀ · B · M_in1`, a
quadratic form whose 1×1 result equals its own transpose, so `B` and `Bᵀ`
produce the *identical* scalar — verified bit-identical, while the same
mutation moves the cross path by 9.99e3. A step-2 migration that transposed
the body would therefore sail past a bar-only check.

Note the shipped configuration runs `CALIBRATE_BODY=True`, which makes the body
**diagonal** (hence symmetric), so the bug would be invisible in production
today — but the uncalibrated body is not symmetric, so this is a live hazard
for step 2 rather than a curiosity. Three other mutations (cross-for-bar,
conjugated input, one mode dropped) are caught by both halves of the bar; that
mutation check is what established the transpose gap is real and not an
artifact of the test data.

That existing test's docstring claimed the match was "bit-for-bit (to floating-
point noise)" while asserting `atol=1e-12`; corrected in place to say what it
actually checks and to point at the new module.

### The real device's matrices — ✅ CAPTURED, job 12070199

`examples/verification/python/capture_mmi_2x2_cascade_golden.py`
(+ `batch/verification/run_capture_mmi_2x2_cascade_golden.batch`) →
`example_output/verification/mmi_2x2_cascade_golden/cascade_golden.npz`.
**COMPLETED 2026-08-20 20:00:17 in 38 s, exit 0:0** — the whole shipped
configuration was cache-resident, so this cost essentially nothing. It writes
only to its own directory and never touches
`example_output/mmi_2x2/mmi_2x2_hybrid/`, which the GDS acceptance run (job
12069544) reads at its own step 7.

Contents (read from the `.npz`, not the job log): `M_in1` (11, 14, 1), `M_in0`
+ `M_in0_derived` (11, 14, 1), `M_body` (201, 11, 14, 14), `coeffs_bar`/
`coeffs_cross` and `circuit_bar`/`circuit_cross` (201, 11), `T_bar`/`T_cross`
(201, 11), plus `git_sha` 29519d0 and the working-tree diffstat, so the golden
records which library state produced it.

**The synthetic result reproduces on the real device matrices.** einsum vs
`Circuit`, recomputed from the saved arrays: max relative difference
**8.885e-16 (4.0 eps)** on the bar path and **8.844e-16 (4.0 eps)** on the
cross path — the same few-eps rounding the unit tests predicted from synthetic
data, on the actual 201 × 11 grid.

Two notes for whoever picks up step 2:

- A plain `lookup()` probe of the three shipped calls **missed** on all three
  beforehand, which was *not* proof the cache was cold: `solve_section()` falls
  back to a full-matrix `lookup()` and then `lookup_superset()`
  (`fdtd.py:2924-2943`). The 38 s runtime settled it — the supersets hit. Do not
  size a job off a hand-built `lookup()` probe.
- `CAPTURE_REAL_M_IN0=1` adds the second input-side FDTD run that
  `SYMMETRIC_MODE_REUSE` exists to avoid, and is **still not run**. **Step 3(b)
  needs it**, because deriving or validating `MODE_SIGNS` requires an
  independently measured `M_in0` to check against; the golden's current `M_in0`
  is the `MODE_SIGNS`-derived one, identical to `M_in0_derived`. The script
  already prints the implied sign vector beside the hard-coded one under that
  flag. Budget it as a real ~1.5–2.5 h 3-D run, not a cache hit.

### Side benefit — this closed an open issue in the GDS validation

The same capture recomputed the parametric 2×2's `T_bar`/`T_cross` on the
current library, which
[`gds_geometry_hybrid_cascade_validation.md`](gds_geometry_hybrid_cascade_validation.md)
step 5 needed in order to separate a geometry-backend difference from the seven
uncommitted library edits behind its 2026-08-18 reference. Both grids came back
**bit-identical** to that reference (max |Δ| 0.000e+00). Caveat recorded there:
the capture was cache-served, so this re-verifies everything downstream of the
cached section matrices but not the meep/FDTD layer itself.

---

## Step 2 — Migrate composition onto `beef.network.Circuit`

**The primitive already exists and has never been used for this.**
`notes.md:583-587`:

> "`beef.network.Circuit` can register both excitations' transfer matrices and
> extract every port-combination coefficient from one `reduce()` +
> `to_smatrix()` call — see Network below. Neither example has been migrated to
> it; this item is about a dedicated `HybridSolver.solve()` API, which is still
> not built."

So the cheap route is not an API extension at all:

1. Call `HybridSolver.solve()` (or `solve_length_sweep()`) **once per excitation
   configuration** — each call has a single excitation, satisfying the one-call
   contract, so **`calibrate_body_neff=True` works unmodified**.
2. Share one `cache=` dict across the calls. Sections are memoized by
   `(solver identity, x0, x1, n_in, n_out, kwargs)` (`hybrid.py:1045-1049`), so
   every section whose kwargs do not change between calls is simulated once.
3. Compose the per-excitation results with `Circuit` —
   `from_transfer_matrix()` (`network/convert.py:17-32`) wraps each
   `[M, n_out, n_in]` matrix as a degenerate 2-port `SMatrixBlock`, then
   `reduce()` + `to_smatrix()`.

Note what `Circuit` does **not** buy here: `from_transfer_matrix` zero-fills
reflection and reverse transmission (`convert.py:29-31`), so a chain built purely
from `HybridSolver` matrices *"reduces exactly to ordinary transfer-matrix
multiplication"* (`ops.py:148-152`). It is a bookkeeping win — port combinations
handled in one place instead of a hand-rolled einsum — not new physics. Say that
plainly in the doc so nobody reads it as adding reflections.

**Cost, ~~ARGUED~~ MEASURED 2026-08-20:** this route recovers saving (3) and,
with `mirror_pairs` declared, saving (1) for the *bar* path where the pairing is
genuinely mirror-symmetric. It does **not** recover (2), nor (1) for the cross
path.

Counted at fake-solver level (`TestStep2RunCount`, the right level for a
"how many solves happen" question):

| | input-side solves | output-side solves | total |
|---|---|---|---|
| hand-rolled cascade (today) | 1 | 0 | **1** |
| step 2 via two `solve()` calls + shared cache | 1 | 1 | **2** |

**This REFINES the estimate in "Why `solve()` cannot express it today" above,
downward.** That section predicts porting would push "input-side FDTD runs from
1 to 2+ **and** add a real output-section simulation" — i.e. 3 or more. Measured,
the shared `cache=` already recovers the input-side duplicate on its own: both
excitations drive the *same* top arm, so section 0's `geometry_kwargs` are
identical between the two calls and it is simulated **once**
(`test_a_shared_cache_reuses_the_input_section_across_excitations`). The real
gap is exactly **one output-side solve**, and that is all step 3 has to recover.

Also confirmed, rather than assumed: the mirror guardrail **correctly rejects**
the cross pairing rather than silently transposing the wrong section —
`mirror_pairs[2]=0: center_at mismatch at x=13.75 vs mirrored x=-13.75: -0.75 vs
0.75 µm`, exactly the failure the plan predicts.

> **Methodology warning for anyone extending these tests.**
> `tests/test_hybrid.py`'s `FakeGeometry` ignores `**kwargs` entirely, so
> `input_idx`/`output_idx` do not move `width_at`/`center_at` and the mirror
> guardrail cannot tell the top arm from the bottom. A first pass at this
> measurement used it and reported that the non-mirror cross pairing raised no
> error and that step 2 cost nothing — both artifacts. Every cross-excitation
> question needs the branching `TwoArmGeometry` in
> `tests/test_multi_excitation_cascade.py`, which is guarded by its own
> `test_the_two_arm_geometry_is_not_flat`.

**Done when:** the migrated 2×2 **closely aligns** with step 1's captured
matrices — max relative difference < 1e-13, step 1's gate — with
`calibrate_body_neff=True` and no hand-rolled einsum, and the run count is
recorded honestly even if it went up.

> Relaxed 2026-08-20: this previously said "bit-identically", which step 1
> measured to be both unreachable (`connect()` contracts with BLAS, the example
> with `np.einsum`) and pointless (every realistic bug lands at O(1) relative
> error). See step 1's bar block, and the body-transpose blind spot any
> bar-only check would miss.

### ✅ DONE 2026-08-21 — but NOT via `Circuit`, and the run count went DOWN

**`Circuit` is not part of the migration.** This step assumed the two
excitations would have to be composed after the fact. They do not: once step 3's
`MirrorSpec` can express the cross pairing, each excitation is an ordinary
single-excitation `solve_length_sweep()` call that already returns the finished
`[n_L, n_wl, 1, 1]` coefficients. There is no leftover composition for `Circuit`
to do. The migrated example imports it not at all.

That is a **better** outcome than this step planned, on two counts:

- **Accuracy.** `solve_length_sweep()` cascades with the same `np.einsum` the
  hand-rolled path used, so the migration is **bit-identical** rather than
  4-eps-close. Routing through `Circuit` would have introduced the
  BLAS-vs-einsum residual step 1 measured, in exchange for nothing.
- **Cost.** This step's own measured estimate was 2 solves (1 input + 1
  output) against the hand-rolled 1. With `MirrorSpec` the output-side solve
  is recovered too, so the migrated total is **1**, not 2. The table at the
  top of this doc is the measurement.

**`Circuit` keeps its place** — `tests/test_multi_excitation_cascade.py`'s
step-1 classes still pin that `Circuit` reproduces the same algebra, and it
remains the right primitive for a genuinely multi-port network (one where
ports actually interconnect rather than chain). It simply is not what this
device needed.

**Measured, `cache_only=True` against the live cache:**

| migrated vs. | `T_bar` | `T_cross` | `coeffs_bar` | `coeffs_cross` |
|---|---|---|---|---|
| shipped run's `mm_length_sweep.npz` | 0.000e+00 | 0.000e+00 | — | — |
| step 1's golden (job 12070199) | 0.000e+00 | 0.000e+00 | **0.000e+00** | **0.000e+00** |

All bit-identical. The complex-coefficient rows matter: `|·|²` would hide a
global phase error, and the golden is the only place the pre-migration complex
coefficients were preserved.

> **A methodology trap this measurement walked into first, recorded because it
> cost a real login-node FDTD run.** The verification harness counted solves by
> **subclassing** `FDTDSolver`. That breaks two things at once:
> `cache/config.py:259` fingerprints on `type(solver).__name__`, so every
> lookup missed; and `hybrid.py:_supports_cache_only()` decides whether to
> forward `cache_only` by **introspecting the method signature**, which a
> `(*a, **k)` override erases — so the guard was silently dropped and a
> genuine 3-D solve started time-stepping. Both are the library behaving
> correctly. **Instrument with `functools.wraps` (which `inspect.signature`
> follows), never a subclass**, and assert `cache_only` arrived rather than
> trusting it. The `_supports_*` introspection is a general sharp edge: any
> wrapper around `solve_section` also silently disables `save_dir` and
> `neff_override` forwarding.

---

## Step 3 — Generalise the mirror relation to recover the run count

This is what makes the manuscript's single-FDTD-simulation claim true through the
supported API. Two sub-problems:

**(a) Cross-excitation mirror pairs.** Today `mirror_pairs` maps section → section
under one excitation. It needs to express "section `j` under excitation
configuration `A` supplies section `k` under configuration `B`", with
`_validate_mirror_pair` checking the geometry with each side's own kwargs — which
it already does (`hybrid.py:411-429`); what is missing is the *declaration*, not
the check.

**(b) The parity sign vector.** `MODE_SIGNS` (`mmi_2x2_hybrid.py:140`) is currently
a hard-coded array, measured once against a real two-arm FDTD run and reused
thereafter. Options:

- **Accept it as a declaration** — a `mode_signs=` argument alongside the pair.
  Simple, but it is unvalidated input; a stale vector after a geometry change
  silently corrupts the cross path. The example already warns about this: set
  `SYMMETRIC_MODE_REUSE=False` "to re-derive MODE_SIGNS after a geometry change"
  (`:29-33`).
- **Derive it from geometry parity.** Each body mode is even or odd about the
  mirror axis, so the sign is determined by the mode's own symmetry — derivable
  from the mode fields `section_mode_fields()` already returns (`fdtd.py:1977`),
  the same fields `beef.mode_identity` consumes.
- **Validate a declared vector cheaply.** Keep the declaration but verify it
  against mode parity at solve time, warning or raising on mismatch.

Recommendation: **derive it**, falling back to a declaration with validation.
Deriving removes the stale-vector failure mode entirely and is the kind of check
this repo already favours — the parity is a property of the modes, not of the
script.

### ✅ 3(a) DONE 2026-08-21 — `MirrorSpec`

`beef.MirrorSpec` (`solvers/hybrid.py`) is the declaration this step asked
for. A `mirror_pairs` value may now be a `MirrorSpec` instead of a bare int:

```python
mirror_pairs={2: MirrorSpec(section=0,
                            partner_kwargs={"input_idx": 0},
                            mode_signs=MODE_SIGNS)}
```

giving `M_2 = ((M_0 * mode_signs)ᵀ)`. It says exactly what the plan wanted:
*section 2 under this call's excitation is section 0 mirrored, but section 0
as it would have been solved on the other arm.* The bare-int form is
untouched — `mirror_pairs` stays a plain `{int: int}` map (which
`plot/geometry.py:569` reads), with the richer objects on a parallel
`mirror_specs`.

**The declaration is checked, not trusted — and this is where the design
work was.** `partner_kwargs` is dangerous in a way a bare int is not: the
existing x-mirror check compares section k against section j *under
partner_kwargs*, so it never looks at what section j was actually solved
with. On its own that makes `partner_kwargs` an unchecked way to pair any
two unrelated branches. So a **second** check was added,
`_validate_two_arm_symmetry()`: the partner's two excitations must
themselves be transverse mirror images — same width at every sampled x,
centres summing to `2 * transverse_axis`. That is the precondition that makes
a parity vector meaningful at all, and it is what closes the hole.

Both checks run **on every solve() call, including a cache hit** — the
requirement this plan's step 3 note preserved from the (already-fixed)
ordering trap, now covered for the new path too by
`test_the_two_arm_check_runs_on_a_cache_hit_too`.

Declaration-time validation also pins: `mode_signs` entries are ±1 only (a
non-unit value would silently *rescale* rather than reinterpret);
`mode_signs` without `partner_kwargs` is refused rather than ignored; and the
vector's length must match the partner's output-interface mode count.

**15 new tests** in `TestStep3CrossExcitationMirror`, meep/meow-free. The
load-bearing ones: the exact pairing the bare-int guardrail *correctly*
rejects is accepted as a spec; both excitations cost **one** input solve and
**zero** output solves; the composed result equals the hand-rolled cascade on
**both** paths (the bar path alone cannot see a transposed body — step 1's
blind spot); an all-`+1` parity vector gives a measurably different answer, so
nothing passes with `mode_signs` ignored; and a `partner_kwargs` naming the
same arm is rejected.

### ✅ 3(b) PREREQUISITE MEASURED 2026-08-21 (job 12159086) — the declared vector is CORRECT

`CAPTURE_REAL_M_IN0=1` ran and landed: **COMPLETED 45:01, exit 0:0**, well
under the ~1.5–2.5 h budget. It adds the second input-side FDTD arm, giving an
independently measured `M_in0` to check the hard-coded vector against.

**All 14 signs agree.** Implied by the measured pair, and hard-coded, are the
identical vector `[1,-1,1,-1,1,-1,1,-1,-1,-1,-1,1,-1,1]`; the script's own
`agree:` line prints `True`, and a per-mode re-derivation from the saved `.npz`
confirms the sign for **every** mode, not just the dominant ones.

**But the parity relation is a statement about signs, not amplitudes — and the
amplitude half only holds where it matters.** Per-mode, measured `M_in0`
against the `MODE_SIGNS`-derived one at λ=1.530 µm:

| modes | \|M_in1\| | max rel err vs derived |
|---|---|---|
| 0-2 (dominant) | 0.42 – 0.70 | 1.7e-05 – 2.6e-05 |
| 3-6, 10, 13 | 9.1e-03 – 7.9e-02 | 9.5e-05 – 1.8e-03 |
| **7, 9, 11, 12** | 2.2e-05 – 4.8e-05 | **1.0e-01 – 2.5e-01** |
| **8** | **1.3e-06** | **2.8e+00** |

The headline `max rel = 9.576e-01` in the job log is this effect: it is
mode 8, whose amplitude is 1.3e-06 against the dominant mode's 0.70 — a
factor of 5e5 below it. These are noise-floor modes, and the relation being
"wrong" by 280 % on a 1e-06 coefficient is not a defect in the parity
argument, it is the amplitude of a mode that carries no power.

**The cost of the assumption, quantified end-to-end.** Rebuilding the cross
path from the *measured* `M_in0` instead of the derived one moves the final
answer by:

| 201 × 11 grid | max abs Δ |
|---|---|
| `T_cross` | **1.585e-05** |
| `coeffs_cross` | 2.128e-05 |
| `T_bar`, `coeffs_bar`, `M_in1`, `M_body` | 0.000e+00 (bit-identical — the bar path never uses `M_in0`) |

So `SYMMETRIC_MODE_REUSE=True` costs **~1.6e-05 on T_cross** and saves a
~45-minute 3-D FDTD run. That is the number the example's own comment
estimated at "≤1.1e-4 relative" — now measured directly rather than inferred.

> **Two goldens now, deliberately.** The run rewrites `cascade_golden.npz` in
> place, and its cross columns are *not* the same quantity as before (they now
> come from the measured `M_in0`). Rather than let one filename mean two
> configurations, the measured-arm capture is kept as
> **`cascade_golden_real_m_in0.npz`** and `cascade_golden.npz` was restored to
> the `SYMMETRIC_MODE_REUSE=True` capture that is the regression target for
> steps 2-4. Re-verified after restoring: the migration is still bit-identical
> against it on all four arrays. This is why the file was backed up before
> submitting.

### ⏸️ 3(b) ITSELF STILL NOT DONE — `MODE_SIGNS` is declared, not derived

`MirrorSpec` takes the vector as an argument, which is this step's
**"accept it as a declaration"** option, hardened with the geometric
precondition above but *not* with mode-parity derivation. The plan
recommends deriving it, and that recommendation still stands.

**Why it was not done here, honestly:** deriving parity from
`section_mode_fields()` can be written, but it cannot be *verified* without
an independently measured `M_in0` to check the derived vector against — and
that is exactly the `CAPTURE_REAL_M_IN0=1` run this doc's step 1 flags as
**still not run** and budgets at a real ~1.5–2.5 h 3-D job. Shipping a
derivation whose output nothing had been checked against would put an
unmeasured claim into the library, which is the thing this repo's validation
standard exists to prevent.

What did land instead, and reduces the stale-vector risk meaningfully:

- The **geometric precondition is now enforced** every call. A geometry
  change that breaks the two arms' symmetry — the scenario the example's own
  warning is about — now raises instead of silently corrupting the cross
  path. It does not catch a stale vector on a *still-symmetric* geometry.
- The **re-derivation path is now first-class**. Setting
  `SYMMETRIC_MODE_REUSE = False` in `mmi_2x2_hybrid.py` no longer just
  "simulates both arms": it configures the cross solver as a plain
  `mirror_pairs={2: 0}` on an `input_idx=0` excitation, and then reads the
  two input-section matrices straight off
  `HybridSolver.last_section_matrices` (a new attribute, added for exactly
  this) to print the implied parity vector beside the hard-coded one, plus
  the ratio's imaginary/real magnitude as a sanity check that the relation
  really is a real ±1. Previously this required re-issuing the
  `solve_section()` calls by hand — i.e. a diagnostic that reconstructs the
  call it is diagnosing, and can therefore diagnose a different one.

**Remaining work — now UNBLOCKED.** The regression target exists:
`cascade_golden_real_m_in0.npz` holds a measured `M_in0` (job 12159086), so
field-parity derivation can finally be written *and verified* rather than
written and hoped for.

**But the priority dropped, and that is a real finding, not an excuse.** This
step's whole motivation was the stale-vector failure mode — "a stale vector
after a geometry change silently corrupts the cross path." Two things now
mitigate that:

1. The declared vector is **measured correct on all 14 modes**, so there is no
   live defect to fix.
2. `MirrorSpec` **enforces the geometric precondition every call**, so the
   geometry change that would strand the vector now raises instead of passing
   silently. That was the actual hazard.

What derivation would still buy: removing the last hand-maintained constant,
and covering the residual case of a geometry that stays transversely
symmetric while its *mode ordering* changes (e.g. a mode crossing — cf.
[[mode-crossing-correspondence-repair-plan]]), which the precondition check
cannot see. Worth doing; no longer urgent.

Implementation note for whoever does: derive from parity of the mode fields
about the transverse axis (`section_mode_fields()`, the same fields
`beef.mode_identity` consumes), and gate acceptance on the strong modes only
— the table above shows modes at the 1e-06 noise floor disagree by up to 280 %
in amplitude, so a derivation validated on *all* modes equally would fail for
reasons that have nothing to do with parity.

**~~Guard against a real trap~~ — ✅ ALREADY FIXED, this sub-item is closed.**
This plan said `_solve_sections()` `continue`d on an in-memory cache hit
*before* reaching `_validate_mirror_pair`, so a cross-excitation pair served
from `cache=` would skip its geometric check — and step 2 deliberately shares
one cache across excitations, making that the common path.

Checked 2026-08-20: the ordering is already correct in the working tree **and
in `HEAD`**. `_validate_mirror_pair` is called at `hybrid.py:1396-1400`, ahead
of the cache short-circuit at `:1413-1418`, under a comment that opens "Before
the cache short-circuit, not after". `git log -S` dates it to commit
**29519d0** — *"Audit the manuscript against the library: one correctness fix,
three build plans"* — i.e. the very commit that created this document. The
"one correctness fix" in that message **is** this fix; the plan text was
written from the pre-fix state and never updated. It is covered by
`tests/test_hybrid.py::test_mirror_guardrail_still_runs_on_a_cache_hit`.

Nothing to do here. The requirement itself still stands for anything step 3(a)
adds: a cross-excitation pair must re-validate per call, not once at
construction.

---

## Step 4 — Retire the hand-rolled path

Once step 3 lands: delete the manual cascade (`mmi_2x2_hybrid.py:249-302`), the
`MODE_SIGNS` constant if step 3 derives it, and the vestigial `output_solver` that
exists only so the section map labels correctly and is never run (`:170-173`).
Same treatment for `directional_coupler_hybrid.py`.

Also retires a side wart `notes.md:576-582` records: because `solve()` derives each
section's length from `boundaries`, it cannot get section lengths wrong, whereas
the hand-rolled path can — so each example currently asserts its own length budget
in script-level code that a `solve()`-based version would not need.

Keep the module docstrings honest afterwards: `mmi_2x2_hybrid.py:41-43` currently
explains *why* the einsum is hand-rolled and will be stale.

### ✅ DONE 2026-08-21 for the 2×2

Deleted from `mmi_2x2_hybrid.py`: the two manual `solve_section()` calls, the
`M_in0 = M_in1 * MODE_SIGNS` derivation, the `M_out_bar`/`M_out_cross`
transposes, the manual `section_mode_beta_details()` +
`beta * wl / 2π` conversion, the direct
`solve_section_length_sweep(neff_override=...)` call, and the three-line
`einsum` cascade. In their place: two `solve_length_sweep()` calls sharing one
`cache=`, and `calibration = hybrid_bar.last_calibrations[1]` for the
diagnostic plots.

Kept, deliberately:

- **`MODE_SIGNS`** — step 3(b) did not derive it, so it survives as the
  `MirrorSpec` argument. Retiring it is 3(b)'s job, not this step's.
- **`output_solver`** — no longer "vestigial". It is a real member of the
  declared chain (it sets section 2's mode counts and labels the section map);
  it simply never reaches its `solve_section()` because `mirror_pairs`
  short-circuits it on every call. The comment now says that rather than
  calling it an artifact.
- The module docstring's hand-roll rationale is replaced with the `MirrorSpec`
  declaration it became.

### ✅ ALSO DONE 2026-08-21 — `directional_coupler_hybrid.py`, after clearing the blocker

**Migrated and bit-identical.** Verified `cache_only=True` against the shipped
`example_output/directional_coupler/directional_coupler_hybrid/length_sweep.npz`
(2026-08-18):

| 201 × 11 grid | max abs Δ | max rel Δ | bit-identical |
|---|---|---|---|
| `T_through` | 0.000e+00 | 0.000e+00 | **yes** |
| `T_cross` | 0.000e+00 | 0.000e+00 | **yes** |

at **1** `solve_section()` call and **1** `section_mode_beta_details()` call
across both excitations — the hand-rolled path's own count. The declarations:

| path | excited | measured | declaration |
|---|---|---|---|
| through | `input_idx=0` | `output_idx=0` | `mirror_pairs={2: 0}` |
| cross | `input_idx=0` | `output_idx=1` | `MirrorSpec(section=0, partner_kwargs={"input_idx": 1}, mode_signs=MODE_SIGNS)` |

**The blocker below was real and is now fixed in the library**:
`HybridSolver` gained **`calibrate_kwargs`**, extra kwargs passed only to the
calibration source's `section_mode_beta()` call and never to the section's own
solver. `None` (default) reproduces the historical auto-extract exactly;
explicit entries win per key. Three tests pin it — reaches the calibration and
not the section, defaults unchanged, explicit beats inferred. The DC passes
`calibrate_kwargs=_iface_mon` while keeping section 1's `geometry_kwargs`
empty, which is what preserves the body's cache hit.

`SYMMETRIC_MODE_REUSE=False` gets the same treatment as the 2×2's: a plain
`mirror_pairs={2: 0}` on an `input_idx=1` excitation, with the implied parity
vector read off `last_section_matrices`.

<details>
<summary>The original 2026-08-21 investigation that identified the blocker</summary>

**Not done: `directional_coupler_hybrid.py` — and it is NOT a drop-in.**
Investigated 2026-08-21 rather than left as a guess. The good news first: the
DC's composition is structurally identical to the 2×2's — `M_in0` solved,
`M_in1 = M_in0 * MODE_SIGNS` derived (`MODE_SIGNS = [1,-1,-1,1,1]`, `N_MODES=5`),
`M_out_through`/`M_out_cross` as the two transposes, the same three-line
einsum — and `DirectionalCoupler` implements the same `input_idx`/`output_idx`
convention on `width_at`/`center_at`/`section_*`. So `MirrorSpec` maps onto it
directly:

| path | excited | measured | declaration |
|---|---|---|---|
| through | `input_idx=0` | `output_idx=0` | `mirror_pairs={2: 0}` |
| cross | `input_idx=0` | `output_idx=1` | `MirrorSpec(section=0, partner_kwargs={"input_idx": 1}, mode_signs=MODE_SIGNS)` |

A reference for the comparison already exists too —
`example_output/directional_coupler/directional_coupler_hybrid/length_sweep.npz`
(2026-08-18, `T_through`/`T_cross` on a 201 × 11 grid), the exact analogue of
the 2×2's `mm_length_sweep.npz`. So no capture job is needed; the DC can be
verified `cache_only` the same way this migration was.

**The actual blocker is a library gap, measured:** the DC threads
`_iface_mon = {"monitor_center": 0.0, "monitor_size": 6.0}` into its
**calibration** call but deliberately **not** into its EME body sweep.
`HybridSolver` has only ONE kwargs dict per section, serving both consumers:
`_calibrated_neff()` extracts `monitor_center`/`monitor_size` as a passthrough
(`hybrid.py:1525`), but the same keys also reach
`solve_section_length_sweep()`, whose `**kwargs` funnels them into
`geometry_kwargs` — and that is in the cache fingerprint. Measured on the DC's
own configuration:

| EME body length-sweep fingerprint | hash |
|---|---|
| shipped call (no `_iface_mon`) | `0c243abea80cb883` |
| migrated call (`_iface_mon` in section 1's kwargs) | `76f7367470c61021` |

**Different** — so the body would miss cache and re-solve, and bit-identity
with the 2026-08-18 reference is unreachable. Omitting `_iface_mon` instead
just moves the problem: then the *calibration* call stops matching its own
cached entry.

**The fix this needs** is a way to give the calibration source its own kwargs,
separate from the section solver's — e.g. a `calibrate_kwargs` argument
alongside `calibrate_lengths`. Filtering by signature introspection is the
tempting alternative and should be resisted: it would silently drop a
misspelled kwarg, and `_supports_*`-style introspection has already produced
one expensive failure this pass (see step 2's methodology trap). That is a
small, explicit API addition, but it is a **shared-library change affecting
every calibrated device**, so it is scoped as its own piece of work rather
than smuggled in under a single example's migration.

*(That fix was then implemented the same day — see the block above.)*

</details>

> **One harness lesson from the DC verification, worth the same shelf as step
> 2's.** The first two attempts missed cache on every lookup, and neither miss
> was a migration bug: the reconstruction took `L_mid` from the *final* 201-point
> sweep grid when the example takes it from the **200-point base grid before
> `PLOT_COUPLING_LENGTH` is spliced in** (`:195-196`) — which shifts
> `dc.coupling_length` and therefore every boundary — and then omitted
> `half_domain_aperture=True` from `_fdtd_common`, which lives in
> `solver_config`. Both were found by **diffing the built config dict against
> the cached `config.json` key-by-key**, not by guessing. That diff is the tool
> to reach for whenever a `cache_only` verification misses: it names the field
> in one step, and it distinguishes "my harness is wrong" from "the code under
> test is wrong" — which is the distinction that matters.

The `notes.md:576-582` length-budget wart does not apply to
`mmi_2x2_hybrid.py` — it carries no script-level length assertion to retire
(checked; the file has none). It presumably refers to the coupler, and travels
with that migration.

### Survey of the remaining hand-rolled cascades — and one API gap it found

Every other example still carrying a `transpose(0, 2, 1)` cascade was checked
2026-08-21 to see whether `MirrorSpec` covers it:

| example | migratable today? |
|---|---|
| `mmi_2x2_hybrid.py`, `directional_coupler_hybrid.py` | ✅ **done** |
| `directional_coupler_soi160_hybrid{,_dense,_dc_test,_length_sweep}.py` | ❌ **no — see below** |
| `mmi_4x4_90hybrid_*`, `mmi_1x4_hybrid.py` | untouched — other sessions' active work |
| `ring_resonator_hybrid_cascaded.py`, `hcg_fabry_perot_*` | N/A — loop/cavity topology, not a mirror pair |
| `capture_mmi_2x2_cascade_golden.py` | ❌ **deliberately not** — it must stay independent |

**The gap: `MirrorSpec` cannot reuse a matrix solved in a *different* call.**
The four SOI-160 coupler variants do **not** use a parity vector — they solve
**both** arms for real (`input_idx=0` and `input_idx=1`), then compose the
cross path as `M_in1ᵀ · B · M_in0`: excited on one arm, measured on the other,
with **both matrices independently simulated**. `MirrorSpec` reinterprets
`mats[j]` *from the current call* via a ±1 parity vector; it has no way to say
"section 2's matrix is the transpose of a section 0 that a different `solve()`
call produced." Running the cross excitation as its own call with
`mirror_pairs={2: 0}` gives `M_in1ᵀ · B · M_in1`, which is the wrong quantity.

This is a real limitation, not a defect: those scripts are correct as they
stand, and nothing about this plan's changes degrades them. Closing it would
mean a `MirrorSpec` variant that names a **cache key** rather than a section
index — worth considering if a third device needs it, over-engineering for two.

**One example DID need updating as a consequence of the GDS plan**, and was:
`mmi_2x2_hybrid_gds.py`'s step 7 told the reader "solver noise … a
percent-level difference means a real bug in the GDS path." The measured
answer is 1.6e-03 with a known, benign cause, so a future run would have been
read as a failure. It now prints the reference numbers and the prism-order
explanation alongside them.

---

## Manuscript dependency

Two claims in `beef_manuscript/sections/test_cases.tex` are written as though
this plan is done:

- `:36` — "The hybrid simulation for this example uses the `calibrate_body_neff`
  flag". Today the mechanism runs but the flag does not; needs step 2.
- `:34` — the two-fold mirror-symmetry reduction, ending "the entire hybrid
  response is obtained for any multimode section length from a single FDTD
  simulation". True today via the hand-rolled path; needs **step 3** to stay true
  through the supported API.

Tracked in `beef_manuscript/unimplemented_changes.md`.

**✅ BOTH RESOLVED 2026-08-21 — clear them from
`beef_manuscript/unimplemented_changes.md`.**

- `:36` is now literally true: `mmi_2x2_hybrid.py` passes
  `calibrate_body_neff=CALIBRATE_BODY` to both `HybridSolver`s and the flag is
  what drives the body's `n_eff`. Measured to move nothing (bit-identical to
  the pre-migration golden), exactly as this doc's headline finding predicted —
  it was an API question, not a physics one.
- `:34` is now true **through the supported API**, which is what step 3 was
  for. Measured: 1 `solve_section()` call across both excitations, 0
  output-side solves. The sentence no longer depends on a hand-rolled einsum
  that the library could not express.

One wording caution for the manuscript: the single-FDTD-simulation reduction
now rests on a **declared** parity vector (`MODE_SIGNS`), geometrically
precondition-checked but not derived from mode parity — see step 3(b). If the
text claims the reduction is automatic or derived, that is still ahead of the
code; "declared and validated" is what was built.
