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

## Summary

| | status |
|---|---|
| `calibrate_body_neff` flag on `HybridSolver` | ✅ implemented, tested, default `False` |
| The 2×2 propagating its body on FDTD-measured β | ✅ **already happening** — hand-rolled, not via the flag |
| `mirror_pairs` reuse within one excitation | ✅ implemented (`hybrid.py:263`, validated per-call at `:1380`) |
| Composing two excitations' matrices in one `solve()` call | ❌ **not supported** |
| A primitive that can compose them after the fact | ✅ **`beef.network.Circuit` — exists since 2026-07-29, never used for this** |
| Declaring one section as another's transpose **across** excitations | ❌ does not exist — this is what preserves the single-FDTD-run count |

Ordered work: **establish numerical equivalence (step 1) → migrate onto
`Circuit` (step 2) → generalise the mirror relation to recover the run count
(step 3) → retire the hand-rolled cascade (step 4).**

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
- These become the regression target for steps 2–4: **bit-identical**, not
  "within tolerance". Any drift is a defect, per the finding above.
- Follows the precedent of
  [`resolved/calibrate_body_neff_mode_solve_skip_validation.md`](resolved/calibrate_body_neff_mode_solve_skip_validation.md),
  which verified the mode-solve skip bit-exact against a reconstructed old path.
- Fake-solver level is enough for the composition algebra, so most of this
  belongs beside `tests/test_hybrid.py` (35 tests) and `tests/test_network.py`
  (29) and needs neither meep nor meow.

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

**Cost, ARGUED:** this route recovers saving (3) and, with `mirror_pairs`
declared, saving (1) for the *bar* path where the pairing is genuinely
mirror-symmetric. It does **not** recover (2), nor (1) for the cross path. Measure
the actual run count before and after rather than assuming.

**Done when:** the migrated 2×2 reproduces step 1's matrices bit-identically with
`calibrate_body_neff=True` and no hand-rolled einsum, and the run count is
recorded honestly even if it went up.

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

**Guard against a real trap:** whatever is built must re-validate per call, and
`hybrid.py:1370-1376` currently `continue`s on an in-memory cache hit **before**
reaching `_validate_mirror_pair` at `:1380`. A cross-excitation pair served from
`cache=` would skip its geometric check — and step 2 deliberately shares one
cache across excitations, which makes that path the common one rather than the
exception. Fix the ordering as part of this step.

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
