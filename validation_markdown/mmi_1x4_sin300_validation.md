# 1×4 MMI on 300 nm Si₃N₄ — build plan

Opened **2026-08-03**, at the user's request: a 1×4 MMI example with **FDTD, hybrid
and EME** versions, on a **300 nm silicon nitride** core.

> **This doc is a PLAN, not a gap list, and it is the first of its kind in this
> directory.** Every other file here records what has to hold before an
> *existing* device's hybrid result can be trusted. Nothing of this device
> exists yet, so this file instead orders the work required to build it and
> states, per step, what would count as evidence that the step is done. It
> converts to a normal gap list against the [four checks](README.md#the-four-checks)
> once the three scripts run.
>
> **Nothing below has been simulated.** No meep or meow solve has been run for a
> nitride stack. The measurements quoted are model-level (analytic index,
> Lorentzian fits) and were produced on the login node in this session; they are
> labelled **MEASURED**. Everything about what a nitride *device* does is
> labelled **ARGUED** and is a prediction to be checked, not a result.

---

## Summary — where the work actually is

The 1×4 geometry is nearly free. **The material model is the whole job.**

| | status |
|---|---|
| 1×4 MMI geometry | ✅ **already supported** — `MMI(n_inputs=1, n_outputs=4, …)`, no new class |
| N-output FDTD / EME / hybrid plumbing | ✅ already supported — monitors and port lists are already N-ary |
| Silicon nitride as a material | ✅ **DONE 2026-08-03 (step 1)** — `"si3n4"` in `_SELLMEIER`, from meep's `Si3N4_NIR`, tested |
| A non-Si core in `MaterialStack` | ✅ **DONE 2026-08-03 (step 2)** — `core_kind`/`clad_kind`, `SIN_300` singleton |
| 300 nm core height | ✅ **DONE 2026-08-03 (step 4)** — `SIN_300 = MaterialStack(wg_height=0.30, sz=5.0, core_kind="si3n4")`, `sz` measured-converged at the design width |
| Device dimensions for a 1×4 on nitride | ❌ unknown — must be derived, not copied from the 1×2 |

So the ordered work is: **material model → stack generalization → cache
fingerprint → cross-section design → body design → three scripts → the four
checks.** Steps 1–3 are library changes shared by every device in the repo and
are where the risk is; steps 4–5 are device design; step 6 is largely copy-work.
**Steps 1–3 are now done** (2026-08-03) — the library work is fully closed;
everything left is device-specific design and measurement.

---

## Part A — findings that set the plan

### A1. `MMI` already does 1×4 — MEASURED (read from source)

`beef/geometry/mmi.py:27` — the class is already N-input/M-output, and its own
docstring names this case:

> `n_inputs=1, n_outputs=2` gives a 1x2 splitter; `n_inputs=2, n_outputs=2` a 2x2;
> **`n_inputs=1, n_outputs=4` a 1x4**; etc.

`_port_ys(n, spacing)` centres *n* ports symmetrically about y=0 for any *n*, and
`polygons()`/`to_gds_eme()`/`section_prisms()` already accept
`output_indices=[…]` to build one branch at a time. `eme_lengths()` returns the
5-region split (`straight | taper_in | mm body | taper_out | straight`)
independent of port count.

Nothing in the geometry layer needs to change for 1×4. Confirm with a cheap
`polygons()` call and a geometry plot before trusting it (step 6a), because
"supported by signature" is not "produces a sane 4-branch layout" — the four
output tapers at `waveguide_offset` spacing must not overlap, which is a real
constraint at 1×4 that never binds at 1×2.

### A2. The geometry layer is already material-agnostic — MEASURED (read from source)

Every `section_prisms()` in `beef/geometry/` takes the core medium **as an
argument** (`section_prisms(self, x0, x1, si_medium, …)`). The parameter is
*named* `si_medium`, but nothing in the geometry inspects it. Renaming is
cosmetic; no geometry change is required to put a nitride medium in there.

`beef/geometry/base.py:16::polygons_to_prisms(vertex_lists, si_medium, wg_height=None)`
is the single funnel. Leave the names alone or rename in one pass — do not do it
halfway.

### A3. `MaterialStack` is hardwired Si-core / SiO2-clad — MEASURED (read from source)

This is the blocking finding. `beef/materials/stack.py`:

- `_SELLMEIER` (line 35) has exactly two entries, `"si"` and `"sio2"`.
- The public accessors are *named for the material*: `si_index()`, `sio2_index()`,
  `si_medium()`, `sio2_medium()`.
- `extrusion_rules_at()` (line 324) hardcodes `meow_material("si")` for the core
  layer and `meow_material("sio2")` for both cladding layers.
- `_medium()` (line 444) dispatches with a literal:
  `index_of = self.si_index if kind == "si" else self.sio2_index`.
- The class docstring is *"Si waveguide on SiO2 substrate with SiO2 upper
  cladding."*

### A4. ⚠ Aliasing the method names is NOT sufficient — ARGUED (code reading; needs a test)

The tempting minimal change is a `SiN300` stack that keeps the method names
`si_index`/`si_medium` and just returns nitride values, since the consumers reach
them by `getattr`:

- `hybrid.py:843` — `getattr(stack, "si_index" if kind == "si" else "sio2_index", None)`
- `hybrid.py:856` — `getattr(stack, "si_medium" if kind == "si" else "sio2_medium", None)`
- `hybrid.py:965` — `stack.si_index(float(np.mean(wl)))`

**That would silently simulate silicon under dispersion.** The dispersive path
does not go through the stack's accessors at all:

```
si_medium(wl) → _medium("si", wl) → meep_medium_reparam("si", wl_arr)
              → reparam_lorentz("si", …) → _SELLMEIER["si"]      ← literal lookup
```

`meep_medium_reparam()` and `reparam_lorentz()` take the **`kind` string** and
read `_SELLMEIER` directly, bypassing the instance. A stack that overrode only
`si_index` would produce the correct *frozen* index and a **silicon** Lorentzian
whenever `dispersive=True` — which is the default since 2026-08-01. The failure
is worse than a wrong number: the frozen and dispersive branches would disagree
about what material the device is made of, which is exactly the engine-asymmetry
class of bug `material_model_validation.md` step 3 exists to prevent.

**Therefore the `kind` string itself must become a stack attribute.** This is
step 2 and it is not optional.

### A5. The nitride model exists in meep and maps exactly onto this codebase's form — MEASURED

`_SELLMEIER` stores `(eps_inf, ((λ_i, σ_i), …))` for
`ε(λ) = ε∞ + Σ σ_i λ²/(λ² − λ_i²)`, lifted from meep's library with γ = 0.
Three nitride models ship with meep; **only one is usable here**:

| meep model | valid λ (µm) | γ | usable at telecom? |
|---|---|---|---|
| `Si3N4` | 0.230 – 0.830 | **1.4413** | ❌ visible-range fit, and lossy — γ≠0 is not representable in the lossless `_SELLMEIER` form |
| `Si3N4_VISNIR` | 0.207 – 1.240 | 0 | ❌ upper edge is below the C-band |
| **`Si3N4_NIR`** | **0.310 – 5.504** | **0** | ✅ **this one** |

Picking `meep.materials.Si3N4` because of the name would be a silent error: it
raises out-of-range at 1.31 and 1.55 µm (verified — `epsilon(1/1.55)` throws
`ValueError: frequency 0.7634 is below the Medium's limit, 1.2048`), and its γ≠0
means it could not be expressed here even if the range were right.

`Si3N4_NIR` has ε∞ = 1.0 and two γ=0 Lorentzians, the same convention as the
existing `"si"`/`"sio2"` entries, so the candidate entry is a direct transcription:

```python
"si3n4": (1.0, ((0.1353406,    3.0249),
                (1239.842,    40314.0))),
```

**MEASURED — this reproduces meep's own `Si3N4_NIR` to `max|Δn| = 4.441e-16`
across 1.10–2.10 µm** (the full `_MODEL_RANGE_UM`), i.e. machine precision. That
is the same standard the existing entries were held to — the module comment
claims 4e-16 for `meep.materials.Si` / `fused_quartz` over 1.4–1.7 µm — so
nitride enters on equal footing rather than as a weaker addition.

Resulting indices: **n(1.31 µm) = 2.003130**, **n(1.55 µm) = 1.996280**.

The second term (λ₂ = 1239.842 µm, σ₂ = 40314) is a far-IR pole with an enormous
strength; that is how the source parameterises it and it is numerically harmless
at telecom, but it is why the entry has two terms where `"sio2"` has three and
`"si"` has three. Do not "simplify" it.

### A6. The re-parameterised Lorentzian fits nitride BETTER than silicon — MEASURED

`dispersive=True` does not use the Sellmeier model directly for the FDTD; it
re-fits it into the single `LorentzianSusceptibility` meep can smooth
(`reparam_lorentz()`, pinned pole `lam0 = 1.15 × max(band)`). Fit quality, run
this session:

| band | kind | ε∞ | σ | max\|Δn\| |
|---|---|---|---|---|
| C 1.50–1.60 | **si3n4** | 4.009902 | 0.010002 | **1.689e-04** |
| C 1.50–1.60 | si | 12.222049 | 0.051282 | 6.231e-04 |
| C 1.50–1.60 | sio2 | 2.092887 | 0.003105 | 6.878e-05 |
| O 1.26–1.36 | **si3n4** | 4.035428 | 0.009561 | **2.172e-04** |
| O 1.26–1.36 | si | 12.474377 | 0.081810 | 1.349e-03 |
| O 1.26–1.36 | sio2 | 2.099539 | 0.002633 | 7.764e-05 |

Nitride's fit residual is **3.7× smaller than silicon's in the C-band and 6.2×
smaller in the O-band**, which is expected — its telecom dispersion is far
weaker (n moves 0.0069 over 1.31→1.55 where silicon moves ~0.010 across the
C-band alone). The pinned-pole machinery needs no new parameter for nitride.

⚠ **This does NOT mean the strip-test gate passes.** The gate is not about fit
quality; it is about whether ε∞ carries enough of the **core/cladding contrast**
for meep's subpixel averaging to still act (`si_medium()`'s docstring: putting
ε∞=1 and the whole index in susceptibilities cost **0.0496 in n_eff** at res 40).
Nitride's contrast is much weaker than silicon's — see A7 — so the gate has to
be re-run, and its result is genuinely unpredictable from the fit table above.
`_REPARAM_LAM0_FACTOR = 1.15` was validated on **four** (stack, band) pairs, all
silicon; a nitride stack adds new pairs that have never been gated.

### A7. Weak confinement is the design risk — MEASURED index, ARGUED consequences

At 1.55 µm: **n(SiN) = 1.9963**, n(SiO2) = 1.4440, n(Si) = 3.4777.

| | core–clad Δn | Δ(n²) |
|---|---|---|
| Si / SiO2 | 2.0337 | 10.0094 |
| **SiN / SiO2** | **0.5523** | **1.8999** |

Nitride's index-squared contrast is **5.3× weaker**. Consequences, all **ARGUED**
and all needing a mode solve to settle:

- **The mode is much larger and reaches further into the cladding.** `sz = 2.0`
  and the examples' `oxide_width` were chosen for silicon. They may clip the
  nitride mode. This must be measured (step 4a), not assumed either way.
- **Single-mode strip width is much wider.** The 1×2 uses `waveguide_width=0.35`
  on 160 nm Si; a 300 nm nitride strip will want something closer to ~1 µm. Copying
  0.35 would give a barely-guided mode.
- **The MMI body must be wider and longer.** Self-imaging length
  `L_π ≈ 4 n_r W_e²/(3λ)` grows with the effective width, and W_e exceeds the
  drawn width by more when contrast is low (larger Goos–Hänchen penetration).
- **Vertically, 300 nm nitride is still single-mode** — for a symmetric slab the
  TE count is `⌈2d√(n₁²−n₂²)/λ⌉`, and `2(0.3)(1.3784)/1.55 = 0.534 → 1 mode`
  (against 0.898 for 220 nm Si). Analytic, not solved: treat as a sanity figure.
- **More PML clearance and a longer decay tail.** A weakly-guided mode radiates
  more at every taper and junction, so `stop_when_dft_decayed` will take longer
  and `decay_tol` may need revisiting.

---

## Part B — the ordered steps

Each step names its **deliverable** and its **evidence** (what makes it done).
Steps 1–3 are library changes that touch every device; do them first, in order,
and do not start step 6 before step 3 is verified.

### Step 1 — nitride into the shared Sellmeier model ✅ **DONE 2026-08-03**

**Deliverable.** The `"si3n4"` entry from A5 added to `_SELLMEIER` in
`beef/materials/stack.py`, with a provenance comment naming
`meep.materials.Si3N4_NIR`, its 0.31–5.504 µm validity, and **why the
similarly-named `Si3N4` is wrong here** (A5's table — this is the kind of trap
that gets re-introduced by someone "fixing" the name later).

**Evidence.** A test in `BEEF/tests/` asserting the entry reproduces
`meep.materials.Si3N4_NIR` to ≤1e-12 across `_MODEL_RANGE_UM`. The 4.441e-16
measurement above is the expected value; the test is what keeps it true. Mirror
whatever the existing si/sio2 model test does — check for one first.

**Cost.** Minutes. No cluster time.

#### What landed

```python
"si3n4": (1.0, ((0.1353406,    3.0249),
                (1239.842,    40314.0))),
```

Provenance in the block comment above `_SELLMEIER`: K. Luke et al., *Optics
Letters* **40**, 4823–26 (2015), via meep's `Si3N4_NIR`, plus the three-model
rejection table from A5 so the `Si3N4`-by-name trap is recorded rather than
left to be rediscovered.

Tests in `tests/test_material_reparam.py`:

- `TestSellmeierMatchesMeepLibrary` — all three materials against meep's own
  library, per-material wavelength ranges (meep's `Si` refuses below 1.36 µm),
  plus a lossless/ε∞=1 invariant, nitride index anchors, an ordering check
  (n_SiO2 < n_SiN < n_Si), and a test that the two **rejected** nitride models
  really do raise at 1.31/1.55 µm — so the choice of `Si3N4_NIR` is checked,
  not merely asserted in a comment.
- `TestNitrideReparam` — nitride through the dispersive path: fit better than
  silicon's, derived pole clears its own guard, ε∞ carries the SiN/SiO2
  contrast, and the built `Medium` tracks Sellmeier across the C-band.

**580 tests pass, 18 skipped** (the `--runslow` set); the material file's own
suite went 52 → 54 → **54 passing**. The addition is inert for every existing
device — no existing test changed behaviour.

#### ⚠ Finding: the SILICON entry has a transcription truncation

Step 1 said *"mirror whatever the existing si/sio2 model test does — check for
one first."* **There was none.** The module comment had claimed since the model
was created that the entries reproduce `meep.materials.Si` / `fused_quartz` to
4e-16, and nothing tested it. Writing the test falsified it:

> `_SELLMEIER["si"]` term 2 carries `sigma = 0.003043475` — 7 significant
> figures where meep has 8 (`0.0030434748`). **MEASURED: n_Si is up to
> 8.35e-11 off `meep.materials.Si` over 1.4–1.7 µm.** SiO2 and Si3N4 both hold
> at <1e-14.

**Not corrected, deliberately — this is a decision, not an oversight.** The
error is ~7 orders of magnitude below the ~1e-3 effects every check here is
built to resolve, i.e. physically nothing. But correcting it **moves the cache
fingerprint**: measured, three entries shift in their 9th decimal —
`resolved_si_index` at 1.70 µm, and `resolved_si_reparam`'s ε∞ on **both** the
C- and O-band — which would restale every dispersive cache entry and all four
goldens repo-wide. In this project that is a scheduled event, not a side effect
of adding a material.

**Carried as an open item:** fix `0.003043475` → `0.0030434748` during the next
scheduled cache/golden regeneration, and tighten
`TestSellmeierMatchesMeepLibrary.TOL["si"]` from 1e-9 to 1e-12 in the same
commit. The test pins the current drift meanwhile, so a *real* coefficient
change still trips it. The stale module comment has been corrected in place.

### Step 2 — generalize `MaterialStack` to a core/clad pair ✅ **DONE 2026-08-03**

**Deliverable.** `MaterialStack(wg_height, sz, dispersive, core_kind="si",
clad_kind="sio2")`, with the `kind` strings threaded through **every** place that
currently hardcodes them:

- `extrusion_rules_at()` — `meow_material(self.core_kind)` / `meow_material(self.clad_kind)`
  for the three extrusion layers, and the same for the `_flat_meow_material` branch.
- `_medium()` — dispatch on `kind == self.core_kind` rather than `kind == "si"`.
- `si_index`/`si_medium`/`sio2_index`/`sio2_medium` — **keep them** as aliases for
  core/clad. Renaming them means touching `hybrid.py` (3 sites), the cache config,
  and every example; keeping them is a one-line-each delegation and the
  `getattr`-based consumers keep working unchanged. Add `core_index`/`core_medium`
  as the new preferred names and leave the old ones documented as aliases.

Defaults must preserve today's behaviour exactly: `core_kind="si"`,
`clad_kind="sio2"`.

Then add the platform singleton next to `SOI_160`/`SOI_220`:

```python
SIN_300 = MaterialStack(wg_height=0.30, sz=..., core_kind="si3n4")   # sz per step 4a
```

and export it from `beef/materials/__init__.py` and `beef/__init__.py` (both the
import list, the `__all__`, and the module docstring's materials block, which
currently reads *"SOI_160, SOI_220 — pre-built stacks for common platforms"* and
*"MaterialStack — Si/SiO2 layer definitions"*).

**Evidence.**
1. A test that `SOI_220.si_index(1.55)` and `SOI_220.extrusion_rules` are
   unchanged by the refactor — this is a shared-code change, and per this
   project's standard an "inert" refactor is an argument until a run says
   otherwise.
2. A test that `SIN_300.si_medium(band, dispersive=True)` really carries
   **nitride's** (ε∞, σ) from A6 and not silicon's. This is the A4 trap; it is
   the single most important assertion in the whole plan.
3. ⚠ **A re-run of one existing calibrated device** to confirm the refactor moved
   nothing. The directional coupler is the cheapest with a known answer.

**Cost.** Library work is hours. The confirming device re-run is one array job.

#### What landed

`core_kind="si"`/`clad_kind="sio2"` on `MaterialStack.__init__`, validated
against `_SELLMEIER` at construction time (`ValueError` immediately, not a
`KeyError` three calls deep). Threaded through:

- `extrusion_rules_at()` — `meow_material(self.core_kind)`/`meow_material(self.clad_kind)`,
  both the dispersive and `_flat_meow_material` branches.
- `_medium()` — takes `kind` from its caller (`self.core_kind` or
  `self.clad_kind`, never a literal) and calls `sellmeier_index(kind, …)`
  directly for the frozen fallback, `meep_medium_reparam(kind, …)` for the
  dispersive path.
- New preferred names `core_index`/`clad_index`/`core_medium`/`clad_medium`.
  `si_index`/`sio2_index`/`si_medium`/`sio2_medium` kept as one-line aliases —
  every existing caller, including `hybrid.py`'s three `getattr(stack,
  "si_index"/"si_medium", …)` sites, is untouched and now correctly resolves to
  whatever the stack's real core/clad is.

`SIN_300 = MaterialStack(wg_height=0.30, sz=2.0, core_kind="si3n4")` added and
exported from both `beef/materials/__init__.py` and `beef/__init__.py`
(import list, `__all__`, module docstring). **`sz=2.0` is explicitly flagged in
its own comment as inherited, not validated** — step 4a's job, not this one.

Tests: new `TestNonDefaultCoreClad` (11 tests) and
`TestExistingStacksUnchangedByGeneralization` (5 tests) in
`tests/test_material_reparam.py`. Evidence #1 and #2 both covered:

- **Evidence #1 (inertness).** `SOI_220.si_index`/`sio2_index`/`si_medium`
  (dispersive and frozen)/`extrusion_rules_at` all asserted equal to values
  computed independently via `sellmeier_index`/`reparam_lorentz`, not by
  comparison against a copy of the old code. `SIN_300` existence and shape
  checked directly.
- **Evidence #2 (the A4 trap, load-bearing).**
  `test_dispersive_medium_carries_the_real_cores_fit_not_silicons` builds a
  `MaterialStack(core_kind="si3n4", dispersive=True)`, calls the OLD alias name
  `si_medium(C_BAND)`, and asserts the returned `mp.Medium`'s (ε∞, σ) match
  `reparam_lorentz("si3n4", …)` and explicitly do **not** match
  `reparam_lorentz("si", …)`. This is a real `meep.Medium` built through the
  real dispersive path, not a read of the source.

**70/70 pass in `test_material_reparam.py`; full suite 596 passed, 18 skipped**
(up from 580 — 16 new tests, no existing test's outcome changed).

**Evidence #3 (confirming device re-run) — ✅ DONE, byte-identical.** Job
**11659362**, `directional_coupler_mode` (the coupler example that reaches
`SOI_220.extrusion_rules` directly through meow, and the cheapest of the
three — 3:33 on a warm cache as part of 11655592; this run took 3:50).
Non-destructive script (`batch/verification/run_confirm_step2_dc_mode.batch`,
no `rm -rf`, correct account) so as not to touch that device's
`*_FROZEN_pre20260801.npz` baseline.

**Measured: `supermode_sweep.npz` is byte-for-byte identical before and after
the step-2 refactor** — all 10 arrays (`wavelengths`, `n_even_wl`, `n_odd_wl`,
`Lc_wl`, `gaps`, `n_even_gap`, `n_odd_gap`, `Lc_gap`, `waveguide_width`,
`gap_nominal`), `np.array_equal` true on every one, diffed against the pre-edit
copy saved to the scratchpad before any step-2 code changed. This is a real
meow mode solve through `SOI_220.extrusion_rules_at()` — exactly the
generalized code path — run twice, once before and once after, with identical
output. The refactor is confirmed inert for the default Si/SiO2 stacks, not
just argued to be.

#### ⚠ Step 2 sharpens, but does not close, the step-3 cache-fingerprint risk

Reading `resolved_material_config()` closely while wiring `_medium()` turned the
step-3 risk from a code-reading argument into something concrete enough to pin
in a test now, even though the fix itself is still step 3's job:

```python
for kind in ("si", "sio2"):                 # ← still a literal, unchanged
    eps_inf, sigma, _ = reparam_lorentz(kind, wl)
    config[f"resolved_{kind}_reparam"] = [...]
```

After step 2, `core_medium()`/`si_medium()` on a nitride stack correctly builds
**nitride's** Lorentzian. But `resolved_material_config()`'s dispersive-branch
loop is still hardcoded to `("si", "sio2")`, so for `dispersive=True` it fits
and records **silicon's** coefficients into `resolved_si_reparam` regardless of
what the stack's `core_kind` is. The medium actually simulated and the medium
recorded in the cache key are now a genuine, demonstrated mismatch — not merely
different names for the same thing.

**Consequence, stated plainly: do not run a non-default-core stack under
`dispersive=True` through the caching solvers (`FDTDSolver`/`EMESolver`/
`HybridSolver`) until step 3 lands.** The frozen (`dispersive=False`) path is
fine — `resolved_si_index`/`resolved_sio2_index` already read the real
core/clad index via the alias chain. Pinned as a live check,
`test_resolved_material_config_still_has_the_step3_gap`, deliberately written
to assert **today's wrong behaviour** and marked to start failing — on purpose
— once step 3's loop is fixed. That failure is step 3's own regression signal;
it is not a bug in this test.

### Step 3 — cache fingerprint ✅ **DONE 2026-08-03, item 2 replaced — see below**

**Deliverable (as originally planned).** Two hardcoded spots fixed:

1. `resolved_material_config()` (`stack.py:462`) hardcodes the key *names* and
   the fit *kinds*:
   ```python
   "resolved_si_index":   [... material_stack.si_index(w) ...],
   "resolved_sio2_index": [... material_stack.sio2_index(w) ...],
   ...
   for kind in ("si", "sio2"):
       eps_inf, sigma, _ = reparam_lorentz(kind, wl)
   ```
   That loop must iterate `(stack.core_kind, stack.clad_kind)`. **As written it
   would fingerprint silicon's Lorentzian coefficients for a nitride run** —
   ARGUED from reading, and it should be turned into a failing test before it is
   fixed.
2. `beef/cache/config.py:119` records a stack as `{"wg_height", "sz"}` only.
   `core_kind`/`clad_kind` must join it, or a nitride stack and an SOI stack at
   the same `wg_height` would share cache keys through that field.

**Evidence.** A test that two stacks differing **only** in `core_kind` produce
different fingerprints, both in `resolved_material_config()` and in
`_solver_config()`. This is the exact defect class that
`resolved_material_config`'s own docstring was written for — a per-solver
override that never reached the key and served old physics from a byte-identical
entry. Do not hand-wave it.

⚠ Note whether adding a key to these dicts invalidates **every existing cache
entry and all four goldens** repo-wide. `resolved_material_config` only adds
`resolved_dispersive`/`resolved_*_reparam` when the band is broadband, so
renaming `resolved_si_index` → `resolved_core_index` would be a global break.
**Prefer keeping the existing key names** and adding `core_kind`/`clad_kind`
alongside, so SOI runs keep their fingerprints and only genuinely-new nitride
runs get new ones. Confirm by fingerprinting an SOI config before and after.

#### What landed

**Item 1 — done exactly as planned, and measured INERT for every existing
stack.** `resolved_material_config()`'s reparam loop now reads:

```python
for label, kind in (("si", material_stack.core_kind),
                    ("sio2", material_stack.clad_kind)):
    eps_inf, sigma, _ = reparam_lorentz(kind, wl)
    config[f"resolved_{label}_reparam"] = [...]
```

`label` stays the fixed `("si", "sio2")` pair — same reasoning as
`resolved_si_index`/`resolved_sio2_index` since step 2: these are now fixed
**slot names** (core/clad), not literal material names. `kind` is what
actually gets fit. For a default stack `kind == label` always, so both the key
names AND the computed values are byte-identical to before — **measured**,
not argued: `resolved_material_config(SOI_220, C_BAND)` and
`resolved_material_config(SOI_160, O_BAND)` reproduce step 1's own A6 table
exactly (`resolved_si_reparam` eps_inf `12.222048631`/`12.474376861`). Zero
cache/golden invalidation from this half of the fix.

For a nitride stack the same call now correctly returns nitride's fit —
confirmed by flipping the test that had pinned the bug
(`test_resolved_material_config_still_has_the_step3_gap` →
`test_resolved_material_config_now_records_the_real_cores_reparam_fit`, plus a
symmetric clad-side test and an explicit default-stack inertness test, all in
`tests/test_material_reparam.py`).

**Item 2 — NOT done as originally planned, replaced by a better-scoped fix
found while implementing.** Investigating exactly how `material_stack` reaches
`beef/cache/config.py`'s fingerprint turned up two things the original plan
text didn't have:

1. **The `{wg_height, sz}` gap named in the plan doesn't actually collide
   today.** `_solver_config()`'s bare dict is only ONE field in a much larger
   fingerprint; `FDTDSolver.solve`/`solve_section` and `EMESolver.solve_section`/
   `solve_section_length_sweep` ALL also merge `resolved_material_config(...)`
   into the SAME outer dict, and that half is core-aware since step 2's
   `si_index`/`sio2_index` aliases. Two stacks differing only in `core_kind`
   already get different `resolved_si_index` values, hence different overall
   fingerprints, for these four methods — the plan's stated collision doesn't
   reach them.
2. **`EMESolver.solve()` and `EMESolver.solve_length_sweep()` call
   `resolved_material_config()` NOT AT ALL** — a gap the original plan missed
   entirely because it was written from `resolved_material_config`'s call
   sites in isolation, not by checking which cache-relevant methods lack a
   call site. For these two methods, `_solver_config()`'s bare `{wg_height,
   sz}` genuinely is the WHOLE material fingerprint. Worse than the planned
   item 2 suggested: it's not just `core_kind` that's invisible there, so is
   `dispersive` — and `dispersive` demonstrably changes what these methods
   simulate (`EMESolver._extrusion_rules()`'s frozen-vs-exact branch). This
   predates nitride entirely; it could in principle already have bitten a
   pure-Si/SiO2 `EMESolver(dispersive=False)` A/B run cached alongside a
   dispersive one.

Given that, adding `core_kind`/`clad_kind` to `_solver_config()`'s
`MaterialStack` branch (the original item 2) would have (a) invalidated
**every cache entry and every golden repo-wide** — the SAME magnitude of break
as the 2026-08-01 dispersion-default flip, since `material_stack` is in the
`_METHOD_ATTRS` allowlist for every cached method — while (b) still not
fixing `dispersive`-blindness for `EMESolver.solve()`/`solve_length_sweep()`,
since `core_kind`/`clad_kind` alone doesn't cover that flag.

**Replacement fix, scoped to the two methods that actually lack material
fingerprinting:** `EMESolver.solve()` and `EMESolver.solve_length_sweep()`
(`beef/solvers/eme.py`) now merge `**resolved_material_config(self.material_stack,
wavelengths, self._resolved_dispersive())` into their cache configs, the exact
pattern `solve_section()`/`solve_section_length_sweep()` already used. This:

- Closes the `core_kind` gap (item 1's fix propagates here for free).
- **Also closes the pre-existing `dispersive` gap** for these two methods,
  which the original item 2 would not have.
- Invalidates cache entries for `EMESolver.solve()`/`solve_length_sweep()`
  callers ONLY — i.e. `mmi_1x2_eme.py`, `mmi_2x2_eme.py`, and the future
  `mmi_1x4_eme.py` — not the hybrid or FDTD examples, which were never at
  risk (item 1 alone, and it's inert for their default stacks).
- `_solver_config()`'s `MaterialStack` branch is **untouched** — still
  `{"wg_height", "sz"}` — deliberately, now that every cache-relevant method
  carries material fingerprinting through `resolved_material_config()` instead.

**Evidence.**
- Fast, meep/meow-free: `TestEMESolveFingerprintsMaterial` in
  `tests/test_cache_config.py` — confirms via `inspect.getsource()` that both
  methods' source actually calls `resolved_material_config(...)` with
  `self._resolved_dispersive()` (not the stack's own flag, matching
  `solve_section`'s established pattern), plus a check that `solve_sweep()`
  (the uncached, memory-batched variant) has no `save_dir` and is correctly
  out of scope.
- Real, `--runslow`-marked, in `tests/test_eme_length_sweep.py::
  TestSolveFingerprintsMaterial`: two ACTUAL `EMESolver.solve(save_dir=...)`
  calls, differing only in `core_kind` at identical `wg_height`/`sz`, produce
  **two distinct cache folders** (measured: `len(run_dirs) == 2`) — the exact
  collision this closes, confirmed end to end rather than argued from source.
  A companion test confirms the fix didn't overshoot into cache-missing on
  every call (two identical solves still hit the same folder).
- Regression: **601 passed, 20 skipped** in the fast suite (up from 598 — 3 new
  fast tests, `test_eme_length_sweep.py`'s 2 new slow tests correctly
  skipped without `--runslow`), plus the full `--runslow` `test_eme_length_sweep.py`
  suite (11/11) and the targeted cache/material/hybrid suites, all green.
- Confirmed **zero** effect on `_solver_config()` or on any solve_section-family
  cache key: those were never touched.

### Step 4 — cross-section design (cheap mode solves, no FDTD) ✅ **DONE 2026-08-03**

**4a. Vertical span.** Sweep `sz` (2.0, 3.0, 4.0 µm) for a 300 nm nitride strip
and find where the fundamental n_eff stops moving. Deliverable: the `sz` that
goes into `SIN_300`. Evidence: an n_eff-vs-`sz` table with a stated tolerance.
Do **not** inherit `sz=2.0` silently — A7 says it is the most likely thing to be
wrong, and this is a minutes-long mode solve, not a simulation.

**4b. Single-mode width.** Sweep strip width at 300 nm nitride and find the
largest width that is still single-mode at the chosen band (count guided modes,
not just n_eff). Deliverable: `waveguide_width`. Evidence: a mode-count-vs-width
table.

**4c. Band.** Pick C-band (1.50–1.60) or O-band (1.26–1.36) explicitly and say
why. The 1×2 uses O-band with SOI_160; the coupler and 2×2 use C-band. Nitride's
`Si3N4_NIR` model covers both. **Recommendation: C-band**, to reuse the
band-specific reparam validation the C-band already has and to sit in the middle
of the nitride model's range — but state the choice in the scripts, do not let it
be implicit.

#### 4c is settled — no mode solve needed

**Decision: C-band (1.50–1.60 µm), λ=1.55 µm centre.** Three reasons, all
already on the record rather than newly argued:

1. **Nitride's own reparam fit is measurably better in C-band** — step 1's A6
   table: max|Δn| = 1.689e-4 (C-band) vs 2.172e-4 (O-band), a real property of
   *this* material, not just infrastructure reuse.
2. Comfortably inside `Si3N4_NIR`'s 0.310–5.504 µm validity with margin on both
   sides — no extrapolation risk.
3. Matches the directional coupler and 2×2 MMI, the two devices with the most
   complete four-checks records, so a future cross-device comparison needs no
   band translation.

Used as the mode-solve wavelength for both 4a and 4b below — stated here, not
left implicit in the script, per this step's own instruction.

#### 4a/4b — script written, hit a real cross-engine bug, fixed, real sweep in flight

Combined into one script (both sweeps need the same cross-section builder and
guided-mode counter, and the plan's own note that these are two 1D sweeps, not
a joint optimisation, makes one file natural):
`examples/verification/python/mmi_1x4_sin_cross_section.py`.

**4a** sweeps `sz` at a fixed reference width. That width is **not a
result** — a scaling-argument seed, stated so it's checkable rather than
silently assumed: matching normalised frequency `V ~ width·√(Δn²)` to
silicon's known single-mode point (450 nm at Δn²=10.0093, A7) under nitride's
Δn²=1.8999 gives `width ≈ 0.45·√(10.0093/1.8999) ≈ 1.03 µm` → seeded at
1.0 µm. Convergence criterion: successive `|Δn_eff| < 1e-4`.

**4b** sweeps width at the sz 4a converges to, counting GUIDED modes
(`Re(n_eff) > n_clad`, the same convention `mode_spectrum_convergence.py` and
`HybridSolver`'s mode-correspondence filter use — README's Finding 1 is why
this matters: every interface basis measured in this repo so far carries at
least one unguided returned mode, so a raw returned-mode count would
over-count). Reports the widest single-mode width.

A confirmatory pass re-runs 4a's sz sweep at the width 4b actually lands on
(not just the 1.0 µm seed), closing the loop between the two 1D sweeps without
a full joint optimisation.

Three mechanical bugs were caught and fixed by smoke-testing before spending
real cluster time: a mesh-construction pattern copied from
`directional_coupler_mode.py` (`linspace(-L/2, L/2, int(res*L)+1)` per axis)
only gives equal x/y pixel spacing when `res*L` lands on an integer for both
axes — true at the fixed dimensions that script uses, not guaranteed for a
*swept* `sz`; fixed with a shared `_axis_grid()` helper building both axes
from the same explicit `1/res` step. Two index bugs in the convergence-report
printing (wrong bracket in the "CONVERGED at sz=…" message; `if converged_idx`
used as a truthy check, wrong when convergence happens on the very first step
since `0` is falsy in Python) were fixed the same pass.

#### ⚠ MEASURED 2026-08-03 — `compute_modes_meep` silently returns the WRONG answer for this structure

The first real submission, **job 11659619**, ran to its 4b phase and reported
**zero guided modes at every width tested (up to 2.4 µm) and every sz tested
(up to 5 µm)**. Cancelled once it was clear this wasn't going to produce a
useful result (the same "0 guided" outcome repeating across widths already
sampled) — this section records what was actually wrong rather than what the
job printed.

**This was not real physics.** An analytic check — the exact symmetric
3-layer slab dispersion relation, infinite-width limit, solved numerically —
predicts a genuine, **tightly-confined** bound TE0 mode for 300 nm Si₃N₄ in
SiO₂ at 1.55 µm: **n_eff = 1.673911**, cladding evanescent decay length
**0.291 µm**. That decay length is comfortably inside every box tested (2–5 µm
half-spans); a finite-width channel confines somewhat less than the infinite
slab but nothing about this margin suggested near-cutoff behaviour.

Direct inspection ruled out the geometry as the cause: the built structure's
index map was verified pixel-by-pixel exact (core `n=1.996280` spanning
precisely ±0.15 µm × ±0.5 µm, cladding `n=1.444024` everywhere else, at
`width=1.0, sz=2.0`). A **field plot of the returned "fundamental" mode**
(`example_output/verification/mmi_1x4_sin_cross_section/fields/cross_section_w1.0_sz2.0.png`)
shows why: the intensity is smooth and nearly uniform across the *entire*
8 µm domain, barely perturbed by the tiny high-index core — a spurious
**box/continuum mode**, not a real guided field. Widening the strip
(1.6 µm, 2.4 µm) at the same `sz` returned the **identical** n_eff regardless
of width — impossible for a real guided mode (which must depend on core
width) and conclusive proof the eigensolver wasn't tracking the waveguide at
all.

**Root cause, isolated by elimination, not guessed:**

1. **Not the geometry or box size** — see above.
2. **Not a general solver/box bug** — a control case (identical box, `core_kind="si"`
   instead of `"si3n4"`) found **3 clean guided modes immediately**
   (n_eff 2.50–2.94) via the exact same `compute_modes_meep` call.
3. **Not "not enough requested modes"** — a 25-mode request was in flight when
   this line of investigation overtook it (cancelled once the real cause was
   found; see next point for why more modes wouldn't have helped anyway).
4. **Not a bad initial k-guess.** `compute_modes_meep` (`meow_MEEP/src/meow/
   fde/meep.py`) finds eigenmodes via meep's own band-structure machinery —
   `sim.get_eigenmode(direction=NO_DIRECTION, kpoint=Vector3(z=1), …)`, a
   Newton search on the propagation constant starting from a **hardcoded**
   guess (`kpoint=1` → guessed n_eff = 1×wl = 1.55, coincidentally close to
   the true 1.674). Explicitly forcing three different guesses directly into
   meep's `get_eigenmode` on the identical structure — `kpoint.z ∈ {1.0, 1.65,
   1.9}`, spanning from the hardcoded default to essentially the core index —
   **converged to the exact same wrong value (n_eff=1.440790) every time.**
   That rules out the guess's accuracy as the variable and points at
   something more fundamental in how that Newton search behaves for
   weakly-index-contrasted structures — not something to patch blind.

**The fix: switch mode-solver backends, not patch the search.** meow ships a
second backend, `compute_modes_tidy3d` — a proper finite-difference
eigenmode (FDE) solver rather than a k-vs-frequency band search, and already
a first-class option elsewhere in this codebase (`EMESolver(mode_solver=
"tidy3d")`, exercised in `tests/test_eme_length_sweep.py`). On the identical
failing cross-section it **immediately** found 2 guided modes
(n_eff = 1.5727, 1.5109) — consistent with the analytic estimate (a finite
1.0 µm-wide channel confining somewhat less than the infinite-slab limit, as
expected). Cross-validated against `compute_modes_meep` on the Si control
*first*, so as not to trade a known-good result for an unvalidated one: the
two backends agree to **~1–3%** there (2.938/2.630/2.501 vs 2.953/2.651/2.535),
in line with this repo's *other* documented engine-agreement gaps (MPB-vs-meep
phase ~2.3%, MMI hybrid-vs-FDTD ~0.5–1%) — not a red flag, and clear evidence
`compute_modes_tidy3d` is a genuine working alternative, not a different or
broken calculation.

**Why not patch `compute_modes_meep` itself**, despite the user's explicit
go-ahead to modify `meow_MEEP`: its search behaviour is shared by every
existing (silicon) device in this repo, which already depend on its current
behaviour working — and "guess-independent convergence to the same wrong
root" is a symptom I don't have the MPB/meep internals visibility to safely
re-derive a fix for without real risk of silently changing results for
devices that currently work. Switching backends for this one script achieves
the goal (correct nitride modes) with no blast radius on anything else.
**This is a real, generally-applicable limitation of `compute_modes_meep`
worth flagging for any future low-index-contrast material** (returns a wrong
answer with **no error or warning** — the dangerous kind of bug) — recorded
here rather than silently worked around, and worth a defensive check inside
`meow_MEEP` itself at some point (e.g. warn when varying the initial k-guess
doesn't change the result, or when the top returned mode sits suspiciously
close to the cladding index) rather than leaving the next low-contrast device
to rediscover this from scratch.

`mmi_1x4_sin_cross_section.py` now imports `compute_modes_tidy3d` throughout
(no `compute_modes_meep` calls remain). Its module docstring carries the full
chain of evidence above.

**Candidate ranges also revised**, from real (not assumed) numbers measured
with the fixed backend on a reduced smoke test:
- `SZ_CANDIDATES`: `2.0→2.5` measured `|Δn_eff|=1.84e-3`, an order of
  magnitude above the `1e-4` tolerance — extended from `{2.0,…,5.0}` to
  `{2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0}` µm.
- `WIDTH_CANDIDATES`: `width=0.6` µm already measured **2 guided modes**
  (n_eff 1.498/1.477) — the original `{0.6,…,2.4}` range never brackets a
  single-mode point at all. Extended down to `{0.2, 0.3, 0.4, 0.5, 0.6, 0.8,
  1.0, 1.4, 1.8, 2.2}` µm.
- `MAX_MODES_REQUESTED`: `2.4` µm measured 5 of 6 requested modes guided —
  bumped `6 → 8` so the widest candidates have headroom.

**Job 11660092** submitted (`batch/verification/run_mmi_1x4_sin_cross_section.batch`,
account `gts-sr108-ece`, single-node/8-cpu/serial, `-t 240`; supersedes
11659619, cancelled, and its diagnostic follow-up 11659878, cancelled once
the root cause was isolated by direct testing rather than needing its
remaining test cases). **Completed.**

#### 4a/4b results — MEASURED, and a real design tension surfaced

**4a, at the reference width (1.0 µm): converged cleanly at `sz = 4.0 µm`.**
A textbook convergence curve — successive `|Δn_eff|` shrinking monotonically
from `2.29e-3` (sz 3→4) to `1.45e-4` (4→5) to `9.93e-6` (5→6) down to
`2.87e-10` by sz 8→10. Fundamental n_eff plateaus at **1.575107**.

**4b, at sz=4.0 µm: a clean, monotonic guided-mode count vs. width.**

| width (µm) | guided modes | n_eff (highest → lowest) |
|---|---|---|
| 0.20 | 1 | 1.4447 |
| 0.30 | 1 | 1.4488 |
| 0.40 | 2 | 1.4569, 1.4552 |
| 0.50 – 1.00 | 2 | (1.4635–1.5751 across this range) |
| 1.40 | 3 | 1.6136, 1.5246, 1.4632 |
| 1.80 | 4 | 1.6337, 1.5366, 1.5207, 1.4572 |
| 2.20 | 5 | 1.6453, 1.5620, 1.5442, 1.4804, 1.4465 |

**Single-mode width bracketed to [0.30, 0.40] µm** — far narrower than the
1.03 µm scaling-argument seed (A6/4a), which sits solidly in 2-mode
territory. Confirms A7's prediction that this platform confines much more
weakly than Si/SiO2, just far more sharply than the seed estimate suggested.

**⚠ Confirmatory check (sz re-swept at width=0.30 µm, the actual single-mode
point) did NOT converge to the 1e-4 tolerance, even out to sz=10 µm:**

```
sz:     2.0      3.0      4.0      5.0      6.0      7.0      8.0      10.0
n_eff:  1.4544   1.4505   1.4488   1.4479   1.4474   1.4471   1.4469   1.4472
|Δn|:      3.88e-3  1.68e-3  8.87e-4  5.14e-4  3.10e-4  1.90e-4  (3.27e-4, 2-unit step)
```

**Why:** the single-mode window is weakly confined by construction — margin
above `n_clad` (1.444024) is only **~0.003–0.004** at width=0.30 µm (heading
toward ~1.4465–1.447 as sz→∞), and **~0.0007** at width=0.20 µm, i.e. barely
guided at all. A weakly-confined mode's evanescent tail decays slowly, so it
keeps "seeing" the domain boundary well past the box sizes that converge the
2-mode reference case instantly. **`sz=4.0 µm`, chosen from the wrong
reference width, is NOT validated at the true single-mode design point** —
exactly the failure mode this confirmatory step exists to catch, and it
caught it.

**This is a genuine property of 300 nm Si₃N₄, not a numerical artifact**: at
this thickness the single-mode width window is narrow and inherently weakly
confined. Real SiN platforms at comparable thickness often lean into
wider/multimode designs for exactly this reason.

**User's call, given this finding:** reconsider the design width — try
**1.1 µm** instead of chasing the narrow single-mode window. Measured
directly (not assumed) in a follow-up job rather than reasoned about from the
1.0/1.4 µm bracketing points already on hand — see below.

#### Follow-up — width = 1.1 µm ✅ MEASURED, clean result

**Job 11660166** (`examples/verification/python/mmi_1x4_sin_width1p1_check.py`,
`batch/verification/run_mmi_1x4_sin_width1p1_check.batch`).

```
sz (µm)   mode 1 n_eff (margin)   mode 2 n_eff (margin)
2.0       1.585233 (+0.14121)     1.516602 (+0.07258)
3.0       1.587228 (+0.14320)     1.511519 (+0.06750)
4.0       1.587341 (+0.14332)     1.510750 (+0.06673)
5.0       1.587347 (+0.14332)     1.510636 (+0.06661)
6.0       1.587348 (+0.14332)     1.510619 (+0.06660)
8.0       1.587348 (+0.14332)     1.510616 (+0.06659)
10.0      1.587348 (+0.14332)     1.510616 (+0.06659)
```

**Converges cleanly at `sz = 5.0 µm`** (`|Δn_eff| < 1e-4` going 4.0→5.0; by
sz=8 the fundamental has settled to 6 decimal places). A completely different
regime from the 0.3 µm case: **2 guided modes, both comfortably confined**
(margins +0.143 and +0.067 above `n_clad`, vs. ~0.003 and ~0.0007 at the
narrow single-mode width) and converging in a handful of µm rather than
never converging within 10.

**`SIN_300.sz` updated 2.0 → 5.0 µm** in `beef/materials/stack.py`, with this
measurement recorded in place of the old placeholder comment. Full suite
re-run clean (72/72 in `test_material_reparam.py`; `SIN_300`'s own test
doesn't pin a specific `sz` value, so it needed no change).

**⚠ Design-scope change, flagged for step 5/6:** 1.1 µm is deliberately a
**2-mode width, not single-mode** — the user's call, made after measuring
that the true single-mode window (0.2–0.35 µm) is impractically weakly
confined for this 300 nm platform. This changes an assumption every
downstream step was written against:

- **`mmi_1x2_eme.py`/`mmi_1x2_hybrid.py`-style per-branch EME solves assume a
  single-mode access waveguide** (`n_modes_in`/`n_modes_out` sized for one
  mode per port). At 1.1 µm there are two real guided modes at the access
  waveguides; step 6's hybrid/EME scripts need to either request 2 modes at
  every port genuinely, or make an explicit, stated decision to launch only
  the fundamental and treat the second mode as unwanted higher-order content
  (common practice, but a decision to write down, not assume).
- **Step 5's MMI body design** (next section) needs the taper from a 1.1 µm
  access width up to `mm_width` to be built and solved with this in mind —
  the taper itself may act as a mode filter/converter depending on its
  length, which is exactly the kind of behaviour `HybridSolver`'s
  `check_mode_correspondence` check exists to catch.
- **Check 2 (mode identity, README's four checks)** becomes more load-bearing
  here than on any existing device: a wider, genuinely 2-mode access
  waveguide is closer to the modal multiplexer's territory (13+ guided
  modes, where check 2 already failed once) than to the 1×2/2×2 MMIs'
  single-mode access waveguides (where it also failed once, at N_MODES=12,
  before being fixed). Run it early once real devices exist, not as an
  afterthought.

None of this blocks step 5 — the numbers needed (`sz=5.0`, `waveguide_width=
1.1`) are now in hand — but it changes what step 5/6 have to solve for, and
is recorded here so it isn't silently assumed away.

### Step 5 — MMI body design (the one genuinely new physics step)

**Deliverable.** `mm_width`, `mm_length`, `taper_width`, `taper_length`,
`waveguide_offset` for a 4-way split.

**Method — do not use a formula alone.** Use the same trick the coupler uses for
L₅₀: `EMESolver.solve_section_length_sweep()` computes the body's modes **once
per wavelength** and evaluates every length analytically, so an `mm_length` sweep
of 200 points costs one mode solve. The MMI body is a uniform cross-section of
width `mm_width`, which is exactly what that method requires
(`eme.py:537`; its docstring: *"[x0, x1] must be a uniform-cross-section
region"*). So:

1. Fix a candidate `mm_width` wide enough to support the modes a 4-way split
   needs and to fit four output tapers at `waveguide_offset` spacing without
   overlap — that geometric constraint is what makes 1×4 different from 1×2 and
   it should be written down as an inequality, then asserted in the script.
2. Sweep `mm_length` and read off the length where the four outputs are equal.
   Use `4 n_r W_e²/(3λ) · (3/N)` only as the **starting bracket** for the sweep.
3. Iterate `mm_width` if the equal-split point is unreachable or the excess loss
   at it is poor.

**Evidence.** An `mm_length` sweep plot with the chosen length marked, and the
four output powers at it. This is the 1×4 analogue of the coupler's
`coupling_length_sweep.png`.

⚠ **Pin the chosen length exactly.** The coupler shipped a bug for weeks where a
pinned length was `argmin`'d onto a `linspace` grid that did not contain it, so
the "same" length differed by 0.0208 µm between two examples (fixed 2026-08-03,
job 11655592). If the 1×4 hybrid pins an `mm_length` for its plots, splice it
into the sweep grid and assert the hit is exact — do not repeat that.

### Step 6 — the three example scripts

Directory `BEEF/examples/mmi_1x4/python/`, mirroring `mmi_1x2/`:

**6a. Geometry first, alone.** Build the `MMI(n_inputs=1, n_outputs=4, …)`,
print `total_length`/`output_ys`, and run `plot_*_geometry`. Confirms A1 and
catches taper overlap before any solve. Minutes.

**6b. `mmi_1x4_eme.py`.** Model on `mmi_1x2_eme.py`. Note its structure: the EME
solver measures **one output branch per solve** (`output_idx=…`, restricting the
cross-sections to that branch so the mode is the isolated waveguide mode rather
than a supermode), then merges the results into one `SMatrix` via
`smat._data.update(...)`. **1×4 means four solves, not two** — and the
"comfortably decoupled" justification in the 1×2's comment
(`waveguide_offset - waveguide_width = 0.65 µm` edge-to-edge) must be
**recomputed for nitride**, where a weakly-confined mode couples across a gap
that would be safe in silicon. If the four branches are not decoupled at the
chosen offset, that invalidates the per-branch method, not just the number.

**6c. `mmi_1x4_fdtd.py`.** Model on `mmi_1x2_fdtd.py`. `output_monitor_centers`
and `output_monitor_widths` are already lists — supply four of each. The 1×2 uses
half-domain windows (`[2.5/2, -2.5/2]`, widths `[2.5, 2.5]`); **four monitors
cannot each take half the domain**, so the windows must be sized to the actual
port pitch without overlapping. This is the one place the 1×4 FDTD script is not
a mechanical copy of the 1×2.

**6d. `mmi_1x4_hybrid.py`.** Model on `mmi_1x2_hybrid.py` — FDTD sections on the
tapered input/output, EME on the uniform body, per-output solves as in 6b.
Decide `N_MODES` and `RESOLUTION` **jointly**, not by copying: README check 3 is
explicit that a shared resolution knob is how a 7% error survived for months, and
the 1×2's check 2 actually *failed* at `N_MODES=12` until it was dropped to 8. A
1×4 body is wider and supports more modes than either existing MMI, so neither
device's numbers transfer.

### Step 7 — wiring

- `BEEF/examples/Makefile`: add `MMI_1X4_EXAMPLES = mmi_1x4_eme mmi_1x4_fdtd
  mmi_1x4_hybrid`, add to `EXAMPLES` and `DEVICES`, add the `mmi_1x4` aggregate
  target and the three per-script targets.
- `batch/run_mmi_1x4_examples.batch`: copy an existing device runner. **Use
  `--account=gts-sr108-ece`** (the balance account; most existing scripts still
  hardcode `gts-sr108-gedc_main`). Index order in the `EXAMPLES` array must match
  the Makefile — `SLURM_ARRAY_TASK_ID` selects positionally.
- ⚠ **Do not copy the `rm -rf ../../example_output/<device>/$NAME` line** without
  thinking. It is why `run_directional_coupler_L151.batch` had to be written from
  scratch on 2026-08-03: the stock runner would have deleted
  `*_FROZEN_pre20260801.*` baselines that cannot be regenerated. A brand-new
  device has nothing to protect yet, so the line is safe *here* — but size the
  walltime off a real cold run, not off the stock 300 min, which has already
  timed out one coupler job at 05:00:26.
- `beef/__init__.py` docstring + `__all__` for `SIN_300` (step 2).

### Step 8 — the four checks

Only now does this device get a normal gap list. From
[README.md](README.md#the-four-checks), in the order that makes each one
cheap:

| check | 1×4-specific note |
|---|---|
| **2. Mode identity** | Do this **first** — it is mode solves only, no time-stepping, and it is the check that actually failed on the 1×2. A wider nitride body supports more near-degenerate modes than either existing MMI, so this is *more* likely to fail here, not less. |
| **3. Resolution convergence, per engine** | Sweep FDTD and EME resolution **independently**. Sweep a Δβ-like quantity, not just output power. Expect nitride to need different resolution than silicon — lower contrast means the index step is easier to resolve, but the mode is larger, so the two effects push opposite ways and it must be measured. |
| **1. Engine agreement on the shared cross-section** | `verify_item_6_engine_agreement.py` on the MMI body. Both MMIs measured a ~0.5% engine bias (0.9946 / 0.9953); whether that is material-independent is an open and interesting question this device can answer. |
| **4. Full-device FDTD reference** | `mmi_1x4_fdtd.py` vs `mmi_1x4_hybrid.py` at the identical `mm_length`. Record provenance in the doc — config, solver fingerprint, geometry — per the coupler's item 5 lesson. |

Plus the **material-model gate** inherited from step 1–3, which is not one of the
four and cannot be skipped:

- **Strip test / subpixel-averaging gate for nitride** (A6). Re-run whatever
  `material_model_validation.md`'s gate script is for the new (stack, band)
  pairs, against its 0.002 tolerance. Nitride's weak contrast is the reason this
  cannot be inherited from silicon's passes.

---

## Measured vs argued — read this before quoting anything above

**MEASURED this session** (login node, `jmh` env, model-level only — no device
simulated):

- meep ships `Si3N4`, `Si3N4_NIR`, `Si3N4_VISNIR`; their validity ranges and γ
  values as tabulated in A5; `Si3N4` and `Si3N4_VISNIR` both raise out-of-range
  at 1.31 and 1.55 µm.
- The `"si3n4"` entry in A5 reproduces `Si3N4_NIR` to **4.441e-16** over
  1.10–2.10 µm. n(1.31) = 2.003130, n(1.55) = 1.996280.
- The `reparam_lorentz` table in A6, all six rows.
- Index contrast figures in A7 (from `sellmeier_index`).

**MEASURED (read from source, not run):** A1, A2, A3 — the `MMI` signature and
docstring, the geometry layer's medium-as-argument pattern, and the four
hardcoded-material sites in `MaterialStack`.

**MEASURED as of 2026-08-03 (were ARGUED when this doc opened):** A4's claim
that method-name aliasing would silently simulate silicon — demonstrated with
a real `meep.Medium` before the fix, then flipped to confirm the fix, in
`tests/test_material_reparam.py::TestNonDefaultCoreClad`. The step-3 claim that
`resolved_material_config` would fingerprint silicon's coefficients for a
nitride run — same treatment, in `TestExistingStacksUnchangedByGeneralization`.
A third, related claim found (not predicted) while implementing step 3:
`EMESolver.solve()`/`solve_length_sweep()` had NO material fingerprinting at
all, not even the wrong-kind version — confirmed with a real, `.solve()`-
executing test in `test_eme_length_sweep.py::TestSolveFingerprintsMaterial`
before being fixed.

**MEASURED 2026-08-03 (step 4 — the first real device-level simulation in
this whole plan, not model arithmetic):** `sz=5.0 µm` (converged, jobs
11660092 → 11660166); `waveguide_width=1.1 µm` is a genuinely 2-mode design
point with comfortable confinement (margins +0.143/+0.067 above `n_clad`);
the *true* single-mode window (0.2–0.35 µm) exists but is impractically
weakly confined (margin ~0.003 or less, doesn't converge even at sz=10 µm).
Band choice (C-band) was ARGUED from step 1's numbers, not solved for, and
stays that way — no mode solve was needed for it. Also measured: a real,
general bug in `compute_modes_meep` (silently wrong for weakly-confined
structures — see step 4's writeup for the full elimination chain), fixed by
switching to `compute_modes_tidy3d`.

**ARGUED — predictions, not results:** the rest of A7's consequence list not
resolved by the above (decoupling of the four MMI output branches at a given
`waveguide_offset`, taper behaviour); `mm_width`/`mm_length`; every predicted
check outcome (steps 6–8); and the design implication itself — that a 2-mode
access waveguide needs `N_MODES≥2` stated through the body/taper design and
step 6's scripts — is reasoned from the step 4 measurement, not yet verified
by actually building and solving that geometry.

**NOT MEASURED AT ALL:** whether the nitride strip-test gate passes (check 3
material-model gate, step 8); `mm_width`/`mm_length`; whether the four MMI
output branches are decoupled at any given `waveguide_offset`; runtime cost of
any of the three scripts; mode identity/correspondence between the FDTD and
EME sides of this 2-mode access waveguide (check 2 — flagged in step 4's
writeup as *more* load-bearing here than on any single-mode existing device).

## Order and cost

| step | cluster time | blocking? |
|---|---|---|
| ~~1 — Sellmeier entry + test~~ | none | ✅ **done 2026-08-03** |
| ~~2 — stack generalization + regression test~~ | 1 array job (confirming re-run, job 11659362) | ✅ **done 2026-08-03** |
| ~~3 — fingerprint + test~~ | none (item 2 replaced — see step 3's writeup — costing zero cache/golden invalidation instead of the originally-planned global one) | ✅ **done 2026-08-03** |
| ~~4 — cross-section (`sz`, width, band)~~ | 2 short jobs (11660092, 11660166) + 2 cancelled false starts | ✅ **done 2026-08-03** — `sz=5.0`, `width=1.1` (2-mode, user's revised design point) |
| 5 — body design (`mm_length` sweep) | 1 short job | blocks 6 |
| 6a — geometry plot | none | blocks 6b–d |
| 6b–d — three scripts | 1 array job (3 tasks) | blocks 8 |
| 7 — Makefile/batch/exports | none | — |
| 8 — four checks + material gate | several jobs | — |

The critical path is **1 → 2 → 3 → 4 → 5 → 6**. Steps 1 and 3 need no cluster
time at all, and step 2's library work can be written while step 4's mode solves
queue.

~~**Recommended first action:** step 1 plus step 2's evidence #2 written as a
*failing* test.~~ **Steps 1–4 all landed 2026-08-03.** Every library change
the plan called for shared code is done: nitride is a real material,
`MaterialStack` builds it correctly under both frozen and dispersive paths,
every cache-relevant method fingerprints the real core/clad material, and the
first real nitride mode solves are in: `sz=5.0 µm`, `waveguide_width=1.1 µm`
(a deliberately 2-mode design point — the true single-mode window, 0.2–0.35 µm,
turned out to be too weakly confined to be practical). Step 4 also surfaced and
fixed a real, general bug: `compute_modes_meep` silently returns wrong
(box-mode) answers for weakly-confined structures; switched to
`compute_modes_tidy3d`, cross-validated against a known-good silicon case
first.

**Next action: step 5** — MMI body design. The numbers it needs (`sz`,
`waveguide_width`) are now in hand. Carry forward the design-scope note from
step 4's follow-up: 1.1 µm is 2-mode, not single-mode, so the body/taper
design and step 6's hybrid/EME scripts need `N_MODES` sized and stated for 2
modes at the access waveguides, not assumed single-mode by default the way
every existing MMI example is.
