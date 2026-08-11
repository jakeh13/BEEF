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
| Device dimensions for a 1×4 on nitride | ✅ **DONE 2026-08-04 (steps 4-5)** — `waveguide_width=1.1`, `taper_width=2.2`, `taper_length=5.0`, `waveguide_offset=2.4`, `mm_width=12.0`, `mm_length=206.0` (user-selected — a genuine footprint/throughput/balance tradeoff, not the pure argmin) |
| Example scripts | ⚠ **drafted, structurally validated, not yet solved for real (step 6)** — `mmi_1x4_eme.py`, `mmi_1x4_hybrid.py`. `mmi_1x4_fdtd.py` dropped: `mm_length=206` makes whole-device single-domain FDTD impractical |

So the ordered work is: **material model → stack generalization → cache
fingerprint → cross-section design → body design → two scripts (not three)
→ the four checks.** Steps 1–3 are library changes shared by every device in
the repo and are where the risk is; steps 4–5 are device design; step 6 is
largely copy-work. **Steps 1–5 are now done** (2026-08-03/04) — all of the
design work is closed; what's left is running the two example scripts for
real, wiring (step 7), and the four validation checks (step 8).

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

**⚠ v1 of the design-tool script was WRONG, not just imprecise — MEASURED
2026-08-03/04.** The first attempt built `M_in`/`M_out[k]` by calling
`EMESolver.solve_section()` directly on the taper+access bands (mirroring the
body's own `solve_section_length_sweep()` call, on the theory that an
"all-EME, cheaper than hybrid" design tool could reuse the same convenience
method everywhere). Smoke-test job 11660322 returned **T = 1.0000 exactly**
for every output, every `mm_length`, every wavelength (`total T = 4.0000`,
violating energy conservation on its face).

Root cause, confirmed both empirically and by reading meow's own source
(`create_cells`/`Cell.structures_2d`/`compute_s_matrix`, job 11661894):
`solve_section()` builds exactly **one** `Ls=[length]` cell, samples the
geometry at **one** z location (`Cell.structures_2d`: `z =
0.5*(z_min+z_max)`), and computes **one** mode basis shared by both the
"left" and "right" ports of that cell. `mw.compute_s_matrix()` with only one
cell has no second, distinct cross-section to build a mode-overlap interface
against — the returned block is the **identity** (mode *i* → mode *i*, unit
magnitude, phase only) **by construction**, regardless of how physically
complex the sampled cross-section is. Verified directly: calling
`solve_section()` on this device's own output taper+access band with
`n_modes_in=n_modes_out=6` returned exactly `np.eye(6)`.

This means `solve_section()` is only ever valid over a genuinely **uniform**
region end to end — exactly what `solve_section_length_sweep()`'s own
docstring already requires for the body, and exactly why `mmi_1x2_hybrid.py`/
`mmi_2x2_hybrid.py` use **FDTD**, not EME `solve_section()`, for every taper
band (`INPUT_METHOD = OUTPUT_METHOD = "fdtd"`, only `BODY_METHOD = "eme"`).
There is no all-EME shortcut cheaper than that for a device with real tapers.
The established all-EME device examples don't take one either:
`mmi_1x2_eme.py` calls `EMESolver.solve()`/`solve_length_sweep()` — the
**full multi-cell** device solve (`num_cells=[1,8,1,8,1]`, each taper
genuinely staircased into 8 sub-cells with real mode-mixing interfaces
between them) — never `solve_section()` on a taper span directly.

**Fix (v2), MEASURED working 2026-08-04 (job 11662670):** rewrote
`mmi_1x4_sin_body_length_sweep.py` on `mmi_1x2_eme.py`'s actual pattern — one
`EMESolver` with 5 regions (`wg_in | taper_in | MMI body | taper_out |
wg_out`), tapers subdivided into real cells, one `output_idx` call per output
branch (four calls, no mirror-symmetry shortcut — this device's one input
side is structurally different from its four-branch output side, same
reasoning v1 already had right). The length sweep rides
`solve_length_sweep(sweep_region=2)` — the multi-cell generalization of the
single-cell trick — so it stays free regardless of grid size. Reduced-scope
smoke test (`N_BODY=6`, `N_TAPER=4`, `C_TAPER=2`, `N_WL=3`, 5 lengths) now
returns genuinely length-dependent, non-trivial T: per-output band-mean T =
[0.0924, 0.1210, 0.1210, 0.0924] at `mm_length=40` — note the exact
left-right mirror pairing (output 0 ↔ 3, 1 ↔ 2), which is the geometry's own
y-symmetry showing up correctly, not a bug. `total T = 0.4267` at this
severely truncated mode count is expected to rise substantially at the real
`N_BODY=20` (most of an MMI's self-imaging power lives across many
interfering modes, not the top 6).

**Compute cost, MEASURED:** the smoke config took ≈15–17 min/branch even at
this reduced scope; scaling to the real config (`N_BODY=20`, `N_TAPER=10`,
`C_TAPER=8`, mode-solve cost dominated) projects to roughly a full day if all
four branches ran sequentially in one job. Restructured into one script that
solves a single branch selected by an `MMI1X4_OUTPUT_IDX` env var, submitted
as a 4-task SLURM job array (`run_mmi1x4_body_sweep.batch`, one task per
`output_idx`, `-t 08:00:00` each, run in parallel) — turning a ~1-day
sequential job into a few-hour wall-clock one at the same total CPU cost.
`mmi_1x4_sin_body_length_sweep_merge.py` combines the four saved
`mm_length_sweep_branch{k}.npz` files into the final `mm_length` pick +
evidence plot once all four tasks finish. `N_WL` also reduced from 11 to 5
for this exploratory design-tool pass (final validation is step 6's job, not
this script's).

**Real run completed 2026-08-04 — job array 11667114, all 4 tasks succeeded
(~5:49–6:29 hrs each, within the 8-hr budget).**
`mmi_1x4_sin_body_length_sweep_merge.py` combined the four
`mm_length_sweep_branch{0,1,2,3}.npz` files: pure `argmin(spread)` lands at
`mm_length=280.0 µm` (`spread=0.0013`), a genuine, sharp interior minimum —
verified NOT a sweep-boundary artifact by inspecting the 240–300 µm tail:
spread falls smoothly from 0.13 at both 270 and 290 µm down to 0.0013 exactly
at 280, a clean symmetric dip, not a value still trending toward the 300 µm
edge.

**But the pure spread-minimum is not obviously the best design point.**
Every local minimum in the swept range trades footprint, split balance, and
total throughput against each other differently:

| `mm_length` | spread | relative imbalance | `T_total` | IL |
|---|---|---|---|---|
| 45.0 µm  | 0.0123 | 8.2%  | 0.6024 | −2.20 dB |
| 115.0 µm | 0.0144 | 9.5%  | 0.6095 | −2.15 dB |
| 206.0 µm | 0.0052 | 3.7%  | 0.5696 | −2.44 dB |
| 280.0 µm | 0.0013 | **1.0%** | 0.5211 | −2.83 dB |

(relative imbalance = `spread / (T_total/4)`; the other local minima in the
sweep — 24, 31, 63, 97, 132, 161, 188, 224 µm — all carry ≥40% relative
imbalance and are not competitive.) Longer body → better split precision but
*worse* total throughput and a larger device — the opposite of what a naive
"longer is safer" intuition would suggest, because points further from a
self-imaging condition bleed more power into modes this design tool's
`N_BODY=20` (of 23 guided) truncation doesn't carry cleanly.

**User-selected 2026-08-04: `mm_length = 206.0 µm`** — the balanced middle
option, not the pure `argmin(spread)` winner (280.0 µm). Reasoning: 3.7%
relative imbalance is already good; 280 µm buys the last ~3x reduction in
imbalance at the cost of 0.4 dB more loss and 1.36x the length. Confirmed
exactly on the sweep's own grid (`linspace(20,300,281)`, integer spacing) —
no splice needed, asserted in both final example scripts per the coupler-bug
lesson.

**Design consequence, also decided by the user:** at `mm_length=206 µm` the
device's total length (`≈216 µm` with tapers) makes a **single-domain,
whole-device FDTD simulation impractical** — 3-D time-domain FDTD scales
with volume, and 200+ µm of propagation length in full 3-D is a different
cost class than the ~5 µm devices `mmi_1x2_fdtd.py`/`mmi_2x2_fdtd.py`
simulate whole. **`mmi_1x4_fdtd.py` was dropped from the plan** — this
device ships as EME + hybrid only (`mmi_1x4_eme.py`, `mmi_1x4_hybrid.py`),
not all three methods step 6's original framing assumed. The hybrid script
still gets FDTD validation on the short input/output taper sections; only
the whole-device single-domain FDTD reference is out of scope for this
`mm_length`.

Evidence: `mm_length_sweep.png` (both panels — band-mean transmission per
output vs. length, and log-scale spread vs. length) in
`example_output/verification/mmi_1x4_sin_body_length_sweep/`.

### Step 6 — the example scripts (two, not three — see 6c) ✅ **DONE 2026-08-04**

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

**6c. `mmi_1x4_fdtd.py` — ❌ REMOVED FROM THE PLAN, 2026-08-04 (user decision).**
Originally planned to model on `mmi_1x2_fdtd.py` (`output_monitor_centers`/
`output_monitor_widths` as four-element lists, sized to the actual port pitch
instead of the 1×2's half-domain trick — see the aperture reasoning kept
below, since it still applies to 6d's FDTD sections). Dropped once step 5
settled `mm_length=206 µm`: total device length ≈216 µm makes a
**single-domain, whole-device FDTD simulation impractical** — 3-D
time-domain FDTD cost scales with volume, and 200+ µm of propagation is a
different cost class than the ~5 µm devices `mmi_1x2_fdtd.py`/
`mmi_2x2_fdtd.py` simulate whole. This device ships as **EME + hybrid only**.
The aperture design below still matters for 6d's FDTD input/output
*sections* (short, ~5.5 µm — not the whole device), just not for a
whole-device FDTD reference. **Use `use_neighbor_aperture=True`** on the
`MMI` geometry (not `half_domain_aperture`) — `beef/geometry/mmi.py`'s own
docstring is explicit that the fixed half-domain rule only works when a
simple two-way symmetry split isolates each port, which fails here
(`output_ys` = [-3.6, -1.2, +1.2, +3.6] µm — a y=0 split still groups ports
0&1 together and 2&3 together). `use_neighbor_aperture` instead sizes each
port's aperture to the distance to its *nearest* neighbor specifically,
generalizing past two ports.

**6b/6c/6d crosstalk check — MEASURED 2026-08-04, done in parallel with step
5 (job 11667138), since it does not depend on `mm_length`.** The plan's own
concern here ("the 1×2's 0.65 µm edge-to-edge justification must be
recomputed for nitride, where a weakly-confined mode couples across a gap
that would be safe in silicon") turned out NOT to be a problem at this
design's actual spacing. Even/odd supermode analysis (Ansys/Lumerical
"Evanescent waveguide couplers" method, same as `directional_coupler_mode.py`
— `compute_modes_tidy3d`, not `compute_modes_meep`, per step 4's fix) on two
1.1 µm SIN_300 strips at `waveguide_offset − waveguide_width = 2.4 − 1.1 =
1.3 µm` edge-to-edge gap: `n_eff,even=1.58798`, `n_eff,odd=1.58671`,
`L_c = 614.2 µm`. Over the actual propagation length these four waveguides
sit this close together in the modeled device (`taper_length + stub ≈
5.5 µm`), `T_cross = sin²(πL/2L_c) = 0.000198` — negligible by three orders
of magnitude. Gap sweep confirms comfortable margin even far tighter than
this device's actual spacing (0.5 µm gap: still only 2.3%). **Conclusion:
`waveguide_offset=2.4` does not need revision for crosstalk, and no
mode-basis contamination is expected in either the EME per-branch solves
(already geometrically isolated via `output_idx`/`to_gds_eme`) or the FDTD
monitors (`use_neighbor_aperture=True` is a numerical-isolation convenience,
not compensating for a real physical coupling problem).**

**6d. `mmi_1x4_hybrid.py`.** Model on `mmi_1x2_hybrid.py` — FDTD sections on the
tapered input/output, EME on the uniform body, per-output solves as in 6b.
Decide `N_MODES` and `RESOLUTION` **jointly**, not by copying: README check 3 is
explicit that a shared resolution knob is how a 7% error survived for months, and
the 1×2's check 2 actually *failed* at `N_MODES=12` until it was dropped to 8. A
1×4 body is wider and supports more modes than either existing MMI, so neither
device's numbers transfer.

**Both scripts drafted 2026-08-04, in parallel with step 5's mm_length sweep
(job array 11667114), since nothing about their structure depends on the
exact `mm_length` value; `mmi_1x4_fdtd.py` was drafted the same way and then
removed per 6c above.** `mmi_1x4_eme.py` (`solver.solve()` once per
`output_idx`, `N_BODY=20`/`N_TAPER=10`/`C_TAPER=8` carried over from the
now-validated step 5 design tool) and `mmi_1x4_hybrid.py` (FDTD input/output
+ EME body, `use_neighbor_aperture=True`, four `hybrid.solve()` calls
sharing one `cache` dict, `N_MODES=20`) initially carried `MM_LENGTH = None`
with a guard that raised `NotImplementedError` until step 5's real value was
filled in — deliberately not a silent placeholder. **Structurally validated
at that stage** (not solved): a temporary copy of each script with a
stand-in `mm_length=150.0` and an early exit right after its geometry plot
ran clean through MMI construction, solver instantiation
(`FDTDSolver`/`EMESolver`/`HybridSolver`), and `plot_*_geometry` (including
`mmi_1x4_hybrid.py`'s `sub_region_configs`, one per output branch) with no
errors. **`MM_LENGTH=206.0` filled in 2026-08-04** once step 5 landed, spliced
onto the sweep's own grid and asserted exact (coupler-bug lesson) in both
scripts; re-ran the same structural check with the real value — still clean.
`mode_solver="tidy3d"` throughout (never `"meep"`) per step 4's fix,
including in the hybrid's EME body sub-solver, where the established
precedent (`mmi_1x2_hybrid.py`/`mmi_2x2_hybrid.py`) uses `"meep"` — flagged
in `mmi_1x4_hybrid.py`'s own docstring as a `calibrate_body_neff` risk to
revisit at step 8 if the engine-agreement check finds it matters.

**`mmi_1x4_eme.py` run for real 2026-08-04 — MEASURED, job 11678090.**
Sequential per-branch execution would have cost ~1 day+ (each `solve()`
call redoes all 19 mode-solve cells across 11 wavelengths — see the
session's own cost breakdown), so the run was split the same way step 5
was: `_prewarm_mmi_1x4_eme.py` (temporary, `MMI1X4_OUTPUT_IDX` env var)
solved each of the 4 output branches into the persistent `save_dir` cache
via a 4-task SLURM array (job 11672738) — **all 4 tasks succeeded in
~2:45 hrs each**, well under the pessimistic ~13 hr/branch extrapolation
from step 5's own design-tool timing. The real, unmodified `mmi_1x4_eme.py`
then ran against the warm cache in 38 seconds (job 11678090), producing the
actual deliverable (`smatrix.pkl`, `performance.png`):

| output | y (µm) | T | IL (dB) |
|---|---|---|---|
| 0 | −3.60 | 0.1415 | −8.49 |
| 1 | −1.20 | 0.1598 | −7.97 |
| 2 | +1.20 | 0.1598 | −7.97 |
| 3 | +3.60 | 0.1415 | −8.49 |
| **total** | | **0.6025** | **−2.20** |

(at λ=1.550 µm, band center of the 11-point C-band sweep). Exact left-right
mirror symmetry (0↔3, 1↔2) holds, matching both the prewarm results and the
design tool's own step-5 measurements at this `mm_length` (band-mean
`T_total=0.5696`, `IL=−2.44` dB there vs. this single-center-wavelength
`T_total=0.6025`/`IL=−2.20` dB here — a mean-over-band vs. single-wavelength
comparison, not a discrepancy). `mmi_1x4_eme.py` is now genuinely done, not
just structurally validated. **Uses `res=40`, the still-open resolution
question from step 8 check 3's preliminary finding above — a real caveat on
these numbers until that check lands.**

**`mmi_1x4_hybrid.py` — first prewarm attempt (job array 11672739,
`N_MODES=20`, `-t 24:00:00`) TIMED OUT 2026-08-05, MEASURED.** All 4 tasks
hit the wall-time limit having only reached MPB band 14 of 20 on the FIRST
of the two body-adjacent interfaces (input↔body) — the second interface
(body↔output) and all FDTD time-stepping hadn't even started. Extrapolating
the observed rate (~14 bands in 24h) puts the true per-branch cost at
60-80+ hours, not the ~24h the walltime was sized against; the original
4-parallel-job structure also redundantly repeats the shared input+body
work (identical across all 4 branches) once per task instead of once
total, since nothing enforces an ordering between the 4 racing jobs.

**Fix, user-directed: reduce `N_MODES` (`interface_modes`), not just add
walltime.** In `HybridSolver`, `interface_modes` directly sets
`solve_section_length_sweep()`'s `n_modes_in`/`n_modes_out` for the body
section — i.e. it IS the body's resolution in this cascade, not an
independent interface-only knob (the `EMESolver(num_modes=...)` constructor
value for `body_solver` is an unused placeholder for this call path, exactly
as `mmi_1x2_hybrid.py`'s/`mmi_2x2_hybrid.py`'s own comment already said).
So there is no way to keep the body at the standalone EME script's 20 modes
while only shrinking the FDTD interface count — reducing `N_MODES` to 12
genuinely trades some self-imaging accuracy for tractability, the same
tradeoff `mmi_1x2_hybrid.py` (`N_MODES=8`) and `mmi_2x2_hybrid.py`
(`N_MODES=14`) already accept relative to their own fully-guided mode
counts. **Preliminary, unvalidated choice — revisit at step 8's
resolution/mode-count convergence check**, alongside the still-open
`res=40` question. MPB's own log from the timed-out run showed near-cutoff
modes taking far more iterations to converge than well-confined ones (up to
191 vs. 8-10), so dropping the marginal top modes should cut cost by more
than the naive 12/20 ratio — not yet confirmed.

**Second prewarm attempt submitted 2026-08-05: job array 11696233**
(`N_MODES=12`, `-t 48:00:00` — still an estimate, no data yet on the second
interface or actual FDTD time-stepping cost, chosen with real margin over a
naive linear scaling of the first attempt's partial-progress rate). Will be
run for real the same way `mmi_1x4_eme.py` was once that array lands.

**⚠ MEASURED, ~15.5 hrs in: `N_MODES=12` fixed the MPB-band-count bottleneck,
but exposed a SECOND, larger one — `decay_tol=1e-11`'s time-stepping cost on
a weakly-confined material.** All 4 tasks: input section's FDTD run
(`Computational cell is 6.75 x 21 x 6`, the short taper+stub, NOT the whole
216 µm device) took **`run 0 finished at t = 213.275 (17062 timesteps)`
— roughly 11.85 hrs of real wall time** (17062 steps at the observed
~2.3-3.0 s/step, accelerating slightly near the end) just to satisfy
`decay_tol=1e-11` on this one short section. That is a real, direct
confirmation of A7's original prediction ("more PML clearance and a longer
decay tail... `decay_tol` may need revisiting") — not the MPB band count
diagnosed from the first (`N_MODES=20`) timeout, which is a separate cost
that `N_MODES` reduction did successfully fix (this run reached time-stepping
within ~20 min, vs. never in 24h before). Both bottlenecks were real; fixing
one exposed the other.

After `run 0` finished, all 4 tasks moved into the input section's mode
decomposition (`get_eigenmode_coefficients`-driven `MPB solved for
frequency_N(...)` calls, N_MODES=12 bands × 11 wavelengths = 132 solves per
task) — reached band ~9-10 of 12 by ~15.5 hrs, i.e. ~3.5 hrs for ~9 bands,
roughly 25-30 min/band. **Extrapolated remaining budget, not yet confirmed**:
~1-1.5 hrs to finish input decomposition, then the body's EME solve (tidy3d,
no time-stepping — expected fast), then the OUTPUT section's own FDTD run
(same cell class, so plausibly another ~12 hrs if similarly weakly-confined)
plus its own ~3.5 hr decomposition. Rough total estimate **~30-33 hrs**,
inside the 48 hr budget with some margin — but this assumes the output
section's decay behaves like the input section's, which has NOT been
measured and is the next thing to check once it starts. All 4 array tasks
redundantly repeat the input+body work independently (flagged as a known
inefficiency, not fixed here — see the first-attempt diagnosis above), so
this ~12 hr cost is paid 4 times in parallel, not once.

**Update, ~26 hrs in: the output section's own FDTD run has started —
confirms the "two comparable time-stepping phases" extrapolation above.** A
second `Computational cell is 6.75 x 21 x 6 with resolution 40` line appears
in all 4 logs (input decomposition + the body's tidy3d EME solve both
completed in between, no separate log markers for either but both are fast
by construction). Tasks 1-3 are at t≈108-111 of what was 213.275 for input's
own full decay (steps ~8700-8940 of the ~17062 input needed) — roughly
halfway, and at the same ~2.3 s/step rate, consistent with the output
section decaying similarly to the input section as hypothesized. **Task 0 is
notably behind** (t≈41.7, step 3339) — same task that was slowest through
input decomposition too; worth watching in case it is heading for the 48 hr
wall, though nothing indicates a stall (steady step rate, not stuck).
Revised total-time picture: if output decays like input, tasks 1-3 finish
output time-stepping around ~32-35 hrs elapsed, then need output
decomposition (~4-8 hrs by input's own pattern) — landing somewhere around
~38-43 hrs, inside the 48 hr budget with real but not huge margin. Task 0's
slower pace is the actual risk to watch, not the other three.

**Update, ~38 hrs in: output section's FDTD run finished, all 4 tasks —
`run 0 finished at t = 223.8625 (17909 timesteps)`, essentially matching
input's own `t = 213.275 (17062 timesteps)` (~5% longer).** Confirms the
"output decays like input" hypothesis directly rather than by extrapolation.
Task 0, previously the laggard, caught back up enough to finish its own
output run at the same elapsed time as the other three. All 4 tasks are now
in output-section mode decomposition (bands 8-11 of 12 as of this check) —
the last phase before each task's own cache write completes. Comfortably
inside the 48 hr budget at this point.

**All 4 tasks COMPLETED 2026-08-07** (task 1: 1d16h57m; task 2: 1d16h59m;
task 3: 1d16h49m; task 0, the persistent laggard: 1d21h13m) — well inside the
48 hr budget. `mmi_1x4_hybrid.py` run for real (job **11734486**) against the
now fully-warm cache: **23.5 s**, hit cache for every section of every
branch.

| output | y (µm) | T | IL (dB) |
|---|---|---|---|
| 0 | −3.60 | 0.1541 | −8.12 |
| 1 | −1.20 | 0.1409 | −8.51 |
| 2 | +1.20 | 0.1409 | −8.51 |
| 3 | +3.60 | 0.1541 | −8.12 |
| **total** | | **0.5898** | **−2.29** |

(at λ=1.550 µm, band center.) Mirror symmetry (0↔3, 1↔2) holds, and total
throughput is close to EME's own 0.6025 (−2.20 dB) — within ~2%, comparable
to the ~0.5% engine bias measured on the existing MMIs, though on the higher
side.

**⚠ But the SPLIT SHAPE disagrees, not just the magnitude — flagged here,
not yet explained.** EME sends MORE power to the inner ports (1, 2: 0.1598)
than the outer ports (0, 3: 0.1415). Hybrid does the OPPOSITE: outer ports
(0, 3: 0.1541) get MORE than inner (1, 2: 0.1409). Both preserve the
device's own mirror symmetry and both sum to a similar total, but which pair
of ports is favored flips between engines — a qualitative disagreement, not
just a quantitative one. Leading unconfirmed hypothesis: `mmi_1x4_hybrid.py`
uses `N_MODES=12` for the body's self-imaging interference (reduced from the
EME script's own `N_BODY=20` purely for the prewarm's tractability, per the
step-6 narrative above) — a coarser modal basis could plausibly redistribute
power between the near-imaging output positions differently, not just lose
some to radiation. **Not yet verified.** This is exactly the kind of
disagreement step 8 check 1 (engine agreement) exists to characterize;
record it there rather than resolving it here from reasoning alone.

Temporary prewarm files (`_prewarm_mmi_1x4_hybrid.py`,
`batch/verification/run_prewarm_mmi1x4_hybrid.batch`) deleted 2026-08-07,
same cleanup EME's own prewarm files got once its real run succeeded.

#### Field visualization + a real, quantified port-alignment finding — MEASURED 2026-08-08

`beef.plot_eme_fields()` existed in the codebase (`beef/plot/fields.py`) but
had **zero call sites** anywhere in examples/tests before this session.
`examples/verification/python/plot_mmi1x4_eme_fields.py` (new) rebuilds the
same `cells`/`modes` `EMESolver.solve()` builds internally
(`beef/solvers/eme.py:765-771`, discarded there once the S-matrix is
extracted) for one wavelength/branch and calls the plotting utility — job
**11736246**, ~18 min, 19 tidy3d eigensolves, no time-stepping. Produces a
textbook self-imaging field pattern for output branch 0 at λ=1.55 µm: input
focus at z=0, interference fan through the body, and a **single, tight
refocus at z≈165 µm, y≈0**.

**That interior refocus is NOT a 4-way split image and is not evidence a
shorter body would help — confirmed, not assumed, with a follow-up
annotated analysis (`analyze_mmi1x4_self_imaging.py`, job 11750449, ~16
min).** The four real output ports sit at y = ±1.2, ±3.6 µm — none at y=0 —
so a bright lobe centered on y=0 partway through the body is the classical
"restored single image" every MMI shows at a fraction of its beat length,
independent of the intended N-way split. Confirmed directly: annotating the
same field with the real region boundaries (computed from the solver's own
region lengths, not eyeballed: body spans z=[5.5, 211.5] µm, taper_out
[211.5, 216.5]) shows the interference pattern re-fanning back out to the
body's full ±6 µm width after the z≈165 spot, well before reaching the
actual taper mouth at z=211.5. This is also consistent with — not merely
argued from — the already-measured 281-point `mm_length` sweep (step 5):
lengths near this interior spot (161, 188 µm) were explicitly measured as
**not competitive** (≥40% relative imbalance), while 206 µm (landing the
cut at the body's real end, not at the interior single-image point) is one
of the genuinely competitive candidates. The field plot, read correctly,
agrees with the sweep rather than contradicting it.

**A real, previously unquantified finding, from the same job: a transverse
cross-cut of the field taken exactly at the taper mouth (z=211.5 µm, not an
arbitrary interior z) shows the self-imaging pattern's four peaks sitting
systematically ~0.3–0.4 µm OUTWARD of the four ports' actual centers**
(±1.2, ±3.6 µm), symmetric on both the inner and outer pairs. Against a
1.1 µm-wide access waveguide and a 2.4 µm port pitch, an offset of this size
is not negligible — plausibly a real contributor to the measured 3.7%
relative imbalance / 0.5696 total T at this `mm_length`, though **not
confirmed as such**: this is one wavelength, one output branch, no sweep.
**Not yet actioned or further investigated — flagged here as a genuine,
evidence-based candidate for a future `waveguide_offset` sweep** (analogous
to step 5's own `mm_length` sweep), which has never been run for this
device. If real, closing a ~0.3–0.4 µm mismatch could plausibly recover some
of the gap between the chosen 206 µm point (3.7% imbalance) and 280 µm's
1.0% — cheaper than lengthening the device, if the two effects are actually
independent (also not yet checked).

#### Hybrid variant of the step-5 `mm_length` sweep — MEASURED 2026-08-09

Step 5's own sweep (`mmi_1x4_sin_body_length_sweep.py`, job array 11667114)
used the pure-EME chain, at real cost (~5:49-6:29 hrs/branch, no FDTD section
to lean a cache on) — that sweep's config is already bit-for-bit
`mmi_1x4_eme.py`'s own, so it already **is** the EME variant; nothing new
needed there.

The hybrid chain is different: `HybridSolver.solve_length_sweep()`'s own
docstring promises the two FDTD sections dedupe against `solve()`'s own
cache entries when configs match exactly, so a hybrid `mm_length` sweep
pinned bit-for-bit to `mmi_1x4_hybrid.py`'s real config (same `mm_length=
206.0`, same `wavelengths`, same `N_MODES=12`) should reuse the cache jobs
11696233/11734486 already warmed, paying only for the body's own cheap
tidy3d mode solve. New script `mmi_1x4_hybrid_body_length_sweep.py` written
on that basis.

**A real bug, caught and fixed before it did lasting damage.** First
submission (job 11752875) computed `save_dir` via `os.path.join(_HERE, "..",
"..", "..", "example_sim_cache")` — **three** `".."`, copied from the
pre-existing EME sweep tool's own path arithmetic without checking it. Three
dots from `examples/verification/python/` lands at `duty/BEEF/
example_sim_cache/` — a real, separate, mostly-empty directory (the EME
tool's own past runs had silently been writing there too, harmlessly for
that script since it never expected cache reuse against the real device
cache anyway). For THIS script the same bug was not harmless: pointed away
from the actually-warm `duty/example_sim_cache/`, every section — including
the two ~12-hr FDTD sections — cache-missed, and the job timed out at 1 hr
mid-time-step, having made no useful progress. Root-caused via a sequence of
`cache_only=True` probes (safe by construction — a miss raises immediately
rather than silently falling through to a real solve) that isolated it to
the path, not to anything about looping over branches or the 281-point grid
(both tested independently and shown fine in isolation once the path was
right). Fixed to four `".."` (matching `mmi_1x4_hybrid.py`'s own depth from
`examples/mmi_1x4/python/` — the same depth as `examples/verification/
python/`, four levels either way). **Not fixed**: the EME sweep tool's own
three-dot path — genuinely inert there, out of scope of what this fix was
for.

**Re-run clean (job 11770532): 25.3 s for all 4 branches, 281 lengths.**

**⚠ SECOND real bug, found 2026-08-09/10 — this run's own numbers were
WRONG, not just the path.** `HybridSolver._stubs(k)` returns `(0.0, 0.0)`
for any section that is neither the very first nor the very last in the
chain — the body (section 1 of [input, body, output]) is neither. So
`solve_section()`'s own `length = (x1-x0) + stub_left + stub_right`
(what a plain `hybrid.solve()` call uses for the body) reduces to just
`x1-x0` = `mm_length - 2*body_stub`, **not** `mm_length` — `body_stub` of
MMI body is folded into each *adjacent* FDTD section's own domain instead
(module docstring's "Build HybridSolver" comment), so
`boundaries[1]`/`boundaries[2]` are already built narrower to account for
that. `solve_section_length_sweep()` (what the swept path uses for the
body) applies whatever `lengths` it is given directly as the propagation
length, with **no** such stub adjustment — so the caller must pre-subtract
`2*body_stub` to mean the same physical length `solve_section()` would at a
given `mm_length`. `mmi_2x2_hybrid.py`'s own hand-rolled cascade already
does exactly this (`body_lengths = mm_lengths - 2 * body_stub`); this
session's new scripts, built on the generic `HybridSolver.solve_length_sweep()`
API instead of hand-rolling, missed it the first time through.

**MEASURED impact**: passing raw `mm_lengths` (not `mm_lengths -
2*body_stub`) gave a real, reproducible ~6% relative error on individual
port T at the chosen point (0.1635 vs the correct 0.1541 at output 0,
λ=1.55, mm_length=206) — large because T(L) oscillates rapidly with length
here (visible in the sweep plot's own many closely-spaced peaks/dips), so
even a 0.5 µm propagation-length error lands on a materially different
point of the curve. Root-caused by direct comparison: `hybrid.solve()`,
`hybrid.solve_length_sweep()` with a single-element `body_lengths=[205.5]`,
and the full 281-point sweep's own slice at `mm_length=206` all agree
exactly (max diff 0.0 between the two swept variants, confirming the grid
itself introduced no error) once the fix landed; before the fix, `solve()`
and the swept path disagreed by 0.0094 in T at the same nominal length.
Fixed in both `mmi_1x4_hybrid_body_length_sweep.py` and
`mmi_1x4/python/mmi_1x4_hybrid.py` (the production script, see step 6's own
account below); **re-verified after the fix: the production script's
single-length results are now byte-identical to the pre-sweep-integration
baseline** (job 11783136 reproduces job 11734486's T=[0.1541, 0.1409,
0.1409, 0.1541] exactly).

**Corrected re-run (job 11782082, 6:56 — slower than the buggy version's
25.3 s because `body_lengths` is a genuinely new cache key, needing one
fresh body mode solve; still cheap, not the ~12 hr/section FDTD class):**

| | EME (job 11667114/step 5) | Hybrid, corrected (job 11782082) |
|---|---|---|
| sweep optimum (min spread) | 280.0 µm (spread 0.0013) | **44.0 µm** (spread 0.0041) |
| spread at chosen 206 µm | 0.0052 | 0.0064 |
| total T (band-mean) at 206 µm | 0.5696 | 0.6014 |
| per-output band-mean T at 206 µm | [0.1398, 0.1450, 0.1450, 0.1398] | [0.1535, 0.1471, 0.1471, 0.1535] |
| outer vs. inner at 206 µm | inner **>** outer | outer **>** inner |

**The single-argmin optima now look 236 µm apart — but the overlay plot
shows this is NOT the disagreement it looks like.** Re-ran
`mmi_1x4_length_sweep_compare.py` (no new solving, just re-reads both
`.npz` files) against the corrected hybrid data:
`example_output/verification/mmi_1x4_length_sweep_compare/
eme_vs_hybrid_length_sweep.png`. **Both total-throughput and spread curves
track almost exactly across the ENTIRE 20–300 µm range** — same
oscillation period, same dip depths, at essentially every point, including
near 44 µm: EME's OWN spread curve also dips sharply there (visibly
comparable in depth to its own 280 µm dip), it just isn't EME's single
*global* minimum. The "236 µm apart" framing describes two engines' argmin
picks landing on two different (both genuinely excellent) local minima of a
curve they otherwise agree on almost everywhere — not a physics
disagreement. **This is a real, reassuring confirmation that both engines
see 44 µm as a genuinely strong candidate**, not an artifact of one engine
or one bug.

**The split-SHAPE discrepancy, however, survives the fix and is real**:
at the corrected numbers, EME still sends more band-mean power to the inner
ports at 206 µm (0.1450 > 0.1398), hybrid still sends more to the outer
ports (0.1535 > 0.1471) — opposite signs, not a magnitude difference, and
not explained by either bug found so far (both were length/path bugs, not
anything that would flip a sign). Well-evidenced, reproducible, and still
unexplained — exactly what step 8 check 1 (engine agreement) exists to
investigate, not yet run for this device.

**New direction opened by the corrected data**: mm_length=44 µm is a
genuine, dual-engine-confirmed candidate for a **>4x shorter** multimode
section than the chosen 206 µm, at comparable (even slightly better, on the
hybrid side) split quality. But the current `waveguide_offset=2.4 µm` port
spacing was tuned against 206 µm's own self-imaging pattern, which is a
physically different point on the beat curve than 44 µm's — nothing has
yet checked whether a *different* offset would align better with 44 µm's
own image spacing. `mmi_1x4_sin_taper_offset_sweep.py` (new) fixes
`mm_length=44.0` and sweeps `waveguide_offset` across the geometrically
feasible range `(2.2, 3.2667)` µm (uniform port pitch, per-offset genuinely
fresh EME solve — no analytic trick available here, unlike length) to
answer that directly.

**Targeted, not blind: `analyze_mmi1x4_self_imaging_L44.py` (job 11783781,
17 min, same field-snapshot method as the 206 µm version) read the field's
actual lobe positions at mm_length=44's own body end first.** Found **7
peaks**, not a clean 4: y = ∓5.01 (h132), ∓3.69 (h83), ∓1.74 (h159,
tallest), 0.01 (h58, weak, no port there anyway). Best 4-port match: outer
pair → ±5.01 peak, inner pair → ±1.74 (tallest) peak → implied
`waveguide_offset≈3.34` µm — just past the geometric feasibility ceiling
(3.2667 µm, set by `mm_width=12`/`taper_width=2.2` fixed) by ~0.07 µm. Real
per-point cost checked first (job 11783183: 2:48:51 for ONE offset, ONE
branch, full config) — far too expensive for a grid scan, so the actual
sweep used the field-plot estimate to center a **5-point** sweep near the
feasible ceiling (`[3.0, 3.217]` µm) at **coarse `C_TAPER=3`** (down from
8) and **center wavelength only** (`N_WL=1`, not the full 11-point band) —
~38 min/branch, 4 branches parallel (job 11789245).

**Result — a real, substantial throughput gain, at a real balance cost:**

| offset (µm) | outer T | inner T | total T | spread |
|---|---|---|---|---|
| 3.000 | 0.195 | 0.181 | 0.751 | 0.028 |
| 3.054 | 0.197 | 0.182 | 0.758 | 0.031 |
| 3.108 | 0.198 | 0.183 | **0.762** | 0.032 |
| 3.163 | 0.198 | 0.183 | 0.762 | 0.030 |
| 3.217 | 0.196 | 0.183 | 0.759 | **0.027** |

Total T is **0.75–0.76 across the whole tested range** — vs. the current
production point's (mm_length=206, offset=2.4) 0.6014 band-mean —
**~27% more light delivered**, a large, consistent effect across all 5
points, not a fluke of one. But balance is worse: spread 0.027–0.032 here
(relative imbalance ≈15%, since per-port mean ≈0.19) vs. production's
0.0064 (≈3.7% relative) — a genuine throughput/balance tradeoff, not a
strictly-better point. **Caveats, stated plainly**: single wavelength
(no band-average/dispersion check), coarse `C_TAPER=3` (taper-staircase
discretization error unmeasured), and the tested range didn't reach the
peak-matched ~3.34 µm estimate (infeasible at the current `mm_width`).
**Not yet a validated design change** — flagged as a promising, real,
measured direction for this device (short body + widened offset, a
markedly more compact and higher-throughput geometry than the current
206 µm/2.4 µm point), worth a real-config confirmation (`N_WL=11`,
`C_TAPER=8`) and possibly a `mm_width` increase (to reach the true
peak-matched 3.34 µm offset) before treating it as a candidate replacement
design.

#### `mmi_1x4_eme.py`'s own length-sweep plot — code landed, real run NOT YET COMPLETE ⚠

Mirroring what landed in `mmi_1x4_hybrid.py` (below), `mmi_1x4_eme.py`
Section 4 was rewritten to call `EMESolver.solve_length_sweep(...,
sweep_region=BODY_REGION=2)` instead of plain `solve()`, deriving the
existing single-`MM_LENGTH=206.0` report AND a new `mm_length_sweep.png`
(same per-output/spread format as the hybrid version) from one computation
— code is written, syntax-checked, **not proven against a real run yet.**

**Real cost is materially higher than the single-length prewarm, and the
first attempt TIMED OUT.** `solve_length_sweep()` caches under its own
method name — a genuinely different cache entry from `solve()`'s, per
`EMESolver.solve_length_sweep()`'s own docstring — so this needs a fresh
prewarm regardless of the already-warm `solve()`-based cache from step 6's
own real run. `_prewarm_mmi_1x4_eme_sweep.py` (new, temporary, mirrors the
original single-length prewarm's per-branch/array pattern) submitted as job
**11783204** (4-task array, `-t 480` / 8h). **All 4 tasks TIMED OUT at 8h,
zero progress markers visible in any log** (only the one "Prewarming
branch…" print line — `mode_solver="tidy3d"` gives no per-step output the
way meep's MPB does for the hybrid script's FDTD sections, so there is
no way to tell from the log alone how close it was). Not a stall/bug by
itself — the pre-existing EME design tool (`mmi_1x4_sin_body_length_sweep.py`,
step 5) measured **5:49–6:29 hrs at `N_WL=5`** for the equivalent full
19-cell-chain `solve_length_sweep()` call; linearly extrapolating to this
script's real `N_WL=11` (est. cost scales with wavelength count, since the
"free" trick only removes the length-count dependency, not the
per-wavelength mode-solve cost) suggests **~13–14 hrs/branch** — i.e. 8h
was always going to be short, this is consistent with the design tool's own
prior measurement, not a new anomaly.

**Left in this state for the next session to pick up.** Options, not yet
decided between:
1. Resubmit the same array with a longer walltime (`-t 900` / 15h, or split
   further) — most direct, costs real cluster time, no code change.
2. Reduce `N_WL` for this specific sweep plot (it is a design/diagnostic
   artifact, not a load-bearing result the way the single-`MM_LENGTH=206`
   report is) — e.g. `N_WL=5`, matching the design tool's own already-paid
   cost, cutting the fresh solve to the already-measured ~6 hr/branch class.
3. Skip re-deriving from scratch and instead have `mmi_1x4_eme.py` read the
   design tool's own already-computed `mm_length_sweep.npz` (step 5, job
   11667114 — genuinely the same solver config, confirmed identical
   N_BODY/N_TAPER/C_TAPER/oxide dims/res/mode_solver) for the PLOT only,
   keeping the existing fast `solve()`-based single-length path for the
   printed/saved production results — zero new cost, at the price of the
   plot being sourced from a companion script's run rather than this one's
   own, same asymmetry the hybrid/EME scripts already have in a few other
   places (documented, not hidden).

No cluster time is currently committed to this (the timed-out array is
finished, not paused) — next step is picking one of the three options
above and running it, not resuming anything in flight.

Companion field-plot scripts written the same session for the other
existing EME examples, confirming `plot_eme_fields()` generalizes cleanly:
`plot_mmi1x2_eme_fields.py` (job 11745434, 19 cells, ~1 min against a
lightly-loaded node — cost varies more with node contention than the
~18 min mmi_1x4 measurement suggested), `plot_mmi2x2_eme_fields.py` (job
11745435, 33 cells, ~2 min) — both produced clean, textbook self-imaging
patterns, visually confirmed. `plot_mmi4x4_90hybrid_eme_fields.py` (a
different, concurrently-worked device — see its own
`mmi_4x4_90_hybrid_validation.md`, not this doc) also written and run at
that device's real config (job 11750455, 63 cells, ~4 min) — notably faster
than that device's own docstring warned (~416 s/mode measured for its body
cross-section alone, Step 3a); not investigated further here since it
belongs to that other validation thread.

### Step 7 — wiring ✅ **DONE 2026-08-05**

- `BEEF/examples/Makefile`: add `MMI_1X4_EXAMPLES = mmi_1x4_eme mmi_1x4_hybrid`
  (two, not three — `mmi_1x4_fdtd` was dropped, step 6c), add to `EXAMPLES` and
  `DEVICES`, add the `mmi_1x4` aggregate target and the two per-script targets.
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

#### What landed

- `BEEF/examples/Makefile`: `MMI_1X4_EXAMPLES = mmi_1x4_eme mmi_1x4_hybrid`
  (two, not three — `mmi_1x4_fdtd` dropped per step 6c), added to `EXAMPLES`
  and `DEVICES`, plus the `mmi_1x4` aggregate target and the two per-script
  targets (`mmi_1x4_eme`, `mmi_1x4_hybrid`). Verified with
  `make -n mmi_1x4_eme mmi_1x4_hybrid mmi_1x4` — prints the expected two
  `srun`/`mpirun` command lines, nothing else.
- `batch/run_mmi_1x4_examples.batch` — new, copied from
  `run_mmi_2x2_examples.batch`'s array-of-`make`-targets pattern
  (`--account=gts-sr108-ece`, `EXAMPLES=(mmi_1x4_eme mmi_1x4_hybrid)` in
  Makefile index order, per-task `rm -rf .../example_output/mmi_1x4/${NAME}`
  — safe here since this device has no `*_FROZEN` baseline to protect,
  unlike the coupler). **Resourced differently from the 1x2/2x2 examples
  batches on purpose**: those request `-N 2 --ntasks-per-node=24` because
  their array includes a whole-device FDTD run; this device ships EME +
  hybrid only, both of which are cache-hit-dependent solves once prewarmed,
  so it mirrors `run_mmi1x4_eme_real.batch`'s much smaller footprint
  (`-N 1 --cpus-per-task=8`) instead. `-t 60` is sized off the EME side's
  real measured warm-cache cost (38 s, job 11678090); the hybrid side's
  warm-cache cost is still unmeasured (prewarm job array 11696233 was still
  running at the time this batch script was written) — flagged in the
  script's own comment to revisit if that task times out.
- `beef/__init__.py`/`beef/materials/__init__.py` exports for `SIN_300` —
  already done as part of step 2/4's work (`si3n4` platform singleton and
  its docstring/`__all__` entries), confirmed present, no further change
  needed.

### Step 8 — the four checks

Only now does this device get a normal gap list. From
[README.md](README.md#the-four-checks), in the order that makes each one
cheap:

| check | 1×4-specific note |
|---|---|
| **2. Mode identity** | Do this **first** — it is mode solves only, no time-stepping, and it is the check that actually failed on the 1×2. A wider nitride body supports more near-degenerate modes than either existing MMI, so this is *more* likely to fail here, not less. **✅ MEASURED 2026-08-05, and the prediction landed.** `audit_interface_modes.py` extended for this device (job **11697897**, res=40, N_MODES=12, wl=1.55 — pinned exactly to `mmi_1x4_hybrid.py`'s real config, on the body cross-section at `mm_width=12.0`). One real bug caught while wiring it up: `meow_neff()` hardcoded `mode_solver="meep"`, but this device's own hybrid script deliberately uses `"tidy3d"` for its EME body sub-solver (step 4's MPB-wrong-for-nitride finding) — auditing with meep would have checked a basis the real `HybridSolver` chain never builds. Fixed by threading `cfg.get("mode_solver", "meep")` through (default unchanged, inert for the other three devices). **Result: 12/12 guided; consistent meep-vs-meow n_eff bias of −0.0022 to −0.0024 (bands 0–8), narrowing to −0.0008 (bands 9–11) — comparable in magnitude to the ~0.5% engine bias measured on both existing MMIs. But `ORDER PRESERVED: False`**: the naive nearest-neighbour-in-n_eff correspondence maps BOTH meep bands 0 and 1 onto the same meow band (1) — the top two guided modes are only 0.0036 apart in n_eff, smaller than the systematic engine bias itself, so a 1-D scalar match cannot resolve them. This is exactly the crowding failure `mode_identity.py`'s own docstring names (observed before on the modal multiplexer's bands 8–10) and exactly what this row predicted — the failure just shows up at the very TOP of the spectrum here rather than the middle. **Resolved, ✅ MEASURED 2026-08-06, job 11706990 — a scalar-check artifact, not a real mismatch.** `verify_mmi1x4_mode_overlap.py`'s field-overlap check (same `mode_identity.py` machinery, `cmp_res=20` per its own cost-driven default — deliberately a different resolution than job 11697897's res=40, since a shape-based check does not need matched-resolution agreement the way nearest-neighbour n_eff does) found **`ORDER PRESERVED: True`, worst overlap on the matched diagonal 0.9979** — every one of the 12 guided bands corresponds 1:1 (identity permutation) between meep and meow, with overlap magnitudes 0.990–0.998 across the board, no ambiguous pairing anywhere. The top two modes that fooled the scalar check (0.0036 apart in n_eff, smaller than the ~0.0023 engine bias) are, by mode SHAPE, unambiguously distinct and correctly ordered — exactly the failure `mode_identity.py`'s docstring predicts nearest-neighbour-in-n_eff cannot resolve but a field-overlap match can. **Check 2 passes for this device.** The interface basis `HybridSolver` actually contracts is sound; the earlier `ORDER PRESERVED: False` finding is now understood and closed, not an open risk. |
| **3. Resolution convergence, per engine** | Sweep FDTD and EME resolution **independently**. Sweep a Δβ-like quantity, not just output power. Expect nitride to need different resolution than silicon — lower contrast means the index step is easier to resolve, but the mode is larger, so the two effects push opposite ways and it must be measured. **⚠ Preliminary attempt made 2026-08-04 (job 11675082), inconclusive, deferred — not the real check.** Swept `res=[15,20,25,30,35,40]` on the body cross-section (mm_width=12.0, 23 guided modes), tracking modes 0/10/19/22 (fundamental through near-cutoff). Result: `n_eff` does NOT converge monotonically anywhere in that range — successive `\|dn_eff\|` bounce around ~1-2e-3 instead of shrinking (e.g. mode 0: 1.6720→1.6715→1.6725→1.6737→1.6726→1.6724), unlike step 4's clean monotonic `sz` convergence. This means `res=40` (what step 6's real example runs currently use) has **not actually been shown to be converged either** — the scatter spans the whole tested range, not just the low end. Leading hypothesis, NOT yet verified: mode-index correspondence swapping between near-degenerate higher-order modes in this densely-multimoded cross-section (23 guided modes in a fairly narrow index range) rather than a genuine solver-convergence failure — this codebase's own `mode_identity.py` field-overlap tools could confirm/refute this directly but have not been run against it. **User decision 2026-08-04: proceed with `res=40` as a working assumption for step 6's real runs (already in flight, jobs 11672738/11672739), run the real resolution-convergence check afterward.** This is a real open risk on the step-6 results until that check lands, not a closed question. |
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
| ~~7 — Makefile/batch/exports~~ | none | ✅ **done 2026-08-05** |
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
