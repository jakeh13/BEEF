# 1×4 MMI on 300 nm Si₃N₄ — build plan

> ## ✅ RESOLVED 2026-08-31 — all four checks measured; three residuals carried
>
> Moved to `resolved/` under this directory's own bar: *"whole docs
> sufficiently addressed (open to-dos may remain)"*. **Read the three
> residuals below before quoting any number from this device.** They are
> stated here rather than left at line ~1869 so they travel with the doc.
>
> **Update 2026-09-05 → 2026-09-13, see "✅ 9d RESULT" below.** The
> `interface_modes` sweep this header's own to-do list called a multi-week
> commitment was run, after a library fix (`monitor_size` pin) made it a
> 25–32 s/count sweep instead. Result: the **hybrid's own error against the
> FDTD reference** (a different question from Residual 1's EME-vs-hybrid
> disagreement, which is untouched and still open) splits into two —
> **blue half (1.50–1.55 µm) is now solved** by raising `interface_modes`
> to 22–24 (worst per-port error 5.99%→0.62%, mechanism confirmed: the
> single TE input excites only even lateral orders, wasting most of a
> 12-mode basis); **red half (1.56–1.60 µm) is untouched**, physically
> insensitive to mode count. `body_stub` cut position was independently
> ruled out as a lever (0.07% effect, all 4 ports, gated at 0.12% against
> the shipped result). A C-band resim at 21 points confirms the FDTD
> reference has no hidden pathology and resolves real beat structure the
> coarse grid was aliasing by luck; the matching hybrid resim is still
> running.
>
> **What is closed.** Steps 1–7 (nitride Sellmeier entry, `MaterialStack`
> core/clad generalization, cache fingerprint, cross-section design,
> body design, two example scripts, wiring) and **all four checks**:
> check 1 measured across the full band, check 2 (mode identity on a 2-mode
> access waveguide), check 3 (resolution convergence, job 12275175),
> check 4 (whole-device FDTD reference at band centre *and* full band,
> jobs 12029541 / 12260631). Plus step 9a (β-calibration — hypothesis
> **refuted**, it is not the split-shape explanation) and the
> `decay_method="energy"` study.
>
> ### ⚠ Residual 1 — check 1's engine disagreement is REAL and unexplained
>
> Not a constant offset: FDTD and EME have **opposite spectral slopes**
> (−0.749/µm vs +0.087/µm), and band centre is the crossover, so every
> band-centre-only comparison flattered the agreement. Concretely, EME vs
> hybrid: split shape disagrees **20× with opposite sign** (EME favours the
> outer pair by 0.0008, hybrid the inner by ~0.017) and **total T differs by
> 2.5 points** (0.9707 vs 0.9454) — larger than either engine's own split
> asymmetry.
>
> Mode truncation was tested and **REFUTED** as the cause (a free 3-engine
> ladder: *more* truncation correlated with *less* disagreement — the wrong
> direction). The confound check was also done properly first: the gap
> survives matching EME's resolution and domain to FDTD's own (jobs
> 12323220 / 12331012).
>
> This residual is **larger than the precedents already sitting in
> `resolved/`** (the 2×2's 0.47% Δβ, the coupler's ~5–10%). It is recorded,
> not closed.
>
> ### ⚠ Residual 2 — the leading suspect was DROPPED, not investigated
>
> The remaining candidate for residual 1 is the **EME taper staircase**
> (`C_TAPER=8` is the only path approximating the port tapers; FDTD meshes
> them exactly). The convergence sweep that would test it was **never run,
> at the user's instruction**. "Dropped" is not "excluded" — if this device
> is ever picked back up, that sweep is the first thing to run.
>
> ### ⚠ Residual 3 — res=30 is NOT established for the EME side
>
> Check 3 found the real mechanism (a common-mode n_eff shift that does not
> propagate to the differential FDTD depends on for self-imaging) and
> concluded **res=30 defensible for FDTD, NOT established for EME** —
> EME's differential drift 30→40 is 2.17e-03, an order of magnitude worse
> than FDTD's 1.99e-04. Note `mmi_1x4_eme.py` still runs at
> **`RESOLUTION = 40.0`** while `_fdtd.py`/`_hybrid.py` moved to 30 on
> 2026-08-18. That split may well be correct given the above, but it is
> **undocumented** — decide and record it before relying on either script.
>
> ### Open to-dos (do not block resolution)
>
> - **Nitride strip-test gate — substantively PASSES, harness needs a call.**
>   Run 2026-08-31 (`verify_item_8_strip_gate.py`, extended to non-silicon
>   cores; SIN_300 C-band added at the device's own 1.1 × 0.30 µm). Arm B
>   (the re-parameterisation under test) lands **Δn_eff = +0.000207**, ten
>   times inside the 0.002 tolerance, with ε∞ smoothed to 4.0099 against a
>   true 3.9851 — so ε∞ *is* carrying the contrast. What did not clear is the
>   script's **harness self-check**, which demands arm C (the known-bad ε∞=1
>   configuration) move by more than an absolute `10 × TOL = 0.02`; nitride's
>   C moved −0.017311. That bar was calibrated on silicon's 0.04–0.11
>   failures, and a low-contrast core *must* fail smaller — the error is the
>   contrast destroyed, and nitride loses 3.99→1 against a 2.085 cladding
>   where silicon loses 12.09→1. By ratio |ΔC|/|ΔB| nitride is **84×**,
>   mid-pack among the silicon rows (73× / 99× / 151× / 395×). **Threshold
>   deliberately NOT changed** — re-calibrating a gate after seeing it fail is
>   the wrong reflex; it needs a decision, not a patch.
> - ~~Step 9d (`interface_modes` sweep) — scoped, not run~~ **RUN
>   2026-09-05 → 2026-09-13, see "✅ 9d RESULT" below.** The "~3-day prewarm
>   per mode count" cost was real for the EME-side truncation question this
>   header refuted, but for the hybrid's own error against FDTD, a
>   `monitor_size` pin made the sweep 25–32 s/count once a 24-mode ceiling
>   was prewarmed once.
> - Step 9b (output-side mirror reuse) — proposed; needs a cache-keyed
>   `MirrorSpec` variant (`beef.MirrorSpec` is x-mirror; this device needs
>   transverse-across-`solve()`).
> - Four scripts carrying pre-2026-08-14 constants were refreshed
>   2026-08-31 (commit `ca42ea9`). Outputs in `example_output/` predating
>   that date were produced with the old geometry; each script now says so.

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
| Example scripts | ✅ `mmi_1x4_eme.py` (job 11887395) and `mmi_1x4_hybrid.py` (job 12227591, calibrated) both **solved for real**. `mmi_1x4_fdtd.py` was dropped 2026-08-04 because `mm_length=206` made whole-device FDTD impractical — **that reason expired** with step 9c (216 µm → 49.90 µm, 4.3×); script revived 2026-08-19, and the whole-device reference is now **MEASURED at band centre** (job 12029541, 41 hrs, production `decay_tol`). Full-band run still pending. See "CHECK 4" under step 9. |

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

> **⚠ This justification is stale as of 2026-08-20 and the decision should be
> re-taken, not inherited.** `MM_LENGTH` became **39.90 µm** on 2026-08-14
> (`examples/mmi_1x4/python/mmi_1x4_hybrid.py:157`), so the total device length is
> now ≈50 µm rather than ≈216 µm — the reasoning above is off by more than 4×, and
> a whole-device FDTD reference (check 4) may well be affordable. The 2×2's own
> full-3-D runs at comparable scale take 1.7–2.5 h each
> (`running_jobs.md:116-121`). Re-ask the cost question before concluding this
> device ships EME + hybrid only.
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

### Step 9 — two cost/correctness items opened 2026-08-13, neither run yet

Both found while reviewing the 6d split-shape discrepancy (EME favors inner
ports 1,2; hybrid favors outer ports 0,3 — still unexplained, step 6d/check 1).

**9a. `calibrate_body_neff` — ✅ DONE 2026-08-23. Ran to completion, and the
hypothesis it was testing is REFUTED — see "9a RESULT" below. The narrative
between here and there is the three-attempt history, kept because each
cancellation was for a real reason.** `mmi_1x4_hybrid.py`
carried the body on tidy3d's own mode index the whole way through the
11734486 run recorded in step 6d; `calibrate_body_neff=True` (HybridSolver's
built-in β-calibration, `dc_validation.md` item 6) was flagged in the
script's own docstring as "available, unused by default here... if step 8's
engine-agreement check finds it matters" but never turned on. Turned on this
session: `CALIBRATE_BODY = True` module constant, threaded into
`HybridSolver(calibrate_body_neff=CALIBRATE_BODY)`.

Why now, not deferred further: this device pairs tidy3d (body) with meep
(input/output FDTD) — a strictly bigger cross-engine basis gap than
`directional_coupler_hybrid.py`'s own meep-vs-MPB case that
`calibrate_body_neff` was built to fix, and the coupler's version of this
exact mismatch was a measured 10.4% error on the coupling-sensitive quantity
(Δβ) before calibration. A plausible, not yet confirmed, explanation for the
6d split-shape sign flip.

**Cost, not yet paid.** `calibrate_body_neff` adds `len(calibrate_lengths)`
(default 2) extra FDTD section solves on the body span, sourced from the
input section's `section_mode_beta()`, cached and shared across all four
`output_idx` calls — same mechanism `directional_coupler_hybrid.py` calls
"three short FDTD section solves, cheap" in its own docstring. **Not cheap
here**: this device's own measured FDTD section cost is ~12–24 hrs/section
(`decay_tol=1e-11` on weakly-confined 300 nm nitride, step 6d), so budget for
that class of job, not the coupler's. `calibrate_body_neff` is part of
`HybridSolver`'s own section cache key, so turning it on invalidates none of
the *existing* (uncalibrated) cache entries — it adds new ones alongside
them — but every number in step 6d's table (T, IL, the split-shape finding
itself) was produced WITHOUT this flag and must be re-measured under it, not
assumed to still hold or to explain the discrepancy until an actual run says
so.

**9a, continued — first real attempt cancelled 2026-08-17, 3 days in, NOT a
failure.** The real 4-task prewarm array (job 11856080, then resubmitted
under `gts-sr108` as job 11859048 after the earlier account-headroom issue)
was running the OLD `section_mode_beta()` calibration path the whole time —
`len(calibrate_lengths)` (2) × `N_MODES` (12) = 24 fully independent FDTD
domain builds, each decomposing the FULL 12×12 transfer matrix (not just the
diagonal `T_ii` the calibration actually needs). Per
`validation_markdown/resolved/section_mode_beta_multi_monitor_validation.md`
step 9 (landed in `beef/solvers/fdtd.py` 2026-08-17 12:00:02, confirmed by
file `mtime`, three days after job 11859048 started) — `measure_diagonal_only`
(~4×) and a cache-key fix that finally let the merged-domain path
(~1.9× more) activate for real production callers, combined ~8.3× on a
comparable device (MMI 2×2: 30.6h calibration → ~3.7h). A running Python
process keeps whatever code it imported at start, so 11859048 could not pick
this up mid-flight — **cancelled** (user's explicit call after confirming
the timeline) rather than let it run several more days on the slow path.

**Resubmitted, 2026-08-17, same job otherwise unchanged**
(`HybridSolver(calibrate_body_neff=True)` needs no code change per that
doc's own "how to apply to another device" section) — smoke-tested first
(job 11938126, reduced scale, 49m41s, clean) before the real array
(job **11939350**, `-t 2-12:00:00`/60h, down from the cancelled run's blind
96h guess — new cost estimate ~43-48h: input/output FDTD sections unchanged
at ~41-45h since `measure_diagonal_only` only touches the calibration's own
diagonal-only need, plus ~1.5-3h for the now-fast calibration itself).

**Correction: default calibration lengths is 3, not 2.** Repeated
throughout this doc's own step 9a writeup above and in several batch-file
comments — `section_mode_beta()`'s own docstring says "Default is three,
evenly spaced," and job 11982321's own real run confirmed it:
`lengths=[1.789, 3.577, 5.366] µm` auto-selected. Not corrected retroactively
everywhere above; noted here so a future reader doesn't propagate the error
further.

**9a, second cancellation, 2026-08-17 — job 11939350's own estimate was
also wrong, for a different reason: `measure_diagonal_only` doesn't touch
the FDTD *stepping* cost, only decomposition, and stepping is what
dominates on this device.** At 21h23m in, task 0's own log showed
calibration mode 1 of 12 needing 12,706 timesteps (~7.3 hrs) to decay —
extrapolating linearly across all 12 modes projected ~85-100h for
calibration alone, far past the 60h budget. Compounding the risk:
`solve_section_multi_length()` (`beef/solvers/fdtd.py:3260`) has **no
incremental caching** — `save_run()` sits after the full per-mode loop
(line 3548, loop starts line 3473) — so a wall-time kill mid-calibration
discards the entire attempt, not just the unfinished modes.
`scontrol update JobId=... TimeLimit=...` on the running job was refused
(permission denied) — no way to extend in place. **Cancelled again**
(user's explicit call after reviewing the projection), rebudgeted to
`-t 6-00:00:00` (144h) from the real per-mode number, and resubmitted as
job **11981595** — immediately hit `AssocGrpBillingMinutes` on `gts-sr108`
(climbed to 0.706 usage from this session's own accumulated work in the
few days since the last account switch), moved to `gts-sr108-gedc_main`
(0.065 usage, the account this exact calibration doc's own step 8 needed
for an analogous large MMI 2×2 job) as job **11981595** for real this time.

**Resolution dropped 40→30 before letting 11981595 run to completion,
measured first.** Rule-of-thumb check: `n_SiN(1.55µm)=1.996` →
23.3 points/wavelength in the core at res=30, comfortably above the 15-20
PPW "good accuracy" bar (vs. 31.1 PPW at res=40). The only prior in-repo
resolution check for this device (`mmi_1x4_sin_body_res_convergence.py`,
step 8 check 3) tested the EME/tidy3d side only, inconclusively — meep/FDTD's
own resolution convergence had never been checked at all. Cancelled
11981595 (only ~23s in, negligible loss) and ran two parallel timing probes
before committing the real run again: `investigate_mmi1x4_res30_input_
timing.py` (job **11982320**, real production input-section solve, cached
and reusable) and `investigate_mmi1x4_res30_calibration_timing.py` (job
**11982321**, one calibration mode, n_modes=1 probe, not cached).

**MEASURED, and the two section types respond differently to the same
resolution drop:**
- Calibration (uniform straight-waveguide domain): genuine **~2.1-2.5×
  speedup** — 3.39 hrs/mode at res=30 vs. 7.3 hrs/mode at res=40.
  `lengths=[1.789, 3.577, 5.366] µm` auto-selected; fitted β ranged
  6.55–7.06 rad/µm across the C-band, a sane dispersion trend.
- Input/output taper+stub sections: **no measured benefit** — 11.49 hrs at
  res=30, essentially unchanged from the ~11-12 hr res=40 baseline. Real
  cause, not assumed: res=30 needed ~1.4× MORE timesteps to reach
  `decay_tol=1e-11` (23,388 steps/t=389.8 vs. res=40's own ~17,000
  steps/t≈213) than the naive Courant-limited-Δt argument predicted, which
  roughly cancels the ~2× per-step speedup the coarser grid does deliver.

Net: still a real, worthwhile overall win (calibration was already the
majority of total cost), just not the ~3.2× a naive `(30/40)^4` grid-cell
argument would predict — closer to **~2×** once the input/output sections'
own non-response is folded in. `RESOLUTION` changed to `30.0` in both
`mmi_1x4_hybrid.py` and `_prewarm_mmi_1x4_hybrid.py` (must match exactly —
same cache key). Revised per-branch projection: input (11.49h, measured) +
calibration (12×3.39h = 40.7h, still a 1-mode extrapolation — per-mode
variance across the other 11 modes unmeasured) + output (~11.5h, ASSUMED
similar to input, same geometry class, itself unmeasured) ≈ **64h**.

**Real run resubmitted at res=30, job 12014403** (`-t 5-00:00:00`/120h, ~1.9×
margin over the 64h projection), account `gts-sr108-gedc_main`. No fresh
smoke test — the two timing probes themselves already exercised this exact
config (real res=30, real `decay_tol=1e-11`) to clean completion, a more
rigorous check than a reduced-scale smoke test would have been.

**Cancelled AGAIN at ~8h in, 2026-08-19 — not a failure, a real further
optimization found while it ran.** User asked directly whether the
transverse (y) and vertical (z) sim domain could shrink — `sim_width =
MM_WIDTH + 8.0 = 20.0` and `sim_height = SIN_300.sz = 5.0` were both
inherited conventions (shared with every other FDTD section on this
device), never validated for necessity. Checked cheaply (mode-solve only,
no FDTD): `mmi_1x4_sin_body_domain_size_convergence.py` (job **12029224**,
4m20s) swept both on the body cross-section, tracking the same modes
`mmi_1x4_sin_body_res_convergence.py` already used. Key realization: modes
19/22 in that tracked set never converge within the swept range, but they
are **irrelevant** — production only keeps `N_MODES=12` (indices 0–11), and
mode 10 (the second-to-last kept mode) converges cleanly by
`mesh_width=13–14` and `sz=4.0`. `sim_width` dropped to `MM_WIDTH+4.0=16.0`
(mode 19's own convergence point, real margin over what mode 10/11 need)
and `sim_height` to `4.0`, decoupled from `SIN_300.sz` (a shared
`MaterialStack` constant other solvers/devices also read — not itself
changed). ~36% fewer grid cells (`16/20 × 4/5 = 0.64×`) on top of the
res=40→30 change, applied to `mmi_1x4_hybrid.py`'s `_fdtd_common` AND its
EME `body_solver`'s own `oxide_width`/`oxide_height` (same tidy3d mesh the
convergence check directly measured), and to `_prewarm_mmi_1x4_hybrid.py`
(must match exactly for the cache). Re-smoke-tested (job **12029535**)
before resubmitting the real array. Only ~8h of the ~64h projection was
sunk at cancellation — a worthwhile trade against a further ~36%+ cut to
the remaining ~56h.

**9a, the run that is actually going: job 12031838 (4-task array),
submitted 2026-08-19. ✅ ALL 12 CALIBRATION MODES COMPLETE on all four
branches as of 2026-08-22 08:26** (2d 14h elapsed, 2d 9h walltime left).
This is the first time this device's β calibration has run to completion —
the two earlier attempts were cancelled at 3d and ~8h for the reasons above,
neither a failure.

Progress was read from the logs' **domain-size fingerprints** rather than
extrapolated from step counts, because the per-mode cost varies too much for
extrapolation to be trustworthy (see the spread below):

| domain (`Computational cell is …`) | count | what it is |
|---|---|---|
| `6.36667 x 17 x 5` | **12 / 12** | the calibration body-section domains, one per excited band — all built and time-stepped |
| `6.76667 x 17 x 5` | 2 | input / output sections (vertices name y = 5.2 µm = the outer output port at 4.65 + 1.1/2) |
| `3 x 17 x 5` | 1 | short section |

Task 0's completed FDTD run lengths, in order: **23387** (input section —
matching the 23,388 measured independently for this section) then
11465 / 11464 / 11461 / 11466 / 25932 / 9555 / 10099 / 24201 / 15832 /
24751 / 10553 for the calibration modes. **That is a 2.4× spread (9.5k–25.9k
steps) across modes of the same device at the same config** — worth
recording, because it is exactly why the walltime projections for this step
kept missing: an average taken over the first few modes does not predict the
rest. A mid-run estimate of "15–20 h of calibration remaining", made from 7
completed runs averaging 3.4 h, came in pessimistic — the remainder finished
inside ~17 h.

**✅ PREWARM COMPLETE 2026-08-22.** All four tasks `COMPLETED`, exit `0:0`,
elapsed 2d23h49m – 3d02h54m. `Branch {0,1,2,3} done -- cache warmed`.

**✅ REAL RUN 2026-08-23, job 12227591 — 49.8 s against the warm cache**
(3 days of prewarm → a 50-second production run, mirroring
`mmi_1x4_eme.py`'s own 38 s). Config was verified field-by-field against
`_prewarm_mmi_1x4_hybrid.py` before submission (res=30, decay_tol=1e-11,
dpml=0.5, sim 16×4, N_MODES=12, CALIBRATE_BODY=True, body_stub=0.25,
identical geometry) — a mismatch on any cache-relevant field would have
silently turned this into a multi-day run instead of a cache hit.

> ⏱ **14 h lost to an account problem, not a technical one.**
> `batch/run_mmi1x4_hybrid_real.batch` hardcoded `--account=gts-sr108-ece`,
> which is **99.94% through its `GrpTRESMins` billing cap** (203,719 left of
> 354,561,109), so Slurm parked the job on `AssocGrpBillingMinutes` overnight.
> **`sshare` does not show this** — it reports fairshare (`EffectvUsage`),
> a different mechanism; `gedc_main` reads 0.064 there and *is* the right
> answer, but for an unrelated reason. The billing cap is only visible via
> `scontrol show assoc_mgr flags=assoc account=<acct>`. Batch file switched to
> `gts-sr108-gedc_main` (no billing cap) with the numbers recorded inline.

### 9a RESULT — calibration is NOT the split-shape explanation ❌ hypothesis refuted

At λ=1.55 µm, all three at the **same current geometry** (39.90 / 3.1):

| run | outer (0,3) | inner (1,2) | favours | spread | total T |
|---|---|---|---|---|---|
| EME, job 11887395 | 0.2431 | 0.2423 | **outer** | 0.0008 | 0.9707 |
| Hybrid **uncalibrated**, job 12227682 | 0.2287 | 0.2460 | **inner** | 0.0173 | 0.9495 |
| Hybrid **calibrated**, job 12227591 | 0.2280 | 0.2447 | **inner** | 0.0167 | 0.9454 |

**Calibration moves essentially nothing**: per-port ΔT ≈ 0.0007–0.0013
(0.3–0.5%), total T 0.9495 → 0.9454, and the inner-favouring split shape is
**already there with calibration off**. So the step-6d hypothesis — that
calibration explains the EME-vs-hybrid split-shape disagreement — is
**refuted**, not confirmed.

**Why this needed a deliberate A/B rather than reading it off the calibrated
run.** Calibration was enabled 2026-08-13 and the geometry replaced
2026-08-14, one day apart, so **no uncalibrated hybrid run existed at the
current geometry** (verified: only two runs at this geometry exist at all).
Comparing the new calibrated result against step 6d's old-geometry numbers
would have confounded the two changes, and would have made calibration look
like it flipped the split shape when it does not. The isolating run was cheap
— **9 minutes** — because the FDTD input/output section cache does not depend
on calibration at all (only the EME body's `neff_override` does), so it reused
the 3-day prewarm. `CALIBRATE_BODY` is now env-overridable
(`MMI1X4_CALIBRATE_BODY=0`), **default unchanged**.

**What remains unexplained, and is now the real open item:** EME and the
hybrid still disagree, and not subtly — EME favours the OUTER pair by 0.0008
while the hybrid favours the INNER pair by ~0.017, a **20× larger spread of
the opposite sign**. Separately, **total T differs by 2.5 points**
(0.9707 vs 0.9454), which is larger than either engine's own split asymmetry
and is a check-1 engine-agreement finding in its own right. Neither is
attributable to β calibration any more.

**A second, unexpected calibration effect worth keeping.** Calibration barely
moves the answer *at* MM_LENGTH=39.9 but substantially changes the
**length-sweep landscape**: the best-balance optimum reads 246.9 µm
(spread 0.0001) uncalibrated vs 129.4 µm (spread 0.0000) calibrated. So
calibration is not inert in general — it is inert *at this short body*. A
plausible (NOT verified) reading is that accumulated β error scales with body
length, so the 39.90 µm body simply has too little length for the correction
to matter, where the superseded 206.0 µm body would have had ~5×
more. That would also predict calibration mattering much more for the old
geometry than the new one — testable, not tested.

**Artifacts.** Both result sets are preserved side by side rather than
overwritten (the two runs share an output directory):
`example_output/mmi_1x4/mmi_1x4_hybrid_CALIBRATED_12227591/` and
`..._UNCALIBRATED_12227682/`; the canonical
`example_output/mmi_1x4/mmi_1x4_hybrid/` holds the calibrated set.

### 9d — the engine gap is a MODE-COUNT mismatch, and sweeping it is NOT cheap

**The two engines are not being asked the same question.**
`mmi_1x4_eme.py` carries `N_BODY=20` body modes; `mmi_1x4_hybrid.py`
truncates every interface to `interface_modes=12` (cut 20→12 for
tractability, with this doc's own N_MODES comment already conceding it
"genuinely trades some self-imaging accuracy"). The hybrid — the truncated
one — reports the LOWER total T (0.9454 vs 0.9707), which is the direction
truncation predicts. So the check-1 disagreement is very likely a check-3
interface-completeness question, not an engine bug.

**Cost probe, MEASURED 2026-08-23 (job 12260083): a mode-count sweep is
NOT free.** `FDTDSolver.solve_section()` does implement superset cache reuse
(a cached 12-mode solve slicing down to serve a smaller request), so a
DOWNWARD sweep looked like it might cost nothing. It does not apply here:
run at `N_MODES=10` against the warm 12-mode cache, the job **started
time-stepping within minutes** (57 timestep lines at 5m39s) rather than
hitting the cache. Killed at that point — the probe was deliberately capped
at `-t 12` so a miss cost nothing.

The likely reason is the documented precondition (`fdtd.py`'s own
family_config comment): superset reuse requires a **pinned `monitor_size`**,
because the default monitor width depends on whether `n_measure > 1` — so
without pinning, the *simulated monitor geometry* changes with mode count,
not merely the number of bands decoded from it, and the two runs are not
physically the same simulation. This chain does not pin it. **Not
root-caused further** — the probe answered the cost question, which is what
it was for.

**Consequence: each mode count needs its own ~3-day prewarm.** An
`interface_modes` convergence sweep — the thing
`beef/interface_modes_convergence.py` exists for, and the only way to
establish basis completeness (layer 2) — is therefore a multi-week
commitment at this device's FDTD cost, not an afternoon. That is a
scheduling decision, not a technical one.

`N_MODES` is now env-overridable (`MMI1X4_N_MODES`), **default unchanged**,
so the sweep can be driven without editing the script when it is run.

### ✅ 9d RESULT — the sweep was finally run, 2026-09-05 → 2026-09-13

**⚠ Scope, read first.** This settles the **hybrid's own error against the
whole-device FDTD reference**. It does **not** touch Residual 1 (the
EME-vs-hybrid split-shape/total-T disagreement in the RESOLVED header
above) — that comparison's own truncation test (the free 3-engine ladder)
stands unchanged and unrefuted. The two questions are orthogonal: this
section is about interface completeness on the **hybrid** side only.

**The blocker was the `monitor_size` pin, not compute.** Added
`monitor_center=0.0`/`monitor_size=FDTD_SIM_WIDTH` to sections 0 and 2 (never
the body — `_calibrated_neff` forwards those exact keys out of the body's own
kwargs into the β calibration, so pinning it there would rehash and re-pay
that calibration too) in both `mmi_1x4_hybrid.py` and
`_prewarm_mmi_1x4_hybrid.py` (`MMI1X4_PIN_MONITOR`, default on). This is what
`solve_section()`'s own `family_config` gate needs (`monitor_size is not
None`) to build the mode-count "superset" cache family at all — without it,
job 12260083 and a repeat this session (12797780) both confirmed a downward
probe starts time-stepping instead of slicing. With the pin, a 24-mode
prewarm (jobs 12963023/13019844/13020152/13019542, 4 branches, ~1 day each
after retries — see Lustre note below) made every smaller count from 8 up
free: the full ladder 8/10/12/14/16/18/20/22/24 (job 13058948) ran in
**25–32 seconds per mode count**, verified bit-identical (0.000% deviation)
against an independently-computed 12-mode run.

**Mechanism, confirmed directly from cached per-mode power, not inferred.**
The single centred input excites TE only, and within TE only **even**
lateral orders — odd orders are antisymmetric, zero by parity; TM modes are
zero by polarization mismatch. At 16 modes the power-carrying indices are
exactly `[0, 2, 4, 6, 8, 15]` (0.357/0.291/0.191/0.099/0.040/0.013);
`interface_modes=12` was spending 7 of 12 slots on modes carrying ~1e-11.
This also refutes the standing prediction it inherited: item 9d's own
"more truncation ⇒ less disagreement" framing (from the EME-side ladder)
does NOT hold for hybrid-vs-FDTD — here, less truncation strictly helps.

**Convergence — two real plateaus, at irregular counts, not every 4.**
A two-family (TE + TM) n_eff ladder fit off the 16-mode data correctly
predicted TE10 at index 15 and TE12 at index 21 (power 0.0034) — confirmed
directly, not just fit. 18–20 and 22–24 are bit-identical to each other
because no new mode exists in that gap:

| interface_modes | mean per-port err vs FDTD | worst | blue-band worst (≤1.55 µm) |
|---|---|---|---|
| 12 (production, calibrated) | 2.78% | 5.99% | 5.99% |
| 12, uncalibrated | 2.98% | 5.89% | 5.89% |
| 16 | 1.65% | 4.40% | 1.59% |
| 18–20 | 0.99% | 2.64% | 1.59% |
| **22–24** | **0.69%** | **2.39%** | **0.62%** |

**The error is now two separable problems, one solved.** Blue half
(1.50–1.55 µm): fixed by mode count, under the 1% target at 22–24 modes
(0.62% worst, down from 5.89–5.99%). Red half (1.56–1.60 µm): **completely
insensitive to mode count** — 16/18/20/22/24 are bit-identical there because
the modes that would help (TE10 at 1.56 µm) are physically cut off past that
wavelength. That residual (3.1–4.4%) needs a different mechanism, plausibly
connected to Residual 1's own opposite-spectral-slope finding — not
re-investigated here.

**Calibration re-examined across the full band, not just band centre —
result sharpens, doesn't reverse, 9a.** Uncalibrated has the smaller
worst-per-port error at 8 of 11 wavelengths; calibrated wins only at the red
end (1.58–1.60) by a wide enough margin to still win on the *mean* (2.78%
vs 2.98%). Split shape is essentially unmoved by calibration either way
(~0.0166–0.0173 vs FDTD's ~0.0025 at band centre, calibrated or not). Net:
calibration is not earning its keep over most of the band; kept as the
default pending a decision, not because the data supports it.

**The `body_stub` cut-position axis — ruled out, all four ports, gate
passed.** New script `mmi_1x4_body_stub_sweep.py`
(`examples/verification/python/`) walks the interface cut 0.25→1.00 µm using
`measure_margins` (one FDTD run per section puts all 6 depths on separate
monitors in the same simulation — 5 runs total, not 5×6) and reconstructs
the cascade with the calibrated β reused across all depths (the body's
per-mode β is a property of the uniform 12 µm slab, identical at every cut
depth — only the propagated length changes). **Gate**: production
`body_stub=0.25` reconstructed to `[0.22832, 0.24497, 0.24496, 0.22833]`
vs the shipped `[0.22804, 0.24468, 0.24467, 0.22805]` (job 12227591), worst
deviation **0.12%** against a 2% bar — validates both the reconstruction and,
independently, `fdtd.py`'s measure-side `PortSpec` LOW/HIGH convention,
which the library's own docstring flagged as unverified until this run.
**Result: dead lever.** Across the full 0.25→1.00 µm range, total T varies
0.0097% and split shape 0.0686% — against the then-5.9% error being chased.
Not the mechanism.

**A C-band (1.530–1.565 µm, 21 pt, 1.75 nm spacing) resim was run to pin the
TE10 cutoff precisely** — the original 11-point 1.50–1.60 grid (10 nm
spacing) straddled it by luck, not resolution. `MMI1X4_FDTD_WL_MIN/MAX/NWL`
and matching `MMI1X4_HPREWARM_WL_MIN/MAX` added (both default to the
original 1.50/1.60/11, additive). **Whole-device FDTD reference: ✅
COMPLETE (job 13092669, 2d13h)** — smooth, monotonic total T through the
whole window (no hidden pathology the coarse grid was aliasing), but real
oscillatory structure in split shape: a genuine bump–dip–bump beat pattern
(local peak ≈1.545 µm, dip ≈1.55 µm, global peak ≈1.5597 µm) that the
11-point grid happened to sample close enough to both peaks and the
intervening dip to look smooth — a different phase offset could have missed
it. **Hybrid-side resim: still running** (job 13169522, 24 modes,
uncalibrated, resubmitted with a 6-day budget after a first attempt
(13092670) TIMED OUT at 2d12h with **nothing cached** — confirmed via a
direct `cache_only=True` probe, `CacheMissError` on all 4 branches, not
assumed; the true per-branch cost at this mode count/wavelength combination
is still unmeasured). Update this row once it lands.

**Two operational findings worth carrying forward, independent of this
device.**

1. `solve_section()`'s cache write is all-or-nothing — `save_run()` fires
   once, only after the full result array is assembled and returned. There
   is no incremental persistence. A crash or walltime kill at 99% loses the
   **entire** call, not a resumable fraction. This cost two real multi-day
   losses this session (see next point, and the C-band hybrid TIMEOUT
   above) and should be assumed by default for any long `solve_section()`
   call, not rediscovered per device.
2. A recurring `os.getcwd() FileNotFoundError` inside meow's import
   (`meow_MEEP/src/meow/materials.py`'s `_validate_path`), confirmed via
   `dmesg` as a real Lustre client error, hit this campaign **3 separate
   times** across 2 different job families, on **7 different compute
   nodes** — not a single bad node. The 24-mode uncalibrated prewarm
   recovered via per-branch resubmits (no compute lost — each failure hit
   before that branch's own expensive calibration/output stage). The
   16-mode **calibrated** prewarm was not so lucky: it failed **twice**,
   both times after the entire ~16-run diagonal-only calibration had
   completed and was mid-decomposition on the very last band — i.e., at
   effectively 99% done — losing ~3 days of compute each time, with
   nothing recoverable per finding 1. **That confirmation was dropped**
   (not re-attempted a third time) rather than accept a third multi-day
   exposure to the same all-or-nothing failure mode for a result that was
   already trending "mildly harmful," consistent with 9a's own finding at
   12 modes. Neither instance was escalated to PACE support; if this
   recurs on a future long single-rank meep/meow job, escalate rather than
   retry a third time blind.

### ✅ CHECK 3 — resolution convergence, RESOLVED 2026-08-24 (job 12275175)

**Both the original finding and its leading explanation were wrong, in
different ways.** Full sweep res=[15,20,25,30,35,40], 16 modes solved,
12 reported (= production `interface_modes`), BOTH engines, mode identity
established by **field overlap** rather than by index. 3h29m. Script:
`examples/verification/python/mmi_1x4_check3_res_convergence_corresponded.py`.

**1. The mode-swapping hypothesis is REFUTED.** The doc's standing explanation
for the 2026-08-04 inconclusive result was "mode-index correspondence swapping
between near-degenerate higher-order modes." It is not that: **correspondence
is the identity permutation at every resolution, in both engines**, worst
diagonal overlap 0.9972-1.0000 on the reported modes. The naive index-tracked
and field-overlap-corrected |dn_eff| series are consequently **identical to
every printed digit**. The original attempt tracked modes correctly after all;
its result was real.

**2. The non-convergence is a COMMON-MODE shift, not a per-mode or
near-cutoff problem.** Per-mode |dn_eff| (FDTD, modes 0->11) is essentially the
same number for every mode at each step -- 15->20 is 3.4e-03 for mode 0 and
3.2e-03 for mode 11; 30->35 is 1.1e-03 vs 9.9e-04. **The fundamental mode is
just as non-monotonic as the near-cutoff ones.** That rules out the
"weakly-confined high-order modes need more mesh" reading that motivated
targeting the body cross-section in the first place, and points instead at a
mode-independent discretisation bias on the effective index.

**3. THE QUANTITY THE DEVICE DEPENDS ON DOES CONVERGE -- on the FDTD side.**
This is what step 8's own instruction ("sweep a **Δβ-like quantity**, not just
output power") was pointing at, and what the first attempt and my own first
pass both missed by tracking absolute n_eff. Because the error is common-mode,
it largely cancels in differences:

| | absolute n_eff drift 30->40 | **differential** Δn_eff drift 30->40 |
|---|---|---|
| FDTD / meep | 1.96e-03 | **1.99e-04 max, 8.14e-05 mean — 10x smaller** |
| EME / tidy3d | 1.35e-03 | 2.17e-03 max, 6.10e-04 mean — *not* smaller |

FDTD's differential also shows a clean asymptotic tail once past the coarse
points: max |ΔΔn_eff| runs 3.2e-03 -> 6.2e-04 -> 1.5e-04 -> **4.6e-05**,
shrinking ~3-4x per step from res 25 up. And the mode-0<->1 spacing, the most
self-imaging-relevant single number, agrees to ~0.1% between res 30 and 40 in
**both** engines (FDTD -0.003563 vs -0.003559; EME -0.003553 vs -0.003555).

**Verdict: res=30 is defensible for the FDTD side** on the differential
quantity that governs MMI self-imaging, residual ~2e-04 max / 8e-05 mean
against res=40. The absolute-n_eff offset is a common-mode discretisation bias
that does not propagate to Δβ.

> ⚠ **Two caveats, both real.**
> **(a) The EME/tidy3d side does NOT show the same cancellation** -- its
> differential drift (2.17e-03) is no smaller than its absolute drift, and its
> per-step sequence stays noisy. res=30 is justified for the meep side but
> **not** established for the tidy3d body solver. Worth connecting to check 1:
> EME is the engine that got the split-shape *sign* wrong against the
> whole-device FDTD reference, and split shape is a Δβ-across-many-modes
> quantity. **Hypothesis, NOT established** -- an under-converged EME
> differential could contribute, and it is testable by re-running the EME body
> at res 40+ to see whether the split shape moves.
> **(b) My script's `monotonically shrinking` flag is too crude** and prints
> `False` for both engines. It demands *every* consecutive step shrink, which
> a sequence with one coarse-grid spike fails even when its tail converges
> cleanly -- exactly FDTD's case. The per-mode/per-step tables are the
> evidence; do not read that boolean as the finding.

---

### ✅ EME converges — to a DIFFERENT answer than FDTD (jobs 12323220, 12331012)

Full (res, domain) grid at band centre. **The split shape never flips.**

| case | res | domain | total T | favours | \|outer−inner\| | dev vs FDTD |
|---|---|---|---|---|---|---|
| A | 40 | 20x5 | 0.9707 | OUTER | 7.93e-04 | 3.09e-03 |
| B | 30 | 16x4 | 0.9703 | OUTER | 2.38e-03 | 3.79e-03 |
| C | 40 | 16x4 | 0.9706 | OUTER | 2.31e-03 | 3.84e-03 |
| D | 50 | 20x5 | 0.9705 | OUTER | 7.62e-04 | 3.04e-03 |
| E | **60** | 20x5 | 0.9707 | **OUTER** | **7.48e-04** | **3.07e-03** |

**1. EME IS resolution-converged.** At the 20x5 domain, |outer−inner| settles
7.93e-04 → 7.62e-04 → **7.48e-04** across res 40 → 50 → 60, and total T holds
0.9705–0.9707. This is the important negative result for check 3's own
hypothesis: EME's *cross-section* differential non-convergence (2.17e-03 drift,
30→40) **does NOT propagate to the device-level split shape**, which converges
cleanly. That hypothesis is therefore **not supported** — dropped, not carried
forward.

**2. Deviation from FDTD does NOT shrink with resolution.** dev vs the
whole-device FDTD reference sits at 3.04–3.84e-03 across *every* case, and is
3.07e-03 at res=60 versus 3.09e-03 at res=40. **EME is converged and
disagrees.** So check 1's engine disagreement is real physics/modelling, not a
mesh artefact — the strongest statement available on it so far.

**3. DOMAIN, not resolution, drives the split magnitude.** 20x5 gives
|outer−inner| ≈ 7.5e-04; 16x4 gives ≈ 2.3e-03 — the smaller domain roughly
**triples** it (B and C agree at 2.3–2.4e-03 across res 30 and 40, so this is
a domain effect, not a resolution one). This confirms the suspicion that the
16x4 domain — validated for the hybrid's **12** interface modes — is **not
adequate for `mmi_1x4_eme.py`'s N_BODY=20**. The published EME run's 20x5 is
the right choice, and is the converged one.

---

### ✅ decay_method="energy" is usable at this device — 3.1x cheaper (jobs 12276203, 12309040)

Whole-device, band centre, scored against the trusted dft/1e-11 reference
(job 12029541: T=[0.239975, 0.242448, 0.242448, 0.239975], total 0.9648).

| method / decay_tol | hrs | speedup | total T | max per-port dev |
|---|---|---|---|---|
| dft, 1e-11 (reference) | 41.03 | 1.0x | 0.9648 | — |
| energy, 1e-3 | 8.87 | 4.6x | 0.9611 | 1.27e-03 |
| **energy, 1e-5** | **13.33** | **3.1x** | **0.9649** | **1.87e-04** |

At 1e-3 the run stops before the fields settle — every port reads low and
total T is down 0.0038, the systematic signature of truncating a still-decaying
field (not random noise), which is why tightening fixed it rather than
reshuffling it.

At **1e-5**: total T agrees with the reference to **1e-4** (0.9649 vs 0.9648)
and max per-port deviation is **1.87e-04**.

> **Stated honestly: 1.87e-04 is ~1.9x the ≤1e-4 bar I set BEFORE running, so
> it does not meet my own pre-registered criterion.** It is nonetheless usable,
> and here is the sizing rather than a moved goalpost: 1.87e-04 is ~13x below
> the FDTD split asymmetry being measured (2.5e-03) and ~16x below the
> EME-vs-FDTD engine gap (3.07e-03). It is comfortably beneath every physical
> effect currently under study, but it is NOT good enough to resolve effects at
> the 1e-4 level, and should not be used for one.

**`decay_tol` means different things to the two methods** — a per-DFT-component
relative threshold for `"dft"`, a BULK-ENERGY decay factor for `"energy"`. So
1e-11 is a sensible dft threshold and an absurd energy target; these are not
"the same tolerance, different method".

**Caveat on what this establishes.** The reference is itself a dft/1e-11 run,
so this measures **agreement between two stopping criteria at this device**,
not an absolute accuracy bound. Enough to justify using `"energy"` where the
cost matters; not an independent validation of either.

**Recommended use: the tapered-body device's own whole-device FDTD reference**,
which does not exist yet and would otherwise cost ~41 hrs at dft. Nothing
changes for check 4 (job 12260631), which is already ~56% through under dft --
switching now would discard ~34 hrs to save ~28. Library defaults unchanged.

---

### ✅ CHECK 4 COMPLETE — full band (job 12260631, 2d18h, 43,225 steps)

**Self-check passes.** This run is the N_WL=11 counterpart of the band-centre
probe (12029541) -- same method, res, domain -- and reproduces it:
band centre T=[0.2400, 0.2425, 0.2425, 0.2400], total **0.9649** vs the probe's
[0.23998, 0.24245, 0.24245, 0.23998] / **0.9648**, agreeing to ~1e-4 despite
43,225 vs 40,420 steps and ~44 vs ~4 DFT components. `decay_method="dft"`
did not stall.

> ⚠ This run printed `Comparison -- EME (job 11840722 ...) total=0.9620`.
> **Disregard it.** 11840722 is `beef_mmi1x4_toc`, a step-9c DESIGN-SWEEP array
> from 2026-08-14 (the day the geometry changed), not a production EME run. The
> script was corrected 2026-08-26, but this job had already imported the old
> version. Authoritative EME reference is 11887395.

#### THE BAND-CENTRE COMPARISON WAS THE MOST FLATTERING ONE POSSIBLE

Every check-1 comparison to date was made at 1.55 µm. Across the band the two
engines have **opposite spectral slopes**, and band centre is almost exactly
where they cross:

| wl | FDTD total T | EME total T | EME−FDTD | FDTD shape | EME shape |
|---|---|---|---|---|---|
| 1.50 | 0.9846 | 0.9445 | **−0.0400** | OUTER | OUTER |
| 1.52 | 0.9821 | 0.9605 | −0.0216 | INNER | OUTER |
| 1.54 | 0.9723 | 0.9690 | −0.0033 | INNER | OUTER |
| **1.55** | **0.9649** | **0.9707** | **+0.0058** | INNER | OUTER |
| 1.56 | 0.9564 | 0.9705 | +0.0141 | INNER | OUTER |
| 1.58 | 0.9356 | 0.9652 | +0.0296 | INNER | OUTER |
| 1.60 | 0.9112 | 0.9538 | **+0.0426** | OUTER | OUTER |

- **FDTD decreases monotonically** across the band: 0.9846 → 0.9112, a **7.3
  point** fall. Physically sensible for an MMI — L_pi ∝ W²/λ, so a fixed
  mm_length drifts off the self-imaging condition as λ grows.
- **EME peaks at band centre** and is roughly symmetric (0.9445 → 0.9707 →
  0.9538). A qualitatively DIFFERENT spectral shape, not an offset.
- **The gap swings sign**, −0.0400 to +0.0426 (range 0.083), passing through
  ~0 near 1.55. Max per-port deviation across the band is **1.25e-02** (at
  1.50) — **4x larger** than the 3.09e-03 measured at band centre.
- **Split shape:** EME favours OUTER at **11/11** wavelengths; FDTD favours
  INNER at **9/11** (flipping to OUTER only at the extreme edges 1.50/1.60).
  So the sign disagreement is not a band-centre quirk either.

**Consequence for how check 1 must be reported.** Previous statements
("engines disagree by ~2.5 points / ~3e-03 per port") were measured at the one
wavelength where the disagreement is smallest, and understate it by ~4x.
Quoting band-centre agreement in the manuscript would be misleading. The
honest statement is: **the two engines disagree on the SPECTRAL SHAPE of this
device, not by a constant offset**, and the whole-device FDTD (no interface
truncation) is the reference.

Plot: `example_output/mmi_1x4/mmi_1x4_fdtd/check1_engine_gap_vs_wavelength.png`

#### Truncation is NOT the explanation — the TAPER staircase is the new suspect

Free test (all three spectra already on disk, no new compute). If mode
truncation drove the spectral-shape difference, agreement with FDTD should
degrade with truncation: FDTD (none) → EME (N_BODY=20) → hybrid (12). It does
the opposite.

| engine | truncation | total-T slope | mean \|dev vs FDTD\| |
|---|---|---|---|
| FDTD | none | **−0.749 /µm** | — |
| hybrid | interface_modes=**12** | **−0.338 /µm** | **0.0184** |
| EME | N_BODY=**20** | **+0.087 /µm** | 0.0235 |

The **more**-truncated hybrid is **closer** to FDTD on both measures. So
`N_BODY`/`interface_modes` does not explain the shape, and item 9d's
interface-modes sweep — ~3 days of prewarm per mode count — would **not** have
answered this question. Worth knowing before spending it.

**What actually separates them is the TAPER treatment**, and it lines up
exactly with the slope ordering:
- FDTD: tapers meshed exactly (slope −0.749)
- hybrid: tapers solved in **FDTD**, body in EME (slope −0.338)
- EME: tapers as an **8-cell staircase** (`C_TAPER=8`), body in EME (+0.087)

EME is the only path approximating the tapers, and the only one whose spectrum
is flat/symmetric rather than monotonically falling. Each staircase step jumps
0.1375 µm in width (1.1 → 2.2 µm over 5 µm in 8 cells), and the mode mismatch
at each step is wavelength-dependent — a plausible mechanism for a
spectral-shape error specifically.

**HYPOTHESIS, NOT ESTABLISHED.** The test is a `C_TAPER` convergence sweep
(8 → 16 → 32) on the whole-device EME: if EME's slope steepens toward FDTD's
−0.749 as the staircase refines, the taper discretisation is the cause. Cheap
if run at 3 wavelengths (1.50/1.55/1.60) to capture slope rather than the full
11-point band.


**What this does NOT establish.** Which engine is right at the band edges. FDTD
is the better-grounded reference (whole-device, no mode truncation, and its
monotonic trend has a physical mechanism), but nothing here independently
confirms its band-edge values. The EME path's fixed N_BODY=20 truncation is a
candidate cause for its flatter/symmetric shape — untested.

---

### ✅ CHECK 4 — whole-device FDTD reference, MEASURED at band centre (job 12029541)

**The check-4 reference exists and has since 2026-08-21.** Job 12029541
(`beef_inv_wholedevice_fdtd`, `investigate_mmi1x4_res30_fdtd_wholedevice_timing.py`)
COMPLETED, exit 0, **41.03 hrs / 40,420 timesteps**, and it ran at
**`decay_tol=1e-11` — the PRODUCTION decay tolerance**, res=30, the real
whole-device geometry. The only reduction from production is `N_WL=1`
(band centre only), which the script's own docstring notes does **not**
reduce decay cost. So these numbers are production-quality at λ=1.55 µm.

**This settles the check-1 engine disagreement.** All four runs, same
geometry, λ=1.55 µm:

| run | outer (0,3) | inner (1,2) | favours | spread | total T |
|---|---|---|---|---|---|
| EME, 11887395 | 0.2431 | 0.2423 | **outer** | 0.0008 | 0.9707 |
| Hybrid uncal, 12227682 | 0.2287 | 0.2460 | inner | 0.0173 | 0.9495 |
| Hybrid cal, 12227591 | 0.2280 | 0.2447 | inner | 0.0167 | 0.9454 |
| **FDTD whole-device, 12029541** | **0.23998** | **0.24245** | **inner** | **0.0025** | **0.9648** |

Three findings, none of them what the earlier framing assumed:

1. **The HYBRID has the split-shape sign RIGHT; EME has it WRONG.** FDTD
   favours the INNER pair, as both hybrid runs do. Step 6d's long-standing
   framing — that the hybrid's outer-favouring was the anomaly needing
   explanation — had it backwards for the current geometry: EME is the
   outlier here.
2. **But the hybrid badly overstates the asymmetry**: spread 0.0167 vs the
   reference's 0.0025, a **6.7× inflation**. Right sign, wrong magnitude.
3. **Total T**: reference 0.9648. EME overestimates by +0.6%; the hybrid
   **underestimates by −2.0%** — the direction and rough size mode
   truncation predicts (N_BODY=20 vs interface_modes=12, item 9d), and the
   best evidence yet that 12 modes is costing real accuracy.

> ### ⚠ The check-1 comparison was CONFOUNDED — retested 2026-08-25, conclusion survives
>
> The comparison above was first reported as "the whole-device FDTD arbitrated
> the engine disagreement." That claim was **not supported by those two runs**,
> because they differ in three ways at once, not one:
>
> | | EME (11887395) | FDTD ref (12029541) |
> |---|---|---|
> | resolution | **40** | **30** |
> | domain (y x z) | **20 x 5** | **16 x 4** |
> | body modes | 20 | none (whole-device) |
>
> The device geometry is identical; the SIMULATION parameters are not. The
> domain gap is not an oversight either: 16x4 was adopted 2026-08-19 on a
> convergence check whose own conclusion was "modes 19/22 never converge but
> are IRRELEVANT because production keeps N_MODES=12" -- true for the hybrid's
> 12-mode interfaces, but `mmi_1x4_eme.py` carries **N_BODY=20**, so it is the
> one caller that check does not cover.
>
> **Tested rather than argued** (job 12323220,
> `mmi_1x4_eme_split_shape_confound.py`, band centre only):
>
> | case | res | domain | total T | favours | dev vs FDTD | dev vs EME pub |
> |---|---|---|---|---|---|---|
> | A | 40 | 20x5 | 0.9707 | **OUTER** | 3.09e-03 | **3.46e-05** |
> | B | 30 | 16x4 | 0.9703 | **OUTER** | 3.79e-03 | 9.10e-04 |
>
> Case A reproduces the published EME result to 3.46e-05, so the harness is
> faithful. Case B then matches the FDTD reference's OWN resolution and domain
> -- and EME **still favours OUTER**, with the outer/inner gap WIDENING
> (7.93e-04 -> 2.38e-03) and agreement with FDTD getting slightly worse.
>
> **So the split-shape disagreement is NOT a resolution/domain artefact.** It
> survives matching parameters, and the check-1 conclusion stands -- now on a
> controlled comparison. The correction was still warranted: the conclusion
> happened to survive, but it was not supported by the evidence originally
> cited for it. Cases C/D/E (domain isolated; res 50/60) are running to test
> whether EME converges toward the FDTD answer at all, which check 3's finding
> (EME's differential does not converge like meep's) predicts it may not.

**Cost question from `mmi_1x4_fdtd.py`'s docstring is also resolved.** It
projected 52 hrs if timesteps stayed near the input probe's 23,388, or
~390 hrs if they scaled with the ~7.5× longer propagation path, and refused
to commit blind. Measured: **40,420 steps — only 1.73× the input section**,
not 7.5×. The optimistic branch was right; whole-device FDTD on this device
is firmly tractable.

**What remains for check 4:** the full 11-wavelength band. Because meep runs
ONE time-stepping simulation with DFT monitors at all frequencies, adding the
other 10 wavelengths should cost roughly the same ~41 hrs, not 11×.
`mmi_1x4_fdtd.py` is already revived and pinned to the current geometry
(mtime 2026-08-19); it has never been run.

> ⚠ **Doc-hygiene:** the summary table at the top of this file still reads
> "`mmi_1x4_fdtd.py` dropped: `mm_length=206` makes whole-device single-domain
> FDTD impractical". That justification expired with step 9c's geometry change
> (216 µm → 49.90 µm, 4.3×), the script was revived 2026-08-19, and the
> reference has now been measured. Fixed in the table.

**9b. Symmetric mode reuse on the OUTPUT side — proposed, not implemented,
no code written.** `mmi_2x2_hybrid.py`/`directional_coupler_hybrid.py` each
cut their large-FDTD count in half by noticing their two *input* arms are
mirror images of each other (`SYMMETRIC_MODE_REUSE`: solve one arm, derive
the other via a measured per-mode sign vector, `MODE_SIGNS`). This device has
one input port (nothing to mirror there) but its **four output branches**
carry the same symmetry: `output_ys` is centered on the same y=0 axis the
single input sits on (`MMI._port_ys()`),
[⚠ this paragraph originally read `output_ys = [-3.6, -1.2, +1.2, +3.6]` µm,
which is the **superseded** `WAVEGUIDE_OFFSET=2.4` geometry; at the current
`WAVEGUIDE_OFFSET=3.1` it is `[-4.65, -1.55, +1.55, +4.65]` µm, verified
directly from `MMI.output_ys` 2026-08-22. The mirror-symmetry argument below
is unaffected — it depends only on the ports being centred on y=0, which
holds for both — but the numbers were wrong and are the kind of thing that
gets quoted.]
so output_idx=0 is the y-mirror of output_idx=3, and output_idx=1 of
output_idx=2 — structurally the same trick, just on the other side of the
device.

**Sizing.** Step 6d's real run solves 1 input + 4 outputs = **5** large FDTD
sections per full sweep. Simulating only output_idx∈{0,1} directly and
deriving output_idx∈{3,2} via mirroring cuts that to 1 input + 2 output =
**3** — a 40% reduction in the dominant cost (each section is the
~12–24 hr class measured in step 6d, not the cheap calibration solves of
9a), i.e. plausibly saving one full ~1.7-day array-task run out of the
~1.7 days a 4-branch prewarm currently costs per branch pair.

**What landing this actually needs — none of it done yet:**

1. **Measure the sign vector for real**, the same way `mmi_2x2_hybrid.py`'s
   `MODE_SIGNS` (line 137) and `directional_coupler_hybrid.py`'s (line 263)
   were measured: run both members of one mirror pair (e.g. output_idx=0 AND
   output_idx=3) for real, `cache_only`-style comparison of their per-mode
   coefficients, confirm each body mode's contribution either matches or
   flips sign between the two — do not assume the sign vector is `[1,-1,1,-1,...]`
   by analogy; both existing examples measured theirs independently and got
   different patterns (`[1,-1,1,-1,1,-1,1,-1,-1,-1,-1,1,-1,1]` for the 2×2,
   `[1,-1,-1,1,1]` for the coupler).
2. **Rewrite the per-output loop as a hand-rolled cascade**, not a
   `hybrid.solve_length_sweep()` call per branch. `HybridSolver`'s one-call
   contract solves a section or derives it from `mirror_pairs` *within one
   chain*; deriving output_idx=3's transfer matrix from output_idx=2's
   already-computed result, across two otherwise-identical
   `solve_length_sweep()` calls, isn't expressible through that API —
   exactly the reason `mmi_2x2_hybrid.py`'s own module docstring gives for
   hand-rolling its cross-path einsum instead of using `HybridSolver.solve()`
   for it ("HybridSolver's one-call contract can't express [this]").
3. **Re-verify the geometry is exactly y-mirror-symmetric** before trusting
   step 1's derivation — `use_neighbor_aperture=True` and the per-branch FDTD
   monitors should already respect this by construction (A1/6a), but no test
   has explicitly asserted `width_at`/`center_at` symmetry across y=0 for
   this device the way `HybridSolver`'s own `mirror_check_tol` guardrail does
   for `mirror_pairs`-based chains — worth the same check here, cheaply,
   before spending cluster time on a 4-branch reference run.

**Interacts with 9a**: if `calibrate_body_neff` changes the numbers, the
sign vector must be measured AFTER 9a lands, on the calibrated cascade — not
before, or 9b would freeze a sign pattern against soon-to-be-stale physics.
Land 9a first.

**9c. Joint taper-offset x mm_length design sweep — user-requested
2026-08-13, IN FLIGHT.** Neither 5's `mm_length` sweep (fixed offset=2.4) nor
the earlier `mmi_1x4_sin_taper_offset_sweep.py` (fixed mm_length=44.0, offset
swept) ever asked the joint question: for EACH `waveguide_offset`, what is
the best `mm_length`? New script
`examples/verification/python/mmi_1x4_sin_taper_offset_length_sweep.py`
(EME, `mode_solver="tidy3d"`, all-guided per-branch `solve_length_sweep()`,
same free-length-sweep trick as step 5) answers it directly: one fresh
`EMESolver.solve_length_sweep()` per `(offset, branch)` pair (offset is NOT
free like length — it's a genuinely different cross-section at every one of
the 19 cells), sweeping the full `mm_length` grid at each.

**Offset grid, decided with the user after two feasibility rounds:**
`[2.2, 2.5, 2.8, 3.1, 3.2167] µm`. `mm_width` held FIXED at today's
production `12.0 µm` (user's explicit choice over widening it — keeps
`N_BODY=20`'s mode-count basis, sized off 23 guided modes at `mm_width=12.0`
in job 11660226, valid at every point without re-deriving it for a wider
body). This pins the offset ceiling at this device's own geometric
containment limit, `offset < (mm_width − taper_width)/3 = 3.2667 µm`
(`mmi_1x4_sin_taper_offset_sweep.py`'s own derivation) — the swept 3.2167 is
that ceiling minus the same 0.05 µm safety margin that script's own
`_FEASIBLE_MAX` already uses (the literal 3.2667 the user first named rounds
UP across the true repeating-decimal boundary and fails the containment
assert by ~5e-5 µm). The low end, 2.2 µm, sits exactly ON the
adjacent-taper-overlap floor (`offset > taper_width = 2.2 µm`, else
neighboring tapers overlap at the body face) — the user's own first choice
of 2.0 µm was rejected here for exactly that reason (0.2 µm of taper
overlap), and 2.2 is the boundary case: 0.0 µm edge-to-edge clearance,
touching but not overlapping. Both endpoints verified by direct arithmetic
before submission (not just asserted in-script): all 5 offsets checked
against both the containment ceiling and the adjacent floor, all pass.

**Config**: `N_WL=5, N_BODY=20, N_TAPER=10, C_TAPER=8`, `res=40` — pinned to
step 5's own real-run config (job array 11667114), not production's
`N_WL=11`, since this is an exploratory design sweep. `mm_length` grid
identical to step 5's (`linspace(20, 300, 281)`), so results compare
directly against that sweep's own numbers.

**Cost**: 5 offsets × 4 branches = 20 independent tasks, each the same cost
class as step 5's own real run (measured there: 5:49–6:29 hrs/branch).

**Status: ✅ DONE 2026-08-14 — all 20 tasks completed, MEASURED, a real
design finding.** Smoke test (job 11817341, reduced scope) confirmed the
script first, per this session's two-step discipline; real array (job
**11817364**, `--array=0-19`) then ran all 20 `(offset, branch)` pairs to
completion, 05:46–06:33 each, zero failures (`sacct`: 20/20 `COMPLETED`,
exit `0:0`; no `Traceback`/`AssertionError`/`OOM` in any of the 20 log
files). `mmi_1x4_sin_taper_offset_length_sweep_merge.py` combined them.

**⚠ Output directory correction — same off-by-one this session's own
`mmi_1x4_sin_body_length_sweep.py` already has (that script's own module
history, and the `mmi_1x4_hybrid_body_length_sweep.py` fix note earlier in
this doc, for the identical `save_dir` version of this exact mistake).**
Copying that script's own `os.path.join(_HERE, "..", "..", "..", …)`
(3 levels up from `examples/verification/python/`) lands at
**`BEEF/example_output/verification/mmi_1x4_sin_taper_offset_length_sweep/`**,
NOT `duty/example_output/...` as this doc's own 9c write-up first said
(4 levels up, matching the production scripts' own depth from
`examples/mmi_1x4/python/`). Not fixed — flagged only, same call this repo
already made for the sibling script: moving it now would just relocate a
complete, verified result set, and every file this session produced is
self-consistently under the `BEEF/`-nested path (the merge script's own
`output_dir` uses the identical 3-dot convention, so it found all 20 files
correctly). **Real path**: `BEEF/example_output/verification/
mmi_1x4_sin_taper_offset_length_sweep/` (relative to repo root).
`save_dir` (the EME mode-solve cache) landed the same way, at
`BEEF/example_sim_cache/` — a fresh, disconnected cache pool, not the real
device cache at `duty/example_sim_cache/`. Harmless here for the same
reason it was harmless for the sibling script: this sweep's cross-sections
(5 new offsets) were never in the real cache to begin with, so nothing was
missed by not sharing it.

**Per-offset best (min spread, naive global argmin — what the merge script
reports by default):**

| offset (µm) | best mm_length (µm) | spread | total T |
|---|---|---|---|
| 2.2000 | 113.00 | 0.0001 | 0.5120 |
| 2.5000 | 205.00 | 0.0043 | 0.6112 |
| 2.8000 | 286.00 | 0.0005 | 0.6302 |
| 3.1000 | 72.00  | 0.0001 | 0.4426 |
| 3.2167 | 197.00 | 0.0002 | 0.7444 |

**But the real finding is NOT the naive argmin — it's a short-device family
near L≈40 µm that the per-offset table above misses entirely** (argmin picks
the single tightest-balance point regardless of throughput, and at several
offsets the tightest-balance point sits at low throughput — e.g. 3.1000's
own "best," L=72, T=0.44). Scanning the full `[offset, mm_length]` grid for
points combining low spread AND high throughput surfaces a **new standout
candidate, dominating every previously known point on this device**:

> **offset = 3.1000 µm, mm_length = 40.0 µm → total T = 0.9598 (IL = −0.18
> dB), spread = 0.00185 (0.77% relative imbalance).**

Neighbors confirm this is a real, smooth peak, not a single-pixel grid
artifact: L=39 → T=0.9496/spread 0.0041 (1.72%), L=41 → T=0.9459/spread
0.0024 (1.01%). A second, slightly weaker peak sits at offset=3.2167,
L=38–42 (T=0.90–0.94, 0.68–3.4% imbalance).

**Comparison against every other candidate on record for this device:**

| design point | mm_length | total T | rel. imbalance |
|---|---|---|---|
| current production | 206.0 µm | 0.6014 | 4.26% (step 6d/9) |
| earlier 44 µm + offset≈3.0–3.2 exploration (coarse, N_WL=1) | 44 µm | 0.75–0.76 | ~15% |
| **this sweep: offset=3.1, L=40** | **40.0 µm** | **0.9598** | **0.77%** |

That is **~60% more throughput than production, ~5.5x shorter than
production, and better balanced than either prior candidate** — not a
marginal improvement, a different regime. Physically plausible, not an
obvious artifact: T=0.96 sits close to the lossless ideal (1.0) with all
four outputs landing within a few percent of each other
(`[0.2405, 0.2394, 0.2394, 0.2405]`-scale per-output values), the signature
of a genuine near-perfect 4-way self-imaging condition rather than basis
truncation noise.

**Caveats, stated plainly — NOT yet a validated design change:**
- Design-tool fidelity only: `N_WL=5` (not production's 11-point band),
  `N_BODY=20/N_TAPER=10/C_TAPER=8` (same preliminary truncation step 6
  already flagged as unvalidated for the production device).
- 1 µm grid spacing — the true optimum could sit slightly off L=40; a finer
  local sweep (e.g. 35–45 µm at 0.1 µm steps, still cheap under
  `solve_length_sweep()`'s free-length trick) would sharpen this before
  treating 40.0 as pinned.
- Single reference `mm_width=12.0` throughout (this sweep's own user-decided
  scope) — untested whether a wider body pushes this peak even higher, the
  same open question step 9c's own module docstring already flags for the
  offset ceiling.
- No hybrid/FDTD cross-check yet (check 1, still open) and no resolution
  convergence check (check 3, still open) — both apply to this candidate
  exactly as much as to the production point.

**Next step, not yet done**: a real-config confirmation run (`N_WL=11`,
finer local `mm_length` grid around 40 µm, offset fixed at 3.1) at this
candidate point, the same "promote a design-tool finding to a validated
result" step every prior candidate on this device has needed.

#### Local refinement — ✅ DONE 2026-08-14, peak confirmed and sharpened

User asked for a follow-up local refinement (offset 2.90–3.20 µm in 0.05 µm
steps × mm_length 30–50 µm in 0.1 µm steps, `N_WL=1` for a fast localization
pass — full config otherwise unchanged). New scripts
`mmi_1x4_sin_taper_offset_length_refine.py` /
`..._refine_merge.py`, same two-step (smoke, job 11839141, then real 28-task
array, job **11839158**) discipline. **All 28 tasks completed cleanly**
(sacct: 28/28 `COMPLETED`, exit `0:0`, 55–58 min/task at `N_WL=1` — the ~5x
speedup vs. the coarse sweep's `N_WL=5` cost was a correct estimate).
Output: `BEEF/example_output/verification/mmi_1x4_sin_taper_offset_length_refine/`.

**Refined peak: offset=3.1000 µm (unchanged from the coarse grid),
mm_length≈39.9–40.1 µm → total T=0.9707, spread=0.00079 (0.33% relative
imbalance)** — better than the coarse grid's own 40.0 µm estimate
(T=0.9598, 0.77%) simply because the finer grid sits closer to the true
continuous optimum. The offset did NOT move between the 0.3 µm-coarse and
0.05 µm-fine grids — real evidence this is a genuine peak, not a grid
artifact.

**The peak is broad, not a knife-edge — a favorable property for a real
device.** Every `mm_length` in [39.5, 40.3] µm gives T=0.968–0.971 with
<0.7% imbalance (top-10 table, all at offset=3.1 except one point at
offset=3.05/L=39.9 also landing at T=0.968) — a ~0.8 µm-wide plateau, not a
single sharp point, meaning the design tolerates realistic length variation
without falling off this optimum.

**Still not yet a validated design change** — same caveats as before
(design-tool `N_BODY/N_TAPER/C_TAPER` fidelity, `N_WL=1` now instead of 5
for this pass specifically, single `mm_width=12.0`, no check 1/3 yet). Next
step unchanged: a real `N_WL=11` confirmation pinned at offset=3.1,
mm_length=40.0 (or 39.9–40.1, functionally equivalent given the plateau).

**Output**: `BEEF/example_output/verification/mmi_1x4_sin_taper_offset_length_refine/refine_offset{0-6}_branch{0-3}.npz`
(28 raw files) and `taper_offset_length_refine.npz`/`.png` (merged result +
plot).

#### Full-fidelity confirmation — ✅ DONE 2026-08-14, holds up at N_WL=11

New scripts `mmi_1x4_sin_taper_offset_confirm.py` / `..._confirm_merge.py`:
offset FIXED at the confirmed 3.1 µm candidate (not swept — this is a
confirmation, not another search), `mm_length` grid narrowed to
[38.0, 42.0] µm (0.1 µm steps, bracketing the refined plateau), **N_WL=11**
— the full production C-band, matching `mmi_1x4_eme.py`/
`mmi_1x4_hybrid.py`'s own real config for the first time on this candidate.
Two-step discipline again: smoke (job 11840705) then real 4-task array (job
**11840722**, one task per output branch — no offset axis this time).

**All 4 tasks completed cleanly** (sacct: 4/4 `COMPLETED`, exit `0:0`,
4:17–4:41 hrs/task — notably faster than the ~10–14 hr estimate that
assumed a steeper-than-actual N_WL cost scaling; real cost did not scale as
steeply with N_WL as the coarse sweep's own N_WL=1→5 step suggested).

**CONFIRMED: offset=3.1000 µm, mm_length=39.90 µm → total T=0.9620
(IL=−0.17 dB), spread=0.00180 (0.75% relative imbalance)**, per-output
band-mean T = `[0.2414, 0.2396, 0.2396, 0.2414]` (outputs at
y = ∓4.65, ∓1.55 µm at this offset). This lands squarely inside the range
both design-tool passes already predicted (T=0.960–0.971, imbalance
0.33–0.77%) — the reduced-`N_WL` shortcuts were NOT hiding a
band-averaging artifact; the full 11-point C-band average is consistent
with both the 5-point and 1-point estimates to within a few tenths of a
percent.

**Final comparison, this device's full candidate history:**

| design point | mm_length | offset | total T | rel. imbalance |
|---|---|---|---|---|
| current production | 206.0 µm | 2.4 µm | 0.6014 | 4.26% |
| 44 µm + offset exploration (coarse, N_WL=1) | 44.0 µm | ~3.0–3.2 µm | ~0.755 | ~15% |
| **this candidate — CONFIRMED (N_WL=11)** | **39.90 µm** | **3.1 µm** | **0.9620** | **0.75%** |

**~60% more throughput, ~5.2× shorter multimode section, and better
balance than production — measured at full production fidelity, not just
design-tool shortcuts.** This is now a real, evidenced candidate to replace
the production `MM_LENGTH=206.0`/`WAVEGUIDE_OFFSET=2.4` pair, pending the
items below.

**What this confirmation does NOT cover — still genuinely open before
calling this a shipped design change:**
- No cross-engine check (check 1) — this candidate has only ever been run
  on EME (`tidy3d`); no hybrid/FDTD cross-check exists for it the way
  production's own split-shape discrepancy (EME vs. hybrid, step 6d) is
  still being chased.
- No resolution convergence check (check 3) — `res=40`/`N_BODY=20`/
  `N_TAPER=10`/`C_TAPER=8` throughout, same preliminary truncation flagged
  as unvalidated for production.
- `mm_width=12.0` held fixed the whole investigation — never tested whether
  a wider body pushes this peak even higher (the offset ceiling this
  candidate sits near, 3.2167 µm, is a `mm_width=12.0`-specific limit).
  **User's explicit call, 2026-08-14: not blocking on this — ship with the
  current `mm_width=12.0`.**
- `mmi_1x4_hybrid.py` not yet updated — still carries the old
  `MM_LENGTH=206.0`/`WAVEGUIDE_OFFSET=2.4` pair. This is also where the
  cross-engine (check 1) gap above would actually get closed, since that
  script is the only one of the two production scripts that touches FDTD.

**`mmi_1x4_eme.py` updated 2026-08-14** — `MM_LENGTH`/`WAVEGUIDE_OFFSET`
changed to `39.90`/`3.1` in place, full provenance in the module docstring
and the constants' own comments. One mechanical bug caught and fixed while
editing: the old `mm_lengths` grid (`linspace(20, 300, 281)`, integer
spacing) doesn't contain `39.90`, and switching to a finer
`linspace(20, 300, 2801)` (0.1 µm spacing) lands `39.90` only to within
~7e-15 of exact due to `linspace`'s own float accumulation — the script's
old `==` grid-alignment assertions (both of them) would have failed on a
real floating-point non-issue. Fixed by switching to
`np.isclose(..., rtol=0, atol=1e-9)`, the same convention
`directional_coupler_hybrid.py` already uses for this identical
pinned-value-on-grid check — tight enough to still catch a real mis-snap
(the coupler's own 0.0208 µm discrepancy) while tolerating genuine
floating-point noise. Confirmation run's own cache (`mmi_1x4_sin_taper_
offset_confirm.py`, job 11840722) lived under `BEEF/example_sim_cache/`
(off-by-one path, same as flagged for the sweep/refine scripts) while
`mmi_1x4_eme.py` reads from `duty/example_sim_cache/` (its own correct
depth) — so this script's own real run needed a fresh prewarm, not a reuse
of that cache.

**✅ DONE 2026-08-15 — prewarmed and run for real, production artifacts
produced.** Prewarm (`_prewarm_mmi_1x4_eme_sweep.py`) hit the length-sweep
cost bug documented in the next section below before landing cleanly; once
fixed (mm_lengths narrowed to the same `[35,45]`/101-point window as the
production script), the real 4-task prewarm array (job 11871657) completed
in 6:37–6:50 hrs/branch, and `mmi_1x4_eme.py` itself then ran against the
warm cache in **38 seconds** (job 11887395, matching the historical ~38s
precedent for the old geometry exactly):

| output | y (µm) | T | IL (dB) |
|---|---|---|---|
| 0 | −4.65 | 0.2431 | −6.14 |
| 1 | −1.55 | 0.2423 | −6.16 |
| 2 | +1.55 | 0.2423 | −6.16 |
| 3 | +4.65 | 0.2431 | −6.14 |
| **total** | | **0.9707** | **−0.13** |

Consistent with the earlier full-band confirmation (T=0.9620, job
11840722) — small difference from the finer local grid, not a
discrepancy. The script's own local length-sweep optimum within `[35,45]`
(40.3 µm, spread 0.0016) is essentially indistinguishable from the chosen
39.9 µm (spread 0.0018), reconfirming the broad plateau found by the
refine sweep. All production artifacts saved: `geometry.png`,
`performance.png`, `smatrix.pkl`, `mm_length_sweep.png`/`.npz` in
`example_output/mmi_1x4/mmi_1x4_eme/`.

**Output**: `BEEF/example_output/verification/mmi_1x4_sin_taper_offset_confirm/confirm_branch{0-3}.npz`
(4 raw files) and `taper_offset_confirm.npz`/`.png` (merged result + plot).

**Output** (paths corrected per the note above):
`BEEF/example_output/verification/mmi_1x4_sin_taper_offset_length_sweep/sweep_offset{0-4}_branch{0-3}.npz`
(20 raw files) and `taper_offset_length_sweep.npz`/`.png` (merged result +
plot).

#### ⚠ Real bug found 2026-08-14: `EMESolver.solve_length_sweep()`'s "free length sweep" is NOT free at large point counts

While prewarming `mmi_1x4_eme.py`'s own cache for the new geometry
(`MM_LENGTH=39.90`/`WAVEGUIDE_OFFSET=3.1`), the real 4-task prewarm array
(job **11848462**, full production config `N_WL=11`, `mm_lengths =
linspace(20, 300, 2801)`) **timed out at 10 hours with zero progress
markers** (`mode_solver="tidy3d"` gives none, same diagnostic blindness
noted for this device's own earlier length-sweep-plot attempt below).

**Root cause, found by reading `EMESolver.solve_length_sweep()`'s own
source** (`beef/solvers/eme.py`), not guessed:

```python
for wl_idx, wl in enumerate(wavelengths):
    ...
    modes = [self.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]  # shared across lengths
    for k, L in enumerate(lengths):
        ...
        S, _ = mw.compute_s_matrix(modes, cells=cells, cell_lengths=cell_lengths)  # NOT shared
```

The eigenmode **solve** genuinely is computed once per wavelength and reused
across every swept length, exactly as the docstring promises. But
`mw.compute_s_matrix()` — the interface-stitching + propagation **cascade**
across this device's full 19-cell chain (`wg_in | 8-cell taper_in | body |
8-cell taper_out | wg_out`) — is called **once per length inside that
loop**, and its own per-call cost is not negligible at this chain size.
2801 lengths × 11 wavelengths = **30,811** such calls, vs. the CONFIRMED
design point's own successful run at 41 lengths × 11 wavelengths = **451**
calls (`mmi_1x4_sin_taper_offset_confirm.py`, job 11840722, 4:17–4:41 hrs)
— a ~68× jump in call count for identical per-wavelength mode-solve cost.

**This also retroactively explains an earlier unresolved mystery**: this
same script's own N_WL=11 length-sweep-plot attempt (job 11783204, only
**281** lengths) timed out at 8h with the identical zero-progress symptom,
attributed at the time only to "N_WL scaling" (see below). The real driver
was ALSO this per-length cascade cost compounding on top of the mode-solve
cost — 281×11 = 3,091 calls, already enough to matter.

**`mmi_1x4_hybrid.py` is NOT affected** — its body-length sweep goes
through `EMESolver.solve_section_length_sweep()` (the **single-cell**
variant: just the uniform body cross-section, no taper staircase in this
call at all), where each `compute_s_matrix()` call is trivial. MEASURED:
25.3s for all 4 branches × 281 lengths (job 11770532). The cost blowup is
specific to the **multi-cell, full-device** `solve_length_sweep()` method
`mmi_1x4_eme.py` alone uses.

**Fix**: `mmi_1x4_eme.py`'s (and its prewarm companion's) `mm_lengths`
narrowed from the full `[20, 300]` (2801 pts) to `[35, 45]` (101 pts,
`np.arange(35.0, 45.001, 0.1)`), bracketing the confirmed `MM_LENGTH=39.90`
with real margin (the refine sweep's own plateau was 39.5–40.3) at a call
count (101×11=1,111) close to the confirmed-safe scale above. The broader
design-space context this grid used to show is not lost — it already lives
in `mmi_1x4_sin_taper_offset_length_sweep.py`'s/`..._refine.py`'s own plots
(this section), just not duplicated here at prohibitive cost.

**Status**: fix applied to both scripts; re-smoke-tested; real prewarm
resubmitted (job 11871620 smoke → real array to follow) under `gts-sr108`
(also switched off `gts-sr108-ece`, which is separately at 0.998 fair-share
usage — see [[slurm-account-gts-sr108-ece]]).

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
| ~~5 — body design (`mm_length` sweep)~~ | job array 11667114 | ✅ **done 2026-08-04**; re-chosen at step 9c (2026-08-14) → `mm_length=39.90` |
| ~~6a — geometry plot~~ | none | ✅ **done 2026-08-04** |
| ~~6b–d — three scripts~~ | 1 array job | ✅ **done 2026-08-04** — two scripts, not three (see 6c); `mmi_1x4_fdtd.py` revived 2026-08-19 once step 9c made it affordable |
| ~~7 — Makefile/batch/exports~~ | none | ✅ **done 2026-08-05** |
| ~~8 — four checks + material gate~~ | jobs 12275175, 12029541, 12260631, 12323220/12331012 | ✅ **all four MEASURED**; material gate substantively passes, harness call open — see the RESOLVED header |

~~The critical path is **1 → 2 → 3 → 4 → 5 → 6**.~~ **Every step above is done.** This table described the un-started plan long after the work landed; corrected 2026-08-31 on resolution. What remains is not a step but the three residuals in the RESOLVED header at the top of this file. Steps 1 and 3 need no cluster
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
