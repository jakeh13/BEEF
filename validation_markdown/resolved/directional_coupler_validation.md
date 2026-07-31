# Directional coupler validation

Status, open work, and outstanding verification. **Everything already closed — items 1–7,
9–12 and 14, with their measurements, jobs and reasoning — lives in
[`archive/dc_validation_20260727.md`](../archive/dc_validation_20260727.md).** (Item 14 closed
2026-07-28, after that file's own date; it was updated in place rather than left stale.)

> **This file was `dc_validation.md` until 2026-07-27.** Roughly 76 references across the
> library, tests and examples cite it by that name, usually as *"dc_validation.md item N"*
> — and almost always a *completed* item. Every numbered item and its evidence is
> preserved in the archive above, so those citations still resolve in one hop. Item
> numbers are stable IDs and are never reused.

## Status — 2026-07-28

**The hybrid chain converges, and the coupler's hybrid results are trustworthy.**

Item 6, the last obstacle and the reason the original document existed, is closed.
(Naming collision: `fdtd_vs_EME_mode_validation.md` has its own, unrelated, still-open
item 6 — `interface_modes` convergence. `verify_item_6_engine_agreement.py` and
`run_item_6_res_sweep.batch` are THIS item 6, closed, not that one.) The
two engines disagreed by ~10% on Δβ through an identical uniform cross-section, which made
the chain's answer depend on where the FDTD/EME cut was made — the one thing a hybrid
decomposition must not do.

The fix was **not** to reconcile the engines. Both mode solvers are MPB and already agreed
to 0.57%; the gap was between MPB's mode index and the phase meep's time-stepper
accumulates (~2.3% per mode, same sign for both supermodes, levered into ~10% on Δβ
because `β/Δβ ≈ 99`). So `HybridSolver(calibrate_body_neff=True)` measures β from the
FDTD's own accumulated phase and propagates the EME body on that.

Job 11517245, against the uncalibrated job 11500412:

| | drift with cut position | L₅₀ spread over 9 stubs | mean L₅₀ | vs 15.13 µm reference |
|---|---|---|---|---|
| uncalibrated | +0.108543 µm/µm | 0.220 µm | 14.429 | **−5.37%** (at stub 0) |
| **calibrated** | **−0.000651 µm/µm** | **0.005 µm** | **15.1468** | **+0.111%** |

167× less drift, below the run's own 0.0015 µm rms scatter — the cut point no longer
measurably matters.

**The flagship tutorial now uses it too** (item 14, job 11526547): its auto-picked 50:50
length moved 14.5 → **15.15 µm**, landing 0.14% from the full-device reference — 0.07 of
its own sweep-grid step, i.e. at the measurement's resolution floor.

**Read that as self-consistency, not proven physical accuracy.** The reference is meep's
own discretization, and MPB's Δβ disagrees with it by ~5% with *both* engines looking
converged. Which is closer to Maxwell is open as **item 13**, and that caveat belongs
wherever the 0.11% is quoted.

### Where things stand per check

The four checks are defined in [`README.md`](../README.md).

| check | coupler |
|---|---|
| 1 engine agreement | ✅ closed — β-calibrated, L₅₀ 0.11% off the reference |
| 2 mode identity | ✅ identity permutation, worst guided disagreement 0.0023 |
| 3 resolution convergence | ✅ engines split — FDTD 40, EME 80 |
| 4 full-device FDTD reference | ✅ done, with provenance in the file |

One caveat carried forward from check 2: the coupler's interface basis asks for `N_MODES =
6` and only **5 are guided** — band 5 sits at n_eff 1.4367, below the 1.4443 cladding.
Dropping to 5 is free and has not been done.

## To-do

**Open: 8 (step 2 only), 13, 15**, plus one regeneration run. Item 14 closed 2026-07-28
(job 11526547 — the tutorial's 50:50 length moved 14.5 → 15.15 µm, 0.14% from the
full-device reference) and its full record has moved to the archive. Item 8's step 1 is
done: both engines now share one index model and agree to 4e-16, and the C-band kink is
gone. Its step 2 (a dispersive FDTD medium) was implemented, measured 0.0496 *worse* in
n_eff, and reverted — see the item for why and what would make it work; the checked plan,
with a strip-test gate, is [`material_model_validation.md`](material_model_validation.md).

Nothing from either is in the published example outputs or the goldens yet. All four
goldens are `@pytest.mark.slow`, were never run during this work, and are **stale** —
item 8 shifted n_Si by ~0.002 for every device. One cold full-device pass covers the
goldens, the `COUPLING_LENGTH` change and the material change together; the emptied cache
means it is one pass rather than three.

13. **Which engine's Δβ is physically right?** *(opened 2026-07-27.)*
    MPB plateaus at Δβ = 0.0894 by res 80–100; Δβ_FDTD sits flat at ~0.085 across res
    20–50. **Two independently converged answers ~5% apart** means at least one carries a
    systematic that a resolution sweep cannot see, and nothing measured so far says which.

    It matters because L₅₀ ∝ 1/Δβ almost exactly here — the ratios agree to **0.057
    percentage points**:

    | ratio | value |
    |---|---|
    | L₅₀: full-device FDTD 15.13 ÷ uncalibrated hybrid 14.429 | 1.04858 (+4.86%) |
    | Δβ: converged MPB 0.0894 ÷ FDTD 0.085212 | 1.04915 (+4.91%) |

    So the entire L₅₀ disagreement between the two methods is the Δβ disagreement and
    nothing else — no bend phase, no seam mismatch, no basis truncation left over. "Which
    L₅₀ is right" reduces completely to this question. If MPB is right, the *uncalibrated*
    hybrid was closer to physical truth all along and the calibration is faithfully
    reproducing meep's error.

    **The test that would settle it** is scoped and not run — see Verification.

8. **`notes.md` Solver #9's material half.** EME's Δn_eff has a +0.00022 kink at exactly
   1.57 µm where meow's Si table switches 3.4784 → 3.4710, while the FDTD half stays
   frozen at 3.4784 for every wavelength. Real, ~0.8% of Δn_eff — worth fixing for
   correctness, not a driver of anything. (Solver #8, wavelength interpolation, is not
   involved: the scripts do a full mode solve at all 11 wavelengths.)

   **Step 0 done, 2026-07-28 — the resolved indices are now in the cache fingerprint.**
   This had to come first: the cache recorded a `MaterialStack` as `{wg_height, sz}`,
   i.e. its *geometry* only, so nothing optical reached the key. Any fix below would have
   kept matching every existing entry and been served the old physics — visible only on a
   cold cache. Same defect class as item 3, surviving in this one corner.

   `resolved_material_config()` (`beef/materials/stack.py`) records the indices a run
   actually simulates with, splatted into `call_kwargs` at all seven config sites across
   both solvers. It follows the existing `_resolved_spectral_config()` /
   `resolved_layout` convention for resolved-vs-ingredient values. Proven both ways: a
   stored config still hashes to its own folder name, adding the field changes the hash
   (so old entries miss), and swapping snap→interpolate changes it *again* (so step 1
   will be seen).

   **Two things this turned up.**

   *It cannot be a `MaterialStack.cache_key()`,* which was the obvious design and is
   wrong. Resolving an index imports meow, and **`import meow` imports meep** — measured;
   meep lands in `sys.modules` on a bare `import meow`. `beef/cache/store.py::_is_master()`
   dispatches on whether meep is loaded while promising to import nothing itself, so a
   `cache_key()` would flip rank detection off its launcher-env fallback as a side effect
   of building a cache config, electing *every* rank of an EME-only job as master — the
   item 4 bug. Fingerprinting from the solver keeps the promise, since both are loaded by
   then.

   *`test_cache_multirank.py::TestRankDetection` was order-dependent.* It exercises
   `_is_master()`'s third-priority launcher-env fallback but never established that the
   first two were absent, so any earlier test importing meow silently loaded meep and
   `am_really_master()` answered instead — True in a serial pytest process, making the
   assertions vacuous. Its `serial_env` fixture now hides meep and mpi4py from
   `sys.modules`, so the fallback is genuinely covered regardless of import order. Note
   the docstring claim that "the EME-only path imports neither meep nor mpi4py" is
   **false** in this environment; production is unaffected because meep initialises MPI
   at import under srun and answers correctly, but the fallback is dead there.

   **Step 1 done, 2026-07-28. Step 2 attempted, measured, and reverted.**

   Step 1 was **not** a missing-interpolation bug: meow's `SampledMaterial.__call__`
   already interpolates correctly, and `meow/fde/meep.py::meep_material()` was *bypassing*
   it to re-implement `argmin` on the raw table, with `MaterialStack.si_index` doing the
   same. Both copies are gone.

   Rather than patch the two lookups separately, both engines are now driven from one
   Sellmeier model in `beef/materials/stack.py::_SELLMEIER` — Salzberg 1957 for Si,
   Malitson 1965 for SiO2, lifted from meep's own library and verified to reproduce it to
   4e-16. The EME gets a densely sampled `SampledMaterial` (0.5 nm, interpolation error
   ~1e-12); the FDTD evaluates the same model at its band centre. **Measured agreement
   across 1.26–1.70 µm: 4e-16** — machine precision, where before one engine interpolated
   a table and the other froze a snapped index. The C-band step in n(λ) is gone.

   **Step 2 — a dispersive FDTD medium — is still open, now for a measured reason.**
   Building the media from the Sellmeier susceptibilities was implemented and then backed
   out. meep subpixel-averages only the *instantaneous* ε and leaves susceptibility terms
   unsmoothed, and the Sellmeier form (like meep's own library materials) puts ε∞ = 1 with
   the entire index in susceptibilities — so meep builds a structure with
   `epsilon diagonal = (1,1,1)`, no dielectric contrast at the Si/SiO2 interface, and
   subpixel smoothing does nothing at all:

   | 450×220 nm Si strip in oxide, res 40, 1.55 µm | n_eff |
   |---|---|
   | non-dispersive `mp.Medium(index=n)` | 2.313398 |
   | dispersive, ε∞ = 1 + Lorentzians | 2.363038 |
   | | **+0.049640** |

   That trades the ~0.002 the dispersion was worth for a 0.0496 discretization loss on the
   same quantity — **25× the wrong way**. `notes.md` recommendation #9 predicted exactly
   this before the attempt ("budget a re-run of the resolution sweeps as part of #5, not
   as a follow-up"), and it held.

   The fix is not to abandon dispersion but to re-parameterise so ε∞ carries the bulk of
   the permittivity and the susceptibilities carry only its wavelength dependence — then
   smoothing acts on the real contrast and the residual dispersive terms are small enough
   to leave unsmoothed. That needs its own fit *and* its own resolution sweep, so it is
   not bundled here. `meep_medium()` stays in `materials/stack.py`, unused by the solvers,
   as the correct model to start from.

   **That prescription has now been checked against meep's actual `LorentzianSusceptibility`
   shape, and the checked version — the fit band, the λ₀ degeneracy, the strip-test gate —
   is recorded in [`material_model_validation.md`](material_model_validation.md), since it
   reaches every device, not just this one.**

   **2026-07-31 — landed, gate-passed, and the resolution question closed on this device's
   own cross-section.** The re-parameterisation shipped as `dispersive=True` (opt-in, off by
   default) and passed the strip-test gate on both stacks/bands. Job 11588997 then swept
   res 20→120 on **this coupler's own body cross-section** under both settings: at band
   centre dispersion changes nothing (identical convergence resolution throughout); at the
   band edge the 1e-3 threshold moves from res 60 to res 80 under dispersion, which is
   exactly the resolution item 11 above already pins production at — so nothing currently
   shipped needs to change. Full table in `material_model_validation.md` work item 4. Step 2
   is therefore no longer simply "open" — it exists, is gated, and its one open question
   (this device's own resolution sensitivity) has an answer; turning it on by default for
   any example is still a separate, unstarted decision.

   **Three findings worth keeping.**

   *meep's `SiO2` is the wrong material.* It is crystalline α-quartz (Ghosh 1999) —
   birefringent, "ordinary axes in X and Y, extraordinary in Z" — reading 1.4578 at
   1.55 µm against fused silica's 1.4440. Nothing on an integrated chip is crystalline
   quartz: BOX is thermal oxide, cladding is deposited oxide, both amorphous. Using that
   name would have introduced a 0.0138 cladding error, four times larger than the effect
   being fixed. `fused_quartz` is the right one, and it matches meow's oxide table to
   **0.00000** in band — they are the same Malitson model, so the cladding did not move.

   *The media are built from raw susceptibilities, not `meep.materials.Si`.* That object
   carries `valid_freq_range` 1.36–11 µm and examples here run the O-band to 1.26 µm, so
   it would raise. Salzberg's 1.36 floor is the 1957 paper's *measurement* range, not a
   pole — the nearest resonance in the fit is at 0.3015 µm — and at 1.26 µm the model is
   0.001 from meow's own table, against a 0.0028 disagreement in the C-band where both
   are nominally valid. The extrapolation is explicit and documented.

   *Fitting meow's own Si table instead was tried and rejected.* Over 1.2–2.0 µm it offers
   8 noisy points (1.372 → 1.400 µm drops 0.0131 while 1.400 → 1.532 drops only 0.0092),
   and least squares returns a near-cancelling A = −1085.9 / B = +1097.5 with a **0.0061**
   max residual — worse than the 0.0028 discrepancy it would have replaced. Hence
   adopting Salzberg, which does shift n_Si by ~0.002: a deliberate source change, not a
   refactor.

   **The cache fingerprint was extended to match**, since one midband index no longer
   describes a dispersive run: `resolved_material_config()` now records n at the band
   edges *and* centre, so a model with the right midpoint but the wrong slope still moves
   the key.

   **Sequencing.** Step 0 alone invalidates the whole on-disk cache without changing any
   number, so it should land together with step 1 (and ideally step 2) and be re-run once
   — bundled with the `COUPLING_LENGTH` regeneration already outstanding under item 14,
   which needs fresh full-device FDTD anyway.

15. **`L = 4.0 µm` hangs in the EME half of `aggregate()`.** *(opened 2026-07-27.)*
    Reproducible: res 20/30/50 (job 11502074 tasks 1/5/9, all killed at the 3h12m wall)
    and previously job `11490488_1`. The FDTD halves of those same tasks completed in
    1.6/3.4/22.5 min, so the hang is specifically in the EME half that `aggregate()` runs
    redundantly on all 48 ranks. Cancelling the three redundant tasks saved ~259
    core-hours. Not diagnosed. Nothing depends on it, but it is a real hang in a code path
    every array sweep executes and it has now appeared in two separate jobs.

### Unblocked cleanups, deliberately not bundled

- **Delete the items 1–2 harness.** `verify_items_1_2.py`, its batch file, and the
  Makefile's `verify_items_1_2` / `clean_verification` targets were to go once item 6
  closed. It has. Deserves its own commit.
- ~~**Drop `N_MODES` 6 → 5** on the coupler — one of the six interface modes is a radiation
  mode (check 2 above).~~ **DONE 2026-07-30** as part of
  `fdtd_vs_EME_mode_validation.md` item 9's cross-device pass (band 5, n_eff
  1.4367, below the 1.4443 cladding), across `directional_coupler_hybrid.py`,
  `_stub_sweep.py` and `_resolution_sweep.py`. This entry was left standing on
  this list until 2026-07-31; the code had already moved.
- **Rank-gate `aggregate()`.** The calibrated re-run wrote *zero* new cache entries and
  still took 34 min/task, because all 48 ranks redundantly unpickle ~55 MB and redo the
  `[200 × 11 × 6 × 6]` einsum and the plotting. Cached is not the same as cheap. Same
  rank-gating pattern item 4 applied to cache writes.
- **Repoint the ~76 `dc_validation.md` citations** at this file or the archive.

## Verification

### Outstanding

**Regenerating the example outputs — not yet run.** `directional_coupler_fdtd.py`, both
resolution sweeps and both notebooks still carry results computed at 13.5 µm. Fresh
full-device 3-D FDTD at 15.13 µm, so nothing hits the cache.

**Item 13 — not scheduled.** A bare 450×220 nm Si strip: one guided mode, no even/odd
splitting, so none of the `β/Δβ ≈ 99` amplification that makes the coupler hypersensitive
— compare `section_mode_beta` against `section_mode_neff` on a geometry with a published
n_eff. Half the harness exists: `section_mode_neff` is already pinned there at 2.3124.
Add a third discretization via `EMESolver(mode_solver="tidy3d")`, a one-line diagnostic on
an already-supported path; if it lands with MPB, that is two operators against one.

**Items 8 and 15 — no harness yet.**

### Standing

- **Fast suite: 271 passed, 1 skipped as of this file's 2026-07-28 date; 489 passed, 1
  skipped as of 2026-07-31** (growth is unrelated library work across other devices, not a
  regression here). Includes `test_beta_calibration.py` (26 tests on
  the β extraction and its two guards), `test_atomic_write.py` (9, on the concurrent-write
  fix), `test_mode_identity.py`, `test_engine_agreement.py`, `test_eme_convergence.py`.
- **Goldens: 4, and all four are now STALE — item 8 invalidated them.** Every golden is
  `@pytest.mark.slow`, so the fast suite deselects them and has *not* been checking them
  through any of this work. Item 14 could not move them by construction: the runner is a
  replica calling `HybridSolver.solve()` on its **defaults**, and `calibrate_body_neff`
  defaults off. Item 8 is different — it changed the silicon index model for everything,
  shifting n_Si by ~0.002, so every golden's numbers move and **all four must be
  regenerated**. They are unverified until that happens. Bundle it with the
  `COUPLING_LENGTH` regeneration under item 14; the cache is empty, so it is one cold
  pass. (If `calibrate_body_neff`'s default ever flips too, replica and json regenerate
  together.)
- **Full-device FDTD reference:** `example_output/verification/ground_truth_sweep/length_sweep.npz`
  — 0.5 µm grid across the full band, provenance in the file. L₅₀ = 15.13 µm (band mean).
- The engine-agreement measurement itself is re-runnable:
  `examples/verification/python/verify_item_6_engine_agreement.py`, with `--res` for the
  resolution sweep.
