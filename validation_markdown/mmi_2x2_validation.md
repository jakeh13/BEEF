# 2×2 MMI validation — to-do

> **2026-07-28 — the material model changed under every device here.** Both engines
> now evaluate one analytic Sellmeier model (Salzberg Si, Malitson SiO2) instead of
> each doing its own nearest-tabulated lookup, which shifts n_Si by ~0.002 and
> removes a 0.0063 step that sat in the middle of the C-band. Any golden or stored
> result for this device predates that and needs regenerating; see
> `directional_coupler_validation.md` item 8 and, for the still-open dispersive-FDTD
> half, `material_model_validation.md`.



Opened 2026-07-26. Nothing is known to be wrong; this is the gap list against the
four checks in [README.md](README.md).

**Status 2026-07-26:** check 2 (mode identity) is **done and passed** — job
11500536, see below. It was expected to be the riskiest device in the repo and
turned out to be the cleanest. Checks 1, 3 and 4 remain unrun, and there are
still no goldens.

**Status 2026-07-30:** two config changes and one new doubt. `N_MODES` is now
**14**, not 15 (gap 3, done — the unguided column is dropped). Gap 5 is in
flight as jobs 11577985/11577986. And gap 6 is **new**: the passing check-2
result above was measured before either duplicate-mode guard existed, on a
mode count this device no longer requests, so it needs a cheap re-run before it
can be quoted as current. Checks 1 and 3 remain unrun; still no goldens.

`examples/mmi_2x2/python/` — `mmi_2x2_eme.py`, `mmi_2x2_fdtd.py`,
`mmi_2x2_hybrid.py`.

## Where it sits

Structurally the same shape as the 1×2 MMI — three scripts, EME/FDTD/hybrid —
but with **no golden fixtures at all**. `tests/golden/` holds only
`directional_coupler_hybrid.json` and the three `mmi_1x2_*.json`. So unlike the
1×2, nothing here is pinned: any change to the solvers moves this device's
numbers silently.

That is the one thing this device has that the 1×2 does not have a version of:
the 1×2's problem is *stale* references, which at least fail loudly under
`--runslow`. Having none is quieter and worse.

## Gaps

1. **No goldens.** Add `mmi_2x2_eme`, `mmi_2x2_fdtd`, `mmi_2x2_hybrid` targets to
   `tests/test_golden_examples.py`'s `TARGETS` and generate them with
   `scripts/generate_golden_values.py`. Do this *after* gap 2, so the recorded
   values are not immediately invalidated by a resolution change — recording a
   golden against a configuration you already intend to change is how
   `directional_coupler_hybrid.json` ended up regenerated twice.

2. **Shared resolution knob — and the EME side is not converged at any
   resolution tried.** `mmi_2x2_hybrid.py:69` passes one `RESOLUTION = 40.0` to
   both sub-solvers, the defect `directional_coupler_validation.md` item 11 fixed for the
   coupler.

   The audit's mesh sweep (same job) gives `n_eff[0]` on this cross-section:

   | res | 40 | 60 | 80 | 100 |
   |---|---|---|---|---|
   | n_eff[0] | 2.806984 | 2.819049 | 2.825888 | 2.828554 |
   | step | — | +0.0121 | +0.0068 | +0.0027 |

   Still moving 0.0027 at 80→100, and the steps shrink by only ~2× each time,
   so this is a slow tail rather than a converged value. **The coupler's res 80
   does not transfer** — that pin was justified in item 11 on Δβ, the even/odd
   *difference*, which converges far faster than absolute n_eff. All three
   audited devices are still moving in absolute n_eff at res 100.

   Whether that matters depends on which quantity this device's result actually
   depends on. Establish that first — the coupler's lesson was that sweeping a
   quantity the physics is not sensitive to (T_cross) hid a 7% error in one it
   is (Δβ). Do not simply raise the resolution and assume it is fixed.

3. ~~**Mode identity — this is the most exposed device in the repo.**~~
   **Done 2026-07-26 — passed, and the prediction behind this gap was wrong.**

   Measured on the uniform body cross-section (`x ∈ [−1, +1]`, 13 µm clear of
   either taper seam) at 1.55 µm, both engines meshed at res 40, job 11500536:

   | | |
   |---|---|
   | modes requested | 15 |
   | **guided** (n_eff > 1.4443) | **14** |
   | correspondence | identity, `[0 … 13]` |
   | worst guided \|Δn_eff\| | **0.0007** |

   The reasoning in this gap — that 15 modes crowd together and so this device
   is the most likely to reorder — **did not hold**. It is the cleanest of the
   three devices audited, better than the coupler's 0.0023 on 6 modes and an
   order of magnitude better than the 1×2's 0.0067 on 12. Mode *count* is not
   what drives disagreement; confinement is. This device is SOI_220 on a 2.8 µm
   body, and even band 13 sits at n_eff 1.508, comfortably above the 1.4443
   cladding. The 1×2 is SOI_160, and a 160 nm core holds its top bands so
   weakly that they pile up at cutoff — that is where its failure is
   (`mmi_1x2_validation.md` gap 3).

   **Band 14 (n_eff 1.4078) is below the cladding index** — a radiation mode.
   So one column of every 2×2 interface basis is unguided, the same finding as
   the coupler's band 5, and its +0.0021 disagreement is 3× worse than any
   guided band here. `N_MODES = 14` would drop it at no cost.

   **DONE 2026-07-30** — `mmi_2x2_hybrid.py:124` now reads `N_MODES = 14`,
   applied as part of `fdtd_vs_EME_mode_validation.md` item 9 (the cross-device
   pass that drops each device's one unguided column). The configuration
   change this gap asked for has therefore already happened, which moves the
   goldens question in gap 1 / `## Order` forward one step — see there.

4. **Engine agreement not measured.** Adapt
   `verify_item_6_engine_agreement.py`: solve a uniform stretch of the MMI body
   at two lengths and check the two engines agree on each mode's accumulated
   phase. The extraction helper
   (`tests/test_engine_agreement.py::dbeta_from_phases`) is device-independent
   and already tested — it takes lengths and phases, nothing coupler-specific.

5. **No hybrid-vs-full-device comparison on record.** `mmi_2x2_fdtd.py` exists,
   so the reference run is available; nothing asserts the hybrid reproduces it.
   This is the only check that tests the cascade as a whole.

   **In progress 2026-07-30** — `examples/verification/python/
   verify_mmi_2x2_hybrid_vs_fdtd.py` added: loads both scripts' saved
   `SMatrix` pickles, compares bar/cross transmission per wavelength
   (reconciling `mmi_2x2_fdtd.py`'s `out_0`/`out_1` port names against
   `mmi_2x2_hybrid.py`'s `bar`/`cross`), and writes a provenance-carrying
   report (unlike the coupler's original check-4 reference — see
   `README.md`). Reports the numbers only, no pass/fail threshold: gap 4
   (engine agreement) is not measured for this device yet, so there is no
   basis to say how much of any observed gap is calibratable
   (`resolved/directional_coupler_validation.md` items 6/14 precedent: ~5%
   uncalibrated, 0.11% after `calibrate_body_neff`) vs. structural.
   SLURM job 11577985 (array indices 1-2, `mmi_2x2_fdtd`/`mmi_2x2_hybrid`)
   generates the two SMatrix files; job 11577986
   (`--dependency=afterany:11577985`) runs the comparison once they land.
   Note this runs *before* gap 2 (resolution) despite the Order section
   below — both examples still share `RESOLUTION = 40.0`, so this first
   pass measures the hybrid at the same not-yet-converged resolution gap 2
   already flags, not a clean result.

6. **New, opened 2026-07-30 — this device's mode-identity pass predates both
   duplicate-mode guards.** Gap 3's result (job 11500536, 2026-07-26) is the
   evidence for "cleanest device audited", and it is quoted that way in
   `README.md`'s Finding 1/2 and in three other files. It was measured before
   `FDTDSolver._iter_distinct_bands()` and `EMESolver.compute_modes()`'s
   `guard_duplicate_modes` existed (both 2026-07-30), so:

   - **The audited basis is not the basis the device now builds.** With the
     guards on, a duplicate at any band makes the solver reach past it, so the
     15 (now 14) modes returned need not be the 15 that were audited.
   - **This is the repo's highest-mode-count SOI_220 device**, and
     `fdtd_vs_EME_mode_validation.md` item 11 shows the meep-side guard's
     duplicate check has no discriminating margin left at these mode counts on
     a 2.8 µm body — which is exactly this device's body width. If that guard
     mis-fires anywhere, here is a likely place.
   - `N_MODES` also moved 15→14 since the audit (gap 3), so the recorded table
     describes a request this device no longer makes.

   Cheap to close: re-run `audit_interface_modes.py mmi_2x2` at `N_MODES = 14`
   with the guards on and confirm the identity permutation and worst
   |Δn_eff| 0.0007 still hold. Any guard firing in that log is itself the
   finding. Until then, treat "the cleanest device audited" as a 2026-07-26
   measurement of a superseded configuration, not a current claim.

## Order

~~**3** first~~ — done 2026-07-26, passed. Remaining:

1. **Gap 6** (re-audit at `N_MODES = 14` with the duplicate guards on) — new
   2026-07-30 and now first, because it is a smoke run and because everything
   below quotes gap 3's numbers as though they still describe this device.
2. **Gap 2**, the resolution question — it gates gap 1, and the audit has
   already shown res 40 is not converged.
3. **Gap 5** (hybrid vs full-device FDTD), which needs no new scripts.
   *In flight as of 2026-07-30 — jobs 11577985 / 11577986, see gap 5.*
4. **Gap 1** (goldens) once the configuration has stopped moving. ~~Dropping
   `N_MODES` to 14 per gap 3 counts as the configuration moving, so do it
   before recording anything.~~ **That drop is done** (gap 3, 2026-07-30), so
   the remaining blocker on recording goldens here is gap 2's resolution
   question alone — plus the repo-wide material-model regeneration
   (`README.md` Finding 6).
5. **Gap 4** last, as on the coupler.
