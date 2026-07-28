# Directional coupler validation

Status, open work, and outstanding verification. **Everything already closed — items 1–7
and 9–12, with their measurements, jobs and reasoning — lives in
[`archive/dc_validation_20260727.md`](archive/dc_validation_20260727.md).**

> **This file was `dc_validation.md` until 2026-07-27.** Roughly 76 references across the
> library, tests and examples cite it by that name, usually as *"dc_validation.md item N"*
> — and almost always a *completed* item. Every numbered item and its evidence is
> preserved in the archive above, so those citations still resolve in one hop. Item
> numbers are stable IDs and are never reused.

## Status — 2026-07-27

**The hybrid chain converges, and the coupler's hybrid results are trustworthy.**

Item 6, the last obstacle and the reason the original document existed, is closed. The
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

**Read that as self-consistency, not proven physical accuracy.** The reference is meep's
own discretization, and MPB's Δβ disagrees with it by ~5% with *both* engines looking
converged. Which is closer to Maxwell is open as **item 13**, and that caveat belongs
wherever the 0.11% is quoted.

### Where things stand per check

The four checks are defined in [`README.md`](README.md).

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

**Open: 8, 13, 14, 15.** Item 14 is the one that still affects a user of this repo.

14. **The flagship tutorial does not use the fix.** *(opened 2026-07-27 — in progress.)*
    `examples/directional_coupler/python/directional_coupler_hybrid.py` auto-picks and
    prints a 50:50 split point, and until now hand-rolled the *uncalibrated* chain — so it
    reported ≈14.5 µm where the full-device FDTD says 15.13 and the calibrated sweep says
    15.147. **A ~4.5% error in the most-read example in the repo.**

    Item 6 closed on the *method*, proven on the stub sweep, but nothing under `examples/`
    consumed it. Adoption was deliberately deferred for the MMIs (they need `N_MODES`
    12 → 8 first); the coupler tutorial fell through the same gap despite being the device
    the fix was validated on. Its check 2 already passes, so nothing blocked it.

    **Change made, verification queued** — `CALIBRATE_BODY = True` plus `neff_override`
    into `solve_section_length_sweep`, and a log line contrasting MPB's mode index with
    the FDTD's own phase. Cheap: `section_mode_beta` cuts at `mid ± L/2` and section
    polygons inside the uniform coupling region do not depend on the device's total
    `coupling_length`, so the three calibration solves hit the cache the stub sweep
    already filled, despite this script's `half_body` differing (15.125 vs 15.0).

    Not yet confirmed — see Verification below.

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
- **Drop `N_MODES` 6 → 5** on the coupler — one of the six interface modes is a radiation
  mode (check 2 above).
- **Rank-gate `aggregate()`.** The calibrated re-run wrote *zero* new cache entries and
  still took 34 min/task, because all 48 ranks redundantly unpickle ~55 MB and redo the
  `[200 × 11 × 6 × 6]` einsum and the plotting. Cached is not the same as cheap. Same
  rank-gating pattern item 4 applied to cache writes.
- **Repoint the ~76 `dc_validation.md` citations** at this file or the archive.

## Verification

### Outstanding

**Item 14 — job 11526547, queued.** Runs the tutorial with `CALIBRATE_BODY = True`. The
prediction is falsifiable and recorded in the batch header: the auto-picked 50:50 length
moves **≈14.5 → ≈15.15 µm**, i.e. from ~4.5% off the full-device reference to well under
1%. If it does not move, `neff_override` is not reaching the propagation.

**Item 13 — not scheduled.** A bare 450×220 nm Si strip: one guided mode, no even/odd
splitting, so none of the `β/Δβ ≈ 99` amplification that makes the coupler hypersensitive
— compare `section_mode_beta` against `section_mode_neff` on a geometry with a published
n_eff. Half the harness exists: `section_mode_neff` is already pinned there at 2.3124.
Add a third discretization via `EMESolver(mode_solver="tidy3d")`, a one-line diagnostic on
an already-supported path; if it lands with MPB, that is two operators against one.

**Items 8 and 15 — no harness yet.**

### Standing

- **Fast suite: 271 passed, 1 skipped.** Includes `test_beta_calibration.py` (26 tests on
  the β extraction and its two guards), `test_atomic_write.py` (9, on the concurrent-write
  fix), `test_mode_identity.py`, `test_engine_agreement.py`, `test_eme_convergence.py`.
- **Goldens: 4, unchanged.** `tests/golden/directional_coupler_hybrid.json` pins
  `HybridSolver.solve()` on its **defaults**, and `calibrate_body_neff` defaults off. The
  golden runner is a replica, not a call into the tutorial, so the item 14 change cannot
  move it — deliberately. If that default ever flips, replica and json regenerate
  together.
- **Full-device FDTD reference:** `example_output/verification/ground_truth_sweep/length_sweep.npz`
  — 0.5 µm grid across the full band, provenance in the file. L₅₀ = 15.13 µm (band mean).
- The engine-agreement measurement itself is re-runnable:
  `examples/verification/python/verify_item_6_engine_agreement.py`, with `--res` for the
  resolution sweep.
