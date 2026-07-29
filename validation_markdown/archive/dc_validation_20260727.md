# Directional coupler validation — to-do

Opened 2026-07-25, from investigating why `directional_coupler_hybrid_stub_sweep.py`
never converges. Root cause found; this is the work to fix it and re-establish trust
in the directional coupler results.

**Status 2026-07-27: the stub sweep converges.** Items 1–7 and 9–12 done and verified;
**8, 13, 14, 15 open**. Item 6 — the last obstacle, and the reason this document stayed
open — is **closed**. The one that still bites a user is **item 14**: the fix is proven
and the flagship tutorial does not use it, so `directional_coupler_hybrid.py` still
reports a 50:50 length ~4.5% off the full-device reference.

The fix was not to reconcile the two engines but to stop asking the EME body to
propagate on a β the FDTD reference does not have. `HybridSolver(calibrate_body_neff=True)`
measures β from the FDTD's own accumulated phase and propagates the body on that. Job
11517245:

| | drift slope (µm/µm) | L₅₀ spread over 9 stubs | mean L₅₀ | vs 15.13 µm reference |
|---|---|---|---|---|
| uncalibrated (job 11500412) | +0.108543 | 0.220 µm | 14.429 | **−5.37%** (at stub 0) |
| **calibrated** | **−0.000651** | **0.005 µm** | **15.1468** | **+0.111%** |

Both predictions the test was registered against held. The drift is **167× smaller** and
below the run's own 0.0015 µm rms scatter, so the cut point no longer measurably matters
— which was the entire design goal of a hybrid decomposition. And L₅₀ landed at 15.147
against a *prediction* of 15.02, i.e. closer to the independent full-device FDTD than the
algebra forecast.

**Read this as self-consistency, not proven physical accuracy.** The hybrid now reproduces
the full-device FDTD to 0.11%, but that reference is meep's own discretization, and MPB
disagrees with it by ~5% on Δβ with *both* engines looking converged. Which one is closer
to Maxwell is genuinely unsettled and is now **item 13**.

Getting there ran through the EME side first. Item 6's headline gap was mostly a
**measurement artifact**: at the res=40 every example inherited, Δβ_EME was overstated by
7%. Item 11 split the two engines' resolutions and closed that. Item 12 checked the
assumption underneath it all — that both engines return interface modes in the same order
— and for this device **they do** (identity permutation, worst guided-band disagreement
0.0023), so `HybridSolver`'s contraction is sound here. It also found that one of the
coupler's six interface modes is not guided at all.

**Open: item 8** (the material half, worth ~0.8% of Δn_eff) and **item 13** (which
engine's Δβ is physically right). Neither blocks using the coupler's hybrid results.

## What's wrong

> Item 1 is **fixed** as of 2026-07-25 — see Done at the end. This section is kept
> as the description of the defect the rest of the list is scoped against.

**`FDTDSolver` mislabelled its own wavelength axis.** `solvers/fdtd.py:352` and `:654`
computed

```python
fwidth = (1.0/wavelengths.min() - 1.0/wavelengths.max()) * 0.5
```

and pass it as meep's `df` to `add_mode_monitor(fcen, df, nfreq)` (`:411`, `:424`,
`:703`, `:708`). Meep spaces `nfreq` frequencies evenly over `[fcen - df/2,
fcen + df/2]`, **ascending in frequency**. For the standard 1.50–1.60 µm, 11-point
C-band sweep that is:

| index | code labels it | meep actually computed | error |
|---|---|---|---|
| 0 | 1.5000 µm | 1.5754 µm | +75.4 nm |
| 5 | 1.5500 µm | 1.5500 µm | 0 |
| 10 | 1.6000 µm | 1.5254 µm | −74.6 nm |

Two independent defects in one line:

1. **Half the span.** `df` is set to *half* the requested frequency width, so the
   monitor covers 1.5254–1.5754 µm, not 1.50–1.60 µm.
2. **Reversed.** Meep's grid ascends in frequency, i.e. *descends* in wavelength, but
   the result is indexed against the caller's ascending-wavelength array.

Only the centre point (index 5, `fcen`) is correctly labelled. Verified directly
against `mp.get_flux_freqs()`.

Note a linear-in-λ array can never be reproduced by a `(fcen, df, nfreq)` triple —
meep spaces linearly in frequency. The fix must pass an explicit frequency list.

### What this explains

- **Every FDTD directional-coupler result shows backwards wavelength dispersion.**
  Archived full-device sweep: 50:50 length 13.16 µm at "1.50" rising to 17.65 µm at
  "1.60". Coupling appears to *weaken* with wavelength, which is backwards for a
  directional coupler. Unflip the axis and it is the physically correct direction.
  Reproduced in the archived FDTD *resolution* sweep at all 8 resolutions, so it is
  not a mesh artifact.
- **Why the stub sweep can't converge.** `np.einsum("Lwij,wjk->Lwik", M_body, M_in0)`
  contracts the EME body at true λ_w against an FDTD section at a *different*
  wavelength — at w=0 it multiplies EME@1.500 µm by FDTD@1.5754 µm. Ten of eleven
  columns are physically meaningless products. Each stub is its own FDTD run, so the
  wavelength-mean the sweep plots drifts smoothly with stub with no reason to ever
  settle: `L_50 = 0.2427*BODY_STUB + 13.4067 µm`, max residual to that line 2.2e-3 µm
  over 0→2 µm. Zero curvature.

### What it does *not* explain

Index 5 is correctly labelled in both engines, and the hybrid still disagrees with the
full-device FDTD there: 13.65 µm (stub 0) → 14.09 µm (stub 2) versus **15.29 µm**.
Restricting the stub sweep to λ=1.55 alone still gives a dead-straight
+0.2185 µm/µm drift. In the two-mode picture,
`L_50(s) = [π/2 − 2φ_bend]/Δβ_EME + 2s(1 − Δβ_FDTD/Δβ_EME)`, so that slope measures
`Δβ_FDTD/Δβ_EME = 0.891` — the FDTD sections accumulate ~11% less even/odd phase per
µm of the *same* straight coupled region than the EME body does. That is a second,
independent problem (item 6) — **since closed, 2026-07-27**, by propagating the EME
body on the FDTD's own β rather than by reconciling the two engines.

Correction to an earlier read of this: the narrow `fwidth` does **not** starve the
source. Source and monitor share the same variable, so the pulse is well matched to
the band meep actually sampled — power conservation is 0.9987–0.9998 at every point.
The bug is the labelling, not the signal level.

## Is the archived FDTD run a valid ground truth?

`archive/example_output_20260717/directional_coupler/directional_coupler_fdtd_sweep/length_sweep.npz`

**Partly — as a single-wavelength reference only.**

- Usable: power-conserving (0.9987–0.9998); resolution-converged — the 2026-07-22
  FDTD resolution sweep gives T_cross = 0.423 ± 0.008 at λ=1.55 across res 30–90, so
  res=40 is fine; and **its index 5 is correctly labelled**, giving a genuine
  full-device 50:50 length of **15.29 µm at 1.55 µm**.
- Not usable as-is for spectra: the same axis bug mislabels its other ten wavelength
  points. **They are recoverable, not garbage** — job 11477911's check 8 validated the
  relabelling recipe (map index *i* to `1/(fcen + df*(i/(n-1) − 1/2))`, i.e. reversed
  onto 1.5754→1.5254 µm) on a sibling run from the same code, agreeing with a fresh
  post-fix solve to 2.5e-4. So this file can be re-read on its true axis if needed.
  The catch is that its true axis spans only 1.5254–1.5754 µm, half the nominal band.
- Provenance is unrecorded. The npz stores only `coupling_lengths`, `wavelengths`,
  `T_through`, `T_cross` — no config, no solver fingerprint, no geometry, and the
  whole-device `solve()` path had no `save_dir` then (the archived cache holds only
  `solve_section`/`solve_section_length_sweep` entries). Its parameters match the
  current script's constants, but nothing in the file proves it.
- Coarse: 2 µm length steps, so 15.29 µm is interpolated across a 2 µm bracket
  (~±0.1 µm), against a stub sweep resolving 0.06 µm steps.

Verdict: **superseded — item 5 is done.** Ground truth is now
`example_output/verification/ground_truth_sweep/length_sweep.npz`, on a 0.5 µm grid
across the full band with provenance in the file. This archive is kept only as an
independent cross-check, and in that role it did well: re-read on its true axis, all 11
of its points land on the new L₅₀(λ) curve **within 0.027 µm** (a uniform +0.023 µm
bias, 0.003 µm scatter). Its 15.29 µm at 1.55 against the new 15.31 µm is the sharpest
single confirmation that provenance and the 2 µm step — not the axis bug — were what
kept it from being truth.

To re-read it, the axis is `1/(fcen + df·(i/(n−1) − 1/2))` with **`fcen = 1/mean(λ)`**
and `df = (1/λ_min − 1/λ_max)/2` — giving 1.5754 µm at index 0 down to 1.5254 µm at
index 10. Using the frequency-mean `0.5(1/λ_min + 1/λ_max)` instead shifts every point
~1.6 nm and biases L₅₀ by ~0.17 µm; see the correction in Done item 5.

## To-do

Numbers are stable IDs — completed items move to Done and keep their number, so
cross-references stay valid.

**Open: 8, 13, 14, 15.** Item 6 — for a while the only thing left — is closed as of
2026-07-27. **Item 14 is the one that matters**: the fix is proven but the tutorial does
not use it, so the repo's headline coupler example still reports a number 4.5% off. 13 is
the honest unknown behind item 6's closure, 8 and 15 are minor. The ordering below is kept
because it records *why* the sequence was chosen and how it turned out.

The plan was **11 → 12 → 7**, then **8**. Both prerequisites paid off, and the payoff is
measurable:

- **11 before 7** because the sweep's drift slope measures `2·(1 − Δβ_FDTD/Δβ_EME)`, and
  the EME half of that ratio was 7% off at res 40. Running 7 first would have spent nine
  array tasks measuring a drift that was more than half discretization error. **Borne
  out**: the slope fell 0.2427 → 0.1085, and the post-fix value matches the item 6 +
  item 11 prediction to 0.82%. On the pre-fix mesh it would not have.
- **12 before 7** because if the two engines ordered their six interface modes
  differently, the chain the sweep tests would be contracting mismatched modes and no
  stub could ever fix it. **Excluded** — identity permutation, worst guided disagreement
  0.0023. Worth having done: the same audit later found this failure *is* real on
  `mmi_1x2`, where meep returns two degenerate bands that meow splits. Had the coupler
  behaved that way, item 7's straight line would have had two candidate explanations
  instead of one.
- **8 after**, worth ~0.8% of Δn_eff (measured, Done item 6) — real, but a minor
  contributor to the ~5% that remained.

One thing the sequence got wrong, recorded because the mistake is instructive: **11 and
12 were both scoped as ways to make the two engines agree**, and neither could have.
Item 12's own measurement is what eventually showed why — both engines run MPB, so
aligning them was never worth more than 0.57% against a 10% gap. The fix that worked
came from the opposite direction: stop trying to make EME match FDTD's *mode solve* and
make it match FDTD's *propagation*. See the comparison under Done item 12.

**10** is done; the ring's *re-run* is tracked in `ring_resonator_validation.md`, not
here.

6. *(done — see Done item 6. **Closed 2026-07-27** by β calibration; the stub sweep
   converges.)* Two follow-on items came out of the original measurement — 11 and 12
   below — and a third, **item 13**, out of the fix.

11. *(done — see Done item 11.)*

12. *(done — see Done item 12.)*

7. *(done — see Done item 7. It did **not** converge, and that answer moved the work to
   item 6 — where it was fixed. Re-running the same sweep with
   `calibrate_body_neff=True` flattens it: slope +0.1085 → −0.0007. The original
   sentence is left standing above because it records what was true, and correctly,
   when the sequence was chosen.)*
13. **Which engine's Δβ is physically right?** *(opened 2026-07-27 by item 6's fix.)*
    Item 6 achieved self-consistency: the hybrid now reproduces the full-device FDTD to
    0.11%. But that reference is meep's own discretization. MPB plateaus at Δβ = 0.0894
    by res 80–100; Δβ_FDTD sits flat at ~0.085 across res 20–50. **Two independently
    converged answers ~5% apart** means at least one carries a systematic that a
    resolution sweep cannot see, and nothing measured so far says which.

    This matters because **L₅₀ ∝ 1/Δβ almost exactly here**, and the numbers say so
    rather than the theory:

    | ratio | value |
    |---|---|
    | L₅₀: full-device FDTD 15.13 ÷ uncalibrated hybrid 14.429 | 1.04858 (+4.86%) |
    | Δβ: converged MPB 0.0894 ÷ FDTD 0.085212 | 1.04915 (+4.91%) |

    They agree to **0.057 percentage points**. The entire L₅₀ disagreement between the
    two methods is the Δβ disagreement and nothing else — no bend phase, no seam
    mismatch, no basis truncation left over. So "which L₅₀ is right" is not a separate
    question; it reduces completely to this one. If MPB is right, the *uncalibrated*
    hybrid was closer to physical truth all along and the calibration is faithfully
    reproducing meep's error.

    **The test that would settle it.** A bare 450×220 nm Si strip: one guided mode, no
    even/odd splitting, so none of the `β/Δβ ≈ 99` amplification that makes the coupler
    hypersensitive — compare `section_mode_beta` against `section_mode_neff` on a
    geometry with a published n_eff. Half the harness exists: `section_mode_neff` is
    already pinned there at 2.3124 (Done item 12). Add a third discretization via
    `EMESolver(mode_solver="tidy3d")`, a one-line diagnostic on an already-supported
    path — if it lands with MPB, that is two operators against one.

    Not scheduled. Nothing downstream is blocked on it as long as the framing above is
    stated wherever the 0.11% agreement is quoted.
8. **Minor, while in there — `notes.md` Solver #9's material half.** Confirmed present
   in the data: EME's Δn_eff has a +0.00022 kink at exactly 1.57 µm where meow's Si
   table switches 3.4784 → 3.4710, while the FDTD half stays frozen at 3.4784 for
   every wavelength. Real, but ~0.8% of Δn_eff — worth fixing for correctness, not a
   driver of anything above. (Solver #8, wavelength interpolation, is not involved:
   this script does a full mode solve at all 11 wavelengths.)
14. **The flagship tutorial does not use the fix.** *(opened 2026-07-27 — **DONE
    2026-07-28**, job 11526547, 2h04m, exit 0. Closed and archived here 2026-07-28; this
    entry was updated in place rather than appended, since the completion supersedes the
    open one below.)*

    **Confirmed, and the prediction was exact.** The tutorial's auto-picked 50:50 length
    moved **14.5 → 15.15 µm** (ratio 0.500), against a registered prediction of ≈15.15.

    | | L₅₀ | vs full-device FDTD (15.129 µm) |
    |---|---|---|
    | stale constant 13.5 | 13.500 | −10.77% |
    | uncalibrated tutorial | 14.429 | −4.63% |
    | **calibrated tutorial** | **15.151** | **+0.14%** |

    That +0.14% is 0.022 µm — **0.07 of the tutorial's own 0.299 µm sweep-grid step**, so
    the agreement is at the resolution floor of the measurement and cannot be improved
    without a finer length grid. The calibration numbers it logged (n_eff 2.375691 /
    2.354670, Δn_eff 0.021021) match the stub sweep's to every printed digit, confirming
    the cache hit and that both scripts calibrate the same cross-section.

    `COUPLING_LENGTH` **has been set to 15.13 µm** in
    `directional_coupler_fdtd.py:80`, `directional_coupler_fdtd_resolution_sweep.py:71`,
    `directional_coupler_hybrid_resolution_sweep.py:85` and
    `notebooks/directional_coupler_fdtd.ipynb`, plus the two module docstrings that quoted
    13.5. The value comes from item 5's ground truth (crossing at 15.129 µm) rather than
    the hybrid's grid-snapped 15.151, since three of the four users are full-device FDTD
    runs and that reference is the finer measurement.

    Deliberately **not** changed: `examples/verification/python/verify_items_1_2.py:131`,
    which also pins 13.5. That harness is a historical artifact — its recorded 9/9 pass
    (job 11477911) *was* run at 13.5 µm, and rewriting the constant would misrepresent
    what it executed. It is due for deletion now that item 6 has closed.

    **Closure is on the code fix and its verification job, not on regenerated outputs.**
    The example outputs for the three scripts above and both notebooks' stored outputs
    were still from 13.5 µm / pre-calibration runs at archival time and need a fresh
    full-device 3-D FDTD pass to update — no cache hits. That regeneration is tracked as a
    standalone chore in `directional_coupler_validation.md`'s To-do (shared with item 8,
    since the emptied cache makes it one cold pass rather than two), not as part of this
    item.

    *(Original entry, kept as the record of why this was opened:)*
    `examples/directional_coupler/python/directional_coupler_hybrid.py` — the repo's
    flagship coupler example, and the one the other devices were written against — still
    hand-rolls the uncalibrated chain. It auto-picks and prints a 50:50 split point
    (`:356`), so it reports **≈14.5 µm** where the full-device FDTD says 15.13 and the
    calibrated sweep says 15.147: a **~4.5% error in the most-read example in the repo**.

    Item 6 closed on the *method* — proven on the stub sweep — but nothing under
    `examples/` consumes it. Adoption was deliberately deferred for the **MMIs** (they
    need `N_MODES` 12 → 8 first, item 12), and the coupler tutorial fell through the same
    gap despite being the device the fix was validated on. Its check 2 already passes
    (identity permutation, worst guided disagreement 0.0023), so nothing blocks it.

    Cheap, too: `section_mode_beta` cuts at `mid ± L/2` and section polygons inside the
    uniform coupling region do not depend on the device's total `coupling_length`, so the
    three calibration solves hit the cache the stub sweep already filled — despite this
    script's `half_body` differing (L_mid ≈ 30.25 µm vs the sweep's 30.0).

    Carry item 13's caveat into whatever the tutorial prints: the calibrated number is
    self-consistent with the full-device FDTD, not independently proven correct.

    **The notebook had drifted on two items, not one.**
    `examples/directional_coupler/notebooks/directional_coupler_hybrid.ipynb` is a
    hand-maintained mirror of the script — `notes.md` Solver #6 flags exactly this risk,
    and notes it already cost a double-application on 2026-07-24. It had missed **item
    11** as well: its `body_solver` still used a single shared `RESOLUTION = 40`, so its
    EME half carried the same 7% Δβ overstatement item 11 fixed in the script on
    2026-07-26. Both are now applied, and every numeric knob in the two files agrees
    (`FDTD_RESOLUTION`, `EME_RESOLUTION`, `WL_MIN/MAX`, `N_WL`, `BODY_STUB`, `N_MODES`,
    `CALIBRATE_BODY`, `PLOT_COUPLING_LENGTH`). The notebook deliberately omits the
    script's MPB-comparison print, because `section_mode_neff()` runs a fresh MPB solve
    and does not read the cache, which would break its `CACHE_ONLY` contract.

    **Why 13.5 µm was wrong**, recorded because the constant sat unchallenged through two
    fixes that should each have caught it: at that length the full-device FDTD gives
    cross/(cross+through) = **0.4318**, a 43:57 splitter, against a reference crossing of
    **15.129 µm** — 10.8% short. It came from the pre-item-1 sweep, which reported
    L₅₀ = 13.65 µm *with the wavelength axis reversed*, and it survived both the item 1
    axis fix and item 5's ground-truth rebuild because nothing rechecked a hard-coded
    number that no test asserted on.

    The two resolution sweeps stayed valid *as convergence studies* throughout — they
    compare one device across meshes, and its length does not decide whether res 40 is
    converged. What was wrong was the claim in their comments, and
    `directional_coupler_fdtd.py` presenting a 43:57 device as a 3 dB coupler.

15. **`L = 4.0 µm` hangs in the EME half of `aggregate()`.** *(opened 2026-07-27.)*
    Reproducible: res 20/30/50 (job 11502074 tasks 1/5/9, all killed at the 3h12m wall)
    and previously job `11490488_1`. The FDTD halves of those same tasks completed in
    1.6/3.4/22.5 min, so the hang is specifically in the EME half that `aggregate()` runs
    redundantly on all 48 ranks.

    Cost so far: the three res-sweep tasks were cancelled after confirming from their
    aggregates (`lengths_completed=4`) that they were redundant, saving ~259 core-hours.
    Not diagnosed. Nothing depends on it — L = 4.0 µm is one point of a convergence
    sweep whose other points all landed — but it is a real hang in a code path every
    array sweep executes, and it has now appeared in two separate jobs.

10. *(done — see Done item 10.)* The **re-run** it unblocks is tracked in
    `ring_resonator_validation.md`, not here.

## Verification

**Items 1–2 have a dedicated harness** (temporary — delete with the items):

    BEEF/examples/verification/python/verify_items_1_2.py
    batch/verification/run_verify_items_1_2.batch   →  make verify_items_1_2

Nine checks, output to `example_output/verification/items_1_2/`
(`verification.json` with per-check PASS/FAIL plus full provenance,
`verification.npz`, and four plots). Exits non-zero on any failure, so the job's exit
status means something. Checks 1–6 are static and need no meep run:

    cd BEEF/examples && OMP_NUM_THREADS=1 python verification/python/verify_items_1_2.py --static-only

### Status: 9/9 pass — job 11477911, 2026-07-25

48 ranks, 56 min, exit 0; meep 1.33.0, beef `c243a46+dirty`, uncached.

| check | result |
|---|---|
| 2. grid matches requested λ | 2.2e-13 nm error (pre-fix: 75 nm) |
| 4. source covers the band | weakest λ at 0.269 of peak (floor 0.10) |
| 7. T_cross rises with λ | dT/dλ = **+3.378 /µm**, 0.276 → 0.613 (pre-fix: backwards) |
| 8. pre-fix data on its true axis | **max \|ΔT_cross\| = 0.00025, rms 0.00014** (tol 0.05) |
| 9. power conserved | T_total 0.9992–0.9997 (tol 0.02); `check_passive()` OK |

Check 8 is the decisive one and it passed by a factor of 200, not by a hair — worth
stating, because a crude straight line would have scored 0.043 against that 0.05
tolerance. Two independent 3-D FDTD runs, one pre-fix and one post-fix with a
different source bandwidth (fwidth 0.0208 → 0.0833), agree to 2.5e-4 in cross
transmission once the old run is placed on the wavelengths meep really sampled. This
is the claim item 1 makes, confirmed at the level of the solver's own convergence:
**the FDTD physics was never wrong, only the labels were.** `dispersion_before_after.png`
shows it directly — the re-plotted pre-fix points sit on the post-fix curve, and the
as-labelled pre-fix curve crosses it at exactly one point, λ=1.55, the single index
both axes agreed on.

Check 8 rebuilds the archived pre-fix config exactly
(`archive/example_output_20260724/.../directional_coupler_fdtd_resolution_sweep/resolution_sweep.npz`,
res-40 row — same device, `coupling_length` 13.5 µm, same wavelengths). The run is
deliberately **uncached** — a verification must not be able to pass on a stored
result, and `save_run()` is not multi-rank safe (item 4).

**Do not misread check 5.** It reports the pre-fix source leaving the band edges at
7.3e-10 of peak, which reads like the old run was noise. It was not, and check 8's
2e-4 agreement is the proof. Check 5 measures the *old* fwidth against the *full*
1.50–1.60 µm request — a counterfactual. Pre-fix, source and monitor shared that narrow
fwidth, so the pulse covered everything meep actually sampled. What check 5 shows is
that item 2 was required *because* item 1 widened the grid: fixing the grid alone would
have starved the new edges. Both were needed; neither alone was sufficient.

Keep the harness until item 7 closes — items 5 and 7 re-run FDTD sweeps and this is
what re-checks them. Then delete the script, the batch file, and the Makefile's
`verify_items_1_2` / `clean_verification` targets.

Remaining items:

- Items 1–2 together: ✅ confirmed by check 7 — the full-device FDTD 50:50 length now
  *decreases* with wavelength (T_cross rises). Any future result where it increases
  means the axis has flipped again.
- Item 4: ✅ done — see Done below. The multi-node check has its own script, since
  pytest cannot reproduce 48 ranks on 2 nodes sharing one cache:

      srun -A <account> -q inferno -N 2 --ntasks-per-node=24 -t 10 --mem-per-cpu=2G \
           --mpi=pmi2 make verify_cache_multirank

  Two nodes minimum and a save_dir on shared storage, or it proves nothing.
- Item 6: the new test fails at ~11% before the fix, passes after.
- Item 7: ❌ **the stated criterion was not met, and the criterion was the wrong one.**
  `L_50` vs `BODY_STUB` did *not* curve or flatten — it stayed straight to 0.0024 µm —
  and the λ=1.55 value did not land on the reference. That is not a failed run: a
  straight line is exactly what `L_50(s) = … + 2s·(1 − Δβ_FDTD/Δβ_EME)` predicts for a
  constant ratio ≠ 1, so "does it flatten?" was only ever a proxy for "do the engines
  agree?". The measurement to keep is the **slope**, which is the ratio, and it now
  agrees with item 6's independent value to 0.82%. See Done item 7.

## Done

### 1. Monitor frequency grid — 2026-07-25

Every DFT monitor in `beef/solvers/fdtd.py` now takes an explicit frequency list
instead of meep's `(fcen, df, nfreq)` triple, and results are reordered back onto the
caller's wavelength array on unpack.

**Changed**

- **`beef/solvers/fdtd.py` — new `monitor_freq_grid(wavelengths)` helper** (module
  level, next to `section_fdtd_layout`). Returns `(mon_freqs, wl_order)`:
  `mon_freqs` is `1/wavelengths` sorted strictly ascending, which is what meep
  requires; `wl_order` is the inverse permutation that lifts a monitor-ordered result
  back onto `wavelengths` as the caller gave it. Works for unsorted input, not just
  ascending λ. Docstring records the defect so the triple doesn't come back.
- **`solve()`** — `nfreq` dropped; `mon_freqs, wl_order = monitor_freq_grid(...)`.
  Both `add_mode_monitor` calls (input normalization + the per-output-waveguide list)
  now pass `mon_freqs`. The normalized amplitude is reordered once at the end:
  `alpha = (coeff.alpha[0, :, 0] / src_alpha)[wl_order]`.
- **`solve_section()`** — same three changes: `nfreq` dropped, both `norm_flux` /
  `mea_flux` monitors take `mon_freqs`, and
  `coeff = (c.alpha[0, :, alpha_idx] / src_alpha)[wl_order]`. Ratio taken in monitor
  order first, reordered once — `src_alpha` and the measured amplitude are both
  monitor-ordered, so they stay aligned.
- `fwidth` is still computed at both sites but now feeds **only** the
  `GaussianSource`. Its value is wrong for that job — that is item 2.

**Added** — `BEEF/tests/test_fdtd_freq_grid.py`, 5 tests:

- monitor frequencies ascending and exactly `1/wavelengths`, with the endpoints
  asserted against the old grid's 1.5254–1.5754 µm band;
- `wl_order` round-trips a monitor-ordered result back onto the caller's array;
- the ascending-λ case is a plain flip, stated explicitly so a refactor can't drop it;
- unsorted input round-trips too;
- end-to-end against meep: `mp.get_flux_freqs()` returns exactly the requested
  wavelengths. Needs a source in the sim (meep allocates field components from
  sources, and rejects a monitor on a fieldless simulation) but no FDTD run —
  the whole file runs in ~2 s.

**Verified** — `pytest tests/test_fdtd_freq_grid.py` 5 passed; full fast suite
`pytest tests/` 130 passed, 6 skipped (the `--runslow` golden tests), no regressions.

**Confirmed physically** — job 11477911 (2026-07-25), checks 7–8: cross transmission
now rises with wavelength, and the archived pre-fix run re-plotted on its true axis
lands on the post-fix curve to `max |ΔT_cross| = 0.00025`. The relabelling account of
this bug is correct; the old data was recoverable, not garbage.

**Not addressed here** — every cached FDTD result predates this fix and is still
mislabelled; the fingerprint does not record the frequency grid, so those caches will
be reused silently. That is item 3, and it blocks trusting any re-run.
*(Superseded by item 2, which changed the fingerprint and stranded them all.)*

### 2. Source bandwidth — 2026-07-25

The Gaussian pulse now covers the band the monitors actually sample, and the width is
a real solver parameter rather than a magic 0.5.

**Changed**

- **`beef/solvers/fdtd.py` — new `source_fwidth(wavelengths, factor)`.** `factor`
  multiplies the **full** span `1/λmin − 1/λmax`; the old code used `span * 0.5`.
  Meep's Gaussian falls as `exp(-2π²(Δf/fwidth)²)`, so measured against a 1.50–1.60 µm
  band the edges sit at 0.005 of peak at factor 1.0, 0.097 at 1.5, and **0.27 at the
  new default 2.0**. The ratio is what governs it, so this holds for any band or point
  count — confirmed identical for an O-band sweep and for a 51-point grid.
- **Single-wavelength requests no longer produce `fwidth = 0`.** Span is zero there,
  and meep's pulse width is `1/fwidth`, so the old expression yielded an unusable
  source. Falls back to `factor * 0.1 * fcen`, i.e. meep's usual `0.2 * fcen` at the
  default. Latent bug, found while here.
- **New `check_source_coverage(fcen, fwidth, wavelengths, min_fraction=0.1)`**, called
  by both `solve()` and `solve_section()`. Evaluates the real
  `mp.GaussianSource.fourier_transform` at each requested wavelength and logs a
  warning naming the worst point if any falls below 10% of peak. A starved point still
  yields a well-conditioned transmission — the source spectrum cancels in the
  measured/normalization ratio — but it burns dynamic range, and `decay_tol` is
  applied to DFT fields orders of magnitude smaller there.
- **`FDTDSolver(source_bandwidth_factor=2.0)`** is a real constructor parameter with a
  docstring entry, not a module constant — `notes.md` Solver #9's sub-note complains
  precisely about physics knobs that live in defaults nothing can see.
- **`beef/cache/config.py`** — `source_bandwidth_factor` added to `_METHOD_ATTRS` for
  both `("FDTDSolver", "solve")` and `("FDTDSolver", "solve_section")`. That list is
  an explicit allowlist, so a new attribute is *not* picked up automatically; without
  this the knob would move results while leaving the hash unchanged.

**Added** — 5 more tests in `BEEF/tests/test_fdtd_freq_grid.py`: the factor scales the
full span; a single wavelength still gets a usable pulse; the shipped default clears
the coverage floor with margin (>0.25) and peaks at fcen; the pre-fix `span * 0.5`
leaves the edges below 1e-6 **and now warns** (regression lock); and the knob is
reachable from the constructor.

**Verified** — `pytest tests/test_fdtd_freq_grid.py` 10 passed; full fast suite
`pytest tests/` 135 passed, 6 skipped, no regressions. Fingerprint check:
`source_bandwidth_factor` appears in the built config and two factors produce
different configs. Cache audit of `BEEF/examples/example_sim_cache`: 119 cached runs,
54 of them `FDTDSolver`, **0** carrying the new key — so every stale FDTD result now
misses instead of being silently reused. The 65 EME entries are unaffected and remain
valid.

**Confirmed physically** — job 11477911 (2026-07-25), check 4: the shipped default puts
the weakest requested wavelength at 0.269 of peak in a real 3-D solve, and check 9 gives
power conservation 0.9992–0.9997 with `check_passive()` clean, so the wider pulse cost
nothing in accuracy. See the check-5 caveat in Verification — the wider bandwidth was
needed *because* item 1 widened the grid, not because the pre-fix run was invalid.

**Consequence** — items 1+2 together made job 11477911 the first trustworthy FDTD output
in this project, and it is now the reference for the coupling-vs-wavelength shape at
`coupling_length` 13.5 µm. Numbers in the diagnostic sections *above* still come from
the old, mislabelled data and are kept as the record of the bug; item 5 replaces the
ground-truth length sweep.

### 3. Cache fingerprint completeness — 2026-07-25

Everything that changes an FDTD result is now in the hash, and the unreachable entries
are off disk.

**(a) Subpixel smoothing is a real parameter now.** `eps_averaging`, `subpixel_tol` and
`subpixel_maxeval` were on meep's implicit defaults — never passed, never fingerprinted
— while moving `n_eff` at the second decimal (`notes.md` Solver #9). They are
`FDTDSolver` constructor arguments, passed to **both** `mp.Simulation` calls
(`fdtd.py:501`, `:798`), and in `_METHOD_ATTRS` for `solve` and `solve_section`.
Defaults are meep's own, so no existing result changes meaning; a test asserts that,
since a wrong default here would silently alter every run.

**(c) The monitor grid does not need its own schema bump — it needs its resolved value
hashed.** New `_resolved_spectral_config()` puts `resolved_monitor_freqs` and
`resolved_source_fwidth` into the `call_kwargs` fingerprint at all three `FDTDSolver`
config sites. Both are derived from `wavelengths`, which is already hashed, so the
fingerprint was *correct* today — but that is exactly the state item 1's bug lived in
for months: the grid was derived wrongly while every cache entry looked identical to a
right one. Now a change to the derivation changes the hash even when the inputs do not.

This is the pattern `build_run_config()`'s own docstring prescribes for
resolved-vs-ingredient values (Solver recommendation #4), already used for
`resolved_layout`. Chosen over bumping `_CACHE_SCHEMA_VERSION`, which would also have
stranded the 20 valid EME entries the bug never touched, and which depends on a human
remembering to bump it.

**(b) Dead entries deleted.** 54 unreachable `FDTDSolver` folders (4.4 MB) removed from
`BEEF/examples/example_sim_cache`; 120 entries → 66, then 65 after the audit below.
"Unreachable" was tested, not assumed: an entry whose `config.json` lacks a key current
code always emits can never be matched by `lookup()`.

`index.jsonl` was also rebuilt, 5760 lines → 65. It had grown 48× because every rank
appended (item 4), and 2640 of its lines pointed at folders just deleted. Regenerated
from the surviving `config.json` files through the same `_manifest_entry()` the writer
uses, so it is now exactly one line per entry, 0 orphans.

**Added** — 15 tests in `tests/test_cache_config.py`: the three subpixel knobs
discriminating on both methods, their defaults matching meep's, and a
`TestResolvedSpectralConfig` class covering the grid being ascending in frequency and
equal to the requested wavelengths, `fwidth` tracking the factor, and — the one that
matters — **a reversed grid derivation changing the hash** while an identical call
still shares one.

**Verified** — full fast suite **175 passed, 6 skipped**. Post-deletion audit: 65
entries, **0 unreachable**, 0 orphaned index lines.

**Found while doing this, not fixed here** — the same monitor-grid bug is live in both
ring resonator example scripts, and their cache entries are still reachable. Raised as
item 10 rather than folded in silently; it is a different code path (raw meep, its own
config builder) and deserves its own change.


### 4. Cache write race — 2026-07-25

Only one rank writes the cache now, and a damaged cache entry can no longer kill the
job that reads it.

**The shape of the bug.** Every rank of an `srun -n 48` job runs the whole solve
redundantly, so every rank also reached every write path in `store.py` — 48 processes
across 2 nodes writing one directory on shared storage. `index.jsonl` is the receipt:
**5760 lines for 120 saves, exactly 48× over-counted.**

**Changed — `beef/cache/store.py`**

- **New `_is_master()`.** Deliberately returns a bool, not a rank: nothing here needs
  to know *which* rank it is, and inventing a number later callers could misuse is how
  this kind of bug starts. Sources in order of authority, importing nothing not already
  loaded (importing mpi4py into a serial process would call `MPI_Init` behind meep's
  back): meep's `am_really_master()` → mpi4py → the launcher's environment.
  `am_really_master()`, **not** `am_master()` — the latter is per-group under
  `divide_parallel_processes` and would elect one writer per group.
- **The environment fallback is load-bearing, not a nicety.** An EME-only run imports
  neither meep nor mpi4py (meow/MPB don't pull meep in) yet is still launched under
  `srun -n 48`, so a meep-only check would have called all 48 ranks master and left the
  bug fully intact on that path. `SLURM_PROCID` / the OMPI/PMI variables answer without
  importing anything.
- **`save_run()` and `update_family_pointer()` write only on master**, then `_barrier()`.
  The barrier is the important half: without it a later `lookup()` could hit on one rank
  and miss on another, and ranks disagreeing about whether to simulate leaves meep's
  collective calls unmatched — the job hangs until it burns its wall clock, which is
  what actually consumed the 5 h. Both are collective, called from all ranks or none;
  every call site is guarded by conditions derived from the solver config, never from
  per-rank data, so that holds.
- **`save_run()` still returns the folder path on every rank**, since callers feed it
  straight into `update_family_pointer()`.
- **New `_read_family_pointer()`.** A pointer is a reuse hint, never data — the result
  it names is also reachable through the plain `lookup()` path — so an empty, truncated
  or schema-drifted one is now a cache miss with a warning, exactly like `lookup()`'s
  config-mismatch guard, instead of an exception. This is the actual crash site: tasks 4
  and 5 raised out of `json.load` here and took a 48-rank job down. Used by both
  `lookup_superset()` and `update_family_pointer()`, so damaged pointers heal on the
  next write.
- **`_atomic_write_json()` temp names carry host as well as PID.** PIDs are allocated
  per host, so a PID-only name is not unique on a filesystem spanning nodes.
- **Docstring corrected.** The claim that concurrent writers cost only "a missed reuse
  opportunity, never an incorrect one" was true across separate processes (a SLURM
  array) but false for ranks within one job. Now says which case it means.

**Added** — `BEEF/tests/test_cache_multirank.py`, 25 tests: rank detection through all
four sources including the `am_master`/`am_really_master` distinction; nothing written
off master while the folder path still returns correctly; `index.jsonl` getting one line
per job rather than per rank; host-qualified temp names; and five flavours of damaged
pointer (`""`, truncated JSON, whitespace, `null`, `[]`) reading as a miss, healing on
the next save, without turning every update into a clobber.

**Verified** — full fast suite `pytest tests/` **160 passed, 6 skipped** (135 baseline +
25). Then the part unit tests cannot reach, via the new
`examples/verification/python/verify_cache_multirank.py` at **48 ranks on 2 nodes**,
production shape, meep imported, no simulation, ~1 s:

    ranks           48 on 2 hosts
    index.jsonl     5 lines (expected 5; pre-fix would be 240)
    stray tmp files 0
    pointer dims    [6, 1]

Audit of the live cache found **no lasting damage**: 0 corrupt pointers, 0 stray temp
files, all 5760 `index.jsonl` lines parseable.

**Correction to the original diagnosis.** This document previously stated as fact that
two ranks on different nodes drew the same PID. That is *plausible but unproven* — the
48-rank probe found no PID reused across hosts, and concurrent `os.replace`/read against
shared storage explains an empty file just as well. The rank gate removes every version
of the race, so the fix does not depend on which mechanism it was; the docstrings now
say so rather than asserting the PID story.

**Re-scoped, not skipped.** Item 4 originally included re-running the two dead stub
tasks. Those seven surviving `S*.npz` files are dated 2026-07-24 — **pre-fix**, since
items 1–2 landed 2026-07-25 11:56 — so re-running only tasks 4 and 5 would leave one
folder holding a silent mixture of pre- and post-fix data. A visible gap beats an
invisible mixture, so the re-run moved to item 7, which now covers all nine tasks after
deleting the folder.


### 5. Ground truth re-established — 2026-07-25

The directional coupler has a trustworthy, self-describing L₅₀(λ) across the full
1.50–1.60 µm band again. **Job 11480725, 39/39 tasks COMPLETED, 40/40 lengths on disk**
(the 40th, L=13.5, came from the validation task 11480098\_12).

**The result.**

| λ (µm) | 1.500 | 1.510 | 1.520 | 1.530 | 1.540 | 1.550 | 1.560 | 1.570 | 1.580 | 1.590 | 1.600 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L₅₀ (µm) | 20.46 | 19.31 | 18.22 | 17.20 | 16.23 | 15.31 | 14.44 | 13.61 | 12.84 | 12.09 | 11.39 |

Monotone decreasing, dL₅₀/dλ = **−90.3 µm/µm**. Coupling strengthens with wavelength,
which is the physically right direction and the sign that items 1–2 hold at the *device*
level, not just at the monitor. Every one of the 11 crossings is bracketed within
0.5 µm, so the fine window was placed correctly and no value is interpolated across a
wider gap than claimed.

**How it was built.** `examples/verification/python/verify_item_5_ground_truth.py`, run
as `--array=0-39` via `batch/verification/run_verify_item_5_ground_truth.batch`. Grid is
`union(0:4:60, 9:0.5:22)` — 40 lengths, the fine window placed from data rather than
guessed: the archived sweep re-read on its true axis gives dL₅₀/dλ ≈ −90 µm/µm, putting
the crossing between ~10.8 and ~19.8 µm across the band. Device and solver kwargs are
pinned to `directional_coupler_fdtd_sweep.py` exactly, so the result is directly
comparable to both the archive and job 11477911. Set up as a verification example rather
than by editing that script, which keeps its deliberate 0:2:60 grid and its teaching
role.

Three things it does that the tutorial sweep does not, each traceable to something this
investigation turned up:

- **Every write is rank-0 gated** with a host-qualified temp name — item 4's defect, still
  live in all five `directional_coupler/python/*.py` sweep scripts, and exactly the setup
  that killed stub-sweep tasks 4 and 5.
- **Provenance is written into the file**, both as a `provenance` key inside the npz and a
  greppable `provenance.json` sidecar: git SHA, meep version, rank count, SLURM ids, and
  the full device/solver kwargs. This is the specific defect that made the archive
  unusable as truth even though its numbers were fine.
- **L₅₀ is stored with its interpolation bracket** (`L50_bracket_lo`/`_hi`). 15.30 µm
  quoted from a 2 µm bracket and from a 0.5 µm bracket are not the same claim, and the
  archived file gave no way to tell which it was.

Caching is **on** here, unlike the items-1–2 harness — that one had to be uncached
because a verification must not pass on a stored result; this one *produces* reference
data at ~0.4–3 h per task, so a re-run after a node failure should be free. Only safe
because item 3 completed the fingerprint and item 4 made the writes rank-safe.
`--no-cache` forces re-derivation.

**Validated before committing 39 tasks — job 11480098\_12**, L=13.5 µm, 56:20 against
11477911's 56:30 on the same device. That length was chosen because 11477911 solved it
*uncached*, making it an independent cross-check rather than a cache echo. Result was
**bit-identical**: max |ΔT_cross| = 0.000e+00 across all 11 wavelengths,
`solve_seconds = 3362` proving it genuinely simulated.

**Verified — four independent checks.**

- **Power conservation:** max |1 − (T_through + T_cross)| = **1.1e-3** across all 440
  points, no NaNs, `check_passive()` clean.
- **Reproduces the verified run:** L=13.5 at 1.55 µm gives T_cross = 0.423022 against
  11477911's 0.423022, agreeing to **5e-7**.
- **Continuity with the archive** — the strongest one. Re-reading the archived pre-fix
  sweep on its true axis and overlaying all 11 of its points on the new curve: **every
  point agrees to within 0.027 µm**, and the residual is a uniform **+0.023 µm bias with
  only 0.003 µm of scatter** about it. Two runs with different source bandwidths,
  different length grids, and one of them relabelled, tracking each other to 3 nm of
  shape. This confirms the relabelling diagnosis independently of the solver's own
  convergence, which is what check 8 of the items-1–2 verification had established
  internally.

  The bias is systematic, not noise, and is not fully explained by check 8's 2.5e-4
  agreement in T — at dT_cross/dL = 0.042 /µm near the crossing that buys only 0.006 µm.
  The remaining ~0.017 µm is most likely the archive's own 2 µm interpolation bracket
  against 0.5 µm here. Not chased further: it is 0.15% of L₅₀ and the new file is the one
  being trusted going forward.

  **Correction — the figure first recorded here (rms 0.009 µm with a 0.18 µm outlier)
  was wrong.** It came from reconstructing the archive's axis with
  `fcen = 0.5*(1/λ_min + 1/λ_max)`, the frequency-mean. The pre-fix code used
  `fcen = 1/mean(λ)` (`_build_domain()` → `wl_center = np.mean(wavelengths)`), 0.7 mHz
  lower, which shifts the whole reconstructed axis by ~1.6 nm. The "outlier" was that
  offset showing up as a uniform bias, misread as scatter. The verification script's own
  `_archived_l50()` had it right all along — the error was in my ad-hoc check, so the
  plot and the shipped numbers were never affected. Worth recording because the same
  mistake would quietly degrade any future re-read of a pre-fix file.
- **The 1.55 µm anchor**, the one index both engines labelled correctly before the fix:
  new 15.31 µm against the archive's 15.29 µm.

Outputs in `example_output/verification/ground_truth_sweep/`: `length_sweep.npz`,
`provenance.json`, `length_sweep.png`, `l50_vs_wavelength.png` (which overlays the
archived true-axis points, the plot behind the continuity check above).

**Cost:** 2387 core-hours, **$16.23** for the array and **$16.84** including both
single-task validations, priced with `pace_spending_analysis/plot_usage.py`'s model
(cpu-small at $0.0068/core-hour × 48 cores). Runtime ran 0:21 at L=0 to 3:10 at L=56
against a 300 min wall — no timeouts, and no correlation tight enough to be worth
modelling: L=56 took longer than L=60.

**Also fixed while here — 48× duplicated log lines.** Job 11480098's log was 48 of 107
lines the *same* message: `logger.info("solve: cache saved at ...")` fires on every rank,
and there are 15 such cache log sites across `solvers/fdtd.py` and `solvers/eme.py`.
Those describe one shared cache, not per-rank state, and with 40 task logs to read it
buries the results. New `beef/_log.py`: off the master rank, keep WARNING and above only,
rank-tagged. Warnings survive on *every* rank deliberately — they are the records that
can genuinely be rank-specific, so suppressing them to de-duplicate INFO would trade
noise for blindness.

It has to be installed **per module**, not once on the parent `beef` logger: a logger's
filters are consulted only for records logged on that logger, never for records
propagating up from a child. The first version got this wrong — it would have filtered
nothing — and a test caught it. Six modules now take their logger from
`beef._log.get_logger()`; the root logger is untouched.

**Added** — 28 tests in `tests/test_log_rank_filter.py`, including a parametrized check
that every one of the six module loggers actually carries the filter (the test that
caught the parent-logger bug). Full fast suite **203 passed, 6 skipped**. Verified at 48
ranks on 2 nodes: INFO 48 → 1, warnings 48/48 with 47 rank-tagged.

**Found while doing this, not part of the original scope** — finishing item 5 does *not*
make the stub sweep's reference line resolve, as had been assumed. `_FDTD_SWEEP_NPZ`
(`directional_coupler_hybrid_stub_sweep.py:93`) points at the tutorial sweep's output
path, which still does not exist in the live tree; the new truth is under
`example_output/verification/ground_truth_sweep/`. Repointed — the consumer at `:352`
already reads `coupling_lengths`/`T_cross`/`T_through`, which are exactly the new file's
key names, so nothing else changed. Item 7 depends on this line drawing.


### 9. Golden reference regenerated — 2026-07-25

`tests/golden/directional_coupler_hybrid.json` is back on the fixed wavelength axis.
**Job 11486066**, serial, 1:04:31, `check_passive() = True`.

|  | 1.50 | 1.55 | 1.60 |
|---|---|---|---|
| through, old | 0.0098 | **0.0601** | 0.4289 |
| through, new | 0.0247 | **0.0601** | 0.4948 |
| cross, old | 0.9786 | **0.9220** | 0.5587 |
| cross, new | 0.9667 | **0.9220** | 0.4957 |

**The centre point is bit-identical in both ports.** With 3 wavelengths, index 1 is
`fcen` — the one point the pre-fix code labelled correctly — so it is precisely the
value that must *not* move, and it does not, to four decimals. Anything that had
perturbed the geometry, the mode count or the physics would have moved it too. That is
the single strongest piece of evidence that this is the axis fix and nothing else.

**Why the endpoints cannot be checked by relabelling, and what was done instead.** For a
pure FDTD result the pre-fix values are recoverable on their true axis — that is how the
archived length sweep was validated to 0.027 µm in item 5. It does not work here,
because this is a *hybrid* result: `np.einsum` contracts an EME body at the true λ
against FDTD sections at the shifted λ, so the old endpoints are unphysical products,
not correct values wearing the wrong label. Attempting the reconstruction anyway gives
0.39 disagreement, which is the expected outcome, not a failure. This is the same
mechanism described under "Why the stub sweep can't converge" above; worth stating
explicitly because the obvious check is the wrong one here.

Checked instead against item 5's **independent** ground truth. Fitting
`cross = sin²(π/4 · L/L₅₀(λ) + φ_bend)` with L = 30 µm and L₅₀(λ) from job 11480725:

| fit to | φ_bend | rms |
|---|---|---|
| new golden, on its own axis | +0.293 rad | **0.012** |
| old golden, labels at face value | +0.245 rad | 0.022 |
| old golden, on the FDTD sections' true axis | +2.739 rad | 0.076 |

The new values are the most consistent with the independently measured L₅₀(λ) under any
axis interpretation. One free parameter over three points, so this is supporting
evidence rather than proof — the bit-identical centre carries the argument.

**Falling `cross` is correct here** and is *not* the old backwards dispersion returning.
At L = 30 µm with L₅₀ running 20.5 → 11.4 µm across the band, the device is past its
first full crossover, so strengthening coupling pushes `sin²` back down. An early draft
of the check asserted `cross` must rise; that assumed a device below the crossover and
was wrong.

**Correction to this item's premise.** It said the pre-fix golden was uncommitted in the
working tree. It was not — it had been committed in `86622d1` (2026-07-25 14:09) as part
of the item-5 harness commit, i.e. the quiet overwrite this item warns against had
already happened, one commit earlier than anyone noticed. Nothing was lost (the history
holds both), but the regeneration was made its own commit so the move is explained where
someone will look for it.

No stale-cache risk: the generator uses an in-memory `cache: dict`, not `save_dir`.

**Also, while in there** — `batch/run_golden_directional_coupler_hybrid.batch` loads the
system `openmpi` module, which is incompatible with the MPICH-built meep in the `jmh`
env. Harmless today because the job is `--ntasks=1` and never calls `srun`, but it would
fail immediately with `PMI_Init returned -1` if this target were ever parallelized.
Commented rather than removed, since removing it changes a working script for no present
benefit.


### 6. Engine agreement measured directly — 2026-07-26

The residual gap is real, it lives in the straight coupled-pair section, and **most of
it is not engine disagreement**. Job 11489846 (4 array tasks) + job 11490488 (re-run).

> **Two later results, read alongside this.**
>
> **(a) The residual gap is ~5%, and it is not MPB-vs-meow.** The "~4.7% genuine engine
> disagreement" below was obtained by refining the EME mesh against a *fixed* res-40
> FDTD. Refining both together (job 11502074, res 20–50) shows Δβ_MPB still falling at
> res 50 while Δβ_FDTD is flat. Taking the res-80/100 plateau (Δβ_MPB = 0.0894) against
> Δβ_FDTD = 0.0852 gives **4.7–4.9%** depending on the denominator — so the number below
> happens to survive, by a different route than it was derived. *A Richardson
> extrapolation of the res 20–50 sequence gave ~7%; that is superseded and should not be
> quoted, because the sequence's convergence exponent was degrading (−1.78 → −1.15 →
> −0.61) and the extrapolation was therefore unsound.* The split is also not
> MPB-vs-meow: both engines' mode solves are MPB and agree to 0.57%; the gap is MPB's
> mode index versus the phase meep's time-stepper accumulates. Full comparison, including
> why editing `meow/fde/meep.py` can close at most 4.6% of it, under item 12.
>
> **(b) This item is closed.** Not by reconciling the engines — which (a) shows is worth
> at most 0.57% — but by propagating the EME body on the FDTD's own β. See
> **item 6 (closure), 2026-07-27**, at the end of this file. The ~5% gap itself is not
> explained and is now **item 13**.

**What was measured.** The stub sweep could only *infer* Δβ_FDTD/Δβ_EME = 0.891, from a
linear fit to L₅₀ through a whole hybrid chain with bend phase and seam mode-matching
folded in. This measures the same ratio on nothing but a uniform coupled pair — sections
cut at x ∈ [−L/2, +L/2] inside the reference device's 30 µm coupling region, so every
cut sits 15 µm clear of the nearest S-bend and both engines are being asked an identical
question.

FDTD side: `solve_section()` with 2 supermodes returns a transfer matrix whose diagonal
carries exp(iβ_even L) and exp(iβ_odd L), so Δφ(L) = arg(T₀₀) − arg(T₁₁) = Δβ·L + const.
The constant is mode-matching and monitor-plane phase — which is exactly why **one length
cannot answer this**: Δφ/L would carry it as a 1/L error. Four lengths (2, 4, 6, 8 µm),
fit the slope, discard the intercept, keep the residual as evidence the model held.

EME side: the same `section_gds()` polygons → the same meow cross-section →
`compute_modes()`, so Δβ_EME comes off the eigenvalues with no propagation involved.

**Result — it reproduces the stub sweep exactly.**

| | ratio at 1.55 µm |
|---|---|
| stub sweep, inferred through the hybrid chain | 0.891 |
| **this measurement, uniform region only** | **0.8912** |

The data is unusually clean: |diagonal| = 1.000 (lossless), max off-diagonal **6.7e-6**
(a uniform region cannot convert even↔odd, and it doesn't — so the supermode basis is
right, which is checked rather than assumed), and Δφ linear in L to ≤1e-4 rad, with
successive differences 0.1705 / 0.1705 / 0.1705. The two-mode model holds essentially
exactly, so this is not noise-limited.

**But the EME side is not converged.** Sweeping the EME mesh — which the 2026-07-22
resolution study never did, it varied only FDTD:

| EME res | 20 | 40 | 60 | 80 | 100 |
|---|---|---|---|---|---|
| Δβ_EME (rad/µm) | 0.1125 | **0.0957** | 0.0926 | 0.0894 | 0.0895 |
| FDTD/EME | 0.7576 | 0.8912 | 0.9208 | 0.9535 | 0.9528 |

Converged by res 80. The production res=40 **overstates Δβ_EME by 7%**, which accounts
for ~6 of the 11 percentage points. The genuine engine disagreement is **~4.7%** — still
larger than FDTD's own resolution scatter, so still real, but less than half the headline
number. Raised as item 11; it reaches every EME result in the repo, not just this one.

**FDTD is not the culprit.** Re-reading the archived FDTD resolution sweep at λ=1.55 (its
one correctly-labelled index) and converting T_cross back to phase, the implied Δβ moves
only **2.74%** across res 30–90, and non-monotonically — noise, not a trend. So the
earlier "mesh discretization is ruled out" claim holds for the FDTD half. meep's subpixel
smoothing is the likely reason it converges so much earlier on the same nominal grid.

**Item 8 now has direct evidence.** The ratio steps 0.8905 (1.56 µm) → 0.8730 (1.57 µm)
and stays flat after — exactly where meow's Si table switches 3.4784 → 3.4710 while the
FDTD side stays frozen at the centre-wavelength index for every wavelength. So the
material half is worth ~2% of the ratio: real, not cosmetic, but not a driver.

**A useful by-product.** `sim.get_eigenmode_coefficients()` already returns `kdom` —
each mode's dominant wavevector per band per frequency — next to the `alpha` beef uses
(`fdtd.py:869`), and beef discards it. Capturing it would give n_eff per FDTD mode from a
**single** run, replacing this four-length sweep entirely and enabling the mode-identity
check that is now item 12.

**Added** — `tests/test_engine_agreement.py`. Six fast tests with no simulation, guarding
the extraction method itself: rate recovery, immunity to the constant offset, phase
unwrapping past π, a residual that reports nonlinearity, a regression against the four
measured phases, and one that documents why a **single length would have been wrong**
(the naive Δφ/L is off 160% at L=2 and still 33% at L=8). Plus the real check,
`@pytest.mark.slow` and `xfail(strict=True)` — it encodes the requirement, not today's
behaviour, so it fails loudly *as a pass* the moment someone closes the gap.

**Verified** — fast suite **209 passed, 7 skipped**. Slow suite in job 11491895: 6
passed, 1 xfailed, 21 min. The slow test meshes FDTD at res 20 and EME at 80 on purpose:
EME must be converged or the check measures discretization, while FDTD is within 0.5% of
res 40 already, which is what makes it affordable serially at all (the res-40
measurement needed 48 MPI ranks).

**Two things left deliberately unfixed.**

Array task `11490488_1` hit TIMEOUT at 25 min. `aggregate()` runs the EME half on all 48
ranks redundantly, and that task drew a slow node. Harmless — the other tasks wrote the
same aggregate, and the data file it produced is intact. Not "fixed" by gating on master:
the EME call goes through meep, and a rank-divergent collective is precisely the hang
item 4 documented. Redundant work on 47 ranks is the cheaper failure.

`np.savez` appends `.npz` to a path that lacks it, so the host-qualified temp names lost
their `os.replace` and all four tasks died *after* simulating. `verify_item_5_ground_truth.py`
avoids this by handing `np.savez` an already-open file object; this script copied the
helper but not the calling convention. The FDTD results were already in the cache, so the
re-run cost seconds rather than 40 minutes — a direct dividend of item 3 and 4 having
made caching trustworthy. `_atomic_replace_write()`'s docstring now states the
requirement.


### 11. EME resolution split from FDTD's — 2026-07-26

Every hybrid example fed a single `RESOLUTION = 40.0` to both sub-solvers. meep is
converged there; meow/MPB on its own `Mesh2D` is not.

| EME res | 20 | 40 | 60 | 80 | 100 |
|---|---|---|---|---|---|
| Δβ_EME at 1.55 µm (rad/µm) | 0.1125 | **0.0957** | 0.0926 | 0.0894 | 0.0895 |

res 40 **overstates the even/odd splitting by 7%**, and 80 is the first converged
point. That was ~6 of the 11 percentage points of item 6's FDTD/EME gap, and it is
present in every EME result in the repo.

The three hybrid scripts now declare `FDTD_RESOLUTION = 40` and `EME_RESOLUTION = 80`
separately. **`EMESolver`'s own default is deliberately unchanged**: raising it would
re-invalidate all four goldens (including item 9's, regenerated the day before) and
strand the 20 valid EME cache entries — for something that should be a per-device
choice made with a convergence check in hand, not a library-wide default.

**Why this hid for months**, both reasons worth keeping because both will recur:

- `directional_coupler_hybrid_resolution_sweep.py` meshes **both** engines at the same
  `R` (`:170`, `:183`), so it measures the two convergences convolved. Its curve looks
  flat because FDTD dominates the observable it plots (`T_cross`) while EME dominates
  Δβ. Now commented in place: a flat *joint* curve is not evidence that either engine
  converged.
- Δβ is a **difference** of two n_eff near 2.3, so a 0.1% error in each mode is ~10% in
  their difference. The quantity the coupler depends on is far more mesh-sensitive than
  the quantity anyone looks at.

**Added** — `tests/test_eme_convergence.py`: two slow tests asserting res 80 is
converged (against 100) and that the old default was *not* (so it reports if a
meow/MPB change ever closes it), plus nine fast tests with no simulation checking the
scripts still declare separate knobs, keep EME ≥ 80, and never reintroduce a shared
`res=RESOLUTION`. The physics tests are `--runslow` only, so on a normal run those
nine are the only guard against an edit silently restoring the 7%.

**Verified** — `11 passed` in 16:19 (9 fast + both slow).


### 12. Mode identity across the FDTD/EME join — 2026-07-26

`HybridSolver` contracts index *i* of an FDTD interface against index *i* of an EME
body, assuming they are the same mode. Nothing tested that. The two engines solve MPB
on **different grids** — meep's Yee mesh vs meow's `Mesh2D` — so a reorder was possible
and would have mixed modes silently.

**They agree.** At 1.55 µm on the coupled-pair cross-section, res 40 both sides:

| band | meep | meow | diff |
|---|---|---|---|
| 0 | 2.324450 | 2.325116 | −0.0007 |
| 1 | 2.300984 | 2.301519 | −0.0005 |
| 2 | 1.744418 | 1.746722 | −0.0023 |
| 3 | 1.631969 | 1.634080 | −0.0021 |
| 4 | 1.466926 | 1.467593 | −0.0007 |
| 5 | 1.436725 | 1.430456 | **+0.0063** |

Ordering preserved, correspondence is the identity permutation, worst guided-band
disagreement 0.0023. **The contraction is sound for this device.**

**New: the interface basis contains a mode that is not guided.** The cladding index is
1.444. Band 5 sits at **1.4367 — below it**, and band 4 at 1.4669 is marginal. The
coupler carries `N_MODES = 6` at every interface, so at least one column is a radiation
mode, whose n_eff and ordering are not well defined — and note its engine disagreement
(+0.0063) is an order of magnitude worse than any guided band, which is exactly what an
unbound mode living on the discretization should look like. The tests are strict only
about guided bands, and assert the unguided one exists so the finding cannot rot.

This bears on `N_MODES` everywhere, and the coupler is the *easiest* case: `mmi_1x2`
carries 12, `mmi_2x2` 15, `modal_multiplexer` 16 across four interfaces. None has been
checked. Recorded in each device's validation doc.

**Added — `FDTDSolver.section_mode_neff()`**, the FDTD-side counterpart of
`EMESolver.compute_modes()`. Uses `sim.get_eigenmode()`, which runs MPB on the
initialized grid and returns with **no time-stepping**, so a mode-identity check costs a
`set_epsilon` plus one eigensolve per mode — not the minutes a `solve_section()` needs.
That is what makes checking this routinely affordable. n_eff is `k/f`: meep gives the
Bloch wavevector in 2π/a units and the frequency in c/a, so their ratio is the index
with no 2π anywhere — pinned on a bare 450×220 nm strip at 2.3124, where the competing
normalizations give 0.37 and 14.5.

**A trap worth recording — never compare across resolutions.** The first version of this
check compared meep@40 against meow@**80** and reported a reordering. There was none:
the res 40→80 shift is ~0.027, *larger* than the 0.022 even/odd splitting, so
nearest-neighbour matching pairs band 0 with band 1 purely from the offset. Only
same-resolution comparison is valid. `assert_modes_correspond()` takes both arrays from
one resolution for that reason, and a test pins the reason so nobody "fixes" the checker
to accept the cross-resolution case.

The same trap bit the test suite once more: a *uniform* offset larger than the tightest
mode spacing is indistinguishable from a reorder by nearest-neighbour matching, so the
tolerance test had to perturb a well-separated band to isolate the magnitude check. Both
failures are the same fact from different directions, and both are now tests.

**Added** — `tests/test_mode_identity.py`: five fast tests with no simulation guarding
the checker itself (it detects an adjacent swap, a missing mode, an over-tolerance
disagreement, and rejects the cross-resolution comparison), plus two slow ones running
both engines for real.

**Verified** — full fast suite **223 passed, 11 skipped**.

#### How the two mode solves actually differ — and why it is *not* item 6's cause

Recorded because it is the obvious first suspect and it has now been checked directly.
"meep's FDTD modes vs meow's FDE modes" is a false framing: **both are MPB**.
`meow/fde/meep.py:114` calls `sim.get_eigenmode()` and so does
`FDTDSolver.section_mode_neff()` (`fdtd.py:823`). They differ only in setup:

| | `compute_modes_meep` | `section_mode_neff` |
|---|---|---|
| propagation axis | **z** — `direction=NO_DIRECTION`, `kpoint=(0,0,1)`, n_eff from `k[2]·λ` | **x** — `direction=mp.X`, n_eff from `k.x/f` |
| cell | `(x_span, y_span, 1)`, **no PML** | `(Sx, Sy, cell_sz)` with `mp.PML(dpml)` |
| grid | `1/dx` from meow's `Mesh2D` | `self.res`, meep's own Yee mesh |
| ε model | `meep_material` — nearest table λ, **non-dispersive** | `material_stack.si_medium()`, dispersive |
| eigensolver tol | explicitly `1e-12` | meep's default (not set) |
| subpixel | `eps_averaging=True`, library defaults | `eps_averaging` / `subpixel_tol` / `subpixel_maxeval` configurable |
| fields | sampled on `mesh.x_/y_`, then `np.real()` — imaginary part discarded | not sampled; only `k` is used |

Every one of those is a genuine difference, and together they are worth almost nothing
here. The table above puts the worst guided-band n_eff disagreement at **0.0023**; on Δβ
at res 40, meow gives 0.09566 and meep's own MPB 0.09512 — **0.57% apart**.

**So this is the ceiling on editing `meow/fde/meep.py`.** Making meow's solve bit-identical
to meep's would move the FDTD/MPB Δβ gap from 10.88% to 10.38% — it closes **4.6%** of it
and leaves 95% standing. The route is not worth taking for item 6's sake. (The defects in
that file — the `np.real()` truncation on field profiles, the non-dispersive
nearest-λ material, the missing PML — are still real and still worth fixing. They land on
*interface overlap integrals* and on item 8's material half, not on Δβ.)

**The difference that does matter is not between the two mode solves.** It is between
MPB's mode index and the phase meep's time-stepper actually accumulates — an eigenproblem
in a planewave basis versus a leapfrog on a staggered Yee grid, two discretizations of the
same continuum operator with different numerical dispersion:

| | β even (rad/µm) | β odd | Δβ |
|---|---|---|---|
| MPB — either engine | 9.42255 | 9.32743 | 0.09512 |
| meep FDTD — fitted from `T(L)` over 4 lengths | 9.63292 | 9.54767 | 0.08525 |
| excess | **+2.233%** | **+2.361%** | **−10.38%** |

Both supermodes are biased the *same direction* by ~2.3%; the differential between them is
only **0.129 percentage points**. But `β/Δβ ≈ 99`, and that lever turns it into the whole
10.4% Δβ error — `Δβ_MPB + β₀e₀ − β₁e₁` reproduces the measured Δβ_FDTD exactly. Matching
meow to meep's MPB cannot help, because meep's MPB carries the same 2.3% bias.

**Nor does resolution.** Sweeping res 20/30/40/50 on the same cross-section (job 11502074),
the two sides move in opposite ways:

| res | β₀ excess | β₁ excess | differential | Δβ_MPB | Δβ_FDTD | gap |
|---|---|---|---|---|---|---|
| 20 | +5.850% | +6.234% | 0.3848% | 0.11252 | 0.08452 | 24.9% |
| 30 | +3.195% | +3.381% | 0.1867% | 0.09926 | 0.08519 | 14.2% |
| 40 | +2.203% | +2.338% | 0.1342% | 0.09566 | 0.08525 | 10.9% |
| 50 | +1.862% | +1.979% | 0.1171% | 0.09221 | 0.08299 | 10.0% |

**Δβ_FDTD is essentially resolution-independent** — 2.7% across the whole range and
non-monotone, i.e. scatter. Δβ_MPB moves 18% and is still falling at res 50, with
decelerating exponents (−1.78 / −1.15 / −0.61), so it is not on a clean Δx² path either.
**A gap survives in the limit, and the honest number is ~5%.** Item 6's own EME sweep runs
to res 100 and plateaus — 0.1125 / 0.0957 / 0.0926 / 0.0894 / 0.0895 at res 20/40/60/80/100
— so Δβ_MPB converges to **0.0894** against Δβ_FDTD's **0.0852**: a **4.7–4.9%**
disagreement depending on which you divide by. Refining the mesh does not remove it.

*(Richardson extrapolation of the res 20–50 sequence alone suggested ~7%, with a
relative-gap fit at 6.6%. Both are superseded by the plateau above and should not be
quoted: that sequence's exponents were **decelerating** — −1.78 → −1.15 → −0.61 — so it
was not on a power law and the extrapolation had no basis. Recorded because the mistake is
easy to repeat: a monotone sequence is not a converging one, and Richardson on a
sequence whose exponent is drifting will confidently return a wrong limit.)*

That is what motivates `HybridSolver(calibrate_body_neff=True)`: not matching the two mode
solvers to each other, but propagating the EME body on the β the FDTD reference actually
has.

### 10. Ring resonator monitor grid — 2026-07-26

The second place the item 1 defect appeared. `ring_resonator_hybrid.py` and
`ring_resonator_hybrid_cascaded.py` each hand-built their meep simulations — the ring's
coupling ports are y-normal and nothing in `beef/solvers/fdtd.py` draws those — so each
carried its own copy of

```python
fwidth = (1/WL_MIN - 1/WL_MAX) * 0.5
sim.add_mode_monitor(fcen, fwidth, N_WL_FDTD, ...)
```

in raw meep, where the solver-side fix could not reach it. A request for 1.50–1.60 µm
sampled 1.5254–1.5754 µm: half the span, reversed.

**Changed** — both scripts now take an explicit ascending frequency list from
`monitor_freq_grid()`, size the pulse with `source_fwidth(wavelengths, 2.0)` (exactly 4×
the old width), reorder every measured amplitude onto the caller's axis with
`[wl_order]` — ratio taken in monitor order first, reordered once, as `solve_section()`
does — and call `check_source_coverage()`. `_sim_params()` now takes `wavelengths` and
emits `resolved_monitor_freqs` / `resolved_source_fwidth`, mirroring
`_resolved_spectral_config()`.

**The stranding was verified, not assumed.** All 6 pre-fix FDTD cache entries now miss
(`lookup()` → `None`); each folder's stored `config.json` was re-hashed first to confirm
it matched its own folder name, so the check ran against the real recorded config rather
than a reconstruction of it. The 39 `solve_neff` / `solve_body_supermodes` entries are
untouched and all still hit — they are MPB solves the monitor grid never affected, and
stranding them would have forced pointless re-simulation. Per-entry hashes are in
`ring_resonator_validation.md`.

**Added** — `tests/test_ring_monitor_grid.py`, 14 fast tests, no meep. Source-level,
because importing either script *runs* it. Docstrings and comments are stripped before
matching, so the prose describing the old line cannot pass for live code.

**Verified** — full fast suite **237 passed, 11 skipped** (was 223).

**Not done — nothing has been re-simulated.** The fix guarantees a re-run cannot
silently serve pre-fix data; it does not produce post-fix data. Every ring plot in
`example_output/` is still pre-fix. The re-run and its centre-wavelength check are
tracked in `ring_resonator_validation.md`.

**`notes.md` is worth checking before assuming these were the last two** — this was the
second place the same line appeared, which is the argument for the library-path
recommendation in `ring_resonator_validation.md` rather than a third hand-port.

### 7. Stub sweep re-run — 2026-07-26

> **Superseded 2026-07-27 — the sweep does converge.** Read this section as the record of
> the **uncalibrated** method; its measurements and its reasoning both stand. Re-running
> the identical sweep with `calibrate_body_neff=True` flattens the drift from +0.1085 to
> −0.0007 µm/µm. The conclusion below — that no amount of stub tuning can fix an engine
> disagreement — was right; what changed is that the engine disagreement itself was
> removed. See item 6 (closure) at the end of this file.

**Job 11500412**, all nine array tasks, 2 nodes × 48 ranks each, ~1h40m per task. The
folder was cleaned first and the seven pre-fix `S*.npz` archived to
`archive/example_output_20260724/`, so nothing on disk is a mixture of pre- and
post-fix data.

**The sweep still does not converge. That is the result, and it is worth more than a
convergence would have been.**

| BODY_STUB (µm) | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 | 1.25 | 1.50 | 1.75 | 2.00 |
|---|---|---|---|---|---|---|---|---|---|
| L₅₀ (µm) | 14.318 | 14.350 | 14.377 | 14.402 | 14.428 | 14.455 | 14.483 | 14.511 | 14.538 |

Fit: **L₅₀ = 0.1085·s + 14.3205 µm**, max residual to that line **0.0024 µm** across the
whole 0→2 µm range. The pre-fix sweep's residual was 0.0022 µm. **The curve is exactly
as straight as it was before.** No curvature, no plateau, no sign of settling.

**What did change: the slope more than halved, and its size is now fully explained.**
From `L₅₀(s) = [π/2 − 2φ_bend]/Δβ_EME + 2s·(1 − Δβ_FDTD/Δβ_EME)`, the slope carries a
factor of 2 for the two FDTD sections:

| | slope | 1 − ratio | Δβ_FDTD/Δβ_EME |
|---|---|---|---|
| pre-fix sweep | 0.2427 | 0.1213 | 0.8787 |
| item 6, measured directly (EME res 40) | 0.2176 | 0.1088 | 0.8912 |
| **predicted** = item 6 × item 11's 7% EME correction | 0.0928 | 0.0464 | **0.9536** |
| **post-fix sweep, measured** | **0.1085** | 0.0543 | **0.9457** |

Predicted 0.9536 against measured 0.9457 — **agreement to 0.82%**. Two completely
independent measurements of the same quantity: item 6 solved four section lengths and
extracted Δβ from the phase slope; item 7 reads it off the drift of a 200-length,
11-wavelength hybrid sweep. They agree.

**So the stub sweep was measuring the engine disagreement all along.** Its drift is not
a discretization artifact, not a stub-length artifact, and not something more stub will
fix — it is `2·(1 − Δβ_FDTD/Δβ_EME)`, and a straight line is precisely what that formula
predicts when the ratio is a constant ≠ 1. Items 1–2 and 11 removed 55% of it by fixing
the wavelength axis and the EME mesh. The remaining 5.4% is item 6's open, genuine
engine gap.

**How far off it still is.** At the largest stub tried, L₅₀ = 14.538 µm against the
full-device reference of **15.13 µm** (item 5's band-mean crossing — not the
per-wavelength L₅₀(1.55) = 15.31 µm; the mean of the ratios crosses lower than the mean
of the crossings because L₅₀ falls with λ). Still 0.59 µm short. Extrapolating the fit,
the sweep would reach the reference at a **7.5 µm** stub — a quarter of the entire 30 µm
coupling region, which is not a split between two engines in any meaningful sense.

**What this closes and what it opens.**

- The question this document was opened to answer — *why does the stub sweep never
  converge?* — is **answered**. Two causes, both now quantified: a mislabelled
  wavelength axis (items 1–2) and an unconverged EME mesh (item 11) accounted for 55%,
  and the residual is a real disagreement between meep and meow about how fast phase
  accumulates through an identical uniform cross-section.
- **`notes.md` Solver #10(b)** — promoting the stub sweep into a reusable convergence
  procedure — should **not** be acted on. The sweep cannot converge while the engines
  disagree, so as a general-purpose tool it would report a failure that has nothing to
  do with the device under test. It is a good *diagnostic* for the engine gap and a bad
  *procedure* for choosing a cut point.
- **Item 6 is now the only thing standing between this device and a trustworthy hybrid
  result.** It was ranked below item 7; that ordering is now inverted by measurement.

**Verification harness** — items 1–2's `verify_items_1_2` was kept open specifically to
re-check this re-run and has now done its job. It, the batch file, and the Makefile's
`verify_items_1_2` / `clean_verification` targets were to be deleted once item 6 closed.
**Item 6 closed 2026-07-27, so that deletion is now unblocked** — not done here, and it
should be a deliberate commit of its own rather than a side effect of the calibration
work.


### 6. Engine agreement — closed by β calibration — 2026-07-27

Job 11517245. The decisive test, registered in advance with two falsifiable predictions,
and both held.

**What was changed.** One thing, and only one: the EME body propagates on the β the FDTD
sections actually accumulate rather than on meow's mode index. Device, resolutions, mode
count, stub values and interface monitors are all identical to job 11500412, so the two
sweeps differ in exactly the quantity under test.

| | drift slope (µm/µm) | L₅₀ spread | rms residual | mean L₅₀ | vs 15.13 µm reference |
|---|---|---|---|---|---|
| uncalibrated | +0.108543 | 0.220 µm | 0.0014 | 14.429 | **−5.37%** (at stub 0) |
| **calibrated** | **−0.000651** | **0.005 µm** | 0.0015 | **15.1468** | **+0.111%** |

**The drift is gone.** 167× smaller, and −0.00065 µm/µm over the 2 µm stub range is
0.0013 µm — below the run's own 0.0015 µm rms scatter, so it is not distinguishable from
zero. The prediction was that the slope collapses *identically*, because it equals
`2s·(1 − Δβ_FDTD/Δβ_EME)` and the change sets that ratio to 1. It did.

**And it moved toward the independent reference, not away.** L₅₀ = 15.147 against the
full-device FDTD's 15.13 — a **0.11%** disagreement where the uncalibrated sweep sat at
5.37%. The forecast was 15.02 (0.74%); the result beat it.

**The calibration cross-validates, which is the part that was not guaranteed.** At
1.55 µm the run logged:

| | n_eff even | n_eff odd | Δn_eff | Δβ (rad/µm) |
|---|---|---|---|---|
| MPB | 2.324450 | 2.300984 | 0.023465 | 0.095119 |
| FDTD, from its own phase | 2.375691 | 2.354670 | 0.021021 | 0.085212 |
| excess | +2.204% | +2.333% | | ratio **0.8958** |

Done item 6 above measured 0.08525 and a ratio of 0.8962 by a completely different route:
lengths 2/4/6/8 µm instead of 1.26/2.52/3.78, two modes instead of six, and an ad-hoc
refit instead of `section_mode_beta`'s seeded branch selection. They agree to **0.05% on
Δβ and 0.026% on the individual β**.

That is the first real check on the branch selection, and it was the one loose end in the
implementation. A wrong-branch pick is invisible to every guard the code has — it produces
a clean straight fit with a residual around 1e-14 and a β that is 33% wrong. No residual
test can catch it. Reproducing a known answer from *different* length spacings is the only
evidence available that the seeding lands on the right branch, and it does.

**The 0.128-percentage-point differential reproduces too** (0.129 previously), which is
the quantity `β/Δβ ≈ 99` levers into the whole 10% Δβ error.

**What this does and does not establish.** It establishes that the hybrid decomposition is
now *self-consistent*: the answer no longer depends on where the chain is cut, which is
the property that makes a hybrid method meaningful at all. It does **not** establish that
15.147 µm is physically correct — that rests on the full-device FDTD being right, and the
MPB/FDTD disagreement behind item 6 is still unresolved. Opened as **item 13**. Quote the
0.11% with that caveat attached.

**A concurrency bug, found the hard way.** 7 of the 9 array tasks died with
`zipfile.BadZipFile` from `np.load` inside `aggregate()`. Each task writes its own
per-stub `.npz` and then reads *every* task's file; the per-stub write used a plain
`np.savez`, which is not atomic, so a task could read a peer's half-written file. The
script already had `_atomic_replace_write` and already used it for the aggregate
outputs — the one file that other tasks read was the one left out, and the helper's own
docstring stated the hazard.

The result survived only by luck: two tasks happened to aggregate after every write had
finished, and `stub_sweep.npz` is intact. The identical code passed in job 11500412 purely
because those tasks finished 12 minutes apart instead of 6. A latent race, not a new one.

**Fixed by consolidation rather than a local patch.** `beef._utils.atomic_replace_write`
is now the single implementation, exported from `beef`, and `beef/cache/store.py`'s
`_atomic_write_json` delegates to it. Seven example scripts had their own copies, all
using a **PID-only** temp name — weaker than the cache's, which carries the hostname
because PIDs are per-host and these jobs span 2 nodes (item 4). The same unguarded
per-task write existed in four more array sweeps (`hybrid_modes_sweep`,
`hybrid_resolution_sweep`, `fdtd_sweep`, `fdtd_resolution_sweep`); it had simply never
fired. All five now write atomically.

The helper also **raises on the `np.savez` suffix trap** rather than failing obscurely:
given a path with no `.npz`, `np.savez` appends one, so handing it the temp name writes
`<tmp>.npz`, leaves `<tmp>` missing, and the rename fails *after* the simulation work is
done. That is how four tasks of an earlier job died (Done item 6 above). Callers must pass
an open file object; if they don't, the error now names the cause.

**A second trap, one layer down, caught in review of this very fix.** The obvious way to
satisfy the above is `lambda tmp: np.savez(open(tmp, "wb"), ...)`. That is wrong in a way
that would have quietly reintroduced the original bug: **`zipfile` does not close a file
object it was handed** (`ZipFile.close()` skips `self.fp` when `_filePassed` is set), so
nothing flushes the handle before `atomic_replace_write` calls `os.replace`. It works on
CPython only because refcounting finalizes the orphaned file object the instant `np.savez`
returns. Anything that delays that — a lingering traceback holding the frame, a
non-refcounted runtime — and `os.replace` moves a **truncated** file into place, which is
precisely the failure being fixed, now with the atomicity machinery lending it false
credibility. All five sweeps use a named writer with `with open(tmp, "wb") as fh:`
instead, and the helper's docstring shows that form rather than the lambda.

**Added** — `tests/test_atomic_write.py`, nine fast tests with no simulation: completeness
and no temp residue, atomic overwrite of an existing target, the `.npz` trap rejected with
a diagnosis and the target left untouched, the hostname requirement pinned so nobody
reverts it, and a threaded reader polling `np.load` while a writer runs. That last one
carries its own control — the same probe against a plain `np.savez`, which *does* catch
`BadZipFile`, so the guard test cannot silently become vacuous.

**Verified** — fast suite **271 passed, 1 skipped** (up from 262). One pre-existing test,
`test_cache_multirank.py::test_temp_name_is_host_qualified`, monkeypatched `store.socket`
and was repointed at `_utils`; it still guards the cache's own call site, and the new file
guards the helper directly. The four goldens are untouched — `calibrate_body_neff`
defaults off and no solver default changed.

**And verified where it actually matters — job 11522231, a real array re-run: 9/9
COMPLETED**, against 2/9 for the same sweep before the fix, with zero `BadZipFile` in any
of the nine logs. Unit tests cannot prove a concurrency fix; a job that previously failed
7 times out of 9 and now doesn't is the evidence. `L_50` is unchanged to 5e-6 (the
precision the pre-fix values were recorded at) and the drift slope reads −0.000648 against
−0.000651 — the change is to the write path only, and the physics did not move.

**One cost worth recording: cached is not the same as cheap.** The re-run wrote **zero**
new cache entries — every FDTD section and every EME body solve hit — and still took 34
minutes per task, against 1h51m uncached. The remaining time is not simulation: all 48
ranks redundantly unpickle the same ~55 MB of cache and redo the `[200 × 11 × 6 × 6]`
einsum chain and the plotting. Anyone estimating the cost of a re-run from "it's cached"
will be wrong by about 30 minutes. Gating the aggregate and the deserialization on the
master rank would remove most of it, and is the same rank-gating pattern item 4 applied to
cache writes — not done here, and not urgent.
