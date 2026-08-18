# 2×2 MMI validation — to-do

> **★ 2026-08-17 — this device is the reference case for the new β-calibration
> cost model, and its calibration was re-run on the merged multi-monitor path.**
> Two library changes landed and both apply automatically to every device
> (`mmi_2x2_hybrid.py` itself is unchanged): `measure_diagonal_only=True` stops
> each run decomposing 14 bands it discards 13 of (~4x, exact — verified
> bit-identical against this device's own cached matrices), and
> `solve_section_multi_length()` gained a cache key, making 14 sims × 3 monitors
> the default instead of 42 sims × 1 monitor (~1.9x more). This device's
> calibration went 30.6 h → **4.14 h measured** (job 11937894, COMPLETED, 14
> sims and 770 eigensolves exactly as predicted); the next run should reach
> ~3.4 h now that the seed is cached and anchored. The device answer moved
> max 5.5e-03 in T_bar / 3.6e-03 in T_cross across the 201×11 grid — step 6's
> domain-length effect at production scale, ~0.1% of the split ratio and an
> order of magnitude inside this device's own gap-5 discrepancy. Its 42 per-length entries in `example_sim_cache/` remain valid and are
> NOT reused by the merged path, deliberately. A third change landed the same
> day: the MPB seed is now cached and solved at `seed_anchors=3` wavelengths
> instead of all 11 (154 eigensolves → 42, ~31.5 min → ~8.6 min, free on a
> re-run) — verified against this device's real 11-wavelength seed by job
> 11940556. **Read step 9 of
> `section_mode_beta_multi_monitor_validation.md` before re-running any other
> device's calibration** — it is self-contained and carries the cost model to
> budget with, since a calibration is ~97.6% of a hybrid example's wall time.

> **2026-07-28 — the material model changed under every device here.** Both engines
> now evaluate one analytic Sellmeier model (Salzberg Si, Malitson SiO2) instead of
> each doing its own nearest-tabulated lookup, which shifts n_Si by ~0.002 and
> removes a 0.0063 step that sat in the middle of the C-band. Any golden or stored
> result for this device predates that and needs regenerating; see
> `directional_coupler_validation.md` item 8 and, for the dispersive-FDTD
> re-parameterisation (built and gate-passed 2026-07-31, opt-in via
> `dispersive=True`, off by default), `material_model_validation.md`.



Opened 2026-07-26. Nothing is known to be wrong; this is the gap list against the
four checks in [README.md](../README.md).

**Status 2026-07-26:** check 2 (mode identity) is **done and passed** — job
11500536, see below. It was expected to be the riskiest device in the repo and
turned out to be the cleanest. Checks 1, 3 and 4 remain unrun, and there are
still no goldens.

**Status 2026-07-30:** two config changes and one new doubt. `N_MODES` is now
**14**, not 15 (gap 3, done — the unguided column is dropped). Gap 5 is in
flight as jobs 11585844/11585849 (re-queued 2026-07-31 under
`gts-sr108-gedc_main` — the original 07-30 submission never started, blocked
on `AssocGrpBillingMinutes`). And gap 6 is **new**: the passing check-2
result above was measured before either duplicate-mode guard existed, on a
mode count this device no longer requests, so it needs a cheap re-run before it
can be quoted as current. Checks 1 and 3 remain unrun; still no goldens.

**Status 2026-07-31:** three of the four remaining gaps closed today. Gap 2
(check 3, resolution) closes at **res 100** on Δβ_01 (job 11590785). Gap 5
(check 4) measured: 0.82%/0.62% max\|ΔT\|, **sign flips across the band**
(not the 1×2's uniform under-prediction) — `calibrate_body_neff` not
obviously indicated the same way here, untried regardless. Gap 6 (guard
reconfirm) closed (job 11588569): "cleanest device audited" is current
again, neither guard fired. **Only check 1 (gap 4, engine agreement) and
gap 1 (goldens — never existed) remain.**

**Update, same day, 19:46 EDT: gap 4's attempt (job 11591249) TIMED OUT**
at its full 5h budget rather than completing — `sacct` confirms
`State=TIMEOUT`, `ExitCode=0:0`, not a crash. All 4 lengths' first FDTD run
finished, but a 5th run was only ~5,300 of a likely ~30,000 timesteps in when
killed — roughly halfway through the job's total FDTD work, per-length. **Gap
4 remains open**; this device's engine agreement is still unmeasured, the
same status as the 1×2. See gap 4 below and `running_jobs.md`'s own 11591249
entry for the full timing breakdown. Goldens (gap 1) are still blocked behind
gap 4 plus the coupler's own cold regeneration pass
(`directional_coupler_validation.md` Verification → Outstanding),
launched today as 11592462/11592463 (11592463 completed 17:42; 11592462
still running as of this update) — see `running_jobs.md`.

**Status 2026-08-01 03:07 EDT — ✅ GAP 4 IS CLOSED. Every check on this
device has now been run.** The reshape held: jobs 11596845 `[0-3]` (array,
one length each) all COMPLETED, and merge job 11596851 assembled the answer
from cache in 25:29. **Mean ratio Δβ_FDTD/Δβ_EME = 0.9953** — the two engines
agree on this body's phase accumulation to **0.47%**, against the coupler's
~10% on the same measurement. So the coupler's discretization mechanism is
**not** present here at a comparable size and `calibrate_body_neff` is not
indicated for this device — a real negative result, and the reason it is
worth having: gap 5's 0.82% now has no known calibratable component behind
it. Numbers, diagonality check and the argument against gap 5 in gap 4 below.
**Only gap 1 (goldens — never existed) remains open on this device**, and it
is blocked on the repo-wide material-model regeneration, not on anything
2×2-specific.

**Status 2026-08-01 (evening) — gap 5 re-measured dispersive, and P4's
candidate mechanism is REFUTED** (jobs 11609211_1 + 11609680 + 11609681).
Gap 5 is now **0.76% / 0.59%** max|ΔT| (frozen-era: 0.82% / 0.62%), and **the
sign flip survives in the same grid interval** — the material model was not
producing it. The FDTD reference moved by at most 0.00033 under the flip,
~20× less than the 1×2's 0.0065, so the scale argument that had looked like
support was measuring the 1×2's much larger material sensitivity rather than
this device's. **P4 stays open and structural — tapers or the cascade.** Full
analysis in the P4 box under *Parity with the coupler*.

Same evening, the 1×2's check 1 landed for the first time at **0.9946**
against this device's **0.9953** — 0.07 points apart, different geometry,
different band, same direction. **The 0.47% Δβ bias is an engine-level
property, not a 2×2 artifact**, which strengthens gap 4's negative result
about `calibrate_body_neff` rather than weakening it.

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

   > **🟡 Half done 2026-08-01 — the runners now exist.**
   > `run_mmi_2x2_eme` / `run_mmi_2x2_fdtd` / `run_mmi_2x2_hybrid` are written
   > and registered in `TARGETS` (4 golden targets → 7). This gap was never
   > merely "not generated": there was no runner to generate from, and that
   > half is code, not compute, so it did **not** have to wait on the cold
   > regeneration pass. Doing it now means that pass covers all three devices
   > in one round rather than needing a second round for this one.
   >
   > Each replica mirrors this device's own example: C-band 1.530–1.565 (not
   > the 1×2's O-band), top input excited (`input_idx=1`), bar = `output_idx=1`,
   > cross = `output_idx=0`, `half_domain_aperture` on **both** end sections
   > since this device has two ports 1.5 µm apart at each end.
   >
   > Two deliberate departures, both recorded in the code:
   > - **`n_modes = 6` against production's 14.** Chosen conservatively, not
   >   scaled: since 2026-07-31 `_iter_distinct_bands()` *raises* when a
   >   cross-section cannot supply the requested guided modes rather than
   >   silently padding, so an over-request at the reduced res 20 would fail
   >   generation outright. 6 is safe by comparison — `run_mmi_1x2_hybrid`
   >   already runs 6 at res 20 on a narrower, thinner core (2.0 µm × 160 nm vs
   >   2.8 × 220), which supports strictly fewer guided modes.
   > - **The hybrid replica does not use the example's mirror-symmetry
   >   shortcut.** `mmi_2x2_hybrid.py` solves the input section once per arm and
   >   transposes; that is a real, verified saving but it is the *example's*
   >   assembly. A golden pins `HybridSolver.solve()` on its defaults — the same
   >   reasoning the coupler's golden records for `calibrate_body_neff` — so
   >   every section is actually simulated here. Costs more and gives slightly
   >   different absolute numbers **on purpose**: the example is the physical
   >   answer, this is the regression pin.
   >
   > Still blocked: **generating** the three files, on the coupler's cold
   > regeneration pass (`README.md` parity item P7).

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

   **✅ MEASURED 2026-07-31 — job 11590785 (1:05:48).** This device's own
   self-imaging length (`mm_length=28.0`) is set by the beat length between
   its two lowest guided modes, Δβ_01 = 2π/λ·(n_eff[0]−n_eff[1]) — the same
   logic as the coupler's Δβ, applied to this device's own relevant mode pair
   rather than assumed to transfer. `verify_mmi_2x2_gap2_dbeta.py` swept res
   40/60/80/100/120:

   | res | Δβ_01 (rad/µm) | step | relative step |
   |---|---|---|---|
   | 40 | 0.16481 | — | — |
   | 60 | 0.16389 | −0.00092 | 0.56% |
   | 80 | 0.16340 | −0.00049 | 0.30% |
   | 100 | 0.16319 | −0.00020 | 0.12% |
   | 120 | 0.16311 | −0.00008 | 0.049% |

   **Correction to this entry's own earlier partial reading.** The 40→60/
   60→80 readout in flight compared Δβ_01's raw step size directly against
   n_eff[0]'s own raw step (0.0121, 0.0068 — a *different quantity*, ~2.8 in
   absolute value against Δβ_01's ~0.16) and called it "an order of magnitude
   tighter." That comparison mixed units and was wrong. Read relative to each
   quantity's own value, the two converge at a **similar rate** — n_eff[0]'s
   steps were 0.43%/0.24%/0.096% (40→60/60→80/80→100) against Δβ_01's
   0.56%/0.30%/0.12% over the same three intervals — not the dramatic
   divergence the coupler's T_cross-vs-Δβ story showed.

   **What the full sweep actually establishes: Δβ_01 is converged by res
   100–120 in every practical sense**, regardless of how its convergence
   *rate* compares to n_eff[0]'s. The step drops to 0.049% (0.05 in a
   thousand) by res 120, and `L_π = π/Δβ_01` — the quantity that actually
   sets the self-imaging length — inherits that same tiny relative error.
   That is the answer this item asked for: the physically load-bearing
   quantity plateaus at a resolution well within reach, so gap 2 is closed
   for this device at **res 100** (the last point with a >0.1% step;
   generous headroom over res 120's 0.05%), not blocked on an unconverged
   absolute n_eff the physics was never sensitive to in the first place.

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
   applied as part of `resolved/fdtd_vs_EME_mode_validation.md` item 9 (the cross-device
   pass that drops each device's one unguided column). The configuration
   change this gap asked for has therefore already happened, which moves the
   goldens question in gap 1 / `## Order` forward one step — see there.

4. **Engine agreement not measured.** Adapt
   `verify_item_6_engine_agreement.py`: solve a uniform stretch of the MMI body
   at two lengths and check the two engines agree on each mode's accumulated
   phase. The extraction helper
   (`tests/test_engine_agreement.py::dbeta_from_phases`) is device-independent
   and already tested — it takes lengths and phases, nothing coupler-specific.

   **Attempted 2026-07-31, job 11591249 — ❌ TIMED OUT, not yet answered.**
   `verify_mmi_2x2_engine_agreement.py` (new), 4 lengths (2/4/6/8 µm, deeper
   margin than the minimum two since the body has room), same
   `dbeta_from_phases` helper, same 2-mode pair gap 2's own Δβ_01 measures.
   Meant to answer whether this device's own hybrid-vs-FDTD gap (0.82%, sign
   flips across the band) has the coupler's discretization mechanism
   underneath it or something else — `calibrate_body_neff` has never been
   tried here. Real FDTD time-domain runs, not eigensolves — 2 nodes, 5h
   budget, same shape as this device's own production jobs. **That budget was
   not enough**: `sacct` confirms `State=TIMEOUT` (`ExitCode=0:0`, killed by
   SLURM, not a crash) at 19:46, the full 5h. All 4 lengths' first FDTD run
   completed (`run 0 finished` at 29,398/30,997/33,334/30,505 timesteps), but
   a 5th run — a second pass, per length or per comparison — was only at
   timestep ~5,311 when killed, well under half the length of the completed
   runs. **Gap 4 is still open**; no engine-agreement number exists yet for
   this device. This script has no `save_dir`/cache reuse the way the modal
   multiplexer's hybrid script does, so a plain resubmit likely restarts all
   4 lengths from scratch rather than resuming — worth a longer wall-clock
   budget or a checkpointing rework before trying again, not a blind rerun.
   Full accounting in `running_jobs.md`'s 11591249 entry.

   **Reshaped and resubmitted, 2026-07-31 ~21:20.** `verify_mmi_2x2_engine_
   agreement.py` now supports running one length per SLURM array task
   (`SLURM_ARRAY_TASK_ID` 0-3) plus a separate merge step that assembles the
   fit/ratio from the on-disk cache (`cache_only=True`, raises loudly if a
   task hasn't finished). Resubmitted as **11596845** (array, 8h/task —
   generous given L=4 alone took 182.7 min and cost grows faster than
   linearly with L) with merge job **11596851** dependent on it. See
   `running_jobs.md`'s reshape note and the script's own docstring for the
   mechanics.

   ---

   > **✅ CLOSED 2026-08-01 — the engines agree to 0.47%; `calibrate_body_neff`
   > is NOT indicated on this device.** Jobs 11596845 `[0-3]` (array) +
   > 11596851 (merge, 25:29). Result in
   > `example_output/verification/mmi_2x2_engine_agreement/engine_agreement.npz`
   > (`beef_git` 9e52e9d, res 40, 2 modes, 11 wavelengths, lengths 2/4/6/8 µm).
   >
   > | quantity | value |
   > |---|---|
   > | mean ratio Δβ_FDTD / Δβ_EME | **0.9953** |
   > | ratio range across 1.530–1.565 µm | 0.9962 → 0.9945, **monotone decreasing** |
   > | Δβ_FDTD | 0.16141 → 0.16597 |
   > | Δβ_EME | 0.16204 → 0.16690 |
   > | worst \|off-diagonal\| of T, all lengths/wavelengths | 1.21e-05 |
   > | worst FDTD phase-fit residual | 4.02e-06 rad |
   >
   > Both diagnostics print as `0.0000` in the job log at the script's `:.4f`
   > — they are genuinely that small, read off the `.npz`, not defaulted or
   > skipped. The section is diagonal as the method requires (**checked, not
   > assumed** — 1.21e-05 is 4 orders below the diagonal's 1.000), and the
   > 4-length linear fit of dφ(L) is clean to 4e-06 rad, i.e. **25,000× inside
   > the 0.1 rad guard** whose failure is what broke the 1×2's calibration
   > attempt (gap 5 there). This is the well-behaved counterpart to that
   > failure, on the same helper.
   >
   > **Reading, per the script's own stated criterion** (`ratio ~ 1.00` ⇒ the
   > mechanism is absent): 0.47% against the coupler's ~10% is **21× tighter**,
   > and it is consistent with the coupler doc's own prediction — the coupler's
   > figure is `β/Δβ ≈ 99` amplification of a ~2.3% per-mode n_eff error
   > between two near-degenerate supermodes (`directional_coupler_
   > validation.md` item 6 / item 13 note), and this device's two modes are
   > split far more strongly, so there is no comparable lever here.
   >
   > **This does NOT explain gap 5's 0.82%, and that is the useful part.** The
   > residual here is a *monotone* bias — 0.9962 → 0.9945 in one direction with
   > no zero crossing — whereas gap 5's Δ T **flips sign across the same band**.
   > A monotone Δβ bias cannot produce a sign-flipping transmission error, so
   > gap 5's discrepancy has to come from somewhere other than body phase
   > accumulation (taper transitions or the cascade itself are the obvious
   > candidates; neither is measured). Flagged as an **argument from the two
   > signatures, not a measurement** — nothing has isolated the taper
   > contribution on this device.
   >
   > **Caveat carried:** res 40, which gap 2 established is *not* converged in
   > absolute n_eff (Δβ_01 needs res 100–120 for 0.05%). Both engines were read
   > at the same res 40 here, so a shared mesh bias partially cancels in the
   > ratio, but the 0.47% is not certified mesh-independent. Re-measuring at
   > res 100 was not done and is not queued.
   >
   > **Cost, for reuse:** per-length wall clock came out 96 min (L=2) / 183
   > (L=4) / 215 (L=6) / **307 (L=8)** — 13.4 h of serial work, which is why
   > no single 5 h slot was ever going to finish it, and L=8 alone would have
   > blown that limit. As an array the whole thing closed in 5:07 + 25 min.
   > `_0`/`_1` returned in **36 s each** on 11591249's cached lengths, exactly
   > as the reshape predicted. The array-plus-`cache_only`-merge shape is worth
   > copying for any other multi-length `solve_section()` sweep.

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
   (`directional_coupler_validation.md` items 6/14 precedent: ~5%
   uncalibrated, 0.11% after `calibrate_body_neff`) vs. structural.
   SLURM job 11585844 (array indices 1-2, `mmi_2x2_fdtd`/`mmi_2x2_hybrid`)
   generates the two SMatrix files; job 11585849
   (`--dependency=afterany:11585844`) runs the comparison once they land.
   Note this runs *before* gap 2 (resolution) despite the Order section
   below — both examples still share `RESOLUTION = 40.0`, so this first
   pass measures the hybrid at the same not-yet-converged resolution gap 2
   already flags, not a clean result.

   **✅ MEASURED 2026-07-31 — job 11585849 (18 s), and the signature is
   different from the 1×2's.** `check_passive()` true both sides.

   | port | max\|ΔT\| | mean\|ΔT\| | % of mean T_fdtd | worst λ |
   |---|---|---|---|---|
   | bar | 0.0028 | 0.0012 | 0.62% | 1.5615 µm |
   | cross | 0.0043 | 0.0017 | 0.82% | 1.5650 µm |

   > ### ✅ RE-MEASURED DISPERSIVE 2026-08-01 — jobs 11609211_1 + 11609680 + 11609681
   >
   > **This is the number the shipped configuration produces**, and it
   > supersedes the frozen-era table above as the current gap 5 result. The
   > frozen numbers are kept because they are the other half of the P4 test.
   >
   > | port | max\|ΔT\| frozen | max\|ΔT\| **dispersive** | % of mean T_fdtd |
   > |---|---|---|---|
   > | bar | 0.002803 | **0.002689** | 0.62% → **0.59%** |
   > | cross | 0.004279 | **0.003951** | 0.82% → **0.76%** |
   >
   > **Everything qualitative in the paragraphs below survives the flip
   > unchanged** — still under 1% on both ports, still one sign change per port
   > in the same grid interval, still an order of magnitude tighter than the
   > 1×2. The material change moved the FDTD reference by at most 0.00033 and
   > the hybrid by at most 0.00011.
   >
   > It also **answers the "not yet tried" note two paragraphs down, in the
   > negative**: the sign-changing slope is *not* a dispersion-sensitive phase
   > term. Full analysis, and the refutation of the P4 mechanism it was the
   > decisive test for, in the **P4** box under *Parity with the coupler*.

   Under 1% on both ports — an order of magnitude tighter than the 1×2's 9.11%
   — but **the sign is not uniform**: both ports start positive (hybrid over
   full-device FDTD) at the blue edge, cross zero near 1.544 µm, and end
   negative at the red edge. That is *not* the 1×2's signature (`mmi_1x2_
   validation.md` gap 5 — every ΔT negative, systematic under-prediction, the
   Δβ-discretization pattern items 6/14 describe on the coupler). A slope that
   changes sign across the band looks more like a resolution- or
   dispersion-sensitive phase term than a single constant offset, so
   `calibrate_body_neff=True` — which corrects one constant bias in the body's
   propagated phase — is not obviously the fix here the way it is for the 1×2.
   **Not yet tried on this device regardless**; worth one run to see whether it
   moves the crossing point or the magnitude at all, but expect a smaller or
   different effect than the 1×2's, not assume the same win.

   Also unresolved and probably entangled with the above: this measurement is
   at the shipped, not-yet-converged `RESOLUTION = 40` (gap 2), so part of the
   sign flip could be a resolution artifact rather than physics — cannot be
   separated from calibration's effect until gap 2 is measured on this device.
   Output under `example_output/verification/mmi_2x2_hybrid_vs_fdtd/`.

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
     `resolved/fdtd_vs_EME_mode_validation.md` item 11 shows the meep-side guard's
     duplicate check has no discriminating margin left at these mode counts on
     a 2.8 µm body — which is exactly this device's body width. If that guard
     mis-fires anywhere, here is a likely place.
   - `N_MODES` also moved 15→14 since the audit (gap 3), so the recorded table
     describes a request this device no longer makes.

   Cheap to close: re-run `audit_interface_modes.py mmi_2x2` at `N_MODES = 14`
   with the guards on and confirm the identity permutation and worst
   |Δn_eff| 0.0007 still hold. Any guard firing in that log is itself the
   finding. ~~Until then, treat "the cleanest device audited" as a 2026-07-26
   measurement of a superseded configuration, not a current claim.~~

   **✅ CLOSED 2026-07-31 — job 11588569 (20:18), and it holds up exactly.**
   `batch/verification/run_audit_mmi_2x2_gap6.batch`, `N_MODES = 14`, res 40,
   λ = 1.55 µm, both duplicate guards live. The audit's own config had still been
   asking for 15 while claiming to be "pinned to `mmi_2x2_hybrid.py`"; corrected
   to 14 as part of this run.

   | | 2026-07-26 (job 11500536) | 2026-07-31 (job 11588569) |
   |---|---|---|
   | modes requested | 15 | **14** |
   | guided | 14 of 15 | **14 of 14** |
   | correspondence | identity `[0 … 13]` | **identity `[0 … 13]`** |
   | worst guided \|Δn_eff\| | 0.0007 | **0.0007** |
   | duplicate guards | did not exist | **both live, neither fired** |

   Three things this settles:

   - **"The cleanest device audited" is a current claim again**, not a dated one.
     `README.md` Findings 1 and 2 and the three other files that quote 0.0007 can
     be read at face value.
   - **Neither guard fired — zero firings in the log.** That is the specific risk
     this gap raised: `resolved/fdtd_vs_EME_mode_validation.md` item 11 identifies a 2.8 µm
     body at high mode count as exactly where the meep-side guard has the least
     discriminating margin, and this device is the repo's highest-mode-count
     SOI_220 geometry on precisely that body width. It had nothing to fire on.
     Combined with item 11 closing the same day, the "guard might mis-fire here"
     concern is answered rather than merely unobserved.
   - **Dropping the unguided band cost nothing measurable.** Worst guided
     |Δn_eff| is unchanged at 0.0007 with band 15 gone, which is what item 9
     predicted ("at no cost") and had not previously been checked on this device.

   Per-band table in `reports/11588569.out`; every band 0–13 is guided, worst
   disagreement 0.0007 at band 12 (meep 1.570726 / meow 1.571412).

## Order

~~**3** first~~ — done 2026-07-26, passed. Remaining:

1. ~~**Gap 6** (re-audit at `N_MODES = 14` with the duplicate guards on) — new
   2026-07-30 and now first, because it is a smoke run and because everything
   below quotes gap 3's numbers as though they still describe this device.~~
   **✅ Done 2026-07-31, job 11588569 — passed unchanged (identity, 0.0007,
   14/14 guided, neither guard fired).** See gap 6.
2. ~~**Gap 2**, the resolution question — **now first.** It gates gap 1, and
   the audit has already shown res 40 is not converged.~~ **✅ Done
   2026-07-31, job 11590785 — Δβ_01 (the quantity the physics actually
   depends on) converges to 0.05% by res 120; gap 2 closes at res 100.**
   See gap 2.
3. ~~**Gap 5** (hybrid vs full-device FDTD), which needs no new scripts.~~
   **✅ Done — 0.82%, sign flips across band.** See gap 5.
4. **Gap 1** (goldens) once the configuration has stopped moving. ~~Dropping
   `N_MODES` to 14 per gap 3 counts as the configuration moving, so do it
   before recording anything.~~ **That drop is done** (gap 3, 2026-07-30), and
   now so is gap 2 above. The remaining blocker is the repo-wide
   material-model regeneration alone (`README.md` Finding 6) — no goldens
   anywhere in the repo have been regenerated since the 2026-07-28 Sellmeier
   change, and recording one here first would need doing twice.
5. ~~**Gap 4** last, as on the coupler. **Attempted 2026-07-31 — job 11591249
   TIMED OUT at the full 5h budget, roughly halfway through its FDTD work.
   Still open**, needs a reshaped (longer or checkpointed) resubmit.~~
   **✅ Done 2026-08-01, jobs 11596845 `[0-3]` + 11596851 — ratio 0.9953, the
   engines agree to 0.47%, `calibrate_body_neff` not indicated here.** See
   gap 4.

**Gap 1 is now the only open item on this device**, and it is not gated on
anything 2×2-specific any more — the blocker is the repo-wide material-model
regeneration (`README.md` Finding 6), shared with every other device. Note
gap 4's answer *removes* a reason to delay: the configuration will not change
on account of calibration, because calibration is not indicated here.

## Parity with the coupler — two items beyond the four checks

Added 2026-08-01. All four checks are ticked on this device and it is still
not at coupler level. Both items are tracked as **P3** and **P4** in
[README.md](../README.md)'s parity backlog; recorded here because this doc owns
them.

> ### 🔬 P4 — a candidate mechanism, 2026-08-01, with a decisive test in flight
>
> **⚠ READ FIRST: this block is the pre-test prediction, preserved unedited.
> The test has since run and REFUTED the mechanism — see the ❌ section
> immediately following.** Kept in full because a prediction written before
> the result is the only kind that can be checked afterwards.
>
> Recorded **before** the test lands, so the prediction is on record rather
> than fitted afterwards.
>
> **The observation.** Reading gap 5's per-wavelength ΔT out of
> `comparison.npz` (preserved as `comparison_FROZEN_pre20260801.npz`), the sign
> change is not scattered — it sits **at band centre**:
>
> ```
>    λ        ΔT bar     ΔT cross
>  1.5300    +0.00069    +0.00249
>  1.5405    +0.00066    +0.00103
>  1.5440    -0.00001    +0.00053
>  1.5475    -0.00066    -0.00048   <-- band centre
>  1.5650    -0.00210    -0.00428
> ```
>
> Bar crosses between 1.5405/1.5440, cross between 1.5440/1.5475. Band centre
> is **1.5475**.
>
> **The mechanism.** Until 2026-08-01 a broadband FDTD run froze one index at
> band centre. That error is **exactly zero at centre and antisymmetric about
> it**. Gap 5 compares a hybrid whose body is EME (exact per wavelength)
> against a full-device FDTD reference whose body is frozen — so their
> difference contains precisely such a term, and it must change sign where the
> frozen index crosses the true one. It does.
>
> **Why gap 4 stays monotone under the same artifact — this is the part that
> dissolves the tension.** Gap 4 measures Δβ, a *difference* of two modes
> solved in the same frozen medium, so the frozen error largely cancels. Gap 5
> depends on *absolute* accumulated phase through the body, where it does not.
> The two checks are not in conflict; they have different sensitivity to a term
> neither was designed to measure.
>
> **Tested against the 1×2, which should NOT show it.** Confirmed: that
> device's ΔT is uniform-negative with no crossing, and its magnitude
> (0.017–0.036) is **~10× the 2×2's ~0.003 swing** — a real offset that swamps
> any antisymmetric term. Consistent with the mechanism rather than a
> counter-example, and the asymmetry is itself predicted: the term is only
> visible where the genuine gap is small, which is exactly the 2×2.
>
> **The decisive test, submitted 2026-08-01: jobs 11609211 `[1-2]` +
> 11609212** re-measure gap 5 with dispersion on, where both sides are exact
> and the term must vanish.
>
> | outcome | reading |
> |---|---|
> | sign flip **disappears**, \|ΔT\| shrinks | mechanism confirmed; **P4 closes**, and checks 1 and 4 reconcile after all |
> | sign flip **survives** | mechanism refuted; the residual is structural (tapers or the cascade) and P4 stands as originally written |
>
> Either way this is worth more than the "not isolated, no cheap first step"
> the text below records — the step was cheap and it existed. Flagged as a
> **hypothesis with evidence**, not a result: one zero crossing at the right
> place is suggestive, and bar's crossing is a grid step early.
>
> **First hard number, and it cuts both ways — 1×2 job 11609215_1, landed
> before this device's.** That device's full-device FDTD reference was
> re-measured dispersive against its preserved frozen baseline
> (`../mmi_1x2_validation.md` gap 5 for the full table). Two findings:
>
> - **Scale supports the mechanism.** The reference moves by up to **0.0065**
>   (1.65% of mean T) purely from the material change. That is ~20% of the
>   1×2's own 0.021–0.036 discrepancy — hidden there, as the falsification
>   control predicted — but it is **2× this device's entire ~0.003 ΔT swing.**
>   A term of this class cannot be a minor correction here; it either dominates
>   or the devices differ in a way not yet identified.
> - **Shape does not.** The 1×2's shift is near-zero at 1.3001 (band centre
>   1.310, same slightly-below-centre offset seen here), but it crosses zero a
>   **second** time at 1.3514 and returns positive at the band edge. A single
>   frozen-index term cannot produce two crossings. So the simple
>   "antisymmetric about the freeze point" picture is **incomplete at best**,
>   and the prediction table above should be read as testing *whether the
>   material term dominates*, not as testing a fully specified shape.
>
> Recorded before 11609212 landed, and deliberately including the half that
> weakens the hypothesis.

### ❌ P4's candidate mechanism is REFUTED — 2026-08-01, jobs 11609680 + 11609681

**The test ran and the prediction table's second row is what happened: the
sign flip survives.** Read the table above before this paragraph — it was
written before the job landed and it is unedited.

Jobs 11609211_1 (dispersive full-device FDTD, 1:51:57) + **11609680**
(dispersive hybrid, 3:50:41, the `-t 720` resubmission of the cancelled
11609211_2) + **11609681** (comparison, 0:15). Both sides dispersive and
mtime-verified against `provenance.json` — `fdtd_pkl` 19:13 UTC, `hybrid_pkl`
21:59 UTC, beef `f8b0d47+dirty`, against the frozen baseline's `28d65f8+dirty`.
No repeat of the 11608692 pickle race.

Per-wavelength ΔT, read from the `.npz` (not the `:.4f` log lines):

| λ (µm) | ΔT bar frozen | ΔT bar disp. | ΔT cross frozen | ΔT cross disp. |
|---|---|---|---|---|
| 1.5300 | +0.000686 | +0.001011 | +0.002492 | +0.002659 |
| 1.5405 | +0.000658 | +0.000743 | +0.001025 | +0.001058 |
| 1.5440 | **−0.000008** | **+0.000029** | +0.000527 | +0.000549 |
| 1.5475 | −0.000662 | −0.000646 | **−0.000479** | **−0.000491** |
| 1.5650 | −0.002105 | −0.001942 | −0.004279 | −0.003951 |

**Exactly one sign change per port, before and after, in the same grid
interval.** The mechanism required the crossing to *vanish*; it did not move.

| | max \|ΔT\| frozen | max \|ΔT\| dispersive | rel to mean T |
|---|---|---|---|
| bar (`out_0`) | 0.002803 | **0.002689** | 0.616% → **0.591%** |
| cross (`out_1`) | 0.004279 | **0.003951** | 0.824% → **0.761%** |

|ΔT| shrank by 4% and 8% — a rounding-level improvement, not the collapse a
dominant material term would have produced.

**The decisive number is how little the material change moved either side.**
Differencing the dispersive and frozen runs directly:

| | max shift, FDTD reference | max shift, hybrid | max change in ΔT | as fraction of the ΔT swing |
|---|---|---|---|---|
| bar | 0.000329 | 0.000114 | 0.000325 | **8.6%** |
| cross | 0.000298 | 0.000030 | 0.000328 | **4.8%** |

**This also resolves the "cuts both ways" tension recorded above, against the
hypothesis.** The 1×2 evidence said the reference moves 0.0065 — *2× this
device's entire ΔT swing* — so a term of that class had to either dominate here
or the devices had to differ. **They differ, by a factor of ~20:** the 2×2's
reference moves at most **0.00033**. So the scale argument never applied to
this device, and the one measurement that looked like support for P4 was in
fact measuring the 1×2's much larger material sensitivity.

That difference is not mysterious — the 1×2 spans 1.26–1.36 µm (±3.8% about
centre) against the 2×2's 1.53–1.565 µm (±1.13%), a 3.4× narrower relative
band, and freezing an index costs roughly in proportion to how far the band
runs from the freeze point. Offered as the plausible reading, not a measured
decomposition.

**And the frozen term's *shape* was wrong too.** The mechanism's premise was
that the reference's frozen error is zero at centre and antisymmetric about
it. Measured, `d(FDTD ref)` on bar is negative at **every** wavelength
(−0.00033 … −0.00009, no crossing at all), and on cross it wanders sign
several times with no structure about 1.5475. It is a small noisy offset, not
an antisymmetric term. The 1×2's two-crossing shape, already flagged above as
weakening the hypothesis, was the first sign of this.

**Conclusion: the residual is structural — tapers or the cascade — exactly as
P4 stood before the hypothesis was raised.** P4 does **not** close. What the
test bought is the elimination of the cheapest candidate, and it cost one
already-required re-measurement, since gap 5 needed its dispersive number
regardless. The hypothesis was worth recording and worth testing; it was
wrong, and the prediction being on record beforehand is what makes that
readable now.

> **Corrected in the text above:** the pre-test note "*Scale supports the
> mechanism*" should be read as refuted — it extrapolated the 1×2's material
> sensitivity to this device, and that extrapolation is what failed.

**P4 — checks 1 and 4 do not reconcile, and this is the substantive one.**
Gap 4 (check 1) measures a **monotone** 0.47% Δβ bias, 0.9962 → 0.9945 across
the band with no zero crossing. Gap 5 (check 4) measures a 0.82% ΔT whose
**sign flips** across that same band. A monotone Δβ bias cannot produce a
sign-flipping transmission error, so **gap 4 does not account for gap 5's
residual.** On the coupler these two reconciled — item 7's measured ratio
0.9457 against the 0.9536 predicted from item 6, agreeing to 0.82% from two
completely independent routes — and that mutual confirmation is a large part
of why the coupler's number is trusted. The 2×2 has no second route and its
two checks actively disagree in signature. The cause lies somewhere the four
checks between them do not cover: the taper transitions, or the cascade
itself. **Neither has been isolated, and there is no cheap first step**;
gap 4's own uniform-body cut is by construction blind to both.

**P3 — the resolution knob is still shared.** `mmi_2x2_hybrid.py` ships one
`RESOLUTION = 40.0` used by both sub-solvers, while the coupler splits its own
(`FDTD_RESOLUTION = 40` / `EME_RESOLUTION = 80`) with the comment *"shared
RESOLUTION was how a 7% error hid for months (dc_validation.md item 11)"*.
Gap 2 measured this device's Δβ_01 as needing **res 100**. Keeping 40 is a
defensible cost call and is recorded as a caveat — but it means check 3 is not
closed here the way it is on the coupler, and gap 4's 0.47% inherits the
caveat (both engines read at the same mesh, so a shared bias partially
cancels in the ratio, but the number is not certified mesh-independent).
**Worth re-deciding explicitly rather than leaving implicit.**
