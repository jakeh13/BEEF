# Active SLURM jobs — snapshot

**Snapshot taken 2026-07-31, 13:05 EDT** (seventh snapshot of the day: 07-30
21:46 → 07-31 08:22 re-queue → 08:46 → 10:35 → 10:53 → 12:05 → now). This is a point-in-time list, not a
live view — job states, elapsed times and queue positions all move. Check
`squeue -u jhiesener3` (or `sacct -j <id>` once finished) for current status;
this file exists so the *purpose* of each job is recorded somewhere durable,
not to be a dashboard. Update or delete entries as jobs complete and their
results get folded into the relevant validation doc.

**Re-queue, 2026-07-31 08:22 EDT:** the five jobs below never started —
`scontrol show job` confirmed `Account=gts-sr108-ece`,
`Reason=AssocGrpBillingMinutes` on the three parents (that association is out
of billing minutes; `gts-sr108-gedc_main` had headroom). Cancelled all five
(old IDs 11577985/11577986/11578018/11578019/11578026) and resubmitted under
`gts-sr108-gedc_main`. Also switched all other batch scripts hardcoding
`gts-sr108-ece` to `gedc_main` so this doesn't recur.

| job | name | state (at snapshot) | for |
|---|---|---|---|
| 11575516 | `beef_mode_spec_converge` | ✅ **COMPLETED** 21:17, 3:36:13 — results folded in, entry retained for provenance | [fdtd_vs_EME_mode_validation.md](fdtd_vs_EME_mode_validation.md) item 1 + [modal_multiplexer_validation.md](modal_multiplexer_validation.md) gap 2 |
| ~~11577985~~ → **11585844** `[1-2]` | `beef_mmi_2x2_examples` | ✅ **COMPLETED** — `_1` 1:16:45, `_2` 4:11:55; both SMatrix files on disk | [mmi_2x2_validation.md](mmi_2x2_validation.md) gap 5 (data) |
| ~~11577986~~ → **11585849** | `beef_mmi_2x2_gap5` | ✅ **COMPLETED** 12:46, 0:18 — **worst max\|ΔT\| 0.0043 (0.82%), sign FLIPS across the band (not the 1×2's uniform signature); results folded in** | [mmi_2x2_validation.md](mmi_2x2_validation.md) gap 5 (comparison) |
| ~~11578018~~ → **11585845** `[1-2]` | `beef_mmi_1x2_examples` | ✅ **COMPLETED** — `_1` 8:57, `_2` 1:46:59; both SMatrix files on disk | [mmi_1x2_validation.md](mmi_1x2_validation.md) gap 5 (data) |
| ~~11578019~~ → **11585850** | `beef_mmi_1x2_gap5` | ✅ **COMPLETED** 10:28, 0:09 — **max\|ΔT\| 0.0359 (9.11%), results folded in** | [mmi_1x2_validation.md](mmi_1x2_validation.md) gap 5 (comparison) |
| ~~11578026~~ → **11585847** `[0]` | `beef_modal_mux_examples` | 🟢 **RUNNING** (4:22 of 5h, 4 nodes) — meow-side guard now firing **6×/2 pairs** on the production path | [modal_multiplexer_validation.md](modal_multiplexer_validation.md) gap 6 (exercise run) |
| ~~**11586046**~~ → **11588309** | `beef_diagnose_items_11_12` | 11586046 ❌ **CANCELLED** 11:12 (hung, not slow); refixed → **11588309 ✅ COMPLETED 11:50, 34:55 — BOTH ITEMS RESOLVED:** item 11 closed (all four grids agree at res 40 too), item 12 decisive (laddering recovers a 14th guided mode at equal cost) | [fdtd_vs_EME_mode_validation.md](fdtd_vs_EME_mode_validation.md) items 11 + 12; it also opened item 14 |
| **11586098** | `beef_diagnose_item_13` | ✅ **COMPLETED** 10:50, 1:54:50 — **C vs A inert (max \|Δn_eff\| 0.000e+00), cost 3.2×; results folded in** | [fdtd_vs_EME_mode_validation.md](fdtd_vs_EME_mode_validation.md) item 13 |
| ~~**11586986**~~ → **11588997** | `beef_dispersive_convergence` | 11586986 ❌ **CANCELLED** 12:00 — **it was measuring nothing**, `dispersive=` never reached the mode solve (library bug, fixed); refixed → **11588997 🟢 RUNNING** (1:03 of 8h), liveness check passed at 7.87e-03 | [material_model_validation.md](material_model_validation.md) work item 4 (cheap half) |
| **11588569** | `beef_audit_mmi_2x2_gap6` | ✅ **COMPLETED** 11:58, 20:18 — **gap 6 closed:** identity `[0…13]`, worst \|Δn_eff\| 0.0007, 14/14 guided, **neither guard fired** | [mmi_2x2_validation.md](mmi_2x2_validation.md) gap 6 |
| **11589770** | `beef_mmi_1x2_hybrid_calibrated` | 🟢 **RUNNING** (21:24 of 5h) | [mmi_1x2_validation.md](mmi_1x2_validation.md) gap 5 follow-up (calibration) |
| **11589779** | `beef_mmi_1x2_gap5_calibrated` | PENDING (dep on 11589770) | [mmi_1x2_validation.md](mmi_1x2_validation.md) gap 5 follow-up (comparison) |
| ~~11589771~~ ~~11589899~~ ~~11589936~~ | `beef_modal_mux_fdtd` | ❌ **BLOCKED, do not resubmit** — three cheap failures (~15 s each) converged on a real library gap: `CompoundGeometry` cannot back `FDTDSolver.solve()` at all | [modal_multiplexer_validation.md](modal_multiplexer_validation.md) gap 3 |

**Six items resolved today; three jobs cancelled as unproductive, two refixed
and rerun, one left blocked.** Since the 12:05 snapshot: the **2×2 MMI's whole
gap-5 chain finished** (11585844_2, then 11585849) at 0.82% worst |ΔT| — an
order of magnitude tighter than the 1×2's 9.11%, but with a sign that flips
across the band rather than the 1×2's uniform under-prediction, so
`calibrate_body_neff` is not obviously indicated the same way here. A
follow-up calibration run for the **1×2** (which does show the uniform
signature) is in flight — 11589770/11589779. And **modal_multiplexer_fdtd.py**,
written to give the mux its first full-device FDTD reference (gap 3), hit a
real, unfinished library gap: `CompoundGeometry` cannot back
`FDTDSolver.solve()`'s single-domain path at all (only `solve_section()`,
what `HybridSolver` uses). One sub-bug in the same area WAS fixed and tested
(`CompoundGeometry.source_aperture_at()` used to raise past its own footprint;
now clamps, matching how every leaf geometry already behaves) — necessary but
not sufficient. Three cheap failed submissions (~15 s each) before landing on
the real blocker; not resubmitting until the library gap is actually closed.
See the mux's own detail section below and `modal_multiplexer_validation.md`
gap 3 / `notes.md` Geometry item 8.

Before that: **items 11 and 12** both resolved by 11588309 (34:55,
after its predecessor hung), **2×2 MMI gap 6** closed by 11588569, and a library
bug found by reading 11586986's output — `dispersive=` never reached
`section_mode_neff()`, so that job's two columns were byte-identical and work
item 4 was never actually started. Fixed, 9 tests added, re-running as 11588997
with a liveness check that aborts if the flag is ever inert again.

Item 13 is closed by 11586098 (the
guard's over-request is a pure filter, and `_DUP_MODE_MAX_EXTRA` is now a
measured 3.2× of pure cost) — folded into
`fdtd_vs_EME_mode_validation.md`. The 1×2 MMI's whole gap-5 chain completed end
to end earlier and its numbers are in `mmi_1x2_validation.md`. `11585847` has
produced a real finding before finishing: the item-0 duplicate guard fires on
`EMESolver`'s **production** path on this device at res 40 (folded into
`modal_multiplexer_validation.md` gap 6), though the item-8 `T`-column scan it
was submitted to exercise is still silent. `11585849` (2×2 comparison) is
pending only on `11585844_2`, which is correct.

**11586046 HUNG and was cancelled at 11:12; resubmitted as 11588309 after a
script fix.** It was not slow — it was stuck, and an earlier revision of this
file diagnosed it wrongly as a wall-clock budget miss. The real cause, from its
own MPB trace: `N_TARGET = 20` distinct modes is unattainable at that
cross-section (today's fixed seed finds 13 distinct guided modes there, and a
laddered seed 14 — job 11588309), so the search
walked up to `band_num 23`, where **no k reaches the target frequency at all**
and MPB's k-search limit-cycles between k = 0.2346 and k = 1.7480 indefinitely.
It burned 161 MPB solves on that one band — 59% of the job's entire eigensolve
work — while bands 1–22 had each converged in 3–8. Bands 1–22 took ~58 min for
22 attempts, so the original ~25 s/solve budget was about right: wall clock was
never the problem, and the 4-task split this file recommended at 11:03 would
have produced four hung tasks.

Two results came out of the failure, both folded into
`fdtd_vs_EME_mode_validation.md`: item 12 gains the measured band table showing
this device's shipped **16-mode basis is unachievable at res 40 by any seeding
policy** — since superseded by 11588309, which found laddering recovers a 14th;
the seeding half IS fixable, the count is 14 not 13 — and the production guard's
identical exposure to the same non-convergent solve is filed as new **item 14**.
Full detail under 11586046 and 11588309 below.

The detail sections below carry the current job IDs — for the five re-queued
jobs only the ID and billing account changed, so each entry's script, resources
and rationale are unchanged from the 07-30 submission.

## Details

### 11575516 — `beef_mode_spec_converge` — ✅ COMPLETED
- **Script:** `batch/verification/run_mode_spectrum_convergence.batch` →
  `examples/verification/python/mode_spectrum_convergence.py`
- **Resources:** 1 node, 8 cpus, serial (no MPI), 8h limit. Used 3:36:13.
- **Output:** `reports/11575516.out`.
- **What it did:** Swept meep and meow resolution *independently* (meep
  20/30/40/60, meow 20/30/40/60/80) on the modal multiplexer's interface 0
  and counted how many of the 16 requested modes each engine calls guided,
  and how many modes the two engines agree on within tolerance
  (`n_agree`), at every (meep_res, meow_res) pair.
- **Why:** discriminate the leading hypothesis for item 1's meep↔meow orphan
  mode — `n_agree` rising toward 16 with resolution would mean a
  band-window/discretization artifact (fixable with a per-engine converged
  resolution); a plateau would mean something structural.
- **Result — both, and the structural half is the one that matters.**
  `n_agree` does rise (12 → 15 at matched res 40 → 60) and the guided orphan
  count falls from 3 to 1, so most of the mode-8+ divergence is
  discretization. But the surviving orphan is present at **every** resolution
  20–60: meow finds two members of a near-degenerate pair (converging to a
  3.1e-4 splitting by res 80) where meep's per-band ladder finds one. That is
  a driver defect, not a mesh one. Full reading under item 1; the new to-dos
  it opened are items 11–13.
- **Second result, unplanned:** the first per-engine converged-resolution
  measurement on SOI_160 — meow plateaus between res 60 and 80 (16/16), meep
  is still moving at 60. The shipped shared `RESOLUTION = 40` is below both
  (`modal_multiplexer_validation.md` gap 2).
- **Correction to this entry's original filing:** it was listed as serving
  item 6 (`interface_modes`/layer-2 convergence) as well. It does not — that
  item sweeps the requested mode *count* and requires the *S-matrix* to
  plateau, neither of which this job did. See item 6's note.

### 11585844 `[1-2]` — `beef_mmi_2x2_examples` (array indices 1,2 only) — ✅ COMPLETED
- **Script:** `batch/run_mmi_2x2_examples.batch` → `make mmi_2x2_fdtd
  mmi_2x2_hybrid`
- **Resources:** 2 nodes × 24 ranks per array task, 5h limit each.
- **What it does:** Runs the 2×2 MMI's full-device FDTD reference
  (`mmi_2x2_fdtd.py`) and the FDTD/EME hybrid cascade
  (`mmi_2x2_hybrid.py`), each saving an `SMatrix` pickle. Array index 0
  (`mmi_2x2_eme`) deliberately skipped — not needed for gap 5.
- **Why:** gap 5 ("no hybrid-vs-full-device comparison on record") needs
  both scripts' saved output to exist before anything can compare them.
- **Result:** `_1` (fdtd) 1:16:45, `_2` (hybrid) 4:11:55. Both SMatrix pickles
  on disk; see 11585849 for the comparison.

### 11585849 — `beef_mmi_2x2_gap5` — ✅ COMPLETED
- **Script:** `batch/verification/run_verify_mmi_2x2_hybrid_vs_fdtd.batch` →
  `examples/verification/python/verify_mmi_2x2_hybrid_vs_fdtd.py`
- **Resources:** 1 node, 1 cpu, 20 min. `--dependency=afterany:11585844`.
- **What it does:** Loads both SMatrix pickles from 11585844, reconciles
  port names (`out_0`/`out_1` vs `bar`/`cross`), reports per-wavelength
  transmission difference with provenance (git sha, config, mtimes). No
  meep/meow import cost of its own beyond loading `beef`.
- **Why:** closes gap 5 as measurement-and-report (no threshold enforced —
  gap 4, engine agreement, isn't measured for this device yet, so there's
  no basis to call a given disagreement expected vs. a bug).
- **Result (12:46, 18 s):** `check_passive()` true both sides. Worst max|ΔT|
  **0.0043** (cross port, λ=1.5650 µm, 0.82% of mean T_fdtd); bar port 0.0028
  (0.62%). Under 1% on both — far tighter than the 1×2's 9.11% — but the sign
  is **not** uniform: both ports start positive, cross zero near 1.544 µm, end
  negative. Not the 1×2's Δβ-discretization signature (uniform negative), so
  `calibrate_body_neff` is not obviously indicated here the way it is for the
  1×2 — a sign flip looks more like a resolution/dispersion-sensitive phase
  term. Also measured at the shipped unconverged `RESOLUTION = 40` (gap 2), so
  part of the flip could be a mesh artifact rather than physics; the two
  cannot be separated until gap 2 is measured on this device. Folded into
  `mmi_2x2_validation.md` gap 5. Output under
  `example_output/verification/mmi_2x2_hybrid_vs_fdtd/`.

### 11585845 `[1-2]` — `beef_mmi_1x2_examples` (array indices 1,2 only)
- Same shape as 11585844, for the 1×2 MMI: `make mmi_1x2_fdtd
  mmi_1x2_hybrid`.
- **Why:** `mmi_1x2_validation.md` gap 5, same reasoning as the 2×2's. This
  device's basis mismatch (gap 3, `N_MODES` 12→8) is already fixed, so —
  unlike the 2×2 — there's no open mode-identity concern underneath this
  comparison; only gap 2 (resolution) remains unsettled.

### 11585850 — `beef_mmi_1x2_gap5` — ✅ COMPLETED
- **Script:** `batch/verification/run_verify_mmi_1x2_hybrid_vs_fdtd.batch` →
  `examples/verification/python/verify_mmi_1x2_hybrid_vs_fdtd.py`
- Same shape as 11585849, ported to `mmi_1x2`'s port names (`out_0`/`out_1`
  vs `out_top`/`out_bot`). `--dependency=afterany:11585845`.
- **Why:** `mmi_1x2_validation.md` gap 5.
- **Result (10:28, 9 s):** max |ΔT| **0.0359** at λ = 1.33 µm, mean |ΔT|
  0.0253, **9.11%** of mean T_fdtd. `check_passive()` true both sides; the two
  output ports agree to 4 decimals, so splitter symmetry survives the cascade.
  **Every ΔT is negative** — a systematic under-prediction, not scatter, which
  is the Δβ-discretization signature from
  `resolved/directional_coupler_validation.md` items 6/14 (coupler: ~5%
  uncalibrated → 0.11% with `calibrate_body_neff=True`, never tried here).
  Folded into `mmi_1x2_validation.md` gap 5; output under
  `example_output/verification/mmi_1x2_hybrid_vs_fdtd/`. Provenance records
  beef `1cc8edc`.

### 11585847 `[0]` — `beef_modal_mux_examples`
- **Script:** `batch/run_modal_multiplexer_examples.batch` → `make
  modal_multiplexer_hybrid`
- **Resources:** 4 nodes × 24 ranks, 5h limit — the most expensive single
  job this session; this device is the repo's "most demanding hybrid
  device" (16 modes × 4 interfaces, SOI_160).
- **What it does:** This device's **first real run** in this tree — no
  `example_output/modal_multiplexer/` existed before this job was
  submitted.
- **Why:** (a) gives the least-supported device in the repo (no golden, no
  FDTD reference, shared resolution knob) its first real output; (b)
  exercises `beef.mode_identity.transfer_matrix_duplicate_scan()` /
  `warn_on_transfer_duplicate()` (`fdtd_vs_EME_mode_validation.md` item 8)
  against this device's real `solve_section()` interfaces for the first
  time — `modal_multiplexer_validation.md` gap 6 flags this as "not yet
  exercised... needs an actual meep/meow run." Watch
  `reports/11585847_0.out` for `warn_on_transfer_duplicate` firing on
  either solver.
- **Partial result at 10:53 (2h12 in, still running):** two FDTD subregions
  done; the first EME body solve after them fired the **item-0 guard** —
  `meow mode 9 duplicates already-accepted mode 8 (self-overlap 1.000000)`,
  logged once per rank across all 96 ranks. By 12:16 this had grown to **six
  distinct firings from two pairs** — 4x mode 9/8 and 2x mode 11/10, the latter
  being item 0's own original pair. First
  catch of the meow duplicate on the *production* `compute_modes()` path rather
  than in an audit script. `warn_on_transfer_duplicate` — the item-8 scan this
  job was actually submitted to exercise — has **not** fired; no `transfer`
  line appears in the log at all, which is what you would expect if the
  upstream guard rejected the duplicate before it reached `T`, but is not yet
  evidence the scan works. Full reading in `modal_multiplexer_validation.md`
  gap 6.

### 11586046 — `beef_diagnose_items_11_12` — ❌ CANCELLED (hung), superseded by 11588309
- **Script:** `batch/verification/run_diagnose_items_11_12.batch` → two python
  scripts in sequence, `examples/verification/python/`:
  `diagnose_kpoint_ladder_v2.py` then `diagnose_dup_band_grid.py`.
- **Resources:** 1 node, 8 cpus, serial (no MPI), 240 min budget. Eigensolve
  only — `get_eigenmode()`, no time-stepping, same affordability argument as
  `run_diagnose_mode_correspondence.batch`.
- **Output:** `reports/11586046.out`.
- **Why one job for two items:** both solve the *same* cross-section (modal
  multiplexer interface 0, res 40), so they share setup cost.
- **Item 12 runs first, deliberately** — it is the root-cause item and the
  more expensive of the two, so a wall-clock overrun loses the cheaper result,
  not the more valuable one. Runs a **real** guarded accept/skip/retry search
  (not v1's flat 20-solve sweep) under four seeding policies: A fixed
  (today's code), B naive ladder (raw previous solve — what v1 tested), C step
  ladder (previous − 0.005, aimed *past* the previous band), D accepted ladder
  (last accepted band, so a skip cannot poison the next seed). C and D are
  exactly the two follow-ups item 12 lists as "not yet tried". Reports attempts
  needed, whether any policy gets stuck re-finding one duplicate, and whether
  all four converge on the same final physical mode set.
- **Item 12's open risk it is testing:** v1 showed naive laddering *moves* the
  duplicate rather than removing it (recovers 2 real modes, but re-finds the
  previous band at a new position). C and D are the two candidate fixes.
- **Item 11 second:** replays the duplicate-guard grid comparison at **res 40**
  — the single data point item 11 names as missing before it can close, since
  its clear result so far was measured only at res 20 and a finer mesh packs
  more oscillation into the same core width. `diagnose_dup_band_grid.py`'s
  `RES` is env-overridable (`BEEF_DIAG_RES`, default 20.0) so the res-20 result
  on record stays reproducible as the script's default.
- **❌ Cancelled 11:12 after 2h30. It was HUNG, not slow** — and the 10:53
  revision of this entry got that wrong, calling it a wall-clock budget miss.
  Zero of the four `===== policy =====` blocks ever printed. Reconstructing the
  converged `k_x` per band from the MPB trace explains it: at this cross-section
  bands 1–15 are guided (n_eff 2.7064 → 1.4848, cladding 1.446804) with
  duplicate pairs at 10/11 and 14/15, bands 16–22 are *below* the cladding
  index, and **band 23 has no k that reaches the target frequency at all**.
  `N_TARGET = 20` distinct modes therefore needs band_num 23, and MPB's k-search
  there limit-cycles between k = 0.2346 (f = 0.4775) and k = 1.7480 (f = 1.0283)
  forever — **161 solves on that one band, 59% of the job's entire eigensolve
  work**, against 3–8 solves for each of bands 1–22. No wall clock would have
  been enough, and the 4-task split recommended at 11:03 would have produced
  four hung tasks. Bands 1–22 took ~58 min for 22 attempts, so the original
  ~25 s/solve budget was about right.
- **What it produced anyway, before dying.** The band table above: 15 guided
  band slots at res 40 containing only **13 distinct** modes under today's fixed
  seed, while the device requests `N_MODES = 16`. Read at the time as "the basis
  is unachievable by *any* seeding policy"; the re-run showed a laddered seed
  reaches **14**, so the ceiling is 14 and part of the shortfall is a software
  defect. Folded into
  `fdtd_vs_EME_mode_validation.md` item 12, with the caveat that whether (10,11)
  is a solver duplicate or a real pair meep cannot separate is still item 1's
  open question.
- **And it opened item 14.** `_iter_distinct_bands()` (`beef/solvers/fdtd.py:1049`)
  walks `band_num` up to `n_modes + 8` and its `RuntimeError` can only fire
  *after* a solve returns, so the same non-convergent eigensolve would hang a
  production `HybridSolver` run with no exception and no timeout. Latent, not
  observed; recorded rather than fixed, per the same
  scope argument gap 6 used for this method.

### 11588309 — `beef_diagnose_items_11_12` (replaces 11586046) — ✅ COMPLETED
- **Submitted** 2026-07-31 11:15, started immediately. Same script as 11586046
  with three changes, all of them consequences of that job's hang:
  - **`diagnose_kpoint_ladder_v2.py` stops at the guided edge.** The search now
    ends when a raw solve returns `n_eff ≤ cladding index` — nothing above is
    guided, and the next band_num is the non-convergent one. A *physical* bound,
    not a timeout: MPB runs in C, so a Python signal cannot interrupt a solve
    already entered.
  - **`N_TARGET` 20 → 16**, what the device actually requests. The 20 was for
    comparability with v1's raw sweep, whose bands 16–20 were unguided anyway.
  - **Item 11 now runs FIRST.** Putting item 12 first to protect the
    higher-value result is exactly what cost both items last time: item 12 hung
    and printed nothing, and item 11 never started. Item 11 is ~50 min and every
    band it solves is known to converge, so it banks a guaranteed result; item 12
    prints per policy, so a late overrun there degrades gracefully.
- **Resources:** 1 node, 8 cpus, serial, **360 min** (was 240). Budget re-derived
  from 11586046's own numbers rather than v1's: ~2.6 min/attempt at res 40,
  ~16–18 attempts per policy × 4 ≈ 3 h for item 12, ~50 min for item 11.
- **Output:** `reports/11588309.out`.
- **Termination verified before submitting**, offline and meep-free: replayed
  11586046's measured n_eff sequence (including all three duplicate pairs)
  through `_guarded_search` with fake `sim`/`mp` objects and the real
  `field_overlap_matrix`. All four policies stop at band 16 on the cladding
  test, none requests the non-convergent band, and policy A reproduces the
  measured 13 accepted / 16 attempts / skips at band 11 and 15.
- **What to watch:** whether any policy reports `stopped because: attempts`
  (would mean duplicates ate the whole +8 headroom production also relies on),
  and whether C or D recovers more than policy A's 13 distinct guided modes.
  **Answered: all of B, C and D recover a 14th (n_eff 1.512589) at the same 16
  attempts; none hit the attempts backstop. See item 12.**
- **✅ Result (11:50, 34:55 — item 11 took 14:56, item 12 19:42). Both items
  resolved, and the reshape worked exactly as intended:** item 11 ran first and
  banked its result, item 12 terminated on the cladding test at band_num 16 under
  every policy rather than walking into the non-convergent region.
  - **Item 11 → CLOSED.** All four diagnostic grids (16×8 shipped through
    128×16) skip the same two band_nums (11, 15), accept the same 16, and give
    identical n_eff. Both skips are exact n_eff ties, so they are confirmed
    genuine without any field sampling. Worst off-diagonal overlap among accepted
    bands on the shipped grid **0.3877** (res 20 was 0.3392) against a 0.9
    threshold — the margin narrowed with resolution as the item predicted, and is
    still less than half way to firing.
  - **Item 12 → decisive.** Policy A (today's fixed seed) finds **13** distinct
    guided modes; B, C and D each find **14**, the extra being n_eff 1.512589,
    *absent from A entirely*. Same 16 attempts, one fewer skip — a strict
    improvement, so the fix direction is settled. C and D returned byte-identical
    results to B, so their extra machinery buys nothing here.
  - **Correction this forced:** the earlier reading that "only 13 distinct guided
    modes exist at res 40, so no seeding policy can help" was inferred from the
    fixed-seed trace — the very enumeration under test — and is wrong. Corrected
    in six places.

### 11588569 — `beef_audit_mmi_2x2_gap6` — ✅ COMPLETED
- **Script:** `batch/verification/run_audit_mmi_2x2_gap6.batch` (new 2026-07-31)
  → `examples/verification/python/audit_interface_modes.py mmi_2x2`.
- **Resources:** 1 node, 8 cpus, serial, 120 min budget. Used **20:18**.
- **Output:** `reports/11588569.out`, plus
  `example_output/verification/interface_mode_audit/`.
- **Why its own batch file** rather than the existing
  `run_audit_interface_modes.batch`: that one runs all three devices *with*
  `--convergence` (an EME mesh sweep, ~4× the mode solves). Gap 6 asks only
  whether the identity permutation and worst |Δn_eff| still hold at res 40 with
  the guards live — the mesh sweep is gap 2's separate question.
- **Also fixed on the way in:** the audit's own `mmi_2x2` config claimed to be
  "pinned to `mmi_2x2_hybrid.py`" while still requesting `n_modes=15`. The device
  moved to 14 on 2026-07-30 (item 9), so the audit had been describing a request
  the device no longer makes. Corrected to 14.
- **✅ Result — passed unchanged, and nothing fired.** 14/14 guided (band 15, the
  sub-cladding one, is gone as intended), correspondence identity `[0…13]`,
  **worst guided |Δn_eff| 0.0007** — the same number as the 2026-07-26 audit at
  15 modes. **Zero duplicate-guard firings** in the log from either engine, which
  answers the specific concern gap 6 raised: this is the repo's highest-mode-count
  SOI_220 device on the 2.8 µm body that `fdtd_vs_EME_mode_validation.md` item 11
  names as the guard's least-margin case, and it had nothing to fire on. "The
  cleanest device audited" is a current claim again.

### 11586098 — `beef_diagnose_item_13` — ✅ COMPLETED
- **Script:** `batch/verification/run_diagnose_guard_over_request.batch` →
  `examples/verification/python/diagnose_guard_over_request.py`.
- **Resources:** 1 node, 8 cpus, serial (no MPI), 360 min budget. meow/MPB
  cross-section solves only — no propagation, no time-stepping.
- **Output:** `reports/11586098.out`.
- **Cross-section:** modal multiplexer interface 0 again — same as 11586046 and
  the item-0 duplicate — but at **res 60**, per item 13: above the shipped 40,
  and where 11575516 found meow converged (16/16 guided, plateau 60→80), so a
  difference seen here is not an unconverged mesh talking.
- **Why three arms, not "guard on vs off":** item 13's own point is that
  flipping the flag is *not* a clean A/B — with the guard on,
  `compute_modes()` over-requests `num_modes + _DUP_MODE_MAX_EXTRA` (16 → 24)
  in one `compute_modes_meep()` call, so the flag changes both the
  eigenproblem and whether filtering happens. The arms separate them:
  **A** guard off / 16 (baseline), **B** guard on / 16 (production default),
  **C** guard off / 24 (B's eigenproblem, unfiltered).
- **What each comparison answers:**
  - **C vs A** — does over-requesting *alone* move the low bands? This is the
    artifact question. Non-zero ⇒ every guarded EME number in the repo carries
    an over-request artifact, whether or not a duplicate was ever filtered.
    Inert ⇒ the guard is a pure filter and `_DUP_MODE_MAX_EXTRA` becomes a free
    cost knob.
  - **B vs C** — what the filter actually drops: which of C's 24 the guard
    rejected as duplicates, and what it reached further up the spectrum to get.
- **Cost is half the point:** wall time is recorded per arm, so if C vs A comes
  back inert the `C/A` ratio is the whole argument for tuning
  `_DUP_MODE_MAX_EXTRA` down — pure cost, no effect on the answer. Item 13
  notes job 11575516's 3h36m was dominated by exactly this over-request.
- **Result (10:50, 1:54:50) — inert, and expensive.** **C vs A: max |Δn_eff| =
  0.000e+00** on all 16 low bands, so the over-request does not perturb the
  eigenproblem at all and the guard is a pure filter. **B vs C:** B's 16 are
  exactly C's first 16 (worst residual 0.000e+00, no index skipped) — no
  duplicate fired at res 60, itself a new data point given the res-40 firings.
  **Cost:** A 934.0 s, B 2948.5 s, C 2997.7 s ⇒ `B/A = 3.16×`, `C/A = 3.21×`;
  the two ratios agreeing confirms the expense is the eigenproblem, not the
  filtering. Item 13 closed; folded into `fdtd_vs_EME_mode_validation.md` item
  13, which also records why `_DUP_MODE_MAX_EXTRA` was **not** cut yet.

### 11586986 — `beef_dispersive_convergence`
- **Script:** `batch/verification/run_dispersive_convergence.batch` →
  `examples/verification/python/verify_item_8_dispersive_convergence.py`
- **Resources:** 1 node, 8 cpus, serial, 480 min budget. Eigensolve only.
- **Output:** `reports/11586986.out`.
- **What it does:** Sweeps res 20→120 on the coupler body cross-section under
  both material settings, reporting Δn_eff of the even/odd supermode pair and
  the *step* between consecutive resolutions.
- **Why:** `material_model_validation.md` work item 4 asks for a re-run of the
  full resolution sweeps, because smoothing moves the mesh-convergence point
  and item 11's res-80 pin is not inherited once dispersion is on. That sweep
  is 11 array tasks × 2 nodes × 24 ranks at up to 24 h. But the pin was
  justified on **Δβ**, a mode-solve quantity — so this answers "does the
  convergence point move" on one serial task first, to decide whether the
  node-days are warranted and over what range.
- ~~**Caveat, recorded 2026-07-31:** submitted *before* λ₀ became band-derived
  (commit `4296ee7`), so its **dispersive column is stale** — it used the old
  fixed λ₀ = 2.0.~~ **Wrong diagnosis, corrected 11:45.** λ₀ cannot matter here,
  because the dispersive medium was never built at all.
- **⚠️ This job measured nothing.** Its "dispersive" column is **byte-identical**
  to its non-dispersive one — same n_eff, Δn_eff and Δβ to six decimals at res
  20/30/40/60/80. `FDTDSolver._build_section_mode_sim()` built media with
  `si_medium(wl_center)` directly, bypassing `_media()`, so it never read
  `self.dispersive`; `section_mode_neff()` takes a whole band and loops it on
  that one frozen-centre structure. Confirmed three ways: this log, a 2.7 s
  direct check (max |Δn_eff| = 0.000e+00 with `resolve_dispersive → True`), and
  the media themselves (ε 12.0630 vs 12.0946 at 1.60 µm — the mode solve sees
  the frozen one under both settings).
- **Fixed the same day**; a real solve now moves 5.29e-03 across the band where
  it moved by exactly 0. Full reading in `material_model_validation.md` work
  item 4, which is back to **not started**. The non-dispersive column (res
  20→120, all seven rows) is valid and is preserved in that entry, so nothing is
  lost by superseding this job.

### 11589770 — `beef_mmi_1x2_hybrid_calibrated`
- **Script:** `batch/run_mmi_1x2_hybrid_calibrated.batch` (new) →
  `MMI1X2_CALIBRATE_BODY=1 make mmi_1x2_hybrid`.
- **Resources:** 2 nodes × 24 ranks, 300 min — same shape as the uncalibrated
  hybrid array task (11585845_2, 1:46:59), since `calibrate_body_neff` changes
  how the body section's phase is propagated, not how many FDTD sections run.
- **Why:** `mmi_1x2_validation.md` gap 5 measured the hybrid 9.11% **under**
  the full-device FDTD reference, every ΔT negative — the Δβ-discretization
  signature `resolved/directional_coupler_validation.md` items 6/14 describe
  on the coupler, where `calibrate_body_neff=True` took it from ~5% to 0.11%.
  Never tried on this device. `mmi_1x2_hybrid.py` gained a
  `MMI1X2_CALIBRATE_BODY` env flag (redirects output to
  `mmi_1x2_hybrid_calibrated/` so the uncalibrated baseline survives; same
  pattern as `directional_coupler_hybrid_stub_sweep.py`'s
  `DC_CALIBRATE_BODY`).

### 11589779 — `beef_mmi_1x2_gap5_calibrated`
- **Script:** `batch/verification/run_verify_mmi_1x2_calibrated.batch` (new) →
  `MMI1X2_HYBRID_CALIBRATED=1 python verify_mmi_1x2_hybrid_vs_fdtd.py`.
- **Resources:** 1 node, 1 cpu, 20 min. `--dependency=afterany:11589770`.
- **What it does:** Same comparison as 11585850, pointed at the calibrated
  hybrid's SMatrix instead. `verify_mmi_1x2_hybrid_vs_fdtd.py` gained a
  `MMI1X2_HYBRID_CALIBRATED` env override for this; `FDTD_PKL` (the shared
  full-device reference) is deliberately **not** redirected.
- **Why:** answers whether `calibrate_body_neff` closes this device's 9.11%
  the way it closed the coupler's ~5%. Compare its worst max|ΔT| against
  11585850's 0.0359 once both have run.

### 11589771 / 11589899 / 11589936 — `beef_modal_mux_fdtd` — ❌ BLOCKED, do not resubmit
- **Script:** `batch/run_modal_multiplexer_fdtd.batch` (new) → `make
  modal_multiplexer_fdtd` → the new (blocked) `modal_multiplexer_fdtd.py`.
- **Why attempted:** `modal_multiplexer_validation.md` gap 3 — this device is
  the only one in the repo with no full-device FDTD reference at all. The
  script mirrors `mmi_2x2_fdtd.py`'s shape: the same `CompoundGeometry` the
  hybrid script builds, one single-domain `FDTDSolver.solve()` per input mode
  (TE0, TE1) instead of five cascaded sections. Resolution deliberately
  reduced to **20** (env-overridable `MUX_FDTD_RES`) per gap 3's own
  recommendation to scope cheap first — this device is the longest
  single-domain run in the repo (62.16 µm) and 3-D FDTD cost is steep in
  resolution.
- **Attempt 1 (11589771, 22 s):** crashed in `plot_fdtd_geometry()` →
  `CompoundGeometry.source_aperture_at()` → `_locate()`, which raised for any
  x outside `[-total_length/2, total_length/2]`. `FDTDSolver.solve()`'s own
  source-position query (`src_x = -Sx/2 + dpml + source_margin`, landing in
  the solver's external PML-clearance stub) sits at exactly this kind of x, so
  the real simulation would have hit the identical crash had the plot not
  caught it first — this was going to fail regardless, just later and after a
  full node allocation.
  - **Real bug, fixed:** every leaf geometry (`MMI.width_at()`, etc.) already
    clamps to its terminal value for arbitrarily large overshoot past its own
    taper — that's why plain-MMI full-device scripts (`mmi_1x2_fdtd.py`,
    `mmi_2x2_fdtd.py`) never hit this. `CompoundGeometry._locate()` was the
    only one that hard-failed, purely because delegating to a segment needs a
    segment to exist "out there". `beef/geometry/compound.py`:
    `_locate(x, tol, clamp=False)` gained an opt-in `clamp=True` that
    delegates out-of-bounds x to the nearest end segment instead of raising —
    wired into `source_aperture_at()` only; `width_at()`/`center_at()`/
    `section_prisms()`/etc. keep raising, since an out-of-bounds cut *there*
    is a genuine caller error building physical geometry, not a solver's
    legitimate external-stub query. 5 new tests,
    `tests/test_geometry.py::TestCompoundGeometrySourceApertureClamps`,
    mutation-checked (three fail against the pre-fix code, two pass either
    way). Fast suite 482 passed, 1 skipped.
- **Attempt 2 (11589899, 15 s):** cleared the plot crash, hit a second one
  inside `solver.solve()`'s cache-fingerprinting path:
  `geometry.polygons(stub_length=...)`, and `CompoundGeometry.polygons
  (stub_length>0)` raises `NotImplementedError` by design (its own docstring:
  "no generic way to know per-port stub width for an arbitrary segment
  type"). **Misdiagnosed at the time** as a caching-only issue — dropping
  `save_dir` (disabling the cache) does skip *that* call site.
- **Attempt 3 (11589936, 13 s), same script with `save_dir` dropped:** still
  crashed, same `NotImplementedError`, from a **second, unconditional** call
  site: `solve()` → `to_meep_prisms()` → `geometry.polygons(stub_length=...)`,
  which runs regardless of caching — it is how the geometry is actually built
  for the sim. **This is the real blocker.** `FDTDSolver.solve()` (the full
  single-domain path) cannot run on a `CompoundGeometry` top-level device at
  all today; only `solve_section()` (what `HybridSolver` always uses,
  including this same device's own hybrid script) avoids it, via
  `section_gds()`/`to_gds_eme()` instead of `polygons()`.
- **Not fixed.** Unlike the aperture clamp, this needs a real design decision
  (which segment's stub-port width wins at an ambiguous end — e.g. a
  multi-port `WaveguideRun`), not a mechanical patch, and the author who wrote
  `polygons()`'s sibling stub-handling explicitly declined to guess at it.
  Same caution bar applied to item 14 and the k-seed guard fix this session:
  not doing production-adjacent, under-verified geometry changes blind. Filed
  as `notes.md` Geometry item 8. `modal_multiplexer_fdtd.py` and its batch
  script are both marked BLOCKED at the top — do not resubmit until this is
  actually implemented.
- **Cost:** three attempts, ~15 s / minimal node-seconds each. No wasted
  budget, no partial/corrupt output (only `geometry.png` was ever written;
  crashes happen before any `smat.save()`).

## Not listed here

An `elevator_*` array job (11585814 at the 07-31 snapshot; the 07-30
snapshot's 11568634 / 11569468 have since cleared) is also queued, under
`/storage/scratch1/8/jhiesener3/esc/elevator_sim_tower_from_gds/` — a
different project entirely, unrelated to BEEF or anything in this
directory. Omitted above; see that project's own tree if they matter.

A second set of these same three device jobs was briefly submitted on
2026-07-31 08:26 to `cpu-medium` at `--mem-per-cpu=16G` (11585878/79/80,
`sbatch` overrides only, no script changes) and **cancelled before any of
them started** — no output, no billing, nothing to fold in. Recorded only so
the gap in the job-ID sequence isn't mistaken for a lost run.
