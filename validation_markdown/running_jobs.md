# Active SLURM jobs — snapshot

**Snapshot taken 2026-07-31, 10:35 EDT** (fourth snapshot of the day: 07-30
21:46 → 07-31 08:22 re-queue → 08:46 → now). This is a point-in-time list, not a
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
| ~~11577985~~ → **11585844** `[1-2]` | `beef_mmi_2x2_examples` | `_1` ✅ **COMPLETED** 1:16:45; `_2` 🟢 **RUNNING** | [mmi_2x2_validation.md](mmi_2x2_validation.md) gap 5 (data) |
| ~~11577986~~ → **11585849** | `beef_mmi_2x2_gap5` | PENDING (dep on 11585844) | [mmi_2x2_validation.md](mmi_2x2_validation.md) gap 5 (comparison) |
| ~~11578018~~ → **11585845** `[1-2]` | `beef_mmi_1x2_examples` | ✅ **COMPLETED** — `_1` 8:57, `_2` 1:46:59; both SMatrix files on disk | [mmi_1x2_validation.md](mmi_1x2_validation.md) gap 5 (data) |
| ~~11578019~~ → **11585850** | `beef_mmi_1x2_gap5` | ✅ **COMPLETED** 10:28, 0:09 — **max\|ΔT\| 0.0359 (9.11%), results folded in** | [mmi_1x2_validation.md](mmi_1x2_validation.md) gap 5 (comparison) |
| ~~11578026~~ → **11585847** `[0]` | `beef_modal_mux_examples` | 🟢 **RUNNING** (08:38, 4 nodes) | [modal_multiplexer_validation.md](modal_multiplexer_validation.md) gap 6 (exercise run) |
| **11586046** | `beef_diagnose_items_11_12` | 🟢 **RUNNING** (08:42) — new, submitted 2026-07-31 | [fdtd_vs_EME_mode_validation.md](fdtd_vs_EME_mode_validation.md) items 11 + 12 |
| **11586098** | `beef_diagnose_item_13` | 🟢 **RUNNING** (08:56) — new, submitted 2026-07-31 | [fdtd_vs_EME_mode_validation.md](fdtd_vs_EME_mode_validation.md) item 13 |
| **11586986** | `beef_dispersive_convergence` | 🟢 **RUNNING** (10:25) — new, submitted 2026-07-31 | [material_model_validation.md](material_model_validation.md) work item 4 (cheap half) |

**First real results are in.** The 1×2 MMI's whole gap-5 chain has completed
end to end (both data tasks plus the comparison) and its numbers are folded
into `mmi_1x2_validation.md`. `11585849` (2×2 comparison) is still pending only
on `11585844_2`, which is correct. Everything else is running.

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

### 11585844 `[1-2]` — `beef_mmi_2x2_examples` (array indices 1,2 only)
- **Script:** `batch/run_mmi_2x2_examples.batch` → `make mmi_2x2_fdtd
  mmi_2x2_hybrid`
- **Resources:** 2 nodes × 24 ranks per array task, 5h limit each.
- **What it does:** Runs the 2×2 MMI's full-device FDTD reference
  (`mmi_2x2_fdtd.py`) and the FDTD/EME hybrid cascade
  (`mmi_2x2_hybrid.py`), each saving an `SMatrix` pickle. Array index 0
  (`mmi_2x2_eme`) deliberately skipped — not needed for gap 5.
- **Why:** gap 5 ("no hybrid-vs-full-device comparison on record") needs
  both scripts' saved output to exist before anything can compare them.

### 11585849 — `beef_mmi_2x2_gap5`
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

### 11586046 — `beef_diagnose_items_11_12`
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

### 11586098 — `beef_diagnose_item_13`
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
- **Caveat, recorded 2026-07-31:** submitted *before* λ₀ became band-derived
  (commit `4296ee7`), so its **dispersive column is stale** — it used the old
  fixed λ₀ = 2.0. The non-dispersive column and the qualitative
  "does-the-point-move" reading are unaffected, but the dispersive half needs
  re-running before any number here is recorded as current.

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
