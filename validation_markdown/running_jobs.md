# Active SLURM jobs — snapshot

**Snapshot taken 2026-07-31, 08:22 EDT** (updated from the 2026-07-30 21:46
snapshot — see re-queue note below). This is a point-in-time list, not a
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
| ~~11577985~~ → **11585844** `[1-2]` | `beef_mmi_2x2_examples` | 🟢 **RUNNING** (task `_1` started 08:24; `gedc_main`) | [mmi_2x2_validation.md](mmi_2x2_validation.md) gap 5 (data) |
| ~~11577986~~ → **11585849** | `beef_mmi_2x2_gap5` | PENDING (dep on 11585844) | [mmi_2x2_validation.md](mmi_2x2_validation.md) gap 5 (comparison) |
| ~~11578018~~ → **11585845** `[1-2]` | `beef_mmi_1x2_examples` | PENDING (`gedc_main`; `ReqNodeNotAvail`, not a billing block) | [mmi_1x2_validation.md](mmi_1x2_validation.md) gap 5 (data) |
| ~~11578019~~ → **11585850** | `beef_mmi_1x2_gap5` | PENDING (dep on 11585845) | [mmi_1x2_validation.md](mmi_1x2_validation.md) gap 5 (comparison) |
| ~~11578026~~ → **11585847** `[0]` | `beef_modal_mux_examples` | PENDING (`gedc_main`; `ReqNodeNotAvail`, not a billing block) | [modal_multiplexer_validation.md](modal_multiplexer_validation.md) gap 6 (exercise run) |

As of 08:24 EDT, `11585844_1` (2×2 MMI) has **started running**; the other
four are still queued. The detail sections below carry the current job IDs —
only the ID and billing account changed in the re-queue, so each entry's
script, resources and rationale are unchanged from the 07-30 submission.

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

### 11585850 — `beef_mmi_1x2_gap5`
- **Script:** `batch/verification/run_verify_mmi_1x2_hybrid_vs_fdtd.batch` →
  `examples/verification/python/verify_mmi_1x2_hybrid_vs_fdtd.py`
- Same shape as 11585849, ported to `mmi_1x2`'s port names (`out_0`/`out_1`
  vs `out_top`/`out_bot`). `--dependency=afterany:11585845`.
- **Why:** `mmi_1x2_validation.md` gap 5.

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
