# Active SLURM jobs — snapshot

> **THREE independent handoffs are live in this file, written by different
> sessions on 2026-09-14. Read the one for your device:**
> * **1×4 SiN300** — immediately below. One array job RUNNING (13169522).
> * **4×4 90° hybrid** — under the THIRTY-FIRST SNAPSHOT further down.
>   **One job RUNNING (13197102)**, the reopened knob-4 convergence curve;
>   a timeout on it is expected, not a failure. Campaign otherwise landed.
> * **2×2, straight + waisted** — immediately after this block. The two
>   shipped examples are DONE and verified; the ALT cascade (FDTD ramps)
>   array 13200430 COMPLETED but its result is **INVALID** (non-passive ramp matrix) — see below.
>
> Each says "no other job of mine is running" — true from that session's own
> scope, not a statement about the whole queue.

> ## ⛔ 2×2 ALT CASCADE — 2026-09-14, array 13200430 COMPLETED but the
> ## RESULT IS INVALID. Do not quote its numbers.
>
> All 11 tasks COMPLETED cleanly (~3 h each, 110884 s of solve summed) and
> the merge produced a full comparison. **The comparison is not a
> measurement**: the FDTD ramp section's transfer matrix is numerically
> broken, and the cascade output still looks plausible, which is exactly what
> makes it dangerous.
>
> **The evidence, from the cached ramp matrix at band centre** (10 modes,
> section x=[-9.0, -2.0]; `sum|M|^2 = 8.8714`):
>
> | band | power out when excited | power in when measured | |
> |---|---|---|---|
> | 1-5 | 0.97-1.00 | 0.99-1.00 | ok |
> | **6** | **0.0095** | 0.9529 | flux-calibrated fallback |
> | **7** | **0.0146** | 0.9516 | flux-calibrated fallback |
> | 8 | 0.9532 | 0.9953 | ok |
> | **9** | **1.9505** | 0.9952 | ⛔ **>1: NON-PASSIVE** |
> | **10** | 0.9949 | **0.0036** | flux-calibrated fallback |
>
> Exciting band 9 returns ~2x the input power, which is impossible for a
> passive section. Bands 6 and 7 lose essentially all of it. Band 10 receives
> essentially nothing.
>
> **These are NOT negligible modes.** The three flagged bands carry **11.5 %
> of the excited power and 21.5 % of the measured power**. The usual excuse
> for a large relative mismatch — that it is a ~1e-6 mode, as recorded in
> `mode-signs-gate-on-strong-modes` — does not apply.
>
> **Where it comes from.** Every task logged the same warning:
> `solve_section: excite band N normalization mismatch between
> get_eigenmode_coefficients() and get_fluxes() (max 100.0%) for section
> x=[-9.0000, -2.0000] — likely a near-degenerate excited mode; falling back
> to flux-calibrated normalization`. 792 instances / 24 ranks = 3 per task,
> on bands 6, 7 and 10, in every one of the 11 wavelengths. (Divide per-rank
> log counts by the rank count before reading them — see
> `divide-log-counts-by-mpi-rank-count`.) The flux-calibrated fallback
> recovers a magnitude but evidently not a consistent one here.
>
> **Suspected cause, NOT confirmed.** The ramp spans 2.80 µm to 2.24 µm. The
> 2.80 µm face supports ~14 modes and the 2.24 µm waist exactly 10, so
> requesting 10 at BOTH ends asks the face for modes 6-10 in a region where
> they are near-degenerate with each other — precisely the condition the
> warning names. `N_FACE` is already a separate knob from `N_MODES` in the
> script for this reason; it was simply left equal.
>
> **What to do next, in order:**
> 1. Diagnose at ONE wavelength, not 11 — the ramp section alone, sweeping
>    `MMI2X2T_ALT_NFACE` (try 8, 12, 14) and checking column power for
>    passivity. Cheap relative to the array.
> 2. Only once every column sums to <= 1 is the alt-vs-sibling comparison
>    worth re-running.
> 3. Do NOT sweep FDTD_STUB first. It was the pre-registered suspect, but
>    passivity violation is a mode-basis failure, not a geometry one.
>
> **For the record, the numbers the invalid run produced** (so nobody
> re-derives them thinking they are new): alt vs sibling vs full FDTD —
> IL 0.3500 / 0.2446 / 0.2010; imbalance 0.0317 / 0.0259 / 0.0124;
> phase -90.3345 / -88.5220 / -89.4822; max|dT| bar 0.0364 / 0.0087;
> max|dT| cross 0.0106 / 0.0201. The cross-port improvement is the only
> thing that would have been interesting, and it is exactly the kind of
> thing a broken normalization can manufacture.
>
> *(original block, kept for provenance)*
>
> ## 🏃 2×2 ALT CASCADE — 2026-09-14, array 13200430 RUNNING (11 tasks)
>
> **`batch/run_mmi_2x2_tapered_alt.batch`**, array 0-10, one WAVELENGTH per
> task, 24 ranks each, `-t 1440`, `gts-sr108-gedc_main`. All 11 started
> immediately. **13198943 was an earlier submission of the same thing at
> `RAMP_DECAY_TOL=1e-05`; it was CANCELLED 30 min in (no shards written) once
> the tolerance sweep came back — see the decay_tol block below.** Merge when the array is COMPLETED:
> `make -C BEEF/examples mmi_2x2_tapered_hybrid_alt_merge`.
>
> **What it is:** `examples/mmi_2x2_tapered/python/mmi_2x2_hybrid_alt.py` —
> the same waisted 2×2 at the same design point, cascaded
> `[fdtd, eme, fdtd, eme, fdtd, eme, fdtd]` instead of `[fdtd, eme×5, fdtd]`.
> The two width RAMPS are solved by FDTD rather than an EME staircase;
> everything else is held fixed, so the difference against the sibling
> `mmi_2x2_hybrid` IS the measurement. Motivation: the staircase is the
> leading suspect for the sibling's method error (its IL overshoots
> full-device FDTD and CLIMBS with cell count, so C_TRANS is an optimum
> rather than a floor — not how a converging discretisation behaves).
>
> ⚠ **~11 FDTD runs per wavelength, against the sibling's 1.**
> `solve_section()` needs one run per excited mode and a ramp carries
> N_MODES on both faces, so the ramp alone is 10 runs. Ramp 2 is solved;
> ramp 4 is a plain `mirror_pairs={4: 2}` (legal — both are on-axis body
> sections) and the output port stays `MirrorSpec`-derived, so neither is
> simulated. 24 ranks is deliberate and is the OPPOSITE of the sibling's
> resourcing: this job is FDTD-stepping-bound (~98 % scaling), not
> MPB-bound (rank-redundant).
>
> ⚠ **`decay_tol`: the ramp keeps 1e-11, and the FIRST ANSWER HERE WAS WRONG.**
> A stalled smoke run was attributed to this device's known decay_tol hard
> floor (1e-11 unreachable on an 84 µm 4×4: floors, then diverges) and the
> ramps were dropped to 1e-05 on that reasoning. The sweep then measured the
> opposite, on the exact ramp section [-9.0, -2.0] (res 12, 2 modes, 1 wl),
> diffing the transfer matrix between consecutive tolerances:
>
> | decay_tol | seconds | max\|dM\| vs previous |
> |---|---|---|
> | 1e-04 | 424 | — |
> | 1e-05 | 632 | 8.4e-04 |
> | 1e-06 | 869 | 8.5e-06 |
> | 1e-07 | 436 | 5.7e-07 |
> | 1e-09 | 476 | 7.3e-07 |
> | **1e-11** | **650** | **3.8e-10** |
>
> 1e-11 IS reachable on a 7 µm section — nowhere near the 84 µm regime — and
> **cost is FLAT in the tolerance** (424 s to 869 s, non-monotonic), because
> the stop is dominated by the field physically leaving the section rather
> than by the threshold. A looser tolerance buys nothing and costs 8.4e-04.
> Both ports and ramps therefore run at 1e-11; `RAMP_DECAY_TOL` stays a
> separate knob because ramp and port are different problems, NOT because
> this one needed a different value. The table is reproduced in the script.
>
> ⚠ **UNRESOLVED — near-degenerate normalization on the ramp.** The sweep
> logged 12 instances of
> `solve_section: excite band N normalization mismatch between
> get_eigenmode_coefficients() and get_fluxes() (max 42.0%) ... falling back
> to flux-calibrated normalization`, on bands 1 and 2 of section
> x=[-9.0, -2.0]. That is the near-degenerate path from the 2026-08-09
> shared-library fix. It fired at res 12 with 2 modes; whether it fires at
> res 40 with 10 modes, and how much it moves the answer, is NOT known. It is
> a plausible confounder for the whole alt-vs-sibling comparison, so grep the
> task logs for "normalization mismatch" BEFORE reading the merged numbers at
> face value.
>
> ⚠ **A concurrent session edited `beef/solvers/fdtd.py` four minutes before
> this was submitted** (phase columns added to the `snapshot_every` CSV
> writer). Judged inert for this job — it touches `solve()`'s snapshot path,
> not `solve_section()`, which is all this script uses, and snapshot args
> are explicitly outside the cache fingerprint. The live tree imported clean
> at submit time. But a PENDING job imports the library when it STARTS, so
> if these tasks fail on an import or an attribute error, suspect that file
> first.
>
> **On completion:** `make -C BEEF/examples mmi_2x2_tapered_hybrid_alt_merge`
> (seconds, no MPI). It refuses a partial merge.
>
> ⚠ **CHECK `squeue` FOR AN EXISTING `beef_mmi_2x2_tapered_alt` BEFORE
> SUBMITTING THIS ARRAY.** Two sessions have now each submitted it
> independently (a duplicate task 5, job 13198990, was submitted by one and
> cancelled). Two jobs writing the same shard path is wasted compute rather
> than corruption — the shard write is temp-file plus atomic rename — but it
> is still two jobs.
>
> *(A second session wrote its own copy of this block while the array was
> being set up. The copies were consolidated here on 2026-09-14; the only
> substantive difference was the decay_tol paragraph above, which that copy
> still carried in its pre-sweep form. Nothing else was dropped.)*
>
> **Open:** the FDTD_STUB = 1.0 µm choice is defensible (seams land in
> uniform material with dpml clearance, four bounds checked, sections tile
> the device to 1e-9) but it was NOT swept. If the alt result disagrees with
> full-device FDTD, sweep `MMI2X2T_ALT_STUB` before blaming the cascade.


> ## ✅ 2×2 (straight + WAISTED, shipped examples) — 2026-09-14, DONE.
>
> **BOTH JOBS COMPLETED, and BOTH reproduced their previous output exactly.**
> `13189623_1` COMPLETED in 01:52:56, `13189859_2` in 00:05:03. Every file
> byte-identical to the pre-run copy; `max|delta| = 0.000e+00` on T_bar and
> T_cross across all 51 FDTD wavelengths and across the hybrid's 201 × 11
> grid. Band centre, before and after: `T_bar 0.455993`, `T_cross 0.520670`,
> `total 0.976664`.
>
> **What the FDTD run actually bought, confirmed by probe:** the same
> `cache_only=True` check that read MISS before it now reads HIT. So the
> prediction from the fingerprint diff held — the two keys added on
> 2026-09-12 are inert here, and the 48-core re-run purchased a cache entry
> reachable under the new fingerprint, not a correction. Anyone re-running
> `examples/mmi_2x2/python/mmi_2x2_fdtd.py` today gets it from cache.
>
> ⚠ **Other devices whose stored `FDTDSolver.solve()` entries predate
> 2026-09-12 are still in the pre-re-run state** — they will miss and
> resimulate on next use. This was only fixed for the straight 2×2 because
> that is what was asked for; it is not a repo-wide sweep.
>
> ---
>
> *(as submitted)*
>
> **Job 13189623_1** — `sbatch --array=1 --account=gts-sr108-gedc_main
> batch/run_mmi_2x2_examples.batch`. Array index 1 is `mmi_2x2_fdtd` in that
> script's own `EXAMPLES=(mmi_2x2_eme mmi_2x2_fdtd mmi_2x2_hybrid)`, so this
> re-runs ONLY the STRAIGHT device's whole-device FDTD. 2 nodes × 24, `-t 300`.
>
> ⚠ **Account overridden to `gts-sr108-gedc_main`.** The batch file hard-codes
> `gts-sr108-ece`, which is currently out of balance — its jobs sit at
> `AssocGrpBillingMinutes`. Check `sshare` before reusing that file; do not
> assume the hard-coded account is spendable.
>
> ⚠ **That script `rm -rf`s `example_output/mmi_2x2/mmi_2x2_fdtd` before
> running.** The pre-run copy was saved outside the repo first. Baseline to
> compare the re-run against, at band centre: `T_bar 0.455993`,
> `T_cross 0.520670`, `total 0.976664`.
>
> **Why it is being re-run at all, and why this is NOT a correctness worry.**
> Its cached entry is a MEASURED cache miss: `beef/cache/config.py` gained
> `geometry_mode` and `design_grid_resolution_factor` in the
> `("FDTDSolver","solve")` allowlist on 2026-09-12, and the stored entry
> predates them. A field-by-field diff of the fingerprint the current code
> builds against the stored one shows those two keys are the ONLY difference —
> geometry polygons, call kwargs and every other solver attribute are
> identical, and both new keys are inert here (`geometry_mode="prism"` is the
> old behaviour). So the numbers are expected to come back unchanged; the run
> buys a cache entry reachable under the new fingerprint, not a correction.
>
> **Job 13189859_2** — the same batch file at array index 2, i.e. the STRAIGHT
> `mmi_2x2_hybrid`, same account override.
>
> Its CASCADE needed no job and is already confirmed: re-run on the login node
> 2026-09-14, it reproduced the Aug-31 output BIT-IDENTICALLY — byte-identical
> `smatrix_hybrid.pkl`, `hybrid_results.png`, `mm_length_sweep.{npz,png}` and
> `geometry_hybrid.png`, with `max|delta| = 0.000e+00` across the whole
> 201 × 11 grid. `solve_section`'s allowlist was not touched by the 2026-09-12
> commit, so its cache still hits.
>
> ⚠ **What that login-node run did NOT cover: step 8**, the
> `plot_hybrid_body_fields()` body-field reconstruction. It was still in that
> step at 30 min and was killed by the wrapper's own timeout (exit 124). No
> damage — it never reached the write, so `fields/hybrid_body_mmi2x2_body_L28.0.png`
> is untouched Aug-31 output — but it is therefore the one part of this device
> NOT re-verified against the current library. That step is a fresh mode solve
> plus a diagnostic plot; it does not feed the S-matrix. 13189859_2 exists to
> close it on a compute node rather than repeat a 30-minute login-node solve.
>
> Note that this job, like 13189623_1, `rm -rf`s its own output subfolder
> first. The pre-run copy was saved outside the repo.
>
> **The WAISTED 2×2 is DONE and needs no job**, despite being the more
> expensive device: both examples reproduce the study's stored numbers
> bit-identically from cache (whole-device FDTD 0.7 s; all 11 hybrid
> wavelengths 109 s against the ~15 h they originally cost). See
> `mmi_2x2_waist_body_validation.md`, "The device now ships as a normal
> example".

> ## 🏃 HANDOFF — ONE JOB CURRENTLY RUNNING, 2026-09-14
>
> **Job 13169522** (array, 4 tasks: `13169522_0..3`), submitted 2026-09-14
> ~00:00, `-t 6-00:00:00`. As of this write-up: tasks 0/1/2 `RUNNING` at
> 12h00m, task 3 still `PENDING`. Check with:
> `sacct -j 13169522 -n -P -X --format=JobID,State,Elapsed,ExitCode`
>
> **What it is:** `beef_prewarm_1x4_cband21_n24` —
> `_prewarm_mmi_1x4_hybrid.py`, one output branch per task, N_MODES=24,
> uncalibrated, `monitor_size` pinned, on a C-band sub-grid (1.530–1.565 µm,
> 21 pts, 1.75 nm spacing). Batch file:
> `batch/verification/run_prewarm_mmi1x4_hybrid_cband21_n24.batch`. Full
> context: `mmi_1x4_sin300_validation.md`'s "✅ 9d RESULT" section, last
> paragraph before "Two operational findings."
>
> **This is attempt 2.** Attempt 1 (job 13092670) TIMEOUT'd at 2d12h with
> **zero cached progress** on any of the 4 branches (confirmed via a direct
> `cache_only=True` probe, not assumed) — no duplication bug found (two
> "Computational cell" log blocks that looked suspicious are legitimately
> the input and output sections, confirmed via differing prism geometry),
> just genuinely more expensive than the 2d12h budget allowed for. Root
> cause of *why* N_MODES=24 costs this much more than 16 did (which
> finished in 14–16h on the coarse 1.50–1.60 grid) was **not** confirmed —
> plausibly more MPB iterations near cutoff at higher mode counts, not
> verified further. **Attempt 2 uses `-t 6-00:00:00` as a real budget
> increase, not a confident estimate — the true cost is still unmeasured.**
>
> **On completion (COMPLETED, all 4 tasks):** cascade the 4 branches'
> SMatrices against the already-completed whole-device FDTD C-band
> reference (`example_output/mmi_1x4/mmi_1x4_fdtd_CBAND21/smatrix.pkl`,
> job 13092669, done) and update the "still running" line in
> `mmi_1x4_sin300_validation.md`'s 9d-RESULT section with the real numbers.
>
> **On failure (FAILED, `os.getcwd()`/meow-import signature):** this would
> be the **4th** occurrence of the same dmesg-confirmed Lustre client error
> across this campaign (previously hit 3 times, 7 different nodes, not one
> bad node). **Escalate to PACE support rather than retrying a 3rd time on
> this job** — that threshold was already crossed once on a sibling job
> (the 16-mode calibrated prewarm, dropped after 2 losses) and a 4th
> instance overall is enough evidence for a ticket, not another blind
> resubmit. If it TIMEOUTs again instead, re-check with a `cache_only=True`
> probe before deciding whether to extend further — do not assume progress
> without checking, both prior assumptions on this exact question turned
> out wrong before direct verification.
>
> No other job of mine is running as of this snapshot — this is the only
> one to track. (Other users'/sessions' jobs sharing this cluster and cache
> — `beef_bo_*`, `beef_mux_*`, `beef_mmi4x4_*`, `beef_mmi2x2_*`, `sbend-*` —
> are unrelated to this device and not mine.)

**★ THIRTY-FIRST SNAPSHOT — 2026-09-14, verified against `squeue -u jhiesener3`.
This is the current one.** Everything below the next rule is the thirtieth
snapshot (2026-09-01), not reverified this pass.

> **Convention:** every snapshot in this file says "This is the current one"
> — each was written when it WAS. Only the topmost is current; the rest are
> a historical log, newest first. Do not read an older snapshot's own claim
> as a live statement.

## 4×4 90° hybrid — ONE JOB RUNNING (13197102), submitted 2026-09-14

> **Supersedes the "nothing running" line below**, which was true when this
> section was written earlier the same day.

**Job 13197102** (`beef_mmi4x4_tolcurve`, 15 h, 24 ranks, 1 node, account
`gts-sr108-gedc_main`). Reopens knob 4: job 13110659 bounded the FDTD
reference's phase noise from **two points** (`decay_tol` 1e-4 vs 1e-5), which
is a delta between consecutive rungs, not a convergence test. This is the
convergence curve — L=15 top-port at `decay_tol=1e-9`, snapshotting every 50
meep time units, which because DFT monitors accumulate yields the whole
phase-vs-`decay_tol` curve from one run.

⚠️ **A TIMEOUT IS THE EXPECTED OUTCOME, NOT A FAILURE.** `1e-9` is set below
what the run can reach deliberately; the stopping point is not what is being
measured. Samples append as they are taken, so a walled run keeps everything it
wrote and the CSV is readable *while the job is in flight*:

    example_output/mmi_4x4_90hybrid/mmi_4x4_90hybrid_fdtd_L15.00_p4_tol1e-09/tolcurve_p4_port4.csv

If `energy_ratio` floors and then GROWS (this device's measured res=20
signature: 5.16e-3 → 51.9 → 2.14e8), that turnaround is itself the second
wanted result — the floor at this domain size, which is unmeasured; the ~5e-6
figure everything assumes was measured on the 84 µm full device, not at L=15.
Discard every sample after the turnaround. **Do not auto-resubmit on timeout —
report and wait.**

Library note: `beef/solvers/fdtd.py`'s snapshot writer now also emits
`ph_out{i}_wl{j}` (`arg(S)`, degrees), appended after the `T` columns. Both
tolcurve plotters select by name, so old and new CSVs stay mutually readable.
Verified against the returned `SMatrix` to <0.01°. Fires only when
`snapshot_every` is set, so it cannot reach any other queued job.

### The rest of the campaign — all COMPLETED, none failed

| job | what | elapsed | outcome |
|---|---|---|---|
| 13102582–589 | calibration probes, steps 3/4/5 | 0:40–3:41 | grid dispersion confirmed; group index + bandwidth ruled out |
| 13110658 (×16) | 16-mode calibration at **res = 60** | 1:02–2:44 | β fit 44× better, agreement 18× WORSE — resolutions must match |
| 13110659 | L=15 FDTD at **decay_tol = 1e-5** | 4:00:49 | reference phase noise 0.029° — not the limit |
| 13152193 (×4) | input sections at **24-mode** basis, L=15 | 2:43:52 | power leak 0.7% → 0.25% |
| 13152194 (×24) | **24-mode** calibration at res = 40 | 2:33–2:45 | mode 24 offset is an outlier (+0.0796) |
| 13162199 | 24-mode cascade vs FDTD, L=15 | 4:04:01 | **null** — out8 unchanged at 2.2° |

Results and full reasoning: `mmi_4x4_90_hybrid_validation.md`, the dated
entries at the end of the file, and the CURRENT STATE banner at its top.

### ⚠️ Two zombie jobs in the queue — will NEVER run, safe to cancel

`13092912` (`beef_mmi4x4_calib2cmp`, PENDING/Dependency) and `13092882`
(`beef_mmi4x4_calib2agg`, PENDING/DependencyNeverSatisfied) are left over from
job **13092881**, the 16-task calibration array that FAILED on the duplicate-
eigenmode guard before that guard was fixed. Their `afterok` dependency can
never be satisfied. The work they were waiting on was redone and completed by
jobs 13098452/453, so nothing is lost by cancelling them. They have been left
in place rather than cancelled unilaterally — they cost nothing while pending.

### Not this device — but note the billing block

`13169522` (`beef_prewarm_1x4_cband21_n24`, 6-day limit) is the 1×4 C-band
prewarm from **concurrent session activity, not this campaign** — recorded
only so it is not mistaken for it. State as of 2026-09-14: tasks 0/1/2
RUNNING (~12:00 in), **task 3 PENDING on `AssocGrpBillingMinutes`**. That is
the recurring account-quota block on `gts-sr108-ece`, not a fault in the job;
the usual remedy is resubmitting the blocked task under a different account
(`sshare` first — headroom shifts, do not hardcode which account has it).
**Whoever owns that 1×4 session should decide — do not move another session's
task without asking.**

---

## 🤝 HANDOFF — 4×4 90° hybrid, state as of 2026-09-14

**Nothing of this device is running, and nothing needs babysitting.** A new
session can pick this up cold. Verified against `squeue` and a `pgrep` for
stray login-node processes: both clean.

**Where the device stands:** it MEETS the CMRR spec — `mm_length = 158.80 µm`,
min CMRR 20.326 dBe (≥ 20), max |Δφ| 3.664° (≤ 5), excess loss 2.444 dB. That
is the calibrated re-score; it supersedes every earlier CMRR number in
`mmi_4x4_90_hybrid_validation.md`. Read that file's **CURRENT STATE banner at
the top** first — the file is ~6000 chronological lines and most of it is
superseded.

**What is closed, by measurement, so nobody re-runs it:**

| closed | verdict |
|---|---|
| calibration at higher resolution (knob 1) | ❌ NEGATIVE — a res=60 calibration against a res=40 reference is 18× worse. Resolutions must MATCH. |
| interface mode count 16 → 24 (knob 3) | ❌ NULL — recovered ⅔ of the leaked power, dominant arm moved 0.0° |
| FDTD reference `decay_tol` (knob 4) | ✅ NOT the limit — reference phase noise is 0.029° |
| group index / pulse bandwidth | ❌ ruled out — offset identical to 6 decimals across bw 1.5/2.0/4.0 |
| `waveguide_offset` | ❌ closed earlier — W/4 is near-optimal independent of ΔL |

**The one live thread, and the cheapest next move.** Power-weighted
hybrid-vs-FDTD phase agreement is 2.12°; the target (1.5% of the 90°
quadrature spacing) is 1.35°. The residual is a systematic Δn_eff ≈ 3.0e-4,
and the leading candidate is now the **number and span of the CALIBRATION
LENGTHS** — currently three, spaced by `max_safe_spacing/2` from the
branch-safety formula. The arithmetic: the β-fit residual at res=40 is
5.58e-3, which over the fit's own 0.84–2.53 µm baseline implies Δn ≈ 4.7e-4 —
the same order as the 3.0e-4 observed. A longer/denser baseline determines β
more precisely, and β precision is what the residual is made of. This knob was
NOT on the original list (it is not in the knobs table in the device doc,
which is marked superseded) and has not been tried.

**Do NOT just raise `res` to fix this** — that is knob 1, already closed
negative, unless the full-device FDTD reference moves with it (~5× cost in 3D).

**Two loose ends:**
1. **Mode 24's calibrated offset is +0.0796** against a smooth +0.0457→+0.0486
   trend across modes 1–23 — the familiar "topmost band of N is unreliable"
   signature. It was NOT excluded from `aggregated_n24res40_n24.npz`. The
   24-mode result was a clean null so it probably did not matter, but drop or
   re-measure mode 24 before reusing that file.
2. **The 42 cached campaign points under `wg_ygap_opt/` cannot be re-scored
   as-is** — calibration is per-cross-section and none of them have one.
   Applying the single measured offset as a uniform shift is a plausible
   shortcut but is UNTESTED; validate on two geometries bracketing the current
   `mm_width` before trusting it, and expect the offset to grow for narrower
   bodies (fewer cells across the field).

**Standing directives that still apply:** if a long run crashes or times out,
report and wait — do not auto-resubmit. Never `git checkout`/`stash` the
`beef/` tree; a concurrent agent may be editing shared library files.

---

**★ THIRTIETH SNAPSHOT — 2026-09-01, verified against `squeue -u jhiesener3`
immediately after submitting. This is the current one.** It covers only the
4×4 90° hybrid BO round; two other jobs in the queue at snapshot time
(`beef_mmi2x2_free` 12674569, `beef_1x2_step0_MPI` 12674181) belong to
concurrent session activity and are recorded here only so they are not
mistaken for this round's. Everything below the next rule is the
twenty-ninth snapshot (2026-08-21), not reverified this pass.

## 4×4 90° hybrid — BO round 7, 4 jobs, all RUNNING from submit

The `(mm_width, y_gap)` campaign's round 7, after the optimizer was
reparameterized to `(mm_width, taper_width)` — full reasoning in
`mmi_4x4_90_hybrid_validation.md`, section "The (mm_width, y_gap)
optimization campaign". 4 nodes / 96 ranks each, `gts-sr108-gedc_main`
(`sshare` at submit: `-gedc_main` Infinity, `-ece` 159.40, `gts-sr108`
71.43), `-t 720`. Both bodies (W=6.5, W=7.0) are cached from Phase A, so
~3-3.7 h/point is the expectation, not the 4:43:45 cold cost.

| job | state | what |
|---|---|---|
| **12677821** | 🏃 RUNNING (started immediately) | `beef_bo_W7.00_y0.65` — taper_width 1.100, the ridge point at W=7.0. **Decides whether the monotone 5.0 → 6.5 objective improvement continues past 6.5.** |
| **12677822** | 🏃 RUNNING | `beef_bo_W7.00_y0.40` — taper 1.350, wide-taper side at W=7.0 |
| **12677823** | 🏃 RUNNING | `beef_bo_W7.00_y0.85` — taper 0.900, narrow-taper side; with the two above and the existing y_gap=1.11 point, W=7.0's interior is bracketed |
| **12677824** | 🏃 RUNNING | `beef_bo_W6.50_y0.44` — taper 1.185, fills the gap at the best width between taper 1.115 (obj 0.677) and 1.305 (obj 0.856) |

All four write `eval_W*_y*.npz` into
`example_output/mmi_4x4_90hybrid/wg_ygap_opt/`. Harvest with
`python mmi_4x4_90hybrid_bo.py status`; refresh the table/plot with
`python mmi_4x4_90hybrid_ygap_collect.py`. **Use
`propose --dry-run` before every new round** — round 6's acquisition
collapse was invisible until a job had already been submitted.

One job was cancelled by hand this pass: **12674307** (`beef_bo_W6.00_y0.52`)
— round 6's single surviving candidate, 0.01-0.02 µm from two points already
measured at objective ≈1.006. Not resubmitted, deliberately.

## 1×4 MMI SiN — hybrid-vs-FDTD campaign, 2026-09-05 → 2026-09-14 (1 job still running)

See `mmi_1x4_sin300_validation.md`'s "✅ 9d RESULT" section for the full
writeup; this is the job ledger.

| job | what | outcome |
|---|---|---|
| 12797780 (×2) | `N_MODES` downward probe, no `monitor_size` pin | TIMEOUT, both tasks — confirmed the pin was missing, motivated the fix below |
| 12824996 (×5) + stage B | `body_stub` cut-position sweep, `measure_margins`, all 4 ports | ✅ COMPLETE — gate passed 0.12% vs shipped; cut position ruled out (0.07% effect) |
| 12853535 (×4) | 16-mode **calibrated** prewarm, attempt 1 | CANCELLED by hand at 1d09h (4/16 calib runs done, projected 5+ days) — switched to uncalibrated |
| 12893611 (×4) | 16-mode **uncalibrated** prewarm, `monitor_size` pinned | ✅ COMPLETE, 14h35m–16h15m/branch — the pin's first real test, worked |
| 12909417 (×5) | `interface_modes` sweep, 16/14/12/10/8 | ✅ COMPLETE, 25–32 s/task — superset slicing verified working, 0.000% vs an independent 12-mode run |
| 12963023 (×4) | 24-mode uncalibrated prewarm | branches 0 and 3 **FAILED** (`os.getcwd()` inside meow's import, dmesg-confirmed Lustre error, same node `atl1-1-02-002-26-2` both times) |
| 13019844, 13020152, 13019542 | resubmits of the 3 failed 24-mode branches above (one more failure mid-sequence, branch 2, on a *different* node — 3 distinct nodes total across the retries) | ✅ all COMPLETE — no compute lost, each failure hit before that branch's own expensive stage |
| 13058948 (×8) | extended `interface_modes` sweep, 24/22/20/18/16/14/12/10 | ✅ COMPLETE, 23–26 s/task — confirmed TE12 at basis index 21 exactly where a two-family fit predicted |
| 12962990 (×4) | 16-mode **calibrated** prewarm, attempt 2 (resumed after the cancel above) | **FAILED**, all 4 branches, 2d21h — same `os.getcwd()`/Lustre signature, 4 *different* nodes. All 16 calibration time-steppings had completed and it was decomposing the last band when it died; `solve_section()` only calls `save_run()` once the full result returns, so **nothing was recoverable**. **Dropped, not retried a 3rd time** — see the doc's operational-findings note |
| 13092669 | whole-device FDTD reference, C-band 1.530–1.565 µm, 21 pt | ✅ COMPLETE, 2d13h — real beat structure in split shape the coarse 11-pt grid was aliasing by luck, no hidden pathology |
| 13092670 (×4) | hybrid 24-mode uncalibrated prewarm, same C-band grid, attempt 1 | TIMEOUT, all 4 branches, 2d12h — initially misread as a duplicate-run bug (it wasn't: two distinct sections, input then output, confirmed via differing prism geometry); direct `cache_only=True` probe confirmed **zero** cached progress on any branch |
| 13169522 (×4) | same, attempt 2, resubmitted at `-t 6-00:00:00` | 🏃 **RUNNING** — true cost at N_MODES=24/N_WL=21 still unmeasured, this is a real budget increase not a confident estimate |

**Two library changes, both additive/default-unchanged, that made the fast
part of this possible:** `MMI1X4_PIN_MONITOR` (default on) on
`mmi_1x4_hybrid.py`/`_prewarm_mmi_1x4_hybrid.py`, pinning
`monitor_center`/`monitor_size` on sections 0 and 2 only — this is what
unlocks `solve_section()`'s mode-count superset cache family; and
`MMI1X4_FDTD_WL_MIN/MAX/NWL` + matching `HPREWARM_` knobs +
`MMI1X4_OUTPUT_SUFFIX` on both scripts, for the C-band sub-sweep without
touching the production 1.50–1.60/11-pt outputs. New verification script:
`mmi_1x4_body_stub_sweep.py`.

---

**★ TWENTY-NINTH SNAPSHOT — 2026-08-21, verified against
`squeue -u jhiesener3` immediately after submitting. This is the current one.
Everything below the next rule is the twenty-eighth snapshot (2026-08-18),
not reverified this pass.**

## C-band measurement-grid re-runs — 4 jobs, all SUBMITTED, none landed

Re-simulating the two 4C10-die "DC_test" couplers (0.5 µm = 1921/1922,
3.0 µm = 1911/1912) and the HCG reflector on the **measurement grid
1530–1565 nm, 1 nm steps, 36 pts**, replacing the 1.50–1.60 µm / 11-pt
(DC) and 21-pt (HCG) production windows. `CALIBRATE_BODY` is on for both
hybrid jobs — that is the point of the re-run, not incidental.

All four on **`gts-sr108-gedc_main`** (unlimited). `pace-quota` at submit
time: `gts-sr108-ece` **$0.31**, `gts-sr108` $6.29, `-gedc_main` Infinity.
The two existing dc_test batch files hardcode `-ece` and would have been
near-unfundable — the new batch files do not reuse them.

| job | state | what |
|---|---|---|
| **12093321** | 🏃 RUNNING (started immediately, 1 node/24) | `beef_dc160_hyb_cband36` — `directional_coupler_soi160_hybrid.py` at the new window. The **3.0 µm calibrated hybrid**, plus its β-calibration plots. `-t 1440`. → `directional_coupler_soi160_hybrid/windows/1.5300-1.5650_n36/` |
| **12093324** | ⏸ PENDING (**Dependency**, `afterok:12093321`) | `beef_dc160_lsweep_cband36` — `..._hybrid_length_sweep.py`. **The 0.5 µm hybrid result.** Held deliberately: it solves the same 2 FDTD sections + β calibration as 12093321, so running it concurrently would duplicate the cold work and put two writers on one set of cache keys. `-t 720`. → `..._length_sweep/windows/1.5300-1.5650_n36/` |
| **12093328_0/_1** | 🏃 RUNNING (both, 2 nodes/24 each) | `beef_dc160_fdtd_cband36` — `..._fdtd_dc_test.py` array over `(0.5, 3.0)`. Single-domain full-FDTD control on both lengths. `-t 720`. → `..._fdtd_dc_test/L{0.500,3.000}/windows/1.5300-1.5650_n36/` |
| **12093330** | ⏸ PENDING (Resources, 2 nodes/24) | `beef_hcg_res80_cband36` — `high_contrast_grating_fdtd.py` at **res=80** on the new window. `-t 1440`. → `high_contrast_grating_fdtd/res_sweep/res80/windows/1.5300-1.5650_n36/` |

**Why 0.5 µm needed the length-sweep script.**
`directional_coupler_soi160_hybrid_dc_test.py` **cannot** produce a 0.5 µm
hybrid result and no parameter tweak fixes it — re-confirmed 2026-08-21 by
executing it, which raises its `_body_stub_cap` assertion. That assert is only
the first wall; the real blocker is the neff calibration. `section_mode_beta()`
needs `body_span - 2·dpml > 0` and then fits β over three sub-lengths spaced
`usable_span/3` apart (`fdtd.py:2320`, `:2410`). At `coupling_length=0.5` the
body is ≤ 0.5 µm total, so even at an unusably thin `dpml=0.15` the fit lever
arm collapses to ~0.07 µm against β ≈ 9.5 rad/µm — the calibration would be
measuring monitor artifacts. **Structural, not a margin.**
`..._hybrid_length_sweep.py` reaches it legitimately: reference geometry at
`L_REF=3.0` (where sections and β are well-conditioned), EME body swept
analytically, floor `SWEEP_MIN = 2·BODY_STUB = 0.5` with the device lengths
injected onto the grid. **Open question for when it lands**: at 0.5 µm the two
S-bends are only 0.5 µm apart, so the independent-scatterer + 2-mode-interface
assumption degrades. 12093328_0 is the direct check. Agreement at 3.0 with
disagreement at 0.5 is the *expected* signature, not a bug.

**Cost note — a narrower band is a LONGER run.** `source_fwidth` scales with
the frequency span `1/λ_min − 1/λ_max`, which drops **2.85×**
(0.041667 → 0.014617), so the Gaussian pulse is ~2.85× longer in time. Every
`-t` above accounts for this. Budgets are margins, not measured bounds; the
precedents for erring generous are job **11726827** (N_WL=601, TIMED OUT at
12 h inside its MPB loop) and `dc_soi160_hybrid` (timed out at 5 h and 4 h
before 12 h worked).

**Cache**: all four are genuinely cold. The wavelength array is hashed into
every FDTD fingerprint (`_resolved_spectral_config`), so nothing in
`example_sim_cache/` is reusable at this window. 12093324 is the one
exception *by construction* — it should be near-all hits off 12093321.

### Landed — 12093328_0, full-FDTD @ 0.5 µm ✅ COMPLETED 00:54:06

First of the four. `-t 720` budget, used **7.5 %** of it — the ~2.85× pulse-
length penalty from the narrower band was real but far smaller than the 2–3 h
estimate (job 11785063_1 was 00:40:42 at 11 pts, so this is ×1.33, not ×2.85;
the ringdown, not the source, dominates this geometry). **Re-budget the
remaining full-FDTD work down** rather than repeating the 12 h ask.

Clean run: no passivity warning, no `check_source_coverage` starvation warning.
`T_total` 0.99864–0.99912 across the band.

**Band-overlap regression — the run is verified, not just finished.** The old
11-pt 1.50–1.60 µm result interpolated onto the new 36-pt grid agrees to:

| | max abs diff | mean abs diff |
|---|---|---|
| `T_through` | 1.2e-04 | 5e-05 |
| `T_cross` | 6e-05 | 2e-05 |
| `T_total` | 1.2e-04 | 6e-05 |

So the narrower source did **not** move the physics — the change is sampling
density only, which is exactly what had to be true before any of these four
runs could be read. This also validates the `windows/` override mechanism
end to end.

**Against measurement**: simulated AIN0/AIN1 split at 1.550 µm = **14.31 dB**
vs the **~14.6 dB** measured for device 1921/1922. Independent corroboration
of the measurement-based id→length assignment in
`directional_coupler_soi160_hybrid_dc_test.py`'s docstring, which had been
flagged there as "a good working hypothesis, not a verified GDS number".
Still not a GDS re-derivation — but the hypothesis now survives a real check.

Values at λ = 1.548 µm: through `T` 0.9637 (−0.16 dB), cross `T` 0.0350
(−14.56 dB), total 0.9987. → `..._fdtd_dc_test/L0.500/windows/1.5300-1.5650_n36/`

**No overwrites**: all six at-risk production files re-checked against their
pre-submit mtimes after this landed — all UNCHANGED.

### Landed — 12093328_1, full-FDTD @ 3.0 µm ✅ COMPLETED 01:00:45

Job C complete (both array tasks). Clean: no passivity or source-coverage
warning, `T_total` 0.9984–0.9988. Band-overlap regression against the old
11-pt run: max abs diff **1.66e-04** (`T_total`), 1.1e-04 (`T_through`),
1.3e-04 (`T_cross`) — same verdict as `_0`, the narrower source did not move
the physics.

Values at λ = 1.548 µm: through 0.9022 (−0.45 dB), cross 0.0962 (−10.17 dB),
total 0.9984, coupling ratio 0.096.

**Against measurement**: split at 1.550 µm = **9.62 dB** vs **~10.0 dB**
measured for 1911/1912. Both dc_test devices now corroborate their
measurement-based length assignment (0.5 µm: 14.31 vs ~14.6 dB).
→ `..._fdtd_dc_test/L3.000/windows/1.5300-1.5650_n36/`

**Full-FDTD re-budget, now on two measured points**: 00:54:06 and 01:00:45
against a `-t 720` ask. Future full-FDTD points on this band are ~1 h jobs.

### 2026-08-22 — ER finding CORRECTED, and an overwrite bug fixed

**Correction to the design finding recorded below.** The earlier claim of a
"peak 33.28 dB at Lc = 26.2 µm" was **over-read**: it was the argmax of a
coarse 40-point Lc axis, and a 400-point rescan shows that structure is
metric ripple, not physics.

**Mechanism, measured**: ER is a MEDIAN over whichever resonances fall inside
the fixed 35 nm band. As coupling_length shifts the comb, dips enter and
leave that window — over the 400-point scan the in-band count moved between
16.2 and 23.6 and **changed 392 times**. The resulting ER jitter is **std
0.223 dB, peak-to-peak 2.561 dB** about a smoothed envelope. The apparent
"global optimum" (32.30 dB at Lc = 26.60) sits only **+1.39 dB above the
envelope** — inside the ripple. It is not a resolvable design point.

**What is actually true** (smoothed envelope, window ≫ the 0.62 µm ripple):

| | |
|---|---|
| envelope peak | **~31.3 dB near Lc ≈ 21.7 µm** |
| within 1 dB | Lc ∈ [17.8, 36.6] µm — **18.8 µm wide** |
| within 2 dB | Lc ∈ [16.2, 38.8] µm |
| within 3 dB | Lc ∈ [14.7, 40.3] µm |
| fabricated Lc = 3.0 µm | 12.35 dB |
| gain from reaching the plateau | **~18.9 dB** |

The ~19 dB headline headroom **survives** the correction — it was never
sensitive to the ripple. What changes is the guidance: this is a **broad,
flat plateau, not a sharp optimum**, so coupling-length fabrication tolerance
is generous (tens of µm within 1 dB), not the ±447 nm the naive peak-width
read suggested. That is better design news, not worse.

Read the ENVELOPE of the ER map, never its argmax. Documented in the script.

**Overwrite bug, found by hitting it.** The sweep's output path encoded band
and all three resolutions — everything that changes what was *simulated* —
but NOT the sweep axes, which change only what was *evaluated*. So the
400×5 rescan silently clobbered the 40×40 production heatmaps at the same
path. Fixed with the same rule used everywhere else in this repo (**the guard
is the env var, not the value**): unset reproduces the documented production
path exactly; setting any sweep-axis knob redirects to
`grids/Lc<..>_n<..>__La<..>_n<..>/`. Both planes now coexist, verified. The
production deliverable was regenerated and is identical to the original
(ER 12.380 / FSR 2.2355 / FWHM 0.0938).

### Landed — 12093330 (HCG res 80) ✅ 03:46:33 and 12117825 (arm β res 80) ✅ 02:32:12

**ALL FIVE JOBS COMPLETE.** Both clean (only meep's grid-rounding notices).
12093330 used 16 % of its 24 h budget, 12117825 21 % of its 12 h — the ~1.3×
narrow-band pulse factor held, confirming the 2.85× projection was wrong in
the safe direction.

**HCG mirror, res 60 → 80 in the overlap band**: R rises from [0.9524, 0.9693]
to [0.9615, 0.9777], mean shift **+0.0088**. A RISE is the expected direction
— `high_contrast_grating_validation.md`'s own sweep climbs monotonically
0.4468 → 0.9697 over res 10 → 60 and is explicitly not converged there. This
adds one point above that; **still not a convergence claim**.

**Arm n_eff, res 60 → 80**: FDTD 2.169312 → 2.172613 (+0.15 %), and the
MPB-to-FDTD excess *shrinks* 1.4627 % → **1.3547 %**. The excess narrowing
with resolution is what you would expect if it is a discretization artifact
rather than a fixed engine offset — a useful datum for the arm-resolution
convergence question that `HCGFP_ARM_RES` exists to answer.

### ★ DELIVERABLE — `hcg_fabry_perot_length_sweep.py` at res 80 on the C band

Run with `HCGFP_CACHE_ONLY=1`, so all four pieces had to hit cache or raise:
**1600 cells in 6 s, zero simulation.** →
`hcg_fabry_perot_length_sweep/1.5300-1.5650_n36/res80_dc40_arm80/`

The 1 nm assembly grid pays off exactly as predicted: per-factor unwrap
margins are **0.026–0.116 rad/step** against the 0.5 rad limit (they were
0.26–1.20 on the old 10 nm grid), and fit residuals ~15× smaller. The
dense-grid adequacy guard clears at **32.6 samples** across the narrowest
resonance in the plane. All 1600 cells carry ≥2 resonances — no NaN holes.

**Fabricated device** (Lc 3.0, L_arm 53.08), res60/old-band → res80/C-band:

| | res 60, 1.50–1.60 | res 80, 1.530–1.565 |
|---|---|---|
| ER | 10.9023 dB | **12.3795 dB** |
| FSR | 2.2455 nm | **2.2355 nm** |
| FWHM | 0.1008 nm | **0.0938 nm** |

All three move coherently with the mirror's R going 0.9693 → 0.9777: higher
R → higher finesse → deeper ER, narrower FWHM, slightly smaller FSR. That
internal consistency is the check that the res-80 rebuild is sound.

**Design findings**
- ER is **essentially 1-D in coupling length**: max std across the whole arm
  axis at fixed Lc is **0.673 dB**. FSR carries the arm dependence (4.11×
  across the arm range). The two design knobs are near-orthogonal — set FSR
  with the arms, ER with the coupler.
- The fabricated coupling length is far off the ER optimum: **11.82 dB**
  (arm-averaged) at Lc = 3.0 vs **31.33 dB** at Lc ≈ 27.0 µm — **+19.5 dB of
  headroom from coupling length alone**. Peak single cell 33.28 dB at
  Lc = 26.2, L_arm = 16.9 µm.
- Caveat on that: the coupler's own S-bend FDTD is res 40, ASSUMED and never
  converged on this platform. The ER optimum sits at Lc ≈ 27 µm, an order of
  magnitude from the 3.0 µm the coupler pieces were validated around. Confirm
  before designing to it.

### Landed — 12093321 (hybrid @ 3.0) ✅ 02:01:32 and 12093324 (length sweep) ✅ 00:00:50

**The DC half of the 2026-08-21 request is COMPLETE**: both requested lengths
(0.5 and 3.0 µm), by both methods, on 1530–1565 nm / 1 nm steps, with the
neff calibration on for both hybrid runs.

12093324 finishing in **50 seconds** is the dependency chain working exactly
as designed — every expensive piece was a cache hit off 12093321, and only
the analytic EME body sweep was new. Had the two run concurrently they would
each have paid the full cold cost.

Job A clean (only meep's benign grid-rounding notices). Band-overlap
regression vs the old 11-pt run: max abs diff **1.67e-04**. Body β
calibration at 1.548 µm: MPB 2.129795/2.101037 → FDTD 2.181470/2.156059,
Δn_eff 0.028758 → 0.025410 — consistent with the pre-existing cached run.
Job B measured the 50:50 length at **12.2316 µm**.

#### ★ METHOD CROSS-CHECK — and a prediction that was WRONG

The expectation recorded above was: *"agreement at 3.0 with disagreement at
0.5 is the expected signature"*, on the grounds that at 0.5 µm the two
S-bends sit only 0.5 µm apart and the independent-scatterer + N_MODES=2
interface assumption should degrade. **That did not happen.** Hybrid vs
single-domain full-FDTD, at λ = 1.5475 µm:

| | L = 0.5 µm | | L = 3.0 µm | |
|---|---|---|---|---|
| | hybrid / full-FDTD | abs diff | hybrid / full-FDTD | abs diff |
| `T_through` | 0.9621 / 0.9637 | 0.0017 | 0.9009 / 0.9022 | 0.0013 |
| `T_cross` | 0.0349 / 0.0350 | 0.0001 | 0.0961 / 0.0962 | 0.0002 |
| coupling ratio | 0.0350 / 0.0350 | **0.0000** | 0.0964 / 0.0964 | **0.0000** |

Band-wide max abs diff: through 0.0018 / cross 0.0001 at 0.5 µm, and
0.0014 / 0.0002 at 3.0 µm. **The 0.5 µm agreement is as good as the 3.0 µm
agreement**, on two independent methods sharing no cache entries.

The residual is a small `T_total` deficit in the hybrid, and it is
**length-INDEPENDENT**, which is what identifies its cause:

| coupling length | band-mean `T_total` deficit (full-FDTD − hybrid) |
|---|---|
| 0.5 µm | **+0.00170** |
| 3.0 µm | **+0.00145** |

The bends are 6× closer at 0.5 µm, yet the deficit grows only 17 %. A
bend-proximity failure would scale strongly with that separation; this does
not. It is consistent with a fixed cost at the two modal interfaces — the
N_MODES=2 basis truncation, paid once per interface regardless of what sits
between them.

**Consequence**: the left edge of the coupling-length axis is no longer
provisional. `directional_coupler_soi160_hybrid_length_sweep.py`'s
shared-FDTD construction is validated down to its 0.5 µm floor against
independent full-FDTD, which also removes the caveat carried in
`hcg_fabry_perot_length_sweep.py`'s docstring for the same reason. Update
that docstring when the queue is clear.

**No overwrites**: all six at-risk production files re-checked — UNCHANGED.

### New — 12117825, HCG-FP connecting-waveguide β @ res 80 🏃 RUNNING

Submitted 2026-08-21 alongside the four above, on `gts-sr108-gedc_main`,
`-N 1 --ntasks-per-node=24`, `-t 720`. Started immediately.

The one piece `hcg_fabry_perot_length_sweep.py` (new, see below) needs that
is neither cached nor already queued. **Confirmed a real cache MISS before
submitting** — a `cache_only=True` probe raised `CacheMissError` from
`section_mode_neff()` — so this is not redundant work.

Run through `examples/verification/python/_prewarm_hcgfp_arm_beta.py`, not
through the sweep script, *specifically* so it cannot touch 12093321's or
12093330's cache keys: the sweep script with `cache_only=0` would not WAIT
for those two jobs, it would re-simulate them. The arm is a bare 0.5 µm strip
and shares no key with either, so it is safe in parallel.

Budget basis: job 11995681's own arm calibration was 47,604 steps / 48.4 min
at res 60 on the 11-pt 1.50–1.60 µm grid. res 60→80 is (80/60)⁴ = 3.16×; the
band change adds ~1.3× (**measured** on today's two full-FDTD landings, not
the naive 2.85× — ringdown dominates these geometries). Estimate 3–4 h.

## New script — `hcg_fabry_perot_length_sweep.py`

2-D design sweep of the HCG Fabry-Perot over the coupler's **coupling length**
and the connecting waveguide's **arm length** → extinction-ratio / FSR / FWHM
heatmaps, plus a dip-count diagnostic panel.

**The whole plane is analytic.** In `hcg_fabry_perot_hybrid.py`'s closed form,
`L_ARM` enters ONLY through `phase = exp(2j·π·n_eff_arm/λ·L_ARM)` and
`coupling_length` ONLY through the EME body's diagonal phase
`exp(1j·β_body·(L_c − 2·BODY_STUB))`. The two FDTD S-bend solves and both β
calibrations are properties of the bend and the cross-section, not of how much
straight coupler sits between them. So a 40×40 plane is **4 cached
simulations + ~9 s of numpy** — the same construction
`directional_coupler_soi160_hybrid_length_sweep.py` uses, extended one axis.

**Validated end-to-end 2026-08-21 against the OLD cached band** (1.50–1.60 µm,
res 60/40/60), `HCGFP_CACHE_ONLY=1`, so all four lookups had to hit or raise:

| check | result |
|---|---|
| mean FSR at the fabricated point | **2.2455 nm** vs job 11995681's **2.245 nm** |
| resonance count | **45** vs job 11995681's **45** |
| arm n_eff | **2.138038 → 2.169312 (+1.4627 %)** vs 11995681's +1.463 % |
| 40×40 wall time | **9 s**, zero simulation |

Two defects that validation exposed, both fixed:
1. **FWHM was grid-quantised** — `_find_dips` read half-power crossings at
   sample positions, a 2–3 cell staircase at coarse density. Now interpolated
   sub-grid (`_cross`).
2. **Dip detection was O(N_dense) Python per cell** — 32 M iterations for a
   40×40 map, dominating the sweep. Vectorised; the anchor values are
   bit-identical before and after.

And one real finding from the adequacy guard added alongside: at the original
default of 20001 dense points the **narrowest cell in the plane** (0.035 nm)
got only 7 samples, which **under-reports extinction ratio** (ER samples the
dip bottom; measured 8.35 → 10.90 → 10.88 dB over 2001 → 8001 → 20001 while
FSR held at 2.2455 throughout). Default raised to **40001**, and the script
now prints samples-per-narrowest-resonance and warns below 20.

**Design finding, provisional pending the res-80 run**: the fabricated device
(L_c = 3.0, L_arm = 53.08) sits at ER ≈ 10.9 dB, while the plane peaks near
**29.6 dB around L_c ≈ 24 µm**. ER is banded almost purely in coupling length
(critical-coupling driven, near-independent of arm length); FSR carries the
1/L_arm dependence. The two axes are close to orthogonal.

**Carried-forward caveat**: short coupling lengths reuse the L_REF=3.0 S-bend
solves. Exact for the bend's own scattering, but as L_c → 2·BODY_STUB the two
bends sit ~0.5 µm apart and the independent-scatterer + N_MODES=2 interface
assumption degrades — the same issue 12093324 carries. Treat the left edge of
the coupling-length axis as provisional until the hybrid/full-FDTD comparison
at 0.5 µm is done.

**Library fix shipped with 12093330**: `HCG_FDTD_RES` scaled the resolution
but did **not** redirect `output_dir`, so any res-sweep point silently
overwrote job 11673560's res=60 production `smatrix.pkl`/`reflection.png`.
`high_contrast_grating_fdtd.py` now redirects to `res_sweep/res<R>/` whenever
that var is set. Pre-existing outputs were mtime-baselined before submission
and must be re-checked when these land.

---


**★ TWENTY-EIGHTH SNAPSHOT — 2026-08-18, verified against
`squeue -u jhiesener3`. This is the current one. Everything below the next
rule is the twenty-seventh snapshot (2026-08-10), not reverified this pass.**

Two new jobs only — a resolution replica of the 2x2 MMI dense stub sweep.
Both are SUBMITTED AND RUNNING, neither has landed, so nothing below is a
result yet.

| job | state | what |
|---|---|---|
| **11995441** | ✅ **COMPLETED** 19:06:00 (2026-08-18 17:50 → 2026-08-19 12:58), exit 0:0 — **26.5% of its 72h budget** | `beef_mmi2x2_bodystub_sweep_alt_2nm_res30` — `mmi_2x2_body_stub_sweep_alt.py` at **res=30**, a replica of the manuscript's dense 2 nm grid (51 monitors, 0–0.10 µm past the taper end, `x1_fixed=-13.85`, 11 wl pts, 14 interface modes, top-arm excitation) with ONLY `self.res` changed. **Verified against the `.npz` itself, not the job-log print**: `M_in1` (11, 51, 14, 1) complex128, all finite, no zeros, `RESOLUTION=30.0` stored, and `body_stubs`/`wavelengths`/`x1_fixed` all bit-identical to the res=40 file — so the two datasets are directly differenceable. Result in the section below. → `data/mode_coeffs_body_stub_alt_dense_res30_2nm_0-100nm.npz` |
| **11995442** | ✅ **COMPLETED** 2-09:30:58 (57.5 h), exit 0:0 — **34% of its 7-day budget**; npz written 2026-08-21 03:23 | `beef_mmi2x2_bodystub_sweep_alt_2nm_res50` — same replica at **res=50**. `Elapsed run time = 207,049 s`. Steady-state rate held at ~3,550 s/position for the whole run; the pre-submission estimate was ~56 h at res³ scaling, actual 57.5 h — **the res³ decomposition scaling was correct**. → `data/mode_coeffs_body_stub_alt_dense_res50_2nm_0-100nm.npz` |
| **11982320** | 🏃 RUNNING (6h58m elapsed at snapshot time) | `beef_inv_res30_input` — another session's work, not independently verified this snapshot |

### Landed this pass — job 11995681, `beef_hcgfp_selfcontained` ✅ COMPLETED

| job | state | what |
|---|---|---|
| **11995681** | ✅ COMPLETED 18:25→19:15, **49:34**, exit 0 (1 node/24 cores, `gts-sr108-gedc_main`) | `hcg_fabry_perot_hybrid` rewritten to solve all three of its pieces itself against `example_sim_cache/` instead of loading two other scripts' pickles, **and to calibrate the connecting arms on the FDTD's own β** — the last uncalibrated segment in the device (2 × 53.08 µm against the coupler's 23 µm). Only the arm calibration time-stepped (47,604 steps, 48.4 min, res 60); mirror, both S-bends, coupler body calibration and EME body all hit cache. Task 2 (the `HCGFP_CALIBRATE_ARM=0` baseline) was 21.6 s, all hits. → `example_output/high_contrast_grating_fabry_perot/hcg_fabry_perot_hybrid{,_uncalibrated_arm}/` |

**Result**: arm n_eff 2.138038 (MPB) → **2.169312** (FDTD), **+1.463 %**, worth
**+4.28 FSR** of comb position over the 4 × L_ARM = 212.32 µm round trip.
45 resonances, mean FSR **2.245 nm** (baseline 2.268 nm), Q 9,998–16,071.
Rebuilt mirror `r(λ)` matched its published pickle to **0.000e+00** and the
coupler's `t`/`k` to 1.8e-15, so the rebuild is a refactor, not a new device.
Full write-up: `resolved/high_contrast_grating_fabry_perot_validation.md`
step 7.

**Account note**: submitted under `gts-sr108-gedc_main`, not the
`gts-sr108` in most batch files — `sshare -U -u jhiesener3` at submit time
showed EffectvUsage 0.979 (`-ece`) / 0.716 (`gts-sr108`) / **0.064**
(`-gedc_main`). Started immediately. Re-check before reusing; which of the
three has headroom has now moved three times.

**Library fix that this job depended on**: `beef/cache/config.py` gained
`_DEFAULT_OMITTED_ATTRS`, because `FDTDSolver.minimum_run_time` had stranded
all 48 `solve_reflection` cache entries (the HCG mirror's ~55 min res-60 3-D
production run among them) while simultaneously being invisible to
`solve`/`solve_section`/`solve_section_multi_length`, where the modal
multiplexer actually sets it non-zero. See
`resolved/section_mode_beta_multi_monitor_validation.md` step 11.

### Landed — job 12260555, `beef_migrated_cascade_examples` — ✅ COMPLETED 2026-08-24, **2:23**, exit 0:0

End-to-end re-run of the two examples migrated off the hand-rolled einsum
(`mmi_2x2_hybrid.py`, `directional_coupler_hybrid.py`) via
`batch/run_migrated_cascade_examples.batch`.

**Why it existed, given the numbers were already verified:** the `cache_only`
harness that proved the migration bit-identical *replicated* each script's
composition rather than importing it, so it never exercised the scripts
themselves. Specifically untested by it: `plot_beta_calibration()` fed from the
new `hybrid.last_calibrations[1]` (previously the examples held the
`BetaCalibration` object directly), `plot_hybrid_geometry()` against a solver
that now carries `mirror_pairs`, the SMatrix save/`check_passive()` path, and
the DC's `section_mode_neff()` log line.

**Result: all four grids bit-identical to the pre-migration backups**
(`T_bar`/`T_cross` for the 2×2, `T_through`/`T_cross` for the DC), every plot
regenerated (calibration `fits/` + `geometry/` file mtimes confirmed current —
the *directory* mtimes stay old because files are overwritten, not created), no
warnings or tracebacks. **2:23** total: both configurations were cache-resident,
as expected.

The migration is now verified at the script level, not just the composition
level. The `*.pre_migration_bak` copies of both reference `.npz` files are kept
beside them.

### Archive snapshot — 2026-08-24

`archive/example_output_20260824/` and `archive/example_sim_cache_20260824/`
(`cp -a`, matching the existing `*_YYYYMMDD` convention — the archives carry no
manifest file, so none was added). Taken **after** job 12260555, so this is the
first snapshot holding the post-migration state of both cascade examples.

Verified rather than assumed, because seven other sessions' jobs were running at
the time: regular-file bytes identical (243,198,070 for the cache; 382,762,112
for the outputs), 2,093 cache files and 789 output files on both sides, `find`
trees `diff`-clean, and **nothing written to either source tree during the 91 s
copy** — so no torn writes despite the concurrent jobs. (`du -sh` reports
240M vs 238M for the cache; that is directory-inode block allocation on a fresh
copy, not missing data — regular-file bytes match exactly.)

### Landed — job 12159086, `beef_mmi2x2_cascade_golden` (**`CAPTURE_REAL_M_IN0=1`**) — ✅ COMPLETED 2026-08-21 17:14:55, **45:01**, exit 0:0

**Well under budget** (45 min vs. the ~1.5-2.5 h estimated) — the second
input-side FDTD arm was the only thing that had to time-step; everything else
was cache-resident.

**Result: the hard-coded `MODE_SIGNS` is correct on all 14 modes.** Implied by
the measured pair == hard-coded ==
`[1,-1,1,-1,1,-1,1,-1,-1,-1,-1,1,-1,1]`, `agree: True`, and a per-mode
re-derivation from the `.npz` confirms every sign, not just the dominant ones.

The log's `max rel = 9.576e-01` is **not** a failure: the parity relation
constrains *signs*, and its amplitude half degrades only on noise-floor modes.
Modes 0-2 agree to 1.7e-05; modes 3-6/10/13 to ≤1.8e-03; modes 7/9/11/12
(|M| ~1e-05) to 0.1-0.25; mode 8 (|M| = **1.3e-06**, i.e. 5e5 below the
dominant mode) to 2.8. **Cost of the whole assumption, end to end: 1.585e-05
on `T_cross`**, with `T_bar`/`coeffs_bar`/`M_in1`/`M_body` bit-identical
(the bar path never touches `M_in0`). That saves a 45-minute 3-D run.

**Two goldens now, on purpose.** The script rewrites `cascade_golden.npz` in
place and its cross columns become a *different quantity* (measured rather than
derived `M_in0`). Kept apart rather than letting one filename mean two
configurations: **`cascade_golden_real_m_in0.npz`** is this run;
`cascade_golden.npz` was restored to the `SYMMETRIC_MODE_REUSE=True` capture
that steps 2-4 regress against, and the migration was **re-verified
bit-identical against the restored file** on all four arrays. The pre-submit
backup is what made that separation possible — take it.

**Step 3(b) is unblocked but de-prioritised**, and the reason is a finding: its
motivating hazard (a stale vector silently corrupting the cross path) is now
covered by `MirrorSpec`'s per-call geometric precondition check, and the vector
itself is measured correct. Full write-up in
`resolved/multi_excitation_cascade_validation.md` step 3(b).

<details>
<summary>Submission entry (2026-08-21) — kept for the sizing record</summary>

The one job `resolved/multi_excitation_cascade_validation.md` leaves outstanding after
steps 2-4 closed. It adds the **second input-side FDTD arm** that
`SYMMETRIC_MODE_REUSE` exists to avoid (~1.5-2.5 h of real 3-D solve; the rest
of the capture is cache-resident and cost 38 s in job 12070199). That measured
`M_in0` is the regression target **step 3(b)** needs before `MODE_SIGNS` can be
*derived* from mode parity instead of declared — today `MirrorSpec` takes it as
a hard-coded argument, geometrically precondition-checked but not derived.

Submitted as `sbatch --export=ALL,CAPTURE_REAL_M_IN0=1
batch/verification/run_capture_mmi_2x2_cascade_golden.batch` (the batch file
already reads the variable through, defaulting to 0). Verified with `scontrol
show job` that `CAPTURE_REAL_M_IN0=1` actually reached the job environment —
the whole point of the run is that one flag, and a silently-defaulted 0 would
produce a second 38 s cache-hit that looks like success. `gts-sr108-gedc_main`
(`sshare` at submit time: **0.071** `-gedc_main` / 0.733 `gts-sr108` / 0.977
`-ece` — re-checked per [[slurm-account-gts-sr108-ece]], not assumed).
`PENDING (Resources)`, 2 nodes × 24 ranks, `-t 420`.

**The script rewrites `cascade_golden.npz` in place**, which is the same file
this session's steps 2-4 verified against. Backed up first to
`cascade_golden_preM_IN0_20260821.npz.bak`. The rewrite should be a strict
superset (everything else recomputes from cache, measured bit-identical
already) — if it is not, that is itself a finding, and the backup is what makes
it visible instead of silently overwritten.

**What to do when it lands:** read the printed `signs implied by the measured
pair` against the hard-coded `MODE_SIGNS` (the script already prints both plus
`max |measured - derived|`), record it in the validation doc's step 3(b), and
only then implement field-parity derivation with the measured `M_in0` as its
regression target.

</details>

### Landed — job 12069544, `beef_mmi_2x2_hybrid_gds` — ✅ COMPLETED 2026-08-21 00:06:32, **4:31:02**, exit 0:0

**Result** (read from the two `.npz` files, not the job log; `mm_lengths` and
`wavelengths` bit-identical between runs, so the grids are directly
differenceable):

| 201 × 11 grid | max abs Δ | max rel Δ | RMS Δ | signed mean Δ |
|---|---|---|---|---|
| `T_bar` | 6.556e-04 | **1.603e-03** | 2.684e-04 | +1.872e-04 |
| `T_cross` | 5.218e-04 | **1.346e-03** | 1.973e-04 | +5.252e-05 |

Looser than step 5's "solver noise" bar, smooth in λ and positively biased —
a systematic offset, not stochastic noise. **Localized rather than accepted**,
in three cache-only measurements:

1. **`beta`/`body_neff` bit-identical (0.000e+00)** — and computed *fresh* on
   the GDS side, so this is the independent re-verification of the meep/FDTD
   layer that this job's cache-served sibling 12070199 could not give. The
   seven uncommitted library edits move nothing there either.
2. **`M_in1` carries 100 % of it** — swapping `neff` alone moves the grids by
   exactly 0.000e+00.
3. **Cause: prism ORDER.** Both backends emit the *same five polygons* for the
   input section (same areas, canonicalized vertex-for-vertex identical) in a
   different order, and meep's subpixel averaging is order-sensitive where
   prisms share an edge. Rebuilding the ε grid with GDS's polygons **reordered
   to `MMI`'s** gives a **bit-identical** structure (0 of 2.52 M voxels differ,
   vs. 872 differing natively, all at the two seams).

**Verdict: acceptance PASS, no code change** (user's call 2026-08-21). The GDS
path is proven correct by the reordered-ε test; the residual is a meep property,
not a GDS defect. Full write-up in
`resolved/gds_geometry_hybrid_cascade_validation.md` step 5 — **that plan is now CLOSED**.

**Spun out of it — a library-wide question, not GDS business.** Merging the
input section's five abutting polygons into the one 14-vertex polygon they
geometrically are moves **41,677** ε voxels — **48× the ordering effect**. meep
infers an interface normal at every seam between geometry objects even when both
carry the identical medium, so *every* BEEF device drawn as abutting prisms
carries this artifact. Deliberately **not** fixed here (would move all six
devices' numbers, invalidate the 566-entry FDTD cache, and change the
manuscript's published 2×2 results); wants its own plan.

<details>
<summary>Original submission entry (2026-08-20) — kept for the sizing record</summary>

`resolved/gds_geometry_hybrid_cascade_validation.md` step 5's full production-scale
run: `mmi_2x2_hybrid_gds.py` at `RESOLUTION=40`, `N_WL=11`, the 201-point
`mm_length` sweep — identical settings to `mmi_2x2_hybrid.py`'s (job
11585844_2, which took 4:11:55 on this same 2-node/24-rank shape), with
`GdsGeometry` reading `mmi_2x2_twin.gds` in place of the parametric `MMI`.
This is the run steps 1-4 and step 5's cheap-first half (both already
MEASURED) were building toward; the doc's transfer-matrix agreement number
stays ARGUED until this lands.

New batch script `batch/run_mmi_2x2_hybrid_gds.batch` (+ Makefile target
`mmi_2x2_hybrid_gds`), same 2-node/24-rank/`-t 420` shape as the parametric
run's array task, `-o reports/%j.out`. No `--dependency` needed — the
comparison in the script's own step 7 reads
`example_output/mmi_2x2/mmi_2x2_hybrid/mm_length_sweep.npz`, which already
exists on disk (2026-08-18, same 201×11 grid, confirmed via `np.load` before
submitting) rather than a sibling array task's fresh output.

Submitted under `gts-sr108-gedc_main` (`sshare -U -u jhiesener3` at submit
time: EffectvUsage 0.979 `-ece` / 0.720 `gts-sr108` / **0.067**
`-gedc_main` — re-checked per [[slurm-account-gts-sr108-ece]], not assumed
from a prior snapshot). `squeue` initially showed `PENDING
(ReqNodeNotAvail, May be reserved for other job)`; the only active
reservation (`md-26-10`) doesn't start until 2026-10-05, so this is
ordinary queueing for 2 free nodes, not a real block. Not yet landed —
update this entry and the validation doc's step 5 with the real
`T_bar`/`T_cross` max-|Δ| numbers once it completes.

The earlier login-node smoke run (RES=12/N_MODES=4/N_WL=2/3-pt sweep,
mentioned in the validation doc as "in flight") is no longer running on this
node and left only a partial `example_output/mmi_2x2/mmi_2x2_hybrid_gds/`
(geometry plot + section plots, no `mm_length_sweep.npz`/SMatrix) — it ended
without the doc being updated with an outcome. Its cache entries (RES=12)
don't collide with this RES=40 production submission's fingerprints, so no
cleanup was needed before submitting.

</details>

### Landed this pass — job 12070199, `beef_mmi2x2_cascade_golden` — ✅ COMPLETED 20:00:17, **38 s**, exit 0:0

`resolved/multi_excitation_cascade_validation.md` step 1's device-capture half:
`examples/verification/python/capture_mmi_2x2_cascade_golden.py` via
`batch/verification/run_capture_mmi_2x2_cascade_golden.batch`. Captures the
2×2's hand-rolled `M_in1`/`M_in0`/`M_body` plus the composed bar/cross grids
as the bit-level regression target for steps 2–4, and measures the
einsum-vs-`Circuit` agreement on the **real** device matrices rather than only
the synthetic ones the new unit tests use.

Writes to `example_output/verification/mmi_2x2_cascade_golden/` **only** — it
deliberately does not touch `example_output/mmi_2x2/mmi_2x2_hybrid/`, which
job 12069544 (running) reads at its step 7.

**Cost: 38 s** — the supersets hit, no meep time-stepping in the log at all.
Worth recording why the estimate was open-ended: a plain `lookup()` probe of
the three shipped calls missed on all three, but `solve_section()` also falls
back to a full-matrix `lookup()` and `lookup_superset()`
(`fdtd.py:2924-2943`), which a hand-built probe cannot safely reproduce, so
the job was sized for the miss case (2 nodes × 24 ranks, `-t 420`, cf.
11585844_2's 4:11:55 on the same shape). **Do not size a job off a hand-built
`lookup()` probe.** `gts-sr108-gedc_main` (EffectvUsage 0.067 at submit time).

**Result**, read from `cascade_golden.npz` rather than the log: einsum vs
`Circuit` agrees to **4.0 eps** on both the bar and cross paths over the full
201 × 11 grid, reproducing on real device matrices what
`tests/test_multi_excitation_cascade.py` establishes synthetically. Full
write-up in `resolved/multi_excitation_cascade_validation.md` step 1.

**Follow-on, 2026-08-21 — the plan this golden served is now CLOSED, with no
further job needed.** `resolved/multi_excitation_cascade_validation.md` steps 2-4 landed
against this capture: the 2×2's hand-rolled einsum cascade is retired in favour
of two `HybridSolver.solve_length_sweep()` calls plus the new `beef.MirrorSpec`
cross-excitation mirror declaration. Verified `cache_only=True` against this
`.npz` — **bit-identical** (0.000e+00) on `T_bar`, `T_cross` *and* the complex
coefficients, at **1** `solve_section()` call across both excitations (the
hand-rolled path's own count; step 2 without `MirrorSpec` would have cost 2).
`Circuit` turned out not to be needed at all. 757 tests green.

**The one job this leaves outstanding — now SUBMITTED, see below.**

**Bonus result — it closed an open issue in the GDS validation.** It also
recomputed the parametric 2×2's `T_bar`/`T_cross` on the current library, which
`resolved/gds_geometry_hybrid_cascade_validation.md` step 5 needed to separate a
geometry-backend difference from the seven uncommitted library edits behind its
2026-08-18 reference. Both 201 × 11 grids came back **bit-identical**
(max |Δ| 0.000e+00) — the library edits moved nothing on this device, and the
08-18 reference is a valid comparison target after all. Caveat recorded in that
doc: the capture was cache-served, so it re-verifies everything downstream of
the cached section matrices but not the meep/FDTD layer, which did not re-run.

### ⛔ FINAL VERDICT — the transition region is NOT mesh-converged at ANY of res 30/40/50

All three resolutions have landed. Adding res=50 settles the question the
res=30 run could only pose, and the answer is the unfavourable one.

| region | max \|Δ\|c\|\| 30→40 | max \|Δ\|c\|\| 40→50 | ratio | verdict |
|---|---|---|---|---|
| `body_stub` ≥ 0.030 µm (plateau) | 0.0283 | **0.0014** | **0.16** | ✅ **converged** — refinement shrinks the difference ~6× |
| `body_stub` < 0.028 µm (transition) | 0.0817 | **0.0877** | **1.07** | ⛔ **NOT converged** — refinement does not shrink the difference at all |

**Going from res 40 to res 50 changes the transition region by MORE than going
from 30 to 40 did.** A converging sequence would give a ratio well under 1;
this gives 1.07. The apparent order of convergence there is ~0.

The mechanism is a monotone transfer of power from mode 0 into mode 2 that has
not saturated by res 50. At `body_stub` = 0.010 µm (the 40→50 worst case):

| mode | res 30 | res 40 | res 50 |
|---|---|---|---|
| 0 | 0.700 | 0.688 | 0.663 |
| 1 | 0.705 | 0.704 | 0.705 |
| **2** | **0.074** | **0.143** | **0.231** |

Mode 1 (the symmetric mode) is completely insensitive to resolution. Mode 2
roughly triples across the range and its increments are still *growing*
(+0.069 then +0.088), so it is not approaching a limit from either side.

**This is NOT the known mode-identity/correspondence artifact.** Tested
explicitly: if the 14-mode basis were permuting between resolutions, sorting
each position's \|c\| spectrum would collapse the difference. It does not —
at `body_stub` = 0.010 µm the per-index max difference is 0.0877 and the
**sorted** max difference is also 0.0877. The power really is moving between
physically distinct modes, not being relabelled. (Two minor exceptions worth
knowing: at `body_stub` = 0 and 0.004 µm the sorted difference IS smaller than
the per-index one — 0.0416 vs 0.0501 and 0.0122 vs 0.0290 — so there is some
partial reordering right at the taper face; and modes 10/11 do swap a tiny
amount of amplitude at res 50, but both sit at \|c\| ≲ 0.016, i.e. \|c\|² ≲
3e-4, far below the 0.05 dominance threshold and irrelevant to the verdict.)

Power conservation is healthy and itself converged at every resolution
(Σ\|c\|² = 0.9988 / 0.9992 / 0.9992), which is what makes the redistribution
interpretation sound rather than a leak.

**What the manuscript may and may not say.**
- ✅ Plateau values at `body_stub` ≳ 0.03 µm are mesh-converged to ~0.0014 in
  \|c\|, demonstrated with a proper three-point sequence rather than two runs
  agreeing. Quote them.
- ⛔ **Nothing about the shape, slope, or onset of the transition below
  ~0.028 µm is supportable from res 30/40/50 data.** In particular mode 2's
  amplitude anywhere in 0–0.024 µm is not a converged number at res 40, and
  the res=40 dense sweep alone would have given no hint of that.
- ⛔ The `body_stub` = 0 point (monitor sitting exactly on the taper
  discontinuity) is the worst-behaved of all: non-monotone in resolution
  (0.0409 → 0.0856 → 0.0356 for mode 2) and partially reordered. Treat it as
  unusable, not merely uncertain.

**Recommended follow-up if the transition shape actually matters.** Do NOT
re-run the full 51-position grid at higher resolution — the measured res³
decomposition scaling makes that punishing (a res=60 51-position job projects
to ~100 h, res=80 to ~250 h). Instead run a *narrow* grid of ~11 monitors over
0–0.020 µm only, where the non-convergence lives. Time-stepping is a
whole-domain cost paid once regardless of monitor count, so at res=60 that is
~12.5 h time-stepping + 11 × ~6,100 s decomposition + ~1.5 h setup ≈ **33 h**
(fits a 72 h budget); res=80 lands near ~85 h and needs the 7-day budget. Two
such jobs would establish whether the mode-0 → mode-2 transfer saturates at
all, which is the actual open physics question.

### RESULT (superseded by the three-resolution verdict above) — res=30 vs res=40

First actual answer this convergence study produced. Comparing
`|mode coefficient|` between the two 51-position grids (identical monitor
positions, so this is a straight difference, no interpolation), at the
sweep's centre λ = 1.5475 µm and also across all 11 wavelengths:

| region | max \|Δ\|c\|\| (dominant modes) | verdict |
|---|---|---|
| `body_stub` ≥ 0.040 µm | **≤ 0.005** (flat floor at 0.0024, in weak mode 5) | mesh-converged |
| `body_stub` ≥ 0.020 µm | 0.078 | not converged |
| `body_stub` 0.012–0.024 µm | **0.074–0.082** (peak 0.0817 at 0.018 µm) | **not converged** |
| `body_stub` = 0 | 0.045 | not converged |

The disagreement is almost entirely **mode 2**, and it is large in relative
terms, not just absolute: at `body_stub` = 0, mode 2's `|c|` is **0.0409 at
res 30 vs 0.0856 at res 40 — a factor of ~2**. Modes 0 and 1 agree to
≤ 0.005 and ≤ 0.003 respectively everywhere. Total captured power is
essentially identical and healthy on both (Σ\|c\|² = 0.9988 res 30 vs
0.9992 res 40, spread < 0.0005 across all 51 positions), so this is NOT a
basis-completeness or power-conservation problem — it is specifically the
*distribution* of power among the low-order modes in the near-taper region.

Physically consistent: `body_stub` = 0 puts the monitor exactly at the taper
end, on the geometric discontinuity, where the field has the most rapidly
varying / evanescent content and is therefore the most mesh-sensitive. As
the monitor retreats into uniform `mm_width` material that content decays
and the two resolutions converge. Notably the crossover (~0.04 µm) coincides
with the point `mmi_2x2_body_stub_sweep.py`'s own 10-point sweep already
identified as where modes 0/2 stop changing physically (0.03–0.05 µm) — so
the physical settling and the mesh settling happen at the same place.

**What this licenses and forbids.** Any quantity read off the dense sweep at
`body_stub` ≳ 0.04 µm is mesh-converged at the 0.005 level and safe to
quote. Any claim about the *shape of the transition* below ~0.024 µm is NOT
resolution-independent between res 30 and res 40 and must not be quoted from
a single resolution.

**res=50 (11995442) is the tiebreaker** and has not landed. If res 50 agrees
with res 40 in the sub-0.024 µm window, res 40 is converged there and res 30
is simply too coarse. If res 50 keeps moving, none of the three is converged
in that window and the transition region needs a finer mesh than 50 before
anything about its shape can be claimed.

### Corrected cost model — there is a ~1h fixed setup offset, so don't use (elapsed − time-stepping)/N

The res=30 run's own numbers close cleanly and give the right way to size
these jobs:

| term | res=30 measured |
|---|---|
| time-stepping | ~3,500 s (0.97 h) — 26,589 steps; mean 0.1355 s/step over 828 progress lines, 0.1258 s/step steady-state (the 7.76 s/step max is the first line, i.e. startup) |
| one-time setup | **~3,500 s (0.97 h)** — structure init, input-side source mode solve, DFT monitor allocation |
| decomposition | 51 × **1,210 s/position** = 61,710 s (17.14 h), slope fitted from a 10-min-cadence sample series (R² = 0.9984) |
| **total** | **68,750 s (19.10 h)**, vs `Elapsed 19:06:00` from `sacct` |

The naive `(elapsed − time-stepping)/N` estimate gives 1,278 s/position —
**68 s/position too high**, and that gap is exactly the ~3,500 s fixed setup
cost smeared over 51 positions. This matters when N is small: an in-flight
job only a few positions in will look far more expensive per position than
it really is. **Fit the slope; don't divide the total.**

**Decomposition cost ACCELERATES with resolution — it is not a single power
law.** Measured per-position slopes, now all from fitted series or
well-diluted totals:

| resolution | s/position | exponent vs. previous |
|---|---|---|
| 30 | **1,210** (fitted, R² = 0.9984, 51 positions) | — |
| 40 | 1,790 (total ÷ 51; setup-inflated by ≲ 5%, so true slope ~1,720–1,790) | res^1.36 |
| 50 | **3,550** (fitted, R² = 0.9995, 73 samples over 12 h) | **res^3.07** |

End-to-end 30 → 50 is res^2.11, but that average hides the shape: the step
from 40 to 50 costs roughly res³ while 30 to 40 costs roughly res^1.4. **Any
res ≥ 60 run sized from a res² assumption will be badly under-budgeted.**
Likely mechanism is the cross-section working set outgrowing a cache level
between res 40 and 50 (res 50's cell is 282 × 300 × 125 = 10.6M pixels), but
that was inferred from the timing curve, not measured with counters.

Hardware is NOT the explanation: all four nodes involved are Intel Gold 6226
Cascade Lake / avx512, and the only difference is that res=40's node ran at
2900 MHz vs 2700 MHz for res=30/res=50 — a ~7% effect, and in the direction
that makes res=40 look *cheaper*, i.e. it slightly *understates* the 40 → 50
jump rather than creating it.

**Corrections to earlier versions of this snapshot** (both were wrong, and
both were wrong in the direction of under-stating res-50 cost):

- An earlier revision called the `res^3.01` exponent "overstated, computed
  the naive way" and replaced it with ~res^1.36 overall. That was wrong. The
  setup offset is real (see the table above) but it is only ~5% of res=40's
  per-position figure and does not touch res=50's *fitted* slope at all.
  res=50's fitted slope (3,550 s/pos) is within 2% of the naive average
  computed when it was only 9 positions in (3,465 s/pos), so the naive
  estimate happened to be fine there. **~res^3 for 40 → 50 stands.**
- An earlier revision claimed res=50's timing was "contaminated as a scaling
  datum" by node contention from a 4,280-CPU co-tenant (job 12013145), citing
  a deceleration from +31 to +15 lines per 10-minute tick and a settled rate
  of 5,775 s/position. **That does not hold.** Those were sub-hour windows
  covering well under one position each — sampling noise aliased against the
  ~1 h position cycle. Over 73 samples spanning 12 h the rate is
  3,550 s/position with R² = 0.9995, and every sub-window (first half, second
  half, last 6 h, last 2 h) agrees to within 2%. That span brackets the
  4,280-CPU job's departure from the node with no detectable change in
  throughput, which is direct evidence against the contention hypothesis.
  **res=50's timing is usable as a scaling datum; no `--exclusive` re-run is
  needed.**

Methodological lesson worth keeping: a per-position rate for these jobs needs
a sample window covering **several positions** (≳ 3–4 h at res 50), not
minutes. Short windows produced three mutually inconsistent rate estimates
here (2,552 / 4,615 / 5,775 s/pos) before the long series settled at 3,550.

### Where the res=30/res=50 time budgets came from (measured, not guessed)

The res=40 2 nm production run (**job 11880200**,
`batch/verification/reports/11880200.out`) is the only real cost datum for
this grid, and it splits cleanly:

| phase | cost | how it was counted |
|---|---|---|
| time-stepping | 8,874 s (2.5 h) | 35,425 timesteps × 0.2505 s/step, from the log's own `on time step` progress lines and `run 0 finished at t = 442.8125 (35425 timesteps)` |
| eigenmode decomposition | 91,293 s (25.4 h) | total − time-stepping, i.e. **1,790 s per monitor position** over 51 positions (consistent with the ~1,712 s/position measured mid-run on the earlier 126-position attempt, job 11847954) |
| **total** | **100,167 s (27.8 h)** | the log's own `Elapsed run time` |

Scaled: time-stepping ~res⁴ in 3-D (three spatial dims + the
Courant-limited step count); decomposition is dominated by the per-monitor
MPB cross-section solves (11 wl × 14 modes each) so it scales somewhere
between res² and res³. **The decomposition exponent is EXTRAPOLATED, not
measured** — when these two jobs land, size any future run from their real
rates instead. Margins are deliberately fat because `solve_section()` saves
its result ONCE at the very end (no incremental checkpoint), so a timeout
loses the whole run including the already-completed time-stepping.

### Two things changed about the job shape, both evidence-backed

1. **1 node/4 cores, down from 2 nodes/48 cores.** Job 11880200's own log
   records `Using MPI version 4.1, 1 processes` and `mem=3361476K` (~3.4 GB)
   against 384 GB requested — there is no `mpirun`/`srun` in these batch
   scripts, so 47 of the 48 requested cores sat idle for 27.8 h. res=50's
   cell is 10.6M pixels (1.95x res=40's 5.4M) and its DFT monitor bank grows
   1.56x, so ~6–8 GB peak; 4 × 8G = 32 GB keeps ~4x headroom on the largest
   leg.
2. **`OMP_NUM_THREADS=1` now explicit.** Previously left to meep's implicit
   default. Threads measurably HURT this workload (`OMP_NUM_THREADS=8` ran
   ~3.7x SLOWER than `=1` on a 101-position smoke case), so pinning removes
   a variable the res=40 baseline left unpinned. Physics is unaffected.

### Script change enabling this

`mmi_2x2_body_stub_sweep_alt.py` gained a `BODY_STUB_RES` env var (default
40.0 = `LEGACY_RES`). Output filenames carry a `res{N}_` tag for every
resolution EXCEPT 40, whose two already-published files
(`mode_coeffs_body_stub_alt_dense_{2nm_0-100nm,5nm_0-250nm}.npz`, read by
`mmi_2x2_body_stub_dense_analysis.py` by exact name) keep their untagged
names — so the new runs are self-describing with no manuscript data renamed.
`RESOLUTION` is also stored inside every `.npz`. `example_sim_cache/`
separates the three resolutions automatically: `res` is in
`solve_section()`'s cache fingerprint allowlist (`beef/cache/config.py`,
`_METHOD_ATTRS`), so neither job can false-positive onto the res=40 cached
result. Verified before submitting: default (unset) still resolves to
res=40 and the untagged legacy filename, and a sub-minute real-meep smoke
(res=10 / 2 modes / 1 wl / 6 monitors) ran `solve_section()` with
`measure_port`/`measure_margins` end to end and wrote a correctly
`res10_`-tagged `.npz`.

**Analysis script — WRITTEN 2026-08-20.**
`beef_manuscript_figures/2x2_MMI_figures/mmi_2x2_body_stub_res_convergence.py`
(new; `mmi_2x2_body_stub_dense_analysis.py` is untouched and still does the
2 nm-vs-5 nm *sampling* comparison at res=40). It reads the 2 nm grid at
every resolution that has landed, hard-errors if `(body_stubs, wavelengths,
x1_fixed, N_MODES)` don't match across files rather than silently
resampling, and runs fine on a subset — res=50 slots in automatically when
11995442 finishes, no edit needed. Emits three figure pairs (png+svg) into
`plots_stub_sweep/`: mode amplitude per resolution, **max \|Δ\|c\|\| vs
`body_stub` per adjacent resolution pair** (the actual convergence plot), and
total captured power per resolution. Re-derives the res30-vs-res40 numbers
above through an independent code path — 0.0817 worst at `body_stub` = 0.018
µm in mode 2, converged below 0.005 from 0.040 µm — which is a real
cross-check, not a restatement.

---

**★ TWENTY-SEVENTH SNAPSHOT — 2026-08-10, verified against
`squeue -u jhiesener3`. This is the current one. Everything below the rule
is the twenty-sixth snapshot (2026-08-05), not reverified this pass except
where noted.**

| job | state | what |
|---|---|---|
| **11791990_1** | ✅ **COMPLETED** 01:52:03 (16:06:45→17:58:48), exit 0:0 | `beef_mmi_2x2_examples`, array index 1 only (`mmi_2x2_fdtd`) — rerun at `N_WL` 11→51 across the same C-band 1.530–1.565 µm span (user request: more wavelength resolution). `check_passive()` held (script reached its post-solve print block and saved `smatrix.pkl` without raising). At λ=1.547 µm (center of 51-pt sweep): bar 0.4560, cross 0.5207, total 0.9767, IL −0.10 dB — close to ideal 50/50 (T_each=0.5) at the nominal self-imaging length. Verified against `example_output/mmi_2x2/mmi_2x2_fdtd/smatrix.pkl` (written 17:58) and `transmission.png`, not just the job-log `:.4f` lines |
| **11791991_0** | ✅ **COMPLETED** 01:44:41 (16:06:45→17:51:26), exit 0:0 | `beef_mmi_2x2_fdtd_alt_length` array task 0, L24 (mm_length=24, total 32.00 µm) — same `N_WL`→51 change. `check_passive()` held. At λ=1.547 µm: bar 0.2786, cross 0.3930, total 0.6716, IL −1.73 dB — expected off-nominal degradation away from the 28 µm self-imaging length. Verified against `L24.00/smatrix.pkl` (written 17:51) |
| **11791991_1** | ✅ **COMPLETED** 02:29:09 (16:06:45→18:35:54), exit 0:0 | `beef_mmi_2x2_fdtd_alt_length` array task 1, L32 (mm_length=32, total 40.00 µm — longest body of the three, hence longest runtime). `check_passive()` held. At λ=1.547 µm: bar 0.3864, cross 0.3468, total 0.7331, IL −1.35 dB. Verified against `L32.00/smatrix.pkl` (written 18:35) and `transmission.png`, not just the job-log `:.4f` lines |

**Rerun complete — all three `N_WL` 11→51 mmi_2x2 FDTD jobs finished clean, no errors, `check_passive()` held on every one.** Rollup at λ=1.547 µm (center of the new 51-point C-band sweep, 1.530–1.565 µm):

| mm_length | bar | cross | total | IL (dB) | wall time |
|---|---|---|---|---|---|
| 28.00 (nominal) | 0.4560 | 0.5207 | 0.9767 | −0.10 | 1:52:03 |
| 24.00 | 0.2786 | 0.3930 | 0.6716 | −1.73 | 1:44:41 |
| 32.00 | 0.3864 | 0.3468 | 0.7331 | −1.35 | 2:29:09 |

Matches the expected physical picture: the nominal 28 µm self-imaging length is the only one near ideal 50/50 (IL ~ 0 dB); both off-nominal lengths show real insertion-loss penalties. All figures/`smatrix.pkl` regenerated at the finer 51-point resolution the user asked for; nothing else about the device configuration changed from the 11-point baseline.

---

**★ TWENTY-SIXTH SNAPSHOT — 2026-08-05, verified against
`squeue -u jhiesener3`.**

| job | state | what |
|---|---|---|
| ~~11688645~~ → ~~11693107~~ → **11707932** | ❌→❌(killed)→🏃 | `beef_modal_mux_fdtd` — the real (non-smoke) res=20 full-device FDTD run (gap 3, README.md check 4). **11688645 TIMED OUT at 4h; 11693107 (24h budget) killed by the user's own diagnosis at 18h53m**, both stuck in just the first input mode (TE0), never converging. User found a real comparison point — an archived TO campaign for this exact device (`04_weary_paw/`) completed a full forward+adjoint cycle at HIGHER resolution in ~70-80 min/iteration — and directed a core-count increase. **Rescaled 1 node/24 ranks → 4 nodes/96 ranks, resubmitted as 11707932**, confirmed launching on all 96 ranks. **Update 2026-08-06, ~9h06m elapsed: TE0 fully decayed (t=77,197.45, 3,087,898 timesteps, ~7.95h wall) — the first time this campaign has ever finished even the first input mode, far past both single-node attempts' kill points. Now into TE1, ~1h15m in, healthy. Still not landed, but this is real evidence core count (not `decay_tol`) was the dominant lever.** See the detail section below for the full account |
| **11688701** `_[..]` | 🏃 RUNNING | `beef_hcg` — unrelated to this thread, another session's work, not independently verified this snapshot |
| **11691379** | 🏃 RUNNING | `beef_dc_...` — unrelated to this thread, another session's work, not independently verified this snapshot |
| **11693043** | 🏃 RUNNING | `beef_mmi...` — unrelated to this thread, another session's work, not independently verified this snapshot |
| **11672738_0..3** | ✅ **COMPLETED** (all 4 tasks, ~2:45 hrs each) | `beef_prewarm_mmi1x4_eme` — cache-prewarm for `mmi_1x4_eme.py` (step 6b), much faster than the ~13h/branch extrapolation from step 5's design tool. Real run **11678090 COMPLETED in 38s** (cache hit for all 4 branches): T=[0.1415, 0.1598, 0.1598, 0.1415], total T=0.6025 (IL=-2.20 dB) at λ=1.550 µm — exact mirror symmetry, consistent with step 5. `mmi_1x4_eme.py` is genuinely done |
| ~~11672739~~ → **11696233** | ❌→🏃 | `beef_prewarm_mmi1x4_hybrid` — **11672739 TIMED OUT at 24h** (all 4 tasks), only reached MPB band 14/20 on the FIRST of two body-adjacent interfaces (input↔body); true per-branch cost extrapolates to 60-80+ hrs at `N_MODES=20`. In `HybridSolver`, `interface_modes` IS the body's resolution for this cascade (not an independent interface-only knob — `EMESolver(num_modes=...)` for `body_solver` is an unused placeholder for the `solve_section_length_sweep()` call path), so there's no way to keep body accuracy at 20 while only shrinking the FDTD interface count. User-directed fix: reduced `N_MODES` 20→12 (preliminary, unvalidated — revisit at step 8), matching the same tradeoff `mmi_1x2_hybrid.py`/`mmi_2x2_hybrid.py` already accept at 8/14 modes. **Resubmitted as 11696233, `-t 48:00:00`**, running |

---

**★ TWENTY-FIFTH SNAPSHOT — 2026-08-04, verified against
`squeue -u jhiesener3`. Superseded by the snapshot above for job states;
still current for full context on the modal multiplexer's earlier jobs.**

| job | state | what |
|---|---|---|
| **11670540** `[0-3]` + **11678100** | ✅ **CLOSED** — mean ratio **0.9958** | `beef_modal_mux_engine_agreement` — `modal_multiplexer_validation.md` gap 5 (README.md check 1), first measurement ever. All 4 array tasks completed (3:07–5:33 hrs); merge **11670544 TIMED OUT** at 30 min (all 4 cached lengths loaded instantly, then hung on the ~20-min EME half with no margin left); re-budgeted to 90 min + thread-pinned, resubmitted as **11678100 — COMPLETED 30:37**. Result: Δβ_FDTD/Δβ_EME = 0.9958 (0.42% bias, monotone 0.9978→0.9948 across the O-band), worst off-diagonal/fit residual both 0.0000 — within 0.05–0.12 points of both MMIs' own engine bias (0.9946/0.9953). `calibrate_body_neff` NOT indicated. Check 1 CLOSED for this device |
| ~~11673180~~ → **11673526** | ✅ **COMPLETED** 4:51:09 | `beef_mode_spec_converge` — extended `MEEP_RESOLUTIONS` to include 80/100 (gap 2). 11673180 FAILED at res 20 (duplicate-band guard, fixed in the script, see detail section below). **Result: clean 14/16 guided at every resolution 20–100 (an improvement), but meep's absolute n_eff does NOT plateau by 100** — 60→80 looked converged (14/14 agree) but 80→100 breaks completely (0/14, shifts 2–3× larger). Full account in the validation doc's gap 2 and the detail section below |
| **11673498** | ✅ **COMPLETED** 11:34 | `beef_modal_mux_fdtd_SMOKE` — first-ever real run of `modal_multiplexer_fdtd.py`, at `MUX_FDTD_RES=10` (gap 3, README.md check 4). Pipeline ran end to end: both broadband input-mode solves completed, `check_passive()` held, real `SMatrix` saved. TE0 in → 56% total transmission, TE1 in → 1% (expected at this deliberately crude smoke resolution, not a device property). Real res=20 run not yet submitted |
| **11670282** | 🏃 RUNNING | `beef_hcg` — high contrast grating full-band production run, unrelated to this thread, unchanged this snapshot |
| **11667114_0..3** | ✅ **COMPLETED** (all 4 tasks, ~5:49-6:29 hrs each) | `beef_mmi1x4_body_sweep` — 1x4 MMI on 300nm Si3N4, step 5's `mm_length` design sweep. Merged: `argmin(spread)`=280.0 µm, but **user picked `mm_length=206.0 µm`** instead (footprint/throughput/balance tradeoff — `mmi_1x4_sin300_validation.md` step 5) and asked to drop `mmi_1x4_fdtd.py` (infeasible at this length). |
| **11672738_0..3** | 🏃 RUNNING (job array, 4 tasks) | `beef_prewarm_mmi1x4_eme` — temporary cache-prewarm for `mmi_1x4_eme.py`'s real deliverable run (step 6b), one task per output branch, `-t 20:00:00`/task. See detail section below. |
| **11672739_0..3** | 🏃 RUNNING (job array, 4 tasks) | `beef_prewarm_mmi1x4_hybrid` — temporary cache-prewarm for `mmi_1x4_hybrid.py`'s real deliverable run (step 6d), one task per output branch, `-t 24:00:00`/task (no timing precedent for this device's FDTD sections, generous margin). See detail section below. |

### `modal_multiplexer_hybrid.py` N_MODES changed 16→14 — a judgment call, not a job

Direct consequence of gap 2's resolution-sweep result: the item-12/14 guard
means `N_MODES=16` will raise on this device's very next run (14 reachable
at the shipped `RESOLUTION=40`, measured cleanly and consistently at every
tested resolution 20–100). Changed `N_MODES` 16→14, kept `RESOLUTION=40`
(no resolution 20–100 is demonstrably converged in absolute n_eff, so a
finer mesh buys no demonstrated accuracy; this device's own check 1
measured the engine ratio — not either engine's absolute value — at
`RESOLUTION=40` and found it excellent, 0.9958). **Made without explicit
user sign-off** — flagged loudly in the script's own comment and in
`modal_multiplexer_validation.md`'s gap 2 boxed note for review. Not
independently verified: whether 14 still carries several supermodes of
each parity at every interface, the constraint the script's own comment
states — predicted yes, not checked against a real run (gap 6's territory,
not attempted this session). Job 11688645 (`run_modal_multiplexer_fdtd.
batch`, res=20) submitted the same round, using the now-fixed
`CompoundGeometry.polygons()` path — unrelated script, no N_MODES
dependency, this is check 4's real (non-smoke) run.

### 11670540 `[0-3]` + 11670544 — `beef_modal_mux_engine_agreement`, this device resumed out of `deferred/`

At the user's direction, `modal_multiplexer_validation.md` moved out of
`deferred/` back to the top level of `validation_markdown/` today, and work
restarted with gap 5 (check 1, engine agreement) rather than picking up where
gap 6's exercise run left off (the 25 h/three-timeout campaign) — gap 5 is
unblocked, cheap relative to this device's other campaigns, and answers a
localized question (does the body's own uniform cross-section accumulate
phase the same way in both engines) independent of everything else open on
this device. See the validation doc's "Re-ordered 2026-08-04" note under
Order for the full reasoning. Gap 3's `CompoundGeometry.polygons()` library
gap was picked up the same day, once gap 5 and the gap-2 extension were both
already in flight — see the detail section below.

`examples/verification/python/verify_modal_multiplexer_engine_agreement.py`
is a direct adaptation of `verify_mmi_2x2_engine_agreement.py` (same
mechanism as `resolved/directional_coupler_validation.md` item 6): cut MMI
1's body (mm_length=16.92, gap 1 already established both of this device's
bodies share an identical cross-section) at 2/4/6/8 µm, fit
`Δβ_FDTD` from FDTD's accumulated phase via `dbeta_from_phases()`, compare
against `Δβ_EME` from the same cross-section's eigenvalues — at the shipped
`RESOLUTION = 40` and this device's own O-band grid, unlike the coupler/2×2
scripts' C-band. Shipped already in the array-per-length + merge shape from
the start, learning from the 2×2 sibling's own history (its first attempt,
11591249, timed out at 5h with 2/4 lengths done) rather than rediscovering
it. Deliberately bare `python`, not `srun --mpi=pmi2` — no rank guards in
this script, and both MMI siblings already proved a single meep process
finishes a problem this size inside an 8h budget.

**✅ CLOSED — mean ratio 0.9958.** All 4 array tasks completed cleanly
(3:07–5:33 hrs). The merge (**11670544**) TIMED OUT at its original 30-min
budget: all 4 cached lengths loaded instantly (confirmed cache hits), then
it hung on "running the EME half..." until the wall clock killed it —
despite a login-node pre-submission measurement of that exact call
finishing in 1237 s (~20.6 min), there wasn't enough margin left once SLURM
startup overhead was added. Re-budgeted `-t` 30→90, pinned
`OMP_NUM_THREADS=1`/`OPENBLAS_NUM_THREADS=1` to match the successful
login-node conditions, resubmitted as **11678100 — COMPLETED in 30:37**:

```
     wl   dbeta_fdtd    dbeta_eme     ratio
 1.2600      0.13674      0.13704    0.9978
 1.3100      0.14379      0.14443    0.9955
 1.3600      0.15115      0.15193    0.9948

mean ratio (FDTD/EME): 0.9958
```

Worst off-diagonal 0.0000, worst FDTD fit residual 0.0000. This lands
within 0.05–0.12 points of both MMIs' own engine-level bias (1×2: 0.9946,
2×2: 0.9953) — a third device confirming the same ~0.4–0.6% MPB/meep
discretization bias on an MMI-body cross-section, strengthening the
"engine-level, not device-specific" reading. `calibrate_body_neff` is NOT
indicated here, same conclusion as the 2×2 MMI. **This is the first of the
four checks to close for this device.**

### 11673526 (was 11673180) — gap 2 extension, and a real guard-interaction finding

`mode_spectrum_convergence.py`'s `MEEP_RESOLUTIONS` grew from `[20, 30, 40,
60]` to `[20, 30, 40, 60, 80, 100]` — gap 2's own flagged "cheap next step,"
answering whether meep/FDTD's mode spectrum plateaus by res 80/100 the way
meow/EME already does at 60→80. First submission (**11673180**) FAILED after
66 s, not from the extension itself but at res 20, the FIRST point in the
(unchanged) low end of the sweep: `FDTDSolver: this cross-section supplies
only 14 guided mode(s) of the 16 requested`. This script's original
successful run (job 11575516) was 2026-07-30; items 12/14
(`resolved/fdtd_vs_EME_mode_validation.md`) — which make the duplicate-band
guard raise the instant it reaches an unguided band, instead of padding —
were fixed the day after, 2026-07-31. Re-running the same script today hit
that guard for the first time, at exactly the resolution where this
cross-section has fewer than 16 guided modes. A real, useful confirmation
(the guard is live here, corroborating gap 2's own account from a second
angle), but the script's fixed `N_MODES=16` request across every resolution
no longer works unmodified. Fixed `_meep_spectrum()` to catch that specific
`RuntimeError`, parse the guard's own reported guided count, and retry at
exactly that count — verified the regex against the real captured error
text before resubmitting as **11673526** (same budget/account, unchanged
otherwise).

**✅ 11673526 COMPLETED, 4:51:09 — the retry-degrade fix worked at every
resolution.** Clean, consistent **14 guided meep modes at every resolution
20/30/40/60/80/100** — better than the pre-fix 13/13/13/14 record, since it
comes from the same robust walk at every point rather than depending on
where a duplicate happened to fall. But the actual question ("does it
plateau?") came back **NO**: self-convergence (|Δn_eff|≤0.005, of 14)
went 20→30: 1/14, 30→40: 0/14, 40→60: 1/14, **60→80: 14/14** (looked like a
plateau!), **80→100: 0/14** — every one of the 14 modes shifted by
0.0083–0.0127 between 80 and 100, roughly 2–3× the 60→80 step (max 0.0048),
the wrong direction for ordinary mesh convergence. Checked by hand against
the raw spectra, not read off the summary table: the shift is small,
uniform in sign and monotonically increasing across all 14 modes — a real
effect, not noise. No resolution from 20–100 is demonstrably converged in
absolute n_eff at this cross-section. Matched-resolution meep-vs-meow
agreement stayed good throughout (13–14/14) — the one residual is the
already-tracked near-degenerate-pair defect (items 11–13), now confirmed to
persist at res 80 too (previously only documented through res 60). Full
account, including why the "orphans" column reads empty despite this,
in the validation doc's gap 2.

### 11673498 — `beef_modal_mux_fdtd_SMOKE`, gap 3/check 4 unblocked

`CompoundGeometry.polygons(stub_length>0)` had raised `NotImplementedError`
since 2026-07-31, blocking `modal_multiplexer_fdtd.py` (this device's only
possible full-device FDTD reference, README.md check 4) entirely — the
stated reason was "no generic way to know per-port stub width for an
arbitrary segment type." That premise was wrong, not just unimplemented:
`_end_ports()` already answers the question for `to_gds_eme()`'s own stub
handling, by returning EVERY port at an end rather than picking a winner,
and a multi-port leaf geometry already draws one stub per port for its own
`polygons()`. Fixed by wiring `_end_ports()` into `polygons()` the same way
`to_gds_eme()` already uses it (`beef/geometry/compound.py`,
`notes.md` Geometry item 8 part 2). Confirmed against this exact device: the
mux's `CompoundGeometry` now builds the correct 4 stub polygons (1 input + 3
output, each at its own port's true width/position), and `to_meep_prisms()`
— the actual call `FDTDSolver.solve()` makes — now succeeds. 7 new tests,
`tests/test_geometry.py::TestCompoundGeometryPolygonsStub`; full fast suite
still green (618 passed, 19 skipped).

**This fixes geometry construction only — the script had never actually
been run.** Per this repo's own two-step discipline for expensive, never-run
scripts, submitted a smoke test first rather than the real run:
`run_modal_multiplexer_fdtd_SMOKE.batch`, `MUX_FDTD_RES=10` (half the
script's own already-reduced res=20 default, ~16× cheaper per its own res⁴
cost scaling), `-t 60`.

**✅ COMPLETED, 11:34 — the pipeline runs end to end, for the first time
ever.** Both broadband input-mode runs finished (58,215 timesteps for TE0),
`check_passive()` did not raise, a real `SMatrix` saved: TE0 in → top/mid/
bot T=0.0902/0.1430/0.3284 (56% total), TE1 in → 0.0072/0.0027/0.0025 (1%
total). Read the TE1 number as a smoke-resolution artifact, not a device
property — `res=10` gives the 0.7 µm multi-mode input waveguide only ~7
pixels across, not enough to launch a higher-order mode cleanly. One
pre-existing, unfixed hazard noted: like `mmi_1x2_fdtd.py`/
`mmi_2x2_fdtd.py`, this script has no `is_master()` gating, so all 24 MPI
ranks wrote the same output files independently (~24 duplicate "saved →"
lines) — harmless here (every rank computes identical results), not new,
not in scope this session.

**The real run — submitted as 11688645, ❌ TIMED OUT at 4h, a real cost
finding, resubmitted at 24h as 11693107.** Still mid-way through just the
FIRST input mode (TE0) when killed: 416,806 timesteps, t=10,420.2, flat
~0.0336 s/step with no sign of approaching its decay stop — already 3.4×
further in simulated time than the smoke test's ENTIRE TE0 run needed to
decay (t≈2910–3217 at res=10), which is more than the ~2× step-count CFL
scaling alone would predict for the same physical duration. Not merely a
budgeting miss: something about decay behavior itself appears
resolution-dependent here, mechanism unconfirmed (open question, not
investigated further). No caching in this script by design, so all 96
core-hours were a total loss. Presented to the user with 4 options (bigger
budget / loosen `decay_tol` / refactor into an array job first / hold off);
**chose bigger budget, unchanged otherwise** — `-t` 240→1440 (24h,
matching this repo's largest existing precedent) in
`run_modal_multiplexer_fdtd.batch`, resubmitted as **11693107**.

**11693107 also did not converge, and this time the user found a real
external comparison point instead of authorizing another blind wait.**
Killed at 18h53m, still stuck in TE0, t=49,375.9 — 4.7× further than
11688645's own stopping point, same flat step rate, no sign of leveling
off. The user pointed at `04_weary_paw/` — an archived output tree (config,
logs, geometries; no driver code) from an actual seeded topology-
optimization run of this exact device (`gds_name: "mux_MMI_2"`,
`res=40`/`dr_res=80` — HIGHER resolution than this script). Its
`logs/master.log` shows a full forward+adjoint gradient cycle (both TE0/TE1
groups, plus constraints — much more work than this script's forward-only
sweep) completing in **~70–80 minutes per iteration**, 56 iterations over
~67.3h total. This repo's own simpler, lower-resolution run had already
spent ~19h without finishing one mode — a ~15–20× gap too large to pin on
core count alone (no core/node count was even recoverable from the
archived logs to confirm), but core count is the one lever definitively
under this repo's control, so at the user's direction: killed 11693107,
rescaled `run_modal_multiplexer_fdtd.batch` 1 node/24 ranks → **4 nodes/96
ranks** (matching `modal_multiplexer_hybrid.py`'s own scale for this
device's FDTD sections), resubmitted as **11707932**. Confirmed at launch
(`Using MPI version 4.1, 96 processes` / `Splitting into 96 chunks by
voxels`) — the full 4-node allocation is actually being used. `decay_tol`
(the more likely dominant lever, per the archived campaign's own much
shorter per-iteration time despite doing more work at higher resolution)
was NOT changed — still a methodology decision needing explicit sign-off.
Not yet landed.

**Live check 2026-08-06, ~1h02m elapsed — measured, not projected from the
worry above.** `reports/11707932.out`: step 387,803, t=9695.08, steady
~0.0089–0.0090 s/step. Average sim-time rate since launch ≈ **156.1
units/min**, vs. 11693107's own flat whole-run rate ≈ **43.6 units/min** (1
node/24 ranks) — a **~3.6× speedup for 4× the ranks** (~88–94% parallel
efficiency), better than "parallel efficiency drops well before 4× more
ranks buys 4× speedup" anticipated a paragraph above. At this rate, reaching
t=49,375.9 (11693107's own kill point, still not decayed there) projects to
≈5.3h elapsed, leaving ≈18.7h of the 24h budget as margin beyond that — a
rate projection, not a completion; still no confirmation `stop_when_
dft_decayed(tol=1e-11)` actually triggers before budget runs out.

Also went back to `04_weary_paw/` for the rank-count question flagged above
as unrecoverable from its own logs. Re-confirmed: zero hits grepping
`logs/master.log`/`logs/group{0,1}.log` for process/mpi/chunk/rank/core/
node/thread/parallel — the custom logger never captured meep's own startup
banner. Found indirect corroboration instead: the only two surviving SLURM
scripts from this repo's older Makefile-driven framework generation
(`archive/FECES_Legacy/batch/run.batch`, `rum.batch` —
`srun --mpi=pmi2`-launched, matching `04_weary_paw`'s implied launcher)
both request `-N 4 --ntasks-per-node=24` (96 ranks) for their own full 3-D
campaigns — but for sibling devices (`mmi_2x2_bend`/`mmi_2x2_gds`), not the
mux, and the driver modules `04_weary_paw`'s own logs name (`Mux_MMI.py`,
`engine.py`, `utils_s_mux.py`) match zero files anywhere in
`archive/FECES_Legacy` — `masters/mmi_mux_master.py`, the only mux-related
survivor, is a different, later rewrite (`class MMI_MUX`, direct `meow`/
`meep` imports), not what produced `04_weary_paw`. Recorded as circumstantial
corroboration of the framework generation's standard job scale, not proof of
what the archived campaign itself used — full account in
`modal_multiplexer_validation.md` gap 3.

**✅ 2026-08-06, ~9h06m elapsed — TE0 fully decayed, first time ever in this
campaign.** `run 0 finished at t = 77197.45 (3087898 timesteps)` — far past
where either single-node attempt died (11688645 at t=10,420/4h;
11693107 at t=49,375.9/18h53m). Wall time for TE0 alone ≈ **7.95h**
(7,048 progress lines × meep's ~4.06 s/line default report cadence, measured
directly from this run's own line-count-vs-elapsed-time, not assumed). This
also explains the odd-looking near-stall between the two live checks above:
not a stall — the per-`mp.Simulation` step counter resets to 0 at the start
of each input mode, and TE0's finish plus ~55 eigenmode-coefficient MPB
solves for its flux monitors landed in that window. Now into
**input_mode=2 (TE1)**, ~1h15m in, same healthy ~0.0087 s/step. No S-matrix
yet (script prints/saves only after both modes complete). Worth noting for
the `decay_tol` debate above: TE0 hit `tol=1e-11` in well under half the 24h
budget once given 4× the ranks, real evidence that core count — not
necessarily `decay_tol` — was the dominant lever here, at least for this
mode; TE1's own decay time is still unmeasured. Full account in
`modal_multiplexer_validation.md` gap 3.

### 11667114 — a real bug found and fixed first, not a routine sweep

The design-tool script for this sweep went through two versions before this
run. v1 (all-EME via `EMESolver.solve_section()` on the taper bands) was
**wrong, not imprecise** — smoke test 11660322 returned T=1.0000 for every
output/length/wavelength (`total T=4.0000`, violating energy conservation).
Root cause, confirmed both empirically (job 11661894: `solve_section()` on
this device's own output band, 6 modes, returned exactly `np.eye(6)`) and by
reading meow's source: a single `Ls=[length]` cell samples one cross-section
and shares one mode basis between both its "ports", so `compute_s_matrix()`
has no adjacent interface to mix modes across — the result is the identity
matrix **by construction**, regardless of the geometry's real complexity.
`solve_section()` is only valid on genuinely uniform regions (as
`solve_section_length_sweep()`'s own docstring already required for the
body) — exactly why `mmi_1x2_hybrid.py`/`mmi_2x2_hybrid.py` use FDTD, not EME
`solve_section()`, on every taper band. v2 rewrote the script on
`mmi_1x2_eme.py`'s actual all-EME pattern (`EMESolver.solve_length_sweep()`,
tapers genuinely subdivided into real cells via `num_cells`) — smoke-tested
working (job 11662670: T now varies with length, shows the expected
left-right mirror symmetry). Smoke timing also surfaced a real cost problem
(~1 day if all 4 output branches ran sequentially in one job), fixed by
splitting into this 4-task job array instead. Full chain of evidence in
`mmi_1x4_sin300_validation.md` step 5.

### 11672738 `[0-3]` + 11672739 `[0-3]` — cache-prewarm for the two real 1x4 example scripts

`mmi_1x4_eme.py`/`mmi_1x4_hybrid.py` (step 6, `MM_LENGTH=206.0` now filled
in, `mmi_1x4_fdtd.py` deleted) each loop over all 4 output branches
SEQUENTIALLY in one process at full `N_WL=11` — the same pattern that made
the step-5 design tool take ~1 day naively, except worse here (design tool
used a reduced `N_WL=5`; this projects to 2+ days for the EME script alone
by extrapolation). Rather than permanently restructure the final example
scripts, temporary scripts `_prewarm_mmi_1x4_eme.py`/
`_prewarm_mmi_1x4_hybrid.py` (in `examples/verification/python/`, NOT
`examples/mmi_1x4/`) solve one branch each (`MMI1X4_OUTPUT_IDX` env var)
into the same persistent `save_dir` the real scripts read. Once all 8 tasks
land, the real, unmodified example scripts run normally and hit the warm
cache for every solve call — near-instant, producing the actual
`smatrix.pkl`/plots. Hybrid's per-branch input/body sections will be
recomputed redundantly a few times (4 parallel jobs racing, no
cross-process lock on which one "wins" the shared cache entry first) — a
compute-time inefficiency, not a correctness issue; this codebase's
`save_dir` cache is already used safely by other parallel job patterns.
Prewarm scripts + their batch files deleted once the real scripts have run.

---

**TWENTY-THIRD SNAPSHOT — 2026-08-03 (even later), verified against
`squeue -u jhiesener3`.**

| job | state | what |
|---|---|---|
| **11655592_0** | ✅ **COMPLETED** 06:38:10 | `beef_dc_L151` hybrid task — closes out the 15.1 µm coupler re-run. `Plot length: user-set 15.100000 µm -> sampled length 15.100000 µm (exact)` confirms the splice fix holds under real load. **Cross T=0.4907, Through T=0.5075, ratio=0.492** — within 0.3% of the FDTD task's ratio=0.492 (Cross T=0.4920/Through T=0.5075), consistent with this device's documented ~0.1% engine-agreement level |
| **11660092** | 🏃 RUNNING | `beef_mmi1x4_sin_xsec` — nitride cross-section sweep (steps 4a/4b), corrected candidate ranges and `compute_modes_tidy3d` backend (see `mmi_1x4_sin300_validation.md`) |
| **11658729_3** | 🏃 RUNNING, ~1:54 | `beef_ring_gap3_fdtd_res` — unrelated, unchanged this snapshot |

### 11660092 — what it's for, and what came before it

Two prior submissions for the same goal were **cancelled**, not completed:

- **11659619** — the first real cross-section sweep, using `compute_modes_meep`
  (this repo's usual mode-solve backend). Cancelled after it reported ZERO
  guided modes at every tested width/sz — later confirmed to be a real,
  general limitation of that backend for weakly-index-contrasted structures
  (silicon nitride's core/clad contrast is ~5.3x weaker than Si/SiO2's), not
  real physics. Full chain of evidence in the validation doc's step 4 section.
- **11659878** — a 4-case diagnostic to isolate the cause. Cancelled once the
  root cause was found via a faster, more direct test (three explicit initial
  k-guesses fed straight into meep's `get_eigenmode`, all three converging to
  the identical wrong answer) before the diagnostic's slower cases finished.

**11660092** uses `compute_modes_tidy3d` (meow's other, FDE-based mode
solver) instead, cross-validated against `compute_modes_meep` on a silicon
control first (agrees ~1-3%, consistent with this repo's other documented
engine gaps). Candidate ranges also widened based on real smoke-test numbers
with the fixed backend, since the original ranges didn't bracket either the
sz convergence point or the single-mode width cutoff.

---

| job | state | what |
|---|---|---|
| **11655592_[0-2]** | ⏳ PENDING, 2 nodes × 24, `-t 720` | `beef_dc_L151` — all three coupler examples at **15.1 µm**. Task 0 = `directional_coupler_hybrid`, 1 = `directional_coupler_mode`, 2 = `directional_coupler_fdtd` |
| **11655405_[0-2]** | ❌ CANCELLED ~1 min in | `beef_dc_L1513` — same run at 15.13 µm. Superseded by 11655592; script deleted |

### What 11655592 is for

Re-run of all three directional coupler examples with the hybrid and the
full-device FDTD **at the same coupling length, 15.1 µm** — user request,
2026-08-03. Script `batch/run_directional_coupler_L151.batch`, account
`gts-sr108-ece`.

**Why 15.1 and not 15.13.** 15.13 is the FROZEN-era 3 dB length (item 5's
ground-truth crossing, 15.129 µm). Dispersion defaulted 2026-08-01 and moved
L₅₀ down about 0.05 µm: **15.089958 µm** measured for the hybrid (job 11609658)
and **~15.074 µm** inferred for the full-device FDTD (job 11609659 — band-mean
ratio 0.500675 at 15.089958, over the hybrid's local slope 0.0418/µm). 15.1 is
the round tenth nearest both. ⚠ It is a *rounded convenience value*, not a
measured crossing for either engine — the item-5 sweep has still never been run
dispersively. Do not quote 15.1 as a measured 3 dB length.

**The two examples were NOT at the same length before this fix, despite both
source files reading `15.13`.** `directional_coupler_fdtd.py` builds its device
at a true 15.13. `directional_coupler_hybrid.py` only ever `argmin`'d its pinned
`PLOT_COUPLING_LENGTH` against `np.linspace(2*BODY_STUB, 60.0, 200)` — a
0.298995 µm grid that contains neither 15.13 nor 15.1 — and so **actually
simulated 15.150754 µm, 0.0208 µm away** (measured: `argmin` index 49 on that
grid). The old `.2f` print rounded this to `"15.13 -> 15.15"`, which is why it
survived. Fixed 2026-08-03 by splicing the pinned length into the grid (201 pts)
and asserting the hit is exact; `L_mid` is still taken from the base grid so the
FDTD section cache keys do not move (confirmed unchanged at 30.399497 µm).

Verified pre-submit by executing each script's geometry header: both now build
`dc.coupling_length = 15.1`, hybrid `argmin` sample `= 15.1` exactly.

⚠ Even at 15.1 this is an **engine-vs-engine comparison at one identical
length**, not a 50:50 search. The two engines' crossings are ~0.1% apart
(15.129/15.151 frozen, 15.074/15.090 dispersive) and that gap is a *result* this
device measures. Expect a ratio near but not exactly 0.5 from both; do not read
the offset as disagreement.

This script does **not** `rm -rf` the output dirs, unlike
`run_directional_coupler_examples.batch`: they hold `*_FROZEN_pre20260801.*`
baselines that cannot be regenerated. Copies parked in
`example_output/_frozen_baseline_backup/` as well. `-t 720`, not that script's
300 — 300 is what killed cache-cold hybrid job 11592462_0 (parity item P7).

---

**TWENTY-FIRST SNAPSHOT — 2026-08-03 (later), verified against
`squeue -u jhiesener3`. Superseded; kept unedited.**

| job | state | what |
|---|---|---|
| **11653509_[0-3]** | ⏳ PENDING (`Resources`), 2 nodes × 24, `-t 600` | `beef_ring_gap3_fdtd_res` — **gap 3's missing half: FDTD-side resolution convergence** (check 3). Array sweeps `RING_FDTD_RES` = **20 / 30 / 40 / 60** |

### What 11653509 is for, and how to read it

Check 3 has been **half** measured since 2026-08-01: the mode side converges at
`MODE_RESOLUTION=100` (job 11616610), but **`FDTD_RESOLUTION=40` has never been
measured at all**, so every resonance position this device reports rests on an
unchecked assumption. This array is that measurement.

**Read the comb POSITION, not the transmission magnitudes.** Gap 3's mode-side
result moved `n_eff` by only 0.028 and slid the comb an eighth of an FSR;
magnitudes settling proves nothing about resonance position. Convergence means
the comb stops moving.

⚠ **Task `_2` (res 40) is a deliberate control, not a spare.** It reproduces the
production config, so it must come back a **pure cache hit**. If it instead runs
a simulation, something has moved in the fingerprint and the sweep is not
comparable to 11634900_1.

Outputs land in
`example_output/.../ring_resonator_hybrid_cascaded/res_sweep/res_<N>/`, one
subfolder per task — added today via a `RING_FDTD_RES` guard, because an
unguarded sweep would have had all four tasks overwrite the same
`ring_data.npz`. That is not hypothetical: it is how a smoke run destroyed the
11634900_1 arrays earlier the same day.

Expected wall clock, from res^4 scaling against res 40's 44:18: res 20 ~3 min,
res 30 ~14 min, res 40 instant (cache), res 60 **~3.7 h**.

---

**TWENTIETH SNAPSHOT — 2026-08-03, verified against `squeue -u jhiesener3`
(empty) and `sacct -u jhiesener3 -S 2026-08-02`. Superseded; kept unedited.**

## NOTHING RUNNING. Queue is empty.

| job | state | what |
|---|---|---|
| **11634859** | ✅ **COMPLETED 3:56:16**, exit `0:0`, ended 2026-08-02T23:10:59 | `beef_ring_xc_clampfix`. Ran to the end well inside `-t 480`. **Right section and composed cross-check both landed** — see below |
| **11634900_0** | ✅ COMPLETED 58:06 | `ring_resonator_hybrid`, 3-D. Mean FSR **12.975 nm (0.8%)** |
| **11634900_1** | ✅ COMPLETED 44:18 | `ring_resonator_hybrid_cascaded`, 3-D. Mean FSR 13.350 nm (3.7%) |

### ✅ 11634859 finished — the whole cross-check passes, not just the left section

The nineteenth snapshot recorded this job's left-section numbers while it was
still at 1:57:04. It has since completed. New material:

```
right bus<-iface:  summed |T|^2 over seam modes = 0.9997 0.9992 0.9970 0.9953 0.9933
right ring<-iface: summed |T|^2 over seam modes = 0.9989 0.9985 0.9989 0.9986 0.9980
composed ports: ('bus_L', 'ring_L', 'bus_R', 'ring_R')
bus-side power sum |t|^2+|k1|^2: min 1.0050 max 1.0067   (lossless coupler -> ~1)
reciprocity max ||k1|-|k2||: 0.0021   (-> 0)
physics checks PASS (power sum ~1, reciprocity holds)
worst |E_legacy - E_network| = 1.735e-18
L_outer = 7.0000 um
VERDICT: routes AGREE to machine precision
```

**The right section is clean too** (0.993–1.000, no clamp pathology), so the
fix holds on both sides at production config. Power sum and reciprocity pass on
the *composed* 4-port. Per the standing entry, the `AGREE to machine precision`
verdict is the weakest line here — it printed identically pre-fix with ~88% of
the power gone.

✅ **`crosscheck.npz` is restored.** This job re-wrote
`example_output/verification/ring_network_crosscheck/crosscheck.npz` at 23:10,
replacing the array the login-node smoke runs destroyed. Contents are now
post-clamp-fix and valid. The `OUT`-redirect warning still stands for future
login-node runs.

⚠ **This job printed `L_outer = 7.0000 um`, which corrects gap 3.** Gap 3's
convergence script hardcodes `L_OUTER = 2*2.0 + 10.0` = 14.0 µm on a comment
asserting `L_C=10.0` is "the shipped scale"; the shipped `L_C` is **4.0**, and
both device scripts compute `L_outer = 2*(L_SIDE - S_MON) + L_FAR` = **7.0 µm**.
Gap 3's lever arm is exactly **2× too large** and its Δφ column halves. Details
and the corrected table in `ring_resonator_validation.md` gap 3.

### The half-FSR route offset — candidate space collapsed by measurement

The nineteenth snapshot left this undiagnosed with two named candidates. **Both
are now eliminated, by source comparison plus the saved arrays — no new FDTD.**
`ring_data.npz` retains all coefficients as `complex128`, so this was
post-processing. Full write-up in `ring_resonator_validation.md`. Short version:

- `L_outer`, `phi`, `phi_bend_correction`, `BEND_D_IN/OUT`, `neff_fit` — the
  whole analytic phase chain is **byte-identical code fed identical constants**
  in both scripts, and `neff` fit residual prints identically (7.29e-06).
- `b_raw` is **identical to all printed digits in magnitude AND phase**
  (`|b|=0.998828`, `arg=-1.275855` rad, difference exactly `0.000000`).
- So the offset enters through the **coupler coefficients**, not the ring
  phase: `arg(t_ring)` differs by **0.883π**, `arg(t)` by 0.404π.

**Measured. Not yet attributed to a defect** — see the validation doc for what
that does and does not license.

---

**NINETEENTH SNAPSHOT — 2026-08-02 ~20:16, verified against `squeue -u
jhiesener3` and `sacct -X`. Superseded; kept unedited.**

| job | state | what |
|---|---|---|
| **11634859** | ⏳ RUNNING **1:57:04** @ 21:11 (`-t 480`) | `beef_ring_xc_clampfix`. **Left-section numbers in and they close item 2's question** — see below |
| **11634900_0** | ✅ COMPLETED **58:06** | `ring_resonator_hybrid`, **3-D**. **Mean FSR 12.975 nm, 0.8% from the analytic 12.87.** Reciprocity `\|k1\|=0.2530 vs \|k2\|=0.2530`, independently normalized |
| **11634900_1** | ✅ COMPLETED **44:18** @ 20:13 | `ring_resonator_hybrid_cascaded`, 3-D. Mean FSR 13.350 nm (3.7%) |

**Order item 1 is CLOSED** — first trustworthy resonance positions this device
has produced. **But the two routes disagree by ~half an FSR**, undiagnosed;
full write-up in `ring_resonator_validation.md`'s new
`## Order item 1 — the re-run` section. Do not quote one route's comb without
saying which.

### ✅ 11634859 — the clamp fix is confirmed at production config

```
left iface<-bus:  1.0056 1.0073 1.0084 1.0098 1.0117
left iface<-ring: 0.9903 0.9906 0.9907 0.9916 0.9920
left ring<-bus:   0.0000 0.0000 0.0000 0.0000 0.0000
left TOTAL out of bus = 1.0056 1.0073 1.0084 1.0098 1.0117
left deficit vs 1.0   = -0.0056 -0.0073 -0.0084 -0.0098 -0.0117
```

Against **11627401's +0.7065 … +0.9307**, at the same 3-D res-20 config. The
~80% deficit that item 2 chased for two days was
`RingResonatorCouplerSection._section_polygons()` drawing bare cladding where
the caller's stub padding should have continued the bus. Right section and
composition still pending.

**Runtime is itself consistent with the fix:** 1:01+ against pre-fix runs of
7:04 and 8:21. Pre-fix the source radiated straight into the PML and
`stop_when_dft_decayed` tripped early — fast *because* the power was
vanishing. Suggestive, not proof.

### ✅ 11634900_1 — cascaded ring, 3-D, first trustworthy resonance comb

```
bus-excited path : |t|^2+|k1|^2 = 0.9998
ring-excited path: |t_ring|^2+|k2|^2 = 0.9928
reciprocity: max|k1-k2| = 1.78e-16   (symmetric BY CONSTRUCTION -- weak evidence)
supermode neff fit residual: even=1.45e-05, odd=7.30e-07
Dense-grid assembly (2001 points): 8 resonance dip(s) found.
Measured mean FSR: 13.350 nm (predicted 12.87)
Results at 1.5500 um: |t|=0.9713 |k1|=0.2374 |k2|=0.2374 |t_ring|=0.9677 |b|=0.9988
```

FSR is **3.7%** off prediction, against **14%** in the 2-D `RING_TEST_MODE`
smoke — consistent with gap 3's measurement that res 20 sits ~75% of an FSR
from res 40.

⚠ **Q is grid-limited, and this run makes that visible rather than suspected.**
Every FWHM is an exact multiple of the 2001-point grid spacing (0.05 nm):

| FWHM (nm) | 0.1500 | 0.2000 | 0.2500 | 0.3000 |
|---|---|---|---|---|
| grid steps | 3 | 4 | 5 | 6 |

So the reported Q (5322–10193) is set by **`N_WL_DENSE`**, a knob independent
of both `FDTD_RESOLUTION` and `MODE_RESOLUTION`. In the 501-point smoke every
FWHM was pinned at exactly 0.4000 nm. **If Q is wanted, raise `N_WL_DENSE`** —
it is post-processing on an interpolated grid and costs no new FDTD.

---

**EIGHTEENTH SNAPSHOT — 2026-08-02 ~19:15, verified against `squeue -u
jhiesener3`. Superseded; kept unedited.**

| job | state | what |
|---|---|---|
| **11634859** | ⏳ PENDING @ 19:15 (`cpu-small`, reason `None` — not a billing hold) | `beef_ring_xc_clampfix` — `verify_ring_network_crosscheck.py` re-run after the `_section_polygons` clamp fix. 1 node × 24, `-t 480`, `--account=gts-sr108-ece` |

**What it is expected to answer.** Every bus-excited number in
`ring_resonator_validation.md` item 2 was launched from a source plane sitting
on bare cladding (see that file's `🛠 DEFECT FOUND AND FIXED` block). A
controlled login-node A/B at 2-D res 8 already puts `TOTAL out of bus` at
**0.2190 → 1.0087** across the fix alone. This job is the 3-D res-20
production confirmation of that, at the config the recorded numbers were taken
at. Expect `iface<-bus` ≈ 1 and the +0.71…+0.93 deficit to vanish.

**When it lands:** fold into `ring_resonator_validation.md` item 2, replacing
the 11627160/11627401 numbers rather than appending to them — those runs are
void, not merely superseded.

⚠ The script `np.savez`es to
`example_output/verification/ring_network_crosscheck/crosscheck.npz`
unconditionally. It currently holds **pre-fix-clamp 2-D res-8 A/B data** from a
login-node run, not a production result. This job overwrites it correctly.

---

**SEVENTEENTH SNAPSHOT — 2026-08-02 ~13:10, verified against `squeue -u
jhiesener3` (empty) and `sacct -X -S 2026-08-02`. Superseded; kept unedited.**

**THE QUEUE IS EMPTY. Nothing is running, nothing is pending.** Every job
this project submitted today has terminated.

| job | final state | answer |
|---|---|---|
| **11627401** | ✅ COMPLETED **8:21** @ 12:00 | `beef_ring_xcC_busbus` — the `bus<-bus` control. Reflection **0.1–2%**, forward **0.049–0.326**. Both surviving readings from item 2 stay alive; ~80% still unaccounted. Folded into `ring_resonator_validation.md` item 2 |
| **11627585** | ❌ **FAILED** (exit 1) **3:21** @ 12:11 | `beef_ring_xcD_launch` — the straight-bus launch control. Died in geometry, **before any FDTD physics**: `ValueError: Cut/stub at x=-3.5000 falls strictly inside the near-bend arc` |
| **11628207** | ❌ **FAILED** (exit 1) **3:20** @ 13:05 | `beef_ring_xcD_launch2` — same control after an attempted fix. **Failed identically**, `x=-6.5000`. The fix was wrong |

### ⛔ The launch control has never run, and the script as committed will fail a third time

**Measured 2026-08-02, not argued.** Both failures are the same defect, and
the correction applied between them did not remove it.

`RingResonatorCouplerSection._check_cut` (`beef/geometry/ring_coupler.py:168`)
rejects any cut **strictly inside** the near-bend arc, which for this device
spans `(stub_x, hc) = (-7.0, -2.0)`. The quantity it validates is not `x1`,
and it is not `x1 + stub_right` either:

```
FDTDSolver.solve_section()  beef/solvers/fdtd.py:1978
    prisms = geometry.section_prisms(
        x0, x1, si_med,
        stub_left=stub_left + self.dpml, stub_right=stub_right + self.dpml, ...)
```

**`solve_section` pads both stubs by `dpml` before the geometry ever sees
them** — deliberately, so the drawn waveguide continues through the full PML
instead of stopping at its inner edge. The cut actually tested is therefore

```
hi = X_STRAIGHT + stub_right + dpml = X_STRAIGHT + 1.0 + 0.5
```

Evaluated directly against the geometry class (conda `jmh`, login node, pure
geometry — no meep, no time-stepping):

| `X_STRAIGHT` | cut `hi` | guard | matches |
|---|---|---|---|
| −5.0 | −3.50 | ILLEGAL | job 11627585's error string exactly |
| −7.5 | −6.00 | ILLEGAL | |
| **−8.0 (committed default)** | **−6.50** | **ILLEGAL** | job 11628207's error string exactly |
| −8.5 | −7.00 | LEGAL | lands on `stub_x`, the boundary the guard allows |
| −9.0 | −7.50 | LEGAL | |

The reproduction is exact in both directions — the two failed jobs' error
strings are recovered from the geometry alone, so this is the whole mechanism
and not a coincidence.

**Why the second attempt failed.** The fix comment now standing at
`verify_ring_network_crosscheck.py:56-67` reasons: *"section_polygons()
validates x1 + stub_right … With common's stub_right=1.0, x1=-8.0 puts the
stub end exactly on stub_x=-7.0 — the boundary, which the guard allows."*
That arithmetic is right and its premise is wrong: it omits the `+ dpml`
term, so the real cut was −6.5, not −7.0. **The comment is still in the file
and still wrong**, and `X_STRAIGHT` still defaults to −8.0. Resubmitting the
batch script unchanged buys a third identical 3-minute failure.

**Cost so far:** 2 × ~3:20 wall on 24 cores. Cheap in absolute terms, but both
failures are pure geometry — reachable in under a second on a login node
without allocating a node at all.

**The open design call, which is why this is not patched here.** `X_STRAIGHT
= -8.5` is the largest legal value, and it leaves `[-9.0, -8.5]` — **0.5 µm**
of straight bus. That may be too short to be the control item 2 asked for
(the question is whether a clipped launch *sheds over a propagation
distance*). Dropping `stub_right` to 0.0 for this one call instead allows
`X_STRAIGHT ≤ -7.5` and a 1.5 µm span, but changes the port padding away from
what every other measurement in the script used. **Pick one deliberately** —
this is exactly the "sweep the wrong thing" trap the coupler's item 11 and the
MMI's gap 2 both recorded.

---

**SIXTEENTH SNAPSHOT — 2026-08-02 ~11:40, verified against `squeue -u
jhiesener3` and `sacct -X`.** Superseded by the seventeenth above; its two
live jobs have both since landed and 11627160 is unchanged. The fifteenth is
kept below it, unedited, for the two lessons it carries.

**Both jobs the fifteenth snapshot listed as live have landed, successfully.**

| job | final state | answer |
|---|---|---|
| **11609658** | ✅ COMPLETED **10:25:31** @ 01:42 (of 12 h) | `beef_dc_hybrid_dispersive`. Dispersive L₅₀ **15.089958 µm** (grid-snapped 15.150754, 0.299 µm grid). The 12 h budget was right — the 3 h that killed 11608703 was never going to be enough |
| **11609659** | ✅ COMPLETED **1:37:19** @ 03:20 | `beef_dc_fdtd_at_hybrid_l50`, `afterok` fired normally. At L=15.089958 the **band-mean** split ratio is **0.500675** |
| **11627160** | ✅ COMPLETED **7:04** @ 11:32 | `beef_ring_xcB_ringbus` — the ring's `ring<-bus` measurement. **Hypothesis refuted:** 0.0000–0.0004. See `ring_resonator_validation.md` item 2 |

### The coupler result, because it is easy to misread

The log's headline line reports **λ=1.550 only** (ratio 0.4919). The constant
`COUPLING_LENGTH` is defined on the **band mean**, which is 0.500675 — 0.13%
off a true 50:50. Reading the centre-λ line as the answer inverts the
conclusion. Computed from `smatrix_fdtd.pkl` directly.

| | frozen | dispersive |
|---|---|---|
| hybrid band-mean L₅₀ | 15.150754 | **15.089958** |
| FDTD band-mean L₅₀ | 15.129 | **~15.074** |

Dispersion moves **both** engines down ~0.05 µm and they stay ~0.1% apart.
**The hybrid-vs-FDTD gap survives the material model** — it is not a dispersion
artifact. This closes the `COUPLING_LENGTH` contradiction (2026-08-01 context
item 3): **15.13 was correct all along**, and the defect was the stale comment
block above it claiming the constant was pinned to the hybrid's 15.15.
Corrected 2026-08-02; value untouched. The ~15.074 figure rests on one length
plus the *hybrid's* slope, so it is **not** adopted as a constant — regenerate
item 5's ground-truth sweep dispersively first.

### ⚠ Not tracked here, and it should have been

**`11609123_7` `elevator_corner_round`** — running since 2026-08-01 13:21,
22 h+ elapsed of a 48 h limit, res=100, the last of an 8-point sweep. It lives
in `/storage/scratch1/8/jhiesener3/esc/`, not in BEEF, which is why the
2026-08-01 23:42 context save asserted "nothing else is running." **This file
is titled "Active SLURM jobs", not "active BEEF jobs".** It is also on
**`gts-sr108-gedc_main`** — the unfunded account, one of the 33 scripts still
pointing at it.

---

**★ FIFTEENTH SNAPSHOT — 2026-08-01 ~18:10, verified against `squeue -u
jhiesener3` and `sacct -X`.** Superseded by the sixteenth above.

**Every BEEF job the fourteenth snapshot tracked has landed.** All results are
folded into the owning device docs; this snapshot records the states and the
two things worth carrying forward.

| job | final state | answer |
|---|---|---|
| **11608875** `[0-3]` + **11609686** | ✅ all COMPLETED (55:15 / 41:22 / 49:33 / **2:23:09**), merge 25:45 @ 15:34 | **P1 CLOSED.** 1×2 gap 4, first ever: mean ratio **0.9946**, 0.54% monotone bias, worst residual 3.01e-06 rad |
| **11609215_1 / _2** + **11609687** | ✅ 13:44 / 2:00:01, comparison 0:12 @ 15:27 | 1×2 gap 5 dispersive: **9.11% → 10.18%** — marginally *worse*, mean offset unmoved |
| **11609211_1** + **11609680** + **11609681** | ✅ 1:51:57 / **3:50:41** / 0:15 @ 17:59 | **the P4 test → mechanism REFUTED.** 2×2 gap 5 dispersive **0.82% → 0.76%**, sign flip survives |
| **11608690** | ✅ COMPLETED 1:28:14 | **P2 CLOSED** — first calibrated hybrid to survive since the `beta_from_phases` fix |
| **11608692** | ✅ COMPLETED 00:07, **number contaminated** | gap 5 calibrated — valid number (7.30%) recovered offline; see the race below |
| **11608941** / **11608940** | ✅ 1:03:21 / 07:13 | gaps 2 and 3 dispersive — **byte-identical**, as predicted |
| **11608703** | ❌ **TIMEOUT 3:00:29** @ 15:17 | coupler L₅₀ re-verification — **did not finish**, see below |

~~**Two live jobs remain**~~ — **both landed 2026-08-02, see the sixteenth
snapshot at the top.** Left here as written:

| job | state | answer |
|---|---|---|
| **11609658** | ~~🟢 RUNNING 2:55 of **12:00**~~ ✅ 10:25:31 | `beef_dc_hybrid_dispersive` — the 11608703 resubmission, dispersive and with 4× the wall clock |
| **11609659** | ~~⏸ PENDING `afterok:11609658`~~ ✅ 1:37:19 | `beef_dc_fdtd_at_hybrid_l50`, 5 h |

### ❌ 11608703 TIMED OUT — two separate lessons

`beef_dc_hybrid_item14` hit its **3 h** limit at 15:17 with no output written.
Two things went wrong, and only one of them was the wall clock:

1. **Budget.** 3 h was sized against frozen-era coupler runs. The dispersive
   1×2 came in at **1.53×** its frozen cost (11609215_1) and the dispersive
   gap-4 lengths at ~2.4×; the coupler is the longest device in the repo. The
   resubmission carries **12 h**.
2. **⚠ It ran under `gts-sr108-gedc_main`, not `gts-sr108-ece`.** That is the
   account whose balance problems already blocked jobs on 2026-07-30
   (`AssocGrpBillingMinutes`). Everything submitted after it — 11609658,
   11609659, 11609686 — is explicitly on **`gts-sr108-ece`**. Check
   `#SBATCH -A` when reviving an older batch script; the timeout made this one
   visible, but a script sitting on the wrong account is silent until it is not.

**The P4 test was rescued from the same failure mode.** 11609211_2 was
cancelled 46 min in on a live cost projection rather than left to time out,
and its resubmission 11609680 took **3:50:41** — well past the original
budget, so the cancellation was correct. Two jobs hit the same
dispersive-cost-underestimate this session; one was caught in flight and one
was not.

---

## Fourteenth snapshot — 2026-08-01 ~13:35 (superseded, kept for its reasoning)

It adds the four jobs
submitted after the thirteenth was written plus the two diagnostics it omitted;
the thirteenth's *reasoning* about the dispersion generational break still
stands and is not repeated.

**Ten BEEF jobs are live.** Nothing has landed since 11608940. Grouped by what
they answer:

| job | state at snapshot | answers |
|---|---|---|
| **11608875** `[0-3]` + **11608876** | `_0` ✅ 55:15, `_1` ✅ 41:22, `_2` ✅ 49:33; `_3` 🟢 running; merge ⏸ `afterok` | **gap 4 / check 1 / P1** re-measured dispersive, fully cache-cold |
| **11608690** | ✅ **COMPLETED 1:28:14** | **P2 CLOSED** — first calibrated hybrid to survive since the `beta_from_phases` fix. No guard fire |
| **11608692** | ✅ COMPLETED 00:07, but its **number is contaminated** | gap 5 calibrated comparison — see the race below |
| **11608703** | 🟢 running 1:08:27 | coupler L₅₀ re-verification through the shared fixed code path |
| **11608941** | ✅ **COMPLETED 1:03:21** @ 14:00 | **gap 2** re-audited dispersive — **byte-identical**, as predicted. Folded in |
| **11609211_1** + ~~`_2`~~ → **11609680** + **11609681** | `_1` 🟢 running; `_2` ⛔ cancelled 46 min in, resubmitted as 11609680 (🟢, `-t 720`); comparison resubmitted as 11609681 | **the P4 test** — 2×2 gap 5 dispersive. **Rescued from a certain timeout** — see below |
| **11609215** `[1-2]` + ~~11609216~~ → **11609687** | `_1` ✅ **13:44**, `_2` 🟢 running; comparison resubmitted `afterok` | 1×2 gap 5 dispersive. **`_1`'s FDTD reference already folded in** — see below |
| **11608876** → **11609686** | ⏸ resubmitted `afterok:11608875` | gap 4 merge, moved to the funded account |
| **11608940** | ✅ COMPLETED 07:13 @ 13:04 | **gap 3** re-audited dispersive — folded in |

**Gaps 2 and 3 are both re-measured dispersive and both came back
byte-identical** (11608940 @ 13:04, 11608941 @ 14:00). Two independent
confirmations, so the rule is now evidence rather than analogy: **a quantity
evaluated at a scalar wavelength is invariant across the dispersion break;
only genuinely broadband quantities move** — gaps 4 and 5, and nothing else.
Spending an hour of compute on a result expected to be unchanged was the point.

**⚠ Do not reach for the ε grep on every job.** It works only where meep builds
a broadband medium and prints `dielectric constant epsilon diagonal`.
11608941 runs through meow/MPB and prints no such line at all — grepping
returns nothing, which is *not* evidence of frozen material. There its log
header carries the discriminator instead: `wl=1.31`, a single scalar
wavelength. Match the check to the code path.

**11608940's result, and the correction it forced.** Gap 3 came back
**byte-identical** to its frozen-era numbers. I had recorded beforehand that
gaps 2 and 3 were stale and that gap 2 was "most likely to move" — that was
**wrong**, and the job disproved it rather than confirming it. A dispersive
medium only differs from a frozen one *across a band*; `section_mode_neff()`
builds its medium at a scalar wavelength deliberately, and meow's
`Environment(wl=…)` is exact Sellmeier either way. Only the genuinely broadband
gaps (4 and 5) move. 11608941 is left running to measure that for gap 2 rather
than infer it by analogy. Full account in `resolved/mmi_1x2_validation.md`.

**Both gap-5 pairs preserved their frozen baseline before launching** —
`comparison_FROZEN_pre20260801.npz` in each of
`example_output/verification/mmi_{1x2,2x2}_hybrid_vs_fdtd/`. The comparison
script overwrites `comparison.npz`; without the copy the pre-flip half of the
P4 test would have been unrecoverable.

**Watch item on the 2×2 hybrid.** It runs at `N_MODES = 14` with zero margin
over the guided count — a duplicate-mode guard trip there fails the job rather
than degrading it.

**⚠ 11609211_2 was cancelled 46 min in because it was going to time out, and
the way that was caught generalises.** `run_mmi_2x2_examples.batch` carries one
`-t 300` shared by all three of its examples, sized when the slowest fitted.
The frozen-era hybrid (11585844_2) took **04:11:55** against that 300 — **48
minutes of headroom**. Dispersion costs ~1.5× (measured: 1×2 FDTD 08:57 → 13:44,
jobs 11585845_1 → 11609215_1), putting the dispersive hybrid near **06:25**.
It was never going to finish. Resubmitted as **11609680** with `-t 720`.

- **The limit cannot be raised in flight.** `scontrol update TimeLimit` returns
  `Access/permission denied` for a non-admin, so the only remedy is cancel and
  resubmit. That makes the check worth doing *at submit time*, not after.
- **The check itself is cheap and should be routine:** before re-running
  anything under the new material model, look up the frozen-era elapsed with
  `sacct -X -j <old> --format=Elapsed,Timelimit` and multiply by ~1.5. A shared
  array time limit is the usual trap, because it is sized for whichever task
  happened to be slowest at the time.
- **The 46 min are not lost** — `example_sim_cache` persists per solve, so
  completed section solves are reused on the rerun.

**⚠ The second hazard was worse, and it was invisible.** The gap-5 comparison
jobs depended on their producers with **`afterany`**, which a TIMEOUT satisfies
just as well as success. Had the hybrid timed out, the comparison would have
fired anyway and compared the **new dispersive FDTD** against the **old frozen
hybrid pickle still on disk** — yielding a plausible, wrong, unlabelled gap 5
number. That is precisely what corrupted **11608692** earlier today, caught only
because `provenance.json` happened to record an mtime.

**Both comparisons are now `afterok`** — 11609681 (2×2) and 11609687 (1×2).
Use `afterok` for any job that consumes a producer's output artifact;
`afterany` is only correct when the dependent job's job is to *report on*
whatever happened, not to read a result.

**Account.** New submissions go to **`gts-sr108-ece`**, which has the balance;
`gts-sr108-gedc_main` is what the older batch scripts hardcode. Applies to
newly run jobs only — jobs already queued or running are left alone.
`scontrol update JobId=<id> Account=...` does not work here regardless: it
fails with a misleading `Requested time limit is invalid`.

**11609215_1's number, folded in without waiting for its merge.** The
full-device FDTD reference alone is enough to measure how far the *reference*
moved across the dispersion break, which is worth having separately from gap 5
itself: max |ΔT| **0.00651** (1.65% of mean T), `out_0` and `out_1` agreeing to
five decimals. Read against the P4 hypothesis it is **mixed** — the scale is
2× the 2×2's entire ΔT swing, which supports it, but the curve crosses zero
**twice**, which a single frozen-index term cannot do. Both halves recorded, in
`resolved/mmi_1x2_validation.md` gap 5 (table) and `resolved/mmi_2x2_validation.md` gap 5
(P4 box), **before 11609212 lands**.

**⚠ A pickle-path race corrupted 11608692, and it will recur.** The comparison
scripts load their FDTD reference **by path**, so they read whatever last wrote
there. 11609215_1 overwrote `example_output/mmi_1x2/mmi_1x2_fdtd/smatrix.pkl`
with the *dispersive* reference at 17:34:48 UTC; 11608692 ran at 17:43 and
compared it against 11608690's *frozen* calibrated hybrid. Its printed
`max|ΔT| = 0.0313` mixes two material eras and must not be quoted. The valid
frozen-era number (**0.02876 / 7.30%**) was recovered offline from the
preserved baseline at zero compute cost — see `resolved/mmi_1x2_validation.md` gap 5.

Three lessons, all cheap:

1. **Two jobs writing and reading one pickle path need an ordering constraint,
   and `afterany:` on the *producer* is not it** — 11608692 depended on
   11608690 (its hybrid), not on the unrelated job that clobbered its FDTD
   input. Nothing in the dependency graph could have caught this.
2. **It was caught only because `provenance.json` records `fdtd_pkl_mtime`.**
   That is P6's argument arriving from a second direction: extend provenance
   from goldens to comparison outputs, and record the **resolved material
   config**, not merely an mtime — an mtime only helps if someone already
   suspects a problem.
3. **The meep log settles material era in one grep**, independent of any
   provenance — see below.

**How to check a run actually got the dispersive material, in one line.** Worth
recording because "the default flipped" and "this job ran dispersive" are
different claims, and the second is the one that matters when folding a number
in. Every meep log prints `dielectric constant epsilon diagonal`. For `SOI_160`
silicon that is **12.4744** dispersive (ε∞ of the fit, plus one Lorentzian) and
**12.2774** frozen (exactly Sellmeier n² at band centre, 3.5039² ). The two are
1.6% apart and unmistakable:

```
python -c "import numpy as np; from beef.materials.stack import SOI_160 as s; \
  wl=np.linspace(1.26,1.36,11); \
  print([s.si_medium(wl,dispersive=d).epsilon_diag.x for d in (True,False)])"
```

**All three completed 11608875 tasks print 12.4744**, so gap 4 is genuinely
being re-measured dispersive. This does not depend on trusting the cache
fingerprint or the wall clock.

**Wall clock, by contrast, is not usable here.** Frozen-era per-length was a
clean monotone 18:30 / 25:33 / 27:01 for L = 0.5/1.0/1.5. Dispersive came back
**55:15 / 41:22 / 49:33** — surcharges of 2.99× / 1.62× / 1.83×, and
*inverted* in L, with the shortest section slowest. Compute cost cannot depend
on L that way, so these numbers are contention-dominated (1 rank inside a
2-node allocation, four tasks placed wherever the scheduler had room). Do not
build a cost model on them; the 1.53× from 11609215_1 below is the cleaner
figure and even that is one sample.

**Cost note.** 11609215_1 took **13:44** against the frozen-era 11585845_1's
**8:57** for the same FDTD example — a **1.53×** dispersion surcharge, from the
Lorentzian susceptibility meep now has to time-step. Real but modest. The
hybrid half is the expensive one either way (11585845_2 was **1:46:59**), so
size the eventual golden regeneration (P7) off that, not off this.

---

**★ THIRTEENTH SNAPSHOT — 2026-08-01, after material dispersion was turned on
by default. Superseded by the fourteenth above for job states; still current
for the generational-break reasoning.**

**Everything measured before this point is on the old material model.** The
`dispersive` default flipped to True (`material_model_validation.md` step 3),
which adds `resolved_dispersive` + the fitted `resolved_{si,sio2}_reparam` to
the cache fingerprint — so **every cache entry misses and every golden is stale
again.** There is a clean generational break here; do not mix numbers across it
without saying which side they came from.

| job | state | note |
|---|---|---|
| **11608448** `[0-3]` | ✅ all COMPLETED — 18:30 / 25:33 / 27:01 / **52:25** (≈2 h 04 total) | **frozen-material era.** Its cache entries are still on disk but unreachable under the new key |
| ~~11608450~~ | ⛔ **CANCELLED by me, deliberately** | the gap-4 merge runs `cache_only=True`; under the new key it would have raised `CacheMissError` through no fault of its own |
| **11608875** `[0-3]` + **11608876** | 🟢 re-submitted | gap 4 re-measured under the configuration that actually ships. Fully cache-cold, so expect ~2 h again |
| **11608690** | 🟢 running | 1×2 calibrated hybrid — **frozen-material era**, and that is fine: it exists to test the `beta_from_phases` fix, which is orthogonal |
| **11608692** | 🟢 released | gap 5 calibrated comparison. Safe across the break — it only loads two saved `SMatrix` pickles, never touching the cache. Verified before releasing rather than assumed |
| **11608703** | 🟢 running | coupler L₅₀ re-verification, also frozen-era, also orthogonal (shared `beta_from_phases` path) |

**Sequencing note worth keeping.** The two pending merges were put on
`scontrol hold` *before* the library edit, not after — a `cache_only=True` job
firing mid-edit would have failed spuriously and looked like a real cache miss.
One was then released (pickles only) and one cancelled (cache-dependent). That
distinction is worth making explicitly whenever a fingerprint changes under a
running queue.

---

**★ TWELFTH SNAPSHOT — 2026-08-01, earlier the same day. Superseded by the
thirteenth above for job states; still current for the parity-work rationale.**

Five jobs submitted, all for the **1×2 MMI parity work** (`README.md`'s new
parity backlog, items P1 and P2). One inert job cancelled (`11597199`, the
orphaned mux merge — done).

| job | name | for | state at snapshot |
|---|---|---|---|
| **11608448** `[0-3]` | `beef_mmi_1x2_engine_agreement` | **P1 / gap 4 / check 1** — this device's last unrun check | `_0` ✅ 18:30, `_1` ✅ 25:33, `_2` ✅ 27:01, `_3` 🟢 running |
| **11608450** | …`_merge` | assembles the ratio from cache | ⏸ `afterok:11608448` |
| **11608690** | `beef_mmi_1x2_hybrid_calibrated` | **P2** — successor to 11589770, first calibrated run since the `beta_from_phases` fix | 🟢 running |
| **11608692** | `beef_mmi_1x2_gap5_calibrated` | gap 5's calibrated comparison, successor to 11589779 | ⏸ `afterany:11608690` |
| **11608703** | `beef_dc_hybrid_item14` | **coupler re-verification** — the shared `beta_from_phases` fix touches the coupler's path too | 🟢 running |

**P1 — gap 4 cost, and the reshape paid off differently than on the 2×2.**
Per-length wall clock 18:30 / 25:33 / 27:01 for L = 0.5/1.0/1.5 — about **2 h
of serial work total** against the 2×2's 13.4 h, as the 6.2×-shorter body
(`mm_length` 4.5 vs 28) predicted. A single serial job would have fitted in one
slot here; the array shape was adopted anyway because it is free and a
per-length timeout then costs one length instead of four. Launched with bare
`python` (1 rank) deliberately — the script has no rank guards, so `srun` would
race 48 ranks on one `np.savez`.

**P2 — why these two reruns exist.** The 2026-07-31 diagnosis of 11589770's
failure (a π sign flip in `T_ii`) was **wrong**; the real cause was a 2π
turn-count bug in `beta_from_phases`, now fixed. See that job's row below for
the full account. 11608690 is the first calibrated run on the fixed code: if it
completes without the residual guard firing, `calibrate_body_neff` works on
this device for the first time and 11608692 can finally measure gap 5
calibrated.

**11608703 is the backward-compatibility check, and it is deliberately a real
run rather than an argument.** The fix is in shared library code that the
coupler's own closed item 6 depends on. Offline evidence is strong — at the
coupler's production seed the fix changes **0 of 66 cached cells**, and 28 of
29 existing tests were unchanged — but that used a *proxy* seed rather than the
run's own MPB output. This job re-derives the coupler's calibrated 50:50 length
through the real path; it should reproduce ~15.15 µm against the 15.13 µm
full-device reference. Chose `run_dc_hybrid_item14.batch` over
`run_golden_directional_coupler_hybrid.batch` on purpose: the golden generator
would **overwrite** `tests/golden/directional_coupler_hybrid.json`, destroying
the pre-fix reference this comparison needs (and goldens are deferred work
besides). The archived stub-sweep script would have been the closest match to
the original 0.11% measurement, but its `make` target no longer exists.

---

**★ ELEVENTH SNAPSHOT — 2026-08-01, verified against `squeue -u jhiesener3`
and `sacct -u jhiesener3 -X`. Superseded by the twelfth above for job states;
still current for the 2×2 gap-4 closure and the mux diagnosis.**

**The queue is now empty of BEEF work.** Nothing is running. One inert entry
remains: **11597199** `beef_modal_mux_examples_merge`, `PENDING` with
`Reason=DependencyNeverSatisfied` — its parent array timed out, so it can
never start. It holds no resources; cancel it as hygiene (`scancel 11597199`).
Not done automatically. (`11585814_7` `elevator_resolution` COMPLETED 07:28
after 22:47:08 — **unrelated project, ignore**.)

Both reshapes submitted in the 21:20 window have now landed, with **opposite
outcomes**:

| job | name | outcome |
|---|---|---|
| **11596845** `[0-3]` | `beef_mmi_2x2_engine_agreement` | ✅ **ALL FOUR COMPLETED.** `_0`/`_1` in **36 s each** (predicted cache hit on 11591249's L=2/L=4); `_2` (L=6) 3:35:10; `_3` (L=8) **5:07:14** |
| **11596851** | …`_merge` | ✅ **COMPLETED** 03:07, 25:29 — assembled the answer from cache, `cache_only=True` raised nothing |
| **11597177** `[0-2]` | `beef_modal_mux_examples` | ❌ **ALL THREE TIMEOUT** (02:44/02:52/03:01), full 5 h × 4 nodes each — **third consecutive timeout**, cause now understood |
| **11597199** | …`_merge` | ⛔ **orphaned**, `DependencyNeverSatisfied`, cancel it |

**✅ The 2×2 MMI's gap 4 is CLOSED — its last open check.** Mean ratio
Δβ_FDTD/Δβ_EME = **0.9953**: the engines agree on body phase accumulation to
**0.47%**, versus the coupler's ~10% on the identical measurement, so
`calibrate_body_neff` is **not** indicated on this device. T diagonal to
1.21e-05, worst phase-fit residual 4.02e-06 rad (25,000× inside the 0.1 rad
guard that broke the 1×2's calibration). The job log prints both as `0.0000`
at `:.4f` — they are genuinely that small; read off
`example_output/verification/mmi_2x2_engine_agreement/engine_agreement.npz`,
**not** defaulted or skipped, which was checked before recording this.
Folded into `resolved/mmi_2x2_validation.md` gap 4 (boxed) and `README.md`'s table.
Caveat carried there: measured at res 40, which gap 2 says is not converged in
absolute n_eff; both engines read at the same mesh so the ratio partially
cancels it, but 0.47% is not certified mesh-independent.

**The reshape's cost model was right.** Per-length wall clock 96 / 183 / 215 /
**307** min for L = 2/4/6/8 — **13.4 h of serial work**, so the original 5 h
single-slot job could never have finished, and L=8 alone would have blown it.
As an array it closed in 5:07 + 25 min. Worth copying for any other
multi-length `solve_section()` sweep.

**❌ The mux reshape did the opposite, and the reason is a mistake worth not
repeating: it silently dropped MPI.** The new batch script calls `python
modal_multiplexer/python/modal_multiplexer_hybrid.py` directly to avoid
`make` running all 3 ports in one process — but the Makefile target is
`$(runner) $(my_python) …` and on Phoenix `runner` = **`srun --mpi=pmi2`**
(`BEEF/examples/Makefile:5`). Bypassing `make` threw away the launcher with
the target:

| job | launcher | meep ranks | s/step |
|---|---|---|---|
| 11585847_0 / 11590663_0 (via `make`) | `srun --mpi=pmi2` | **96** | 0.0064 |
| 11597177_0/_1/_2 (reshaped, direct) | none | **1** | 0.27–0.36 |

`Using MPI version 4.1, 1 processes` in all three new logs against
`96 processes` in both old ones — **~42× slower per timestep**, while holding
4 nodes × 24 ranks idle (`cput=20-00:35:12` for `walltime=05:00:22` on 96
cpus). **Second, independent flaw:** the reshape assumed sections 0-3 were
already warm on disk; they were not — each task re-ran the same three
`Initializing structure` + MPB + time-stepping blocks at near-identical log
offsets. Whether the earlier timeouts never persisted them, or the cache key
differs between 96-rank and 1-rank runs, is **not established**.

**The mux is deferred, so no fourth attempt was made** (per instruction, and
consistent with `deferred/modal_multiplexer_validation.md`). Fix recipe
recorded in that doc's gap 6 for whenever it resumes: restore
`srun --mpi=pmi2`, *then* verify the cache actually hits, then resubmit — and
do not merely raise the wall clock, since at 1 rank no affordable limit
finishes. **Gap 6 has now consumed ~25 h of compute across three attempts and
is still unanswered either way.** The array-plus-merge *pattern* is not at
fault — the same design worked first try on the 2×2.

**Still open from earlier, untouched by this snapshot:** `beef_dc_examples`
`_0` (11592462_0) TIMEOUT, cause still undiagnosed — the standing
recommendation is a short-wall-clock resubmit watched live, **not** submitted
here. The mmi_1x2 π-sign-flip diagnosis and the ring y-normal excite blocker
are both analysis work with no job in the queue.

---

**Snapshot taken 2026-07-31, 16:35 EDT** (ninth snapshot of the day: 07-30
21:46 → 07-31 08:22 re-queue → 08:46 → 10:35 → 10:53 → 12:05 → 13:05 → 13:44 →
now). This update: mmi_1x2 gap 2 (11591224) completed and closed (res 120,
`RESOLUTION` held at 40 by decision — see `resolved/mmi_1x2_validation.md`); the
coupler's cold regeneration pass launched, 11592462 (example outputs,
array 0-2) + 11592463 (golden) — see the flagged anomaly on array task 2
below before treating that output as a clean fresh result. This is a point-in-time list, not a
live view — job states, elapsed times and queue positions all move. Check
`squeue -u jhiesener3` (or `sacct -j <id>` once finished) for current status;
this file exists so the *purpose* of each job is recorded somewhere durable,
not to be a dashboard. Update or delete entries as jobs complete and their
results get folded into the relevant validation doc.

**Reshaped and resubmitted, 2026-07-31 ~21:20 EDT.** Both `beef_modal_mux_
examples` and `beef_mmi_2x2_engine_agreement` timed out with an understood
cost problem (see their entries below), so both scripts were reshaped into
an array-job-plus-merge-job pair instead of one long serial process, and
resubmitted:

| new job | replaces | shape |
|---|---|---|
| **11596845** `[0-3]` | 11591249 | `verify_mmi_2x2_engine_agreement.py`, now one array task per length (2/4/6/8 µm), 8h/task |
| **11596851** | (new) | merge step for 11596845 — `--dependency=afterok:11596845`, assembles the ratio/fit from cache, 30 min |
| ~~11596846~~ → **11597177** `[0-2]` | 11590663 | `modal_multiplexer_hybrid.py`, one array task per output port (top/mid/bot), 5h/task |
| ~~11596852~~ → **11597199** | (new) | merge step — `--dependency=afterok:11597177`, assembles the SMatrix from cache, 30 min |

See `verify_mmi_2x2_engine_agreement.py`'s and `modal_multiplexer_hybrid.py`'s
own "Reshaped 2026-07-31" docstring sections, and the four
`batch/.../*.batch` files, for the mechanics (`SLURM_ARRAY_TASK_ID` selects
one length/port; `cache_only=True` in the merge step turns a missing result
into a loud `CacheMissError` instead of a silent multi-hour re-simulation).
**`beef_dc_examples` idx 0 (11592462_0) also timed out** during this same
window (21:34) but its cause is not understood yet — **not** reshaped or
resubmitted; see its own entry below before touching it.

**Correction, 21:34-21:41 EDT:** 11596846 `[0-2]` (the modal-mux array,
first attempt) **FAILED instantly** (~8s each task) — the new batch script
had the wrong relative path to `modal_multiplexer_hybrid.py`
(`python/modal_multiplexer_hybrid.py` instead of
`modal_multiplexer/python/modal_multiplexer_hybrid.py`, a mistake made while
writing this reshape, not a cluster/library issue). Its merge job
(11596852) was left stuck on `DependencyNeverSatisfied` and was cancelled.
Fixed and resubmitted as **11597177** (array) / **11597199** (merge,
`--dependency=afterok:11597177`). **11596845** (mmi_2x2 engine agreement
array) had no such bug: tasks `_0`/`_1` (L=2, L=4) both **completed in 36s**
each, confirming the predicted cache hit on 11591249's already-computed
lengths; `_2`/`_3` (L=6, L=8 — the real new work) are running for real as of
this snapshot.

**Re-verified against `squeue`/`sacct` 2026-07-31, 21:12 EDT (tenth snapshot).**
Two jobs that were "🟢 RUNNING" in the 16:35 snapshot both **hit their 5h wall
and were killed by SLURM, not by an error**: `beef_modal_mux_examples`
(11590663_0) at 18:58, and `beef_mmi_2x2_engine_agreement` (11591249) at
19:46. Neither is a crash — both show `ExitCode=0:0`, `State=TIMEOUT` — but
neither produced its intended result: the 2×2 MMI's gap 4 (engine agreement,
its last open check) did **not** close, and the mux's gap 6 exercise run
(`warn_on_transfer_duplicate` against real interfaces) is **still** unproven
either way. `beef_golden_dc_hybrid` (11592463) completed cleanly (1:09:17).
`beef_dc_examples` array task 0 (11592462_0, hybrid) is still running at
4:38:25 of its own 5:00:00 limit as of this snapshot — **~21 min of budget
left; watch this one, it may time out the same way** unless it finishes
first. See the three detail sections below (11590663, 11591249, 11592462) for
what each actually reached before stopping.

**Re-queue, 2026-07-31 08:22 EDT:** the five jobs below never started —
`scontrol show job` confirmed `Account=gts-sr108-ece`,
`Reason=AssocGrpBillingMinutes` on the three parents (that association is out
of billing minutes; `gts-sr108-gedc_main` had headroom). Cancelled all five
(old IDs 11577985/11577986/11578018/11578019/11578026) and resubmitted under
`gts-sr108-gedc_main`. Also switched all other batch scripts hardcoding
`gts-sr108-ece` to `gedc_main` so this doesn't recur.

| job | name | state (at snapshot) | for |
|---|---|---|---|
| 11575516 | `beef_mode_spec_converge` | ✅ **COMPLETED** 21:17, 3:36:13 — results folded in, entry retained for provenance | [resolved/fdtd_vs_EME_mode_validation.md](resolved/fdtd_vs_EME_mode_validation.md) item 1 + [modal_multiplexer_validation.md](modal_multiplexer_validation.md) gap 2 |
| ~~11577985~~ → **11585844** `[1-2]` | `beef_mmi_2x2_examples` | ✅ **COMPLETED** — `_1` 1:16:45, `_2` 4:11:55; both SMatrix files on disk | [resolved/mmi_2x2_validation.md](resolved/mmi_2x2_validation.md) gap 5 (data) |
| ~~11577986~~ → **11585849** | `beef_mmi_2x2_gap5` | ✅ **COMPLETED** 12:46, 0:18 — **worst max\|ΔT\| 0.0043 (0.82%), sign FLIPS across the band (not the 1×2's uniform signature); results folded in** | [resolved/mmi_2x2_validation.md](resolved/mmi_2x2_validation.md) gap 5 (comparison) |
| ~~11578018~~ → **11585845** `[1-2]` | `beef_mmi_1x2_examples` | ✅ **COMPLETED** — `_1` 8:57, `_2` 1:46:59; both SMatrix files on disk | [resolved/mmi_1x2_validation.md](resolved/mmi_1x2_validation.md) gap 5 (data) |
| ~~11578019~~ → **11585850** | `beef_mmi_1x2_gap5` | ✅ **COMPLETED** 10:28, 0:09 — **max\|ΔT\| 0.0359 (9.11%), results folded in** | [resolved/mmi_1x2_validation.md](resolved/mmi_1x2_validation.md) gap 5 (comparison) |
| ~~11578026~~ → ~~11585847~~ → **11590663** `[0]` | `beef_modal_mux_examples` | 11585847 ❌ **TIMEOUT** 13:42 (full 5h, 4 nodes) — sections 0-3 cached to disk before it died mid-section-4; resubmitted unchanged as 11590663 → ❌ **TIMEOUT again**, 18:58, full 5h — got further (deep into a per-port section-4 eigenmode-overlap monitor, band 13/14) but still did not finish or reach the item-8 `warn_on_transfer_duplicate` scan this run exists to exercise → reshaped as **11597177** `[0-2]` → ❌ **THIRD TIMEOUT 2026-08-01**, all 3 tasks, full 5h × 4 nodes each. **Cause found: the reshaped batch script dropped `srun --mpi=pmi2`, so meep ran on 1 rank instead of 96 (~42× slower per step) with 96 cores idle**; and the "sections 0-3 are cached" premise it rested on is false — every task recomputed them. **NOT resubmitted, device is deferred**; merge 11597199 orphaned on `DependencyNeverSatisfied`, cancel it. ~25 h of compute spent, gap 6 still unanswered | [modal_multiplexer_validation.md](modal_multiplexer_validation.md) gap 6 (exercise run) |
| ~~**11586046**~~ → **11588309** | `beef_diagnose_items_11_12` | 11586046 ❌ **CANCELLED** 11:12 (hung, not slow); refixed → **11588309 ✅ COMPLETED 11:50, 34:55 — BOTH ITEMS RESOLVED:** item 11 closed (all four grids agree at res 40 too), item 12 decisive (laddering recovers a 14th guided mode at equal cost) | [resolved/fdtd_vs_EME_mode_validation.md](resolved/fdtd_vs_EME_mode_validation.md) items 11 + 12; it also opened item 14 |
| **11586098** | `beef_diagnose_item_13` | ✅ **COMPLETED** 10:50, 1:54:50 — **C vs A inert (max \|Δn_eff\| 0.000e+00), cost 3.2×; results folded in** | [resolved/fdtd_vs_EME_mode_validation.md](resolved/fdtd_vs_EME_mode_validation.md) item 13 |
| ~~11586986~~ → ~~11588997~~ | `beef_dispersive_convergence` | 11586986 ❌ **CANCELLED** 12:00 — **it was measuring nothing**, `dispersive=` never reached the mode solve (library bug, fixed); refixed → **11588997 ✅ COMPLETED** 13:56, 1:56:53 — **pin holds at centre, holds at the edge only up to the 1e-3 tolerance production actually uses (needs res 80, not 60, under dispersion — which is what's already shipped); does not fully converge at 3e-5 at the edge within res 20–120** | [material_model_validation.md](resolved/material_model_validation.md) work item 4 (cheap half) |
| **11588569** | `beef_audit_mmi_2x2_gap6` | ✅ **COMPLETED** 11:58, 20:18 — **gap 6 closed:** identity `[0…13]`, worst \|Δn_eff\| 0.0007, 14/14 guided, **neither guard fired** | [resolved/mmi_2x2_validation.md](resolved/mmi_2x2_validation.md) gap 6 |
| **11589770** → **11608690** | `beef_mmi_1x2_hybrid_calibrated` | ❌ **FAILED** 16:59, 4:17:26 — `section_mode_beta`'s phase-fit residual (2.101 rad) blows its own 0.1 rad sanity guard on the `[-2.0, 2.0]` section, 3 lengths. ~~DIAGNOSED 2026-07-31 as a **π sign flip in `T_ii`** at 4 of 88 cells; narrow upstream mode-extraction bug; the guard caught a real defect~~ → **that diagnosis was WRONG. ✅ REAL CAUSE FOUND AND FIXED 2026-08-01: a 2π turn-count bug in `beta_from_phases` itself.** It rounded wrap counts independently per length against one seed slope, so the ~2% MPB-vs-FDTD seed error — *the quantity `calibrate_body_neff` exists to measure* — drifted with L until a length straddled a `.5` boundary. The two readings are algebraically indistinguishable on the residual (both give `δ·(+⅓,−⅔,+⅓)`); physics separates them: turn-count → n_eff **2.783, +1.96%** over MPB, in family with the coupler's +2.233%/+2.361%; sign-flip → **2.101, −23%**, wrong sign and size. Seed-free proof: second difference `(φ₃−2φ₂+φ₁)/2π` is an exact integer to **0.00308 in all 88 cells** — nothing needs a sign flip. Not a corner case either: 4025-trial sweep, old code wrong 13.7%, new 0%. The doc's own "unresolved" branch margin 0.922 **was** the tell → now 0.0785. **The guard did catch a real defect — its own unwrapping.** Re-run in flight as **11608690** | [resolved/mmi_1x2_validation.md](resolved/mmi_1x2_validation.md) gap 5 (boxed, rewritten) |
| **11589779** | `beef_mmi_1x2_gap5_calibrated` | ❌ **FAILED** 16:59, 0:00:09 — dependent on 11589770, nothing to load, expected given the above | [resolved/mmi_1x2_validation.md](resolved/mmi_1x2_validation.md) gap 5 follow-up (comparison) |
| **11590785** | `beef_mmi_2x2_gap2_dbeta` | ✅ **COMPLETED** 1:05:48 — **gap 2 closed at res 100**, Δβ_01 converges to 0.05% by res 120 | [resolved/mmi_2x2_validation.md](resolved/mmi_2x2_validation.md) gap 2 |
| **11591224** | `beef_mmi_1x2_gap2_dbeta` | ✅ **COMPLETED** 15:45, 1:05:15 — **gap 2 closed at res 120** (non-monotone step, unlike the 2×2's clean decay); `RESOLUTION` held at 40 by decision, gap recorded as a known mesh caveat instead | [resolved/mmi_1x2_validation.md](resolved/mmi_1x2_validation.md) gap 2 |
| **11591249** → **11596845** `[0-3]` + **11596851** | `beef_mmi_2x2_engine_agreement` | ✅ **CLOSED 2026-08-01 by the reshape** — array all 4 COMPLETED (`_0`/`_1` 36 s cache hits, `_2` 3:35:10, `_3` 5:07:14) + merge 25:29. **Ratio 0.9953 — engines agree to 0.47%, `calibrate_body_neff` NOT indicated on this device.** T diagonal to 1.21e-05, fit residual 4.02e-06 rad (log shows both as `0.0000` at `:.4f`; real values read from the `.npz`). Gap 4 was this device's last open check. Original attempt below: ❌ **TIMEOUT** 19:46, full 5h — 4 of the 4 lengths' first run completed (2/4/6/8 µm all show a `run 0 finished` line), a 5th run (second pass, presumably for the reverse-direction or FDTD-vs-EME comparison half) was ~5,300 of a likely ~30,000 timestep run when killed — well under half done on that leg. Gap 4 (this device's last open check) did **not** close | [resolved/mmi_2x2_validation.md](resolved/mmi_2x2_validation.md) gap 4 / [README.md](README.md) check 1 — never measured for either MMI; 4 real FDTD solve_section() runs, adapted from the coupler's item 6 script |
| **11592462** `[0-2]` | `beef_dc_examples` | `_0` (hybrid) ❌ **TIMEOUT** 21:34, full 5h — see detail section below, not yet diagnosed; `_1` (mode) ✅ **COMPLETED 16:37, 3:49** — legit, not flagged: serial MPB even/odd supermode sweep (1 process), inherently cheap, no BEEF solver/cache involved; L_c = 32.813 µm at λ=1.55, gap sweep 0.15–0.5 µm recorded, `_2` (fdtd) ✅ **COMPLETED 16:33, 0:00:17** — **flagged, not yet trusted as a clean result:** 17s for a full 3-D broadband FDTD sweep is a cache hit, not a fresh compute. Confirmed NOT the item-8 fingerprint bug (`FDTDSolver.solve()` does splat `resolved_material_config()` into its cache key, `fdtd.py:707`, inside `solve()` itself) — so whatever it hit was already computed post-item-8, most likely the ground-truth length_sweep's own 15.13 µm grid point. Printed numbers (T=0.5061/0.4933, IL≈-3 dB) look physically sane, but this has not been cross-checked against that sweep's stored value yet — do that before citing this run's output as the "outstanding — not yet run" regeneration `resolved/directional_coupler_validation.md`'s Verification section calls for | `resolved/directional_coupler_validation.md` Verification → Outstanding (example-output regeneration) |
| **11592463** | `beef_golden_dc_hybrid` | ✅ **COMPLETED** 17:42, 1:09:17 | `resolved/directional_coupler_validation.md` Verification → Standing (goldens stale since item 8) |
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
`resolved/fdtd_vs_EME_mode_validation.md`. The 1×2 MMI's whole gap-5 chain completed end
to end earlier and its numbers are in `resolved/mmi_1x2_validation.md`. `11585847` has
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
`resolved/fdtd_vs_EME_mode_validation.md`: item 12 gains the measured band table showing
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
  `resolved/mmi_2x2_validation.md` gap 5. Output under
  `example_output/verification/mmi_2x2_hybrid_vs_fdtd/`.

### 11585845 `[1-2]` — `beef_mmi_1x2_examples` (array indices 1,2 only)
- Same shape as 11585844, for the 1×2 MMI: `make mmi_1x2_fdtd
  mmi_1x2_hybrid`.
- **Why:** `resolved/mmi_1x2_validation.md` gap 5, same reasoning as the 2×2's. This
  device's basis mismatch (gap 3, `N_MODES` 12→8) is already fixed, so —
  unlike the 2×2 — there's no open mode-identity concern underneath this
  comparison; only gap 2 (resolution) remains unsettled.

### 11585850 — `beef_mmi_1x2_gap5` — ✅ COMPLETED
- **Script:** `batch/verification/run_verify_mmi_1x2_hybrid_vs_fdtd.batch` →
  `examples/verification/python/verify_mmi_1x2_hybrid_vs_fdtd.py`
- Same shape as 11585849, ported to `mmi_1x2`'s port names (`out_0`/`out_1`
  vs `out_top`/`out_bot`). `--dependency=afterany:11585845`.
- **Why:** `resolved/mmi_1x2_validation.md` gap 5.
- **Result (10:28, 9 s):** max |ΔT| **0.0359** at λ = 1.33 µm, mean |ΔT|
  0.0253, **9.11%** of mean T_fdtd. `check_passive()` true both sides; the two
  output ports agree to 4 decimals, so splitter symmetry survives the cascade.
  **Every ΔT is negative** — a systematic under-prediction, not scatter, which
  is the Δβ-discretization signature from
  `resolved/directional_coupler_validation.md` items 6/14 (coupler: ~5%
  uncalibrated → 0.11% with `calibrate_body_neff=True`, never tried here).
  Folded into `resolved/mmi_1x2_validation.md` gap 5; output under
  `example_output/verification/mmi_1x2_hybrid_vs_fdtd/`. Provenance records
  beef `1cc8edc`.

### 11585847 `[0]` → 11590663 — `beef_modal_mux_examples`
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
  `warn_on_transfer_duplicate()` (`resolved/fdtd_vs_EME_mode_validation.md` item 8)
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
- **❌ TIMEOUT 2026-07-31 13:42, full 5h.** `sacct` confirms `State=TIMEOUT`,
  `ExitCode=0:0` — killed by SLURM mid-run, not a crash. The batch script's `-t
  300` was sized off the earlier device estimates, not this device's own cost:
  this is the repo's most expensive hybrid chain (16 modes × 4 interfaces) and
  its own `save_dir`/`cache=` mean a re-run doesn't restart from zero, so the
  fix is a plain resubmit rather than a reshape (contrast the 11586046 case,
  where the batch script's *cost model* was wrong and needed reshaping, not
  just more wall clock).
  - **What finished, read off `example_sim_cache/` mtimes against the
    script's fixed section order (0 input FDTD, 1 EME MMI1 body, 2 FDTD
    transition, 3 EME MMI2 body, 4 output FDTD, once per port):** section 0 by
    08:50, section 1 (EME) by 10:10, section 2 (FDTD transition,
    `CompoundGeometry`) by 10:23, then a **2h17m** gap to 12:45 for section 3's
    EME solve (a `solve_section_length_sweep` call — by far the single most
    expensive piece of this device), then one more EME
    `CompoundGeometry`-keyed cache write at 13:40:31 — 1m41s before the kill.
    The tail of `reports/11585847_0.out` is mid-way through building the first
    section-4 (output) FDTD grid when SLURM cancelled it. **Net: sections 0-3
    are fully cached; zero of the 3 per-port section-4 solves completed.**
  - **Resubmitted unchanged as 11590663** (same `batch/run_modal_multiplexer_
    examples.batch`, no script edits — nothing touched since 08:41 changes
    `modal_multiplexer_hybrid.py`'s or the solvers' cache fingerprints, so a
    clean disk-cache hit on sections 0-3 is expected). Remaining real work is
    the 3 section-4 FDTD solves; comparable single-section FDTD solves earlier
    in this same log took 5-33 min, so the 5h budget should now be generous
    rather than tight. Watch `reports/11590663_0.out` — if it doesn't reach
    `example_output/modal_multiplexer/modal_multiplexer_hybrid/` output within
    the first ~15 min, the cache did not hit as expected and this entry needs
    another look.
  - **❌ TIMEOUT again, 2026-07-31 18:58, full 5h.** `sacct`: `State=TIMEOUT`,
    `ExitCode=0:0`, same as the first attempt. The "generous rather than tight"
    estimate above was wrong — the 3 section-4 FDTD solves were not a short
    tail. `example_output/modal_multiplexer/modal_multiplexer_hybrid/sections/`
    confirms all 3 per-port section-4 geometry plots exist (`port_top`,
    `port_bot`, `port_mid`), so the disk-cache hit on sections 0-3 worked as
    predicted and this run did get into section 4 — but the log's tail is deep
    inside a per-port eigenmode-overlap monitor computation (repeated `MPB
    solved for frequency_13(...)` lines sweeping many k-points, then starting
    band 14) when SLURM killed it, not writing final S-matrix output. No new
    `example_sim_cache/` entries appear for this device in the run's own
    13:58–18:58 window (the cache dirs that changed in that window all belong
    to other jobs' MMI runs sharing the login node's cache directory), so
    section 4's per-port solves are not yet caching intermediate progress the
    way sections 0-3 apparently do — a plain resubmit will very likely restart
    section 4 from scratch rather than resume mid-monitor-sweep. **Gap 6's
    original ask — get `warn_on_transfer_duplicate` to fire or not fire on
    this device's real interfaces — is still unanswered**, now after two full
    5h attempts (10h of 4-node time total). Worth reshaping before a third
    blind resubmit: either raise the wall-clock limit substantially (this
    device's section-4 monitor cost was underestimated twice now), or split
    section 4's 3 per-port solves into their own array job so a future timeout
    loses at most one port's progress instead of all three. Not resubmitted as
    of this snapshot — see `modal_multiplexer_validation.md` gap 6 for the
    validation-doc-side reasoning.

### 11592462 `[0]` — `beef_dc_examples` (directional_coupler_hybrid) — ❌ TIMEOUT, diagnosed 2026-07-31 late
- **Script:** `batch/run_directional_coupler_examples.batch` → `make
  directional_coupler_hybrid` (array index 0 of the three coupler examples;
  `_1` mode and `_2` fdtd both completed separately, see table above).
- **❌ TIMEOUT 2026-07-31 21:34, full 5h.** `sacct`: `State=TIMEOUT`,
  `ExitCode=0:0`.
- **Diagnosed, full read of `reports/11592462_0.out` (all 512 lines).** The
  two `Initializing structure...` blocks (lines 22 and 330) are **not**
  `M_in0`/`M_in1` (the script's own "2 FDTD solves" cost estimate) — they
  don't match that shape at all. A real `solve_section()` call (confirmed
  against 11591249's successful log) prints ONE small mode-source block
  then goes straight into `on time step` lines; these two blocks instead
  each solve **every wavelength** (11, then a 1-wavelength repeat) × 5 bands
  before anything else — that is `section_mode_neff()`'s own documented
  shape ("one structure, every wavelength... seconds, not minutes"). Cross-
  referencing `directional_coupler_hybrid.py`'s source: `CALIBRATE_BODY =
  True` (hardcoded, not env-gated) runs `section_mode_beta()` first, whose
  *seed* step is exactly one such call (N_MODES=5, all 11 wavelengths) —
  that's block 1 (lines 22-328). The script's separate single-wavelength
  `_mpb = section_mode_neff(...)` diagnostic (used only for a comparison
  print line) is block 2 (lines 330-376, 1 wavelength since
  `wavelengths[_i:_i+1]`) — its first seed value (`frequency_1(1.6129,...)`)
  exactly matches block 1's 6th (middle, λ=1.55) repeat, confirming it's the
  same cross-section at the band-centre wavelength, not a new mystery.
  **Both blocks completed successfully** (block 2 reaches "Dominant
  planewave for band 5" at line 376). **What's missing: the 15 real
  `solve_section()` calibration runs `section_mode_beta()` is supposed to
  do in between** (3 lengths × min(5,5) excited modes) — there is a bare
  `-----------` divider between blocks 1 and 2 with zero evidence any of
  them ran (a real `solve_section()` call always leaves a trace, per the
  11591249 comparison above). Two readings: either all 15 hit the on-disk
  cache (`save_dir=save_dir` is passed) and returned near-instantly, or the
  job never got that far before block 2 started — either way, **the
  documented-cheap `section_mode_neff()` seed/diagnostic calls, not the
  documented-expensive `solve_section()` calibration runs, consumed
  essentially the entire 5h budget.** That contradicts
  `section_mode_beta()`'s own docstring by orders of magnitude at this
  device's production N_MODES=5 / 11-wavelength / `Computational cell is 3 x
  7 x 2.5` configuration — not diagnosed further (whether it's a full-3D-grid
  MPB cost rather than a cheap 2D cross-section, or something else). The
  final 120 `Warning: grid volume` lines (right after block 2 finishes, no
  MPB solves interleaved) most likely mark the start of whatever comes next
  in the script (the real calibration solves and/or `M_body`'s 200-length
  EME sweep) being killed before any could report output — not independently
  explained.
- **Recommendation: do not resubmit this whole script blind a third time.**
  Before resubmitting, benchmark `section_mode_neff()`'s cost in isolation
  at this device's real N_MODES=5/11-wavelength/domain-size configuration
  (a short, cheap diagnostic job, the same pattern as this session's other
  `diagnose_*` scripts) to find out whether that step is really the
  multi-hour cost here — if so, `CALIBRATE_BODY`'s seed computation itself
  needs reshaping (e.g., fewer wavelengths for the seed, since it only needs
  to pick the right 2π branch, not a converged answer) before this job can
  ever finish in one slot. Not resubmitted as of this writing.
- **Benchmarked, 2026-07-31 22:05 — theory REFUTED.**
  `diagnose_dc_calibration_seed_cost.py` (new) measures
  `section_mode_neff()` directly at this exact configuration (N_MODES=5,
  same cross-section, same `Computational cell is 3 x 7 x 2.5`). Login-node
  serial smoke (1 process): **28.7 s/wavelength.** Same script under the
  production resource shape, job **11597622** (2 nodes, 48 ranks, 3
  wavelengths): **44.4 s/wavelength** — MPI rank overhead makes it *slower*
  per-wavelength, as suspected, but only ~1.5×, not the 30-40× that would be
  needed to explain a multi-hour cost. **Extrapolated full cost (11
  wavelengths, production MPI config): ~8 min.** That directly refutes the
  "cheap eigensolve is secretly expensive here" theory above — whatever ate
  the ~4.5h between the two `section_mode_neff`-shaped log blocks in
  11592462_0, it was NOT this call itself.
  **Root cause is back to unknown.** Given seed+comparison together cost
  ~9 min and the 15 real `solve_section()` calibration runs left literally
  zero log trace (no `Initializing structure` between the two blocks),
  those 15 runs most likely all hit the on-disk cache (near-instant, which
  also prints nothing) — meaning **essentially the entire accounted-for
  workload in this script should complete in well under 15 minutes**, yet
  the real job ran 5 hours. Something consumed the other ~4h45m with no
  log signature at all (no print, no crash, no partial structure-init).
  Live open questions, none checked yet: whether disk I/O contention
  (concurrent jobs hitting the same `example_sim_cache`/`save_dir`) stalled
  a cache lookup or write; whether `M_body`'s 200-point
  `solve_section_length_sweep()` (never confirmed to have started) is
  silently expensive in a way that also prints nothing; or whether the job
  was simply stuck/hung (cf. job 11586046's unrelated hang earlier this
  session, also silent). **Next step: resubmit with a shorter time limit
  and watch it interactively** (`squeue`/log tail) rather than diagnosing
  from a dead job's log alone — a live job can be interrupted and inspected
  mid-run, which the two theories above cannot be told apart without.

### 11591249 — `beef_mmi_2x2_engine_agreement` — ❌ TIMEOUT
- **Script:** `batch/verification/run_mmi_2x2_engine_agreement.batch` →
  `examples/verification/python/verify_mmi_2x2_engine_agreement.py` (new).
- **Resources:** 2 nodes × 24 ranks, 5h limit — same shape as this device's
  own production hybrid/FDTD jobs, since this is real FDTD time-domain work,
  not an eigensolve.
- **What it does:** Solves a uniform stretch of the 2×2 MMI body at 4 lengths
  (2/4/6/8 µm — deeper margin than the coupler's minimum two, since this
  body has more room) and checks FDTD and EME agree on each mode's
  accumulated phase via the shared `dbeta_from_phases` helper — the same
  2-mode pair gap 2's own Δβ_01 measures. This is check 1 (engine agreement)
  for the 2×2, the device's last unmeasured check.
- **Why:** answers whether this device's own hybrid-vs-FDTD gap (gap 5:
  0.82% max\|ΔT\|, sign flips across the band) has the coupler's
  Δβ-discretization mechanism underneath it or something else —
  `calibrate_body_neff` has never been tried on this device.
- **❌ TIMEOUT 2026-07-31 19:46, full 5h.** `sacct`: `State=TIMEOUT`,
  `ExitCode=0:0`. Log shows all 4 lengths' first FDTD run completed (`run 0
  finished` appears 4 times, at 29,398 / 30,997 / 33,334 / 30,505
  timesteps), then a 5th run began and was only at timestep ~5,311 (of a
  likely ~30,000-step run, by the pattern of the first four) when killed —
  under a fifth of the way into whatever this run is (a second pass per
  length, or a reference run; the script's own structure was not re-read to
  confirm which). **Gap 4 — this device's last open check — did not close.**
  Net: the job got roughly halfway through its total FDTD work in one 5h
  slot; a straight resubmit (no cache to fall back on here — this script has
  no `save_dir` reuse the way the modal mux's does) would need to redo all of
  it, or the script needs checkpointing / a longer wall-clock budget before
  trying again. Not resubmitted as of this snapshot — see
  `resolved/mmi_2x2_validation.md` gap 4 for the validation-doc-side reasoning.

### 11596845 `[0-3]` + 11596851 — `beef_mmi_2x2_engine_agreement` (reshape of 11591249) — ✅ COMPLETED, GAP 4 CLOSED
- **Scripts:** `batch/verification/run_mmi_2x2_engine_agreement.batch`
  (`--array=0-3`, `-t 480`) and `…_merge.batch`, both →
  `verify_mmi_2x2_engine_agreement.py`. `SLURM_ARRAY_TASK_ID` picks one
  length; the merge step runs with none set and re-reads all four from the
  on-disk cache under `cache_only=True`.
- **Outcome, 2026-08-01:** `_0` (L=2) and `_1` (L=4) **COMPLETED in 36 s
  each** — the predicted cache hit on 11591249's already-computed lengths,
  which is the one thing that attempt did leave behind. `_2` (L=6) 3:35:10,
  `_3` (L=8) **5:07:14**. Merge **11596851** COMPLETED 03:07 in 25:29 on 1
  rank; no `CacheMissError`, so every piece was where it was supposed to be.
- **✅ Result — gap 4 closed:**

  | quantity | value |
  |---|---|
  | mean ratio Δβ_FDTD / Δβ_EME | **0.9953** (0.47% low) |
  | across 1.530–1.565 µm | 0.9962 → 0.9945, **monotone** |
  | worst \|off-diagonal\| of T | 1.21e-05 |
  | worst phase-fit residual | 4.02e-06 rad |

  Against the coupler's ~10% on the identical measurement this is **21×
  tighter**, so `calibrate_body_neff` is **not indicated** for the 2×2 — a
  negative result, and the one gap 4 existed to get. Saved to
  `example_output/verification/mmi_2x2_engine_agreement/engine_agreement.npz`
  (`beef_git` 9e52e9d, res 40).
- **⚠️ Read the log carefully here:** it prints both diagnostics as `0.0000`
  because the script formats them `:.4f`. They are *not* zero and *not*
  defaulted — the `.npz` was opened and the real values (1.21e-05, 4.02e-06)
  read out before this was recorded. An exact-looking zero in a merge step
  that could in principle skip work is worth checking rather than quoting.
- **Why it does not explain gap 5:** the bias here is monotone across the
  band; gap 5's ΔT flips sign across that same band. A monotone Δβ bias
  cannot produce a sign-flipping ΔT, so gap 5's 0.82% comes from elsewhere
  (tapers or the cascade — neither isolated). Argument from signatures, not a
  measurement.
- **Cost datum worth keeping:** 96 / 183 / 215 / 307 min for L = 2/4/6/8 —
  13.4 h serial, i.e. the original single 5 h slot was never going to make it
  and L=8 alone exceeds it. Folded into `resolved/mmi_2x2_validation.md` gap 4 and
  `README.md`.

### 11597177 `[0-2]` + 11597199 — `beef_modal_mux_examples` (reshape of 11590663) — ❌ THIRD TIMEOUT, cause found, NOT resubmitted
- **❌ All three array tasks TIMEOUT** 2026-08-01 02:44 / 02:52 / 03:01, each
  the full 5 h on 4 nodes. **11597199** (merge) never started —
  `DependencyNeverSatisfied`; inert, cancel it.
- **Root cause — the reshape dropped `srun`.** The batch script invokes
  `python modal_multiplexer/python/modal_multiplexer_hybrid.py` directly to
  stop `make` from running all 3 ports in one process, but the Makefile target
  is `$(runner) $(my_python) …` with `runner = srun --mpi=pmi2` on Phoenix
  (`BEEF/examples/Makefile:5`). Losing `make` lost the launcher:
  `Using MPI version 4.1, **1 processes**` in all three new logs vs
  `96 processes` in 11585847_0 and 11590663_0 — **0.27–0.36 s/step against
  0.0064, ~42× slower**, with 96 cores held idle (`cput=20-00:35:12` for
  `walltime=05:00:22`).
- **Second flaw, independent of the first:** the premise that sections 0-3
  were already warm on disk is **false** — every task re-ran the same three
  `Initializing structure`/MPB/time-stepping blocks at near-identical log
  offsets (18/2189/4528, 18/2184/4523, 18/2224/…). Not established whether
  the earlier timeouts never persisted them or the cache key differs between
  96-rank and 1-rank runs. Tasks died at unrelated points (`_0` timestep 289,
  `_1` timestep 2945, `_2` inside MPB band 10), consistent with all three
  grinding through that duplicated prefix.
- **Not resubmitted — the mux is deferred.** Recipe if resumed, in order:
  restore `srun --mpi=pmi2`; verify sections 0-3 actually hit cache (a task
  that skips them should be time-stepping within a minute or two, not an
  hour); then resubmit. **Do not just raise the wall clock** — at 1 rank no
  affordable limit finishes. Full write-up in
  `deferred/modal_multiplexer_validation.md` gap 6.
- **Note the pattern is not at fault:** the same array-plus-`cache_only`-merge
  design submitted in the same window worked first try on the 2×2 (11596845
  above). This broke on the launcher.

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
  `resolved/fdtd_vs_EME_mode_validation.md` item 12, with the caveat that whether (10,11)
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
  SOI_220 device on the 2.8 µm body that `resolved/fdtd_vs_EME_mode_validation.md` item 11
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
  filtering. Item 13 closed; folded into `resolved/fdtd_vs_EME_mode_validation.md` item
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

### 11588997 — `beef_dispersive_convergence` (replaces 11586986) — ✅ COMPLETED
- **Same script, post-fix.** `reports/11588997.out`, 1:56:53.
- **Result:** at band centre (1.55 µm) dispersion changes nothing — identical
  convergence resolution to the non-dispersive column at every threshold
  (1e-3/3e-4/1e-4/3e-5 all first clear at res 60/100/100/100 either way). At
  the band edge (1.60 µm, where dispersion is not near-inert by construction),
  the 1e-3 threshold moves from res 60 to res 80 under dispersion — but item
  11 already pins production's EME resolution at 80, so the dispersive case's
  actual requirement lands exactly on what already ships; nothing currently
  deployed was ever run at res 60. The tightest threshold (3e-5) is the one
  gap: non-dispersive clears it everywhere by res 100, dispersive never drops
  under it at the edge inside the tested range (20–120) — though the raw
  dispersive-vs-non-dispersive Δn_eff difference oscillates at ±(3–5)e-4 with
  no resolution trend, reading like noise floor rather than a growing effect.
- **Verdict: the proposed node-day full-device sweep (11 array tasks × 2 nodes
  × 24 ranks, up to 24h) is not warranted** at the 1e-3-class tolerance item 11
  actually used to set the shipped resolution. Full table in
  `material_model_validation.md` work item 4, now closed for the mode-solve
  half.

### 11589770 — `beef_mmi_1x2_hybrid_calibrated`
- **Script:** `batch/run_mmi_1x2_hybrid_calibrated.batch` (new) →
  `MMI1X2_CALIBRATE_BODY=1 make mmi_1x2_hybrid`.
- **Resources:** 2 nodes × 24 ranks, 300 min — same shape as the uncalibrated
  hybrid array task (11585845_2, 1:46:59), since `calibrate_body_neff` changes
  how the body section's phase is propagated, not how many FDTD sections run.
- **Why:** `resolved/mmi_1x2_validation.md` gap 5 measured the hybrid 9.11% **under**
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

### 11590785 — `beef_mmi_2x2_gap2_dbeta`
- **Script:** `batch/verification/run_mmi_2x2_gap2_dbeta.batch` (new) →
  `examples/verification/python/verify_mmi_2x2_gap2_dbeta.py` (new).
- **Resources:** 1 node, 8 cpus, serial, 180 min budget. Eigensolve only —
  `EMESolver.compute_modes()` at `num_modes=2`, no FDTD, no time-stepping.
- **Why:** `resolved/mmi_2x2_validation.md` gap 2's existing measurement (job 11500536)
  shows the MMI body's fundamental n_eff still moving at res 100, but that
  doc explicitly declines to conclude anything from an arbitrary quantity —
  same "sweep the wrong thing, hide the real error" lesson item 11 taught on
  the coupler (there: T_cross vs Δβ). For a 2×2 MMI the physically load-bearing
  quantity is Δβ_01 = 2π/λ·(n_eff[0]−n_eff[1]) between the two lowest guided
  modes — the beat length `mm_length=28.0` was chosen against — not either
  mode's own n_eff. This sweeps res 40/60/80/100/120 measuring Δβ_01 directly,
  the same pattern item 11 used for the coupler's even/odd pair.
- **Cost note:** a login-node smoke check (`OMP_NUM_THREADS=1`, res 40 only)
  measured 205.8 s for one point — slower than expected for a 2-mode solve,
  hence the batch submission (8 threads) rather than running the full sweep
  on the login node. 180 min budgets generously for the higher-resolution
  points, whose mesh cost does not scale linearly with `res`.
- **✅ Result (1:05:48).** Δβ_01 steps: −0.00092/−0.00049/−0.00020/−0.00008
  rad/µm at 40→60/60→80/80→100/100→120 — 0.56%/0.30%/0.12%/0.049% relative.
  **Gap 2 closes at res 100** for this device: the physically load-bearing
  quantity is converged to well under 0.1% there, regardless of absolute
  n_eff's own slower-looking tail. One correction made to this session's own
  earlier partial reading: comparing Δβ_01's raw step size directly against
  n_eff[0]'s raw step (different units/magnitudes, ~0.16 vs ~2.8) and calling
  it "an order of magnitude tighter" was a mixed-units comparison — read
  relative to each quantity's own value, the two converge at a *similar*
  rate; what actually matters is that Δβ_01 itself is small enough to trust
  by res 100, not how its convergence rate compares to another quantity's.
  Full table in `resolved/mmi_2x2_validation.md` gap 2.

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
