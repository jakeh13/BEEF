# 4×4 MMI 90° optical hybrid — SOI 220 nm, O-band — build plan

Opened **2026-08-05**, at the user's request: an optical 90° hybrid based on a
4×4 MMI, on **220 nm silicon-on-insulator** (`SOI_220`), with real **input and
output access tapers**, reproducing the device topology of

> T. Yasui, J. Sugisaka, and K. Hirayama, "Structural Optimization of an
> Optical 90 Degree Hybrid Based on a Weakly-guided 4×4 Multimode Interference
> Coupler Using a Parallelized Real-coded Micro-genetic Algorithm," *ACES
> Journal*, vol. 36, no. 5, pp. 526–532, May 2021.
> https://doi.org/10.47037/2020.ACES.J.360506

Full citation and the two papers it in turn cites that matter here are in
[References](#references) at the bottom.

> ## ⚠️ CURRENT STATE — read this before any older section (2026-09-13)
>
> This document is ~6000 lines and largely chronological, so **most CMRR
> numbers in it are superseded.** The short version:
>
> * **The device MEETS the CMRR spec.** `mm_length = 158.80 µm`, min CMRR
>   **20.326 dBe** (spec ≥ 20), max |Δφ| 3.664° (spec ≤ 5), excess loss
>   2.444 dB. Every earlier figure — 18.482 dBe, 15.615 dBe, and the long
>   "fails both specs" era before that — was computed with a body propagation
>   constant 0.0458 too low, at a length chosen for that wrong constant.
> * **n_eff calibration is MANDATORY here, not a refinement.** The mode solver
>   and meep's time-stepper realise different β on the same cross-section, and
>   the gap only closes as ~res^-1.45 — reaching 0.005 would need res ≳ 180,
>   infeasible in 3D.
> * **A calibration's resolution must MATCH the simulation it reproduces.**
>   Calibrating at res=60 against a res=40 reference made agreement 18× WORSE,
>   by exactly the predicted amount. Never mix.
> * **Hybrid vs full-device FDTD now agrees to 2.12° power-weighted**, and the
>   FDTD reference's own phase noise is 0.029°, so that is a real number and
>   not a floor artifact.
> * Closed by measurement: calibration resolution (negative), interface mode
>   count 16→24 (negative), reference `decay_tol` (not the limit). Leading open
>   candidate: the number and span of the **calibration lengths**.
>
> Jump to the dated entries at the END of this file for the full trail.

> **This doc is a PLAN, not a gap list** — same convention as
> `mmi_1x4_sin300_validation.md` and
> `high_contrast_grating_fabry_perot_validation.md`, the two other build
> plans in this directory: nothing below has been simulated, no geometry
> exists yet, and every device *number* in Part B is either taken directly
> from the paper's own general design equations or explicitly marked
> **TBD — measure, don't assume**. Findings about the *codebase* (what
> already supports this device, what doesn't) are labelled **MEASURED (read
> from source)**, exactly as the other two plans do. This converts to a
> normal gap list against the [four checks](README.md#the-four-checks) once
> the example scripts exist and run for the first time.

**⚠ Terminology collision, flagged once so it doesn't cause confusion for
the rest of this file.** "Hybrid" means two unrelated things in this
project:

1. **The device under test** — an *optical 90° hybrid*, the coherent-receiver
   component this whole paper is about.
2. **`HybridSolver` / `*_hybrid.py`** — this codebase's FDTD+EME cascade
   *simulation method*, used for every device in this directory regardless
   of what the device itself is called.

Every other file here is implicitly "the hybrid-method validation for
device X." This file is "the hybrid-method validation for a device that is
*itself* called a hybrid." Read "hybrid" as the simulation method unless a
sentence explicitly says "optical hybrid" or "90° hybrid."

---

**CURRENT STATUS (2026-09-01): the active work is the (mm_width, y_gap)
optimization campaign — jump to
[The (mm_width, y_gap) optimization campaign](#the-mm_width-y_gap-optimization-campaign--2026-08-26--2026-09-01)
at the bottom of this file. Best geometry so far: `mm_width=6.50`,
`y_gap=0.56`, `mm_length=91.0` — phase spec MET at 4.764°, CMRR 12.109 dBe
still short of the 20 dBe bar. The 2026-08-09 status block below is retained
as history.**

**CURRENT STATUS (2026-08-09) — read this first, history follows below.**
Steps 1-7 done. Step 6 (all 3 example scripts) done and RE-RUN at a
properly-converged body mode basis -- see "Body/taper mode-count
convergence investigation" (before Step 8) for the full trail. **The
numbers that matter right now:**

- `N_MODES_TAPER=12` (mmi_4x4_90hybrid_eme.py), `N_MODES_BODY=16`
  (body_sweep.py/hybrid.py **and, since 2026-08-21, the EME script too**)
  -- both converged, both measured, not guessed. `num_cells=30` for the
  EME-only script's taper (job 11736983).
  > ⚠ **The "Step 6 (all 3 example scripts) done and RE-RUN at a
  > properly-converged body mode basis" claim above was FALSE for 12 days,
  > and the parenthetical "(body_sweep.py/hybrid.py)" on this very line is
  > what gave it away.** The 2026-08-09 redo reached exactly those two
  > scripts. `mmi_4x4_90hybrid_eme.py` and
  > `verification/python/plot_mmi4x4_90hybrid_eme_fields.py` were both
  > missed and sat at the superseded `N_MODES_BODY=8` / `ΔL=-1.50` /
  > `mm_length=80.1896` point -- the basis this doc's own convergence
  > check measured as 50-64% under-converged, at a design point this
  > project had already rejected. Found 2026-08-21 by checking artifact
  > mtimes against the redo date, not by reading the scripts. Both
  > corrected; see "EME script staleness" below.
- Settled geometry: `mm_width=6.0`, `waveguide_width=0.35`,
  `waveguide_offset=1.5`, `taper_width=0.8`, `taper_length=1.0`,
  **`ΔL=-2.75` (mm_length=78.9396)** -- re-derived from the old
  `ΔL=-1.50`/`mm_length=80.1896` once the body basis converged; the old
  point is stale, do not reuse it.
- **Real Step 6c (hybrid, FDTD+EME) numbers at the converged point (job
  11778128): worse than previously reported.** CMRR fails the cited
  ≥20 dBe bar at most wavelengths (was: comfortable at band centre);
  phase error breaches the cited ≤5° bar at multiple ports (was: only
  one); excess loss 6.5-8.4 dB peak (was 2.1-4.3 dB). The earlier
  "close to spec" read was an artifact of the under-converged N=8 body
  basis, not real performance -- corrected, not regressed.
- Steps 8 (four checks) and 9 (golden fixture) not started.
- **Open, not yet started: whether the geometry itself (mm_width, taper
  dimensions, access waveguide width) needs redesigning to actually meet
  the paper's cited spec, now that the real performance is known.**
  Paused here at the user's direction (2026-08-09) -- pick this up
  next if/when resumed.
- **2026-08-14, mm_width Stage-1 scoping — job 11840371 TIMED OUT at -t
  480 (8h), partial real results; continuation submitted as job 11855837.**
  Cheap L_π-vs-width + guided-mode-count probes (9 candidate widths,
  3.5–10.0 µm, `mmi_4x4_90hybrid_mm_width_scope.py`) to pick 1-2 promising
  widths before committing to the expensive Stage-2 full treatment (own
  body solve + baseline 4-solve per width, no cache reuse possible).
  **L_π-vs-width fully completed and clean**: `L_π/W²` stays 3.023–3.026
  µm⁻¹ across the ENTIRE 3.5–10.0 µm range — the same self-imaging scaling
  already confirmed at 6–16 µm now holds down to the narrowest
  geometrically-valid width too, not just near the current design point.
  **Mode-count probe (the expensive half) only reached 3 of 9 widths
  before the time limit**: 3.5 µm → 23 guided modes (6434 s), 4.0 µm → 24
  (7051 s), 4.5 µm → 24 (4774 s) — **~1.3–2h per width, not cheap despite
  the "Stage 1" framing, and going narrower does NOT reduce mode count
  substantially** (23-24 modes at 3.5-4.5 µm vs. the current 6.0 µm
  design's own already-known >16). Real duplicate-mode collisions logged
  at every width (`meow mode N duplicates already-accepted mode N-1`,
  `fdtd_vs_EME_mode_validation.md` item 0's own guard firing, not a new
  bug) — this cross-section is richly multimode even at the narrow end.
  **Continuation submitted as job 11855837 — ✅ COMPLETED 2026-08-15
  (14:41h wall).** Full mode-count picture, all 9 widths:

  | width (µm) | guided | wall time |
  |---|---|---|
  | 3.5 | 23 | 6434 s |
  | 4.0 | 24 | 7051 s |
  | 4.5 | 24 | 4774 s |
  | 5.0 | 24 | 5078 s |
  | 5.5 | 24 | 4931 s |
  | 6.0 (current) | 24 | 7427 s |
  | 7.0 | 24 | 10351 s |
  | 8.0 | 24 | 10407 s |
  | 10.0 | 24 | 14651 s |

  **Sobering, load-bearing result: every width from 4.0 to 10.0 µm hits
  the full `n_request=24` CAP** — meaning the true guided-mode count is
  `≥24` almost everywhere in the tested range, this probe only ever
  measured a floor, never found a ceiling. Wall time trends upward with
  width (roughly 5-7k s at 4.5-6.0 µm to ~10-15k s at 7-10 µm), consistent
  with more/costlier duplicate-mode collisions to search past at wider
  cross-sections. **There is no evidence, across the entire geometrically-
  valid range (3.2 < W ≤ 10), that any candidate `mm_width` needs fewer
  modes — hence is any cheaper to fully converge — than the current 6.0 µm
  design.** Combined with `waveguide_offset` already being closed off
  (above), the geometry-redesign search inside this (offset, mm_width)
  family is not finding a cheap lever; the remaining untried candidates
  from the 2026-08-09 open question are the still-unmeasured `taper_width`/
  `taper_length` placeholders, or accepting the current CMRR/phase gap as
  a real design limitation rather than a bug to fix.
- **2026-08-14, `calibrate_body_neff`/`section_mode_beta()` investigation
  — real crash found, HANDED OFF to a separate, already-in-progress
  effort, not pursued further in this thread.** This device has never
  used FDTD-consistent body-neff calibration (the mechanism `dc_
  validation.md` item 6 introduced — MPB's body n_eff and what meep's
  FDTD actually accumulates differ by ~2.3%/mode, amplified by Δβ into a
  ~10% coupling error on the coupler). A direct attempt to replicate
  `HybridSolver`'s own calibration for this device's hand-rolled cascade
  hit a real `RuntimeError: meep: simulation fields are NaN or Inf` inside
  `section_mode_beta()`'s own internal FDTD call — traced to the method's
  branch-safety length formula picking an absurdly short calibration
  domain (~1.9 µm total, 53% of it PML) for this body's high `β_max`
  (`n_eff≈2.96`), a real edge case for wide/high-index-contrast multimode
  bodies the method has apparently never been exercised against (the only
  working precedent, `mmi_2x2_hybrid.py`, has a much narrower `mm_width=
  2.8` body). Not a usage mistake — the call matches that working
  precedent exactly. **Not resolved here**: the user has a separate agent
  already working this implementation; this thread stood down to avoid
  duplicating cluster time/edits.
- **2026-08-17, `section_mode_beta()` calibration RESUMED at the user's
  direction — real fix found and confirmed, one new flag surfaced.** A
  separate, unrelated shared-library effort completed 2026-08-12→16
  (`resolved/section_mode_beta_multi_monitor_validation.md`) added
  `use_merged_domain=True` (now `section_mode_beta()`'s default): ONE FDTD
  domain sized to `max(lengths)` with multiple monitors, instead of
  `len(lengths)` fully separate domains. Reasoned root cause of the
  original crash: it almost certainly happened in the SHORTEST individual
  per-length domain (0.885 µm + 2×0.5 dpml = 1.885 µm, 53% PML) — the
  merged path builds to `max(lengths)`=2.655 µm instead (27% PML) and so
  never constructs that domain at all. Caveat: the merged path silently
  falls back to the old per-length path whenever `save_dir` is passed (no
  cache-key design yet) — must be called without it to actually engage.
  **Retested directly (same body cross-section, `save_dir` dropped,
  reduced test settings) — ✅ SUCCESS.** FDTD decayed completely normally
  (`t=111.5`, vs. the original crash's non-decaying `t≈4228`), all 16
  modes' β computed, and at LOWER cost than the crashed attempt (5488 s vs.
  13853 s). Root-cause diagnosis confirmed correct, not just plausible.
  **New flag, not yet resolved**: the code's own built-in sanity check
  fired — `"FDTD β exceeds the MPB seed by +13.69% at λ=1.26 µm mode 15 —
  out of family with the +2%-ish discretization gap this calibration
  exists to measure"`. Mode 15 (the highest of 16 requested, likely the
  most weakly-guided/near-cutoff mode) shows a gap ~7× the size this whole
  mechanism is designed to correct for. Not yet distinguished whether this
  is real physics at this device's dense, near-degenerate mode spectrum,
  or an artifact of the reduced test settings (`res=15`, not production
  `res=40`). **Next step (not yet run): full production-settings
  calibration (`res=40`/`decay_tol=1e-11`/`N_WL=3`)**, then use the
  calibrated neff to re-check the current baseline design's CMRR/phase
  numbers against the uncalibrated production run (job 11778128) — the
  original motivation for this whole thread. **HOLDING, then UNBLOCKED,
  2026-08-17.** The original 3-9+ day estimate (naive res-scaling) turned
  out to be wrong: a real, measured shared-library fix
  (`resolved/section_mode_beta_multi_monitor_validation.md` step 9, same
  day) found `section_mode_beta()`'s dominant cost (89%, measured on job
  11872191) is MPB eigenmode solves, which do **not** scale with meep
  resolution at all — only the FDTD-stepping minority (11%) does. Two
  fixes now apply automatically, no code change needed:
  `measure_diagonal_only` (~4×) and the merged-domain path no longer
  falling back to per-length just because `save_dir` is set (~1.9×
  further, confirmed directly in `beef/solvers/fdtd.py`:
  `merged_active = use_merged_domain`, no `save_dir`/`cache_only` gate
  anymore). Real cost model: `n_modes×(n_lengths+1)×n_wl×14s +
  n_modes×n_wl×14s + n_modes×timesteps×0.0089s` — for this device
  (N=16, 3 lengths, N_WL=3), **~90 min for calibration alone**, not days.
  **Submitted at full production settings, job 11938100**
  (`run_mmi_4x4_90hybrid_calibrated_baseline.batch`, `-t 360`/6h budget,
  real margin over the ~3h total estimate — calibration + the body's own
  now-calibrated EME re-solve, reusing job 11778128's 4 cached input
  solves unchanged).

  **❌ TIMED OUT at the full 6h — never reached the calibration
  measurement.** Root cause diagnosed: the script's own `cache_only=True`
  fetch of the 4 baseline input matrices MISSED (job's own log: "input
  matrices: cache MISS (expected at non-production settings) — solving
  fresh", printed despite this being a production-settings run — a real
  signal, not a benign message) — the entire 6h budget went to redoing
  those 4 already-cached (job 11778128) solves from scratch, with the
  actual `section_mode_beta()` calibration only starting in the run's
  final moments (1 MPB solve logged) before the time limit hit. **The
  underlying cache-fingerprint mismatch has since been fixed** by the
  other agent's ongoing `fdtd.py` work — re-ran the identical
  `cache_only=True` fetch immediately after with current code: all 4 HIT
  instantly (0.25 s total, vs. the 6h the stale-code job burned). A fresh
  resubmission now would skip straight to the real ~90 min calibration
  measurement. **Resubmitted at the user's direction, job 11949712**
  (same batch script, no changes needed — the script itself was correct,
  the fix landed in the shared library). **Job 11949712 stuck PENDING on
  `AssocGrpBillingMinutes` again** (same class of issue as Step 6b's own
  task-1 delay, 2026-08-14) — cancelled and resubmitted under the plain
  `gts-sr108` account as **job 11981048** (`sbatch --account=gts-sr108`),
  confirmed picked up (reason `Priority`, normal queue wait), started
  2026-08-18 10:25:56.

  **❌ TIMED OUT AGAIN at the full 6h — a different, more concerning
  failure this time.** The cache-fingerprint bug is confirmed fixed
  (`input matrices: cache HIT (all 4, no new FDTD)` in this run's own
  log). But the calibration itself hit a severe, unestimated cost: only
  ONE `Initializing structure...` block appears after the MPB seed phase
  — after the full 6h, the job was still stuck inside the FIRST of 16
  excited-mode calibration runs, at meep time≈8303 (660,245+ timesteps),
  still stepping cleanly (no NaN, no crash) but nowhere near decaying.
  The reduced-settings smoke test of this exact same domain (res=15,
  decay_tol=1e-9) decayed cleanly by t=111.5 — this production run is
  ~74× further in simulated time alone and still not converged. **This
  is the same failure signature as Step 6b's full-device timeouts**: a
  slow, non-decaying (not crashing) run, plausibly tied to this device's
  dense, closely-spaced mode spectrum interacting badly with
  `decay_tol=1e-11`'s strictness — not simply "needs more wall-clock
  budget," the same lesson Step 6b already taught twice. The earlier
  ~90 min calibration estimate (step 9's cost formula) assumed a
  `timesteps` value that does not hold for THIS device's first excited
  mode specifically — the formula itself may be fine for a well-behaved
  mode; this device's mode 1 is evidently not one. **Not resubmitted** —
  treating this the same as the user's standing Step 6b policy rather
  than assuming a bigger `-t` fixes it.
- **2026-08-18, investigating WHY mode 1 stalls — real findings, not yet
  a root cause.** Two checks, both at zero/near-zero new cluster cost:
  1. **Mode spacing is NOT anomalously close.** Pulled the MPB seed's own
     band-1..16 values from job 11981048's log: the gap between band 1
     and band 2 (0.004370, in the same units as the printed k-vector) IS
     the smallest gap in the whole spectrum, but the full sequence of
     gaps grows smoothly and monotonically (ratio of successive gaps
     drops smoothly from 1.67 down to 1.09 by band 16) — textbook
     dispersion for a wide multimode slab where the lowest-order modes
     are always most closely spaced near the bulk index. **This is the
     expected shape, not evidence of true pathological near-degeneracy**
     — retracting the initial hypothesis that mode 1 has an anomalously
     close partner.
  2. **Directly compared bands 1, 2, 3's own decay behavior at matched
     reduced settings** (res=15, decay_tol=1e-9, same merged-domain
     mechanism, explicit solver construction — not via env var, after a
     first attempt accidentally ran at res=40 because
     `mmi_4x4_90hybrid_offset_sweep.py`'s env var name
     (`MMI4X4_OFFSET_RES`) doesn't match `MMI4X4_CAL_RES`, caught before
     wasting real time on it). **Band 1 decayed FASTEST of the three**
     (t=108.0, 3241 steps), not slowest (band 2: t=114.7/3441 steps;
     band 3: t=150.7/4521 steps) — the opposite of what an intrinsic
     mode-1 pathology would predict. All three also triggered the same
     "near-degenerate excited mode" normalization-mismatch guard at
     essentially identical magnitude (36.4-36.8%), so that warning isn't
     singling out mode 1 either — it fires generically on this basis.
  **Conclusion so far: the production-scale stall on mode 1 is very
  unlikely to be an intrinsic property of that physical mode** — reduced
  settings show it decaying cleanly and fastest of the three tested.
  The failure most likely depends specifically on production resolution/
  `decay_tol`, not on which mode is excited per se. **Two live,
  undistinguished candidates**: (a) a genuine numerical effect (marginal
  reflection/resonance) that only appears at finer discretization or
  stricter tolerance and gets damped out at reduced settings; (b) mode
  reordering — "band 1" at res=15 is not guaranteed to be the same
  physical mode as "band 1" at res=40, and this repo has hit real
  mode-identity/ordering issues across resolutions/solvers before
  (`fdtd_vs_EME_mode_validation.md`, the modal multiplexer's own gaps).
  **Checked at res=25 (decay_tol still 1e-9, same 3-band comparison):**

  | band | res=15 decay (t) | res=25 decay (t) |
  |---|---|---|
  | 1 | 108.0 | 99.0 |
  | 2 | 114.7 | 130.7 |
  | 3 | 150.7 | 93.7 |

  All three still decay in a modest, similar range — no sign of building
  toward the catastrophic res=40 stall. The near-degenerate normalization
  mismatch actually SHRANK from ~36% (res=15) to ~15.6% (res=25) — moving
  in the healthy direction, not the unhealthy one. **This makes the
  res=40 failure look like a sharp threshold, not a gradual slowdown** —
  shifts the balance toward candidate (b), mode reordering/crossing
  somewhere between res=25 and res=40 (or a resolution-specific numerical
  resonance coincident with res≈40 for this cross-section), over a smooth
  numerical degradation.

  **✅ Mode-identity check done — candidate (b), reordering, is RULED
  OUT.** Cheap `section_mode_neff()` calls (no time-stepping, ~seconds,
  short window at the body's own midpoint per that method's own PML
  convention) at res=15/25/40 for bands 0-5:

  | band | res=15 | res=25 | res=40 |
  |---|---|---|---|
  | 0 | 2.893802 | 2.940424 | 2.962350 |
  | 1 | 2.887617 | 2.934357 | 2.956336 |
  | 2 | 2.877280 | 2.924218 | 2.946286 |
  | 3 | 2.862749 | 2.909967 | 2.932160 |
  | 4 | 2.843960 | 2.891543 | 2.913899 |
  | 5 | 2.820833 | 2.868868 | 2.891428 |

  Every band increases smoothly and monotonically with resolution
  (expected discretization convergence), the ORDERING never changes at
  any resolution, and the band-0-to-1 gap stays essentially constant
  (~0.006) across all three — no resolution-dependent near-degeneracy
  emerging either. **Band 1 at res=40 is unambiguously the same physical
  mode as band 1 at res=15/25.** This rules out mode-crossing/reordering
  as the explanation, leaving only candidate (a): a genuine,
  resolution-specific numerical effect in the FDTD time-stepping/
  decay-detection itself (not which mode is targeted) — e.g. a
  coincidental numerical resonance specific to res≈40 for this exact
  domain.

  **New shared-library option added 2026-08-18 to test this directly:
  `FDTDSolver(decay_method=...)`.** Investigated which meep early-stop
  condition is actually in use — confirmed directly against meep's own
  API (`dir(meep)` lists six: `stop_when_dft_decayed`,
  `stop_when_fields_decayed`, `stop_when_energy_decayed`,
  `stop_after_walltime`, `stop_at_direction`, `stop_on_interrupt`) — that
  **every** `sim.run()` call in `beef/solvers/fdtd.py` uses only
  `stop_when_dft_decayed(tol=decay_tol, minimum_run_time=minimum_run_time)`,
  with `maximum_run_time` never set. Two properties of that criterion
  matter here: it requires **every** field component and frequency on
  **every** registered DFT monitor to individually decay (confirmed via
  the code's own comment: "watches the whole field, not any one
  monitor"), and it has **no built-in ceiling** — a single anomalously
  slow-decaying component has nothing to stop it early, consistent with
  the observed stall (no error, just running with no exit).

  Added `decay_method: str = "dft"` (default, byte-identical behavior —
  every existing device/golden unaffected) and `energy_decay_dt: float =
  20.0` to `FDTDSolver.__init__`, plus a `_stop_condition()` helper every
  `sim.run()` call site now goes through. `decay_method="energy"` switches
  to `mp.stop_when_energy_decayed(dt=energy_decay_dt, decay_by=decay_tol)`
  — one bulk-energy-over-the-whole-cell check instead of a per-component
  requirement. Targeted test suite (`-k "fdtd or beta"`, 117 passed/2
  skipped) confirms no regression from the new default.

  **Sanity check (res=15, decay_tol=1e-9, mode 1): energy method gives
  t=100.2 (3005 steps) vs. the dft method's own t=108.0 (3241 steps) for
  the identical setup** — consistent, sensible agreement, confirms the
  new option itself is implemented correctly before trusting it at the
  resolution that actually matters.

  **✅ Real test (res=40, decay_tol=1e-11, mode 1) — CONVERGED, decisively.**
  The exact configuration that stalled past 660,000+ timesteps under
  `stop_when_dft_decayed` (never converging) reached decay at **t=140.1
  (11,207 timesteps)** under `decay_method="energy"` — a completely
  normal decay time, in line with the lower-resolution results (res=15:
  t=108.0, res=25: t=99.0), not an outlier. **This directly confirms the
  diagnosis**: the DFT criterion's requirement that *every* individual
  field component/frequency on *every* monitor decay was getting stuck on
  something specific to res=40, while the bulk energy -- the physically
  meaningful quantity -- was decaying completely normally the whole time.
  `stop_when_energy_decayed` sidesteps that per-component requirement
  entirely.

  **Not yet done, before trusting this for a real reported calibration
  number**: verify the beta/phase extracted under the energy criterion is
  actually as accurate as a fully-DFT-converged run would give -- stopping
  on bulk energy rather than per-component DFT decay is a different
  (looser) criterion, and this single-mode test only confirms it
  *terminates*, not that the extracted phase is fully settled. The
  natural next check: run all 16 modes with `decay_method="energy"` at
  production settings and see whether the result changes the CMRR/
  phase-error picture in a way consistent with everything already known
  about this design point, or compare mode 1's own extracted beta against
  a loosened-decay_tol DFT run (if one can be found that terminates) as a
  cross-check.

  **Full 16-mode production calibration submitted, 2026-08-18.** Script
  updated (`mmi_4x4_90hybrid_calibrated_baseline.py`): `make_solvers()`
  now returns a SEPARATE `calibration_solver` (`decay_method="energy"`)
  distinct from `input_solver` (left at the "dft" default) — setting
  `decay_method="energy"` on `input_solver` directly would also change
  the fingerprint for the already-cached (job 11778128) baseline
  input-matrix fetches, turning free cache hits back into wasted
  re-solves, the exact bug hit twice already. Verified directly before
  submitting: default-`"dft"` solvers still cache-hit all 4 baseline
  input matrices (0.17s), and an `"energy"` solver produces a genuinely
  different cache config (`beef/cache/config.py` updated the same way
  `minimum_run_time` was — `decay_method`/`energy_decay_dt` added to
  `_METHOD_ATTRS`'s three FDTD allowlists AND to `_DEFAULT_OMITTED_ATTRS`,
  so a default-valued run's fingerprint is untouched but a non-default
  one is never silently served by/serves a `"dft"` entry). Targeted
  suite (`-k "fdtd or beta or cache or config"`) re-run clean: 237
  passed, 4 skipped. **Submitted as job 12052383 → hit
  `AssocGrpBillingMinutes` again → cancelled, resubmitted under
  `gts-sr108` as job 12052385** (now `ReqNodeNotAvail`, a normal resource
  wait, not a block).

  **❌ TIMED OUT AGAIN at the full 6h — but this time a real, different,
  and much more encouraging failure.** Modes 1-15 of 16 completed
  CLEANLY under `decay_method="energy"`, with sensible, gently-climbing
  decay times (160.1 -> 260.2 meep-time-units across the 15 completions)
  — the energy fix is confirmed working broadly, not a fluke on mode 1
  alone. **Mode 16 (the last, weakest-confined/highest-order of the 16
  requested — closest to the cladding cutoff) is a NEW stuck case**: its
  own energy-decay ratio sat at exactly `1.0` at two separate checkpoints
  (t=1260.8 and t=1280.8, log: `energy decay(t=...): X/X = 1.0`) — no
  decay progress at all, not even slow. Ran to its own t≈1280+ (5-8× the
  other modes' own ~160-260 range) before the wall hit. **All 15
  completed modes' work is lost again** — `solve_section_multi_length()`
  only writes to cache after the FULL 16-mode loop finishes, no partial
  credit. Physically plausible explanation, not yet confirmed: the
  weakest-guided mode near cutoff genuinely radiates/couples to the
  continuum differently than the other 15, and may need its own
  treatment (a real `maximum_run_time`-style ceiling — which
  `stop_when_energy_decayed()` does not natively support, unlike
  `stop_when_dft_decayed()` — or a custom combined stop condition) rather
  than assuming energy-decay alone fixes every mode uniformly.

  **Checked directly, not assumed: mode 16 is NOT near cutoff.** Real
  n_eff/margin data for all 16 modes at res=40: mode 16's margin above
  n_clad is 0.883 (vs. mode 1's 1.516) — the weakest-confined of the 16
  requested, but still solidly guided, not a marginal/near-cutoff case.
  This rules out "near-cutoff weak guiding" as the explanation.

  **Testing mode 16 in TRUE isolation** (a runtime-patched copy of
  `solve_section_multi_length()` exciting ONLY band 16, first and only
  simulation in a fresh process — `mmi_4x4_90hybrid_mode16_isolated.py`)
  to distinguish an intrinsic mode-16 property from a "16th-in-a-long-
  sequence" resource-accumulation artifact. First attempt on the login
  node showed anomalous per-step cost (0.3-1.0 s/step, 10-30× every other
  run at this resolution) — traced to genuine contention on a shared
  login node (load avg ~10-13, 37 other users), not the simulation
  itself; killed and moved to batch. **Submitted as job 12066640 → hit
  `AssocGrpBillingMinutes` → resubmitted under `gts-sr108` as job
  12066668** (`-t 240`, normal queue wait).

  **✅ DECISIVE RESULT — mode 16 converges completely normally in
  isolation, ruling out an intrinsic mode-16 problem.** Job 12066668
  COMPLETED in 23m28s: mode 16 alone decayed at t=220.1 (17,611
  timesteps), right in line with modes 1-15's own 160-260 range. The
  energy-decay log shows genuine, clean, monotonic decay — ratio=1.0 only
  at the earliest two checkpoints (t=20, t=40, the normal pulse-buildup
  transient, not a stall), then dropping smoothly (8.3e-4 → 3.6e-7 → ...
  → 7.5e-12) well past `decay_tol`. This is qualitatively different from
  the failing production run, where the ratio stayed stuck at exactly
  `1.0` at LATE times (t=1260.8, t=1280.8) with no transient excuse.

  **Conclusion: mode 16 has no intrinsic problem — the bug is specific to
  running it as the 16th simulation in one long sequential loop within a
  single process**, some form of resource accumulation across 15 prior
  `mp.Simulation` builds/`reset_meep()` calls, not a property of the mode
  itself. **Dropping mode 16 (N_MODES=15) is NOT the right fix** — it
  would lose real accuracy (the extended body convergence check already
  showed the basis isn't fully flat even at N=16-36) to work around a
  sequence-position bug, not a genuine mode-16 limitation. The better fix
  is isolating each excited mode's own simulation (e.g. one process/job
  per mode, matching how this isolation test itself was structured).

  **Built and submitted, 2026-08-20.** Two new scripts implement the
  one-process-per-mode architecture:
  - `mmi_4x4_90hybrid_calibration_mode_worker.py` — takes a mode index
    (1-16) via CLI arg, computes the SAME seed/lengths `section_mode_beta_
    details()` would (unpatched — `section_mode_neff()` is cheap/cached,
    not the stalling part), then calls a runtime-patched
    `solve_section_multi_length` that excites ONLY that one band (same
    technique as the mode-16 isolation test), and saves that mode's own
    phase-vs-length data + seed β to a per-mode `.npz` file.
  - `mmi_4x4_90hybrid_calibration_aggregate.py` — combines all 16
    per-mode files by calling the SAME imported `beta_fit_from_phases()`
    and the SAME guard thresholds (`_BETA_FIT_RESIDUAL_MAX`,
    `_BETA_BRANCH_MARGIN_RAISE`) the library itself uses, so a mode that
    would fail those guards in a single-process run fails them here too.

  **Correctness, not just architecture, checked directly against the
  library's own code**: `section_mode_beta_details()`'s post-fit step
  calls `beta_fit_from_phases(lengths, phases[:,i,m], seed_beta[i,m])`
  once per (wavelength, mode) cell, using ONLY that mode's own data — no
  cross-mode mixing. Splitting modes across processes and merging
  afterward is therefore mathematically identical to a single-process
  run, not an approximation.

  **Verified before spending any real cluster time** (all free,
  cache-only): the worker's own seed computation hits job 12052385's
  already-cached `section_mode_neff()` result exactly (`cache_only=True`
  succeeded), and its derived `lengths` (`[0.842, 1.683, 2.525]` µm, using
  all 3 production wavelengths' own β_max — slightly different from the
  earlier single-λ diagnostic tests' `[0.885, 1.77, 2.655]` for exactly
  the expected reason) match the library's own formula exactly.

  **Submitted as a 16-task SLURM array, job 12069540** (`-t 90`/1.5h per
  task, real margin over the ~20-25 min the mode-16 isolation test itself
  took), all 16 modes running in PARALLEL rather than the previous
  sequential 6h-timeout approach. Results pending.
  **Tasks 2-16 hit `AssocGrpBillingMinutes` on `gts-sr108-ece`** (task 1
  started before the account's balance was exhausted and kept running);
  cancelled just the pending 15 (`scancel 12069540_[2-16]`, left task 1
  alone) and resubmitted them as **job 12069773** (`--account=gts-sr108
  --array=2-16`, same batch file). All 16 tasks (`12069540_1` +
  `12069773_[2-16]`) reached RUNNING by 2026-08-20. Results pending —
  run `mmi_4x4_90hybrid_calibration_aggregate.py` once all 16
  `mode_XX.npz` files exist.

  **Cross-device relevance checked, 2026-08-20**: is either the
  `decay_method="energy"` fix or the sequence-position stall relevant to
  any OTHER device's own calibration? Every device in this repo currently
  running `calibrate_body_neff=True` was checked for its own body
  `N_MODES` — MMI 2×2 (14, the highest), 1×4 (12), 1×2 (8), directional
  coupler (5, SOI160 variant 2), HCG Fabry-Perot's coupler body (2) + arm
  (1), ring resonator's coupled-pair body (2). None are close to 16, and
  none have been observed to stall — but the trigger condition for the
  sequence-position bug is confirmed only at "16th of 16 on this specific
  device," not tied to a known numeric threshold, so this isn't proof
  those devices are safe by construction, just that nothing needs to
  change on them today. Full note, with the file:line table, added to
  `resolved/section_mode_beta_multi_monitor_validation.md` (the central
  cross-device calibration doc) and to the
  [[beta-calibration-cost-model-step9]] memory, so future work on any
  calibrating device sees this before assuming a bigger `-t` fixes an
  unexpectedly long calibration run.

  **Correction, same day — running the full production array reopens this,
  it isn't closed.** Modes 1-15 (jobs `12069540_1` + `12069773_[2-15]`) all
  completed cleanly. **Mode 16 (`12069773_16`) did NOT converge in
  isolation either — it hit `TIMEOUT` at the full 90-min budget**, and its
  own log shows why: not a stall, **exponential energy growth** (ratio
  bottoms at 1.95e-11 near t≈260 — just above `decay_tol=1e-11`, so the
  criterion correctly never fires — then climbs monotonically to 1.7e53 by
  t≈2081, still climbing at timeout). Re-reading job 12066668's own energy
  trace with this in mind: it was *also* still falling with no turnaround
  visible through its last checkpoint, and crossed below `decay_tol` only
  because its domain (`max(lengths)=2.655µm`, single-λ diagnostic lengths)
  was slightly longer than production's (`2.525µm`) — a near-miss, not
  evidence of stability. **The "mode 16 has no intrinsic problem, only a
  sequence-position bug" conclusion is retracted** — isolating the mode
  changes nothing about whether it converges; the instability is intrinsic
  to mode 16's own configuration.

  **Root cause, from the repo's own code**: `beef/materials/stack.py`
  builds Si's dispersion as `mp.LorentzianSusceptibility(..., gamma=0,
  ...)` — an undamped pole (same convention as meep's own `materials.Si`),
  a known marginal-stability risk over long FDTD integration times.
  Calibration's ~2000+ time-unit runs at `decay_tol=1e-11` are far more
  demanding than any ordinary S-parameter solve in this repo, which is
  plausibly why this has never surfaced before. Full account, including why
  this also changes the "which other devices are at risk" answer (mode
  count was never the right proxy — integration time and decay_tol are),
  in `resolved/section_mode_beta_multi_monitor_validation.md`'s own
  same-day correction section. **Not pursued further without direction —
  a materially different, more open question than what was believed a few
  hours earlier today.**

  **Decision, 2026-08-20 (user's call, option 3 of 3 offered): drop mode 16,
  proceed on N_MODES=15.** The alternatives offered and not taken were (1) a
  small nonzero `gamma` on the dispersion model to damp the unstable pole,
  and (2) retrying mode 16 under `decay_method="dft"`. Chosen over both
  because the instability is a real numerical property of this
  configuration, not a stopping-criterion artifact — neither alternative is
  known to fix it, and both cost more cluster time to find out.

  **Aggregation of the 15 modes that exist PASSES every library guard**
  (`MMI4X4_CAL_MODES=15`, free, login-node): worst fit residual 5.617e-03
  rad (guard `_BETA_FIT_RESIDUAL_MAX`=0.1), worst branch margin 0.161
  (guard `_BETA_BRANCH_MARGIN_RAISE`=0.9), worst seed excess +5.65% at
  λ=1.36 µm mode 15 (warn threshold 10%). That last number is worth
  noting: the earlier reduced-settings (res=15) test flagged **+13.69%**
  at this same highest mode, and at real production settings it comes in
  at less than half that — so that anomaly was substantially a
  reduced-settings artifact, though +5.65% is still above the ~2%/mode
  discretization gap this calibration exists to measure. Calibrated
  n_eff at λ=1.26 µm: 3.0403 (mode 1) monotonically down to 2.5982
  (mode 15). Saved as `aggregated_beta_n15.npz`.

  **The 4 input FDTD matrices are NOT re-solved at 15 modes.** They are
  fetched from job 11778128's cache at 16 modes (verified free before
  submitting: all 4 HIT in 1.42 s, shape `(3, 16, 1)`) and **truncated** to
  the first 15 (`M[:, :15, :]`). This rests on the same per-band
  independence the repo already relies on for `measure_diagonal_only` and
  for the per-mode calibration split — meep decomposes each band
  separately, so bands 1-15 do not depend on whether band 16 was also
  requested. Re-solving at `n_modes_out=15` would be a genuine cache miss
  costing ~10 h of FDTD for a result expected to be bit-comparable; not
  spent. **Noted as a deliberate assumption, not a measured equivalence** —
  it has not been checked directly at production scale here.

  **Submitted as job 12078959** (`-t 360`,
  `run_mmi_4x4_90hybrid_calibrated_baseline_n15.batch`,
  `mmi_4x4_90hybrid_calibrated_baseline_n15.py`), running two fresh
  81-length EME body sweeps at 15 modes — uncalibrated AND calibrated —
  **deliberately both**, because an N=15 calibrated result differs from the
  published baseline in two ways at once:
  - **[A] calibrated-15 vs uncalibrated-15** — identical basis both sides,
    so this isolates the *calibration* effect. This is the number that
    answers the original question this whole thread exists to answer.
  - **[B] calibrated-15 vs the published uncalibrated-16** (job 11778128:
    max|Δφ|=12.3°) — confounds calibration with the dropped mode, reported
    as context only, NOT as the calibration effect. The script also prints
    **[B'] uncal-16 → uncal-15**, isolating the dropped mode's own cost, so
    the two effects can be read apart rather than guessed at.

  **✅ RESULT, job 12078959 COMPLETED — and it ANSWERS this whole thread's
  motivating question: no, calibration does NOT rescue this device.**

  | metric | uncal n=16 (published) | uncal n=15 | **cal n=15** | requirement |
  |---|---|---|---|---|
  | max\|Δφ\| | 12.320° | 12.322° | **12.472°** | ≤ 5° |
  | min CMRR | 1.531 dBe | 1.532 dBe | **4.078 dBe** | ≥ 20 dBe |
  | max excess loss | — | 8.366 dB | **7.556 dB** | (not a cited bar) |
  | max imbalance | — | 14.615 dB | **9.053 dB** | (not a cited bar) |

  **[A] the calibration effect, isolated** (identical 15-mode basis both
  sides): max|Δφ| **12.322 → 12.472°** (+0.150, slightly *worse*), min CMRR
  **1.532 → 4.078 dBe** (+2.547), excess loss −0.810 dB, imbalance
  **−5.562 dB** (the largest single move). So the FDTD-consistent β is a
  real, measurable effect — it substantially improves imbalance and lifts
  the CMRR floor — but it leaves the headline phase error **unchanged**, and
  **both cited requirements still fail by a wide margin** (CMRR 4.1 vs 20
  dBe needed; |Δφ| 12.5 vs 5° needed).

  **This closes the hypothesis that motivated the entire calibration
  thread.** The premise (script docstring, and open-decisions item 5) was
  that because Δφ_j is driven by Δβ over the body length, an uncalibrated
  ~2.3%/mode β error could be inflating the phase error — i.e. that this
  device's failure might be substantially a *modeling artifact* rather than
  a real limitation. Measured: it is not. Phase error moves by +0.15° when
  the β's are corrected, against a 7.3° shortfall. **The device's CMRR/phase
  failure is real physics of this design point, not a calibration artifact.**

  **[B'] dropping mode 16 costs essentially nothing**: uncal 16 → uncal 15
  moves max|Δφ| by **+0.002°** and min CMRR by +0.0002 dBe. This
  retroactively validates the decision to drop it — the mode that could not
  be simulated turns out to be the one that barely matters. (Consistent
  with the earlier N_MODES_TAPER finding: fundamental-mode launch couples
  very weakly into the highest-order modes of this geometry.)

  **Calibration does not shift the optimum length either** — scanning the
  full 81-point ΔL grid on the calibrated body returns **ΔL = −2.75**, the
  same optimum as uncalibrated, at the same metrics. So there is no
  "re-tune ΔL for the calibrated β" recovery available.

  **Cost note, and a real mechanism worth knowing**: the job finished in
  **22 s**, not the ~3 h budgeted for two fresh 81-length sweeps. The
  calibrated sweep was genuinely computed fresh (its own new cache entry,
  `neff_override` in the key) but takes `neff_override`'s documented
  diagonal fast path; the *uncalibrated* n=15 sweep was served by
  `lookup_superset()` (`beef/cache/store.py:348`) — a **designed** cache
  tier that answers an n-mode request from a stored ≥n-mode result by
  slicing its trailing axes down. Verified this was that tier and not a
  wrong-key hit: returned shape is `(81, 3, 15, 15)`, and the only n=15
  `length_sweep` entry on disk is the calibrated one written by this job.
  **This makes comparison [A] cleaner than the script was designed for** —
  both sides share the *same* 16-mode-solved mode shapes sliced to 15, so
  the only thing differing between them is the propagation constants,
  which is exactly the quantity being tested.
- **2026-08-14, Step 6b (full-device FDTD reference) RE-QUEUED at the
  user's direction, job 11847078** — see Step 6b's own updated section
  below for the full account (stale-geometry fix + budget raised to 3
  days given job 11733865's inconclusive 6h timeout). Results pending.
- **2026-08-13, queued as the first concrete lever on that open question:
  a coarse `waveguide_offset` sweep** (open decisions item 5 below) --
  `waveguide_offset` is currently *derived* (`=W/4`, exactly reproduces
  paper Eq. 8) rather than independently validated for this platform, and
  directly controls where each access waveguide couples into the body's
  self-imaging pattern -- a plausible lever for the CMRR/phase gap. Made
  cheap by two findings, both checked, not assumed (see "Top-bottom
  mirror-symmetry solve reuse" and "waveguide_offset body-cache
  independence" sections above/below): the body mode solve doesn't need
  re-running per offset point (confirmed via `cache_only=True` against
  job 11768847's cache), and each offset point only needs 2 new FDTD
  input solves, not 4, via the same mirror-reuse trick as
  `mmi_2x2_hybrid.py`'s own `MODE_SIGNS`. **Script:
  `examples/verification/python/mmi_4x4_90hybrid_offset_sweep.py`,
  smoke-tested end to end at reduced settings first** (caught and fixed one
  real bug: the body-mode cache lookup also needs a production-vs-test
  fallback, same class of issue as the sign-vector derivation). **Submitted
  2026-08-13, job 11817429** (4-task array, `run_mmi_4x4_90hybrid_offset_
  sweep.batch`, `gts-sr108-ece`, offsets {1.2, 1.35, 1.65, 1.8} µm — the
  baseline `offset=1.5` needs no new run, it is job 11778128's own cached
  result, already fetched and confirmed to reproduce those exact numbers).
  **✅ RESULTS LANDED 2026-08-13, job 11817429 (all 4 tasks COMPLETED,
  2:10h–3:02h each) — the offset lever makes things WORSE, not better, at
  every point tested, and the derivation mechanism itself was independently
  confirmed correct, not just theoretically argued:**

  | offset (µm) | max\|Δφ\| (deg) | min CMRR (dBe) | max excess loss (dB) | max imbalance (dB) |
  |---|---|---|---|---|
  | 1.20 | 157.348 | 0.713 | 5.553 | 22.293 |
  | 1.35 | 107.460 | 0.044 | 6.901 | 35.416 |
  | **1.50 (baseline, job 11778128)** | **12.320** | **1.531** | **8.366** | **14.615** |
  | 1.65 | 51.365 | 2.908 | 5.668 | 10.674 |
  | 1.80 | 162.849 | 1.135 | 5.101 | 15.620 |

  Every offset away from the baseline is DRAMATICALLY worse on phase error
  (51°–163° vs. 12°) — not a smooth, mild degradation. Before trusting
  numbers this surprising (they echo A11's own old phase-pinning bug
  superficially), ran a direct, independent confirmation rather than
  accepting `check_passive()`'s clean pass as sufficient (it only bounds
  power, not phase — a sign-derivation bug on the phase-sensitive derived
  ports would pass it too): at offset=1.35, solved ports 2 and 3 for REAL
  (not derived) and compared against the sign-vector-derived values used
  in the sweep — **agreement to ≤0.84° angle, ≤0.2% magnitude**, at a
  genuinely new offset never used to derive the sign vector. The reuse
  mechanism is confirmed correct away from the baseline, not just at it;
  the phase blowup is real physics, not a screening artifact.

  **Interpretation: `mm_length`/`ΔL=-2.75` was tuned specifically FOR
  `offset=1.5` (Step 3's own ΔL sweep, run only at that offset) — moving
  the ports without re-tuning the body length breaks the self-imaging
  phase alignment badly, consistent with MMI self-imaging being a discrete,
  position-sensitive interference pattern, not a smoothly-varying one.**
  `waveguide_offset` is NOT a free, independent lever on its own; it is
  coupled to `ΔL` — **checked directly the same day, not left as a
  hypothesis, see "Joint (offset, ΔL) grid" below: re-optimizing ΔL per
  offset does not rescue it either.** offset is now a closed, ruled-out
  lever.

### Joint `(offset, ΔL)` grid — 2026-08-13, zero new cluster time

The single-ΔL screen above only ever evaluated each offset against the ONE
`ΔL=-2.75` slice tuned for `offset=1.5`. Two things needed for a real
answer were already sitting on disk, needing no new solves: job
11817429's own 2 real + 2 sign-derived input matrices per new offset (job
11778128's 4 real for the baseline), and the FULL 81-point `ΔL` body
length-sweep grid (job 11768847) — confirmed offset-independent (the
"`waveguide_offset` body-cache independence" finding above), so the SAME
grid applies at every offset, not just 1.5. `mmi_4x4_90hybrid_offset_
deltaL_grid.py` re-combines the already-cached input matrices against
every one of the 81 `ΔL` points per offset (pure numpy, no meep/meow call)
— 21 s on the login node for the full 5×81 = 405-point grid.

**Result: even each offset's own best achievable `ΔL` stays far worse than
the baseline's best point.**

| offset (µm) | best ΔL in [-10,+10] µm | max\|Δφ\| at that ΔL (deg) | min CMRR there (dBe) |
|---|---|---|---|
| 1.20 | +7.50 | 125.563 | 1.790 |
| 1.35 | +0.25 | 64.802 | 2.723 |
| **1.50 (baseline)** | **−2.75** | **12.320** | **1.531** |
| 1.65 | −1.00 | 47.638 | 5.170 |
| 1.80 | +6.50 | 75.121 | 0.420 |

The global best point over the entire explored 5×81 grid remains
`offset=1.5, ΔL=-2.75` — the paper's own Eq. 8 port-spacing formula
(`offset=mm_width/4`) is not merely "whatever `ΔL` happened to be tuned
for," it is genuinely close to optimal across this range, independent of
`ΔL`. `waveguide_offset` is a closed, ruled-out lever now, not an open
one — the next candidate levers (per the 2026-08-09 open question) are
`mm_width`, or the still-unmeasured `taper_width`/`taper_length`
placeholders (Step 2's own flagged, never-independently-optimized
defaults), both of which would need genuinely new FDTD/EME solves, not
free reuse of what is already cached.

---

**Status 2026-08-05:** Step 1's cross-section sweep (1a single-mode width,
1b sz convergence, 1c L_π vs. body width, all at O-band on `SOI_220`) is
running as SLURM job **11692410** (`batch/verification/
run_mmi_4x4_90hybrid_cross_section.batch`, account `gts-sr108-ece`, 8 cpus,
`-t 240`) — moved to batch after a login-node attempt at
`OMP_NUM_THREADS=1` was still running past 12 minutes on a single 2-mode
point with no output yet; `run_mmi_2x2_gap2_dbeta.batch`'s own comment
records 206 s for one comparable point on that device's (smaller) body
cross-section, so this is genuine batch-scale work, not a hung process.
While that runs, **Step 7's math primitives landed** (see A8) and **A5's
mirror-symmetry claim is now confirmed directly**, not just argued from
source (see A5) — both meep-free, done on the login node.

**Update, same day — 11692410 cancelled, resubmitted corrected as
11693043 (see A10 below).** Partial output before cancellation: even the
sweep's *narrowest* width (0.30 µm) already showed 2 guided modes at O-band
on `SOI_220`, not the single mode a 220 nm Si strip normally supports —
enough of a surprise to investigate before trusting the range, and enough
to trigger a real, recorded dead end (A10): trying to separate the count by
polarization (the paper excites only a TE-polarized fundamental mode) using
meow's `Mode.te_fraction`, then finding that metric disagrees with this
repo's own already-shipped precedent. Reverted; the width range is now
widened downward (0.10–0.70 µm, was 0.30–0.70) and the job resubmitted on
the plain, repo-standard guided-mode count.

**Update, same day — 11693043 COMPLETED (3:19:14). Step 1 is closed.**
`waveguide_width` is DECIDED: 0.35 µm, the user's own choice, not derived
from the search. 1c's `L_π/W²` came out essentially constant (3.026 →
3.020 µm⁻¹ over 6–16 µm) — textbook self-imaging scaling, a real
cross-check that these are genuine body modes. Full results, the
2×2-precedent parallel that reframed Step 1a's own success criterion, and
the caveat on 1b's box-size convergence (measured at 0.45 µm, not yet
re-checked at the decided 0.35 µm) are all in the "Step 1 MEASURED" box
right after A10 below.

**Step 2 is now also done.** `W=6.0 µm` (the most compact measured point),
`waveguide_offset=W/4=1.5 µm` (checked directly to reproduce paper Eq. 8
exactly via `beef`'s own default port generator — no custom port lists
needed), `L=(3/4)L_π=81.69 µm` — **5.6× shorter than the paper's own
459 µm silica body**, which is the compactness A9 predicted, now measured
rather than just expected. `taper_width`/`taper_length` are explicit,
flagged placeholders (`MMI`'s own class defaults), not yet measured. Full
parameter table and the geometry-overlap check (0.7 µm clearance, checked
not assumed) are in Step 2's own box.

**Step 3 is done (first pass) as of 2026-08-06.** Real production sweep
(job 11708624, 48 ranks, 3:54h) hit a real bug on first read — max|Δφ|
pinned near 180° across the entire ΔL sweep while imbalance behaved
normally — diagnosed with ZERO new cluster time (every solve already
cached) as a global S-parameter phase-sign-convention mismatch between
meep/meow and the paper's own convention, not a device or code indexing
bug (A8's own predicted cross-call phase-reference risk was checked
directly and ruled out first). Fixed and reverified: **best point
ΔL=-1.50 µm, max|Δφ|=8.28°, max imbalance=4.90 dB** — close to but not
quite meeting the cited 5°/CMRR≥20dBe requirements, on a first pass with
only 8 body modes and unmeasured placeholder tapers. Full diagnostic
trail and corrected sweep table in **A11** below. **Step 4 (geometry with
real numbers) is next**, though refining N_MODES/taper dimensions first
may close the remaining phase gap — an open call, see A11's own closing
paragraph.

> ### A10. `Mode.te_fraction` does not match this repo's own "TE" convention — a dead end, recorded not hidden
>
> **The observation.** A first pass at 1a filtered guided modes by
> `te_fraction > 0.5` (0.0 = pure TM, 1.0 = pure TE, per meow's own
> docstring), reasoning that the paper's device excites only a TE-polarized
> fundamental mode (Sec. IV) and that these narrow, near-square access
> waveguides (width ≈ height ≈ 0.2–0.4 µm) are exactly the aspect ratio
> where a TE0/TM0 pair sits close in n_eff.
>
> **The sanity check, and why it matters.** Before trusting a brand-new
> filter, it was checked against a cross-section this repo already trusts:
> `mmi_2x2_hybrid.py`'s own access waveguide (`width=0.45`, C-band, SOI_220
> — documented throughout this repo, including that script's own module
> docstring reasoning, as running with a TE-polarized fundamental mode).
> Its dominant mode (`n_eff=2.311740`, exactly the index this repo's own
> devices already report and depend on) came back **`te_fraction=0.000199`
> — i.e. meow calls this repo's own working "TE" mode "TM."**
>
> **What this does and doesn't establish.** It does NOT prove meow's
> `te_fraction` is buggy — it may be using an internal axis/field-component
> convention that doesn't line up with this codebase's informal "TE"
> labelling, or the real vector-mode character may be more subtly
> hybridized than the paper's own 2-D scalar TE/TM split (Eq. 2) assumes —
> that split is an effective-index-method reduction this codebase's full
> 3-D vector solve (meow/meep) does not make, so the paper's clean
> TE-vs-TM decoupling does not automatically transfer. What it DOES
> establish: **this metric cannot currently be trusted to separate TE from
> TM in this codebase**, because it disagrees with a case this whole repo
> already treats as settled.
>
> **Why this is not blocking, and why no device has ever needed it.** Every
> existing device in this repo (2×2, 1×2, 1×4, the coupler, the modal
> multiplexer) counts guided modes as plain `n_eff > n_clad` bands, with
> zero polarization filtering, and none of their validation docs raise TE/TM
> conflation as an open risk. The reason is structural, not incidental: a
> plain rib/strip cross-section with no lateral or vertical asymmetry keeps
> TE and TM as exact, non-interacting symmetry classes through any
> symmetric perturbation (taper, MMI body, …) — so an extra wrong-
> polarization band sitting unused in an interface basis costs a wasted
> mode slot, never a real coupling error. This device's access waveguides
> and MMI body have exactly that same symmetry, so the same reasoning
> applies without needing a runtime filter.
>
> **Decision: reverted to the plain, repo-standard guided-mode count.** Not
> because the TE/TM question is unimportant, but because (a) the tool tried
> for answering it does not currently agree with this repo's own ground
> truth, and (b) the repo's own structural argument already covers why it
> doesn't matter for this class of device. If a genuinely asymmetric
> section is ever added to this device (it is not planned to be), this
> reasoning would need revisiting.

### ✅ Step 1 MEASURED 2026-08-05 — job 11693043, and `waveguide_width` DECIDED by the user (0.35 µm), not derived from the search

**1a (single-mode width, plain guided-mode count, no polarization filter):
empty across the whole swept range.** Every width from 0.10 to 0.70 µm shows
≥2 guided modes at O-band on `SOI_220` — even 0.10 µm (second band margin
only 0.0007 above cladding). The trend is smooth and monotonically
increasing with width (not the flat, width-independent signature that
flagged a real `compute_modes_meep` bug in the SiN 1×4 saga — see that doc's
own step 4a/4b), so this is read as genuine physics, not a solver artifact.

| width (µm) | guided | n_eff (highest → lowest) |
|---|---|---|
| 0.10 | 2 | 1.474411, 1.447551 |
| 0.15 | 2 | 1.543006, 1.471209 |
| 0.20 | 2 | 1.648716, 1.600483 |
| 0.25 | 2 | 1.897532, 1.760734 |
| 0.30 | 2 | 2.185559, 1.859633 |
| 0.35 | 3 | 2.385726, 1.940421, 1.457579 |
| 0.40 | 3 | 2.520250, 2.004680, 1.512177 |
| 0.45 | 4 | 2.613334, 2.055597, 1.615018, 1.459231 |
| 0.50 – 0.70 | 4 | (climbing further) |

**This sharpened, not just failed to close, A10's structural argument.**
Re-checking `mmi_2x2_hybrid.py`'s own shipped access waveguide (0.45 µm,
C-band — used for A10's `te_fraction` sanity check) shows the *same*
signature: 2 clearly-guided modes (n_eff 2.311740, 1.696198), not 1. **This
repo's own already-validated, in-production SOI_220 device does not meet
the "1 guided mode of any polarization" bar either** — it works because of
A10's structural argument (exact symmetry decouples TE/TM regardless of
what else is in the basis; the FDTD source explicitly launches one target
eigenmode), not because its cross-section is literally single-band-guided.
So Step 1a's own original success criterion was stricter than this repo's
own working precedent actually requires — the goalpost, not just the
search range, needed correcting.

**Put to the user directly rather than decided unilaterally**, given the
direct parallel to the SiN 1×4 plan's own "reconsider the design width —
user's call" moment. **Decision: `waveguide_width = 0.35 µm`** — the user's
own choice, not derived from the mode-count search above (0.35 is not the
narrowest 2-guided-mode point measured — 0.30 is — and is not framed as a
deliberate N-mode design either; it is stated as a direct design value).
At 0.35 µm the real measured cross-section is **3 guided bands**: n_eff
2.385726 (fundamental), 1.940421 (second), 1.457579 (third, weakly guided,
margin 0.011 above cladding). Carried forward into Step 2 as-is.

**1b (box convergence at 0.45 µm — the trial width `main()` picked before
the user's decision above; not re-run at 0.35, see caveat below):**
converges cleanly, settling to 6 decimal places by `sz=2.0` — already
`SOI_220`'s shipped default, no change needed on this axis.

| sz (µm) | n_eff[0] | step |
|---|---|---|
| 1.0 | 2.613445 | — |
| 1.5 | 2.613335 | −0.000110 |
| 2.0 | 2.613334 | −0.000000 |
| 2.5 | 2.613334 | −0.000000 |
| 3.0 | 2.613334 | −0.000000 |

⚠ **Caveat, explicitly carried, not silently papered over:** this
convergence was measured at 0.45 µm (the script's own `main()` picks the
widest single-mode width from 1a, which was empty, so it fell back to its
0.45 µm default) — *before* the user fixed `waveguide_width=0.35`. Given how
fast and clean the 0.45 µm convergence was, and that 0.35 µm is narrower
(a smaller, more tightly bound mode is not expected to need a *bigger* box),
`sz=2.0` is very likely fine at 0.35 µm too — but this is an argued
expectation, not a re-measurement at the real, decided width. Cheap to
close if it ever matters (one more mode solve).

**1c (L_π vs. candidate MMI body width):**

| mm_width (µm) | n_eff[0] | n_eff[1] | Δβ (rad/µm) | L_π (µm) |
|---|---|---|---|---|
| 6.0 | 2.962350 | 2.956337 | 0.02884 | 108.92 |
| 8.0 | 2.963225 | 2.959841 | 0.01623 | 193.56 |
| 10.0 | 2.963630 | 2.961464 | 0.01039 | 302.34 |
| 12.0 | 2.963851 | 2.962346 | 0.00722 | 435.25 |
| 14.0 | 2.963984 | 2.962878 | 0.00530 | 592.30 |
| 16.0 | 2.964070 | 2.963223 | 0.00406 | 773.21 |

**✅ Job 11693043 COMPLETED 2026-08-05, 3:19:14 wall.** All 6 points landed.
Both `n_eff[0]` and `n_eff[1]` climb toward the slab limit as the body
widens (their difference, hence Δβ, shrinking correspondingly), and
`L_π/W²` is essentially constant across the whole range —
**3.026 → 3.020 µm⁻¹, a 0.2% drift over 6→16 µm** — exactly the `L_π ∝ W²`
scaling self-imaging theory (Eq. 9, Soldano & Pennings) predicts. That
clean a match is itself evidence these are genuine, well-behaved guided
body modes, not an artifact — a real, independent cross-check the earlier
`te_fraction` detour (A10) never had. **These are TE0/TE1-style body modes
by the same reasoning A10 gives for why polarization mixing doesn't matter
structurally here** — not independently re-verified per point, consistent
with treating this as a non-issue for symmetric cross-sections throughout
this device.

**Step 1 is now closed.** Feeds directly into Step 2: pick a candidate
`mm_width` and read `L_π` either from this table or via
`L_π ≈ 3.02 × mm_width²` (µm), good to ~0.2% over 6–16 µm.

## Summary — what's free, what's new

| | status |
|---|---|
| 4-input / 4-output MMI geometry | ✅ **already supported, no new class** — `MMI(n_inputs=4, n_outputs=4, …)`, same constructor every other MMI example uses |
| Input **and** output tapers | ✅ **already supported** — `MMI`'s default port constructor gives every port its own `taper_length`/`taper_width`, on both the input and output side, for free |
| 220 nm SOI at O-band | ⚠ **new (platform, band) *pair*, not new physics** — `SOI_220` and O-band each already exist and are separately validated elsewhere in this repo, but never together (`SOI_220` has only ever been run at C-band; O-band has only ever been run on `SOI_160`) |
| Si/SiO2 material model at O-band | ✅ **already validated, zero new work** — unlike the SiN 1×4 plan, where "the material model was the whole job," this device's material system is the same Sellmeier model the 1×2 MMI already runs at O-band on `SOI_160`. Nothing in `beef/materials/stack.py` needs to change. |
| Mirror symmetry (paper Eq. 8, `y_i = y_{i+4}`) | ✅ **CONFIRMED 2026-08-05** — not just the default behaviour, checked directly on a real `MMI(4,4)` instance — see A5 |
| 4-port aperture isolation (both ends) | ⚠ **needs `use_neighbor_aperture=True`**, not `half_domain_aperture` — see A6, precedent already exists in this repo (`mmi_1x4_hybrid.py`) |
| CMRR / phase-error / imbalance / excess-loss metrics (paper Eqs. 10–15) | ✅ **Step 7 landed 2026-08-05** — `SMatrix.imbalance_db`/`excess_loss_db`/`cmrr_db`/`cross_excitation_phase_diff_deg`, 13 new tests, synthetic data only. The *math* primitive for phase exists now; the *physical* risk A8 flags (do two independently-run FDTD excitations share a consistent phase reference?) is **not** resolved by this and still needs a real check once Step 6 produces actual solves. |
| Self-imaging design equations (paper Eqs. 7–9) | ✅ **reusable directly** — they are general self-imaging relations, not specific to weakly-guided silica — see A3 |
| The paper's micro-GA optimizer | ❌ **not ported** — deliberate substitution, see A3/Step 3. `beef/sweep/` has no optimization loop, only `parallel.py` (a map-over-workers helper) |

So the ordered work is: **settle the (platform, band) cross-section → derive
starting dimensions from the paper's own equations → refine with a cheap EME
sweep (not a GA) → geometry/aperture → interface mode audit → three example
scripts → a new phase-aware metrics helper → the four checks.** Steps 1–2 are
cheap mode solves; step 3 is the one genuine simulation campaign before any
example script is written; nothing here is blocked on a library gap the way
the SiN 1×4 plan was blocked on the material model, or the mux was blocked on
`CompoundGeometry`.

---

## Part A — findings that set the plan

### A1. The geometry is already fully general — MEASURED (read from source)

`beef/geometry/mmi.py:27` — `MMI`'s own docstring: *"n_inputs=1, n_outputs=2
gives a 1x2 splitter; n_inputs=2, n_outputs=2 a 2x2; n_inputs=1, n_outputs=4 a
1x4; etc."* `n_inputs=4, n_outputs=4` is the same constructor path, no special
case anywhere in the class. `_port_ys(n, spacing)` (`mmi.py:74`) centres *n*
ports symmetrically about y=0 for any *n*; `MMIPort` (`mmi.py:16`) already
carries a per-port `taper_length`/`taper_width`, so **both** input and output
tapers are drawn by the default constructor with no extra plumbing — the
"input and output tapers" the user asked for are not a feature to add, they
are what `MMI` already does for every port on every device in this directory.

### A2. `SOI_220` × O-band is a new *pair*, not new physics — MEASURED (read from source)

`beef/materials/stack.py:711-712`:

```python
SOI_160 = MaterialStack(wg_height=0.16, sz=2.0)   # 160 nm Si (default in this codebase)
SOI_220 = MaterialStack(wg_height=0.22, sz=2.0)   # 220 nm Si (standard commercial SOI)
```

Every existing `SOI_220` device (2×2 MMI, directional coupler) runs at
C-band. Every existing O-band device (1×2 MMI, modal multiplexer) runs on
`SOI_160`. **No example in this repo has ever combined 220 nm Si with
O-band.** This repo's own precedent (`modal_multiplexer_validation.md` gap 2:
*"It had to be measured fresh on SOI_160 rather than copied from the
coupler's res 80, since a 160 nm core confines more weakly than 220 nm"*)
treats a new (platform, band) pair as needing its own cross-section check
even when both halves are individually well understood — the same standard
applies here, in the *other* direction: `sz=2.0` and the existing SOI_220
examples' `waveguide_width=0.45` were chosen (or inherited) for **C-band**;
neither is validated at O-band's shorter wavelength, where the same physical
core confines more tightly. **Do not reuse `waveguide_width=0.45` from
`mmi_2x2_hybrid.py` without checking it is still single-mode at O-band** —
tighter confinement at shorter λ makes a mode *more* guided, which is exactly
the direction that can turn a single-mode C-band width multi-mode at O-band.

The Si/SiO2 Sellmeier model itself needs nothing new: `_SELLMEIER["si"]`/
`["sio2"]` already cover O-band (the 1×2 MMI proves it in production), so
this device's material work is zero, unlike the SiN 1×4 plan where the
material model was blocking everything.

### A3. This device's physics regime does NOT need the paper's micro-GA — ARGUED, from the paper's own text

The paper's central problem statement (its own §I): self-imaging theory
"cannot be assumed in weakly guided waveguides," so it needed a real-coded
micro-GA to hit CMRR/phase-error targets on **weakly-guided ZrO2-doped
silica** (Δ = 5.5%, Table 1's optimized `w = 2.98 µm`, `L = 459 µm`). Our
platform is the opposite regime: 220 nm Si core in SiO2 clad is
**strongly-guided, high-index-contrast** — the same Δn² ≈ 10 vs. silica's
much smaller contrast already established in this repo's own SiN-vs-Si
comparison (`mmi_1x4_sin300_validation.md` A7: "Si/SiO2 Δ(n²) = 10.0094" vs.
nitride's 1.8999, and silica's contrast is smaller again). Self-imaging
theory is expected to hold well here — which is exactly the regime the
paper's own reference [4] (Halir et al., *Opt. Lett.* 2011, "High-performance
90º hybrid based on a silicon-on-insulator multimode interference coupler")
addresses: a real SOI 90° hybrid, designed from self-imaging without a GA.

**Decision for this plan: use the paper's own general self-imaging equations
(Eqs. 7–9, general — not weakly-guided-specific) to get starting dimensions,
then refine with a coarse EME width/length sweep (this repo's
`solve_section_length_sweep()` pattern, already used by `mmi_2x2_hybrid.py`
and `mmi_1x4_hybrid.py` for exactly this kind of cheap length refinement) —
not a ported micro-GA.** `beef/sweep/` has no optimizer (`parallel.py` is a
worker-pool map, not a GA), so porting one would be new infrastructure, not
device work. **Flagged as a judgment call, not a silent one** — if the
self-imaging + sweep approach does not hit the paper's own quoted system
requirements (CMRR ≥ 20 dBe, |Δφ| ≤ 5°, from its refs [4,5]) within a
reasonable sweep budget, porting a real optimizer becomes the next honest
step rather than hand-tuning forever.

### A4. Only ports 1 and 3 are ever driven — MEASURED, read from the paper

Paper §III: *"input waveguides 2 and 4 are not used for optical 90º hybrids
[3]"* — `I = {1,3}`, `O = {5,6,7,8}` in the paper's own notation. Mapped to
this codebase's 0-indexed `input_idx`/`output_idx`: only `input_idx ∈ {0,2}`
of the four input ports are ever excited; all four `output_idx ∈ {0,1,2,3}`
are measured for both excitations. A5 below means the other two input ports
(`input_idx ∈ {1,3}`) get solved anyway as a side effect of the mirror-symmetry
reuse — not wasted, just not separately required by the paper's own figures
of merit.

### A5. Mirror symmetry (paper Eq. 8) is already the default, not a feature to build — ARGUED from source + the paper

Paper Eq. (8) plus its own follow-up sentence: *"we assume y_i = y_{i+4} for
i = 1,…,4"* — each input port sits at the **same** y-position as its
same-indexed output port, and (implicitly, since both sides share `w`) the
same access-waveguide width. That is precisely what `MMI`'s default
constructor already gives for free: `_default_ports(n)` builds `self.inputs`
and `self.outputs` from the *same* `waveguide_width`/`taper_length`/
`taper_width`/`waveguide_offset` call (`mmi.py:58-66`), so as long as this
device is built via the plain `MMI(n_inputs=4, n_outputs=4, waveguide_width=…,
taper_length=…, taper_width=…, waveguide_offset=…)` call — no custom
`inputs=`/`outputs=` port lists — the input and output faces are
*structurally* identical, exactly the symmetry the paper's Eq. 8 states.

**This is the direct answer to "check the paper for symmetry between input
and output ports."** It also means this device is a strong candidate for the
same left-right mirror-symmetry reuse `mmi_2x2_hybrid.py` already exploits
(module docstring: *"the output section with output_idx=k is exactly the
x-mirror of the input section with input_idx=k"*) — generalized from 2 ports
to 4:

- Solve the **input** section once per port index (`input_idx = 0,1,2,3` —
  **4 FDTD solves**, not 8), plus **1 EME body solve** (length-swept
  analytically, same pattern as `mmi_2x2_hybrid.py`/`mmi_1x4_hybrid.py`).
- Every output-side transfer matrix is the corresponding input matrix's
  transpose — the output section is **never separately simulated**, same as
  the 2×2's own doc.
- To get the paper's full `S_ji(λ)` table for `j ∈ {5,6,7,8}` (all four
  outputs) and `i ∈ {1,3}` (the two driven inputs), cascade: `M_body @
  M_in_i` then `M_out_j.T @ (that)` for each of the 4 output indices — 4
  input solves cover **all 4×2 = 8** entries the paper's fitness function
  needs, at the cost of 4 FDTD solves total (vs. up to 8 if input/output
  were solved independently, or 2 if only the driven ports were solved and
  outputs were separately simulated).

**Not yet confirmed, flagged as a real risk rather than assumed to hold:**
the reuse is only valid if the *body* is genuinely left-right symmetric too
(a plain rectangular MMI slab, `ΔL` folded symmetrically) and if
`body_stub`/`half_domain_aperture` settings on both ends match exactly —
worth a direct check once geometry exists (Step 4), not before.

> **✅ CONFIRMED 2026-08-05** — checked directly, not just read from source:
> `MMI(n_inputs=4, n_outputs=4, mm_length=20.0, mm_width=10.0,
> waveguide_width=0.45, taper_length=3.0, taper_width=1.2,
> waveguide_offset=1.5, waveguide_height=0.22, use_neighbor_aperture=True)`
> gives `input_ys == output_ys == [-2.25, -0.75, 0.75, 2.25]`, and every
> input/output port pair matches on `y`, `waveguide_width`, `taper_length`,
> and `taper_width`. `frame_offset == 0.0` as well (equal `taper_length` on
> both sides), confirming the body-centred frame used for the mirror
> reasoning above is also bbox-centred here — no hidden offset between the
> "mirror" picture and the actual coordinate frame `section_prisms()`/
> `section_gds()` use. Not yet checked: the *body itself* being
> left-right symmetric once `ΔL`/taper lengths are chosen for real (Step 2/4)
> — this only confirms the port-level claim.

### A6. Aperture: needs `use_neighbor_aperture=True` on BOTH ends — MEASURED (read from source), by direct analogy

`mmi_2x2_hybrid.py`'s `half_domain_aperture` works only because its two
ports per end split cleanly across a shared y=0 axis. This device has **four**
ports per end, on **both** the input and the output side (unlike the 1×4
SiN device, which only had the problem on its output side) — a y=0 cut
still groups ports {1,2} and {3,4} together on each end, exactly the failure
mode `mmi_1x4_hybrid.py`'s own module docstring identifies for its four
*output* ports:

> *"This device's four ports … don't split that way: a y=0 cut still groups
> ports 0&1 … and 2&3 … together. `use_neighbor_aperture` instead sizes each
> port's aperture to the distance to its OWN nearest neighbor."*

**Decision: `MMI(..., use_neighbor_aperture=True)`, and no
`half_domain_aperture=True` calls on either FDTD section** — same mechanism,
now needed on both faces rather than one. `mmi.py:332-355`
(`source_aperture_at`) already implements this generally; nothing new to
write.

### A7. Interface mode count is unknown — must be measured, not assumed

Every MMI hybrid example in this repo picked its own `N_MODES` from a real
mode-count measurement at its own cross-section, never copied across
devices at face value: 2×2 MMI → 14 (of 15, one unguided dropped), 1×4 SiN →
20 (of 23 guided), modal multiplexer → 14 (of 16 requested, cut down after
measuring only 14 are actually reachable). This device's body is wider
still (four access waveguides plus tapers on each end, vs. two or one), so
**expect more supermodes than any MMI example so far**, closer to the modal
multiplexer's territory than to either existing MMI's. `run
audit_interface_modes.py`-style — a real mode solve at the settled
cross-section (Step 5) — before picking a number, and budget for the
duplicate-band risk the modal multiplexer hit at high mode count (the
library-level guards, `_iter_distinct_bands()`/`compute_modes()`'s
`guard_duplicate_modes`, are on by default repo-wide since 2026-07-31, so
this device inherits the fix for free — but the *measurement* of how many
distinct modes actually exist here is still device-specific and unmeasured).

### A8. This is the first device whose own metric needs S-parameter phase, not just power — MEASURED (read from source)

`beef/results/smatrix.py` stores `_data[(port, in_mode, out_mode)]` as
**complex** amplitudes (`SMatrix.set`/`coeff` docstring: *"All coefficients
are complex; power = |coeff|^2"*), but every existing helper —
`transmission()`, `insertion_loss_db()`, `splitting_ratio()`,
`check_passive()`, `total_power()` — reads only `|coeff|^2`. Nothing in
`beef` currently computes a phase difference between two `SMatrix` entries.

That is exactly what this device's headline metrics need. From the paper's
own equations, in this codebase's notation (`i ∈ {0,2}` = driven inputs,
`j ∈ {0,1,2,3}` = the four outputs, `S = smatrix.coeff(port=f"out{j}",
in_mode=..., out_mode=...)` per whatever port-naming convention Step 6
settles on):

- **Imbalance** (paper Eq. 10): `IB_i(λ) = −10·log10[min_j|S_ji|² /
  max_j|S_ji|²]` — power-only, `transmission()` already gives this per port;
  only the min/max-over-outputs reduction is new.
- **Excess loss** (Eq. 11): `XL_i(λ) = −10·log10[Σ_j|S_ji|² / P_in(λ)]` —
  power-only, a sum of `transmission()` calls.
- **CMRR_I / CMRR_Q** (Eqs. 13–14): `−20·log10[||S_5i|²−|S_8i|²| /
  (|S_5i|²+|S_8i|²)]` and the 6/7 analogue — power-only, same shape as
  `splitting_ratio()` but as a *difference* over a *sum*, not a ratio; a
  small new reduction, not a new physical quantity.
- **Phase error** (Eq. 15): `Δφ_j(λ) = ∠S_j3(λ) − ∠S_j1(λ) − (φ̃_j − φ̃_5)`
  with target angles `φ̃_5=−45°, φ̃_6=225°, φ̃_7=45°, φ̃_8=135°` — **this one
  is genuinely new**: it needs `np.angle()` of two *different-excitation*
  `SMatrix` entries, on a shared, physically meaningful phase reference.

**Open risk, worth its own check (Step 8a) before trusting any phase-error
number**: A5's mirror-symmetry reuse builds `M_in_1` and `M_in_3` from two
*independent* `solve_section()` calls (different `input_idx`) sharing one
body solve. Every prior device in this repo only ever needed the *magnitude*
of its transfer matrices to agree across separate excitations
(`check_passive()`, splitting ratios) — nothing has ever checked that two
independently-run FDTD eigenmode decompositions preserve a consistent
**absolute phase reference** relative to each other. If they do not (e.g. an
arbitrary per-call phase convention in `get_eigenmode_coefficients()`),
`Δφ_j` would be corrupted while every power-based check
(`check_passive()`, imbalance, excess loss, CMRR) still looks perfectly
healthy — exactly the "each individual piece looks fine" failure mode the
four checks exist to catch, just for a quantity none of them currently
touch. **This needs a dedicated, direct check**, e.g. verifying `Δφ_j` from
a full single-domain FDTD reference (check 4) reproduces the hybrid's own
`Δφ_j` before trusting the hybrid number at all — a stronger requirement
than check 4 has ever had to meet on a magnitude-only device.

**Decision needed, not made here:** whether these four helpers become
script-local functions in the example (mirroring
`beef.solvers.fdtd.dbeta_from_phases()`'s precedent of a small, reusable,
device-independent helper) or a proper addition to `SMatrix`
(`SMatrix.phase()`, `SMatrix.cmrr()`, …). Given `beef.mode_identity`'s and
`dbeta_from_phases`' own precedent of factoring out device-independent
numeric helpers rather than leaving them buried in example scripts, and
given CMRR/phase-error are generic coherent-receiver metrics that could
plausibly matter to a future device, **leaning toward small, tested
additions to `beef/results/smatrix.py`** — but flagged here for a decision
before Step 7, not decided unilaterally.

> **✅ DECIDED AND LANDED 2026-08-05 — small, tested additions to
> `SMatrix`.** `beef/results/smatrix.py` gained four methods:
>
> - `imbalance_db(ports, in_mode=1, out_mode=1)` — Eq. 10, generalized from a
>   fixed pair to any number of output ports.
> - `excess_loss_db(ports, in_mode=1, out_mode=1, p_in=1.0)` — Eq. 11.
>   **Positive = loss**, the paper's own convention and the *opposite* sign
>   of this class's existing `insertion_loss_db()` — documented in the
>   method's own docstring so the two are never mixed in one report.
> - `cmrr_db(port_a, port_b, in_mode=1, out_mode=1)` — Eqs. 13/14's formula,
>   generic to any two-port balanced pair (the caller picks which two ports
>   make it CMRR_I vs. CMRR_Q).
> - `cross_excitation_phase_diff_deg(other, port, in_mode=1, out_mode=1)` —
>   the genuinely new primitive. Takes a *second* `SMatrix` explicitly
>   (`self`/`other` built from different input-port excitations of the same
>   device) rather than being a single-instance method, because **every
>   existing SMatrix-producing example in this repo holds exactly one input
>   excitation's results** (confirmed by reading `mmi_2x2_hybrid.py`'s own
>   single `smat = SMatrix(...)` for its one driven arm) — there is no
>   existing convention for "one SMatrix, two input ports" to extend, so this
>   had to be designed, not copied. Reproduces the paper's own
>   `phi_j(λ) = ∠S_j3(λ) − ∠S_j1(λ)` exactly when called as
>   `smat_in1.cross_excitation_phase_diff_deg(smat_in3, port=...)`.
>   **Deliberately stops short of the paper's full Eq. 15** — the reference-
>   port subtraction and hardcoded target angles (`φ̃_5=−45°`, etc.) are
>   device-specific, and building them from this primitive is Step 6/7's
>   example-script job, not a library concern. Output is wrapped to
>   `(-180, 180]`, with an explicit docstring warning that a target stated in
>   `[0, 360)` (`φ̃_6=225°`, `φ̃_8=135°`) must be compared modulo 360, not by
>   direct subtraction.
>
> **13 new tests** in `tests/test_smatrix.py::TestCoherentReceiverMetrics`,
> synthetic data only, hand-computable (e.g. `cmrr_db` on a 0.6/0.4 W split
> asserted against `-20·log10(0.2) = 13.9794` dBe by hand) — including a
> wraparound case (`170°` vs. `-170°` excitations correctly returning `+20°`,
> not `-340°`) and an in-range replica of two of the paper's own Eq. 15
> targets (`φ̃_5=-45°`, `φ̃_7=45°`). Full fast suite re-run clean: **639
> passed, 19 skipped** (up from 626 passed pre-change — 13 new, zero
> regressions).
>
> **What this does NOT close**, restated because it is easy to misread "the
> metrics are implemented" as "the metrics are trustworthy": the *physical*
> risk this finding opened with — whether two independently-run
> `solve_section()` FDTD excitations preserve a consistent absolute phase
> reference relative to each other — is untouched by this. These four
> methods are correct arithmetic over whatever `SMatrix` coefficients they
> are given; they cannot tell a real, physically-meaningful phase from an
> arbitrary per-call convention. That check still needs Step 6's real solves
> and is still Step 8a's job.

### A9. No published ground truth for THIS device — ARGUED, by construction

The paper's own Table 2/Fig. 3 numbers (CMRR ≥ 28.9 dBe, phase error ≤
2.52°, excess loss ≤ 0.68 dB) describe a **weakly-guided silica, C-band**
device with dimensions in the hundreds of µm (`L = 459 µm`). This device
will be **strongly-guided silicon, O-band**, with dimensions expected to be
an order of magnitude smaller (shorter `L_π` from much larger index
contrast). **Those numbers are a design-methodology reference (the
equations, the port topology, the ports-2-and-4-unused convention, the
system-requirement thresholds it itself cites), not a numerical target** —
the same caution the modal multiplexer doc already states about its own
"Seed vs. Optimized" reference figure (`modal_multiplexer_validation.md`,
"Reference design" section): don't quote a paper's numbers as this device's
ground truth when the underlying physics differs. What DOES carry over
directly and should be checked against: the **typical system requirements**
the paper itself cites from refs [4,5] — **CMRR ≥ 20 dBe, |Δφ_j| ≤ 5°** —
since those are receiver-system specs, not properties of the silica design.

---

## Part B — the ordered steps

Each step names its **deliverable** and its **evidence** (what makes it
done), same convention as `mmi_1x4_sin300_validation.md`.

### Step 0 — settle the band explicitly

**Deliverable.** `WL_MIN = 1.26`, `WL_MAX = 1.36 µm` (O-band), matching the
existing O-band devices in this repo (1×2 MMI, modal multiplexer) rather
than deriving new bounds — stated explicitly per this repo's own standard
("state the choice … do not let it be implicit," `mmi_1x4_sin300_validation.md`
step 4c). Centre wavelength 1.31 µm for any single-λ mode solve (`β0`, `β1`
for `L_π`, per paper Eq. 9).

**Evidence.** This line, plus the constant in whatever script Step 6
produces.

### Step 1 — cross-section design (cheap mode solves, no FDTD)

**1a. Single-mode access-waveguide width, at O-band, on `SOI_220`.** Sweep
strip width and count guided modes (same convention as every mode-count
table in this repo: `Re(n_eff) > n_clad`) at λ=1.31 µm. **Do not inherit
`waveguide_width=0.45`** from `mmi_2x2_hybrid.py` (A2) — that number was
never independently validated even at its own C-band, and O-band's shorter
wavelength confines more tightly. Deliverable: `waveguide_width`. Evidence:
a mode-count-vs-width table with a stated single-mode window.

**1b. Vertical/lateral box convergence (`sz`, `sim_height`).** `SOI_220`
ships `sz=2.0`, inherited from C-band devices. Sweep `sz` at the chosen
width, O-band, and confirm `|Δn_eff|` convergence — likely converges fast
given Si/SiO2's strong confinement (contrast with the SiN saga's slow
convergence at its own near-cutoff width), but state the number measured,
don't assume `sz=2.0` transfers unchanged.

**1c. `L_π` and the self-imaging length (paper Eqs. 7, 9).** Solve
`β0`, `β1` (fundamental and first-order MMI-body modes) at a candidate
`mm_width`, λ=1.31 µm: `L_π = π/(β0−β1)`. This number is expected to be far
shorter than the paper's own 459 µm silica body, by roughly the ratio of
index contrasts — compute it, don't guess the scale from the paper's number.

### Step 2 — starting dimensions from the paper's own equations

**Deliverable.** A first-pass `(W, ΔL, w, y_1..y_4)` set from Eqs. 7–9,
Table 1's parameter *roles* (not its silica-tuned *values*): `mm_width (W)`,
`ΔL` (length deviation), `waveguide_width (w)`, per-port `Δy_i` (kept at 0
initially — a reasonable start given A5's mirror symmetry already handles
the topology-level placement `y_i = 2i-5)/8 · W`, per Eq. 8, with `Δy_i` as
a later fine-tuning knob only if Step 3's sweep needs it). `L = (3/4)L_π +
ΔL` (Eq. 7) with `ΔL` initially 0, refined in Step 3.

**Evidence.** The computed starting values, with `L_π` from Step 1c, not
copied from Table 2.

### ✅ Step 2 DONE 2026-08-05 — first-pass dimensions computed and geometry-checked

**Candidate `mm_width`: `W = 6.0 µm`** — the narrowest (most compact) point
actually measured in Step 1c, not interpolated or extrapolated. Chosen
because a strongly-guided platform's whole advantage over the paper's own
weakly-guided silica (A3/A9) is compactness, so the starting point should
lean into that rather than default to the middle of the swept range.
Step 3's sweep may move this — this is a starting point, not a final value.

**A clean derivation, checked directly rather than assumed:** paper Eq. 8,
`y_i = ((2i−5)/8)·W` for `i=1..4`, is *exactly* what `beef`'s own default
`MMI` port generator (`_port_ys`) produces when
**`waveguide_offset = W/4`**. Algebra: `_port_ys` gives
`y_j = offset·(j−1.5)` for 0-indexed `j=0..3`; matching term-by-term against
Eq. 8 (`i=j+1`) gives `offset = W/4` for every port, not just one. Verified
on a real `MMI(4,4)` instance, not just by hand:

```
paper Eq.8 y_i        : [-2.25, -0.75, 0.75, 2.25]
beef default input_ys : [-2.25, -0.75, 0.75, 2.25]   # W=6.0, offset=1.5=W/4
```

This means Step 4's geometry needs no custom `inputs=`/`outputs=` port
lists to satisfy Eq. 8 — the plain constructor with `waveguide_offset=W/4`
does it exactly, on both faces (A5).

**Computed first-pass parameter set:**

| parameter | value | source |
|---|---|---|
| `w` (`waveguide_width`) | 0.35 µm | user's direct decision (Step 1) |
| `W` (`mm_width`) | 6.0 µm | candidate, measured point (Step 1c) |
| `waveguide_offset` | 1.5 µm (`= W/4`) | derived, reproduces Eq. 8 exactly |
| `y_1..y_4` | −2.25, −0.75, 0.75, 2.25 µm | derived, Eq. 8, Δy_i = 0 |
| `L_π` | 108.9195 µm | **measured**, job 11693043, at `W=6.0` |
| `ΔL` | 0 µm | initial (Eq. 7); Step 3's refinement knob |
| `L` (`mm_length`) | 81.6896 µm | derived, Eq. 7: `(3/4)·L_π + ΔL` |
| `taper_width` | 0.8 µm | **placeholder** — `MMI`'s own class default, not measured |
| `taper_length` | 1.0 µm | **placeholder** — `MMI`'s own class default, not measured |

**taper_width/taper_length are flagged, not silently assumed.** Unlike
every other row, these two are not derived from the paper or from a
measurement — they are `MMI.__init__`'s own neutral defaults, carried
forward exactly because Step 2's own plan named them "TBD — measure clean-
transition length, not a default copy." Real values need either a quick
adiabatic-transition check or fall out of Step 3's sweep; do not treat 0.8/
1.0 as validated.

**Geometry feasibility, checked directly (no meep — pure `MMI` construction):**

```
mm_length (Eq.7, ΔL=0) = 81.689625 um
total_length           = 83.689625 um
bbox                    = (83.689625, 6.0)
widest port width = 0.8 um vs spacing (offset) = 1.5 um -> clearance = 0.700 um per side
OK -- no taper overlap, checked directly
```

0.7 µm of clearance per side between adjacent taper edges at their widest
point — comfortable margin, not a near-miss. **Step 2 is closed; Step 3
(EME width/length refinement sweep) is next**, starting from this exact
parameter set.

### Step 3 — EME width/length refinement (replaces the paper's μGA — see A3)

**Deliverable.** A coarse grid or 1D sweep over `mm_width`/`ΔL` (and, only
if needed, `w`), using `EMESolver.solve_section_length_sweep()` on the
uniform body slab — the same cheap length-sweep machinery
`mmi_2x2_hybrid.py`/`mmi_1x4_hybrid.py` already use for their own body
refinement (a length sweep is one mode solve plus analytic phase
propagation per length, not a re-solve per point). Target: equal splitting
among the four outputs for each driven input, and the four relative phases
(Eq. 15's `φ̃_j`) within the cited system requirement (`|Δφ_j| ≤ 5°`) across
the O-band.

**Evidence.** A splitting-ratio-vs-length and phase-error-vs-length plot,
with the chosen `(mm_width, ΔL)` stated and justified (same "state the
tradeoff, don't silently pick the argmin" discipline
`mmi_1x4_sin300_validation.md`'s own `MM_LENGTH=206.0` decision used — that
device explicitly chose a *non*-optimal point for footprint/robustness
reasons and said so). **If this sweep cannot reach the cited system
requirement**, that is the trigger for A3's fallback (a real optimizer),
not a reason to loosen the requirement quietly.

#### In progress 2026-08-06

Concretely, Step 3 needs: (3a) a real body mode-count audit at `W=6.0`
(Step 1c only ever solved 2 modes there, just enough for `L_π`) to pin
`N_MODES`; (3b) 4 real FDTD input-section solves, one per port index — not
2, even though only ports 1/3 are driven (A4), because all 4 are needed as
mirror-reused output matrices too (A5); (3c) the analytic body length
sweep around Eq. 7's point; (3d) evaluate every length via the new
`SMatrix` methods (Step 7/A8) and pick a `ΔL`. Script:
`examples/verification/python/mmi_4x4_90hybrid_body_sweep.py` (drafted,
not yet run — waiting on 3a for `N_MODES`). Its own `sub_region_configs`
mirror-pairs argument was caught and fixed before ever running (`{2: k}`,
wrong — section indices only go 0–2 in a 3-section `HybridSolver`; ports go
0–3 — corrected to `{2: 0}`, section 2 always mirrors section 0 regardless
of which port a config excited).

**3a hit a real, informative cost surprise.** First attempt, job
**11706986**, requested 16 modes (guessing "a few times" the 2-mode
point's 874.6 s) — **TIMED OUT at the full 60 min with ZERO output**, not
even a partial result. Whatever the true per-mode cost curve is here, "a
few times more for 8× the modes" was wrong. Resubmitted as **11707835**
requesting only 8 modes (a calibration point, not the final answer) at a
90-min budget, specifically to measure the real scaling before guessing
again. `W=6.0` is wider than the 2×2's own `mm_width=2.8`, and O-band's
shorter wavelength packs in more modes than C-band at the same width, so
this cross-section plausibly supports considerably more guided modes than
any existing device's own `N_MODES` — worth knowing before picking a
number, not assuming a small one transfers.

**✅ 3a MEASURED 2026-08-06 — job 11707835 COMPLETED (55:58).** All 8
requested modes guided, n_eff 2.962350 → 2.833452, cladding at 1.446804 —
healthy, still-widening margin, no sign of approaching cutoff. Real cost:
**3324.7 s for 8 modes vs. 874.6 s for 2 (job 11693043) — ~416 s/mode,
close to linear**, not the "few times more" originally guessed. This
cross-section evidently supports considerably more than 16 guided modes
(the band-to-band spacing is still *growing*, not shrinking, at mode 8).

**Consequence: `solve_section_length_sweep()` costs one full mode solve
PER WAVELENGTH, so `N_MODES=16 × N_WL=11`** (this repo's usual convention,
which the draft script had copied uncritically) **would have cost roughly
20+ hours for the body alone.** Rescoped before submitting anything real:

- **`N_MODES=8`**, reusing job 11707835 directly — a deliberate coarse
  first pass, not an exhaustive audit (MMI self-imaging is typically
  dominated by low-order modes). Flagged for revisiting with a real,
  higher-mode-count audit later if this pass needs refinement.
- **`N_WL=3`**, at `{1.26, 1.31, 1.36} µm`** — this is not an ad hoc cut;
  it mirrors the **paper's own fitness convention exactly** (Sec. IV / Eq.
  12: `Λ = {1520, 1550, 1580 nm}`, 3 points). The other MMI examples in
  this repo used 11 points because their own cross-sections were cheap
  enough to afford it; this one is not, and 3 has real precedent here.

Brings the body sweep to ~3×55 min serial (or ~55 min array-parallelized,
one wavelength per task — this repo's established pattern for this exact
cost shape). `mmi_4x4_90hybrid_body_sweep.py` updated accordingly, plus
`RESOLUTION`/`DECAY_TOL`/`N_WL` made env-var overridable
(`MMI4X4_RES`/`MMI4X4_DECAY_TOL`/`MMI4X4_NWL`) for a smoke pass — same
pattern `modal_multiplexer_fdtd.py`'s own `MUX_FDTD_RES` smoke test used.

**Smoke test submitted, job 11708372** (`res=15`, `decay_tol=1e-9`,
`N_WL=1`) — the input-section FDTD path (`FDTDSolver.solve_section()`,
real time-stepping on this 4-port tapered access geometry) has never been
exercised at all; only the pure EME mode-count path has (11693043/
11707835). Per this repo's own two-step discipline for expensive,
never-run scripts.

**✅ Smoke test PASSED 2026-08-06, job 11708613 (18.4 s).** Geometry
correct (`mm_length=81.6896`, matches Step 2 exactly), all 4 input-section
solves and the body length-sweep produced correctly-shaped output
(`(1,8,1)` per input, `(81,1,8,8)` for the body), the full evaluate/save/
plot pipeline ran without error. The printed numbers themselves (best ΔL
at the search-range edge, `max|Δφ|=86.8°`) are not meaningful — res=15,
decay_tol=1e-9, N_WL=1 is a pure plumbing check, exactly as designed, and
say nothing about the real device.

**A real, measured confirmation the `srun --mpi=pmi2` fix mattered, not
just a precaution.** Job **11708372** — the batch script's pre-fix version
(`module load openmpi` + bare `python`, no `srun --mpi=pmi2`, caught and
fixed before submission per the mux's own 42x-slowdown precedent — see the
batch script's own comment) — happened to run too (submitted moments
before the fix, from the same background drafting process) and returned
**byte-identical results** in **592.6 s, a measured 32× slower** than the
fixed version's 18.4 s. Correct answer either way at this tiny smoke
scale (meep runs fine single-process, just far slower without proper
launching) — but the penalty compounds with real settings (res 15→40 is
already ~51× costlier per the res⁴ 3-D FDTD scaling; a further ~32×
launcher penalty on top of that would plausibly turn a feasible 12h budget
into an infeasible one). Worth having a real number for this, not just the
mux's own precedent to point to.

**Real production sweep submitted 2026-08-06, job 11708624**
(`run_mmi_4x4_90hybrid_body_sweep.batch`, 48 ranks / 2 nodes, `-t 720`,
`gts-sr108-ece`, `srun --mpi=pmi2`) — production settings, no env
overrides: `RESOLUTION=40`, `DECAY_TOL=1e-11`, `N_WL=3` (the paper's own
Eq. 12 wavelength convention). Budgeted generously rather than tightly,
matching `run_mmi_2x2_hybrid_dispersive.batch`'s own real-hybrid-production
scale — 18.4 s at res=15 does not extrapolate to a real-settings estimate.
Every `solve_section()`/`solve_section_length_sweep()` call already passes
`save_dir=save_dir`, so a timeout partway through the 4 sequential
input-port solves would leave completed ones cached — a plain resubmission
resumes rather than restarts; array/merge reshaping (this repo's own
established fallback for exactly this cost shape, per the 2×2's own gap-4
sweep) is the next move only if a plain resubmit turns out not to be
enough, not the first one.

**✅ Job 11708624 COMPLETED 2026-08-06 (3:53:59) — and the first read of
the result triggered a real, diagnosed, and fixed bug (A11).**

> ### A11. A global phase-sign-convention mismatch, found, diagnosed with cached data (no new FDTD), and fixed
>
> **The observation.** The raw result: `Best ΔL = -1.50 µm: max|Δφ|=
> 179.857°, max imbalance=4.9044 dB`. Loading the full 81-point sweep from
> `body_sweep.npz`: `max_dphi` sat between **172.4° and 180.0°across the
> ENTIRE ±10 µm ΔL range** (mean 179.1°) — essentially the worst possible
> value, nearly invariant to a parameter (body length) that should strongly
> move real relative phase. `max_imbalance`, by contrast, varied sensibly
> over the same sweep (4.9–39.6 dB) — a magnitude-only metric behaving
> normally while a phase-only metric sits pinned near its theoretical
> extreme is the signature of a **global sign-convention mismatch**, not a
> per-length physical result.
>
> **Diagnosed using only cached data — zero new FDTD/cluster time.** Every
> `solve_section()`/`solve_section_length_sweep()` call in the sweep script
> already passes `save_dir`, so the 4 real FDTD input-port matrices and the
> body length-sweep matrix were all pulled back with `cache_only=True` in
> under a second each, on the login node.
>
> First checked A8's own predicted risk directly: do independently-run
> `solve_section()` calls preserve a consistent absolute phase reference?
> Compared the mirror-image port pair (`input_idx=1`/`2`, at y=∓0.75) across
> all 8 body modes — **exactly** the physics a mirror-symmetric slab
> predicts: even modes (0,2,4,6) show 0.00° phase difference between the
> pair, odd modes (1,3,5,7) show ±180.00°, magnitude ratios 1.0000, all to
> 2–3 decimal places. **This ruled out A8's cross-call phase-reference risk
> directly** — separate `solve_section()` invocations do preserve a
> consistent absolute phase reference; that was not the bug.
>
> **Then compared the actual measured `phi_j = ∠S_j3 − ∠S_j1` against the
> paper's own target angles directly, not just the wrapped residual:**
>
> | j (port) | measured `phi_j` (wl=1.26) | paper target `φ̃_j` | `-φ̃_j` |
> |---|---|---|---|
> | 0 (out5) | +46.25° | −45° | **+45°** |
> | 1 (out6) | +128.91° | 225° | **135°** (225 mod 360, sign-flipped) |
> | 2 (out7) | −45.97° | 45° | **−45°** |
> | 3 (out8) | −132.23° | 135° | **−135°** |
>
> Every measured value matches **the *negative* of its target** to within
> 1–6° — not coincidence at this consistency across all four independently-
> measured ports. **Conclusion: `meep`/`meow`'s own S-parameter phase
> convention is conjugated relative to the paper's implicit one** (plausibly
> a `exp(-iωt)` vs `exp(+iωt)` time-harmonic convention difference, a known,
> standard source of exactly this kind of uniform phase-difference sign
> flip) — a property of the solver pairing, not a device design failure and
> not a bug in `cross_excitation_phase_diff_deg` itself (which correctly
> computes `angle(other)-angle(self)` exactly as documented; it was called
> with the paper's literal `self=port1, other=port3` ordering, which is
> exactly the ordering that turns out to need negating here).
>
> **Fix, verified against the same cached data (recomputed in 5.5 s, no new
> FDTD):** swap `self`/`other` (equivalently, `φ_j = ∠S_j1 − ∠S_j3`, the
> negative of the paper's literal Eq. 15 as written — necessary because of
> the convention mismatch just diagnosed, not a contradiction of the
> paper). `mmi_4x4_90hybrid_body_sweep.py` corrected accordingly, with the
> full diagnostic trail recorded in the code comment at the fix site so a
> future reader sees the *why*, not just the *what*.
>
> **Corrected result — a real, clean local optimum, not noise:**
>
> | ΔL (µm) | max|Δφ| | max imbalance |
> |---|---|---|
> | −4.00 | 68.9° | 9.79 dB |
> | −2.00 | 11.9° | 6.21 dB |
> | **−1.50** | **8.28°** | **4.90 dB** |
> | −0.75 | 8.46° | 7.08 dB |
> | 0.00 (pure Eq. 7, no correction) | **7.35° (global min)** | 10.44 dB |
> | +1.00 | 12.4° | 20.27 dB |
> | +2.00 | 138.2° | 15.60 dB |
>
> A smooth trough, not scatter — `max|Δφ|` falls monotonically from 68.9°
> at ΔL=−4 to a 7.3–8.5° minimum region between −1.75 and 0, then rises
> again past +1. **This is real physics: self-imaging theory (A3's
> hypothesis) gets this first-pass design — 8 body modes, unoptimized
> `taper_width=0.8`/`taper_length=1.0` placeholders, a pure analytic Eq. 7
> starting length — to within 2.3–3.3° of the paper's own cited 5°
> requirement, without any optimizer.** Neither point in this ±10 µm sweep
> quite clears 5° (global min 7.345° at ΔL=0.00, but with worse imbalance
> there), and there's a genuine trade-off between the two metrics near the
> optimum — worth refining (wider ΔL search, more body modes past the
> coarse `N_MODES=8`, or real taper dimensions replacing the placeholders)
> before concluding A3's fallback (a real optimizer) is or isn't needed.
> **Best combined point: `ΔL = -1.50 µm` (`mm_length = 80.1896 µm`),
> `max|Δφ|=8.28°`, `max imbalance=4.90 dB`.**
>
> Full corrected grids in `example_output/mmi_4x4_90hybrid/body_sweep/
> body_sweep_FIXED.npz`; plot regenerated at `body_sweep.png`.

**Carry forward `mmi_1x4_hybrid.py`'s grid-snap lesson explicitly:** that
script's own comment records a **real, previously-shipped bug** — a
"pinned" length argmin'd onto a sweep grid that didn't actually contain it,
0.0208 µm off for weeks. Whatever length this step picks, assert it is an
exact point on the sweep's own grid before shipping it in a script, the way
`mmi_1x4_hybrid.py:_nearest` now does.

### Step 4 — geometry

**Deliverable.**

```python
mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=...,        # Step 2/3
    mm_width=...,         # Step 2/3
    waveguide_width=...,  # Step 1a
    taper_length=...,     # TBD — measure clean-transition length, not a default copy
    taper_width=...,      # TBD
    waveguide_offset=...,
    waveguide_height=0.22,  # SOI_220
    use_neighbor_aperture=True,   # A6 — needed on both 4-port ends
)
```

Built via the **plain** constructor (no custom `inputs=`/`outputs=` lists) —
required for A5's mirror-symmetry claim to actually hold, not just be
plausible.

**Evidence.** A geometry plot (`plot_hybrid_geometry` or a bare `polygons()`
render) confirming: four input tapers, four output tapers, no taper overlap
at the chosen `waveguide_offset` (the real constraint A1 in the SiN plan's
own 1×4 analogue flags — "the four output tapers at waveguide_offset spacing
must not overlap, which is a real constraint at 1×4 that never binds at
1×2" — doubly true here with tapers on *both* ends), and that the input and
output faces are visually identical (A5's symmetry, checked, not assumed).

**✅ Step 4 DONE 2026-08-06 — already satisfied by Step 3's own output, not
duplicated here.** `mmi_4x4_90hybrid_body_sweep.py`'s `plot_hybrid_geometry`
call already renders exactly this evidence: `example_output/
mmi_4x4_90hybrid/body_sweep/geometry_hybrid.png` shows four tapered input
ports fanning into the body on the left, an identical mirrored fan on the
right (visually confirming A5, not just the port-position check from Step
2), and the FDTD|EME|FDTD section boundaries. At the user's direction,
proceeding with the current best design point (`ΔL=-1.50`, `N_MODES=8`,
placeholder tapers) rather than spending more cluster time chasing the
band-edge CMRR/phase gap (A11) before continuing — that gap is carried
forward as a known, documented limitation, not silently dropped.

### Step 5 — interface mode audit

**Deliverable.** A real mode solve at the settled input/output-taper ↔ body
interface (both ends, confirmed identical per A5 rather than assumed and
audited only once), reporting guided-mode count and the FDTD/EME
correspondence (n_eff nearest-neighbor, then field-overlap η — the two-step
process `modal_multiplexer_validation.md` gap 1 established as necessary,
not the n_eff-only check that gap's own Step 1 showed can be misleading).
`N_MODES` is whatever this measures, not assumed from A7's cross-device
table.

**Evidence.** The audit table (band, n_eff both engines, Δ, guided/not),
plus explicit confirmation of whether either duplicate-band guard fires —
if it does, that is itself the finding (per `mmi_2x2_validation.md` gap 6's
own framing).

### Step 6 — three example scripts

Following every other MMI device's naming: `examples/mmi_4x4_90hybrid/python/`
— `mmi_4x4_90hybrid_eme.py`, `mmi_4x4_90hybrid_fdtd.py`,
`mmi_4x4_90hybrid_hybrid.py`.

**6a. EME script.** Full-device EME solve across the O-band, both driven
inputs — the cheapest of the three, and the first place A3's design (does
the self-imaging + sweep result actually hit the cited CMRR/phase-error
requirements) gets checked end-to-end rather than on a bare uniform-body
cut.

**6b. FDTD script (full single-domain reference, for check 4).** Flag the
cost concern up front, same as the modal multiplexer's own doc: a 4-port,
double-tapered device is a long box, and 3-D FDTD cost scales steeply with
domain size. The modal multiplexer's own `CompoundGeometry.polygons()`
stub-per-port fix (`notes.md` Geometry item 8 part 2, landed 2026-08-04 —
"`_end_ports()` already answers exactly this 'which port's width wins'
question … by *not* picking a single winner: it returns every port at that
end") is **directly relevant and already landed** — this device needs the
exact same multi-port stub behaviour at both ends that fix enables, so this
script is not blocked on it the way the mux's once was. Budget a smoke test
at reduced resolution first (this repo's own "two-step discipline for
expensive, never-run scripts" — every prior full-device FDTD script in this
directory did this before its real run).

**6c. Hybrid script.** Input/body/output split per A5/A6: `INPUT_METHOD =
OUTPUT_METHOD = "fdtd"`, `BODY_METHOD = "eme"` (this repo's standard split
for every MMI hybrid so far), `use_neighbor_aperture=True` on the geometry
(not `half_domain_aperture` on the solvers — A6), 4 FDTD input-arm solves +
1 EME body length-sweep, output section declared but never separately
simulated (mirroring `mmi_2x2_hybrid.py`'s own `sub_region_configs`
convention for documenting *which* sub-region plots correspond to a real
simulation vs. a mirrored reuse).

**⚠ In progress 2026-08-06/07 — 6c only, at the user's direction (Step 4's
"proceed now, refine later" decision).** `examples/mmi_4x4_90hybrid/python/
mmi_4x4_90hybrid_hybrid.py` written: Step 1-3's settled parameters
(`waveguide_width=0.35`, `mm_width=6.0`, `waveguide_offset=1.5`,
`mm_length=80.1896` — the ΔL=-1.50 point), `N_MODES=8`, A11's phase-sign
fix applied at the same call site with the diagnostic trail recorded
in-line, and a headline CMRR/phase-error-vs-λ plot (Step 7's own
deliverable) using the new `SMatrix` methods for real, not synthetic, data.

**First submission (job 11716285) TIMED OUT at the full 4h budget** — not
a bug, a genuinely new cache fingerprint. The script builds its own `mmi`
object with `mm_length=80.1896` (the *settled* design point, matching
`mmi_2x2_hybrid.py`'s own convention of the geometry object representing
the reported point, e.g. its `mm_length=28.0`), which shifts the input/
output section boundaries (`x0`/`x1`, part of `solve_section()`'s cache
key) by 0.75 µm relative to Step 3's own exploration script (which used
the un-corrected `mm_length=81.6896` as its geometry reference throughout
the sweep) — a real, different fingerprint, not a mismatch to fix. All 4
`input_idx=k solved` lines DID print before the timeout, meaning all 4 real
FDTD input-section solves completed and were cached (cache writes are
atomic — `beef._utils.atomic_replace_write`) before time ran out; only the
(cheap, EME-only) body sweep never started. **Resubmitted unchanged as job
11723505**, expecting cache hits on all 4 inputs this time and a fast
finish.

**✅ Job 11723505 COMPLETED 2026-08-07 (1:13:00)** — all 4 inputs hit
cache as predicted (near-instant), the body sweep (reused Step 3's exact
81-point grid, sliced at ΔL=-1.50) was the only real new work. Full,
per-port results, not just the earlier scalar summary:

| | 1.26 µm | 1.31 µm (centre) | 1.36 µm |
|---|---|---|---|
| CMRR_I, port 1 | 15.68 dBe | 33.32 dBe | 18.64 dBe |
| CMRR_Q, port 1 | 12.94 dBe | 20.61 dBe | 12.65 dBe |
| CMRR_I, port 3 | 12.94 dBe | 20.61 dBe | 12.65 dBe |
| CMRR_Q, port 3 | 12.45 dBe | 29.46 dBe | 17.66 dBe |
| Δφ_6 | 7.34° | 4.75° | 8.28° |
| Δφ_7 | 2.22° | 0.14° | −0.63° |
| Δφ_8 | −1.52° | −5.03° | −1.20° |
| imbalance, port 1 | 1.99 dB | 0.81 dB | 2.06 dB |
| imbalance, port 3 | 4.90 dB | 0.91 dB | 4.89 dB |
| excess loss, port 1 | 3.09 dB | 2.06 dB | 2.66 dB |
| excess loss, port 3 | 4.27 dB | 2.21 dB | 3.98 dB |

**A more granular reading than the earlier scalar "max|Δφ|=8.28°"
summary: only Δφ_6 actually breaches the cited 5° bar** (up to 8.28° at
the band edges) — Δφ_7 stays under 2.2° throughout, Δφ_8 stays within
~5.0° (touching the boundary at centre, −5.03°). CMRR clears 20 dBe
comfortably at centre (up to 33.3 dBe) but falls to 12.5–19 dBe at both
edges — the same wavelength-dependent pattern the phase error shows,
consistent with a self-imaging condition that is well-tuned near centre
and drifts at the band edges, exactly as expected for a first pass with
no dispersion-aware refinement.

**Excess loss (2.1–4.3 dB) is real and worth flagging even though it is
not one of the paper's own cited pass/fail requirements** — the paper's
own device achieved <0.68 dB, but that number belongs to a completely
different platform/geometry (A9) and is not a fair comparison point; still,
2-4 dB is a meaningful insertion-loss cost for a real receiver and a
natural target for the same taper-refinement pass that would address the
phase gap.

Results plot: `example_output/mmi_4x4_90hybrid/mmi_4x4_90hybrid_hybrid/
hybrid_results.png`. SMatrix saved for both driven ports. **Step 6c is
done.**

**6a and 6b started 2026-08-07, at the user's direction ("let's complete
Step 6").**

**6a (EME-only).** `mmi_4x4_90hybrid_eme.py` written: 8 independent full
5-region-chain `EMESolver.solve()` calls (2 driven inputs x 4 outputs),
matching `mmi_2x2_eme.py`'s own per-port-pair pattern generalized to 4
ports — flagged in the script's own docstring that this device's port gap
(0.7 µm, Step 2) is narrower than the 2×2's own 1.05 µm, so the "ports
comfortably decoupled" approximation every EME-only script in this repo
makes is less obviously safe here, not independently checked. Smoke-tested
first: one `(input_idx=0, output_idx=0)` call at `RES=15`/
`N_MODES_TAPER=2`/`N_MODES_BODY=4`/1 wavelength took 29.4 s on the login
node — much cheaper than the hybrid script's own FDTD input-section
solves, since `EMESolver.solve()` has no time-stepping. Extrapolated via
res² (2-D eigensolve) × ~2× (mode count) × 3× (wavelengths, no per-
wavelength reuse in `solve()`'s own loop) to ~15–25 min/call at production
settings, ~2–4 h for all 8 serially. Real run submitted as **job
11726819** (`-t 360`).

**✅ Job 11726819 COMPLETED 2026-08-07 (0:40:46)** — faster than the
conservative estimate. All 8 calls succeeded:

| | 1.26 µm | 1.31 µm | 1.36 µm |
|---|---|---|---|
| CMRR_I, port 1 | 12.15 dBe | 37.15 dBe | 33.16 dBe |
| CMRR_Q, port 1 | 8.59 dBe | 20.61 dBe | 15.23 dBe |
| CMRR_I, port 3 | 8.59 dBe | 20.61 dBe | 15.23 dBe |
| CMRR_Q, port 3 | 9.04 dBe | 41.96 dBe | 35.30 dBe |
| Δφ_6 | 26.82° | 3.18° | 7.19° |
| Δφ_7 | 4.10° | −0.10° | 1.44° |
| Δφ_8 | −5.69° | −3.95° | −4.15° |
| imbalance, port 1 | 3.39 dB | 0.81 dB | 1.52 dB |
| imbalance, port 3 | 8.00 dB | 1.01 dB | 2.32 dB |

**A real, meaningful divergence from 6c's own FDTD-hybrid results — not
noise, a genuine cross-method finding.** At 1.26 µm, Δφ_6 is **26.82° in
the EME script vs. 7.34° in the hybrid script** — a large gap. CMRR values
also disagree by several dB at multiple points (e.g. port 1 CMRR_I at
1.36 µm: 33.16 dBe EME vs. 18.64 dBe hybrid). This is consistent with the
EME script's own flagged, unverified approximation: each of its 8 calls
draws only ONE input port and ONE output port at a time, treating the
other 6 as physically absent — a decoupling assumption this repo's own
2×2 MMI example states outright ("comfortably decoupled" at its 1.05 µm
gap), but this device's gap is narrower (0.7 µm, Step 2) and was
explicitly flagged in the EME script's own docstring as *not*
independently checked. **The size of this divergence is itself evidence
the approximation is not free here** — the hybrid script's real FDTD
input sections, which include all 4 ports' actual geometry simultaneously,
are the more physically trustworthy of the two for anything beyond a
quick sanity check. Not yet formally isolated (would need e.g. a direct
crosstalk measurement between adjacent ports, mirroring `mmi_1x4_hybrid.py`'s
own T_cross=0.0002 check for its output ports) — flagged as a real,
open finding rather than either dismissed or over-interpreted.

**6b (full-device FDTD).** Hit a real, genuine library gap, not a device
design issue: `FDTDSolver.solve()`'s `port` parameter only ever accepted
`"top"`/`"bot"` (`geometry.input_ys`' *last*/*first* entry) — fine for
every existing 2-port MMI/coupler example, but this device needs the
paper's own ports 1 **and** 3 (finding A4), and port 3 (`input_ys[2]` of
4) is an *interior* port neither string can reach. Fixed properly, not
worked around: `FDTDSolver.solve()` gained an explicit `input_idx`
parameter that overrides `port` when given, plus a `TypeError` guard
against passing it both as the new named argument and inside `**kwargs`
(where it would silently size the source aperture for the wrong port).

**An incidental, real bug caught and fixed in the same change:** before
this, when a caller didn't pass `input_idx` via `**kwargs` at all —
*every* existing example script, including this repo's own
`mmi_2x2_fdtd.py` — the source **aperture width** was computed using
`geometry.source_aperture_at()`'s own default `input_idx=0`, regardless
of which port `port="top"`/`"bot"` actually excited. Invisible for the
2×2 MMI specifically only because its aperture is symmetric between its
two ports (both give the same width by construction) — a real latent
mismatch for any asymmetric device, now closed for every caller, not just
this one. Verified via source inspection of all 5
`source_aperture_at()` implementations in the repo (`mmi.py`,
`coupler.py`, `base.py`, `compound.py`, `uniform.py`) — all accept
`input_idx` either explicitly or via a `**kwargs`/`**_kwargs` catch-all,
so passing it explicitly cannot break any of them. Fast suite re-run
clean (639 passed, 19 skipped) after the change. **Not yet independently
re-verified with a real meep run against a known-good existing device**
(the symmetric-aperture argument is why that's safe to defer, not a
substitute for eventually doing it — flagged, not silently skipped).

`mmi_4x4_90hybrid_fdtd.py` written: one output monitor per port
(`width=waveguide_offset=1.5 µm`, touching neighbor midpoints exactly),
`port="bot"` for port 1 (already reachable) and the new `input_idx=2` for
port 3. Smoke-tested (`RES=10`, `DECAY_TOL=1e-6`, `N_WL=1`, `srun
--mpi=pmi2`, job **11726835**) before any real cluster time, per this
repo's own two-step discipline — the first-ever full single-domain FDTD
run on this geometry AND the first real exercise of the new `input_idx`
path.

**Real production run submitted 2026-08-07 as array job 11733865**
(`-t 360`, `-N 2 --ntasks-per-node=24` per task, `--array=0-1`, one task
per driven input port: task 0 → `port="bot"` (port 1), task 1 →
`input_idx=2` (port 3)). Two lower-resolution calibration attempts first
(job 11726835, `RES=10`/`decay_tol=1e-6`: ~0 transmission; job 11726854,
`RES=20`/`decay_tol=1e-8`, 24 ranks, 6:44: **still ~1e-17 transmission —
essentially the noise floor**, not a smoothly improving trend to
extrapolate from) showed this device's 0.35 µm waveguide genuinely does
not guide at coarse resolution over an 84 µm propagation length — a real
threshold effect, not evidence the full RES=40 run will also fail (the
hybrid script's own short input-section smoke test showed real coupling
already at `RES=15`, confirming this is about propagation *length*
amplifying under-resolution loss over 84 µm, not this cross-section being
fundamentally unguidable). Rather than keep calibrating at resolutions
that aren't informative, sized directly from this repo's own best
precedent for a comparable-length device: **the directional coupler's own
~78 µm device took ~3 h on 48 ranks at full production settings** — this
device is 84 µm, closely comparable. Budgeted at `-t 360` (2× that
precedent) per task.

**⚠️ Both array tasks TIMED OUT 2026-08-07/08 at the full `-t 360` (6h)
budget** (task 0/port 1: 06:00:10 elapsed; task 1/port 3: 06:00:29) — the
directional-coupler precedent underestimated this device's real cost by
at least 2×. Not a crash or a bug: both tasks were still stepping cleanly
at a steady, unremarkable ~0.29 s/step right up to the kill (step ~70120,
meep time ≈876), no NaNs, no slowdown, no error — genuinely still
converging toward `decay_tol=1e-11` when the clock ran out. `FDTDSolver`
has no checkpoint/restart path (checked directly — no `checkpoint`/
`restart`/`resume` in `beef/solvers/fdtd.py`), so any resubmission restarts
from t=0 and re-pays the full cost, not just the remaining delta. No
decay-percentage diagnostic is printed by this run to extrapolate
remaining time from, so the real distance to convergence is unknown —
could be a modest overrun or a much larger one. Partition (`cpu-small`)
has no meaningful wall-time ceiling for this account (21 days), so a
bigger `-t` budget is available if wanted; the open question is how much
more compute is worth spending chasing this cross-check, not whether the
cluster allows it. **Paused here, not resubmitted, pending the user's
call** — this repo's whole discipline (Step 3a's 16-mode timeout,
resolution calibration here already) has been to measure before spending
more, not to double the budget and hope.

**Decision 2026-08-08 (user's direct call): drop Step 6b for now.** Not
resubmitted at a bigger budget or a loosened `decay_tol` — the hybrid
script (Step 6c) stands as the trusted, already-completed result for this
device; the full-device FDTD cross-check (check 4, Step 8) remains
formally open, not closed, if picked up again later.

**Re-queued 2026-08-14 at the user's direction.** Two real issues fixed
before resubmitting, not just a retry at the same config:

1. **The script and notebook were stale.** `mmi_4x4_90hybrid_fdtd.py`/
   `.ipynb` still carried `ΔL=-1.50`/`mm_length=80.1896` — the point from
   *before* the 2026-08-09 `N_MODES_BODY=8→16` body-basis convergence fix
   moved the settled point to `ΔL=-2.75`/`mm_length=78.9396` (job
   11768847). Job 11733865's own timeout was therefore already against a
   superseded geometry, not just an incomplete one. Both files updated to
   the current point, with the same grid-snap `assert` `mmi_4x4_90hybrid_
   hybrid.py` uses (do not silently let a "settled" length drift off the
   Step 3 sweep grid).
2. **Budget raised, not just resubmitted at `-t 360`.** Job 11733865
   reached meep time≈876 in the full 6h with no sign of trouble (steady
   ~0.29 s/step, no NaNs) — genuinely still converging, not stuck. New
   budget: **`-t 4320` (3 days, 12× the previous attempt)**, `run_
   mmi_4x4_90hybrid_fdtd.batch` (updated in place, same array-of-2-tasks
   shape). Still a real, uncertain cost, stated plainly: `FDTDSolver` has
   no checkpoint/restart (confirmed directly), so if 3 days is not enough
   the entire budget is lost, not partially credited, and a fresh, larger
   resubmission would be the only path forward — this is not a confident
   estimate, it is the most defensible next step given the one real data
   point available (876 meep-time units reached in 6h, unknown distance
   to actual decay).

A cheap reduced-resolution smoke pass was NOT re-run first — this exact
device already showed (job 11726835/11726854, above) that reduced
resolution/decay_tol sits below the noise floor over this 84 µm
propagation length, so a smoke pass here would be uninformative by this
device's own established precedent, not skipped carelessly.

**Submitted 2026-08-14, job 11847078** (array 0-1, same task→port mapping
as before: task 0 → `port="bot"` (port 1), task 1 → `input_idx=2`
(port 3)). Task 0 (port 1) started running normally. **Task 1 (port 3)
stayed `PENDING` on `AssocGrpBillingMinutes`** — the `gts-sr108-ece`
account's concurrent usage allocation was saturated (this session's own
jobs plus a separate agent's concurrent work on other devices sharing the
same account). Cancelled just the still-pending task 1 (`scancel
11847078_1`, task 0 left running untouched) and **resubmitted it under
the plain `gts-sr108` account (a separate allocation this user also has)
as job 11859047** (`sbatch --account=gts-sr108 --array=1 batch/run_
mmi_4x4_90hybrid_fdtd.batch`) — confirmed picked up correctly
(`Account=gts-sr108`, reason changed from `AssocGrpBillingMinutes` to
plain `Priority`, i.e. now just a normal queue wait). **Port 1 = job
11847078 task 0; port 3 = job 11859047 task 1** — two different job IDs
for the two halves of this one check, going forward.

**❌ BOTH TIMED OUT 2026-08-17 at the full `-t 4320` (3-day) budget —
still no successful full-device FDTD reference for this device, on the
second real attempt.** Port 1: meep time≈11007 (step 880580). Port 3:
meep time≈10991 (step 879252) — both essentially the same distance
reached, both still stepping cleanly at a steady ~0.286 s/step, no NaNs,
no slowdown, genuinely still converging toward `decay_tol=1e-11` when
killed — the identical signature the ORIGINAL 6h attempt (job 11733865)
showed at meep time≈876, just ~12.5× further along this time for
~12× the budget (6h → 3 days). **Not resubmitted, per the user's
explicit direction** — holding here rather than escalating the budget
again without a decay-percentage diagnostic to bound how much further
would actually be needed. `FDTDSolver` still has no checkpoint/restart,
so both 3-day attempts are a full loss, not partial progress toward a
future resubmission. Check 4 (Step 8) remains formally open — this
device's full-device FDTD cross-check has now failed to complete at 6h,
then at 3 days, and whether it is tractable at all within a reasonable
budget is itself an open question, not just "needs more time.

**2026-08-21 — the "still converging" reading was never MEASURED, and
today's instability finding says it may be wrong. Diagnostic submitted,
job 12079104.**

Both 3-day runs were recorded above as "genuinely still converging toward
`decay_tol`." That was an *inference* from a steady ~0.286 s/step and no
NaNs — not a measurement. Checked directly today: `grep -c "decay("` on
both logs returns **0**. `stop_when_dft_decayed` prints no decay
diagnostic at all, so **those two runs burned 144 CPU-days each and
produced zero information about how far from convergence they were.**
That is why "how much more budget?" has never been answerable.

Meanwhile the body calibration on this same device found a real late-time
**instability** (see this doc's own correction section): mode 16's energy
bottomed just above `decay_tol` near t≈260 and then *grew* monotonically
to 1.7e53 by t≈2081, traced to the undamped (`gamma=0`) Lorentzian poles
every material in this repo is built from (`beef/materials/stack.py:284`).
**The 6b logs carry that same pole on their own face** —
`lorentzian susceptibility: frequency=0.639386, gamma=0`. With no energy
readout, a diverging run and a slowly-converging one are
**indistinguishable**: both step steadily, neither NaNs, neither stops.
So the recorded "still converging" and "no NaNs" observations do not
actually discriminate between the two.

*Counter-evidence, stated so this is not overclaimed*: both runs reached
t≈11000 without NaN/Inf. Unbounded exponential growth would plausibly
have overflowed float64 (~1e308) by then, which argues against divergence
as fast as mode 16's. The point is that nobody knows, and the previous
spend could not have found out.

**Job 12079104** (`run_mmi_4x4_90hybrid_fdtd_DECAYDIAG.batch`) settles it
cheaply: `decay_method="energy"` prints the ratio every 20 time units.
Deliberately a *measurement*, not an attempt to finish check 4 — **one
port** (port 1) not two, **`-t 720` (12 h)** not 3 days, but otherwise
**production physics** (res=40, `decay_tol=1e-11`, 3 λ, same ΔL=−2.75
geometry), since the instability is itself resolution-dependent and a
reduced-resolution probe would not answer the question. At the previously
measured ~153 meep-time-units/hour this reaches t≈1800 — well past where
mode 16's divergence was unmistakable (t≈580). Three outcomes, all
actionable, none of which the prior 6 h + 3 d + 3 d could produce:
  - **ratio falling steadily** → genuinely a budget problem, *and* the
    printed sequence finally gives a real extrapolation to completion;
  - **ratio plateaued or rising** → the instability; check 4 is not
    fixable by a bigger budget and that question closes for good;
  - **ratio oscillating above `decay_tol`** → mode 16's own near-miss
    case, i.e. wants a looser tol rather than more time.

`decay_method` is part of the cache fingerprint by design, so this run
neither hits nor pollutes the existing default-`"dft"` entries. The
user's standing no-auto-resubmit policy for Step 6b still applies to it.

#### ✅ RESULT — job 12079104, and it answers a question open since 2026-08-08

**TIMEOUT at 12 h, exactly as designed** (this was scoped as a
measurement, not an attempt to finish check 4). Reached **t=1801, ratio
1.441e-3, still falling monotonically, no divergence.** Outcome **(a)**:
Step 6b is a genuine budget problem — and for the first time there is a
real extrapolation instead of a guess.

**The decay curve, 70 ringdown checkpoints.** Energy peaked at 29.06
(t≈60-80), stayed flat while the pulse crossed the 84 µm device (one-way
transit ≈340 time units — a flat plateau there is correct for a lossless
structure, not a stall), then dumped ~93% into the PML at t≈340-420. What
remains is a **long-lived trapped component**, ~5.7% of peak at t=420,
decaying steadily since:

| window | k |
|---|---|
| t=420→821 | 3.24e-3 |
| t=821→1221 | 2.94e-3 |
| t=1221→1621 | 2.32e-3 |
| t=1621→1801 | 2.54e-3 |

k decelerates mildly then **stabilizes around 2.3-2.5e-3** — it does not
collapse toward a floor. One earlier feature is worth recording because it
nearly caused a wrong call in real time: the ratio **plateaued at ≈0.0566
from t≈520-660** (k dropped to 6.9e-5), which looked exactly like an
asymptote; decay then resumed and *accelerated*. Consistent with beating
between trapped modes. **Do not call a floor from a plateau on this
device** — this is the second time today that pattern appeared (mode 16's
own near-miss was the first).

**Extrapolation to `decay_tol=1e-11`:**

| basis | k | converges at | wall time |
|---|---|---|---|
| average since ringdown | 2.794e-3 | t≈8,524 | **56 h** |
| recent ~300 t-units (conservative) | 2.063e-3 | t≈10,907 | **72 h** |

**The finding that matters.** The two 3-day DFT runs reached **t≈11,007**
— *at or past where energy convergence lands on either estimate* — and
**still did not stop.** So the physics converges; what never converged is
`stop_when_dft_decayed`, which requires every DFT component at every
monitor to individually reach tol. That is the same over-strict criterion
that stalled mode 1 in the calibration. **This one stopping criterion has
now blocked this device three separate times** (mode 1, mode 16, Step 6b),
and in 6b's case it burned 288 CPU-days behind a diagnostic that printed
nothing at all.

**Recommended, NOT submitted** (standing no-auto-resubmit policy): re-run
Step 6b with `decay_method="energy"`, **`-N 4 --ntasks-per-node=24`
(96 ranks)** and **`-t 2880` (48 h)**. Both ports, as the original array did.

> ⚠ **SUPERSEDED 2026-08-24 — that re-run was tried (job 12227584, both
> ports, 96 ranks) and cancelled at t≈6,200 after 21.5 h, because the
> premise was wrong. `decay_tol=1e-11` is not merely expensive on this
> device, it is UNREACHABLE. See "decay_tol has a floor" below.**

#### `decay_tol` has a hard floor, and 1e-11 is below it (2026-08-24)

Prompted by the user asking whether `decay_tol=1e-11` is necessary at all,
and whether a shorter or lower-resolution run could settle it. A ladder of
tolerances was run at res=20 (jobs 12260633, 12262138-141) on the reasoning
that the question is RELATIVE -- at what tolerance do the S-parameters stop
moving -- so the repo's "reduced resolution sits below the noise floor"
caveat does not bind. The ladder answered a different and more important
question instead.

**Every rung below ~1e-4 diverges.** Absolute energy from the 1e-5 rung:

| t | energy | ratio |
|---|---|---|
| 260 | 29.58 | 1.00 (peak, source done) |
| 3,624 | 5.16e-3 | 1.74e-4 |
| **3,705** | — | **1.682e-4 ← FLOOR** |
| 5,787 | 51.9 | 1.00 (already growing) |
| 8,030 | 2.14e8 | 1.00 (runaway) |

The 1e-7, 1e-9 and 1e-11 rungs bottom at **exactly 1.682e-4 at exactly
t=3,705** -- bit-identical across three independent jobs, i.e. deterministic
numerics, not noise. Energy decays normally, floors, then grows
exponentially: the `gamma=0` undamped-Lorentzian instability
(`beef/materials/stack.py:284`). **Third independent observation** after
calibration mode 16 (bottomed 1.95e-11, grew to 1.7e53) and its isolated
re-run -- and the first on a completely different configuration (res=20,
full device, no calibration), which is what makes the root cause solid
rather than inferred.

**The consequence is not "expensive", it is "impossible".** `decay_tol`
below the instability floor can never be satisfied, because the field
starts growing before the energy gets there. **This explains all five 6b
failures at once** -- job 11733865 (6 h), 11847078 + 11859047 (3 days
each), and 12227584 (two ports, 48 h) -- none were under-budgeted; all were
chasing a target the simulation cannot reach. Every extrapolation in the
sections above, including my own 28-37 h and 56-72 h estimates, assumed a
single exponential continuing indefinitely and was therefore wrong in
principle, not just in magnitude.

**The floor is resolution-dependent**, which is itself the signature of a
numerical rather than physical effect: 1.68e-4 at res=20, versus res=40
reaching 5.27e-6 by t=6,224 with no turnaround yet observed. res=40's floor
is lower but unmeasured, and there is no reason to think it is below 1e-11.

**Open, and now the load-bearing question**: (a) measure res=40's actual
floor by running until it turns around, then set `decay_tol` above it with
margin; or (b) fix the instability at source by giving the Lorentzian poles
a small nonzero `gamma`, which would remove the floor and unblock BOTH 6b
and the mode-16 calibration. (b) touches shared material code used by every
device; (a) is cheap and device-local. Neither started.

#### The slow tail is PHYSICAL — `dpml` ruled out (job 12260446, 2026-08-24)

Two effects had been tangled together all day and needed separating. The
first was whether the long tail exists at all for a legitimate reason.

`dpml=0.5` µm is only ~0.4 wavelengths at 1.31 µm, and **`dpml` has never
been convergence-tested anywhere in this repo** (audited: the PML-adjacent
sweeps all vary domain *clearance*, and `dpml` is only ever raised 0.5→1.0
to preserve pixel count at low resolution). A too-thin PML reflects
low-angle and evanescent components, and reflected energy rattling between
domain walls mimics a trapped mode exactly. Worth testing before accepting
the tail as real.

Re-ran the identical configuration at **`dpml=1.5`** (cell 12.0 × 5.0 µm vs
the baseline's 10.0 × 3.0, ~2× the voxels) and compared at matched meep
time:

| t | dpml=1.5 | dpml=0.5 baseline | ratio |
|---|---|---|---|
| 440 | 6.533e-2 | 6.519e-2 | 1.002 |
| 540 | 5.673e-2 | 5.673e-2 | **1.000** |

**Agreement to four significant figures.** Tripling the PML changes the
trapped fraction by nothing measurable. **The ~5.7% trapped component is
real device physics — genuinely long-lived resonant modes in the MMI body —
not a boundary artifact.** `dpml=0.5` is vindicated for this device; the
untested-`dpml` gap is real but is not what causes slow decay here.

**This separates the two phenomena, which is why "more budget" failed five
times — it addressed neither:**

| | cause | consequence | fix |
|---|---|---|---|
| slow tail | **physical** trapped modes (confirmed above) | decay takes t~19,000+ to reach 1e-11 | pick a tolerance the physics justifies |
| floor + reversal | numerical (hypothesis: `gamma=0` poles) | tolerances below the floor are UNREACHABLE | damping, or a reachable target |

Job 12273630 (Courant halved, res=20, everything else identical) is the
cheap discriminator for the second: it changes only the timestep, so a floor
that moves down/later confirms a numerical instability and justifies editing
shared material code, while an unchanged floor would remove most of the
support for the `gamma=0` hypothesis. Deliberately run BEFORE touching
`beef/materials/stack.py`, which every device draws on.

#### ❌ RESULT — the `gamma=0` timestep-stability hypothesis is FALSIFIED (job 12273630)

| | floor ratio | floor time |
|---|---|---|
| Courant=0.5 (reference) | 1.6824e-04 | t=3,705 |
| Courant=0.25 | 1.6826e-04 | t=3,742 |
| change | **0.01%** | +1.0% |

**Halving the timestep did not move the floor.** The proposed mechanism —
undamped Lorentzian poles sitting on the imaginary axis, with time
discretization nudging them into exponential growth — predicts a
substantially lower and later floor at half the timestep. It did not happen,
so that mechanism is not what is occurring, and **the `gamma` fix is NOT
justified on this evidence.** A shared edit to `beef/materials/stack.py`,
which every device in the repo draws on, was one step away from being
recommended on a diagnosis that then failed its own test. The ~1 h
discriminator was worth it.

**Control value**: the two runs agreed to 0.1-0.5% at every matched time
BEFORE the floor (t=500/900/1300/1700), confirming that halving the timestep
did not perturb the physics — so the null result is a real null, not a
botched comparison.

**What is now established about the floor, all measured:**

| property | result | implication |
|---|---|---|
| resolution | 1.68e-4 @ res=20, ~5e-6 @ res=40 | **spatial** discretization matters |
| timestep | invariant (0.01%) | **not** temporal / not Courant |
| PML thickness | invariant (0.05%, job 12260446) | not a boundary artifact |
| repeatability | bit-identical across 3 jobs | deterministic, not noise |

Resolution-dependent but timestep-independent points to a **spatial**
discretization effect. One candidate worth recording, NOT a conclusion:
`resolved/gds_geometry_hybrid_cascade_validation.md` documents that meep's
subpixel averaging treats abutting same-material prisms as having an
internal interface, moving 41,677 epsilon voxels at res=40 — a spatial
discretization error of the right character and the right resolution
scaling. Untested against this floor.

**The practical answer is unchanged**: the floor is real, `decay_tol=1e-11`
remains unreachable, and picking a reachable tolerance (option 1, jobs
12272571/12272573) is still the path. Only the mechanism reopens.

---

## ✅ STEP 6b COMPLETED — first time ever, job 12272571 (2026-08-25)

**Port 1 finished: COMPLETED in 20:56:44, 107 curve points, reached
`decay_tol=1e-5` (final ratio 1.020e-05 at t=5,350).** Every previous attempt
failed — 11733865 (6 h), 11847078 + 11859047 (3 days each), 12227584 (two
ports, 48 h). None were under-budgeted; all were chasing `decay_tol=1e-11`,
which this device's own numerical floor makes unreachable.

**Final |T|² at the reachable tolerance** (rows = wavelength, cols = out0..3):

| λ | out0 | out1 | out2 | out3 |
|---|---|---|---|---|
| 1.26 µm | 0.088976 | 0.098147 | 0.064560 | 0.042439 |
| 1.31 µm | 0.141433 | 0.163326 | 0.185942 | 0.200832 |
| 1.36 µm | 0.230228 | 0.187837 | 0.062458 | 0.131395 |

### What `decay_tol` this device actually needs — measured against the TRUE endpoint

The run snapshotted its S-matrix every 50 time units (`solve(snapshot_every=)`,
added the same day), each sample tagged with the energy ratio that would have
stopped the run there. So the cost of any tolerance is read directly off ONE
simulation rather than inferred from a ladder of runs:

| `decay_tol` | t | max ΔT | ≈dB | verdict |
|---|---|---|---|---|
| 1e-2 | 1,000 | 2.078e-2 | 0.564 | too loose |
| 5e-3 | 1,300 | 4.859e-3 | 0.138 | too loose |
| 1e-3 | 2,000 | 2.466e-3 | 0.071 | marginal |
| 5e-4 | 2,450 | 1.837e-3 | 0.053 | marginal |
| **1e-4** | **3,500** | **5.334e-4** | **0.015** | **recommended** |
| 5e-5 | 4,050 | 3.153e-4 | 0.009 | OK |

**`decay_tol=1e-4` is the setting for this device**: 0.015 dB against the
0.1 dB the conclusions rest on, at t=3,500 ≈ **11.4 h** on 96 ranks — versus
the 63.6 h that `1e-11` would have needed had it been reachable at all.
Provisional until port 3 (job 12272573) confirms independently.

**The worst wavelength is 1.36 µm at EVERY tolerance — never band centre.**
This matters methodologically: an earlier reading of this same data at 1.31 µm
alone produced a confident "converged by ratio 5.7e-2" claim that was wrong by
two orders of magnitude, because band centre converges fastest while the band
EDGES are both slowest AND where this device's CMRR/phase failures concentrate.
Always evaluate convergence at the worst wavelength.

**Two retractions this measurement settles**, recorded because both were stated
confidently before being checked:
- "10⁻²–10⁻³ is ample" — WRONG. 1e-2 costs 0.564 dB, ~6x the budget.
- "converged by ratio 5.7e-2, 0.002 dB" — WRONG, an artifact of comparing
  against a moving reference (the then-latest sample at t=800, with 4,500 time
  units still to come) and of reading band centre only.

### ⚠ CHECK 4 — FIRST EVER RESULT, and the two methods DISAGREE (port 1, preliminary)

The full-device FDTD S-matrix above compared against the hybrid result (job
11778128, 2026-08-10, same corrected design point, **identical 3-wavelength
grid so no interpolation is involved**). Port labels map in order
(`out5..out8` = `out_0..out_3`).

| port | hybrid 1.26 / 1.31 / 1.36 | FDTD 1.26 / 1.31 / 1.36 | max abs | max rel |
|---|---|---|---|---|
| out5 | 0.081456 / 0.131037 / 0.233193 | 0.088976 / 0.141433 / 0.230228 | 0.0104 | 9.2% |
| out6 | 0.074165 / 0.154555 / 0.194348 | 0.098147 / 0.163326 / 0.187837 | 0.0240 | 32.3% |
| out7 | 0.039811 / 0.174411 / 0.060343 | 0.064560 / 0.185942 / 0.062458 | 0.0247 | 62.2% |
| out8 | 0.030496 / 0.204279 / 0.125779 | 0.042439 / 0.200832 / 0.131395 | 0.0119 | 39.2% |

**Totals per wavelength** — the disagreement is strongly wavelength-dependent:

| | 1.26 µm | 1.31 µm | 1.36 µm |
|---|---|---|---|
| hybrid | 0.225929 | 0.664282 | 0.613663 |
| full-device FDTD | 0.294121 | 0.691534 | 0.611918 |
| difference | **+30%** | +4% | −0.3% |

Agreement is essentially exact at 1.36 µm and degrades toward the SHORT
wavelength edge, where the FDTD finds 30% more transmitted power.

**This is not convergence error.** The FDTD side carries ~0.015 dB
(5.3e-4 in |T|²) of tolerance error at `decay_tol=1e-5`; the disagreements are
~0.024 in |T|², roughly 50x larger. The large *relative* figures on out7/out8
are partly small-denominator artifacts at 1.26 µm, but out6/out7's absolute
0.024 is real.

**Candidate causes, none tested:**
1. **The hybrid's port-isolation approximation.** Each `solve()` draws only ONE
   input and ONE output branch, treating the other 6 ports as physically absent
   — no crosstalk modelled at all. This device's port clearance is 0.7 µm
   (offset=1.5, taper_width=0.8), and `mmi_4x4_90hybrid_eme.py`'s own docstring
   already flags this as a stronger approximation than the 2×2's 1.05 µm gap.
   Crosstalk growing toward the band edge would match the observed pattern.
2. **Interface-basis truncation at N_MODES=16.** The hybrid contracts an FDTD
   interface basis against an EME body basis; the full-device FDTD has no such
   contraction. The extended convergence sweep (job 11936081) already showed the
   basis is not fully flat even at N=16-36.
3. **Subpixel averaging** (`resolved/gds_geometry_hybrid_cascade_validation.md`)
   — affects both methods but not necessarily equally.

### ✅ CONFIRMED BY BOTH PORTS — job 12272573 COMPLETED 1-03:39:25, 127 points

Port 3 finished (final ratio 1.03e-05 at t=6,350) and reproduces the pattern
independently:

| λ | port 1 driven | port 3 driven |
|---|---|---|
| 1.26 µm | **+30.2%** | **+49.6%** |
| 1.31 µm | +4.1% | +3.9% |
| 1.36 µm | −0.3% | −0.3% |

Two independent simulations, different driven ports, agreeing to 0.2% at
1.36 and 1.31 µm while both showing a large POSITIVE discrepancy at 1.26 µm.
Max |T|² difference over all ports/wavelengths/inputs: **0.0255**, versus a
convergence error of 0.0005 at `decay_tol=1e-5` — **48×**. This is a real
systematic effect, not noise and not under-convergence.

**The sign is diagnostic.** The full-device FDTD finds MORE transmitted power
than the hybrid, and the gap grows toward SHORT wavelength. That is precisely
what the hybrid's port-isolation approximation predicts: each `solve()` draws
only ONE input and ONE output branch and treats the other six ports as
physically absent, so power coupling into an undrawn port is never collected.
As λ shortens the beat length shifts and self-imaging degrades, sending more
power to the "wrong" outputs — lost by the hybrid, captured by the full-device
FDTD. At 1.36 µm self-imaging is well formed and the methods agree essentially
exactly.

**Implication, stated carefully.** The device's headline failures (CMRR, Δφ)
are worst at the band edges, and the hybrid-vs-FDTD disagreement is ALSO worst
at a band edge, in the same direction. Every design conclusion for this device
— the `waveguide_offset` sweep, the `mm_width` scoping, the CMRR/phase
assessment — was computed through the hybrid method. If the hybrid
systematically under-collects at short wavelength, those conclusions inherit
that error. CMRR and Δφ are ratios, so a UNIFORM loss would cancel; this loss
is wavelength- and port-dependent, so it does not obviously cancel.

**NOT established**: that the port-isolation approximation IS the cause
(three candidates were listed; only the sign and wavelength-trend favour it).

### ✅ CHECK 4 CLOSED — the device metrics, both methods (2026-08-25)

`verification/python/mmi_4x4_90hybrid_check4_metrics.py` recomputes the paper's
metrics from BOTH S-matrices using the same `SMatrix` methods and the same
construction `mmi_4x4_90hybrid_hybrid.py` uses (including A11's corrected
self/other order on `cross_excitation_phase_diff_deg`). No new simulation —
the FDTD side comes back from cache via `cache_only=True`, the hybrid side
from its published pickles.

| | hybrid | full-device FDTD | spec |
|---|---|---|---|
| max \|Δφ\| | 12.32° | **20.46°** | ≤ 5° |
| min CMRR | 1.53 dBe | 5.96 dBe | ≥ 20 dBe |
| max excess loss | 8.37 dB | 6.62 dB | — |
| max imbalance | 14.62 dB | 8.03 dB | — |

**THE DEVICE'S VERDICT IS CONFIRMED, NOT OVERTURNED.** Both methods fail both
cited specs by large margins — phase error 2.5-4x over budget, CMRR 3-13x
under. The sectioning error does not rescue the device. **For the first time
the design conclusion rests on two independent methods.**

**But the raw-|T|² disagreement DID carry into the metrics — it did not
cancel**, which was the open question (ratios cancel a *uniform* loss; this one
is wavelength- and port-dependent). The damage is concentrated at 1.26 µm:

| metric @ 1.26 µm | hybrid | FDTD | diff |
|---|---|---|---|
| dphi_7 | 10.45° | 20.46° | **+10.0°** |
| imbalance_port3 | 14.62 dB | 8.03 dB | **−6.6 dB** |
| cmrr_Q_port3 | 1.53 dBe | 7.07 dBe | +5.5 dBe |

At 1.31 and 1.36 µm the two methods agree within ~1-2 units on almost every
metric.

**The hybrid was OPTIMISTIC on the headline metric.** Full-device FDTD puts
max|Δφ| at 20.46° against the hybrid's 12.32° — the device is FURTHER from
spec than the sectioned method reported, not closer.

**Do the earlier design conclusions survive? Argued, not re-run.** The
`waveguide_offset` sweep found best-achievable max|Δφ| of 125.6/64.8/47.6/75.1°
at alternative offsets against 12.3° at baseline, and the `mm_width` scoping
was similarly lopsided. Those gaps are an order of magnitude larger than this
~10° systematic discrepancy, so a band-edge error of this size cannot reverse
them. That is an argument from magnitude, not a re-run through the FDTD path —
which would be prohibitive (each offset point is a full-device FDTD run).

#### MPI scaling — measured 2026-08-22, jobs 12217184 / 12217185

Prompted by the user asking whether the re-run could be parallelized. It
can, and **no 6b attempt has ever tried**: the 6 h attempt, both 3-day
runs and the 12 h diagnostic all used the identical `-N 2
--ntasks-per-node=24` = 48 ranks. The domain is 82.95 × 10 × 3 µm at
res=40 ≈ **159 M voxels**, i.e. 3.3 M voxels/rank at 48 — far above the
~1e5-1e6 where meep's scaling normally rolls off.

| ranks | nodes | s/step | speedup | efficiency |
|---|---|---|---|---|
| 48 | 2 | 0.287 | 1.00× | — |
| **96** | **4** | **0.1469** | **1.95×** | **98%** |
| 192 | 8 | never reached time-stepping in 30 min | — | — |

**96 ranks is the operating point.** 192 did not merely scale poorly — it
never printed `set_epsilon` (21.8 s at 48 ranks) in a full 30 minutes, so
it hung in collective-heavy structure setup. SLURM had scattered those 8
nodes over four racks (`02-002`, `02-003`, `02-004`, `03-004`); inter-rack
latency is the likely cause. Retrying 192 with `--switches=1` and a longer
budget might work, but 96 already halves the runtime at 98% efficiency and
is not worth risking a multi-day job on.

**Only one axis was actually available**, checked rather than assumed:
- *ports* — already parallel, the original 6b was a 2-task array;
- *wavelengths* — already optimal, and this was VERIFIED, not reasoned:
  the 12 h diagnostic log contains exactly **one** `Initializing
  structure` and one time-stepping loop for all 3 wavelengths (single
  broadband pulse + DFT monitors). Splitting λ would triple the work;
- *time* — impossible, FDTD is sequential in time and `FDTDSolver` has no
  checkpoint/restart (which is also why all three timeouts were total
  losses);
- *MPI ranks* — the one real lever, and it was sitting unused.

**Revised estimate at 96 ranks: 29 h (average-k) to 37 h (conservative
recent-k)**, from 56-72 h at 48. `-t 2880` covers the conservative case
with ~30% margin.

**✅ CONFIRMED IN PRODUCTION — job 12227584.** The probe number transferred
essentially exactly: port 1 steady **0.1468 s/step**, port 3 **0.1466**,
against the probe's 0.1469 and the 48-rank baseline's 0.2865 — a measured
**1.95×** on the real both-ports run, not just on a 30-minute probe.
Completion estimate at that rate: **27.8 h** (average-k target t≈8,524 =
681,920 steps) to **35.6 h** (conservative-k t≈10,907 = 872,560 steps),
inside the 48 h budget with ~12 h margin in the worst case.

**Billing note.** The first submission (job 12218117) was blocked
`AssocGrpBillingMinutes` on `gts-sr108` — the fallback account this
project has leaned on all session — after the day's probe work drew it
down. Checked all three rather than guessing (see
[[slurm-account-gts-sr108-ece]]): `gts-sr108-ece` billing cap 354,561,109
and `gts-sr108` 174,062,774 both exhausted, while
**`gts-sr108-gedc_main` has no `GrpTRESMins` cap set at all** and was
already running work. Resubmitted there as 12227584 (precedent: the 1×4
sin300 doc's own 64 h job used the same account). Worth remembering that
the cheap diagnostic probes are what exhausted the fallback — individually
minor, collectively enough to block a headline job.

**Caveat to settle before trusting the S-matrix it produces**: energy-decay
and DFT-decay are *different* convergence criteria. Stopping when bulk
energy is 1e-11 of peak means the residual field is negligible, so the DFT
monitors should have accumulated essentially the full result — that is the
argument, and it is the one already implicitly accepted for every
calibration run on this device, but it is an assumption to verify against a
converged reference, not to assert.


### EME script staleness — Step 6a was never re-run after the body-mode redo (found 2026-08-21)

**Prompted by the user asking a plain question: "are the EME field plots
up-to-date?"** They are not, and the cause was not stale images — it was
two stale *scripts*.

`mmi_4x4_90hybrid_eme.py` (mtime 2026-08-08 10:44, untouched since) and
`verification/python/plot_mmi4x4_90hybrid_eme_fields.py` both still
carried `DELTA_L=-1.50`, `MM_LENGTH=80.1896`, `N_MODES_BODY=8` — the point
superseded on 2026-08-09 when the body basis went 8→16 and moved the
optimum to `ΔL=-2.75`/`mm_length=78.9396`. `mmi_4x4_90hybrid_hybrid.py`
(mtime 2026-08-10 02:13) and `body_sweep.py` did get the change.

**Every EME artifact on disk therefore shows a rejected design point at a
50-64% under-converged basis:**

| artifact | mtime | reflects |
|---|---|---|
| `mmi_4x4_90hybrid_eme/geometry.png` | 2026-08-08 10:45 | ΔL=−1.50 geometry |
| `mmi_4x4_90hybrid_eme/eme_results.png` | 2026-08-08 13:11 | N_MODES_BODY=8 |
| `fields/eme_{Ex,Ey,Ez,Hx,Hy,Hz}_out5_in0_1.310um.png` | 2026-08-08 17:01 | N_MODES_BODY=8 |
| (hybrid's own outputs, for contrast) | 2026-08-10 02:26–12:36 | correctly regenerated |

**How it hid for 12 days**: the doc's own CURRENT STATUS header asserted
"Step 6 (all 3 example scripts) done and RE-RUN at a properly-converged
body mode basis," while the line immediately below it scoped the fix to
"(body_sweep.py/hybrid.py)". The summary and its own supporting detail
contradicted each other on adjacent lines. Nothing anywhere recorded
leaving the EME script behind deliberately, so this was an oversight, not
a decision. Header now corrected in place.

**Also a real contributing cause, worth fixing properly later**: the
field-plot script *duplicates* the design-point constants rather than
importing them from `mmi_4x4_90hybrid_eme.py`, which is exactly how the
two drifted apart silently. Noted in both files; not refactored today.

**Corrected 2026-08-21**: both scripts moved to `DELTA_L=-2.75` /
`N_MODES_BODY=16`, each with a comment recording why. All EME artifacts
need regeneration, plus `notebooks/mmi_4x4_90hybrid_eme.ipynb` (which this
repo generates mechanically from the `.py`, never hand-copied).

**✅ FIELD PLOTS REGENERATED — job 12212996, 8m26s, 2026-08-22.** All six
`eme_{Ex,Ey,Ez,Hx,Hy,Hz}_out5_in0_1.310um.png` now at the settled
`ΔL=-2.75` / `N_MODES_BODY=16` / `N_MODES_TAPER=12` point. Verified by
LOOKING at the output, not by exit code: source at x≈-2.2 µm (port 1),
beam diverging into the body, rich self-imaging interference consistent
with a 16-mode basis, and total length ≈81 µm matching the corrected
`total_length=80.9396` (the superseded ΔL=-1.50 geometry was 82.19 µm).
Prior plots preserved in `fields_backup_nbody8_dL150_nmt12_20260822_084608/`.

**It took two attempts, and the first one is the more instructive.** Job
12212950 ran to COMPLETED in 4m33s and produced six plausible-looking
plots — with **`N_MODES_TAPER=4`**. That is a THIRD stale constant in
`plot_mmi4x4_90hybrid_eme_fields.py`, missed on 2026-08-21 when `DELTA_L`
and `N_MODES_BODY` were corrected there. The 4→12 bump landed in
`mmi_4x4_90hybrid_eme.py` on 2026-08-08 (job 11745928 rendered the
comparison; at 4 the taper's 0.8 µm end is missing 2 of its 6 real guided
modes) and this file never received it. **Caught only because the script
prints its own config and the log said `N_MODES_TAPER=4`** — and because
the pre-existing `fields_backup_nmt4_20260808_102153/` directory is named
for precisely that value, i.e. the run had silently regressed to the plots
the 2026-08-08 bump existed to supersede. A green exit code and a
reasonable-looking image were both worthless as verification here.

**Root cause is now three-for-three, and the obvious fix does NOT work.**
`plot_mmi4x4_90hybrid_eme_fields.py` duplicates `mmi_4x4_90hybrid_eme.py`'s
design-point constants (`DELTA_L`, `MM_LENGTH`, `N_MODES_BODY`,
`N_MODES_TAPER`) instead of sharing them, and all four drifted. The
tempting fix — have the plot script `import` them from `eme.py` — is
**unsafe**: `eme.py` has no `def main()` and no `__main__` guard, and
carries 6 module-level statements that call `solve()`/`savefig()`, so
importing it would execute the entire >96 h 8-call sweep at import time
(checked via `ast`, not assumed). A correct fix needs a shared
constants/params module both scripts import, which means editing a shipped
example whose notebook is generated mechanically from it — a real refactor
deserving its own decision, not a drive-by. **Left as a documented wart** in the sense
that no shared module was created -- but the drift is now CHECKED rather
than merely noted; see below.

#### The drift is now gated, and gating it found a stale published notebook (2026-08-22)

The obvious fix -- one params module every script imports -- was scoped and
**rejected on evidence**, not on effort:
- `mmi_4x4_90hybrid_eme.py` has no `main()` and no `__main__` guard and runs
  8 EME `solve()` calls at module level (>96 h at the current basis), so
  nothing can import it and making it importable means restructuring a
  shipped example (verified via `ast`, not assumed);
- the three notebooks carry their **own** copies of these constants, so a
  `.py`-only refactor leaves them drifted and forces regenerating all three.

**What already existed, with a gap.** `scripts/check_notebook_sync.py` is
built for exactly this, and its own docstring records the same failure in
July -- *"2026-07-31: four mismatched constants, one of them reproducing a
known-bad mode basis in a published tutorial."* It missed this drift only
because `DELTA_L`/`N_MODES_BODY`/`N_MODES_TAPER` were absent from its
`CONSTANTS` list. **Added.** Unlike the 2026-08-18 addition (which was pure
coverage, 0 drift found), this one immediately surfaced a real one:

> **`mmi_4x4_90hybrid_eme.ipynb` was still at `DELTA_L=-1.50` /
> `N_MODES_BODY=8`.** The `.py` was corrected 2026-08-21 but the hand-mirrored
> notebook was not — a *published tutorial* teaching the rejected design point
> at the 50-64% under-converged basis. Exactly the failure that checker exists
> to prevent, recurring because its constant list did not reach this device.

Fixed, along with a stale derived comment (`# = 80.1896` → `78.9396`).

**Plus `tests/test_mmi_4x4_90hybrid_design_point.py`**, covering the axis the
notebook checker does not: `.py`↔`.py` agreement across all 9 scripts that
hardcode this device's design point. `mode_convergence.py` and
`taper_sweep.py` are explicitly EXEMPT (their recorded results genuinely are
at `ΔL=-1.50`; rewriting them would misrepresent what ran) and pinned as such,
so re-using either -- `taper_sweep.py` being the natural base for a
`taper_width` sweep -- forces a conscious decision instead of silently
inheriting stale geometry.

**Two bugs in that test, both caught by running it, both worth recording**
because each produced a confident false report of drift on files that were
correct: (1) a line-anchored regex over raw text matched this device's own
*prose* (`"At the converged N_MODES_BODY=16 design point"`,
`"N_MODES_BODY=4 proxy"`); (2) after switching to AST, regexing the
assignment's source segment pulled the `4` out of the env-var **name**
`MMI4X4_EME_NMODES_TAPER` instead of its default `12`. Now reads real AST
literal nodes. A checker that cries wolf is worse than none, and both
failures looked exactly like genuine drift.

**Verified**: `check_notebook_sync.py` clean for this device (18 pairs), new
test 11 passed, full suite **777 passed / 21 skipped**.

**Pre-existing drift deliberately NOT touched**: `directional_coupler_hybrid`
and `mmi_2x2_hybrid` notebooks are both missing a `calibrate_body_neff()`
call. Different devices, concurrent work by another agent -- reported, not
fixed.

**Budget being MEASURED, not extrapolated — probe job 12130641.** Job
11736983 did all 8 `solve()` calls in 1:13:27 at `N_MODES_BODY=8`;
doubling the basis scales by no factor this repo has measured, and this
device has already been burned twice by guessed EME costs (the
`num_cells` 0.3-0.6 dB extrapolation that measured 0.08-0.10 dB, and the
isolated-`solve()` 8× reading that was mostly one-time JAX JIT).
`mmi_4x4_90hybrid_eme_timing_probe.py` therefore times **two** calls —
call 1 carries the JIT, call 2 is steady state, estimate = `call1 +
7*call2`. Two deliberate safety properties: it writes **nothing** to the
EME script's hardcoded `output_dir` (two coarse smoke runs overwrote job
11726819's official results that exact way on 2026-08-07), and it shares
`save_dir` at the real config so its 2 completed solves become cache hits
for the full re-run rather than throwaway work.

### Step 7 — metrics helper (A8) — ✅ math primitives landed 2026-08-05; the plot/report itself still pending Step 6

**Deliverable.** Imbalance, excess loss, CMRR_I/Q, and phase-error
functions per A8's equations, plus the phase-consistency check A8 flags as
a new risk class (verify the hybrid script's cross-excitation `Δφ_j`
against the full-device FDTD reference from Step 6b, not just against
itself). Location (script-local vs. `SMatrix` method) per A8's open
decision.

**✅ Landed:** `SMatrix.imbalance_db`/`excess_loss_db`/`cmrr_db`/
`cross_excitation_phase_diff_deg` — see A8's own boxed update for the full
account. **Still pending, blocked on Step 6:** there is no real device yet
to plot, and the phase-consistency risk these primitives cannot answer on
their own.

**Evidence.** A CMRR-vs-λ / phase-error-vs-λ / excess-loss-vs-λ plot in the
shape of the paper's own Fig. 3, annotated against the **cited system
requirements** (CMRR ≥ 20 dBe, |Δφ_j| ≤ 5°) — not against the paper's own
silica numbers (A9). Not yet produced — needs Step 6's real S-matrices.

### Plotting/reporting improvements — added 2026-08-07, at the user's request

Every one of the three scripts' results plots (`hybrid_results.png`,
`eme_results.png`, `fdtd_results.png`) was a 1×2 (CMRR, phase error) —
extended to a 2×2, adding:

- **Excess loss vs. wavelength** (Eq. 11), one curve per driven input port.
- **Quadrature phase vs. wavelength** — the *actual* measured
  `φ_j = ∠S_j3 − ∠S_j1` (A11's sign-corrected direction) for all four
  output ports, each against its own paper-cited target angle (dotted
  reference line, wrapped to the same `(-180, 180]` range
  `cross_excitation_phase_diff_deg()` returns) — a more direct view of the
  I/Q structure than the existing Δφ (residual-only) panel: shows the
  actual ~90°-spaced port phases, not just their error from target.

Regenerated for the two scripts that already have real data
(`mmi_4x4_90hybrid_hybrid.py`/`_eme.py`, jobs 11723505/11726819) **from
their already-saved `SMatrix` `.pkl` files — zero new cluster time**, via
a one-off loader script; the printed metrics matched the previously-
recorded values exactly, confirming no drift. Visually, both plots show
the four ports sitting close to their expected quadrature positions
(hybrid: φ_5≈−45°, φ_6≈−130°, φ_7≈+46°, φ_8≈+134°, each within a few
degrees of its target dotted line). `mmi_4x4_90hybrid_fdtd.py` got the
same plotting code, ready for when job 11733865 lands — nothing to
regenerate there yet.

**Notebooks added for all three scripts**, `examples/mmi_4x4_90hybrid/
notebooks/*.ipynb`, matching this repo's own established per-device
convention (`notes.md`'s own item 6 flags this hand-copy pattern as a real
source of past drift across 12 other examples — mitigated here, not
solved repo-wide, by generating all three mechanically from a small
converter script rather than hand-copying: module docstring → markdown
title cell, code split at the scripts' own `# ---...---` section
delimiters, `matplotlib.use("Agg")` dropped. Verified byte-identical
code content against each script directly, not just visually).

### Excess-loss investigation — started 2026-08-07, at the user's request ("is there any way to reduce the loss?")

**Hypothesis going in:** `taper_width=0.8`/`taper_length=1.0` are unmeasured
placeholders (Step 2) — over just 1.0 µm the access waveguide widens
0.35→0.8 µm, a ~12° per-side taper half-angle, aggressive for an adiabatic
transition in strongly-guided 220 nm SOI. Proposed a taper_length scoping
sweep, EME-only, to check cheaply before committing FDTD time. User: "sure."

**Built** `examples/verification/python/mmi_4x4_90hybrid_taper_sweep.py` —
same method/geometry as `mmi_4x4_90hybrid_eme.py`, sweeping `taper_length`
only (taper_width fixed), ONE driven input (port 1) to keep it cheap (taper
loss is local physics near the access waveguide, not expected to depend on
edge- vs. interior-port position the way phase/CMRR do — flagged, not
independently checked for this side investigation).

**Smoke test (login node, RES=15, N_MODES_TAPER=2, N_MODES_BODY=4, N_WL=1,
taper_length∈{1.0,3.0}): cheap (29-39s per 4-call point) but the excess
loss AND imbalance came back IDENTICAL to 3 decimals between the two
lengths** — not just similar, exactly equal. Suspicious rather than
reassuring; investigated with zero further cluster cost (isolated
single-port-pair `solve()` calls on the login node) before trusting or
discarding the sweep:

- **At fixed `num_cells=15` (this device's taper cell count, inherited
  unexamined from `mmi_2x2_eme.py`'s own convention), `|S|²` for one
  representative port pair is invariant to `taper_length` to 3-4 decimal
  places** — confirmed at both `N_MODES_TAPER=2` (0.010538 at tl=1, 3, 10 µm)
  and `N_MODES_TAPER=4` (0.063994-0.063997 at tl=1, 3, 10 µm). Only the
  *phase* moves substantially with length (as expected — propagation phase
  accumulates over distance); the *transmitted power* does not.
- **Root cause, confirmed directly: it's a discretization artifact, not
  real taper-length physics.** At `N_MODES_TAPER=4`, fixing `taper_length`
  and sweeping `num_cells` instead (5→15→30→60→100) gave |S|² = 0.05910 →
  0.06400 → 0.06559 → 0.06645 → 0.06679 — a clear, still-climbing trend
  with diminishing increments (+0.0049, +0.0016, +0.0009, +0.0003), and
  **identical values (to 4 decimals) between taper_length=1.0 and 3.0 at
  every one of those `num_cells` settings.** The taper's linear width
  profile is discretized into `num_cells` uniform steps of FIXED width-jump
  size `(taper_width−waveguide_width)/num_cells`, independent of physical
  length — so with `num_cells` fixed, lengthening the taper changes nothing
  about the per-step mismatch loss this model computes; only a finer
  staircase (more `num_cells`) does.
- **`num_cells=15` is measurably under-converged.** The 15→100 climb
  (+0.00279, diminishing but not flat) is roughly a 0.18 dB improvement in
  transmitted power at just one taper interface, not yet plateaued at 100
  (extrapolating the shrinking increments, perhaps another ~0.1 dB left
  before a true asymptote).
- **Correction, checked directly before acting on this (2026-08-07):
  this artifact does NOT touch Step 6c's reported 2.1-4.3 dB excess loss.**
  Grepped `num_cells` across all three official scripts + the body sweep:
  only `mmi_4x4_90hybrid_eme.py` (Step 6a, `num_cells=[1,15,1,15,1]`)
  subdivides its tapers via this EME staircase at all. The hybrid script
  (`mmi_4x4_90hybrid_hybrid.py`, Step 6c — the more physically complete,
  already-reported result) and `mmi_4x4_90hybrid_body_sweep.py` both use
  `INPUT_METHOD="fdtd"` for their input/taper section (`num_cells=[1]` only
  ever appears for their own uniform, taper-free EME *body* sub-solver) —
  real time-stepped FDTD over the continuous taper geometry, resolved by
  grid `RES` the normal way, no discretization staircase of this kind.
  **This finding is scoped to Step 6a's own EME-only numbers only** — an
  initial framing of this as touching the settled 2.1-4.3 dB hybrid figure
  was wrong and has been corrected here before it propagated further.
  Each signal path in the EME-only script crosses two tapers (in + out),
  so **a rough 0.3-0.6 dB slice of Step 6a's own reported numbers**
  (already flagged as the less-trusted of the two methods, per the
  EME-vs-hybrid divergence finding above) **is plausibly this numerical
  artifact, not physical radiation loss** — not measured to a firm number,
  extrapolated from a clear,
  monotonic, diminishing-returns trend on a coarse (RES=15, N_MODES_BODY=4)
  proxy, not the production settings.
- Cost escalates fast with `num_cells`: 100 cells took 759s (12.7 min)
  serial on the login node for ONE port pair at these coarse settings; 200
  was killed after >28 min still running (not needed — the diminishing
  increments already say enough without chasing the last ~0.05 dB).

**Net finding: taper_length is not (yet) a demonstrated lever in this
model — num_cells is, and only in `mmi_4x4_90hybrid_eme.py` (Step 6a).**
The originally-planned taper_length sweep can't show a real physical
trend until `num_cells` is bumped to something more converged first (the
effect was being masked entirely in this EME-only script).

**Actioned 2026-08-07, user's direction ("err on the lower side as 60
cells is a lot"): bumped `mmi_4x4_90hybrid_eme.py`'s taper `num_cells`
from 15 to 30** — a 2× step, not the 60-100 the convergence check used,
trading some residual under-convergence for cost (the diminishing-returns
data above says 30 already captures roughly half of the measured 15→100
gain). Only this one script changed — the hybrid (6c) and body-sweep
(Step 3) scripts don't use this parameter for their tapers at all (they
use FDTD there), so nothing about their own already-reported numbers is
affected or needs re-running.

**Cost check before committing production time — isolated single calls
were misleading, measured directly rather than trusted.** An isolated
`solve()` call at `num_cells=30` (coarse settings) took ~240s vs.
`num_cells=15`'s ~29s — an alarming 8x that would extrapolate to hours for
the full 8-call production run. Before accepting that estimate, ran the
ACTUAL script (not isolated calls) at the same coarse settings for both
`num_cells` values: the full 8-call batch took 74s at 15 and 133s at 30 —
only ~1.8x, not 8x. The isolated-call number was dominated by a one-time
JAX JIT-compilation cost for that specific array-shape configuration, paid
once and reused across the other 7 same-shape calls within one script run
— invisible in an isolated single call, which pays that cost with nothing
to amortize it against. Applying the real ~1.8x to the original
production run's own measured time (job 11726819, 40:46 for all 8 calls)
gives an estimated **~73 min** for the `num_cells=30` re-run — comfortably
inside the existing `-t 360` (6h) SLURM budget, left unchanged.

**Real near-miss, caught and fixed, not shipped:** the coarse-setting
smoke invocations above wrote to the SAME `example_output/.../
mmi_4x4_90hybrid_eme/` directory the real script always uses (env-var
overrides change solver settings, not the output path) — so two smoke
runs briefly overwrote the official job 11726819 results (`eme_results.
npz`/`.png`, `smatrix_eme_port{1,3}.pkl`) with coarse garbage. Caught
immediately (the file mtimes/values didn't match what should be there),
and recovered with zero data loss: re-ran the script with
`MMI4X4_EME_NCELLS_TAPER=15` and every other setting at its production
default, which hit the on-disk solve cache instantly (11s for all 8 calls)
and reproduced job 11726819's own recorded numbers exactly (imbalance,
CMRR, Δφ_6/7/8 all matched to the displayed precision) — confirming the
underlying cache was never touched, only the output directory's own
files were transiently clobbered. Lesson for next time: point smoke-test
env-var experiments at a scratch output path, not the script's own
official `output_dir`, even when only solver *parameters* are overridden.

**✅ Job 11736983 COMPLETED 2026-08-07 (1:13:27 — matched the ~73 min
estimate closely).** Real, measured before/after comparison (num_cells=15,
job 11726819, vs. num_cells=30, job 11736983), all 3 wavelengths:

| | 1.26 µm | 1.31 µm | 1.36 µm |
|---|---|---|---|
| excess loss, port 1 (15→30) | 3.570→3.475 dB | 2.325→2.237 dB | 2.440→2.359 dB |
| excess loss, port 3 (15→30) | 5.299→5.200 dB | 2.582→2.496 dB | 3.099→3.018 dB |

**A real, consistent, but modest improvement: ~0.08-0.10 dB per port at
every wavelength** — smaller than the ~0.3-0.6 dB extrapolated earlier
from the coarse single-port proxy (that estimate was from a very different,
unconverged settings regime — RES=15/N_MODES_BODY=4 — not a substitute
for measuring it at the real production settings, which is what this is).
Imbalance/CMRR/phase error all shifted by only hundredths of a dB/degree
— consistent with excess loss being the metric this specific numerical
artifact should affect (it changes total transmitted power), while
imbalance/CMRR/phase are ratios between output ports that a
roughly-uniform-across-ports discretization effect should mostly cancel
out of, and did.

**Net: num_cells was a real lever, but a small one at the 15→30 step
taken.** `num_cells=30` is itself still short of the ~15→100 convergence
climb measured earlier (60-100 kept climbing on the coarse proxy) — more
of this specific correction is available at higher `num_cells`, at a
steeper, non-trivially-scaling cost (see the cost-check above). Not
pursued further without the user's direction, given the returns here were
modest relative to the still-unaddressed EME-vs-hybrid divergence and the
still-unresolved taper_length question (masked by `num_cells=15` before,
now only partially unmasked at `num_cells=30` — the original
`mmi_4x4_90hybrid_taper_sweep.py` scoping sweep has not been re-run at the
new `num_cells=30`).

### Taper-length re-check at num_cells=30 — 2026-08-07

Re-ran `mmi_4x4_90hybrid_taper_sweep.py` (updated to default `num_cells=30`,
matching `mmi_4x4_90hybrid_eme.py`'s now-settled value) at a coarse smoke
setting (`N_MODES_TAPER=2`, `N_MODES_BODY=4`, `RES=15`, taper_length ∈
{1.0, 3.0}). **Still exactly flat** (excess_loss=10.802 dB, imbalance=7.41
dB, identical at both lengths) — bumping `num_cells` alone was NOT enough
to unmask a length trend at this mode count. Cross-referencing the earlier
diagnostic data: `N_MODES_TAPER=4` (this script's own real-sweep default)
also showed flatness to displayed precision at every `num_cells` tested
(5/15/30/60/100, all identical between taper_length=1 and 3). The ONLY
mode count that showed ANY resolvable taper_length dependence at all was
`N_MODES_TAPER=8`, and even there the effect was tiny (~1.3% relative
`|S|²` change from tl=1→10 µm) and pointed the "wrong" way (loss
*increasing* with length) — more consistent with a residual numerical
artifact than a real adiabatic-taper benefit, not independently confirmed
at a higher `num_cells` either.

**Not yet run: the full production-parameter sweep** (6 taper lengths ×
4 outputs, `N_MODES_TAPER=4`/`N_MODES_BODY=8`/`RES=40`) — estimated at
roughly an hour of real cluster time by analogy to the EME production
run's own per-call cost, likely to reconfirm "no resolvable trend," given
every mode count tested below production settings already shows one.
Paused here to check with the user before spending that time chasing a
signal that looks null at every setting actually usable in the settled
scripts.

**Decision 2026-08-08 (user's direct call): stop here.** Every setting
tried (num_cells 5-100 x N_MODES_TAPER 2/4, and now num_cells=30 x
N_MODES_TAPER=2) shows a flat or negligible taper_length effect — not
pursuing the full production-parameter sweep to chase a signal already
this consistently null. `taper_length=1.0` (MMI's own class default,
Step 2's original placeholder) stands unchanged; the only real,
production-validated loss-reduction lever from this whole investigation
is the `num_cells` 15→30 bump already landed in `mmi_4x4_90hybrid_eme.py`
(job 11736983, ~0.08-0.10 dB real reduction).

### Body/taper mode-count convergence investigation — 2026-08-09

User request: field plots (examples/verification/python/
plot_mmi4x4_90hybrid_eme_fields.py) showed taper section under-resolved.
Asked N_MODES_TAPER>10.

**Taper side.** N_MODES_TAPER 4->12 in mmi_4x4_90hybrid_eme.py. Field
plot dramatically smoother at 12 (fans into clean self-imaging pattern
vs staying confined/speckly at 4) -- BUT production S-matrix metrics
barely moved (job 11746259 vs original: agreement to 3-4 decimals).
Convergence sweep (job 11747369, single port pair, num_cells=30): flat
4-16, anomaly at 20. Root cause found: taper's wide end (0.8um) has
**6 real guided modes** (checked directly, meep backend, after catching
a strip-tuple-order bug in first attempt) -- N=4 missed 2 of 6, N=12
covers all + margin. Flatness explained: fundamental-mode launch barely
couples to the missing higher-order guided modes for this geometry --
real physics, not luck. N=20 anomaly: likely instability among
non-guided continuum pseudo-modes (modes 7+), not a guided threshold
(guided count is 6, nowhere near 20). Net: N_MODES_TAPER=12 kept,
well-justified now, not just "field plot looked better."

**Body side -- the big one.** S-parameter convergence sweep (job
11752350, tidy3d backend, single port pair, N=[4,8,12,16,20,24]):
**N_MODES_BODY=8 (used in EVERY script for this device) under-converged
by 50-64% relative to true value**, which converges cleanly by N=16-20
(agreement <1% through 24). This is not the taper's small effect --
N=8 gives |S|^2=[0.139,0.150,0.147], true value ~[0.085,0.149,0.223].

Body sub-solver deliberately locked to "meep" backend (not tidy3d) for
FDTD/EME mode-basis consistency across hybrid cascade (body_sweep.py's
own comment). Bare N=16/20 guided-ceiling probes via meep at the WIDE
diagnostic domain (oxide_width=MM_WIDTH+14.0, body_modes.py's own
convention) TIMED OUT at both 3h (N=16, job 11753042) and 4h (N=20, job
11749479) -- ZERO output both times, same signature as the original
Step 3a 16-mode failure. User direction: skip ceiling-finding, go
straight to S-param convergence (done above); then investigate
meep/tidy3d interchangeability before accepting the wall.

**meep vs tidy3d comparison** (N=8, both handle this fine): n_eff
agrees to ~1% (consistent small offset, all 8 modes). Overlap matrix
between the two mode sets: perfectly diagonal (zero cross-terms,
same ordering/content) but diagonal magnitude only ~0.05 even with
power-normalized inner product (each self-normalizes to 1.0
individually). Dug deeper: sign-change count along raw field profile
--meep 2/7/11 for modes 1/4/8 (smooth), tidy3d 101/182/111 (noisy).
Envelope/n_eff/width/position all agree closely; tidy3d's raw field
carries real numerical noise that would corrupt a strict interface
overlap. **Conclusion: NOT clearly safe to swap body backend to
tidy3d.**

Bisected meep's own ceiling instead: N=12 timed out at 2h (job
11764161, zero output), N=10 probe in flight when superseded (see
below). **User asked: have we checked the SIMULATION SPAN (domain
size) for the mode calc, not just mode count?** Real catch -- checked
oxide_width/oxide_height actually used: body_modes.py/cross_section.py
(and every prior "meep wall" test) used oxide_width=MM_WIDTH+14.0=20um;
body_sweep.py/hybrid.py's ACTUAL production domain is oxide_width=
MM_WIDTH+3.0=9um. Never the same span.

Domain convergence screen (tidy3d, cheap/fast, N=8): n_eff BIT-FOR-BIT
IDENTICAL across oxide_width margin 3.0 through 20.0, and oxide_height
2.0 through 6.0. Domain was never too small for ACCURACY.

Re-tested meep at the REAL production width (9um): **N=16 completed in
77 min (job 11767748), N=8 in 23 min (job 11767747)** -- both faster
than anything seen at the wide domain, no wall at all. Root cause of
the "meep wall": the WIDE diagnostic-only domain gave the eigensolver
far more numerical room for near-degenerate/spurious modes at higher
orders (MPB log: mode 16 needed 160 iterations at the wide domain vs
few iterations for low-order modes) -- a property of the wrong test
config, not of meep or N=16 itself. Cancelled the now-moot wide-domain
N=10 bisection job (11767698).

**Full redo, as originally approved, now unblocked.** N_MODES_BODY
8->16 in body_sweep.py (production oxide_width unchanged, already
correct). Step 3 ΔL re-sweep: job 11768847, COMPLETED 8:58h (budget
15h). **Result: optimal ΔL shifted -1.50 -> -2.75 um** (mm_length
78.9396 vs old 80.1896) -- best-point metrics WORSE than the N=8
result: max|Δφ|=12.320° (was 8.28°), max imbalance=14.6155dB (was
4.90dB). The old "close to spec" result was itself an artifact of the
under-converged basis, not real.

hybrid.py updated to N_MODES=16, DELTA_L=-2.75 -- NOT a cache hit
(this script's own boundaries derive from mmi.mm_length, which shifted;
same mechanism as the original mm_length near-miss). Job 11778128,
COMPLETED 10:10h (budget 15h), all 4 FDTD input sections + body solve
genuinely fresh.

**Final, converged Step 6c numbers (N=16, ΔL=-2.75), 1.26/1.31/1.36um:**

| | port 1 | port 3 |
|---|---|---|
| imbalance (dB) | 4.267 / 1.928 / 5.871 | 14.615 / 1.744 / 6.495 |
| excess loss (dB) | 6.460 / 1.776 / 2.121 | 8.366 / 1.452 / 3.959 |
| CMRR_I (dBe) | 6.836 / 13.214 / 10.480 | 10.415 / 24.389 / 5.578 |
| CMRR_Q (dBe) | 10.417 / 24.385 / 5.578 | 1.531 / 14.241 / 7.784 |
| Δφ_6/7/8 (deg) | 2.291 / 10.449 / 9.971 | -4.324 / 7.514 / 12.320 | (third col 9.973/9.347 at 1.36um, see report)

**Net: this device's real, converged performance is meaningfully worse
against the paper's cited requirements (CMRR>=20dBe, |Δφ|<=5°) than
the earlier under-converged (N=8) numbers suggested** -- CMRR fails at
most wavelengths now (was: comfortable at centre), phase error breaches
5° at multiple ports/wavelengths now (was: only Δφ_6). Excess loss also
higher (6.5-8.4dB peak vs old 2.1-4.3dB). This is the correct, trustworthy
number -- not a regression to fix, a measurement error to stop trusting.
Whether the geometry itself (mm_width, taper params, access waveguide
width) needs redesigning to actually meet spec is now an open question
this convergence work exposed, not one it answers.

### Top-bottom mirror-symmetry solve reuse — checked 2026-08-13, zero new cluster time

Question: `mmi_2x2_hybrid.py`'s own `SYMMETRIC_MODE_REUSE`/`MODE_SIGNS`
mechanism (top vs. bottom arm, a ±1 sign vector over body modes, empirically
measured once) already drops that device's FDTD count from 2 to 1, on top
of the left-right (input↔output) mirror reuse every MMI hybrid example
here already uses (A5). This device has the same y=0 mirror symmetry
(`_port_ys` places all 4 ports symmetrically, `_default_ports` gives every
port identical `waveguide_width`/`taper_length`/`taper_width` —
`beef/geometry/mmi.py:58-77`), across TWO mirror pairs (`input_idx` 0↔3 at
y=∓2.25, 1↔2 at y=∓0.75) — so the same trick would in principle drop the 4
input-section solves to 2, the natural floor for a 4-port device (2 distinct
`|y|` magnitudes).

**Checked directly against job 11778128's already-cached results
(`solve_section(..., cache_only=True)`, all 4 `input_idx` cache HIT, 0.17 s,
no new FDTD run)**, not just architecturally argued:

- `M_in[3]/M_in[0]` and `M_in[2]/M_in[1]` (both mirror pairs, all 3
  wavelengths, N_MODES_BODY=16): **all 16 of 16 body modes give a clean,
  real, wavelength-stable ±1 sign** — checked with proper circular
  statistics (`R = |mean(unit ratio)|`), not a naive arithmetic mean/std of
  the angle (a first pass using the latter wrongly flagged 2 modes as
  "unstable" — a branch-cut artifact: angles landing at alternating
  +179.9°/−180.0° across wavelengths are the *same* physical sign,
  literally 180° apart from 0 either way, but a plain `.mean()`/`.std()` on
  wrapped angle values doesn't know that and reports a spuriously huge
  swing. Caught before it propagated into anything real by re-checking
  with `R` — `R=1.0000` (min 0.99999995) for every one of the 16 modes,
  both pairs, meaning perfect wavelength-to-wavelength consistency, no
  exceptions). Angular deviation from an exact 0°/180° sign: ≤0.04° (pair
  0↔3), ≤0.15° (pair 1↔2) — every mode. Magnitude residual `||ratio|-1|`:
  ≤0.17% (pair 0↔3), ≤0.29% (pair 1↔2) — every mode. **The same alternating
  sign pattern `[+1,-1,+1,-1,...]` holds for BOTH pairs**, i.e. one shared,
  single sign vector (not two independently-fit ones) explains both 0↔3
  and 1↔2 — real evidence the sign is an intrinsic y-parity property of
  each body eigenmode, not a per-pair coincidence.

**Verdict: the reduction is real, clean on all 16 modes (no flagged
exceptions), and at a precision comparable to the 2×2's own `MODE_SIGNS`
validation (~1e-4-level angle agreement, sub-0.3% magnitude residual) — not
just a screening approximation.** Still worth one caveat before treating it
as equal-status to the 2×2's own precedent: this is one reference geometry
(N_MODES_BODY=16, `mm_length=78.9396`, the current design point) checked
once, not re-derived across a length/wavelength range the way the 2×2's
`MODE_SIGNS` comment describes being re-validated after design changes —
worth a quick re-check (same cache-only technique, ~free) if `mm_length`/
`ΔL` or `N_MODES_BODY` ever change again.

### `waveguide_offset` body-cache independence — checked 2026-08-13, zero new cluster time

Before spending real solves on an offset sweep: does `EMESolver`'s cached
body-mode result (job 11768847's 81-point `ΔL` grid, N_MODES_BODY=16) stay
valid if `waveguide_offset` changes? `waveguide_offset` only moves the four
ports' y-position, entirely inside the *input* section's own x-range
(`[boundaries[0], boundaries[1]]`); the body section's cache fingerprint
(`beef/cache/config.py`'s `build_run_config()`) is built from
`geometry.section_gds(x0, x1, ...)` — the polygons **sliced to
`[boundaries[1], boundaries[2]]` only** — so an offset change cannot appear
in it. **Confirmed directly, not just argued from the fingerprint code**:
called `solve_section_length_sweep(..., cache_only=True)` with
`waveguide_offset=1.3` (a candidate sweep point, never solved before)
against the same body boundaries/lengths grid — cache HIT, identical
`(81, 3, 16, 16)` shape, no new run. **The body solve never needs
re-running across an offset sweep** — only the 2 (via the mirror-reuse
trick above) or 4 input-section solves per offset point do.

### Step 8 — the four checks

Once Step 6's three scripts exist and run for the first time, this doc
converts to a normal gap list against
[README.md's four checks](README.md#the-four-checks):

- **Check 1 (engine agreement).** Uniform-body-stretch phase fit
  (`verify_item_6_engine_agreement.py`-style), on this device's own
  cross-section — cannot borrow any existing device's ratio (every prior
  measurement is at a different (platform, band) pair).
- **Check 2 (mode identity).** Elevated risk here per A7 — four
  closely-spaced ports per end is closer to the modal multiplexer's
  16-mode regime (where check 2 first failed, on a real bug) than to either
  existing MMI's two-or-one-port ends.
- **Check 3 (resolution convergence, per engine).** Sweep independently at
  the new (SOI_220, O-band) cross-section — the 2×2's res 100 (C-band) does
  not transfer, same logic as every prior "new pair" in this repo.
- **Check 4 (full-device FDTD reference).** Per A8, this device's check 4
  carries a strictly stronger requirement than any prior device's: it must
  reproduce `Δφ_j`, not just total transmitted power.

### Step 9 — golden fixture

Add to `tests/test_golden_examples.py`'s `TARGETS` **after** Step 8's
resolution question settles — every other device's doc repeats this
ordering for the same reason (recording a golden against a configuration
still expected to move means regenerating it twice).

---

## Open decisions carried into Part B, not made unilaterally here

1. **A3** — self-imaging + EME sweep vs. porting a real optimizer, if the
   sweep can't hit the cited requirements.
2. **A8** — script-local metric helpers vs. new `SMatrix` methods.
3. **Exact taper length/width, `waveguide_offset`, and whether any per-port
   `Δy_i` fine-tuning (paper Eq. 8's own knob) turns out to be needed** —
   all explicitly TBD in Step 2/3, not guessed here.
4. **Directory/file naming** — used `mmi_4x4_90hybrid` above for
   consistency with `mmi_2x2`/`mmi_1x4`'s own `mmi_<NxM>` convention plus a
   disambiguating suffix for the terminology collision flagged at the top;
5. ~~Coarse `waveguide_offset` sweep~~ ✅ **DONE 2026-08-13, job 11817429 —
   then extended into a full joint `(offset, ΔL)` grid, zero new cost,
   same day — see "Joint (offset, ΔL) grid" section below. Answered
   firmly: offset alone is not the fix, and re-optimizing ΔL per offset
   does not rescue it either.** Every tested offset (1.2, 1.35, 1.65, 1.8)
   was dramatically worse on phase error at the shared ΔL (51°–163° vs.
   baseline 1.5's 12.3°); the sign-vector derivation was independently
   confirmed correct at a new offset (≤0.84° angle / ≤0.2% magnitude vs. a
   real 4-solve check at offset=1.35). Extending to the FULL 81-point ΔL
   grid per offset (already-cached data, `mmi_4x4_90hybrid_offset_deltaL_
   grid.py`, 21 s on the login node) shows even each offset's OWN best
   achievable ΔL stays far worse than baseline (best max\|Δφ\| = 125.6° /
   64.8° / 47.6° / 75.1° at offset 1.2/1.35/1.65/1.8, vs. baseline's
   12.3°) — the paper's own Eq. 8 port-spacing formula (`offset=W/4`) is
   not merely "what ΔL happened to be tuned for," it is genuinely close to
   optimal across the explored range, independent of ΔL. `waveguide_offset`
   is now a closed, ruled-out lever, not an open one.
6. ~~Top-bottom mirror-symmetry reuse for production reporting, not just
   screening~~ ✅ **No longer blocked** — the original "2 flagged
   near-degenerate modes" finding was itself an analysis artifact (a
   branch-cut bug in how the check averaged wrapped angles, not a real
   mode-basis instability); corrected check shows all 16 modes clean. The
   sign vector is trusted for this design point; re-validate (cheap,
   cache-only) if `mm_length`/`N_MODES_BODY` change again.
   open to a different name before Step 6 creates the directory.

---

## References

Primary reference for this device's topology, equations, and terminology:

> T. Yasui, J. Sugisaka, and K. Hirayama, "Structural Optimization of an
> Optical 90 Degree Hybrid Based on a Weakly-guided 4×4 Multimode
> Interference Coupler Using a Parallelized Real-coded Micro-genetic
> Algorithm," *ACES Journal*, vol. 36, no. 5, pp. 526–532, May 2021.
> https://doi.org/10.47037/2020.ACES.J.360506

Two papers it cites that matter directly to this plan's decisions (A3, A9):

> R. Halir, G. Roelkens, A. Ortega-Moñux, and J. G. Wangüemert-Pérez,
> "High-performance 90º hybrid based on a silicon-on-insulator multimode
> interference coupler," *Opt. Lett.*, vol. 36, no. 2, pp. 178–180, Jan.
> 2011. — the SOI-based precedent motivating A3's self-imaging-without-GA
> approach, and one of the two sources (with [5] below) for the "typical
> system requirements" (CMRR ≥ 20 dBe, |Δφ| ≤ 5°) A9 checks against instead
> of the paper's own silica numbers.

> L. B. Soldano and E. C. M. Pennings, "Optical multi-mode interference
> devices based on self-imaging: Principles and applications," *J. Lightw.
> Technol.*, vol. 13, no. 4, pp. 615–627, Apr. 1995. — the general
> self-imaging theory Eqs. 7–9 are drawn from, referenced in A3 for why it
> is expected to hold well on a strongly-guided platform even though it
> fails on the paper's own weakly-guided one.

> J. S. Fandino and P. Munoz, "Manufacturing tolerance analysis of an
> MMI-based 90º optical hybrid for InP integrated coherent receivers,"
> *IEEE Photon. J.*, vol. 5, no. 2, pp. 7900512, Apr. 2013. — the second
> source (with Halir above) for the CMRR/phase-error system requirements.


---

## Toward a (mm_width, y_gap) Bayesian optimization — the bias study first (2026-08-25)

User's plan: treat the device as unoptimized in length, multimode width and
taper geometry, and run a **2D Bayesian optimization over (mm_width, y_gap)**
in EME with tidy3d modes, using `solve_length_sweep()` to find the best
`mm_length` at each candidate. Objective: a weighted combination of phase
error and CMRR, at 1.31 ± 0.02 µm.

### The parameterization

Port centres are at `waveguide_offset * (i - 1.5)`, so adjacent port spacing IS
`waveguide_offset` and the tip-to-tip gap between neighbouring tapers is
`y_gap = waveguide_offset - taper_width`. Holding the paper's Eq. 8 relation:

    waveguide_offset = mm_width / 4
    taper_width      = mm_width / 4 - y_gap

Two free parameters fix the whole geometry. **Sanity check: the current design
(mm_width=6.0, taper_width=0.8) is exactly y_gap=0.70**, matching the 0.7 µm
port clearance `mmi_4x4_90hybrid_eme.py`'s own docstring already quotes.
Validity: `0 < y_gap < mm_width/4 - 0.35`.

### Why a bias study has to come first

`EMESolver.solve()` draws ONE input and ONE output branch per call and treats
the other six ports as physically absent — no crosstalk modelled at all. That
restriction is deliberate, not an oversight: `mmi_1x2_eme.py`'s own comment
explains it gives "the isolated waveguide mode rather than a supermode of both
outputs". Removing it means drawing all ports, raising end-cell `num_modes`
from 1 to ≥4, and projecting supermodes back onto per-port modes — squarely in
this repo's most fragile area (`MODE_SIGNS`, `check_mode_correspondence`, the
modal-mux basis gap).

Normally a systematic approximation is tolerable in a BO, because a constant
bias cancels in a *relative* ranking. **That argument fails here**: `y_gap` IS
the tip-to-tip separation that sets the neglected crosstalk, and coupling is
roughly exponential in gap. A bias that varies with the optimized parameter
does not cancel — it steers the optimizer. Check 4 already measured this
approximation costing up to 10° of phase and 6.6 dB of imbalance at y_gap=0.70.

### The study — 7 jobs, submitted 2026-08-25

EME (`mmi_4x4_90hybrid_ygap_eme.py`, jobs 12352460/61/62) and full-device FDTD
(`mmi_4x4_90hybrid_ygap_fdtd.py`, jobs 12352475-78) at **y_gap = 0.40 / 0.70 /
1.00**, `mm_width=6.0`, `mm_length` FIXED at 78.9396 so the two methods are
compared on identical geometry. y_gap=0.70's FDTD is REUSED from jobs
12272571/12272573 rather than re-run, saving two ~14 h jobs. Wavelengths
1.26-1.36 to match that reference — measuring the bias on the WIDER band is
conservative, since band-edge error is the worst case.

FDTD settings come from this session's own measurements: `decay_tol=1e-4`
(measured sufficient at 0.015 dB), `decay_method="energy"`, 96 ranks. At the
inherited `1e-11` these four runs would have been impossible.

**Verdict rule, fixed in advance:**
- error roughly FLAT in y_gap → largely cancels → the 2D BO on plain EME is sound
- error VARIES strongly → the BO would be biased → fix the solver, narrow the
  y_gap range, or optimize mm_width alone

The EME script IS the BO's evaluation function, so this validates exactly what
the optimizer will call, objective included.
---

## The (mm_width, y_gap) optimization campaign — 2026-08-26 → 2026-09-01

**Documentation gap closed 2026-09-01.** Everything in this section happened
after this file's last edit (2026-08-26) and was undocumented here: the bias
study's close-out, the hybrid-path evaluator, Phase A bodies, the 106-point
grid's replacement by a Bayesian optimizer, and 24 real BO evaluations. The
crash-recovery snapshot `context/context_2026-08-28_0913_4x4_hybrid.md` covers
Aug 26-28 only, and calls itself non-durable; this is the durable record.

### ✅ Pure-EME BO DISQUALIFIED — bias study CLOSED, all three points (2026-09-01)

The verdict rule fixed in advance above ("error VARIES strongly → the BO would
be biased") fired. Full table, `max|Δφ|` over the 1.26-1.36 µm band,
`mm_width=6.0`, `mm_length=78.9396` on both methods:

| y_gap | EME max\|Δφ\| | FDTD max\|Δφ\| | EME error | FDTD min CMRR |
|---|---|---|---|---|
| 0.40 | 66.293° | 22.436° | **+43.86°** | 9.746 dBe |
| 0.70 | 82.478° | 20.46° | **+62.02°** | 9.19 dBe |
| 1.00 | 89.203° | 53.707° | **+35.50°** | 1.264 dBe |

- The y_gap=1.00 row was the one still missing on 2026-08-28 (job 12352478,
  port 3, was still running at 22h12m). It landed; **computed 2026-09-01 from
  the npz already on disk, zero new cluster time**, by rebuilding `SMatrix`
  objects from `ygap_study/fdtd_ygap1.00_{bot,idx2}.npz`'s `c_out_*` and
  running `mmi_4x4_90hybrid_check4_metrics.py`'s own `metrics()` on them.
  Pipeline verified first by reproducing the already-published y_gap=0.40
  FDTD number (22.436° vs. 22.44° documented) before trusting the new row.
  The y_gap=0.70 FDTD entry is Check 4's cached full-device run, not an
  `ygap_study` npz — that point was deliberately reused, not re-run.
- **The bias is not a constant offset and not even monotonic**: +43.9°,
  +62.0°, +35.5° as y_gap goes 0.40 → 0.70 → 1.00, a 26° swing against a true
  (FDTD) signal that itself only moves 33° over the same range. Nothing
  cancels in a relative ranking, and no fitted offset repairs it.
- **The two methods rank the candidates differently**: FDTD orders
  0.70 < 0.40 << 1.00 (best to worst); EME orders 0.40 < 0.70 < 1.00. EME gets
  the worst point right and flips the two that matter.
- Cause, as predicted: `EMESolver.solve()` draws one input and one output
  branch per call and models no port crosstalk, and y_gap *is* the tip-to-tip
  separation that sets crosstalk.

### The evaluator that replaced it — `mmi_4x4_90hybrid_wg_eval.py`

Hybrid path, not pure EME: 4 real FDTD input sections with **all four ports
present simultaneously** (so crosstalk is real), cascaded with an EME body at
`mode_solver="meep"` (same-engine bases, per `body_sweep.py`'s convention).
Per point it also sweeps 81 lengths (`DELTA_L_GRID = linspace(-10, 10, 81)`,
centred on `0.75·L_π` scaled as `W²` from the `mm_width=6.0` measurement) and
returns the best — cheap, since tier-3 cache reconstructs `T` from cached
`n_eff` with no new mode solve per length. Scalar objective, lower-is-better,
0 = both specs met:

    objective = 0.5·(max|Δφ| / 5°) + 0.5·max(0, (20 − min CMRR) / 20)

**Measured per-point cost: 4:43:45 on 4 nodes** (job 12383245, the nominal
point with no pre-cached body) — roughly 1.0-2.25 h body + 3-3.7 h FDTD input
sections. That is 3-4× the earlier extrapolation from the offset sweep; budget
from this number, not that one.

### Why BO, not the 106-point grid

At the measured 3.5-4.7 h/point, 106 points at the batch script's own 8-way
concurrency cap is 45+ wall-clock hours *with nodes free on demand*, and this
campaign has seen single 4-node requests pend 3.5-22 h. The grid
(`mmi_4x4_90hybrid_ygap_grid.py`, `POINTS`) was never submitted as Phase B;
`mmi_4x4_90hybrid_bo.py` replaced it over the **same** per-point evaluator and
objective. Phase A's 5 body solves (job 12387804 + 12383245 backfilling
W=6.0) are why `mm_width` is **Categorical on {5.0, 5.5, 6.0, 6.5, 7.0}**: on
that lattice every candidate's body solve is a cache hit, and only the FDTD
input sections are paid per point. `y_gap` stays continuous (it only changes
which input sections are solved).

### The 24 evaluated points (2026-09-01, `ygap_campaign_summary.txt`)

| W | y_gap | taper | ΔL | mm_len | max\|Δφ\| | min CMRR | obj |
|---|---|---|---|---|---|---|---|
| 5.0 | 0.10 | 1.150 | −3.23 | 53.500 | 20.058 | 14.516 | 2.1429 |
| 5.0 | 0.43 | 0.820 | −3.58 | 53.150 | 19.933 | 5.600 | 2.3533 |
| 5.0 | 0.89 | 0.360 | −3.33 | 53.400 | 29.197 | 0.423 | 3.4091 |
| 5.5 | 0.34 | 1.035 | −4.14 | 64.500 | 15.977 | 9.710 | 1.8549 |
| 5.5 | 0.47 | 0.905 | −3.99 | 64.650 | 13.459 | 7.368 | 1.6617 |
| 6.0 | 0.10 | 1.400 | −5.24 | 76.450 | 10.278 | 8.990 | 1.3030 |
| 6.0 | 0.30 | 1.200 | −4.74 | 76.950 | 10.776 | 11.471 | 1.2908 |
| 6.0 | 0.40 | 1.100 | −4.44 | 77.250 | 9.771 | 11.576 | 1.1877 |
| 6.0 | 0.42 | 1.080 | −4.44 | 77.250 | 9.676 | 11.426 | 1.1819 |
| 6.0 | 0.44 | 1.060 | −4.39 | 77.300 | 9.239 | 11.186 | 1.1442 |
| 6.0 | 0.45 | 1.050 | −4.34 | 77.350 | 8.809 | 11.034 | 1.1050 |
| 6.0 | 0.50 | 1.000 | −4.24 | 77.450 | 7.666 | 10.477 | 1.0047 |
| 6.0 | 0.51 | 0.990 | −4.24 | 77.450 | 7.663 | 10.415 | 1.0060 |
| 6.0 | 0.65 | 0.850 | −4.29 | 77.400 | 7.473 | 9.606 | 1.0071 |
| 6.0 | 0.80 | 0.700 | −4.09 | 77.600 | 11.213 | 7.924 | 1.4232 |
| 6.5 | 0.10 | 1.525 | −5.57 | 90.300 | 6.709 | 7.641 | 0.9799 |
| 6.5 | 0.32 | 1.305 | −5.07 | 90.800 | 5.900 | 9.351 | 0.8562 |
| 6.5 | 0.51 | 1.115 | −4.87 | 91.000 | 4.788 | 12.081 | 0.6767 |
| **6.5** | **0.56** | **1.065** | **−4.87** | **91.000** | **4.764** | **12.109** | **0.6737** |
| 6.5 | 0.65 | 0.975 | −4.67 | 91.200 | 6.093 | 11.672 | 0.8175 |
| 6.5 | 0.85 | 0.775 | −4.12 | 91.750 | 10.824 | 7.939 | 1.3839 |
| 7.0 | 0.21 | 1.540 | −5.34 | 105.850 | 9.475 | 6.221 | 1.2920 |
| 7.0 | 1.11 | 0.640 | −3.44 | 107.750 | 17.026 | 3.419 | 2.1171 |
| 7.0 | 1.36 | 0.390 | −3.29 | 107.900 | 24.179 | 1.399 | 2.8829 |

**Best point: `mm_width=6.50`, `y_gap=0.56` (`taper_width=1.065`),
`mm_length=91.0000` (ΔL=−4.87) — max|Δφ| = 4.764°, min CMRR = 12.109 dBe.**

- **The phase spec is now MET** (4.764° ≤ 5°) — it was 12.32° at the settled
  geometry in Check 4, and 8.09° at the nominal point once length was
  optimized. This is the first geometry in this device's history to clear it.
- **The CMRR spec is not close, anywhere.** Best CMRR over all 24 points is
  14.516 dBe (W=5.0, y_gap=0.10) and that point's phase is 20.058° — 4× over
  spec. Along the good-phase ridge CMRR saturates at 11-12 dBe. Nothing in
  this 2D space approaches 20 dBe, and the two metrics trade off against each
  other: the objective's own CMRR term never reaches 0.
- Per-width picture: the objective improves monotonically 5.0 → 6.5 (2.14 →
  1.66 → 1.00 → 0.674), and each width has an interior optimum near
  `taper_width ≈ 0.85-1.12`, i.e. the good region tracks roughly constant
  *taper width*, not constant *y_gap*.

### ⚠ The GP has collapsed — MEASURED 2026-09-01, and the driver is now fixed

Round 6's `propose` submitted **one** job instead of four: three of four
candidates were `W=6.00` at `y_gap` 0.50/0.51/0.52, all inside an
already-measured neighbourhood (W=6.00 y=0.50 and 0.51 were in hand at
objective 1.0047/1.0060), and the driver's only guard was an exact-tuple
`(mm_width, y_gap)` match, so it silently skipped 3 and shrank the round.

Then measured directly: oversampling `opt.ask(n_points=16, strategy="cl_min")`
returns **16 of 16** candidates inside 0.03 µm of an evaluated point — 15 at
`W=6.00, y_gap=0.50-0.53` and one at `W=6.50, y_gap=0.55`. EI is effectively
zero everywhere else; the acquisition is fully exploitative and the campaign
cannot make progress by re-running `propose` as-is. (Each `ask` costs ~2-6 min
of login-node GP fitting at 24 observations — `cl_min` refits per point.)

`mmi_4x4_90hybrid_bo.py` changed, 2026-09-01:
- `MIN_YGAP_SEP = 0.03` — a candidate within 0.03 µm of an evaluated point at
  the **same** width is rejected, not just an exact duplicate. 3× the npz
  naming precision, and well inside the objective differences the table shows
  at 0.01 µm spacing (y_gap 0.50 vs 0.51: 1.0047 vs 1.0060).
- `ASK_OVERSAMPLE = 4` — rejections are refilled from the GP's own ranking
  instead of shrinking the round.
- `explore_fill()` — farthest-point fallback in (width-index, y_gap-frac)
  space when the GP cannot supply distinct candidates, so a round always
  submits its full batch.
- `propose --dry-run` — prints the round's candidates and submits nothing.
  Use it before every round now; the GP's collapse was invisible until a job
  had already been submitted.

**But the fallback's own answer is not good enough either, and this is the
load-bearing finding.** With the GP exhausted, `explore_fill` proposes
`(5.5, 0.85)`, `(6.0, 1.10)`, `(5.0, 0.65)`, `(6.5, 1.15)` — all four in
regions the table already shows are bad (large y_gap, or the two worst
widths). Pure exploration would spend 4 × ~4 h on known-bad territory.

### The real gap the GP cannot see: W=7.0's interior is unmeasured

`W=7.0` has been sampled at `y_gap` 0.21, 1.11, 1.36 — `taper_width` 1.540,
0.640, 0.390. **Every other width's optimum sits at `taper_width` ≈
0.85-1.12; at W=7.0 that band is `y_gap` ≈ 0.63-0.90 — entirely unmeasured.**
The GP will not ask for it because it works in `(mm_width, y_gap_frac)` space,
where the good ridge is diagonal, and because EI has already collapsed onto
W=6.00.

Two things follow, both cheap (W=7.0's body is cached from Phase A):

1. **Highest-information points available: `W=7.0, y_gap≈0.65`
   (`taper_width=1.10`) and `W=7.0, y_gap≈0.40` (`taper_width=1.35`)** — they
   decide whether the monotone 5.0 → 6.5 improvement continues past 6.5 or
   turns over, which is the whole question of whether the width lattice needs
   extending to 7.5/8.0 (each new width = one new ~1-2 h body solve).
2. **Reparameterize the GP to `(mm_width, taper_width)`** instead of
   `(mm_width, y_gap_frac)`. `y_gap = mm_width/4 − taper_width` is exact, so
   **all 24 existing points convert for free, zero cluster time** — and the
   ridge becomes axis-aligned, which an isotropic kernel can actually learn.
   This is a coordinate change on the search only; the device, the evaluator
   and the objective are untouched.

### Open, needs a decision (2026-09-01)

- **Job 12674307 (`beef_bo_W6.00_y0.52`, 4 nodes) is in flight and
  redundant** — it is the one job round 6 managed to submit, 0.01-0.02 µm from
  two points already measured at objective ≈1.006, and 0.33 in objective worse
  than the best point. `scancel` on it was blocked by this session's tool
  permissions; it needs to be cancelled by hand. Until it clears, `propose`
  refuses to open a new round by design (one round outstanding at a time).
- **The CMRR wall.** No point in the 2D space gets past ~12 dBe on the
  good-phase ridge (spec: 20 dBe). Options, none yet chosen: extend the width
  lattice upward (7.5/8.0, one body solve each); free a third parameter
  (`taper_length`, `waveguide_width`, or decoupling `taper_width` from the
  paper's Eq. 8 relation entirely); or accept the measured frontier and report
  the device as phase-compliant / CMRR-limited.

### ✅ Reparameterization landed + round 7 submitted (2026-09-01)

`mmi_4x4_90hybrid_bo.py` now searches **`(mm_width, taper_width)`** —
`Categorical(MM_WIDTHS)` × `Real(TAPER_MIN=0.36, TAPER_MAX=1.65)`, with
`y_gap = mm_width/4 − taper_width`. Exact coordinate change, all 24 results
converted with no new simulation. Immediate effect, visible in the first
`--dry-run`: the GP's own top candidate became `W=7.00, taper=1.060`, i.e. it
now asks for the ridge it could not see in `y_gap`-fraction space.

Two further fixes found by that dry run:
- **`TAPER_CLIP_SLACK = 0.02`, out-of-range candidates rejected not clipped.**
  The `Real` dimension spans the *widest* width's taper range, so at
  W=5.0-6.0 a large fraction of it is invalid; clipping collapsed those
  candidates onto the same `y_gap = 0.10` floor point (3 of 4 in one round).
- `propose --seed W,YGAP` (repeatable) — hand-picked points enter the round
  ahead of the GP's candidates, for exactly the case where a documented
  physics read beats an exhausted acquisition.

**Round 7, submitted 2026-09-01, all four started immediately (4 nodes each,
`gts-sr108-gedc_main`; W=6.5 and W=7.0 bodies both cached from Phase A, so
~3-3.7 h/point expected):**

| job | mm_width | y_gap | taper_width | what it decides |
|---|---|---|---|---|
| 12677821 | 7.0 | 0.65 | 1.100 | the ridge point at W=7.0 — same taper as W=6.0's 1.100 (obj 1.188) and near W=6.5's best 1.065/1.115 (obj 0.674/0.677). **Does the monotone 5.0 → 6.5 improvement continue past 6.5, or turn over?** |
| 12677822 | 7.0 | 0.40 | 1.350 | the ridge's wide-taper side at W=7.0, against W=6.5's 1.305 (obj 0.856) |
| 12677823 | 7.0 | 0.85 | 0.900 | the narrow-taper side, connecting to W=7.0's already-measured 0.640 (obj 2.117) — brackets 7.0's interior optimum from both sides |
| 12677824 | 6.5 | 0.44 | 1.185 | fills the one real gap at the best width, between taper 1.115 (obj 0.677) and 1.305 (obj 0.856) |

If W=7.0 at taper ≈ 1.10 beats 6.5's 0.674, the width lattice needs extending
to 7.5/8.0 (one ~1-2 h Phase A body solve per new width). If it does not,
`mm_width = 6.5` is the optimum on this lattice and the CMRR wall is the only
remaining question — see the open decision above.

### ✅ ROUND 7 LANDED — the width trend TURNS OVER at 6.5 (2026-09-01)

All four jobs COMPLETED, **04:02:31 - 04:07:20** each (4 nodes). Per-section
timing from the logs: ~3620-3720 s × 4 input sections, i.e. **the body was NOT
a cache hit at these points** — the 3-3.7 h/point estimate assumed it was, so
budget ~4.1 h/point for W=6.5-7.0 from now on.

| job | W | y_gap | taper | max\|Δφ\| | min CMRR | obj |
|---|---|---|---|---|---|---|
| 12677821 | 7.0 | 0.65 | 1.100 | 6.534° | **13.587** | 0.8138 |
| 12677822 | 7.0 | 0.40 | 1.350 | 6.552° | 8.909 | 0.9325 |
| 12677823 | 7.0 | 0.85 | 0.900 | 10.047° | 10.928 | 1.2315 |
| 12677824 | 6.5 | 0.44 | 1.185 | 4.980° | 10.895 | 0.7256 |

**The question round 7 was built to answer is answered: the objective turns
over at `mm_width = 6.5`.** W=7.0's own ridge point (taper 1.100, the same
taper band where 6.0 and 6.5 optimize) lands at **0.8138 vs 6.5's 0.6737** —
20% worse, and worse than three separate W=6.5 points. The monotone 5.0 → 6.5
improvement (2.14 → 1.66 → 1.00 → 0.674) does not continue.

**Consequence: do NOT extend the width lattice to 7.5/8.0.** That was the
contingency in round 7's own table; it is now ruled out by measurement rather
than left open, and the two Phase A body solves it would have needed are not
worth spending.

W=7.0's interior is now bracketed on both sides (taper 0.900 → 1.2315,
1.100 → 0.8138, 1.350 → 0.9325), so 7.0's own best is ≈0.81 near taper 1.10 —
still above 6.5's. W=6.5's optimum is bracketed to **taper 1.065-1.115**
(y_gap 0.51-0.56): the new taper-1.185 point is 0.7256, worse than both.

**New physics datum, and the most useful thing round 7 produced beyond the
turn-over: CMRR on the good-phase ridge is still climbing with width.**
Best min-CMRR at each width's ridge point: 6.0 → 11.576, 6.5 → 12.109,
7.0 → **13.587** dBe. Phase moves the opposite way (9.771° → 4.764° →
6.534°). So the two spec metrics are *anti-correlated in width* along the
ridge, and the 20 dBe CMRR bar is not reachable by widening — 7.0 buys
1.5 dBe of CMRR for 1.8° of phase, and the objective gets worse. **The CMRR
wall is structural in this 2-parameter family, not a search failure.** That
strengthens the case for freeing a third parameter over any further
(mm_width, y_gap) search.

**28 points evaluated. Best unchanged: `mm_width=6.50`, `y_gap=0.56`
(taper 1.065), `mm_length=91.0` — max|Δφ| 4.764°, min CMRR 12.109 dBe,
objective 0.6737.**

### Round 8 — `taper_length` opened as a third axis (2026-09-01)

**Why no round 8 in `(mm_width, y_gap)`, decided from the 28 points at zero
cluster cost:**
- The basin is flat: the top two points are 0.4% apart in objective (0.6737
  vs 0.6767) and 0.05 µm apart in `y_gap`. Not resolvable at ~4.1 h/point.
- **Re-scored every point with phase as a hard CONSTRAINT** —
  `0.5·max(0, (max|Δφ| − 5)/5) + 0.5·max(0, (20 − CMRR)/20)`, giving no credit
  for driving phase below spec, unlike the live objective's graded phase term
  — and the ranking is **identical**, same top three in the same order
  (W=6.5 at taper 1.065 / 1.115 / 1.185). The objective's shape is not hiding
  a better region, so there is no "fix the objective and re-search" move.
- Only **3 of 28** points meet the phase spec, all at W=6.5, and their CMRR
  spans just 10.895-12.109 dBe. Best CMRR anywhere is 14.516 dBe at
  (5.0, 0.10) — whose phase is 20.058°, 4× over spec.

**The lever that aims at the failing metric.** The taper is 0.35 µm →
`taper_width` over `TAPER_LENGTH`, and `TAPER_LENGTH` has been **1.0 µm** the
whole campaign (the `MMI` class default, Step 2's placeholder). At the settled
`taper_width=1.065` that is a 0.7 µm width change in 1.0 µm of length — steep
enough to be non-adiabatic, which puts power into unwanted body modes. That
shows up as the 6.5-8.4 dB excess loss and as output imbalance, and imbalance
is exactly what caps CMRR.

**⚠ This axis was closed once, on 2026-08-08, by the user's direct call — and
that close does not survive the method change.** The null result came from
`mmi_4x4_90hybrid_taper_sweep.py`, i.e. **pure EME**, the path this campaign
later disqualified for this device (no port crosstalk; see the bias study
above). Its "flat" answer was identical to 4 displayed digits between
`taper_length` 1.0 and 3.0 µm — a 3× geometry change producing no digit of
difference reads as solver insensitivity, not physics. Under this evaluator's
FDTD input sections, with all four ports present, `taper_length` has never
been measured. Reopened on that basis, not because the earlier decision was
wrong on its own evidence.

Code changes, both minimal and backward-compatible:
- `mmi_4x4_90hybrid_wg_eval.py`: `TAPER_LENGTH` now reads
  `MMI4X4_TAPER_LENGTH` (default 1.0, so every existing point is unchanged),
  is echoed in the header log line, and is stored in the npz. The output
  filename keeps the **exact** historical form at the default and gains a
  `_tl<len>` tag otherwise.
- `mmi_4x4_90hybrid_bo.py`: `absorb_disk()` now requires the filename to
  equal `eval_W{w}_y{y}.npz` exactly. Without that guard a `_tl` file would
  have parsed as the default-taper point at the same (W, y_gap) and
  overwritten its objective with a different geometry's — silently.

**Round 8, submitted 2026-09-01, 4 nodes each, `-t 720`:**

| job | W | y_gap | taper_length | what it decides |
|---|---|---|---|---|
| 12699318 | 6.5 | 0.56 | 2.0 | best geometry, 2× taper |
| 12699319 | 6.5 | 0.56 | 4.0 | 4× |
| 12699320 | 6.5 | 0.56 | 8.0 | 8× — where adiabaticity should show if it is the mechanism |
| 12699321 | 7.0 | 0.65 | 4.0 | does the lever's sign hold at the width with the best measured CMRR (13.587 dBe)? |

`taper_length=1.0` at (6.5, 0.56) is the existing 0.6737 point, so no job
spends time re-measuring the baseline. **If this comes back null under FDTD
too, the honest conclusion is that the device is phase-compliant and
CMRR-limited in this topology, and the campaign closes** — there is no further
2-parameter search worth running.

#### ✅ FIRST ROUND-8 RESULT — `taper_length` is REAL, not null (job 12699318, 2026-09-01)

`W=6.5, y_gap=0.56, taper_length=2.0` — COMPLETED 04:06:17.

| | taper_length=1.0 (baseline) | taper_length=2.0 | change |
|---|---|---|---|
| max\|Δφ\| | 4.764° | **4.154°** | −0.61° (−12.8%) |
| min CMRR | 12.109 dBe | **13.021 dBe** | **+0.91 dBe** |
| objective | 0.6737 | **0.5898** | −12.5% |
| best mm_length | 91.00 | 91.50 | +0.5 |

**Both metrics improved at once.** That is the important part, and it is
qualitatively different from the `mm_width` axis, where phase and CMRR are
anti-correlated (round 7: 7.0 bought +1.5 dBe of CMRR for +1.8° of phase and a
worse objective). `taper_length` moves them together — consistent with the
non-adiabaticity hypothesis, where a steeper taper injects power into unwanted
body modes and degrades *both* the phase relationship and the output balance.
Best objective on the campaign is now **0.5898**, beating the 28-point
(mm_width, y_gap) optimum.

**This directly contradicts the 2026-08-08 pure-EME null** (`taper_sweep.py`:
"still exactly flat", identical to 4 displayed digits between tl=1.0 and 3.0).
Under FDTD input sections with all four ports present, doubling the taper is
worth 12.5% of the objective. The earlier flatness was the disqualified
solver's insensitivity, exactly as suspected when the axis was reopened — the
close was sound on its own evidence and wrong on the physics.

Cost note: 4 FDTD input sections in 245.6 min (~3680 s each, same as tl=1.0 —
taper length does not change per-section cost), and the **body solve + 7601-
length sweep took 0.2 min**, i.e. a clean cache hit. The body depends on
`mm_width` only, and that holds across `taper_length` too.

Three jobs still running: tl=4.0 and tl=8.0 at the same geometry, and tl=4.0
at (7.0, 0.65). **The open question is whether this keeps improving with
length or turns over** — 2.0 is only the first step, and a monotone trend to
8.0 would mean the taper should be longer still.

#### tl=4.0 — the lever SATURATES fast (job 12699319, COMPLETED 04:20:12)

| taper_length | max\|Δφ\| | min CMRR | objective | Δobj vs previous |
|---|---|---|---|---|
| 1.0 (baseline) | 4.764° | 12.109 dBe | 0.6737 | — |
| 2.0 | 4.154° | 13.021 dBe | 0.5898 | **−0.0839** |
| 4.0 | **4.133°** | **13.280 dBe** | **0.5813** | −0.0085 |

**Doubling 1→2 bought 10× what doubling 2→4 bought.** Phase is essentially
converged (4.154 → 4.133°, a 0.02° move), CMRR is still creeping (+0.26 dBe).
The shape is a saturating curve, not a linear trend — the taper is already
close to adiabatic by 2 µm, and the residual is something else.

**Consequence for the CMRR wall: this lever does not open it.** Extrapolating
the CMRR sequence (12.109 → 13.021 → 13.280, steps of +0.91 then +0.26) gives
an asymptote near **13.4-13.5 dBe**, nowhere near the 20 dBe bar. Round 8 has
found a real ~13% objective improvement and a better device, but it has not
changed the campaign's headline conclusion.

`mm_length` keeps drifting with taper length (91.00 → 91.50 → 91.80), as
expected — a longer taper moves the effective self-imaging plane, and the
per-point 81-length sweep re-finds it each time. Body solve was again a cache
hit (0.3 min).

tl=8.0 (12699320) will say whether the curve is flat or still creeping past 4.

#### ⚠ W=7.0 at tl=4.0 — the width ordering FLIPS, round 7's turn-over was an artifact of the short taper (job 12699321, COMPLETED 04:23:35)

`W=7.0, y_gap=0.65, taper_length=4.0` → **max|Δφ| = 3.732°, min CMRR = 13.672
dBe, objective = 0.5314** — the best point of the entire campaign, beating
W=6.5 at the same taper length (0.5813).

| geometry | tl=1.0 | tl=4.0 | Δobj from the taper |
|---|---|---|---|
| W=6.5, y=0.56 | 0.6737 | 0.5813 | −0.092 |
| W=7.0, y=0.65 | 0.8138 | **0.5314** | **−0.282** |

**The two axes interact, and round 7's headline is now void.** Round 7
concluded "the objective turns over at mm_width = 6.5" and, on that basis,
ruled out extending the width lattice to 7.5/8.0. That conclusion was measured
entirely at `taper_length = 1.0`. At tl=4.0 the ordering reverses: 7.0 beats
6.5. The turn-over was not a property of the width axis — it was the wider
device being hurt more by a non-adiabatic taper, which is exactly what one
would expect, since a wider body supports more modes for a steep taper to
scatter into.

Per-metric, the mechanism is clear: **W=7.0's gain is almost entirely phase**
(6.534° → 3.732°, −2.8°), while its CMRR barely moves (13.587 → 13.672 dBe).
W=6.5's phase was already near-converged at 4.7° and had little to give.

**Two consequences:**
1. **`mm_width` 7.5 / 8.0 must be re-opened**, at tl≈4.0, not at tl=1.0. Each
   new width costs one Phase A body solve (~1-2 h) plus ~4.4 h/point. The
   "ruled out by measurement" note in the round 7 section above applies ONLY
   at tl=1.0 and should not be cited without that qualifier.
2. **Any conclusion in this document drawn from the 28-point (mm_width, y_gap)
   campaign is conditional on `taper_length = 1.0`** — including the shape of
   the ridge, the taper_width optimum, and the width trend. The 2D surface was
   mapped with a systematically non-adiabatic taper.

**The CMRR wall survives all of it.** Best CMRR anywhere is now 13.672 dBe
against a 20 dBe spec, and the taper lever moved it +0.09 dBe at W=7.0. Phase
is comfortably met (3.732° vs 5°); CMRR is not close, and nothing measured so
far — width, y_gap, or taper length — moves it more than ~1.5 dBe.

#### ✅ ROUND 8 COMPLETE — tl=8.0 landed (job 12699320, 04:35:12), and the saturation read was WRONG

Full `taper_length` series at `W=6.5, y_gap=0.56`:

| taper_length | max\|Δφ\| | min CMRR | objective | ΔCMRR per doubling |
|---|---|---|---|---|
| 1.0 | 4.764° | 12.109 dBe | 0.6737 | — |
| 2.0 | 4.154° | 13.021 dBe | 0.5898 | +0.91 |
| 4.0 | 4.133° | 13.280 dBe | 0.5813 | +0.26 |
| 8.0 | 4.159° | **13.615 dBe** | **0.5755** | +0.34 |

**Correction to the tl=4.0 entry above: CMRR is not saturating to 13.4-13.5.**
That extrapolation was made from two steps (+0.91, +0.26) and the third step
came in *larger* (+0.34), not smaller. The real behaviour is two different
curves in one:
- **Phase saturates hard.** 4.154° → 4.133° → 4.159° from tl=2 to tl=8 — a
  0.026° spread, at or below this evaluator's own noise. The taper is
  adiabatic for phase purposes by 2 µm, full stop.
- **CMRR does not saturate; it creeps ~+0.3 dBe per doubling** and shows no
  sign of stopping over 3 doublings. That is a logarithmic trend, not an
  approach to an asymptote.

**The logarithmic creep makes the CMRR verdict stronger, not weaker.** At
+0.3 dBe per doubling from 13.6, reaching the 20 dBe bar needs of order
**twenty more doublings** — a taper millions of µm long. `taper_length` cannot
deliver the spec even in principle, and the honest statement is now
quantitative rather than a guess at an asymptote.

**Round 8's verdict: the lever is REAL but it is a phase/loss lever, not a
CMRR lever.** It bought a genuine 12-15% objective improvement at W=6.5, a 35%
improvement at W=7.0, overturned round 7's width conclusion, and is worth
keeping in the design (tl=4.0 is the sensible setting — tl=8.0 buys 0.006 of
objective for 2× the taper footprint). It does not move the device toward
CMRR compliance.

### Campaign state after round 8 (2026-09-01)

**Best geometry measured: `mm_width=7.0`, `y_gap=0.65` (taper_width=1.100),
`taper_length=4.0`, `mm_length=106.85` — max|Δφ| = 3.732°, min CMRR = 13.672
dBe, objective = 0.5314.** Phase spec MET with 25% margin. CMRR spec MISSED by
6.3 dBe.

32 evaluated points (28 at tl=1.0, 4 at tl>1.0). Round 8 cost 4 jobs at
04:06-04:35 each; per-section FDTD cost is independent of taper length
(~3680-4130 s), and the body remained a cache hit throughout (0.2-0.3 min),
confirming body independence from `taper_length` as well as from `y_gap`.

**The two open items, in priority order:**
1. **`mm_width` 7.5 / 8.0 at tl=4.0** — round 7 ruled these out at tl=1.0 and
   that ruling does not transfer (see the flip above). This is the only
   remaining unexplored direction with evidence behind it: W=7.0 improved by
   0.282 when given an adiabatic taper, and nothing says 7.5 will not.
   Cost: one Phase A body solve per width (~1-2 h) + ~4.4 h/point.
2. **The CMRR wall, now characterized rather than merely observed.** Three
   independent levers have been measured against it: `mm_width` (±1.5 dBe,
   anti-correlated with phase), `y_gap` (±1.5 dBe across the whole valid
   range), `taper_length` (+0.3 dBe per doubling, logarithmic). None reaches
   20 dBe. If the spec is to be met, it needs a change of kind — a different
   port/aperture topology, or a re-read of what CMRR bar this device actually
   has to hit — not another parameter sweep.

### Round 9 — re-opening `mm_width` 7.5 / 8.0 at an adiabatic taper (2026-09-02)

Item 1 from the round-8 close-out. Round 7 ruled these widths out; that ruling
was measured at `taper_length=1.0` and does not transfer (see the tl=4.0 flip).

**Phase A first — jobs 12735983_0 (W=7.5) and _1 (W=8.0), 1 node/8 ranks each,
RUNNING from submit.** Body-only solves, so the four evaluation points that
follow are cache hits on the body rather than each paying its own ~1-2 h solve
— and, just as important, two points at the same width do not become two
concurrent writers on one set of cache keys.

`run_mmi_4x4_90hybrid_ygap_bodies.batch` now takes its width list from
`MMI4X4_WIDTHS` (default: the original Phase A five, so `--array=0-4` is
unchanged), and echoes the width per array task so the mapping is in the log
rather than inferred. No `taper_length` is set — round 8 measured the body to
be independent of it (tl=2/4/8 all hit the tl=1-era cache in 0.2-0.3 min).

**Planned Phase B points, all at `taper_length=4.0`** (round 8's sensible
setting: tl=8.0 bought only 0.006 of objective for 2× the footprint):

| W | taper_width | y_gap | why |
|---|---|---|---|
| 7.5 | 1.10 | 0.775 | the constant-taper line through 6.5/7.0 — the direct width comparison |
| 7.5 | 1.20 | 0.675 | the measured optimum taper drifts up with width (6.0 ≈ 0.85-1.00, 6.5 ≈ 1.065-1.115, 7.0 ≈ 1.10), so probe above the line too |
| 8.0 | 1.10 | 0.900 | constant-taper line at the far end |
| 8.0 | 1.20 | 0.800 | shifted-taper probe at the far end |

**⚠ Methodological caveat to settle before these numbers are load-bearing.**
Every width in this campaign is evaluated at `N_MODES = 16`, and that basis
size was converged at `mm_width = 6.0` (see "Body/taper mode-count convergence
investigation", 2026-08-09). The Step-1c width scope measured 24 guided modes
at every width from 4.0 to 10.0 µm — the `n_request=24` cap — so a fixed
16-mode basis is not obviously adequate at 8.0 just because it was at 6.0. If
a wide point wins, re-run it at `N_MODES=20` before treating the win as real;
if the wide points lose, the caveat does not matter, since a truncated basis
would if anything flatter the wider device's self-imaging, not penalise it.

Expect >4.4 h/point at these widths: `sim_width = mm_width + 3`, so the FDTD
domain grows with the device (11 µm at W=8.0 vs 9.5 at 6.5). `-t 720` has
ample headroom.

#### ✅ Round 9 Phase A DONE, Phase B submitted (2026-09-02)

Both body solves COMPLETED: **W=8.0 in 1:05:07** (body + 7601-length sweep in
64.9 min), **W=7.5 in 1:10:09** (70.0 min) — 1 node/8 ranks each, in line with
the original Phase A's 1.0-2.25 h range. `L_pi = 193.635 µm` at W=8.0, centre
145.2 µm, comfortably inside the shared 20-400 µm length grid, so no grid
widening is needed at these widths.

The two W=8.0 points were submitted as soon as *its* body landed rather than
waiting for W=7.5 — the bodies are per-width cache keys, so the still-running
W=7.5 solve could not collide with them, and it saved an hour of wall clock.

| job | W | y_gap | taper_width |
|---|---|---|---|
| 12737440 | 8.0 | 0.900 | 1.10 |
| 12737441 | 8.0 | 0.800 | 1.20 |
| 12737705 | 7.5 | 0.775 | 1.10 |
| 12737706 | 7.5 | 0.675 | 1.20 |

All at `taper_length=4.0`, 4 nodes/96 ranks, `-t 720`. Baselines to beat:
**W=7.0 y=0.65 tl=4.0 → 3.732° / 13.672 dBe / 0.5314** (campaign best) and
W=6.5 y=0.56 tl=4.0 → 4.133° / 13.280 dBe / 0.5813.

#### ⚠⚠ W=8.0 BREAKS THE CMRR WALL — 17.664 dBe (job 12737440, COMPLETED 04:52:22)

`W=8.0, y_gap=0.900 (taper 1.10), taper_length=4.0, mm_length=140.75` →
**max|Δφ| = 6.794°, min CMRR = 17.664 dBe, objective = 0.7378.**

CMRR by width along the constant-taper≈1.10 line, all at tl=4.0:

| W | taper | max\|Δφ\| | min CMRR | objective |
|---|---|---|---|---|
| 6.5 | 1.065 | 4.133° | 13.280 dBe | 0.5813 |
| 7.0 | 1.100 | 3.732° | 13.672 dBe | 0.5314 |
| 8.0 | 1.100 | 6.794° | **17.664 dBe** | 0.7378 |

**Correction, and it is a big one: my "the CMRR wall is structural, no lever
moves it more than ~1.5 dBe" conclusion was wrong.** It was generalized from
widths ≤ 7.0, where CMRR really does sit at 12-13.7 dBe. Going to 8.0 moves it
**+4.0 dBe in one step** — 17.664 dBe, within 2.3 dBe of the 20 dBe spec, and
by far the best CMRR this device has ever produced. The earlier statement
should be read as "no lever moves CMRR within 5.0 ≤ mm_width ≤ 7.0", which is
what was actually measured.

**Why the objective still looks worse (0.7378 vs 0.5314):** phase degraded to
6.794°, over the 5° bar, and the live objective grades phase linearly with no
credit for meeting it. Under a phase-as-constraint scoring this point is
`0.5·(6.794−5)/5 + 0.5·(20−17.664)/20 = 0.238`, which would make it the best
point of the campaign by a wide margin, versus 0.197 for the old W=6.5
optimum... i.e. **the two scorings now disagree about which geometry wins.**
They agreed on all 28 earlier points; they do not agree here. That has to be
resolved before anything is called "the" optimum — see below.

**This is not a converged result yet, and the pre-registered caveat now
matters.** Round 9's plan said: *"if a wide point wins, re-run it at
`N_MODES=20` before treating the win as real."* A wide point has now won on
the metric that was the whole problem, at `N_MODES = 16` — a basis converged at
`mm_width = 6.0`, against a cross-section the Step-1c scope measured at 24
guided modes. **The N_MODES=20 re-run is now required, not optional.**

Cost note: 4 input sections in 291.7 min (~4380 s each vs ~3680 s at W=6.5-7.0
— the `sim_width = mm_width + 3` domain growth, +19%/section as expected).
`mm_length = 140.75 µm` at this width, still well inside the 20-400 grid.

Three jobs still running, and two of them are now much more interesting than
when they were submitted: `W=8.0 y=0.800` (taper 1.20 — does the wider taper
recover the phase at this width?) and the two W=7.5 points, which bracket
where the CMRR climb begins.

#### ✅✅ W=8.0 AT taper 1.20 — BEST POINT OF THE CAMPAIGN, phase MET and CMRR 17.4 dBe (job 12737441, COMPLETED 04:53:23)

`W=8.0, y_gap=0.800 (taper 1.20), taper_length=4.0, mm_length=140.70` →
**max|Δφ| = 4.306° (MEETS the 5° spec), min CMRR = 17.393 dBe,
objective = 0.4957.**

The wider taper recovered the phase at this width without giving back the CMRR
gain, which is exactly what the "optimum taper drifts up with width" trend
predicted (6.0 ≈ 0.85-1.00, 6.5 ≈ 1.065, 7.0 ≈ 1.10, **8.0 ≈ 1.20**):

| W | taper | max\|Δφ\| | min CMRR | obj (live) | obj (phase-as-constraint) |
|---|---|---|---|---|---|
| 6.5 | 1.065 | 4.133° | 13.280 | 0.5813 | 0.1973 |
| 7.0 | 1.100 | 3.732° | 13.672 | 0.5314 | 0.1582 |
| 8.0 | 1.100 | 6.794° | 17.664 | 0.7378 | 0.2379 |
| **8.0** | **1.200** | **4.306°** | **17.393** | **0.4957** | **0.0652** |

**The scoring disagreement flagged one entry above is resolved: this point
dominates under BOTH scorings** — best live objective (0.4957 vs 0.5314) and
best constraint objective (0.0652 vs 0.1582) in the campaign. No trade-off to
adjudicate.

**Where the device now stands: phase MET at 4.306°, CMRR 2.6 dBe short of the
20 dBe bar** — against 12.109 dBe when this campaign's optimizer started, and
against the "structurally capped ~13.6" I wrote yesterday. The campaign moved
CMRR +5.3 dBe while keeping phase compliant.

**Pre-registered confirmation now in flight: job 12748566** — body-only solve
at `W=8.0, N_MODES=20`, 1 node/8 ranks. The evaluation point at N_MODES=20
follows once that body lands. This is the round-9 plan's own gate on a wide
win, and until it clears, the 17.4 dBe number is provisional.

*(Note for anyone re-running this: `run_mmi_4x4_90hybrid_ygap_bodies.batch`
takes its width from `MMI4X4_WIDTHS`, and the in-script `WIDTHS` assignment
OVERRIDES a plain `MMI4X4_W` export — a first attempt at this job with
`MMI4X4_W=8.0` would silently have solved W=5.0 instead. Caught before it ran;
cancelled and resubmitted with `MMI4X4_WIDTHS="8.0"`, and the per-task echo
in the log is what makes this checkable.)*

#### W=7.5 at taper 1.10 — the CMRR climb is monotone AND accelerating (job 12737705, COMPLETED 04:58:40)

`W=7.5, y_gap=0.775 (taper 1.10), taper_length=4.0, mm_length=123.45` →
max|Δφ| = 5.229°, min CMRR = 15.251 dBe, objective = 0.6417.

The constant-taper-1.10 line is now four points, all at tl=4.0:

| W | max\|Δφ\| | min CMRR | ΔCMRR per +0.5 µm |
|---|---|---|---|
| 6.5 | 4.133° | 13.280 | — |
| 7.0 | 3.732° | 13.672 | +0.39 |
| 7.5 | 5.229° | 15.251 | +1.58 |
| 8.0 | 6.794° | 17.664 | +2.41 |

**CMRR rises monotonically with width and the steps are getting BIGGER, not
smaller** (+0.39, +1.58, +2.41 dBe). Nothing in this sequence looks like an
approach to a ceiling — which is the opposite of what the 5.0-7.0 range
suggested, and the reason the "structural wall" reading was wrong.

Phase along the same line degrades past 7.0 (3.732° → 5.229° → 6.794°), but
that is recoverable: at W=8.0 moving the taper from 1.10 to 1.20 brought phase
back to 4.306° while holding CMRR at 17.393 dBe. The per-width recipe is
"widen the taper along with the body," and the un-recovered taper-1.10 points
should be read as the *floor* of what each width can do on phase, not its best.

**Extrapolating the trend, CMRR ≥ 20 dBe looks reachable near W ≈ 9** — with
the same caveat as everything else here: it needs measuring, not extrapolating,
and it needs the `N_MODES` question settled first, since basis adequacy is
exactly what a wider body strains.

*Minor bookkeeping wart, noted so it does not confuse a later reader: the npz
for this point is named `eval_W7.50_y0.78_tl4.00.npz` — the filename rounds
y_gap to 2 decimals, while the value actually simulated (and stored inside the
npz) is 0.775. The BO driver's exact-filename guard means such a point is not
absorbed into `bo_state.json` as a y=0.78 result, which is the safe failure
mode, but the name should not be read as the geometry.*

#### ✅ ROUND 9 PHASE B COMPLETE — and the two scorings now genuinely disagree (job 12737706, COMPLETED 04:59:49)

`W=7.5, y_gap=0.675 (taper 1.20), taper_length=4.0, mm_length=123.15` →
**max|Δφ| = 3.261° (best phase of the entire campaign), min CMRR = 15.495 dBe,
objective = 0.4387 (best LIVE objective of the campaign).**

All four round-9 points, tl=4.0:

| W | taper | max\|Δφ\| | min CMRR | obj (live) | obj (constraint) |
|---|---|---|---|---|---|
| 7.5 | 1.10 | 5.229° | 15.251 | 0.6417 | 0.1417 |
| 7.5 | 1.20 | **3.261°** | 15.495 | **0.4387** | 0.1126 |
| 8.0 | 1.10 | 6.794° | 17.664 | 0.7378 | 0.2379 |
| 8.0 | 1.20 | 4.306° | **17.393** | 0.4957 | **0.0652** |

Two clean sub-results first:
- **Taper 1.20 beats taper 1.10 at both new widths**, confirming the
  taper-drifts-up-with-width recipe as a rule rather than a one-off at 8.0.
- **Round 7's width turn-over is now fully dead.** With an adiabatic taper and
  a width-matched taper_width, the objective improves monotonically
  6.5 → 7.0 → 7.5, and 8.0 is better than everything on CMRR.

**The scoring disagreement is now decision-relevant, and it needs settling.**
The live objective ranks W=7.5/1.20 first (0.4387); phase-as-constraint ranks
W=8.0/1.20 first (0.0652). The disagreement is real, not cosmetic: the live
objective is still paying for phase improvements *below* the 5° spec — at
3.261° it banks 0.326 of credit for margin the device cannot use — while
treating a 4.5 dBe CMRR shortfall as the smaller problem.

**Recommendation, and how this document will score from here: phase as a
constraint, CMRR as the objective.** Both candidates meet the phase spec
(3.261° and 4.306° against a 5° bar), so the only difference that reaches the
finished device is CMRR, where 17.393 beats 15.495 by 1.9 dBe. On that basis
**`W=8.0, y_gap=0.800, taper_width=1.20, taper_length=4.0, mm_length=140.70`
is the best geometry this campaign has produced.** The earlier check that
found the two scorings agreeing was run on the 28 tl=1.0 points, where nothing
met the phase spec and the question could not arise; that check should not be
cited as evidence they agree in general.

**Round 9 answered item 1 in full: the width lattice DID need extending.**
Round 7's "ruled out by measurement" was an artifact of `taper_length = 1.0`,
and the widths it excluded contain the best points found so far — including
the only ones above 15 dBe CMRR.

**Still open, in flight:** the `N_MODES = 20` confirmation at W=8.0
(body 12748566 running). Until it lands, every number in the table above
carries the 16-mode-basis caveat, and the wide points carry it most.

#### N_MODES=20 confirmation in flight — and a silent-overwrite bug caught on the way (2026-09-02)

`W=8.0, N_MODES=20` body-only COMPLETED in **1:59:52** (119.4 min vs 64.9 min
for the same body at N_MODES=16 — the larger basis costs ~1.8×). Evaluation
point submitted as **job 12752148** (`W=8.0, y_gap=0.800, taper_length=4.0,
N_MODES=20`).

**Bug caught before it did damage: `wg_eval.py`'s npz filename encoded
`mm_width`, `y_gap` and `taper_length` but NOT `N_MODES`.** The confirmation
run would therefore have written to `eval_W8.00_y0.80_tl4.00.npz` — silently
overwriting the N_MODES=16 result that it exists to be compared against, and
destroying the very baseline the check needs. Two fixes, both applied while
12752148 was still PENDING:

1. `eval_W8.00_y0.80_tl4.00.npz` copied to
   `eval_W8.00_y0.80_tl4.00_nmodes16.npz.bak` immediately (the running job
   imports its code at start, so a patch alone could not be relied on to
   protect it).
2. `wg_eval.py` now appends `_n{N_MODES}` to the filename whenever
   `N_MODES != 16`, and stores `n_modes` inside the npz. Same reasoning as the
   `_tl` tag: **a basis-convergence re-run is a different measurement of the
   same geometry, not a replacement for it.**

This is the second naming-collision hazard in this campaign (the first was the
`_tl` tag itself, guarded in `mmi_4x4_90hybrid_bo.py`'s `absorb_disk()`). The
pattern to watch: every time an axis is added to this evaluator, its output
name has to grow with it, or an old result gets overwritten by a new one that
is not comparable to it.

#### ✅✅ N_MODES=20 CONFIRMS THE WIDE WIN — basis is converged at W=8.0 (job 12752148, COMPLETED 06:54:36)

Same geometry, same everything, basis 16 → 20:

| | N_MODES=16 | N_MODES=20 | difference |
|---|---|---|---|
| max\|Δφ\| | 4.306° | 4.288° | 0.018° (0.4%) |
| min CMRR | 17.393 dBe | 17.368 dBe | 0.025 dBe (0.14%) |
| objective | 0.4957 | 0.4946 | 0.0011 (0.2%) |
| wall time | 04:53:23 | 06:54:36 | +41% |

**The 16-mode basis is adequate at `mm_width = 8.0`.** Every metric moves by
well under 1%, far below the differences the campaign is ranking geometries
on (the 6.5→8.0 CMRR gain is 4.1 dBe; this check's disagreement is 0.025 dBe,
160× smaller). **The pre-registered gate is cleared and the wide-point result
stands as measured** — 17.4 dBe CMRR with the phase spec met is real, not a
truncation artifact.

That also retires the standing worry that the whole width comparison was
biased by a basis converged only at `mm_width = 6.0`: it now has a direct
check at the far end of the range, where the strain is worst. Points between
6.0 and 8.0 are bracketed by two confirmed widths, and no re-run of the
intermediate widths is needed.

Both npz files survive side by side — `eval_W8.00_y0.80_tl4.00.npz` (n=16) and
`eval_W8.00_y0.80_tl4.00_n20.npz` — which is exactly what the `_n` filename
tag added an hour earlier exists to guarantee.

**Round 9 is COMPLETE.** Final campaign state:

- **Best geometry: `mm_width=8.0`, `y_gap=0.800`, `taper_width=1.20`,
  `taper_length=4.0`, `mm_length=140.70` — max|Δφ| = 4.31° (spec ≤5°, MET),
  min CMRR = 17.39 dBe (spec ≥20 dBe, 2.6 dBe short).** Basis-converged.
- Campaign progress on CMRR: 12.109 dBe (28-point 2D optimum, tl=1.0) →
  13.672 (tl=4.0 at W=7.0) → **17.393** (W=8.0, width-matched taper). Phase
  has stayed inside spec throughout the last three of those.
- Three of this document's own conclusions were overturned by measurement in
  the last two days, each recorded above rather than quietly edited: the
  width turn-over at 6.5 (round 7 — artifact of `taper_length=1.0`), the
  `taper_length`-is-null close of 2026-08-08 (artifact of pure EME), and the
  "CMRR wall is structural at ~13.6 dBe" reading (artifact of never having
  measured above `mm_width = 7.0`).

**The one open direction: `mm_width` 8.5 / 9.0.** The CMRR-vs-width steps are
+0.39, +1.58, +2.41 dBe per +0.5 µm and still growing; a linear read of the
last two steps puts 20 dBe near W ≈ 8.5-9.0. Each new width needs one Phase A
body (~2 h at N_MODES=16) plus ~5 h/point, and taper_width should start near
1.25-1.30 following the width-matched-taper rule.

### ✅ Excess-loss audit — loss does NOT degrade with width (job 12780374, 2026-09-03)

**The gap this closes:** `wg_eval.py` recorded only `max|Δφ|` and `min CMRR`
for all 33 campaign points, while the optimum was being driven from
`mm_width` 6.5 → 8.0 and `mm_length` 91 → 141 µm. Neither recorded metric can
see a device getting lossier as it grows, so the ranking that named W=8.0 best
rested on two of the three numbers that decide whether the device is usable.

`metrics_at()` now also returns `max_excess_loss_db` and `max_imbalance_db`
(max over both driven ports), printed and stored in the npz. **The objective is
deliberately unchanged** — these are recorded, not optimised, so every ranking
already in this document stays comparable. A new `MMI4X4_CACHE_ONLY` flag
threads `cache_only=True` into both the FDTD sections and the body length
sweep, so an already-simulated geometry can be re-scored with **zero new
simulation** and a cache miss raises instead of silently starting a 5 h solve.

**Audit: 7 geometries, all cache hits, COMPLETED in 02:26** — the entire
recovery cost 2.5 minutes of one node, against ~35 h of FDTD already banked.

| W | y_gap | taper | tl | max\|Δφ\| | min CMRR | **max XL** | **max imb** |
|---|---|---|---|---|---|---|---|
| 6.5 | 0.56 | 1.065 | 1.0 | 4.764° | 12.109 | **4.076 dB** | **5.202 dB** |
| 6.5 | 0.56 | 1.065 | 4.0 | 4.133° | 13.280 | 2.940 | 4.535 |
| 7.0 | 0.65 | 1.100 | 4.0 | 3.732° | 13.672 | 3.340 | 3.976 |
| 7.5 | 0.775 | 1.100 | 4.0 | 5.229° | 15.251 | 2.980 | 2.698 |
| 7.5 | 0.675 | 1.200 | 4.0 | 3.261° | 15.495 | 3.036 | 3.117 |
| 8.0 | 0.900 | 1.100 | 4.0 | 6.794° | 17.664 | 3.314 | 2.409 |
| **8.0** | **0.800** | **1.200** | **4.0** | **4.306°** | **17.393** | **2.956** | **2.129** |

**Every phase and CMRR value reproduces the original run exactly** — same
digits, from cache, through a code path that now computes four metrics instead
of two. That is a free end-to-end check on the section cache and on the audit
itself, and it passed. Zero `AUDIT FAILED` lines.

**Three findings:**

1. **Excess loss is FLAT in `mm_width`: 2.94-3.34 dB across 6.5 → 8.0, with no
   trend and no monotonicity** (2.94, 3.34, 2.98/3.04, 3.31/2.96). The concern
   that drove this audit — that the campaign was optimising toward a longer,
   lossier device without measuring it — **is not borne out.** Doubling the
   body length from 91 to 141 µm costs nothing measurable in excess loss.
2. **Imbalance IMPROVES monotonically with width: 4.535 → 3.976 → 2.698 →
   2.129 dB.** This is the mechanism behind the CMRR climb rather than an
   independent result — CMRR is set by the amplitude balance of the paired
   outputs, so imbalance falling and CMRR rising are two views of one effect.
   It is the first *physical* explanation the campaign has for why wider is
   better, as opposed to an empirical trend.
3. **`taper_length` is also a loss lever, retro-confirming its mechanism.**
   At W=6.5, going tl=1.0 → 4.0 cut excess loss **4.076 → 2.940 dB (−1.14 dB)**
   and imbalance 5.202 → 4.535 dB. Round 8 argued the taper was non-adiabatic
   and scattering power into unwanted modes, but could not show it because loss
   was not recorded then. It is now: the lengthened taper recovers over a dB.

**The ranking is unchanged, and is now better supported than before.**
`W=8.0, y_gap=0.800, taper 1.20, tl=4.0` has the best CMRR of any
phase-compliant point, the **lowest imbalance (2.129 dB)** and the
second-lowest excess loss (2.956 dB) of all seven. It wins on every recorded
metric except phase, where W=7.5/1.20's 3.261° beats its 4.306° — and both
meet the 5° spec, so that margin is not usable.

*Absolute caveat, unchanged by this audit: ~3 dB of excess loss is still
substantial (half the power) for a coherent-receiver front end. It is simply
not a function of width, so it is a separate problem from the one this
campaign has been solving, and it should not be read as "fixed."*

### Round 10 — `mm_width` 8.5 / 9.0, chasing the 20 dBe CMRR bar (2026-09-03)

Motivated by round 9's CMRR-vs-width steps (+0.39, +1.58, +2.41 dBe per
+0.5 µm, still growing) and cleared to proceed by the excess-loss audit, which
showed loss flat in width so a wider device is not being bought with loss.

**Phase A done:** W=8.5 COMPLETED 01:11:29 (body + 7601-length sweep in
71.2 min), W=9.0 COMPLETED 01:27:23 (87.2 min) — 1 node/8 ranks each. Body cost
grows with width (64.9 → 71.2 → 87.2 min at 8.0 → 8.5 → 9.0) but stays cheap
next to a ~5 h evaluation point.

**Phase B submitted, all at `taper_length=4.0`,** following the
width-matched-taper rule (measured optima: 6.5 → 1.065, 7.0 → 1.10,
7.5 → 1.20, 8.0 → 1.20), so each width is probed at taper 1.20 and 1.30:

| job | W | y_gap | taper_width | expected mm_length |
|---|---|---|---|---|
| 12783703 | 8.5 | 0.925 | 1.20 | ~164 µm (0.75·L_π, L_π ≈ 218) |
| 12783704 | 8.5 | 0.825 | 1.30 | ~164 µm |
| 12783705 | 9.0 | 1.05 | 1.20 | ~184 µm (L_π ≈ 245) |
| 12783706 | 9.0 | 0.95 | 1.30 | ~184 µm |

Both widths are well inside the validity bound (`y_gap < mm_width/4 − 0.35`,
i.e. < 1.775 at 8.5 and < 1.90 at 9.0) and inside the 20-400 µm length grid.
Expect >5 h/point — the FDTD domain is `mm_width + 3`, so 12 µm wide at W=9.0
against 11 at 8.0, on top of the longer body.

**What decides the outcome:** whether `min CMRR` crosses **20 dBe** while
`max|Δφ|` stays under 5°, with excess loss still near 3 dB. Beating the
current best needs CMRR > 17.393 dBe at a phase-compliant point.

#### W=8.5 at taper 1.20 — best CMRR yet (18.126 dBe), but phase breaks the spec (job 12783703, COMPLETED 05:30:53)

`W=8.5, y_gap=0.925 (taper 1.20), taper_length=4.0, mm_length=159.40` →
**max|Δφ| = 5.438° (FAILS the 5° bar by 0.44°), min CMRR = 18.126 dBe (best of
the campaign), max XL = 2.949 dB, max imbalance = 2.413 dB, objective 0.5906.**

Against the current best (W=8.0, taper 1.20 — 4.306°, 17.393 dBe, XL 2.956,
imb 2.129):
- **CMRR: +0.73 dBe, a new campaign high**, now 1.87 dBe from the 20 dBe spec.
- **Excess loss: FLAT — 2.949 vs 2.956 dB.** The audit's finding holds at the
  new width; loss is still not a function of `mm_width`.
- **Phase: FAILS.** 5.438° vs 4.306°.
- Imbalance ticked up slightly (2.413 vs 2.129 dB), tracking the phase.

**This is the same pattern W=8.0 showed at taper 1.10** — a wider body at a
taper that has not widened with it overshoots on phase while still gaining
CMRR. There, moving taper 1.10 → 1.20 recovered phase (6.794° → 4.306°) at
essentially no CMRR cost. **The companion point `W=8.5, taper 1.30`
(job 12783704) is exactly that probe at this width**, and it is the one that
decides whether 8.5 beats 8.0.

**Caution on the extrapolation that motivated round 10.** Along the *taper
1.20* line, CMRR now reads 15.495 (7.5) → 17.393 (8.0) → 18.126 (8.5): steps
of **+1.90 then +0.73 — decelerating**, unlike the accelerating steps seen
along the taper-1.10 line. If that deceleration is real, 20 dBe is further
away than the linear read suggested, and W=9.0 (already queued) may land near
18.5-19 rather than clearing the bar.

#### ✅✅ W=8.5 at taper 1.30 — NEW BEST ON EVERY METRIC (job 12783704, COMPLETED 05:28:43)

`W=8.5, y_gap=0.825 (taper 1.30), taper_length=4.0, mm_length=159.30` →
**max|Δφ| = 3.933° (MEETS 5°), min CMRR = 18.482 dBe (campaign high),
max XL = 2.706 dB (campaign LOW), max imbalance = 1.972 dB (campaign low),
objective = 0.4312 (best live objective of the campaign).**

vs the previous best (W=8.0, taper 1.20):

| | W=8.0 / 1.20 | **W=8.5 / 1.30** | change |
|---|---|---|---|
| max\|Δφ\| | 4.306° | **3.933°** | −0.373° |
| min CMRR | 17.393 dBe | **18.482 dBe** | **+1.09 dBe** |
| max excess loss | 2.956 dB | **2.706 dB** | −0.25 dB |
| max imbalance | 2.129 dB | **1.972 dB** | −0.16 dB |
| objective (live) | 0.4957 | **0.4312** | better |
| objective (constraint) | 0.0652 | **0.0380** | better |

**No trade-off anywhere: it wins on all four measured quantities and under
both scorings.** The scoring ambiguity that had to be adjudicated at W=7.5 vs
8.0 does not arise here.

**The taper-widening recipe is confirmed a third time.** At this width,
taper 1.20 → 1.30 moved phase 5.438° → 3.933° (recovering the spec) *and*
raised CMRR 18.126 → 18.482 dBe — the same "widen the taper with the body"
pattern seen at W=8.0 (1.10 → 1.20) and W=7.5 (1.10 → 1.20). The per-width
optimum taper now reads: 6.5 → 1.065, 7.0 → 1.10, 7.5 → 1.20, 8.0 → 1.20,
**8.5 → 1.30**.

**Excess loss did not just stay flat — it improved to the campaign minimum**
(2.706 dB), and imbalance to 1.972 dB. The loss audit's conclusion (loss is
not a function of width) survives at 8.5, and the imbalance-falls-with-width
trend that explains the CMRR climb continues: 4.535 → 3.976 → 2.698 → 2.129
→ **1.972 dB**.

**Still 1.52 dBe short of the 20 dBe CMRR spec.** The two W=9.0 points
(12783705 taper 1.20, 12783706 taper 1.30) are running and will say whether
the remaining gap closes at the next width.

#### ✅ ROUND 10 COMPLETE — W=9.0 lands, and the CMRR climb has PLATEAUED (jobs 12783705/12783706, COMPLETED 05:52)

| W | taper | max\|Δφ\| | min CMRR | max XL | max imb | obj (live) |
|---|---|---|---|---|---|---|
| 9.0 | 1.20 | 5.965° (FAILS) | 16.593 | 2.651 | 2.932 | 0.6817 |
| 9.0 | 1.30 | 4.499° (meets) | **18.488** | 2.671 | 2.363 | 0.4877 |
| **8.5** | **1.30** | **3.933°** | **18.482** | **2.706** | **1.972** | **0.4312** |

**The decisive number: W=9.0 at taper 1.30 gives 18.488 dBe against W=8.5's
18.482 dBe — a difference of 0.006 dBe, i.e. nothing.** The best-at-width CMRR
sequence is now 17.393 (8.0) → 18.482 (8.5) → **18.488 (9.0)**: steps of
+1.09 then **+0.006**. The climb that motivated rounds 9 and 10 has flattened
into a plateau at ≈18.5 dBe.

**Round 10's question is answered: the 20 dBe CMRR spec is NOT reachable by
widening `mm_width`.** The linear extrapolation that put 20 dBe "near W ≈ 9"
was wrong — the deceleration first visible along the taper-1.20 line
(+1.90 → +0.73) continued to zero. This supersedes that projection.

**W=8.5 / taper 1.30 stays the best geometry**, and W=9.0 is worse on every
metric that separates them: phase 4.499° vs 3.933°, imbalance 2.363 vs 1.972
dB, objective 0.4877 vs 0.4312. Only excess loss is a wash (2.671 vs 2.706 dB,
both at the campaign floor). The width axis has genuinely turned over this
time — at *matched* taper, not at a stale taper as in round 7.

**Caveat, stated because it is the one loose end.** W=9.0 was probed at tapers
1.20 and 1.30 only, while the per-width optimum taper has been drifting up
(8.0 → 1.20, 8.5 → 1.30), so 9.0's own optimum may be ≈1.40 and was not
measured. A single point at `W=9.0, y_gap=0.85 (taper 1.40)` would close it.
It is unlikely to change the conclusion — 9.0/1.30 already matches 8.5/1.30 on
CMRR while losing on phase and imbalance — but "the plateau" rests on one
matched-taper comparison until that point exists.

### Campaign summary — where the device stands (2026-09-04)

**Best geometry: `mm_width=8.5`, `y_gap=0.825`, `taper_width=1.30`,
`taper_length=4.0`, `mm_length=159.30` — max|Δφ| = 3.933° (spec ≤5°, MET with
21% margin), min CMRR = 18.482 dBe (spec ≥20 dBe, 1.52 dBe short),
excess loss = 2.706 dB, imbalance = 1.972 dB.**

Progress from the campaign's own starting point (the 28-point 2D optimum at
`mm_width=6.5`, `taper_length=1.0`):

| | start | **now** | change |
|---|---|---|---|
| max\|Δφ\| | 4.764° | 3.933° | −0.83° |
| min CMRR | 12.109 dBe | **18.482 dBe** | **+6.37 dBe** |
| excess loss | 4.076 dB | **2.706 dB** | **−1.37 dB** |
| imbalance | 5.202 dB | **1.972 dB** | **−3.23 dB** |

Every measured quantity improved, and phase stayed inside spec throughout.

**Levers measured against the CMRR spec, and what each is worth:**
- `taper_length` (rounds 8): +1.5 dBe from tl=1→8, logarithmic, ~+0.3/doubling
- `mm_width` (rounds 9-10): +5.1 dBe from 6.5→8.5, then **plateaus at 8.5**
- `y_gap` / `taper_width`: sets whether a width meets phase at all; the
  width-matched-taper rule (6.5→1.065, 7.0→1.10, 7.5→1.20, 8.0→1.20,
  8.5→1.30) is the single most reliable finding of the campaign
- `mm_length`: re-optimized per point, worth ~1-2 dBe within a width

**The remaining 1.5 dBe is not available from any parameter this campaign has
searched.** Closing it needs a change of kind — a different port/aperture
topology, a non-linear taper profile, or a re-read of what CMRR bar this
device actually has to hit.

### ⚠⚠ THE DEVICE IS BANDWIDTH-LIMITED, NOT LOSS- OR CMRR-LIMITED (job 12831929, 2026-09-04)

Five cache-only probes, **01:30 total, zero new simulation**, printing the
per-wavelength breakdown the campaign never recorded. New script
`mmi_4x4_90hybrid_band_probe.py` reuses wg_eval's own solvers and cascade but
keeps the wavelength axis instead of collapsing it, and scores each geometry at
the `best_mm_length` stored in that geometry's own npz. Zero cache misses.

**Best geometry (`W=8.5, y_gap=0.825, taper 1.30, tl=4.0`), per wavelength:**

| λ (µm) | max\|Δφ\| | min CMRR | max XL | max imb |
|---|---|---|---|---|
| 1.290 | 3.933° | 18.734 | 2.706 | 1.619 |
| **1.310** | **3.490°** | **20.761** | **0.654** | **0.940** |
| 1.330 | 3.930° | 18.482 | 1.696 | 1.972 |
| **BAND (reported)** | 3.933° | 18.482 | 2.706 | 1.972 |
| band-edge penalty | 0.44° | **2.28 dBe** | **2.05 dB** | **1.03 dB** |

**At 1.310 µm the device MEETS the 20 dBe CMRR spec (20.761 dBe) and its excess
loss is 0.654 dB.** The 1.52 dBe CMRR shortfall and three-quarters of the
2.7 dB "excess loss" are **band-edge effects, not a lossy or unbalanced
device.** Every conclusion this document has drawn about "the CMRR wall" and
"~3 dB of excess loss" describes a ±20 nm *bandwidth* limit, not the device at
its design wavelength.

W=9.0 is even starker and asymmetric: CMRR reads 21.044 / 21.159 / 18.488 at
1.29 / 1.31 / 1.33 — **it clears 20 dBe at two of three wavelengths and fails
only at the red edge.** An asymmetric band response means the operating point
is not centred, which is a length-tuning question, not a geometry question.

**Answering the question that prompted this — would optimising over
1310 ± 10 nm instead of ± 20 nm help?** Yes, substantially, and it should be
called what it is:
- It would recover a large share of the 2.28 dBe CMRR and 2.05 dB loss penalty,
  very likely putting CMRR over the bar. **That is a specification change, not
  a design improvement** — the same fabricated device, judged over a narrower
  band.
- **It is not yet quantified.** The cached solves are at 1.29/1.31/1.33 only;
  1.30 and 1.32 are new cache keys and need a real ~5.5 h FDTD run at the
  chosen geometry. The ±10 nm worst case lies somewhere between the centre and
  ±20 nm values, and the sampling here (3 points) cannot say where.
- The honest framing for a write-up: report both bands, or state the band with
  the numbers rather than quoting a band-worst figure without it.

**Excess loss IS basis-converged — the 2.7 dB is physical, not truncation.**
The N_MODES=20 confirmation predated loss recording, so this pair closes it at
`W=8.0, y_gap=0.800, tl=4.0`:

| | N_MODES=16 | N_MODES=20 | difference |
|---|---|---|---|
| max XL | 2.956 dB | 2.975 dB | +0.019 (0.6%) |
| max XL at 1.31 | 0.755 dB | 0.755 dB | 0.000 |
| max imb | 2.129 dB | 2.266 dB | +0.137 (6.4%) |
| min CMRR | 17.393 dBe | 17.368 dBe | −0.025 |

Loss agrees to 0.6% and at band centre to three decimals. Imbalance moves 6.4%
at the 1.33 edge (1.4% at centre) — same order, not a truncation artifact.
**So a `waveguide_width` sweep aimed at "fixing 2.7 dB of loss" would be
chasing a band-edge number, and the centre-wavelength loss it would target is
already 0.654 dB.**

**Follow-up launched, also free: job 12832166,
`mmi_4x4_90hybrid_length_rescore.py`.** Every published length in this campaign
was chosen by minimising the live objective, which grades phase linearly and so
keeps buying phase margin *below* the 5° spec at CMRR's expense. The cached
7601-length sweep is re-scored under two further rules — phase-as-constraint,
and maximise band-worst CMRR subject to |Δφ| ≤ 5° — across the five best
geometries. If a different length balances the band, **the 20 dBe spec may be
reachable at ±20 nm with no change to the fabricated geometry at all.**

#### ❌ Length re-selection buys NOTHING — hypothesis refuted (job 12832166, COMPLETED 01:33)

The idea tested: every published length was chosen by minimising the live
objective, which grades phase linearly and so keeps buying phase margin below
the 5° spec at CMRR's expense; a length that balanced the band instead might
lift band-worst CMRR over 20 dBe with no geometry change. **It does not.**

All three selection rules — live objective, phase-as-constraint, and
maximise-band-worst-CMRR-subject-to-|Δφ|≤5° — pick the **identical** length at
every geometry:

| W | y_gap | taper | mm_length (all 3 rules) | max\|Δφ\| | min CMRR | XL | imb | phase-compliant lengths |
|---|---|---|---|---|---|---|---|---|
| 8.5 | 0.825 | 1.30 | 159.30 | 3.933° | 18.482 | 2.706 | 1.972 | 19 / 7601 |
| 9.0 | 0.95 | 1.30 | 179.35 | 4.499° | 18.488 | 2.671 | 2.363 | 15 / 7601 |
| 8.0 | 0.800 | 1.20 | 140.70 | 4.306° | 17.393 | 2.956 | 2.129 | 20 / 7601 |
| 8.5 | 0.925 | 1.20 | 159.40 | 5.438° | 18.126 | 2.949 | 2.413 | **0 / 7601** |
| 9.0 | 1.05 | 1.20 | 179.65 | 5.965° | 16.593 | 2.651 | 2.932 | **0 / 7601** |

**The length axis is exhausted.** The published lengths were already optimal
under every rule, so the campaign's phase-weighted objective never cost it any
CMRR — a nice retro-validation of ~40 points' worth of length selection, and a
clean refutation of the hypothesis that motivated this job. **The 20 dBe spec
is NOT reachable at 1310 ± 20 nm by length tuning with the geometry fixed.**

Two side findings worth keeping:
- **Only 15-20 of 7601 lengths (≈0.25%) meet the phase spec at all** — about
  1 µm of usable length in narrow islands. Comfortable for lithography
  (±0.1 µm is routine), but it means length is a knife-edge parameter that
  cannot be traded for anything else, which is exactly why all three rules
  agree.
- **The two taper-mismatched geometries have ZERO phase-compliant lengths out
  of 7601.** A taper that has not widened with the body cannot be rescued by
  any length whatsoever — the strongest form yet of the width-matched-taper
  rule (6.5→1.065, 7.0→1.10, 7.5→1.20, 8.0→1.20, 8.5→1.30).

**Where that leaves the 1.5 dBe gap.** Two routes remain, and they are
different in kind:
1. **Narrow the judging band.** At 1.310 µm the device already makes
   20.761 dBe. A ±10 nm worst case is bounded between that and the ±20 nm
   figure, but is unmeasured — 1.30/1.32 are new cache keys and cost one
   ~5.5 h job at the best geometry. This is a specification decision.
2. **Change of kind** — different port/aperture topology, or a non-linear
   taper profile (noting the 1×4 tapered-body work in this repo closed both
   profiles after a doubled-resolution, 2×-length search, with a
   good-balance/bad-total-T signature pointing at radiative loss).

A `waveguide_width` sweep is NOT on that list any more: the loss it would
target is a band-edge number (0.654 dB at centre), and loss was just shown
basis-converged.

#### ❌ 1310 ± 10 nm does NOT reach spec — narrowing the band is not worth it (job 12852577, COMPLETED 07:00:22)

Same geometry (`W=8.5, y_gap=0.825, taper 1.30, tl=4.0`), re-optimised over
**1.30 / 1.31 / 1.32 µm**. Nothing was cached at this band, so this was a full
fresh run: 4 FDTD input sections (346.3 min) plus a fresh body solve.

| | ±20 nm (1.29-1.33) | **±10 nm (1.30-1.32)** | 1.310 µm alone |
|---|---|---|---|
| max\|Δφ\| | 3.933° | **3.691°** | 3.490° |
| min CMRR | 18.482 dBe | **18.902 dBe** | 20.761 dBe |
| max excess loss | 2.706 dB | **2.658 dB** | 0.654 dB |
| max imbalance | 1.972 dB | **1.652 dB** | 0.940 dB |
| objective | 0.4312 | **0.3966** | — |
| best mm_length | 159.30 | **157.70** | — |

**Halving the band buys +0.42 dBe of CMRR and 0.05 dB of loss. It does not
reach the 20 dBe spec (18.902), and it does not fix the loss (2.658 dB).**

**Correction to this document's own prediction.** The previous section argued
narrowing "would recover a large share of the 2.28 dBe CMRR penalty, very
likely putting CMRR over the bar." **That was wrong — it recovers 18% of it
(0.42 of 2.28 dBe).** The error was assuming a smooth, roughly quadratic
wavelength dependence, so that halving the band would remove most of the
edge penalty. The measurement says the dependence is far steeper near centre:
by ±10 nm the device has already given up 1.86 dBe of CMRR and 2.00 dB of
excess loss relative to 1.310 µm. **The 20.761 dBe / 0.654 dB figures at
1.310 µm are a narrow peak, not a plateau** — which is consistent with the
length sensitivity found by the re-score (only ~0.25% of lengths meet the
phase spec).

Note this is the *generous* reading of "narrow the spec": mm_length was
re-optimised for the narrow band (159.30 → 157.70 µm), so this is a device
redesigned for ±10 nm, not merely the ±20 nm device judged over a smaller
window. The plain re-labelling would score no better than this.

**Verdict: do not adopt 1310 ± 10 nm.** It halves the operating bandwidth,
still misses CMRR, and leaves excess loss essentially unchanged — the worst
kind of trade. **Report the device at ±20 nm.**

**What this settles about the device.** The CMRR ceiling is not a
bandwidth-accounting artifact after all. Both bands land at 18.5-18.9 dBe, and
the only place 20 dBe appears is a single wavelength. The remaining route to
spec is a change of kind — port/aperture topology or a non-linear taper
profile — and that is now the only route left on the table, since `mm_width`
plateaued (round 10), `taper_length` saturated (round 8), `mm_length` is
exhausted under every selection rule (job 12832166), and the band is exhausted
here.

*Bookkeeping: the `_wl` filename tag added before this run did its job — the
±20 nm baseline `eval_W8.50_y0.82_tl4.00.npz` is intact (wavelengths
[1.29 1.31 1.33], CMRR 18.482, mm_length 159.30) alongside the new
`eval_W8.50_y0.82_tl4.00_wl1.300-1.320.npz` (wavelengths [1.30 1.31 1.32],
CMRR 18.902, mm_length 157.70). Verified by reading both npz files.*

## ✅ Step 6 re-landed at the OPTIMISED geometry — the deliverable gap closed (2026-09-06)

**The gap.** All three published example scripts still carried Step 2's
first-pass point — `mm_width=6.0`, `taper_width=0.8`, `taper_length=1.0`,
`mm_length≈78.94`, with both taper values labelled "PLACEHOLDER — MMI's own
class default, not measured". Everything rounds 8-10 produced lived in
`wg_ygap_opt/*.npz` and in this document, not in the examples the build plan
exists to deliver. Ported now.

**Geometry, all three scripts** (`mmi_4x4_90hybrid_{eme,fdtd,hybrid}.py`):

| | was (first pass) | **now (optimised)** |
|---|---|---|
| `mm_width` | 6.0 | **8.5** |
| `waveguide_offset` | 1.5 | **2.125** |
| `taper_width` | 0.8 (placeholder) | **1.30** (measured, job 12783704) |
| `taper_length` | 1.0 (placeholder) | **4.0** (measured, round 8) |
| `L_PI` | 108.9195 (measured @ 6.0) | **218.5954** (scaled by W², documented as scaled) |
| `mm_length` | 78.9396 | **159.30** |
| measured performance | 4.764° / 12.109 dBe / 4.076 dB / 5.202 dB | **3.933° / 18.482 dBe / 2.706 dB / 1.972 dB** |

**Length pinning changed form.** The old scripts derived `MM_LENGTH` from a
`DELTA_L` on a ±10 µm / 0.25 µm grid around 0.75·L_π, asserted against that
grid (a guard added after a real grid-snap bug). The campaign selected lengths
from the shared **absolute 20-400 µm / 0.05 µm** `MM_LENGTH_GRID`, so the
scripts now pin `MM_LENGTH = 159.30` and assert *that* grid instead, keeping
`DELTA_L` as an informational derived value. The guard matters more here, not
less: only 19 of 7601 lengths meet the phase spec at this geometry.

**Two stale settings fixed at the same time, both measured wrong rather than
merely old:**
1. **Band: 1.26-1.36 → 1.29-1.33 µm** (env-overridable, so the old ±50 nm
   window is still reachable). The geometry was optimised over ±20 nm, and
   this device is steeply wavelength-dependent — min CMRR falls 20.761 →
   18.482 dBe from 1.310 µm alone to ±20 nm. Publishing the optimised
   geometry judged over a 2.5× wider window would report it far outside its
   design band.
2. **`decay_tol`: 1e-11 → 1e-4 with `decay_method="energy"`.** 1e-11 is
   MEASURED UNREACHABLE for this project's solver (the residual floors, then
   diverges — `mmi_1x2_decay_tol_measurement.md`); 1e-4 with energy decay was
   measured sufficient for this device at 0.015 dB and is what all ~45
   campaign points used. Left at 1e-11 these runs could not have converged.

**Geometry smoke-checked before submitting** (no solve): total_length 167.30 µm,
`mm_length` 159.30, ports at ±3.1875 / ±1.0625 µm (spacing 2.125 = the offset),
`input_ys == output_ys` so paper Eq. 8 mirror symmetry holds, and
0.35 < taper_width 1.30 < offset 2.125 so the taper is neither inverse nor
overlapping.

**Runs submitted: 12877997 (`mmi_4x4_90hybrid_eme.py`, 1 node/8 cpu, -t 360)
and 12877998 (`mmi_4x4_90hybrid_hybrid.py`, 1 node/24, -t 900).** The
full-device FDTD (`mmi_4x4_90hybrid_fdtd.py`) is deliberately held: it is the
long pole (the device is now 167 µm long and 11.5 µm wide vs 84 µm before, and
Check 4's runs at the *old* geometry took 14-21 h), and it is only needed for
Step 8's check 4. Submit it once the two above land and confirm the geometry
behaves.

**Steps 8 and 9 remain open**, and Step 8's four checks are all still to run at
this (SOI_220, O-band) cross-section — none of them transfer from another
device or another geometry.

### ✅ Step 6a re-landed at the optimised geometry — EME example (job 12877997, COMPLETED 01:42:38)

Geometry confirmed in the log: `MMI total length: 167.3000 µm
(mm_length=159.3000)` — the port took effect. Ran well inside `-t 360`.

Per-wavelength output at 1.29 / 1.31 / 1.33 µm:

| | 1.29 | 1.31 | 1.33 |
|---|---|---|---|
| imbalance, port 1 / port 3 (dB) | 2.608 / 4.646 | 0.558 / 1.384 | 0.685 / 0.834 |
| excess loss, port 1 / port 3 (dB) | 4.862 / 5.843 | 1.227 / 1.623 | 0.951 / 0.880 |
| min CMRR (dBe) | 8.396 | 24.150 | 21.888 |
| max \|Δφ\| (deg) | 13.366 | 0.836 | 5.962 |

Band-worst: **max|Δφ| = 13.366°, min CMRR = 8.396 dBe, max XL = 5.843 dB,
max imbalance = 4.646 dB** — all set by the 1.29 µm edge.

**These differ substantially from the hybrid-path campaign values for the same
geometry (3.933° / 18.482 dBe / 2.706 dB / 1.972 dB), and that is EXPECTED, not
a sign of a broken port.** This script is **pure EME**, which this campaign
disqualified for this device by direct measurement (the 3-point bias study,
jobs 12352460-78): `EMESolver.solve()` draws one input and one output branch
per call and models no port crosstalk, and at the *old* geometry its phase
error against full-device FDTD ran +43.9° to +62.0°. A large EME-vs-hybrid gap
here is that same documented limitation, re-appearing at the new geometry.

Two things the EME run does corroborate independently:
- **The steep wavelength dependence.** At 1.31 µm alone this script reports
  0.836° and 24.150 dBe; at 1.29 µm, 13.366° and 8.396 dBe. Same
  narrow-peak-not-plateau shape the hybrid path showed (job 12831929), from a
  different solver.
- **The direction of the EME bias is optimistic at centre** (24.15 dBe vs the
  hybrid path's 20.76 dBe at 1.310), consistent with neglected crosstalk
  flattering the device.

The like-for-like comparison against the campaign numbers is
`mmi_4x4_90hybrid_hybrid.py` (job 12877998), still running.

### ❌ Hybrid example FAILED at 05:51:25 — my port bug, fixed, NOT auto-resubmitted (job 12877998)

```
NameError: name '_DELTA_L_GRID' is not defined
  File ".../mmi_4x4_90hybrid_hybrid.py", line 290, in <module>
    _body_lengths_grid = (0.75 * L_PI + _DELTA_L_GRID) - 2 * body_stub
```

**Cause: the geometry port removed `_DELTA_L_GRID`'s definition but not its
use.** The length block was rewritten to pin `MM_LENGTH = 159.30` on the shared
absolute grid, which deleted the `_DELTA_L_GRID = np.linspace(-10, 10, 81)`
definition — but 140 lines further down, the body length sweep still built its
`lengths` array from that grid. `mmi_4x4_90hybrid_{eme,fdtd}.py` do not do a
length sweep and were unaffected (the EME run completed fine).

**The failure came AFTER the expensive part, not before it.** The traceback is
at 5:51:25, i.e. the FDTD interface solves had already completed and are now in
`example_sim_cache/`. A resubmission should skip straight past them.

**Fix — and it makes the script cheaper than it was.** The old grid was Step 3's
±10 µm / 0.25 µm `DELTA_L` window, whose cache entries all live at the
superseded `mm_width=6.0` geometry and would have missed anyway. It now uses
the **shared absolute `MM_LENGTH_GRID` (20-400 µm @ 0.05 µm)** — the same grid
every campaign point and the `mm_width=8.5` Phase A body solve (job
12780661_0) used, reproduced inline rather than imported so a published example
does not reach into `examples/verification/`. `body_stub = 0.25` and the
boundary construction are byte-identical to `wg_eval.boundaries_for()`, which
is what makes that a **cache hit** instead of a fresh ~70 min body solve.
Guarded with an assert that `MM_LENGTH` is actually on the grid (it is: index
2786 of 7601 → 159.30).

**Not resubmitted.** Standing directive of 2026-08-17: if these long runs
crash, report and wait for direction rather than auto-resubmitting.

### ✅✅ Step 6c re-landed — the example REPRODUCES the campaign numbers exactly (job 12889625, COMPLETED 01:11:34)

Resubmitted after the `_DELTA_L_GRID` fix. **01:11:34 vs the 05:51:25 the failed
attempt burned** — the FDTD interface solves came from cache, and the body sweep
hit the `mm_width=8.5` Phase A entry exactly as the fix intended:
`body solved (shared MM_LENGTH_GRID reused, sliced at mm_length=159.30 um):
(1, 3, 16, 16)`. Geometry confirmed: `MMI total length: 167.3000 µm
(mm_length=159.3000)`, ports at ±3.1875 / ±1.0625 µm, `input_ys == output_ys`.

Per wavelength (1.29 / 1.31 / 1.33 µm):

| | 1.29 | 1.31 | 1.33 |
|---|---|---|---|
| imbalance p1 / p3 (dB) | 0.915 / 1.619 | 0.725 / 0.940 | 1.039 / 1.972 |
| excess loss p1 / p3 (dB) | 2.282 / 2.706 | 0.395 / 0.654 | 1.696 / 1.036 |
| min CMRR (dBe) | 18.734 | 20.761 | 18.482 |
| Δφ_6 / Δφ_7 / Δφ_8 (deg) | −2.642 / 1.682 / 3.933 | −3.490 / −0.153 / 3.391 | 2.467 / 3.930 / 2.192 |

**Band-worst, against the campaign's hybrid-path values for this geometry:**

| | campaign (wg_eval) | example (this run) | difference |
|---|---|---|---|
| max\|Δφ\| | 3.933° | **3.933°** | 0.000 |
| min CMRR | 18.482 dBe | **18.482 dBe** | 0.000 |
| max excess loss | 2.706 dB | **2.706 dB** | 0.000 |
| max imbalance | 1.972 dB | **1.972 dB** | 0.000 |

**Exact to three decimals on all four metrics** — and the per-wavelength values
match the band probe (job 12831929) exactly too (1.31 µm: 0.654 dB / 0.940 dB /
20.761 dBe). This is a reproduction, not an approximation, which is the
expected outcome: the published example's hybrid path and `wg_eval` build the
same cascade at the same `N_MODES=16`, off the same cached body, at the same
length. **The port is verified end to end** — the deliverable example now
reports the optimised device, and reports it identically to the optimiser that
found it.

Artifacts written to
`example_output/mmi_4x4_90hybrid/mmi_4x4_90hybrid_hybrid/`:
`smatrix_hybrid_port{1,3}.pkl`, `hybrid_results.npz`, `hybrid_results.png`.

**Step 6 status: 6a (EME) ✅ and 6c (hybrid) ✅ both re-landed at the optimised
geometry. 6b (full-device FDTD) is the remaining one**, and it is also Step 8's
check 4. It is the long pole: the device is now 167.3 µm long and 11.5 µm wide
against 84 µm before, and check 4's runs at the *old* geometry took 14-21 h.

### Dense wavelength sweep — EME at N_WL=11 submitted (job 12904783, 2026-09-07)

**Why.** Every plot in this campaign interpolates between exactly three
wavelengths (1.29 / 1.31 / 1.33), and the device's dominant behaviour *is* its
wavelength dependence — the sharp 54.8 dBe CMRR "peak" and the V-shaped loss
minimum in `hybrid_results.png` are shapes imposed by 3-point sampling, not
measured curves. 11 points over the same design band gives 4 mm resolution in
wavelength instead of 20.

**The cost asymmetry, read from the solvers rather than assumed:**
- **FDTD is broadband** — `fdtd.py:170` builds DFT monitors from
  `monitor_freq_grid(wavelengths)` and drives one wideband source
  (`source_fwidth`). One time-domain run yields every wavelength, so extra
  wavelengths are nearly free in stepping time.
- **EME is per-wavelength** — `eme.py` loops
  `for wl_idx, wl in enumerate(wavelengths)` with a mode solve each, at six
  separate sites. Cost is linear in `N_WL`.
- The hybrid script is a mix: flat FDTD stepping plus per-wavelength interface
  mode bases and a per-wavelength EME body solve (those "Dominant planewave for
  band N" lines that filled its log).
- **No partial cache reuse.** `cache/config.py` hashes `wavelengths` as an
  array, so changing `N_WL` invalidates the 3-point entries wholesale — the
  same thing that made the ±10 nm run re-solve from scratch.

Measured 1:42:38 at N_WL=3 → **~6.3 h expected at N_WL=11**, so
`run_mmi_4x4_90hybrid_eme_nwl11.batch` raises `-t` from 360 to 480.

**Fourth instance of the overwrite hazard, guarded before launch.** These
examples write fixed filenames (`eme_results.npz`, `hybrid_results.png`,
`smatrix_*_port{1,3}.pkl`) with no `N_WL` in the path, so a denser run would
have silently destroyed the 3-point results verified on 2026-09-06 to reproduce
the optimiser exactly. Both scripts now tag `output_dir` with `_nwl{N}` whenever
`N_WL != 3` (default stays untagged, so every existing path keeps working), and
`mmi_4x4_90hybrid_hybrid.py`'s hardcoded `N_WL = 3` became env-overridable via
`MMI4X4_HYB_NWL` to match the EME script's `MMI4X4_EME_NWL`. Verified: N_WL=3 →
untagged dir; 7 → `_nwl7`, grid 1.2967/1.3033/… ; 11 → `_nwl11`, grid at 4 nm
spacing.

**EME first deliberately** — besides being the cheap one, it measures the
per-wavelength mode-solve cost at this cross-section directly, which is the
unknown term in the hybrid estimate (~8-21 h at N_WL=7, the spread being how
much of its cold 5:51 was FDTD stepping vs mode solves). The hybrid dense run
follows once that number is real.

### ⚠⚠ DENSE SWEEP OVERTURNS THE "NARROW PEAK" READING — there is a ~24 nm dual-spec window (job 12904783, COMPLETED 06:07:19)

`mmi_4x4_90hybrid_eme.py` at **N_WL=11** (4 nm spacing over 1.29-1.33), same
optimised geometry. Output went to the tagged directory
`mmi_4x4_90hybrid_eme_nwl11/` and **the 3-point results are intact** —
`mmi_4x4_90hybrid_eme/eme_results.npz` still dated Sep 6 11:11, the new one
Sep 7 15:33. The `_nwl` tag did its job.

| λ (µm) | min CMRR | max XL | max imb | max\|Δφ\| |
|---|---|---|---|---|
| 1.290 | 8.396 | 5.843 | 4.646 | 13.366 |
| 1.294 | 12.035 | 4.891 | 2.849 | 9.552 |
| 1.298 | 15.888 | 3.860 | 1.910 | 6.195 |
| 1.302 | 19.659 | 2.929 | 1.509 | 3.387 |
| 1.306 | 22.735 | 2.177 | 1.389 | 1.582 |
| **1.310** | **24.150** | 1.623 | 1.384 | **0.836** |
| 1.314 | 23.824 | 1.248 | 1.382 | 1.014 |
| 1.318 | 22.272 | 1.022 | 1.287 | 1.912 |
| 1.322 | 20.905 | 0.905 | 1.020 | 3.424 |
| 1.326 | 20.551 | **0.865** | **0.818** | 4.766 |
| 1.330 | 21.888 | 0.951 | 0.834 | 5.962 |

**Both "shapes" the 3-point sampling implied are wrong.**

1. **CMRR is not a triangular peak — it is a broad plateau with a steep blue
   edge.** It climbs 8.4 → 24.2 dBe from 1.290 to 1.310, then decays *gently*:
   23.8, 22.3, 20.9, 20.6, and **turns back up to 21.9 at 1.330**. CMRR stays
   **above 20 dBe from ≈1.304 all the way to 1.330** — a ≥26 nm window. The
   3-point samples (8.4 / 24.2 / 21.9) linearly interpolated into a spike and
   hid the entire plateau.
2. **Excess loss is not a V centred at 1.310.** It decays monotonically across
   almost the whole band, 5.843 → 0.865 dB, with its **minimum at 1.326, not
   1.310**, and only a token rise to 0.951 at 1.330. Imbalance behaves the same
   way (minimum 0.818 dB at 1.326).

**Consequence — the device has a real operating window, not a single good
wavelength.** Taking both specs together: CMRR crosses 20 dBe between 1.302 and
1.306, and |Δφ| crosses 5° between 1.326 and 1.330, giving a dual-compliant
band of roughly **1.304-1.328 µm (~24 nm)** — and it is **centred near 1.316,
not 1.310.**

**This corrects two statements written earlier in this document**, both of
which were 3-point artifacts: "the 20.761 dBe / 0.654 dB figures at 1.310 µm
are a narrow peak, not a plateau", and "the only place 20 dBe appears is a
single wavelength". Neither survives 11-point sampling.

**Caveat, and it is the important one: this is the EME path.** Pure EME is
disqualified for *absolute* numbers on this device (no port crosstalk; the
3-point comparison shows it optimistic by ≈+3.4 dBe on CMRR at both 1.310 and
1.330 against the hybrid path). What transfers is the **shape**, and the shape
says plateau. Applying the ≈3.4 dBe offset to the plateau puts the hybrid-path
CMRR somewhere near 17-21 dBe across 1.31-1.33 — which is exactly the range
where "does it meet spec over a usable band" is decided, and cannot be settled
without the hybrid dense sweep.

**Cost, now measured rather than extrapolated.** 06:07:19 for 11 wavelengths
vs 01:42:38 for 3 → **2003 s/wavelength vs 2053 s/wavelength**. EME is cleanly
linear in `N_WL`, within 2.5%, so the earlier 6.3 h estimate was right.

**Hybrid dense-sweep estimate.** Its cold 3-wavelength run was 5:51:25, of
which the post-interface part (body solve + 7601-length sweep + metrics)
measured 4274 s when interfaces were cached. FDTD *stepping* is flat in `N_WL`
but the interface mode bases are not, so the honest bracket is:
- fully linear → **N_WL=7 ≈ 13.7 h, N_WL=11 ≈ 21.5 h**
- stepping-dominated → **N_WL=7 ≈ 9 h, N_WL=11 ≈ 13 h**

`-t` must go to 1440 either way; `-t 900` would risk a timeout at N_WL=11.

### ⚠ HYBRID DENSE SWEEP: THE EME PLATEAU DOES NOT SURVIVE — narrower window, shifted lower, peaked not flat (job 12963168, COMPLETED 2026-09-09 01:01:56, 13:44:28 elapsed)

`mmi_4x4_90hybrid_hybrid.py` at **N_WL=7** (1.29, 1.2967, 1.3033, 1.31, 1.3167,
1.3233, 1.33), same optimised geometry (mm_width=8.5, taper_width=1.30,
taper_length=4.0, mm_length=159.30). Output went to the tagged directory
`mmi_4x4_90hybrid_hybrid_nwl7/` and **the 3-point directory is intact** —
`mmi_4x4_90hybrid_hybrid/hybrid_results.npz` still dated Sep 6 19:39, the new
one Sep 9 01:01. All three previously-verified 3-point values reproduce
exactly at their matching indices (18.734/2.706/3.933 at 1.29; 20.761/0.654/
3.490 at 1.31; 18.482/1.696/3.930 at 1.33) — confirms this is the same
evaluator, not a re-derivation.

| λ (µm) | min CMRR (dBe) | max XL (dB) | max imb (dB) | max\|Δφ\| (°) | CMRR≥20 | \|Δφ\|≤5 | both |
|---|---|---|---|---|---|---|---|
| 1.2900 | 18.734 | 2.706 | 1.619 | 3.933 | ✗ | ✓ | ✗ |
| 1.2967 | 23.085 | 1.573 | 1.408 | 1.367 | ✓ | ✓ | ✓ |
| 1.3033 | 22.745 | 0.930 | 1.334 | 1.488 | ✓ | ✓ | ✓ |
| **1.3100** | **20.761** | 0.654 | 0.940 | 3.490 | ✓ | ✓ | ✓ |
| 1.3167 | 21.577 | 0.640 | 0.726 | **5.619** | ✓ | ✗ | ✗ |
| 1.3233 | 19.616 | 1.135 | 1.341 | 6.041 | ✗ | ✗ | ✗ |
| 1.3300 | 18.482 | 1.696 | 1.972 | 3.930 | ✗ | ✓ | ✗ |

(`min CMRR` = min over {CMRR_I, CMRR_Q} × {port 1, port 3}; `max XL`/`max imb`
= max over the two driven ports, matching the evaluator's existing metric
definitions.)

**Answer to the open question: the hybrid path does *not* show the EME
sweep's broad flat plateau.** The EME 11-point sweep found CMRR climbing to a
plateau ≥20 dBe from ≈1.304 to 1.330, decaying gently and turning back up at
the top edge. The hybrid path instead shows a **single peaked/decaying curve**
— it jumps from 18.734 (1.290, below spec) up to a peak of 23.085 at 1.2967,
then declines with a small wobble (22.745, 20.761, 21.577) before dropping
further to 19.616 and 18.482 at the two longest wavelengths. There is no
turn-back-up at 1.330 in the hybrid path — it keeps falling, mirroring the
already-known 3-point value there rather than the EME plateau's recovery.

**CMRR-only window** (ignoring phase): 1.2967-1.3167 µm, ≈20 nm, contiguous.
**Dual-compliant window** (CMRR≥20 dBe **and** max\|Δφ\|≤5°): **1.2967-1.3100
µm, ≈13.3 nm**, centred near **1.303 µm** — narrower than the EME-based
estimate of ~24 nm, and shifted to a *lower* wavelength than EME's ~1.316 µm
centre. The phase constraint, not CMRR, is what cuts the window short at the
long-wavelength end: 1.3167 clears CMRR (21.577) but fails phase (5.619° >
5°), a failure mode the EME-only view (no phase-degrading crosstalk paths
modelled) could not have shown.

**This is consistent with, not contradicted by, the previously-documented
+3.4 dBe EME optimism** — the offset was measured only at 1.310 and 1.330 in
the 3-point comparison, both of which are exactly where the two paths diverge
most in *shape*, not just magnitude. The hybrid path's real port crosstalk
evidently reshapes the curve rather than uniformly translating it.

**Net finding for the device:** a genuine dual-compliant operating window
exists, but it is materially narrower (~13 nm vs the ~24 nm EME suggested) and
sits at the *blue* edge of the design band (1.2967-1.3100 µm) rather than
centred on it. 1.310 µm — the original single-wavelength design target — is
the *edge* of this window, not its centre; 1.290 and 1.330 µm, and everything
at 1.3167 µm and above, fail at least one spec.

**Cost.** 13:44:28 (49445 s) for N_WL=7, inside the predicted 9-13.7 h bracket
(near its upper end) — confirms the interface mode bases are not purely
stepping-dominated; the linear extrapolation is closer to reality than the
optimistic one. Per-wavelength marginal cost above the 3-point baseline:
(49445 s - 21085 s [5:51:25 cold 3-pt run]) / 4 new wavelengths ≈ 7090
s/wavelength — consistent with N_WL=11 costing roughly 21-22 h if run next.

### Step 8, Check 1 (engine agreement) — short-body FDTD-vs-hybrid probe at THIS device's own cross-section (jobs 13024684 / 13024689, both COMPLETED 2026-09-09)

Full-device FDTD at the production optimum (mm_length=159.30) was excluded as
too big (~100-140 h/port). Shrinking `mm_length` to <20 µm shrinks the
single-domain FDTD box ~6-12x by volume and makes a real Check 1 (engine
agreement, this device's own SOI_220/O-band/mm_width=8.5 cross-section)
affordable. Both scripts gained an `MMI4X4_MM_LENGTH` env override for this
(default 159.30 unchanged; output dirs get a `_L{length}` suffix whenever
used, so nothing collides with the production 159.30 µm results).

**mm_length=5.00 µm is NOT a device operating point** — L_pi=218.5954, so
3/4·L_pi≈164 µm is where self-imaging actually happens; at 5 µm the device
is nowhere near it. Both methods agree on this at face value: CMRR is
essentially 0 dBe, imbalance 33-74 dB, phase errors up to 171°, all far
outside spec. This run exists purely to compare the two methods' raw
S-matrices at a length cheap enough for full-device FDTD to reach true decay.

**Cost, measured, not extrapolated:**
- FDTD probe (job 13024684): single node, 24 ranks (deliberately NOT the
  96-rank/4-node config the full-length device measured as efficient — that
  was tuned for a 159M-voxel domain, 3.3M/rank at 48; this ~20.6M-voxel
  domain would put only ~215k voxels/rank at 96, likely communication-bound).
  Steady **0.134 s/step**, comparable to the OLD full-length domain's
  96-rank rate (0.147 s/step) despite ~7.7x fewer voxels — confirms this
  small box is already communication-bound at only 24 ranks, not
  compute-bound; more ranks would likely not have helped here. Both ports
  (1 and 3) fully decayed to `decay_tol=1e-4` well inside budget: **1:49:41
  total wall**, vs the 6h budget and the full-length device's ~100-140h
  estimate — a ~55-75x speedup, roughly matching the ~7x volume shrink times
  a shorter physical ring-down time.
- Hybrid probe (job 13024689): single node, 24 ranks, **6:58:47 total wall**
  — NOT proportionally cheaper than the full-length production run's own
  ~9h fresh-solve precedent, confirming the prediction in both scripts'
  updated docstrings: the 4 FDTD input-section solves span a region whose
  physical size (taper_length + body_stub) is independent of `mm_length`,
  and the body EME mode basis is a per-cross-section, not a per-length,
  solve. Shortening the body bought no discount here, as expected.

**Results at mm_length=5.00 µm, both methods** (Δφ_6/7/8, wavelengths
1.29/1.31/1.33 µm — all values expected to be far outside spec, this is not
a device-performance measurement):

| | full-device FDTD | hybrid |
|---|---|---|
| Input 1 imbalance (dB) | 49.517 / 38.116 / 42.526 | 73.529 / 69.179 / 71.547 |
| Input 1 excess loss (dB) | 3.444 / 3.302 / 3.470 | 3.378 / 3.376 / 3.362 |
| Input 1 CMRR_I (dBe) | 0.00117 / 0.00314 / 0.00097 | ~0 (7.7e-7 / 2.1e-6 / 1.2e-6) |
| Input 1 CMRR_Q (dBe) | 0.0165 / 0.174 / 0.625 | 0.0249 / 0.0364 / 0.0460 |
| Input 3 imbalance (dB) | 47.689 / 37.857 / 33.486 | 48.407 / 46.120 / 44.521 |
| Input 3 excess loss (dB) | 3.705 / 3.192 / 3.641 | 3.367 / 3.423 / 3.477 |
| Δφ_6 (deg) | −120.2 / −117.6 / −21.9 | −70.0 / −78.9 / −88.5 |
| Δφ_7 (deg) | −154.4 / −143.0 / +61.4 | −47.9 / −65.5 / −84.7 |
| Δφ_8 (deg) | −46.0 / −0.9 / +171.4 | +84.8 / +119.1 / +98.2 |

The derived spec metrics themselves disagree substantially between methods
at this length — expected: at an arbitrary off-self-imaging length, small
differences in how each method distributes power among nominally-equal
output ports get amplified by the CMRR/phase-error ratios, which are not a
meaningful comparison here. **The comparison that matters is the raw
per-port |T|² transmission**, computed directly from both `SMatrix` objects:

**Input port 1 driven** (out5 = the dominant, "correct" output at this port):

| port | hybrid 1.29/1.31/1.33 | FDTD 1.29/1.31/1.33 | max abs | max rel |
|---|---|---|---|---|
| out5 | 0.4548/0.4544/0.4552 | 0.4472/0.4603/0.4444 | 0.0108 | 2.4% |
| out6 | 0.00454/0.00518/0.00583 | 0.00527/0.00707/0.00514 | 0.0019 | 26.8% |
| out7 | 7.0e-6/1.1e-5/1.5e-5 | 5.0e-6/7.1e-5/1.85e-4 | 0.0002 | 91.6% |
| out8 | 0/0/0 | 3.0e-5/8.3e-5/2.5e-5 | 0.0001 | 99.9% |

Totals: 1.29 µm hybrid=0.4594 FDTD=0.4525 (**+1.5%**); 1.31 µm
hybrid=0.4596 FDTD=0.4675 (**−1.7%**); 1.33 µm hybrid=0.4611 FDTD=0.4497
(**+2.5%**).

**Input port 3 driven** (out7 = the dominant output at this port):

| port | hybrid 1.29/1.31/1.33 | FDTD 1.29/1.31/1.33 | max abs | max rel |
|---|---|---|---|---|
| out5 | 7.0e-6/1.1e-5/1.5e-5 | 7.0e-6/7.6e-5/1.88e-4 | 0.0002 | 91.8% |
| out6 | 0.00467/0.00530/0.00595 | 0.00564/0.00639/0.00747 | 0.0015 | 20.3% |
| out7 | 0.4514/0.4442/0.4373 | 0.4152/0.4660/0.4198 | 0.0362 | 8.7% |
| out8 | 0.00454/0.00518/0.00583 | 0.00522/0.00702/0.00506 | 0.0018 | 26.2% |

Totals: 1.29 µm hybrid=0.4606 FDTD=0.4261 (**+8.1%**); 1.31 µm
hybrid=0.4547 FDTD=0.4795 (**−5.2%**); 1.33 µm hybrid=0.4491 FDTD=0.4325
(**+3.8%**).

**This does NOT reproduce the old geometry's engine-agreement pattern.** The
old-geometry comparison (job 11778128/12272573) found a clean, monotonic
trend: FDTD found MORE power than hybrid at every wavelength, by an amount
that grew steadily toward short wavelength (+30.2%/+49.6% at 1.26 µm down to
−0.3% at 1.36 µm — one consistent sign, one consistent direction). Here the
sign **oscillates** with wavelength for both driven ports (port 1:
+1.5%/−1.7%/+2.5%; port 3: +8.1%/−5.2%/+3.8%) — no monotonic trend toward
either band edge, and no single dominant sign. The dominant "correct" output
port agrees far better for port 1 (out5, max 2.4%) than for port 3 (out7,
max 8.7%), an asymmetry the old geometry's I/O structure didn't show either.

**Interpretation.** The old geometry's finding was measured at a length near
its own self-imaging condition, where the beat pattern between the excited
modes is smooth and slowly varying with wavelength — exactly the regime
where a systematic, monotonic bias (e.g. from the hybrid's port-isolation
crosstalk approximation) would show up cleanly. At mm_length=5 µm the device
is far from any self-imaging condition; many modes are simultaneously
mid-beat, and small differences in mode content and phase between the two
methods land on different parts of an oscillatory pattern at each
wavelength, which can flip sign non-monotonically without implying a
different *magnitude* of underlying error. **The magnitude of disagreement
here (max port-total 8.1%, single dominant-port max 8.7%) is comparable to,
not larger than, the old geometry's 1.31/1.36 µm points (+4.1%/−0.3%,
+3.9%/−0.3%)** — this short-length probe does not show the two engines
agreeing WORSE at this device's own cross-section than they did at the old
one; it shows a *different-shaped* disagreement, consistent with the
length-regime argument above rather than a cross-section-specific problem.

**What this settles and what it doesn't.** It gives real, this-cross-section
engine-agreement data (Check 1) that did not exist before, at a cost 55-75x
cheaper than the full production length would require. It does NOT settle
whether the OLD geometry's specific "hybrid under-collects, worst at short
wavelength" mechanism (candidate: the port-isolation crosstalk
approximation, per the original Check 4 section above) is present or absent
at mm_width=8.5's wider/closer-packed ports — the short length here changes
too many things at once (length regime, not just cross-section) to isolate
that. A clean test of the port-isolation hypothesis specifically would need
either a length near this device's OWN self-imaging point (159.30 µm,
prohibitively expensive) or a different diagnostic that varies port spacing
independent of length.

### Step 8, Check 2 (mode identity) — audit FAILED on a duplicate-mode guard, diagnostic dump does NOT reproduce it (jobs 13060190 FAILED, 13062141 COMPLETED, 2026-09-10)

Added a `_mmi_4x4_90hybrid()` entry to `examples/verification/python/
audit_interface_modes.py`'s `DEVICES` (this device had no Check 2 entry
before now). Config: single-cut (not `MULTI_DEVICES`), section=(-1, 1) at
body center — the body is one uniform mm_width=8.5 slab with straight
sidewalls over its whole length, so a short interior cut audits the same
cross-section `HybridSolver` actually contracts at (near-edge, body_stub
inside the taper seam) by the same translation-invariance argument every
other MMI device's own Check 2 config already relies on. n_modes=16 (this
device's real `N_MODES`), wl=1.31, stack=SOI_220, sim_width=11.5,
sim_height=2.0, res=40 — pinned to `mmi_4x4_90hybrid_hybrid.py` exactly.

**Job 13060190 FAILED at 41:33 elapsed** (`AUDIT_ARGS=mmi_4x4_90hybrid`, no
`--convergence`, single node/1 task/8 cpus, serial):

```
RuntimeError: FDTDSolver: could not find 16 distinct eigenmodes within 24
band_num attempts at this cross-section -- meep/MPB returned a duplicate
more times than this guard tolerates. See validation_markdown/
modal_multiplexer_validation.md gap 6.
```

This is `FDTDSolver._iter_distinct_bands()`'s own duplicate-eigenmode guard
(`beef/solvers/fdtd.py`) — it silently skips any `band_num` whose field
overlaps an already-accepted one above 0.9 (using its own small, fixed
16×8-point diagnostic grid, independent of the caller's real resolution)
and raises once it exhausts `n_modes + 8 = 24` attempts. This device's own
Step 8 plan flagged exactly this risk in advance (A7): four closely-spaced
ports per end puts it closer to the modal multiplexer's 16-mode regime
(where this guard was originally built after a real duplicate-mode bug)
than to either existing MMI's simpler port count.

**Diagnostic follow-up: `mmi_4x4_90hybrid_mode_field_dump.py` (NEW script),
job 13062141 COMPLETED 1:15:05.** The guard discards the fields it skips,
so nothing in job 13060190's log said *which* bands collided or what they
looked like. This new script re-walks the identical cross-section, same
k-seed-laddering algorithm, band-by-band via `sim.get_eigenmode()` directly
— but keeps **every** raw band's field (no skip/raise), plots each one's
`|Ey|²+|Ez|²` intensity in the y-z plane, and computes the **full 30×30
self-overlap matrix** on a much finer 231×41-point (0.05 µm) grid, rather
than inferring anything from where an exception fired.

**Result: all 30 requested bands came back guided** (this 8.5 µm slab
supports far more than 16 guided modes — n_eff ranged 2.963 down to 2.287,
cladding is 1.4468, nowhere close to the guided edge even at band 30), and
**the self-overlap table found NO duplicates at all** — max off-diagonal
overlap across every one of the 30×29/2 pairs is **0.017** (bands 22/27),
orders of magnitude below the 0.9 threshold that trips the guard. Full
table and 30 field plots at `example_output/verification/
mmi_4x4_90hybrid_mode_field_dump/`.

**This does not reproduce the audit's failure.** Bands 1-24 (the exact range
the guard exhausted its attempts within) show clean, well-separated n_eff
(spacing 0.02-0.05/band) and no field-shape overlap anywhere near duplicate
territory. Bands 25-30 show n_eff spacing collapsing sharply (2.331 → 2.328
→ 2.321 → 2.320 → 2.313 → 2.301 → 2.287, vs ~0.02-0.05/band through 24) —
modes genuinely crowding as expected approaching a dense higher-order
region on a wide multimode slab — but even there, the finest-spaced pair
(bands 24/26, Δn_eff=0.0018) overlaps at only 0.0003 in field shape. Not
one pair anywhere in the walk is ambiguous, let alone a duplicate.

**Leading hypothesis, not yet confirmed: the guard's own coarse diagnostic
grid is the likely source of a false positive, not a real degeneracy.**
`_iter_distinct_bands()`'s duplicate check uses a FIXED `_DUP_BAND_DIAG_NY
=16, _DUP_BAND_DIAG_NZ=8` grid regardless of the caller's real resolution
or geometry — 8 points spanning the *whole* `sim_height=2.0` µm gives
~0.25 µm z-spacing, against a **0.22 µm-thick guided core**. That is under
one grid cell resolving the entire vertical mode profile of every mode in
the basis, for a device whose sole distinguishing structure in z is that
0.22 µm slab. The 16×8 diagnostic grid was sized generically (its own
docstring: "it only needs to tell 'same eigenvector' from 'different
eigenvector'... stays cheap regardless of how finely a caller wants its own
returned fields sampled") and was validated against the modal multiplexer's
own cross-section (also SOI_220-ish core height, but not audited against a
wide mm_width=8.5 slab with this many closely-packed higher-order modes
sharing similar vertical structure) — a resolution adequate for a handful
of well-separated low-order modes can plausibly alias two genuinely
different but similarly-shaped high-order lobes into an apparent overlap
above 0.9 once enough modes are requested on a wide enough slab.

**Not yet established as fact** — this is the strongest available
explanation given the data (clean fine-grid result vs. a real coarse-grid
failure on the identical cross-section/algorithm), but confirming it
requires directly comparing the coarse-grid overlap against the fine-grid
one on the SAME live eigenmode objects (the field-dump job did not save raw
field arrays, only the fine-grid overlap matrix, so this needs one more
short run rather than post-hoc analysis of what's already on disk).

**Status: Check 2 is OPEN, not closed and not failed-as-a-device-problem.**
The evidence so far argues the device's real 16-mode interface basis is
very likely sound (no genuine near-degeneracy found at any resolution
tested), but the production audit tool itself needs either a resolution fix
(scale `_DUP_BAND_DIAG_NY`/`_DUP_BAND_DIAG_NZ` to the caller's real z-extent
relative to core thickness) or a documented exception for wide slabs before
this device's Check 2 can be run to actual completion through the shared
`audit_interface_modes.py` harness.

### ✅✅ Check 2 coarse-grid-aliasing hypothesis CONFIRMED — real shared-library bug, not a device problem (job 13064925, COMPLETED 59:16)

`mmi_4x4_90hybrid_mode_dup_grid_compare.py` (NEW script) computed BOTH the
guard's exact 16×8 diagnostic grid and the 231×41 fine grid from the SAME
live eigenmode objects, one pass, 27 bands. Full table at
`example_output/verification/mmi_4x4_90hybrid_mode_dup_grid_compare/
grid_compare_table.md`.

**Confirmed, and far more extensively than expected.** 22 of 27 bands flag
a coarse-grid-only "duplicate" (coarse overlap 0.96–0.996, fine overlap
0.0006–0.017 for the identical pair). This is not scattered noise — it is a
clean **spatial-aliasing pattern**: band *i* coarse-duplicates band
*(24 − i)* almost exactly across the whole range (1↔23, 2↔22, 3↔21, 4↔20,
5↔19, 6↔18, 7↔17, 8↔16, 9↔15, 10↔14, 11↔13 — the last pair peaking at
0.996, the tightest "duplicate" in the table). Band 12, sitting at this
series' center, is the one clean exception (coarse overlap only 0.140).
Past the mode-crowding point found in the earlier dump (bands 24-27), a new
false pair appears (27↔22, coarse 0.909) with a different mechanism (genuine
n_eff crowding there, not the same mirror pattern).

**Mechanism, now understood, not just hypothesized:** `_DUP_BAND_DIAG_NY=16`
samples across the full `sim_width` gives roughly one sample per
`mode_index/8` lobes for this cross-section's mode family. Once a mode's
transverse lobe count exceeds what 16 points can resolve (Nyquist), the
coarse sampling **aliases** higher-order spatial oscillation back onto a
lower apparent pattern — and modes with complementary/mirrored lobe counts
alias onto nearly-identical sampled values, which is exactly the *(i, 24−i)*
mirror pattern observed. This is a real, generic failure mode: the guard's
diagnostic grid is a small FIXED constant (16×8, chosen once, validated
against the coupler/mmi_1x2/mmi_2x2/modal_multiplexer's narrower/lower-mode
cross-sections) with no relationship to how many modes or how wide a slab a
given caller requests — it will alias on ANY device whose requested mode
count approaches or exceeds what 16 transverse samples can resolve, not
just this one.

**This closes the device-side half of Check 2.** The device's real 16-mode
interface basis is sound — no genuine near-degeneracy exists among the
first 27 guided modes at any resolution actually capable of resolving mode
shape (confirmed twice now: the fine-grid dump and this same-run
comparison). **The audit tool itself has a real bug**: `FDTDSolver.
_iter_distinct_bands()`'s duplicate guard produces false positives on wide,
many-mode cross-sections, and would have silently (or, as here, loudly)
mis-served this device's Check 2 regardless of how many times it was
re-run, because the failure is deterministic given the same k-seed-laddering
algorithm and the same undersized diagnostic grid.

**Fix needed, not yet made:** `_DUP_BAND_DIAG_NY`/`_DUP_BAND_DIAG_NZ` should
scale with the cross-section actually being queried (e.g. tied to
`n_modes` and `sim_width`/core thickness, Nyquist-safe for the highest
mode index requested) rather than being a fixed 16×8 constant. This is a
`beef/solvers/fdtd.py` change, shared across every device that calls
`section_mode_neff()`/`section_mode_fields()`/`audit_interface_modes.py`,
not a `mmi_4x4_90hybrid`-local one — filed here as an open item, not fixed
inline, since it touches the shared library's own duplicate-guard contract
and the four already-audited devices' recorded results.

**Step 8 Check 2 status: DEVICE SIDE CLOSED, TOOL SIDE OPEN.** This
device's interface mode identity/ordering is not blocked by anything found
here; the `audit_interface_modes.py` harness cannot currently complete a
clean run for this device without either the diagnostic-grid fix above or a
one-off override (e.g. passing a wider diagnostic grid just for this
device's config) — a decision for whoever picks up the shared-library fix,
not made unilaterally here.

### Step 8, Check 1 follow-up — top-port (port 4) FDTD-vs-hybrid comparison at 3 lengths (jobs 13069782-787, all COMPLETED 2026-09-10)

The earlier short-body Check 1 probe (mm_length=5.00) only ever excited ports
1 and 3 (the paper's own I={1,3}). At the user's request, added same-monitor-
locations FDTD-vs-hybrid engine-agreement data for the **top port (paper
port 4, `input_ys[-1]`)** — never driven through either script before this —
at mm_length = 5.00 / 10.00 / 15.00 µm. Both `mmi_4x4_90hybrid_fdtd.py` and
`mmi_4x4_90hybrid_hybrid.py` gained an `MMI4X4_DRIVEN_PORTS` env override
(comma list of 1-based paper ports, default `"1,3"` unchanged) — when the
driven set isn't exactly `{1,3}`, the paper's CMRR/Δφ metrics are skipped
(they need both driven inputs) and plain per-arm `|T|²`/phase(S) is reported
instead. Output directories carry a `_p{ports}` tag so nothing collides with
the canonical ports-1&3 files. The hybrid script's 4 input-section solves
always run regardless of the driven set (mirror reuse needs all 4), so the
L=5.00 case (job 13069785) was mostly a `save_dir` cache hit from job
13024689 — **26 s wall time**; L=10.00/15.00 (jobs 13069786/787) were fresh,
**6:55:30 / 6:58:45**. FDTD probes (13069782/783/784) were each fresh:
**54:32 / 1:35:00 / 2:52:02**.

**Per-arm `|T|²`, FDTD vs hybrid, at 1.29/1.31/1.33 µm:**

| L (µm) | port | FDTD `\|T\|²` | hybrid `\|T\|²` | max abs diff | max rel diff |
|---|---|---|---|---|---|
| 5 | out5 | 3.02e-5/8.31e-5/2.48e-5 | 2.02e-8/5.49e-8/3.19e-8 | 8.30e-5 | ~100% (both near noise floor) |
| 5 | out6 | 5.02e-6/7.11e-5/1.85e-4 | 6.52e-6/1.08e-5/1.54e-5 | 1.69e-4 | 30-92% (near noise floor) |
| 5 | out7 | 0.00527/0.00707/0.00514 | 0.00454/0.00518/0.00583 | 0.00189 | 27% |
| **5** | **out8 (dominant)** | 0.4472/0.4602/0.4444 | 0.4548/0.4544/0.4552 | **0.0108** | **2.4%** |
| 10 | out5 | 3.80e-6/4.84e-5/1.89e-4 | 3.22e-6/3.84e-6/5.52e-6 | 1.84e-4 | 15-97% (near noise floor) |
| 10 | out6 | 4.89e-4/8.47e-4/1.378e-3 | 7.86e-4/9.56e-4/1.176e-3 | 2.97e-4 | 61% |
| 10 | out7 | 0.05854/0.05839/0.06177 | 0.06425/0.06490/0.06478 | 0.00651 | 11.2% |
| **10** | **out8 (dominant)** | 0.5521/0.5560/0.5522 | 0.5356/0.5338/0.5319 | **0.0222** | **4.0%** |
| 15 | out5 | 6.39e-5/3.37e-4/6.07e-4 | 3.04e-4/3.93e-4/5.04e-4 | 2.41e-4 | 17-376% (near noise floor) |
| 15 | out6 | 0.01312/0.01882/0.02060 | 0.01509/0.01674/0.01842 | 0.00219 | 15% |
| 15 | out7 | 0.03560/0.04184/0.03510 | 0.03806/0.03836/0.03923 | 0.00412 | 11.7% |
| **15** | **out8 (dominant)** | 0.3771/0.3897/0.3203 | 0.3738/0.3619/0.3508 | **0.0306** | **9.5%** |

**The dominant output arm (out8) agrees well at every length** — max
absolute difference 0.011-0.031, relative 2.4-9.5%, growing mildly with
length (more mode content accumulated by L=15 than L=5). This is comparable
to or better than the earlier ports-1/3 probe's own agreement level at
these same short, off-self-imaging lengths, and is independent confirmation
across a THIRD physical port (this device's engine agreement does not look
port-specific). The weaker arms (out5/out6/out7) show larger relative
disagreement, but out5 in particular is dominated by numerical noise floor
(`|T|²` ~1e-5 to 1e-4) where a percent-level absolute difference between two
independently-converged near-zero numbers reads as tens to hundreds of
percent — not a meaningful disagreement.

**Phase(S) does NOT usefully compare at these lengths, and this is
expected, not a red flag.** E.g. out8's phase at L=15: FDTD
[-78.4°,-120.4°,-146.8°] vs hybrid [94.8°,60.9°,38.3°] — a ~170-185°
offset, and at L=5/L=10 the offset is neither constant nor small (tens to
~120°, varying with wavelength). None of these lengths are anywhere near
this device's self-imaging condition (3/4·L_π≈164 µm), so the ABSOLUTE
propagation phase accumulated between input and output is highly sensitive
to small path-length/dispersion differences between the two methods'
cascades — a difference of a fraction of a degree in effective index
compounds into tens of degrees of absolute phase over even a few µm at
these off-resonance lengths. The paper's own metrics (Δφ, Eq. 15) are
built entirely from *differences* between two excitations' phases for
exactly this reason (an absolute reference cancels out); a single-port
absolute phase was never a meaningful comparison target to begin with.
**Magnitude agreement (the table above) is the real signal here.**

**Unrelated to this numeric comparison**: the same-session attempt to build
a rigorous stitched-field-plot capability for the output taper FAILED its
own correctness check — see the next section for the full writeup. That
failure blocks only the field-plot/stitching work, not the numeric per-arm
comparison documented here.

### Step 8 field-plot work — output-taper capture mechanism FAILS its own correctness check (job 13069673, COMPLETED 2:25:24, HYPOTHESIS NOT CONFIRMED)

At the user's request for stitched hybrid + full-device FDTD field plots at
L=5/10/15 for the top-port excitation, added `FDTDSolver.solve()` a new
`capture_fields_dir` parameter (mirrors `solve_section()`'s own, reusing its
`_write_section_fields()` writer verbatim) so the full-device FDTD side can
dump a single z=0 in-plane field slice — straightforward, and unblocks the
FDTD-side field picture on its own.

The hybrid side is harder: the output taper is driven by a superposition of
up to 16 body modes, not one physical mode, so a real output-taper field
picture needs each body eigenmode launched INDIVIDUALLY from the interior
(body) side into the output taper — the reverse of every existing call on
this device (which always excite from the physical-port side, the "fewer
modes" side `solve_section()` picks by default). Forcing the reversed
direction needs `solve_section()`'s `excite_port`/`measure_port` `PortSpec`
override (both given bypasses the "fewer side" tie-break entirely) — a code
path exercised only once before in this repo
(`verify_ring_bend_reciprocity_3d.py`, a completely different port
topology: a bent ring coupler's x/y-normal ports, not two x-normal ports on
a straight MMI slab).

**Before committing to the full 16-mode capture across 3 lengths (a real
~40-70 cluster-hour cost), `mmi_4x4_90hybrid_output_capture_validate.py`
(NEW script) checked the mechanism cheaply**: launch just `N_TEST_MODES=3`
body modes from the interior at mm_length=5.00, and compare the resulting
transfer coefficients against the one free, exact ground truth already
available — reciprocity says this reversed call's own coefficients must
equal `M_in[TOP_PORT].transpose()`, already computed (and cached from job
13024689) by the production hybrid script's own mirror-reuse trick.

**Result: HYPOTHESIS NOT CONFIRMED.** At all 3 wavelengths × 3 test modes,
relative difference was **198.4-202.1%**, not scattered noise — e.g. at
1.31 µm mode 0: test=0.068235+0.191649j vs predicted=0.024942-0.190566j
(max abs diff 0.894). The two values' **magnitudes are comparable** (same
order, 10-20% apart) but their **phase/direction disagrees sharply** — the
signature of a sign/direction convention error (a wrong `PortSpec.sign`,
the wrong forward/backward `alpha_idx`, or a conjugate mismatch somewhere
in the reversed-excitation construction), not a magnitude/normalization
bug. Not yet root-caused.

**A real, unrelated bug was found and fixed along the way**: the
validation script's own `save_dir` pointed at the wrong, near-empty
`BEEF/example_sim_cache` (5 entries) instead of the real, actively-used
`example_sim_cache` at the repo root (1076 entries) — a copy-paste path
error, not a cache-mechanism problem. Fixed in the script for future runs;
the already-running job was left to finish rather than restarted, since it
still produced a valid (if slower than necessary) result.

**Status: the output-taper field-plot capability is BLOCKED**, pending a
root cause for the sign/direction mismatch above. The full-device FDTD
field capture (the `solve()` `capture_fields_dir` addition) is unaffected
and ready to use independently. Not attempted further without new
direction, per the user's own explicit instruction after this result.

### ❌ n_eff re-calibration at the CURRENT geometry BLOCKED by the same Check 2 guard bug (job 13092881, ALL 16 TASKS FAILED, 2026-09-11)

Attempted to re-measure this device's only prior n_eff calibration (job
12078959, "does not rescue this device") at the CURRENT optimised design
point (mm_width=8.5, N_MODES=16, mm_length=159.30) — the old result was
measured at the superseded mm_width=6.0 geometry and has never been
checked against the point actually being validated. Reused the proven
one-process-per-mode architecture from the 2026-08-20 campaign
(`mmi_4x4_90hybrid_calib2_mode_worker.py`, a 16-task SLURM array +
`..._calib2_aggregate.py`, adapted only for the new geometry constants and
`decay_tol=1e-4`).

**All 16 array tasks FAILED identically, ~50 min each**, before reaching
any real FDTD calibration work at all — in `section_mode_neff()`'s own
16-mode seed step (called once per task, `solver.section_mode_neff(mmi,
wl_seed, ..., n_modes=16, ...)`):

```
RuntimeError: FDTDSolver: could not find 16 distinct eigenmodes within 24
band_num attempts at this cross-section -- meep/MPB returned a duplicate
more times than this guard tolerates.
```

This is the **identical failure mode** diagnosed under Step 8 Check 2
above (job 13060190) — and the identical walk shape: band_num 13 through
23 (11 consecutive attempts) were rejected as duplicates in every one of
the 16 tasks, deterministically (all 16 solve the same seed cross-section),
before band_num=24 (the last allowed attempt) also failed to close out 16
accepted bands. This is strong independent confirmation that the false
positives found in the Check 2 diagnostic (the guard's fixed 16×8
diagnostic grid aliasing this wide 8.5 µm slab's higher-order mode lobes)
are not an audit-tool-only problem — **the same guard blocks any
`section_mode_neff()`/`section_mode_beta()` call requesting 16 modes at
this cross-section, including calibration**, which calls it internally via
`section_mode_beta_details()`'s own seed step.

**Consequence: the n_eff re-calibration cannot proceed at the current
geometry until the shared-library fix already flagged under Check 2 is
made** (scaling `_DUP_BAND_DIAG_NY`/`_DUP_BAND_DIAG_NZ` to the requested
mode count / cross-section width, rather than a fixed 16×8 constant). The
downstream aggregator (job 13092882) and with/without comparison table
(job 13092912) never ran — both were queued on `--dependency=afterok` and
that condition cannot be satisfied by a failed array.

**Status.** The only calibration data this device has is still the
2026-08-20 result at the superseded mm_width=6.0 geometry ("calibration
does not rescue this device": phase +0.150° worse, CMRR +2.547 dBe better
but 16 dBe short of spec either way). Whether that conclusion holds at the
current mm_width=8.5 design remains genuinely unknown, and re-checking it
is now blocked on the same shared-library item Check 2 already surfaced —
fixing that guard would unblock both.

### ✅✅ Shared-library fix VERIFIED — Check 2 now passes cleanly, calibration pipeline unblocked (job 13096902, COMPLETED 38:47, 2026-09-11)

The `beef/solvers/fdtd.py` `_iter_distinct_bands()` fix (fixed-spacing
0.05 µm diagnostic grid, replacing the fixed 16×8 count that produced the
false positives documented above) was verified against the exact audit
that originally failed as job 13060190.

**Result: `audit_interface_modes.py mmi_4x4_90hybrid` now completes
cleanly.**

```
guided modes: 16/16
correspondence (guided, matched res): [0, 1, 2, ..., 15]   -- identity
ORDER PRESERVED: True
worst guided |dn_eff|: 0.0000
```

**Zero duplicate-skip warnings occurred anywhere in the run** — not just
"below threshold," but none at all, across all 16 meep bands. This is a
stronger result than merely "the guard no longer trips": the underlying
n_eff agreement between meep and meow at this cross-section is essentially
exact (worst |dn_eff| rounds to 0.0000 at the printed precision), matching
the clean picture the fine-grid field dump already established (Step 8
Check 2's earlier sections) — this device's real interface basis was never
the problem; the audit tool's own coarse grid was.

**Step 8 Check 2: fully CLOSED**, both device side and tool side. No open
items remain for mode identity/ordering on this device.

**Consequence for calibration**: the n_eff re-calibration pipeline (jobs
13092881/13092882/13092912, blocked by this exact bug in the previous
section) is now unblocked — the same `section_mode_neff()` 16-mode seed
call that failed there uses the identical, now-fixed guard. Not
relaunched yet; ready whenever wanted.

**One earlier data point revised in light of this**: job 13096233 (an
interim resubmit attempt) completed in 8 s with zero output — traced to a
one-off scheduling/node anomaly (a second resubmit, job 13096902, above,
ran normally), not a code problem. Noted here only so a future reader
doesn't mistake it for a real regression.

### ✅✅ n_eff calibration at the CURRENT geometry: RE-MEASURED, and it makes CMRR WORSE (jobs 13098452/453/454, all COMPLETED, 2026-09-11)

With the shared-library duplicate-guard fix verified (previous section),
relaunched the identical calibration pipeline that failed pre-fix (jobs
13092881/882/912). All 16 array tasks completed cleanly this time
(~1:17-1:38 h each for modes 1-15, 2:07 h for mode 16 — slower, more MPB
iterations per solve as expected for the mode closest to the guided-mode
cutoff, but no instability recurrence: its own FDTD decay converged
cleanly to ratio 1.3e-9, nothing like the old undamped-Lorentzian blow-up).
Zero duplicate-mode warnings across all 16 tasks.

**Aggregator guards passed cleanly**: worst residual 5.923e-3 rad (guard
0.1), worst branch margin 0.063 (guard 0.9), worst seed excess +1.71%
(warn 10%) — the calibration is real and usable, not a guard-bypassed
artifact.

**With-vs-without comparison, all 4 lengths this session has probed:**

| L (µm) | max\|Δφ\| uncal | max\|Δφ\| cal | min CMRR uncal | min CMRR cal | max XL uncal | max XL cal | max imb uncal | max imb cal |
|---|---|---|---|---|---|---|---|---|
| 5.00 | 119.052° | 118.968° | 0.000 dBe | 0.000 dBe | 3.477 dB | 3.488 dB | 73.529 dB | 73.509 dB |
| 10.00 | 176.296° | 179.098° | 0.000 dBe | 0.000 dBe | 4.293 dB | 4.295 dB | 52.214 dB | 52.308 dB |
| 15.00 | 143.328° | 143.943° | 0.014 dBe | 0.014 dBe | 5.095 dB | 5.091 dB | 30.890 dB | 30.811 dB |
| **159.30 (production)** | **3.933°** | **3.958°** | **18.482 dBe** | **15.615 dBe** | **2.706 dB** | **1.924 dB** | **1.972 dB** | **2.064 dB** |

The three short lengths (5/10/15 µm) are off-self-imaging probe points, not
device operating conditions (established earlier this session) — their
huge phase/imbalance numbers either way are expected and calibration moves
them by only a few tenths of a degree/dB, noise-level next to numbers that
large.

**The number that matters is L=159.30, the production design point.**
Calibration there:
- max|Δφ| 3.933° → 3.958° (+0.025°, essentially unchanged, still meets
  the ≤5° spec either way)
- **min CMRR 18.482 → 15.615 dBe (−2.867 dBe, WORSE)** — moves further
  from the ≥20 dBe spec, not closer
- max excess loss 2.706 → 1.924 dB (better)
- max imbalance 1.972 → 2.064 dB (slightly worse)

**This both confirms and sharpens the 2026-08-20 old-geometry finding.**
At the superseded mm_width=6.0 point, calibration *improved* CMRR (+2.547
dBe, still 16 dBe short) while making phase slightly worse. At the
CURRENT mm_width=8.5 optimized point, calibration makes CMRR **worse**
while phase stays flat — the opposite CMRR direction from before. Taken
together, the two results agree on the one thing that matters for this
device's verdict: **calibration does not rescue it, at either geometry,
and the CMRR/phase shortfall is real physics of the design, not a
calibration artifact of the uncalibrated MPB-derived n_eff.** The
CMRR-direction reversal between the two geometries is itself informative:
it shows the calibration correction's sign/magnitude is geometry-specific
and not a fixed, predictable nudge — reinforcing that it should be
re-measured per design point rather than assumed to transfer, exactly as
was done here.

**Status: this device's n_eff calibration question is now CLOSED at the
current design point**, superseding the 2026-08-20 old-geometry-only
result as the operative answer for anything referencing "has this device
been calibrated."

### ⚠⚠⚠ ROOT CAUSE FOUND — the hybrid-vs-FDTD phase drift is the body's propagation constant, and it is fully quantified (2026-09-11)

**The observation.** Comparing per-arm `phase(S)` between the hybrid
cascade and the full-device FDTD for the top port (paper port 4) at
mm_length = 5/10/15 µm — data from jobs 13069782-787, already tabulated
above — the phases do NOT match, and the mismatch is not noise:

| L (µm) | hybrid − FDTD phase offset, out8 (dominant arm) | all four arms? |
|---|---|---|
| 5 | −56.4° / −57.8° / −51.7° | scattered, out8 tightest |
| 10 | −120.0° / −118.0° / −114.7° | **yes — all 4 arms cluster −100° to −121°** |
| 15 | +173.2° / −178.7° / −174.9° | **yes — all 4 arms cluster near ±180°** |

Two features make this diagnostic rather than merely noisy: the offset is
nearly **wavelength-independent** at fixed length (spread ~5° across
1.29-1.33 µm), and at L=10/15 it appears **coherently on all four output
arms at once**. Random per-arm error cannot do that.

**The mechanism, confirmed quantitatively.** A wavelength-independent
offset growing linearly with length is the signature of a constant
propagation-constant mismatch, `Δφ = Δβ·L`. Back-solving the ~−13.6°/µm
slope gave Δn_eff ≈ 0.05. That prediction was then confirmed against data
already on disk, with no new simulation: **the calibration (jobs
13098452/453) measures n_eff +0.0458 higher than the MPB/meow mode
solver**, and — critically — that offset is essentially **mode-independent**
(+0.04558 to +0.04604 across all 16 modes; spread 0.00046, under 1%):

| mode | mode solver n_eff (1.31) | calibrated n_eff (1.31) | Δ |
|---|---|---|---|
| 1 | 2.963353 | 3.009390 | +0.04604 |
| 8 | 2.899712 | 2.945469 | +0.04576 |
| 16 | 2.696700 | 2.742437 | +0.04574 |

Closing the loop arithmetically, with `L_body = L − 2·body_stub`:

| L (µm) | L_body | predicted Δφ = (2π/λ)·0.0458·L_body | observed | ratio |
|---|---|---|---|---|
| 5.00 | 4.50 | 56.6° | 55.3° | **0.977** |
| 10.00 | 9.50 | 119.5° | 117.6° | **0.984** |
| 15.00 | 14.50 | 182.4° | 184.4° | **1.011** |

**Agreement within 2.3% at every length.** The phase discrepancy is not
partially explained or plausibly attributed — it is *fully accounted for*,
to the precision of the measurement, by the difference between the n_eff
the mode solver reports and the n_eff meep's own time-stepper actually
realizes.

**Why this does NOT contradict Check 2.** Check 2 found meep-vs-meow
*mode-solver* agreement to be essentially exact (worst |dn_eff| 0.0000).
That is a statement about two EIGENSOLVERS agreeing with each other. The
+0.0458 here is a different comparison entirely: eigensolver-vs-
time-stepper. Both can be true at once, and both are. The uncalibrated
hybrid propagates its body at the eigensolver's index; the full-device
FDTD propagates at the time-stepper's. Nothing in Check 2 could have
caught this.

**Leading mechanism for the +0.0458 itself: Yee-grid numerical
dispersion.** Sign is right (discrete phase velocity runs slow, so the
realized index reads HIGH), mode-independence is right (a grid property,
not a mode property), and the magnitude is the right order: 0.0458/2.96 =
1.55% at ~15 points per wavelength in-material (res=40, λ/n ≈ 0.374 µm).
Second-order Yee dispersion predicts ~0.7% there — same order, factor ~2
short, with the gap plausibly from subpixel averaging and the real
3-D dielectric structure. NOT yet confirmed; see the debug plan below.

**Implication.** `calibrate_body_neff` exists to measure exactly this
quantity, and the numbers say it should very nearly eliminate the drift.
That is being tested directly (job 13102478, cascade arithmetic on cached
inputs only). Note the tension with the calibrated production-length
result recorded above (CMRR 18.482 → 15.615 dBe, worse): if calibration is
the physically correct propagation constant, the calibrated number is the
more trustworthy one, and the *uncalibrated* 18.482 dBe would be the
artifact — which would make the device's real CMRR worse than the headline
figure this campaign has been quoting. That question is now the single
highest-value open item on this device.

### ✅✅✅ CONFIRMED — calibration collapses the phase drift, and that forces a downgrade of the device's headline CMRR (job 13102478, 2026-09-11)

Debug-plan step 1, run on cached inputs only (cascade arithmetic, 5.7 s of
actual work). Re-ran the top-port cascade at L=5/10/15 with
`neff_override=neff_calibrated` and compared per-arm `phase(S)` against the
full-device FDTD reference:

| L (µm) | mean \|Δφ\| uncalibrated | mean \|Δφ\| **calibrated** |
|---|---|---|
| 5 | 46.71° | 32.40° |
| 10 | 111.01° | **15.02°** |
| 15 | 167.02° | **11.03°** |

Per-arm, the arms that actually carry power are fixed almost exactly:

| arm | L=5 | L=10 | L=15 |
|---|---|---|---|
| **out8** (dominant) | 55.3° → **1.9°** | 117.6° → **1.8°** | 175.6° → **2.2°** |
| **out7** | 51.0° → 11.8° | 117.0° → **1.9°** | 177.0° → **2.3°** |
| out6 | 48.0° → 55.6° | 102.6° → 15.5° | 171.8° → 5.8° |
| out5 | 32.6° → 60.2° | 106.8° → 40.8° | 143.8° → 33.7° |

**Up to an 80x improvement on the dominant arm** (175.6° → 2.2° at L=15).

The residual tracks arm POWER, not arm index. Pairing each phase error with
the full-device FDTD |T|² actually measured in that arm (mean over
1.29/1.31/1.33 µm; jobs 13069782/783/784) makes the pattern unambiguous:

| L (µm) | arm | FDTD \|T\|² | \|Δφ\| uncal | \|Δφ\| **cal** |
|---|---|---|---|---|
| 5 | out5 | 4.60e-05 | 32.6° | 60.2° |
| 5 | out6 | 8.70e-05 | 48.0° | 55.6° |
| 5 | out7 | 5.83e-03 | 51.0° | 11.8° |
| 5 | out8 | 4.51e-01 | 55.3° | **1.9°** |
| 10 | out5 | 8.05e-05 | 106.8° | 40.8° |
| 10 | out6 | 9.05e-04 | 102.6° | 15.5° |
| 10 | out7 | 5.96e-02 | 117.0° | **1.9°** |
| 10 | out8 | 5.53e-01 | 117.6° | **1.8°** |
| 15 | out5 | 3.36e-04 | 143.8° | 33.7° |
| 15 | out6 | 1.75e-02 | 171.8° | **5.8°** |
| 15 | out7 | 3.75e-02 | 177.0° | **2.3°** |
| 15 | out8 | 3.62e-01 | 175.6° | **2.2°** |

Sorted by power the residual is monotone: every arm above |T|² ≈ 1e-2
lands at 1.8-5.8°, every arm below 1e-3 stays at 33-60°, and the one
in-between arm (out7 at L=5, 5.8e-3) lands in between at 11.8°. out5 at
L=5 getting "worse" (32.6° → 60.2°) is not a counterexample: at |T|²
= 4.6e-5, ~4 orders below the dominant arm, both numbers are the phase of
a near-zero field and neither is a physical quantity. **Every arm carrying
real power agrees to ~2°.** Note this corrects a looser earlier reading of
this data — out6 at L=15 is at 1.75e-2, NOT down in the noise, and it is
fixed to 5.8°, so "the residual is all in out5/out6" was the wrong framing;
"the residual is all in the low-power arms" is the right one.

Figure: `example_output/verification/mmi_4x4_90hybrid_calib2_phase_check/calib_phase_vs_power.png`
(script: `examples/verification/python/plot_mmi4x4_90hybrid_calib_phase_vs_power.py`)
— per-length grouped bars plus the error-vs-power scatter on a log power axis.

**This closes the mechanism completely.** The hybrid-vs-FDTD phase
disagreement at reduced lengths was never a cascade/interface/sectioning
problem, and never a group-index problem. It was one number: the body's
propagation constant, where the uncalibrated hybrid used the mode solver's
n_eff and the full-device FDTD used what meep's time-stepper actually
realizes on its grid, 0.0458 higher. Correct that one number and the two
methods agree to ~2° over 15 µm of propagation.

**The consequence is NOT good news for the device.** Calibration is now
demonstrated to be the physically correct propagation model on this
cross-section — it is the setting that reproduces the full-device FDTD,
which is the most complete simulation available. That makes the
*calibrated* production-length numbers the trustworthy ones:

| metric | uncalibrated (previously reported) | **calibrated (now the defensible figure)** | spec |
|---|---|---|---|
| min CMRR | 18.482 dBe | **15.615 dBe** | ≥ 20 dBe |
| max \|Δφ\| | 3.933° | 3.958° | ≤ 5° |

**The device's real CMRR shortfall is ~4.4 dBe, not ~1.5 dBe.** Every
CMRR number this campaign has quoted at the optimised design point —
including the headline "18.482 dBe, 1.52 short of spec" that the whole
rounds-8-10 optimisation was steered by — was computed with a body index
0.0458 too low. The phase spec still passes either way (3.933° vs 3.958°,
both under 5°), so the *phase* conclusions survive; the *CMRR* conclusions
do not.

**What this does NOT yet establish**: whether the +0.0458 is genuine
physics of a correctly-converged simulation or an artifact of running at
res=40 (leading hypothesis: Yee-grid numerical dispersion, which would
shrink as ~1/res² and vanish in a converged run). If it is grid
dispersion, then BOTH methods are wrong at res=40 and the true device CMRR
is neither 18.482 nor 15.615 — it is whatever the res→∞ limit gives, and
the calibrated number is merely the one that correctly reproduces the
res=40 FDTD. Debug-plan step 3 (resolution sweep, jobs 13102582/583/584)
is what settles this, and it is now the single highest-value open item on
this device.

**Optimisation-campaign caveat, stated explicitly.** The rounds 8-10
(mm_width, taper_width, taper_length, mm_length) campaign that produced
the current design point evaluated ~45 points through the UNCALIBRATED
hybrid path. If the calibration offset shifts CMRR by ~-2.9 dBe
non-uniformly across geometry, the campaign's ranking of candidate points
is not guaranteed to survive re-scoring. Whether the chosen optimum is
still the optimum under the corrected propagation constant is an open
question this finding raises and does not answer.

---

### ✅ ANSWERED — that open question, answered the same day: the optimum DOES move, and the device MEETS the CMRR spec at the moved optimum (2026-09-11)

The paragraph immediately above asks whether the campaign's chosen optimum
survives re-scoring under the corrected propagation constant. It does not —
and the correction is in the device's **favour**.

`examples/verification/python/mmi_4x4_90hybrid_calibrated_length_rescore.py`
re-scores the entire `mm_length` axis (7601 points, 20.00–400.00 µm @ 0.05 µm,
3 wavelengths) through both the uncalibrated and the calibrated body index,
with the scoring code replicating `mmi_4x4_90hybrid_hybrid.py` exactly (same
`SMatrix` methods, same A11 phase convention, same Eq. 15 targets).

**This cost zero new simulations.** The calibration is a property of the body
CROSS-SECTION, not its LENGTH, so everything needed was already cached: `M_in`
for all four ports (length-independent — the taper/stub geometry does not
change when the body grows) and the calibrated body length-sweep grid solved
by job 13098454. Total runtime 13.7 s, of which 11.2 s was loading the two
93 MB body grids off disk.

| | uncalibrated | **calibrated** |
|---|---|---|
| optimal `mm_length` | 159.30 µm | **158.80 µm** |
| min CMRR | 18.482 dBe | **20.326 dBe** |
| max \|Δφ\| | 3.933° | **3.664°** |
| excess loss | 2.706 dB | **2.444 dB** |
| meets CMRR ≥ 20 dBe | ✗ | **✓** |

**The script validates itself against production on two independent points:**
the uncalibrated path reproduces the campaign's headline 18.482 dBe / 3.933° /
2.706 dB at 159.30 µm exactly, and the calibrated path reproduces the 15.615
dBe / 3.958° reported in the section above at that same length. Identical to
three decimals on both, so the re-scoring machinery is provably the same
scoring the campaign used.

**The 4.4 dBe "shortfall" was a wrong LENGTH, not a worse device.** Scoring the
calibrated device at the *uncalibrated* optimum (159.30 µm) gives 15.615 dBe;
moving 0.50 µm to the calibrated optimum gives 20.326 dBe. The earlier
"calibration downgrades the device" reading in the section above was correct
about the number but wrong about the cause — it compared a calibrated device
to a length chosen for an uncalibrated one. **Every headline CMRR figure in
this document that predates this entry should be read as superseded.**

The shift direction was predictable but its MAGNITUDE was not. To first order
the optimum should move to hold `β·L` invariant, giving
`159.30 × 2.9632/3.0094 ≈ 156.85 µm`, a −2.4 µm shift. The measured shift is
−0.50 µm, about **5× smaller** than that estimate. The naive scaling treats
the body as a single-mode phase accumulator; the self-imaging condition
actually depends on *beat lengths between mode pairs*, i.e. on
`Δβ_mn`, and the calibration offset is nearly uniform across modes (spread
0.00046 over 16 modes) so it largely cancels in the differences. That near-
cancellation is why the device survives at all — a mode-DEPENDENT offset of
the same size would have moved the optimum ~5× further and degraded it.

Figure: `example_output/verification/mmi_4x4_90hybrid_calibrated_length_rescore/calibrated_length_rescore.png`

**Caveat that governs this entire entry: this is res = 40, and res = 40 is not
converged** — see the resolution entry below. The 20.326 dBe figure is "the
number that reproduces the res = 40 full-device FDTD", which is the best
available reference but not yet the converged answer.

---

### ✅ POWER-WEIGHTED phase agreement: 2.12°, and out8 alone owns ~85% of it (2026-09-11)

`examples/verification/python/mmi_4x4_90hybrid_weighted_phase_metric.py`.
The unweighted mean over the four arms is the wrong figure of merit at the
reduced probe lengths, because two of the four arms sit 3–4 decades down in
power and contribute the phase of numerical noise. Weighting each arm by its
share of output power, `Σ_j (|T_j|²/Σ|T|²)·|Δφ_j|`:

| L (µm) | uncal | **calibrated** | % of 360° | % of 90° |
|---|---|---|---|---|
| 5 | 55.25° | **2.08°** | 0.578 | 2.31 |
| 10 | 117.49° | **1.86°** | 0.517 | 2.07 |
| 15 | 175.49° | **2.43°** | 0.674 | 2.69 |
| **all** | **116.08°** | **2.12°** | **0.590** | **2.36** |

Against a "below 1.5%" target the answer depends entirely on the denominator,
which must therefore be stated: **1.5% of 360° = 5.40°, already met with 2.5×
margin; 1.5% of 90° = 1.35°, NOT met (2.12°).** The 90° reading is the
meaningful one for a quadrature hybrid — it is the scale the ≤ 5° spec is
implicitly set against — so the target is read here as **1.35°**, and the gap
to close is 2.12° → 1.35°, a factor of 1.6.

**Where the residual lives — the residual budget, after weighting:**

| L (µm) | out5 | out6 | out7 | **out8** |
|---|---|---|---|---|
| 5 | 0.3% | 0.6% | 7.1% | **92.0%** |
| 10 | 0.2% | 1.2% | 10.3% | **88.3%** |
| 15 | 1.0% | 11.0% | 8.6% | **79.3%** |

This inverts the priority implied by the raw per-arm table. The 33–60°
residuals in out5/out6 are **irrelevant to the weighted metric** — they
contribute under 1.5% combined at two of three lengths. **Reaching 1.35°
means improving out8 from ~1.9° to ~1.4°, and nothing else matters.** Any
effort spent on the low-power arms is wasted against this target.

---

### ✅ STEPS 3/4/5 COMPLETE — grid dispersion CONFIRMED as the mechanism, group index and bandwidth RULED OUT (jobs 13102582–589, 2026-09-11)

All eight parameterised probes (`mmi_4x4_90hybrid_calib_probe.py`) landed.
Mode 1, n_eff at 1.31 µm.

**Step 5 — pulse bandwidth: RULED OUT, definitively.**

| bandwidth factor | 1.5 | 2.0 | 4.0 |
|---|---|---|---|
| offset | +0.046181 | +0.046181 | +0.046181 |

Identical to **six decimal places**. The fitted β does not depend on pulse
bandwidth at all, so envelope/group-velocity contamination of the phase fit
is not merely small — it is absent. Combined with the magnitude argument
(n_g − n_eff = 0.89, 20× the offset), **the group index is definitively not
involved**. This closes the "do we need to incorporate the group index"
question with a measured null.

**Step 4 — source band and centre: a real but SECOND-ORDER effect.**

| band (µm) | 1.29–1.33 | 1.26–1.36 | 1.25–1.31 | 1.31–1.37 |
|---|---|---|---|---|
| centre | 1.31 (on) | 1.31 (on) | 1.28 (off) | 1.34 (off) |
| offset | +0.046181 | +0.044989 | +0.044507 | +0.044964 |

Spread 0.0017, i.e. **3.6% of the offset** — so band choice is not the
mechanism. But it is not negligible in absolute terms either: 0.0017 of n_eff
over a 159 µm body is ≈ 74° of accumulated phase at 1.31 µm. **The full-device
FDTD and the hybrid's calibration therefore have to agree on their source
band, and nothing in the code enforces that** — they happen to both use
1.29–1.33 today. This is an unguarded assumption and should become an assert
(see open items). Notably, moving the band centre OFF the monitored
wavelength changes the offset by less than widening the band does, so the
`dispersive=False` band-centre material freeze is not the dominant term.

**Step 3 — resolution: the mechanism, CONFIRMED, and it is not cheap to escape.**

| res | n_eff mode solver | n_eff calibrated (FDTD) | offset | β-fit resid |
|---|---|---|---|---|
| 20 | 2.916967 | 3.037510 | +0.120543 | 5.72e-02 |
| 40 | 2.963209 | 3.009390 | +0.046181 | 5.58e-03 |
| 60 | 2.974730 | 2.998948 | +0.024218 | 1.43e-04 |

The offset falls monotonically and fast: ratio 2.61 over 20→40 and 1.91 over
40→60, giving pairwise exponents p = 1.38 and 1.59 and a three-point power-law
fit of **p ≈ 1.45** for `offset ∝ res^-p`. Not the clean p = 2 a pure Yee
second-order argument predicts, but unambiguously a discretisation error
converging to zero.

**The most important detail is that BOTH engines move, in OPPOSITE
directions.** The mode solver climbs (2.9170 → 2.9632 → 2.9747) while the FDTD
time-stepper falls (3.0375 → 3.0094 → 2.9989). This is **not** a case of one
engine being right and the other wrong, and the earlier framing in this
document — "the FDTD realises a different index than the mode solver" — was
half the story. They are both under-resolved at res = 40 and are converging
toward each other from opposite sides, somewhere near n_eff ≈ 2.98–2.99.

Richardson extrapolation is offered with explicit low confidence: the mode
solver's three points fit p = 1.50 almost exactly and extrapolate to
n_∞ ≈ 2.9885; the FDTD's fit a much softer p ≈ 0.82 and extrapolate to
n_∞ ≈ 2.9727. **These do not meet** (0.016 apart), but a three-parameter fit
to three points cannot distinguish a genuine exponent from a stalling ratio —
exactly the failure mode recorded for the HCG work. Treat the extrapolation
as a hypothesis needing a res = 80 point, not a result.

**The practical consequence is decisive: brute-force resolution convergence is
not reachable for this device.** Extrapolating `offset ≈ 9.7·res^-1.45`,
driving the offset below 0.005 needs **res ≳ 180** — in a 3D simulation of an
8.5 × 159 µm body that is far beyond feasible. **Calibration is therefore not
a stopgap to be removed once resolution improves; it is the mechanism by which
this device can be simulated at all at achievable resolution.** That upgrades
n_eff calibration from "a correction we discovered we needed" to "a required
part of the method for any wide multimode body", and it applies to every other
MMI in this repo, not just this one.

---

### 📋 PLAN — driving the power-weighted phase error from 2.12° to below 1.35°

**The residual is systematic, not noise — measured, not assumed.** Signed
calibrated residual on out8 (the arm that owns ~85% of the weighted metric):

| L (µm) | 1.29 | 1.31 | 1.33 | mean | std |
|---|---|---|---|---|---|
| 5 | +1.68 | −0.89 | +3.26 | **+1.35** | 1.71 |
| 10 | +2.28 | +1.96 | +1.22 | **+1.82** | 0.44 |
| 15 | −0.09 | +4.49 | +2.12 | **+2.18** | 1.87 |

The mean is positive at all three lengths and **grows monotonically with
length** — the signature of a residual Δβ, not of random error. Slope
≈ 0.083°/µm, i.e. a **residual Δn_eff ≈ 3.0e-4**. Calibration removed
0.0458 − 0.0003, or **99.3% of the original offset**; what is left is 0.66% of
it. Superimposed is a per-wavelength scatter of ~0.4–1.9° std, comparable in
size to the mean and estimated from only three wavelengths, so it is itself
poorly determined.

**Both components must be attacked, and only out8 matters.** The residual
budget shows out5/out6 contribute under 1.5% combined at two of three lengths;
their 33–60° errors are irrelevant to this target.

**Knobs, ranked by measured leverage** — ⚠️ **this ranking has since been
OVERTAKEN BY MEASUREMENT; read the dated entries below before acting on it.**
Knobs 1, 3 and 4 have all been run and closed (1 and 3 negative, 4 clearing
the way), and the leading candidate is now one that is not in this table at
all — the number and span of the CALIBRATION LENGTHS. The table is kept as the
record of what was believed before the runs, not as current guidance:

| # | knob | expected effect | evidence | cost |
|---|---|---|---|---|
| 1 | **calibration at res = 60** | attacks the systematic 3e-4 | β-fit residual falls 5.58e-3 → 1.43e-4, **39×**, measured in the res probe | 16 modes × the res-60 probe |
| 2 | **more calibration wavelengths (3 → 5/7)** | attacks the scatter; 3 points fit a quadratic with zero residual, so any real dispersion curvature is invisible | scatter std 0.4–1.9° estimated from 3 samples — statistically meaningless | 16 modes, ~2× the 3-λ cost |
| 3 | **N_MODES 16 → 20/24** | basis truncation at the taper→body interface; never tested at this geometry | untested — this is the largest *unmeasured* term | new EME body grid + new calibration |
| 4 | **FDTD reference `decay_tol` 1e-4 → tighter** | lowers the REFERENCE's own noise floor, which bounds what any hybrid-side fix can achieve | not yet bounded; note `decay_tol` has a measured hard floor well above 1e-11 | re-run 3 full-device FDTDs |
| 5 | **`BODY_STUB` 0.25 µm** | interface placement between FDTD and EME | untested sensitivity | cheap sweep on cached inputs |
| 6 | **assert band match** | correctness guard, not accuracy | band spread 0.0017 ≈ 74° over 159 µm; nothing currently enforces that calibration and production share a band | trivial |

Knob 1 is first because it is the only one whose benefit is already measured
rather than hypothesised, and it targets the component proven to be
systematic. Knob 4 must be done early regardless: **if the reference's own
phase noise is ~1°, then 1.35° is at or below the measurement floor and the
target is unachievable as stated** — that has to be known before spending on
knobs 1–3.

**Re-optimisation status.** `mm_length` is the ONE geometry axis that
re-optimises for free under calibration (done — see the rescore entry above:
158.80 µm, 20.326 dBe, spec met). Every other axis — `mm_width`,
`taper_width`, `taper_length`, `waveguide_width`, `waveguide_offset` — changes
the body CROSS-SECTION and therefore invalidates the calibration, so each
candidate point needs its own 16-mode calibration run before it can be scored
honestly. 42 evaluated campaign points exist under `wg_ygap_opt/`, but they
cannot simply be re-scored: the calibration is per-cross-section and none of
them have one.

The cheap path worth testing, stated as the approximation it is: the offset is
a **discretisation artifact** (confirmed by the res sweep) and is nearly
uniform across modes (spread 0.00046 over 16). If it is also weakly dependent
on cross-section geometry, a **single measured offset could be applied as a
uniform shift** to the mode-solver n_eff of other geometries, making the whole
campaign re-scorable for the cost of validating that assumption on two points.
**This is an untested hypothesis, and a plausible failure mode is specific:**
the offset should scale with how many grid cells span the guided field, so a
narrower body would show a larger offset. Validate by running the full 16-mode
calibration at two campaign geometries bracketing the current `mm_width` and
comparing against the uniform-shift prediction; only adopt the shortcut if
those agree.

**Queued 2026-09-11 — knobs 1 and 4, the two that gate the rest:**

| job | what | budget | batch file |
|---|---|---|---|
| **13110658** (array 1-16) | knob 1 — full 16-mode calibration at **res = 60** | 12 h/task | `run_mmi_4x4_90hybrid_calib_res60.batch` |
| **13110659** | knob 4 — L=15 top-port full-device FDTD at **decay_tol = 1e-5** | 15 h | `run_mmi_4x4_90hybrid_fdtd_tol1e5.batch` |

Knob 4 is queued alongside knob 1 rather than after it because it can
*invalidate the target*: `decay_tol = 1e-4` was selected on a `|T|²`
AMPLITUDE criterion (0.015 dB), which says nothing about phase. If the
reference's own phase noise turns out to be ~1°, then 1.35° is at or below the
measurement floor and knobs 1–3 cannot reach it no matter what they cost. L=15
is the probe case deliberately — largest calibrated residual (2.43° weighted)
and the longest propagation, so reference noise has the most room to show up.

**One real hazard fixed before submitting.** `MMI4X4_FDTD_DECAY_TOL` was
already an env override in `mmi_4x4_90hybrid_fdtd.py`, but `output_dir` was
tagged only by `mm_length` and the driven-port set — so a tolerance probe
would have written straight over job 13069784's reference results, the exact
failure this device already had a near-miss on when coarse smoke runs landed
in the official `output_dir`. Added `_decay_tol_tag` (empty at the 1e-4
default, so no existing path changes); the probe now lands in
`mmi_4x4_90hybrid_fdtd_L15.00_p4_tol1e-05`.

Aggregation for the res=60 family is
`examples/verification/python/mmi_4x4_90hybrid_calib_probe_aggregate.py`
(`MMI4X4_CALIB_TAG=res60cal`). It refuses to run on a partial set, and writes
to `calib_probe/aggregated_res60cal_n16.npz` rather than into
`calib2_per_mode/` — the res=40 calibration is what every current result was
produced with, and replacing it has to be a deliberate act, not a file landing
in the right place. Per-mode independence makes the split exact: the library's
own post-fit step fits each mode from only that mode's data.

Standing directive applies to both: **if either times out, report and wait —
do not auto-resubmit.**

---

### ✅ KNOB 4 ANSWERED — the FDTD reference's own phase noise is 0.029°, so the 1.35° target is NOT floor-limited (job 13110659, 2026-09-12)

L=15, top port, `decay_tol` 1e-4 → 1e-5 (4:00:49). Per-arm phase shift:

| arm | \|T\|² | Δphase 1.29 / 1.31 / 1.33 |
|---|---|---|
| **out8** | 3.6e-01 | −0.043 / −0.011 / −0.005 |
| out7 | 3.8e-02 | +0.123 / +0.041 / −0.076 |
| out6 | 1.7e-02 | −0.307 / +0.054 / −0.005 |
| out5 | 3.4e-04 | +1.237 / +0.028 / −0.728 |

**Power-weighted reference phase uncertainty: 0.029°.** Worst among arms
carrying real power: 0.307°. The dominant arm out8 — which owns ~85% of the
weighted metric — moves by at most 0.043°.

**The 1.35° target sits ~46× above the reference's own noise floor, so it is
not floor-limited and knobs 1–3 are worth spending on.** This also extends the
earlier `decay_tol` finding: 1e-4 was selected on a `|T|²` AMPLITUDE criterion
(0.015 dB) and is now confirmed adequate for PHASE as well, which had never
been checked. (out5's 1.2° swing at 1.29 µm is a 6e-5-power arm — noise, and
weighted to irrelevance.)

---

### ❌ KNOB 1 ANSWERED, NEGATIVELY — a BETTER-converged calibration makes agreement WORSE, and the amount is exactly predictable (job 13110658, 2026-09-12)

All 16 tasks COMPLETED (1:02–2:44 each, far under the 12 h budget — the res=60
mode-solver seed was already cached from the single-mode probe). The
calibration itself is excellent, exactly as the probe predicted:

| | res = 40 | **res = 60** |
|---|---|---|
| offset vs mode solver @1.31 | +0.0462 | **+0.0244** |
| offset spread over 16 modes | 0.00046 | 0.00057 |
| β-fit residual (mean) | 5.58e-03 | **1.27e-04** |

**A 44× better β fit — and the measured hybrid-vs-FDTD agreement got 18×
WORSE:**

| arm | L=5 | L=10 | L=15 |
|---|---|---|---|
| out8, res=40 calibration | 1.9° | 1.8° | 2.2° |
| out8, **res=60 calibration** | **11.6°** | **25.0°** | **39.0°** |

**This is not a failure, it is the mechanism confirming itself.** The
full-device FDTD reference is at **res = 40**. A res=60 calibration measures
what a res=60 time-stepper would propagate — a *different, better* simulation
than the one being compared against. The residual is therefore the res=40 →
res=60 difference in calibrated n_eff, which the aggregate puts at 0.009594:

| L (µm) | predicted 360·Δn·L/λ | observed (out8) | ratio |
|---|---|---|---|
| 5 | 13.2° | 11.6° | 0.880 |
| 10 | 26.4° | 25.0° | 0.948 |
| 15 | 39.5° | 39.0° | 0.986 |

Predicted from an independently-measured n_eff difference, with no fitting.
The agreement is as good as the original +0.0458 prediction was.

**THE RULE THIS ESTABLISHES: the calibration's resolution must MATCH the
resolution of the simulation it is meant to reproduce.** Calibration is not a
move toward absolute truth — it is a measurement of *what a particular
discretisation propagates*. Improving the calibration's own convergence while
leaving the reference at res=40 is not a refinement; it is a mismatch, and it
degrades agreement in direct proportion to how much better-converged it is.
Benefiting from a res=60 calibration requires a res=60 full-device FDTD too,
which in 3D costs roughly (60/40)⁴ ≈ 5× the reference run.

**Knob 1 is therefore CLOSED for the stated goal.** The 2.12° weighted error
cannot be reduced by calibrating harder. The res=40 calibration's residual
against its own matched res=40 reference is Δn ≈ 3.0e-4 — 32× smaller than the
0.0096 res-convergence difference — and *that* is the real remaining term.
Remaining live knobs: 2 (more calibration wavelengths), 3 (`N_MODES`), 5
(`BODY_STUB`). Knob 3 is now the largest unmeasured one.

Guard added while running this: `MMI4X4_CALPATH` selects a non-default
calibration, and when set the output npz is tagged with that calibration's
name, so a probe can never overwrite the res=40 result the published figures
are built from.

---

### 🔬 KNOB 3 QUEUED — the 16-mode interface basis leaks 0.7–0.8% of the input power (2026-09-13)

Free check first, from the cached 16-mode `M_in`, before committing to any
new solves. Top port, modal power fractions at 1.31 µm:

| mode | 1 | 2 | 3 | 4 | 5 | 6 | 7 | **8** |
|---|---|---|---|---|---|---|---|---|
| frac | 3.72% | 12.27% | 19.70% | 21.09% | 15.99% | 8.08% | 1.95% | **0.0014%** |

| mode | 9 | 10 | 11 | 12 | 13 | 14 | 15 | **16** |
|---|---|---|---|---|---|---|---|---|
| frac | 1.45% | 3.76% | 4.75% | 3.94% | 2.28% | 0.86% | 0.15% | **0.0010%** |

**Total captured power: 0.9932 / 0.9926 / 0.9919 at 1.29/1.31/1.33 µm** — so
**0.68–0.81% of the input power is not represented in the basis at all.** As a
field amplitude that is √0.007 ≈ 8.4e-2, i.e. of order **4.8°** of phase if it
re-enters coherently — the same order as the ~2° residual being chased. Knob 3
is therefore worth the cluster time, unlike knob 1.

**The tail is NOT simply converged.** Modes 8 and 16 sit at ~1e-5 while their
neighbours carry 1.5–4.8% each — a clean every-8th-mode symmetry suppression,
not a decaying tail. Reading "mode 16 ≈ 0, so 16 is enough" off the last entry
would be exactly wrong: it is a suppression node, and modes 17+ above it can
carry real power again. This is also why the test goes to **24, not 20** — 24
clears the next node.

| job | what | budget |
|---|---|---|
| **13152193** (array 0–3) | input sections at a 24-mode basis, L=15, one port per task | 10 h/task |
| **13152194** (array 1–24) | full 24-mode calibration **at res = 40** | 8 h/task |

**res = 40 is load-bearing here, not incidental.** Job 13110658 established
that a calibration's resolution must match the reference's; this knob varies
the MODE COUNT only, holding resolution at the reference's own 40. And a
24-mode calibration is not the 16-mode one plus 8 more — raising `N_MODES`
re-seeds the whole band solve, so all 24 are re-measured.

L=15 is the single probe length used: largest calibrated residual (2.43°
weighted) and longest propagation, so the most sensitive test per CPU-hour.

Plumbing added, with the production path regression-checked afterwards
(re-ran the default configuration and reproduced 2.2 / 2.3 / 5.8 / 33.7° at
L=15 exactly, writing to the untagged filename):
- `MMI4X4_CALIB_NMODES` on `mmi_4x4_90hybrid_calib_probe.py`, and `n_modes`
  now saved into each probe npz.
- `MMI4X4_NMODES` on `mmi_4x4_90hybrid_calib2_phase_check.py`, plus a
  `_PRODUCTION` flag (default calibration **and** 16 modes) that gates
  `cache_only`, the full-length grids, and the output filename. A probe at a
  different basis can therefore neither fail on a cache miss it was always
  going to have, nor overwrite the production result.
- New `mmi_4x4_90hybrid_input_presolve.py` — going 16 → 24 modes is a genuine
  re-solve, since the cache's `lookup_superset()` tier slices a stored result
  DOWN to fewer modes but never up.

---

### ❌ KNOB 3 ANSWERED, NEGATIVELY — a 24-mode basis recovers the missing power and changes the phase agreement by NOTHING (jobs 13152193/194/13162199, 2026-09-13)

All 28 tasks COMPLETED clean (4 input pre-solves at 2:43:52; 24 calibration
tasks at 2:33–2:45; final cascade 4:04:01).

**The premise was confirmed — the basis really was leaking power.** Captured
input power at L=15, top port:

| basis | 1.29 µm | 1.31 µm | 1.33 µm | leak |
|---|---|---|---|---|
| 16 modes | 0.9932 | 0.9926 | 0.9919 | 0.68–0.81% |
| **24 modes** | 0.9977 | 0.9973 | 0.9969 | **0.23–0.31%** |

Modes 17–24 carry 0.48–0.64%, so the every-8th-mode suppression argument was
right: the tail was NOT converged at 16, and going past the mode-16 node found
real power. **Roughly two thirds of the missing power was recovered.**

**And it bought nothing.** Same res=40 reference, same L=15, same resolution,
only the basis changed:

| arm | 16 modes | 24 modes |
|---|---|---|
| **out8** (79% of the weighted metric) | **2.2°** | **2.2°** |
| out7 | 2.3° | 2.4° |
| out6 | 5.8° | 6.3° |
| out5 | 33.7° | 35.9° |
| **power-weighted** | **2.425°** | **2.468°** |

The dominant arm did not move at all, and the aggregate got marginally worse
(+1.8%, within the per-wavelength scatter). **Basis truncation is RULED OUT as
the source of the ~2° residual.** Recovering two thirds of the leaked power
moved the dominant arm by 0.0°, which is only consistent with that power never
having been the mechanism — it leaves as radiation/continuum rather than
re-entering the guided cascade coherently.

**Knob 3 is CLOSED. So are knobs 1 and 4.** Of the original list only knob 2
(calibration wavelengths) and knob 5 (`BODY_STUB`) remain — and the arithmetic
below suggests neither is the real target either.

**WHERE THE RESIDUAL ACTUALLY LIVES — a new candidate, from the numbers now in
hand.** The residual is a systematic Δn_eff ≈ 3.0e-4 (it grows linearly with
body length). The calibration's own β-fit residual at res=40 is 5.58e-3. Over
the fit's own baseline (lengths 0.84–2.53 µm), a phase residual of that size
implies Δβ ≈ 5.58e-3/2.53 ≈ 2.2e-3 rad/µm, i.e.

    Δn_eff ≈ 2.2e-3 × 1.31 / 2π ≈ 4.7e-4

**the same order as the 3.0e-4 actually observed.** That points at the
**number and SPAN of the calibration lengths** — currently three, with spacing
set by the branch-safety formula `max_safe_spacing/2` — rather than at
anything in the cascade. A longer baseline determines β more precisely, and β
precision is what the residual is made of. This is a knob that was not on the
original list and is now the leading candidate; it is also cheap relative to
knobs 1 and 3, since it changes only the calibration's own geometry.

Corroborating: at res=60 the β-fit residual fell to 1.27e-4 (44×), which by the
same arithmetic would give Δn_eff ≈ 1.1e-5 — but that gain is unusable because
a res=60 calibration must be paired with a res=60 reference (job 13110658).
Improving the FIT at fixed resolution is the way to get the same benefit
without moving the reference.

**One anomaly recorded, not yet resolved.** Mode 24's calibrated offset is
**+0.0796** against a smooth +0.0457 → +0.0486 trend across modes 1–23 — a 1.7×
outlier that single-handedly takes the mode-spread from 0.00046 (16 modes) to
0.0338. This is the familiar "topmost requested band of N is unreliable"
signature this device already hit at mode 16 of 16. It was NOT excluded from
the 24-mode calibration used above. Given the result was a clean null on the
dominant arm, it is unlikely to have changed the conclusion, but any future
use of `aggregated_n24res40_n24.npz` should either drop mode 24 or re-measure
it at a larger `N_MODES` so it is no longer the top band.

**Method note worth carrying:** at 16 modes this comparison is cached
arithmetic (~19 s); at 24 it is a fresh EME mode solve and took 4:04:01. A
first attempt on a login node produced no output in ~4 minutes and was killed
rather than left to contend with other users. `MMI4X4_LENGTHS` was added so a
non-16 basis cannot request a length whose input sections were never
pre-solved — `M_in` is fetched `cache_only`, so that would fail rather than
silently start a multi-hour FDTD in the wrong place.

---

### 🔬 KNOB 4 REOPENED — it was a two-point delta, not a convergence test (job 13197102, 2026-09-14)

**The gap.** Job 13110659 measured the reference's phase shift between
`decay_tol` 1e-4 and 1e-5 and reported 0.029° power-weighted, concluding the
1.35° target is not floor-limited. That is **two points**. A single pair of
rungs bounds the step between them and cannot distinguish "the phase has
converged" from "the phase is still drifting and both rungs are early" — and
the entire "knobs 1–3 are worth spending on" decision rests on it. The
amplitude side of this question got a real convergence curve on 2026-08-25
(1e-2 → 5e-5, six points); the phase side never did.

**Can the tolerance go tighter than 1e-5? Unknown, and the assumption in force
is borrowed from a different geometry.** The practical ceiling of ~1e-5 comes
from an estimated floor of ~5e-6 at res=40 — but that floor was measured on
the **84 µm full device**, not on this L=15 probe. What is actually established
about the floor is that it is **spatial**: resolution-dependent (1.68e-4 at
res=20), invariant under a 2× timestep change (0.01%), and invariant under
`dpml` 0.5→1.5 (≤1%). A much smaller domain is therefore not bound by the
84 µm number, and the 1x2 MMI at 6.5 µm reached 1.06e-13 with **no turnaround
at all**. So tolerances below 1e-5 may well be reachable at L=15. Nobody has
looked.

**One run answers both, because DFT monitors accumulate.** Decomposing the
monitors at time *t* returns exactly what a run stopped at *t* would have
returned, so a single densely-sampled run IS the phase-vs-`decay_tol` curve —
the same mechanism that produced the amplitude table, collapsing a ladder of
runs to one job.

**The library was throwing the phase away.** `solve(snapshot_every=)` recorded
`{t, energy, energy_ratio, T}` and the complex coefficient never reached the
CSV, even though `snap` held it the whole time. `beef/solvers/fdtd.py` now also
writes `ph_out{i}_wl{j}` = `arg(S)` in degrees, **appended after** every `T`
column — both tolcurve plotters select by name on `startswith("T_out")`, and a
positional reader keeps the same leading fields, so old CSVs stay readable by
new code and new CSVs by old. Verified on a 2-D smoke solve against the
returned `SMatrix`: phase agrees to **<0.01°** on every arm and wavelength, `T`
columns unchanged. Not in the cache fingerprint (a method argument, never
stored on `self`), so a snapshotted run still shares a cache entry with a plain
one. `MMI4X4_FDTD_SNAP_EVERY` exposes it on `mmi_4x4_90hybrid_fdtd.py`, one CSV
per driven port (the writer appends, so a shared path would interleave two
runs into one unreadable file). Production path re-checked: with the variable
unset, `snapshot_every=None` and `output_dir` is byte-identical to before.

| job | what | budget | batch file |
|---|---|---|---|
| **13197102** | L=15 top-port, `decay_tol=1e-9`, snapshots every 50 t-units | 15 h | `run_mmi_4x4_90hybrid_fdtd_tolcurve.batch` |

**`decay_tol=1e-9` is set below what the run can reach ON PURPOSE**, and a
**timeout is the expected outcome, not a failure** — the stopping point is not
what is being measured. Samples append as they are taken, so a walled run keeps
every sample it wrote. If the energy instead floors and turns around (this
device's measured res=20 signature: 5.16e-3 → 51.9 → 2.14e8), **that turnaround
is the floor measurement** — read it off `energy_ratio` and discard every
sample after it. Standing directive applies: if it times out, report and wait,
do not auto-resubmit.

Output lands in `mmi_4x4_90hybrid_fdtd_L15.00_p4_tol1e-09/tolcurve_p4_port4.csv`
— distinct from the reference (13069784, `_L15.00_p4`) and from knob 4
(13110659, `_L15.00_p4_tol1e-05`). `decay_tol` is in the cache fingerprint, so
this both cannot be served to a production request and cannot be served *from*
one, which is what guarantees it time-steps rather than returning instantly.

**What this can change.** If the phase curve is flat well before 1e-5, knob 4's
conclusion stands on real evidence instead of one interval and the matter is
closed. If it is still drifting at 1e-5, the 0.029° bound is an underestimate,
the reference is less settled than assumed, and the ~2° hybrid-vs-FDTD residual
is partly the reference's own — which would reframe knobs 2 and 5 and the
calibration-length candidate alike, since all three are being scored against
that reference.

Submitted against `gts-sr108-gedc_main`: `gts-sr108-ece` is down to 0.82
available (it is what is holding task 3 of the 1x4 prewarm array in
`AssocGrpBillingMinutes`) and `gts-sr108` to 10.50.

### ✅ MATCHED res=60 comparison — 2026-09-15, jobs 13193430/31/32

The res=60 calibration (job 13110658) was unusable because no res=60 reference
existed to pair it with. Those three references now exist:

| L (µm) | job | wall | stopped at | res=40 twin stopped at |
|---|---|---|---|---|
| 5 | 13193430 | 4:11:16 | t=280.08 (33610 steps) | 280.18 |
| 10 | 13193431 | 3:44:12 | t=360.12 (43214) | 380.24 |
| 15 | 13193432 | 6:44:49 | t=540.19 (64823) | 560.35 |

All three decayed slightly EARLIER in physical time than their res=40 twins
(−0.03 / −5.3 / −3.6%), so finer resolution did not lengthen propagation here.
Reference saved as `example_output/mmi_4x4_90hybrid/fdtd_reference_res60.npz`.
`mmi_4x4_90hybrid_calib2_phase_check.py` gained an **`MMI4X4_FDTDREF`**
override (additive; the default reproduces the hardcoded res=40 values
exactly) because its existing `CALPATH` override deliberately holds the
reference fixed — right for isolating the calibration, wrong for a matched
test.

#### The rule is CONFIRMED — on the length-SCALING term

| configuration | out8 L=5 / 10 / 15 | slope |
|---|---|---|
| res60 cal + res40 ref (mismatched) | 11.6 / 25.0 / 39.0 | **+2.74 °/µm** |
| **res60 cal + res60 ref (MATCHED)** | **20.0 / 20.0 / 18.4** | **−0.16 °/µm** |
| res40 cal + res40 ref (baseline) | 1.9 / 1.8 / 2.2 | +0.03 °/µm |

Matching the resolutions flattened the slope **17×**. The β-mismatch term the
previous section predicted and measured is gone.

#### 🔴 But it exposed a LENGTH-INDEPENDENT ~19.4° offset

Power-weighted, matched res=60: **19.41°** pooled (5.39% of 360°, 21.6% of
90°) versus **2.12°** at matched res=40.

**Quantitatively explained by the taper sections still being res=40.** The
hybrid's M_in/M_out are cached at res 40 while the reference is now res 60, and
a fixed path length with a Δn error gives a constant phase:

> 360 × Δn(res40→res60) × L_taper/λ = 360 × 0.009594 × (2 × 4.0)/1.31 = **21.1°**

Observed **19.41°**, ratio 0.92, no fitting — the same predictive form the
original mismatch used.

**EACH UNMATCHED COMPONENT HAS ITS OWN SIGNATURE**, which makes this
diagnosable rather than mysterious:

| unmatched component | signature |
|---|---|
| body calibration | error **scales with length** (Δβ·L) |
| taper sections | error **constant with length** (fixed path) |

#### What this does and does not answer

**Does NOT yet answer** whether the 2.12° res=40 residual is
resolution-limited. A fully matched res=60 comparison needs res=60 M_in/M_out
as well; until then the 19.41° is dominated by a known, predicted taper
mismatch and says nothing about convergence.

**Next step, costed**: re-solve the input/output sections at res 60
(`mmi_4x4_90hybrid_input_presolve.py` at `MMI4X4_FDTD_RES=60`). Only then does
"res=60 everywhere vs res=40 everywhere" become a real convergence datum.
res=80 remains blocked behind a res=80 calibration that does not exist.

### 🟢 RESOLUTION CONVERGENCE ANSWERED — and the absolute metric was the wrong one

Prompted by the question "is the phase relative to the input or to the other
outputs?" — the weighted metric differences **absolute per-arm phase(S)**. A
90° hybrid is specified on the **quadrature relationship between outputs**, so
a phase common to all four arms is physically meaningless (it depends on
arbitrary reference-plane placement and cancels in every spec quantity).

Recomputing both configurations with each arm referenced to out8:

| configuration | absolute (wtd) | **relative / quadrature (wtd)** |
|---|---|---|
| res=40 matched (baseline) | 2.12° | **0.21°** |
| res=60 matched | 19.41° | **0.21°** |

**THE RELATIVE ERROR IS IDENTICAL AT BOTH RESOLUTIONS — 0.21°.** In the metric
the device actually operates on, the hybrid-vs-FDTD phase agreement is
**already converged at res=40**; res=60 moves it by nothing measurable.

0.21° is **0.23% of the 90° quadrature spacing** — it passes even the strict
1.5%-of-90° reading (1.35°) by 6×, the target the absolute metric reported as
NOT met (2.36%).

**The 19.41° at res=60 is pure common-mode**: the per-wavelength offsets are
19–22° and uniform across all four arms at every length. The earlier
attribution to res=40 tapers against a res=60 reference was CORRECT about
origin (predicted 21.1°, observed 19.41°, ratio 0.92) and WRONG about
importance — it is a global phase.

**Consequences:**
- The proposed res=60 input-section presolve is **withdrawn as unnecessary**;
  it would chase a term that cannot affect device function.
- **res=80 is likewise unnecessary** for this metric, which removes the
  blocker that no res=80 calibration exists.
- The length-scaling term still matters and calibration is still mandatory —
  that term is NOT common-mode (it is a genuine Δβ·L per-arm error, and at the
  159.30 µm production length the uncalibrated Δn≈0.0458 is >4500° of drift).

⚠ **METHODOLOGICAL: `mmi_4x4_90hybrid_weighted_phase_metric.py` should report
the relative-phase reading alongside the absolute one.** Its docstring is
already careful to power-weight so noise arms cannot dominate — the same care
is needed for common-mode phase, which at res=60 inflated the score **90×**
(19.41° vs 0.21°) without any change in device behaviour. Compare
[[criterion-text-vs-code-phase-mean-max]]: this repo has already been bitten
once by a phase criterion whose stated and computed quantities differed.
