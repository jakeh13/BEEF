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

**CURRENT STATUS (2026-08-09) — read this first, history follows below.**
Steps 1-7 done. Step 6 (all 3 example scripts) done and RE-RUN at a
properly-converged body mode basis -- see "Body/taper mode-count
convergence investigation" (before Step 8) for the full trail. **The
numbers that matter right now:**

- `N_MODES_TAPER=12` (mmi_4x4_90hybrid_eme.py), `N_MODES_BODY=16`
  (body_sweep.py/hybrid.py) -- both converged, both measured, not
  guessed. `num_cells=30` for the EME-only script's taper (job 11736983).
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
- **2026-08-14, mm_width Stage-1 scoping in progress, job 11840371** —
  cheap L_π-vs-width + guided-mode-count probes (9 candidate widths,
  3.5–10.0 µm, `mmi_4x4_90hybrid_mm_width_scope.py`) to pick 1-2 promising
  widths before committing to the expensive Stage-2 full treatment (own
  body solve + baseline 4-solve per width, no cache reuse possible).
  Still `RUNNING` within its 8h budget as of this update; results pending.
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
(port 3)). Results pending.

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
