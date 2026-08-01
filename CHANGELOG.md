# Changelog

Dated record of library changes. Companion to `notes.md` (which holds the
standing recommendation lists and the reasoning behind each item) and
`validation_markdown/` (per-device validation status). Where an entry has a
fuller write-up, it points there rather than repeating it.

---

## 2026-08-01

Two library bugs, both found by disbelieving something already written down as
settled, plus the material-dispersion deployment.

### Material dispersion is ON by default

`MaterialStack(dispersive=True)` is now the default, so `SOI_160`/`SOI_220` —
and therefore every example and both solvers — track the band instead of
freezing one index at band centre. The gate that blocked this passed 2026-07-31
on all four stack/band combinations (worst Δn_eff −0.00141 vs a 0.002
tolerance).

**This closes the engine asymmetry, which was the actual gap.** A broadband
`FDTDSolver` froze its index while an `EMESolver` re-solved per wavelength, so
every hybrid chain contracted modes solved in two different materials:

| | worst n_Si gap, O-band |
|---|---|
| before (frozen FDTD + exact EME) | **0.007862** — over the 2e-3 guard |
| after (both dispersive) | **0.001349** — under it |

`HybridSolver`'s material-consistency guard was warn-only *because the shipped
default failed it*; `material_consistency_action="raise"` is now usable on a
default-constructed chain. The residual 0.001349 is fit error, not disagreement
— the FDTD runs the re-parameterised single-Lorentzian fit, the EME stays exact
Sellmeier.

**A cache-fingerprint bug had to be fixed first, or this would have poisoned the
cache silently.** `resolved_material_config()` keyed off
`material_stack.dispersive` while all seven call sites passed the stack, so a
per-solver override never reached the key. Measured:
`FDTDSolver(material_stack=SOI_220, dispersive=True)` built a medium with a real
Lorentzian susceptibility against a flat ε = 12.096 and produced a
**byte-identical cache key** — a dispersive run colliding with frozen-physics
entries. Same defect that function exists to kill, one level up; it bit in both
directions (an explicit `EMESolver(dispersive=False)` built frozen materials
under a key that said dispersive). The function now takes the caller's resolved
value; `EMESolver` gained `_resolved_dispersive()` so the fingerprint and the
extrusion rules cannot disagree about what `None` means. `notes.md` Solvers
item 21.

**Cost: every cache entry misses and every golden is stale again** (third
invalidation — Sellmeier n_Si, then λ₀, now this). The EME half was invalidated
by the correctness fix alone, independent of the default flip: those entries
were dispersive all along while keyed as frozen.

**Not deployed to the ring.** `ring_resonator_hybrid{,_cascaded}.py` hand-roll
raw meep and pass a scalar `wl_center` to `si_medium()`, which correctly stays
frozen — leaving them the only devices on frozen media while their meow half is
already dispersive, i.e. carrying the exact asymmetry this flip closed
everywhere else. Their hand-rolled `_sim_params()` also records no material
fingerprint at all, so neither the Sellmeier change nor this one would move a
ring cache key. Itemised as M1–M4 in `ring_resonator_validation.md`.

Tests **508 → 514 passed**, 15 skipped. Nine encoded the old default and were
rewritten rather than relaxed; the clearest is
`test_default_configuration_is_flagged` → `test_the_default_configuration_is_now_quiet`,
whose inversion is the deployment's acceptance criterion.

### `beta_from_phases()`: wrap counts were chosen per length, not consistently

`beef/solvers/fdtd.py`. `turns = round((seed_beta*lengths - phases) / 2π)` ran
**independently at every length against one seed slope**, so a seed error ε
drifted by `ε·L/(2π)` and eventually straddled a `.5` boundary, giving one
length a turn count inconsistent with its neighbours. ε is the FDTD-vs-MPB gap
that `calibrate_body_neff` exists to measure (~2%), so **the unwrapping was
broken by exactly the quantity it was built to measure.** It killed job
11589770 (mmi_1x2 calibrated hybrid) with a 2.101 rad residual against its own
0.1 guard.

That failure had been recorded on 2026-07-31 as a **π sign flip in `T_ii`**.
It is not one. The residual pattern `δ·(+⅓, −⅔, +⅓)` is **degenerate** — for
three equally spaced points, "middle point off by −π" and "points 2 and 3 each
off by −2π" give it identically — so the arithmetic that looked like
confirmation fits both readings. Physics separates them: turn-count gives
n_eff **2.783** (+1.96% over MPB, in family with the coupler's +2.233%/+2.361%),
sign-flip gives **2.101** (−23%, wrong sign *and* size). A seed-free invariant
settles it without any fit at all: the second difference
`(φ₃ − 2φ₂ + φ₁)/2π` must be an exact integer for a linear ramp, and it is —
within **0.00308 in all 88** (wavelength, mode) cells, 31 needing a ±1 turn.
Nothing needs a sign flip.

- **Added `beta_fit_from_phases()`** returning a `BetaFit` namedtuple
  (`beta, resid, margin, seed_excess, turns, seed_turns, turns_shifted,
  branch_ambiguity, turn_ambiguity`). It decides exactly **one** integer
  against the seed — the branch — and rounds every later turn against a slope
  **refit from the shorter lengths**, so drift cannot accumulate. O(n).
  `beta_from_phases()` stays as a 3-tuple facade, so no existing unpack site
  changed.
- **Rejected** a ±1 brute-force turn search: it is residual-degenerate between
  neighbouring branches (differing in the 9th significant figure) and
  min-residual picks β = 7.072, 49% wrong, with a *clean* residual — trading a
  loud failure for a silent one.
- **`_BETA_FIT_RESIDUAL_MAX` (0.1, unchanged) now means one thing** — the phase
  is not linear in L. Its comment's "clean fits give 1.4e-5" was the 2-mode
  configuration only; the real 6-mode coupler cache reaches **2.2e-2**, 4.5×
  under, so the threshold was deliberately *not* tightened.
- **New `_BETA_BRANCH_MARGIN_RAISE = 0.9`** makes `section_mode_beta`'s
  long-standing docstring promise of a margin `ValueError` true for the first
  time. New `_BETA_SEED_EXCESS_WARN = 0.10` (warn only — the excess is the
  measurement, not an error) and `_BETA_TURN_POLISH_MAX = 8` (raises rather
  than returning an unconverged turn vector).
- **`section_mode_beta()` keeps full `(n_wl, n_modes)` grids** instead of
  running maxima, so raises and logs **name the offending cell**. The old code
  could report "2.101 rad" without saying where — and where (modes 0–1, the
  *shortest* wavelengths, i.e. largest β and tightest branch spacing) was the
  most informative fact about the failure.

**Backward compatibility, measured on cached production data.** At the
coupler's own seed the fix changes **0 of 66 cells**, so its L₅₀ 0.11% result
is untouched. On mmi_1x2 it changes **exactly the 4 of 88 cells** that blew the
guard. Seed-free sweep over ±3% seed error: mmi_1x2 4025 trials, old wrong
**13.7%**, new **0%**; coupler 3650 trials, old **5.4%**, new **0%** — the
coupler's working calibration was luck, not margin.

**Tests: 29 → 42** in `tests/test_beta_calibration.py`. Exactly one existing
test broke — `test_a_boundary_seed_shows_up_as_a_broken_fit`, which encoded the
*old* contract (a boundary seed no longer shows up as a blown residual, and
never should have been the tell) — rewritten as
`test_a_boundary_seed_is_undecidable_and_the_margin_is_what_says_so`. Widened
`_wire_fakes` to multi-mode/multi-wavelength with `offset`/`seed_scale`, which
buys a meep-free end-to-end reproduction of the production failure. Added
`TestTheMmi1x2TurnCountRegression`, pinned to the real cached numbers.
Fast suite **508 passed, 15 skipped**.

`examples/verification/python/diagnose_mmi1x2_calibration_phase_fit.py`
rewritten to carry the whole proof: old-vs-new grids side by side, the
seed-free second-difference check, the robustness sweep, and a printed seed
literal so re-validation needs no MPB call. Full write-up in `notes.md`
Solvers item 20 and `validation_markdown/mmi_1x2_validation.md` gap 5.

### The 2×2 MMI has golden runners for the first time

`run_mmi_2x2_{eme,fdtd,hybrid}` added to `tests/test_golden_examples.py` and
registered in `TARGETS` — 4 golden targets → 7. That device had no golden of any
kind and, more to the point, **no runner to generate one from**, so
`mmi_2x2_validation.md` gap 1 was never merely "not generated". The runner half
is code rather than compute and did not have to wait on the outstanding cold
regeneration pass; doing it now means that pass covers all three devices in one
round.

Each replica follows its own example's conventions — C-band 1.530–1.565, top
input excited (`input_idx=1`), bar = `output_idx=1`, cross = `output_idx=0`, and
`half_domain_aperture` on both end sections (this device has two ports 1.5 µm
apart at *each* end, unlike the 1×2). Two deliberate departures, both commented
in place: `n_modes = 6` against production's 14, chosen conservatively because
`_iter_distinct_bands()` now *raises* on an unsatisfiable mode request rather
than padding; and the hybrid replica simulates every section instead of reusing
the example's mirror-symmetry transpose, since a golden pins
`HybridSolver.solve()` on its defaults.

### Goldens now record what they were generated against

`scripts/generate_golden_values.py` writes a `provenance` block: `beef_git`,
`generated_utc`, `dispersive_default`, and the Sellmeier si/sio2 indices at the
run's own band edges and centre.

Until now a golden held `wavelengths` + `transmission` and nothing else, so
**staleness could not be read off the file** — only inferred from the validation
docs. These four fixtures have been invalidated four separate times (Sellmeier
n_Si shift; `N_MODES` 12→8 on the 1×2; the λ₀ re-parameterisation; today's
dispersion flip) and on none of those occasions did anything in the file change
to say so. Same lesson `README.md` check 4 already recorded for the full-device
FDTD reference, never applied to the `.json` fixtures.

Records the resolved material model rather than a version tag, so a change
underneath moves it on its own — the reasoning `resolved_material_config()`
already uses, and the same three probe points. Descriptive only:
`test_golden_examples.py` still compares `transmission` alone, and one of the
new tests pins that it stays that way (the git sha changes every commit, so
asserting on provenance would fail every golden immediately).

5 new tests, `tests/test_golden_provenance.py`, meep-free. Needed no compute and
did not wait on the outstanding cold regeneration pass; the existing four files
stay bare until they are next regenerated.

### Examples

- `verify_mmi_1x2_engine_agreement.py` reshaped to the array+merge shape the
  2×2's sibling proved (`SLURM_ARRAY_TASK_ID` picks one section length; the
  merge step re-reads all four under `cache_only=True`). The script already
  existed, fully ported to this device, and had simply never been submitted —
  this is `mmi_1x2_validation.md` gap 4 / README check 1, the device's last
  unrun check. Batch pair `run_mmi_1x2_engine_agreement{,_merge}.batch`;
  deliberately launched with bare `python`, not `srun`, since the script has no
  rank guards.
- `run_audit_interface_modes.batch` honours an `AUDIT_ARGS` override instead of
  hardcoding `--convergence` across all devices, so a single device can be
  re-audited without editing the file. Added because the dispersion flip makes
  every recorded identity number frozen-era and worth refreshing per device —
  and because an all-device run is not currently safe: the modal multiplexer
  requests 16 modes where ~14 are reachable at res 40 and now raises on the
  guided-spectrum stop.

  **Correction, later the same day.** The first use of that override (the 1×2)
  came back **byte-identical**, so "frozen-era" overstates it for this class of
  number: these are per-wavelength mode solves, and a dispersive medium only
  differs from a frozen one *across a band*. Only genuinely broadband
  quantities move. The override is still worth having — that was measured, not
  assumed — but expect re-audits to confirm rather than overturn.

---

## 2026-07-31

Started as infrastructure (a billing-account fix and a re-queue), then became a
material-model session: dispersion is now implemented, gated on real solves
across every stack/band the examples use, and wired to per-solver flags. Sections
below are in the order they were written, so later ones supersede earlier ones
where they overlap — notably λ₀, which moved twice (10.0 → 2.0 → derived per
band) as each measurement replaced the previous reasoning. Fast suite reaches
**477 passed, 1 skipped** by the `section_mode_neff` dispersion fix below (468
before its 9 new tests), then **489 passed, 1 skipped** by day's end (the
`interface_modes_convergence()` and `CompoundGeometry` clamp work adds 12
more).

The re-queued jobs then started landing, and the last sections are mostly their
results rather than planned code changes: the 1×2 MMI's first hybrid-vs-FDTD
number (9.11% low, systematic), **item 13 closed** (the EME duplicate guard's
over-request is exactly inert and costs 3.2×), and the meow guard caught firing
on the *production* path for the first time.

Then a job **hung** and a job **measured nothing**, and between them they
produced the day's largest findings — none of which were what those jobs were
sent to measure:

- the production duplicate guard can hang a run indefinitely on an eigensolve
  that has no solution (new **item 14**), and short of that it pads a basis with
  unguided modes;
- `FDTDSolver.section_mode_neff()` ignored `dispersive=` outright, so a material
  model work item had been recorded as *started* when it had never run;
- **items 11 and 12 both resolved** once the hung job was reshaped — the guard's
  diagnostic grid is clear at res 40 too, and k-seed laddering recovers a guided
  channel today's fixed seed loses silently, at no extra cost;
- **2×2 MMI gap 6 closed** — "the cleanest device audited" holds up under both
  duplicate guards, neither of which fired.

Two jobs were cancelled mid-flight and resubmitted after fixes, both called from
reading job output rather than waiting on exit status. In both cases an earlier
diagnosis recorded in these docs had to be corrected: a wall-clock budget miss
that was really a non-terminating eigensolve, and a stale λ₀ that was really a
flag never reaching the physics. A third correction followed from the first — the
mux's reachable guided-mode count is 14, not the 13 inferred from the hung job's
own fixed-seed trace.

The 2×2 MMI's own gap-5 chain then finished: 0.82% worst |ΔT| against the 1×2's
9.11%, but with a sign that flips across the band instead of the 1×2's uniform
under-prediction — different mechanism, calibration not obviously indicated the
same way. A calibration run for the 1×2 (where the signature *does* match the
coupler's) went out. And writing the mux's first full-device FDTD reference hit
a real, unfinished library gap — `CompoundGeometry` cannot back
`FDTDSolver.solve()`'s single-domain path — with one genuine sub-bug found and
fixed along the way and the larger design question correctly left alone.

Item 6 (layer 2, completeness) got its reusable tool, the same way item 7
(layer 3) got one on 2026-07-30: generic function built and tested, an actual
device-level exercise run explicitly deferred rather than rushed. Then, at the
user's direction, all five directional-coupler sweep drivers and their batch
scripts were archived pending new, more general sweeping functionality. Last,
the mux's own hybrid exercise run (job 11585847) timed out after its full 5h
budget — a plain resubmit once the disk cache was confirmed to have carried
sections 0-3 forward, not a reshape, since the batch script's cost estimate
wasn't wrong, the device was just this expensive.

### Found — job 11585847 timed out after its full 5h; disk cache carried sections 0-3 forward, resubmitted as 11590663

Gap 6's mux hybrid exercise run (`modal_multiplexer_hybrid.py`, 4 nodes × 24
ranks, `-t 300`) hit `State=TIMEOUT` at 13:42, `ExitCode=0:0` — SLURM killed
it, it did not crash. Reading `example_sim_cache/` mtimes against the script's
fixed section order (0 input FDTD, 1 EME MMI1 body, 2 FDTD transition, 3 EME
MMI2 body — a 2h17m `solve_section_length_sweep` call, by far the single most
expensive piece — 4 output FDTD, once per port) shows sections 0-3 fully
cached by 13:40:31, 1m41s before the kill, with the log tail mid-way through
building the first section-4 (output) FDTD grid. **Net: zero of the 3 per-port
section-4 solves completed**, but the entire shared prefix is now on disk.

Resubmitted unchanged as **11590663** — nothing touched since the original
08:41 submission changes `modal_multiplexer_hybrid.py`'s or the solvers'
cache fingerprints, so a clean hit on sections 0-3 is expected, leaving only
the 3 section-4 solves as real new work (comparable single-section FDTD
solves earlier in the same log took 5-33 min, well inside the unchanged 5h
budget). `running_jobs.md` and `modal_multiplexer_validation.md` gap 6
corrected — both had the job recorded as still `RUNNING`, and one line in
gap 3 read "has run fine (job 11585847) the whole time," which was true of
`solve_section()` never hitting the `CompoundGeometry.polygons()` bug but
overstated the job's own outcome.

### Measured — mmi_2x2_validation.md gap 2 closed at res 100 (job 11590785)

Δβ_01 (beat length between the MMI body's two lowest guided modes — the
quantity `mm_length=28.0` was actually chosen against, not either mode's own
n_eff) converges to 0.05% relative step by res 120, closing at res 100 with
generous headroom. Corrected an imprecise comparison in this session's own
earlier partial reading of this job: Δβ_01's raw step size (~0.16 scale) was
compared directly against absolute n_eff's raw step (~2.8 scale) and called
"an order of magnitude tighter" — a mixed-units comparison. Read relative to
each quantity's own value, the two converge at a similar rate; what actually
settles gap 2 is that Δβ_01 itself is small enough to trust by res 100, not
a dramatic difference in convergence rate. `mmi_2x2_validation.md`'s own
`## Order` updated — the goldens (gap 1) are now blocked only on the
repo-wide material-model regeneration, not this device's own resolution
question.

### Fixed — `FDTDSolver._iter_distinct_bands()`: k-seed laddering (item 12) and the guided-edge stop (item 14)

Both items this session's earlier diagnostic jobs (11588309, and the hang that
produced it) opened, closed together in the one method they share. The seed
now ladders — attempt *k* starts from attempt *(k−1)*'s own raw converged
n_eff instead of a fixed guess every time — recovering a guided mode the
fixed seed silently drops, at no extra cost (policy B, measured
byte-identical to the fancier step/accepted-only variants, so shipped as the
simplest one). And the walk now stops the instant a raw solve's n_eff drops
to or below the cladding index, raising immediately if fewer than `n_modes`
were found by then, instead of either padding the basis with unguided modes
or continuing toward a `band_num` where MPB's k-search can limit-cycle
forever (measured: 161 solves on one such band before SLURM killed the job
that found this). Landed together deliberately — laddering alone would still
be exposed to the hang, and the stop condition is what makes it safe to ship.

Both fixes were diagnosed and verified offline first
(`diagnose_kpoint_ladder_v2.py`, job 11588309) before touching the production
method, per this session's own standing caution about changing this
particular code path under-verified. 6 new meep-free tests
(`tests/test_duplicate_band_guard.py::TestSeedLaddering`,
`TestGuidedEdgeStop`), mutation-checked — all 11 tests in the file fail
against the pre-fix signature. Fast suite 489 → 495.

**Known, deliberate consequence:** any device whose shipped `N_MODES` exceeds
what a cross-section actually supports will now raise instead of silently
padding. The modal multiplexer is exactly this case (16 requested, 14
reachable at res 40) — its currently-running hybrid job (11590663) started
before this change and is unaffected, but the *next* resubmission of
`modal_multiplexer_hybrid.py` will raise until `N_MODES`/`RESOLUTION` are
decided jointly (`modal_multiplexer_validation.md` gap 2). Not fixed as a
side effect of this change — that doc is explicit the two need a real joint
decision, not a quick unblock to whatever today's ceiling happens to be.

### Measured — material model work item 4 closed for the mode-solve half (job 11588997)

The re-run of `verify_item_8_dispersive_convergence.py` after its `dispersive=`
plumbing fix (below) completed in 1:56:53. At the coupler body's band centre,
dispersion changes nothing: identical convergence resolution to the
non-dispersive column at every threshold tested. At the band edge — the column
that actually matters, since centre is near-inert by construction — the 1e-3
threshold moves from res 60 to res 80 under dispersion, but item 11 already
pins production's EME resolution at 80 for an unrelated reason (Δβ on this
same cross-section), so the dispersive requirement lands exactly on what
already ships. The one real gap: at the tightest threshold (3e-5) the edge
never converges under dispersion within the tested range (20-120), though the
dispersive-vs-non-dispersive difference itself oscillates at noise-floor
levels with no resolution trend. **Verdict: the proposed 11-array-task,
up-to-24h full-device sweep is not warranted** at the tolerance item 11
actually used. Full table in `material_model_validation.md` work item 4;
`running_jobs.md` updated.

### Added — `verify_mmi_2x2_gap2_dbeta.py`, testing the right quantity for mmi_2x2's own resolution question (job 11590785, in flight)

`mmi_2x2_validation.md` gap 2 measured the MMI body's fundamental n_eff still
moving at res 100 but explicitly declined to conclude anything from it — the
doc names item 11's coupler lesson (sweep Δβ, not T_cross) as the reason.
Applied here: a 2×2 MMI's self-imaging length depends on the beat length
between its two lowest guided modes, Δβ_01 = 2π/λ·(n_eff[0]-n_eff[1]), not
either mode's own n_eff — same logic, different device, different mode pair.
New script sweeps res 40/60/80/100/120 measuring Δβ_01 directly via
`EMESolver.compute_modes(num_modes=2)` (cheap by construction — MPB solves
the lowest N bands in one call, so 2 modes instead of the device's real 14
costs little extra). A login-node smoke check at res 40 alone measured 206 s
(`OMP_NUM_THREADS=1`) — slower than expected for 2 modes — so the full sweep
runs as a proper batch job instead
(`batch/verification/run_mmi_2x2_gap2_dbeta.batch`, 1 node/8 cpus/180 min).
Chosen over an mmi_2x2 golden test (notes.md High-level #1) at the user's
direction: `mmi_2x2_validation.md`'s own `## Order` section sequences golden
tests after gap 2 closes, not just after gaps 5/6 — a sequencing this session
had initially missed.

### Added — `interface_modes_convergence()`, the reusable layer-2 completeness check (item 6)

`beef/interface_modes_convergence.py`, direct sibling of
`beef/cut_position_invariance.py` (item 7): takes one already-solved `SMatrix`
per `interface_modes` count plus the counts, reports per-`(port, in_mode,
out_mode)` wavelength-mean transmission and the **step** between consecutive
counts. Deliberately not a fitted slope like layer 3's tool — cut-position
invariance's correct answer is a horizontal line (an unmoving S-matrix, so a
near-zero slope is the criterion); a mode-count sweep is expected to move
while the basis is incomplete and plateau once it is not, with no reason to
be linear in between, so a step-based read (matching the convention already
used by `mode_spectrum_convergence.py` and `verify_item_8_dispersive_
convergence.py`) is the honest fit. `last_step` — has it stopped moving at
the *highest* count tried — is the headline metric; `max_step` is
informational. Raises on mismatched key sets, identical reasoning to layer
3's tool: SMatrices claiming to be the same device at different mode counts
must carry the same ports. 7 new tests
(`tests/test_interface_modes_convergence.py`), meep/meow-free. Exported from
`beef/__init__.py` alongside `cut_position_invariance`. Fast suite 489
passed, 1 skipped.

**Not yet exercised on a real device, and now explicitly deferred rather than
merely unstarted** — see the archive entry immediately below. Same split item
7 landed in: generic comparison done, wiring an actual sweep onto it needs a
real meep/meow campaign. `directional_coupler_hybrid_modes_sweep.py` already
swept `N_MODES = 2..12` on the coupler and would have been the closest
existing candidate (the same relationship its `_stub_sweep.py` sibling had to
item 7 before migration), but it saves raw `.npz` arrays rather than
`SMatrix` objects and grids over two axes (mode count and coupling length)
where this tool wants only mode count varied at one fixed configuration — not
a drop-in fit regardless. Whichever device is picked when this is taken up
must already be resolution-converged (item 6's own precondition): the coupler
(converged at res 80) is the safer first target than the mux, whose own gap 2
is still open and whose res-40 basis is already known incomplete (item 12).

### Archived — all five directional-coupler sweep drivers, pending new sweeping functionality

User call, not a code-quality finding: `examples/directional_coupler/python/
directional_coupler_{fdtd_sweep,fdtd_resolution_sweep,hybrid_modes_sweep,
hybrid_resolution_sweep,hybrid_stub_sweep}.py` moved to
`examples/directional_coupler/archive/`, and their six batch launchers
(`batch/sweep/*.batch`, including `run_dc_stub_sweep_calibrated.batch`) moved
to `batch/archive/` — `batch/sweep/` is now empty and removed. New, more
general sweeping functionality is planned once every device is verified;
these five hand-rolled, per-device drivers (each duplicating the same
harness — `_tag()`, `atomic_replace_write()`, an array-index CLI, per-value
`data/*.npz` + `plots/*.png`, re-aggregation — see `notes.md` Misc #5, which
named exactly this duplication before today) are retired rather than carried
forward into that redesign.

**Not a silent removal.** `tests/test_eme_convergence.py`'s
`TestExampleScriptsPinTheConvergedResolution` reads
`directional_coupler_hybrid_stub_sweep.py` and `..._hybrid_modes_sweep.py`
directly by path to guard the item-11 shared-resolution-knob regression (a
7% error that went unnoticed for months the first time). Updated to check
`examples/directional_coupler/archive/` as a fallback location rather than
dropping the two archived scripts from the parametrized list — an archived
script can still bit-rot on exactly this bug, and it might be revived, so the
guard stays. All 9 tests still pass. `directional_coupler_hybrid.py` itself
did not move and needed no change.

`examples/Makefile` lost the five sweep targets, their `clean_*` counterparts,
`SWEEP_INDEX`, and their names from `.PHONY` — dead references to files that
no longer live at those paths. `make directional_coupler_hybrid` and every
other live target still resolve correctly (checked with `make -n`); requesting
an archived target now fails with Make's own "No rule to make target" rather
than a Python `FileNotFoundError` three layers down.

Fast suite unaffected by the move itself (still 489 passed, 1 skipped) —
the `test_eme_convergence.py` fix is what keeps it that way, not luck.

### Measured — 2×2 MMI gap 5: 0.82% worst |ΔT|, a different signature from the 1×2 (job 11585849)

Coupler body reasoning aside, the two MMIs now both have a hybrid-vs-full-device
number, and they disagree in *kind*, not just magnitude:

| device | worst \|ΔT\| | sign |
|---|---|---|
| 1×2 (11585850) | 0.0359 (9.11%) | uniformly negative |
| 2×2 (11585849) | 0.0043 (0.82%) | positive at blue edge, crosses zero ~1.544 µm, negative at red edge |

The 1×2's number is `resolved/directional_coupler_validation.md` items 6/14's
Δβ-discretization signature — a single roughly-constant bias, which
`calibrate_body_neff=True` corrects and took the coupler from ~5% to 0.11%. The
2×2's sign flip is not that shape: a slope that changes sign across the band
looks like a resolution- or dispersion-sensitive phase term rather than one
constant offset, so the same fix is not obviously indicated here — worth trying
regardless, but expect a smaller or qualitatively different effect. Also
unresolved: this measurement is at the shipped unconverged `RESOLUTION = 40`
(gap 2), so part of the sign flip could be a mesh artifact rather than physics,
and the two cannot be separated until gap 2 is measured on this device.

### Added — testing `calibrate_body_neff` on the 1×2 MMI for the first time (jobs 11589770/11589779, in flight)

`mmi_1x2_hybrid.py` gained a `MMI1X2_CALIBRATE_BODY` env flag
(`HybridSolver(calibrate_body_neff=CALIBRATE_BODY)`, output redirected to
`mmi_1x2_hybrid_calibrated/` so the uncalibrated gap-5 baseline survives — same
pattern as `directional_coupler_hybrid_stub_sweep.py`'s `DC_CALIBRATE_BODY`).
`verify_mmi_1x2_hybrid_vs_fdtd.py` gained the matching
`MMI1X2_HYBRID_CALIBRATED` override, pointed at the calibrated pickle and a
separate output dir; the shared `mmi_1x2_fdtd.py` reference is **not**
redirected, same reasoning as the coupler's stub sweep not redirecting its own
shared reference. `calibrate_body_neff` has never been exercised on this device
before, despite fixing the exact discrepancy shape (`~5% → 0.11%` on the
coupler) that gap 5 just measured here (`9.11%`, uniformly negative).

### Found — `CompoundGeometry` cannot back `FDTDSolver.solve()`'s single-domain path (writing modal_multiplexer_fdtd.py, gap 3)

Wrote `modal_multiplexer_fdtd.py` — the mux's first attempted full-device FDTD
reference, same `CompoundGeometry` the hybrid script already runs successfully,
one single-domain `solve()` call per input mode (TE0, TE1) instead of five
cascaded sections. It is the first script in this repo to call `.solve()`
rather than `.solve_section()` (the only thing `HybridSolver` ever calls) with
a `CompoundGeometry` as the top-level device — every existing full-device
script passes a plain `MMI`. That untested combination broke twice, on two
independent bugs.

**Bug 1, FIXED.** `solve()` queries `geometry.source_aperture_at(x)` at the
solver's external PML-clearance stub, past the device's own physical extent.
Every leaf geometry already tolerates this by construction —
`MMI.width_at()`/`source_aperture_at()` clamp to the terminal value for
arbitrary overshoot, which is exactly why the plain-`MMI` full-device scripts
never hit it. `CompoundGeometry._locate()` had no segment "out there" to
delegate to and raised instead — confirmed as load-bearing, not cosmetic:
`solve()`'s own source-position query lands in precisely that region, so the
*real simulation* would have crashed here too, not only the geometry-preview
plot that caught it first (attempt 1, 22 s).

`beef/geometry/compound.py`'s `_locate(x, tol, clamp=False)` gained an opt-in
`clamp=True` that delegates out-of-bounds x to the nearest end segment instead
of raising, wired into `source_aperture_at()` only. `width_at()`/`center_at()`/
`section_prisms()`/etc. are unaffected — still raise by default, since an
out-of-bounds cut there is a genuine caller error building physical geometry,
not a solver's legitimate external-stub query. 5 new tests,
`tests/test_geometry.py::TestCompoundGeometrySourceApertureClamps`,
mutation-checked: three fail against the pre-fix code (proving they exercise
the regression), the existing `test_x_outside_bounds_raises` still passes
either way (proving the opt-in changed nothing for any existing caller). Fast
suite 482 passed, 1 skipped.

**Bug 2, NOT FIXED, and the real blocker.** Clearing bug 1 exposed a second,
unconditional crash one layer in: `solve() → to_meep_prisms() →
geometry.polygons(stub_length=self.waveguide_length + self.dpml)`, and
`CompoundGeometry.polygons(stub_length>0)` raises `NotImplementedError` by
design — its own docstring: "no generic way to know per-port stub width for an
arbitrary segment type" (attempt 3, 13 s, after attempt 2's 15 s wrongly
attributed the same error to caching alone — `save_dir=None` only clears a
*different*, cache-fingerprinting call to the same method). This is real,
unfinished library work: `to_gds_eme()`/`section_gds()` already solve the
equivalent stub-drawing problem via `_end_ports()`, but generalising that to
plain `polygons()` needs an actual design decision for an ambiguous end (a
multi-port `WaveguideRun`'s stub inherits which port's width?), which the
existing code explicitly declined to guess at.

**Left alone**, on the same caution basis as item 14 and item 12's guard fix
this session: this is production-adjacent geometry logic (`polygons()`/
`to_meep_prisms()` back every FDTD example in the repo) and deserves a real
design pass plus device-level verification, not a mechanical patch under time
pressure. `modal_multiplexer_fdtd.py` and `batch/run_modal_multiplexer_fdtd
.batch` are both marked BLOCKED at the top of the file — do not run either
until `notes.md` Geometry item 8 is actually closed. Three attempts total,
~15 s and negligible node-seconds each; no wasted budget, no partial output
(only `geometry.png` was ever written — every crash happens before
`smat.save()`).

### Fixed — batch scripts billed an account that was out of billing minutes

The five jobs submitted the evening of 2026-07-30 (`mmi_1x2`/`mmi_2x2` gap-5
data + comparisons, modal-multiplexer gap-6 exercise run) never started.
`scontrol show job` gave the reason: the three parent jobs carried
`#SBATCH --account=gts-sr108-ece`, and that association is out of allocated
billing minutes — `Reason=AssocGrpBillingMinutes`, indefinitely pending.
`gts-sr108-gedc_main` had headroom (`sshare` `EffectvUsage` 0.045 vs. ece's
0.40).

Switched **every** batch script hardcoding `gts-sr108-ece` to
`gts-sr108-gedc_main` — seven files, not just the three behind the stuck jobs,
so the next submission of any of them doesn't hit the same wall:
`batch/run_mmi_1x2_examples.batch`, `batch/run_mmi_2x2_examples.batch`,
`batch/run_modal_multiplexer_examples.batch`,
`batch/run_golden_directional_coupler_hybrid.batch`, `batch/run_examples.batch`,
`batch/sweep/run_directional_coupler_fdtd_sweep.batch`,
`batch/sweep/run_directional_coupler_hybrid_modes_sweep.batch`.
(`batch/run_ring_resonator_examples.batch` bills plain `gts-sr108` and was left
alone.) Cancelled and resubmitted all five jobs; the two dependent comparison
jobs needed fresh `--dependency=afterany:<new_id>` values, since that flag was
passed on the command line, not stored in the script.

Job IDs moved 11577985→11585844, 11577986→11585849, 11578018→11585845,
11578019→11585850, 11578026→11585847. The 2×2 MMI's first array task started
within minutes of the re-queue — the first of these five to run at all.

**Not addressed:** SLURM takes one account per job, so no batch script can
express a multi-account fallback on its own. Rotating a stuck job across
accounts would need `scontrol update JobId=<id> Account=<next>` driven by a
watcher polling for `Reason=AssocGrpBillingMinutes`; noted as an option, not
built. A submit-time "pick the account with headroom" check is the weaker
alternative — `sacctmgr show assoc ... format=GrpTRESMins` returns empty for
all three accounts here, so remaining balance isn't reliably queryable before
submitting.

### Added — one batch job for items 11 and 12 (job 11586046)

`batch/verification/run_diagnose_items_11_12.batch` runs the two follow-ups the
07-30 login-node smoke runs left open, in one job because both solve the same
cross-section (modal multiplexer interface 0, res 40) and share setup cost.
Serial, eigensolve-only, 240 min budget.

Item 12 first, deliberately — it is the root-cause item and the more expensive
of the two, so a wall-clock overrun costs the cheaper result rather than the
more valuable one. It runs `diagnose_kpoint_ladder_v2.py`, written 07-30 but
never executed: a **real** guarded accept/skip/retry search (v1 was a flat
20-solve sweep with no retry) under four seeding policies — fixed, naive
ladder, step ladder, accepted-band ladder. The last two are exactly what item
12 lists as "not yet tried".

Item 11 second, at **res 40** — the one data point it names as missing, its
clear result so far being res-20 only. `diagnose_dup_band_grid.py`'s `RES` is
now read from `BEEF_DIAG_RES` (default `20.0`), so the res-20 result on record
stays the script's default and the re-check is a launch-time override rather
than an edit that would strand the earlier number.

### Added — item 13's first measurement, three arms not two (job 11586098)

`examples/verification/python/diagnose_guard_over_request.py` +
`batch/verification/run_diagnose_guard_over_request.batch`. Item 13 asks for
"16 modes at res 60, guard on vs off, compare the low bands and the wall time"
— but that comparison is precisely the one the item itself calls unclean:
`EMESolver.compute_modes()` over-requests `num_modes + _DUP_MODE_MAX_EXTRA`
(16 → 24) in a single `compute_modes_meep()` call when the guard is on, so the
flag moves both the eigenproblem and the filtering at once.

Three arms separate them — **A** guard off/16, **B** guard on/16, **C** guard
off/24 (B's eigenproblem, unfiltered):

- **C vs A** isolates the artifact question: does over-requesting *alone* move
  the low bands? Non-zero means every guarded EME number carries an
  over-request artifact independent of any filtering; inert means the guard is
  a pure filter.
- **B vs C** shows what the filter drops — which of C's 24 were rejected as
  duplicates and what the guard reached further up the spectrum to replace
  them with.

Wall time is recorded per arm, because if C vs A is inert the `C/A` ratio is
the entire argument for cutting `_DUP_MODE_MAX_EXTRA` on cost grounds — item 13
notes job 11575516's 3h36m was dominated by this over-request. Run at res 60:
above the shipped 40, and where 11575516 measured meow as converged, so a
difference cannot be blamed on an unconverged mesh.

### Added — dispersive FDTD media, re-parameterised and opt-in (material model step 2, work items 1–3)

`material_model_validation.md` step 2 was implemented once and reverted: the Sellmeier
form puts ε∞ = 1 with the whole index in susceptibilities, meep leaves those unsmoothed,
and the resulting `epsilon diagonal = (1,1,1)` structure lost +0.0496 in n_eff to
disabled subpixel averaging — 25× the ~0.002 the dispersion was worth. Work items 1–3
build the re-parameterised version that keeps smoothing. **Off by default** (the gate had
not run when this was written; it has since passed — see below, and note λ₀ has moved twice
since this section).

- **`reparam_lorentz(kind, band, lam0=None)`** (`beef/materials/stack.py`; λ₀ was a
  fixed 10.0 when first written — see the two later sections) — least-squares
  fit of `ε(λ) ≈ ε∞ + σ·λ²/(λ² − λ₀²)`, linear in both parameters at fixed λ₀, so a
  deterministic `lstsq` rather than an optimiser in the medium-build path. The fitted form
  *is* one meep `LorentzianSusceptibility` by construction — which is why re-fitting works
  where peeling ε∞ off analytically cannot (the remainder carries a λ₀² numerator meep has
  no shape for). Verified end to end: meep's own `Medium.epsilon()` tracks Sellmeier across
  the C-band to 1.96e-4, exactly the fit error.
- **`meep_medium_reparam()`** — ε∞ carries 104% of the true Si/SiO2 ε contrast (against 0%
  today), so smoothing acts again; the single Lorentzian carries only the dispersion.
- **`MaterialStack(dispersive=False)`** — new opt-in flag. `si_medium()`/`sio2_medium()`
  now take a scalar *or* a band and only disperse when the flag is on and a real band was
  passed, so single-wavelength mode solves are unaffected by construction.
- **Cache fingerprint** — `resolved_material_config()` records `[ε∞, σ, λ₀]` per material.
  Without this, flipping the flag would leave the key **identical** (the probe indices come
  from `sellmeier_index()`, which the re-parameterisation does not change) and every stale
  entry would keep matching — the defect step 0 was created to kill, one level up.
- `examples/verification/python/verify_item_8_reparameterisation.py` reproduces every table
  in the validation doc. Pure numpy, no meep/meow import.
- 14 new tests, `tests/test_material_reparam.py` (28 by end of day). Fast suite 447
  passed at this point, 468 by end of day.

**Solver-level flags (follow-up, same day).** `dispersive=` now exists on `FDTDSolver` and
`EMESolver`, so the setting no longer has to be baked into a shared `MaterialStack` — the
`SOI_160`/`SOI_220` singletons are process-wide, so mutating one to A/B a run would leak
into every other solver.

The two defaults are **deliberately asymmetric**, because the engines start from opposite
places. `FDTDSolver(dispersive=None)` inherits the stack (off by default) — that engine is
the one with the frozen-band problem. `EMESolver(dispersive=None)` means **dispersive**,
and only an explicit `False` freezes it: EME has always been exact per wavelength via
`SampledMaterial` + one `Environment` per wavelength, so inheriting an off-by-default stack
would have silently downgraded every existing EME run. `EMESolver(dispersive=False)` builds
a flat `SampledMaterial` pinned at band centre — strictly a downgrade, and it exists only so
a hybrid A/B can hold *both* engines frozen the same way.

`HybridSolver` gets no flag of its own: it composes pre-built solvers, so it inherits
theirs, and a hybrid-level setting could only work by mutating the caller's objects. It
gets a guard instead — `_check_material_consistency()` warns when sections disagree about
dispersion, naming the sections and the worst n_Si gap. That is the cheap half of work
item 5; the hard assert waits on the strip test to fix a tolerance.

**The guard fires on the default configuration, and that is the point.** A broadband
`FDTDSolver` frozen at band centre chained to an `EMESolver` solving per wavelength have
always been contracting modes computed in two different materials — 0.00786 apart in n_Si
at the O-band edges. That mismatch predates all of this work; nothing had ever checked it.
It also sharpens the doc's "engine agreement regresses 4e-16 → 2e-4": that 4e-16 measures
agreement between the two *index models*, not between the media actually simulated. On the
latter, enabling dispersion improves FDTD-vs-EME agreement ~16× (0.00786 → 4.9e-4). Work
item 5's guard must therefore compare built media, not `si_index()` lookups — the latter
reads 4e-16 in every configuration, including today's broken one.

### Measured — 1×2 MMI check 4: hybrid vs full-device FDTD (job 11585850)

**The repo's first hybrid-vs-full-device comparison on an MMI**, closing
`mmi_1x2_validation.md` gap 5 as a measurement. Both data tasks (11585845_1/_2)
plus the comparison completed end to end.

| | max \|ΔT\| | mean \|ΔT\| | rel. to mean T_fdtd | worst λ |
|---|---|---|---|---|
| `out_top` | 0.03588 | 0.02530 | **9.11%** | 1.33 µm |
| `out_bot` | 0.03590 | 0.02532 | 9.11% | 1.33 µm |

`check_passive()` true both sides. The two ports agree to four decimals, so
**splitter symmetry survives the cascade** — worth recording as a negative
result, since interface mode-mixing would not be expected to preserve it.

**The disagreement is systematic, not scatter: every ΔT is negative**, growing
0.021 → 0.036 from 1.26 to 1.33 µm. The hybrid under-predicts at all 11
wavelengths. That is the Δβ-discretization signature from
`resolved/directional_coupler_validation.md` items 6/14, where the coupler sat
~5% off uncalibrated and **0.11%** with `HybridSolver(calibrate_body_neff=True)`
— a flag this device has never been run with. Trying it is the obvious next
step and costs one hybrid run, not a new script.

Two caveats against reading 9.11% as the device's true error: both examples
still share `RESOLUTION = 40.0`, which gap 2 shows is unconverged here (and this
is the worst-converging device audited), and gap 4 (engine agreement) is still
unmeasured, so there is no basis to split the number into calibratable versus
structural. No threshold enforced, deliberately.

### ✅ Strip-test gate PASSED — and it moved λ₀ from 10.0 to 2.0

`examples/verification/python/verify_item_8_strip_gate.py` (new). Three arms on a
450×220 nm Si strip at res 40, 1.55 µm — A non-dispersive, B re-parameterised, **C the
reverted ε∞=1 attempt**. C exists so the harness can prove itself: it reproduces the
recorded +0.049640 failure to 3e-4, which is what makes B's number evidence about the fix
rather than about the test. Runs in ~4 s on a login node, no queue slot.

**The gate initially FAILED at the shipped λ₀ = 10** (Δn_eff −0.00222 C-band, −0.00258
O-band, tolerance 0.002). That is the doc's own risk firing — *"the ε∞ contrast is 104% of
the true contrast… whether that biases n_eff, and in which direction, is not predictable
from the fit"* — and the answer is: it biases **low**.

λ₀ does not really control fit quality; it controls how the permittivity splits between
the **smoothed** ε∞ and the **unsmoothed** susceptibility, and only that split decides
whether subpixel averaging still works:

| λ₀ | Δn_eff C / O band | σ (C/O) | fit, % of band swing |
|---|---|---|---|
| 10.0 | −0.00222 ❌ / −0.00258 ❌ | 17.6 / 37.7 | 2.4 / 3.3 |
| 3.0 | −0.00160 ✅ / −0.00203 ❌ | 0.89 / 2.30 | 3.1 / 3.9 |
| **2.0** | **−0.00069 ✅ / −0.00123 ✅** | **0.12 / 0.51** | **5.5 / 5.2** |

2.0 is the only value clearing both bands, so `_REPARAM_LAM0_UM` is now 2.0. Fit quality
is the *wrong* criterion — it improves monotonically with λ₀ while the gate result worsens.
Cost is ~3% more of the band's dispersion given up, still keeping ~95%.

Two further consequences. σ drops two orders of magnitude, which **retires** the standing
"meep stability is unverified at σ = 17.6" risk (no stability warnings logged). And λ₀ = 2.0
sits inside `_MODEL_RANGE_UM` (1.10–2.10), so `reparam_lorentz()` now **raises** if a band
comes within 10% of the pole instead of returning a quietly near-resonant model; the widest
band the examples span (1.26–1.70) clears it by 18%.

Also recorded: arm A lands 0.0017 off the doc's 2.313398 baseline, because this harness is
a bare strip in a 2×2 µm cell rather than the original ad-hoc setup. Δ must be read within
a run, not as absolute n_eff across runs — which is exactly what arm C validates.

**Items 4 (resolution sweeps) and 5's hard assert are unblocked.** 3 more tests
(28 total in `tests/test_material_reparam.py`); fast suite 459 passed, 14 skipped.

**Three doc corrections found by implementing it.** The prescribed call sites
(`fdtd.py` 526/776/1175) were stale — the real ones are 615/616, 930/931, 1581/1582, and
only **two** should change: 930 is a single-wavelength mode solve where the frozen index is
already correct. The O-band and full-range fit errors are 4.9e-4 and 5.6e-3, not 4.7e-4 and
5.5e-3. And λ₀ is not strictly degenerate — fit quality improves slightly and converges as
λ₀ grows (2.11e-4 → 1.92e-4 over λ₀ 5 → 50 µm) rather than being identical; σ still swings
3.78 → 461 across that range, so the pinning argument stands unchanged.

### Docs — item 6 naming collision, v1 ladder script marked superseded

Two docs both number an item 6, unrelated to each other:
`fdtd_vs_EME_mode_validation.md` item 6 (`interface_modes` convergence, open,
no tooling) and `resolved/directional_coupler_validation.md` item 6 (engine
agreement, closed). `verify_item_6_engine_agreement.py` and
`run_item_6_res_sweep.batch` are the latter only. Added a cross-reference note
at each doc's item 6 so "the item 6 script" can't resolve to the wrong one.

`diagnose_kpoint_ladder.py`'s docstring now states up front that it is
superseded by `diagnose_kpoint_ladder_v2.py` — kept, not deleted, because its
numbers (3 ties/17 distinct fixed, 1 tie/19 distinct naive ladder) are quoted
verbatim elsewhere in this file and in notes.md.

### Docs — job IDs refreshed after the re-queue

`validation_markdown/running_jobs.md` re-snapshotted (07-31 08:22) with the
old→new ID mapping, the billing diagnosis, and current states; the per-job
detail sections and `mmi_1x2_validation.md` gap 5, `mmi_2x2_validation.md`
gap 5 + `Order`, `modal_multiplexer_validation.md` gap 6, and
`validation_markdown/README.md`'s status table all updated to the new IDs
(including the `reports/11585847_0.out` path to watch). Dates that read
"in flight as of 2026-07-30" corrected to 07-31 where they now name a 07-31
job.

### Measured — item 13 closed: the duplicate guard's over-request is inert, and costs 3.2× (job 11586098)

Job 11586098 completed in 1:54:50. Modal multiplexer interface 0, λ = 1.31 µm,
res 60, 16 modes, `_DUP_MODE_MAX_EXTRA = 8`.

| comparison | result |
|---|---|
| **C vs A** (24 requested vs 16, both unguarded) | max \|Δn_eff\| = **0.000e+00** across all 16 low bands |
| **B vs C** (guarded-16 vs unguarded-24) | B = C's first 16 exactly; worst residual 0.000e+00; **no index skipped** |
| **cost** | A 934.0 s, B 2948.5 s, C 2997.7 s → `B/A = 3.16×`, `C/A = 3.21×` |

**The over-request does not move the eigenproblem.** Not "to within 1e-6" —
exactly zero on every band. So `EMESolver.compute_modes(guard_duplicate_modes=
True)` is a pure filter at production scale, item 0's res-12/3-mode inertness
check does generalise to 16 modes on a res-60 mesh, and **no guarded EME number
recorded anywhere in `validation_markdown/` carries an over-request artifact**.
That was the open risk item 13 was raised for; it is closed.

**What follows is a cost decision, deliberately not taken yet.** `B/A` and
`C/A` agreeing to within 1.6% is the proof that B's expense is the extra 8
bands and not the filtering, so `_DUP_MODE_MAX_EXTRA = 8` is ~2/3 of the wall
time of every guarded `compute_modes()` call — including the 3h36m of job
11575516 — in exchange for headroom the observed skip counts (≤ 3 on this
device) never use. Left at 8 for now on two grounds: the skip-count evidence is
from one device, and item 12 may replace the guard with a k-seeding fix that
makes the constant moot. Recorded so the change is a one-line decision when
those resolve, not a re-measurement.

**Unplanned third finding — no duplicate fires at res 60.** Arm B filtered
nothing, on the same cross-section where item 0 found meow's modes 10/11
identical at res 40 and where the production run below fires the guard at res
40. Consistent with job 11575516's meow-converged-by-60 plateau: the
near-degenerate pair looks resolved into two genuine eigenvectors by res 60.
Suggestive only — one cross-section, one wavelength, and the 24-band list still
holds tight clusters above band 16 that a 16-mode request never reaches. It
does turn "does the duplicate survive mesh refinement" into a real question for
item 12.

### Found — the meow duplicate guard fires on the *production* path (job 11585847, partial)

Job 11585847 (modal multiplexer hybrid, this device's first real run in this
tree) is still running, but its log already answers something the audit scripts
could not. After the first two FDTD subregions, the first EME body solve logged
— once per rank on all 96 ranks, one distinct event —

```
EMESolver.compute_modes: meow mode 9 duplicates already-accepted mode 8
(self-overlap 1.000000 > 0.90) -- skipping it and reaching further up the spectrum.
```

This is the item-0 guard on `EMESolver`'s production `compute_modes()`, at the
shipped `RESOLUTION = 40`, doing real `HybridSolver` work rather than replaying
gap 1's field-overlap check (which is how item 0 was verified on 07-30). Absent
the guard, this run's 16-channel interface basis would have been rank-deficient
by one channel, on the device whose entire purpose is keeping 16 channels apart.

`warn_on_transfer_duplicate` — the item-8 post-hoc `T`-column scan this job was
submitted to exercise — has **not** fired; there is no `transfer` line in the
log at all. That is the expected outcome if the upstream guard removes the
duplicate before it reaches `T`, but a detector that has never fired and a
detector that cannot fire are indistinguishable, so
`modal_multiplexer_validation.md` gap 6's original ask stays open. Three
interfaces and ~2h45 of wall remain.

### Fixed — job 11586046 hung on an eigensolve with no solution; script now stops at the guided edge (resubmitted as 11588309)

Job 11586046 (items 11 + 12) printed none of `diagnose_kpoint_ladder_v2.py`'s
four per-policy blocks in 2h30 and was cancelled. **It was not slow, it was
stuck** — an earlier revision of `running_jobs.md` diagnosed it as a wall-clock
budget miss, which was wrong. Reconstructing the converged `k_x` per `band_num`
from its MPB trace (target f = 0.763359, cladding n = 1.446804):

| band_num | n_eff | guided? | |
|---|---|---|---|
| 1–9 | 2.7064 → 1.7637 | yes | |
| **10, 11** | **1.7477, 1.7477** | yes | duplicate |
| 12, 13 | 1.6836, 1.5936 | yes | |
| **14, 15** | **1.4848, 1.4848** | yes | duplicate |
| 16–18 | 1.4337 → 1.3877 | **no** | below cladding |
| **19, 20** | **1.3491, 1.3491** | no | duplicate |
| 21, 22 | 1.3226, 1.3075 | no | |
| **23** | — | — | **no k reaches f** |

`N_TARGET = 20` *distinct* modes needs band_num 23 because of the three
duplicate pairs, and band 23 has no solution at the target frequency: MPB's
k-search limit-cycles between k = 0.2346 (f = 0.4775) and k = 1.7480
(f = 1.0283) indefinitely. It burned **161 MPB solves on that one band — 59% of
the job's entire eigensolve work** — against 3–8 solves for each of bands 1–22.
Bands 1–22 took ~58 min for 22 attempts, so the original ~25 s/solve budget was
about right; no wall clock would have helped, and the 4-task split briefly
recommended would have produced four hung tasks.

**Script fix.** `N_TARGET` 20 → 16 (what `modal_multiplexer_hybrid.py` actually
requests; the 20 was for comparability with v1's raw sweep, whose bands 16–20
were unguided anyway), and the search now **stops when a raw solve returns
`n_eff ≤ the cladding index`** — past there nothing is guided and the next
band_num is the non-convergent one. This is a *physical* bound, not a timeout:
MPB runs in C, so a Python signal or alarm cannot interrupt a solve already
entered. The per-policy report now prints an explicit stop reason
(`target` / `unguided` / `attempts`), because "found 13 of 16 because only 13
exist" and "found 13 because it ran out of attempts" are different findings and
must not read alike. Termination checked offline before resubmitting, meep-free:
11586046's own measured n_eff sequence replayed through `_guarded_search` with
fake `sim`/`mp` objects and the real `field_overlap_matrix` — all four policies
stop at band 16, none requests the non-convergent band, and policy A reproduces
the measured 13 accepted / 16 attempts / skips at 11 and 15.

**Batch fix.** `run_diagnose_items_11_12.batch` now runs **item 11 first** and
budgets 360 min (was 240). Ordering item 12 first was deliberate — protect the
higher-value result from a wall-clock overrun — and it is precisely what cost
both items: item 12 hung so it printed nothing, *and* item 11 never started.
Item 11 is ~50 min with every band known to converge, so it banks a guaranteed
result first; item 12 prints per policy and so degrades gracefully.

### Fixed — four constants had drifted between notebooks and scripts, and one drifted inside an unrelated commit

Found by diffing every notebook's constants against its script — an audit, not a
failing test, because nothing checks this.

| file | notebook | script | which was wrong |
|---|---|---|---|
| `directional_coupler_hybrid` | `N_MODES = 6` | **5** | notebook |
| `mmi_1x2_hybrid` | `N_MODES = 12` | **8** | notebook |
| `mmi_2x2_hybrid` | `N_MODES = 15` | **14** | notebook |
| `directional_coupler_fdtd` | `COUPLING_LENGTH = 15.13` | 15.15 | **script** |

**The first three: item 9 was not the cross-device pass it claimed to be.** It
dropped each device's unguided mode in every `.py` and touched no notebook. The
1×2's stale 12 matters most — that is precisely the configuration
`mmi_1x2_validation.md` gap 3 identified as the failure mode (SOI_160's top bands
piling up at cutoff), so anyone running that tutorial reproduced a known-bad
basis. Notebooks synced to their scripts; one line changed in each, no
reformatting churn.

**The fourth is the reverse, and worse.** `COUPLING_LENGTH` was deliberately set
to **15.13** in four places, with reasoning recorded in
`archive/dc_validation_20260727.md` item 14: it comes from item 5's ground-truth
crossing (15.129 µm) rather than the hybrid's grid-snapped 15.151, "since three
of the four users are full-device FDTD runs and that reference is the finer
measurement". Commit `e85b350` — message *"Remove cached example simulation files
for ring resonator hybrid solve_neff BareStrip configurations"* — moved
`directional_coupler_fdtd.py` to 15.15 and said nothing about it. The other three
sites never moved, so this was a stray edit, not a decision. Restored to 15.13
with the reasoning inlined at the constant so the next drift has to argue with
it. No stored output was generated at 15.15 (this script's outputs are still the
owed-regeneration ones at 13.5 µm), so nothing downstream needs redoing.

Also cleared two stale doc entries: `resolved/directional_coupler_validation.md`
still listed "Drop `N_MODES` 6 → 5" as an unbundled cleanup that item 9 had
already done on 07-30, and four coupler example docstrings still described
`N_MODES=6` in prose while their code said 5.

### Measured — items 11 and 12 both resolved; k-seed laddering recovers a guided mode the shipped seed loses (job 11588309)

The reshaped job completed in **34:55** — less than the hung predecessor had
burned on a single non-convergent band. Both items it carries are answered.

**Item 11 → CLOSED. The guard's diagnostic grid is clear at res 40 as well.**
Same four grids, at the shipped `RESOLUTION` this time:

| grid | points in core | skipped band_nums | worst overlap among accepted |
|---|---|---|---|
| 16×8 (shipped) | 8/16 | 11, 15 | **0.3877** |
| 32×8 | 18/32 | 11, 15 | 0.1924 |
| 64×16 | 36/64 | 11, 15 | 0.1201 |
| 128×16 | 72/128 | 11, 15 | 0.1209 |

All four agree exactly — same skips, same 16 accepted band_nums, identical n_eff
— and both skips are exact n_eff ties, so they are confirmed genuine by a
threshold-free test that needs no field sampling at all. The margin *did* narrow
with resolution as the item predicted (0.3392 at res 20 → 0.3877 at res 40, on a
0.9 threshold), which is now recorded as the quantity to re-check if this device
ever moves to res 60+. Two resolutions is enough to close it as not live.

Incidental: the third duplicate pair moves band position with the mesh — res 20
skipped 11 and 16, res 40 skips 11 and 15. Both are exact ties at their own
resolution, so both skips are correct; pair positions are mesh-dependent, not
fixed properties of the device.

**Item 12 → the fixed k-seed is silently losing a guided channel.** Four seeding
policies under a real accept/skip/retry search, modal multiplexer interface 0,
res 40:

| policy | distinct guided found | attempts | skips | stopped because |
|---|---|---|---|---|
| **A** fixed (today's code) | **13** of 16 | 16 | 2 — bands 11, 15 | guided spectrum exhausted |
| **B** naive ladder | **14** of 16 | 16 | 1 — band 10 | guided spectrum exhausted |
| **C** step ladder (−0.005) | **14** of 16 | 16 | 1 — band 10 | guided spectrum exhausted |
| **D** accepted ladder | **14** of 16 | 16 | 1 — band 10 | guided spectrum exhausted |

B, C and D return A's set **plus** n_eff 1.512589 — a mode absent from A
entirely, not reordered and not duplicated. It reproduces v1's recovered
1.512589 exactly, now with the duplicate guard in the loop rather than a flat
unguarded sweep.

- **The loss is silent today.** A's two skips are both genuine duplicates,
  correctly rejected; it still ends a mode short because bands 14 *and* 15 both
  converge to 1.484797 from a constant seed. Nothing in the log distinguishes a
  retry that found something new from one that re-found a neighbour — on the
  device whose whole purpose is keeping 16 channels apart.
- **Laddering costs nothing**: same 16 attempts, one fewer skip. The "strict
  improvement" case the script was written to detect, so the fix direction is
  settled — seed band *k* from band *k−1*, not from a constant.
- **C and D are unnecessary here.** All three laddered policies gave
  byte-identical sets, attempts and skips. The self-reconvergence risk that
  motivated them did fire (all three skip band 10 at overlap 1.0, raw n_eff
  1.763693 = band 9's answer), but the existing guard absorbed it and the retry
  landed on a real mode. Not adopting the extra complexity on one cross-section's
  evidence.
- **Every policy stopped on the new cladding test**, none on the attempts
  backstop — so no policy had duplicates eat the +8 headroom production also
  relies on, and "found 13 because only 13 were reachable" is distinguishable in
  the log from "found 13 because it ran out of attempts". That distinction is what
  makes these counts readable as physics.

**A correction this forced, in six files.** Yesterday's reading — "only 13
distinct guided modes exist at res 40, so the 16-mode basis is unachievable by
any seeding policy" — was inferred from the hung job's fixed-seed trace, i.e.
from the enumeration under test. At least 14 are reachable. The shortfall against
16 is real and the guard still pads it with unguided bands (item 14), but part of
it is a software defect, which is the opposite of the conclusion first recorded.

**Also noted against the script, not fixed mid-run:** the cross-policy summary
reports "DIFFERENT modes, max |Δn_eff| 0.027793" for B/C/D against A. That is an
artefact of comparing sorted 13- and 14-element prefixes when an inserted mode
shifts the alignment — B ⊃ A exactly. The comparison should align on set
membership, not position, when the counts differ.

### Measured — 2×2 MMI gap 6 closed: "the cleanest device audited" survives both duplicate guards (job 11588569)

`batch/verification/run_audit_mmi_2x2_gap6.batch` (new), 20:18. The 2026-07-26
result behind that phrase was measured at `N_MODES = 15` before either duplicate
guard existed, and is quoted in `README.md` Findings 1 and 2 plus three other
files. Re-run at 14 with both guards live:

| | 2026-07-26 (11500536) | 2026-07-31 (11588569) |
|---|---|---|
| modes requested | 15 | **14** |
| guided | 14 of 15 | **14 of 14** |
| correspondence | identity `[0 … 13]` | **identity `[0 … 13]`** |
| worst guided \|Δn_eff\| | 0.0007 | **0.0007** |
| duplicate guards | did not exist | **both live, neither fired** |

Three things settled. The claim is current again rather than dated. **Neither
guard fired at all** — the specific worry was that this is the repo's
highest-mode-count SOI_220 device on the 2.8 µm body
`fdtd_vs_EME_mode_validation.md` item 11 identifies as the meep guard's
least-margin case, and it turns out to have nothing to fire on (item 11 closing
the same day makes that a positive result, not merely an absence). And dropping
the unguided band cost nothing measurable: worst guided |Δn_eff| is unchanged,
which is what item 9 asserted and had never been checked here.

Fixed on the way in: the audit's own `mmi_2x2` config still requested
`n_modes=15` while its comment claimed to be pinned to `mmi_2x2_hybrid.py`, which
moved to 14 on 07-30. Corrected.

### Fixed — `section_mode_neff()` ignored `dispersive=` entirely, so job 11586986 measured nothing

Found by reviewing the running jobs' output rather than by a test. Job 11586986's
"dispersive" column came back **byte-identical** to its non-dispersive one — same
n_eff, Δn_eff and Δβ to six decimals at every resolution it reached.

`FDTDSolver._build_section_mode_sim()` built its media with
`si_medium(wl_center)` **directly, bypassing `_media()`**, and so never read
`self.dispersive`. That was a deliberate choice, recorded in
`material_model_validation.md` work item 2 as *"930 is a single-wavelength mode
solve, where the frozen index at that exact wavelength is already the correct
medium"*. True for `section_mode_fields()`, which takes one wavelength. **False
for `section_mode_neff()`**, which accepts a whole band and loops it on one
structure — for that caller the frozen centre index is precisely the defect work
item 2 fixed for `_build_domain()`/`solve_section()`, left in place by an
argument that did not apply to it.

`_build_section_mode_sim()` now takes `wavelengths=` and routes media through
`_media()`; `section_mode_neff()` forwards its band. Scalar callers
(`section_mode_fields()`, the three diagnostic scripts, the pad tests) are
byte-unchanged, and structurally rather than by convention:
`MaterialStack._medium()` requires **both** the flag and a real band
(`size > 1 and min < max`), so a scalar cannot reach the dispersive path at all.

Measured on a real solve, coupler body, SOI_220, 1.50–1.60 µm:

| | max \|Δn_eff\| over the band | at centre 1.55 | at edge 1.60 |
|---|---|---|---|
| before | **0.000e+00** | 0 | 0 |
| after | **5.29e-03** | 5.6e-05 | 4.6e-03 |

Agreement at band centre with divergence at the edges *in opposite directions*
is the signature of unfreezing a centre-pinned index — that is what says the fix
does the intended thing and not merely something. meep also now reports
`lorentzian susceptibility: frequency=0.543478` (λ₀ = 1.84 µm = 1.15 × band max,
the band-derived pole from `4296ee7`) and `Dispersive medium detected on
eigenmode monitor`, neither of which appeared before.

**9 new tests**, `tests/test_fdtd_mode_sim_dispersive.py`, meep-free. They pin
both halves — band caller gets per-wavelength media, scalar callers unchanged in
either configuration, degenerate "bands" (one point, or min == max) stay inert —
plus the one that would have caught this: `section_mode_neff()` must forward its
whole band, not `mean(wavelengths)`. Verified by mutation: the two key tests fail
against the pre-fix code while the seven unchanged-behaviour tests still pass,
which is what establishes the change is confined to the band path. Fast suite
**477 passed, 1 skipped**.

`verify_item_8_dispersive_convergence.py` also no longer runs vacuously: it
solves once per setting at the lowest resolution **first** and aborts non-zero if
the two agree exactly, so an inert independent variable cannot be recorded as a
measurement again. It additionally reports Δn_eff at a band **edge**, because
dispersion is near-inert at band centre by construction — a smoke run gives
+1.7e-5 at centre against +5.9e-4 at the edge, so the centre-only reading this
script used understates the effect ~34×.

Material-model work item 4 is therefore back to **not started**; its record has
been corrected in three places, including an earlier note of mine that blamed a
stale λ₀ (irrelevant — λ₀ cannot matter when the dispersive medium is never
built).

**Salvaged:** 11586986's non-dispersive column is valid and preserved in work
item 4. It carries an unrelated finding worth having — Δβ's convergence step
collapses to +1.5e-05 at res 100, then reopens to −4.08e-04 at 120. Res 100
looks converged and is not, on the exact quantity item 11's res-80 pin was set
on, in a sweep run to check that pin (`README.md` Finding 4).

### Found — this device's 16-mode basis is unachievable at res 40, and the production guard can hang (new item 14)

Two findings fell out of diagnosing the hang above, both larger than what the
job was sent to measure.

**The requested basis cannot be supplied at res 40 — and part of the shortfall
is a software defect, not the spectrum.** At the mux's interface 0, res 40,
today's fixed k-seed resolves **13 distinct guided modes** of the 15 guided band
slots; a laddered seed resolves **14** (job 11588309, below).
`modal_multiplexer_hybrid.py` requests `N_MODES = 16`. So the shipped basis is
short either way and the guard makes up the difference with bands below the
cladding index — but one of the missing channels is recoverable in software,
which is what item 12 was opened to determine. (This section first claimed the
ceiling was 13 and that no seed policy could beat it; that was inferred from the
fixed-seed trace and the re-run falsified it.) This sharpens
`modal_multiplexer_validation.md` gap 2 from "the shipped resolution is below
both engines' converged point" to "the shipped resolution cannot supply the
requested basis", and it is item 9's concern arriving through the guard rather
than through a device's `N_MODES`. Caveat kept with the number: whether the
(10,11) pair is a solver duplicate or a genuine near-degenerate pair meep cannot
separate is still `fdtd_vs_EME_mode_validation.md` item 1's open question, so 13
is "distinct as meep resolves them at res 40", not a physics claim.

**New item 14 — `_iter_distinct_bands()` can hang a production run.**
`beef/solvers/fdtd.py:1049` walks `band_num` up to
`n_modes + _DUP_BAND_MAX_EXTRA`, and its `RuntimeError` for "could not find
`n_modes` distinct eigenmodes" can only fire *after* a solve returns — so a
non-converging solve never reaches it and the loop hangs inside meep with no
exception and no timeout. For this device at res 40 the guard already reaches
band_num 19 (16 requested + 3 skips) against a failure point of 23: four
band_nums of margin, on a device whose duplicate count is known at only two
resolutions. Latent, not observed in production. **Recorded, not fixed** — the
eventual fix (the same cladding stop) changes which modes a basis *contains*,
not just how many attempts it makes, so it needs its own tests and a real
before/after on a device known to skip; the same scope argument
`modal_multiplexer_validation.md` gap 6 used when it declined to modify this
method under-verified.

---

## 2026-07-30

A single session, working down `validation_markdown/`'s gap lists. Three real
bugs fixed, **four** more found and recorded but not fixed (two of them by the
day's last measurement — see `Found and recorded, not fixed — two more` near
the end), one new opt-in validation feature added — and two bugs in that new
feature caught and fixed during its own verification. Sections below are in the
order they were written, so the later ones supersede earlier ones where they
overlap.

### Fixed — `FDTDSolver` returned duplicate eigenmodes for distinct bands

`sim.get_eigenmode()` can return the **same eigenvector** for two different
requested `band_num` values when bands sit close together, silently making an
`N_MODES`-mode interface basis rank-deficient — strictly worse than a
reordering, since the displaced mode is missing from the basis entirely.
Confirmed by self-overlap 0.9999999999 between two "different" bands on a
cross-section where every other pair overlapped ≤ 0.126.

- **Fix**: new `FDTDSolver._iter_distinct_bands()` — checks each candidate
  band's transverse-E overlap against every already-accepted band in the same
  call, skips and retries on a near-1 match (threshold 0.9, up to 8 extra
  attempts), raises if it cannot find enough distinct bands. Wired into
  `section_mode_neff()` and `section_mode_fields()`.
- **Verified**: inert on the directional coupler (identical results, runtime
  unchanged within noise); on the modal multiplexer the guard fires exactly
  twice (`band_num` 11 and 15) and the resulting n_eff sequence has no exact
  ties, against a pre-fix pair of exact ties at 1.747750 and 1.484797.
- **Not changed, deliberately**: `solve_section()`, the production
  interface-basis builder, reaches the same meep behaviour through a different
  API (`get_eigenmode_coefficients()` on DFT-decayed flux). Deferred with a
  concrete zero-cost follow-up recorded — a post-hoc check of its
  already-computed `T` coefficient columns for near-duplicates.
- Tests: `tests/test_duplicate_band_guard.py` (5, meep-free via a mocked
  `sys.modules["meep"]`).
- See `validation_markdown/modal_multiplexer_validation.md` gap 6, `notes.md`
  Solvers #12.

### Fixed — `CompoundGeometry` double-applied `stub_left`/`stub_right`

`_overlaps()` already expands the query window by the stubs, then also handed
the raw stub values to the first/last overlapping segment, which applied them
again. Whenever the padded window landed fully inside one segment (nothing
clipping the expansion back), geometry came out 2× the stub too long, or the
doubled coordinate raised a spurious "cut falls strictly inside a taper/bend".

- **Fix**: `_overlaps()` now returns only the part of each stub falling
  *outside* its segment. Both call sites already assumed the resulting
  invariant (`ov_lo - stub_left == lo`), so neither needed changing.
- **Reachability corrected**: first assessed as unreachable in production.
  That was wrong — `beef/plot/geometry.py` passes `stub_left + dpml` for FDTD
  sections, which for an *interior* section is `0.0 + dpml`. The modal
  multiplexer's section plots were drawing 10.5 µm where 9.5 µm is correct
  (0.5 µm of extra MMI body at each end). Plots only; no solver result was
  affected.
- Tests: `tests/test_geometry.py::TestCompoundGeometryStubIsAppliedExactlyOnce`
  (4 — fully-inside, reaching-outside, and the stub=0 path production uses).
- See `notes.md` Geometry #6.

### Added — `HybridSolver(check_mode_correspondence=True)`

`solve()`'s matrix-chain contraction assumes index *i* of one section's output
basis is the same physical mode as index *i* of the next section's input
basis, with **no runtime check at all**. This wires `beef/mode_identity.py`'s
field-overlap machinery into `solve()`/`solve_length_sweep()` as an opt-in
guard.

- Constructor params: `check_mode_correspondence`,
  `mode_correspondence_action` (`"raise"`/`"warn"`),
  `mode_correspondence_min_overlap`, `mode_correspondence_res`.
- Duck-typed and best-effort throughout: skipped with a log where a solver
  lacks a field-sampling hook or a domain size, and where the field query
  raises a geometry error (which is "could not check here", never "checked and
  found a mismatch").
- Supporting addition: **`EMESolver.section_modes()`**, exposing the raw
  `meow.Mode` list `solve_section()` already computes internally but never
  returned.
- Tests: `tests/test_hybrid_mode_correspondence.py` (17, meep/meow-free) +
  `tests/test_eme_length_sweep.py::TestSectionModes` (2, real meow).
- See `notes.md` Solvers #14.

#### Two bugs in the above, found by running it against real devices

Both produced confident-looking but entirely spurious "correspondence
failures". Recording them because the *pattern* is the lesson: a validation
tool that has never been run against a known-good reference is an untested
tool, and in both cases the diagnostic tell was the **shape** of the failure,
not its magnitude — the kind of thing a summary statistic hides.

1. **Diagnostic mesh was a fixed pixel count**, reasoning by analogy from
   `_iter_distinct_bands()`'s deliberately-coarse grid. That analogy is
   invalid: that function compares meep against *itself*, where identical
   eigenvectors overlap to exactly 1.0 on any grid. Cross-engine comparison
   has no such immunity — and an EME side **solves** on whatever mesh it is
   handed, so 24 pixels across 5 µm (4.8 px/µm, leaving SOI_160's 0.16 µm core
   *0.77 cells* tall) made one side's modes physically wrong. Tell: modes 0–6
   clean, 7+ scrambled. **Fixed** — mesh is now specified as a resolution,
   defaulting to the coarser of the two solvers' own `res`.
2. **No guided-mode filter.** Radiation modes (n_eff below the cladding index)
   have no engine-independent shape — their profile is set by the box and each
   engine's discretization — so they overlap ~0 however well the real basis
   corresponds. Per `validation_markdown/README.md` Finding 1, *every*
   interface basis in this repo carries at least one, so the check would have
   false-alarmed on every device. Tell: exactly the three sub-cladding bands at
   fault. **Fixed** — both sides filtered to `n_eff > sio2_index(wl)` before
   comparison, excluded count logged.

### Found and recorded, not fixed (at the time)

Both latent (not reachable by any in-tree configuration today), both
independently reproduced before recording. **Both fixed later the same
day — see the two "Fixed" entries below.**

- **`CompoundGeometry` assumes every segment's local frame is bbox-centred;
  `MMI`'s is body-centred.** With asymmetric per-port `taper_length` (a
  supported, documented `MMIPort` field), `compound.polygons()` hangs 1 µm past
  its own declared edge and `width_at()` returns a mid-taper width. Note
  `to_gds_eme()` is correct, so two paths in the same class disagree.
  → `notes.md` Geometry #7.
- **`section_mode_beta()` clamps for `dpml` draw padding in its seed path but
  not its fit path.** For short bodies the fit window overruns the uniform
  region by `dpml` at each end; reproduced as a real `ValueError` against a
  `DirectionalCoupler` with `coupling_length=3.0`. Unreachable today because
  `calibrate_body_neff` defaults off and both live callers use a ~29.5 µm body.
  → `notes.md` Solvers #15.

### Fixed — `CompoundGeometry` assumed every segment's frame is bbox-centered; `MMI`'s is body-centered

Follow-up to the gap recorded just above, same day. Every `CompoundGeometry`
device→local-frame conversion (`polygons()`, `section_prisms()`,
`section_gds()`, `width_at()`, `center_at()`, `source_aperture_at()`,
`_end_ports()`) assumed a segment's local-frame origin sits at its bbox
midpoint `(seg_lo + seg_hi) / 2`. True for every geometry except `MMI`,
whose frame is centered on the MM body — the two only coincide when
`max(input taper lengths) == max(output taper lengths)`.

- **Fix**: new `DeviceGeometry.frame_offset` property (default `0.0`,
  bbox-centred — unchanged behavior for every existing geometry),
  overridden in `MMI` as `(in_max - out_max) / 2`. Every `CompoundGeometry`
  conversion now adds `seg.frame_offset` to the bbox-midpoint `center` it
  already computed. `to_gds_eme()` needed no change (already correct).
- **Verified** against the notes.md repro case (`mm_length=5`, input taper
  3.0, output taper 1.0): `MMI.frame_offset == 1.0`; `polygons()`/
  `section_gds()` now stay within the declared bbox; `width_at()` at the true
  left edge returns the access-waveguide width 0.35, not the old mid-taper
  0.467.
- Tests: `tests/test_geometry.py::TestCompoundGeometryMMIFrameOffset` (4).
- See `notes.md` Geometry #7.

### Fixed — `section_mode_beta()`'s default fit-length spacing missed the `-2*dpml` clearance its seed span already applies

Also a same-day follow-up. `fdtd.py`'s `seed_span = min(2.0, span -
2.0*self.dpml)` exists because the mode solve pads its query by `dpml` on
each side — but the fit-length spacing a few lines below carried no such
clamp, so for a short uniform body the longest fit length's own
dpml-padded `solve_section()` window overran the region into whatever
taper/bend flanks it.

- **Fix**: `spacing` now derives from `usable_span = span -
  2.0*self.dpml` instead of the raw `span`, mirroring `seed_span` exactly.
  No new raise needed — `seed_span`'s existing check already guards the
  identical zero-clearance condition and fires first.
- **Verified** against the `coupling_length=3.0` repro: default lengths are
  now `[0.5, 1.0, 1.5]`, a drawn half-window of exactly 1.250 — flush with,
  not past, the 1.25 µm uniform half (was 1.75 before the fix).
- Tests: `tests/test_beta_calibration.py::TestSectionModeBetaDpmlClearance`
  (4, meep-free — `section_mode_neff()`/`solve_section()` swapped for
  numpy-only fakes), confirmed to fail against the pre-fix code (caught the
  exact 1.3333-vs-1.2500 overrun) before passing after.
  `TestTheSpacingPrecondition.test_the_actual_default_is_span_independent`
  updated to check the `dpml`-aware invariant.
- See `notes.md` Solvers #15.

### Example / doc changes

- `examples/mmi_1x2/python/mmi_1x2_hybrid.py`: `N_MODES` 12 → 8. At 12 the
  basis carried near-cutoff modes whose correspondence broke down (gap 3); at 8
  correspondence is identity with worst |Δn_eff| 0.0003.
- `validation_markdown/resolved/` created; `directional_coupler_validation.md`
  moved there.
- `notes.md`: Solvers #12–#15, Geometry #6–#7 added/updated.
- `README.md`, `validation_markdown/README.md`,
  `validation_markdown/modal_multiplexer_validation.md` updated to match.
- `examples/verification/python/audit_interface_modes.py`: extended to handle
  multi-interface `CompoundGeometry` devices; its `CompoundGeometry` workaround
  comment updated to note the underlying bug is now fixed.

### Test suite

362 → 376 passing (389 collected, 14 skipped). New: 5 duplicate-band-guard,
17 mode-correspondence, 4 CompoundGeometry stub, 2 `section_modes`.

**Further same-day additions (the two follow-up fixes above):** 376 → 403
passing, 14 still skipped. New: 4 `TestCompoundGeometryMMIFrameOffset`, 4
`TestSectionModeBetaDpmlClearance`, plus edits to one pre-existing
`test_beta_calibration.py` assertion (no new test) — a net +8 over the
figure above, since collection also picked up whatever else had landed by
that point.

### Open at end of session — ~~superseded~~, see `Measured — job 11575516` below

**The mode-correspondence check is not yet trustworthy on the modal
multiplexer.** After both fixes above it still reports worst-overlap 0.0000 at
interface 0, where `audit_interface_modes.py` measures 0.9735 on what should be
the same comparison. That discrepancy is unexplained; a meow-vs-meow isolation
run was in flight at end of session. Until it is resolved, treat this check's
output on that device as unverified — it is expected to flag *something* there
(gap 1 documents a genuine one-index shift), but not with a 0.0000 overlap.

*Resolved later the same day.* The two numbers were different measurements (the
audit calls a different meow mode solver), not contradicting ones; the
meow-vs-meow isolation landed as `self_overlap_diagnostics()`; and job 11575516
localized what remains to a single missing near-degenerate pair member. Kept
above as the state it was recorded in.

### Added — items 5 and 8 (`fdtd_vs_EME_mode_validation.md`)

- **Item 5 fixed** — `body_stub` (0.25) < `dpml` (0.5) was blocking
  `check_mode_correspondence` at the interior interfaces of 3 of 4 hybrid
  examples, since `FDTDSolver.section_mode_fields()` always self-pads its
  query by `dpml` on each side. `_build_section_mode_sim()`/
  `section_mode_neff()`/`section_mode_fields()` now accept a `pad` override
  (default `None` → `self.dpml`, behavior-preserving); `HybridSolver` gained
  `mode_correspondence_pad`, threaded to whichever solver's
  `section_mode_fields()` declares a `pad` parameter. Safe because
  `get_eigenmode()` never time-steps, so the PML built around a mode-field
  query is never actually exercised — a smaller clearance is exactly as valid
  as `dpml` for this call, just thinner. See the doc entry for the full
  reasoning.
- **Item 8's scoped zero-cost step done** — `beef.mode_identity` gained
  `transfer_matrix_duplicate_scan()`/`excited_mode_axis()`/
  `warn_on_transfer_duplicate()`: a post-hoc scan of `solve_section()`'s own
  `T` for near-parallel excited-side columns, wired into both
  `FDTDSolver.solve_section()` and `EMESolver.solve_section()` right before
  they return. Logged as a WARNING (symptom, not proof), never raised.
- Tests: 22 new (`tests/test_fdtd_mode_sim_pad.py`,
  `tests/test_hybrid_mode_correspondence.py::TestModeCorrespondencePad`,
  `tests/test_transfer_duplicate_scan.py`), all meep/meow-free. 403 → 425
  passing, 14 still skipped.

### Added — item 7 (`fdtd_vs_EME_mode_validation.md`): reusable cut-position invariance

- New top-level module `beef.cut_position_invariance` —
  `cut_position_invariance(smats, stubs)` fits transmission drift vs stub for
  every `(port, in_mode, out_mode)` key across a sequence of already-solved
  `SMatrix` objects, reporting the fitted slope as the headline metric.
  Generalizes `directional_coupler_hybrid_stub_sweep.py`'s one-off
  50:50-crossing sweep so any device's stub sweep can reuse the same
  comparison; per-device boundary/geometry bookkeeping stays with the caller.
- `SMatrix` gained `keys()` (all stored `(port, in_mode, out_mode)` tuples),
  needed to compare two SMatrix instances entry-by-entry.
- Exported from the top-level `beef` package.
- Tests: 8 new (`tests/test_cut_position_invariance.py`), meep/meow-free
  (pure `SMatrix`/numpy). 425 → 433 passing, 14 still skipped.
- Not yet done: migrating an actual device sweep script onto this function to
  confirm it reproduces the coupler's known 167x drift-collapse result —
  needs a real meep/meow run, not a code change.

### Docs — item 10 (`fdtd_vs_EME_mode_validation.md`): re-scope n_eff nearest-neighbour matching

- `tests/test_mode_identity.py::assert_modes_correspond()` re-scoped, not
  retired: module docstring, function docstring, and both real-device call
  sites now state plainly it is a fast, resolution-locked, shape-blind smoke
  check, superseded by `beef.mode_identity`'s field-overlap machinery (also
  reachable at runtime via `HybridSolver(check_mode_correspondence=True)`)
  as the authoritative mechanism. `beef/mode_identity.py`'s own module
  docstring updated to point back.
- Docs/docstrings only, no behavior change — `tests/test_mode_identity.py -m
  "not slow"` still passes unchanged (5 passed, 2 deselected).
- Also fixed staleness in three docs still describing item 5 as unfixed:
  `validation_markdown/README.md` (check 2's "still outstanding" paragraph),
  `validation_markdown/modal_multiplexer_validation.md` (gap 6, `Order`),
  `notes.md` (Solver #10(b), #12, #14).

### Measured — job 11575516: mode-spectrum resolution convergence (item 1)

No library change. `examples/verification/python/mode_spectrum_convergence.py`
ran to completion (3:36:13, 8 cpus) and settled item 1's open question, which
supersedes the `Open at end of session` section above.

- **The 0.0000 is explained** and was never a cross-engine failure: the audit
  script and the runtime check were measuring different things (different meow
  mode solver). `FDTDSolver`'s `get_eigenmode()` and
  `EMESolver(mode_solver="meep")` are both MPB on the same mesh and agree to
  ≤1.6e-4 on modes 0–7 at every resolution — identical to five decimals at res
  20/30/60.
- **The band-window hypothesis was half right, in the wrong direction.** The
  guard does push meep's window upward as predicted, but the surviving orphan
  is a *meow* mode meep's basis lacks, mid-basis, not a meep top mode without a
  partner. Over-requesting on the meow side — the proposed fix — would not have
  closed it.
- **What survives is narrower and structural.** At res 60 the entire
  disagreement is one mode: meep's 14 guided modes are exactly meow's 15 minus
  meow's index 8, matching pairwise to ≤4.9e-4. That mode is the upper member
  of a genuine near-degenerate pair (meow: 1.79727/1.79197 at res 60,
  1.79725/1.79694 at res 80 — converging together, and meow's own duplicate
  guard fired at neither resolution). meep's per-band ladder returns one member
  at **every** resolution 20–60, so a finer mesh does not reach it.
- **First per-engine converged resolution on SOI_160**: the EME side plateaus
  between res 60 and 80 (16/16 modes within 5e-3); the FDTD side is still
  moving at 60. The modal multiplexer's shipped shared `RESOLUTION = 40` is
  below both — `modal_multiplexer_validation.md` gap 2, now measured.

### Found and recorded, not fixed — two more (items 11–13)

Both opened by the run above, both about the duplicate-mode guards added
earlier the same day now being load-bearing on the production path.

- **The meep-side guard's duplicate check runs on a grid too coarse to judge
  the modes it is judging** (`_DUP_BAND_DIAG_NY = 16` × `_DUP_BAND_DIAG_NZ = 8`
  over the *monitor* width — 8 of 16 points inside a 2.8 µm core, 0.3125 µm
  pitch against a ~0.19 µm lobe pitch). Idealized slab profiles on that exact
  grid alias mode *m* onto mode *18−m* at η ≈ 0.9988, past the 0.9 threshold.
  Latent, not demonstrated: on the mux every band it skipped also has a genuine
  n_eff tie in the pre-guard audit, so its verdicts there look right — but not
  by construction, and a wrong skip logs identically to a right one. Item 11 /
  `notes.md` Solver #17.
- **Every band's eigensolve starts from the same k guess** (`_KPOINT_SEED_NEFF
  = 2.5` for band 1 and band 18 alike) — the likely root cause of the
  duplicates both guards work around, and the reason the meow side differs:
  `compute_modes_meep(num_modes=N+8)` is one simultaneous solve. Item 12 /
  `notes.md` Solver #18. Smaller companion (item 13): `guard_duplicate_modes`
  over-requests +8 modes, so guarded vs unguarded are different eigenproblems,
  and inertness was only checked at 3 modes/res 12.

### Docs — staleness sweep against current code (2026-07-30, late)

Verified doc claims against the shipped scripts and fixed what had drifted. No
behavior change; fast suite 433 passed, 1 skipped.

- `validation_markdown/README.md` — status table's 1×2 MMI check 2 still read
  `❌ FAILS — bands 8/9` (fixed that day via `N_MODES` 12→8) and the modal
  multiplexer's check 2 still read `not run`. Finding 1's `requested` column
  labelled as dated evidence, since all three counts have since moved
  (6→5, 12→8, 15→14). Finding 3 extended with the per-engine convergence
  measurement.
- `validation_markdown/mmi_2x2_validation.md` — gap 3's `N_MODES = 14`
  recommendation is **done** (`mmi_2x2_hybrid.py:124`), so `Order` step 3 no
  longer blocks goldens on it. New **gap 6**: this device's passing
  mode-identity result predates both duplicate guards and its own `N_MODES`
  change, so it describes a superseded configuration.
- `validation_markdown/mmi_1x2_validation.md` — recorded that the 2026-07-30
  re-audit's guard status is not in the record (guards landed the same day).
- `validation_markdown/running_jobs.md` — 11575516 marked completed with its
  result; the other five confirmed still queued, none started.
- `mode_spectrum_convergence.py` docstring and `running_jobs.md` both filed
  that job partly under item 6 (`interface_modes` convergence). It does not
  serve item 6 — different swept quantity, different observable. Corrected in
  both places and noted under item 6 itself.

### Measured — items 11 and 12, on real fields (2026-07-30, later)

Two login-node smoke scripts, following up the two open findings above. No
library change.

- **Item 11 (guard grid resolution) — checked, not live at res 20.**
  `examples/verification/python/diagnose_dup_band_grid.py`: replayed
  `_iter_distinct_bands()`'s accept/skip logic on the shipped 16×8 diagnostic
  grid and three finer ones (up to 128×16) against the same 20 real
  eigenmodes. All four grids agree exactly — same two skips, same n_eff — and
  the worst overlap among accepted bands (0.34 at the shipped grid) sits well
  under the 0.9 threshold. The idealized-slab aliasing estimate that opened
  this item does not transfer to real fields. Not yet re-checked at res 40.
- **Item 12 (fixed vs. laddered k seed) — tested, net positive, not closed.**
  `examples/verification/python/diagnose_kpoint_ladder.py`: 20 bands
  unguarded at res 40, fixed seed vs. seeding band *k* from band *(k−1)*'s
  own converged n_eff. Laddering cuts exact ties 3→1 and recovers 2
  genuinely distinct modes the fixed seed silently drops — real basis rank,
  not just reordering. It also reproduces the predicted risk: seeding
  directly at the previous band's answer can re-find that same eigenvector,
  and did, at a different band pair. Net win, not a fix. Also surfaced a
  third raw duplicate on this device's basis (band_num 19/20) beyond the two
  already on record.
- See `validation_markdown/fdtd_vs_EME_mode_validation.md` items 11–12 and
  `notes.md` Solver #17–18 for full tables and next steps (a smaller ladder
  step, or laddering off the last accepted band instead of the raw previous
  solve).
