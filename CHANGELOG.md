# Changelog

Dated record of library changes. Companion to `notes.md` (which holds the
standing recommendation lists and the reasoning behind each item) and
`validation_markdown/` (per-device validation status). Where an entry has a
fuller write-up, it points there rather than repeating it.

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
