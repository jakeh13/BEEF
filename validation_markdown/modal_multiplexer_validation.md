# Modal multiplexer validation — to-do

> **2026-07-28 — the material model changed under every device here.** Both engines
> now evaluate one analytic Sellmeier model (Salzberg Si, Malitson SiO2) instead of
> each doing its own nearest-tabulated lookup, which shifts n_Si by ~0.002 and
> removes a 0.0063 step that sat in the middle of the C-band. Any golden or stored
> result for this device predates that and needs regenerating; see
> `directional_coupler_validation.md` item 8 and, for the still-open dispersive-FDTD
> half, `material_model_validation.md`.



Opened 2026-07-26. This is the gap list against the four checks in
[README.md](README.md).

**Status 2026-07-30 (late):** gap 2 (resolution) is **measured** — job
11575516 swept both engines' `res` independently at interface 0 and found the
EME side plateaus between res 60 and 80 while the FDTD side is still moving at
60, so the shipped shared `RESOLUTION = 40` is below both. Details in gap 2.
That job also left a residue gap 1's fix did not remove: the two engines
disagree, at *every* resolution from 20 to 60, about one member of a
near-degenerate mode pair — meow finds both, meep's per-band ladder finds one.
That is a driver-level defect, not a mesh one, and it is tracked as
`fdtd_vs_EME_mode_validation.md` items 11–13 rather than here, since it is not
specific to this device. **Items 11, 12 and 13 all closed 2026-07-31** (jobs
11588309, 11586098) — the guard's diagnostic grid is clear at res 40, the EME
guard's over-request is inert, and k-seed laddering recovers one guided mode the
shipped fixed seed loses. Gap 6's exercise run (job 11585847) is **running** and
has already fired the meow-side guard six times on the production path; see
gap 6.

**Status 2026-07-30:** gap 1 (mode identity) was **run and FAILED — and the
root cause was not what it looked like.** All four interfaces disagree between
engines (the cross-section is the same body slab at all four, so the result
is identical at each — see gap 1 below), but the failure was dominated by a
**meep `get_eigenmode()` duplicate-band bug** (gap 6), not a genuine
mode-crowding or correspondence problem: two of meep's requested bands were,
provably, the exact same eigenvector returned twice. **Gap 6 is now fixed**
(for the audit-only entry points; `solve_section()` itself is deliberately
deferred — see gap 6) and the fix is confirmed to resolve this device's own
reproduction: no more exact neff ties, no more catastrophic η mismatches.
Details and evidence are in gaps 1 and 6.

`examples/modal_multiplexer/python/modal_multiplexer_hybrid.py` — the only script
for this device. SOI_160, not SOI_220 like the coupler and MMIs.

## Where it sits

The **most demanding** hybrid device here and the **least supported**:

- `N_MODES = 16` carried at **four** section interfaces (`:153`, `:208`) — more
  modes and more joins than any other example. Every interface is another place
  the FDTD and EME bases have to agree.
- **No full-device FDTD script.** Unlike both MMIs and the coupler, there is
  nothing to compare the cascade against as a whole; check 4 has no starting
  point at all.
- **No golden fixture.** Nothing pins its output.
- Different material stack (SOI_160), so no measurement from the other devices
  transfers numerically — not the coupler's converged EME res 80, not its Δβ
  ratio.

So it has the most ways to be wrong and the fewest ways to find out.

## Gaps

1. **Mode identity, at four interfaces with 16 modes each — RUN, FAILED, root
   cause found. DONE 2026-07-30.** `HybridSolver` contracts index *i* of an
   FDTD interface against index *i* of the EME body at every join, assuming
   they are the same mode. With 16 modes the higher-order ones are closely
   spaced and the two engines solve MPB on different grids, so reordering is
   more likely here than anywhere else — and a swap at *any* of the four
   interfaces corrupts the chain silently.

   This is the single highest-value check for this device. A modal
   multiplexer's entire purpose is keeping modes distinct, so a solver that
   mixes them is failing at exactly the thing being demonstrated.

   **Step 1 — n_eff nearest-neighbor audit** (`examples/verification/python/
   audit_interface_modes.py modal_multiplexer`, extended this session to
   handle `CompoundGeometry`'s four interfaces — see the script's own
   docstring). All four interfaces (`body1_start+BODY_STUB`,
   `body1_end-BODY_STUB`, `body2_start+BODY_STUB`, `body2_end-BODY_STUB`) sit
   on cross-sections of identical width/material (both MMI bodies are
   2.8 µm × 0.16 µm SOI_160), so all four give **identical** results — this
   was verified directly, not assumed, before treating one interface as
   representative of all four.

   | band | meep | meow | Δ | guided |
   |---|---|---|---|---|
   | 0–7 | 2.706–1.992 | 2.746–2.050 | ≤ small | yes |
   | 8 | 1.763693 | — | | yes |
   | 9 | 1.747750 | — | | yes |
   | 10 | 1.747750 | — | | yes |
   | 11 | 1.683602 | — | | yes |
   | 12 | 1.593585 | — | | yes |
   | 13 | 1.484797 | — | | yes |
   | 14 | 1.484797 | — | | yes |
   | 15 | 1.433726 | — | | **NO** (15/16 guided) |

   Nearest-neighbour n_eff correspondence: `[0,1,2,3,4,5,6,7,9,9,9,11,12,14,14]`
   — **not** the identity. **Worst guided |Δn_eff| 0.0640**, an order of
   magnitude worse than any other device audited (coupler 0.0023, mmi_1x2
   0.0067 pre-fix, mmi_2x2 0.0007). At face value this reads as "the worst
   basis mismatch in the repo" and matches this gap's original prediction
   that 16-mode crowding would be the failure mode.

   **Step 2 — field-overlap correspondence (the real diagnosis).** The n_eff
   read above turned out to be misleading about *why* it fails.
   `beef.mode_identity.field_overlap_matrix()`/`best_correspondence()` (new
   this session — see `notes.md`/`README.md`) compare mode *shape* (a
   normalized transverse E-field overlap) instead of the scalar n_eff, and
   are far more discriminating. Run at both a coarse comparison resolution
   (`cmp_res=20`, deliberately different from the audit's own res 40, to
   demonstrate the check needs no shared resolution) and at production
   res 40:

   | i (meep) | meep n_eff | → j (meow) | meow n_eff | η (overlap) | margin |
   |---|---|---|---|---|---|
   | 0–7 | 2.706→1.992 | 0–7 | 2.746→2.050 | 0.996–0.999 | 0.81–0.99 |
   | 8 | 1.763693 | **9** | 1.828169 | 0.9945 | 0.89 |
   | 9 | 1.747750 | **10** | 1.793915 | 0.9982 | 0.98 |
   | 10 | 1.747750 | **8** | 1.831949 | **0.0109** | **−0.99** |
   | 11 | 1.683602 | 11 | 1.729909 | 0.9972 | 0.97 |
   | 12 | 1.593585 | 12 | 1.639543 | 0.9956 | 0.97 |
   | 13 | 1.484797 | **14** | 1.525086 | 0.9735 | 0.93 |
   | 14 | 1.484797 | **13** | 1.571648 | **0.0035** | **−0.97** |

   (res 40 table shown; `cmp_res=20` agrees qualitatively — same shift-by-one
   from band 8 on, same two catastrophic mismatches — confirming the finding
   is not a resolution artifact.)

   Reading this row by row: bands 0–7 are clean, high-confidence identity
   matches (η ≈ 0.996–0.999) — genuinely fine, contrary to what the n_eff
   table's blanket "worst 0.064" suggested. From band 8 on, meep's index
   runs **one behind** meow's — but every individual match is still a
   **confident** one (η ≥ 0.97, large positive margin) **except** bands 10
   and 14, which are **not** confident at all (η ≈ 0.01–0.04, and — telling —
   a *negative* margin: the Hungarian assignment could not find any good
   column for them because the good one was already taken by their neighbor).

   **Root cause, confirmed directly: meep's `get_eigenmode()` returned the
   same eigenvector twice**, at bands 9≡10 and again at 13≡14. Checked by
   computing `field_overlap_matrix()` of meep's own fields against
   themselves (no meow involved): band 9 vs band 10 self-overlap =
   **0.99999999997**; band 13 vs band 14 self-overlap = **0.99999999990** —
   not "close," but numerically identical to 10 decimal places, while every
   other off-diagonal self-overlap in the same matrix is ≤ 0.126. meow's own
   n_eff sequence at this cross-section, by contrast, is strictly, uniquely
   descending with no ties anywhere (2.746, 2.717, 2.667, …, 1.525) — a
   genuine physical degeneracy would show up in *both* engines' spectra, or
   at least as an unusually small (not exactly zero) gap in both; seeing an
   *exact* tie in only one engine, at a specific pair of consecutive
   requested band indices, is the signature of an iterative eigensolver
   converging to the same eigenvector twice rather than finding the next
   distinct one — not a property of the physical cross-section.

   **What this means:** the interface basis this device requests at
   `N_MODES=16` is not actually 16 linearly independent FDTD modes — it is
   14 distinct modes plus two exact duplicates. `HybridSolver`'s contraction
   against those two duplicated columns is not just *mis-indexed*, it is
   *rank-deficient*: two of the 16 "channels" it thinks it has are the same
   channel counted twice, and whatever mode the duplicate silently displaced
   is missing from the basis entirely. This is worse than a simple
   reordering (which the n_eff table's framing suggested) and better in one
   respect (12 of the mismatched-looking bands, 8/9/11–13, are actually fine,
   just index-shifted by the two duplicates ahead of them) — the field-shape
   check is what separates these two very different failure modes; the
   n_eff-only check cannot.

   **Partially checked as of 2026-07-30** — see gap 6's "Recommended
   follow-up", now implemented as `fdtd_vs_EME_mode_validation.md` item 8:
   `FDTDSolver.solve_section()` (and `EMESolver.solve_section()`) now run a
   post-hoc scan of their own already-computed `T` for near-duplicate
   excited-side columns (`beef.mode_identity.transfer_matrix_duplicate_scan()`
   /`warn_on_transfer_duplicate()`), logging a WARNING if found. This closes
   the *detection* gap cheaply (no extra eigensolve) but not the *proof*
   gap — it is explicitly a symptom, not confirmation, so `solve_section()`'s
   production basis on this device has still never been run through it for
   real (that needs an actual meep run of this device's interfaces, not yet
   done). Note that `HybridSolver(check_mode_correspondence=True)` (see
   `notes.md` Solvers #14) does NOT answer this question either, even though
   it can now run against this device: it compares the two engines' *mode
   bases* at an interface via their own eigensolves, which is a different
   code path from `solve_section()`'s `get_eigenmode_coefficients()`-on-
   DFT-flux basis construction.

   **Recommended fix — since IMPLEMENTED, see gap 6 (done 2026-07-30):** a
   duplicate-band guard in
   `FDTDSolver`'s eigenmode-solving methods — after solving each band, check
   its self-overlap against previously-solved bands in the same call; a
   near-1 match means the same eigenvector came back and the solver needs to
   be nudged (e.g. a different seed `kpoint`, or explicitly skip ahead) to
   find the genuinely next one. Reducing `N_MODES` to dodge these two
   specific indices is not a real fix — the duplicate is a property of *this
   cross-section's exact band spacing at whatever resolution is used*, not a
   fixed pair of indices, so a different `N_MODES` or resolution could just
   as easily hit a duplicate somewhere else.

2. **Shared resolution knob — now MEASURED (2026-07-30, job 11575516), and
   res 40 is wrong for both engines.** `:86` passes one `RESOLUTION = 40.0` to
   both sub-solvers (`:167` FDTD, `:186` EME) —
   `directional_coupler_validation.md` item 11's defect. It had to be measured
   fresh on SOI_160 rather than copied from the coupler's res 80, since a
   160 nm core confines more weakly than 220 nm and its modes reach further
   into the cladding. That measurement now exists:
   `examples/verification/python/mode_spectrum_convergence.py` swept the two
   engines' `res` **independently** at this device's interface 0, 16 modes,
   1.31 µm.

   Modes agreeing between consecutive resolutions of the *same* engine, out of
   16, within |Δn_eff| ≤ 0.005:

   | | 20→30 | 30→40 | 40→60 | 60→80 |
   |---|---|---|---|---|
   | meow (EME) | 1 | 1 | 2 | **16** |
   | meep (FDTD) | 1 | 1 | 3 | *(not run)* |

   Guided-mode count at each resolution, of 16 requested: meep 13/13/13/14 at
   res 20/30/40/60; meow 14/14/14/15/15 at res 20/30/40/60/80.

   **Read it as three separate statements:**

   - **The EME side plateaus between res 60 and 80** — all 16 modes agree
     across that step. That is the first converged-resolution result for
     SOI_160 in this repo, and it is *not* the coupler's 80 arrived at by
     coincidence: 60 already suffices here, 40 does not.
   - **The FDTD side is still moving at res 60** (40→60 agrees on only 3 of
     16), so its converged resolution is above 60 and remains unmeasured.
     Extending the sweep to meep res 80/100 is the open action.
   - **res 40 is below both**, so every number this device has ever produced —
     including the gap 1 audit and gap 6's before/after — was taken on an
     unconverged mesh for both engines. That does not invalidate gap 1's
     *identity* findings (mode correspondence is a shape question, and the
     mode-8+ story is `fdtd_vs_EME_mode_validation.md` item 1), but it does
     mean no absolute n_eff or S-matrix number from this device is a converged
     one yet.

   **Sharpened 2026-07-31 — res 40 does not merely under-resolve this basis, it
   cannot supply it.** The "meep 13/16 guided at res 40" entry above was read as
   a resolution shortfall. Job 11586046's MPB trace shows the mechanism, and it
   is worse than that. Reconstructing the converged `k_x` band by band at this
   same cross-section (cladding n = 1.446804):

   - bands 1–15 are guided, but **10/11 and 14/15 are duplicate solves** — so
     **10/11 and 14/15 are duplicate solves under today's fixed k-seed**, so
     it finds only **13 distinct guided modes** at res 40 — while a laddered
     seed finds **14** on the same mesh (job 11588309), i.e. at least 14 exist
     and today's enumeration misses one;
   - bands 16–22 sit *below* the cladding index — unguided;
   - band 23 has no k that reaches the target frequency at all.

   13 is exactly the guided count the table above recorded, which corroborates
   it from a second job by a different route. What is new is what the other 3 of
   the 16 are: the duplicate guard reaches past band 15 and **pads the basis
   with unguided modes** to reach the requested count. So `N_MODES = 16` at
   `RESOLUTION = 40` does not describe a basis that exists — three of its
   channels are not guided modes of this cross-section, whatever the resolution
   sweep says about the other thirteen. Two consequences: the fix is a *joint*
   choice of `RESOLUTION` and `N_MODES`, not resolution alone, and the padding
   behaviour itself is a library defect
   (`fdtd_vs_EME_mode_validation.md` **item 14**, which also covers the hang one
   band further up).

   Per check 3, the remaining question is which quantity *this device's*
   result depends on — the coupler's lesson was that sweeping an insensitive
   observable hid a 7% error in a sensitive one. A mode spectrum is not that
   observable; it is the precondition for choosing one. So this gap is
   measured, not closed.

3. **No full-device FDTD reference.** Needs a `modal_multiplexer_fdtd.py` written
   before check 4 is possible. Cost is the concern — a multiplexer is long, and a
   single-domain 3-D run at res 40 over that span is expensive (the coupler's
   78 µm device took ~3 h on 48 ranks). Worth scoping at reduced resolution
   first: even a res-20 reference that only resolves the qualitative
   mode-sorting would catch a cascade that is mixing modes, which is the failure
   mode gap 1 is about.

   **⚠️ BLOCKED 2026-07-31 — this is a real library gap, not a cost problem.**
   `modal_multiplexer_fdtd.py` was written (same `CompoundGeometry` the hybrid
   script builds, one single-domain `FDTDSolver.solve()` per input mode instead
   of five cascaded sections, res 20 default per the note above). It cannot run:
   `FDTDSolver.solve()` (the full single-domain path) is unconditionally
   incompatible with `CompoundGeometry` as the top-level device.
   `solve() → to_meep_prisms() → geometry.polygons(stub_length=...)`, and
   `CompoundGeometry.polygons(stub_length>0)` raises `NotImplementedError` by
   design — "no generic way to know per-port stub width for an arbitrary
   segment type" (its own docstring). `HybridSolver` never hits this because
   it only ever calls `solve_section()`, which uses `section_gds()`/
   `to_gds_eme()` instead — this device's own `modal_multiplexer_hybrid.py`
   has run fine (job 11585847) the whole time.

   Confirmed by three cheap failed submissions (11589771/11589899/11589936,
   ~15 s each — see `running_jobs.md`), not by inspection alone. One real,
   smaller bug found and fixed along the way: `CompoundGeometry.
   source_aperture_at()` used to raise for any x past its own footprint, where
   every leaf geometry (e.g. `MMI.width_at()`) already clamps to its terminal
   value for arbitrary overshoot — `FDTDSolver.solve()`'s own source-position
   query lands exactly in that region for PML clearance, so the *actual
   simulation* would have hit this too, not just the geometry-preview plot
   that caught it first. Fixed and tested (`beef/geometry/compound.py`,
   `tests/test_geometry.py::TestCompoundGeometrySourceApertureClamps`) — real,
   but not sufficient on its own.

   **Not fixed further.** The `polygons()` gap needs a genuine design decision
   (which segment's stub-port width wins at an ambiguous end, e.g. a
   multi-port `WaveguideRun`), not a mechanical patch — the same caution this
   doc's own gap 6 already applied when declining to touch `solve_section()`
   under-verified. See `notes.md` Geometry item 8. `modal_multiplexer_fdtd.py`
   and its batch script are both marked BLOCKED; do not run either until this
   is implemented.

4. **No golden.** Add a `modal_multiplexer_hybrid` target to
   `tests/test_golden_examples.py`'s `TARGETS`. As with the 2×2 MMI, do it after
   the resolution question settles, so the recorded values are not immediately
   invalidated.

5. **Engine agreement not measured.** Adapt
   `verify_item_6_engine_agreement.py` onto a uniform stretch of the body. The
   extraction helper (`tests/test_engine_agreement.py::dbeta_from_phases`) is
   device-independent.

6. **New, opened 2026-07-30 — `FDTDSolver`'s duplicate-eigenmode-band bug.
   DONE 2026-07-30 for `section_mode_neff()`/`section_mode_fields()`;
   `solve_section()` deliberately deferred (see below).**

   Gap 1's field-overlap check found that `sim.get_eigenmode()`, called in a
   `for band in range(1, n_modes+1)` loop, can return the SAME eigenvector
   for two different requested band numbers when they are closely spaced
   (confirmed: self-overlap 0.9999999999 between two "different" bands, on a
   cross-section where every other band pair overlaps ≤ 0.126). This silently
   makes an `N_MODES`-mode basis rank-deficient by one mode per duplicate,
   which is a strictly worse failure than reordering.

   **Fix implemented**: `FDTDSolver._iter_distinct_bands()` (`beef/solvers/
   fdtd.py`) wraps the band loop. It solves a cheap diagnostic-resolution
   eigenmode per `band_num`, computes its transverse-E-field overlap (via
   `beef.mode_identity.field_overlap_matrix`) against every already-accepted
   band in the same call, and if the overlap exceeds `_DUP_BAND_OVERLAP_
   THRESHOLD = 0.9` it logs a warning and moves on to the next `band_num`
   instead of accepting the duplicate — up to `_DUP_BAND_MAX_EXTRA = 8` extra
   attempts before raising `RuntimeError`. Wired into both
   `section_mode_neff()` and `section_mode_fields()`.

   **Sanity check — inert on a known-good device**: re-ran the directional
   coupler's audit before/after the fix. Identical results, no bands skipped,
   runtime unchanged within noise (108.5s vs 111.0s).

   **Confirmation — fixes the real reproduction**: re-ran this device's
   interface-0 field-overlap check (gap 1's Step 2) at production resolution
   (cmp_res=40) with the fix applied.

   *Before* (from gap 1, Step 2): meep neff had exact ties at
   1.747750/1.747750 (bands 9/10) and 1.484797/1.484797 (bands 13/14), with
   catastrophic η ≈ 0.01–0.04 mismatches against meow at those rows.

   *After*: the guard fired twice —
   ```
   FDTDSolver: meep band_num=11 duplicates an already-accepted eigenmode
   (overlap 1.000000 > 0.90 threshold) -- skipping and trying the next band_num.
   FDTDSolver: meep band_num=15 duplicates an already-accepted eigenmode
   (overlap 1.000000 > 0.90 threshold) -- skipping and trying the next band_num.
   ```
   The resulting neff sequence has no exact ties anywhere:
   ```
   [2.706419 2.676302 2.625392 2.552536 2.455959 2.333052 2.180088 1.991978
    1.763693 1.747750 1.683602 1.593585 1.484797 1.433726 1.413238 1.387737]
   ```
   and every guided row's diagonal η is clean (0.973–0.998, all with margin
   ≥ 0.80 except one 0.9852-vs-lower-margin row) — no more ~0.01–0.04
   catastrophic mismatches. The one remaining oddity (`order_preserved:
   False`, meep index 8 → meow index 9) is the same one-index shift already
   identified in gap 1's Step 1 as a genuine, confidently-matched
   correspondence gap (not a solver bug) — the duplicate-band bug is what's
   fixed here, not that separate finding.

   **Unit tests**: 5 new meep-free tests in `tests/test_duplicate_band_guard.py`
   (fakes `sys.modules["meep"]` with a minimal stand-in module) covering the
   no-duplicates path, one-duplicate skip-and-retry (plus its warning log),
   exhaustion → `RuntimeError`, and — the important negative case — that two
   genuinely different mode shapes with near-identical `k_x` are NOT flagged
   as duplicates (this check is field-shape-based, not neff-proximity-based).
   All 5 pass, alongside the existing mode-identity/overlap/portspec suites
   (56 passed, 2 skipped) and a clean full-repo collection (366 tests, no
   import-time regressions).

   **Scope decision — `solve_section()` deliberately NOT changed.** This is
   the production interface-basis builder used by `HybridSolver`, and it was
   confirmed to hit the *same* underlying meep behavior through a different
   API: it calls `sim.get_eigenmode_coefficients()` on a DFT-decayed flux
   after a full `sim.run()`, not the cheap no-time-stepping `get_eigenmode()`
   that `_iter_distinct_bands()` wraps. Extending the same guard there would
   require either (a) an extra cheap eigensolve per band purely for the
   duplicate check, which is easy to get wrong silently since it's a
   different code path than the one actually producing the basis, or (b) a
   real time-domain verification run to prove it's correct, which this
   session didn't have budget for. Given `solve_section()` is a heavily used,
   delicately tuned, already well-tested production method, changing it
   under-verified was judged too risky.

   **Recommended follow-up — DONE 2026-07-30, see
   `fdtd_vs_EME_mode_validation.md` item 8.** A *post-hoc* diagnostic inside
   `solve_section()` that checks its already-computed `T` coefficient-array
   columns for near-duplicate values, once the flux-decayed run has already
   completed. Costs zero extra simulation time (no new eigensolves), reuses
   data `solve_section()` already produces, and turns today's silent
   rank-deficiency into a loud (logged, not raised) warning. Wired into both
   `FDTDSolver.solve_section()` and `EMESolver.solve_section()`, not yet
   exercised against this device's own real interfaces (needs an actual
   meep/meow run, not a code change).

   **Exercise run submitted 2026-07-30, re-queued 2026-07-31** — SLURM job
   **11585847** (the 07-30 submission never started, blocked on
   `AssocGrpBillingMinutes` under `gts-sr108-ece`; resubmitted under
   `gts-sr108-gedc_main`)
   (`batch/run_modal_multiplexer_examples.batch`), this device's first real
   run of `modal_multiplexer_hybrid.py` in this tree (no
   `example_output/modal_multiplexer/` existed before this). 4 nodes × 24
   ranks, up to 5h — the most expensive single job run this session, per
   "most demanding hybrid device" above. Exercises the item-8 `T`-column
   scan against this device's real four interfaces for the first time;
   watch `reports/11585847_0.out` for `warn_on_transfer_duplicate` firing on
   either solver.

   **Partial reading, 2026-07-31 10:53 (job still RUNNING, 2h12 of 5h) — the
   meow-side guard fires on the production path.** Two of the FDTD subregions
   have completed (`run 0 finished`), and the first EME body solve after them
   logged, once per rank across all 96 ranks:

   ```
   EMESolver.compute_modes: meow mode 9 duplicates already-accepted mode 8
   (self-overlap 1.000000 > 0.90) -- skipping it and reaching further up the
   spectrum.
   ```

   That is the item-0 guard, on `EMESolver`'s production `compute_modes()`, on
   this device's real interfaces at the shipped `RESOLUTION = 40` — the first
   time the meow-side duplicate has been caught outside an audit script.

   **Updated 12:16 (3h35 in): six firings now, from two distinct pairs.**
   576 log lines / 96 ranks = 6 events — **four** of `mode 9 duplicates mode 8`
   and **two** of `mode 11 duplicates mode 10`, both at self-overlap 1.000000.
   The second pair is the one item 0 originally found (meow's modes 10 and 11
   identical at res 40 on this cross-section), so **both** of item 0's pairs are
   now confirmed firing on the production path, repeatedly, across this device's
   interfaces and wavelength points rather than once. Three
   things follow:

   - **It confirms the defect is live in production, not an artifact of the
     audit harness.** Item 0 was verified against a re-run of gap 1's own
     field-overlap check (job 11573880); this is `HybridSolver` doing real work.
     Exactly one distinct event so far (mode 9 → mode 8), replicated across
     ranks, not 96 separate firings.
   - **Had the guard not been there, this run's interface basis would be
     rank-deficient by one of the 16 requested channels** — on the device whose
     entire purpose is to keep 16 channels separate. The guard's own log line
     says as much; this is the first device where that sentence is not
     hypothetical.
   - **`warn_on_transfer_duplicate` has NOT fired** (no `transfer` line in the
     log at all) — so the item-8 post-hoc `T`-column scan, the thing this job
     was submitted to exercise, is still silent. That is the expected outcome
     *if* the upstream guard is doing its job — the scan looks for duplicates
     that survive into `T`, and the duplicate was rejected before it got there.
     It is **not yet** evidence the scan works, since a detector that has never
     fired and a detector that cannot fire look identical. Gap 6's original ask
     (exercise the scan) is therefore still open at the time of writing; the run
     has ~2h45 left and three more interfaces to reach.

## Order

**1** first and by a wide margin — cheap, and it tests the property this device
exists to demonstrate. **Done 2026-07-30** — see gap 1. It surfaced **6**,
which blocked everything below more directly than gap 1's original
mismatched-basis concern did: a rank-deficient interface basis, from a
solver-level bug, would corrupt gaps 2 through 5's results regardless of
resolution or `N_MODES` choice, and is not specific to this device — any
example with closely-spaced higher-order bands (the coupler, both MMIs) could
be silently affected too, just not yet checked there. **6 is now done**
(2026-07-30) for the audit-only entry points, confirmed by re-running this
device's own reproduction with the fix applied — see gap 6. `solve_section()`
itself is deliberately deferred (also gap 6) — its own post-hoc `T`-column
duplicate scan is now implemented (`fdtd_vs_EME_mode_validation.md` item 8)
but not yet exercised against this device's real interfaces, so treat it as a
still-open risk for **2** through **5** below, though a lower-probability one
now that both the mechanism is understood and a detector (not yet proof) is
wired in. ~~Next up: **2**.~~ **2 is now measured** (2026-07-30, job 11575516)
— the EME side is converged at res 60, the FDTD side is not converged at 60 and
needs its sweep extended, and neither is at the shipped 40. Two things follow
for the order below:

- **Do not record a golden (4) at `RESOLUTION = 40`.** That was already the
  stated reason for sequencing 4 after 2; the measurement makes it concrete
  rather than precautionary. Pinning values on a mesh both engines are known to
  be short of would need a third regeneration.
- **Extending 2 to meep res 80/100 is cheap and comes first** — same script,
  same cross-section, one more row. Everything below it is a simulation
  campaign; this is not.

Then **4**, to stop the numbers drifting unnoticed while the rest is worked on.
**3** and **5** last: both are real simulation campaigns.
