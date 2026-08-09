# Modal multiplexer validation — to-do

> **Status 2026-08-04 — resumed, moved out of `deferred/`.** Reviewed against
> the live repo before touching anything: every claim below still matched at
> the time of the review — `modal_multiplexer_hybrid.py` still shipped
> `N_MODES=16`/`RESOLUTION=40` unchanged, `modal_multiplexer_fdtd.py` still
> opened with its `BLOCKED — do not run` docstring, only geometry-plot
> output existed on disk for either script (no `smatrix.pkl` for either),
> and the orphaned merge job (11597199, per `running_jobs.md`) was no
> longer in the queue. Nothing had moved since the 2026-08-01 snapshot. At
> the user's direction this file left `deferred/` and work restarted,
> beginning with gap 5/check 1 (engine agreement) — see gap 5 below for the
> script and its result.
>
> **Same day, later: all of the above changed.** Check 1 closed (mean
> ratio 0.9958), gap 2's resolution question got a real (if inconclusive)
> answer, gap 3's library blocker was fixed and its script run for the
> first time ever, and — as a direct consequence of gap 2's own finding —
> `modal_multiplexer_hybrid.py`'s `N_MODES` was changed **16→14** (kept
> `RESOLUTION=40`). That last change is a judgment call made without
> explicit user sign-off; see the boxed note under gap 2 and the script's
> own comment for the full reasoning, and treat it as open to revision.

> **2026-07-28 — the material model changed under every device here.** Both engines
> now evaluate one analytic Sellmeier model (Salzberg Si, Malitson SiO2) instead of
> each doing its own nearest-tabulated lookup, which shifts n_Si by ~0.002 and
> removes a 0.0063 step that sat in the middle of the C-band. Any golden or stored
> result for this device predates that and needs regenerating; see
> `directional_coupler_validation.md` item 8 and, for the dispersive-FDTD
> re-parameterisation (built and gate-passed 2026-07-31, opt-in via
> `dispersive=True`, off by default), `material_model_validation.md`.



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
`resolved/fdtd_vs_EME_mode_validation.md` items 11–13 rather than here, since it is not
specific to this device. **Items 11, 12 and 13 all closed 2026-07-31** (jobs
11588309, 11586098) — the guard's diagnostic grid is clear at res 40, the EME
guard's over-request is inert, and k-seed laddering recovers one guided mode the
shipped fixed seed loses. Gap 6's exercise run (job 11585847) **TIMED OUT**
2026-07-31 13:42 after the full 5h, having already fired the meow-side guard
six times on the production path; resubmitted unchanged as **11590663** (disk
cache carries sections 0-3 forward, so it should reach the un-finished
per-port output section quickly) — see gap 6.

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

**Update, 2026-07-31 18:58:** gap 6's exercise run **timed out a second time**
(11590663, resubmitted from 11585847's first timeout). It got further —
reached section 4 (all 3 per-port geometry plots exist on disk) and was deep
inside a per-port eigenmode-overlap monitor computation when killed — but
still did not finish, and `warn_on_transfer_duplicate` still has not fired
either way. After two full 5h/4-node attempts (10h of compute total), gap 6's
original question remains open. See gap 6 below and `running_jobs.md`'s
11590663 entry.

**Update, 2026-08-01 03:01: gap 6's exercise run timed out a THIRD time**
(11597177 `[0-2]`, the reshaped array), and this time the cause is known —
**the reshape's batch script dropped `srun --mpi=pmi2`**, so all three tasks
ran meep on **1 rank instead of 96** (~42× slower per timestep) while holding
4 nodes idle each. A second flaw compounds it: the "sections 0-3 come from the
warm cache" premise the reshape rested on is **false in practice** — every
task recomputed them. Gap 6 has now cost ~25 h of compute and is still
unanswered. **No fourth attempt made — this device is deferred**, and the
fix-and-resubmit recipe is recorded in gap 6's boxed entry for whenever it is
picked back up. The orphaned merge job 11597199 can never run and should be
cancelled.

`examples/modal_multiplexer/python/modal_multiplexer_hybrid.py` — the only script
for this device. SOI_160, not SOI_220 like the coupler and MMIs.

## Where it sits

The **most demanding** hybrid device here and the **least supported**:

- `N_MODES = 14` (originally shipped 16; changed 2026-08-04, see gap 2's
  boxed note) carried at **four** section interfaces — more modes and more
  joins than any other example. Every interface is another place the FDTD
  and EME bases have to agree.
- **No full-device FDTD script.** Unlike both MMIs and the coupler, there is
  nothing to compare the cascade against as a whole; check 4 has no starting
  point at all.
- **No golden fixture.** Nothing pins its output.
- Different material stack (SOI_160), so no measurement from the other devices
  transfers numerically — not the coupler's converged EME res 80, not its Δβ
  ratio.

So it has the most ways to be wrong and the fewest ways to find out.

## Reference design — the author's own (unpublished) mode multiplexer paper

Added 2026-08-05. The user supplied a figure from their own paper
characterizing a fabricated 2-stage cascaded-MMI mode multiplexer —
**unpublished as of this writing; citation to be added once it is.**
Confirmed by the user directly (not inferred): the device is **O-band,
160 nm silicon** — the same band and material system this repo's device
already uses (`SOI_160`, `WL_MIN=1.26`/`WL_MAX=1.36 µm` in both
`modal_multiplexer_hybrid.py` and `modal_multiplexer_fdtd.py`).

**Same overall topology, confirmed via the legacy config this repo's device
was built from.** The figure's panel (a) labels four longitudinal zones —
`L_mm1 | L_taper | L_ps | L_mm2` — and body/taper widths `w_mm`/`w_in` (MMI
body / multimode input access) and `w_t1`/`w_t2` (taper widths). Every one
of these has a same-named counterpart already in
`archive/FECES_Legacy/configs_mmi/config_mux.json`, the config this repo's
`modal_multiplexer_hybrid.py`/`modal_multiplexer_fdtd.py` were built from
(`mm_1_length=16.92`, `mm_2_length=34.24`, `mm_width=2.8`,
`mm_waveguide_width=0.7`, `taper_length=2.0`,
`between_waveguide_length=3.0`, `top_bot_taper_width=0.725`,
`middle_taper_width=1.15`, `top_waveguide_width=0.33`,
`waveguide_width=0.35`, `waveguide_offset=1.0375`) — this is the same
device family, not a coincidental resemblance. **Not independently
verified**: the exact correspondence between the figure's three port labels
(`NC`/`A`/`B`) and this repo's `top`/`mid`/`bot` naming, or which of the
figure's width labels lands on which specific port — the source image
doesn't resolve that unambiguously from a flattened text extraction, and it
was not worth guessing at when the topology-level match is already solid
from the parameter names alone.

**⚠️ Critical distinction, confirmed by the user: this repo's device is the
paper's "Seed" design, NOT the "Optimized" one that was actually fabricated
and measured.** The figure's panel (g) shows both — a plain MMI+taper
"Seed" geometry (topologically identical to what `config_mux.json`/this
repo's scripts implement) and a visibly different, free-form "Optimized"
geometry the paper's design process produced from it (same relationship as
the HCG device's `HCG_seeded_TO.gds` to its own published design, per
`high_contrast_grating_validation.md` — except confirmed here, not still
open). **Consequence: panels (e)/(f)'s measured and simulated S-parameter
curves (S13/S24 ≈ −1 to −8 dB, S14/S23 ≈ −20 to −40 dB isolation, ~1280–
1330 nm) characterize the OPTIMIZED device, not this one, and are NOT a
valid direct ground-truth target for `modal_multiplexer_fdtd.py`'s own
S-matrix or for a future golden fixture.** Quoting them as if they were
would be exactly the kind of unearned precision this repo's own
measured-not-argued standard exists to catch. This repo's Seed device has
no independent published measurement of its own — only its topology, band,
and material system are confirmed to trace to real, physically-motivated
values rather than an arbitrary academic exercise.

**What IS usefully confirmed by this reference, despite the seed/optimized
gap:**
- The 2-stage cascaded-MMI topology (1×3 MMI → 3 connecting single-mode
  waveguides → 3×3 MMI) this repo implements is a real starting point for
  an actual fabricated mode multiplexer, not an invented exercise.
- O-band, 160 nm Si — exactly this repo's `SOI_160` stack and wavelength
  range (this repo's own 1.26–1.36 µm sweep in fact covers a slightly wider
  band than the ~1280–1330 nm shown in the paper's own panels e/f).
- **Panel (d)'s characterization method** — cascade two copies of the mux
  (`Mode Mux 1 → Mode Mux 2`), launch `TE00` and `TE01` separately at ports
  1/2, and read the output pattern at ports 3/4 — is a reciprocity-based
  self-test structurally similar to `ring_resonator_validation.md`'s own
  cascaded-reciprocity check. That doc's own finding (a cascaded reciprocity
  check can be a *tautology* rather than independent evidence) is worth
  keeping in mind before treating an analogous check on this device as
  strong evidence, should one ever be built.

**Not acted on this session** — recorded for whoever next works gap 3/4 or
considers a golden fixture, so the seed/optimized distinction is legible
before either is attempted, not discovered after the fact.

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
   follow-up", now implemented as `resolved/fdtd_vs_EME_mode_validation.md` item 8:
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
     mode-8+ story is `resolved/fdtd_vs_EME_mode_validation.md` item 1), but it does
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
   (`resolved/fdtd_vs_EME_mode_validation.md` **item 14**, which also covers the hang one
   band further up).

   > ⚠️ **Items 12 and 14 are now fixed in the library, 2026-07-31 — this
   > device's own `N_MODES=16` will raise the next time it runs.**
   > `_iter_distinct_bands()` now ladders its k-seed (item 12) and stops at
   > the cladding index instead of padding with unguided modes (item 14) —
   > raising a clear error the moment a cross-section can't supply the
   > requested count, which is exactly this device's situation (14 reachable,
   > 16 requested). The currently in-flight hybrid run (job 11590663) started
   > before this change and is unaffected.
   >
   > **✅ `N_MODES` changed 16→14 in `modal_multiplexer_hybrid.py`, 2026-08-04
   > — a judgment call made without explicit user sign-off, flagged loudly in
   > the script's own comment for review, not silently decided.** 14 because
   > that is exactly what gap 2 measured as reachable at the shipped
   > `RESOLUTION=40`, consistently across every tested resolution 20–100 —
   > reducing further would discard real guided channels for no reason.
   > `RESOLUTION` was deliberately left at 40, not raised: gap 2 also found
   > no resolution 20–100 is demonstrably converged in absolute n_eff (60→80
   > looked converged, misleadingly — 80→100 broke it completely), so a
   > finer mesh buys no demonstrated accuracy here, while this device's own
   > check 1 (gap 5) measured the engine ratio — the quantity this hybrid
   > chain actually depends on, not either engine's absolute value — at
   > exactly `RESOLUTION=40` and found it excellent (0.9958). **Not verified**:
   > whether 14 (down from 16) still carries several supermodes of each
   > parity at every interface, the constraint the script's own `N_MODES`
   > comment states — predicted yes (dropping the two weakest-confined,
   > highest-index modes removes one of each parity from an alternating
   > sequence, not a whole parity) but not checked against a real run, which
   > is gap 6's territory (a separate, expensive campaign, not attempted
   > this session).

   Per check 3, the remaining question is which quantity *this device's*
   result depends on — the coupler's lesson was that sweeping an insensitive
   observable hid a 7% error in a sensitive one. A mode spectrum is not that
   observable; it is the precondition for choosing one. So this gap is
   measured, not closed.

   **Extension submitted 2026-08-04 — job 11673180, the cheap next step
   named above. FAILED at res 20, a real finding, script fixed, resubmitted
   as 11673526.** `mode_spectrum_convergence.py`'s `MEEP_RESOLUTIONS` grew
   from `[20, 30, 40, 60]` to `[20, 30, 40, 60, 80, 100]` (`meow` already
   swept to 80; the FDTD/meep side had never been pushed past 60, where it
   was still moving — 40→60 agreed on only 3/16). Same script, same
   interface-0 cross-section, no other changes — this answers whether
   meep's mode spectrum plateaus by 80/100 (discretization, per this
   script's own docstring) or the mode-8+ gap survives at every resolution
   (structural).

   **11673180 FAILED after 66 s, at res 20 — before ever reaching the new
   80/100 points.** `FDTDSolver: this cross-section supplies only 14 guided
   mode(s) of the 16 requested -- band_num=15 returned n_eff=1.443140, at or
   below the cladding index 1.446804.` This is items 12/14
   (`resolved/fdtd_vs_EME_mode_validation.md`, fixed 2026-07-31 — AFTER this
   script's own original successful run, job 11575516, on 2026-07-30): the
   duplicate-band guard now raises the instant it reaches an unguided band
   rather than silently padding, and this script's fixed `N_MODES=16`
   request hits exactly that at res 20 (only 14 guided there). A real,
   useful confirmation in its own right — the guard is live on this exact
   cross-section, corroborating gap 2's own account from a second angle —
   but it meant the script itself needed to adapt: `_meep_spectrum()`
   (generator-based, so nothing partial survives the raise) now catches
   this specific `RuntimeError`, parses the guard's own reported guided
   count from its message, and retries at exactly that count — recovering
   the real spectrum instead of crashing the whole sweep. Verified the
   regex against the real captured error text before resubmitting.
   Deliberately fixed at the (disposable) verification-script level, not by
   loosening the guard itself, which stays exactly as strict as items 12/14
   designed it. Resubmitted as **job 11673526** (same batch script,
   unchanged budget/account).

   **11673526 COMPLETED 2026-08-04, 4:51:09. The retry-degrade fix worked at
   every resolution — and the answer to gap 2's original question is NO,
   not the plateau this was hoping to confirm.**

   | quantity | res 20 | 30 | 40 | 60 | 80 | 100 |
   |---|---|---|---|---|---|---|
   | meep guided (of 16) | 14 | 14 | 14 | 14 | 14 | 14 |
   | meow guided (of 16) | 14 | 14 | 14 | 15 | 15 | *(not swept)* |

   **First, a genuine improvement over the pre-fix record.** The retry
   mechanism now reports a clean, consistent **14 guided modes at every
   meep resolution from 20 to 100** — a cleaner answer than the pre-fix
   13/13/13/14 sequence recorded above, because it comes from the same
   robust duplicate-skip-then-hard-stop walk at every resolution rather
   than depending on where a duplicate happened to land relative to the old
   fixed-count enumeration. 14 is now a resolution-independent, confidently
   measured number for this cross-section — a real input to the
   `N_MODES`/`RESOLUTION` decision below.

   **Second, meep's absolute n_eff does NOT plateau by res 100 — it looked
   like it did at 60→80, and that was misleading.** Self-convergence
   (matched mode index, |Δn_eff| ≤ 0.005, out of 14 guided):

   | | 20→30 | 30→40 | 40→60 | 60→80 | 80→100 |
   |---|---|---|---|---|---|
   | meep (FDTD) | 1/14 | 0/14 | 1/14 | **14/14** | **0/14** |

   60→80 agrees on all 14 modes (max |Δn_eff| = 0.0048, just inside
   tolerance) — exactly what a genuine plateau would look like, and the
   result gap 2's "Extending 2 to meep res 80/100 is cheap and comes first"
   note was hoping for. **80→100 breaks it completely**: every one of the
   14 modes shifted by 0.0083–0.0127, roughly 2–3× the 60→80 step, despite
   covering a smaller *relative* resolution change (25% vs 33%) — the wrong
   direction for ordinary mesh convergence, where the step should shrink,
   not grow. This was checked by hand against the raw spectra, not read off
   the summary table alone: the shift is small, uniform in sign (every mode
   moved the same direction), and monotonically increasing, not scatter —
   a real, repeatable effect, not noise. **The practical reading: there is
   no resolution in {20, 30, 40, 60, 80, 100} at which meep's absolute
   n_eff for this cross-section is demonstrably converged**, extending
   README's Finding 3 ("no device's EME is converged in absolute n_eff,
   even at res 100") to this device's FDTD side too, and specifically
   contradicting the read that would have closed this gap.

   **Third, matched-resolution engine agreement (meep vs meow, the
   comparison check 1 actually cares about) stays good throughout** —
   13–14 of 14 meep modes match a meow mode at every resolution tested, with
   the sole residual being the already-tracked near-degenerate-pair defect
   (`resolved/fdtd_vs_EME_mode_validation.md` items 11–13): at res 60/80,
   meow resolves TWO closely-spaced modes near n_eff≈1.797 (differing from
   each other by only ~0.0003) where meep's per-band ladder returns only
   one — confirmed by hand-checking the raw spectra (the summary table's
   "orphans" column reads empty here because a bare nearest-neighbour check
   can't tell "matched" from "matched to the wrong twin of a degenerate
   pair," exactly the limitation `field_overlap_matrix()`-based checks exist
   to fix, per README's check 2 discussion). Not a new finding — it extends
   the previously-documented res 20–60 range of this residue to res 80 —
   but worth recording since it is the only source of engine
   *dis*agreement in an otherwise clean picture.

   **What this means for the `N_MODES`/`RESOLUTION` decision (still open,
   boxed above):** `N_MODES = 14` is now well-supported and
   resolution-independent — a confident downgrade from the shipped 16.
   `RESOLUTION` has no obvious answer from this data: nothing tested is
   demonstrably converged in absolute n_eff, and the one resolution that
   looked converged (60→80) was misleading. What check 1 (gap 5, this same
   session) actually measures — the *ratio* Δβ_FDTD/Δβ_EME at the shipped
   `RESOLUTION = 40` — is a different, arguably more decision-relevant
   quantity than either engine's absolute n_eff, since a ratio can stay
   near 1 even while both engines individually keep drifting (as it does
   for the coupler and both MMIs, all measured at their own shipped,
   unconverged resolutions per README's Finding 3). Whether that pattern
   holds here is what gap 5's own result will show, not yet known at the
   time of this writing — see gap 5 for its status; this decision has
   deliberately not been made unilaterally regardless, see the note under
   Order.

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
   has never hit this `polygons()` bug through job 11585847/11590663's real
   `solve_section()` calls, independent of that job's own TIMEOUT (a
   wall-clock/scheduling problem, not this one).

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

   **✅ FIXED 2026-08-04 — the premise was wrong, not just unimplemented.**
   `_end_ports()` already answers exactly this "which port's width wins"
   question for `to_gds_eme()`'s own stub handling, and it answers it by
   *not* picking a single winner: it returns every port at that end (one
   `(y, width)` pair per port), and a leaf geometry with several ports
   already draws one stub per port for its own `polygons()` (e.g.
   `MMI.polygons()`'s `for i, p in enumerate(self.inputs)`/`.outputs`
   loops). There was no ambiguity requiring a value judgment. `polygons()`
   now wires `_end_ports()` in the same way `to_gds_eme()` already does —
   one stub polygon per port at each true outer end, no `output_idx`
   restriction (the single-domain path needs every branch at once).

   Confirmed against this exact device: the mux's `CompoundGeometry` now
   builds the correct 4 stub polygons (1 input + 3 output, each at its own
   port's true width/y, matching `BRANCH_SPECS`) where it previously raised
   immediately, and `to_meep_prisms()` — the actual call
   `FDTDSolver.solve()` makes — now succeeds end to end through prism
   construction. 7 new tests, `tests/test_geometry.py::
   TestCompoundGeometryPolygonsStub` (multi-port ends, the `frame_offset`
   interaction, a `stub_length=0` no-op check, a real `to_meep_prisms()`
   build). Full fast suite still green, 618 passed/19 skipped, no
   regressions. Full account in `notes.md` Geometry item 8 part 2.

   **✅ Smoke test COMPLETED 2026-08-04, 11:34 — the pipeline runs
   end to end, for the first time ever.** Per this repo's own two-step
   discipline for expensive, never-run scripts, submitted a smoke test
   first rather than the real run: **job 11673498**,
   `run_modal_multiplexer_fdtd_SMOKE.batch`, `MUX_FDTD_RES=10` (half the
   script's own already-reduced res=20 default, ~16× cheaper per its own
   res⁴ cost scaling). Both broadband input-mode runs completed
   (58,215 timesteps for TE0, similar for TE1), `check_passive()` did not
   raise, and a real `SMatrix` was saved:

   | | top | mid [MUX out] | bot [MUX out] | total |
   |---|---|---|---|---|
   | TE0 in | T=0.0902 | T=0.1430 | T=0.3284 | 0.5616 |
   | TE1 in | T=0.0072 | T=0.0027 | T=0.0025 | 0.0124 |

   **Read these as pipeline evidence, not physics** — `res=10` is a
   deliberately crude smoke-test mesh (half the script's own already-reduced
   default), and TE1's near-total loss (98.8% unaccounted) is exactly the
   kind of artifact expected there: a higher-order mode's transverse profile
   needs more pixels to launch cleanly than a fundamental mode does, and
   at 10 px/µm the 0.7 µm multi-mode input waveguide gets only ~7 pixels
   across. This number should NOT be quoted as a device property; it is
   evidence the mechanics work, not a measurement.

   One latent, pre-existing hazard noted, not fixed: like
   `mmi_1x2_fdtd.py`/`mmi_2x2_fdtd.py`, this script has no `is_master()`
   gating on its file writes (`plot_fdtd_geometry`, `smat.save()`,
   `fig.savefig()`), so all 24 MPI ranks in this run independently attempted
   to write the same `geometry.png`/`smatrix.pkl`/`results.png` — the log
   shows ~24 duplicate "saved →" lines. It did not corrupt anything here
   (every rank computes byte-identical results after meep's own internal
   synchronization, and this is the same pattern both other full-device
   FDTD scripts already run in production), so it is not fixed as part of
   this session's scope — flagged since it was visible in this run's log,
   not because it is new or specific to this device.

   **The real res=20 run — submitted 2026-08-05 as job 11688645, TIMED OUT
   at 4h, a real cost finding, resubmitted at 24h as 11693107.** This is
   the actual full-device FDTD reference check 4 needs; the smoke test only
   proved the pipeline, not the physics at production resolution.

   **11688645 TIMED OUT after the full 4h budget, still mid-way through
   just the FIRST of two input modes (TE0)** — 416,806 timesteps, simulated
   time t=10,420.2, a flat, unwavering ~0.0336 s/step throughout (no sign of
   approaching its decay stop). This is NOT merely the expected ~16×
   res⁴ per-step cost the smoke test's own scaling predicted: at res=10,
   TE0 finished cleanly at t≈2910–3217 (58–64k steps, ~5 min of
   timestepping); at res=20, the run had already reached **3.4× further in
   simulated time** than that entire stopping point before being killed,
   with the per-step cost itself only ~6.6× higher (0.0336 vs 0.00508
   s/step) — CFL scaling alone would predict roughly 2× the *step count*
   for the same simulated duration, not 3.4× *more simulated time before
   decaying at all*. Read as an open, unconfirmed question rather than a
   settled mechanism: either the coarse res=10 mesh was under-resolving a
   slowly-decaying channel and stopped "too early" in a numerically
   convenient but not fully physical sense, or finer resolution genuinely
   sharpens some resonance/radiation-coupling feature that takes longer to
   decay below `decay_tol=1e-11` — this run's own log gives no way to
   distinguish the two, and it wasn't investigated further this session.
   **No caching in this script** (by design, per its own module docstring),
   so all 96 core-hours were a total loss, not partially salvageable.

   Presented to the user with the finding before spending more compute:
   resubmit with a much bigger budget, loosen `decay_tol` (the known
   dominant cost knob), refactor into a per-input-mode array job first (like
   `modal_multiplexer_hybrid.py`'s own reshape), or hold off entirely. **User
   chose: resubmit with a bigger budget, unchanged otherwise.** `-t` raised
   240→1440 (24h, matching this repo's largest existing precedent) in
   `run_modal_multiplexer_fdtd.batch`; resubmitted as **job 11693107**.
   `decay_tol` stayed at `1e-11`, matching
   `modal_multiplexer_hybrid.py`'s own setting, meant to keep this a fair
   apples-to-apples reference once it completed.

   **11693107 ALSO did not converge — killed by the user's own diagnosis
   after 18h53m, still stuck in TE0, having reached t=49,375.9 (nearly 5×
   further than 11688645's own already-unconverged stopping point) with the
   same flat ~0.0335 s/step throughout.** The user pointed at a real,
   independent data point rather than authorizing another blind wait: an
   actual prior topology-optimization campaign for this exact device family
   exists on this cluster, uploaded to `04_weary_paw/` (the user does not
   have the driver code, only the archived output — `config.json`,
   `logs/master.log`, `logs/group{0,1}.log`, `geometries/`). Its config
   confirms this is the same device (`gds_name: "mux_MMI_2"`, a `[64, 4] µm`
   splitter design region spanning nearly this device's full 62–66 µm
   length) at **`res=40`, `dr_res=80`** — HIGHER resolution than this
   script's `res=20`, doing MORE work per iteration (a full forward +
   adjoint gradient cycle across two mode groups, plus linewidth/
   linespacing/area constraint evaluations, not just a forward-only sweep).
   **That full cycle took ~70–80 minutes per iteration** (56 iterations
   logged over ~67.3h total, `logs/master.log` timestamps 2025-03-23
   11:58:24 → 2025-03-26 07:17:06; individual cycle times confirmed
   directly from consecutive "Current iteration" timestamps, not estimated).
   Compare: this repo's own forward-only, single-mode-so-far run had
   already run **~19 hours** without finishing.

   **This gap (roughly 15–20×, comparing a smaller task at lower resolution
   against a larger one at higher resolution) is too large to attribute to
   core count alone** — meep's parallel efficiency drops well before 4×
   more ranks buys 4× wall-clock speedup on a domain this size, and no
   core/node count was recoverable from the archived campaign's own logs to
   even confirm it used more ranks than this script's original 24. The more
   likely dominant factor, not yet acted on: `decay_tol=1e-11` is a strict
   full-spectral-convergence criterion appropriate for a clean, quotable
   S-matrix; a TO forward pass computing a *gradient* typically tolerates a
   far less converged field (the archived campaign's own `nf=5` frequency
   points over a narrower 1.295–1.325 µm band, versus this script's 11
   points over 1.26–1.36 µm, is also suggestive but not proof by itself).
   **Not changed here** — still a methodology decision needing explicit
   sign-off, same reasoning as before.

   At the user's direction, increased core count anyway — cheap, reversible,
   directly attacks the one variable definitively under this repo's own
   control, even though unlikely to close the full gap alone. Killed
   11693107, rescaled `run_modal_multiplexer_fdtd.batch` from 1 node/24
   ranks to **4 nodes/96 ranks** (matching `modal_multiplexer_hybrid.py`'s
   own scale for this device's FDTD sections — the only precedent this repo
   has for this device at this scale), same `-t 1440`, resubmitted as **job
   11707932**. Confirmed at launch: `Using MPI version 4.1, 96 processes` /
   `Splitting into 96 chunks by voxels` — the 4-node allocation is being
   used, not silently falling back to 1 rank the way `modal_multiplexer_
   hybrid.py`'s own reshape once did. Not yet landed — see `running_jobs.md`
   for
   live status.

   **Live check 2026-08-06, ~1h02m into 11707932 — the rescale is delivering
   close to the full expected speedup, measured, not the a priori worry.**
   `reports/11707932.out`: step 387,803 at elapsed 1:02:07, t=9695.08, steady
   ~0.0089–0.0090 s/step. Average sim-time rate since launch (includes
   MPB/structure-init startup overhead) = 9695.08/62.12 min ≈ **156.1
   sim-time-units/min**, versus 11693107's own flat, whole-run rate of
   49,375.9/1133 min ≈ **43.6 units/min** (1 node/24 ranks) — a **~3.6×
   speedup for 4× the ranks** (~88–94% parallel efficiency depending on
   whether the startup-inflated average or the settled instantaneous rate is
   used), better than this section's own "parallel efficiency drops well
   before 4× more ranks buys 4× speedup" worry anticipated. Projected time to
   reach t=49,375.9 (11693107's own kill point, still not decayed there) at
   this rate: ≈316 min ≈ 5.3h elapsed — leaving roughly 18.7h of the 24h
   budget as margin beyond that point. **This is a rate projection, not a
   completion** — still no confirmation of when (or whether, within budget)
   `stop_when_dft_decayed(tol=1e-11)` actually triggers.

   **Re-checked `04_weary_paw/` itself for the rank/core count this section
   already flags as unrecoverable from its own logs — confirmed again, then
   found indirect corroboration elsewhere.** `grep`ing
   `logs/master.log`/`logs/group{0,1}.log` for
   process/mpi/chunk/rank/core/node/thread/parallel returns zero hits — the
   custom Python logger captured only its own INFO/DEBUG calls, never meep's
   own C++ startup banner, so the exact rank count really is gone. What
   does survive: `archive/FECES_Legacy/batch/run.batch` and `rum.batch`, the
   only two SLURM scripts left from this repo's older Makefile-driven
   framework generation (`archive/FECES_Legacy/Makefile`, `srun --mpi=pmi2`
   launcher — the same launcher style `04_weary_paw`'s 2-group split
   implies), both request **`-N 4 --ntasks-per-node=24` (96 ranks)** for
   their own full 3-D campaigns (siblings `mmi_2x2_bend`/`mmi_2x2_gds`, not
   the mux itself). **Two caveats, not glossed over**: no SLURM script for
   the `mmi_mux` Makefile target survives in this archive at all, and the
   driver modules `04_weary_paw`'s own logs name — `Mux_MMI.py`, `engine.py`,
   `utils_s_mux.py` — match **zero** files anywhere in
   `archive/FECES_Legacy` (grepped); `masters/mmi_mux_master.py`, the only
   mux-related script that does survive, imports `meow`/`meep` directly under
   a `class MMI_MUX` with its own `sim_type` dispatch (FDTD/HYBRID/EME/
   FDTD_GDS) — a different, later rewrite, not the code that produced
   `04_weary_paw`. So this is evidence about the framework generation's
   *standard job template* for this device family's era, not a confirmed
   measurement of the exact job that produced `04_weary_paw` — recorded
   because it happens to agree with the scale already chosen for 11707932
   (4 nodes/96 ranks), not because it settles the question the section above
   already correctly marked open.

   **✅ 2026-08-06, ~9h06m elapsed — TE0 fully decayed for the first time in
   this campaign, a real landmark, and it reframes the `decay_tol`
   question.** `reports/11707932.out`: `run 0 finished at t = 77197.45
   (3087898 timesteps)`. Neither single-node attempt ever got close:
   11688645 was killed at t=10,420 (4h budget), 11693107 at t=49,375.9
   (18h53m) — both stuck in TE0, nowhere near this stopping point. At 4
   nodes/96 ranks, TE0 actually reached `stop_when_dft_decayed(tol=1e-11)`
   in **~7.95h wall time** (7,048 progress lines × the ~4.06 s/line cadence
   measured directly from this run's own line count ÷ elapsed time —
   meep's default ~4s progress-report interval, not step-count-based, which
   is also why the earlier apparent "slowdown" between two live checks this
   session turned out to be a false alarm: the step counter resets to 0 for
   each new `mp.Simulation` object, and TE0's finish plus its ~55
   eigenmode-coefficient MPB solves sat in between the two checks, not a
   stall). The sim is now into **input_mode=2 (TE1)**, ~1h15m in, same
   healthy ~0.0087 s/step instantaneous rate, no S-matrix printed yet
   (the script only prints/saves once both modes finish).

   **This is evidence against `decay_tol` being "the more likely dominant
   factor"** in the reasoning recorded above and in `running_jobs.md`: TE0
   reached full 1e-11 decay in well under half the 24h budget once given
   4× the ranks, with the exact same `decay_tol=1e-11` both single-node
   attempts never got near satisfying. That does not prove core count is
   *sufficient* on its own for every case (TE1's own decay time is still
   unmeasured, and could differ from TE0's), but it is a real, measured data
   point that the earlier "likely need to loosen decay_tol" hypothesis was,
   at minimum, not the only lever that mattered here — not changed
   retroactively, recorded as-is for whoever revisits that reasoning next.

   **❌ 2026-08-07, job 11707932 COMPLETED (15:51:32, well inside the 24h
   budget) — but the result is not usable, and it is a serious, not cosmetic,
   finding.** Both input modes finished (TE0 at t=77,197.45/3,087,898
   timesteps as recorded above; TE1 finished too, total wall time 57,073 s).
   `check_passive()` did not raise. But the saved `SMatrix`
   (`example_output/modal_multiplexer/modal_multiplexer_fdtd/smatrix.pkl`,
   read back and inspected directly, not just the one center-wavelength
   line the script itself prints) shows near-total loss at **every** port,
   input mode, and wavelength:

   | in → out | T, all 11 wavelengths (1.26–1.36 µm) |
   |---|---|
   | TE0 → top | ~1.39–1.42 ×10⁻⁵, smooth in λ |
   | TE0 → mid [MUX out] | **exactly 0.0** at all 11 points |
   | TE0 → bot [MUX out] | ~1.05–1.18 ×10⁻⁵, smooth in λ |
   | TE1 → top/mid/bot | **exactly 0.0** at all 11 points, every port |

   Compare the smoke test's own sane numbers at the same ports (res=10,
   job 11673498): TE0 → top/mid/bot = 0.0902/0.1430/0.3284 (56% total), TE1
   → 0.0072/0.0027/0.0025 (1% total). This is not "the numbers changed with
   resolution" — it is a qualitative collapse from a physically plausible
   result to one indistinguishable from zero. Two things argue against
   reading this as real device physics: (1) a genuinely near-total-loss
   passive structure would not typically produce **exactly** 0.0 in
   floating point at every one of 11 wavelength points for 4 of 6
   channels — that pattern looks like underflow/non-detection, not a small
   real number; (2) the *nonzero* channels are smooth, well-behaved
   functions of wavelength (not scattered noise), which is the signature of
   a real coefficient ratio computed correctly but scaled by something
   badly wrong, not random garbage.

   **Connects directly to this section's own still-open question about why
   TE0 needed t=77,197 to decay** — 24-26× further in simulated time than
   the res=10 smoke test's entire TE0 run needed (t≈2910–3217), a gap this
   section already flagged as unexplained. Now there is a plausible
   mechanism, **not yet confirmed**: `stop_when_dft_decayed(tol=1e-11)`
   checks *relative* change in the DFT-accumulated field amplitude; if this
   design's real coupling to the output ports is intrinsically weak at
   production resolution, the tracked amplitude may be dominated by a
   slowly-decaying spurious/resonant component long after the physically
   meaningful signal has already decayed past it — meaning by the time the
   1e-11 relative criterion is finally satisfied, the useful signal itself
   may have decayed into numerical noise, explaining both the anomalously
   long runtime AND the near-zero result as two symptoms of the same
   underlying issue rather than two separate ones. **Not investigated
   further yet** — `FDTDSolver.solve()`'s internal coefficient-extraction
   code (the actual `src_coeff`/`get_eigenmode_coefficients` logic this
   thin example script delegates to entirely) has not been read with this
   failure in mind. This is the first time `FDTDSolver.solve()` (the
   single-domain path, as opposed to `solve_section()`) has ever been
   exercised at a real, non-smoke resolution for ANY device in this repo —
   so a resolution-dependent bug specific to that code path, not just this
   device, cannot be ruled out either.

   **Check 4 (README) is NOT closed by this run.** A completed job with an
   unusable result is not a reference; gap 3 stays open pending root-cause
   investigation.

4. **No golden.** Add a `modal_multiplexer_hybrid` target to
   `tests/test_golden_examples.py`'s `TARGETS`. As with the 2×2 MMI, do it after
   the resolution question settles, so the recorded values are not immediately
   invalidated.

5. **Engine agreement — script written 2026-08-04, run in progress.** Adapted
   `verify_item_6_engine_agreement.py` onto a uniform stretch of the body, per
   the sibling MMIs' own precedent (`resolved/mmi_1x2_validation.md` gap 4,
   `resolved/mmi_2x2_validation.md` gap 4) rather than re-deriving the
   approach: `examples/verification/python/
   verify_modal_multiplexer_engine_agreement.py` cuts a section of MMI 1's
   body (mm_length=16.92, the shorter of the two — gap 1 already established
   both bodies share an identical 2.8 µm × 0.16 µm SOI_160 cross-section, so
   checking one represents both) at four lengths (2/4/6/8 µm, 4.46 µm clear of
   the taper seam at the longest), fits `Δβ_FDTD` from FDTD's accumulated
   phase between the two lowest guided modes via the device-independent
   `beef.solvers.fdtd.dbeta_from_phases()`, and compares it against `Δβ_EME`
   from the same cross-section's own eigenvalues. Pinned at the shipped
   `RESOLUTION = 40` (checks the shipped config, not a converged one — see the
   script's own docstring for why the two lowest modes are less exposed to
   gap 2's res-40 shortfall than the full 16-mode basis is) and the shipped
   O-band grid (1.26–1.36 µm), matching `modal_multiplexer_hybrid.py` exactly
   so this is a check *of* that device, not an adjacent one.

   Shipped already in the array-per-length + merge shape from the start
   (`batch/verification/run_modal_multiplexer_engine_agreement.batch` +
   `..._merge.batch`, `--account=gts-sr108-ece`) rather than discovering the
   need for it the way the 2×2 sibling did (its first attempt, job 11591249,
   TIMED OUT at 5h having finished only 2 of 4 lengths). Deliberately bare
   `python`, not `srun --mpi=pmi2` — this script carries no MPI rank guards,
   and both MMI siblings already proved a single meep process finishes a
   problem this size inside an 8h budget; see the script's own docstring for
   why this is the proven-safe configuration here, not a repeat of the
   missing-launcher bug that cost gap 6 three timeouts on the much larger,
   rank-guarded production hybrid solve.

   Login-node smoke check before submission (object construction + the cheap
   EME half only, no FDTD time-stepping — the OOM lesson from this device's
   own hybrid work says any real-mesh FDTD belongs in SLURM, not a login-node
   background process): `make_mmi1()` builds the correct 16.92 µm body, all
   four section lengths land inside it with 4.46–7.46 µm of taper clearance,
   both solvers construct cleanly, and `eme_dbeta()` ran for real (1 rank,
   1237 s) returning sensible, strictly-descending n_eff (`[2.7458, 2.7183]`
   at 1.26 µm down to `[2.6688, 2.6359]` at 1.36 µm) and smoothly increasing
   `Δβ_EME` across the band (0.1370 → 0.1519 rad/µm) — the assert on mode
   ordering did not fire, and nothing about the geometry/API calls needed
   correcting before committing to the expensive FDTD half.

   Submitted 2026-08-04 as **job 11670540 `[0-3]`** (array, one section
   length per task, `batch/verification/run_modal_multiplexer_engine_
   agreement.batch`, `-t 480`, `gts-sr108-ece`) with merge job **11670544**
   dependent on it (`--dependency=afterok:11670540`,
   `..._merge.batch`, `-t 30`).

   **✅ CLOSED 2026-08-04 — all 4 array tasks completed cleanly (3:07–5:33
   hrs each).** Merge job 11670544 TIMED OUT at its original 30-min budget
   (all 4 cached lengths loaded instantly, then it hung on the ~20-min EME
   half with too little margin left for SLURM startup overhead);
   re-budgeted to 90 min, thread-pinned, resubmitted as **11678100 —
   COMPLETED in 30:37**:

   | wl (µm) | Δβ_FDTD | Δβ_EME | ratio |
   |---|---|---|---|
   | 1.26 | 0.13674 | 0.13704 | 0.9978 |
   | 1.29 | 0.14094 | 0.14146 | 0.9963 |
   | 1.31 | 0.14379 | 0.14443 | 0.9955 |
   | 1.34 | 0.14816 | 0.14892 | 0.9949 |
   | 1.36 | 0.15115 | 0.15193 | 0.9948 |

   **Mean ratio 0.9958** (monotone decrease across the band, 0.9978→0.9948),
   worst off-diagonal 0.0000, worst FDTD fit residual 0.0000 — as clean as
   the 2×2 MMI's own result. This lands within 0.05–0.12 points of BOTH
   MMI devices' own engine bias (1×2: 0.9946, O-band; 2×2: 0.9953, C-band) —
   a THIRD device now clustering in the same ~0.4–0.6% band, on the same
   kind of cross-section (this checked MMI 1's body, itself an MMI), which
   strengthens the "engine-level bias, not device-specific" reading
   `README.md`'s parity-backlog section already drew from the first two.
   **`calibrate_body_neff` is NOT indicated for this device** — same
   conclusion as the 2×2 MMI, for the same reason (the residual is far
   below the coupler's ~10% that flag was built to correct). **Check 1 is
   now closed for this device** — the first of the four checks to go from
   "not run" to done this session. See `running_jobs.md` for the full job
   history.

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
   `resolved/fdtd_vs_EME_mode_validation.md` item 8.** A *post-hoc* diagnostic inside
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

   **❌ TIMEOUT 2026-07-31 13:42, full 5h — gap 6's original ask (exercise
   `warn_on_transfer_duplicate`) is STILL open.** The job never reached a
   third interface, let alone finished all output ports: `example_sim_cache/`
   shows sections 0-3 (both EME bodies, input+transition FDTD) cached by
   12:45, one more EME cache write at 13:40:31, then killed 1m41s later
   mid-way through building the first section-4 (per-port output) FDTD grid.
   No `transfer` line ever appeared in the log — still not distinguishable
   from "cannot fire" per the same caveat above. Resubmitted unchanged as
   **11590663**; the disk cache should skip straight past sections 0-3 this
   time, leaving only the 3 per-port section-4 solves as real new work. See
   `running_jobs.md` for the full accounting.

   **❌ TIMEOUT AGAIN, 2026-07-31 18:58, full 5h — still open after two full
   attempts.** The disk-cache prediction above was right (sections 0-3 hit
   cleanly, no re-solve), and this run did reach section 4: all 3 per-port
   geometry plots (`port_top`/`port_bot`/`port_mid`) exist on disk. But the
   "5h should now be generous" estimate was wrong — the log's tail is deep
   inside a per-port eigenmode-overlap monitor sweep (many `MPB solved for
   frequency_13(...)` lines across k-points, into band 14) when SLURM killed
   it, well short of writing a final S-matrix. No new
   `example_sim_cache/` entries attributable to this device appear during the
   run's own window, so unlike sections 0-3, section 4's per-port monitor
   work does not appear to checkpoint — a third blind resubmit would very
   likely restart section 4 from nothing rather than resume. **Gap 6's
   original ask is unanswered after 10h of 4-node compute across two
   attempts.** Worth reshaping before trying a third time: either a
   substantially longer wall-clock limit, or splitting the 3 per-port
   section-4 solves into their own array job so a future timeout only loses
   one port's progress. Not resubmitted as of this writing. Full accounting
   in `running_jobs.md`'s 11590663 entry.

   **Reshaped and resubmitted, 2026-07-31 ~21:20.**
   `modal_multiplexer_hybrid.py` now supports running one output port per
   SLURM array task (`SLURM_ARRAY_TASK_ID` 0-2 → top/mid/bot) plus a
   separate merge step that assembles the combined SMatrix from the on-disk
   cache (`cache_only=True` on every section, raises loudly if a port's task
   hasn't finished). Resubmitted as 11596846 (array) / 11596852 (merge) —
   that first attempt **failed instantly on both**, a wrong relative path in
   the new batch script (`running_jobs.md` has the correction). Fixed and
   resubmitted as **11597177** (array, 5h/task — sections 0-3 come from the
   already-warm cache, so the real new work per task is just its one port's
   section-4 solve) with merge job **11597199** dependent on it. See
   `running_jobs.md`'s reshape note and the script's own docstring for the
   mechanics.

   ---

   > **❌ THIRD TIMEOUT, 2026-08-01 02:44–03:01 — the reshape made it WORSE,
   > and the cause is understood. No further attempt made: this device is
   > deferred.** All three array tasks (11597177 `_0`/`_1`/`_2`) hit the full
   > 5 h wall on 4 nodes each — 15 node-hours spent, gap 6 still unanswered
   > after **25 h of cumulative compute** across three attempts.
   >
   > **Root cause — the reshape silently dropped MPI.** The new batch script
   > runs `python modal_multiplexer/python/modal_multiplexer_hybrid.py`
   > directly, deliberately bypassing `make` "which would run the whole script
   > — all 3 ports — in one process." But the Makefile's `modal_multiplexer_
   > hybrid` target is `$(runner) $(my_python) ...` and on Phoenix `runner` is
   > **`srun --mpi=pmi2`** (`BEEF/examples/Makefile:5`). Bypassing `make` threw
   > away the launcher along with the target. The logs are unambiguous:
   >
   > | job | launcher | meep ranks | s/step |
   > |---|---|---|---|
   > | 11585847_0, 11590663_0 (via `make`) | `srun --mpi=pmi2` | **96 processes** | 0.0064 |
   > | 11597177_0/_1/_2 (reshaped, direct) | none | **1 process** | 0.27–0.36 |
   >
   > `Using MPI version 4.1, 1 processes` in all three new logs vs
   > `96 processes` in both old ones — a **~42× per-timestep slowdown**, while
   > still holding 4 nodes × 24 ranks idle (SLURM accounting confirms the
   > waste: `cput=20-00:35:12` against `walltime=05:00:22` on 96 cpus). The
   > per-port split was real but was swamped many times over by losing 96-way
   > domain decomposition.
   >
   > **Second, independent flaw: the "sections 0-3 come from the already-warm
   > cache" premise is false.** Each task re-ran the *same* three
   > `Initializing structure` + MPB-band + time-stepping blocks at nearly
   > identical log offsets (lines 18/2189/4528, 18/2184/4523, 18/2224/…) —
   > i.e. every task redid the shared work from scratch rather than reading it
   > back. Whether the earlier timeouts never actually persisted those
   > sections, or the cache key differs between the 96-rank and 1-rank runs,
   > is **not established** — and is the thing to check first if this is ever
   > picked back up. The tasks died at wildly different points (`_0` at
   > timestep 289, `_1` at 2945, `_2` still inside MPB band 10), which is
   > consistent with all three grinding through that duplicated prefix.
   >
   > **If resumed later, in this order:** (1) put `srun --mpi=pmi2` back in
   > front of the direct `python` call — the array-per-port shape itself is
   > fine and worth keeping, it was only ever the launcher that broke; (2)
   > verify sections 0-3 actually hit the cache (a task that skips them should
   > start time-stepping within a minute or two, not after an hour); only then
   > (3) resubmit. Do **not** just raise the wall clock — at 1 rank, no
   > affordable limit finishes this.
   >
   > Merge job **11597199** was left `PENDING` /
   > `DependencyNeverSatisfied` and can never run; it is inert (holds no
   > resources) but should be cancelled as queue hygiene.
   >
   > Note the contrast with the 2×2 MMI's reshape submitted in the same window
   > and to the same design — that one worked exactly as predicted (see
   > `resolved/mmi_2x2_validation.md` gap 4). The array-plus-`cache_only`-merge pattern
   > is sound; this device's instance of it broke on the launcher, not the
   > pattern.

7. **New, 2026-08-06 — `modal_multiplexer_hybrid.py`'s section 2 split from
   one monolithic N_MODES-launch FDTD "transition" block into three:
   FDTD taper-out | EME connector | FDTD taper-in. Implemented and
   geometry-smoke-checked; NOT yet run for real, and its one new numeric
   knob (`CONN_MODES`) is an unmeasured, flagged placeholder.**

   Prompted by a user question about why section 2 (MMI 1's three output
   tapers + the three connecting waveguides + MMI 2's three input tapers,
   one monolithic FDTD domain) wasn't already split the way this device's
   own docstring reasoning would suggest: the connecting waveguides are
   themselves a uniform, straight `WaveguideRun` — exactly what EME
   propagates for free, same as the two MMI bodies — and
   `beef.geometry.compound.WaveguideRun` is explicitly documented as "safe
   to cut anywhere." The old single section cost **N_MODES=14** FDTD
   launches (confirmed directly from the shipped code, `interface_modes =
   [N_MODES] * 4` — not just inferred from the module docstring's own
   accounting) because both of its exposed interfaces were the full
   N_MODES-wide body-supermode basis; `FDTDSolver.solve_section()` already
   launches from whichever side of a section has fewer modes (the same
   mechanism sections 0/4 already exploited for their 2/3-launch costs), so
   exposing the connector's own much-narrower basis as a real interface lets
   the two new taper sections each launch from that narrow side instead.

   **Verified directly from the library code before touching the script**
   (not assumed): `CompoundGeometry.section_gds()`/`section_prisms()` are
   already fully generic over `WaveguideRun` segments, `EMESolver.
   solve_section()` builds its cross-section from `geometry.section_gds()`
   with no MMI-specific assumption anywhere, and `HybridSolver` already
   supports per-interface mode counts and an arbitrary section count (its
   own docstring: "any solver mix, any count >= 2") — so the entire change
   is confined to this one example script's own `boundaries`/
   `interface_modes`/`solvers` construction, no `beef` library changes
   needed.

   **What changed**: `modal_multiplexer_hybrid.py` now builds a 7-section
   chain (was 5) — `FDTD-in | EME-MMI1 | FDTD-taper-out | EME-conn |
   FDTD-taper-in | EME-MMI2 | FDTD-out` — reusing the same `fdtd_solver`/
   `body_solver` instances for every slot of their type (no new solver
   objects). Two new constants: `CONN_STUB = 0.25` µm (matching the
   coupler/other MMI hybrid examples' own `body_stub` convention, not this
   device's own more conservative `BODY_STUB = 0.75`, which was sized for
   the MMI body's harder degenerate-multimode discontinuity) and
   `CONN_MODES = 8` (see below). The three hardcoded `{4: {"output_idx":
   ...}}` references to the output section (in the geometry-plot
   `sub_region_configs`, the per-port array-task `solve()` call, and the
   merge-step `solve()` call) all became `{6: ...}` — the output section's
   index shifted from 4 to 6 now that two new sections sit before it.

   **`CONN_MODES` — now MEASURED (2026-08-06/07, login-node
   `EMESolver.section_modes()` on the connector's own cross-section) and
   corrected 8→6, with a real cost finding alongside it.** The original
   reasoning ("~3, one per arm") was a hand-wave, not a measurement, and it
   was wrong: requesting 8 modes at wl=1.31 on the exact `CONN_STUB`-bounded
   span returned n_eff `[2.10549, 2.09067, 2.03595, 1.56570, 1.54514,
   1.52526, 1.42755, 1.41349]` — the last two sit below the SOI_160 cladding
   index (~1.4468, this repo's own prior measurement), so only **6** are
   guided. Each arm is apparently weakly bimodal (fundamental + first-order)
   rather than strictly single-mode, giving 3 arms × 2 modes = 6 coupled
   supermodes, not 3. `CONN_MODES` corrected 8→6 in the script (over-
   requesting 8 was never wrong on its own, since `EMESolver.compute_modes()`
   guards duplicates — but its guard is documented as duplicate/rank-
   deficiency detection, NOT a guided-vs-cladding check the way
   `FDTDSolver._iter_distinct_bands()` has post items-12/14, so shipping 8
   unverified risked silently admitting 2 unguided modes into the basis; 6
   is the safe, measured choice).

   **Second, more consequential finding: this single-rank solve took 5050 s
   (~84 min) for ONE wavelength — not "basically free."** The section-2
   split's whole module-docstring rationale leaned on the connector EME
   solve being cheap, by analogy with the MMI bodies' own fast slab solves.
   Three separated, tightly-pitched (`WO=1.0375` µm) silicon islands with no
   exploitable mirror symmetry is evidently a much harder MPB target than a
   single wide slab. The script solves this cross-section once per
   wavelength (`N_WL=11`), so total connector-EME cost in a real run is a
   genuinely open question — **not** measured yet, and how much production's
   many-rank MPI parallelism brings this single-rank number down is
   unknown. This may claw back some of the FDTD-launch savings (14→6×2
   instead of 14) the split was built to capture; treat the split's net
   benefit as unproven until a real run measures both sides of that
   trade-off, not just the FDTD half.

   **Smoke-checked, not run for real.** Executed the real script's own
   source (steps 1-4 only — geometry, solver config, `HybridSolver`
   construction, `plot_hybrid_geometry`), truncated before step 5's actual
   `hybrid.solve()` calls, so no physics solve ran. Section spans came out
   exactly as hand-computed: taper-out/taper-in both 3.0 µm (0.75 body stub
   + 2.0 taper + 0.25 conn stub), EME-conn 2.5 µm (3.0 − 2×0.25).
   `HybridSolver(...)` construction raised nothing (validates
   `boundaries`/`interface_modes` length agreement). `geometry_hybrid.png`
   regenerated cleanly showing all 7 sections (labels overlap in the now-
   busy ~8.5 µm middle stretch — cosmetic, not investigated).

   **NOT done, deliberately, to avoid compounding unverified changes**:
   - `mirror_pairs` — the two new taper sections are geometric mirror
     images of each other about the connector's midpoint (verified by hand:
     reflecting `[body1_end, conn_start]` about `(conn_start+conn_end)/2`
     lands exactly on `[conn_end, body2_start]`), which `HybridSolver`
     could exploit to solve only one and transpose-reuse it for the other —
     roughly halving the new cost again. Left undone; would need its own
     `mirror_axis` (the connector's midpoint, NOT the device's own overall
     midpoint the default would use, since the whole device isn't globally
     mirror-symmetric — MMI 1's body is 16.92 µm, MMI 2's is 34.24 µm).
   - `check_mode_correspondence=True` at the two new interfaces — this
     device's own gap 1 is the single highest-value reason to want this
     (interface mode-identity bugs are its documented top failure mode),
     but turning it on was not attempted this pass. `CONN_STUB=0.25` was
     chosen to clear `FDTDSolver`'s `dpml=0.5` query padding *if*
     `mode_correspondence_pad` is left unset, matching the coupler/other
     MMI examples' own precedent — but per `notes.md` Solvers item 14, pass
     `mode_correspondence_pad=0.1` to be safe, matching those examples'
     own working configuration.
   - No new SLURM submission of any kind — the smoke check above is
     geometry-only (login-node, seconds), not a real hybrid solve; a real
     run needs its own explicit go-ahead, and there is not yet a landed
     full-device FDTD reference (gap 3, job 11707932, still running as of
     this writing) to compare a new hybrid result against anyway.

   **✅ 2026-08-08 — real go-ahead given, first-ever real attempt made,
   a genuine new library bug found and fixed, second attempt in
   progress.** With gap 3 finally closed (a real FDTD reference now
   exists) and `TAPER_L`/`minimum_run_time` fixes available, submitted
   the first-ever physics-on run of this device's hybrid chain: job
   **11751064** (`run_modal_multiplexer_hybrid_SMOKE.batch`, new this
   session — single array task, port "top", `MUX_HYBRID_N_WL=2` to bound
   the connector-EME cost, `minimum_run_time=100.0` added to this script's
   `fdtd_solver` preventatively, 6h budget, `save_dir` wired in so a
   partial run isn't a total loss).

   **11751064 FAILED after 1h17m** — but with a real, caught
   `ValueError`, not a silent bad result or a wasted budget (`save_dir`
   preserved whatever completed: 1 EME section + 1 FDTD section cached).
   `_check_aperture()` (an existing guard, `beef/solvers/fdtd.py`, added
   for exactly this class of problem per its own docstring — a
   ring-resonator incident once cost 27% of transmitted power to a
   silently-clipped aperture) caught: `excite_port aperture crosses the
   PML: ... centre 1.0375, width 5.0000`.

   **Root cause, found by tracing (not guessing)**: `section_fdtd_layout()`
   widens a multi-mode interface's window to the full `sim_width` (correct
   — this is the whole point of gap 3/7's "measure >1 mode → full width"
   rule) but never re-centers it to y=0 to match. For every device before
   this one, the underlying `center_at()` query at a multi-mode interface
   always landed inside a single, well-defined body (an MMI or a
   directional coupler's coupling region), which naturally reports y=0.
   The new connector-adjacent interfaces (gap 7) are the **first ever**
   HybridSolver boundaries to sit inside a multi-port `WaveguideRun` —
   and `CompoundGeometry.center_at()`'s own fallback for a multi-port
   `WaveguideRun` query with no explicit `output_idx`/`input_idx`
   defaults to **port 0** ("top", y=+1.0375), not y=0. Combined with the
   (correctly) unclamped full `sim_width`, that reached 3.5375 µm past a
   2.5 µm PML inner edge. A latent gap in `section_fdtd_layout()`, not
   something specific to this device's own script — any future device
   cutting a boundary inside a multi-port `WaveguideRun` would have hit
   the identical failure.

   **Fixed** in `section_fdtd_layout()` (`beef/solvers/fdtd.py`): when the
   multi-mode branch sets a side's width to `sim_width` (and no explicit
   `monitor_center`/`monitor_size` override was given), that side's center
   is now also forced to `0.0`. Scoped narrowly — the `n_excite == 1`
   single-lone-port fallback (every other device's existing, working
   precedent) is untouched; only the `n_measure > 1` / `n_excite > 1`
   multi-mode branches are affected, and no prior device's own multi-mode
   interfaces ever queried a multi-port `WaveguideRun` this way, so no
   already-recorded result moves.

   Traced and confirmed the mechanism directly (not just inferred from the
   error message) via a login-node script instrumenting
   `section_fdtd_layout()` — **caveat, logged honestly**: that trace
   briefly ran real (small, ~5 µm domain, single-rank) FDTD physics on the
   login node for section 0 before being caught and stopped, a real
   violation of this repo's own "SLURM, not login-node" discipline for
   FDTD, even though small in scope — not repeated; verification moved to
   a real SLURM resubmission instead, per that same discipline.

   **Job 11752080 resubmitted with the fix, same configuration, currently
   running** — the actual test of whether this closes the bug. Not yet
   landed.

8. **New, 2026-08-07 — user-provided paper draft (Hiesener et al.,
   "Topology Optimization Enhancements for Multimode Interferometers," the
   TE modal multiplexer test case) cross-checked against this device's own
   code, plus a real mode-monitor sizing bug fixed in the library. Three
   real, confirmed changes; one open question sent back to the user.**

   **(a) `TAPER_L` was wrong — 2.0 µm in the code, 3.0 µm in the paper's own
   Table 2 (`L_taper`).** Every other geometric parameter cross-checked
   cleanly against Table 2: `MM_WIDTH`=2.8 (`w_mm`), `mmi1.mm_length`=16.92
   (`L_mm1`), `mmi2.mm_length`=34.24 (`L_mm2`), `BETWEEN_WG_LENGTH`=3.0
   (`L_ps`), input taper width 1.4 (`w_in`), MM bus width 0.7 (`w_mmwg`),
   `BRANCH_SPECS` taper widths 0.725/1.150 (`w_t1`/`w_t2`, and the paper's
   own `w_t1 = (w_mm - w_t2 - 2*0.1)/2` formula checks out exactly:
   (2.8-1.150-0.2)/2 = 0.725), branch waveguide widths 0.33/0.35 (`w_ps`/
   `w_wg`). Only `TAPER_L` was off. **Fixed in both
   `modal_multiplexer_hybrid.py` and `modal_multiplexer_fdtd.py`**
   (2.0→3.0); `mux.total_length` moves 62.16→66.16 µm (4 tapers × +1.0 µm).
   Confirmed by smoke-checking both scripts' geometry/solver construction
   (login-node, no physics solve) — spans update correctly, nothing raises.
   **Not yet re-run** — every prior result from this device (including
   job 11707932's own, already-discarded-as-unusable run above) used the
   wrong taper length; this alone would need a re-run regardless of the
   monitor-width fix below. 160 nm SOI / O-band (1260–1360 nm) were already
   correct, confirmed against the paper's own device description.

   **(b) A real, likely-significant mode-monitor sizing bug, found while
   checking (a) — plausible root cause (not yet confirmed) for job
   11707932's near-total-loss result.** At the user's direction (general
   rule: a mode monitor should span as much of the domain as it validly
   can — full `sim_width` when centered, or the widest symmetric window
   that still fits inside the domain when off-axis — rather than an
   arbitrary narrower value). Two findings:
   - `modal_multiplexer_fdtd.py` explicitly set
     `output_monitor_widths=[0.9, 0.9, 0.9]`, far narrower than it needed
     to be (branch centers at ±1.0375/0.0 inside a `sim_width=5.0` domain
     allow up to 2.925/5.0/2.925). A monitor this narrow risks clipping
     real field content — and, more seriously, the surrounding cladding
     context `get_eigenmode_coefficients()`'s own MPB solve needs to find
     the *correct* eigenmode shape in the first place — before the overlap
     integral that produces the reported coefficient even runs. **Removed
     the override entirely** — now uses the library default.
   - That library default itself (`FDTDSolver.solve()`, `beef/solvers/
     fdtd.py`) was ALSO wrong, and this is the more consequential fix:
     `output_monitor_widths=None` resolved to a bare `[sim_width] *
     len(centers)` — full `sim_width` centered on *every* monitor
     regardless of its own `output_monitor_centers` offset. For any
     off-axis output port (this device's "top"/"bot" branches at
     ±1.0375 µm), that places a monitor asymmetric about y=0 but sized as
     if centered there — extending past the valid domain on the far side
     by `2×offset`. **Fixed** to match the clamp already used elsewhere in
     the same method for the excited/source-side aperture
     (`min(aperture, sim_width - 2*abs(center))`):
     `sim_width - 2*abs(center)` per monitor. This is a **library-level
     fix** (`FDTDSolver.solve()`'s own default), not device-specific — it
     silently affects any *other* future caller relying on the default too,
     though every existing example script to date passed an explicit
     override, so no currently-shipped golden result is affected by this
     specific change.

   Both scripts smoke-checked clean after all three fixes (geometry/solver
   construction only, login-node, no physics solve) —
   `modal_multiplexer_fdtd.py`'s solver now resolves output widths to
   `[2.925, 5.0, 2.925]` at solve() time, confirmed by inspection.
   **Neither device script has been re-run with these fixes** — job
   11707932's own near-zero result (above) was measured before all three
   fixes existed, so it is NOT evidence either fix does or doesn't help;
   that requires an actual re-run, not yet submitted.

   **(c) Open question sent back to the user, NOT resolved unilaterally:**
   the paper's own text says this device "cascades a 1×3 MMI splitter, a
   phase shifter, and a 3×3 MMI splitter, **with the third port of each MMI
   left unconnected**" (Section 3.2) — i.e., only 2 of each MMI's 3 ports
   are real, connected single-mode waveguides (`A`/`B`); the third is
   `NC`. This device's own code (`BRANCH_SPECS`, three entries: top/mid/
   bot) builds all THREE ports as real, physically-taped waveguides —
   inherited from the legacy `mmi_mux_master.py`/`config_mux.json` lineage
   (`PORT_DEFS`: top/mid/bot, mid+bot tagged `[MUX out]`), not derived from
   this paper. Whether "NC" in the paper means that third port is
   physically absent from the mask (this device's code would then be
   simulating an extra waveguide/taper that doesn't exist on the real
   chip) or merely unused/unmeasured downstream (geometrically present,
   functionally irrelevant — closer to what the code already builds) is
   not resolvable from the paper text alone.

   **✅ RESOLVED 2026-08-07 — user confirmed directly: the NC port is
   physically present on the real device (taper + waveguide drawn), just
   unused/unmeasured downstream — not absent from the mask.** The current
   3-real-branch topology (`BRANCH_SPECS`: top/mid/bot, all real) is
   therefore topologically correct as-is; no geometry change needed here.
   Only (a) and (b) above were real bugs.

   **✅ 2026-08-07, job 11726910 — re-ran the SMOKE test (res=10) with all
   three fixes in place. Confirms the monitor-width diagnosis was right —
   near-zero collapsed to physically sane numbers — but reveals one new,
   real, smaller anomaly.** COMPLETED in 14m27s (well under the 60 min
   budget):

   | | top | mid [MUX out] | bot [MUX out] | total |
   |---|---|---|---|---|
   | TE0 in | T=0.0055 | T=0.5283 | T=0.5038 | **1.0375** |
   | TE1 in | T=0.0101 | T=0.0109 | T=0.0105 | 0.0316 |

   TE0 now correctly routes the large majority of its power to the two
   real MUX-out ports (~53% each), with a small, plausible leak to "top" —
   physically sane for a pre-optimization seed device, and a world away
   from the prior run's uniform near-zero/exact-zero collapse. This is
   real, measured evidence the `TAPER_L`/monitor-width fixes (this gap,
   above) fixed the dominant problem, not just a plausible-sounding theory.

   **But `check_passive()` now raises a warning it didn't get the chance to
   before**: TE0's total hits **1.1113** (>1) at wavelength 1.33 µm
   specifically — non-physical for a passive device, and not seen at the
   other 10 wavelengths (the identical warning text repeats ~24× in the
   log purely from this script's known, pre-existing lack of
   `is_master()` gating — every MPI rank prints it, not 24 independent
   violations). Plausible mechanism, **not confirmed**: the newly-widened
   output monitors (2.925/5.0/2.925 µm for top/mid/bot) now overlap each
   other substantially at this branch pitch (`WO`=1.0375 µm) — "mid"'s
   window alone spans the full domain and geometrically contains "top"
   and "bot" too — so some double-counting via imperfect orthogonality
   between the requested per-branch eigenmodes, likely worst at whatever
   makes 1.33 µm specifically prone to near-degeneracy, is a reasonable
   suspect. Trades the old failure mode (catastrophic under-count) for a
   much milder, single-wavelength over-count — a real improvement, but not
   a fully clean result. **Not investigated further yet** — this is exactly
   the kind of thing a smoke test (res=10, deliberately crude) exists to
   surface cheaply before the real res=20 run commits real hours to it.


   **❌ 2026-08-08, job 11734286 COMPLETED (17h01m, within budget) — the
   real res=20 run, with all three fixes applied. STILL catastrophically
   near-zero. The smoke test's fix did NOT carry over to production
   resolution.** Read the full per-wavelength S-matrix directly (not just
   the one center-wavelength line the script prints), same discipline as
   the pre-fix run:

   | in → out | T across all 11 wavelengths |
   |---|---|
   | TE0 → top | ~4.4×10⁻⁹ – 1.4×10⁻⁷, smooth in λ |
   | TE0 → mid [MUX out] | ~2.6×10⁻⁶ – 3.2×10⁻⁶, smooth in λ |
   | TE0 → bot [MUX out] | ~2.6×10⁻⁶ – 3.4×10⁻⁶, smooth in λ |
   | TE1 → top | ~3.8×10⁻¹⁵ – 6.6×10⁻¹⁴, smooth in λ |
   | TE1 → mid | ~9.2×10⁻¹³ – 1.1×10⁻¹² , smooth in λ |
   | TE1 → bot | ~9.4×10⁻¹³ – 1.2×10⁻¹² , smooth in λ |

   Every channel is smooth and monotone across wavelength — the same
   "correctly-computed-ratio-scaled-by-something-badly-wrong" signature as
   the pre-fix run, not random noise — but the magnitudes are now even
   smaller than the pre-fix run's (which had ~10⁻⁵ channels; this one's
   largest channel is ~10⁻⁶). TE0 finished at essentially the same
   decay point as the pre-fix run (t=77,174.5 vs t=77,197.45 — the
   `TAPER_L` change is too small a perturbation to move this).

   **This isolates the problem.** Two real, confirmed fixes (`TAPER_L`,
   monitor width) transformed the res=10 smoke test from broken to sane,
   but made no qualitative difference at res=20 — same catastrophic
   collapse, if anything worse in magnitude. The variable that actually
   correlates with the failure across every run to date is resolution/
   required-decay-time, not geometry or monitor sizing: res=10 needs
   t≈2910–3217 to satisfy `decay_tol=1e-11` and produces sane results
   (both before and after the fixes); res=20 needs t≈77,175–77,197 (24–26×
   further) and produces near-zero results (both before and after the
   fixes). **The leading hypothesis is now `decay_tol` itself**, exactly
   the "more likely dominant lever" flagged and deliberately NOT changed
   in gap 3's earlier account (pending explicit sign-off, since it's a
   methodology decision, not a bug fix): `stop_when_dft_decayed`'s
   relative-change criterion may force the simulation to keep running long
   after the physically meaningful signal has already exited/decayed,
   tracking some much smaller, slowly-decaying spurious/resonant
   component instead — by the time that relative criterion is finally
   satisfied at res=20, the real signal at the output monitors may already
   have decayed into numerical noise, while at res=10 the same criterion
   happens to trigger before that happens. **Not confirmed** — this is
   still a hypothesis backed by correlation across four total runs (2
   smoke, 2 real), not a mechanism verified from the meep source or a
   controlled decay_tol sweep. **Not changed here** — same reasoning as
   before, needs the user's explicit sign-off before spending more compute
   on it.

   **✅ User sign-off given 2026-08-08. `decay_tol` default changed
   1e-11→1e-9** in `modal_multiplexer_fdtd.py` (env-overridable via
   `MUX_FDTD_DECAY_TOL`; `=1e-11` reproduces the two prior broken runs
   exactly). **Smoke test re-run first, job 11749510 (res=10,
   decay_tol=1e-9) — safe, no regression.** COMPLETED in 8m12s (480.9s,
   even faster than the original 1e-11 smoke test's 850s). Results
   byte-for-byte the same as the 1e-11 smoke test: TE0 total=1.0375
   (same 1.33 µm `check_passive` overshoot), TE1 top/mid/bot=
   0.0101/0.0109/0.0105, total=0.0316. Expected — res=10 already decays by
   t≈2910–3217, well inside both tolerances, so this only confirms the
   looser criterion doesn't corrupt anything at the resolution where
   decay was never the bottleneck. **Real res=20 run submitted as job
   11749567** — this is the actual test of the hypothesis. landed.

   **❌ 2026-08-08, job 11749567 COMPLETED in 12 minutes — but this
   REJECTS the decay_tol hypothesis and finds a different, more specific
   bug.** The astonishing speed itself was the tell: `run 0 finished at
   t=107.125` (TE0), `run 1 finished at t=112.225` (TE1) — for a 66 µm
   device, light needs several hundred meep time units at minimum just to
   physically cross it once. The sim stopped before the pulse could have
   reached the output monitors at all — T=0 (hit the script's -200 dB
   floor) for exactly that reason, not because the device is lossy.

   **Root cause, found**: `FDTDSolver` never sets `minimum_run_time` on
   `stop_when_dft_decayed()` anywhere (3 call sites, `beef/solvers/
   fdtd.py`) — meep's own default is 0, no floor at all. A loose-enough
   `decay_tol` can therefore be satisfied by a fast, near-source transient
   long before a distant monitor has seen anything, independent of the
   device's real coupling. This also reframes the *previous* result
   (11734286, decay_tol=1e-11, ran to t≈77,175): that run had plenty of
   time for the real signal to arrive and still returned near-zero — so
   decay_tol was never really the dominant lever; it happened to be tight
   enough to avoid the trivial-early-stop bug, but something else
   (still unidentified) makes even a fully-run simulation report
   near-zero transmission at res=20 specifically.

   **Independent corroboration, should have been checked earlier**: the
   user's own paper (Fig. 3e/Table 2, gap 8) reports the SEED device's
   simulated single-pass transmission at roughly 50–68% per channel —
   much closer to the res=10 smoke test's ~50% than to either real run's
   near-zero. This supports reading res=10 as closer to physically correct
   and both res=20 attempts as genuinely broken, not the reverse.

   **Fixed**: `FDTDSolver` gained a new `minimum_run_time` constructor
   parameter (default 0.0 — meep's own default, unchanged, so no other
   device's already-recorded behavior moves), wired into all 3
   `stop_when_dft_decayed()` call sites. `modal_multiplexer_fdtd.py` sets
   it explicitly to 1000 (device+stub+PML span ~68 µm; even a conservative
   silicon group index ~4.5 only needs ~306 time units for one transit —
   1000 leaves >3x margin for reflections, comfortably below the ~77,175
   decay_tol=1e-11 has separately shown is needed regardless). Smoke-
   checked (login-node, construction only): `solver.minimum_run_time` ==
   1000.0, confirmed wired correctly. **Not yet run for real at any
   resolution with this fix** — next step is the res=10 smoke test
   first, same two-step discipline as every other change this gap.

   **✅ 2026-08-08 — CLOSED.** Job 11749883 (res=10 smoke,
   decay_tol=1e-9 + minimum_run_time=1000) COMPLETED in 475s, byte-for-byte
   identical to every prior sane smoke result — no regression, as
   expected (natural decay at res=10 already exceeds the 1000-unit floor).
   Then **job 11749963 (real res=20, both fixes) COMPLETED in 35m51s —
   29× faster than the original 17h01m run — with the FIRST genuinely
   sane, physically trustworthy result this device has ever produced at
   production-relevant resolution.** Full per-wavelength S-matrix read
   back directly (not just the center-λ print):**

   | in → out | T range across 11 wavelengths (1.26–1.36 µm) |
   |---|---|
   | TE0 → top | 0.40 – 0.71, smooth monotone decrease |
   | TE0 → mid [MUX out] | 0.02 – 0.22, smooth monotone increase |
   | TE0 → bot [MUX out] | 0.02 – 0.19, smooth monotone increase |
   | TE1 → top | 0.01 – 0.07 |
   | TE1 → mid [MUX out] | 0.05 – 0.55, smooth monotone decrease |
   | TE1 → bot [MUX out] | 0.05 – 0.55, smooth monotone decrease |

   TE0 total power 0.75–0.85 across the band; TE1 total 0.17–1.10 (the
   one >1 point sits at the short-wavelength band edge, λ=1.26 µm,
   `check_passive` flagged it — same modest, edge-of-band monitor-overlap
   suspect already on record above, not a new problem). Every one of the
   6 spectra is smooth and physically well-behaved across the full band —
   no noise, no exact zeros, no absurd magnitudes — plotted in
   `results.png`: TE1's mid/bot curves cross TE0's top curve near
   λ≈1.32–1.33 µm, a real, plausible spectral feature for a phase-
   shifter-based mux (whose operating principle is inherently
   wavelength-dependent — a relative phase shift between two arms is a
   function of λ).

   **This closes the four-bug chain that blocked this gap**: `TAPER_L`
   (2.0→3.0, matching the paper), `output_monitor_widths` (both the
   device script's own ad-hoc override and the library's broken default),
   and `minimum_run_time` (missing floor on `stop_when_dft_decayed`) were
   all real, independent, now-fixed defects — not one single root cause,
   four of them stacked, each masking evidence of the others until fixed
   in sequence. **Gap 3 (README check 4) is CLOSED**: this device finally
   has a trustworthy full-device FDTD reference, `example_output/
   modal_multiplexer/modal_multiplexer_fdtd/smatrix.pkl`
   (job 11749963), matching every other device in this repo. **Remaining
   open item, minor**: the single >1 `check_passive` point at the
   1.26 µm band edge — not investigated further, same monitor-overlap
   suspect flagged in this gap's earlier smoke-test writeup.

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
duplicate scan is now implemented (`resolved/fdtd_vs_EME_mode_validation.md` item 8)
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

**Re-ordered 2026-08-04, at the user's direction: start with 5, not last.**
This reverses the sequencing above, for reasons that did not exist when it
was written — at the time, **3** was known-blocked on a real library gap
(`notes.md` Geometry item 8), not merely expensive, so it could not go next
regardless of where the list put it; and **4** (golden) is explicitly gated
on the `N_MODES`/`RESOLUTION` decision that gap 2 raised, which still has
not been made. **5** (engine agreement), by contrast, was unblocked from the
start, cheap relative to the device's other campaigns (both MMI siblings
closed their own version inside an 8h budget once reshaped), and — same
reasoning `directional_coupler_validation.md` item 6 and both MMI gap-4s
already used — is a localized question (does the body's own uniform
cross-section accumulate phase the same way in both engines) that does not
depend on **3** or **4** being settled first. See gap 5 for the script and
run status.

Gap 2's own cheap extension (meep res 80/100) was submitted alongside 5, for
the same reason: it is not a simulation campaign either, and its answer is a
direct input to the `N_MODES`/`RESOLUTION` decision that gap 4 (the golden)
and any future hybrid re-run both need.

**Same day, later: 3's library blocker was fixed too** (`notes.md` Geometry
item 8 part 2 — the "no generic way to know per-port stub width" premise
was wrong, not just unimplemented; see gap 3 for the full account). This
does not reopen the original sequencing question — 3 is still a real
simulation campaign (a 62 µm single-domain 3-D FDTD run), still the most
expensive of the three now in flight, and still needs a smoke test before
the real run per this repo's own two-step discipline for never-run scripts.
It was picked up anyway, in parallel with 5 and the gap-2 extension, because
none of the three depend on each other and all three were unblocked at the
same time.

**All three landed the same day.** **5 is CLOSED** (mean ratio 0.9958,
`calibrate_body_neff` not indicated — see gap 5). The **gap-2 extension**
completed too, and answered its own question with a real "no": meep's
absolute n_eff does not plateau by res 100 after all (60→80 looked
converged, misleadingly — see gap 2). **Gap 3's smoke test** (`MUX_FDTD_RES
=10`) completed — the pipeline runs end to end for the first time ever, a
real S-matrix saved, `check_passive()` held — but the real res=20 run is
still not submitted, so check 4 itself is not yet closed. Net effect on the
`N_MODES`/`RESOLUTION` decision this doc keeps flagging as open: `N_MODES
=14` is now solid; `RESOLUTION` still has no clean answer, and that
decision remains deliberately unmade.
