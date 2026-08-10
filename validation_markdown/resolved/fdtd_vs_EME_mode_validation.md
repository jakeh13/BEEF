# FDTD ↔ EME mode validation — mechanism and to-do

Opened 2026-07-30. Cross-device, like `material_model_validation.md`: this is
about the *mechanism* by which BEEF establishes that an FDTD interface basis
and an EME body basis name the same modes, not about any one device's status
(those live in the per-device files).

> **Note 2026-08-01 — the material-dispersion default flip does NOT invalidate
> this document's numbers, and that is measured.** `dispersive` defaulted to
> True on 2026-08-01 (`material_model_validation.md`), changing the
> cache fingerprint and staling every cached entry and golden. The reflex is to
> mark every recorded number here frozen-era; that reflex is **wrong here**.
> Job **11608940** re-ran the 1×2's interface-mode audit under dispersion and
> returned **byte-identical** overlaps and n_eff. Reason: everything in this
> file is a *per-wavelength* mode solve. `section_mode_neff()` builds its
> medium at a scalar wavelength deliberately, and meow's `Environment(wl=…)`
> evaluates exact Sellmeier either way — a dispersive medium only differs from
> a frozen one *across a band*. Only genuinely broadband quantities (per-device
> gap 4 and gap 5) move. Re-audits are cheap and worth running per device
> (`AUDIT_ARGS` on `run_audit_interface_modes.batch` exists for exactly that),
> but expect them to confirm, not overturn. **Caveat:** do not launch the
> all-device default — the modal multiplexer requests 16 modes where ~14 are
> reachable at res 40 and will now RAISE on the guided-spectrum stop (item 14).
> Name the device.

## The requirement

BEEF must work for **any device geometry that supports a mode**, at **any
resolution pair**. That rules out most of what is in place today, so it is
worth being precise about what the requirement actually implies.

`HybridSolver.solve()` computes `M_K @ … @ M_1`, contracting over the mode
index shared at each interior interface. Three *separate* properties must hold
for that to be correct, and they are usually conflated:

- **(a) Identity / ordering** — index *i* of section *k*'s output basis is the
  same physical mode as index *i* of section *k+1*'s input basis.
- **(b) Completeness** — the truncated `interface_modes[k]`-mode basis captures
  the field crossing that plane; power in neglected modes is negligible.
- **(c) Distinctness / rank** — the *N* requested modes really are *N*
  independent modes.

Today only (a) is checked, and only partially. (c) is guarded on the
audit-only entry points as of 2026-07-30 (`modal_multiplexer_validation.md`
gap 6). **(b) is not checked anywhere at all** — and for a novel geometry it
cannot be assumed, because nobody knows in advance how many modes that
geometry needs.

## Why the existing mechanisms are not enough

| mechanism | where | addresses | why it fails the "any geometry" bar |
|---|---|---|---|
| n_eff nearest-neighbour | `tests/test_mode_identity.py::assert_modes_correspond()` | (a) | 1-D scalar heuristic: picks the wrong pairing outright once modes crowd in n_eff (three modal-multiplexer bands matched the same meow band). Requires **matched resolution** — and production deliberately uses mismatched (coupler: FDTD 40, EME 80). |
| transverse-E field overlap + Hungarian | `beef/mode_identity.py` | (a) | Right primitive, but only its **diagonal** is tested. Cannot see a mode present in one engine and missing from the other, so it is blind to (b). |
| runtime check | `HybridSolver(check_mode_correspondence=True)` | (a) | Same diagonal-only limitation, and **currently unverified** — see Open questions. |
| duplicate-band guard | `FDTDSolver._iter_distinct_bands()` | (c) | Sound, but only on `section_mode_neff()`/`section_mode_fields()`, never on `solve_section()`, the basis production actually uses. |

Two hard constraints any cross-engine field comparison must respect, both
established empirically and both easy to get wrong (each cost a debugging
cycle this session — see `notes.md` Solvers #14):

1. **No cross-engine E×H overlap is available.** The textbook mode overlap
   `∫(Eₐ × H_b*)·x̂ dA` mixes engine A's E with engine B's H, and meep and meow
   report E:H ~375× apart. Per-engine power normalization does **not** rescue
   it: with `E_A = αE`, `H_A = βH`, `E_B = γE`, `H_B = δH`, the normalized cross
   term comes to `√(αδ/βγ)`, which is 1 only if `αδ = βγ`. The symmetrized form
   fails the same way. So a **field-type-pure** overlap (E against E, or H
   against H) is not a convenience — it is the only cross-engine form available
   without an explicit E:H calibration per engine pair.
2. **The comparison mesh must resolve the physics.** An EME side *solves* on
   whatever mesh it is handed, so an under-resolved comparison grid makes one
   side's modes physically wrong rather than merely coarsely sampled.

## Recommended mechanism — layered, not a single test

No single test covers (a)+(b)+(c) at acceptable cost. The recommendation is
four layers, cheapest first, each catching what the one above cannot.

### Layer 1 (primary, per interface, cheap) — treat η as a *matrix*, not a diagonal

Keep the normalized transverse-field overlap matrix `η` and the Hungarian
assignment, but test the whole matrix. If both bases span the same subspace and
correspond, **η is close to a permutation matrix**. That single statement
covers far more than the current diagonal check, and costs nothing extra — it
is the same matrix already being computed:

- **Diagonal** `η[i, corr[i]]` — the matched pair's similarity. *(today's check)*
- **Margin** best-minus-second-best per row — flags genuine degeneracy rather
  than silently resolving it. *(already computed, currently unused)*
- **Row completeness** `Σⱼ |η[i,j]|²` — this is the (b) test, and it is free.
  For an orthonormal basis `{b_j}` and unit `a_i`, that sum is 1 **iff** `a_i`
  lies in `span{b_j}`. A row summing to noticeably less than 1 means engine A
  has a mode that engine B's retained basis cannot represent — i.e. truncation
  or basis mismatch, exactly what a novel geometry will hit and what a
  diagonal-only check cannot see.
- **Column completeness** — the same in reverse, catching a mode engine B
  found that engine A missed.
- **Off-diagonal mass** — two rows peaking on one column means degeneracy or a
  duplicate that survived layer (c).

Preconditions, all established: guided-mode filtering (radiation modes have no
engine-independent shape); an adequately resolved shared mesh; a field-type-pure
overlap. Computing the **H-only overlap independently** and requiring it to
agree with the E-only one is a cheap, genuinely independent corroboration —
currently the H components are fetched and discarded.

### Layer 2 (necessary companion) — `interface_modes` convergence

Sweep `interface_modes[k]` and require the final S-matrix to plateau. This is
the only way to establish (b) for a geometry whose modal content is not known
in advance, and it is device-agnostic. Relatively cheap on the EME side (modes
are reused); the FDTD side costs one run per mode. Non-negotiable for the
"any geometry" claim.

### Layer 3 (authoritative, per device) — cut-position invariance

Move the FDTD/EME boundary along a uniform stretch and require the final
S-matrix not to change. This is the strongest available test and the only one
that needs **no** mode comparison, no shared grid, and no cross-engine field
calibration — so it is immune to every constraint above and applies to any
geometry at any resolution pair.

It is also already validated as *the* meaningful property, in this repo, on
real data: the coupler's stub-sweep drift collapsed from +0.1085 to −0.0007
µm/µm (167×, below run scatter) once the body was β-calibrated, and
`validation_markdown/README.md` puts it plainly — the answer ceasing to depend
on where the cut is made "is the property that makes a hybrid decomposition
meaningful at all." A broken correspondence *cannot* be cut-invariant, so this
subsumes layers 1–2 as a pass/fail criterion. Its weakness is diagnostic
resolution, not validity: it tells you something is wrong, not which mode.
Hence layer 1 stays, as the cheap localizer.

### Layer 4 (standing control) — same-engine self-comparison

Always run meep-vs-meep and meow-vs-meow at the same interface as controls.
Both must give identity with η ≈ 1; anything else is a **checker** bug, not a
device finding. This session produced two spurious "correspondence failures"
that a control would have caught immediately, and the ad-hoc attempt to run
one after the fact timed out. This belongs in the tooling, not in a scratch
script.

### What this means for the current implementation

`HybridSolver(check_mode_correspondence=True)` is layer 1, restricted to the
diagonal. The upgrade path is additive: the matrix statistics are computable
from the `η` it already builds.

**2026-07-31 — layers 2 and 3 are now both reusable, device-agnostic
functions**, `beef.interface_modes_convergence.interface_modes_convergence()`
and `beef.cut_position_invariance.cut_position_invariance()`. Both take
already-solved `SMatrix` sequences and leave the device-specific sweep
construction to the caller; both are unexercised against a real device sweep
as of this writing (item 6's tool is brand new, item 7's is still awaiting
its own migration — see each item).

## To-do

Ordered. Item 0 was discovered by item 2's own machinery on its first
full-scale run and now leads, because it invalidates the EME side of every
comparison here.

- [x] **0. NEW 2026-07-30 — `EMESolver`/meow returns duplicate eigenmodes too.
      The `FDTDSolver` bug, on the other engine, unguarded.** **VERIFIED
      2026-07-30 against the real failure** — job **11573880** (guard ON)
      completed and, at the same 16-mode/res-40 modal-multiplexer
      cross-section as the pre-fix baseline (job 11573125, guard N/A), meow's
      guided-mode count drops from 15 (14 distinct + mode 10/11 counted
      twice) to 14 genuinely distinct modes — the duplicate is gone, exactly
      as designed. The layer-4
      self-overlap control (item 2) failed on its first real run: at the modal
      multiplexer's interface 0, 16 modes at res 40, **meow's own modes 10 and
      11 are the same eigenvector** — normalized overlap 1.0000 exactly, row
      completeness 2.0003 (i.e. that mode counted twice). Several other modes
      also show off-diagonal self-overlap of 0.10–0.19 in apparently-pairing
      sets (7↔13 at 0.1886, 9↔15 at 0.1345, 5↔6 at 0.103), which on a
      symmetric MMI cross-section suggests near-degenerate pairs the solver is
      mixing rather than resolving.

      This is the same failure class `FDTDSolver._iter_distinct_bands()` was
      built for (`modal_multiplexer_validation.md` gap 6) — a rank-deficient
      basis, where one of the *N* requested channels is a copy of another and
      whatever mode it displaced is absent entirely. It was fixed on the meep
      side only; nothing guards the meow side, and `EMESolver.solve_section()`
      (the production path) calls the same `compute_modes()`.

      **Consequences, in order of importance:**
      - `EMESolver`'s body bases are rank-deficient wherever this fires. This
        is a **production** defect, not an audit-tooling one.
      - It explains the 0.0000 that item 1 was chasing without needing the
        band-window hypothesis: meep supplies 16 genuinely distinct modes, meow
        supplies 15 distinct + 1 duplicate, so at least one meep mode has no
        partner and its row completeness collapses. The band-window effect may
        also contribute; both are now measurable separately.
      - Every meow-side number in `validation_markdown/` that assumed 16
        independent modes needs re-reading.

      **Fix implemented 2026-07-30; verified on an audit re-run the same day and
      on the production path 2026-07-31 (job 11585847).**
      `EMESolver.compute_modes()` now over-requests by
      `_DUP_MODE_MAX_EXTRA = 8`, rejects any mode whose self-overlap against an
      already-accepted one exceeds `_DUP_MODE_OVERLAP_THRESHOLD = 0.9`, and
      returns `num_modes` genuinely independent modes (raising if it cannot).
      Placed in `compute_modes()` because that is the single chokepoint every
      call path goes through — `solve_section()`,
      `solve_section_length_sweep()`, `solve()`, `section_modes()` — so unlike
      the meep-side guard this covers the **production** path from the start.
      Switchable via `guard_duplicate_modes` (default **True**) so before/after
      can be measured directly.

      The 0.9 threshold is deliberately far above the 0.10–0.19 mutual overlap
      genuinely near-degenerate *pairs* show here: those are real, independent
      modes and must not be discarded. A dedicated test pins that
      (`TestNearDegenerateModesSurvive`) — identical n_eff must never trigger a
      drop, since only field shape decides.

      Status:
      - 7 meow-free unit tests, `tests/test_eme_duplicate_mode_guard.py`
        (skip-and-reach-further, warning log, near-degenerate survival,
        exhaustion → `RuntimeError`, and both switchable-off behaviours). Pass.
      - Confirmed **inert when there is nothing to fix**: 3 modes at res 12
        gives byte-identical n_eff guarded and unguarded, the same sanity
        property the meep-side guard has on the coupler.
      - **Verified against the real failure** (16 modes, res 40, where modes
        10/11 actually duplicate) — job **11573880** (guard ON) vs job
        11573125 (pre-fix baseline, guard N/A): guided count 14 vs 15, the
        duplicate is gone, nothing else about the comparison moved (see
        item 1 below — the pre-existing orphan at meep mode 9 is byte-for-byte
        the same in both runs, confirming it is a *separate* defect from this
        one).
      - **Fires on the production path, confirmed 2026-07-31** — job
        **11585847**, the modal multiplexer's real `modal_multiplexer_hybrid.py`
        run at the shipped res 40, logged `meow mode 9 duplicates
        already-accepted mode 8 (self-overlap 1.000000)` from a body solve
        inside `HybridSolver`, not from an audit script. The "covers the
        production path from the start" claim above is now observed rather than
        argued. See `modal_multiplexer_validation.md` gap 6.
      - **Does not fire at res 60** on that same cross-section — job
        **11586098** (item 13) requested 24 and the guard filtered none of them,
        so the res-40 duplicate looks like a mesh-resolvable near-degeneracy
        rather than an unconditional driver defect. One cross-section only; see
        item 13's third reading.
      - **Results changed**: any EME result where this fires shifts. Goldens
        are already stale from the 2026-07-28 material-model change
        (`README.md` Finding 6), so they need regeneration regardless — but
        that regeneration must now happen *after* this, not before.

      Re-ran item 1's diagnostic (that is what jobs 11573125/11573880 both
      are): the guard fix does **not** close item 1 — the orphan persists
      identically before and after, so it is not a duplicate-eigenmode
      artifact. Still open: whether gap 6's `solve_section()` deferral
      conclusion (`modal_multiplexer_validation.md`) changes now that the EME
      side is guarded at its chokepoint (item 8).

Item 1 is measured; item 2 is done. Both were prerequisites for trusting
anything else here. What item 1's measurement left behind is items 11–13, and
those are now the prerequisites. **Item 13 is closed as of 2026-07-31** (job
11586098 — the guard's over-request is exactly inert, so nothing measured under
it needs re-reading, and what remains there is a cost decision). Items 11 and 12
are still the live prerequisites; the job carrying both hung on its first attempt
(11586046 → refixed and resubmitted as **11588309**), and diagnosing that hang
opened **item 14** — the same non-convergent eigensolve can hang the production
guard, and item 12's band table shows this device's shipped 16-mode basis is
unachievable at res 40 in the first place.

- [~] **1. Resolve the outstanding meep↔meow discrepancy on the modal
      multiplexer.** *(**Measured at scale 2026-07-30 — job 11575516**
      (`mode_spectrum_convergence.py`, 3:36:13 on 8 cpus, COMPLETED). The
      band-window hypothesis below is **half right, and the wrong half is the
      important one**: the mode-8+ disagreement does shrink with resolution,
      but it does not close, and what is left is not a window artifact and not
      fixable with a finer mesh. The residue is re-homed as items 11–13.
      Earlier stages: localized by `diagnose_mode_correspondence.py`, then jobs
      11573125 (pre-fix) / 11573880 (post item-0-fix) which ruled out the
      duplicate-eigenmode bug as its cause.)*

      **The decisive table** — matched resolution only, 16 modes requested at
      the mux's interface 0, cladding index 1.446804, `MATCH_TOL = 0.005`:

      | res | n_agree | guided meep | guided meow | guided orphans (meep \| meow) |
      |---|---|---|---|---|
      | 20 | 13 | 13 | 14 | — \| 1.71134 |
      | 30 | 12 | 13 | 14 | 1.71622 \| 1.75382, 1.48944 |
      | 40 | 12 | 13 | 14 | 1.74775 \| 1.78589, 1.51288 |
      | 60 | **15** | 14 | 15 | — \| 1.79727 |

      Per-engine self-convergence, consecutive resolutions, same metric:

      | | 20→30 | 30→40 | 40→60 | 60→80 |
      |---|---|---|---|---|
      | meep | 1 | 1 | 3 | *(not run)* |
      | meow | 1 | 1 | 2 | **16** |

      **What it establishes, in order of importance:**

      1. **Modes 0–7 were never in question and are not now.** They agree to
         **≤1.6e-4** at every resolution, and are identical to all five printed
         decimals at res 20, 30 and 60. Expected, and worth stating plainly:
         `EMESolver(mode_solver="meep")` and `FDTDSolver` are both MPB on the
         same mesh, so these are not two answers that happen to be close, they
         are the same solve.
      2. **At res 60 the entire disagreement is ONE mode.** meep's 14 guided
         modes are exactly meow's 15 with meow's index 8 (**1.79727**)
         deleted, matching pairwise to ≤4.9e-4. Not a re-ordering, not a
         window offset — a mode present in one engine's basis and absent from
         the other's.
      3. **That mode is the upper member of a genuine near-degenerate pair.**
         meow gives 1.79727 / 1.79197 at res 60 (split 5.3e-3) and 1.79725 /
         1.79694 at res 80 (split **3.1e-4**) — the pair converges *together*,
         and meow's own duplicate guard fired at **neither** resolution, so
         these are two distinct eigenvectors and not one counted twice.
         meep's per-band `get_eigenmode()` returns only one of them, **at
         every resolution tried**. On a modal multiplexer — a device whose
         whole purpose is keeping modes distinct — silently dropping one
         member of a degenerate pair is the worst-shaped defect available.
      4. **The hypothesis's own criterion, applied honestly.** The bulk *does*
         behave like discretization: guided orphans fall 3 (res 30/40) → 1
         (res 20, 60) and `n_agree` rises 12 → 15, so most of the mode-8+
         divergence is degeneracy-resolution improving with mesh. But the
         missing pair member is present at **all four** resolutions, so by
         "n_agree flat with the same orphan persisting => structural" the
         residue is structural. Note the orphan's *n_eff* moves with
         resolution (1.71134 → 1.79727) while its *identity* — meow's index 8,
         upper member of that pair — does not; reading the moving number as a
         moving defect is the trap here.
      5. **First per-engine converged-resolution measurement on SOI_160.**
         meow plateaus between res 60 and 80 (**16/16** modes within 5e-3);
         meep is still moving at 60 (40→60 agrees on only 3/16). So the
         shared `RESOLUTION = 40` is wrong for *both* engines on this device,
         and they do not want the same number — `modal_multiplexer_
         validation.md` gap 2, now measured rather than asserted.
      6. **The two engines' duplicate guards disagree with each other on the
         same nominal mesh.** Both log overlap `1.000000`; neither agrees on
         where. meep skipped 2 bands at res 20/30/40 (11 and 16, then 11 and
         15 twice) and 1 at res 60; meow skipped 1 at res 30 and 40 and
         **none** at 20, 60 or 80. Two guards on the same eigensolver at the
         same cross-section reaching different conclusions is a defect in the
         guards, and it is the direct mechanical cause of the guided-count
         gap in the table above. Items 11–13.

      **What this does not settle:** which driver is right about that pair.
      A simultaneous solve for the lowest *N* bands is the more trustworthy
      instrument for a *count*, which points at meow — but that is an argument
      from construction, not a measurement, and it is item 12's job to make it
      one.

      *Earlier stages of this item, kept because each one narrowed the search
      and the eliminations still hold.* A reduced smoke run (2 modes, res 12)
      **cleared two suspects and convicted a third**:

      - **Geometry is innocent.** `meow(CompoundGeometry mux)` vs
        `meow(mmi1 local frame)` is an exact identity, n_eff equal to 5
        decimals. The `CompoundGeometry` routing the check uses is fine (as
        expected, now that Geometry #6 is fixed).
      - **Half-span is innocent.** 0.05 vs 0.10 is an exact identity.
      - **`meep` and `EMESolver`'s meow agree EXACTLY** — η = 1.0000 in both E
        and H, n_eff identical to 5 decimals (2.49879 both). That is the
        expected result and worth stating plainly: `FDTDSolver`'s
        `get_eigenmode()` and `EMESolver(mode_solver="meep")`'s
        `compute_modes_meep()` are both MPB on the same mesh, so at matched
        resolution they are not merely *similar*, they are the same solve. **So
        the runtime check's cross-engine comparison is not the broken part.**
      - **The audit script is the outlier.** `mw.compute_modes` (its call)
        returns n_eff 2.68730 where `EMESolver` returns 2.49879 on the same
        mesh and cross-section — a different mode solver, hence a different
        answer. Its 0.9735 is a *different measurement* than the check's, not a
        contradicting one, which dissolves the original "one of the two must be
        wrong" framing.

      **The band-window hypothesis, and its verdict (job 11575516).** The
      hypothesis was: a band-window mismatch created by the duplicate-band
      guard. `_iter_distinct_bands()` skips duplicate meep bands, so to return
      *N* distinct modes it reaches past `band_num N` (on this device it
      skipped 11 and 15, so it went to 18), while meow asked for *N* returns
      its own *N* lowest. The two sets are then genuinely different, and meep's
      top modes can sit above meow's window with no partner at all — which
      presents as **row completeness ≈ 0, not as a mis-ordering**.

      **Verdict: the mechanism is real but it is the wrong sign.** The guard
      does move meep's window upward exactly as described, and that is
      confirmed (meep reaches band 18 for 16 modes). But the surviving orphan
      is **not** one of meep's extra top modes lacking a partner — it is a
      *meow* mode that meep's basis does not contain at all (point 2 above), at
      an index in the middle of the guided set rather than at the top. So the
      proposed fix — over-requesting on the meow side to cover meep's band
      reach — would not have closed it: meow already *has* the mode; meep is
      the side that is short. The (b)-completeness framing was right, the
      direction was backwards, and the remedy is on the meep driver (items
      11–13), not on the comparison window.
- [x] **2. Add layer 4 (same-engine controls) as a first-class helper.**
      **DONE 2026-07-30** — `beef.mode_identity.self_overlap_diagnostics()`,
      and wired into `HybridSolver`'s check as a *gate*: both sides' controls
      run first, and a failing one logs `SELF-OVERLAP CONTROL FAILED` and
      **skips** the cross-engine comparison rather than reporting a number
      from broken data. Costs nothing (reuses the modes already solved).
      Already earned its keep — it correctly caught a duplicate eigenmode
      inside one basis while a test of item 3 was being written, which is a
      measurement fault rather than a basis mismatch.

      One honest limitation, pinned by
      `tests/test_mode_identity_overlap.py::test_coarse_sampling_of_orthogonal_modes_still_passes`:
      coarse *sampling* of an already-orthogonal set does **not** fail the
      control (three midpoint samples of `sin(kπy)` preserve the orthogonality
      integral exactly — a quadrature coincidence, verified). Its power against
      an under-resolved mesh comes from the *solve* degrading into
      non-eigenmodes, not from the grid being coarse per se.
- [x] **3. Upgrade layer 1 from diagonal to matrix.** **DONE 2026-07-30** —
      `beef.mode_identity.overlap_diagnostics()` returns per-mode diagonal,
      row/column completeness, off-diagonal mass and margin; `HybridSolver`'s
      failure message now **names the dominant symptom** rather than printing
      one scalar, because the remedies differ:
      `BASIS INCOMPLETE` (a mode unspanned by the other side — truncation or
      band-window mismatch) vs `MIS-ORDERED` (all modes present, wrong
      indices — the contraction is mixing channels) vs `LOW OVERLAP`. It also
      prints the per-mode arrays, since both of this session's false alarms
      were diagnosed from the failure's *shape* across modes. 12 new tests.
- [x] **4. Add the H-only overlap as an independent corroboration** of the
      E-only result. **DONE 2026-07-30** — `beef.mode_identity.field_overlap_matrix()`
      (and `self_overlap_diagnostics()`) gained a `components` parameter,
      default `("Ey", "Ez")`; passing `("Hy", "Hz")` runs the identical
      normalized-overlap machinery on H instead (`_transverse_e_overlap()`
      renamed `_transverse_overlap()` since it is no longer E-specific — the
      module docstring's argument for why a field-type-pure overlap works
      applies verbatim to H against H). Replaces the ad-hoc key-swap dict
      `diagnose_mode_correspondence.py` used to reuse the E-only code path;
      that script now calls `field_overlap_matrix(..., components=("Hy",
      "Hz"))` directly.

      Wired into `HybridSolver._check_interface_correspondence()`: right
      after the E-only `overlap_diagnostics()` call, computes the H-only
      matrix on the same already-sampled fields (`section_mode_fields()`/
      `stack_meow_mode_fields()` always return all six components, so this
      is free — no extra eigensolve) and logs its worst diagonal, worst row
      completeness, and correspondence assignment alongside the E-only
      numbers, plus a WARNING if the two independent measurements disagree
      on which mode pairs with which.

      **Corroboration only, not a second gate** — it does not raise and does
      not affect `mode_correspondence_action`. The E-only comparison is the
      one validated against real devices to date (Layer 1); disagreement
      between the two is a signal that the E-only result — or this
      interface's H-field sampling — may not be trustworthy, not proof
      either way, per the module's still-open E:H-calibration question.
      57 existing `mode_identity`/`hybrid`-correspondence tests pass
      unchanged against the refactor (`_transverse_e_overlap` was private
      and had no external callers).
- [x] **5. Fix the `body_stub` < `dpml` blocker.** **DONE 2026-07-30** — took
      the second option: gave the mode-field entry points a way to query
      without the full PML pad, rather than widening the example constants
      (which would have meant re-recording every affected golden for no
      physics reason). `FDTDSolver._build_section_mode_sim()` (and
      `section_mode_neff()`/`section_mode_fields()` above it) now accept a
      `pad` override, defaulting to `None` → `self.dpml` (behavior-preserving
      for every existing caller). This is safe specifically because
      `get_eigenmode()` never time-steps — the PML is built but never
      exercised for a mode-field query — so the clearance's only real job is
      keeping `geometry.section_prisms()` inside uniform material, which a
      value well under `self.dpml` still does. `HybridSolver` gained
      `mode_correspondence_pad`, threaded to `section_mode_fields()` only for
      solvers whose signature declares a `pad` parameter (duck-typed like
      `neff_override`/`save_dir`), so `EMESolver.section_modes()` (which needs
      no such override at all) is unaffected. A value like 0.1 comfortably
      fits inside the coupler/MMI examples' `body_stub = 0.25`. 12 new tests
      (`tests/test_fdtd_mode_sim_pad.py`,
      `tests/test_hybrid_mode_correspondence.py::TestModeCorrespondencePad`),
      meep/meow-free. Does not itself turn `check_mode_correspondence=True` on
      for the three affected examples — that is a separate decision requiring
      real-device verification, not bundled here.
- [~] **6. Address (b) — `interface_modes` convergence (layer 2).** No device
      has this measured. Do it after item 5, since the sweep needs the check to
      be runnable per device.

      **The reusable check itself: DONE 2026-07-31** —
      `beef.interface_modes_convergence.interface_modes_convergence()`, built
      the same way item 7's `cut_position_invariance()` was: hand it one
      already-solved `SMatrix` per `interface_modes` count plus the counts
      themselves, and it reports, for every `(port, in_mode, out_mode)` key,
      the wavelength-mean transmission at each count and the **step** between
      consecutive counts. Deliberately NOT a fitted slope like layer 3's tool
      — cut-position invariance's correct answer is a horizontal line (the
      S-matrix must not move at all), so a near-zero slope is the right
      criterion; a mode-count sweep is expected to move while the basis is
      incomplete and plateau once it is not, with no reason to be linear in
      between. `last_step` (has it stopped moving at the HIGHEST count tried
      — the literal "does the final S-matrix plateau" question this item
      asks) is the headline metric, `max_step` is informational. Raises if
      the supplied SMatrices don't all carry the same key set, same
      reasoning as layer 3's tool. 7 new tests
      (`tests/test_interface_modes_convergence.py`), meep/meow-free.

      **Not yet done — exercising it on a real device, and now explicitly
      deferred.** Same split item 7 landed in: the generic comparison exists,
      wiring a real sweep onto it needs an actual meep/meow campaign, not a
      code change. `directional_coupler_hybrid_modes_sweep.py`
      (`examples/directional_coupler/archive/` as of 2026-07-31 — see below)
      already sweeps `N_MODES = 2..12` on the coupler (across the full
      `coupling_length` range, saving per-N `.npz` grids) and would have been
      the closest existing candidate, the same relationship
      `directional_coupler_hybrid_stub_sweep.py` had to item 7 before it was
      migrated — not a drop-in fit either way, since it saves raw
      `T_bar`/`T_cross` arrays rather than `SMatrix` objects and grids over
      *two* axes (mode count and coupling length) where this tool wants only
      the mode-count axis varied at one fixed device configuration.

      **Deliberately not picked up now.** All five directional-coupler sweep
      drivers, including this one, moved to
      `examples/directional_coupler/archive/` the same day this tool was
      built — retired pending new, more general sweeping functionality once
      every device is verified (`notes.md` Misc #5). Migrating a real device
      sweep onto this tool is real future work, not blocked by anything here;
      it is simply not this session's task. Whichever device is chosen when
      it is picked up must be resolution-converged first (item 6's own
      precondition, echoing item 1 point 5) — the mux is a poor first choice
      on that count alone, since gap 2 there is still open and res 40 does
      not even supply the requested basis (item 12); the coupler, already
      converged at res 80, is the safer first target, archived driver or a
      fresh one.

      **Naming collision, not a status update:** `verify_item_6_engine_agreement.py`
      and `batch/verification/run_item_6_res_sweep.batch` serve `resolved/
      directional_coupler_validation.md`'s item 6 (engine agreement, closed),
      not this one — that doc's own numbering, unrelated question, unrelated
      device. ~~This item has no tooling yet.~~ The reusable comparison now
      exists (above); no *device* has been exercised through it yet. Reaching
      for "the item 6 script" by name alone still gets the wrong item.

      **Still open, and job 11575516 did NOT close it** — recorded because both
      `mode_spectrum_convergence.py`'s docstring and `running_jobs.md`'s
      snapshot filed that job partly under this item, which overstates it. That
      job swept **resolution** at a fixed `N_MODES = 16` and measured the *mode
      spectrum*. Layer 2 is the other axis entirely: sweep `interface_modes[k]`
      itself and require the final **S-matrix** to plateau. Neither the swept
      quantity nor the observable matches.
      What 11575516 does contribute is a precondition: it says how many of a
      requested 16 are even guided at a given mesh on this device (13–14 meep,
      14–15 meow across res 20–60), so an `interface_modes` sweep here has a
      sane upper end instead of guessing. It also, per item 1 point 5, means
      such a sweep must not be run at res 40 — neither engine is converged
      there, so a plateau found at that mesh would not be a plateau in the
      physics.
- [x] **7. Promote cut-position invariance (layer 3) to a reusable check.**
      **DONE 2026-07-30** — `beef.cut_position_invariance.cut_position_invariance()`.
      `directional_coupler_hybrid_stub_sweep.py`'s own boundary/geometry
      bookkeeping and its single device-specific scalar (the 50:50
      coupling-length crossing) stay exactly where they were — genuinely
      device-specific, not something a generic function should own — but the
      *comparison* is now generic: hand it one already-solved `SMatrix` per
      stub value plus the stub values themselves, and it fits transmission
      (`|coeff|^2`, wavelength-mean) against stub for every
      `(port, in_mode, out_mode)` key, reporting the fitted **slope** — drift
      per µm of stub — as the headline metric, not just an endpoint
      difference. Works for any device/port structure (1 port, 4, N) with no
      per-device code, since it operates purely on `SMatrix`'s own key/coeff
      interface (`SMatrix` gained a `keys()` accessor for this). Raises if
      the supplied SMatrices don't all carry the same key set — a caller
      building inconsistent stubs, not a property of the sweep.
      12 new tests (`tests/test_cut_position_invariance.py`), meep/meow-free
      (pure `SMatrix`/numpy). **Not yet done**: migrating
      `directional_coupler_hybrid_stub_sweep.py` itself (and a similar sweep
      for another device) onto this function, to confirm it reproduces the
      167x drift-collapse result — deferred since that needs a real
      meep/meow run, not a code change.

      **2026-07-31 — explicitly deferred further, not just unstarted.**
      `directional_coupler_hybrid_stub_sweep.py` and every other
      directional-coupler sweep driver moved to
      `examples/directional_coupler/archive/` (`notes.md` Misc #5): new,
      more general sweeping functionality is planned once every device is
      verified, and migrating this one script onto `cut_position_invariance()`
      is no longer the obvious next step ahead of that redesign. The function
      itself is unaffected — it never depended on the script, only on
      already-solved `SMatrix` objects — and the archived script is still
      exactly where a future migration (of it or its replacement) would start
      from.
- [x] **8. Close the `solve_section()` gap for (a) and (c) — the scoped
      zero-cost step.** **DONE 2026-07-30** — added
      `beef.mode_identity.transfer_matrix_duplicate_scan()` (plus
      `excited_mode_axis()` and the shared `warn_on_transfer_duplicate()`
      logging wrapper), wired into both `FDTDSolver.solve_section()` and
      `EMESolver.solve_section()` right before they return `T`. Genuinely
      zero-cost as scoped: `T` is already computed by the time either method
      returns, so this is `O(n_excited^2)` vector arithmetic on data already
      in memory, no extra eigensolve or simulation.

      `solve_section()` excites from whichever side has FEWER modes (its own
      reciprocity rule); a duplicate excited mode shows up as two nearly
      PARALLEL vectors along that side's axis, each spanning the
      wavelength-resolved response measured on the other side — the
      transfer-matrix analogue of `_iter_distinct_bands()`/
      `EMESolver.compute_modes()`'s duplicate-eigenvector guard
      (`modal_multiplexer_validation.md` gap 6), but observed on the matrix
      `HybridSolver` actually contracts rather than a side-channel eigensolve.
      Similarity is complex cosine similarity maxed over wavelength, same
      normalization/threshold (`_DUPLICATE_OVERLAP` = 0.9) as the field-based
      guards.

      **Logged as a WARNING, not raised or gated on.** Unlike the field-shape
      same-engine controls, a high similarity here is a symptom, not proof —
      two excited modes could in principle produce correlated responses for a
      reason other than being the same eigenvector. Left for investigation at
      the call site rather than an automatic hard failure.

      **Scope note, not yet closed**: only the freshly-simulated path is
      scanned — a `save_dir` cache-hit early return in either `solve_section()`
      skips the scan (that `T` was scanned, if at all, when it was first
      computed). `EMESolver.solve_section()`'s scan is downstream of
      `compute_modes()`'s own `guard_duplicate_modes` (item 0), which already
      rejects a duplicate WITHIN the mode basis a section solves on — this is
      a different (if related) object, the *transfer matrix*, and can still
      flag independently. 13 new tests
      (`tests/test_transfer_duplicate_scan.py`), meep/meow-free — exercising
      the scan/threshold/logging logic directly rather than through a full
      FDTD/EME simulation fake, which would be disproportionate to this
      "scoped zero-cost step."
- [x] **9. Drop the unguided mode from every device's `N_MODES`.** **DONE
      2026-07-30, partially** — `validation_markdown/README.md` Finding 1 says
      every interface basis carries at least one radiation mode "at no cost"
      to remove. Today's work sharpened why it matters: those columns have no
      engine-independent shape, cannot be validated by any mechanism here, and
      are still contracted by `solve()`.

      Dropped where the measurement is current and the doc already names it
      free: coupler `N_MODES` 6→5 (band 5, n_eff 1.4367, below the 1.4443
      cladding — `directional_coupler_validation.md` check 2) across
      `directional_coupler_hybrid.py`, `_stub_sweep.py`, `_resolution_sweep.py`
      (plus the stale `N_MODES=6` label in `_modes_sweep.py`'s plot legend);
      2×2 MMI `N_MODES` 15→14 (band 14, n_eff 1.4078 — `mmi_2x2_validation.md`
      gap 3).

      **Not touched, deliberately:**
      - 1×2 MMI — already at `N_MODES = 8`, down from 12 via a *different* fix
        (`mmi_1x2_validation.md` gap 3, degenerate-band ambiguity near cutoff,
        not this finding). Its own comment states bands 0-7 are "comfortably
        guided" — i.e. 8/8 guided at the current setting — so Finding 1's
        stale 12-requested/11-guided row no longer describes it and there is
        no unguided column left to drop.
      - modal multiplexer — the resolution-convergence job (item 1, job
        11575516, now **complete**) confirms the guided count moves with
        resolution rather than sitting at a fixed always-unguided column:
        meep 13/13/13/14 at res 20/30/40/60, meow 14/14/14/15/15 at res
        20/30/40/60/80. So 1–2 of the 16 are unguided depending on mesh, and
        which ones is mesh-dependent. Still premature to pick a number to
        drop — and now for a sharper reason than when this was written: at
        res 40 (the shipped setting) the two engines do not even agree on the
        guided *count*, so a drop chosen there would be chosen against a
        basis that changes as soon as the resolution is fixed
        (`modal_multiplexer_validation.md` gap 2).
      - ring resonator — check 2 is "not applicable" (no `FDTDSolver`/
        `EMESolver` objects built).

      No golden regeneration triggered: `tests/test_golden_examples.py`'s
      directional-coupler replica pins its own self-contained `n_modes = 3`
      at different resolution/calibration settings by design, and no
      mmi_2x2 entry exists in `TARGETS`. Goldens are stale regardless
      (`README.md` Finding 6) and still need a full re-run.
- [x] **10. Retire or re-scope n_eff nearest-neighbour matching.** **DONE
      2026-07-30** — re-scoped, not retired: `assert_modes_correspond()`
      (`tests/test_mode_identity.py`) is still used (its own
      `TestModeCorrespondenceLogic` suite is a genuinely useful fast
      correspondence-logic regression, and its two `@pytest.mark.slow`
      real-device tests are still cheap same-resolution smoke checks worth
      keeping), but its module docstring, its own docstring, and both
      real-device call sites now state plainly that it is a fast smoke check
      only — resolution-locked and shape-blind (it cannot tell two genuinely
      different modes with coincidentally close n_eff from an actual index
      swap) — superseded by `beef.mode_identity`'s field-overlap machinery
      (`field_overlap_matrix()`/`best_correspondence()`/
      `overlap_diagnostics()`, also reachable at runtime via
      `HybridSolver(check_mode_correspondence=True)`) as the authoritative
      mechanism. `beef/mode_identity.py`'s own module docstring points back
      the other way, so the two read consistently regardless of which one a
      reader opens first. No behavior change — docs/docstrings only; the
      existing fast test suite (`tests/test_mode_identity.py -m "not slow"`)
      still passes unchanged.

Items 11–13 are all **opened 2026-07-30 by item 1's measurement** (job
11575516) and are the same subject from three angles: the duplicate-mode guards
added earlier that day are now load-bearing on the production path, and none of
the three properties they need has been checked.

- [x] **11. The meep-side guard's duplicate check runs on a grid that cannot
      resolve the modes it is judging — CLOSED 2026-07-31, not live at res 20
      OR res 40.** `_iter_distinct_bands()`
      (`beef/solvers/fdtd.py`) decides "same eigenvector or not" from a fixed
      `_DUP_BAND_DIAG_NY = 16` × `_DUP_BAND_DIAG_NZ = 8` sample spanning the
      whole **monitor** width. On job 11575516's cross-section that monitor is
      5.0 µm wide, so the transverse pitch is 0.3125 µm and only **8 of the 16
      points fall inside the 2.8 µm core** — against a lobe pitch of ~0.19 µm
      for the 15th mode of that core. The grid is coarser than the structure it
      is asked to distinguish.

      The constant's own comment justifies staying coarse: it "only needs to
      distinguish same eigenvector from different eigenvector (0.9999... vs
      ≤ 0.1257)". That 0.1257 was measured on the **6-mode** coupler. It does
      not carry to 16 modes, and the arithmetic says so: idealized slab
      profiles `sin(mπy/W)` sampled on that exact grid alias mode *m* onto mode
      *18−m* at **η ≈ 0.9988** for every *m* from 1 to 8 — comfortably above
      the 0.9 threshold, i.e. a *genuine, distinct* mode would be discarded as
      a duplicate and silently skipped.

      **Latent, not demonstrated.** Real modes here extend into the cladding
      rather than vanishing at the core edge, so they alias less than the
      idealized model, and the guard did accept bands 12, 13, 14, 17 and 18 on
      this device. Its three skips (bands 11, 15, 16) each coincide with a
      genuine n_eff tie in the pre-guard audit
      (`modal_multiplexer_validation.md` gap 1), so on this cross-section its
      verdicts look correct. The defect is that they are not correct *by
      construction* — there is no discriminating margin left at these mode
      counts, and the failure direction (drop a real mode, reach past it, log
      it as a duplicate) is indistinguishable in the log from a correct skip.
      This is the same constraint the "Two hard constraints" section above
      already states for the cross-engine comparison mesh, applied to the
      guard's own grid; it was not applied there.

      **Measured 2026-07-30**, `examples/verification/python/
      diagnose_dup_band_grid.py`: solved band_num 1–20 unguarded at the mux's
      interface 0, res 20 (one of the two resolutions where job 11575516
      confirmed the shipped guard actually fires), then replayed the guard's
      own accept/skip logic on the shipped 16×8 grid and three finer ones
      (32×8, 64×16, 128×16, up to 72/128 points inside the core vs the shipped
      8/16).

      | grid | skipped band_nums | worst overlap among accepted |
      |---|---|---|
      | 16×8 (shipped) | 11, 16 | 0.3392 |
      | 32×8 | 11, 16 | 0.2265 |
      | 64×16 | 11, 16 | 0.1815 |
      | 128×16 | 11, 16 | 0.1940 |

      **All four grids agree exactly** — same two skips, same accepted set,
      same resulting n_eff to the printed precision — and both skips are
      independently confirmed as genuine: `(10,11)` and `(15,16)` are exact
      n_eff ties (`|Δn_eff| < 1e-9`), a threshold-free check that needs no
      field sampling at all. The idealized-slab aliasing argument (mode *m*
      onto mode *18−m* at η ≈ 0.9988 for `sin(mπy/W)` profiles) does not
      transfer to this device's real fields: real modes extend into the
      cladding rather than vanishing at the core edge, and the actual worst
      margin (0.34 at the shipped grid) sits well below the 0.9 threshold, not
      approaching it.

      **Verdict: not live at res 20.** The shipped grid is not merely
      "correct here by luck of which bands happened to be true duplicates" —
      it agrees with a grid 8× finer on every accepted band's overlap
      structure, not just the two skip decisions. Downgraded from "latent,
      not demonstrated" to "checked and clear" at this resolution. **Not
      re-checked at res 40** (the shipped `RESOLUTION`, and where job
      11575516 separately confirmed guard firings at band_num 11/15) — worth
      a second data point before closing this item outright, since higher
      resolution packs more oscillation into the same core width and could in
      principle narrow the margin further. Re-run the same script with
      `RES = 40.0` to check.

      **CLOSED 2026-07-31 — clear at res 40 too, job 11588309.** The missing
      data point is in. Same script (`BEEF_DIAG_RES=40`), same four grids, at the
      shipped `RESOLUTION`:

      | grid | points in core | skipped band_nums | worst overlap among accepted |
      |---|---|---|---|
      | 16×8 (shipped) | 8/16 | 11, 15 | **0.3877** |
      | 32×8 | 18/32 | 11, 15 | 0.1924 |
      | 64×16 | 36/64 | 11, 15 | 0.1201 |
      | 128×16 | 72/128 | 11, 15 | 0.1209 |

      **All four grids agree exactly again** — same two skips, same 16 accepted
      band_nums, identical resulting n_eff to the printed precision — and both
      skips are independently confirmed as genuine by the threshold-free test:
      `(10,11)` and `(14,15)` are exact n_eff ties (|Δn_eff| < 1e-9). The
      idealized-slab aliasing argument that opened this item does not transfer at
      res 40 either.

      **The margin did narrow with resolution, as predicted, and it is still not
      close.** Worst off-diagonal overlap among accepted bands on the shipped
      grid: **0.3392 at res 20 → 0.3877 at res 40**. That is the direction this
      item warned about — a finer mesh packs more oscillation into the same
      2.8 µm core — but it sits less than half way to the 0.9 threshold, and the
      finer grids' 0.12–0.19 bracket the coupler's 0.1257 reference. Two
      resolutions is enough to close the item as *not live*, while recording that
      the quantity to re-check if this device ever moves to res 60+ is the
      accepted-band margin, not the skip decisions.

      **Note the skip list changed between resolutions:** res 20 skipped 11 and
      **16**, res 40 skips 11 and **15**. Both are exact ties at their own
      resolution, so both are correct — the third duplicate pair simply sits at a
      different band_num on a different mesh. Read the pair positions as
      mesh-dependent, not as fixed properties of the device.

      **Blast radius, now bounded rather than open.** Every recorded
      mode-identity number measured before 2026-07-30 was taken with no guard
      at all; every one taken since (mux, res 20) has this grid underneath it
      and it is now shown to agree with a resolved measurement there. The
      2×2 MMI's "cleanest device audited" result (job 11500536, 15 modes)
      still predates both guards and has never been re-measured with them on
      — `mmi_2x2_validation.md` gap 6 — but that gap is now about the guards
      firing at all on that device, not about this grid being untrustworthy
      when they do.

- [x] **12. DONE 2026-07-31 (library fix landed). Every band's eigensolve starts from the same k guess — the likely
      root cause of the duplicates both guards are working around. TESTED
      2026-07-30, net positive but not a clean fix.**
      `_iter_distinct_bands()` calls `sim.get_eigenmode(..., band_num,
      mp.Vector3(f * _KPOINT_SEED_NEFF, 0, 0), parity=mp.NO_PARITY)` with
      `_KPOINT_SEED_NEFF = 2.5` — the *same* starting vector for band 1 and
      band 18 alike, N independent targeted solves from one guess. Two
      neighbouring targeted solves started from an identical seed converging to
      the same eigenpair is the textbook shape of the failure item 0 and gap 6
      both document, and it explains why the meow side does not have it in the
      same form: `compute_modes_meep(cross_section, num_modes=N+8)` is **one
      simultaneous solve** for the lowest N+8 bands, not N targeted ones.

      That asymmetry is not incidental — the simultaneous side is exactly the
      side that finds the near-degenerate pair member meep misses (item 1
      point 3). Enumeration method, not physics, is the leading explanation for
      the one disagreement resolution does not close.

      **Two things to try, cheapest first.** (a) Ladder the seed: start band
      *k* from band *k−1*'s converged `k_x` instead of the fixed 2.5, so
      successive solves are not launched from the same point. (b) Use one
      simultaneous `compute_modes_meep()` call on the same mesh as the
      reference for the *count*, and check whether meep's targeted ladder
      reproduces it band for band. If (a) removes the duplicates, both guards
      become a backstop rather than load-bearing machinery, and item 11's
      resolution problem stops mattering as much.

      **(a) tested 2026-07-30**, `examples/verification/python/
      diagnose_kpoint_ladder.py`: 20 bands, unguarded, at the SAME
      cross-section and resolution (res 40) as the confirmed pre-guard
      reproduction, seeding fixed (`_KPOINT_SEED_NEFF` every band, today's
      code) against laddered (band *k* seeded from band *k−1*'s own converged
      n_eff, naive — chains off whatever *k−1* converged to, duplicate or
      not).

      | | exact ties | distinct values / 20 | worst overlap |
      |---|---|---|---|
      | fixed (today) | (10,11), (14,15), (19,20) | 17 | 1.0000 (10 vs 11) |
      | laddered | (9,10) | 19 | 1.0000 (9 vs 10) |

      **Net positive, and it recovers real basis rank**: 2 of the 3 raw
      duplicates resolve into genuinely distinct modes under laddering —
      1.512589 at position 14 and 1.357849 at position 19, both **entirely
      absent** from the fixed-seed run (silently replaced by a copy of a
      neighbour). 19 of 20 solves land on distinct eigenvectors, against 17 of
      20 today. Confirms the "different eigenproblem per band" framing this
      item opened with: the fixed seed is not merely inefficient, it is
      **losing modes the ladder recovers**, on the exact cross-section gap 6
      first found this bug on.

      **Also confirms the anticipated risk, exactly where predicted.** The
      docstring's own risk note — seeding band *k* directly at band
      *(k−1)*'s converged value could instead re-find band *(k−1)*'s own
      eigenvector — happened: band 9 (1.763693, unique in the fixed run) and
      band 10 (1.74775 in the fixed run) both converge to **1.763693** under
      laddering, a new duplicate the fixed seed did not have at that position.
      Naive laddering trades duplicate location, and mostly wins, but does not
      close the guard's job.

      **Unplanned finding: the pre-guard reproduction has a third duplicate,
      not two.** `modal_multiplexer_validation.md` gap 6 documents ties at
      (what this script's 1-indexing reports as) positions 10/11 and 14/15 —
      confirmed present again here, same n_eff to the printed precision. The
      (19,20) tie at 1.349097 is new to the record: earlier work solved only
      up to band_num 18 for this device (16 requested + 2 skips observed), one
      short of where this third duplicate sits. Both guards' `_DUP_BAND_MAX_
      EXTRA` / `_DUP_MODE_MAX_EXTRA = 8` headroom easily covers it, so it is
      not a guard-exhaustion risk — but it means the guard is doing more work
      on this device's basis than previously measured.

      **Not yet tried**: a less naive ladder (seed from *(k−1)*'s neff minus a
      small step, since this basis is monotonically descending in n_eff, so as
      to search *past* the previous band rather than land on it) or ladder
      band *k* from the last **accepted** band rather than the raw previous
      solve. Either is a plausible fix for the self-reconvergence failure mode
      just demonstrated. (b) (one simultaneous `compute_modes_meep()` call as
      a count reference) not yet run.

      **Attempt at C and D, 2026-07-31 — job 11586046 HUNG and was cancelled.
      Resubmitted as 11588309 after a script fix. But it produced one result
      before dying, and it is a bigger one than what it was sent to measure.**

      The job asked `diagnose_kpoint_ladder_v2.py` for `N_TARGET = 20` distinct
      modes and never printed a single policy block in 2h21. Reconstructing the
      converged `k_x` per `band_num` from its MPB trace (target f = 0.763359,
      cladding n = 1.446804) shows why:

      | band_num | n_eff | guided? | |
      |---|---|---|---|
      | 1–9 | 2.7064 → 1.7637 | yes | |
      | **10, 11** | **1.7477, 1.7477** | yes | duplicate |
      | 12, 13 | 1.6836, 1.5936 | yes | |
      | **14, 15** | **1.4848, 1.4848** | yes | duplicate |
      | 16–18 | 1.4337 → 1.3877 | **no** | below cladding |
      | **19, 20** | **1.3491, 1.3491** | no | duplicate |
      | 21, 22 | 1.3226, 1.3075 | no | |
      | **23** | — | — | **no k reaches f at all** |

      Bands 1–22 each converged in 3–8 MPB solves (~58 min for 22 attempts, so
      the original ~25 s/solve budget was about right). **Band 23 burned 161
      solves — 59% of the job's entire eigensolve work — without converging**,
      limit-cycling between k = 0.2346 (f = 0.4775) and k = 1.7480 (f = 1.0283)
      either side of the target. Twenty *distinct* modes needs band_num 23
      because of the three duplicate pairs, and band 23 does not exist at this
      frequency. No wall clock would ever have been enough, and MPB runs in C so
      no Python timeout can interrupt it — the only defence is not asking for
      the band.

      **Superseded by the re-run's own measurement — the ceiling is 14, not 13,
      and the difference is the whole point of this item.** The paragraph that
      stood here inferred from the hung job's fixed-seed trace that "only 13
      distinct guided modes exist at res 40", and concluded the 16-mode basis was
      "unachievable by any seeding policy". Job 11588309 falsified the premise:
      three of the four seeding policies find **14**. The inference was drawn
      from the enumeration under test, which is exactly the thing this item says
      cannot be trusted to enumerate.

      Corrected reading: at the mux's interface 0, res 40, **at least 14 distinct
      guided modes exist, and today's fixed seed finds 13.** The device requests
      16, so a shortfall remains — but part of it is a seeding defect that is
      fixable in software, not a property of the spectrum, and that part is
      exactly what item 12 was opened to find.

      **Caveat on "duplicate" for the (10,11) pair specifically.** meep returns
      byte-identical `k_x` for both, but 11575516 found meow resolves that same
      region into a genuine pair with a 3.1e-4 splitting by res 80. Whether
      (10,11) is a solver duplicate or a real near-degenerate pair meep fails to
      separate is precisely item 1's open question, still open. Every count here
      is "distinct as meep resolves them at res 40", not a claim about the
      physics.

      **Script fixed and resubmitted as job 11588309.**
      `diagnose_kpoint_ladder_v2.py` now targets 16 (what the device requests,
      not 20 chosen for comparability with v1's raw sweep — whose bands 16–20
      were unguided anyway) and **stops the search when a raw solve returns
      n_eff ≤ the cladding index**: past that point nothing above is guided, and
      the next band_num is where the k-search stops converging. A physical
      bound, not a timeout. Termination verified offline against this job's own
      measured n_eff sequence with fake `sim`/`mp` objects (all four policies
      stop at band 16, none reaches the non-convergent region, policy A
      reproduces the measured 13 accepted / 16 attempts / skips at 11 and 15).
      The batch script also **runs item 11 first now** — putting item 12 first
      to protect the higher-value result is what cost both items when item 12
      hung.

      **RESULT — job 11588309, COMPLETED in 34:55 (item 12's half took 19:42).
      Laddering is a strict improvement, and the fixed seed is losing a real
      channel.** Modal multiplexer interface 0, res 40, λ = 1.31 µm, target 16
      distinct guided modes, cladding 1.446804:

      | policy | distinct guided found | attempts | skips | stopped |
      |---|---|---|---|---|
      | **A** fixed (today's code) | **13** of 16 | 16 | 2 — bands 11, 15 | guided spectrum exhausted |
      | **B** naive ladder (raw previous) | **14** of 16 | 16 | 1 — band 10 | guided spectrum exhausted |
      | **C** step ladder (raw − 0.005) | **14** of 16 | 16 | 1 — band 10 | guided spectrum exhausted |
      | **D** accepted ladder | **14** of 16 | 16 | 1 — band 10 | guided spectrum exhausted |

      Every policy terminated on the cladding test at band_num 16
      (n_eff 1.433726 ≤ 1.446804), none on the attempts backstop — so no policy
      had duplicates eat the +8 headroom production also relies on, and the
      script's own stop-reason line distinguishes "found 13 because only 13 were
      there" from "found 13 because it ran out of attempts". That distinction is
      why the counts below can be read as physics.

      **B, C and D return A's set plus one mode: n_eff 1.512589**, sitting
      between 1.593585 and 1.484797. It is *absent entirely* from A — not
      reordered, not duplicated, never found. This reproduces v1's recovered
      1.512589 exactly, now under a real accept/skip/retry search rather than a
      flat sweep, so the recovery survives the guard being in the loop.

      **Four readings, in order of what they change:**

      - **The fixed seed loses a guided channel, and the loss is silent.** On the
        device whose entire purpose is keeping 16 channels separate, today's code
        builds an interface basis missing one of them, with no warning: A's two
        skips are both genuine duplicates, correctly rejected, and it still ends
        up a mode short because band 14 *and* band 15 both converge to 1.484797
        under a fixed seed. Nothing in the log says a mode was missed — a skip
        looks the same whether the retry finds something new or re-finds a
        neighbour.
      - **Laddering costs nothing.** 16 attempts either way, and one *fewer*
        skip. This is the "strict improvement: cheaper AND correct" case the
        script was written to look for, so item 12's fix direction is settled:
        seed band *k* from band *k−1*, do not re-seed from a constant.
      - **C and D are unnecessary refinements — B is enough here.** All three
        laddered policies returned byte-identical sets, attempt counts and skip
        lists. The self-reconvergence risk that motivated C and D did fire (all
        three skip band 10, overlap 1.0, raw n_eff 1.763693 = band 9's answer, so
        band 10 re-found band 9), but the existing duplicate guard absorbed it
        and the retry landed on a genuine mode. So the guard and the naive ladder
        together are sufficient on this cross-section; C's step and D's
        accepted-only seed buy nothing measurable. Do not adopt the extra
        complexity on this evidence — one cross-section, and the argument for C/D
        was always about robustness in cases not seen here.
      - **The cross-policy "DIFFERENT modes" flag is an artefact of my own
        comparison, not a finding.** It reports max |Δn_eff| = 0.027793 vs A,
        which is what you get comparing sorted 13- and 14-element prefixes when
        the extra mode shifts the alignment. B ⊃ A exactly. The comparison should
        align on set membership, not position, when the counts differ — noted
        against the script rather than fixed mid-run.

      **What is still open for this item.** The count is 14 of 16 requested, so
      two channels remain unaccounted for even with laddering, and whether those
      are (a) genuinely absent from the spectrum at res 40, (b) more fixed-seed
      losses a different policy would recover, or (c) the (10,11) pair being one
      mode rather than two (item 1) is not settled by this run. Sub-item (b) is
      now the cheap one to test: the same harness at res 60, where 11575516 found
      meow reaches 16/16 and item 13 found no duplicate fires at all.

      **Library fix landed 2026-07-31.** `FDTDSolver._iter_distinct_bands()` now
      seeds attempt *k* from attempt *(k−1)*'s own raw (accepted-or-skipped)
      converged n_eff instead of the fixed `_KPOINT_SEED_NEFF` every time —
      policy B, the naive ladder, since B/C/D were shown byte-identical on the
      cross-section this was measured on and the extra machinery in C/D buys
      nothing. Landed together with item 14's cladding stop below, since both
      touch the same method and the stop condition is what makes the fix safe
      to ship (laddering alone, with no guided-edge bound, would still be
      exposed to the same non-convergent-k hang). 6 new meep-free tests in
      `tests/test_duplicate_band_guard.py::TestSeedLaddering` (+3) and
      `TestGuidedEdgeStop` (+3, item 14's own), pinning the seed mechanism
      directly rather than re-deriving the physics claim (that's what job
      11588309 already measured). Mutation-checked: all 11 tests in the file
      fail against the pre-fix signature. Fast suite 489 → 495.

      **Known consequence, not a regression:** any device whose shipped
      `N_MODES` exceeds what a cross-section can actually supply at its
      resolution will now raise instead of silently padding with unguided
      modes (see item 14 below) — the modal multiplexer is exactly this case
      today (16 requested, 14 reachable) and its own example script will need
      `N_MODES` addressed before it can run again. Deliberately not fixed as
      part of this change: `modal_multiplexer_validation.md` gap 2 is explicit
      that `N_MODES`/`RESOLUTION` need to be decided jointly, not patched to
      whatever today's ceiling happens to be.

- [x] **13. `guard_duplicate_modes=True` changes the solve, not just the
      filter.** `EMESolver.compute_modes()` over-requests
      `num_modes + _DUP_MODE_MAX_EXTRA` (= +8) in a single
      `compute_modes_meep()` call whenever the guard is on — so the guarded and
      unguarded paths hand MPB *different eigenproblems*, and "guard on vs off"
      is not a clean A/B. Item 0 checked inertness at **3 modes, res 12**
      (byte-identical n_eff); that is the regime where it is least likely to
      matter. At the production end — 16 modes becoming 24 on a res-80 mesh —
      neither the inertness nor the cost has been measured, and job 11575516's
      3h36m runtime was dominated by exactly that. Measure once: 16 modes at
      res 60, guard on vs off, compare the low bands and the wall time. If the
      low bands move, every guarded EME number carries an over-request
      artifact; if they don't, `_DUP_MODE_MAX_EXTRA` can be cut on cost
      grounds.

      **MEASURED 2026-07-31 — job 11586098 COMPLETED (1:54:50)**,
      `examples/verification/python/diagnose_guard_over_request.py`. Measured
      with **three** arms rather than the two this item asked for, because
      "guard on vs off" is exactly the unclean comparison described above —
      the flag moves the eigenproblem *and* the filtering together:

      | arm | guard | requests | what it is |
      |---|---|---|---|
      | A | off | 16 | unguarded baseline |
      | B | on | 16 | production default (solves 24, filters to 16) |
      | C | off | 24 | B's eigenproblem, unfiltered |

      **C vs A** isolates the question this item actually turns on — whether
      the over-request *by itself* perturbs the low bands, with no filtering
      involved. **B vs C** then reads off what the guard rejected and what it
      reached further up the spectrum to replace it with. Wall time is taken
      per arm, so an inert C-vs-A result immediately licenses cutting
      `_DUP_MODE_MAX_EXTRA` on the `C/A` cost ratio alone. Run at res 60 for
      the reason item 1 point 5 gives: res 40 is not converged on this device,
      so a difference measured there would not be attributable.

      **Result — the over-request is inert, and it costs 3.2×.** Modal
      multiplexer interface 0 (`mmi1` local x = −7.7100 µm), λ = 1.31 µm,
      res 60, `N_MODES = 16`, `_DUP_MODE_MAX_EXTRA = 8` (so B and C both solve
      24).

      | comparison | result |
      |---|---|
      | **C vs A** — first 16 bands, positionally | max \|Δn_eff\| = **0.000e+00**, every per-band entry 0.000e+00 |
      | **B vs C** — what the filter dropped | B's 16 map to C indices 0…15, worst match residual 0.000e+00, **no C index skipped** |

      | arm | guard | requests | wall |
      |---|---|---|---|
      | A | off | 16 | **934.0 s** |
      | B | on | 16 (solves 24) | **2948.5 s** |
      | C | off | 24 | **2997.7 s** |

      Two readings, and the cost one is the actionable half:

      - **The guard is a pure filter at production scale.** Asking MPB for 24
        instead of 16 does not perturb the low bands *at all* — not "to 1e-6",
        exactly 0.000e+00 across all 16. So no guarded EME number anywhere in
        `validation_markdown/` carries an over-request artifact, and item 0's
        res-12/3-mode inertness check does generalise to 16 modes on a res-60
        mesh. That was the open risk this item was raised for, and it is
        closed.
      - **`_DUP_MODE_MAX_EXTRA = 8` is therefore pure cost.** `B/A = 3.16×`,
        `C/A = 3.21×` — the two agree, which is itself the confirmation that
        B's expense is the eigenproblem and not the filtering (the filter's own
        cost is the 1.6% gap between B and C). Since the extra 8 bands changed
        no answer here, the constant is a free knob: cutting it to the smallest
        headroom that still covers the observed skip counts (≤ 3 on this
        device, per item 12's ladder runs) would take ~2/3 of the wall time off
        every guarded `compute_modes()` call, including the 3h36m of job
        11575516. **Not yet changed** — the skip-count evidence is from one
        device, and item 12 may yet replace the guard with a seeding fix that
        makes the question moot.

      **Third reading, unplanned and relevant elsewhere: at res 60 no
      duplicate fired on this cross-section.** Arm B filtered nothing — B is
      literally C's first 16 — where item 0 found meow's modes 10/11 identical
      at res **40** on this device, and where the production hybrid run (job
      11585847, res 40) does fire the guard (mode 9 duplicating mode 8, see
      `modal_multiplexer_validation.md` gap 6). Consistent with 11575516's
      finding that meow reaches 16/16 guided and plateaus by res 60: the
      near-degenerate pair the guard trips over at res 40 appears to be
      resolved into two distinct eigenvectors by res 60. **Suggestive, not
      established** — this is one cross-section at one wavelength, and B/C's
      24-band list still shows tight non-degenerate clusters higher up
      (1.261133 / 1.260187 / 1.259385 at indices 16–18; 1.158071 / 1.157491 at
      22–23) which the 16-mode production request never reaches. It does mean
      "does the duplicate survive mesh refinement" is now a real question for
      item 12 rather than an assumption, and that the res-40 firings are the
      regime that needs the fix.

- [x] **14. NEW 2026-07-31, FIXED 2026-07-31 — `_iter_distinct_bands()` can hang a production run
      indefinitely, with no exception and no timeout.** Found by job 11586046
      hanging (item 12): when `sim.get_eigenmode()` is asked for a `band_num`
      that has no solution at the target frequency, MPB's k-search does not
      error — it limit-cycles between two k values forever. 161 solves and
      counting on one band, until SLURM killed it.

      The production guard has the same exposure. `_iter_distinct_bands()`
      (`beef/solvers/fdtd.py:1049`) walks `band_num` upward to
      `n_modes + _DUP_BAND_MAX_EXTRA`, and its `RuntimeError` for
      "could not find `n_modes` distinct eigenmodes" can only fire **after** a
      solve returns. A non-converging solve never reaches it, so the bound is
      not a bound: the loop hangs inside meep instead of raising. `parity=
      mp.NO_PARITY`, `_KPOINT_SEED_NEFF` and the retry count are all irrelevant
      to this — it is the `band_num` request itself that has no root.

      **How close this is to firing today.** For the modal multiplexer at res
      40: 16 requested + 3 skips ⇒ the guard already reaches band_num 19, and
      band 23 is where convergence fails. Four band_nums of margin, on a device
      whose duplicate count is only known at two resolutions and whose shipped
      `RESOLUTION` is already below both engines' converged point (gap 2). Any
      device with a shallower guided spectrum, or one more duplicate pair, hangs.
      It has **not** been observed in production — job 11585847 is the only real
      hybrid run of this device on record and it fires the meow-side guard, not
      this one — so this is latent, not demonstrated.

      **Fix, when it is taken:** the same physical bound the diagnostic now
      uses — stop when a raw solve returns `n_eff ≤ cladding index`, since
      nothing above it is guided and the next band_num is the non-convergent
      one. That is strictly better than a larger retry budget, and it fixes a
      second thing at the same time: as of today the guard will happily accept
      **unguided** bands as basis members to reach its count (bands 16–18 at the
      mux's interface 0, n_eff 1.4337–1.3877 against a 1.446804 cladding), which
      is item 9's concern arriving through the guard rather than through a
      device's `N_MODES`. A raise-instead-of-pad decision belongs with it: if
      the cross-section cannot supply `n_modes` guided modes, a loud failure is
      more useful than a silently unguided basis.

      **Deliberately not fixed 2026-07-31 [morning] ... FIXED 2026-07-31
      [same day, later].** The caution above was right at the time — this
      needed its own tests plus a real before/after, not a blind patch — and
      both now exist. `_iter_distinct_bands()` checks every raw solve's
      n_eff against the cross-section's cladding index (computed once per
      caller from `self.material_stack.sio2_index(wl)`, threaded through
      both call sites) *before* the duplicate-overlap sampling, and raises
      immediately if the guided spectrum is exhausted with fewer than
      `n_modes` accepted — the raise-instead-of-pad direction this item
      itself recommended, since a basis silently padded with unguided modes
      (item 9's concern) has a rank that is a lie. Landed together with item
      12's seed laddering, since both touch the same method and the stop
      condition is what makes laddering safe to ship without reopening the
      hang risk. 3 new tests
      (`tests/test_duplicate_band_guard.py::TestGuidedEdgeStop`), covering:
      raises when a raw solve crosses cladding before `n_modes` accepted;
      stays inert when every requested mode is guided; and the cladding
      check fires *before* the duplicate-overlap check would have silently
      skipped the same band and continued past it. All meep-free, mutation-
      checked against the pre-fix signature (see item 12's entry for the
      shared verification).

      **Consequence, not a regression, and worth restating plainly:** this
      raises the bar from "never observed in production" to "will fire the
      next time a cross-section can't supply what its caller requests" —
      which is exactly the modal multiplexer's own shipped `N_MODES=16`
      against the 14 reachable at res 40 (item 12 above). The currently
      in-flight hybrid run (job 11590663) started before this change and is
      unaffected; any *future* resubmission of `modal_multiplexer_hybrid.py`
      will raise until its `N_MODES`/`RESOLUTION` pair is actually decided
      (`modal_multiplexer_validation.md` gap 2) — deliberately not patched
      to 14 as a side effect of this fix, since that doc is explicit the two
      need a joint decision, not a quick unblock.

## Open questions

- **Which is right on the modal multiplexer — the audit or the runtime check?**
  Item 1. Note the audit script is itself unverified in the layer-4 sense.
- **Which driver is right about the near-degenerate pair at the mux's
  interface 0?** meow's simultaneous solve reports two modes (converging to a
  3.1e-4 splitting at res 80); meep's per-band ladder reports one, at every
  resolution from 20 to 60. A count from a simultaneous solve is the more
  trustworthy instrument in principle, but that is an argument from
  construction and item 12 is what would turn it into a measurement.
- **Is an explicit per-engine-pair E:H calibration worth having?** It would
  unlock the true reciprocity overlap (the quantity that actually governs the
  coupling coefficient) instead of a shape proxy. Measurable once on a bare
  strip fundamental. Unclear whether it is stable across geometry/wavelength —
  if it is not, field-type-pure overlaps are the permanent answer, and that is
  worth knowing definitively rather than by default.
- **Does layer 3 have enough dynamic range on a strongly multimode device?**
  It is proven on the coupler, where the observable (L₅₀) is exquisitely
  sensitive to Δβ. On a device whose output is dominated by one strong channel,
  a mis-ordered weak mode may not move the S-matrix measurably.
- **2026-08-09 — the near-degenerate pair at the mux's interface 0, above,
  is real and reproduces at a different wavelength/mode-count.** The modal
  multiplexer's `check_mode_correspondence` (finally exercised at the
  actual swept wavelengths, not just the band midpoint) failed at interface
  0, λ=1.26 µm, N_MODES=14: mode 8's row completeness only 0.009, unmoved
  by raising N_MODES to 15 — inconsistent with truncation, consistent with
  this same meep-ladder-vs-meow-simultaneous-solve resolution gap. See
  `mode_crossing_correspondence_repair_validation.md` (new, this date) for
  the full write-up, a regime taxonomy for what a wavelength-dependent
  index change can mean, and a staged plan — the repair question this
  section's own open items motivate but do not yet answer.
