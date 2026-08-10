# Mode-crossing / near-degenerate correspondence repair — plan

**BUILD PLAN — the mechanism this document specs does not exist yet.**
Cross-device (any `HybridSolver` chain with a multi-mode interface is
exposed), prompted by a real, still-open failure on the modal multiplexer
(`modal_multiplexer_validation.md` gap 7, 2026-08-09 entries) and a direct
question about the general case: *if the correspondence check can detect a
clean index swap between two wavelengths, why doesn't it just fix it?*

Builds on `resolved/fdtd_vs_EME_mode_validation.md` (the correspondence
**checking** mechanism — layers 1-4, `η` matrix, `best_correspondence()`,
row/column completeness, margin — all already built and in active use).
This document is about the **repair** question that mechanism's own
diagnostics make possible but that nothing currently acts on, plus a
classification question that has to be answered *before* any repair
mechanism can be trusted not to make things worse.

## The question that prompted this

`beef/solvers/hybrid.py`'s `check_mode_correspondence` computes, per
interface, an optimal one-to-one assignment (`best_correspondence()`,
Hungarian algorithm on the overlap matrix `η`) between the two adjoining
solvers' mode bases, along with per-mode diagonal, row/column completeness,
and a margin (best overlap minus second-best, in the same row). Confirmed
directly by grep (2026-08-09): **the computed `correspondence` permutation
is never applied anywhere in `solve()`'s actual matrix contraction.** The
check only ever `raise`s (abort the run) or `warn`s (log and proceed with
the still-misaligned contraction) — it is a gate, not a repair.

So: if mode 11 and mode 12 swap rank between two wavelengths — both engines
still see two well-formed, mutually orthogonal eigenvectors, just sorted
differently — today's mechanism can *detect* that (a clean off-diagonal
`η` with high margin, not a collapsed row) but has no way to *fix* it short
of aborting the whole run. This document is about whether and how to close
that gap safely.

## Why "just auto-permute on a clean swap" is not automatically safe

A wavelength-dependent index change can mean three genuinely different
things, and the correct response differs for each. Conflating them — in
particular, building an auto-repair mechanism before knowing which one a
given failure is — risks turning an honest, loud failure into a silent,
wrong one, which is a strictly worse outcome than today's status quo.

### Regime A — true (allowed) crossing

Two modes belong to different symmetry classes of the actual physical
structure (e.g. one strictly even, one strictly odd about a mirror plane
the whole structure respects). Their eigenvectors never hybridize — as a
swept parameter (wavelength here) moves their eigenvalues through each
other, each mode's own shape stays perfectly well-defined the entire time,
including *at* the crossing wavelength itself. What breaks is only a
**sorting convention** (e.g. "list modes by n_eff, descending") that
arbitrarily relabels which index each mode occupies before vs. after the
crossing. This is the case the user's question describes, and it is the
one case where a correspondence-based relabeling is unambiguous at every
wavelength, including exactly at the crossing.

### Regime B — avoided crossing (anticrossing)

Two modes share a symmetry class and *can* couple. As their eigenvalues
approach, the modes hybridize — each becomes a nontrivial superposition of
the two "bare" characters — and the eigenvalue traces repel rather than
cross. In a window around the point of closest approach, there is **no
correct single-mode assignment**: `best_correspondence()`'s own `margin`
output is exactly the signal for this (small margin = "genuinely ambiguous
pairing... a true physical degeneracy," per its own docstring). Forcing a
hard 1:1 relabel inside that window is wrong regardless of which assignment
is chosen, because the physically correct description at that point is
"the light is coupling continuously between both channels," not "channel A
is now channel B." Away from the window on either side, each mode rapidly
reacquires its bare character and a relabeling is safe again — the
practical question is how wide the ambiguous window is and how to detect
its edges (the margin threshold) reliably.

### Regime C — cross-engine resolution failure (no clean pair to relabel at all)

Independent of whether the underlying physical situation is A or B, one
engine's specific eigenmode-extraction *method* may simply fail to produce
two independent eigenvectors for a close pair in the first place. **This
has direct, already-recorded prior evidence on this exact device, at this
exact interface**, predating today's investigation:
`resolved/fdtd_vs_EME_mode_validation.md`'s own "Open questions" section
and item 12: at the modal multiplexer's interface 0, meep's per-band ladder
search (`sim.get_eigenmode()` called sequentially per `band_num`) returned
a **byte-identical `k_x`** for a close pair of bands, while meow's
simultaneous MPB solve resolved the same region as a genuine pair with only
a 3.1e-4 n_eff splitting (converged at res 80) — an open question the file
itself states is "still open." A per-band ladder is not the right tool for
eigenvalues this close together; there is no permutation that fixes this,
because meep never produced a second, independent eigenvector to permute
onto — the ladder converged onto (effectively) the same vector twice, or
skipped one entirely.

## What today's actual failure (mode 8, interface 0, λ=1.26 µm) most likely is

Three pieces of evidence, all measured today (`modal_multiplexer_validation.md`
gap 7, 2026-08-09 entries), point at regime C rather than A:

1. **The failure is surgically isolated, not a gradient.** Per-mode diagonal
   at N_MODES=14: `[1,1,1,1,1,1,1,1,0.0,1,1,1,0.9999,0.9999]` — mode 8 is
   the *only* zero, sitting in the middle of an otherwise perfectly clean
   spectrum, including the higher-index (weaker-confined, more
   marginal/near-cutoff) modes 12-13, which are fine. Generic
   under-resolution or near-cutoff solver disagreement would be expected to
   degrade gradually from the weakly-guided end inward — this instead looks
   exactly like a genuine, localized degeneracy specific to one mode's
   neighborhood, which is exactly item 12's already-flagged concern.
2. **Raising the basis size did not help.** N_MODES=15: row completeness for
   mode 8 moved from 0.0090 to only 0.0391 — still catastrophic, while
   nothing else changed. If this were basis truncation (a mode present but
   outside the requested window), adding one more slot should have helped
   substantially. It did not, which is inconsistent with a truncation
   explanation and consistent with "the mode meep calls index 8 does not
   correspond to a real, independent eigenvector meow can match at all,"
   i.e. a resolution failure, not a window-size failure.
3. **This interface already has a documented, unresolved near-degenerate
   pair on record**, independent of today's work — `fdtd_vs_EME_mode_validation.md`'s
   open question about interface 0's (10,11)-type pair (numbering there was
   at N_MODES=16, wl=1.31; today's failure is at N_MODES=14, wl=1.26 — not
   proven to be the literal same pair, since which two bands sit closest
   can itself shift with wavelength, but it is the same interface and the
   same general phenomenon: meep's ladder losing resolution against a
   near-degenerate pair that meow resolves cleanly).

**Working hypothesis, not yet independently confirmed: today's failure is a
wavelength-shifted instance of the same already-open meep-ladder resolution
problem, i.e. regime C.** This matters directly for the auto-permute
proposal from the prior message in this conversation: if that hypothesis is
right, **a permutation-based repair would not fix today's actual problem at
all**, because meep never produces a valid second eigenvector to permute
onto in the first place. Building that mechanism first, without this
classification step, risks either wasted effort (it doesn't address this
failure) or worse — silently forcing some permutation onto data that
doesn't support one, converting today's honest `raise` into a confidently
wrong number.

## Staged plan

### Stage 0 — classify before repairing (cheap, do first)

Determine which regime today's failure actually is, reusing tooling that
already exists rather than building new instruments:

- Re-run (or adapt) `examples/verification/python/diagnose_kpoint_ladder_v2.py`
  (built for item 12, currently only ever exercised at wl=1.31) at **wl=1.26
  µm** specifically, at this device's interface 0 — check directly whether
  meep's per-band ladder returns a byte-identical (or near-identical) `k_x`
  for two adjacent bands near index 8/9, the same signature item 12 already
  found at 1.31.
- Run meow's own mode solve at the same interface/wavelength and check
  whether it reports a genuinely resolved close pair near that index (small
  but nonzero n_eff splitting, two independent, mutually orthogonal
  eigenvectors) — the same comparison item 1/item 12 already made at 1.31,
  just at the actual failing wavelength this time.
- **If meep shows a duplicate/near-duplicate and meow shows a real
  resolved pair → regime C confirmed.** The fix is not permutation-based;
  see Stage 2.
- **If meep instead shows two genuinely distinct, orthogonal eigenvectors
  that are simply sorted differently from meow's own ordering → regime A.**
  The permutation-based repair (Stage 1) is the right and sufficient fix.
- Either way, run the same check at a few more wavelengths across the swept
  band (not just 1.26/1.31/1.36) to bound how wide the affected region is —
  relevant to Stage 3 regardless of which regime this turns out to be.

**Confound to rule out first, or at least alongside the above: transverse
domain size.** Checked directly (2026-08-09) — the mode-solve's transverse
(y/z) domain (`sim_width`/`sim_height` for FDTD, `oxide_width`/
`oxide_height` for EME, both 5.0/1.5 µm for this device) has **never been
convergence-tested for any mode-solving context in this repo**, for this
device or any other. The only domain-size sweep that exists anywhere here
is `examples/verification/python/hcg_reflection_convergence.py`'s
`kind="domain_y"` points, applied to HCG's own full-device reflection
measurement, never to an interface mode-solve. `sim_width=5.0` for this
device traces back to "matches modal_multiplexer_hybrid.py's own
oxide_width" (a comment in `modal_multiplexer_fdtd.py`) — i.e. carried by
convention/consistency between the two engines, not independently
justified. An under-sized domain would truncate/reflect a mode's evanescent
tail off the nearby PML, perturbing n_eff and shape for *both* engines
(possibly by different amounts, since each has its own padding — see
`section_mode_fields()`'s own `pad` parameter, defaulting to `self.dpml`)
— which could distort the very A/B/C classification above, not just add a
fourth candidate explanation. **Run this check before trusting a Stage 0
regime conclusion, not after.**

- `examples/verification/python/diag_mode8_domain_sweep.py` (new, this
  date, TEMPORARY): reuses the real script's construction, then mutates
  the already-built `fdtd_solver`/`body_solver` instances'
  `sim_width`/`sim_height`/`oxide_width`/`oxide_height` directly between
  checks (confirmed by reading `_build_section_mode_sim()`/`EMESolver`'s
  own methods that these are read at call time, not baked in at
  construction — no need to rebuild solvers or geometry per domain size).
  Sweeps interface 0, wl=1.26 µm only (the one already-confirmed failure)
  across `(sim_width, sim_height)` = (5.0, 1.5) [shipped baseline], (7.0,
  2.0), (9.0, 2.5) — a few points, matching HCG's own sweep style, not a
  dense grid; scope kept to one interface/wavelength since cost grows with
  domain area at fixed resolution.
- **If mode 8's row completeness improves with a bigger domain**: real
  evidence of transverse truncation as a contributing or causal factor —
  re-run the meep-ladder/meow-simultaneous comparison above at the
  converged domain size before concluding anything about regime A vs C,
  since the "byte-identical k_x" signature itself could conceivably be a
  domain-truncation artifact rather than a genuine degeneracy.
- **If it stays flat across all three domain sizes**: rules domain
  truncation out, and strengthens the regime-C (cross-engine resolution
  failure) hypothesis on its own terms.
- **✅ MEASURED 2026-08-09 (jobs 11769172, then 11769720 after fixing a
  guided-mode-count crash — see below) — domain size matters, and the
  earlier "surgically isolated, therefore a real degeneracy" argument was
  wrong, or at least premature.**

  First attempt (job 11769172, shipped N_MODES=15) crashed going from
  sim_width=5.0→7.0: `FDTDSolver`'s own guided-mode guard found only **14**
  guided modes at 7.0/2.0, fewer than the 15 found at the 5.0/1.5 baseline
  — a bigger domain finding *fewer* guided modes is itself the signature
  of a finite-domain confinement artifact (a marginal near-cutoff "mode"
  that only looked bound because the nearby PML was propping it up).
  Resubmitted (job 11769720) at a lower, domain-size-robust n_modes=12 to
  get past that guard and actually test correspondence.

  | domain (sim_width, sim_height) | mode 8 | new failure |
  |---|---|---|
  | (5.0, 1.5) — shipped baseline | **diagonal 0.0, row completeness 0.009 (FAILS)** | — |
  | (7.0, 2.0) | **diagonal 1.0, row completeness 1.0 (CLEAN)** | mode 11 (the highest/last in the 12-mode basis): row completeness 0.030, margin −0.14 — same basis-incomplete signature, different mode |
  | (9.0, 2.5) | **diagonal 1.0, row completeness 1.0 (CLEAN)** | **none — interface FULLY clean, exact identity correspondence [0..11], worst diagonal 1.0000, worst row completeness 1.0000** |

  Third point completed 2026-08-09 (job 11770757, resubmitted at -t 240 /
  8 threads / trimmed to just this domain size after 11769720 timed out —
  8 CPU threads offered since this script is a single un-MPI'd process and
  more ranks would not apply, genuinely untested whether the build honors
  OpenMP/BLAS threading here). **Result: completely clean at the largest
  domain tested — no failure anywhere in the 12-mode basis.**

  **Reading the full three-point trend**: this is convergent, not noise.
  Mode 8's failure disappears at 7.0 and stays gone at 9.0; the new
  failure that appeared at 7.0 (mode 11, the marginal top-of-basis one)
  is also gone by 9.0. **The interface-0 correspondence failure this
  entire investigation was chasing is very likely a domain-truncation
  artifact, not a real physical mode crossing/degeneracy (regime A/B) or
  a cross-engine resolution limitation (regime C).** The shipped
  `sim_width=5.0`/`sim_height=1.5` — carried over by convention, never
  independently justified (see above) — appears to simply be too small
  for this cross-section's more weakly-confined modes.

  **This substantially reframes Stage 0.** The regime-A/B/C taxonomy
  above remains valid in general (a real crossing/degeneracy could still
  exist and still matters for other devices/interfaces), but for THIS
  specific, already-investigated failure, domain size now looks like the
  primary and possibly sole driver — not the meep-ladder-vs-meow
  resolution question the taxonomy was originally built to chase.
  Likely implication: the earlier N_MODES=16 crash (only 15 guided at
  5.0/1.5) and the N_MODES=15 result (mode 8 still failing) may
  themselves have been domain-truncation artifacts, not genuine ceilings
  on this cross-section's real guided-mode count — a larger, adequate
  domain may support 15 modes, or more, cleanly. **Not yet confirmed**:
  this sweep used a reduced n_modes=12 to dodge the earlier crash: the
  natural next check is the shipped N_MODES (15) at the largest domain
  (9.0, 2.5) — if that is ALSO clean, domain size is confirmed as the
  fix and the regime classification work above becomes unnecessary for
  this device's actual production configuration, not just abandoned
  half-finished. Also still unconfirmed: whether 9.0/2.5 itself is truly
  converged (only one point beyond where things went clean, not a proper
  plateau the way `hcg_reflection_convergence.py`'s multi-point sweep
  established one) and whether interfaces 1-5 / other wavelengths show
  the same domain sensitivity.

  **✅ MEASURED 2026-08-09, same day — confirmed at the shipped N_MODES,
  and the investigation resolves cleanly.** Job 11772198 (1h42m): shipped
  `N_MODES=15`, domain (9.0, 2.5), interface 0, wl=1.26 µm. **Modes 0-13
  (14 of the 15 requested) are essentially perfect** — diagonal 1.0 (two
  at 0.9999/0.9997), row completeness ~1.0 across the board, including
  the ORIGINAL mode 8 this whole investigation started from. **Only mode
  14 — the highest-index, most marginal mode in the 15-requested basis —
  fails**, same signature (row completeness 0.02), matching the
  "top-of-basis" pattern already seen at 7.0/2.0 with a 12-mode request
  (mode 11 failed there). This is not a new mystery: the mode nearest
  cutoff is always the hardest to resolve, at any domain size — it is not
  evidence against the domain-size conclusion, it is the expected residual
  at the edge of whatever count is requested.

  **Conclusion: domain-size truncation was the real, essentially complete
  explanation for this whole investigation.** The original mode-8 failure,
  the N_MODES=16 crash, and N_MODES=15's unmoved mode-8 result earlier
  today were all measured on a too-small domain
  (`sim_width=5.0`/`sim_height=1.5`, never independently justified) and
  should be considered unreliable, not genuine ceilings on this
  cross-section's guided-mode count. **The concrete, well-evidenced fix
  for the production script**: enlarge the domain (9.0/2.5 tested clean;
  a real multi-point convergence sweep, not just this single confirming
  point, should still pick the final shipped value — see the open item
  below) and revert `N_MODES` 15→**14** (not 15) — 14 is exactly the
  count that comes out perfectly clean at this domain size, matching the
  device's own pre-2026-08-09 shipped value, now for a *correct* reason
  (a converged domain) rather than the earlier, since-superseded reason
  (gap 2's small-domain guided-count measurement).

  **Not yet done, deliberately, pending explicit sign-off** (this changes
  the actual production script, not just a diagnostic, and invalidates
  every `solve_section()` cache entry computed at the old 5.0/1.5 domain
  — a real, nontrivial recomputation cost, not a free edit):
  - Change `modal_multiplexer_hybrid.py`'s `sim_width`/`sim_height` (both
    `fdtd_solver` and `body_solver`/`oxide_width`/`oxide_height`
    constructor args) from 5.0/1.5 to the converged value, and `N_MODES`
    back to 14.
  - A real multi-point domain-size convergence sweep (matching
    `hcg_reflection_convergence.py`'s own methodology) before picking the
    final shipped value — this investigation confirms 5.0 is too small
    and 9.0 is at least good enough for interface 0/wl=1.26/N_MODES≤14,
    but does not establish where the true plateau is, nor whether an even
    smaller value than 9.0 (e.g. 7.0, already shown clean up to mode 10)
    might suffice more cheaply.
  - Re-check interfaces 1-5 and the other swept wavelengths (1.31, 1.36)
    at the new domain size — this entire sweep only ever tested interface
    0 at wl=1.26 µm; nothing about the other 5 interfaces or 2 other
    wavelengths has been re-verified at a larger domain.
  - Once a domain size is settled, invalidate the relevant stale
    `example_sim_cache/` entries (config-hash-keyed on `sim_width`/
    `sim_height` themselves, so a real config change will naturally miss
    cache and recompute — no manual deletion needed here, unlike the
    normalization-fix cache invalidation earlier this session) and re-run
    the real hybrid smoke test to confirm the S-matrix itself changes
    sensibly (or does not materially change, if the earlier bad modes
    carried negligible power).

### Stage 1 — margin-gated auto-permutation in `HybridSolver` (build only if Stage 0 finds a real regime-A case)

- Extend `HybridSolver` so that, when `check_mode_correspondence` finds a
  non-identity `correspondence` that is a **clean permutation** — every
  row's diagonal-after-permutation near 1.0, every row/column completeness
  near 1.0 (no collapse), and every margin above a confidence threshold
  (reuse the `mode_correspondence_min_overlap`-style knob already on the
  class) — it permutes the affected section's matrix rows/columns to
  realign indices before contracting, instead of raising, and logs what it
  did (never silent).
- When any of those conditions fail (row/col collapse → regime C-like
  truncation; low margin → regime B-like ambiguity), fall back to today's
  behavior (`raise`/`warn`) — auto-repair must never attempt to paper over
  a case its own diagnostics say is not a clean swap.
- Needs unit tests with **synthetic** overlap matrices for both the
  positive case (clean permutation, should auto-fix) and both negative
  cases (a collapsed row; a low-margin near-degenerate pair) confirming the
  fallback path triggers correctly — this is exactly the kind of case a
  real device is expensive to reproduce on demand, so the safety logic
  needs to be provable without a SLURM job.

### Stage 2 — coupled-subspace handling (only if Stage 0 finds this device's case is regime B or C for real)

A genuinely bigger, riskier architecture change, out of scope to attempt
without its own separate design discussion:

- For regime C specifically: give `FDTDSolver`'s eigenmode extraction a
  fallback for closely-spaced pairs — e.g. a seeded/simultaneous solve
  analogous to what item 12 already built for the *counting* problem
  (`diagnose_kpoint_ladder_v2.py`'s laddering fix), adapted to actually
  *extract* both eigenvectors of a close pair rather than just detect that
  the ladder undercounts.
- For regime B (and as a more general answer for either regime, if
  per-mode extraction turns out not to be reliably fixable): stop
  insisting that every interface decomposes into independent scalar
  channels. Propagate a small coupled 2×2 (or n×n) block through the chain
  for any pair/group whose margin stays low across the wavelengths that
  matter, rather than forcing two single-index channels that do not
  correspond to well-defined physical entities there.
- Not attempted this session. Flagging it here so it is not lost, and so a
  future push at N_MODES or at the permutation mechanism does not get
  mistaken for having addressed this.

### Stage 3 — exercise `check_mode_correspondence` across the actual swept band, not just its midpoint (cheap, valuable regardless of regime)

- `_check_all_interfaces()` currently checks **only** at
  `wl_center = mean(wavelengths)` — for the modal multiplexer's swept
  `[1.26, 1.36]`, that is 1.31 µm, neither endpoint, and today's entire
  investigation only happened because a one-off diagnostic script checked
  1.26 µm by hand. Under normal use (`check_mode_correspondence=True` in a
  real production script) this failure would have been **invisible** —
  the built-in check would have reported everything fine at 1.31 and the
  actual swept wavelengths would never have been examined at all.
- Fix: when `check_mode_correspondence=True`, check at (at minimum) every
  wavelength actually being solved, or at a caller-configurable subsample
  — not a single fixed center point. This is cheap (`section_mode_fields`/
  `section_modes` are "seconds"-scale, no time-stepping, per
  `audit_interface_modes.py`'s own docstring) relative to the real
  `solve_section()` cost it gates, and it is the only way this class of
  band-edge-specific failure is ever caught by default rather than by
  accident.
- This stage is independently worth doing regardless of how Stages 0-2
  resolve, and probably belongs in the library's own default behavior
  going forward, not just as a one-off diagnostic pattern.

## Explicit non-goals / what not to do

- **Do not build Stage 1 before Stage 0.** The whole point of this
  document is that today's evidence favors regime C, where a permutation
  fix is a no-op at best and a silent-corruption risk at worst.
- **Do not treat `mode_correspondence_action="warn"` as safe for a real
  production run.** It exists for diagnostics like today's (report
  everything instead of stopping at the first failure) — used in a real
  `solve()` call, it lets a known-broken contraction proceed. `"raise"`
  remains the only safe default.
- **Do not assume raising `N_MODES` is a general answer to a future
  correspondence failure.** Measured twice, directly, on this device
  today: N_MODES=16 is physically impossible at 1.26 µm (only 15 modes
  guided) and N_MODES=15 did not move mode 8's correspondence at all.
  Basis size and mode-identity correspondence are different failure
  classes with different fixes.
- **Do not assume interfaces or wavelengths not yet checked are clean.**
  Interfaces 1-5 and wl=1.31/1.36 remain completely unchecked as of this
  writing (`modal_multiplexer_validation.md` gap 7's 2026-08-09 entries).

## Cross-device relevance

Any `HybridSolver`-based example with a multi-mode interface is exposed to
regimes A-C in principle — this is a library-level mechanism gap, not a
modal-multiplexer-specific one. The modal multiplexer is simply the first
device where a real instance has been caught, because it is the first
device in this repo with (a) a deliberately near-degenerate design element
(the top/bot branch width asymmetry, `modal_multiplexer_fdtd.py`
`BRANCH_SPECS`) and (b) the first real exercise of
`check_mode_correspondence` away from a single fixed wavelength. Any other
device with closely-spaced output ports, intentional near-degenerate
routing, or a design wavelength range spanning a mode-count change (a new
mode turning on/off within the swept band) is a candidate for the same
class of failure and should not be assumed safe by default.
