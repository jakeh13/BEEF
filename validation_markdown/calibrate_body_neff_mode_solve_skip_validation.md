# Skipping the EME mode solve when `calibrate_body_neff` makes it moot — plan

**BUILD PLAN — the mechanism this document specs does not exist yet.**
Cross-device (any `HybridSolver` body solved with `calibrate_body_neff=True`
and a single-cell EME section — the shape every shipped example's body
actually is). Prompted by a direct question during manuscript-writing,
2026-08-11: if `calibrate_body_neff` overrides the body's propagation
constant with an FDTD-measured β, is the EME mode solve that produced the
(discarded) MPB β doing anything the final transfer matrix still needs?

## The question that prompted this

`HybridSolver(calibrate_body_neff=True)` measures β from a neighbouring
`FDTDSolver.section_mode_beta()` and hands it to `EMESolver.solve_section()`
as `neff_override` (`hybrid.py:1250-1282`). Inside `solve_section()`
(`eme.py:387-535`), the body's modes are still fully solved
(`self.compute_modes(cs, n_modes)`, `eme.py:504`) before `_override_neff()`
replaces each mode's `.neff` and `mw.compute_s_matrix()` builds the transfer
matrix from the (now-recalibrated) mode list. So the expensive part —
solving the eigenmode problem for the field profile — happens regardless of
whether calibration is on, and its one physically load-bearing output
(β/`neff`) gets thrown away and replaced immediately after.

## Why the field profile is not needed, for a single-cell body specifically

Checked directly against `mw.compute_s_matrix()`'s call site
(`eme.py:497,507`): with `Ls=[length]` (one cell — every shipped example's
body), left and right ports are the *same* cross-section's mode basis.
Distinct eigenmodes of one cross-section are orthogonal by construction, so
the returned `T` is diagonal, `T[i,i] = exp(2j*pi*neff_i/wl*L)`, with no
off-diagonal terms — a value determined entirely by `neff`, never by the
field's spatial shape. `neff_override`'s own docstring already says this in
passing ("leaving field profiles ... untouched",
`eme.py:409-411`) but the codebase has never acted on it: the field is
still fully solved and then never read for anything the returned `T`
depends on.

## Checked, not assumed: does anything else quietly depend on that solve?

- **`check_mode_correspondence`'s field fetch is independent.**
  `EMESolver.section_modes()` (`eme.py:337-385`) — what
  `check_mode_correspondence` actually calls for field data
  (`hybrid.py`'s own integration doc) — does its own `mw.create_cells()`/
  `mw.CrossSection.from_cell()`/`compute_modes()` from scratch. It does not
  read anything `solve_section()` computed. **Confirmed by reading both
  methods side by side, 2026-08-11: skipping fields inside `solve_section()`
  cannot break the correspondence check.** This was the main risk worth
  ruling out before writing the rest of this plan, and it rules out clean.
- **Mode count is not discovered by the solve.** `n_modes` is passed in by
  the caller (from `interface_modes`, fixed per device — e.g. the 2×2 MMI's
  `N_MODES=14`), not derived from the solve result. A fast path does not
  lose any "how many modes are there" auto-detection, because that
  detection is not currently happening per-call anyway.
- **Not yet checked**: `warn_on_transfer_duplicate()` (`eme.py:524-526`,
  called on the `T` `solve_section()` just built) currently flags
  near-identical *rows/columns* of `T` as a proxy for near-degenerate
  eigenvectors. On a directly-built diagonal `T`, two entries can look
  "duplicate" only if two modes' calibrated phases happen to coincide —
  a different, weaker signal than the eigenvector check it currently
  stands in for. Needs a decision (see step 5 below), not yet made.

## Scope — where this plan does *not* claim to apply

- `num_cells > 1` (a body actually split into several EME cells with
  differing cross-sections, e.g. a real taper modeled as multiple cells) —
  `mw.compute_s_matrix()` genuinely needs inter-cell field overlaps there.
  No shipped example's body currently does this, but the fast path must not
  silently apply to one that does.
- `calibrate_body_neff=False` — without an override, `neff` itself is the
  thing being solved for; there is nothing to substitute.
- Any future caller that reads `solve_section()`'s modes for something
  besides `T` and the correspondence check (none currently does, per the
  independence finding above — worth re-checking if that changes).

## Implementation steps

1. **Measure first — this whole plan is only worth it if the split is
   real.** Time `compute_modes()` versus `mw.compute_s_matrix()` separately
   for a representative single-cell body (e.g. the 2×2 MMI's `N_MODES=14`
   body at its production `res`). Not measured yet. If `compute_modes()` is
   not the dominant cost, this optimization is not worth the complexity
   below and this plan should stop here.
2. Add a fast path in `EMESolver.solve_section()`, gated on `len(cells) == 1
   and neff_override is not None`: build `T` directly as
   `diag(exp(2j*pi*neff_override[wl_idx]/wl*length))` per wavelength,
   skipping `mw.create_cells()`/`CrossSection.from_cell()`/
   `compute_modes()`/`mw.compute_s_matrix()` entirely on that path.
3. Decide `warn_on_transfer_duplicate()`'s fate on the fast path (see the
   "not yet checked" item above) — keep it running against the
   directly-built diagonal (weaker signal, still free), skip it explicitly
   with a one-line reason, or repurpose it to compare calibrated β values
   instead of matrix rows. Pick one and record why.
4. Cache-fingerprint check: `build_run_config()`'s `_METHOD_ATTRS` allowlist
   for `("EMESolver", "solve_section")` (`cache/config.py`) must still
   produce the same fingerprint for a fast-path run and a full-solve run of
   the *same* physical config, since they must be numerically identical and
   therefore cache-interchangeable — confirm this holds rather than
   assuming it, since a cache key that accidentally depends on which path
   ran would silently double the cache's footprint for no reason.
5. **Numerical equivalence test, before this ships as a default.** Run both
   paths on a real body (2×2 MMI, `N_MODES=14`) and assert the two `T`
   matrices agree to floating-point precision — same standard this repo
   already holds itself to elsewhere (cf. the HCG Fabry-Pérot `Circuit`
   assembly vs. closed-form cross-check, `max|delta| ~1e-16`). Only after
   this passes does the fast path become the default; keep the full-solve
   path reachable (e.g. a debug/force flag) in case a future finding needs
   real fields here after all.

## Open questions

- Does this shortcut generalize to a `num_cells > 1` body where every cell
  happens to share the *same* cross-section (split only for reporting, not
  because the geometry varies)? Not assumed — out of scope until checked
  directly against `mw.compute_s_matrix()`'s behavior in that case.
- Is there a cheaper meow/MPB call that returns eigenvalues without the
  full eigenvector solve, in case a future use of this path still wants
  *some* sanity confirmation of mode count/ordering rather than trusting
  `n_modes` blindly? Not needed for the plan above (mode count is already
  externally supplied, not discovered), but worth knowing if scope grows.
