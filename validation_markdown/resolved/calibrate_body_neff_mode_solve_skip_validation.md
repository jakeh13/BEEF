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
   body at its production `res`). **In progress, 2026-08-11** — `cProfile`
   run on a cache-disabled `EMESolver.solve_section()` call, real MMI body
   geometry (`mm_length=28.0`, `N_MODES=14`, `res=40`), single wavelength
   (1.55 µm). Local run via the `jmh` conda env (has `meow`; the repo's
   ambient `python3` does not). Script:
   `/private/tmp/claude-501/-Users-jhiesener-Documents-Research-duty/a6038aa9-a350-4c1d-a396-b023e0aad6c8/scratchpad/time_solve_section.py`.
   Backgrounded (job `bv08mvjr2`). **Result, 2026-08-11:**
   `compute_modes()` (→ `_solve_modes_raw` → meow's `compute_modes_meep` →
   `meep._get_eigenmode`) = 651.9 s cumulative; `solve_section()` total =
   653.9 s. `compute_s_matrix` never showed up as its own line in the
   profile — everything besides the mode solve (mesh/cell build, S-matrix
   build, overhead) is ≈2 s combined, i.e. ~0.3% of the call. **Split is
   real, dominance confirmed — proceed with step 2 (already written
   below) and step 5.** Ran fully local, no cluster (1 MPI process, `jmh`
   conda env) — just slow, ~11 min for one wavelength × 14 modes at
   `res=40`, so budget for that when running step 5's equivalence check.
2. Add a fast path in `EMESolver.solve_section()`, gated on `len(cells) == 1
   and neff_override is not None`: build `T` directly as
   `diag(exp(2j*pi*neff_override[wl_idx]/wl*length))` per wavelength,
   skipping `mw.create_cells()`/`CrossSection.from_cell()`/
   `compute_modes()`/`mw.compute_s_matrix()` entirely on that path.
   **Written, 2026-08-11** — `eme.py` `solve_section()`, the
   `neff_override is not None` branch right after the cache-miss check
   (before, that block was unconditional). `len(cells) == 1` isn't a
   runtime check: `solve_section()` only ever builds `Ls=[length]`, a
   singleton, so the condition collapses to `neff_override is not None` and
   is noted as an invariant, not tested per-call. Diagonal built as
   `T[wl_idx, k, k]` for `k in range(min(n_modes_out, n_modes_in))` —
   `n_modes_out`/`n_modes_in` can differ, and off-diagonal terms are zero
   either way (per the orthogonality argument above), so this covers the
   unequal case too, not just the square one the `diag(...)` shorthand
   above implies. **Gated, not yet the default's only path** — the
   `else` branch (original full solve) is kept intact and still runs
   whenever `neff_override is None`; this is not yet validated against
   step 5, so treat as unshipped until that passes.
3. Decide `warn_on_transfer_duplicate()`'s fate on the fast path (see the
   "not yet checked" item above) — keep it running against the
   directly-built diagonal (weaker signal, still free), skip it explicitly
   with a one-line reason, or repurpose it to compare calibrated β values
   instead of matrix rows. Pick one and record why. **Decided, 2026-08-11:
   keep running as-is.** No code change needed — the call site
   (`warn_on_transfer_duplicate(T, logger, x0, x1)`, `eme.py:574`) sits
   after both the fast and full-solve branches and was untouched by step
   2's edit, so it already runs against whichever `T` came out. Accepted
   as weaker-but-free per the "not yet checked" analysis above; revisit
   only if it turns out to miss something the eigenvector check used to
   catch in practice.
4. Cache-fingerprint check: `build_run_config()`'s `_METHOD_ATTRS` allowlist
   for `("EMESolver", "solve_section")` (`cache/config.py`) must still
   produce the same fingerprint for a fast-path run and a full-solve run of
   the *same* physical config, since they must be numerically identical and
   therefore cache-interchangeable — confirm this holds rather than
   assuming it, since a cache key that accidentally depends on which path
   ran would silently double the cache's footprint for no reason.
   **Checked, 2026-08-11 — holds, by construction.** `config`/
   `family_config` (`eme.py:493-521`) are built entirely before the
   fast/slow branch point introduced in step 2; the fingerprint depends
   only on `_METHOD_ATTRS[("EMESolver","solve_section")]`
   (`material_stack, oxide_width, oxide_height, res, temperature,
   mode_solver`) plus call args (`x0, x1, wavelengths, n_modes_in/out,
   stub_left/right, geometry_kwargs`, and `neff_override`'s own
   fingerprint). None of these are branch-dependent — step 2 added no new
   solver attribute or config field, so there is nothing for the two
   paths to disagree on. (`neff_override`'s fingerprint already
   distinguishes calibrated-vs-uncalibrated results, correctly, since
   those are physically different outputs, not a case this check is
   about.)
5. **Numerical equivalence test, before this ships as a default.** Run both
   paths on a real body (2×2 MMI, `N_MODES=14`) and assert the two `T`
   matrices agree to floating-point precision — same standard this repo
   already holds itself to elsewhere (cf. the HCG Fabry-Pérot `Circuit`
   assembly vs. closed-form cross-check, `max|delta| ~1e-16`). Only after
   this passes does the fast path become the default; keep the full-solve
   path reachable (e.g. a debug/force flag) in case a future finding needs
   real fields here after all.
   **PASSED, 2026-08-11.** One real mode solve (2×2 MMI body, `N_MODES=14`,
   `res=40`, single wavelength) reused for both sides of the comparison:
   the OLD path reconstructed by hand from the pre-edit source
   (`compute_modes()` → `_override_neff()` → `mw.compute_s_matrix()` →
   port-map lookup) versus `EMESolver.solve_section()` called for real
   with the same synthetic `neff_override` (deliberately different from
   the modes' own solved neff — first 5 real values were
   `2.806, 2.766, 2.697, 2.597, 2.464`, override used `2.0, 2.013, 2.026, ...`
   — so the test can't pass by the override silently doing nothing).
   Result: `max|T_full - T_fast| = 0.0` — exact bit-for-bit match, not
   merely within `1e-16`. `T_full`'s off-diagonal also came back exactly
   `0.0`, empirically confirming the single-cross-section orthogonality
   claim (## Why the field profile is not needed) directly, not just the
   override arithmetic. Script:
   `/private/tmp/claude-501/-Users-jhiesener-Documents-Research-duty/a6038aa9-a350-4c1d-a396-b023e0aad6c8/scratchpad/step5_equivalence.py`,
   run in the `jmh` conda env, ~638s wall time, fully local (no cluster).
   **No debug/force flag was added** — `neff_override is not None` is the
   only gate, so the full-solve path stays reachable simply by calling
   `solve_section()` without `neff_override` (a different, physically
   distinct call). Adding a same-input force-full-solve flag remains a
   possible future addition if a real-fields re-check is ever needed, but
   isn't required for the fast path to be correct/default now.

## Status: closed

All 5 steps done, 2026-08-11, entirely local (no cluster) via the `jmh`
conda env. `EMESolver.solve_section()`'s `neff_override is not None`
branch (`eme.py`) is live and is the only path taken on that condition —
no flag gates it further. The header's "BUILD PLAN — mechanism does not
exist yet" no longer applies; leaving the rest of this document as the
record of how that was established.

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
