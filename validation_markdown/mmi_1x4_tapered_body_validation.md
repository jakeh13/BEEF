# 1×4 MMI with a tapered multimode body, 300 nm Si₃N₄ — build plan

Opened **2026-08-20**, at the user's request: a 1×4 MMI splitter on nitride whose
multimode section **tapers in and out** rather than running at constant width, to
shorten the device. **This becomes the manuscript's `sec:1x4_mmi_splitter` test
case**; the existing straight-body `mmi_1x4` stays in the repo as the untapered
baseline to validate against.

> **This doc is a PLAN, not a gap list** — same category as
> [`mmi_1x4_sin300_validation.md`](mmi_1x4_sin300_validation.md) and
> [`mmi_4x4_90_hybrid_validation.md`](mmi_4x4_90_hybrid_validation.md), and it
> depends on the first of those for its material model and cross-section design.
> It converts to a normal [four-checks](README.md#the-four-checks) gap list once
> the scripts run.
>
> **Nothing below has been simulated.** No tapered-body device exists, and the
> library cannot currently propagate one correctly (see *Why this does not work
> today*). Every device-level number here is **ARGUED**. The library-behaviour
> findings are **MEASURED** in the sense that they were read directly off the
> current source, with line references; they are not simulation results.

---

## Summary — where the work actually is

The geometry is the cheap part again, as with the straight 1×4. **The blocker is
that `EMESolver.solve_section()` is single-cell by construction, so no tapered
body can be propagated correctly through a hybrid cascade at all.**

| | status |
|---|---|
| Nitride material model, `SIN_300`, cross-section design | ✅ inherited from `mmi_1x4_sin300_validation.md` steps 1–5 |
| 1×4 port plumbing (N-ary monitors, port lists) | ✅ already supported |
| A *uniform*-body MMI in a hybrid cascade | ✅ already supported — the straight 1×4 and the 2×2 both do it |
| A **tapered-width body** in the `MMI` geometry | ✅ **done 2026-08-20** — `body_taper_length`/`body_end_width` on `MMI`; `width_at`/`_body_polys`/`_check_cut` all updated, untapered limit byte-identical to the pre-2026-08-20 straight body |
| **Multi-cell EME inside a cascade section** | ✅ **done 2026-08-20** — `solve_section()`'s real-solve branch now builds `self.num_cells[0]` cells via the existing `_flat_cell_lengths()`/region-list machinery; per-cell mode count and the `neff_override` diagonal path (step 1) both unchanged |
| `calibrate_body_neff` semantics for a non-uniform body | ✅ **decided 2026-08-21** — option (a) **refuse**, recorded in step 4 and pinned by a test on the real `MMI(body_taper_length=...)`; per-cell (b) / middle-only (c) deliberately deferred. Cost: this device's body carries the ~2.3%/mode MPB gap uncorrected |

Ordered work: **guard the silent-wrong-answer path (step 1) → tapered-body
geometry (step 2) → real multi-cell `solve_section` (step 3) → decide calibration
semantics (step 4) → convergence + untapered-limit check (step 5) → the device and
its four checks (steps 6–7).**

Step 1 is small and lands first on purpose: until it does, anyone pointing a
calibrated hybrid at a tapered body gets a plausible wrong answer with no warning.

> **Progress, 2026-08-21.** Steps **1–4 are done** — the library can now build
> and propagate a tapered body, and refuses to mis-calibrate one. Everything
> remaining is device-and-measurement work, not library work: **step 5**
> (convergence: cell count, per-cell `num_modes`, untapered limit), **step 6**
> (the two device scripts + `Makefile` target), **step 7** (the four checks).
> The *Length payoff* section below is still **ARGUED** and is what step 6
> finally measures.

---

## Why this does not work today

`EMESolver.solve_section()` builds **exactly one meow cell spanning the whole
section**, on both of its branches:

- **`neff_override` branch** (`beef/solvers/eme.py:541-557`) — the fast path added
  by [`calibrate_body_neff_mode_solve_skip_validation.md`](resolved/calibrate_body_neff_mode_solve_skip_validation.md).
  It skips the eigenmode solve entirely and writes a diagonal
  `T[i,i] = exp(2j·π·neff_override/wl·L)`. Its own comment states the assumption:

  > "solve_section always builds a single cell (Ls=[length]), so left/right ports
  > share one cross-section's mode basis -- modes are orthogonal by construction,
  > making the transfer matrix diagonal"

  That reasoning is sound *for a uniform section* and false for a tapered one. The
  gate is `if neff_override is not None` — **not** a uniformity test.

- **`else` branch** (`eme.py:558-581`) — a real meow solve, but still
  `cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])`
  (`eme.py:565`). One cell over the whole span, so a varying cross-section is
  collapsed to a single sampled cross-section rather than propagated through a
  sequence of them.

So there is no multi-cross-section branch in the cascade path to fall back on.
The multi-cell machinery exists only on the whole-device methods —
`_build_cells` (`eme.py:233`), `_modes_per_cell` (`eme.py:254`),
`_flat_cell_lengths` (`eme.py:262`), consumed by `EMESolver.solve()` and
`solve_length_sweep()` (`eme.py:946`), which describes itself as *"the multi-cell
generalization of solve_section_length_sweep()'s single-uniform-cell trick"*.

And nothing guards it. The only uniformity check in the library,
`HybridSolver._validate_uniform_sweep_section` (`hybrid.py:1229-1246`), is called
solely for `solve_length_sweep()`'s swept section. `_calibrated_neff`
(`hybrid.py:1272-1321`) applies `neff_override` to **any** section whose solver
lacks `section_mode_beta()` — i.e. every EME body — with no cross-section test.

*source: `beef/solvers/eme.py:541-581`, `:565`, `:233-270`, `:946-971`;
`beef/solvers/hybrid.py:1229-1246`, `:1272-1321`*

---

## Step 1 — Refuse a non-uniform span instead of silently diagonalising it ✅ **DONE 2026-08-20**

**Smallest change that removes a wrong-answer path.** Reuse the existing checker
rather than writing a second one: `_validate_uniform_sweep_section`
(`hybrid.py:1229-1246`) already samples `width_at`/`center_at` across a section at
`mirror_check_samples` points against `mirror_check_tol`, and already produces a
good error message.

Two candidate placements — the doc should pick one and say why:

- **(a) In `HybridSolver._calibrated_neff`**, before returning a calibration for
  section `k`. Cheap, and keeps the geometry knowledge in the class that already
  owns `boundaries` and the per-section `kwargs`. Does not protect a caller who
  passes `neff_override=` to `EMESolver.solve_section()` directly — which is
  exactly what `mmi_2x2_hybrid.py:298-302` does today.
- **(b) In `EMESolver.solve_section`**, on the `neff_override is not None` branch.
  Catches the direct callers too, but `EMESolver` would need to sample the
  geometry itself (it has `geometry` and `x0`/`x1` in hand, so this is feasible).

Recommendation: **both** — (b) as the real guard, (a) as an earlier, clearer
failure with the section index named. Generalise the checker out of
`HybridSolver` into a shared helper first so there is still one implementation;
the repo has already been bitten by duplicated geometry logic drifting
(`section_fdtd_layout()` is deliberately shared with
`beef/plot/geometry.py::_draw_section_fdtd_markers()` for this reason).

**Done when:** a tapered EME body with `calibrate_body_neff=True` raises with the
offending x named, and the straight 1×4 / 2×2 / coupler all still pass
`pytest tests/ -m "not slow"` bit-identically (their bodies are uniform, so the
guard must be inert for them — a cache-hit run is the cheap evidence).

### What landed

Plan followed as recommended — shared helper first, then **both** placements.
Every source claim in *Why this does not work today* was re-verified against the
freshly-pulled tree before building on it, and all held exactly (both
`solve_section` branches really do build `Ls=[length]`; the fast-path gate really
is `if neff_override is not None`, not a uniformity test).

- **`assert_uniform_cross_section()`** — new shared helper in
  `beef/geometry/base.py`, the single implementation. Takes `what`/`why` clauses
  so each caller states its own requirement while the sampling logic stays in one
  place. Kept the module **stdlib-only** (plain arithmetic instead of
  `np.linspace`) rather than adding a numpy import to a geometry predicate.
- **`HybridSolver._validate_uniform_sweep_section`** — refactored to a thin
  wrapper. Same test, same message, no behaviour change.
- **(a) `HybridSolver._calibrated_neff`** — guards before calibrating, naming the
  section (`"Section 2/3 (calibrate_body_neff=True)"`).
- **(b) `EMESolver.solve_section`**, `neff_override` branch — the real guard,
  catching direct callers `HybridSolver` cannot see. Checks
  `[x0 - stub_left, x1 + stub_right]`, not `[x0, x1]`: `length` includes both
  stubs, so both must be uniform too.

**Evidence.** Full fast suite **674 passed, 0 failed** (job 12066434). Guard
inertness was additionally confirmed by sampling the **real** production
geometries directly — 1×4 (both the 39.90/3.1 and old 206/2.4 geometries), 2×2,
1×2, and the directional coupler all pass over their actual EME body spans —
and by confirming it *fires* on a span crossing the 1×4's input taper band and
on the coupler's own S-bend.

New tests: `TestCalibrationRefusesANonUniformBody` — parametrized over a
width-varying body **and** a constant-width/walking-centre (S-bend-shaped) body,
plus an inertness case. The second parametrization earns its place: a
width-only check would wrongly pass it, which is exactly why the helper tests
`center_at` as well as `width_at`.

> ⚠ **One deviation, and the reasoning, since it is the kind of choice that
> looks settled and is not.** The first run came back **5 failed** — all in
> `tests/test_beta_calibration.py`, all passing `geometry=None` and so hitting
> `AttributeError` inside the new guard. (Same failure shape
> `resolved/section_mode_beta_multi_monitor_validation.md` step 7 already
> recorded for that same file.) There were two candidate fixes: soften the guard
> to skip when a geometry lacks `width_at`/`center_at` — there is precedent, as
> `check_mode_correspondence` is explicitly documented as best-effort — or fix
> the tests. **Chose the tests.** `DeviceGeometry.width_at()`'s own
> `NotImplementedError` says it is *"required for use with HybridSolver"*, so
> only stubs and `None` can lack it, never a real geometry; and a silent-skip
> path inside a guard whose entire purpose is to prevent a silent wrong answer
> would be self-defeating. Added a `UniformStubGeometry` and used it **only** in
> the two call sites that reach this path, leaving the other four `geometry=None`
> uses untouched so no passing test was disturbed.

---

## Step 2 — Tapered-body geometry ✅ **DONE 2026-08-20**

**Design decision to make first.** Three routes:

1. **Extend `MMI` with body-taper bands.** The body is currently one constant-width
   rectangle (`mmi.py:_section_polygons`, `_rect(max(lo, in1), min(hi, out0), W/2)`).
   Adding e.g. `body_taper_length` and a start/end body width means extending
   `_taper_bounds()` to report the two new bands, drawing them in
   `_section_polygons`, teaching `width_at` to interpolate inside them, and
   extending `_check_cut` to forbid cuts strictly inside them.
2. **`CompoundGeometry` chain** (`geometry/compound.py:37`) — concatenate
   segments along x with each query translated into the owning segment's local
   frame. Already backs a 7-section chain. But it composes whole geometries, and
   the 1×4's four output branches live inside one `MMI`, so the branching
   structure does not decompose cleanly.
3. **New `TaperedMMI` class.** Cleanest conceptually, most duplicated code.

Recommendation: **route 1.** All the machinery that matters — N-ary ports
(`mmi.py:29-34`, `_port_ys`), `_check_cut`, `_section_polygons`,
`section_prisms`/`section_gds`, `width_at`/`center_at`/`source_aperture_at` — is
already generic and would be inherited. Route 3 re-implements all of it.

**Note the interaction with `_check_cut`.** Once the body has taper bands, the
hybrid's cuts must land at band endpoints or in the uniform middle, exactly as
they do today for the input/output taper bands. That is *not* a workaround for
step 3 — the FDTD/EME cut still lands on a constant cross-section, and it is the
**EME section itself** that spans a varying body. Those are different
requirements and both must hold.

**Done when:** `width_at` traces the intended profile across the whole device,
`_check_cut` rejects a mid-body-taper cut, `plot_hybrid_geometry()` draws the
device correctly, and the untapered limit (body taper length → 0) reproduces the
straight `MMI` polygons exactly. Geometry-only, so it belongs in
`tests/test_geometry.py` (67 tests, imports neither meep nor meow).

### What landed

Route 1 as recommended — `MMI` extended in place, nothing new duplicated.

- **`MMI.__init__`** — two new keyword args, `body_taper_length: float = 0.0` and
  `body_end_width: Optional[float] = None` (defaults to `mm_width`, so an
  unspecified end width is a no-op, not an error). Validated: negative length
  rejected; `2 * body_taper_length > mm_length` rejected (the two bands would
  overlap, leaving no uniform middle — message names both numbers and the max
  allowed value); non-positive end width rejected.
- **`is_body_tapered`** — `True` only when *both* `body_taper_length > 0` and
  `body_end_width != mm_width`; a zero-length taper with a end-width set is still
  untapered by construction (there is no band for it to apply over), so the
  untapered fast path stays reachable from either default.
- **`_body_taper_bounds()`/`body_width_at()`/`_body_polys()`** — the single shared
  implementation, called from both `polygons()` (whole device) and
  `_section_polygons()` (a cascade section), same reasoning as
  `assert_uniform_cross_section()` in step 1: one place that can disagree with
  itself is one place, not two. `_body_polys` returns the untapered case as
  **one** rectangle, byte-identical to the pre-2026-08-20 code path, and the
  tapered case as three trapezoids (taper-in, uniform middle, taper-out) — all
  three all-or-nothing per `_check_cut`, matching the port-taper convention.
- **`width_at`** — the body branch (`-half <= x <= half`) now calls
  `body_width_at(x)` instead of returning `mm_width` unconditionally.
- **`center_at`** — confirmed to need **no change**: the body taper is
  width-only and symmetric about y=0, so its existing unconditional `return 0.0`
  inside the body span already covers it. Checked explicitly rather than assumed,
  since step 1's own `TestCalibrationRefusesANonUniformBody` precedent is exactly
  that a width-only check can miss a centre-only problem — here it's the
  reverse (a width-only *change* with no centre effect), confirmed by direct
  sampling, not by symmetry argument alone.
- **`_check_cut`** — new block rejecting a cut strictly inside either body taper
  band, same message shape as the existing port-taper-band rejection, naming
  both band bounds.
- Two now-dead `W = self.mm_width` locals removed (`polygons()`, then found and
  removed a second one in `_section_polygons()` on the same pass) once both
  methods' body-drawing lines were replaced by the shared `_body_polys()` call —
  caught by grepping `\bW\b` after each edit rather than assuming the first
  removal was the only one.

**Evidence.** New `TestMMIBodyTaper` class in `tests/test_geometry.py`, 17 tests:
untapered-limit regression (`polygons()` byte-identical to a plain `MMI()` with
no body-taper kwargs), validation errors, tapered profile correctness (band
endpoints, uniform-middle interior, linear-ramp midpoint), `center_at`
invariance, `_body_polys` trapezoid count, and `_check_cut` firing inside a body
taper band while still passing at its boundary. **17 passed**, and the full fast
suite regressed clean at **725 passed, 0 failed** (job 12066598) — same
zero-failure bar as step 1, now against a tree that also carries some tests
added by other concurrent work on this shared repo (674 → 725 baseline shift is
not this step's tests alone; the 17 new ones are a strict subset).
`plot_hybrid_geometry()` was confirmed by direct render (`to_gds_eme()` →
`polygons()` → `_body_polys()` is the same call chain the plot function already
used, so no plot-code change was needed) — a 1×4 with `body_taper_length=6.0`,
`body_end_width=8.0` against `mm_width=12.0` draws the expected shape: narrow at
the port-taper interface, widening to the full body width over 6 µm, holding
uniform, then symmetrically narrowing before the output fan-out.

---

## Step 3 — Real multi-cell EME inside a cascade section ✅ **DONE 2026-08-20**

The substantive library change. `solve_section` must be able to subdivide
`[x0, x1]` into cells, solve modes per cell, and let
`mw.compute_s_matrix(modes, cells)` handle the interface overlaps between
successive cross-sections — which is precisely the mode-field computation the
manuscript's Solvers subsection says a tapered body requires.

Points to settle in the doc:

- **Where the cell count comes from.** `EMESolver.__init__` already carries
  `num_modes` / `num_cells` / `lengths` as per-region lists for the whole-device
  path, but `solve_section` currently ignores them and uses `Ls=[length]`. Decide
  whether a cascade section reuses those lists, takes a new
  `section_cells`-style argument, or derives a count from a width-change
  tolerance. Deriving it is attractive (the caller should not have to hand-tune
  per section) but adds a convergence knob that step 5 then has to pin down.
- **Reuse, do not duplicate.** `_flat_cell_lengths` (`eme.py:262`) and
  `_modes_per_cell` (`eme.py:254`) already do the region → flat-cell expansion;
  the section path should call them rather than growing a parallel copy.
- **Cache fingerprint.** Any new argument that changes the answer must reach
  `beef/cache/config.py`'s `_METHOD_ATTRS` allowlist for
  `("EMESolver", "solve_section")`. This has bitten the repo twice — the
  `minimum_run_time` false-positive hit (`cache/config.py:98-127`) and the
  `dispersive` byte-identical key (`materials/stack.py:619-626`). A new cell-count
  argument silently omitted would serve a single-cell result for a multi-cell
  request.
- **`solve_section_length_sweep` interaction.** Its analytic phase trick is
  uniform-section-only by construction; a tapered body simply cannot use it. So
  the tapered 1×4 has no length sweep, which is fine — but it means the
  manuscript's `solve_length_sweep()` claim must not be attached to this device.

**Done when:** a tapered body's transfer matrix is built from a real multi-cell
meow solve; the untapered limit reproduces the single-cell result to a stated
tolerance (the natural regression test, cheap because it needs only a small
section); and the duplicate-mode guards still run — `compute_modes()`'s
`guard_duplicate_modes` per cell, plus `warn_on_transfer_duplicate` on the
assembled `T` (`eme.py:583-592`).

### Recommendation (resolving the three points above)

> **Superseded in part — read the ⚠ correction at the end of *What landed*
> before acting on this paragraph.** Its conclusion for `num_cells` holds; its
> conclusion for `num_modes` was wrong, and the grep it rests on is the reason
> why.

**Cell count and per-cell mode count: reuse `num_modes`/`num_cells` exactly as
constructed today, not a new argument.** Every existing `*_hybrid.py` script
already builds its body `EMESolver` with single-element lists —
`num_modes=[N_MODES], num_cells=[1], lengths=[1.0]` — and the `num_modes=[N_MODES]`
value already, if accidentally, equals both `n_modes_in` and `n_modes_out` at
every call site (checked by grep across every `*_hybrid.py`/`*_hybrid_*.py`
device script and every test that calls `.solve_section()` — none passes a
multi-region list there; the multi-region examples, e.g.
`mmi_4x4_90hybrid_eme.py`'s `num_modes=[1, N_MODES_TAPER, N_MODES_BODY,
N_MODES_TAPER, 1]`, only ever feed the whole-device `solve()`/`solve_length_sweep()`
methods, never `solve_section()`). So `solve_section()` can validate
`len(self.num_modes) == len(self.num_cells) == 1` (a HybridSolver cascade
already carries one dedicated `EMESolver` per section — a multi-region instance
reaching this path is a caller error, not a case to silently half-honor) and
then read `self.num_cells[0]` as this section's cell count and
`self.num_modes[0]` as every cell's mode count, via the exact
`_flat_cell_lengths()`/`_modes_per_cell()` helpers the whole-device path already
uses — **zero new arguments, zero new duplicated region→cell logic**, and
`num_cells[0]=1` (today's value everywhere) reproduces today's single cell
exactly, `Ls=[length]` included (stubs folded into the one length, unchanged).
A caller that wants the body to carry more internal resolution than its
FDTD-facing interfaces expose (e.g. `num_modes=[20]` internally against
`n_modes_in=n_modes_out=12`) gets it for free via the same over-request-then-slice
`port_map` read-out `solve_section()` already does — validated by requiring
`self.num_modes[0] >= max(n_modes_in, n_modes_out)`.

**Port-map mechanics need no new code.** Read meow's own netlist assembly
(`meow/eme/cascade.py:_get_netlist`, in the `meow_MEEP` checkout this repo's
editable install points at): external ports are built *only* from the first
cell's `left@i` faces and the last cell's `right@i` faces — every interior
cell-to-cell interface is wired internally and never exposed. So
`solve_section()`'s existing `port_map[f"right@{i}"]`/`port_map[f"left@{j}"]`
read-out loop is already correct for N cells, not just one; the only things
that actually change are the `cells`/`modes` lists fed into
`mw.compute_s_matrix()`.

**Cache fingerprint.** `num_modes`/`num_cells` move from unread (today) to read
by `solve_section()`'s real-solve branch, so both must join
`_METHOD_ATTRS[("EMESolver", "solve_section")]`. Every existing cache entry
predates these keys entirely and so will not match post-change — a full,
one-time miss-and-resimulate for that method, same acceptable cost this file
already accepts for a tracked-field addition (`cache/config.py`'s own
`minimum_run_time` precedent), not a silent false-positive risk.
`solve_section_length_sweep()` is untouched — it stays single-cell/diagonal by
construction (this step doesn't reach it) — so its allowlist entry does not
change.

### What landed

Implemented exactly as recommended above — no design deviation this time.

- **`EMESolver.solve_section()`**, real-solve (`neff_override is None`) branch:
  validates `len(self.num_cells) == 1` (raises otherwise), then builds
  `cell_lengths = self._flat_cell_lengths(region_lengths=[length])` and
  `modes_per_cell = [n_modes] * len(cell_lengths)` — reusing the existing
  helper exactly, not a parallel copy — and calls `mw.compute_s_matrix()` over
  the full per-cell `cells`/`modes` lists. The `port_map` read-out loop
  (`right@i`/`left@j`) is untouched: confirmed against meow's own
  `_get_netlist()` (`meow_MEEP/src/meow/eme/cascade.py`, the editable-install
  checkout this repo's `meow` points at) that external ports come only from
  the first/last cell's faces regardless of cell count, so the existing
  lookup was already correct for N cells. `num_cells[0]=1` (every existing
  device script's value, unchanged) reproduces the old `Ls=[length]` call
  exactly. The `neff_override` branch (diagonal fast path, guarded by step
  1) is untouched.
- **`beef/cache/config.py`** — `num_cells` added to
  `_METHOD_ATTRS[("EMESolver", "solve_section")]`. `num_modes` deliberately
  **not** added (see the correction below).
  `solve_section_length_sweep()`'s entry is unchanged.

> ⚠ **Correction, made 2026-08-21 before step 5 — the recommendation above was
> wrong on one point, and only checking the real callers caught it.** The plan
> said to read *both* `num_modes[0]` (per-cell mode count) and `num_cells[0]`,
> justified by "every call site already has `num_modes[0]` equal to both
> interface counts" — which was verified by grepping the `EMESolver(...)`
> constructor lines, and **that grep was the wrong check**. It confirmed the
> two numbers match where a script defines them *together*, but not that they
> match at every *call*. `modal_multiplexer_hybrid.py` fills **three** EME
> slots with **one** `body_solver(num_modes=[15])`, and one of those slots is a
> connector section whose interfaces carry `CONN_MODES=6`. Reading
> `num_modes[0]` there silently solves **15 modes where 6 were solved before**:
> a behaviour change (the top-6 slice of a 15-mode solve is not guaranteed
> identical to a direct 6-mode solve — the same caveat `solve_section()`'s own
> superset-cache comment already records) *and* a cost increase, on a device
> nobody was touching. The script even says so in its own comment: *"solve_section()
> takes n_modes from interface_modes above, not from body_solver's own
> constructor num_modes, so no new EMESolver is needed for the connector
> section."*
>
> **Fixed by reading only `num_cells`.** The per-cell mode count stays
> `max(n_modes_in, n_modes_out)`, exactly as before this step. Step 3's actual
> goal — a real multi-cell body — is unaffected, since that is what `num_cells`
> delivers; `num_modes` was never needed for it. A caller wanting a richer
> internal body basis than its interfaces expose should raise
> `n_modes_in`/`n_modes_out`, which it controls directly. This also removes the
> `num_modes[0] >= max(...)` validation (nothing to validate now) and keeps
> `num_modes` out of the cache key, so existing entries are not stranded on an
> attribute the method does not read. Pinned by
> `test_solve_section_does_not_read_num_modes` (two solvers differing only in
> `num_modes` must return a bit-identical `T`, `atol=0`) and by
> `test_cache_config.py::...::test_eme_solve_section_still_ignores_num_modes`.
>
> The lesson worth keeping: *"every call site already does X"* is a claim about
> **calls**, and a grep over **constructors** cannot establish it — one shared
> instance reused across several sections breaks the correspondence silently.

**Evidence.** New `tests/test_eme_multicell_section.py`, real meow (`--runslow`):
- `test_multicell_matches_single_cell_on_a_uniform_span` — **the untapered-limit
  regression this step's own "done when" requires.** `num_cells=[1]` vs
  `num_cells=[3]` over the *same* uniform span agree to `atol=5e-3` on `|T|`:
  measured cross-solve noise between independent eigensolves of an identical
  cross-section, not a body-taper effect (this span has none). **Passed.**
- `test_multicell_body_taper_changes_the_transfer_matrix` — the same
  multi-cell machinery, tapered body vs. straight body (step 2's
  `body_taper_length=1.0, body_end_width=1.5` vs. the default). Measured:
  `|T|` diagonal is exactly 1.0 for the straight body (a uniform body is
  lossless per mode, as expected) vs. **0.9802/0.9884 for the tapered one** —
  a real, mode-dependent ~1–2% amplitude effect an order of magnitude past the
  ~5e-3 noise floor the previous test establishes, i.e. the taper genuinely
  reaches the eigensolve rather than every cell silently sampling the same
  cross-section. **First deviation, caught and fixed by the numbers, not
  argued around them:** the test's original `atol=0.02` sat right at the
  0.0198 measured difference (a guess made before running, not derived from
  it) and so returned a false "indistinguishable" failure on the first run
  (job 12069743); tightened to `atol=5e-3` — between the noise floor and the
  measured effect — against the real numbers, then re-run clean (job
  12069809).
- Two fast (no real eigensolve reached) validation tests for the new guards.
- **Regression caught in the existing suite, not just added tests:**
  `test_cache_config.py::TestRealSolverAllowlists::test_eme_solve_section_ignores_num_modes_and_num_cells`
  encoded the pre-step-3 invariant this step deliberately reverses. Renamed to
  `..._is_sensitive_to_num_modes_and_num_cells` and its assertion flipped
  (`!=` in place of `==`), rather than left failing or the allowlist change
  reverted — the test was asserting the *old*, now-intentionally-wrong
  behaviour, the same category of update step 1 made to
  `test_beta_calibration.py` for a guard's own new failure mode.
- Full fast suite: **727 passed, 0 failed** (job 12069809) — up from step 2's
  725 baseline by exactly the 2 new fast tests in this step's own file, net of
  the fixed pre-existing one.

---

## Step 4 — What `calibrate_body_neff` means for a multi-cell body ✅ **DECIDED 2026-08-21**

**The crux, and the reason step 1 exists.** For a uniform single-cell body, one
per-mode β fully determines `T` — that is the whole basis of the fast path. For a
tapered body it does not: the transfer matrix depends on the mode *fields* at
every internal interface, not just on how fast each mode accumulates phase.

Options:

- **(a) Refuse.** `calibrate_body_neff=True` on a non-uniform body raises. Honest,
  trivially correct, and step 1 already implements it. Cost: the tapered 1×4 gets
  no calibration, so its body propagates on MPB's n_eff and inherits the ~2.3%/mode
  MPB-vs-time-stepper gap documented in
  [`resolved/directional_coupler_validation.md`](resolved/directional_coupler_validation.md)
  item 6.
- **(b) Per-cell calibration.** Measure β per cell cross-section and override each
  cell's modes. Physically the right generalisation, but the cost scales with cell
  count — and β calibration is already the expensive step (89% MPB modal
  decomposition, `resolved/section_mode_beta_multi_monitor_validation.md` step 9).
  An N-cell body means N calibrations unless cells share cross-sections.
- **(c) Calibrate the uniform middle only**, leave the taper cells on MPB. A
  defensible middle ground if the middle carries most of the phase.

Recommendation: ship **(a)** with step 1, and treat (b)/(c) as a separate
follow-on with its own measurement. **Do not** let the tapered 1×4 block on
solving calibration semantics — the device's value in the manuscript is that it
exercises multi-cell EME with real mode fields, which is orthogonal.

Whichever is chosen, record it explicitly: this is the kind of question that
looks settled and is not.

### Decision: **(a) refuse**, recorded explicitly

**A tapered body gets no β calibration.** It propagates on the mode solver's own
n_eff and inherits the ~2.3%/mode MPB-vs-time-stepper gap
(`resolved/directional_coupler_validation.md` item 6). Options (b) per-cell and
(c) middle-only are **deferred**, each needing its own measurement — see *What
this costs* below for why deferring is defensible rather than merely convenient.

No new library code: step 1's guard already implements exactly this, in both
placements. What step 4 adds is the **decision, recorded**, plus the test that
pins it to the real geometry.

**What landed.** `TestCalibrationRefusesARealTaperedMMIBody` in
`tests/test_beta_calibration.py`, two tests:

- `test_it_refuses_a_real_tapered_mmi_body` — a real
  `MMI(body_taper_length=4.0, body_end_width=8.0)`, body section spanning the
  whole body with both taper bands included (exactly what the tapered 1×4's own
  cascade will ask for), raises through `HybridSolver._calibrated_neff()`.
- `test_it_still_calibrates_the_straight_mmi_body` — the untapered default
  still calibrates and returns the right shape. This is the existing straight
  1×4's own production path, and step 2's whole premise is that it is
  undisturbed.

**Why this test exists on top of step 1's.** Step 1's own guard tests
(`TestCalibrationRefusesANonUniformBody`) predate step 2, so they could only use
hand-written stub geometries (`TaperedStubGeometry`, `WalkingCentreStubGeometry`).
Those prove the *guard* works; they cannot prove the guard can *see a real
`MMI`'s taper*, because no such geometry existed when they were written. A
future change to `MMI.width_at`/`body_width_at` that hid the taper from the
guard would leave every stub test green while silently restoring the
wrong-answer path. This closes that loop on the actual class the device will be
built from. **5 passed** (3 stub + 2 real, job 12078941); full fast suite
**746 passed, 0 failed**. Note the suite total moved 727 → 746 while this step
added only 2 tests — the other +17 are another concurrent session's work on this
shared tree, so the absolute count is not a per-step signal here; the
zero-failure result is.

**What this costs, stated plainly.** The tapered 1×4 is the one device in this
repo whose body cannot be calibrated, so it is also the one whose EME body
carries the MPB-vs-time-stepper gap uncorrected. That is a real accuracy
concession, not a free choice — worth re-reading before comparing this device's
absolute numbers against the straight 1×4's calibrated ones. It is defensible
because the device's value in the manuscript is that it exercises multi-cell EME
with real mode fields (step 3), which is orthogonal to calibration; and because
(b) would multiply the already-dominant calibration cost by the cell count (89%
MPB modal decomposition, `resolved/section_mode_beta_multi_monitor_validation.md`
step 9) for an accuracy gain nobody has measured yet on this geometry. Revisit
if step 5's convergence or step 7's check 1 shows the uncalibrated body is the
binding error term.

---

## ⛔ Blocker found 2026-08-21, before step 5 — the symmetric taper has almost no room

**MEASURED (pure geometry, no solver).** The 1×4's own output fan-out consumes
essentially the entire body face, so a **symmetric** body taper — which is what
`MMI.body_end_width` implements (step 2) — can narrow this device by at most
**4.2%**:

| face | ports there | required face width | headroom vs `mm_width=12.0` |
|---|---|---|---|
| **output** (`x=+L/2`) | 4 tapers at y = ±1.55, ±4.65, each `taper_width=2.2` wide | 2×(4.65 + 1.1) = **11.50 µm** | **0.50 µm total (4.2%)** |
| **input** (`x=−L/2`) | 1 taper on axis, `taper_width=2.2` | **2.20 µm** | **9.80 µm (81.7%)** |

The outermost output port sits at |y| = 4.65 and is 2.2 µm wide at the body
face, so it reaches |y| = 5.75 against the body's own 6.0 half-width. **A
symmetric taper is capped by the output face**, and 4.2% of width will not
produce a length saving worth reporting — `L_π ∝ W²`, so 4.2% off the *faces
only* (the uniform middle is unchanged) is a fraction of a percent of device
length.

The input face, by contrast, is 81.7% free.

**This was not anticipated by step 2**, which added a single `body_end_width`
applying to both faces without checking whether the port layout left room for
it. The geometry code is correct; the *design* it enables is the constrained
one. Verified against the real production constants
(`MM_WIDTH=12.0`, `WAVEGUIDE_OFFSET=3.1`, `TAPER_WIDTH=2.2`) — and note
`WAVEGUIDE_OFFSET=3.1` was itself only just optimised in
`mmi_1x4_sin300_validation.md` step 9c, so moving the ports to free up room
re-opens a settled question rather than being free.

Options, none yet chosen — **this changes the device and the manuscript's
footprint claim, so it is not mine to settle**:

- **(a) Asymmetric taper** — narrow at the input face, full width at the
  output face. Uses the 81.7% headroom that actually exists. Costs a step 2
  revision: `body_end_width` splits into a start/end pair (the geometry
  helpers `_body_taper_bounds`/`body_width_at`/`_body_polys` already handle
  two independent bands, so this is a parameter change, not a rewrite).
- **(b) Move the ports inward** to make room for a symmetric taper — re-opens
  step 9c's optimised `WAVEGUIDE_OFFSET`.
- **(c) Accept the 4.2% cap** and report the payoff as negligible. The doc's
  own *Length payoff* section pre-authorises exactly this: *"If the saving is
  not material, the device is still worth building — it is what makes the
  manuscript's tapered-body claim real — but say so plainly rather than
  overselling the footprint."*

**Step 5's convergence sweeps are deliberately not started until this is
settled** — cell-count and per-cell-mode convergence are properties *of a
specific taper profile*, so running them against a profile that may be
replaced would be measuring the wrong device.

### Resolved 2026-08-21 — option (a), asymmetric taper

**Chosen by the user: narrow the input face, leave the output face at
`mm_width`.** Step 2's geometry was revised accordingly, same day.

- `MMI.body_end_width` split into **`body_start_width`** (input face) and
  **`body_end_width`** (output face), each defaulting to `None` → `mm_width`.
  The two bands were already independent inside
  `_body_taper_bounds()`/`body_width_at()`/`_body_polys()`, so this was a
  parameter change, as predicted — no restructuring.
- `is_body_tapered` now fires when **either** face differs from the middle
  (previously it tested one shared value), so the asymmetric profile is
  correctly recognised as tapered while both-faces-equal stays straight.
- **`_check_cut` was the one non-obvious consequence.** With an asymmetric
  profile the untapered side's band is geometrically a plain *rectangle*, so
  the original guard would have refused cuts on a genuinely uniform
  cross-section — forbidding valid section boundaries, the opposite of its
  purpose. It now skips any band whose face width equals `mm_width`, and only
  refuses bands where the width actually varies. Pinned by
  `test_a_face_at_mm_width_does_not_block_cuts`, which asserts both halves:
  the flat band allows a cut, the varying band still refuses.

Verified on the real production constants — `body_taper_length=8.0,
body_start_width=6.0` against `mm_width=12.0` gives 6.0 µm at the input face,
ramping to 12.0 over 8 µm, then flat through the output face (which clears the
11.5 µm the fan-out needs). Rendered and visually confirmed. New tests:
`test_asymmetric_faces_are_independent`,
`test_asymmetric_ramps_only_on_the_varying_side`,
`test_a_face_at_mm_width_does_not_block_cuts`,
`test_the_real_1x4_asymmetric_profile` (which re-derives the 11.5 µm fan-out
requirement from `output_ys` rather than hardcoding it, so a future port-layout
change fails the test instead of silently invalidating the profile), plus
`test_both_faces_at_mm_width_is_not_body_tapered` and
`test_default_body_start_width_is_mm_width`. The nine pre-existing symmetric
cases were made *explicitly* symmetric (`body_start_width=` and
`body_end_width=` both passed) rather than deleted — the symmetric profile is
still supported and still worth covering.

### ❌ Face-taper (body_start_width/body_end_width) profile: DEAD END, closed 2026-08-27

Full grid (jobs 12383125, 12387482), 4 (Ws, Lt) combinations x 5 lengths each,
whole-device EME at design settings (res=30, 16x4):

| Ws | Lt | L_best | totalT | balance | vs straight (L=39.90, T=0.9703, bal=2.38e-03) |
|---|---|---|---|---|---|
| 6.0 | 8.0 | 42.00 | 0.9136 | 4.44e-03 | **+2.10 µm LONGER**, worse T, worse balance |
| 6.0 | 4.0 | 42.00 | 0.9140 | 5.98e-03 | same |
| 9.0 | 8.0 | 42.00 | 0.9136 | 5.17e-03 | same |
| 9.0 | 4.0 | 42.00 | 0.9130 | 5.29e-03 | same |

**Every combination converges to the same L=42.00, regardless of how
aggressive the taper is** (Ws 6→9, a 50% change in narrowing; Lt 4→8, a 2x
change in taper length) — totalT and balance barely move (0.9130-0.9140,
4.4-6.0e-03). This means the face-taper amount is not the variable governing
self-imaging length here; the taper is too small a perturbation on the
dominant uniform-mm_width middle to matter, and mm_length alone sets the
optimum. **The device gets LONGER, not shorter, and worse on both metrics, at
every profile tried.** Earlier (jobs 12383125 at the ARGUED-predicted L=31.90)
found outright imaging failure (totalT=0.554); this grid shows the true
optimum sits past the straight device's own length, not before it.

**Verdict: face-taper does not deliver the length payoff, full stop, for
Ws in [6,9] / Lt in [4,8].** Superseded by the WAIST profile (below), which
sidesteps the port-fan-out constraint that motivated the face-taper in the
first place rather than fighting it.

### ⚠ Waist profile, first real result: Ww=6.0 is too aggressive (job 12402940)

Corrected smoke (post N_LEN=1 bug fix, see the script's own inline comment),
real 3-length window centred on L_straight=39.90:

| L | totalT | balance |
|---|---|---|
| 35.91 | 0.2798 | 1.26e-01 |
| 39.90 | 0.2338 | 9.79e-02 |
| 43.89 | 0.2720 | 1.13e-01 |
| control (39.90, straight) | **0.9703** | **2.38e-03** |

**Catastrophic across the whole window, not an edge artefact this time** —
totalT 0.23-0.28 at all three lengths (70-77% power loss). Same qualitative
failure the face-taper profile showed at its aggressive settings: halving
the width (Ww=6 vs mm_width=12, a 50% narrowing) breaks self-imaging rather
than just costing some radiation loss. Mode-reuse mechanism itself worked
correctly throughout (104 solves/40 reused per point, confirming the width-
keyed cache functions) -- this is a real profile result, not a script bug.

**Natural next step, NOT yet run:** gentler Ww (9-10 range, i.e. 15-25%
narrowing rather than 50%) — the same lesson the face-taper sweep already
taught (Ws=6/9 both failed; never tried Ws=10-11). Untested whether a mild
waist gives ANY net length saving once it stops breaking imaging, or whether
this device's self-imaging condition is simply too sensitive to width
changes anywhere along its length to benefit from tapering at all.

### ⚠ Waist profile, Ww=9.0 smoke: still fails, less catastrophically (job 12435720)

Single-point smoke at design centre length only (`N_LEN=1`, reduced cells
`C_TRANS=4/C_WAIST=2` for speed), Ww=9.0 (25% narrowing, vs. Ww=6.0's 50%):

| L | totalT | balance |
|---|---|---|
| 39.90 | 0.6609 | 2.05e-02 |
| control (39.90, straight) | **0.9703** | **2.38e-03** |

**Better than Ww=6.0 (0.66 vs. 0.23-0.28) but still a real failure** — 34%
power loss and balance ~8.6x worse than control, not just a minor cost. Not
the outright imaging collapse Ww=6.0 showed, so narrowing amount clearly
matters (unlike the face-taper's Ws=6/9 insensitivity) — but 25% narrowing at
this single length is still well outside "holds against the straight
control." Flagged by the script itself as landing ON WINDOW EDGE (a
single-point window, so not diagnostic on its own — the real question is
whether a length *other* than 39.90 recovers totalT for this Ww, since a
narrower interior may shift the self-imaging length rather than just costing
loss at the straight device's own length).

### ❌ Waist profile: gentle Ww also fails, CLOSED 2026-08-28 (job 12436352)

Full grid, Ww in {9.0, 10.0} x 5 lengths around centre (33.91-45.88), design
settings (res=30, 16x4):

| Ww | L_best | totalT | balance | vs control (L=39.90, T=0.9703, bal=2.38e-03) |
|---|---|---|---|---|
| 9.0 | 36.91 | 0.6186 | 1.26e-02 | 2.99 µm SHORTER, but 36% power loss, 5.3x worse balance |
| 10.0 | 45.88 | 0.6858 | 6.71e-02 | 5.98 µm LONGER, 29% power loss, 28x worse balance, **on window edge** |

**No (Ww, L) point in this grid holds totalT/balance against the straight
control.** Ww=9.0's best is genuinely shorter but still loses over a third
of the power; Ww=10.0's best sits at the edge of the tested length window
(true optimum may lie further out) but is already longer than the control
there and no better on either metric, so chasing the edge further cannot
produce a win against this device's own reason for existing (a real length
saving at matched performance).

**Verdict: the waist profile fails at every Ww tried (6.0, 9.0, 10.0), same
conclusion as the face-taper.** Two independent taper-profile shapes
(face-taper: narrow one/both faces; waist: narrow the interior only) have
now both failed to produce ANY design that is shorter than the straight
1x4 while holding its totalT/balance. This is evidence, not just an
unlucky parameter search: this device's self-imaging condition appears too
sensitive to width perturbation ANYWHERE along its length to benefit from
tapering at all, at 1x4's own aspect ratio (mm_width=12.0, mm_length=39.90).
**Tapered-body premise for this device is closed pending a new idea** —
not a parameter to keep sweeping.

**Still to pick, in step 6:** the actual `body_taper_length` /
`body_start_width` values. The 6.0/8.0 pair above is a **placeholder used for
verification only**, not an optimised design — no sweep has been run.

### ⚠ Waist profile REOPENED 2026-08-28: field-plot diagnostic finds cause for doubt on the closure

New script `plot_mmi1x4_waist_profile_eme_fields.py` (copies the sweep
script's own body-segment cell-building rather than importing it, since
`mmi_1x4_waist_profile_sweep.py` has no `__main__` guard and would run its
whole grid as an import side effect) renders straight-vs-waisted EME fields
at a chosen design point, to answer the closure note's own open diagnostic
question: is the totalT loss localized at the two width-transition
interfaces, or spread through the whole waist region?

**Ww=9.0, L=36.91 (job 12445061, real design settings, res=30, matching the
grid's own C_TRANS=6):** confirms the smoke-resolution result exactly. The
field narrows correctly into the waist, but a **spurious bright central
lobe reappears around z=35-40**, right at the up-transition, before finally
spreading into 4 output lobes over a much shorter remaining distance than
the straight control needs. The refocus sits exactly where the interior's
accumulated phase first re-expresses in the wide-body mode basis, so this
picture does not cleanly separate "scattering happens at that step" from
"phase error built up gradually through the waist and only becomes visible
once re-projected there" — a per-mode power trace through the waist would
be needed to fully disambiguate, not run.

**Ww=10.0, L=45.88, C_TRANS DOUBLED to 12 (job 12446756, at the user's
request):** qualitatively DIFFERENT field character, not just a smoother
version of the Ww=9.0 picture. Instead of one spurious refocus, the field
shows **genuine multi-lobe interference developing progressively** through
z=28-45 — several bright/dark band revivals near the eventual output
positions, closer to real Talbot self-imaging structure than the Ww=9.0
case, just with more revival cycles than the straight control's single
clean image. This is suggestive that:
1. The original grid's `C_TRANS=6` may have been under-resolving the
   transition itself (a numerical confound on TOP of the physical
   sensitivity finding, not instead of it) — untested until now.
2. L=45.88 (the grid's own window-edge result for Ww=10.0) may have simply
   caught this profile mid-revival rather than at a clean image plane;
   the multi-lobe pattern visible here suggests a longer length might land
   on one.

**Follow-up submitted, not yet returned (job 12449008):** Ww=10.0,
`C_TRANS=12` (matching the field-plot's doubled resolution), 9-length sweep
from 34.8 to 81.2 µm — both a rough C_TRANS convergence check (does
doubling the transition slice count move totalT/balance materially at a
similar length to the original grid's 45.88 point?) and a search for a
clean revival further out than the original grid's window covered. **This
does not retract the closure verdict above** — both taper-profile shapes
still fail at every point measured so far — but it means the closure should
be read as "no working point found in the grids run so far," not
"exhaustively shown impossible," until this job lands.

**Plan, contingent on 12449008 finding a real design point:** promote the
waist profile from a `examples/verification/python/` sweep script into a
proper independent example, mirroring step 6's original scope for the
straight 1x4 (`mmi_1x4_eme.py`/`mmi_1x4_hybrid.py`) —
`examples/mmi_1x4_tapered/python/mmi_1x4_tapered_eme.py` and
`mmi_1x4_tapered_hybrid.py`, with their own `Makefile` target (see step 6's
own note that `mmi_4x4_90hybrid` has none, which is how that device stayed
unrun). Not started; this is a plan, not a commitment, until a design point
actually survives the length-extension check.

### EME field plots, both geometries ✅ **MEASURED 2026-08-21** (job 12121984)

`examples/verification/python/plot_mmi1x4_straight_vs_tapered_eme_fields.py` —
new, and it renders **both** geometries from one config so the only difference
between the two runs is the body taper itself. Output under
`example_output/mmi_1x4/mmi_1x4_straight_vs_tapered_eme/fields/`.

**Why a new script.** The existing `plot_mmi1x4_eme_fields.py` (2026-08-07) is
pinned to the **superseded** geometry — `MM_LENGTH=206.0`,
`WAVEGUIDE_OFFSET=2.4`, `res=40`, all three replaced by step 9c (2026-08-14)
and the res 40→30 change (2026-08-18) — and had **no saved output anywhere**,
so there was no straight-body field plot of the *current* device to compare
against either. Two separately-pinned scripts also could not guarantee that the
only difference is the taper, which is the entire point of the comparison.

**Result.** Both geometries reach clean 4-lobe self-imaging at the output, with
the lobes landing on the real port positions (x ≈ ±1.55, ±4.65 — matching
`output_ys`). The straight body steps to full 12.0 µm at the body face
(z ≈ 5.5 µm); the tapered body flares gradually across its 8 µm taper band
(z ≈ 5.5 → 13.5 µm) before reaching full width. This is the visual counterpart
to step 3's numeric evidence that a tapered body genuinely propagates
differently.

**Cost, measured:** 559 s (straight) + 557 s (tapered) at res=30, 12 body
modes, 24 body cells = 450 mode-solves per geometry. A smoke pass at res=15 /
4 modes / 6 body cells ran both in **73 s** total and is what the real pass was
scaled from — the script keeps the env-var overrides for it
(`MMI1X4_FIELD_RES` etc.).

> **Component naming — worth knowing before reading these plots.** The
> filenames use **meow's** component names, not the device frame's. meow's
> `Ex` is the device's transverse-horizontal field, i.e. the dominant TE
> component: **measured |Ex| ≈ 34 vs |Ey| ≈ 0.63**, a factor of 54. So
> `eme_Ex_*.png` is the propagating field; `eme_Ey_*.png` is the
> sidewall-discontinuity component, which is not the field of interest but
> *does* trace the body **outline** — a genuinely useful independent check
> that the taper band sits where intended (the width step at z ≈ 12 µm is
> plainly visible in it). The script's own closing message originally pointed
> readers at `Ey`; corrected, with the measured magnitudes recorded in it.

**Two caveats, since these images are manuscript-adjacent.** The taper profile
is the placeholder above, *not* a designed device — step 6 picks the real
values. And `BODY_CELLS=24` was chosen for a legible picture, **not** from a
convergence study; step 5 is what establishes the converged count, and this
script deliberately does not pre-empt it. Both geometries use the *same* cell
count so the comparison isolates geometry rather than discretisation.

---

## Step 5 — Convergence

Two axes, both new for this device:

- **Cell count** across the taper. Sweep it and show the S-matrix plateaus.
- **`num_modes` per cell.** The body widens, so the guided-mode count changes
  along the taper; the widest cell sets the requirement. Note the 4×4's warning
  here (`mmi_4x4_90_hybrid_validation.md`): *"every width 4.0–10.0 µm hits the
  full n_request=24 cap — this probe only ever measured a floor, never found a
  ceiling."* Do not repeat that mistake; raise the cap until a real ceiling
  appears.
- **Untapered limit**, the cheapest and most diagnostic check: body taper → 0
  must reproduce the straight 1×4's answer.

`beef/interface_modes_convergence.py` is the tool for the mode-count axis
(README's "Layer 2" completeness check) and `beef/cut_position_invariance.py`
for cut placement ("Layer 3"). Use them rather than hand-rolling — that is what
they are for.

---

## Step 6 — The device and its scripts

- `examples/mmi_1x4_tapered/python/mmi_1x4_tapered_eme.py` and
  `_hybrid.py`, following `mmi_1x4_hybrid.py`'s structure.
- Whether a whole-device FDTD reference is affordable depends on the final
  length. **Note the straight 1×4's stale justification:**
  `mmi_1x4_sin300_validation.md:1034` drops `mmi_1x4_fdtd.py` because *"total
  device length ≈216 µm makes a single-domain, whole-device FDTD simulation
  impractical"* — but `MM_LENGTH` became **39.90 µm** on 2026-08-14
  (`mmi_1x4_hybrid.py:157`), so that reasoning is stale by 4× and a full-FDTD
  check 4 may now be affordable for both the straight and tapered devices. Re-ask
  the question rather than inheriting the answer.
- Add a `Makefile` target; note that `mmi_4x4_90hybrid` has none, which is how
  that device stayed unrun.

---

## Step 7 — The four checks

Standard, per [README](README.md#the-four-checks). The tapered body makes check 1
(engine agreement on β through the shared cross-section) more interesting than
usual: there is no single body cross-section, so decide which one the check is
asked about — probably the uniform middle, with the taper cells covered by step 5's
convergence instead.

---

## Length payoff — the point of the device

**ARGUED, to be measured in step 6.** The straight 1×4 sits at `mm_length = 39.90`
µm (`mmi_1x4_hybrid.py:157`) with `mm_width = 12.0`. A tapered body reaches the
same modal beat length in less distance by running narrower where the field does
not need the full width. Quantify: tapered total length vs. 39.90 + 2×5.0 taper,
at matched imbalance and excess loss. If the saving is not material, the device is
still worth building — it is what makes the manuscript's tapered-body claim real —
but say so plainly rather than overselling the footprint.

---

## Manuscript dependency

`beef_manuscript/sections/simulation_paradigm.tex:53`'s final sentence —

> "A body spanning more than one internal cross-section (such as a tapered
> multimode region of a MMI) still requires the mode field profiles for EME
> propagation, and computes them as before."

— is written as though steps 1–3 are done. It is currently **false**: there is no
multi-cross-section path in `solve_section`. Tracked in
`beef_manuscript/unimplemented_changes.md`. `sec:1x4_mmi_splitter` also depends on
this device existing.
