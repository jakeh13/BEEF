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
| A **tapered-width body** in the `MMI` geometry | ❌ **does not exist** — `MMI`'s body is a single constant `mm_width` (`geometry/mmi.py:_section_polygons`, `_rect(max(lo, in1), min(hi, out0), W / 2)` with `W = self.mm_width`) |
| **Multi-cell EME inside a cascade section** | ❌ **does not exist** — `solve_section()` builds one cell, both branches |
| `calibrate_body_neff` semantics for a non-uniform body | ❌ **undefined**, and currently unguarded |

Ordered work: **guard the silent-wrong-answer path (step 1) → tapered-body
geometry (step 2) → real multi-cell `solve_section` (step 3) → decide calibration
semantics (step 4) → convergence + untapered-limit check (step 5) → the device and
its four checks (steps 6–7).**

Step 1 is small and lands first on purpose: until it does, anyone pointing a
calibrated hybrid at a tapered body gets a plausible wrong answer with no warning.

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

## Step 1 — Refuse a non-uniform span instead of silently diagonalising it

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

---

## Step 2 — Tapered-body geometry

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

---

## Step 3 — Real multi-cell EME inside a cascade section

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

---

## Step 4 — What `calibrate_body_neff` means for a multi-cell body

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
