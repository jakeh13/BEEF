# Merging `section_mode_beta()`'s per-length FDTD sims into one multi-monitor simulation — plan

**BUILD PLAN — the mechanism this document specs does not exist yet.**
Prompted by a direct question during manuscript-figure work on the 2×2 MMI's
`calibrate_body_neff` geometry (`calibrate_body_neff_mmi_2x2_geometry.ipynb`),
2026-08-12: `FDTDSolver.section_mode_beta()` fits body β from `len(lengths)`
independently-simulated FDTD sections — for the MMI 2×2 body specifically,
the default is 3 lengths × `N_MODES=14` excited bands each = up to **42
separate `mp.Simulation` builds** for one calibration. Three questions
prompted this: (1) is that really the count — confirmed, yes; (2) has the
choice of 3 lengths ever been convergence-checked — checked directly against
`validation_markdown/`, no; (3) could the per-length sims be merged into one
longer domain with multiple monitors instead — this document specs how, and
makes the previously-prohibitive convergence check (2) nearly free as a side
effect of doing (3).

## Status — paused 2026-08-12, resume here on the cluster

Steps 1-3 done locally (M3 MacBook, no cluster needed for these — pure
design/implementation/smoke-test, not real production-scale compute).
Steps 4-8 need real FDTD time at a scale this machine can't reasonably
absorb — paused here deliberately, not blocked on anything unresolved.

**Done:**
- Step 1: real per-sim cost measured (3609.53 s, one band). Old-path
  (42.1 hr) / merged-path (14.0 hr, projected) numbers below.
- Step 2: multi-monitor layout designed — formula, boundary check, and the
  source-proximity concern resolved with real numbers, all below.
- Step 3: `FDTDSolver.solve_section_multi_length()` implemented
  (`beef/solvers/fdtd.py`) and smoke-tested (structural correctness only).

**Not done — pick up here:**
- Step 4 (numerical equivalence test) is next. Real cost estimate at a
  *reduced* scale (`n_modes=2`, 1 wavelength, real `res=40`/
  `decay_tol=1e-11`/real lengths — NOT full `N_MODES=14`/11-wavelength
  production scale, which would cost tens of hours just for this check):
  old path (3 lengths × 2 bands ≈ 6 sims) **≈ 3.5-4 hours**; merged path
  (2 bands × 1 longer domain) **≈ 1.2-1.5 hours**; **≈ 5-6 hours
  combined**. The time-stepping/decay cost doesn't shrink with fewer
  modes — only step 1's measured output-decomposition tail (~26 of the 60
  min) does, which is why this is still hours, not minutes, even reduced.
  Steps 5-8 unstarted, all real-FDTD-time-gated, all described in full
  below.
- Both new scripts are permanent, git-tracked, NOT in any session
  scratchpad: `examples/verification/python/
  time_section_mode_beta_single_sim.py` (step 1) and `examples/
  verification/python/smoke_test_solve_section_multi_length.py` (step 3).
  Re-run the smoke test first on the cluster env before step 4, as a cheap
  sanity check that nothing broke in transit.

## The current mechanism, checked directly against the source

`section_mode_beta()` (`beef/solvers/fdtd.py:2049` on) picks 3 lengths
(`[spacing, 2·spacing, 3·spacing]`, spacing set by a cheap MPB-only seed —
see its own docstring) and loops:

```python
for j, L in enumerate(lengths):
    T = self.solve_section(geometry, wavelengths, mid - L/2, mid + L/2,
                            n_modes_in=n_modes, n_modes_out=n_modes,
                            save_dir=save_dir, cache_only=cache_only, **iface, **geometry_kwargs)
    phases[j] = np.angle(np.diagonal(T, axis1=1, axis2=2))
```

Each `solve_section()` call independently: builds a fresh domain sized to
`[x0, x1]` (`fdtd.py:2570-2580`, `mp.Simulation(...)`), excites each of
`n_modes_in` bands in its own separate `mp.Simulation` (one build + one
`sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol,
...))` per band, `fdtd.py:2561-2657`), and reads exactly **two** mode
monitors per run: `norm_flux` near the source, `mea_flux` at the far
(interface) end (`fdtd.py:2581-2590`). Three lengths means three full domain
rebuilds and three independent decay-to-tolerance waits per excited band,
even though the underlying cross-section is identical across all three (only
the propagation distance changes) — the same "solve something that doesn't
need re-solving" shape the `calibrate_body_neff` mode-solve-skip work
(`resolved/calibrate_body_neff_mode_solve_skip_validation.md`) already fixed
one layer up, here one layer lower.

## Why merging is plausible: the dominant cost doesn't scale with monitor count

`sim.run(until_after_sources=mp.stop_when_dft_decayed(...))` is gated on DFT
convergence — the field has to decay to `decay_tol` at the monitors present,
which is governed by domain length, PML absorption, and the structure's own
group velocity/loss, not by how many DFT-accumulating mode monitors are
attached. `add_mode_monitor()` is a passive DFT probe (no back-action on the
simulated field) — meep already runs multiple monitors in one `Simulation`
elsewhere in this codebase (`norm_flux` + `mea_flux` together, every call).
Nothing here is new meep capability; the change is architectural: build ONE
domain sized to the *longest* length, place a `mea_flux` monitor at each of
the (currently 3, to be convergence-checked) length positions as *interior*
points rather than at the domain's own edge, excite each band once, and read
`get_eigenmode_coefficients()` off every monitor after the single
decay-to-tolerance wait. Projected: 3 domain rebuilds × 14 bands (42 sims)
→ 1 domain × 14 bands (14 sims) for the MMI body specifically — not
assumed, to be measured for real in step 1 below.

## Checked, not assumed: what could break this

- **PML/domain sizing currently keys off a single `x1` per call.** The
  domain must be built for the *longest* requested length, with shorter
  lengths' monitors sitting inside that domain rather than at its edge —
  `section_fdtd_layout()` (`fdtd.py:613`) computes source/monitor
  positions from a single `x0,x1` pair today; a multi-monitor variant needs
  its own layout logic, not a bolt-on. **Not designed yet** — step 2 below.
- **`solve_section()`'s public contract must not change shape/semantics
  for its many other callers.** `T: [n_wl, n_out, n_in]`, cache-keyed on a
  single `x0,x1` — used throughout this codebase (`HybridSolver`,
  `verify_*` scripts, every hybrid example). A multi-monitor length sweep
  returns `T: [n_L, n_wl, n_out, n_in]` instead, so this has to be a
  **new, dedicated method** (or an internal helper `section_mode_beta()`
  calls directly) — not an edit to `solve_section()` itself. Confirmed by
  reading every `solve_section()` call site's expected return shape; none
  tolerate an extra leading axis silently.
- **The shortest monitor's proximity to the source.** Today's shortest
  length already has to clear `dpml` + whatever near-field settling the
  source needs (implicit in each independently-built domain's own margins).
  Merged into one domain, the shortest monitor sits closer to the source in
  absolute terms than its own independently-built domain would have placed
  it relative to the (now-longer, shared) domain's own PML — needs an
  explicit minimum source-to-monitor clearance check, analogous to the
  existing `_check_aperture()` PML-clipping guard (`fdtd.py:2541-2556`), not
  assumed safe by extrapolation from the single-length case.
- **The linearity/residual self-check's statistical independence.**
  Today's 3 separate sims give 3 fully independent noise realizations to
  fit a line through and check residuals against
  (`_BETA_FIT_RESIDUAL_MAX`/branch-margin guards). One shared simulation
  correlates that noise across the 3 (now merged) monitors — the guards
  might read artificially cleaner than they should. Given `decay_tol=1e-11`
  is already tight, likely a second-order effect, but **not measured**,
  flagged for step 4's equivalence check to actually look at, not assume
  away.
- **Cache-key fingerprinting.** A merged-sim run's natural cache
  granularity (one run, N lengths) doesn't match today's per-length
  `solve_section()` cache entries (`build_run_config()` keyed on a single
  `x0,x1`) — old per-length cache entries and a new merged run would not be
  cache-interchangeable, and a fresh design is needed for what the merged
  method's own cache key covers (see `beef/cache/config.py`'s
  `_METHOD_ATTRS` docstring for the safe-default-direction reasoning to
  follow here too).

## Scope — where this plan does *not* claim to apply

- Any caller of `solve_section()` directly (not through
  `section_mode_beta()`) — unaffected; `solve_section()` itself is
  untouched by this plan.
- Devices/examples that don't use `calibrate_body_neff` — nothing to speed
  up there; this plan only touches `section_mode_beta()`'s internal cost.
- The EME-side mode-solve-skip work
  (`resolved/calibrate_body_neff_mode_solve_skip_validation.md`) — a
  separate, already-closed optimization one layer up in the same pipeline;
  this plan is about the FDTD measurement `calibrate_body_neff` depends on,
  not the EME propagation it feeds.

## Implementation steps

1. **Measure first — same discipline as the mode-solve-skip plan.** Time
   one real `solve_section()` call (one excited band, one length, the MMI
   2×2 body's actual calibration window) to get a real per-sim cost number,
   then project the 42→14 (or whatever the final monitor count turns out to
   be) savings against it for real, not by assumption.

   **Measured, 2026-08-12.** `n_modes_in=1` (one excited band), the real
   calibration window from `calibrate_body_neff_mmi_2x2_geometry.ipynb`
   (`mid=0.0`, `L=2.2046 µm`, `n_modes_out=14`, all 11 wavelengths,
   `decay_tol=1e-11` — identical to production), local M3 MacBook, no
   cache. **3609.53 s (≈60.2 min) for this single sim** — the earlier,
   pre-measurement guess a few turns back ("tens of minutes" for the whole
   calibration) badly understated this; one band alone is an hour. Profile
   breakdown: the time-stepping phase (decay to `1e-11`) took ~34 min
   (steps 1-30063, ~0.06-0.09 s/step, no deceleration visible — genuinely
   decay-limited, not stalled), and the remaining ~26 min was the
   *output-side* decomposition — `get_eigenmode_coefficients()` for each of
   the 14 output bands is its own MPB eigenmode solve × 11 wavelengths
   (154 small solves total), not free as the "cheap, MPB-only" language
   elsewhere in this codebase might suggest for a *single* such call — it
   adds up at `n_modes_out=14`.

   Projected from this one real data point (script:
   `examples/verification/python/time_section_mode_beta_single_sim.py`):
   - **Old path** (3 lengths × 14 bands = 42 sims): **151,600 s ≈ 42.1
     hours**.
   - **Merged path, projected** (14 bands, 1 domain each): **50,533 s ≈
     14.0 hours** — assumes per-band cost at the longest merged domain
     equals this single-length measurement; **not yet verified** (the
     merged domain is longer than this L=2.2046 µm window, so real
     per-band cost there could be higher, e.g. from a longer decay time —
     step 3's real implementation is what actually measures this, this
     projection is a first-order estimate only).
   - **Projected speedup: 3.00x** — before accounting for whether the
     longer merged domain's own per-band cost changes that ratio.

   This resets the scale of the whole plan: this isn't a "nice speedup on a
   minutes-long calibration," it's the difference between ~42 hours and
   ~14 hours of real FDTD time for one device's calibration — both numbers
   far beyond anything reasonable to run outside a cluster, which also
   raises the bar on how carefully steps 4-6 need to verify correctness
   before this becomes the default (a 3x-wrong number that takes 14 hours
   to produce is much worse to discover late than a fast one).
2. **Design the multi-monitor domain/layout.**

   **Designed, 2026-08-12** (not yet implemented — that's step 3), by
   working through `section_fdtd_layout()`'s actual formulas
   (`fdtd.py:613-719`) rather than assuming a shape.

   **The key simplification**: the merged design does *not* need to
   reproduce the old per-length `[x0, x1]` windows exactly. It only needs
   each monitor's *distance from the shared source* to equal the old
   `lengths[i]` values — the phase fit only ever consumed
   `(distance, phase)` pairs, never the windows themselves. That
   reformulation is what makes one shared domain work at all.

   For `section_mode_beta()`'s actual calibration calls, `n_in == n_out ==
   N_MODES` always, so `section_fdtd_layout()`'s own `reverse = n_out <=
   n_in` is always `True` — source sits at the domain's *right* edge
   (`src_x = L - source_margin`), interface monitor at the *left* edge
   (`mea_x = 0`), confirmed against the actual rendered geometry plots in
   `calibrate_body_neff_mmi_2x2_geometry.ipynb` ("source (reverse)" label
   on the right, "interface mon" on the left).

   Formula, `reverse=True` case (the case that actually matters here;
   `reverse=False` is the mirror image — swap left/right):
   - `L_domain = max(lengths)` — the merged domain's local length, playing
     the role today's single `L` plays.
   - `src_x = L_domain - source_margin`, `norm_x = L_domain -
     monitor_margin` — **both fixed, computed once**, same formulas as
     today with `L_domain` substituted for `L`.
   - `mea_positions[i] = L_domain - lengths[i]` — one interface monitor
     per requested length.

   **Why this is exactly right, not approximately**: distance from
   `norm_x` to `mea_positions[i]` = `(L_domain - monitor_margin) -
   (L_domain - lengths[i])` = `lengths[i] - monitor_margin` — the *same*
   `length - monitor_margin` distance today's single-length case already
   produces (its `mea_x=0` gives distance `norm_x - 0 = L -
   monitor_margin`). The merged design reproduces today's per-length
   measurement geometry exactly, just embedded in one shared domain,
   rather than approximating it.

   **Boundary check**: at `lengths[i] = L_domain` (the longest length),
   `mea_positions[i] = 0` — exactly matches today's single-longest-length
   case. Consistent by construction, not by coincidence.

   **Resolves the "shortest monitor's proximity to source" concern from
   above, with real numbers**: for the MMI 2×2 body's actual calibration
   (`lengths = [1.1023, 2.2046, 3.3069]` µm, `L_domain = 3.3069`,
   `source_margin = 0.0` for this device), the *shortest* length's monitor
   sits at `mea_positions = 3.3069 - 1.1023 = 2.2046` µm — distance from
   the source = `2.2046 - 0.0 = 2.2046` µm, not close to the source at
   all. The monitor closest to the source is actually the shortest
   length's (largest `mea_positions[i]`, since it's `L_domain -
   lengths[i]` and shorter lengths give larger values) — still `spacing -
   source_margin` away, not zero. Reasoned through with this device's real
   numbers, not just algebraically; still needs an explicit runtime
   assertion in the step 3 implementation (an `_check_aperture()`-style
   guard per monitor, not just per source/single-monitor as today), not
   left as "should be fine."
   - **No new data structure needed, found during implementation (step 3):**
     the plan here originally called for a `MultiSectionFDTDLayout`
     dataclass. Turned out unnecessary — `section_fdtd_layout()` itself is
     reused completely unchanged (called once, with `x0,x1` spanning
     `max(lengths)`); only `mea_positions` is computed locally in the new
     method from `layout.L`/`layout.reverse`, per the formula above. A
     simpler design than planned, corrected here rather than left stale.

3. **Implement the merged sim.**

   **Implemented, 2026-08-12** — `FDTDSolver.solve_section_multi_length()`
   (`beef/solvers/fdtd.py`, right after `solve_section()`). One
   `mp.Simulation` per excited band (matching today's per-band loop,
   `solve_section()`'s own `fdtd.py:2561`-equivalent), one `norm_flux`
   monitor plus N `mea_flux` monitors (one per length) via
   `add_mode_monitor()`, one shared `sim.run()`, then N × `n_measure`
   `get_eigenmode_coefficients()` reads afterward. New method — does not
   modify `solve_section()`'s existing contract, confirmed by construction
   (a separate method, `solve_section()` untouched). Deliberately narrower
   than `solve_section()`: no `PortSpec`/`y0`/`y1` support (unneeded by
   `section_mode_beta()`'s calls), no `save_dir`/`cache_only` (cache-key
   shape is still an open question per this doc — deferred, not silently
   decided). Returns `T: [len(lengths), len(wavelengths), n_modes,
   n_modes]`.

   **Smoke-tested, 2026-08-12** (structural correctness only — real
   equivalence is step 4, not this): low `res=10`, loose `decay_tol=1e-3`,
   `n_modes=2`, 3 short lengths, MMI body. Ran to completion in 3.2s, no
   Python errors, correct output shape, no NaNs, and — the important
   structural check — all 3 length slices returned genuinely distinct
   complex values (confirms the 3 monitors are actually reading different
   positions, not aliasing onto the same one). Script:
   `examples/verification/python/smoke_test_solve_section_multi_length.py`
   (re-run from this permanent location, 2026-08-12, confirmed identical
   pass: shape `(3, 1, 2, 2)`, no NaNs, all 3 slices distinct, 3.1 s).
4. **Numerical equivalence test, before this ships as a default.** Same
   standard as the mode-solve-skip plan: run both the merged approach and
   today's 3-independent-sims approach on the real MMI 2×2 body, compare
   fitted β (and the raw per-length phases) to floating-point precision or
   document the real residual if merging's correlated-noise concern (above)
   turns out to matter. Also explicitly check the linearity/branch-margin
   guards still fire correctly on a case constructed to need them (don't
   just check the happy path).
5. **Convergence study #1 — number of monitors/lengths.** The gap this
   whole document exists to close: with merging, adding a monitor is cheap
   (no new domain build/decay wait), so sweep N ∈ {3, 4, 5, 6, ...} monitors
   in one merged run and check whether the fitted β stabilizes, and how the
   linearity-residual guard's sensitivity changes with N. This was
   prohibitively expensive to check under the old N-separate-sims
   architecture (each extra length was another full sim) — merging is what
   makes this checkable at all, not just cheaper.
6. **Convergence study #2 — total simulation/domain length.** Does the
   longest monitor's distance need to keep growing — e.g., do the *shorter*
   lengths' fitted phases change once they're measured inside a longer
   shared domain (near-field/launch-transient contamination) versus their
   own independently-sized domain today? Sweep the total domain length
   (independent of monitor count) and check β-fit stability against it.
7. **Wire in with a fallback, same pattern as `SYMMETRIC_MODE_REUSE`.**
   Once steps 4-6 pass, make the merged approach `section_mode_beta()`'s
   default, keep the old per-length path reachable behind a flag (mirrors
   `mmi_2x2_hybrid.py`/`directional_coupler_hybrid.py`'s own
   `SYMMETRIC_MODE_REUSE` toggle from the same session) for re-validation if
   a future finding needs the fully-independent-sims version again.
8. **Re-run the shipped examples' calibration for real.** `mmi_2x2_hybrid.py`
   and `directional_coupler_hybrid.py`/`directional_coupler_soi160_hybrid*`
   all call `calibrate_body_neff` through `section_mode_beta()` — their
   calibration numbers should be regenerated (cluster, given this is still
   real FDTD compute even after merging) once the merged path is validated,
   to confirm the shipped β values don't move outside noise.

## Open questions

- Does the merged approach's cache-key design (flagged above, not designed)
  want its own tier-style reuse (cf.
  `EMESolver.solve_section_length_sweep()`'s 3-tier cache,
  `resolved/calibrate_body_neff_mode_solve_skip_validation.md`'s extension
  session) — e.g. a length-independent fingerprint so a later call asking
  for a different/longer set of monitor positions on the same cross-section
  doesn't force a full re-run? Not designed, worth considering once step 2
  lands.
- Is 3 monitors (today's default) actually going to turn out to be enough
  once convergence study #1 has real data, or does this reveal the shipped
  calibration numbers (DC's `job 11517245`, etc.) were under-sampled all
  along? Genuinely unknown until measured — this document takes no
  position on the answer, only that the question was never asked before.
- Does the correlated-noise concern on the linearity/residual guard (flagged
  above) turn out to matter in practice, or is `decay_tol=1e-11` tight
  enough that it's inert? Step 4 answers this for real; not assumed either
  way here.
