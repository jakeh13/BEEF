# Merging `section_mode_beta()`'s per-length FDTD sims into one multi-monitor simulation — plan

**COMPLETE — the mechanism this document specs is built, validated, and wired in (steps 1-8, 2026-08-12 → 2026-08-16).** It is not yet the path production actually takes: every real caller passes `save_dir`, and the merged path has no cache-key design, so all of them still fall back to the per-length path. See the Status section and the open questions at the end.
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

## Status — COMPLETE: steps 1-8 all done and measured

Steps 1-3 done locally (M3 MacBook, no cluster needed for these — pure
design/implementation/smoke-test, not real production-scale compute).
Steps 4-6 done for real on the cluster, 2026-08-13/14. Step 7 (wiring) done
locally, 2026-08-14 — a code change plus tests, no cluster needed. Step 8
(jobs 11872191/11872184/11872185) landed 2026-08-15/16, all COMPLETED,
exit 0 — full production scale, ~42 cluster-hours across the three.

**Done:**
- Step 1: real per-sim cost measured (3609.53 s, one band). Old-path
  (42.1 hr) / merged-path (14.0 hr, projected) numbers below.
- Step 2: multi-monitor layout designed — formula, boundary check, and the
  source-proximity concern resolved with real numbers, all below.
- Step 3: `FDTDSolver.solve_section_multi_length()` implemented
  (`beef/solvers/fdtd.py`) and smoke-tested (structural correctness only).
- Step 4: numerical equivalence test **PASSED for real, 2026-08-13**, job
  11813378 (`gts-sr108-ece`, 07:08:01 wall, exit 0). Reduced scale
  (`n_modes=2`, 1 wavelength, real `res=40`/`decay_tol=1e-11`, real
  auto-picked lengths `[1.1023, 2.2046, 3.3069]` µm — NOT full
  `N_MODES=14`/11-wavelength production scale). Full numbers, guard checks,
  and script/batch record below at step 4 itself.
- Step 5: convergence study #1 (monitor count) **done, 2026-08-14**, job
  11838953 — β converged at N=3, residual grows with N but the writeup's
  original "confounded with domain length" explanation was wrong (caught
  and corrected in place) and got resolved for real by step 6.
- Step 6: convergence study #2 (domain length) **done, 2026-08-14**, job
  11844627 — found a real, non-monotonic domain-length sensitivity on the
  same 3 fit lengths (β range `5.78e-04` rad/µm across 3×-8×spacing
  domains, ~15x step 4's equivalence noise floor), still 3+ orders of
  magnitude under both guard thresholds. Full numbers at step 6 itself.
  - Scripts: `examples/verification/python/
    equivalence_test_solve_section_multi_length.py` (step 4),
    `convergence_study_monitor_count.py` (step 5),
    `convergence_study_domain_length.py` (step 6). Batches:
    `batch/verification/run_equivalence_section_multi_length.batch`,
    `run_convergence_monitor_count.batch`, `run_convergence_domain_length.batch`.
- Step 7: wired in with a fallback, **done, 2026-08-14** — new
  `use_merged_domain=True` default on `section_mode_beta()`
  (`beef/solvers/fdtd.py`), automatic fallback to the per-length path when
  caching is requested (today, every real production caller), a regression
  found and fixed in `tests/test_beta_calibration.py`, and Guard 2 now
  covers the merged path for free. Full record at step 7 itself.
- Step 8: shipped examples re-run at full production scale, **done,
  2026-08-15/16** — jobs 11872185 (DC SOI160, 1h51m), 11872184 (DC hybrid,
  3h31m), 11872191 (MMI 2×2, 31h23m), all exit 0. Both DC scripts came back
  bit-for-bit identical to their prior runs; the MMI 2×2 had no valid
  pre-step-7 baseline and was confirmed inert by a build-count measurement
  instead (44 builds = the per-length path exactly). Full record at step 8.
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

   **PASSED, 2026-08-13.** Job 11813378 (`gts-sr108-ece`, 07:08:01 wall,
   exit 0). Reduced scale: `n_modes=2`, 1 wavelength (`1.5475` µm, band
   centre), real `res=40`, real `decay_tol=1e-11`, lengths auto-picked the
   same way `section_mode_beta()` picks them internally (not hand-copied):
   `[1.1023, 2.2046, 3.3069]` µm (`beta_max=11.4001` rad/µm,
   `max_safe_spacing=2.2046` µm). Old path 3 independent `solve_section()`
   calls, 18297.8 s (5.08 h) total (4995.7 / 5724.4 / 7577.6 s per length,
   growing with length as expected). Merged path one
   `solve_section_multi_length()` call, 7339.3 s (2.04 h). **Measured
   speedup at this reduced scale: 2.49x** (doc step 1's full-production-scale
   projection was 3.00x — in the same range, not identical, since the merged
   domain here is longer than any single old-path domain and per-band cost
   isn't perfectly flat with domain length, exactly the caveat step 1 itself
   flagged as unverified).

   Agreement, old vs. merged:
   - **Raw phases**: `max|phase_old - phase_merged| = 7.72e-05 rad` — about
     1300x tighter than `_BETA_FIT_RESIDUAL_MAX` (0.1 rad), the guard's own
     straight-line tolerance, used here only as a scale reference.
   - **Raw `|T_ii|`**: `max||T_ii|_old - |T_ii|_merged|| = 2.64e-03`.
   - **Fitted β** (`beta_fit_from_phases()`, same seed fed to both):
     `max|beta_old - beta_merged| = 3.50e-05 rad/µm` (**0.0003%**), both
     wl/mode cells' residuals ≤ 8.00e-05 rad (guard fires above 0.1) and
     margins 0.060 (guard fires above 0.9) — nowhere near either threshold
     on the happy path, for both paths.
   - **The correlated-noise concern this doc flagged above turned out not to
     matter** at this decay_tol: agreement is 3-4 orders of magnitude
     tighter than the guards themselves, not merely "inside tolerance."

   Guard checks, all real, not assumed:
   - **B1** (pure-math sanity, no FDTD): `beta_fit_from_phases()` recovers a
     known synthetic β=9.4 rad/µm at these real lengths to `1.78e-15`
     — regression baseline for the shared fit function both paths rely on.
   - **B2** (residual guard, reusing this run's real phases, no extra FDTD):
     a deliberate 1.0 rad outlier injected into one (wl, mode) cell of each
     path's real phases pushes `resid` to 0.6667 rad on **both** — guard
     fires correctly, consistently, on data from either path.
   - **B3** (branch-margin guard / Guard 2, cheap — seed solve only, no
     time-stepping): `section_mode_beta()` called for real with deliberately
     too-widely-spaced lengths (`[4.4092, 8.8185]` µm against
     `max_safe_spacing=2.2046` µm) raised `ValueError` exactly as designed
     ("branches only 1.425 rad/µm apart on a β of ~11.400 ... Space them no
     more than 2.205 µm apart"). This only exercises the *old* path's wired
     guard — `solve_section_multi_length()` has no guard wrapper of its own
     yet (step 7 wiring, not this step) — flagged in the script's own
     output, not silently assumed equivalent.

   Script: `examples/verification/python/
   equivalence_test_solve_section_multi_length.py`. Batch:
   `batch/verification/run_equivalence_section_multi_length.batch`
   (`--account=gts-sr108-ece`, `-t 600`). Full log:
   `reports/11813378.out`.
5. **Convergence study #1 — number of monitors/lengths.** The gap this
   whole document exists to close: with merging, adding a monitor is cheap
   (no new domain build/decay wait), so sweep N ∈ {3, 4, 5, 6, ...} monitors
   in one merged run and check whether the fitted β stabilizes, and how the
   linearity-residual guard's sensitivity changes with N. This was
   prohibitively expensive to check under the old N-separate-sims
   architecture (each extra length was another full sim) — merging is what
   makes this checkable at all, not just cheaper.

   **PASSED (β side), residual trend flagged, 2026-08-14.** Job 11838953
   (`gts-sr108-ece`, 03:11:17 wall, exit 0). Same reduced scale as step 4
   (`n_modes=2`, 1 wavelength, real `res=40`/`decay_tol=1e-11`), N sweep
   {3,4,5,6} monitors from ONE `solve_section_multi_length()` call
   (spacing=1.1023 µm shared across all N, domain = 6×spacing = 6.6138 µm,
   well inside the 26.5 µm usable body span — nowhere near the tapers).
   Every N < 6 is a free re-fit of a prefix of that one run's phase data
   via `beta_fit_from_phases()` — no extra simulation.

   Fitted β (rad/µm, both (wl, mode) cells):

   | N | β mode 0 | β mode 1 | worst resid (rad) | worst margin |
   |---|----------|----------|--------------------|--------------|
   | 3 | 11.571207 | 11.407468 | 1.069e-04 | 0.0603 |
   | 4 | 11.571265 | 11.407524 | 1.283e-04 | 0.0603 |
   | 5 | 11.571236 | 11.407497 | 1.289e-04 | 0.0603 |
   | 6 | 11.571113 | 11.407379 | 4.503e-04 | 0.0603 |

   **β is already converged at N=3** — `max|β(N) - β(N=6)|` across all four
   N values stays at `9.4e-05`-`1.5e-04` rad/µm regardless of N (not a
   monotonic decay toward N=6, just noise-floor fluctuation around the same
   value): today's production default of 3 monitors gives the same β, to
   this precision, as 6 monitors would. Successive-N steps
   (`|β(N)-β(N-1)|`): `5.8e-05` (4v3), `2.9e-05` (5v4), `1.2e-04` (6v5) —
   also non-monotonic, consistent with noise rather than a real trend
   needing more points to average out.

   **The residual, unlike β, DOES show a real, monotonically-growing trend
   with N** — `1.07e-04 → 1.28e-04 → 1.29e-04 → 4.50e-04` rad from N=3 to
   N=6, a ~4.2x jump overall (biggest single step N=5→6). Still 2-3 orders
   of magnitude under `_BETA_FIT_RESIDUAL_MAX` (0.1 rad) — the guard would
   not have fired at any N tested — but this is a measured trend, not
   noise.

   **Correction, 2026-08-14** (caught while designing step 6, not left
   stale): this entry originally attributed the residual growth to N and
   domain length being confounded in this run (`domain = N×spacing`
   growing alongside N). That's wrong — job 11838953 made exactly ONE
   `solve_section_multi_length()` call with all 6 lengths at once, so the
   domain was fixed at `6×spacing = 6.6138 µm` for *every* row of the table
   above, including the N=3 row. What actually varies across N here is
   only how many of that one run's 6 already-measured monitors get fed
   into the fit — domain length never moved. The real domain-length
   question (does embedding the *same* 3 lengths in a bigger shared domain
   shift their own measured phase/β) turned out to be answerable directly
   by comparing this job's N=3 row against step 4's separate merged-path
   run (job 11813378, domain = `3×spacing = 3.3069` µm, the same 3
   lengths) — see step 6 below, which resolves this for real using that
   comparison plus new data.

   Script: `examples/verification/python/
   convergence_study_monitor_count.py`. Batch: `batch/verification/
   run_convergence_monitor_count.batch`. Full log: `reports/11838953.out`.
   **Caveat carried from step 4**: reduced scale only (`n_modes=2`, 1
   wavelength) — not yet checked at full `N_MODES=14`/11-wavelength
   production scale.
6. **Convergence study #2 — total simulation/domain length.** Does the
   longest monitor's distance need to keep growing — e.g., do the *shorter*
   lengths' fitted phases change once they're measured inside a longer
   shared domain (near-field/launch-transient contamination) versus their
   own independently-sized domain today? Sweep the total domain length
   (independent of monitor count) and check β-fit stability against it.

   **PASSED (real domain-length sensitivity found, still guard-safe),
   2026-08-14.** Job 11844627 (`gts-sr108-ece`, 09:15:11 wall, exit 0).
   Holds the 3 real production/step-4/step-5 lengths fixed and varies ONLY
   the shared domain size via one extra, ignored padding length (forces
   `L_domain` bigger, never enters the β fit) — cleanly separating "domain
   length" from "monitor count," which step 5's own run (corrected above)
   did not actually vary. New domains: pad=4×spacing (6579.0 s, 1.83 h) and
   pad=8×spacing (26704.7 s, 7.42 h — nearly 2x the ~3.9 h projected from a
   length-only cost model; a bigger domain costs more per FDTD *step* too,
   not just more steps, exactly as observed live before this job finished).
   Combined with the two free points (`3×spacing`, job 11813378;
   `6×spacing`, job 11838953's N=3 row):

   | domain (µm) | β mode 0 | β mode 1 | worst resid (rad) | worst margin | source |
   |---|---|---|---|---|---|
   | 3.3069 (3×) | 11.571281 | 11.407538 | 8.000e-05 | 0.0600 | free (job 11813378) |
   | 4.4092 (4×) | 11.571353 | 11.407607 | 2.321e-07 | 0.0604 | new |
   | 6.6138 (6×) | 11.571207 | 11.407468 | 1.069e-04 | 0.0603 | free (job 11838953) |
   | 8.8185 (8×) | 11.570775 | 11.407054 | 2.103e-04 | 0.0602 | new |

   **Real, measurable, non-monotonic domain-length sensitivity** on the
   *same* 3 fit lengths — not the noise-floor-sized wiggle steps 4-5 saw.
   β range across the sweep: `5.78e-04` rad/µm (mode 0), `5.53e-04` rad/µm
   (mode 1) — about **15x step 4's old-vs-merged equivalence agreement**
   (`3.50e-05` rad/µm) and clearly outside steps 4-5's noise floor, though
   still tiny in relative terms (`5.06e-4/11.57 ≈ 0.0044%` peak-to-peak).
   Not a clean "farther domain = worse" trend: residual is *smallest* at
   4×spacing (`2.321e-07` rad, an almost exact fit) — smaller than the
   no-padding baseline itself (`8.00e-05` at 3×) — then grows through 6×
   and 8×. β itself peaks at 4×spacing too, not at either extreme. **This
   answers step 2's original source-proximity/near-field worry for real:
   yes, embedding the same near monitors in a bigger shared domain does
   measurably shift their fitted β/residual, and the shift does not
   monotonically track domain size** — plausibly a domain-size-dependent
   interference pattern between the source's near-field and the PML's own
   (non-zero, if small) residual reflection, though this run cannot isolate
   the mechanism, only its size. Margin stayed flat (~0.060) throughout,
   same as steps 4-5 — unaffected by domain length.

   **Still deeply guard-safe**: worst resid `2.103e-04` rad is 3 orders of
   magnitude under `_BETA_FIT_RESIDUAL_MAX` (0.1 rad); worst margin `0.0604`
   is over an order of magnitude under `_BETA_BRANCH_MARGIN_RAISE` (0.9).
   Nothing here would have tripped either guard at any domain size tested.
   The practical read: today's production 3-monitor, no-padding domain
   (3×spacing) sits well within this sensitivity band, not at either
   extreme of it — no evidence the shipped calibration numbers are an
   outlier case.

   Script: `examples/verification/python/
   convergence_study_domain_length.py`. Batch: `batch/verification/
   run_convergence_domain_length.batch`. Full log: `reports/11844627.out`.
   **Caveat carried from steps 4-5**: reduced scale only (`n_modes=2`, 1
   wavelength) — not yet checked at full `N_MODES=14`/11-wavelength
   production scale, where a real (if small) domain-length effect like this
   one is worth re-checking rather than assuming it scales away.
7. **Wire in with a fallback, same pattern as `SYMMETRIC_MODE_REUSE`.**
   Once steps 4-6 pass, make the merged approach `section_mode_beta()`'s
   default, keep the old per-length path reachable behind a flag (mirrors
   `mmi_2x2_hybrid.py`/`directional_coupler_hybrid.py`'s own
   `SYMMETRIC_MODE_REUSE` toggle from the same session) for re-validation if
   a future finding needs the fully-independent-sims version again.

   **Done, 2026-08-14.** New `use_merged_domain: bool = True` parameter on
   `section_mode_beta()` (`beef/solvers/fdtd.py`). When `True` (default) and
   no caching is requested, the phase-collection block calls the merged
   `solve_section_multi_length()` once instead of looping `solve_section()`
   over `len(lengths)`; everything downstream (both post-fit guards, and
   Guard 2 above them) reads `phases` unchanged either way, so **the merged
   path now gets Guard 2's pre-simulation branch-margin check for free**,
   closing one of step 4's own open questions ("does the merged path need
   its own guard wrapper?") — yes, and wiring it through this entry point
   is how, with no new guard code needed.

   **Caching gate, not a blanket default.** Checked every real caller before
   writing this (`beef/solvers/hybrid.py`, `mmi_2x2_hybrid.py`,
   `directional_coupler_hybrid.py`, `directional_coupler_soi160_hybrid*.py`,
   `mmi_4x4_90hybrid_calibrated_baseline.py`) — every one of them passes
   `save_dir` for caching. The merged path has no cache-key design yet (open
   question above, still open), so `merged_active = use_merged_domain and
   save_dir is None and not cache_only` — any caller requesting caching
   automatically falls back to the per-length path, logged at `info` level,
   regardless of the flag. Practical consequence: **today, every real
   production caller keeps using the proven per-length path unaffected**,
   because they all pass `save_dir`. The flag is real and takes effect the
   moment a caller either drops `save_dir` or the cache-key design (open
   question above) lands — it is not dead code, but it is not yet the path
   production actually takes either. That gap is real and open, not
   papered over.

   **Regression found and fixed while wiring this in.**
   `tests/test_beta_calibration.py::TestSectionModeBetaDpmlClearance` (5
   tests) monkeypatches `solve_section()`/`section_mode_neff()` to test the
   per-length orchestration meep-free — running the suite confirmed exactly
   the predicted break: with the new default active, those tests never
   touch the fake `solve_section()` at all, and the real
   `solve_section_multi_length()` crashed instead
   (`AttributeError: 'NoneType' object has no attribute 'center_at'` on the
   test's `geometry=None`). Fixed by adding `use_merged_domain=False` to all
   5 call sites — they test the per-length path specifically, so pinning it
   preserves exactly what they were designed to check. Full suite after the
   fix: **639 passed, 19 skipped** (pre-existing skips, unrelated).

   **Smoke-tested end to end, real meep, low res=10/decay_tol=1e-3** (same
   deliberately-cheap scale as step 3's own smoke test — structural
   confirmation only, not a re-check of equivalence, which steps 4-6 already
   did at production-representative settings): default call → **3**
   `mp.Simulation` builds (1 seed + 2 excited-band merged sims); explicit
   `use_merged_domain=False` → **7** builds (1 seed + 3 lengths × 2 bands);
   default **with** `save_dir` set → **7** builds, identical count and
   bit-identical β to the explicit-False case, confirming the caching
   fallback takes the per-length path, not the merged one. (Merged-vs-old β
   disagreement at this loose smoke scale was `0.24` rad/µm — expected noise
   from `decay_tol=1e-3`/`res=10`/sub-µm lengths, not a real equivalence
   concern; the real number is step 4's `3.50e-05` rad/µm at production-
   representative settings.)
8. **Re-run the shipped examples' calibration for real.** `mmi_2x2_hybrid.py`
   and `directional_coupler_hybrid.py`/`directional_coupler_soi160_hybrid*`
   all call `calibrate_body_neff` through `section_mode_beta()` — their
   calibration numbers should be regenerated (cluster, given this is still
   real FDTD compute even after merging) once the merged path is validated,
   to confirm the shipped β values don't move outside noise.

   **ALL 3 PASSED, 2026-08-15/16.** Both DC scripts came back **bit-for-bit
   identical** to their prior runs; the MMI 2×2 had no valid pre-step-7
   baseline to diff against (an unrelated `SYMMETRIC_MODE_REUSE` commit
   landed in between), and was confirmed inert instead by a direct
   build-count measurement off its own run — details under each below. All
   three pass `save_dir`, so per step 7's fallback all three take the
   per-length path: this is a sanity re-run confirming nothing about these
   shipped numbers moved from step 7's code change, not a test of the
   merged path itself.

   Full production scale, not steps 4-6's reduced scale: `N_MODES=14`/11
   wavelengths (MMI 2×2), `N_MODES=5` (DC), `N_MODES=2` (DC SOI160), all at
   `res=40`. Real cost is far larger than anything else in this document —
   checked before submitting, not assumed: no direct cache-cold precedent
   existed for the MMI 2×2 script under real 48-way MPI parallelism
   (`LMOD_SYSHOST=phoenix` → `srun --mpi=pmi2` confirmed live), so its
   walltime is anchored on `directional_coupler_hybrid.py`'s own REAL
   cache-cold measurement under that same parallel regime (job 11609658,
   10:25:31 at `N_MODES=5`, see `batch/run_directional_coupler_L151.batch`'s
   own comment history), scaled by `N_MODES` (14/5 ≈ 2.8x) with a large
   margin on top, not a naive serial extrapolation (an earlier, serial-only
   projection put this at 42.1h — irrelevant once real MPI parallelism is
   accounted for).

   - MMI 2×2 hybrid: job 11872191, `gts-sr108-gedc_main` (see account note
     below), `-N 2 --ntasks-per-node=24`, `-t 2880` (48h). Batch:
     `batch/verification/run_step8_mmi_2x2_hybrid_recalibration.batch`.

     **PASSED on the question step 8 actually asks, but by a different
     route than the other two — no bit-for-bit comparison was possible.**
     1-07:22:45 wall (31h23m), exit 0, inside the 48h cap (and close to the
     ~29h pre-submit estimate). Results at λ=1.547 µm, `mm_length`=28.0 µm:
     bar T=0.4562 (−3.41 dB), cross T=0.5187 (−2.85 dB), total T=0.9749.

     **Why there is no before/after table here.** The only prior recorded
     run of this script is job 11791990_1 (2026-08-10), and it **predates
     commit `d7910e8` (2026-08-12)**, which introduced this example's own
     `SYMMETRIC_MODE_REUSE` — deriving `M_in0` from `M_in1` by mirror
     symmetry instead of simulating the second arm. That is a different
     change by a different mechanism, and it moves these numbers
     (11791990_1 reported bar 0.4560 / cross 0.5207 / total 0.9767). No run
     of `mmi_2x2_hybrid.py` exists between that commit and step 7 (checked
     `sacct` across the whole window — the several `beef_mmi2x2_*` jobs
     there are body-stub sweeps and body-mode field scripts, not this
     example). So the ~0.4% cross-port difference against the August 10 log
     is attributable to `SYMMETRIC_MODE_REUSE`, **not** to step 7, and this
     device simply has no clean pre-step-7 baseline to diff against. Stated
     rather than papered over with a comparison that would not mean what it
     appears to.

     **What does confirm step 7 was inert here, directly and from this very
     run**: the job made **44 meep domain builds**. The per-length path
     predicts exactly 1 (`M_in1`) + 1 (seed) + 3 lengths × 14 excited bands
     = 44; the merged path would have been 1 + 1 + 14 = 16. The production
     run therefore took the unchanged per-length path, exactly as step 7's
     `save_dir` fallback designs it to — a direct structural measurement,
     stronger than a numeric diff, and consistent with the bit-for-bit
     matches the two DC scripts (which *do* have valid baselines) already
     showed.
   - DC hybrid: job 11872184, `gts-sr108`, `-N 2 --ntasks-per-node=24`,
     `-t 900` (15h, ~44% margin over the measured 10:25:31 precedent).
     Batch: `batch/verification/
     run_step8_directional_coupler_hybrid_recalibration.batch`.

     **PASSED, 2026-08-15 — bit-for-bit identical.** 03:30:40 wall, exit 0
     (well under the 15h cap; faster than the 10:25:31 cache-cold precedent,
     consistent with some cache reuse — not chased further, since the
     comparison below is what actually matters). Compared directly against
     job 11655592_0's own printed output (`reports/11655592_0.out`, the
     most recent prior run pinned to the same exact 15.10 µm coupling
     length, pre-step-7):
     | quantity | job 11655592_0 (before) | job 11872184 (after) |
     |---|---|---|
     | MPB n_eff (even/odd) | 2.323677 / 2.300177 | 2.323677 / 2.300177 |
     | FDTD n_eff (even/odd) | 2.375225 / 2.354162 | 2.375225 / 2.354162 |
     | Δn_eff (FDTD) | 0.021063 | 0.021063 |
     | Through T / IL | 0.5075 / −2.95 dB | 0.5075 / −2.95 dB |
     | Cross T / IL | 0.4907 / −3.09 dB | 0.4907 / −3.09 dB |
     | Coupling ratio | 0.492 | 0.492 |

     Second confirmation of the same result as DC SOI160 above — expected,
     same reason: `save_dir` set, step 7's fallback took the identical
     per-length path, deterministic FDTD.
   - DC SOI160 hybrid: job 11872185, `gts-sr108`, `-N 1
     --ntasks-per-node=24`, `-t 1440` (24h — double the 720/12h that
     succeeded in `batch/dc_soi160_hybrid.batch`, but that completion
     (job 11716260) reused 5 cached entries from an earlier attempt;
     zero cache exists now, so this run is genuinely colder). Batch:
     `batch/verification/
     run_step8_directional_coupler_soi160_hybrid_recalibration.batch`.

     **PASSED, 2026-08-15 — bit-for-bit identical, not merely within
     noise.** 01:51:23 wall, exit 0 (well under the 24h cap). 10
     `Initializing structure` domain builds in this run's log vs. only 2 in
     job 11716260's (which reused 5 cached entries from an even earlier
     attempt) — real, direct evidence this run did more fresh computation,
     i.e. was colder, though not the fully-uncached case checked for
     before submitting (`example_sim_cache` inspection was of the wrong
     cache dir; this script's own cache may partially persist elsewhere —
     not chased further since it doesn't change the result below). Compared
     directly against job 11716260's own printed output
     (`reports/11716260.out`, the prior real run, pre-step-7):
     | quantity | job 11716260 (before) | job 11872185 (after) |
     |---|---|---|
     | MPB n_eff (even/odd) | 2.129795 / 2.101037 | 2.129795 / 2.101037 |
     | FDTD n_eff (even/odd) | 2.181661 / 2.156256 | 2.181661 / 2.156256 |
     | Δn_eff (FDTD) | 0.025405 | 0.025405 |
     | Through T / IL | 0.8989 / −0.46 dB | 0.8989 / −0.46 dB |
     | Cross T / IL | 0.0979 / −10.09 dB | 0.0979 / −10.09 dB |
     | Coupling ratio | 0.098 | 0.098 |

     Every printed digit matches — expected and correct, not a coincidence:
     since this script passes `save_dir`, step 7's fallback routes it
     through the exact same per-length code path as before the change, and
     FDTD is deterministic given identical inputs. Confirms step 7's wiring
     is a genuine no-op for this caller, exactly as designed.

   **Account note**: `gts-sr108-ece` hit `AssocGrpBillingMinutes` on all
   three submissions (confirmed broadly exhausted — even unrelated other
   users'/sessions' pending jobs were stuck on it at the same time, not
   specific to these jobs). DC/DC-SOI160 cleared under `gts-sr108`; the
   MMI 2×2 job's much larger request (`-N 2 × 24 ranks × 48h`) still hit
   the same gate there and only cleared under `gts-sr108-gedc_main` — the
   three accounts' available headroom shifts and needs checking per
   submission, not assumed from a prior session (see
   [[slurm-account-gts-sr108-ece]]).

   Results to be filled in here once each job completes.

## Open questions

- Does the merged approach's cache-key design (flagged above, not designed)
  want its own tier-style reuse (cf.
  `EMESolver.solve_section_length_sweep()`'s 3-tier cache,
  `resolved/calibrate_body_neff_mode_solve_skip_validation.md`'s extension
  session) — e.g. a length-independent fingerprint so a later call asking
  for a different/longer set of monitor positions on the same cross-section
  doesn't force a full re-run? Not designed, worth considering once step 2
  lands.
- Is 3 monitors (today's default) actually going to turn out to be enough?
  **Answered, 2026-08-14 (steps 5-6): yes for β, and the residual growth
  step 5 saw is now explained.** β at N=3 matches β at N=6 to `9.4e-05`
  rad/µm (non-monotonic fluctuation, not a trend converging toward more
  monitors) — no sign the shipped calibration numbers (DC's
  `job 11517245`, etc.) are under-sampled on the β axis. Step 5's residual
  growth turned out NOT to be a monitor-count effect (step 5's domain was
  fixed the whole time, corrected above) — step 6 isolated it as a genuine
  domain-length effect instead: the same 3 fit lengths' β/residual shift
  measurably (β range `5.78e-04` rad/µm) as the *shared domain around them*
  grows, independent of how many monitors are in it. Both effects stay 3+
  orders of magnitude under the guard thresholds throughout. Reduced scale
  only (`n_modes=2`, 1 wavelength); full production scale unchecked.
- Does the merged domain's size need its own minimum-clearance recommendation
  (beyond the existing `_check_aperture()` PML-crossing guard), given step 6
  found the domain-length effect is real but non-monotonic rather than
  simply "bigger is worse"? **Opened, 2026-08-14 (step 6).** Only 4 domain
  sizes tested (3×-8×spacing); not enough to fit a model or recommend a safe
  range, only to confirm the effect exists and stays guard-safe at every
  size tried. Worth a denser sweep if this ever needs a hard bound rather
  than an existence proof.
- Does the correlated-noise concern on the linearity/residual guard (flagged
  above) turn out to matter in practice, or is `decay_tol=1e-11` tight
  enough that it's inert? **Answered, 2026-08-13 (step 4): inert.** Old-vs-
  merged agreement (`max|beta diff| = 3.50e-05 rad/µm`) is 3-4 orders of
  magnitude tighter than the residual/margin guard thresholds — checked at
  reduced scale (`n_modes=2`, 1 wavelength), not yet at full
  `N_MODES=14`/11-wavelength production scale, so treat as strong evidence,
  not a production-scale guarantee.
- Does the merged path need its own guard wrapper mirroring
  `section_mode_beta()`'s pre-simulation branch-margin guard (Guard 2)?
  **Resolved, 2026-08-14 (step 7): no separate wrapper needed.**
  `solve_section_multi_length()` itself still has no guard of its own, but
  step 7 wired it in *behind* `section_mode_beta()`'s existing Guard 2 (and
  both post-fit guards), so calling through that entry point — which is how
  every real caller reaches it — gets identical guard coverage on either
  path. A caller that bypasses `section_mode_beta()` and calls
  `solve_section_multi_length()` directly still gets no guard at all, same
  as before; not a real gap since no production caller does that.
- Given the merged path won't actually activate for any current production
  caller until its cache-key design lands (all of them pass `save_dir`,
  step 7's own finding) — is designing that cache key worth doing before
  step 8's re-run, or should step 8 run once on the per-length path as
  today and revisit once caching is designed? Not decided; step 8 itself
  doesn't need an answer to proceed, since it's re-validating shipped
  numbers either way.
