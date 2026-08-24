# Merging `section_mode_beta()`'s per-length FDTD sims into one multi-monitor simulation — plan

**COMPLETE — the mechanism this document specs is built, validated, and wired in (steps 1-8, 2026-08-12 → 2026-08-16).**

**★ SUPERSEDED IN PART, 2026-08-17 — step 9 below. The merged path is now the
DEFAULT path for every caller, and the reason it is worth taking turned out not
to be the reason this document gives.** The cost model underneath steps 1-8
("the dominant cost doesn't scale with monitor count") was measured wrong: the
dominant cost is the modal decomposition, which scales *exactly* with monitor
count, and merging alone is worth ~1.09x at production scale, not the 2.49x-3x
this document projects. The fix that actually pays is `measure_diagonal_only`,
worth ~4x, and merging is worth a further ~1.9x on top of it. Both landed
2026-08-17. **If you are here to make another device's calibration cheaper, read
step 9 and nothing else — it is self-contained and tells you what to do.**
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

---

# Step 9 — the real cost model, and the standard way to calibrate (2026-08-17)

**Read this section on its own if you are here to fix another device's β
calibration.** It supersedes this document's cost reasoning; the mechanism and
convergence work in steps 1-8 still stand.

## What every calibration was actually spending its time on

Measured on job **11872191** (`mmi_2x2_hybrid.py`, cold cache, 48 ranks,
31.4 h), by tallying its own log rather than by reasoning:

| | amount | cost | share |
|---|---|---|---|
| MPB eigenmode solves | 7,249 @ ~13.9 s | **27.9 h** | **89 %** |
| FDTD time-stepping | 1,428,073 steps @ 0.0089 s | 3.5 h | 11 % |
| structure builds (`set_epsilon`) | 44 | 21 s | ~0 % |

The per-solve figure is confirmed independently at two rank counts — 13.8 s at
48 ranks (11872191) and 13.9 s at 24 ranks (job 11935291: 154 solves, **zero
timesteps**, 2,142 s). **The dominant cost does not scale with rank count.**
Adding nodes buys speed on the 11 % and nothing on the 89 %.

Stage breakdown of the same job, from cache-folder timestamps (the log prints
only one total):

| stage | wall | share |
|---|---|---|
| startup + geometry plots + input-arm FDTD (`M_in1`) | 44 min | 2.3 % |
| **β calibration** | **30 h 38 min** | **97.6 %** |
| EME body length sweep (+ its n_eff) | 23 s | 0.02 % |
| cascade + result plots | 13 s | 0.01 % |

The EME is free — `neff_override` skips the mode solve and propagation is
analytic. **Calibration is essentially the entire cost of a hybrid example.**

## Why: 42 simulations, each computing 14× what it needs

Two multipliers, and the second one was pure waste:

1. `EigenModeSource` launches **one band per simulation** (`eig_band=k`), so a
   14×14 transfer matrix costs 14 runs — × 3 fit lengths = **42**. The
   wavelength axis is *not* one of the multipliers: each run is one broadband
   Gaussian pulse with DFT monitors at all 11 wavelengths.
2. Each run then decomposed its monitor into **all** `n_measure` bands at
   **all** wavelengths — 14 × 11 = 154 MPB solves, ~36 min — and
   `section_mode_beta_details()` kept one of them, via `np.diagonal(T)`.

Split the log into per-`mp.Simulation` blocks and every calibration run is
**exactly 165 eigensolves**, with the same composition each time:

| what | band | solves |
|---|---|---|
| measurement monitor, the 13 bands the caller discards | ≠ k | 13 × 11 = **143** |
| measurement monitor, the one band the caller keeps | k | 11 |
| normalization monitor near the source | k | 11 |
| the `EigenModeSource` itself | k | **0** (no separate eigensolve logged) |
| | | **165** |

Whole-job accounting, which closes to the line:

    42 calibration sims × 165   = 6,930
    + MPB seed (14 bands × 11λ) =   154
    + M_in1 (14×11 + 11 norm)   =   165
                                  -----
                                  7,249   ← exactly the log's count

Per-solve cost follows from the same job two ways: (30.6 h − 3.4 h stepping) /
6,930 = **14.1 s**, and job 11935291's seed-only run 2,142 s / 154 = 13.9 s.
Take it as ~14 s.

## The fix: `measure_diagonal_only` (~4×, exact)

`FDTDSolver.solve_section(measure_diagonal_only=True)` decomposes band `k` only
in run `k`. **Exact, not an approximation** — the removed entries are the ones
the caller already discards. Off-diagonals return **NaN, not 0**, so a consumer
reading a never-measured entry fails loudly instead of taking "not computed"
for "no coupling".

Verified against job 11872191's cached full matrices, all three fit lengths:
diagonals **bit-identical** (`max|Δ| = 0.0e+00`), fitted phases bit-identical,
off-diagonals all NaN.

Two design points that are load-bearing:

- **The flag is written into the cache fingerprint ONLY when True.** Writing it
  unconditionally changes every config's hash and invalidated all 566 entries in
  `example_sim_cache/` at a stroke — caught during implementation, and now
  pinned by `tests/test_beta_calibration.py::TestMeasureDiagonalOnly::
  test_the_flag_is_omitted_from_the_fingerprint_when_false`.
- **A cached full matrix serves a diagonal request** by extracting its diagonal,
  so turning the flag on never discards a cache an earlier full run paid for.
  The reverse (a diagonal serving a full request) is refused, and mode-count
  superset reuse is disabled for diagonal runs.

## The merged path is now the default (a further ~1.9×)

`solve_section_multi_length()` has a cache key of its own as of 2026-08-17,
which was **the only thing** keeping it out of production. `merged_active` no
longer depends on `save_dir`/`cache_only`; the per-length path is still
reachable with `use_merged_domain=False`.

Post-fix, at 14 modes × 11 λ:

| | per sim | measure solves | normalization | stepping | total |
|---|---|---|---|---|---|
| 42 sims × 1 monitor | 11 + 11 = **22** | 42 × 11 = **462** | 42 × 11 = 462 | 3.4 h | **~7.0 h** |
| 14 sims × 3 monitors | 33 + 11 = **44** | 14 × 33 = **462** | 14 × 11 = 154 | 1.3 h | **~3.7 h** |

Note what merging does and does not do: **the measurement work is identical
either way** — 462 solves, because you decompose band `k` at three planes
regardless of whether those planes live in one domain or three. Merging removes
28 of the 42 *normalization* monitors (one per simulation, and there are 28
fewer simulations) and two-thirds of the decay waits. That is a real ~1.9×, but it is not the ~3× this document projected, and
step 4's measured 2.49× came from a configuration doing 36 MPB solves against
production's 7,249: **reducing `n_modes` 14→2 and wavelengths 11→1 removed
precisely the axis that dominates a real run.** Reduced-scale speedups from this
mechanism do not transfer; reduced-scale *equivalence* numbers still do.

Cumulative on the MMI 2×2: **30.6 h → ~3.7 h, ~8.3×.**

## What was considered and rejected

- **One simulation for all 14 modes.** Superposing sources makes the monitor
  read `c_i = Σ_j T_ij a_j` with no way to separate `T_ii`. Measured from the
  cached matrices, off-diagonal leakage reaches **5.6 %** of the diagonal
  (worst at modes 9 and 12), i.e. ~0.055 rad of phase error — over half the
  0.1 rad residual guard and ~10 % of the 1.55-4.2 % β excess being measured.
  It would also undo the diagonal fix, since with all bands present every band
  must be decomposed at every monitor.
- **More ranks.** See the two-rank-count measurement above.

## The MPB seed: cached, and solved at 3 wavelengths (2026-08-17, later same day)

Cutting the decode exposed a third cost that had always been there and had
always been invisible. **`section_mode_neff` — the MPB seed — is 154
eigensolves, ~31.5 min, and became ~39 % of the calibration** (job 11937894's
own breakdown at 4 of 14 modes). It was 1.9 % of the 30.6 h version. Nothing
about it changed; everything around it got ~10× cheaper.

The seed exists to do two things, and neither is produce the answer:

1. **Pick the 2π branch.** `arg(T_ii)` comes back wrapped into (−π, π], and with
   β ≈ 6-11.7 rad/µm over 1.086 µm spacings the true phase advances 7-13 rad
   between fit points. The wrapped values fit infinitely many lines, one per turn
   count; MPB's β — within 1.55-4.2 % of the answer — says which.
2. **Space the fit lengths**, via `beta_max` → `max_safe_spacing` (Guard 2).

Both are coarse. Branch selection needs the estimate within ~25 %
(`_BETA_MIN_BRANCH_FRACTION`), and production margins sit at 0.05-0.095 of the
half-distance to the neighbouring branch. Two changes followed:

### Cached

`section_mode_neff()` now takes `save_dir`/`cache_only`, looked up **before the
simulation is built**. Measured: cold 38.4 s → warm 0.001 s on a scaled-down
seed, identical `n_eff`. Mode-count superset reuse applies, legitimately here in
a way it does not for `solve_section()`: this method's monitor plane is pinned by
`monitor_center`/`monitor_size` and does not widen with the mode count, and band
*i* is MPB's *i*-th distinct band however many were requested — so a stored
14-mode seed serves a 6-mode request by slicing.

Its `_METHOD_ATTRS` allowlist is deliberately tighter than `solve_section()`'s:
a mode solve never time-steps and has no source, so `decay_tol`,
`source_margin`, `source_bandwidth_factor` and `monitor_margin` cannot reach it
and must not be fingerprinted, or two physically identical seeds miss each other.
`pad` **is** fingerprinted — it changes the solved domain without appearing in
any other kwarg. Pinned by `tests/test_beta_calibration.py::TestSeedCache`.

This does not help a first cold run. What it kills is the repeat: job 11935291
was 35:42 of which essentially all was the seed, and `cache_only=True` now
covers the seed too, so a plot/documentation re-run can no longer silently spend
half an hour there.

### Solved at 3 wavelengths and interpolated — `seed_anchors=3`, on by default

`n_eff(λ)` across a device band is smooth and very nearly quadratic, so the seed
does not need the caller's whole grid. `seed_anchors` (default **3**, active for
any caller with more than 3 wavelengths) solves `linspace(min, max, n)` and fits
each band with a degree-(n−1) polynomial on band-centred x — exact interpolation
through those points, not a least-squares compromise. `None` restores the old
behaviour; subsampling is skipped automatically for a single-point band or ≤ 3
wavelengths, so no degenerate fit is reachable.

**Only the seed is subsampled.** The FDTD time-stepping and the modal decode
always run on every wavelength passed in. Confirmed on the live call path: the
log reads `seed solved at 3 of 11 wavelengths and interpolated (degree 2); FDTD
and decode still use all 11`, and the run then demanded merged simulations on
all 11.

**MEASURED against the real 11-wavelength MPB seed — job 11940556, COMPLETED
46:31, exit 0** (`verify_seed_anchor_interpolation.py`, production scale: 14
modes, res 40, the same 2 µm seed window the calibration uses). The pre-
implementation budget had used this device's *calibrated* n_eff curve as a
stand-in, because the MPB seed had never been cached and so its own 11-point
curve did not exist anywhere. It is now cached, and these are the real numbers:

| | worst \|Δn_eff\| | as Δβ_seed | of branch half-spacing |
|---|---|---|---|
| modes 0-11 | ≤ 5.5e-06 | ≤ 2.2e-05 rad/µm | ≤ 7.7e-06 |
| mode 12 | 2.05e-05 | 8.3e-05 rad/µm | 2.9e-05 |
| mode 13 (nearest cutoff, most curved) | **1.53e-04** | 6.2e-04 rad/µm | **2.15e-04** |

Against production branch margins of 0.05-0.095, the interpolation consumes
~0.2 % of the available headroom at the worst mode. (The stand-in curve had
predicted 1.10e-04; the real seed is 1.53e-04 — same order, slightly worse.)

**The fit lengths do not move at all**, which is the result that matters for
every cache entry on disk:

    beta_max  truth 11.571417845   interpolated 11.571417845   d = 0.00e+00
    lengths   truth  [1.085983652, 2.171967305, 3.257950957]
              interp [1.085983652, 2.171967305, 3.257950957]

Exactly zero, for a structural reason rather than luck: `beta_max` is
`max(2πn/λ)`, which sits at mode 0 at λ_min — and λ_min is an anchor, so it is
solved, not interpolated. The anchors also reproduce the truth **exactly** where
they overlap (`max|Δ| = 0.00e+00`): the anchored solve is the same solve, not an
approximation of one. So enabling this default invalidates nothing — not the
merged entries, not the 30.6 h of per-length entries from job 11872191.

Cost: **154 eigensolves → 42, 36.3 min → 9.8 min (3.7×)**, and zero on any
re-run now that the seed is cached.

### Convergence in `seed_anchors` — measured, free

On-grid subsets of the same 11-point truth, so no extra solves:

| anchors | worst \|Δn_eff\| | of branch half-spacing |
|---|---|---|
| 3 (the default) | 1.53e-04 | 2.15e-04 |
| 5 | 3.01e-06 | 4.23e-06 |
| 7 | 3.60e-08 | 5.05e-08 |

Roughly two orders per two anchors. `seed_anchors=5` costs ~3 extra minutes if a
device's band ever needs it. **This is the answer to "was the default
convergence-tested"** — yes, at production scale, on this device.

One caveat stands: **`seed_excess` (β/β_seed − 1) is now computed against an
interpolated β_seed**, and that quantity is *the* reported measurement of this
calibration, not an internal. Measured range on this device 1.535 %-4.227 %,
worst shift from interpolating **0.0104 percentage points** (mode 13). The
docstring says so rather than leaving it implied.

### Raise `seed_anchors` when

a band is wide, `n_eff(λ)` is not well described by a degree-(n−1) polynomial
across it, or a mode crossing sits inside the band. Cost is linear in the
number. `tests/.../test_more_anchors_help_the_curved_mode` asserts the knob
works in the direction you would want.

## Still on the table: the modal basis is re-derived ~47 times

The eigenbasis depends only on (cross-section, band, frequency), and a uniform
body section presents the same cross-section to every run. Job 11872191's log:
**7,249 mode solves producing 165 distinct results** — 38 minutes of distinct
work stretched into 27.9 h, because each `mp.Simulation` starts fresh and
`get_eigenmode_coefficients()` has no cross-simulation basis cache.
`section_mode_neff()` computes that exact basis for the seed and discards it.

Post-fix, merged + diagonal still executes 616 solves against those 165 — a
remaining ~3.7×. Reaching it means computing the overlap integrals by hand from
`get_dft_array()` instead of calling `get_eigenmode_coefficients()`, which means
re-deriving meep's normalization convention. **Treat that as high-risk**: this
repo already carries `_EXCITE_NORM_MISMATCH_TOL` from being bitten by exactly
that (transfer-matrix columns 10-20× too large, `fdtd.py:51`). Not attempted.

## How to apply this to another device

`mmi_2x2_hybrid.py` needed **no change** — both fixes are in the library, so any
script calling `section_mode_beta()` / `section_mode_beta_details()` /
`HybridSolver(calibrate_body_neff=True)` gets them automatically. For another
device:

1. **Change nothing in the example.** Just re-run it.
2. **Budget the re-run from the cost model**, not from its last wall time.
   Three terms, all measured on the MMI 2×2:

       merged calibration   n_modes x (n_lengths + 1) x n_wl x 14 s
                            one decode per fit-length monitor plus one for
                            normalization, excited band only
       MPB seed             n_modes x seed_anchors x 14 s   (0 if cached)
                            seed_anchors is 3 by default, NOT n_wl
       time-stepping        n_modes x timesteps x 0.0089 s

   Rank count does not help the first two — they are MPB, measured at the same
   ~14 s per solve at both 24 and 48 ranks.
3. **Expect the calibration to re-simulate once.** Its per-length entries do not
   serve the merged path, deliberately — see step 6's domain-length sensitivity.
   Everything else in the script (input-section FDTD, EME sweeps) still hits
   cache. Devices with a calibration in `example_sim_cache/` today:
   `directional_coupler_hybrid.py` (last cold: 3 h 31 m, job 11872184),
   `directional_coupler_soi160_hybrid.py` (1 h 51 m, job 11872185).
4. **Verify against the previous β, don't assume.** Two things changed under
   the calibration, and both are guard-safe but real:
   - the merged domain is a different simulation (step 6: `5.78e-04` rad/µm of
     domain-length sensitivity, reduced scale, **production scale still
     unchecked**);
   - `β_seed` is interpolated from `seed_anchors` wavelengths, which shifts the
     reported `seed_excess` by ~0.007 pp at the worst mode.

   Compare the new `calibration/fits/calibration_summary.png` against the old
   one, and check the residual and branch-margin grids against their 0.1 rad and
   0.9 guards. **If the branch margin moved materially, raise `seed_anchors`
   before anything else** — that is the one failure mode the interpolation can
   introduce, and it is visible in exactly that panel.
5. **`cache_only=True` now covers the seed as well as the simulations**, so a
   plot-only re-run is genuinely free and can no longer surprise you with half
   an hour of MPB.
6. **Do not reduce `n_modes` or `n_wl` to make a convergence study affordable**
   without saying so — that is the exact move that produced this document's
   wrong cost model.

## Production result — job 11937894, COMPLETED 04:14:12, exit 0

The first real run of the merged path at production scale (`mmi_2x2_hybrid.py`,
14 modes, 11 wavelengths, 48 ranks). Both fixes confirmed by the log's own
arithmetic, not by inspection:

    FDTD sims  :  14   (42 -> 14: merged domains)
    eigensolves: 770   = 154 seed + 14 x 44   (165 -> 44 per sim: diagonal decode)

**Cost 30.6 h -> 4.14 h, 7.4x.** Short of the ~3.7 h projection only because
this job predates the seed cache and the 3-anchor default by a few hours and so
still paid a full 36-minute 154-solve seed; the next run of this script pays
neither and should land near 3.4 h.

### The device answer moves — measured, and it is step 6's effect

Against the per-length baseline from job 11872191 (same script, same cache, only
the calibration path differing), over the full 201-length x 11-wavelength grid:

| | max \|Δ\| | mean \|Δ\| |
|---|---|---|
| `T_bar` | **5.5e-03** | 1.7e-03 |
| `T_cross` | **3.6e-03** | 1.5e-03 |

At the nominal 28.0 µm, λ=1.5475 µm: bar 0.456172 -> 0.455807, cross 0.518707 ->
0.518836.

**This is not noise and not a regression.** It is step 6's domain-length
sensitivity (`5.78e-04` rad/µm of β, reduced scale, non-monotonic) arriving at
production scale for the first time, propagated through a 27.5 µm body and the
full cascade. Three snug domains and one shared long one are different
simulations; that was the stated reason per-length cache entries deliberately do
not serve the merged path, and this is the size of the thing that reasoning was
protecting against.

Scale: ~0.1 % of the split ratio, an order of magnitude inside this device's own
gap-5 discrepancy (0.76 %/0.59 % max|ΔT| vs the full-device FDTD reference), and
both post-fit guards stayed clean throughout.

**So the correct claim is narrower than "the two paths agree":** they agree
exactly on the fit lengths (β_max sits on an anchor, verified `d = 0.00e+00` by
job 11940556), and they differ at the 1e-3 level on the device answer. Any
device re-run onto this path should expect a shift of that order and check it
against its own error budget rather than assuming bit-identity. Production-scale
merged-vs-per-length β equivalence — the per-mode number behind this
transmission shift — is still unmeasured; see the open questions.

## Rollout, 2026-08-17 — every calibrating example moved to the new path

`mmi_2x2_hybrid.py` first (job 11937894), then the rest. **Only one script
needed an edit**, which is the point: all three changes live in the library, so
any script reaching the calibration through `section_mode_beta()`,
`section_mode_beta_details()` or `HybridSolver(calibrate_body_neff=True)` picks
them up by re-running.

| example | script change | job |
|---|---|---|
| `mmi_2x2_hybrid` | none | 11937894 |
| `directional_coupler_hybrid` | none — already `CALIBRATE_BODY = True` | 11942633_0 |
| `directional_coupler_soi160_hybrid` | none — already `CALIBRATE_BODY = True` | 11942633_1 |
| `mmi_1x2_hybrid` | **default flipped off → on** | 11942633_2 |
| `directional_coupler_soi160_hybrid_length_sweep` | none | 11942662 |
| `hcg_fabry_perot_hybrid` | none — calibrates nothing itself | 11942662 |

Two ordering facts that matter:

- **`hcg_fabry_perot_hybrid` is downstream of the SOI160 coupler.** It is a
  `beef.network.Circuit` assembly that consumes that coupler's saved
  `smatrix_hybrid.pkl` rather than solving anything itself, so it must run
  *after* 11942633_1 or it republishes a cavity spectrum built on the old
  coupler. Submitted with `--dependency=afterok:11942633_1`. The HCG mirror half
  is unaffected — no body section, nothing to calibrate.
- **Each calibration re-simulates once**, deliberately: per-length entries do
  not serve the merged path (step 6's domain-length sensitivity is real).
  Everything else in every script still hits cache.

### The one script edit: `mmi_1x2_hybrid.py`

`CALIBRATE_BODY` defaulted **off** (`MMI1X2_CALIBRATE_BODY=0`) on measured
grounds that still stand: job 11608690 (2026-08-01) found calibration takes this
device's gap from 9.11 % to 7.30 % — ~20 % relative against the coupler's 98 % —
with the sign structure untouched. **That has not changed, and gap 5 stays open
on exactly those grounds.** What changed is the tradeoff: a calibration used to
be ~97.6 % of a hybrid example's wall time, so "off unless it earns its keep"
described a real cost. At ~3.7 h on the 2×2 it no longer does, and propagating
the body on the FDTD's own β is the self-consistent choice regardless of whether
it closes a given device's residual error — `directional_coupler_validation.md`
calls it "a correctness fix rather than a fudge" for that reason.

The escape hatch is still one env var (`MMI1X2_CALIBRATE_BODY=0`), the output
suffix inverted (`_uncalibrated` now, `_calibrated` before), and gap 5's
baseline was copied to `example_output/mmi_1x2/mmi_1x2_hybrid_uncalibrated/`
before the flip. **`tests/golden/mmi_1x2_hybrid.json` is unaffected**: that
target is an independent reduced-scale re-implementation in
`tests/test_golden_examples.py` (res 20, `n_modes=6`, no `calibrate_body_neff`),
not a run of this script — checked, not assumed.

### Notebooks: a checker, and the drift it found

Notebooks here are hand-mirrored copies of the scripts (`notes.md` Solvers item
6), and its predicted drift had already happened once by audit. The calibration
work would have made it worse — three notebooks were still calling the old
`section_mode_beta()` and had never received the `plot_beta_calibration` block.

`scripts/check_notebook_sync.py` is the cheap half of that item's
recommendation: it diffs a named set of constants and API calls between every
`.ipynb` and its paired `.py`, normalises the deliberate differences (a script's
`os.environ.get(NAME, "11")` default equals a notebook's literal `11`) and
exits 1 on drift, so it can gate CI. It does not generate notebooks from
scripts; that remains the better fix and remains unbuilt.

Run on the tree it found **six** out of sync — three from this session's work
and three pre-existing, both of the latter being exactly the failure mode item 6
warned about:

| notebook | drift | age |
|---|---|---|
| `directional_coupler_hybrid` | missing `plot_beta_calibration` / `section_mode_beta_details` | this session |
| `directional_coupler_soi160_hybrid` | same | this session |
| `mmi_2x2_hybrid` | same | this session |
| `mmi_1x2_hybrid` | `CALIBRATE_BODY` default | this session |
| **`mmi_2x2_fdtd`** | **`N_WL` 11 vs the script's 51** | since 2026-08-10 (job 11791990_1) |
| **`modal_multiplexer_hybrid`** | **`N_MODES` 16 vs the script's 15** | pre-existing |

All six fixed; the checker now exits 0 on 18 pairs. Every notebook re-validated
as parseable JSON with every code cell compiling — worth doing, because the
first attempt at injecting the new import produced `from beef import (,
plot_beta_calibration` in two of them and only that check caught it.

### The geometry plot now draws the domain that ran (2026-08-17, later)

The rollout above left one thing stale, found by looking at job 11937894's own
output rather than at the code: `calibration/geometry/` still held **three**
PNGs, one per fit length, each drawing a snug window with its own private
source. The merged run is one domain with three monitors. The plot was
merged-*aware* — `plot_beta_calibration_geometry()` appends "(merged run —
per-length equivalent shown)" to each title when `cal.merged` — so it was never
false, but the module's own stated purpose is "what was actually simulated," and
under the new default it no longer showed that.

`plot_beta_calibration_merged_geometry()` draws the real thing: one domain at
`max(lengths)`, one source, one interface monitor per length.
`plot_beta_calibration()` now picks off `cal.merged`, so **no example script
changed** — same rule as the rest of step 9. Per-length windows still render for
a per-length calibration, and `per_length=True` asks for both on a merged one
(useful against a pre-merge run's figures); asking for a merged figure from a
per-length calibration raises rather than inventing a domain that never existed.

The monitor positions are **not** re-derived in the plot. The formula moved out
of `solve_section_multi_length()` into `merged_mea_positions()`, which both now
call — the same anti-duplication rule `section_fdtd_layout()` and
`BetaCalibration` each already record, and the third time this file has had to
apply it.

Verified against job 11937894's cached calibration (`cache_only=True`, no
simulation): `merged=True`, lengths 1.086/2.172/3.258, guards clean. The drawn
domain is `[-1.629, 1.629]`, source at local `x=3.258` (`source_margin=0`), norm
monitor at `3.008`, interface monitors at `2.172 / 1.086 / 0.000` — each exactly
`L_i - source_margin` from the source, i.e. its per-length twin's distance.
Every affected script passes the same `monitor_center`/`monitor_size` to the
calibration and to the plot (`_iface_mon` in the three couplers, nothing in
either for the 1×2), so the drawn window is the simulated one in all of them.
Landed before 11942633 started, so all four pending devices produce it on this
run.

**What the picture is for:** the shortest length's monitor sits 1.086 µm from
the source in the 2×2 — interior, with two more monitors and the source's own
near field between it and the domain edge. That clearance is the open question
below, and it is not visible at all in a plot that redraws each length with its
own source.

## Open questions this step opens or changes

- **Production-scale merged-vs-per-length equivalence is still unmeasured.**
  A script exists for it —
  `examples/verification/python/equivalence_test_solve_section_multi_length_production.py`
  — comparing merged against the cached per-length β at 14 modes × 11 λ, with
  the per-length half served free from cache. Submitted as job 11935871 and
  **cancelled**: it was budgeted at 20 h from the old 3× cost model when the
  real figure is ~28 h. Re-budget before resubmitting.
- **Does the merged path need a minimum source-to-monitor clearance rule?**
  Unchanged from step 6 — still open, still only 4 domain sizes tested. Now at
  least *visible* per device, in the merged geometry plot above (1.086 µm on
  the 2×2); a plot is not a threshold, and none is asserted.
- **Is the basis-reuse 4× worth the normalization risk?** New, undecided.

---

# Step 10 — the audit that closes the rollout (2026-08-18)

Step 9's rollout table was written as "every calibrating example moved to the
new path." A direct audit of that claim, prompted by the user asking whether the
2×2, the 1×2, the coupler, the HCG Fabry-Perot and the **ring resonator** were
all on the current calibration with a run behind them, found it true for the
five it listed and **incomplete in two ways it could not see**: one device was
never in the table at all, and one of the devices that was in it had run with
its calibration diagnostics silently disabled.

## Every shipped example, re-run on one library state — jobs 11981703, 11982196

The five step-9 devices were re-run together rather than trusted, because
`mmi_2x2_hybrid`'s result (job 11937894, finished 16:19) predated the last three
library edits of 2026-08-17 (16:48–16:50: `merged_mea_positions()` extraction,
`plot_beta_calibration_merged_geometry()`, and the seed cache / `seed_anchors=3`
default). Its calibration *figures* had been regenerated afterwards by a
`cache_only` pass, but its device answer had not been produced under the current
library, and "nothing should move" was an argument, not a measurement.

Nothing moved. Every device came back bit-identical to its step-9 run, in
**under 90 seconds each** — the merged calibration entries all hit cache, which
is itself the evidence that none of the three late edits perturbed a cache key:

| example | job | wall | result vs step 9 |
|---|---|---|---|
| `mmi_1x2_hybrid` | 11981703_0 | 0:53 | T 0.3749 / 0.3749 / 0.7498 — identical |
| `mmi_2x2_hybrid` | 11981703_1 | 0:56 | bar 0.4558, cross 0.5188, total 0.9746 — identical |
| `directional_coupler_hybrid` | 11981703_2 | 1:21 | L₅₀ 15.080479 µm, T 0.5071 / 0.4911 — identical |
| `directional_coupler_soi160_hybrid` | 11981703_3 | 0:46 | T 0.8989 / 0.0979 — identical |
| `..._soi160_hybrid_length_sweep` | 11981703_4 | 0:45 | L₅₀ 11.9016 µm — identical |
| `hcg_fabry_perot_hybrid` | 11982196 | 7:50 | FSR 2.255 nm — identical |

So the seed cache and `seed_anchors=3` really are inert on an already-calibrated
device, at production scale, on four devices rather than the one job 11940556
measured. `hcg_fabry_perot_hybrid` ran after the coupler again, for the reason
step 9 gives.

Two stale artifacts cleaned while there: both couplers'
`calibration/geometry/` still held the three `length_*.png` windows from their
last *pre-merged* run (2026-08-16) sitting beside the current
`merged_domain.png`. Neither was regenerated by the merged path and neither
describes the simulation that ran. Deleted; the length sweep's folder, which
never had a pre-merge run, already showed the correct clean state and is the
reference for what these should look like.

## Finding 1 — `mmi_1x2_hybrid` was calibrating and reporting nothing

Job 11942633_2 is recorded in step 9's table as the 1×2's rollout run, and it
was: `calibrate_body_neff = True`, printed by the script, body propagated on the
FDTD's β. But its output directory has **no `calibration_section1/`**, and its
log never printed the "Calibration plots saved" line. A calibrated run with no
calibration diagnostics, exit 0, no warning.

The cause is an interaction between two things that are individually correct.
`HybridSolver._solve_sections()` cleared `self.last_calibrations` at the top of
every pass — deliberately, so a stale entry from an earlier pass can never be
plotted as this one's. And `mmi_1x2_hybrid.py` calls `solve()` **once per output
port with a shared `cache=` dict**, deliberately, so the input and body sections
are simulated once rather than twice. On the second call the body section hit
that cache and `continue`d before `_calibrated_neff()` ran, so the dict the
script then loops over was empty and its `plot_beta_calibration()` block did
nothing at all.

**Fixed in `beef/solvers/hybrid.py`:** a section served from the in-memory cache
now carries its previous entry forward. That is not a weakening of the
no-stale-entries rule — the cache key includes `calibrate_body_neff` along with
everything else that changes the matrix, so a hit means this pass is *reusing
the matrix that very calibration produced*. The entry describes what is being
used, which is exactly the contract.

Pinned by `tests/test_beta_calibration.py::
TestLastCalibrationsSurvivesTheInMemoryCache`, four tests: the carry-forward
itself, that a genuinely re-solving pass still reports its own measurement and
not the previous one, that the no-cache path is unchanged, and that turning the
flag off still empties the dict. Verified to fail without the fix
(`assert set() == {1}`), not just to pass with it.

Re-run as job 11981703_0: `calibration_section1/` now holds 8 per-mode fits, the
summary, and `merged_domain.png`, and the device answer is unchanged — which is
the expected shape, since only the *reporting* was broken, never the β.

**Scope: this affects any script that calls `solve()` more than once with a
shared cache and then reads `last_calibrations`.** Two scripts in this repo do,
not one:

- `mmi_1x2_hybrid.py` — two `solve()` calls, one per output port. Fixed and
  re-run here.
- **`mmi_1x4_hybrid.py` — four `solve()` calls sharing one `section_cache`
  (line 333), `CALIBRATE_BODY = True` since 2026-08-13, and the same
  `last_calibrations` plotting loop copied from the 1×2 (line 482).** It has the
  bug identically, and its shipped output directory has no `calibration_section*`
  folder either. The library fix repairs it, but **this device is not re-run
  here**: it is mid-campaign (`mmi_1x4_sin300_validation.md`, prewarm jobs still
  running), and starting a run against that would collide with work in flight.
  Its next run picks the diagnostics up for free — no script edit needed. Flagged
  rather than done.

The 2×2 and both couplers call `section_mode_beta_details()` directly and were
never affected. Fixing this in the library rather than in the script is what
makes the 1×4 a free repair instead of a second edit — and what keeps it from
being a trap for the next multi-port device.

## Finding 2 — the ring resonator was never audited for this at all

**`ring_resonator_hybrid.py` and `ring_resonator_hybrid_cascaded.py` do not
appear in step 9's rollout table, and the reason is that they had never been
considered for calibration at any point.** `ring_resonator_validation.md` — a
3,000-line document — does not contain the strings `calibrate_body_neff` or
`section_mode_beta` anywhere.

The reason is structural, not an oversight in the rollout: **neither ring script
goes through `HybridSolver`.** They are raw meep/meow, because `FDTDSolver` only
builds x-normal sources and monitors and the two ring-side ports are y-normal
(each script's module docstring says so). So `calibrate_body_neff=True` had
nothing to attach to, and a rollout that swept "every script reaching the
calibration through `section_mode_beta*()` or `HybridSolver`" swept past the one
device that reached it through neither.

### It applies, and on the cascaded route it applies harder than on a coupler

Three consumers of an MPB index, all of them inside a loop whose every other
term is a meep measurement:

1. **`phi_bend_correction`, both scripts.** `run_bend_fdtd()`'s coefficient is
   measured between monitor planes that sit `BEND_D_IN`/`BEND_D_OUT` of straight
   waveguide from the arc ends, and the assembly de-embeds that straight with
   `exp(-1j·2π·n_eff·(D_IN+D_OUT)/λ)`. The phase being removed was put there by
   meep's time-stepper. Removing it with MPB's index subtracts the wrong number
   from a measured quantity — **there is no self-consistency argument on the
   other side of this one**; it is simply the wrong β.
2. **`phi = 2π·n_eff·L_outer/λ`, both scripts.** The analytic straights. Same
   argument as an EME body spliced between FDTD sections: a round trip assembled
   from `b`, `t`, `k1`, `k2`, `t_ring` (all meep) plus `phi` (MPB) changes answer
   depending on which engine carried which segment.
3. **`P_body`, cascaded only — and this is the textbook `calibrate_body_neff`
   case, not an analogy to it.** `S_coupler = M_left.T @ P_body @ M_left`, with
   `M_left` measured by meep at the interface plane and
   `P_body = diag(exp(1j·φ_even), exp(1j·φ_odd))` propagated over `BODY_LENGTH`
   on meow's supermode indices. That is an EME body between two FDTD sections,
   with item 6's amplification intact: both supermodes biased the same way, the
   differential tiny, and coupling riding the `β/Δβ ≈ 99` lever on it.

Scale, on the monolithic route's own numbers: the res-40 MPB-vs-meep gap
measured on the coupler is ~2.2–2.4 % in β; on `L_OUTER = 7.00 µm` at 1.55 µm
that is ~1.5 rad of round-trip phase, **roughly a quarter of an FSR**. That is
the same order as this device's gap-3 mode-resolution finding, which its own
validation document calls "not a caveat, it is a wrong answer for the one
quantity this device exists to report." The cascaded route's analytic arm is
2.00 µm rather than 7.00 µm (the DC-style bend accounting moved 5.00 µm into
`b`), so term 2 is ~3.5× smaller there — but term 3 is unique to it.

### What was built

`UniformCrossSection` — which exists for exactly this and whose module docstring
names **both** of these cross-sections as its motivating examples, verbatim —
plus `FDTDSolver.section_mode_beta_details()`. Nothing was re-derived: going
through the library inherits the merged multi-monitor domain, the diagonal-only
decode, the cached 3-anchor seed, the 2π branch guard, the residual guard and
`plot_beta_calibration()`'s figures. Hand-rolling a β fit against this script's
own conventions would have duplicated tested code and produced a β with none of
the guards attached.

The solvers are configured from **each script's own module-level constants** —
resolution, PML, source/monitor margins, decay tolerance, source bandwidth,
2-D/3-D, and the resolved dispersion flag — because the β being measured is only
the β the ring's bend and coupler runs accumulate if the simulation around it
matches. `sim_width` is `MESH_WIDTH` for the strip and `BODY_MESH_WIDTH` for the
pair, matching each meow counterpart's transverse window, so the excess between
the two engines is a statement about the engines and not about two domain sizes.

`CALIBRATE_SPAN = 6.0 µm` is a budget, not a device dimension: `spacing =
min((span − 2·dpml)/3, max_safe_spacing/2)` and at β ≈ 9.5 rad/µm the physics cap
is 1.32 µm, needing span ≥ 5.0 for the span not to bind. It is deliberately
**not** `BODY_LENGTH` — tying fit lengths to the body's real length is the
mistake `section_mode_beta()`'s own docstring warns about.

Both scripts keep the meow/MPB curves and report the calibrated ones against
them (per-mode excess, and for the pair the Δn_eff excess, which is the number
item 6 is actually about). Both write `calibrated` and the uncalibrated curves
into `ring_data.npz`, so the comparison never needs a second run to reconstruct.
`RING_CALIBRATE_BODY=0` restores the old behaviour and inverts the output suffix
to `_uncalibrated`, the same escape hatch `mmi_1x2_hybrid.py` uses.

**Skipped in 2-D, loudly.** A `RING_THREE_D=0` run has no (y, z) cross-section
for `UniformCrossSection` to present, so the calibration is skipped with a
printed notice rather than quietly measuring a different structure.

Notebooks: both ring notebooks were brought back into sync and
`scripts/check_notebook_sync.py` is at **0 out of sync across 18 pairs**, with
every code cell in both re-verified as compiling — worth doing for the reason
step 9 records, since a previous hand-edit of notebook JSON produced
`from beef import (, plot_beta_calibration`.

### Status: submitted, NOT yet landed

`batch/verification/run_ring_calibrated.batch` — four runs sequentially in one
job (monolithic calibrated, cascaded calibrated, then both uncalibrated
baselines). Sequential deliberately: the two scripts request an **identical**
bare-strip calibration and so share a cache key, which concurrent tasks would
race to write.

Two things were checked before submitting, both cheaply and both worth the
minutes:

- **The calibration path itself runs on a `UniformCrossSection`**, measured
  rather than assumed — a reduced-scale 3-D probe (res 20, 3 λ, 22 min). It
  came back `merged=True`, fit lengths `[1.3102, 2.6204, 3.9307]` µm — i.e.
  capped by `max_safe_spacing/2` as intended and **not** by `CALIBRATE_SPAN`,
  which is the check that 6.0 µm is large enough — `n_eff` 2.4287 / 2.3635 /
  2.2962 across 1.50–1.60 µm, worst residual **3.62e-03 rad** against the 0.1
  guard, worst branch margin **0.243** against 0.9, turn ambiguity 0.0017, and
  no wrap count moved off the seed's own rounding. Production res 40 should
  tighten all of these.
- **`plot_beta_calibration()` renders both cross-sections.** This runs *after*
  the expensive part in both scripts, so a raise there would waste the
  simulation it had just paid for. Exercised with a synthetic
  `BetaCalibration` of the real shape — no meep, seconds — for the strip and
  the pair, at their respective `sim_width`s.

That second check paid for itself immediately, on something no assertion would
have caught: the drawn domain showed **the coupled pair sitting 0.35 µm
off-axis**. The first version handed `UniformCrossSection` the ring's own
`(0, Y_RING)` centres, on the reasoning that the monitors span the full
`sim_width` about y=0 so shifting would move the pair off its window. That is
backwards — FDTDSolver's domain *is* centred on y=0, so not shifting is what
puts the pair off-centre, in the simulation and in the monitor that decomposes
it. `verify_ring_eme_engine_agreement.py` had already settled this convention
and states it explicitly: the ring centres this cross-section's oxide on
`BODY_Y_CENTER = Y_RING/2` while every `DeviceGeometry.section_gds()` centres
its oxide on y=0, so the strips go in at `(-Y_RING/2, +Y_RING/2)` and the
comparison is between the same structure rather than two different ones.
Corrected to match; clearance goes from 1.43/2.13 µm to 2.13/2.13. The first
submission (job 11982549) was **cancelled 2 min in** to make the change — it
was still on the strip calibration, which this does not touch. Resubmitted as
**job 11983058**.

### The cost model was wrong, and the reason is a finding of its own

The batch script budgeted this from step 9's calibration cost model "not from
these scripts' last wall times (their FDTD halves are all cached; only the
calibrations are new)". **The second half of that was false**, and job 11983058
showed it within the hour: `run_coupler_fdtd('bus')` took **~55 min to
simulate** rather than loading.

Checked against the cache rather than guessed at — every stored entry's
`config.json` compared field-by-field against the one this run wrote:

| family | entries usable at 3-D / res 40 / 101 λ |
|---|---|
| `RH_run_coupler_fdtd_bus` | 1 of 8 — and that one is this run's own, written 11:56 |
| `RH_run_coupler_fdtd_ring` | **0 of 7** |
| `RH_run_bend_fdtd` | **0 of 7** |
| `RH_cascaded_run_left_fdtd_bus` | 1 of 16 |
| `RH_cascaded_run_left_fdtd_ring` | 1 of 16 |
| `RH_cascaded_run_bend_fdtd` | 1 of 16 |

The older monolithic entries miss on `three_d` (eight of them are 2-D runs) and,
for the one 3-D entry, on `resolved_dispersive`/`resolved_si_reparam`/
`resolved_sio2_reparam` being `None` where today's run resolves them — i.e. it
predates the 2026-08-03 **M3/M4 material fingerprint**. So:

**`ring_resonator_hybrid.py`'s FDTD half has been stranded since 2026-08-03 and
this is its first 3-D production run under the current material model.** The
ring's own validation document predicted exactly this ("⚠️ This strands the 39
mode-solve cache entries… The 6 FDTD entries… remain stranded by the earlier axis
fix, and are now *additionally* stranded by the M3 material fingerprint") — the
prediction simply was never cashed in, because the device was not re-run again.
**The cascaded route was NOT cached either — retracting a claim made one
paragraph too early.** The table above was built from a filter over four fields
(`three_d`, `resolution`, a non-null `resolved_dispersive`, 101 wavelengths),
which is *necessary but not sufficient*: it reported 1 usable entry per cascaded
family, and stage 2 then re-simulated all three anyway, with **zero cache hits**.
A full field-by-field diff — which is what should have been done from the start,
and was done for the monolithic — finds the real reason, and it is a different
one:

    closest pre-existing entry (2026-08-03 18:30), differs in 2 fields:
      L_C                : now None   old 4.0
      geometry_polygons  : now x in [-4.5, 4.5]   old x in [-9.5, -0.5]

`run_left_fdtd()` acquired a **translation-invariance normalisation** after that
run: it expresses its polygons in a canonical centred frame and drops `L_C` from
the key, precisely so that editing `L_C` no longer invalidates an entry whose
physics does not depend on it. `ring_resonator_hybrid.py`'s own `_sim_params()`
docstring refers to it ("no translation-invariance trick available the way there
is in `ring_resonator_hybrid_cascaded.py`'s `run_left_fdtd()`"). Like any
re-keying it turns its old entries over exactly once, and this is that once.

So **all six ring FDTD families were cold, for two unrelated reasons**: the
monolithic on the M3/M4 material fingerprint plus `three_d`, the cascaded on the
`run_left_fdtd` re-keying. The methodological lesson is the one this document
keeps re-learning about duplicated formulas: **diff the whole config, never a
hand-picked subset of its fields.**

Revised estimate, from measured wall times rather than from a model: monolithic
came in at **1 h 55 m** (bus + ring-excited + bend + strip calibration);
cascaded needs its own three FDTD domains too (see the retraction above) plus
the coupled-pair calibration, running ~**1 h 30 m**; the two uncalibrated stages
are a few minutes each once everything above is cached. **≈ 3.5 h against an
8 h wall** — the total happens to be what was projected, for different reasons
than projected. A timeout would be recoverable rather
than wasteful — everything computed is cached, so a resubmit resumes.

Worth stating plainly: **the ring's comb will move for two reasons at once on
this run**, the calibration and the first post-M3/M4 3-D FDTD, and the
uncalibrated baseline in the same job separates them — it is being produced under
the same material model, so the calibrated-vs-uncalibrated difference isolates
the calibration alone.

### Stage 1 landed — and it reproduces item 6's number on a different stack

`ring_resonator_hybrid.py`, calibrated (job 11983058):

    Bare-strip n_eff at lambda = 1.550 um:
        MPB mode index  2.339523
        FDTD own phase  2.363938   <- the straights propagate on this
        excess +1.044%  (band range +0.998% .. +1.044%)
        -> round-trip phase shift +0.6928 rad = +11.0% of a full FSR
    Measured mean FSR: 12.560 nm (predicted 12.87)

Guards clean with room to spare, read off `calibration/fits/
calibration_summary.png`: fit residual **0.0018–0.0026 rad** against the 0.1
guard (~40×), branch margin **0.084–0.089** against 0.9 — and that margin range
sits exactly where step 9 says production margins sit (0.05–0.095). Fit lengths
`[1.264, 2.528, 3.792]` µm inside the 6.0 µm span, so `CALIBRATE_SPAN` was not
the binding constraint, as intended.

**The +1.044% is smaller than the ~2.2–2.4 % this document's estimate assumed,
and the reason is not a discrepancy — it is gap 3.** The script reports the
excess against the ring's own **res-100** meow curve, because that is what it
would otherwise have propagated on. The calibration's own `seed_excess` panel,
which measures against `FDTDSolver`'s res-40 MPB seed, reads **2.23 % → 2.37 %
across the band** — the like-for-like comparison. By hand against gap 3's own
res-40 figure for this recipe (2.311706): `2.363938 / 2.311706 − 1 = +2.259 %`.

Set beside `directional_coupler_validation.md` item 6, measured on **SOI_160**,
different waveguide dimensions, different device:

| | meep(res 40) vs MPB(res 40) |
|---|---|
| coupler, even supermode | +2.233 % |
| **ring, bare strip** | **+2.259 %** |
| coupler, odd supermode | +2.361 % |

The ring lands between the coupler's two modes. **This is the first independent
reproduction of item 6's discretization gap on another stack**, and it was not
sought — it fell out of a device audit. It also settles what the two fixes are
each doing: gap 3's move to `MODE_RESOLUTION = 100` already closed a bit over
half the res-40 gap, and the residual +1.044 % is the part **mesh refinement of
MPB cannot reach**, which is exactly item 6's claim about why this is a
correctness fix and not a convergence problem.

## LANDED — job 11983058 COMPLETED, exit 0, all four stages, ~3 h 20 m

| stage | route | FSR (nm) | dips |
|---|---|---|---|
| 1 | monolithic, **calibrated** | 12.560 | 7 |
| 2 | cascaded, **calibrated** | 12.753 | 8 |
| 3 | monolithic, uncalibrated | 12.560 | 7 |
| 4 | cascaded, uncalibrated | 12.752 | 8 |

**Control, and a good one: stage 4's 12.752 nm is exactly the FSR this device's
own validation document records for the post-M1 3-D cascaded route (job
11658568, 2026-08-03).** Every cascaded FDTD domain was re-simulated on this run
under the new translation-normalised cache key, and it reproduced the old physics
to the printed precision — the re-keying turned entries over without changing an
answer, which is what a re-keying is supposed to do.

**FSR barely moves; the COMB does.** Judging the calibration by FSR would have
found nothing (12.560 → 12.560, 12.752 → 12.753). That is correct behaviour, not
a null result: the excess is nearly wavelength-independent (+0.998 % … +1.044 %
across the band), so the comb translates almost rigidly while the *spacing*,
set by the group index, barely changes. The measurement is the dip positions:

| route | comb shift | as a fraction of an FSR |
|---|---|---|
| monolithic | **+0.392 nm** | **+3.12 %** |
| cascaded | **+0.951 nm** | **+7.45 %** |

### Two predictions in this document were wrong, and both are instructive

**1. The "roughly a quarter of an FSR" estimate was wrong, and so was the
script's own printed +11.0 %.** Both used `L_OUTER = 7.00 µm` as the lever arm.
**n_eff enters the round trip twice with opposite signs:**
`P = b_corrected² · exp(1j·phi)` with
`b_corrected = b_d · exp(-1j·phi_bend_correction)`, so the analytic phase is
`phi − 2·phi_bend_correction` and the **net** arm is
`L_OUTER − BEND_STRAIGHT_TOTAL = 7.00 − 5.00 = 2.00 µm`. The net arm predicts
**+3.15 %** against **+3.12 %** measured; the discrepancy was exactly
`L_OUTER / (L_OUTER − BEND_STRAIGHT_TOTAL) = 3.5`. Fixed in
`ring_resonator_hybrid.py` to the same expression the cascaded script's print
already used correctly. **The printed line in job 11983058's log is wrong; the
comb it produced is not** — this only ever fed a diagnostic print, never the
assembly, so no re-run is owed.

**2. "The two routes' analytic arms differ by 3.5×, so the calibration will not
shift them equally — a real prediction this run tests."** Right conclusion,
wrong at every step of the reasoning. The arms do **not** differ: both are
2.00 µm net, which this document already knew from the other direction
(`max |P_phase-strip − P_DC-style| = 7.2e-15`, "not a physics difference"). And
the ordering came out **backwards** — the cascaded route moved *more* (+7.45 %)
than the monolithic (+3.12 %), not 3.5× less.

The real asymmetry is that **the cascaded route has a second calibration.** Its
`P_body` supermode phase is calibrated too, and that term is the difference:
`7.45 − 3.12 = +4.33 %` of an FSR on top of the strip term the routes share. The
strip term alone reproduces on both at +3.12 / +3.15 %.

**So the half-FSR route offset is NOT closed.** Monolithic 12.560 nm against
cascaded 12.753 nm; the calibration *changed* the offset by ~4.3 % of an FSR
rather than removing it. `ring_resonator_validation.md`'s `## Order` item 1 stays
open, with one more term now accounted for.

### `P_body` reproduces item 6's amplification, and independently validates the new code

    Body supermode n_eff at lambda = 1.550 um (even/odd):
        MPB mode index  2.349085 2.326731   delta_n_eff 0.022353
        FDTD own phase  2.375605 2.354529   delta_n_eff 0.021076   <- P_body propagates on this
        per-mode excess +1.129% / +1.195%, but delta_n_eff -5.714%

Per-mode excesses differing by 0.066 percentage points produce a **−5.714 %**
change in `Δn_eff` — the `β/Δβ` lever, on this device, measured.

And a check that was not designed, which is worth more than one that was.
`directional_coupler_hybrid.py` is **SOI_220, `waveguide_width = 0.45`,
`gap = 0.25`, res 40, C-band** — the *same physical cross-section* as the ring's
coupled bus+ring pair (`W = 0.45`, `GAP = 0.25`, `SOI_220`, res 40). Two
independent code paths: one through `HybridSolver(calibrate_body_neff=True)` on a
`DirectionalCoupler` geometry at `oxide_width = 6.0`, the other through this
session's new `UniformCrossSection` calibration inside a raw-meep script at
`BODY_MESH_WIDTH = 4.70`:

| | even | odd | Δn_eff |
|---|---|---|---|
| coupler (job 11981703_2) | 2.375702 | 2.354625 | **0.021077** |
| ring pair (job 11983058) | 2.375605 | 2.354529 | **0.021076** |
| \|Δ\| | 9.7e-05 | 9.6e-05 | **1e-06** |

**Δn_eff agrees to 1e-06 (0.005 %)** while the individual modes differ by ~1e-04
— the expected signature rather than a coincidence: the residual per-mode offset
is a common-mode effect of the different transverse window (6.0 vs 4.70 µm) and
cancels in the differential, which is the quantity that matters. **This is the
strongest available evidence that the new ring calibration measures the right
thing**: it reproduces a separately validated device's FDTD β on a shared
cross-section, through a different geometry class, a different solver
configuration, and a different script.

---

# Step 11 — the last uncalibrated segment, and a cache-key bug found on the way there (2026-08-18)

Step 9's rollout table recorded `hcg_fabry_perot_hybrid` as needing **"none —
calibrates nothing itself"**, and step 10 re-ran it unchanged. Both were true of
the script as it then stood and both missed the same thing, for the same reason
the ring was missed: *a device that reaches no calibration cannot be audited by
looking at which calibrations it reaches.*

That script was a `beef.network.Circuit` assembly of two saved pickles plus one
`meow`/MPB mode solve. The mode solve was the **connecting waveguide** —
2 × 53.08 µm of straight strip, against the coupler's 23 µm — i.e. **~106 µm of
the ~130 µm one-way optical path was propagating on MPB's index while every
block it cascaded with was referenced to meep's β.** Same argument as the ring's
`phi`, with a much longer lever: the resonant denominator is `1 - ρ²t²` with
`ρ = r·phase²`, so the arm index acts over a **4 × L_ARM = 212.32 µm** round trip.

## What changed

`examples/high_contrast_grating_fabry_perot/python/hcg_fabry_perot_hybrid.py`
now solves all three pieces itself against `example_sim_cache/` — the mirror's
`solve_reflection()`, the coupler's two `solve_section()` S-bends and its
`section_mode_beta_details()` body calibration, and a **new
`section_mode_beta_details()` calibration of the arm**, reached through
`UniformCrossSection` exactly as the ring's straights are (step 10, finding 2).
`CALIBRATE_ARM` defaults on; `HCGFP_CALIBRATE_ARM=0` inverts the output suffix.
It no longer reads `smatrix.pkl`/`smatrix_hybrid.pkl` from either constituent
script's output folder.

**The rebuild is a refactor, and that is measured, not argued.** The script
compares its own rebuilt coefficients against those two published pickles
whenever they are still on disk, and both agreed on the first run:

| rebuilt in `hcg_fabry_perot_hybrid.py` | vs published | max abs difference |
|---|---|---|
| HCG mirror `r(λ)`, 21 pts | `high_contrast_grating_fdtd/smatrix.pkl` | **0.000e+00** |
| coupler `t(λ)`, 11 pts | `directional_coupler_soi160_hybrid/smatrix_hybrid.pkl` | **1.85e-15** |
| coupler `k(λ)`, 11 pts | same | **1.88e-15** |

The coupler's are float noise rather than exact because the body is now reached
through `EMESolver.solve_section()` at the one real length instead of
`solve_section_length_sweep()` at 200; with `neff_override` both take the same
diagonal fast path, so the difference is summation order.

## The cache-key bug — `minimum_run_time`, found because a cache hit was expected and missed

Building the above, the mirror's `solve_reflection()` call **missed** a cache
entry that visibly matched it (`FDTDSolver_solve_reflection_GdsGeometry_
993cdeb82669a277`: 21 pts, 1.50–1.60 µm, res 60, 3-D, `extend_input_waveguide`).
Diffing the built fingerprint against the stored `config.json` left exactly one
difference — `solver_config[minimum_run_time]: mine=0.0, cached=<absent>` — and
deleting that one attribute off a freshly built solver made the same call hit,
returning the real result (peak R 0.9714 @ 1.580 µm, matching job 11673551).

`FDTDSolver.minimum_run_time` was added in commit `d7910e8` with the modal
multiplexer's early-stop fix. It has **two independent faults**, in opposite
directions, and the second is the dangerous one:

| | mechanism | consequence |
|---|---|---|
| `solve_reflection` | no `_METHOD_ATTRS` entry → full `vars()` dump → the new key has no counterpart in older configs | **48 stranded entries**, including a ~55 min res-60 3-D production run that was still perfectly valid — the runs it described had no floor, which is exactly what the new default `0.0` means |
| `solve`, `solve_section`, `solve_section_multi_length` | all pass it to `stop_when_dft_decayed()`, all **omitted it from their allowlists** | a **silent false-positive hit**: a floored run reusing an unfloored result. `_METHOD_ATTRS`' own comment calls this "strictly worse than the spurious-miss cost of over-including" |

The second was **live, not latent**: `modal_multiplexer_hybrid.py` sets
`minimum_run_time=100.0` and `modal_multiplexer_fdtd.py` sets `1000.0`, both
added 2026-08-08 *because* a loose `decay_tol` stopped a 66 µm device at
t≈107–112, before the pulse could physically cross it. A key that cannot see the
floor can serve the pre-fix result to the post-fix request — defeating the fix.
(`decay_tol` **is** in the allowlists and changed at the same time, which
narrows but does not close the exposure: `MUX_FDTD_DECAY_TOL=1e-11`, documented
as the way to "reproduce the two prior (broken) runs exactly", lands on exactly
the ambiguous combination, and `modal_multiplexer_hybrid.py`'s own
`decay_tol=1e-11` never changed at all.)

### The fix — `_DEFAULT_OMITTED_ATTRS` (`beef/cache/config.py`)

An attribute in that set is **left out of the fingerprint while it holds its
constructor default, and tracked normally when it does not.** Defaults are read
off `inspect.signature`, not duplicated. Both faults close at once:

- entries written before the attribute existed still hash equal to a
  default-valued solver today → the 48 `solve_reflection` entries un-strand, and
  the HCG mirror did **not** have to be re-simulated;
- `minimum_run_time` is now listed in all three allowlists that read it, so a
  non-zero floor is a different cache entry — the modal multiplexer's 100.0 and
  1000.0 stop being invisible.

**Deliberately not "skip every default-valued attribute".** Most tracked
attributes *are* at their defaults in real runs (`material_stack`, `dpml`,
`decay_tol`, `sim_height`, `eps_averaging`, the subpixel pair) and every existing
entry records them; dropping those would strand the whole cache — the exact
damage being undone. The set is explicit and additive: an attribute belongs there
only because it postdates entries that are still wanted.

**Also deliberately not done**: giving `solve_reflection` its own `_METHOD_ATTRS`
entry. It genuinely over-includes today (`num_input_modes`, `num_output_modes`,
`output_monitor_centers`, `output_monitor_widths`, `half_domain_aperture`,
`waveguide_length` are all unread by it), but narrowing the key would strand
those same 48 entries a second time. Over-inclusion only costs spurious misses,
which is the safe direction; left as a known, deliberate wart.

**Verified**: `tests/test_cache_config.py::TestDefaultOmittedAttrs` (7 tests,
4 of them parameterised across all four affected methods) pins both directions,
plus the `False == 0` / `0 == 0.0` type-exactness the default check needs. Full
suite 670 passed / 19 skipped. Measured, not argued: the same read-only
`cache_only=True` probe that found the miss was re-run after the fix and the
mirror, both S-bends, the coupler body calibration and the EME body all hit.


---

# `decay_method` and a sequence-position stall — found on the 4×4 hybrid, checked against every other calibrating device (2026-08-20)

Calibrating the (not yet shipped) 4×4 90° hybrid's body — the first device to
exercise `section_mode_beta()` at as many as 16 modes — turned up **two
separate problems**, only one of which this doc's existing machinery already
covered. Recorded here because the fix for the first is now a general
`FDTDSolver` option any device can reach for, and the second is an unresolved
risk this doc's own "is it safe" audits (step 10 above) should be able to
check against going forward.

## Problem 1 — mode 1 stalled under `stop_when_dft_decayed()`, fixed

`stop_when_dft_decayed()` (meep's default, and until now this library's only,
early-stop criterion) requires *every* tracked DFT component to have decayed
to `decay_tol` before it will end a run — with a floor from `minimum_run_time`
but no ceiling. Mode 1's own body-calibration run hit a component that would
not decay to tol in any reasonable time (job 12052383/12052385's stall).

**Fix**: `FDTDSolver` now takes `decay_method: str = "dft"` and
`energy_decay_dt: float = 20.0` (`beef/solvers/fdtd.py`, `__init__`), with
`decay_method="energy"` switching every `sim.run()` call in the solver to
`mp.stop_when_energy_decayed(energy_decay_dt, decay_tol)` via a new
`_stop_condition()` helper. Both new attributes are in `_METHOD_ATTRS` for
`solve`/`solve_section`/`solve_section_multi_length`, and in
`_DEFAULT_OMITTED_ATTRS` (`beef/cache/config.py`) — same two-part pattern as
`minimum_run_time` above, checked the same way (default-valued solver still
hits the pre-existing cache; a `decay_method="energy"` solver's config key
differs and misses correctly). Mode 1 decays cleanly under this criterion.

## Problem 2 — the 16th sequential mode stalled too, NOT fixed by the above, root cause NOT found

With problem 1 fixed, the same 16-mode run reached mode 16 — last in
`solve_section_multi_length()`'s own `for k in range(1, n_excite + 1):` loop,
i.e. the 16th `mp.Simulation` build + run + `reset_meep()` in one long-lived
process — and **stalled again**, under `decay_method="energy"` too. Isolating
mode 16 alone in a fresh process (patching the loop to `for k in [16]:` via
`inspect.getsource`/`exec`, job 12066668) converged normally: t=220.1375,
17,611 timesteps, energy ratio 1.0 → 7.5e-12. Mode 16 has no intrinsic
problem — margin from cutoff is 0.883, nothing near the guard thresholds.
Something accumulates across repeated `mp.Simulation`/`reset_meep()` cycles
within one process and eventually stalls the last one; **what, specifically,
is not identified**, only worked around.

**Workaround, not a fix**: one process per excited mode.
`examples/verification/python/mmi_4x4_90hybrid_calibration_mode_worker.py`
(one mode per SLURM array task, each a fresh process) +
`_calibration_aggregate.py` (combines the 16 `mode_XX.npz` files exactly as
`section_mode_beta_details()` would in one process — `beta_fit_from_phases()`
depends only on that mode's own phase-vs-length data and `seed_beta`, so
per-mode isolation changes nothing about what gets computed, applying the same
`_BETA_FIT_RESIDUAL_MAX`/`_BETA_BRANCH_MARGIN_RAISE` guards on the far side).
Reusable as a pattern if another device hits the same wall.

## Checked: is any other calibrating device close enough to be at risk?

Every device currently calibrating its own body (`calibrate_body_neff=True`,
same call path this doc covers) was checked for its own sequential mode count
— the dimension that matters, since the stall only appears deep into one
process's own sequence of excited-mode simulations:

| device | own-body N_MODES | where set |
|---|---|---|
| MMI 2×2 | 14 | `examples/mmi_2x2/python/mmi_2x2_hybrid.py:126` |
| MMI 1×4 | 12 | `examples/mmi_1x4/python/mmi_1x4_hybrid.py:249` |
| MMI 1×2 | 8 | `examples/mmi_1x2/python/mmi_1x2_hybrid.py:148` |
| directional coupler | 5 | `examples/directional_coupler/python/directional_coupler_hybrid.py:249` |
| directional coupler, SOI160 | 2 | `examples/directional_coupler/python/directional_coupler_soi160_hybrid.py:211` |
| HCG Fabry-Perot, coupler body | 2 | `examples/high_contrast_grating_fabry_perot/python/hcg_fabry_perot_hybrid.py:210` |
| HCG Fabry-Perot, arm | 1 | `hcg_fabry_perot_hybrid.py:493` |
| ring resonator (cascaded), coupled-pair body | 2 | `examples/ring_resonator/python/ring_resonator_hybrid_cascaded.py:296` |
| ring resonator, bare strip (both scripts) | 1 | `ring_resonator_hybrid.py:936`, `ring_resonator_hybrid_cascaded.py:1863` |

None of these have been observed to stall, and 14 (the highest, MMI 2×2,
step 10's own re-run — jobs 11981703/11982196) is not close to 16. **That is
not the same as "safe below some known threshold"** — the trigger condition
for problem 2 is confirmed only at position 16 of 16, on this one device's
geometry and resolution, and nothing here pins down whether it depends on
mode count, process lifetime, resolution, or something else entirely. Nothing
above needs to change today. If any future device's own `N_MODES` grows past
~14, or any script runs many sequential `mp.Simulation`+`reset_meep()` cycles
in one process for a different reason (a length sweep, a wavelength sweep),
and that job runs far longer than its own per-mode/per-length cost model
predicts, treat that as this same failure mode rather than a slow machine —
`decay_method="energy"` and the one-process-per-mode pattern above are both
now available as the first two things to try, in that order.

## Correction, same day — mode 16's own isolated run diverges too; this is not a sequence-position bug

The paragraphs above were written right after job 12066668 (mode 16 alone,
in its own process) appeared to converge cleanly (t=220.1, ratio 7.5e-12).
Running the full one-process-per-mode production array to actually get a
usable β for all 16 modes (job 12069540_1 + 12069773_[2-16], all at the
real 3-wavelength production seed/lengths) falsifies that conclusion.

**Modes 1-15 completed cleanly** (12-30 min each). **Mode 16
(`12069773_16`) did not converge — it hit `TIMEOUT` at the full 90 min
budget, exit `0:15` (SIGTERM, not a crash)**, still running at meep-time
~2100 when killed. Its own energy trace, read directly from the job's log
rather than assumed, is not a stall — it is **exponential growth**: the
ratio bottoms out around t≈260 at 1.95e-11 (just *above* `decay_tol=1e-11`,
so `stop_when_energy_decayed()` correctly did not fire), then climbs
monotonically — 1.0 again by t≈580, 1.7e53 by t≈2081, still climbing at
the last checkpoint before timeout. This is a genuine numerical
instability developing at late simulated time, not a criterion that is
merely slow to satisfy.

**That reframes job 12066668's own "clean" result as a near-miss, not a
resolution.** Its own energy trace (checked directly, same log) was *also*
still monotonically falling with no turnaround visible through its own
last checkpoint (t=220, ratio 7.5e-12) — it happened to cross below
`decay_tol` on a very slightly longer domain (`max(lengths)=2.655 µm`, the
single-wavelength diagnostic lengths) than the real production run
(`2.525 µm`, this device's 3-wavelength β_max), which was enough to change
exactly when the noise floor bottoms out relative to when this instability
overtakes it. **Isolating mode 16 into its own process changes nothing
about whether it converges — the earlier "no intrinsic problem, only a
sequence-position bug" conclusion is retracted.** The one-process-per-mode
architecture itself is still sound and was not wasted: it correctly
isolated mode 16's own data from any cross-mode interaction, and that
isolation is precisely what proves the instability is intrinsic to mode
16's own configuration, not an artifact of running 16th in a sequence.

**Root cause, grounded in the repo's own code, not guessed**:
`beef/materials/stack.py:284/337` builds every Sellmeier resonance term as
`mp.LorentzianSusceptibility(frequency=1/lam_i, gamma=0, sigma=...)` — a
**lossless (zero-damping) pole**, the same convention meep's own
`materials.Si` uses. Undamped Lorentzian susceptibilities are a
documented marginal-stability risk in FDTD over long integration times —
accumulated discretization noise can resonantly grow on an undamped pole
given enough timesteps. This has apparently never surfaced elsewhere in
this repo because nothing else runs this long at this tight a
`decay_tol`: ordinary S-parameter solves neither integrate for ~2000+
time units nor demand 1e-11 decay. Calibration's own demands are what
expose it.

**This changes the "which other devices should worry" answer from the
section above.** That section reasoned about risk in terms of sequential
mode count (N_MODES), because the believed mechanism was resource
accumulation across repeated `mp.Simulation` calls. The real mechanism —
an undamped dispersive pole slowly resonating over long simulated time —
has nothing to do with mode count or process lifetime; it depends on how
long a run has to integrate and how tight its decay tolerance is, which
calibration's own extreme requirements make far more likely to trigger
than a mode-count comparison alone would suggest. No other device's
calibration is known to have hit it, but "N_MODES well under 14" is not
a reason to assume it's safe.

**Open, not resolved**: whether a small nonzero `gamma` on the dispersion
model suppresses this without a material accuracy cost; whether
`decay_method="dft"` behaves differently for mode 16 specifically (the
original all-in-one-process run under `dft` stalled on mode 1 before ever
reaching mode 16, so this is untested); whether a shorter calibration
domain or looser `decay_tol` sidesteps it at an acceptable accuracy cost.
**Not pursued further without direction** — this is a materially
different, more open question than "isolate each mode and it'll work,"
and deserves a decision on how far to chase it before more cluster time
goes into it.



---

# Rank count: calibration and full-device FDTD want OPPOSITE resourcing (measured 2026-08-22)

Prompted by the user asking whether all cores are actually used at every
stage of a simulation. They are not, and the answer differs so sharply by
job type that using one rule for both wastes either time or allocation.

## Phase breakdown, from this repo's own logs

Wall time split into time-stepping (parallel, domain-decomposed) vs
everything else (structure setup, the MPB source solve, eigenmode
decomposition, IO):

| job | non-stepping | parallel fraction | Amdahl ceiling |
|---|---|---|---|
| 4×4 hybrid full-device FDTD (12079104) | **4.2%** | 95.8% | ~24× |
| calibration worker, mode 2 (12069773_2) | 47.3% | 52.7% | 2.1× |
| calibration worker, mode 8 | 65.1% | 34.9% | 1.5× |
| calibration worker, mode 14 | 68.3% | 31.7% | 1.5× |
| mode-16 isolation (12066668) | 77.6% | 22.4% | 1.3× |

This independently reproduces step 9's "MPB eigenmode decomposition is
~89% of calibration cost" from a different measurement. Note the trend:
the non-stepping share **rises with mode index**, consistent with
higher-order modes needing more MPB iterations in a dense near-degenerate
spectrum.

## The MPB phase does not scale AT ALL — measured, jobs 12218118 / 12218119

`section_mode_neff()` (pure MPB, no time-stepping) on the 4×4 body
cross-section at production `n_modes`, deliberately with **no `save_dir`**
so it could not return a cache hit and measure nothing:

| ranks | wall |
|---|---|
| 8 | 3110.7 s |
| 24 | 3055.2 s |

**3× the cores, 1.8% difference.** meep solves the MPB eigenproblem
redundantly on every rank; N−1 of N cores duplicate work rather than
contribute.

## What follows

- **Full-device FDTD: add ranks.** 96 ranks measured 1.95× over 48 at 98%
  efficiency on the 4×4 (job 12217184). 192 across scattered racks never
  reached time-stepping in 30 min (12217185) — 96 is the sweet spot.
- **Calibration: do NOT add ranks.** The 16-task array used 24 ranks/task
  (384 cores) and could never have exceeded ~1.67× regardless. At 8
  ranks/task, mode 8 would run ~1.7× slower in wall time but consume
  ~1.8× fewer core-seconds (33.4k → 18.9k). Given how often this account
  hits `AssocGrpBillingMinutes`, that trade is usually worth taking.
- **Optimize calibration by reducing MPB work, not by adding cores** —
  which is exactly what the two fixes that worked did
  (`measure_diagonal_only` ~4×, `seed_anchors=3` ~10× on the seed step).

## A load-bearing comment that is wrong when uncached

`mmi_4x4_90hybrid_calibration_mode_worker.py` calls `section_mode_neff()`
"unpatched -- cheap and already cache-safe" / "cheap/cached, not the
stalling part". Uncached it took **51 minutes** here. The claim is true
only *because* it is normally a cache hit; it is not cheap in itself, and
nothing in that comment says so.
