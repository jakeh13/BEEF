# Ring resonator validation — to-do

> ## 📁 FILED INTO `resolved/` — 2026-08-03
>
> Per this directory's README, `resolved/` means *"whole docs sufficiently
> addressed (open to-dos may remain)"*. It is **not** a claim that everything
> below is closed. Read this box before quoting anything from this file.
>
> **The validated route is the CASCADED one** (`ring_resonator_hybrid_cascaded.py`).
> Its basis, all measured:
>
> | check | result |
> |---|---|
> | 1 — engine agreement | bare strip bit-identical (`0.000e+00`); coupled pair `1.3–1.6e-10` = the translation floor exactly |
> | 2 — mode identity | field overlap at the translation floor |
> | unitarity `\|t\|²+\|k1\|²` | **0.999640 – 0.999926**, deviation **3.60e-04**, on an FDTD-measured `M_left` — post-M1 3-D, job 11658568 |
> | FSR convergence | p ≈ 1.4–2.1 |
> | mode-solve resolution | converged at res 100; res 40 → 100 moves the comb **3.5–3.7%** of an FSR (job 11658491) |
> | `interp_complex`, degree-2 `n_eff` polyfit | both measured, both fine |
>
> **✅ The basis is stable under M1** — the material fix did not move it, which is
> the outcome that makes it usable. Post-M1 3-D (11658568, 1h05) against the
> pre-M1 numbers: unitarity deviation **3.595e-04** vs 3.5e-04; `|t_ring|` vs `|t|`
> **3.774e-03** vs 3.768e-03; reciprocity **1.665e-16** vs 2.220e-16. 8 resonances,
> mean FSR **12.752 nm** against 12.87 predicted from the 43.42 µm perimeter
> (0.9%). Supermode `n_eff` fit residuals 1.45e-05 (even) / 7.30e-07 (odd),
> bare-strip 7.29e-06. The assembly ledger reports the corrected DC-style spans
> (0.25 / 1.50 / 0.25, each matching its physical span to `+0.0000`).
>
> **⛔ OUT OF SCOPE BY DECISION, not by omission:**
> * **Absolute resonance wavelength.** Finesse `F = Q·FSR/λ ≈ 120`, so any
>   per-round-trip phase error is amplified ~120×. Neither model positions the
>   comb to better than a linewidth. FSR and coefficient-level quantities are
>   validatable; absolute comb position is not.
> * **Comparison against a full-device FDTD.** Retired 2026-08-03 at the user's
>   direction and independently disqualified by its own output: it reports
>   transmission **> 1** for a passive device over **46%** of the band, and the
>   excess *grows* with refinement (+0.067 at res 30 → +0.084 at res 40), so it is
>   not a coarse-grid artifact. See the ⛔ RETIRED section below.
>
> **🔴 OPEN, carried into `resolved/` deliberately:**
> 1. **The MONOLITHIC route (`ring_resonator_hybrid.py`) has a real, unexplained
>    defect** — a k1/k2 residual of `1.320e-04` with `arg` drift `0.0068 π` that no
>    plane correction removes and **six tested mechanisms do not explain**. The
>    cascaded route does not have it (`2.220e-16`). Do not quote monolithic comb
>    positions. Its 3-D output is also still pre-M1, deliberately not re-run.
> 2. **⚠ The cascaded route's reciprocity check is a TAUTOLOGY, not evidence.**
>    `S_coupler = M_left.T @ P_body @ M_left` makes `k1` and `k2` the same matrix
>    element, so `1.665e-16` is machine epsilon — the route *enforces* reciprocity
>    rather than verifying it. Unitarity above is the real coefficient-level test.
>    The script itself used to print this as *"a stronger guarantee than the
>    monolithic version's independently-normalized `|k1|~=|k2|`"*, which is
>    backwards — an identity is the weakest possible evidence, and that wording is
>    why it was briefly recorded as part of the validation basis. Corrected in
>    source 2026-08-03; the entry is withdrawn.
> 3. **Check 3's resolution sweep is pre-M1.** The `p ≈ 0.70` comb-position and
>    `p ≈ 1.4–2.1` FSR numbers come from job 11653509 (ended 12:33–14:40); M1
>    landed ~15:00, and those cache configs lack the `resolved_dispersive` key
>    entirely. Re-run **RUNNING as 11658729** (released when 11658568 completed);
>    pre-M1 data archived at `archive/ring_res_sweep_pre_m1_20260803/`. **Update
>    the numbers here when it lands.**
> 4. ~~**Post-M1 3-D production**~~ — **✅ CLOSED, job 11658568, 2026-08-03,
>    1h05.** The production directory is now internally consistent: post-M1
>    assembly on post-M1 FDTD. Numbers folded into the table above.
> 5. **Bookkeeping:** `_CACHE_SCHEMA_VERSION` is still 1 and the monitor-plane
>    convention change wants the bump; **no ring golden exists at all**
>    (`tests/golden/` has none, so nothing in `test_golden_examples.py` /
>    `test_golden_provenance.py` covers this device); and the two routes now ship
>    different bend conventions — phase-strip in the monolithic
>    (`ring_resonator_hybrid.py:1148`), DC-style in the cascaded
>    (`ring_resonator_hybrid_cascaded.py:1657`). Those two are algebraically
>    identical (`max |P₁−P₂| = 7.2e-15`), so this is an inconsistency between two
>    files meant to be compared, not a physics difference.

> **2026-07-28 — the material model changed under every device here.** Both engines
> now evaluate one analytic Sellmeier model (Salzberg Si, Malitson SiO2) instead of
> each doing its own nearest-tabulated lookup, which shifts n_Si by ~0.002 and
> removes a 0.0063 step that sat in the middle of the C-band. Any golden or stored
> result for this device predates that and needs regenerating; see
> `directional_coupler_validation.md` item 8 and, for the dispersive-FDTD
> re-parameterisation (built and gate-passed 2026-07-31, opt-in via
> `dispersive=True`, off by default), `material_model_validation.md`.



Opened 2026-07-26. **The least validated device in the repo.**

**Status 2026-07-26:** the axis bug is **fixed in code** (`directional_coupler_validation.md`
item 10 — see Done below). Nothing has been re-simulated yet, so every ring
plot and every stored ring FDTD result on disk is still pre-fix. The fix
guarantees the old data can no longer be served silently, not that new data
exists.

`examples/ring_resonator/python/` — `ring_resonator_hybrid.py`,
`ring_resonator_hybrid_cascaded.py`. Despite the `_hybrid` names, neither uses
`HybridSolver`, `FDTDSolver` or `EMESolver`: both hand-build their meep
simulations and their meow cross-section directly. That single fact is the root
of everything below.

## Done — the axis bug (`directional_coupler_validation.md` item 10), 2026-07-26

Both scripts computed `fwidth = (1/WL_MIN - 1/WL_MAX) * 0.5` and handed it to
`sim.add_mode_monitor(fcen, fwidth, N_WL_FDTD)` — character for character the
defect item 1 fixed, in raw meep rather than through `FDTDSolver`, which is
exactly why the solver fix never reached it. A request for 1.50–1.60 µm
actually sampled **1.5254–1.5754 µm**: half the span, reversed.

For a resonator this is worse than for a coupler. A reversed, half-span axis
does not merely mislabel a smooth curve — it puts every **resonance
wavelength** in the wrong place and compresses the free spectral range by 2×.
No FSR or Q read off a pre-fix ring plot is meaningful.

**Changed**, in both `ring_resonator_hybrid.py` and
`ring_resonator_hybrid_cascaded.py`:

- `add_x_monitor()` / `add_y_monitor()` now take an explicit ascending
  frequency list from `monitor_freq_grid()` instead of `(fcen, fwidth, nfreq)`.
- `fwidth` comes from `source_fwidth(wavelengths, SOURCE_BANDWIDTH_FACTOR)`
  with the factor pinned at `FDTDSolver`'s own default of 2.0 — **exactly 4× the
  old pulse** (2× for covering the full span rather than half, 2× again for the
  bandwidth factor). Widening the monitors without widening the source would
  have starved the new band edges, which is why item 2 was needed alongside
  item 1.
- Every measured amplitude is now reordered onto the caller's axis with
  `[wl_order]`. Each ratio is taken in *monitor* order first — numerator and
  denominator are both monitor-ordered, so they stay aligned — then reordered
  once, as `FDTDSolver.solve_section()` does. In the cascaded script the
  reorder is along axis 0 of the `[n_wl, N_MODES]` interface column.
- `check_source_coverage()` now runs on every FDTD entry point.
- `_sim_params()` takes `wavelengths` and emits `resolved_monitor_freqs` /
  `resolved_source_fwidth`, mirroring
  `FDTDSolver._resolved_spectral_config()`.

**The cache stranding was verified against the real cache, not assumed.** All
**6** pre-fix FDTD entries now miss — each folder's stored `config.json` was
re-hashed to confirm it matched its own folder name first, so the comparison is
against the real recorded config rather than a reconstruction:

| entry | stored hash | after fix |
|---|---|---|
| `ring_resonator_hybrid` `run_coupler_fdtd` bus | `6ec67e004cc0d073` | `50ca0d12f91a1f76` |
| `ring_resonator_hybrid` `run_coupler_fdtd` ring | `7dc2d9ae7706f4b0` | `3dd5aceabde0d1ee` |
| `ring_resonator_hybrid` `run_bend_fdtd` | `9f79cce29bef6ee2` | `88188f3c1262fe16` |
| `..._cascaded` `run_left_fdtd` bus | `f90d6cc8c4f64816` | `ac3912345252c8c4` |
| `..._cascaded` `run_left_fdtd` ring | `c5566a9725a64929` | `e9402786ee01387b` |
| `..._cascaded` `run_bend_fdtd` | `6c2f7d6c5669e69c` | `d347ce2dde611894` |

`lookup()` returns `None` for all six. They are not deleted — deleting only
forces re-simulation of equally wrong results, and the stranded folders are
still readable if the pre-fix numbers are ever wanted for a relabelling check
like item 5's.

**Equally important: the 39 mode-solve entries are untouched** — all
`solve_neff` and `solve_body_supermodes` folders still hit. Those are MPB
solves that the monitor grid never affected, and neither calls `_sim_params()`,
so nothing in this change can reach their hash. Confirmed 39/39 reachable.

**Tests** — `tests/test_ring_monitor_grid.py`, 14 fast tests, no meep. They are
source-level because importing either script *runs* it. Five guard the scripts
(no `(fcen, df, nfreq)` triple, no half-span `fwidth`, helpers used, `wl_order`
actually applied, `_sim_params()` fingerprints the resolved grid), with
docstrings and comments stripped first so the prose above describing the old
line cannot pass for live code. Four pin the arithmetic — including a
reconstruction of the pre-fix band from meep's documented rule, so the
1.5254–1.5754 figure is derived in the test rather than asserted.

**Not done: nothing has been re-simulated.** The fix guarantees a re-run cannot
silently serve pre-fix data; it does not produce post-fix data. Every ring plot
in `example_output/` is still pre-fix and should be regarded as wrong until
regenerated. When that re-run happens, the check to apply is item 9's: **the
centre wavelength must not move**, since `fcen` is the one point the pre-fix
code labelled correctly. `fcen` itself is unchanged by this work — `wl_center`
switched from `(WL_MIN + WL_MAX)/2` to `np.mean(wavelengths)` to match
`FDTDSolver._build_domain()`'s convention by construction rather than by
coincidence, and for the uniform `linspace` both scripts pass these agree
exactly at N=101.

## Recommendation — put the ring on `HybridSolver`

Item 10 exists *because* this device is off the library path. Items 1–4, 11 and
12 all landed in `beef/` and reached every other example automatically; the ring
needed each ported by hand, and the axis fix was not. That will keep happening.
The durable fix is to make the ring a first-class device, which needs three
things — the first two mechanical, the third a real design change.

**(a) Arbitrary-orientation ports in `FDTDSolver` — done, 2026-07-29.** Every
source and monitor in `beef/solvers/fdtd.py` was x-normal — `mp.Vector3(0, w,
h)`, 8 occurrences — and the ring's two coupling ports, sitting on top of the
near bends, are **y-normal** (x–z planes). `solve_section()` now takes
`y0`/`y1` (a y-axis domain bound, mirroring `x0`/`x1`) and `excite_port`/
`measure_port` (`PortSpec(axis="x"|"y", sign=±1, center, width)`), additive
and default-`None` — every existing caller is unaffected, verified by the
full fast suite staying at 338 passed after the change.
`section_fdtd_layout()`/`SectionFDTDLayout` are untouched; `PortSpec`
bypasses them per-side rather than extending their scalar contract.

Two real bugs found while building this, not just anticipated: (1)
`section_fdtd_layout()`'s `reverse = n_out <= n_in` tie-break silently
excites from the wrong physical side whenever `n_in == n_out` — exactly the
standalone-bend case (`n_modes_in=n_modes_out=1`) — now bypassed via a forced
`reverse=False` whenever either `PortSpec` is given, with the excited/
measured role taken literally from `n_modes_in`/`n_modes_out`. (2) a
measure port's edge convention is the *opposite* of an excite port's for the
same `sign` value — a source's launch direction physically forces it to a
specific edge, but a passive monitor has no such constraint, so `sign=+1`
means "low edge, launch +axis" for `excite_port` but "high edge, read
+axis-forward light" for `measure_port`. Caught mid-implementation, before
either shipped silently wrong.

**Verified against real physics, not just unit tests.** A
`@pytest.mark.slow` check (`tests/test_ring_orientation.py`, run with
`--runslow`) builds the new `RingResonatorBendSection` (below), excites
x-normal, measures y-normal through the 90° bend at reduced resolution
(res 15, 2-D — a login-node smoke run, not production), and passes: meep's
own log shows the y-normal monitor resolving a mode whose dominant
planewave is genuinely y-directed (`(0, 2.01, 0)`) against the x-normal
source's x-directed one (`(2.04, 0, 0)`) — physical confirmation, not a
numeric coincidence. **Not yet run at production resolution/3-D.**

> ### ⚠️ NEW BLOCKER, 2026-07-31 — the y-normal *excite* path fails reciprocity
>
> The y-normal **excite** + flux-calibration path (needed for the coupler
> section's ring input; not exercised by the x-excite/y-measure test above)
> has now been tested for the first time, and **it does not satisfy
> reciprocity.** New test:
> `tests/test_ring_orientation.py::test_bend_y_normal_excite_matches_reciprocity`
> — same bend, same mesh, same apertures, only which port excites changes.
> A passive reciprocal structure must give S₂₁ = S₁₂, so |T| must match.
> Measured at res 15, 2-D, λ = 1.53/1.55/1.57:
>
> | | 1.53 | 1.55 | 1.57 |
> |---|---|---|---|
> | \|T\| x→y (validated direction) | 0.71439 | 0.71114 | 0.70866 |
> | \|T\| y→x (**y-excite, flux-calibrated**) | 0.51147 | 0.52247 | 0.51708 |
> | ratio | 0.7160 | 0.7347 | 0.7297 |
>
> **The ratio is geometry-invariant, which is what makes this a real finding
> rather than mesh noise.** First hypothesis was that the test itself was
> unfair — the x port has a 1.0 µm straight lead-in (`stub_left=1.0`) while
> the y port had only `output_margin=0.5`, so the two monitors sit at
> different distances from the arc and see differently-contaminated
> near-fields. **Refuted by direct sweep** (`output_margin` 0.5 → 1.0 → 1.5,
> λ=1.55): |T| x→y itself moves 27% (0.71114 → 0.60250 → 0.51622) while the
> ratio barely budges (0.7347 → 0.7297 → 0.7242). A near-field/geometry
> effect would track the geometry; a normalization error would not.
>
> **Suspected cause, NOT yet proven.** A y-normal `excite_port` is the only
> case that triggers flux-calibrated normalization
> (`src_alpha = sqrt(|get_fluxes(norm_flux)|) · exp(i·arg(norm_coeff.alpha))`,
> `fdtd.py` ~1880) instead of using the eigenmode coefficient alone. If
> meep's `|alpha|²` and its flux at that plane do not agree to a constant,
> mixing the two — flux for the source magnitude, coefficient for the
> measured response — introduces exactly a constant multiplicative error,
> which is the signature observed. `PortSpec`'s own docstring already flags
> this path as empirically derived from `ring_resonator_hybrid.py` and notes
> `sign=+1` was applied "conservatively" and never tested; this result says
> `sign=-1` (the one that script *did* use) is also not right, or at least
> not consistent with the coefficient-based convention on the other side.
>
> **Control run — done, and it localizes the fault.** The failing comparison
> changed two things at once: a y-normal excite AND `sign=-1` on both ports.
> A uniform straight (`Taper` with `width_start == width_end`,
> mirror-symmetric, both ports x-normal so flux calibration never triggers),
> forward (`sign=+1`) vs reverse (`sign=-1`), λ=1.55, same solver settings:
>
> | | \|T\| |
> |---|---|
> | forward, `sign=+1` | 0.993016 |
> | reverse, `sign=-1` | 0.991270 |
> | **ratio** | **0.998242** |
>
> `sign=-1`, the measure-side edge convention, and the forced
> `reverse=False` bypass are therefore **all clean to 0.18%** — and |T| ≈
> 0.993 on a uniform straight is itself the right answer physically
> (near-lossless), so the harness as a whole is sound. **That leaves the
> y-normal excite / flux-calibrated normalization as the only remaining
> explanation for the bend's 0.73.** Scratch probes:
> `probe_recip_symmetry.py`, `probe_sign_control.py`.
>
> **Next step to close it:** instrument the y-excite path to print
> `sqrt(|get_fluxes(norm_flux)|)` against `|norm_coeff.alpha[...]|` at the
> same plane. If their ratio is the observed ~0.73, the fix is to use one
> convention consistently rather than mixing them, and the correct choice
> has to be pinned by which one makes a *known* structure (the uniform
> straight above, driven y-normal) return |T| ≈ 1. Do NOT simply divide the
> discrepancy out — `ring_resonator_hybrid.py` adopted this calibration
> because the uncalibrated coefficient was empirically wrong for a
> backward off-axis launch, so the real question is which quantity is the
> correct source normalization, not what factor patches this one case.
>
> ### 🔬 Investigated 2026-08-01 — the localization above is WRONG, and the "validated direction" was itself 27% off
>
> Ran the doc's own prescribed next step (instrument the y-excite path, print
> `sqrt(|get_fluxes|)` against `|norm_coeff.alpha|` at the same plane) plus
> four controls. **Blocker 0 is NOT closed**, but three things are now
> established that change what to do next. Diagnostic is in the shipped code,
> env-gated: `BEEF_DEBUG_NORM=1` (pure print, no behaviour change).
>
> **Finding 1 — the two conventions disagree only on y, and only in the
> source's near field.** At the norm plane, `sqrt(flux)/|alpha|`:
>
> | | ratio |
> |---|---|
> | x-excite (validated direction) | **0.992–0.995** — agree |
> | y-excite, `monitor_margin=0.25` (shipped) | **1.516–1.580** — disagree |
> | y-excite, `monitor_margin=0.75` | **0.991** — agree again |
>
> The coefficient is what moves (5.21 → 8.07); `sqrt(flux)` is flat (8.00 →
> 8.00). So the norm monitor at 0.25 µm sits inside the y-source's near field
> where the eigenmode has not formed, and **flux is the robust quantity there,
> which is exactly why `ring_resonator_hybrid.py` adopted it empirically.**
>
> **Finding 2 — but that is NOT what breaks reciprocity.** Moving the monitor
> out of the near field changes `|T|` by 0.3% (0.737037 → 0.739323). The
> normalization convention is a red herring for the 0.73: both conventions
> leave reciprocity violated, and the bare-`alpha` convention is actively
> worse — it returns **|T| = 1.13–1.15, non-physical** (energy gain) once the
> geometry below is fixed. **Do not "fix" this by switching to bare `alpha`.**
>
> **Finding 3 — the real defect found so far is geometric, and it invalidates
> the reference direction.** `Sx = span + stub_left + stub_right + 2*dpml`, so
> with `solve_section(g, wl, 0.0, bend_radius)` the y-normal port sits at
> device x = 5.0 = **exactly the PML inner edge**, and its `width=3.0`
> aperture spans 1.5 µm into PML and 1.0 µm past the cell entirely. Giving it
> x headroom (`stub_right=2.0`, which for `RingResonatorBendSection` draws
> **no** extra geometry — `_section_polygons` emits nothing past `x=R`, so it
> is pure domain padding):
>
> | | \|T\| x→y |
> |---|---|
> | shipped test geometry | 0.711 |
> | with x headroom | **0.978** |
>
> **0.978 is the physically right answer** — a 5 µm-radius 90° bend in 450×220
> SOI is near-lossless. So test 1's 0.711, and the `>0.3` sanity range that
> blessed it, were passing on a 27%-clipped result. **The doc's control logic
> could never have caught this: the `output_margin` sweep varies *y* only, so
> it left the x-overrun constant — which is precisely why the ratio looked
> "geometry-invariant".** That invariance was evidence of an uncontrolled
> variable, not of a normalization constant.
>
> **What is still open.** With the x-overrun fixed and the monitor out of the
> near field, reciprocity is still violated at **~0.75** (0.977 vs 0.734).
> Giving the x-normal measure port y headroom (`y0=-2.0`, aperture clear of
> PML) did not move it either. So there is a real remaining asymmetry that is
> neither the flux calibration nor either aperture overrun.
>
> ## ✅✅ BLOCKER 0 IS CLOSED — 2026-08-01, jobs 11616920 + 11616921
>
> **The reciprocity "violation" was second-order discretization error the
> whole time. There is no orientation asymmetry, and no normalization bug.**
>
> Measured in 3-D on clean geometry at two resolutions:
>
> | res | \|T\| x→y | \|T\| y→x | ratio | error `1−ratio` |
> |---|---|---|---|---|
> | 15 | 1.003338 | 0.859508 | 0.855346 | 0.144654 |
> | 25 | 0.998865 | 0.946531 | **0.947608** | **0.052392** |
>
> **The convergence order is the proof, not the trend.**
>
> ```
> measured error ratio   e(15)/e(25) = 2.7610
> second-order prediction  (25/15)^2 = 2.7778     -> agree to 0.60%
> fitted order             p = 1.988              -> meep is 2nd-order
> ```
>
> A normalization error, a wrong `alpha_idx`, a polarization mismatch or a
> genuine geometric asymmetry would all give a **resolution-independent**
> ratio. This one converges as *h²*, hitting meep's own formal accuracy to
> within 0.6%. That is not consistent with anything except discretization.
>
> Richardson-extrapolating the same fit:
>
> | res | extrapolated ratio |
> |---|---|
> | 40 (`FDTD_RESOLUTION`) | **0.979** |
> | 100 | 0.997 |
>
> The x→y direction independently confirms it, converging to unity **from
> above** (1.0053 → 0.9989) exactly as a near-lossless bend should.
>
> **The y-normal excite path and its flux-calibrated normalization are
> correct.** `ring_resonator_hybrid.py` adopted that calibration empirically
> and it was right to; the earlier suspicion against it is fully withdrawn.
>
> ### The history is worth keeping, because the failure mode repeated
>
> Blocker 0 was diagnosed wrongly **three times** before this, each time by
> reasoning from a control rather than measuring the thing itself:
>
> | # | claim | why it fell |
> |---|---|---|
> | 1 | π-era mode-extraction bug | superseded before testing |
> | 2 | flux-calibrated normalization is mis-scaled | conventions agree once the norm plane leaves the source near field; switching them changes \|T\| by 0.3% |
> | 3 | y-port aperture buried in PML | real bug, real 29% error — but **not this one**; fixing it left the ratio at 0.855 |
>
> The doc asserted the ratio was *"geometry-invariant, which is what makes
> this a real finding rather than mesh noise."* It was invariant under
> `output_margin`, which varies **y only** — the one axis that could not move
> either the x-overrun or the mesh. **Resolution was never varied until
> today, and resolution was the answer.**
>
> Two of those three wrong diagnoses were nonetheless worth having: they found
> the PML-aperture defect (now guarded, and it immediately caught a fourth
> live instance) and the near-field norm-plane placement. But the lesson is
> the repo's own: *vary the thing you are claiming invariance over.*
>
> ---
>
> ### ⚠️ Job 11616611 was WITHDRAWN as contaminated, then CONFIRMED by a clean re-run (11616920)
>
> **Do not quote the table below.** The run predates the aperture guard and
> tripped exactly the defect the guard was written for: the script hardcoded
> `y0=-1.0`, which leaves the **x-normal** port's own aperture spanning
> y=[−3.75, −0.75] against a PML inner edge at ±3.25.
>
> That is not symmetric between the two directions, which is what makes it
> fatal here: in the **y→x** run the x port is the **measure** side, so the
> clip removes collected power and reads as an orientation asymmetry —
> plausibly the entire residual 14.5%. In the x→y run the same port is the
> *source*, which largely self-normalizes against its own norm monitor at the
> same aperture, so that direction looked fine (1.005) while its partner was
> depressed.
>
> The guard caught this by failing the res-25 resubmission (11616692, 9 s)
> before it burned a node — the first thing it did in anger was catch a fresh
> script written by the same person who wrote the guard.
>
> Re-running clean at res 15 (**11616920**) and res 25 (**11616921**) with
> `y0=-2.0`. ~~**Blocker 0's status is UNKNOWN until those land** — it may close.~~
>
> **✅ RESOLVED — 11616920 (17:08) reproduces the withdrawn numbers almost
> exactly, so the conclusion stands unchanged:**
>
> | λ | contaminated (11616611) | **clean (11616920)** |
> |---|---|---|
> | 1.53 | 0.849053 | 0.848557 |
> | 1.55 | 0.856094 | 0.855835 |
> | 1.57 | 0.860587 | 0.861646 |
> | **mean** | 0.855245 | **0.855346** |
>
> Agreement to 4 decimals. **The clipped x-port aperture was not material
> here**, and the physics says why: that port's aperture only overran by
> 0.5 µm at its outer *edges*, where a 0.45 µm-wide guided mode carries
> essentially no field. Contrast the y-port case, where the port's own
> **centre** sat on the PML edge — that one cost 29%. Aperture clipping
> matters in proportion to how much mode power is in the clipped region, not
> to how much geometry is.
>
> So: the withdrawal was the right procedure on a number computed against
> known-bad geometry, and re-measuring rather than arguing settled it in 17
> minutes. **Blocker 0 genuinely survives into 3-D at ~0.855, now confirmed
> on clean geometry.** ~~res 25 (11616921) is in flight to separate
> discretization from a real constant.~~ **✅ LANDED 2026-08-01 — it was
> discretization: fitted order p = 1.988, extrapolating to 0.979 at
> `FDTD_RESOLUTION = 40`. Blocker 0 CLOSED; see `## Order` item 0.**
>
> ### ❌ 3-D verdict — blocker 0 SURVIVES, but ~40% of it was a 2-D artifact
>
> (Table below is 11616611's; reproduced by the clean 11616920 above.)
>
> Every earlier probe was `three_d=False`; production is 3-D, where TE/TM no
> longer decouple. Re-ran the reciprocity pair in 3-D with **both** of the
> fixes above applied (`stub_right=2.0`, `monitor_margin=0.75`), res 15,
> `verify_ring_bend_reciprocity_3d.py`:
>
> | λ | \|T\| x→y | \|T\| y→x | ratio |
> |---|---|---|---|
> | 1.53 | 1.005730 | 0.853918 | 0.849053 |
> | 1.55 | 1.004699 | 0.860117 | 0.856094 |
> | 1.57 | 1.005504 | 0.865324 | 0.860587 |
> | | | **mean** | **0.855245** |
>
> **Two readings, and both matter.**
>
> 1. **The aperture fix is confirmed independently.** |T| x→y is now
>    **1.005** — a 5 µm-radius bend in 450×220 SOI is near-lossless, so this
>    is the physically correct answer arrived at from a second direction. The
>    0.5% overshoot is res-15 discretization. The original 0.711 was 29% low
>    purely from clipped apertures.
> 2. **The reciprocity violation is real and survives into production
>    geometry.** 2-D 0.751 → 3-D 0.855, so roughly 40% of the original gap
>    was a 2-D polarization artifact and **~14.5% is not.** The hypothesis
>    that this was purely 2-D is refuted.
>
> Mild monotone λ-dependence (0.849 → 0.861 across the band) — a constant
> normalization error would be flat, so there is some dispersive character to
> the residual too.
>
> Raw data: `example_output/verification/ring_bend_reciprocity/reciprocity.npz`.
>
> **Blocker 0 therefore still blocks the wiring**, for the third distinct
> reason in its history: not the π-era guess, not the flux calibration, not
> the aperture clipping, but a residual y-excite asymmetry that none of those
> explain.
>
> **Next steps, in order:**
> 1. **Re-baseline test 1 with `stub_right`** and tighten its assertion from
>    `>0.3` to something near 0.98. The current test cannot fail usefully.
> 2. **Audit every `PortSpec` aperture against the PML inner edge** — add a
>    guard in `solve_section()` that raises when a port's transverse extent
>    crosses PML, rather than silently clipping. This class of bug is invisible
>    and cost a 27% error that looked plausible.
> 3. Only then re-examine the residual 0.75, on geometry known to be clean.
>
> **Do not build the wiring on this yet** — the conclusion below still holds,
> for a partly different reason than it was written for.
>
> **Why this blocks the wiring below, not just item (a):** every route from
> the coupler section to a `beef.network` block needs the ring-side port's
> coefficient (`k1`, and `t_ring`/`k2` from the ring-excited run). A
> constant ~0.73 on those relative to the bus-side `t` would leave the
> closed-form `E_thru = t + k1·k2·P/(1−t_ring·P)` with mis-scaled coupling —
> corrupting resonance depth and Q while still producing a plausible-looking
> curve. Gap 4 below (no full-device reference) is exactly why that could
> have gone unnoticed. **This same suspicion applies to
> `ring_resonator_hybrid.py`'s existing `k2`/`t_ring`**, which use the
> hand-rolled equivalent of the same calibration — so it is not only a
> new-code concern.

**(b) A `RingResonator` DeviceGeometry — done, 2026-07-29.**
`RingResonatorCouplerSection` (bus + one near-bend + ring stub, `side="left"`
or `"right"`) and `RingResonatorBendSection` (the standalone 90° bend,
x-normal in / y-normal out) both implement `section_prisms`/`section_gds`/
`width_at`/`center_at`, in `beef/geometry/ring_coupler.py` and `ring_bend.py`.
Breakpoints verified to match the reference scripts' own constants exactly
(`stub_x=-7.0`, `hc=-2.0`, `iface_x=-1.0`, `ring_port().y_boundary=6.20` for
`side="left"`; left/right confirmed exact x-mirrors). Since `width_at`/
`center_at`'s x-in/y-out contract can't honestly describe a fixed-x slice
through the curved ring stub (not a function of x there, the same reason
`SBend`/`DirectionalCoupler` forbid an interior bend cut at all),
`ring_port()`/`out_port()` are the documented escape hatch, returning the
perpendicular port's own geometry as plain floats for a caller to build a
`PortSpec` + `y0`/`y1` from. 46 new fast, meep-free geometry tests, following
`tests/test_geometry.py`'s existing style.

The body (straight bus+ring pair) gets no new geometry code — it stays a
raw-polygon `meow` call, matching how `solve_body_supermodes()` already
works in the reference scripts.

**(c) `HybridSolver` must express feedback.** This was the substantive one. The
solver still contracts a straight cascade — threading light through the
sections once, not a loop. Two call sites, both unchanged in substance since
this item opened, only their line numbers have drifted again (checked
2026-07-30): `solve()`'s `np.einsum("mij,mjk->mik", mat, full)` at
`hybrid.py:874`, and `solve_length_sweep()`'s batched counterpart —
`np.einsum("Lmij,Lmjk->Lmik", mat, full)`, an extra leading length axis, same
cascade — at `hybrid.py:973`. A ring is a **closed loop**: light re-enters
the coupling region every round trip, which is a closed-form all-pass sum,
not a matrix product. Two ways out were scoped:

- *Narrow:* keep the cascade for the open path (coupler section + bends +
  straights) and add an explicit loop-closure step that applies the resonator sum
  to a designated feedback port pair. Less work, and it matches how the script
  already assembles the answer — but it hard-codes "ring" into the solver.
- *General:* let `HybridSolver` describe an S-matrix **network** with named ports
  rather than an ordered chain, and solve it by Redheffer star products /
  node elimination. Cycles then fall out for free, and it would also cover
  Mach–Zehnders, cascaded rings, and any device where a section feeds back. This
  is the architecturally right answer and considerably more work.

**The general form's composition primitive is done, 2026-07-29 — `beef.network`
(`connect()`/`close_loop()`/`Circuit`, see `notes.md`'s Network section).**
`close_loop()` reproduces this device's own proven closed-form formula
(`E_thru = t + k1·k2·P/(1-t_ring·P)`, `ring_resonator_hybrid.py:762-774`) term
for term when applied to a 2-port `[bus, ring]` block — the concrete evidence
the general primitive is physically correct on exactly this device's math, not
just internally consistent. **Not yet wired to this example or to
`HybridSolver`**: the module depends only on `SMatrix`/`numpy` and has no
caller in `beef/solvers/` or any example — (a) and (b) above are now also
done, so nothing structural blocks the wiring anymore, but doing it (and
rewiring `HybridSolver`'s own cascade onto `beef.network`) is separate,
undone work. `close_loop()`
also only implements the same-port self-loop case (this device's case); the
two-different-ports case some `Circuit`-discovered cycles would need (e.g.
certain cascaded-ring topologies) is deferred.

**All three of (a), (b) and (c) now exist as standalone, tested library
capabilities — none are wired to each other or to this device yet.** That
wiring is the actual remaining work: build a real `RingResonatorCouplerSection`
+ `EMESolver` body + mirrored right section via `HybridSolver`/`solve_section()`
(exercising the still-unverified y-normal excite/flux-calibration path),
compose with `beef.network.close_loop()`, and check the result against this
device's own proven closed-form/full-device numbers. A new example script
would naturally follow from that, replacing the current hand-rolled scripts.

### Wiring to-do — itemized 2026-07-31

Scoped by reading the three modules' actual signatures (not just this doc's
prose). In dependency order:

0. ~~**⚠️ Fix the y-normal excite normalization**~~ — **✅ CLOSED 2026-08-01
   (11616920/11616921). There was nothing to fix.** The apparent reciprocity
   violation is second-order discretization: in 3-D on clean geometry the
   ratio converges as h^2 (0.8553 at res 15 -> 0.9476 at res 25, fitted order
   **p = 1.988** against meep's formal 2nd order), extrapolating to **0.979 at
   `FDTD_RESOLUTION = 40`**. A normalization or orientation bug would be
   resolution-independent. `k1`/`k2`/`t_ring` can be trusted from this path.
   **The wiring below is UNBLOCKED.** Full evidence in the boxed section under
   prerequisite (a).

1. **Decide the port topology.** `RingResonatorCouplerSection` is
   **per-side** (`side="left"` or `"right"`), but `ring_resonator_hybrid.py`
   simulates *both* sides in one meep run (its bus spans both stub columns,
   `±INNER_X_HALF`, with ring stubs at `±STUB_X`). So the sectioned form is
   naturally a **4-port** network — `bus_left`, `bus_right`, `ring_left`,
   `ring_right` — not the 2-port `[bus, ring]` block `close_loop()`'s
   docstring worked through. Two options, unchosen: compose two per-side
   2-port-ish blocks with `connect()` across the bus, or assemble the pair
   as a single 4-port block and close two loops. **This is a real design
   decision and the first thing to settle after item 0.**

   > ### 🔬 Scoped 2026-08-01 — the physical 4-port route is BLOCKED on `close_loop()` Branch B
   >
   > Read against the three modules' actual signatures. **The two options are
   > not equal-cost, and the difference was not previously identified.**
   >
   > **Port count, from the geometry rather than from prose.** Each
   > `RingResonatorCouplerSection` exposes **three** ports: `outer_x` (the
   > true physical bus port), `iface_x` (the interior interface cut past the
   > near bend), and `ring_port()` (the y-normal stub). Composing the pair
   > with `connect(left, "iface", right, "iface")` leaves exactly the
   > 4-port set this item predicted: `bus_L`, `ring_L`, `bus_R`, `ring_R`.
   >
   > **The blocker.** The physical ring loop runs `ring_R` → around the outer
   > ring → `ring_L`. That is `port_out != port_in`, and
   > `close_loop()` **raises `NotImplementedError`** for exactly that case.
   > `ops.py`'s own module docstring is explicit that this is *not* a
   > degenerate case of the implemented one: *"Branch B starts from two
   > independent free variables (a_p, a_q) and removes only one … it's two
   > different formulas."* It records **"Branch B has no required caller in
   > this pass (deferred)."**
   >
   > **This device is that caller.** The deferral and this to-do were written
   > independently and never joined up.
   >
   > **Why the legacy form escapes it.** `E_thru = t + k1·k2·P/(1−t_ring·P)`
   > treats `bus` and `ring` as **aggregate** ports (`S_bus,bus = t` is a
   > through-path, not a reflection), which makes the loop a same-port
   > self-loop — Branch A, implemented, and already cross-checked term for
   > term by `tests/test_network_ring_crosscheck.py`. The cost is that those
   > ports are not physical ends, so they do **not** compose with
   > `connect()`, and the shape hard-codes single-coupler racetrack topology
   > (no add-drop, no cascade).
   >
   > | | Branch A / 2-port aggregate | Branch B / physical 4-port |
   > |---|---|---|
   > | works today | **yes** | no — `NotImplementedError` |
   > | reproduces legacy formula | exactly, already tested | yes, but needs new algebra |
   > | composes via `connect()` | no | yes |
   > | extends to add-drop / cascaded rings | no | yes |
   >
   > **Recommendation — do BOTH, in this order, and the order is the point.**
   >
   > 1. **Build the 2-port aggregate form first**, and use it purely as the
   >    numerical cross-check item 4 already calls for. It needs no new
   >    library work, and because `close_loop()`'s algebra is *already*
   >    proven against the legacy formula, **any disagreement is necessarily
   >    in items 2/3** — the hand-assembled `SMatrixBlock` or the propagator
   >    `P`. That is the whole diagnostic value: it isolates the inputs.
   > 2. **Then implement Branch B and move to the physical 4-port**, which is
   >    the form that composes and the only one that survives contact with
   >    add-drop or cascaded topologies.
   >
   > Doing 2 first means debugging new solver inputs *and* new network
   > algebra simultaneously, with no reference — the same mistake pattern
   > this repo has now paid for twice (the π-sign-flip misdiagnosis, and
   > blocker 0's own localization).
   >
   > **Schedule Branch B explicitly.** It is a real work item with no owner
   > today, it gates the general route, and its absence is currently
   > invisible — nothing fails until someone calls `close_loop()` with two
   > different ports.
   >
   > ### ✅ DECIDED 2026-08-01 — 2-port aggregate first, then Branch B + 4-port
   >
   > The recommended sequencing is adopted. **Both forms get built; the order
   > is the deliverable.**
   >
   > 1. **2-port `[bus, ring]` aggregate, Branch A.** No new library work.
   >    Its only job is wiring item 4's cross-check — and because
   >    `close_loop()`'s algebra is *already* proven against the legacy
   >    formula by `tests/test_network_ring_crosscheck.py`, **any disagreement
   >    is necessarily in items 2/3** (the hand-assembled `SMatrixBlock` or
   >    the propagator `P`). One unknown at a time.
   > 2. **Then `close_loop()` Branch B + the physical 4-port**, which is the
   >    form that composes via `connect()` and the only one that survives
   >    add-drop or cascaded topologies — `ring_resonator_hybrid_cascaded.py`
   >    has no path forward without it.
   >
   > **Branch B is hereby a scheduled work item**, blocking step 2 only. It
   > should be tracked in `beef/network/ops.py`'s own docstring (which
   > currently says "no required caller in this pass") so the deferral and its
   > caller stop being recorded in two places that never reference each other.

> ### 🔬 Item 2 in progress, 2026-08-01 — three defects found, two fixed, one open
>
> `examples/verification/python/verify_ring_network_crosscheck.py` builds the
> 2-port aggregate block from `RingResonatorCouplerSection` through
> `FDTDSolver.solve_section()` and runs both routes. Status after three 2-D
> res-20 iterations (~2:40 each):
>
> | | \|t\|²+\|k1\|² (→1) | max \|\|k1\|−\|k2\|\| (→0) |
> |---|---|---|
> | first run | 0.008 – 0.198 | 0.380 |
> | + port-role fix | 0.008 – 0.198 | 0.380 |
> | + half-domain bus aperture | **0.049 – 0.326** | **0.097** ✅ |
>
> **Defect 1 — the block index convention (FIXED, no compute).** Verified
> `data[:,0,0]=t, [0,1]=k2, [1,0]=k1, [1,1]=t_ring` against the closed form
> with synthetic coefficients: agreement **5.6e-17**. An `S[out,in]` vs
> `S[in,out]` slip would have silently swapped `k1`/`k2` and still produced a
> plausible resonance curve.
>
> **Defect 2 — excite/measure edge convention (FIXED).** A single `PortSpec`
> cannot serve both roles: `ring_port().sign = -1` is the *excite* convention,
> and reusing it on the measure side puts the ring monitor at the LOW y edge,
> where there is no ring. `t_ring` went **0.006 → 0.54** on this fix alone.
>
> **Defect 3 — the bus aperture contained the ring (FIXED).** Bus centreline
> y=0, ring y=0.70; a `PORT_APERTURE=3.0` window spans y=[−1.5,+1.5] and so
> launches a two-guide **supermode**, not the bus mode. Note
> **`half_domain_aperture` does NOT reach a `PortSpec`** — it is consumed only
> by `section_fdtd_layout()`, which `PortSpec` bypasses by design, so setting
> it on the solver is silently ignored. The split is computed explicitly
> instead: y=[−1.750, 0.350], which excludes the ring and clears PML
> (reach 3.850 vs inner 4.100). **This fixed reciprocity** (0.380 → 0.097).
>
> **⚠ STILL OPEN — power conservation.** `|t|²+|k1|²` is 0.049–0.326 against
> ~1. Improved 4× but nowhere near closed: roughly 70% of the excited power
> is leaving through something this 2-port model does not measure, and `k1`
> stays tiny (~0.01) so it is not going into the ring either.
>
> ### 🛑 STOPPED ITERATING 2026-08-01 — a control exonerates everything generic
>
> Seven successive hypotheses were tried and **all seven were wrong**: port
> role, bus aperture width, "move the cut", two-section composition, seam mode
> count, seam aperture width, aperture offset. Each was inferred from the
> composed result on a complicated geometry. The pattern, not any one guess,
> is the finding.
>
> **The control that should have come first**
> (`examples/verification/python/verify_solve_section_power_control.py`, job
> 11620125): a uniform straight waveguide, no ring/bend/seam/coupling, at the
> **same solver settings** the ring cross-check uses. `|T|²` must be 1.
>
> | case | \|T\|² |
> |---|---|
> | layout path (no PortSpec) | 0.8874 |
> | PortSpec, aperture 1.0 / 2.0 / 3.0 | 0.99988 / 0.99983 / 0.99982 |
> | PortSpec + tall y=[−2.0, 6.2] | 0.99982 |
> | PortSpec, aperture offset −0.35 | 0.99980 |
> | PortSpec, aperture offset −0.70 | **0.99861** |
>
> **`solve_section` + `PortSpec` conserve power to ≤0.14% across every knob
> the ring work spent seven rounds adjusting.** Solver, resolution,
> `decay_tol`, `PortSpec` itself, aperture width, the tall forced y-domain and
> the half-domain offset are **all exonerated**. The loss is specific to the
> ring coupler geometry.
>
> **Incidental finding worth its own item:** the *layout* path loses **11%**
> on the same trivial guide where `PortSpec` loses 0.02%. That is the path
> every other device uses. It is not this device's problem and was not chased
> here, but it should be.
>
> **The concrete next step, and it is a measurement not a guess.** The left
> section is a genuine **3-port** (bus, ring, iface), and the audit only ever
> measured two paths out of it — `iface←bus` and `iface←ring`. It never
> measured **`ring←bus`**. Power coupling from the bus up the ring stub and
> out through the top boundary is a legitimate third exit that the audit
> counts as loss. So `capture(iface←bus)` alone was never required to be 1;
> the quantity that must sum to ~1 is `capture(iface←bus) + capture(ring←bus)`.
> Measure that before forming any further hypothesis.
>
> ### ❌ MEASURED 2026-08-02 (job 11627160, 7:04) — the third exit is NOT there
>
> `ring<-bus` was added to the left-section audit as a diagnostic (not folded
> into `data`, so the composition under test is unchanged). It is **zero**:
>
> ```
> left iface<-bus:  0.1422 0.2458 0.1908 0.0668 0.0367
> left iface<-ring: 0.5813 0.5643 0.4914 0.5756 0.4167
> left ring<-bus:   0.0004 0.0002 0.0001 0.0001 0.0000   <-- the hypothesis
> left TOTAL:       0.1426 0.2460 0.1909 0.0669 0.0367
> ```
>
> Adding the ring stub moves the total by **0.0004**. Hypothesis eight joins
> the other seven. **It was still the right thing to run** — it was cheap, it
> was a measurement rather than a guess, and it eliminates the whole
> "legitimate third exit" class rather than one variant of it.
>
> **What the run did establish, from all four seam paths:**
>
> | path | captured |
> |---|---|
> | `iface<-bus` | 0.037 – 0.246 |
> | `bus<-iface` | 0.033 – 0.216 |
> | `ring<-iface` | 0.229 – 0.253 |
> | `iface<-ring` | 0.417 – 0.581 |
>
> The bus↔seam pair agrees in both directions (0.033–0.246 either way), so the
> bus path is **reciprocal and genuinely low** — this is not an edge/sign slip
> on one side. Every path through the seam is lossy; only `iface<-ring` is
> even half-healthy. Composed, `|t|^2+|k1|^2` = **0.0022 – 0.0656**, which is
> about what multiplying two ~0.15 sections predicts.
>
> ⚠ **That composed range is worse than the 0.049–0.326 in the table above.**
> That row was measured before the narrow seam aperture landed; the committed
> configuration is not the one that table records. Do not read the two as a
> trend.
>
> **The routes still agreed to 1.388e-17.** Third run in a row where the
> cross-check passes at machine precision while the inputs are ~85% wrong —
> the method note below is now confirmed three times over.
>
> **Next measurement, and it is a control, not a hypothesis: `bus<-bus`.**
> Excite the bus at `outer_x`, measure at `iface_x` with the *bus* aperture
> instead of the seam aperture — one more `solve_section` at the same cost.
> `bus_measure` already exists for it. This splits the remaining space in two:
>
> - **`bus<-bus` ≈ 1** → the power arrives at the seam plane fine and the
>   seam's 4 supermodes fail to represent it. A **basis** problem: fix the
>   seam window and mode count, no geometry change.
> - **`bus<-bus` ≈ 0.15** → the section really loses it before the seam.
>   Radiation or reflection — and note this script models forward
>   transmission only, so **reflection is currently not measured at all** and
>   is an unexamined exit.
>
> A standing hypothesis, recorded but *not* acted on: the bus excite window is
> y=[−1.750, **0.350**] — 1.75 µm below the guide at y=0 but only 0.35 µm
> above it, to clear the ring at y=0.70. A 0.45 µm strip mode's evanescent
> tail runs past 0.35 µm, so the launch may be a clipped bus mode. Job
> 11620125 does **not** exonerate this: its offset-aperture cases used the
> same aperture at both ends, where a clipped mode still transmits perfectly
> into itself. Here the excite and measure windows are different shapes.
> **Run `bus<-bus` before touching this.**
>
> ### ✅ RAN 2026-08-02 (job 11627401, 8:21) — not reflection, not the basis
>
> ```
> T bus<-bus (seam plane, bus basis): 0.2757 0.3258 0.2120 0.0486 0.0681
> R bus<-bus (reflection at outer_x): 0.0175 0.0020 0.0026 0.0206 0.0012
> BUDGET fwd + refl + ring:           0.2935 0.3279 0.2147 0.0693 0.0694
> unaccounted:                       +0.7065 +0.6721 +0.7853 +0.9307 +0.9306
> ```
>
> **Reflection is 0.1–2%.** Ruled out — and with it the `iface_x` cut as the
> loss mechanism, since a cut inside the coupled region would reflect.
>
> **A second, independent basis at the same plane recovers the same order.**
> The bus-shaped 1-mode window gives 0.049–0.326 against the 4-supermode seam
> window's 0.037–0.246. Two different bases agreeing kills the "the seam's
> supermodes can't represent it" reading. ~80% is unaccounted.
>
> **Two readings survive, and they need different fixes.** Either the bend
> genuinely radiates, or the source launches a clipped non-eigenmode field
> that sheds immediately — the bus window is y=[−1.750, **0.350**], only
> 0.35 µm above a centreline at y=0. Control queued (job 11627585): the same
> launch and the same measure window over `[outer_x, −5.0]`, a straight bus
> with no bend and no coupling. ~1 means the launch is sound; ~0.2 means every
> bus-excited number since defect 3 inherits an aperture defect.
>
> ### ⛔ NOT MEASURED 2026-08-02 (jobs 11627585, 11628207) — the control never ran
>
> **Both attempts died in geometry construction, before any FDTD physics.**
> Neither produced a number. The question above is still open, and the two
> readings that need different fixes are both still alive.
>
> ```
> 11627585  X_STRAIGHT=-5.0  ValueError: Cut/stub at x=-3.5000 ... strictly inside the near-bend arc
> 11628207  X_STRAIGHT=-8.0  ValueError: Cut/stub at x=-6.5000 ... strictly inside the near-bend arc
> ```
>
> **The mechanism, measured directly off the geometry class rather than read
> off the logs.** `_check_cut` (`beef/geometry/ring_coupler.py:168`) rejects
> cuts strictly inside the arc, which here spans `(stub_x, hc) = (-7.0,
> -2.0)`. The cut it validates is **not** `x1 + stub_right`:
> `FDTDSolver.solve_section()` (`beef/solvers/fdtd.py:1978`) pads both stubs
> by `dpml` before handing them to `section_prisms()`, so the tested cut is
> `X_STRAIGHT + stub_right + dpml = X_STRAIGHT + 1.5`. Evaluating that against
> the class reproduces **both** jobs' error strings exactly (−5.0 → −3.50;
> −8.0 → −6.50), which is what makes this the mechanism and not a story about
> it. Legal requires `X_STRAIGHT ≤ −8.5`.
>
> ⚠ **The committed script is still broken.** `X_STRAIGHT` still defaults to
> −8.0 and the explanatory comment at
> `verify_ring_network_crosscheck.py:56-67` still asserts the `+ dpml`-free
> arithmetic that produced the second failure. Resubmitting as-is fails a
> third time in ~3 minutes. **Fix the constant and the comment together** —
> the comment is what will otherwise re-derive the same wrong value.
>
> **The design call blocking a one-line fix.** −8.5 is the largest legal
> value and leaves only **0.5 µm** of straight bus. The control's whole
> premise is that a clipped launch sheds *over a propagation distance*, so a
> 0.5 µm span may not answer it. The alternative — `stub_right=0.0` for this
> call alone, allowing `X_STRAIGHT ≤ −7.5` and a 1.5 µm span — buys length by
> changing the port padding away from what every other measurement in this
> script used, which is its own confound. Choose deliberately; this is the
> same "sweep the wrong thing and hide the real error" shape as coupler item
> 11 and MMI gap 2.
>
> **What this costs.** Two ~3:20 failures on 24 cores. Small, but avoidable:
> both are pure geometry and reproduce on a login node in under a second with
> no allocation at all. **Cut-legality for this device is a login-node check,
> not a batch job** — run it before submitting anything that moves a section
> boundary.
>
> ### 🛠 DEFECT FOUND AND FIXED 2026-08-02 — the section padding drew no waveguide
>
> **This displaces both surviving readings above. Neither was tested against
> a third possibility that turns out to be a real library bug.**
>
> `RingResonatorCouplerSection._section_polygons()` clamped the bus rect to
> the class's own `[outer_x, iface_x]` range:
>
> ```python
> _rect(max(lo, outer_lo), min(hi, outer_hi), 0.0, ww)      # before
> _rect(lo, hi, 0.0, ww)                                    # after
> ```
>
> So the caller's `stub_left`/`stub_right` padding drew **nothing** — the
> padded domain past `outer_x` and past `iface_x` was bare cladding.
>
> **Measured, on the login node, before the fix** (`_section_polygons(x0=-9.0,
> x1=-1.0, stub_left=stub_right=1.5)`, i.e. the script's 1.0 plus
> `solve_section()`'s own `dpml`):
>
> | | padded domain | bus actually drawn |
> |---|---|---|
> | before | `[-10.500, 0.500]` | `[-9.000, -1.000]` |
> | after | `[-10.500, 0.500]` | `[-10.500, 0.500]` |
>
> **Where the source and monitor sit is what makes this fatal.**
> `solve_section()` places them at `lo + source_margin` and
> `hi - monitor_margin` — off the **padded** boundary, not off `x0`/`x1`
> (`beef/solvers/fdtd.py:2055,2092`). With this script's `source_margin=0.0`,
> `monitor_margin=0.75`, `stub_left=stub_right=1.0`:
>
> - source plane **x = −10.0** — bus started at −9.0, so **1.0 µm out in bare
>   cladding, and exactly on the PML inner edge**;
> - seam monitor plane **x = −0.75** — geometry ended at −1.0, so **0.25 µm
>   past the end of the waveguide**, also cladding.
>
> An eigenmode source resolves the eigenmode of the cross-section *at its own
> plane*. On bare cladding there is no guided mode to launch. This is the same
> failure class as blocker 0's y-port-in-PML (which cost 27%), and the
> `_check_aperture` guard at `fdtd.py:2103` does not catch it — that guard
> checks the *transverse* aperture against the PML, never the along-axis
> position against drawn material.
>
> **This class was the outlier, which is the strongest evidence it was
> wrong.** `DirectionalCoupler._section_polygons()` draws its access arms out
> to `lo`/`hi` (`beef/geometry/coupler.py:221,238`) and
> `RingResonatorBendSection` draws its input straight from `lo`
> (`ring_bend.py:110`). The ring coupler's own class docstring already
> promised the DC convention — *"outer_x is the true physical bus port, where
> an external caller's own stub_left/waveguide_length padding would be
> added"* — and the code contradicted it.
>
> **Why it survived.** No test exercised this class with non-zero stubs:
> `test_section_polygons_full_span` passed `stub_left=stub_right=0`, where the
> clamped and unclamped forms are identical. Regression test added
> (`tests/test_geometry.py::test_stub_padding_draws_waveguide_not_bare_cladding`,
> both sides): the bus must span the full padded domain, and the sliver must
> follow the padding out through the seam end while still stopping at `hc`.
> Full suite green after the change — **527 passed, 18 skipped** (13.7 s),
> including all three ring test files.
>
> **What the fix does NOT do.** The `hc` end of the sliver still clamps —
> past it is the bend, not more straight guide — and the arc/stub column are
> untouched. **No cache is stranded:** `verify_ring_network_crosscheck.py`
> passes no `save_dir`, and it is the only non-test caller of this class, so
> the polygon fingerprint change reaches nothing on disk.
>
> **Status of the numbers above: every bus-excited result in item 2 is void,
> not just the control.** `T bus<-bus`, `R bus<-bus`, the BUDGET row and the
> 11627160 `ring<-bus` result all launched from x = −10.0, on bare cladding.
>
> ### ✅ MEASURED 2026-08-02 — controlled A/B, and the clamp IS the missing 80%
>
> Two login-node runs of `verify_ring_network_crosscheck.py`, **2-D, res 8,
> 1 wavelength, N_IFACE=2**, differing in exactly one thing: whether
> `_section_polygons()` clamps. The pre-fix form was monkeypatched back in
> (`_section_polygons_PREFIX`, with an assertion that the clamp really was
> restored before the run was spent) so that resolution, dimensionality,
> aperture, wavelength count and every solver setting are held.
>
> | left section | pre-fix clamp | post-fix |
> |---|---|---|
> | `iface<-bus` | 0.2065 | **1.0039** |
> | `iface<-ring` | 0.3074 | **0.6449** |
> | `ring<-bus` | 0.0125 | 0.0048 |
> | **TOTAL out of bus** | **0.2190** | **1.0087** |
> | deficit vs 1.0 | **+0.7810** | **−0.0087** |
> | `T bus<-bus` (seam plane, bus basis) | 0.1059 | **0.9941** |
> | BUDGET fwd + refl + ring | 0.1216 | 1.1735 |
> | unaccounted | **+0.8784** | −0.1735 |
> | `CONTROL bus<-bus` straight, 0.50 µm | **0.1353** | **1.0580** |
>
> **+0.7810 reproduces the production +0.7065…+0.9307 range** at a completely
> different resolution and dimensionality. That is what makes this the cause
> and not a correlate: the defect's signature survives the config change, and
> removing the defect removes the signature.
>
> Post-fix the script also reports `T bus<-bus` (seam plane, bus basis) =
> **0.9941** against 0.049–0.326 before, and `physics checks PASS (power sum
> ~1, reciprocity holds)` — the assertions that had been the only thing
> catching defects here (the route comparison never caught any of them).
>
> **The launch control answered itself: 1.0580.** Post-fix, `CONTROL bus<-bus
> straight` over the derived 0.50 µm span comes back ≈1. The question it was
> written to settle — "~1 means the launch is sound, ~0.2 means the aperture
> clips" — landed unambiguously on the ~1 side, so **the half-domain aperture
> is exonerated and the 0.5-µm-span design call is moot.** The span knob
> (`bus_clearance`) is documented in the script for anyone who wants a longer
> control later; it is no longer blocking anything.
>
> ### ⚠ Methodological warning — the control would have MISLED us
>
> Run under the pre-fix clamp, the same control returns **0.1353**. The
> hypothesis space written on 2026-08-02 read that as *"~0.2 → the
> half-domain aperture clips the bus mode"*. **That conclusion would have
> been wrong**, drawn from a correct measurement: the aperture was fine, and
> what the control was actually detecting was a source plane sitting on bare
> cladding. Had jobs 11627585/11628207 succeeded, this device would very
> likely have spent the next day re-shaping an aperture that was never the
> problem. The two geometry crashes, which cost ~6:40 of queue time and read
> as pure waste at the time, bought that day back.
>
> The transferable lesson is the one already on record and violated anyway:
> **an enumerated hypothesis space is not a closed one.** Both surviving
> readings ("the bend radiates" / "the aperture clips") were about the
> *physics being modelled*; the actual cause was in the *model being built*,
> a category neither branch covered. Before running a discriminator, ask what
> it returns if the geometry itself is wrong — here, both branches predicted
> the same ~0.2 the real defect produced.
>
> **Note also what did NOT catch it.** The pre-fix run still printed
> `VERDICT: routes AGREE to machine precision`, with ~88% of the power
> missing. This is the fourth independent confirmation of the standing entry
> below: both routes consume the same inputs, so the route comparison agrees
> regardless of how wrong those inputs are. Every defect on this device has
> been caught by a physics assertion (power sum, reciprocity) or by reading
> the geometry — never once by the cross-check the script is named after.
>
> Values slightly above 1.0 (1.0087, 1.0580) are res-8 normalisation error,
> not gain — they are the reason this pair is a *diagnostic* A/B and the 3-D
> res-20 production re-run was still owed.
>
> ### ✅✅ CONFIRMED AT PRODUCTION CONFIG — job 11634859, 2026-08-02
>
> **The A/B is no longer the only evidence.** Re-run at the exact config the
> void numbers were taken at — **3-D, res 20, 5 wavelengths, N_IFACE=4**:
>
> ```
> left iface<-bus:  1.0056 1.0073 1.0084 1.0098 1.0117
> left iface<-ring: 0.9903 0.9906 0.9907 0.9916 0.9920
> left ring<-bus:   0.0000 0.0000 0.0000 0.0000 0.0000
> left TOTAL out of bus = 1.0056 1.0073 1.0084 1.0098 1.0117
> left deficit vs 1.0   = -0.0056 -0.0073 -0.0084 -0.0098 -0.0117
> ```
>
> | | 11627401 (pre-fix) | **11634859 (post-fix)** |
> |---|---|---|
> | `iface<-bus` | 0.037–0.246 | **1.0056–1.0117** |
> | deficit vs 1.0 | **+0.7065 … +0.9307** | **−0.0056 … −0.0117** |
>
> **Item 2's central question is closed.** The ~80% was never radiation, never
> a basis deficiency, and never the aperture — it was
> `_section_polygons()` drawing bare cladding where the caller's stub padding
> should have continued the bus. All seven-plus hypotheses raised on
> 2026-08-01 and 2026-08-02 were about the physics; the defect was in the
> geometry construction.
>
> `ring<-bus = 0.0000` reproduces 11627160 exactly, so **hypothesis eight
> stays refuted** — that result was correct all along and is unaffected by the
> clamp (the ring stub is drawn by the all-or-nothing bend branch, which never
> clamped).
>
> **The left section's full audit closed at 21:16, and the last two lines are
> the ones that end the argument:**
>
> ```
> left unaccounted (-> radiation) = +0.0040 +0.0064 +0.0090 +0.0098 +0.0117
> left CONTROL bus<-bus straight [-9.00,-8.50] (0.50 um, no bend)
>                                 = 1.0020 1.0029 1.0038 1.0047 1.0056
> ```
>
> | | 11627401 (pre-fix) | **11634859 (post-fix)** |
> |---|---|---|
> | unaccounted → radiation | **+0.7065 … +0.9307** | **+0.0040 … +0.0117** |
> | `CONTROL` straight bus | *never ran* | **1.0020 … 1.0056** |
>
> **The budget closes to within ~1%.** "The bend genuinely radiates" — the
> first of the two surviving readings — is refuted: there is no longer any
> meaningful radiation to attribute to it.
>
> **The launch control finally ran, at production config, and returns ≈1.**
> Its own printed rubric reads *"~1 → launch sound; ~0.2 → the half-domain
> aperture clips the bus mode."* It landed at 1.002–1.006. The second
> surviving reading — the aperture — is refuted too. **Both branches of the
> 2026-08-02 hypothesis space are now closed by measurement, and neither was
> the cause.**
>
> Note what this control cost to obtain: two failed jobs, a derived-cut fix,
> and a library defect found on the way. Had it run on the first attempt it
> would have returned **0.1353** (see the A/B) and pointed at the aperture —
> the wrong answer. **Runtime corroborates:** 1:01+ against pre-fix runs of 7:04
> and 8:21 — pre-fix, the source radiated into the PML and
> `stop_when_dft_decayed` tripped early, i.e. it ran fast *because* the power
> was vanishing. Suggestive, not proof.
>
> ⚠ **Cost of this measurement, recorded because it was avoidable.** The
> script hardcodes `OUT` and `np.savez`es to it unconditionally, so both smoke
> runs overwrote
> `example_output/verification/ring_network_crosscheck/crosscheck.npz` —
> destroying the array from the last successful production run (almost
> certainly 11627401). Those numbers were pre-clamp-fix and therefore void
> anyway, and the printed diagnostics survive in `reports/`, but the `.npz`
> was the artifact this doc's own standard says to quote. **Redirect `OUT`
> before running this script outside the queue.**
>
> ✅ **RESOLVED 2026-08-02T23:10 — the `.npz` is back.** 11634859 ran to
> completion and re-wrote `crosscheck.npz` on its way out. The restored array is
> post-clamp-fix, i.e. strictly better than what was lost. The redirect warning
> above still applies to future login-node runs; the data loss itself cost
> nothing in the end.
>
> ### ✅ Right section and composition — job 11634859 completed, 3:56:16
>
> **The job finished** (exit `0:0`, ended 2026-08-02T23:10:59, well inside
> `-t 480`). The block above was written from its left-section output while it
> was still at 1:57:04. The remainder:
>
> ```
> right bus<-iface:  summed |T|^2 over seam modes = 0.9997 0.9992 0.9970 0.9953 0.9933
> right ring<-iface: summed |T|^2 over seam modes = 0.9989 0.9985 0.9989 0.9986 0.9980
> composed ports: ('bus_L', 'ring_L', 'bus_R', 'ring_R')
> bus-side power sum |t|^2+|k1|^2: min 1.0050 max 1.0067   (lossless coupler -> ~1)
> reciprocity max ||k1|-|k2||: 0.0021   (-> 0)
> physics checks PASS (power sum ~1, reciprocity holds)
> worst |E_legacy - E_network| = 1.735e-18
> L_outer = 7.0000 um
> VERDICT: routes AGREE to machine precision
> ```
>
> **The fix holds on both sides.** The right section shows no clamp pathology
> (0.993–1.000 across the band); it was never expected to, since the defect was
> in a range clamp that both sides share, but this is the measurement rather
> than the inference. Power sum and reciprocity now pass on the **composed
> 4-port**, not merely per-section.
>
> ⚠ **Weight these three lines correctly.** `VERDICT: routes AGREE to machine
> precision` and `worst |E_legacy - E_network| = 1.735e-18` are the *weakest*
> lines in the block — the pre-fix run printed the same verdict with ~88% of the
> power missing (fourth confirmation, standing entry below). The load-bearing
> lines are the power sum and the per-section budgets, which are physics
> assertions and did move when the defect was fixed.
>
> 📌 **`L_outer = 7.0000 um` corrects gap 3.** This is the first time the true
> lever arm has been printed by a running job. Gap 3's convergence script
> assumed **14.0 µm**. See gap 3 below — its Δφ column is 2× too large.
>
> ### ⚠ STRUCTURAL LIMIT — this audit can only ever see modal content
>
> **New item, 2026-08-02.** Every number in item 2, in all four runs, is a
> **modal projection coefficient** from `get_eigenmode_coefficients()`. None
> is a basis-free power measurement. That is why the audit cannot presently
> distinguish "the power radiated" from "the power is there and no basis we
> tried represents it" — it can only accumulate bases and watch them agree,
> which is evidence but never proof.
>
> `solve_section()` returns modal coefficients only. The machinery for the
> basis-free answer already exists inside it — `mp.get_fluxes()` is called on
> the y-excite calibration path (`beef/solvers/fdtd.py:2176,2187`) — but there
> is **no caller-facing way to read raw Poynting flux through a plane**.
>
> **The work item:** expose total flux from `solve_section()` (a debug return,
> an optional out-parameter, or a sibling method — the shape is a design call,
> the gap is not). Then the power budget closes against a number that owes
> nothing to a mode basis, and questions of this kind stop needing three jobs
> and a process of elimination.
>
> This is not ring-specific. Any device whose section audit fails to conserve
> power hits the same wall; the ring is simply the first to fail badly enough
> to expose it. **Applies to the 11% layout-path loss too** — that finding is
> also a modal number.
>
> ~~**Leading hypothesis — the section cut.**~~ **Refuted** — see above; the
> cut is legal and the composition is correct. The section spans
> `outer_x = -9.0` → `iface_x = -1.0`, and `iface_x` sits beside the bend
> centre `hc = -2.0`, so the measure plane may still be inside the coupled
> region. That is the doc's own "move the cut" option and it is a **design
> decision, not a bug fix** — `iface_x` is fixed by the geometry class's
> breakpoints, chosen to match the reference script.
>
> **Method note worth keeping.** The route comparison passed at
> **machine precision on every one of these runs**, including while
> `|k1| = 0.03` against `|k2| = 0.38`. It has to: both routes consume the same
> inputs and `close_loop()`'s algebra is already proven, so it agrees no
> matter how wrong the coefficients are. **Every defect above was caught by
> the physics assertions** (power sum, reciprocity), which were added only
> after the first run passed while being wrong. A cross-check between two
> routes over shared inputs validates the routes, never the inputs.

2. **Hand-assemble the `SMatrixBlock`.** `beef.network.convert.
   from_transfer_matrix()` is *not* sufficient here: it produces a
   degenerate 2-port with **zero self-terms**, and the coupler's whole
   physics lives in its self-terms (`t` = bus→bus, `t_ring` = ring→ring).
   The block's `data[m]` must be populated directly at the four
   (port, port) index blocks via `port_slice()`, from two
   `FDTDSolver.solve_section()` runs (bus-excited → `t`, `k1`;
   ring-excited → `t_ring`, `k2`), matching
   `run_coupler_fdtd()`'s existing two-run structure.

3. **Build the loop propagator `P`.** Legacy form (`ring_resonator_hybrid.py`
   :762-774) is `P = b_corrected**2 * exp(1j*phi)`, where `b_corrected =
   b_d * exp(-1j*phi_bend_correction)` strips the straight-run phase off the
   raw bend coefficient using `BEND_D_IN = BEND_D_OUT = STUB_LEN -
   MONITOR_MARGIN`, and `phi = 2π·neff·L_outer/λ` with
   `L_outer = 2*(L_SIDE - S_MON) + L_FAR`. **That reference-plane correction
   is load-bearing for resonance position** — reproduce it exactly or every
   resonance shifts. Note both ring geometry classes' `to_gds_eme()` raise
   `NotImplementedError` (they are FDTD-only, mixed-orientation sections),
   so an `EMESolver` cell cannot cover the bend; the body stays either
   analytic-`neff` (as today) or a raw-polygon meow call.

4. **Close the loop and cross-check numerically.**
   `close_loop(block, port_out="ring", port_in="ring", loop=P)` collapses to
   the surviving bus port; `to_smatrix()` back out gives `E_thru`. Compare
   against the legacy closed form **term by term**, not just visually —
   `tests/test_network_ring_crosscheck.py` already proves `close_loop()`
   reproduces that algebra symbolically, so any disagreement here is in the
   inputs (items 2/3), which is exactly what makes this a useful check.
   `close_loop()` implements only the same-port self-loop, which is this
   device's case; fine as-is.

5. **New example script** replacing both hand-rolled ones, plus a golden.

**Note on `HybridSolver`'s role:** it has **no `solve_section()` of its
own** — that method lives only on `FDTDSolver`/`EMESolver`, and
`HybridSolver.solve()` cascades sections by plain matrix product
(`hybrid.py:1063-1066`) with no feedback. So the ring's loop is closed by
`beef.network`, not by `HybridSolver`; `HybridSolver` can only ever cover
the open path. This doc's item (c) framing ("`HybridSolver` must express
feedback") is therefore better read as "the *feedback* mechanism now exists
in `beef.network`" — rewiring `HybridSolver` itself onto `beef.network` is
listed there as separate, still-undone work and is **not** a prerequisite
for the ring.

Until that wiring lands the ring cannot be validated the way the other
devices are: checks 1–3 in [README.md](../README.md) all assume
`FDTDSolver`/`EMESolver` objects that this device never constructs — the
existing example scripts are unchanged and still don't call any of (a)/(b)/(c).

## Material model — the ring is now the only device left frozen

Added 2026-08-01, when `dispersive` was flipped on by default across the
library and every example. **The ring did not follow, and it fails silently
rather than loudly.** All four items below are consequences of the same root
cause this doc is already about — the ring hand-rolls raw meep/meow instead of
going through `FDTDSolver`/`EMESolver` — so they are cheap once the
`HybridSolver` recommendation above lands, and awkward before it.

~~Nothing here is currently *corrupt*: `example_sim_cache/` holds no ring entries
at all, so items M3/M4 were latent mechanism gaps, not live staleness. They will
become live the first time either ring script is run.~~

> **✅ SUPERSEDED 2026-08-03 — they became live, and they landed on the right
> side of the fix.** `example_sim_cache/` now holds **66 ring entries**, all
> dated 2026-08-02, all written by job 11634900:
>
> | slug | count |
> |---|---|
> | `ring_resonator_hybrid_solve_neff_..._BareStrip` | 18 |
> | `ring_resonator_hybrid_cascaded_solve_neff_..._BareStrip` | 18 |
> | `ring_resonator_hybrid_cascaded_solve_body_supermodes_..._CoupledPairCrossSection` | 18 |
> | `..._run_coupler_fdtd_{bus,ring}_RingResonatorCoupler` | 2 + 2 |
> | `..._cascaded_run_left_fdtd_{bus,ring}_RingResonatorLeftSection` | 2 + 2 |
> | `..._run_bend_fdtd_RingResonatorBend` (both scripts) | 2 + 2 |
>
> **They postdate the M3/M4 fingerprint fix (2026-08-01), so they are correctly
> keyed.** Verified by reading a stored `config.json`: it carries
> `resolved_si_index`, `resolved_sio2_index` and `resolved_probe_wl` alongside
> `geometry_polygons` and every `_sim_params()` constant. Had these entries been
> written a day earlier they would have carried no material fingerprint at all
> and would now be silent staleness traps. **Nothing is corrupt — but that was
> timing, not design.**
>
> Note what the fingerprint honestly reveals: `resolved_si_index` holds **three
> frozen values** (one per probe wavelength), which is M1 visible in the cache
> key. The fingerprint is right; the physics is still half-frozen.

> **Status 2026-08-02 — M3 and M4 are DONE, M1 and M2 are NOT.** The cache
> fingerprint is fixed; the physics is not. The ring's FDTD half still freezes
> one index at band centre while its meow half tracks the band. Read M1/M2 as
> live, M3/M4 as history.

**M1 — the ring's FDTD half is still frozen, and nothing warns.**
`ring_resonator_hybrid.py:258` and `_cascaded.py:284` call
`MATERIAL_STACK.si_medium(wl_center)` with a **scalar**. That is not affected by
the default flip by design — a scalar carries no band to fit, so `si_medium`
correctly returns the frozen index at that exact wavelength. The consequence is
that these two scripts keep freezing one index at band centre while every other
device now tracks the band. Fix: pass the run's `wavelengths` array, not
`wl_center`, at both call sites (and the matching `sio2_medium` at `:265`/`:291`).

**M2 — and its EME half is already dispersive, so the ring now carries the exact
asymmetry the flip just closed everywhere else.** `ring_resonator_hybrid.py:582`
uses `MATERIAL_STACK.extrusion_rules`, i.e. meow `SampledMaterial`s, which
re-evaluate per wavelength. So this device contracts a frozen FDTD result
against an exact-per-wavelength EME mode basis — the 0.00786-in-n_Si mismatch
(SOI_220 numbers; the ring's own band should be re-measured) that
`material_model_validation.md` work item 5 exists to catch. **It is not caught
here**, because that guard lives in `HybridSolver._check_material_index_agreement()`
and the ring never constructs a `HybridSolver`. This is the first concrete cost
of the bypass that shows up as a *wrong number* rather than as missing
functionality.

**M3 — ✅ FIXED 2026-08-01.** All 6 fingerprint sites across both ring scripts
now splat `resolved_material_config()` (`ring_resonator_hybrid.py:380,600`;
`_cascaded.py:384,567,631`). Verified before the fix that the geometry-only key
was **byte-identical** frozen-vs-dispersive, which is the staleness this item
predicted. The original text is kept below because the *mechanism* is the
lesson, and because M1/M2 are still open — see the deliberate `False`/`True`
asymmetry noted under M4.

*Original item:* the ring's hand-rolled cache config omits the material
fingerprint entirely. `_sim_params()` (`ring_resonator_hybrid.py:333-353`) records
`"material_stack": {"wg_height": ..., "sz": ...}` — geometry only. The library's
`resolved_material_config()` additionally records `resolved_si_index` /
`resolved_sio2_index` at three probe wavelengths, plus `resolved_dispersive` and
the fitted `resolved_{si,sio2}_reparam`. **None of that reaches the ring's key.**
So neither the 2026-07-28 Sellmeier change (n_Si shifted ~0.002) nor the
2026-08-01 dispersion flip would move a ring cache entry: a stale run would keep
matching and be served with the old physics.

That is precisely the failure mode `_sim_params`' own docstring describes for
the *spectral* config ("Nothing in dc_validation.md items 1–3 touched this
script's hash, which is exactly why those entries survived (item 10)") —
recurring on the material axis instead. Fix: call
`resolved_material_config(MATERIAL_STACK, wavelengths, <resolved dispersive>)`
and splat it in, rather than hand-copying two geometry fields.

**M4 — ✅ FIXED 2026-08-01, and the two halves deliberately disagree.** The
resolved value is passed explicitly as the third argument at every site — and
it is **`False` on the FDTD side** (still `si_medium(wl_center)`, i.e. M1 is
genuinely unfixed) and **`True` on the meow side**. That asymmetry *is* M2.
Recording it truthfully in the key is the point: the fingerprint now says the
device is half-frozen instead of silently claiming it is not. When M1 lands,
the FDTD side's `False` becomes `True` and the keys turn over correctly.

*Original item:* do not hand-copy the resolved flag, either. When M1 is fixed, the
dispersive setting must reach the key as the *resolved* value, not the stack's
raw flag. The library had exactly this bug until 2026-08-01: every
`resolved_material_config()` call site passed the stack, so
`FDTDSolver(dispersive=True)` on a frozen stack built a genuinely different
medium under a byte-identical key. The function now takes the resolved value as
a third argument; the ring should pass it too. See `notes.md` Solvers item 21.

**Sequencing note.** M1–M4 are all subsumed by the "put the ring on
`HybridSolver`" recommendation above — routing through the solvers gets the
media, the guard and the fingerprint for free, and is the reason none of this is
worth patching twice. If the ring is going to be rewired anyway, do M1–M4 as
part of that rather than before it. What is *not* optional either way: **do not
publish ring numbers computed on the frozen media as if they were current**,
since every other device's numbers moved on 2026-08-01 and the ring's did not.

## ✅ Order item 1 — the re-run, 2026-08-02 (job 11634900, both tasks)

3-D, `FDTD_RESOLUTION = 40`, `MODE_RESOLUTION = 100`, `N_WL_DENSE = 2001`.
`_0` = `ring_resonator_hybrid` (58:06), `_1` = `..._cascaded` (44:18).

| | `_0` monolithic | `_1` cascaded |
|---|---|---|
| power conservation | 0.9998 / 0.9932 | 0.9998 / 0.9928 |
| reciprocity | **\|k1\|=0.2530 vs \|k2\|=0.2530** | max\|k1−k2\|=1.78e-16 *(symmetric by construction)* |
| `neff` fit residual | 7.29e-06 | even 1.45e-05 / odd 7.30e-07 |
| dips found | 7 | 8 |
| **mean FSR** | **12.975 nm** (0.8% from the analytic 12.87) | 13.350 nm (3.7%) |
| at 1.5500 µm | \|t\|=0.9674 \|k1\|=**0.2530** \|t_ring\|=0.9639 \|b\|=0.9988 | \|t\|=0.9713 \|k1\|=**0.2374** \|t_ring\|=0.9677 \|b\|=0.9988 |

**The monolithic reciprocity result is the load-bearing one.** `|k1|` and
`|k2|` are *independently normalized* there and agree to four decimals at
production 3-D. That is the y-normal excite / flux-calibration path blocker 0
accused three separate times, and it is the number that *could* have come out
wrong. The cascaded 1.78e-16 cannot — `S_coupler = M_left.T @ P_body @ M_left`
is symmetric by construction, so it tests the algebra, not the physics. Two
different guarantees, not one result confirmed twice.

### ⚠ NEW, UNDIAGNOSED — the two routes disagree by ~half an FSR

Same device, same 43.42 µm perimeter, same analytic prediction.

```
monolithic: 1.5101 1.5226 1.5352 1.5481 1.5612 1.5744 1.5879
cascaded:   1.5034 1.5161 1.5291 1.5422 1.5555 1.5691 1.5829 1.5969
```

Near band centre the cascaded dip at **1.5555** sits almost exactly midway
between the monolithic **1.5481** and **1.5612** — the two combs are **~half
an FSR out of phase**, i.e. ~π of round-trip phase. The offset drifts from
~0.52 FSR at the short-wavelength end to ~0.37 FSR at the long end, consistent
with the 2.9% FSR difference accumulating across the band. The differing dip
counts (7 vs 8) follow from the offset comb catching one more dip inside the
same 1.50–1.60 µm window, not from a different resonance density.

`|k1|` also differs by **6.6%** (0.2530 vs 0.2374). That accounts for
different dip *depths*; it does **not** account for a comb offset — position
is set by round-trip phase and both routes use the same physical length.

~~**Deliberately not diagnosed here.** The candidate space is phase bookkeeping
— the reference-plane correction on `b(λ)`, or the `L_outer` the cascaded
route assumes — and none of it has been tested.~~

### ✅ CANDIDATE SPACE COLLAPSED, 2026-08-03 — phase bookkeeping is ELIMINATED

**No new FDTD was needed.** Both `ring_data.npz` files retain `t`, `k1`, `k2`,
`t_ring`, `b_raw` as `complex128` (101 wavelengths), so this was post-processing
on saved arrays plus a line-by-line source comparison.

**1. The analytic phase chain is byte-identical in both scripts.** Not
"equivalent" — identical source fed identical constants:

| quantity | monolithic | cascaded | equal? |
|---|---|---|---|
| `L_outer = 2*(L_SIDE - S_MON) + L_FAR` | `:805` | `:864` | **identical, = 7.0 µm** |
| `phi = 2π·neff_dense·L_outer/wl_dense` | `:806` | `:865` | **identical** |
| `phi_bend_correction` over `BEND_D_IN + BEND_D_OUT` | `:803` | `:862` | **identical** |
| `P = b_corrected**2 * exp(1j*phi)` | `:807` | `:866` | **identical** |
| `L_C`, `L_SIDE`, `L_FAR`, `S_MON`, `R`, `STUB_LEN`, `MONITOR_MARGIN` | — | — | **all identical** |
| `neff` fit residual (max abs), printed | 7.29e-06 | 7.29e-06 | **identical** |

The cascaded script's own section header says it: *"Dense-grid assembly
(unchanged from `ring_resonator_hybrid.py`)"*.

**2. `b_raw` is identical to all printed digits — magnitude AND phase.** At
λ=1.5500 µm, difference exactly `0.000000`:

| | monolithic | cascaded | diff |
|---|---|---|---|
| `\|b_raw\|` | 0.998828 | 0.998828 | **0.00e+00** |
| `arg(b_raw)` | −1.275855 rad | −1.275855 rad | **0.000000** |

So the far-bend coefficient and the entire ring-arm round-trip phase are common
to both routes. **Neither named candidate can produce the offset.**

**3. The difference is in the coupler coefficients, and it is in their PHASE:**

| | monolithic | cascaded | diff | diff/π |
|---|---|---|---|---|
| `\|t\|` | 0.967365 | 0.971299 | −3.9e-03 | — |
| `arg(t)` | −1.898810 | −0.629439 | −1.269371 | **−0.404** |
| `\|k1\|` | 0.253006 | 0.237425 | +1.6e-02 | — |
| `arg(k1)` | +2.491699 | +0.621226 | +1.870474 | **+0.595** |
| `\|t_ring\|` | 0.963932 | 0.967696 | −3.8e-03 | — |
| `arg(t_ring)` | +1.502908 | −1.272172 | +2.775080 | **+0.883** |

Magnitudes agree to ~0.4–6%. **Phases disagree by fractions of π of order one.**
Resonance position is set by the round-trip argument, and `arg(t_ring)` — the
coefficient that closes the loop — differs by **0.883π**. That is the same order
as the observed ~0.5 FSR (≈π) comb offset.

**What this licenses, and what it does not.** Established by measurement: the
offset originates in the complex coupler coefficients, not the ring phase
bookkeeping; the two routes' magnitudes broadly agree while their phases do not.
**Not established:** which route's phases are right, or that any single sign or
convention is responsible. This doc's one recorded diagnosis ("π sign flip in
`T_ii`") was confidently wrong and cost a day, and `arg(t_ring)` differing by
0.883π rather than exactly 1.000π is precisely the kind of near-miss that
invites the same error twice. **Measured: where the offset lives. Still not
argued: why.**

⚠ **Separate observation, and it weakens a claim made above.** The reciprocity
check compares **magnitudes only** (`max ||k1| − |k2||`). In the monolithic run
the magnitudes match to 2.4e-05 while the *phases* do not:
`arg(k1) = +2.4917` vs `arg(k2) = +0.2551` rad, a gap of **0.712π**. The
cascaded route has `arg(k1) = arg(k2)` exactly, as its `M.T @ P @ M` construction
forces. Full reciprocity for this device requires `k1 = k2` as complex numbers,
so the monolithic's "reciprocity holds" — called the load-bearing number in
`## Order item 1` — has only ever been tested on half the quantity.
**Caveat before treating this as a defect:** if `k1` and `k2` are read at
different reference planes, a phase offset between them is expected and
physical. ~~That has not been checked.~~ **CHECKED — see below. The planes DO
differ, so a phase offset is expected; but the arithmetic does not close.**
Still worth printing `arg` alongside `abs` in both scripts' summary blocks.

### 🔍 Port planes DO differ between the routes — but they do NOT account for the phases

**The two scripts do not measure their coefficients at the same reference
planes, and the cascaded script's own docstring says so** (`ring_resonator_hybrid
_cascaded.py`, "Note on port-plane consistency — a design constraint the
monolithic version didn't need"). The `M_right = M_left.T` reciprocity trick
requires the *same* ring port for both directions, so the cascaded route puts
its ring port at `Y_RING + R + S_MON` in both directions. The monolithic route
is **asymmetric**:

| monolithic ring-excited run | plane | height above arc end |
|---|---|---|
| source | `Y_TOP - SOURCE_MARGIN` | — |
| **input normalization** | `Y_TOP - MONITOR_MARGIN` (`:471`) | **+1.25 µm** |
| **output** `ring_r_mon` | `Y_RING + R + S_MON` (`:476`) | **+0.50 µm** |

`Y_TOP = Y_RING + R + STUB_LEN`, so the asymmetry is
`(STUB_LEN − MONITOR_MARGIN) − S_MON` = **0.75 µm** of ring stub that the
monolithic's ring-excited coefficients traverse un-compensated. The docstring
notes the stub top was *"only ever used there as an INPUT launch point, with no
matching output use on that same side"* — i.e. the asymmetry is by construction,
not by accident.

**This is a real inconsistency, and it is the tempting explanation. It is also
not sufficient.** Converting the measured `arg` differences into an equivalent
propagation length (`Δarg·λ / 2π·n_eff`, `n_eff@1.55 = 2.339522`):

| λ (µm) | Δarg(t_ring) → path | Δarg(k1) → path | Δarg(t) → path |
|---|---|---|---|
| 1.5000 | −2.8344 → **−0.282** | +2.1213 → **+0.211** | −1.0196 → **−0.101** |
| 1.5500 | +2.7751 → **+0.293** | +1.8705 → **+0.197** | −1.2694 → **−0.134** |
| 1.6000 | +2.1383 → **+0.239** | +1.6318 → **+0.183** | −1.5078 → **−0.169** |

Candidate lengths from the geometry: `d = 0.75`, `2d = 1.5`, `STUB_LEN −
MONITOR_MARGIN = 1.25`. **None of them match any column.** Worse, the three
coefficients require **three different** equivalent lengths (≈0.29, ≈0.20,
≈−0.13 µm), so no single common path offset can reconcile them — and each drifts
with wavelength rather than holding a fixed length.

**What that means.** The 0.75 µm asymmetry is real and should be fixed on its own
merits, but **"the monolithic forgot 0.75 µm of ring stub" is refuted as the
explanation of the comb offset.** Because the cascaded route splits the device
into different FDTD domains, *every* port plane differs between the routes — bus
as well as ring, which is why `arg(t)` (a bus→bus quantity that never touches the
ring port) also disagrees. The phases are not comparable between routes at all
without a full per-coefficient reference-plane accounting, which neither script
does and which `L_outer` does not supply.

### 🛠 DEFECT FOUND 2026-08-03 — `P_body` under-propagates the body by 0.5 µm

**Found by the new assembly ledger, not by inspection.**
`ring_resonator_hybrid_cascaded.py` gained `device_assembly_ledger()` /
`plot_device_assembly()` (section 10b): the round trip drawn and tallied as the
performance algebra assembles it, each section placed between the reference
planes its own coefficients are normalized to.

**The ledger validates itself first.** Summing the *physical* spans reproduces the
`PERIMETER` constant exactly — `43.4159` vs `43.4159`, delta `+0.00e+00`. That
cross-check caught a sign error in the ledger's own right-section term during
development, so it is load-bearing, not decoration.

With that established, exactly one row disagrees:

| span | modelled by | assigned | physical | mismatch |
|---|---|---|---|---|
| left section: ring port → interface | `M_left` (FDTD) | *measured* | 9.1040 | — |
| **body: interface → interface** | **`P_body`** | **2.0000** | **2.5000** | **−0.5000** |
| right section: interface → ring port | `M_left.T` | *measured* | 9.1040 | — |
| side straights (×2), far straight | `L_outer` | 1.5 / 4.0 | 1.5 / 4.0 | +0.0000 |
| far bend arcs (×2) | `b_corrected` | 7.8540 | 7.8540 | +0.0000 |

**The mechanism.** `BODY_LENGTH = L_C − 2*IFACE_BUFFER` = 2.0 µm assumes the body
spans `IFACE_X` to `−IFACE_X`. But the supermode monitor is placed at
**`IFACE_X − MONITOR_MARGIN`** (`:462`), not `IFACE_X` — so the actual span
between the two decomposition planes is

```
(−IFACE_X + MONITOR_MARGIN) − (IFACE_X − MONITOR_MARGIN)
    = BODY_LENGTH + 2*MONITOR_MARGIN = 2.5 µm
```

The shipped comment says the body is shortened by `2*IFACE_BUFFER` "since each
section's own FDTD domain already simulates `IFACE_BUFFER` of the physical `L_C`."
It simulates it but **measures at `IFACE_BUFFER − MONITOR_MARGIN` = 0.75 µm**, so
the shortening should have been `2*0.75`, not `2*1.0`. **The model's round trip is
0.5 µm shorter than the device's actual perimeter.**

### The decisive test — real defect, PARTIAL explanation

Run read-only from the cached `M_left` columns plus the script's own `find_dips`;
no FDTD, nothing written. The shipped case reproduces production exactly
(`|t|=0.9713`, `arg(t_ring)=−1.2722`, FSR 13.350 nm, the 11634900_1 comb to four
decimals), so the harness is faithful.

| | mean FSR | error vs analytic 12.87 | dips near band centre |
|---|---|---|---|
| cascaded, **as shipped** (2.00 µm) | 13.350 nm | **3.7%** | 1.5422 · **1.5555** · 1.5691 |
| cascaded, **corrected** (2.50 µm) | **13.143 nm** | **2.1%** | 1.5391 · **1.5522** · 1.5655 |
| monolithic (reference) | 12.975 nm | 0.8% | 1.5481 · 1.5612 |

**The correction improves both metrics in the right direction.** FSR error halves
(3.7% → 2.1%); the comb offset against the monolithic near band centre drops from
**~0.55 FSR to ~0.31 FSR**. **It does not close the gap.** Something else remains.

**So: a real, geometrically proven bookkeeping defect that accounts for part of
the disagreement and not all of it.** Recorded that way deliberately — the scale
was suggestive enough (0.5 µm ≈ 0.75 of an FSR of bare-strip round-trip phase) to
have been written up as *the* cause on the strength of the estimate alone, which
is the error this doc has already paid for once.

### 📌 ROOT CAUSE — the ring applied the excite-side margin to a measure plane

**The library already does this correctly, and the DC proves it.** In
`_section_layout()` (`beef/solvers/fdtd.py:602`) the two sides use deliberately
different conventions:

```python
src_x, norm_x, mea_x = source_margin, monitor_margin, L
```

- the **excite/norm** side is inset by `source_margin`/`monitor_margin` — forced,
  since a source must sit inside the domain and its normalization monitor is
  co-located with it;
- the **measure** side sits at `mea_x = L`, i.e. **exactly on the section
  boundary**, with `dpml` added outside it (`L = (x1-x0) + stub_left + stub_right`,
  "PML not" included; `Sx = span + stubs + 2*dpml`).

Measure planes landing exactly on the seam is what makes section lengths compose
exactly. `directional_coupler_hybrid.py` depends on it and **audits it**:

```python
body_lengths = coupling_lengths - 2 * BODY_STUB
assert np.allclose(body_lengths + 2 * BODY_STUB, coupling_lengths)
```

That audit is *correct* because its interface monitor really is at
`boundaries[1] = -half_body + BODY_STUB`, exactly as the script prints.

**The ring hand-rolls raw meep and put its interface monitor at
`IFACE_X - MONITOR_MARGIN` — the excite-side convention on a measure plane.** The
library solver would have placed it at `IFACE_X`, and `BODY_LENGTH = L_C -
2*IFACE_BUFFER` would then have been exactly right. So this is not an arithmetic
slip: it is one specific consequence of the ring bypassing `FDTDSolver`, the same
root cause this whole document is about. It is also a concrete, previously
unlisted cost of that bypass — add it to the M1–M4 list's rationale.

### ✅ FIXED 2026-08-03 — the length is now derived from the plane, and the prediction held

**What changed.** The interface plane had been computed in two independent places;
it is now defined once and the length derived from it:

```python
IFACE_PLANE_L = IFACE_X - MONITOR_MARGIN     # the plane actually measured at
IFACE_PLANE_R = -IFACE_PLANE_L
BODY_LENGTH   = IFACE_PLANE_R - IFACE_PLANE_L
```

Both monitor placements (`run_left_fdtd()`'s `iface_mon` and the section preview)
now name `IFACE_PLANE_L` instead of recomputing `IFACE_X - MONITOR_MARGIN`, so a
future change to the plane cannot desynchronize the length again. **This repair is
correct regardless of which plane is chosen** — that is why it was done first.

**Free to apply.** `BODY_LENGTH` only reaches `phi_even`/`phi_odd` → `P_body`,
i.e. post-processing. Supermodes are length-independent, so it appears in no cache
fingerprint. Verified: all three 3-D FDTD entries still hit, and the cache stayed
at 97 entries with none added.

**The prediction was registered before the change and matched exactly.** From the
read-only cache rebuild: mean FSR 13.143 nm and the comb at
`1.5009 1.5134 1.5261 1.5391 1.5522 1.5655 1.5791 1.5929`. The patched script
produced **13.143 nm** and precisely that comb.

**Ledger after the fix:**

```
body: interface -> interface   P_body (supermode phase)   2.5000   2.5000   +0.0000
physical round trip (sum of spans)                                43.4159
PERIMETER constant                                                43.4159   (delta +0.00e+00)
Every assigned length matches its physical span. Bookkeeping closes.
```

**Every metric moved toward the monolithic and the analytic prediction:**

| at λ=1.5500 µm | monolithic | cascaded, before | cascaded, **after** |
|---|---|---|---|
| mean FSR | 12.975 nm | 13.350 nm (3.7% vs analytic) | **13.143 nm (2.1%)** |
| \|t\| | 0.9674 | 0.9713 (0.40%) | **0.9657 (0.18%)** |
| \|k1\| | 0.2530 | 0.2374 (**6.2%**) | **0.2593 (2.5%)** |
| \|t_ring\| | 0.9639 | 0.9677 (0.39%) | **0.9621 (0.19%)** |
| comb offset near centre | — | ~0.55 FSR | **~0.31 FSR** |

⚠ **The route disagreement is NOT closed.** ~0.31 FSR of comb offset survives, dip
counts are still 8 vs 7, and per the earlier per-coefficient analysis the residual
lives in the coefficient *phases*, which this fix does not address. This was a real
defect worth fixing on its own merits; it was never established to be the whole
cause and it is not.

~~⚠ **Q is still grid-limited**~~ — **✅ RESOLVED later the same day, see
`## Q is no longer grid-limited` below.**

### Still open — should the plane move to the library convention?

Deliberately NOT done. `FDTDSolver` puts a measure plane exactly on the section
boundary (`mea_x = L`), insetting only the excite/norm side; the ring keeps a
`MONITOR_MARGIN` inset. Adopting the library convention would:

1. change what the FDTD measures → **all three 3-D runs must be redone**; and
2. ⚠ **hit a cache trap.** `monitor_margin` is in `_sim_params()`, but the
   *derived* monitor position is not. Changing `IFACE_X - MONITOR_MARGIN` →
   `IFACE_X` with `MONITOR_MARGIN` unchanged leaves the fingerprint **identical**,
   so the cache would silently serve the old measurement. Needs a fingerprint
   field or a `_CACHE_SCHEMA_VERSION` bump **first**. This is the standing
   schema-version exposure biting in a concrete place.

📌 **Reframing worth keeping.** `L_outer = 2*(L_SIDE − S_MON) + L_FAR` subtracts
`S_MON`, which matches the *cascaded* convention. The monolithic's ring input
plane sits at `STUB_LEN − MONITOR_MARGIN`, not `S_MON`. So the shared assembly
formula is consistent with at most one of the two routes. **This makes it likelier
that the offset is a reference-plane bookkeeping mismatch than a physics error in
either FDTD** — but that sentence is a hypothesis, and the paragraph above is
exactly why it must not be promoted before it is measured. **Measured: the planes
differ, and the naive length correction fails. Still not established: what the
correct per-coefficient correction is.**

**Why this is good news rather than bad.** It is the first time the ring's two
routes have disagreed *at all*. Per the standing entry below, two routes
agreeing proves little because they consume the same inputs; two routes
disagreeing means at least one is wrong and is cheap to adjudicate — the
analytic FSR already favours the monolithic (0.8% vs 3.7%).

⚠ **[RESOLVED 2026-08-03 — `N_WL_DENSE` raised to 20001, see
`## Q is no longer grid-limited`. Left in place because the reasoning is still
the reason the knob matters.]** Q in both runs is grid-limited, not physical. Every FWHM is an exact
multiple of the 2001-point grid spacing (0.05 nm) — 0.15/0.20/0.25/0.30 =
3/4/5/6 steps. The reported Q (5293–10193) is set by **`N_WL_DENSE`**, a knob
independent of both resolutions. In the 501-point `RING_TEST_MODE` smoke every
FWHM was pinned at exactly 0.4000 nm. **If Q is wanted, raise `N_WL_DENSE`** —
it is post-processing on an interpolated grid and costs no new FDTD.

## Filing status — this doc stays at top level, NOT in `resolved/`

**Asked 2026-08-02, re-asked and re-answered 2026-08-03. Recorded here so it
does not have to be re-litigated a third time.** This analysis previously lived
only in a session context snapshot, which is explicitly not a durable record —
that was the gap.

The bar taken from [README.md](../README.md): `resolved/` means *"whole docs
sufficiently addressed (open to-dos may remain)"*, and the directional coupler is
described there as **the only device where all four checks have actually been
run**. Against that:

| check | ring | is this an "open to-do"? |
|---|---|---|
| 1. Engine agreement | ~~**not run — not runnable.** No shared cross-section object exists~~ **✅ RUN AND PASSES 2026-08-03** — `UniformCrossSection` supplied the shared object. Job 11654797 at the shipped res 100: bare strip **bit-identical** (0.000e+00), coupled pair 1.3–1.6e-10 = its translation floor exactly at all 3 wavelengths | Done |
| 2. Mode identity | ~~**not run — not applicable.** No interface mode basis is exchanged~~ **✅ RUN AND PASSES 2026-08-03** — matched overlap 1.000000, identity permutation, unit completeness, clean self-overlap control | Done |
| 3. Resolution convergence | **MEASURED, and the FDTD side does not converge.** Mode side converged at 100. Job 11653509 (20/30/40/60): FSR converges (0.16% of an FSR at res 60) but **absolute comb position converges at p ≈ 0.70** and is still 29% of an FSR out at res 60 — the signature of a wavelength-flat phase, probably the 0.976π | Yes |
| 4. Full-device FDTD reference | ~~**not run.** No `ring_resonator_fdtd.py` exists~~ **WRITTEN 2026-08-03**, job 11655143 (2-D res 20/30/40). 3-D is 167 h at res 40 — beyond the wall limit | Yes |

~~Checks 1 and 2 are not unfinished work; they are **unattemptable** until the ring
moves onto `HybridSolver`.~~
**↑ Superseded 2026-08-03.** They turned out not to need `HybridSolver` at all —
only a `DeviceGeometry` both engines could consume, which is a much smaller thing
than this paragraph claimed. `UniformCrossSection` (75 lines) unblocked both, and
both now pass. The standing recommendation below is still worth doing for
composability, but it was **never the blocker** these two checks were waiting on.
Worth noting as a case where "structurally blocked" was an over-reading of the
obstacle, and where re-deriving the actual dependency was cheaper than the
migration it seemed to require. Prerequisites (a) and (b) landed 2026-07-29;
**(c) "HybridSolver must express feedback" is the substantive one.** Two further
live items as of 2026-08-03: the **half-FSR route disagreement is undiagnosed**
(see `## Order item 1`), and **M1/M2 leave the device half-frozen** — FDTD side
frozen at band centre while the meow side tracks the band.

**`deferred/` is also wrong**, and for a different reason: that folder means work
postponed *by choice* (the mux, at the user's direction). The ring is under
active work with an unmade gating decision in front of it — **port topology,
4-port `bus_L/bus_R/ring_L/ring_R` vs composed 2-port.** "Blocked and active"
reads differently from both "addressed" and "postponed," and top level is the
only one of the three that says it.

**Nothing in the 2026-08-02 or 2026-08-03 results moved this.** The clamp fix
closed a defect and job 11634859 confirmed it at production config, but neither
touched checks 1–4; the same session *opened* the route disagreement. Progress on
correctness is not progress against the filing bar.

## 📌 CORRECTION — wiring items 2–4 are already DONE (Branch A), 2026-08-03

**A "what is left before `resolved/`" answer given earlier today listed wiring
items 2–5 as outstanding. That was wrong, and the error was not checking the
verification script before summarising the to-do list.**

`examples/verification/python/verify_ring_network_crosscheck.py` already
implements the whole Branch A path this document recommends doing first:

| wiring item | where | state |
|---|---|---|
| 2. hand-assemble the `SMatrixBlock` | `:292`, `:311`, `:369` | **done** |
| 3. build the loop propagator `P` | `loop_propagator()` | **done** |
| 4. close the loop and cross-check | `close_loop(block, port_out="ring", port_in="ring", loop=P)` at `:393` | **done, passes** |

`port_out == port_in` is the same-port self-loop, i.e. **Branch A** — the
implemented branch. Its own module docstring says so at line 4. Job 11634859
measured `worst |E_legacy - E_network| = 1.735e-18` at production 3-D config.

**What actually remains** is wiring item **1** (port topology → needs
`close_loop()` **Branch B**, which raises `NotImplementedError` for
`port_out != port_in`) and item **5** (a production example script replacing both
hand-rolled ones, plus a golden). Checks 1 and 2 stay blocked behind item 1.

So the ring is further along the library path than the gap list implied: Branch A
runs *through* `beef.network` on `RingResonatorCouplerSection` geometry today.

## 🔍 The monolithic has its OWN reference-plane inconsistency — partial result

> **SUPERSEDED the same day — see `## The route disagreement is essentially RESOLVED`
> below. The mixed result recorded here is what motivated building the ledger instead of
> patching further, and that is why it is kept.**

**Tested 2026-08-03, read-only from the monolithic `ring_data.npz`. Applied to
nothing.** After fixing the cascaded body length, the natural next hypothesis was
that the monolithic carries the mirror-image defect — and geometrically it does:

- monolithic ring **input** normalization plane: arc end **+1.25 µm**
  (`Y_TOP - MONITOR_MARGIN`)
- monolithic ring **output** plane: arc end **+0.50 µm** (`Y_RING + R + S_MON`)
- but `L_outer = 2*(L_SIDE - S_MON) + L_FAR` subtracts `S_MON` = 0.5 on **both**
  sides

The input side should subtract 1.25, so the correct lever arm is
`(L_SIDE - 0.5) + (L_SIDE - 1.25) + L_FAR` = **6.25 µm**, not 7.0 — the monolithic
**over**-counts by 0.75 µm where the cascaded **under**-counted by 0.5.

**The test result is genuinely mixed, which is why nothing was changed:**

| | mean FSR | vs cascaded post-fix (13.141) | comb offset near centre |
|---|---|---|---|
| monolithic, shipped `L_outer=7.00` | 12.973 nm | 1.3% apart | ~0.31 FSR |
| monolithic, corrected `L_outer=6.25` | **13.175 nm** | **0.24% apart** | **~0.45 FSR** |

**FSR agreement between the two routes improves ~5× (2.9% → 0.24%). Comb position
agreement gets WORSE (0.31 → 0.45 FSR).**

**And the analytic FSR cannot adjudicate.** `FSR_GUESS` uses `NG_GUESS = 4.3`, a
hand-entered group-index guess the script itself labels *"just for the printed
FSR estimate"*. Treating 12.87 nm as a reference at the 1–2% level is not sound,
so "the shipped monolithic is closer to analytic" is not evidence.

**Reading, held deliberately loose.** The 0.75 µm inconsistency is real and
geometrically demonstrable. Correcting it alone does not improve the thing that
matters, which means **at least one further uncompensated offset exists** in the
monolithic. The leading candidate is already measured: `arg(k1)` and `arg(k2)`
differ by **0.712π** there, because `k1` is normalized at the bus plane and `k2`
at the ring input plane — and `k1·k2` enters `E_thru` directly.

**Next step is the principled one, not another single-term patch:** build a
monolithic round-trip ledger the way `device_assembly_ledger()` does for the
cascaded, covering **all four** coefficients' reference planes at once. Patching
terms one at a time on a device with several interacting offsets is how the
"π sign flip in `T_ii`" day was lost.

## ✅ The route disagreement is essentially RESOLVED — 2026-08-03

**From ~0.55 FSR to ~0.05 FSR, with matching dip counts.** Two monolithic defects,
both found by building the ledger rather than by patching terms.

### The monolithic ledger, and it self-validates

`ring_resonator_hybrid.py` section 10b now carries `device_assembly_ledger()` /
`coefficient_reference_planes()` / `print_assembly_ledger()`, the counterpart to
the cascaded's. Physical spans sum to **43.4159 = `PERIMETER`, delta `+0.00e+00`**,
so the tally is trustworthy before any conclusion is drawn from it. One row
disagrees:

| span | modelled by | assigned | physical | mismatch |
|---|---|---|---|---|
| coupler section: ring_in → ring_out | one FDTD | *measured* | 21.4580 | — |
| right side straight | `L_outer` | 1.5000 | 1.5000 | +0.0000 |
| far bend arcs (×2) | `b_corrected` | 7.8540 | 7.8540 | +0.0000 |
| far straight | `L_outer` | 4.0000 | 4.0000 | +0.0000 |
| **left side straight** | **`L_outer`** | **1.5000** | **0.7500** | **+0.7500** |

`L_outer` as shipped **7.0000**, plane-consistent value **6.2500**. The monolithic's
modelled round trip is 0.75 µm **longer** than the device — the opposite sign to
the cascaded's 0.5 µm shortfall.

### Defect 2: a wavelength-flat π in the ring-excited coefficients

The reference-plane table makes the second defect visible: `k1` and `k2` are **not
referenced to the same ring plane** (`k1` measured at arc end +0.50, `k2`
normalized at +1.25), so `max ||k1| − |k2||` was never a like-for-like test.

Referencing `k2` onto `k1`'s plane and asking whether the residual is flat
(a convention constant) or drifting (another length) — **the discriminator was
chosen before looking**:

| λ (µm) | \|k1\| | \|k2\| | raw Δarg/π | after +βd | **after −βd** |
|---|---|---|---|---|---|
| 1.5000 | 0.1974 | 0.1973 | 0.5775 | 0.1752 | **0.9798** |
| 1.5250 | 0.2238 | 0.2237 | 0.6457 | 0.3137 | **0.9778** |
| 1.5500 | 0.2530 | 0.2530 | 0.7119 | 0.4479 | **0.9760** |
| 1.5750 | 0.2852 | 0.2852 | 0.7762 | 0.5781 | **0.9744** |
| 1.6000 | 0.3206 | 0.3205 | 0.8386 | 0.7043 | **0.9730** |

**The −βd column is flat: 0.976π ± 0.003 across 100 nm**, while the raw gap drifts
0.26π. Magnitudes agree to four decimals. So after the plane correction the two
differ by a **wavelength-independent constant near π** — a sign, not a length.

### Applying both collapses the disagreement

Read-only from the monolithic `ring_data.npz`, using the script's own
`interp_complex()`/`find_dips()`:

| case | dips | FSR | offset vs cascaded |
|---|---|---|---|
| monolithic as shipped | 7 | 12.973 nm | ~0.31 FSR |
| + `L_outer` → 6.25 only | 7 | 13.175 | ~0.45 FSR |
| + ring-excited sign only | 8 | 12.976 | ~0.23 FSR |
| **+ BOTH** | **8** | **13.177** | **~0.053 FSR** |
| cascaded (post-fix, reference) | 8 | 13.141 | — |

Comb, both corrections: `1.5015 1.5140 1.5268 1.5397 1.5529 1.5663 1.5799 1.5937`
against the cascaded's `1.5009 1.5134 1.5261 1.5390 1.5522 1.5655 1.5791 1.5929`
— **0.7 nm apart, uniformly, with the dip count finally matching (8 vs 8)** and FSR
within 0.27%. Both corrections are needed; neither suffices alone.

### Mechanism localized — and why it was never caught

`get_coeffs(..., forward=False)` returns `alpha[0, :, 1]`, the opposite-traveling
index, and its own docstring says it is needed **"exactly once in this script"** —
the ring-excited run's norm monitor, whose source launches −y. That single call is
the normalization denominator for exactly the two coefficients carrying the π
(`k2` and `t_ring`), and for neither of the two that don't.

⚠ **This is the gap in blocker 0's closure.** Blocker 0 was the y-normal excite
normalization; it was closed 2026-08-01 (11616920/11616921) on the finding that the
ratio converges as h² with fitted **p = 1.988**, i.e. discretization rather than a
bug. **That test was magnitude-only** — `max ||k1| − |k2||`. A π in the phase is
invisible to it. The closure was right about what it measured and silent about
this.

### ❌ MECHANISM REFUTED — the backward index is sign-consistent, 2026-08-03

**`verify_ring_ynormal_sign.py`, and it says no.** Two mirror-image runs on a
uniform y-oriented strip (2-D, res 30, 5 wavelengths): source at y=−S launching
+y read via `alpha[...,0]`, versus source at y=+S launching −y read via
`alpha[...,1]`, each monitor `MONITOR_MARGIN` downstream of its own source. Mirror
symmetry demands the two coefficients be identical.

```
|ratio| deviation from 1 : 1.289e-06
arg(ratio)/pi mean       : +0.00000
arg(ratio)/pi spread     : 2.509e-07
VERDICT: CONVENTIONS AGREE (ratio ~ +1)
```

**`alpha[..., 1]` is NOT sign-flipped relative to `alpha[..., 0]`.** The
`get_coeffs(forward=False)` call is not the source of the ring's flat π, and
**`get_coeffs()` must not be changed.** Straight guide, no bend, no coupler — the
one geometry where the expected answer is known exactly, so this is not a
near-miss to be argued with.

**This is the second mechanism proposed for that π and the second one refuted**
(the first, in an older session, was "π sign flip in `T_ii`"). The measurement —
0.976π flat across 100 nm with magnitudes equal to four decimals, and a 10×
collapse in route disagreement when applied — is unaffected. Only the explanation
died. **Measured: the π and what it fixes. Still not established: where it comes
from.**

⚠ **What the test does and does not cover.** It compares
`source(−y) + alpha[...,1]` against `source(+y) + alpha[...,0]` as a *pairing*, so
it would not detect a sign that cancels within that pairing. More importantly it
compares **y-normal against y-normal**. The ring's π is between `k1` (bus-excited,
**x-normal** source) and `k2` (ring-excited, **y-normal** source) — a
**cross-orientation** comparison this test never makes.

**Next candidate, and the next test:** does an x-normal port carry the same phase
convention as a y-normal one? Same uniform strip, run once along x and once along
y (a 90° rotation of an identical structure), compare the extracted coefficients.
That is the only remaining difference between how `k1` and `k2` are obtained, and
it is squarely in the territory blocker 0 examined with magnitudes only.

### ❌ SECOND MECHANISM REFUTED — port orientation is not it either, 2026-08-03

`verify_ring_ynormal_sign.py` test 2: **cross-orientation reciprocity on the ring's
own standalone bend**, the simplest structure with one x-normal and one y-normal
port. Excite the x stub and measure the y stub (`c_A`); excite the y stub and
measure the x stub (`c_B`). The bend is reciprocal and both runs use the same two
planes at the same two distances, so reciprocity demands `c_A == c_B` as **complex**
numbers — this is the ring's `k1`/`k2` pairing with no coupler, no feedback and no
second bend.

```
arg(ratio)/pi mean   : +0.00040
arg(ratio)/pi spread : 5.086e-03
VERDICT: ORIENTATIONS AGREE (ratio ~ +1)
```

**An x-normal port and a y-normal port carry the same phase convention.** Combined
with test 1 (backward vs forward index, `+0.00000π`), **both port-convention
explanations for the ring's flat π are now dead.**

**Where that leaves it.** The measurement is unchanged: 0.976π flat across 100 nm,
magnitudes equal to four decimals, and applying it (with the plane fix) collapses
the route disagreement 10× and fixes the dip count. But after **two controlled
refutations** the origin is unknown, and that raises the weight on a possibility
worth stating plainly: **0.976π is not π.** The 0.024π offset and its 0.007π drift
may mean the correction is not a sign at all, but something else that happens to
sit near π. **Do not apply it to the scripts.** An empirical constant that works
without a mechanism is exactly what the "π sign flip in `T_ii`" entry warns about.

📌 **A separate finding from the same test, and it is blocker 0's signature.**
Magnitudes did *not* agree: `|ratio|` deviates by **4.9%**, and `|c_B|` comes out
**above 1.0** (1.00559–1.00884) — unphysical transmission for a lossy bend, i.e.
the y-excited run's normalization is a few percent high. That is precisely blocker
0's symptom, at res 30, and consistent with its h² closure (p = 1.988,
extrapolating to 0.979 at res 40). **So blocker 0's magnitude story reproduces on
a bare bend, while its phase was never examined by either investigation.**

## ✅ CHECK 4 RAN, AND IT LOCALIZED THE DEFECT — 2026-08-03

Job 11655143 (2-D, res 20/30/40) against the 2-D hybrid partner (11655349), both
pre-M1 so matched. Comparison in `verify_ring_check4_comparison.py`.

| res | FSR difference | comb offset vs cascaded hybrid |
|---|---|---|
| 20 | −0.603 nm | +3.610 nm (24.8% of an FSR) |
| 30 | −0.143 nm | +5.123 nm (34.4%) |
| 40 | −0.192 nm | +4.639 nm (30.9%) |

**FSR agrees to ~1.2% and converges. Absolute comb position disagrees by ~30% of an
FSR and does NOT shrink with resolution.**

That was the pre-registered discriminator and it came out on the bookkeeping side:

> the flat phase lives in the hybrid's reference-plane arithmetic, not in meep's
> treatment of the staircased arcs.

The full-device run shares **none** of the hybrid's phase bookkeeping — no reference
planes, no decomposition, no network algebra — so a resolution-independent
disagreement cannot be discretization the hybrid is faithfully reproducing. This is
what rules out the explanation the p ≈ 0.70 comb drift had made plausible.

⚠ **At 34% the offset approaches the half-FSR aliasing limit**, so the line pairing
could be off by one. The script now warns on this itself. The monolithic route is the
cross-check (data from 11656190, not yet compared).

## 🔴 THE k1/k2 RESIDUAL IS MONOLITHIC-ONLY — and the cascaded route's reciprocity is a TAUTOLOGY, 2026-08-03

**This should have been checked six mechanisms ago.** The residual was measured on the
**monolithic** route's `ring_data.npz` and treated throughout as "the ring's" defect.
It is not. Measured on both:

| route | plane offset | max \|\|k1\|−\|k2\|\| | arg(k2)−arg(k1), corrected | drift |
|---|---|---|---|---|
| monolithic | 0.75 µm | 1.320e-04 | −2.9798 .. −2.9730 π | 0.0068 π |
| **cascaded** | **0.00 µm** | **2.220e-16** | **0.0000 π exactly** | **0.0000 π** |

**The cascaded route has no residual at all.** Its ring port sits at `LEFT_RING_Y` in
*both* directions by construction — that symmetry is precisely what makes
`M_right = M_left.T` a valid reuse — so there is no plane offset to correct and
nothing to disagree.

### ⚠ But that zero is a TAUTOLOGY, not a measurement

`S_coupler = M_left.T @ P_body @ M_left`, so `k1` and `k2` are **the same matrix
element** read from either side. `2.220e-16` is machine epsilon: the cascaded route
does not *verify* reciprocity, it **enforces** it. Same for `arg(k2/k1)` at 6.3e-16.

So one of the entries I listed as a validation basis is not evidence. Sorting the
cascaded route's coefficient checks honestly:

| check | value | real measurement? |
|---|---|---|
| reciprocity \|k1\| = \|k2\| | 2.220e-16 | ❌ **tautology** — same matrix element |
| reciprocity arg(k1) = arg(k2) | 6.341e-16 rad | ❌ **tautology** |
| **unitarity \|t\|²+\|k1\|²** | **0.999647 – 0.999941** | ✅ **real** — nothing forces the FDTD-measured `M_left` to conserve power; deviation **3.5e-04** |
| \|t_ring\| vs \|t\| | 3.768e-03 | ⚠ partly structural (flows from the `M_left.T` reuse) |

**Unitarity to 3.5e-04 on an FDTD-measured transfer matrix is the substantive
coefficient-level result**, and it is genuinely independent: `M_left` comes from two
meep runs and nothing in the algebra makes `|t|²+|k1|²` come out at 1.

### What this means for the two routes

* **The cascaded example** — the one being validated — **does not have this defect.**
  The six-mechanism hunt was chasing a monolithic bookkeeping error.
* **The monolithic route has a real defect**: its asymmetric ring port (+1.25 in,
  +0.50 out) leaves a residual that no plane correction removes, and which six tested
  mechanisms do not explain. It is not a blocker for the cascaded example; it is a
  reason to prefer the cascaded route, and a reason the monolithic's comb positions
  should not be quoted.
* **The ~0.05 FSR route disagreement** is therefore between a route with a known
  plane-bookkeeping defect and a route without one. That is no longer a mystery
  requiring a tiebreaker — it is an expected consequence.

### ❌ SIXTH MECHANISM REFUTED — `get_norm_calibrated()` is phase-transparent

Tested before the above was noticed, and still worth recording since the monolithic
defect is open. `verify_ring_norm_calibrated_phase.py`:

* `arg(get_norm_calibrated) − arg(coeff)` = **5.35e-17 π**. The helper takes its phase
  from the modal coefficient and changes nothing — it cannot inject a π.
* `mp.get_fluxes()` is negative at **11/11** points and `np.abs()` does discard that
  sign — but correctly: the sign belonged to a *magnitude*, and direction is carried by
  the separately-taken phase.
* The calibration changes the normalizer **magnitude** by a factor **1.2662**. That is
  real, and it is what "blocker 0" was actually about.

⚠ **This script's first verdict was wrong and I nearly reported it as a find.** It
compared against the *forward* coefficient at that monitor and got 0.995 π. But the
source launches −y, so the forward coefficient there is a **weak reflection**:
`|a_fwd|/|a_bwd| = 4.4%`. Two tells were already on screen and I read past both — the
magnitude ratio, and the fact that the "π" varied **0.73–0.99 π** across the band where
a genuine sign convention would sit at exactly π. The verdict now checks whether the
comparison quantity is even meaningful before ruling on it.

### ⛔ RETIRED — a full-device FDTD is NOT a valid reference for this device, 2026-08-03

**Decision: stop validating the ring by comparison against a full-device FDTD.**
Not deferred, not caveated — the reference is disqualified by its own output.

**The transmission of a passive device exceeds 1:**

| run | T max | excess | points T>1 | % of band |
|---|---|---|---|---|
| res 20 | 0.9829 | −0.0171 | 0 | 0% (the run that resolved only 2 lines) |
| res 30 | 1.0671 | **+0.0671** | 186 | **46.4%** |
| res 40 | 1.0837 | **+0.0837** | 187 | **46.6%** |

Up to **8.4% gain, over 46% of the band, and WORSENING with refinement**
(0.067 → 0.084). That last fact is what makes it a disqualification rather than a
coarse-grid artifact: a discretization error refines away, and this does not. The
likely mechanism is the bus-only normalization — the ring perturbs the through-port
mode near the coupler, so `|a_dev|² / |a_ref|²` is not a true transmission — but the
cause hardly matters once the reference reports gain in a lossless passive device.

**Two independent reasons, and either alone is sufficient:**

1. **Unphysical output** (above). The reference violates energy conservation.
2. **The reference does not converge.** Between res 30 and 40 its own comb moves
   **−2.357 nm = 18.1% of an FSR**, per-line shifts −4.0 to −0.3 nm (so the FSR moves
   too, +4.3%), against a hybrid-vs-reference offset of ~30% of an FSR. The arbiter
   carried uncertainty of the same order as the thing it was arbitrating.

⚠ **Everything the "check 4 verdict" section below claims is therefore withdrawn as a
validation result.** The measurement happened and the numbers are real, but they do
not adjudicate the hybrid. Kept in place because the reasoning — why the observable
had to be group delay, the 3-D cost table, the artifact at res 20 — is worth keeping,
and because a withdrawn conclusion should be visible rather than deleted.

Job 11656189_2 (check 4, res 40, post-M1) was **cancelled at 1:58:19** rather than
finished: more output from a disqualified reference has no value.

### 📌 WHY a full-device FDTD cannot referee a high-finesse resonator

Finesse `F = Q·FSR/λ = 14294 × 13/1550 ≈ 120`. The field circulates ~120 times, so
**any per-round-trip phase error is amplified ~120×** — a 0.5% round-trip error moves
a resonance a full linewidth. FDTD's staircased arcs give O(h) boundary-position
error, which *is* a round-trip path-length error. The reference is therefore subject
to exactly the amplification it is supposed to be checking, and at the resolutions a
full ring can afford it has no accuracy advantage over the model.

⚠ The device is **linear** throughout — no Kerr, no TPA. It is *resonant*, which for
validation is a harder problem than nonlinearity would be, because the amplification
is structural rather than intensity-dependent.

### ✅ WHAT VALIDATES THIS DEVICE INSTEAD

Absolute resonance wavelength is **out of scope** — not validatable at any affordable
resolution, by the hybrid or by an FDTD reference. The validation rests on quantities
that carry no finesse amplification:

| basis | status |
|---|---|
| **Engine agreement** (check 1) | ✅ closed — bare strip bit-identical, pair = translation floor |
| **Mode identity** (check 2) | ✅ closed — overlap 1.000000, identity permutation |
| **Network algebra vs the closed form** | ✅ routes A/B/C agree to machine precision |
| **Coefficient reciprocity** \|k1\| = \|k2\| | ✅ agree to 1.3e-04 |
| **Coefficient unitarity** \|t\|²+\|k1\|² | ✅ ~1 |
| **FSR** | ✅ converges p ≈ 1.4–2.1, 0.16% of an FSR at res 60 |
| **Mode-solve resolution** (check 3, EME half) | ✅ converged at 100 |
| **FDTD resolution** (check 3, FDTD half) | ⚠️ measured; sets the absolute-position limit |
| `k1`/`k2` **phase residual** | ⚠️ **OPEN** — the one remaining defect |

Every entry above is a **single-pass or algebraic** quantity. That is the point: they
are measurable to 1e-4 or better and they do not depend on accumulated round-trip
phase. The `k1`/`k2` phase residual stays open precisely because it is also
single-pass, so it can be chased cleanly — unlike the comb.

### 🔴 CAVEAT — the reference is NOT converged either, and the verdict above is overstated

**Measured 2026-08-03, after the verdict was written.** The full-device FDTD was
treated as the arbiter without checking whether it converges. It does not:

| res | lines | FSR | comb shift from previous |
|---|---|---|---|
| 20 | **2** | 80.25 nm | — (unusable post-M1) |
| 30 | 7 | 14.500 nm | — |
| 40 | 7 | 15.125 nm | **−2.357 nm = 18.1% of an FSR** |

Per-line shifts run −4.0 to −0.3 nm, so it is not even a rigid translation — the FSR
moves +4.3% too. **The reference carries ~18% of an FSR of its own uncertainty at the
only two usable resolutions**, while the hybrid-vs-reference offset varies ~11% of an
FSR across those same resolutions.

**What survives:** the offset does not trend toward zero with refinement, which still
argues against pure discretization.
**What does not:** "~30% of an FSR, attributable to the decomposition's bookkeeping"
as a clean number. That is overstated and should be read as ~30% ± something of the
same order.

### 📌 WHY a resonator resists this kind of validation — the finesse amplifies everything

This is the deeper reason, and it reframes the whole check:

    F = Q·FSR/lambda = 14294 x 13/1550 ~ 120

The field circulates ~120 times, so **any per-round-trip phase error is amplified by
~120×**; a 0.5% round-trip error moves a resonance by a full linewidth. FDTD's
staircased arcs give O(h) boundary-position error, which *is* a round-trip
path-length error. That is the mechanism behind the hybrid's p ≈ 0.70 comb-position
convergence, and behind the reference's own wobble above. **It is not a bug in either
model — it is what a high-finesse resonator does to a discretization error.**

⚠ Note the ring's response is **not** nonlinear in the optical sense — no Kerr, no
TPA, linear passive device throughout. It is *resonant*, which for validation purposes
is a harder problem than nonlinearity would be, because the amplification is
structural rather than intensity-dependent.

**Consequence for what this device can claim.** Absolute resonance wavelength is the
quantity the finesse amplifies and the grid corrupts, in the hybrid *and* in a
full-device FDTD. It is **not validatable at any affordable resolution by either
model.** Two quantities are well behaved and should carry the validation instead:

* **FSR** — converges at p ≈ 1.4–2.1, within 0.16% of an FSR at res 60. Set by n_g and
  perimeter, not by accumulated round-trip phase.
* **The coupling coefficients** `t`, `k1`, `k2`, `t_ring` — single-pass quantities with
  no finesse amplification. `|k1|` vs `|k2|` already agree to `1.3e-04`.

So check 4's useful output is the **FSR comparison** (−0.14 to −0.19 nm at res 30/40,
~1.2%, and converging) and not the comb offset. The `k1`/`k2` phase residual remains
worth chasing precisely because it is single-pass: measurable cleanly, unlike the comb.

### ⚠ THE OBSERVABLE HAD TO CHANGE — amplitude dips were the wrong measurement

The first version of `ring_resonator_fdtd.py` copied the hybrid scripts'
`find_dips()`, looking for `|T|²` minima. For this device that is wrong, and the job
showed it:

| res | gap in pixels | ringdown | amplitude dips | group-delay peaks |
|---|---|---|---|---|
| 20 | 5.0 | 3.6 τ | 5 | **7** |
| 30 | 7.5 | 0.3 τ | **0** | **7** |
| 40 | 10.0 | 0.4 τ | 3 | **7** |

The same device cannot resonate at res 20 and not at res 30. This is an **all-pass**
ring — one bus, no drop port — and a lossless all-pass ring has **|T| ≡ 1
identically**; its resonances live entirely in the PHASE. Amplitude dips therefore
exist only in proportion to loss and **shrink as the discretization improves**: the
better the simulation, the less there is to see.

Res 20's apparent resonance was a **5-pixel-gap subpixel-smoothing artifact** (the
0.25 µm gap smeared into artificially strong coupling) plus staircase scatter, and its
comb was junk on its own terms — spacings 14.0 / 14.5 / 14.75 / **30.25** nm (a missed
line) and a mean FSR of 18.4 against 12.87 predicted. **It must not be used.**

`τ_g = −dφ/dω` peaks at every resonance for any coupling and any loss including zero,
because the phase winds 2π per FSR regardless. Amplitude dips are still reported as a
secondary diagnostic whose depth measures loss.

⚠ **Peak-finding needs a half-FSR minimum separation** or each resonance is reported
two or three times — 22 "peaks" for 7 at res 40, in triplets 1.25 nm apart between
real 11.5–12.25 nm gaps. `ring_resonator_fdtd.py` and
`verify_ring_check4_comparison.py` briefly disagreed about how many resonances the
same spectrum held because only the latter had the filter.

### 3-D full-ring is out of reach, and that is measured

Q = 14294 → τ = Q/(2πf) = **3526 meep time units**, and
`stop_when_dft_decayed(1e-9)` stops near **26 τ**. Scaling job 11653509_3's measured
0.060 s/step by the 4.00× cell-volume ratio and by res³:

| res | 2-D | 3-D |
|---|---|---|
| 20 | 0.2 h | 10 h |
| 30 | 0.7 h | 53 h |
| 40 | **1.7 h** | **167 h** |

**3-D exceeds the 10 h wall at every usable resolution.** ⚠ An earlier estimate said
45 h, from assuming T = 8 τ; 8 τ leaves 1.8% residual amplitude and the criterion
actually used stops at 26 τ. Both numbers are kept in the script's docstring because
the 45 h figure had already been used to justify defaulting to 2-D.

## ✅ M1 FIXED — the ring's FDTD half now tracks the band, 2026-08-03

`si_medium(wl_center)` with a **scalar** froze one index at band centre while the
meow half re-evaluated per wavelength: n_Si off by **+0.004338 @1.50** and
**−0.003913 @1.60**, i.e. **+1.42% / −1.20% of an FSR** in round-trip phase.

**Why the ring alone was left behind.** The library fixed this on 2026-08-01 inside
`FDTDSolver._media()` — six lines: resolve the flag, pass the whole band. Every
example that calls `FDTDSolver` inherited it with no edit. Grepping every direct
caller of `si_medium`/`sio2_medium` outside the library returns **only the three ring
scripts**. They hand-roll `mp.Simulation` and never construct an `FDTDSolver`, so
they could not inherit a solver fix. **M1 is not a separate bug from the hand-rolled
bypass — it is that bypass, surfacing as a wrong number.**

`_media()` is now replicated in all three scripts (deliberately line-for-line, not a
second opinion about what "dispersive" means) and `make_sim()` takes the band.
Verified:

* `eps_inf` **12.0946 → 12.2220** with one Lorentzian — that is
  `meep_medium_reparam()`, and it matters: a naive dispersive medium *loses* 0.0496 in
  n_eff at res 40 because meep subpixel-averages only the instantaneous ε and meep's
  Sellmeier form puts ε∞ = 1, leaving no Si/SiO₂ contrast to smooth.
* `gamma = 0`, `Im(n) = 0` — lossless. Index accurate to 6e-4, contrast to 0.03%.
* The fingerprint gains `resolved_dispersive: true` and the first run wrote **new**
  cache entries rather than being served frozen media — M3/M4's machinery working on
  its first real test.

**Not a regression.** Res 30 post-M1 resonates (spacings 1.25, 1.25, **11.50**, …).
Res 20's collapse from 3.6 τ to 0.1 τ was that resolution losing its 5-pixel-gap
artifact, which is correct behaviour.

⚠ **Every 3-D result on disk predates M1** — 0 of the 3-D cascaded FDTD cache entries
carry `resolved_dispersive=True`. So the top-level production `ring_data.npz`,
`smatrix_hybrid.pkl` and every figure beside them describe the frozen model. Job
11657558 re-runs the cascaded 3-D production; **the monolithic is still stale.**

## ✅ Bend-straight accounting — DC-style is now the default, 2026-08-03

`phi_bend_correction` strips `b`'s internal 2×1.25 µm and lets `L_OUTER` carry the
full 7.00 µm. The DC-style alternative never touches `b` and shortens the analytic
length to `L_OUTER − BEND_STRAIGHT_TOTAL = 2.00 µm`, i.e. **0.25 µm per side straight
and 1.50 µm across the far straight**.

**They are algebraically identical, and stay so even when `b`'s internal phase
corresponds to a different index than the analytic β.** Tested with
`n_true = n_analytic + 0.025` — exactly the error under investigation:
`max |P1 − P2| = 7.2e-15`. Both reduce to `b² · exp(iβₐ · 2.00)`; β only ever
multiplies the 2.00 µm remainder, and neither form re-derives the FDTD's internal
straight phase. **So this cannot be the phase error**, and switching moved no number:
dips **bit-identical**, `max |diff| = 0.000e+00`.

Switched anyway, for robustness: `b_corrected` was labelled *"the arc alone"*, which
is only true if the analytic `n_eff` describes the bend's 2.5 µm — the exact
assumption under investigation; and `n_eff` entered in **two** places that had to stay
in sync, which would silently stop cancelling the moment a bend-corrected index is
used for one and not the other. Phase-strip is auto-selected when
`BEND_STRAIGHT_HEADROOM < 0`, where it is mandatory. A **standing assertion** now
verifies equivalence every run (`2.5e-14 in P`), replacing a one-off measurement.

### ⚠ Three separate figure/ledger defects, and I reported the fix three times before it was fixed

1. **Ledger and labels quoted 1.50 / 4.00** while the code propagated 2.00 µm total.
   Now **0.25 / 1.50**, with `assigned` *and* `physical` moving together — the DC form
   draws the span boundaries in different places — so every mismatch stays `+0.0000`
   and PERIMETER still closes at **43.4159, delta +0.00e+00**.
2. **`bend_strip` spanned the full 2.00 µm stub-through-PML**, not `b`'s 1.25 µm
   reference plane. Drawn solid (it is *kept* under the DC form) it covered the side
   straight. Now clipped to `BEND_D_IN/OUT`.
3. **The analytic spans were hollow**, so the 1.50 µm far-straight gap read as
   continuous blue — the pale EME wash sits inside it and thick green dashes frame it.
   Now white-filled under the DC form; the old *"don't fill, it would hide the
   hatch"* reason no longer applies because the hatch is clipped away.

Everything keys off one flag, so a negative headroom reverts arithmetic, spans, stub
styling, fill and labels together.

**Each of the first two was checked by eye on a downscaled figure where 1.5 µm is a
few pixels wide, and reported as done while still wrong.** A probe on the drawing
code settled it in one step — `x −7.225..−2.000 arc | −2.000..−0.750 stub |
+0.750..+2.000 stub | +2.000..+7.225 arc`, gap **1.50 µm**. Read the geometry from
the code, not from the PNG.

⚠ **Also: only one of seven `geometry_assembly.png` files was regenerated** at first.
The `res_sweep/` + `two_d/` output isolation added the same day — specifically so runs
could not overwrite each other — multiplied the figure into seven copies. All three
2-D ones are now current; the four 3-D ones are stale pending the re-run.

## ⚠ HAZARD — `run_ring_resonator_examples.batch` is destructive, 2026-08-03

It runs `rm -rf "../../example_output/ring_resonator/${NAME}"` before `make`, to give
each device a fresh tree. That was safe when one device owned that path. **It is not
safe now:** the same tree holds `res_sweep/res_{20,30,40,60}/` — the gap-3 FDTD
convergence sweep behind the p ≈ 0.70 finding, including a res-60 run that cost
2h07 — and `res_sweep/res_*/two_d/`, check 4's comparison partner. The subfolder
redirects added the same day silently widened what that `rm -rf` destroys.

`batch/run_ring_cascaded_3d_postm1.batch` (job 11657558) is a **non-destructive**
replacement for the cascaded re-run, and everything under
`example_output/ring_resonator/` was copied to **`archive/ring_pre_m1_20260803/`**
(31 MB, 13 `ring_data.npz`) first. **Do not submit the old batch until that `rm -rf`
is scoped.**

Separately: the monolithic script **never read `RING_FDTD_RES`** and had no
`res_sweep` redirect — only the cascaded did. Job 11655349 therefore ran it **three
times at res 40 into one file**; two of its six results never existed. Both scripts
are symmetric now.

## 🔴 THE RECORDED "flat 0.976π" READING IS WRONG — re-measured on 101 points, 2026-08-03

The two refutations above were aimed by a misreading, and this is the correction.

The table under *"Defect 2: a wavelength-flat π"* has **5 wavelengths** and was read
as `0.976π ± 0.003` — flat, therefore *"a sign, not a length"*. Recomputed from the
production `ring_data.npz` at all **101 FDTD wavelengths**, with the same 0.75 µm
plane correction:

| model fitted to the corrected residual | rms |
|---|---|
| pure constant (**the recorded reading**) | 1.98e-03 |
| constant + propagation length | **8.59e-05** |
| exactly π + a length | 7.28e-02 |
| quadratic in 1/λ | **8.99e-06** |

* **It is not flat.** The residual drifts **0.0068π monotonically** across the band.
  A free propagation term beats a pure constant by **23×**. What was read as
  `±0.003` scatter is a trend.
* **It is not a clean sign either.** Forcing the constant to exactly π and fitting
  only a length is **850× worse** than leaving it free — so "π plus a small path
  length" is also refuted.

**Why this matters more than the arithmetic.** The discriminator was sound and was
pre-registered correctly (*flat → convention; drifting → another length*). Its
**verdict** was misread, and that misreading is what pointed mechanisms 1 and 2 at
**sign conventions** — where both died. Every hypothesis tested against this residual
so far has been wavelength-INDEPENDENT, and the residual is wavelength-DEPENDENT.

The free fit gives `const ≈ 0.87π`, `d ≈ −35 nm`, but 1/λ spans only 6.5% across this
band so that split is poorly constrained. **Quote the existence of the drift, not
those two numbers.**

### ❌ THIRD MECHANISM REFUTED — the stub cross-section is the bare strip, 2026-08-03

The natural reading of a dispersive residual is that `β = 2πn_eff/λ` uses the wrong
`n_eff`: the transport planes sit in the ring stub, not in an isolated straight. If
the stub's cross-section differed, the transport index would differ, and the error
would drift. Tested directly, at the ring's own `MODE_RESOLUTION = 100`:

| λ | stub cross-section | bare strip |
|---|---|---|
| 1.50 | 2.402279 | 2.402279 |
| 1.55 | 2.339523 | 2.339523 |

**Bit-identical.** And it has to be: the stub is the vertical guide at
`x = ±STUB_X = ±7.0`, and its nearest neighbour — the bus, the near straight, the
other stub — is **≥ 5 µm away**. It is isolated straight guide, so its mode *is* the
bare strip's. The transport `n_eff` is not wrong for cross-section reasons.

⚠ **Also corrected here:** `RING_IN_PLANE` (+1.25) is on the **LEFT** stub and
`RING_OUT_PLANE` (+0.50) on the **RIGHT** — opposite sides of the device, not two
heights on one stub. The 0.75 µm `exp(−iβd)` correction is still the right
correction, because the two stubs are identical mirrored guide, but the earlier
description of it as a single 0.75 µm span was wrong. See
`plot_reference_planes()` in `ring_resonator_hybrid.py`, added 2026-08-03, which
draws both planes on their actual stubs against the arc-end datum.

### ❌ FIFTH MECHANISM REFUTED — the field past the arc IS modal, 2026-08-03

`examples/verification/python/verify_ring_stub_phase_evolution.py`. Mode monitors up
the bend's output stub at 13 distances from the arc end, with a **straight-guide
control** in the same cell so a positive result could not be a harness artifact.
Local index from the phase gradient between adjacent planes; 2-D, res 20.

| quantity, over the ring's own 0.50-1.25 um window | value |
|---|---|
| worst \|n_local - relaxed\| , BEND | **0.000207** |
| worst \|n_local - relaxed\| , straight CONTROL | 0.000016 (harness floor) |
| what the transport residual demands | **~0.025** |
| shortfall | **121x too small** |
| \|a\| spread across all 13 planes | **2.2e-04** relative |

**Refuted on two independent grounds.** The local index departs from its relaxed
value by 121x less than the residual needs, in the very window where `k1` and `k2`
are referenced. And `|a|` is **flat**, not rising — a relaxing superposition projects
less onto the fundamental near the arc and more further out, and there is no such
trend. The field is modally pure essentially immediately.

⚠ **2-D caveat.** A 2-D slab confines more tightly than the real 450x220 nm strip, so
it may relax faster. The 121x margin is wide enough that 3-D is unlikely to reverse
it; `RING_SPE_3D=1` closes that if wanted.

⚠ **The script's first verdict rule was wrong** and would have sent a decided question
back for compute: it asked "is the excess > 0.005?" and called a 121x-too-small
result *inconclusive, raise the resolution*. The bar is not an absolute threshold, it
is the 0.025 the residual actually demands. Fixed.

### 🔬 The hypothesis this refuted (kept for the record)
*Below was the surviving hypothesis before the measurement above. It is now dead.*

### 🔬 THE SURVIVING HYPOTHESIS — the field near the arc is not the fundamental mode

Assume the residual is a π sign flip and solve for the `n_eff` the transport would
have needed:

| λ | needed | bare strip | Δ |
|---|---|---|---|
| 1.500 | 2.42249 | 2.40228 | **+0.02021** |
| 1.550 | 2.36433 | 2.33952 | **+0.02480** |
| 1.600 | 2.30544 | 2.27662 | **+0.02883** |

A smooth **+0.025 (1.06%)** offset, and a physically sensible guided Si mode. Since
the cross-section is provably right, the remaining explanation is that **the field
0.50 µm past an arc is not the pure fundamental mode.** A bend excites a
superposition; over the transport distance that mixture's phase does not evolve at
any single `n_eff`. That is consistent with every constraint:

* dispersive residual ✓ (a modal mismatch is wavelength-dependent)
* correct cross-section ✓ (measured above)
* correct bend-stub accounting ✓ (see the DC-style section)
* defect in the bookkeeping, not the physics ✓ (check 4)

**None of the four refuted mechanisms touched this.** Δ is not perfectly constant
either (it drifts +0.0086 across the band, ~7% of the bare index's own dispersion),
so a single wrong constant `n_eff` is not the whole story — which is itself what a
non-modal field would look like.

**Next step, and it is a measurement:** monitor the phase along the stub at several
distances from the arc end and infer the *local* effective index from the phase
gradient. If it starts ~+0.025 high near the arc and relaxes toward the bare strip
with distance, the hypothesis is confirmed and the fix is to reference `b` and the
coupler coefficients at a plane far enough out for the field to have relaxed — or to
carry the measured phase rather than an assumed modal one.

## ✅ Checks 1 and 2 are now RUNNABLE — `UniformCrossSection`, 2026-08-03

**What "not runnable" actually meant.** This doc recorded check 1 as *"not
measurable today — the mode solve is a direct meow call, and the FDTD is raw meep,
so there is no shared cross-section object to compare."* Chasing that down: the
blocker is that **both ring geometry classes raise on `to_gds_eme()`**, by design —

> `RingResonatorCouplerSection.to_gds_eme()`: *"covers only the mixed-orientation
> coupler section, never a whole-device EME cell layout — the straight bus+ring
> 'body' cross-section is built directly (raw polygons + `mw.compute_modes_meep`),
> not through a DeviceGeometry"*

So there was no object an `EMESolver` and an `FDTDSolver` could both consume. That
is a much smaller gap than "not runnable" implies, and it is now closed.

**`beef/geometry/uniform.py` — `UniformCrossSection`.** N parallel straight strips,
constant along x. Covers both of the ring's cross-sections:

```python
UniformCrossSection([(0.0, W)], length=...)                 # solve_neff()'s bare strip
UniformCrossSection([(0.0, W), (Y_RING, W)], length=...)     # solve_body_supermodes()' pair
```

**Verified end to end, not just unit-tested:**

```
EMESolver.solve_section on UniformCrossSection(coupled pair) -> (1, 2, 2)
  |S| = [[1. 0.]
         [0. 1.]]
bare strip -> (1, 1, 1) |S| = 1.0
```

`|S|` being the identity is the **correct** answer, not a null result: supermodes
are the eigenbasis of a uniform section, so its transfer matrix is diagonal with
unit magnitude and all the content is in the *phase*. That is precisely the
structure the ring builds by hand as `P_body` — so the EME path reproduces its
shape from a shared geometry object.

**15 tests** (`tests/test_geometry_uniform.py`). The load-bearing one pins the
drawn geometry against the ring's own constants — both core y-spans, the `GAP`
between them, the oxide span — because this only unblocks the checks if it draws
*exactly* what the ring builds by hand. It also asserts, rather than absorbs, that
`to_gds_eme()` centres its oxide on y=0 while the ring centres on `BODY_Y_CENTER`,
so a caller must shift it.

Every cut is legal here (no `_check_cut`), which makes it the one geometry where a
cut-position invariance sweep has a known-exact answer — useful for check 2.

⚠ **Runnable is not run.** The demonstration above is res 14 against the ring's
res-100 cached supermodes, so it shows the *path* works, not that the engines
agree. **Check 1 proper is now a matched-resolution comparison** of `EMESolver`'s
supermode `n_eff` against `solve_body_supermodes()`'s — a real run, and the next
concrete step. Check 2 (mode identity) needs the interface basis exchanged across
the same object, which this makes possible for the first time.

## ✅ Checks 1 and 2 RUN — and PASS at the shipped res 100, 2026-08-03

`examples/verification/python/verify_ring_eme_engine_agreement.py`. Mode solve
only, no FDTD. **Job 11654797, res 100 = the ring's shipped `MODE_RESOLUTION`,
three wavelengths, 17.6 min.** (A res-20 login-node smoke run passed first at
`9.116e-08`.)

| cross-section | λ | ring `n_eff` | beef `n_eff` | \|Δ\| | translation floor |
|---|---|---|---|---|---|
| bare strip | 1.500 | 2.402279484 | 2.402279484 | **0.000e+00** | 0.000e+00 |
| bare strip | 1.550 | 2.339522768 | 2.339522768 | **0.000e+00** | 0.000e+00 |
| bare strip | 1.600 | 2.276618695 | 2.276618695 | **0.000e+00** | 0.000e+00 |
| pair, even | 1.500 | 2.408939184 | 2.408939184 | 1.520e-11 | — |
| pair, odd | 1.500 | 2.392196768 | 2.392196768 | 1.285e-10 | **1.285e-10** |
| pair, even | 1.550 | 2.349084582 | 2.349084582 | 1.977e-11 | — |
| pair, odd | 1.550 | 2.326731377 | 2.326731377 | 1.347e-10 | **1.347e-10** |
| pair, even | 1.600 | 2.289906149 | 2.289906149 | 3.564e-11 | — |
| pair, odd | 1.600 | 2.260436363 | 2.260436363 | 1.615e-10 | **1.615e-10** |

Check 2 at every one of the six mode/wavelength combinations: matched overlap
`1.000000`, identity correspondence, row/col completeness `1.000000`, clean
same-engine control (max off-diagonal `1e-06`).

**Two different results, worth not conflating.** The **bare strip is
bit-identical** — `0.000e+00`, not "small" — because it needs no y-shift at all
(only the coupled pair's oxide is off-axis), so that comparison is completely
free of the translation. The **coupled pair matches its own translation floor
to every digit at all three wavelengths**, so its `~1e-10` is the rigid shift
and nothing else. Neither number contains an engine disagreement.

**The number that matters is that TEST 1 equals TEST 0 to every printed digit.**
The two are independent measurements. TEST 0 is a *within-engine* control: the
ring's own recipe, run twice, once shifted by `-BODY_Y_CENTER`. That shift is
forced — `section_gds()` centres its oxide on `y=0` while the ring centres on
`BODY_Y_CENTER`, so an exact reproduction is impossible without it, and pretending
otherwise would compare two different structures. TEST 1 is the *cross-engine*
comparison. Their agreeing at `9.116e-08` says the entire residual is the
deliberate translation and **none of it is the engines**. Had the script only
printed TEST 1, `9e-08` would have looked like a small engine disagreement rather
than what it is: zero.

Check 2 goes further than "same `n_eff`". Two different structures can share an
`n_eff`, so the script takes the field overlap matrix between the two engines'
mode sets, and reports:

* **matched diagonal** 1.000000 — same modes, not just same numbers;
* **identity permutation** — the ring's `(even, odd)` band index carries across
  engines. This is the failure that would corrupt `P_body`'s diagonal, and it is
  now excluded rather than assumed;
* **row/col completeness** 1.000000 — neither engine found a mode the other's
  basis cannot represent. The diagonal alone is blind to that;
* **a same-engine control** — the ring's own basis against itself, which must be
  the identity by construction. Two confident-looking "correspondence failures"
  elsewhere in this project turned out to be artifacts of the *checker*; this
  control is what distinguishes them from real findings.

**Resolutions are not free parameters here.** `compute_modes_meep` needs `dy ==
dz`, and the coupled pair's `BODY_MESH_WIDTH = 4.7` does not divide evenly at
every `res`. Legal: 20, 30, 40, 50, 60, 80, 100. **Illegal: 14 and 25** — the
first smoke attempt was at res 14 and died inside meow three minutes in. The
script now refuses an illegal resolution up front with the legal list, and does
*not* silently round to a legal grid, which would change the structure under test.

**Scope, stated in the script's own verdict.** This is the EME/mode-solve half
only. It does not compare against FDTD, and it does not test propagation —
`solve_section()`'s S-matrix, reflection and reverse transmission stay out of
scope. Checks 1 and 2 close for the mode solve; check 4 is still the only thing
that puts the whole device against a reference.

## ✅ The physical 4-port route now RUNS — and what it does and doesn't buy

**Added to `verify_ring_network_crosscheck.py` 2026-08-03 as route C.** The 4-port
block already existed there (`connect(left, "iface", right, "iface")` →
`('bus_L','ring_L','bus_R','ring_R')`) but was collapsed to the `[bus, ring]`
aggregate before closure. It is now also closed *as a 4-port* with the new
`close_loop_pair()`, keeping `bus_L`/`bus_R` distinct.

**C1 is a control, and it passes exactly.** With the counter direction dead
(`loop_ab=0`) the paired closure must reproduce the aggregate's Branch A result,
because the aggregate's "ring" port *is* the one-directional `ring_L`-in /
`ring_R`-out identification:

```
CONTROL  worst |E_4port(1-way) - E_aggregate| = 0.000e+00   -> PASS
route C ports after pair closure: ('bus_L', 'bus_R')
```

Exactly zero. **The 4-port port bookkeeping is verified**, which is what had to be
true before any 4-port number meant anything.

### ⚠ C2 agreeing is a TAUTOLOGY, and the script now says so itself

Closing both circulation directions also gives `0.000e+00` — and that is **not** a
physics result. The script now measures the terms the second loop would act on:

```
backward terms: S[ring_L,bus_L]=0.00e+00  S[ring_L,ring_R]=0.00e+00
                S[bus_L,ring_R]=0.00e+00  S[bus_L,bus_R]=0.00e+00
-> TAUTOLOGY, not a physics result.
```

All structurally zero, because `_side_block` fills **forward transmission only** —
its own docstring: *"reflection and reverse transmission are zero, matching what
`solve_section()` actually models."* The counter-circulation loop has nothing to act
on, so 2-way *must* equal 1-way.

I had first written this verdict as "the counter path is dark, as a co-propagating
coupler implies," which reads as a measurement of the device. It is not. The
detector was added so the script cannot make that claim on anyone's behalf again.

### What this changes for the gap list

**Wiring item 1 is no longer blocked.** The 4-port route runs, its bookkeeping is
verified against the proven aggregate, and `close_loop_pair()` has independent
mutation-checked tests. What remains under item 1 is a decision, not a blocker:
whether `P_LR` and `P_RL` are one propagator (equal for a symmetric racetrack).

**But a NEW gap is now named:** the 4-port route currently buys **composability,
not new physics**. To differ from the aggregate at all it needs side blocks
carrying reflection and reverse transmission, and `solve_section()` does not model
them. Anything requiring backward paths — add-drop, cascaded rings, backscatter,
counter-propagating resonances — is blocked on *that*, not on the network layer.
This is the same class as the item-2 structural limit (no basis-free flux from
`solve_section()`): the network algebra is ahead of what the solver can measure.

📌 **`RING_XC_OUT` added.** The script `np.savez()`es unconditionally, and on
2026-08-02 two login-node smoke runs used that to destroy the production
`crosscheck.npz`. This doc's standing advice was "redirect `OUT` before running
outside the queue" — which required a code change to be possible. Now it is:
`RING_XC_OUT` redirects the output tree, and today's smoke runs used it.

## 📌 Wiring item 1 rescoped — Branch B as specified is not what the ring needs

**Derived 2026-08-03 from `close_loop()`'s actual algebra.** Branch B is defined as
`a_{port_in} = loop @ b_{port_out}` with `port_in = p != q = port_out`. Working it
through: with `E` the remaining ports and `a_p` forced by `b_q`,

```
b_q = (I - S_qp L)^-1 (S_qq a_q + S_qE a_E)
b_E = S_Ep L b_q + S_Eq a_q + S_EE a_E
b_p = S_pp L b_q + S_pq a_q + S_pE a_E
```

**The free inputs are `a_q, a_E`; the outputs are `b_p, b_q, b_E`.** That is not a
valid S-matrix reduction: `b_q` is consumed by the loop yet still appears as an
output, and `a_q` and `b_p` are left dangling. A one-way `p != q` closure is
**under-determined for the ring**, because a real ring circulates both ways —
closing only `ring_R -> ring_L` leaves `a_{ring_R}` and `b_{ring_L}` unconnected.

**What the ring actually needs is a PAIRED closure**, both circulation directions
at once:

```
a_L = P_RL b_R          and          a_R = P_LR b_L
```

Then both ring ports are fully eliminated and `bus_L, bus_R` remain — a well-posed
2-port. With `x = [b_L; b_R]`:

```
[ I - S_LR P_LR      -S_LL P_RL  ] [b_L]   [S_LE]
[   -S_RR P_LR     I - S_RL P_RL ] [b_R] = [S_RE]  a_E

b_E = S_EL P_RL b_R + S_ER P_LR b_L + S_EE a_E
```

A single 2n×2n solve per wavelength, structurally the same shape as Branch A's
`(I - S_pp P)^-1` but on the port pair. **This is implementable and well-posed,
where the one-way Branch B is not.**

**Consequence for the to-do list.** Wiring item 1 was recorded as "decide the port
topology, then implement Branch B." The decision is narrower than that: the 4-port
route needs a **`close_loop_pair()`-style paired elimination**, not the one-way
Branch B that `ops.py` currently raises on. The existing `NotImplementedError`
message is accurate about Branch A vs B being different formulas, but names an
operation that does not close the ring by itself. **Not implemented here** — this
is scoping, and it needs its own tests plus a decision on whether `P_LR` and `P_RL`
are independent (they are equal for a symmetric racetrack, but an add-drop or a
cascade would not be).

### ⚠ NOT patched, and deliberately so

**What is established:** the plane inconsistency (geometric, ledger-verified); the
π being flat across the band (measured); that applying both collapses the route
disagreement 10× and fixes the dip count; and that the suspect path is a single
documented `forward=False` call.

**What is NOT established:** that meep's `alpha[..., 1]` actually carries a −1
relative to `alpha[..., 0]` for a y-normal port in this geometry. That is a claim
about meep's convention and it has not been tested. **This document already
records one confidently-wrong "π sign flip in `T_ii`" diagnosis that cost a day**,
and inserting a bare `-1` on the strength of a curve fit would be that mistake
again with better arithmetic.

**Next step, and it is a measurement:** a controlled 2-D straight-waveguide test —
launch −y, read the same coefficient via the forward and backward indices at two
planes, and compare against the +y case. Cheap, decisive, and it settles whether
the fix belongs in `get_coeffs()` (a convention bug affecting anything using
`forward=False`) or in the ring script alone.

## Q is no longer grid-limited — `N_WL_DENSE` 2001 → 20001, 2026-08-03

**The knob was raised and the convergence was measured, not assumed.** Both ring
scripts now use `N_WL_DENSE = 20001` (2001 in `RING_TEST_MODE`), giving 0.005 nm
steps across 1.50–1.60 µm.

Rebuilding the dense assembly from the *same* cached coefficients at increasing
density — no new FDTD, and using the scripts' own `interp_complex()`/`find_dips()`
rather than a reimplementation:

| `N_WL_DENSE` | step (nm) | first-dip FWHM (nm) | first-dip Q |
|---|---|---|---|
| 2001 | 0.0500 | 0.1500 | 10006 |
| 5001 | 0.0200 | 0.1400 | 10720 |
| 10001 | 0.0100 | 0.1100 | 13644 |
| **20001** | **0.0050** | **0.1050** | **14294** |
| 40001 | 0.0025 | 0.1050 | 14294 |
| 80001 | 0.0013 | 0.1050 | 14294 |

**Flat from 20001 on** — 4× and 8× further refinement change nothing. So 20001 is
converged, and Q is now a property of the device rather than of the grid.

**Every previously reported Q was an under-estimate**, by up to ~43% (10006 →
14294 on the first dip). Production values after the change:

| | dips | Q range |
|---|---|---|
| monolithic | 7 | 6757 – 13131 |
| cascaded | 8 | 5792 – 14294 |

against the old 5293–10193. `T_min` also deepened (e.g. 0.4145 → 0.3512): at
0.05 nm steps the sampled minimum was simply missing the true one.

⚠ **A first attempt at this table was wrong and is worth recording.** Rebuilding
with plain real/imag interpolation gave Q still climbing at 80001 and never
converging. The scripts interpolate **magnitude + unwrapped phase**
(`interp_complex()`), which is not the same thing for a rapidly rotating
coefficient — real/imag chords across the circle and distorts the lineshape. Dip
*positions* were unaffected (those are set by the analytic `exp(iφ)` on the dense
grid), which is why the earlier `BODY_LENGTH` test still stands, but every FWHM
was wrong. **Extract the script's own helpers rather than re-deriving them.**

This closes the "if Q is wanted, raise `N_WL_DENSE`" item outstanding since
2026-08-02. It cost no FDTD and no cache invalidation — the dense grid is
post-processing and appears in no fingerprint.

## Gaps against the four checks

1. ~~**Engine agreement** — not measurable today. The mode solve is a direct meow
   call, and the FDTD is raw meep, so there is no shared cross-section object to
   compare. Blocked on the recommendation above.~~
   **✅ MEASURED AND PASSES 2026-08-03.** `UniformCrossSection` supplied the
   shared object; `verify_ring_eme_engine_agreement.py` runs the comparison.
   Worst `|Δn_eff|` at res 20 is `9.116e-08`, **equal to the script's own
   within-engine translation-floor control to every printed digit** — so the
   residual is the forced `-BODY_Y_CENTER` y-shift and none of it is the
   engines. Res-100 run: job 11654797. See the dedicated section above.
   ⚠ Scope: this is EME-vs-EME through two different call paths (hand-rolled
   meow vs `EMESolver`), **not** EME-vs-FDTD. It closes the "no shared object"
   blocker and shows the library path reproduces the ring's own numbers; an
   FDTD comparison of the same cross-section is still unwritten.
2. ~~**Mode identity** — not applicable in the current form (no interface mode
   basis is exchanged; the scripts thread scalar coefficients).~~
   **✅ MEASURED AND PASSES 2026-08-03**, same script and same run. Field
   overlap matrix between the two engines' mode sets: matched diagonal
   `1.000000`, **identity permutation** (so the ring's `(even, odd)` band index
   carries across — the failure that would corrupt `P_body`'s diagonal is now
   excluded rather than assumed), unit row/col completeness (neither engine
   found a mode the other's basis cannot represent), and a same-engine control
   that comes back clean. The scripts still thread scalar coefficients in
   production; what changed is that the bases can now be compared at all.
3. ~~**Resolution convergence** — never measured, for either engine.~~
   **✅ MEASURED 2026-08-01 — job 11616610 (4:34). The shipped `RESOLUTION = 40`
   slides the resonance comb by a QUARTER of an FSR.**

   `verify_ring_gap3_neff_convergence.py` sweeps `solve_neff()`'s exact recipe
   (same meow `Mesh2D`/`CrossSection`/`compute_modes_meep`, same `SOI_220`,
   same `MESH_WIDTH/HEIGHT`) and converts each mesh step into the quantity
   that actually sets resonance position — the round-trip phase
   `φ = 2π·n_eff·L_outer/λ`:

   ⚠ **The Δφ columns below are superseded — the lever arm was 2× too large.
   See "Lever arm CORRECTED" immediately after the table.** The as-run numbers
   are kept so the correction is auditable.

   | res | n_eff @1.55 | Δφ @1.55 (rad) *as run, 2× too large* | fraction of 2π *as run* | **Δφ CORRECTED** | **fraction of 2π CORRECTED** |
   |---|---|---|---|---|---|
   | 20 | 2.229033 | — | — | — | — |
   | 30 | 2.283999 | +3.1194 | 50% | **+1.5597** | **24.8%** |
   | 40 (**shipped**) | 2.311706 | +1.5724 | 25% | **+0.7862** | **12.5%** |
   | 60 | 2.323453 | +0.6666 | 11% | **+0.3333** | **5.3%** |
   | 80 | 2.338308 | +0.8430 | 13% | **+0.4215** | **6.7%** |
   | 100 | 2.339523 | +0.0690 | 1.1% | **+0.0345** | **0.55%** |

   ### 🛠 Lever arm CORRECTED, 2026-08-03 — every Δφ above halves

   `verify_ring_gap3_neff_convergence.py:31` hardcodes

   ```python
   L_OUTER = 2 * 2.0 + 10.0        # = 14.0
   ```

   on a comment asserting `L_C=10.0` is "the shipped scale". **The shipped
   `L_C` is 4.0** (`ring_resonator_hybrid.py:125`), and both device scripts
   compute the lever arm with a reference-plane pullback the convergence script
   omits entirely:

   ```python
   L_outer = 2 * (L_SIDE - S_MON) + L_FAR    # 2*(2.0-0.5) + 4.0 = 7.0
   ```

   (`ring_resonator_hybrid.py:805`, `ring_resonator_hybrid_cascaded.py:864` —
   identical.) **Job 11634859 printed `L_outer = 7.0000 um` directly**, so this
   is measured, not re-derived. The convergence script's 14.0 is wrong twice
   over: wrong `L_C`, and no `−2·S_MON` term.

   Because Δφ scales linearly in `L_outer`, **every Δφ is exactly half what was
   reported.** The n_eff column, the mesh sweep, and the res-100 convergence
   conclusion are untouched — as gap 3's own caveat correctly predicted they
   would be.

   **Reading, corrected and still serious.** Going from the shipped res 40 to
   res 100 moves n_eff by **0.0278**, i.e. a round-trip phase shift of
   **0.789 rad = 12.6% of a full FSR** — an **eighth** of an FSR, not the
   quarter previously recorded. The direction of the conclusion does not change:
   *resonance wavelengths at the shipped resolution are still positioned
   materially by the mesh.* But the headline figure was overstated 2×, and
   anything downstream quoting "a quarter of an FSR" or "50%" needs the halved
   number. `RING_TEST_MODE`'s res 20 is ~37% of an FSR from res 40, not 75%.

   **This does NOT explain the half-FSR route offset** in `## Order item 1`.
   ~~Both device scripts compute `L_outer` identically, so no `L_outer` error can
   separate them~~ — and a 2× lever-arm error would move FSR by 2×, while the two
   measured FSRs differ by 2.9%.

   ⚠ **CORRECTED AGAIN 2026-08-03 — the lever arm was wrong a THIRD time, and the
   two scripts no longer compute `L_outer` identically.** The struck sentence above
   was true when written and stopped being true the same day, when the DC-style
   bend accounting moved 5.00 µm of straight out of the analytic phase and into the
   FDTD-measured `b` — **for the cascaded route only**:

   | route | source | analytic straight = lever arm |
   |---|---|---|
   | monolithic (phase-strip) | `ring_resonator_hybrid.py:1148` | 7.00 µm |
   | **cascaded (DC-style, shipped default)** | `ring_resonator_hybrid_cascaded.py:1657` | **2.00 µm** |

   Nothing in the verification script changed to cause this; a hardcoded
   `L_OUTER = 7.0` simply stopped describing the shipped default.

   **Re-measured, job 11658491** (4m35s, mode solve only). Shipped res 40 → res 100,
   Δn_eff = +0.027384 … +0.028121:

   | route | Δφ, res 40 → 100 | as a fraction of an FSR |
   |---|---|---|
   | monolithic (7.00 µm arm) | 0.773 – 0.803 rad | **12.3 – 12.8%** — confirms the recorded 12.6% |
   | **cascaded (2.00 µm arm)** | **0.221 – 0.229 rad** | **3.5 – 3.7%** |

   ✅ **Independent cross-check inside the same run:** res-100 `n_eff@1.50` came out
   **2.402279**, against the **2.402279484** control that checks 1 and 2 measured as
   bit-identical across resolutions. Agreement to the printed precision (4.8e-07)
   confirms this script's mode-solve recipe really is the one the ring uses, which
   is the assumption the whole gap-3 study rests on.

   The per-step table also reproduces the recorded **non-monotonicity at 60→80**
   (worst-case 0.0550 → 0.0672 of an FSR, monolithic arm), so a two-point trend
   still must not be extrapolated. Converged at res 100: 0.61% of an FSR
   (monolithic) / 0.17% (cascaded).

   ⚠ **It is NOT a 3.5× improvement in accuracy.** The 5.00 µm did not stop
   mattering — it moved into the bend coefficient, whose own resolution error is
   the p ≈ 0.70 problem. The sensitivity was reallocated between two studies, not
   reduced. Anyone quoting 3.6% must quote the FDTD resolution result alongside it.

   The script now derives both arms from `L_OUTER` and `BEND_STRAIGHT_TOTAL`,
   prints one Δφ table per route with the worst value as a fraction of an FSR, and
   stores `lever_arms` in the npz, so the next accounting change surfaces as a
   number instead of a stale comment. Re-run: **job 11658491**. The npz on disk
   still held `L_outer = 14.0`, the *first* wrong value, because the previous fix
   was source-only and the artifact was never regenerated.

   📌 **The two routes now ship different (algebraically identical) bend
   conventions.** `max |P_phase-strip − P_DC-style| = 7.2e-15`, so this is not a
   physics difference and not a candidate for the route offset. It is still an
   inconsistency between two files that exist to be compared — open decision
   whether to move the monolithic to DC-style too.

   **Converges at res 100** (last step 1.1% of 2π). That is the same tier as
   the 2×2 and one below the 1×2's 120, so the ring is not an outlier — it was
   simply never measured. **Non-monotone at 60→80** (0.667 → 0.843 rad), like
   the 1×2's gap 2, so do not extrapolate from a two-point trend.

   ~~**Caveat on the lever arm:** `L_outer` was taken as `2·L_SIDE + L_FAR` =
   14.0 µm from the shipped constants. The Δφ column scales linearly with it,
   so a different true `L_outer` rescales the phase but not the n_eff column
   or the res-100 convergence conclusion.~~
   **↑ This caveat cashed in on 2026-08-03: the true `L_outer` is 7.0 µm and the
   phase column has been rescaled. The caveat's scoping was exactly right — the
   n_eff column and the res-100 conclusion did survive untouched.** Worth noting
   as a case where flagging an un-derived constant at the time it was used made a
   later correction cheap and bounded instead of invalidating the run.

   **Consequence for order:** this gate is now measured and it is *ahead of*
   the re-run, exactly as step 2 of `## Order` predicted. Re-running either
   script at res 40 would produce resonance positions already known to be a
   quarter-FSR off. Decide the resolution first. Unlike the MMIs — where
   holding res 40 was a defensible cost call against a few-percent
   transmission error — a quarter-FSR comb slide is not a caveat, it is a
   wrong answer for the one quantity this device exists to report.

   Raw data: `example_output/verification/ring_gap3_neff/neff_convergence.npz`.

   ### 🔴 FDTD SIDE MEASURED 2026-08-03 — job 11653509, and it does NOT converge

   The sweep the split above was made to enable. `FDTD_RESOLUTION` 20/30/40/60
   on the cascaded route, mode side pinned. **Control: `neff_solve[0]` is
   `2.402279484` bit-identical at all four resolutions**, so everything that
   moved is the FDTD half and nothing leaked in from the mode solve.

   | res | dips | comb (first 4, µm) | FSR over shared lines |
   |---|---|---|---|
   | 20 | 8 | 1.5043 1.5167 1.5293 1.5421 | 13.124 nm |
   | 30 | 8 | 1.5021 1.5146 1.5273 1.5402 | 13.212 nm |
   | 40 (**shipped**) | 8 | 1.5009 1.5134 1.5261 1.5390 | 13.243 nm |
   | 60 | **7** | 1.5121 1.5248 1.5378 1.5509 | 13.268 nm |

   ⚠ **The as-printed "mean FSR" is not comparable across rows and must not be
   read.** Res 60 reports 7 dips, not 8: the comb shifted down far enough that
   the first line slid below the 1.500 µm window edge (it sat only 0.9 nm
   inside at res 40). Averaging over a different set of intervals made the
   printed mean jump 13.141 → 13.268 nm, which looks like divergence and is
   pure bookkeeping. Compared over the **7 shared lines** the FSRs are 13.243
   vs 13.268. Every number in the table above is measured that way.

   **Two quantities, two completely different convergence orders.** Fitting
   `f(h) = A - B·h^p` on the comb:

   | quantity | fitted order | res 40 | res 60 |
   |---|---|---|---|
   | FSR | p ≈ 1.4–2.1 | 0.35% of an FSR out | **0.16% out** |
   | absolute comb position | **p ≈ 0.70** | 39% out | **29% out** |

   The position order is `0.680 / 0.719 / 0.703` from the three independent
   triplets (20,30,40), (30,40,60), (20,40,60) — consistent, and **sub-linear**.
   Per-step comb shift: 20→30 **−1.86 nm**, 30→40 **−1.11**, 40→60 **−1.24**.
   Not shrinking.

   ### 📌 What that signature means — this is probably the flat 0.976π

   A *difference* of positions converging at ~h² while the *absolute* position
   crawls at h^0.7 has one natural cause: a **wavelength-flat phase offset**. It
   shifts every dip equally and **cancels exactly in the FSR**. That is the same
   signature as the unexplained flat 0.976π in the ring-excited coefficients
   (two mechanism hypotheses already refuted by direct test), so these are
   likely **one defect, not two** — and this is the first evidence that ties
   them together.

   ### ⛔ Consequence — do NOT set `FDTD_RESOLUTION` from this sweep

   Picking a value here would be reading a number off a curve that is not
   converging. Refinement cannot fix it: at p ≈ 0.70 the extrapolated
   requirement to reach 1% of an FSR is far outside anything runnable, which is
   the finding rather than a prescription.

   **What can be trusted today:** FSR, which is converged at res 40+ (0.35% of
   an FSR) and barely improves at 60. **What cannot:** absolute resonance
   wavelength, at *any* resolution, until the flat phase is found.

   This also retires the standing "decide the resolution first, then re-run"
   ordering. That ordering assumed the FDTD side had a convergent answer to
   decide on. It does not, so M1/M2 and the monitor-plane convention are no
   longer gated behind a resolution choice — they are gated behind the flat
   phase.

   ### 🔬 Check 4 is now the discriminator, and it is running

   The full-device FDTD (`ring_resonator_fdtd.py`, job 11655143) shares **none**
   of the hybrid's phase bookkeeping — no reference planes, no decomposition, no
   network algebra.

   * Full-device comb *also* drifts at p ≈ 0.7 → the flat phase is
     discretization (meep resolving the staircased arcs) and the hybrid is
     faithfully reproducing it.
   * Full-device comb converges cleanly while the hybrid's does not → the flat
     phase is in the **decomposition's reference-plane arithmetic**, which
     localizes 0.976π to code rather than physics.

   Comparison harness written and waiting: `verify_ring_check4_comparison.py`.
   Its 2-D hybrid partner is job 11655349, run with `RING_THREE_D=0` (one knob)
   rather than `RING_TEST_MODE=1` (five knobs, one of them a grid-limited
   `N_WL_DENSE`). The harness aligns combs by best offset rather than comparing
   `mean(diff(·))`, precisely because of the 8→7 dip trap above.

   ### ✅ DECIDED 2026-08-01 — split the knob, don't bump it

   `RESOLUTION` was **one** constant driving both the meep FDTD
   (`mp.Simulation(resolution=...)`) and the meow/MPB mode-solve `Mesh2D` —
   the identical defect the coupler split its own knob to kill (*"shared
   RESOLUTION was how a 7% error hid for months"*, `dc_validation.md` item
   11), carried by this device unnoticed. **Gap 3 above measured only the
   mode-solve half.**

   ```
   FDTD_RESOLUTION = 40.0    # meep -- NOT converged-tested, unchanged
   MODE_RESOLUTION = 100.0   # meow/MPB -- gap 3 converged
   ```

   **Why not simply set both to 100.** The mode solve is what sets resonance
   position, it converges at 100, and it is cheap — the entire six-resolution
   sweep took 4:34. The FDTD side has never been measured at *any* resolution,
   and 3-D FDTD cost scales roughly as res⁴, so bumping it would pay ~39× for
   an unmeasured claim. Splitting is what makes those two independent
   statements expressible at all; a single knob forces one answer onto two
   different questions.

   13 call sites updated across both scripts (5 in `ring_resonator_hybrid.py`,
   8 in `_cascaded.py`, the latter including `solve_body_supermodes`' own body
   mesh).

   > ⚠️ **This strands the 39 mode-solve cache entries**, which is correct and
   > intended: they were computed at res 40 and `"resolution"` is in their
   > key, so they now miss rather than being served for a run that solves a
   > different mesh. The 6 FDTD entries are unaffected by this change (they
   > key off `FDTD_RESOLUTION`, still 40) — though they remain stranded by the
   > earlier axis fix, and are now *additionally* stranded by the M3 material
   > fingerprint added the same day.
   >
   > **`FDTD_RESOLUTION` is now the open half of gap 3.** Nothing has measured
   > the FDTD side's convergence on this device. That is a real remaining
   > to-do, not a closed question — it is simply not the one that was
   > producing a quarter-FSR error.

   Original gap text, for the record: never measured, for either engine. The mode solve
   feeds `neff(λ)` into three straight-propagation phases *and* the reference-plane
   correction on `b(λ)`, so an n_eff error propagates into resonance position
   directly. Given `directional_coupler_validation.md` item 11 found the shared res=40 overstates
   Δn_eff by 7% on a coupled pair, the bare-strip n_eff here needs its own
   convergence check before any resonance wavelength is trusted.
4. ~~**Full-device FDTD reference** — **none exists**, and there is no
   `ring_resonator_fdtd.py` to build one from.~~
   **✅ WRITTEN 2026-08-03 — `examples/ring_resonator/python/ring_resonator_fdtd.py`,
   running as job 11655143 (2-D, res 20/30/40).**

   Bus + complete racetrack in one meep run, normalized against the same cell
   with the racetrack deleted, so source spectrum, mode profile, bus
   propagation and both PML interfaces divide out. No decomposition, no EME, no
   network algebra, no reference planes — which is the point: the monolithic
   and cascaded routes are two ways of *assembling* this device, and no
   comparison between them can say which assembly is right.

   **The cost was measured, and it kills 3-D.** Driver is Q, not mesh. Q=14294
   → τ = Q/(2πf) = **3526 meep time units**, and `stop_when_dft_decayed(1e-9)`
   stops near **26 τ**. Scaling job 11653509_3's measured 0.060 s/step by the
   4.00× cell-volume ratio and by res³:

   | res | 2-D @10 τ | 2-D @30 τ | 3-D @10 τ | 3-D @30 τ |
   |---|---|---|---|---|
   | 20 | 0.07 h | 0.21 h | 3 h | 10 h |
   | 30 | 0.24 h | 0.71 h | 18 h | 53 h |
   | 40 | 0.56 h | 1.67 h | 56 h | **167 h** |

   So **3-D full-ring is out of reach of the 10 h wall at every usable
   resolution** — it needs a multi-day allocation or checkpoint/restart. That
   is a finding about this check, not a way around it, and the script's verdict
   prints the 2-D scope rather than implying 3-D coverage.

   ⚠ **My first cost estimate was wrong and is recorded rather than replaced.**
   It said 45 h for 3-D res 40, from assuming T = 8 τ. 8 τ leaves 1.8% residual
   amplitude and the decay criterion actually used stops near 26 τ. The real
   figure is 167 h — 3.7× worse. The 45 h number had already been used to
   justify defaulting to 2-D; the corrected one justifies it far more strongly.

   **2-D still answers the question in play**, because both hybrid scripts have
   a 2-D path (`RING_TEST_MODE=1`), so the comparison runs at *matched* physics
   instead of conflating decomposition error with a 2-D-vs-3-D material
   difference. It does **not** validate the 3-D numbers the device ships, and
   the verdict says so.

   The array sweeps resolution for the same reason check 3 existed: a
   full-device reference that is itself mesh-limited would "disagree" with the
   hybrid for a reason unrelated to the decomposition.

   Two guards, both exercised on a deliberately-too-short smoke run (res 8, cap
   600 = 0.2 τ): a **truncation detector** that refuses to let linewidths be
   quoted from a capped run, and a **τ check** that fires when the device run
   decays too fast for the assumed Q. The smoke run produced 32 spurious dips
   and FSR 6.5 nm against the predicted 12.87 — and was correctly flagged
   rather than reported. `assert_geometry_closes()` re-derives the racetrack
   perimeter from the polygons the file actually draws (43.4159 µm) and checks
   the far-bend corners meet their straights, so a drifted mirrored constant
   fails loudly instead of shifting the FSR.

~~Also unchecked: the interpolation strategy. FDTD coefficients are interpolated
magnitude + unwrapped-phase, and n_eff by a degree-2 polynomial fit, onto a dense
grid. Both choices are argued in the docstring and both are load-bearing for
resonance position — a piecewise-linear n_eff would put an artificial kink at
every solve node — but neither has a test.~~

**✅ Half-closed 2026-08-03.** `interp_complex` is now `beef/spectral.py` with 11
tests (`tests/test_spectral.py`), and both ring scripts import it — a test
asserts bit-identical output against the removed inline definition, so promoting
it changed no result. The mag/phase-vs-real/imag claim is now *measured* on an
analytic ramp (exact to `<1e-12` vs 4+ orders worse) rather than argued. The
strategy's own failure mode is pinned too, including a **false negative I did not
expect**: `unwrap_is_safe()` measures the post-unwrap step, which `np.unwrap`
forces into `(-π, π]`, so on already-aliased data it reports safe while the
interpolant is badly wrong with a perfectly flat magnitude. That limit is
fundamental — above Nyquist the aliased ramp is self-consistent even under
decimation — and is now documented and tested rather than assumed away.

**✅ Fully closed 2026-08-03 — the degree-2 `n_eff` polyfit is measured too.**
`examples/verification/python/verify_ring_neff_polyfit.py`. No new simulation:
it reads the 26 bare-strip and 13 coupled-pair `n_eff` values already in
`example_sim_cache/` at res 100, and does **leave-one-out** — drop an interior
wavelength, fit the rest, predict the dropped point, compare against what the
mode solver returned there. A held-out prediction, not a goodness-of-fit on the
training points, because a low-order fit can hug its own nodes and still be
wrong between them, which is the failure in question.

Every number is converted into the currency check 3 uses, `Δφ = 2π·Δn_eff·L_outer/λ`
as a fraction of one FSR, so it is directly comparable rather than a residual in
isolation:

| degree | LOO max \|Δn_eff\| | % of one FSR (strip) | % of one FSR (even supermode) |
|---|---|---|---|
| 1 | 5.99e-05 | 0.0271% | 0.0642% |
| **2 (shipped)** | **9.83e-06** | **0.0044%** | **0.0055%** |
| 3 | 5.28e-07 | 0.0002% | 0.0002% |
| 4 | 6.27e-09 | 0.0000% | 0.0000% |

Production-like — fit degree 2 on 13 points, residual against all 26 including
the ones the fit never saw — gives **0.0033% of an FSR**.

**Verdict: the fit is not a limiting error source.** 0.0044% against the 12.5%
of an FSR the shipped FDTD resolution contributes is a factor of **2817**.
Raising the degree would change no resonance position anyone could measure.
*Even degree 1* comes in at 0.0271%, so degree 2 has real headroom — the choice
was safe, and is now known to be safe rather than argued to be.

⚠ Scope: this validates the **fit**, not the `n_eff` values fitted to. Those are
check 3's business.

With this, **both halves of the interpolation gap are closed** — the coefficient
resampler by measurement against an analytic function, the `n_eff` fit by
held-out prediction against cached solves.

## Order

~~**Fix item 10 first**~~ — done 2026-07-26, see above. The code is fixed; the
**re-run is not done** and is now the top item.

0. ~~**NEW, 2026-07-31 — fix the y-normal excite normalization first.**~~
   **✅ CLOSED 2026-08-01 (11616920/11616921) — there was nothing to fix.**
   The suspected ~27% coupling error was second-order discretization, not a
   normalization bug: fitted order p = 1.988, extrapolating to 0.979 at
   `FDTD_RESOLUTION = 40`. `ring_resonator_hybrid.py`'s hand-rolled `k2`/
   `t_ring` can be trusted from this path, so this no longer sits ahead of
   the re-run. Full evidence in the boxed section under prerequisite (a) and
   in the prerequisite list's own item 0 — **this line contradicted both
   until 2026-08-02.**
1. ~~**Re-run both scripts** and check the centre wavelength does not move.~~
   **✅ DONE 2026-08-02 — job 11634900, both tasks COMPLETED (3-D, res 40 /
   mode 100).** First trustworthy resonance positions this device has
   produced. Full results in `running_jobs.md`'s nineteenth snapshot; the
   headline is that **the two scripts now disagree**, see the boxed entry
   below.
2. ~~**Gap 3 (resolution convergence)** — cheap, and it gates every resonance
   number.~~ **✅ DONE 2026-08-01 (11616610).** ~~shipped res 40 is a
   QUARTER-FSR off res 100~~ → **CORRECTED 2026-08-03: an EIGHTH of an FSR
   (12.6%), not a quarter.** The convergence script's lever arm was 2× too
   large; job 11634859 printed the true `L_outer = 7.0000 µm`. Full corrected
   table in gap 3. **The decision this forces is unchanged in kind but smaller
   in size**: re-running at res 40 still regenerates resonance positions
   materially set by the mesh, so decide the resolution first. The MMIs' "hold
   40 and carry it as a caveat" precedent is now **closer to transferable than
   it looked** — an eighth-FSR slide is still not a caveat for the one quantity
   this device exists to report, but the gap between "caveat" and "wrong answer"
   narrowed by 2×. Re-read that trade-off before spending res-100 node-hours.
3. Then the (a)/(b)/(c) work below, after which gaps 1, 2 and 4 become possible
   in the same form as the other devices.
