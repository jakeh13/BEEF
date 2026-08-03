# Ring resonator validation — to-do

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
devices are: checks 1–3 in [README.md](README.md) all assume
`FDTDSolver`/`EMESolver` objects that this device never constructs — the
existing example scripts are unchanged and still don't call any of (a)/(b)/(c).

## Material model — the ring is now the only device left frozen

Added 2026-08-01, when `dispersive` was flipped on by default across the
library and every example. **The ring did not follow, and it fails silently
rather than loudly.** All four items below are consequences of the same root
cause this doc is already about — the ring hand-rolls raw meep/meow instead of
going through `FDTDSolver`/`EMESolver` — so they are cheap once the
`HybridSolver` recommendation above lands, and awkward before it.

Nothing here is currently *corrupt*: `example_sim_cache/` holds no ring entries
at all, so items M3/M4 were latent mechanism gaps, not live staleness. They will
become live the first time either ring script is run.

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

⚠ **Q in both runs is grid-limited, not physical.** Every FWHM is an exact
multiple of the 2001-point grid spacing (0.05 nm) — 0.15/0.20/0.25/0.30 =
3/4/5/6 steps. The reported Q (5293–10193) is set by **`N_WL_DENSE`**, a knob
independent of both resolutions. In the 501-point `RING_TEST_MODE` smoke every
FWHM was pinned at exactly 0.4000 nm. **If Q is wanted, raise `N_WL_DENSE`** —
it is post-processing on an interpolated grid and costs no new FDTD.

## Gaps against the four checks

1. **Engine agreement** — not measurable today. The mode solve is a direct meow
   call, and the FDTD is raw meep, so there is no shared cross-section object to
   compare. Blocked on the recommendation above.
2. **Mode identity** — not applicable in the current form (no interface mode
   basis is exchanged; the scripts thread scalar coefficients). Becomes check 2
   proper once (c) lands.
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
   Both device scripts compute `L_outer` identically, so no `L_outer` error can
   separate them — and a 2× lever-arm error would move FSR by 2×, while the two
   measured FSRs differ by 2.9%.

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
4. **Full-device FDTD reference** — **none exists**, and there is no
   `ring_resonator_fdtd.py` to build one from. A whole-ring 3-D FDTD is expensive
   but is the only thing that would test the feedback-loop assembly end to end;
   without it the closed-form sum is unchecked against anything.

Also unchecked: the interpolation strategy. FDTD coefficients are interpolated
magnitude + unwrapped-phase, and n_eff by a degree-2 polynomial fit, onto a dense
grid. Both choices are argued in the docstring and both are load-bearing for
resonance position — a piecewise-linear n_eff would put an artificial kink at
every solve node — but neither has a test.

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
