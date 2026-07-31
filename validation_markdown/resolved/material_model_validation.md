# Material model validation

**Cross-device.** The other files in this directory are one-per-example; this one is not,
because the material model is shared by every device and by both engines. It is the
validation record for `directional_coupler_validation.md` **item 8** and `notes.md`
recommendations **#5** and **#9**, which are three views of the same work.

Scope: what the FDTD and the EME each believe silicon and oxide are at a given wavelength,
whether they agree, and whether either is right.

## Status — 2026-07-28

| step | what it does | state |
|---|---|---|
| 0 | resolved indices reach the cache fingerprint | ✅ done |
| 1 | one shared analytic model, both engines | ✅ done — agree to 4e-16 |
| 2 | FDTD medium actually dispersive across the band | ✅ **done 2026-07-31** — implemented, gate PASSED on all 4 stack/band combinations (worst Δn_eff −0.00141 on SOI_160/O-band, tol 0.002); item 4's resolution sweeps closed the same day (job 11588997 — the pin holds at the tolerance production uses). **Opt-in** (`dispersive=True`, default off) is a deliberate design choice now, not a pending blocker — turning it on for any example is separate, unstarted work |

**Steps 0 and 1 are closed.** Both engines evaluate one Sellmeier model
(`beef/materials/stack.py::_SELLMEIER` — Salzberg 1957 Si, Malitson 1965 fused silica,
lifted from meep's own library and verified against it to 4e-16). The two independent
`np.argmin` nearest-sample lookups that preceded it are gone, along with the 0.0063 step
they put in the middle of the C-band. Full reasoning, including the three rejected
alternatives, is in `directional_coupler_validation.md` item 8.

Adopting Salzberg shifts n_Si by ~0.002 against meow's old table. **All four goldens are
stale because of it**, and being `@pytest.mark.slow` they have not been re-run.

**Step 2 is closed as of 2026-07-31** — the re-parameterisation was implemented, gate-passed,
and its resolution-sensitivity question (work item 4) answered, all measured rather than
assumed; see below.

## Why step 2 was reverted

Building both media from the Sellmeier susceptibilities was implemented and backed out.
meep subpixel-averages only the *instantaneous* ε and assigns susceptibility terms per
pixel unsmoothed. The Sellmeier form — like meep's own library materials — puts ε∞ = 1
with the entire index in susceptibilities, so meep builds a structure whose
`epsilon diagonal = (1,1,1)`: no dielectric contrast at the Si/SiO2 interface, and
subpixel smoothing does nothing at all.

| 450×220 nm Si strip in oxide, res 40, 1.55 µm | n_eff |
|---|---|
| non-dispersive `mp.Medium(index=n)` | 2.313398 |
| dispersive, ε∞ = 1 + Lorentzians | 2.363038 |
| | **+0.049640** |

That is a 0.0496 discretization loss in exchange for the ~0.002 the dispersion was worth —
**25× the wrong way**, on the same quantity. `notes.md` #9 predicted it before the attempt
and is kept as the standing reason the media are still `mp.Medium(index=...)`.

## Is the documented fix actually realizable?

Item 8 and `notes.md` #5 both prescribe the same remedy: *re-parameterise so ε∞ carries the
bulk of the permittivity and the susceptibilities carry only its wavelength dependence.*
That was written as an assumption. It has now been checked, and it holds — but not in the
form the prescription implies.

**The exact split is not expressible in meep.** The identity is clean,

    σ·λ²/(λ² − λ₀²)  =  σ  +  σ·λ₀²/(λ² − λ₀²)

and the constant `σ` is exactly what belongs in ε∞. But the remainder is `σω²/(ω₀² − ω²)`,
and meep's `LorentzianSusceptibility` is `σω₀²/(ω₀² − ω²)` — an ω₀² numerator where an ω²
is needed. `DrudeSusceptibility` does not have the right shape either. So ε∞ cannot be
peeled off analytically; the model has to be **re-fit**.

Two constraints fall out of doing that, neither anticipated in the existing docs.

**Finding 1 — it cannot be fit once, globally.** Silicon's telecom-band dispersion comes
from a UV resonance at 0.3015 µm, and no combination of far-IR Lorentzians reproduces that
shape across the full range the examples use. Best achievable, fitting ε over the band and
reporting the resulting index error:

| fit band | max abs Δn vs Salzberg (Si) |
|---|---|
| 1.50–1.60 µm | **1.9e-4** |
| 1.26–1.36 µm | 4.7e-4 |
| 1.26–1.70 µm | 5.5e-3 |

The full-range fit is *worse than the error being fixed* and is not usable. Per-band it is
comfortable. This is compatible with how the solver already works — `_build_domain` builds
the media once per run from that run's own wavelengths — but it does change the media from
a global constant into a per-run quantity.

**Finding 2 — λ₀ is degenerate and must be pinned.** For λ₀ ≫ λ only the ratio σ/λ₀²
matters, so a free fit walks a ridge to its bounds (σ = 500 at λ₀ = 52 µm) with no gain.
Pinning λ₀ and solving for ε∞ and σ gives *identical* fit quality at λ₀ = 5, 10, 20 and
50 µm, and makes the model **linear in both remaining parameters** — a deterministic
least-squares solve, no optimiser in the build path.

> ⚠️ **λ₀ = 10 below is superseded. The shipped pole is now DERIVED per band —
> 1.15 × the band's upper edge — not a constant at all.** Everything in the rest of this
> section is correct as *fit* analysis, and the "pin it, don't fit it" conclusion stands.
> But λ₀ = 10 **fails the strip gate**: fit quality is not what λ₀ actually controls — it
> controls how the permittivity splits between the smoothed ε∞ and the unsmoothed
> susceptibility. The 104% contrast figure quoted below is precisely the symptom that
> failed. A single global value cannot work at all; see **GATE PASSED** below.

With λ₀ = 10 µm over the coupler's own 1.50–1.60 µm sweep:

| | ε∞ | σ | max abs Δn over band |
|---|---|---|---|
| Si | 12.5288 | +17.611 | 1.9e-4 |
| SiO2 | 2.1114 | +1.063 | 3.6e-6 |

The number that decides it: the ε∞ contrast is 10.417 against a true contrast of 10.009 at
1.55 µm — **104% of the contrast is available to subpixel averaging, against 0% today.**
The residual unsmoothed susceptibility is 3.7% of ε_Si. Fit error is 2.3% of the 0.00825
that n_Si genuinely moves across that band, so the parameterisation keeps ~98% of the
dispersion it is there to capture.

*Reproducing this:* evaluate `_SELLMEIER` on a grid across the band, then solve
`ε(λ) ≈ ε∞ + σ·λ²/(λ² − λ₀²)` for `[ε∞, σ]` by linear least squares at fixed λ₀.
~~Not yet committed as a script — see work item 1.~~ **Committed 2026-07-31** as
`beef.materials.stack.reparam_lorentz()`, with
`examples/verification/python/verify_item_8_reparameterisation.py` emitting every table
above. All figures reproduce (ε∞ 12.5287 / σ +17.6106 / max Δn 1.96e-4 for Si; 2.1114 /
+1.0635 / 3.71e-6 for SiO2; contrast 10.4174 vs true 10.0094 = 104.1%).

**Two corrections from running it.** The O-band fit error is **4.9e-4**, not 4.7e-4, and
the full-range is **5.6e-3**, not 5.5e-3 — grid-density artifacts of the original
hand-computation, conclusions unchanged. And λ₀ is *not* strictly degenerate: max |Δn|
runs 2.11e-4 / 1.96e-4 / 1.93e-4 / 1.92e-4 at λ₀ = 5 / 10 / 20 / 50 µm, improving
slightly and converging as λ₀ grows rather than being identical. The operative point is
unaffected — σ swings 3.78 → 461 across that range, so σ/λ₀² is what is determined and λ₀
must still be pinned, not fitted.

**Why the re-fit works where the algebraic split cannot** (worth stating plainly, since
the doc above only says the split fails): meep's `LorentzianSusceptibility`
`σ·f₀²/(f₀² − f²)` is, in wavelength, exactly `σ·λ²/(λ² − λ₀²)` — the same form as a
Sellmeier term, which is why `meep_medium()` can represent the model exactly today. What
is *not* expressible is the **remainder** after peeling ε∞ off analytically, which carries
a λ₀² numerator. Re-fitting directly into `ε∞ + σ·λ²/(λ² − λ₀²)` sidesteps that: the
fitted form is one meep Lorentzian by construction. Verified end to end — meep's own
`Medium.epsilon()` tracks Sellmeier across the C-band to 1.96e-4, exactly the fit error.

## Work

1. ~~**Commit the fit.**~~ **DONE 2026-07-31.** `reparam_lorentz()` in
   `beef/materials/stack.py` (shipped, so the medium builder and the verification script
   share one implementation), with
   `examples/verification/python/verify_item_8_reparameterisation.py` emitting the tables.
   Pure numpy — imports neither meep nor meow, runs in milliseconds on a login node.

2. ~~**Per-band media.**~~ **DONE 2026-07-31**, opt-in and inert by default.
   `MaterialStack(dispersive=False)` is the new constructor flag;
   `si_medium(wl)`/`sio2_medium(wl)` now accept a scalar *or* the run's whole wavelength
   array, and build `meep_medium_reparam()` only when the flag is on **and** a real band
   was passed.

   **The doc's "three call sites at lines 526, 776, 1175" was wrong on both counts.** The
   real sites are `beef/solvers/fdtd.py` 615/616, 930/931, 1581/1582 — and only **two**
   should change. 615 (`_build_domain`) and 1581 (`solve_section`) are broadband and now
   receive `wavelengths`; **930 is a single-wavelength mode solve**, where the frozen index
   at that exact wavelength is already the correct medium and dispersion is meaningless.
   Left scalar deliberately.

   Two further callers are off the library path and unchanged:
   `examples/ring_resonator/python/ring_resonator_hybrid{,_cascaded}.py` call
   `si_medium(wl_center)` directly — the same "ring bypasses the solvers" problem
   `ring_resonator_validation.md` is about.

   Default is **False** because the strip-test gate below has not run. Switchable so the
   gate is a clean A/B, matching `guard_duplicate_modes` / `calibrate_body_neff`.

3. ~~**Fingerprint the parameterisation — this one bites.**~~ **DONE 2026-07-31.**
   `resolved_material_config()` now adds `resolved_dispersive` and
   `resolved_{si,sio2}_reparam = [ε∞, σ, λ₀]` whenever the stack is dispersive over a real
   band. Confirmed by test that flipping the flag moves the key, and that two different
   bands give different keys (the media are per-run now, so they must not collide). The
   fitted numbers are recorded **directly** rather than via a model-version tag, so any
   change to band, pinned λ₀, sample count, or the Sellmeier coefficients underneath moves
   the key on its own with nothing to remember to bump. The non-dispersive key is
   byte-unchanged, so everything already cached keeps matching.

4. **Re-run the resolution sweeps.** `notes.md` #9 already budgets this as part of the
   change, not as a follow-up: smoothing moves the mesh-convergence point, so item 11's
   res-80 EME pin is not inherited across this change.
   (`batch/archive/run_directional_coupler_*_resolution_sweep.batch` -- moved
   to archive 2026-07-31, retired pending new sweeping functionality, but
   still the literal scripts this item names; revive or replace them when
   this item is actually run.)

   **Started 2026-07-31 with the cheap half first — job 11586986.**
   `examples/verification/python/verify_item_8_dispersive_convergence.py`. The literal
   instruction is 11 array tasks × 2 nodes × 24 ranks at up to 24 h each, which answers the
   question by brute force. But **item 11's pin was justified on Δβ**, the even/odd
   supermode difference, and Δβ is a *mode-solve* quantity — no time-stepping, no
   propagation, no S-matrix. So "does the convergence point move" is answerable on **one
   serial task** (~2 h, measured ~40 s per point at res 20, growing ~res²), and that runs
   first to decide whether the node-days are warranted and over what range.

   > ### ⚠️ Job 11586986 measured nothing. This item is NOT started.
   >
   > **Corrected 2026-07-31.** That job's "dispersive" column came back
   > **byte-identical** to its non-dispersive one — same n_eff, same Δn_eff, same Δβ to
   > six decimals at res 20/30/40/60/80. Not a small effect: the flag never reached the
   > calculation.
   >
   > `FDTDSolver._build_section_mode_sim()` built its media with `si_medium(wl_center)`
   > **directly, bypassing `_media()`**, so it never read `self.dispersive`. The
   > justification was work item 2's own — *"930 is a single-wavelength mode solve, where
   > the frozen index at that exact wavelength is already the correct medium"* — and that
   > is true for `section_mode_fields()`, which takes one wavelength. It is **false for
   > `section_mode_neff()`**, which accepts a whole band and loops it on one structure. For
   > a band caller the frozen centre index is the exact defect work item 2 fixed for
   > `_build_domain()`/`solve_section()`, left in place by an argument that did not apply
   > to it.
   >
   > An earlier revision of `running_jobs.md` blamed a stale λ₀ (the job predates λ₀
   > becoming band-derived in `4296ee7`). That was also wrong: λ₀ cannot matter when the
   > dispersive medium is never built.
   >
   > **Fixed the same day.** `_build_section_mode_sim()` gained a `wavelengths=` parameter
   > and routes media through `_media()`; `section_mode_neff()` forwards its band. Scalar
   > callers are byte-unchanged and not by convention — `MaterialStack._medium()` requires
   > **both** the flag and a real band (`size > 1 and min < max`), so a scalar cannot reach
   > the dispersive path. Measured on a real solve after the fix:
   >
   > | | max \|Δn_eff\| over 1.50–1.60 µm | at centre 1.55 | at edge 1.60 |
   > |---|---|---|---|
   > | before | **0.000e+00** | 0 | 0 |
   > | after | **5.29e-03** | 5.6e-05 | 4.6e-03 |
   >
   > Agreement at band centre and divergence at the edges with opposite signs is the
   > signature of unfreezing a centre-pinned index, which is what confirms the fix does the
   > intended thing rather than merely *something*.
   >
   > **9 new tests** (`tests/test_fdtd_mode_sim_dispersive.py`), meep-free, pinning both
   > halves: a band caller under `dispersive=True` gets per-wavelength media, every scalar
   > caller is unchanged in either configuration, and — the one that would have caught
   > this — `section_mode_neff()` forwards its whole band rather than just the centre.
   > Both key tests fail against the pre-fix code; the seven unchanged-behaviour tests pass
   > against it, which is what shows the fix is confined to the band path.
   >
   > **The script now refuses to run vacuously.** It solves once per setting at the lowest
   > resolution first and **aborts with a non-zero exit** if the two agree exactly, so an
   > inert independent variable can never again be recorded as a measurement. It also
   > reports Δn_eff at a band **edge** as well as at centre, because dispersion is
   > near-inert at centre by construction — a smoke run shows +1.7e-5 at centre against
   > +5.9e-4 at the edge, so a centre-only reading understates it ~34×.
   >
   > **Salvaged from 11586986: its non-dispersive column is valid**, and it is the only
   > half that ever meant anything. Coupler body, SOI_220, band centre 1.55 µm, Δn_eff of
   > the even/odd pair:
   >
   > | res | n_eff[0] | n_eff[1] | Δn_eff | step in Δn | Δβ | sec |
   > |---|---|---|---|---|---|---|
   > | 20 | 2.243923 | 2.216122 | 0.027802 | — | 0.11270 | 35.9 |
   > | 30 | 2.301684 | 2.277158 | 0.024526 | −3.28e-03 | 0.09942 | 85.1 |
   > | 40 | 2.323677 | 2.300177 | 0.023501 | −1.03e-03 | 0.09526 | 137.1 |
   > | 60 | 2.335676 | 2.312802 | 0.022874 | −6.27e-04 | 0.09272 | 334.1 |
   > | 80 | 2.350117 | 2.328029 | 0.022089 | −7.85e-04 | 0.08954 | 618.2 |
   > | 100 | 2.351338 | 2.329235 | 0.022104 | **+1.50e-05** | 0.08960 | 886.1 |
   > | 120 | 2.355171 | 2.333476 | 0.021695 | **−4.08e-04** | 0.08795 | 1392.2 |
   >
   > **Worth reading on its own, independent of dispersion: Δβ is not converged at res
   > 120, and res 100 looks converged only by accident.** The step collapses to +1.5e-05 at
   > 100 — a factor of 50 below its neighbours, and the kind of number a sweep stopping
   > there would report as a plateau — then reopens to −4.08e-04 at 120, an order of
   > magnitude larger. That is `README.md` Finding 4 ("a converged-looking number can still
   > be wrong, and a resolution sweep will not tell you") appearing on the very quantity
   > item 11's res-80 pin was set on, in a sweep run to check that pin. Not acted on here —
   > it is `resolved/directional_coupler_validation.md` item 11's territory, and one
   > non-monotone step is not yet a case — but any re-derivation of that pin should sweep
   > past 120 rather than stopping at the first small step.

   Sweeps res 20→120 on the coupler body cross-section, both material settings, reporting
   the *step* between consecutive resolutions rather than the absolute value — the same
   criterion item 11 and `mode_spectrum_convergence.py` use.

   > ### ✅ MEASURED 2026-07-31 — job 11588997 (1:56:53). The pin holds where production
   > uses it; it does not hold at every threshold.
   >
   > First resolution whose step drops under each threshold, non-dispersive vs dispersive:
   >
   > | threshold | centre — non-disp | centre — disp | edge — non-disp | edge — disp |
   > |---|---|---|---|---|
   > | 1e-3 | 60 | 60 | 60 | **80** |
   > | 3e-4 | 100 | 100 | 100 | 100 |
   > | 1e-4 | 100 | 100 | 100 | 100 |
   > | 3e-5 | 100 | 100 | 100 | **not reached by 120** |
   >
   > **At band centre, dispersion changes nothing** — identical convergence resolution at
   > every threshold, the clean "unchanged" case this item's setup anticipated. The edge is
   > where it has to be judged instead (dispersion is near-inert at centre by construction,
   > per the script's own read-me), and there it is not free: the 1e-3 threshold moves from
   > res 60 to res 80 under dispersion.
   >
   > **That move does not cost anything production doesn't already have.** Item 11 already
   > pins the shipped EME resolution at 80, for an unrelated reason (Δβ on this exact
   > cross-section). The dispersive case's actual requirement, measured here, lands exactly
   > on the resolution already shipped — so while a literal reading of the two-branch
   > criterion above ("moved ⇒ full sweep required") is triggered by the 60→80 change, there
   > is no gap between what this measures and what is already running: nothing currently
   > deployed was ever simulated at res 60.
   >
   > **What does not close:** at the tightest threshold (3e-5), non-dispersive converges
   > everywhere by res 100; dispersive converges at centre but **never drops under 3e-5 at
   > the edge inside the tested range (20–120)**. The raw dispersive-minus-non-dispersive
   > Δn_eff table oscillates at ±(3–5)e-4 with no visible trend across resolution — reading
   > like noise-floor disagreement rather than a growing one — but "fully converged at any
   > tolerance" cannot be claimed at the edge under dispersion from this data.
   >
   > **Verdict on the node-day sweep this item originally proposed: not warranted**, at the
   > 1e-3-class tolerance item 11 actually used to set the shipped resolution. Downgrading a
   > global material-model check to "did the specific number production ships change" is
   > exactly what this cheap serial task was for. This closes the mode-solve half of work
   > item 4; it still says nothing about the cascade's *S-matrix* convergence (only the full
   > sweep would), and dispersion itself remains opt-in (`dispersive=False` default).

5. **Add the `HybridSolver` guard from `notes.md` #9** — assert the FDTD and EME sections
   resolve the same index at the same wavelength — as the regression test for this work.
   It needs a **tolerance, not equality**: see the last risk below.

   **DONE 2026-07-31.** `HybridSolver._check_material_index_agreement()`, controlled by
   the new `material_consistency_action` ("warn" default / "raise" / "ignore"), tolerance
   `_MATERIAL_INDEX_TOL = 2e-3`.

   Both of the doc's requirements turned out to be load-bearing:

   - **Compare built media, not index lookups.** `si_index()` is the shared Sellmeier
     model both engines call, so it agrees to machine precision in *every* configuration,
     including the broken one. A guard reading it would be vacuous. `_effective_index()`
     instead evaluates the `mp.Medium` the FDTD actually builds and the per-wavelength
     value the EME actually solves at. A test pins this specifically.
   - **Tolerance, not equality.** With dispersion on the FDTD runs a *fit* while the EME
     stays exact, so they legitimately differ by ~1.3e-3. Equality would fail the correct
     configuration. 2e-3 sits above that and far below the 7.9e-3 a frozen medium is wrong
     by at the O-band edge.

   Measured on a real chain (FDTD + EME, O-band):

   | configuration | worst index gap | verdict |
   |---|---|---|
   | **today's default** (frozen FDTD, exact EME) | **0.007862** | over tolerance |
   | both dispersive | 0.001349 | passes |
   | both frozen | 0.000000 | passes |

   Default is "warn" rather than "raise" precisely because the shipped default does not
   satisfy it — raising would break every existing hybrid chain over a pre-existing
   condition. Best-effort: solvers or stacks that cannot be interrogated are skipped, not
   guessed at (caught by the existing `_FakeStack` tests, which this initially crashed).

   ~~**Cheap half done 2026-07-31, warning-only.**~~ `HybridSolver._check_material_consistency()`
   runs on every `solve()`/`solve_length_sweep()` and logs a WARNING when the chain's
   sections disagree about dispersion, naming which sections are on which side and the
   worst n_Si gap across the band. The hard assert still waits on the gate, which is what
   decides the tolerance.

   **It fires on today's default configuration, and that is correct.** A broadband
   `FDTDSolver` freezes its index at band centre while an `EMESolver` body re-solves per
   wavelength, so the chain has *always* contracted modes solved in two different
   materials — 0.00786 apart in n_Si at the O-band edges. That mismatch predates this work
   entirely; it was simply never checked.

## Verification

**The gate is the strip test**, and it is cheap: re-measure the 450×220 nm strip at res 40,
1.55 µm against the 2.313398 non-dispersive baseline already recorded above. Acceptance is
**|Δn_eff| ≲ 0.002** — below what the dispersion buys. The reverted attempt scored 0.0496
on exactly this test, so it is a like-for-like comparison against a known failure.

Work items 4 and 5 should not be queued until that gate passes; item 4 is the only
expensive part of this change.

### ✅ GATE PASSED — 2026-07-31, and it moved λ₀

`examples/verification/python/verify_item_8_strip_gate.py`, run on a login node in ~4 s.
Three arms, because two cannot establish that the harness reproduces the original problem:

| arm | n_eff | Δ vs A | ε∞(Si) actually smoothed |
|---|---|---|---|
| A non-dispersive (shipped) | 2.311707 | — | 12.0946 |
| C legacy ε∞=1 (the reverted attempt) | 2.361645 | **+0.049938** | **1.0000** |
| B re-parameterised | 2.311019 | **−0.000688** ✅ | 12.2714 |

**C reproduces the recorded +0.049640 to 3e-4**, so the harness measures what the original
attempt measured and B's number is evidence about the fix rather than about the test. The
ε∞ column is the mechanism, stated numerically for the first time: the reverted attempt
handed meep a structure with *no dielectric contrast at all* to average.

Arm A lands 0.0017 below the doc's 2.313398 baseline — this harness is a bare strip in a
2×2 µm cell, not whatever the original ad-hoc measurement used. Read Δ *within* a run, not
absolute n_eff across runs; the C-vs-A reproduction is what validates that.

**λ₀ is now derived per band (1.15 × band max), not a constant — and the gate is what
forced that.** The original gate covered only SOI_220; extending it to SOI_160 (2026-07-31,
same day) **failed** on the O-band at −0.00375, the combination `mmi_1x2` and the modal
multiplexer both use. A 160 nm core confines weakly, so more of the mode sits at the
Si/SiO2 interface where smoothing acts — arm C's failure is ~2× larger there (+0.108 vs
+0.050), and a SOI_220 pass does not transfer.

The conflict is **structural, not tuning**: SOI_160/O-band needs λ₀ ≤ 1.60, while the
C-band's own near-resonance floor is 1.10 × 1.60 = 1.76. No single global λ₀ satisfies
both. Scaling λ₀ to the band removes the conflict and clears the guard by construction:

| case | λ₀ = 1.15·max | λ₀ = 1.25·max |
|---|---|---|
| SOI_220 C-band | −0.00033 ✅ | −0.00069 ✅ |
| SOI_220 O-band | −0.00010 ✅ | −0.00060 ✅ |
| SOI_160 C-band | −0.00109 ✅ | −0.00172 ✅ |
| SOI_160 O-band | **−0.00141 ✅** | **−0.00239 ❌** |

**Historical note — λ₀ 10.0 → 2.0 was the first correction, since superseded by the rule
above.** The reasoning that produced it still stands: At λ₀ = 10 the gate **failed**
(−0.00222 C-band, −0.00258 O-band) — the doc's own risk *"ε∞ contrast is 104% of the true
contrast… whether that biases n_eff, and in which direction, is not predictable from the
fit"* materialising, biasing **low**. λ₀ sets how the permittivity splits between the
smoothed ε∞ and the unsmoothed susceptibility, and only that split decides whether
smoothing still works:

| λ₀ | Δn_eff C-band | Δn_eff O-band | σ (C/O) | fit, % of band swing |
|---|---|---|---|---|
| 10.0 | −0.00222 ❌ | −0.00258 ❌ | 17.6 / 37.7 | 2.4 / 3.3 |
| 3.0 | −0.00160 ✅ | −0.00203 ❌ | 0.89 / 2.30 | 3.1 / 3.9 |
| **2.0** | **−0.00069 ✅** | **−0.00123 ✅** | **0.12 / 0.51** | **5.5 / 5.2** |

2.0 is the only value clearing both bands. **Fit quality is the wrong criterion** — it
improves monotonically with λ₀ over this range while the gate result gets worse. The cost
is ~3% more of the band's dispersion given up (still keeping ~95%).

**This also retires the σ risk below.** σ drops from 17.6 to 0.12; "meep stability is
unverified at σ = 17.6" no longer applies, and no stability warnings appeared.

λ₀ = 2.0 sits inside `_MODEL_RANGE_UM` (1.10–2.10), so `reparam_lorentz()` now **raises**
if a band comes within 10% of the pole rather than returning a near-resonant model. The
widest band the examples span (1.26–1.70) clears it by 18%.

**Items 4 and 5 are unblocked.**

**Sequencing.** This invalidates the four goldens a second time. It must land **before**
the outstanding cold regeneration pass — which already owes the goldens, the
`COUPLING_LENGTH` change and step 1's n_Si shift — or that pass gets paid for twice.

## Risks

- ~~**The ε∞ contrast is 104% of the true contrast, not 100%.**~~ **Measured 2026-07-31 —
  the risk was real.** It biases n_eff **low**, and at λ₀ = 10 that bias (−0.0022 C-band,
  −0.0026 O-band) failed the gate outright. Fixed by moving λ₀ to 2.0, which brings the
  contrast to 101.7% and the bias to −0.0007. This is the one risk on this list that
  actually fired.
- ~~**meep stability is unverified** at σ = 17.6 with f₀ = 0.1.~~ **Retired 2026-07-31.**
  The strip gate moved λ₀ to 2.0, which drops σ to 0.12 (C-band) / 0.51 (O-band) — two
  orders of magnitude down, squarely unremarkable. The gate run logged no meep stability
  warnings.
- **Engine agreement regresses from 4e-16 to ~2e-4.** The EME keeps the exact Sellmeier
  model via `SampledMaterial` and solves per wavelength, so it needs no re-parameterisation;
  only the FDTD is fit. That is a real cost and the reason work item 5 needs a tolerance.
  It is worth paying: 2e-4 of disagreement *across the band* replaces a band that is
  currently frozen at one wavelength and wrong by up to 0.008 at its edges.

  **Sharpened 2026-07-31 — the two numbers measure different things, and on the one that
  decides correctness this is an improvement, not a regression.** The 4e-16 is agreement
  between the two engines' *index models*: both call `sellmeier_index()`, so they return
  the same number at any wavelength you ask. It is not agreement between the media they
  actually simulate with. The FDTD then collapses its model to a single wavelength when it
  builds `mp.Medium`, while the EME evaluates per wavelength — so the *simulated* materials
  already differ today, by up to 0.00786 in n_Si across the O-band. Measured:

  | | FDTD medium vs EME's exact n(λ), worst over O-band |
  |---|---|
  | today (frozen at centre) | **0.00786** |
  | `dispersive=True` | **4.9e-4** |

  So enabling dispersion improves real FDTD-vs-EME material agreement ~16×. What genuinely
  regresses is model-level agreement — the FDTD's fit differs from exact Sellmeier by
  ~2e-4 where before it differed by 0. Work item 5's guard must therefore compare **built
  media, not index lookups**: a guard written against `si_index()` would read 4e-16 in
  every configuration and be vacuous, including the mismatched one shipping today.
