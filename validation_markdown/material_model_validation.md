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
| 2 | FDTD medium actually dispersive across the band | ❌ **open** — attempted, measured, reverted |

**Steps 0 and 1 are closed.** Both engines evaluate one Sellmeier model
(`beef/materials/stack.py::_SELLMEIER` — Salzberg 1957 Si, Malitson 1965 fused silica,
lifted from meep's own library and verified against it to 4e-16). The two independent
`np.argmin` nearest-sample lookups that preceded it are gone, along with the 0.0063 step
they put in the middle of the C-band. Full reasoning, including the three rejected
alternatives, is in `directional_coupler_validation.md` item 8.

Adopting Salzberg shifts n_Si by ~0.002 against meow's old table. **All four goldens are
stale because of it**, and being `@pytest.mark.slow` they have not been re-run.

**Step 2 is open, and the reason is measured rather than assumed** — see below.

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
`ε(λ) ≈ ε∞ + σ·λ²/(λ² − λ₀²)` for `[ε∞, σ]` by linear least squares at fixed λ₀. Not yet
committed as a script — see work item 1.

## Work

1. **Commit the fit.** `examples/verification/python/verify_item_8_reparameterisation.py`,
   emitting the tables above so the choice of λ₀ and the per-band restriction are
   reproducible rather than asserted here.

2. **Per-band media.** `MaterialStack.si_medium(wl)` / `sio2_medium(wl)` take a scalar
   today and must take the band. Three call sites, all of the form
   `si_medium(wl_center)`: `beef/solvers/fdtd.py` lines 526, 776, 1175. `meep_medium()`
   already exists in `materials/stack.py` as the model to start from.

3. **Fingerprint the parameterisation — this one bites.** `resolved_material_config()`
   records three probe *indices* from `sellmeier_index()`, and step 2 does not change that
   function. Turning dispersion on would therefore leave the cache key **identical**, and
   every existing entry would keep matching and be served the old physics — precisely the
   defect step 0 was created to kill, recurring one level up. The band itself is already
   pinned (`resolved_probe_wl` records min/mean/max), so what is missing is ε∞ and σ.
   Record them directly rather than relying on the band plus a model-version tag.

4. **Re-run the resolution sweeps.** `notes.md` #9 already budgets this as part of the
   change, not as a follow-up: smoothing moves the mesh-convergence point, so item 11's
   res-80 EME pin is not inherited across this change.
   (`batch/sweep/run_directional_coupler_*_resolution_sweep.batch`.)

5. **Add the `HybridSolver` guard from `notes.md` #9** — assert the FDTD and EME sections
   resolve the same index at the same wavelength — as the regression test for this work.
   It needs a **tolerance, not equality**: see the last risk below.

## Verification

**The gate is the strip test**, and it is cheap: re-measure the 450×220 nm strip at res 40,
1.55 µm against the 2.313398 non-dispersive baseline already recorded above. Acceptance is
**|Δn_eff| ≲ 0.002** — below what the dispersion buys. The reverted attempt scored 0.0496
on exactly this test, so it is a like-for-like comparison against a known failure.

Work items 4 and 5 should not be queued until that gate passes; item 4 is the only
expensive part of this change.

**Sequencing.** This invalidates the four goldens a second time. It must land **before**
the outstanding cold regeneration pass — which already owes the goldens, the
`COUPLING_LENGTH` change and step 1's n_Si shift — or that pass gets paid for twice.

## Risks

- **The ε∞ contrast is 104% of the true contrast, not 100%.** Smoothing will act on
  slightly too much. Whether that biases n_eff, and in which direction, is not predictable
  from the fit and is what the strip test measures.
- **meep stability is unverified** at σ = 17.6 with f₀ = 0.1. Expected to be fine — ε∞ is
  large and the resonance is far below the band — but meep's stability warnings must be
  read, not assumed silent.
- **Engine agreement regresses from 4e-16 to ~2e-4.** The EME keeps the exact Sellmeier
  model via `SampledMaterial` and solves per wavelength, so it needs no re-parameterisation;
  only the FDTD is fit. That is a real cost and the reason work item 5 needs a tolerance.
  It is worth paying: 2e-4 of disagreement *across the band* replaces a band that is
  currently frozen at one wavelength and wrong by up to 0.008 at its edges.
