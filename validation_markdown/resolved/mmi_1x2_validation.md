# 1×2 MMI validation — to-do

> **2026-07-28 — the material model changed under every device here.** Both engines
> now evaluate one analytic Sellmeier model (Salzberg Si, Malitson SiO2) instead of
> each doing its own nearest-tabulated lookup, which shifts n_Si by ~0.002 and
> removes a 0.0063 step that sat in the middle of the C-band. Any golden or stored
> result for this device predates that and needs regenerating; see
> `directional_coupler_validation.md` item 8 and, for the dispersive-FDTD
> re-parameterisation (built and gate-passed 2026-07-31, opt-in via
> `dispersive=True`, off by default), `material_model_validation.md`.



Opened 2026-07-26. Gap list against the four checks in [README.md](../README.md).

**Status 2026-07-26:** gap 1 (stale goldens) is **closed** — both were
regenerated and independently confirmed (jobs 11500524, 11500663). Gap 3 (mode
identity) is **run and FAILED at N_MODES=12**: the two engines do not agree on
modes 8 and 9, so `HybridSolver`'s index-*i*-to-index-*i* contraction was not
valid on this device as configured. **This was the only device in the repo
where check 2 failed.**

**Update 2026-07-30:** gap 3 is **fixed and reconfirmed** — `N_MODES` dropped
12→8 in `mmi_1x2_hybrid.py`, re-audited, identity permutation holds
(`[0..7]`, worst guided |Δn_eff| 0.0003). Mode identity no longer blocks this
device; gap 1's goldens are stale again as a result (this is the third
invalidation, see `## Order` step 0) and gap 2 (resolution) is the next open
item. Details in gap 3.

~~**Update 2026-07-31 (late):** gap 5's calibration failure is **diagnosed** —
`calibrate_body_neff=True` fails on this device because of a **π sign flip in
`T_ii`** at 4 of 88 (wavelength, mode) cells (modes 0-1, λ ≤ 1.28); the other
84 fit to ≤0.0059 rad, well inside the 0.1 rad guard.~~ **Superseded
2026-08-01 — there is no sign flip.** The real cause is a **2π turn-count bug
in `beta_from_phases` itself**: it rounded the wrap count independently per
length against one seed slope, so the ~2% MPB-vs-FDTD seed error drifted with
L until a length straddled a rounding boundary. The two readings fit the same
residual pattern exactly and only physics separates them (n_eff 2.783, +1.96%
over MPB, versus the sign-flip reading's 2.101, −23%). **Fixed** —
`beta_fit_from_phases`; validated seed-free across all 88 cells (second
difference an exact integer to 0.00308 everywhere) and on the coupler's cache
(0 of 66 cells change at its production seed). The two *original* hypotheses
(short body / weak SOI_160 confinement) stay **refuted**, now on stronger
grounds. Full evidence in gap 5's boxed section. Calibration reruns in flight
as 11608690/11608692.

**Update 2026-08-01 (late) — the dispersion flip touches FEWER of this device's
results than it first appears, and that is now measured rather than assumed.**

The initial reading of today's material-dispersion default flip
(`material_model_validation.md` step 3) was that every number on this
device had gone stale. **That was wrong for half of them**, and the reason is
worth stating because it applies to every device:

> **A dispersive medium only differs from a frozen one across a *band*.** The
> flip changes `_build_domain()` and `solve_section()`, which are broadband and
> previously froze one index at band centre. It deliberately does **not** reach
> single-wavelength mode solves — `section_mode_neff()` builds its media at a
> scalar wavelength (`fdtd.py`'s own note: *"a single-wavelength mode solve,
> where the frozen index at that exact wavelength is already the correct medium
> and dispersion is meaningless. Left scalar deliberately."*), and a meow
> `Environment(wl=…)` is exact Sellmeier at that wavelength either way.

| gap | what it solves | affected? |
|---|---|---|
| 2 (resolution) | EME `compute_modes()` at a single `WL = 1.31` | **No** |
| 3 (mode identity) | `section_mode_neff()` / `compute_modes()` at 1.31 µm | **No — confirmed by measurement, below** |
| 4 (engine agreement) | `solve_section()` over an 11-point band | **Yes** — ✅ measured dispersive, 11608875 + 11609686 |
| 5 (hybrid vs FDTD) | full-device broadband FDTD + cascade | **Yes** — ✅ re-measured dispersive, 11609215 + 11609687 |

**Gap 3 re-audited under dispersion — job 11608940 (7:13), and every number is
byte-identical to the 2026-07-30 frozen-era table**: bands 0–7 at 2.696912 …
1.705834, correspondence `[0..7]`, 8/8 guided, worst guided |Δn_eff| **0.0003**.
Not "close" — identical to all six printed decimals. That is the deliberate
scalar at `fdtd.py`'s mode-solve path doing exactly what its comment claims,
confirmed on a real device rather than argued from the code.

It also **discharges this gap's own parked note** — the re-audit was to be
re-confirmed *"opportunistically, with the next run of this audit"* to settle
whether the duplicate-mode guards were live. This run is that: guards are in
the shipped code, 8/8 guided, order preserved, and no guard warning appears.

Job 11608941 re-runs gap 2 the same way. On the reasoning above it should also
come back unchanged; it is cheap and left running so that is measured rather
than inferred by analogy from gap 3.

**✅ It did — 11608941 completed 1:03:21 and reproduced all five Δβ values to
every printed digit.** Both gaps are now measured. The general rule, stated
from evidence rather than argument: **a quantity evaluated at a scalar
wavelength is invariant across the dispersion break; only genuinely broadband
quantities move.** That is gaps 4 and 5, and nothing else on this device.

**What this does *not* excuse.** Gaps 4 and 5 genuinely are broadband and do
move. Gap 5's calibrated re-run in flight (11608690/11608692) started **before**
the flip, so it lands frozen-era: valid as a test of the `beta_from_phases` fix,
which is what it exists for, but **gap 5 still needs re-measuring dispersive
afterwards.** Say which side of the flip a broadband number came from.

**✅ Both broadband gaps have now been measured dispersive, 2026-08-01
(evening).** Gap 4 landed for the **first time on this device** (11608875 +
11609686): mean Δβ ratio **0.9946**, a 0.54% monotone bias — all four checks
are now ticked on the 1×2. Gap 5 re-measured (11609215 + 11609687):
**9.11% → 10.18%**, i.e. it got marginally *worse* under the more correct
material model, with the mean offset essentially unmoved (−0.025304 →
−0.025413). The one number still owed on the frozen side of the break is the
**calibrated** gap 5 (7.30%, frozen, offline-derived). Full results in gaps 4
and 5 below.

**That dispersive re-measure is submitted — jobs 11609215 `[1-2]` +
11609216**, the 1×2's half of the pair whose 2×2 twin (11609211/11609212) is
the P4 test (`mmi_2x2_validation.md` gap 5). Here it is *not* a test of a
hypothesis, only the number the shipped configuration produces. The frozen-era
result was preserved first as
`example_output/verification/mmi_1x2_hybrid_vs_fdtd/comparison_FROZEN_pre20260801.npz`
— the comparison script overwrites `comparison.npz`, so without that copy the
before/after pair would not exist.

**Update 2026-07-31:** gap 2 closed (job 11591224) — converges at **res 120**,
one tier above the 2×2's res 100, non-monotone unlike the 2×2's clean decay.
Decision: shipped `RESOLUTION` stays **40**; the gap is carried as a known
mesh-convergence caveat rather than closed by a config bump (real cost jump,
existing budgets sized for 40 — see gap 2 for the reasoning). ~~**Only gap 4
(check 1, engine agreement) is now fully unstarted** for this device — unlike
the 2×2, which at least has it in flight (`mmi_2x2_validation.md` gap 4,
job 11591249).~~ **Superseded 2026-08-01: gap 4 is measured and closed
(11608875 + 11609686, mean ratio 0.9946). No check is unstarted on this device
any more.** Goldens (gap 1) stay blocked behind gap 4 plus the coupler's
own cold regeneration pass, launched today as 11592462/11592463 — see
`running_jobs.md`.

`examples/mmi_1x2/python/` — `mmi_1x2_eme.py`, `mmi_1x2_fdtd.py`,
`mmi_1x2_hybrid.py`. Best-covered example after the directional coupler: it is
the only device with all three golden fixtures (`mmi_1x2_eme.json`,
`mmi_1x2_fdtd.json`, `mmi_1x2_hybrid.json`).

## What it already has

- **A full-device FDTD script** (`mmi_1x2_fdtd.py`) — check 4's prerequisite, and
  the only other example besides the coupler and the 2×2 MMI to have one.
- **Three goldens**, so EME, FDTD and hybrid each have a recorded reference. That
  is stronger than any other device here.

## Gaps

1. **Two of the three goldens are stale — confirmed, not suspected.**

   | golden | last regenerated | affected by the axis bug? |
   |---|---|---|
   | `mmi_1x2_eme.json` | 2026-07-16 (`2220ee1`) | no — pure EME path |
   | `mmi_1x2_fdtd.json` | 2026-07-16 (`2220ee1`) | **yes** — `FDTDSolver` |
   | `mmi_1x2_hybrid.json` | 2026-07-16 (`ed70802`) | **yes** — FDTD sections |

   `directional_coupler_validation.md` items 1–2 landed 2026-07-25, so both FDTD-path goldens
   were recorded with the reversed, half-span wavelength axis — exactly the
   state `directional_coupler_hybrid.json` was in before item 9 regenerated it.
   The `--runslow` golden tests for those two should be failing right now.

   **The stored values prove it without needing to re-run anything.** On the
   same device, at the same three wavelengths:

   | | 1.26 | 1.31 | 1.36 | direction |
   |---|---|---|---|---|
   | `mmi_1x2_eme` (unaffected) | 0.4264 | 0.3863 | 0.3524 | **falling** |
   | `mmi_1x2_fdtd` (as stored) | 0.3752 | 0.3943 | 0.4295 | **rising** |

   The two engines disagree about the *sign* of the dispersion on one device,
   which no amount of engine difference explains. Three things confirm the
   diagnosis rather than merely suggesting it:

   - Reverse the FDTD row and it falls, matching EME: 0.4295 / 0.3943 / 0.3752.
   - Its true axis is 1.3355 / 1.3100 / 1.2854 — a narrower band than
     1.26–1.36 — and its swing is correspondingly smaller: 0.054 against EME's
     0.074 over the wider span.
   - The **centre point does not move**: 0.3943 sits at 1.31 = `fcen`, the one
     index the pre-fix code labelled correctly. That is the same signature item
     9 used to tell the axis fix from a regression.

   `mmi_1x2_hybrid.json` falls with wavelength like EME, because the EME body
   dominates the trend — but it is contaminated all the same (0.4086 against
   EME's 0.4264 at 1.26). Direction is not a sufficient test for the hybrid;
   only the centre-point check is.

   **DONE 2026-07-26 — job 11500524, confirmed.** `mmi_1x2_eme.json` was
   deliberately not regenerated; it was never affected, and it is the one
   uncontaminated reference this device has.

   | `out_0` / `out_top` | 1.26 | 1.31 | 1.36 | direction |
   |---|---|---|---|---|
   | `mmi_1x2_fdtd`, stored | 0.37521 | **0.3942762** | 0.42954 | rising |
   | `mmi_1x2_fdtd`, new | 0.43283 | **0.3942750** | 0.33287 | falling |
   | `mmi_1x2_hybrid`, stored | 0.40863 | **0.3763402** | 0.34574 | falling |
   | `mmi_1x2_hybrid`, new | 0.43130 | **0.3763395** | 0.32767 | falling |
   | `mmi_1x2_eme` (untouched) | 0.42638 | 0.38632 | 0.35238 | falling |

   Both centre points are unmoved — **1.1e-6** for FDTD, **7.2e-7** for hybrid —
   at 1.31 = `fcen`, the one index the pre-fix code labelled correctly. That is
   item 9's signature, and it is what distinguishes the axis fix from a
   regression: anything that had perturbed the geometry, mode count or physics
   would have moved the centre too. The FDTD golden's direction now matches EME.
   `check_passive()` true for both. `out_1`/`out_bot` track their partner to 4
   decimals, so splitter symmetry survived.

   **Verified end to end:** `pytest tests/test_golden_examples.py --runslow -k
   mmi_1x2` — **3 passed** (job 11500663, 27:48). Two of the three were failing
   before.

   *One caveat, recorded so it is not mistaken for a discrepancy later.* The
   relabelling cross-check that validated the coupler to 0.00025 in item 5 is
   much weaker here. The pre-fix true axis was 1.2854 / 1.3100 / 1.3355;
   interpolating the new curve onto those points gives 0.416 / 0.394 / 0.366
   against the pre-fix 0.4295 / 0.3943 / 0.3752 — agreement to ~0.013. But the
   coupler's figure came from a 40-length sweep and this is a 3-point quadratic
   interpolation, which is a weak instrument. Treat it as consistent with the
   fix, not as confirmation of it; the centre-point check is what carries the
   argument.

2. **Shared resolution knob — and this is the worst-converging device audited.**
   `mmi_1x2_hybrid.py:41` passes one `RESOLUTION = 40.0` to both sub-solvers,
   the defect `directional_coupler_validation.md` item 11 fixed for the coupler.

   The audit's mesh sweep (same job) gives `n_eff[0]`:

   | res | 40 | 60 | 80 | 100 |
   |---|---|---|---|---|
   | n_eff[0] | 2.693784 | 2.713348 | 2.716740 | 2.725059 |
   | step | — | +0.0196 | +0.0034 | +0.0083 |

   **Non-monotone** — the 80→100 step is larger than 60→80 — and still moving
   0.0083 at res 100, three times the 2×2's 0.0027 over the same interval. The
   weak 160 nm confinement that drives gap 3 shows up here too: modes reaching
   further into the cladding need a wider, finer mesh to settle.

   The coupler's res 80 is **not** transferable. That pin was justified in
   item 11 on Δβ, the even/odd *difference*, which converges much faster than
   absolute n_eff — all three audited devices are still moving at res 100.
   Establish which quantity *this* device's result depends on before picking a
   number; the coupler's lesson was that sweeping an insensitive observable
   (T_cross) hid a 7% error in a sensitive one (Δβ).

   Do this **after** gap 3's `N_MODES` reduction — dropping the near-cutoff
   bands removes the modes that are hardest to converge, so the answer will
   change.

   **✅ MEASURED 2026-07-31 — job 11591224 (1:05:15).** Same recipe as the
   2×2's own gap 2 (job 11590785, `mmi_2x2_validation.md`): this device's
   self-imaging length (`mm_length=4.5`) is set by Δβ_01 between its two
   lowest guided modes, not either mode's own n_eff. `verify_mmi_1x2_gap2_
   dbeta.py` swept res 40/60/80/100/120:

   | res | Δβ_01 (rad/µm) | step | relative step |
   |---|---|---|---|
   | 40 | 0.28396 | — | — |
   | 60 | 0.28149 | −0.00247 | 0.870% |
   | 80 | 0.28095 | −0.00054 | 0.192% |
   | 100 | 0.27998 | −0.00097 | 0.345% |
   | 120 | 0.27987 | −0.00011 | 0.039% |

   > ✅ **Re-measured dispersive 2026-08-01, job 11608941 (1:03:21): every
   > value above is reproduced to all five printed digits.** 0.28396 / 0.28149
   > / 0.28095 / 0.27998 / 0.27987, and both n_eff columns match too. The
   > res-120 closure stands unchanged across the material break.
   >
   > **This was predicted before the job ran, and the reason is visible in its
   > own log header:** `[gap2] mmi_1x2 body, section=(-1.0, 1.0), wl=1.31,
   > N_MODES=2` — a **single scalar wavelength**. `resolved_material_config()`
   > only sets `resolved_dispersive` when `wl.min() < wl.max()`, and a medium
   > built at one wavelength is exact Sellmeier whether or not the dispersive
   > flag is on. Same mechanism as gap 3's byte-identical re-audit (11608940).
   >
   > **Gaps 2 and 3 are now both measured rather than one measured and one
   > inferred by analogy** — which was the point of spending an hour of compute
   > on a result expected to be unchanged. The rule can now be stated from
   > evidence: **a quantity evaluated at a scalar wavelength is invariant
   > across the dispersion break; only genuinely broadband quantities (gaps 4
   > and 5) move.**
   >
   > ⚠ **The ε log test does not apply to this job.** `run_mmi_1x2_gap2_dbeta`
   > goes through meow/MPB, which never prints meep's `dielectric constant
   > epsilon diagonal` line — grepping for it returns nothing, which is not
   > evidence of anything. Here the discriminator is the `wl=` field in the
   > header. Check the mechanism that actually applies to the code path rather
   > than reaching for the same grep every time.

   **Non-monotone, unlike the 2×2's clean decay.** The 80→100 relative step
   (0.345%) is larger than 60→80's (0.192%) — the same wobble already noted
   in this device's raw n_eff[0], now showing up in the physically-relevant
   quantity too, not just the one gap 2 originally flagged as possibly
   insensitive. Applying the 2×2's own rule verbatim — pin at the last
   resolution whose incoming step still exceeds 0.1% — does **not** land on
   res 100 here: that step is 0.345%, still 3.5× over threshold. res 120 is
   the first point under it (0.039%). **Gap 2 closes for this device at res
   120**, one tier above the 2×2's res 100 — consistent with this being the
   worst-converging of the three devices audited (weak SOI_160 confinement,
   `mmi_1x2_validation.md` gap 3's own diagnosis).

   **Decision, 2026-07-31: hold shipped `RESOLUTION` at 40, do not bump to
   120.** Both `mmi_1x2_hybrid.py:61` and `mmi_1x2_fdtd.py:38` keep
   `RESOLUTION = 40.0`. Moving to 120 is a ~3× mesh-density increase per
   spatial dimension plus a smaller Courant timestep — full 3-D FDTD cost does
   not scale linearly with that, nothing in this repo has run this device's
   full-device or hybrid script anywhere near res 120, and existing batch
   budgets (5h/job) were sized around res 40. Gap 2's res-120 finding is
   recorded as a **known, unclosed mesh error** instead — the same treatment
   check 4's 9.11% already carries a resolution caveat under. Downstream
   steps (goldens, gap 4, any re-run of gap 5) proceed at `RESOLUTION = 40`
   with this gap cited as the reason the number isn't mesh-converged.

3. **Mode identity — RUN, FAILED at N_MODES=12, and the fix is DONE (2026-07-30).**

   Measured on the uniform body cross-section (`x ∈ [−1, +1]`, 1.25 µm clear of
   each taper seam) at 1.31 µm, both engines meshed at res 40, job 11500536.
   Cladding index 1.4462.

   | band | meep | meow | Δ | guided |
   |---|---|---|---|---|
   | 0 | 2.693780 | 2.693784 | −0.00000 | yes |
   | 1 | 2.634501 | 2.634514 | −0.00001 | yes |
   | 2 | 2.532938 | 2.532970 | −0.00003 | yes |
   | 3 | 2.384403 | 2.384470 | −0.00007 | yes |
   | 4 | 2.181192 | 2.181325 | −0.00013 | yes |
   | 5 | 1.912386 | 1.912664 | −0.00028 | yes |
   | 6 | 1.773581 | 1.773560 | +0.00002 | yes |
   | 7 | 1.703956 | 1.704134 | −0.00018 | yes |
   | **8** | **1.5830239553117** | 1.589765 | **−0.00674** | yes |
   | **9** | **1.5830239553119** | 1.583575 | −0.00055 | yes |
   | 10 | 1.450410 | 1.450113 | +0.00030 | yes (barely) |
   | 11 | 1.432746 | 1.424674 | +0.00807 | **NO** |

   **meep returns bands 8 and 9 degenerate to 2.5e-13** — machine precision —
   while meow resolves them **6.19e-3 apart**. The nearest-neighbour
   correspondence is `[0,1,2,3,4,5,6,7,9,9,10]`, not the identity: both meep
   bands match meow's band 9, and **meow's band 8 has no meep counterpart at
   all** (nearest is 6.7e-3 away, an order of magnitude worse than any other
   guided band, all of which agree to ≤3e-4).

   Whichever engine is right, the consequence is the same: **index 8 of the FDTD
   interface basis is not the same mode as index 8 of the EME body basis**, and
   `HybridSolver` contracts them as though they were. Every 1×2 hybrid result
   currently mixes those two channels.

   **The obvious explanation is wrong, which is worth recording.** This gap
   originally predicted trouble because 12 modes crowd more than the coupler's
   6. But `mmi_2x2` asks for **15** and is the cleanest device audited —
   identity permutation, worst guided |Δn_eff| 0.0007. Mode *count* is not the
   driver. What distinguishes this device is **SOI_160**: a 160 nm core confines
   far more weakly than 220 nm, so the top of its basis piles up at cutoff.
   Band 10 sits only 0.0042 above the cladding index and band 11 is below it —
   the last four bands are all near or past cutoff, and that is exactly where
   two different grids stop agreeing about what the modes are.

   **Fix: drop `N_MODES` from 12 to 8** (`mmi_1x2_hybrid.py:95`). Bands 0–7
   agree to ≤1.8e-4 and are all comfortably guided; 8–11 are the entire problem
   region. This changes the device's configuration, so it invalidates the
   goldens just regenerated in gap 1 — they will need regenerating again, and
   that is the right order, not a reason to leave a broken basis in place.
   Re-run this audit afterwards to confirm the identity permutation holds.

   **DONE 2026-07-30.** `N_MODES` dropped to 8 in `mmi_1x2_hybrid.py:95`, and
   the audit (`examples/verification/python/audit_interface_modes.py mmi_1x2`,
   updated to match) re-run at the new count on the same cross-section (res
   40, 1.31 µm, `x ∈ [−1, +1]`):

   | band | meep | meow | Δ | guided |
   |---|---|---|---|---|
   | 0 | 2.696912 | 2.696915 | −0.00000 | yes |
   | 1 | 2.637697 | 2.637710 | −0.00000 | yes |
   | 2 | 2.536250 | 2.536282 | −0.00003 | yes |
   | 3 | 2.387898 | 2.387965 | −0.00007 | yes |
   | 4 | 2.184959 | 2.185091 | −0.00013 | yes |
   | 5 | 1.916513 | 1.916790 | −0.00028 | yes |
   | 6 | 1.775455 | 1.775436 | +0.00004 | yes |
   | 7 | 1.705834 | 1.706013 | −0.00018 | yes |

   All 8 requested modes are guided (no unbound column at this count, unlike
   the coupler and the 2×2), correspondence is the identity `[0,1,2,3,4,5,6,
   7]`, and worst guided |Δn_eff| is **0.0003** — the same numbers as bands
   0–7 measured at N_MODES=12 above, as expected since dropping unused high
   bands from the request doesn't change how the low ones are solved. The
   mismatched-basis problem this gap opened with is closed.

   **This invalidates the goldens confirmed in gap 1 a second time** — see
   `## Order` below; they must not be regenerated again until gap 2's
   resolution question is also settled, per the existing sequencing.

   *One thing this record does not say, noted 2026-07-30 so it is not assumed
   later:* whether the re-audit above ran with the duplicate-mode guards on.
   `FDTDSolver._iter_distinct_bands()` and `EMESolver.compute_modes()`'s
   `guard_duplicate_modes` both landed the **same day** as this re-audit, and
   the log evidence quoted here (a clean identity permutation, no skipped
   bands mentioned) is consistent with either. It matters less here than on the
   2×2 (`mmi_2x2_validation.md` gap 6) because 8 modes on a 2.8 µm body is well
   inside where the meep-side guard's diagnostic grid still discriminates
   (`resolved/fdtd_vs_EME_mode_validation.md` item 11), and because no guard warning
   appears in this device's audit output. Re-confirm opportunistically — with
   the next run of this audit, not as a job of its own.

4. ~~**Engine agreement has not been measured.**~~ **✅ MEASURED AND CLOSED
   2026-08-01 — jobs 11608875 `[0-3]` + merge 11609686.** The coupler's test
   uses even/odd supermode phase, which has no direct analogue in a 1×2 MMI.
   The transferable form is: solve a uniform stretch of the MMI body at two
   lengths, and check that FDTD and EME agree on each mode's accumulated phase.
   Adapted from `verify_item_6_engine_agreement.py` — the extraction helper
   (`tests/test_engine_agreement.py::dbeta_from_phases`) is device-independent
   and already tested.

   ### ✅ Check 1 result — the last of the four checks to land on this device

   Measured **dispersive** (`SECTION_LENGTHS = [0.5, 1.0, 1.5, 2.0]`,
   `N_MODES = 2`, `RESOLUTION = 40`, fully cache-cold). Costs 55:15 / 41:22 /
   49:33 / 2:23:09 plus a 25:45 merge. **All four checks are now ticked on the
   1×2**, matching the 2×2.

   | λ (µm) | Δβ FDTD | Δβ EME | ratio FDTD/EME |
   |---|---|---|---|
   | 1.26 | 0.26836 | 0.26937 | 0.99627 |
   | 1.30 | 0.27949 | 0.28103 | 0.99453 |
   | 1.31 | 0.28233 | 0.28396 | 0.99424 |
   | 1.35 | 0.29394 | 0.29579 | **0.99375** |
   | 1.36 | 0.29692 | 0.29877 | 0.99382 |

   **Mean ratio 0.9946 — a 0.54% bias, monotone across ten of the eleven
   points** (shallow minimum at 1.35, a +0.00007 tick back at the band edge).
   **No sign change.** Fit quality is excellent: worst FDTD phase-fit residual
   **3.01e-06 rad** and worst |off-diagonal| **1.12e-05**, so the two-mode
   basis is clean and the ratio is a real physical disagreement, not fit noise.

   > ⚠ **Both of those diagnostics print as `0.0000` in the job log** — the
   > `:.4f` trap that already burned the 2×2's merge. The values above come
   > from `engine_agreement.npz`. **Read the `.npz`, not the log.**

   **This lands almost exactly on the 2×2's number** (0.9962 → 0.9945, mean
   bias 0.47%), same direction — EME over-predicts Δβ relative to FDTD — and
   the same near-monotone shape. Two independent MMI geometries in different
   bands agreeing to within 0.07 points is a much stronger statement about the
   *engines* than either device made alone: this is an engine-level bias, not a
   device artifact.

   **Does it account for gap 5's 10.18%? No, and that is now settled from two
   directions.** The sign is consistent — a systematic Δβ under-prediction by
   FDTD is the pattern that produces gap 5's uniform-negative ΔT — but the
   magnitude is not: `calibrate_body_neff` corrects exactly this class of
   constant body-phase bias, and the P2 experiment above measured it buying
   only ~20% of the discrepancy (9.11% → 7.30%). So gap 4 is **consistent with
   but insufficient for** gap 5.

   Note this check does **not** depend on the `beta_from_phases` bug fixed the
   same day: it uses `dbeta_from_phases`, which is seed-free (Δβ ≈ 0.284
   rad/µm gives 0.142 rad per 0.5 µm step, nowhere near a wrap) and so has
   never had a branch to choose.

   ⚠ **Inherits the res-40 caveat** (see *Parity with the coupler* / P3): both
   engines are read at the same unconverged mesh, so a shared bias partially
   cancels in the ratio. The number is not certified mesh-independent.

5. **No hybrid-vs-full-device comparison on record.** Both scripts exist and both
   have goldens, but nothing asserts the hybrid reproduces `mmi_1x2_fdtd.py` on
   the same device. That is check 4 and it is the only one that tests the cascade
   as a whole. Cheapest real validation available here, since neither run has to
   be built first.

   **In progress 2026-07-30** — `examples/verification/python/
   verify_mmi_1x2_hybrid_vs_fdtd.py` added, ported directly from the 2×2's
   version (`mmi_2x2_validation.md` gap 5): loads both scripts' saved
   `SMatrix` pickles, reconciles `mmi_1x2_fdtd.py`'s `out_0`/`out_1` against
   `mmi_1x2_hybrid.py`'s `out_top`/`out_bot`, reports per-wavelength ΔT with
   provenance, no pass/fail threshold (same reasoning as the 2×2: gap 4 not
   measured yet). SLURM job 11585845 (array indices 1-2, re-queued
   2026-07-31 under `gts-sr108-gedc_main` — the original 07-30 submission
   never started, blocked on `AssocGrpBillingMinutes`) generates the two
   SMatrix files; job 11585850 (`--dependency=afterany:11585845`) runs the
   comparison once they land.

   **✅ MEASURED 2026-07-31 — job 11585850** (beef `1cc8edc`; FDTD data from
   11585845_1, hybrid from 11585845_2). This is **this device's first
   hybrid-vs-full-device number**, and check 4's first result on any MMI.

   > ⚠ **Frozen-material era.** Every number in this subsection was measured
   > before `dispersive` defaulted to True (2026-08-01). Gap 5 is broadband,
   > so it is one of the two gaps that genuinely moves across that break.
   > Preserved as `comparison_FROZEN_pre20260801.npz`; **the dispersive
   > re-measure has now landed — see the boxed result immediately below, which
   > is the current number.**

   > ### ✅ GAP 5 RE-MEASURED DISPERSIVE — 2026-08-01, jobs 11609215 `[1-2]` + 11609687
   >
   > **This is the number the shipped configuration produces.** Both sides
   > dispersive and mtime-verified in `provenance.json` (`fdtd_pkl` 17:34 UTC
   > from 11609215_1, `hybrid_pkl` 19:26 UTC from 11609215_2, beef
   > `f8b0d47+dirty` against the frozen baseline's `1cc8edc`) — no repeat of
   > the 11608692 pickle race. Comparison job 11609687, 0:12.
   >
   > | | frozen | **dispersive** |
   > |---|---|---|
   > | max \|ΔT\| (`out_0`) | 0.03588 | **0.04006** |
   > | mean \|ΔT\| | 0.02530 | **0.02541** |
   > | rel to mean T | 9.11% | **10.18%** |
   > | worst λ | 1.33 µm | **1.36 µm** (band edge) |
   > | all ΔT negative | yes | **yes** |
   >
   > **The gap got slightly worse, not better — 9.11% → 10.18%.** Say this
   > plainly: turning on the more physically correct material model moved this
   > device's check-4 discrepancy in the wrong direction. It is not a
   > regression in the material work (the O-band n_Si gap it fixed is real and
   > measured, `material_model_validation.md`), but it does mean the
   > material model was never what stood between this device and agreement.
   >
   > **What actually moved is only the band edge.** The mean offset is
   > essentially untouched — **−0.025304 → −0.025413, a 0.4% change** — while
   > the worst point migrated from 1.33 µm to the 1.36 µm edge, where the FDTD
   > reference itself shifted most (−0.0065 at 1.34, +0.0040 at 1.36). The
   > large uniform negative offset that *is* this device's gap 5 result is
   > **material-model-independent**, which is the strongest evidence yet that
   > it is structural.
   >
   > Splitter symmetry survives: `out_1` tracks `out_0` to five decimals
   > throughout (0.040060 vs 0.040076), as it did frozen.
   >
   > **Combined with the calibration result above, both cheap explanations for
   > this device's ~9–10% are now eliminated by measurement** — calibration
   > buys ~20% relative and leaves the sign structure intact; dispersion buys
   > nothing and costs 1 point. What remains is the tapers or the cascade, the
   > same conclusion the 2×2's P4 test reached independently on the same day.
   >
   > ⚠ **The calibrated number is still frozen-era.** The 7.30% in the P2 box
   > above has not been re-measured dispersive; a calibrated+dispersive run is
   > still owed before the two can be quoted side by side.

   ### ✅ P2 CLOSED, and calibration does NOT close this device's gap — 2026-08-01, jobs 11608690 + 11608692

   **`calibrate_body_neff=True` ran to completion on the 1×2 for the first
   time** — 11608690, 1:28:14, no residual guard fire, successor to the failed
   11589770. The `beta_from_phases` → `beta_fit_from_phases` fix works on the
   real device, not just on cached cells and synthetic seeds. That is P2's
   deliverable and it is met.

   > ⚠ **The job's own reported number, `max|ΔT| = 0.0313` / 7.95%, is
   > CONTAMINATED — do not quote it.** 11608692 loads the FDTD reference *by
   > path*, and job 11609215_1 had overwritten
   > `example_output/mmi_1x2/mmi_1x2_fdtd/smatrix.pkl` with the **dispersive**
   > reference nine minutes earlier (17:34:48 UTC vs the comparison's 17:43).
   > So it compared a **dispersive FDTD** against a **frozen calibrated
   > hybrid**. Verified from the meep logs, which are unambiguous: 11609215_1
   > prints ε = **12.4744** (dispersive), 11608690 prints ε = **12.2774**
   > (frozen). Caught only because `provenance.json` records `fdtd_pkl_mtime`
   > — a direct argument for extending P6-style provenance from goldens to
   > these comparison outputs, and for recording the *resolved material config*
   > rather than just an mtime.

   **The valid number was recovered offline, zero new compute**, by comparing
   11608690's saved calibrated hybrid against the preserved frozen reference
   `comparison_FROZEN_pre20260801.npz` — both sides frozen, like for like. The
   method self-checks: it reproduces the recorded uncalibrated 0.03588 / 9.11%
   exactly.

   | | max \|ΔT\| | mean \|ΔT\| | rel to mean T | all ΔT negative |
   |---|---|---|---|---|
   | uncalibrated | 0.03588 | 0.02530 | **9.11%** | yes |
   | **calibrated** | **0.02876** | **0.01916** | **7.30%** | yes |

   **Calibration buys ~20% relative and leaves the sign structure untouched.**
   Set against the coupler — **5% → 0.11%**, a 98% reduction — this is close to
   nothing. The conclusion is that **this device's 9.11% is mostly NOT the Δβ
   discretization gap `calibrate_body_neff` corrects.** The residual is
   structural, and gap 5's long-standing "trying it is the obvious next step"
   has now been tried and answered: it is not the explanation.

   Both remaining numbers are frozen-era; the dispersive calibrated re-run is
   still owed, and until it exists the calibrated result should be cited as
   7.30% (frozen, offline-derived) rather than the job's printed 7.95%.

   ---

   **Partial result 2026-08-01 — the FDTD side of that re-measure has landed
   (11609215_1, 13:44).** This is the full-device reference only; the hybrid
   half (`_2`) and the comparison (11609216) are still running, so **this is
   not yet a gap 5 number.** But the reference's own movement across the
   dispersion break is now measured directly against the preserved baseline:

   > **⚠ Superseded — the complete gap 5 result is the box above (10.18%).**
   > This subsection is kept because it isolates the *reference's* own
   > movement, which the finished comparison folds together with the hybrid's.
   > Its closing P4 argument was **written before the 2×2's test landed and is
   > now refuted** — see the correction at the end of it. (The comparison
   > actually ran as **11609687**, not 11609216.)

   | λ (µm) | T frozen | T dispersive | Δ |
   |---|---|---|---|
   | 1.260 | 0.44678 | 0.44716 | +0.00037 |
   | 1.270 | 0.43895 | 0.44068 | +0.00173 |
   | 1.280 | 0.42232 | 0.42384 | +0.00152 |
   | 1.290 | 0.41242 | 0.41259 | +0.00018 |
   | 1.300 | 0.40821 | 0.40821 | **+0.00001** |
   | 1.310 | 0.39585 | 0.39527 | −0.00058 |
   | 1.320 | 0.38730 | 0.38649 | −0.00080 |
   | 1.330 | 0.38053 | 0.37839 | −0.00214 |
   | 1.340 | 0.35853 | 0.35202 | **−0.00651** |
   | 1.350 | 0.33810 | 0.33743 | −0.00067 |
   | 1.360 | 0.34441 | 0.34845 | +0.00404 |

   Max |Δ| **0.00651**, i.e. **1.65%** of mean T. `out_1` tracks `out_0` to
   five decimals throughout, so **splitter symmetry survives the material
   change** as well as the cascade.

   **Read against P4, this is mixed evidence and should not be reported as
   support.** The mechanism (`mmi_2x2_validation.md` gap 5) predicts a frozen
   medium's error is zero at the freeze point and antisymmetric about it. There
   *is* a near-exact zero at **1.3001** (interpolated), one sample below band
   centre 1.310 — and the 2×2's crossing sits just below its centre too, the
   same direction. But the curve is **not** cleanly antisymmetric: it crosses a
   **second** time at 1.3514 and returns positive (+0.00404) at the band edge.
   A single frozen-index term cannot produce two crossings on its own, so
   either something else is superposed or the shape is not the simple thing the
   mechanism describes.

   **What it does settle is scale, and that is the useful part.** 0.0065 is
   **~20%** of this device's 0.021–0.036 gap 5 discrepancy — consistent with
   what was recorded when the 1×2 was used as P4's falsification control: here
   a large uniform-negative offset swamps any band-antisymmetric term. On the
   2×2, whose whole ΔT swing is ~0.003, a term of this class would **dominate**
   rather than hide. That is the quantitative version of the P4 prediction, and
   it was written down before 11609212 landed.

   > ### ❌ That scale argument was wrong, and the 2×2's own measurement is why
   >
   > It assumed a material term of **this device's** size would appear on the
   > 2×2. When the 2×2 was measured (11609680/11609681), its FDTD reference
   > moved by at most **0.00033** — not the 0.0065 seen here, but **~20×
   > less**, and only 8.6% of its ΔT swing. The term did not dominate there
   > because it is not that size there.
   >
   > **The unstated assumption was that material sensitivity transfers between
   > devices.** It does not: this device spans 1.26–1.36 µm (±3.8% about
   > centre) against the 2×2's ±1.13%, a 3.4× narrower relative band, and the
   > cost of freezing an index grows with distance from the freeze point. The
   > paragraph above extrapolated from the device with the *largest* material
   > sensitivity in the repo to the one with the smallest.
   >
   > The two-crossing shape flagged just above — the half that weakened the
   > hypothesis — was the better signal of the two and should have carried more
   > weight. Full refutation in `mmi_2x2_validation.md`'s P4 box.

   | | `out_0`/`out_top` | `out_1`/`out_bot` |
   |---|---|---|
   | max \|ΔT\| | **0.03588** | 0.03590 |
   | mean \|ΔT\| | 0.02530 | 0.02532 |
   | max \|ΔT\| / mean(T_fdtd) | **9.11%** | 9.11% |
   | worst λ | 1.33 µm | 1.33 µm |

   `check_passive()` true for both S-matrices. The two ports agree to four
   decimals, so **splitter symmetry survives the cascade** — a useful negative
   result, since a mode-mixing failure at an interface would not be expected to
   preserve it.

   **The disagreement is systematic, not scatter: every ΔT is negative** — the
   hybrid under-predicts transmission at all 11 wavelengths, by a margin that
   grows from 0.021 at 1.26 µm to 0.036 at 1.33 µm. A random basis or meshing
   error would not sign itself so consistently. That is the signature of the Δβ
   discretization gap `directional_coupler_validation.md` items 6/14
   diagnosed: MPB's index and the phase meep's time-stepper accumulates differ
   ~2.3% per mode, which the `β/Δβ` lever amplifies.

   **For scale:** the coupler's uncalibrated hybrid sat ~5% off its own
   full-device reference and dropped to **0.11%** with
   `HybridSolver(calibrate_body_neff=True)`. This device is at 9.11%
   *uncalibrated* and has never been run with that flag. Trying it is the
   obvious next step and costs one more hybrid run, not a new script.

   **❌ TRIED, FAILED — 2026-07-31, job 11589770 (4:17:26, then errored) →
   11589779 (dependent, failed in 9s with nothing to load).** Not a cluster
   fluke: `HybridSolver._calibrated_neff()` → `section_mode_beta()` raises
   outright —

   > `section_mode_beta: phase-fit residual 2.101 rad exceeds 0.1 at lengths
   > [0.923, 1.846, 2.770]. The measured phase is not linear in length, so
   > either the wrap counts are inconsistent between lengths or [-2.0, 2.0]
   > is not a uniform cross-section.`

   Calibration measures β by fitting FDTD's accumulated phase across three
   lengths of the body cross-section; here the fit itself is bad enough to
   trip its own 0.1 rad sanity guard, not just noisy. Two live hypotheses,
   neither checked yet: (1) the `[-2.0, 2.0]` section isn't actually uniform
   for this device the way the coupler's calibration cut was — plausible
   given `mm_length=4.5` is a lot less body to find 2 µm of clear uniform
   span in than the coupler's much longer straight; or (2) this is the same
   weak-SOI_160-confinement signature already driving gaps 2 and 3 — modes
   that piling up at cutoff also beat/leak enough over a few µm to break the
   linear-phase assumption `section_mode_beta` depends on. ~~**Not
   diagnosed.**~~ **Diagnosed 2026-07-31 late — see the box below; both of
   these two hypotheses are refuted.** ~~This closes off the obvious path to
   explaining 9.11% via the coupler's Δβ-discretization mechanism —
   calibration cannot currently be tried on this device at all, let alone
   shown to help.~~ **Revised: the path is not closed off, it is narrowed.**
   84 of 88 (wavelength, mode) cells calibrate cleanly; the blocker is a π
   sign flip in 4 cells (modes 0-1, λ ≤ 1.28), which is an upstream
   mode-extraction defect rather than a statement about whether this device
   can be calibrated.

   > ### ✅ FIXED 2026-08-01 — it was a **2π turn-count bug in
   > ### `beta_from_phases` itself**, not a π sign flip. Read this first;
   > ### everything below it, including the π-flip box, is superseded.
   >
   > **The π-sign-flip diagnosis recorded on 2026-07-31 (boxed below) was
   > wrong, and the way it was wrong is the useful part: the two readings are
   > algebraically indistinguishable.** For three equally spaced points,
   > "middle point displaced by −π" and "points 2 and 3 each displaced by
   > −2π" both give exactly `δ·(+⅓, −⅔, +⅓)`. The residual arithmetic that
   > looked like confirmation could not discriminate — it fits both.
   >
   > **Physics discriminates, decisively:**
   >
   > | reading | turns | β | n_eff | vs MPB seed 2.7295 |
   > |---|---|---|---|---|
   > | π sign flip in `T_ii` | `[0,2,4]` | 10.4752 | **2.1007** | **−23%** — wrong sign *and* size |
   > | 2π turn-count error | `[0,3,5]` | 13.8780 | **2.7830** | **+1.96%** — in family with the coupler's +2.233%/+2.361% |
   >
   > The FDTD time-stepper is *documented* to run a couple of percent above
   > MPB — that excess is the whole quantity `calibrate_body_neff` exists to
   > measure. 2.783 is it. 2.101 is below any plausible n_eff for a
   > 2.0 µm × 160 nm Si core at 1.26 µm.
   >
   > **Seed-free proof, which needs no MPB call and no branch assumption.**
   > For a linear phase ramp on equally spaced lengths the second difference
   > `(φ₃ − 2φ₂ + φ₁)/2π` must be an *exact integer*. Measured across the
   > whole grid: **within 0.00308 of an exact integer in all 88 cells**, 31 of
   > which need a ±1 turn correction. Every cell is a clean linear ramp;
   > **nothing anywhere needs a sign flip to become linear.** A π displacement
   > would sit at a *half* integer. This also re-refutes hypotheses (1) and
   > (2) below on stronger grounds than before — a non-uniform cut or an
   > unsettled near field would leave a *non*-integer second difference.
   >
   > Independent cross-check: that 0.00308 must equal the worst
   > `turn_ambiguity` the new fit reports, computed a completely different
   > way. It does, to 5 decimals.
   >
   > **The actual bug.** `beta_from_phases` rounded the wrap count
   > *independently at every length* against one seed slope, with no
   > cross-length consistency constraint. A seed error ε therefore produced a
   > fractional-turn offset `ε·L/(2π)` that **grew with L** until some length
   > straddled a `.5` boundary. The whole failure is visible in one row:
   >
   > ```
   > lengths                  [0.9232, 1.8465, 2.7697]
   > raw arg(T_ii)            [3.0618, -2.9846, -2.7285]
   > fractional turns vs seed [1.5127, 4.4750, 6.4343]
   >                                    ^^^^^^ 0.025 below the 4.5 boundary;
   >                                           the truth is 4.5539
   > ```
   >
   > ε here was **1.96% — the very FDTD-vs-MPB gap the calibration exists to
   > measure. The unwrapping was broken by exactly the quantity it was built
   > to measure.** That is why it hit modes 0–1 at the shortest wavelengths:
   > that is where β is *largest*, so branch spacing is tightest and the
   > fractional turn drifts furthest. The 2026-07-31 box read that corner as
   > evidence *against* a mesh/confinement story and therefore *for* a sign
   > flip; it is really the fingerprint of the branch arithmetic.
   >
   > **And it was never a corner case.** Seed-free sweep over every cell ×
   > every physically plausible branch × seed error ±3%: **4025 trials, the
   > old code wrong in 550 (13.7%), the new code wrong in 0.** On the
   > *coupler's* own cached data: 3650 trials, old wrong in 197 (5.4%), new
   > wrong in 0. The coupler's success was luck, not margin.
   >
   > **The doc's own "unresolved secondary finding" was the real tell.** That
   > cell's branch margin of **0.922** (against 0.5 warn / 1.0 useless) was
   > flagged below as possibly-caused-by-the-flip and left open. It was the
   > symptom. It is now **0.0785**.
   >
   > **Fix:** `beta_fit_from_phases` (`beef/solvers/fdtd.py`) decides exactly
   > **one** integer against the seed — the branch — and derives every later
   > turn from a slope refit on the shorter lengths, so ε·L drift cannot
   > accumulate. That *narrows* the seed's authority from n decisions to 1,
   > and 1 is precisely what `_BETA_MIN_BRANCH_FRACTION` already licensed.
   > `beta_from_phases` survives as a 3-tuple facade. New guard
   > `_BETA_BRANCH_MARGIN_RAISE = 0.9` makes `section_mode_beta`'s
   > long-standing docstring promise of a margin `ValueError` true for the
   > first time; the residual guard stays at 0.1 but now means exactly one
   > thing (the phase is not linear in L), since turn ambiguity no longer
   > contributes to it. Raise and log messages now name the offending
   > (wavelength, mode) cell — the old code kept only running maxima, which is
   > why "2.101 rad" could be reported without saying *where*.
   >
   > **Backward compatibility, measured not argued.** At the coupler's own
   > production seed the fix changes **0 of 66 cells** (`turns_shifted` 0
   > everywhere), so L₅₀ 0.11% is untouched. On this device it changes
   > **exactly the 4 of 88 cells** that blew the guard — the same 4 the box
   > below called sign-flipped. Existing suite: 28 of 29 tests unchanged; the
   > one that broke,
   > `test_a_boundary_seed_shows_up_as_a_broken_fit`, encoded the *old*
   > contract (it asserted a boundary seed shows up as a blown residual, which
   > is no longer true and was never a property of the physics). Rewritten as
   > `test_a_boundary_seed_is_undecidable_and_the_margin_is_what_says_so`.
   > Suite is now 42 tests, including
   > `TestTheMmi1x2TurnCountRegression` pinned to this cell's real cached
   > numbers.
   >
   > **In flight:** job **11608690** re-runs this device's calibrated hybrid
   > (`MMI1X2_CALIBRATE_BODY=1`, successor to the failed 11589770) with
   > **11608692** dependent on it for gap 5's calibrated comparison; job
   > **11608703** re-verifies the coupler's own calibrated L₅₀ through the
   > same shared code path.
   >
   > **Withdrawn, kept for the reasoning — the 2026-07-31 π-flip box:**
   >
   > `examples/verification/python/diagnose_mmi1x2_calibration_phase_fit.py`
   > (new) recovers the raw phase-vs-length data the failing fit was built
   > from — job 11589770's three `solve_section()` runs all persisted to the
   > on-disk cache, so this is a disk read, not a re-simulation
   > (`cache_only=True`, so a miss would have raised rather than silently
   > recomputed). It **reproduces the failure exactly: worst residual
   > 2.1008 rad against the job's reported 2.101.**
   >
   > **The failure is confined to 4 of 88 (wavelength, mode) cells:**
   >
   > | λ (µm) | offending modes | worst residual |
   > |---|---|---|
   > | 1.26 | 0 **and** 1 | 2.1008 |
   > | 1.27 | 0 | 2.101 |
   > | 1.28 | 0 | 2.100 |
   > | 1.29 – 1.36 | none | ≤ 0.0059 |
   >
   > Every other cell fits to **0.0014 – 0.0059 rad**, i.e. 17–70× *under*
   > the 0.1 rad guard. So calibration is not broadly broken on this device;
   > it is broken in one corner.
   >
   > **The mechanism, from the worst cell's shape** (λ=1.26, mode 0):
   >
   > ```
   > lengths        [0.9232, 1.8465, 2.7697]
   > unwrapped      [15.6281, 22.1481, 34.9706]
   > fitted line    [14.5777, 24.2490, 33.9202]
   > deviation      [+1.0504, -2.1008, +1.0504]
   > ```
   >
   > That deviation pattern is the least-squares fingerprint of **exactly one
   > displaced point**: for three equally spaced points with the middle one
   > offset by δ, the residuals are `δ·(−⅓, +⅔, −⅓)`, which for
   > δ = −3.1512 rad gives `(+1.0504, −2.1008, +1.0504)` — the observed row
   > to four decimals. Confirmed a second, independent way: the unwrapped
   > gaps are 6.52 and 12.82 (≈ 1:2), and adding 3.1512 back to the middle
   > point makes both gaps 9.6712 — exactly linear, same fitted slope.
   >
   > ~~**δ ≈ −π (3.1512 vs 3.14159, 0.3% off) means a sign flip of `T_ii`**,
   > since `arg(−z) = arg(z) ± π`. It is *not* a 2π branch error (that would
   > be a multiple of 6.2832 and, per `beta_from_phases`' own docstring,
   > would leave a *clean* residual) and *not* a smooth non-linearity (that
   > would deviate monotonically across all three lengths, not displace one).~~
   >
   > **↑ THIS IS THE WRONG STEP.** It *was* a 2π branch error. The
   > "multiple of 6.2832 ⇒ clean residual" reading of the docstring only holds
   > when *every* length is off by the same whole number of turns; here two
   > lengths were off and one was not, which is neither a uniform branch error
   > nor a non-linearity. The `δ·(+⅓, −⅔, +⅓)` pattern does not distinguish
   > "one point off by π" from "two points off by 2π" — it fits both exactly.
   >
   > **This actively refutes hypothesis (2) below**, rather than merely
   > failing to confirm it: weak SOI_160 confinement and near-cutoff
   > behaviour would strike the **highest** mode indices at the **longest**
   > wavelengths. The observed failure is the exact opposite corner — modes
   > **0 and 1** (the most strongly confined) at the **shortest** wavelengths.
   > Hypothesis (1) is refuted separately by the span arithmetic already
   > recorded below, and now also by the 84 clean cells on that same span.
   >
   > **Secondary finding, unresolved:** that cell's branch margin is
   > **0.922**, against `_BETA_BRANCH_MARGIN_MAX = 0.5` for a warning and
   > 1.0 for "the seed could not distinguish the branches at all". Fitted
   > β = 10.475 rad/µm against MPB's seed 13.611 — 23% apart. So even setting
   > the π flip aside, this corner is close to the branch-ambiguity limit.
   > Whether the sign flip *causes* the margin blow-up (a flipped `T_ii`
   > misleads the turn selection) or both are symptoms of one upstream mode
   > problem is **not** established. Note modes 0 *and* 1 both fail at
   > λ=1.26 while only mode 0 fails at 1.27/1.28, which is at least
   > consistent with a 0↔1 mode-ordering swap at that corner — a hypothesis,
   > not a result.
   >
   > **Suggested next step:** this is now a narrow, cheap question — dump
   > `T` itself (not just `arg(T_ii)`) for the three lengths at λ=1.26 and
   > check whether mode 0's column has swapped with mode 1's, or whether the
   > diagonal entry is simply negated. ~~Either way the fix is upstream of
   > `section_mode_beta`, which is behaving correctly here: **its guard
   > caught a real defect rather than a false alarm**, which is worth
   > recording as a positive result for the guard.~~
   >
   > **↑ Half right, and the wrong half mattered.** The `T` dump was the right
   > call and was done (the diagnostic script now retains `T`): the diagonal is
   > **not** negated and mode 0's column has **not** swapped. But the fix was
   > *not* upstream — it was inside `section_mode_beta`'s own helper. The guard
   > did catch a real defect, just not the one assumed: **it caught a defect in
   > its own unwrapping.** Recording it as "a positive result for the guard"
   > is what made the bug look like someone else's problem for a day.

   **Older, superseded analysis (2026-07-31, geometry-only) — kept because
   its span arithmetic is still correct and still rules out the strict
   "no room for 3 lengths" reading:** `[-2.0, 2.0]` is exactly
   `boundaries[1]`/`boundaries[2]` — i.e. `±half_body ∓ body_stub` =
   `±2.25 ∓ 0.25` — the **entire** EME-body span this device hands to
   calibration, the same pattern the coupler uses (`section_mode_beta`
   called on `boundaries[1]`/`boundaries[2]` there too), just on a body
   90%+ shorter. Checking `section_mode_beta`'s own length-selection logic
   against the failed run's reported lengths (`[0.923, 1.846, 2.770]`):
   `usable_span = span − 2·dpml = 4.0 − 1.0 = 3.0` µm, and the chosen
   spacing (0.923) is *below* the span-based cap (`usable_span/3 = 1.0`) —
   so the β-branch-safety rule, not the physical span, was what actually
   bound the fit here, and the longest fitted length (2.770) still fits
   inside the 3.0 µm usable span with ~0.23 µm to spare. **That rules out
   the fit hitting a hard span limit** — hypothesis (1) in its strictest
   "not enough room to even construct 3 valid lengths" form doesn't hold.

   It does **not** rule out a softer version of (1) combined with (2):
   `body_stub = 0.25` µm (this device's clearance past the taper seam) is
   *less than* `dpml = 0.5`, a pattern shared with every other device in
   this repo (the standard "monitor sits body_stub past the seam"
   convention), so that alone isn't unique to the 1×2. What **is** unique is
   that this device's *total* calibration span (4.0 µm) is the shortest of
   any device calibration has been tried on by a wide margin — the coupler's
   equivalent span is its whole coupling body, tens of µm. If SOI_160's weak
   confinement (independently established by gaps 2 and 3 on this exact
   device) means the near-field genuinely needs more than ~0.25-0.5 µm past
   the taper to settle into a pure uniform-body mode, a 4 µm-wide device
   simply has nowhere for that settling to finish before the fitted lengths
   run out of room — a coupler-scale body would absorb the same settling
   distance as a rounding error. **Best-supported reading: this is a
   combination of hypotheses (1) and (2), not a clean either/or** — the
   short body (a geometric fact about this device) and weak confinement (a
   physical fact about SOI_160, already twice implicated here) compound
   rather than either alone being sufficient. The residual's size (2.101 rad
   against a 0.1 rad guard, 21×) is too large to read as a marginal
   near-miss, consistent with a genuine near-field/non-uniformity effect
   rather than pure numerical noise.

   **Concrete next diagnostic step, not yet run:** plot the raw fitted phase
   vs. length points from job 11589770's own data (if still recoverable) or
   re-run just the 3-length phase measurement standalone — the *shape* of
   the non-linearity (a smooth curve vs. a clean jump between two lengths)
   would distinguish "genuine near-field settling" from "inconsistent 2π
   wrap count," which the residual number alone cannot.

   **Two caveats on reading 9.11% as the device's real error.** Both examples
   still share `RESOLUTION = 40.0`, which gap 2 shows is not converged here
   (and this is the worst-converging device audited), so part of this is mesh.
   And no threshold is enforced, deliberately: gap 4 (engine agreement) is
   still unmeasured for this device, so there is no basis to split this number
   into calibratable versus structural. Output:
   `example_output/verification/mmi_1x2_hybrid_vs_fdtd/`
   (`hybrid_vs_fdtd.png`, `comparison.npz`, `provenance.json`). Also run *before* gap 2 here, same as the
   2×2 — this device's basis is already fixed (gap 3, `N_MODES` 12→8) so the
   ordering concern the `## Order` section raises below (comparing a
   mismatched basis against a reference with "no such notion") no longer
   applies; only the resolution question (gap 2) is still open underneath
   this first pass.

## Order

~~**1** first~~ — done 2026-07-26, confirmed. The order below is **changed by
what gap 3 found**: the goldens are correct for the device as currently
configured, but the configuration itself is wrong.

0. **Wait on the coupler's cold regeneration pass before touching this device's
   goldens.** *(added 2026-07-28.)* The material-model change at the top of this
   file (`directional_coupler_validation.md` item 8) already invalidates
   `mmi_1x2_hybrid.json`/`mmi_1x2_fdtd.json` a third time, on top of gaps 3 and
   2 below — and that pass is not yet run. Folding it in here is exactly the
   same reasoning step 1 below already states for gaps 3+2 (don't regenerate
   against a configuration you're about to change again), just extended to a
   repo-wide change instead of a device-local one. Once that pass lands,
   `N_MODES` 12→8 and the resolution fix are config changes on top of an
   already-current material model, not underneath a stale one — one cold run
   covers all three instead of staggering three.
~~**1** `N_MODES` 12 → 8~~ — **done 2026-07-30**, confirmed: identity
permutation `[0..7]`, worst guided |Δn_eff| 0.0003, all 8 requested modes
guided. See gap 3.
~~**2** Gap 2, the resolution question, on the reduced basis~~ — **done
   2026-07-31** (job 11591224): converges at res 120, not the 2×2's res 100 —
   see gap 2. **Decision, 2026-07-31: shipped `RESOLUTION` stays 40** (cost of
   120 not worth it, existing budgets sized for 40); the gap is recorded as a
   known unconverged-mesh caveat rather than closed by config change. Steps
   3–5 below proceed at `RESOLUTION = 40`.
3. **Regenerate the goldens again.** Gaps 3 and 2 both change the configuration,
   so the values confirmed above will no longer apply. Recording a golden
   against a configuration you intend to change is how
   `directional_coupler_hybrid.json` ended up regenerated twice — the fix is to
   let the configuration settle first, not to skip the regeneration.
4. **Gap 5** (hybrid vs full-device FDTD). Both runs already exist, so this is
   the cheapest real validation available here — but it is only meaningful once
   the basis is right, since a cascade contracting mismatched modes would be
   compared against a reference that has no such notion.
5. ~~**Gap 4** last: the most work, and on the coupler the check that has
   stayed open longest.~~ **🟢 In flight 2026-08-01 — and it was much less
   work than this line assumed.** `verify_mmi_1x2_engine_agreement.py` already
   existed, fully ported to this device (`SECTION_LENGTHS = [0.5, 1.0, 1.5,
   2.0]`, ≥1.25 µm clear of the ±2.25 taper seams, `N_MODES = 2`); it had
   simply never been submitted. Reshaped to the 2×2's proven array+merge shape
   (`SLURM_ARRAY_TASK_ID` picks one length; the merge step re-reads all four
   with `cache_only=True`) and submitted as **11608448 `[0-3]`** +
   **11608450**. Costs came in at 18:30 / 25:33 / 27:01 for the first three
   lengths — roughly 2 h of serial work against the 2×2's 13.4 h, as the
   6.2×-shorter body predicted.

   **✅ CLOSED 2026-08-01.** The 11608448 pass was frozen-era and its merge
   (11608450) was deliberately cancelled rather than run against a cache it
   could no longer reach. Re-run cache-cold and dispersive as **11608875
   `[0-3]`** + merge **11609686**: mean Δβ ratio **0.9946**, 0.54% monotone
   bias, worst fit residual 3.01e-06 rad. Full table in gap 4. Dispersive
   costs ran 55:15 / 41:22 / 49:33 / 2:23:09 — roughly 2.4× the frozen-era
   pass, consistent with the 1.53× seen on 11609215_1.

   Note this check does **not** depend on the `beta_from_phases` bug fixed the
   same day: it uses `dbeta_from_phases`, which is seed-free (Δβ ≈ 0.284
   rad/µm gives 0.142 rad per 0.5 µm step, nowhere near a wrap) and so has
   never had a branch to choose.

## Parity with the coupler — one item beyond the four checks

Added 2026-08-01, tracked as **P3** in [README.md](../README.md)'s parity
backlog and recorded here because gap 2 owns it.

**The resolution knob is still shared.** `mmi_1x2_hybrid.py` ships one
`RESOLUTION = 40.0` for both sub-solvers, while the coupler splits its own
(`FDTD_RESOLUTION = 40` / `EME_RESOLUTION = 80`) with the comment *"shared
RESOLUTION was how a 7% error hid for months (dc_validation.md item 11)"* —
i.e. the exact pattern this device still carries. Gap 2 measured the need here
as **res 120**, one tier above the 2×2's 100 and the worst of the three
audited devices, and the decision recorded above was to hold 40 on cost
grounds and carry the gap as a caveat.

That decision is defensible and is not being reopened here. What was *not*
recorded is its consequence: **check 3 is not closed on this device the way it
is on the coupler**, because the coupler's closure was the knob split, not a
convergence measurement. Anything measured at res 40 on this device — gap 5's
9.11% (now 10.18% dispersive), gap 4's **measured 0.9946 ratio** — inherits
that caveat. Worth re-deciding
explicitly rather than leaving implicit, especially now that gap 4 gives a
second res-40 number to interpret.
