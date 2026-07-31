# 1×2 MMI validation — to-do

> **2026-07-28 — the material model changed under every device here.** Both engines
> now evaluate one analytic Sellmeier model (Salzberg Si, Malitson SiO2) instead of
> each doing its own nearest-tabulated lookup, which shifts n_Si by ~0.002 and
> removes a 0.0063 step that sat in the middle of the C-band. Any golden or stored
> result for this device predates that and needs regenerating; see
> `directional_coupler_validation.md` item 8 and, for the still-open dispersive-FDTD
> half, `material_model_validation.md`.



Opened 2026-07-26. Gap list against the four checks in [README.md](README.md).

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
   (`fdtd_vs_EME_mode_validation.md` item 11), and because no guard warning
   appears in this device's audit output. Re-confirm opportunistically — with
   the next run of this audit, not as a job of its own.

4. **Engine agreement has not been measured.** The coupler's test uses even/odd
   supermode phase, which has no direct analogue in a 1×2 MMI. The transferable
   form is: solve a uniform stretch of the MMI body at two lengths, and check
   that FDTD and EME agree on each mode's accumulated phase. Adapt
   `verify_item_6_engine_agreement.py` — the extraction helper
   (`tests/test_engine_agreement.py::dbeta_from_phases`) is device-independent
   and already tested.

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
   comparison once they land. Also run *before* gap 2 here, same as the
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
2. **Gap 2**, the resolution question, on the reduced basis. This is now the
   next open step — nothing else in this file is blocked on anything but this
   and step 0 above.
3. **Regenerate the goldens again.** Gaps 3 and 2 both change the configuration,
   so the values confirmed above will no longer apply. Recording a golden
   against a configuration you intend to change is how
   `directional_coupler_hybrid.json` ended up regenerated twice — the fix is to
   let the configuration settle first, not to skip the regeneration.
4. **Gap 5** (hybrid vs full-device FDTD). Both runs already exist, so this is
   the cheapest real validation available here — but it is only meaningful once
   the basis is right, since a cascade contracting mismatched modes would be
   compared against a reference that has no such notion.
5. **Gap 4** last: the most work, and on the coupler the check that has stayed
   open longest.
