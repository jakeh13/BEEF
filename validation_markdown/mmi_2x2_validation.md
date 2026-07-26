# 2×2 MMI validation — to-do

Opened 2026-07-26. Nothing is known to be wrong; this is the gap list against the
four checks in [README.md](README.md).

**Status 2026-07-26:** check 2 (mode identity) is **done and passed** — job
11500536, see below. It was expected to be the riskiest device in the repo and
turned out to be the cleanest. Checks 1, 3 and 4 remain unrun, and there are
still no goldens.

`examples/mmi_2x2/python/` — `mmi_2x2_eme.py`, `mmi_2x2_fdtd.py`,
`mmi_2x2_hybrid.py`.

## Where it sits

Structurally the same shape as the 1×2 MMI — three scripts, EME/FDTD/hybrid —
but with **no golden fixtures at all**. `tests/golden/` holds only
`directional_coupler_hybrid.json` and the three `mmi_1x2_*.json`. So unlike the
1×2, nothing here is pinned: any change to the solvers moves this device's
numbers silently.

That is the one thing this device has that the 1×2 does not have a version of:
the 1×2's problem is *stale* references, which at least fail loudly under
`--runslow`. Having none is quieter and worse.

## Gaps

1. **No goldens.** Add `mmi_2x2_eme`, `mmi_2x2_fdtd`, `mmi_2x2_hybrid` targets to
   `tests/test_golden_examples.py`'s `TARGETS` and generate them with
   `scripts/generate_golden_values.py`. Do this *after* gap 2, so the recorded
   values are not immediately invalidated by a resolution change — recording a
   golden against a configuration you already intend to change is how
   `directional_coupler_hybrid.json` ended up regenerated twice.

2. **Shared resolution knob — and the EME side is not converged at any
   resolution tried.** `mmi_2x2_hybrid.py:69` passes one `RESOLUTION = 40.0` to
   both sub-solvers, the defect `dc_validation.md` item 11 fixed for the
   coupler.

   The audit's mesh sweep (same job) gives `n_eff[0]` on this cross-section:

   | res | 40 | 60 | 80 | 100 |
   |---|---|---|---|---|
   | n_eff[0] | 2.806984 | 2.819049 | 2.825888 | 2.828554 |
   | step | — | +0.0121 | +0.0068 | +0.0027 |

   Still moving 0.0027 at 80→100, and the steps shrink by only ~2× each time,
   so this is a slow tail rather than a converged value. **The coupler's res 80
   does not transfer** — that pin was justified in item 11 on Δβ, the even/odd
   *difference*, which converges far faster than absolute n_eff. All three
   audited devices are still moving in absolute n_eff at res 100.

   Whether that matters depends on which quantity this device's result actually
   depends on. Establish that first — the coupler's lesson was that sweeping a
   quantity the physics is not sensitive to (T_cross) hid a 7% error in one it
   is (Δβ). Do not simply raise the resolution and assume it is fixed.

3. ~~**Mode identity — this is the most exposed device in the repo.**~~
   **Done 2026-07-26 — passed, and the prediction behind this gap was wrong.**

   Measured on the uniform body cross-section (`x ∈ [−1, +1]`, 13 µm clear of
   either taper seam) at 1.55 µm, both engines meshed at res 40, job 11500536:

   | | |
   |---|---|
   | modes requested | 15 |
   | **guided** (n_eff > 1.4443) | **14** |
   | correspondence | identity, `[0 … 13]` |
   | worst guided \|Δn_eff\| | **0.0007** |

   The reasoning in this gap — that 15 modes crowd together and so this device
   is the most likely to reorder — **did not hold**. It is the cleanest of the
   three devices audited, better than the coupler's 0.0023 on 6 modes and an
   order of magnitude better than the 1×2's 0.0067 on 12. Mode *count* is not
   what drives disagreement; confinement is. This device is SOI_220 on a 2.8 µm
   body, and even band 13 sits at n_eff 1.508, comfortably above the 1.4443
   cladding. The 1×2 is SOI_160, and a 160 nm core holds its top bands so
   weakly that they pile up at cutoff — that is where its failure is
   (`mmi_1x2_validation.md` gap 3).

   **Band 14 (n_eff 1.4078) is below the cladding index** — a radiation mode.
   So one column of every 2×2 interface basis is unguided, the same finding as
   the coupler's band 5, and its +0.0021 disagreement is 3× worse than any
   guided band here. `N_MODES = 14` would drop it at no cost.

4. **Engine agreement not measured.** Adapt
   `verify_item_6_engine_agreement.py`: solve a uniform stretch of the MMI body
   at two lengths and check the two engines agree on each mode's accumulated
   phase. The extraction helper
   (`tests/test_engine_agreement.py::dbeta_from_phases`) is device-independent
   and already tested — it takes lengths and phases, nothing coupler-specific.

5. **No hybrid-vs-full-device comparison on record.** `mmi_2x2_fdtd.py` exists,
   so the reference run is available; nothing asserts the hybrid reproduces it.
   This is the only check that tests the cascade as a whole.

## Order

~~**3** first~~ — done 2026-07-26, passed. Remaining:

1. **Gap 2**, the resolution question — it gates gap 1, and the audit has
   already shown res 40 is not converged.
2. **Gap 5** (hybrid vs full-device FDTD), which needs no new scripts.
3. **Gap 1** (goldens) once the configuration has stopped moving. Dropping
   `N_MODES` to 14 per gap 3 counts as the configuration moving, so do it
   before recording anything.
4. **Gap 4** last, as on the coupler.
