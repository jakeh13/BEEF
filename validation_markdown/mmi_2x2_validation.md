# 2×2 MMI validation — to-do

Opened 2026-07-26. Nothing is known to be wrong; this is the gap list against the
four checks in [README.md](README.md), none of which have been run here.

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

2. **Shared resolution knob.** `mmi_2x2_hybrid.py:69` passes one
   `RESOLUTION = 40.0` to both sub-solvers — the defect `dc_validation.md` item
   11 fixed for the coupler, where res 40 overstated Δn_eff by 7%. Split into
   `FDTD_RESOLUTION` / `EME_RESOLUTION` and measure convergence on *this*
   cross-section. The coupler's res 80 does not transfer: an MMI body is a wide
   multimode slab, and its higher-order modes converge differently from a
   two-waveguide coupled pair.

3. **Mode identity — this is the most exposed device in the repo.**
   `N_MODES = 15` at the interfaces (`mmi_2x2_hybrid.py:124`), against the
   coupler's 6 and the 1×2's 12. Higher-order slab modes crowd together, so the
   near-degeneracies that make FDTD and EME disagree about mode *order* are most
   likely here of anywhere. `HybridSolver` contracts index *i* against index *i*
   assuming they match; if they do not, the chain mixes modes silently.
   `FDTDSolver.section_mode_neff()` gives the FDTD side without time-stepping,
   so checking all 15 costs seconds. **Run this before trusting any 2×2 result.**

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

**3** first — it is cheap, and it is the one that could invalidate everything
else this device produces. Then **2**, then **5** (both runs already exist), then
**1** once the configuration has stopped moving. **4** last, as on the coupler.
