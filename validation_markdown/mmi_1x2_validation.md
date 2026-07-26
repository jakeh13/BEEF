# 1×2 MMI validation — to-do

Opened 2026-07-26. Nothing is known to be wrong with this device; this is the
gap list against the four checks in [README.md](README.md), none of which have
been run here.

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

   `dc_validation.md` items 1–2 landed 2026-07-25, so both FDTD-path goldens
   were recorded with the reversed, half-span wavelength axis — exactly the
   state `directional_coupler_hybrid.json` was in before item 9 regenerated it.
   The `--runslow` golden tests for those two should be failing right now.

   Regenerate via `scripts/generate_golden_values.py mmi_1x2_fdtd` /
   `mmi_1x2_hybrid`, and verify the move is the axis fix rather than a
   regression the same way item 9 did: **the centre wavelength must not move**,
   since `fcen` is the one index the pre-fix code labelled correctly. Do not
   regenerate `mmi_1x2_eme.json` — it was never affected, and re-recording it
   would destroy the one uncontaminated reference this device has.

2. **Shared resolution knob.** `mmi_1x2_hybrid.py:41` passes one
   `RESOLUTION = 40.0` to both sub-solvers, the defect `dc_validation.md` item 11
   fixed for the coupler. The EME half is presumptively carrying the same 7%
   Δn_eff error. Split into `FDTD_RESOLUTION` / `EME_RESOLUTION` and measure the
   convergence on *this* cross-section — the coupler's res 80 is not
   transferable, since an MMI body is a much wider multimode slab and its higher
   modes will converge differently.

3. **Mode identity is unchecked, and this device is more exposed than the
   coupler.** `N_MODES = 12` at the interfaces (`mmi_1x2_hybrid.py:95`) against
   the coupler's 6. A wide MMI slab's higher-order modes are closely spaced, so
   the near-degeneracies that make FDTD and EME disagree about mode *order* are
   far more likely here. Use `FDTDSolver.section_mode_neff()` — no time-stepping,
   so it costs seconds — and compare all 12 pairwise.

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

## Order

**1** (are the goldens stale?) is the only one that could mean current recorded
values are wrong, so it comes first and is cheap — a timestamp check, then one
batch job if it fails. Then **5**, which uses what already exists. Then **2** and
**3** together, since both are cross-section measurements on the same geometry.
**4** last: it is the most work and, on the coupler, the check that has stayed
open longest.
