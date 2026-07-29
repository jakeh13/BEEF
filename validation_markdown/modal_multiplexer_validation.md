# Modal multiplexer validation — to-do

> **2026-07-28 — the material model changed under every device here.** Both engines
> now evaluate one analytic Sellmeier model (Salzberg Si, Malitson SiO2) instead of
> each doing its own nearest-tabulated lookup, which shifts n_Si by ~0.002 and
> removes a 0.0063 step that sat in the middle of the C-band. Any golden or stored
> result for this device predates that and needs regenerating; see
> `directional_coupler_validation.md` item 8 and, for the still-open dispersive-FDTD
> half, `material_model_validation.md`.



Opened 2026-07-26. Nothing is known to be wrong; this is the gap list against the
four checks in [README.md](README.md), none of which have been run here.

`examples/modal_multiplexer/python/modal_multiplexer_hybrid.py` — the only script
for this device. SOI_160, not SOI_220 like the coupler and MMIs.

## Where it sits

The **most demanding** hybrid device here and the **least supported**:

- `N_MODES = 16` carried at **four** section interfaces (`:153`, `:208`) — more
  modes and more joins than any other example. Every interface is another place
  the FDTD and EME bases have to agree.
- **No full-device FDTD script.** Unlike both MMIs and the coupler, there is
  nothing to compare the cascade against as a whole; check 4 has no starting
  point at all.
- **No golden fixture.** Nothing pins its output.
- Different material stack (SOI_160), so no measurement from the other devices
  transfers numerically — not the coupler's converged EME res 80, not its Δβ
  ratio.

So it has the most ways to be wrong and the fewest ways to find out.

## Gaps

1. **Mode identity, at four interfaces with 16 modes each.** `HybridSolver`
   contracts index *i* of an FDTD interface against index *i* of the EME body at
   every join, assuming they are the same mode. With 16 modes the higher-order
   ones are closely spaced and the two engines solve MPB on different grids, so
   reordering is more likely here than anywhere else — and a swap at *any* of the
   four interfaces corrupts the chain silently. `FDTDSolver.section_mode_neff()`
   costs seconds per cross-section; check every interface, not just one.

   This is the single highest-value check for this device. A modal multiplexer's
   entire purpose is keeping modes distinct, so a solver that mixes them is
   failing at exactly the thing being demonstrated.

2. **Shared resolution knob.** `:86` passes one `RESOLUTION = 40.0` to both
   sub-solvers (`:167` FDTD, `:186` EME) — `directional_coupler_validation.md` item 11's defect.
   Must be measured fresh on SOI_160: a 160 nm core confines more weakly than
   220 nm, so its modes extend further into the cladding and may well need a
   *different* mesh than the coupler's res 80 to converge. Do not copy the number
   across.

3. **No full-device FDTD reference.** Needs a `modal_multiplexer_fdtd.py` written
   before check 4 is possible. Cost is the concern — a multiplexer is long, and a
   single-domain 3-D run at res 40 over that span is expensive (the coupler's
   78 µm device took ~3 h on 48 ranks). Worth scoping at reduced resolution
   first: even a res-20 reference that only resolves the qualitative
   mode-sorting would catch a cascade that is mixing modes, which is the failure
   mode gap 1 is about.

4. **No golden.** Add a `modal_multiplexer_hybrid` target to
   `tests/test_golden_examples.py`'s `TARGETS`. As with the 2×2 MMI, do it after
   the resolution question settles, so the recorded values are not immediately
   invalidated.

5. **Engine agreement not measured.** Adapt
   `verify_item_6_engine_agreement.py` onto a uniform stretch of the body. The
   extraction helper (`tests/test_engine_agreement.py::dbeta_from_phases`) is
   device-independent.

## Order

**1** first and by a wide margin — cheap, and it tests the property this device
exists to demonstrate. Then **2**. Then **4**, to stop the numbers drifting
unnoticed while the rest is worked on. **3** and **5** last: both are real
simulation campaigns, and neither is worth starting until 1 has ruled out a
mismatched basis, which would invalidate their results anyway.
