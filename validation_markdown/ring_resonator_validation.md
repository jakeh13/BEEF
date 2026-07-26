# Ring resonator validation — to-do

Opened 2026-07-26. **The least validated device in the repo**, and the only one
with a known-live bug.

`examples/ring_resonator/python/` — `ring_resonator_hybrid.py`,
`ring_resonator_hybrid_cascaded.py`. Despite the `_hybrid` names, neither uses
`HybridSolver`, `FDTDSolver` or `EMESolver`: both hand-build their meep
simulations and their meow cross-section directly. That single fact is the root
of everything below.

## Known defect — the axis bug is still live here

`ring_resonator_hybrid.py` (lines 385, 469) and
`ring_resonator_hybrid_cascaded.py` (lines 392, 463) each compute

```python
fwidth = (1/WL_MIN - 1/WL_MAX) * 0.5
```

and hand it to `sim.add_mode_monitor(fcen, fwidth, N_WL_FDTD)` — character for
character the defect `dc_validation.md` item 1 fixed, in raw meep rather than
through `FDTDSolver`, so the solver fix does not reach it. **Every ring FDTD
result carries the same reversed, half-span wavelength axis.**

Worse, their 6 cache entries are still *reachable*: these scripts build their own
flat config dicts via `_sim_params()` (`ring_resonator_hybrid.py:297`), so
nothing in `dc_validation.md` items 1–3 changed their hash and a re-run silently
reuses pre-fix data. Entries were left on disk deliberately — deleting them only
forces re-simulation of equally wrong results.

This is `dc_validation.md` item 10. Fix: port both scripts to
`monitor_freq_grid()` / `source_fwidth()`, add the resolved grid to
`_sim_params()` so the hash strands the old entries automatically (the mechanism
item 2 used), then re-run.

For a resonator the consequence is worse than for a coupler. A reversed,
half-span axis does not just mislabel a smooth curve — it puts every **resonance
wavelength** in the wrong place and compresses the free spectral range by 2×.
Any FSR or Q read off these plots is wrong.

## Recommendation — put the ring on `HybridSolver`

Item 10 exists *because* this device is off the library path. Items 1–4, 11 and
12 all landed in `beef/` and reached every other example automatically; the ring
needed each ported by hand, and the axis fix was not. That will keep happening.
The durable fix is to make the ring a first-class device, which needs three
things — the first two mechanical, the third a real design change.

**(a) Arbitrary-orientation ports in `FDTDSolver`.** Every source and monitor in
`beef/solvers/fdtd.py` is x-normal — `mp.Vector3(0, w, h)`, 8 occurrences. The
ring's two coupling ports sit on top of the near bends and are **y-normal**
(x–z planes). Nothing in BEEF draws or measures those, and
`beef/geometry/_compose.py` has translation only, no rotation. Needs a port
orientation carried through `section_fdtd_layout()` and the monitor builders.

**(b) A `RingResonator` DeviceGeometry**, so `section_gds()` / `section_prisms()`
work and the device stops hand-building meep objects.

**(c) `HybridSolver` must express feedback.** This is the substantive one. The
solver today contracts a straight cascade — `np.einsum("mij,mjk->mik", ...)` at
`hybrid.py:318` and `:416` — threading light through the sections once. A ring
is a **closed loop**: light re-enters the coupling region every round trip, which
is a closed-form all-pass sum, not a matrix product. Two ways out:

- *Narrow:* keep the cascade for the open path (coupler section + bends +
  straights) and add an explicit loop-closure step that applies the resonator sum
  to a designated feedback port pair. Less work, and it matches how the script
  already assembles the answer — but it hard-codes "ring" into the solver.
- *General:* let `HybridSolver` describe an S-matrix **network** with named ports
  rather than an ordered chain, and solve it by Redheffer star products /
  node elimination. Cycles then fall out for free, and it would also cover
  Mach–Zehnders, cascaded rings, and any device where a section feeds back. This
  is the architecturally right answer and considerably more work.

Recommend the general form if `ring_resonator_hybrid_cascaded.py` is a sign of
where this is heading — cascaded rings are exactly the case the narrow option
starts fighting. Either way, (a) and (b) are prerequisites and can land first,
independently, and are useful on their own.

Until then the ring cannot be validated the way the other devices are: checks 1–3
in [README.md](README.md) all assume `FDTDSolver`/`EMESolver` objects that this
device never constructs.

## Gaps against the four checks

1. **Engine agreement** — not measurable today. The mode solve is a direct meow
   call, and the FDTD is raw meep, so there is no shared cross-section object to
   compare. Blocked on the recommendation above.
2. **Mode identity** — not applicable in the current form (no interface mode
   basis is exchanged; the scripts thread scalar coefficients). Becomes check 2
   proper once (c) lands.
3. **Resolution convergence** — never measured, for either engine. The mode solve
   feeds `neff(λ)` into three straight-propagation phases *and* the reference-plane
   correction on `b(λ)`, so an n_eff error propagates into resonance position
   directly. Given `dc_validation.md` item 11 found the shared res=40 overstates
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

**Fix item 10 first**, on its own, without waiting for the recommendation: it is
a known-wrong result on disk, the fix is mechanical, and it is worth doing even
if the scripts are later replaced. Then **3**, which is cheap and gates every
resonance number. Then the (a)/(b)/(c) work, after which 1, 2 and 4 become
possible in the same form as the other devices.
