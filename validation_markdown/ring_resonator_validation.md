# Ring resonator validation — to-do

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

1. **Re-run both scripts** and check the centre wavelength does not move. This
   is what turns a code fix into a trustworthy result, and nothing below is
   worth reading off a plot until it happens.
2. **Gap 3 (resolution convergence)** — cheap, and it gates every resonance
   number. Do it before or alongside the re-run rather than after, so the
   re-run is not immediately invalidated by a resolution change. `directional_coupler_validation.md`
   item 11 found the shared res=40 overstates Δn_eff by 7% on a coupled pair,
   and the audit's coupler sweep shows absolute `n_eff` still moving 0.0044
   between res 80 and 100 — the ring's bare-strip n_eff feeds three propagation
   phases *and* the reference-plane correction, so it needs its own number.
3. Then the (a)/(b)/(c) work below, after which gaps 1, 2 and 4 become possible
   in the same form as the other devices.
