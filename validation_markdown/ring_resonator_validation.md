# Ring resonator validation — to-do

> **2026-07-28 — the material model changed under every device here.** Both engines
> now evaluate one analytic Sellmeier model (Salzberg Si, Malitson SiO2) instead of
> each doing its own nearest-tabulated lookup, which shifts n_Si by ~0.002 and
> removes a 0.0063 step that sat in the middle of the C-band. Any golden or stored
> result for this device predates that and needs regenerating; see
> `directional_coupler_validation.md` item 8 and, for the still-open dispersive-FDTD
> half, `material_model_validation.md`.



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

**(a) Arbitrary-orientation ports in `FDTDSolver` — done, 2026-07-29.** Every
source and monitor in `beef/solvers/fdtd.py` was x-normal — `mp.Vector3(0, w,
h)`, 8 occurrences — and the ring's two coupling ports, sitting on top of the
near bends, are **y-normal** (x–z planes). `solve_section()` now takes
`y0`/`y1` (a y-axis domain bound, mirroring `x0`/`x1`) and `excite_port`/
`measure_port` (`PortSpec(axis="x"|"y", sign=±1, center, width)`), additive
and default-`None` — every existing caller is unaffected, verified by the
full fast suite staying at 338 passed after the change.
`section_fdtd_layout()`/`SectionFDTDLayout` are untouched; `PortSpec`
bypasses them per-side rather than extending their scalar contract.

Two real bugs found while building this, not just anticipated: (1)
`section_fdtd_layout()`'s `reverse = n_out <= n_in` tie-break silently
excites from the wrong physical side whenever `n_in == n_out` — exactly the
standalone-bend case (`n_modes_in=n_modes_out=1`) — now bypassed via a forced
`reverse=False` whenever either `PortSpec` is given, with the excited/
measured role taken literally from `n_modes_in`/`n_modes_out`. (2) a
measure port's edge convention is the *opposite* of an excite port's for the
same `sign` value — a source's launch direction physically forces it to a
specific edge, but a passive monitor has no such constraint, so `sign=+1`
means "low edge, launch +axis" for `excite_port` but "high edge, read
+axis-forward light" for `measure_port`. Caught mid-implementation, before
either shipped silently wrong.

**Verified against real physics, not just unit tests.** A
`@pytest.mark.slow` check (`tests/test_ring_orientation.py`, run with
`--runslow`) builds the new `RingResonatorBendSection` (below), excites
x-normal, measures y-normal through the 90° bend at reduced resolution
(res 15, 2-D — a login-node smoke run, not production), and passes: meep's
own log shows the y-normal monitor resolving a mode whose dominant
planewave is genuinely y-directed (`(0, 2.01, 0)`) against the x-normal
source's x-directed one (`(2.04, 0, 0)`) — physical confirmation, not a
numeric coincidence. **Not yet run at production resolution/3-D**, and the
y-normal *excite* + flux-calibration path (needed for the coupler section's
actual ring input, not exercised by the bend alone) is still unverified by
any real meep run.

**(b) A `RingResonator` DeviceGeometry — done, 2026-07-29.**
`RingResonatorCouplerSection` (bus + one near-bend + ring stub, `side="left"`
or `"right"`) and `RingResonatorBendSection` (the standalone 90° bend,
x-normal in / y-normal out) both implement `section_prisms`/`section_gds`/
`width_at`/`center_at`, in `beef/geometry/ring_coupler.py` and `ring_bend.py`.
Breakpoints verified to match the reference scripts' own constants exactly
(`stub_x=-7.0`, `hc=-2.0`, `iface_x=-1.0`, `ring_port().y_boundary=6.20` for
`side="left"`; left/right confirmed exact x-mirrors). Since `width_at`/
`center_at`'s x-in/y-out contract can't honestly describe a fixed-x slice
through the curved ring stub (not a function of x there, the same reason
`SBend`/`DirectionalCoupler` forbid an interior bend cut at all),
`ring_port()`/`out_port()` are the documented escape hatch, returning the
perpendicular port's own geometry as plain floats for a caller to build a
`PortSpec` + `y0`/`y1` from. 46 new fast, meep-free geometry tests, following
`tests/test_geometry.py`'s existing style.

The body (straight bus+ring pair) gets no new geometry code — it stays a
raw-polygon `meow` call, matching how `solve_body_supermodes()` already
works in the reference scripts.

**(c) `HybridSolver` must express feedback.** This was the substantive one. The
solver still contracts a straight cascade — threading light through the
sections once, not a loop. Two call sites, both unchanged in substance since
this item opened, only their line numbers have drifted again (checked
2026-07-30): `solve()`'s `np.einsum("mij,mjk->mik", mat, full)` at
`hybrid.py:874`, and `solve_length_sweep()`'s batched counterpart —
`np.einsum("Lmij,Lmjk->Lmik", mat, full)`, an extra leading length axis, same
cascade — at `hybrid.py:973`. A ring is a **closed loop**: light re-enters
the coupling region every round trip, which is a closed-form all-pass sum,
not a matrix product. Two ways out were scoped:

- *Narrow:* keep the cascade for the open path (coupler section + bends +
  straights) and add an explicit loop-closure step that applies the resonator sum
  to a designated feedback port pair. Less work, and it matches how the script
  already assembles the answer — but it hard-codes "ring" into the solver.
- *General:* let `HybridSolver` describe an S-matrix **network** with named ports
  rather than an ordered chain, and solve it by Redheffer star products /
  node elimination. Cycles then fall out for free, and it would also cover
  Mach–Zehnders, cascaded rings, and any device where a section feeds back. This
  is the architecturally right answer and considerably more work.

**The general form's composition primitive is done, 2026-07-29 — `beef.network`
(`connect()`/`close_loop()`/`Circuit`, see `notes.md`'s Network section).**
`close_loop()` reproduces this device's own proven closed-form formula
(`E_thru = t + k1·k2·P/(1-t_ring·P)`, `ring_resonator_hybrid.py:762-774`) term
for term when applied to a 2-port `[bus, ring]` block — the concrete evidence
the general primitive is physically correct on exactly this device's math, not
just internally consistent. **Not yet wired to this example or to
`HybridSolver`**: the module depends only on `SMatrix`/`numpy` and has no
caller in `beef/solvers/` or any example — (a) and (b) above are now also
done, so nothing structural blocks the wiring anymore, but doing it (and
rewiring `HybridSolver`'s own cascade onto `beef.network`) is separate,
undone work. `close_loop()`
also only implements the same-port self-loop case (this device's case); the
two-different-ports case some `Circuit`-discovered cycles would need (e.g.
certain cascaded-ring topologies) is deferred.

**All three of (a), (b) and (c) now exist as standalone, tested library
capabilities — none are wired to each other or to this device yet.** That
wiring is the actual remaining work: build a real `RingResonatorCouplerSection`
+ `EMESolver` body + mirrored right section via `HybridSolver`/`solve_section()`
(exercising the still-unverified y-normal excite/flux-calibration path),
compose with `beef.network.close_loop()`, and check the result against this
device's own proven closed-form/full-device numbers. A new example script
would naturally follow from that, replacing the current hand-rolled scripts.

Until that wiring lands the ring cannot be validated the way the other
devices are: checks 1–3 in [README.md](README.md) all assume
`FDTDSolver`/`EMESolver` objects that this device never constructs — the
existing example scripts are unchanged and still don't call any of (a)/(b)/(c).

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
