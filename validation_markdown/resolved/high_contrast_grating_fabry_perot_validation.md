# HCG Fabry-Pérot resonator — build plan

> ## 📁 FILED INTO `resolved/` — 2026-08-08
>
> Per this directory's README, `resolved/` means *"whole docs sufficiently
> addressed (open to-dos may remain)"*. For this device that is close to
> literal: Steps 0-5 all landed with real data, and all four checks have a
> real, resolved status.
>
> | | status |
> |---|---|
> | Step 0 — platform decision | ✅ decided, `SOI_160` both pieces |
> | Step 1 — cavity length | ✅ MEASURED — `L_ARM=53.08 µm`, traced from `duty/Mux_lpbk_9$4.gds`'s real layout |
> | Step 2 — DC coupling + resolution | ✅ real coupling landed; resolution convergence MEASURED (= check 3) |
> | Step 3 — waveguide phase blocks | ✅ real 21-pt `n_eff(λ)`, bidirectional block (real bug found+fixed) |
> | Step 4 — `Circuit` assembly | ✅ matches the closed form to floating-point precision (`max\|ΔR\|=6.9e-17`) |
> | Step 5 — resolved spectrum | ✅ resolved by per-factor interpolation, not brute-force — **superseded by Step 7's numbers** |
> | Step 7 — self-contained + arms calibrated | ✅ **2026-08-18, job 11995681** — every piece solved in-script; arms moved off MPB onto the FDTD's own β (+1.463 %), comb **+4.28 FSR**, **45 resonances**, FSR **2.245 nm**, Q 9,998–16,071 |
> | Check 1 — engine agreement | ✅ CLOSED — full-band ratio stable (0.28%); stub/cut-point drift fully explained by `CALIBRATE_BODY` |
> | Check 2 — mode identity/ordering | ✅ CLOSED — identity correspondence, worst `\|dn_eff\|=0.00048` |
> | Check 3 — resolution convergence | ✅ CLOSED — converged at production `res=40/80` |
> | Check 4 — full-device FDTD | ⛔ decided OUT OF SCOPE (see A4) — a documented choice, not a gap |
>
> **⚠️ Open and NOT merely cosmetic (2026-08-18)**: the real device's measured
> FSR (2.2608 nm, `right_chip_N_LUNA`) sits ABOVE both the calibrated and
> the uncalibrated simulation, and calibrating the arms moved the simulation
> ~2x FURTHER from it. Not evidence the calibration is wrong — `L_ARM`'s own
> approximations (two real arms averaged 3.0 µm apart; bends modelled as
> straight) are of comparable size, and the arm β is unconverged — but the
> two explanations are not yet separated. See Step 7's own ⚠️ section for
> what would settle it. Manuscript FSR/FWHM figures are stale until re-run.
>
> **Open, not blocking**: the arm β is measured at ONE resolution per run
> (res 60) and nothing here converges it — `HCGFP_ARM_RES` is the hook for
> that sweep, and it is the largest remaining uncertainty on absolute comb
> position (Step 7). The coupler body (res 40) and the arms (res 60) also
> now sit on different discretization errors, a stated cost of that choice.
> Otherwise none for this device's own core physics. Two
> minor, optional items remain — dedicated pytest coverage for the newer
> interpolation scripts, and formally reformatting this doc's status into
> README's own four-checks gap-list table (this box is the informal
> equivalent). Companion notebooks
> (`examples/directional_coupler/notebooks/directional_coupler_soi160_
> hybrid.ipynb`, `examples/high_contrast_grating_fabry_perot/notebooks/
> hcg_fabry_perot_hybrid.ipynb`) were added 2026-08-08, executed end-to-end,
> and confirmed to reproduce these exact real numbers.

Opened **2026-08-05**, at the user's request: an FDTD/EME **hybrid** example
of a Fabry-Pérot cavity terminated by two high-contrast-grating (HCG)
mirrors, tapped in the middle by a directional coupler (DC), per the
diagram supplied with the request. Reproduced in text form since the
source is a user-supplied PDF, not a repo file:

```
 "input/        ______________                    ______________
  reflection"  /              \                  /              \   "output"
      --------|                \________________/                |--------
               \  (DC lower arm, coupled straight body)          /
                \______________/                  \______________/
   FDTD    EME          FDTD          EME (DC)         FDTD          EME    FDTD
  (HCG)  (waveguide)   (DC bend)     (coupled pair)   (DC bend)   (waveguide) (HCG)
   |<---------------------------- Fabry-Pérot cavity (upper rail) --------------------->|
```

The **upper rail** is the resonant cavity itself: HCG mirror — straight
single-mode waveguide — DC's "through" arm (bending down into the coupling
region and back up) — straight single-mode waveguide — HCG mirror, with no
break in the guided path. The **lower rail** is the DC's other two ports
("input/reflection" and "output") — physically disconnected from the
cavity, coupled to it only evanescently across the DC's straight coupled
region. Light injected at "input/reflection" partially cross-couples into
the cavity, bounces between the two HCG mirrors (a genuine Fabry-Pérot
resonance, since both ends are high-but-not-unity reflectors), and leaks
back out at both "input/reflection" (net cavity reflection) and "output"
(net cavity transmission, tapped off evanescently on every round trip).

> **This doc is a PLAN, not a gap list**, same convention as
> `mmi_1x4_sin300_validation.md`: nothing below has been simulated for this
> combined device. Every quoted number is either read directly from
> existing code (labeled **MEASURED (read from source)** — these are
> facts about what already exists, not about this new device) or is a
> **PROPOSED** design choice / **ARGUED** prediction to be checked once the
> pieces actually run. Converts to a normal gap list against the
> [four checks](../README.md#the-four-checks) once a hybrid result exists.
>
> **Two pieces already exist and already work**:
> `examples/high_contrast_grating/python/high_contrast_grating_fdtd.py`
> (the HCG mirror, `FDTDSolver.solve_reflection()`, see
> `high_contrast_grating_validation.md`) and
> `examples/directional_coupler/python/directional_coupler_hybrid.py` (the
> DC, `FDTDSolver.solve_section()` + `EMESolver.solve_section_length_sweep()`).
> **Nothing currently cascades independently-solved devices together** —
> that composition layer (`beef.network`) exists and is unused by any
> example script today. This doc's real job is wiring, not new physics.
>
> **Updated 2026-08-05 — Step 0 decided, Step 1 pinned.** The platform
> mismatch (A1) is resolved: reuse the HCG mirror on `SOI_160` as-is,
> redesign the DC for `SOI_160` at a REAL fabricated coupler's own
> dimensions (`duty/dc.gds`, read directly with `gdstk`, fed into the
> existing parametric `DirectionalCoupler` class — not a `GdsGeometry`
> import). **That real device's cross-section supports 4 supermodes, not
> 2** — a measured finding — **but a follow-up measurement (exciting the
> real S-bend with mode 1, decomposed at the coupled-region interface)
> shows only 0.0085% of the launched power reaches the 2nd-order-mode
> pair.** User's call: treat the device as single-mode, `N_MODES=2`, a
> measured simplification, not an assumed one. See Step 0. Step 1 (cavity
> length/FSR) is pinned at 20 µm/arm (~9.5 nm FSR estimated from a real
> measured group index) — provisionally, per the user's own framing:
> revisit once Step 4 actually simulates something and the real resonance
> spacing can be judged directly. **Step 2 is running, on its 3rd
> submission after two timeouts**: measures the real coupling at the
> GDS's fixed 3.0 µm `coupling_length`, `N_MODES=2`, at
> `FDTD_RESOLUTION=40`/`EME_RESOLUTION=80` ASSUMED from the `SOI_220`
> script (user's explicit instruction — not re-derived for this
> platform/device this pass). **Job 11691379 (`-t 300`) timed out** — 5 of
> its 6 expensive FDTD sub-solves had already completed and cached
> cleanly (confirmed from the log and `example_sim_cache/`, both arm
> excitations plus all 3 of `section_mode_beta()`'s calibration-length
> sub-solves). **Resubmit job 11706987 (`-t 240`) ALSO timed out** —
> despite all 5 expensive pieces correctly cache-hitting this time (no
> matching re-computation in its log), the whole 4 h budget was still
> consumed by something inside two nominally-cheap `section_mode_neff()`
> calls, whose real cost driver is not yet fully diagnosed. **Shown this,
> the user's call was to resubmit with a much larger margin rather than
> keep investigating — job 11716260, `-t 720` (12 h).**
>
> **While that runs: Steps 3-5's DESIGN is now done and verified with
> synthetic data** (`tests/test_fabry_perot_
> network_crosscheck.py`) — including a real bug caught before any real
> device numbers touched it: the connecting-waveguide block cannot be
> built via `beef.network.convert.from_transfer_matrix()` (one-directional
> by design; this cavity needs the waveguide to carry the mirror's
> reflection back through itself) — fixed with a direct bidirectional
> `SMatrixBlock`. The full 5-block `Circuit` wiring and a hand-derived
> Fabry-Pérot-with-DC-tap closed form are cross-checked against each other
> and against a brute-force solve, all to floating-point precision.
> **Step 3 then landed for real too**: a real 21-point fundamental-mode
> `n_eff(λ)` sweep for the connecting waveguide (matching the HCG's own
> production grid), connected to the HCG's own real mirror block to give
> `rho(λ)` — the actual input Step 5's closed form needs. Only Step 2's
> real `t(λ)`/`k(λ)` are still missing before a fully real cavity spectrum
> can be evaluated.
>
> **Updated 2026-08-07 — job 11716260 landed (4h56m real walltime), and
> the first fully real cavity spectrum has been evaluated.** DC coupling
> at the fixed `coupling_length=3.0 µm`: `T_through=0.8989`,
> `T_cross=0.0979` at λ=1.55 µm (≈10% coupling ratio). Combined with the
> HCG's real mirror spectrum and Step 3's real waveguide phase through the
> verified closed form (Step 5): **a genuine Fabry-Pérot resonance shows
> up at λ=1.570 µm** — transmission drops from a flat ~0.97–0.99 baseline
> to 0.7395, reflection spikes from ~0.001–0.02 up to 0.1710, both moving
> together as real interference, not noise. **But the 11-point/10 nm-step
> wavelength grid is comparable to the ~9.5 nm estimated FSR, so this
> resonance is almost certainly under-sampled** — this result establishes
> the assembled model DOES produce real resonant behavior at a physically
> plausible location, not that the resonance lineshape (depth/width/Q/
> exact position) is actually characterized yet. See Step 5.
>
> **A linewidth sweep was submitted to actually resolve it**: sized from
> the real numbers above (finesse ~23, FWHM ~0.41 nm), job 11726827 swept
> 1.555–1.585 µm at 0.05 nm spacing (601 points, ~8 across the estimated
> linewidth) — not a guess, a calculation from job 11716260's own measured
> `|rho|²`/`|t|²`. **TIMED OUT 2026-08-08 at the full `-t 720` (12 h)
> budget** — the wider window doesn't cache-hit job 11716260's saved
> pieces, so every expensive FDTD sub-solve had to redo, each now also
> carrying a 601-frequency (not 11) decomposition tail; only 3 of the
> needed pieces finished in 12 h (two ~2.6 h arm excitations, one ~4.9 h
> calibration-length probe). A real, measured cost blowup, not the
> "~4h56m plus a cheap tail" originally guessed. Not yet resubmitted —
> genuine cost/design trade-off (bigger budget vs. fewer points vs.
> narrower window), put to the user. See Step 5.
>
> **Worked while that ran: A4/Step 6 updated with the ring resonator's own
> precedent** (its full-device FDTD check was retired 2026-08-03 for
> reporting T>1 in a passive device and a non-converging comb — a
> disqualification, not just a cost problem — for a mechanistic reason
> that generalizes: round-trip phase error is amplified ~F× by finesse).
> At this device's own F≈23 (vs. the ring's F≈120) that threshold is 5×
> more forgiving, but the HCG's grating mirrors staircase more per bounce
> than the ring's smooth bends — net direction unmeasured either way.
>
> **User's call, 2026-08-07: skip it outright — no full-device FDTD control
> run for this device, planned or attempted.** Not the cheap-pilot-first
> option floated initially — decided directly, given the likely
> disqualification (by the argument above) stacked on the likely cost (the
> ring's 167 h is a floor for a smaller, lower-finesse device, not a
> ceiling for this one). This device's validation rests entirely on the
> single-pass/algebraic checks already real: Step 4's closed-form-vs-
> `Circuit` match to floating-point precision, the `k'=k`/`t'=t`
> reciprocity check, and passivity at every evaluated point. See A4,
> Step 6.
>
> **Updated 2026-08-08 — the resonance is now FULLY RESOLVED, and NOT via
> job 11726827.** That job timed out at its full 12 h budget (see Step 5)
> after completing only 3 of the ~7+ needed expensive FDTD sub-solves —
> brute-force dense re-sampling was the wrong shape for this problem.
> Fixed instead with interpolation: per-factor phase analysis of the REAL
> cached 11-point data showed the mode-basis pieces (`M_in0`/`M_in1`,
> the body's calibrated `beta(λ)`) are each comfortably safe to
> interpolate (~65-69°/10nm and ~15-16°/10nm respectively) — only their
> COMBINATION looked dangerous (~153°/10nm), because this device's own
> mirror-symmetry reuse counts one S-bend's phase twice. This is standard
> practice, not a special-case trick (confirmed against Lumerical's own
> EME documentation and coupled-mode theory, both describing exactly this
> "solve modes sparsely, handle propagation phase analytically" split —
> see Step 5). **Real validation, not just an argument**: 2 of job
> 11726827's own pieces (both 601-point arm excitations) finished before
> its timeout, and interpolating the sparse 11-point data reproduces that
> real measurement to ≤1.3% error. **Result: 22 real resonances resolved
> across the C-band, mean FSR 4.500 nm** (Step 1's original ~9.5 nm
> estimate had a missing round-trip factor of 2 — corrected, see Step 1).
> The resonance near the original 11-point sample (T=0.7395 at 1.570 µm)
> turned out to be a real dip, but ~18× shallower than its true depth
> (T_min=0.0415 at λ=1.5698 µm, FWHM=0.235 nm, Q≈6680) — confirming the
> under-sampling concern was real. See Step 5 for the full table.
>
> **Updated 2026-08-08 — Step 1's cavity length is now PINNED FOR REAL,
> not provisional.** The 20 µm/arm placeholder is replaced with 53.08 µm
> (average of two real arms), traced instance-by-instance from
> `duty/Mux_lpbk_9$4.gds` — a real fabricated test structure sharing this
> device's own mirror-DC-mirror topology — at the user's direction, bends
> treated as straight waveguide (bend radius `R=10.0 µm`, confirmed two
> independent ways: from the bend cell's own polygon extent and from its
> port-to-port chord length). Re-running the fully-resolved dense
> reconstruction at this real length (the interpolation scheme itself is
> unaffected — only the connecting-waveguide phase term changes) gives
> **44 real resonances across the C-band, mean FSR 2.255 nm, Q up to
> ~16000** — roughly double the resonance count and Q of the 20 µm case,
> exactly as expected for a ~2.05× longer real cavity. See Step 1 and
> Step 5 for the full derivation and tables.
>
> **Updated 2026-08-08 — Step 6, check 2 (mode identity/ordering) is now
> CLOSED, measured.** Extended `examples/verification/python/
> audit_interface_modes.py` (the same cheap, no-time-stepping tool already
> used for the SOI_220 coupler and both MMIs) with this device's own real
> geometry. Result: both engines return 2/2 guided modes, **identity
> correspondence, worst `|dn_eff|=0.00048`** — no mode-mixing risk found.
>
> **Same day — check 1 (engine agreement) partially closed too, at zero
> new compute.** `section_mode_beta()` already computes a per-wavelength
> FDTD/MPB comparison internally but only returns the aggregate result,
> discarding the rest. Recovered the full grid by reusing already-cached
> calibration phases with the library's own `beta_fit_from_phases()`:
> the `Δβ_FDTD/Δβ_MPB` ratio is **0.8801 ± 0.00084 across the entire
> 100 nm band** (0.28% relative variation) and every phase-fit residual is
> `≤9.3×10⁻⁴ rad` (>100× under the guard threshold) — the ~12%
> discrepancy `CALIBRATE_BODY` compensates for is a smooth, stable,
> physically real effect, not a one-wavelength fluke. This does NOT
> replicate the SOI_220 coupler's own stronger stub/cut-point invariance
> test, which would need genuinely new FDTD (real cost, not attempted).
> Check 3 (resolution convergence) remains open, explicitly deferred;
> check 4 stays out of scope by decision. See Step 6.
>
> **User's call, 2026-08-08: fund both remaining real-cost checks.**
> `directional_coupler_soi160_hybrid.py` gained `DC160_BODY_STUB` (check 1's
> cut-point sweep) and `DC160_RES_SCALE` (check 3's resolution sweep)
> env-var overrides — both force `N_WL=1` at λ=1.55 µm (these test a
> numerical property, not the spectrum, so paying for the full 11-point
> decomposition would buy nothing new). Submitted as a 5-task SLURM array,
> **job 11747764** (`batch/dc_soi160_hybrid_checks_1_3.batch`, `-t 720`
> each): `BODY_STUB ∈ {0.1, 0.5, 0.75}` (check 1, alongside the already-real
> `0.25` production point) and `RES_SCALE ∈ {0.75, 1.25}` (check 3,
> bracketing the production `1.0`). Cost is unmeasured at `N_WL=1` — the
> underlying 5-real-FDTD-sub-solve cost this device's own history shows
> (~5 h at `N_WL=11`) is not known to shrink just because fewer
> frequencies are monitored, so `-t 720` reuses the same working margin as
> every other real measurement here, not a measured bound. See Step 6.
>
> **Job 11747764 landed, 2026-08-08 (real costs: 37 min to 4h55m per task,
> rising steeply with resolution).** **Check 3: ✅ MEASURED, converged at
> the production resolution** — `res=40→50` moves `T_cross` by 0.00%
> (`0.0979` both), `res=30→40` moves it +7.35% — the production
> `FDTD_RESOLUTION=40`/`EME_RESOLUTION=80`, carried over unverified from
> the SOI_220 device since 2026-08-06, is independently confirmed correct.
> **Check 1: ✅ CLOSED, and resolved for real, 2026-08-08** — `T_cross` is
> stable to <1% across `BODY_STUB=0.1–0.75 µm`. The naive adapted formula
> predicted a ~23× larger drift, but that formula was borrowed from the
> SOI_220 coupler's own PRE-calibration stub sweep (derived to expose the
> Δβ gap using MPB's uncalibrated body rate) — this device already runs
> `CALIBRATE_BODY=True`, so the comparison was apples-to-oranges. Checked
> directly instead: the arm sections' own real FDTD phase-accumulation
> rate (branch-safe fit via the library's own `beta_fit_from_phases()`,
> avoiding the wrap-ambiguity a naive unwrap risks) matches the isolated
> calibration's rate to 0.3% per mode — meaning the stub-dependent terms
> in the combined coefficient cancel almost exactly BY CONSTRUCTION,
> because `CALIBRATE_BODY` makes the body propagate on that same measured
> rate. The small residual drift is fully explained by that 0.3% gap —
> nothing left open. All four checks now have a real, resolved status.
> See Step 6.

---

## Summary — where the work actually is

| | status |
|---|---|
| HCG mirror, FDTD reflection (complex, with phase) | ✅ already exists — `FDTDSolver.solve_reflection()`, `SOI_160`, C-band |
| Directional coupler, FDTD/EME hybrid | ✅ **redesigned for `SOI_160` and re-solved for real (Step 2)** — `T_through=0.8989`, `T_cross=0.0979` at λ=1.55 µm, coupling ratio 0.098 |
| Both on the same material platform | ✅ **decided 2026-08-05** — both `SOI_160`, DC using `duty/dc.gds`'s real dimensions |
| DC cross-section is single-mode | ⚠ **cross-section supports 4 supermodes, but MEASURED <0.01% of excited power reaches the 2nd-order pair** — treated as `N_MODES=2` by user decision (Step 0), a measured simplification |
| Cavity geometry (mirror spacing, DC coupling ratio, target finesse) | ✅ **pinned FOR REAL (Step 1)** — coupling_length fixed by the real GDS (3.0 µm); waveguide length 53.08 µm/arm, traced from `duty/Mux_lpbk_9$4.gds`'s real layout; FSR MEASURED at 2.255 nm (44 real resonances, Step 5) |
| Cross-device S-matrix composition (mirror + waveguide + DC + waveguide + mirror) | ✅ **DONE with real device data (Step 4)** — `examples/high_contrast_grating_fabry_perot/python/hcg_fabry_perot_hybrid.py`, matches the closed form to floating-point precision |
| A 1-port self-reflecting block (the HCG, in network terms) | ✅ **DONE with real data (A3/Step 4)** — a direct `SMatrixBlock` construction from the real production `SMatrix`, not via `from_smatrix()` (which cannot represent it) |
| The connecting-waveguide block | ✅ **real bug found and fixed (Step 3)** — `from_transfer_matrix()` is one-directional by design and silently zeroes the mirror's reflection; fixed with a direct bidirectional `SMatrixBlock`. Solved for the real cross-section, both sparse (11 pt) and reused (fit) for the dense grid |
| Closed-form sanity check (Airy/Fabry-Pérot formula) | ✅ **FULLY RESOLVED with real numbers, at the real cavity length (Step 5)** — 44 real resonances, mean FSR 2.255 nm, real FWHM/Q/depth per resonance (Q≈9990-16000), via interpolation of the safe intermediate transfer-matrix pieces (validated to ≤1.3% against real partial ground truth), not a brute-force dense FDTD re-solve |
| Full-device FDTD control run | ⛔ **decided 2026-08-07: out of scope, not attempted** — see A4/Step 6 |

So the ordered work is: **resolve the platform mismatch (Step 0, ✅ done)
→ pin cavity geometry (Step 1, ✅ pinned, FSR now measured) → re-derive
the DC for `SOI_160` at `N_MODES=2` (Step 2, ✅ all items done, real
coupling AND resolution convergence both measured) →
build the two EME waveguide phase blocks (Step 3, ✅ real data landed) →
assemble with `beef.network.Circuit` (Step 4, ✅ done with real data,
matches the closed form to floating-point precision) → closed-form
cross-check (Step 5, ✅ FULLY RESOLVED — 44 real resonances at the real
cavity length, via interpolation, not brute-force re-solving) → the four
checks (Step 6, ✅ checks 1-3 all CLOSED with real measurements, check 4
⛔ decided out of scope — no full-device FDTD control run for this
device).**

---

## Part A — findings that set the plan

### A1. The two existing devices are on incompatible material platforms — MEASURED (read from source), blocking

`examples/high_contrast_grating/python/high_contrast_grating_fdtd.py` uses
`SOI_160` (`beef.materials.SOI_160`, 160 nm Si core) with a 500 nm-wide
access waveguide (`GdsGeometry(..., input_width=0.5)`, matching the GDS's
own drawn cross-section). `examples/directional_coupler/python/
directional_coupler_hybrid.py` uses `SOI_220` (220 nm Si core) with a
450 nm single-mode strip (`waveguide_width=0.45`).

A single-mode strip's cross-section is not a free parameter once the core
thickness is fixed — 450 nm is single-mode at 220 nm thickness, not
necessarily at 160 nm (and vice versa for 500 nm at 220 nm). **There is no
version of this combined device today where the mirror and the coupler
agree on what "the waveguide" physically is.** Three ways to resolve this,
none simulated yet:

1. **Redesign the DC for `SOI_160`** at a new single-mode width, and re-run
   `directional_coupler_hybrid.py`'s whole beat-length characterization on
   the new platform (this repo's own convergence/resolution numbers for the
   DC are `SOI_220`-specific — `FDTD_RESOLUTION=40`/`EME_RESOLUTION=80`
   were measured on that stack, per that script's own comments referencing
   `dc_validation.md item 11`; nothing says they transfer).
2. **Redesign the HCG mirror for `SOI_220`** at the DC's existing 450 nm
   width, and re-run the HCG's own convergence sweep on the new platform
   (`high_contrast_grating_validation.md`'s `res=60`/PML-clearance/etc.
   numbers are `SOI_160`-specific to this exact device's ~7 µm domain —
   nothing says a 220 nm core changes them, but nothing says it doesn't
   either).
3. **Keep both platforms and add a thin, lossless-in-principle mode
   converter** between the HCG's 160 nm/500 nm access waveguide and the
   DC's 220 nm/450 nm strip — physically realistic (real photonic chips do
   this at layer transitions) but adds a THIRD FDTD simulation type this
   doc's Part B does not currently budget for, and a new loss/reflection
   source at exactly the interface the whole point of this exercise is to
   characterize cleanly.

**This doc does not pick one.** Option 1 is the smallest number of new
FDTD/EME solves (the HCG mirror is reused as-is); option 2 keeps the DC's
already-measured convergence numbers but requires re-running the HCG's own
resolution/PML sweep (recent, expensive — see
`high_contrast_grating_validation.md` step 8, ~75-90 min/point at
`res=60` alone) on a platform nobody has tested it on. **Recommend option
1** on cost grounds alone, but this is the user's call, not something to
assume silently and discover wrong after both hybrid halves are built.

### A2. `beef.network` is the right composition layer, and no example uses it across devices yet — MEASURED (read from source)

`beef/network/` (`SMatrixBlock`, `connect()`, `close_loop()`,
`close_loop_pair()`, `Circuit`) is a general dense-multi-port S-parameter
composition layer, distinct from `beef.results.SMatrix` (a flat, sparse
per-excitation dict with no notion of a port's own reflection —
`block.py`'s own docstring is explicit about this). Every existing hybrid
example (`directional_coupler_hybrid.py`, `ring_resonator_hybrid.py`,
`mmi_2x2_hybrid.py`, …) composes sections **within one device** via either
`HybridSolver`'s own matrix-chain multiplication or (for the ring's
cascaded variant) a hand-rolled closed-form assembly — never via
`beef.network.Circuit` connecting the *outputs of independently-solved
devices*. `tests/test_network.py`/`test_network_ring_crosscheck.py` prove
the primitives are individually correct (including a cross-check that
`close_loop()` reproduces the ring's own closed-form all-pass formula
exactly). One non-test script,
`examples/verification/python/verify_ring_network_crosscheck.py`, already
feeds `close_loop()`/`connect()` real `FDTDSolver.solve_section()`
coefficients — but only for the ring resonator's own coupler section, a
single physical device, and only via the free functions directly, never
the `Circuit` class. **Nothing anywhere uses `Circuit.add_block()` to
compose more than one independently-solved physical device together, and
nothing feeds this layer a genuinely 1-port self-reflecting block (A3).
This device is the first real user of both.**

`Circuit` is exactly this doc's target shape: `add_block(name, block)` for
each of the five physical elements (HCG-left, waveguide-left, DC,
waveguide-right, HCG-right), `connect(name_a, port_a, name_b, port_b)` for
each of the four internal seams, and `reduce(external_ports=...)` to get
back a single `SMatrixBlock` exposing exactly "input/reflection" and
"output" — no loop-closing needed, since the cavity here is a **finite,
mirror-terminated line**, not a closed ring (contrast
`ring_resonator_hybrid_cascaded.py`, which *does* need `close_loop()`
because its ring genuinely closes on itself).

### A3. The HCG mirror is a 1-port self-reflecting block — not what `from_smatrix()` assumes — MEASURED (read from source)

`beef.network.convert.from_smatrix()`'s own docstring: *"Self-terms
(input_port's own reflection, and cross-terms between two different output
ports) are zero-filled — an honest limitation of what SMatrix records, not
an approximation this function introduces."* It requires `input_port` and
`output_ports` to be given as **distinct** labels and always zero-fills the
excited port's own self-term.

`FDTDSolver.solve_reflection()`'s `SMatrix` has exactly one entry: `port="in",
in_mode=1, out_mode=1` — the *complex* reflection coefficient
`r(λ) = bwd/fwd` (phase included, confirmed by reading
`beef/solvers/fdtd.py`'s `solve_reflection()` directly: it stores the raw
ratio, not `|r|²`). This is precisely a port reflecting into **itself** —
the thing `from_smatrix()` is explicitly not built to carry. Routing it
through `from_smatrix()` with a different `input_port` label would
misrepresent the mirror as a 2-port transmissive device (power going in
one labeled port and out a different one, with **zero** self-reflection at
either) — exactly backwards from what a mirror does.

**Not a library bug** — `from_smatrix()` is honest about the limitation —
but a real integration step: build the HCG's block directly,
bypassing `from_smatrix()`:

```python
r = hcg_smat.coeff(port="in", in_mode=1, out_mode=1)   # [n_wl] complex
hcg_block = SMatrixBlock(
    wavelengths=hcg_smat.wavelengths, ports=("p",),
    modes_per_port={"p": 1}, data=r[:, None, None],
)
```

`SMatrixBlock`'s own validation (a dense `[M, N, N]` array, `N =
sum(modes_per_port.values())`) places no lower bound on port *count* — a
genuine 1-port block is structurally valid, just untested by
`tests/test_network.py` today (worth a small new test alongside whatever
implements this, per Step 4/6 below).

### A4. A full-device FDTD control run — ⛔ DECIDED 2026-08-07: not running one for this device

Every other hybrid device in this repo (DC, ring, MMIs) has, or plans, a
full-device FDTD run as an independent cross-check on the hybrid
decomposition (check 4 of the [four checks](../README.md#the-four-checks)).
For a Fabry-Pérot cavity, the whole point of the device is a HIGH-Q
resonance — `ring_resonator_fdtd.py`'s own validation doc already reports
**167 h estimated for a 3-D `res=40` full-device FDTD** on a racetrack ring
with a much smaller total path length and far lower finesse than a
double-HCG-mirror cavity is likely to have. A Fabry-Pérot cavity's `Q`
scales with mirror reflectivity — the HCG's own measured R≈97% per mirror
(`high_contrast_grating_validation.md` step 9) implies a finesse in the
tens, meaning the cavity's `decay_tol`-driven ringdown time in a full-device
FDTD run would be tens of round trips long, at whatever cavity length Step
1 picks. **This should be measured, not assumed impossible outright** —
state it as a real open cost question in Step 6, not silently skip the
check.

**Updated 2026-08-07 — ARGUED, by direct precedent, not yet measured for
this device: cost is very likely NOT the only problem.**
`validation_markdown/resolved/ring_resonator_validation.md` retired its own
full-device FDTD check outright (2026-08-03), for a reason that has nothing
to do with wall-clock budget: the reference itself reported **transmission
above 1 in a lossless passive device** (up to +8.4%, over 46% of the band,
**worsening** from res 30→40 rather than refining away), and its own
resonance comb shifted **18% of an FSR** between those same two
resolutions — i.e. the "ground truth" carried more uncertainty than the
thing it was supposed to arbitrate. Their diagnosis, stated as a general
mechanism, not specific to the ring: *"Finesse `F = Q·FSR/λ`. The field
circulates ~`F` times, so any per-round-trip phase error is amplified
~`F`×... FDTD's staircased arcs give O(h) boundary-position error, which
*is* a round-trip path-length error. The reference is therefore subject to
exactly the amplification it is supposed to be checking."* Their threshold:
at their `F≈120`, a **0.5% round-trip error moves a resonance a full
linewidth**.

This device's own measured finesse (Step 5's follow-up calculation) is
**F≈23** — about 5× more tolerant (threshold ≈1/F≈4.3% round-trip error
before a resonance moves a full linewidth, vs. the ring's ≈0.83%), so this
is NOT a slam-dunk "definitely also disqualified" the way it was for the
ring's much higher finesse. But the underlying mechanism (round-trip phase
error scales with mirror count × per-bounce staircasing error, amplified by
`F`) applies to any resonant cavity validated by time-domain FDTD,
regardless of device — and this cavity's mirrors are HCG **gratings**, each
presenting many periodic staircased tooth-edges per round trip, plausibly
*more* discretization surface per bounce than the ring's two smooth bus/ring
coupler arcs, not less. Which direction dominates (lower F helping vs. more
staircased edges per bounce hurting) is a real open question, **not
resolved by analogy alone**.

**Decision (user), 2026-08-07: skip check 4 for this device — no full-device
FDTD control run, planned or attempted.** Not deferred pending the cheap
2-D pilot floated above — decided outright, the same way the ring
resonator's own doc eventually landed (its check 4 was RETIRED after
actually running and failing; this one is scoped out before spending any
compute on it at all, since the mechanism argued above already predicts
the same class of failure and the 167 h ring number is a floor, not a
ceiling, for what a full double-HCG-mirror cavity would cost). Consequence
for Step 6: **this device's validation rests entirely on the single-pass/
algebraic checks that carry no round-trip amplification** — Step 4's
closed-form-vs-`Circuit` match to floating-point precision, the
`k'=k`/`t'=t` reciprocity check, passivity at every evaluated point, and
each individual device's own existing checks (HCG mirror, DC coupling).
No full-device time-domain reference will exist for this cavity, by
choice, not by oversight — record it as a deliberate scope decision in
Step 6, mirroring how `ring_resonator_validation.md` documents its own.

---

## Part B — the ordered steps

### Step 0 — resolve the platform mismatch (A1) — ✅ DECIDED 2026-08-05

**Decision (user): option 1 from A1** — reuse the HCG mirror as-is on
`SOI_160`, redesign the DC for `SOI_160` too, and use a REAL fabricated
coupler's own dimensions rather than a freshly-derived cross-section:
`duty/dc.gds` (top-level `duty/`, not `BEEF/examples/gds_files/` —
different from the HCG's own GDS). **Dimensions only, not a `GdsGeometry`
import** — read directly from the GDS's layer-(2,0) polygons with
`gdstk`, then fed into the existing parametric `DirectionalCoupler` class
(matching every other hybrid example's approach; nothing about this class
requires or benefits from a raw-polygon import path):

| parameter | value | how it was read |
|---|---|---|
| `waveguide_width` | 0.5 µm | coupled-region arm width, both arms |
| `gap` | 0.25 µm | coupled-region edge-to-edge separation |
| `coupling_length` | 3.0 µm | coupled-region x-span |
| `bend_length` | 10.0 µm | each S-bend's x-span (left and right, confirmed mirror-symmetric) |
| `waveguide_offset` | 4.0 µm | far-port centre-to-centre pitch (2.375 − (−1.625), confirmed symmetric bend shapes) |
| `bend_type` | `"euler"` (assumed) | not re-derived from the GDS's curve family — a modeling choice, not a measured one; matches the existing `SOI_220` script's own default |
| `waveguide_height` | 0.16 µm (`SOI_160`) | the GDS itself carries no thickness — supplied externally, same as the HCG needed `input_width`/`input_y` supplied externally |

Verified against `beef/geometry/coupler.py`'s exact constructor semantics
(`gap`=edge-to-edge, `waveguide_offset`=far centre-to-centre pitch,
`coupling_pitch`=derived=`waveguide_width+gap`) before reading these
numbers — confirmed `DirectionalCoupler(waveguide_width=0.5, gap=0.25,
coupling_length=3.0, bend_length=10.0, waveguide_offset=4.0)` builds
without error (`total_length=23.0 µm`, `coupling_pitch=0.75 µm`).

**Real finding that changes the plan: this waveguide is NOT single-mode at
160 nm in C-band — MEASURED** (cheap MPB mode-count check, no FDTD, ~9-11
min on the login node at `res=100`, `OMP_NUM_THREADS=1`):

Bare 0.5 µm-wide, 0.16 µm-thick strip, guided modes (real n_eff above the
SiO2 cladding index):

| λ (µm) | mode 1 (fundamental) | mode 2 (2nd order) | cladding n |
|---|---|---|---|
| 1.50 | 2.2113 | 1.5188 | 1.4446 |
| 1.55 | 2.1507 | 1.5039 | 1.4440 |
| 1.60 | 2.0912 | 1.4926 | 1.4434 |

Both modes sit with real margin above cladding at every C-band point
checked — this is genuinely 2-mode, not a marginal/near-cutoff case.
Coupled-pair cross-section (both arms present, `gap=0.25`) at λ=1.55 µm:

```
n_eff = [2.1628, 2.1366, 1.5396, 1.4500, 1.4430, 1.4327]   (cladding n=1.4440)
```

4 bands sit clearly above cladding (the fundamental-mode even/odd
supermode pair near 2.15-2.16, and the 2nd-order-mode even/odd supermode
pair near 1.45-1.54); band 5 (1.4430) sits at/just below the cladding
index (1.4440) — the same "marginal, drop it" call the `SOI_220` script's
own N_MODES comment makes for its own band 5. **So the cross-section
itself supports 4 supermodes, not 2.**

**But a follow-up real measurement shows the S-bend barely excites the
2nd-order-mode pair at all — MEASURED, 2026-08-05, `res=25`,
`decay_tol=1e-10`, one real `FDTDSolver.solve_section()` run exciting the
far access port with mode 1 only, decomposed onto all 4 supermodes at the
coupled-region interface**:

| output supermode | power fraction |
|---|---|
| 1, 2 (fundamental-mode pair) | 0.4965, 0.4953 (99.18% combined) |
| 3, 4 (2nd-order-mode pair) | 5.5×10⁻⁵, 2.9×10⁻⁵ (0.0085% combined) |

(total 0.9919 — the small shortfall from 1 is ordinary radiation/PML
loss, not a red flag.) **Less than a hundredth of a percent of the
launched power reaches the 2nd-order-mode supermodes** — the large
index gap between the two modes (Δn≈0.65 at 1.55 µm) suppresses
bend-induced conversion at this bend length exactly as expected, and now
this is measured, not argued.

**User's call on this finding (2026-08-05): treat the device as
single-mode, `N_MODES=2`, not 4.** The cross-section is genuinely
multimode (bands exist), but the physics this device's own S-bend
actually excites is, to measured precision, single-mode — so the
simpler, cheaper 2-supermode treatment (matching the `SOI_220` script's
own structure, though NOT its specific numbers — different device, do not
reuse `N_MODES=5` or its resolution/coupling-length values) is used from
here forward. This is a measured simplification, not an assumed one: the
0.0085% figure above is the actual justification, re-derived for this
exact geometry rather than inherited from the general "large Δn suppresses
conversion" argument alone.

### Step 1 — pin cavity geometry — ✅ PINNED FOR REAL 2026-08-08 (was provisional 2026-08-05, now from a real fabricated test structure)

`coupling_length` is no longer a free knob (Step 0: fixed at the GDS's own
3.0 µm, not re-derived off a 50:50 target — the user's choice to use the
real device's own dimensions applies here too, not just to the
cross-section). What Step 0 does NOT fix:

- **Mirror-to-DC waveguide length: PINNED at 53.08 µm each side (average),
  2026-08-08 — MEASURED from a real fabricated test structure, not
  argued.** Originally a placeholder (20 µm, 2026-08-05, explicitly
  provisional — "use 20 µm for now, revisit once the notebook actually
  simulates something real"). Once the resonance was fully resolved
  (Step 5), the user pointed at `duty/Mux_lpbk_9$4.gds` — a real
  fabricated test/loopback structure with this exact topology (two
  reflector elements either side of a directional coupler) — and asked
  for the real waveguide length it uses between the DC and the mirrors,
  treating bends as straight waveguide (matching this doc's own existing
  simplification: the connecting arm is already modeled as a single
  straight bare-strip mode, Step 3, so a bend's own arc length slots in
  as ordinary added length rather than requiring a new bend-mode model).

  **Traced instance-by-instance** (`gdstk`, reading each reference's own
  `origin`/`rotation`/`x_reflection` and chaining port-to-port, the same
  method Step 0 used to read `dc.gds`): each arm is `straight → bend →
  straight`. Bend radius `R=10.0 µm`, confirmed two independent ways —
  from the bend cell's own core-layer polygon extent, and from each
  bend's own port-to-port chord length (`chord = R√2` for a 90° arc,
  measured chord `14.142 µm` ⟹ `R=10.000`) — both agree exactly. Arc
  length (quarter circle, per the "treat bends as straight" instruction):
  `(π/2)·10.0 = 15.708 µm`.

  | arm | straight 1 | bend | straight 2 | total |
  |---|---|---|---|---|
  | left | 25.375 µm | 15.708 µm | 13.5 µm | **54.583 µm** |
  | right | 25.375 µm | 15.708 µm | 10.5 µm | **51.583 µm** |

  The two real arms differ by exactly 3.0 µm (a real layout asymmetry —
  the two reflectors sit at the same distance from the coupler's own
  centerline, but route to different coupler ports whose own internal
  arms differ by that much). This device's own symmetric-cavity
  assumption (one `EMESolver` solve serving both arms) does not capture
  that asymmetry — **`L_ARM = 53.08 µm`, the average of the two real
  arms**, is used as the single representative length. Implemented in
  `hcg_fabry_perot_hybrid.py`'s own `L_ARM` constant (see its comment for
  the full derivation); no longer provisional.
- **DC coupling ratio at this device's real, fixed geometry** — NOT a
  design knob here (contrast the `SOI_220` script, which sweeps
  `coupling_length` to find its own 50:50 point). With `coupling_length`
  fixed at 3.0 µm, the actual power-coupling ratio (fundamental-mode
  even/odd split, per Step 0's `N_MODES=2` decision) at this length is an
  OUTCOME of Step 2's FDTD/EME solves, not a target to hit — measure it,
  don't assume 50:50 or any other specific split.
- **Target resonance wavelength / FSR** — not independently pinned; a
  first estimate below, from a REAL measured group index (not a guess).

**FSR estimate — ARGUED, not yet a real EME solve of the full cavity**,
built from the Step 0 bare-strip fundamental-mode data (2-point-plus-one
finite difference across the 3 measured wavelengths, `n_g = n_eff −
λ·dn_eff/dλ`):

| λ (µm) | n_eff (fundamental) | n_g |
|---|---|---|
| 1.50 | 2.2113 | 4.0293 |
| 1.55 | 2.1507 | 4.0122 |
| 1.60 | 2.0912 | 3.9952 |

At the proposed 20 µm/arm: `L_cavity = 2×20 + dc.total_length (23.0) = 63
µm`. `FSR ≈ λ²/(n_g·L_cavity)` at λ=1.55 µm, n_g=4.0122: **FSR ≈ 9.5 nm**,
~10-11 resonances across the 1.50-1.60 µm C-band sweep — a similar density
to the ring resonator's own C-band characterization (`ring_resonator_
hybrid.py` predicts ~10 resonances at its own shipped `PERIMETER`).
Reasonable to characterize with, not so dense that individual resonances
would be hard to resolve at this repo's usual `N_WL_DENSE` interpolation
grid (see `ring_resonator_hybrid_cascaded.py`'s own grid-convergence note).

**This FSR estimate uses only the FUNDAMENTAL mode's n_g** — and per
Step 0's follow-up measurement (<0.01% of excited power reaches the
2nd-order-mode pair), that is now a well-justified simplification rather
than a known gap: the real cavity's resonances should be dominated by the
single fundamental-mode family, matching the single-mode Airy-function
picture Step 5 plans to use as a closed-form check. Not a proof at the
full-cavity level (the 0.0085% figure was measured at one bend, one
excitation, one wavelength) — worth a similar spot-check once Step 4
assembles the real network, but no longer treated as an open, possibly-
invalidating gap the way it was before the measurement above.

**Corrected 2026-08-08 — this estimate was off by ~2x, MEASURED, not
argued.** Step 5's real, fully-resolved dense spectrum (22 real resonances
across the C-band) gives a **measured mean FSR of 4.500 nm** — about half
the 9.5 nm estimated above. Two real effects, both missing from the rough
estimate above, account for it:

1. **A Fabry-Pérot cavity's FSR uses the ROUND-TRIP length, not the
   one-way length** — `FSR ≈ λ²/(n_g · L_roundtrip)` with
   `L_roundtrip = 2 · L_cavity = 126 µm`, not `L_cavity = 63 µm` as used
   above. That alone predicts **FSR ≈ 4.75 nm** — already within ~5% of
   the real measured 4.5 nm. This was a plain arithmetic omission in the
   estimate above (an ARGUED quantity, exactly the kind of thing this
   doc's own convention expects a later real measurement to catch),
   not a new physical finding.
2. **The remaining ~5% gap** is consistent with the rough estimate's
   other simplification: it applied the BARE STRAIGHT WAVEGUIDE's own
   n_g uniformly across all 63 µm, including the 23 µm spent inside the
   DC's own coupled-pair cross-section (narrower gap, different
   dispersion) — Step 5's real dense evaluation instead uses the DC's
   own FDTD-calibrated `beta(λ)` for that stretch, which is not required
   to share the bare strip's exact group index.

Net effect: **twice as many resonances (22, not ~10-11) across the same
C-band**, all real and independently confirmed via 22 separately-measured
peak-to-peak spacings, not a single formula. See Step 5 for the full
resolved table.

### Step 2 — re-derive whichever existing example moves platforms — ✅ ALL ITEMS DONE — item 3 (real coupling) 2026-08-07, item 2 (resolution convergence, = check 3) 2026-08-08

Under Step 0's decision: redesign `directional_coupler_hybrid.py`'s device
for `SOI_160`, at `dc.gds`'s own real dimensions, **treated as single-mode
(`N_MODES=2`)** per Step 0's follow-up measurement (<0.01% of excited
power reaches the 2nd-order-mode pair). Step 0's mode-count check replaces
item 1 below (no single-mode width to find — the real cross-section is
4-supermode, but the excited physics is measured single-mode):

1. ~~Find a single-mode strip width at 160 nm core thickness~~ — **N/A,
   superseded by Step 0**: the real device's own width (0.5 µm) supports 2
   guided modes per arm, but the S-bend's own measured mode conversion
   into the 2nd-order pair is negligible (0.0085%) — `N_MODES=2` is used
   as a measured simplification, not a redesigned single-mode
   cross-section.
2. ✅ **`FDTD_RESOLUTION`/`EME_RESOLUTION` convergence study on `SOI_160`
   — MEASURED 2026-08-08 (= check 3), not just assumed.** Originally
   deferred per the user's explicit instruction to first assume the
   `SOI_220` script's own numbers (40/80) — later picked back up (see
   Step 6, check 3). `RES_SCALE ∈ {0.75, 1.25}` (`FDTD_RES/EME_RES =
   30/60` and `50/100`) bracketing the production `40/80`: `res=40→50`
   moves the coupling ratio by 0.00%, `res=30→40` moves it +7.35% — real,
   direct confirmation the production resolution sits at or past the
   convergence knee, not a silently-wrong platform-specific guess the way
   `directional_coupler_validation.md`'s own item 11 warned this class of
   assumption can be.
3. **Real coupling measurement.** Script:
   `examples/directional_coupler/python/directional_coupler_soi160_hybrid.py`
   — 2 FDTD arm excitations (`in`/`isolated`) + 1 EME body mode-solve at
   the GDS's fixed 3.0 µm `coupling_length` (not a 50:50 search — see
   Step 1), the same 2×2 transfer-matrix shape the `SOI_220` script's own
   core physics uses, at `res=40`/`80` (item 2's assumption). A login-node
   plumbing smoke test (`res=8`, `decay_tol=1e-3`) confirmed the real code
   path runs correctly before failing an expected accuracy gate at that
   deliberately crude resolution.

   **Timed out twice before landing — a real, measured cost finding, not
   an inconvenience to smooth over:**

   - **Job 11691379** (`-t 300`, 5 h) TIMEOUT. Diagnosed from the job's
     own log + `example_sim_cache/` (not guessed): both FDTD arm
     excitations (`M_in0`, `M_in1`) and all 3 of `section_mode_beta()`'s
     own internal calibration-length sub-solves (lengths 0.5/1.0/1.5 µm)
     completed and cached to disk cleanly before the timeout (5 distinct
     cache entries, valid `config.json`+`result.pkl`, exact
     `geometry_polygons` x-spans confirming which piece each one is; no
     leftover `.tmp` files).
   - **Job 11706987** (`-t 240`, 4 h, same shared cache) ALSO TIMEOUT.
     All 5 expensive pieces from job 11691379 correctly cache-hit this
     time (confirmed: no matching "Initializing structure" domain or new
     cache entry for any of them appears in this run's own log) — but the
     full 4 h budget was still consumed, apparently inside two calls to
     `section_mode_neff()` (one internal to `section_mode_beta()`'s own
     seed step, one from this script's own MPB-comparison print), a
     method whose own docstring calls it "cheap... seconds, not minutes."
     **The real cost driver inside this step is not yet fully
     diagnosed** — flagged honestly as an open gap, not resolved by
     further guessing.
   - **Shown both timeouts, the user's explicit call (2026-08-06) was to
     resubmit with a much larger time margin rather than keep
     investigating the exact cost driver first.** Job **11716260**
     (`batch/dc_soi160_hybrid.batch`, `-t 720`, 12 h, same shared cache)
     — not a measured bound, a deliberately large margin given two prior
     estimates (5 h, 4 h) both proved insufficient.
   - **Job 11716260 COMPLETED, 2026-08-07 — 4h56m real walltime end to
     end** (17774 s). So the true cost of this measurement, uninterrupted,
     is just under 5 h — close to job 11691379's original `-t 300`
     estimate, meaning that first timeout was a near-miss, not wildly
     under-budgeted; the second attempt's own timeout (with 5 of 6
     expensive pieces cache-hit) is still the genuinely unexplained
     result. **Real result, at the fixed `coupling_length=3.0 µm`, λ=1.55
     µm**: `T_through=0.8989` (IL −0.46 dB), `T_cross=0.0979` (IL
     −10.09 dB), `T_total=0.9969`, coupling ratio 0.098 (≈10% cross,
     90% through) — passive (`T_total<1`, small real loss, not a red
     flag). Body calibration (`CALIBRATE_BODY=True`) ran successfully:
     MPB's own Δn_eff (0.028758) vs. the FDTD-measured Δn_eff (0.025405)
     differ by ~12%, the same class of MPB/FDTD phase disagreement
     `directional_coupler_hybrid.py`'s own docstring documents for its
     `SOI_220` device — confirms this calibration step is doing real,
     necessary work here too, not a formality. Full complex `t(λ)`/`k(λ)`
     at 11 wavelengths saved to
     `example_output/directional_coupler/directional_coupler_soi160_hybrid/smatrix_hybrid.pkl`.

Cheaper than originally scoped (Step 0's own quick check, `res=25` on the
login node, ran in ~24 min for one excitation — the full convergence
sweep is still real SLURM-scale work across several resolutions, but at
`N_MODES=2` rather than 4) — should still be scoped as its own SLURM-batch
effort, not folded silently into Step 4's assembly work. **These two
timeouts are themselves evidence this whole cost assumption needs
revisiting if Step 2 is ever repeated at other resolutions** — a single
fixed-length coupling measurement at `res=40` has now needed close to
9 hours of wall-clock across two attempts and still hasn't landed; a
genuine resolution CONVERGENCE study (item 2, not attempted this pass)
would cost several times that, and the specific cost driver behind the
second timeout is worth understanding properly before committing to one.
several times that.

### Step 3 — the two EME waveguide phase blocks — ✅ REAL BUG FOUND+FIXED, REAL DATA COMPUTED 2026-08-05

Per the user's own count (confirmed against A2): **one** `EMESolver`
mode-solve for the waveguide connecting each HCG to the DC, reused for
both the left and right arms under the Step 1 symmetric-cavity assumption
(mirroring how the HCG's own convergence sweep and the DC's own S-bend
both get reused left/right by geometric symmetry rather than re-solved).

**Treated as single-mode, per Step 0's decision.** The bare strip itself
carries 2 real guided modes (Step 0's own check), but this is a plain
STRAIGHT section — no bend, no discontinuity, nothing to convert power
between modes — and both of its neighbors only ever carry mode-1 content
in this device's assembly: the HCG mirror's own `solve_reflection()` only
tracks `input_mode=1`, and the DC is itself now treated at `N_MODES=2`
(Step 2, the fundamental-mode supermode pair only). So the connecting
waveguide only ever needs to propagate the fundamental mode's phase — a
genuinely single-entry `exp(2πi·n_eff,1/λ·L)` phase block, not a
multimode diagonal.

**This plan originally said to build that block via
`beef.network.convert.from_transfer_matrix()` — WRONG, caught with
synthetic data before any real device numbers were involved
(`tests/test_fabry_perot_network_crosscheck.py`).**
`from_transfer_matrix()` is documented as one-directional by design (its
own docstring: *"Self-terms — reflection at either port, or reverse
transmission — are zero-filled, matching what those solvers actually
model: pure forward transmission, no reflection"*) — correct for
`HybridSolver`'s own straight-through cascades, where nothing downstream
ever reflects back through an earlier section. **This cavity's connecting
waveguides must carry light BOTH forward (toward the HCG mirror) and
backward (the mirror's own reflection, through the SAME physical
waveguide)** — using `from_transfer_matrix()` here builds a waveguide
block that cannot carry anything back through itself, and `connect()`ing
a mirror to it reduces to EXACTLY ZERO, not a small error (verified
directly: a synthetic `connect(mirror_block, "p", from_transfer_matrix(...)
, "in")` gives `0+0j`, when the physical answer is `r*p` — one mirror
reflection seen through one waveguide pass).

**Fix**: build the waveguide's `SMatrixBlock` directly, populating BOTH
directions with the same phase (`S[out,in] = S[in,out] = exp(2πi·n_eff,1/λ·L)`,
zero self-reflection at either end) — ordinary EM reciprocity for a
uniform, lossless section, not a new physics assumption. Verified against
a full `Circuit` assembly (5 blocks, the exact wiring Step 4 needs) AND an
independent brute-force 4-port linear solve, both to floating-point
precision — see Step 5.

**Now solved for real, 2026-08-05**: fundamental-mode `n_eff(λ)` for the
real 0.5 µm/160 nm bare strip, at all 21 wavelengths of the HCG's own
production grid (1.50–1.60 µm) — a real `meow`/MPB mode-solve (`res=100`,
num_modes=1, ~8.8 min total on the login node, no FDTD needed since this
is a straight, uniform section). `n_eff` runs 2.2113 (1.50 µm) down to
2.0912 (1.60 µm), monotonic across the band.

**The real HCG mirror block (A3) and this real waveguide block have been
connected** (`connect(hcg_block, "p", wg_block, "out")`, the exact
Step 4 wiring pattern, using REAL numbers on both sides) to produce
`rho(λ) = r(λ)·p(λ)²` — the HCG mirror as seen from the DC's side of one
20 µm arm, at all 21 wavelengths. `|rho|²` reproduces the HCG's own
measured `R(λ)` spectrum exactly (as it must — propagation phase doesn't
touch magnitude), peaking at 0.9714 at 1.580 µm; `arg(rho)` wraps rapidly
across the band (the 20 µm path length accumulates ~170 rad of phase at
these wavelengths, many multiples of 2π between adjacent 5 nm-spaced grid
points — expected for a path this many wavelengths long, not a bug).
**This `rho(λ)` is the actual input Step 5's closed form needs** — only
Step 2's real `t(λ)`/`k(λ)` remain before the full cavity spectrum can be
evaluated with entirely real numbers.

### Step 4 — assemble with `beef.network.Circuit` — ✅ DONE WITH REAL DATA 2026-08-07

```python
from beef.network import Circuit

circuit = Circuit()
circuit.add_block("hcg_left", hcg_block)          # A3's direct 1-port construction
circuit.add_block("wg_left", wg_block)             # Step 3's FIXED bidirectional construction
circuit.add_block("dc", dc_block)                  # Step 2's re-solved DC, as a 4-port block
circuit.add_block("wg_right", wg_block)            # same block reused, per Step 3
circuit.add_block("hcg_right", hcg_block)          # same block reused, per HCG's own left-right symmetry

circuit.connect("hcg_left", "p", "wg_left", "in")
circuit.connect("wg_left", "out", "dc", "through_left")
circuit.connect("dc", "through_right", "wg_right", "in")
circuit.connect("wg_right", "out", "hcg_right", "p")

cavity = circuit.reduce(external_ports={
    "dc.input_reflection": "input_reflection",
    "dc.output": "output",
})
```

**This exact 5-block, 4-`connect()` wiring was first verified**
(`tests/test_fabry_perot_network_crosscheck.py::
test_circuit_matches_closed_form_symmetric_cavity`) against the
independent closed form derived in Step 5, using synthetic mirror/DC/
waveguide blocks — the first real exercise of `Circuit.add_block()`
composing more than one independently-defined physical device (A2) and of
a 1-port self-reflecting block (A3), both previously unexercised anywhere
in this repo.

**Now run with entirely real device data, 2026-08-07 — new permanent
script**: `examples/high_contrast_grating_fabry_perot/python/
hcg_fabry_perot_hybrid.py`. Loads the HCG's real production `SMatrix`
and the DC's real `smatrix_hybrid.pkl` directly (does not re-solve
either — both are their own expensive SLURM measurements, already
landed), runs a fresh real fundamental-mode-only `meow`/MPB sweep for the
connecting waveguide (~5 min, no FDTD), and assembles all five real
blocks through the literal `Circuit` wiring above. Runs end-to-end on a
login node in ~5 minutes (dominated entirely by that one mode-solve;
everything else is either loading existing files or instant network
algebra).

**Updated 2026-08-08 — the output folder is now a self-contained picture
of the whole device, not just the combined spectrum.** The script now
copies the HCG mirror's own `geometry.png`/`reflection.png` (into
`hcg_mirror/`) and the DC's own `geometry_hybrid.png`/`hybrid_results.png`/
`hybrid_results_dense.png`/`sections/*.png` (into `directional_coupler/`)
from each constituent device's own real output — copies, not moves, so
each device's own folder stays the authoritative source. It also builds a
new `overall_geometry.png`: a labeled schematic of the assembled cavity to
real x-scale (HCG mirror 5.389 µm — waveguide 53.08 µm — DC 23.00 µm —
waveguide 53.08 µm — HCG mirror, total 139.94 µm, with the DC's tap ports
marked) — explicitly NOT a shared FDTD domain (each piece is
independently solved in its own coordinate frame, A2), a visualization
aid, not a new simulation.

**The DC's own 4-port block port-mapping question is resolved, and
checked, not just asserted**: `directional_coupler_soi160_hybrid.py`
only ever computes and stores `through`/`cross` given excitation from
`in` — the `isolated`-excited branch (`k'`, `t'`) was never separately
saved. Rather than assume the standard symmetric-coupler equality
(`t'==t`, `k'==k`), it was checked directly against this device's own
cached `M_in0`/`M_in1`/`M_body` pieces (a cache-hit reload, no new FDTD):
**`k'` matched `k` to machine precision, `t'` matched `t` to ~1e-5** — so
building the 4-port block from the saved `t`/`k` alone (used twice, for
both the `L↔R` and `I↔O` pairs) is a verified equality for this real
device, not an unchecked assumption. Port mapping used: `in`→
`through_left`, its mirrored/output-side partner→`through_right`,
`isolated`→`input_reflection`, its partner→`output`.

**Result: the real `Circuit` assembly matches Step 5's real closed-form
spectrum to floating-point precision** (`max|ΔR|=6.9e-17`,
`max|ΔT|=1.1e-15` at Step 1's real `L_ARM=53.08 µm`, re-confirmed after
that length was corrected from its original 20 µm placeholder — asserted
inside the script itself, not just checked by hand) — two independent
code paths (general network composition vs. a hand-derived scalar
formula), same real device data, identical answer, unaffected by the
arm-length correction.
This is the strongest evidence yet that both the wiring and the physics
are self-consistent, though — same caveat as Step 5 — the underlying
10 nm-spaced DC grid is still far too coarse to resolve the resonance
this spectrum itself reveals near 1.57 µm.

### Step 5 — closed-form sanity check — ✅ FULLY RESOLVED WITH REAL NUMBERS 2026-08-08 — 44 real resonances at the real cavity length, real FWHM/Q/depth per resonance

Every ring-resonator example in this repo has a closed-form check to catch
a wiring/sign/port-order mistake in the network assembly before trusting
it: the ring's own `t + k1·k2·P/(1 - t_ring·P)` all-pass formula,
independently cross-checked against `close_loop()` in
`tests/test_network_ring_crosscheck.py`. This device's equivalent, for a
symmetric cavity (both mirrors reflection `r`, both connecting waveguides
single-pass phase `p`, DC through/cross amplitudes `t`/`k`, one
wavelength):

```
rho = r * p**2                                # mirror + ROUND-TRIP waveguide phase
D = 1 - rho**2 * t**2
self_reflection = k**2 * rho / D              # tap port -> itself (input_reflection or output)
transmission    = t * (1 + rho**2*(k**2 - t**2)) / D   # input_reflection <-> output
```

Hand-derived from the DC's 4-port S-matrix (`L,R` = the two cavity-rail
ports, `I,O` = the two tap ports) plus the two mirror boundary conditions
(`a_L = rho·b_L`, `a_R = rho·b_R`), then **verified two independent ways**
(`tests/test_fabry_perot_network_crosscheck.py`):

1. Against a brute-force direct linear solve of the same 4-port system
   (200 random trials, `atol=1e-9`) — catches an error shared between the
   closed form and any Step-4-style reduction, which a single cross-check
   cannot.
2. Against the actual `Circuit` assembly (Step 4's exact 5-block wiring),
   to `atol=1e-10`, including the reciprocity check
   `S[output, input_reflection] == S[input_reflection, output]`.

Also confirmed passive (`|self_reflection|² + |transmission|² ≤ 1`) for a
physically valid, lossy-mirror case using a genuinely unitary DC
(`t=cos θ`, `k=i·sin θ` — the standard symmetric-coupler quadrature form;
an earlier attempt with arbitrary, non-unitary `t`/`k` phases spuriously
"violated" passivity, which was the test's own flaw, not the formula's —
see that test's docstring for the derived 90°-phase-relationship
condition this device's real, physically-lossless DC will satisfy
automatically).

**Now evaluated with entirely real numbers — MEASURED, 2026-08-07.** Job
11716260 (Step 2) finally landed the DC's real `t(λ)`/`k(λ)` (11 pts,
1.50–1.60 µm, step 0.01 µm — coarser than the HCG's own 21-pt/0.005 µm
grid, but an exact subset of it). Combined with the HCG's real `r(λ)` (its
own production `smatrix.pkl`) and Step 3's real waveguide `n_eff(λ)`
(giving `rho(λ) = r(λ)·p(λ)²`, sliced onto the DC's coarser 11-pt grid),
the closed form above gives the first real predicted cavity spectrum:

| λ (µm) | \|rho\|² | R (input_reflection) | T (output) |
|---|---|---|---|
| 1.500 | 0.9284 | 0.0032 | 0.9854 |
| 1.510 | 0.9368 | 0.0011 | 0.9937 |
| 1.520 | 0.9449 | 0.0074 | 0.9776 |
| 1.530 | 0.9524 | 0.0079 | 0.9794 |
| 1.540 | 0.9588 | 0.0023 | 0.9925 |
| 1.550 | 0.9640 | 0.0034 | 0.9908 |
| 1.560 | 0.9679 | 0.0183 | 0.9677 |
| **1.570** | **0.9704** | **0.1710** | **0.7395** |
| 1.580 | 0.9714 | 0.0150 | 0.9758 |
| 1.590 | 0.9708 | 0.0087 | 0.9858 |
| 1.600 | 0.9686 | 0.0081 | 0.9872 |

Passivity holds everywhere (`R+T` max 0.9953, ≤1 as required). **A real
resonance feature shows up at λ=1.570 µm**: transmission drops from its
otherwise-flat ~0.97–0.99 to 0.7395, and reflection spikes from
~0.001–0.02 up to 0.1710 — exactly the interference signature a
Fabry-Pérot cavity is supposed to produce, not noise (both channels move
together, consistently with energy redirected from the output port back
toward the input/reflection port at resonance).

**But this is very likely under-sampled, not a resolved resonance
lineshape — stated honestly, not oversold.** Step 1's own FSR estimate
was ~9.5 nm; this grid's point spacing is 10 nm — comparable to, not much
finer than, the estimated resonance spacing. A single grid point landing
near a dip is consistent with a real, much narrower resonance being
mostly invisible between samples elsewhere in the band (nothing rules out
additional resonances between the other 10 points, and the true peak
depth/width at 1.570 µm itself could differ substantially from this one
sample). **A denser wavelength sweep is needed before this spectrum
supports any quantitative claim about resonance depth, width, Q, or exact
position** — this result establishes that the assembled model DOES
produce genuine resonant behavior at a physically plausible location, not
that the resonance is now characterized.

**Follow-up: linewidth sweep started, 2026-08-07 — sized from real
numbers, not guessed.** `directional_coupler_soi160_hybrid.py` gained
`DC160_WL_MIN`/`DC160_WL_MAX`/`DC160_N_WL` env-var overrides (unset
reproduces job 11716260's original full-band run exactly; any override
redirects output to its own `windows/` subfolder so it can never
overwrite that production result). Sizing, from job 11716260's own real
measured numbers:

```
round-trip power |rho^2 t^2| = |rho|^2 * |t|^2 = 0.9714 * 0.8989 = 0.8732
finesse F = pi*sqrt(0.8732)/(1-0.8732) ~ 23
FWHM ~ FSR/F = 9.5 nm / 23 ~ 0.41 nm   (Step 1's own FSR estimate)
```

— roughly 24× narrower than the original 10 nm grid, not just "a bit
coarse." Submitted job **11726827** (`batch/dc_soi160_hybrid_zoom.batch`,
`-t 720`): window 1.555–1.585 µm (30 nm, ~3 FSRs of context around the
found dip), `N_WL=601` → 0.05 nm spacing, ~8 points across the estimated
FWHM.

**Job 11726827 TIMED OUT, 2026-08-08 — a real, measured cost blowup, much
larger than the `-t 720` margin assumed.** Confirmed via `sacct`
(`TIMEOUT`, `12:00:25` elapsed, ended 2026-08-08T00:41:45) and the job's
own log (`reports/11726827.out`). Root cause, now measured directly
instead of argued: **the narrower/denser wavelength window does not
cache-hit against job 11716260's saved pieces** (the cache key includes
the wavelength grid, and this is a different grid), so this run had to
redo essentially every expensive FDTD sub-solve from scratch — AND each
one now carries a much heavier end-of-run eigenmode decomposition (601
frequencies instead of 11).

**Only 3 pieces fully completed and cached in the full 12 h**: both arm
excitations (`input_idx=0` at 15:17, `input_idx=1` at 17:55 — each ~2.6 h
including decomposition) and ONE calibration-length sub-solve (`L=0.5 µm`,
`reverse=True`, completed 22:49 — ~4.9 h by itself, longer than either arm
excitation, likely because it decomposes onto 2 modes × 601 frequencies
instead of 1 × 601). The job's own log never printed `section_mode_beta`'s
own "Body beta calibration at lambda=..." line, confirming calibration
itself never finished. **The log shows at least 7 distinct
`Initializing structure...` blocks (2 arm + up to 5 more)** — more than
the 3 calibration lengths alone would need at one sub-solve each, hinting
`section_mode_beta()` may run each calibration length in BOTH directions
(forward/reverse) rather than once — a real structural detail this doc
had not previously pinned down precisely (the earlier "5 expensive FDTD
pieces" count may itself have been an undercount, just never stressed
enough by N_WL to show it before). **Not yet resubmitted — this is a real
allocation/design decision (bigger time budget vs. fewer wavelength
points vs. a narrower window) put to the user rather than assumed.**

**Superseded 2026-08-08 — solved with interpolation instead, no
resubmission needed.** The brute-force approach (resubmitting
`directional_coupler_soi160_hybrid.py` at higher `N_WL`) was the wrong
shape for this problem: it re-derives every expensive FDTD sub-solve from
scratch just to get more monitored frequencies, when the underlying
physics doesn't need that. **Per-factor phase analysis of the REAL cached
11-point data** (not assumed):

| quantity | phase step per 10 nm | `unwrap_is_safe(margin=0.5)` |
|---|---|---|
| `M_in0`/`M_in1` (arm→body-mode overlap) | ~65-69° | ✅ safe |
| `beta` (FDTD-calibrated body propagation const.) | ~15-16° | ✅ safe (degree-2 poly fit, residual 0.02% of range) |
| `coeff_through`/`coeff_cross` (the FINAL combined coefficient) | ~153° | ❌ unsafe — this is what made the original 11-point grid look "too aliased to interpolate" |

The final coefficient's large phase step is not new physics needing dense
resampling — it's this device's own mirror-symmetry reuse
(`M_out_through = M_in0ᵀ`, Step 2/A2) counting the S-bend's phase TWICE.
Each individual factor is comfortably safe to interpolate; only their
combination looked dangerous. This matches standard practice, not a
special-case trick: Lumerical's own EME documentation describes exactly
this as a "perturbative method... appropriate for applications where the
wavelength dependence is mainly due to phase change rather than mode
profile change," and coupled-mode theory describes directional-coupler
coupling the same way — governed by `(β_even − β_odd)·L`, a smooth
propagation-constant difference times a known, fixed length, exactly
`M_body`'s own measured structure (diagonal, `|M_body|=1`, pure phase).

**Real validation, not just an argument**: job 11726827's timeout wasn't
wasted — 2 of its expensive pieces (both 601-point arm excitations) fully
completed and cached before the timeout, at the exact same dense window.
Interpolating the sparse 11-point `M_in0`/`M_in1` (magnitude + unwrapped
phase, `beef.spectral.interp_complex`) reproduces that REAL 601-point
measured data to **≤1.3% relative error at every point** — genuine ground
truth this scheme was not designed with, not a self-consistency check.

**New permanent scripts**: `directional_coupler_soi160_hybrid_dense.py`
(loads the sparse run's own `coupling_pieces.npz`, interpolates the safe
factors, fits `beta(λ)` with a degree-2 polynomial matching
`ring_resonator_hybrid.py`'s own `neff_fit` convention, recombines via the
same contraction the sparse script uses, saves `smatrix_hybrid_dense.pkl`)
and a new "section 8" appended to `hcg_fabry_perot_hybrid.py` (interpolates
the HCG's own `r(λ)` — already safe, ~7.5°/5nm — fits the connecting
waveguide's `neff(λ)` by reusing the ALREADY-COMPUTED 11 points above, no
new MPB solves, and evaluates the closed form at the dense grid; does NOT
re-run `Circuit` at thousands of points, since Circuit and the closed form
already agree to floating-point precision on the sparse grid, section 6).

**The real, fully resolved spectrum — 20001 points, 1.50–1.60 µm, 0.005 nm
spacing — MEASURED 2026-08-08**:

| λ (µm) | T_min | FWHM (nm) | Q |
|---|---|---|---|
| 1.5035 | 0.2926 | 0.185 | 8127 |
| 1.5077 | 0.2666 | 0.185 | 8150 |
| 1.5120 | 0.2412 | 0.185 | 8173 |
| 1.5163 | 0.2168 | 0.185 | 8196 |
| 1.5206 | 0.1932 | 0.190 | 8003 |
| 1.5250 | 0.1701 | 0.190 | 8026 |
| 1.5293 | 0.1492 | 0.190 | 8049 |
| 1.5337 | 0.1301 | 0.195 | 7865 |
| 1.5381 | 0.1131 | 0.195 | 7888 |
| 1.5426 | 0.0981 | 0.200 | 7713 |
| 1.5471 | 0.0845 | 0.205 | 7547 |
| 1.5516 | 0.0729 | 0.210 | 7388 |
| 1.5561 | 0.0626 | 0.215 | 7238 |
| 1.5606 | 0.0542 | 0.225 | 6936 |
| 1.5652 | 0.0472 | 0.230 | 6805 |
| **1.5698** | **0.0415** | **0.235** | **6680** |
| 1.5744 | 0.0369 | 0.250 | 6298 |
| 1.5791 | 0.0335 | 0.260 | 6073 |
| 1.5838 | 0.0307 | 0.270 | 5866 |
| 1.5885 | 0.0284 | 0.285 | 5574 |
| 1.5932 | 0.0268 | 0.300 | 5311 |
| 1.5980 | 0.0259 | 0.310 | 5155 |

**22 resonances, mean FSR 4.500 nm** (see Step 1's correction above —
roughly half the original ~9.5 nm estimate, mostly a missing round-trip
factor of 2). Passivity holds everywhere (`R+T` max 0.9955). Depth and Q
both trend smoothly across the band (T_min 0.29→0.026, Q 8200→5150) as
the HCG mirror's own reflectivity and the DC's coupling ratio each drift
with wavelength — a smooth, physically expected trend, not a numerical
artifact. **The resonance nearest the original 11-point grid's own sample
(λ≈1.570 µm) is real, but was badly under-sampled exactly as suspected**:
the coarse grid read `T=0.7395` at that point; the actual resolved dip is
**T_min=0.0415 at λ=1.5698 µm, FWHM=0.235 nm, Q≈6680** — about 18× deeper
and far narrower than the single sparse sample suggested. Full 22-point
table, `smatrix_fabry_perot_dense.pkl`, and `spectrum_dense.png` (a clean,
evenly-spaced comb, visually confirming the reconstruction is physical)
saved under `example_output/high_contrast_grating_fabry_perot/
hcg_fabry_perot_hybrid/`.

**SUPERSEDED 2026-08-08 — Step 1's real arm length (53.08 µm, not the
20 µm placeholder above) changes every number in this step.** Everything
above this note (the sparse 11-point table, the 22-resonance dense table,
the specific 1.570 µm anecdote) was computed at the OLD provisional
20 µm/arm length — kept in place as a real record of what was found at
that length, not deleted, but no longer the current result. Re-run with
`L_ARM=53.08 µm` (`hcg_fabry_perot_hybrid.py`'s own `L_ARM` constant,
derived in Step 1 from `duty/Mux_lpbk_9$4.gds`'s real layout): the
interpolation scheme, the per-factor safety checks, and the ≤1.3% real
validation are all UNCHANGED (none of that depends on `L_ARM`) — only the
connecting-waveguide phase term changes, which is exactly what a longer
real cavity should do: shrink the FSR and raise the Q, not invalidate the
method.

**Circuit vs. closed-form cross-check, at the real length**:
`max|ΔR|=6.9e-17`, `max|ΔT|=1.1e-15` — still floating-point precision.
**Sparse (11-point) spectrum, at the real length**:

| λ (µm) | R (input_reflection) | T (output) |
|---|---|---|
| 1.500 | 0.0072 | 0.9713 |
| 1.510 | 0.0042 | 0.9844 |
| 1.520 | 0.0014 | 0.9932 |
| 1.530 | 0.0379 | 0.9122 |
| 1.540 | 0.0023 | 0.9924 |
| 1.550 | 0.0862 | 0.8444 |
| 1.560 | 0.0034 | 0.9919 |
| 1.570 | 0.0112 | 0.9805 |
| 1.580 | 0.0718 | 0.8937 |
| 1.590 | 0.0091 | 0.9852 |
| 1.600 | 0.0077 | 0.9879 |

Passivity holds (`R+T` max 0.9955). With the longer, real cavity the
under-sampling is worse, not better — the sparse grid now lands its
deepest reads at 1.550/1.580 µm rather than 1.570 µm, itself more evidence
that a 10 nm grid cannot be trusted to even locate a resonance, let alone
characterize it.

**Dense (20001-point) resolved spectrum, at the real length — MEASURED
2026-08-08: 44 real resonances, mean FSR 2.255 nm** (roughly half of the
20 µm case's 4.5 nm — consistent with the real cavity being ~2.05× longer,
`L_roundtrip` scaling from 126 to 258.32 µm). Representative rows (full
44-point table in the script's own stdout / `smatrix_fabry_perot_dense.pkl`):

| λ (µm) | T_min | FWHM (nm) | Q |
|---|---|---|---|
| 1.5019 | 0.3029 | 0.095 | 15809 |
| 1.5104 | 0.2513 | 0.095 | 15899 |
| 1.5233 | 0.1800 | 0.095 | 16035 |
| 1.5365 | 0.1211 | 0.105 | 14633 |
| 1.5499 | 0.0773 | 0.105 | 14761 |
| 1.5613 | 0.0533 | 0.115 | 13576 |
| **1.5705** | **0.0406** | **0.120** | **13087** |
| 1.5798 | 0.0340 | 0.135 | 11702 |
| 1.5893 | 0.0287 | 0.145 | 10961 |
| 1.5988 | 0.0261 | 0.160 | 9993 |

Passivity holds everywhere (`R+T` max 0.9955). Depth and FWHM both trend
smoothly across the band (T_min 0.30→0.026, FWHM 0.095→0.16 nm, Q
16000→9990) — the same physically-expected drift as before, at roughly
double the Q (a longer cavity, same round-trip loss fraction, is a higher-
finesse-per-unit-bandwidth resonator, exactly as expected). This is now
the current, real, fully-resolved spectrum for this device.

### Step 6 — the four checks

Convert to a normal gap list against
[the four checks](../README.md#the-four-checks) once Steps 0-5 land.

**Check 2 (mode identity/ordering) — ✅ CLOSED 2026-08-08, MEASURED.**
Never run for this device before now. Extended
`examples/verification/python/audit_interface_modes.py` (the same
"cheap by construction" tool already used for the SOI_220 coupler and
both MMIs — `FDTDSolver.section_mode_neff()` runs meep's own eigenmode
solver with no time-stepping, minutes not hours) with a new
`directional_coupler_soi160` config, pinned to this device's own real
geometry (`waveguide_width=0.5, gap=0.25, coupling_length=3.0,
bend_length=10.0, waveguide_offset=4.0`, `SOI_160`, `N_MODES=2`), audited
at a short cut safely inside the real 2.5 µm body (`section=(-0.25,0.25)`,
padding to ±0.75 after `section_mode_neff()`'s own `dpml=0.5` clearance —
0.5 µm of real margin from the actual ±1.25 µm body edge). **Result**:
both engines return 2/2 guided modes, **identity correspondence
(`[0, 1]`), order preserved, worst guided `|dn_eff|=0.00048`** — meep's
and meow's own eigenmode solvers agree to ~0.02%, tighter than the
~0.57% this repo's own check-1 discussion cites as typical MPB-vs-MPB
agreement. No reordering risk found — this device's `N_MODES=2` basis is
simple enough (vs. e.g. the modal multiplexer's own 16-mode, 4-interface
case, where this exact check caught a real bug) that a clean result here
is not surprising, but it is now measured, not assumed. Saved to
`example_output/verification/interface_mode_audit/
directional_coupler_soi160.json`.

**Check 1 (engine agreement) — ✅ CLOSED 2026-08-08.** `CALIBRATE_BODY=True`
already measures and compensates for the real FDTD-time-stepped-phase
vs. MPB-mode-index
disagreement (Step 2: MPB Δn_eff=0.028758 vs. FDTD-measured Δn_eff=0.025405
at λ=1.55 µm, ~12% apart, same class of gap the SOI_220 coupler's own
check 1 found and fixed the same way) — but that was only ever checked at
ONE wavelength.

**Extended to all 11 wavelengths, at ZERO new compute** — `section_mode_
beta()` already computes a per-(wavelength, mode) phase-fit residual and
FDTD/MPB "seed excess" internally for every calibration it runs, but only
returns the aggregate `beta` array, discarding the rest after logging just
the single worst cell. Reloaded the SAME cached 3-length calibration
phases (`solve_section()` at lengths `[0.5, 1.0, 1.5] µm`, a pure
cache-hit) and called `beef.solvers.fdtd.beta_fit_from_phases()` directly
— the same public function the library itself uses — to recover the full
grid (~55 s total, no new FDTD, no new EME beyond one cheap MPB reseed):

| quantity | result |
|---|---|
| `Δβ_FDTD/Δβ_MPB` ratio, all 11 λ | **0.8801 ± 0.00084** (0.28% relative variation across the full 100 nm band) |
| Worst phase-fit residual, all 22 cells | `9.3×10⁻⁴ rad` (`_BETA_FIT_RESIDUAL_MAX=0.1` — >100× under threshold) |
| Worst FDTD/MPB seed excess, all 22 cells | `2.72%` (same class as the SOI_220 coupler's own +2.23%/+2.36%) |

The ~12% `Δβ` gap `CALIBRATE_BODY` compensates for is not a one-wavelength
fluke — it is a smooth, remarkably stable (<0.3% relative drift) function
across the entire operating band, and every single calibration phase-fit
is a genuinely clean straight line (residuals two orders of magnitude
under the guard threshold). That is real, if partial, independent
evidence the compensation is measuring something physically real, not an
artifact.

**What this does NOT fully replicate**: the SOI_220 coupler's own STRONGER
check-1 confirmation — a stub/cut-point sweep testing whether a
PHYSICALLY MEANINGLESS modeling choice (where the FDTD/EME boundary sits)
leaves the predicted L₅₀ invariant, with the observed drift slope matching
the value `2·(1−Δβ_FDTD/Δβ_EME)` predicts from the direct measurement
(agreement to 0.82% there). This device has no L₅₀ (`coupling_length` is
fixed, not swept for a 50:50 point) — the adapted analogue tested is
whether the FULL hybrid pipeline's own predicted coupling ANGLE
(`θ=asin(√T_cross)`) at the FIXED `coupling_length=3.0 µm` drifts with
`BODY_STUB` by the amount `dθ/d(stub) = Δβ_EME·(ratio−1)` predicts from
the already-measured `ratio=0.8801`, `Δβ_EME=0.11704` rad/µm at λ=1.55 µm.

**Job 11747764 landed, 2026-08-08 — real result, but inconclusive at this
sensitivity, reported honestly rather than oversold.**
`DC160_BODY_STUB ∈ {0.1, 0.5, 0.75}` (single λ=1.55 µm point each,
`-t 720`, actual cost 1h40m/1h40m/2h25m) plus the already-real `0.25`
production point:

| `BODY_STUB` (µm) | `T_through` | `T_cross` | θ (rad) |
|---|---|---|---|
| 0.10 | 0.8995 | 0.0973 | 0.31722 |
| 0.25 (production) | 0.8989 | 0.0979 | 0.31823 |
| 0.50 | 0.8990 | 0.0980 | 0.31840 |
| 0.75 | 0.9000 | 0.0971 | 0.31689 |

`T_cross` (and θ) is stable to <1% across the whole tested range — a real,
reassuring result for the production `BODY_STUB=0.25` choice, whichever
way the rest of this reads. Fitting a line to θ vs. stub: **measured
`dθ/d(stub) = −0.00062` rad/µm**, with point-to-point residuals
(`≈5-8×10⁻⁴` rad) of the same order as the fitted slope's own signal over
this range — i.e. the fit is not well-resolved above scatter. That
measured slope is **~23× smaller than the `−0.01403` rad/µm the adapted
formula predicts** (predicted total drift over this stub range: `−0.0091`
rad; measured endpoint-to-endpoint: `−0.00034` rad).

**RESOLVED 2026-08-08 — the mismatch was the adapted formula, not this
device, and the real explanation is a POSITIVE result, not an open
question.** The `2·(1−Δβ_FDTD/Δβ_EME)`-style formula was derived for the
SOI_220 coupler's own PRE-calibration stub sweep — measured before
`calibrate_body_neff` existed, specifically to EXPOSE the Δβ gap by
propagating the EME body on MPB's own uncalibrated rate. This device's
production run already has `CALIBRATE_BODY=True` from the start, so
applying that same formula compared an uncalibrated prediction against a
calibrated measurement.

**Checked directly, branch-safely, with data already in hand (zero new
compute)**: reloaded each stub's own cached `M_in0` (the real arm-section
FDTD solve) and fit its phase-vs-stub rate with the library's own
`beta_fit_from_phases()` (seeded on the isolated calibration's own β, to
avoid exactly the wrap-branch ambiguity a naive `np.unwrap` across
coarse, unevenly-spaced stub steps risks — the same trap that function
exists to avoid, per its own docstring):

| mode | arm's real FDTD phase rate (branch-safe fit) | isolated calibration rate | agreement | fit residual | margin |
|---|---|---|---|---|---|
| 0 | 8.8720 rad/µm | 8.8437 rad/µm | 0.32% | 0.026 rad | 0.002 (confident) |
| 1 | 8.7696 rad/µm | 8.7408 rad/µm | 0.33% | 0.026 rad | 0.002 (confident) |

Both fits are clean and unambiguous (`turns_shifted=False` — the seed's
own naive branch guess was already right). **This is the mechanism**:
`coeff_through = Σ_m c_m²·e^{i2φ_bend,m}·e^{iβ_m·(2·stub)}·
e^{iβ_body,m·(coupling_length−2·stub)}`. With `CALIBRATE_BODY=True`,
`β_body,m` IS the FDTD-measured rate — the *same* rate the real arm
sections independently arrive at (confirmed above, to 0.3%) — so the
`stub`-dependent terms cancel almost exactly, BY CONSTRUCTION, regardless
of where the cut sits. **The ~23× smaller-than-naively-predicted drift is
not a gap — it is direct, positive confirmation that `CALIBRATE_BODY` is
doing exactly what it is for**: making the hybrid construction
self-consistent and robust to the arbitrary FDTD/EME cut-point choice.
The small residual drift that WAS measured is fully explained by the
0.3%-per-mode gap between the arm's own real accumulation rate and the
isolated calibration's rate (a genuine, small, second-order effect,
plausibly near-field/mode-transition physics right at the S-bend-to-
straight seam) — nothing left unexplained. **Check 1 is now considered
CLOSED**, on the same footing as checks 2 and 3.

**Check 3 (resolution convergence, = Step 2 item 2) — ✅ MEASURED
2026-08-08, converged at the production resolution.**
`DC160_RES_SCALE ∈ {0.75, 1.25}` (`FDTD_RES/EME_RES = 30/60` and `50/100`,
actual cost 36m40s and 4h55m respectively — cost rises steeply with
resolution, consistent with FDTD's own ~res⁴ 3-D time-stepping scaling),
bracketing the production `scale=1.0` (`40/80`):

| `RES_SCALE` | `FDTD_RES`/`EME_RES` | `T_through` | `T_cross` | coupling ratio |
|---|---|---|---|---|
| 0.75 | 30/60 | 0.9059 | 0.0912 | 0.0915 |
| **1.00 (production)** | **40/80** | **0.8989** | **0.0979** | **0.0982** |
| 1.25 | 50/100 | 0.9009 | 0.0979 | 0.0980 |

**0.75→1.00: `T_cross` moves +7.35%. 1.00→1.25: `T_cross` moves +0.00%**
(`0.0979` at both). This is real, direct evidence the production
resolution (`40/80`) sits AT OR PAST the convergence knee — `res=30` is
clearly under-resolved (the outlier), while `res=40` and `res=50` agree to
4 decimal places on `T_cross`. The underlying FDTD-measured `Δn_eff`
wobbles slightly more (`0.024480 → 0.025405 → 0.025113` across the three
resolutions — a −1.15% step from 40→50 against a +3.78% step from 30→40),
consistent with a small residual calibration-level scatter that does not
propagate through to the device-level observable. **Step 2's own
`FDTD_RESOLUTION=40`/`EME_RESOLUTION=80` assumption, carried over
unverified from the SOI_220 device since 2026-08-06, turns out to be
correct — not by luck, now independently confirmed.**

**Check 4 (full-device FDTD control) — ⛔ DECIDED 2026-08-07: OUT OF SCOPE,
not attempted, by user decision (see A4).** Not "deferred pending a cheap
pilot" — skipped outright, before spending any compute on it. Reasoning
(A4): the ring resonator's own full-device FDTD check wasn't just slow, it
was disqualified after actually running — it reported T>1 in a passive
device (up to +8.4%) and a resonance comb that shifted 18% of an FSR
between resolutions, worsening with refinement — with a stated mechanism
(round-trip phase error amplified ~F× by finesse, and FDTD's staircased
boundaries ARE a round-trip path error) that isn't specific to the ring.
This device's own lower finesse (F≈23 vs. the ring's F≈120) helps, but its
HCG mirrors are gratings — many staircased tooth-edges per round trip,
plausibly more discretization surface per bounce than the ring's smooth
bends — which hurts. Net direction unmeasured, and the ring's own 167 h
number is a floor (a smaller, lower-finesse device) not a ceiling for what
a double-HCG-mirror cavity would cost. Given the likely disqualification
and the likely cost, the call is to not run it at all, rather than spend
SLURM time discovering that firsthand.

**This device's validation instead rests entirely on the single-pass/
algebraic checks that carry no round-trip amplification** — the same
resolution the ring resonator's own doc settled on for the same reason
(see its own "what validates this device instead" section):

| basis | status |
|---|---|
| HCG mirror's own FDTD reflection spectrum | ✅ exists (`high_contrast_grating_validation.md`) |
| DC's own FDTD/EME coupling measurement | ✅ real, landed (Step 2) |
| DC's own mode identity/ordering (check 2) | ✅ MEASURED — identity correspondence, worst `\|dn_eff\|=0.00048` |
| Closed-form Fabry-Pérot-with-DC-tap formula vs. brute-force 4-port solve | ✅ machine precision (Step 5) |
| Closed-form vs. the actual `Circuit` assembly, real device data | ✅ floating-point precision, `max\|ΔR\|=6.9e-17` (Step 4, at the real `L_ARM=53.08 µm`) |
| DC coupler reciprocity (`k'=k`, `t'=t`) | ✅ machine/near-machine precision (Step 4) |
| Passivity (`R+T≤1`) at every evaluated point | ✅ holds everywhere measured (Step 5) |
| Engine agreement, independently cross-validated (check 1) | ✅ CLOSED — full-band ratio stable (0.28%); stub sweep's small residual drift fully explained by `CALIBRATE_BODY` making the arm and body propagate on the same measured rate, confirmed branch-safe to 0.3% per mode |
| Resolution convergence (check 3 / Step 2 item 2) | ✅ MEASURED — converged at production `res=40/80` (`res=40→50` moves `T_cross` 0.00%; `res=30→40` moves it +7.35%) |
| Full-device FDTD reference (check 4) | ⛔ out of scope, by decision — not a gap, a documented choice |

No full-device time-domain reference will exist for this cavity, by
choice, not by oversight.

---

### Step 7 — self-contained rebuild, and the last uncalibrated segment — ✅ LANDED 2026-08-18, job 11995681

Two problems, one rewrite of
`examples/high_contrast_grating_fabry_perot/python/hcg_fabry_perot_hybrid.py`.

**Problem 1 — provenance.** The script consumed
`high_contrast_grating_fdtd/smatrix.pkl` and
`directional_coupler_soi160_hybrid/smatrix_hybrid.pkl` from two other
scripts' output folders, plus a dense pickle from a third. So this cavity's
answer depended on which revision of three other scripts had last written
those files, and the run had to be ordered behind them by SLURM dependency
(step 9 of `section_mode_beta_multi_monitor_validation.md` records exactly
that `--dependency=afterok` reasoning).

**Problem 2 — the arms were the one segment never calibrated.** Step 3 built
them from a `meow`/MPB mode index. Everything they cascade with is
referenced to meep's own β: the HCG's `r(λ)` is a meep mode-monitor
measurement, the coupler's S-bends are meep, and its EME body propagates on
its own FDTD-calibrated β. At 2 × 53.08 µm against the coupler's 23 µm,
**~106 µm of the ~130 µm one-way path was carried by the other engine's
index** — and the closed round trip is exactly where mixing engines is not
allowed. The β-calibration audit could not see this, because it audited
which examples *reach* a calibration, and this one reached none.

Both are fixed: the script now solves all three pieces itself against
`example_sim_cache/`, and calibrates the arms.

**The rebuild is a refactor — measured, not asserted.** The script compares
its own rebuilt coefficients against the published pickles whenever those are
still on disk:

| rebuilt here | vs published | max abs difference |
|---|---|---|
| HCG mirror `r(λ)`, 21 pts | `high_contrast_grating_fdtd/smatrix.pkl` | **0.000e+00** |
| coupler `t(λ)`, 11 pts | `directional_coupler_soi160_hybrid/smatrix_hybrid.pkl` | **1.848e-15** |
| coupler `k(λ)`, 11 pts | same | **1.878e-15** |

The coupler's are float noise, not exactness, because the body is now reached
through `EMESolver.solve_section()` at the one real length rather than
`solve_section_length_sweep()` at 200 — with `neff_override` both take the
same diagonal fast path, so the residual is summation order.

#### The arm calibration — MEASURED, job 11995681 task 1

`UniformCrossSection([(0.0, 0.5)], length=6.0)` through
`FDTDSolver.section_mode_beta_details()`, the same route
`ring_resonator_hybrid.py` uses for its own straights. Merged domain
5.1 × 5.0 × 2.5 µm, fit lengths [1.365, 2.731, 4.096] µm, 47,604 timesteps to
t = 396.7, **48.4 min on 24 ranks** — the only thing in the job that
time-stepped; every other piece hit cache.

**`ARM_RESOLUTION = 60`, matching the HCG mirror rather than the coupler's 40
— the user's explicit call.** There is no single device resolution to inherit
here, and a bare strip is cheap enough at 60 that the dominant lever on the
comb is measured on the finer of the two meshes this device already uses.

| at λ = 1.55 µm | n_eff | excess over MPB |
|---|---|---|
| arm, MPB mode index | 2.138038 | — |
| arm, FDTD's own accumulated phase (res 60) | **2.169312** | **+1.463 %** (band +1.442 % … +1.475 %) |
| coupler body even, MPB / FDTD (res 40) | 2.128829 / 2.181649 | +2.481 % |
| coupler body odd, MPB / FDTD (res 40) | 2.099957 / 2.156245 | +2.680 % |

**An unplanned result worth more than the planned one.** The arm's excess at
res 60 (+1.46 %) is a little over *half* the coupler body's at res 40
(+2.48 %/+2.68 %) — on a simpler cross-section, but the same MPB-vs-meep
comparison. That is direct evidence the gap is discretization rather than
anything physical, which is what the whole calibration argument assumes and
what `section_mode_beta_details()`' own docstring asserts from the coupler
alone. **It is also the stated cost of the res-60 choice**: the two blocks in
this cavity now sit on *different* discretization errors, so a milder version
of the inconsistency this step removes survives — over the coupler's 23 µm
instead of the arms' 106 µm. `HCGFP_ARM_RES` redirects output to its own
folder so a real convergence sweep can settle it without touching this.

#### What it does to the spectrum — MEASURED, both arms of job 11995681

**The lever arm is 4 × L_ARM = 212.32 µm, not L_ARM.** The resonant
denominator is `1 − ρ²t²` with `ρ = r·phase²`, so the arm index enters the
round trip through `phase⁴`. (The ring hit the analogous trap from the other
direction — job 11983058, `L_OUTER` vs the net arm.)

| | arms on MPB index (task 2) | arms on FDTD β (task 1, production) |
|---|---|---|
| resonances resolved, 1.50–1.60 µm | 44 | **45** |
| mean FSR | 2.268 nm | **2.245 nm** (−1.01 %) |
| Q range | — | **9,998 – 16,071** |
| T_min range | — | 0.026 – 0.309 |
| comb position | — | **shifted +26.9163 rad = +4.28 FSR** (≈ 9.6 nm) |
| `Circuit` vs closed form | agrees | max\|ΔR\| 9.71e-17, max\|ΔT\| 6.66e-16 |
| dense passivity, max R+T | 0.9955 | 0.9955 |

**This is not a refinement; it re-registers the comb.** The +1.46 % index
change is worth more than four free spectral ranges of position, while FSR
moves 1 % and Q/depth barely at all — the expected signature, since the
former is an absolute phase over 212 µm and the latter are ratios.

**Stated honestly, and stated in the script's own docstring and printout:
calibrating makes the loop self-consistent, it does NOT make the absolute
resonance positions converged.** Both MPB's index and meep's β are
discretizations of the same continuum problem and neither is the truth. What
changes is that one engine now carries every segment of a closed loop, which
is the thing a round trip actually requires. The residual is a real,
unmeasured error on absolute comb position; FSR, Q and depth are
differences/ratios and are insensitive to it at the same order.

**The baseline is preserved, and is not a bit-identical replay.**
`HCGFP_CALIBRATE_ARM=0` writes to
`example_output/high_contrast_grating_fabry_perot/hcg_fabry_perot_hybrid_uncalibrated_arm/`
(21.6 s, all cache hits). It puts the arms back on an MPB index but asks
*meep's* eigensolver on the FDTD grid, not `meow`'s `Mesh2D` at res 100 as
the pre-rewrite script did — deliberately, so that holding the grid fixed
makes the two curves a statement about the two ENGINES rather than about two
meshes. That is why its FSR (2.268 nm) is near but not equal to the
2026-08-18 published 2.255 nm. **The original is archived byte-identical at
`archive/example_output_20260818/`** (verified by md5 before anything was
overwritten).

#### ⚠️ OPEN, and it cuts against this change — the measured FSR sits ABOVE both simulations

Found 2026-08-18 while auditing what consumes this device's outputs, **not**
part of the work above and not resolved by it.
`beef_manuscript_figures/fp_resonator_figures/HCG_plotting_LUNA.ipynb` is a
measured-vs-simulated manuscript figure, and its stored outputs (last executed
elsewhere, against the PRE-rewrite dense pickle) record real ground truth this
device's simulations have never been read against in this doc:

| | mean FSR, in the figure's own window | vs measured |
|---|---|---|
| **measured** — `right_chip_N_LUNA`, real fabricated device, 15 points | **2.2608 nm** | — |
| simulated, pre-rewrite (arms on `meow`/MPB res 100), 16 points | 2.2497 nm | −0.49 % |

Both are *windowed* statistics from that notebook's own peak finder, so they
are not directly comparable to the band-wide means above — but the offset is
calibratable: the same pre-rewrite data gives 2.255 nm band-wide against
2.2497 nm windowed, ≈0.005 nm. Carrying that across:

| band-wide (this doc) | → windowed estimate | vs measured 2.2608 nm |
|---|---|---|
| uncalibrated arms, 2.268 nm | ≈2.263 nm | **−0.1 %** |
| **calibrated arms (production), 2.245 nm** | ≈2.240 nm | **−0.9 %** |

**So calibrating the arms moved the simulated FSR AWAY from the measured
value, by roughly a factor of two.** Stated plainly rather than buried: the
real device's round-trip optical path is *shorter* than either model's, and
raising `n_eff` by +1.46 % lengthens it further. Consistent with meep's Yee
operator biasing `n_eff` high — which is exactly what the +1.46 % excess is —
and therefore consistent with MPB's index happening to sit closer to reality
on this cross-section.

**Why this does NOT by itself say the calibration is wrong**, and why it is
recorded as open rather than acted on:

1. **The cavity-length model has errors of comparable or larger size.**
   `L_ARM = 53.08 µm` is the *average* of two real arms that differ by
   **3.0 µm** — 2.8 % of the 106 µm arm total — and each arm's 90° bend
   (R = 10 µm, 15.708 µm of arc) is **treated as straight waveguide** per
   Step 1's own framing, though a bend's effective index is not a straight's.
   Either approximation can move FSR by ~1 %, i.e. as much as or more than
   the calibration does. Absolute FSR on this device is currently limited by
   cavity-length modeling, not by which β the arms carry.
2. **The self-consistency argument is about the loop, not about absolute
   accuracy.** A closed round trip through blocks measured by two different
   discretizations is not measuring the device, whatever the FSR comes out
   at. That argument stands unchanged.
3. **The arm β is not converged** (one resolution, res 60 — see the open item
   above). A res sweep could move the +1.46 % excess materially, and the gap
   already shrinks with resolution (+2.48 % at res 40 on the coupler body).

**What would actually settle it**, cheapest first, all with every other piece
already cached: (a) an `HCGFP_ARM_RES` sweep, to see whether the excess
extrapolates toward zero — if it does, the converged answer may land near the
MPB curve and the tension dissolves; (b) modelling the bends as bends rather
than straights, and/or running the two real arm lengths asymmetrically, to
size approximation 1 against the calibration; (c) only then re-reading the
measured comparison. **Nothing here should be read as validating either
convention against the LUNA data until (a) and (b) are done** — the two
candidate explanations are not yet separated.

**Downstream consequence, unactioned**: that manuscript notebook's
`fp_resonator_fsr.{svg,png}`, its FWHM/extinction-ratio tables, and any
figure built on `smatrix_fabry_perot_dense.pkl` will change when re-run —
the comb moved **+4.28 FSR** and the FSR by −1.01 %. It has NOT been re-run
here; its stored outputs still describe the pre-rewrite simulation.

#### A cache-key bug found on the way, and fixed

The mirror's `solve_reflection()` call missed a cache entry that visibly
matched it. Cause: `FDTDSolver.minimum_run_time` (commit `d7910e8`) is
tracked inconsistently — over-included for `solve_reflection` (no
`_METHOD_ATTRS` entry, so the full `vars()` dump), stranding all 48
`solve_reflection` entries including this mirror's ~55 min res-60 3-D
production run; and **under**-included for `solve`/`solve_section`/
`solve_section_multi_length`, which all pass it to `stop_when_dft_decayed()`
— a silent false-positive hit, live on the modal multiplexer (100.0 and
1000.0, both added 2026-08-08 to fix a real early-stop bug). Fixed in
`beef/cache/config.py` (`_DEFAULT_OMITTED_ATTRS`); the mirror hit cache again
with **no re-simulation**. Full write-up:
`section_mode_beta_multi_monitor_validation.md` step 11.

#### Consequences for how this example is run

- **No `--dependency` any more.** It is not downstream of either constituent
  script; `batch/run_hcg_fabry_perot_hybrid.batch` submits standalone.
- **The Makefile target moved from `$(my_python)` to `$(runner)`.** It used to
  be login-node safe (two pickle loads plus one cheap mode solve); a cold run
  is now real 3-D FDTD.
- **The output folder is self-contained.** Both calibrations' merged-domain
  and per-mode fit figures, the mirror's domain/reflection plots, and the
  coupler's hybrid-section plots are all *generated* here from the same
  geometry/solver objects that were simulated, rather than copied out of
  another script's folder — the picture cannot drift from what ran.
- **The companion notebook was regenerated from the script** and executed
  end-to-end under `CACHE_ONLY = True`, reproducing these exact numbers
  (+1.463 %, +4.28 FSR, 45 dips, FSR 2.245 nm).
  `scripts/check_notebook_sync.py`: 18 pairs, 0 out of sync.

## Measured vs argued — read this before quoting anything above

**MEASURED (read from source or real MPB mode solves, no FDTD)**: A1's
platform mismatch (both scripts' own module-level constants read
directly); A2's characterization of `beef.network`'s existing primitives
and their current lack of a cross-device, `Circuit`-based caller (read
`beef/network/*.py` and grepped every `examples/**/*.py` for
`Circuit`/`connect(`/`close_loop(` — one non-test hit,
`verify_ring_network_crosscheck.py`, using the free functions directly on
one device's own coupler section, not `Circuit`, and not composing
separate devices); A3's `from_smatrix()` limitation (its own docstring)
and `solve_reflection()`'s complex (not just power) coefficient (read
directly from `beef/solvers/fdtd.py`); **Step 0's real dimensions from
`duty/dc.gds`** (`gdstk`-read polygon extents, cross-checked against
`beef/geometry/coupler.py`'s exact parameter semantics); **Step 0's mode
counts** (real `meow`/MPB mode solves, both the bare strip and the coupled
pair, at 3 C-band wavelengths — 4 supermodes exist at the cross-section);
**Step 0's mode-conversion measurement** (a real `res=25`,
`decay_tol=1e-10` `FDTDSolver.solve_section()` run of the actual S-bend,
mode 1 excited, decomposed onto all 4 supermodes — <0.01% reaches the
2nd-order pair); **Step 1's group-index table** (finite-difference from
the same 3 measured `n_eff` points); **Steps 3-5's synthetic-data
verification** (`tests/test_fabry_perot_network_crosscheck.py`, 3 tests,
all passing) — the `from_transfer_matrix()` bug and its fix, the 5-block
`Circuit` wiring, and the closed-form formula are each checked against
either a brute-force linear solve or each other, to floating-point
precision; **Step 2's real DC coupling** (job 11716260, `T_through=0.8989`,
`T_cross=0.0979` at λ=1.55 µm, complex `t(λ)`/`k(λ)` at 11 wavelengths,
`decay_tol=1e-11` not loosened); **Step 3's real waveguide `n_eff(λ)`**
(21-pt `meow`/MPB sweep — SUPERSEDED 2026-08-18 by Step 7's FDTD β
calibration, job 11995681: the arms propagate on `n_eff = 2.169312` at
λ=1.55 µm, not MPB's `2.138038`) and the real `rho(λ) = r(λ)·p(λ)²` built from it
plus the HCG's own real production spectrum; **Step 5's real evaluated
sparse spectrum at the 20 µm placeholder length (2026-08-07, superseded)**
— the closed form fed entirely real numbers, giving a real resonance at
λ=1.570 µm (transmission 0.7395 vs. an elsewhere-flat ~0.97–0.99,
reflection 0.1710 vs. ~0.001–0.02), passive everywhere (`R+T≤0.9953`);
**Step 5's real DENSE spectrum at that same placeholder length** — 22 real
resonances across the C-band via interpolation of the safe intermediate
transfer-matrix pieces (per-factor phase steps directly measured: `M_in0`/
`M_in1` ~65-69°/10nm, `beta` ~15-16°/10nm, both comfortably safe; only
their combination, ~153°/10nm, looked dangerous), validated against 2 real
601-point pieces that finished before job 11726827's timeout (≤1.3% error,
genuine ground truth); **Step 1's arm length, 2026-08-08** — traced
instance-by-instance from `duty/Mux_lpbk_9$4.gds`'s real fabricated
mirror-DC-mirror test structure (bend radius `R=10.0 µm`, confirmed two
independent ways), giving `L_ARM=53.08 µm` (average of two real arms,
54.583/51.583 µm, a real 3.0 µm asymmetry); **Step 5's real DENSE
spectrum at the real length** — re-run at `L_ARM=53.08 µm`, giving **44
real resonances, mean FSR 2.255 nm, Q≈9990-16000**, the interpolation
scheme and its validation unaffected (only the connecting-waveguide phase
term depends on `L_ARM`); **Step 1's FSR correction** — the original
~9.5 nm estimate had a missing round-trip factor of 2, caught by the
first (20 µm placeholder) real dense measurement before the length itself
was corrected.

**PROPOSED / ARGUED, not measured — NONE remaining for this device's own
core physics.** Step 2's `FDTD_RESOLUTION=40`/`EME_RESOLUTION=80`, once
ARGUED (assumed from a different device on a different platform), is now
MEASURED (check 3, 2026-08-08) — converged at the production value. The
only remaining ARGUED item in this whole doc is A4's own directional
question (whether HCG grating staircasing vs. this cavity's lower finesse
nets out favorably or not for a full-device FDTD check) — moot, since
check 4 is out of scope by decision regardless of which way that would
resolve.

**OUT OF SCOPE BY DECISION, not a measurement gap**: check 4, the
full-device FDTD control run — the user's call (2026-08-07, see A4/Step 6)
was to skip it outright, not attempt even a cheap pilot, given the ring
resonator's own precedent (its own full-device reference was disqualified
outright, not just slow, for reporting T>1 in a passive device and a
non-converging comb) stacked on the likely cost (this cavity's estimated
finesse and mirror-plus-grating geometry make the ring's own 167 h a floor,
not a ceiling). This device's validation rests on the single-pass/
algebraic checks only — same resolution the ring resonator's own doc
reached, for the same reason.

**What IS now a real physical result for the combined device**: the
closed-form Fabry-Pérot-with-DC-tap spectrum in Step 5, built entirely
from real measurements (HCG mirror, connecting waveguide, DC coupling),
independently reproduced by the literal `Circuit`/port-relabeling
assembly (Step 4, `hcg_fabry_perot_hybrid.py`) to floating-point precision
— AND now fully resolved at a dense (20001-point) grid via interpolation,
AND at the real cavity length (Step 1, traced from a real fabricated test
structure, not a placeholder) — giving 44 real resonances with real
FWHM/Q/depth each, validated to ≤1.3% against genuine partial ground
truth; check 2 (mode identity/ordering) is closed — measured identity
correspondence, worst `|dn_eff|=0.00048`; check 1 (engine agreement) is
now CLOSED — the ~12% `Δβ_FDTD`/`Δβ_MPB` discrepancy `CALIBRATE_BODY`
compensates for is confirmed smooth and stable (0.28% relative variation)
across the entire band, AND the stub/cut-point invariance test (job
11747764) confirmed `CALIBRATE_BODY` makes the arm and body sections
propagate on the same real FDTD-measured rate (0.3% agreement, branch-
safe) — the small residual drift measured is fully explained, not an open
gap; check 3 (resolution convergence) is CLOSED — `res=40→50` moves
`T_cross` 0% while `res=30→40` moves it +7.35%, real confirmation the
production `res=40/80` sits at or past the convergence knee. **All four
checks now have a real, resolved status** (check 4 out of scope by
documented decision, the other three closed with real measurements) —
nothing left open in this device's own validation.
