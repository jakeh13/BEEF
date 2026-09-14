# Validation docs

One file per example device, recording what has to hold before that device's
**hybrid** result can be trusted, and what is currently missing.

    resolved/mmi_1x2_validation.md              1×2 MMI
    resolved/mmi_2x2_validation.md              2×2 MMI
    resolved/ring_resonator_validation.md       ring resonator
    resolved/mmi_1x4_sin300_validation.md       1×4 MMI on 300 nm Si3N4 — RESOLVED 2026-08-31. All four checks MEASURED (1: full-band engine agreement; 2: mode identity on a 2-mode access waveguide; 3: job 12275175; 4: jobs 12029541 band-centre + 12260631 full-band). Steps 1-7 done. Step 9a's beta-calibration hypothesis REFUTED. ⚠ THREE RESIDUALS CARRIED, stated in the doc's own RESOLVED header -- read them before quoting any number: (1) check 1's engine disagreement is REAL and unexplained, and LARGER than the precedents already in resolved/ -- opposite spectral slopes (-0.749 vs +0.087 /um), EME-vs-hybrid split shape 20x apart WITH OPPOSITE SIGN, total T differing 2.5 points; mode truncation tested and REFUTED as the cause, confound check done first; (2) the leading remaining suspect, the EME taper staircase, was DROPPED at the user's instruction and never tested -- dropped is not excluded; (3) res=30 is NOT established for the EME side (its differential drift 30->40 is 2.17e-03 vs FDTD's 1.99e-04), and mmi_1x4_eme.py still runs at 40 while _fdtd/_hybrid moved to 30 -- possibly correct, undocumented either way. Open to-dos: the nitride strip-test gate substantively PASSES (arm B dn_eff=+0.000207, 10x inside tolerance) but its harness self-check uses an absolute bar calibrated on silicon that a low-contrast core cannot clear -- threshold deliberately NOT changed, needs a decision; step 9b scoped not run. ★ UPDATE 2026-09-05 to -09-14: step 9d (interface_modes sweep), previously scoped as a multi-week cost, WAS RUN after a monitor_size cache-pin fix cut it to 25-32s/mode-count -- this is the HYBRID's own error vs the FDTD reference, a DIFFERENT question from residual 1's EME-vs-hybrid disagreement (still open, untouched). Blue half (1.50-1.55um) SOLVED at interface_modes=22-24 (worst per-port err 5.99%->0.62%; mechanism confirmed: single TE input excites only even lateral orders, wasting most of a 12-mode basis on ~1e-11-power slots). Red half (1.56-1.60um) untouched -- physically insensitive to mode count, the modes that would help are cut off there. body_stub cut position independently ruled out as a lever (0.07% effect, all 4 ports, gated at 0.12% vs shipped). C-band (1.530-1.565um, 21pt) FDTD reference resim COMPLETE, no hidden pathology, resolves real beat structure the coarse grid aliased by luck; matching hybrid resim still running. Two ~3-day compute losses to a recurring dmesg-confirmed Lustre os.getcwd() failure (7 different nodes across 3 instances) -- the 16-mode CALIBRATED confirmation was dropped after its 2nd loss rather than retried a 3rd time. See mmi_1x4_sin300_validation.md's "9d RESULT" section and running_jobs.md's 1x4 ledger

    resolved/fdtd_vs_EME_mode_validation.md     cross-device — the mode-correspondence MECHANISM itself
    resolved/calibrate_body_neff_mode_solve_skip_validation.md  cross-device — CLOSED 2026-08-11; EMESolver.solve_section()'s neff_override fast path (skips compute_modes()/compute_s_matrix() when calibrate_body_neff makes them moot) built, measured (mode solve was 99.7% of solve_section's cost), and verified bit-exact against the reconstructed old path on the real 2x2 MMI body
    resolved/gds_geometry_hybrid_cascade_validation.md  cross-device — CLOSED 2026-08-21, all 5 steps MEASURED. GdsGeometry could back a whole-device single-solver run but could NOT occupy a cascade section (section_prisms/section_gds/width_at inherited unimplemented, base.py:155-192; center_at silently returned 0.0, wrong for any off-axis device). All four hooks plus source_aperture_at and a sampled _check_cut now exist (gds_import.py), verified against a GDS twin generated from MMI.polygons() at the shipped 2x2 dimensions so the comparison isolates the QUERY path rather than testing two hand-drawn layouts. Geometry agreement is float noise (width_at/center_at 4.44e-16/1.11e-16 um over 1401 x; section polygons 0 mismatched pixels on a 1600x400 raster; all 13 cut-legality verdicts matching MMI's). 44 tests. ★ ACCEPTANCE (job 12069544, 4:31:02): T_bar/T_cross agree to 1.6e-03 relative — LOOSER than the "solver noise" bar this plan set, and the whole residual is EXPLAINED, not a GDS defect. Localized in three cache-only measurements: (1) beta/body_neff bit-identical AND computed fresh on the GDS side, which independently re-verifies the meep/FDTD layer and closes the library-drift caveat job 12070199 could not; (2) M_in1 carries 100% of it (swapping neff alone moves 0.000e+00); (3) cause is PRISM ORDER — both backends emit the SAME five polygons for the input section in a DIFFERENT order, and meep's subpixel averaging is order-sensitive at abutting-prism seams. Reordering GDS's polygons to MMI's makes the epsilon grid BIT-IDENTICAL (0 of 2.52M voxels; 872 differ natively, all at the two seams). Recorded, not fixed, at the user's direction. ★ SPUN OUT, LIBRARY-WIDE: merging the five abutting polygons into the one they geometrically are moves 41,677 epsilon voxels — 48x the ordering effect — so EVERY BEEF device drawn as abutting prisms carries a subpixel-averaging artifact of that size. Wants its own plan; would move all six devices' numbers, invalidate the 566-entry FDTD cache and change the manuscript's published 2x2 results
    resolved/multi_excitation_cascade_validation.md  cross-device — CLOSED 2026-08-21, all 4 steps done; the library-integration half of notes.md Solver recommendation #2. mmi_2x2_hybrid.py and directional_coupler_hybrid.py both bypassed HybridSolver.solve() for a hand-rolled einsum, because their cross path pairs one arm's transposed matrix with the OTHER arm's excitation — so neither could use calibrate_body_neff. MEASURED finding that framed the whole plan: routing the 2x2 through the flag changes not a single number (_calibrated_neff and the hand-rolled path compute the same quantity from the same source, boundaries, n_modes and conversion). An API question, not a physics one. ★ BOTH DEVICES NOW MIGRATED AND BIT-IDENTICAL to what they replaced (2x2: 0.000e+00 on T_bar/T_cross AND the phase-sensitive complex coefficients vs both the shipped run and step 1's golden; DC: 0.000e+00 on T_through/T_cross vs its shipped length_sweep.npz), each at the SAME 1 FDTD simulation the hand-rolled path cost. Two library additions made it possible: MirrorSpec (declares a cross-excitation mirror — reuse section j reinterpreted as the other arm via a parity vector — backed by a SECOND geometric check that the partner's two excitations really are transverse mirror images, so partner_kwargs cannot become an unchecked escape hatch) and calibrate_kwargs (kwargs for the calibration source only, never the section solver — the DC needs an interface monitor on one and not the other, and EMESolver funnels unknown **kwargs into geometry_kwargs, which is in the cache fingerprint). ★ TWO PLAN ASSUMPTIONS WERE WRONG, both favourably: beef.network.Circuit turned out NOT to be needed (once MirrorSpec exists each excitation is an ordinary single-excitation solve_length_sweep call with nothing left to compose) and would have COST accuracy (BLAS-vs-einsum, ~4 eps); and the "bit-identical" bar step 1 relaxed as unreachable was reachable after all, via solve_length_sweep rather than Circuit. Also measured: the BAR PATH CANNOT DETECT A TRANSPOSED BODY (a quadratic form — B and B^T give the identical scalar), so a bar-only check like the pre-existing test_network.py one would miss it. 30 tests. ★ MODE_SIGNS VALIDATED (job 12159086, 45:01): an independently measured M_in0 confirms all 14 signs; the parity relation constrains SIGNS, and its amplitude half degrades only on noise-floor modes (modes 0-2 agree to 1.7e-05, mode 8 — |M|=1.3e-06, 5e5 below the dominant mode — to 2.8). Cost of the whole SYMMETRIC_MODE_REUSE assumption end to end: 1.585e-05 on T_cross, saving a 45-min 3-D run. OPEN TO-DO: step 3(b), deriving MODE_SIGNS from mode-field parity instead of declaring it — unblocked (cascade_golden_real_m_in0.npz is the regression target) but DE-PRIORITISED, since its motivating hazard (a stale vector silently corrupting the cross path) is now covered by MirrorSpec's per-call geometric precondition and the vector is measured correct
    running_jobs.md                    snapshot — active SLURM jobs and which gap/item each is for
    mmi_4x4_90_hybrid_validation.md    4×4 MMI 90° optical hybrid, SOI 220nm/O-band — BUILD PLAN, device does not exist yet; ref: Yasui et al., ACES J. 36(5), 2021
    mmi_1x4_tapered_body_validation.md 1×4 MMI with a TAPERED multimode body, 300 nm Si3N4 — BUILD PLAN, device does not exist yet and the library cannot propagate it. Opened 2026-08-20. Depends on mmi_1x4_sin300_validation.md for material model + cross-section design; becomes the manuscript's 1x4 test case, with the straight-body 1x4 kept as the untapered baseline. Headline: EMESolver.solve_section() is single-cell by construction on BOTH branches (neff_override fast path writes a diagonal T, eme.py:541-557; the else branch still calls create_cells(..., Ls=[length]), eme.py:565), so no tapered body can be propagated correctly through a cascade — and nothing guards it, since the fast path gates on `neff_override is not None` rather than uniformity and _calibrated_neff (hybrid.py:1272-1321) applies it to any EME body untested. Step 1 (reuse _validate_uniform_sweep_section as a shared guard) removes a silent-wrong-answer path and should land regardless of the device
    mmi_2x2_waist_body_validation.md   2×2 MMI with an interior-WAIST multimode body, SOI 220nm/C-band -- ✅ VALIDATED 2026-09-14. DESIGN POINT Ww2.24 Lw6.0 Lt5.0 L22.0 = 31.00um total, 16.2% shorter than the 37.00um straight control, confirmed by full-device FDTD: IL 0.2010 dB, imbalance 0.0124, phase mean-dev 0.267 deg / MAX-dev 0.447 deg -- it passes all four criteria under BOTH the mean and max phase readings, so it does not depend on that ambiguity. The 1x4's waist failure does NOT reproduce: this device works, which REFUTES the radiative-loss-at-transitions hypothesis as a general claim. ❌ The earlier 29.50um/20.3% headline is RETRACTED (fails on both settled taper bases AND full-device FDTD; also never hybrid-buildable at s=0.25um). Grid re-scored at a settled taper basis and verified stable across two independent bases (660/668 verdicts agree, 70 passers). ⚠ C_TRANS is an OPTIMUM not a floor on a tapered body -- step 5a's 'converged at 8' does not hold for the hybrid, where more cells ADD spurious loss (~2% |E|^2 discontinuity per ramp cell boundary). METHOD (separate from the design claim): hybrid matches full-device FDTD to 0.76% on the bar port; a ~1-2% CROSS-port deficit persists on both geometries tested, with 5 hypotheses measured and 4 falsified. Per-cell beta calibration for TAPERED bodies added to the library (cell_neff_override/cell_neff_fn)
    resolved/high_contrast_grating_fabry_perot_validation.md HCG Fabry-Perot cavity + DC tap — Steps 0-5 + checks 1-3 all closed with real data; 44 real resonances at the real (GDS-derived) cavity length, FSR 2.255nm; check 4 out of scope by decision
    modal_multiplexer_validation.md    modal multiplexer — resumed 2026-08-04, gap list against the four checks
    resolved/section_mode_beta_multi_monitor_validation.md  cross-device — steps 1-8 CLOSED 2026-08-16, then SUPERSEDED IN PART by step 9 on 2026-08-17. ★ READ STEP 9 BEFORE TOUCHING ANY DEVICE'S β CALIBRATION — it is self-contained and is the standard way calibration is now done. Headline: a β calibration's cost is 89% MPB modal decomposition and 11% FDTD time-stepping (measured, job 11872191: 7,249 solves @ ~13.9 s vs 1.4M timesteps @ 0.0089 s), and ~14x of that decomposition was computed and discarded. measure_diagonal_only=True fixes it exactly (~4x, verified bit-identical against cached matrices); solve_section_multi_length() gained a cache key and is now the DEFAULT path (a further ~1.9x, not the 3x steps 1-8 projected — merging leaves the measurement work identical and only removes decay waits and norm/source solves). MMI 2×2: 30.6 h → ~3.7 h. THIRD change same day: the MPB seed (section_mode_neff) is now cached AND solved at seed_anchors=3 wavelengths with the rest interpolated — it had become ~39% of the shrunken calibration at 154 eigensolves/~31.5 min; now 42 solves/~8.6 min, and free on any re-run. Only the seed is subsampled (FDTD and decode keep every wavelength); it perturbs the reported seed_excess by 0.0104 pp. Seed anchoring VERIFIED against the real 11-wavelength MPB seed (job 11940556): worst error 1.53e-04 in n_eff = 2.15e-04 of the branch half-spacing against margins of 0.05-0.095, fit lengths bit-identical so no cache entry is invalidated, and a measured convergence curve (3/5/7 anchors -> 1.53e-04/3.01e-06/3.60e-08). Steps 1-8's cost model is wrong; its mechanism and convergence work still stand. Per-solve cost does NOT improve with rank count. Verification jobs: 11937894 (merged calibration, production), 11940556 (seed anchors vs the real 11-wavelength seed). ★ STEP 10 (2026-08-18) is the audit that closes the rollout: all six step-9 devices re-run on ONE library state and bit-identical (jobs 11981703/11982196, <90 s each — the cache hits are the evidence the seed cache/seed_anchors/merged-plot edits perturbed no key); a HybridSolver bug found where .last_calibrations was wiped by a second solve() served from a shared in-memory cache, so mmi_1x2_hybrid (and mmi_1x4_hybrid, not yet re-run) produced NO calibration figures on a calibrated run — fixed + 4 tests; and the ring resonator brought onto the calibration for the first time ever (it had never been asked the question — neither script goes through HybridSolver, so the rollout swept past it), via UniformCrossSection + section_mode_beta_details(), covering the bare strip in both scripts and the coupled-pair P_body in the cascaded one. Ring stage 1 landed: excess +1.044% against the res-100 meow curve, +2.259% against a res-40 MPB seed — which reproduces item 6's +2.233%/+2.361% on a DIFFERENT stack (SOI_220 vs SOI_160), the first independent confirmation of that gap. Job 11983058 still running
    archive/                           closed items, kept whole
    resolved/                          whole docs sufficiently addressed (open to-dos may remain)
    resolved/directional_coupler_validation.md  directional coupler — the worked example
    resolved/material_model_validation.md       cross-device — shared Si/SiO2 model, both engines
    deferred/                          work intentionally postponed, not blocked or closed
    deferred/high_contrast_grating_validation.md HCG reflector, FDTD-only (no EME/hybrid) — moved here 2026-08-13; corrected peak R 97.1% vs paper's 98.9%, gap widened after the waveguide-extension fix; res=70/80 probe shows resolution residual has NOT settled by res=60; postponed pending a decision to spend the multi-hour cluster time on a res=70/80 full-band re-run
    deferred/mode_crossing_correspondence_repair_validation.md  cross-device — BUILD PLAN, moved here 2026-08-13; check_mode_correspondence can DETECT a mode-identity mismatch but nothing REPAIRS one, prompted by the modal multiplexer's interface 0 (mode 8, wl=1.26 um); Stage 0 classified the failure as regime C with a real mechanism (meep's per-band-count ladder re-diagonalizes independently each call, no shared state — confirmed via meep's own C++ source); Stage 2 feasibility CONFIRMED (a single MPB simultaneous num_bands solve keeps the known duplicate pair orthogonal, overlap 0.34 vs meep's 1.0) but production integration (frequency-matching per band, field-API adapter) not started — postponed by user choice, not blocked

Every per-device file for an **existing** device lives in `resolved/`.

**`mmi_1x4_sin300_validation.md` (2026-08-03) is the first doc here that is a
plan rather than a gap list**, and it sits at the top level deliberately: the
device does not exist, so `resolved/` (sufficiently addressed) and `deferred/`
(postponed by choice) both misdescribe it. It orders the work to build a 1×4 MMI
example — FDTD, hybrid and EME — on a 300 nm silicon nitride core, and converts
to a normal four-checks gap list once those three scripts run. Its headline
finding is that the 1×4 *geometry* is already supported (`MMI(n_inputs=1,
n_outputs=4, …)`) and the real work is in the shared material model:
`beef/materials/stack.py` is hardwired Si-core/SiO2-clad, and `_SELLMEIER` has no
nitride. Read its measured-vs-argued section before quoting any number from it —
the material model numbers are measured, every device dimension is not.

**`high_contrast_grating_validation.md` (2026-08-04) is also a plan, not a
gap list, but for a different reason than the 1×4's** — the device (loaded
from `examples/gds_files/HCG_seeded_TO.gds`) exists and runs; nothing in
the [four checks](#the-four-checks) below applies to it, because it is
FDTD-only by design (no EME, no `HybridSolver`, per the user's request) and
because the GDS turned out to be a **one-port** device (a single access
waveguide feeding a freeform topology-optimized mirror region, not a
two-port through structure) — `FDTDSolver.solve()` cannot report anything
physical for it at all. Two real library gaps blocked using it with FDTD
in the first place (`gdstk` missing from the `jmh` env; `GdsGeometry` had
no `polygons()`, only the EME-only `to_gds_eme()`), both fixed this
session, plus a new `FDTDSolver.solve_reflection()` method for the
one-port reflection measurement itself.

**Same day, later: a real 3-D resolution + x/y/z PML convergence sweep
landed** (19 SLURM array tasks, `gts-sr108-ece`), after catching and fixing
two real mistakes — the whole sweep had first been run at the **wrong
band** (O-band, copying every other `SOI_160` example, when the paper this
device is from reports a **C-band** design peaked at 1569 nm), and the
x-domain (PML clearance) sweep was briefly self-corrupting (`source_margin`/
`monitor_margin` are PML-edge-relative, not geometry-edge-relative, so
sweeping `edge_clearance` alone moved the source outside the physical
waveguide). A full-band production run at the settled `res=60` then landed
(job 11670282) with peak R=98.10% at 1.550 µm against the paper's 98.9% at
1.569 nm — agreement within ~1 point.

**Evening, same day: a THIRD real issue, caught by the user directly** —
the input waveguide ended abruptly instead of continuing through the PML,
a genuine physical defect (the clean-cut facet reflects part of the guided
mode on its own, contaminating exactly this reflection measurement).
Fixed (`GdsGeometry` gained `input_width`/one-sided `stub_length`
extension; `FDTDSolver.solve_reflection()` gained
`extend_input_waveguide=True`), verified visually via a rendered ε-map, and
the full 17-point convergence sweep re-run: **systematically lower R at
every resolution** than the pre-fix numbers (the old numbers were partly
inflated by the facet artifact), same clean convergence shape, `res=60`
now at R=0.9697 (vs. the pre-fix 0.9789 — which, in hindsight, matched the
paper's 97.9% partly by the artifact coincidentally cancelling against the
resolution residual, not because the pre-fix pipeline was more accurate).

**The re-run production run has landed (job 11673560, 54 min) — and the
honest number is worse, not better.** Corrected: **peak R=97.14% at
1.580 µm, a 1.76-point gap from the paper's 98.9% — larger than the
(flawed) pre-fix run's 0.8-point gap.** The old, better-looking number was
partly an artifact: removing the facet reflection also shifted the
resonance peak by ~30 nm, from 19 nm blue of the paper's 1569 nm to 11 nm
red of it — a real position change, not just a magnitude correction, that
the pre-fix comparison could not have revealed. Reported as the honest
current result rather than smoothed over. See that file's step 9 for the
full spectrum and step 10 for what it opens — chiefly, whether resolution
above `res=60` keeps moving the resonance position (this session's
strongest open question), plus whether this GDS is the paper's final
design or a staged seed, where the non-reflected power goes, and
generalizing the domain plot.

**Later still: two of the three caught mistakes now have an automated
guard.** `FDTDSolver.solve_reflection()` gained a pre-flight
`check_reflection_layout()` call (before `mp.Simulation` is built, so a bad
layout fails in milliseconds) that verifies the source and monitor land
inside the geometry's real silicon, and — when
`extend_input_waveguide=True` — that the silicon actually reaches the
domain boundary. Both checks are reconstructed against this device's own
historical bug-2/bug-3 numbers in `tests/test_fdtd_reflection_layout_guard.py`
(7 new tests, no meep import needed), and confirmed clean against the real
production layout. Full suite: 625 passed, 19 skipped. The wrong-band
mistake (1 of 3) still has no guard — nothing in the codebase represents a
device's intended design band to check a request against — see that file's
step 10 item 6.

**2026-08-05: the peak-position open question was partially checked, and
the answer is not reassuring.** A single-wavelength `res=70`/`80` probe
(job 11688701) added to the resolution series shows R at the fixed
1.569 µm probe wavelength still climbing by a non-trivial amount through
res=80 (res=60: 0.9697 → res=70: 0.9741 → res=80: 0.9780), and the
deceleration ratio that looked clean through res=60 stalled rather than
continuing to shrink — the resolution residual has evidently NOT settled.
Since this is only a single-wavelength probe, it can't directly confirm a
peak-*position* shift, but a magnitude change this size (comparable to
what the waveguide-extension fix itself caused, which produced a
documented ~30 nm peak shift) is real evidence the full-band lineshape
likely hasn't settled either — making the open question MORE open, not
closed. A real res=70/80 full-band production re-run is what would
actually settle it; not attempted (cost scales from the single-wavelength
point's own 262 min at res=80).

**`high_contrast_grating_fabry_perot_validation.md` (2026-08-05) is a build
plan for a new, hybrid (FDTD+EME) device — an HCG-mirror Fabry-Pérot cavity
tapped by a directional coupler** — reusing the HCG mirror above, composed
via `beef.network.Circuit` (a cross-device composition layer no example
has used before) with a directional coupler redesigned for `SOI_160`
(the two existing devices sat on incompatible platforms — HCG `SOI_160`,
the shipped DC example `SOI_220` — a real decision the plan flagged rather
than picking silently). **Step 0's platform decision landed same-day**:
reuse a REAL fabricated coupler's own dimensions, `duty/dc.gds` (read
directly with `gdstk`, fed into the existing parametric `DirectionalCoupler`
class — not a `GdsGeometry` import), rather than a freshly-derived
cross-section. **That real device's cross-section supports 4 supermodes,
not 2, at 160 nm in C-band** (measured via real MPB mode solves) — but a
follow-up measurement (exciting the real S-bend with mode 1, a real
`res=25` FDTD run) found **<0.01% of the launched power reaches the
2nd-order-mode pair**, so the user's call was to treat the whole device as
single-mode (`N_MODES=2`) — a measured simplification, not an assumed one.
**Step 1 (cavity length) is now pinned for real, 2026-08-08** — 53.08 µm/arm
(average of two real arms, 54.583/51.583 µm), traced instance-by-instance
from `duty/Mux_lpbk_9$4.gds`, a real fabricated test structure with this
exact mirror-DC-mirror topology (bend radius `R=10.0 µm`, confirmed two
independent ways; user's instruction: treat bends as straight waveguide).
Not the original 20 µm placeholder.
**Step 2 (re-deriving the DC for `SOI_160` at `N_MODES=2`) landed for
real 2026-08-07, after timing out twice.** Job 11691379 (`-t 300`)
TIMEOUT — 5 of 6 expensive FDTD sub-solves had already cached cleanly.
Resubmit job 11706987 (`-t 240`) correctly cache-hit all 5 of those and
ALSO timed out anyway, apparently inside two nominally "seconds, not
minutes" `section_mode_neff()` calls whose real cost was never fully
diagnosed. Shown both timeouts, the user's call was to resubmit with a
much larger margin — job **11716260** (`-t 720`) **COMPLETED in 4h56m
real walltime**: `T_through=0.8989`, `T_cross=0.0979` at λ=1.55 µm (≈10%
coupling ratio), passive, body-calibration confirmed doing real work
(MPB vs. FDTD-measured Δn_eff differ ~12%, same class of disagreement the
`SOI_220` script's own docstring documents).

**Steps 3-5 also landed with real data.** A real bug was caught before
any real numbers were involved: `beef.network.convert.from_transfer_matrix()`
(planned for the connecting waveguide) is one-directional BY DESIGN and
silently zeroes a mirror's reflection when connected to it — fixed with a
direct bidirectional block. A hand-derived Fabry-Pérot-with-DC-tap closed
form was verified against a brute-force 4-port solve and the real 5-block
`Circuit` wiring, both to floating-point precision
(`tests/test_fabry_perot_network_crosscheck.py`, 3 tests; full suite: 628
passed, 19 skipped). **Fed entirely real numbers (HCG mirror + real
waveguide `n_eff` + real DC coupling), that closed form produces the
first real predicted cavity spectrum — and finds a genuine resonance at
λ=1.570 µm** (transmission drops to 0.7395 from an elsewhere-flat
~0.97–0.99; reflection spikes to 0.1710 from ~0.001–0.02; passive
everywhere). **Read this honestly, not triumphantly: the 11-point/10 nm
grid is comparable to the ~9.5 nm estimated FSR, so the resonance is
almost certainly under-sampled** — this establishes the model produces
real resonant behavior, not that the resonance is characterized.

**Now fully resolved, 2026-08-08 — via interpolation, not more FDTD
compute.** A first attempt to densify (job 11726827, a brute-force
601-point re-run) TIMED OUT at its full 12 h budget after completing only
3 of ~7+ needed expensive sub-solves — going 11→601 wavelength points
multiplies each one's own eigenmode-decomposition cost by ~55x. Fixed
instead by measuring which quantities are actually safe to interpolate:
the real cached data shows the mode-basis pieces (`M_in0`/`M_in1`, the
DC body's calibrated `beta(λ)`) each have modest phase steps (~65-69°/10nm,
~15-16°/10nm) safe to interpolate directly, and only their COMBINATION
(~153°/10nm) looked dangerous — because this device's own mirror-symmetry
reuse counts one S-bend's phase twice. Standard practice, confirmed
against the literature (Lumerical's own EME docs and coupled-mode theory
both describe exactly this "solve modes sparsely, handle phase
analytically" split), and validated for real: 2 of job 11726827's own
pieces finished before its timeout, and interpolating the sparse data
reproduces that real measurement to ≤1.3% error. **Result: 22 real
resonances resolved across the C-band, mean FSR 4.500 nm** — the original
~9.5 nm FSR estimate turns out to have had a missing round-trip factor of
2. The resonance near the original sparse sample (read as T=0.7395 at
1.570 µm) is real but was ~18× shallower than its true depth
(T_min=0.0415 at 1.5698 µm, FWHM=0.235 nm, Q≈6680) — confirming the
under-sampling concern. New permanent scripts:
`directional_coupler_soi160_hybrid_dense.py` and a new "section 8"
appended to `hcg_fabry_perot_hybrid.py`.

**Superseded again, same day — Step 1's cavity length is now real, not a
20 µm placeholder.** The paragraph above (22 resonances, FSR 4.5 nm) was
computed at that placeholder length. Traced instead from
`duty/Mux_lpbk_9$4.gds` (a real fabricated test structure sharing this
device's own mirror-DC-mirror topology), `L_ARM=53.08 µm` — average of
two real arms (54.583/51.583 µm), each `straight → bend(R=10.0 µm,
confirmed two independent ways) → straight`, bends treated as straight
waveguide per the user's instruction. Re-running the same interpolation
scheme (unaffected by the length — only the phase term changes) at this
real length gives **44 real resonances, mean FSR 2.255 nm, Q≈9990-16000**
— roughly double the resonance count and Q of the placeholder-length
case, consistent with a ~2.05× longer real cavity.

**Step 4 (the `Circuit` assembly) is now also done with entirely real
data, not just synthetic-verified** — new permanent script
`examples/high_contrast_grating_fabry_perot/python/
hcg_fabry_perot_hybrid.py` loads the HCG's and DC's real production
`SMatrix` files directly, runs a fresh ~5 min real waveguide mode-solve,
and assembles all five real blocks through the actual `Circuit` wiring.
Along the way, the DC's own port-mapping question was resolved by
checking rather than assuming: the `isolated`-excited branch (`k'`/`t'`,
never separately saved by `directional_coupler_soi160_hybrid.py`) was
reconstructed from cached solver pieces and found to match the saved
`in`-excited branch (`t`/`k`) to machine precision — confirming the
standard symmetric-coupler equality holds for this real device rather
than assuming it. **Result: the real `Circuit` assembly matches Step 5's
real closed-form spectrum to floating-point precision**
(`max|ΔR|=5.6e-17`, `max|ΔT|=4.4e-16`, asserted in the script itself) —
two independent code paths, same real data, identical answer.

**Separately: check 4 (full-device FDTD control run) was decided OUT OF
SCOPE outright, 2026-08-07 — not deferred, not piloted, skipped.** The reasoning (A4/Step 6) picked up a direct precedent from
`ring_resonator_validation.md`: that doc didn't just find its own
full-device FDTD check slow, it found the reference itself reported T>1 in
a passive device (up to +8.4%) and a comb that shifted 18% of an FSR
between two resolutions, and retired the check outright, with a mechanism
that generalizes — round-trip phase error is amplified ~F× by finesse, and
FDTD's staircased boundaries ARE a round-trip path-length error, so the
reference is subject to exactly what it's supposed to be checking. This
device's own measured finesse (F≈23) is ~5× more tolerant than the ring's
(F≈120), but its HCG mirrors are gratings — many staircased tooth-edges per
round trip, plausibly more discretization surface per bounce than the
ring's smooth bends. Net direction unmeasured either way — but between the
likely disqualification and the likely cost (the ring's 167 h is a floor,
not a ceiling, for a smaller/lower-finesse device), the user's call was to
skip it rather than spend SLURM time confirming either. This device's
validation now rests entirely on the single-pass/algebraic checks that
carry no round-trip amplification (Step 4's closed-form-vs-`Circuit` match
to floating-point precision, DC reciprocity, passivity) — the same
resolution the ring resonator's own doc reached, for the same reason.

**Check 2 (mode identity/ordering) is now closed, 2026-08-08 — measured,
not assumed.** Extended `examples/verification/python/
audit_interface_modes.py` (the same cheap, no-time-stepping tool already
used for the SOI_220 coupler and both MMIs, one of which it caught a real
bug on) with this device's own real geometry. Result: identity
correspondence, worst `|dn_eff|=0.00048` — no mode-mixing risk.

**Check 1 (engine agreement) is now partially closed too, same day, at
zero new compute.** `section_mode_beta()` already computes a
per-wavelength FDTD/MPB comparison internally but only ever returned the
aggregate result. Recovered the full grid by reusing already-cached
calibration phases with the library's own `beta_fit_from_phases()`: the
`Δβ_FDTD/Δβ_MPB` ratio is **0.8801 ± 0.00084 across the entire 100 nm
band** — the ~12% discrepancy `CALIBRATE_BODY` compensates for is a
smooth, stable, physically real effect confirmed over the whole
operating range, not a one-wavelength fluke.

**Both remaining real-cost checks landed, 2026-08-08 (job 11747764, 5-task
array, 37 min–4h55m per task).** `directional_coupler_soi160_hybrid.py`
gained `DC160_BODY_STUB` and `DC160_RES_SCALE` env-var overrides (both
force a single-wavelength `N_WL=1` mode — these test numerical
properties, not the spectrum). **Check 3 (resolution convergence): ✅
MEASURED, converged.** `res=40→50` moves the coupling ratio by 0.00%
while `res=30→40` moves it +7.35% — real, direct confirmation the
production `FDTD_RESOLUTION=40`/`EME_RESOLUTION=80` (carried over
unverified from the SOI_220 device since 2026-08-06) sits at or past the
convergence knee. **Check 1's stub/cut-point sweep: ✅ CLOSED, resolved for real.**
`T_cross` is stable to <1% across `BODY_STUB=0.1–0.75 µm`. The naive
SOI_220-style formula predicted ~23× more drift than measured — but that
formula was derived for the SOI_220 coupler's own PRE-calibration stub
sweep (before `calibrate_body_neff` existed, specifically to expose the
Δβ gap via MPB's uncalibrated body rate); this device already runs
`CALIBRATE_BODY=True`, so the comparison was apples-to-oranges. Checked
directly instead, with data already in hand: the arm sections' own real
FDTD phase-accumulation rate (branch-safe fit, avoiding the wrap
ambiguity a naive unwrap risks) matches the isolated calibration's rate
to 0.3% per mode — meaning the stub-dependent terms in the combined
coefficient cancel almost exactly BY CONSTRUCTION, because
`CALIBRATE_BODY` makes the arm and body propagate on the same measured
rate. The small residual drift is fully explained by that 0.3% gap. All
four checks now have a real, resolved status (check 4 out of scope by
documented decision, the other three closed with real measurements).

**`mmi_4x4_90_hybrid_validation.md` (2026-08-05) is a build plan for a new
device — a 4×4 MMI optical 90° hybrid on 220 nm SOI at O-band, with real
input and output tapers** — at the user's request, reproducing the topology
of Yasui, Sugisaka & Hirayama's ACES J. 36(5) 2021 paper (full citation in
the doc), but on a different (strongly-guided Si, not weakly-guided silica)
platform and band, so none of that paper's own numbers are a target (see the
doc's A9). Headline findings: the 4×4 geometry, and both input **and**
output tapers, need no new code (`MMI(n_inputs=4, n_outputs=4, …)` is the
same constructor every MMI example already uses); the paper's own
`y_i = y_{i+4}` port symmetry (Eq. 8) is already this codebase's default
port-construction behaviour, not a feature to build, and sets up the same
mirror-symmetry transfer-matrix reuse `mmi_2x2_hybrid.py` already exploits,
generalized from 2 to 4 ports; and — the one genuinely new gap — **nothing
in `beef.results.SMatrix` reads S-parameter *phase*, and this is the first
device in the repo whose own headline metric (the paper's phase-error, Eq.
15) needs it**, which also surfaces a new, unmeasured risk class no
power-only check here has ever had to cover (whether two independently-run
FDTD excitations share a consistent absolute phase reference). Deliberately
does **not** port the paper's own micro-GA optimizer (`beef/sweep/` has no
optimization loop) — self-imaging design equations (general, not
weakly-guided-specific) plus a coarse EME sweep are used instead, flagged as
a decision to revisit only if that approach can't hit the cited system
requirements. Nothing has been simulated yet; `SOI_220` and O-band have each
been separately validated elsewhere in this repo but never together.

**`ring_resonator_validation.md` moved into `resolved/` on 2026-08-03.** It carries
a `📁 FILED INTO resolved/` box at the top listing exactly what is closed, what is
out of scope **by decision** (absolute resonance wavelength, on a finesse argument;
comparison against a full-device FDTD, which reports transmission > 1 for a passive
device), and five open items deliberately carried in — most importantly that the
**monolithic** route keeps an unexplained k1/k2 residual and that the cascaded
route's reciprocity check is a *tautology* rather than evidence. Read that box
before quoting the file. As with the directional coupler below, citations elsewhere
name it bare; there are no markdown links to repoint (checked).

Check 2 below (mode identity) has grown its own cross-device file:
`resolved/fdtd_vs_EME_mode_validation.md` covers the *mechanism* — what must hold for any
geometry at any resolution pair, why n_eff matching and a diagonal-only overlap
check are both insufficient, the layered scheme recommended instead, and the
ordered to-do list. Per-device mode-identity status stays in the per-device
files.

`material_model_validation.md` is the one file here that is not per-device: the material
model in `beef/materials/stack.py` is shared by every example and by both engines, so a gap
in it (or a fix to it) reaches all five at once. It is the validation record for
`directional_coupler_validation.md` item 8 and `notes.md` #5/#9.

The coupler came first, from investigating why
`directional_coupler_hybrid_stub_sweep.py` never converged. It is the only
device where all four checks below have actually been run, so it doubles as the
template and the others are written as gap lists against it.

**Its closed work — items 1–7 and 9–12, with every measurement, job number and
piece of reasoning — is in `archive/dc_validation_20260727.md`.** That file was
named `dc_validation.md` until 2026-07-27, and ~76 citations across the library,
tests and examples still say *"dc_validation.md item N"*. Item numbers are stable
IDs and are never reused, so those citations resolve against the archive
directly; `directional_coupler_validation.md` carries only current status, open
to-dos and outstanding verification.

That file moved into `resolved/` on 2026-07-30, and ~23 citations across these
docs still name it bare, as `directional_coupler_validation.md`. They are
backticked prose references rather than markdown links, so they resolve by name;
read a bare citation as `resolved/directional_coupler_validation.md`.

⚠ **"No link in this directory is broken — checked" was true for prose citations
and wrong about relative links.** Moving a file down one level silently breaks
every `](README.md)` inside it, because that now resolves to
`<subdir>/README.md`. Found 2026-08-03 while filing the ring doc: two in
`resolved/ring_resonator_validation.md` (mine, same day) and **one pre-existing in
`deferred/modal_multiplexer_validation.md`**, which had been broken since its own
2026-07-31 move. All three now point at `../README.md`, and every link under this
directory is verified resolvable by walking them, not by reading. **When moving a
doc between these folders, re-check its relative links.**

**`modal_multiplexer_validation.md` moved into the new `deferred/` on
2026-07-31**, at the user's direction: verify the other devices first, come
back to the mux after. Unlike `resolved/`, this was not a claim that its work
is addressed — check 4 is flat-out blocked (`CompoundGeometry` can't back
`FDTDSolver.solve()`, `notes.md` Geometry item 8), check 3 is the worst of
any device (shipped res 40 can't even supply the requested 16-mode basis),
and the item-12/14 library fix means its own hybrid script will now raise on
its next run until `N_MODES`/`RESOLUTION` are decided jointly. `deferred/`
existed so that was legible at a glance — "not started/blocked" reads
differently from "postponed by choice," and this device was genuinely the
latter, not just the least-finished one. Same citation convention as both
moves below: the 3 actual markdown-link references (all in `running_jobs.md`)
were repointed, bare-backtick citations elsewhere were left as-is.

**Moved back out of `deferred/` to the top level on 2026-08-04**, at the
user's direction, and work resumed starting with check 1 (gap 5, engine
agreement — see the file's own status box and `running_jobs.md`). Reviewed
against the live repo before the move: nothing changed underneath it since
the 2026-08-01 snapshot above — same blocked check 4, same check-3 basis
shortfall, same shipped `N_MODES=16`/`RESOLUTION=40`. It sits at the top
level rather than back in `resolved/` for the same reason `mmi_1x4_sin300_
validation.md`/`high_contrast_grating_validation.md` do: real open gaps
remain, so neither `resolved/` (sufficiently addressed) nor `deferred/`
(postponed by choice) describes it — it is simply active, like the coupler
and both MMIs were before their own checks closed. `deferred/` is now empty;
nothing else in this tree is currently postponed. Same link-repointing
discipline as the move in: its one internal relative link
(`[README.md](../README.md)` → `[README.md](README.md)`) and the 3
markdown-link references in `running_jobs.md` were fixed, checked by walking
them, not by reading.

**2026-08-05 — a reference design landed, and it comes with a caveat worth
reading before quoting it.** The user supplied a figure from their own
(unpublished) mode-multiplexer paper — same topology (confirmed via
`archive/FECES_Legacy/configs_mmi/config_mux.json`, the legacy config this
repo's device was built from), same O-band/160 nm-Si system. But this
repo's device is the paper's pre-optimization **"Seed"** design, not the
**"Optimized"** one that was actually fabricated and measured — confirmed
by the user directly, not inferred. The paper's measured S-parameters are
therefore NOT a valid ground-truth target for this repo's own
`modal_multiplexer_fdtd.py` S-matrix. See the validation doc's own
"Reference design" section for the full account, including what the
reference DOES usefully confirm despite that gap.

**`material_model_validation.md` moved into `resolved/` on 2026-07-31**, once
its five work items (`beef.materials.stack.reparam_lorentz()` committed, the
per-band media, the fingerprint, `HybridSolver`'s index-agreement guard, and
the resolution-sensitivity sweep) all closed and its strip-test gate passed on
every stack/band combination. Same convention as the coupler's own move: the
handful of bare-backtick citations across these docs are left as-is (they
resolve by name), and the three actual markdown-link references were
repointed. ~~`dispersive=True` remains opt-in, default off, for every solver —
that is a design choice now, not an open gap.~~ **Superseded 2026-08-01:
`dispersive` now defaults to True**, so every example and both solvers track
the band. That closes the engine asymmetry — the hybrid chain's worst n_Si gap
over the O-band drops **0.007862 → 0.001349**, from over `HybridSolver`'s 2e-3
material guard to under it. A cache-fingerprint bug had to be fixed first (a
per-solver override never reached the key); see `material_model_validation.md`
step 3 and `notes.md` Solvers item 21. **Consequence for everything below:
every cache entry and all four goldens are stale again**, and the ring is the
one device not carried over (`ring_resonator_validation.md` M1–M4).

**Both MMI docs moved into `resolved/` on 2026-08-01**, once every one of the
four checks had actually been run on both devices — the 1×2's check 1 (P1)
landed that evening and was the last unrun check on either. Same convention as
the three moves above: the markdown-link references (all in `running_jobs.md`)
were repointed, and bare-backtick citations elsewhere resolve by name.

**This is the `resolved/` sense of the word — "sufficiently addressed", not
"nothing left".** Both still carry open items, and they are named here so the
move does not bury them:

- **Gap 1 / P5 — goldens.** The 1×2's three are stale (4th invalidation); the
  2×2 has never had any and its runners exist but have not been run. Both are
  blocked on **P7**, the coupler's cold regeneration pass.
- **P4 — the 2×2's checks 1 and 4 still do not reconcile.** Its candidate
  mechanism was tested and refuted the same day; the residual is structural
  (tapers or the cascade) and has no cheap first step.
- **P3 — both ship a shared `RESOLUTION = 40`** against measured needs of 120
  (1×2) and 100 (2×2). A recorded cost decision, not a bug, but it means
  neither meets check 3 the way the coupler does.
- **The 1×2's calibrated gap 5 (7.30%) is still frozen-era** and owes a
  dispersive re-run before it can be quoted beside the 10.18%.

**`fdtd_vs_EME_mode_validation.md` moved into `resolved/` on 2026-08-01** as
well. Check 2's *mechanism* is settled — the layered scheme is specified, the
guards are in the shipped code, and every device that has been audited passes
through them. Same citation convention as every move above.

Its two open to-dos are tooling, not mechanism, and both survive the move:

- **Item 6 — `interface_modes_convergence` has zero devices run through it.**
  The tool exists and is tested; nothing has been put through it. Coupler at
  res 80 is the safe first target. The archived sweep driver is **not** a
  drop-in replacement.
- **Item 1 — the laddering harness at res 60.** Not blocked by the mux
  deferral: it is mode solves only.

## The four checks

`HybridSolver` cascades FDTD sections onto EME bodies. Each check below closes a
way that cascade can be wrong while every individual piece looks fine.

**1. Engine agreement on the shared cross-section.** The two engines must agree
on how fast phase accumulates through the *same* uniform region. If they do not,
the chain's answer depends on where the cut between them is made, and no amount
of tuning the cut will make it converge. Measured by
`examples/verification/python/verify_item_6_engine_agreement.py`; ~5% off for
the coupler even with EME converged (`directional_coupler_validation.md` item 6).

That gap was confirmed by a second, independent route. The stub sweep's L₅₀
drift slope is `2·(1 − Δβ_FDTD/Δβ_EME)`, and item 7's re-run measures a ratio of
0.9457 against the 0.9536 predicted from item 6 — agreement to 0.82%, from a
200-length hybrid sweep versus a four-length phase fit. It also settles what
that sweep was for: **it never converges while the engines disagree**, because a
straight line is exactly what the formula predicts for a constant ratio ≠ 1. A
non-converging cut-point sweep is a symptom of check 1 failing, not a tuning
problem.

**Closed for the coupler, 2026-07-27 — and the fix is not what it looks like.**
The instinct is to make the two engines' mode solvers agree. That cannot work:
*both* are MPB (`meow/fde/meep.py` and `FDTDSolver.section_mode_neff()` each
call `sim.get_eigenmode()`), and they already agree to 0.57%. The gap is between
MPB's mode index and the phase meep's *time-stepper* accumulates — a planewave
eigenproblem versus a leapfrog on a staggered Yee grid, ~2.3% apart per mode. So
`HybridSolver(calibrate_body_neff=True)` measures β from the FDTD's own phase and
propagates the EME body on that instead. The stub drift collapses from +0.1085 to
−0.0007 µm/µm (167×, below the run's own scatter) and L₅₀ moves from 5.37% off
the full-device reference to **0.11%**.

Read that as **self-consistency, not proven accuracy**: the reference is meep's
discretization, and which engine's Δβ is closer to Maxwell is open as
`directional_coupler_validation.md` item 13. Adopt the flag per device only
after check 2 passes there — it is off by default.

**2. Mode identity and ordering.** The cascade contracts index *i* of the FDTD
interface against index *i* of the EME body, assuming they are the same mode.
The two solve MPB on different grids, so near-degenerate higher-order modes can
come back reordered — and if they do, the chain silently mixes modes with
nothing else in the codebase noticing. `FDTDSolver.section_mode_neff()` gives
the FDTD side's n_eff without any time-stepping, so this is cheap to check
(`directional_coupler_validation.md` item 12).

**This is not hypothetical — it failed on the 1×2 MMI.** Run across three
devices by `examples/verification/python/audit_interface_modes.py`, which puts
the coupler first as a self-check since item 12 already fixed its answer.
Results, and the two things they overturn, are in Status below. (**Fixed
2026-07-30** — `N_MODES` 12→8 on that device; the failure was in bands at
cutoff, which the reduced request no longer asks for. The *mechanism* it
demonstrated is why everything in this section exists, so it stays stated in the
present tense as a hazard.)

**2026-07-30 — the same check is now available at runtime, not just via the
audit script.** `HybridSolver(check_mode_correspondence=True)` wires
`beef.mode_identity`'s field-overlap machinery directly into `solve()`/
`solve_length_sweep()`, opt-in, so this check can run as part of an actual
device solve instead of a separate manual pass. It is best-effort and
duck-typed (see `notes.md` Solvers #14) — confirmed against real devices to
skip gracefully with a warning, not crash, whenever either adjoining
solver's own automatic padding (`FDTDSolver`'s `dpml`) reaches past a
section's straight run-out. It compares only *guided* modes, for the reason
recorded under Finding 1 below, on a shared mesh whose resolution defaults to
the coarser of the two adjoining solvers' own `res`.

Two caveats worth carrying over from building it, both recorded in full in
`notes.md` Solvers #14: the comparison mesh must genuinely resolve the
waveguide layer (an EME side *solves* on whatever mesh it is handed, so
under-resolving makes one side's modes wrong rather than merely coarse), and
radiation modes must be excluded. Both produced confident-looking but entirely
spurious "correspondence failures" before being found — worth knowing if this
check ever reports something surprising: check the *shape* of the reported
failure against those two modes before trusting it.

Running it turned up a real, previously-unrecorded bug of its own:
`CompoundGeometry` was double-applying `stub_left`/`stub_right`, which is
why the original audit script had to route around `CompoundGeometry`
entirely. **Fixed 2026-07-30** (`notes.md` Geometry item 6), so the modal
multiplexer's interfaces can now be checked through the composed geometry.

The coupler and both MMI examples use a 0.25 µm `body_stub`, thinner than
`FDTDSolver`'s automatic `dpml=0.5` query padding, which used to mean the
check could not run at those interfaces at all. **Fixed 2026-07-30**
(`resolved/fdtd_vs_EME_mode_validation.md` item 5): `HybridSolver` gained
`mode_correspondence_pad`, letting the query ask for less clearance than a
full PML footprint (safe because the mode-field solve never time-steps, so
the PML built around it is never actually exercised) — a value like 0.1 µm
fits inside `body_stub = 0.25` where `dpml = 0.5` alone would not. Fixing
the mechanism did not itself turn `check_mode_correspondence=True` on for
these three examples; that is a separate decision needing real-device
verification, still not done.

**3. Resolution convergence, per engine separately.** A single shared
`RESOLUTION` is how a 7% error survived for months: meep is converged at 40 and
meow/MPB is not, but a sweep that moves both together cannot tell you that — the
joint curve looks flat because FDTD dominates the plotted observable while EME
dominates Δβ. Sweep them independently, and sweep a quantity that is sensitive:
Δβ is a *difference* of two n_eff near 2.3, so a 0.1% error in each mode is ~10%
in their difference (`directional_coupler_validation.md` item 11).

**4. A full-device FDTD reference.** The hybrid has to be checked against a
single-domain brute-force run of the same device — that is the only check that
tests the cascade *as a whole* rather than its parts. It must carry provenance
in the file: the coupler's old reference was numerically fine but had no config,
no solver fingerprint and no geometry recorded, so its validity had to be argued
from circumstantial evidence rather than read off
(`directional_coupler_validation.md` item 5).

Power conservation is assumed throughout and is necessary but weak — every
pre-fix directional coupler result conserved power to 0.999 while carrying a
reversed wavelength axis. It rules out gross errors, not wrong physics.

## Status

Updated **2026-08-01** after the material-dispersion default flip (`dispersive`
now True — `resolved/material_model_validation.md`), which adds
`resolved_dispersive` plus the fitted Lorentz reparameterisation to the cache
fingerprint. **Every cache entry and every golden is stale across that break.**
Not everything moves, though, and that is measured **twice**, not assumed: jobs
11608940 (gap 3) and 11608941 (gap 2) both re-ran dispersive on the 1×2 and
both returned **byte-identical** numbers. The rule, from evidence: **a quantity
evaluated at a scalar wavelength is invariant across the dispersion break** —
`resolved_dispersive` is only set when `wl.min() < wl.max()`, and a medium
built at one wavelength is exact Sellmeier either way. Only the genuinely
broadband gaps (**4 and 5**) move. Same day: the 1×2's `beta_from_phases` turn-count bug found and
fixed (P2), gap 4 submitted for the first time on that device (P1), golden
provenance shipped (P6), 2×2 golden runners written (P5 half), and P4's
candidate mechanism identified with its decisive test in flight.

**Updated 2026-08-01 (evening) — all of that has now landed, and the two MMIs
are in a materially different state than the paragraph above describes.**

- **P1 CLOSED.** The 1×2's check 1 ran for the first time ever (11608875 +
  11609686): **mean Δβ ratio 0.9946**, 0.54% monotone bias. **All four checks
  are now ticked on both MMIs.**
- **The two MMIs agree on the engine bias.** 0.9946 (1×2, O-band) against
  0.9953 (2×2, C-band) — 0.07 points apart on different geometries in
  different bands, same direction. That is an **engine-level** bias, a
  stronger claim than either device supported alone.
- **P4's candidate mechanism is REFUTED** (11609680 + 11609681). The 2×2's ΔT
  sign flip survives dispersion in the same grid interval; the material change
  moves its FDTD reference by at most 0.00033. **P4 stays open, structural.**
- **Both MMIs' gap 5 re-measured dispersive.** 2×2 **0.82% → 0.76%**; 1×2
  **9.11% → 10.18%** — the 1×2 got marginally *worse* under the more correct
  material model, with its mean offset essentially unmoved, which is itself
  evidence the residual is material-model-independent.
- **Converging conclusion from two independent routes on the same day:** on
  both MMIs the hybrid-vs-FDTD residual is **structural — tapers or the
  cascade** — not calibration (1×2, P2: ~20% relative) and not the material
  model (2×2, P4: 4–8%).

**Updated 2026-08-04 — a THIRD device now confirms the engine-level bias.**
The modal multiplexer's own check 1 (its gap 5, run for the first time this
day) measured **0.9958** on its own MMI-body cross-section (O-band) —
0.05–0.12 points from both MMIs' 0.9946/0.9953, on a different material
stack (SOI_160 vs SOI_220) and a different device family entirely. Three
independent measurements now cluster in the same ~0.4–0.6% band, which is a
materially stronger claim than two: this reads as a property of MPB/meep's
discretization on an MMI-body-style cross-section specifically, largely
independent of which device the body sits inside. `calibrate_body_neff` is
not indicated on any of the three MMI-body measurements taken so far — only
the coupler (a different cross-section shape entirely, ~10%) needs it.

Remaining live: **11609658** (coupler hybrid dispersive, 12 h) + **11609659**
(dependent FDTD at L₅₀), the resubmission after 11608703 TIMED OUT at its 3 h
budget. See `running_jobs.md`.

Previously updated 2026-07-31 after jobs 11585850 (1×2 MMI hybrid vs full-device FDTD,
9.11%, uniform sign), 11585849 (2×2 MMI hybrid vs full-device FDTD, 0.82%, sign
FLIPS across the band — a different mechanism from the 1×2's), 11586098
(`resolved/fdtd_vs_EME_mode_validation.md` item 13 — the EME duplicate guard's
over-request measured inert at 16 modes/res 60, at a 3.2× cost), 11588309
(items 11 and 12 — the guard's diagnostic grid cleared at res 40 as well, and
k-seed laddering recovers a guided mode the shipped fixed seed loses) and
11588569 (2×2 MMI gap 6 — the "cleanest device audited" result re-confirmed
under both duplicate guards). A `calibrate_body_neff` test for the 1×2 was in
flight then as 11589770/11589779 — it **failed**, and the failure turned out to
be a library bug rather than a device limit; fixed 2026-08-01, see P2. Writing
the mux's first full-device FDTD reference
hit a real library gap — `CompoundGeometry` cannot back `FDTDSolver.solve()`'s
single-domain path (`notes.md` Geometry item 8) — blocked, not run.
2026-07-30
after job 11575516 (mode-spectrum resolution convergence on
the modal multiplexer) and the day's mode-identity work
(`resolved/fdtd_vs_EME_mode_validation.md` items 0–13); 2026-07-28 after job 11526547
(tutorial calibration) and the item 8 material-model change; 2026-07-27 after
job 11517245 (calibrated stub sweep); 2026-07-26 after jobs 11500412 (stub
sweep), 11500524/11500663 (`mmi_1x2` goldens) and 11500536 (interface-mode
audit).

| device | 1 engine agreement | 2 mode identity | 3 convergence | 4 FDTD reference |
|---|---|---|---|---|
| directional coupler | ✅ **closed** — β-calibrated, L₅₀ 0.11% off reference | ✅ identity, worst 0.0023 | done, engines split | done, with provenance |
| 1×2 MMI | ✅ **CLOSED 2026-08-01** — first run ever (11608875 `[0-3]` + merge 11609686), dispersive, cache-cold. Mean Δβ_FDTD/Δβ_EME = **0.9946** (0.54% bias, monotone over 10 of 11 points, no sign change), worst fit residual **3.01e-06 rad**, worst off-diagonal **1.12e-05**. Lands within 0.07 points of the 2×2's 0.9953 in a different band — **an engine-level bias, not a device artifact.** Res-40 caveat applies (P3) | ✅ **fixed** — `N_MODES` 12→8, identity, worst 0.0003 | measured, worst of the three | ✅ **RE-MEASURED DISPERSIVE 2026-08-01** (11609215 `[1-2]` + 11609687) — hybrid **10.18%** under full-device FDTD, uniform negative sign. Frozen-era was 9.11% (11585850, preserved): **the gap got marginally worse under the more correct material model**, while the mean offset barely moved (−0.025304 → −0.025413), so the residual is material-model-independent. Calibration attempt FAILED (11589770); root cause **fixed 2026-08-01** — a 2π turn-count bug in `beta_from_phases`, **not** the π sign flip this cell previously asserted (that reading is withdrawn — see `resolved/mmi_1x2_validation.md` gap 5 for the degeneracy that made it look confirmed). Calibrated rerun ✅ **succeeded 2026-08-01** (11608690) — **7.30%** (frozen-era, still owed dispersive), only ~20% better than uncalibrated, so the residual here is structural rather than calibratable |
| 2×2 MMI | ✅ **closed 2026-08-01** — Δβ_FDTD/Δβ_EME = **0.9953** (0.47%, monotone across band), T diagonal to 1.2e-05, fit residual 4e-06 rad; **21× tighter than the coupler's ~10%, so `calibrate_body_neff` is NOT indicated here** (11596845 `[0-3]` + 11596851, gap 4). Measured at res 40, which gap 2 says is not converged in absolute n_eff — both engines read at the same mesh, so the ratio partially cancels it | ✅ identity, worst 0.0007 — **re-confirmed 2026-07-31 at `N_MODES` 14 with both guards live, neither fired** (11588569, gap 6) | ✅ **closed 2026-07-31** — absolute n_eff still moves at res 100, but Δβ_01 (the quantity the self-imaging length actually depends on) converges to 0.05% by res 120 (11590785, gap 2) | ✅ **RE-MEASURED DISPERSIVE 2026-08-01** (11609211_1 + 11609680 + 11609681) — hybrid **0.76%** / 0.59% under full-device FDTD, **sign still FLIPS across band, in the same grid interval**. Frozen-era was 0.82%/0.62% (11585849, preserved). The material change moved the FDTD reference by at most 0.00033 and the hybrid by 0.00011 — **P4's candidate mechanism is REFUTED**; the residual is structural |
| modal multiplexer | ✅ **CLOSED 2026-08-04** — job 11670540 `[0-3]` + merge 11678100. Mean ratio Δβ_FDTD/Δβ_EME = **0.9958** (0.42% bias, monotone across the O-band), worst off-diagonal/fit residual both 0.0000. Within 0.05–0.12 points of both MMIs' own bias (0.9946/0.9953) — a third device in the same engine-level band. `calibrate_body_neff` NOT indicated | run, failed, root cause fixed — **one pair-member disagreement open**, now confirmed persisting to res 80 too (gap 2) | ✅ **N_MODES=14 now confidently measured, resolution-independent** (job 11673526) — clean 14/16 guided at EVERY meep resolution 20–100, a cleaner answer than the pre-fix 13/13/13/14. But **RESOLUTION does NOT plateau by 100**: 60→80 looked converged (14/14 agree) — misleadingly, since 80→100 breaks completely (0/14 agree, shifts 2–3× larger than the 60→80 step). No resolution 20–100 is demonstrably converged in absolute n_eff | ✅ **library gap FIXED** (`notes.md` Geometry item 8 part 2) AND ✅ **smoke test COMPLETED** (job 11673498, `MUX_FDTD_RES=10`) — first-ever successful run of this script. ⏳ **Real res=20 run: two attempts timed out/killed** (11688645 at 4h, 11693107 killed at 18h53m — both stuck in TE0; an archived TO campaign for this device at higher resolution did far more work in ~70-80min/iteration, see gap 3). Rescaled 1→4 nodes (96 ranks) at user's direction, resubmitted as **11707932**, running |
| ring resonator | not applicable | not applicable | shared knob | **no FDTD script** |

## Parity backlog — what the MMIs still lack versus the coupler

Opened 2026-08-01, from an audit asking what remains to bring both MMIs to
the coupler's level. **Every item here was untracked before that audit** —
the four-check table above can go all-ticks while these stay open, which is
exactly why they are listed separately. Ticking all four checks is necessary,
not sufficient.

Read the bar honestly: the coupler itself is not pristine. Its own items 8
(step 2), 13 and 15 are open, its cold regeneration pass is unfinished, and
**all four of its goldens are stale.** "Coupler level" means the four checks
are closed *and the gaps are written down*, not that nothing is left.

| # | device | item | state |
|---|---|---|---|
| **P1** | 1×2 | **Check 1 never run.** The script existed (`verify_mmi_1x2_engine_agreement.py`) but had never been submitted. | ✅ **CLOSED 2026-08-01.** Reshaped to the 2×2's array+merge pattern. First pass (11608448 `[0-3]`) completed in ~2 h but landed **frozen-material era**; its merge 11608450 was cancelled deliberately rather than let it raise `CacheMissError` under the new fingerprint. Re-run cache-cold and dispersive as **11608875 `[0-3]` + merge 11609686** (55:15 / 41:22 / 49:33 / 2:23:09 + 25:45): **mean ratio 0.9946**, 0.54% monotone bias. **All four checks now ticked on the 1×2.** Agrees with the 2×2's 0.9953 to 0.07 points across a different band, which upgrades the reading from "this device's bias" to an engine-level one. See `resolved/mmi_1x2_validation.md` gap 4 |
| **P2** | 1×2 | **`calibrate_body_neff` blocked.** The coupler ships `CALIBRATE_BODY = True`; the 1×2 defaults it *off* (`MMI1X2_CALIBRATE_BODY`) and the one attempt to enable it failed (11589770). | ✅ **CLOSED 2026-08-01.** Root cause was a 2π turn-count bug in `beta_from_phases`, *not* the π sign flip previously recorded. Fixed, and **job 11608690 ran `calibrate_body_neff=True` to completion on the real device** (1:28:14, no guard fire) — the first time it has ever succeeded here. **But calibration does not close this device's gap:** 9.11% → **7.30%**, ~20% relative, against the coupler's 5% → 0.11%. So the 1×2's residual is mostly **structural**, not the Δβ discretization error calibration corrects. Gap 5's standing "obvious next step" is now tried and answered. See `resolved/mmi_1x2_validation.md` gap 5 |
| **P3** | both | **Resolution knob still shared.** Both MMIs ship one `RESOLUTION = 40` for EME *and* FDTD. The coupler splits it — `FDTD_RESOLUTION = 40` / `EME_RESOLUTION = 80` — with the comment *"shared RESOLUTION was how a 7% error hid for months"*. Measured need: **res 120** (1×2), **res 100** (2×2). | ⚠️ **open decision, not a bug.** Both docs record the gap as a known caveat and keep 40 on cost grounds. Defensible, but it means neither MMI meets check 3 *as the coupler defines it*. Worth re-deciding explicitly |
| **P4** | 2×2 | **Checks 1 and 4 do not reconcile.** Check 1: 0.47% Δβ error, **monotone** across the band. Check 4: 0.82% ΔT, **sign-flipping** across the same band. A monotone bias cannot produce a sign-flipping error, so check 1 does not account for check 4. On the coupler the two *did* reconcile (item 7's ratio 0.9457 against 0.9536 predicted from item 6, agreeing to 0.82%) — and that mutual confirmation is much of why the coupler is trusted. | ⚠️ **STILL OPEN — the candidate mechanism was tested and REFUTED 2026-08-01** (11609211_1 + 11609680 + 11609681). The hypothesis: gap 5's ΔT zero crossing sits at band centre, and a medium frozen at band centre has error exactly zero there and antisymmetric about it, so the sign flip might be the material model rather than the device. **Re-measured with dispersion on, the sign flip survives in the same grid interval** and \|ΔT\| shrinks only 4–8%. The material change moves the FDTD reference by at most **0.00033** — 8.6% of the ΔT swing, and ~20× *less* than the 1×2's 0.0065, so the scale argument that looked like support was measuring the 1×2's much larger material sensitivity (3.4× wider relative band), not this device's. The premise fails too: the reference's frozen error is a small near-uniform offset here, **not** antisymmetric about centre. **P4 stands as originally written — the residual is structural (tapers or the cascade).** The cheapest candidate is eliminated at no extra cost, since gap 5 needed its dispersive number anyway. See `resolved/mmi_2x2_validation.md`'s P4 box |
| **P5** | 2×2 | **No goldens at all** — and there was no *runner* to generate one from, so this was never merely "not run". | 🟡 **half done 2026-08-01.** `run_mmi_2x2_{eme,fdtd,hybrid}` written and registered in `TARGETS` (7 golden targets now, up from 4), each mirroring its example's own port convention — top input excited (`input_idx=1`), bar = `output_idx=1`, cross = `output_idx=0`. Writing them needed no compute and did **not** wait on P7, so the outstanding regeneration pass now covers all three devices in one round instead of needing a second for this one. **Generating** the files still waits on P7 |
| **P6** | both | ~~**Goldens carry no provenance.**~~ ✅ **DONE 2026-08-01.** `scripts/generate_golden_values.py` now writes a `provenance` block: `beef_git`, `generated_utc`, `dispersive_default`, and the Sellmeier `si`/`sio2` indices at the run's own band edges and centre. Records the *material model* rather than a version tag, so any change underneath moves it with nothing to remember to bump — same reasoning as `resolved_material_config()`. Descriptive only: `test_golden_examples.py` still compares `transmission` alone, and a test pins that it stays that way (the git sha changes every commit). 5 new tests in `tests/test_golden_provenance.py`, meep-free. **Needed no compute and did not wait on P7** — the next regeneration will simply carry it. Existing four files stay bare until then |
| **P7** | both | **The coupler's cold regeneration pass is unfinished** — 11592462_0 TIMEOUT, still undiagnosed. Both MMIs' Order step 0 says wait for it before touching goldens, so this one job gates goldens on **both** devices. | ⏸ deferred, but it is the cheapest item on this list and the top of the critical path |

Ordering, if picked up cold: ~~**P7 → P2 → P1 → P3 → P5/P6.**~~ **P1 and P2
are now closed; the live order is P7 → P3 → P5/P6.** P4's cheap first step was
taken and came back negative, so it is no longer a candidate-mechanism
question — it needs the taper/cascade investigation it always did, and there is
still no cheap first step for that. **P7 is now unambiguously the top of the
critical path.**

Four entries need their qualifier read, not just their tick:

- **2×2 MMI, checks 1 and 4 together** — both are ticked, and read together
  they leave a real question open. Check 1 (2026-08-01) says the body's phase
  accumulation agrees between engines to 0.47%, biased **monotonically** across
  the band. Check 4 says the hybrid sits 0.82% under full-device FDTD with a
  ΔT that **flips sign** across that same band. A monotone Δβ bias cannot
  produce a sign-flipping ΔT, so check 1 does **not** account for check 4's
  residual — the cause is somewhere the two checks between them do not cover
  (taper transitions, or the cascade). Neither has been isolated. This is the
  one device where all four checks pass and the passes still do not add up;
  see `resolved/mmi_2x2_validation.md` gap 4's boxed result.
- **1×2 MMI, check 2** — the `❌ FAILS — bands 8/9` this table carried until
  2026-07-30 is closed: dropping `N_MODES` 12→8 removed the near-cutoff bands
  that were the whole problem region (`resolved/mmi_1x2_validation.md` gap 3). Finding 1
  and Finding 2 below still describe the pre-fix 12-mode measurement, and are
  labelled as such.
- **2×2 MMI, check 2** — passed, but measured 2026-07-26 against a
  configuration this device no longer uses (15 modes, and before either
  duplicate-mode guard existed). `resolved/mmi_2x2_validation.md` gap 6 is the cheap
  re-run that would make it current.
- **modal multiplexer, check 2** — the catastrophic η ≈ 0.01 mismatches were a
  duplicate-eigenmode bug and are fixed. What survives is narrower and
  sharper: at every resolution from 20 to 60 the two engines disagree about one
  member of a near-degenerate pair (meow finds both, meep's per-band ladder
  finds one). That is a driver defect, not a device finding —
  `resolved/fdtd_vs_EME_mode_validation.md` items 11–13.

`shared knob` means the example still passes one `RESOLUTION` to both
sub-solvers, so its EME half is presumptively carrying the same 7% error the
coupler had. The modal multiplexer still ships that way; the difference as of
2026-07-30 is that the number it should carry is now measured rather than
guessed. The ring's two `not applicable` entries are structural: it never
constructs `FDTDSolver`/`EMESolver` objects, so checks 1 and 2 have nothing to
compare — see `ring_resonator_validation.md`.

### Six findings that generalize

**1. Every interface basis in the repo contains at least one unguided mode.**

| device | requested | guided | worst guided \|Δn_eff\| | order |
|---|---|---|---|---|
| directional coupler | 6 | 5 | 0.0023 | ✅ identity |
| 1×2 MMI | 12 | 11 | 0.0067 | ❌ `[…,9,9,10]` |
| 2×2 MMI | 15 | 14 | 0.0007 | ✅ identity |

A mode below the cladding index is a radiation mode: its n_eff and ordering are
not well defined, and in all three devices it is the worst-agreeing column in
the basis. Each `N_MODES` should drop by one at no cost.

**The `requested` column is the state this was measured in, not the state the
examples ship in — all three have since moved** (`resolved/fdtd_vs_EME_mode_validation.md`
item 9, and `resolved/mmi_1x2_validation.md` gap 3 for the 1×2's separate fix):
coupler 6→**5**, 1×2 MMI 12→**8**, 2×2 MMI 15→**14**. The `guided`/`worst`/
`order` columns are still the correct record of what was found at the requested
count shown, and the 1×2's ❌ is fixed at its new count (worst 0.0003, identity).
Read the table as dated evidence for the finding, not as a description of
current configuration.

**2026-07-30 — this finding is stronger than "worst-agreeing".** A radiation
mode has no engine-independent *shape* either: its profile is set by the
simulation box and each engine's own discretization, so two engines' radiation
modes are essentially unrelated even when the guided basis corresponds
perfectly. Measured on the modal multiplexer with the field-overlap check:
including its three sub-cladding bands (n_eff 1.434/1.413/1.388 against a
1.4468 cladding index) drove the worst mode-to-mode overlap to **0.0000**,
while the 13 guided bands were unaffected. That is why
`HybridSolver(check_mode_correspondence=True)` filters to guided modes before
comparing — without it the check false-alarms on every device in this table.
It does **not** make the unguided columns harmless: `solve()` still contracts
them, which is exactly what "each `N_MODES` should drop by one" above is
asking for and what a correspondence check cannot fix on its own.

**2. Mode count is not what drives disagreement — confinement is.** The
expectation was that the 2×2's 15 modes would be the riskiest case and the
coupler's 6 the safest. The opposite happened: the 2×2 is the cleanest device
audited, and the failure is on the 1×2, which asks for fewer modes but is
**SOI_160**. A 160 nm core confines far more weakly than 220 nm, so the top of
its basis piles up at cutoff — its band 10 sits 0.0042 above the cladding index
and band 11 below it. That is where meep returns two bands degenerate to 2.5e-13
that meow splits by 6.19e-3. Judge a basis by where its modes sit relative to
cutoff, not by how many there are.

**3. No device's EME is converged in absolute n_eff, even at res 100.**

| device | res 40 → 100 | 40→60 | 60→80 | 80→100 |
|---|---|---|---|---|
| directional coupler | 2.3251 → 2.3521 | +0.0113 | +0.0144 | +0.0012 |
| 1×2 MMI | 2.6938 → 2.7251 | +0.0196 | +0.0034 | +0.0083 |
| 2×2 MMI | 2.8070 → 2.8286 | +0.0121 | +0.0068 | +0.0027 |

None is monotone. **Item 11's res-80 pin is not a general "converged"
resolution** and must not be copied to other devices — it was justified on Δβ,
the even/odd *difference*, which converges far faster than absolute n_eff. Per
check 3, sweep the quantity the device's result actually depends on; sweeping an
insensitive one is what hid the coupler's 7% error for months.

**2026-07-30 — and the two engines do not converge at the same resolution
either.** Job 11575516 swept meep's and meow's `res` *independently* at the
modal multiplexer's interface 0 (SOI_160, 16 modes) and counted how many of the
16 modes each engine reproduces across a resolution step:

| | 20→30 | 30→40 | 40→60 | 60→80 |
|---|---|---|---|---|
| meow (EME) | 1 | 1 | 2 | **16** |
| meep (FDTD) | 1 | 1 | 3 | *(not run)* |

The EME side plateaus between 60 and 80; the FDTD side is still moving at 60.
Two consequences. First, this is the **first converged-resolution number
measured on SOI_160** — and it is 60, not the coupler's 80, which is one more
reason the coupler's pin does not transfer. Second, `directional_coupler_
validation.md` item 11's defect is worse than "one knob for two engines that
each need a different number": the *right* numbers differ by engine **and** by
material stack, so there is no single resolution to converge onto and a shared
knob cannot be correct for both sides of any hybrid. See
`modal_multiplexer_validation.md` gap 2.

**4. A converged-looking number can still be wrong, and a resolution sweep will
not tell you.** Chasing check 1 on the coupler produced two answers that each
look converged and disagree by ~5%: meow/MPB plateaus at Δβ = 0.0894 by res
80–100, while the FDTD's Δβ sits flat at ~0.085 across res 20–50. Flat is not
the same as correct — it only means the error is not dominated by Δx. Settling
which is right needs an external reference or a third discretization, not a
finer mesh (`directional_coupler_validation.md` item 13). Worth remembering
before quoting any convergence study in this repo as proof of accuracy.

**5. Per-task outputs in an array sweep must be written atomically.** Every
sweep here ends by aggregating *all* tasks' `.npz` files, including ones a peer
is writing that instant. A plain `np.savez` is not atomic, and job 11517245 lost
**7 of 9 tasks** to `BadZipFile` because of it — the same code had passed
earlier only because those tasks finished 12 minutes apart instead of 6. Use
`beef.atomic_replace_write`, and pass `np.savez` an open **file object**: given a
suffix-less path it appends `.npz`, which defeats the rename after the
simulation work is already done. Both traps are pinned by
`tests/test_atomic_write.py`.

**6. One material model now backs both engines — and every golden is stale until
re-run.** Until 2026-07-28 the EME and the FDTD each did their own
nearest-tabulated `argmin` over meow's index tables, so at the same wavelength
they could sit on different silicon: meow's Si table has only two samples
between 1.45 and 1.65 µm, putting a 0.0063 step in the middle of the C-band.
Both now evaluate one analytic Sellmeier model (Salzberg Si, Malitson SiO2) and
agree to 4e-16. Two consequences that reach every device here:

- n_Si shifts by ~0.002, because meow's silicon table was a different, noisier
  source. **All four goldens must be regenerated**; they are `@pytest.mark.slow`,
  so the fast suite has not been checking them.
- Making the FDTD genuinely *dispersive* was tried the naive way first and
  reverted: meep subpixel-averages only the instantaneous ε, and the Sellmeier
  form puts ε∞ = 1 with the whole index in susceptibilities, so the structure
  has no contrast to smooth. On a 450x220 nm strip at res 40 that moved n_eff
  by **+0.0496** — 25x the ~0.002 the dispersion was worth. See `notes.md` #5.
  **2026-07-31 — the re-parameterised fix landed, passed its strip-test gate on
  every stack/band combination, and (job 11588997) its resolution sensitivity
  on the coupler's own body cross-section matches what production already
  ships at.** ~~It is opt-in (`dispersive=True`, default off) for every solver —
  turning it on by default anywhere is a separate, unstarted decision.~~
  **Deployed 2026-08-01: `dispersive` defaults to True everywhere**, closing the
  frozen-FDTD-vs-exact-EME asymmetry (worst n_Si gap over the O-band 0.007862 →
  0.001349, i.e. from over `HybridSolver`'s 2e-3 guard to under it). The
  remaining 0.001349 is fit error, not disagreement — the FDTD runs the
  single-Lorentzian re-parameterisation while the EME stays exact Sellmeier,
  which is why that guard needed a tolerance rather than equality. Full record
  in `material_model_validation.md`.
