# Validation docs

One file per example device, recording what has to hold before that device's
**hybrid** result can be trusted, and what is currently missing.

    resolved/mmi_1x2_validation.md              1×2 MMI
    resolved/mmi_2x2_validation.md              2×2 MMI
    resolved/ring_resonator_validation.md       ring resonator

    resolved/fdtd_vs_EME_mode_validation.md     cross-device — the mode-correspondence MECHANISM itself
    running_jobs.md                    snapshot — active SLURM jobs and which gap/item each is for
    mmi_1x4_sin300_validation.md       1×4 MMI on 300 nm Si3N4 — BUILD PLAN, device does not exist yet
    high_contrast_grating_validation.md HCG reflector, FDTD-only (no EME/hybrid) — convergence measured, production run in flight
    modal_multiplexer_validation.md    modal multiplexer — resumed 2026-08-04, gap list against the four checks
    archive/                           closed items, kept whole
    resolved/                          whole docs sufficiently addressed (open to-dos may remain)
    resolved/directional_coupler_validation.md  directional coupler — the worked example
    resolved/material_model_validation.md       cross-device — shared Si/SiO2 model, both engines
    deferred/                          work intentionally postponed, not blocked or closed (currently empty)

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
waveguide). Both fixed; see the doc's own top-of-file box for the full
story. Corrected, the resolution series (`res` 10→60) converges
geometrically to **R=0.9789 at `res=60`** — matching the paper's own
reported C-band minimum reflectance (97.9%) to three decimals — while x/y/z
domain-size sensitivity is negligible by comparison. The full-band (21 λ,
1.50–1.60 µm) production run at these settled settings is **submitted, job
11670282, not yet landed**. Comparison against the design's own published
spectrum (Hiesener, Pesch, Ralph, *IEEE Photonics Technology Letters*
**36**(9), 597–600, 2024) is single-wavelength-confirmed already and
full-band once that job lands — see that file's step 10 for what's still
open (the production run itself, where the non-reflected power goes, and
generalizing the domain plot).

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
| modal multiplexer | ⏳ **submitted 2026-08-04** — job 11670540 `[0-3]` + merge 11670544, uniform section of MMI 1's body (gap 5), smoke-tested (EME half real, sensible descending n_eff) before submission | run, failed, root cause fixed — **one pair-member disagreement open**; guard confirmed firing on the *production* path 2026-07-31 (11585847) | **measured** — EME converged res 60, FDTD not at 60, and as of 2026-07-31 shipped res 40 **cannot supply** the requested 16-mode basis — at most 14 distinct guided modes are reachable there, and today's fixed k-seed finds only 13 (gap 2) | ⚠️ **BLOCKED 2026-07-31** — script written, cannot run: `CompoundGeometry` can't back `FDTDSolver.solve()` (`notes.md` Geometry item 8) |
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
