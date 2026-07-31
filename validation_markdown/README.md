# Validation docs

One file per example device, recording what has to hold before that device's
**hybrid** result can be trusted, and what is currently missing.

    mmi_1x2_validation.md              1×2 MMI
    mmi_2x2_validation.md              2×2 MMI
    modal_multiplexer_validation.md    modal multiplexer
    ring_resonator_validation.md       ring resonator

    fdtd_vs_EME_mode_validation.md     cross-device — the mode-correspondence MECHANISM itself
    running_jobs.md                    snapshot — active SLURM jobs and which gap/item each is for
    archive/                           closed items, kept whole
    resolved/                          whole docs sufficiently addressed (open to-dos may remain)
    resolved/directional_coupler_validation.md  directional coupler — the worked example
    resolved/material_model_validation.md       cross-device — shared Si/SiO2 model, both engines

Check 2 below (mode identity) has grown its own cross-device file:
`fdtd_vs_EME_mode_validation.md` covers the *mechanism* — what must hold for any
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
backticked prose references rather than markdown links (no link in this
directory is broken — checked), so they resolve by name; read a bare citation as
`resolved/directional_coupler_validation.md`.

**`material_model_validation.md` moved into `resolved/` on 2026-07-31**, once
its five work items (`beef.materials.stack.reparam_lorentz()` committed, the
per-band media, the fingerprint, `HybridSolver`'s index-agreement guard, and
the resolution-sensitivity sweep) all closed and its strip-test gate passed on
every stack/band combination. Same convention as the coupler's own move: the
handful of bare-backtick citations across these docs are left as-is (they
resolve by name), and the three actual markdown-link references were
repointed. `dispersive=True` remains opt-in, default off, for every solver —
that is a design choice now, not an open gap.

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
(`fdtd_vs_EME_mode_validation.md` item 5): `HybridSolver` gained
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

Updated 2026-07-31 after jobs 11585850 (1×2 MMI hybrid vs full-device FDTD,
9.11%, uniform sign), 11585849 (2×2 MMI hybrid vs full-device FDTD, 0.82%, sign
FLIPS across the band — a different mechanism from the 1×2's), 11586098
(`fdtd_vs_EME_mode_validation.md` item 13 — the EME duplicate guard's
over-request measured inert at 16 modes/res 60, at a 3.2× cost), 11588309
(items 11 and 12 — the guard's diagnostic grid cleared at res 40 as well, and
k-seed laddering recovers a guided mode the shipped fixed seed loses) and
11588569 (2×2 MMI gap 6 — the "cleanest device audited" result re-confirmed
under both duplicate guards). A `calibrate_body_neff` test for the 1×2 is in
flight (11589770/11589779); writing the mux's first full-device FDTD reference
hit a real library gap — `CompoundGeometry` cannot back `FDTDSolver.solve()`'s
single-domain path (`notes.md` Geometry item 8) — blocked, not run.
2026-07-30
after job 11575516 (mode-spectrum resolution convergence on
the modal multiplexer) and the day's mode-identity work
(`fdtd_vs_EME_mode_validation.md` items 0–13); 2026-07-28 after job 11526547
(tutorial calibration) and the item 8 material-model change; 2026-07-27 after
job 11517245 (calibrated stub sweep); 2026-07-26 after jobs 11500412 (stub
sweep), 11500524/11500663 (`mmi_1x2` goldens) and 11500536 (interface-mode
audit).

| device | 1 engine agreement | 2 mode identity | 3 convergence | 4 FDTD reference |
|---|---|---|---|---|
| directional coupler | ✅ **closed** — β-calibrated, L₅₀ 0.11% off reference | ✅ identity, worst 0.0023 | done, engines split | done, with provenance |
| 1×2 MMI | not run | ✅ **fixed** — `N_MODES` 12→8, identity, worst 0.0003 | measured, worst of the three | ✅ **measured 2026-07-31** — hybrid 9.11% under full-device FDTD, uniform sign, uncalibrated (11585850); calibration test in flight (11589770/11589779) |
| 2×2 MMI | not run | ✅ identity, worst 0.0007 — **re-confirmed 2026-07-31 at `N_MODES` 14 with both guards live, neither fired** (11588569, gap 6) | measured, not converged | ✅ **measured 2026-07-31** — hybrid 0.82% under full-device FDTD, sign FLIPS across band (11585849) |
| modal multiplexer | not run | run, failed, root cause fixed — **one pair-member disagreement open**; guard confirmed firing on the *production* path 2026-07-31 (11585847) | **measured** — EME converged res 60, FDTD not at 60, and as of 2026-07-31 shipped res 40 **cannot supply** the requested 16-mode basis — at most 14 distinct guided modes are reachable there, and today's fixed k-seed finds only 13 (gap 2) | ⚠️ **BLOCKED 2026-07-31** — script written, cannot run: `CompoundGeometry` can't back `FDTDSolver.solve()` (`notes.md` Geometry item 8) |
| ring resonator | not applicable | not applicable | shared knob | **no FDTD script** |

Three entries need their qualifier read, not just their tick:

- **1×2 MMI, check 2** — the `❌ FAILS — bands 8/9` this table carried until
  2026-07-30 is closed: dropping `N_MODES` 12→8 removed the near-cutoff bands
  that were the whole problem region (`mmi_1x2_validation.md` gap 3). Finding 1
  and Finding 2 below still describe the pre-fix 12-mode measurement, and are
  labelled as such.
- **2×2 MMI, check 2** — passed, but measured 2026-07-26 against a
  configuration this device no longer uses (15 modes, and before either
  duplicate-mode guard existed). `mmi_2x2_validation.md` gap 6 is the cheap
  re-run that would make it current.
- **modal multiplexer, check 2** — the catastrophic η ≈ 0.01 mismatches were a
  duplicate-eigenmode bug and are fixed. What survives is narrower and
  sharper: at every resolution from 20 to 60 the two engines disagree about one
  member of a near-degenerate pair (meow finds both, meep's per-band ladder
  finds one). That is a driver defect, not a device finding —
  `fdtd_vs_EME_mode_validation.md` items 11–13.

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
examples ship in — all three have since moved** (`fdtd_vs_EME_mode_validation.md`
item 9, and `mmi_1x2_validation.md` gap 3 for the 1×2's separate fix):
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
  ships at.** It is opt-in (`dispersive=True`, default off) for every solver —
  turning it on by default anywhere is a separate, unstarted decision. Full
  record in `material_model_validation.md`.
