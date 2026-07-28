# Validation docs

One file per example device, recording what has to hold before that device's
**hybrid** result can be trusted, and what is currently missing.

    dc_validation.md                 directional coupler   — the worked example
    mmi_1x2_validation.md            1×2 MMI
    mmi_2x2_validation.md            2×2 MMI
    modal_multiplexer_validation.md  modal multiplexer
    ring_resonator_validation.md     ring resonator

`dc_validation.md` came first, from investigating why
`directional_coupler_hybrid_stub_sweep.py` never converged. It is much longer
than the others because it is the only device where all four checks below have
actually been run — so it doubles as the template, and the others are written as
gap lists against it.

## The four checks

`HybridSolver` cascades FDTD sections onto EME bodies. Each check below closes a
way that cascade can be wrong while every individual piece looks fine.

**1. Engine agreement on the shared cross-section.** The two engines must agree
on how fast phase accumulates through the *same* uniform region. If they do not,
the chain's answer depends on where the cut between them is made, and no amount
of tuning the cut will make it converge. Measured by
`examples/verification/python/verify_item_6_engine_agreement.py`; ~5% off for
the coupler even with EME converged (`dc_validation.md` item 6).

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
`dc_validation.md` item 13. Adopt the flag per device only after check 2 passes
there — it is off by default.

**2. Mode identity and ordering.** The cascade contracts index *i* of the FDTD
interface against index *i* of the EME body, assuming they are the same mode.
The two solve MPB on different grids, so near-degenerate higher-order modes can
come back reordered — and if they do, the chain silently mixes modes with
nothing else in the codebase noticing. `FDTDSolver.section_mode_neff()` gives
the FDTD side's n_eff without any time-stepping, so this is cheap to check
(`dc_validation.md` item 12).

**This is not hypothetical — it fails on the 1×2 MMI.** Run across three devices
by `examples/verification/python/audit_interface_modes.py`, which puts the
coupler first as a self-check since item 12 already fixed its answer. Results,
and the two things they overturn, are in Status below.

**3. Resolution convergence, per engine separately.** A single shared
`RESOLUTION` is how a 7% error survived for months: meep is converged at 40 and
meow/MPB is not, but a sweep that moves both together cannot tell you that — the
joint curve looks flat because FDTD dominates the plotted observable while EME
dominates Δβ. Sweep them independently, and sweep a quantity that is sensitive:
Δβ is a *difference* of two n_eff near 2.3, so a 0.1% error in each mode is ~10%
in their difference (`dc_validation.md` item 11).

**4. A full-device FDTD reference.** The hybrid has to be checked against a
single-domain brute-force run of the same device — that is the only check that
tests the cascade *as a whole* rather than its parts. It must carry provenance
in the file: the coupler's old reference was numerically fine but had no config,
no solver fingerprint and no geometry recorded, so its validity had to be argued
from circumstantial evidence rather than read off (`dc_validation.md` item 5).

Power conservation is assumed throughout and is necessary but weak — every
pre-fix directional coupler result conserved power to 0.999 while carrying a
reversed wavelength axis. It rules out gross errors, not wrong physics.

## Status

Updated 2026-07-27 after job 11517245 (calibrated stub sweep), and 2026-07-26
after jobs 11500412 (stub sweep), 11500524/11500663 (`mmi_1x2` goldens) and
11500536 (interface-mode audit).

| device | 1 engine agreement | 2 mode identity | 3 convergence | 4 FDTD reference |
|---|---|---|---|---|
| directional coupler | ✅ **closed** — β-calibrated, L₅₀ 0.11% off reference | ✅ identity, worst 0.0023 | done, engines split | done, with provenance |
| 1×2 MMI | not run | ❌ **FAILS** — bands 8/9 | measured, worst of the three | script exists, unverified |
| 2×2 MMI | not run | ✅ identity, worst 0.0007 | measured, not converged | script exists, unverified |
| modal multiplexer | not run | not run | shared knob | **no FDTD script** |
| ring resonator | not applicable | not applicable | shared knob | **no FDTD script** |

`shared knob` means the example still passes one `RESOLUTION` to both
sub-solvers, so its EME half is presumptively carrying the same 7% error the
coupler had. The ring's two `not applicable` entries are structural: it never
constructs `FDTDSolver`/`EMESolver` objects, so checks 1 and 2 have nothing to
compare — see `ring_resonator_validation.md`.

### Three findings that generalize

**1. Every interface basis in the repo contains at least one unguided mode.**

| device | requested | guided | worst guided \|Δn_eff\| | order |
|---|---|---|---|---|
| directional coupler | 6 | 5 | 0.0023 | ✅ identity |
| 1×2 MMI | 12 | 11 | 0.0067 | ❌ `[…,9,9,10]` |
| 2×2 MMI | 15 | 14 | 0.0007 | ✅ identity |

A mode below the cladding index is a radiation mode: its n_eff and ordering are
not well defined, and in all three devices it is the worst-agreeing column in
the basis. Each `N_MODES` should drop by one at no cost.

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

**4. A converged-looking number can still be wrong, and a resolution sweep will
not tell you.** Chasing check 1 on the coupler produced two answers that each
look converged and disagree by ~5%: meow/MPB plateaus at Δβ = 0.0894 by res
80–100, while the FDTD's Δβ sits flat at ~0.085 across res 20–50. Flat is not
the same as correct — it only means the error is not dominated by Δx. Settling
which is right needs an external reference or a third discretization, not a
finer mesh (`dc_validation.md` item 13). Worth remembering before quoting any
convergence study in this repo as proof of accuracy.

**5. Per-task outputs in an array sweep must be written atomically.** Every
sweep here ends by aggregating *all* tasks' `.npz` files, including ones a peer
is writing that instant. A plain `np.savez` is not atomic, and job 11517245 lost
**7 of 9 tasks** to `BadZipFile` because of it — the same code had passed
earlier only because those tasks finished 12 minutes apart instead of 6. Use
`beef.atomic_replace_write`, and pass `np.savez` an open **file object**: given a
suffix-less path it appends `.npz`, which defeats the rename after the
simulation work is already done. Both traps are pinned by
`tests/test_atomic_write.py`.
