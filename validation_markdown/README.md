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
`examples/verification/python/verify_item_6_engine_agreement.py`; currently
~5% off for the coupler even with EME converged, and unexplained
(`dc_validation.md` item 6).

**2. Mode identity and ordering.** The cascade contracts index *i* of the FDTD
interface against index *i* of the EME body, assuming they are the same mode.
The two solve MPB on different grids, so near-degenerate higher-order modes can
come back reordered — and if they do, the chain silently mixes modes with
nothing else in the codebase noticing. `FDTDSolver.section_mode_neff()` gives
the FDTD side's n_eff without any time-stepping, so this is cheap to check
(`dc_validation.md` item 12).

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

| device | 1 engine agreement | 2 mode identity | 3 convergence | 4 FDTD reference |
|---|---|---|---|---|
| directional coupler | measured, ~5% open | in progress | done, engines split | done, with provenance |
| 1×2 MMI | not run | not run | shared knob | script exists, unverified |
| 2×2 MMI | not run | not run | shared knob | script exists, unverified |
| modal multiplexer | not run | not run | shared knob | **no FDTD script** |
| ring resonator | not run | not run | shared knob | **no FDTD script** |

The `shared knob` entries all mean the same thing: that example still passes one
`RESOLUTION` to both sub-solvers, so its EME half is presumptively carrying the
same 7% error the coupler had. None of them have been measured.
