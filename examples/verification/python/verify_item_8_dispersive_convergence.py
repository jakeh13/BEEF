"""Does enabling dispersion move the mesh-convergence point?
(validation_markdown/material_model_validation.md, work item 4 — the cheap
half, run before committing to the expensive full-device sweep.)

What item 4 actually asks
-------------------------
`notes.md` #9 budgets a re-run of the resolution sweeps as part of turning
dispersion on, because **subpixel smoothing moves the mesh-convergence point**
and `directional_coupler_validation.md` item 11's res-80 EME pin therefore is
not automatically inherited across this change.

The literal instruction is to re-run
`batch/sweep/run_directional_coupler_*_resolution_sweep.batch`. That is 11
array tasks x 2 nodes x 24 ranks, up to 24 h each — the single most expensive
thing in this change, and it answers the question by brute force.

But **the res-80 pin was justified on Δβ**, the even/odd supermode difference
(`resolved/directional_coupler_validation.md` item 11), and Δβ is a *mode
solve* quantity: no time-stepping, no propagation, no S-matrix. So whether the
convergence point moves can be measured directly, in minutes, on the same
cross-section the pin was set on. That is what this script does.

Read it as a **decision input, not a replacement**:

    convergence point does NOT move  ->  item 11's res-80 pin is inherited;
        the full sweep is confirmation, not discovery, and can be scoped to a
        few resolutions rather than 11.
    convergence point DOES move      ->  the full sweep is required, and this
        says roughly where to look, so the swept range can be chosen instead
        of guessed.

Either way this costs minutes and the sweep costs node-days, so it runs first.

Method
------
Coupler body cross-section (the uniform straight coupling region, the same
place item 11 measured), SOI_220, at each resolution in RES_SWEEP, under both
material settings:

    non-dispersive  — today's shipped behaviour, one frozen index at band centre
    dispersive      — `FDTDSolver(dispersive=True)`, the re-parameterised media

For each: n_eff of the first two modes (the even and odd supermodes) and their
difference, from which Δβ = 2π·Δn_eff/λ. Convergence is judged on the *step
between consecutive resolutions*, the same way item 11 and
`mode_spectrum_convergence.py` judge it — an absolute value is not converged
just because it is stable to the eye.

Δn_eff is the headline rather than n_eff because that is what the coupling
length depends on, and it is the quantity item 11's pin was chosen for: the
``β/Δβ ≈ 99`` lever means a small common-mode bias in n_eff barely matters
while the difference matters enormously.
"""
from __future__ import annotations

import time

import numpy as np

from beef import DirectionalCoupler, FDTDSolver, SOI_220

WL_MIN, WL_MAX = 1.50, 1.60
WL = 0.5 * (WL_MIN + WL_MAX)
BAND = np.linspace(WL_MIN, WL_MAX, 21)

# Same range the full sweep uses, so the two are directly comparable. Trimmed
# at the top only if wall time demands it — see the printout's timing column.
RES_SWEEP = (20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 120.0)

N_MODES = 2          # even + odd supermode; the pair Δβ is defined on
HALF_SPAN = 0.05     # thin slice of the uniform coupling region
SIM_W, SIM_H = 6.0, 1.5
DPML = 0.5


def _coupler():
    """The same device the resolution sweeps and item 11 use."""
    return DirectionalCoupler(
        waveguide_width=0.45,
        gap=0.25,
        coupling_length=15.13,
        bend_length=9.0,
        waveguide_offset=3.0,
        bend_type="euler",
        waveguide_height=0.22,
    )


def _neff_pair(dc, res, dispersive):
    """(n_eff[0], n_eff[1], seconds) at the coupler body, one resolution."""
    solver = FDTDSolver(
        res=res, sim_width=SIM_W, sim_height=SIM_H, dpml=DPML,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_220, dispersive=dispersive, decay_tol=1e-11,
    )
    t0 = time.perf_counter()
    neff = solver.section_mode_neff(
        dc, BAND, -HALF_SPAN, +HALF_SPAN, N_MODES,
        monitor_center=0.0, monitor_size=SIM_W)
    dt = time.perf_counter() - t0
    # section_mode_neff returns [n_wl, n_modes] or [n_modes]; take band centre.
    arr = np.atleast_2d(np.real(neff))
    row = arr[arr.shape[0] // 2] if arr.shape[0] > 1 else arr[0]
    return float(row[0]), float(row[1]), dt


def main():
    dc = _coupler()
    print("=" * 78, flush=True)
    print("Does dispersion move the mesh-convergence point? "
          "(material model work item 4)", flush=True)
    print(f"coupler body cross-section, SOI_220, wl {WL:.3f} um, "
          f"{N_MODES} modes, band {WL_MIN}-{WL_MAX}", flush=True)
    print("=" * 78, flush=True)

    results = {}
    for dispersive in (False, True):
        label = "dispersive" if dispersive else "non-dispersive"
        print(f"\n===== {label} =====", flush=True)
        print(f"  {'res':>5} {'n_eff[0]':>11} {'n_eff[1]':>11} {'dn_eff':>11} "
              f"{'step in dn':>11} {'d(beta)':>10} {'sec':>7}", flush=True)
        print(f"  {'-'*5} {'-'*11} {'-'*11} {'-'*11} {'-'*11} {'-'*10} {'-'*7}",
              flush=True)
        rows, prev = [], None
        for res in RES_SWEEP:
            try:
                n0, n1, dt = _neff_pair(dc, res, dispersive)
            except Exception as exc:
                print(f"  {res:5.0f}  FAILED: {type(exc).__name__}: {exc}",
                      flush=True)
                continue
            dn = n0 - n1
            dbeta = 2 * np.pi * dn / WL
            step = "" if prev is None else f"{dn - prev:+.2e}"
            print(f"  {res:5.0f} {n0:11.6f} {n1:11.6f} {dn:11.6f} "
                  f"{step:>11} {dbeta:10.5f} {dt:7.1f}", flush=True)
            rows.append((res, dn))
            prev = dn
        results[dispersive] = rows

    # ---- where does each settle? ---------------------------------------
    print("\n" + "=" * 78, flush=True)
    print("CONVERGENCE — first resolution whose step falls under each "
          "threshold", flush=True)
    print("=" * 78, flush=True)
    print(f"  {'threshold':>12} {'non-dispersive':>16} {'dispersive':>16}",
          flush=True)
    print(f"  {'-'*12} {'-'*16} {'-'*16}", flush=True)
    for tol in (1e-3, 3e-4, 1e-4, 3e-5):
        cell = {}
        for dispersive in (False, True):
            hit = next((r for (r, d), (rp, dp)
                        in zip(results[dispersive][1:], results[dispersive][:-1])
                        if abs(d - dp) < tol), None)
            cell[dispersive] = f"{hit:.0f}" if hit else "not reached"
        print(f"  {tol:12.0e} {cell[False]:>16} {cell[True]:>16}", flush=True)

    # ---- direct comparison at matched resolution -----------------------
    print("\n  Δn_eff, dispersive minus non-dispersive, at matched resolution:",
          flush=True)
    d_off = dict(results[False])
    d_on = dict(results[True])
    for res in RES_SWEEP:
        if res in d_off and res in d_on:
            print(f"    res {res:5.0f}: {d_on[res] - d_off[res]:+.3e}", flush=True)

    print("\n" + "=" * 78, flush=True)
    print("HOW TO READ THIS", flush=True)
    print("=" * 78, flush=True)
    print("  Same convergence resolution in both columns  ->  item 11's res-80", flush=True)
    print("  pin survives the change. The full 11-point device sweep becomes", flush=True)
    print("  confirmation rather than discovery and can be scoped down.", flush=True)
    print("  Different convergence resolution  ->  the pin does NOT carry, the", flush=True)
    print("  full sweep is required, and the rows above say which range it", flush=True)
    print("  needs to cover.", flush=True)
    print("\n  Note this is the MODE SOLVE only. It cannot speak for the", flush=True)
    print("  cascade's S-matrix convergence, which is what the full sweep", flush=True)
    print("  measures — only for the quantity item 11's pin was set on.", flush=True)


if __name__ == "__main__":
    main()
