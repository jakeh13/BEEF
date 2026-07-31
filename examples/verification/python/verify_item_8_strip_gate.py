"""The strip test: does the re-parameterisation get subpixel smoothing back?
(validation_markdown/material_model_validation.md, Verification — the gate.)

What this decides
-----------------
Step 2 (dispersive FDTD media) was implemented once and reverted. The Sellmeier
form puts eps_inf = 1 with the whole index in susceptibilities; meep smooths
only the *instantaneous* eps, so the structure came out with
`epsilon diagonal = (1,1,1)`, no dielectric contrast for subpixel averaging to
act on, and a 450x220 nm strip at res 40 moved n_eff 2.313398 -> 2.363038,
**+0.0496** -- 25x the ~0.002 the dispersion was worth.

`reparam_lorentz()` re-fits into `eps_inf + one Lorentzian` so eps_inf carries
~104% of the Si/SiO2 contrast. The fit itself is already verified as a *model*
(`verify_item_8_reparameterisation.py`: tracks Sellmeier to 1.96e-4). What no
amount of fitting can tell us is whether meep then actually smooths it. Only a
real solve does, and this is that solve.

**Acceptance: |dn_eff| <= 0.002 against the non-dispersive baseline** -- below
what the dispersion is worth, per the validation doc. Items 4 (resolution
sweeps) and 5's hard assert are gated on this passing.

Three arms, because two cannot prove the harness works
------------------------------------------------------
    A. non-dispersive   mp.Medium(index=n)      -- today's shipped behaviour,
                                                   the baseline (expect ~2.3134)
    B. re-parameterised eps_inf + 1 Lorentzian  -- what is under test
    C. legacy Sellmeier eps_inf=1 + 3 Lorentzians -- the REVERTED attempt

C is not redundant. If B comes back clean but the harness simply fails to
reproduce the original problem, B's result means nothing. C must reproduce the
known +0.0496 failure; only then is B's number evidence about the fix rather
than about the test. Run A and C first for exactly that reason -- if C does not
fail, stop reading B.

Cheap on purpose: a 450x220 nm strip in a 2x2 um cell at res 40 is a ~80x80
cross-section and one eigensolve per arm, seconds each. No time-stepping.
"""
from __future__ import annotations

import numpy as np

from beef.materials.stack import (
    meep_medium, meep_medium_reparam, sellmeier_index,
)

RES = 40.0                # the baseline's resolution
W_SI = 0.45               # 450 nm strip width
CELL_Y, CELL_Z = 2.0, 2.0 # cladding span around it
X_SPAN = 0.2              # short run of uniform waveguide; no propagation here
BASELINE = 2.313398       # recorded non-dispersive n_eff (the doc's number)
REVERTED = 2.363038       # what the eps_inf=1 attempt scored
TOL = 0.002               # acceptance, per the validation doc

# (label, core height um, band) -- every (stack, band) combination the examples
# actually run. The original gate covered only the first row; SOI_160 is what
# mmi_1x2 and the modal multiplexer use, and a 160 nm core confines far more
# weakly, so proportionally more of the mode sits at the Si/SiO2 interface --
# which is exactly where subpixel smoothing acts and therefore exactly what
# this gate measures. A 220 nm pass does not transfer to it for free
# (cf. mmi_1x2 gap 3 and modal_multiplexer gap 2, where SOI_160 behaved
# differently from SOI_220 on unrelated axes).
CASES = [
    ("SOI_220, C-band (coupler, mmi_2x2, ring)", 0.22, (1.50, 1.60)),
    ("SOI_220, O-band", 0.22, (1.26, 1.36)),
    ("SOI_160, C-band", 0.16, (1.50, 1.60)),
    ("SOI_160, O-band (mmi_1x2, modal_mux)", 0.16, (1.26, 1.36)),
]


def _solve_neff(si_med, sio2_med, label, h_si=0.22, wl=1.55):
    """Fundamental-mode n_eff of the strip under one pair of media.

    Subpixel smoothing is left at meep's default (`eps_averaging=True`) --
    that is the entire quantity under test, so it must not be disabled.
    """
    import meep as mp

    cell = mp.Vector3(X_SPAN, CELL_Y, CELL_Z)
    geometry = [mp.Block(size=mp.Vector3(mp.inf, W_SI, h_si),
                         center=mp.Vector3(), material=si_med)]
    sim = mp.Simulation(
        cell_size=cell,
        geometry=geometry,
        default_material=sio2_med,
        resolution=int(RES),
        dimensions=3,
        eps_averaging=True,
    )
    sim.init_sim()
    f = 1.0 / wl
    where = mp.Volume(center=mp.Vector3(), size=mp.Vector3(0, CELL_Y, CELL_Z))
    em = sim.get_eigenmode(f, mp.X, where, 1,
                           mp.Vector3(f * 2.4, 0, 0), parity=mp.NO_PARITY)
    neff = float(np.real(em.k.x)) / f
    sim.reset_meep()
    # The INSTANTANEOUS permittivity — the only part meep subpixel-averages.
    # Deliberately not `get_epsilon_point(pt, f)`, which evaluates the full
    # frequency-dependent eps (susceptibilities included) and therefore reads
    # ~12.09 for every arm, hiding the very mechanism under test.
    eps_inst = float(np.real(si_med.epsilon_diag.x))
    print(f"  {label:<34} n_eff = {neff:.6f}   "
          f"eps_inf(Si, smoothed) = {eps_inst:.4f}", flush=True)
    return neff, eps_inst


def _run_case(label, h_si, band_lo, band_hi):
    """One (core height, band) combination. Returns (delta_B, delta_C, ok)."""
    import meep as mp

    wl = 0.5 * (band_lo + band_hi)
    band = np.linspace(band_lo, band_hi, 21)
    n_si = float(sellmeier_index("si", wl))
    n_ox = float(sellmeier_index("sio2", wl))

    print("\n" + "=" * 74, flush=True)
    print(f"CASE: {label}", flush=True)
    print(f"  {W_SI*1000:.0f} x {h_si*1000:.0f} nm Si in SiO2, res {RES:.0f}, "
          f"band {band_lo}-{band_hi} um, solved at {wl:.3f}", flush=True)
    print(f"  exact n_Si = {n_si:.6f}, n_SiO2 = {n_ox:.6f}, "
          f"true eps contrast = {n_si**2 - n_ox**2:.4f}", flush=True)
    print("=" * 74, flush=True)

    n_a, eps_a = _solve_neff(mp.Medium(index=n_si), mp.Medium(index=n_ox),
                             "A non-dispersive", h_si, wl)
    n_c, eps_c = _solve_neff(meep_medium("si"), meep_medium("sio2"),
                             "C legacy eps_inf=1 (reverted)", h_si, wl)
    n_b, eps_b = _solve_neff(meep_medium_reparam("si", band),
                             meep_medium_reparam("sio2", band),
                             "B re-parameterised", h_si, wl)

    d_b, d_c = n_b - n_a, n_c - n_a
    harness_ok = abs(d_c) > 10 * TOL
    passed = harness_ok and abs(d_b) <= TOL
    print(f"  -> C reproduces the failure: "
          f"{'YES' if harness_ok else 'NO -- result inconclusive'} "
          f"({d_c:+.6f})", flush=True)
    print(f"  -> B dn_eff vs A: {d_b:+.6f}   "
          f"{'PASS' if passed else ('INCONCLUSIVE' if not harness_ok else 'FAIL')}"
          f"  (tol {TOL})", flush=True)
    print(f"  -> eps_inf(Si) smoothed: A {eps_a:.4f} | C {eps_c:.4f} "
          f"| B {eps_b:.4f}", flush=True)
    return d_b, d_c, passed


def main():
    print("=" * 74, flush=True)
    print("Strip gate -- every (material stack, band) the examples run",
          flush=True)
    print(f"acceptance: |dn_eff| <= {TOL} vs the non-dispersive arm, and arm C",
          flush=True)
    print("must reproduce the reverted attempt's failure or the case is void.",
          flush=True)
    print("=" * 74, flush=True)

    results = {}
    for label, h_si, (lo, hi) in CASES:
        results[label] = _run_case(label, h_si, lo, hi)

    print("\n" + "=" * 74, flush=True)
    print("SUMMARY", flush=True)
    print("=" * 74, flush=True)
    print(f"  {'case':<40} {'dn_eff(B)':>11} {'dn_eff(C)':>11} {'gate':>8}",
          flush=True)
    print(f"  {'-'*40} {'-'*11} {'-'*11} {'-'*8}", flush=True)
    for label, (d_b, d_c, ok) in results.items():
        print(f"  {label:<40} {d_b:+11.6f} {d_c:+11.6f} "
              f"{'PASS' if ok else 'FAIL':>8}", flush=True)

    every = all(ok for _, _, ok in results.values())
    print("\n  OVERALL: " + ("PASS -- every stack/band combination the examples "
                             "use clears the\n  gate, so dispersion is safe to "
                             "enable across the device set."
                             if every else
                             "FAIL -- at least one combination does not clear "
                             "the gate.\n  Dispersion must not be enabled by "
                             "default for the affected devices."), flush=True)
    print("=" * 74, flush=True)


if __name__ == "__main__":
    main()
