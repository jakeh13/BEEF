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

WL = 1.55                 # where the baseline was measured
BAND = (1.50, 1.60)       # C-band, the fit window (see the reparam script)
RES = 40.0                # the baseline's resolution
W_SI, H_SI = 0.45, 0.22   # 450 x 220 nm strip
CELL_Y, CELL_Z = 2.0, 2.0 # cladding span around it
X_SPAN = 0.2              # short run of uniform waveguide; no propagation here
BASELINE = 2.313398       # recorded non-dispersive n_eff (the doc's number)
REVERTED = 2.363038       # what the eps_inf=1 attempt scored
TOL = 0.002               # acceptance, per the validation doc


def _solve_neff(si_med, sio2_med, label):
    """Fundamental-mode n_eff of the strip under one pair of media.

    Subpixel smoothing is left at meep's default (`eps_averaging=True`) --
    that is the entire quantity under test, so it must not be disabled.
    """
    import meep as mp

    cell = mp.Vector3(X_SPAN, CELL_Y, CELL_Z)
    geometry = [mp.Block(size=mp.Vector3(mp.inf, W_SI, H_SI),
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
    f = 1.0 / WL
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


def main():
    import meep as mp

    n_si = float(sellmeier_index("si", WL))
    n_ox = float(sellmeier_index("sio2", WL))
    print("=" * 74, flush=True)
    print(f"Strip gate: {W_SI*1000:.0f} x {H_SI*1000:.0f} nm Si in SiO2, "
          f"res {RES:.0f}, wl {WL} um", flush=True)
    print(f"exact n_Si({WL}) = {n_si:.6f}   n_SiO2 = {n_ox:.6f}   "
          f"true eps contrast = {n_si**2 - n_ox**2:.4f}", flush=True)
    print(f"acceptance: |dn_eff| <= {TOL} vs arm A", flush=True)
    print("=" * 74, flush=True)

    band = np.linspace(*BAND, 21)

    print("\n--- A. non-dispersive (today's shipped behaviour) ---", flush=True)
    n_a, eps_a = _solve_neff(mp.Medium(index=n_si), mp.Medium(index=n_ox),
                             "mp.Medium(index=n)")

    print("\n--- C. legacy Sellmeier, eps_inf=1 (the REVERTED attempt) ---",
          flush=True)
    print("  If this does not reproduce the failure, the harness is not "
          "measuring\n  what the original attempt measured and B proves "
          "nothing.", flush=True)
    n_c, eps_c = _solve_neff(meep_medium("si"), meep_medium("sio2"),
                             "eps_inf=1 + Lorentzians")

    print("\n--- B. re-parameterised (UNDER TEST) ---", flush=True)
    n_b, eps_b = _solve_neff(meep_medium_reparam("si", band),
                             meep_medium_reparam("sio2", band),
                             "eps_inf + 1 Lorentzian")

    # ---- verdict ------------------------------------------------------
    d_b, d_c = n_b - n_a, n_c - n_a
    print("\n" + "=" * 74, flush=True)
    print("RESULTS", flush=True)
    print("=" * 74, flush=True)
    print(f"  A non-dispersive : {n_a:.6f}   (doc baseline {BASELINE:.6f}, "
          f"delta {n_a - BASELINE:+.6f})", flush=True)
    print(f"  C legacy eps_inf=1: {n_c:.6f}   dn vs A = {d_c:+.6f}   "
          f"(doc recorded {REVERTED - BASELINE:+.6f})", flush=True)
    print(f"  B re-parameterised: {n_b:.6f}   dn vs A = {d_b:+.6f}", flush=True)

    print("\n  Harness validity (does C reproduce the known failure?)",
          flush=True)
    if abs(d_c) > 10 * TOL:
        print(f"    YES -- C is off by {d_c:+.4f}, the smoothing loss is real "
              "here.\n    B's number is therefore evidence about the fix.",
              flush=True)
        harness_ok = True
    else:
        print(f"    NO -- C only moved {d_c:+.4f}. This harness does not "
              "reproduce the\n    original failure, so B tells us nothing. "
              "Do NOT pass the gate on it.", flush=True)
        harness_ok = False

    print(f"\n  eps_inf(Si) actually smoothed: A {eps_a:.4f} | C {eps_c:.4f} "
          f"| B {eps_b:.4f}", flush=True)
    print("    C at 1.0 IS the mechanism: nothing for subpixel averaging to\n"
          "    act on. B restores it, which is the whole point of the re-fit.",
          flush=True)

    print("\n  GATE", flush=True)
    if not harness_ok:
        print("    INCONCLUSIVE -- fix the harness before reading the gate.",
              flush=True)
    elif abs(d_b) <= TOL:
        print(f"    PASS -- |dn_eff| = {abs(d_b):.6f} <= {TOL}. The "
              "re-parameterisation keeps\n    subpixel smoothing. Work items 4 "
              "(resolution sweeps) and 5's hard\n    assert are unblocked.",
              flush=True)
    else:
        print(f"    FAIL -- |dn_eff| = {abs(d_b):.6f} > {TOL}. Smoothing is "
              "still degraded;\n    do not enable dispersion by default.",
              flush=True)
    print("=" * 74, flush=True)


if __name__ == "__main__":
    main()
