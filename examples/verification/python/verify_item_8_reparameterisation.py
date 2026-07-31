"""Is the ε∞ re-parameterisation of the Sellmeier model actually realizable?
(validation_markdown/material_model_validation.md, work item 1.)

Why this script exists
----------------------
`material_model_validation.md` step 2 — make the FDTD medium dispersive across
the band — was implemented and **reverted**. meep subpixel-averages only the
*instantaneous* ε and leaves susceptibility terms unsmoothed per pixel. The
Sellmeier form (like meep's own library materials) puts ε∞ = 1 and carries the
whole index in susceptibilities, so meep builds a structure with
`epsilon diagonal = (1,1,1)`: no dielectric contrast at the Si/SiO2 interface,
and subpixel smoothing does nothing at all. On a 450x220 nm strip at res 40
that cost +0.049640 in n_eff against the ~0.002 the dispersion was worth --
25x the wrong way.

The prescribed remedy (item 8, `notes.md` #5/#9) is to re-parameterise so ε∞
carries the bulk of the permittivity and the susceptibilities carry only its
wavelength dependence. That was written as an *assumption*. The doc's analysis
says it holds, but the numbers behind it were never committed as runnable code
-- so this script exists to make λ₀, the per-band restriction, and the
resulting fit quality reproducible rather than asserted.

The exact split is not expressible in meep
-------------------------------------------
The algebra is clean:

    σ·λ²/(λ² - λ₀²)  =  σ  +  σ·λ₀²/(λ² - λ₀²)

and that leading constant σ is exactly what belongs in ε∞. But the remainder is
`σω²/(ω₀² - ω²)`, while meep's `LorentzianSusceptibility` is
`σω₀²/(ω₀² - ω²)` -- an ω₀² numerator where an ω² is needed, and
`DrudeSusceptibility` is not the right shape either. So ε∞ cannot be peeled off
analytically. The model has to be **re-fit**, which is what this script does.

What it measures
----------------
Fits `ε(λ) ≈ ε∞ + σ·λ²/(λ² - λ₀²)` by linear least squares at fixed λ₀, and
reports three things the doc asserts:

    1. It cannot be fit once, globally. Silicon's telecom dispersion comes from
       a UV resonance at 0.3015 µm and no far-IR Lorentzian reproduces it
       across the full range the examples span. Per-band it is comfortable;
       full-range is worse than the error being fixed.
    2. λ₀ is weakly determined and must be pinned. For λ₀ >> λ only the ratio
       σ/λ₀² matters, so a free fit walks a ridge to its bounds. Pinning λ₀
       makes the model linear in both remaining parameters -- a deterministic
       lstsq, no optimiser in the build path. **Which λ₀ to pin is NOT decided
       here**: fit quality improves monotonically with λ₀ over the usable
       range, while the smoothing split gets worse. `verify_item_8_strip_gate.py`
       is what chose the shipped value, on the opposite criterion.
    3. The number that decides it: how much of the true Si/SiO2 ε contrast the
       re-parameterised ε∞ makes available to subpixel averaging. Today it is
       0% (ε∞ = 1 for both). Anything near 100% means smoothing is back.

Read the output as: per-band max |Δn| small (≲2e-4) and ε∞ contrast ≈ 100% of
true contrast  ->  the re-parameterisation is realizable and step 2 can proceed
to its strip-test gate. Full-range max |Δn| larger than the ~0.002 dispersion
is worth  ->  a single global fit is not an option, per-band is mandatory.

This script solves no eigenproblem and imports neither meep nor meow -- it is
pure numpy over `beef.materials.stack.sellmeier_index`, and runs in
milliseconds on a login node.
"""
from __future__ import annotations

import numpy as np

from beef.materials.stack import (
    reparam_lam0, reparam_lorentz, sellmeier_index,
)

# Bands the examples actually run. The C-band is the coupler's sweep; the
# O-band is what both MMIs and the modal multiplexer use; the union is what a
# single global fit would have to cover.
BANDS = {
    "C-band  1.50-1.60 (coupler)": (1.50, 1.60),
    "O-band  1.26-1.36 (MMIs, mux)": (1.26, 1.36),
    "full    1.26-1.70 (union)": (1.26, 1.70),
}

# The shipped pole is derived per band (1.15 x band max), so there is no
# single default -- this is the C-band's, used for the comparison tables.
LAM0_DEFAULT = reparam_lam0((1.50, 1.60))
LAM0_TRIALS = (1.84, 2.0, 3.0, 10.0)
N_GRID = 201              # sample density across a band for the lstsq
CONTRAST_WL = 1.55        # where the doc quotes the true Si/SiO2 ε contrast


def fit_reparam(kind: str, lo: float, hi: float, lam0: float = LAM0_DEFAULT):
    """Thin shim onto `beef.materials.stack.reparam_lorentz`.

    Deliberately NOT a second implementation: this script's whole purpose is to
    make the shipped fit reproducible, which it cannot do by re-deriving it.
    The band is passed as its two endpoints because `reparam_lorentz` resamples
    onto its own grid internally.
    """
    return reparam_lorentz(kind, (lo, hi), lam0=lam0)


def main():
    print("=" * 72, flush=True)
    print("Sellmeier -> (eps_inf + one Lorentzian) re-parameterisation",
          flush=True)
    print(f"lam0 pinned at {LAM0_DEFAULT} um, {N_GRID} samples per band",
          flush=True)
    print("=" * 72, flush=True)

    # ---- 1. per-band vs global ------------------------------------------
    print("\n1. Fit quality per band -- can this be fit once, globally?\n",
          flush=True)
    print(f"  {'band':<32} {'Si max|dn|':>12} {'SiO2 max|dn|':>13}", flush=True)
    print(f"  {'-'*32} {'-'*12} {'-'*13}", flush=True)
    for label, (lo, hi) in BANDS.items():
        _, _, dn_si = fit_reparam("si", lo, hi)
        _, _, dn_ox = fit_reparam("sio2", lo, hi)
        print(f"  {label:<32} {dn_si:12.2e} {dn_ox:13.2e}", flush=True)
    print("\n  The full-range row is the verdict: if it exceeds the ~2e-3 the",
          flush=True)
    print("  dispersion itself is worth, one global fit is not usable and the",
          flush=True)
    print("  media must be built per run from that run's own band.", flush=True)

    # ---- 2. lam0 degeneracy ---------------------------------------------
    print("\n2. Is lam0 degenerate? (Si, C-band)\n", flush=True)
    lo, hi = BANDS["C-band  1.50-1.60 (coupler)"]
    print(f"  {'lam0 (um)':>10} {'eps_inf':>12} {'sigma':>12} {'max|dn|':>12}",
          flush=True)
    print(f"  {'-'*10} {'-'*12} {'-'*12} {'-'*12}", flush=True)
    for lam0 in LAM0_TRIALS:
        e, s, dn = fit_reparam("si", lo, hi, lam0=lam0)
        print(f"  {lam0:10.1f} {e:12.4f} {s:12.4f} {dn:12.2e}", flush=True)
    print("\n  sigma swings ~150x while max|dn| moves only ~2x -> the fit is",
          flush=True)
    print("  weakly determined in lam0: it must be PINNED, not fitted, or a free",
          flush=True)
    print("  optimiser walks the sigma/lam0^2 ridge to its bounds.", flush=True)
    print("\n  Note the trend runs OPPOSITE to how lam0 was actually chosen.",
          flush=True)
    print("  Larger lam0 fits better here, but pushes eps_inf past the true",
          flush=True)
    print("  contrast and leaves more permittivity in the UNSMOOTHED", flush=True)
    print("  susceptibility -- which is what the strip gate measures, and what",
          flush=True)
    print(f"  actually decided lam0 = {LAM0_DEFAULT}. Fit quality is not the",
          flush=True)
    print("  criterion; see verify_item_8_strip_gate.py.", flush=True)

    # ---- 3. the number that decides it ----------------------------------
    print(f"\n3. How much eps contrast becomes available to subpixel smoothing?"
          f"\n   (at {CONTRAST_WL} um, fitted over the C-band)\n", flush=True)
    e_si, s_si, dn_si = fit_reparam("si", lo, hi)
    e_ox, s_ox, dn_ox = fit_reparam("sio2", lo, hi)
    print(f"  {'':6} {'eps_inf':>12} {'sigma':>12} {'max|dn|':>12}", flush=True)
    print(f"  {'-'*6} {'-'*12} {'-'*12} {'-'*12}", flush=True)
    print(f"  {'Si':<6} {e_si:12.4f} {s_si:12.4f} {dn_si:12.2e}", flush=True)
    print(f"  {'SiO2':<6} {e_ox:12.4f} {s_ox:12.4f} {dn_ox:12.2e}", flush=True)

    eps_si = float(sellmeier_index("si", CONTRAST_WL) ** 2)
    eps_ox = float(sellmeier_index("sio2", CONTRAST_WL) ** 2)
    true_contrast = eps_si - eps_ox
    fit_contrast = e_si - e_ox
    print(f"\n  eps_inf contrast (smoothable) = {fit_contrast:.4f}", flush=True)
    print(f"  true eps contrast             = {true_contrast:.4f}", flush=True)
    print(f"  available to smoothing        = "
          f"{100 * fit_contrast / true_contrast:.1f}%   (today: 0%)", flush=True)

    resid_si = abs(s_si) * CONTRAST_WL ** 2 / abs(CONTRAST_WL ** 2 - LAM0_DEFAULT ** 2)
    print(f"  residual unsmoothed susceptibility = "
          f"{100 * resid_si / eps_si:.1f}% of eps_Si", flush=True)

    # How much of the real dispersion the fit keeps.
    swing = float(sellmeier_index("si", lo) - sellmeier_index("si", hi))
    print(f"\n  n_Si genuinely moves {swing:.5f} across this band;", flush=True)
    print(f"  fit error is {100 * dn_si / swing:.1f}% of that, so the "
          f"parameterisation keeps", flush=True)
    print(f"  ~{100 - 100 * dn_si / swing:.0f}% of the dispersion it exists to "
          f"capture.", flush=True)

    print("\n" + "=" * 72, flush=True)
    print("Gate this does NOT clear: the strip test. Fit quality says the", flush=True)
    print("model is right; only a real meep solve says whether restoring", flush=True)
    print("smoothing recovers the 0.0496 the reverted attempt lost.", flush=True)
    print("=" * 72, flush=True)


if __name__ == "__main__":
    main()
