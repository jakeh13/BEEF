# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Georgia Tech Research Corporation
"""MaterialStack: single source of truth for waveguide core/cladding layer
definitions.

Replaces the per-script layer-extrusion and meow-to-meep material-conversion
helpers that were duplicated across every master and sweep script.

Originally Si-core/SiO2-clad only; generalized 2026-08-03 to any core/clad
pair in `_SELLMEIER` via `core_kind`/`clad_kind`. `SOI_160`/
`SOI_220` are unaffected (still Si/SiO2, still the codebase's default); `si_*`
method names on `MaterialStack` are now aliases for the generalized `core_*`
names, kept so no existing caller had to change.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

# ---------------------------------------------------------------------------
# One dispersion model per material, shared by BOTH engines.
#
# Sellmeier coefficients as (eps_inf, ((resonance_wl_um, strength), ...)), so
#
#     eps(l) = eps_inf + SUM_i  sigma_i * l^2 / (l^2 - lambda_i^2)
#
# These are lifted from meep's own materials library, converting its Lorentzian
# form (frequency = 1/lambda_i, gamma = 0) into wavelength terms.
#
# ⚠ "verified to reproduce `meep.materials.Si` / `fused_quartz` to 4e-16 over
# 1.4-1.7 um" is what this comment claimed until 2026-08-03, and it was NEVER
# TESTED. It is now, and it is **false for silicon**: term 2's sigma is written
# here to 7 significant figures (0.003043475) where meep carries 8
# (0.0030434748), which puts n_Si up to **8.35e-11** off the library over
# 1.4-1.7 um. SiO2 and Si3N4 do hold at machine precision (<1e-14, measured).
#
# NOT corrected, deliberately: 8.35e-11 is ~7 orders below anything modelled
# here, but the fix moves `resolved_si_index`/`resolved_si_reparam` in their 9th
# decimal and would restale every dispersive cache entry and all four goldens.
# Do it during a scheduled regeneration, not on its own. Guarded either way by
# `tests/test_material_reparam.py::TestSellmeierMatchesMeepLibrary`, which pins
# the current drift so a real coefficient change still trips.
#
# Holding the numbers here rather than importing them means meow
# and meep are driven from one source instead of two, which is the whole point:
# before 2026-07-28 the EME interpolated meow's sampled tables while the FDTD
# froze a single snapped index, and nothing kept the two consistent.
#
#   Si    - Salzberg & Villa, JOSA 47, 244 (1957)
#   SiO2  - Malitson, JOSA 55, 1205 (1965), fused silica
#   Si3N4 - K. Luke et al., Optics Letters 40, 4823-26 (2015), via meep's
#           `Si3N4_NIR`; see the block comment below before touching it.
#
# meow's own `silicon_oxide` table *is* Malitson: they agree to 0.00000 across
# 1.50-1.60 um. Its `silicon` table is a different, noisier source sitting
# ~0.002 above Salzberg, so adopting this model shifts n_Si by that much —
# a deliberate change, not a refactor. It buys a smooth analytic n(lambda)
# where the table had only two samples between 1.45 and 1.65 um.
#
# ---------------------------------------------------------------------------
# SILICON NITRIDE, added 2026-08-03 for the 300 nm SiN platform.
#
# ⚠ THREE nitride models ship with meep and ONLY `Si3N4_NIR` is usable here.
# Reaching for the one called plainly `Si3N4` because of its name is a silent
# error, so the alternatives are recorded rather than left to be rediscovered:
#
#   meep model      valid wl (um)    gamma      verdict
#   Si3N4           0.230 - 0.830    1.4413     NO - visible-range fit, AND
#                                               lossy: a nonzero gamma is not
#                                               representable in this lossless
#                                               (eps_inf, ((lam_i, sigma_i),..))
#                                               form at all. Raises
#                                               out-of-range at 1.31/1.55 um.
#   Si3N4_VISNIR    0.207 - 1.240    0          NO - upper edge is below the
#                                               C-band; raises at 1.31/1.55 um.
#   Si3N4_NIR       0.310 - 5.504    0          YES - covers O- and C-band and
#                                               all of _MODEL_RANGE_UM below.
#
# `Si3N4_NIR` has eps_inf = 1 and two gamma=0 Lorentzians, the same convention
# as the two entries above, so this is a direct transcription rather than a
# re-fit: meep writes the poles literally as 0.1353406 um and 1239.842 um.
# MEASURED: reproduces `meep.materials.Si3N4_NIR` to max|dn| = 4.44e-16 across
# _MODEL_RANGE_UM, i.e. machine precision — the same standard the Si/SiO2
# entries are held to. Guarded by
# `tests/test_material_reparam.py::TestSellmeierMatchesMeepLibrary`.
#
# The second term is a far-IR pole with an enormous strength (lam_2 = 1239.842
# um, sigma_2 = 40314). That is how Luke parameterises it; it is numerically
# harmless at telecom and is why this entry has two terms where the other two
# have three. Do NOT "simplify" it away.
#
# n(1.31 um) = 2.003130, n(1.55 um) = 1.996280 — roughly half silicon's index,
# so a SiN/SiO2 stack has ~5.3x less eps contrast than Si/SiO2. Consequences
# for confinement, mode size and the subpixel-averaging gate are NOT settled by
# adding this entry; see the validation doc's part A7.
# ---------------------------------------------------------------------------
_SELLMEIER: Dict[str, Any] = {
    "si": (1.0, ((0.301516485, 10.6684293),
                 (1.13475115, 0.003043475),
                 (1104.0, 1.54133408))),
    "sio2": (1.0, ((0.0684043, 0.6961663),
                   (0.1162414, 0.4079426),
                   (9.896161, 0.8974794))),
    "si3n4": (1.0, ((0.1353406, 3.0249),
                    (1239.842, 40314.0))),
}

# Range the meow-side sampled copies span, and the range these models are used
# over. Deliberately wider than meep's declared validity for Si (1.36-11 um):
# examples here run the O-band down to 1.26 um, and Salzberg's stated floor is
# the *measurement* range of the 1957 paper, not a pole — the nearest resonance
# in the fit sits at 0.3015 um, so evaluating at 1.26 um is a mild
# extrapolation of a smooth function, 0.001 from meow's own table there
# (against a 0.0028 disagreement in the C-band where both are "valid").
# `meep_medium()` therefore builds the susceptibilities directly rather than
# importing `meep.materials.Si`, whose `valid_freq_range` would raise at
# 1.26 um. NOTE it is currently unused by the solvers -- see `si_medium()` for
# the measurement that sent the FDTD back to a non-dispersive Medium, and what
# has to change before dispersion can be turned on.
_MODEL_RANGE_UM = (1.10, 2.10)
_MEOW_SAMPLES = 2001          # 0.5 nm spacing; interpolation error ~1e-12

_material_cache: Dict[str, Any] = {}


def sellmeier_index(kind: str, wl) -> np.ndarray:
    """n(wl) for "si" or "sio2" from the shared model above. wl in µm."""
    eps_inf, terms = _SELLMEIER[kind]
    wl = np.asarray(wl, dtype=float)
    eps = eps_inf + sum(s * wl**2 / (wl**2 - li**2) for li, s in terms)
    return np.sqrt(eps)


def meow_material(kind: str):
    """meow material for `kind`, sampled densely from the shared model.

    A `SampledMaterial` rather than a new analytic `MaterialBase` subclass
    because meow types `Structure3D.material` as a closed union
    (`IndexMaterial | SampledMaterial | TidyMaterial`) — a custom class would
    be rejected or silently coerced. At 0.5 nm spacing the linear
    interpolation meow applies on top is accurate to ~1e-12, far below every
    effect being modelled, so nothing is lost by going through samples.
    """
    key = f"meow:{kind}"
    if key not in _material_cache:
        import meow as mw

        wls = np.linspace(*_MODEL_RANGE_UM, _MEOW_SAMPLES)
        _material_cache[key] = mw.SampledMaterial(
            name=f"beef_{kind}",
            n=sellmeier_index(kind, wls).astype(complex),
            params={"wl": wls},
        )
    return _material_cache[key]


# Pole for the re-parameterised single-Lorentzian fit, as a MULTIPLE of the
# band's upper edge rather than a fixed wavelength.
#
# It cannot be *fitted*: for lam0 >> lambda only the ratio sigma/lam0^2 is
# determined, so a free fit walks that ridge to its bounds. Pinning it makes the
# model linear in the two remaining parameters -- a deterministic lstsq, no
# optimiser in the medium-build path.
#
# **Scaled to the band, because no single global value works.** lam0 sets how
# the permittivity splits between the SMOOTHED eps_inf and the UNSMOOTHED
# susceptibility, and only that split decides whether subpixel averaging still
# works. Larger lam0 fits better but leaves more unsmoothed, biasing n_eff low.
# The conflict is structural, not a preference: SOI_160 over the O-band needs
# lam0 <= 1.60 um to clear the strip gate, while the C-band's own
# near-resonance floor is 1.10 * 1.60 = 1.76 um. Incompatible. Deriving lam0
# from the band removes the conflict, and by construction always clears
# `_REPARAM_MIN_LAM0_RATIO`.
#
# 1.15 measured on every (stack, band) combination the examples use, against a
# 0.002 tolerance -- worst case -0.00141, ~30% margin:
#
#   factor   SOI_220 C   SOI_220 O   SOI_160 C   SOI_160 O
#    1.15     -0.00033    -0.00010    -0.00109    -0.00141   all PASS
#    1.25     -0.00069    -0.00060    -0.00172    -0.00239   SOI_160 O FAILS
#
# SOI_160 is the binding constraint: a 160 nm core confines weakly, so more of
# the mode sits at the Si/SiO2 interface where smoothing acts -- the reverted
# attempt's failure is ~2x larger there (+0.108 vs +0.050). A tolerance met on
# SOI_220 alone does not transfer. Measured by
# `examples/verification/python/verify_item_8_strip_gate.py`.
_REPARAM_LAM0_FACTOR = 1.15
_REPARAM_SAMPLES = 201
# The fitted term blows up as the band approaches lam0, so refuse rather than
# silently return a near-resonant model. The derived lam0 above always clears
# this; the check exists for callers passing lam0 explicitly.
_REPARAM_MIN_LAM0_RATIO = 1.10


def reparam_lam0(wl_band) -> float:
    """The pinned Lorentzian pole (µm) this band's fit should use.

    Derived rather than constant so one rule covers every band — see
    `_REPARAM_LAM0_FACTOR` for why a single global value cannot work.
    """
    wl = np.atleast_1d(np.asarray(wl_band, dtype=float))
    return float(wl.max()) * _REPARAM_LAM0_FACTOR


def reparam_lorentz(kind: str, wl_band, lam0: Optional[float] = None):
    """Fit ``eps(wl) ~= eps_inf + sigma*wl^2/(wl^2 - lam0^2)`` over `wl_band`.

    Returns ``(eps_inf, sigma, max_abs_dn)`` -- the last being the worst
    refractive-index error against the exact Sellmeier model across the band.

    **Why a re-fit rather than an algebraic split.** The obvious approach is
    to "move the bulk of the permittivity into eps_inf". Algebraically
    ``s*wl^2/(wl^2 - li^2) = s + s*li^2/(wl^2 - li^2)``, and that leading `s`
    is exactly what belongs in eps_inf -- but the remainder has an ``li^2``
    numerator, while meep's `LorentzianSusceptibility` is
    ``sigma*f0^2/(f0^2 - f^2)``, which in wavelength is
    ``sigma*wl^2/(wl^2 - lam0^2)`` -- a ``wl^2`` numerator. The exact split is
    therefore not expressible in meep and the model has to be re-fit into the
    form meep *can* represent. That form is what this function returns, so its
    output drops straight into a single `LorentzianSusceptibility`.

    **Per band, not global.** Silicon's telecom dispersion comes from a UV
    resonance at 0.3015 µm and no far-IR Lorentzian reproduces it across the
    whole range the examples span: fitting 1.26-1.70 µm at once gives max |dn|
    ~5.6e-3, worse than the ~2e-3 the dispersion is worth. Per band it is
    ~2e-4 (C-band) to ~4.9e-4 (O-band). This is why the media become a
    per-run quantity -- see `MaterialStack.si_medium()`.

    Reproduced end to end, with the tables, by
    `examples/verification/python/verify_item_8_reparameterisation.py`.
    """
    wl = np.atleast_1d(np.asarray(wl_band, dtype=float))
    if lam0 is None:
        lam0 = reparam_lam0(wl)
    if wl.size == 1:
        raise ValueError(
            "reparam_lorentz() needs a wavelength band, not a single point -- "
            "there is no dispersion to fit across one wavelength."
        )
    if float(wl.max()) * _REPARAM_MIN_LAM0_RATIO > lam0:
        raise ValueError(
            f"band reaches {float(wl.max()):.3f} um, too close to the fitted "
            f"pole lam0={lam0} um: the Lorentzian is near-resonant there and "
            f"the fit would not describe a smooth medium. Need "
            f"lam0 >= {_REPARAM_MIN_LAM0_RATIO} * max(band)."
        )
    grid = np.linspace(float(wl.min()), float(wl.max()), _REPARAM_SAMPLES)
    eps_true = sellmeier_index(kind, grid) ** 2
    basis = np.column_stack(
        [np.ones_like(grid), grid**2 / (grid**2 - lam0**2)])
    (eps_inf, sigma), *_ = np.linalg.lstsq(basis, eps_true, rcond=None)
    n_fit = np.sqrt(basis @ np.array([eps_inf, sigma]))
    max_dn = float(np.max(np.abs(n_fit - sellmeier_index(kind, grid))))
    return float(eps_inf), float(sigma), max_dn


def meep_medium_reparam(kind: str, wl_band, lam0: Optional[float] = None):
    """meep Medium for `kind`, dispersive across `wl_band` AND smoothable.

    The point of the re-parameterisation: `meep_medium()` below puts eps_inf=1
    and the entire index in susceptibilities, which meep leaves unsmoothed --
    the structure has `epsilon diagonal = (1,1,1)`, no contrast for subpixel
    averaging to act on, and that cost +0.0496 in n_eff when it was tried. Here
    eps_inf carries ~100% of the Si/SiO2 contrast, so smoothing is back, and
    the single Lorentzian carries only the wavelength dependence.

    Not cached: unlike `meep_medium()`/`meow_material()` this depends on the
    caller's band, so a global key would collide across runs.
    """
    import meep as mp

    # Resolve the pole BEFORE fitting, so the Lorentzian below is built at the
    # same lam0 the fit used. Deriving it twice, or letting `reparam_lorentz`
    # resolve it internally while this function kept `None`, would silently
    # pair one band's coefficients with another band's pole.
    if lam0 is None:
        lam0 = reparam_lam0(wl_band)
    eps_inf, sigma, _ = reparam_lorentz(kind, wl_band, lam0=lam0)
    return mp.Medium(
        epsilon=eps_inf,
        E_susceptibilities=[
            mp.LorentzianSusceptibility(
                frequency=1.0 / lam0, gamma=0, sigma=sigma)
        ],
    )


def _flat_meow_material(kind: str, wl: float):
    """meow material with a CONSTANT index, pinned at `wl`.

    The EME-side counterpart of the FDTD's frozen `mp.Medium(index=n)`. Built
    as a two-sample `SampledMaterial` with equal endpoints rather than an
    `IndexMaterial` so it stays the same meow type as the dispersive path —
    meow types `Structure3D.material` as a closed union, and keeping one type
    means the frozen and dispersive branches cannot diverge in how meow
    handles them.

    Spans `_MODEL_RANGE_UM` so any `Environment` wavelength the solvers use
    interpolates to exactly this index rather than falling off the end.
    """
    key = f"meow_flat:{kind}:{wl:.9f}"
    if key not in _material_cache:
        import meow as mw

        n = float(sellmeier_index(kind, wl))
        _material_cache[key] = mw.SampledMaterial(
            name=f"beef_{kind}_flat_{wl:.4f}",
            n=np.array([n, n], dtype=complex),
            params={"wl": np.asarray(_MODEL_RANGE_UM, dtype=float)},
        )
    return _material_cache[key]


def meep_medium(kind: str):
    """Dispersive meep Medium for `kind`, from the same coefficients.

    **Not currently used by `MaterialStack.si_medium()`/`sio2_medium()`** —
    turning it on cost 0.0496 in n_eff by disabling subpixel smoothing; see
    `MaterialStack.si_medium()`. Kept because it is correct as a model and is
    the starting point for the ε∞ re-parameterisation that would make
    dispersion usable.

    Built from raw susceptibilities rather than `meep.materials.Si` /
    `fused_quartz` so no `valid_freq_range` is attached — see
    `_MODEL_RANGE_UM` for why the O-band examples need that.
    """
    key = f"meep:{kind}"
    if key not in _material_cache:
        import meep as mp

        eps_inf, terms = _SELLMEIER[kind]
        _material_cache[key] = mp.Medium(
            epsilon=eps_inf,
            E_susceptibilities=[
                mp.LorentzianSusceptibility(frequency=1.0 / li, gamma=0, sigma=s)
                for li, s in terms
            ],
        )
    return _material_cache[key]


class MaterialStack:
    """A core waveguide material clad in a second material, both drawn from
    `_SELLMEIER`. Defaults to Si core / SiO2 clad, this codebase's original and
    still most common platform.

    Args:
        wg_height: core thickness (µm).
        sz: total vertical simulation span for meow cross-sections (µm).
        dispersive: build the FDTD media from the re-parameterised
            single-Lorentzian fit (`reparam_lorentz`) over the run's own band
            instead of freezing one index at band centre.

            **Default True as of 2026-08-01.** It was False while the
            re-parameterisation was verified only as a *model*; the strip-test
            gate has since passed on real meep solves across all four
            stack/band combinations (worst Δn_eff −0.00141 on SOI_160/O-band
            against a 0.002 tolerance), and the resolution question is answered
            — under dispersion the 1e-3 threshold wants res 80 at the band
            edge, which is what production already ships.

            Turning it on is what closes the engine-asymmetry gap: a broadband
            FDTD run froze one index at band centre while the EME re-solved per
            wavelength, so a hybrid chain contracted modes solved in two
            different materials — **0.00786 apart in n_Si at the O-band edges**,
            against a 2e-3 tolerance. Both-dispersive drops that to 0.001349.
            `HybridSolver`'s material-consistency guard was warn-only precisely
            because the old default failed it.

            Still switchable, and `False` on both engines remains the way to
            run a frozen A/B — see `EMESolver._extrusion_rules`.

            ⚠ **The strip-test gate above was measured on Si/SiO2 only.**
            `_REPARAM_LAM0_FACTOR = 1.15` (the pinned Lorentzian pole) has been
            validated on four (stack, band) combinations, all Si/SiO2. A
            non-default `core_kind` inherits `dispersive=True` without having
            been through that gate — see `core_kind` below.

        core_kind: key into `_SELLMEIER` for the core material. Default `"si"`
            preserves every existing stack's behaviour exactly.

            Added 2026-08-03 so a non-silicon core — nitride, so far — can reuse this
            class instead of a parallel one. Threaded through
            `extrusion_rules_at()`, `si_index()`/`si_medium()` (now aliases,
            see below) and `_medium()`'s dispersive dispatch.

            ⚠ **This does NOT mean a non-default core is validated.** It makes
            the class *able* to build one; A6/A7 in the validation doc list what
            still has to be measured — the subpixel-averaging strip-test gate
            above chief among them, since it has only ever been run at Si's
            much stronger core/clad contrast.
        clad_kind: key into `_SELLMEIER` for the cladding, both upper and
            substrate. Default `"sio2"`, same preservation rationale as
            `core_kind`.
    """

    def __init__(self, wg_height: float = 0.16, sz: float = 2.0,
                 dispersive: bool = True,
                 core_kind: str = "si", clad_kind: str = "sio2"):
        if core_kind not in _SELLMEIER:
            raise ValueError(
                f"core_kind={core_kind!r} is not in _SELLMEIER "
                f"({sorted(_SELLMEIER)}) -- add it there first"
            )
        if clad_kind not in _SELLMEIER:
            raise ValueError(
                f"clad_kind={clad_kind!r} is not in _SELLMEIER "
                f"({sorted(_SELLMEIER)}) -- add it there first"
            )
        self.wg_height = wg_height
        self.sz = sz
        self.dispersive = bool(dispersive)
        self.core_kind = core_kind
        self.clad_kind = clad_kind

    def resolve_dispersive(self, override: Optional[bool] = None) -> bool:
        """This stack's `dispersive` setting, or `override` when not None.

        The single place the three solvers agree on what "dispersive" means,
        so `FDTDSolver(dispersive=...)`, `EMESolver(dispersive=...)` and the
        stack cannot drift apart in how they interpret the flag.
        """
        return self.dispersive if override is None else bool(override)

    def extrusion_rules_at(self, wl: Optional[float] = None) -> dict:
        """meow GDS extrusion rules, optionally frozen at one wavelength.

        `wl=None` (the default, and what `extrusion_rules` gives) uses the
        dispersive `SampledMaterial`s — meow then re-evaluates the index at
        whatever wavelength each `Environment` names, which is why the EME
        side has never had the FDTD's frozen-band problem.

        Passing a wavelength substitutes **flat** materials pinned at that
        index. That is strictly a downgrade of the EME side and exists for one
        reason: `EMESolver(dispersive=False)` lets both engines be frozen the
        same way, so a hybrid A/B can separate "dispersion changed the answer"
        from "the two engines were modelling different materials all along"
        (`material_model_validation.md`, engine-agreement risk).
        """
        import meow as mw

        h = self.wg_height
        sz = self.sz
        if wl is None:
            core, clad = meow_material(self.core_kind), meow_material(self.clad_kind)
        else:
            core, clad = (_flat_meow_material(self.core_kind, float(wl)),
                          _flat_meow_material(self.clad_kind, float(wl)))
        return {
            (2, 0): [
                mw.GdsExtrusionRule(material=clad, h_min=-sz/2, h_max=-h/2, mesh_order=2),
                mw.GdsExtrusionRule(material=core, h_min=-h/2,  h_max=h/2,  mesh_order=1),
                mw.GdsExtrusionRule(material=clad, h_min=h/2,   h_max=sz/2, mesh_order=2),
            ],
            (1, 0): [
                mw.GdsExtrusionRule(material=clad, h_min=-sz/2, h_max=sz/2, mesh_order=2),
            ],
        }

    @property
    def extrusion_rules(self) -> dict:
        """meow GDS extrusion rules (layer → list[GdsExtrusionRule]).

        Dispersive. Kept as a property for the many existing callers; see
        `extrusion_rules_at()` for the frozen variant.
        """
        return self.extrusion_rules_at(None)

    def core_index(self, wl: float) -> float:
        """Refractive index of this stack's CORE material at wavelength wl
        (µm), from the shared Sellmeier model both engines use.

        Until 2026-07-28 the silicon case snapped to the nearest sample of
        meow's table (`argmin(|wls - wl|)`), which has only *two* entries
        between 1.45 and 1.65 µm — so the whole band below their midpoint got
        n = 3.4784 and everything above got 3.4710, a 0.0063 step mid-C-band
        plus a +0.00196 bias at 1.55 µm. That snap, duplicated in
        `meow/fde/meep.py::meep_material()`, is the Δn_eff kink recorded as
        directional_coupler_validation.md item 8. Both copies are gone.
        """
        return float(sellmeier_index(self.core_kind, wl))

    def clad_index(self, wl: float) -> float:
        """Refractive index of this stack's CLADDING material at wavelength wl
        (µm). See `core_index()`."""
        return float(sellmeier_index(self.clad_kind, wl))

    def si_index(self, wl: float) -> float:
        """Alias for `core_index()`, kept for every caller written before
        `core_kind` existed (2026-08-03). Despite the name this is the
        stack's CORE index, not necessarily silicon's — on a stack built with
        `core_kind="si3n4"` it returns nitride's index. See `core_kind` on the
        class docstring for why the alias, rather than a rename, was kept."""
        return self.core_index(wl)

    def sio2_index(self, wl: float) -> float:
        """Alias for `clad_index()`. See `si_index()`."""
        return self.clad_index(wl)

    def core_medium(self, wl, dispersive: Optional[bool] = None):
        """meep Medium for this stack's CORE material at wavelength wl (µm).

        **Non-dispersive, and that is now a measured choice rather than an
        oversight.** A broadband FDTD run gets one Medium for its whole band,
        so this frozen index is wrong away from `wl` — silicon really moves
        0.28% across 1.50–1.60 µm. Item 8 step 2 made both media dispersive to
        fix that, and it was **reverted after measurement**:

            450x220 nm Si strip in oxide, res 40, 1.55 µm
              non-dispersive  mp.Medium(index=n)        n_eff = 2.313398
              dispersive      eps_inf=1 + Lorentzians   n_eff = 2.363038
                                                        diff  = +0.049640

        meep subpixel-averages only the *instantaneous* ε and assigns
        susceptibility terms per pixel unsmoothed. meep's own library materials
        (and therefore the Sellmeier form lifted from them) put ε∞ = 1 and carry
        the entire index in susceptibilities — so the structure meep builds has
        `epsilon diagonal = (1,1,1)`, no dielectric contrast at the Si/SiO2
        interface, and subpixel smoothing does nothing. Losing it costs 0.0496
        in n_eff at res 40, against the ~0.002 the dispersion was worth: a 25x
        adverse trade on the same quantity, predicted before the attempt.

        The fix is not to abandon dispersion but to re-parameterise it so ε∞
        carries the bulk of the permittivity and the susceptibilities carry only
        its wavelength dependence — then smoothing acts on the contrast and the
        residual dispersive terms are small enough that leaving them unsmoothed
        does not matter. That needs its own fit and its own resolution sweep,
        which is why it is not bundled here.

        What *is* fixed: `core_index()` now comes from the shared analytic
        model, so the EME and FDTD agree on n(wl) to machine precision at any
        given wavelength — the frozen-band limitation remains, the silent
        engine-to-engine disagreement does not.

        **`dispersive=True` changes this** (the default since 2026-08-01; pass
        `dispersive=False` for the frozen-band behaviour described above). `wl` may
        then be the run's whole wavelength array, and the medium is built by
        `meep_medium_reparam()` — eps_inf carrying ~100% of the core/clad
        contrast so subpixel smoothing still acts, one Lorentzian carrying the
        band's wavelength dependence. A scalar `wl` still gives the plain
        frozen medium under either setting: there is no band to disperse over,
        and a single-wavelength mode solve wants exactly the index at that
        wavelength.

        `dispersive=None` (the default) uses this stack's own setting; pass
        True/False to override it for one call, which is how
        `FDTDSolver(dispersive=...)` overrides a shared stack without mutating
        it — the module-level `SOI_160`/`SOI_220` singletons are shared by
        every solver in the process.

        Reaches `meep_medium_reparam(self.core_kind, ...)` — **not** a literal
        `"si"` — since 2026-08-03. Before that, `si_medium()` passed the
        literal string `"si"` straight through, so any stack would have built
        a *silicon* Lorentzian under dispersion regardless of what its own
        `si_index()` returned. That silent mismatch was traced but never
        exercised (no non-Si stack existed yet); this is the fix, and
        `TestNonDefaultCoreClad` in `tests/test_material_reparam.py` is the
        regression guard.
        """
        return self._medium(self.core_kind, wl, dispersive)

    def clad_medium(self, wl, dispersive: Optional[bool] = None):
        """meep Medium for this stack's CLADDING material at wavelength wl
        (µm), or across a band when dispersion is enabled. See
        `core_medium()`."""
        return self._medium(self.clad_kind, wl, dispersive)

    def si_medium(self, wl, dispersive: Optional[bool] = None):
        """Alias for `core_medium()`. See `si_index()` for why this stays an
        alias rather than being renamed."""
        return self.core_medium(wl, dispersive)

    def sio2_medium(self, wl, dispersive: Optional[bool] = None):
        """Alias for `clad_medium()`. See `si_index()`."""
        return self.clad_medium(wl, dispersive)

    def _medium(self, kind: str, wl, dispersive: Optional[bool] = None):
        """Shared body of `core_medium()`/`clad_medium()`.

        `kind` is always `self.core_kind` or `self.clad_kind` by construction —
        every caller goes through one of those two, never a bare literal. Before
        2026-08-03 the frozen-index fallback dispatched on a literal
        `kind == "si"` (`self.si_index if kind == "si" else self.sio2_index`),
        which was harmless only because `kind` itself was always the literal
        `"si"`/`"sio2"` passed in by `si_medium()`/`sio2_medium()` — the same
        two strings `si_index`/`sio2_index` were hardcoded to. It would have
        broken the moment either side stopped being literally silicon or
        silica, which is exactly what `core_kind` now makes possible. Calling
        `sellmeier_index(kind, ...)` directly removes the redundant dispatch
        rather than generalizing it, since `kind` already names the right
        entry.

        Dispersion requires BOTH the flag and a real band — a scalar `wl` falls
        through to the frozen index under either setting, so single-wavelength
        callers (mode solves) are unaffected by the flag by construction.
        `dispersive=None` inherits the stack's own setting.
        """
        import meep as mp

        wl_arr = np.atleast_1d(np.asarray(wl, dtype=float))
        if (self.resolve_dispersive(dispersive)
                and wl_arr.size > 1 and wl_arr.min() < wl_arr.max()):
            return meep_medium_reparam(kind, wl_arr)
        return mp.Medium(index=float(sellmeier_index(kind, float(np.mean(wl_arr)))))


def resolved_material_config(material_stack, wavelengths: np.ndarray,
                             dispersive: Optional[bool] = None) -> Dict[str, Any]:
    """The *resolved* refractive indices this run will actually simulate with.

    `dispersive` is the CALLER's resolved setting — pass
    `stack.resolve_dispersive(solver.dispersive)`, not the stack's own flag.
    Omitting it falls back to the stack, which is correct only when the caller
    has no override.

    **Why that parameter exists (added 2026-08-01, measured):** every call site
    used to pass only the stack, so a per-solver override never reached the
    fingerprint. `FDTDSolver(material_stack=SOI_220, dispersive=True)` built a
    genuinely different medium — one Lorentzian susceptibility against a flat
    ε = 12.096 — and produced a **byte-identical cache key**, so a dispersive
    run collided with frozen-physics entries and was served the old physics.
    That is exactly the defect described below, recurring one level up: the
    stack-level hole was plugged, the solver-override hole was not.

    Same reasoning as `beef/solvers/fdtd.py::_resolved_spectral_config()`,
    applied to the material half — and here the gap it closes is real rather
    than precautionary. The cache records a `MaterialStack` as
    `{wg_height, sz}` (see `beef/cache/config.py::_solver_config()`), i.e. its
    *geometry* only, so nothing about the optical data reaches the fingerprint.
    Change the dispersion model, or how an index is derived from it, and every
    existing entry keeps matching and is served with the old physics.

    That is why `directional_coupler_validation.md` item 8 could not be fixed
    before this existed: the fix would have been silently ignored on every
    cached run and surfaced only on a cold cache.

    Sampled at the band edges *and* centre, not just the centre: both media are
    dispersive now, so one number no longer describes what the run uses. Three
    points pin the slope as well as the offset, so swapping in a different
    model with a coincidentally equal midpoint still moves the key.

    Deliberately **not** a `MaterialStack.cache_key()`: resolving an index used
    to import meow, meow imports meep (measured — `import meow` alone puts meep
    in `sys.modules`), and `beef/cache/store.py::_is_master()` dispatches on
    whether meep is loaded while promising to import nothing itself. A
    `cache_key()` would therefore flip rank detection off its launcher-env
    fallback as a side effect of building a config, electing every rank of an
    EME-only job as master — the item 4 bug. Called from a solver instead, both
    are already loaded and the promise holds.
    """
    wl = np.asarray(wavelengths, dtype=float)
    probes = (float(wl.min()), float(np.mean(wl)), float(wl.max()))
    # `si_index()`/`sio2_index()` are the 2026-08-03 core/clad aliases (see
    # `MaterialStack.si_index()`) -- they already resolve to the real core/clad
    # index for a non-Si stack, e.g. nitride's ~2.00 rather than silicon's
    # ~3.48. That is what makes the FROZEN half of this key correctly
    # stack-agnostic already.
    config = {
        "resolved_probe_wl": [round(w, 9) for w in probes],
        "resolved_si_index": [round(material_stack.si_index(w), 9) for w in probes],
        "resolved_sio2_index": [round(material_stack.sio2_index(w), 9) for w in probes],
    }

    # The probe indices above come from `sellmeier_index()`, which the
    # dispersive re-parameterisation does NOT change -- so without the
    # following, flipping `dispersive` on would leave this key byte-identical
    # and every existing cache entry would keep matching and be served the old
    # non-dispersive physics. That is precisely the defect this function was
    # created to kill, recurring one level up
    # (`material_model_validation.md` work item 3).
    #
    # Recorded directly as the fitted (eps_inf, sigma) rather than as a
    # model-version tag, so any change to the fit -- band, pinned lam0, sample
    # count, or the Sellmeier coefficients underneath it -- moves the key on
    # its own, with nothing to remember to bump.
    #
    # FIXED 2026-08-03 (`mmi_1x4_sin300_validation.md` step 3). `label` is a
    # fixed ("si", "sio2") pair purely for the KEY NAMES below
    # (`resolved_si_reparam`/`resolved_sio2_reparam`) -- kept literal
    # deliberately, not renamed to `resolved_core_reparam`/`resolved_clad_
    # reparam`, so a default (core_kind="si", clad_kind="sio2") stack's cache
    # key is byte-for-byte what it always was; only the newly-possible
    # non-default case gets a new key. `kind` is what actually gets fit --
    # `material_stack.core_kind`/`clad_kind` -- which before this fix was the
    # SAME hardcoded literal as `label`, so `reparam_lorentz("si", wl)` ran
    # for every stack regardless of what its own `core_medium()` built. On a
    # nitride stack that meant the cache key recorded SILICON's (eps_inf,
    # sigma) while the medium actually simulated was nitride's -- the
    # `resolved_si_index` half of this config, one field up, had already been
    # stack-correct since step 2; this was the still-hardcoded other half.
    # `TestNonDefaultCoreClad.test_resolved_material_config_still_has_the_
    # step3_gap` in `tests/test_material_reparam.py` demonstrated this on real
    # code before the fix and is updated alongside it to assert the fix
    # instead (see that test for the historical record).
    resolved = (getattr(material_stack, "dispersive", False) if dispersive is None
                else bool(dispersive))
    if resolved and wl.min() < wl.max():
        config["resolved_dispersive"] = True
        for label, kind in (("si", material_stack.core_kind),
                            ("sio2", material_stack.clad_kind)):
            eps_inf, sigma, _ = reparam_lorentz(kind, wl)
            config[f"resolved_{label}_reparam"] = [
                round(eps_inf, 9), round(sigma, 9), round(reparam_lam0(wl), 9)]
    return config


# Pre-built stacks for common platforms.
SOI_160 = MaterialStack(wg_height=0.16, sz=2.0)   # 160 nm Si (default in this codebase)
SOI_220 = MaterialStack(wg_height=0.22, sz=2.0)   # 220 nm Si (standard commercial SOI)

# 300 nm Si3N4 on SiO2. Added 2026-08-03 for the 1x4 MMI build.
#
# sz=5.0 MEASURED 2026-08-03 (step 4, job 11660166): converged (|dn_eff| < 1e-4
# going sz 4.0->5.0) at the actual design width, 1.1 um -- NOT the placeholder
# 2.0 inherited from the SOI stacks, which A7 correctly flagged as suspect
# (nitride's core/clad eps contrast is ~5.3x weaker than Si/SiO2's, so its
# mode reaches further into the cladding). Fundamental n_eff plateaus at
# 1.587348 (margin +0.14332 above n_clad=1.444024); a second guided mode at
# 1.510616 (margin +0.06659) -- 1.1 um is a genuinely 2-mode width, not
# single-mode (the true single-mode window, 0.2-0.35 um, turned out to be so
# weakly confined it would not converge even at sz=10 um -- see the doc's
# step 4 writeup for the full measurement and why 1.1 um was chosen instead).
#
# `dispersive` is still the class default, True -- see `core_kind` on
# MaterialStack's docstring for why that inherits the flag without having
# been through the Si/SiO2-only strip-test gate.
SIN_300 = MaterialStack(wg_height=0.30, sz=5.0, core_kind="si3n4")
