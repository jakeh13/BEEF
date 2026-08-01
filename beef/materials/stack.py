"""MaterialStack: single source of truth for Si/SiO2 layer definitions.

Replaces the duplicated stack_extrusions_GF() and meep_material_from_meow()
calls scattered across every master and sweep script.
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
# These are lifted verbatim from meep's own materials library, converting its
# Lorentzian form (frequency = 1/lambda_i, gamma = 0) into wavelength terms —
# verified to reproduce `meep.materials.Si` / `fused_quartz` to 4e-16 over
# 1.4-1.7 um. Holding the numbers here rather than importing them means meow
# and meep are driven from one source instead of two, which is the whole point:
# before 2026-07-28 the EME interpolated meow's sampled tables while the FDTD
# froze a single snapped index, and nothing kept the two consistent.
#
#   Si   - Salzberg & Villa, JOSA 47, 244 (1957)
#   SiO2 - Malitson, JOSA 55, 1205 (1965), fused silica
#
# meow's own `silicon_oxide` table *is* Malitson: they agree to 0.00000 across
# 1.50-1.60 um. Its `silicon` table is a different, noisier source sitting
# ~0.002 above Salzberg, so adopting this model shifts n_Si by that much —
# a deliberate change, not a refactor. It buys a smooth analytic n(lambda)
# where the table had only two samples between 1.45 and 1.65 um.
_SELLMEIER: Dict[str, Any] = {
    "si": (1.0, ((0.301516485, 10.6684293),
                 (1.13475115, 0.003043475),
                 (1104.0, 1.54133408))),
    "sio2": (1.0, ((0.0684043, 0.6961663),
                   (0.1162414, 0.4079426),
                   (9.896161, 0.8974794))),
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

    **Why a re-fit rather than an algebraic split.** Item 8 / `notes.md` #5
    prescribe "move the bulk of the permittivity into eps_inf". Algebraically
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
    `examples/verification/python/verify_item_8_reparameterisation.py`
    (`validation_markdown/material_model_validation.md` work item 1).
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
    """Si waveguide on SiO2 substrate with SiO2 upper cladding.

    Args:
        wg_height: Si core thickness (µm).
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
    """

    def __init__(self, wg_height: float = 0.16, sz: float = 2.0,
                 dispersive: bool = True):
        self.wg_height = wg_height
        self.sz = sz
        self.dispersive = bool(dispersive)

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
            si, ox = meow_material("si"), meow_material("sio2")
        else:
            si, ox = (_flat_meow_material("si", float(wl)),
                      _flat_meow_material("sio2", float(wl)))
        return {
            (2, 0): [
                mw.GdsExtrusionRule(material=ox, h_min=-sz/2, h_max=-h/2, mesh_order=2),
                mw.GdsExtrusionRule(material=si, h_min=-h/2,  h_max=h/2,  mesh_order=1),
                mw.GdsExtrusionRule(material=ox, h_min=h/2,   h_max=sz/2, mesh_order=2),
            ],
            (1, 0): [
                mw.GdsExtrusionRule(material=ox, h_min=-sz/2, h_max=sz/2, mesh_order=2),
            ],
        }

    @property
    def extrusion_rules(self) -> dict:
        """meow GDS extrusion rules (layer → list[GdsExtrusionRule]).

        Dispersive. Kept as a property for the many existing callers; see
        `extrusion_rules_at()` for the frozen variant.
        """
        return self.extrusion_rules_at(None)

    def si_index(self, wl: float) -> float:
        """Refractive index of silicon at wavelength wl (µm), from the shared
        Sellmeier model both engines now use.

        Until 2026-07-28 this snapped to the nearest sample of meow's table
        (`argmin(|wls - wl|)`), which has only *two* entries between 1.45 and
        1.65 µm — so the whole band below their midpoint got n = 3.4784 and
        everything above got 3.4710, a 0.0063 step mid-C-band plus a +0.00196
        bias at 1.55 µm. That snap, duplicated in
        `meow/fde/meep.py::meep_material()`, is the Δn_eff kink recorded as
        directional_coupler_validation.md item 8. Both copies are gone.
        """
        return float(sellmeier_index("si", wl))

    def sio2_index(self, wl: float) -> float:
        """Refractive index of SiO2 at wavelength wl (µm). See `si_index()`."""
        return float(sellmeier_index("sio2", wl))

    def si_medium(self, wl, dispersive: Optional[bool] = None):
        """meep Medium for silicon at wavelength wl (µm).

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
        adverse trade on the same quantity. `notes.md` recommendation #9
        predicted exactly this before the attempt.

        The fix is not to abandon dispersion but to re-parameterise it so ε∞
        carries the bulk of the permittivity and the susceptibilities carry only
        its wavelength dependence — then smoothing acts on the contrast and the
        residual dispersive terms are small enough that leaving them unsmoothed
        does not matter. That needs its own fit and its own resolution sweep,
        which is why it is not bundled here.

        What *is* fixed: `si_index()` now comes from the shared analytic model,
        so the EME and FDTD agree on n(wl) to machine precision at any given
        wavelength — the frozen-band limitation remains, the silent
        engine-to-engine disagreement does not.

        **`dispersive=True` changes this** (opt-in, off by default). `wl` may
        then be the run's whole wavelength array, and the medium is built by
        `meep_medium_reparam()` — eps_inf carrying ~100% of the Si/SiO2
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
        """
        return self._medium("si", wl, dispersive)

    def sio2_medium(self, wl, dispersive: Optional[bool] = None):
        """meep Medium for SiO2 at wavelength wl (µm), or across a band when
        dispersion is enabled. See `si_medium()`."""
        return self._medium("sio2", wl, dispersive)

    def _medium(self, kind: str, wl, dispersive: Optional[bool] = None):
        """Shared body of `si_medium()`/`sio2_medium()`.

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
        index_of = self.si_index if kind == "si" else self.sio2_index
        return mp.Medium(index=index_of(float(np.mean(wl_arr))))


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
    resolved = (getattr(material_stack, "dispersive", False) if dispersive is None
                else bool(dispersive))
    if resolved and wl.min() < wl.max():
        config["resolved_dispersive"] = True
        for kind in ("si", "sio2"):
            eps_inf, sigma, _ = reparam_lorentz(kind, wl)
            config[f"resolved_{kind}_reparam"] = [
                round(eps_inf, 9), round(sigma, 9), round(reparam_lam0(wl), 9)]
    return config


# Pre-built stacks for common SOI platforms.
SOI_160 = MaterialStack(wg_height=0.16, sz=2.0)   # 160 nm Si (default in this codebase)
SOI_220 = MaterialStack(wg_height=0.22, sz=2.0)   # 220 nm Si (standard commercial SOI)
