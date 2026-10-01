"""Tests for the eps_inf re-parameterisation of the Sellmeier model.

These are meep-free except where a `Medium` is actually built.

The property under test throughout: turning dispersion on must be an *opt-in*
change that (a) leaves every existing caller byte-identical while off, (b)
restores subpixel-smoothable dielectric contrast when on -- the specific thing
whose absence caused the first attempt to be reverted -- and (c) cannot be
turned on without moving the cache fingerprint.
"""
import numpy as np
import pytest

from beef.materials.stack import (
    _MODEL_RANGE_UM,
    _SELLMEIER,
    MaterialStack,
    SOI_220,
    reparam_lam0,
    meep_medium_reparam,
    reparam_lorentz,
    resolved_material_config,
    sellmeier_index,
)

C_BAND = np.linspace(1.50, 1.60, 21)
O_BAND = np.linspace(1.26, 1.36, 21)


def _si_material(rules):
    """The silicon material out of a meow extrusion-rule dict.

    Layer (2, 0) is [oxide below, silicon core, oxide above]; the core is the
    mesh_order=1 entry.
    """
    return next(r.material for r in rules[(2, 0)] if r.mesh_order == 1)


class TestSellmeierMatchesMeepLibrary:
    """`_SELLMEIER` must reproduce meep's own library materials.

    The module comment above `_SELLMEIER` has always claimed the entries were
    "lifted verbatim from meep's own materials library ... verified to reproduce
    `meep.materials.Si` / `fused_quartz` to 4e-16" — but until 2026-08-03
    **nothing tested it**. The whole point of holding the coefficients locally
    is that meow and meep are driven from one source; if a transcription is off,
    both engines are consistently wrong together and every validation check
    is measuring the wrong material in agreement.

    Added with the `si3n4` entry (nitride step 1) and deliberately written to
    cover all three, not just the new one: the guard the new material needed was
    the guard the old ones never had.

    Ranges are per material because meep's own `valid_freq_range` differs.
    `meep.materials.Si` refuses below 1.36 um — which is exactly why
    `meep_medium()` builds susceptibilities by hand instead of importing it, and
    why the examples can run the O-band down to 1.26 um at all.

    ⚠ **The tolerances are NOT all machine precision, and that is a finding, not
    a slack bound.** Writing this test is what turned up a transcription
    truncation in the SILICON entry that had been there since the model was
    created: `_SELLMEIER["si"]` term 2 carries `sigma = 0.003043475` where meep
    has `0.0030434748` — 7 significant figures against meep's 8. MEASURED
    2026-08-03: that makes n_Si differ from `meep.materials.Si` by up to
    **8.35e-11** over 1.4-1.7 um, so the module comment's claim of "4e-16" is
    **false for silicon** (it holds for SiO2 and Si3N4, both of which pass at
    <1e-12 here).

    **Deliberately not fixed, and the reason is cache economics, not doubt.**
    8.35e-11 in n is ~7 orders of magnitude below the ~1e-3 effects every
    validation check is built to resolve — it is physically nothing.
    But correcting it MOVES THE CACHE FINGERPRINT: measured, three entries shift
    in their 9th decimal (`resolved_si_index` at 1.70 um, and
    `resolved_si_reparam`'s eps_inf on BOTH the C- and O-band), which would
    invalidate every dispersive cache entry and all four goldens repo-wide. That
    is a deliberate, scheduled event in this project, not a side effect of
    adding a material. Fix it when a regeneration is happening anyway, and tighten
    `TOL["si"]` to 1e-12 in the same commit.
    """

    # (kind, meep attribute, wl range to compare over)
    CASES = [
        ("si", "Si", (1.40, 1.70)),                 # meep valid 1.36-11.0
        ("sio2", "fused_quartz", _MODEL_RANGE_UM),  # meep valid 0.21-6.7
        ("si3n4", "Si3N4_NIR", _MODEL_RANGE_UM),    # meep valid 0.31-5.504
    ]

    # Bounds sit just above each material's MEASURED drift, so a real
    # coefficient change trips them. See the class docstring for why si's is not
    # 1e-12 — it is pinning a known truncation, not tolerating unknown slop.
    TOL = {"si": 1e-9, "sio2": 1e-12, "si3n4": 1e-12}

    @pytest.mark.parametrize("kind,attr,wl_range", CASES)
    def test_reproduces_meep_library(self, kind, attr, wl_range):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        import meep.materials as mm

        med = getattr(mm, attr)
        wls = np.linspace(*wl_range, 201)
        ours = sellmeier_index(kind, wls)
        theirs = np.array(
            [np.sqrt(float(np.real(med.epsilon(1.0 / w)[0][0]))) for w in wls])
        worst = float(np.max(np.abs(ours - theirs)))
        assert worst < self.TOL[kind], (
            f"{kind} drifted from meep.materials.{attr} by {worst:.3e} over "
            f"{wl_range[0]}-{wl_range[1]} um -- the local coefficients no "
            "longer match the library they were lifted from"
        )

    def test_the_new_entries_really_are_at_machine_precision(self):
        """Separated from the parametrised case so the si truncation above can
        never be read as licence for a sloppy transcription. Anything added to
        `_SELLMEIER` from meep's library should land here, not in the 1e-9 tier.
        """
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        import meep.materials as mm

        for kind, attr in (("sio2", "fused_quartz"), ("si3n4", "Si3N4_NIR")):
            med = getattr(mm, attr)
            wls = np.linspace(*_MODEL_RANGE_UM, 201)
            worst = float(np.max(np.abs(
                sellmeier_index(kind, wls)
                - np.array([np.sqrt(float(np.real(med.epsilon(1.0 / w)[0][0])))
                            for w in wls]))))
            assert worst < 1e-14, f"{kind} at {worst:.3e}, expected ~1e-16"

    def test_every_entry_is_lossless_and_eps_inf_is_one(self):
        """The stored form has no gamma, so a lossy source material cannot be
        represented — it would be silently transcribed as lossless.

        This is not hypothetical: meep's plain `Si3N4` has gamma = 1.4413, and
        it is the model someone reaching by name would pick. Dropping its gamma
        would look like a working nitride and be a different material.
        """
        for kind, (eps_inf, terms) in _SELLMEIER.items():
            assert eps_inf == 1.0, (
                f"{kind} has eps_inf={eps_inf}; meep's library materials put "
                "eps_inf=1 and carry the index in susceptibilities, and "
                "`reparam_lorentz` assumes that convention"
            )
            assert len(terms) >= 1
            for lam_i, sigma in terms:
                assert lam_i > 0 and sigma != 0

    def test_nitride_index_is_where_it_should_be(self):
        """Anchor values, so a plausible-looking wrong transcription (e.g. the
        VISNIR single-term fit) does not pass the range test by luck."""
        assert float(sellmeier_index("si3n4", 1.31)) == pytest.approx(2.003130, abs=5e-7)
        assert float(sellmeier_index("si3n4", 1.55)) == pytest.approx(1.996280, abs=5e-7)

    def test_nitride_sits_between_oxide_and_silicon(self):
        """Sanity on the platform this was added for: a SiN core must still
        guide in an SiO2 cladding, and must not be mistaken for silicon."""
        n_ox = float(sellmeier_index("sio2", 1.55))
        n_sin = float(sellmeier_index("si3n4", 1.55))
        n_si = float(sellmeier_index("si", 1.55))
        assert n_ox < n_sin < n_si
        assert n_sin - n_ox == pytest.approx(0.5523, abs=1e-3)

    @pytest.mark.parametrize("attr", ["Si3N4", "Si3N4_VISNIR"])
    def test_the_rejected_nitride_models_really_are_unusable_at_telecom(self, attr):
        """Pins the reason `Si3N4_NIR` was chosen, so the block comment in
        `stack.py` is checked rather than merely asserted. If a future meep
        widens these ranges this fails loudly and the choice can be revisited
        on purpose instead of by accident."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        import meep.materials as mm

        med = getattr(mm, attr)
        for wl in (1.31, 1.55):
            with pytest.raises(ValueError):
                med.epsilon(1.0 / wl)


class TestNitrideReparam:
    """The new material has to survive the same dispersive path as Si/SiO2.

    Passing these does NOT mean nitride is ready to simulate — the
    subpixel-averaging strip gate is a real-meep measurement and has never been
    run for this stack (`mmi_1x4_sin300_validation.md` A6). These only pin that
    the fit machinery accepts nitride and produces a sane, smoothable medium.
    """

    def test_fit_is_at_least_as_good_as_silicons(self):
        """Nitride's telecom dispersion is far weaker than silicon's, so its
        single-Lorentzian residual should be smaller, not larger. A regression
        here means the pinned-pole rule stopped suiting this material."""
        for band in (C_BAND, O_BAND):
            _, _, dn_sin = reparam_lorentz("si3n4", band)
            _, _, dn_si = reparam_lorentz("si", band)
            assert dn_sin < dn_si, (
                f"nitride fit {dn_sin:.2e} is worse than silicon's {dn_si:.2e}")
            assert dn_sin < 1e-3

    def test_derived_pole_clears_the_guard(self):
        for band in (C_BAND, O_BAND):
            reparam_lorentz("si3n4", band)      # must not raise

    def test_eps_inf_carries_the_contrast_against_oxide(self):
        """The reverted-attempt guard (see TestSmoothableContrastIsRestored),
        repeated for the nitride/oxide pair — the contrast is 5.3x weaker here,
        which is precisely why it cannot be inherited from silicon's pass."""
        e_sin, _, _ = reparam_lorentz("si3n4", C_BAND)
        e_ox, _, _ = reparam_lorentz("sio2", C_BAND)
        true_contrast = float(sellmeier_index("si3n4", 1.55) ** 2
                              - sellmeier_index("sio2", 1.55) ** 2)
        frac = (e_sin - e_ox) / true_contrast
        assert 0.95 < frac < 1.10, (
            f"only {frac:.1%} of the SiN/SiO2 eps contrast is in eps_inf")

    def test_built_medium_tracks_sellmeier_across_the_band(self):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        med = meep_medium_reparam("si3n4", C_BAND)
        for wl in (1.50, 1.55, 1.60):
            n_meep = float(np.sqrt(med.epsilon(1.0 / wl)[0][0].real))
            assert abs(n_meep - float(sellmeier_index("si3n4", wl))) < 1e-3
        assert len(med.E_susceptibilities) == 1
        assert med.epsilon_diag.x > 3.5, "eps_inf collapsed toward 1"


class TestNonDefaultCoreClad:
    """`MaterialStack(core_kind=..., clad_kind=...)`.

    The class was hardwired Si-core/SiO2-clad in four places before
    2026-08-03: `_SELLMEIER` lookups in `extrusion_rules_at()`, the dispatch in
    `_medium()`, and the literal `"si"`/`"sio2"` strings `si_medium()`/
    `sio2_medium()` passed through. `si_index`/`sio2_index`/`si_medium`/
    `sio2_medium` are kept as aliases for the new `core_*`/`clad_*` names
    (deliberately, not renamed — see `MaterialStack.si_index()`), so every
    existing caller (`hybrid.py`'s `getattr(stack, "si_index", ...)` among
    them) keeps working unchanged while now resolving to the real core/clad
    material instead of a hardcoded one.
    """

    def test_invalid_core_kind_raises_immediately(self):
        with pytest.raises(ValueError, match="core_kind"):
            MaterialStack(wg_height=0.30, core_kind="gaas")

    def test_invalid_clad_kind_raises_immediately(self):
        with pytest.raises(ValueError, match="clad_kind"):
            MaterialStack(wg_height=0.30, clad_kind="gaas")

    def test_defaults_are_silicon_on_oxide(self):
        """The whole point of keeping core_kind/clad_kind as keyword-only
        parameters with these specific defaults: every stack built the old way
        (`MaterialStack(wg_height=...)`, no mention of core/clad) must behave
        exactly as it did before this existed."""
        stack = MaterialStack(wg_height=0.22)
        assert stack.core_kind == "si"
        assert stack.clad_kind == "sio2"

    def test_index_aliases_resolve_to_the_real_core_not_hardcoded_silicon(self):
        """The A4 trap, made concrete. A naive alias-only implementation could
        satisfy every OLD test while still returning silicon's index for a
        nitride stack if `si_index()` stayed hardcoded to `sellmeier_index("si",
        ...)` rather than `sellmeier_index(self.core_kind, ...)`."""
        sin = MaterialStack(wg_height=0.30, core_kind="si3n4")
        assert sin.si_index(1.55) == pytest.approx(sellmeier_index("si3n4", 1.55))
        assert sin.si_index(1.55) != pytest.approx(sellmeier_index("si", 1.55))
        assert sin.core_index(1.55) == sin.si_index(1.55)

    def test_clad_aliases_still_work_symmetrically(self):
        """Exercises clad_kind on the material that is USUALLY the default —
        a stack with an unusual core but the standard oxide clad, the expected
        shape for the 1x4 device (nitride core, SiO2 clad)."""
        sin = MaterialStack(wg_height=0.30, core_kind="si3n4")
        assert sin.sio2_index(1.55) == pytest.approx(sellmeier_index("sio2", 1.55))
        assert sin.clad_index(1.55) == sin.sio2_index(1.55)

    def test_dispersive_medium_carries_the_real_cores_fit_not_silicons(self):
        """THE load-bearing test in this class, and the one the validation doc
        calls "the single most important assertion in the whole plan" (A4 /
        step 2 evidence #2).

        Before 2026-08-03, `si_medium()` called `self._medium("si", wl, ...)`
        with a LITERAL `"si"`, so under `dispersive=True` ANY stack would have
        built `meep_medium_reparam("si", wl)` — silicon's Lorentzian —
        regardless of what `core_kind` said. That mismatch was traced in the
        validation doc but never exercised, because no non-Si stack existed to
        exercise it with. This is that exercise.
        """
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        sin = MaterialStack(wg_height=0.30, core_kind="si3n4", dispersive=True)
        med = sin.si_medium(C_BAND)               # the pre-existing alias name
        expect_eps_inf, expect_sigma, _ = reparam_lorentz("si3n4", C_BAND)
        wrong_eps_inf, wrong_sigma, _ = reparam_lorentz("si", C_BAND)

        assert med.epsilon_diag.x == pytest.approx(expect_eps_inf, abs=1e-6)
        assert med.E_susceptibilities[0].sigma_diag.x == pytest.approx(
            expect_sigma, abs=1e-6)
        assert med.epsilon_diag.x != pytest.approx(wrong_eps_inf, abs=1e-3), (
            "medium carries SILICON's eps_inf on a nitride-core stack -- the "
            "A4 regression"
        )
        assert med.E_susceptibilities[0].sigma_diag.x != pytest.approx(
            wrong_sigma, abs=1e-3
        ), "medium carries SILICON's sigma on a nitride-core stack -- the A4 regression"

    def test_core_medium_and_si_medium_alias_agree(self):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        sin = MaterialStack(wg_height=0.30, core_kind="si3n4", dispersive=True)
        a, b = sin.core_medium(C_BAND), sin.si_medium(C_BAND)
        assert a.epsilon_diag.x == b.epsilon_diag.x
        assert a.E_susceptibilities[0].sigma_diag.x == b.E_susceptibilities[0].sigma_diag.x

    def test_frozen_scalar_medium_is_also_the_real_core(self):
        """The non-dispersive / single-wavelength path goes through a
        different branch of `_medium()` -- pin it separately from the
        dispersive case above."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        sin = MaterialStack(wg_height=0.30, core_kind="si3n4", dispersive=False)
        med = sin.si_medium(1.55)
        assert not med.E_susceptibilities
        assert float(med.epsilon_diag.x) == pytest.approx(
            sellmeier_index("si3n4", 1.55) ** 2, abs=1e-9)

    def test_extrusion_rules_use_the_real_core_and_clad_materials(self):
        """`extrusion_rules_at()` -- the EME-side path, independent of
        `_medium()`. Checks the actual `SampledMaterial.n` meow will mesh with,
        not just that the call succeeds."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        sin = MaterialStack(wg_height=0.30, sz=3.0, core_kind="si3n4")
        rules = sin.extrusion_rules_at(1.55)
        core = next(r.material for r in rules[(2, 0)] if r.mesh_order == 1)
        clads = [r.material for r in rules[(2, 0)] if r.mesh_order == 2]
        assert float(np.real(core.n[0])) == pytest.approx(
            sellmeier_index("si3n4", 1.55), abs=1e-9)
        for clad in clads:
            assert float(np.real(clad.n[0])) == pytest.approx(
                sellmeier_index("sio2", 1.55), abs=1e-9)
        # (1, 0) is the etch/background layer -- clad on both stacks, unrelated
        # to core_kind. Confirms the generalization didn't leak into it.
        bg = rules[(1, 0)][0].material
        assert float(np.real(bg.n[0])) == pytest.approx(
            sellmeier_index("sio2", 1.55), abs=1e-9)

    def test_sz_and_wg_height_still_reach_the_extrusion_geometry(self):
        """Regression: core/clad generalization must not disturb the
        geometric (h_min/h_max) side of extrusion_rules_at()."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        sin = MaterialStack(wg_height=0.30, sz=4.0, core_kind="si3n4")
        rules = sin.extrusion_rules_at(1.55)
        core_rule = next(r for r in rules[(2, 0)] if r.mesh_order == 1)
        assert core_rule.h_min == pytest.approx(-0.15)
        assert core_rule.h_max == pytest.approx(0.15)
        # Two clad slabs -- substrate [-sz/2, -h/2] and upper [h/2, sz/2] -- so
        # h_min/h_max differ between them; compare as the pair of intervals.
        clad_rules = [r for r in rules[(2, 0)] if r.mesh_order == 2]
        intervals = {(round(r.h_min, 6), round(r.h_max, 6)) for r in clad_rules}
        assert intervals == {(-2.0, -0.15), (0.15, 2.0)}


class TestExistingStacksUnchangedByGeneralization:
    """Evidence #1 for step 2: SOI_160/SOI_220 (and any default-constructed
    MaterialStack) must be byte-for-byte identical to what they were before
    `core_kind`/`clad_kind` existed. This is a regression suite, not new
    coverage -- every value here is what the pre-refactor hardcoded methods
    would have returned, computed independently via `sellmeier_index()` rather
    than by importing the old code.
    """

    def test_si_index_matches_the_shared_model_directly(self):
        for wl in (1.26, 1.31, 1.50, 1.55, 1.60):
            assert SOI_220.si_index(wl) == pytest.approx(sellmeier_index("si", wl))
            assert SOI_220.sio2_index(wl) == pytest.approx(sellmeier_index("sio2", wl))

    def test_core_and_clad_aliases_agree_with_si_sio2_aliases(self):
        for wl in (1.31, 1.55):
            assert SOI_220.core_index(wl) == SOI_220.si_index(wl)
            assert SOI_220.clad_index(wl) == SOI_220.sio2_index(wl)

    def test_dispersive_medium_unchanged(self):
        """The exact assertion `TestSolverFlags.test_fdtd_default_stack_gives_a_
        dispersive_solver` already made, repeated here at the MaterialStack
        level so this class stands alone as the step-2 regression record."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        med = SOI_220.si_medium(C_BAND, dispersive=True)
        eps_inf, sigma, _ = reparam_lorentz("si", C_BAND)
        assert med.epsilon_diag.x == pytest.approx(eps_inf, abs=1e-9)
        assert med.E_susceptibilities[0].sigma_diag.x == pytest.approx(sigma, abs=1e-9)

    def test_extrusion_rules_unchanged(self):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        rules = SOI_220.extrusion_rules_at(1.55)
        core = next(r.material for r in rules[(2, 0)] if r.mesh_order == 1)
        assert float(np.real(core.n[0])) == pytest.approx(
            sellmeier_index("si", 1.55), abs=1e-9)
        assert core.name == "beef_si_flat_1.5500"

    def test_sin_300_singleton_exists_and_is_nitride(self):
        from beef.materials.stack import SIN_300

        assert SIN_300.core_kind == "si3n4"
        assert SIN_300.clad_kind == "sio2"
        assert SIN_300.wg_height == pytest.approx(0.30)
        assert SIN_300.dispersive is True

    def test_resolved_material_config_now_records_the_real_cores_reparam_fit(self):
        """FIXED 2026-08-03 (step 3). Until this fix landed, this test
        (formerly `..._still_has_the_step3_gap`) pinned the opposite: on real
        code, `resolved_material_config()`'s dispersive branch recorded
        SILICON's (eps_inf, sigma) for a nitride stack's cache key, while
        `core_medium()`/`si_medium()` (fixed in step 2) had already been
        correctly building NITRIDE's medium -- the key had stopped tracking
        the thing it was supposed to fingerprint. That failure mode is why
        this class exists at all: recorded here, not just in a comment, so a
        future regression on either half is caught the same way this one was.

        `resolved_si_reparam` is a fixed KEY NAME (same reasoning as
        `resolved_si_index`/`resolved_sio2_index` since step 2) -- it is now
        the CORE slot's reparam fit regardless of what the core material
        actually is, not literally silicon's.
        """
        sin = MaterialStack(wg_height=0.30, core_kind="si3n4", dispersive=True)
        cfg = resolved_material_config(sin, C_BAND)
        nitride_fit = reparam_lorentz("si3n4", C_BAND)
        silicon_fit = reparam_lorentz("si", C_BAND)
        assert cfg["resolved_si_reparam"][:2] == [
            round(nitride_fit[0], 9), round(nitride_fit[1], 9)]
        assert cfg["resolved_si_reparam"][:2] != [
            round(silicon_fit[0], 9), round(silicon_fit[1], 9)], (
            "cache key still carries SILICON's reparam fit on a nitride-core "
            "stack -- the step-3 regression"
        )

    def test_resolved_material_config_clad_reparam_also_uses_the_real_clad(self):
        """Symmetric coverage on the clad side, at a stack whose clad is NOT
        the default oxide -- catches a fix that only threads core_kind through
        and leaves clad_kind hardcoded."""
        stack = MaterialStack(wg_height=0.30, core_kind="sio2", clad_kind="si3n4",
                              dispersive=True)
        cfg = resolved_material_config(stack, C_BAND)
        nitride_fit = reparam_lorentz("si3n4", C_BAND)
        assert cfg["resolved_sio2_reparam"][:2] == [
            round(nitride_fit[0], 9), round(nitride_fit[1], 9)]

    def test_resolved_material_config_is_unchanged_for_default_stacks(self):
        """Evidence the fix is INERT for every existing device: for
        core_kind="si"/clad_kind="sio2" (the class default), the key names
        (`resolved_si_reparam`/`resolved_sio2_reparam`) AND the values fitted
        into them must be byte-identical to what a literal `("si", "sio2")`
        loop always produced -- no existing cache entry or golden should move
        because of this fix."""
        default = MaterialStack(wg_height=0.22, dispersive=True)
        cfg = resolved_material_config(default, C_BAND)
        si_fit = reparam_lorentz("si", C_BAND)
        sio2_fit = reparam_lorentz("sio2", C_BAND)
        assert cfg["resolved_si_reparam"] == [
            round(si_fit[0], 9), round(si_fit[1], 9), round(reparam_lam0(C_BAND), 9)]
        assert cfg["resolved_sio2_reparam"] == [
            round(sio2_fit[0], 9), round(sio2_fit[1], 9), round(reparam_lam0(C_BAND), 9)]


class TestReparamFit:
    def test_per_band_fit_is_far_better_than_the_dispersion_is_worth(self):
        """Fit error must stay small against the ~0.008-0.015 that n_Si
        actually moves across a band.

        The bound is 10%, not tighter, because the shipped lam0 deliberately
        trades fit quality for a better smoothed/unsmoothed split — lam0=10
        fits to ~2.4% of swing but fails the strip gate, lam0=2 fits to ~5.5%
        and passes it. Fit quality is the wrong thing to optimise here; this
        test only guards against it becoming badly wrong.
        """
        for band in (C_BAND, O_BAND):
            _, _, max_dn = reparam_lorentz("si", band)
            swing = abs(float(sellmeier_index("si", band.min())
                              - sellmeier_index("si", band.max())))
            assert max_dn < 0.10 * swing, (
                f"fit error {max_dn:.2e} is more than 10% of the {swing:.5f} "
                "dispersion it exists to capture"
            )

    def test_shipped_lam0_keeps_sigma_small_enough_to_be_unremarkable(self):
        """lam0=10 needed sigma ~17.6 (Si, C-band) and ~37.7 (O-band), and
        meep stability was never verified there. The shipped pole keeps sigma
        O(1), which retires that risk rather than arguing about it."""
        for band in (C_BAND, O_BAND):
            _, sigma, _ = reparam_lorentz("si", band)
            assert abs(sigma) < 5.0, f"sigma {sigma:.3f} is back in the risky range"

    def test_band_approaching_the_pole_is_refused(self):
        """The fitted term is near-resonant as the band approaches lam0, so it
        must raise rather than return a model that is quietly meaningless.

        Only reachable with an EXPLICIT lam0: the derived pole is 1.15x the
        band's top edge and the guard trips at 1.10x, so auto-derived fits
        clear their own guard by construction. That is the intended behaviour,
        and this pins the guard for callers who override it.
        """
        with pytest.raises(ValueError, match="too close to the fitted pole"):
            reparam_lorentz("si", C_BAND, lam0=1.65)

    def test_derived_lam0_always_clears_its_own_guard(self):
        """The auto rule must never produce a fit it would itself refuse —
        including on the widest band the examples span (1.26-1.70)."""
        for band in (C_BAND, O_BAND, np.linspace(1.26, 1.70, 21)):
            reparam_lorentz("si", band)          # must not raise
            assert reparam_lam0(band) > float(band.max()) * 1.10

    def test_lam0_tracks_the_band(self):
        """A single global pole cannot serve both bands — SOI_160/O-band needs
        <= 1.60 to clear the strip gate while the C-band's own guard floor is
        1.76. The pole must therefore move with the band."""
        assert reparam_lam0(O_BAND) < reparam_lam0(C_BAND)
        assert reparam_lam0(O_BAND) < 1.60

    def test_global_fit_is_worse_than_the_error_being_fixed(self):
        """Finding 1: this cannot be fit once across the full example range,
        which is *why* the media became a per-run quantity."""
        _, _, full = reparam_lorentz("si", np.linspace(1.26, 1.70, 201))
        _, _, cband = reparam_lorentz("si", C_BAND)
        assert full > 1e-3, "full-range fit unexpectedly good -- re-check Finding 1"
        assert full > 10 * cband

    def test_lam0_is_degenerate_so_it_must_be_pinned(self):
        """Finding 2: fit quality barely moves with lam0 while sigma swings by
        two orders of magnitude -- only sigma/lam0^2 is determined, so a free
        optimiser would walk the ridge to its bounds."""
        out = [reparam_lorentz("si", C_BAND, lam0=l) for l in (5.0, 10.0, 20.0, 50.0)]
        dns = [o[2] for o in out]
        sigmas = [o[1] for o in out]
        assert max(dns) - min(dns) < 5e-5, "lam0 changes fit quality more than expected"
        assert max(sigmas) / min(sigmas) > 50, "sigma should trade off strongly against lam0"

    def test_single_wavelength_is_rejected(self):
        with pytest.raises(ValueError, match="band, not a single point"):
            reparam_lorentz("si", 1.55)


class TestSmoothableContrastIsRestored:
    """The reverted attempt failed because eps_inf = 1 left meep nothing to
    subpixel-average. This is the regression guard for that specific failure."""

    def test_eps_inf_carries_essentially_all_the_contrast(self):
        e_si, _, _ = reparam_lorentz("si", C_BAND)
        e_ox, _, _ = reparam_lorentz("sio2", C_BAND)
        true_contrast = float(sellmeier_index("si", 1.55) ** 2
                              - sellmeier_index("sio2", 1.55) ** 2)
        frac = (e_si - e_ox) / true_contrast
        assert 0.95 < frac < 1.10, (
            f"only {frac:.1%} of the eps contrast is in eps_inf -- subpixel "
            "smoothing would again act on nearly nothing"
        )

    def test_built_medium_disperses_to_the_fit_accuracy(self):
        """meep's own evaluation of the Medium must track Sellmeier across the
        band, not just at one point -- otherwise nothing was gained.

        Bound is 1e-3 rather than tighter because the shipped pole trades fit
        accuracy for a better smoothed/unsmoothed split; what matters is that
        this stays far below the ~0.008 the frozen medium is wrong by at the
        band edges, which the next test pins directly.
        """
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        med = meep_medium_reparam("si", C_BAND)
        for wl in (1.50, 1.55, 1.60):
            n_meep = float(np.sqrt(med.epsilon(1.0 / wl)[0][0].real))
            assert abs(n_meep - float(sellmeier_index("si", wl))) < 1e-3

    def test_dispersive_beats_the_frozen_medium_by_a_wide_margin(self):
        """The actual value proposition, stated as a test rather than a claim:
        across the band the fitted medium must track n(λ) far better than a
        single index frozen at band centre — which is what ships today."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        for band in (C_BAND, O_BAND):
            med = meep_medium_reparam("si", band)
            frozen = float(sellmeier_index("si", float(np.mean(band))))
            worst_fit = worst_frozen = 0.0
            for wl in (float(band.min()), float(band.max())):
                true = float(sellmeier_index("si", wl))
                fit = float(np.sqrt(med.epsilon(1.0 / wl)[0][0].real))
                worst_fit = max(worst_fit, abs(fit - true))
                worst_frozen = max(worst_frozen, abs(frozen - true))
            assert worst_fit < 0.25 * worst_frozen, (
                f"dispersive error {worst_fit:.2e} is not clearly better than "
                f"the frozen medium's {worst_frozen:.2e}"
            )

    def test_medium_has_real_dielectric_contrast(self):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        med = meep_medium_reparam("si", C_BAND)
        assert med.epsilon_diag.x > 10.0, (
            "eps_inf collapsed toward 1 -- this is the reverted attempt's bug"
        )
        assert len(med.E_susceptibilities) == 1


class TestOptInIsInert:
    """Renamed in spirit 2026-08-01: dispersion is no longer opt-in. What these
    still pin is that the flag reaches exactly what it should and nothing
    else — in particular that a scalar wavelength stays frozen either way."""

    def test_the_default_stack_is_dispersive_over_a_band(self):
        """THE DEFAULT, pinned. Flipped 2026-08-01 (was
        `test_default_stack_is_not_dispersive_even_given_a_band`) once the
        strip-test gate passed on all four stack/band combinations. This is the
        one test that should fail if someone flips it back without meaning to."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        stack = MaterialStack(wg_height=0.22)
        assert stack.dispersive is True
        med = stack.si_medium(C_BAND)
        assert len(med.E_susceptibilities) == 1
        assert med.epsilon_diag.x > 10.0, "eps_inf must still carry the contrast"

    def test_an_explicitly_frozen_stack_is_still_available(self):
        """The A/B path has to survive the default flip, or the frozen arm of
        every comparison in material_model_validation.md becomes unreachable."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        stack = MaterialStack(wg_height=0.22, dispersive=False)
        med = stack.si_medium(C_BAND)
        assert not med.E_susceptibilities
        assert abs(med.epsilon_diag.x - stack.si_index(float(np.mean(C_BAND))) ** 2) < 1e-9

    def test_scalar_stays_frozen_even_when_dispersive(self):
        """Single-wavelength callers (mode solves) want exactly the index at
        that wavelength; the flag must not reach them."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        stack = MaterialStack(wg_height=0.22, dispersive=True)
        assert not stack.si_medium(1.55).E_susceptibilities

    def test_band_is_dispersive_when_enabled(self):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        stack = MaterialStack(wg_height=0.22, dispersive=True)
        assert len(stack.si_medium(C_BAND).E_susceptibilities) == 1
        assert len(stack.sio2_medium(C_BAND).E_susceptibilities) == 1


class TestCacheFingerprint:
    """Work item 3 -- the sharp one. The probe indices come from
    `sellmeier_index()`, which the re-parameterisation does not touch, so
    without an explicit record the key would be identical with dispersion on
    and every stale entry would keep matching."""

    def test_flipping_dispersion_moves_the_key(self):
        off = resolved_material_config(
            MaterialStack(wg_height=0.22, dispersive=False), C_BAND)
        on = resolved_material_config(
            MaterialStack(wg_height=0.22, dispersive=True), C_BAND)
        assert off != on, (
            "cache key unchanged by enabling dispersion -- every existing "
            "entry would be served the old physics"
        )

    def test_key_records_the_fit_not_a_version_tag(self):
        cfg = resolved_material_config(
            MaterialStack(wg_height=0.22, dispersive=True), C_BAND)
        eps_inf, sigma, lam0 = cfg["resolved_si_reparam"]
        fit = reparam_lorentz("si", C_BAND)
        assert eps_inf == pytest.approx(fit[0], abs=1e-9)
        assert sigma == pytest.approx(fit[1], abs=1e-9)
        assert lam0 == pytest.approx(reparam_lam0(C_BAND))

    def test_different_bands_give_different_keys(self):
        """The media are per-run now, so two bands must not collide."""
        stack = MaterialStack(wg_height=0.22, dispersive=True)
        assert (resolved_material_config(stack, C_BAND)
                != resolved_material_config(stack, O_BAND))

    def test_non_dispersive_key_is_unchanged_by_this_work(self):
        """A frozen run's key must stay exactly what it always was, so entries
        cached before any of this work keep matching when frozen is asked for
        explicitly. (Before 2026-08-01 this was also the default; it is not any
        more, hence the explicit flag.)"""
        cfg = resolved_material_config(
            MaterialStack(wg_height=0.22, dispersive=False), C_BAND)
        assert set(cfg) == {
            "resolved_probe_wl", "resolved_si_index", "resolved_sio2_index"}

    def test_a_solver_override_moves_the_key_not_just_the_stack_flag(self):
        """REGRESSION, 2026-08-01. `resolved_material_config` read the stack's
        own flag, so a per-solver override never reached the fingerprint:
        `FDTDSolver(material_stack=<frozen>, dispersive=True)` built a medium
        with a real Lorentzian susceptibility and produced a byte-identical key,
        colliding with frozen-physics entries. Same defect the class docstring
        above describes, one level up."""
        frozen = MaterialStack(wg_height=0.22, dispersive=False)
        as_stack_says = resolved_material_config(frozen, C_BAND)
        as_solver_resolved = resolved_material_config(
            frozen, C_BAND, frozen.resolve_dispersive(True))
        assert as_stack_says != as_solver_resolved
        assert "resolved_dispersive" in as_solver_resolved


class TestSolverFlags:
    """`dispersive=` on FDTDSolver / EMESolver, and the deliberate asymmetry
    in what `None` means on each."""

    def test_fdtd_none_inherits_the_stack(self):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        from beef import FDTDSolver

        off = FDTDSolver(res=10, material_stack=MaterialStack(wg_height=0.22,
                                                             dispersive=False))
        on = FDTDSolver(res=10, material_stack=MaterialStack(wg_height=0.22,
                                                            dispersive=True))
        assert not off._media(C_BAND)[0].E_susceptibilities
        assert on._media(C_BAND)[0].E_susceptibilities

    def test_fdtd_default_stack_gives_a_dispersive_solver(self):
        """The deployment, end to end: build an FDTDSolver the way the examples
        do and get dispersive media without asking for them."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        from beef import FDTDSolver

        assert FDTDSolver(res=10)._media(C_BAND)[0].E_susceptibilities

    def test_fdtd_flag_overrides_the_stack_without_mutating_it(self):
        """The SOI_* stacks are process-wide singletons; a solver-level
        override must not leak into other solvers sharing one."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        from beef import FDTDSolver

        stack = MaterialStack(wg_height=0.22, dispersive=False)
        solver = FDTDSolver(res=10, material_stack=stack, dispersive=True)
        assert solver._media(C_BAND)[0].E_susceptibilities
        assert stack.dispersive is False
        assert not FDTDSolver(res=10, material_stack=stack)._media(
            C_BAND)[0].E_susceptibilities

    def test_eme_defaults_to_dispersive_even_on_a_non_dispersive_stack(self):
        """The regression the EME's non-inheritance exists to prevent: EME has
        always been exact per wavelength, so inheriting a frozen stack would
        silently downgrade every existing EME run. Still explicitly a *frozen*
        stack here — since 2026-08-01 the default stack is dispersive, so
        leaving it implicit would stop testing anything."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        from beef import EMESolver

        from beef.materials.stack import meow_material

        eme = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0],
                        material_stack=MaterialStack(wg_height=0.22,
                                                     dispersive=False))
        assert _si_material(eme._extrusion_rules(C_BAND)) is meow_material("si")

    def test_eme_explicit_false_freezes(self):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        from beef import EMESolver
        from beef.materials.stack import meow_material

        eme = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0],
                        dispersive=False)
        si = _si_material(eme._extrusion_rules(C_BAND))
        assert si is not meow_material("si")
        assert np.allclose(si.n, si.n[0]), "frozen material still varies with wl"

    def test_frozen_eme_material_is_constant_across_the_band(self):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        from beef.materials.stack import _flat_meow_material

        mat = _flat_meow_material("si", 1.31)
        assert float(np.real(mat.n[0])) == pytest.approx(
            float(sellmeier_index("si", 1.31)))
        assert np.allclose(mat.n, mat.n[0])


class TestHybridConsistencyGuard:
    def _hybrid(self, fdtd_disp=None, eme_disp=None):
        from beef import EMESolver, FDTDSolver, HybridSolver

        f = FDTDSolver(res=10, dispersive=fdtd_disp)
        e = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0],
                      dispersive=eme_disp)
        return HybridSolver(solvers=[f, e], boundaries=[0.0, 1.0, 2.0],
                            interface_modes=[2])

    def test_the_default_configuration_is_now_quiet(self, caplog):
        """INVERTED 2026-08-01, and this inversion is the whole point of the
        change. This test used to be `test_default_configuration_is_flagged`
        and asserted the guard fired on the shipped default, because it did: a
        frozen FDTD chained to an exact-per-wavelength EME, 0.00786 apart in
        n_Si at the O-band edges. Both engines now default dispersive, so the
        chain the examples build no longer trips its own guard."""
        with caplog.at_level("WARNING"):
            self._hybrid()._check_material_consistency(O_BAND)
        assert "disagree on material dispersion" not in caplog.text

    def test_a_deliberately_mixed_chain_is_still_flagged(self, caplog):
        """The guard must keep working — it just needs an actually-mixed chain
        now that the default is not one. This is the configuration that used to
        ship, constructed explicitly."""
        with caplog.at_level("WARNING"):
            self._hybrid(fdtd_disp=False, eme_disp=True)._check_material_consistency(O_BAND)
        assert "disagree on material dispersion" in caplog.text
        assert "0.00786" in caplog.text

    def test_both_dispersive_is_quiet(self, caplog):
        with caplog.at_level("WARNING"):
            self._hybrid(fdtd_disp=True, eme_disp=True)._check_material_consistency(O_BAND)
        assert "disagree on material dispersion" not in caplog.text

    def test_both_frozen_is_quiet(self, caplog):
        """The A/B configuration the EME flag exists for."""
        with caplog.at_level("WARNING"):
            self._hybrid(fdtd_disp=False, eme_disp=False)._check_material_consistency(O_BAND)
        assert "disagree on material dispersion" not in caplog.text

    def test_single_wavelength_is_not_flagged(self, caplog):
        """No band, nothing to disagree about."""
        with caplog.at_level("WARNING"):
            self._hybrid()._check_material_consistency(np.array([1.31]))
        assert "disagree on material dispersion" not in caplog.text


class TestMaterialIndexAgreement:
    """Work item 5's hard check. Two properties decide whether it is worth
    anything: it must compare the media actually built (not the shared index
    model, which agrees in every configuration), and it must use a tolerance
    (the dispersive FDTD is a *fit*; the EME stays exact)."""

    def _hybrid(self, fdtd_disp=None, eme_disp=None, action="warn"):
        from beef import EMESolver, FDTDSolver, HybridSolver

        stack = MaterialStack(wg_height=0.22)
        f = FDTDSolver(res=10, material_stack=stack, dispersive=fdtd_disp)
        e = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0],
                      material_stack=stack, dispersive=eme_disp)
        return HybridSolver(solvers=[f, e], boundaries=[0.0, 1.0, 2.0],
                            interface_modes=[2],
                            material_consistency_action=action)

    def test_the_configuration_that_used_to_ship_really_does_disagree(self):
        """The measured size of the gap the 2026-08-01 default flip closed. A
        frozen FDTD against a per-wavelength EME differs by ~0.008 in n_Si at
        the O-band edges. Constructed explicitly now — it was the default when
        this test was written (as `test_todays_default_really_does_disagree`)
        and is not any more."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        gap, wl, kind, a, b = self._hybrid(
            fdtd_disp=False, eme_disp=True)._material_disagreement(O_BAND)
        assert gap > 5e-3, f"expected the known ~0.008 mismatch, got {gap:.2e}"

    def test_the_shipped_default_now_passes_its_own_guard(self):
        """The deployment's acceptance criterion, stated as a test. Build the
        chain the way the examples do and the worst index gap must sit under
        the guard's own 2e-3 tolerance — it was 0.00786, over it, before."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        from beef.solvers.hybrid import _MATERIAL_INDEX_TOL

        gap, *_ = self._hybrid()._material_disagreement(O_BAND)
        assert gap < _MATERIAL_INDEX_TOL, (
            f"the shipped default still fails the material guard at {gap:.2e}")

    def test_enabling_dispersion_closes_it_to_the_fit_error(self):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        gap, *_ = self._hybrid(fdtd_disp=True, eme_disp=True)._material_disagreement(O_BAND)
        assert gap < 2e-3, f"dispersive chain still disagrees by {gap:.2e}"

    def test_both_frozen_agree_exactly(self):
        """Same frozen index on both sides -- the A/B configuration."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        gap, *_ = self._hybrid(fdtd_disp=False, eme_disp=False)._material_disagreement(O_BAND)
        assert gap < 1e-12

    def test_raise_action_rejects_the_broken_chain(self):
        """Explicitly mixed since 2026-08-01 — the default chain is no longer
        broken, so `self._hybrid(action="raise")` would not raise."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        with pytest.raises(ValueError, match="resolve different si indices"):
            self._hybrid(fdtd_disp=False, eme_disp=True, action="raise") \
                ._check_material_index_agreement(O_BAND)

    def test_raise_action_accepts_the_shipped_default(self):
        """Corollary worth its own name: `material_consistency_action="raise"`
        is now a usable setting for a default-constructed chain. It could not
        be before — that is why the default is "warn"."""
        self._hybrid(action="raise")._check_material_index_agreement(O_BAND)

    def test_raise_action_accepts_the_fixed_chain(self):
        self._hybrid(fdtd_disp=True, eme_disp=True, action="raise") \
            ._check_material_index_agreement(O_BAND)

    def test_ignore_action_is_silent(self, caplog):
        with caplog.at_level("WARNING"):
            self._hybrid(action="ignore")._check_material_index_agreement(O_BAND)
        assert "resolve different" not in caplog.text

    def test_it_does_not_compare_index_lookups(self):
        """The trap this check exists to avoid: `si_index()` is the shared
        Sellmeier model, so it agrees to machine precision in EVERY
        configuration, including the broken one. If the guard ever starts
        reading that instead of the built media, this fails."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        stack = MaterialStack(wg_height=0.22)
        assert abs(stack.si_index(1.26) - stack.si_index(1.26)) < 1e-15
        gap, *_ = self._hybrid()._material_disagreement(O_BAND)
        assert gap > 1e-6, (
            "guard is reading the shared index model, not the built media"
        )
