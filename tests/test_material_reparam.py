"""Tests for the eps_inf re-parameterisation of the Sellmeier model.

`validation_markdown/material_model_validation.md` step 2 / work items 1-3.
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
    MaterialStack,
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
