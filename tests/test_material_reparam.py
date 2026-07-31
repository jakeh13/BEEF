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
    _REPARAM_LAM0_UM,
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
        """The whole justification: fit error must be small against the ~0.008
        that n_Si actually moves across a band."""
        for band in (C_BAND, O_BAND):
            _, _, max_dn = reparam_lorentz("si", band)
            swing = abs(float(sellmeier_index("si", band.min())
                              - sellmeier_index("si", band.max())))
            assert max_dn < 0.05 * swing, (
                f"fit error {max_dn:.2e} is more than 5% of the {swing:.5f} "
                "dispersion it exists to capture"
            )

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
        band, not just at one point -- otherwise nothing was gained."""
        med = meep_medium_reparam("si", C_BAND)
        for wl in (1.50, 1.55, 1.60):
            n_meep = float(np.sqrt(med.epsilon(1.0 / wl)[0][0].real))
            assert abs(n_meep - float(sellmeier_index("si", wl))) < 5e-4

    def test_medium_has_real_dielectric_contrast(self):
        med = meep_medium_reparam("si", C_BAND)
        assert med.epsilon_diag.x > 10.0, (
            "eps_inf collapsed toward 1 -- this is the reverted attempt's bug"
        )
        assert len(med.E_susceptibilities) == 1


class TestOptInIsInert:
    def test_default_stack_is_not_dispersive_even_given_a_band(self):
        stack = MaterialStack(wg_height=0.22)
        med = stack.si_medium(C_BAND)
        assert not med.E_susceptibilities
        assert abs(med.epsilon_diag.x - stack.si_index(float(np.mean(C_BAND))) ** 2) < 1e-9

    def test_scalar_stays_frozen_even_when_dispersive(self):
        """Single-wavelength callers (mode solves) want exactly the index at
        that wavelength; the flag must not reach them."""
        stack = MaterialStack(wg_height=0.22, dispersive=True)
        assert not stack.si_medium(1.55).E_susceptibilities

    def test_band_is_dispersive_when_enabled(self):
        stack = MaterialStack(wg_height=0.22, dispersive=True)
        assert len(stack.si_medium(C_BAND).E_susceptibilities) == 1
        assert len(stack.sio2_medium(C_BAND).E_susceptibilities) == 1


class TestCacheFingerprint:
    """Work item 3 -- the sharp one. The probe indices come from
    `sellmeier_index()`, which the re-parameterisation does not touch, so
    without an explicit record the key would be identical with dispersion on
    and every stale entry would keep matching."""

    def test_flipping_dispersion_moves_the_key(self):
        off = resolved_material_config(MaterialStack(wg_height=0.22), C_BAND)
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
        assert lam0 == pytest.approx(_REPARAM_LAM0_UM)

    def test_different_bands_give_different_keys(self):
        """The media are per-run now, so two bands must not collide."""
        stack = MaterialStack(wg_height=0.22, dispersive=True)
        assert (resolved_material_config(stack, C_BAND)
                != resolved_material_config(stack, O_BAND))

    def test_non_dispersive_key_is_unchanged_by_this_work(self):
        """Everything already cached must keep matching while the flag is off."""
        cfg = resolved_material_config(MaterialStack(wg_height=0.22), C_BAND)
        assert set(cfg) == {
            "resolved_probe_wl", "resolved_si_index", "resolved_sio2_index"}


class TestSolverFlags:
    """`dispersive=` on FDTDSolver / EMESolver, and the deliberate asymmetry
    in what `None` means on each."""

    def test_fdtd_none_inherits_the_stack(self):
        from beef import FDTDSolver

        off = FDTDSolver(res=10, material_stack=MaterialStack(wg_height=0.22))
        on = FDTDSolver(res=10, material_stack=MaterialStack(wg_height=0.22,
                                                            dispersive=True))
        assert not off._media(C_BAND)[0].E_susceptibilities
        assert on._media(C_BAND)[0].E_susceptibilities

    def test_fdtd_flag_overrides_the_stack_without_mutating_it(self):
        """The SOI_* stacks are process-wide singletons; a solver-level
        override must not leak into other solvers sharing one."""
        from beef import FDTDSolver

        stack = MaterialStack(wg_height=0.22)          # dispersive=False
        solver = FDTDSolver(res=10, material_stack=stack, dispersive=True)
        assert solver._media(C_BAND)[0].E_susceptibilities
        assert stack.dispersive is False
        assert not FDTDSolver(res=10, material_stack=stack)._media(
            C_BAND)[0].E_susceptibilities

    def test_eme_defaults_to_dispersive_even_on_a_non_dispersive_stack(self):
        """The regression this asymmetry exists to prevent: EME has always been
        exact per wavelength, so inheriting an off-by-default stack would
        silently downgrade every existing EME run."""
        from beef import EMESolver

        from beef.materials.stack import meow_material

        eme = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0],
                        material_stack=MaterialStack(wg_height=0.22))
        assert _si_material(eme._extrusion_rules(C_BAND)) is meow_material("si")

    def test_eme_explicit_false_freezes(self):
        from beef import EMESolver
        from beef.materials.stack import meow_material

        eme = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0],
                        dispersive=False)
        si = _si_material(eme._extrusion_rules(C_BAND))
        assert si is not meow_material("si")
        assert np.allclose(si.n, si.n[0]), "frozen material still varies with wl"

    def test_frozen_eme_material_is_constant_across_the_band(self):
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

    def test_default_configuration_is_flagged(self, caplog):
        """Today's default really is mismatched -- frozen FDTD against an
        exact-per-wavelength EME. The guard must say so rather than pass."""
        with caplog.at_level("WARNING"):
            self._hybrid()._check_material_consistency(O_BAND)
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
