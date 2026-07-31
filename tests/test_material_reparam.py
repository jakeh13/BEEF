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
