"""Fast, meep-free tests for which MEDIUM `_build_section_mode_sim()` builds.

The regression these exist for (found 2026-07-31): the method built its media
with `si_medium(wl_center)` directly, bypassing `FDTDSolver._media()` and so
never reading `self.dispersive`. `section_mode_neff()` accepts a whole band and
loops it on one structure, so under `dispersive=True` it silently returned
results frozen at band centre -- measurably identical to `dispersive=False`
(max |dn_eff| = 0.000e+00), which made job 11586986's dispersive column a
byte-for-byte copy of its non-dispersive one. Nothing in the suite noticed,
because nothing asserted on the medium.

Both halves of the contract are pinned here:

    a band caller under dispersive=True gets per-wavelength media;
    every scalar caller stays byte-identical, in either configuration.

meep is faked out entirely, same pattern as test_fdtd_mode_sim_pad.py -- this
is about which Medium is constructed, not about an eigensolve.
"""
from __future__ import annotations

import sys
import types

import numpy as np
import pytest

from beef.materials import SOI_220
from beef.solvers.fdtd import FDTDSolver


class _FakeVector3:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x, self.y, self.z = x, y, z


class _FakePrism:
    def __init__(self, vertices, height, material):
        self.vertices = vertices
        self.height = height
        self.material = material


class _FakeMedium:
    """Frozen media arrive as `index=`; re-parameterised dispersive media as
    `epsilon=` plus a Lorentzian, which is how the tests below tell them
    apart."""

    def __init__(self, index=None, epsilon=None, E_susceptibilities=()):
        self.index = index
        self.epsilon_inf = epsilon
        self.E_susceptibilities = list(E_susceptibilities)

    @property
    def is_dispersive(self):
        return bool(self.E_susceptibilities)


class _FakeLorentzian:
    def __init__(self, frequency, gamma, sigma):
        self.frequency = frequency
        self.gamma = gamma
        self.sigma = sigma


class _FakePML:
    def __init__(self, thickness):
        self.thickness = thickness


class _FakeVolume:
    def __init__(self, center, size):
        self.center = center
        self.size = size


class _FakeSimulation:
    instances: list = []

    def __init__(self, resolution, cell_size, boundary_layers, geometry,
                 default_material, eps_averaging, subpixel_tol, subpixel_maxeval):
        self.geometry = geometry
        self.default_material = default_material
        _FakeSimulation.instances.append(self)

    def init_sim(self):
        pass

    def reset_meep(self):
        pass


@pytest.fixture(autouse=True)
def fake_meep(monkeypatch):
    fake = types.ModuleType("meep")
    fake.Vector3 = _FakeVector3
    fake.Prism = _FakePrism
    fake.Medium = _FakeMedium
    fake.LorentzianSusceptibility = _FakeLorentzian
    fake.PML = _FakePML
    fake.Volume = _FakeVolume
    fake.Simulation = _FakeSimulation
    fake.X = "X"
    fake.NO_PARITY = "NO_PARITY"
    fake.Ex, fake.Ey, fake.Ez = "Ex", "Ey", "Ez"
    fake.Hx, fake.Hy, fake.Hz = "Hx", "Hy", "Hz"
    monkeypatch.setitem(sys.modules, "meep", fake)
    _FakeSimulation.instances.clear()
    yield


class FakeGeometry:
    def section_prisms(self, x0, x1, si_med, stub_left, stub_right, wg_height, **kwargs):
        v0 = _FakeVector3(x0 - stub_left, -0.1, -0.1)
        v1 = _FakeVector3(x1 + stub_right, 0.1, 0.1)
        return [_FakePrism([v0, v1], height=0.22, material=si_med)]


BAND = np.array([1.50, 1.55, 1.60])


def _solver(dispersive):
    return FDTDSolver(res=20.0, sim_width=2.0, sim_height=1.0, dpml=0.5,
                      material_stack=SOI_220, dispersive=dispersive)


def _media_of(solver, **kwargs):
    """(si_medium, sio2_medium) that reached the built Simulation."""
    solver._build_section_mode_sim(
        FakeGeometry(), 1.55, -0.05, 0.05, None, None, **kwargs)
    sim = _FakeSimulation.instances[-1]
    return sim.geometry[0].material, sim.default_material


class TestScalarCallersUnchanged:
    """`section_mode_fields()` and every diagnostic script pass one wavelength.
    Those must be untouched by this change -- and not by convention:
    `MaterialStack._medium()` requires BOTH the flag and a real band
    (`size > 1 and min < max`), so a scalar cannot reach the dispersive path."""

    @pytest.mark.parametrize("dispersive", [False, True])
    def test_scalar_stays_frozen(self, dispersive):
        si, sio2 = _media_of(_solver(dispersive))
        assert not si.is_dispersive
        assert not sio2.is_dispersive
        assert si.index == pytest.approx(SOI_220.si_index(1.55))
        assert sio2.index == pytest.approx(SOI_220.sio2_index(1.55))

    def test_scalar_identical_across_the_flag(self):
        off_si, off_sio2 = _media_of(_solver(False))
        _FakeSimulation.instances.clear()
        on_si, on_sio2 = _media_of(_solver(True))
        assert on_si.index == off_si.index
        assert on_sio2.index == off_sio2.index


class TestBandCallers:
    def test_band_with_dispersion_on_is_dispersive(self):
        si, sio2 = _media_of(_solver(True), wavelengths=BAND)
        assert si.is_dispersive, (
            "a band caller under dispersive=True must get per-wavelength media; "
            "this is the exact regression that made job 11586986 vacuous")
        assert sio2.is_dispersive
        # eps_inf must carry the bulk of the permittivity -- that is the whole
        # point of the re-parameterisation (material_model_validation.md).
        assert si.epsilon_inf > 10.0
        assert si.E_susceptibilities[0].sigma > 0.0

    def test_band_with_dispersion_off_freezes_at_mean(self):
        si, sio2 = _media_of(_solver(False), wavelengths=BAND)
        assert not si.is_dispersive
        assert si.index == pytest.approx(SOI_220.si_index(float(np.mean(BAND))))
        assert sio2.index == pytest.approx(SOI_220.sio2_index(float(np.mean(BAND))))

    @pytest.mark.parametrize("band", [
        np.array([1.55]),          # one point is not a band
        np.array([1.55, 1.55]),    # nor is a degenerate one
    ])
    def test_degenerate_band_is_inert(self, band):
        si, _ = _media_of(_solver(True), wavelengths=band)
        assert not si.is_dispersive
        assert si.index == pytest.approx(SOI_220.si_index(1.55))


class _Bail(Exception):
    pass


def _stop_at_eigensolve(monkeypatch):
    def _boom(*a, **kw):
        raise _Bail
    monkeypatch.setattr(_FakeSimulation, "get_eigenmode", _boom, raising=False)


class TestEntryPointsForwardTheBand:
    """The regression test proper: `section_mode_neff()` must hand its whole
    band down, not just `mean(wavelengths)`. Everything above can be correct
    while this one line is wrong -- which is what happened."""

    def test_section_mode_neff_forwards_the_band(self, monkeypatch):
        solver = _solver(True)
        seen = {}
        orig = solver._build_section_mode_sim

        def spy(self_, geometry, wl_center, x0, x1, monitor_center, monitor_size,
                pad=None, wavelengths=None, **kwargs):
            seen["wavelengths"] = wavelengths
            seen["wl_center"] = wl_center
            return orig(geometry, wl_center, x0, x1, monitor_center, monitor_size,
                        pad=pad, wavelengths=wavelengths, **kwargs)

        monkeypatch.setattr(FDTDSolver, "_build_section_mode_sim", spy)
        _stop_at_eigensolve(monkeypatch)

        with pytest.raises(_Bail):
            solver.section_mode_neff(FakeGeometry(), BAND, -0.05, 0.05, 1)

        assert seen["wavelengths"] is not None, "band was dropped on the way down"
        np.testing.assert_allclose(seen["wavelengths"], BAND)
        assert seen["wl_center"] == pytest.approx(float(np.mean(BAND)))

    def test_section_mode_fields_stays_scalar(self, monkeypatch):
        """Deliberate, not an oversight: this entry point takes ONE wavelength,
        where the frozen index at that exact wavelength is the correct medium.
        Pinned so a future 'consistency' change does not invent a band here."""
        solver = _solver(True)
        seen = {}
        orig = solver._build_section_mode_sim

        def spy(self_, geometry, wl_center, x0, x1, monitor_center, monitor_size,
                pad=None, wavelengths=None, **kwargs):
            seen["wavelengths"] = wavelengths
            return orig(geometry, wl_center, x0, x1, monitor_center, monitor_size,
                        pad=pad, wavelengths=wavelengths, **kwargs)

        monkeypatch.setattr(FDTDSolver, "_build_section_mode_sim", spy)
        _stop_at_eigensolve(monkeypatch)

        with pytest.raises(_Bail):
            solver.section_mode_fields(
                FakeGeometry(), 1.55, -0.05, 0.05, 1,
                y_grid=[0.0], z_grid=[0.0])

        assert seen["wavelengths"] is None
