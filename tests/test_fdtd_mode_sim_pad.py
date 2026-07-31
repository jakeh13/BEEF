"""Fast, meep-free tests for FDTDSolver._build_section_mode_sim()'s `pad`
override (fdtd_vs_EME_mode_validation.md item 5: the body_stub < dpml
blocker).

meep is faked out entirely, same pattern as tests/test_duplicate_band_guard.py
-- this exercises the pad plumbing itself (what stub_left/stub_right reach
geometry.section_prisms(), what PML thickness/cell size reach mp.Simulation),
not a real eigensolve.
"""
from __future__ import annotations

import sys
import types

import pytest

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
    def __init__(self, index):
        self.index = index


class _FakePML:
    def __init__(self, thickness):
        self.thickness = thickness


class _FakeVolume:
    def __init__(self, center, size):
        self.center = center
        self.size = size


class _FakeSimulation:
    """Records the cell_size/boundary_layers it was built with; init_sim() is
    a no-op since this test never calls get_eigenmode()."""

    instances: list = []

    def __init__(self, resolution, cell_size, boundary_layers, geometry,
                 default_material, eps_averaging, subpixel_tol, subpixel_maxeval):
        self.resolution = resolution
        self.cell_size = cell_size
        self.boundary_layers = boundary_layers
        self.geometry = geometry
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
    """Records the stub_left/stub_right it was asked to draw."""

    def __init__(self):
        self.calls = []

    def section_prisms(self, x0, x1, si_med, stub_left, stub_right, wg_height, **kwargs):
        self.calls.append(dict(x0=x0, x1=x1, stub_left=stub_left, stub_right=stub_right))
        v0 = _FakeVector3(x0 - stub_left, -0.1, -0.1)
        v1 = _FakeVector3(x1 + stub_right, 0.1, 0.1)
        return [_FakePrism([v0, v1], height=0.22, material=si_med)]


class ThinStubGeometry(FakeGeometry):
    """A geometry whose straight run-out past x0/x1 is only 0.25 um on each
    side -- the coupler/MMI hybrid examples' real body_stub -- raising if
    asked to draw further than that, the same failure
    section_mode_fields()'s self-padded query hits in production."""

    def section_prisms(self, x0, x1, si_med, stub_left, stub_right, wg_height, **kwargs):
        if stub_left > 0.25 or stub_right > 0.25:
            raise ValueError(
                "cut reaches past this section's straight run-out into a bend/taper")
        return super().section_prisms(x0, x1, si_med, stub_left, stub_right, wg_height, **kwargs)


@pytest.fixture
def solver():
    return FDTDSolver(res=20.0, sim_width=2.0, sim_height=1.0, dpml=0.5)


class TestPadDefaultsToDpml:
    def test_default_pad_matches_self_dpml(self, solver):
        geometry = FakeGeometry()
        solver._build_section_mode_sim(
            geometry, wl_center=1.55, x0=-0.05, x1=0.05,
            monitor_center=None, monitor_size=None)
        call = geometry.calls[0]
        assert call["stub_left"] == pytest.approx(0.5)
        assert call["stub_right"] == pytest.approx(0.5)
        sim = _FakeSimulation.instances[-1]
        assert sim.cell_size.x == pytest.approx(0.1 + 2 * 0.5)
        assert sim.boundary_layers[0].thickness == pytest.approx(0.5)


class TestExplicitPadOverrides:
    def test_pad_overrides_stub_and_pml_thickness(self, solver):
        geometry = FakeGeometry()
        solver._build_section_mode_sim(
            geometry, wl_center=1.55, x0=-0.05, x1=0.05,
            monitor_center=None, monitor_size=None, pad=0.1)
        call = geometry.calls[0]
        assert call["stub_left"] == pytest.approx(0.1)
        assert call["stub_right"] == pytest.approx(0.1)
        sim = _FakeSimulation.instances[-1]
        assert sim.cell_size.x == pytest.approx(0.1 + 2 * 0.1)
        assert sim.boundary_layers[0].thickness == pytest.approx(0.1)

    def test_pad_fits_a_thin_interior_stub_that_dpml_would_not(self, solver):
        """The actual item-5 blocker: body_stub=0.25 against dpml=0.5 means
        the default-padded query needs 0.55 um of straight run each side but
        only 0.25 is available -- reproduced here via ThinStubGeometry. A
        smaller explicit pad (0.1, still comfortably under 0.25) succeeds
        where the default does not."""
        geometry = ThinStubGeometry()
        with pytest.raises(ValueError, match="straight run-out"):
            solver._build_section_mode_sim(
                geometry, wl_center=1.55, x0=-0.05, x1=0.05,
                monitor_center=None, monitor_size=None)  # default pad = dpml = 0.5
        solver._build_section_mode_sim(
            geometry, wl_center=1.55, x0=-0.05, x1=0.05,
            monitor_center=None, monitor_size=None, pad=0.1)  # fits


class TestPadThreadsThroughEntryPoints:
    """section_mode_neff()/section_mode_fields() must forward their own
    `pad` argument to _build_section_mode_sim() rather than silently keeping
    the old dpml-only behavior."""

    def test_section_mode_neff_forwards_pad(self, solver, monkeypatch):
        seen = {}
        orig = solver._build_section_mode_sim

        def spy(self_, geometry, wl_center, x0, x1, monitor_center, monitor_size,
                pad=None, **kwargs):
            seen["pad"] = pad
            return orig(geometry, wl_center, x0, x1, monitor_center, monitor_size,
                        pad=pad, **kwargs)

        monkeypatch.setattr(FDTDSolver, "_build_section_mode_sim", spy)
        geometry = FakeGeometry()

        class _StopIteration(Exception):
            pass

        # We only care that pad was forwarded, not about the eigensolve loop
        # itself (untouched by this change) -- bail out right after the sim
        # is built by making get_eigenmode explode.
        def _boom(*a, **kw):
            raise _StopIteration

        monkeypatch.setattr(_FakeSimulation, "get_eigenmode", _boom, raising=False)
        with pytest.raises(_StopIteration):
            solver.section_mode_neff(geometry, [1.55], -0.05, 0.05, 1, pad=0.1)
        assert seen["pad"] == 0.1

    def test_section_mode_fields_forwards_pad(self, solver, monkeypatch):
        seen = {}
        orig = solver._build_section_mode_sim

        def spy(self_, geometry, wl_center, x0, x1, monitor_center, monitor_size,
                pad=None, **kwargs):
            seen["pad"] = pad
            return orig(geometry, wl_center, x0, x1, monitor_center, monitor_size,
                        pad=pad, **kwargs)

        monkeypatch.setattr(FDTDSolver, "_build_section_mode_sim", spy)
        geometry = FakeGeometry()

        class _StopIteration(Exception):
            pass

        def _boom(*a, **kw):
            raise _StopIteration

        monkeypatch.setattr(_FakeSimulation, "get_eigenmode", _boom, raising=False)
        with pytest.raises(_StopIteration):
            solver.section_mode_fields(
                geometry, 1.55, -0.05, 0.05, 1, y_grid=[0.0], z_grid=[0.0], pad=0.1)
        assert seen["pad"] == 0.1
