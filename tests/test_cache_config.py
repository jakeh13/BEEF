"""Fast unit tests for beef.cache.config — no meep/meow import anywhere.
FDTDSolver()/EMESolver() are constructed (cheap, pure Python) but
solve_section()/solve() are never called — those do `import meep`/`import
meow` unconditionally, even down the cache_only=True short-circuit path, so
calling them here would break CI (pip-installed, no meep — see notes.md's
Testing section). The FDTDSolver.solve_section()-shaped resolved_layout
behavior is instead exercised via section_fdtd_layout() + build_run_config()
directly — the same two functions solve_section() itself calls, just without
going through the method that imports meep."""
from __future__ import annotations

import numpy as np
import pytest

import beef.cache.config as config
from beef.cache.config import _fingerprint, _solver_config
from beef.cache.store import _hash_config
from beef.solvers.eme import EMESolver
from beef.solvers.fdtd import (
    FDTDSolver, _resolved_spectral_config, section_fdtd_layout,
)


class _Fake:
    """Bare object standing in for a solver — vars(obj) is whatever's passed."""

    def __init__(self, **attrs):
        for key, value in attrs.items():
            setattr(self, key, value)


class TestMethodAttrsScoping:
    def test_unknown_pair_falls_back_to_everything(self):
        a = _Fake(x=1, y=2)
        b = _Fake(x=1, y=3)
        assert _solver_config(a, "some_unknown_method") != _solver_config(b, "some_unknown_method")

    def test_attribute_in_allowlist_changes_config(self, monkeypatch):
        monkeypatch.setitem(config._METHOD_ATTRS, ("_Fake", "solve_section"), ("x",))
        a = _Fake(x=1, y=100)
        b = _Fake(x=2, y=100)
        assert _solver_config(a, "solve_section") != _solver_config(b, "solve_section")

    def test_attribute_not_in_allowlist_is_ignored(self, monkeypatch):
        monkeypatch.setitem(config._METHOD_ATTRS, ("_Fake", "solve_section"), ("x",))
        a = _Fake(x=1, y=100)
        b = _Fake(x=1, y=999)  # y differs but isn't in the allowlist
        assert _solver_config(a, "solve_section") == _solver_config(b, "solve_section")

    def test_allowlist_is_per_method_not_per_class(self, monkeypatch):
        monkeypatch.setitem(config._METHOD_ATTRS, ("_Fake", "solve"), ("x",))
        # No entry for ("_Fake", "solve_section") -> falls back to everything,
        # so a "y" difference DOES show up there even though "solve" ignores it.
        a = _Fake(x=1, y=100)
        b = _Fake(x=1, y=999)
        assert _solver_config(a, "solve") == _solver_config(b, "solve")
        assert _solver_config(a, "solve_section") != _solver_config(b, "solve_section")


class TestRealSolverAllowlists:
    """Spot-checks the hand-traced FDTDSolver/EMESolver allowlists against
    the real classes — the thing most likely to silently drift as the
    solvers evolve."""

    def test_fdtd_solve_section_ignores_num_input_modes(self):
        a = FDTDSolver(num_input_modes=1)
        b = FDTDSolver(num_input_modes=99)  # solve_section() never reads this
        assert _solver_config(a, "solve_section") == _solver_config(b, "solve_section")

    def test_fdtd_solve_section_is_sensitive_to_res(self):
        a = FDTDSolver(res=20.0)
        b = FDTDSolver(res=40.0)  # solve_section() does read this
        assert _solver_config(a, "solve_section") != _solver_config(b, "solve_section")

    def test_fdtd_solve_is_sensitive_to_num_output_modes(self):
        a = FDTDSolver(num_output_modes=1)
        b = FDTDSolver(num_output_modes=4)  # solve() reads this, solve_section() doesn't
        assert _solver_config(a, "solve") != _solver_config(b, "solve")

    def test_eme_solve_section_ignores_num_modes_and_num_cells(self):
        a = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0])
        b = EMESolver(num_modes=[12], num_cells=[3], lengths=[1.0])  # only solve() reads these
        assert _solver_config(a, "solve_section") == _solver_config(b, "solve_section")

    def test_eme_solve_section_is_sensitive_to_oxide_width(self):
        a = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0], oxide_width=5.0)
        b = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0], oxide_width=6.0)
        assert _solver_config(a, "solve_section") != _solver_config(b, "solve_section")

    def test_eme_solve_is_sensitive_to_num_modes(self):
        a = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0])
        b = EMESolver(num_modes=[6], num_cells=[1], lengths=[1.0])
        assert _solver_config(a, "solve") != _solver_config(b, "solve")

    @pytest.mark.parametrize("method", ["solve", "solve_section"])
    @pytest.mark.parametrize("kwargs", [
        {"eps_averaging": False},
        {"subpixel_tol": 1e-6},
        {"subpixel_maxeval": 50},
    ])
    def test_fdtd_is_sensitive_to_subpixel_smoothing(self, method, kwargs):
        """These reach mp.Simulation and move n_eff at the second decimal
        (notes.md Solver #9). They were on meep's implicit defaults and
        invisible to the hash, so two solvers that would simulate differently
        shared one cache entry."""
        a = FDTDSolver()
        b = FDTDSolver(**kwargs)
        assert _solver_config(a, method) != _solver_config(b, method)

    @pytest.mark.parametrize("method", ["solve", "solve_section"])
    def test_subpixel_defaults_match_meep(self, method):
        """Defaults must be meep's own, or adding these knobs would silently
        change every existing result."""
        cfg = _solver_config(FDTDSolver(), method)
        assert cfg["eps_averaging"] is True
        assert cfg["subpixel_tol"] == 1e-4
        assert cfg["subpixel_maxeval"] == 100000


class _StackWithCacheKey:
    """Stand-in for a custom (non-MaterialStack) material stack that opts into
    cache fingerprinting via cache_key() — deliberately NOT JSON-safe on its
    own (has no __dict__ entry _json_safe could handle) to prove the hook,
    not a lucky fallback, is what includes it."""

    __slots__ = ("_params",)

    def __init__(self, **params):
        self._params = params

    def cache_key(self):
        return dict(self._params)


class _StackNoCacheKey:
    __slots__ = ("_params",)

    def __init__(self, **params):
        self._params = params


class TestCacheKeyHook:
    def test_cache_key_object_is_included_and_sensitive(self):
        a = _Fake(material_stack=_StackWithCacheKey(t=0.16))
        b = _Fake(material_stack=_StackWithCacheKey(t=0.30))
        assert _solver_config(a, "m") != _solver_config(b, "m")

    def test_cache_key_object_equal_when_params_equal(self):
        a = _Fake(material_stack=_StackWithCacheKey(t=0.16, s=4.5))
        b = _Fake(material_stack=_StackWithCacheKey(t=0.16, s=4.5))
        assert _solver_config(a, "m") == _solver_config(b, "m")

    def test_object_without_cache_key_is_dropped_with_warning(self, caplog):
        import logging
        a = _Fake(material_stack=_StackNoCacheKey(t=0.16))
        with caplog.at_level(logging.WARNING):
            cfg = _solver_config(a, "m")
        assert "material_stack" not in cfg
        assert any("cache_key" in r.message for r in caplog.records)


class TestSchemaVersion:
    def test_fingerprint_includes_schema_version(self):
        solver = _Fake(x=1)
        fp = _fingerprint(solver, "m", {}, "SomeGeometry", {})
        assert fp["cache_schema_version"] == config._CACHE_SCHEMA_VERSION

    def test_bumping_schema_version_changes_the_hash(self, monkeypatch):
        solver = _Fake(x=1)
        fp = _fingerprint(solver, "m", {}, "SomeGeometry", {})
        before = _hash_config(fp)

        monkeypatch.setattr(config, "_CACHE_SCHEMA_VERSION", config._CACHE_SCHEMA_VERSION + 1)
        fp2 = _fingerprint(solver, "m", {}, "SomeGeometry", {})
        after = _hash_config(fp2)

        assert before != after


class TestFDTDResolvedLayout:
    """FDTDSolver.solve_section() must fingerprint the fully-resolved
    SectionFDTDLayout, not just the monitor_center/monitor_size ingredients
    that derive it — otherwise a future change to that derivation logic
    can't be detected by the cache. Exercises section_fdtd_layout() +
    build_run_config() directly (the same two calls solve_section() itself
    makes before ever touching meep) instead of calling solve_section(),
    since that method's unconditional `import meep` would break CI here."""

    def _build_config(self, geometry, solver, monitor_center=None, monitor_size=None):
        layout = section_fdtd_layout(
            geometry, -2.0, 2.0, 0.0, 0.0,
            1, 2, solver.sim_width,
            solver.source_margin, solver.monitor_margin, solver.half_domain_aperture,
            {}, monitor_center=monitor_center, monitor_size=monitor_size,
        )
        from dataclasses import asdict
        return config.build_run_config(
            solver, "solve_section",
            dict(x0=-2.0, x1=2.0, wavelengths=[1.55],
                 n_modes_in=1, n_modes_out=2,
                 stub_left=0.0, stub_right=0.0,
                 monitor_center=monitor_center, monitor_size=monitor_size,
                 geometry_kwargs={}, resolved_layout=asdict(layout)),
            geometry,
        )

    def test_resolved_layout_is_included_in_the_cache_config(self):
        from beef.geometry.taper import Taper

        geometry = Taper(length=4.0, width_start=0.45, width_end=0.45)
        solver = FDTDSolver(res=20.0, sim_width=4.0)

        cfg = self._build_config(geometry, solver)
        layout = cfg["call_kwargs"]["resolved_layout"]
        for field in ("src_x", "norm_x", "mea_x", "src_c", "src_width", "mea_c", "mea_w", "kpoint"):
            assert field in layout

    def test_resolved_layout_differs_when_monitor_placement_differs(self):
        from beef.geometry.taper import Taper

        geometry = Taper(length=4.0, width_start=0.45, width_end=0.45)
        solver = FDTDSolver(res=20.0, sim_width=4.0)

        cfg_a = self._build_config(geometry, solver, monitor_center=0.0, monitor_size=4.0)
        cfg_b = self._build_config(geometry, solver, monitor_center=1.0, monitor_size=4.0)

        assert cfg_a["call_kwargs"]["resolved_layout"]["mea_c"] != \
            cfg_b["call_kwargs"]["resolved_layout"]["mea_c"]
        assert cfg_a != cfg_b


class TestResolvedSpectralConfig:
    """FDTDSolver fingerprints the *resolved* monitor frequency grid and
    source bandwidth, not just the `wavelengths` they derive from.

    Both are correct today without this — `wavelengths` is already hashed —
    so the tests below are about the failure mode that actually happened: the
    monitor grid was derived wrongly for months while every cache entry looked
    identical to a correct one. A derivation change with unchanged inputs must
    change the hash.
    """

    WLS = np.array([1.50, 1.55, 1.60])

    def test_grid_is_the_requested_wavelengths(self):
        r = _resolved_spectral_config(self.WLS, 2.0)
        assert np.allclose(sorted(1.0 / np.asarray(r["resolved_monitor_freqs"])),
                           sorted(self.WLS))

    def test_grid_is_ascending_in_frequency(self):
        freqs = np.asarray(_resolved_spectral_config(self.WLS, 2.0)["resolved_monitor_freqs"])
        assert np.all(np.diff(freqs) > 0), "meep requires an ascending frequency list"

    def test_fwidth_tracks_the_bandwidth_factor(self):
        a = _resolved_spectral_config(self.WLS, 2.0)["resolved_source_fwidth"]
        b = _resolved_spectral_config(self.WLS, 3.0)["resolved_source_fwidth"]
        assert b == pytest.approx(1.5 * a)

    def _config(self, extra=None):
        from beef.geometry.taper import Taper

        geometry = Taper(length=4.0, width_start=0.45, width_end=0.45)
        solver = FDTDSolver(res=20.0, sim_width=4.0)
        call_kwargs = dict(x0=-2.0, x1=2.0, wavelengths=list(self.WLS),
                           n_modes_in=1, n_modes_out=1)
        call_kwargs.update(_resolved_spectral_config(self.WLS, solver.source_bandwidth_factor))
        if extra:
            call_kwargs.update(extra)
        return config.build_run_config(solver, "solve_section", call_kwargs, geometry)

    def test_resolved_fields_reach_the_cache_config(self):
        ck = self._config()["call_kwargs"]
        assert "resolved_monitor_freqs" in ck
        assert "resolved_source_fwidth" in ck

    def test_a_reversed_grid_changes_the_hash(self):
        """The item-1 bug in miniature: same wavelengths, grid derived
        backwards. Pre-fix this was indistinguishable from a correct entry,
        which is why stale results were reused silently."""
        good = self._config()
        bad = self._config({"resolved_monitor_freqs":
                            list(reversed(good["call_kwargs"]["resolved_monitor_freqs"]))})
        assert _hash_config(good) != _hash_config(bad)

    def test_a_changed_fwidth_derivation_changes_the_hash(self):
        good = self._config()
        bad = self._config({"resolved_source_fwidth":
                            0.5 * good["call_kwargs"]["resolved_source_fwidth"]})
        assert _hash_config(good) != _hash_config(bad)

    def test_identical_calls_still_share_a_hash(self):
        """The point is to catch derivation changes, not to defeat caching."""
        assert _hash_config(self._config()) == _hash_config(self._config())
