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


class TestMaterialIndicesAreFingerprinted:
    """Until 2026-07-28 a MaterialStack reached the cache as {wg_height, sz}
    and nothing else, so the *optical* half of a material stack was invisible.
    Change the index table, or change how an index is looked up from it, and
    every cached entry kept matching and served the old physics.

    That is why directional_coupler_validation.md item 8 (the nearest-neighbour
    material lookup) could not be fixed first: the fix would have been silently
    ignored on every cached run and shown up only on a cold cache.
    """

    def test_resolved_indices_are_recorded(self):
        import numpy as np
        from beef.materials.stack import resolved_material_config
        from beef.materials import SOI_220

        cfg = resolved_material_config(SOI_220, np.linspace(1.5, 1.6, 11))
        # SOI_220 is dispersive by default since 2026-08-01, so the key also
        # carries the fitted re-parameterisation. The three index probes below
        # are what this test is actually about and are unchanged by that.
        assert {"resolved_probe_wl", "resolved_si_index",
                "resolved_sio2_index"} <= set(cfg)
        assert cfg["resolved_dispersive"] is True
        assert len(cfg["resolved_si_reparam"]) == 3
        # Band edges and centre, so the fingerprint sees dispersion (both media
        # are dispersive since item 8 step 2), not just one midband number.
        assert cfg["resolved_probe_wl"] == [1.5, 1.55, 1.6]
        assert all(3.0 < n < 4.0 for n in cfg["resolved_si_index"])
        assert all(1.3 < n < 1.6 for n in cfg["resolved_sio2_index"])

    def test_a_lookup_change_moves_the_fingerprint(self):
        """Snap -> interpolate is the change item 8 will actually make. It is
        invisible to a hash of the raw table (the data is untouched), which is
        why this records what the stack *returns*, not what it stores."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        import numpy as np
        import meow as mw
        from beef.materials.stack import resolved_material_config
        from beef.materials import SOI_220

        class _Interpolating(type(SOI_220)):
            def si_index(self, wl):
                return float(np.interp(wl, np.asarray(mw.silicon.params["wl"]),
                                       np.asarray(mw.silicon.n).real))

        wls = np.linspace(1.5, 1.6, 11)
        a = resolved_material_config(SOI_220, wls)
        b = resolved_material_config(_Interpolating(wg_height=0.22, sz=2.0), wls)
        assert a["resolved_si_index"] != b["resolved_si_index"], (
            "a snap-vs-interpolate change must move the cache fingerprint"
        )

    def test_geometry_alone_no_longer_determines_the_key(self):
        """The exact regression: two stacks agreeing on wg_height/sz but
        disagreeing on n must not produce the same fingerprint."""
        import numpy as np
        from beef.materials.stack import resolved_material_config
        from beef.materials import SOI_220

        class _Shifted(type(SOI_220)):
            def si_index(self, wl):
                return super().si_index(wl) + 1e-3

        a, b = SOI_220, _Shifted(wg_height=0.22, sz=2.0)
        assert (a.wg_height, a.sz) == (b.wg_height, b.sz)
        wls = np.linspace(1.5, 1.6, 11)
        assert resolved_material_config(a, wls) != resolved_material_config(b, wls)

    def test_fingerprinting_does_not_import_meep(self):
        """beef/cache/store.py::_is_master() dispatches on whether meep is in
        sys.modules and promises to import nothing itself. A MaterialStack
        cache_key() would have imported meow -> meep while a config was being
        built, flipping rank detection off its launcher-env fallback and
        electing every rank of an EME-only job as master -- the item 4 bug.
        Keeping the lookup out of the cache-config path is what avoids that,
        so no import of beef.cache may pull meep in."""
        import subprocess, sys
        out = subprocess.run(
            [sys.executable, "-c",
             "import sys; import beef.cache.config as c; "
             "print('meep' in sys.modules)"],
            capture_output=True, text=True,
        )
        assert out.stdout.strip().endswith("False"), out.stdout + out.stderr


class TestEMESolveFingerprintsMaterial:
    """`EMESolver.solve()`/`solve_length_sweep()` -- found and fixed 2026-08-03
    while implementing `mmi_1x4_sin300_validation.md` step 3.

    Every OTHER cache-relevant method (`FDTDSolver.solve`/`solve_section`,
    `EMESolver.solve_section`/`solve_section_length_sweep`) has merged
    `resolved_material_config(...)` into its cache config since 2026-08-01.
    `solve()`/`solve_length_sweep()` never did -- their fingerprint carried
    only `_solver_config()`'s bare `material_stack: {wg_height, sz}`. Two real
    consequences, not a hypothetical: (1) a stack differing only in
    `core_kind` at the same wg_height/sz -- e.g. `SIN_300` vs a same-thickness
    Si stack -- would have shared a cache entry, wrong material and all;
    (2) `dispersive` genuinely changes what these methods draw
    (`EMESolver._extrusion_rules()` frozen-vs-exact branch), and was not
    fingerprinted AT ALL for these two methods, meep/Si/SiO2 included -- this
    part of the bug predates nitride and could in principle already have bitten
    an `EMESolver(dispersive=False)` A/B run cached alongside a dispersive one
    at the same save_dir.

    Meep/meow-free per this file's own convention (see module docstring): does
    NOT call `.solve()` (its first line is `import meow as mw`, unconditional,
    even down the `cache_only=True` short-circuit). Instead confirms the fix
    is wired into the real method source, via `inspect` rather than a
    reimplementation. A real end-to-end confirmation (an actual `.solve()`
    call under `--runslow`) lives in `test_eme_length_sweep.py`.
    """

    def test_solve_source_calls_resolved_material_config(self):
        import inspect

        src = inspect.getsource(EMESolver.solve)
        assert "resolved_material_config(" in src
        assert "self._resolved_dispersive()" in src, (
            "must use the EME-specific dispersive resolution (None means "
            "dispersive here), not the stack's own flag -- see "
            "_resolved_dispersive()'s docstring"
        )

    def test_solve_length_sweep_source_calls_resolved_material_config(self):
        import inspect

        src = inspect.getsource(EMESolver.solve_length_sweep)
        assert "resolved_material_config(" in src
        assert "self._resolved_dispersive()" in src

    def test_solve_sweep_is_uncached_so_needs_no_fix(self):
        """`solve_sweep()` (the memory-batched variant) has no save_dir
        parameter at all -- confirms it is out of scope for this fix rather
        than silently unfixed."""
        import inspect

        sig = inspect.signature(EMESolver.solve_sweep)
        assert "save_dir" not in sig.parameters


class TestOneDispersionModelForBothEngines:
    """directional_coupler_validation.md item 8, steps 1 and 2.

    Two independent defects shared one root cause. The EME snapped n to the
    nearest sample of meow's table -- which has only *two* entries between 1.45
    and 1.65 um -- putting a 0.0063 step in the middle of the C-band. The FDTD
    froze one index at the band centre for a 100 nm broadband run, so it had no
    dispersion at all. Both engines now evaluate the same Sellmeier model.
    """

    WLS = (1.26, 1.31, 1.50, 1.55, 1.60, 1.70)

    def test_both_engines_agree_to_machine_precision(self):
        """The point of the exercise. EME goes through meow's material call,
        FDTD through a dispersive meep Medium; both must be the same function."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        import numpy as np
        from meow.environment import Environment
        from meow.fde.meep import meep_material
        from beef.materials.stack import sellmeier_index, meow_material, meep_medium

        for kind in ("si", "sio2"):
            for wl in self.WLS:
                want = float(sellmeier_index(kind, wl))
                eme, _ = meep_material(meow_material(kind), Environment(wl=wl))
                fdtd = float(np.sqrt(meep_medium(kind).epsilon(1.0 / wl)[0, 0].real))
                assert abs(eme - want) < 1e-12, (kind, wl, eme, want)
                assert abs(fdtd - want) < 1e-12, (kind, wl, fdtd, want)

    def test_the_c_band_kink_is_gone(self):
        """The specific symptom: n stepped 3.4784 -> 3.4710 between 1.56 and
        1.57 um. n(wl) must now be smooth and monotonically decreasing."""
        import numpy as np
        from beef.materials.stack import sellmeier_index

        wl = np.linspace(1.50, 1.60, 101)
        n = sellmeier_index("si", wl)
        d = np.diff(n)
        assert (d < 0).all(), "n(wl) must decrease across the C-band"
        # No step: the largest single increment is within 2x the mean one.
        assert abs(d).max() < 2 * abs(d).mean(), f"discontinuity in n(wl): {abs(d).max():.2e}"

    def test_the_fdtd_medium_stays_non_dispersive_on_purpose(self):
        """Item 8 step 2 made both media dispersive and it was reverted after
        measurement: meep subpixel-averages only the instantaneous eps, and the
        Sellmeier form (like meep's own library materials) puts eps_inf = 1 with
        the whole index in susceptibilities. The structure then has no
        dielectric contrast to smooth -- meep logs `epsilon diagonal = (1,1,1)`
        -- and on a 450x220 nm strip at res 40 that moved n_eff 2.313398 ->
        2.363038, +0.0496, against the ~0.002 the dispersion was worth.

        So `si_medium()` returns a plain index Medium, and this pins that. If
        someone re-enables dispersion, they must re-parameterise so eps_inf
        carries the bulk of the permittivity, and re-run the resolution sweeps
        -- notes.md recommendation #9."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        from beef.materials import SOI_220

        med = SOI_220.si_medium(1.55)
        assert not med.E_susceptibilities, (
            "si_medium() is dispersive again -- subpixel smoothing is now "
            "disabled at every material interface; see the docstring"
        )
        assert abs(med.epsilon_diag.x - SOI_220.si_index(1.55) ** 2) < 1e-9

    def test_the_dispersive_model_is_still_available_and_correct(self):
        """`meep_medium()` is unused by the solvers but kept as the starting
        point for the eps_inf re-parameterisation, so it must stay right."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        import numpy as np
        from beef.materials.stack import meep_medium, sellmeier_index

        for wl in (1.26, 1.55, 1.60):
            got = float(np.sqrt(meep_medium("si").epsilon(1.0 / wl)[0, 0].real))
            assert abs(got - float(sellmeier_index("si", wl))) < 1e-12

    def test_oxide_still_matches_meows_table_in_band(self):
        """Malitson either way, so step 2 must not have moved the cladding.
        (meep's `SiO2` is crystalline quartz and would have shifted it 0.0138.)"""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        import numpy as np
        import meow as mw
        from beef.materials.stack import sellmeier_index

        w = np.asarray(mw.silicon_oxide.params["wl"])
        n = np.asarray(mw.silicon_oxide.n).real
        for wl in (1.50, 1.55, 1.60):
            assert abs(float(sellmeier_index("sio2", wl)) - float(np.interp(wl, w, n))) < 1e-4

    def test_model_covers_every_example_band(self):
        """Salzberg is declared valid from 1.36 um but examples run the O-band
        to 1.26 um, so the media are built from raw susceptibilities with no
        valid_freq_range. Evaluating below 1.36 um must not raise."""
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        import numpy as np
        from beef.materials.stack import meep_medium

        for kind in ("si", "sio2"):
            for wl in (1.26, 1.36, 1.60):
                assert np.isfinite(meep_medium(kind).epsilon(1.0 / wl)[0, 0].real)

    def test_fingerprint_tracks_the_dispersion_not_just_one_point(self):
        """Both media are dispersive now, so a single centre index no longer
        describes the run. Two models agreeing at midband must still differ."""
        import numpy as np
        from beef.materials.stack import resolved_material_config
        from beef.materials import SOI_220

        class _FlatAtMidband(type(SOI_220)):
            def si_index(self, wl):
                return SOI_220.si_index(1.55)      # same centre, no slope

        wls = np.linspace(1.5, 1.6, 11)
        a = resolved_material_config(SOI_220, wls)
        b = resolved_material_config(_FlatAtMidband(wg_height=0.22, sz=2.0), wls)
        assert a["resolved_si_index"][1] == b["resolved_si_index"][1]   # agree at centre
        assert a != b, "fingerprint must see the dispersion, not just midband"
