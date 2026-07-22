"""EMESolver: meow-based eigenmode expansion solver."""
from __future__ import annotations

import logging
import os
from typing import List, Optional

import numpy as np

from .._utils import extrude_polygons
from ..cache import (
    CacheMissError, build_run_config, build_whole_device_config, load_result, lookup,
    lookup_superset, save_run, update_family_pointer,
)
from ..geometry.base import DeviceGeometry
from ..materials.stack import MaterialStack, SOI_160
from ..results.smatrix import SMatrix

logger = logging.getLogger(__name__)


class EMESolver:
    """Eigenmode expansion solver wrapping meow.

    Args:
        num_modes: modes per EME region, left to right.
                   e.g. [1, 7, 15, 7, 1] for (wg, taper, MMI, taper, wg).
        num_cells: cells per region (same length as num_modes).
                   e.g. [1, 12, 1, 12, 1].  Regions with num_cells > 1 are
                   subdivided equally — use this for tapers.
        oxide_width: y-span of the cross-section mesh (µm).
        oxide_height: z-span of the cross-section mesh (µm).
        res: mesh resolution in pts/µm.
        material_stack: MaterialStack instance for extrusion rules.
        straight_length: access waveguide length on each side (µm).
        temperature: simulation temperature for meow Environment (°C).
        mode_solver: which meow mode solver compute_modes() dispatches to —
            "meep" (meow.compute_modes_meep, the default) or "tidy3d"
            (meow.compute_modes_tidy3d). "meep" matches the same underlying
            eigenmode solver/conventions FDTDSolver's own
            get_eigenmode_coefficients() calls use in HybridSolver chains,
            avoiding basis mismatches between EME and FDTD sections — use
            that for hybrid tutorials. A standalone EME-only solve (no FDTD
            sections to stay consistent with) is free to use "tidy3d"
            instead. Does not affect solve_sweep(), which always uses
            tidy3d's batched mode solver regardless of this setting.
    """

    def __init__(
        self,
        num_modes: List[int],
        num_cells: List[int],
        lengths: List[float],
        oxide_width: float = 5.0,
        oxide_height: float = 1.5,
        res: float = 40.0,
        material_stack: MaterialStack = SOI_160,
        temperature: float = 25.0,
        mode_solver: str = "meep",
    ):
        if len(num_modes) != len(num_cells):
            raise ValueError(
                f"num_modes ({len(num_modes)}) and num_cells ({len(num_cells)}) must be same length"
            )
        if len(lengths) != len(num_modes):
            raise ValueError(
                f"lengths ({len(lengths)}) must match num_modes ({len(num_modes)})"
            )
        if mode_solver not in ("meep", "tidy3d"):
            raise ValueError(f"unknown mode_solver {mode_solver!r}, expected 'meep' or 'tidy3d'")
        self.num_modes = num_modes
        self.num_cells = num_cells
        self.lengths = lengths
        self.straight_length = lengths[0]  # input access waveguide length
        self.oxide_width = oxide_width
        self.oxide_height = oxide_height
        self.res = res
        self.material_stack = material_stack
        self.temperature = temperature
        self.mode_solver = mode_solver

    def _build_mesh(self):
        import meow as mw

        return mw.Mesh2D(
            x=np.linspace(-self.oxide_width/2, self.oxide_width/2,
                          int(self.res * self.oxide_width) + 1),
            y=np.linspace(-self.oxide_height/2, self.oxide_height/2,
                          int(self.res * self.oxide_height) + 1),
        )

    def _build_cells(self, geometry: DeviceGeometry, side: str = "top", **geometry_kwargs):
        import meow as mw

        gds_polys = geometry.to_gds_eme(
            straight_length=self.straight_length,
            oxide_width=self.oxide_width,
            x_start=0.0,
            side=side,
            **geometry_kwargs,
        )
        structs = extrude_polygons(gds_polys, self.material_stack.extrusion_rules)

        Ls = []
        for region_len, n_cells in zip(self.lengths, self.num_cells):
            cell_len = region_len / n_cells
            Ls.extend([cell_len] * n_cells)

        mesh = self._build_mesh()
        return mw.create_cells(structures=structs, mesh=mesh, Ls=Ls)

    def _modes_per_cell(self) -> List[int]:
        """Flat list of mode counts, one entry per cell (not per region)."""
        out = []
        for n_modes, n_cells in zip(self.num_modes, self.num_cells):
            out.extend([n_modes] * n_cells)
        return out

    def _flat_cell_lengths(self, region_lengths: Optional[List[float]] = None) -> List[float]:
        """Flat per-cell length list (one entry per cell), from region lengths
        divided evenly across each region's num_cells. Defaults to
        self.lengths — the same Ls _build_cells() hands to mw.create_cells()."""
        lengths = self.lengths if region_lengths is None else region_lengths
        out: List[float] = []
        for region_len, n_cells in zip(lengths, self.num_cells):
            out.extend([region_len / n_cells] * n_cells)
        return out

    def _region_cell_span(self, region_idx: int) -> tuple:
        """(start, stop) indices of region_idx's cells within the flat cell list."""
        start = sum(self.num_cells[:region_idx])
        return start, start + self.num_cells[region_idx]

    def _infer_sweep_region(self) -> int:
        """The unique region index with num_cells > 1 — i.e. the subdivided
        taper region a length sweep varies. Raises if that's ambiguous."""
        variable = [i for i, n in enumerate(self.num_cells) if n > 1]
        if len(variable) != 1:
            raise ValueError(
                f"cannot infer sweep_region: expected exactly one region with "
                f"num_cells > 1, found indices {variable} (num_cells={self.num_cells}). "
                "Pass sweep_region explicitly."
            )
        return variable[0]

    def compute_modes(self, cross_section, num_modes: int) -> list:
        """Compute modes for a single cross-section (exposed for caching).

        Dispatches to meow's MEEP/MPB-backed solver or its tidy3d-backed one
        per self.mode_solver (see the class docstring) — "meep" matches the
        same underlying eigenmode solver/conventions as FDTDSolver's own
        get_eigenmode_coefficients() calls in HybridSolver chains, avoiding
        basis mismatches between EME and FDTD sections.
        """
        if self.mode_solver == "meep":
            from meow import compute_modes_meep
            return compute_modes_meep(cross_section, num_modes=num_modes)
        from meow import compute_modes_tidy3d
        return compute_modes_tidy3d(cross_section, num_modes=num_modes)

    def solve_section(
        self,
        geometry,
        wavelengths: np.ndarray,
        x0: float,
        x1: float,
        n_modes_in: int = 1,
        n_modes_out: int = 1,
        stub_left: float = 0.0,
        stub_right: float = 0.0,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        **geometry_kwargs,
    ) -> np.ndarray:
        """Modal transfer matrix for the device section [x0, x1] (device frame).

        Single-cell EME over the section — mode expansion naturally gives the
        full n_modes_out x n_modes_in block in one solve (no per-mode source
        runs needed), unlike FDTD.

        Args:
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache — see
                FDTDSolver.solve_section()'s save_dir docstring for the full
                semantics (identical here).
            cache_only: if True, never run the solve — a cache hit returns as
                usual, a miss raises CacheMissError. See
                FDTDSolver.solve_section()'s cache_only docstring for the
                full semantics (identical here).

        Returns:
            Complex ndarray [len(wavelengths), n_modes_out, n_modes_in].
        """
        import meow as mw

        wavelengths = np.asarray(wavelengths)
        n_modes = max(n_modes_in, n_modes_out)
        length = (x1 - x0) + stub_left + stub_right

        # A request for fewer modes than a prior call already computed can
        # reuse that larger result by slicing its mode axes — unlike FDTD,
        # EME has no excited/measured-side asymmetry and no monitor geometry
        # to pin, since compute_modes(cs, n_modes) solves the full
        # cross-section once regardless of how the n_modes_in/n_modes_out
        # split is drawn from it, so no extra precondition/invariant is
        # needed beyond dropping the mode counts from the family config
        # (Solver recommendation #2). Caveat: this assumes the top-k slice
        # of an N-mode solve matches a direct k-mode solve, which holds as
        # long as meow's mode solver doesn't reorder near-degenerate modes
        # differently between differently-sized requests — see
        # tests/test_cache_config.py's EME superset-reuse coverage.
        family_config = None
        if save_dir:
            family_config = build_run_config(
                self, "solve_section",
                dict(x0=x0, x1=x1, wavelengths=wavelengths,
                     stub_left=stub_left, stub_right=stub_right,
                     geometry_kwargs=geometry_kwargs),
                geometry,
            )
            config = build_run_config(
                self, "solve_section",
                dict(x0=x0, x1=x1, wavelengths=wavelengths,
                     n_modes_in=n_modes_in, n_modes_out=n_modes_out,
                     stub_left=stub_left, stub_right=stub_right,
                     geometry_kwargs=geometry_kwargs),
                geometry,
            )
            hit = lookup(save_dir, config)
            if hit is not None:
                logger.info("solve_section: cache hit at %s", hit)
                return load_result(hit)
            superset = lookup_superset(save_dir, family_config, (n_modes_out, n_modes_in))
            if superset is not None:
                folder, _stored_dims = superset
                logger.info("solve_section: superset cache hit at %s (sliced)", folder)
                return load_result(folder)[:, :n_modes_out, :n_modes_in]
        if cache_only:
            raise CacheMissError(
                f"solve_section(): cache_only=True but no cached run found for "
                f"this config (save_dir={save_dir!r}). Run once with "
                f"cache_only=False to populate the cache."
            )

        gds_polys = geometry.section_gds(
            x0, x1, self.oxide_width,
            stub_left=stub_left, stub_right=stub_right, **geometry_kwargs
        )
        structs = extrude_polygons(gds_polys, self.material_stack.extrusion_rules)
        mesh = self._build_mesh()
        cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])

        T = np.zeros((len(wavelengths), n_modes_out, n_modes_in), dtype=complex)

        for wl_idx, wl in enumerate(wavelengths):
            env = mw.Environment(wl=float(wl), T=self.temperature)
            css = [mw.CrossSection.from_cell(cell=cell, env=env) for cell in cells]
            modes = [self.compute_modes(cs, n_modes) for cs in css]
            S, port_map = mw.compute_s_matrix(modes, cells=cells)
            # S is a flat dense matrix keyed by string port labels ("left@n",
            # "right@n"), not an NxN block per port pair — look up indices
            # explicitly rather than assuming S[0,1] is the transmission block.
            for i in range(n_modes_out):
                for j in range(n_modes_in):
                    i_out = port_map[f"right@{i}"]
                    j_in = port_map[f"left@{j}"]
                    T[wl_idx, i, j] = S[i_out, j_in]

        if save_dir:
            folder = save_run(save_dir, config, T)
            logger.info("solve_section: cache saved at %s", folder)
            update_family_pointer(
                save_dir, family_config, (n_modes_out, n_modes_in), os.path.basename(folder),
            )

        return T

    def solve_section_length_sweep(
        self,
        geometry,
        wavelengths: np.ndarray,
        x0: float,
        x1: float,
        lengths: np.ndarray,
        n_modes_in: int = 1,
        n_modes_out: int = 1,
        stub_left: float = 0.0,
        stub_right: float = 0.0,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        **geometry_kwargs,
    ) -> np.ndarray:
        """Modal transfer matrices for [x0, x1] swept over many section
        lengths, from a single eigenmode solve per wavelength.

        [x0, x1] must be a uniform-cross-section region (the actual cell
        length used for meshing/mode-solving is irrelevant here — only the
        cross-section shape matters, and that doesn't change with length).
        Mode profiles and effective indices depend only on the cross-section
        and wavelength (meow's CrossSection.from_cell() never reads cell
        length), so for a single uniform cell the propagation S-matrix is
        exactly diagonal in the mode basis: S_ii(L) = exp(2j*pi*neff_i/wl*L).
        meow's compute_s_matrix() accepts an explicit cell_lengths override
        that swaps in a different L for the same already-computed modes, with
        no new mode solve — so this method costs exactly the same as one
        solve_section() call (one mode solve per wavelength), regardless of
        how many lengths are swept.

        Args:
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache — see
                FDTDSolver.solve_section()'s save_dir docstring for the full
                semantics (identical here; lengths is included in the config).
            cache_only: if True, never run the solve — a cache hit returns as
                usual, a miss raises CacheMissError. See
                FDTDSolver.solve_section()'s cache_only docstring for the
                full semantics (identical here).

        Returns:
            Complex ndarray [len(lengths), len(wavelengths), n_modes_out, n_modes_in].
        """
        import meow as mw

        wavelengths = np.asarray(wavelengths)
        lengths = np.asarray(lengths)
        n_modes = max(n_modes_in, n_modes_out)

        # See solve_section()'s comment on family_config/lookup_superset —
        # same reuse mechanism, same no-extra-precondition reasoning; lengths
        # stays part of the family key (unchanged, only mode counts dropped)
        # since it's a real physical input, not something mode-count reuse
        # can vary (Solver recommendation #2).
        family_config = None
        if save_dir:
            family_config = build_run_config(
                self, "solve_section_length_sweep",
                dict(x0=x0, x1=x1, wavelengths=wavelengths, lengths=lengths,
                     stub_left=stub_left, stub_right=stub_right,
                     geometry_kwargs=geometry_kwargs),
                geometry,
            )
            config = build_run_config(
                self, "solve_section_length_sweep",
                dict(x0=x0, x1=x1, wavelengths=wavelengths, lengths=lengths,
                     n_modes_in=n_modes_in, n_modes_out=n_modes_out,
                     stub_left=stub_left, stub_right=stub_right,
                     geometry_kwargs=geometry_kwargs),
                geometry,
            )
            hit = lookup(save_dir, config)
            if hit is not None:
                logger.info("solve_section_length_sweep: cache hit at %s", hit)
                return load_result(hit)
            superset = lookup_superset(save_dir, family_config, (n_modes_out, n_modes_in))
            if superset is not None:
                folder, _stored_dims = superset
                logger.info("solve_section_length_sweep: superset cache hit at %s (sliced)", folder)
                return load_result(folder)[:, :, :n_modes_out, :n_modes_in]
        if cache_only:
            raise CacheMissError(
                f"solve_section_length_sweep(): cache_only=True but no cached "
                f"run found for this config (save_dir={save_dir!r}). Run once "
                f"with cache_only=False to populate the cache."
            )

        gds_polys = geometry.section_gds(
            x0, x1, self.oxide_width,
            stub_left=stub_left, stub_right=stub_right, **geometry_kwargs
        )
        structs = extrude_polygons(gds_polys, self.material_stack.extrusion_rules)
        mesh = self._build_mesh()
        cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[1.0])  # length is a placeholder

        T = np.zeros((len(lengths), len(wavelengths), n_modes_out, n_modes_in), dtype=complex)

        for wl_idx, wl in enumerate(wavelengths):
            env = mw.Environment(wl=float(wl), T=self.temperature)
            css = [mw.CrossSection.from_cell(cell=cell, env=env) for cell in cells]
            modes = [self.compute_modes(cs, n_modes) for cs in css]  # one solve, reused below

            if wl_idx == len(wavelengths) // 2:
                # Diagnostic: confirm this cross-section actually resolves
                # n_modes distinct, non-degenerate supermodes (neff split
                # apart) rather than a near-degenerate/duplicate pair, which
                # would explain a single-arm launch failing to couple into
                # one of them.
                neffs = [m.neff for m in modes[0]]
                logger.debug("EME cross-section neff @ λ=%.3f µm: %s", wl, neffs)

            for L_idx, L in enumerate(lengths):
                S, port_map = mw.compute_s_matrix(modes, cells=cells, cell_lengths=[float(L)])
                for i in range(n_modes_out):
                    for j in range(n_modes_in):
                        i_out = port_map[f"right@{i}"]
                        j_in = port_map[f"left@{j}"]
                        T[L_idx, wl_idx, i, j] = S[i_out, j_in]

        if save_dir:
            folder = save_run(save_dir, config, T)
            logger.info("solve_section_length_sweep: cache saved at %s", folder)
            update_family_pointer(
                save_dir, family_config, (n_modes_out, n_modes_in), os.path.basename(folder),
            )

        return T

    def solve(self, geometry: DeviceGeometry, wavelengths: np.ndarray,
              ports: Optional[List[str]] = None, save_dir: Optional[str] = None,
              cache_only: bool = False,
              **kwargs) -> SMatrix:
        """Run EME simulation across all wavelengths.

        Args:
            geometry: device geometry implementing eme_lengths().
            wavelengths: [M] wavelength array (µm).
            ports: list of port sides to simulate.
                   Defaults to ["top", "bot"] for 2-output devices.
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache — see
                solve_section()'s save_dir docstring for the full semantics
                (identical here, fingerprinted via geometry.to_gds_eme()
                instead of section_gds()). Cached per port side, since each
                side builds its own geometry/cells.
            cache_only: if True, never run the solve — a cache hit for a
                given port returns as usual, a miss raises CacheMissError
                immediately (fail-fast — does not check the remaining
                ports first). See solve_section()'s cache_only docstring for
                the full semantics.
            **kwargs: forwarded to geometry.to_gds_eme() (e.g. output_idx to
                      restrict the drawn output branch to a single port).

        Returns:
            SMatrix with entries keyed by port name.
        """
        import meow as mw

        if ports is None:
            ports = ["top", "bot"]

        wavelengths = np.asarray(wavelengths)
        smat = SMatrix(wavelengths=wavelengths)
        modes_per_cell = self._modes_per_cell()

        for side in ports:
            config = None
            if save_dir:
                gds_polys = geometry.to_gds_eme(
                    straight_length=self.straight_length, oxide_width=self.oxide_width,
                    x_start=0.0, side=side, **kwargs,
                )
                config = build_whole_device_config(
                    self, "solve",
                    dict(side=side, wavelengths=wavelengths, kwargs=kwargs),
                    geometry, gds_polys,
                )
                hit = lookup(save_dir, config)
                if hit is not None:
                    logger.info("solve: cache hit at %s", hit)
                    smat.set(load_result(hit), port=side)
                    continue
            if cache_only:
                raise CacheMissError(
                    f"solve(): cache_only=True but no cached run found for "
                    f"port={side!r} (save_dir={save_dir!r}). Run once with "
                    f"cache_only=False to populate the cache."
                )

            cells = self._build_cells(geometry, side=side, **kwargs)
            amplitudes = np.zeros(len(wavelengths), dtype=complex)

            for wl_idx, wl in enumerate(wavelengths):
                env = mw.Environment(wl=float(wl), T=self.temperature)
                css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
                modes = [self.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]
                S, _ = mw.compute_s_matrix(modes, cells=cells)
                amplitudes[wl_idx] = S[0, 1]

            if save_dir:
                folder = save_run(save_dir, config, amplitudes)
                logger.info("solve: cache saved at %s", folder)

            smat.set(amplitudes, port=side)

        return smat

    def solve_length_sweep(
        self,
        geometry: DeviceGeometry,
        wavelengths: np.ndarray,
        lengths: np.ndarray,
        ports: Optional[List[str]] = None,
        sweep_region: Optional[int] = None,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        **kwargs,
    ) -> List[SMatrix]:
        """Whole-device S-matrices for many device lengths, from ONE mode
        solve per wavelength — the mode basis is reused across every length.

        The cross-section *shape* an EME cell samples depends only on the
        taper fraction (i+0.5)/num_cells the solver evaluates, never on a
        region's physical length: a linear taper at fraction t has the same
        width regardless of how long the taper is (only Δx-per-cell scales).
        meow's S-matrix factors cleanly into an interface part (depends only
        on the modes) and a propagation part (the only place a cell's
        physical length enters — as exp(2j·π·neff/wl·L_cell)); see
        meow.eme.compute_propagation_s_matrix. So sweeping one region's
        length changes only a cell_lengths list handed to the cascade, never
        the modes or interfaces — turning an N-length sweep from N mode
        solves per wavelength into one. This is the multi-cell generalization
        of solve_section_length_sweep()'s single-uniform-cell trick.

        Only a region whose cross-section scales with its length (a taper or a
        straight) may be swept. sweep_region selects which self.lengths entry
        varies; it defaults to the unique region with num_cells > 1 (the
        subdivided taper). num_cells for that region is held FIXED across the
        sweep — it's the cross-section-staircase resolution (an accuracy dial,
        checkable by raising it), independent of the physical length, which
        only sets each cell's propagation length lengths[k]/num_cells.

        The geometry is drawn once, at whatever reference length it was built
        with, purely to sample those (length-independent) cross-section
        shapes — its swept-region length does not affect the result; only the
        values in `lengths` do. Every other (non-swept) region keeps its
        self.lengths value, so those must already match the real device
        (e.g. fixed access-waveguide stubs).

        Args:
            geometry: device geometry, built at any reference length for the
                swept region.
            wavelengths: [M] wavelength array (µm).
            lengths: [K] swept-region lengths to evaluate (µm).
            ports: port sides to simulate (as in solve()); defaults to
                ["top", "bot"].
            sweep_region: index into self.lengths that varies (default: the
                unique num_cells > 1 region — see _infer_sweep_region()).
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache. The full [K, M] result
                array is cached per port side, keyed on the solver config, the
                reference geometry, `lengths`, and sweep_region.
            cache_only: if True, never run the solve — a cache hit returns as
                usual, a miss raises CacheMissError (fail-fast per port).
            **kwargs: forwarded to geometry.to_gds_eme() (e.g. output_idx).

        Returns:
            list[SMatrix] aligned with `lengths`: result[k] is the SMatrix for
            device length lengths[k], keyed by port name (same convention as
            solve()).
        """
        import meow as mw

        if ports is None:
            ports = ["top", "bot"]
        if sweep_region is None:
            sweep_region = self._infer_sweep_region()

        wavelengths = np.asarray(wavelengths)
        lengths = np.asarray(lengths)
        modes_per_cell = self._modes_per_cell()
        ref_cell_lengths = self._flat_cell_lengths()
        c_start, c_stop = self._region_cell_span(sweep_region)
        n_region_cells = c_stop - c_start

        results = [SMatrix(wavelengths=wavelengths) for _ in lengths]

        for side in ports:
            config = None
            if save_dir:
                gds_polys = geometry.to_gds_eme(
                    straight_length=self.straight_length, oxide_width=self.oxide_width,
                    x_start=0.0, side=side, **kwargs,
                )
                config = build_whole_device_config(
                    self, "solve_length_sweep",
                    dict(side=side, wavelengths=wavelengths, lengths=lengths,
                         sweep_region=sweep_region, kwargs=kwargs),
                    geometry, gds_polys,
                )
                hit = lookup(save_dir, config)
                if hit is not None:
                    logger.info("solve_length_sweep: cache hit at %s", hit)
                    amps = load_result(hit)
                    for k in range(len(lengths)):
                        results[k].set(amps[k], port=side)
                    continue
            if cache_only:
                raise CacheMissError(
                    f"solve_length_sweep(): cache_only=True but no cached run "
                    f"found for port={side!r} (save_dir={save_dir!r}). Run once "
                    f"with cache_only=False to populate the cache."
                )

            cells = self._build_cells(geometry, side=side, **kwargs)
            amps = np.zeros((len(lengths), len(wavelengths)), dtype=complex)

            for wl_idx, wl in enumerate(wavelengths):
                env = mw.Environment(wl=float(wl), T=self.temperature)
                css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
                modes = [self.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]
                for k, L in enumerate(lengths):
                    cell_lengths = list(ref_cell_lengths)
                    for j in range(c_start, c_stop):
                        cell_lengths[j] = float(L) / n_region_cells
                    S, _ = mw.compute_s_matrix(modes, cells=cells, cell_lengths=cell_lengths)
                    amps[k, wl_idx] = S[0, 1]

            if save_dir:
                folder = save_run(save_dir, config, amps)
                logger.info("solve_length_sweep: cache saved at %s", folder)

            for k in range(len(lengths)):
                results[k].set(amps[k], port=side)

        return results

    def solve_sweep(self, geometry: DeviceGeometry, wavelengths: np.ndarray,
                    ports: Optional[List[str]] = None,
                    wl_batch_size: int = 1, **kwargs) -> SMatrix:
        """Memory-efficient EME sweep using batched wavelength processing.

        The center region (index len(num_modes)//2) is treated as the
        variable region — modes are freed after each batch while left/right
        cell modes are pre-computed and reused.
        """
        import meow as mw
        from meow import compute_modes_tidy3d as _compute_modes

        if ports is None:
            ports = ["top", "bot"]

        wavelengths = np.asarray(wavelengths)
        M = len(wavelengths)
        smat = SMatrix(wavelengths=wavelengths)
        modes_per_cell = self._modes_per_cell()

        # Identify center region boundaries in the flat cell list
        center_region = len(self.num_modes) // 2
        cumcells = [0] + list(np.cumsum(self.num_cells))
        split_left = cumcells[center_region]
        split_right = cumcells[center_region + 1]

        cells_per_port = {side: self._build_cells(geometry, side=side) for side in ports}
        left_mode_counts = modes_per_cell[:split_left]
        center_mode_count = modes_per_cell[split_left]
        right_mode_counts = modes_per_cell[split_right:]

        left_cells = cells_per_port[ports[0]][:split_left]
        center_cell = cells_per_port[ports[0]][split_left]
        right_cells_per_port = {side: cells_per_port[side][split_right:] for side in ports}

        # Pre-compute left and right modes (constant across sweep)
        left_modes: dict = {}
        right_modes_per_port: dict = {side: {} for side in ports}

        for wl_idx, wl in enumerate(wavelengths):
            env = mw.Environment(wl=float(wl), T=self.temperature)
            left_modes[wl_idx] = [
                _compute_modes(mw.CrossSection.from_cell(cell=c, env=env), nm)
                for c, nm in zip(left_cells, left_mode_counts)
            ]
            for side in ports:
                right_modes_per_port[side][wl_idx] = [
                    _compute_modes(mw.CrossSection.from_cell(cell=c, env=env), nm)
                    for c, nm in zip(right_cells_per_port[side], right_mode_counts)
                ]

        all_wl = list(range(M))
        batches = [all_wl[i:i+wl_batch_size] for i in range(0, M, wl_batch_size)]
        acc: dict = {side: np.zeros(M, dtype=complex) for side in ports}

        for batch in batches:
            center_modes_batch = {}
            for wl_idx in batch:
                env = mw.Environment(wl=float(wavelengths[wl_idx]), T=self.temperature)
                cs = mw.CrossSection.from_cell(cell=center_cell, env=env)
                center_modes_batch[wl_idx] = [_compute_modes(cs, center_mode_count)]

            for wl_idx in batch:
                for side in ports:
                    full_modes = (left_modes[wl_idx]
                                  + center_modes_batch[wl_idx]
                                  + right_modes_per_port[side][wl_idx])
                    S, _ = mw.compute_s_matrix(full_modes, cells=cells_per_port[side])
                    acc[side][wl_idx] = S[0, 1]

            del center_modes_batch

        for side in ports:
            smat.set(acc[side], port=side)

        return smat
