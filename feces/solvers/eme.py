"""EMESolver: meow-based eigenmode expansion solver."""
from __future__ import annotations

from typing import List, Optional

import numpy as np

from ..geometry.base import DeviceGeometry
from ..materials.stack import MaterialStack, SOI_160
from ..results.smatrix import SMatrix


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
    ):
        if len(num_modes) != len(num_cells):
            raise ValueError(
                f"num_modes ({len(num_modes)}) and num_cells ({len(num_cells)}) must be same length"
            )
        if len(lengths) != len(num_modes):
            raise ValueError(
                f"lengths ({len(lengths)}) must match num_modes ({len(num_modes)})"
            )
        self.num_modes = num_modes
        self.num_cells = num_cells
        self.lengths = lengths
        self.straight_length = lengths[0]  # input access waveguide length
        self.oxide_width = oxide_width
        self.oxide_height = oxide_height
        self.res = res
        self.material_stack = material_stack
        self.temperature = temperature

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

        gds_comp = geometry.to_gds_eme(
            straight_length=self.straight_length,
            oxide_width=self.oxide_width,
            x_start=0.0,
            side=side,
            **geometry_kwargs,
        )
        structs = mw.extrude_gds(gds_comp, self.material_stack.extrusion_rules)

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

    def compute_modes(self, cross_section, num_modes: int) -> list:
        """Compute modes for a single cross-section (exposed for caching).

        Uses meow's MEEP/MPB-backed mode solver (not the tidy3d default) so
        EME sections use the same underlying eigenmode solver/conventions as
        FDTDSolver's own get_eigenmode_coefficients() calls in HybridSolver
        chains — avoids basis mismatches between EME and FDTD sections.
        """
        from meow import compute_modes_meep
        return compute_modes_meep(cross_section, num_modes=num_modes)

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
        **geometry_kwargs,
    ) -> np.ndarray:
        """Modal transfer matrix for the device section [x0, x1] (device frame).

        Single-cell EME over the section — mode expansion naturally gives the
        full n_modes_out x n_modes_in block in one solve (no per-mode source
        runs needed), unlike FDTD.

        Returns:
            Complex ndarray [len(wavelengths), n_modes_out, n_modes_in].
        """
        import meow as mw

        wavelengths = np.asarray(wavelengths)
        n_modes = max(n_modes_in, n_modes_out)
        length = (x1 - x0) + stub_left + stub_right

        gds_comp = geometry.section_gds(
            x0, x1, self.oxide_width,
            stub_left=stub_left, stub_right=stub_right, **geometry_kwargs
        )
        structs = mw.extrude_gds(gds_comp, self.material_stack.extrusion_rules)
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

        Returns:
            Complex ndarray [len(lengths), len(wavelengths), n_modes_out, n_modes_in].
        """
        import meow as mw

        wavelengths = np.asarray(wavelengths)
        lengths = np.asarray(lengths)
        n_modes = max(n_modes_in, n_modes_out)

        gds_comp = geometry.section_gds(
            x0, x1, self.oxide_width,
            stub_left=stub_left, stub_right=stub_right, **geometry_kwargs
        )
        structs = mw.extrude_gds(gds_comp, self.material_stack.extrusion_rules)
        mesh = self._build_mesh()
        cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[1.0])  # length is a placeholder

        T = np.zeros((len(lengths), len(wavelengths), n_modes_out, n_modes_in), dtype=complex)

        for wl_idx, wl in enumerate(wavelengths):
            env = mw.Environment(wl=float(wl), T=self.temperature)
            css = [mw.CrossSection.from_cell(cell=cell, env=env) for cell in cells]
            modes = [self.compute_modes(cs, n_modes) for cs in css]  # one solve, reused below

            for L_idx, L in enumerate(lengths):
                S, port_map = mw.compute_s_matrix(modes, cells=cells, cell_lengths=[float(L)])
                for i in range(n_modes_out):
                    for j in range(n_modes_in):
                        i_out = port_map[f"right@{i}"]
                        j_in = port_map[f"left@{j}"]
                        T[L_idx, wl_idx, i, j] = S[i_out, j_in]

        return T

    def solve(self, geometry: DeviceGeometry, wavelengths: np.ndarray,
              ports: Optional[List[str]] = None, **kwargs) -> SMatrix:
        """Run EME simulation across all wavelengths.

        Args:
            geometry: device geometry implementing eme_lengths().
            wavelengths: [M] wavelength array (µm).
            ports: list of port sides to simulate.
                   Defaults to ["top", "bot"] for 2-output devices.
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
            cells = self._build_cells(geometry, side=side, **kwargs)
            amplitudes = np.zeros(len(wavelengths), dtype=complex)

            for wl_idx, wl in enumerate(wavelengths):
                env = mw.Environment(wl=float(wl), T=self.temperature)
                css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
                modes = [self.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]
                S, _ = mw.compute_s_matrix(modes, cells=cells)
                amplitudes[wl_idx] = S[0, 1]

            smat.set(amplitudes, port=side)

        return smat

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
