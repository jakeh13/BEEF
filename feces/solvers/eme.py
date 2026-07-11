"""EMESolver: meow-based eigenmode expansion solver.

Core loop migrated from masters/mmi_2x2_master.py::eme_simulation() and
sweeps/mmi_2x2_sweep_parallel_efficient.py::_run_width().
"""
from __future__ import annotations

from typing import List, Optional

import numpy as np

from ..geometry.base import DeviceGeometry
from ..geometry.mmi import MMI2x2
from ..materials.stack import MaterialStack, SOI_160
from ..results.smatrix import SMatrix


class EMESolver:
    """Eigenmode expansion solver wrapping meow.

    Args:
        num_modes_mmi: number of modes to compute in the MMI body.
        num_modes_taper: number of modes in the taper transition cells.
        num_cells_taper: number of cells to sub-divide each taper into.
        oxide_width: y-span of the cross-section mesh (µm).
        oxide_height: z-span of the cross-section mesh (µm).
        res: mesh resolution in pts/µm.
        material_stack: MaterialStack instance for extrusion rules.
        straight_length: access waveguide length on each side (µm).
        temperature: simulation temperature for meow Environment (°C).
    """

    def __init__(
        self,
        num_modes_mmi: int = 15,
        num_modes_taper: int = 7,
        num_cells_taper: int = 15,
        oxide_width: float = 5.0,
        oxide_height: float = 1.5,
        res: float = 40.0,
        material_stack: MaterialStack = SOI_160,
        straight_length: float = 0.5,
        temperature: float = 25.0,
    ):
        self.num_modes_mmi = num_modes_mmi
        self.num_modes_taper = num_modes_taper
        self.num_cells_taper = num_cells_taper
        self.oxide_width = oxide_width
        self.oxide_height = oxide_height
        self.res = res
        self.material_stack = material_stack
        self.straight_length = straight_length
        self.temperature = temperature

    def _build_mesh(self):
        import meow as mw

        return mw.Mesh2D(
            x=np.linspace(-self.oxide_width/2, self.oxide_width/2,
                          int(self.res * self.oxide_width) + 1),
            y=np.linspace(-self.oxide_height/2, self.oxide_height/2,
                          int(self.res * self.oxide_height) + 1),
        )

    def _build_cells(self, geometry: MMI2x2, side: str = "top"):
        """Build meow cell list for given geometry and port side."""
        import meow as mw

        gds_comp = geometry.to_gds_eme(
            straight_length=self.straight_length,
            oxide_width=self.oxide_width,
            x_start=0.0,
            side=side,
        )
        structs = mw.extrude_gds(gds_comp, self.material_stack.extrusion_rules)
        taper_Ls = np.array(
            [geometry.taper_length / self.num_cells_taper] * self.num_cells_taper
        )
        Ls = np.concatenate((
            [self.straight_length],
            taper_Ls,
            [geometry.mm_length],
            taper_Ls,
            [self.straight_length],
        ))
        mesh = self._build_mesh()
        return mw.create_cells(structures=structs, mesh=mesh, Ls=Ls)

    def compute_modes(self, cross_section, num_modes: int) -> list:
        """Compute modes for a single cross-section (exposed for caching)."""
        from meow import compute_modes_tidy3d
        return compute_modes_tidy3d(cross_section, num_modes=num_modes)

    def _num_modes_list(self, geometry: DeviceGeometry) -> List[int]:
        N = self.num_cells_taper
        return [1] + [self.num_modes_taper]*N + [self.num_modes_mmi] + [self.num_modes_taper]*N + [1]

    def solve(self, geometry: DeviceGeometry, wavelengths: np.ndarray,
              ports: Optional[List[str]] = None, **kwargs) -> SMatrix:
        """Run EME simulation across all wavelengths.

        Args:
            geometry: MMI2x2 (or compatible) geometry object.
            wavelengths: [M] wavelength array (µm).
            ports: list of port sides to simulate, e.g. ["top", "bot"].
                   Defaults to ["top", "bot"] for MMI2x2.

        Returns:
            SMatrix with entries keyed by port name.
        """
        import meow as mw

        if ports is None:
            ports = ["top", "bot"]

        wavelengths = np.asarray(wavelengths)
        smat = SMatrix(wavelengths=wavelengths)

        for side in ports:
            cells = self._build_cells(geometry, side=side)
            num_modes_list = self._num_modes_list(geometry)
            # split indices: left wg + tapers | MMI | tapers + right wg
            split_left = 1 + self.num_cells_taper
            split_right = split_left + 1

            amplitudes = np.zeros(len(wavelengths), dtype=complex)

            for wl_idx, wl in enumerate(wavelengths):
                env = mw.Environment(wl=float(wl), T=self.temperature)
                css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
                modes = []
                for cs, nm in zip(css, num_modes_list):
                    modes.append(self.compute_modes(cs, nm))
                S, _ = mw.compute_s_matrix(modes, cells=cells)
                amplitudes[wl_idx] = S[0, 1]

            smat.set(amplitudes, port=side)

        return smat

    def solve_sweep(self, geometry: DeviceGeometry, wavelengths: np.ndarray,
                    ports: Optional[List[str]] = None,
                    wl_batch_size: int = 1, **kwargs) -> SMatrix:
        """Memory-efficient EME sweep using batched wavelength processing.

        Equivalent to the batched logic in _run_width() from the 2x2 sweep:
        MMI body modes are solved in batches of wl_batch_size to limit peak
        memory usage, while left/right cell modes are reused across batches.

        Args:
            wl_batch_size: max wavelengths held in memory simultaneously.
        """
        import meow as mw
        from meow import compute_modes_tidy3d as _compute_modes

        if ports is None:
            ports = ["top", "bot"]

        wavelengths = np.asarray(wavelengths)
        M = len(wavelengths)
        smat = SMatrix(wavelengths=wavelengths)

        # Build cells for first port (left cells are port-independent).
        cells_per_port = {side: self._build_cells(geometry, side=side) for side in ports}
        split_left = 1 + self.num_cells_taper
        split_right = split_left + 1

        left_cells = cells_per_port[ports[0]][:split_left]
        left_num_modes = [1] + [self.num_modes_taper]*self.num_cells_taper
        mmi_cell = cells_per_port[ports[0]][split_left]
        right_cells_per_port = {side: cells_per_port[side][split_right:] for side in ports}
        right_num_modes = [self.num_modes_taper]*self.num_cells_taper + [1]

        # Pre-compute left and right cell modes (constant across wavelength sweep).
        left_modes = {}
        right_modes_per_port: dict = {side: {} for side in ports}

        for wl_idx, wl in enumerate(wavelengths):
            env = mw.Environment(wl=float(wl), T=self.temperature)
            lm = []
            for cell, nm in zip(left_cells, left_num_modes):
                cs = mw.CrossSection.from_cell(cell=cell, env=env)
                lm.append(_compute_modes(cs, nm))
            left_modes[wl_idx] = lm

            for side in ports:
                rm = []
                for cell, nm in zip(right_cells_per_port[side], right_num_modes):
                    cs = mw.CrossSection.from_cell(cell=cell, env=env)
                    rm.append(_compute_modes(cs, nm))
                right_modes_per_port[side][wl_idx] = rm

        # Process wavelengths in batches — MMI modes freed after each batch.
        all_wl = list(range(M))
        batches = [all_wl[i:i+wl_batch_size] for i in range(0, M, wl_batch_size)]
        acc: dict = {side: np.zeros(M, dtype=complex) for side in ports}

        for batch in batches:
            mmi_modes_batch = {}
            for wl_idx in batch:
                env = mw.Environment(wl=float(wavelengths[wl_idx]), T=self.temperature)
                cs = mw.CrossSection.from_cell(cell=mmi_cell, env=env)
                mmi_modes_batch[wl_idx] = [_compute_modes(cs, self.num_modes_mmi)]

            for wl_idx in batch:
                for side in ports:
                    full_modes = (left_modes[wl_idx]
                                  + mmi_modes_batch[wl_idx]
                                  + right_modes_per_port[side][wl_idx])
                    S, _ = mw.compute_s_matrix(full_modes, cells=cells_per_port[side])
                    acc[side][wl_idx] = S[0, 1]

            del mmi_modes_batch

        for side in ports:
            smat.set(acc[side], port=side)

        return smat
