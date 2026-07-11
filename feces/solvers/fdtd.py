"""FDTDSolver: meep-based FDTD solver.

Wraps the FDTD simulation logic from masters/mmi_2x2_master.py::_port_simulation()
and _run_fdtd().  MPI safety guards are consolidated here so calling code does not
need to repeat them.
"""
from __future__ import annotations

from typing import List, Optional

import numpy as np

from ..geometry.base import DeviceGeometry
from ..materials.stack import MaterialStack, SOI_160
from ..results.smatrix import SMatrix


class FDTDSolver:
    """MEEP-based FDTD solver for one port of an MMI device.

    Args:
        res: simulation resolution (pixels/µm).
        dpml: PML thickness (µm).
        waveguide_length: access waveguide length inside domain (µm).
        source_position: x offset from PML where the source is placed (µm).
        sm_monitor_position: x offset from PML for the single-mode monitor (µm).
        mm_monitor_position: x offset from PML for the multimode monitor (µm).
        oxide_width: simulation y-span (µm).
        oxide_height: simulation z-span / height for 2D slabs (µm).
        mm_sim_length: length of the MMI section captured in the domain (µm).
        num_modes_mmi: number of monitored modes in the MMI cross-section.
        material_stack: MaterialStack for material index lookup.
        decay_tol: DFT decay tolerance for early stopping.
        mpi_safe: wrap writes with mp.am_master() / mp.all_wait() guards.
    """

    def __init__(
        self,
        res: float = 40.0,
        dpml: float = 0.5,
        waveguide_length: float = 0.5,
        source_position: float = 0.0,
        sm_monitor_position: float = 0.25,
        mm_monitor_position: float = 0.25,
        oxide_width: float = 5.0,
        oxide_height: float = 1.5,
        mm_sim_length: float = 1.0,
        num_modes_mmi: int = 15,
        material_stack: MaterialStack = SOI_160,
        decay_tol: float = 1e-11,
        mpi_safe: bool = True,
    ):
        self.res = res
        self.dpml = dpml
        self.waveguide_length = waveguide_length
        self.source_position = source_position
        self.sm_monitor_position = sm_monitor_position
        self.mm_monitor_position = mm_monitor_position
        self.oxide_width = oxide_width
        self.oxide_height = oxide_height
        self.mm_sim_length = mm_sim_length
        self.num_modes_mmi = num_modes_mmi
        self.material_stack = material_stack
        self.decay_tol = decay_tol
        self.mpi_safe = mpi_safe

    def _am_master(self) -> bool:
        if not self.mpi_safe:
            return True
        try:
            import meep as mp
            return mp.am_master()
        except ImportError:
            return True

    def _all_wait(self) -> None:
        if not self.mpi_safe:
            return
        try:
            import meep as mp
            mp.all_wait()
        except ImportError:
            pass

    def _build_domain(self, geometry: DeviceGeometry, wavelengths: np.ndarray):
        """Compute domain size and return (Sx, Sy, wl_center, si_medium, sio2_medium)."""
        wl_center = float(np.mean(wavelengths))
        total_len = geometry.total_length
        Sx = total_len + 2 * self.waveguide_length + 2 * self.dpml
        Sy = self.oxide_width
        si_med = self.material_stack.si_medium(wl_center)
        sio2_med = self.material_stack.sio2_medium(wl_center)
        return Sx, Sy, wl_center, si_med, sio2_med

    def solve(self, geometry: DeviceGeometry, wavelengths: np.ndarray,
              port: str = "top", input_mode: int = 1,
              gds_dir: str = "gds", **kwargs) -> SMatrix:
        """Run FDTD simulation for one input port.

        Simulates the input-side sub-domain: access waveguide → taper → partial MMI body.
        By reciprocity, the output-side coefficients (α_out) are obtained by sourcing
        from each output port and measuring at the MMI monitor (same geometry, reversed).

        Args:
            port: "top" or "bot" — which input port to source.
            input_mode: 1-based waveguide mode index to excite.
            gds_dir: directory to write temporary GDS files.

        Returns:
            SMatrix with complex amplitudes at the MMI monitor cross-section.
        """
        import meep as mp
        import os

        wavelengths = np.asarray(wavelengths)
        Sx, Sy, wl_center, si_med, sio2_med = self._build_domain(geometry, wavelengths)
        fcen = 1.0 / wl_center
        fwidth = (1.0/float(wavelengths.min()) - 1.0/float(wavelengths.max())) * 0.5
        nfreq = len(wavelengths)

        os.makedirs(gds_dir, exist_ok=True)

        if self._am_master():
            gds_comp = geometry.to_gds_fdtd(
                domain_sx=Sx, dpml=self.dpml,
                waveguide_length=self.waveguide_length,
                side="input", name=f"fdtd_{port}",
            )
            gds_path = os.path.join(gds_dir, f"fdtd_{port}.gds")
            gds_comp.write_gds(gds_path)
        self._all_wait()

        gds_path = os.path.join(gds_dir, f"fdtd_{port}.gds")
        geometry_meep = mp.get_GDSII_prisms(si_med, gds_path, layer=2, zmin=0, zmax=0)

        wg_offset = getattr(geometry, "waveguide_offset", 0.0) / 2
        src_yo = wg_offset if port == "top" else -wg_offset
        src_x = -Sx/2 + self.dpml + self.source_position
        src_kpt = mp.Vector3(1, 0, 0)

        src = mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=fwidth),
            center=mp.Vector3(src_x, src_yo),
            size=mp.Vector3(0, Sy),
            eig_band=input_mode,
            eig_kpoint=src_kpt,
            eig_match_freq=True,
        )

        sm_mon_x = -Sx/2 + self.dpml + self.sm_monitor_position
        mm_mon_x = Sx/2 - self.dpml - self.mm_monitor_position

        sim = mp.Simulation(
            resolution=self.res,
            cell_size=mp.Vector3(Sx, Sy),
            boundary_layers=[mp.PML(self.dpml)],
            geometry=geometry_meep,
            sources=[src],
            default_material=sio2_med,
        )

        sm_flux = sim.add_mode_monitor(
            fcen, fwidth, nfreq,
            mp.FluxRegion(center=mp.Vector3(sm_mon_x, src_yo), size=mp.Vector3(0, Sy)),
        )
        mm_flux = sim.add_mode_monitor(
            fcen, fwidth, nfreq,
            mp.FluxRegion(center=mp.Vector3(mm_mon_x), size=mp.Vector3(0, Sy)),
        )

        sim.run(
            mp.dft_ldos(),
            until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol),
        )

        smat = SMatrix(wavelengths=wavelengths)
        for mode_n in range(1, self.num_modes_mmi + 1):
            coeff = sim.get_eigenmode_coefficients(mm_flux, [mode_n], eig_parity=mp.NO_PARITY)
            alpha = coeff.alpha[:, :, 0].flatten()  # [nfreq] forward amplitude
            smat.set(alpha, port=port, in_mode=input_mode, out_mode=mode_n)

        sim.reset_meep()
        return smat
