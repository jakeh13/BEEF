"""HybridSolver: FDTD/EME combined solver.

Implements the formula from slide 4 of the FECES planning deck:

    α_full[m] = Σ_n  α_in[m,n] · α_mm[m,n] · α_out[m,n]
    T[m]       = |α_full[m]|²

where:
    α_in[m,n]  = input FDTD coupling coefficient at wavelength m, MMI mode n
    α_mm[m,n]  = EME propagation coefficient through multimode section
    α_out[m,n] = output FDTD coupling coefficient (by reciprocity)
"""
from __future__ import annotations

from typing import List, Optional

import numpy as np

from ..geometry.base import DeviceGeometry
from ..results.smatrix import SMatrix
from .eme import EMESolver
from .fdtd import FDTDSolver


class HybridSolver:
    """FDTD/EME hybrid solver.

    The input and output FDTD simulations run on partial domains (taper + partial MMI),
    while the EME solver handles the full multimode section.  The three sets of mode
    overlap coefficients are combined analytically.

    Args:
        input_fdtd: FDTDSolver configured for the input-side sub-domain.
        eme: EMESolver for the multimode section.
        output_fdtd: FDTDSolver configured for the output-side sub-domain.
    """

    def __init__(
        self,
        input_fdtd: FDTDSolver,
        eme: EMESolver,
        output_fdtd: FDTDSolver,
    ):
        self.input_fdtd = input_fdtd
        self.eme = eme
        self.output_fdtd = output_fdtd

    def solve(
        self,
        geometry: DeviceGeometry,
        wavelengths: np.ndarray,
        input_port: str = "top",
        input_mode: int = 1,
        output_ports: Optional[List[str]] = None,
        eme_ports: Optional[List[str]] = None,
        **kwargs,
    ) -> SMatrix:
        """Run FDTD/EME hybrid simulation.

        Steps:
          1. Input FDTD: source at input_port → measure MMI mode overlap coefficients α_in[m,n]
          2. EME: propagate through multimode section → α_mm[m,n]
          3. Output FDTD: source at each output port (reciprocity) → α_out[m,n]
          4. Combine: α_full[m] = Σ_n α_in[m,n] · α_mm[m,n] · α_out[m,n]

        Args:
            geometry: device geometry (must support to_gds_fdtd and to_gds_eme).
            wavelengths: [M] wavelength array (µm).
            input_port: "top" or "bot" — which input port to excite.
            input_mode: 1-based waveguide mode index at input.
            output_ports: list of output port labels to compute, e.g. ["top", "bot"].
            eme_ports: port sides for EME simulation; defaults to output_ports.

        Returns:
            SMatrix with T[m] = |α_full[m]|² for each output port.
        """
        wavelengths = np.asarray(wavelengths)
        if output_ports is None:
            output_ports = ["top", "bot"]
        if eme_ports is None:
            eme_ports = output_ports

        center_region = len(self.eme.num_modes) // 2
        N = self.eme.num_modes[center_region]

        # Step 1: input FDTD — α_in[m, n] for n = 1..N
        in_smat = self.input_fdtd.solve(
            geometry, wavelengths, port=input_port, input_mode=input_mode, **kwargs
        )
        # α_in shape: [M, N] complex
        alpha_in = np.stack([
            in_smat.coeff(port=input_port, in_mode=input_mode, out_mode=n)
            for n in range(1, N+1)
        ], axis=1)  # [M, N]

        # Step 2: EME — propagation coefficients α_mm[m, n] for each port side
        eme_smat = self.eme.solve(geometry, wavelengths, ports=eme_ports, **kwargs)
        # EME SMatrix: amplitude(port=side) = S[0,1] = the dominant coupling term.
        # For the hybrid formula we need the full modal decomposition.
        # Here we use a simplified form: α_mm[m,n] = S_n[m] from the n-th mode's contribution.
        # Full modal hybrid requires per-mode EME output; current EME returns S[0,1] only.
        # TODO: extend EMESolver.solve() to return per-mode amplitudes for full hybrid support.

        result = SMatrix(wavelengths=wavelengths)

        for out_port in output_ports:
            # Step 3: output FDTD — α_out[m, n] by reciprocity
            out_smat = self.output_fdtd.solve(
                geometry, wavelengths, port=out_port, input_mode=1, **kwargs
            )
            alpha_out = np.stack([
                out_smat.coeff(port=out_port, in_mode=1, out_mode=n)
                for n in range(1, N+1)
            ], axis=1)  # [M, N]

            # Use EME S[0,1] as the single-mode propagation approximation.
            # Full multi-mode hybrid: α_mm[m,n] per mode would come from EME per-mode output.
            alpha_mm_approx = eme_smat.coeff(port=out_port).reshape(-1, 1)  # [M, 1] broadcast

            # Step 4: combine
            alpha_full = np.sum(alpha_in * alpha_mm_approx * alpha_out, axis=1)  # [M]
            result.set(alpha_full, port=out_port, in_mode=input_mode, out_mode=1)

        return result
