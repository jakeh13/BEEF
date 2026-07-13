"""HybridSolver: generic N-section cascade combining any mix of solvers.

Any solver type — FDTDSolver, EMESolver, or a custom class — can occupy any
position in the chain, in any order, as long as it implements:

    solve_section(geometry, wavelengths, x0, x1,
                  n_modes_in, n_modes_out,
                  stub_left=0.0, stub_right=0.0,
                  **geometry_kwargs) -> complex ndarray [M, n_modes_out, n_modes_in]

Each section's matrix is the modal transfer function between the section's
left and right interfaces. The full-device response is the matrix-chain
product of the per-section matrices, contracted over the mode index shared
at each internal interface:

    alpha_full[m] = M_K[m] @ M_{K-1}[m] @ ... @ M_1[m]

This is the direct generalization of the slide-4 formula
    alpha_full[m] = sum_n alpha_in[m,n] * alpha_mm[m,n] * alpha_out[m,n]
to an arbitrary number of sections (>= 2) and arbitrary (possibly
non-diagonal) per-section transfer matrices — the sum over n above is exactly
what the matrix product contracts over.

Nothing in this module knows about any particular device type; the geometry
and its section boundaries are supplied by the caller.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from ..results.smatrix import SMatrix


class HybridSolver:
    """Chains an arbitrary sequence of solvers across contiguous device sections.

    Args:
        solvers: ordered list of solver instances (length >= 2), left to
            right. Any mix of solver types/order is supported as long as
            each implements solve_section() (see module docstring).
        boundaries: list of len(solvers)+1 x-positions (device frame, µm)
            partitioning the device. boundaries[0] / boundaries[-1] are the
            true physical port locations; the len(solvers)-1 interior points
            are the cuts between consecutive solvers.
        interface_modes: list of len(solvers)-1 ints — number of modes
            carried at each internal interface between consecutive solvers.
    """

    def __init__(
        self,
        solvers: List[Any],
        boundaries: List[float],
        interface_modes: List[int],
    ):
        if len(solvers) < 2:
            raise ValueError("HybridSolver needs at least 2 solvers.")
        if len(boundaries) != len(solvers) + 1:
            raise ValueError(
                f"boundaries must have len(solvers)+1 = {len(solvers) + 1} "
                f"entries, got {len(boundaries)}."
            )
        if len(interface_modes) != len(solvers) - 1:
            raise ValueError(
                f"interface_modes must have len(solvers)-1 = {len(solvers) - 1} "
                f"entries, got {len(interface_modes)}."
            )
        self.solvers = solvers
        self.boundaries = boundaries
        self.interface_modes = interface_modes

    def solve(
        self,
        geometry,
        wavelengths: np.ndarray,
        in_modes: int = 1,
        out_modes: int = 1,
        port: str = "out",
        geometry_kwargs: Optional[Any] = None,
        cache: Optional[Dict[Any, np.ndarray]] = None,
    ) -> SMatrix:
        """Run every section solver and cascade the resulting transfer matrices.

        Args:
            geometry: device geometry implementing section_prisms/section_gds/
                width_at/center_at (see feces.geometry.base.DeviceGeometry).
            wavelengths: [M] wavelength array (µm).
            in_modes: mode count at the true left (input) port. Usually 1.
            out_modes: mode count at the true right (output) port. Usually 1.
            port: port label to store the result under in the returned SMatrix.
            geometry_kwargs: forwarded to solve_section() calls. Either a
                single dict applied to every section, or a list of
                len(solvers) dicts giving different kwargs per section (e.g.
                port_idx only on the section that actually branches — see
                cache below for why this distinction matters).
            cache: optional dict shared across multiple solve() calls (e.g.
                once per output port of a multi-output device). Sections are
                memoized by (solver identity, x0, x1, n_in, n_out, kwargs) —
                so a section whose own kwargs don't change between calls
                (typically every section upstream of the one that branches)
                is only ever actually simulated once.

        Returns:
            SMatrix with coeff(port=port, in_mode=i, out_mode=o) = alpha_full[m].
        """
        wavelengths = np.asarray(wavelengths)
        K = len(self.solvers)

        if geometry_kwargs is None:
            kwargs_per_section = [{}] * K
        elif isinstance(geometry_kwargs, dict):
            kwargs_per_section = [geometry_kwargs] * K
        else:
            kwargs_per_section = list(geometry_kwargs)
            if len(kwargs_per_section) != K:
                raise ValueError(
                    f"geometry_kwargs list must have len(solvers) = {K} "
                    f"entries, got {len(kwargs_per_section)}."
                )

        n_in_list = [in_modes] + list(self.interface_modes)
        n_out_list = list(self.interface_modes) + [out_modes]

        mats = []
        for k, solver in enumerate(self.solvers):
            x0, x1 = self.boundaries[k], self.boundaries[k + 1]
            stub_left = solver.waveguide_length if (k == 0 and hasattr(solver, "waveguide_length")) else 0.0
            stub_right = solver.waveguide_length if (k == K - 1 and hasattr(solver, "waveguide_length")) else 0.0
            kwargs = kwargs_per_section[k]

            cache_key = None
            if cache is not None:
                cache_key = (id(solver), x0, x1, n_in_list[k], n_out_list[k],
                             tuple(sorted(kwargs.items())))
                if cache_key in cache:
                    print(f"[Hybrid] Section {k + 1}/{K}: {type(solver).__name__} (cached)")
                    mats.append(cache[cache_key])
                    continue

            print(
                f"[Hybrid] Section {k + 1}/{K}: {type(solver).__name__} "
                f"x=[{x0:.3f}, {x1:.3f}]  modes {n_in_list[k]}->{n_out_list[k]}"
            )
            mat = solver.solve_section(
                geometry, wavelengths, x0, x1,
                n_modes_in=n_in_list[k], n_modes_out=n_out_list[k],
                stub_left=stub_left, stub_right=stub_right,
                **kwargs,
            )
            mat = np.asarray(mat, dtype=complex)
            if cache is not None:
                cache[cache_key] = mat
            mats.append(mat)

        # Cascade: full[m] = mats[-1][m] @ ... @ mats[0][m]
        full = mats[0]
        for mat in mats[1:]:
            full = np.einsum("mij,mjk->mik", mat, full)  # [M, out, in]

        result = SMatrix(wavelengths=wavelengths)
        for o in range(out_modes):
            for i in range(in_modes):
                result.set(full[:, o, i], port=port, in_mode=i + 1, out_mode=o + 1)
        return result
