# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Georgia Tech Research Corporation
"""Solver protocol: common interface for EME, FDTD, and Hybrid solvers."""
from __future__ import annotations

from typing import Protocol, runtime_checkable

import numpy as np

from ..geometry.base import DeviceGeometry
from ..results.smatrix import SMatrix


@runtime_checkable
class Solver(Protocol):
    """All BEEF solvers implement this interface."""

    def solve(self, geometry: DeviceGeometry, wavelengths: np.ndarray, **kwargs) -> SMatrix:
        """Run simulation and return S-parameter results.

        Args:
            geometry: device geometry object (MMI, DirectionalCoupler, CompoundGeometry, etc.).
            wavelengths: [M] array of wavelengths (µm).

        Returns:
            SMatrix containing complex amplitudes indexed by (port, in_mode, out_mode).
        """
        ...
