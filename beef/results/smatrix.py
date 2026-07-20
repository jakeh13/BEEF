"""SMatrix: unified S-parameter result type.

Replaces the inconsistent dict-keyed results (different tuple arity per sim type)
and ad-hoc pickle dumps across master and sweep scripts.
"""
from __future__ import annotations

import pickle
import warnings
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np


@dataclass
class SMatrix:
    """S-parameter result indexed by (port, in_mode, out_mode).

    All coefficients are complex; power = |coeff|^2.

    Attributes:
        wavelengths: [M] wavelength array (µm).
        _data: dict[(port, in_mode, out_mode)] -> np.ndarray[M] complex amplitudes.
    """

    wavelengths: np.ndarray
    _data: Dict[Tuple[str, int, int], np.ndarray] = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Building
    # ------------------------------------------------------------------

    def set(self, coeff: np.ndarray, port: str = "out",
            in_mode: int = 1, out_mode: int = 1) -> None:
        """Store complex mode overlap coefficient for a (port, in_mode, out_mode) key."""
        self._data[(port, in_mode, out_mode)] = np.asarray(coeff, dtype=complex)

    def set_power(self, power: np.ndarray, port: str = "out",
                  in_mode: int = 1, out_mode: int = 1) -> None:
        """Store real power (coeff = sqrt(power), phase unknown → stored as real)."""
        self._data[(port, in_mode, out_mode)] = np.sqrt(np.asarray(power, dtype=float)).astype(complex)

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def coeff(self, port: str = "out", in_mode: int = 1, out_mode: int = 1) -> np.ndarray:
        """Complex mode overlap coefficient array [M]."""
        return self._data[(port, in_mode, out_mode)]

    def transmission(self, port: str = "out", in_mode: int = 1, out_mode: int = 1) -> np.ndarray:
        """Power transmission [M] = |coeff|^2."""
        return np.abs(self._data[(port, in_mode, out_mode)]) ** 2

    def insertion_loss_db(self, port: str = "out", in_mode: int = 1, out_mode: int = 1) -> np.ndarray:
        """Insertion loss in dB [M] (negative = loss, e.g. −3 dB for 50/50 split)."""
        T = self.transmission(port, in_mode, out_mode)
        T = np.clip(T, 1e-20, None)
        return 10 * np.log10(T)

    def splitting_ratio(self, port_a: str, port_b: str,
                        in_mode: int = 1, out_mode: int = 1) -> np.ndarray:
        """Power splitting ratio: T_a / (T_a + T_b) [M]."""
        Ta = self.transmission(port_a, in_mode, out_mode)
        Tb = self.transmission(port_b, in_mode, out_mode)
        total = Ta + Tb
        total = np.where(total == 0, 1e-20, total)
        return Ta / total

    def ports(self) -> list:
        """All port labels present in this SMatrix."""
        return sorted({k[0] for k in self._data})

    # ------------------------------------------------------------------
    # Sanity checks
    # ------------------------------------------------------------------

    def total_power(self, in_mode: int = 1) -> np.ndarray:
        """Total output power [M] for a given input mode, summed over every
        (port, out_mode) entry present with that in_mode — i.e. all the
        power this input mode couples to anywhere in the device."""
        keys = [k for k in self._data if k[1] == in_mode]
        if not keys:
            raise KeyError(f"no SMatrix entries with in_mode={in_mode}")
        total = np.zeros(len(self.wavelengths))
        for k in keys:
            total += np.abs(self._data[k]) ** 2
        return total

    def check_passive(self, tol: float = 1e-3, warn: bool = True) -> bool:
        """Sanity check: for every input mode present, total output power
        (summed over every port/out_mode via total_power()) must not exceed
        1+tol at any wavelength — a passive device can only lose power
        (to radiation, reflection, or a truncated mode basis), never gain
        it. A violation almost always means a bug (mode-count mismatch, a
        double-counted output branch, wrong normalization) rather than real
        gain, and this check catches it immediately instead of requiring
        someone to notice an anomalous plot first.

        Args:
            tol: allowed slack above 1.0 before a wavelength/mode counts as
                a violation.
            warn: emit a warnings.warn() identifying the worst violation
                (input mode, wavelength, total power) when the check fails.

        Returns:
            True if every input mode passes at every wavelength.
        """
        in_modes = sorted({k[1] for k in self._data})
        ok = True
        for m in in_modes:
            total = self.total_power(in_mode=m)
            worst = np.argmax(total)
            if total[worst] > 1.0 + tol:
                ok = False
                if warn:
                    warnings.warn(
                        f"SMatrix.check_passive: in_mode={m} totals "
                        f"{total[worst]:.4f} (> 1+{tol:g}) at wavelength "
                        f"{self.wavelengths[worst]:.4f} — check for a mode-"
                        "count mismatch, double-counted output branch, or "
                        "normalization bug.",
                        stacklevel=2,
                    )
        return ok

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def save(self, path: str) -> None:
        with open(path, "wb") as f:
            pickle.dump({"wavelengths": self.wavelengths, "data": self._data}, f)

    @classmethod
    def load(cls, path: str) -> "SMatrix":
        with open(path, "rb") as f:
            d = pickle.load(f)
        obj = cls(wavelengths=d["wavelengths"])
        obj._data = d["data"]
        return obj

    # ------------------------------------------------------------------
    # Interop with legacy sweep dicts
    # ------------------------------------------------------------------

    @classmethod
    def from_legacy_powers(cls, wavelengths: np.ndarray,
                            powers_ports: list,
                            port_names: list) -> "SMatrix":
        """Build SMatrix from the legacy (powers_ports[port_idx][wl_idx]) format."""
        obj = cls(wavelengths=np.asarray(wavelengths))
        for port_idx, port_name in enumerate(port_names):
            obj.set_power(np.asarray(powers_ports[port_idx]), port=port_name)
        return obj

    def __repr__(self) -> str:
        ports = self.ports()
        return f"SMatrix(wavelengths={len(self.wavelengths)} pts, ports={ports})"
