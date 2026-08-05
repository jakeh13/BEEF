"""SMatrix: unified S-parameter result type.

Replaces the inconsistent dict-keyed results (different tuple arity per sim type)
and ad-hoc pickle dumps across master and sweep scripts.
"""
from __future__ import annotations

import pickle
import warnings
from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence, Tuple

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

    def keys(self) -> list:
        """All ``(port, in_mode, out_mode)`` keys present in this SMatrix,
        each retrievable via ``coeff()``/``transmission()``. Unlike ``ports()``
        (port labels only), this is every full key a caller comparing two
        SMatrix instances entry-by-entry needs — see
        ``beef.cut_position_invariance.cut_position_invariance()``."""
        return list(self._data.keys())

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
    # Coherent-receiver metrics.
    #
    # Added for the 4x4 MMI 90-degree hybrid (validation_markdown/
    # mmi_4x4_90_hybrid_validation.md, finding A8): every existing helper
    # above reads only |coeff|^2 (power). imbalance_db/excess_loss_db/
    # cmrr_db are still power-only and generalize cleanly to a single
    # SMatrix (one input excitation, several output ports). Phase error
    # (paper Eq. 15) is different in kind: it needs the phase DIFFERENCE
    # between two SEPARATELY-excited inputs, and every SMatrix in this repo
    # so far holds exactly one input excitation's results (e.g. mmi_2x2_
    # hybrid.py's single `smat` for its one driven arm) -- so that one takes
    # a second SMatrix explicitly rather than pretending to be a single-
    # instance method. The paper's own reference-port subtraction and
    # hardcoded target angles (phi~_5=-45 deg, etc.) are device-specific, not
    # generalized here -- cross_excitation_phase_diff_deg is the reusable
    # primitive; building Eq. 15 from it is the example script's job.
    # ------------------------------------------------------------------

    def imbalance_db(self, ports: Sequence[str], in_mode: int = 1, out_mode: int = 1) -> np.ndarray:
        """Power imbalance across several output ports of ONE input excitation
        (paper Eq. 10, generalized from a fixed pair to any number of ports):
        ``-10*log10(min_j T_j / max_j T_j)`` [M]. 0 dB = perfectly balanced;
        larger is worse."""
        if len(ports) < 2:
            raise ValueError(f"imbalance_db needs >= 2 ports, got {list(ports)!r}")
        Ts = np.stack([self.transmission(p, in_mode, out_mode) for p in ports])
        tmax = np.clip(Ts.max(axis=0), 1e-20, None)
        tmin = np.clip(Ts.min(axis=0), 1e-20, None)
        return -10 * np.log10(tmin / tmax)

    def excess_loss_db(self, ports: Sequence[str], in_mode: int = 1, out_mode: int = 1,
                        p_in: float = 1.0) -> np.ndarray:
        """Excess loss summed over several output ports of ONE input
        excitation (paper Eq. 11): ``-10*log10(sum_j T_j / p_in)`` [M].

        Positive = loss -- the paper's own sign convention, and the OPPOSITE
        of this class's own insertion_loss_db() (negative = loss). Kept as
        written so a value can be compared directly against a paper table
        without a sign flip; do not mix the two conventions in one report.
        """
        total = np.clip(sum(self.transmission(p, in_mode, out_mode) for p in ports), 1e-20, None)
        return -10 * np.log10(total / p_in)

    def cmrr_db(self, port_a: str, port_b: str, in_mode: int = 1, out_mode: int = 1) -> np.ndarray:
        """Common-mode rejection ratio between two output ports of ONE input
        excitation (paper Eqs. 13-14 -- CMRR_I compares a 4x4 90-degree
        hybrid's ports 5/8, CMRR_Q compares 6/7, but the formula itself is
        generic to any balanced two-port pair):
        ``-20*log10(|T_a - T_b| / (T_a + T_b))`` [M]. Larger (dBe) is better
        (more common-mode rejection)."""
        Ta = self.transmission(port_a, in_mode, out_mode)
        Tb = self.transmission(port_b, in_mode, out_mode)
        total = np.clip(Ta + Tb, 1e-20, None)
        diff = np.clip(np.abs(Ta - Tb), 1e-20, None)
        return -20 * np.log10(diff / total)

    def cross_excitation_phase_diff_deg(self, other: "SMatrix", port: str,
                                         in_mode: int = 1, out_mode: int = 1) -> np.ndarray:
        """Phase difference [M, degrees, wrapped to (-180, 180]] between two
        SMatrix instances built from DIFFERENT input excitations of the same
        multi-input device, at the same (port, in_mode, out_mode):
        ``angle(other.coeff) - angle(self.coeff)``.

        This is the ingredient paper Eq. 15's ``phi_j(lambda) = angle(S_j3) -
        angle(S_j1)`` needs -- pass ``self`` as the port-1-excited SMatrix and
        ``other`` as the port-3-excited one to reproduce phi_j exactly. The
        paper's own reference-port subtraction (``- (phi~_j - phi~_5)``) and
        target angles are device-specific and are NOT applied here; compute
        them in the caller from this primitive plus the device's own target
        angles.

        Output is wrapped to (-180, 180]. A target stated in [0, 360)
        (e.g. this paper's own phi~_6 = 225 deg) will not compare equal by
        direct subtraction against this range -- compare modulo 360, not by
        assuming a shared quadrant.
        """
        a = self.coeff(port, in_mode, out_mode)
        b = other.coeff(port, in_mode, out_mode)
        diff = np.degrees(np.angle(b) - np.angle(a))
        return ((diff + 180.0) % 360.0) - 180.0

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
