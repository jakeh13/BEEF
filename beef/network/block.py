"""SMatrixBlock: dense multi-port S-parameter block for network composition.

Data model for beef.network's Redheffer-star-product composition layer
(see ops.py). Distinct from beef.results.smatrix.SMatrix, which stores a
flat, sparse dict of (port, in_mode, out_mode) -> amplitude and has no
notion of a port's own reflection or of two different output ports'
cross-coupling; SMatrixBlock is a dense [M, N, N] matrix over every port's
modes so that connect()/close_loop() have something to do real linear
algebra on. convert.py bridges the two.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Tuple

import numpy as np


@dataclass
class SMatrixBlock:
    """Dense S-parameter block over named, multi-mode ports.

    data[m] @ a[m] = b[m]: for wavelength index m, `data[m]` maps a vector
    of input amplitudes (one entry per mode of every port, concatenated in
    `ports` order) to a vector of output amplitudes in the same layout.
    Every port is bidirectional — its in-mode count and out-mode count are
    both `modes_per_port[port]` — matching standard S-parameter convention
    (as opposed to beef.solvers.hybrid's one-directional transfer matrices,
    which have independent in/out mode counts; see convert.from_transfer_matrix
    for how those become a degenerate 2-port block here).

    Attributes:
        wavelengths: [M] real wavelength array (µm).
        ports: ordered port names. Order defines the flat-index layout used
            by `data` — port `ports[0]` occupies indices [0, modes_per_port[ports[0]]),
            `ports[1]` the next block, etc.
        modes_per_port: mode count for each name in `ports`.
        data: [M, N, N] complex, N = sum(modes_per_port.values()).
    """

    wavelengths: np.ndarray
    ports: Tuple[str, ...]
    modes_per_port: Dict[str, int]
    data: np.ndarray
    _port_slices: Dict[str, slice] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.wavelengths = np.asarray(self.wavelengths, dtype=float)
        if self.wavelengths.ndim != 1:
            raise ValueError(
                f"wavelengths must be 1-D, got shape {self.wavelengths.shape}."
            )
        self.ports = tuple(self.ports)
        if len(set(self.ports)) != len(self.ports):
            raise ValueError(f"ports must be unique, got {self.ports}.")
        missing = [p for p in self.ports if p not in self.modes_per_port]
        if missing:
            raise ValueError(
                f"modes_per_port is missing entries for ports {missing}."
            )
        extra = [p for p in self.modes_per_port if p not in self.ports]
        if extra:
            raise ValueError(
                f"modes_per_port has entries for ports not in `ports`: {extra}."
            )
        bad_counts = {p: n for p, n in self.modes_per_port.items() if n < 1}
        if bad_counts:
            raise ValueError(
                f"modes_per_port values must be >= 1, got {bad_counts}."
            )

        slices: Dict[str, slice] = {}
        offset = 0
        for p in self.ports:
            n = self.modes_per_port[p]
            slices[p] = slice(offset, offset + n)
            offset += n
        self._port_slices = slices
        N = offset

        self.data = np.asarray(self.data, dtype=complex)
        M = len(self.wavelengths)
        if self.data.shape != (M, N, N):
            raise ValueError(
                f"data must have shape (M, N, N) = ({M}, {N}, {N}) given "
                f"{len(self.wavelengths)} wavelengths and total mode count "
                f"{N} (ports={self.ports}, modes_per_port={self.modes_per_port}); "
                f"got {self.data.shape}."
            )

    def port_slice(self, port: str) -> slice:
        """Flat-index slice occupied by `port` on either axis of `data`."""
        if port not in self._port_slices:
            raise KeyError(f"unknown port {port!r}; have {self.ports}.")
        return self._port_slices[port]

    def n_modes(self, port: str) -> int:
        """Mode count carried by `port`."""
        return self.modes_per_port[port]

    def submatrix(self, port_out: str, port_in: str) -> np.ndarray:
        """[M, n_out, n_in] block relating `port_in`'s input to `port_out`'s output."""
        return self.data[:, self.port_slice(port_out), self.port_slice(port_in)]

    def flat_index(self, port: str, mode: int) -> int:
        """Flat row/column index for `port`'s 1-indexed `mode`."""
        sl = self.port_slice(port)
        if not (1 <= mode <= sl.stop - sl.start):
            raise ValueError(
                f"mode={mode} out of range for port {port!r} "
                f"(has {sl.stop - sl.start} modes)."
            )
        return sl.start + (mode - 1)

    def __repr__(self) -> str:
        return (
            f"SMatrixBlock(wavelengths={len(self.wavelengths)} pts, "
            f"ports={ {p: self.modes_per_port[p] for p in self.ports} })"
        )


def rename_ports(block: SMatrixBlock, mapping: Dict[str, str]) -> SMatrixBlock:
    """Return a new SMatrixBlock with `mapping`'s named ports relabeled.

    Ports not present in `mapping` keep their existing name. Does not copy
    or touch `data` — only `ports`/`modes_per_port` are rebuilt, so this is
    cheap regardless of block size.
    """
    new_ports = tuple(mapping.get(p, p) for p in block.ports)
    if len(set(new_ports)) != len(new_ports):
        collisions = {p for p in new_ports if new_ports.count(p) > 1}
        raise ValueError(
            f"rename_ports: mapping {mapping} produces duplicate port "
            f"name(s) {collisions} from {block.ports}."
        )
    new_modes = {mapping.get(p, p): n for p, n in block.modes_per_port.items()}
    return SMatrixBlock(
        wavelengths=block.wavelengths,
        ports=new_ports,
        modes_per_port=new_modes,
        data=block.data,
    )


def check_passive(block: SMatrixBlock, tol: float = 1e-3, warn: bool = True) -> bool:
    """Sanity check: for every input mode, total output power (summed over
    every port/out-mode, including the exciting port's own reflection) must
    not exceed 1+tol at any wavelength.

    Same contract as beef.results.smatrix.SMatrix.check_passive(), operating
    directly on the dense block instead of a sparse dict — not auto-invoked
    by connect()/close_loop(), callers invoke it explicitly, matching how
    every example calls smat.check_passive() itself today.
    """
    import warnings

    total = np.sum(np.abs(block.data) ** 2, axis=1)  # [M, N] summed over output rows
    worst_idx = np.unravel_index(np.argmax(total), total.shape)
    worst = total[worst_idx]
    ok = bool(worst <= 1.0 + tol)
    if not ok and warn:
        m_idx, col = worst_idx
        port_in = next(p for p in block.ports if block.port_slice(p).start <= col < block.port_slice(p).stop)
        mode_in = col - block.port_slice(port_in).start + 1
        warnings.warn(
            f"network.check_passive: input (port={port_in!r}, mode={mode_in}) "
            f"totals {worst:.4f} (> 1+{tol:g}) at wavelength "
            f"{block.wavelengths[m_idx]:.4f} — check for a mode-count "
            "mismatch or a composition bug.",
            stacklevel=2,
        )
    return ok
