# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Georgia Tech Research Corporation
"""Conversions between beef.network's SMatrixBlock and the rest of BEEF.

The only file in beef.network that imports beef.results.smatrix.SMatrix —
block.py/ops.py/circuit.py depend on numpy only, so nothing outside this
file couples the network layer to SMatrix's schema.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from ..results.smatrix import SMatrix
from .block import SMatrixBlock


def from_transfer_matrix(
    matrix: np.ndarray,
    wavelengths: np.ndarray,
    port_in: str = "left",
    port_out: str = "right",
) -> SMatrixBlock:
    """Wrap a HybridSolver-style one-directional transfer matrix as a
    degenerate 2-port SMatrixBlock.

    `matrix` is [M, n_out, n_in] with alpha_out = matrix @ alpha_in — the
    exact contract of beef.solvers.hybrid's solve_section()/
    solve_section_length_sweep() (for one length slice). Self-terms
    (reflection at either port, or reverse transmission) are zero-filled,
    matching what those solvers actually model: pure forward transmission,
    no reflection.
    """
    matrix = np.asarray(matrix, dtype=complex)
    wavelengths = np.asarray(wavelengths, dtype=float)
    M = len(wavelengths)
    if matrix.shape[0] != M:
        raise ValueError(
            f"from_transfer_matrix: matrix has {matrix.shape[0]} wavelength "
            f"entries, wavelengths has {M}."
        )
    n_out, n_in = matrix.shape[1], matrix.shape[2]
    N = n_in + n_out
    data = np.zeros((M, N, N), dtype=complex)
    data[:, n_in:, :n_in] = matrix  # rows=port_out, cols=port_in
    return SMatrixBlock(
        wavelengths=wavelengths,
        ports=(port_in, port_out),
        modes_per_port={port_in: n_in, port_out: n_out},
        data=data,
    )


def from_smatrix(
    smat: SMatrix,
    input_port: str,
    output_ports: Optional[Sequence[str]] = None,
) -> SMatrixBlock:
    """Wrap one SMatrix (the response to a single excitation config) as a block.

    `input_port` names whichever port was actually excited to produce
    `smat` — SMatrix has no port-direction metadata of its own (a
    (port, in_mode, out_mode) entry never records which port in_mode
    belongs to), so the caller must supply it. Self-terms (input_port's
    own reflection, and cross-terms between two different output ports)
    are zero-filled — an honest limitation of what SMatrix records, not an
    approximation this function introduces.

    Raises:
        ValueError: no entries found for the requested output_ports, or
            the in_mode/out_mode values present are not a contiguous
            1..n range (SMatrixBlock requires a dense mode-index layout).
    """
    if output_ports is None:
        output_ports = smat.ports()
    output_ports = list(output_ports)

    keys = [k for k in smat._data if k[0] in output_ports]
    if not keys:
        raise ValueError(
            f"from_smatrix: no entries found for output_ports={output_ports}."
        )
    in_modes = sorted({k[1] for k in keys})
    n_in = max(in_modes)
    if in_modes != list(range(1, n_in + 1)):
        raise ValueError(
            f"from_smatrix: in_mode values {in_modes} for "
            f"output_ports={output_ports} are not contiguous 1..{n_in}."
        )

    modes_per_out = {}
    for port in output_ports:
        out_modes = sorted({k[2] for k in keys if k[0] == port})
        if not out_modes:
            raise ValueError(f"from_smatrix: no entries found for output port {port!r}.")
        n_out = max(out_modes)
        if out_modes != list(range(1, n_out + 1)):
            raise ValueError(
                f"from_smatrix: out_mode values {out_modes} for port "
                f"{port!r} are not contiguous 1..{n_out}."
            )
        modes_per_out[port] = n_out

    ports = (input_port,) + tuple(output_ports)
    modes_per_port = {input_port: n_in, **modes_per_out}
    M = len(smat.wavelengths)
    N = n_in + sum(modes_per_out.values())
    block = SMatrixBlock(
        wavelengths=smat.wavelengths, ports=ports, modes_per_port=modes_per_port,
        data=np.zeros((M, N, N), dtype=complex),
    )
    for out_port, n_out in modes_per_out.items():
        for i in range(1, n_in + 1):
            for o in range(1, n_out + 1):
                try:
                    coeff = smat.coeff(port=out_port, in_mode=i, out_mode=o)
                except KeyError:
                    continue
                block.data[:, block.flat_index(out_port, o), block.flat_index(input_port, i)] = coeff
    return block


def to_smatrix(block: SMatrixBlock, input_port: str) -> SMatrix:
    """Flatten a (typically post-composition) block back to SMatrix's
    (port, in_mode, out_mode) schema.

    Includes input_port itself among the output ports if
    block.submatrix(input_port, input_port) is nonzero — i.e. this can
    express a reflection coefficient at the excited port, which no
    existing example's SMatrix has ever populated (a strict superset of
    today's usage, not a behavior change to anything that reads SMatrix).
    """
    if input_port not in block.ports:
        raise ValueError(
            f"to_smatrix: input_port={input_port!r} not in block (have {block.ports})."
        )
    smat = SMatrix(wavelengths=block.wavelengths)
    n_in = block.modes_per_port[input_port]
    col = block.port_slice(input_port).start
    for port in block.ports:
        n_out = block.modes_per_port[port]
        row0 = block.port_slice(port).start
        for i in range(1, n_in + 1):
            for o in range(1, n_out + 1):
                coeff = block.data[:, row0 + (o - 1), col + (i - 1)]
                smat.set(coeff, port=port, in_mode=i, out_mode=o)
    return smat
