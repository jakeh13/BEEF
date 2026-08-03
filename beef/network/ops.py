"""connect() / close_loop(): Redheffer star product and self-loop elimination.

Free functions rather than methods — connect() in particular takes two peer
blocks with no natural "self" (mirrors resolved_material_config() operating
on a MaterialStack in beef/materials/stack.py, the existing precedent for a
derived-object transform living outside its input dataclass).

Both operate on the S-parameter convention: for a block with data[m] (see
SMatrixBlock), b[m] = data[m] @ a[m] where a/b are input/output amplitude
vectors over every port's modes.

connect(): standard Redheffer star product joining one port of block_a to
one port of block_b (arbitrary, possibly different, blocks).

close_loop(): eliminates a block's OWN port against itself — the general
form of the ring resonator's proven closed-form all-pass formula
(examples/ring_resonator/python/ring_resonator_hybrid.py:762-774):

    E_thru = t + k1*k2*P / (1 - t_ring*P)

which is exactly close_loop() applied to a 2-port [bus, ring] block with
S_bus,bus=t, S_bus,ring=k2, S_ring,bus=k1, S_ring,ring=t_ring, loop=P — see
tests/test_network_ring_crosscheck.py.

Only the port_out == port_in case ("Branch A": one port eliminated against
itself) is implemented in this pass. The port_out != port_in case ("Branch
B": two DIFFERENT ports of the same block closed against each other, e.g.
when Circuit.reduce() discovers a cycle joining two distinctly-named ports)
is a materially different elimination — Branch A removes one free variable
entirely; Branch B starts from two independent free variables (a_p, a_q)
and removes only one of them (a_q, forced by the loop), leaving a_p and its
own port still live. Naively taking the p->q limit of a Branch-B-shaped
formula does not reduce to Branch A's, so this isn't a single formula with
a degenerate case — it's two different formulas. Branch B raises
NotImplementedError with this same explanation.

Branch B DOES now have a required caller, identified 2026-08-01: the ring
resonator's physical 4-port form. Its loop runs ring_R -> around the outer
ring -> ring_L, i.e. port_out != port_in, so the composable/add-drop/cascaded
route is blocked on this function. See validation_markdown/
ring_resonator_validation.md's "Wiring to-do" item 1. That doc's decision is
to ship the 2-port aggregate form (Branch A) first as a numerical
cross-check, then implement Branch B -- so this is scheduled work with a
named dependent, no longer an open-ended deferral.
"""
from __future__ import annotations

import warnings
from typing import List, Optional, Sequence, Union

import numpy as np

from .block import SMatrixBlock

Number = Union[int, float, complex]


def _check_same_wavelengths(block_a: SMatrixBlock, block_b: SMatrixBlock, fname: str) -> None:
    if not np.array_equal(block_a.wavelengths, block_b.wavelengths):
        raise ValueError(
            f"{fname}: block_a and block_b have different wavelength arrays "
            f"({len(block_a.wavelengths)} vs {len(block_b.wavelengths)} points, "
            "or matching lengths with different values) — beef.network "
            "performs no wavelength interpolation; resample one side to "
            "match the other before composing (see e.g. how "
            "ring_resonator_hybrid.py's interp_complex() does this)."
        )


def _group_indices(block: SMatrixBlock, group_ports: Sequence[str]) -> np.ndarray:
    """Flat row/column indices spanning every port in `group_ports`, in order."""
    if not group_ports:
        return np.array([], dtype=int)
    return np.concatenate(
        [np.arange(block.port_slice(p).start, block.port_slice(p).stop) for p in group_ports]
    )


def _submatrix_by_index(data: np.ndarray, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    """[M, len(rows), len(cols)] selection, robust to either index array being empty."""
    M = data.shape[0]
    if rows.size == 0 or cols.size == 0:
        return np.zeros((M, rows.size, cols.size), dtype=complex)
    return data[:, rows[:, None], cols[None, :]]


def _check_conditioning(
    matrix: np.ndarray, wavelengths: np.ndarray, context: str,
    cond_warn: float, cond_error: float,
) -> None:
    """Flag a near-singular star-product/loop-closure system.

    Deliberately NOT np.linalg.cond() (the standard relative condition
    number ||M|| * ||M^-1||): for a 1x1 matrix — exactly the case for a
    single-mode loop port, i.e. the ring resonator's real use case —
    cond() is identically 1.0 no matter how close M is to singular
    (verified: cond([[1e-9]]) == 1.0), because the *relative* sensitivity
    of solving Mx=b to a relative perturbation in M is genuinely O(1) even
    as |M| -> 0. What actually matters here is the *absolute* nearness to
    singular — how large np.linalg.solve's answer gets for a unit-scale
    RHS — which is the reciprocal smallest singular value, 1/sigma_min(M).
    For a 1x1 M=[[m]] this is exactly 1/|m|, so it correctly blows up as a
    ring approaches critical coupling (t_ring*P -> 1, m = 1 - t_ring*P -> 0).
    """
    sigma_min = np.linalg.svd(matrix, compute_uv=False)[..., -1]
    with np.errstate(divide="ignore"):
        inv_sigma_min = 1.0 / sigma_min  # sigma_min==0 (exact singular) -> inf, caught below
    worst = int(np.argmax(inv_sigma_min))
    if inv_sigma_min[worst] > cond_error:
        raise ValueError(
            f"{context}: near-singular system at wavelength "
            f"{wavelengths[worst]:.6g} (1/sigma_min = {inv_sigma_min[worst]:.3e} "
            f"> cond_error={cond_error:.1e}) — this is usually a real "
            "resonance/pole (e.g. a closed loop at or past critical "
            "coupling), not numerical noise. Check the physical "
            "configuration before raising cond_error to silence this."
        )
    if inv_sigma_min[worst] > cond_warn:
        warnings.warn(
            f"{context}: poorly conditioned system at wavelength "
            f"{wavelengths[worst]:.6g} (1/sigma_min = {inv_sigma_min[worst]:.3e} "
            f"> cond_warn={cond_warn:.1e}) — results may be numerically "
            "unreliable near this wavelength; this often indicates "
            "proximity to a real resonance (e.g. near-critical coupling).",
            stacklevel=3,
        )


def connect(
    block_a: SMatrixBlock, port_a: str,
    block_b: SMatrixBlock, port_b: str,
    *,
    cond_warn: float = 1e6,
    cond_error: float = 1e12,
) -> SMatrixBlock:
    """Redheffer star product joining block_a's port_a to block_b's port_b.

    Every other port of both blocks survives in the result, named
    block_a's surviving ports (in block_a's order) followed by block_b's
    surviving ports (in block_b's order). block_a/block_b are not mutated.

    When both blocks come from convert.from_transfer_matrix() (pure
    one-directional transmission, zero self-terms), this reduces exactly
    to ordinary transfer-matrix multiplication — see
    tests/test_network.py::TestConnect for the cross-check against
    HybridSolver-style manual einsum chains.

    Raises:
        ValueError: mismatched wavelength arrays, unknown port name,
            mismatched mode counts between port_a/port_b, a port-name
            collision among the surviving ports, or a near/exactly-singular
            star-product system (see cond_warn/cond_error).
    """
    _check_same_wavelengths(block_a, block_b, "connect")
    if port_a not in block_a.ports:
        raise ValueError(f"connect: port_a={port_a!r} not in block_a (have {block_a.ports}).")
    if port_b not in block_b.ports:
        raise ValueError(f"connect: port_b={port_b!r} not in block_b (have {block_b.ports}).")

    n = block_a.modes_per_port[port_a]
    n_b = block_b.modes_per_port[port_b]
    if n != n_b:
        raise ValueError(
            f"connect: port_a={port_a!r} has {n} modes, port_b={port_b!r} "
            f"has {n_b} modes — mode counts must match to connect."
        )

    E_a = [p for p in block_a.ports if p != port_a]
    E_b = [p for p in block_b.ports if p != port_b]
    collide = set(E_a) & set(E_b)
    if collide:
        raise ValueError(
            f"connect: port name collision between block_a and block_b: "
            f"{sorted(collide)} — rename_ports() one side first."
        )

    wavelengths = block_a.wavelengths
    M = len(wavelengths)
    idx_a_E = _group_indices(block_a, E_a)
    idx_a_x = np.arange(block_a.port_slice(port_a).start, block_a.port_slice(port_a).stop)
    idx_b_E = _group_indices(block_b, E_b)
    idx_b_x = np.arange(block_b.port_slice(port_b).start, block_b.port_slice(port_b).stop)

    S_a_EE = _submatrix_by_index(block_a.data, idx_a_E, idx_a_E)
    S_a_Ex = _submatrix_by_index(block_a.data, idx_a_E, idx_a_x)
    S_a_xE = _submatrix_by_index(block_a.data, idx_a_x, idx_a_E)
    S_a_xx = _submatrix_by_index(block_a.data, idx_a_x, idx_a_x)

    S_b_EE = _submatrix_by_index(block_b.data, idx_b_E, idx_b_E)
    S_b_Ex = _submatrix_by_index(block_b.data, idx_b_E, idx_b_x)
    S_b_xE = _submatrix_by_index(block_b.data, idx_b_x, idx_b_E)
    S_b_xx = _submatrix_by_index(block_b.data, idx_b_x, idx_b_x)

    I_n = np.broadcast_to(np.eye(n, dtype=complex), (M, n, n))
    M1 = I_n - S_b_xx @ S_a_xx
    M2 = I_n - S_a_xx @ S_b_xx
    _check_conditioning(M1, wavelengths, "connect (block_a-facing system)", cond_warn, cond_error)
    _check_conditioning(M2, wavelengths, "connect (block_b-facing system)", cond_warn, cond_error)

    solveM1_Sbxx = np.linalg.solve(M1, S_b_xx)
    solveM1_SbxE = np.linalg.solve(M1, S_b_xE) if idx_b_E.size else np.zeros((M, n, 0), dtype=complex)
    solveM2_Saxx = np.linalg.solve(M2, S_a_xx)
    solveM2_SaxE = np.linalg.solve(M2, S_a_xE) if idx_a_E.size else np.zeros((M, n, 0), dtype=complex)

    S_EE_a_new = S_a_EE + S_a_Ex @ solveM1_Sbxx @ S_a_xE
    S_Ea_Eb_new = S_a_Ex @ solveM1_SbxE
    S_Eb_Ea_new = S_b_Ex @ solveM2_SaxE
    S_EE_b_new = S_b_EE + S_b_Ex @ solveM2_Saxx @ S_b_xE

    N_Ea, N_Eb = idx_a_E.size, idx_b_E.size
    new_data = np.zeros((M, N_Ea + N_Eb, N_Ea + N_Eb), dtype=complex)
    new_data[:, :N_Ea, :N_Ea] = S_EE_a_new
    new_data[:, :N_Ea, N_Ea:] = S_Ea_Eb_new
    new_data[:, N_Ea:, :N_Ea] = S_Eb_Ea_new
    new_data[:, N_Ea:, N_Ea:] = S_EE_b_new

    new_ports = tuple(E_a + E_b)
    new_modes = {p: block_a.modes_per_port[p] for p in E_a}
    new_modes.update({p: block_b.modes_per_port[p] for p in E_b})
    return SMatrixBlock(wavelengths=wavelengths, ports=new_ports, modes_per_port=new_modes, data=new_data)


def _as_loop_operator(loop, n: int, wavelengths: np.ndarray) -> np.ndarray:
    """Normalize `loop` to a [M, n, n] complex array.

    Accepts None (identity), a Python/NumPy scalar (scalar * I, broadcast
    over every wavelength), a [M] array (per-wavelength scalar * I — the
    ring script's actual `P`), or a full [M, n, n] array.
    """
    M = len(wavelengths)
    eye = np.eye(n, dtype=complex)
    if loop is None:
        return np.broadcast_to(eye, (M, n, n)).copy()
    arr = np.asarray(loop)
    if arr.ndim == 0:
        return complex(arr) * np.broadcast_to(eye, (M, n, n))
    if arr.ndim == 1:
        if arr.shape[0] != M:
            raise ValueError(
                f"close_loop: loop array of shape {arr.shape} does not "
                f"match {M} wavelengths."
            )
        return arr.astype(complex)[:, None, None] * eye[None, :, :]
    if arr.ndim == 3:
        if arr.shape != (M, n, n):
            raise ValueError(
                f"close_loop: loop array shape {arr.shape} != expected "
                f"({M}, {n}, {n})."
            )
        return arr.astype(complex)
    raise ValueError(
        f"close_loop: loop must be None, a scalar, a [M] array, or a "
        f"[M, n, n] array; got ndim={arr.ndim}."
    )


def close_loop(
    block: SMatrixBlock, port_out: str, port_in: str,
    loop: Optional[Union[Number, np.ndarray]] = None,
    *,
    cond_warn: float = 1e6,
    cond_error: float = 1e12,
) -> SMatrixBlock:
    """Eliminate block's own port(s) by feeding port_out's output through
    `loop` into port_in's input: a_{port_in} = loop @ b_{port_out}.

    Only port_out == port_in ("Branch A") is implemented in this pass —
    see the module docstring for why the port_out != port_in case
    ("Branch B") is a different formula, not a degenerate case of this one.

    For port_out == port_in == p, with the block's remaining ports E:
        b_p = (I - S_pp @ loop)^-1 @ S_pE @ a_E
        b_E = [S_EE + S_Ep @ loop @ (I - S_pp @ loop)^-1 @ S_pE] @ a_E
    Substituting a 2-port [bus, ring] block's S_bus,bus=t, S_bus,ring=k2,
    S_ring,bus=k1, S_ring,ring=t_ring and loop=P reproduces
    t + k1*k2*P/(1 - t_ring*P) exactly — the ring resonator's proven
    closed-form formula (see tests/test_network_ring_crosscheck.py).

    Raises:
        ValueError: unknown port, or a near/exactly-singular loop system
            (see cond_warn/cond_error) — a near-critical-coupling ring
            (t_ring*P -> 1) genuinely produces this condition; don't raise
            cond_error to silence it without checking whether the device
            is actually near critical coupling.
        NotImplementedError: port_out != port_in.
    """
    if port_out != port_in:
        raise NotImplementedError(
            "close_loop: port_out != port_in (closing two DIFFERENT ports "
            "of the same block against each other) is not implemented in "
            "this pass — see this module's docstring for the derivation "
            "and why it is a different formula from the port_out == "
            "port_in case, not a degenerate case of it. Only self-loop "
            "closure (port_out == port_in) is supported today."
        )
    p = port_out
    if p not in block.ports:
        raise ValueError(f"close_loop: port {p!r} not in block (have {block.ports}).")

    n = block.modes_per_port[p]
    wavelengths = block.wavelengths
    M = len(wavelengths)
    P = _as_loop_operator(loop, n, wavelengths)

    E = [q for q in block.ports if q != p]
    idx_E = _group_indices(block, E)
    idx_p = np.arange(block.port_slice(p).start, block.port_slice(p).stop)

    S_EE = _submatrix_by_index(block.data, idx_E, idx_E)
    S_Ep = _submatrix_by_index(block.data, idx_E, idx_p)
    S_pE = _submatrix_by_index(block.data, idx_p, idx_E)
    S_pp = _submatrix_by_index(block.data, idx_p, idx_p)

    I_n = np.broadcast_to(np.eye(n, dtype=complex), (M, n, n))
    Minv_sys = I_n - S_pp @ P
    _check_conditioning(Minv_sys, wavelengths, "close_loop", cond_warn, cond_error)

    solved = np.linalg.solve(Minv_sys, S_pE) if idx_E.size else np.zeros((M, n, 0), dtype=complex)
    S_EE_new = S_EE + S_Ep @ P @ solved

    new_modes = {q: block.modes_per_port[q] for q in E}
    return SMatrixBlock(wavelengths=wavelengths, ports=tuple(E), modes_per_port=new_modes, data=S_EE_new)
