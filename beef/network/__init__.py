"""beef.network: general S-matrix network composition (Redheffer star products).

See beef/network/ops.py's module docstring for the core math, and the
implementation plan (issued 2026-07-29) for the motivating gap this closes:
HybridSolver only cascades sections in a straight line, which cannot express
a ring resonator's feedback loop (validation_markdown/ring_resonator_validation.md)
or composing several independently-solved devices' SMatrix results into one
circuit (notes.md Solver recommendation #3).

Depends only on beef.results.smatrix.SMatrix and numpy — nothing in
beef.solvers or any example imports this yet (see the plan's Deferred section).
"""
from __future__ import annotations

from .block import SMatrixBlock, check_passive, rename_ports
from .circuit import Circuit
from .convert import from_smatrix, from_transfer_matrix, to_smatrix
from .ops import close_loop, connect

__all__ = [
    "SMatrixBlock",
    "rename_ports",
    "check_passive",
    "connect",
    "close_loop",
    "Circuit",
    "from_transfer_matrix",
    "from_smatrix",
    "to_smatrix",
]
