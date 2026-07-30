"""Circuit: stateful multi-block network builder on top of connect()/close_loop().

Matches ParallelSweep (beef/sweep/parallel.py) as this repo's existing
precedent for a stateful multi-step orchestrator implemented as a class
rather than a free-function pipeline.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from .block import SMatrixBlock, rename_ports
from .ops import Number, close_loop, connect

_ConnectOp = Tuple[str, str, str, str, str]          # ("connect", name_a, port_a, name_b, port_b)
_CloseLoopOp = Tuple[str, str, str, str, object]      # ("close_loop", name, port_out, port_in, loop)


class Circuit:
    """Accumulates connect()/close_loop() calls across several named blocks
    and resolves them all in one reduce() call.

    Each block registered via add_block() has its ports namespaced as
    f"{name}.{port}" immediately, so two instances of the same device (e.g.
    two rings) can never collide on port names. connect()/close_loop() here
    take the block's ORIGINAL (unnamespaced) port names — the namespacing
    is an implementation detail resolved internally.
    """

    def __init__(self) -> None:
        self._blocks: Dict[str, SMatrixBlock] = {}
        self._ops: List[Union[_ConnectOp, _CloseLoopOp]] = []

    def add_block(self, name: str, block: SMatrixBlock) -> None:
        if name in self._blocks:
            raise ValueError(f"Circuit.add_block: name {name!r} is already registered.")
        self._blocks[name] = rename_ports(block, {p: f"{name}.{p}" for p in block.ports})

    def _check_registered_port(self, name: str, port: str) -> None:
        if name not in self._blocks:
            raise ValueError(f"Circuit: block name {name!r} is not registered (call add_block first).")
        namespaced = f"{name}.{port}"
        if namespaced not in self._blocks[name].ports:
            original_ports = [p.split(".", 1)[1] for p in self._blocks[name].ports]
            raise ValueError(
                f"Circuit: port {port!r} not found on block {name!r} "
                f"(has ports {original_ports})."
            )

    def connect(self, name_a: str, port_a: str, name_b: str, port_b: str) -> None:
        """Queue a connection between two registered blocks' own port names.

        Not resolved until reduce(). If name_a and name_b have already
        been merged into the same component by an earlier op (i.e. this
        connection closes a cycle), reduce() dispatches to close_loop()
        instead of connect() automatically.
        """
        self._check_registered_port(name_a, port_a)
        self._check_registered_port(name_b, port_b)
        self._ops.append(("connect", name_a, port_a, name_b, port_b))

    def close_loop(
        self, name: str, port_out: str, port_in: str,
        loop: Optional[Union[Number, np.ndarray]] = None,
    ) -> None:
        """Queue a self-loop closure on one registered block's own ports."""
        self._check_registered_port(name, port_out)
        self._check_registered_port(name, port_in)
        self._ops.append(("close_loop", name, port_out, port_in, loop))

    def reduce(self, external_ports: Optional[Dict[str, str]] = None) -> SMatrixBlock:
        """Resolve every queued op in declaration order and return the
        single fully-reduced block.

        Args:
            external_ports: optional {"name.port": "final_name"} renaming
                applied to the result. Ports not listed keep their internal
                "name.port" label.

        Raises:
            ValueError: no blocks/ops registered, or more than one
                disconnected component remains after every queued op is
                applied (message lists the disjoint name groups).
        """
        if not self._blocks:
            raise ValueError("Circuit.reduce: no blocks registered (call add_block first).")

        component: Dict[str, SMatrixBlock] = dict(self._blocks)

        for op in self._ops:
            if op[0] == "connect":
                _, name_a, port_a, name_b, port_b = op
                block_x, block_y = component[name_a], component[name_b]
                ns_a, ns_b = f"{name_a}.{port_a}", f"{name_b}.{port_b}"
                if block_x is block_y:
                    # name_a and name_b already share a component: this
                    # "connect" actually closes a cycle within one block.
                    merged = close_loop(block_x, ns_a, ns_b, loop=None)
                else:
                    merged = connect(block_x, ns_a, block_y, ns_b)
                for nm, blk in component.items():
                    if blk is block_x or blk is block_y:
                        component[nm] = merged
            else:  # "close_loop"
                _, name, port_out, port_in, loop = op
                block_x = component[name]
                ns_out, ns_in = f"{name}.{port_out}", f"{name}.{port_in}"
                merged = close_loop(block_x, ns_out, ns_in, loop=loop)
                for nm, blk in component.items():
                    if blk is block_x:
                        component[nm] = merged

        distinct = {id(blk): blk for blk in component.values()}
        if len(distinct) > 1:
            groups: Dict[int, List[str]] = {}
            for nm, blk in component.items():
                groups.setdefault(id(blk), []).append(nm)
            group_strs = [f"{{{', '.join(names)}}}" for names in groups.values()]
            raise ValueError(
                f"Circuit.reduce: {len(distinct)} disconnected components "
                f"remain — connect them all before reducing: {group_strs}."
            )

        final = next(iter(distinct.values()))
        if external_ports:
            final = rename_ports(final, external_ports)
        return final
