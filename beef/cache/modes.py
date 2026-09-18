"""The length-independent factorization of an EME cascade.

meow's S-matrix factors cleanly into two parts (see
``meow.eme.cascade.compute_s_matrix_sax``, which is what ``mw.compute_s_matrix``
resolves to)::

    propagations = propagations_fn(modes, cells, cell_lengths=cell_lengths)
    interfaces   = interfaces_fn(modes)
    net          = _get_netlist(propagations, interfaces)

* ``interfaces`` depends ONLY on the modes -- overlap integrals between adjacent
  cross-sections. No length enters.
* ``propagations`` is where a cell's physical length enters, and
  ``compute_propagation_s_matrix()`` reads nothing from a mode but ``neff`` and
  ``env.wl``.

So the pair (interface S-matrices, per-cell neff) is a complete, length-
independent summary of a solved mode basis: it reproduces the cascade at ANY
``cell_lengths`` without re-solving. That pair is also small -- interface
matrices are [2n, 2n] per boundary and neffs are [n] per cell -- where the
expanded result grows with the length grid.

Both ``modes`` and ``cells`` are threaded to the injected callables and used
nowhere else in the cascade, so a cached factorization can be replayed with
placeholders in their place.
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)


def extract_factorization(modes: List[Sequence[Any]]) -> Dict[str, Any]:
    """Summarize a solved mode basis by what the cascade actually consumes.

    Args:
        modes: one mode list per cell, exactly as handed to
            ``mw.compute_s_matrix``.

    Returns:
        A picklable dict with ``interfaces`` (per-boundary ``(array, port_map)``
        as numpy), ``neffs`` (one complex array per cell) and ``wls`` (each
        cell's wavelength, so propagation can be rebuilt without a mode object).
    """
    import meow as mw

    interfaces = mw.compute_interface_s_matrices(modes)
    return {
        "interfaces": {
            key: (np.asarray(arr), dict(port_map))
            for key, (arr, port_map) in interfaces.items()
        },
        "neffs": [np.array([m.neff for m in cell_modes], dtype=complex)
                  for cell_modes in modes],
        "wls": [float(cell_modes[0].env.wl) for cell_modes in modes],
    }


def _propagation_s_dict(neff: np.ndarray, wl: float, cell_length: float) -> dict:
    """``compute_propagation_s_matrix()`` for one cell, from neff alone."""
    phases = np.exp(2j * np.pi * neff / wl * cell_length)
    s = {(f"left@{i}", f"right@{i}"): phases[i] for i in range(len(neff))}
    s.update({(p2, p1): v for (p1, p2), v in s.items()})
    return s


def s_matrix_from_factorization(
    fac: Dict[str, Any],
    cell_lengths: Sequence[float],
    sax_backend: str = "klu",
) -> Tuple[np.ndarray, dict]:
    """Replay the cascade at arbitrary ``cell_lengths`` from a factorization.

    No eigenmode solve and no meow CrossSection: the interface matrices come
    straight from the cache and the propagation part is rebuilt from neff.
    """
    import meow as mw

    neffs, wls = fac["neffs"], fac["wls"]
    if len(cell_lengths) != len(neffs):
        raise ValueError(
            f"len(cell_lengths) != cached cells: {len(cell_lengths)} != {len(neffs)}"
        )

    def interfaces_fn(_modes, **_kw):
        return fac["interfaces"]

    def propagations_fn(_modes, _cells=None, *, cell_lengths=None, **_kw):
        return {
            f"p_{i}": _propagation_s_dict(neffs[i], wls[i], float(cell_lengths[i]))
            for i in range(len(neffs))
        }

    # modes/cells are consumed only by the two callables above, both replaced.
    placeholder = [()] * len(neffs)
    return mw.compute_s_matrix(
        placeholder, cells=None, cell_lengths=list(cell_lengths),
        sax_backend=sax_backend,
        interfaces_fn=interfaces_fn, propagations_fn=propagations_fn,
    )
