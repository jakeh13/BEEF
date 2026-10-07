# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Georgia Tech Research Corporation
"""`interface_modes` convergence: the device-agnostic completeness check.

`HybridSolver.solve()` contracts index *i* of one section's output basis
against index *i* of the next section's input basis, truncated to
``interface_modes[k]`` modes at interface *k*. Three separate properties have
to hold for that contraction to be correct (see the doc's "The requirement"
section): (a) identity/ordering, (b) completeness — the truncated basis
actually captures the power crossing that plane — and (c) distinctness/rank.
Layer 1 (``beef.mode_identity``) and layer 4 (self-overlap controls) address
(a) and (c) by comparing modes that were ALREADY requested; neither can see a
mode that was never asked for in the first place. Layer 3
(``beef.cut_position_invariance``) is correctness-given-the-count, not
sufficiency-of-the-count. **Only a direct sweep of the mode count itself,
with the final S-matrix as the observable, establishes (b)** — which is the
only property that generalizes to a novel geometry, since nobody knows in
advance how many modes it needs.

This module is the reusable generalization item 6 asks for, in the same
spirit as ``beef.cut_position_invariance`` (item 7): it operates directly on
``HybridSolver.solve()``'s own ``SMatrix`` output for a sequence of
already-solved mode counts, so it works for ANY device/port structure with no
per-device code. Constructing each solve — which interface(s) get how many
modes, on which device, at which resolution — is genuinely device-specific
and stays out of this module, exactly as cut-position construction stays out
of ``cut_position_invariance``.

**Caller responsibility, not enforced here:** the mesh must already be
converged at the resolution the sweep is run at. A plateau found on an
unconverged mesh is a mesh artifact, not a statement about the physics (item
6's own text, echoing item 1 point 5) — this module has no way to know
whether that precondition holds and does not check it.
"""
from __future__ import annotations

from typing import Any, Dict, Sequence, Tuple

import numpy as np

from .results.smatrix import SMatrix

SMatrixKey = Tuple[str, int, int]


def interface_modes_convergence(
    smats: Sequence[SMatrix], mode_counts: Sequence[int],
) -> Dict[str, Any]:
    """Convergence of every ``(port, in_mode, out_mode)`` transmission
    coefficient across a sequence of already-solved interface-mode-count
    values.

    Args:
        smats: one ``SMatrix`` per mode count (already solved by the
            caller — this function runs no simulation). Must all carry the
            exact same set of ``(port, in_mode, out_mode)`` keys: they are
            claimed to be the same physical device solved with a different
            ``interface_modes`` truncation, so a key present in one but not
            another means the caller built an inconsistent sweep, not a
            property of the mode count.
        mode_counts: the interface mode count each ``smats[i]`` was solved
            with, same length as ``smats``. At least 2 values are required
            to compute a step at all.

    Unlike ``cut_position_invariance`` (whose correct answer is a horizontal
    line — the S-matrix must not depend on cut position at all, so a fitted
    slope near 0 is the criterion), a mode-count sweep is expected to MOVE
    while the basis is still incomplete and then plateau once it is not —
    there is no reason to expect the approach to be linear. So this reports
    the wavelength-mean transmission at each mode count and the **step**
    between consecutive counts (matching the step-based convergence
    convention already used elsewhere in this repo — e.g.
    ``mode_spectrum_convergence.py``, ``verify_item_8_dispersive_
    convergence.py`` — rather than fitting a trend line to a quantity with no
    reason to be linear).

    Returns:
        dict with:
            ``keys``: ``{(port, in_mode, out_mode): {"mode_counts", "values"
                (wavelength-mean transmission per mode count, sorted by
                mode_counts), "steps" (diff of values, length len-1),
                "last_step" (steps[-1] — has it stopped moving at the
                HIGHEST count tried, the diagnostic question this item asks),
                "max_step" (max |step| over the whole sweep, informational —
                shows whether early low-count points moved much more than
                recent ones, the expected shape of real convergence)}}``.
            ``worst_key``: the key with the largest-magnitude ``last_step``.
            ``worst_last_step``: that key's last step — the single number to
                watch as the sweep is extended to a higher mode count; it
                should shrink toward 0, the same criterion item 6 states
                ("require the final S-matrix to plateau").

    Raises:
        ValueError: fewer than 2 (smat, mode_count) pairs, mismatched
            lengths, or the supplied SMatrices do not all carry the same key
            set.
    """
    smats = list(smats)
    counts_arr = np.asarray(mode_counts, dtype=int)

    if len(smats) != len(counts_arr):
        raise ValueError(
            f"interface_modes_convergence: smats and mode_counts must have "
            f"the same length, got {len(smats)} and {len(counts_arr)}."
        )
    if len(smats) < 2:
        raise ValueError(
            f"interface_modes_convergence: need at least 2 (smat, "
            f"mode_count) pairs to compute a step, got {len(smats)}."
        )

    key_sets = [set(s.keys()) for s in smats]
    if any(ks != key_sets[0] for ks in key_sets[1:]):
        raise ValueError(
            "interface_modes_convergence: every supplied SMatrix must carry "
            "the same (port, in_mode, out_mode) keys -- they are claimed to "
            "be the same physical device solved at different interface_modes "
            f"counts. Got key sets: {[sorted(ks) for ks in key_sets]}."
        )

    order = np.argsort(counts_arr)
    counts_sorted = counts_arr[order]

    keys: Dict[SMatrixKey, Dict[str, Any]] = {}
    worst_key = None
    worst_last_step = 0.0
    for key in sorted(key_sets[0]):
        values = np.array(
            [float(np.mean(np.abs(s.coeff(*key)) ** 2)) for s in smats]
        )[order]
        steps = np.diff(values)
        last_step = float(steps[-1])
        max_step = float(np.max(np.abs(steps))) if len(steps) else 0.0
        keys[key] = dict(
            mode_counts=counts_sorted, values=values, steps=steps,
            last_step=last_step, max_step=max_step,
        )
        if worst_key is None or abs(last_step) > abs(worst_last_step):
            worst_key, worst_last_step = key, last_step

    return dict(keys=keys, worst_key=worst_key, worst_last_step=worst_last_step)
