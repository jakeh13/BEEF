"""Cut-position invariance: the authoritative, device-agnostic hybrid
correspondence check (``validation_markdown/fdtd_vs_EME_mode_validation.md``
layer 3, item 7).

A correctly-decomposed hybrid chain's S-matrix cannot depend on WHERE along a
uniform stretch of geometry the FDTD/EME interface is cut — the physical
device is identical for every cut position, so any movement in the answer as
the cut moves is method error, not a different device. This is the strongest
available correspondence test: unlike ``beef.mode_identity``'s field-overlap
machinery (layers 1/2/4), it needs no mode comparison, no shared diagnostic
grid, and no cross-engine field calibration, so it applies to any geometry at
any resolution pair. A broken correspondence *cannot* be cut-invariant, which
is why it subsumes those layers as a pass/fail criterion (its weakness is
diagnostic resolution, not validity — it says something is wrong, not which
mode). It is already validated as the meaningful property on real data: the
directional coupler's stub-sweep drift collapsed from +0.1085 to −0.0007
µm/µm (167x, below run scatter) once its EME body was β-calibrated.

That result came from
``examples/directional_coupler/python/directional_coupler_hybrid_stub_sweep.py``,
a one-off script hard-coded to one device's geometry, boundary bookkeeping,
and a single scalar observable (the 50:50 coupling-length crossing point).
This module is the reusable generalization item 7 asks for: it operates
directly on ``HybridSolver.solve()``'s own ``SMatrix`` output for a sequence
of already-solved stub values, so it works for ANY device/port structure with
no per-device code. Constructing each stub's ``HybridSolver``/geometry/
boundaries is left to the caller — how a stub trades length between FDTD and
EME sections, at how many interior interfaces, for which physical device, is
genuinely device-specific and is exactly what stays out of this module.
"""
from __future__ import annotations

from typing import Any, Dict, Sequence, Tuple

import numpy as np

from .results.smatrix import SMatrix

SMatrixKey = Tuple[str, int, int]


def cut_position_invariance(
    smats: Sequence[SMatrix], stubs: Sequence[float],
) -> Dict[str, Any]:
    """Drift of every ``(port, in_mode, out_mode)`` transmission coefficient
    across a sequence of already-solved cut-position (stub) values.

    Args:
        smats: one ``SMatrix`` per stub value (already solved by the
            caller — this function runs no simulation). Must all carry the
            exact same set of ``(port, in_mode, out_mode)`` keys: they are
            claimed to be the same physical device cut at different
            positions, so a key present in one but not another means the
            caller built inconsistent stubs, not a property of the sweep.
        stubs: the µm stub value each ``smats[i]`` was solved at, same
            length as ``smats``. At least 2 values are required to fit a
            drift slope at all.

    For each key, transmission (``|coeff|^2``, averaged over wavelength —
    the same reduction ``directional_coupler_hybrid_stub_sweep.py``'s own
    50:50-crossing curve uses) is fit against `stubs` with a straight line.
    The fitted **slope** — drift per µm of stub — is the headline metric
    this item asks for, since a converged interface placement is one where
    the answer stops moving as the cut moves further from the discontinuity
    entirely, not merely one where two particular stubs happen to agree
    (which the fitted `intercept`/`max_drift` alone cannot distinguish from
    a slope that is merely small over the sampled range). Only the
    wavelength-MEAN transmission is fit; a caller wanting the wavelength-
    resolved drift (e.g. to check whether it is worse at the band edges)
    should slice ``smats``/re-call per wavelength — deliberately out of
    scope here, matching how the one-off script itself only ever plotted
    the mean-over-wavelength curve.

    Returns:
        dict with:
            ``keys``: ``{(port, in_mode, out_mode): {"stubs", "values"
                (wavelength-mean transmission per stub), "slope",
                "intercept", "max_drift" (max - min over the swept stubs)}}``.
            ``worst_key``: the key with the largest-magnitude `slope`.
            ``worst_slope``: that key's slope — the single number to watch
                as the sweep is extended; it should shrink toward 0 as the
                stub grows, the same convergence
                ``directional_coupler_hybrid_stub_sweep.py`` demonstrated for
                the 50:50 split length once the EME body was β-calibrated.

    Raises:
        ValueError: fewer than 2 stubs/smats, mismatched lengths, or the
            supplied SMatrices do not all carry the same key set.
    """
    smats = list(smats)
    stubs_arr = np.asarray(stubs, dtype=float)

    if len(smats) != len(stubs_arr):
        raise ValueError(
            f"cut_position_invariance: smats and stubs must have the same "
            f"length, got {len(smats)} and {len(stubs_arr)}."
        )
    if len(smats) < 2:
        raise ValueError(
            f"cut_position_invariance: need at least 2 (smat, stub) pairs to "
            f"fit a drift slope, got {len(smats)}."
        )

    key_sets = [set(s.keys()) for s in smats]
    if any(ks != key_sets[0] for ks in key_sets[1:]):
        raise ValueError(
            "cut_position_invariance: every supplied SMatrix must carry the "
            "same (port, in_mode, out_mode) keys -- they are claimed to be "
            "the same physical device cut at different stub positions. Got "
            f"key sets: {[sorted(ks) for ks in key_sets]}."
        )

    order = np.argsort(stubs_arr)
    stubs_sorted = stubs_arr[order]

    keys: Dict[SMatrixKey, Dict[str, Any]] = {}
    worst_key = None
    worst_slope = 0.0
    for key in sorted(key_sets[0]):
        values = np.array(
            [float(np.mean(np.abs(s.coeff(*key)) ** 2)) for s in smats]
        )[order]
        slope, intercept = np.polyfit(stubs_sorted, values, 1)
        max_drift = float(values.max() - values.min())
        keys[key] = dict(
            stubs=stubs_sorted, values=values,
            slope=float(slope), intercept=float(intercept), max_drift=max_drift,
        )
        if worst_key is None or abs(slope) > abs(worst_slope):
            worst_key, worst_slope = key, float(slope)

    return dict(keys=keys, worst_key=worst_key, worst_slope=worst_slope)
