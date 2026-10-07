# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Georgia Tech Research Corporation
"""Spectral resampling of complex coefficients onto a dense wavelength grid.

Why this is a library function and not three lines inline
---------------------------------------------------------
A hybrid device is simulated at a modest number of wavelengths (the ring's
`N_WL_FDTD = 101`) and then evaluated on a dense grid to resolve resonances
(`N_WL_DENSE = 20001`). Everything between those two grids is interpolation,
and the choice of *what* to interpolate is not cosmetic:

    * **magnitude and unwrapped phase** (what `interp_complex` does) — a
      coefficient that is mostly a propagation phase, `A·exp(i·β·λ)`, has a
      slowly-varying magnitude and a phase that is nearly linear in λ. Linear
      interpolation of those two is near-exact.

    * **real and imaginary parts** — the same coefficient oscillates in Re and
      Im at the full phase rate. Linear interpolation chords across each
      oscillation, which *shrinks the magnitude* wherever the phase advances
      appreciably between samples.

The second failure is not hypothetical and not small. It was measured on the
ring on 2026-08-03: a Q-convergence table built with real/imag interpolation
showed Q never converging with grid refinement, because the interpolant itself
was damping the coefficients. Rebuilt with magnitude/unwrapped-phase, the same
data converges cleanly (FWHM flat at 0.1050 nm from 20001 points on). Dip
*positions* were unaffected, which is exactly what makes this dangerous — the
resonance wavelengths looked right the whole time.

`test_spectral.py` measures that difference on an analytic function with a
known answer rather than restating this paragraph as an assertion.

The unwrap caveat
-----------------
`np.unwrap` can only remove 2π jumps it can *see*: it assumes the true phase
advances by less than π between consecutive SOURCE samples. If the source grid
is too coarse for the coefficient's phase rate, unwrapping picks the wrong
branch and the interpolant is confidently wrong — worse than the real/imag
version, which merely blurs. `unwrap_is_safe()` measures that margin so a
caller can check it instead of hoping; `interp_complex(check_unwrap=...)` wires
it in directly.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

__all__ = ["interp_complex", "unwrap_is_safe"]


def unwrap_is_safe(c: np.ndarray, margin: float = 0.5) -> Tuple[bool, float]:
    """Flag a source grid that is MARGINAL for `np.unwrap`.

    Returns `(safe, worst_step)` where `worst_step` is the largest absolute
    phase advance between consecutive samples after unwrapping, in radians.

    **What this cannot do, stated plainly.** It cannot detect aliasing that has
    already happened. `np.unwrap` maps every step into `(-π, π]` by
    construction, so `worst_step <= π` *always*, no matter how badly
    undersampled the input is. A coefficient whose true advance is 5.0 rad per
    sample folds to an observed 1.28 rad and this function reports it safe —
    measured, not hypothesised: see
    `test_reports_safe_on_already_aliased_data` in tests/test_spectral.py,
    which pins that false negative rather than letting a future reader assume
    the guard covers it.

    That limit is fundamental, not an oversight. Recovering a phase rate above
    Nyquist from a single sampled sequence is undecidable — the aliased ramp is
    a genuine, self-consistent solution, and it stays self-consistent under
    decimation, so no data-only refinement test rescues it either. Deciding it
    needs information from outside the array: a physical phase rate (`β·L`
    from the geometry) or a run on a finer source grid.

    What it *is* good for: catching the regime where the sampling is close
    enough to the branch boundary that a modest change in the device — a longer
    section, a band edge, a higher-index mode — would push it over. That is a
    real and common way to break a working script, and it is cheap to check.

    Args:
        c: complex coefficient on the source grid.
        margin: fraction of π to treat as the safe ceiling. The default 0.5
            (i.e. π/2) leaves a factor-of-two headroom before the true break.
    """
    c = np.asarray(c)
    if c.size < 2:
        return True, 0.0
    steps = np.abs(np.diff(np.unwrap(np.angle(c))))
    worst = float(np.max(steps))
    return bool(worst < margin * np.pi), worst


def interp_complex(
    wl_src: np.ndarray,
    c: np.ndarray,
    wl_dense: np.ndarray,
    *,
    check_unwrap: Optional[float] = None,
) -> np.ndarray:
    """Resample complex `c(wl_src)` onto `wl_dense` via magnitude and
    unwrapped phase, each interpolated linearly.

    This is the resampling both ring resonator examples use for their coupler
    coefficients (`t`, `k1`, `k2`, `t_ring`, `b`). See the module docstring for
    why it is not real/imag.

    Args:
        wl_src: source wavelengths, increasing. Length must match `c`.
        c: complex coefficient at `wl_src`.
        wl_dense: target wavelengths. Values outside `wl_src`'s span are
            clamped by `np.interp`'s edge behaviour, NOT extrapolated — a
            deliberately conservative choice, since extrapolating an unwrapped
            phase ramp compounds error linearly in distance.
        check_unwrap: if given, the `margin` for `unwrap_is_safe()`; raises
            `ValueError` when the source grid is marginal. Off by default so
            this stays a drop-in for the plain three-line form, but worth
            enabling on any new call site. Read `unwrap_is_safe()`'s docstring
            before trusting it: it catches near-misses, NOT already-aliased
            input, which no data-only check can catch.

    Returns:
        Complex array shaped like `wl_dense`.
    """
    wl_src = np.asarray(wl_src, dtype=float)
    c = np.asarray(c)
    wl_dense = np.asarray(wl_dense, dtype=float)
    if wl_src.shape != c.shape:
        raise ValueError(
            f"interp_complex: wl_src has shape {wl_src.shape} but c has "
            f"{c.shape}; they must match.")
    if wl_src.size < 2:
        raise ValueError(
            f"interp_complex needs at least 2 source points, got {wl_src.size}.")
    if np.any(np.diff(wl_src) <= 0):
        raise ValueError(
            "interp_complex: wl_src must be strictly increasing (np.interp and "
            "np.unwrap both assume it; an unsorted grid produces a silently "
            "wrong interpolant rather than an error).")

    if check_unwrap is not None:
        safe, worst = unwrap_is_safe(c, margin=check_unwrap)
        if not safe:
            raise ValueError(
                f"interp_complex: phase advances by up to {worst:.3f} rad "
                f"({worst / np.pi:.2f}pi) between source samples, above the "
                f"{check_unwrap:.2f}pi margin. The unwrap branch is close "
                f"enough to the boundary that a modest change in the device "
                f"would flip it — refine wl_src rather than interpolating "
                f"this coefficient. (Note this guard cannot see aliasing that "
                f"has ALREADY folded; see unwrap_is_safe().)")

    mag = np.interp(wl_dense, wl_src, np.abs(c))
    phase = np.interp(wl_dense, wl_src, np.unwrap(np.angle(c)))
    return mag * np.exp(1j * phase)
