"""FDTDSolver: meep-based FDTD solver for MMI and related photonic devices."""
from __future__ import annotations

import logging
import os
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Literal, NamedTuple, Optional, Sequence, Tuple

import numpy as np

from ..cache import (
    CacheMissError, build_run_config, build_whole_device_config, load_result, lookup,
    lookup_superset, save_run, update_family_pointer,
)
from ..geometry.base import DeviceGeometry
from ..materials.stack import MaterialStack, resolved_material_config, SOI_160
from ..results.smatrix import SMatrix
from .._log import get_logger

logger = get_logger(__name__)

# Starting guess for section_mode_neff()'s eigensolves, as an effective index.
# Only a seed — meep's match_frequency iterates from here to the true k — so it
# just has to sit in the guided range for a Si-core/oxide-clad waveguide
# (1.44 < n_eff < 3.48). Roughly the middle of that window.
_KPOINT_SEED_NEFF = 2.5

# --- _iter_distinct_bands()'s duplicate-eigenmode guard ---------------------
#
# See validation_markdown/modal_multiplexer_validation.md gap 6: sim.get_eigenmode()
# can return the SAME eigenvector for two different requested band_num values.
# Measured directly on the reproduction case: true duplicates overlap at
# 0.9999999999+; every genuinely distinct band pair on the same cross-section
# overlapped <= 0.1257. 0.9 sits with wide margin on both sides of that gap.
_DUP_BAND_OVERLAP_THRESHOLD = 0.9
# Diagnostic-only grid for the duplicate check -- deliberately small and fixed,
# independent of self.res or any caller-requested field grid, since it only
# needs to distinguish "same eigenvector" from "different eigenvector".
_DUP_BAND_DIAG_NY = 16
_DUP_BAND_DIAG_NZ = 8
# Extra get_eigenmode() calls tolerated beyond n_modes before giving up.
# get_eigenmode() is a "seconds" operation (see section_mode_neff()'s
# docstring) and every duplicate seen in practice needed only one retry, so
# this stays cheap even at its ceiling.
_DUP_BAND_MAX_EXTRA = 8

# --- solve_section() excite-normalization mismatch guard -------------------
#
# See validation_markdown/modal_multiplexer_validation.md gap 7: the connector
# section's excited bands 1-2 (Δn_eff=0.0148, near-degenerate) produced
# transfer-matrix columns 10-20x too large -- get_eigenmode_coefficients()'s
# alpha at the norm-flux plane silently undercounted launched power for those
# two bands specifically, while mp.get_fluxes() (a direct field integral, not
# a modal decomposition) stayed correct. Measured mismatch for the broken
# bands was 10-20x; every other band on the same cross-section agreed with
# flux to within numerical noise. 0.1 sits with wide margin below the former
# and above the latter.
_EXCITE_NORM_MISMATCH_TOL = 0.1

# --- section_mode_beta() safety thresholds ---------------------------------
#
# REWRITTEN 2026-08-01. This block used to describe "failure mode 1" as the
# residual guard's job: a seed near a branch boundary made the per-length wrap
# counts inconsistent, so the points stopped lying on a line and the residual
# blew up. That was an accurate description of a REAL BUG IN THIS MODULE, not
# of a hazard callers had to live with. `beta_from_phases` rounded the wrap
# count independently at every length against one seed slope, with no
# cross-length consistency constraint, so a seed error ε produced a
# fractional-turn offset ε·L/(2π) that GREW with L until some length straddled
# a .5 boundary and got a turn count inconsistent with its neighbours.
#
# It bit for real: job 11589770 (mmi_1x2, calibrate_body_neff=True) died on a
# 2.101 rad residual, and the failure was recorded for a day as "a π sign flip
# in T_ii" -- a reading that fits the same residual pattern exactly, because
# for three equally spaced points "middle point displaced by −π" and "points 2
# and 3 each displaced by −2π" are algebraically indistinguishable (both give
# δ·(+1/3, −2/3, +1/3)). The turn-count reading is the true one: it returns
# n_eff 2.783, +1.96% over MPB and in family with the +2.233%/+2.361% measured
# on the coupler, where the sign-flip reading returned 2.101, −23% off, wrong
# in sign as well as size. Every one of that job's 88 (wavelength, mode) cells
# is a clean linear ramp once whole turns are assigned consistently: the second
# difference (φ₃ − 2φ₂ + φ₁)/2π is within 0.0031 of an exact integer in all 88.
#
# The sting: ε there was 1.96% -- the very MPB-vs-time-stepper gap that
# `calibrate_body_neff` exists to measure. **The unwrapping was broken by
# exactly the quantity it was built to measure.** Measured failure rate over
# seed ±3%: 550/4025 (13.7%) on that device's cached data, and 104/2475 (4.2%)
# on the coupler's -- so the coupler's success was luck, not margin. At
# LENGTHS=[2,4,6,8] the old code's tolerable seed error was 0.000 for an
# adverse intercept.
#
# `beta_fit_from_phases` now decides ONE integer against the seed -- the branch
# -- and derives every later turn from a slope refit on the shorter lengths, so
# the ε·L drift cannot accumulate. That NARROWS the seed's authority from n
# decisions to 1, and 1 is exactly what guard 2 below already licenses.
#
# What remains are two INDEPENDENT failure modes, neither catching the other:
#
#   1. Seed is a WHOLE branch off -> every length agrees on that wrong branch,
#      so the fit is perfectly clean (residual 1e-14) at a β 33% wrong, and the
#      margin-to-seed is 0.000. Nothing in the fitted numbers can detect this.
#      The only defence is a precondition on the LENGTHS: keep the branch
#      spacing 2π/ΔL large compared with any physically possible FDTD-vs-MPB
#      disagreement, so a whole-branch error cannot happen in the first place.
#      Guard: _BETA_MIN_BRANCH_FRACTION.
#   2. Seed lands within a hair of a branch BOUNDARY -> the one seed-dependent
#      decision is a coin flip, and no arithmetic can rescue it because the
#      seed genuinely carries no information there. The residual is now CLEAN
#      in this case (the fit is self-consistent on whichever branch it lands),
#      so the residual cannot be the tell any more. Guard:
#      _BETA_BRANCH_MARGIN_RAISE.
#
# Minimum branch spacing as a fraction of β, for guard 1. The measured
# disagreement is ~2%; requiring the branches to be >=25% of β apart leaves an
# order of magnitude. At the coupler's β = 9.6 rad/µm this caps the length
# spacing at 2.6 µm, which the item 6 measurement's 2.0 µm satisfies.
_BETA_MIN_BRANCH_FRACTION = 0.25

# Max phase-fit residual (rad). Post-2026-08-01 this means exactly ONE thing:
# the phase is not linear in L -- a non-uniform cut, a near field that has not
# settled, mode mixing, or band contamination. Turn-count ambiguity no longer
# contributes, so this no longer double-duties as a seed-quality check.
#
# Clean values measured: 1e-14 synthetic, 1.4e-5 on the 2-mode item-6 data,
# 6.45e-3 worst of 88 cells on mmi_1x2's res-40 data, and 2.2e-2 worst of 66
# cells on the REAL 6-mode coupler production cache. That last one is only
# 4.5x under the threshold -- thinner than the "1.4e-5" this comment used to
# quote, which was the 2-mode configuration only. Kept at 0.1 deliberately:
# tightening it toward the observed worst clean value would start failing real
# runs.
_BETA_FIT_RESIDUAL_MAX = 0.1

# How far the answer moved from the seed, as a fraction of the half-distance to
# the neighbouring branch (see beta_fit_from_phases). Two levels:
#
#   MAX (warn)   the coupler sits at 0.13, mmi_1x2 at 0.078. Above this,
#                say so in the log -- worth a human look, not wrong per se.
#   RAISE (hard) the seed sat within 10% of a branch boundary, so the single
#                seed-dependent decision was near-arbitrary. This makes
#                section_mode_beta's own docstring (">= 1.0 means the seed
#                could not distinguish the branches") enforced for the first
#                time -- it claimed a ValueError it never raised. 0.9 rather
#                than 1.0 to catch the case just before it degenerates, and
#                rather than 0.5 because promoting a heuristic look-at-this
#                level to a hard failure would break live callers for no
#                physics reason. Note that for equally spaced lengths the new
#                algorithm guarantees margin <= 1 by construction, so this
#                fires only for a genuinely near-boundary seed or for unequal
#                spacing whose extrapolation is bad -- both real caller errors.
_BETA_BRANCH_MARGIN_MAX = 0.5
_BETA_BRANCH_MARGIN_RAISE = 0.9

# Warn if FDTD's β exceeds the MPB seed by more than this fraction. The whole
# point of the calibration is that this excess is nonzero and physical:
# +2.233%/+2.361% on the coupler, +1.96% on mmi_1x2. An excess past 10% is out
# of family and means something other than discretization is going on. Warn
# only -- there is no principled hard limit, and this quantity is the
# measurement, not an error.
_BETA_SEED_EXCESS_WARN = 0.10

# Cap on beta_fit_from_phases' fixed-point polish iterations. The sequential
# pass is already self-consistent in every case measured (~7000 trials
# converged in 0 extra iterations); the loop exists only so a noisy two-point
# bootstrap cannot leave a stale turn vector behind. Exhausting it RAISES
# rather than returning the last iterate -- silently returning an unconverged
# turn vector is precisely the class of bug this rewrite removes.
_BETA_TURN_POLISH_MAX = 8


def monitor_freq_grid(wavelengths: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """DFT-monitor frequencies for `wavelengths`, plus the index map that puts
    meep's results back into the caller's wavelength order.

    Every monitor in this module must be built from this explicit frequency
    list, never from meep's ``(fcen, df, nfreq)`` triple. That triple spreads
    nfreq points evenly over ``[fcen - df/2, fcen + df/2]`` *in frequency*, so
    it cannot reproduce a wavelength array that is linear in λ, and it orders
    them ascending in frequency — i.e. **descending in wavelength**. Indexing
    such a result against an ascending-λ array mislabels every point but the
    centre one, silently: with the ``df = (1/λ_min - 1/λ_max) * 0.5`` this
    module used to pass, a 1.50–1.60 µm 11-point sweep actually sampled
    1.5254–1.5754 µm — half the requested span, reversed, up to 75 nm off. It
    made every FDTD result here report a directional coupler's coupling as
    *weakening* with wavelength (backwards), and made HybridSolver chains
    contract an EME section at one wavelength against an FDTD section at a
    different one. See dc_validation.md.

    Returns:
        mon_freqs: strictly ascending frequencies (1/µm), for
            ``sim.add_mode_monitor(mon_freqs, ...)``.
        wl_order: index array such that ``alpha[wl_order]`` reorders a
            monitor-ordered result back onto `wavelengths` as given.
    """
    freqs = 1.0 / np.asarray(wavelengths, dtype=float)
    order = np.argsort(freqs)              # ascending-freq position -> caller index
    return freqs[order], np.argsort(order)  # caller index -> ascending-freq position


# Relative bandwidth used when a caller asks for a single wavelength, where the
# requested span is zero and there is nothing to scale: 0.1 * fcen before the
# factor, i.e. meep's usual fwidth = 0.2 * fcen at the default factor of 2.
_SINGLE_WL_REL_BANDWIDTH = 0.1


def source_fwidth(wavelengths: np.ndarray, factor: float) -> float:
    """Gaussian-pulse `fwidth` covering the requested wavelength span.

    `factor` multiplies the *full* frequency span ``1/λ_min - 1/λ_max``. It has
    to exceed 1: meep's Gaussian falls as ``exp(-2π²(Δf/fwidth)²)``, so at
    ``factor = 1`` the band edges sit at 0.005 of peak amplitude and at
    ``factor = 1.5`` at 0.097 — the default 2.0 puts them at 0.27, comfortably
    above the 10% floor `check_source_coverage()` warns below. The ratio is
    what matters, so this holds for any band or point count.

    Before the frequency-grid fix (see `monitor_freq_grid`) this was
    ``span * 0.5``, which looked adequate only because the monitors shared the
    same value and therefore sampled the same too-narrow band.
    """
    wl = np.asarray(wavelengths, dtype=float)
    span = 1.0 / wl.min() - 1.0 / wl.max()
    if span <= 0.0:                        # single wavelength — nothing to span
        return factor * _SINGLE_WL_REL_BANDWIDTH / float(wl.mean())
    return factor * span


def dbeta_from_phases(lengths, dphi) -> Tuple[float, float, float]:
    """Differential phase rate (rad/µm) from ``dphi`` measured at several lengths.

    ``dphi(L) = arg(T_ii) - arg(T_jj) = Δβ·L + const``, where the constant is
    mode-matching and monitor-plane phase, identical at every length. The SLOPE
    is the physics; the intercept is discarded. One length cannot do this —
    ``dphi(L)/L`` would carry the offset as a 1/L error.

    Well-conditioned and needs no prior, because Δβ is small: at the coupler's
    0.095 rad/µm a 2 µm step accrues 0.19 rad, nowhere near a wrap. That is
    exactly what `beta_from_phases()` below cannot assume, and why the two are
    separate functions.

    Returns ``(dbeta, intercept, rms_residual)``. The residual is the evidence
    that dphi really is linear in L, i.e. that this is one clean mode pair.
    """
    lengths = np.asarray(lengths, dtype=float)
    y = np.unwrap(np.asarray(dphi, dtype=float))
    A = np.vstack([lengths, np.ones_like(lengths)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = float(np.sqrt(np.mean((y - (slope * lengths + intercept)) ** 2)))
    return float(slope), float(intercept), resid


class BetaFit(NamedTuple):
    """Everything `beta_fit_from_phases` learned, not just the answer.

    `beta`/`resid`/`margin` are the historical 3-tuple, in order, so
    ``beta_from_phases`` can stay a thin facade over this.
    """

    beta: float             # fitted phase rate, rad/µm
    resid: float            # max |deviation| from the fitted line, rad
    margin: float           # |beta - seed| / half-distance to next branch
    seed_excess: float      # beta/seed_beta - 1. THE PHYSICS: what calibration measures.
    turns: np.ndarray       # chosen wrap counts, gauge-normalised to turns[0] == 0
    seed_turns: np.ndarray  # what per-length rounding would have given, same gauge
    turns_shifted: bool     # turns != seed_turns, i.e. the old code would have differed
    branch_ambiguity: float # |x - round(x)| of the ONE seed-dependent decision, [0, 0.5]
    turn_ambiguity: float   # worst of the DATA-determined decisions, [0, 0.5]


class BetaCalibration(NamedTuple):
    """Everything one `section_mode_beta()` call measured, not just its answer.

    `section_mode_beta()` returns only `beta`, which is all its callers need
    to propagate a body on. Plotting and diagnostics need the *inputs* to
    that answer too — the fit lengths, the raw per-length phases, the MPB
    seed each branch was chosen against, and the per-cell `BetaFit` — so
    `section_mode_beta_details()` returns this and `section_mode_beta()` is
    a one-line facade over it.

    Exists to keep that data in ONE place. Before it, every consumer
    (`diagnose_mmi1x2_calibration_phase_fit.py`, the three
    `resolved/section_mode_beta_multi_monitor_validation.md` step 4-6 scripts)
    re-derived the length selection by hand against the source; five copies
    of a formula that has to agree exactly is how `section_fdtd_layout()`'s
    own duplicate drifted (see that function's docstring).

    Attributes:
        beta: ``[n_wl, n_modes]`` rad/µm — exactly what `section_mode_beta()`
            returns.
        lengths: the fit lengths actually used (µm), in the order simulated.
        phases: ``[n_L, n_wl, n_modes]`` wrapped ``arg(T_ii)``, i.e. the raw
            measurement before any wrap counting.
        wavelengths: ``[n_wl]`` µm, as passed in.
        seed_wavelengths: the wavelengths the seed was actually solved at.
            Equal to `wavelengths` unless `seed_anchors` subsampled it, in
            which case `seed_beta` is interpolated from these.
        seed_beta: ``[n_wl, n_modes]`` rad/µm, MPB's answer — the branch
            prior, never the result.
        fits: ``fits[i][m]`` is the full `BetaFit` for wavelength i, mode m,
            including the chosen `turns` a plot needs to unwrap `phases`.
        x0, x1, mid: the device-frame span this calibrated, and its midpoint.
        merged: True if the phases came from one shared-domain
            `solve_section_multi_length()` call rather than `len(lengths)`
            independent `solve_section()` runs (see `use_merged_domain`).
        guard_failures: post-fit guard messages that would have raised under
            `strict=True`. **Empty on any calibration usable for
            propagation** — non-empty means the caller passed `strict=False`
            and `beta` is for inspection only. A tuple, not a list, so the
            default cannot be a shared mutable.
    """

    beta: np.ndarray
    lengths: List[float]
    phases: np.ndarray
    wavelengths: np.ndarray
    seed_beta: np.ndarray
    fits: List[List[BetaFit]]
    x0: float
    x1: float
    mid: float
    merged: bool
    guard_failures: Tuple[str, ...] = ()
    seed_wavelengths: Optional[np.ndarray] = None


def _fit_line(lengths, values) -> Tuple[float, float]:
    A = np.vstack([lengths, np.ones_like(lengths)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, values, rcond=None)
    return float(slope), float(intercept)


def _diagonal_of(T: np.ndarray) -> np.ndarray:
    """`T` with every off-diagonal entry replaced by NaN, shape unchanged.

    The return shape of a `solve_section(measure_diagonal_only=True)` run,
    and the conversion that lets a cached FULL matrix serve one. NaN rather
    than 0 so a consumer reading an entry this run never measured fails
    loudly instead of taking "not computed" for "measured zero" — the whole
    point of the flag is that those entries are unknown, not small.

    Handles both `solve_section()` shapes: `[n_wl, n_out, n_in]` and
    `measure_margins`' `[n_wl, n_pos, n_out, n_in]`.
    """
    out = np.full_like(T, np.nan + 0j)
    n = min(T.shape[-2], T.shape[-1])
    idx = np.arange(n)
    out[..., idx, idx] = T[..., idx, idx]
    return out


def beta_fit_from_phases(lengths, phases, seed_beta) -> BetaFit:
    """Absolute phase rate (rad/µm) for ONE mode, from ``arg(T_ii)`` vs length.

    Unlike `dbeta_from_phases()` this is a *large* rate — β ≈ 9.4 rad/µm for a
    silicon strip — so the phase wraps many times between lengths and
    ``np.unwrap`` is useless: at the coupler's 2 µm spacing, ``β·ΔL`` is 2.9993
    full turns, and all that survives is the 0.002 rad remainder. Getting this
    wrong is not subtle but it is silent — a naive fit here returned Δβ = −0.331
    against a true 0.08525 during the item 6 analysis.

    **The seed decides exactly one integer, and it is the branch.** Everything
    else is read off the data. Concretely: the first gap's wrap increment is
    rounded against `seed_beta`; every later length's increment is rounded
    against a slope *refit from all shorter lengths*. That refit is the whole
    point — see the block comment above `_BETA_MIN_BRANCH_FRACTION`. Until
    2026-08-01 this function rounded every length independently against the one
    seed slope, so a seed error ε drifted by ε·L/(2π) and eventually straddled a
    .5 boundary. A 1.96% seed error — the FDTD-vs-MPB gap that
    `calibrate_body_neff` exists to *measure* — was enough to do it, which is
    how job 11589770 died on a 2.101 rad residual.

    Choosing the branch against the seed is safe rather than circular: with
    equally spaced lengths the candidate slopes are ``2π/ΔL`` apart — 3.14 rad/µm
    at ΔL = 2 µm, a 33% step — while the disagreement being measured is of order
    2%. `_BETA_MIN_BRANCH_FRACTION` is the precondition that keeps that true.

    Args:
        lengths: section lengths (µm). Need not be equally spaced or sorted;
            sorted internally, but `turns`/`seed_turns` come back in the order
            given. The wrap-safety margin is governed by the largest gap.
        phases: ``arg(T_ii)`` at each length, each wrapped into (−π, π].
        seed_beta: prior for branch selection (rad/µm), normally MPB's answer
            from `FDTDSolver.section_mode_neff()`.

    Returns:
        `BetaFit`. Read `resid` as "is this a straight line" and
        `branch_ambiguity`/`margin` as "was the seed good enough to pick the
        branch" — post-2026-08-01 those are separate questions, and the residual
        answers only the first.

    Raises:
        ValueError: if fewer than two lengths are given, or if the turn vector
            fails to reach a fixed point within `_BETA_TURN_POLISH_MAX`.
    """
    lengths_in = np.asarray(lengths, dtype=float)
    phases_in = np.asarray(phases, dtype=float)
    seed_beta = float(seed_beta)
    n = lengths_in.size
    if n < 2:
        raise ValueError(
            f"beta_from_phases: need at least 2 lengths to fit a rate, got {n}"
        )

    order = np.argsort(lengths_in)          # work ascending; report as given
    L = lengths_in[order]
    phi = phases_in[order]
    two_pi = 2.0 * np.pi

    # What the pre-2026-08-01 code would have done, kept for the record. Gauge
    # both this and the chosen vector to [0] == 0: only DIFFERENCES of turns are
    # physical (a constant is absorbed by the fitted intercept), so comparing
    # them un-normalised would report spurious shifts.
    seed_turns = np.round((seed_beta * L - phi) / two_pi)
    seed_turns = seed_turns - seed_turns[0]

    turns = np.zeros(n)
    # THE one seed-dependent integer: the first gap's increment, i.e. the branch.
    x = (seed_beta * (L[1] - L[0]) - (phi[1] - phi[0])) / two_pi
    turns[1] = np.round(x)
    branch_ambiguity = float(abs(x - turns[1]))

    # Every later turn is decided against a slope REFIT from the shorter
    # lengths, so the seed's error never gets multiplied by a longer L.
    turn_ambiguity = 0.0
    for j in range(2, n):
        slope, intercept = _fit_line(L[:j], phi[:j] + two_pi * turns[:j])
        y = (slope * L[j] + intercept - phi[j]) / two_pi
        turns[j] = np.round(y)
        turn_ambiguity = max(turn_ambiguity, float(abs(y - turns[j])))

    # Safety net against a noisy two-point bootstrap: re-round every turn
    # against the full-set slope until nothing moves. Converged immediately in
    # every case measured; exhausting the cap is a bug, not a soft outcome.
    for _ in range(_BETA_TURN_POLISH_MAX):
        slope, intercept = _fit_line(L, phi + two_pi * turns)
        polished = np.round((slope * L + intercept - phi) / two_pi)
        polished = polished - polished[0]
        if np.array_equal(polished, turns):
            break
        turns = polished
    else:
        raise ValueError(
            "beta_from_phases: wrap counts did not reach a fixed point in "
            f"{_BETA_TURN_POLISH_MAX} iterations at lengths {L.tolist()} — "
            "refusing to return an unconverged turn vector"
        )

    slope, intercept = _fit_line(L, phi + two_pi * turns)
    resid = float(np.max(np.abs((phi + two_pi * turns) - (slope * L + intercept))))

    # Neighbouring branches differ by 2π/ΔL in slope, so the seed is decisive
    # only within ±π/ΔL of the truth. Use the largest gap: it is the one that
    # wraps first.
    max_gap = float(np.max(np.diff(L))) if n > 1 else np.inf
    half_branch = np.pi / max_gap if np.isfinite(max_gap) and max_gap > 0 else np.inf
    margin = float(abs(slope - seed_beta) / half_branch) if np.isfinite(half_branch) else 0.0
    seed_excess = float(slope / seed_beta - 1.0) if seed_beta != 0.0 else float("nan")

    # Back to the caller's ordering before anyone compares these to `lengths`.
    unsort = np.argsort(order)
    return BetaFit(
        beta=float(slope), resid=resid, margin=margin, seed_excess=seed_excess,
        turns=turns[unsort], seed_turns=seed_turns[unsort],
        turns_shifted=not np.array_equal(turns, seed_turns),
        branch_ambiguity=branch_ambiguity, turn_ambiguity=turn_ambiguity,
    )


def beta_from_phases(lengths, phases, seed_beta) -> Tuple[float, float, float]:
    """``(beta, resid, margin)`` — the first three fields of `beta_fit_from_phases`.

    Kept as-is so the dozen existing 3-tuple unpack sites (tests, the mmi_1x2
    diagnostic script) do not have to change. New callers that want to know
    *why* a fit came out the way it did should call `beta_fit_from_phases` and
    read `seed_excess` / `turns_shifted` / `branch_ambiguity`.
    """
    fit = beta_fit_from_phases(lengths, phases, seed_beta)
    return fit.beta, fit.resid, fit.margin


def _resolved_spectral_config(wavelengths: np.ndarray, factor: float) -> Dict[str, Any]:
    """The *resolved* spectral quantities, for the cache fingerprint.

    Both are derived from `wavelengths`, which is already hashed, so the
    fingerprint is correct today without them. They are recorded anyway
    because a future change to how they are derived would otherwise be
    invisible to the cache — old entries would keep matching new calls that
    simulate something different. That is exactly how the monitor-grid bug
    survived as long as it did, and it is the pattern
    `build_run_config()`'s docstring prescribes for resolved-vs-ingredient
    values (Solver recommendation #4), already used for `resolved_layout`.

    Preferred over bumping `_CACHE_SCHEMA_VERSION`, which would also strand
    the 65 valid EME entries the bug never touched, and which relies on a
    human remembering to bump it.
    """
    return {
        "resolved_monitor_freqs": monitor_freq_grid(wavelengths)[0],
        "resolved_source_fwidth": source_fwidth(wavelengths, factor),
    }



def check_source_coverage(fcen: float, fwidth: float, wavelengths: np.ndarray,
                          min_fraction: float = 0.1) -> np.ndarray:
    """Warn if the pulse barely reaches some requested wavelength.

    A starved point still yields a well-conditioned transmission — the source
    spectrum cancels in the measured/normalization amplitude ratio — but it
    spends the simulation's dynamic range, and `decay_tol` is applied to DFT
    fields that are orders of magnitude smaller there. Cheap to check, and the
    failure it guards against is silent.

    Returns the per-wavelength amplitude as a fraction of peak.
    """
    import meep as mp

    src = mp.GaussianSource(fcen, fwidth=fwidth)
    amp = np.array([abs(src.fourier_transform(1.0 / float(w))) for w in wavelengths])
    peak = amp.max()
    if peak <= 0.0:
        logger.warning("source carries no power at any requested wavelength "
                       "(fcen=%.5f, fwidth=%.5f)", fcen, fwidth)
        return amp
    frac = amp / peak
    if frac.min() < min_fraction:
        worst = int(np.argmin(frac))
        logger.warning(
            "source is starved at %d/%d requested wavelengths (fcen=%.5f, "
            "fwidth=%.5f): worst is λ=%.4f µm at %.2g of peak amplitude, below "
            "the %.0f%% floor. Raise FDTDSolver(source_bandwidth_factor=...).",
            int((frac < min_fraction).sum()), len(frac), fcen, fwidth,
            float(wavelengths[worst]), frac[worst], 100 * min_fraction,
        )
    return frac


def _point_in_polygon(x: float, y: float, poly) -> bool:
    """Even-odd ray-casting test, boundary-exclusive. `poly`: [(x, y), ...]."""
    inside = False
    x1, y1 = poly[-1]
    for x2, y2 in poly:
        if (y1 > y) != (y2 > y):
            x_at_y = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < x_at_y:
                inside = not inside
        x1, y1 = x2, y2
    return inside


def check_reflection_layout(core_polygons, src_x: float, mon_x: float, y0: float,
                            extend_input_waveguide: bool, domain_x_edge: float) -> None:
    """Guard against the two real, silent bugs `solve_reflection()` has
    already shipped with once each (validation_markdown/
    high_contrast_grating_validation.md step 10 item 6): `source_margin`/
    `monitor_margin` are measured from the PML inner edge, not from the
    geometry's own edge, so a source/monitor position that looks reasonable
    in isolation can still land in bare oxide once `edge_clearance` (or any
    other domain-size parameter) changes; and a geometry that claims to
    honor `extend_input_waveguide` can silently ignore `stub_length` for its
    input side, leaving the cut facet the flag was meant to remove.

    Both checks are on the same `(2, 0)`-layer polygons `to_meep_prisms()`
    is about to extrude, so a pass here means the *actual* simulated
    silicon — not some other approximation of it — covers what the caller
    assumed it would. Raises before `mp.Simulation` is built, not after, so
    a bad layout fails in milliseconds instead of after a full solve.
    """
    for label, x in (("source", src_x), ("monitor", mon_x)):
        if not any(_point_in_polygon(x, y0, poly) for poly in core_polygons):
            raise ValueError(
                f"solve_reflection(): {label} at (x={x:.4f}, y={y0:.4f}) falls "
                f"outside the geometry's own silicon (layer (2, 0)) -- "
                f"source_margin/monitor_margin are measured from the PML inner "
                f"edge, not from the geometry's edge, so changing edge_clearance "
                f"(or sim_width/dpml) without re-checking them can silently move "
                f"the {label} into bare oxide. Check the geometry with a plot, "
                f"or adjust source_margin/monitor_margin, before re-running."
            )
    if extend_input_waveguide:
        probe_x = domain_x_edge + 1e-6
        if not any(_point_in_polygon(probe_x, y0, poly) for poly in core_polygons):
            raise ValueError(
                f"solve_reflection(): extend_input_waveguide=True but the "
                f"geometry's own silicon does not reach the domain boundary "
                f"(x={domain_x_edge:.4f}) at y={y0:.4f} after "
                f"geometry.polygons(stub_length=...) -- this geometry may be "
                f"ignoring stub_length for its input side, which defeats the "
                f"point of extend_input_waveguide: the cut facet the flag is "
                f"meant to remove would still be there, just moved deeper into "
                f"the PML instead of eliminated."
            )


@dataclass
class SectionFDTDLayout:
    """Source/monitor position + size for one solve_section() call, in the
    section's own local frame: x_local=0 at x0-stub_left, x_local=L at
    x1+stub_right — the same frame geometry.section_gds()/section_prisms()
    build their geometry in. Returned by section_fdtd_layout(), the single
    place this math lives; solve_section() shifts it into its own
    meep-centered frame (x_meep = x_local - L/2), and
    beef/plot/geometry.py's _draw_section_fdtd_markers() draws directly in
    this local frame — see section_fdtd_layout()'s docstring."""
    L: float
    n_excite: int
    n_measure: int
    reverse: bool  # True: excite from the right (n_out <= n_in)
    src_x: float
    norm_x: float
    mea_x: float
    src_c: float
    src_width: float
    mea_c: float
    mea_w: float
    kpoint: Tuple[float, float, float]
    alpha_idx: int  # 0 forward / 1 backward


@dataclass(frozen=True)
class PortSpec:
    """Explicit orientation + aperture for one solve_section() port,
    overriding section_fdtd_layout()'s geometry-derived default for EITHER
    the excited (`excite_port`) or measured (`measure_port`) side
    independently — see solve_section()'s own docstring for how the two
    compose. Unlike section_fdtd_layout()'s geometry-driven aperture
    (center_at/width_at/source_aperture_at/half_domain_aperture), a
    PortSpec always fully and explicitly specifies its own aperture: no
    DeviceGeometry today has an accessor for a y-normal port's local
    geometry (an x-in/y-out width_at() can't describe a fixed-y slice's
    x-extent) — see e.g. beef.geometry.ring_coupler.RingResonatorCouplerSection
    .ring_port(), which hands the caller exactly the plain floats a PortSpec
    needs rather than trying to satisfy width_at()'s existing contract.

    Attributes:
        axis: "x" (default — today's only orientation, a source/monitor
            plane at fixed x spanning y) or "y" (a plane at fixed y
            spanning x).
        sign: +1 (default) selects get_eigenmode_coefficients() band index 0
            (forward, +axis); -1 selects index 1 (backward, -axis). Which
            EDGE of `axis`'s local domain this port sits at depends on
            whether it's used as solve_section()'s excite_port or
            measure_port — see that method's docstring: an excite_port's
            edge is physically forced by its launch direction (LOW edge
            for +1, HIGH for -1), a measure_port's is the opposite
            convention (HIGH edge for +1, LOW for -1), since a monitor has
            no such constraint.
        center: device-frame coordinate TRANSVERSE to axis (y if axis="x",
            x if axis="y") of the source/monitor plane's own center.
        width: aperture span transverse to axis (µm). Must be > 0.
    """
    axis: Literal["x", "y"] = "x"
    sign: int = 1
    center: float = 0.0
    width: float = 0.0

    def __post_init__(self) -> None:
        if self.axis not in ("x", "y"):
            raise ValueError(f"PortSpec.axis must be 'x' or 'y', got {self.axis!r}.")
        if self.sign not in (1, -1):
            raise ValueError(f"PortSpec.sign must be +1 or -1, got {self.sign!r}.")
        if self.width <= 0.0:
            raise ValueError(f"PortSpec.width must be > 0, got {self.width!r}.")

    @property
    def alpha_idx(self) -> int:
        """get_eigenmode_coefficients() band index this port reads: 0
        (forward) for sign=+1, 1 (backward) for sign=-1 — verified against
        ring_resonator_hybrid_cascaded.py's get_coeffs()/get_norm_calibrated()
        (forward=True/False maps to alpha[...,0]/alpha[...,1] there too)."""
        return 0 if self.sign == 1 else 1


def _port_point(axis: str, along: float, transverse: float) -> Tuple[float, float]:
    """(x, y) for a source/monitor plane at position `along` its own axis,
    centered at `transverse` on the perpendicular axis. axis="x": along is
    x, transverse is y (today's only case). axis="y": along is y,
    transverse is x — verified against the reference scripts' own
    add_y_monitor()/y-normal EigenModeSource construction, mp.Vector3(xc, y)."""
    return (along, transverse) if axis == "x" else (transverse, along)


def _port_size(axis: str, width: float, mon_sz: float) -> Tuple[float, float, float]:
    """meep Vector3 size for a source/monitor plane: zero-thickness along
    its own normal, `width` transverse, `mon_sz` vertical (3-D core height
    or 0 for a 2-D sim)."""
    return (0.0, width, mon_sz) if axis == "x" else (width, 0.0, mon_sz)


def section_fdtd_layout(
    geometry,
    x0: float,
    x1: float,
    stub_left: float,
    stub_right: float,
    n_in: int,
    n_out: int,
    sim_width: float,
    source_margin: float,
    monitor_margin: float,
    half_domain_aperture: bool,
    geometry_kwargs: dict,
    monitor_center: Optional[float] = None,
    monitor_size: Optional[float] = None,
) -> SectionFDTDLayout:
    """Shared position/size math for one FDTDSolver.solve_section() call —
    which boundary excites vs. measures, the source/monitor y-centre and
    y-span on each side, and where they sit along x — factored out so
    solve_section() and the plotting code that visualizes it
    (_draw_section_fdtd_markers in beef/plot/geometry.py) can't drift
    apart; each previously reimplemented this independently, and an earlier
    bug-hunt found real mismatches between the two copies.

    Excites from whichever boundary has fewer modes (reciprocity makes the
    resulting transfer matrix identical either way — see solve_section()'s
    own docstring). The excited/normalization side's aperture follows a
    priority order: a geometry-provided single-mode launch window
    (geometry.source_aperture_at(), e.g. an MMI's closely-spaced output
    branches) wins over half_domain_aperture (span only sim_width/2, from
    the shared y=0 axis to the near PML edge, for a genuinely single-mode
    port with an unexcited neighbor nearby), which wins over the full
    sim_width. The measured/interface side is sized to sim_width whenever
    it must resolve >1 mode, else to the bare local waveguide width;
    monitor_center/monitor_size override that side's centre/span outright.

    Returns everything in the section's local frame (0..L, stub included,
    PML not) — callers needing a different frame (e.g. solve_section()'s
    own meep-centered one) shift the returned x positions themselves.
    """
    L = (x1 - x0) + stub_left + stub_right
    eps = min(1e-3, (x1 - x0) * 0.01)
    c_left = geometry.center_at(x0 + eps, **geometry_kwargs)
    c_right = geometry.center_at(x1 - eps, **geometry_kwargs)
    w_left = geometry.width_at(x0 + eps, **geometry_kwargs)
    w_right = geometry.width_at(x1 - eps, **geometry_kwargs)

    reverse = n_out <= n_in
    n_excite = n_out if reverse else n_in
    n_measure = n_in if reverse else n_out

    if reverse:
        src_x, norm_x, mea_x = L - source_margin, L - monitor_margin, 0.0
        src_c, mea_c, mea_w0 = c_right, c_left, w_left
        src_ap_x = x1 - eps
        kpoint = (-1.0, 0.0, 0.0)
        alpha_idx = 1
    else:
        src_x, norm_x, mea_x = source_margin, monitor_margin, L
        src_c, mea_c, mea_w0 = c_left, c_right, w_right
        src_ap_x = x0 + eps
        kpoint = (1.0, 0.0, 0.0)
        alpha_idx = 0

    mea_w = (sim_width if n_measure > 1 else mea_w0) if monitor_size is None else monitor_size
    if monitor_center is not None:
        mea_c = monitor_center
    elif n_measure > 1 and monitor_size is None:
        # A multi-mode interface's window is the FULL domain width, which is
        # only geometrically valid centered at y=0 -- center_at() may return
        # an arbitrary single PORT's position here, meaningless for "the
        # whole shared cross-section" this window is meant to represent
        # (e.g. a multi-port WaveguideRun with no explicit output_idx/
        # input_idx in geometry_kwargs silently defaults to port 0's
        # position), and combined with sim_width can push the window past
        # the PML on the far side. FOUND 2026-08-08 on the modal
        # multiplexer's new connector-adjacent interfaces (gap 7):
        # CompoundGeometry.center_at() returned +1.0375 ("top" branch, port
        # 0's default) for a query into the connector WaveguideRun with no
        # port index given, caught by _check_aperture() rather than
        # silently clipping (see that function's own docstring/history).
        mea_c = 0.0

    src_aperture = None
    if n_excite == 1:
        _aperture_at = getattr(geometry, "source_aperture_at", None)
        if _aperture_at is not None:
            src_aperture = _aperture_at(src_ap_x, **geometry_kwargs)
    if src_aperture is not None:
        src_width = min(src_aperture, sim_width - 2.0 * abs(src_c))
    elif n_excite == 1 and half_domain_aperture:
        src_width = sim_width / 2.0
        src_c = sim_width / 4.0 if src_c > 0 else (-sim_width / 4.0 if src_c < 0 else 0.0)
    else:
        src_width = sim_width
        if n_excite > 1:
            # Same reasoning/fix as mea_c above, mirrored for the excite
            # side's own multi-mode case. n_excite == 1 (a genuine lone
            # single-mode port) is deliberately left untouched here --
            # every existing device's own precedent for that fallback
            # keeps center_at()'s real position, and changing it would be
            # an unrelated, unverified behavior change for devices that
            # already work.
            src_c = 0.0

    return SectionFDTDLayout(
        L=L, n_excite=n_excite, n_measure=n_measure, reverse=reverse,
        src_x=src_x, norm_x=norm_x, mea_x=mea_x,
        src_c=src_c, src_width=src_width, mea_c=mea_c, mea_w=mea_w,
        kpoint=kpoint, alpha_idx=alpha_idx,
    )


def merged_mea_positions(layout: SectionFDTDLayout, lengths: Sequence[float]) -> List[float]:
    """Interior interface-monitor positions for a merged multi-length domain,
    in the same section-local frame `section_fdtd_layout()` returns.

    One monitor per requested length, all reading a single shared source in
    ONE domain sized to ``max(lengths)``:

        position[i] = (layout.L - lengths[i]) if layout.reverse else lengths[i]

    which reproduces each length's own single-domain source-to-monitor
    distance exactly -- derived algebraically in
    `resolved/section_mode_beta_multi_monitor_validation.md` step 2, not
    approximated. The longest length lands on `layout.mea_x` (the domain
    boundary); every shorter one sits interior, nearer the source.

    Factored out of `solve_section_multi_length()` for the same reason
    `section_fdtd_layout()` itself was: `beef/plot/calibration.py` draws
    these monitors and must not re-derive the formula to do it. A silently
    drifting second copy of exactly this kind of position math is the
    failure `section_fdtd_layout()`'s and `BetaCalibration`'s own docstrings
    both record.

    Args:
        layout: the merged domain's resolved layout -- i.e.
            `section_fdtd_layout()` called once with
            ``x0, x1 = mid -+ max(lengths)/2``.
        lengths: the propagation distances being measured (µm), in the
            order they are simulated/returned.

    Returns:
        One local-frame x position per entry of `lengths`, same order.
    """
    return [(layout.L - float(L_i)) if layout.reverse else float(L_i) for L_i in lengths]


class FDTDSolver:
    """MEEP-based FDTD solver.

    Device geometry is converted to meep Prism objects via
    DeviceGeometry.to_meep_prisms() — no GDS files are written.

    Args:
        res: mesh resolution (pixels/µm).
        sim_width: simulation y-span (µm).
        sim_height: simulation z-span (µm) when three_d=True — inner span,
            excluding PML, analogous to sim_width; the silicon core
            (material_stack.wg_height thick, centered on z=0) sits in oxide
            within it. In 2-D mode it is only a nominal slab reference and
            never enters the simulation.
        three_d: False (default) runs the legacy 2-D planar (x, y)
            simulation — waveguides extruded infinitely in z with bulk-index
            silicon, cheap but quantitatively a different structure (no
            vertical confinement; e.g. it under-predicts a directional
            coupler's coupling rate by ~2×). True runs full 3-D: the cell
            gains a z axis (sim_height + 2*dpml tall, PML on all six faces),
            the core is a real material_stack.wg_height-thick silicon layer
            in oxide, and every source/monitor spans sim_height in z, so
            their MPB eigenmode solves become true (y, z) cross-section
            solves — the same mode basis EMESolver's meow cross-sections
            use, which is what makes 3-D FDTD sections quantitatively
            consistent with EME sections in a hybrid chain.
        dpml: PML thickness on all four sides (µm).
        source_margin: x-offset of EigenModeSource from left PML inner edge (µm).
            Must satisfy source_margin < waveguide_length (see below).
        monitor_margin: x-offset of monitors from PML inner edges (µm).
            Input normalization monitor: left PML inner + monitor_margin.
            Output monitors: right PML inner - monitor_margin.
            Must satisfy monitor_margin < waveguide_length.
        num_input_modes: number of input modes captured by normalization monitor.
        num_output_modes: number of modes extracted per output waveguide monitor.
        output_monitor_centers: y-centers of the output waveguide monitors (µm).
            Required; must match the output waveguide positions of the device.
        output_monitor_widths: y-spans of the output monitors (µm).
            None (default) → as wide as each monitor can be while staying
            centered on its own output_monitor_centers entry and inside the
            non-PML domain: sim_width - 2*abs(center) — full sim_width for
            a centered (y=0) port, narrowing symmetrically for an off-axis
            one so it never reaches into the PML on its far side. Matches
            the clamp already used for the excited/source-side aperture
            (`min(src_aperture, sim_width - 2*abs(src_yo))`, this method's
            body) rather than an arbitrary fixed width — a monitor narrower
            than this default risks clipping real field content (and, for
            get_eigenmode_coefficients()'s own MPB solve, the surrounding
            cladding context that eigenmode solve needs) before the overlap
            integral runs. FIXED 2026-08-07 from a bare `[sim_width] *
            len(...)` default, which was wrong for any off-axis port (it
            would place a monitor of full sim_width off-center, reaching
            past the PML on the far side) — not previously hit because
            every caller to date passed an explicit override.
        half_domain_aperture: when True, solve_section()'s source+norm
            monitor on a side being excited with exactly 1 mode spans only
            *half* the domain — from the shared symmetry axis (y=0) out to
            the near PML edge on whichever side the excited arm sits (i.e.
            width sim_width/2, centered at +-sim_width/4) — instead of the
            default full sim_width. Matches the convention validated by the
            legacy (pre-BEEF) FDTD/EME hybrid reference implementation
            (BEEF_Legacy/masters/mmi_2x2_master.py's _port_simulation:
            aperture width = oxide_width/2 centered at +-oxide_width/4,
            i.e. this codebase's sim_width/sim_width equivalents) — a full
            sim_width aperture centered on the arm would also reach a
            second, physically-present-but-unexcited waveguide sitting
            nearby (e.g. a directional coupler's other arm, or an MMI's
            other output branch), corrupting the "single-mode" eigenmode
            solve into seeing the near-degenerate even/odd pair of both
            waveguides instead of isolating the intended one. False
            (default) keeps the full sim_width — correct for a lone port
            with no nearby companion (e.g. a single-input MMI's input side).
            Ignored (sim_width used) whenever more than 1 mode is being
            excited at that boundary, where the full aperture is still
            needed/correct. Also superseded whenever the geometry itself
            provides a launch window (source_aperture_at() returns non-None
            at the excited cross-section — see DeviceGeometry.
            source_aperture_at): that geometry-specific window wins, since
            it can isolate ports this fixed half-domain rule can't (e.g. a
            port sitting exactly on the y=0 axis). Not used by solve() — that method's source is
            always the one true input port, with no neighboring port nearby
            to avoid, so it always spans the full sim_width regardless of
            this setting.
        material_stack: material index lookup for Si / SiO₂.
        dispersive: build the media from the re-parameterised single-Lorentzian
            fit over this run's band instead of freezing one index at band
            centre. `None` (default) inherits `material_stack.dispersive`;
            True/False override it for this solver only, without mutating the
            stack — the `SOI_160`/`SOI_220` singletons are shared process-wide.

            A broadband run gets one Medium for its whole span, so the frozen
            index is wrong at the band edges by up to 0.008 in n_Si on the
            O-band, while the EME side re-solves per wavelength. Turning this
            on cuts that mismatch ~16×.

            **On by default as of 2026-08-01**, via the stack: `SOI_160` and
            `SOI_220` now carry `dispersive=True`, so `None` inherits True. The
            strip-test gate that gated it has passed on all four stack/band
            combinations. Pass `False` here to freeze this solver alone without
            mutating the shared stack singletons.
            Single-wavelength calls (mode solves) are unaffected either way —
            a scalar wavelength has no band to fit, so the frozen index at that
            exact wavelength is already the correct medium.
        decay_tol: DFT field decay tolerance for early-stop criterion.
            Lower = more accurate but slower; 1e-9 is faster, 1e-11 is tight.
        minimum_run_time: floor (meep time units) before
            stop_when_dft_decayed()'s check is even evaluated. meep's own
            default is 0 -- no floor at all -- so a loose-enough decay_tol
            can be satisfied by a fast, near-source transient before light
            has physically had time to reach a distant monitor, silently
            reporting near-zero transmission for a reason that has nothing
            to do with the device's real coupling. FOUND 2026-08-08 on the
            modal multiplexer (66 um device): decay_tol=1e-9 with no floor
            stopped the sim at t=107-112 -- not even one transit's worth of
            time for a device this long -- while decay_tol=1e-11 (no floor
            either, just tight enough that the trivial-early-stop can't
            satisfy it) ran to t~77,175 and STILL returned near-zero
            transmission, an unrelated, still-open problem (see
            modal_multiplexer_validation.md gap 3). Default 0 (meep's own
            default, unchanged) so no other device's already-recorded
            behavior/golden moves; every caller to date has relied on this
            default implicitly, same reasoning as output_monitor_widths'
            own history. Set explicitly, per device, to comfortably exceed
            one full pulse transit of the actual structure (device length /
            group velocity, in meep's c=1 units) before trusting decay_tol
            alone to know when to stop.
        source_bandwidth_factor: Gaussian pulse width as a multiple of the
            requested frequency span 1/λ_min - 1/λ_max (see source_fwidth()).
            Must exceed 1 or the band edges are starved — the default 2.0
            leaves them at 0.27 of peak amplitude, and every solve warns
            (check_source_coverage()) if any requested wavelength drops below
            10%. For a single-wavelength request, where there is no span to
            scale, it multiplies a 10%-of-fcen relative bandwidth instead.
        eps_averaging, subpixel_tol, subpixel_maxeval: meep's subpixel
            smoothing controls, passed straight through to mp.Simulation. The
            defaults are meep's own, so behaviour is unchanged; they are
            constructor parameters rather than implicit defaults because they
            move n_eff at the second decimal (notes.md Solver #9) and anything
            that moves the physical result has to be visible to the cache
            fingerprint. Turning eps_averaging off is the usual first check
            when a result looks resolution-sensitive.

    Derived attribute:
        waveguide_length = max(source_margin, monitor_margin) + 0.25 (µm)
            Access waveguide stub between PML and device edge. Auto-sized to
            contain source and monitors with a 250 nm safety margin.
    """

    def __init__(
        self,
        res: float = 20.0,
        sim_width: float = 5.0,
        sim_height: float = 1.5,
        dpml: float = 0.5,
        source_margin: float = 0.0,
        monitor_margin: float = 0.25,
        num_input_modes: int = 1,
        num_output_modes: int = 1,
        output_monitor_centers: Optional[List[float]] = None,
        output_monitor_widths: Optional[List[float]] = None,
        half_domain_aperture: bool = False,
        three_d: bool = False,
        material_stack: MaterialStack = SOI_160,
        dispersive: Optional[bool] = None,
        decay_tol: float = 1e-11,
        minimum_run_time: float = 0.0,
        source_bandwidth_factor: float = 2.0,
        eps_averaging: bool = True,
        subpixel_tol: float = 1e-4,
        subpixel_maxeval: int = 100000,
    ):
        self.res = res
        self.sim_width = sim_width
        self.sim_height = sim_height
        self.three_d = three_d
        self.dpml = dpml
        self.source_margin = source_margin
        self.monitor_margin = monitor_margin
        self.num_input_modes = num_input_modes
        self.num_output_modes = num_output_modes
        self.output_monitor_centers = list(output_monitor_centers or [])
        self.output_monitor_widths = output_monitor_widths  # None → sim_width per monitor
        self.half_domain_aperture = half_domain_aperture     # True → half-domain aperture for single-mode ports
        self.material_stack = material_stack
        self.dispersive = dispersive
        self.decay_tol = decay_tol
        self.minimum_run_time = minimum_run_time
        self.source_bandwidth_factor = source_bandwidth_factor
        self.eps_averaging = eps_averaging
        self.subpixel_tol = subpixel_tol
        self.subpixel_maxeval = subpixel_maxeval
        # Derived: stub long enough that source + monitors sit in physical waveguide
        self.waveguide_length = max(source_margin, monitor_margin) + 0.25

    # 2-D/3-D switches — the three quantities that differ between the modes,
    # so solve()/solve_section() share one code path (a z extent of 0 makes
    # meep run the cell/source/monitor as 2-D planar objects).
    @property
    def _cell_sz(self) -> float:
        """Total z cell span incl. PML (0 → meep runs 2-D)."""
        return self.sim_height + 2 * self.dpml if self.three_d else 0.0

    @property
    def _mon_sz(self) -> float:
        """z-span of every source/monitor (inner sim only, not PML; 0 in 2-D)."""
        return self.sim_height if self.three_d else 0.0

    @property
    def _wg_height(self) -> Optional[float]:
        """Core extrusion thickness handed to the geometry's prism builders
        (None → legacy infinite 2-D extrusion)."""
        return self.material_stack.wg_height if self.three_d else None

    def _media(self, wavelengths):
        """(si_medium, sio2_medium) for this run's band.

        One place so `_build_domain()` and `solve_section()` cannot drift.
        `self.dispersive` overrides `self.material_stack.dispersive` when set;
        None inherits the stack, which is what every existing caller gets.

        A broadband FDTD run gets ONE Medium for its whole span, so a frozen
        centre index is wrong at the band edges — up to 0.008 in n_Si on the
        O-band, against an EME side that re-solves per wavelength. Passing the
        band lets the dispersive path fix that; the non-dispersive path ignores
        it and returns the same frozen medium as before.
        """
        dispersive = self.material_stack.resolve_dispersive(self.dispersive)
        wl = wavelengths if dispersive else float(np.mean(wavelengths))
        return (self.material_stack.si_medium(wl, dispersive=dispersive),
                self.material_stack.sio2_medium(wl, dispersive=dispersive))

    def _build_domain(self, geometry: DeviceGeometry, wavelengths: np.ndarray):
        """Return (Sx, Sy, wl_center, si_medium, sio2_medium).

        sim_width/sim_height are inner dimensions (excluding PML).
        Sx and Sy include 2*dpml on each axis (as does the z axis in 3-D —
        see _cell_sz).
        """
        wl_center = float(np.mean(wavelengths))
        Sx = geometry.total_length + 2 * self.waveguide_length + 2 * self.dpml
        Sy = self.sim_width + 2 * self.dpml   # full cell: inner + PML on both sides
        si_med, sio2_med = self._media(wavelengths)
        return Sx, Sy, wl_center, si_med, sio2_med

    def solve(
        self,
        geometry: DeviceGeometry,
        wavelengths: np.ndarray,
        port: str = "top",
        input_mode: int = 1,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        input_idx: Optional[int] = None,
        **kwargs,
    ) -> SMatrix:
        """Run a broadband FDTD simulation for one input port.

        A single GaussianSource pulse excites ``input_mode``; DFT mode monitors
        capture complex amplitudes at all wavelengths in one run.

        Args:
            geometry: device geometry (implements to_meep_prisms()).
            wavelengths: wavelengths to evaluate (µm), determines DFT frequencies.
            port: input port label — "top" selects geometry.input_ys' last
                  (highest-y) entry, "bot" its first (lowest-y) entry; both
                  collapse to the same position for a single-input device
                  (1×2 MMI), since input_ys has only one entry there. Ignored
                  when ``input_idx`` is given.
            input_mode: 1-based waveguide mode index to excite and normalize against.
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache — see
                solve_section()'s save_dir docstring for the full semantics
                (identical here, fingerprinted via
                geometry.polygons(stub_length=...) instead of section_gds()).
            cache_only: if True, never run the simulation — return the
                cached result or raise CacheMissError. See solve_section()'s
                cache_only docstring for the full semantics (identical here).
            input_idx: explicit 0-based index into ``geometry.input_ys`` to
                excite, overriding ``port``. Added for devices with more than
                two access ports (e.g. a 4x4 MMI, `mmi_4x4_90hybrid_fdtd.py`),
                where "top"/"bot" can only ever reach the first/last entry —
                an interior port (e.g. input_ys[2] of 4) was previously
                unreachable through this method at all. None (default)
                preserves the exact original port="top"/"bot" behaviour.
                Also forwarded to geometry.source_aperture_at() (replacing
                whatever the caller might separately pass via **kwargs under
                that name — passing both raises immediately, so the aperture
                width is never computed for a different port than the one
                actually excited).
            **kwargs: forwarded to geometry.source_aperture_at() (e.g.
                  input_idx for DirectionalCoupler) when the geometry
                  provides one — see the source-width note in the body.

        Returns:
            SMatrix with ports ``"out_0"``, ``"out_1"``, … matching the order of
            output_monitor_centers.  Entry (out_i, input_mode, mode_n) holds the
            complex amplitude coefficient, normalised so |coeff|² = transmission.
        """
        import meep as mp

        if not self.output_monitor_centers:
            raise ValueError(
                "FDTDSolver.output_monitor_centers is empty — "
                "set it to the y-positions of the output waveguides (µm)."
            )

        wavelengths = np.asarray(wavelengths)

        config = None
        if input_idx is not None and "input_idx" in kwargs:
            raise TypeError(
                "solve() got 'input_idx' both as a named argument and inside "
                "**kwargs -- pass it once; the named argument selects the "
                "excited port AND is what gets forwarded to "
                "source_aperture_at(), so a second, different value in "
                "kwargs would silently size the aperture for the wrong port."
            )

        if save_dir:
            stub_length = self.waveguide_length + self.dpml
            config = build_whole_device_config(
                self, "solve",
                dict(port=port, input_idx=input_idx, input_mode=input_mode,
                     wavelengths=wavelengths, kwargs=kwargs,
                     **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor),
                     **resolved_material_config(
                         self.material_stack, wavelengths,
                         self.material_stack.resolve_dispersive(self.dispersive))),
                geometry, geometry.polygons(stub_length=stub_length),
            )
            hit = lookup(save_dir, config)
            if hit is not None:
                logger.info("solve: cache hit at %s", hit)
                return load_result(hit)
        if cache_only:
            raise CacheMissError(
                f"solve(): cache_only=True but no cached run found for this "
                f"config (save_dir={save_dir!r}). Run once with cache_only=False "
                f"to populate the cache."
            )

        Sx, Sy, wl_center, si_med, sio2_med = self._build_domain(geometry, wavelengths)
        fcen = 1.0 / wl_center
        fwidth = source_fwidth(wavelengths, self.source_bandwidth_factor)
        mon_freqs, wl_order = monitor_freq_grid(wavelengths)
        check_source_coverage(fcen, fwidth, wavelengths)

        geometry_meep = geometry.to_meep_prisms(
            si_med, stub_length=self.waveguide_length + self.dpml,
            wg_height=self._wg_height,
        )

        # Source position: centered on the actual input waveguide via
        # geometry.input_ys, not a hand-rolled +-waveguide_offset/2 offset —
        # that formula ignored the geometry's real port positions (e.g. it
        # placed the source off-axis for a single-input device, whose one
        # true input sits at input_ys[0]=0.0, not +-waveguide_offset/2).
        input_ys = getattr(geometry, "input_ys", [0.0])
        # input_idx param (if given) overrides the port="top"/"bot" string --
        # that string can only ever reach input_ys' first/last entry, which
        # is a real, previously-unfixable gap for any device with more than
        # two access ports (an interior port, e.g. input_ys[2] of 4, was
        # simply unreachable through this method before). See this method's
        # own docstring for the full reasoning.
        resolved_input_idx = input_idx if input_idx is not None else (
            len(input_ys) - 1 if port == "top" else 0
        )
        src_yo = input_ys[resolved_input_idx]
        src_x = -Sx / 2 + self.dpml + self.source_margin
        # Default: the full non-PML domain — PML edge to PML edge, exclusive
        # of PML — correct for the usual case of one lone input waveguide
        # with no neighboring port nearby (so half_domain_aperture, which is
        # for narrow *output* ports, never applies here). But when the
        # geometry itself reports a launch window at the source cross-
        # section (source_aperture_at() non-None — e.g. DirectionalCoupler,
        # whose unexcited second arm is physically present on the input side
        # too), the source + norm monitor shrink to that window, centered on
        # the input arm and clamped so an off-axis window never reaches into
        # the PML — same priority rule as solve_section()'s.
        # src_x is in the device frame too: the meep cell is centered on the
        # device (Sx = total_length + 2*(waveguide_length + dpml)), so
        # meep x == device x here.
        src_width = self.sim_width
        _aperture_at = getattr(geometry, "source_aperture_at", None)
        if _aperture_at is not None:
            # input_idx explicitly passed through (not left to whatever's in
            # kwargs) so the aperture is always sized for the port actually
            # excited above -- the mismatch the TypeError guard earlier in
            # this method exists to prevent.
            src_aperture = _aperture_at(src_x, input_idx=resolved_input_idx, **kwargs)
            if src_aperture is not None:
                src_width = min(src_aperture, self.sim_width - 2.0 * abs(src_yo))

        src = mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=fwidth),
            center=mp.Vector3(src_x, src_yo),
            size=mp.Vector3(0, src_width, self._mon_sz),   # span inner sim only, not PML
            eig_band=input_mode,
            eig_kpoint=mp.Vector3(1, 0, 0),
            eig_match_freq=True,
        )

        in_mon_x = -Sx / 2 + self.dpml + self.monitor_margin
        out_mon_x = Sx / 2 - self.dpml - self.monitor_margin

        sim = mp.Simulation(
            resolution=self.res,
            cell_size=mp.Vector3(Sx, Sy, self._cell_sz),
            boundary_layers=[mp.PML(self.dpml)],
            geometry=geometry_meep,
            sources=[src],
            default_material=sio2_med,
            eps_averaging=self.eps_averaging,
            subpixel_tol=self.subpixel_tol,
            subpixel_maxeval=self.subpixel_maxeval,
        )

        # Input normalization monitor — same aperture as the source
        in_flux = sim.add_mode_monitor(
            mon_freqs,
            mp.FluxRegion(center=mp.Vector3(in_mon_x, src_yo),
                          size=mp.Vector3(0, src_width, self._mon_sz)),
        )

        # One output monitor per output waveguide. Default width clamps to
        # TWO limits, not just one: the PML edge (sim_width - 2*abs(yc), as
        # before) AND the gap to the nearest neighboring output port. The
        # PML-only clamp is correct for a single, isolated multi-mode
        # monitor (e.g. an MMI body cross-section spanning every waveguide
        # together on purpose), but for a multi-PORT solve() call with
        # several individually-narrow output_monitor_centers close together
        # (this device's 3 output branches at WO=1.0375 um pitch,
        # sim_width=5.0), the PML-only clamp let the on-axis "mid" monitor
        # reach the full 5.0 um sim_width -- completely engulfing the other
        # two branches at +/-1.0375 um. get_eigenmode_coefficients() solves
        # eigenmodes of whatever dielectric sits inside its window, so that
        # monitor was decomposing onto the fundamental supermode of all
        # three coupled waveguides, not the isolated mid-branch mode --
        # exactly the "double-counted output branch" check_passive() itself
        # warns about. Traced directly to this via modal_multiplexer's own
        # gap-3 reference (TE1 totals 1.1021 at 1.26 um, persisting even
        # after the input-side normalization fix above, which correctly
        # never fired since the corruption is per-output-port, not at the
        # shared source). Clamping to the nearest-neighbor gap makes
        # adjacent monitors meet exactly at their shared midpoint with zero
        # overlap for evenly-spaced ports (this device's case); for a lone
        # output port (no neighbors), or one far enough from its neighbors
        # that the PML edge is already the tighter limit, this clamp is a
        # no-op and prior behavior is unchanged.
        def _neighbor_gap(i: int) -> float:
            yc = self.output_monitor_centers[i]
            others = [abs(yc - oc) for j, oc in enumerate(self.output_monitor_centers) if j != i]
            return min(others) if others else self.sim_width

        out_widths = (
            self.output_monitor_widths
            if self.output_monitor_widths is not None
            else [
                min(self.sim_width - 2.0 * abs(yc), _neighbor_gap(i))
                for i, yc in enumerate(self.output_monitor_centers)
            ]
        )
        out_fluxes = [
            sim.add_mode_monitor(
                mon_freqs,
                mp.FluxRegion(
                    center=mp.Vector3(out_mon_x, yc),
                    size=mp.Vector3(0, yw, self._mon_sz),
                ),
            )
            for yc, yw in zip(self.output_monitor_centers, out_widths)
        ]

        sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol, minimum_run_time=self.minimum_run_time))

        # Normalise by forward input-mode amplitude so |coeff|² = physical T.
        # Same mismatch guard as solve_section()'s x-axis excite branch (see
        # _EXCITE_NORM_MISMATCH_TOL's docstring): get_eigenmode_coefficients()
        # alone can silently undercount launched power for an ill-conditioned
        # decomposition, so cross-check against mp.get_fluxes() (a direct
        # field integral) and only fall back to flux-calibrated normalization
        # when they disagree. Found via modal_multiplexer's own gap-3
        # reference: the ONE check_passive() violation in that "closed" run
        # (TE1 totals 1.1021 at 1.26 um) traces to this exact unconditional
        # alpha here, this device's own input waveguide being the excited
        # geometry rather than a connector this time — the failure mode isn't
        # specific to near-degenerate multimode bases the way gap 7's fix
        # was motivated, just to get_eigenmode_coefficients() itself.
        in_coeff = sim.get_eigenmode_coefficients(
            in_flux, [input_mode], eig_parity=mp.NO_PARITY
        )
        in_alpha_raw = in_coeff.alpha[0, :, 0]  # [nfreq] forward amplitude, monitor order
        in_flux_mag = np.abs(np.asarray(mp.get_fluxes(in_flux)))
        in_mismatch = np.abs(np.sqrt(in_flux_mag) - np.abs(in_alpha_raw)) / np.maximum(np.sqrt(in_flux_mag), 1e-30)
        if np.any(in_mismatch > _EXCITE_NORM_MISMATCH_TOL):
            logger.warning(
                "solve: input mode %d normalization mismatch between "
                "get_eigenmode_coefficients() and get_fluxes() (max %.1f%% at freq "
                "index %d) — falling back to flux-calibrated normalization.",
                input_mode, 100 * float(in_mismatch.max()), int(in_mismatch.argmax()),
            )
            src_alpha = np.sqrt(in_flux_mag) * np.exp(1j * np.angle(in_alpha_raw))
        else:
            src_alpha = in_alpha_raw

        smat = SMatrix(wavelengths=wavelengths)
        for i, flux_mon in enumerate(out_fluxes):
            out_port = f"out_{i}"
            for mode_n in range(1, self.num_output_modes + 1):
                coeff = sim.get_eigenmode_coefficients(
                    flux_mon, [mode_n], eig_parity=mp.NO_PARITY
                )
                # [wl_order] lifts the ratio out of meep's ascending-frequency
                # monitor order back onto `wavelengths` as the caller gave it.
                alpha = (coeff.alpha[0, :, 0] / src_alpha)[wl_order]
                smat.set(alpha, port=out_port, in_mode=input_mode, out_mode=mode_n)

        sim.reset_meep()

        if save_dir:
            folder = save_run(save_dir, config, smat)
            logger.info("solve: cache saved at %s", folder)

        return smat

    def solve_reflection(
        self,
        geometry: DeviceGeometry,
        wavelengths: np.ndarray,
        port: str = "top",
        input_mode: int = 1,
        edge_clearance: float = 0.0,
        extend_input_waveguide: bool = False,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        **kwargs,
    ) -> SMatrix:
        """Broadband FDTD reflection spectrum for a single-port device.

        solve() cannot report anything physical for a device with no far-side
        output port to normalize against — a mirror, a grating coupler, or
        any other structure whose only well-defined waveguide cross-section
        is the one it is excited from (e.g. a GdsGeometry import with one
        access waveguide feeding a freeform region and no matching waveguide
        on the opposite side). It requires output_monitor_centers on the
        opposite domain edge and, at its own normalization monitor, only
        ever reads the forward-traveling mode coefficient
        (``coeff.alpha[0, :, 0]``) — never the backward one.

        This method reads BOTH directions off that same monitor instead:
        ``alpha[0, :, 0]`` (forward, the incident/normalization amplitude)
        and ``alpha[0, :, 1]`` (backward, the amplitude reflected back into
        the excited mode) — meep's standard single-run S11 decomposition,
        and the same forward/backward convention ``PortSpec.alpha_idx``
        documents elsewhere in this module (validated there against
        ``ring_resonator_hybrid_cascaded.py``'s ``get_coeffs()``). No second
        "incident-only" reference run is needed: both amplitudes come off
        one monitor's DFT fields in the single structured run.

        Domain sizing deliberately does NOT reuse solve()'s
        ``Sx = total_length + 2*waveguide_length + 2*dpml``. That formula
        assumes ``geometry.polygons(stub_length=...)`` extends the drawn
        silicon by ``waveguide_length`` on each side to reach the PML —
        true for every parametric geometry (MMI, Taper, ...), false for
        GdsGeometry and any geometry whose docstring says stub_length is
        ignored. Relying on it here would place the source/monitor in bare
        oxide past the geometry's own silicon, not inside a waveguide.
        Instead PML sits ``edge_clearance`` (plain oxide, default 0 — flush)
        outside ``geometry.total_length``'s own bounding box on both sides,
        and ``source_margin``/``monitor_margin`` (measured inward from the
        PML inner edge, same meaning as solve()) must be picked small enough
        to land inside whatever real input waveguide the geometry itself
        drew there — check with a geometry plot before trusting a run; there
        is no automatic guard for this, because the geometry does not
        generally expose its own stub length.

        A geometry whose input waveguide simply ENDS at
        ``geometry.total_length``'s own edge — the default, and the only
        option before 2026-08-04 — leaves a clean Si/oxide facet sitting in
        the ``edge_clearance`` gap (or right at the PML if
        edge_clearance=0). That facet reflects part of the guided mode on
        its own, independent of anything downstream — contaminating this
        exact reflection measurement with an artifact that has nothing to
        do with the device being characterized. ``extend_input_waveguide=
        True`` fixes this by requesting ``geometry.polygons(stub_length=
        edge_clearance + self.dpml)`` — enough to run the input waveguide
        fully through the PML to the true simulation-cell boundary, where
        PML absorbs it as an ordinary waveguide mode instead of reflecting
        off a cut end. Requires the geometry to actually honor
        ``stub_length`` for its input side (e.g. ``GdsGeometry`` with
        ``input_width`` set at construction — see its own ``polygons()``
        docstring for why that extension is one-sided, not symmetric);
        raises whatever that geometry itself raises if it cannot.

        Args:
            geometry: device geometry (implements to_meep_prisms()).
            wavelengths: wavelengths to evaluate (µm), determines DFT frequencies.
            port: input port label — as solve(): "top" selects
                geometry.input_ys' last (highest-y) entry, "bot" its first.
            input_mode: 1-based waveguide mode index to excite and measure.
            edge_clearance: plain-oxide gap (µm) between geometry.total_length's
                bounding box and the PML, both sides. 0 (default) puts PML
                flush against the geometry's own edges.
            extend_input_waveguide: see above. False (default) preserves the
                pre-2026-08-04 behavior (geometry.polygons(stub_length=0)).
            save_dir, cache_only: as solve().
            **kwargs: forwarded to geometry.source_aperture_at(), as solve().

        Returns:
            SMatrix with a single entry, ``port="in"``,
            ``in_mode=out_mode=input_mode``: complex amplitude whose
            ``|.|**2`` is the power reflected back into that mode, i.e.
            ``smat.transmission(port="in", in_mode=input_mode,
            out_mode=input_mode)`` gives R(λ). ``smat.check_passive()`` is a
            necessary-but-not-sufficient sanity check that R never exceeds 1
            — it says nothing about where the rest of the power went
            (radiated, scattered, or genuinely absorbed).
        """
        import meep as mp

        wavelengths = np.asarray(wavelengths)
        # 0.0 (default): pre-2026-08-04 behavior, geometry draws exactly what
        # it imported/parametrized, nothing more. edge_clearance + self.dpml
        # (extend_input_waveguide=True): run the input waveguide fully
        # through the PML to the true cell boundary — see docstring above.
        stub_length = (edge_clearance + self.dpml) if extend_input_waveguide else 0.0

        config = None
        if save_dir:
            config = build_whole_device_config(
                self, "solve_reflection",
                dict(port=port, input_mode=input_mode, edge_clearance=edge_clearance,
                     extend_input_waveguide=extend_input_waveguide,
                     wavelengths=wavelengths, kwargs=kwargs,
                     **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor),
                     **resolved_material_config(
                         self.material_stack, wavelengths,
                         self.material_stack.resolve_dispersive(self.dispersive))),
                geometry, geometry.polygons(stub_length=stub_length),
            )
            hit = lookup(save_dir, config)
            if hit is not None:
                logger.info("solve_reflection: cache hit at %s", hit)
                return load_result(hit)
        if cache_only:
            raise CacheMissError(
                f"solve_reflection(): cache_only=True but no cached run found "
                f"for this config (save_dir={save_dir!r}). Run once with "
                f"cache_only=False to populate the cache."
            )

        wl_center = float(np.mean(wavelengths))
        Sx = geometry.total_length + 2 * (edge_clearance + self.dpml)
        Sy = self.sim_width + 2 * self.dpml
        si_med, sio2_med = self._media(wavelengths)
        fcen = 1.0 / wl_center
        fwidth = source_fwidth(wavelengths, self.source_bandwidth_factor)
        mon_freqs, wl_order = monitor_freq_grid(wavelengths)
        check_source_coverage(fcen, fwidth, wavelengths)

        geometry_meep = geometry.to_meep_prisms(
            si_med, stub_length=stub_length, wg_height=self._wg_height,
        )

        input_ys = getattr(geometry, "input_ys", [0.0])
        input_idx = len(input_ys) - 1 if port == "top" else 0
        src_yo = input_ys[input_idx]
        src_x = -Sx / 2 + self.dpml + self.source_margin
        src_width = self.sim_width
        _aperture_at = getattr(geometry, "source_aperture_at", None)
        if _aperture_at is not None:
            src_aperture = _aperture_at(src_x, **kwargs)
            if src_aperture is not None:
                src_width = min(src_aperture, self.sim_width - 2.0 * abs(src_yo))

        src = mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=fwidth),
            center=mp.Vector3(src_x, src_yo),
            size=mp.Vector3(0, src_width, self._mon_sz),
            eig_band=input_mode,
            eig_kpoint=mp.Vector3(1, 0, 0),
            eig_match_freq=True,
        )

        mon_x = -Sx / 2 + self.dpml + self.monitor_margin

        core_polygons = geometry.polygons(stub_length=stub_length).get((2, 0), [])
        check_reflection_layout(
            core_polygons, src_x, mon_x, src_yo,
            extend_input_waveguide=extend_input_waveguide, domain_x_edge=-Sx / 2,
        )

        sim = mp.Simulation(
            resolution=self.res,
            cell_size=mp.Vector3(Sx, Sy, self._cell_sz),
            boundary_layers=[mp.PML(self.dpml)],
            geometry=geometry_meep,
            sources=[src],
            default_material=sio2_med,
            eps_averaging=self.eps_averaging,
            subpixel_tol=self.subpixel_tol,
            subpixel_maxeval=self.subpixel_maxeval,
        )

        mon = sim.add_mode_monitor(
            mon_freqs,
            mp.FluxRegion(center=mp.Vector3(mon_x, src_yo),
                          size=mp.Vector3(0, src_width, self._mon_sz)),
        )

        sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol, minimum_run_time=self.minimum_run_time))

        coeff = sim.get_eigenmode_coefficients(mon, [input_mode], eig_parity=mp.NO_PARITY)
        fwd = coeff.alpha[0, :, 0]   # [nfreq] incident amplitude, monitor order
        bwd = coeff.alpha[0, :, 1]   # [nfreq] reflected amplitude, monitor order

        smat = SMatrix(wavelengths=wavelengths)
        r = (bwd / fwd)[wl_order]
        smat.set(r, port="in", in_mode=input_mode, out_mode=input_mode)

        sim.reset_meep()

        if save_dir:
            folder = save_run(save_dir, config, smat)
            logger.info("solve_reflection: cache saved at %s", folder)

        return smat

    # ------------------------------------------------------------------
    # Mode inspection (no time-stepping)
    # ------------------------------------------------------------------

    def section_mode_neff(
        self,
        geometry,
        wavelengths: np.ndarray,
        x0: float,
        x1: float,
        n_modes: int,
        monitor_center: Optional[float] = None,
        monitor_size: Optional[float] = None,
        pad: Optional[float] = None,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        **geometry_kwargs,
    ) -> np.ndarray:
        """Effective indices of the first `n_modes` modes at this section's
        cross-section, from meep's own eigenmode solver.

        This is the FDTD-side counterpart of ``EMESolver.compute_modes()``, and
        exists so the two engines' mode bases can be compared directly.
        ``HybridSolver`` contracts index *i* of an FDTD interface against index
        *i* of an EME body assuming they are the same mode, but the two solve
        MPB on *different grids* — meep's own Yee mesh here, meow's ``Mesh2D``
        there — so their n_eff differ slightly and near-degenerate higher-order
        modes can come back in a different order. If they do, the chain
        silently mixes modes and nothing else in the codebase would notice.
        See dc_validation.md item 12.

        No time-stepping: ``get_eigenmode()`` runs MPB on the initialized grid
        and returns, so this costs a ``set_epsilon`` plus one eigensolve per
        mode per wavelength.

        **Cheap only relative to what it used to sit next to.** "Seconds, not
        minutes" was this docstring's claim while a `solve_section()` run took
        ~44 min; measured on the MMI 2x2 body it is 14 modes x 11 wavelengths =
        154 eigensolves at ~12.3 s = **31.5 min** (job 11937894). That was 1.9%
        of the pre-2026-08-17 calibration and is ~39% of one now — nothing
        about it changed, everything around it got ~10x cheaper. Hence
        `save_dir`: the result is a pure function of the cross-section,
        wavelengths, mode count and monitor plane, and re-running a calibrating
        example otherwise pays those 31.5 min again for an answer that cannot
        have moved.

        n_eff is ``k/f``: meep reports the Bloch wavevector in units of 2π/a
        and the frequency in c/a, so their ratio is the index directly, with no
        2π anywhere. Pinned against a bare 450×220 nm Si strip, where the
        answer must fall between the cladding and core indices — the competing
        normalizations land at 0.37 and 14.5 and are not close calls.

        **One structure, every wavelength.** All ``len(wavelengths)`` solves run
        on a single ``mp.Simulation``, because the geometry does not depend on
        wavelength and rebuilding it per point would multiply the cost by
        ``len(wavelengths)`` for nothing. The *material* does depend on
        wavelength, though, so under ``dispersive=True`` the media are built
        across the band (one meep medium whose ε is frequency-dependent) rather
        than frozen at its centre — fixed 2026-07-31, see
        ``_build_section_mode_sim()``. Under ``dispersive=False`` (today's
        default) the index is still frozen at ``mean(wavelengths)``, so results
        away from band centre carry that error: up to 0.008 in n_Si across the
        O-band (``material_model_validation.md``).

        Args:
            x0, x1: device-frame bounds of the section (µm). The cross-section
                is taken at their midpoint, so both must lie in the same
                uniform stretch of geometry for the result to mean anything.
            n_modes: number of bands to solve, ordered as meep returns them.
            monitor_center / monitor_size: y-centre and y-span of the plane the
                modes are solved on, same meaning as in `solve_section()`.
                Defaults to the full non-PML width on the ``y=0`` axis, which
                is the interface convention a coupled pair needs — the
                geometry-driven default would inherit one arm's centreline and
                make "even/odd" meaningless.
            pad: see `_build_section_mode_sim()`. None (default) is `self.dpml`,
                unchanged from before this parameter existed.
            save_dir / cache_only: as `solve_section()`. A hit returns the
                stored n_eff without building the simulation at all. Mode-count
                superset reuse applies: a cached 14-mode solve serves a later
                request for fewer modes by slicing, since band *i* does not
                depend on how many bands were asked for.

        Returns:
            Real ndarray ``[len(wavelengths), n_modes]`` of effective indices,
            in the caller's wavelength order.
        """
        wavelengths = np.asarray(wavelengths, dtype=float)
        wl_center = float(np.mean(wavelengths))

        # Cache before building anything: on a hit there is no simulation to
        # build at all. The fingerprint carries `pad` explicitly because it
        # changes the solved domain (via _build_section_mode_sim) without
        # appearing anywhere else in these kwargs.
        family_config = config = None
        if save_dir:
            call = dict(
                x0=x0, x1=x1, wavelengths=wavelengths,
                monitor_center=monitor_center, monitor_size=monitor_size,
                pad=(self.dpml if pad is None else pad),
                geometry_kwargs=geometry_kwargs,
                **resolved_material_config(
                    self.material_stack, wavelengths,
                    self.material_stack.resolve_dispersive(self.dispersive)),
            )
            # family = everything except the mode count, so a stored 14-mode
            # solve can serve a 6-mode request by slicing (lookup_superset).
            # Legitimate here in a way it is not for solve_section(): this
            # method's monitor plane is pinned by monitor_center/monitor_size
            # and does not widen with n_modes, and band i is whatever MPB's
            # i-th distinct band is regardless of how many were requested.
            family_config = build_run_config(
                self, "section_mode_neff", dict(call), geometry)
            config = build_run_config(
                self, "section_mode_neff", dict(call, n_modes=n_modes), geometry)
            hit = lookup(save_dir, config)
            if hit is not None:
                logger.info("section_mode_neff: cache hit at %s", hit)
                return load_result(hit)
            superset = lookup_superset(save_dir, family_config, (n_modes,))
            if superset is not None:
                folder, _stored = superset
                logger.info("section_mode_neff: superset cache hit at %s (sliced)", folder)
                return load_result(folder)[:, :n_modes]
        if cache_only:
            raise CacheMissError(
                f"section_mode_neff(): cache_only=True but no cached run found "
                f"for this config (save_dir={save_dir!r}). Run once with "
                f"cache_only=False to populate the cache."
            )

        # Hand the whole band down, not just its centre: under
        # `dispersive=True` the media are then built across it rather than
        # frozen at wl_center, which is the point of the flag for a multi-
        # wavelength caller. Inert when the band is one point, or when
        # dispersion is off. See `_build_section_mode_sim()`.
        sim, where, mea_c, mea_w = self._build_section_mode_sim(
            geometry, wl_center, x0, x1, monitor_center, monitor_size, pad=pad,
            wavelengths=wavelengths, **geometry_kwargs)

        neff = np.zeros((len(wavelengths), n_modes))
        for i, wl in enumerate(wavelengths):
            f = 1.0 / float(wl)
            n_clad = self.material_stack.sio2_index(float(wl))
            # kpoint only seeds the solve; match_frequency then iterates to the
            # true k. Any guided mode sits between the cladding and core
            # indices, so the midpoint of that range is a seed that cannot be
            # badly wrong for any waveguide this solver handles.
            # _iter_distinct_bands() guards against get_eigenmode() returning
            # the same eigenvector for two different band numbers (see its
            # docstring / gap 6) -- band as yielded may skip some meep
            # band_nums, but always fills n_modes distinct logical positions.
            for band, em in enumerate(
                self._iter_distinct_bands(sim, where, f, mea_c, mea_w, n_modes, n_clad), start=1
            ):
                neff[i, band - 1] = float(np.real(em.k.x)) / f

        sim.reset_meep()

        if save_dir:
            folder = save_run(save_dir, config, neff)
            update_family_pointer(save_dir, family_config, (n_modes,),
                                  os.path.basename(folder))
            logger.info("section_mode_neff: cache saved at %s", folder)
        return neff

    def _build_section_mode_sim(
        self, geometry, wl_center: float, x0: float, x1: float,
        monitor_center: Optional[float], monitor_size: Optional[float],
        pad: Optional[float] = None,
        wavelengths: Optional[np.ndarray] = None,
        **geometry_kwargs,
    ):
        """Shared preamble for section_mode_neff()/section_mode_fields(): the
        minimal self-contained meep Simulation for this section's
        cross-section, initialized and ready for get_eigenmode() calls, plus
        the eigenmode monitor Volume. Caller owns sim.reset_meep().

        Split out of section_mode_neff() verbatim (no behavior change) so
        section_mode_fields() can reuse the identical sim setup rather than
        duplicating it.

        `wavelengths` (µm, default None -> use `wl_center` alone): the caller's
        whole band, when it has one. Media are built from it via `_media()`, so
        a band caller under `dispersive=True` gets per-wavelength media instead
        of one index frozen at `wl_center`.

        Why this parameter exists at all. Until 2026-07-31 this method built its
        media with `si_medium(wl_center)` directly, bypassing `_media()` and so
        never consulting `self.dispersive`. That was justified as "a mode solve
        is single-wavelength, where the frozen index at that exact wavelength is
        already the correct medium" (`material_model_validation.md` work item 2)
        — true for `section_mode_fields()`, which takes one wavelength, but NOT
        for `section_mode_neff()`, which accepts a whole band and loops it on
        this one structure. For that caller the frozen centre index is the exact
        defect work item 2 fixed for `_build_domain()`/`solve_section()`, left
        in place here by an argument that did not apply. The measured
        consequence: `FDTDSolver(dispersive=True).section_mode_neff(g, band, …)`
        returned n_eff byte-identical to `dispersive=False` — max |Δn_eff| =
        0.000e+00 — while the media differ by 0.031 in ε (12.063 vs 12.095 at
        1.60 µm on SOI_220). Job 11586986 swept both settings for ~1 h and
        produced two identical columns because of it.

        Scalar callers are byte-unchanged, and not by convention — `_medium()`
        needs BOTH the flag and a real band (`size > 1 and min < max`), so a
        scalar falls through to the frozen index under either setting.

        `pad` (µm, default None -> `self.dpml`): straight-material clearance
        drawn past x0/x1 on each side, and the `mp.PML()` thickness of the
        cell built around it. Only `self.dpml` was ever used here before this
        parameter existed, so the default is behavior-preserving.

        Why this can safely be smaller than `self.dpml` for a mode-field
        query (unlike everywhere else `self.dpml` appears in this file):
        `get_eigenmode()` runs MPB on the initialized structure and returns
        without any time-stepping (see `section_mode_neff()`'s docstring), so
        the PML absorbing boundary is never actually exercised for this call
        — it exists only so `mp.Simulation`'s cell/boundary_layers construct
        without raising. The clearance's real job is purely geometric: giving
        `geometry.section_prisms()` enough straight run past x0/x1 that the
        cut sits inside uniform material, not into a taper/bend. `self.dpml`
        is simply the value already in hand for that; it is not a physical
        requirement of the mode solve itself. This is what lets
        `HybridSolver(check_mode_correspondence=True)`'s diagnostic query
        (`mode_correspondence_pad`, threaded through `section_mode_fields()`
        below) fit inside a device's real straight run-out, which can be
        thinner than `self.dpml` — confirmed the actual blocker: the coupler
        and both MMI hybrid examples use `body_stub = 0.25` against
        `dpml = 0.5` (`fdtd_vs_EME_mode_validation.md` item 5).
        """
        import meep as mp

        clearance = self.dpml if pad is None else float(pad)

        si_med, sio2_med = self._media(
            wl_center if wavelengths is None else wavelengths)

        prisms = geometry.section_prisms(
            x0, x1, si_med,
            stub_left=clearance, stub_right=clearance,
            wg_height=self._wg_height,
            **geometry_kwargs,
        )
        span = x1 - x0
        Sx = span + 2 * clearance
        Sy = self.sim_width + 2 * clearance
        lo, hi = x0 - clearance, x1 + clearance
        x_center = (lo + hi) / 2.0
        # Same vertex-wise shift into meep's origin-centred cell that
        # solve_section() does — prisms carry vertices, not a centre.
        shifted = [
            mp.Prism(
                vertices=[mp.Vector3(v.x - x_center, v.y, v.z) for v in prism.vertices],
                height=prism.height, material=prism.material,
            )
            for prism in prisms
        ]

        mea_c = 0.0 if monitor_center is None else float(monitor_center)
        mea_w = self.sim_width if monitor_size is None else float(monitor_size)

        sim = mp.Simulation(
            resolution=self.res,
            cell_size=mp.Vector3(Sx, Sy, self._cell_sz),
            boundary_layers=[mp.PML(clearance)],
            geometry=shifted,
            default_material=sio2_med,
            eps_averaging=self.eps_averaging,
            subpixel_tol=self.subpixel_tol,
            subpixel_maxeval=self.subpixel_maxeval,
        )
        sim.init_sim()

        where = mp.Volume(center=mp.Vector3(0, mea_c),
                          size=mp.Vector3(0, mea_w, self._mon_sz))
        return sim, where, mea_c, mea_w

    def _iter_distinct_bands(self, sim, where, f, mea_c, mea_w, n_modes, n_clad):
        """Yields `n_modes` distinct-eigenvector ``EigenmodeData`` objects
        for consecutive logical band positions, silently skipping any meep
        `band_num` whose field is (numerically) IDENTICAL to an
        already-yielded one.

        Guards a real meep/MPB failure mode, found auditing the modal
        multiplexer's 16-mode interface basis (see
        validation_markdown/modal_multiplexer_validation.md gap 6):
        ``sim.get_eigenmode()`` can return the SAME eigenvector for two
        different requested `band_num` values when they sit close together
        in frequency — confirmed directly (self-overlap 0.9999999999 between
        two "different" bands, on a cross-section where every other band
        pair overlapped ≤ 0.126, and the other engine's own spectrum there
        has no ties at all — the signature of a duplicate solve, not a
        genuine physical degeneracy). Silently accepting both as distinct
        makes an `n_modes`-mode basis rank-deficient by one per duplicate —
        `HybridSolver`'s index-*i*-to-index-*i* contraction loses that
        channel entirely, not just mis-indexes it.

        The duplicate check uses its own small, fixed diagnostic grid
        (`_DUP_BAND_DIAG_NY` x `_DUP_BAND_DIAG_NZ`), independent of
        `self.res` or any caller-requested grid — it only needs to tell
        "same eigenvector" from "different eigenvector" (0.9999... vs
        ≤ 0.13 in the case this was found on), not resolve mode shape
        precisely, so it stays cheap regardless of how finely a caller like
        `section_mode_fields()` wants its own returned fields sampled.
        `field_overlap_matrix()`'s `dy`/`dz` cancel exactly in its
        normalized ratio, so passing 1.0 for both here (rather than this
        grid's real spacing) gives the identical answer.

        **2026-07-31 — two more failure modes closed**
        (`fdtd_vs_EME_mode_validation.md` items 12 and 14), both verified
        offline against a real device's own measured band trace
        (`diagnose_kpoint_ladder_v2.py`, job 11588309) before landing here:

        - **The k-seed now ladders instead of staying fixed.** Every attempt
          used to start `get_eigenmode()`'s k-search from the same
          `_KPOINT_SEED_NEFF` guess — the textbook shape of the duplicate
          this guard works around, since two adjacent targeted solves
          launched from an identical point can converge to the same
          eigenpair. Seeding attempt *k* from attempt *(k-1)*'s own raw
          (accepted-or-skipped) converged n_eff instead recovers a genuinely
          distinct mode the fixed seed silently drops, at no extra attempt
          cost (measured: same attempt count, one fewer skip, on the
          cross-section this was found on — three seeding variants tried,
          all equivalent, so the simplest one is what ships).
        - **The walk now stops at the guided edge, not just at
          `max_band_num`.** A `band_num` with no solution at the target
          frequency does not raise — MPB's k-search limit-cycles between two
          k values forever, hanging the process (measured: 161 solves on one
          band, 59% of a job's entire eigensolve time, before SLURM killed
          it). Every raw solve's n_eff is now checked against `n_clad`
          *before* the duplicate-overlap sampling: once a raw solve returns
          `n_eff <= n_clad`, nothing above it is guided and the *next*
          `band_num` is where a real cross-section's non-convergent region
          begins, so the walk stops there rather than reaching it. If fewer
          than `n_modes` were accepted by then, this raises immediately —
          a basis silently padded with unguided modes has a rank that is a
          lie (see `fdtd_vs_EME_mode_validation.md` item 9 for why an
          unguided column is worthless to `HybridSolver` anyway), and
          failing loudly here beats discovering it downstream.

        Args:
            n_clad: cladding refractive index at this cross-section's
                wavelength (``self.material_stack.sio2_index(wl)``) — the
                guided/unguided boundary the stop condition above checks
                against.
        """
        import meep as mp

        from ..mode_identity import field_overlap_matrix

        diag_y_edges = np.linspace(mea_c - mea_w / 2, mea_c + mea_w / 2, _DUP_BAND_DIAG_NY + 1)
        diag_y = 0.5 * (diag_y_edges[:-1] + diag_y_edges[1:])
        diag_z_edges = np.linspace(-self._mon_sz / 2, self._mon_sz / 2, _DUP_BAND_DIAG_NZ + 1)
        diag_z = 0.5 * (diag_z_edges[:-1] + diag_z_edges[1:])

        accepted_ey: List[np.ndarray] = []
        accepted_ez: List[np.ndarray] = []
        band_num = 0
        n_accepted = 0
        last_raw_neff: Optional[float] = None
        max_band_num = n_modes + _DUP_BAND_MAX_EXTRA
        while n_accepted < n_modes:
            band_num += 1
            if band_num > max_band_num:
                raise RuntimeError(
                    f"FDTDSolver: could not find {n_modes} distinct eigenmodes "
                    f"within {max_band_num} band_num attempts at this cross-"
                    f"section -- meep/MPB returned a duplicate more times than "
                    f"this guard tolerates. See validation_markdown/"
                    f"modal_multiplexer_validation.md gap 6."
                )
            seed_neff = last_raw_neff if last_raw_neff is not None else _KPOINT_SEED_NEFF
            em = sim.get_eigenmode(
                f, mp.X, where, band_num, mp.Vector3(f * seed_neff, 0, 0),
                parity=mp.NO_PARITY,
            )
            n = float(np.real(em.k.x)) / f
            last_raw_neff = n
            if n <= n_clad:
                raise RuntimeError(
                    f"FDTDSolver: this cross-section supplies only {n_accepted} "
                    f"guided mode(s) of the {n_modes} requested -- band_num="
                    f"{band_num} returned n_eff={n:.6f}, at or below the "
                    f"cladding index {n_clad:.6f}. Nothing above the guided "
                    f"edge is a real channel, and walking further risks a "
                    f"non-convergent MPB k-search (see fdtd_vs_EME_mode_"
                    f"validation.md item 14). Reduce n_modes to what this "
                    f"cross-section actually supports, or choose a "
                    f"resolution where it supports more."
                )
            ey = np.array([[em.amplitude(mp.Vector3(0, float(y), float(z)), mp.Ey)
                            for z in diag_z] for y in diag_y])
            ez = np.array([[em.amplitude(mp.Vector3(0, float(y), float(z)), mp.Ez)
                            for z in diag_z] for y in diag_y])
            if n_accepted > 0:
                existing = dict(Ey=np.stack(accepted_ey), Ez=np.stack(accepted_ez))
                candidate = dict(Ey=ey[None], Ez=ez[None])
                eta = field_overlap_matrix(existing, candidate, dy=1.0, dz=1.0)
                worst = float(np.max(eta))
                if worst > _DUP_BAND_OVERLAP_THRESHOLD:
                    logger.warning(
                        "FDTDSolver: meep band_num=%d duplicates an already-"
                        "accepted eigenmode (overlap %.6f > %.2f threshold) -- "
                        "skipping and trying the next band_num. See "
                        "validation_markdown/modal_multiplexer_validation.md gap 6.",
                        band_num, worst, _DUP_BAND_OVERLAP_THRESHOLD,
                    )
                    continue
            accepted_ey.append(ey)
            accepted_ez.append(ez)
            n_accepted += 1
            yield em

    def section_mode_fields(
        self,
        geometry,
        wavelength: float,
        x0: float,
        x1: float,
        n_modes: int,
        y_grid: np.ndarray,
        z_grid: np.ndarray,
        monitor_center: Optional[float] = None,
        monitor_size: Optional[float] = None,
        pad: Optional[float] = None,
        **geometry_kwargs,
    ) -> Dict[str, np.ndarray]:
        """Effective indices AND transverse field profiles of the first
        `n_modes` modes at this section's cross-section, sampled at the
        caller's own ``(y_grid, z_grid)`` points.

        The field-based counterpart of ``section_mode_neff()``. n_eff is a
        single scalar per mode, and on a wide multimode cross-section several
        modes can sit close enough together that nearest-n_eff matching
        against meow picks the wrong pairing outright — observed directly on
        the modal multiplexer's 16-mode basis, where bands 8, 9 and 10 all
        nearest-matched to the same meow band (see
        validation_markdown/modal_multiplexer_validation.md gap 1). Comparing
        field *shape* instead of the propagation constant resolves that
        ambiguity and, unlike n_eff matching, does not require the two
        engines/resolutions to already agree numerically to work — see
        ``beef.mode_identity.field_overlap_matrix()``/``best_correspondence()``,
        the counterparts that consume this method's output.

        No meep object leaks past this method — fields are meep's own
        eigenmode amplitude at each requested point, in ``beef``'s device
        frame. Sample on whatever grid the *other* engine's mode sits on
        (meow's ``Mode`` fields sit at its ``Mesh2D``'s cell centers, not its
        node coordinates — see ``beef.mode_identity.mesh_cell_centers()``) so
        the two are directly comparable point-for-point.

        Args:
            wavelength: single wavelength (µm) — unlike ``section_mode_neff()``,
                this does not sweep: field grids are large, and a field
                comparison only needs one representative wavelength.
            y_grid, z_grid: 1-D device-frame sample points (µm), spanning the
                monitor plane's y (in-plane transverse) and z (vertical) axes.
            x0, x1, n_modes, monitor_center, monitor_size, geometry_kwargs: same
                meaning as ``section_mode_neff()``.
            pad: see ``_build_section_mode_sim()``. None (default) is
                ``self.dpml``, unchanged from before this parameter existed.
                ``HybridSolver``'s ``mode_correspondence_pad`` is the intended
                caller of a non-default value: the correspondence check's own
                query window is only ``2 * _CORR_HALF_SPAN`` wide (0.1 µm),
                so it needs far less clearance than a real PML footprint to
                stay inside a thin interior stub — see this method's calling
                context in ``fdtd_vs_EME_mode_validation.md`` item 5.

        Returns:
            dict with ``'neff'`` (real, shape ``[n_modes]``) and
            ``'Ex','Ey','Ez','Hx','Hy','Hz'`` (complex, each shape
            ``[n_modes, len(y_grid), len(z_grid)]``), band order as meep
            returns it.
        """
        import meep as mp

        wl_center = float(wavelength)
        sim, where, mea_c, mea_w = self._build_section_mode_sim(
            geometry, wl_center, x0, x1, monitor_center, monitor_size, pad=pad, **geometry_kwargs)

        y_grid = np.asarray(y_grid, dtype=float)
        z_grid = np.asarray(z_grid, dtype=float)
        components = dict(Ex=mp.Ex, Ey=mp.Ey, Ez=mp.Ez, Hx=mp.Hx, Hy=mp.Hy, Hz=mp.Hz)
        fields = {name: np.zeros((n_modes, len(y_grid), len(z_grid)), dtype=complex)
                  for name in components}
        neff = np.zeros(n_modes)

        f = 1.0 / wl_center
        n_clad = self.material_stack.sio2_index(wl_center)
        # _iter_distinct_bands() guards against get_eigenmode() returning the
        # same eigenvector for two different band numbers (see its docstring
        # / validation_markdown/modal_multiplexer_validation.md gap 6).
        for band, em in enumerate(
            self._iter_distinct_bands(sim, where, f, mea_c, mea_w, n_modes, n_clad), start=1
        ):
            neff[band - 1] = float(np.real(em.k.x)) / f
            for iy, y in enumerate(y_grid):
                for iz, z in enumerate(z_grid):
                    point = mp.Vector3(0, float(y), float(z))
                    for name, comp in components.items():
                        fields[name][band - 1, iy, iz] = em.amplitude(point, comp)

        sim.reset_meep()
        return dict(neff=neff, **fields)

    def section_mode_beta(
        self,
        geometry,
        wavelengths: np.ndarray,
        x0: float,
        x1: float,
        n_modes: int,
        lengths: Optional[Sequence[float]] = None,
        monitor_center: Optional[float] = None,
        monitor_size: Optional[float] = None,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        use_merged_domain: bool = True,
        seed_anchors: Optional[int] = 3,
        **geometry_kwargs,
    ) -> np.ndarray:
        """Propagation constants (rad/µm) per mode per wavelength — see
        `section_mode_beta_details()`, of which this is the facade.

        Returns only the fitted β, which is all any propagation caller needs.
        Call `section_mode_beta_details()` instead for the measurement behind
        it (fit lengths, raw per-length phases, MPB seed, per-cell `BetaFit`),
        e.g. to plot the fits — same arguments, same cost, same guards.
        """
        return self.section_mode_beta_details(
            geometry, wavelengths, x0, x1, n_modes,
            lengths=lengths,
            monitor_center=monitor_center, monitor_size=monitor_size,
            save_dir=save_dir, cache_only=cache_only,
            use_merged_domain=use_merged_domain, seed_anchors=seed_anchors,
            **geometry_kwargs,
        ).beta

    def section_mode_beta_details(
        self,
        geometry,
        wavelengths: np.ndarray,
        x0: float,
        x1: float,
        n_modes: int,
        lengths: Optional[Sequence[float]] = None,
        monitor_center: Optional[float] = None,
        monitor_size: Optional[float] = None,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        use_merged_domain: bool = True,
        seed_anchors: Optional[int] = 3,
        strict: bool = True,
        **geometry_kwargs,
    ) -> "BetaCalibration":
        """Propagation constants this solver's FDTD *actually accumulates*
        through the section's cross-section (rad/µm), per mode per wavelength.

        The counterpart to `section_mode_neff()`, and deliberately not the same
        number. `section_mode_neff()` asks MPB for the mode index of the
        discretized cross-section; this runs the time-stepper over several
        lengths and fits the phase it really picks up. On the coupler at res 40
        they differ by ~2.3% per mode — MPB 9.42255 vs FDTD 9.63292 rad/µm for
        the even supermode — because MPB's planewave operator and meep's Yee
        operator are different discretizations of the same continuum problem.

        **Why that 2.3% matters far more than it looks.** A hybrid chain's
        coupling is set by Δβ, a *difference* of two β near 9.4. Both modes are
        biased the same way, and the differential between them is only 0.129
        percentage points — but ``β/Δβ ≈ 99``, so it lands on Δβ as a 10.4%
        error. That is the whole of dc_validation.md item 6, and it is why the
        stub sweep never converges: its L₅₀ drift is exactly
        ``2s·(1 - Δβ_FDTD/Δβ_EME)``. Closing it by mesh refinement would need
        res ≈ 143, some 164× the 3-D cost.

        Feeding this into the EME body instead makes the cascade
        self-consistent — see `HybridSolver(calibrate_body_neff=True)`. That is
        not a fudge: the hybrid's job is to reproduce a full-device FDTD run,
        and that reference *is* meep's discretization, so propagating the body
        on MPB's β is the actual mismatch.

        Cost: ``len(lengths)`` real `solve_section()` runs — genuinely
        expensive, unlike `section_mode_neff()`. Pass `save_dir` so a given
        cross-section is paid for once.

        Args:
            x0, x1: device-frame bounds. Their midpoint sets the cross-section,
                and the *whole* span [x0, x1] must be uniform — a section that
                clips a taper or bend makes the phase non-linear in L, which
                shows up in the returned residual rather than silently.
            lengths: section lengths to fit over (µm). Default is three,
                evenly spaced, with the spacing capped by β rather than by the
                section's size — see the body for why both of those matter.
                Two lengths determine a slope but make the residual identically
                zero, disabling the linearity guard; pass two only if you have
                another reason to trust the branch.
            monitor_center / monitor_size: as `solve_section()`. The default
                here pins the full non-PML width on the ``y=0`` axis, which is
                what a coupled pair's even/odd basis needs.
            use_merged_domain: if True (the default) and no caching is
                requested (`save_dir is None` and `cache_only` is False), fit
                β from ONE shared-domain `solve_section_multi_length()` call
                instead of `len(lengths)` independent `solve_section()` runs —
                validated equivalent to float-noise precision (β agrees to
                ~0.0003%, both guards' margins unaffected), see
                `validation_markdown/section_mode_beta_multi_monitor_
                validation.md` steps 4-6. Both post-fit guards below (and
                Guard 2, above) apply identically either way, since they read
                `phases` after it is filled, not how it was filled.

                If `save_dir`/`cache_only` ARE requested, this method
                automatically falls back to the per-length path regardless of
                this flag: the merged path has no cache-key design of its own
                yet (that doc's own open question), so every caller that
                wants caching — which today means every real production
                caller, `HybridSolver` and every hybrid example — keeps using
                the proven, cacheable per-length path, unaffected. Set to
                `False` explicitly to force the per-length path even without
                caching, e.g. to re-validate against a future finding (mirrors
                `mmi_2x2_hybrid.py`/`directional_coupler_hybrid.py`'s own
                `SYMMETRIC_MODE_REUSE` fallback pattern).

            seed_anchors: how many wavelengths the MPB seed is actually
                solved at, the rest being interpolated. 3 (default, and active
                for any caller with more than 3 wavelengths); None or 0 solves
                the seed on the caller's full grid, as before this existed.

                **Only the seed is affected.** The FDTD time-stepping and the
                modal decode always run on every wavelength passed in — this
                changes the *estimate* used to pick the 2π branch and to space
                the fit lengths, never the measurement.

                The seed can afford it because both of its jobs are coarse.
                Branch selection needs the estimate within ~25% of β
                (`_BETA_MIN_BRANCH_FRACTION`), and production margins sit at
                0.05-0.095 of the half-distance to the next branch. A 3-anchor
                quadratic reproduced the MMI 2×2's real 14-mode C-band n_eff
                curve to 1.10e-04 worst case (modes 0-12: ≤ 9.1e-06; mode 13,
                nearest cutoff, is the outlier) = 4.5e-04 rad/µm of β_seed
                = 1.6e-04 of that half-distance.

                Cost: linear in this number. On the MMI 2×2 it takes the seed
                from 154 eigensolves to 42, ~31.5 min to ~8.6 min. Raise it if
                a device's n_eff(λ) is not well described by a degree-(n−1)
                polynomial across its band — a strongly non-parabolic curve, a
                very wide band, or a mode crossing inside it — or set None to
                switch the subsampling off entirely.

                Visible in the returned calibration as `seed_wavelengths`, and
                the one reported quantity it perturbs is `seed_excess`
                (β/β_seed − 1): ~0.007 percentage points on a 2.5% excess at
                the worst mode.
            strict: True (default) raises on either post-fit guard below,
                exactly as `section_mode_beta()` always has. False downgrades
                both to a WARNING and returns the calibration anyway, with the
                messages in `.guard_failures` — **the returned β is then not
                usable for propagation**, and the only reason this exists is
                that a failing fit is precisely the one worth plotting: the
                residual guard reduces the whole measurement to one number and
                aborts, which is what sent `diagnose_mmi1x2_calibration_
                phase_fit.py` to the on-disk cache to reconstruct by hand.
                Guard 2 (pre-simulation, branch spacing) is unaffected and
                still raises either way — it rejects the *request* before
                spending FDTD time, so there is nothing measured to inspect.

        Returns:
            `BetaCalibration` — `.beta` is the real ndarray
            ``[len(wavelengths), n_modes]`` of β in rad/µm that
            `section_mode_beta()` hands back; the rest is the measurement it
            came from (fit lengths, raw phases, seed, per-cell `BetaFit`),
            kept so plotting/diagnostics need not re-derive any of it. Check
            `.guard_failures` is empty before trusting `.beta` when
            `strict=False`.

        Raises:
            ValueError: if the phase-fit residual exceeds
                `_BETA_FIT_RESIDUAL_MAX` for some cell — the phase is not linear
                in L, so the cut is probably not uniform; or if the
                branch-selection margin exceeds `_BETA_BRANCH_MARGIN_RAISE`,
                i.e. the lengths are spaced too widely for MPB's seed to
                identify which wrap the FDTD phase is on (fix by passing
                closer-spaced `lengths`). Both messages name the offending
                (wavelength, mode) cell.

                Until 2026-08-01 this block promised the margin raise and did
                not implement it, while the residual raise it did not mention
                was firing for a reason that was really a bug here — see the
                comment above `_BETA_MIN_BRANCH_FRACTION`.
        """
        span = float(x1 - x0)
        wavelengths = np.asarray(wavelengths, dtype=float)
        mid = 0.5 * (x0 + x1)
        iface = dict(
            monitor_center=0.0 if monitor_center is None else monitor_center,
            monitor_size=self.sim_width if monitor_size is None else monitor_size,
        )

        # MPB's answer, used only to pick which wrap each phase is on. Cheap
        # (no time-stepping), and it is also what sets a safe default spacing
        # below, so it has to come first.
        #
        # Solved on a SHORT window at the midpoint rather than on [x0, x1]:
        # section_mode_neff() extends its cut by dpml on each side to make PML
        # room, so handing it the full span of a body that runs right up to the
        # S-bends pushes the cut into them and DirectionalCoupler rejects it.
        # The cross-section is identical anywhere in a uniform region, so a
        # narrow window is not an approximation — it is the same modes, cheaper.
        seed_span = min(2.0, span - 2.0 * self.dpml)
        if seed_span <= 0.0:
            raise ValueError(
                f"section_mode_beta: section [{x0}, {x1}] is only {span:.3f} µm "
                f"long, leaving no room for the {self.dpml} µm PML clearance the "
                f"mode solve needs on each side."
            )
        # The seed is solved on `seed_anchors` wavelengths and fitted to the
        # rest, rather than on the caller's whole grid. Its two jobs are coarse
        # decisions -- which 2*pi branch each measured phase sits on, and (via
        # beta_max) how far apart the fit lengths may be -- while n_eff(lambda)
        # across a device band is smooth and very nearly quadratic. Measured on
        # this device's own 14-mode C-band curve, a 3-anchor quadratic
        # reproduces the 11-point n_eff to 1.10e-04 worst case (modes 0-12 are
        # <= 9.1e-06; mode 13, nearest cutoff and so most curved, is the
        # outlier). That is 4.5e-04 rad/um on beta_seed, or 1.6e-04 of the
        # half-distance to the neighbouring branch -- against production branch
        # margins of 0.05-0.095, three orders of magnitude of headroom.
        #
        # ONLY the seed is subsampled. The FDTD time-stepping and the modal
        # decode below both run on every wavelength the caller passed; nothing
        # downstream of seed_beta ever sees a reduced grid. The one visible
        # consequence is that `seed_excess` (beta/beta_seed - 1, the
        # discretization gap this calibration exists to report) is computed
        # against an interpolated beta_seed -- shifting a 2.5% excess by
        # ~0.007 percentage points at the worst mode.
        #
        # save_dir/cache_only ride along: the seed is the single largest item
        # in a calibration once the decode is cut down (154 eigensolves,
        # ~31.5 min on the MMI 2x2), so a re-run that hits cache for every
        # simulation should not pay it again, and a documentation/plot re-run
        # under cache_only must not silently spend half an hour here.
        wl_seed = wavelengths
        if (seed_anchors is not None and seed_anchors >= 2
                and len(wavelengths) > seed_anchors
                and float(np.ptp(wavelengths)) > 0.0):
            wl_seed = np.linspace(
                float(np.min(wavelengths)), float(np.max(wavelengths)), seed_anchors)

        seed_at = self.section_mode_neff(
            geometry, wl_seed, mid - seed_span / 2.0, mid + seed_span / 2.0,
            n_modes,
            monitor_center=iface["monitor_center"],
            monitor_size=iface["monitor_size"],
            save_dir=save_dir, cache_only=cache_only,
            **geometry_kwargs,
        )

        if wl_seed is wavelengths:
            seed = seed_at
        else:
            # Degree len-1 through len points is exact interpolation, not a
            # least-squares compromise. x is centred on the band so the fit is
            # well conditioned at lambda ~ 1.5 with a 0.035 um span.
            x_mid = float(np.mean(wl_seed))
            seed = np.empty((len(wavelengths), n_modes))
            for m in range(n_modes):
                coeffs = np.polyfit(wl_seed - x_mid, seed_at[:, m], len(wl_seed) - 1)
                seed[:, m] = np.polyval(coeffs, wavelengths - x_mid)
            logger.info(
                "section_mode_beta: seed solved at %d of %d wavelengths and "
                "interpolated (degree %d); FDTD and decode still use all %d.",
                len(wl_seed), len(wavelengths), len(wl_seed) - 1, len(wavelengths),
            )
        seed_beta = 2.0 * np.pi * seed / wavelengths[:, None]
        beta_max = float(np.max(np.abs(seed_beta)))
        max_safe_spacing = 2.0 * np.pi / (_BETA_MIN_BRANCH_FRACTION * beta_max)

        if lengths is None:
            # NOT span/2 and span: for a 30 µm body that spaces the lengths
            # 15 µm apart, putting the 2π branches 0.42 rad/µm apart on a
            # β of 9.6 — far closer together than the disagreement being
            # measured, so the seed could not identify the branch at all.
            # Spacing is capped by the physics, not by the section's size.
            #
            # THREE lengths, not two: with two points the straight-line fit is
            # exact by construction, so the residual is identically zero and
            # guard 1 below can never fire. A third point is what gives it
            # something to disagree with, and these sections are short (~1.3 µm
            # at the coupler) so the extra solve is cheap.
            #
            # -2*dpml, mirroring seed_span above: solve_section() below draws
            # each length's own window padded by dpml on each side, so the
            # longest length (3*spacing) must leave that much clearance inside
            # the uniform span or it pushes the drawn cut into whatever flanks
            # this section (a taper/bend). usable_span > 0 is already
            # guaranteed here — seed_span's check above raises first on the
            # identical span - 2*dpml <= 0 condition. Bites short bodies
            # first: unreachable at the ~29.5 µm bodies both live callers use,
            # but reproduced directly on a DirectionalCoupler(coupling_length=3.0).
            usable_span = span - 2.0 * self.dpml
            spacing = min(usable_span / 3.0, max_safe_spacing / 2.0)
            lengths = [spacing, 2.0 * spacing, 3.0 * spacing]
        lengths = [float(L) for L in lengths]
        if len(lengths) < 2:
            raise ValueError(
                "section_mode_beta needs at least 2 lengths: one length cannot "
                "separate the propagation phase from the constant mode-matching "
                "phase at the monitor planes."
            )

        # Guard 2 (see the constants): a whole-branch seed error is invisible in
        # the fitted numbers -- clean residual, zero margin, β 33% wrong -- so it
        # has to be excluded before simulating rather than detected after.
        max_gap = float(np.max(np.diff(np.sort(lengths))))
        if max_gap > max_safe_spacing:
            raise ValueError(
                f"section_mode_beta: lengths {lengths} are spaced up to "
                f"{max_gap:.3f} µm apart, which puts the 2π branches only "
                f"{2*np.pi/max_gap:.3f} rad/µm apart on a β of ~{beta_max:.3f}. "
                f"MPB's seed could then be a whole branch out and the fit would "
                f"look perfect at a badly wrong β. Space them no more than "
                f"{max_safe_spacing:.3f} µm apart."
            )

        # The merged path now carries its own cache key (see
        # solve_section_multi_length), so caching no longer forces the
        # per-length fallback and this is the path every caller takes by
        # default. Set use_merged_domain=False for the per-length path.
        merged_active = use_merged_domain

        if merged_active:
            # solve_section_multi_length() -- section_mode_beta_multi_monitor_
            # validation.md steps 3-6: one shared domain, one decay wait per
            # excited band, instead of len(lengths) independent domain builds.
            T_multi = self.solve_section_multi_length(
                geometry, wavelengths, mid, lengths, n_modes,
                measure_diagonal_only=True, save_dir=save_dir, cache_only=cache_only,
                **iface, **geometry_kwargs,
            )
            phases = np.angle(np.diagonal(T_multi, axis1=2, axis2=3))
        else:
            phases = np.zeros((len(lengths), len(wavelengths), n_modes))
            for j, L in enumerate(lengths):
                T = self.solve_section(
                    geometry, wavelengths, mid - L / 2.0, mid + L / 2.0,
                    n_modes_in=n_modes, n_modes_out=n_modes,
                    measure_diagonal_only=True,
                    save_dir=save_dir, cache_only=cache_only,
                    **iface, **geometry_kwargs,
                )
                phases[j] = np.angle(np.diagonal(T, axis1=1, axis2=2))

        # Full grids, not running scalar maxima. The old code kept only
        # max()-so-far, so when job 11589770 blew the residual guard the message
        # could say "2.101 rad" but not WHICH (wavelength, mode) cell — and the
        # answer (4 of 88 cells, all at the short-wavelength end) turned out to be
        # the most informative fact about the failure. Keeping the grid costs
        # n_wl*n_modes floats.
        beta = np.zeros((len(wavelengths), n_modes))
        grid_resid = np.zeros((len(wavelengths), n_modes))
        grid_margin = np.zeros((len(wavelengths), n_modes))
        grid_excess = np.zeros((len(wavelengths), n_modes))
        grid_turn_amb = np.zeros((len(wavelengths), n_modes))
        n_shifted = 0
        fits: List[List[BetaFit]] = []      # [i][m], retained for BetaCalibration
        for i in range(len(wavelengths)):
            fits.append([])
            for m in range(n_modes):
                fit = beta_fit_from_phases(
                    lengths, phases[:, i, m], seed_beta[i, m])
                fits[i].append(fit)
                beta[i, m] = fit.beta
                grid_resid[i, m] = fit.resid
                grid_margin[i, m] = fit.margin
                grid_excess[i, m] = fit.seed_excess
                grid_turn_amb[i, m] = fit.turn_ambiguity
                n_shifted += int(fit.turns_shifted)

        def _cell(grid):
            """(worst value, 'λ=… mode=…') for the largest entry of `grid`."""
            i, m = np.unravel_index(int(np.argmax(np.abs(grid))), grid.shape)
            return float(grid[i, m]), f"λ={wavelengths[i]:.4f} µm mode {m}"

        worst_resid, at_resid = _cell(grid_resid)
        worst_margin, at_margin = _cell(grid_margin)
        worst_excess, at_excess = _cell(grid_excess)
        worst_turn_amb, at_turn_amb = _cell(grid_turn_amb)

        # Guard: the phases must actually lie on a line. Since 2026-08-01 this
        # means exactly one thing — the section is not uniform, the near field
        # has not settled, or the modes are mixing. Inconsistent wrap counts used
        # to land here too and no longer do (see beta_fit_from_phases), so
        # turn_ambiguity is reported alongside to separate "an integer turn was
        # genuinely undecidable" from "the phase is smoothly curved".
        guard_failures: List[str] = []

        def _guard_failed(message: str) -> None:
            """Raise (strict, the default and every propagation caller's path)
            or record-and-warn, so a failing fit can still be returned for
            inspection — see `strict` in the Args."""
            if strict:
                raise ValueError(message)
            guard_failures.append(message)
            logger.warning("section_mode_beta (strict=False, NOT usable for "
                           "propagation): %s", message)

        if worst_resid > _BETA_FIT_RESIDUAL_MAX:
            _guard_failed(
                f"section_mode_beta: phase-fit residual {worst_resid:.3f} rad "
                f"exceeds {_BETA_FIT_RESIDUAL_MAX} at {at_resid} (worst of "
                f"{grid_resid.size} cells; {int((grid_resid > _BETA_FIT_RESIDUAL_MAX).sum())} "
                f"over threshold) with lengths {lengths}. The measured phase is "
                f"not linear in length, so [{x0}, {x1}] is most likely not a "
                f"uniform cross-section. Worst turn ambiguity {worst_turn_amb:.4f} "
                f"({at_turn_amb}) — near 0 means smooth curvature rather than an "
                f"undecidable wrap count."
            )
        # Guard: the seed has to be able to pick the branch. This is the raise
        # section_mode_beta's own docstring has always promised and never made.
        if worst_margin > _BETA_BRANCH_MARGIN_RAISE:
            _guard_failed(
                f"section_mode_beta: branch margin {worst_margin:.3f} exceeds "
                f"{_BETA_BRANCH_MARGIN_RAISE} at {at_margin} with lengths "
                f"{lengths} — the MPB seed sat too close to a branch boundary to "
                f"identify which one, so the returned β would be a coin flip. "
                f"Use more closely spaced lengths (the branches are 2π/ΔL apart, "
                f"so shorter gaps spread them further)."
            )
        if abs(worst_excess) > _BETA_SEED_EXCESS_WARN:
            logger.warning(
                "section_mode_beta: FDTD β exceeds the MPB seed by %+.2f%% at %s — "
                "out of family with the +2%%-ish discretization gap this "
                "calibration exists to measure; check the seed and the cut.",
                100.0 * worst_excess, at_excess,
            )
        logger.info(
            "section_mode_beta: %d lengths %s | worst residual %.2e rad (%s) | "
            "worst margin %.3f (%s)%s | worst seed excess %+.2f%% (%s) | "
            "worst turn ambiguity %.4f | wrap counts moved off the seed's own "
            "rounding in %d of %d cells",
            len(lengths), lengths, worst_resid, at_resid, worst_margin, at_margin,
            "" if worst_margin < _BETA_BRANCH_MARGIN_MAX
            else "  <-- moved far from MPB's seed; check this is real",
            100.0 * worst_excess, at_excess, worst_turn_amb,
            n_shifted, grid_resid.size,
        )
        return BetaCalibration(
            beta=beta, lengths=lengths, phases=phases, wavelengths=wavelengths,
            seed_beta=seed_beta, fits=fits, x0=float(x0), x1=float(x1),
            mid=float(mid), merged=bool(merged_active),
            guard_failures=tuple(guard_failures),
            seed_wavelengths=wl_seed,
        )

    # ------------------------------------------------------------------
    # Generic section solver (used by HybridSolver)
    # ------------------------------------------------------------------

    def solve_section(
        self,
        geometry,
        wavelengths: np.ndarray,
        x0: float,
        x1: float,
        n_modes_in: int = 1,
        n_modes_out: int = 1,
        stub_left: float = 0.0,
        stub_right: float = 0.0,
        monitor_center: Optional[float] = None,
        monitor_size: Optional[float] = None,
        y0: Optional[float] = None,
        y1: Optional[float] = None,
        excite_port: Optional[PortSpec] = None,
        measure_port: Optional[PortSpec] = None,
        measure_margin: Optional[float] = None,
        measure_margins: Optional[Sequence[float]] = None,
        measure_diagonal_only: bool = False,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        **geometry_kwargs,
    ) -> np.ndarray:
        """Modal transfer matrix for the device section [x0, x1] (device frame).

        One FDTD run is needed per excited mode, so the method excites from
        whichever boundary has FEWER modes — min(n_modes_in, n_modes_out) runs
        instead of always n_modes_in. By reciprocity, launching mode k
        backward from the right and measuring mode m at the left gives
        exactly the same coefficient as launching m forward from the left and
        measuring k at the right (same passive, reciprocal structure, same
        mode normalization on both sides) — so the returned matrix is
        identical regardless of which side actually did the exciting. This is
        also position-agnostic: the same method handles a true input port, a
        true output port, or a purely interior interface, since meep does not
        care whether a cross-section is a "real" device port or an arbitrary
        cut through a larger structure.

        Args:
            x0, x1: device-frame bounds of this section (µm).
            n_modes_in / n_modes_out: mode counts at the left / right boundary.
            stub_left / stub_right: straight uniform-material extension beyond
                x0 / x1 (µm) — only needed at the two ends of the full hybrid
                chain, to create PML clearance; 0 for interior sections whose
                boundaries already sit inside a uniform region.
            monitor_center: explicit y-centre (µm) for the measured/interface-
                side monitor, overriding geometry.center_at() there. Useful
                when center_at() sampled just inside the section boundary
                reports an off-axis arm rather than the axis the interface
                basis is symmetric about — e.g. DirectionalCoupler's S-bend/
                coupling-region seam, where the default inherits the excited
                arm's ∓coupling_pitch/2 centreline; pass 0.0 to re-centre the
                supermode monitor on the shared y=0 axis, matching the
                adjacent EME body's own basis window. None (default) keeps
                the geometry-driven placement.
            monitor_size: explicit y-span (µm) for that same monitor,
                overriding the default rule (sim_width when the side must
                resolve >1 mode, geometry.width_at() for a genuine single-
                mode port). Pass sim_width (with monitor_center=0.0) for a
                window spanning the full non-PML domain, PML inner edge to
                PML inner edge. None (default) keeps the rule. Neither
                override touches the source/norm monitor on the excited side
                (see source_aperture_at / half_domain_aperture for those).
            y0, y1: device-frame y-bounds for this section's domain, mirroring
                x0/x1 but for the transverse axis. Both None (default): the
                domain is sim_width wide, centred on y=0 — today's only
                behavior, unchanged for every caller that doesn't pass these.
                Both given: the domain spans exactly [y0, y1] (+ 2*dpml).
                Giving only one raises ValueError. Required whenever
                excite_port or measure_port has axis="y" (there is no
                y-extent to place a y-normal plane in otherwise).
            excite_port, measure_port: explicit PortSpec overrides for the
                excited/measured side, for a section with a non-x-normal
                port (e.g. one of a ring resonator's two coupling ports —
                see beef.geometry.ring_coupler/ring_bend). None (default,
                both): identical to today's behavior, section_fdtd_layout()
                picks the excited side by reciprocity (excite from whichever
                boundary has fewer modes) and both sides are x-normal. If
                EITHER is given, that tie-break is bypassed entirely: the
                excited side is always n_modes_in, the measured side always
                n_modes_out, regardless of counts — section_fdtd_layout()'s
                own `n_out <= n_in` tie-break silently picks the wrong
                physical side whenever n_in == n_out, which is exactly the
                case a single bend/coupling section with one mode per port
                hits. A None side (with the other given) keeps using
                section_fdtd_layout()'s geometry-derived value for that side
                only. A y-normal excite_port additionally triggers
                flux-calibrated normalization (see the method body) — an
                empirically necessary correction for a backward, off-axis-
                launched source, first found in
                examples/ring_resonator/python/ring_resonator_hybrid.py.
                excite_port.sign and measure_port.sign use OPPOSITE edge
                conventions for the same value: sign=+1 puts an
                excite_port at axis's LOW edge (a source launching +axis
                must sit upstream, physically forced) but a measure_port
                at the HIGH edge (a monitor reading +axis-forward light
                naturally sits downstream, where forward-propagating light
                is captured after crossing the section) — sign=-1 is the
                mirror of both. Only the alpha_idx/kpoint meaning (+1=
                forward/index 0, -1=backward/index 1) is shared between the
                two roles.
            measure_margin: explicit along-axis offset (µm) for
                measure_port's own monitor position, overriding
                self.monitor_margin for that one side only. None (default):
                identical to today's behavior — measure_port (when given)
                falls back to self.monitor_margin, the SAME attribute the
                excited side's own normalization monitor uses (ex_along_norm
                above), so without this override the two sides cannot be
                offset independently of one another. Pass an explicit value
                to decouple them — e.g. a fixed simulation domain whose
                measured-side monitor sweeps to an arbitrary interior
                position while the excited side's own monitor_margin stays
                untouched (mmi_2x2_body_stub_sweep_alt.py's own original
                use case: domain size held fixed by pinning x1, the
                measured monitor's distance from x1 swept via
                measure_margin instead). Ignored (and must be None) when
                measure_port is None — raises ValueError otherwise, since
                nothing would read it. Mutually exclusive with
                measure_margins (below) — only one of the two may be given.
            measure_margins: like measure_margin, but a whole list of
                along-axis offsets (µm) for measure_port, each one getting
                its OWN mode monitor within the SAME simulation — one FDTD
                run captures the mode content at every listed position
                instead of one run per position. Since sim.run()'s cost is
                dominated by stop_when_dft_decayed()'s decay wait (a
                property of the whole field, independent of how many DFT
                monitors are recording it), N monitors cost roughly the
                same wall-clock time as 1 — the motivating case is a dense
                position sweep (e.g. every 5 nm over some span) that would
                otherwise need N separate multi-hour runs
                (mmi_2x2_body_stub_sweep_alt.py's rewrite from a
                one-run-per-body_stub sweep to this). Returned array gains
                an axis for it: [len(wavelengths), len(measure_margins),
                n_modes_out, n_modes_in] instead of the usual
                [len(wavelengths), n_modes_out, n_modes_in]. None (default):
                today's single-monitor behavior/shape, completely
                unaffected. Requires measure_port (same rule as
                measure_margin); mutually exclusive with measure_margin.
                The cache "superset" mode-count-slicing reuse
                (lookup_superset(), family_config below) is skipped
                whenever this is given — that mechanism assumes the usual
                3-D return shape.
            measure_diagonal_only: compute ONLY the diagonal T_ii, leaving
                every off-diagonal entry NaN. False (default) is unchanged
                behavior in every respect.

                This is a cost knob, not a physics one, and the cost it
                cuts is the dominant one. Each run already time-steps once
                per excited band; what it then pays for is decomposing the
                measurement monitor into `n_measure` bands at every
                wavelength — `n_measure * len(wavelengths)` MPB solves,
                ~14 s each, against ~5 min of actual time-stepping. On the
                MMI 2×2 body calibration (14 bands, 11 wavelengths) that is
                154 MPB solves per run at ~36 min, i.e. 89% of the job:
                job 11872191 spent 27.9 h of its 31.4 h there and 3.5 h
                time-stepping. Run k needs exactly one of those 154 columns
                — its own, T_kk — so 93% of the dominant cost is computed
                and discarded by the caller's `np.diagonal(T)`.

                Restricting the decomposition to band k is EXACT, not an
                approximation: the discarded entries are the ones the
                caller never reads. Measured 4.0x on the 2×2 calibration
                (30.6 h -> ~7.6 h).

                Off-diagonals come back NaN rather than 0 so that any
                consumer who reads one fails loudly instead of silently
                treating "not computed" as "measured zero". Requires a
                square problem (n_excite == n_measure); raises otherwise,
                since a diagonal is not defined otherwise. The post-hoc
                duplicate-mode scan is skipped (it compares whole columns),
                and the cache's mode-count "superset" reuse is disabled for
                these runs (a diagonal cannot serve a full request). The
                REVERSE reuse does work and is automatic: a diagonal-only
                request whose full-matrix twin is already cached is served
                from it by extracting the diagonal, so turning this on
                never invalidates a cache a previous full run paid for.
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache entirely — every call
                re-runs the simulation, same as before this parameter
                existed. When set, a call whose full config (solver params +
                these args + the exact section geometry) matches a prior
                saved run loads that run's result instead of re-simulating;
                otherwise the simulation runs as usual and its result,
                config, and a geometry plot are written to save_dir for
                future calls to match against.
            cache_only: if True, never run the simulation — a cache hit
                returns as usual, a miss raises CacheMissError instead of
                falling through to meep. Intended for notebooks/docs
                rendering that must reuse a script's already-populated
                cache rather than silently kicking off a fresh (possibly
                multi-hour) run. False (default) is the normal behavior.
            **geometry_kwargs: forwarded to geometry.section_prisms/width_at/
                center_at (e.g. port_idx for a branching output geometry).

        Returns:
            Complex ndarray [len(wavelengths), n_modes_out, n_modes_in], or
            [len(wavelengths), len(measure_margins), n_modes_out, n_modes_in]
            when measure_margins is given.
        """
        if (y0 is None) != (y1 is None):
            raise ValueError("solve_section: y0 and y1 must both be given or both omitted.")
        for _label, _port in (("excite_port", excite_port), ("measure_port", measure_port)):
            if _port is not None and _port.axis == "y" and y0 is None:
                raise ValueError(
                    f"solve_section: {_label}.axis='y' requires y0 and y1 to be given "
                    "(there is no y-extent to place a y-normal plane in otherwise)."
                )
        if measure_margin is not None and measure_port is None:
            raise ValueError(
                "solve_section: measure_margin is only meaningful together with "
                "measure_port (it overrides that port's own along-axis offset) — "
                "got measure_margin but measure_port=None."
            )
        if measure_margins is not None:
            if measure_port is None:
                raise ValueError(
                    "solve_section: measure_margins is only meaningful together with "
                    "measure_port (each entry overrides that port's own along-axis "
                    "offset for one of its many monitors) — got measure_margins but "
                    "measure_port=None."
                )
            if measure_margin is not None:
                raise ValueError(
                    "solve_section: measure_margin and measure_margins are mutually "
                    "exclusive — pass one along-axis offset (measure_margin) for a "
                    "single monitor, or a list (measure_margins) for many monitors "
                    "in one simulation, not both."
                )
            if len(measure_margins) == 0:
                raise ValueError("solve_section: measure_margins must not be empty.")

        import meep as mp

        wavelengths = np.asarray(wavelengths)

        # Position/size math (which side excites vs. measures, source/monitor
        # apertures, the source_aperture_at()/half_domain_aperture priority
        # order) lives in section_fdtd_layout() — shared with the plotting
        # code in beef/plot/geometry.py's _draw_section_fdtd_markers(), so
        # the two can't drift apart (see that function's docstring). Returned
        # in the section's local frame (0..L); shifted into this method's
        # own meep-centered frame below. Computed before the cache-config
        # block (moved up from its previous spot further down) so its
        # fully-resolved positions can be fingerprinted alongside the
        # monitor_center/monitor_size/margin *ingredients* that derive
        # them — a future change to that derivation logic then shows up as
        # a config diff instead of a silent stale-cache hit (Solver
        # recommendation #4).
        # Skip section_fdtd_layout() entirely when BOTH sides are PortSpec-
        # overridden: nothing downstream would use its geometry-derived
        # defaults anyway, and calling it regardless would be actively wrong
        # for a geometry whose entire section IS a curved/mixed-orientation
        # feature with no "uniform region eps inside the boundary" the way
        # every straight-flanked geometry (DirectionalCoupler, MMI, SBend,
        # Taper) has — section_fdtd_layout() unconditionally samples
        # geometry.width_at()/center_at() near x0/x1 to compute values that
        # would be discarded once both PortSpecs are given. E.g.
        # RingResonatorBendSection's entire [0, bend_radius] span is its own
        # bend, so that eps-sampling always lands inside its own _check_cut()
        # guard — found by tests/test_ring_orientation.py's real-meep check,
        # which failed here before this fix.
        both_overridden = excite_port is not None and measure_port is not None
        layout = None if both_overridden else section_fdtd_layout(
            geometry, x0, x1, stub_left, stub_right,
            n_modes_in, n_modes_out, self.sim_width,
            self.source_margin, self.monitor_margin, self.half_domain_aperture,
            geometry_kwargs, monitor_center=monitor_center, monitor_size=monitor_size,
        )

        # A request for fewer modes than a prior call already computed can
        # reuse that larger result by slicing its mode axes, IF the two
        # calls are otherwise physically identical — same excited/measured
        # side (layout.reverse) and a pinned monitor_size (the default
        # monitor width depends on whether n_measure > 1, so without an
        # explicit monitor_size the *simulated* monitor geometry itself
        # would differ across mode counts, not just the number of bands
        # decoded from it — see lookup_superset()'s docstring). family_config
        # is exactly build_run_config()'s usual fingerprint with the mode
        # counts themselves removed and layout.reverse added in their place
        # (Solver recommendation #2).
        # override/reverse: whenever either PortSpec is given, the excited/
        # measured side is taken literally from n_modes_in/n_modes_out
        # rather than section_fdtd_layout()'s own `n_out <= n_in` tie-break
        # — that tie-break silently picks the wrong physical side whenever
        # n_in == n_out (verified: the standalone-bend case, n_in=n_out=1,
        # excites from the wrong end under the default rule). `reverse`
        # here is the value actually used below, not layout.reverse's raw
        # (possibly bypassed) value — fingerprinted as such.
        override = excite_port is not None or measure_port is not None
        reverse = False if override else layout.reverse

        def _port_key(p: Optional[PortSpec]) -> Optional[dict]:
            return asdict(p) if p is not None else None

        if measure_diagonal_only and n_modes_in != n_modes_out:
            raise ValueError(
                f"solve_section(measure_diagonal_only=True) needs a square "
                f"problem, got n_modes_in={n_modes_in}, n_modes_out={n_modes_out}. "
                f"T_ii is only defined when the excited and measured bases are "
                f"the same size."
            )

        # A diagonal-only result cannot serve a later full-matrix request for
        # fewer modes, so it must never be registered as (or matched against)
        # a mode-count family superset. The useful direction -- a full matrix
        # serving a diagonal request -- is handled explicitly below instead.
        family_config = None
        if (save_dir and monitor_size is not None and measure_margins is None
                and not measure_diagonal_only):
            family_config = build_run_config(
                self, "solve_section",
                dict(x0=x0, x1=x1, wavelengths=wavelengths,
                     stub_left=stub_left, stub_right=stub_right,
                     monitor_center=monitor_center, monitor_size=monitor_size,
                     y0=y0, y1=y1, excite_port=_port_key(excite_port),
                     measure_port=_port_key(measure_port), measure_margin=measure_margin,
                     geometry_kwargs=geometry_kwargs, reverse=reverse,
                     **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor),
                     **resolved_material_config(
                         self.material_stack, wavelengths,
                         self.material_stack.resolve_dispersive(self.dispersive))),
                geometry,
            )

        if save_dir:
            def _run_config(diagonal_only: bool) -> dict:
                # `measure_diagonal_only` is added to the fingerprint ONLY when
                # True. A key present on every config would change the hash of
                # every config, invalidating the entire existing save_dir --
                # 566 entries here, including the 30.6 h of MMI 2x2 calibration
                # this flag exists to avoid ever paying again. Keeping it
                # additive means a False call fingerprints exactly as it did
                # before the flag existed, so every prior entry still hits and
                # `_run_config(False)` below is a real lookup of the full-matrix
                # twin rather than of something that was never written.
                call = dict(
                    x0=x0, x1=x1, wavelengths=wavelengths,
                    n_modes_in=n_modes_in, n_modes_out=n_modes_out,
                    stub_left=stub_left, stub_right=stub_right,
                    monitor_center=monitor_center, monitor_size=monitor_size,
                    y0=y0, y1=y1, excite_port=_port_key(excite_port),
                    measure_port=_port_key(measure_port), measure_margin=measure_margin,
                    measure_margins=(list(measure_margins) if measure_margins is not None else None),
                    geometry_kwargs=geometry_kwargs,
                    resolved_layout=(asdict(layout) if layout is not None else None),
                    **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor),
                    **resolved_material_config(
                        self.material_stack, wavelengths,
                        self.material_stack.resolve_dispersive(self.dispersive)),
                )
                if diagonal_only:
                    call["measure_diagonal_only"] = True
                return build_run_config(self, "solve_section", call, geometry)

            config = _run_config(measure_diagonal_only)
            hit = lookup(save_dir, config)
            if hit is not None:
                logger.info("solve_section: cache hit at %s", hit)
                return load_result(hit)
            if measure_diagonal_only:
                # A full matrix already contains the diagonal this call wants.
                # Reusing it is what keeps switching a calibration over to
                # diagonal-only from throwing away the cache an earlier full
                # run paid for -- e.g. the MMI 2x2's own 30.6 h in job
                # 11872191, whose three cached matrices still serve the
                # calibration unchanged after this flag was introduced.
                full_hit = lookup(save_dir, _run_config(False))
                if full_hit is not None:
                    logger.info(
                        "solve_section: diagonal served from the full-matrix "
                        "cache entry at %s", full_hit,
                    )
                    return _diagonal_of(load_result(full_hit))
            if family_config is not None:
                superset = lookup_superset(save_dir, family_config, (n_modes_out, n_modes_in))
                if superset is not None:
                    folder, _stored_dims = superset
                    logger.info("solve_section: superset cache hit at %s (sliced)", folder)
                    return load_result(folder)[:, :n_modes_out, :n_modes_in]
        if cache_only:
            raise CacheMissError(
                f"solve_section(): cache_only=True but no cached run found for "
                f"this config (save_dir={save_dir!r}). Run once with "
                f"cache_only=False to populate the cache."
            )

        wl_center = float(np.mean(wavelengths))
        si_med, sio2_med = self._media(wavelengths)

        # Draw the waveguide dpml further than it functionally needs to reach
        # (on *both* sides, always — every section is its own independent
        # little simulation with its own PML on both sides, interior
        # HybridSolver seams included) so it continues through the full PML
        # thickness instead of stopping dead at the PML's inner edge, which
        # would present the absorbing layer with an abrupt core->cladding
        # step and risk spurious back-reflection. This only pads what's
        # drawn — lo/hi/x_center/Sx below (and hence every monitor/source
        # position) still key off the true stub_left/stub_right, matching
        # solve()'s own geometry.to_meep_prisms(stub_length=waveguide_length
        # + dpml) precedent.
        prisms = geometry.section_prisms(
            x0, x1, si_med,
            stub_left=stub_left + self.dpml, stub_right=stub_right + self.dpml,
            wg_height=self._wg_height,
            **geometry_kwargs,
        )

        span = x1 - x0
        Sx = span + stub_left + stub_right + 2 * self.dpml
        if y0 is not None:
            Sy = (y1 - y0) + 2 * self.dpml
            y_center = (y0 + y1) / 2.0
        else:
            Sy = self.sim_width + 2 * self.dpml
            y_center = 0.0

        # Centre the meep cell on the *padded* domain [x0-stub_left, x1+stub_right],
        # not the bare section [x0, x1] — meep always centres its cell_size box at
        # the origin, so when stub_left != stub_right (every true end-of-chain
        # section, since only one side ever gets the access-waveguide stub) using
        # the section's own midpoint (x0+x1)/2 as the shift leaves the geometry
        # off-centre within the cell by (stub_right-stub_left)/2 — the PML would
        # sit asymmetrically relative to the drawn stub instead of exactly dpml
        # beyond it on both sides. lo/hi (and hence x_center) already match the
        # frame section_gds()/section_prisms() build their geometry in.
        lo, hi = x0 - stub_left, x1 + stub_right
        x_center = (lo + hi) / 2.0

        fcen = 1.0 / wl_center
        fwidth = source_fwidth(wavelengths, self.source_bandwidth_factor)
        mon_freqs, wl_order = monitor_freq_grid(wavelengths)
        check_source_coverage(fcen, fwidth, wavelengths)

        # shift prisms from device frame into the meep-centered frame,
        # preserving each prism's own z placement/extrusion (vertex z and
        # height carry the 2-D-infinite vs 3-D-core convention, see
        # polygons_to_prisms in beef/geometry/base.py). y_center is 0 (a
        # no-op subtraction) whenever y0/y1 are not given, i.e. for every
        # caller before PortSpec existed.
        shifted = [
            mp.Prism(
                vertices=[mp.Vector3(v.x - x_center, v.y - y_center, v.z) for v in p.vertices],
                height=p.height, material=p.material,
            )
            for p in prisms
        ]

        # NaN, not 0, under measure_diagonal_only: the untouched entries are
        # never measured by this run, and must not read as a measured zero.
        _fill = np.nan + 0j if measure_diagonal_only else 0
        if measure_margins is not None:
            T = np.full((len(wavelengths), len(measure_margins), n_modes_out, n_modes_in),
                        _fill, dtype=complex)
        else:
            T = np.full((len(wavelengths), n_modes_out, n_modes_in), _fill, dtype=complex)

        # source_margin/monitor_margin are documented (see class docstring) as
        # offsets from the PML *inner* edge — i.e. from lo/hi (device frame),
        # not from x0/x1 directly; those only coincide when stub_left/right is
        # 0 (a pure interior cut). layout's positions are in the local frame
        # (0 at lo, L at hi) and shift into this method's meep-centered frame
        # by -L/2 == lo - x_center.
        #
        # excite_port/measure_port bypass layout's x-normal-only fields
        # entirely for the side they cover: each resolves its own axis,
        # along-axis position (source and norm-monitor separately, offset
        # by source_margin/monitor_margin from whichever edge `sign` picks),
        # transverse center, aperture width, launch kpoint, and
        # get_eigenmode_coefficients() band index — independently per side,
        # since (verified against ring_resonator_hybrid_cascaded.py's
        # run_left_fdtd()) the excite/norm side's index and the measure
        # side's index are genuinely different quantities, not one shared
        # value. A None side keeps using layout's geometry-derived value,
        # still shifted by -y_center (a no-op unless the OTHER side forced
        # y0/y1 to be non-default).
        n_excite = n_modes_in if override else layout.n_excite
        n_measure = n_modes_out if override else layout.n_measure

        if excite_port is not None:
            p = excite_port
            ex_axis = p.axis
            ex_width = p.width
            lo_ax, hi_ax = (lo, hi) if ex_axis == "x" else (y0, y1)
            along_shift, transverse_shift = (x_center, y_center) if ex_axis == "x" else (y_center, x_center)
            ex_along_src = (lo_ax if p.sign == 1 else hi_ax) + p.sign * self.source_margin - along_shift
            ex_along_norm = (lo_ax if p.sign == 1 else hi_ax) + p.sign * self.monitor_margin - along_shift
            ex_transverse = p.center - transverse_shift
            ex_kpoint = (float(p.sign), 0.0, 0.0) if ex_axis == "x" else (0.0, float(p.sign), 0.0)
            ex_alpha_idx = p.alpha_idx
            ex_sign_dbg = p.sign
        else:
            ex_axis = "x"
            ex_sign_dbg = 0
            ex_width = layout.src_width
            ex_transverse = layout.src_c - y_center
            ex_along_src = layout.src_x + lo - x_center
            ex_along_norm = layout.norm_x + lo - x_center
            ex_kpoint = layout.kpoint
            ex_alpha_idx = layout.alpha_idx

        if measure_port is not None:
            p = measure_port
            mea_axis = p.axis
            mea_width = p.width
            lo_ax, hi_ax = (lo, hi) if mea_axis == "x" else (y0, y1)
            along_shift, transverse_shift = (x_center, y_center) if mea_axis == "x" else (y_center, x_center)
            # Measure ports use the OPPOSITE edge convention from excite
            # ports for the same sign, unlike excite_port's edge choice
            # (which is physically forced: a source launching +axis must
            # sit at the LOW edge or it can't illuminate the domain at
            # all). A monitor has no such constraint — it passively reads
            # whichever direction is present. sign=+1 (reading +axis-
            # forward light, alpha_idx=0) naturally sits at the HIGH edge
            # (downstream, where forward-propagating light is captured
            # after crossing the section, e.g. an interface plane at the
            # far end of a bent section whose input was x-normal); sign=-1
            # at the LOW edge. This is a physically-motivated choice, not
            # verified against the reference scripts — they never override
            # the measure side's position via an explicit PortSpec, only
            # excite_port's LOW/HIGH convention is directly verified
            # (ring_resonator_hybrid_cascaded.py's run_left_fdtd()).
            # measure_margin (None by default) lets this side's offset be
            # set independently of self.monitor_margin, which the excite
            # side's own ex_along_norm above always uses regardless —
            # without this override the two sides would be forced to share
            # one offset (see solve_section()'s own docstring).
            # measure_margins (a list) instead builds one along-axis
            # position PER entry — mea_along becomes a list too, one flux
            # monitor per position added into the SAME simulation below.
            if measure_margins is not None:
                mea_along = [
                    (hi_ax if p.sign == 1 else lo_ax) - p.sign * mm - along_shift
                    for mm in measure_margins
                ]
            else:
                _mea_margin = self.monitor_margin if measure_margin is None else measure_margin
                mea_along = (hi_ax if p.sign == 1 else lo_ax) - p.sign * _mea_margin - along_shift
            mea_transverse = p.center - transverse_shift
            mea_alpha_idx = p.alpha_idx
        else:
            mea_axis = "x"
            mea_width = layout.mea_w
            mea_transverse = layout.mea_c - y_center
            mea_along = layout.mea_x + lo - x_center
            mea_alpha_idx = layout.alpha_idx

        # --- Aperture-vs-PML guard -------------------------------------
        # A port's aperture spans the axis TRANSVERSE to its own normal, and
        # nothing clamps an explicit PortSpec width the way
        # section_fdtd_layout() clamps its geometry-derived one
        # (`min(src_aperture, sim_width - 2*abs(src_c))`). meep silently
        # clips a monitor that runs past the cell, so an overrunning aperture
        # reads as loss: ring_resonator_validation.md's blocker-0 work found
        # a y-normal port sitting exactly at the PML inner edge with a
        # width=3.0 aperture, which cost 27% of the transmitted power
        # (|T| 0.711 vs the correct 0.978) while still looking physically
        # plausible, and passed a `> 0.3` sanity assertion for days. Raise
        # instead -- an aperture in PML is always a caller error.
        def _check_aperture(kind: str, axis: str, transverse: float, width: float) -> None:
            span = Sy if axis == "x" else Sx      # aperture runs along the OTHER axis
            inner = span / 2.0 - self.dpml        # PML inner edge, meep-centred frame
            reach = abs(transverse) + width / 2.0
            if reach > inner + 1e-9:
                other = "y" if axis == "x" else "x"
                raise ValueError(
                    f"{kind} aperture crosses the PML: axis={axis!r} port spans "
                    f"{other}=[{transverse - width/2:.4f}, {transverse + width/2:.4f}] "
                    f"(centre {transverse:.4f}, width {width:.4f}) but the PML inner "
                    f"edge is at {other}=+/-{inner:.4f} (cell {span:.4f}, dpml "
                    f"{self.dpml}). meep would silently CLIP it and the missing "
                    f"power would read as loss. Widen the domain along {other} "
                    f"(stub_left/stub_right for x, y0/y1 for y) or narrow the "
                    f"aperture to <= {2 * (inner - abs(transverse)):.4f}."
                )

        _check_aperture("excite_port", ex_axis, ex_transverse, ex_width)
        _check_aperture("measure_port", mea_axis, mea_transverse, mea_width)

        for k in range(1, n_excite + 1):
            src = mp.EigenModeSource(
                mp.GaussianSource(fcen, fwidth=fwidth),
                center=mp.Vector3(*_port_point(ex_axis, ex_along_src, ex_transverse)),
                size=mp.Vector3(*_port_size(ex_axis, ex_width, self._mon_sz)),
                eig_band=k,
                eig_kpoint=mp.Vector3(*ex_kpoint),
                eig_match_freq=True,
            )
            sim = mp.Simulation(
                resolution=self.res,
                cell_size=mp.Vector3(Sx, Sy, self._cell_sz),
                boundary_layers=[mp.PML(self.dpml)],
                geometry=shifted,
                sources=[src],
                default_material=sio2_med,
                eps_averaging=self.eps_averaging,
                subpixel_tol=self.subpixel_tol,
                subpixel_maxeval=self.subpixel_maxeval,
            )
            norm_flux = sim.add_mode_monitor(
                mon_freqs,
                mp.FluxRegion(center=mp.Vector3(*_port_point(ex_axis, ex_along_norm, ex_transverse)),
                              size=mp.Vector3(*_port_size(ex_axis, ex_width, self._mon_sz))),
            )
            if measure_margins is not None:
                # One flux monitor per swept position, all in this SAME
                # simulation -- sim.run() below times-steps once regardless
                # of monitor count (stop_when_dft_decayed() watches the
                # whole field, not any one monitor), so N positions cost
                # about the same wall-clock time as 1.
                mea_fluxes = [
                    sim.add_mode_monitor(
                        mon_freqs,
                        mp.FluxRegion(center=mp.Vector3(*_port_point(mea_axis, along, mea_transverse)),
                                      size=mp.Vector3(*_port_size(mea_axis, mea_width, self._mon_sz))),
                    )
                    for along in mea_along
                ]
            else:
                mea_flux = sim.add_mode_monitor(
                    mon_freqs,
                    mp.FluxRegion(center=mp.Vector3(*_port_point(mea_axis, mea_along, mea_transverse)),
                                  size=mp.Vector3(*_port_size(mea_axis, mea_width, self._mon_sz))),
                )

            sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol, minimum_run_time=self.minimum_run_time))

            norm_coeff = sim.get_eigenmode_coefficients(norm_flux, [k], eig_parity=mp.NO_PARITY)
            alpha_raw = norm_coeff.alpha[0, :, ex_alpha_idx]
            if ex_axis == "y":
                # Flux-calibrated normalization: get_eigenmode_coefficients()
                # alone undercounts power for an off-axis, backward-launched
                # source — an empirically-found correction, first needed in
                # examples/ring_resonator/python/ring_resonator_hybrid.py's
                # get_norm_calibrated() (magnitude from mp.get_fluxes(),
                # phase from the coefficient). Applied for both signs
                # conservatively; only sign=-1 was ever tested there.
                flux = np.abs(np.asarray(mp.get_fluxes(norm_flux)))
                src_alpha = np.sqrt(flux) * np.exp(1j * np.angle(alpha_raw))
            else:
                # x-axis (in-line) excite: get_eigenmode_coefficients()'s own
                # alpha is normally trustworthy and left untouched here, so
                # every already-validated x-axis section elsewhere (mmi_1x4,
                # mmi_4x4_90hybrid, mmi_2x2, directional_coupler, ...) sees
                # bit-identical output. But it can silently undercount power
                # for a near-degenerate excited band (gap 7's connector,
                # bands 1-2 at Δn_eff=0.0148) — detect that by cross-checking
                # against mp.get_fluxes() (a direct field integral, immune to
                # modal-decomposition ill-conditioning) and fall back to the
                # same flux-calibrated form used above only when they disagree.
                flux = np.abs(np.asarray(mp.get_fluxes(norm_flux)))
                mismatch = np.abs(np.sqrt(flux) - np.abs(alpha_raw)) / np.maximum(np.sqrt(flux), 1e-30)
                if np.any(mismatch > _EXCITE_NORM_MISMATCH_TOL):
                    logger.warning(
                        "solve_section: excite band %d normalization mismatch between "
                        "get_eigenmode_coefficients() and get_fluxes() (max %.1f%% at freq "
                        "index %d) for section x=[%.4f, %.4f] — likely a near-degenerate "
                        "excited mode; falling back to flux-calibrated normalization for "
                        "this band.",
                        k, 100 * float(mismatch.max()), int(mismatch.argmax()), x0, x1,
                    )
                    src_alpha = np.sqrt(flux) * np.exp(1j * np.angle(alpha_raw))
                else:
                    src_alpha = alpha_raw


            if os.environ.get("BEEF_DEBUG_NORM"):
                # Diagnostic for ring_resonator_validation.md's blocker 0: the
                # y-excite path normalizes by sqrt(|get_fluxes|) while every
                # measure side divides by a bare |alpha|. Print both at the
                # SAME plane so the two conventions can be compared directly.
                _f = np.abs(np.asarray(mp.get_fluxes(norm_flux)))
                _a = np.abs(norm_coeff.alpha[0, :, ex_alpha_idx])
                print(f"[norm] ex_axis={ex_axis} sign={ex_sign_dbg} band={k} "
                      f"alpha_idx={ex_alpha_idx}", flush=True)
                for _i in range(len(_a)):
                    print(f"[norm]   f{_i}: |alpha|={_a[_i]:.6g}  sqrt(flux)={np.sqrt(_f[_i]):.6g}  "
                          f"flux={_f[_i]:.6g}  sqrt(flux)/|alpha|={np.sqrt(_f[_i])/_a[_i]:.6g}  "
                          f"flux/|alpha|^2={_f[_i]/_a[_i]**2:.6g}", flush=True)

            # Run k needs only its own column entry when the caller asked for
            # the diagonal -- decomposing the other n_measure-1 bands is the
            # single largest cost in the method (see measure_diagonal_only).
            measure_bands = [k] if measure_diagonal_only else range(1, n_measure + 1)

            if measure_margins is not None:
                for pos_idx, mea_flux_i in enumerate(mea_fluxes):
                    for m in measure_bands:
                        c = sim.get_eigenmode_coefficients(mea_flux_i, [m], eig_parity=mp.NO_PARITY)
                        # [wl_order] lifts the ratio out of meep's ascending-frequency
                        # monitor order back onto `wavelengths` as the caller gave it.
                        coeff = (c.alpha[0, :, mea_alpha_idx] / src_alpha)[wl_order]
                        if reverse:
                            T[:, pos_idx, k - 1, m - 1] = coeff  # excited=out(k), measured=in(m)
                        else:
                            T[:, pos_idx, m - 1, k - 1] = coeff  # excited=in(k), measured=out(m)
            else:
                for m in measure_bands:
                    c = sim.get_eigenmode_coefficients(mea_flux, [m], eig_parity=mp.NO_PARITY)
                    # [wl_order] lifts the ratio out of meep's ascending-frequency
                    # monitor order back onto `wavelengths` as the caller gave it.
                    coeff = (c.alpha[0, :, mea_alpha_idx] / src_alpha)[wl_order]
                    if reverse:
                        T[:, k - 1, m - 1] = coeff  # excited=out(k), measured=in(m)
                    else:
                        T[:, m - 1, k - 1] = coeff  # excited=in(k), measured=out(m)

            sim.reset_meep()

        # Post-hoc duplicate-mode scan of the transfer matrix just computed —
        # fdtd_vs_EME_mode_validation.md item 8. Every other correspondence
        # mechanism in this codebase checks section_mode_fields()'s own
        # separate eigensolve; T here comes from get_eigenmode_coefficients()
        # on DFT-decayed flux, a different code path, never checked before
        # this. Free: T is already computed, this is O(n_excited^2) vector
        # arithmetic on it.
        from ..mode_identity import warn_on_transfer_duplicate

        # The scan compares whole transfer-matrix columns; under
        # measure_diagonal_only there are no columns, only a diagonal.
        if measure_diagonal_only:
            pass
        elif measure_margins is not None:
            for pos_idx, mm in enumerate(measure_margins):
                # hi - mm: this position's actual absolute measurement
                # plane (matches how mea_along itself was derived above),
                # not x1 -- the two only coincide when stub_right == 0.
                warn_on_transfer_duplicate(T[:, pos_idx], logger, x0, hi - mm)
        else:
            warn_on_transfer_duplicate(T, logger, x0, x1)

        if save_dir:
            folder = save_run(save_dir, config, T)
            logger.info("solve_section: cache saved at %s", folder)
            if family_config is not None:
                update_family_pointer(
                    save_dir, family_config, (n_modes_out, n_modes_in), os.path.basename(folder),
                )

        return T

    def solve_section_multi_length(
        self,
        geometry,
        wavelengths: np.ndarray,
        mid: float,
        lengths: Sequence[float],
        n_modes: int,
        monitor_center: Optional[float] = None,
        monitor_size: Optional[float] = None,
        measure_diagonal_only: bool = False,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        **geometry_kwargs,
    ) -> np.ndarray:
        """Modal transfer matrices at MULTIPLE lengths from ONE merged FDTD
        domain — resolved/section_mode_beta_multi_monitor_validation.md's mechanism.

        `section_mode_beta()` fits body beta from `len(lengths)` fully
        independent `solve_section()` calls, each a separate domain build +
        separate decay-to-tolerance wait, even though the underlying
        cross-section never changes (only the propagation distance does).
        This method builds ONE domain sized to `max(lengths)`, with one
        excited-mode source (shared across all lengths) and one interior
        interface monitor per length, and reads all of them off a single
        `sim.run()` — one decay wait instead of `len(lengths)`.

        Reuses `section_fdtd_layout(x0, x1, ...)` completely unchanged
        (called once, with `x0, x1 = mid -+ max(lengths)/2`) rather than a
        new layout function: `layout.src_x`/`norm_x`/`src_c`/`src_width`/
        `mea_c`/`mea_w`/`kpoint`/`alpha_idx` don't depend on length at all,
        and `layout.mea_x` (computed for the `max(lengths)`-length window)
        already gives the boundary case directly. Only `mea_x` needs a
        per-length variant, which lives in `merged_mea_positions()` (shared
        with the plotting code that draws these monitors, so the two cannot
        drift):

            mea_positions[i] = (layout.L - lengths[i]) if layout.reverse
                                else lengths[i]

        which reproduces each length's own single-domain interface-monitor
        distance from the source (`lengths[i] - monitor_margin`) exactly —
        checked algebraically in
        resolved/section_mode_beta_multi_monitor_validation.md step 2, not
        approximated.

        Caching keys on the whole call — `mid`, the full `lengths` list, the
        merged domain's own resolved layout, and the section geometry drawn
        for `[mid ± max(lengths)/2]` — so one entry covers one merged run's
        complete N-length result, which is the granularity the run actually
        has. That deliberately does NOT interoperate with the per-length
        `solve_section()` entries a pre-merge calibration wrote: they are
        different simulations (three snug domains vs one long shared one),
        and step 6 measured that the difference is real, if guard-safe
        (`5.78e-04` rad/µm of domain-length sensitivity). Serving one from
        the other would silently mix two measurements; a merged run whose
        per-length twin is cached therefore re-simulates, once.

        Unlike `solve_section()`: no `PortSpec`/`y0`/`y1` support (not needed
        by `section_mode_beta()`'s calls, which are always plain x-normal
        with `n_modes_in == n_modes_out`), and no mode-count "superset"
        family reuse — narrower than `solve_section()` on purpose, not a
        general replacement for it.

        Args:
            mid: device-frame x-coordinate the lengths are centered on
                (matches `section_mode_beta()`'s own `mid = 0.5*(x0+x1)`).
            lengths: propagation distances to measure (µm) — the same
                values `section_mode_beta()` would otherwise simulate one
                domain each for.
            n_modes: mode count on both sides (matches
                `section_mode_beta()`'s own `n_modes_in=n_modes_out=n_modes`
                calls).
            monitor_center, monitor_size: as `solve_section()` — override
                the (length-independent) interface monitor's y-centre/span.
            **geometry_kwargs: forwarded to `geometry.section_prisms()`/
                `section_fdtd_layout()`.

        Returns:
            Complex ndarray [len(lengths), len(wavelengths), n_modes, n_modes].
        """
        import meep as mp

        wavelengths = np.asarray(wavelengths)
        lengths = [float(L) for L in lengths]
        if not lengths:
            raise ValueError("solve_section_multi_length: lengths must be non-empty.")
        if any(L <= 0.0 for L in lengths):
            raise ValueError(f"solve_section_multi_length: lengths must be > 0, got {lengths}.")

        L_domain = max(lengths)
        x0, x1 = mid - L_domain / 2.0, mid + L_domain / 2.0

        layout = section_fdtd_layout(
            geometry, x0, x1, 0.0, 0.0,
            n_modes, n_modes, self.sim_width,
            self.source_margin, self.monitor_margin, self.half_domain_aperture,
            geometry_kwargs, monitor_center=monitor_center, monitor_size=monitor_size,
        )
        n_excite = layout.n_excite
        n_measure = layout.n_measure

        # Per-length interior monitor positions, local frame (0..L_domain) --
        # see the docstring derivation above. The shortest length's monitor
        # sits CLOSEST to the source (largest mea_positions[i] when reverse,
        # smallest when not) -- how close, for a real device, is exactly
        # what step 6's convergence study (simulation/domain length) is for;
        # not asserted against an invented threshold here.
        mea_positions_local = merged_mea_positions(layout, lengths)

        # Cache lookup sits here, after the monitor positions are resolved, so
        # the fingerprint can carry them: `lengths` alone is the *request*,
        # while mea_positions_local is what actually gets simulated, and
        # Solver recommendation #4 is specifically that a future change to
        # that derivation must not leave an old entry looking identical.
        if save_dir:
            call = dict(
                x0=x0, x1=x1, mid=mid, lengths=lengths, wavelengths=wavelengths,
                n_modes=n_modes, stub_left=0.0, stub_right=0.0,
                monitor_center=monitor_center, monitor_size=monitor_size,
                geometry_kwargs=geometry_kwargs,
                resolved_layout=asdict(layout),
                resolved_mea_positions=list(mea_positions_local),
                **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor),
                **resolved_material_config(
                    self.material_stack, wavelengths,
                    self.material_stack.resolve_dispersive(self.dispersive)),
            )
            if measure_diagonal_only:
                call["measure_diagonal_only"] = True   # additive, as solve_section()
            config = build_run_config(self, "solve_section_multi_length", call, geometry)
            hit = lookup(save_dir, config)
            if hit is not None:
                logger.info("solve_section_multi_length: cache hit at %s", hit)
                return load_result(hit)
        if cache_only:
            raise CacheMissError(
                f"solve_section_multi_length(): cache_only=True but no cached run "
                f"found for this config (save_dir={save_dir!r}). Run once with "
                f"cache_only=False to populate the cache. NOTE: a per-length "
                f"solve_section() cache from a pre-merge calibration does NOT "
                f"satisfy this — the merged domain is a different simulation."
            )

        wl_center = float(np.mean(wavelengths))
        si_med, sio2_med = self._media(wavelengths)

        # stub_left/stub_right always 0 here (section_mode_beta() never
        # passes them) -- see solve_section()'s own comment for why dpml is
        # still added on both sides of the drawn waveguide regardless.
        prisms = geometry.section_prisms(
            x0, x1, si_med,
            stub_left=self.dpml, stub_right=self.dpml,
            wg_height=self._wg_height,
            **geometry_kwargs,
        )

        span = x1 - x0
        Sx = span + 2 * self.dpml
        Sy = self.sim_width + 2 * self.dpml

        # lo == x0, hi == x1 (stub_left=stub_right=0), so x_center == mid
        # exactly -- the meep cell is centred on the same point the caller's
        # lengths are centred on, by construction.
        lo, hi = x0, x1
        x_center = (lo + hi) / 2.0

        fcen = 1.0 / wl_center
        fwidth = source_fwidth(wavelengths, self.source_bandwidth_factor)
        mon_freqs, wl_order = monitor_freq_grid(wavelengths)
        check_source_coverage(fcen, fwidth, wavelengths)

        shifted = [
            mp.Prism(
                vertices=[mp.Vector3(v.x - x_center, v.y, v.z) for v in p.vertices],
                height=p.height, material=p.material,
            )
            for p in prisms
        ]

        T = np.full((len(lengths), len(wavelengths), n_modes, n_modes),
                    np.nan + 0j if measure_diagonal_only else 0, dtype=complex)

        ex_width = layout.src_width
        ex_transverse = layout.src_c
        ex_along_src = layout.src_x + lo - x_center
        ex_along_norm = layout.norm_x + lo - x_center
        ex_kpoint = layout.kpoint
        ex_alpha_idx = layout.alpha_idx

        mea_width = layout.mea_w
        mea_transverse = layout.mea_c
        mea_alpha_idx = layout.alpha_idx
        mea_along_list = [pos + lo - x_center for pos in mea_positions_local]

        # Same PML-crossing guard as solve_section() -- see that method's
        # own docstring for the ring_resonator_validation.md blocker-0
        # history this is guarding against. Transverse aperture is
        # length-independent (checked once for excite, once for measure).
        def _check_aperture(kind: str, transverse: float, width: float) -> None:
            inner = Sy / 2.0 - self.dpml
            reach = abs(transverse) + width / 2.0
            if reach > inner + 1e-9:
                raise ValueError(
                    f"{kind} aperture crosses the PML: spans y="
                    f"[{transverse - width/2:.4f}, {transverse + width/2:.4f}] "
                    f"(centre {transverse:.4f}, width {width:.4f}) but the PML "
                    f"inner edge is at y=+/-{inner:.4f} (cell {Sy:.4f}, dpml "
                    f"{self.dpml})."
                )

        _check_aperture("excite", ex_transverse, ex_width)
        _check_aperture("measure", mea_transverse, mea_width)

        for k in range(1, n_excite + 1):
            src = mp.EigenModeSource(
                mp.GaussianSource(fcen, fwidth=fwidth),
                center=mp.Vector3(*_port_point("x", ex_along_src, ex_transverse)),
                size=mp.Vector3(*_port_size("x", ex_width, self._mon_sz)),
                eig_band=k,
                eig_kpoint=mp.Vector3(*ex_kpoint),
                eig_match_freq=True,
            )
            sim = mp.Simulation(
                resolution=self.res,
                cell_size=mp.Vector3(Sx, Sy, self._cell_sz),
                boundary_layers=[mp.PML(self.dpml)],
                geometry=shifted,
                sources=[src],
                default_material=sio2_med,
                eps_averaging=self.eps_averaging,
                subpixel_tol=self.subpixel_tol,
                subpixel_maxeval=self.subpixel_maxeval,
            )
            norm_flux = sim.add_mode_monitor(
                mon_freqs,
                mp.FluxRegion(center=mp.Vector3(*_port_point("x", ex_along_norm, ex_transverse)),
                              size=mp.Vector3(*_port_size("x", ex_width, self._mon_sz))),
            )
            mea_fluxes = [
                sim.add_mode_monitor(
                    mon_freqs,
                    mp.FluxRegion(center=mp.Vector3(*_port_point("x", mea_along, mea_transverse)),
                                  size=mp.Vector3(*_port_size("x", mea_width, self._mon_sz))),
                )
                for mea_along in mea_along_list
            ]

            sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol, minimum_run_time=self.minimum_run_time))

            norm_coeff = sim.get_eigenmode_coefficients(norm_flux, [k], eig_parity=mp.NO_PARITY)
            alpha_raw = norm_coeff.alpha[0, :, ex_alpha_idx]
            # x-axis excite only (no PortSpec/y-axis support here) -- same
            # flux-calibrated fallback as solve_section()'s own x-axis
            # branch, guarding against a near-degenerate excited band
            # silently undercounting power (modal_multiplexer_validation.md
            # gap 7).
            flux = np.abs(np.asarray(mp.get_fluxes(norm_flux)))
            mismatch = np.abs(np.sqrt(flux) - np.abs(alpha_raw)) / np.maximum(np.sqrt(flux), 1e-30)
            if np.any(mismatch > _EXCITE_NORM_MISMATCH_TOL):
                logger.warning(
                    "solve_section_multi_length: excite band %d normalization "
                    "mismatch between get_eigenmode_coefficients() and "
                    "get_fluxes() (max %.1f%% at freq index %d) for domain "
                    "x=[%.4f, %.4f] -- likely a near-degenerate excited mode; "
                    "falling back to flux-calibrated normalization for this band.",
                    k, 100 * float(mismatch.max()), int(mismatch.argmax()), x0, x1,
                )
                src_alpha = np.sqrt(flux) * np.exp(1j * np.angle(alpha_raw))
            else:
                src_alpha = alpha_raw

            # As solve_section(): the caller wanting only T_ii makes the other
            # n_measure-1 decompositions per monitor pure waste, and here that
            # waste is multiplied by the monitor count.
            measure_bands = [k] if measure_diagonal_only else range(1, n_measure + 1)

            for i, mea_flux in enumerate(mea_fluxes):
                for m in measure_bands:
                    c = sim.get_eigenmode_coefficients(mea_flux, [m], eig_parity=mp.NO_PARITY)
                    coeff = (c.alpha[0, :, mea_alpha_idx] / src_alpha)[wl_order]
                    if layout.reverse:
                        T[i, :, k - 1, m - 1] = coeff
                    else:
                        T[i, :, m - 1, k - 1] = coeff

            sim.reset_meep()

        if save_dir:
            folder = save_run(save_dir, config, T)
            logger.info("solve_section_multi_length: cache saved at %s", folder)

        from ..mode_identity import warn_on_transfer_duplicate

        for i, L_i in enumerate(lengths):
            warn_on_transfer_duplicate(T[i], logger, x0, x1)

        return T
