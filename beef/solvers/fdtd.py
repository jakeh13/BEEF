"""FDTDSolver: meep-based FDTD solver for MMI and related photonic devices."""
from __future__ import annotations

import logging
import os
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple

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

# --- section_mode_beta() safety thresholds ---------------------------------
#
# Two INDEPENDENT failure modes, because neither guard catches the other:
#
#   1. Seed lands near a branch boundary -> the per-length wrap counts come out
#      inconsistent, the points stop lying on a line, and the RESIDUAL blows up
#      (2.5 rad in testing). Caught by _BETA_FIT_RESIDUAL_MAX.
#   2. Seed is a WHOLE branch off -> every length agrees on that wrong branch,
#      so the fit is perfectly clean (residual 1e-14) at a β 33% wrong, and the
#      margin-to-seed is 0.000. Nothing in the fitted numbers can detect this.
#      The only defence is a precondition on the LENGTHS: keep the branch
#      spacing 2π/ΔL large compared with any physically possible FDTD-vs-MPB
#      disagreement, so a whole-branch error cannot happen in the first place.
#
# Minimum branch spacing as a fraction of β, for guard 2. The measured
# disagreement is ~2%; requiring the branches to be >=25% of β apart leaves an
# order of magnitude. At the coupler's β = 9.6 rad/µm this caps the length
# spacing at 2.6 µm, which the item 6 measurement's 2.0 µm satisfies.
_BETA_MIN_BRANCH_FRACTION = 0.25

# Max phase-fit residual (rad) for guard 1. Clean fits give 1e-14 synthetically
# and 1.4e-5 on the real res-40 data; a boundary-straddling seed gives 2.5. Also
# catches a section that is not actually uniform, where the phase is genuinely
# not linear in L.
_BETA_FIT_RESIDUAL_MAX = 0.1

# Reported diagnostic only (see beta_from_phases): how far the answer moved from
# the seed, as a fraction of the half-distance to the next branch. The coupler
# sits at 0.13. Logged rather than enforced, because guard 2 above is what makes
# the branch choice sound.
_BETA_BRANCH_MARGIN_MAX = 0.5


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


def beta_from_phases(lengths, phases, seed_beta) -> Tuple[float, float, float]:
    """Absolute phase rate (rad/µm) for ONE mode, from ``arg(T_ii)`` vs length.

    Unlike `dbeta_from_phases()` this is a *large* rate — β ≈ 9.4 rad/µm for a
    silicon strip — so the phase wraps many times between lengths and
    ``np.unwrap`` is useless: at the coupler's 2 µm spacing, ``β·ΔL`` is 2.9993
    full turns, and all that survives is the 0.002 rad remainder. Getting this
    wrong is not subtle but it is silent — a naive fit here returned Δβ = −0.331
    against a true 0.08525 during the item 6 analysis.

    So the branch is chosen against `seed_beta`, normally MPB's answer from
    `FDTDSolver.section_mode_neff()`. That is safe rather than circular: with
    equally spaced lengths the candidate slopes are ``2π/ΔL`` apart — 3.14 rad/µm
    at ΔL = 2 µm, a 33% step — while the FDTD/MPB disagreement being measured is
    of order 2%. The seed picks which turn; it does not set the answer within
    that turn, and `margin` below reports how much room was left.

    Args:
        lengths: section lengths (µm). Need not be equally spaced, but the
            wrap-safety margin is governed by the largest gap.
        phases: ``arg(T_ii)`` at each length, each wrapped into (−π, π].
        seed_beta: prior for branch selection (rad/µm).

    Returns:
        ``(beta, resid, margin)``. `resid` is the max absolute deviation from
        the fitted line, in radians — small only if the branch is right AND the
        section is uniform. `margin` is ``|beta - seed_beta|`` as a fraction of
        the half-distance to the neighbouring branch: ≪1 is safe, ≥1 means the
        seed could not distinguish the branches and the result is not
        trustworthy.
    """
    lengths = np.asarray(lengths, dtype=float)
    phases = np.asarray(phases, dtype=float)
    seed_beta = float(seed_beta)

    turns = np.round((seed_beta * lengths - phases) / (2.0 * np.pi))
    unwrapped = phases + 2.0 * np.pi * turns
    A = np.vstack([lengths, np.ones_like(lengths)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, unwrapped, rcond=None)
    resid = float(np.max(np.abs(unwrapped - (slope * lengths + intercept))))

    # Neighbouring branches differ by 2π/ΔL in slope, so the seed is decisive
    # only within ±π/ΔL of the truth. Use the largest gap: it is the one that
    # wraps first.
    max_gap = float(np.max(np.diff(np.sort(lengths)))) if lengths.size > 1 else np.inf
    half_branch = np.pi / max_gap if np.isfinite(max_gap) and max_gap > 0 else np.inf
    margin = float(abs(slope - seed_beta) / half_branch) if np.isfinite(half_branch) else 0.0
    return float(slope), resid, margin


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

    if monitor_center is not None:
        mea_c = monitor_center
    mea_w = (sim_width if n_measure > 1 else mea_w0) if monitor_size is None else monitor_size

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

    return SectionFDTDLayout(
        L=L, n_excite=n_excite, n_measure=n_measure, reverse=reverse,
        src_x=src_x, norm_x=norm_x, mea_x=mea_x,
        src_c=src_c, src_width=src_width, mea_c=mea_c, mea_w=mea_w,
        kpoint=kpoint, alpha_idx=alpha_idx,
    )


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
            None → use sim_width for each (full cross-section).
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
        decay_tol: DFT field decay tolerance for early-stop criterion.
            Lower = more accurate but slower; 1e-9 is faster, 1e-11 is tight.
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
        decay_tol: float = 1e-11,
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
        self.decay_tol = decay_tol
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

    def _build_domain(self, geometry: DeviceGeometry, wavelengths: np.ndarray):
        """Return (Sx, Sy, wl_center, si_medium, sio2_medium).

        sim_width/sim_height are inner dimensions (excluding PML).
        Sx and Sy include 2*dpml on each axis (as does the z axis in 3-D —
        see _cell_sz).
        """
        wl_center = float(np.mean(wavelengths))
        Sx = geometry.total_length + 2 * self.waveguide_length + 2 * self.dpml
        Sy = self.sim_width + 2 * self.dpml   # full cell: inner + PML on both sides
        si_med = self.material_stack.si_medium(wl_center)
        sio2_med = self.material_stack.sio2_medium(wl_center)
        return Sx, Sy, wl_center, si_med, sio2_med

    def solve(
        self,
        geometry: DeviceGeometry,
        wavelengths: np.ndarray,
        port: str = "top",
        input_mode: int = 1,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
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
                  (1×2 MMI), since input_ys has only one entry there.
            input_mode: 1-based waveguide mode index to excite and normalize against.
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache — see
                solve_section()'s save_dir docstring for the full semantics
                (identical here, fingerprinted via
                geometry.polygons(stub_length=...) instead of section_gds()).
            cache_only: if True, never run the simulation — return the
                cached result or raise CacheMissError. See solve_section()'s
                cache_only docstring for the full semantics (identical here).
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
        if save_dir:
            stub_length = self.waveguide_length + self.dpml
            config = build_whole_device_config(
                self, "solve",
                dict(port=port, input_mode=input_mode, wavelengths=wavelengths, kwargs=kwargs,
                     **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor),
                     **resolved_material_config(self.material_stack, wavelengths)),
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
        input_idx = len(input_ys) - 1 if port == "top" else 0
        src_yo = input_ys[input_idx]
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
            src_aperture = _aperture_at(src_x, **kwargs)
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

        # One output monitor per output waveguide
        out_widths = (
            self.output_monitor_widths
            if self.output_monitor_widths is not None
            else [self.sim_width] * len(self.output_monitor_centers)
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

        sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol))

        # Normalise by forward input-mode amplitude so |coeff|² = physical T
        in_coeff = sim.get_eigenmode_coefficients(
            in_flux, [input_mode], eig_parity=mp.NO_PARITY
        )
        src_alpha = in_coeff.alpha[0, :, 0]  # [nfreq] forward amplitude, monitor order

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

        Cheap: ``get_eigenmode()`` runs MPB on the initialized grid and returns
        without any time-stepping, so this costs a ``set_epsilon`` plus one
        eigensolve per mode per wavelength — seconds, not the minutes a
        ``solve_section()`` run needs.

        n_eff is ``k/f``: meep reports the Bloch wavevector in units of 2π/a
        and the frequency in c/a, so their ratio is the index directly, with no
        2π anywhere. Pinned against a bare 450×220 nm Si strip, where the
        answer must fall between the cladding and core indices — the competing
        normalizations land at 0.37 and 14.5 and are not close calls.

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

        Returns:
            Real ndarray ``[len(wavelengths), n_modes]`` of effective indices,
            in the caller's wavelength order.
        """
        wavelengths = np.asarray(wavelengths, dtype=float)
        wl_center = float(np.mean(wavelengths))
        sim, where, mea_c, mea_w = self._build_section_mode_sim(
            geometry, wl_center, x0, x1, monitor_center, monitor_size, pad=pad, **geometry_kwargs)

        neff = np.zeros((len(wavelengths), n_modes))
        for i, wl in enumerate(wavelengths):
            f = 1.0 / float(wl)
            # kpoint only seeds the solve; match_frequency then iterates to the
            # true k. Any guided mode sits between the cladding and core
            # indices, so the midpoint of that range is a seed that cannot be
            # badly wrong for any waveguide this solver handles.
            # _iter_distinct_bands() guards against get_eigenmode() returning
            # the same eigenvector for two different band numbers (see its
            # docstring / gap 6) -- band as yielded may skip some meep
            # band_nums, but always fills n_modes distinct logical positions.
            for band, em in enumerate(
                self._iter_distinct_bands(sim, where, f, mea_c, mea_w, n_modes), start=1
            ):
                neff[i, band - 1] = float(np.real(em.k.x)) / f

        sim.reset_meep()
        return neff

    def _build_section_mode_sim(
        self, geometry, wl_center: float, x0: float, x1: float,
        monitor_center: Optional[float], monitor_size: Optional[float],
        pad: Optional[float] = None,
        **geometry_kwargs,
    ):
        """Shared preamble for section_mode_neff()/section_mode_fields(): the
        minimal self-contained meep Simulation for this section's
        cross-section, initialized and ready for get_eigenmode() calls, plus
        the eigenmode monitor Volume. Caller owns sim.reset_meep().

        Split out of section_mode_neff() verbatim (no behavior change) so
        section_mode_fields() can reuse the identical sim setup rather than
        duplicating it.

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

        si_med = self.material_stack.si_medium(wl_center)
        sio2_med = self.material_stack.sio2_medium(wl_center)

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

    def _iter_distinct_bands(self, sim, where, f, mea_c, mea_w, n_modes):
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
            em = sim.get_eigenmode(
                f, mp.X, where, band_num, mp.Vector3(f * _KPOINT_SEED_NEFF, 0, 0),
                parity=mp.NO_PARITY,
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
        # _iter_distinct_bands() guards against get_eigenmode() returning the
        # same eigenvector for two different band numbers (see its docstring
        # / validation_markdown/modal_multiplexer_validation.md gap 6).
        for band, em in enumerate(
            self._iter_distinct_bands(sim, where, f, mea_c, mea_w, n_modes), start=1
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
        **geometry_kwargs,
    ) -> np.ndarray:
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

        Returns:
            Real ndarray ``[len(wavelengths), n_modes]`` of β in rad/µm.

        Raises:
            ValueError: if the branch-selection margin is not comfortably below
                1 for some mode, i.e. the lengths are spaced too widely for
                MPB's seed to identify which wrap the FDTD phase is on. Fix by
                passing closer-spaced `lengths`.
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
        seed = self.section_mode_neff(
            geometry, wavelengths, mid - seed_span / 2.0, mid + seed_span / 2.0,
            n_modes,
            monitor_center=iface["monitor_center"],
            monitor_size=iface["monitor_size"],
            **geometry_kwargs,
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

        phases = np.zeros((len(lengths), len(wavelengths), n_modes))
        for j, L in enumerate(lengths):
            T = self.solve_section(
                geometry, wavelengths, mid - L / 2.0, mid + L / 2.0,
                n_modes_in=n_modes, n_modes_out=n_modes,
                save_dir=save_dir, cache_only=cache_only,
                **iface, **geometry_kwargs,
            )
            phases[j] = np.angle(np.diagonal(T, axis1=1, axis2=2))

        beta = np.zeros((len(wavelengths), n_modes))
        worst_margin = 0.0
        worst_resid = 0.0
        for i in range(len(wavelengths)):
            for m in range(n_modes):
                b, resid, margin = beta_from_phases(
                    lengths, phases[:, i, m], seed_beta[i, m])
                beta[i, m] = b
                worst_margin = max(worst_margin, margin)
                worst_resid = max(worst_resid, resid)

        # Guard 1: the phases must actually lie on a line. Fails if the seed
        # straddled a branch boundary (wrap counts inconsistent between
        # lengths) or if the section is not uniform after all.
        if worst_resid > _BETA_FIT_RESIDUAL_MAX:
            raise ValueError(
                f"section_mode_beta: phase-fit residual {worst_resid:.3f} rad "
                f"exceeds {_BETA_FIT_RESIDUAL_MAX} at lengths {lengths}. The "
                f"measured phase is not linear in length, so either the wrap "
                f"counts are inconsistent between lengths or [{x0}, {x1}] is "
                f"not a uniform cross-section."
            )
        logger.info(
            "section_mode_beta: %d lengths %s, worst phase-fit residual %.2e rad, "
            "worst branch margin %.3f%s",
            len(lengths), lengths, worst_resid, worst_margin,
            "" if worst_margin < _BETA_BRANCH_MARGIN_MAX
            else "  <-- moved far from MPB's seed; check this is real",
        )
        return beta

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
            Complex ndarray [len(wavelengths), n_modes_out, n_modes_in].
        """
        if (y0 is None) != (y1 is None):
            raise ValueError("solve_section: y0 and y1 must both be given or both omitted.")
        for _label, _port in (("excite_port", excite_port), ("measure_port", measure_port)):
            if _port is not None and _port.axis == "y" and y0 is None:
                raise ValueError(
                    f"solve_section: {_label}.axis='y' requires y0 and y1 to be given "
                    "(there is no y-extent to place a y-normal plane in otherwise)."
                )

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

        family_config = None
        if save_dir and monitor_size is not None:
            family_config = build_run_config(
                self, "solve_section",
                dict(x0=x0, x1=x1, wavelengths=wavelengths,
                     stub_left=stub_left, stub_right=stub_right,
                     monitor_center=monitor_center, monitor_size=monitor_size,
                     y0=y0, y1=y1, excite_port=_port_key(excite_port),
                     measure_port=_port_key(measure_port),
                     geometry_kwargs=geometry_kwargs, reverse=reverse,
                     **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor),
                     **resolved_material_config(self.material_stack, wavelengths)),
                geometry,
            )

        if save_dir:
            config = build_run_config(
                self, "solve_section",
                dict(x0=x0, x1=x1, wavelengths=wavelengths,
                     n_modes_in=n_modes_in, n_modes_out=n_modes_out,
                     stub_left=stub_left, stub_right=stub_right,
                     monitor_center=monitor_center, monitor_size=monitor_size,
                     y0=y0, y1=y1, excite_port=_port_key(excite_port),
                     measure_port=_port_key(measure_port),
                     geometry_kwargs=geometry_kwargs,
                     resolved_layout=(asdict(layout) if layout is not None else None),
                     **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor),
                     **resolved_material_config(self.material_stack, wavelengths)),
                geometry,
            )
            hit = lookup(save_dir, config)
            if hit is not None:
                logger.info("solve_section: cache hit at %s", hit)
                return load_result(hit)
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
        si_med = self.material_stack.si_medium(wl_center)
        sio2_med = self.material_stack.sio2_medium(wl_center)

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

        T = np.zeros((len(wavelengths), n_modes_out, n_modes_in), dtype=complex)

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
        else:
            ex_axis = "x"
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
            mea_along = (hi_ax if p.sign == 1 else lo_ax) - p.sign * self.monitor_margin - along_shift
            mea_transverse = p.center - transverse_shift
            mea_alpha_idx = p.alpha_idx
        else:
            mea_axis = "x"
            mea_width = layout.mea_w
            mea_transverse = layout.mea_c - y_center
            mea_along = layout.mea_x + lo - x_center
            mea_alpha_idx = layout.alpha_idx

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
            mea_flux = sim.add_mode_monitor(
                mon_freqs,
                mp.FluxRegion(center=mp.Vector3(*_port_point(mea_axis, mea_along, mea_transverse)),
                              size=mp.Vector3(*_port_size(mea_axis, mea_width, self._mon_sz))),
            )

            sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol))

            norm_coeff = sim.get_eigenmode_coefficients(norm_flux, [k], eig_parity=mp.NO_PARITY)
            if ex_axis == "y":
                # Flux-calibrated normalization: get_eigenmode_coefficients()
                # alone undercounts power for an off-axis, backward-launched
                # source — an empirically-found correction, first needed in
                # examples/ring_resonator/python/ring_resonator_hybrid.py's
                # get_norm_calibrated() (magnitude from mp.get_fluxes(),
                # phase from the coefficient). Applied for both signs
                # conservatively; only sign=-1 was ever tested there.
                flux = np.abs(np.asarray(mp.get_fluxes(norm_flux)))
                src_alpha = np.sqrt(flux) * np.exp(1j * np.angle(norm_coeff.alpha[0, :, ex_alpha_idx]))
            else:
                src_alpha = norm_coeff.alpha[0, :, ex_alpha_idx]

            for m in range(1, n_measure + 1):
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

        warn_on_transfer_duplicate(T, logger, x0, x1)

        if save_dir:
            folder = save_run(save_dir, config, T)
            logger.info("solve_section: cache saved at %s", folder)
            if family_config is not None:
                update_family_pointer(
                    save_dir, family_config, (n_modes_out, n_modes_in), os.path.basename(folder),
                )

        return T
