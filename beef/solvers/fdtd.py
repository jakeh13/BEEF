"""FDTDSolver: meep-based FDTD solver for MMI and related photonic devices."""
from __future__ import annotations

import logging
import os
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..cache import (
    CacheMissError, build_run_config, build_whole_device_config, load_result, lookup,
    lookup_superset, save_run, update_family_pointer,
)
from ..geometry.base import DeviceGeometry
from ..materials.stack import MaterialStack, SOI_160
from ..results.smatrix import SMatrix
from .._log import get_logger

logger = get_logger(__name__)

# Starting guess for section_mode_neff()'s eigensolves, as an effective index.
# Only a seed — meep's match_frequency iterates from here to the true k — so it
# just has to sit in the guided range for a Si-core/oxide-clad waveguide
# (1.44 < n_eff < 3.48). Roughly the middle of that window.
_KPOINT_SEED_NEFF = 2.5


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
                     **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor)),
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

        Returns:
            Real ndarray ``[len(wavelengths), n_modes]`` of effective indices,
            in the caller's wavelength order.
        """
        import meep as mp

        wavelengths = np.asarray(wavelengths, dtype=float)
        wl_center = float(np.mean(wavelengths))
        si_med = self.material_stack.si_medium(wl_center)
        sio2_med = self.material_stack.sio2_medium(wl_center)

        prisms = geometry.section_prisms(
            x0, x1, si_med,
            stub_left=self.dpml, stub_right=self.dpml,
            wg_height=self._wg_height,
            **geometry_kwargs,
        )
        span = x1 - x0
        Sx = span + 2 * self.dpml
        Sy = self.sim_width + 2 * self.dpml
        lo, hi = x0 - self.dpml, x1 + self.dpml
        x_center = (lo + hi) / 2.0
        # Same vertex-wise shift into meep's origin-centred cell that
        # solve_section() does — prisms carry vertices, not a centre.
        shifted = [
            mp.Prism(
                vertices=[mp.Vector3(v.x - x_center, v.y, v.z) for v in p.vertices],
                height=p.height, material=p.material,
            )
            for p in prisms
        ]

        mea_c = 0.0 if monitor_center is None else float(monitor_center)
        mea_w = self.sim_width if monitor_size is None else float(monitor_size)

        sim = mp.Simulation(
            resolution=self.res,
            cell_size=mp.Vector3(Sx, Sy, self._cell_sz),
            boundary_layers=[mp.PML(self.dpml)],
            geometry=shifted,
            default_material=sio2_med,
            eps_averaging=self.eps_averaging,
            subpixel_tol=self.subpixel_tol,
            subpixel_maxeval=self.subpixel_maxeval,
        )
        sim.init_sim()

        where = mp.Volume(center=mp.Vector3(0, mea_c),
                          size=mp.Vector3(0, mea_w, self._mon_sz))
        neff = np.zeros((len(wavelengths), n_modes))
        for i, wl in enumerate(wavelengths):
            f = 1.0 / float(wl)
            for band in range(1, n_modes + 1):
                # kpoint only seeds the solve; match_frequency then iterates to
                # the true k. Any guided mode sits between the cladding and
                # core indices, so the midpoint of that range is a seed that
                # cannot be badly wrong for any waveguide this solver handles.
                em = sim.get_eigenmode(
                    f, mp.X, where, band, mp.Vector3(f * _KPOINT_SEED_NEFF, 0, 0),
                    parity=mp.NO_PARITY,
                )
                neff[i, band - 1] = float(np.real(em.k.x)) / f

        sim.reset_meep()
        return neff

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
        layout = section_fdtd_layout(
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
        family_config = None
        if save_dir and monitor_size is not None:
            family_config = build_run_config(
                self, "solve_section",
                dict(x0=x0, x1=x1, wavelengths=wavelengths,
                     stub_left=stub_left, stub_right=stub_right,
                     monitor_center=monitor_center, monitor_size=monitor_size,
                     geometry_kwargs=geometry_kwargs, reverse=layout.reverse,
                     **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor)),
                geometry,
            )

        if save_dir:
            config = build_run_config(
                self, "solve_section",
                dict(x0=x0, x1=x1, wavelengths=wavelengths,
                     n_modes_in=n_modes_in, n_modes_out=n_modes_out,
                     stub_left=stub_left, stub_right=stub_right,
                     monitor_center=monitor_center, monitor_size=monitor_size,
                     geometry_kwargs=geometry_kwargs,
                     resolved_layout=asdict(layout),
                     **_resolved_spectral_config(wavelengths, self.source_bandwidth_factor)),
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
        Sy = self.sim_width + 2 * self.dpml

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
        # polygons_to_prisms in beef/geometry/base.py)
        shifted = [
            mp.Prism(
                vertices=[mp.Vector3(v.x - x_center, v.y, v.z) for v in p.vertices],
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
        n_excite, n_measure = layout.n_excite, layout.n_measure
        src_x = layout.src_x + lo - x_center
        norm_x = layout.norm_x + lo - x_center
        mea_x = layout.mea_x + lo - x_center
        src_c_ap, src_width = layout.src_c, layout.src_width
        mea_c, w_meas = layout.mea_c, layout.mea_w
        kpoint = mp.Vector3(*layout.kpoint)
        alpha_idx = layout.alpha_idx

        for k in range(1, n_excite + 1):
            src = mp.EigenModeSource(
                mp.GaussianSource(fcen, fwidth=fwidth),
                center=mp.Vector3(src_x, src_c_ap),
                size=mp.Vector3(0, src_width, self._mon_sz),
                eig_band=k,
                eig_kpoint=kpoint,
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
                mp.FluxRegion(center=mp.Vector3(norm_x, src_c_ap),
                              size=mp.Vector3(0, src_width, self._mon_sz)),
            )
            mea_flux = sim.add_mode_monitor(
                mon_freqs,
                mp.FluxRegion(center=mp.Vector3(mea_x, mea_c),
                              size=mp.Vector3(0, w_meas, self._mon_sz)),
            )

            sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol))

            norm_coeff = sim.get_eigenmode_coefficients(norm_flux, [k], eig_parity=mp.NO_PARITY)
            src_alpha = norm_coeff.alpha[0, :, alpha_idx]

            for m in range(1, n_measure + 1):
                c = sim.get_eigenmode_coefficients(mea_flux, [m], eig_parity=mp.NO_PARITY)
                # [wl_order] lifts the ratio out of meep's ascending-frequency
                # monitor order back onto `wavelengths` as the caller gave it.
                coeff = (c.alpha[0, :, alpha_idx] / src_alpha)[wl_order]
                if layout.reverse:
                    T[:, k - 1, m - 1] = coeff  # excited=out(k), measured=in(m)
                else:
                    T[:, m - 1, k - 1] = coeff  # excited=in(k), measured=out(m)

            sim.reset_meep()

        if save_dir:
            folder = save_run(save_dir, config, T)
            logger.info("solve_section: cache saved at %s", folder)
            if family_config is not None:
                update_family_pointer(
                    save_dir, family_config, (n_modes_out, n_modes_in), os.path.basename(folder),
                )

        return T
