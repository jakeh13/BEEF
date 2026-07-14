"""FDTDSolver: meep-based FDTD solver for MMI and related photonic devices."""
from __future__ import annotations

from typing import List, Optional

import numpy as np

from ..geometry.base import DeviceGeometry
from ..materials.stack import MaterialStack, SOI_160
from ..results.smatrix import SMatrix


class FDTDSolver:
    """MEEP-based FDTD solver.

    Device geometry is converted to meep Prism objects via
    DeviceGeometry.to_meep_prisms() — no GDS files are written.

    Args:
        res: mesh resolution (pixels/µm).
        sim_width: simulation y-span (µm).
        sim_height: slab z-height reference for material lookup (µm).
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
            legacy (pre-FECES) FDTD/EME hybrid reference implementation
            (FECES_Legacy/masters/mmi_2x2_master.py's _port_simulation:
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
            needed/correct. Not used by solve() — that method's source is
            always the one true input port, with no neighboring port nearby
            to avoid, so it always spans the full sim_width regardless of
            this setting.
        material_stack: material index lookup for Si / SiO₂.
        decay_tol: DFT field decay tolerance for early-stop criterion.
            Lower = more accurate but slower; 1e-9 is faster, 1e-11 is tight.

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
        material_stack: MaterialStack = SOI_160,
        decay_tol: float = 1e-11,
    ):
        self.res = res
        self.sim_width = sim_width
        self.sim_height = sim_height
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
        # Derived: stub long enough that source + monitors sit in physical waveguide
        self.waveguide_length = max(source_margin, monitor_margin) + 0.25

    def _build_domain(self, geometry: DeviceGeometry, wavelengths: np.ndarray):
        """Return (Sx, Sy, wl_center, si_medium, sio2_medium).

        sim_width/sim_height are inner dimensions (excluding PML).
        Sx and Sy include 2*dpml on each axis.
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
        Sx, Sy, wl_center, si_med, sio2_med = self._build_domain(geometry, wavelengths)
        fcen = 1.0 / wl_center
        fwidth = (1.0 / float(wavelengths.min()) - 1.0 / float(wavelengths.max())) * 0.5
        nfreq = len(wavelengths)

        geometry_meep = geometry.to_meep_prisms(si_med, stub_length=self.waveguide_length + self.dpml)

        # Source position: centered on the actual input waveguide via
        # geometry.input_ys, not a hand-rolled +-waveguide_offset/2 offset —
        # that formula ignored the geometry's real port positions (e.g. it
        # placed the source off-axis for a single-input device, whose one
        # true input sits at input_ys[0]=0.0, not +-waveguide_offset/2).
        input_ys = getattr(geometry, "input_ys", [0.0])
        input_idx = len(input_ys) - 1 if port == "top" else 0
        src_yo = input_ys[input_idx]
        src_x = -Sx / 2 + self.dpml + self.source_margin
        # This is always the one true input port (no neighboring port
        # nearby that a half-domain aperture would need to avoid, unlike a
        # closely-spaced output branch), so it always spans the full
        # non-PML domain — PML edge to PML edge, exclusive of PML — rather
        # than half_domain_aperture (that setting is for narrow output ports).
        src_width = self.sim_width

        src = mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=fwidth),
            center=mp.Vector3(src_x, src_yo),
            size=mp.Vector3(0, src_width),   # span inner sim only, not PML
            eig_band=input_mode,
            eig_kpoint=mp.Vector3(1, 0, 0),
            eig_match_freq=True,
        )

        in_mon_x = -Sx / 2 + self.dpml + self.monitor_margin
        out_mon_x = Sx / 2 - self.dpml - self.monitor_margin

        sim = mp.Simulation(
            resolution=self.res,
            cell_size=mp.Vector3(Sx, Sy),
            boundary_layers=[mp.PML(self.dpml)],
            geometry=geometry_meep,
            sources=[src],
            default_material=sio2_med,
        )

        # Input normalization monitor — same aperture as the source
        in_flux = sim.add_mode_monitor(
            fcen, fwidth, nfreq,
            mp.FluxRegion(center=mp.Vector3(in_mon_x, src_yo), size=mp.Vector3(0, src_width)),
        )

        # One output monitor per output waveguide
        out_widths = (
            self.output_monitor_widths
            if self.output_monitor_widths is not None
            else [self.sim_width] * len(self.output_monitor_centers)
        )
        out_fluxes = [
            sim.add_mode_monitor(
                fcen, fwidth, nfreq,
                mp.FluxRegion(
                    center=mp.Vector3(out_mon_x, yc),
                    size=mp.Vector3(0, yw),
                ),
            )
            for yc, yw in zip(self.output_monitor_centers, out_widths)
        ]

        sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol))

        # Normalise by forward input-mode amplitude so |coeff|² = physical T
        in_coeff = sim.get_eigenmode_coefficients(
            in_flux, [input_mode], eig_parity=mp.NO_PARITY
        )
        src_alpha = in_coeff.alpha[0, :, 0]  # [nfreq] forward amplitude

        smat = SMatrix(wavelengths=wavelengths)
        for i, flux_mon in enumerate(out_fluxes):
            out_port = f"out_{i}"
            for mode_n in range(1, self.num_output_modes + 1):
                coeff = sim.get_eigenmode_coefficients(
                    flux_mon, [mode_n], eig_parity=mp.NO_PARITY
                )
                alpha = coeff.alpha[0, :, 0] / src_alpha
                smat.set(alpha, port=out_port, in_mode=input_mode, out_mode=mode_n)

        sim.reset_meep()
        return smat

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
            **geometry_kwargs: forwarded to geometry.section_prisms/width_at/
                center_at (e.g. port_idx for a branching output geometry).

        Returns:
            Complex ndarray [len(wavelengths), n_modes_out, n_modes_in].
        """
        import meep as mp

        wavelengths = np.asarray(wavelengths)
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
            **geometry_kwargs,
        )

        span = x1 - x0
        eps = min(1e-3, span * 0.01)
        # The launch aperture must span most of the non-PML domain — not just
        # the local waveguide core — so MPB's eigenmode solve for the source
        # includes the mode's evanescent tail; a source/monitor sized to the
        # bare core (e.g. a 0.35 µm single-mode waveguide in a 5 µm domain)
        # truncates that tail and corrupts the mode solve, which then
        # corrupts every normalization downstream. The far (measurement) side
        # is correctly sized to the local cross-section width ONLY when it's
        # a genuine single-mode port (well-isolated waveguide, small tail).
        # Whenever that side must resolve >1 mode (e.g. a directional
        # coupler's coupled-pair interface), width_at() there is just the
        # bare silicon span with zero cladding margin — nowhere near enough
        # room for MPB to resolve closely-spaced supermodes, whose whole
        # distinguishing feature is the evanescent tail shape. Widened to
        # self.sim_width below once n_measure is known.
        c_left = geometry.center_at(x0 + eps, **geometry_kwargs)
        c_right = geometry.center_at(x1 - eps, **geometry_kwargs)
        w_meas_left = geometry.width_at(x0 + eps, **geometry_kwargs)
        w_meas_right = geometry.width_at(x1 - eps, **geometry_kwargs)

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
        fwidth = (1.0 / float(wavelengths.min()) - 1.0 / float(wavelengths.max())) * 0.5
        nfreq = len(wavelengths)

        # shift prisms from device frame into the meep-centered frame
        shifted = [
            mp.Prism(
                vertices=[mp.Vector3(v.x - x_center, v.y, v.z) for v in p.vertices],
                height=mp.inf, material=p.material,
            )
            for p in prisms
        ]

        T = np.zeros((len(wavelengths), n_modes_out, n_modes_in), dtype=complex)

        # Excite from whichever boundary has fewer modes — min(n_in, n_out)
        # simulations instead of always n_modes_in. Reciprocity guarantees
        # the same coefficient matrix either way (see docstring).
        #
        # source_margin/monitor_margin are documented (see class docstring) as
        # offsets from the PML *inner* edge — i.e. from lo/hi (device frame),
        # not from x0/x1 directly; those only coincide when stub_left/right is
        # 0 (a pure interior cut). Shifted into the meep frame: lo - x_center
        # and hi - x_center. monitor_margin only offsets the source-side
        # normalization monitor (norm_x, collocated with the source, kept
        # clear of its near-field launch region). The measured/interface side
        # (mea_x) sits exactly at the PML interface (lo_m/hi_m, with no
        # margin) so that this section's measured plane exactly matches an
        # adjacent section's own modeled boundary at that same x — e.g.
        # EMESolver.solve_section()'s `length = (x1-x0)+stub_left+stub_right`
        # starts/ends exactly at x0/x1, not at x0+margin/x1-margin.
        lo_m, hi_m = lo - x_center, hi - x_center
        if n_modes_out <= n_modes_in:
            n_excite, n_measure = n_modes_out, n_modes_in
            src_x, norm_x, mea_x = hi_m - self.source_margin, hi_m - self.monitor_margin, lo_m
            src_c, mea_c = c_right, c_left
            w_meas = self.sim_width if n_measure > 1 else w_meas_left
            kpoint = mp.Vector3(-1, 0, 0)
            alpha_idx = 1  # backward amplitude
        else:
            n_excite, n_measure = n_modes_in, n_modes_out
            src_x, norm_x, mea_x = lo_m + self.source_margin, lo_m + self.monitor_margin, hi_m
            src_c, mea_c = c_left, c_right
            w_meas = self.sim_width if n_measure > 1 else w_meas_right
            kpoint = mp.Vector3(1, 0, 0)
            alpha_idx = 0  # forward amplitude

        # n_excite == 1 means this boundary is a genuinely single-mode port
        # (e.g. one arm of a directional coupler, well after its S-bend). A
        # full sim_width aperture centered on the arm would also reach a
        # second, physically-present-but-unexcited waveguide sitting nearby;
        # half_domain_aperture instead spans only sim_width/2 — from the
        # shared symmetry axis (y=0) out to the near PML edge on whichever
        # side the arm sits — centered at +-sim_width/4, matching the legacy
        # reference implementation's validated convention (see
        # half_domain_aperture's docstring). A multi-mode excite side (e.g.
        # the coupling-region interface, where the aperture must see every
        # arm) always keeps the full sim_width, centered on its own src_c.
        if n_excite == 1 and self.half_domain_aperture:
            src_width = self.sim_width / 2.0
            src_c_ap = self.sim_width / 4.0 if src_c > 0 else (-self.sim_width / 4.0 if src_c < 0 else 0.0)
        else:
            src_width = self.sim_width
            src_c_ap = src_c

        for k in range(1, n_excite + 1):
            src = mp.EigenModeSource(
                mp.GaussianSource(fcen, fwidth=fwidth),
                center=mp.Vector3(src_x, src_c_ap),
                size=mp.Vector3(0, src_width),
                eig_band=k,
                eig_kpoint=kpoint,
                eig_match_freq=True,
            )
            sim = mp.Simulation(
                resolution=self.res,
                cell_size=mp.Vector3(Sx, Sy),
                boundary_layers=[mp.PML(self.dpml)],
                geometry=shifted,
                sources=[src],
                default_material=sio2_med,
            )
            norm_flux = sim.add_mode_monitor(
                fcen, fwidth, nfreq,
                mp.FluxRegion(center=mp.Vector3(norm_x, src_c_ap), size=mp.Vector3(0, src_width)),
            )
            mea_flux = sim.add_mode_monitor(
                fcen, fwidth, nfreq,
                mp.FluxRegion(center=mp.Vector3(mea_x, mea_c), size=mp.Vector3(0, w_meas)),
            )

            sim.run(until_after_sources=mp.stop_when_dft_decayed(tol=self.decay_tol))

            norm_coeff = sim.get_eigenmode_coefficients(norm_flux, [k], eig_parity=mp.NO_PARITY)
            src_alpha = norm_coeff.alpha[0, :, alpha_idx]

            for m in range(1, n_measure + 1):
                c = sim.get_eigenmode_coefficients(mea_flux, [m], eig_parity=mp.NO_PARITY)
                coeff = c.alpha[0, :, alpha_idx] / src_alpha
                if n_modes_out <= n_modes_in:
                    T[:, k - 1, m - 1] = coeff  # excited=out(k), measured=in(m)
                else:
                    T[:, m - 1, k - 1] = coeff  # excited=in(k), measured=out(m)

            sim.reset_meep()

        return T
