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
            port: input port label — "top" places source at +wg_offset/2,
                  "bot" at −wg_offset/2.  For single-input devices (1×2 MMI)
                  both land at y=0 effectively (full-width monitor).
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

        # Source position (y derived from device port geometry)
        wg_offset = getattr(geometry, "waveguide_offset", 0.0) / 2
        src_yo = wg_offset if port == "top" else -wg_offset
        src_x = -Sx / 2 + self.dpml + self.source_margin

        src = mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=fwidth),
            center=mp.Vector3(src_x, src_yo),
            size=mp.Vector3(0, self.sim_width),   # span inner sim only, not PML
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

        # Input normalization monitor — span inner sim width, not PML
        in_flux = sim.add_mode_monitor(
            fcen, fwidth, nfreq,
            mp.FluxRegion(center=mp.Vector3(in_mon_x, src_yo), size=mp.Vector3(0, self.sim_width)),
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

        prisms = geometry.section_prisms(
            x0, x1, si_med, stub_left=stub_left, stub_right=stub_right, **geometry_kwargs
        )

        span = x1 - x0
        eps = min(1e-3, span * 0.01)
        # The launch aperture must span most of the non-PML domain — not just
        # the local waveguide core — so MPB's eigenmode solve for the source
        # includes the mode's evanescent tail; a source/monitor sized to the
        # bare core (e.g. a 0.35 µm single-mode waveguide in a 5 µm domain)
        # truncates that tail and corrupts the mode solve, which then
        # corrupts every normalization downstream. The far (measurement) side
        # is a genuine per-interface decomposition and is correctly sized to
        # the local cross-section width.
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
        # and hi - x_center.
        lo_m, hi_m = lo - x_center, hi - x_center
        if n_modes_out <= n_modes_in:
            n_excite, n_measure = n_modes_out, n_modes_in
            src_x, norm_x, mea_x = hi_m - self.source_margin, hi_m - self.monitor_margin, lo_m + self.monitor_margin
            src_c, mea_c = c_right, c_left
            w_meas = w_meas_left
            kpoint = mp.Vector3(-1, 0, 0)
            alpha_idx = 1  # backward amplitude
        else:
            n_excite, n_measure = n_modes_in, n_modes_out
            src_x, norm_x, mea_x = lo_m + self.source_margin, lo_m + self.monitor_margin, hi_m - self.monitor_margin
            src_c, mea_c = c_left, c_right
            w_meas = w_meas_right
            kpoint = mp.Vector3(1, 0, 0)
            alpha_idx = 0  # forward amplitude

        for k in range(1, n_excite + 1):
            src = mp.EigenModeSource(
                mp.GaussianSource(fcen, fwidth=fwidth),
                center=mp.Vector3(src_x, src_c),
                size=mp.Vector3(0, self.sim_width),
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
                mp.FluxRegion(center=mp.Vector3(norm_x, src_c), size=mp.Vector3(0, self.sim_width)),
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
