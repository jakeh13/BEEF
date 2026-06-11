import argparse
import dataconf
import numpy as np
import logging
from matplotlib import pyplot as plt
import pickle
import os
import json
from dataclasses import asdict

import meow as mw
import meep as mp

from meow import compute_modes_tidy3d as compute_modes

from engine.hyperparams.MMI_MUX import HyperParamsMMIMUX as HyperParams

from engine.helpers.general import make_output_folders, resolve_output_folder
from engine.helpers.meep import (
    meep_material_from_meow,
    prepare_gds_for_meep,
    create_temp_gds_for_full_mmi_mux,
    create_temp_gds_for_mmi_mux_input,
    create_temp_gds_for_mmi_mux_output,
)
from engine.helpers.meow import (
    load_MMI_MUX_full,
    load_MMI_MUX_cascade,
    load_MMI_MUX_port,
    stack_extrusions_GF,
)
from engine.plotting.meow import plot_eme_fields


# Port definitions: (y_offset, wg_width, label)
PORT_DEFS = [
    ("top", +1,  "top_waveguide_width"),   # port_idx 0
    ("mid",  0,  "waveguide_width"),        # port_idx 1
    ("bot", -1,  "waveguide_width"),        # port_idx 2
]


class MMI_MUX:
    def __init__(self, args):
        self.config = args.config
        self.hp: HyperParams = dataconf.load(self.config, HyperParams)
        self.hp.outputfolder = resolve_output_folder(self.hp.outputfolder)
        make_output_folders(self.hp.outputfolder)

        s = json.dumps(asdict(self.hp), indent=4, default=lambda x: x.name)
        with open(f"{self.hp.outputfolder}/config.json", "w") as f:
            print(s, file=f)

        handler = logging.FileHandler(f"{self.hp.outputfolder}/logs/master.log")
        handler.setLevel("INFO")
        handler.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)-4s | %(filename)s | %(funcName)s:%(lineno)04d | %(message)s"
        ))
        self.master_logger = logging.getLogger(__name__.replace("_", ""))
        self.master_logger.addHandler(handler)
        self.master_logger.setLevel("INFO")

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def _wavelengths(self):
        return np.linspace(self.hp.sim.wavelength_min, self.hp.sim.wavelength_max, self.hp.sim.wavelength_pts)

    def _wl_center(self):
        return (self.hp.sim.wavelength_max + self.hp.sim.wavelength_min) / 2

    def _materials(self, wl_center):
        Si   = meep_material_from_meow(mw.silicon.n,       mw.silicon.params['wl'],       wl_center)
        SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'], wl_center)
        return Si, SiO2

    def _port_y(self, port_idx):
        sign = PORT_DEFS[port_idx][1]
        return sign * self.hp.geo.waveguide_offset

    def _port_width(self, port_idx):
        attr = PORT_DEFS[port_idx][2]
        return getattr(self.hp.geo, attr)

    # ------------------------------------------------------------------
    # FDTD: input coupling simulation
    # ------------------------------------------------------------------

    def input_simulation(self, wavelengths, mode: int = 1):
        """FDTD of [input MM WG → input taper → partial MMI 1 body].

        Returns array [n_wl, num_modes_MMI]: coupling amplitude from WG mode
        to each MMI 1 body mode at the monitor plane inside MMI 1.
        """
        if mp.am_master(): self.master_logger.info(f"input_simulation mode={mode}")
        geo, sim = self.hp.geo, self.hp.sim
        Sx = sim.waveguide_length + geo.taper_length + sim.mm_sim_length + 2*sim.dpml
        Sy = sim.oxide_width  + 2*sim.dpml
        Sz = sim.oxide_height + 2*sim.dpml
        wl_c = self._wl_center()
        freqs = 1 / wavelengths.astype(list)
        Si, SiO2 = self._materials(wl_c)

        # input MM WG stub
        geometry = [mp.Block(
            center=mp.Vector3(-Sx/2 + sim.waveguide_length/2 + sim.dpml/2, 0, 0),
            size=mp.Vector3(sim.waveguide_length + sim.dpml, geo.mm_waveguide_width, geo.waveguide_height),
            material=Si
        )]

        if mp.am_master():
            create_temp_gds_for_mmi_mux_input(geo, sim, Sx, name=f"temp_mux_in_{mode}")
        mp.all_wait()
        mmi = mp.get_GDSII_prisms(Si, f"gds/temp_mux_in_{mode}.gds", 2,
                                   -geo.waveguide_height/2, geo.waveguide_height/2)
        geometry += mmi

        sources = [mp.EigenModeSource(
            mp.GaussianSource(frequency=1/wl_c, fwidth=0.4/wl_c),
            eig_band=mode,
            direction=mp.NO_DIRECTION,
            eig_kpoint=mp.Vector3(x=1),
            size=mp.Vector3(0, sim.oxide_width/2, sim.oxide_height),
            center=mp.Vector3(-Sx/2 + sim.dpml + sim.source_position, 0, 0),
        )]

        sim_obj = mp.Simulation(
            cell_size=mp.Vector3(Sx, Sy, Sz),
            boundary_layers=[mp.PML(sim.dpml)],
            geometry=geometry,
            sources=sources,
            eps_averaging=True,
            resolution=sim.res,
            default_material=SiO2,
        )

        in_flux = sim_obj.add_mode_monitor(freqs, mp.FluxRegion(
            center=mp.Vector3(-Sx/2 + sim.dpml + sim.sm_monitor_position, 0, 0),
            size=mp.Vector3(0, sim.oxide_width/2, Sz - 2*sim.dpml),
        ))
        out_flux = sim_obj.add_mode_monitor(freqs, mp.FluxRegion(
            center=mp.Vector3(-Sx/2 + sim.dpml + sim.waveguide_length + geo.taper_length + sim.mm_monitor_position, 0, 0),
            size=mp.Vector3(0, Sy - 2*sim.dpml, Sz - 2*sim.dpml),
        ))

        if mode == 1:
            if mp.am_master(): self.master_logger.info("Plotting input_simulation geometry")
            plt.figure(); sim_obj.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, 0)))
            plt.savefig(f'{self.hp.outputfolder}/layout/mux_in_xy.png', dpi=300); plt.close()

        if mp.am_master(): self.master_logger.info(f"Running input_simulation mode={mode}")
        sim_obj.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol=1e-11))

        src_coeff = sim_obj.get_eigenmode_coefficients(
            in_flux, [mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
        coeffs = []
        for n in range(1, sim.num_modes_MMI + 1):
            c = sim_obj.get_eigenmode_coefficients(
                out_flux, [n], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
            s = c.alpha[:, :, 0] / src_coeff.alpha[:, :, 0]
            coeffs.append(s[0, :])
            if mp.am_master(): self.master_logger.info(f"  MMI mode {n}: power={np.abs(s[0,:])**2}")
        return np.asarray(coeffs).T.astype(list)  # [n_wl, num_modes_MMI]

    # ------------------------------------------------------------------
    # FDTD: output coupling simulation (by reciprocity)
    # ------------------------------------------------------------------

    def output_simulation(self, wavelengths, port_idx: int, mode: int = 1):
        """FDTD of [partial MMI 2 body → output taper → output WG] with source in WG.

        Returns array [n_wl, num_modes_MMI]: coupling amplitude from MMI 2 mode
        to WG mode at the monitor plane inside MMI 2.
        """
        if mp.am_master(): self.master_logger.info(f"output_simulation port={port_idx} mode={mode}")
        geo, sim = self.hp.geo, self.hp.sim
        Sx = sim.waveguide_length + geo.taper_length + sim.mm_sim_length + 2*sim.dpml
        Sy = sim.oxide_width  + 2*sim.dpml
        Sz = sim.oxide_height + 2*sim.dpml
        wl_c = self._wl_center()
        freqs = 1 / wavelengths.astype(list)
        Si, SiO2 = self._materials(wl_c)

        py = self._port_y(port_idx)
        pw = self._port_width(port_idx)

        # output WG stub
        geometry = [mp.Block(
            center=mp.Vector3(Sx/2 - sim.waveguide_length/2 - sim.dpml/2, py, 0),
            size=mp.Vector3(sim.waveguide_length + sim.dpml, pw, geo.waveguide_height),
            material=Si,
        )]

        name = f"temp_mux_out_{port_idx}_{mode}"
        if mp.am_master():
            create_temp_gds_for_mmi_mux_output(geo, sim, Sx, port_idx, name=name)
        mp.all_wait()
        mmi = mp.get_GDSII_prisms(Si, f"gds/{name}.gds", 2,
                                   -geo.waveguide_height/2, geo.waveguide_height/2)
        geometry += mmi

        sources = [mp.EigenModeSource(
            mp.GaussianSource(frequency=1/wl_c, fwidth=0.4/wl_c),
            eig_band=mode,
            direction=mp.NO_DIRECTION,
            eig_kpoint=mp.Vector3(x=-1),  # propagate left (into MMI 2)
            size=mp.Vector3(0, sim.oxide_width/2, sim.oxide_height),
            center=mp.Vector3(Sx/2 - sim.dpml - sim.source_position, py, 0),
        )]

        sim_obj = mp.Simulation(
            cell_size=mp.Vector3(Sx, Sy, Sz),
            boundary_layers=[mp.PML(sim.dpml)],
            geometry=geometry,
            sources=sources,
            eps_averaging=True,
            resolution=sim.res,
            default_material=SiO2,
        )

        in_flux = sim_obj.add_mode_monitor(freqs, mp.FluxRegion(
            center=mp.Vector3(Sx/2 - sim.dpml - sim.sm_monitor_position, py, 0),
            size=mp.Vector3(0, sim.oxide_width/2, Sz - 2*sim.dpml),
        ))
        out_flux = sim_obj.add_mode_monitor(freqs, mp.FluxRegion(
            center=mp.Vector3(Sx/2 - sim.dpml - sim.waveguide_length - geo.taper_length - sim.mm_monitor_position, 0, 0),
            size=mp.Vector3(0, Sy - 2*sim.dpml, Sz - 2*sim.dpml),
        ))

        if mp.am_master(): self.master_logger.info(f"Running output_simulation port={port_idx} mode={mode}")
        sim_obj.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol=1e-11))

        src_coeff = sim_obj.get_eigenmode_coefficients(
            in_flux, [mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=-1))
        coeffs = []
        for n in range(1, sim.num_modes_MMI + 1):
            c = sim_obj.get_eigenmode_coefficients(
                out_flux, [n], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=-1))
            s = c.alpha[:, :, 0] / src_coeff.alpha[:, :, 0]
            coeffs.append(s[0, :])
            if mp.am_master(): self.master_logger.info(f"  MMI mode {n}: power={np.abs(s[0,:])**2}")
        return np.asarray(coeffs).T.astype(list)  # [n_wl, num_modes_MMI]

    # ------------------------------------------------------------------
    # EME: cascade mm_simulation (for HYBRID)
    # ------------------------------------------------------------------

    def mm_simulation(self, wavelengths):
        """EME of full MUX cascade from mm_monitor_position inside MMI 1
        to mm_monitor_position inside MMI 2 (from output face).

        Returns list of S-matrices [n_wl], each shape [2*num_modes_MMI, 2*num_modes_MMI].
        S[i, num_modes_MMI+j] = coupling from input mode j to output mode i.
        """
        geo, sim = self.hp.geo, self.hp.sim
        N = sim.num_modes_MMI
        Nt = sim.num_cells_taper
        Mp = sim.mm_monitor_position

        cascade_gds = load_MMI_MUX_cascade(self.hp)
        extrusion_rules = stack_extrusions_GF(geo.waveguide_height, sim.oxide_height)
        structs = mw.extrude_gds(cascade_gds, extrusion_rules)

        taper_L = geo.taper_length / Nt
        Ls = (
            [geo.mm_1_length - Mp]
            + [taper_L] * Nt
            + [geo.between_waveguide_length]
            + [taper_L] * Nt
            + [geo.mm_2_length - Mp]
        )
        num_modes = (
            [N]
            + [sim.num_modes_taper] * Nt
            + [sim.num_modes_taper]
            + [sim.num_modes_taper] * Nt
            + [N]
        )

        cells = mw.create_cells(
            structures=structs,
            mesh=mw.Mesh2D(
                x=np.linspace(-sim.oxide_width/2,  sim.oxide_width/2,  int(sim.res * sim.oxide_width)  + 1),
                y=np.linspace(-sim.oxide_height/2, sim.oxide_height/2, int(sim.res * sim.oxide_height) + 1),
            ),
            Ls=Ls,
        )

        s_matrices = []
        for wl in wavelengths.astype(list):
            if mp.am_master(): self.master_logger.info(f"mm_simulation EME wl={wl:.4f}")
            env = mw.Environment(wl=wl, T=25.0)
            css = [mw.CrossSection.from_cell(cell=cell, env=env) for cell in cells]
            modes = [compute_modes(cs, num_modes=nm) for cs, nm in zip(css, num_modes)]
            S, _ = mw.compute_s_matrix(modes, cells=cells)
            s_matrices.append(S)
        return s_matrices

    # ------------------------------------------------------------------
    # HYBRID simulation
    # ------------------------------------------------------------------

    def hybrid_simulation(self):
        wavelengths = self._wavelengths()
        input_modes  = self.hp.sim.input_modes
        output_modes = self.hp.sim.output_modes
        N = self.hp.sim.num_modes_MMI
        n_wl = len(wavelengths)

        # input FDTD: C_in[m] shape [n_wl, N]
        C_in = {}
        for m in set(input_modes):
            C_in[m] = np.array(self.input_simulation(wavelengths, m))

        # cascade EME: list of [2N,2N] S-matrices
        S_cascade = self.mm_simulation(wavelengths)

        # output FDTD by reciprocity: C_out[(port_idx,p)] shape [n_wl, N]
        output_ports = list(range(3))
        C_out = {}
        for port_idx in output_ports:
            for p in set(output_modes):
                C_out[(port_idx, p)] = np.array(self.output_simulation(wavelengths, port_idx, p))

        powers_by_port = {}
        angles_by_port = {}
        for in_m in input_modes:
            for port_idx in output_ports:
                for out_m in output_modes:
                    key = (in_m, port_idx, out_m)
                    pwr_list, ang_list = [], []
                    for pt in range(n_wl):
                        c_in = np.array(C_in[in_m][pt])          # [N]
                        S    = np.array(S_cascade[pt])            # [2N, 2N]
                        # S[i, N+j] = backward coupling (output j → input i)
                        # = forward coupling (input i → output j) by reciprocity
                        # field_out[j] = sum_i S[i, N+j] * c_in[i] = (S[:N,N:].T @ c_in)[j]
                        field_out = S[:N, N:].T @ c_in            # [N]
                        c_out = np.array(C_out[(port_idx, out_m)][pt])  # [N]
                        amplitude = np.dot(c_out, field_out)
                        pwr_list.append(np.abs(amplitude)**2)
                        ang_list.append(np.angle(amplitude, deg=True))
                        if mp.am_master():
                            self.master_logger.info(
                                f"in={in_m} port={PORT_DEFS[port_idx][0]} out={out_m} "
                                f"power={pwr_list[-1]:.4f} angle={ang_list[-1]:.2f}")
                    powers_by_port[key] = pwr_list
                    angles_by_port[key] = ang_list

        data = {
            'C_in': C_in, 'S_cascade': S_cascade, 'C_out': C_out,
            'powers_by_port': powers_by_port, 'angles_by_port': angles_by_port,
            'wavelengths': wavelengths,
        }
        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data, f)
        return powers_by_port, wavelengths

    # ------------------------------------------------------------------
    # Full FDTD simulation
    # ------------------------------------------------------------------

    def fdtd_simulation(self):
        if mp.am_master(): self.master_logger.info("Creating full FDTD simulation")
        geo, sim = self.hp.geo, self.hp.sim
        device_len = 4*geo.taper_length + geo.mm_1_length + geo.between_waveguide_length + geo.mm_2_length
        Sx = 2*sim.dpml + 2*sim.waveguide_length + device_len
        Sy = sim.oxide_width  + 2*sim.dpml
        Sz = sim.oxide_height + 2*sim.dpml
        wavelengths = self._wavelengths()
        wl_c = self._wl_center()
        freqs = 1 / wavelengths.astype(list)
        Si, SiO2 = self._materials(wl_c)

        # Build geometry once — shared across all input mode simulations
        geometry_base = [mp.Block(
            center=mp.Vector3(-Sx/2 + sim.waveguide_length/2 + sim.dpml/2, 0, 0),
            size=mp.Vector3(sim.waveguide_length + sim.dpml, geo.mm_waveguide_width, geo.waveguide_height),
            material=Si,
        )]
        for port_idx in range(3):
            py = self._port_y(port_idx)
            pw = self._port_width(port_idx)
            geometry_base.append(mp.Block(
                center=mp.Vector3(Sx/2 - sim.waveguide_length/2 - sim.dpml/2, py, 0),
                size=mp.Vector3(sim.waveguide_length + sim.dpml, pw, geo.waveguide_height),
                material=Si,
            ))

        if mp.am_master():
            create_temp_gds_for_full_mmi_mux(geo)
        mp.all_wait()
        mmi = mp.get_GDSII_prisms(Si, "gds/temp_mux_full.gds", 2,
                                   -geo.waveguide_height/2, geo.waveguide_height/2)
        geometry_base += mmi

        powers_by_port = {}
        angles_by_port = {}

        for in_m_idx, fdtd_in_mode in enumerate(sim.input_modes):
            if mp.am_master(): self.master_logger.info(f"Full FDTD input mode {fdtd_in_mode}")
            sources = [mp.EigenModeSource(
                mp.GaussianSource(frequency=1/wl_c, fwidth=0.4/wl_c),
                eig_band=fdtd_in_mode,
                direction=mp.NO_DIRECTION,
                eig_kpoint=mp.Vector3(x=1),
                size=mp.Vector3(0, sim.oxide_width/2, sim.oxide_height),
                center=mp.Vector3(-Sx/2 + sim.dpml + sim.source_position, 0, 0),
            )]

            sim_obj = mp.Simulation(
                cell_size=mp.Vector3(Sx, Sy, Sz),
                boundary_layers=[mp.PML(sim.dpml)],
                geometry=geometry_base,
                sources=sources,
                eps_averaging=True,
                resolution=sim.res,
                default_material=SiO2,
            )

            in_flux = sim_obj.add_mode_monitor(freqs, mp.FluxRegion(
                center=mp.Vector3(-Sx/2 + sim.dpml + sim.sm_monitor_position, 0, 0),
                size=mp.Vector3(0, sim.oxide_width/2, Sz - 2*sim.dpml),
            ))
            out_fluxes = []
            for port_idx in range(3):
                py = self._port_y(port_idx)
                out_fluxes.append(sim_obj.add_mode_monitor(freqs, mp.FluxRegion(
                    center=mp.Vector3(Sx/2 - sim.dpml - sim.sm_monitor_position, py, 0),
                    size=mp.Vector3(0, sim.oxide_width/2, Sz - 2*sim.dpml),
                )))

            if in_m_idx == 0:
                if mp.am_master(): self.master_logger.info("Plotting full FDTD geometry")
                plt.figure(); sim_obj.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, 0)))
                plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_xy.png', dpi=300); plt.close()
                plt.figure(); sim_obj.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, 0, Sz)))
                plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_xz.png', dpi=300); plt.close()
                plt.figure(); sim_obj.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(0, Sy, Sz)))
                plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_yz.png', dpi=300); plt.close()

            if mp.am_master(): self.master_logger.info(f"Running full FDTD (mode {fdtd_in_mode})")
            sim_obj.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol=1e-11))

            src_coeff = sim_obj.get_eigenmode_coefficients(
                in_flux, [fdtd_in_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))

            for port_idx, flux in enumerate(out_fluxes):
                for out_m in sim.output_modes:
                    coeff = sim_obj.get_eigenmode_coefficients(
                        flux, [out_m], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
                    s = coeff.alpha[:, :, 0] / src_coeff.alpha[:, :, 0]
                    s = s[0, :]
                    key = (fdtd_in_mode, port_idx, out_m)
                    pwr, ang = [], []
                    for pt in range(len(wavelengths)):
                        pwr.append(np.abs(s[pt])**2)
                        ang.append(np.angle(s[pt], deg=True))
                        if mp.am_master():
                            self.master_logger.info(
                                f"in={fdtd_in_mode} port={PORT_DEFS[port_idx][0]} out={out_m} "
                                f"power={pwr[-1]:.4f} angle={ang[-1]:.2f}")
                    powers_by_port[key] = pwr
                    angles_by_port[key] = ang

        data = {
            'powers_by_port': powers_by_port, 'angles_by_port': angles_by_port,
            'wavelengths': wavelengths,
        }
        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data, f)
        return powers_by_port, wavelengths

    # ------------------------------------------------------------------
    # FDTD from GDS
    # ------------------------------------------------------------------

    def fdtd_simulation_from_gds(self):
        if mp.am_master(): self.master_logger.info(f"Creating FDTD_GDS simulation from {self.hp.geo.gds_file}")
        geo, sim = self.hp.geo, self.hp.sim
        wl_c = self._wl_center()
        wavelengths = self._wavelengths()
        freqs = 1 / wavelengths.astype(list)
        Si, SiO2 = self._materials(wl_c)

        device_Sx = prepare_gds_for_meep(geo.gds_file, "gds/temp_gds_fdtd.gds")
        Sx = device_Sx + 2*sim.waveguide_length + 2*sim.dpml
        Sy = sim.oxide_width  + 2*sim.dpml
        Sz = sim.oxide_height + 2*sim.dpml

        # Build geometry once — shared across all input mode simulations
        geometry_base = [mp.Block(
            center=mp.Vector3(-Sx/2 + sim.waveguide_length/2 + sim.dpml/2, 0, 0),
            size=mp.Vector3(sim.waveguide_length + sim.dpml, geo.mm_waveguide_width, geo.waveguide_height),
            material=Si,
        )]
        for port_idx in range(3):
            py = self._port_y(port_idx)
            pw = self._port_width(port_idx)
            geometry_base.append(mp.Block(
                center=mp.Vector3(Sx/2 - sim.waveguide_length/2 - sim.dpml/2, py, 0),
                size=mp.Vector3(sim.waveguide_length + sim.dpml, pw, geo.waveguide_height),
                material=Si,
            ))

        mmi = mp.get_GDSII_prisms(Si, "gds/temp_gds_fdtd.gds", 2,
                                   -geo.waveguide_height/2, geo.waveguide_height/2)
        geometry_base += mmi

        powers_by_port = {}
        angles_by_port = {}

        for in_m_idx, fdtd_in_mode in enumerate(sim.input_modes):
            if mp.am_master(): self.master_logger.info(f"FDTD_GDS input mode {fdtd_in_mode}")
            sources = [mp.EigenModeSource(
                mp.GaussianSource(frequency=1/wl_c, fwidth=0.4/wl_c),
                eig_band=fdtd_in_mode,
                direction=mp.NO_DIRECTION,
                eig_kpoint=mp.Vector3(x=1),
                size=mp.Vector3(0, sim.oxide_width/2, sim.oxide_height),
                center=mp.Vector3(-Sx/2 + sim.dpml + sim.source_position, 0, 0),
            )]

            sim_obj = mp.Simulation(
                cell_size=mp.Vector3(Sx, Sy, Sz),
                boundary_layers=[mp.PML(sim.dpml)],
                geometry=geometry_base,
                sources=sources,
                eps_averaging=True,
                resolution=sim.res,
                default_material=SiO2,
            )

            in_flux = sim_obj.add_mode_monitor(freqs, mp.FluxRegion(
                center=mp.Vector3(-Sx/2 + sim.dpml + sim.sm_monitor_position, 0, 0),
                size=mp.Vector3(0, sim.oxide_width/2, Sz - 2*sim.dpml),
            ))
            out_fluxes = []
            for port_idx in range(3):
                py = self._port_y(port_idx)
                out_fluxes.append(sim_obj.add_mode_monitor(freqs, mp.FluxRegion(
                    center=mp.Vector3(Sx/2 - sim.dpml - sim.sm_monitor_position, py, 0),
                    size=mp.Vector3(0, sim.oxide_width/2, Sz - 2*sim.dpml),
                )))

            if in_m_idx == 0:
                if mp.am_master(): self.master_logger.info("Plotting FDTD_GDS geometry")
                plt.figure(); sim_obj.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, 0)))
                plt.savefig(f'{self.hp.outputfolder}/layout/FDTD_gds_xy.png', dpi=300); plt.close()

            if mp.am_master(): self.master_logger.info(f"Running FDTD_GDS (mode {fdtd_in_mode})")
            sim_obj.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol=1e-11))

            src_coeff = sim_obj.get_eigenmode_coefficients(
                in_flux, [fdtd_in_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))

            for port_idx, flux in enumerate(out_fluxes):
                for out_m in sim.output_modes:
                    coeff = sim_obj.get_eigenmode_coefficients(
                        flux, [out_m], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
                    s = coeff.alpha[:, :, 0] / src_coeff.alpha[:, :, 0]
                    s = s[0, :]
                    key = (fdtd_in_mode, port_idx, out_m)
                    pwr, ang = [], []
                    for pt in range(len(wavelengths)):
                        pwr.append(np.abs(s[pt])**2)
                        ang.append(np.angle(s[pt], deg=True))
                    powers_by_port[key] = pwr
                    angles_by_port[key] = ang

        data = {
            'powers_by_port': powers_by_port, 'angles_by_port': angles_by_port,
            'wavelengths': wavelengths,
        }
        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data, f)
        return powers_by_port, wavelengths

    # ------------------------------------------------------------------
    # Full EME simulation
    # ------------------------------------------------------------------

    def eme_simulation(self):
        wavelengths = self._wavelengths()
        geo, sim = self.hp.geo, self.hp.sim
        N = sim.num_modes_MMI
        Nt = sim.num_cells_taper
        n_in  = max(sim.input_modes)
        n_out = max(sim.output_modes)

        extrusion_rules = stack_extrusions_GF(geo.waveguide_height, sim.oxide_height)

        taper_L = geo.taper_length / Nt
        Ls = (
            [sim.waveguide_length]
            + [taper_L] * Nt        # input taper
            + [geo.mm_1_length]
            + [taper_L] * Nt        # MMI 1 output tapers
            + [geo.between_waveguide_length]
            + [taper_L] * Nt        # MMI 2 input tapers
            + [geo.mm_2_length]
            + [taper_L] * Nt        # MMI 2 output tapers
            + [sim.waveguide_length]
        )
        num_modes = (
            [n_in]
            + [sim.num_modes_taper] * Nt
            + [N]
            + [sim.num_modes_taper] * Nt
            + [sim.num_modes_taper]
            + [sim.num_modes_taper] * Nt
            + [N]
            + [sim.num_modes_taper] * Nt
            + [n_out]
        )
        mesh = mw.Mesh2D(
            x=np.linspace(-sim.oxide_width/2,  sim.oxide_width/2,  int(sim.res * sim.oxide_width)  + 1),
            y=np.linspace(-sim.oxide_height/2, sim.oxide_height/2, int(sim.res * sim.oxide_height) + 1),
        )

        # Per-port cells: one device GDS per output port (like 2x2 top/bot), so each
        # output cross-section contains only one waveguide and the S-matrix gives clean
        # port-specific transmission.
        cells_per_port = []
        for port_idx in range(3):
            port_gds = load_MMI_MUX_port(self.hp, port_idx)
            if port_idx == 0:
                plt.figure(); port_gds.plot()
                plt.savefig(f'{self.hp.outputfolder}/layout/eme_gds_full.png', dpi=300); plt.close()
            structs = mw.extrude_gds(port_gds, extrusion_rules)
            cells = mw.create_cells(structures=structs, mesh=mesh, Ls=Ls)
            cells_per_port.append(cells)

        keys = [(in_m, port_idx, out_m)
                for in_m in sim.input_modes
                for port_idx in range(3)
                for out_m in sim.output_modes]
        powers_by_port = {k: [] for k in keys}
        angles_by_port = {k: [] for k in keys}

        for wl in wavelengths:
            if mp.am_master(): self.master_logger.info(f"EME wl={wl:.4f}")
            env = mw.Environment(wl=float(wl), T=25.0)

            # Compute modes once for port 0 (intermediate cells are the same structure for all ports)
            css_base = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells_per_port[0]]
            for cs_num, cs in enumerate(css_base):
                plt.figure(); mw.visualize(cs, show=False)
                plt.savefig(f'{self.hp.outputfolder}/layout/cross_sections/cs_{cs_num}.png', dpi=300)
                plt.close()
            modes_base = [compute_modes(cs, num_modes=nm) for cs, nm in zip(css_base, num_modes)]

            for port_idx in range(3):
                if port_idx == 0:
                    modes_port = modes_base
                    cells_port = cells_per_port[0]
                else:
                    # Only swap the last cell mode (output WG differs per port)
                    cs_last = mw.CrossSection.from_cell(cell=cells_per_port[port_idx][-1], env=env)
                    last_mode = compute_modes(cs_last, num_modes=n_out)
                    modes_port = modes_base[:-1] + [last_mode]
                    cells_port = cells_per_port[port_idx]

                S, _ = mw.compute_s_matrix(modes_port, cells=cells_port)

                for in_m in sim.input_modes:
                    for out_m in sim.output_modes:
                        # S[i, n_in+j] = backward coupling (output j → input i)
                        # = forward transmission (input i → output j) by reciprocity
                        s_val = S[in_m - 1, n_in + (out_m - 1)]
                        key = (in_m, port_idx, out_m)
                        powers_by_port[key].append(np.abs(s_val)**2)
                        angles_by_port[key].append(np.angle(s_val, deg=True))
                        if mp.am_master():
                            self.master_logger.info(
                                f"in={in_m} port={PORT_DEFS[port_idx][0]} out={out_m} "
                                f"power={powers_by_port[key][-1]:.4f}")

            if float(wl) == float(wavelengths[len(wavelengths) // 2]):
                for in_m in sim.input_modes:
                    plot_eme_fields(self.hp, cells_per_port[0], modes_base, tag=in_m, input_mode=in_m)

        data = {
            'powers_by_port': powers_by_port, 'angles_by_port': angles_by_port,
            'wavelengths': wavelengths,
        }
        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data, f)
        return powers_by_port, wavelengths

    # ------------------------------------------------------------------
    # run() dispatcher + plotting
    # ------------------------------------------------------------------

    def run(self):
        if self.hp.sim.sim_type == "HYBRID":
            powers_by_port, wavelengths = self.hybrid_simulation()
        elif self.hp.sim.sim_type == "FDTD":
            powers_by_port, wavelengths = self.fdtd_simulation()
        elif self.hp.sim.sim_type == "FDTD_GDS":
            powers_by_port, wavelengths = self.fdtd_simulation_from_gds()
        elif self.hp.sim.sim_type == "EME":
            powers_by_port, wavelengths = self.eme_simulation()
        else:
            raise ValueError(f"Unknown sim_type: {self.hp.sim.sim_type}")

        for dB, suffix in [(False, ""), (True, "_dB")]:
            plt.figure()
            for (in_m, port_idx, out_m), powers in powers_by_port.items():
                port_label = PORT_DEFS[port_idx][0]
                mux_tag = " [MUX out]" if port_idx in (1, 2) else ""
                label = f"In {in_m} → {port_label} out{out_m}{mux_tag}"
                y = 10 * np.log10(np.array(powers) + 1e-30) if dB else powers
                plt.plot(wavelengths, y, label=label)
            plt.xlabel("Wavelength (µm)")
            plt.ylabel("Transmission (dB)" if dB else "Transmission")
            plt.legend()
            plt.savefig(f'{self.hp.outputfolder}/performance{suffix}.png', dpi=300)
            plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file", required=True)
    args = parser.parse_args()
    mmi = MMI_MUX(args)
    mmi.run()
