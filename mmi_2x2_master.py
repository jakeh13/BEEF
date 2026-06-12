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

from engine.hyperparams.MMI import HyperParamsMMI as HyperParams

from engine.helpers.general import(
    make_output_folders,
    resolve_output_folder,
)
from engine.helpers.meep import(
     meep_material_from_meow,
     create_temp_gds_for_mmi_2x2_input,
     create_temp_gds_for_mmi_2x2_output,
     create_temp_gds_for_full_mmi_2x2,
     prepare_gds_for_meep,
     create_bend_gds_for_meep,
)
from engine.helpers.meow import(
    load_2x2_MMI_top,
    load_2x2_MMI_bot,
    stack_extrusions_GF
)
from engine.plotting.meow import(
    plot_eme_fields
)

class MMI_2x2:
    def __init__(self,args):
        self.config = args.config
        self.hp: HyperParams = dataconf.load(self.config, HyperParams)
        if mp.am_master():
            self.hp.outputfolder = resolve_output_folder(self.hp.outputfolder)
            make_output_folders(self.hp.outputfolder)
        try:
            from mpi4py import MPI
            self.hp.outputfolder = MPI.COMM_WORLD.bcast(self.hp.outputfolder, root=0)
        except ImportError:
            pass

        # save config in output folder
        s = json.dumps(asdict(self.hp), indent=4, default=lambda x: x.name)
        with open(f"{self.hp.outputfolder}/config.json", "w") as f:
            print(s, file=f)

        # Create logger
        handler = logging.FileHandler(f"{self.hp.outputfolder}/logs/master.log")
        handler.setLevel("INFO")

        handler.setFormatter(
            logging.Formatter(
                "%(asctime)s | %(levelname)-4s | %(filename)s | %(funcName)s:%(lineno)04d | %(message)s"
            )
        )
        self.master_logger = logging.getLogger(__name__.replace("_", ""))
        self.master_logger.addHandler(handler)
        self.master_logger.setLevel("INFO")

    def top_simulation(self, wavelengths, mode: int = 1):
        if mp.am_master(): self.master_logger.info(f"Creating top FDTD simulation (mode {mode})")
        # Compute the simulation dimensions
        Sx = self.hp.geo.taper_length + self.hp.sim.waveguide_length + self.hp.sim.mm_sim_length + 2*self.hp.sim.dpml  
        Sy = self.hp.sim.oxide_width + 2*self.hp.sim.dpml  
        Sz = self.hp.sim.oxide_height + 2*self.hp.sim.dpml  
        wavelength_center = (self.hp.sim.wavelength_max - self.hp.sim.wavelength_min)/2 + self.hp.sim.wavelength_min
        freqs = 1/wavelengths.astype(list)

        # materials
        Si = meep_material_from_meow(mw.silicon.n, mw.silicon.params['wl'],wavelength_center)
        SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'],wavelength_center)

        # Set waveguide geometry
        geometry = [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2,self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry += [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2,-self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]

        # Create temp gds file for MMI geometry (master only, then sync)
        if mp.am_master():
            create_temp_gds_for_mmi_2x2_input(self.hp.geo.mm_length, self.hp.geo.taper_length, Sx, self.hp.sim.dpml, self.hp.sim.waveguide_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset, f"top_{mode}")
        mp.all_wait()

        # Get MMI geometry from temp gds file
        mmi = mp.get_GDSII_prisms(Si, f"gds/top_{mode}.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry += mmi

        # Create source MEEP simulation
        source_gaus = [
            mp.EigenModeSource(
                mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                eig_band=mode,
                direction=mp.NO_DIRECTION,
                eig_kpoint=mp.Vector3(x=1),
                size=mp.Vector3(0, self.hp.sim.oxide_width/2, self.hp.sim.oxide_height),
                center=[-Sx/2+self.hp.sim.dpml+self.hp.sim.source_position, self.hp.sim.oxide_width/4, 0],
            )
        ]

        # Create MEEP simulation
        sim = mp.Simulation(
            cell_size=mp.Vector3(Sx, Sy, Sz),
            boundary_layers=[mp.PML(self.hp.sim.dpml)],
            geometry=geometry,
            sources=source_gaus,
            eps_averaging=True,
            resolution=self.hp.sim.res,
            default_material=SiO2
        )

        # Add monitors
        in_flux = sim.add_mode_monitor(
            freqs,
            mp.FluxRegion(
                center = mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position, self.hp.sim.oxide_width/4, 0),
                size = mp.Vector3(0, self.hp.sim.oxide_width/2,Sz - 2*self.hp.sim.dpml)
            )
        )

        out_flux = sim.add_mode_monitor(
            freqs,
            mp.FluxRegion(
                center = mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.waveguide_length + self.hp.geo.taper_length + self.hp.sim.mm_monitor_position, 0, 0),
                size = mp.Vector3(0, Sy - 2*self.hp.sim.dpml,Sz - 2*self.hp.sim.dpml)
            )
        )

        # Plot simulation geometry (only for mode 1 to avoid duplicate layout plots)
        if mode == 1:
            if mp.am_master(): self.master_logger.info("Plotting top FDTD simulation geometry")
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, Sy, 0)))
            plt.savefig(f'{self.hp.outputfolder}/layout/top_xy.png', dpi=600)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/top_xz.png', dpi=600)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/top_yz.png', dpi=600)

        # run simulation
        if mp.am_master(): self.master_logger.info(f"Running top FDTD simulation (mode {mode})")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

        # get eigenmode coefficients
        if mp.am_master(): self.master_logger.info(f"Getting top FDTD simulation eigenmode coefficients (mode {mode}):")
        src_coeff = sim.get_eigenmode_coefficients(in_flux,[mode],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
        top_coeffs = []
        for mode_num in range(1,self.hp.sim.num_modes_MMI+1):
            coeff = sim.get_eigenmode_coefficients(out_flux,[mode_num],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
            s = coeff.alpha[:,:,0]/src_coeff.alpha[:,:,0]
            top_coeffs.append(s[0,:])
            if mp.am_master(): self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(s[0,:])**2}, Angle: {np.angle(s[0,:], deg=True)}")
        top_coeffs = np.asarray(top_coeffs)
        return top_coeffs.T.astype(list)

    def bot_simulation(self, wavelengths, mode: int = 1):
        if mp.am_master(): self.master_logger.info(f"Creating bottom FDTD simulation (mode {mode})")
        # Compute the simulation dimensions
        Sx = self.hp.geo.taper_length + self.hp.sim.waveguide_length + self.hp.sim.mm_sim_length + 2*self.hp.sim.dpml  
        Sy = self.hp.sim.oxide_width + 2*self.hp.sim.dpml  
        Sz = self.hp.sim.oxide_height + 2*self.hp.sim.dpml  
        wavelength_center = (self.hp.sim.wavelength_max - self.hp.sim.wavelength_min)/2 + self.hp.sim.wavelength_min
        freqs = 1/wavelengths.astype(list)

        # materials
        Si = meep_material_from_meow(mw.silicon.n, mw.silicon.params['wl'],wavelength_center)
        SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'],wavelength_center)

        # Set waveguide geometry
        geometry = [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2,self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry += [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2,-self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]

        # Create temp gds file for MMI geometry (master only, then sync)
        if mp.am_master():
            create_temp_gds_for_mmi_2x2_input(self.hp.geo.mm_length, self.hp.geo.taper_length, Sx, self.hp.sim.dpml, self.hp.sim.waveguide_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset, f"bot_{mode}")
        mp.all_wait()

        # Get MMI geometry from temp gds file
        mmi = mp.get_GDSII_prisms(Si, f"gds/bot_{mode}.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry += mmi

        # Create source MEEP simulation
        source_gaus = [
            mp.EigenModeSource(
                mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                eig_band=mode,
                direction=mp.NO_DIRECTION,
                eig_kpoint=mp.Vector3(x=1),
                size=mp.Vector3(0, self.hp.sim.oxide_width/2, self.hp.sim.oxide_height),
                center=[-Sx/2+self.hp.sim.dpml+self.hp.sim.source_position, -self.hp.sim.oxide_width/4, 0],
            )
        ]

        # Create MEEP simulation
        sim = mp.Simulation(
            cell_size=mp.Vector3(Sx, Sy, Sz),
            boundary_layers=[mp.PML(self.hp.sim.dpml)],
            geometry=geometry,
            sources=source_gaus,
            eps_averaging=True,
            resolution=self.hp.sim.res,
            default_material=SiO2
        )

        # Add monitors
        in_flux = sim.add_mode_monitor(
            freqs,
            mp.FluxRegion(
                center = mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position, -self.hp.sim.oxide_width/4, 0),
                size = mp.Vector3(0, self.hp.sim.oxide_width/2,Sz - 2*self.hp.sim.dpml)
            )
        )

        out_flux = sim.add_mode_monitor(
            freqs,
            mp.FluxRegion(
                center = mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.waveguide_length + self.hp.geo.taper_length + self.hp.sim.mm_monitor_position, 0, 0),
                size = mp.Vector3(0, Sy - 2*self.hp.sim.dpml,Sz - 2*self.hp.sim.dpml)
            )
        )

        # Plot simulation geometry (only for mode 1 to avoid duplicate layout plots)
        if mode == 1:
            if mp.am_master(): self.master_logger.info("Plotting bottom FDTD simulation geometry")
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, Sy, 0)))
            plt.savefig(f'{self.hp.outputfolder}/layout/bot_xy.png', dpi=600)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/bot_xz.png', dpi=600)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/bot_yz.png', dpi=600)

        # run simulation
        if mp.am_master(): self.master_logger.info(f"Running bottom FDTD simulation (mode {mode})")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

        # get eigenmode coefficients
        if mp.am_master(): self.master_logger.info(f"Getting bottom FDTD simulation eigenmode coefficients (mode {mode}):")
        src_coeff = sim.get_eigenmode_coefficients(in_flux,[mode],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
        top_coeffs = []
        for mode_num in range(1,self.hp.sim.num_modes_MMI+1):
            coeff = sim.get_eigenmode_coefficients(out_flux,[mode_num],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
            s = coeff.alpha[:,:,0]/src_coeff.alpha[:,:,0]
            top_coeffs.append(s[0,:])
            if mp.am_master(): self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(s[0,:])**2}, Angle: {np.angle(s[0,:], deg=True)}")
        top_coeffs = np.asarray(top_coeffs)
        return top_coeffs.T.astype(list)
    
    def mm_simulation(self,wavelengths):
        mm_coeffs = []
        # Load the GDSII file for the MMI splitter
        mmi_top = load_2x2_MMI_top(self.hp.geo.mm_length,0.0,self.hp.geo.waveguide_width,self.hp.geo.waveguide_offset,self.hp.sim.oxide_width,-self.hp.geo.taper_length-self.hp.sim.mm_monitor_position,self.hp.geo.taper_length,self.hp.geo.mm_width,self.hp.geo.taper_width)
        extrusion_rules = stack_extrusions_GF(self.hp.geo.waveguide_height, self.hp.sim.oxide_height)
        structs = mw.extrude_gds(mmi_top, extrusion_rules)            
        # EME cell for the MMI region
        Ls = [self.hp.geo.mm_length - 2*self.hp.sim.mm_monitor_position]
        num_modes = [self.hp.sim.num_modes_MMI]
        cells = mw.create_cells(
            structures=structs,
            mesh=mw.Mesh2D(
                x=np.linspace(-self.hp.sim.oxide_width/2, self.hp.sim.oxide_width/2, int(self.hp.sim.res*self.hp.sim.oxide_width) + 1),
                y=np.linspace(-self.hp.sim.oxide_height/2, self.hp.sim.oxide_height/2, int(self.hp.sim.res*self.hp.sim.oxide_height) + 1),
            ),
            Ls=Ls,
        )
        for wavelength in wavelengths.astype(list):
            if mp.am_master(): self.master_logger.info(f"Starting EME simulation for multimode section at wavelength {wavelength}")
            
            gf_environment = mw.Environment(wl=wavelength, T=25.0)
            css = [mw.CrossSection.from_cell(cell=cell, env=gf_environment) for cell in cells]

            # Performe mode simulations
            if mp.am_master(): self.master_logger.info(f"Performing mode simulations")
            modes = []
            for cell_num, (cs, num_mode) in enumerate(zip(css,num_modes)):
                if mp.am_master(): self.master_logger.info(f"Computing mode {num_mode} of cell {cell_num}")
                modes.append(compute_modes(cs, num_modes=num_mode))

            S, port_map = mw.compute_s_matrix(modes, cells=cells)
            mm_coeffs_1_wl = []
            for mode_num in range(self.hp.sim.num_modes_MMI):
                mm_coeffs_1_wl.append(S[mode_num,self.hp.sim.num_modes_MMI+mode_num])
                if mp.am_master(): self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(S[mode_num,self.hp.sim.num_modes_MMI+mode_num])**2}, Angle: {np.angle(S[mode_num,self.hp.sim.num_modes_MMI+mode_num], deg=True)}")
            mm_coeffs.append(mm_coeffs_1_wl)
        return mm_coeffs

    def hybrid_simulation(self):
        wavelengths = np.linspace(self.hp.sim.wavelength_min, self.hp.sim.wavelength_max, self.hp.sim.wavelength_pts)
        input_modes = self.hp.sim.input_modes
        output_modes = self.hp.sim.output_modes

        # Minimum set of FDTD simulations needed:
        #   top: used for both input excitation and output-port overlap (by reciprocity)
        #   bot: used only for output-port overlap
        modes_needing_top = set(input_modes) | set(output_modes)
        modes_needing_bot = set(output_modes)

        all_top_coeffs = {}
        for m in modes_needing_top:
            all_top_coeffs[m] = self.top_simulation(wavelengths, m)

        all_bot_coeffs = {}
        for m in modes_needing_bot:
            all_bot_coeffs[m] = self.bot_simulation(wavelengths, m)

        mm_coeffs = self.mm_simulation(wavelengths)

        keys = [(in_m, out_m) for in_m in input_modes for out_m in output_modes]
        powers_top_by_mode = {k: [] for k in keys}
        angles_top_by_mode = {k: [] for k in keys}
        powers_bot_by_mode = {k: [] for k in keys}
        angles_bot_by_mode = {k: [] for k in keys}

        for point in range(len(wavelengths)):
            mm_coeff = mm_coeffs[point]
            for in_m in input_modes:
                input_coeff = all_top_coeffs[in_m][point]
                propagated_coeffs = [x * y for x, y in zip(input_coeff, mm_coeff)]
                for out_m in output_modes:
                    top_coeff = all_top_coeffs[out_m][point]
                    bot_coeff = all_bot_coeffs[out_m][point]
                    out_top = sum(x * y for x, y in zip(propagated_coeffs, top_coeff))
                    out_bot = sum(x * y for x, y in zip(propagated_coeffs, bot_coeff))

                    powers_top_by_mode[(in_m, out_m)].append(np.abs(out_top)**2)
                    angles_top_by_mode[(in_m, out_m)].append(np.angle(out_top, deg=True))
                    powers_bot_by_mode[(in_m, out_m)].append(np.abs(out_bot)**2)
                    angles_bot_by_mode[(in_m, out_m)].append(np.angle(out_bot, deg=True))

                    if mp.am_master(): self.master_logger.info(f"In {in_m} Out {out_m} | Top: Power {powers_top_by_mode[(in_m, out_m)][-1]:.4f}, Angle {angles_top_by_mode[(in_m, out_m)][-1]:.2f}")
                    if mp.am_master(): self.master_logger.info(f"In {in_m} Out {out_m} | Bot: Power {powers_bot_by_mode[(in_m, out_m)][-1]:.4f}, Angle {angles_bot_by_mode[(in_m, out_m)][-1]:.2f}")

        data_to_save = {
            'all_top_coeffs': all_top_coeffs,
            'all_bot_coeffs': all_bot_coeffs,
            'mm_coeffs': mm_coeffs,
            'powers_top_by_mode': powers_top_by_mode,
            'angles_top_by_mode': angles_top_by_mode,
            'powers_bot_by_mode': powers_bot_by_mode,
            'angles_bot_by_mode': angles_bot_by_mode,
            'wavelengths': wavelengths,
        }

        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)

        return powers_top_by_mode, powers_bot_by_mode, wavelengths

    def fdtd_simulation(self):
        if mp.am_master(): self.master_logger.info(f"Creating full FDTD simulation")
        # Compute the simulation dimensions
        Sx = 2*self.hp.geo.taper_length + 2*self.hp.sim.waveguide_length + self.hp.geo.mm_length + 2*self.hp.sim.dpml
        Sy = self.hp.sim.oxide_width + 2*self.hp.sim.dpml
        Sz = self.hp.sim.oxide_height + 2*self.hp.sim.dpml
        wavelength_center = (self.hp.sim.wavelength_max - self.hp.sim.wavelength_min)/2 + self.hp.sim.wavelength_min
        wavelengths = np.linspace(self.hp.sim.wavelength_min,self.hp.sim.wavelength_max,self.hp.sim.wavelength_pts)
        freqs = 1/wavelengths.astype(list)

        # materials
        Si = meep_material_from_meow(mw.silicon.n, mw.silicon.params['wl'],wavelength_center)
        SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'],wavelength_center)

        # Build geometry once — shared across all input mode simulations
        geometry_base = [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2,self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry_base += [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2,-self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry_base += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2,self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry_base += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2,-self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]

        # Create temp gds file for MMI geometry (master only, then sync)
        if mp.am_master():
            create_temp_gds_for_full_mmi_2x2(self.hp.geo.mm_length, self.hp.geo.taper_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset)
        mp.all_wait()

        # Get MMI geometry from temp gds file
        mmi = mp.get_GDSII_prisms(Si, "gds/temp_full.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry_base += mmi

        coeffs_top_by_mode = {}
        coeffs_bot_by_mode = {}
        powers_top_by_mode = {}
        angles_top_by_mode = {}
        powers_bot_by_mode = {}
        angles_bot_by_mode = {}

        _port_y_map = {"top": self.hp.sim.oxide_width / 4, "bot": -self.hp.sim.oxide_width / 4}
        for port_label in self.hp.sim.input_ports:
            y_port = _port_y_map[port_label]
            for in_m_idx, fdtd_input_mode in enumerate(self.hp.sim.input_modes):
                if mp.am_master(): self.master_logger.info(f"Full FDTD port={port_label} mode={fdtd_input_mode}")

                source_gaus = [
                    mp.EigenModeSource(
                        mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                        eig_band=fdtd_input_mode,
                        direction=mp.NO_DIRECTION,
                        eig_kpoint=mp.Vector3(x=1),
                        size=mp.Vector3(0, self.hp.sim.oxide_width/2, self.hp.sim.oxide_height),
                        center=[-Sx/2+self.hp.sim.dpml+self.hp.sim.source_position, y_port, 0],
                    )
                ]

                sim = mp.Simulation(
                    cell_size=mp.Vector3(Sx, Sy, Sz),
                    boundary_layers=[mp.PML(self.hp.sim.dpml)],
                    geometry=geometry_base,
                    sources=source_gaus,
                    eps_averaging=True,
                    resolution=self.hp.sim.res,
                    default_material=SiO2
                )

                in_flux = sim.add_mode_monitor(
                    freqs,
                    mp.FluxRegion(
                        center = mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position, y_port, 0),
                        size = mp.Vector3(0, self.hp.sim.oxide_width/2,Sz - 2*self.hp.sim.dpml)
                    )
                )

                out_flux_top = sim.add_mode_monitor(
                    freqs,
                    mp.FluxRegion(
                        center = mp.Vector3(Sx/2 - self.hp.sim.dpml - self.hp.sim.sm_monitor_position, (Sy - 2*self.hp.sim.dpml)/4, 0),
                        size = mp.Vector3(0, (Sy - 2*self.hp.sim.dpml)/2,Sz - 2*self.hp.sim.dpml)
                    )
                )

                out_flux_bot = sim.add_mode_monitor(
                    freqs,
                    mp.FluxRegion(
                        center = mp.Vector3(Sx/2 - self.hp.sim.dpml - self.hp.sim.sm_monitor_position, -(Sy - 2*self.hp.sim.dpml)/4, 0),
                        size = mp.Vector3(0, (Sy - 2*self.hp.sim.dpml)/2,Sz - 2*self.hp.sim.dpml)
                    )
                )

                if in_m_idx == 0:
                    if mp.am_master(): self.master_logger.info(f"Plotting FDTD simulation geometry (port={port_label})")
                    plt.figure()
                    sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, Sy, 0)))
                    plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_xy_{port_label}.png', dpi=600)
                    plt.figure()
                    sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
                    plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_xz_{port_label}.png', dpi=600)
                    plt.figure()
                    sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
                    plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_yz_{port_label}.png', dpi=600)

                if mp.am_master(): self.master_logger.info(f"Running FDTD simulation (port={port_label} mode={fdtd_input_mode})")
                sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

                if mp.am_master(): self.master_logger.info(f"Getting FDTD eigenmode coefficients (port={port_label} mode={fdtd_input_mode}):")
                src_coeff = sim.get_eigenmode_coefficients(in_flux, [fdtd_input_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))

                for out_mode in self.hp.sim.output_modes:
                    coeff_top = sim.get_eigenmode_coefficients(out_flux_top, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
                    coeff_bot = sim.get_eigenmode_coefficients(out_flux_bot, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
                    s_top = coeff_top.alpha[:,:,0] / src_coeff.alpha[:,:,0]
                    s_top = s_top[0,:]
                    s_bot = coeff_bot.alpha[:,:,0] / src_coeff.alpha[:,:,0]
                    s_bot = s_bot[0,:]
                    key = (port_label, fdtd_input_mode, out_mode)
                    coeffs_top_by_mode[key] = s_top
                    coeffs_bot_by_mode[key] = s_bot
                    pwr_top, ang_top, pwr_bot, ang_bot = [], [], [], []
                    for point in range(len(wavelengths)):
                        pwr_top.append(np.abs(s_top[point])**2)
                        ang_top.append(np.angle(s_top[point], deg=True))
                        pwr_bot.append(np.abs(s_bot[point])**2)
                        ang_bot.append(np.angle(s_bot[point], deg=True))
                        if mp.am_master(): self.master_logger.info(f"Port {port_label} In {fdtd_input_mode} Out {out_mode} | Top: Power {pwr_top[-1]:.4f}, Angle {ang_top[-1]:.2f}")
                        if mp.am_master(): self.master_logger.info(f"Port {port_label} In {fdtd_input_mode} Out {out_mode} | Bot: Power {pwr_bot[-1]:.4f}, Angle {ang_bot[-1]:.2f}")
                    powers_top_by_mode[key] = pwr_top
                    angles_top_by_mode[key] = ang_top
                    powers_bot_by_mode[key] = pwr_bot
                    angles_bot_by_mode[key] = ang_bot

        data_to_save = {
            'coeffs_top_by_mode': coeffs_top_by_mode,
            'coeffs_bot_by_mode': coeffs_bot_by_mode,
            'powers_top_by_mode': powers_top_by_mode,
            'angles_top_by_mode': angles_top_by_mode,
            'powers_bot_by_mode': powers_bot_by_mode,
            'angles_bot_by_mode': angles_bot_by_mode,
            'wavelengths': wavelengths,
        }

        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)

        return powers_top_by_mode, powers_bot_by_mode, wavelengths

    def fdtd_simulation_from_gds(self):
        if mp.am_master(): self.master_logger.info(f"Creating full FDTD simulation from GDS: {self.hp.geo.gds_file}")
        wavelength_center = (self.hp.sim.wavelength_max - self.hp.sim.wavelength_min)/2 + self.hp.sim.wavelength_min
        wavelengths = np.linspace(self.hp.sim.wavelength_min, self.hp.sim.wavelength_max, self.hp.sim.wavelength_pts)
        freqs = 1/wavelengths.astype(list)

        # Load and center the GDS file; derive simulation x-extent from its bounding box
        if mp.am_master():
            device_Sx = prepare_gds_for_meep(self.hp.geo.gds_file, "gds/temp_gds_fdtd.gds")
        else:
            device_Sx = None
        try:
            from mpi4py import MPI
            device_Sx = MPI.COMM_WORLD.bcast(device_Sx, root=0)
        except ImportError:
            pass
        mp.all_wait()
        Sx = device_Sx + 2*self.hp.sim.waveguide_length + 2*self.hp.sim.dpml
        Sy = self.hp.sim.oxide_width + 2*self.hp.sim.dpml
        Sz = self.hp.sim.oxide_height + 2*self.hp.sim.dpml

        # materials
        Si = meep_material_from_meow(mw.silicon.n, mw.silicon.params['wl'], wavelength_center)
        SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'], wavelength_center)

        # Build geometry once — shared across all input mode simulations
        geometry_base = [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2, self.hp.geo.waveguide_offset/2, 0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material=Si)]
        geometry_base += [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2, -self.hp.geo.waveguide_offset/2, 0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material=Si)]
        geometry_base += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2, self.hp.geo.waveguide_offset/2, 0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material=Si)]
        geometry_base += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2, -self.hp.geo.waveguide_offset/2, 0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material=Si)]

        # Get MMI geometry from the centered GDS file
        mmi = mp.get_GDSII_prisms(Si, "gds/temp_gds_fdtd.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry_base += mmi

        coeffs_top_by_mode = {}
        coeffs_bot_by_mode = {}
        powers_top_by_mode = {}
        angles_top_by_mode = {}
        powers_bot_by_mode = {}
        angles_bot_by_mode = {}

        _port_y_map = {"top": self.hp.sim.oxide_width / 4, "bot": -self.hp.sim.oxide_width / 4}
        for port_label in self.hp.sim.input_ports:
            y_port = _port_y_map[port_label]
            for in_m_idx, fdtd_input_mode in enumerate(self.hp.sim.input_modes):
                if mp.am_master(): self.master_logger.info(f"FDTD_GDS port={port_label} mode={fdtd_input_mode}")

                source_gaus = [
                    mp.EigenModeSource(
                        mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                        eig_band=fdtd_input_mode,
                        direction=mp.NO_DIRECTION,
                        eig_kpoint=mp.Vector3(x=1),
                        size=mp.Vector3(0, self.hp.sim.oxide_width/2, self.hp.sim.oxide_height),
                        center=[-Sx/2+self.hp.sim.dpml+self.hp.sim.source_position, y_port, 0],
                    )
                ]

                sim = mp.Simulation(
                    cell_size=mp.Vector3(Sx, Sy, Sz),
                    boundary_layers=[mp.PML(self.hp.sim.dpml)],
                    geometry=geometry_base,
                    sources=source_gaus,
                    eps_averaging=True,
                    resolution=self.hp.sim.res,
                    default_material=SiO2
                )

                in_flux = sim.add_mode_monitor(
                    freqs,
                    mp.FluxRegion(
                        center=mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position, y_port, 0),
                        size=mp.Vector3(0, self.hp.sim.oxide_width/2, Sz - 2*self.hp.sim.dpml)
                    )
                )

                out_flux_top = sim.add_mode_monitor(
                    freqs,
                    mp.FluxRegion(
                        center=mp.Vector3(Sx/2 - self.hp.sim.dpml - self.hp.sim.sm_monitor_position, (Sy - 2*self.hp.sim.dpml)/4, 0),
                        size=mp.Vector3(0, (Sy - 2*self.hp.sim.dpml)/2, Sz - 2*self.hp.sim.dpml)
                    )
                )

                out_flux_bot = sim.add_mode_monitor(
                    freqs,
                    mp.FluxRegion(
                        center=mp.Vector3(Sx/2 - self.hp.sim.dpml - self.hp.sim.sm_monitor_position, -(Sy - 2*self.hp.sim.dpml)/4, 0),
                        size=mp.Vector3(0, (Sy - 2*self.hp.sim.dpml)/2, Sz - 2*self.hp.sim.dpml)
                    )
                )

                if in_m_idx == 0:
                    if mp.am_master(): self.master_logger.info(f"Plotting FDTD GDS simulation geometry (port={port_label})")
                    plt.figure()
                    sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, 0)))
                    plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_gds_xy_{port_label}.png', dpi=600)
                    plt.figure()
                    sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, 0, Sz)))
                    plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_gds_xz_{port_label}.png', dpi=600)
                    plt.figure()
                    sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(0, Sy, Sz)))
                    plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_gds_yz_{port_label}.png', dpi=600)

                if mp.am_master(): self.master_logger.info(f"Running FDTD GDS (port={port_label} mode={fdtd_input_mode})")
                sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol=1e-11))

                if mp.am_master(): self.master_logger.info(f"Getting FDTD GDS eigenmode coefficients (port={port_label} mode={fdtd_input_mode}):")
                src_coeff = sim.get_eigenmode_coefficients(in_flux, [fdtd_input_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))

                for out_mode in self.hp.sim.output_modes:
                    coeff_top = sim.get_eigenmode_coefficients(out_flux_top, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
                    coeff_bot = sim.get_eigenmode_coefficients(out_flux_bot, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
                    s_top = coeff_top.alpha[:, :, 0] / src_coeff.alpha[:, :, 0]
                    s_top = s_top[0, :]
                    s_bot = coeff_bot.alpha[:, :, 0] / src_coeff.alpha[:, :, 0]
                    s_bot = s_bot[0, :]
                    key = (port_label, fdtd_input_mode, out_mode)
                    coeffs_top_by_mode[key] = s_top
                    coeffs_bot_by_mode[key] = s_bot
                    pwr_top, ang_top, pwr_bot, ang_bot = [], [], [], []
                    for point in range(len(wavelengths)):
                        pwr_top.append(np.abs(s_top[point])**2)
                        ang_top.append(np.angle(s_top[point], deg=True))
                        pwr_bot.append(np.abs(s_bot[point])**2)
                        ang_bot.append(np.angle(s_bot[point], deg=True))
                        if mp.am_master(): self.master_logger.info(f"Port {port_label} In {fdtd_input_mode} Out {out_mode} | Top: Power {pwr_top[-1]:.4f}, Angle {ang_top[-1]:.2f}")
                        if mp.am_master(): self.master_logger.info(f"Port {port_label} In {fdtd_input_mode} Out {out_mode} | Bot: Power {pwr_bot[-1]:.4f}, Angle {ang_bot[-1]:.2f}")
                    powers_top_by_mode[key] = pwr_top
                    angles_top_by_mode[key] = ang_top
                    powers_bot_by_mode[key] = pwr_bot
                    angles_bot_by_mode[key] = ang_bot

        data_to_save = {
            'coeffs_top_by_mode': coeffs_top_by_mode,
            'coeffs_bot_by_mode': coeffs_bot_by_mode,
            'powers_top_by_mode': powers_top_by_mode,
            'angles_top_by_mode': angles_top_by_mode,
            'powers_bot_by_mode': powers_bot_by_mode,
            'angles_bot_by_mode': angles_bot_by_mode,
            'wavelengths': wavelengths,
        }

        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)

        return powers_top_by_mode, powers_bot_by_mode, wavelengths

    def fdtd_simulation_from_gds_with_bends(self):
        if mp.am_master(): self.master_logger.info(f"Creating FDTD_GDS_BEND simulation from GDS: {self.hp.geo.gds_file}")
        wavelength_center = (self.hp.sim.wavelength_max - self.hp.sim.wavelength_min)/2 + self.hp.sim.wavelength_min
        wavelengths = np.linspace(self.hp.sim.wavelength_min, self.hp.sim.wavelength_max, self.hp.sim.wavelength_pts)
        freqs = 1/wavelengths.astype(list)

        if mp.am_master():
            device_Sx = prepare_gds_for_meep(self.hp.geo.gds_file, "gds/temp_gds_fdtd.gds")
            bend_info = create_bend_gds_for_meep(
                self.hp.geo.waveguide_width, self.hp.geo.waveguide_offset,
                self.hp.geo.bend_radius, self.hp.geo.bend_type, self.hp.geo.bend_locations,
                device_Sx, self.hp.geo.bend_stub_length, self.hp.geo.bend_wg_length,
                "gds/temp_gds_bends.gds",
                dpml=self.hp.sim.dpml,
                monitor_margin=self.hp.sim.oxide_height/2,
            )
        else:
            device_Sx = None
            bend_info = None
        try:
            from mpi4py import MPI
            device_Sx = MPI.COMM_WORLD.bcast(device_Sx, root=0)
            bend_info = MPI.COMM_WORLD.bcast(bend_info, root=0)
        except ImportError:
            pass
        mp.all_wait()
        x_extra_left  = bend_info["x_extra_left"]
        x_extra_right = bend_info["x_extra_right"]
        y_extra       = bend_info["y_extra"]
        R = self.hp.geo.bend_radius

        Sx = device_Sx + (self.hp.geo.bend_stub_length + x_extra_left) + (self.hp.geo.bend_stub_length + x_extra_right) + 2*self.hp.sim.dpml
        Sy = self.hp.geo.waveguide_offset + 2*y_extra + 2*self.hp.sim.dpml
        Sz = self.hp.sim.oxide_height + 2*self.hp.sim.dpml

        Si  = meep_material_from_meow(mw.silicon.n, mw.silicon.params['wl'], wavelength_center)
        SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'], wavelength_center)

        has_output_bends = "output" in self.hp.geo.bend_locations
        has_input_bends  = "input"  in self.hp.geo.bend_locations

        # Waveguide stubs: stop at the bend start when that side has a bend; extend
        # through x-PML when no bend (source/monitor live in the stub on that side).
        if has_input_bends:
            x_in_center = -Sx/2 + self.hp.sim.dpml + x_extra_left + self.hp.geo.bend_stub_length / 2
            in_stub_sz  = self.hp.geo.bend_stub_length
        else:
            x_in_center = -Sx/2 + (self.hp.geo.bend_stub_length + self.hp.sim.dpml) / 2
            in_stub_sz  = self.hp.geo.bend_stub_length + self.hp.sim.dpml

        if has_output_bends:
            x_out_center = Sx/2 - self.hp.sim.dpml - x_extra_right - self.hp.geo.bend_stub_length / 2
            out_stub_sz  = self.hp.geo.bend_stub_length
        else:
            x_out_center = Sx/2 - (self.hp.geo.bend_stub_length + self.hp.sim.dpml) / 2
            out_stub_sz  = self.hp.geo.bend_stub_length + self.hp.sim.dpml

        stub_size     = mp.Vector3(in_stub_sz,  self.hp.geo.waveguide_width, self.hp.geo.waveguide_height)
        stub_size_out = mp.Vector3(out_stub_sz, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height)

        geometry_base = [mp.Block(center=mp.Vector3(x_in_center,  self.hp.geo.waveguide_offset/2, 0), size=stub_size,     material=Si)]
        geometry_base += [mp.Block(center=mp.Vector3(x_in_center, -self.hp.geo.waveguide_offset/2, 0), size=stub_size,     material=Si)]
        geometry_base += [mp.Block(center=mp.Vector3(x_out_center, self.hp.geo.waveguide_offset/2, 0), size=stub_size_out, material=Si)]
        geometry_base += [mp.Block(center=mp.Vector3(x_out_center,-self.hp.geo.waveguide_offset/2, 0), size=stub_size_out, material=Si)]

        device_prisms = mp.get_GDSII_prisms(Si, "gds/temp_gds_fdtd.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry_base += device_prisms
        if self.hp.geo.bend_locations:
            bend_prisms = mp.get_GDSII_prisms(Si, "gds/temp_gds_bends.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
            geometry_base += bend_prisms

        # Taper monitor x-position: sm_monitor_position from the GDS right edge (middle of output stub)
        x_taper_mon = Sx/2 - self.hp.sim.dpml - x_extra_right - self.hp.geo.bend_stub_length + self.hp.sim.sm_monitor_position
        # Input taper monitor x-position: mirror of x_taper_mon (sm_monitor_position from device edge into input stub)
        x_in_taper_mon = -Sx/2 + self.hp.sim.dpml + x_extra_left + self.hp.geo.bend_stub_length - self.hp.sim.sm_monitor_position

        # Bend exit monitor positions (y-flux, only when output bends present)

        # x-center of bend exit (right-side bends go R past the end of the stub)
        x_bend_exit = Sx/2 - self.hp.sim.dpml - x_extra_right + R
        # Monitor sits sm_monitor_position inside the active region (= just at PML boundary)
        y_top_bend_mon =  Sy/2 - self.hp.sim.dpml
        y_bot_bend_mon = -(Sy/2 - self.hp.sim.dpml)

        # Input bend pre-bend monitor x-center and y positions (y-flux, left side)
        x_in_bend_exit = -Sx/2 + self.hp.sim.dpml + x_extra_left - R
        y_top_in_bend_mon =  Sy/2 - self.hp.sim.dpml
        y_bot_in_bend_mon = -(Sy/2 - self.hp.sim.dpml)

        coeffs_taper_top = {}
        coeffs_taper_bot = {}
        powers_taper_top = {}
        angles_taper_top = {}
        powers_taper_bot = {}
        angles_taper_bot = {}
        powers_bend_top  = {}
        angles_bend_top  = {}
        powers_bend_bot  = {}
        angles_bend_bot  = {}
        powers_in_top_fwd = {}
        powers_in_bot_fwd = {}
        powers_in_top_bwd = {}
        powers_in_bot_bwd = {}

        # Build port configs filtered by configured input_ports
        if has_input_bends:
            _all_port_cfgs = {
                "top": (self.hp.geo.waveguide_offset/2, -1),
                "bot": (-self.hp.geo.waveguide_offset/2, 1),
            }
        else:
            _all_port_cfgs = {
                "top": (self.hp.sim.oxide_width/4, None),
                "bot": (-self.hp.sim.oxide_width/4, None),
            }
        input_port_configs = [(p,) + _all_port_cfgs[p] for p in self.hp.sim.input_ports if p in _all_port_cfgs]

        sim_idx = 0
        for port_label, y_center, y_kpoint_sign in input_port_configs:
            for mode_idx_for_port, fdtd_input_mode in enumerate(self.hp.sim.input_modes):
                if mp.am_master(): self.master_logger.info(f"FDTD_GDS_BEND port={port_label} mode={fdtd_input_mode}")

                if has_input_bends:
                    # Source/monitor in pre-bend vertical waveguide, offset inward from PML face
                    # (mirrors how non-bend FDTD offsets from PML in x).
                    src_x = x_in_bend_exit
                    y_pml_face = (Sy/2 - self.hp.sim.dpml) * (-1 if port_label == "bot" else 1)
                    src_y  = y_pml_face - np.sign(y_pml_face) * self.hp.sim.source_position
                    mon_y_in = y_pml_face - np.sign(y_pml_face) * self.hp.sim.sm_monitor_position
                    source_gaus = [mp.EigenModeSource(
                        mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                        eig_band=fdtd_input_mode,
                        direction=mp.NO_DIRECTION,
                        eig_kpoint=mp.Vector3(y=float(y_kpoint_sign)),
                        size=mp.Vector3(Sz - 2*self.hp.sim.dpml, 0, self.hp.sim.oxide_height),
                        center=[src_x, src_y, 0],
                    )]
                else:
                    source_gaus = [mp.EigenModeSource(
                        mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                        eig_band=fdtd_input_mode,
                        direction=mp.NO_DIRECTION,
                        eig_kpoint=mp.Vector3(x=1),
                        size=mp.Vector3(0, self.hp.sim.oxide_width/2, self.hp.sim.oxide_height),
                        center=[-Sx/2+self.hp.sim.dpml+self.hp.sim.source_position, y_center, 0],
                    )]

                sim = mp.Simulation(
                    cell_size=mp.Vector3(Sx, Sy, Sz),
                    boundary_layers=[mp.PML(self.hp.sim.dpml)],
                    geometry=geometry_base,
                    sources=source_gaus,
                    eps_averaging=True,
                    resolution=self.hp.sim.res,
                    default_material=SiO2
                )

                if has_input_bends:
                    # in_flux: pre-bend waveguide (y-flux); source normalization monitor
                    in_flux = sim.add_mode_monitor(freqs, mp.FluxRegion(
                        center=mp.Vector3(x_in_bend_exit, mon_y_in, 0),
                        size=mp.Vector3(Sz - 2*self.hp.sim.dpml, 0, self.hp.sim.oxide_height)
                    ))
                    # Post-bend input monitors (x-flux) — both stubs, measuring fwd + bwd power
                    taper_in_flux_top = sim.add_mode_monitor(freqs, mp.FluxRegion(
                        center=mp.Vector3(x_in_taper_mon, self.hp.geo.waveguide_offset/2, 0),
                        size=mp.Vector3(0, self.hp.geo.waveguide_offset, Sz - 2*self.hp.sim.dpml)
                    ))
                    taper_in_flux_bot = sim.add_mode_monitor(freqs, mp.FluxRegion(
                        center=mp.Vector3(x_in_taper_mon, -self.hp.geo.waveguide_offset/2, 0),
                        size=mp.Vector3(0, self.hp.geo.waveguide_offset, Sz - 2*self.hp.sim.dpml)
                    ))
                else:
                    in_flux = sim.add_mode_monitor(freqs, mp.FluxRegion(
                        center=mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position, y_center, 0),
                        size=mp.Vector3(0, self.hp.sim.oxide_width/2, Sz - 2*self.hp.sim.dpml)
                    ))

                # Taper output monitors (x-flux) — between stub and bend
                taper_flux_top = sim.add_mode_monitor(freqs, mp.FluxRegion(
                    center=mp.Vector3(x_taper_mon, self.hp.geo.waveguide_offset/2, 0),
                    size=mp.Vector3(0, self.hp.geo.waveguide_offset, Sz - 2*self.hp.sim.dpml)
                ))
                taper_flux_bot = sim.add_mode_monitor(freqs, mp.FluxRegion(
                    center=mp.Vector3(x_taper_mon, -self.hp.geo.waveguide_offset/2, 0),
                    size=mp.Vector3(0, self.hp.geo.waveguide_offset, Sz - 2*self.hp.sim.dpml)
                ))

                # Post-bend monitors (y-flux) — only if output bends present.
                # Monitor x-size = oxide_height so the right edge lands on the right PML inner face.
                _bend_fx_sx = self.hp.sim.oxide_height
                if has_output_bends:
                    bend_flux_top = sim.add_mode_monitor(freqs, mp.FluxRegion(
                        center=mp.Vector3(x_bend_exit, y_top_bend_mon, 0),
                        size=mp.Vector3(_bend_fx_sx, 0, Sz - 2*self.hp.sim.dpml)
                    ))
                    bend_flux_bot = sim.add_mode_monitor(freqs, mp.FluxRegion(
                        center=mp.Vector3(x_bend_exit, y_bot_bend_mon, 0),
                        size=mp.Vector3(_bend_fx_sx, 0, Sz - 2*self.hp.sim.dpml)
                    ))

                if mode_idx_for_port == 0:
                    if mp.am_master(): self.master_logger.info(f"Plotting FDTD_GDS_BEND geometry (port={port_label})")
                    if mp.am_master():
                        plt.figure()
                        sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, 0)))
                        plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_gds_bend_xy_{port_label}.png', dpi=600)
                        plt.figure()
                        sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, 0, Sz)))
                        plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_gds_bend_xz_{port_label}.png', dpi=600)
                        plt.figure()
                        sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(0, Sy, Sz)))
                        plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_gds_bend_yz_{port_label}.png', dpi=600)
                        plt.close("all")

                if mp.am_master():
                    monitors_dir = f'{self.hp.outputfolder}/layout/monitors'
                    xy_plane = mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, 0))

                    def _plot_component(cx, cy, sx_2d, sy_2d, label, color='cyan'):
                        fig, ax = plt.subplots(figsize=(10, 6))
                        sim.plot2D(ax=ax, output_plane=xy_plane,
                                   plot_sources_flag=False, plot_monitors_flag=False, plot_fields_flag=False)
                        if sx_2d == 0:
                            ax.plot([cx, cx], [cy - sy_2d/2, cy + sy_2d/2],
                                    color=color, linewidth=2.5, label=label)
                        elif sy_2d == 0:
                            ax.plot([cx - sx_2d/2, cx + sx_2d/2], [cy, cy],
                                    color=color, linewidth=2.5, label=label)
                        else:
                            from matplotlib.patches import Rectangle as MplRect
                            ax.add_patch(MplRect((cx - sx_2d/2, cy - sy_2d/2), sx_2d, sy_2d,
                                                 linewidth=2, edgecolor=color, facecolor='none', label=label))
                        ax.set_title(label.replace('_', ' '))
                        ax.legend(loc='upper right', fontsize=8)
                        fig.tight_layout()
                        fig.savefig(f'{monitors_dir}/{label}.png', dpi=150, bbox_inches='tight')
                        plt.close(fig)

                    # Source and port-specific input monitors
                    if has_input_bends:
                        _plot_component(src_x, src_y, Sz - 2*self.hp.sim.dpml, 0,
                                        f'{port_label}_input_source', color='red')
                        _plot_component(x_in_bend_exit, mon_y_in,
                                        Sz - 2*self.hp.sim.dpml, 0, f'{port_label}_input_pre_bend')
                        if sim_idx == 0:
                            _plot_component(x_in_taper_mon, self.hp.geo.waveguide_offset/2,
                                            0, self.hp.geo.waveguide_offset, 'top_input_post_bend')
                            _plot_component(x_in_taper_mon, -self.hp.geo.waveguide_offset/2,
                                            0, self.hp.geo.waveguide_offset, 'bot_input_post_bend')
                    else:
                        _plot_component(-Sx/2 + self.hp.sim.dpml + self.hp.sim.source_position,
                                        y_center, 0, self.hp.sim.oxide_width/2,
                                        f'{port_label}_input_source', color='red')
                        _plot_component(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position,
                                        y_center, 0, self.hp.sim.oxide_width/2,
                                        f'{port_label}_input_flux')

                    # Shared output monitors — only plot on first sim iteration
                    if sim_idx == 0:
                        _plot_component(x_taper_mon, self.hp.geo.waveguide_offset/2,
                                        0, self.hp.geo.waveguide_offset, 'top_output_pre_bend')
                        _plot_component(x_taper_mon, -self.hp.geo.waveguide_offset/2,
                                        0, self.hp.geo.waveguide_offset, 'bot_output_pre_bend')
                        if has_output_bends:
                            _plot_component(x_bend_exit, y_top_bend_mon,
                                            _bend_fx_sx, 0, 'top_output_post_bend')
                            _plot_component(x_bend_exit, y_bot_bend_mon,
                                            _bend_fx_sx, 0, 'bot_output_post_bend')

                if mp.am_master(): self.master_logger.info(f"Running FDTD_GDS_BEND (port={port_label} mode={fdtd_input_mode})")
                sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol=1e-11))

                if has_input_bends:
                    src_coeff = sim.get_eigenmode_coefficients(in_flux, [fdtd_input_mode], eig_parity=mp.NO_PARITY,
                                                               kpoint_func=lambda a, b: mp.Vector3(y=float(y_kpoint_sign)))
                    # For a y-flux monitor MEEP defines alpha[:,:,0] as +y ("forward") regardless of kpoint.
                    # When y_kpoint_sign < 0 the source propagates in -y ("backward"), so the injected
                    # power is in alpha[:,:,1]; using alpha[:,:,0] (the small reflection) as the
                    # denominator would inflate all S-parameters by ~100-1000x.
                    src_fwd_idx = 0 if y_kpoint_sign > 0 else 1
                else:
                    src_coeff = sim.get_eigenmode_coefficients(in_flux, [fdtd_input_mode], eig_parity=mp.NO_PARITY,
                                                               kpoint_func=lambda a, b: mp.Vector3(x=1))
                    src_fwd_idx = 0  # x-flux: +x is always "forward"

                if has_input_bends:
                    in_key = (port_label, fdtd_input_mode)
                    c_it = sim.get_eigenmode_coefficients(taper_in_flux_top, [fdtd_input_mode], eig_parity=mp.NO_PARITY,
                                                          kpoint_func=lambda a, b: mp.Vector3(x=1))
                    c_ib = sim.get_eigenmode_coefficients(taper_in_flux_bot, [fdtd_input_mode], eig_parity=mp.NO_PARITY,
                                                          kpoint_func=lambda a, b: mp.Vector3(x=1))
                    s_it_fwd = (c_it.alpha[:, :, 0] / src_coeff.alpha[:, :, src_fwd_idx])[0, :]
                    s_ib_fwd = (c_ib.alpha[:, :, 0] / src_coeff.alpha[:, :, src_fwd_idx])[0, :]
                    s_it_bwd = (c_it.alpha[:, :, 1] / src_coeff.alpha[:, :, src_fwd_idx])[0, :]
                    s_ib_bwd = (c_ib.alpha[:, :, 1] / src_coeff.alpha[:, :, src_fwd_idx])[0, :]
                    powers_in_top_fwd[in_key] = [np.abs(s)**2 for s in s_it_fwd]
                    powers_in_bot_fwd[in_key] = [np.abs(s)**2 for s in s_ib_fwd]
                    powers_in_top_bwd[in_key] = [np.abs(s)**2 for s in s_it_bwd]
                    powers_in_bot_bwd[in_key] = [np.abs(s)**2 for s in s_ib_bwd]
                    for pt in range(len(wavelengths)):
                        if mp.am_master(): self.master_logger.info(
                            f"InTaper | Port {port_label} In {fdtd_input_mode} | "
                            f"TopFwd: {powers_in_top_fwd[in_key][pt]:.4f}  TopBwd: {powers_in_top_bwd[in_key][pt]:.4f}  "
                            f"BotFwd: {powers_in_bot_fwd[in_key][pt]:.4f}  BotBwd: {powers_in_bot_bwd[in_key][pt]:.4f}")

                for out_mode in self.hp.sim.output_modes:
                    key = (port_label, fdtd_input_mode, out_mode)

                    ct = sim.get_eigenmode_coefficients(taper_flux_top, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
                    cb = sim.get_eigenmode_coefficients(taper_flux_bot, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
                    s_top = (ct.alpha[:, :, 0] / src_coeff.alpha[:, :, src_fwd_idx])[0, :]
                    s_bot = (cb.alpha[:, :, 0] / src_coeff.alpha[:, :, src_fwd_idx])[0, :]
                    coeffs_taper_top[key] = s_top
                    coeffs_taper_bot[key] = s_bot
                    pwr_t, ang_t, pwr_b, ang_b = [], [], [], []
                    for pt in range(len(wavelengths)):
                        pwr_t.append(np.abs(s_top[pt])**2); ang_t.append(np.angle(s_top[pt], deg=True))
                        pwr_b.append(np.abs(s_bot[pt])**2); ang_b.append(np.angle(s_bot[pt], deg=True))
                        if mp.am_master(): self.master_logger.info(f"Taper | Port {port_label} In {fdtd_input_mode} Out {out_mode} | Top: {pwr_t[-1]:.4f}, Bot: {pwr_b[-1]:.4f}")
                    powers_taper_top[key] = pwr_t; angles_taper_top[key] = ang_t
                    powers_taper_bot[key] = pwr_b; angles_taper_bot[key] = ang_b

                    if has_output_bends:
                        cbt = sim.get_eigenmode_coefficients(bend_flux_top, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(y=1))
                        cbb = sim.get_eigenmode_coefficients(bend_flux_bot, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(y=-1))
                        sb_top = (cbt.alpha[:, :, 0] / src_coeff.alpha[:, :, src_fwd_idx])[0, :]
                        # Bottom output bend exits in -y; for a y-flux monitor alpha[:,:,1] is the -y component.
                        sb_bot = (cbb.alpha[:, :, 1] / src_coeff.alpha[:, :, src_fwd_idx])[0, :]
                        pwr_bt, ang_bt, pwr_bb, ang_bb = [], [], [], []
                        for pt in range(len(wavelengths)):
                            pwr_bt.append(np.abs(sb_top[pt])**2); ang_bt.append(np.angle(sb_top[pt], deg=True))
                            pwr_bb.append(np.abs(sb_bot[pt])**2); ang_bb.append(np.angle(sb_bot[pt], deg=True))
                            if mp.am_master(): self.master_logger.info(f"Bend  | Port {port_label} In {fdtd_input_mode} Out {out_mode} | Top: {pwr_bt[-1]:.4f}, Bot: {pwr_bb[-1]:.4f}")
                        powers_bend_top[key] = pwr_bt; angles_bend_top[key] = ang_bt
                        powers_bend_bot[key] = pwr_bb; angles_bend_bot[key] = ang_bb

                sim_idx += 1

        data_to_save = {
            'coeffs_taper_top':   coeffs_taper_top,
            'coeffs_taper_bot':   coeffs_taper_bot,
            'powers_taper_top':   powers_taper_top,
            'angles_taper_top':   angles_taper_top,
            'powers_taper_bot':   powers_taper_bot,
            'angles_taper_bot':   angles_taper_bot,
            'powers_bend_top':    powers_bend_top,
            'angles_bend_top':    angles_bend_top,
            'powers_bend_bot':    powers_bend_bot,
            'angles_bend_bot':    angles_bend_bot,
            'powers_in_top_fwd':  powers_in_top_fwd,
            'powers_in_bot_fwd':  powers_in_bot_fwd,
            'powers_in_top_bwd':  powers_in_top_bwd,
            'powers_in_bot_bwd':  powers_in_bot_bwd,
            'wavelengths': wavelengths,
        }
        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)

        return powers_taper_top, powers_taper_bot, powers_bend_top, powers_bend_bot, wavelengths

    def eme_simulation(self):
        wavelengths = np.linspace(self.hp.sim.wavelength_min,self.hp.sim.wavelength_max,self.hp.sim.wavelength_pts)
        # Load the GDSII file for the MMI splitter
        mmi_top = load_2x2_MMI_top(self.hp.geo.mm_length,self.hp.sim.waveguide_length,self.hp.geo.waveguide_width,self.hp.geo.waveguide_offset,self.hp.sim.oxide_width,0.0,self.hp.geo.taper_length,self.hp.geo.mm_width,self.hp.geo.taper_width)
        plt.figure()
        mmi_top.plot()
        plt.savefig(f'{self.hp.outputfolder}/layout/eme_gds_top.png', dpi=600)

        extrusion_rules = stack_extrusions_GF(self.hp.geo.waveguide_height, self.hp.sim.oxide_height)
        structs = mw.extrude_gds(mmi_top, extrusion_rules)            
        # EME cell for the MMI region
        taper_Ls = np.array([self.hp.geo.taper_length / self.hp.sim.num_cells_taper for _ in range(self.hp.sim.num_cells_taper)])
        Ls = np.concatenate(([self.hp.sim.waveguide_length], taper_Ls, [self.hp.geo.mm_length], taper_Ls, [self.hp.sim.waveguide_length]))
        num_output_modes_eme = max(self.hp.sim.output_modes)
        num_input_modes_eme = max(self.hp.sim.input_modes)
        num_modes = [num_input_modes_eme] + [self.hp.sim.num_modes_taper for _ in range(self.hp.sim.num_cells_taper)] + [self.hp.sim.num_modes_MMI] + [self.hp.sim.num_modes_taper for _ in range(self.hp.sim.num_cells_taper)] + [num_output_modes_eme]
        cells = mw.create_cells(
            structures=structs,
            mesh=mw.Mesh2D(
                x=np.linspace(-self.hp.sim.oxide_width/2, self.hp.sim.oxide_width/2, int(self.hp.sim.res*self.hp.sim.oxide_width) + 1),
                y=np.linspace(-self.hp.sim.oxide_height/2, self.hp.sim.oxide_height/2, int(self.hp.sim.res*self.hp.sim.oxide_height) + 1),
            ),
            Ls=Ls,
        )

        # Create seperate cells for output to bottom waveguide
        mmi_bot = load_2x2_MMI_bot(self.hp.geo.mm_length,self.hp.sim.waveguide_length,self.hp.geo.waveguide_width,self.hp.geo.waveguide_offset,self.hp.sim.oxide_width,0.0,self.hp.geo.taper_length,self.hp.geo.mm_width,self.hp.geo.taper_width)
        plt.figure()
        mmi_bot.plot()
        plt.savefig(f'{self.hp.outputfolder}/layout/eme_gds_bot.png', dpi=600)
        structs_bot = mw.extrude_gds(mmi_bot, extrusion_rules) 
        cells_bot = mw.create_cells(
            structures=structs_bot,
            mesh=mw.Mesh2D(
                x=np.linspace(-self.hp.sim.oxide_width/2, self.hp.sim.oxide_width/2, int(self.hp.sim.res*self.hp.sim.oxide_width) + 1),
                y=np.linspace(-self.hp.sim.oxide_height/2, self.hp.sim.oxide_height/2, int(self.hp.sim.res*self.hp.sim.oxide_height) + 1),
            ),
            Ls=Ls,
        )

        input_modes = self.hp.sim.input_modes
        n_in_max = max(input_modes)
        keys = [(in_m, out_m) for in_m in input_modes for out_m in self.hp.sim.output_modes]
        coeffs_top_by_mode = {k: [] for k in keys}
        powers_top_by_mode = {k: [] for k in keys}
        angles_top_by_mode = {k: [] for k in keys}
        coeffs_bot_by_mode = {k: [] for k in keys}
        powers_bot_by_mode = {k: [] for k in keys}
        angles_bot_by_mode = {k: [] for k in keys}

        for wavelength in wavelengths:
            if mp.am_master(): self.master_logger.info(f"Starting EME simulation at wavelength {wavelength}")

            gf_environment = mw.Environment(wl=wavelength, T=25.0)
            css = [mw.CrossSection.from_cell(cell=cell, env=gf_environment) for cell in cells]
            cs_bot = mw.CrossSection.from_cell(cell=cells_bot[-1], env=gf_environment)

            for cs_num, cs in enumerate(css):
                plt.clf()
                plt.figure()
                mw.visualize(cs, show=False)
                plt.savefig(f'{self.hp.outputfolder}/layout/cross_sections/cs_{cs_num}.png', dpi=600)
                plt.close()

            # Perform mode simulations
            if mp.am_master(): self.master_logger.info(f"Performing mode simulations")
            modes = []
            for cell_num, (cs, num_mode) in enumerate(zip(css, num_modes)):
                if mp.am_master(): self.master_logger.info(f"Computing {num_mode} modes of cell {cell_num}")
                modes.append(compute_modes(cs, num_modes=num_mode))
            if mp.am_master(): self.master_logger.info(f"Computing modes for bottom output")
            bot_mode = compute_modes(cs_bot, num_modes=num_output_modes_eme)

            if mp.am_master(): self.master_logger.info(f"Computing s matrices for output to top waveguide")
            S_top, port_map = mw.compute_s_matrix(modes, cells=cells)
            if mp.am_master(): self.master_logger.info(f"Computing s matrices for output to bottom waveguide")
            bot_modes = modes[0:-1]
            bot_modes.append(bot_mode)
            S_bot, port_map = mw.compute_s_matrix(bot_modes, cells=cells_bot)

            for in_m in input_modes:
                for out_mode in self.hp.sim.output_modes:
                    # S-matrix rows: input port modes 0..n_in_max-1; cols: output port modes n_in_max..n_in_max+n_out-1
                    s_top_val = S_top[in_m - 1, n_in_max + out_mode - 1]
                    s_bot_val = S_bot[in_m - 1, n_in_max + out_mode - 1]
                    coeffs_top_by_mode[(in_m, out_mode)].append(s_top_val)
                    powers_top_by_mode[(in_m, out_mode)].append(np.abs(s_top_val)**2)
                    angles_top_by_mode[(in_m, out_mode)].append(np.angle(s_top_val, deg=True))
                    coeffs_bot_by_mode[(in_m, out_mode)].append(s_bot_val)
                    powers_bot_by_mode[(in_m, out_mode)].append(np.abs(s_bot_val)**2)
                    angles_bot_by_mode[(in_m, out_mode)].append(np.angle(s_bot_val, deg=True))
                    if mp.am_master(): self.master_logger.info(f"In {in_m} Out {out_mode} | Top: Power {powers_top_by_mode[(in_m, out_mode)][-1]:.4f}, Bot: Power {powers_bot_by_mode[(in_m, out_mode)][-1]:.4f}")

            if wavelength == wavelengths[int(np.floor(len(wavelengths)/2))]:
                for in_m in input_modes:
                    plot_eme_fields(self.hp, cells, modes, tag=in_m, input_mode=in_m)
                if mp.am_master(): self.master_logger.info(f"Transmission at center wavelength: top={list(powers_top_by_mode.values())}, bot={list(powers_bot_by_mode.values())}")

        data_to_save = {
            'coeffs_top_by_mode': coeffs_top_by_mode,
            'powers_top_by_mode': powers_top_by_mode,
            'angles_top_by_mode': angles_top_by_mode,
            'coeffs_bot_by_mode': coeffs_bot_by_mode,
            'powers_bot_by_mode': powers_bot_by_mode,
            'angles_bot_by_mode': angles_bot_by_mode,
            'wavelengths': wavelengths,
        }

        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)
        return powers_top_by_mode, powers_bot_by_mode, wavelengths

    def run(self):
        if self.hp.sim.sim_type == "HYBRID":
            powers_top_by_mode, powers_bot_by_mode, wavelengths = self.hybrid_simulation()
        elif self.hp.sim.sim_type == "FDTD":
            powers_top_by_mode, powers_bot_by_mode, wavelengths = self.fdtd_simulation()
        elif self.hp.sim.sim_type == "FDTD_GDS":
            powers_top_by_mode, powers_bot_by_mode, wavelengths = self.fdtd_simulation_from_gds()
        elif self.hp.sim.sim_type == "FDTD_GDS_BEND":
            pt_top, pt_bot, pb_top, pb_bot, wavelengths = self.fdtd_simulation_from_gds_with_bends()
            if mp.am_master():
                def _label(key):
                    port, in_m, out_m = key
                    return f"Port {port} In {in_m} → Out {out_m}"
                all_ports = sorted(set(k[0] for k in pt_top.keys()))
                for port in all_ports:
                    port_suffix = f"_{port}" if len(all_ports) > 1 else ""
                    plt.figure()
                    for key, powers in pt_top.items():
                        if key[0] == port: plt.plot(wavelengths, powers, label=f"Taper top {_label(key)}")
                    for key, powers in pt_bot.items():
                        if key[0] == port: plt.plot(wavelengths, powers, label=f"Taper bot {_label(key)}")
                    for key, powers in pb_top.items():
                        if key[0] == port: plt.plot(wavelengths, powers, label=f"Bend top {_label(key)}", linestyle="--")
                    for key, powers in pb_bot.items():
                        if key[0] == port: plt.plot(wavelengths, powers, label=f"Bend bot {_label(key)}", linestyle="--")
                    plt.xlabel("Wavelength (micron)"); plt.ylabel("Transmission"); plt.legend()
                    plt.savefig(f'{self.hp.outputfolder}/performance{port_suffix}.png', dpi=600)
                    plt.figure()
                    for key, powers in pt_top.items():
                        if key[0] == port: plt.plot(wavelengths, 10*np.log10(np.array(powers)), label=f"Taper top {_label(key)}")
                    for key, powers in pt_bot.items():
                        if key[0] == port: plt.plot(wavelengths, 10*np.log10(np.array(powers)), label=f"Taper bot {_label(key)}")
                    for key, powers in pb_top.items():
                        if key[0] == port: plt.plot(wavelengths, 10*np.log10(np.array(powers)), label=f"Bend top {_label(key)}", linestyle="--")
                    for key, powers in pb_bot.items():
                        if key[0] == port: plt.plot(wavelengths, 10*np.log10(np.array(powers)), label=f"Bend bot {_label(key)}", linestyle="--")
                    plt.xlabel("Wavelength (micron)"); plt.ylabel("Transmission (dB)"); plt.legend()
                    plt.savefig(f'{self.hp.outputfolder}/performance_dB{port_suffix}.png', dpi=600)
                plt.close("all")
            return
        elif self.hp.sim.sim_type == "EME":
            powers_top_by_mode, powers_bot_by_mode, wavelengths = self.eme_simulation()

        # EME/HYBRID keys are (in_m, out_m); FDTD/FDTD_GDS keys are (port, in_m, out_m)
        _sample_key = next(iter(powers_top_by_mode))
        _has_port = len(_sample_key) == 3
        if _has_port:
            all_ports = sorted(set(k[0] for k in powers_top_by_mode.keys()))
        else:
            all_ports = [None]
        for port in all_ports:
            port_suffix = f"_{port}" if port is not None and len(all_ports) > 1 else ""
            plt.figure()
            for key, powers in powers_top_by_mode.items():
                if port is None or key[0] == port:
                    in_m, out_m = (key[0], key[1]) if not _has_port else (key[1], key[2])
                    plt.plot(wavelengths, powers, label=f"In {in_m} → Top {out_m}")
            for key, powers in powers_bot_by_mode.items():
                if port is None or key[0] == port:
                    in_m, out_m = (key[0], key[1]) if not _has_port else (key[1], key[2])
                    plt.plot(wavelengths, powers, label=f"In {in_m} → Bot {out_m}")
            plt.xlabel("Wavelength (micron)")
            plt.ylabel("Transmission")
            plt.legend()
            plt.savefig(f'{self.hp.outputfolder}/performance{port_suffix}.png', dpi=600)

            plt.figure()
            for key, powers in powers_top_by_mode.items():
                if port is None or key[0] == port:
                    in_m, out_m = (key[0], key[1]) if not _has_port else (key[1], key[2])
                    plt.plot(wavelengths, 10 * np.log10(np.array(powers)), label=f"In {in_m} → Top {out_m}")
            for key, powers in powers_bot_by_mode.items():
                if port is None or key[0] == port:
                    in_m, out_m = (key[0], key[1]) if not _has_port else (key[1], key[2])
                    plt.plot(wavelengths, 10 * np.log10(np.array(powers)), label=f"In {in_m} → Bot {out_m}")
            plt.xlabel("Wavelength (micron)")
            plt.ylabel("Transmission (dB)")
            plt.legend()
            plt.savefig(f'{self.hp.outputfolder}/performance_dB{port_suffix}.png', dpi=600)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file", required=True)
    args = parser.parse_args()
    mmi = MMI_2x2(args)

    mmi.run()