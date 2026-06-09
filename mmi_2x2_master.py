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

from meow import compute_modes_meep as compute_modes

from engine.hyperparams.MMI import HyperParamsMMI as HyperParams

from engine.helpers.general import(
    make_output_folders
)
from engine.helpers.meep import(
     meep_material_from_meow,
     create_temp_gds_for_mmi_2x2_input,
     create_temp_gds_for_mmi_2x2_output,
     create_temp_gds_for_full_mmi_2x2,
     prepare_gds_for_meep,
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
        make_output_folders(self.hp.outputfolder)

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
            plt.savefig(f'{self.hp.outputfolder}/layout/top_xy.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/top_xz.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/top_yz.png', dpi=300)

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
            plt.savefig(f'{self.hp.outputfolder}/layout/bot_xy.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/bot_xz.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/bot_yz.png', dpi=300)

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

        # Set waveguide geometry
        geometry = [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2,self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry += [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2,-self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2,self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2,-self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]

        # Create temp gds file for MMI geometry (master only, then sync)
        if mp.am_master():
            create_temp_gds_for_full_mmi_2x2(self.hp.geo.mm_length, self.hp.geo.taper_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset)
        mp.all_wait()

        # Get MMI geometry from temp gds file
        mmi = mp.get_GDSII_prisms(Si, "gds/temp_full.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry += mmi

        # FDTD supports a single eigenmode source; use the first requested input mode
        fdtd_input_mode = self.hp.sim.input_modes[0]

        # Create source MEEP simulation
        source_gaus = [
            mp.EigenModeSource(
                mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                eig_band=fdtd_input_mode,
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

        # Plot simulation geometry
        if mp.am_master(): self.master_logger.info("Plotting FDTD simulation geometry")
        plt.figure()
        sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, Sy, 0)))
        plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_xy.png', dpi=300)
        plt.figure()
        sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
        plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_xz.png', dpi=300)
        plt.figure()
        sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
        plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_yz.png', dpi=300)

        # run simulation
        if mp.am_master(): self.master_logger.info("Running FDTD simulation")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

        # get eigenmode coefficients
        if mp.am_master(): self.master_logger.info("Getting FDTD simulation eigenmode coefficients:")
        src_coeff = sim.get_eigenmode_coefficients(in_flux, [fdtd_input_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))

        coeffs_top_by_mode = {}
        coeffs_bot_by_mode = {}
        powers_top_by_mode = {}
        angles_top_by_mode = {}
        powers_bot_by_mode = {}
        angles_bot_by_mode = {}
        for out_mode in self.hp.sim.output_modes:
            coeff_top = sim.get_eigenmode_coefficients(out_flux_top, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
            coeff_bot = sim.get_eigenmode_coefficients(out_flux_bot, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
            s_top = coeff_top.alpha[:,:,0] / src_coeff.alpha[:,:,0]
            s_top = s_top[0,:]
            s_bot = coeff_bot.alpha[:,:,0] / src_coeff.alpha[:,:,0]
            s_bot = s_bot[0,:]
            key = (fdtd_input_mode, out_mode)
            coeffs_top_by_mode[key] = s_top
            coeffs_bot_by_mode[key] = s_bot
            pwr_top, ang_top, pwr_bot, ang_bot = [], [], [], []
            for point in range(len(wavelengths)):
                pwr_top.append(np.abs(s_top[point])**2)
                ang_top.append(np.angle(s_top[point], deg=True))
                pwr_bot.append(np.abs(s_bot[point])**2)
                ang_bot.append(np.angle(s_bot[point], deg=True))
                if mp.am_master(): self.master_logger.info(f"In {fdtd_input_mode} Out {out_mode} | Top: Power {pwr_top[-1]:.4f}, Angle {ang_top[-1]:.2f}")
                if mp.am_master(): self.master_logger.info(f"In {fdtd_input_mode} Out {out_mode} | Bot: Power {pwr_bot[-1]:.4f}, Angle {ang_bot[-1]:.2f}")
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
        device_Sx = prepare_gds_for_meep(self.hp.geo.gds_file, "gds/temp_gds_fdtd.gds")
        Sx = device_Sx + 2*self.hp.sim.waveguide_length + 2*self.hp.sim.dpml
        Sy = self.hp.sim.oxide_width + 2*self.hp.sim.dpml
        Sz = self.hp.sim.oxide_height + 2*self.hp.sim.dpml

        # materials
        Si = meep_material_from_meow(mw.silicon.n, mw.silicon.params['wl'], wavelength_center)
        SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'], wavelength_center)

        # Waveguide stubs bridging PML to device input/output faces
        geometry = [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2, self.hp.geo.waveguide_offset/2, 0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material=Si)]
        geometry += [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2, -self.hp.geo.waveguide_offset/2, 0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material=Si)]
        geometry += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2, self.hp.geo.waveguide_offset/2, 0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material=Si)]
        geometry += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2, -self.hp.geo.waveguide_offset/2, 0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material=Si)]

        # Get MMI geometry from the centered GDS file
        mmi = mp.get_GDSII_prisms(Si, "gds/temp_gds_fdtd.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry += mmi

        # FDTD supports a single eigenmode source; use the first requested input mode
        fdtd_input_mode = self.hp.sim.input_modes[0]

        # Create source MEEP simulation
        source_gaus = [
            mp.EigenModeSource(
                mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                eig_band=fdtd_input_mode,
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
                center=mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position, self.hp.sim.oxide_width/4, 0),
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

        # Plot simulation geometry
        if mp.am_master(): self.master_logger.info("Plotting FDTD GDS simulation geometry")
        plt.figure()
        sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, 0)))
        plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_gds_xy.png', dpi=300)
        plt.figure()
        sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, 0, Sz)))
        plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_gds_xz.png', dpi=300)
        plt.figure()
        sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(0, Sy, Sz)))
        plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_gds_yz.png', dpi=300)

        # Run simulation
        if mp.am_master(): self.master_logger.info("Running FDTD GDS simulation")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol=1e-11))

        # Get eigenmode coefficients
        if mp.am_master(): self.master_logger.info("Getting FDTD GDS simulation eigenmode coefficients:")
        src_coeff = sim.get_eigenmode_coefficients(in_flux, [fdtd_input_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))

        coeffs_top_by_mode = {}
        coeffs_bot_by_mode = {}
        powers_top_by_mode = {}
        angles_top_by_mode = {}
        powers_bot_by_mode = {}
        angles_bot_by_mode = {}
        for out_mode in self.hp.sim.output_modes:
            coeff_top = sim.get_eigenmode_coefficients(out_flux_top, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
            coeff_bot = sim.get_eigenmode_coefficients(out_flux_bot, [out_mode], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
            s_top = coeff_top.alpha[:, :, 0] / src_coeff.alpha[:, :, 0]
            s_top = s_top[0, :]
            s_bot = coeff_bot.alpha[:, :, 0] / src_coeff.alpha[:, :, 0]
            s_bot = s_bot[0, :]
            key = (fdtd_input_mode, out_mode)
            coeffs_top_by_mode[key] = s_top
            coeffs_bot_by_mode[key] = s_bot
            pwr_top, ang_top, pwr_bot, ang_bot = [], [], [], []
            for point in range(len(wavelengths)):
                pwr_top.append(np.abs(s_top[point])**2)
                ang_top.append(np.angle(s_top[point], deg=True))
                pwr_bot.append(np.abs(s_bot[point])**2)
                ang_bot.append(np.angle(s_bot[point], deg=True))
                if mp.am_master(): self.master_logger.info(f"In {fdtd_input_mode} Out {out_mode} | Top: Power {pwr_top[-1]:.4f}, Angle {ang_top[-1]:.2f}")
                if mp.am_master(): self.master_logger.info(f"In {fdtd_input_mode} Out {out_mode} | Bot: Power {pwr_bot[-1]:.4f}, Angle {ang_bot[-1]:.2f}")
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

    def eme_simulation(self):
        wavelengths = np.linspace(self.hp.sim.wavelength_min,self.hp.sim.wavelength_max,self.hp.sim.wavelength_pts)
        # Load the GDSII file for the MMI splitter
        mmi_top = load_2x2_MMI_top(self.hp.geo.mm_length,self.hp.sim.waveguide_length,self.hp.geo.waveguide_width,self.hp.geo.waveguide_offset,self.hp.sim.oxide_width,0.0,self.hp.geo.taper_length,self.hp.geo.mm_width,self.hp.geo.taper_width)
        plt.figure()
        mmi_top.plot()
        plt.savefig(f'{self.hp.outputfolder}/layout/eme_gds_top.png', dpi=300)

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
        plt.savefig(f'{self.hp.outputfolder}/layout/eme_gds_bot.png', dpi=300)
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
                plt.savefig(f'{self.hp.outputfolder}/layout/cross_sections/cs_{cs_num}.png', dpi=300)

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
        elif self.hp.sim.sim_type == "EME":
            powers_top_by_mode, powers_bot_by_mode, wavelengths = self.eme_simulation()

        plt.figure()
        for (in_m, out_m), powers in powers_top_by_mode.items():
            plt.plot(wavelengths, powers, label=f"In {in_m} → Top {out_m}")
        for (in_m, out_m), powers in powers_bot_by_mode.items():
            plt.plot(wavelengths, powers, label=f"In {in_m} → Bot {out_m}")
        plt.xlabel("Wavelength (micron)")
        plt.ylabel("Transmission")
        plt.legend()
        plt.savefig(f'{self.hp.outputfolder}/performance.png', dpi=300)

        plt.figure()
        for (in_m, out_m), powers in powers_top_by_mode.items():
            plt.plot(wavelengths, 10 * np.log10(np.array(powers)), label=f"In {in_m} → Top {out_m}")
        for (in_m, out_m), powers in powers_bot_by_mode.items():
            plt.plot(wavelengths, 10 * np.log10(np.array(powers)), label=f"In {in_m} → Bot {out_m}")
        plt.xlabel("Wavelength (micron)")
        plt.ylabel("Transmission (dB)")
        plt.legend()
        plt.savefig(f'{self.hp.outputfolder}/performance_dB.png', dpi=300)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file", required=True)
    args = parser.parse_args()
    mmi = MMI_2x2(args)

    mmi.run()