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

from engine.hyperparams.MMI_1x2 import HyperParamsMMI as HyperParams

from engine.helpers.general import(
    make_output_folders
)
from engine.helpers.meep import(
     meep_material_from_meow, 
     create_temp_gds_for_mmi_2x2_input, 
     create_temp_gds_for_mmi_2x2_output,
     create_temp_gds_for_full_mmi_2x2
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

    def top_simulation(self,wavelengths):
        self.master_logger.info(f"Creating top FDTD simulation")
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

        # Create temp gds file for MMI geometry
        create_temp_gds_for_mmi_2x2_input(self.hp.geo.mm_length, self.hp.geo.taper_length, Sx, self.hp.sim.dpml, self.hp.sim.waveguide_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset,"top")

        # Get MMI geometry from temp gds file
        mmi = mp.get_GDSII_prisms(Si, "gds/temp.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry += mmi

        # Create source MEEP simulation
        source_gaus = [
            mp.EigenModeSource(
                mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                eig_band=1,
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

        # Plot simulation geometry
        self.master_logger.info("Plotting top FDTD simulation geometry")
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
        self.master_logger.info("Running top FDTD simulation")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

        # get eigenmode coefficients
        self.master_logger.info("Getting top FDTD simulation eigenmode coefficients:")
        src_coeff = sim.get_eigenmode_coefficients(in_flux,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
        top_coeffs = []
        for mode_num in range(1,self.hp.sim.num_modes_MMI+1):
            coeff = sim.get_eigenmode_coefficients(out_flux,[mode_num],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
            s = coeff.alpha[:,:,0]/src_coeff.alpha[:,:,0]
            top_coeffs.append(s[0,:])
            self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(s[0,:])**2}, Angle: {np.angle(s[0,:], deg=True)}")
        top_coeffs = np.asarray(top_coeffs)
        return top_coeffs.T.astype(list)

    def bot_simulation(self,wavelengths):
        self.master_logger.info(f"Creating bottom FDTD simulation")
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

        # Create temp gds file for MMI geometry
        create_temp_gds_for_mmi_2x2_input(self.hp.geo.mm_length, self.hp.geo.taper_length, Sx, self.hp.sim.dpml, self.hp.sim.waveguide_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset,"bot")

        # Get MMI geometry from temp gds file
        mmi = mp.get_GDSII_prisms(Si, "gds/temp.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry += mmi

        # Create source MEEP simulation
        source_gaus = [
            mp.EigenModeSource(
                mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                eig_band=1,
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

        # Plot simulation geometry
        self.master_logger.info("Plotting bottom FDTD simulation geometry")
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
        self.master_logger.info("Running bottom FDTD simulation")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

        # get eigenmode coefficients
        self.master_logger.info("Getting bottom FDTD simulation eigenmode coefficients:")
        src_coeff = sim.get_eigenmode_coefficients(in_flux,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
        top_coeffs = []
        for mode_num in range(1,self.hp.sim.num_modes_MMI+1):
            coeff = sim.get_eigenmode_coefficients(out_flux,[mode_num],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
            s = coeff.alpha[:,:,0]/src_coeff.alpha[:,:,0]
            top_coeffs.append(s[0,:])
            self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(s[0,:])**2}, Angle: {np.angle(s[0,:], deg=True)}")
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
            self.master_logger.info(f"Starting EME simulation for multimode section at wavelength {wavelength}")
            
            gf_environment = mw.Environment(wl=wavelength, T=25.0)
            css = [mw.CrossSection.from_cell(cell=cell, env=gf_environment) for cell in cells]

            # Performe mode simulations
            self.master_logger.info(f"Performing mode simulations")
            modes = []
            for cell_num, (cs, num_mode) in enumerate(zip(css,num_modes)):
                self.master_logger.info(f"Computing mode {num_mode} of cell {cell_num}")
                modes.append(mw.compute_modes_meep(cs, num_modes=num_mode))

            S, port_map = mw.compute_s_matrix(modes, cells=cells)
            mm_coeffs_1_wl = []
            for mode_num in range(self.hp.sim.num_modes_MMI):
                mm_coeffs_1_wl.append(S[mode_num,self.hp.sim.num_modes_MMI+mode_num])
                self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(S[mode_num,self.hp.sim.num_modes_MMI+mode_num])**2}, Angle: {np.angle(S[mode_num,self.hp.sim.num_modes_MMI+mode_num], deg=True)}")
            mm_coeffs.append(mm_coeffs_1_wl)
        return mm_coeffs

    def hybrid_simulation(self):
        wavelengths = np.linspace(self.hp.sim.wavelength_min,self.hp.sim.wavelength_max,self.hp.sim.wavelength_pts)
        top_coeffs = self.top_simulation(wavelengths)
        bot_coeffs = self.bot_simulation(wavelengths)
        mm_coeffs = self.mm_simulation(wavelengths)
        
        powers_top = []
        angles_top = []
        powers_bot = []
        angles_bot = []
        for point in range(len(wavelengths)):
            top_coeff = top_coeffs[point]
            mm_coeff = mm_coeffs[point]
            bot_coeff = bot_coeffs[point]
            propagated_coeffs = [x * y for x, y in zip(top_coeff, mm_coeff)]
            propagated_coeffs_top = [x * y for x, y in zip(propagated_coeffs, top_coeff)]
            propagated_coeffs_bot = [x * y for x, y in zip(propagated_coeffs, bot_coeff)]

            powers_top.append(np.abs(sum(propagated_coeffs_top))**2)
            angles_top.append(np.angle(sum(propagated_coeffs_top), deg=True))
            powers_bot.append(np.abs(sum(propagated_coeffs_bot))**2)
            angles_bot.append(np.angle(sum(propagated_coeffs_bot), deg=True))

            self.master_logger.info(f"Top: Power: {powers_top[point]}, Angle: {angles_top[point]}")
            self.master_logger.info(f"Bottom: Power: {powers_bot[point]}, Angle: {angles_bot[point]}")

        data_to_save = {
            'top_coeffs': top_coeffs,
            'mm_coeffs': mm_coeffs,
            'bot_coeffs': bot_coeffs,
            'powers_top': powers_top,
            'angles_top': angles_top,
            'powers_bot': powers_bot,
            'angles_bot': angles_bot,
            'wavelengths': wavelengths,
        }

        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)
        
        return powers_top, powers_bot, wavelengths

    def fdtd_simulation(self):
        self.master_logger.info(f"Creating full FDTD simulation")
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

        # Create temp gds file for MMI geometry
        create_temp_gds_for_full_mmi_2x2(self.hp.geo.mm_length, self.hp.geo.taper_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset)

        # Get MMI geometry from temp gds file
        mmi = mp.get_GDSII_prisms(Si, "gds/temp_full.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry += mmi

        # Create source MEEP simulation
        source_gaus = [
            mp.EigenModeSource(
                mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                eig_band=1,
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
        self.master_logger.info("Plotting FDTD simulation geometry")
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
        self.master_logger.info("Running FDTD simulation")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

        # get eigenmode coefficients
        self.master_logger.info("Getting input FDTD simulation eigenmode coefficients:")
        src_coeff = sim.get_eigenmode_coefficients(in_flux,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
        coeff_top = sim.get_eigenmode_coefficients(out_flux_top,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
        coeff_bot = sim.get_eigenmode_coefficients(out_flux_bot,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
        
        s_top = coeff_top.alpha[:,:,0]/src_coeff.alpha[:,:,0]
        s_top = s_top[0,:]
        s_bot = coeff_bot.alpha[:,:,0]/src_coeff.alpha[:,:,0]
        s_bot = s_bot[0,:]

        powers_top = []
        angles_top = []
        powers_bot = []
        angles_bot = []
        for point in range(len(wavelengths)):
            powers_top.append(np.abs(s_top[point])**2)
            angles_top.append(np.angle(s_top[point], deg=True))
            powers_bot.append(np.abs(s_bot[point])**2)
            angles_bot.append(np.angle(s_bot[point], deg=True))
            self.master_logger.info(f"Top: Power: {powers_top[point]}, Angle: {angles_top[point]}")
            self.master_logger.info(f"Bottom: Power: {powers_bot[point]}, Angle: {angles_bot[point]}")

        data_to_save = {
            'coeffs_top': s_top,
            'coeffs_bot': s_bot,
            'powers_top': powers_top,
            'powers_bot': powers_bot,
            'angles_top': angles_top,
            'angles_bot': angles_bot,
            'wavelengths': wavelengths,
        }

        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)

        return powers_top, powers_bot, wavelengths
    
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
        num_modes = [1] + [self.hp.sim.num_modes_taper for _ in range(self.hp.sim.num_cells_taper)] + [self.hp.sim.num_modes_MMI] + [self.hp.sim.num_modes_taper for _ in range(self.hp.sim.num_cells_taper)] + [1]
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

        powers_top = []
        coeffs_top = []
        angles_top = []
        powers_bot = []
        coeffs_bot = []
        angles_bot = []
        for wavelength in wavelengths:
            self.master_logger.info(f"Starting EME simulation at wavelength {wavelength}")
            
            gf_environment = mw.Environment(wl=wavelength, T=25.0)
            css = [mw.CrossSection.from_cell(cell=cell, env=gf_environment) for cell in cells]
            cs_bot = mw.CrossSection.from_cell(cell=cells_bot[-1], env=gf_environment)

            for cs_num, cs in enumerate(css):
                plt.clf()
                plt.figure()
                mw.visualize(cs, show=False)
                plt.savefig(f'{self.hp.outputfolder}/layout/cross_sections/cs_{cs_num}.png', dpi=300)


            # Performe mode simulations
            self.master_logger.info(f"Performing mode simulations")
            modes = []
            for cell_num, (cs, num_mode) in enumerate(zip(css,num_modes)):
                self.master_logger.info(f"Computing {num_mode} modes of cell {cell_num}")
                modes.append(mw.compute_modes_meep(cs, num_modes=num_mode))
            self.master_logger.info(f"Computing mode for bottom output")
            bot_mode = mw.compute_modes_meep(cs_bot, num_modes=1)

            self.master_logger.info(f"Computing s matrices for output to top waveguide")
            S_top, port_map = mw.compute_s_matrix(modes, cells=cells)
            self.master_logger.info(f"Computing s matrices for output to bottom waveguide")
            bot_modes = modes[0:-1]
            bot_modes.append(bot_mode)
            S_bot, port_map = mw.compute_s_matrix(bot_modes, cells=cells_bot)

            coeffs_top.append(S_top[0,1])
            powers_top.append(np.abs(S_top[0,1])**2)
            angles_top.append(np.angle(S_top[0,1],deg=True))
            coeffs_bot.append(S_bot[0,1])
            powers_bot.append(np.abs(S_bot[0,1])**2)
            angles_bot.append(np.angle(S_bot[0,1],deg=True))

            if wavelength == wavelengths[int(np.floor(len(wavelengths)/2))]:
                plot_eme_fields(self.hp,cells,modes)
                self.master_logger.info(f"Transmission at center wavelength: {powers_top[-1]}, {powers_bot[-1]}")


        data_to_save = {
            'coeffs_top': coeffs_top,
            'powers_top': powers_top,
            'angles_top': angles_top,
            'coeffs_bot': coeffs_bot,
            'powers_bot': powers_bot,
            'angles_bot': angles_bot,
            'wavelengths': wavelengths,
        }

        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)
        return powers_top, powers_bot, wavelengths

    def run(self):
        if self.hp.sim.sim_type == "HYBRID":
            powers_top, powers_bot, wavelengths = self.hybrid_simulation()
        elif self.hp.sim.sim_type == "FDTD":
            powers_top, powers_bot, wavelengths = self.fdtd_simulation()
        elif self.hp.sim.sim_type == "EME":
            powers_top, powers_bot, wavelengths = self.eme_simulation()

        plt.figure()
        plt.plot(wavelengths,powers_top)
        plt.plot(wavelengths,powers_bot)
        plt.xlabel("Wavelength (micron)")
        plt.ylabel("Transmission")
        plt.savefig(f'{self.hp.outputfolder}/performance.png', dpi=300)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file", required=True)
    args = parser.parse_args()
    mmi = MMI_2x2(args)

    mmi.run()