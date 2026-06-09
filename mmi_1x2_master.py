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

from engine.hyperparams.MMI import HyperParamsMMI as HyperParams

from engine.helpers.general import(
    make_output_folders
)
from engine.helpers.meep import(
     meep_material_from_meow, 
     create_temp_gds_for_mmi_input, 
     create_temp_gds_for_mmi_output,
     create_temp_gds_for_full_mmi
)
from engine.helpers.meow import(
    load_1x2_MMI,
    stack_extrusions_GF
)
from engine.plotting.meow import(
    plot_eme_fields
)

class MMI_1x2:
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

    def input_simulation(self,wavelengths):
        if mp.am_master(): self.master_logger.info(f"Creating input FDTD simulation")
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
        geometry = [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2,0,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]

        # Create temp gds file for MMI geometry
        create_temp_gds_for_mmi_input(self.hp.geo.mm_length, self.hp.geo.taper_length, Sx, self.hp.sim.dpml, self.hp.sim.waveguide_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset)

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
                size=mp.Vector3(0, self.hp.sim.oxide_width, self.hp.sim.oxide_height),
                center=[-Sx/2+self.hp.sim.dpml+self.hp.sim.source_position, 0, 0],
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
                center = mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position, 0, 0), 
                size = mp.Vector3(0, Sy - 2*self.hp.sim.dpml,Sz - 2*self.hp.sim.dpml)
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
        if mp.am_master(): self.master_logger.info("Plotting input FDTD simulation geometry")
        if mp.am_master():
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, Sy, 0)))
            plt.savefig(f'{self.hp.outputfolder}/layout/input_xy.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/input_xz.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/input_yz.png', dpi=300)
            plt.close("all")

        # run simulation
        if mp.am_master(): self.master_logger.info("Running input FDTD simulation")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

        # get eigenmode coefficients
        if mp.am_master(): self.master_logger.info("Getting input FDTD simulation eigenmode coefficients:")
        src_coeff = sim.get_eigenmode_coefficients(in_flux,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
        input_coeffs = []
        for mode_num in range(1,self.hp.sim.num_modes_MMI+1):
            coeff = sim.get_eigenmode_coefficients(out_flux,[mode_num],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
            s = coeff.alpha[:,:,0]/src_coeff.alpha[:,:,0]
            input_coeffs.append(s[0,:])
            if mp.am_master(): self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(s[0,:])**2}, Angle: {np.angle(s[0,:], deg=True)}")
        input_coeffs = np.asarray(input_coeffs)
        return input_coeffs.T.astype(list)

    def mm_simulation(self,wavelengths):
        mm_coeffs = []
        # Load the GDSII file for the MMI splitter
        mmi_top = load_1x2_MMI(self.hp.geo.mm_length,0.0,self.hp.geo.waveguide_width,self.hp.geo.waveguide_offset,self.hp.sim.oxide_width,-self.hp.geo.taper_length-self.hp.sim.mm_monitor_position,self.hp.geo.taper_length,self.hp.geo.mm_width,self.hp.geo.taper_width)
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
                modes.append(mw.compute_modes_meep(cs, num_modes=num_mode))

            S, port_map = mw.compute_s_matrix(modes, cells=cells)
            mm_coeffs_1_wl = []
            for mode_num in range(self.hp.sim.num_modes_MMI):
                mm_coeffs_1_wl.append(S[mode_num,self.hp.sim.num_modes_MMI+mode_num])
                if mp.am_master(): self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(S[mode_num,self.hp.sim.num_modes_MMI+mode_num])**2}, Angle: {np.angle(S[mode_num,self.hp.sim.num_modes_MMI+mode_num], deg=True)}")
            mm_coeffs.append(mm_coeffs_1_wl)
        return mm_coeffs
    
    def output_simulation(self,wavelengths):
        if mp.am_master(): self.master_logger.info(f"Creating output FDTD simulation")
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
        geometry = [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2,self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2,-self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]

        # Create temp gds file for MMI geometry
        create_temp_gds_for_mmi_output(self.hp.geo.mm_length, self.hp.geo.taper_length, Sx, self.hp.sim.dpml, self.hp.sim.waveguide_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset)

        # Get MMI geometry from temp gds file
        mmi = mp.get_GDSII_prisms(Si, "gds/temp_out.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry += mmi

        # Create source MEEP simulation
        source_gaus = [
            mp.EigenModeSource(
                mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                eig_band=1,
                direction=mp.NO_DIRECTION,
                eig_kpoint=mp.Vector3(x=-1),
                size=mp.Vector3(0, self.hp.sim.oxide_width/2, self.hp.sim.oxide_height),
                center=[Sx/2-self.hp.sim.dpml-self.hp.sim.source_position, self.hp.sim.oxide_width/4, 0],
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
                center = mp.Vector3(Sx/2 - self.hp.sim.dpml - self.hp.sim.mm_monitor_position, (Sy - 2*self.hp.sim.dpml)/4, 0), 
                size = mp.Vector3(0, (Sy - 2*self.hp.sim.dpml)/2,Sz - 2*self.hp.sim.dpml)
            )
        )

        out_flux = sim.add_mode_monitor(
            freqs,
            mp.FluxRegion(
                center = mp.Vector3(Sx/2 - self.hp.sim.dpml - self.hp.sim.waveguide_length - self.hp.geo.taper_length - self.hp.sim.sm_monitor_position, 0, 0), 
                size = mp.Vector3(0, Sy - 2*self.hp.sim.dpml,Sz - 2*self.hp.sim.dpml)
            )
        )

        # Plot simulation geometry
        if mp.am_master(): self.master_logger.info("Plotting output FDTD simulation geometry")
        if mp.am_master():
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, Sy, 0)))
            plt.savefig(f'{self.hp.outputfolder}/layout/output_xy.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/output_xz.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/output_yz.png', dpi=300)
            plt.close("all")

        # run simulation
        if mp.am_master(): self.master_logger.info("Running output FDTD simulation")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

        # get eigenmode coefficients
        if mp.am_master(): self.master_logger.info("Getting output FDTD simulation eigenmode coefficients:")
        src_coeff = sim.get_eigenmode_coefficients(in_flux,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=-1))
        output_coeffs = []
        for mode_num in range(1,self.hp.sim.num_modes_MMI+1):
            coeff = sim.get_eigenmode_coefficients(out_flux,[mode_num],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=-1))
            s = coeff.alpha[:,:,1]/src_coeff.alpha[:,:,1]
            output_coeffs.append(s[0,:])
            if mp.am_master(): self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(s[0,:])**2}, Angle: {np.angle(s[0,:], deg=True)}")

        output_coeffs = np.asarray(output_coeffs)
        return output_coeffs.T.astype(list)

    def hybrid_simulation(self):
        wavelengths = np.linspace(self.hp.sim.wavelength_min,self.hp.sim.wavelength_max,self.hp.sim.wavelength_pts)
        input_coeffs = self.input_simulation(wavelengths)
        mm_coeffs = self.mm_simulation(wavelengths)
        output_coeffs = self.output_simulation(wavelengths)
        powers = []
        angles = []
        for point in range(len(wavelengths)):
            input_coeff = input_coeffs[point]
            mm_coeff = mm_coeffs[point]
            output_coeff = output_coeffs[point]
            propagated_coeffs = [x * y for x, y in zip(input_coeff, mm_coeff)]
            propagated_coeffs = [x * y for x, y in zip(propagated_coeffs, output_coeff)]
            powers.append(np.abs(sum(propagated_coeffs))**2)
            angles.append(np.angle(sum(propagated_coeffs), deg=True))
            if mp.am_master(): self.master_logger.info(f"Power: {powers[point]}, Angle: {angles[point]}")

        data_to_save = {
            'input_coeffs': input_coeffs,
            'mm_coeffs': mm_coeffs,
            'output_coeffs': output_coeffs,
            'powers': powers,
            'angles': angles,
            'wavelengths': wavelengths,
        }

        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)
        
        return powers, wavelengths

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
        geometry = [mp.Block(center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2,0,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2,self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]
        geometry += [mp.Block(center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2,-self.hp.geo.waveguide_offset/2,0), size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height), material = Si)]

        # Create temp gds file for MMI geometry
        create_temp_gds_for_full_mmi(self.hp.geo.mm_length, self.hp.geo.taper_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset)

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
                size=mp.Vector3(0, self.hp.sim.oxide_width, self.hp.sim.oxide_height),
                center=[-Sx/2+self.hp.sim.dpml+self.hp.sim.source_position, 0, 0],
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
                center = mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position, 0, 0), 
                size = mp.Vector3(0, Sy - 2*self.hp.sim.dpml,Sz - 2*self.hp.sim.dpml)
            )
        )

        out_flux = sim.add_mode_monitor(
            freqs,
            mp.FluxRegion(
                center = mp.Vector3(Sx/2 - self.hp.sim.dpml - self.hp.sim.sm_monitor_position, (Sy - 2*self.hp.sim.dpml)/4, 0), 
                size = mp.Vector3(0, (Sy - 2*self.hp.sim.dpml)/2,Sz - 2*self.hp.sim.dpml)
            )
        )

        # Plot simulation geometry
        if mp.am_master(): self.master_logger.info("Plotting FDTD simulation geometry")
        if mp.am_master():
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, Sy, 0)))
            plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_xy.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_xz.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_yz.png', dpi=300)
            plt.close("all")

        # run simulation
        if mp.am_master(): self.master_logger.info("Running FDTD simulation")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

        # get eigenmode coefficients
        if mp.am_master(): self.master_logger.info("Getting input FDTD simulation eigenmode coefficients:")
        src_coeff = sim.get_eigenmode_coefficients(in_flux,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
        coeff = sim.get_eigenmode_coefficients(out_flux,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
        s = coeff.alpha[:,:,0]/src_coeff.alpha[:,:,0]
        s = s[0,:]
        powers = []
        angles = []
        for point in range(len(wavelengths)):
            powers.append(np.abs(s[point])**2)
            angles.append(np.angle(s[point], deg=True))
            if mp.am_master(): self.master_logger.info(f"Power: {powers[point]}, Angle: {angles[point]}")

        data_to_save = {
            'coeffs': s,
            'powers': powers,
            'angles': angles,
            'wavelengths': wavelengths,
        }

        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)

        return powers, wavelengths
    
    def eme_simulation(self):
        wavelengths = np.linspace(self.hp.sim.wavelength_min,self.hp.sim.wavelength_max,self.hp.sim.wavelength_pts)
        # Load the GDSII file for the MMI splitter
        mmi_top = load_1x2_MMI(self.hp.geo.mm_length,self.hp.sim.waveguide_length,self.hp.geo.waveguide_width,self.hp.geo.waveguide_offset,self.hp.sim.oxide_width,0.0,self.hp.geo.taper_length,self.hp.geo.mm_width,self.hp.geo.taper_width)
        if mp.am_master():
            plt.figure()
            mmi_top.plot()
            plt.savefig(f'{self.hp.outputfolder}/layout/eme_gds.png', dpi=300)
            plt.close("all")


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
        powers = []
        coeffs = []
        angles = []
        for wavelength in wavelengths:
            if mp.am_master(): self.master_logger.info(f"Starting EME simulation at wavelength {wavelength}")
            
            gf_environment = mw.Environment(wl=wavelength, T=25.0)
            css = [mw.CrossSection.from_cell(cell=cell, env=gf_environment) for cell in cells]

            if mp.am_master():
                for cs_num, cs in enumerate(css):
                    plt.clf()
                    plt.figure()
                    mw.visualize(cs, show=False)
                    plt.savefig(f'{self.hp.outputfolder}/layout/cross_sections/cs_{cs_num}.png', dpi=300)
                    plt.close("all")


            # Performe mode simulations
            if mp.am_master(): self.master_logger.info(f"Performing mode simulations")
            modes = []
            for cell_num, (cs, num_mode) in enumerate(zip(css,num_modes)):
                if mp.am_master(): self.master_logger.info(f"Computing mode {num_mode} of cell {cell_num}")
                modes.append(mw.compute_modes_meep(cs, num_modes=num_mode))
            S, port_map = mw.compute_s_matrix(modes, cells=cells)
            coeffs.append(S[0,1])
            powers.append(np.abs(S[0,1])**2)
            angles.append(np.angle(S[0,1],deg=True))

            if wavelength == wavelengths[int(np.floor(len(wavelengths)/2))]:
                plot_eme_fields(self.hp,cells,modes)
                if mp.am_master(): self.master_logger.info(f"Transmission at center wavelength: {powers[-1]}")

        data_to_save = {
            'coeffs': coeffs,
            'powers': powers,
            'angles': angles,
            'wavelengths': wavelengths,
        }

        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)
        return powers, wavelengths

    def run(self):
        if self.hp.sim.sim_type == "HYBRID":
            powers, wavelengths = self.hybrid_simulation()
        elif self.hp.sim.sim_type == "FDTD":
            powers, wavelengths = self.fdtd_simulation()
        elif self.hp.sim.sim_type == "EME":
            powers, wavelengths = self.eme_simulation()

        if mp.am_master():
            plt.figure()
            plt.plot(wavelengths,powers)
            plt.xlabel("Wavelength (micron)")
            plt.ylabel("Transmission")
            plt.savefig(f'{self.hp.outputfolder}/performance.png', dpi=300)
            plt.close("all")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file", required=True)
    args = parser.parse_args()
    mmi = MMI_1x2(args)

    mmi.run()