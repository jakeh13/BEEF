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
    make_output_folders
)
from engine.helpers.meep import(
     meep_material_from_meow,
     create_temp_gds_for_mmi_1x4_input,
     create_temp_gds_for_mmi_1x4_output,
     create_temp_gds_for_full_mmi_1x4
)
from engine.helpers.meow import(
    load_1x4_MMI,
    stack_extrusions_GF
)
from engine.plotting.meow import(
    plot_eme_fields
)

OUTPUT_OFFSETS_FUNC = lambda wo: [3*wo/2, wo/2, -wo/2, -3*wo/2]

class MMI_1x4:
    def __init__(self, args):
        self.config = args.config
        self.hp: HyperParams = dataconf.load(self.config, HyperParams)
        make_output_folders(self.hp.outputfolder)

        s = json.dumps(asdict(self.hp), indent=4, default=lambda x: x.name)
        with open(f"{self.hp.outputfolder}/config.json", "w") as f:
            print(s, file=f)

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

    def input_simulation(self, wavelengths):
        if mp.am_master(): self.master_logger.info(f"Creating input FDTD simulation")
        Sx = self.hp.geo.taper_length + self.hp.sim.waveguide_length + self.hp.sim.mm_sim_length + 2*self.hp.sim.dpml
        Sy = self.hp.sim.oxide_width + 2*self.hp.sim.dpml
        Sz = self.hp.sim.oxide_height + 2*self.hp.sim.dpml
        wavelength_center = (self.hp.sim.wavelength_max - self.hp.sim.wavelength_min)/2 + self.hp.sim.wavelength_min
        freqs = 1/wavelengths.astype(list)

        Si = meep_material_from_meow(mw.silicon.n, mw.silicon.params['wl'], wavelength_center)
        SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'], wavelength_center)

        geometry = [mp.Block(
            center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2, 0, 0),
            size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height),
            material=Si
        )]

        create_temp_gds_for_mmi_1x4_input(
            self.hp.geo.mm_length, self.hp.geo.taper_length, Sx, self.hp.sim.dpml,
            self.hp.sim.waveguide_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width,
            self.hp.geo.taper_width, self.hp.geo.waveguide_offset
        )
        mmi = mp.get_GDSII_prisms(Si, "gds/temp.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry += mmi

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

        sim = mp.Simulation(
            cell_size=mp.Vector3(Sx, Sy, Sz),
            boundary_layers=[mp.PML(self.hp.sim.dpml)],
            geometry=geometry,
            sources=source_gaus,
            eps_averaging=True,
            resolution=self.hp.sim.res,
            default_material=SiO2
        )

        in_flux = sim.add_mode_monitor(
            freqs,
            mp.FluxRegion(
                center=mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position, 0, 0),
                size=mp.Vector3(0, Sy - 2*self.hp.sim.dpml, Sz - 2*self.hp.sim.dpml)
            )
        )
        out_flux = sim.add_mode_monitor(
            freqs,
            mp.FluxRegion(
                center=mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.waveguide_length + self.hp.geo.taper_length + self.hp.sim.mm_monitor_position, 0, 0),
                size=mp.Vector3(0, Sy - 2*self.hp.sim.dpml, Sz - 2*self.hp.sim.dpml)
            )
        )

        if mp.am_master(): self.master_logger.info("Plotting input FDTD simulation geometry")
        if mp.am_master():
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, 0)))
            plt.savefig(f'{self.hp.outputfolder}/layout/input_xy.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, 0, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/input_xz.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(0, Sy, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/input_yz.png', dpi=300)
            plt.close("all")

        if mp.am_master(): self.master_logger.info("Running input FDTD simulation")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol=1e-11))

        if mp.am_master(): self.master_logger.info("Getting input FDTD simulation eigenmode coefficients:")
        src_coeff = sim.get_eigenmode_coefficients(in_flux, [1], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
        input_coeffs = []
        for mode_num in range(1, self.hp.sim.num_modes_MMI+1):
            coeff = sim.get_eigenmode_coefficients(out_flux, [mode_num], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
            s = coeff.alpha[:,:,0]/src_coeff.alpha[:,:,0]
            input_coeffs.append(s[0,:])
            if mp.am_master(): self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(s[0,:])**2}, Angle: {np.angle(s[0,:], deg=True)}")
        input_coeffs = np.asarray(input_coeffs)
        return input_coeffs.T.astype(list)

    def mm_simulation(self, wavelengths):
        mm_coeffs = []
        mmi_gds = load_1x4_MMI(
            self.hp.geo.mm_length, 0.0, self.hp.geo.waveguide_width, self.hp.geo.waveguide_offset,
            self.hp.sim.oxide_width, -self.hp.geo.taper_length-self.hp.sim.mm_monitor_position,
            self.hp.geo.taper_length, self.hp.geo.mm_width, self.hp.geo.taper_width, 0
        )
        extrusion_rules = stack_extrusions_GF(self.hp.geo.waveguide_height, self.hp.sim.oxide_height)
        structs = mw.extrude_gds(mmi_gds, extrusion_rules)
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
            if mp.am_master(): self.master_logger.info(f"Performing mode simulations")
            modes = []
            for cell_num, (cs, num_mode) in enumerate(zip(css, num_modes)):
                if mp.am_master(): self.master_logger.info(f"Computing mode {num_mode} of cell {cell_num}")
                modes.append(compute_modes(cs, num_modes=num_mode))
            S, port_map = mw.compute_s_matrix(modes, cells=cells)
            mm_coeffs_1_wl = []
            for mode_num in range(self.hp.sim.num_modes_MMI):
                mm_coeffs_1_wl.append(S[mode_num, self.hp.sim.num_modes_MMI+mode_num])
                if mp.am_master(): self.master_logger.info(f"Mode {mode_num}, Power: {np.abs(S[mode_num,self.hp.sim.num_modes_MMI+mode_num])**2}, Angle: {np.angle(S[mode_num,self.hp.sim.num_modes_MMI+mode_num], deg=True)}")
            mm_coeffs.append(mm_coeffs_1_wl)
        return mm_coeffs

    def output_simulation(self, wavelengths):
        output_offsets = OUTPUT_OFFSETS_FUNC(self.hp.geo.waveguide_offset)
        all_output_coeffs = []

        Sx = self.hp.geo.taper_length + self.hp.sim.waveguide_length + self.hp.sim.mm_sim_length + 2*self.hp.sim.dpml
        Sy = self.hp.sim.oxide_width + 2*self.hp.sim.dpml
        Sz = self.hp.sim.oxide_height + 2*self.hp.sim.dpml
        wavelength_center = (self.hp.sim.wavelength_max - self.hp.sim.wavelength_min)/2 + self.hp.sim.wavelength_min
        freqs = 1/wavelengths.astype(list)

        Si = meep_material_from_meow(mw.silicon.n, mw.silicon.params['wl'], wavelength_center)
        SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'], wavelength_center)

        # Build geometry once (same for all 4 port sims)
        create_temp_gds_for_mmi_1x4_output(
            self.hp.geo.mm_length, self.hp.geo.taper_length, Sx, self.hp.sim.dpml,
            self.hp.sim.waveguide_length, self.hp.geo.waveguide_width, self.hp.geo.mm_width,
            self.hp.geo.taper_width, self.hp.geo.waveguide_offset
        )

        for port_idx, yo in enumerate(output_offsets):
            if mp.am_master(): self.master_logger.info(f"Creating output FDTD simulation for port {port_idx}")

            geometry = [mp.Block(
                center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2, 0, 0),
                size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height),
                material=Si
            )]
            for out_yo in output_offsets:
                geometry += [mp.Block(
                    center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2, out_yo, 0),
                    size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height),
                    material=Si
                )]
            mmi = mp.get_GDSII_prisms(Si, "gds/temp_out.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
            geometry += mmi

            source_gaus = [
                mp.EigenModeSource(
                    mp.GaussianSource(frequency=1/wavelength_center, fwidth=0.4/wavelength_center),
                    eig_band=1,
                    direction=mp.NO_DIRECTION,
                    eig_kpoint=mp.Vector3(x=-1),
                    size=mp.Vector3(0, self.hp.geo.waveguide_offset, self.hp.sim.oxide_height),
                    center=[Sx/2 - self.hp.sim.dpml - self.hp.sim.source_position, yo, 0],
                )
            ]

            sim = mp.Simulation(
                cell_size=mp.Vector3(Sx, Sy, Sz),
                boundary_layers=[mp.PML(self.hp.sim.dpml)],
                geometry=geometry,
                sources=source_gaus,
                eps_averaging=True,
                resolution=self.hp.sim.res,
                default_material=SiO2
            )

            in_flux = sim.add_mode_monitor(
                freqs,
                mp.FluxRegion(
                    center=mp.Vector3(Sx/2 - self.hp.sim.dpml - self.hp.sim.mm_monitor_position, yo, 0),
                    size=mp.Vector3(0, self.hp.geo.waveguide_offset, Sz - 2*self.hp.sim.dpml)
                )
            )
            out_flux = sim.add_mode_monitor(
                freqs,
                mp.FluxRegion(
                    center=mp.Vector3(Sx/2 - self.hp.sim.dpml - self.hp.sim.waveguide_length - self.hp.geo.taper_length - self.hp.sim.sm_monitor_position, 0, 0),
                    size=mp.Vector3(0, Sy - 2*self.hp.sim.dpml, Sz - 2*self.hp.sim.dpml)
                )
            )

            if mp.am_master():
                plt.figure()
                sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, 0)))
                plt.savefig(f'{self.hp.outputfolder}/layout/output_{port_idx}_xy.png', dpi=300)
                plt.close("all")

            if mp.am_master(): self.master_logger.info(f"Running output FDTD simulation for port {port_idx}")
            sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol=1e-11))

            if mp.am_master(): self.master_logger.info(f"Getting output FDTD eigenmode coefficients for port {port_idx}:")
            src_coeff = sim.get_eigenmode_coefficients(in_flux, [1], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=-1))
            output_coeffs = []
            for mode_num in range(1, self.hp.sim.num_modes_MMI+1):
                coeff = sim.get_eigenmode_coefficients(out_flux, [mode_num], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=-1))
                s = coeff.alpha[:,:,1]/src_coeff.alpha[:,:,1]
                output_coeffs.append(s[0,:])
                if mp.am_master(): self.master_logger.info(f"Port {port_idx}, Mode {mode_num}, Power: {np.abs(s[0,:])**2}, Angle: {np.angle(s[0,:], deg=True)}")
            output_coeffs = np.asarray(output_coeffs)
            all_output_coeffs.append(output_coeffs.T.astype(list))

            sim.reset_meep()

        return all_output_coeffs

    def hybrid_simulation(self):
        wavelengths = np.linspace(self.hp.sim.wavelength_min, self.hp.sim.wavelength_max, self.hp.sim.wavelength_pts)
        input_coeffs = self.input_simulation(wavelengths)
        mm_coeffs = self.mm_simulation(wavelengths)
        all_output_coeffs = self.output_simulation(wavelengths)

        powers_ports = [[] for _ in range(4)]
        angles_ports = [[] for _ in range(4)]

        for point in range(len(wavelengths)):
            input_coeff = input_coeffs[point]
            mm_coeff = mm_coeffs[point]
            for port_idx in range(4):
                output_coeff = all_output_coeffs[port_idx][point]
                propagated = [x * y * z for x, y, z in zip(input_coeff, mm_coeff, output_coeff)]
                powers_ports[port_idx].append(np.abs(sum(propagated))**2)
                angles_ports[port_idx].append(np.angle(sum(propagated), deg=True))
                if mp.am_master(): self.master_logger.info(f"Port {port_idx}: Power: {powers_ports[port_idx][point]}, Angle: {angles_ports[port_idx][point]}")

        data_to_save = {
            'input_coeffs': input_coeffs,
            'mm_coeffs': mm_coeffs,
            'all_output_coeffs': all_output_coeffs,
            'powers_ports': powers_ports,
            'angles_ports': angles_ports,
            'wavelengths': wavelengths,
        }
        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)

        return powers_ports, wavelengths

    def fdtd_simulation(self):
        if mp.am_master(): self.master_logger.info(f"Creating full FDTD simulation")
        Sx = 2*self.hp.geo.taper_length + 2*self.hp.sim.waveguide_length + self.hp.geo.mm_length + 2*self.hp.sim.dpml
        Sy = self.hp.sim.oxide_width + 2*self.hp.sim.dpml
        Sz = self.hp.sim.oxide_height + 2*self.hp.sim.dpml
        wavelength_center = (self.hp.sim.wavelength_max - self.hp.sim.wavelength_min)/2 + self.hp.sim.wavelength_min
        wavelengths = np.linspace(self.hp.sim.wavelength_min, self.hp.sim.wavelength_max, self.hp.sim.wavelength_pts)
        freqs = 1/wavelengths.astype(list)

        Si = meep_material_from_meow(mw.silicon.n, mw.silicon.params['wl'], wavelength_center)
        SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'], wavelength_center)

        output_offsets = OUTPUT_OFFSETS_FUNC(self.hp.geo.waveguide_offset)

        # Input waveguide
        geometry = [mp.Block(
            center=mp.Vector3(-Sx/2 + self.hp.sim.waveguide_length/2 + self.hp.sim.dpml/2, 0, 0),
            size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height),
            material=Si
        )]
        # 4 output waveguides
        for yo in output_offsets:
            geometry += [mp.Block(
                center=mp.Vector3(Sx/2 - self.hp.sim.waveguide_length/2 - self.hp.sim.dpml/2, yo, 0),
                size=mp.Vector3(self.hp.sim.waveguide_length + self.hp.sim.dpml, self.hp.geo.waveguide_width, self.hp.geo.waveguide_height),
                material=Si
            )]

        create_temp_gds_for_full_mmi_1x4(
            self.hp.geo.mm_length, self.hp.geo.taper_length, self.hp.geo.waveguide_width,
            self.hp.geo.mm_width, self.hp.geo.taper_width, self.hp.geo.waveguide_offset
        )
        mmi = mp.get_GDSII_prisms(Si, "gds/temp_full.gds", 2, -self.hp.geo.waveguide_height/2, self.hp.geo.waveguide_height/2)
        geometry += mmi

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

        sim = mp.Simulation(
            cell_size=mp.Vector3(Sx, Sy, Sz),
            boundary_layers=[mp.PML(self.hp.sim.dpml)],
            geometry=geometry,
            sources=source_gaus,
            eps_averaging=True,
            resolution=self.hp.sim.res,
            default_material=SiO2
        )

        in_flux = sim.add_mode_monitor(
            freqs,
            mp.FluxRegion(
                center=mp.Vector3(-Sx/2 + self.hp.sim.dpml + self.hp.sim.sm_monitor_position, 0, 0),
                size=mp.Vector3(0, Sy - 2*self.hp.sim.dpml, Sz - 2*self.hp.sim.dpml)
            )
        )
        out_fluxes = []
        for yo in output_offsets:
            out_fluxes.append(sim.add_mode_monitor(
                freqs,
                mp.FluxRegion(
                    center=mp.Vector3(Sx/2 - self.hp.sim.dpml - self.hp.sim.sm_monitor_position, yo, 0),
                    size=mp.Vector3(0, self.hp.geo.waveguide_offset, Sz - 2*self.hp.sim.dpml)
                )
            ))

        if mp.am_master(): self.master_logger.info("Plotting FDTD simulation geometry")
        if mp.am_master():
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, Sy, 0)))
            plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_xy.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(Sx, 0, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_xz.png', dpi=300)
            plt.figure()
            sim.plot2D(output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(0, Sy, Sz)))
            plt.savefig(f'{self.hp.outputfolder}/layout/full_FDTD_yz.png', dpi=300)
            plt.close("all")

        if mp.am_master(): self.master_logger.info("Running FDTD simulation")
        sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol=1e-11))

        if mp.am_master(): self.master_logger.info("Getting FDTD eigenmode coefficients:")
        src_coeff = sim.get_eigenmode_coefficients(in_flux, [1], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))

        powers_ports = []
        angles_ports = []
        all_s = []
        for port_idx, out_flux in enumerate(out_fluxes):
            coeff = sim.get_eigenmode_coefficients(out_flux, [1], eig_parity=mp.NO_PARITY, kpoint_func=lambda a, b: mp.Vector3(x=1))
            s = coeff.alpha[:,:,0]/src_coeff.alpha[:,:,0]
            s = s[0,:]
            all_s.append(s)
            powers = [np.abs(s[pt])**2 for pt in range(len(wavelengths))]
            angles = [np.angle(s[pt], deg=True) for pt in range(len(wavelengths))]
            powers_ports.append(powers)
            angles_ports.append(angles)
            for pt in range(len(wavelengths)):
                if mp.am_master(): self.master_logger.info(f"Port {port_idx}: Power: {powers[pt]}, Angle: {angles[pt]}")

        data_to_save = {
            'all_s': all_s,
            'powers_ports': powers_ports,
            'angles_ports': angles_ports,
            'wavelengths': wavelengths,
        }
        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)

        return powers_ports, wavelengths

    def eme_simulation(self):
        wavelengths = np.linspace(self.hp.sim.wavelength_min, self.hp.sim.wavelength_max, self.hp.sim.wavelength_pts)
        extrusion_rules = stack_extrusions_GF(self.hp.geo.waveguide_height, self.hp.sim.oxide_height)
        taper_Ls = np.array([self.hp.geo.taper_length / self.hp.sim.num_cells_taper for _ in range(self.hp.sim.num_cells_taper)])
        Ls = np.concatenate(([self.hp.sim.waveguide_length], taper_Ls, [self.hp.geo.mm_length], taper_Ls, [self.hp.sim.waveguide_length]))
        num_modes = [1] + [self.hp.sim.num_modes_taper for _ in range(self.hp.sim.num_cells_taper)] + [self.hp.sim.num_modes_MMI] + [self.hp.sim.num_modes_taper for _ in range(self.hp.sim.num_cells_taper)] + [1]

        # Build ouput cells for each of the 4 output ports
        output_cells = []
        for output_idx in range(4):
            mmi_gds = load_1x4_MMI(
                self.hp.geo.mm_length, self.hp.sim.waveguide_length, self.hp.geo.waveguide_width,
                self.hp.geo.waveguide_offset, self.hp.sim.oxide_width, 0.0,
                self.hp.geo.taper_length, self.hp.geo.mm_width, self.hp.geo.taper_width, output_idx
            )
            if mp.am_master():
                plt.figure()
                mmi_gds.plot()
                plt.savefig(f'{self.hp.outputfolder}/layout/eme_gds_port{output_idx}.png', dpi=300)
                plt.close("all")
            structs = mw.extrude_gds(mmi_gds, extrusion_rules)
            cells = mw.create_cells(
                structures=structs,
                mesh=mw.Mesh2D(
                    x=np.linspace(-self.hp.sim.oxide_width/2, self.hp.sim.oxide_width/2, int(self.hp.sim.res*self.hp.sim.oxide_width) + 1),
                    y=np.linspace(-self.hp.sim.oxide_height/2, self.hp.sim.oxide_height/2, int(self.hp.sim.res*self.hp.sim.oxide_height) + 1),
                ),
                Ls=Ls,
            )
            if output_idx == 0:
                input_cells = cells[0:-1]
            output_cells.append(cells[-1])

        powers_ports = [[] for _ in range(4)]
        coeffs_ports = [[] for _ in range(4)]
        angles_ports = [[] for _ in range(4)]

        for wavelength in wavelengths:
            if mp.am_master(): self.master_logger.info(f"Starting EME simulation at wavelength {wavelength}")
            gf_environment = mw.Environment(wl=wavelength, T=25.0)
            
            css = [mw.CrossSection.from_cell(cell=cell, env=gf_environment) for cell in input_cells]

            if mp.am_master():
                for cs_num, cs in enumerate(css):
                    plt.clf()
                    plt.figure()
                    mw.visualize(cs, show=False)
                    plt.savefig(f'{self.hp.outputfolder}/layout/cross_sections/cs_input_{cs_num}.png', dpi=300)
                    plt.close("all")

            if mp.am_master(): self.master_logger.info(f"Performing input mode simulations")
            input_modes = []
            for cell_num, (cs, num_mode) in enumerate(zip(css, num_modes)):
                if mp.am_master(): self.master_logger.info(f"Computing mode {num_mode} of cell {cell_num}")
                input_modes.append(compute_modes(cs, num_modes=num_mode))

            output_modes = []
            for output_idx, out_cell in enumerate(output_cells):
                css = [mw.CrossSection.from_cell(cell=out_cell, env=gf_environment)]
                if mp.am_master():
                    for cs_num, cs in enumerate(css):
                        plt.clf()
                        plt.figure()
                        mw.visualize(cs, show=False)
                        plt.savefig(f'{self.hp.outputfolder}/layout/cross_sections/cs_out_{output_idx}_{cs_num}.png', dpi=300)
                        plt.close("all")

                if mp.am_master(): self.master_logger.info(f"Performing output mode simulations for port {output_idx}")
                
                for cell_num, (cs, num_mode) in enumerate(zip(css, num_modes)):
                    if mp.am_master(): self.master_logger.info(f"Computing mode {num_mode} of cell {cell_num}")
                    output_modes.append(compute_modes(cs, num_modes=num_mode))

            for output_idx in range(4):
                full_modes = input_modes + output_modes[output_idx:output_idx+1]
                full_cells = input_cells + [output_cells[output_idx]]   
                S, port_map = mw.compute_s_matrix(full_modes, cells=full_cells)
                coeffs_ports[output_idx].append(S[0,1])
                powers_ports[output_idx].append(np.abs(S[0,1])**2)
                angles_ports[output_idx].append(np.angle(S[0,1], deg=True))

                if wavelength == wavelengths[int(np.floor(len(wavelengths)/2))]:
                    plot_eme_fields(self.hp, full_cells, full_modes, output_idx)
                    if mp.am_master(): self.master_logger.info(f"Port {output_idx} transmission at center wavelength: {powers_ports[output_idx][-1]}")

        data_to_save = {
            'coeffs_ports': coeffs_ports,
            'powers_ports': powers_ports,
            'angles_ports': angles_ports,
            'wavelengths': wavelengths,
        }
        with open(f'{self.hp.outputfolder}/data/sim_results.pkl', 'wb') as f:
            pickle.dump(data_to_save, f)

        return powers_ports, wavelengths

    def run(self):
        if self.hp.sim.sim_type == "HYBRID":
            powers_ports, wavelengths = self.hybrid_simulation()
        elif self.hp.sim.sim_type == "FDTD":
            powers_ports, wavelengths = self.fdtd_simulation()
        elif self.hp.sim.sim_type == "EME":
            powers_ports, wavelengths = self.eme_simulation()

        if mp.am_master():
            plt.figure()
            for port_idx, powers in enumerate(powers_ports):
                plt.plot(wavelengths, powers, label=f'Port {port_idx}')
            plt.xlabel("Wavelength (micron)")
            plt.ylabel("Transmission")
            plt.legend()
            plt.savefig(f'{self.hp.outputfolder}/performance.png', dpi=300)
            plt.close("all")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", help="config file", required=True)
    args = parser.parse_args()
    mmi = MMI_1x4(args)

    mmi.run()
