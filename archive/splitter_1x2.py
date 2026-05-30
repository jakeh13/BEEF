# Import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pickle
import meep as mp
import meow as mw

from engine.helpers.meep import meep_material_from_meow, create_temp_gds_for_mmi_input, create_temp_gds_for_mmi_output
from engine.helpers.general import make_output_folders
from engine.helpers.meow import stack_extrusions_GF, load_1x2_MMI

# Set directory for output folders
output_dir = "out/test_runs"
make_output_folders(output_dir)

# MMI parameters (dependent on gds design)
mm_length = 4.0
taper_length = 1.0
mmi_total_length = mm_length + 2*taper_length
waveguide_width = 0.35
waveguide_offset = 1.0
waveguide_height = 0.16
gds_name = "MMI_splitter"

# Used to set simulation dimensions
oxide_height = 1.5
oxide_width = 5.0

# Simulation parameters
res = 10 # Pixels per micron
wavelength = 1.31 # Wavelength in microns
num_modes_MMI = 15 # Number of modes computed in the MMI section
dpml = 0.5 # Thickness of PML layer in microns
mm_sim_length = 1.0 # Length of MMI section to simulate (should be less than mm_length to avoid taper)
waveguide_length = 0.5 # Length of waveguide between pml and taper
source_position = 0.0 # x_position from the pml where the input source is placed
input_monitor_position = 0.25 # x_position from the pml where the input mode monitor is placed
output_monitor_position = 0.25 # x_position from the end of the taper where the output mode monitor is placed

Sx = taper_length + waveguide_length + mm_sim_length + 2*dpml  
Sy = oxide_width + 2*dpml  
Sz = oxide_height + 2*dpml    

# materials
Si = meep_material_from_meow(mw.silicon.n, mw.silicon.params['wl'],wavelength)
SiO2 = meep_material_from_meow(mw.silicon_oxide.n, mw.silicon_oxide.params['wl'],wavelength)

# Set waveguide geometry
geometry = [mp.Block(center=mp.Vector3(-Sx/2 + waveguide_length/2 + dpml/2,0,0), size=mp.Vector3(waveguide_length + dpml, waveguide_width, waveguide_height), material = Si)]

# Create temp gds file for MMI geometry
create_temp_gds_for_mmi_input(mmi_total_length, Sx, dpml, waveguide_length, gds_name)

# Get MMI geometry from temp gds file
mmi = mp.get_GDSII_prisms(Si, "gds/temp.gds", 2, -waveguide_height/2, waveguide_height/2)
geometry += mmi

# Create source MEEP simulation
source_gaus = [
    mp.EigenModeSource(
        mp.GaussianSource(frequency=1/wavelength, fwidth=0.4/wavelength),
        eig_band=1,
        direction=mp.NO_DIRECTION,
        eig_kpoint=mp.Vector3(x=1),
        size=mp.Vector3(0, oxide_width, oxide_height),
        center=[-Sx/2+dpml+source_position, 0, 0],
    )
]

# Create MEEP simulation
sim = mp.Simulation(
    cell_size=mp.Vector3(Sx, Sy, Sz),
    boundary_layers=[mp.PML(dpml)],
    geometry=geometry,
    sources=source_gaus,
    eps_averaging=True,
    resolution=res,
    default_material=SiO2
)

# Add monitors

in_flux = sim.add_mode_monitor(
    [1/wavelength],
    mp.FluxRegion(
        center = mp.Vector3(-Sx/2 + dpml + input_monitor_position, 0, 0), 
        size = mp.Vector3(0, Sy - 2*dpml,Sz - 2*dpml)
    )
)

out_flux = sim.add_mode_monitor(
    [1/wavelength],
    mp.FluxRegion(
        center = mp.Vector3(-Sx/2 + dpml + waveguide_length + taper_length + output_monitor_position, 0, 0), 
        size = mp.Vector3(0, Sy - 2*dpml,Sz - 2*dpml)
    )
)

# Plot simulation geometry
plt.figure()
sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, Sy, 0)))
plt.savefig(f'{output_dir}/layout/input_xy.png', dpi=300)
plt.figure()
sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
plt.savefig(f'{output_dir}/layout/input_xz.png', dpi=300)
plt.figure()
sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
plt.savefig(f'{output_dir}/layout/input_yz.png', dpi=300)

# run simulation
sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

# get eigenmode coefficients
src_coeff = sim.get_eigenmode_coefficients(in_flux,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
in_coeffs_list = []
for mode_num in range(1,num_modes_MMI+1):
    coeff = sim.get_eigenmode_coefficients(out_flux,[mode_num],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=1))
    s = coeff.alpha[:,:,0]/src_coeff.alpha[:,:,0]
    in_coeffs_list.append(s[0,0])

# Load the GDSII file for the MMI splitter
mmi_top = load_1x2_MMI(mm_length + 2*taper_length,0.0,waveguide_width,waveguide_offset,oxide_width,-taper_length-output_monitor_position,gds_name)
extrusion_rules = stack_extrusions_GF(waveguide_height, oxide_height)
structs = mw.extrude_gds(mmi_top, extrusion_rules)

# EME cell for the MMI region
Ls = [mm_length - 2*output_monitor_position]
num_modes =[num_modes_MMI]
cells = mw.create_cells(
    structures=structs,
    mesh=mw.Mesh2D(
        x=np.linspace(-oxide_width/2, oxide_width/2, int(res*oxide_width) + 1),
        y=np.linspace(-oxide_height/2, oxide_height/2, int(res*oxide_height) + 1),
    ),
    Ls=Ls,
)
gf_environment = mw.Environment(wl=wavelength, T=25.0)
css = [mw.CrossSection.from_cell(cell=cell, env=gf_environment) for cell in cells]

# Performe mode simulations
modes = []
for cell_num, (cs, num_mode) in enumerate(zip(css,num_modes)):
    modes.append(mw.compute_modes_meep(cs, num_modes=num_mode))

S, port_map = mw.compute_s_matrix(modes, cells=cells)
mm_coeffs_list = []
for mode_num in range(num_modes_MMI):
    mm_coeffs_list.append(S[mode_num,num_modes_MMI+mode_num])

# Calculate propagated s-matrices
propagated_coeffs = [x * y for x, y in zip(in_coeffs_list, mm_coeffs_list)]


## Simulate output section

# Set waveguide geometry
geometry = [mp.Block(center=mp.Vector3(Sx/2 - waveguide_length/2 - dpml/2,waveguide_offset/2,0), size=mp.Vector3(waveguide_length + dpml, waveguide_width, waveguide_height), material = Si)]

# Create temp gds file for MMI geometry
create_temp_gds_for_mmi_output(mmi_total_length, Sx, dpml, waveguide_length, gds_name)

# Get MMI geometry from temp gds file
mmi = mp.get_GDSII_prisms(Si, "gds/temp_out.gds", 2, -waveguide_height/2, waveguide_height/2)
geometry += mmi

# Create source MEEP simulation
source_gaus = [
    mp.EigenModeSource(
        mp.GaussianSource(frequency=1/wavelength, fwidth=0.4/wavelength),
        eig_band=1,
        direction=mp.NO_DIRECTION,
        eig_kpoint=mp.Vector3(x=-1),
        size=mp.Vector3(0, oxide_width, oxide_height),
        center=[Sx/2-dpml-source_position, 0, 0],
    )
]

# Create MEEP simulation
sim = mp.Simulation(
    cell_size=mp.Vector3(Sx, Sy, Sz),
    boundary_layers=[mp.PML(dpml)],
    geometry=geometry,
    sources=source_gaus,
    eps_averaging=True,
    resolution=res,
    default_material=SiO2
)

# Add monitors

in_flux = sim.add_mode_monitor(
    [1/wavelength],
    mp.FluxRegion(
        center = mp.Vector3(Sx/2 - dpml - input_monitor_position, 0, 0), 
        size = mp.Vector3(0, Sy - 2*dpml,Sz - 2*dpml)
    )
)

out_flux = sim.add_mode_monitor(
    [1/wavelength],
    mp.FluxRegion(
        center = mp.Vector3(Sx/2 - dpml - waveguide_length - taper_length - output_monitor_position, 0, 0), 
        size = mp.Vector3(0, Sy - 2*dpml,Sz - 2*dpml)
    )
)

# Plot simulation geometry
plt.figure()
sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, Sy, 0)))
plt.savefig(f'{output_dir}/layout/output_xy.png', dpi=300)
plt.figure()
sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(Sx, 0, Sz)))
plt.savefig(f'{output_dir}/layout/output_xz.png', dpi=300)
plt.figure()
sim.plot2D(output_plane=mp.Volume(center = mp.Vector3(), size = mp.Vector3(0, Sy, Sz)))
plt.savefig(f'{output_dir}/layout/output_yz.png', dpi=300)

# run simulation
sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=20, tol = 1e-11))

# get eigenmode coefficients
src_coeffs = sim.get_eigenmode_coefficients(in_flux,[1],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=-1))
out_coeffs_list = []
for mode_num in range(1,num_modes_MMI+1):
    coeff = sim.get_eigenmode_coefficients(out_flux,[mode_num],eig_parity=mp.NO_PARITY, kpoint_func= lambda a, b: mp.Vector3(x=-1))
    s = coeff.alpha[:,:,1]/src_coeffs.alpha[:,:,1]
    out_coeffs_list.append(s[0,0])

# Calculate propagated s-matrices
propagated_coeffs_final = [x * y for x, y in zip(propagated_coeffs, out_coeffs_list)]

final_coeff = sum(propagated_coeffs_final)

mode_power = np.abs(final_coeff)**2
mode_angle = np.angle(final_coeff, deg=True)
print(f"Coeff {final_coeff}, Power: {mode_power}, Angle: {mode_angle}")