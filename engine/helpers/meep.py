import numpy as np
import meep as mp
import gdsfactory as gf

gf.gpdk.PDK.activate()

# Get MEEP material from index vs wavelength lists used in MEOW
def meep_material_from_meow(ns, wls, wl):
    idx = np.argmin(np.abs(wls - wl))
    n = ns[idx]
    return mp.Medium(epsilon=n**2)

def create_MMI_1x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset,name):
    mmi = gf.Component(name)
    # multimode section
    mmi.add_polygon([(-mm_length/2,-mm_width/2),
                     (mm_length/2,-mm_width/2),
                     (mm_length/2,mm_width/2),
                     (-mm_length/2,mm_width/2)], 
                     layer=(2,0))
    # input taper
    mmi.add_polygon([(-mm_length/2-taper_length,-waveguide_width/2),
                     (-mm_length/2,-taper_width/2),
                     (-mm_length/2,taper_width/2),
                     (-mm_length/2-taper_length,waveguide_width/2)], 
                     layer=(2,0))
    # output tapers
    mmi.add_polygon([(mm_length/2,-taper_width/2+waveguide_offset),
                     (mm_length/2+taper_length,-waveguide_width/2+waveguide_offset),
                     (mm_length/2+taper_length,waveguide_width/2+waveguide_offset),
                     (mm_length/2,taper_width/2+waveguide_offset)], 
                     layer=(2,0))
    mmi.add_polygon([(mm_length/2,-taper_width/2-waveguide_offset),
                     (mm_length/2+taper_length,-waveguide_width/2-waveguide_offset),
                     (mm_length/2+taper_length,waveguide_width/2-waveguide_offset),
                     (mm_length/2,taper_width/2-waveguide_offset)], 
                     layer=(2,0))
    return mmi

# Function to generate temp gds file for MMI geometry (to shift the MMI to the correct simulation coordinates)
def create_temp_gds_for_mmi_input(mm_length, taper_length, Sx, dpml, waveguide_length, waveguide_width, mm_width, taper_width, waveguide_offset):
    mmi_total_length = mm_length + 2*taper_length
    temp = create_MMI_1x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset/2,"temp")
    temp.move((mmi_total_length/2 - Sx/2 + dpml + waveguide_length, 0))
    # save as temp.gds
    temp.write_gds(gdsdir='gds')

# Function to generate temp gds file for MMI geometry (to shift the MMI to the correct simulation coordinates)
def create_temp_gds_for_mmi_output(mm_length, taper_length, Sx, dpml, waveguide_length, waveguide_width, mm_width, taper_width, waveguide_offset):
    mmi_total_length = mm_length + 2*taper_length
    temp = create_MMI_1x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset/2,"temp_out")
    temp.move((-mmi_total_length/2 + Sx/2 - dpml - waveguide_length, 0))
    # save as temp.gds
    temp.write_gds(gdsdir='gds')

# Function to generate temp gds file for MMI geometry (to shift the MMI to the correct simulation coordinates)
def create_temp_gds_for_full_mmi(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset):
    temp = create_MMI_1x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset/2,"temp_full")
    # save as temp.gds
    temp.write_gds(gdsdir='gds')