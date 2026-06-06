import meow as mw
import gdsfactory as gf

gf.gpdk.PDK.activate()

# gdsfactory cross section for o-band waveguides on GF
xs_so = gf.cross_section.cross_section(
    width=0.35,
    offset=0,
    layer=(2,0)
)

# gdsfactory cross section for c-band waveguide on .22 Si
xs_sc = gf.cross_section.cross_section(
    width=0.45,
    offset=0,
    layer=(2,0)
)

def stack_extrusions_GF(
    wg_height: float = 0.160,
    Sz: float = 2.0,
):
    """create some simple extrusion rules

    Args:
        wg_height: the Si waveguide thickness
    """
    extrusions = {
        (2, 0): [
            mw.GdsExtrusionRule(
                material=mw.silicon_oxide,
                h_min=-Sz/2,
                h_max=-wg_height/2,
                mesh_order=2,
            ),
            mw.GdsExtrusionRule(
                material=mw.silicon,
                h_min=-wg_height/2,
                h_max=wg_height/2,
                mesh_order=1,
            ),
            mw.GdsExtrusionRule(
                material=mw.silicon_oxide,
                h_min=wg_height/2,
                h_max=Sz/2,
                mesh_order=2,
            ),
        ],
        (1, 0): [
            mw.GdsExtrusionRule(
                material=mw.silicon_oxide,
                h_min=-Sz/2,
                h_max=Sz/2,
                mesh_order=2,
            ),
        ],
    }
    return extrusions

def create_MMI_1x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset):
    mmi = gf.Component('mmi')
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

def load_1x2_MMI(mm_length,straight_length,waveguide_width,waveguide_offset,oxide_width,x_offset,taper_length,mm_width,taper_width):
    mmi_length = mm_length + 2*taper_length
    mmi_top = create_MMI_1x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset/2)
    mmi_top.move((x_offset+mmi_length/2 + straight_length, 0))
    # Add input and output waveguide for measuring transmission to the top waveguide
    input_waveguide = mmi_top << gf.components.straight(length=straight_length,width=waveguide_width,cross_section=xs_so)
    input_waveguide.move((x_offset, 0))
    output_waveguide = mmi_top << gf.components.straight(length=straight_length,width=waveguide_width,cross_section=xs_so)
    output_waveguide.move((x_offset+mmi_length+straight_length, waveguide_offset/2))
    # Add layer (1,0) for oxide cladding
    mmi_top.add_polygon([(x_offset,-oxide_width/2),(x_offset+mmi_length+2*straight_length,-oxide_width/2),(x_offset+mmi_length+2*straight_length,oxide_width/2),(x_offset,oxide_width/2)], layer=(1,0))
    return mmi_top

def create_MMI_2x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset,name):
    mmi = gf.Component(name)
    # multimode section
    mmi.add_polygon([(-mm_length/2,-mm_width/2),
                     (mm_length/2,-mm_width/2),
                     (mm_length/2,mm_width/2),
                     (-mm_length/2,mm_width/2)], 
                     layer=(2,0))
    # input tapers
    mmi.add_polygon([(-mm_length/2-taper_length,-waveguide_width/2+waveguide_offset),
                     (-mm_length/2,-taper_width/2+waveguide_offset),
                     (-mm_length/2,taper_width/2+waveguide_offset),
                     (-mm_length/2-taper_length,waveguide_width/2+waveguide_offset)], 
                     layer=(2,0))
    mmi.add_polygon([(-mm_length/2-taper_length,-waveguide_width/2-waveguide_offset),
                     (-mm_length/2,-taper_width/2-waveguide_offset),
                     (-mm_length/2,taper_width/2-waveguide_offset),
                     (-mm_length/2-taper_length,waveguide_width/2-waveguide_offset)], 
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

def load_2x2_MMI_top(mm_length,straight_length,waveguide_width,waveguide_offset,oxide_width,x_offset,taper_length,mm_width,taper_width):
    mmi_length = mm_length + 2*taper_length
    mmi_top = create_MMI_2x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset/2,"mmi_top")
    mmi_top.move((x_offset+mmi_length/2 + straight_length, 0))
    # Add input and output waveguide for measuring transmission to the top waveguide
    input_waveguide = mmi_top << gf.components.straight(length=straight_length,width=waveguide_width,cross_section=xs_sc)
    input_waveguide.move((x_offset, waveguide_offset/2))
    output_waveguide = mmi_top << gf.components.straight(length=straight_length,width=waveguide_width,cross_section=xs_sc)
    output_waveguide.move((x_offset+mmi_length+straight_length, waveguide_offset/2))
    # Add layer (1,0) for oxide cladding
    mmi_top.add_polygon([(x_offset,-oxide_width/2),(x_offset+mmi_length+2*straight_length,-oxide_width/2),(x_offset+mmi_length+2*straight_length,oxide_width/2),(x_offset,oxide_width/2)], layer=(1,0))
    return mmi_top

def create_MMI_1x4(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset, name):
    mmi = gf.Component(name)
    mmi.add_polygon([(-mm_length/2,-mm_width/2),
                     (mm_length/2,-mm_width/2),
                     (mm_length/2,mm_width/2),
                     (-mm_length/2,mm_width/2)],
                     layer=(2,0))
    mmi.add_polygon([(-mm_length/2-taper_length,-waveguide_width/2),
                     (-mm_length/2,-taper_width/2),
                     (-mm_length/2,taper_width/2),
                     (-mm_length/2-taper_length,waveguide_width/2)],
                     layer=(2,0))
    for yo in [3*waveguide_offset/2, waveguide_offset/2, -waveguide_offset/2, -3*waveguide_offset/2]:
        mmi.add_polygon([(mm_length/2,-taper_width/2+yo),
                         (mm_length/2+taper_length,-waveguide_width/2+yo),
                         (mm_length/2+taper_length,waveguide_width/2+yo),
                         (mm_length/2,taper_width/2+yo)],
                         layer=(2,0))
    return mmi

def load_1x4_MMI(mm_length, straight_length, waveguide_width, waveguide_offset, oxide_width, x_offset, taper_length, mm_width, taper_width, output_idx):
    mmi_length = mm_length + 2*taper_length
    output_offsets = [3*waveguide_offset/2, waveguide_offset/2, -waveguide_offset/2, -3*waveguide_offset/2]
    yo = output_offsets[output_idx]
    mmi = create_MMI_1x4(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset, f'mmi_1x4_{output_idx}')
    mmi.move((x_offset+mmi_length/2 + straight_length, 0))
    input_waveguide = mmi << gf.components.straight(length=straight_length, width=waveguide_width, cross_section=xs_so)
    input_waveguide.move((x_offset, 0))
    output_waveguide = mmi << gf.components.straight(length=straight_length, width=waveguide_width, cross_section=xs_so)
    output_waveguide.move((x_offset+mmi_length+straight_length, yo))
    mmi.add_polygon([(x_offset,-oxide_width/2),(x_offset+mmi_length+2*straight_length,-oxide_width/2),(x_offset+mmi_length+2*straight_length,oxide_width/2),(x_offset,oxide_width/2)], layer=(1,0))
    return mmi

def load_2x2_MMI_bot(mm_length,straight_length,waveguide_width,waveguide_offset,oxide_width,x_offset,taper_length,mm_width,taper_width):
    mmi_length = mm_length + 2*taper_length
    mmi_top = create_MMI_2x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset/2,"mmi_bot")
    mmi_top.move((x_offset+mmi_length/2 + straight_length, 0))
    # Add input and output waveguide for measuring transmission to the top waveguide
    input_waveguide = mmi_top << gf.components.straight(length=straight_length,width=waveguide_width,cross_section=xs_sc)
    input_waveguide.move((x_offset, waveguide_offset/2))
    output_waveguide = mmi_top << gf.components.straight(length=straight_length,width=waveguide_width,cross_section=xs_sc)
    output_waveguide.move((x_offset+mmi_length+straight_length, -waveguide_offset/2))
    # Add layer (1,0) for oxide cladding
    mmi_top.add_polygon([(x_offset,-oxide_width/2),(x_offset+mmi_length+2*straight_length,-oxide_width/2),(x_offset+mmi_length+2*straight_length,oxide_width/2),(x_offset,oxide_width/2)], layer=(1,0))
    return mmi_top