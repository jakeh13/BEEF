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

# Function to generate temp gds file for MMI geometry (to shift the MMI to the correct simulation coordinates)
def create_temp_gds_for_mmi_2x2_input(mm_length, taper_length, Sx, dpml, waveguide_length, waveguide_width, mm_width, taper_width, waveguide_offset,name):
    mmi_total_length = mm_length + 2*taper_length
    temp = create_MMI_2x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset/2,name)
    temp.move((mmi_total_length/2 - Sx/2 + dpml + waveguide_length, 0))
    # save as temp.gds
    temp.write_gds(gdsdir='gds')

# Function to generate temp gds file for MMI geometry (to shift the MMI to the correct simulation coordinates)
def create_temp_gds_for_mmi_2x2_output(mm_length, taper_length, Sx, dpml, waveguide_length, waveguide_width, mm_width, taper_width, waveguide_offset):
    mmi_total_length = mm_length + 2*taper_length
    temp = create_MMI_2x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset/2,"temp_out")
    temp.move((-mmi_total_length/2 + Sx/2 - dpml - waveguide_length, 0))
    # save as temp.gds
    temp.write_gds(gdsdir='gds')

# Function to generate temp gds file for MMI geometry (to shift the MMI to the correct simulation coordinates)
def create_temp_gds_for_full_mmi_2x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset):
    temp = create_MMI_2x2(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset/2,"temp_full")
    # save as temp.gds
    temp.write_gds(gdsdir='gds')

def prepare_gds_for_meep(gds_path, output_path):
    component = gf.import_gds(gds_path)
    cx, cy = component.center
    component.move((-cx, -cy))
    component.write_gds(output_path)
    bbox = component.bbox()
    return float(bbox.right - bbox.left)

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

def create_temp_gds_for_mmi_1x4_input(mm_length, taper_length, Sx, dpml, waveguide_length, waveguide_width, mm_width, taper_width, waveguide_offset):
    mmi_total_length = mm_length + 2*taper_length
    temp = create_MMI_1x4(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset, "temp")
    temp.move((mmi_total_length/2 - Sx/2 + dpml + waveguide_length, 0))
    temp.write_gds(gdsdir='gds')

def create_temp_gds_for_mmi_1x4_output(mm_length, taper_length, Sx, dpml, waveguide_length, waveguide_width, mm_width, taper_width, waveguide_offset):
    mmi_total_length = mm_length + 2*taper_length
    temp = create_MMI_1x4(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset, "temp_out")
    temp.move((-mmi_total_length/2 + Sx/2 - dpml - waveguide_length, 0))
    temp.write_gds(gdsdir='gds')

def create_temp_gds_for_full_mmi_1x4(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset):
    temp = create_MMI_1x4(mm_length, taper_length, waveguide_width, mm_width, taper_width, waveguide_offset, "temp_full")
    temp.write_gds(gdsdir='gds')

# --- MMI MUX (1x3 → 3x3 cascade) GDS helpers ---

def create_MMI_MUX(geo, name="mmi_mux"):
    """Full MUX device polygon geometry, centered at x=0.

    Device x-span: 4*taper_length + mm_1_length + between_wg_length + mm_2_length
    Ports (y): top=+waveguide_offset, mid=0, bot=-waveguide_offset
    """
    c = gf.Component(name)
    L1 = geo.mm_1_length
    L2 = geo.mm_2_length
    Lt = geo.taper_length
    Lc = geo.between_waveguide_length
    mw = geo.mm_width
    wo = geo.waveguide_offset

    total = 4*Lt + L1 + Lc + L2
    x0 = -total / 2  # left edge of input taper

    # --- input taper (mm_waveguide_width → input_taper_width) ---
    iw = geo.mm_waveguide_width
    itw = geo.input_taper_width
    c.add_polygon([
        (x0,              -iw/2),
        (x0 + Lt,         -itw/2),
        (x0 + Lt,          itw/2),
        (x0,               iw/2),
    ], layer=(2, 0))

    # --- MMI 1 body ---
    x1 = x0 + Lt
    c.add_polygon([
        (x1,       -mw/2),
        (x1 + L1,  -mw/2),
        (x1 + L1,   mw/2),
        (x1,        mw/2),
    ], layer=(2, 0))

    # --- output tapers of MMI 1 (3 ports) ---
    # taper from MMI-side width → WG-side width
    x2 = x1 + L1  # right face of MMI 1
    # top port: top_bot_taper_width → top_waveguide_width
    tbw = geo.top_bot_taper_width
    midw = geo.middle_taper_width
    tww = geo.top_waveguide_width
    ww = geo.waveguide_width
    for (yo, mmi_w, wg_w) in [
        ( wo, tbw,  tww),   # top
        (  0, midw, ww),    # mid
        (-wo, tbw,  ww),    # bot
    ]:
        c.add_polygon([
            (x2,        -mmi_w/2 + yo),
            (x2 + Lt,   -wg_w/2  + yo),
            (x2 + Lt,    wg_w/2  + yo),
            (x2,         mmi_w/2 + yo),
        ], layer=(2, 0))

    # --- connecting waveguides ---
    x3 = x2 + Lt
    for (yo, w) in [(wo, tww), (0, ww), (-wo, ww)]:
        c.add_polygon([
            (x3,       -w/2 + yo),
            (x3 + Lc,  -w/2 + yo),
            (x3 + Lc,   w/2 + yo),
            (x3,        w/2 + yo),
        ], layer=(2, 0))

    # --- input tapers of MMI 2 (WG-side → MMI-side, reversed) ---
    x4 = x3 + Lc
    for (yo, wg_w, mmi_w) in [
        ( wo, tww,  tbw),
        (  0, ww,   midw),
        (-wo, ww,   tbw),
    ]:
        c.add_polygon([
            (x4,        -wg_w/2  + yo),
            (x4 + Lt,   -mmi_w/2 + yo),
            (x4 + Lt,    mmi_w/2 + yo),
            (x4,         wg_w/2  + yo),
        ], layer=(2, 0))

    # --- MMI 2 body ---
    x5 = x4 + Lt
    c.add_polygon([
        (x5,       -mw/2),
        (x5 + L2,  -mw/2),
        (x5 + L2,   mw/2),
        (x5,        mw/2),
    ], layer=(2, 0))

    # --- output tapers of MMI 2 (same as MMI 1 output tapers) ---
    x6 = x5 + L2
    for (yo, mmi_w, wg_w) in [
        ( wo, tbw,  tww),
        (  0, midw, ww),
        (-wo, tbw,  ww),
    ]:
        c.add_polygon([
            (x6,        -mmi_w/2 + yo),
            (x6 + Lt,   -wg_w/2  + yo),
            (x6 + Lt,    wg_w/2  + yo),
            (x6,         mmi_w/2 + yo),
        ], layer=(2, 0))

    return c


def create_temp_gds_for_full_mmi_mux(geo):
    """Full device GDS for FDTD simulation, saved as gds/temp_mux_full.gds."""
    c = create_MMI_MUX(geo, name="temp_mux_full")
    c.write_gds(gdsdir='gds')


def create_temp_gds_for_mmi_mux_input(geo, sim, Sx, name="temp_mux_in"):
    """Input-side GDS for FDTD input simulation.

    Contains: input taper + first sim.mm_sim_length of MMI 1 body.
    Positioned so taper left edge = x = -Sx/2 + dpml + waveguide_length.
    """
    Lt = geo.taper_length
    iw = geo.mm_waveguide_width
    itw = geo.input_taper_width
    mw = geo.mm_width
    mm_sim = sim.mm_sim_length
    dpml = sim.dpml
    wl = sim.waveguide_length

    c = gf.Component(name)
    # x=0 is left edge of taper in local coords; shift to MEEP coords below
    # input taper
    c.add_polygon([
        (0,      -iw/2),
        (Lt,     -itw/2),
        (Lt,      itw/2),
        (0,       iw/2),
    ], layer=(2, 0))
    # partial MMI 1 body
    c.add_polygon([
        (Lt,          -mw/2),
        (Lt + mm_sim, -mw/2),
        (Lt + mm_sim,  mw/2),
        (Lt,           mw/2),
    ], layer=(2, 0))

    # shift: left edge of taper → -Sx/2 + dpml + wl
    x_shift = (-Sx/2 + dpml + wl)
    c.move((x_shift, 0))
    c.write_gds(gdsdir='gds')


def create_temp_gds_for_mmi_mux_output(geo, sim, Sx, port_idx, name=None):
    """Output-side GDS for FDTD output simulation (by reciprocity).

    Contains: last sim.mm_sim_length of MMI 2 body + output taper for port_idx.
    port_idx: 0=top (y=+wo), 1=mid (y=0), 2=bot (y=-wo)
    Positioned so taper right edge = x = Sx/2 - dpml - waveguide_length.
    """
    if name is None:
        name = f"temp_mux_out_{port_idx}"
    Lt = geo.taper_length
    mw = geo.mm_width
    wo = geo.waveguide_offset
    mm_sim = sim.mm_sim_length
    dpml = sim.dpml
    wl = sim.waveguide_length

    port_params = [
        ( wo, geo.top_bot_taper_width,  geo.top_waveguide_width),
        (  0, geo.middle_taper_width,   geo.waveguide_width),
        (-wo, geo.top_bot_taper_width,  geo.waveguide_width),
    ]
    yo, mmi_w, wg_w = port_params[port_idx]

    c = gf.Component(name)
    # local coords: x=0 at right edge of MMI 2 body (= left edge of output taper)
    # partial MMI 2 body (extends left)
    c.add_polygon([
        (-mm_sim, -mw/2),
        (0,       -mw/2),
        (0,        mw/2),
        (-mm_sim,  mw/2),
    ], layer=(2, 0))
    # output taper
    c.add_polygon([
        (0,      -mmi_w/2 + yo),
        (Lt,     -wg_w/2  + yo),
        (Lt,      wg_w/2  + yo),
        (0,       mmi_w/2 + yo),
    ], layer=(2, 0))

    # shift: right edge of taper → Sx/2 - dpml - wl
    # right edge of taper in local coords = Lt
    x_shift = (Sx/2 - dpml - wl) - Lt
    c.move((x_shift, 0))
    c.write_gds(gdsdir='gds')