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


def _arc_polygon(cx, cy, r_inner, r_outer, theta0_deg, theta1_deg, n_pts=64):
    """Return polygon vertices for an arc ring (annulus sector) in layer (2,0)."""
    t0 = np.radians(theta0_deg)
    t1 = np.radians(theta1_deg)
    thetas = np.linspace(t0, t1, n_pts)
    outer = [(cx + r_outer*np.cos(t), cy + r_outer*np.sin(t)) for t in thetas]
    inner = [(cx + r_inner*np.cos(t), cy + r_inner*np.sin(t)) for t in reversed(thetas)]
    return outer + inner


def create_bend_gds_for_meep(
    waveguide_width, waveguide_offset,
    bend_radius, bend_type, bend_locations,
    device_Sx, stub_length, bend_wg_length,
    output_path,
    dpml=0.0,
):
    """Generate GDS with 90° bends (and post/pre-bend straight segments) for MEEP.

    Bends placed at ±waveguide_offset/2 on input and/or output side.
    bend_type: "euler" uses gdsfactory bend_euler; "circular" uses hand-built arc polygon
               (euler uses gf path, circular uses explicit arc ring polygon).
    Straight segments extend bend_wg_length + dpml so they reach through the PML.
    Returns dict with domain expansion: x_extra_left, x_extra_right, y_extra.
    """
    R = bend_radius
    hw = waveguide_width / 2
    yo = waveguide_offset / 2
    straight_len = bend_wg_length + dpml

    c = gf.Component("bends")
    xs = gf.cross_section.cross_section(width=waveguide_width, layer=(2, 0))

    x_extra_right = 0.0
    x_extra_left  = 0.0
    y_extra        = 0.0

    def _add_straight_x(component, x0, y_center, length):
        """Horizontal straight waveguide from x0 to x0+length."""
        component.add_polygon([
            (x0,          y_center - hw),
            (x0 + length, y_center - hw),
            (x0 + length, y_center + hw),
            (x0,          y_center + hw),
        ], layer=(2, 0))

    def _add_straight_y(component, x_center, y0, length, upward=True):
        """Vertical straight waveguide from y0 upward (or downward)."""
        dy = length if upward else -length
        y1 = y0 + dy
        y_lo, y_hi = min(y0, y1), max(y0, y1)
        component.add_polygon([
            (x_center - hw, y_lo),
            (x_center + hw, y_lo),
            (x_center + hw, y_hi),
            (x_center - hw, y_hi),
        ], layer=(2, 0))

    def _add_bend_ccw(component, arc_cx, arc_cy, R, hw, theta0_deg, theta1_deg):
        """Counter-clockwise 90° bend arc polygon."""
        verts = _arc_polygon(arc_cx, arc_cy, R - hw, R + hw, theta0_deg, theta1_deg)
        component.add_polygon(verts, layer=(2, 0))

    def _add_bend_euler_gf(component, x_start, y_start, angle_deg):
        """Euler bend via gdsfactory. angle_deg=+90 turns CCW, -90 turns CW.
        Input port is at (x_start, y_start) traveling in +x direction.
        """
        b = gf.components.bend_euler(radius=R, angle=angle_deg, cross_section=xs)
        ref = component.add_ref(b)
        # gf bend_euler: o1 port at (0,0) facing +x (orientation=180 means facing left,
        # so the waveguide enters from left going right). Move o1 to target position.
        o1 = b.ports["o1"]
        ref.move(-np.array(o1.center) + np.array([x_start, y_start]))
        return b

    if "output" in bend_locations:
        x_stub_end = device_Sx / 2 + stub_length  # x where stub ends (bend starts)

        if bend_type == "euler":
            # Top output: CCW (+x → +y). x-mirror of CW bend: (px,py)→(px,-py).
            # Bypass ComponentReference transforms; manipulate polygon vertices directly.
            b_ot = gf.components.bend_euler(radius=R, angle=-90, cross_section=xs)
            o2_ot = np.array(b_ot.ports["o2"].center)
            for pts in [np.array(p) for ps in b_ot.get_polygons_points(by="tuple").values() for p in ps]:
                c.add_polygon(np.column_stack([pts[:, 0] + x_stub_end,
                                               -pts[:, 1] + yo]), layer=(2, 0))
            _add_straight_y(c, o2_ot[0] + x_stub_end, -o2_ot[1] + yo,
                            straight_len, upward=True)
            # Bottom output: CW (+x → -y, angle=-90).
            b_ob = _add_bend_euler_gf(c, x_stub_end, -yo, -90)
            o2_ob = np.array(b_ob.ports["o2"].center)
            _add_straight_y(c, x_stub_end + o2_ob[0], -yo + o2_ob[1],
                            straight_len, upward=False)
        else:
            # Circular arc: top CCW from 0° to 90°, center at (x_stub_end, yo + R)
            # Wait — for top output: start at (x_stub_end, yo) going +x, center is above
            # Arc center for turning +x→+y: center at (x_stub_end, yo + R)...
            # No: for CCW turn from east to north, center is to the left of travel = above
            # Center = (x_stub_end, yo + R), arc goes from 270° to 0° (CCW)
            # Actually: start at (x_stub_end, yo) which is at angle 270° from center (x_stub_end, yo+R)
            # End at (x_stub_end + R, yo + R) which is at angle 0° from center
            _add_bend_ccw(c, x_stub_end, yo + R, R, hw, -90, 0)   # top: CCW -90°→0°
            _add_bend_ccw(c, x_stub_end, -yo - R, R, hw, 0, 90)   # bot: CW (mirror) 0°→90°
            _add_straight_y(c, x_stub_end + R, yo + R, straight_len, upward=True)
            _add_straight_y(c, x_stub_end + R, -yo - R, straight_len, upward=False)

        x_extra_right = R + bend_wg_length
        y_extra = max(y_extra, R + bend_wg_length)

    if "input" in bend_locations:
        x_stub_start = -(device_Sx / 2 + stub_length)  # x where input stub starts (bend exits here going +x)
        # Top input bend: -y → +x; arc center at (x_stub_start, yo + R)
        # Start of arc (pre-bend wg): (x_stub_start - R, yo + R), traveling -y
        # Arc: CCW from 0° → -90° (i.e., from east going to south),
        # Wait: for -y→+x turn (CW), center is to the right of -y = +x direction
        # center = (x_stub_start, yo + R), start = (x_stub_start - R, yo + R) at 180°, end = (x_stub_start, yo) at 270°
        # That's a CW arc from 180° to 270°.

        if bend_type == "euler":
            # Top input: -y → +x. Combined transform (x-mirror then -90° rotate): (px,py)→(-py,-px).
            # Bypass ComponentReference transforms; manipulate polygon vertices directly.
            b_it = gf.components.bend_euler(radius=R, angle=-90, cross_section=xs)
            o1_it = np.array(b_it.ports["o1"].center)
            o2_it = np.array(b_it.ports["o2"].center)
            o2_it_t = np.array([-o2_it[1], -o2_it[0]])        # o2 after transform
            tx_it = np.array([x_stub_start, yo]) - o2_it_t    # translate so o2→(x_stub_start,yo)
            for pts in [np.array(p) for ps in b_it.get_polygons_points(by="tuple").values() for p in ps]:
                c.add_polygon(np.column_stack([-pts[:, 1] + tx_it[0],
                                               -pts[:, 0] + tx_it[1]]), layer=(2, 0))
            o1_it_w = np.array([-o1_it[1], -o1_it[0]]) + tx_it
            _add_straight_y(c, o1_it_w[0], o1_it_w[1], straight_len, upward=True)

            b_in_bot = gf.components.bend_euler(radius=R, angle=-90, cross_section=xs)
            ref_ib = c.add_ref(b_in_bot)
            ref_ib.rotate(+90, center=(0, 0))
            o2_bot = np.array(b_in_bot.ports["o2"].center)
            o2_bot_rotated = np.array([-o2_bot[1], o2_bot[0]])   # CCW: (-y, x)
            ref_ib.move(-o2_bot_rotated + np.array([x_stub_start, -yo]))
            o1_bot_world = np.array(ref_ib.ports["o1"].center)
            _add_straight_y(c, o1_bot_world[0], o1_bot_world[1], straight_len, upward=False)
        else:
            # Circular: top input CW arc from 180°→270°
            _add_bend_ccw(c, x_stub_start, yo + R, R, hw, 180, 270)  # CW = reversed CCW
            _add_straight_y(c, x_stub_start - R, yo + R, straight_len, upward=True)
            # Bottom input CW arc from 90°→180°
            _add_bend_ccw(c, x_stub_start, -yo - R, R, hw, 90, 180)
            _add_straight_y(c, x_stub_start - R, -yo - R, straight_len, upward=False)

        x_extra_left = R + bend_wg_length
        y_extra = max(y_extra, R + bend_wg_length)

    c.write_gds(output_path)
    return {"x_extra_right": x_extra_right, "x_extra_left": x_extra_left, "y_extra": y_extra}

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