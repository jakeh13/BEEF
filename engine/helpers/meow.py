import meow as mw
import gdsfactory as gf

def reset_gf():
    gf.clear_cache()
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
    return xs_so, xs_sc

xs_so, xs_sc = reset_gf()

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
    xs_so, xs_sc = reset_gf()
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
    xs_so, xs_sc = reset_gf()
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
    xs_so, xs_sc = reset_gf()
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

# --- MMI MUX EME geometry helpers ---

def _add_cascade_polygons(c, geo, x_start=0.0):
    """Add cascade polygons to c starting at x=x_start.

    Sequence: MMI1_body | 3×out_tapers | 3×conn_wgs | 3×in_tapers2 | MMI2_body
    (no input taper, no output tapers — used for HYBRID mm_simulation GDS)
    """
    L1 = geo.mm_1_length
    L2 = geo.mm_2_length
    Lt = geo.taper_length
    Lc = geo.between_waveguide_length
    mw_val = geo.mm_width
    wo = geo.waveguide_offset
    tbw = geo.top_bot_taper_width
    midw = geo.middle_taper_width
    tww = geo.top_waveguide_width
    ww = geo.waveguide_width

    x = x_start
    c.add_polygon([(x,-mw_val/2),(x+L1,-mw_val/2),(x+L1,mw_val/2),(x,mw_val/2)], layer=(2,0))
    x += L1
    for (yo, mmi_w, wg_w) in [(wo,tbw,tww),(0,midw,ww),(-wo,tbw,ww)]:
        c.add_polygon([(x,-mmi_w/2+yo),(x+Lt,-wg_w/2+yo),(x+Lt,wg_w/2+yo),(x,mmi_w/2+yo)], layer=(2,0))
    x += Lt
    for (yo, w) in [(wo,tww),(0,ww),(-wo,ww)]:
        c.add_polygon([(x,-w/2+yo),(x+Lc,-w/2+yo),(x+Lc,w/2+yo),(x,w/2+yo)], layer=(2,0))
    x += Lc
    for (yo, wg_w, mmi_w) in [(wo,tww,tbw),(0,ww,midw),(-wo,ww,tbw)]:
        c.add_polygon([(x,-wg_w/2+yo),(x+Lt,-mmi_w/2+yo),(x+Lt,mmi_w/2+yo),(x,wg_w/2+yo)], layer=(2,0))
    x += Lt
    c.add_polygon([(x,-mw_val/2),(x+L2,-mw_val/2),(x+L2,mw_val/2),(x,mw_val/2)], layer=(2,0))


def _add_mux_polygons(c, geo, x_start=0.0):
    """Add all MUX device polygons to component c starting at x=x_start.

    Polygon sequence (left→right):
      input_taper | MMI1_body | 3×out_tapers | 3×conn_wgs | 3×in_tapers2 | MMI2_body | 3×out_tapers2
    """
    L1 = geo.mm_1_length
    L2 = geo.mm_2_length
    Lt = geo.taper_length
    Lc = geo.between_waveguide_length
    mw = geo.mm_width
    wo = geo.waveguide_offset
    iw = geo.mm_waveguide_width
    itw = geo.input_taper_width
    tbw = geo.top_bot_taper_width
    midw = geo.middle_taper_width
    tww = geo.top_waveguide_width
    ww = geo.waveguide_width

    x = x_start
    # input taper
    c.add_polygon([(x,-iw/2),(x+Lt,-itw/2),(x+Lt,itw/2),(x,iw/2)], layer=(2,0))
    x += Lt
    # MMI 1 body
    c.add_polygon([(x,-mw/2),(x+L1,-mw/2),(x+L1,mw/2),(x,mw/2)], layer=(2,0))
    x += L1
    # output tapers of MMI 1
    for (yo, mmi_w, wg_w) in [(wo,tbw,tww),(0,midw,ww),(-wo,tbw,ww)]:
        c.add_polygon([(x,-mmi_w/2+yo),(x+Lt,-wg_w/2+yo),(x+Lt,wg_w/2+yo),(x,mmi_w/2+yo)], layer=(2,0))
    x += Lt
    # connecting waveguides
    for (yo, w) in [(wo,tww),(0,ww),(-wo,ww)]:
        c.add_polygon([(x,-w/2+yo),(x+Lc,-w/2+yo),(x+Lc,w/2+yo),(x,w/2+yo)], layer=(2,0))
    x += Lc
    # input tapers of MMI 2
    for (yo, wg_w, mmi_w) in [(wo,tww,tbw),(0,ww,midw),(-wo,ww,tbw)]:
        c.add_polygon([(x,-wg_w/2+yo),(x+Lt,-mmi_w/2+yo),(x+Lt,mmi_w/2+yo),(x,wg_w/2+yo)], layer=(2,0))
    x += Lt
    # MMI 2 body
    c.add_polygon([(x,-mw/2),(x+L2,-mw/2),(x+L2,mw/2),(x,mw/2)], layer=(2,0))
    x += L2
    # output tapers of MMI 2
    for (yo, mmi_w, wg_w) in [(wo,tbw,tww),(0,midw,ww),(-wo,tbw,ww)]:
        c.add_polygon([(x,-mmi_w/2+yo),(x+Lt,-wg_w/2+yo),(x+Lt,wg_w/2+yo),(x,mmi_w/2+yo)], layer=(2,0))


def load_MMI_MUX_full(hp):
    """Full MUX device for EME simulation.

    Layout (left→right):
      input_MM_wg | input_taper | MMI1 | tapers | conn_wgs | tapers | MMI2 | tapers | 3×output_wgs
    Oxide cladding in layer (1,0).
    """
    xs_so, xs_sc = reset_gf()
    geo = hp.geo
    sim = hp.sim
    sl = sim.waveguide_length
    ow = sim.oxide_width

    device_len = 4*geo.taper_length + geo.mm_1_length + geo.between_waveguide_length + geo.mm_2_length
    total_len = sl + device_len + sl

    c = gf.Component('mmi_mux_full')
    _add_mux_polygons(c, geo, x_start=sl)

    # input MM waveguide
    inp = c << gf.components.straight(length=sl, width=geo.mm_waveguide_width,
                                      cross_section=gf.cross_section.cross_section(width=geo.mm_waveguide_width, offset=0, layer=(2,0)))
    inp.move((0, 0))

    # 3 output waveguides
    for (yo, w) in [(geo.waveguide_offset, geo.top_waveguide_width),
                    (0, geo.waveguide_width),
                    (-geo.waveguide_offset, geo.waveguide_width)]:
        xs = gf.cross_section.cross_section(width=w, offset=0, layer=(2,0))
        out_wg = c << gf.components.straight(length=sl, width=w, cross_section=xs)
        out_wg.move((sl + device_len, yo))

    # oxide cladding
    c.add_polygon([(0,-ow/2),(total_len,-ow/2),(total_len,ow/2),(0,ow/2)], layer=(1,0))
    return c


def load_MMI_MUX_cascade(hp):
    """MUX cascade for HYBRID mm_simulation.

    Geometry: MMI1_body | out_tapers | conn_wgs | in_tapers | MMI2_body
    (no input taper, no output tapers).

    x=0 aligns with mm_monitor_position inside MMI 1 (from left face).
    Cells therefore span:
      (mm_1_length - mp) + taper_length + between_wg_length + taper_length + (mm_2_length - mp)
    Output reference plane: mm_monitor_position inside MMI 2 from right face.
    """
    xs_so, xs_sc = reset_gf()
    geo = hp.geo
    sim = hp.sim
    mp = sim.mm_monitor_position
    ow = sim.oxide_width

    # MMI 1 body left face at x=-mp → x=0 is mp inside MMI 1
    x_start = -mp
    total_gds_len = geo.mm_1_length + 2*geo.taper_length + geo.between_waveguide_length + geo.mm_2_length

    c = gf.Component('mmi_mux_cascade')
    _add_cascade_polygons(c, geo, x_start=x_start)
    c.add_polygon([(x_start, -ow/2), (x_start + total_gds_len, -ow/2),
                   (x_start + total_gds_len,  ow/2), (x_start,  ow/2)], layer=(1,0))
    return c


def load_MMI_MUX_port(hp, port_idx):
    """Full MUX device ending at a single output port, for per-port EME S-matrix.

    port_idx: 0=top (y=+wo, w=top_waveguide_width),
              1=mid (y=0,   w=waveguide_width),
              2=bot (y=-wo, w=waveguide_width)
    Input: MM WG at y=0. Output: single WG for the chosen port.
    """
    xs_so, xs_sc = reset_gf()
    geo = hp.geo
    sim = hp.sim
    sl = sim.waveguide_length
    ow = sim.oxide_width
    device_len = 4*geo.taper_length + geo.mm_1_length + geo.between_waveguide_length + geo.mm_2_length
    total_len = sl + device_len + sl

    port_y = [geo.waveguide_offset, 0.0, -geo.waveguide_offset][port_idx]
    port_w = [geo.top_waveguide_width, geo.waveguide_width, geo.waveguide_width][port_idx]

    c = gf.Component(f'mmi_mux_port{port_idx}')
    _add_mux_polygons(c, geo, x_start=sl)

    inp_xs = gf.cross_section.cross_section(width=geo.mm_waveguide_width, offset=0, layer=(2,0))
    inp = c << gf.components.straight(length=sl, width=geo.mm_waveguide_width, cross_section=inp_xs)
    inp.move((0, 0))

    out_xs = gf.cross_section.cross_section(width=port_w, offset=0, layer=(2,0))
    out_wg = c << gf.components.straight(length=sl, width=port_w, cross_section=out_xs)
    out_wg.move((sl + device_len, port_y))

    c.add_polygon([(0,-ow/2),(total_len,-ow/2),(total_len,ow/2),(0,ow/2)], layer=(1,0))
    return c


def load_2x2_MMI_bot(mm_length,straight_length,waveguide_width,waveguide_offset,oxide_width,x_offset,taper_length,mm_width,taper_width):
    xs_so, xs_sc = reset_gf()
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