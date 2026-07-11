"""BendGeometry: 90° bend GDS for non-collinear FDTD I/O ports.

Migrated from engine/helpers/meep.py::create_bend_gds_for_meep.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np
import gdsfactory as gf

_WG_LAYER = (2, 0)


def _arc_polygon(cx, cy, r_inner, r_outer, theta0_deg, theta1_deg, n_pts=64):
    t0, t1 = np.radians(theta0_deg), np.radians(theta1_deg)
    thetas = np.linspace(t0, t1, n_pts)
    outer = [(cx + r_outer*np.cos(t), cy + r_outer*np.sin(t)) for t in thetas]
    inner = [(cx + r_inner*np.cos(t), cy + r_inner*np.sin(t)) for t in reversed(thetas)]
    return outer + inner


@dataclass
class BendGeometry:
    """90° bend geometry for attaching non-collinear I/O to a MEEP simulation.

    Args:
        waveguide_width: width of the waveguide (µm).
        waveguide_offset: full center-to-center port spacing (µm).
        bend_radius: bend radius (µm).
        bend_type: "euler" (gdsfactory) or "circular" (arc polygon).
        bend_locations: list containing "input", "output", or both.
        stub_length: straight segment between device edge and bend start (µm).
        monitor_margin: sm_monitor_position — how far monitor sits from PML (µm).
    """
    waveguide_width: float = 0.35
    waveguide_offset: float = 1.0
    bend_radius: float = 5.0
    bend_type: str = "euler"
    bend_locations: List[str] = field(default_factory=list)
    stub_length: float = 0.5
    monitor_margin: float = 0.25

    def build(self, device_sx: float, dpml: float, output_path: str) -> Dict[str, float]:
        """Write bend GDS and return domain expansion dict.

        Returns:
            {"x_extra_right", "x_extra_left", "y_extra"}: domain extensions (µm).
        """
        R = self.bend_radius
        hw = self.waveguide_width / 2
        yo = self.waveguide_offset / 2
        straight_len = self.monitor_margin + dpml
        xs = gf.cross_section.cross_section(width=self.waveguide_width, layer=(2, 0))

        c = gf.Component("bends")
        x_extra_right = x_extra_left = y_extra = 0.0

        def _add_straight_x(x0, y_center, length):
            c.add_polygon([(x0, y_center-hw), (x0+length, y_center-hw),
                           (x0+length, y_center+hw), (x0, y_center+hw)], layer=_WG_LAYER)

        def _add_straight_y(x_center, y0, length, upward=True):
            dy = length if upward else -length
            y_lo, y_hi = min(y0, y0+dy), max(y0, y0+dy)
            c.add_polygon([(x_center-hw, y_lo), (x_center+hw, y_lo),
                           (x_center+hw, y_hi), (x_center-hw, y_hi)], layer=_WG_LAYER)

        def _add_arc(arc_cx, arc_cy, t0, t1):
            c.add_polygon(_arc_polygon(arc_cx, arc_cy, R-hw, R+hw, t0, t1), layer=_WG_LAYER)

        if "output" in self.bend_locations:
            x_stub_end = device_sx / 2 + self.stub_length

            if self.bend_type == "euler":
                b_ot = gf.components.bend_euler(radius=R, angle=-90, cross_section=xs)
                o2_ot = np.array(b_ot.ports["o2"].center)
                for pts in [np.array(p) for ps in b_ot.get_polygons_points(by="tuple").values() for p in ps]:
                    c.add_polygon(np.column_stack([pts[:, 0] + x_stub_end,
                                                   -pts[:, 1] + yo]), layer=_WG_LAYER)
                _add_straight_y(o2_ot[0] + x_stub_end, -o2_ot[1] + yo, straight_len, upward=True)

                b_ob = gf.components.bend_euler(radius=R, angle=-90, cross_section=xs)
                ref_ob = c.add_ref(b_ob)
                ref_ob.move((x_stub_end, -yo))
                o2_ob = np.array(b_ob.ports["o2"].center)
                _add_straight_y(x_stub_end + o2_ob[0], -yo + o2_ob[1], straight_len, upward=False)
            else:
                _add_arc(x_stub_end, yo + R, -90, 0)
                _add_arc(x_stub_end, -yo - R, 0, 90)
                _add_straight_y(x_stub_end + R, yo + R, straight_len, upward=True)
                _add_straight_y(x_stub_end + R, -yo - R, straight_len, upward=False)

            x_extra_right = R + self.monitor_margin
            y_extra = max(y_extra, R + self.monitor_margin)

        if "input" in self.bend_locations:
            x_stub_start = -(device_sx / 2 + self.stub_length)

            if self.bend_type == "euler":
                b_it = gf.components.bend_euler(radius=R, angle=-90, cross_section=xs)
                o1_it = np.array(b_it.ports["o1"].center)
                o2_it = np.array(b_it.ports["o2"].center)
                o2_it_t = np.array([-o2_it[1], -o2_it[0]])
                tx_it = np.array([x_stub_start, yo]) - o2_it_t
                for pts in [np.array(p) for ps in b_it.get_polygons_points(by="tuple").values() for p in ps]:
                    c.add_polygon(np.column_stack([-pts[:, 1] + tx_it[0],
                                                   -pts[:, 0] + tx_it[1]]), layer=_WG_LAYER)
                o1_it_w = np.array([-o1_it[1], -o1_it[0]]) + tx_it
                _add_straight_y(o1_it_w[0], o1_it_w[1], straight_len, upward=True)

                b_ib = gf.components.bend_euler(radius=R, angle=-90, cross_section=xs)
                ref_ib = c.add_ref(b_ib)
                ref_ib.rotate(+90, center=(0, 0))
                o2_bot = np.array(b_ib.ports["o2"].center)
                o2_bot_rotated = np.array([-o2_bot[1], o2_bot[0]])
                ref_ib.move(-o2_bot_rotated + np.array([x_stub_start, -yo]))
                o1_bot_world = np.array(ref_ib.ports["o1"].center)
                _add_straight_y(o1_bot_world[0], o1_bot_world[1], straight_len, upward=False)
            else:
                _add_arc(x_stub_start, yo + R, 180, 270)
                _add_straight_y(x_stub_start - R, yo + R, straight_len, upward=True)
                _add_arc(x_stub_start, -yo - R, 90, 180)
                _add_straight_y(x_stub_start - R, -yo - R, straight_len, upward=False)

            x_extra_left = R + self.monitor_margin
            y_extra = max(y_extra, R + self.monitor_margin)

        c.write_gds(output_path)
        return {"x_extra_right": x_extra_right, "x_extra_left": x_extra_left, "y_extra": y_extra}
