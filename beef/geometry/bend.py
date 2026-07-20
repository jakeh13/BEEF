"""BendGeometry: 90° bend GDS for non-collinear FDTD I/O ports.

Migrated from engine/helpers/meep.py::create_bend_gds_for_meep.
"""
from __future__ import annotations

import struct
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np

_WG_LAYER = (2, 0)


def _arc_polygon(cx, cy, r_inner, r_outer, theta0_deg, theta1_deg, n_pts=64):
    t0, t1 = np.radians(theta0_deg), np.radians(theta1_deg)
    thetas = np.linspace(t0, t1, n_pts)
    outer = [(cx + r_outer*np.cos(t), cy + r_outer*np.sin(t)) for t in thetas]
    inner = [(cx + r_inner*np.cos(t), cy + r_inner*np.sin(t)) for t in reversed(thetas)]
    return outer + inner


def _real8(value: float) -> bytes:
    """Encode a float as GDSII's 8-byte excess-64 base-16 REAL8 format."""
    if value == 0.0:
        return b"\x00" * 8
    sign = 0
    if value < 0:
        sign = 0x80
        value = -value
    exponent = 0
    while value >= 1.0:
        value /= 16.0
        exponent += 1
    while value < 1.0 / 16.0:
        value *= 16.0
        exponent -= 1
    mantissa = int(round(value * (1 << 56)))
    if mantissa >= (1 << 56):
        mantissa >>= 4
        exponent += 1
    byte0 = sign | ((exponent + 64) & 0x7F)
    return bytes([byte0]) + mantissa.to_bytes(7, "big")


def _gds_record(rec_type: int, data_type: int, payload: bytes = b"") -> bytes:
    return struct.pack(">HBB", 4 + len(payload), rec_type, data_type) + payload


def _gds_str(s: str) -> bytes:
    b = s.encode("ascii")
    return b if len(b) % 2 == 0 else b + b"\x00"


def _write_gds(
    path: str,
    polys_by_layer: Dict[Tuple[int, int], List[List[Tuple[float, float]]]],
    user_unit_m: float = 1e-6,
    dbu_m: float = 1e-9,
) -> None:
    """Write a minimal single-structure GDSII file with a BOUNDARY element per
    polygon — a self-contained implementation of the plain, documented GDSII
    binary record format (no external GDS library needed): HEADER / BGNLIB /
    LIBNAME / UNITS / BGNSTR / STRNAME, one BOUNDARY+LAYER+DATATYPE+XY+ENDEL
    per polygon, ENDSTR / ENDLIB.

    user_unit_m/dbu_m: metres per user unit / per database unit (defaults:
    user unit = 1 µm, database unit = 1 nm — the same convention klayout/
    gdsfactory use).
    """
    dbu_per_user_unit = user_unit_m / dbu_m
    ts = struct.pack(">6h", 0, 1, 1, 0, 0, 0) * 2  # placeholder mod/access timestamps

    out = bytearray()
    out += _gds_record(0x00, 0x02, struct.pack(">h", 600))                      # HEADER
    out += _gds_record(0x01, 0x02, ts)                                          # BGNLIB
    out += _gds_record(0x02, 0x06, _gds_str("BEEF_LIB"))                       # LIBNAME
    out += _gds_record(0x03, 0x05, _real8(dbu_m / user_unit_m) + _real8(dbu_m))  # UNITS
    out += _gds_record(0x05, 0x02, ts)                                          # BGNSTR
    out += _gds_record(0x06, 0x06, _gds_str("bends"))                           # STRNAME

    for (layer, datatype), polys in polys_by_layer.items():
        for poly in polys:
            pts = list(poly)
            if pts[0] != pts[-1]:
                pts = pts + [pts[0]]  # GDSII boundaries must be closed
            xy = b"".join(
                struct.pack(">ii", round(x * dbu_per_user_unit), round(y * dbu_per_user_unit))
                for x, y in pts
            )
            out += _gds_record(0x08, 0x00)                                      # BOUNDARY
            out += _gds_record(0x0D, 0x02, struct.pack(">h", layer))            # LAYER
            out += _gds_record(0x0E, 0x02, struct.pack(">h", datatype))         # DATATYPE
            out += _gds_record(0x10, 0x03, xy)                                  # XY
            out += _gds_record(0x11, 0x00)                                      # ENDEL

    out += _gds_record(0x07, 0x00)                                              # ENDSTR
    out += _gds_record(0x04, 0x00)                                              # ENDLIB

    with open(path, "wb") as f:
        f.write(bytes(out))


@dataclass
class BendGeometry:
    """90° bend geometry for attaching non-collinear I/O to a MEEP simulation.

    Args:
        waveguide_width: width of the waveguide (µm).
        waveguide_offset: full center-to-center port spacing (µm).
        bend_radius: bend radius (µm).
        bend_type: "euler" or "circular" — both currently resolve to the same
            constant-curvature arc-polygon construction below (see
            coupler.py's _sbend_polygon for the analogous split, where
            "euler" instead gets a distinct Bezier curve; this codebase has
            never implemented a literal Euler-spiral/clothoid bend). Kept as
            two named, validated options for API stability / to describe
            intent at call sites.
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
        if self.bend_type not in ("euler", "circular"):
            raise ValueError(f"unknown bend_type {self.bend_type!r}, expected 'euler' or 'circular'")

        R = self.bend_radius
        hw = self.waveguide_width / 2
        yo = self.waveguide_offset / 2
        straight_len = self.monitor_margin + dpml

        polys: List[List[Tuple[float, float]]] = []
        x_extra_right = x_extra_left = y_extra = 0.0

        def _add_straight_y(x_center, y0, length, upward=True):
            dy = length if upward else -length
            y_lo, y_hi = min(y0, y0+dy), max(y0, y0+dy)
            polys.append([(x_center-hw, y_lo), (x_center+hw, y_lo),
                          (x_center+hw, y_hi), (x_center-hw, y_hi)])

        def _add_arc(arc_cx, arc_cy, t0, t1):
            polys.append(_arc_polygon(arc_cx, arc_cy, R-hw, R+hw, t0, t1))

        if "output" in self.bend_locations:
            x_stub_end = device_sx / 2 + self.stub_length
            _add_arc(x_stub_end, yo + R, -90, 0)
            _add_arc(x_stub_end, -yo - R, 0, 90)
            _add_straight_y(x_stub_end + R, yo + R, straight_len, upward=True)
            _add_straight_y(x_stub_end + R, -yo - R, straight_len, upward=False)

            x_extra_right = R + self.monitor_margin
            y_extra = max(y_extra, R + self.monitor_margin)

        if "input" in self.bend_locations:
            x_stub_start = -(device_sx / 2 + self.stub_length)
            _add_arc(x_stub_start, yo + R, 180, 270)
            _add_straight_y(x_stub_start - R, yo + R, straight_len, upward=True)
            _add_arc(x_stub_start, -yo - R, 90, 180)
            _add_straight_y(x_stub_start - R, -yo - R, straight_len, upward=False)

            x_extra_left = R + self.monitor_margin
            y_extra = max(y_extra, R + self.monitor_margin)

        _write_gds(output_path, {_WG_LAYER: polys})
        return {"x_extra_right": x_extra_right, "x_extra_left": x_extra_left, "y_extra": y_extra}
