"""Fast, meep/meow-free unit tests for beef.geometry.

These exercise the pure-Python parts of the DeviceGeometry contract
(width_at/center_at/eme_lengths/section_gds's cut validation) without ever
building a real mp.Prism, so they run everywhere pytest does.
"""
from __future__ import annotations

import pytest

from beef.geometry._arc import arc_polygon
from beef.geometry.compound import CompoundGeometry, WaveguideRun
from beef.geometry.coupler import DirectionalCoupler
from beef.geometry.mmi import MMI, MMIPort
from beef.geometry.ring_bend import RingResonatorBendSection
from beef.geometry.ring_coupler import RingResonatorCouplerSection
from beef.geometry.sbend import SBend
from beef.geometry.taper import Taper


# ---------------------------------------------------------------------------
# Taper
# ---------------------------------------------------------------------------

class TestTaper:
    def test_width_at_endpoints_and_midpoint(self):
        t = Taper(length=4.0, width_start=0.35, width_end=0.8)
        assert t.width_at(-2.0) == pytest.approx(0.35)
        assert t.width_at(2.0) == pytest.approx(0.8)
        assert t.width_at(0.0) == pytest.approx((0.35 + 0.8) / 2)

    def test_width_at_clamps_beyond_ends(self):
        t = Taper(length=4.0, width_start=0.35, width_end=0.8)
        assert t.width_at(-10.0) == pytest.approx(0.35)
        assert t.width_at(10.0) == pytest.approx(0.8)

    def test_total_length_and_eme_lengths(self):
        t = Taper(length=4.0, width_start=0.35, width_end=0.8)
        assert t.total_length == 4.0
        assert t.eme_lengths(straight_length=1.5) == [1.5, 4.0, 1.5]

    def test_center_at_is_always_zero(self):
        t = Taper(length=4.0, width_start=0.35, width_end=0.8)
        assert t.center_at(-1.0) == 0.0
        assert t.center_at(1.0) == 0.0

    def test_section_gds_cut_inside_taper_raises(self):
        t = Taper(length=4.0, width_start=0.35, width_end=0.8)
        with pytest.raises(ValueError, match="falls strictly inside a Taper"):
            t.section_gds(0.0, 1.0, oxide_width=3.0)

    def test_section_gds_whole_taper_ok(self):
        t = Taper(length=4.0, width_start=0.35, width_end=0.8)
        d = t.section_gds(-2.0, 2.0, oxide_width=3.0)
        assert (2, 0) in d and (1, 0) in d


# ---------------------------------------------------------------------------
# SBend
# ---------------------------------------------------------------------------

class TestSBend:
    def test_width_at_is_constant(self):
        b = SBend(length=10.0, dy=2.0, width=0.4)
        assert b.width_at(-5.0) == pytest.approx(0.4)
        assert b.width_at(0.0) == pytest.approx(0.4)
        assert b.width_at(5.0) == pytest.approx(0.4)

    def test_center_at_endpoints(self):
        b = SBend(length=10.0, dy=2.0, width=0.4)
        assert b.center_at(-5.0) == pytest.approx(-1.0)
        assert b.center_at(5.0) == pytest.approx(1.0)

    def test_eme_lengths(self):
        b = SBend(length=10.0, dy=2.0, width=0.4)
        assert b.eme_lengths(straight_length=3.0) == [3.0, 10.0, 3.0]

    def test_section_gds_cut_inside_bend_raises(self):
        b = SBend(length=10.0, dy=2.0, width=0.4)
        with pytest.raises(ValueError, match="falls strictly inside an SBend"):
            b.section_gds(0.0, 1.0, oxide_width=3.0)

    def test_unknown_bend_type_raises(self):
        b = SBend(length=10.0, dy=2.0, width=0.4, bend_type="nonsense")
        with pytest.raises(ValueError, match="unknown bend_type"):
            b.polygons()


# ---------------------------------------------------------------------------
# MMI
# ---------------------------------------------------------------------------

class TestMMI:
    def _mmi_1x2(self, **kw):
        return MMI(n_inputs=1, n_outputs=2, mm_length=6.0, mm_width=2.5,
                   waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
                   waveguide_offset=1.0, **kw)

    def test_default_port_placement_symmetric(self):
        m = self._mmi_1x2()
        assert m.input_ys == [0.0]
        assert m.output_ys == pytest.approx([-0.5, 0.5])

    def test_port_count_mismatch_raises(self):
        with pytest.raises(ValueError, match="len\\(inputs\\)"):
            MMI(n_inputs=2, n_outputs=1, mm_length=6.0, mm_width=2.5,
                inputs=[MMIPort(y=0.0, waveguide_width=0.35, taper_length=1.0, taper_width=0.8)])

    def test_total_length(self):
        m = self._mmi_1x2()
        assert m.total_length == pytest.approx(6.0 + 1.0 + 1.0)

    def test_width_at_body_is_mm_width(self):
        m = self._mmi_1x2()
        assert m.width_at(0.0) == pytest.approx(2.5)
        assert m.width_at(-3.0) == pytest.approx(2.5)
        assert m.width_at(3.0) == pytest.approx(2.5)

    def test_width_at_input_taper(self):
        m = self._mmi_1x2()
        # far end of the input taper band: bare access-waveguide width
        assert m.width_at(-4.0, input_idx=0) == pytest.approx(0.35)
        # midpoint of the taper band: linear interpolation
        assert m.width_at(-3.5, input_idx=0) == pytest.approx((0.35 + 0.8) / 2)
        # at/inside the body face (x=-half) the body branch wins: mm_width
        assert m.width_at(-3.0, input_idx=0) == pytest.approx(2.5)

    def test_center_at_selects_requested_port(self):
        m = self._mmi_1x2()
        assert m.center_at(4.0, output_idx=0) == pytest.approx(-0.5)
        assert m.center_at(4.0, output_idx=1) == pytest.approx(0.5)

    def test_eme_lengths_shape(self):
        m = self._mmi_1x2()
        lengths = m.eme_lengths(straight_length=2.0)
        assert lengths == [2.0, 1.0, 6.0, 1.0, 2.0]

    def test_section_gds_cut_inside_taper_band_raises(self):
        m = self._mmi_1x2()
        # input taper band is (-4, -3); cut at -3.5 lands strictly inside it
        with pytest.raises(ValueError, match="falls strictly inside a taper band"):
            m.section_gds(-3.5, -3.0, oxide_width=6.0)

    def test_section_gds_band_boundary_ok(self):
        m = self._mmi_1x2()
        d = m.section_gds(-4.0, -3.0, oxide_width=6.0)
        assert (2, 0) in d

    def test_source_aperture_default_off(self):
        m = self._mmi_1x2()
        assert m.source_aperture_at(0.0) is None

    def test_source_aperture_neighbor_isolation(self):
        m = self._mmi_1x2(use_neighbor_aperture=True)
        # two output ports at y=-0.5/+0.5 (pitch 1.0); at x=4.0 (the taper's
        # far/access-waveguide end) the neighbor is at its bare waveguide_width
        # (0.35): half = 1.0 - 0.35/2 = 0.825 -> full aperture = 1.65
        assert m.source_aperture_at(4.0, output_idx=0) == pytest.approx(1.65)


# ---------------------------------------------------------------------------
# DirectionalCoupler
# ---------------------------------------------------------------------------

class TestDirectionalCoupler:
    def test_total_length(self):
        c = DirectionalCoupler(coupling_length=8.0, bend_length=12.0)
        assert c.total_length == pytest.approx(8.0 + 2 * 12.0)

    def test_eme_lengths(self):
        c = DirectionalCoupler(coupling_length=8.0, bend_length=12.0)
        assert c.eme_lengths(straight_length=5.0) == [5.0, 12.0, 8.0, 12.0, 5.0]


# ---------------------------------------------------------------------------
# CompoundGeometry
# ---------------------------------------------------------------------------

class TestCompoundGeometry:
    def _compound(self):
        seg1 = Taper(length=4.0, width_start=0.35, width_end=0.8)
        seg2 = WaveguideRun(length=3.0, ports=[(0.0, 0.8)])
        seg3 = Taper(length=4.0, width_start=0.8, width_end=0.35)
        return CompoundGeometry([seg1, seg2, seg3])

    def test_requires_at_least_one_segment(self):
        with pytest.raises(ValueError, match="at least 1 segment"):
            CompoundGeometry([])

    def test_total_length_is_sum_of_segments(self):
        g = self._compound()
        assert g.total_length == pytest.approx(4.0 + 3.0 + 4.0)

    def test_width_at_delegates_to_owning_segment(self):
        g = self._compound()
        # left edge of seg1 (its own width_start)
        assert g.width_at(-5.5) == pytest.approx(0.35)
        # inside seg2 (WaveguideRun) -> its single port width
        assert g.width_at(0.0) == pytest.approx(0.8)
        # right edge of seg3 (its own width_end)
        assert g.width_at(5.5) == pytest.approx(0.35)

    def test_eme_lengths(self):
        g = self._compound()
        assert g.eme_lengths(straight_length=1.0) == [1.0, 4.0, 3.0, 4.0, 1.0]

    def test_section_gds_cut_inside_a_taper_segment_raises(self):
        g = self._compound()
        # seg1 spans [-5.5, -1.5]; x0=-4 is strictly inside seg1's taper
        with pytest.raises(ValueError, match="falls strictly inside a Taper"):
            g.section_gds(-4.0, -3.0, oxide_width=3.0)

    def test_section_gds_within_waveguide_run_ok(self):
        g = self._compound()
        d = g.section_gds(-1.0, 1.0, oxide_width=3.0)
        assert (2, 0) in d

    def test_x_outside_bounds_raises(self):
        g = self._compound()
        with pytest.raises(ValueError, match="outside CompoundGeometry"):
            g.width_at(100.0)


# ---------------------------------------------------------------------------
# arc_polygon
# ---------------------------------------------------------------------------

class TestArcPolygon:
    def test_endpoints_and_vertex_count(self):
        poly = arc_polygon(0.0, 0.0, 1.0, 2.0, 0.0, 90.0, n_pts=5)
        assert len(poly) == 2 * 5
        outer_start, outer_end = poly[0], poly[4]
        inner_end, inner_start = poly[5], poly[9]
        assert outer_start == pytest.approx((2.0, 0.0))
        assert outer_end == pytest.approx((0.0, 2.0), abs=1e-9)
        assert inner_start == pytest.approx((1.0, 0.0))
        assert inner_end == pytest.approx((0.0, 1.0), abs=1e-9)


# ---------------------------------------------------------------------------
# RingResonatorCouplerSection
# ---------------------------------------------------------------------------

class TestRingResonatorCouplerSection:
    def _left(self):
        return RingResonatorCouplerSection(
            waveguide_width=0.45, gap=0.25, bend_radius=5.0, near_length=4.0,
            ring_stub_margin=0.5, iface_buffer=1.0, bus_clearance=2.0, side="left",
        )

    def _right(self):
        return RingResonatorCouplerSection(
            waveguide_width=0.45, gap=0.25, bend_radius=5.0, near_length=4.0,
            ring_stub_margin=0.5, iface_buffer=1.0, bus_clearance=2.0, side="right",
        )

    def test_breakpoints_match_reference_constants(self):
        c = self._left()
        assert c.stub_x == pytest.approx(-7.0)
        assert c.hc == pytest.approx(-2.0)
        assert c.iface_x == pytest.approx(-1.0)
        assert c.outer_x == pytest.approx(-9.0)
        assert c.total_length == pytest.approx(8.0)

    def test_ring_port_y_boundary_matches_left_ring_y(self):
        c = self._left()
        assert c.ring_port().y_boundary == pytest.approx(6.20)  # Y_RING+R+S_MON = 0.70+5.0+0.5

    def test_left_and_right_are_x_mirrors(self):
        left, right = self._left(), self._right()
        assert right.stub_x == pytest.approx(-left.stub_x)
        assert right.hc == pytest.approx(-left.hc)
        assert right.iface_x == pytest.approx(-left.iface_x)
        assert right.outer_x == pytest.approx(-left.outer_x)

    def test_width_at_and_center_at_bus_alone_vs_coupled_region(self):
        c = self._left()
        assert c.width_at(-8.0) == pytest.approx(0.45)
        assert c.center_at(-8.0) == pytest.approx(0.0)
        assert c.width_at(-1.5) == pytest.approx(2 * 0.45 + 0.25)
        assert c.center_at(-1.5) == pytest.approx((0.45 + 0.25) / 2)

    def test_cut_inside_bend_raises(self):
        c = self._left()
        with pytest.raises(ValueError, match="falls strictly inside the near-bend arc"):
            c.width_at(-4.0)

    def test_invalid_side_raises(self):
        with pytest.raises(ValueError, match="side must be"):
            RingResonatorCouplerSection(side="up")

    def test_invalid_iface_buffer_raises(self):
        with pytest.raises(ValueError, match="iface_buffer"):
            RingResonatorCouplerSection(near_length=4.0, iface_buffer=3.0)

    def test_section_polygons_full_span(self):
        c = self._left()
        polys = c._section_polygons(c.outer_x, c.iface_x, y_pml=0.5)
        assert len(polys) == 4  # bus, near-bend arc, stub column, sliver

    def test_to_gds_eme_raises(self):
        c = self._left()
        with pytest.raises(NotImplementedError):
            c.to_gds_eme(straight_length=1.0, oxide_width=3.0)


# ---------------------------------------------------------------------------
# RingResonatorBendSection
# ---------------------------------------------------------------------------

class TestRingResonatorBendSection:
    def test_total_length_is_bend_radius(self):
        b = RingResonatorBendSection(bend_radius=5.0)
        assert b.total_length == pytest.approx(5.0)

    def test_out_port(self):
        b = RingResonatorBendSection(bend_radius=5.0, output_margin=0.5)
        port = b.out_port()
        assert port.x == pytest.approx(5.0)
        assert port.y_boundary == pytest.approx(5.5)
        assert port.sign == 1

    def test_cut_inside_bend_raises(self):
        b = RingResonatorBendSection(bend_radius=5.0)
        with pytest.raises(ValueError, match="falls strictly inside the bend arc"):
            b.width_at(2.5)

    def test_section_polygons(self):
        b = RingResonatorBendSection(bend_radius=5.0)
        polys = b._section_polygons(0.0, 5.0, stub_left=1.5, y_pml=0.25)
        assert len(polys) == 3  # input straight, arc, output stub

    def test_to_gds_eme_raises(self):
        b = RingResonatorBendSection()
        with pytest.raises(NotImplementedError):
            b.to_gds_eme(straight_length=1.0, oxide_width=3.0)
