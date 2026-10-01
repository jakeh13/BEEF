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
# MMI body taper (mmi_1x4_tapered_body_validation.md, Step 2)
# ---------------------------------------------------------------------------


class TestMMIBodyTaper:
    def _mmi_1x2(self, **kw):
        return MMI(n_inputs=1, n_outputs=2, mm_length=6.0, mm_width=2.5,
                   waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
                   waveguide_offset=1.0, **kw)

    # -- untapered limit: the default must reproduce the pre-taper device
    # exactly, since every existing production body is uniform. --

    def test_default_is_not_body_tapered(self):
        assert not self._mmi_1x2().is_body_tapered

    def test_default_body_end_width_is_mm_width(self):
        assert self._mmi_1x2().body_end_width == pytest.approx(2.5)

    def test_default_body_start_width_is_mm_width(self):
        assert self._mmi_1x2().body_start_width == pytest.approx(2.5)

    def test_zero_length_taper_is_not_body_tapered_even_with_end_width_set(self):
        # body_taper_length=0.0 (the default) means "no taper band exists",
        # regardless of body_end_width -- is_body_tapered requires both a
        # positive length AND a genuinely different face width.
        m = self._mmi_1x2(body_end_width=1.0)
        assert not m.is_body_tapered
        assert m.width_at(0.0) == pytest.approx(2.5)
        assert m.width_at(-3.0) == pytest.approx(2.5)

    def test_both_faces_at_mm_width_is_not_body_tapered(self):
        # A positive length but neither face differing is still a STRAIGHT
        # body -- the untapered fast path must stay reachable.
        m = self._mmi_1x2(body_taper_length=1.0,
                          body_start_width=2.5, body_end_width=2.5)
        assert not m.is_body_tapered

    def test_untapered_body_polys_is_single_rectangle(self):
        m = self._mmi_1x2()
        half = m.mm_length / 2
        polys = m._body_polys(-half, half)
        assert polys == [[(-half, -1.25), (half, -1.25), (half, 1.25), (-half, 1.25)]]

    def test_untapered_polygons_unchanged_by_new_kwargs(self):
        # byte-identical to a plain MMI() built with no body-taper kwargs
        baseline = self._mmi_1x2().polygons()
        m = self._mmi_1x2(body_taper_length=0.0, body_end_width=2.5)
        assert m.polygons() == baseline

    # -- validation --

    def test_negative_body_taper_length_raises(self):
        with pytest.raises(ValueError, match="body_taper_length must be >= 0"):
            self._mmi_1x2(body_taper_length=-1.0)

    def test_oversized_body_taper_length_raises(self):
        # 2 x 3.5 > mm_length=6.0 -- the two bands would overlap
        with pytest.raises(ValueError, match="exceeds mm_length"):
            self._mmi_1x2(body_taper_length=3.5)

    def test_nonpositive_body_end_width_raises(self):
        with pytest.raises(ValueError, match="body_end_width must be > 0"):
            self._mmi_1x2(body_taper_length=1.0, body_end_width=0.0)

    # -- tapered profile --

    def test_is_body_tapered_true(self):
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5, body_end_width=1.5)
        assert m.is_body_tapered

    def test_width_at_body_faces_is_end_width(self):
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5, body_end_width=1.5)
        half = m.mm_length / 2
        assert m.width_at(-half) == pytest.approx(1.5)
        assert m.width_at(half) == pytest.approx(1.5)

    def test_width_at_uniform_middle_is_mm_width(self):
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5, body_end_width=1.5)
        half = m.mm_length / 2
        # taper band is [-half, -half+1.0] and [half-1.0, half]; the middle
        # is uniform at mm_width right up to the band boundary.
        assert m.width_at(0.0) == pytest.approx(2.5)
        assert m.width_at(-half + 1.0) == pytest.approx(2.5)
        assert m.width_at(half - 1.0) == pytest.approx(2.5)

    def test_width_at_taper_band_midpoint_is_linear_average(self):
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5, body_end_width=1.5)
        half = m.mm_length / 2
        mid = -half + 0.5
        assert m.width_at(mid) == pytest.approx((1.5 + 2.5) / 2)

    def test_center_at_unaffected_by_body_taper(self):
        # the body taper only varies width -- centre stays on y=0 throughout,
        # same as an untapered body.
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5, body_end_width=1.5)
        half = m.mm_length / 2
        for x in (-half, -half + 0.5, 0.0, half - 0.5, half):
            assert m.center_at(x) == pytest.approx(0.0)

    def test_body_polys_is_three_trapezoids(self):
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5, body_end_width=1.5)
        half = m.mm_length / 2
        polys = m._body_polys(-half, half)
        assert len(polys) == 3

    def test_polygons_has_two_extra_trapezoids_vs_untapered(self):
        n_straight = len(self._mmi_1x2().polygons()[
            next(iter(self._mmi_1x2().polygons()))])
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5, body_end_width=1.5)
        layer = next(iter(m.polygons()))
        assert len(m.polygons()[layer]) == n_straight + 2

    # -- _check_cut must refuse a cut inside a body taper band --

    def test_section_gds_cut_inside_body_taper_band_raises(self):
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5, body_end_width=1.5)
        half = m.mm_length / 2
        with pytest.raises(ValueError, match="falls strictly inside a BODY taper band"):
            m.section_gds(-half, -half + 0.5, oxide_width=6.0)

    def test_section_gds_cut_at_body_taper_band_boundary_ok(self):
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5, body_end_width=1.5)
        half = m.mm_length / 2
        d = m.section_gds(-half, -half + 1.0, oxide_width=6.0)
        assert (2, 0) in d

    # -- asymmetric profile (2026-08-21): the two faces are independent.
    # The 1x4 needs this -- its output fan-out pins that face at mm_width
    # while its single on-axis input port leaves 81.7% of room. --

    def test_asymmetric_faces_are_independent(self):
        m = self._mmi_1x2(body_taper_length=1.0,
                          body_start_width=1.5)  # end face defaults to mm_width
        half = m.mm_length / 2
        assert m.is_body_tapered
        assert m.width_at(-half) == pytest.approx(1.5)   # narrowed input face
        assert m.width_at(half) == pytest.approx(2.5)    # output face at mm_width
        assert m.width_at(0.0) == pytest.approx(2.5)

    def test_asymmetric_ramps_only_on_the_varying_side(self):
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5)
        half = m.mm_length / 2
        # input band ramps 1.5 -> 2.5
        assert m.width_at(-half + 0.5) == pytest.approx(2.0)
        # output band is a plain rectangle at mm_width, no ramp
        assert m.width_at(half - 0.5) == pytest.approx(2.5)

    def test_a_face_at_mm_width_does_not_block_cuts(self):
        """The output band is geometrically a rectangle when that face sits at
        mm_width, so a cut inside it lands on a genuinely uniform
        cross-section and must be ALLOWED -- refusing it would forbid a valid
        section boundary, the opposite of the guard's purpose."""
        m = self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5)
        half = m.mm_length / 2
        d = m.section_gds(half - 0.5, half, oxide_width=6.0)
        assert (2, 0) in d
        # ...while the genuinely varying input band still refuses.
        with pytest.raises(ValueError, match="falls strictly inside a BODY taper band"):
            m.section_gds(-half, -half + 0.5, oxide_width=6.0)

    def test_the_real_1x4_asymmetric_profile(self):
        """The actual device this parameterisation was added for: narrow the
        input face, leave the output face at mm_width for the fan-out.
        Constants are the production 1x4's (mmi_1x4_hybrid.py)."""
        m = MMI(n_inputs=1, n_outputs=4, mm_length=39.90, mm_width=12.0,
                waveguide_width=1.1, taper_length=5.0, taper_width=2.2,
                waveguide_offset=3.1, waveguide_height=0.30,
                body_taper_length=8.0, body_start_width=6.0)
        half = m.mm_length / 2
        assert m.is_body_tapered
        assert m.width_at(-half) == pytest.approx(6.0)
        assert m.width_at(half) == pytest.approx(12.0)
        # The output face must still clear the outermost output taper:
        # |y|=4.65 + taper_width/2 = 5.75 -> needs >= 11.5 um.
        needed = 2 * (max(abs(y) for y in m.output_ys) + 2.2 / 2)
        assert m.width_at(half) >= needed


class TestMMIBodyWaist:
    """The INTERIOR-waist body profile (body_waist_width/body_waist_length/
    body_transition_length, added 2026-08-27/28).

    Both faces stay at mm_width and the narrowing happens in the interior:

        |<-stub->|<-Lt->|<--- Lw (width Ww) --->|<-Lt->|<-stub->|
           W          W->Ww           Ww          Ww->W      W

    ADDED 2026-08-29, RETROACTIVELY. This profile shipped with NO test
    coverage anywhere in tests/ despite having produced several jobs' worth of
    real physics results (the 1x4 waist sweeps). The claim recorded at the
    time -- "780 passed, 0 failed" -- only established that the feature is
    INERT to the existing suite, which it is by construction: every waist
    parameter defaults to off. It never exercised the profile itself. These
    tests close that gap and pin the contract for the 2x2 waist work.
    """

    def _mmi_1x2(self, **kw):
        kw.setdefault("mm_length", 6.0)
        return MMI(n_inputs=1, n_outputs=2, mm_width=2.5,
                   waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
                   waveguide_offset=1.0, **kw)

    def _waisted(self, **kw):
        kw.setdefault("body_waist_width", 1.5)
        kw.setdefault("body_waist_length", 2.0)
        kw.setdefault("body_transition_length", 1.0)
        return self._mmi_1x2(**kw)

    # -- untapered limit: every default-off path must reproduce a straight
    # body exactly, since every existing production body is uniform. --

    def test_default_is_not_body_waisted(self):
        assert not self._mmi_1x2().is_body_waisted

    def test_default_body_waist_width_is_mm_width(self):
        assert self._mmi_1x2().body_waist_width == pytest.approx(2.5)

    def test_zero_transition_length_is_not_waisted_even_with_waist_width_set(self):
        # body_transition_length=0.0 means "no transition band exists",
        # regardless of body_waist_width -- mirrors the body-taper rule.
        m = self._mmi_1x2(body_waist_width=1.0, body_waist_length=2.0)
        assert not m.is_body_waisted
        assert m.width_at(0.0) == pytest.approx(2.5)

    def test_waist_width_equal_to_mm_width_is_not_waisted(self):
        m = self._mmi_1x2(body_waist_width=2.5, body_waist_length=2.0,
                          body_transition_length=1.0)
        assert not m.is_body_waisted
        assert m.width_at(0.0) == pytest.approx(2.5)

    # -- the profile itself --

    def test_bounds_are_symmetric_about_centre(self):
        (din0, din1), (dout0, dout1) = self._waisted()._body_waist_bounds()
        # waist [-1.0, 1.0]; transitions [-2.0,-1.0] and [1.0, 2.0]
        assert (din0, din1) == pytest.approx((-2.0, -1.0))
        assert (dout0, dout1) == pytest.approx((1.0, 2.0))

    def test_bounds_do_not_depend_on_mm_length(self):
        """The bands are defined about x=0 from Lt/Lw alone -- mm_length only
        changes the leftover stubs. This is what makes a length sweep reuse
        the same cross-sections, and the sweep scripts rely on it."""
        a = self._waisted()._body_waist_bounds()
        b = self._waisted(mm_length=12.0)._body_waist_bounds()
        assert [x for band in a for x in band] == pytest.approx(
            [x for band in b for x in band])

    def test_widths_at_the_five_segment_boundaries(self):
        m = self._waisted()
        half = m.mm_length / 2
        assert m.is_body_waisted
        assert m.width_at(-half) == pytest.approx(2.5)    # in face
        assert m.width_at(-2.0) == pytest.approx(2.5)     # stub -> transition
        assert m.width_at(-1.0) == pytest.approx(1.5)     # transition -> waist
        assert m.width_at(0.0) == pytest.approx(1.5)      # waist centre
        assert m.width_at(1.0) == pytest.approx(1.5)      # waist -> transition
        assert m.width_at(2.0) == pytest.approx(2.5)      # transition -> stub
        assert m.width_at(half) == pytest.approx(2.5)     # out face

    def test_transitions_ramp_linearly_and_mirror_each_other(self):
        m = self._waisted()
        # halfway down the in-transition and up-transition: same width.
        assert m.width_at(-1.5) == pytest.approx(2.0)
        assert m.width_at(1.5) == pytest.approx(2.0)

    def test_faces_stay_at_mm_width_which_is_the_whole_point(self):
        """The waist exists to leave BOTH faces at mm_width, so the port
        fan-out never has to be recompromised -- unlike the face-taper."""
        m = self._waisted()
        half = m.mm_length / 2
        assert m.width_at(-half) == pytest.approx(m.mm_width)
        assert m.width_at(half) == pytest.approx(m.mm_width)

    # -- guards --

    def test_waist_and_face_taper_are_mutually_exclusive(self):
        with pytest.raises(ValueError, match="mutually exclusive"):
            self._mmi_1x2(body_taper_length=1.0, body_start_width=1.5,
                          body_waist_width=1.5, body_waist_length=1.0,
                          body_transition_length=1.0)

    def test_waist_span_cannot_exceed_mm_length(self):
        # 2*Lt + Lw = 2*2.0 + 4.0 = 8.0 > mm_length 6.0
        with pytest.raises(ValueError):
            self._mmi_1x2(body_waist_width=1.5, body_waist_length=4.0,
                          body_transition_length=2.0)

    def test_negative_and_nonpositive_parameters_are_refused(self):
        with pytest.raises(ValueError):
            self._mmi_1x2(body_transition_length=-1.0)
        with pytest.raises(ValueError):
            self._mmi_1x2(body_waist_length=-1.0)
        with pytest.raises(ValueError):
            self._mmi_1x2(body_waist_width=0.0, body_waist_length=1.0,
                          body_transition_length=1.0)

    # -- cut validation: uniform regions allow cuts, varying bands refuse. --

    def test_cut_inside_a_transition_band_is_refused(self):
        m = self._waisted()
        with pytest.raises(ValueError, match="BODY waist transition band"):
            m.section_gds(-1.5, -1.0, oxide_width=6.0)

    def test_cut_inside_the_uniform_waist_is_allowed(self):
        """The waist is a genuine rectangle at body_waist_width, so a section
        boundary there lands on a single cross-section and must be ALLOWED --
        refusing it would forbid a valid cut, the opposite of the guard's
        purpose (same reasoning as the face-taper's own flat-band test)."""
        d = self._waisted().section_gds(-0.5, 0.5, oxide_width=6.0)
        assert (2, 0) in d

    def test_cut_inside_a_uniform_stub_is_allowed(self):
        m = self._waisted()
        half = m.mm_length / 2
        d = m.section_gds(-half, -2.0, oxide_width=6.0)
        assert (2, 0) in d

    # -- the two real devices this profile is used on --

    def test_the_real_1x4_waist_profile(self):
        """The 1x4 SiN device the waist was added for (production constants,
        mmi_1x4_waist_profile_sweep.py)."""
        m = MMI(n_inputs=1, n_outputs=4, mm_length=36.91, mm_width=12.0,
                waveguide_width=1.1, taper_length=5.0, taper_width=2.2,
                waveguide_offset=3.1, waveguide_height=0.30,
                use_neighbor_aperture=True,
                body_waist_width=9.0, body_waist_length=10.0,
                body_transition_length=8.0)
        half = m.mm_length / 2
        assert m.is_body_waisted
        (din0, din1), (dout0, dout1) = m._body_waist_bounds()
        assert [din0, din1, dout0, dout1] == pytest.approx(
            [-13.0, -5.0, 5.0, 13.0])
        assert m.width_at(-half) == pytest.approx(12.0)
        assert m.body_width_at(0.0) == pytest.approx(9.0)
        assert m.width_at(half) == pytest.approx(12.0)

    def test_the_real_2x2_waist_profile_and_its_zero_face_headroom(self):
        """The 2x2 SOI device. Constants are the production
        straight 2x2's (mmi_2x2_hybrid.py), which is already fully resolved.

        The load-bearing assertion is the SECOND one: this device's two port
        tapers consume 100.0% of the body face (2.800 um of 2.800 um, zero
        headroom), so the face-taper profile is not merely unattractive here
        but geometrically IMPOSSIBLE -- the interior waist is the only body
        profile the device can host. Re-derived from output_ys rather than
        hardcoded, so a future port-layout change fails this test instead of
        silently invalidating the premise of that whole validation doc."""
        taper_width = 1.3
        m = MMI(n_inputs=2, n_outputs=2, mm_length=28.0, mm_width=2.8,
                waveguide_width=0.45, taper_length=4.0, taper_width=taper_width,
                waveguide_offset=1.5, waveguide_height=0.22,
                body_waist_width=2.24, body_waist_length=8.0,
                body_transition_length=6.0)
        half = m.mm_length / 2
        assert m.is_body_waisted
        assert m.width_at(-half) == pytest.approx(2.8)
        assert m.body_width_at(0.0) == pytest.approx(2.24)
        assert m.width_at(half) == pytest.approx(2.8)

        # Zero face headroom: the ports span the entire face.
        needed = 2 * (max(abs(y) for y in m.output_ys) + taper_width / 2)
        assert needed == pytest.approx(2.8)
        assert m.mm_width == pytest.approx(needed)


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


class TestCompoundGeometrySourceApertureClamps:
    """`source_aperture_at()`'s out-of-bounds x is deliberately NOT a caller
    error, unlike width_at/center_at/section_prisms/section_gds above.

    Found 2026-07-31: `FDTDSolver.solve()` (the full single-domain path, as
    opposed to solve_section()'s per-section calls, which HybridSolver keeps
    well inside a section's own bounds) queries source_aperture_at() at the
    solver's external stub position -- `waveguide_length` past the device's
    own physical extent, for PML clearance. Every leaf geometry already
    tolerates this by construction (MMI.width_at() clamps to the terminal
    waveguide_width for any x past its taper start, however far past), but
    CompoundGeometry._locate() had no segment "out there" to delegate to and
    raised instead -- the first script to call solve() (not solve_section())
    with a CompoundGeometry as the top-level device (modal_multiplexer_fdtd.py)
    hit it immediately, on its very first FDTD run, before this fix.
    """

    def _compound(self):
        # First segment: a 2-output MMI with use_neighbor_aperture on, so its
        # source_aperture_at() is the interesting (non-None, width_at-driven)
        # case rather than a Taper/WaveguideRun's unconditional None.
        mmi = MMI(
            n_inputs=1, n_outputs=2,
            mm_length=4.0, mm_width=2.0, waveguide_width=0.35,
            taper_length=1.0, taper_width=0.8, waveguide_offset=1.0,
            use_neighbor_aperture=True,
        )
        run = WaveguideRun(length=2.0, ports=[(0.5, 0.35), (-0.5, 0.35)])
        return CompoundGeometry([mmi, run]), mmi

    def test_within_bounds_unaffected_by_the_clamp_flag(self):
        """clamp=True must not change any answer that clamp=False already
        gave -- it only extends behaviour past the edges, never inside."""
        g, _ = self._compound()
        lo, hi = -g.total_length / 2, g.total_length / 2
        for x in (lo, lo / 2, 0.0, hi / 2, hi):
            seg, seg_lo, _ = g._locate(x, clamp=False)
            seg_clamped, seg_lo_clamped, _ = g._locate(x, clamp=True)
            assert seg is seg_clamped and seg_lo == seg_lo_clamped

    def test_past_the_input_end_delegates_instead_of_raising(self):
        g, mmi = self._compound()
        half_total = g.total_length / 2
        x_stub = -half_total - 0.5   # solver's external stub, past the MMI's own input taper
        # Must not raise (this is exactly what crashed before the fix), and
        # must genuinely delegate to the end segment (mmi), not just avoid
        # crashing some other way.
        seg, _, _ = g._locate(x_stub, clamp=True)
        assert seg is mmi
        # Since MMI.width_at()/source_aperture_at() clamp to the terminal
        # value for any x past the input taper start, querying deep in the
        # stub must equal querying right at the segment's own edge -- both
        # land in the same "clamped to terminal" regime.
        aperture = g.source_aperture_at(x_stub)
        at_edge = g.source_aperture_at(-half_total)
        assert aperture == pytest.approx(at_edge)

    def test_past_the_output_end_delegates_instead_of_raising(self):
        g, _ = self._compound()
        half_total = g.total_length / 2
        x_stub = half_total + 0.5   # past the WaveguideRun's own end
        # WaveguideRun.source_aperture_at() is always None (checked via
        # isinstance in CompoundGeometry.source_aperture_at()), so the only
        # claim here is that it does not raise.
        assert g.source_aperture_at(x_stub) is None

    def test_far_outside_still_delegates_to_the_same_end_segment(self):
        """Not just 'doesn't raise at a small overshoot' -- arbitrarily far
        past the edge must keep resolving to the same end segment, matching
        solve()'s stub position however large waveguide_length+dpml is."""
        g, _ = self._compound()
        half_total = g.total_length / 2
        near = g.source_aperture_at(-half_total - 0.5)
        far = g.source_aperture_at(-half_total - 50.0)
        assert near == pytest.approx(far)

    def test_width_at_and_center_at_are_NOT_clamped_by_default(self):
        """The opt-in is scoped to source_aperture_at() only -- width_at()/
        center_at() must keep raising outside bounds, since an out-of-bounds
        cut there is a genuine caller error building physical geometry, not
        a solver's legitimate external-stub query."""
        g, _ = self._compound()
        half_total = g.total_length / 2
        with pytest.raises(ValueError, match="outside CompoundGeometry"):
            g.width_at(-half_total - 0.5)
        with pytest.raises(ValueError, match="outside CompoundGeometry"):
            g.center_at(-half_total - 0.5)


class TestCompoundGeometryStubIsAppliedExactlyOnce:
    """Regression guards for a real stub double-application bug.

    `_overlaps()` expands the query window by the stubs itself
    (`lo, hi = x0-stub_left, x1+stub_right`) before clipping each segment's
    overlap to it, so `ov_lo`/`ov_hi` are already stub-inclusive. It ALSO
    used to hand the raw stub values down to the first/last overlapping
    segment, which applied them a SECOND time -- drawing geometry 2x the
    stub too long whenever the padded window landed fully inside one
    segment (nothing clipped it back). Confirmed twice in the wild before
    being fixed: `audit_interface_modes.py` had to bypass
    CompoundGeometry entirely, and HybridSolver's mode-correspondence check
    tripped a spurious "cut falls strictly inside a taper" from the
    double-subtracted coordinate.

    The invariant these tests pin: the total drawn x-extent must equal
    `(x1 + stub_right) - (x0 - stub_left)` exactly, whether the padded
    window lands inside one segment or spans several.
    """

    def _compound(self):
        # seg2, the WaveguideRun, spans [-1.5, 1.5] in device frame.
        return CompoundGeometry([
            Taper(length=4.0, width_start=0.35, width_end=0.8),
            WaveguideRun(length=3.0, ports=[(0.0, 0.8)]),
            Taper(length=4.0, width_start=0.8, width_end=0.35),
        ])

    @staticmethod
    def _x_extent(polys):
        xs = [v[0] for poly in polys for v in poly]
        return min(xs), max(xs)

    def test_section_gds_stub_not_doubled_when_window_is_inside_one_segment(self):
        g = self._compound()
        # [-0.5, 0.5] padded by 0.5 both sides -> [-1.0, 1.0], entirely
        # within seg2, so nothing clips the stub expansion back.
        d = g.section_gds(-0.5, 0.5, oxide_width=3.0, stub_left=0.5, stub_right=0.5)
        lo_x, hi_x = self._x_extent(d[(2, 0)])
        assert hi_x - lo_x == pytest.approx(2.0)  # was 3.0 with the stub doubled
        # section_gds returns a 0-based section frame, so it must start at 0
        # and line up with the oxide box (which was always correct).
        assert lo_x == pytest.approx(0.0)
        ox_lo, ox_hi = self._x_extent(d[(1, 0)])
        assert (lo_x, hi_x) == pytest.approx((ox_lo, ox_hi))

    def test_section_prisms_stub_not_doubled_when_window_is_inside_one_segment(self):
        pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
        g = self._compound()
        prisms = g.section_prisms(-0.5, 0.5, si_medium=None, stub_left=0.5, stub_right=0.5,
                                  wg_height=0.22)
        xs = [v.x for p in prisms for v in p.vertices]
        # Device frame here (not 0-based): must be exactly [-1.0, 1.0].
        assert min(xs) == pytest.approx(-1.0)
        assert max(xs) == pytest.approx(1.0)

    def test_zero_stub_is_unchanged(self):
        """The stub=0 path is what every production HybridSolver interior
        boundary uses (hybrid.py's _stubs()), so it must be untouched."""
        g = self._compound()
        d = g.section_gds(-0.5, 0.5, oxide_width=3.0)
        lo_x, hi_x = self._x_extent(d[(2, 0)])
        assert hi_x - lo_x == pytest.approx(1.0)
        assert lo_x == pytest.approx(0.0)

    def test_stub_reaching_outside_a_segment_is_still_drawn(self):
        """The stub must NOT simply be dropped: when the padded window
        genuinely reaches past a segment's edge, that part is a real
        access-waveguide extension the delegate segment still has to draw.
        Here the window spans seg1/seg2/seg3, so the outer stubs fall
        outside the compound device entirely."""
        g = self._compound()
        half = g.total_length / 2  # 5.5
        d = g.section_gds(-half, half, oxide_width=3.0, stub_left=0.5, stub_right=0.5)
        lo_x, hi_x = self._x_extent(d[(2, 0)])
        assert hi_x - lo_x == pytest.approx(g.total_length + 1.0)  # 11.0 + both stubs
        assert lo_x == pytest.approx(0.0)


class TestCompoundGeometryMMIFrameOffset:
    """Regression guard: CompoundGeometry assumed
    every segment's local frame is centered on its own bbox
    (`(seg_lo + seg_hi) / 2`), but MMI's frame is centered on the MM body
    instead -- the two only coincide when
    max(input taper lengths) == max(output taper lengths). Reproduction from
    the notes: mm_length=5 with one input taper of 3.0 and one output taper
    of 1.0 (asymmetric). MMI.total_length is 9 either way, but MMI's own
    frame spans [-5.5, 3.5] (bbox-centered would be [-4.5, 4.5]), a genuine
    `frame_offset` of (3.0 - 1.0) / 2 = 1.0. Before the fix, `compound.
    polygons()` hung 1 um past its own declared left edge and `width_at()` at
    that edge returned a mid-taper width instead of the true access-waveguide
    width.
    """

    def _compound(self):
        mmi = MMI(
            n_inputs=1, n_outputs=1, mm_length=5.0, mm_width=2.0,
            inputs=[MMIPort(y=0.0, waveguide_width=0.35, taper_length=3.0, taper_width=0.8)],
            outputs=[MMIPort(y=0.0, waveguide_width=0.35, taper_length=1.0, taper_width=0.8)],
        )
        run = WaveguideRun(length=2.0, ports=[(0.0, 0.35)])
        return CompoundGeometry([mmi, run]), mmi

    def test_mmi_frame_offset_is_half_the_taper_asymmetry(self):
        _, mmi = self._compound()
        assert mmi.total_length == pytest.approx(9.0)
        assert mmi.frame_offset == pytest.approx((3.0 - 1.0) / 2)

    def test_polygons_stay_within_declared_bbox(self):
        g, _ = self._compound()
        half = g.total_length / 2  # 5.5
        polys = g.polygons()[(2, 0)]
        xs = [v[0] for poly in polys for v in poly]
        assert min(xs) == pytest.approx(-half)
        assert max(xs) == pytest.approx(half)

    def test_width_at_left_edge_is_the_access_waveguide_not_mid_taper(self):
        g, _ = self._compound()
        half = g.total_length / 2  # 5.5
        # True left edge of the MMI's own frame: its input access-waveguide
        # width (0.35), not a mid-taper value (0.467) from a bbox-mis-centered
        # lookup.
        assert g.width_at(-half) == pytest.approx(0.35)

    def test_section_gds_stays_within_declared_bbox(self):
        g, _ = self._compound()
        half = g.total_length / 2
        d = g.section_gds(-half, half, oxide_width=3.0)
        xs = [v[0] for poly in d[(2, 0)] for v in poly]
        assert max(xs) - min(xs) == pytest.approx(2 * half)


class TestCompoundGeometryPolygonsStub:
    """Regression guard: `polygons(
    stub_length>0)` used to raise `NotImplementedError` unconditionally,
    on the premise that an arbitrary segment type has no single well-defined
    port width to extend. That premise was wrong -- `_end_ports()` already
    answers the question for `to_gds_eme()`'s own stub handling by returning
    EVERY port at that end rather than picking one, and a leaf geometry with
    several ports already draws one stub per port for its own `polygons()`
    (e.g. `MMI.polygons()`). The fix reuses `_end_ports()` the same way
    `to_gds_eme()` already does. This is also what unblocks
    `modal_multiplexer_fdtd.py` (`modal_multiplexer_validation.md` gap 3):
    `FDTDSolver.solve()`'s single-domain path calls exactly this method.
    """

    def _simple_compound(self):
        # Same fixture as TestCompoundGeometry._compound(): single-port ends.
        return CompoundGeometry([
            Taper(length=4.0, width_start=0.35, width_end=0.8),
            WaveguideRun(length=3.0, ports=[(0.0, 0.8)]),
            Taper(length=4.0, width_start=0.8, width_end=0.35),
        ])

    def test_no_longer_raises(self):
        g = self._simple_compound()
        g.polygons(stub_length=1.5)  # must not raise NotImplementedError

    def test_stub_length_zero_is_unchanged(self):
        g = self._simple_compound()
        d0 = g.polygons(stub_length=0.0)
        d_explicit = g.polygons()
        assert d0[(2, 0)] == d_explicit[(2, 0)]

    def test_single_port_ends_add_exactly_one_stub_polygon_each(self):
        g = self._simple_compound()
        d_bare = g.polygons(stub_length=0.0)
        d_stub = g.polygons(stub_length=1.5)
        assert len(d_stub[(2, 0)]) == len(d_bare[(2, 0)]) + 2  # one per end

    def test_stub_polygons_sit_at_the_true_outer_edges_with_terminal_width(self):
        g = self._simple_compound()
        half = g.total_length / 2  # 5.5 -- seg1 width_start=0.35, seg3 width_end=0.35
        d_bare = g.polygons(stub_length=0.0)
        n_bare = len(d_bare[(2, 0)])
        stubs = g.polygons(stub_length=1.5)[(2, 0)][n_bare:]
        assert len(stubs) == 2
        xs_by_poly = [[v[0] for v in poly] for poly in stubs]
        left, right = sorted(xs_by_poly, key=min)
        assert min(left) == pytest.approx(-half - 1.5)
        assert max(left) == pytest.approx(-half)
        assert min(right) == pytest.approx(half)
        assert max(right) == pytest.approx(half + 1.5)
        for poly in stubs:
            ys = [v[1] for v in poly]
            assert max(ys) - min(ys) == pytest.approx(0.35)  # terminal width both ends

    def test_multi_port_end_adds_one_stub_per_port_not_one(self):
        """The case this fix exists for: a multi-port end (here, an MMI's 3
        outputs) must get one stub polygon per port, at each port's own
        width and y -- not a single arbitrarily-chosen winner."""
        mmi_in = MMI(
            n_inputs=1, n_outputs=1, mm_length=5.0, mm_width=2.0,
            inputs=[MMIPort(y=0.0, waveguide_width=0.5, taper_length=2.0, taper_width=1.0)],
            outputs=[MMIPort(y=0.0, waveguide_width=0.5, taper_length=2.0, taper_width=1.0)],
        )
        branch_specs = [(+1.0, 0.33, 0.7), (0.0, 0.35, 1.1), (-1.0, 0.33, 0.7)]
        mmi_out = MMI(
            n_inputs=1, n_outputs=3, mm_length=6.0, mm_width=3.0,
            inputs=[MMIPort(y=0.0, waveguide_width=0.5, taper_length=2.0, taper_width=1.0)],
            outputs=[MMIPort(y=y, waveguide_width=ww, taper_length=2.0, taper_width=tw)
                     for y, ww, tw in branch_specs],
        )
        g = CompoundGeometry([mmi_in, mmi_out])
        d_bare = g.polygons(stub_length=0.0)
        n_bare = len(d_bare[(2, 0)])
        stubs = g.polygons(stub_length=0.75)[(2, 0)][n_bare:]
        # 1 input stub (mmi_in's single input) + 3 output stubs (mmi_out's
        # three branches) = 4, not 1.
        assert len(stubs) == 4
        half = g.total_length / 2
        widths_at_right_edge = sorted(
            round(max(v[1] for v in poly) - min(v[1] for v in poly), 6)
            for poly in stubs if min(v[0] for v in poly) == pytest.approx(half)
        )
        assert widths_at_right_edge == sorted([0.33, 0.35, 0.33])
        # Every output stub starts exactly at the true right edge and
        # extends outward by stub_length, never inward.
        for poly in stubs:
            xs = [v[0] for v in poly]
            if min(xs) == pytest.approx(half):
                assert max(xs) == pytest.approx(half + 0.75)

    def test_input_stub_respects_mmi_frame_offset(self):
        """Same asymmetric-taper MMI as TestCompoundGeometryMMIFrameOffset:
        frame_offset != 0, so a naive bbox-midpoint stub placement would land
        1 um off the true left edge. The stub must land at the true edge."""
        mmi = MMI(
            n_inputs=1, n_outputs=1, mm_length=5.0, mm_width=2.0,
            inputs=[MMIPort(y=0.0, waveguide_width=0.35, taper_length=3.0, taper_width=0.8)],
            outputs=[MMIPort(y=0.0, waveguide_width=0.35, taper_length=1.0, taper_width=0.8)],
        )
        run = WaveguideRun(length=2.0, ports=[(0.0, 0.35)])
        g = CompoundGeometry([mmi, run])
        half = g.total_length / 2  # 5.5, the TRUE (frame_offset-corrected) edge
        d_bare = g.polygons(stub_length=0.0)
        n_bare = len(d_bare[(2, 0)])
        stubs = g.polygons(stub_length=1.0)[(2, 0)][n_bare:]
        left_stub = min(stubs, key=lambda poly: min(v[0] for v in poly))
        xs = [v[0] for v in left_stub]
        assert min(xs) == pytest.approx(-half - 1.0)
        assert max(xs) == pytest.approx(-half)  # not -half+1.0 (the bbox-only midpoint bug)
        ys = [v[1] for v in left_stub]
        assert max(ys) - min(ys) == pytest.approx(0.35)  # true input access width

    def test_to_meep_prisms_no_longer_raises(self):
        """The actual crash chain this fix clears: FDTDSolver.solve() ->
        to_meep_prisms() -> polygons(stub_length>0). modal_multiplexer_fdtd.py
        was the first script in the repo to hit this."""
        pytest.importorskip("meep")
        import meep as mp
        g = self._simple_compound()
        prisms = g.to_meep_prisms(mp.Medium(index=3.45), stub_length=0.75, wg_height=0.22)
        assert len(prisms) > 0


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

    @pytest.mark.parametrize("side", ["left", "right"])
    def test_stub_padding_draws_waveguide_not_bare_cladding(self, side):
        """The bus and the seam-side sliver must reach the padded boundary.

        solve_section() places the source at lo=x0-stub_left and the monitor
        near hi=x1+stub_right, NOT at x0/x1. Until 2026-08-02 this method
        clamped the bus to [outer_x, iface_x], so both planes sat on bare
        cladding and the eigenmode launch had no guide to couple into.
        """
        c = self._left() if side == "left" else self._right()
        x0, x1 = sorted((c.outer_x, c.iface_x))
        pad = 1.5
        lo, hi = x0 - pad, x1 + pad
        polys = c._section_polygons(x0, x1, stub_left=pad, stub_right=pad, y_pml=0.5)

        def _spans_y(yc):
            # 4-vertex rects only; the near-bend arc is a 128-vertex polygon
            # that also straddles y_ring and is not what this test is about.
            return [(min(v[0] for v in p), max(v[0] for v in p)) for p in polys
                    if len(p) == 4 and min(v[1] for v in p) <= yc <= max(v[1] for v in p)]

        bus = _spans_y(0.0)
        assert len(bus) == 1
        assert bus[0] == pytest.approx((lo, hi))     # bus spans the FULL padded domain

        # the sliver follows the padding out through the seam (iface_x) end
        # only; its other end stops at hc, where the bend takes over.
        sliver = _spans_y(c.y_ring)
        assert len(sliver) == 1
        seam_end = hi if c.iface_x > c.hc else lo
        assert sliver[0] == pytest.approx(tuple(sorted((c.hc, seam_end))))

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
