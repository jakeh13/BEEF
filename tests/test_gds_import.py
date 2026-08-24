"""Fast, meep/meow-free unit tests for beef.geometry.gds_import.

Writes a small synthetic 2-layer GDS (mimicking the elevator device: a
straight + a taper on layer (2,0), a taper + a straight on layer (3738,0),
overlapping in x) with gdstk, then round-trips it through
load_gds_polygons()/GdsGeometry.
"""
from __future__ import annotations

import pytest

gdstk = pytest.importorskip("gdstk")

from beef.geometry.gds_import import GdsGeometry, load_gds_polygons

LAYER_SI = (2, 0)
LAYER_SIN = (3738, 0)

STUB = 2.0
L = 6.0
SPAN = 2 * STUB + L  # 10.0


@pytest.fixture
def elevator_gds(tmp_path):
    lib = gdstk.Library()
    cell = lib.new_cell("elevator_test")

    # layer (2,0): straight [0, L/2+stub] then taper down to [L/2+stub, L+stub]
    cell.add(gdstk.rectangle((0, -0.175), (L / 2 + STUB, 0.175), layer=2, datatype=0))
    cell.add(gdstk.Polygon(
        [(L / 2 + STUB, -0.175), (L + STUB, -0.045), (L + STUB, 0.045), (L / 2 + STUB, 0.175)],
        layer=2, datatype=0,
    ))
    # layer (3738,0): taper [stub, L/2+stub] then straight [L/2+stub, L+2*stub]
    cell.add(gdstk.Polygon(
        [(STUB, -0.06), (L / 2 + STUB, -0.4), (L / 2 + STUB, 0.4), (STUB, 0.06)],
        layer=3738, datatype=0,
    ))
    cell.add(gdstk.rectangle((L / 2 + STUB, -0.4), (L + 2 * STUB, 0.4), layer=3738, datatype=0))

    path = tmp_path / "elevator_test.gds"
    lib.write_gds(str(path))
    return path


def test_load_gds_polygons_groups_by_layer(elevator_gds):
    polys = load_gds_polygons(elevator_gds)
    assert set(polys.keys()) == {LAYER_SI, LAYER_SIN}
    assert len(polys[LAYER_SI]) == 2
    assert len(polys[LAYER_SIN]) == 2
    for layer_polys in polys.values():
        for poly in layer_polys:
            assert all(isinstance(v, tuple) and len(v) == 2 for v in poly)


def test_load_gds_polygons_layer_filter(elevator_gds):
    polys = load_gds_polygons(elevator_gds, layers=[LAYER_SI])
    assert set(polys.keys()) == {LAYER_SI}


def test_load_gds_polygons_missing_cell_name_raises(elevator_gds):
    with pytest.raises(ValueError):
        load_gds_polygons(elevator_gds, cell_name="does_not_exist")


class TestGdsGeometry:
    def test_total_length_matches_bbox_span(self, elevator_gds):
        geom = GdsGeometry(
            elevator_gds, core_layers=[LAYER_SI, LAYER_SIN],
            region_lengths=[STUB, L, STUB],
        )
        assert geom.total_length == pytest.approx(SPAN)

    def test_eme_lengths_returns_region_lengths_unchanged(self, elevator_gds):
        geom = GdsGeometry(
            elevator_gds, core_layers=[LAYER_SI, LAYER_SIN],
            region_lengths=[STUB, L, STUB],
        )
        assert geom.eme_lengths(straight_length=99.0) == [STUB, L, STUB]

    def test_mismatched_region_lengths_raises(self, elevator_gds):
        with pytest.raises(ValueError):
            GdsGeometry(
                elevator_gds, core_layers=[LAYER_SI, LAYER_SIN],
                region_lengths=[STUB, L],  # missing the trailing stub
            )

    def test_bbox_centered_at_x_zero(self, elevator_gds):
        geom = GdsGeometry(
            elevator_gds, core_layers=[LAYER_SI, LAYER_SIN],
            region_lengths=[STUB, L, STUB],
        )
        xs = [x for polys in geom._polys.values() for poly in polys for x, _y in poly]
        assert min(xs) == pytest.approx(-SPAN / 2)
        assert max(xs) == pytest.approx(SPAN / 2)

    def test_to_gds_eme_places_geometry_at_x_start_and_adds_oxide_box(self, elevator_gds):
        geom = GdsGeometry(
            elevator_gds, core_layers=[LAYER_SI, LAYER_SIN],
            region_lengths=[STUB, L, STUB],
        )
        out = geom.to_gds_eme(straight_length=0.0, oxide_width=5.0, x_start=0.0)
        assert set(out.keys()) == {LAYER_SI, LAYER_SIN, (1, 0)}
        xs = [x for poly in out[LAYER_SI] + out[LAYER_SIN] for x, _y in poly]
        assert min(xs) == pytest.approx(0.0, abs=1e-9)
        assert max(xs) == pytest.approx(SPAN, abs=1e-9)
        ox_xs = [x for x, _y in out[(1, 0)][0]]
        assert min(ox_xs) == pytest.approx(0.0)
        assert max(ox_xs) == pytest.approx(SPAN)

    def test_to_gds_eme_ignores_straight_length(self, elevator_gds):
        """Unlike Taper/SBend/MMI, the imported geometry already contains its
        own access waveguides — to_gds_eme() must not append more."""
        geom = GdsGeometry(
            elevator_gds, core_layers=[LAYER_SI, LAYER_SIN],
            region_lengths=[STUB, L, STUB],
        )
        out_zero = geom.to_gds_eme(straight_length=0.0, oxide_width=5.0, x_start=0.0)
        out_nonzero = geom.to_gds_eme(straight_length=3.5, oxide_width=5.0, x_start=0.0)
        assert out_zero[LAYER_SI] == out_nonzero[LAYER_SI]
        assert out_zero[LAYER_SIN] == out_nonzero[LAYER_SIN]
        assert out_zero[(1, 0)] == out_nonzero[(1, 0)]

    def test_no_matching_layers_raises(self, elevator_gds):
        with pytest.raises(ValueError):
            GdsGeometry(elevator_gds, core_layers=[(99, 0)], region_lengths=[SPAN])


# ---------------------------------------------------------------------------
# HybridSolver section hooks (resolved/gds_geometry_hybrid_cascade_validation.md).
#
# The validation vehicle that doc's step 1 specifies is a GDS twin of the
# shipped 2×2 MMI, so these tests build the twin the same way its generator
# script does — from MMI.polygons() itself, exact by construction — and
# compare every hook against MMI's own answer. Comparing against MMI rather
# than against hardcoded numbers is the point: it is the *query* path
# (width_at/center_at/section_*) under test, not whether two independently
# drawn layouts agree.
#
# meep/meow-free like the rest of this file: section_gds() is polygons in,
# polygons out, and section_prisms() is the only hook needing meep (covered
# by asserting it delegates to the same _section_polygons()).
# ---------------------------------------------------------------------------

from beef.geometry.mmi import MMI  # noqa: E402

TWIN_PARAMS = dict(
    n_inputs=2, n_outputs=2, mm_length=28.0, mm_width=2.8,
    waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
    waveguide_offset=1.5, waveguide_height=0.22,
)
TWIN_SPAN = 36.0          # mm_length + 2 * taper_length
TWIN_HALF_BODY = 14.0
TWIN_REGIONS = [4.0, 28.0, 4.0]


def _polygon_area(polys):
    total = 0.0
    for poly in polys:
        s = 0.0
        for i in range(len(poly)):
            x0, y0 = poly[i]
            x1, y1 = poly[(i + 1) % len(poly)]
            s += x0 * y1 - x1 * y0
        total += abs(s) / 2.0
    return total


@pytest.fixture
def mmi_twin_gds(tmp_path):
    """A GDS whose (2,0) polygons ARE MMI.polygons(stub_length=0) — the same
    file examples/gds_files/make_mmi_2x2_twin.py writes."""
    lib = gdstk.Library()
    cell = lib.new_cell("mmi_2x2_twin")
    for verts in MMI(**TWIN_PARAMS).polygons(stub_length=0.0)[LAYER_SI]:
        cell.add(gdstk.Polygon(verts, layer=2, datatype=0))
    path = tmp_path / "mmi_2x2_twin.gds"
    lib.write_gds(str(path))
    return path


@pytest.fixture
def twin(mmi_twin_gds):
    return GdsGeometry(
        mmi_twin_gds, core_layers=[LAYER_SI], region_lengths=TWIN_REGIONS,
        access_extension=True,
    )


@pytest.fixture
def twin_parametric():
    return MMI(**TWIN_PARAMS)


class TestSectionHooks:
    def test_twin_span_matches_mmi_total_length(self, twin, twin_parametric):
        """The two device frames must coincide, or no boundary transfers:
        GdsGeometry centres the imported bbox on x=0, MMI's origin is its
        body centre, and those are the same point only when frame_offset is
        0 (equal in/out taper lengths)."""
        assert twin.total_length == pytest.approx(twin_parametric.total_length)
        assert twin_parametric.frame_offset == pytest.approx(0.0)

    @pytest.mark.parametrize("kwargs", [
        {"input_idx": 1, "output_idx": 1},
        {"input_idx": 0, "output_idx": 0},
        {"input_idx": 1, "output_idx": 0},
    ])
    def test_width_and_center_match_mmi_densely(self, twin, twin_parametric, kwargs):
        """Dense sweep across the access-waveguide region, BOTH taper bands
        (where the width varies continuously) and the body."""
        n = 1401
        x_lo, x_hi = -TWIN_SPAN / 2 - 2.0, TWIN_SPAN / 2 + 2.0
        for i in range(n):
            x = x_lo + (x_hi - x_lo) * i / (n - 1)
            assert twin.width_at(x, **kwargs) == pytest.approx(
                twin_parametric.width_at(x, **kwargs), abs=1e-9), f"width_at({x})"
            assert twin.center_at(x, **kwargs) == pytest.approx(
                twin_parametric.center_at(x, **kwargs), abs=1e-9), f"center_at({x})"

    def test_center_at_is_not_the_base_class_zero(self, twin):
        """The quiet failure mode the plan calls out: DeviceGeometry.center_at
        returns 0.0 unconditionally, which is wrong for every off-axis branch
        and raises nothing."""
        assert twin.center_at(-TWIN_SPAN / 2, input_idx=1) == pytest.approx(0.75)
        assert twin.center_at(-TWIN_SPAN / 2, input_idx=0) == pytest.approx(-0.75)
        assert twin.center_at(+TWIN_SPAN / 2, output_idx=1) == pytest.approx(0.75)

    def test_width_at_reports_the_full_body_face_at_the_taper_seam(self, twin):
        """x exactly on the body/taper seam cuts the body polygon's own
        vertical edge. A half-open scanline rule drops it and reports a
        taper's 1.3 µm; the closed rule reports the body's 2.8 µm."""
        assert twin.width_at(-TWIN_HALF_BODY) == pytest.approx(2.8)
        assert twin.width_at(+TWIN_HALF_BODY) == pytest.approx(2.8)

    def test_body_ignores_branch_index(self, twin, twin_parametric):
        """One branch at this cross-section, so there is nothing to select —
        matching MMI, whose width_at returns mm_width inside the body
        regardless of input_idx/output_idx."""
        for idx in (0, 1):
            assert twin.width_at(0.0, input_idx=idx) == pytest.approx(2.8)
            assert twin.center_at(0.0, input_idx=idx) == pytest.approx(0.0)

    def test_out_of_range_branch_index_raises_with_the_branches_listed(self, twin):
        with pytest.raises(ValueError, match="out of range"):
            twin.width_at(-TWIN_SPAN / 2, input_idx=7)

    def test_beyond_the_bbox_holds_the_facet_cross_section(self, twin, twin_parametric):
        """The (virtual, arbitrary-length) access-waveguide region every
        parametric geometry reports in too."""
        for x in (-TWIN_SPAN / 2 - 0.001, -TWIN_SPAN / 2 - 5.0, -TWIN_SPAN / 2 - 500.0):
            assert twin.width_at(x, input_idx=1) == pytest.approx(0.45)
            assert twin.center_at(x, input_idx=1) == pytest.approx(0.75)
            assert twin.width_at(x, input_idx=1) == pytest.approx(
                twin_parametric.width_at(x, input_idx=1))

    def test_source_aperture_matches_mmi_when_both_opt_in(self, mmi_twin_gds):
        gds = GdsGeometry(mmi_twin_gds, core_layers=[LAYER_SI],
                          region_lengths=TWIN_REGIONS, access_extension=True,
                          neighbor_aperture=True)
        mmi = MMI(**TWIN_PARAMS, use_neighbor_aperture=True)
        # 1e-4, not 1e-9: a probe within the scanline's own vertex-snap
        # tolerance (_X_TOL=1e-9) of a seam is exactly the ambiguous case
        # the plan's step 2 warns about ("prefer sampling just inside the
        # span over exactly at a boundary") -- stay safely clear of it.
        for x in (-TWIN_SPAN / 2, -17.0, -15.0, -TWIN_HALF_BODY - 1e-4,
                  TWIN_HALF_BODY + 1e-4, 15.0, TWIN_SPAN / 2):
            a = mmi.source_aperture_at(x, input_idx=1, output_idx=1)
            b = gds.source_aperture_at(x, input_idx=1, output_idx=1)
            assert (a is None) == (b is None), f"x={x}"
            if a is not None:
                assert b == pytest.approx(a, abs=1e-9), f"x={x}"

    def test_source_aperture_is_opt_out_by_default(self, twin):
        """Must stay None by default, or half_domain_aperture stops winning
        and every existing run's launch window silently changes."""
        assert twin.source_aperture_at(-TWIN_SPAN / 2, input_idx=1) is None


class TestSectionGeometry:
    """section_gds() against MMI's, for the three sections
    examples/mmi_2x2/python/mmi_2x2_hybrid_gds.py actually cascades."""

    BODY_STUB = 0.25
    STUB = 1.0  # FDTDSolver.waveguide_length (0.5) + dpml (0.5)

    def _bounds(self):
        half_total = TWIN_SPAN / 2
        return [-half_total, -TWIN_HALF_BODY + self.BODY_STUB,
                TWIN_HALF_BODY - self.BODY_STUB, half_total]

    @pytest.mark.parametrize("k,stub_left,stub_right,kwargs", [
        (0, STUB, 0.0, {"input_idx": 1}),
        (1, 0.0, 0.0, {}),
        (2, 0.0, STUB, {"output_idx": 1}),
    ])
    def test_section_gds_matches_mmi(self, twin, twin_parametric,
                                     k, stub_left, stub_right, kwargs):
        b = self._bounds()
        x0, x1 = b[k], b[k + 1]
        a = twin_parametric.section_gds(x0, x1, 6.0, stub_left=stub_left,
                                        stub_right=stub_right, **kwargs)
        g = twin.section_gds(x0, x1, 6.0, stub_left=stub_left,
                             stub_right=stub_right, **kwargs)
        assert set(a) == set(g)
        assert a[(1, 0)] == g[(1, 0)]                       # oxide box
        assert len(g[LAYER_SI]) == len(a[LAYER_SI])         # not unioned together
        assert _polygon_area(g[LAYER_SI]) == pytest.approx(
            _polygon_area(a[LAYER_SI]), abs=1e-9)
        # Local frame: left edge of the stubbed span at x=0, same as MMI's.
        xs_g = [x for poly in g[LAYER_SI] for x, _y in poly]
        xs_a = [x for poly in a[LAYER_SI] for x, _y in poly]
        assert min(xs_g) == pytest.approx(min(xs_a), abs=1e-9)
        assert max(xs_g) == pytest.approx(max(xs_a), abs=1e-9)
        assert min(xs_g) == pytest.approx(0.0, abs=1e-9)

    def test_section_draws_every_branch_not_just_the_selected_one(self, twin):
        """input_idx selects which branch width_at/center_at report; the
        drawn section is a real cross-section of the whole device, so the
        unexcited arm is still sitting in it."""
        b = self._bounds()
        g = twin.section_gds(b[0], b[1], 6.0, stub_left=self.STUB, input_idx=1)
        ys = [y for poly in g[LAYER_SI] for _x, y in poly]
        assert min(ys) < -0.5 and max(ys) > 0.5

    def test_access_extension_synthesizes_both_facet_branches(self, twin):
        """The stub past the physical port is two waveguides here, not one —
        the generalization of the one-sided input_width mechanism."""
        b = self._bounds()
        g = twin.section_gds(b[0], b[1], 6.0, stub_left=self.STUB, input_idx=1)
        # In local frame the stub occupies x in [0, STUB]; sample its middle.
        stub_mid_dev = b[0] - self.STUB / 2
        assert twin._intervals_at(stub_mid_dev) == twin._facet_left
        assert len(twin._facet_left) == 2

    def test_without_access_extension_a_cascade_end_section_raises(self, mmi_twin_gds):
        gds = GdsGeometry(mmi_twin_gds, core_layers=[LAYER_SI],
                          region_lengths=TWIN_REGIONS)  # access_extension=False
        b = self._bounds()
        with pytest.raises(ValueError, match="access_extension"):
            gds.section_gds(b[0], b[1], 6.0, stub_left=self.STUB, input_idx=1)

    def test_without_access_extension_an_interior_section_still_works(self, mmi_twin_gds):
        gds = GdsGeometry(mmi_twin_gds, core_layers=[LAYER_SI],
                          region_lengths=TWIN_REGIONS)
        b = self._bounds()
        g = gds.section_gds(b[1], b[2], 6.0)
        assert _polygon_area(g[LAYER_SI]) == pytest.approx(2.8 * (b[2] - b[1]))


class TestCutValidity:
    """The measured stand-in for MMI._check_cut's parametric band lookup."""

    @pytest.mark.parametrize("x", [-18.0, -14.0, -13.75, 0.0, 13.75, 14.0, 18.0])
    def test_legal_cuts_agree_with_mmi(self, twin, twin_parametric, x):
        twin_parametric._check_cut(x)                       # must not raise
        twin._check_cut(x, {"input_idx": 1, "output_idx": 1})

    @pytest.mark.parametrize("x", [-17.0, -16.0, -15.0, 15.0, 16.0, 17.0])
    def test_mid_taper_cuts_raise_like_mmi(self, twin, twin_parametric, x):
        with pytest.raises(ValueError):
            twin_parametric._check_cut(x)
        with pytest.raises(ValueError, match="varying cross-section"):
            twin._check_cut(x, {"input_idx": 1, "output_idx": 1})

    def test_band_endpoint_is_legal_though_one_side_varies(self, twin):
        """The 'uniform on at least ONE side' rule, which is what keeps this
        from being stricter than the parametric geometries it stands in for:
        at the body face a taper varies on one side and the uniform body sits
        on the other, and MMI permits that cut."""
        twin._check_cut(-TWIN_HALF_BODY, {"input_idx": 1})
        twin._check_cut(+TWIN_HALF_BODY, {"output_idx": 1})

    def test_section_with_a_mid_taper_cut_raises(self, twin):
        with pytest.raises(ValueError, match="varying cross-section"):
            twin.section_gds(-16.0, 0.0, 6.0, input_idx=1)


class TestToGdsEmeExtension:
    def test_straight_length_is_honoured_only_with_access_extension(
            self, mmi_twin_gds, twin):
        """The historical contract (ignore straight_length) is the default;
        access_extension=True opts into the parametric geometries' behaviour,
        which is what plot_hybrid_geometry's own
        `total_length/2 + straight_length` frame shift assumes."""
        plain = GdsGeometry(mmi_twin_gds, core_layers=[LAYER_SI],
                            region_lengths=TWIN_REGIONS)
        a = plain.to_gds_eme(straight_length=0.0, oxide_width=6.0)
        b = plain.to_gds_eme(straight_length=3.0, oxide_width=6.0)
        assert a[LAYER_SI] == b[LAYER_SI]

        ext = twin.to_gds_eme(straight_length=3.0, oxide_width=6.0)
        xs = [x for poly in ext[LAYER_SI] for x, _y in poly]
        assert min(xs) == pytest.approx(0.0, abs=1e-9)
        assert max(xs) == pytest.approx(TWIN_SPAN + 6.0, abs=1e-9)
        ox_xs = [x for x, _y in ext[(1, 0)][0]]
        assert max(ox_xs) - min(ox_xs) == pytest.approx(TWIN_SPAN + 6.0)
