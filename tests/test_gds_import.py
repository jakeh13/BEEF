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
