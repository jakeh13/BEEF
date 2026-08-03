"""UniformCrossSection — the shared cross-section object the ring's EME half
lacked (ring_resonator_validation.md checks 1 and 2).

The load-bearing test here is
test_reproduces_the_ring_scripts_own_cross_sections: this class only unblocks
those checks if it draws EXACTLY what the ring scripts build by hand from raw
polygons. If it draws something merely similar, an "engine agreement" comparison
built on it would be comparing two different structures and would be worse than
no check at all.
"""
import numpy as np
import pytest

from beef.geometry import UniformCrossSection

_WG, _OX = (2, 0), (1, 0)

# The ring's own constants (examples/ring_resonator/python/
# ring_resonator_hybrid_cascaded.py): W, GAP -> Y_RING = W + GAP.
W, GAP = 0.45, 0.25
Y_RING = W + GAP
MESH_WIDTH = 4.0
_MODE_MARGIN = MESH_WIDTH / 2 - W / 2
BODY_Y_CENTER = Y_RING / 2
BODY_MESH_WIDTH = Y_RING + W + 2 * _MODE_MARGIN


def _yspans(polys):
    """(y_lo, y_hi) per polygon, rounded to kill float noise."""
    return sorted((round(min(y for _, y in p), 9), round(max(y for _, y in p), 9))
                  for p in polys)


# ---------------------------------------------------------------------------
# The one that decides whether this class is fit for purpose
# ---------------------------------------------------------------------------
def test_reproduces_the_ring_scripts_own_cross_sections():
    """Must match ring_resonator_hybrid_cascaded.py's build_strip_cross_section()
    and build_body_cross_section() exactly, in the same layer layout."""
    length = 1.0

    # --- bare strip: rect_poly(0, -W/2, length, W/2) on layer 2,0;
    #     oxide rect_poly(0, -MESH_WIDTH/2, length, MESH_WIDTH/2) on 1,0
    strip = UniformCrossSection([(0.0, W)], length=length)
    got = strip.to_gds_eme(straight_length=0.0, oxide_width=MESH_WIDTH, x_start=0.0)
    assert _yspans(got[_WG]) == [(-W / 2, W / 2)]
    assert _yspans(got[_OX]) == [(-MESH_WIDTH / 2, MESH_WIDTH / 2)]
    xs = [x for p in got[_WG] for x, _ in p]
    assert (min(xs), max(xs)) == (0.0, length)

    # --- coupled pair: bus at y=0 and ring at y=Y_RING, both width W; oxide
    #     centred on BODY_Y_CENTER with span BODY_MESH_WIDTH
    pair = UniformCrossSection([(0.0, W), (Y_RING, W)], length=length)
    got = pair.to_gds_eme(straight_length=0.0, oxide_width=BODY_MESH_WIDTH, x_start=0.0)
    assert _yspans(got[_WG]) == [
        (round(-W / 2, 9), round(W / 2, 9)),
        (round(Y_RING - W / 2, 9), round(Y_RING + W / 2, 9)),
    ]
    # The ring centres its oxide on BODY_Y_CENTER; _oxide_box centres on 0, so
    # the shared object's oxide must be shifted by the caller. Assert the SPAN
    # matches and record the offset requirement rather than pretending it does
    # not exist.
    (ox_lo, ox_hi), = _yspans(got[_OX])
    assert round(ox_hi - ox_lo, 9) == round(BODY_MESH_WIDTH, 9)
    assert round((ox_lo + ox_hi) / 2, 9) == 0.0, (
        "to_gds_eme() centres the oxide on y=0; a caller wanting the ring's "
        "BODY_Y_CENTER-centred box must shift it (or pass strips shifted so the "
        "pair straddles y=0). Documented, not silently absorbed."
    )

    # The gap between the two cores must be the ring's GAP exactly.
    spans = _yspans(got[_WG])
    assert round(spans[1][0] - spans[0][1], 9) == round(GAP, 9)


# ---------------------------------------------------------------------------
# Uniformity: every cut legal, nothing varies with x
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("x", [-0.49, -0.1, 0.0, 0.1, 0.49])   # interior cuts only
def test_every_cut_is_legal_and_nothing_varies_with_x(x):
    """The defining property: no _check_cut, and width/centre are x-independent.
    Contrast Taper (raises inside), SBend / DirectionalCoupler (no interior bend
    cut), RingResonatorCouplerSection (no cut inside the arc)."""
    g = UniformCrossSection([(0.0, W), (Y_RING, W)], length=1.0)
    left = g.section_gds(-0.5, x, oxide_width=MESH_WIDTH)
    right = g.section_gds(x, 0.5, oxide_width=MESH_WIDTH)
    assert len(left[_WG]) == len(right[_WG]) == 2
    for idx in (0, 1):
        assert g.width_at(x, strip_idx=idx) == g.width_at(-0.4, strip_idx=idx)
        assert g.center_at(x, strip_idx=idx) == g.center_at(0.3, strip_idx=idx)


def test_a_zero_length_section_draws_nothing():
    """Cutting exactly at a boundary gives an empty section rather than a
    degenerate zero-area polygon. Pinned because the cut-legality test above
    deliberately stays off the endpoints, so this is where that edge lives."""
    g = UniformCrossSection([(0.0, W)], length=1.0)
    assert g.section_gds(-0.5, -0.5, oxide_width=MESH_WIDTH)[_WG] == []
    assert g.section_gds(0.5, 0.5, oxide_width=MESH_WIDTH)[_WG] == []
    # ...but a stub on either side revives it, since the span is then nonzero.
    assert len(g.section_gds(-0.5, -0.5, oxide_width=MESH_WIDTH,
                             stub_right=0.2)[_WG]) == 1


def test_section_lengths_add_up_across_an_arbitrary_cut():
    """Two abutting sections must together span exactly the original — the
    property a cut-position invariance sweep rests on."""
    g = UniformCrossSection([(0.0, W)], length=2.0)
    for cut in (-0.7, 0.0, 0.31, 0.9):
        a = g.section_gds(-1.0, cut, oxide_width=MESH_WIDTH)
        b = g.section_gds(cut, 1.0, oxide_width=MESH_WIDTH)
        la = max(x for p in a[_WG] for x, _ in p)
        lb = max(x for p in b[_WG] for x, _ in p)
        assert round(la + lb, 9) == 2.0


def test_stubs_extend_the_section_without_changing_the_cross_section():
    g = UniformCrossSection([(0.0, W), (Y_RING, W)], length=1.0)
    plain = g.section_gds(-0.2, 0.2, oxide_width=MESH_WIDTH)
    padded = g.section_gds(-0.2, 0.2, oxide_width=MESH_WIDTH,
                           stub_left=0.3, stub_right=0.4)
    assert _yspans(plain[_WG]) == _yspans(padded[_WG])
    assert round(max(x for p in padded[_WG] for x, _ in p), 9) == round(0.4 + 0.3 + 0.4, 9)


# ---------------------------------------------------------------------------
# source_aperture_at: the coupled-pair launch window
# ---------------------------------------------------------------------------
def test_aperture_is_none_for_a_single_strip():
    g = UniformCrossSection([(0.0, W)], length=1.0)
    assert g.source_aperture_at(0.0) is None


def test_aperture_excludes_the_neighbouring_core():
    """Halfway to the neighbour's centreline: for equal widths that is the gap
    midpoint, so the window must not reach the neighbour's near edge."""
    g = UniformCrossSection([(0.0, W), (Y_RING, W)], length=1.0)
    for idx in (0, 1):
        span = g.source_aperture_at(0.0, strip_idx=idx)
        yc = g.center_at(0.0, strip_idx=idx)
        near, far = yc - span / 2, yc + span / 2
        other = g.center_at(0.0, strip_idx=1 - idx)
        other_near_edge = other - np.sign(other - yc) * W / 2
        # The window edge on the neighbour's side must stop SHORT of that
        # neighbour's near edge -- the whole point of the aperture.
        window_edge = far if other > yc else near
        assert abs(window_edge - yc) < abs(other_near_edge - yc), (
            f"window edge {window_edge:.4f} reaches past the neighbour's near "
            f"edge {other_near_edge:.4f}")
        assert span == pytest.approx(Y_RING)


def test_aperture_uses_the_NEAREST_neighbour_when_there_are_several():
    g = UniformCrossSection([(0.0, 0.4), (1.0, 0.4), (5.0, 0.4)], length=1.0)
    assert g.source_aperture_at(0.0, strip_idx=0) == pytest.approx(1.0)
    assert g.source_aperture_at(0.0, strip_idx=2) == pytest.approx(4.0)


# ---------------------------------------------------------------------------
# Contract
# ---------------------------------------------------------------------------
def test_total_length_bbox_and_eme_lengths():
    g = UniformCrossSection([(0.0, W), (Y_RING, W)], length=1.5)
    assert g.total_length == 1.5
    length, yspan = g.bbox
    assert length == 1.5
    assert yspan == pytest.approx(Y_RING + W)
    assert g.eme_lengths(2.0) == [2.0, 1.5, 2.0]


def test_rejects_empty_strips_bad_width_and_bad_length():
    with pytest.raises(ValueError, match="at least one strip"):
        UniformCrossSection([], length=1.0)
    with pytest.raises(ValueError, match="widths must be"):
        UniformCrossSection([(0.0, 0.0)], length=1.0)
    with pytest.raises(ValueError, match="length must be"):
        UniformCrossSection([(0.0, W)], length=0.0)


def test_rejects_out_of_range_strip_idx():
    g = UniformCrossSection([(0.0, W)], length=1.0)
    with pytest.raises(ValueError, match="out of range"):
        g.width_at(0.0, strip_idx=3)
