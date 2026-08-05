"""Regression tests for solve_reflection()'s pre-flight layout guard
(validation_markdown/high_contrast_grating_validation.md step 10 item 6).

Locks down the two real, silent bugs this method has already shipped with
once each: source_margin/monitor_margin are measured from the PML inner
edge, not the geometry's own edge, so a source/monitor position that looks
fine in isolation can land in bare oxide once edge_clearance changes (the
"R=0.753, meaningless" bug); and a geometry that claims to honor
extend_input_waveguide can silently ignore stub_length for its input side,
leaving the exact cut facet the flag exists to remove (the HCG
waveguide-extension bug). No meep import needed — both checks are pure
point-in-polygon geometry.
"""
from __future__ import annotations

import pytest

from beef.solvers.fdtd import _point_in_polygon, check_reflection_layout

# A single 1x1 um square waveguide centered at y=0, spanning x in [-1, 1] --
# stands in for geometry.polygons(stub_length=...)[(2, 0)].
WG = [(-1.0, -0.5), (1.0, -0.5), (1.0, 0.5), (-1.0, 0.5)]


def test_point_in_polygon_inside_and_outside():
    assert _point_in_polygon(0.0, 0.0, WG)
    assert _point_in_polygon(-0.99, 0.4, WG)
    assert not _point_in_polygon(0.0, 0.6, WG)      # outside in y
    assert not _point_in_polygon(1.01, 0.0, WG)     # outside in x


def test_passes_when_source_and_monitor_are_inside_silicon():
    check_reflection_layout(
        [WG], src_x=-0.9, mon_x=-0.7, y0=0.0,
        extend_input_waveguide=False, domain_x_edge=-1.0,
    )  # must not raise


def test_raises_when_source_falls_outside_silicon():
    """The exact HCG bug: source_margin held fixed while edge_clearance grew
    moved the source past the real input stub's edge -- R=0.753, a large,
    physically meaningless deviation that was not caught until traced by
    hand."""
    with pytest.raises(ValueError, match="source.*outside the geometry's own silicon"):
        check_reflection_layout(
            [WG], src_x=-1.05, mon_x=-0.7, y0=0.0,
            extend_input_waveguide=False, domain_x_edge=-1.0,
        )


def test_raises_when_monitor_falls_outside_silicon():
    with pytest.raises(ValueError, match="monitor.*outside the geometry's own silicon"):
        check_reflection_layout(
            [WG], src_x=-0.9, mon_x=-1.2, y0=0.0,
            extend_input_waveguide=False, domain_x_edge=-1.0,
        )


def test_passes_when_extension_reaches_the_domain_boundary():
    extended = [(-2.0, -0.5), (1.0, -0.5), (1.0, 0.5), (-2.0, 0.5)]
    check_reflection_layout(
        [extended], src_x=-1.5, mon_x=-1.2, y0=0.0,
        extend_input_waveguide=True, domain_x_edge=-2.0,
    )  # must not raise


def test_raises_when_extend_input_waveguide_true_but_geometry_ignored_stub_length():
    """The exact HCG bug: a geometry that silently ignores stub_length for
    its input side would leave the flag's cut-facet artifact fully in
    place, just relocated -- must be caught, not passed through quietly."""
    with pytest.raises(ValueError, match="does not reach the domain boundary"):
        check_reflection_layout(
            [WG], src_x=-0.9, mon_x=-0.7, y0=0.0,
            extend_input_waveguide=True, domain_x_edge=-2.0,
        )


def test_extend_input_waveguide_false_does_not_check_domain_boundary():
    """False is the documented, still-supported default (facet artifact
    accepted, not fixed) -- the boundary-reach check must only fire when
    the caller actually asked for the extension."""
    check_reflection_layout(
        [WG], src_x=-0.9, mon_x=-0.7, y0=0.0,
        extend_input_waveguide=False, domain_x_edge=-2.0,
    )  # must not raise even though WG stops well short of domain_x_edge
