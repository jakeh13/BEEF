"""Generate `mmi_2x2_twin.gds` — a GDS twin of the shipped 2×2 MMI geometry.

Step 1 of `validation_markdown/resolved/gds_geometry_hybrid_cascade_validation.md`:
the validation vehicle for `GdsGeometry` in a `HybridSolver` cascade. The
file's core-layer polygons are *geometrically identical* to what
``MMI(...)`` draws for the parameters `examples/mmi_2x2/python/
mmi_2x2_hybrid.py` uses, so running that example twice — once on `MMI`,
once on `GdsGeometry(mmi_2x2_twin.gds)` — is a controlled comparison in
which the geometry *backend* is the only thing that differs.

Generated FROM ``MMI.polygons()`` rather than drawn by hand, deliberately:
a hand-drawn twin would test whether two independently drawn layouts agree,
where what we actually want to isolate is the polygon *query* path
(`width_at`/`center_at`/`section_prisms`/`section_gds`).

``stub_length=0``: the file holds the bare device (input tapers + body +
output tapers, x ∈ [-18, +18] with the body centred at x=0), and NOT the
access waveguides a cascade's end sections need. Those are synthesized by
`GdsGeometry(..., access_extension=True)` from the cross-section measured
at each end facet — which is exactly the code path this vehicle exists to
exercise. Keeping them out of the file also makes the two geometries'
device frames coincide exactly: `GdsGeometry` centres the imported bbox on
x=0, `MMI`'s frame origin is its body centre, and for a symmetric 2×2 those
are the same point (`MMI.frame_offset == 0`), so `mmi_2x2_hybrid.py`'s
`boundaries` transfer to the GDS run unchanged.

Run:
    python /path/to/examples/gds_files/make_mmi_2x2_twin.py
"""

import os
import sys

import gdstk

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..", "..")))

from beef import MMI  # noqa: E402

WG_LAYER = (2, 0)
CELL_NAME = "mmi_2x2_twin"
OUT_PATH = os.path.join(_HERE, "mmi_2x2_twin.gds")

# Verbatim from examples/mmi_2x2/python/mmi_2x2_hybrid.py's MMI(...) call.
# Any edit here must be mirrored there (and vice versa) or the twin stops
# being a twin — mmi_2x2_hybrid_gds.py re-asserts the match at import time.
MMI_PARAMS = dict(
    n_inputs=2,
    n_outputs=2,
    mm_length=28.0,
    mm_width=2.8,
    waveguide_width=0.45,
    taper_length=4.0,
    taper_width=1.3,
    waveguide_offset=1.5,
    waveguide_height=0.22,  # SOI_220
)

# gdstk stores vertices as integers of `precision` metres. Every vertex the
# 2×2 produces (±0.225, ±0.65, ±0.75, ±1.4, ±14.0, ±18.0 µm) is an exact
# multiple of 1 nm, so the default 1e-9 precision round-trips losslessly —
# asserted below rather than assumed, since a silent 1 nm snap would show up
# later as an unexplained mismatch in the acceptance comparison.
UNIT = 1e-6
PRECISION = 1e-9


def build_twin_polygons():
    """The (2, 0) vertex lists the twin file holds, straight from MMI."""
    mmi = MMI(**MMI_PARAMS)
    return mmi, mmi.polygons(stub_length=0.0)[WG_LAYER]


def main():
    mmi, wg_polys = build_twin_polygons()

    lib = gdstk.Library(name="beef_validation", unit=UNIT, precision=PRECISION)
    cell = lib.new_cell(CELL_NAME)
    for verts in wg_polys:
        cell.add(gdstk.Polygon(verts, layer=WG_LAYER[0], datatype=WG_LAYER[1]))
    lib.write_gds(OUT_PATH)

    # Round-trip check: the file must give back the same vertices, in the
    # same order, to floating-point equality.
    from beef import load_gds_polygons

    back = load_gds_polygons(OUT_PATH, layers=[WG_LAYER])[WG_LAYER]
    if len(back) != len(wg_polys):
        raise AssertionError(f"round-trip lost polygons: {len(wg_polys)} -> {len(back)}")
    worst = 0.0
    for src, dst in zip(wg_polys, back):
        if len(src) != len(dst):
            raise AssertionError(f"round-trip changed vertex count: {len(src)} -> {len(dst)}")
        for (xs, ys), (xd, yd) in zip(src, dst):
            worst = max(worst, abs(xs - xd), abs(ys - yd))
    if worst > 1e-12:
        raise AssertionError(
            f"round-trip moved a vertex by {worst:.3g} µm — the 2×2's vertices "
            f"are no longer exact multiples of precision={PRECISION:g} m"
        )

    xs = [x for poly in back for x, _y in poly]
    ys = [y for poly in back for _x, y in poly]
    print(f"Wrote {OUT_PATH}")
    print(f"  cell            : {CELL_NAME!r}, layer {WG_LAYER}")
    print(f"  polygons        : {len(back)}")
    print(f"  x span          : [{min(xs):.3f}, {max(xs):.3f}] µm "
          f"(= {max(xs) - min(xs):.3f}; MMI.total_length = {mmi.total_length:.3f})")
    print(f"  y span          : [{min(ys):.3f}, {max(ys):.3f}] µm")
    print(f"  round-trip error: {worst:.3g} µm (exact)")


if __name__ == "__main__":
    main()
