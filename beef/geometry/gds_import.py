"""GDS file import: load raw polygons from an on-disk .gds layout and wrap
them as a DeviceGeometry for EME.

BEEF's other geometries (MMI, Taper, ...) build their polygons themselves,
parametrically, in Python. This module instead reads polygons already drawn
by an external tool (gdspy, gdsfactory, ...) and adapts them to the same
Dict[(layer, datatype), List[polygon]] convention every other geometry uses,
so they flow through extrude_polygons()/EMESolver unchanged.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from .base import DeviceGeometry, _oxide_box

_OX_LAYER = (1, 0)


def load_gds_polygons(
    gds_path: str,
    layers: Optional[Sequence[Tuple[int, int]]] = None,
    cell_name: Optional[str] = None,
) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
    """Read a GDS file and return its polygons as BEEF's standard
    Dict[(layer, datatype), List[polygon]] dict (polygon = list of (x, y)
    vertex tuples, µm — gdstk reports coordinates in µm regardless of the
    file's own unit/precision, so no manual unit conversion is needed).

    Args:
        gds_path: path to the .gds file.
        layers: if given, only these (layer, datatype) tuples are extracted
            (everything else in the file is ignored). Extracts every layer
            present otherwise.
        cell_name: name of the cell to read; defaults to the library's sole
            top-level cell (raises if there is more than one and no name is
            given).

    Returns:
        Dict[(layer, datatype), List[List[(x, y)]]], vertices in the GDS
        file's own coordinate frame (no shift/centering applied).
    """
    import gdstk

    lib = gdstk.read_gds(str(gds_path))

    if cell_name is not None:
        matches = [c for c in lib.cells if c.name == cell_name]
        if not matches:
            raise ValueError(f"cell {cell_name!r} not found in {gds_path}")
        cell = matches[0]
    else:
        top = lib.top_level()
        if len(top) != 1:
            names = [c.name for c in top]
            raise ValueError(
                f"{gds_path} has {len(top)} top-level cells {names} — "
                "pass cell_name= to disambiguate"
            )
        cell = top[0]

    polys = cell.get_polygons(apply_repetitions=True, include_paths=True, depth=None)

    out: Dict[Tuple[int, int], List[List[Tuple[float, float]]]] = {}
    for p in polys:
        layer = (p.layer, p.datatype)
        if layers is not None and layer not in layers:
            continue
        out.setdefault(layer, []).append([(float(x), float(y)) for x, y in p.points])
    return out


class GdsGeometry(DeviceGeometry):
    """DeviceGeometry backed by polygons imported directly from a GDS file.

    A straight-through (single path, no branching) device: the imported
    polygons already contain whatever access waveguides the GDS layout
    drew, so unlike Taper/MMI/SBend (which synthesize their own stub
    extensions from the `straight_length` argument on demand),
    to_gds_eme() places the imported polygons as-is and appends nothing.

    `region_lengths` is a caller-supplied breakdown of the imported
    geometry's own total span (e.g. [stub, taper, stub]), used only as the
    return value of eme_lengths() — a convenience for building the
    `lengths=` list passed to EMESolver. It must sum to the bounding-box
    span of `core_layers`' polygons (checked at construction).
    """

    def __init__(
        self,
        gds_path: str,
        core_layers: Sequence[Tuple[int, int]],
        region_lengths: Sequence[float],
        cell_name: Optional[str] = None,
    ):
        self.gds_path = str(gds_path)
        self.core_layers = list(core_layers)
        self.region_lengths = list(region_lengths)

        raw = load_gds_polygons(self.gds_path, layers=self.core_layers, cell_name=cell_name)

        xs = [x for layer in self.core_layers for poly in raw.get(layer, []) for x, _y in poly]
        if not xs:
            raise ValueError(
                f"no polygons found on layers {self.core_layers} in {gds_path}"
            )
        x_min, x_max = min(xs), max(xs)
        span = x_max - x_min
        lengths_sum = sum(self.region_lengths)
        if abs(span - lengths_sum) > 1e-6:
            raise ValueError(
                f"region_lengths sum to {lengths_sum:.4f} µm but the imported "
                f"geometry spans {span:.4f} µm on layers {self.core_layers} "
                f"({gds_path})"
            )

        self._span = span
        x_off = -(x_min + x_max) / 2.0  # center the bbox at x=0 (device frame)
        self._polys = {
            layer: [[(x + x_off, y) for x, y in poly] for poly in polys]
            for layer, polys in raw.items()
        }

    @property
    def total_length(self) -> float:
        return self._span

    def eme_lengths(self, straight_length: float) -> List[float]:
        return list(self.region_lengths)

    def to_gds_eme(
        self, straight_length: float, oxide_width: float,
        x_start: float = 0.0, **kwargs
    ) -> Dict[Tuple[int, int], List[List[Tuple[float, float]]]]:
        x_off = x_start + self._span / 2.0
        out = {
            layer: [[(x + x_off, y) for x, y in poly] for poly in polys]
            for layer, polys in self._polys.items()
        }
        out[_OX_LAYER] = [_oxide_box(x_start, x_start + self._span, oxide_width)]
        return out
