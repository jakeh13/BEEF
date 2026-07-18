"""build_run_config(): fingerprint one solve_section()-style call into a
JSON-safe dict, for the persistent save-folder cache (see feces/cache/__init__.py).
"""
from __future__ import annotations

import logging
from typing import Any, Dict

import numpy as np

from ..materials.stack import MaterialStack

logger = logging.getLogger(__name__)

_PRECISION = 9


def _json_safe(obj: Any, precision: int = _PRECISION) -> Any:
    if obj is None or isinstance(obj, (bool, str, int)):
        return obj
    if isinstance(obj, float):
        return round(obj, precision)
    if isinstance(obj, np.generic):
        return _json_safe(obj.item(), precision)
    if isinstance(obj, np.ndarray):
        return [_json_safe(v, precision) for v in obj.tolist()]
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v, precision) for v in obj]
    if isinstance(obj, dict):
        return {str(k): _json_safe(v, precision) for k, v in sorted(obj.items(), key=lambda kv: str(kv[0]))}
    raise TypeError(f"cannot serialize {type(obj)!r} into a cache config")


def _solver_config(solver) -> Dict[str, Any]:
    """vars(solver) filtered to JSON-safe scalars/lists/dicts.

    material_stack is an object (not a scalar), so it's serialized
    explicitly; any other non-scalar attribute a solver subclass might carry
    is skipped with a warning rather than silently dropped, since an omitted
    attribute that actually affects the physical result would produce a
    silent false-positive cache hit — see feces/cache/__init__.py.
    """
    out: Dict[str, Any] = {}
    for key, value in vars(solver).items():
        if isinstance(value, MaterialStack):
            out[key] = {"wg_height": value.wg_height, "sz": value.sz}
            continue
        try:
            out[key] = _json_safe(value)
        except TypeError:
            logger.warning(
                "cache: solver attribute %r (%s) is not JSON-safe and was "
                "excluded from the cache config — if it affects the "
                "physical result, cache hits involving this solver may be "
                "false positives.",
                key, type(value).__name__,
            )
    return out


def build_run_config(solver, method_name: str, call_kwargs: Dict[str, Any], geometry) -> Dict[str, Any]:
    """Fingerprint one solve_section()/solve_section_length_sweep() call.

    Args:
        solver: FDTDSolver or EMESolver instance.
        method_name: "solve_section" or "solve_section_length_sweep".
        call_kwargs: every argument the method call received that can affect
            the physical result (x0, x1, wavelengths, n_modes_in/out,
            stub_left/right, monitor_center/size, lengths, geometry_kwargs,
            ...) — save_dir itself must not be included.
        geometry: the DeviceGeometry instance being solved.

    Returns:
        A JSON-safe dict: solver config + call kwargs + the exact section
        geometry (via geometry.section_gds(), not the whole-device
        polygons() — this captures precisely the shape solve_section() drew
        for x0/x1/stub_left/stub_right, with no meep/meow objects involved).
    """
    geometry_kwargs = call_kwargs.get("geometry_kwargs") or {}
    x0 = call_kwargs["x0"]
    x1 = call_kwargs["x1"]
    stub_left = call_kwargs.get("stub_left", 0.0)
    stub_right = call_kwargs.get("stub_right", 0.0)
    oxide_width = getattr(solver, "oxide_width", None)
    if oxide_width is None:
        oxide_width = getattr(solver, "sim_width", 5.0)

    section = geometry.section_gds(
        x0, x1, oxide_width, stub_left=stub_left, stub_right=stub_right, **geometry_kwargs
    )
    geometry_polygons = {
        f"{layer[0]},{layer[1]}": _json_safe(polys)
        for layer, polys in sorted(section.items())
    }

    return {
        "solver_type": type(solver).__name__,
        "solver_config": _solver_config(solver),
        "method": method_name,
        "call_kwargs": _json_safe(call_kwargs),
        "geometry_class": type(geometry).__name__,
        "geometry_polygons": geometry_polygons,
    }
