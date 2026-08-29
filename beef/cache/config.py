"""build_run_config(): fingerprint one solve_section()-style call into a
JSON-safe dict, for the persistent save-folder cache (see beef/cache/__init__.py).
"""
from __future__ import annotations

import inspect
import logging
from typing import Any, Dict, Tuple

import numpy as np

from ..materials.stack import MaterialStack
from .._log import get_logger

logger = get_logger(__name__)

_PRECISION = 9

# Bump this whenever solver-internals derivation logic changes in a way that
# could change a simulated result without changing any tracked config field
# (e.g. a bug fix to how a default monitor/source position is computed from
# already-tracked ingredients) — every existing cache entry's hash changes,
# so old entries are simply never matched again (not read, not deleted) and
# get transparently resimulated/re-saved under the new hash on next use. Not
# tied to beef.__version__ (currently static, see High-level recommendation
# #3 in notes.md) — a plain hand-bumped counter is the cheaper, more robust
# half of Solver recommendation #4's two complementary fixes; the other half
# (fingerprinting fully-resolved positions, not just their ingredients) is
# handled per-method where it applies — see FDTDSolver.solve_section()'s use
# of SectionFDTDLayout in its call_kwargs.
_CACHE_SCHEMA_VERSION = 1

# Per (solver_type, method_name) allowlist of solver attributes that method
# actually reads (including transitively, through helper properties/methods
# it calls) — traced by hand against the current beef/solvers/*.py source.
# Scoping _solver_config() to this instead of dumping all of vars(solver)
# closes the opposite-direction gap from the schema-version bump above: two
# solver instances differing only in an attribute a given method never
# touches (e.g. FDTDSolver.num_input_modes, which solve_section() never
# reads) would otherwise get two disjoint cache entries for what is actually
# the same physical result. A (solver_type, method_name) pair with no entry
# here falls back to the old "dump everything" behavior (see
# _solver_config()) — under-including an attribute that DOES matter would
# risk a silent false-positive cache hit, which is strictly worse than the
# spurious-miss cost of over-including, so "unknown pair -> include
# everything" is the safe default direction for this fallback.
_METHOD_ATTRS: Dict[Tuple[str, str], Tuple[str, ...]] = {
    ("FDTDSolver", "solve"): (
        "res", "sim_width", "sim_height", "three_d", "dpml", "source_margin",
        "monitor_margin", "num_output_modes", "output_monitor_centers",
        "output_monitor_widths", "material_stack", "waveguide_length",
        "decay_tol", "minimum_run_time", "decay_method", "energy_decay_dt",
        "source_bandwidth_factor",
        "eps_averaging", "subpixel_tol", "subpixel_maxeval",
    ),
    ("FDTDSolver", "solve_section"): (
        "res", "sim_width", "sim_height", "three_d", "dpml", "source_margin",
        "monitor_margin", "half_domain_aperture", "material_stack",
        "decay_tol", "minimum_run_time", "decay_method", "energy_decay_dt",
        "source_bandwidth_factor",
        "eps_averaging", "subpixel_tol", "subpixel_maxeval",
    ),
    # section_mode_neff() runs MPB on an initialized grid with no time-stepping
    # and no source, so decay_tol/source_margin/source_bandwidth_factor cannot
    # reach it; monitor_margin cannot either (its monitor plane comes from
    # monitor_center/monitor_size, both fingerprinted as call kwargs).
    ("FDTDSolver", "section_mode_neff"): (
        "res", "sim_width", "sim_height", "three_d", "dpml", "material_stack",
        "eps_averaging", "subpixel_tol", "subpixel_maxeval",
    ),
    # Same attribute set as solve_section(): the merged method builds its
    # domain through the very same section_fdtd_layout() call, differing only
    # in how many monitors it hangs inside it.
    ("FDTDSolver", "solve_section_multi_length"): (
        "res", "sim_width", "sim_height", "three_d", "dpml", "source_margin",
        "monitor_margin", "half_domain_aperture", "material_stack",
        "decay_tol", "minimum_run_time", "decay_method", "energy_decay_dt",
        "source_bandwidth_factor",
        "eps_averaging", "subpixel_tol", "subpixel_maxeval",
    ),
    # num_cells added 2026-08-20 (mmi_1x4_tapered_body_validation.md step 3):
    # solve_section()'s real-solve branch now reads self.num_cells[0] to build
    # a genuine multi-cell EME body instead of always Ls=[length]. Every
    # existing entry predates this key and so misses once (same accepted cost
    # as the minimum_run_time/decay_method precedent above) rather than
    # risking a single-cell result silently serving a multi-cell request.
    #
    # num_modes is deliberately NOT here: solve_section() still takes its mode
    # count from n_modes_in/n_modes_out (call kwargs, already fingerprinted),
    # never from self.num_modes -- see that method's own comment on why
    # reading it would silently change modal_multiplexer_hybrid.py. Adding it
    # here would additionally strand every existing entry for no behavioural
    # reason, splitting the cache on an attribute the method does not read.
    ("EMESolver", "solve_section"): (
        "material_stack", "oxide_width", "oxide_height", "res", "temperature",
        "mode_solver", "num_cells",
    ),
    ("EMESolver", "solve_section_length_sweep"): (
        "material_stack", "oxide_width", "oxide_height", "res", "temperature", "mode_solver",
    ),
    ("EMESolver", "solve"): (
        "num_modes", "num_cells", "lengths", "straight_length", "oxide_width",
        "oxide_height", "res", "material_stack", "temperature", "mode_solver",
        "region_cell_lengths",
    ),
    ("EMESolver", "solve_length_sweep"): (
        "num_modes", "num_cells", "lengths", "straight_length", "oxide_width",
        "oxide_height", "res", "material_stack", "temperature", "mode_solver",
        "region_cell_lengths",
    ),
}


# Attributes omitted from the fingerprint **when they hold their constructor
# default**, and included as normal when they do not.
#
# WHY THIS EXISTS, from the failure that produced it. Adding an attribute to a
# solver's __init__ silently invalidates every cache entry written before it,
# because the entry's stored config simply has no such key and so cannot hash
# equal. That is not hypothetical: `FDTDSolver.minimum_run_time` (added
# 2026-08-08 with the modal multiplexer's early-stop fix) stranded all 48
# `solve_reflection` entries in `example_sim_cache/`, including the HCG
# mirror's own res-60 3-D production run — a ~55 min result that was still
# perfectly valid, since the runs it described had no floor, which is exactly
# what the new default (0.0) means. Found 2026-08-18 while making
# `hcg_fabry_perot_hybrid.py` solve its own pieces; confirmed by deleting the
# one attribute off a freshly built solver and watching the same call hit.
#
# WHY NOT JUST "SKIP EVERY DEFAULT-VALUED ATTRIBUTE". Because most tracked
# attributes ARE at their defaults in real runs — `material_stack`, `dpml`,
# `decay_tol`, `sim_height`, `eps_averaging`, the subpixel pair — and every
# existing entry records them. Dropping those from the key would strand the
# whole cache, i.e. cause precisely the damage this is meant to undo. The set
# is therefore explicit and additive: an attribute belongs here only because
# it postdates entries that are still wanted.
#
# THE OTHER HALF OF THE FIX. Omission alone is not enough — a method whose
# _METHOD_ATTRS allowlist never listed the attribute could not distinguish a
# NON-default value either, which is the strictly worse direction (a silent
# false-positive hit, see _METHOD_ATTRS above). `solve`, `solve_section` and
# `solve_section_multi_length` all pass `minimum_run_time` to
# `stop_when_dft_decayed()` but omitted it from their allowlists, so
# `modal_multiplexer_hybrid.py` (100.0) and `modal_multiplexer_fdtd.py`
# (1000.0) could have reused a no-floor result — defeating the very fix those
# values were added for. It is listed in all three now; the default-omission
# below is what keeps their existing 0.0 entries matching regardless.
# decay_method/energy_decay_dt (added 2026-08-18, same reasoning as
# minimum_run_time above): a "dft"-default run must keep matching every
# entry written before these attributes existed. A non-default
# decay_method="energy" run is a genuinely different simulation (a
# different meep stop condition — see FDTDSolver's own docstring) and
# must NOT be served by, or serve, a "dft" entry, which is exactly what
# listing both in _METHOD_ATTRS above (not just omitting the default
# here) guarantees.
# region_cell_lengths (added 2026-08-29, same reasoning as the three above):
# it defaults to None, meaning "divide every region evenly", which is exactly
# what every entry written before it existed did. Omitting the default keeps
# all of those valid. A non-default value places cell boundaries somewhere an
# even division cannot reach (a waisted MMI body's five segments), which is a
# genuinely different discretisation of the same device and must NOT be served
# by, or serve, an evenly-divided entry -- which is what listing it in
# _METHOD_ATTRS above guarantees.
_DEFAULT_OMITTED_ATTRS = frozenset({"minimum_run_time", "decay_method", "energy_decay_dt",
                                    "region_cell_lengths"})


def _constructor_defaults(solver_cls) -> Dict[str, Any]:
    """`{parameter: default}` for one solver's ``__init__``, cached per class.

    Read off the signature rather than hardcoded, so the omission rule cannot
    drift from the constructor the way a duplicated constant would.
    """
    cached = getattr(solver_cls, "_beef_ctor_defaults", None)
    if cached is None:
        cached = {
            name: p.default
            for name, p in inspect.signature(solver_cls.__init__).parameters.items()
            if p.default is not inspect.Parameter.empty
        }
        try:
            solver_cls._beef_ctor_defaults = cached
        except (AttributeError, TypeError):      # e.g. a slotted/frozen class
            pass
    return cached


def _is_omitted_default(solver, key: str, value: Any) -> bool:
    """True when `key` is default-omitted AND `value` is that default.

    `type(...) is type(...)` guards the one way `==` lies here: `False == 0`
    and `0 == 0.0`, so a bool attribute defaulting to False would otherwise
    also swallow a real 0, and vice versa.
    """
    if key not in _DEFAULT_OMITTED_ATTRS:
        return False
    defaults = _constructor_defaults(type(solver))
    if key not in defaults:
        return False
    default = defaults[key]
    return type(value) is type(default) and value == default


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


def _solver_config(solver, method_name: str) -> Dict[str, Any]:
    """vars(solver) filtered to JSON-safe scalars/lists/dicts, and scoped to
    the attributes _METHOD_ATTRS says method_name actually reads (falling
    back to every attribute when (solver_type, method_name) has no entry —
    see _METHOD_ATTRS' module-level docstring for why that's the safe
    fallback direction).

    A non-scalar attribute (e.g. a material stack) is serialized via its own
    cache_key() method if it has one — the general hook any custom material
    stack or other solver-attribute object should implement to fingerprint
    itself (see MaterialStack's special-case below for the built-in stacks).
    An object with neither cache_key() nor a JSON-safe form is skipped with a
    warning rather than silently dropped, since an omitted attribute that
    actually affects the physical result would produce a silent false-
    positive cache hit — see beef/cache/__init__.py.

    An attribute in `_DEFAULT_OMITTED_ATTRS` sitting at its constructor
    default is left out entirely, so that adding it to a solver does not
    invalidate entries written before it existed — see that set's own
    comment for the stranding it was written to undo.
    """
    allowlist = _METHOD_ATTRS.get((type(solver).__name__, method_name))
    out: Dict[str, Any] = {}
    for key, value in vars(solver).items():
        if allowlist is not None and key not in allowlist:
            continue
        if _is_omitted_default(solver, key, value):
            continue
        cache_key = getattr(value, "cache_key", None)
        if callable(cache_key):
            out[key] = _json_safe(cache_key())
            continue
        if isinstance(value, MaterialStack):
            out[key] = {"wg_height": value.wg_height, "sz": value.sz}
            continue
        try:
            out[key] = _json_safe(value)
        except TypeError:
            logger.warning(
                "cache: solver attribute %r (%s) is not JSON-safe and has no "
                "cache_key() method, so it was excluded from the cache config "
                "— if it affects the physical result, cache hits involving "
                "this solver may be false positives. Give it a cache_key() "
                "returning a JSON-safe fingerprint to fix this.",
                key, type(value).__name__,
            )
    return out


def _fingerprint(
    solver, method_name: str, call_kwargs: Dict[str, Any],
    geometry_class_name: str, geometry_polygons: Dict[Any, Any],
) -> Dict[str, Any]:
    """Assemble the JSON-safe fingerprint dict shared by both config builders."""
    return {
        "cache_schema_version": _CACHE_SCHEMA_VERSION,
        "solver_type": type(solver).__name__,
        "solver_config": _solver_config(solver, method_name),
        "method": method_name,
        "call_kwargs": _json_safe(call_kwargs),
        "geometry_class": geometry_class_name,
        "geometry_polygons": {
            f"{layer[0]},{layer[1]}": _json_safe(polys)
            for layer, polys in sorted(geometry_polygons.items())
        },
    }


def build_run_config(solver, method_name: str, call_kwargs: Dict[str, Any], geometry) -> Dict[str, Any]:
    """Fingerprint one solve_section()/solve_section_length_sweep() call.

    Args:
        solver: FDTDSolver or EMESolver instance.
        method_name: "solve_section" or "solve_section_length_sweep".
        call_kwargs: every argument the method call received that can affect
            the physical result (x0, x1, wavelengths, n_modes_in/out,
            stub_left/right, monitor_center/size, lengths, geometry_kwargs,
            ...) — save_dir itself must not be included. If any of those
            arguments are optional/defaulted and get *resolved* into derived
            absolute positions before the simulation runs (e.g.
            FDTDSolver.solve_section()'s monitor_center/monitor_size ->
            SectionFDTDLayout), include the resolved values too (as a
            separate key, not a replacement) — otherwise a future change to
            that derivation logic leaves an old config.json looking
            identical to a new call's fingerprint despite simulating
            different physical positions (Solver recommendation #4).
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
    return _fingerprint(solver, method_name, call_kwargs, type(geometry).__name__, section)


def build_whole_device_config(
    solver, method_name: str, call_kwargs: Dict[str, Any], geometry, geometry_polygons,
) -> Dict[str, Any]:
    """Fingerprint one whole-device solve() call (FDTDSolver/EMESolver).

    Same shape as build_run_config(), but for the plain (non-section)
    solve() methods, which have no x0/x1 to hand to geometry.section_gds().
    Callers must instead pass the exact polygon dict their own solve() body
    already builds and hands to the underlying meep/meow simulation —
    geometry.polygons(stub_length=...) for FDTDSolver.solve() (the same
    builder to_meep_prisms() calls), geometry.to_gds_eme(...) for
    EMESolver.solve() — so the fingerprint matches precisely what gets
    simulated, with no meep/meow objects involved.

    Args:
        solver: FDTDSolver or EMESolver instance.
        method_name: "solve".
        call_kwargs: every argument the method call received that can affect
            the physical result — save_dir itself must not be included.
        geometry: the DeviceGeometry instance being solved.
        geometry_polygons: Dict[layer_tuple, List[polygon]], as already built
            by the caller (geometry.polygons() or geometry.to_gds_eme()).

    Returns:
        A JSON-safe dict, same shape as build_run_config()'s return value.
    """
    return _fingerprint(solver, method_name, call_kwargs, type(geometry).__name__, geometry_polygons)
