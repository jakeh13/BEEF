"""HybridSolver: generic N-section cascade combining any mix of solvers.

Any solver type — FDTDSolver, EMESolver, or a custom class — can occupy any
position in the chain, in any order, as long as it implements:

    solve_section(geometry, wavelengths, x0, x1,
                  n_modes_in, n_modes_out,
                  stub_left=0.0, stub_right=0.0,
                  **geometry_kwargs) -> complex ndarray [M, n_modes_out, n_modes_in]

Each section's matrix is the modal transfer function between the section's
left and right interfaces. The full-device response is the matrix-chain
product of the per-section matrices, contracted over the mode index shared
at each internal interface:

    alpha_full[m] = M_K[m] @ M_{K-1}[m] @ ... @ M_1[m]

This is the direct generalization of the slide-4 formula
    alpha_full[m] = sum_n alpha_in[m,n] * alpha_mm[m,n] * alpha_out[m,n]
to an arbitrary number of sections (>= 2) and arbitrary (possibly
non-diagonal) per-section transfer matrices — the sum over n above is exactly
what the matrix product contracts over.

Nothing in this module knows about any particular device type; the geometry
and its section boundaries are supplied by the caller.

Mirror-section reuse: for a left-right mirror-symmetric chain (e.g. a
splitter whose output taper is its input taper reversed), the output-side
section's transfer matrix is just the input-side section's transpose (mode
axes swapped) — the same reciprocity argument solve_section() already uses
internally, one level up (see FDTDSolver.solve_section()'s docstring).
`mirror_pairs` declares this at construction time, like meep's mp.Mirror is
declared once on mp.Simulation rather than re-specified per run — it's a
fact about the device's section topology, not something that changes call
to call. Unlike meep's symmetry (a structural property of one simulation),
BEEF's geometry_kwargs (input_idx/output_idx) can change which physical
port is actually live between solve() calls, so the geometric guardrail
that backs a mirror_pairs entry is re-checked every solve() call using that
call's actual kwargs, not just once at construction. If a device needs
different mirror configurations across different solve() calls (e.g. one
output port mirrors the input, another doesn't), construct two HybridSolver
instances sharing the same solvers/boundaries/interface_modes and a single
`cache` dict — cache keys are solver-identity/x-range/kwargs based, not tied
to the HybridSolver instance, so sections still dedupe across both.

Persistent save-folder cache: solve()/solve_length_sweep()'s `save_dir`
parameter is a separate, disk-persistent counterpart to `cache` (which only
dedupes in-memory, within/across calls in one process) — see beef/cache/.
Forwarded to each section solver's solve_section()/
solve_section_length_sweep() call, opt-in per solver (only solvers that
implement a save_dir parameter receive it; see _supports_save_dir()).
`cache_only` rides alongside it the same way (see _supports_cache_only()) —
demands a cache hit per section instead of falling back to simulating.
"""
from __future__ import annotations

import inspect
import logging
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np

from ..mode_identity import (
    best_correspondence, field_overlap_matrix, mesh_cell_centers, overlap_diagnostics,
    self_overlap_diagnostics, stack_meow_mode_fields,
)
from ..results.smatrix import SMatrix
from .._log import get_logger

logger = get_logger(__name__)

# Fallback diagnostic-mesh resolution (pixels/µm) for
# check_mode_correspondence, used only when neither adjoining solver exposes
# a `res` to derive one from. See _interface_diagnostic_mesh() for why this
# is a RESOLUTION and not a fixed pixel count.
_CORR_DEFAULT_RES = 40.0
# Half-width (µm) of the query window straddling an interior interface x —
# small enough to sit inside a straight/uniform stub, per
# CompoundGeometry.section_prisms()'s own cut-inside-a-taper failure mode.
_CORR_HALF_SPAN = 0.05


def _supports_save_dir(method) -> bool:
    return "save_dir" in inspect.signature(method).parameters


def _supports_cache_only(method) -> bool:
    return "cache_only" in inspect.signature(method).parameters


def _supports_neff_override(method) -> bool:
    return "neff_override" in inspect.signature(method).parameters


def _supports_pad(method) -> bool:
    return "pad" in inspect.signature(method).parameters


class HybridSolver:
    """Chains an arbitrary sequence of solvers across contiguous device sections.

    Args:
        solvers: ordered list of solver instances (length >= 2), left to
            right. Any mix of solver types/order is supported as long as
            each implements solve_section() (see module docstring).
        boundaries: list of len(solvers)+1 x-positions (device frame, µm)
            partitioning the device. boundaries[0] / boundaries[-1] are the
            true physical port locations; the len(solvers)-1 interior points
            are the cuts between consecutive solvers.
        interface_modes: list of len(solvers)-1 ints — number of modes
            carried at each internal interface between consecutive solvers.
        mirror_pairs: optional {k: j} map — section k's transfer matrix is
            reused as section j's transposed result instead of being solved
            (see module docstring). j must be < k (the partner must already
            be solved earlier in the left-to-right loop) and must not itself
            be a mirrored section (no chained/transitive mirroring).
        mirror_axis: x-position sections are reflected about, for the
            geometric guardrail. Defaults to the device midpoint,
            (boundaries[0] + boundaries[-1]) / 2.
        mirror_check_samples: number of x points sampled across each
            mirrored section's span when checking width_at/center_at
            symmetry (see _validate_mirror_pair).
        mirror_check_tol: tolerance (µm) for the width_at/center_at
            symmetry check — tight enough to catch a wrong input_idx/
            output_idx pairing (a port-pitch-scale mismatch) while
            tolerant of floating-point noise.
        calibrate_body_neff: propagate each body section on the β its
            neighbouring FDTD sections actually accumulate, instead of on its
            own mode solver's n_eff. **Default False** — this changes results
            and every recorded golden predates it.

            Why it exists (dc_validation.md item 6). The two engines are both
            MPB and agree on mode *index* to 0.56%, but MPB's index and the
            phase meep's time-stepper really picks up differ by ~2.3% per mode
            at res 40 — different discretizations of the same continuum
            problem. Both modes are biased the same way, so the differential is
            only 0.129 percentage points; but coupling depends on Δβ, a
            difference of two β near 9.4 rad/µm, and that ``β/Δβ ≈ 99`` lever
            turns it into a 10.4% error. Mesh refinement cannot reach it: ~res
            143, about 164× the 3-D cost.

            This is a correctness fix rather than a fudge. The hybrid's job is
            to reproduce a full-device FDTD run at lower cost, and that
            reference *is* meep's discretization — so propagating the body on
            MPB's β was the mismatch. On the coupler it removes the stub-sweep
            drift identically (the drift is exactly
            ``2s·(1 - Δβ_FDTD/Δβ_EME)``) and moves L₅₀ from 14.54 µm to
            ~15.02 µm against a 15.13 µm full-device reference.

            Costs `len(calibrate_lengths)` extra FDTD section solves per body,
            cached under `save_dir` like any other. Applies only to sections
            whose solver has no `section_mode_beta()` of its own — in practice
            the EME bodies; an all-FDTD chain is unaffected.
        calibrate_lengths: section lengths (µm) the calibration fits β over.
            Default `[span/2, span]` of the body being calibrated. Widening
            the spacing risks the 2π branch ambiguity
            `FDTDSolver.section_mode_beta()` guards against; it will raise
            rather than guess.
        check_mode_correspondence: verify, at each interior interface, that
            the adjoining solvers' mode bases actually name the same modes in
            the same order — the assumption solve()'s matrix-chain
            contraction otherwise makes with no runtime check at all (see
            module docstring). **Default False** — opt-in, since it costs two
            extra cheap eigenmode solves per interior interface per solve()/
            solve_length_sweep() call, on top of the section solves
            themselves. Compares mode *shape* (transverse E-field overlap,
            via beef.mode_identity — see that module's docstring) rather than
            n_eff, on a shared diagnostic mesh built per interface (see
            `mode_correspondence_res`). Only *guided* modes are compared:
            radiation modes (n_eff below the cladding index) have no
            engine-independent shape to match, and every interface basis in
            this repo carries at least one — see
            `_check_interface_correspondence()`. Silently skipped (with an info log) at
            any interface where either adjoining solver doesn't implement
            `section_mode_fields()` (FDTDSolver-style) or `section_modes()`
            (EMESolver-style), or exposes no domain size
            (`sim_width`/`sim_height` or `oxide_width`/`oxide_height`) — this
            is a best-effort guard, not a hard requirement on custom solver
            types. Also skipped (with a WARNING log, regardless of
            `mode_correspondence_action`) if the field-sampling query itself
            fails with a geometry `ValueError` — this is common in practice,
            not an edge case: `FDTDSolver.section_mode_fields()` always pads
            its own query window by `dpml` on both sides, so any interior
            interface whose straight run-out is thinner than that (verified:
            every `directional_coupler`/`mmi_*` hybrid example's
            `body_stub`/`BODY_STUB` is 0.25 µm against `dpml=0.5`) makes the
            check unable to run there at all; a `CompoundGeometry`-composed
            device can additionally trip its own known stub
            double-application defect the same way (verified: the modal
            multiplexer's real chain hits this). Neither is a correspondence
            finding — always treat the WARNING log as "could not check here,"
            distinct from "checked and found a mismatch." Pass
            `mode_correspondence_pad` (below) to unblock the `body_stub` <
            `dpml` case specifically. Confirmed the
            actual comparison passes cleanly on a real meep/meow FDTD<->EME
            interface where it CAN run (a uniform `WaveguideRun` section).
            Found two real bugs this way on the modal multiplexer
            (validation_markdown/modal_multiplexer_validation.md gaps 1/6)
            when run by hand via
            examples/verification/python/audit_interface_modes.py; this
            wires the same check into HybridSolver itself.
        mode_correspondence_action: "raise" (default) aborts solve() with a
            RuntimeError naming the offending interface; "warn" only logs.
        mode_correspondence_res: resolution (pixels/µm) of the shared
            diagnostic mesh both sides' modes are compared on. Default None
            derives it as the coarser of the two adjoining solvers' own
            `res`, so the check never demands a finer mesh than the device is
            actually simulated at. Raise it if a check reports a mismatch you
            suspect is under-resolution rather than a real basis
            disagreement — an EME side *solves* its modes on this mesh, so
            too coarse a value makes one side's modes physically wrong rather
            than merely coarsely sampled (see
            `_interface_diagnostic_res()` for the measured failure).
        mode_correspondence_min_overlap: minimum per-mode transverse-E-field
            overlap (in [0, 1]) required on the diagonal of the best
            correspondence found; below this, or a non-identity
            correspondence, is treated as a failure.
        mode_correspondence_pad: straight-material clearance (µm) requested
            from each adjoining solver's mode-field query, overriding that
            solver's own default (e.g. `FDTDSolver.dpml`). None (default)
            asks nothing special — each solver pads its query however it
            always has, so an interior interface whose straight run-out is
            thinner than that padding still raises the geometry `ValueError`
            documented above and gets skipped-with-a-WARNING as before.
            Set this to fix that specifically: `FDTDSolver.section_mode_fields()`
            only needs enough clearance to keep `geometry.section_prisms()`
            inside uniform material at the query window (`2 * _CORR_HALF_SPAN`
            = 0.1 µm wide) — it does not need a real PML footprint, since
            `get_eigenmode()` never time-steps (see
            `FDTDSolver._build_section_mode_sim()`'s docstring) — so a value
            well under the device's actual `body_stub`/similar interior
            clearance (e.g. 0.1 for the `body_stub=0.25` the coupler and both
            MMI hybrid examples use) fits where `dpml` alone would not.
            Forwarded only to solvers whose `section_mode_fields()` accepts a
            `pad` parameter (duck-typed like `neff_override`); a solver that
            doesn't (e.g. `EMESolver.section_modes()`, which needs no such
            padding at all — its stubs default to 0.0) is called exactly as
            if this were unset.
    """

    def __init__(
        self,
        solvers: List[Any],
        boundaries: List[float],
        interface_modes: List[int],
        mirror_pairs: Optional[Dict[int, int]] = None,
        mirror_axis: Optional[float] = None,
        mirror_check_samples: int = 25,
        mirror_check_tol: float = 1e-3,
        calibrate_body_neff: bool = False,
        calibrate_lengths: Optional[List[float]] = None,
        check_mode_correspondence: bool = False,
        mode_correspondence_action: Literal["raise", "warn"] = "raise",
        mode_correspondence_min_overlap: float = 0.9,
        mode_correspondence_res: Optional[float] = None,
        mode_correspondence_pad: Optional[float] = None,
    ):
        if len(solvers) < 2:
            raise ValueError("HybridSolver needs at least 2 solvers.")
        if len(boundaries) != len(solvers) + 1:
            raise ValueError(
                f"boundaries must have len(solvers)+1 = {len(solvers) + 1} "
                f"entries, got {len(boundaries)}."
            )
        if len(interface_modes) != len(solvers) - 1:
            raise ValueError(
                f"interface_modes must have len(solvers)-1 = {len(solvers) - 1} "
                f"entries, got {len(interface_modes)}."
            )
        self.solvers = solvers
        self.boundaries = boundaries
        self.interface_modes = interface_modes
        self.mirror_pairs = dict(mirror_pairs) if mirror_pairs else {}
        self.mirror_axis = (
            mirror_axis if mirror_axis is not None else (boundaries[0] + boundaries[-1]) / 2.0
        )
        self.mirror_check_samples = mirror_check_samples
        self.mirror_check_tol = mirror_check_tol
        self.calibrate_body_neff = bool(calibrate_body_neff)
        self.calibrate_lengths = (
            None if calibrate_lengths is None else [float(L) for L in calibrate_lengths]
        )
        if mode_correspondence_action not in ("raise", "warn"):
            raise ValueError(
                f"mode_correspondence_action must be 'raise' or 'warn', got "
                f"{mode_correspondence_action!r}."
            )
        self.check_mode_correspondence = bool(check_mode_correspondence)
        self.mode_correspondence_action = mode_correspondence_action
        self.mode_correspondence_min_overlap = float(mode_correspondence_min_overlap)
        self.mode_correspondence_res = (
            None if mode_correspondence_res is None else float(mode_correspondence_res)
        )
        self.mode_correspondence_pad = (
            None if mode_correspondence_pad is None else float(mode_correspondence_pad)
        )

        K = len(solvers)
        for k, j in self.mirror_pairs.items():
            if not (0 <= k < K) or not (0 <= j < K):
                raise ValueError(
                    f"mirror_pairs entries must be valid section indices in "
                    f"[0, {K - 1}], got {k}: {j}."
                )
            if k == j:
                raise ValueError(f"mirror_pairs[{k}]={j}: a section cannot mirror itself.")
            if j >= k:
                raise ValueError(
                    f"mirror_pairs[{k}]={j}: partner section index must be < {k} "
                    f"(the partner must be solved earlier in the cascade)."
                )
            if j in self.mirror_pairs:
                raise ValueError(
                    f"mirror_pairs[{k}]={j}: partner section {j} is itself a "
                    f"mirrored section (mirror_pairs[{j}]={self.mirror_pairs[j]}); "
                    f"chained/transitive mirroring is not supported — point to "
                    f"an independently-solved section."
                )

        n_in_list, n_out_list = self._mode_lists(1, 1)
        for k, j in self.mirror_pairs.items():
            if n_in_list[k] != n_out_list[j] or n_out_list[k] != n_in_list[j]:
                raise ValueError(
                    f"mirror_pairs[{k}]={j}: mode-count mismatch (section {k}: "
                    f"{n_in_list[k]}->{n_out_list[k]}, section {j}: "
                    f"{n_in_list[j]}->{n_out_list[j]}); transpose reuse requires "
                    f"n_in[{k}]==n_out[{j}] and n_out[{k}]==n_in[{j}]."
                )

    def _mode_lists(self, in_modes: int, out_modes: int) -> Tuple[List[int], List[int]]:
        n_in_list = [in_modes] + list(self.interface_modes)
        n_out_list = list(self.interface_modes) + [out_modes]
        return n_in_list, n_out_list

    def _stubs(self, k: int) -> Tuple[float, float]:
        solver = self.solvers[k]
        K = len(self.solvers)
        stub_left = solver.waveguide_length if (k == 0 and hasattr(solver, "waveguide_length")) else 0.0
        stub_right = solver.waveguide_length if (k == K - 1 and hasattr(solver, "waveguide_length")) else 0.0
        return stub_left, stub_right

    def _validate_mirror_pair(self, geometry, k: int, j: int, kwargs_k: dict, kwargs_j: dict) -> None:
        """Guard against a wrong/stale mirror_pairs entry before trusting it.

        A wrong pairing silently produces wrong physics (transpose of a
        section that isn't actually a mirror of its partner), so this checks
        the two sections' geometry actually IS mirror-symmetric for the
        kwargs in play this solve() call, not just that the pairing was
        declared.
        """
        x0_k, x1_k = self.boundaries[k], self.boundaries[k + 1]
        x0_j, x1_j = self.boundaries[j], self.boundaries[j + 1]
        struct_tol = 1e-6

        if abs((x1_k - x0_k) - (x1_j - x0_j)) > struct_tol:
            raise ValueError(
                f"mirror_pairs[{k}]={j}: section spans differ in length "
                f"({x1_k - x0_k:.6g} vs {x1_j - x0_j:.6g} µm)."
            )

        m = self.mirror_axis
        if abs((2 * m - x1_k) - x0_j) > struct_tol or abs((2 * m - x0_k) - x1_j) > struct_tol:
            raise ValueError(
                f"mirror_pairs[{k}]={j}: sections are not positioned as mirror "
                f"images about mirror_axis={m:.6g} (section {k}="
                f"[{x0_k:.6g},{x1_k:.6g}], section {j}=[{x0_j:.6g},{x1_j:.6g}])."
            )

        stub_left_k, stub_right_k = self._stubs(k)
        stub_left_j, stub_right_j = self._stubs(j)
        if (abs(stub_left_k - stub_right_j) > struct_tol
                or abs(stub_right_k - stub_left_j) > struct_tol):
            raise ValueError(
                f"mirror_pairs[{k}]={j}: stub lengths inconsistent under "
                f"mirroring (section {k} L/R={stub_left_k:.4g}/{stub_right_k:.4g}, "
                f"section {j} L/R={stub_left_j:.4g}/{stub_right_j:.4g}); "
                f"solvers[{k}].waveguide_length must equal solvers[{j}].waveguide_length."
            )

        for x_k in np.linspace(x0_k, x1_k, self.mirror_check_samples):
            x_j = 2 * m - x_k
            w_k = geometry.width_at(x_k, **kwargs_k)
            w_j = geometry.width_at(x_j, **kwargs_j)
            if abs(w_k - w_j) > self.mirror_check_tol:
                raise ValueError(
                    f"mirror_pairs[{k}]={j}: width_at mismatch at x={x_k:.4g} "
                    f"(section {k}, kwargs={kwargs_k}) vs mirrored x={x_j:.4g} "
                    f"(section {j}, kwargs={kwargs_j}): {w_k:.6g} vs {w_j:.6g} µm "
                    f"(tol={self.mirror_check_tol:.2g}). Check input_idx/output_idx "
                    f"and any asymmetric geometry_kwargs for this pair."
                )
            c_k = geometry.center_at(x_k, **kwargs_k)
            c_j = geometry.center_at(x_j, **kwargs_j)
            if abs(c_k - c_j) > self.mirror_check_tol:
                raise ValueError(
                    f"mirror_pairs[{k}]={j}: center_at mismatch at x={x_k:.4g} vs "
                    f"mirrored x={x_j:.4g}: {c_k:.6g} vs {c_j:.6g} µm "
                    f"(tol={self.mirror_check_tol:.2g})."
                )

    @staticmethod
    def _cladding_index(solver_a, solver_b, wl: float) -> Optional[float]:
        """Cladding refractive index at `wl`, for separating guided modes
        from radiation modes. Read off either adjoining solver's own
        `material_stack` (both FDTDSolver and EMESolver carry one);
        `MaterialStack.sio2_index()` is the meep-independent accessor for
        exactly this. None if neither solver exposes one, in which case the
        guided-mode filter is skipped rather than guessed at."""
        for solver in (solver_a, solver_b):
            stack = getattr(solver, "material_stack", None)
            index_fn = getattr(stack, "sio2_index", None)
            if callable(index_fn):
                return float(index_fn(wl))
        return None

    @staticmethod
    def _interface_domain_size(solver) -> Optional[Tuple[float, float]]:
        """(width, height) of a solver's own simulation domain (µm), duck-typed
        across solver types — FDTDSolver's sim_width/sim_height, EMESolver's
        oxide_width/oxide_height. None if neither pair is present."""
        width = getattr(solver, "sim_width", None)
        height = getattr(solver, "sim_height", None)
        if width is None or height is None:
            width = getattr(solver, "oxide_width", None)
            height = getattr(solver, "oxide_height", None)
        if width is None or height is None:
            return None
        return float(width), float(height)

    def _interface_diagnostic_res(self, solver_a, solver_b) -> float:
        """Diagnostic-mesh resolution (pixels/µm) for a correspondence check
        between these two solvers.

        **This must be a resolution, not a fixed pixel count** — that was a
        real bug: a fixed 24-pixel-wide grid works out to only ~4.8 px/µm on
        a 5 µm-wide domain, leaving SOI_160's 0.16 µm core *0.77 cells* tall.
        An `EMESolver` side then SOLVES its modes on the mesh it is handed
        (`section_modes(mesh=...)`), so under-resolving it does not merely
        sample the comparison coarsely, it makes one side's modes physically
        wrong — measured directly on the modal multiplexer: worst diagonal
        overlap 0.0004 and modes 7-15 scrambled, versus 0.973-0.998 and clean
        ordering for the same interface at 40 px/µm. Note this is exactly why
        `FDTDSolver._iter_distinct_bands()`'s much coarser fixed grid IS
        sound: it compares meep against *itself*, where two identical
        eigenvectors overlap to 1.0 on any grid however coarse. A
        cross-engine comparison has no such immunity.

        Default: the coarser of the two solvers' own `res` (each already
        chosen to resolve that engine's physics for this material stack), so
        the check never asks for a finer mesh than the device itself is
        simulated at. `mode_correspondence_res` overrides.
        """
        if self.mode_correspondence_res is not None:
            return float(self.mode_correspondence_res)
        res_values = [float(r) for r in (getattr(solver_a, "res", None),
                                        getattr(solver_b, "res", None)) if r]
        return min(res_values) if res_values else _CORR_DEFAULT_RES

    @staticmethod
    def _interface_diagnostic_mesh(width: float, height: float, res: float):
        """An mw.Mesh2D spanning ~(width, height) at `res` pixels/µm.

        Pixels are exactly square (dz == dx == 1/res): meow's
        compute_modes_meep() (used whenever a chain's EME side runs
        mode_solver="meep", matching the FDTD side's own eigensolve engine)
        requires equal x/y spacing and raises otherwise — confirmed directly
        against the directional coupler's real body_solver. Each span is
        rounded to a whole number of pixels rather than dividing the span
        into a fixed point count, which is what guarantees that equality for
        an arbitrary width/height aspect ratio.
        """
        import meow as mw

        dx = 1.0 / float(res)
        ny = max(2, round(width / dx))
        nz = max(2, round(height / dx))
        w_used, h_used = ny * dx, nz * dx
        return mw.Mesh2D(
            x=np.linspace(-w_used / 2, w_used / 2, ny + 1),
            y=np.linspace(-h_used / 2, h_used / 2, nz + 1),
        )

    @staticmethod
    def _geometry_only_kwargs(solver, kwargs: dict) -> dict:
        """Drop any key in `kwargs` that binds to one of solver.solve_section()'s
        own named parameters (monitor_center, neff_override, save_dir, ...) —
        those are solve_section()-specific tuning, not geometry hooks, and
        may not exist on section_mode_fields()/section_modes()'s signatures."""
        solve_section_params = set(inspect.signature(solver.solve_section).parameters)
        return {key: val for key, val in kwargs.items() if key not in solve_section_params}

    def _interface_fields(
        self, solver, geometry, wavelength: float, x0: float, x1: float,
        n_modes: int, mesh, geometry_kwargs: dict,
    ) -> Tuple[Optional[Dict[str, np.ndarray]], Optional[np.ndarray]]:
        """(FieldSet, neff) for `solver`'s first n_modes modes at [x0, x1],
        sampled on `mesh`'s cell centers -- dispatching on whichever
        field-sampling hook `solver` implements. (None, None) if it
        implements neither (see check_mode_correspondence's docstring)."""
        if hasattr(solver, "section_mode_fields"):
            y_grid, z_grid = mesh_cell_centers(mesh)
            if self.mode_correspondence_pad is not None and _supports_pad(solver.section_mode_fields):
                geometry_kwargs = dict(geometry_kwargs, pad=self.mode_correspondence_pad)
            result = solver.section_mode_fields(
                geometry, wavelength, x0, x1, n_modes, y_grid, z_grid, **geometry_kwargs)
            fields = {name: result[name] for name in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")}
            return fields, np.asarray(result["neff"], dtype=float)
        if hasattr(solver, "section_modes"):
            modes, _ = solver.section_modes(
                geometry, wavelength, x0, x1, n_modes, mesh=mesh, **geometry_kwargs)
            neff = np.array([float(np.real(m.neff)) for m in modes])
            return stack_meow_mode_fields(modes), neff
        return None, None

    def _check_interface_correspondence(
        self, geometry, k: int, kwargs_k: dict, kwargs_k1: dict, wl_center: float,
    ) -> None:
        """Verify solvers[k]'s n_out basis and solvers[k+1]'s n_in basis at
        their shared interface (boundaries[k+1]) name the same physical
        modes, in the same order -- the assumption solve()'s matrix-chain
        contraction otherwise makes with no runtime check. See
        check_mode_correspondence's docstring in __init__.
        """
        solver_a, solver_b = self.solvers[k], self.solvers[k + 1]
        n_modes = self.interface_modes[k]
        x_iface = self.boundaries[k + 1]

        has_hook = lambda s: hasattr(s, "section_mode_fields") or hasattr(s, "section_modes")
        if not (has_hook(solver_a) and has_hook(solver_b)):
            logger.info(
                "check_mode_correspondence: solvers[%d] (%s) / solvers[%d] (%s) "
                "-- at least one has no field-sampling hook, skipping the "
                "interface-%d correspondence check.",
                k, type(solver_a).__name__, k + 1, type(solver_b).__name__, k,
            )
            return

        domain_a = self._interface_domain_size(solver_a)
        domain_b = self._interface_domain_size(solver_b)
        if domain_a is None and domain_b is None:
            logger.info(
                "check_mode_correspondence: neither solvers[%d] nor solvers[%d] "
                "exposes a domain size (sim_width/sim_height or oxide_width/"
                "oxide_height) -- skipping the interface-%d correspondence check.",
                k, k + 1, k,
            )
            return
        if domain_a is not None and domain_b is not None:
            width, height = max(domain_a[0], domain_b[0]), max(domain_a[1], domain_b[1])
        else:
            width, height = domain_a if domain_a is not None else domain_b

        res = self._interface_diagnostic_res(solver_a, solver_b)
        mesh = self._interface_diagnostic_mesh(width, height, res)
        x0, x1 = x_iface - _CORR_HALF_SPAN, x_iface + _CORR_HALF_SPAN

        try:
            fields_a, neff_a = self._interface_fields(
                solver_a, geometry, wl_center, x0, x1, n_modes, mesh,
                self._geometry_only_kwargs(solver_a, kwargs_k))
            fields_b, neff_b = self._interface_fields(
                solver_b, geometry, wl_center, x0, x1, n_modes, mesh,
                self._geometry_only_kwargs(solver_b, kwargs_k1))
        except ValueError as exc:
            # A field-sampling solver pads its own query window (e.g.
            # FDTDSolver adds +/-dpml on both sides automatically) -- on a
            # section whose interior interface has less straight run-out
            # than that padding needs (a real device/config combination,
            # confirmed against the directional coupler's BODY_STUB=0.25 <
            # dpml=0.5), the geometry itself raises rather than silently
            # measuring the wrong cross-section. That is a limitation of
            # THIS check at this interface, not a correspondence finding --
            # always skip-with-a-warning here regardless of
            # mode_correspondence_action.
            logger.warning(
                "check_mode_correspondence: could not sample fields at "
                "interface %d (x=%.4g µm) -- the query window likely reaches "
                "past this interface's straight run-out into a bend/taper. "
                "Skipping this interface's check. Underlying error: %s",
                k, x_iface, exc,
            )
            return

        # Restrict to GUIDED modes before comparing. A mode below the
        # cladding index is a radiation mode: it is not a discrete physical
        # mode at all, its shape is set by the simulation box and each
        # engine's own discretization, so two engines' radiation modes are
        # essentially unrelated and overlap ~0 no matter how well the real
        # basis corresponds. Skipping this filter makes the check
        # false-alarm on every device in the repo -- validation_markdown/
        # README.md's "Finding 1" records that EVERY interface basis here
        # contains at least one unguided mode. Measured on the modal
        # multiplexer: without the filter, worst diagonal overlap 0.0000,
        # driven entirely by bands 13-15 (n_eff 1.434/1.413/1.388 against a
        # 1.4468 cladding index); the 13 guided bands are unaffected.
        # (Those unguided modes ARE still contracted by solve() -- that is
        # its own problem, the one Finding 1 raises, and not something a
        # correspondence check can answer.)
        n_clad = self._cladding_index(solver_a, solver_b, wl_center)
        if n_clad is not None and neff_a is not None and neff_b is not None:
            guided_a, guided_b = neff_a > n_clad, neff_b > n_clad
            n_drop = int((~guided_a).sum() + (~guided_b).sum())
            if n_drop:
                logger.info(
                    "check_mode_correspondence: interface %d — excluding %d "
                    "unguided mode(s) (n_eff below the %.4f cladding index) "
                    "from the comparison: %d/%d guided on the solvers[%d] "
                    "side, %d/%d on the solvers[%d] side. Radiation modes "
                    "have no engine-independent shape to match.",
                    k, n_drop, n_clad, int(guided_a.sum()), len(neff_a), k,
                    int(guided_b.sum()), len(neff_b), k + 1,
                )
            fields_a = {key: val[guided_a] for key, val in fields_a.items()}
            fields_b = {key: val[guided_b] for key, val in fields_b.items()}
            n_modes = int(min(guided_a.sum(), guided_b.sum()))
            if n_modes == 0:
                logger.warning(
                    "check_mode_correspondence: interface %d has no guided "
                    "modes on at least one side — skipping.", k,
                )
                return

        na, nb = fields_a["Ey"].shape[0], fields_b["Ey"].shape[0]
        if min(na, nb) < n_modes:
            message = (
                f"HybridSolver: mode-correspondence check at interface {k} "
                f"(x={x_iface:.4g} µm) needed {n_modes} comparable modes but "
                f"solvers[{k}]={type(solver_a).__name__} supplied {na} and "
                f"solvers[{k + 1}]={type(solver_b).__name__} supplied {nb}."
            )
            if self.mode_correspondence_action == "raise":
                raise RuntimeError(message)
            logger.warning(message)
            return

        dy = float(mesh.x[1] - mesh.x[0])
        dz = float(mesh.y[1] - mesh.y[0])
        eta = field_overlap_matrix(fields_a, fields_b, dy=dy, dz=dz)

        # Same-engine controls first. Each side against itself MUST be the
        # identity -- these are eigenmodes of one operator on one grid. If
        # either fails, the fault is in the MEASUREMENT (comparison grid too
        # coarse for the solve, bad axis remap, a duplicate eigenvector,
        # non-orthogonal radiation modes left in) and any cross-engine number
        # from the same data is meaningless. Reported and bailed on rather
        # than allowed to masquerade as a device finding -- 2026-07-30
        # produced two such false alarms, and this control is what localized
        # both. Free: reuses the modes already solved.
        for side, fields, solver in ((k, fields_a, solver_a), (k + 1, fields_b, solver_b)):
            control = self_overlap_diagnostics(fields, dy=dy, dz=dz)
            if not control["is_orthogonal"] and not control["has_duplicate"]:
                # Mild non-orthogonality WITHOUT a duplicate: real and benign
                # on a symmetric multimode cross-section, where near-degenerate
                # pairs mix. Measured at the modal multiplexer's interface 0 --
                # meep 0.1257, meow 0.1886 -- while the cross-engine comparison
                # matched at diagonal 1.0000 across every mode that had a
                # partner. Report it; do NOT disqualify the basis, or the check
                # switches itself off on exactly the devices it exists for.
                logger.info(
                    "check_mode_correspondence: solvers[%d] (%s) at interface %d "
                    "has mildly non-orthogonal modes (max self-off-diagonal "
                    "%.4f) -- expected where near-degenerate pairs mix, not a "
                    "duplicate. Proceeding with the comparison.",
                    side, type(solver).__name__, k, control["max_off_diagonal"],
                )
            if control["has_duplicate"] or not control["diagonal_is_unit"]:
                logger.warning(
                    "check_mode_correspondence: SELF-OVERLAP CONTROL FAILED for "
                    "solvers[%d] (%s) at interface %d -- its own modes are not "
                    "mutually orthogonal (worst self-diagonal %.4f, max "
                    "off-diagonal %.4f). This is a measurement fault, not a "
                    "device finding: suspect the comparison mesh resolution "
                    "(mode_correspondence_res=%s), a duplicate eigenmode, or "
                    "unfiltered radiation modes. Skipping the cross-engine "
                    "comparison at this interface, whose result would be "
                    "meaningless.",
                    side, type(solver).__name__, k, control["worst_diagonal"],
                    control["max_off_diagonal"], self.mode_correspondence_res,
                )
                return

        # Full-matrix diagnostics, not just the matched diagonal: row
        # completeness is what separates "a mode is MISSING from the other
        # basis" (truncation -- eta has no partner for it at all) from "the
        # modes are all present but mis-ordered". The diagonal alone cannot
        # draw that distinction, and the two have completely different fixes.
        # See beef.mode_identity.overlap_diagnostics() and
        # validation_markdown/fdtd_vs_EME_mode_validation.md layer 1.
        diag = overlap_diagnostics(eta)
        identity = list(range(n_modes))
        checked = diag["correspondence"][:n_modes]
        rows = slice(0, n_modes)
        worst_diag = float(np.nanmin(diag["diagonal"][rows])) if n_modes else 1.0
        worst_complete = float(np.nanmin(diag["row_completeness"][rows])) if n_modes else 1.0
        min_ok = self.mode_correspondence_min_overlap

        # H-only overlap, as an independent corroboration of the E-only
        # result above (fdtd_vs_EME_mode_validation.md item 4). Both fields
        # are already sampled -- section_mode_fields()/stack_meow_mode_fields()
        # always return all six components -- so this is free: no extra
        # eigensolve, just field_overlap_matrix() pointed at Hy/Hz instead of
        # Ey/Ez. Reported alongside the E-only verdict, not gated on: unlike
        # the same-engine self-overlap control above, a disagreement here is
        # not by itself proof of anything (see the module's open question on
        # per-engine E:H calibration), so it does not raise and does not
        # affect mode_correspondence_action -- only the E-only comparison
        # does that, as validated on real devices to date.
        eta_h = field_overlap_matrix(fields_a, fields_b, dy=dy, dz=dz, components=("Hy", "Hz"))
        diag_h = overlap_diagnostics(eta_h)
        checked_h = diag_h["correspondence"][:n_modes]
        worst_diag_h = float(np.nanmin(diag_h["diagonal"][rows])) if n_modes else 1.0
        worst_complete_h = float(np.nanmin(diag_h["row_completeness"][rows])) if n_modes else 1.0
        h_agrees = checked_h == checked
        logger.info(
            "check_mode_correspondence: interface %d — H-only corroboration: "
            "worst diagonal %.4f (E-only %.4f), worst row completeness %.4f "
            "(E-only %.4f), correspondence %s (%s the E-only assignment).",
            k, worst_diag_h, worst_diag, worst_complete_h, worst_complete,
            checked_h, "matches" if h_agrees else "DISAGREES WITH",
        )
        if not h_agrees:
            logger.warning(
                "check_mode_correspondence: interface %d — H-only overlap "
                "disagrees with the E-only assignment (H-only %s vs E-only "
                "%s). Not gated on (see fdtd_vs_EME_mode_validation.md item "
                "4), but worth investigating: either the E-only result above "
                "is not as clean as it looks, or this interface's H-field "
                "sampling has its own fault (comparison mesh resolution, axis "
                "remap, radiation-mode contamination).",
                k, checked_h, checked,
            )

        mis_ordered = checked != identity
        incomplete = worst_complete < min_ok
        low_overlap = worst_diag < min_ok

        if mis_ordered or incomplete or low_overlap:
            # Name the dominant symptom, because the remedies differ: an
            # incomplete row means one side's basis does not span the other's
            # (raise the mode count / align the band windows); a mis-ordering
            # with complete rows means every mode is present but indexed
            # differently (the contraction is mixing channels).
            if incomplete:
                worst_row = int(np.nanargmin(diag["row_completeness"][rows]))
                verdict = (
                    f"BASIS INCOMPLETE -- solvers[{k}]'s mode {worst_row} is only "
                    f"{worst_complete:.4f} representable in solvers[{k + 1}]'s basis "
                    f"(row completeness, 1.0 = fully spanned). One side carries a "
                    f"mode the other does not span at all, so this is a truncation / "
                    f"band-window mismatch rather than a re-ordering: the two sides "
                    f"are not describing the same set of modes."
                )
            elif mis_ordered:
                verdict = (
                    f"MIS-ORDERED -- every mode is present and well matched "
                    f"(worst row completeness {worst_complete:.4f}) but the pairing "
                    f"is not the identity, so solve()'s index-i-to-index-i "
                    f"contraction is mixing channels."
                )
            else:
                verdict = (
                    f"LOW OVERLAP -- ordering is the identity but the worst matched "
                    f"pair only reaches {worst_diag:.4f}; the bases correspond but "
                    f"the mode shapes disagree."
                )
            message = (
                f"HybridSolver: mode-correspondence check failed at interface {k} "
                f"(x={x_iface:.4g} µm, between solvers[{k}]={type(solver_a).__name__} "
                f"and solvers[{k + 1}]={type(solver_b).__name__}). {verdict} "
                f"correspondence={checked} (expected identity {identity}); "
                f"worst diagonal={worst_diag:.4f}, worst row completeness="
                f"{worst_complete:.4f} (min required={min_ok:.2f}). "
                f"Per-mode diagonal={np.round(diag['diagonal'][rows], 4).tolist()}, "
                f"row completeness={np.round(diag['row_completeness'][rows], 4).tolist()}, "
                f"margin={np.round(diag['margin'][rows], 4).tolist()}. "
                f"The per-mode SHAPE of this failure is the diagnostic -- high "
                f"modes degrading means under-resolution, one row collapsing means "
                f"a missing mode, all rows high on wrong columns means a pure "
                f"re-ordering. See validation_markdown/"
                f"fdtd_vs_EME_mode_validation.md."
            )
            if self.mode_correspondence_action == "raise":
                raise RuntimeError(message)
            logger.warning(message)

    def _check_all_interfaces(
        self, geometry, wavelengths: np.ndarray, kwargs_per_section: List[dict],
    ) -> None:
        if not self.check_mode_correspondence:
            return
        wl_center = float(np.mean(wavelengths))
        for k in range(len(self.solvers) - 1):
            self._check_interface_correspondence(
                geometry, k, kwargs_per_section[k], kwargs_per_section[k + 1], wl_center)

    def solve(
        self,
        geometry,
        wavelengths: np.ndarray,
        in_modes: int = 1,
        out_modes: int = 1,
        port: str = "out",
        geometry_kwargs: Optional[Any] = None,
        cache: Optional[Dict[Any, np.ndarray]] = None,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
    ) -> SMatrix:
        """Run every section solver and cascade the resulting transfer matrices.

        Args:
            geometry: device geometry implementing section_prisms/section_gds/
                width_at/center_at (see beef.geometry.base.DeviceGeometry).
            wavelengths: [M] wavelength array (µm).
            in_modes: mode count at the true left (input) port. Usually 1.
            out_modes: mode count at the true right (output) port. Usually 1.
            port: port label to store the result under in the returned SMatrix.
            geometry_kwargs: forwarded to solve_section() calls. Either a
                single dict applied to every section, a list of
                len(solvers) dicts giving different kwargs per section, or a
                sparse {section_idx: kwargs} dict mapping only the sections
                that need kwargs — unnamed sections get {} (e.g.
                {2: {"output_idx": 1}} targets just the section that
                actually branches, with no empty-dict padding — see cache
                below for why the per-section distinction matters; a dict
                whose keys are ints is always read as this sparse form, so
                a plain shared dict must not use int keys). Entries that
                match a solver's own solve_section() named parameters (e.g.
                FDTDSolver's monitor_center/monitor_size overrides) bind to
                those parameters rather than reaching the geometry hooks, so
                per-section solver options ride in the same dicts.
            cache: optional dict shared across multiple solve() calls (e.g.
                once per output port of a multi-output device). Sections are
                memoized by (solver identity, x0, x1, n_in, n_out, kwargs) —
                so a section whose own kwargs don't change between calls
                (typically every section upstream of the one that branches)
                is only ever actually simulated once.
            save_dir: optional persistent on-disk cache directory (see
                beef/cache/), forwarded to every section solver's
                solve_section()/solve_section_length_sweep() call that
                accepts a save_dir parameter (custom solver types that don't
                implement it are simply called without it — save_dir is
                opt-in per solver, not required by the module docstring's
                solve_section() contract). None/"" (default) disables it,
                same as cache=None. This is a separate mechanism from cache:
                cache dedupes within/across solve() calls in one process;
                save_dir persists results to disk across process runs. The
                two compose freely.
            cache_only: forwarded to every section solver's solve_section()/
                solve_section_length_sweep() call that accepts a cache_only
                parameter (same opt-in-per-solver rule as save_dir, via
                _supports_cache_only()). A section already served from the
                in-memory cache: is unaffected (that's already a hit); a
                section that must hit the solver raises CacheMissError on a
                disk-cache miss instead of running. False (default) is the
                normal behavior.

        Returns:
            SMatrix with coeff(port=port, in_mode=i, out_mode=o) = alpha_full[m].
        """
        wavelengths = np.asarray(wavelengths)
        kwargs_per_section = self._kwargs_per_section(geometry_kwargs)
        n_in_list, n_out_list = self._mode_lists(in_modes, out_modes)

        self._check_all_interfaces(geometry, wavelengths, kwargs_per_section)

        mats = self._solve_sections(
            geometry, wavelengths, kwargs_per_section, n_in_list, n_out_list, cache,
            save_dir=save_dir, cache_only=cache_only,
        )

        # Cascade: full[m] = mats[-1][m] @ ... @ mats[0][m]
        full = mats[0]
        for mat in mats[1:]:
            full = np.einsum("mij,mjk->mik", mat, full)  # [M, out, in]

        result = SMatrix(wavelengths=wavelengths)
        for o in range(out_modes):
            for i in range(in_modes):
                result.set(full[:, o, i], port=port, in_mode=i + 1, out_mode=o + 1)
        return result

    def solve_length_sweep(
        self,
        geometry,
        wavelengths: np.ndarray,
        lengths: np.ndarray,
        sweep_section: int,
        in_modes: int = 1,
        out_modes: int = 1,
        geometry_kwargs: Optional[Any] = None,
        cache: Optional[Dict[Any, np.ndarray]] = None,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
    ) -> np.ndarray:
        """Cascade with one section's length swept analytically.

        Section `sweep_section` is solved via its solver's
        solve_section_length_sweep() (EMESolver implements this; FDTDSolver
        does not) — modes are solved once per wavelength and each length is
        applied as an analytic propagation phase, so the whole sweep costs
        exactly the same as a single solve() call regardless of
        len(lengths). Every other section is solved exactly as in solve(),
        with byte-identical cache keys, so sections dedupe freely between
        solve() and solve_length_sweep() calls sharing one cache dict; the
        swept section's cache entry gets the lengths appended to its key so
        it can never collide with a single-length entry.

        The swept section must span a uniform-cross-section region
        (solve_section_length_sweep()'s precondition) — checked by sampling
        geometry.width_at/center_at across the span before solving — and
        must not appear in mirror_pairs on either side.

        Args:
            geometry: as in solve().
            wavelengths: [M] wavelength array (µm).
            lengths: [n_L] physical propagation lengths (µm) for the swept
                section, forwarded verbatim to solve_section_length_sweep().
            sweep_section: index of the section to sweep.
            in_modes/out_modes/geometry_kwargs/cache/save_dir/cache_only: as in solve().

        Returns:
            Complex ndarray [n_L, M, out_modes, in_modes] — the full-device
            transfer coefficients per length and wavelength. No SMatrix is
            built; slice one length and feed SMatrix.set() for that.
        """
        wavelengths = np.asarray(wavelengths)
        lengths = np.asarray(lengths, dtype=float)
        K = len(self.solvers)

        if lengths.ndim != 1 or lengths.size == 0:
            raise ValueError(
                f"lengths must be a non-empty 1-D array, got shape {lengths.shape}."
            )
        if not (0 <= sweep_section < K):
            raise ValueError(
                f"sweep_section must be a valid section index in [0, {K - 1}], "
                f"got {sweep_section}."
            )
        solver = self.solvers[sweep_section]
        if not hasattr(solver, "solve_section_length_sweep"):
            raise TypeError(
                f"solvers[{sweep_section}] ({type(solver).__name__}) does not "
                f"implement solve_section_length_sweep(); analytic length "
                f"sweeping needs an EMESolver-style solver in the swept slot "
                f"(FDTDSolver does not support it)."
            )
        if sweep_section in self.mirror_pairs or sweep_section in self.mirror_pairs.values():
            raise ValueError(
                f"sweep_section={sweep_section} appears in mirror_pairs "
                f"({self.mirror_pairs}); the swept section cannot be mirrored "
                f"or be a mirror partner — sweep an unmirrored section."
            )

        kwargs_per_section = self._kwargs_per_section(geometry_kwargs)
        n_in_list, n_out_list = self._mode_lists(in_modes, out_modes)

        self._validate_uniform_sweep_section(
            geometry, sweep_section, kwargs_per_section[sweep_section]
        )
        self._check_all_interfaces(geometry, wavelengths, kwargs_per_section)

        mats = self._solve_sections(
            geometry, wavelengths, kwargs_per_section, n_in_list, n_out_list, cache,
            sweep_section=sweep_section, lengths=lengths, save_dir=save_dir,
            cache_only=cache_only,
        )

        n_L = len(lengths)
        full = None
        for mat in mats:
            if mat.ndim == 3:
                mat = np.broadcast_to(mat, (n_L,) + mat.shape)  # zero-stride view
            full = mat if full is None else np.einsum("Lmij,Lmjk->Lmik", mat, full)
        return np.ascontiguousarray(full)  # [n_L, M, out, in]

    def _kwargs_per_section(self, geometry_kwargs: Optional[Any]) -> List[dict]:
        K = len(self.solvers)
        if geometry_kwargs is None:
            return [{}] * K
        if isinstance(geometry_kwargs, dict):
            if any(isinstance(key, int) for key in geometry_kwargs):
                # Sparse per-section mapping {section_idx: kwargs} — sections
                # not named get {}. Saves callers from padding a full list
                # with empties just to target one section.
                if not all(isinstance(key, int) for key in geometry_kwargs):
                    raise ValueError(
                        f"geometry_kwargs dict mixes int section indices with "
                        f"other keys ({sorted(map(repr, geometry_kwargs))}); use "
                        f"either a sparse {{section_idx: kwargs}} mapping or a "
                        f"single plain kwargs dict, not both."
                    )
                bad = [key for key in geometry_kwargs if not (0 <= key < K)]
                if bad:
                    raise ValueError(
                        f"geometry_kwargs section indices must be in [0, {K - 1}], "
                        f"got {sorted(bad)}."
                    )
                for key, val in geometry_kwargs.items():
                    if not isinstance(val, dict):
                        raise ValueError(
                            f"geometry_kwargs[{key}] must be a kwargs dict, "
                            f"got {type(val).__name__}."
                        )
                return [dict(geometry_kwargs.get(k, {})) for k in range(K)]
            return [geometry_kwargs] * K
        kwargs_per_section = list(geometry_kwargs)
        if len(kwargs_per_section) != K:
            raise ValueError(
                f"geometry_kwargs list must have len(solvers) = {K} "
                f"entries, got {len(kwargs_per_section)}."
            )
        return kwargs_per_section

    def _validate_uniform_sweep_section(self, geometry, k: int, kwargs: dict) -> None:
        """solve_section_length_sweep() extrapolates an analytic propagation
        phase, which is only physical over a constant cross-section — a
        taper in the swept span would silently produce wrong results."""
        x0, x1 = self.boundaries[k], self.boundaries[k + 1]
        xs = np.linspace(x0, x1, self.mirror_check_samples)
        w0 = geometry.width_at(xs[0], **kwargs)
        c0 = geometry.center_at(xs[0], **kwargs)
        for x in xs[1:]:
            if (abs(geometry.width_at(x, **kwargs) - w0) > self.mirror_check_tol
                    or abs(geometry.center_at(x, **kwargs) - c0) > self.mirror_check_tol):
                raise ValueError(
                    f"sweep_section={k}: cross-section is not uniform over "
                    f"[{x0:.4g}, {x1:.4g}] (width_at/center_at varies at "
                    f"x={x:.4g}); solve_section_length_sweep() requires a "
                    f"uniform region."
                )

    def _calibration_source(self, k: int):
        """The FDTD solver whose β section `k`'s body should be propagated on,
        or None if section `k` needs no calibration.

        Only sections that cannot measure their own propagation are calibrated
        — i.e. those lacking `section_mode_beta()`, which in practice means the
        EME bodies. The β comes from a *neighbouring* FDTD section's solver, so
        the body inherits the discretization of the sections it is spliced to,
        which is the entire point: the chain has to agree with itself.

        Prefers the left neighbour purely for determinism; in every chain in
        this repo the sections flanking a body share solver settings, so the
        choice does not change the answer. If it ever does, that is a
        misconfigured chain and worth the difference being visible.
        """
        if not self.calibrate_body_neff:
            return None
        if hasattr(self.solvers[k], "section_mode_beta"):
            return None
        for j in (k - 1, k + 1):
            if 0 <= j < len(self.solvers) and hasattr(self.solvers[j], "section_mode_beta"):
                return self.solvers[j]
        return None

    def _calibrated_neff(
        self, geometry, wavelengths, k, n_modes, kwargs, save_dir, cache_only,
    ) -> Optional[np.ndarray]:
        """n_eff for section `k`'s body, measured from a neighbour's FDTD phase.

        Returns ``[len(wavelengths), n_modes]`` or None when calibration is off
        or no neighbouring FDTD solver exists.

        The conversion is just ``neff = β·λ/2π`` — meow's propagation applies
        ``exp(2j·π·neff/wl·L)``, so handing it this n_eff makes the body
        accumulate exactly the FDTD's own phase per µm.
        """
        source = self._calibration_source(k)
        if source is None:
            return None
        x0, x1 = self.boundaries[k], self.boundaries[k + 1]
        # Only the parameters section_mode_beta() actually reads — passing the
        # EME section's full kwargs would hand it monitor/aperture settings
        # meant for a different solver.
        passthrough = {key: kwargs[key] for key in ("monitor_center", "monitor_size")
                       if key in kwargs}
        logger.info(
            "Section %d/%d: calibrating body n_eff from %s's own propagation "
            "phase (dc_validation.md item 6)",
            k + 1, len(self.solvers), type(source).__name__,
        )
        beta = source.section_mode_beta(
            geometry, wavelengths, x0, x1, n_modes,
            lengths=self.calibrate_lengths,
            save_dir=save_dir, cache_only=cache_only,
            **passthrough,
        )
        return beta * np.asarray(wavelengths, dtype=float)[:, None] / (2.0 * np.pi)

    def _solve_sections(
        self,
        geometry,
        wavelengths: np.ndarray,
        kwargs_per_section: List[dict],
        n_in_list: List[int],
        n_out_list: List[int],
        cache: Optional[Dict[Any, np.ndarray]],
        sweep_section: Optional[int] = None,
        lengths: Optional[np.ndarray] = None,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
    ) -> List[np.ndarray]:
        K = len(self.solvers)
        mats = []
        for k, solver in enumerate(self.solvers):
            x0, x1 = self.boundaries[k], self.boundaries[k + 1]
            stub_left, stub_right = self._stubs(k)
            kwargs = kwargs_per_section[k]
            is_swept = k == sweep_section

            cache_key = None
            if cache is not None:
                # calibrate_body_neff belongs in the key for the same reason it
                # belongs in the on-disk config: the same section propagated on
                # a different n_eff is a different matrix, and this dict is
                # shared across the port-by-port solves of one device.
                cache_key = (id(solver), x0, x1, n_in_list[k], n_out_list[k],
                             self.calibrate_body_neff,
                             tuple(sorted(kwargs.items())))
                if is_swept:
                    cache_key = cache_key + (("lengths", tuple(float(L) for L in lengths)),)
                if cache_key in cache:
                    logger.info("Section %d/%d: %s (cached)", k + 1, K, type(solver).__name__)
                    mats.append(cache[cache_key])
                    continue

            if k in self.mirror_pairs:
                j = self.mirror_pairs[k]
                self._validate_mirror_pair(geometry, k, j, kwargs, kwargs_per_section[j])
                mat = mats[j].transpose(0, 2, 1).copy()
                logger.info("Section %d/%d: %s (mirror of section %d)", k + 1, K, type(solver).__name__, j + 1)
            elif is_swept:
                logger.info(
                    "Section %d/%d: %s x=[%.3f, %.3f]  modes %s->%s (length sweep, %d lengths)",
                    k + 1, K, type(solver).__name__, x0, x1, n_in_list[k], n_out_list[k], len(lengths),
                )
                sweep_kwargs = dict(kwargs)
                if save_dir and _supports_save_dir(solver.solve_section_length_sweep):
                    sweep_kwargs["save_dir"] = save_dir
                if cache_only and _supports_cache_only(solver.solve_section_length_sweep):
                    sweep_kwargs["cache_only"] = cache_only
                if _supports_neff_override(solver.solve_section_length_sweep):
                    neff = self._calibrated_neff(
                        geometry, wavelengths, k, max(n_in_list[k], n_out_list[k]),
                        kwargs, save_dir, cache_only)
                    if neff is not None:
                        sweep_kwargs["neff_override"] = neff
                mat = solver.solve_section_length_sweep(
                    geometry, wavelengths, x0, x1, lengths,
                    n_modes_in=n_in_list[k], n_modes_out=n_out_list[k],
                    stub_left=stub_left, stub_right=stub_right,
                    **sweep_kwargs,
                )
                mat = np.asarray(mat, dtype=complex)  # [n_L, M, out, in]
            else:
                logger.info(
                    "Section %d/%d: %s x=[%.3f, %.3f]  modes %s->%s",
                    k + 1, K, type(solver).__name__, x0, x1, n_in_list[k], n_out_list[k],
                )
                section_kwargs = dict(kwargs)
                if save_dir and _supports_save_dir(solver.solve_section):
                    section_kwargs["save_dir"] = save_dir
                if cache_only and _supports_cache_only(solver.solve_section):
                    section_kwargs["cache_only"] = cache_only
                if _supports_neff_override(solver.solve_section):
                    neff = self._calibrated_neff(
                        geometry, wavelengths, k, max(n_in_list[k], n_out_list[k]),
                        kwargs, save_dir, cache_only)
                    if neff is not None:
                        section_kwargs["neff_override"] = neff
                mat = solver.solve_section(
                    geometry, wavelengths, x0, x1,
                    n_modes_in=n_in_list[k], n_modes_out=n_out_list[k],
                    stub_left=stub_left, stub_right=stub_right,
                    **section_kwargs,
                )
                mat = np.asarray(mat, dtype=complex)

            if cache is not None:
                cache[cache_key] = mat
            mats.append(mat)
        return mats
