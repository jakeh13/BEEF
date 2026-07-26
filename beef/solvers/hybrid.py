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
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..results.smatrix import SMatrix
from .._log import get_logger

logger = get_logger(__name__)


def _supports_save_dir(method) -> bool:
    return "save_dir" in inspect.signature(method).parameters


def _supports_cache_only(method) -> bool:
    return "cache_only" in inspect.signature(method).parameters


def _supports_neff_override(method) -> bool:
    return "neff_override" in inspect.signature(method).parameters


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
