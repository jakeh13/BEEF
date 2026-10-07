# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Georgia Tech Research Corporation
"""EMESolver: meow-based eigenmode expansion solver."""
from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional, Sequence

import numpy as np

from .._utils import extrude_polygons
from ..cache import (
    CacheMissError, build_run_config, build_whole_device_config, load_result, lookup,
    lookup_superset, save_run, update_family_pointer,
)
from ..cache.modes import extract_factorization, s_matrix_from_factorization
from ..geometry.base import DeviceGeometry, assert_uniform_cross_section
from ..materials.stack import MaterialStack, resolved_material_config, SOI_160
from ..results.smatrix import SMatrix
from .._log import get_logger

logger = get_logger(__name__)

# Duplicate-mode guard for compute_modes(). Mirrors fdtd.py's
# _DUP_BAND_OVERLAP_THRESHOLD / _DUP_BAND_MAX_EXTRA -- same failure, other
# engine. Threshold is deliberately well above the 0.10-0.19 self-overlap
# genuinely near-degenerate mode PAIRS show on a symmetric cross-section
# (measured on the modal multiplexer): those are real, distinct modes and must
# not be discarded. Only a near-exact 1.0 is a duplicate.
_DUP_MODE_OVERLAP_THRESHOLD = 0.9
_DUP_MODE_MAX_EXTRA = 8


def _neff_fingerprint(neff_override: Optional[np.ndarray]):
    """JSON-safe cache-key form of a `neff_override` array, or None.

    Must be part of every EME cache config that accepts the override: the same
    geometry, mesh and wavelengths propagated on a *different* n_eff is a
    different physical result, and without this a calibrated run and an
    uncalibrated one would share a cache entry and silently serve each other.
    Rounded to the same 9 decimals beef.cache.config uses elsewhere so a
    JSON round-trip cannot turn a hit into a miss.
    """
    if neff_override is None:
        return None
    arr = np.asarray(neff_override)
    return [[[round(float(np.real(v)), 9), round(float(np.imag(v)), 9)] for v in row]
            for row in np.atleast_2d(arr)]


def _cell_neff_fingerprint(cell_neff_override: Optional[np.ndarray]):
    """JSON-safe cache-key form of a `cell_neff_override` array, or None.

    Same reasoning as `_neff_fingerprint` — a section propagated on different
    per-cell n_eff is a different physical result and must not share a cache
    entry — but one axis deeper: [wavelength, cell, mode] rather than
    [wavelength, mode]. Same 9-decimal rounding so a JSON round-trip cannot
    turn a hit into a miss.

    NOTE this is only ever ADDED to a config dict when the override is not
    None; see solve_section(). An uncalibrated run must keep the exact
    fingerprint it had before this parameter existed, or every cached FDTD
    section in the repo is invalidated at once.
    """
    if cell_neff_override is None:
        return None
    arr = np.asarray(cell_neff_override)
    if arr.ndim != 3:
        raise ValueError(
            f"cell_neff_override must be [n_wavelengths, n_cells, n_modes]; "
            f"got shape {arr.shape}")
    # NaN is serialized as the STRING "nan", not a float NaN.
    #
    # json.dumps happily writes a bare NaN token and json.load reads it back as
    # float('nan') -- but NaN != NaN, so beef.cache.store.lookup()'s
    # `cached_config != config` guard then NEVER matches and every calibrated
    # run re-solves from scratch. Measured: the first per-cell calibrated run
    # logged "hash match ... but config.json differs" on all 8 ranks, because
    # modes left uncalibrated carry NaN by design. A string sentinel compares
    # equal to itself and keeps the fingerprint faithful.
    def _fp(v):
        re_, im_ = float(np.real(v)), float(np.imag(v))
        if np.isnan(re_):
            return "nan"
        return [round(re_, 9), round(im_, 9)]

    return [[[_fp(v) for v in mode_row] for mode_row in cell_row]
            for cell_row in arr]


# Tolerance for the diagonal-reconstruction self-check in
# solve_section_length_sweep() -- see that method's "neff cache" tier. Loose
# enough to absorb meow's own solver noise, tight enough that a real ordering/
# convention mismatch (not just floating-point jitter) always fails it.
_DIAGONAL_RECONSTRUCTION_TOL = 1e-6
# Lengths per block when self-checking that reconstruction; bounds peak memory
# since T here can be hundreds of MB and the check builds a second copy.
_DIAGONAL_CHECK_CHUNK = 512
# solve_length_sweep() stores the cascade factorization only if replaying it
# reproduces that run's own amplitudes to this tolerance. Measured bit-exact
# (max|dS| = 0) on a real multi-cell solve, so this is a guard against a future
# meow change, not an expected error budget.
_FACTORIZATION_REPLAY_TOL = 1e-12


def _propagate_diagonal(neff_all: np.ndarray, wavelengths: np.ndarray, lengths: np.ndarray,
                         n_modes_out: int, n_modes_in: int) -> np.ndarray:
    """Reconstruct transfer matrices from per-wavelength mode effective
    indices alone, with no meow call at all.

    A uniform-cross-section EME section's propagation S-matrix is exactly
    diagonal in the mode basis (solve_section_length_sweep()'s own docstring):
    ``S_ii(L) = exp(2j*pi*neff_i(wl)/wl * L)``. Given just ``neff_all``
    (``[len(wavelengths), n_modes]``, already carrying any neff_override), the
    transfer matrix at *any* length is this formula, not a new eigenmode
    solve -- the reuse mechanism this function exists to serve, and also the
    ground truth solve_section_length_sweep()'s own real (meow) solve is
    checked against once per fresh solve before its neff is ever cached
    (see _DIAGONAL_RECONSTRUCTION_TOL).

    Args:
        neff_all: [len(wavelengths), n_modes] complex, mode index i in the
            same right@i/left@i basis solve_section_length_sweep()'s own
            mw.compute_s_matrix() calls use.
        wavelengths: [len(wavelengths)] um.
        lengths: [n_L] um.
        n_modes_out/n_modes_in: trailing shape of the returned array; only
            min(n_modes_out, n_modes_in, neff_all.shape[1]) diagonal entries
            are filled, matching a mode-count-truncated diagonal block —
            same truncation convention lookup_superset()'s slicing assumes.

    Returns:
        Complex ndarray [len(lengths), len(wavelengths), n_modes_out, n_modes_in].
    """
    wavelengths = np.asarray(wavelengths)
    lengths = np.asarray(lengths)
    n = min(n_modes_out, n_modes_in, neff_all.shape[1])
    T = np.zeros((len(lengths), len(wavelengths), n_modes_out, n_modes_in), dtype=complex)
    if n == 0:
        return T
    phase = np.exp(
        2j * np.pi * neff_all[None, :, :n] / wavelengths[None, :, None] * lengths[:, None, None]
    )  # [n_L, n_wl, n]
    for i in range(n):
        T[:, :, i, i] = phase[:, :, i]
    return T


class EMESolver:
    """Eigenmode expansion solver wrapping meow.

    Args:
        num_modes: modes per EME region, left to right.
                   e.g. [1, 7, 15, 7, 1] for (wg, taper, MMI, taper, wg).
                   Read by the whole-device methods (solve/solve_length_sweep)
                   ONLY — solve_section() takes its mode count from its own
                   n_modes_in/n_modes_out arguments instead, which is why
                   several hybrid scripts pass a single-element placeholder
                   here and why one shared solver can serve sections with
                   different interface mode counts.
        num_cells: cells per region (same length as num_modes).
                   e.g. [1, 12, 1, 12, 1].  Regions with num_cells > 1 are
                   subdivided equally — use this for tapers. Unlike num_modes,
                   this IS read by solve_section() (as num_cells[0]), which is
                   how a tapered body gets real multi-cell propagation inside a
                   hybrid cascade.
        oxide_width: y-span of the cross-section mesh (µm).
        oxide_height: z-span of the cross-section mesh (µm).
        res: mesh resolution in pts/µm.
        material_stack: MaterialStack instance for extrusion rules.
        straight_length: access waveguide length on each side (µm).
        temperature: simulation temperature for meow Environment (°C).
        mode_solver: which meow mode solver compute_modes() dispatches to —
            "meep" (meow.compute_modes_meep, the default) or "tidy3d"
            (meow.compute_modes_tidy3d). "meep" matches the same underlying
            eigenmode solver/conventions FDTDSolver's own
            get_eigenmode_coefficients() calls use in HybridSolver chains,
            avoiding basis mismatches between EME and FDTD sections — use
            that for hybrid tutorials. A standalone EME-only solve (no FDTD
            sections to stay consistent with) is free to use "tidy3d"
            instead. Does not affect solve_sweep(), which always uses
            tidy3d's batched mode solver regardless of this setting.
        guard_duplicate_modes: reject duplicate eigenvectors from meow, so an
            `n`-mode basis really carries `n` independent modes. **Default
            True** — this is a correctness guard against a measured bug, not a
            tuning knob (see `compute_modes()`), but it does change results
            wherever it fires, so it is switchable for before/after
            comparison. Costs one over-requested mode solve; the duplicate
            check itself is a cheap self-overlap on already-solved fields.
    """

    def __init__(
        self,
        num_modes: List[int],
        num_cells: List[int],
        lengths: List[float],
        oxide_width: float = 5.0,
        oxide_height: float = 1.5,
        res: float = 40.0,
        material_stack: MaterialStack = SOI_160,
        temperature: float = 25.0,
        mode_solver: str = "meep",
        guard_duplicate_modes: bool = True,
        dispersive: Optional[bool] = None,
        region_cell_lengths: Optional[Dict[int, Sequence[float]]] = None,
    ):
        if len(num_modes) != len(num_cells):
            raise ValueError(
                f"num_modes ({len(num_modes)}) and num_cells ({len(num_cells)}) must be same length"
            )
        if len(lengths) != len(num_modes):
            raise ValueError(
                f"lengths ({len(lengths)}) must match num_modes ({len(num_modes)})"
            )
        if mode_solver not in ("meep", "tidy3d"):
            raise ValueError(f"unknown mode_solver {mode_solver!r}, expected 'meep' or 'tidy3d'")
        self.region_cell_lengths = self._validate_region_cell_lengths(
            region_cell_lengths, num_cells, lengths)
        self.guard_duplicate_modes = bool(guard_duplicate_modes)
        self.dispersive = dispersive
        self.num_modes = num_modes
        self.num_cells = num_cells
        self.lengths = lengths
        self.straight_length = lengths[0]  # input access waveguide length
        self.oxide_width = oxide_width
        self.oxide_height = oxide_height
        self.res = res
        self.material_stack = material_stack
        self.temperature = temperature
        self.mode_solver = mode_solver

    @staticmethod
    def _validate_region_cell_lengths(region_cell_lengths, num_cells, lengths):
        """Normalize/validate the explicit per-region cell-length override.

        `region_cell_lengths` maps region index -> that region's own list of
        cell lengths, replacing the default "divide this region evenly by its
        num_cells". It exists for bodies whose cross-section changes at
        INTERNAL boundaries an even division cannot land on — e.g. an MMI with
        a waisted body, whose five segments (stub/transition/waist/transition/
        stub) must each get whole cells so no single cell straddles a width
        discontinuity.

        Guarded hard, because the two ways to get this wrong are both silent:
        a wrong cell COUNT would desync `num_modes`/`_modes_per_cell()` from
        the cells, and a wrong SUM would quietly change the device's physical
        length while every other number still looked plausible.
        """
        if region_cell_lengths is None:
            return None
        out: Dict[int, List[float]] = {}
        for idx, cells in region_cell_lengths.items():
            idx = int(idx)
            if not 0 <= idx < len(num_cells):
                raise ValueError(
                    f"region_cell_lengths key {idx} is out of range for "
                    f"{len(num_cells)} regions")
            cells = [float(c) for c in cells]
            if len(cells) != num_cells[idx]:
                raise ValueError(
                    f"region_cell_lengths[{idx}] has {len(cells)} cells but "
                    f"num_cells[{idx}]={num_cells[idx]} — they must match, or "
                    "the per-cell mode counts would desync from the cells")
            if any(c <= 0 for c in cells):
                raise ValueError(
                    f"region_cell_lengths[{idx}] has a non-positive cell length: {cells}")
            total = sum(cells)
            if abs(total - lengths[idx]) > 1e-9:
                raise ValueError(
                    f"region_cell_lengths[{idx}] sums to {total} but "
                    f"lengths[{idx}]={lengths[idx]} — an explicit subdivision "
                    "must tile its region exactly, not change the device length")
            out[idx] = cells
        return out

    def _resolved_dispersive(self) -> bool:
        """Whether this solver will actually track the band.

        `None` means **dispersive** here — see `_extrusion_rules`' own
        "asymmetric on purpose" note for why this does not defer to the stack.
        Split out 2026-08-01 so the cache fingerprint and the extrusion rules
        cannot disagree about it: passing the stack's flag to
        `resolved_material_config()` instead of this used to let an explicit
        `EMESolver(dispersive=False)` build frozen materials under a key that
        still said dispersive.
        """
        return True if self.dispersive is None else bool(self.dispersive)

    def _extrusion_rules(self, wavelengths):
        """meow extrusion rules for this run, dispersive unless frozen off.

        **The EME side has never had the FDTD's frozen-band problem**: meow
        re-evaluates a `SampledMaterial` at whatever wavelength each
        `Environment` names, and this solver builds one `Environment` per
        wavelength, so its materials are already exact Sellmeier everywhere.
        `dispersive` therefore only ever *removes* accuracy here.

        It exists for one reason: with `dispersive=False` on both engines, a
        hybrid A/B can hold the material model fixed and frozen on both sides,
        separating "dispersion changed the answer" from "the two engines were
        modelling different materials all along" — which today they are, by up
        to 0.008 in n_Si across the O-band.

        **The two engines' defaults now agree, as of 2026-08-01.** They did
        not before: the stack defaulted to False, so `FDTDSolver`'s `None`
        inherited frozen while this solver's `None` meant dispersive. That
        asymmetry *was* the material-model gap — a hybrid chain contracted an
        FDTD section frozen at band centre against an EME body solved exactly
        per wavelength, 0.00786 apart in n_Si at the O-band edges. The stack
        default is now True, so both `None`s resolve to dispersive.

        This solver still does not defer to the stack, deliberately: `None`
        means **dispersive** here regardless, and only an explicit `False`
        freezes. Inheriting would mean a future stack default flipping back
        silently downgrades every EME run from an exact model to a frozen one
        with no caller asking for it — the EME side has never had the FDTD's
        frozen-band problem and should not acquire it by inheritance.
        """
        dispersive = self._resolved_dispersive()
        if dispersive or wavelengths is None:
            # wavelengths=None means a caller that never had a band to freeze
            # at; the dispersive rules are correct there and are also what it
            # got before this flag existed.
            return self.material_stack.extrusion_rules
        return self.material_stack.extrusion_rules_at(
            float(np.mean(np.atleast_1d(np.asarray(wavelengths, dtype=float)))))

    def _build_mesh(self):
        import meow as mw

        return mw.Mesh2D(
            x=np.linspace(-self.oxide_width/2, self.oxide_width/2,
                          int(self.res * self.oxide_width) + 1),
            y=np.linspace(-self.oxide_height/2, self.oxide_height/2,
                          int(self.res * self.oxide_height) + 1),
        )

    def _build_cells(self, geometry: DeviceGeometry, side: str = "top",
                     wavelengths=None, **geometry_kwargs):
        import meow as mw

        gds_polys = geometry.to_gds_eme(
            straight_length=self.straight_length,
            oxide_width=self.oxide_width,
            x_start=0.0,
            side=side,
            **geometry_kwargs,
        )
        structs = extrude_polygons(gds_polys, self._extrusion_rules(wavelengths))

        Ls = self._flat_cell_lengths()

        mesh = self._build_mesh()
        return mw.create_cells(structures=structs, mesh=mesh, Ls=Ls)

    def _modes_per_cell(self) -> List[int]:
        """Flat list of mode counts, one entry per cell (not per region)."""
        out = []
        for n_modes, n_cells in zip(self.num_modes, self.num_cells):
            out.extend([n_modes] * n_cells)
        return out

    def _flat_cell_lengths(self, region_lengths: Optional[List[float]] = None) -> List[float]:
        """Flat per-cell length list (one entry per cell), from region lengths
        divided evenly across each region's num_cells. Defaults to
        self.lengths — the same Ls _build_cells() hands to mw.create_cells().

        A region with an explicit `region_cell_lengths` entry uses that list
        verbatim instead of an even division. Such an entry is validated at
        construction to tile its region exactly, so it is only usable while
        that region still has its configured length — see the guard below.
        """
        lengths = self.lengths if region_lengths is None else region_lengths
        overrides = self.region_cell_lengths or {}
        out: List[float] = []
        for idx, (region_len, n_cells) in enumerate(zip(lengths, self.num_cells)):
            explicit = overrides.get(idx)
            if explicit is not None:
                # The override tiles the region it was validated against. If a
                # caller is now varying that region's length (a length sweep,
                # or solve_section's single-region rebuild), the subdivision no
                # longer tiles it -- refuse rather than silently emit cells
                # that do not sum to the requested length.
                if abs(region_len - self.lengths[idx]) > 1e-9:
                    raise ValueError(
                        f"region {idx} has an explicit region_cell_lengths "
                        f"subdivision tiling length {self.lengths[idx]}, but is "
                        f"being solved at length {region_len}. An explicit "
                        "subdivision cannot be rescaled automatically -- build a "
                        "solver per length, the way a waisted-body sweep does.")
                out.extend(explicit)
            else:
                out.extend([region_len / n_cells] * n_cells)
        return out

    def _region_cell_span(self, region_idx: int) -> tuple:
        """(start, stop) indices of region_idx's cells within the flat cell list."""
        start = sum(self.num_cells[:region_idx])
        return start, start + self.num_cells[region_idx]

    def _infer_sweep_region(self) -> int:
        """The unique region index with num_cells > 1 — i.e. the subdivided
        taper region a length sweep varies. Raises if that's ambiguous."""
        variable = [i for i, n in enumerate(self.num_cells) if n > 1]
        if len(variable) != 1:
            raise ValueError(
                f"cannot infer sweep_region: expected exactly one region with "
                f"num_cells > 1, found indices {variable} (num_cells={self.num_cells}). "
                "Pass sweep_region explicitly."
            )
        return variable[0]

    @staticmethod
    def _override_neff(modes: list, neff_row) -> list:
        """Return `modes` with each mode's n_eff replaced by `neff_row[i]`.

        The only thing meow's propagation S-matrix reads is `Mode.neff` —
        `exp(2j·π·neff/wl·L)` in meow/eme/propagation.py — so substituting it
        here changes how the body propagates and nothing else. Field profiles,
        and therefore every interface overlap integral, are left exactly as the
        mode solver produced them. That split is the point: only the
        *propagation* is being taken from the FDTD (dc_validation.md item 6).

        `Mode` is a frozen pydantic model, so this copies rather than mutates —
        which also means a caller's `modes` list is never disturbed.
        """
        if len(neff_row) != len(modes):
            raise ValueError(
                f"neff_override has {len(neff_row)} entries but this "
                f"cross-section resolved {len(modes)} modes"
            )
        return [m.model_copy(update={"neff": complex(n)})
                for m, n in zip(modes, neff_row)]

    def _solve_modes_raw(self, cross_section, num_modes: int) -> list:
        """Unguarded dispatch to meow's mode solver — see compute_modes()."""
        if self.mode_solver == "meep":
            from meow import compute_modes_meep
            return compute_modes_meep(cross_section, num_modes=num_modes)
        from meow import compute_modes_tidy3d
        return compute_modes_tidy3d(cross_section, num_modes=num_modes)

    def compute_modes(self, cross_section, num_modes: int) -> list:
        """Compute modes for a single cross-section (exposed for caching).

        Dispatches to meow's MEEP/MPB-backed solver or its tidy3d-backed one
        per self.mode_solver (see the class docstring) — "meep" matches the
        same underlying eigenmode solver/conventions as FDTDSolver's own
        get_eigenmode_coefficients() calls in HybridSolver chains, avoiding
        basis mismatches between EME and FDTD sections.

        **Guards against meow returning the same eigenvector twice**, when
        `guard_duplicate_modes` is on (the default). This is a real, measured
        failure and not a hypothetical: at the modal multiplexer's interface 0,
        16 modes at res 40, meow's own modes 10 and 11 came back as the
        identical eigenvector (self-overlap 1.0000 exactly, row completeness
        2.0003 — that one mode counted twice). It is the same failure class
        `FDTDSolver._iter_distinct_bands()` was built for, on the other
        engine, where nothing was checking for it. The cost is a
        rank-deficient basis: one of the `num_modes` channels is a copy of
        another, and whatever mode the duplicate displaced is missing from the
        basis entirely — strictly worse than a re-ordering, and silent.

        Found by `beef.mode_identity.self_overlap_diagnostics()`, the
        same-engine control this now reuses to detect it.
        """
        if not self.guard_duplicate_modes:
            return self._solve_modes_raw(cross_section, num_modes)

        from ..mode_identity import field_overlap_matrix, stack_meow_mode_fields

        # Keep the first `num_modes` DISTINCT modes. Same policy as the
        # meep-side guard: a duplicate means the solver converged onto an
        # already-found eigenvector, so reach further up the spectrum rather
        # than accepting a rank-deficient basis. dy/dz cancel exactly in the
        # normalized ratio, so 1.0/1.0 is not an approximation.
        #
        # LAZY, as of 2026-08-25 -- this used to solve `num_modes +
        # _DUP_MODE_MAX_EXTRA` unconditionally, on every cross-section, whether
        # or not a duplicate was ever found. That is what the comment above the
        # constants means by "mirrors fdtd.py": the CONSTANT matched, but the
        # cost model did not. `_iter_distinct_bands()` treats its own
        # `n_modes + _DUP_BAND_MAX_EXTRA` as a CAP on a generator and stops as
        # soon as it has enough distinct bands; this path paid the cap up front
        # every time, so requesting 2 modes always solved 10.
        #
        # Escalating rather than over-requesting is safe because the
        # over-request is a PURE FILTER and never perturbs the low bands --
        # measured, not assumed: job 11586098 arm C vs arm A gave
        # max |delta n_eff| = 0.000e+00 across all 16 low bands, and recorded
        # the over-request as a 3.2x pure cost (running_jobs.md item 13, folded
        # into resolved/fdtd_vs_EME_mode_validation.md). So the top-k of a
        # k-mode solve and the top-k of a (k+8)-mode solve agree, and solving
        # k first costs nothing in accuracy.
        #
        # Two attempts, not a ladder: duplicates are rare, so the common case
        # must be exactly `num_modes` solves, and the rare case should not pay
        # for a long staircase of re-solves (each `_solve_modes_raw(n)` redoes
        # bands 1..n from scratch). Attempt 2 jumps straight to the old
        # unconditional count, so the fallback behaviour is bit-identical to
        # what this method did before.
        #
        # Cache note: neither `guard_duplicate_modes` nor _DUP_MODE_MAX_EXTRA
        # is part of any EME cache key (see beef/cache/config.py -- the EME
        # entries key on oxide_width/oxide_height/res/material_stack/
        # temperature/mode_solver/num_cells and the call's own arguments), so
        # this change re-keys nothing and every existing entry stays valid.
        attempts = [num_modes]
        if _DUP_MODE_MAX_EXTRA > 0:
            attempts.append(num_modes + _DUP_MODE_MAX_EXTRA)

        kept: List[int] = []
        raw: list = []
        for attempt, n_solve in enumerate(attempts):
            raw = self._solve_modes_raw(cross_section, n_solve)
            fields = stack_meow_mode_fields(raw)
            eta = field_overlap_matrix(fields, fields, dy=1.0, dz=1.0)

            kept = []
            for j in range(len(raw)):
                if len(kept) == num_modes:
                    break
                dup_of = next(
                    (i for i in kept if eta[i, j] > _DUP_MODE_OVERLAP_THRESHOLD), None)
                if dup_of is not None:
                    logger.warning(
                        "EMESolver.compute_modes: meow mode %d duplicates already-"
                        "accepted mode %d (self-overlap %.6f > %.2f) -- skipping it "
                        "and reaching further up the spectrum. A rank-deficient "
                        "basis would silently drop one channel of the %d requested.",
                        j, dup_of, eta[dup_of, j], _DUP_MODE_OVERLAP_THRESHOLD,
                        num_modes,
                    )
                    continue
                kept.append(j)

            if len(kept) == num_modes:
                break
            if attempt + 1 < len(attempts):
                logger.warning(
                    "EMESolver.compute_modes: only %d of %d requested modes were "
                    "distinct among %d solved -- re-solving with %d to reach "
                    "further up the spectrum.",
                    len(kept), num_modes, n_solve, attempts[attempt + 1],
                )

        if len(kept) < num_modes:
            raise RuntimeError(
                f"EMESolver.compute_modes: could not find {num_modes} distinct "
                f"modes within {attempts[-1]} solved modes at "
                f"this cross-section -- meow returned more duplicates than this "
                f"guard tolerates (kept {len(kept)})."
            )
        return [raw[j] for j in kept]

    def section_modes(
        self,
        geometry,
        wavelength: float,
        x0: float,
        x1: float,
        n_modes: int,
        mesh=None,
        stub_left: float = 0.0,
        stub_right: float = 0.0,
        **geometry_kwargs,
    ):
        """Raw meow ``Mode`` objects at this section's cross-section, no
        propagation — the field-sampling counterpart of the intermediate
        ``modes`` list ``solve_section()`` already computes internally but
        never surfaces (it goes straight into ``mw.compute_s_matrix()``).

        Reuses ``solve_section()``'s exact cell-building calls
        (``geometry.section_gds()``, ``extrude_polygons()``,
        ``mw.create_cells()``, ``mw.CrossSection.from_cell()``,
        ``self.compute_modes()``); `solve_section()` itself is untouched.

        Args:
            mesh: optional caller-supplied ``mw.Mesh2D`` (e.g. a shared
                diagnostic mesh independent of ``self.res``/``self.oxide_*``,
                for comparing against another solver's field grid). Defaults
                to ``self._build_mesh()``, matching ``solve_section()``.

        Returns:
            ``(modes, mesh)`` — pass ``mesh`` to
            ``beef.mode_identity.mesh_cell_centers()`` for the ``(y, z)``
            points these modes' field arrays sit at, and the mode list to
            ``beef.mode_identity.stack_meow_mode_fields()``.
        """
        import meow as mw

        length = (x1 - x0) + stub_left + stub_right
        gds_polys = geometry.section_gds(
            x0, x1, self.oxide_width,
            stub_left=stub_left, stub_right=stub_right, **geometry_kwargs
        )
        # Scalar wavelength: freezing at it and evaluating exactly at it are
        # the same material, so this site is inert under the flag either way.
        structs = extrude_polygons(gds_polys, self._extrusion_rules(wavelength))
        mesh = self._build_mesh() if mesh is None else mesh
        cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
        env = mw.Environment(wl=float(wavelength), T=self.temperature)
        cs = mw.CrossSection.from_cell(cell=cells[0], env=env)
        return self.compute_modes(cs, n_modes), mesh

    def solve_section(
        self,
        geometry,
        wavelengths: np.ndarray,
        x0: float,
        x1: float,
        n_modes_in: int = 1,
        n_modes_out: int = 1,
        stub_left: float = 0.0,
        stub_right: float = 0.0,
        neff_override: Optional[np.ndarray] = None,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        *,
        cell_neff_override: Optional[np.ndarray] = None,
        **geometry_kwargs,
    ) -> np.ndarray:
        """Modal transfer matrix for the device section [x0, x1] (device frame).

        EME over the section — mode expansion naturally gives the full
        n_modes_out x n_modes_in block in one solve (no per-mode source runs
        needed), unlike FDTD.

        The section is subdivided into ``self.num_cells[0]`` cells, so a
        **non-uniform** span (e.g. a tapered MMI body) is propagated through
        real interface overlaps between successive cross-sections rather than
        collapsed onto one. ``num_cells=[1]`` — the value every device script
        in this repo uses — is a single cell spanning the whole section,
        exactly what this method always did. The per-cell mode count is
        ``max(n_modes_in, n_modes_out)``; ``self.num_modes`` is **not** read
        here (see the real-solve branch's own comment on why, and
        ``modal_multiplexer_hybrid.py``, which relies on that).

        Note the `neff_override` path below is single-cell by construction and
        guarded to uniform spans — the two are different propagation models,
        not two settings of one.

        Args:
            neff_override: optional array ``[len(wavelengths), n_modes]``
                replacing each mode's n_eff before propagation, leaving field
                profiles — and so every interface overlap — untouched. Lets a
                HybridSolver body propagate on the β its FDTD sections
                actually accumulate rather than on MPB's mode index: the two
                differ by ~2.3% per mode at res 40, and because Δβ is a
                difference of two β near 9.4 rad/µm that lands as a ~10% error
                on the coupling (dc_validation.md item 6). Normally supplied by
                ``HybridSolver(calibrate_body_neff=True)`` rather than by hand;
                see ``FDTDSolver.section_mode_beta()``. None (default) uses the
                mode solver's own n_eff, unchanged.
            cell_neff_override: optional array ``[len(wavelengths), n_cells,
                n_modes]`` replacing each mode's n_eff **per cell**, for the
                REAL multi-cell solve. This is the tapered-body counterpart of
                ``neff_override`` and a different operation, not a variant of
                it:

                * ``neff_override`` replaces the whole section's transfer
                  matrix with a diagonal ``exp(2j·π·n_eff·L/λ)``, skipping the
                  mode solve entirely. Diagonal means NO mode mixing, which is
                  exact for a uniform slab and destroys the physics of a taper,
                  where conversion between modes is the whole point. Hence its
                  uniform-cross-section guard.
                * ``cell_neff_override`` keeps every part of the multi-cell
                  solve — cells, modes, and above all the interface S-matrices
                  — and substitutes n_eff ONLY in the propagation within each
                  cell. The split is exact: interface overlaps depend on the
                  mode FIELDS, propagation depends on n_eff, so calibrating one
                  cannot disturb the other.

                That makes a tapered body calibratable cell by cell, which the
                section-level flag cannot express. Each taper slice IS a
                uniform cell with a well-defined β; what was missing was a way
                to use it without collapsing the section.

                Supply β as n_eff (``β·λ/2π``), e.g. interpolated from a
                β(width) curve measured with
                ``FDTDSolver.section_mode_beta()``. None (default) leaves
                every code path byte-identical to before this parameter
                existed, including the cache fingerprint.
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache — see
                FDTDSolver.solve_section()'s save_dir docstring for the full
                semantics (identical here).
            cache_only: if True, never run the solve — a cache hit returns as
                usual, a miss raises CacheMissError. See
                FDTDSolver.solve_section()'s cache_only docstring for the
                full semantics (identical here).

        Returns:
            Complex ndarray [len(wavelengths), n_modes_out, n_modes_in].
        """
        import meow as mw

        wavelengths = np.asarray(wavelengths)
        n_modes = max(n_modes_in, n_modes_out)
        length = (x1 - x0) + stub_left + stub_right

        # Mutually exclusive by construction: neff_override REPLACES the
        # section with a diagonal and skips the mode solve, while
        # cell_neff_override calibrates propagation INSIDE a real multi-cell
        # solve. Honouring both is not meaningful -- the diagonal path never
        # builds the cells the per-cell array indexes -- so refuse rather than
        # silently let one win.
        if neff_override is not None and cell_neff_override is not None:
            raise ValueError(
                "neff_override and cell_neff_override are different "
                "propagation models, not two settings of one: the former "
                "collapses the section to a diagonal transfer matrix (no mode "
                "mixing, uniform spans only), the latter keeps the multi-cell "
                "solve and its interface S-matrices. Pass exactly one.")

        # A request for fewer modes than a prior call already computed can
        # reuse that larger result by slicing its mode axes — unlike FDTD,
        # EME has no excited/measured-side asymmetry and no monitor geometry
        # to pin, since compute_modes(cs, n_modes) solves the full
        # cross-section once regardless of how the n_modes_in/n_modes_out
        # split is drawn from it, so no extra precondition/invariant is
        # needed beyond dropping the mode counts from the family config
        # (Solver recommendation #2). Caveat: this assumes the top-k slice
        # of an N-mode solve matches a direct k-mode solve, which holds as
        # long as meow's mode solver doesn't reorder near-degenerate modes
        # differently between differently-sized requests — see
        # tests/test_cache_config.py's EME superset-reuse coverage.
        # Included in the fingerprint ONLY when set. Adding a key
        # unconditionally would change EVERY existing solve_section hash and
        # invalidate the whole on-disk cache for every device in this repo --
        # the expensive FDTD sections included. Keeping it conditional means a
        # run that does not use per-cell calibration fingerprints exactly as it
        # did before this parameter existed. tests/test_cache_config.py pins
        # that with a golden hash.
        _cell_fp = ({} if cell_neff_override is None
                    else {"cell_neff_override": _cell_neff_fingerprint(cell_neff_override)})
        family_config = None
        if save_dir:
            family_config = build_run_config(
                self, "solve_section",
                dict(x0=x0, x1=x1, wavelengths=wavelengths,
                     stub_left=stub_left, stub_right=stub_right,
                     geometry_kwargs=geometry_kwargs,
                     neff_override=_neff_fingerprint(neff_override),
                     **_cell_fp,
                     **resolved_material_config(
                         self.material_stack, wavelengths,
                         self._resolved_dispersive())),
                geometry,
            )
            config = build_run_config(
                self, "solve_section",
                dict(x0=x0, x1=x1, wavelengths=wavelengths,
                     n_modes_in=n_modes_in, n_modes_out=n_modes_out,
                     stub_left=stub_left, stub_right=stub_right,
                     geometry_kwargs=geometry_kwargs,
                     neff_override=_neff_fingerprint(neff_override),
                     **_cell_fp,
                     **resolved_material_config(
                         self.material_stack, wavelengths,
                         self._resolved_dispersive())),
                geometry,
            )
            hit = lookup(save_dir, config)
            if hit is not None:
                logger.info("solve_section: cache hit at %s", hit)
                return load_result(hit)
            superset = lookup_superset(save_dir, family_config, (n_modes_out, n_modes_in))
            if superset is not None:
                folder, _stored_dims = superset
                logger.info("solve_section: superset cache hit at %s (sliced)", folder)
                return load_result(folder)[:, :n_modes_out, :n_modes_in]
        if cache_only:
            raise CacheMissError(
                f"solve_section(): cache_only=True but no cached run found for "
                f"this config (save_dir={save_dir!r}). Run once with "
                f"cache_only=False to populate the cache."
            )

        # Fast path: solve_section always builds a single cell (Ls=[length]),
        # so left/right ports share one cross-section's mode basis -- modes
        # are orthogonal by construction, making the transfer matrix diagonal
        # with T[i,i] = exp(2j*pi*neff_i/wl*L). When neff_override replaces
        # neff right after the solve anyway (calibrate_body_neff), that
        # diagonal is determined entirely by neff_override, never by the
        # field profile -- so the eigenmode solve and mw.compute_s_matrix()
        # are skipped outright rather than run and discarded.
        if neff_override is not None:
            # The diagonal form above is a CONSEQUENCE of the single-cell
            # build, not a property of neff_override -- so the gate has to be
            # "is this span uniform?", not "was an override passed?". A
            # tapered span collapsed onto one cross-section's basis returns a
            # confidently wrong, non-diagonal-in-reality answer with no
            # warning. mmi_1x4_tapered_body_validation.md step 1 (option (b));
            # HybridSolver._calibrated_neff() also checks, earlier and with
            # the section index, but this catches direct callers it cannot
            # see (e.g. mmi_2x2_hybrid.py's own hand-rolled cascade).
            # [x0 - stub_left, x1 + stub_right] is the real drawn span --
            # `length` above includes both stubs, so both must be uniform too.
            assert_uniform_cross_section(
                geometry, x0 - stub_left, x1 + stub_right,
                what="EMESolver.solve_section(neff_override=...)",
                why="an overridden-neff section is propagated as a diagonal "
                    "exp(2j*pi*neff/wl*L), which only describes a uniform "
                    "cross-section; a varying one needs a real multi-cell solve",
                geometry_kwargs=geometry_kwargs,
            )
            T = np.zeros((len(wavelengths), n_modes_out, n_modes_in), dtype=complex)
            for wl_idx, wl in enumerate(wavelengths):
                for k in range(min(n_modes_out, n_modes_in)):
                    T[wl_idx, k, k] = np.exp(
                        2j * np.pi * neff_override[wl_idx][k] / wl * length
                    )
        else:
            # Real multi-cell EME (mmi_1x4_tapered_body_validation.md step 3).
            # solve_section() assumes ONE EME region per solver instance -- a
            # HybridSolver cascade already carries one dedicated EMESolver per
            # section, so a multi-region instance reaching here is a caller
            # error (build a separate EMESolver for this section), not a case
            # to silently half-honor by reading only its first region.
            if len(self.num_cells) != 1:
                raise ValueError(
                    f"solve_section() assumes a single-region EMESolver "
                    f"(one EME region per hybrid-cascade section) -- got "
                    f"num_cells={self.num_cells}. Build a dedicated EMESolver "
                    "for this section instead of reusing a whole-device, "
                    "multi-region instance."
                )

            gds_polys = geometry.section_gds(
                x0, x1, self.oxide_width,
                stub_left=stub_left, stub_right=stub_right, **geometry_kwargs
            )
            structs = extrude_polygons(gds_polys, self._extrusion_rules(wavelengths))
            mesh = self._build_mesh()
            # num_cells[0]=1 (every existing device script's value) reproduces
            # today's single cell of Ls=[length] exactly -- stubs folded into
            # the one cell length, unchanged.
            #
            # Per-cell mode count stays `n_modes` (= max(n_modes_in,
            # n_modes_out)), NOT self.num_modes[0], deliberately: `num_modes`
            # is documented as an unused placeholder on this path by several
            # callers, and modal_multiplexer_hybrid.py actively RELIES on that
            # -- it fills three EME slots with ONE body_solver(num_modes=[15]),
            # including a connector section whose interfaces carry CONN_MODES=6,
            # and says so in its own comment. Reading num_modes[0] here would
            # silently solve 15 modes where 6 were solved before: a behaviour
            # change (the top-6 slice of a 15-mode solve is not guaranteed
            # identical to a direct 6-mode solve -- the same caveat this
            # method's own superset-cache comment above already records) and a
            # cost increase, on a device nobody was changing. A caller wanting
            # a richer internal body basis than its interfaces expose should
            # raise n_modes_in/n_modes_out, which it controls directly.
            cell_lengths = self._flat_cell_lengths(region_lengths=[length])
            modes_per_cell = [n_modes] * len(cell_lengths)
            cells = mw.create_cells(structures=structs, mesh=mesh, Ls=cell_lengths)

            if cell_neff_override is not None:
                _cno = np.asarray(cell_neff_override)
                if _cno.shape[:2] != (len(wavelengths), len(cells)):
                    raise ValueError(
                        f"cell_neff_override has shape {_cno.shape}; expected "
                        f"[{len(wavelengths)}, {len(cells)}, >= {n_modes}] "
                        f"(wavelengths, cells, modes). A silent shape mismatch "
                        f"here would calibrate the wrong cell.")
                if _cno.shape[2] < n_modes:
                    raise ValueError(
                        f"cell_neff_override supplies {_cno.shape[2]} modes but "
                        f"this section carries {n_modes}; every propagating "
                        f"mode needs a calibrated n_eff.")

            T = np.zeros((len(wavelengths), n_modes_out, n_modes_in), dtype=complex)

            for wl_idx, wl in enumerate(wavelengths):
                env = mw.Environment(wl=float(wl), T=self.temperature)
                css = [mw.CrossSection.from_cell(cell=cell, env=env) for cell in cells]
                modes = [self.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]
                if cell_neff_override is not None:
                    # Substitute n_eff PER CELL, leaving every field untouched.
                    #
                    # This is what makes a tapered body calibratable. meow
                    # splits cleanly: `compute_interface_s_matrix` integrates
                    # the mode FIELDS, while propagation within a cell is
                    # exp(2j*pi*neff/wl*L). So replacing neff changes only the
                    # propagation and cannot perturb a single interface
                    # overlap -- unlike `neff_override`, which discards the
                    # interfaces altogether by going diagonal.
                    #
                    # Modes are rebuilt rather than mutated because meow's Mode
                    # is immutable; `meow.mode.normalize()` constructs new ones
                    # the same way.
                    # NaN means "leave this mode's solved n_eff alone".
                    #
                    # Needed because mode IDENTITY is not always stable across
                    # the widths a calibration curve is measured on: the
                    # eigensolver sorts by n_eff, so two modes that cross as
                    # the width changes swap indices, and an interpolated beta
                    # for such a mode blends two different physical modes.
                    # Measured on the waisted 2x2's beta(width) curve -- modes
                    # 6/7 swap between w=2.66 and 2.73 um. Those modes carried
                    # 0.16% of the power, so calibrating the well-ordered ones
                    # and leaving the rest on the solver's own n_eff is both
                    # safer and near-lossless. All-or-nothing would have forced
                    # a choice between a known-wrong beta and no calibration.
                    def _sub(m, value):
                        if np.isnan(np.real(value)):
                            return m
                        return type(m)(neff=complex(value), cs=m.cs,
                                       Ex=m.Ex, Ey=m.Ey, Ez=m.Ez,
                                       Hx=m.Hx, Hy=m.Hy, Hz=m.Hz)

                    modes = [
                        [_sub(m, _cno[wl_idx, ci, mi])
                         for mi, m in enumerate(cell_modes)]
                        for ci, cell_modes in enumerate(modes)
                    ]
                S, port_map = mw.compute_s_matrix(modes, cells=cells)
                # S is a flat dense matrix keyed by string port labels ("left@n",
                # "right@n"), not an NxN block per port pair — look up indices
                # explicitly rather than assuming S[0,1] is the transmission block.
                # These labels come only from the FIRST cell's left faces and
                # the LAST cell's right faces regardless of how many cells sit
                # in between (meow/eme/cascade.py:_get_netlist) -- correct for
                # N cells unchanged from the single-cell case.
                for i in range(n_modes_out):
                    for j in range(n_modes_in):
                        i_out = port_map[f"right@{i}"]
                        j_in = port_map[f"left@{j}"]
                        T[wl_idx, i, j] = S[i_out, j_in]

        # Post-hoc duplicate-mode scan of the transfer matrix just computed —
        # fdtd_vs_EME_mode_validation.md item 8. `compute_modes()`'s own
        # guard_duplicate_modes already rejects a duplicate WITHIN the mode
        # basis this section solves on (see its docstring); this instead
        # checks the transfer matrix HybridSolver actually contracts, which
        # is a different (if related) object -- free, since T is already
        # computed.
        from ..mode_identity import warn_on_transfer_duplicate

        warn_on_transfer_duplicate(T, logger, x0, x1)

        if save_dir:
            folder = save_run(save_dir, config, T)
            logger.info("solve_section: cache saved at %s", folder)
            update_family_pointer(
                save_dir, family_config, (n_modes_out, n_modes_in), os.path.basename(folder),
            )

        return T

    def solve_section_length_sweep(
        self,
        geometry,
        wavelengths: np.ndarray,
        x0: float,
        x1: float,
        lengths: np.ndarray,
        n_modes_in: int = 1,
        n_modes_out: int = 1,
        stub_left: float = 0.0,
        stub_right: float = 0.0,
        neff_override: Optional[np.ndarray] = None,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        **geometry_kwargs,
    ) -> np.ndarray:
        """Modal transfer matrices for [x0, x1] swept over many section
        lengths, from a single eigenmode solve per wavelength.

        [x0, x1] must be a uniform-cross-section region (the actual cell
        length used for meshing/mode-solving is irrelevant here — only the
        cross-section shape matters, and that doesn't change with length).
        Mode profiles and effective indices depend only on the cross-section
        and wavelength (meow's CrossSection.from_cell() never reads cell
        length), so for a single uniform cell the propagation S-matrix is
        exactly diagonal in the mode basis: S_ii(L) = exp(2j*pi*neff_i/wl*L).
        meow's compute_s_matrix() accepts an explicit cell_lengths override
        that swaps in a different L for the same already-computed modes, with
        no new mode solve — so this method costs exactly the same as one
        solve_section() call (one mode solve per wavelength), regardless of
        how many lengths are swept.

        Caching has THREE tiers, cheapest first (see beef/cache/ for the
        underlying lookup()/lookup_superset() primitives):
          1. Exact match on every arg including `lengths` — unchanged from
             before, a plain on-disk hit.
          2. Mode-count superset at the SAME `lengths` — unchanged from
             before (`family_config`/lookup_superset()): a cached run with
             >= the requested mode count, at this exact lengths array, gets
             sliced down.
          3. NEW — length-independent neff reuse (`neff_config`): per the
             diagonal-in-mode-basis fact above, the only genuinely expensive
             step is the eigenmode solve itself (one per wavelength), which
             never reads `lengths` at all. So a *previous* call at this same
             geometry/wavelengths/mode-count — with ANY lengths array, swept
             or not — already solved and can hand over its per-wavelength
             neff; this tier reconstructs T for the requested `lengths`
             purely analytically (_propagate_diagonal(), no meow call) and
             is what actually eliminates a re-solve when only the sweep
             range changes.

             Every fresh (non-cached) solve below validates this reuse path
             is safe for the current cross-section BEFORE trusting it for
             future calls: it reconstructs one length's T from the neff it
             just solved and compares against the real mw.compute_s_matrix()
             result for that same length (_DIAGONAL_RECONSTRUCTION_TOL). Only
             a passing check gets its neff saved to tier 3's cache — a
             cross-section this fact somehow doesn't hold for (it always has,
             so far) just never becomes eligible for tier-3 reuse, and still
             returns its own correct (real-solve) T regardless.

        Args:
            neff_override: optional array ``[len(wavelengths), n_modes]``
                replacing each mode's n_eff before propagation, leaving field
                profiles — and so every interface overlap — untouched. Lets a
                HybridSolver body propagate on the β its FDTD sections
                actually accumulate rather than on MPB's mode index: the two
                differ by ~2.3% per mode at res 40, and because Δβ is a
                difference of two β near 9.4 rad/µm that lands as a ~10% error
                on the coupling (dc_validation.md item 6). Normally supplied by
                ``HybridSolver(calibrate_body_neff=True)`` rather than by hand;
                see ``FDTDSolver.section_mode_beta()``. None (default) uses the
                mode solver's own n_eff, unchanged. Included as part of every
                cache tier's key (including tier 3's neff_config) — a
                calibrated and an uncalibrated run never share a cache entry.
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache — see
                FDTDSolver.solve_section()'s save_dir docstring for the full
                semantics (identical here; lengths is included in tiers 1/2's
                config, deliberately excluded from tier 3's neff_config).
            cache_only: if True, never run the solve — a cache hit (any tier)
                returns as usual, a miss at all three raises CacheMissError.
                See FDTDSolver.solve_section()'s cache_only docstring for the
                full semantics (identical here).

        Returns:
            Complex ndarray [len(lengths), len(wavelengths), n_modes_out, n_modes_in].
        """
        import meow as mw

        wavelengths = np.asarray(wavelengths)
        lengths = np.asarray(lengths)
        n_modes = max(n_modes_in, n_modes_out)

        # See solve_section()'s comment on family_config/lookup_superset —
        # same reuse mechanism, same no-extra-precondition reasoning; lengths
        # stays part of the family key (tier 2 only drops mode counts, same
        # as before — Solver recommendation #2). Tier 3's neff_config below
        # additionally drops lengths.
        family_config = None
        neff_config = None
        if save_dir:
            _material_cfg = resolved_material_config(
                self.material_stack, wavelengths, self._resolved_dispersive())
            family_config = build_run_config(
                self, "solve_section_length_sweep",
                dict(x0=x0, x1=x1, wavelengths=wavelengths, lengths=lengths,
                     stub_left=stub_left, stub_right=stub_right,
                     geometry_kwargs=geometry_kwargs,
                     neff_override=_neff_fingerprint(neff_override),
                     **_material_cfg),
                geometry,
            )
            config = build_run_config(
                self, "solve_section_length_sweep",
                dict(x0=x0, x1=x1, wavelengths=wavelengths, lengths=lengths,
                     n_modes_in=n_modes_in, n_modes_out=n_modes_out,
                     stub_left=stub_left, stub_right=stub_right,
                     geometry_kwargs=geometry_kwargs,
                     neff_override=_neff_fingerprint(neff_override),
                     **_material_cfg),
                geometry,
            )
            hit = lookup(save_dir, config)
            if hit is not None:
                logger.info("solve_section_length_sweep: cache hit at %s", hit)
                return load_result(hit)
            superset = lookup_superset(save_dir, family_config, (n_modes_out, n_modes_in))
            if superset is not None:
                folder, _stored_dims = superset
                logger.info("solve_section_length_sweep: superset cache hit at %s (sliced)", folder)
                return load_result(folder)[:, :, :n_modes_out, :n_modes_in]

            # Tier 3: same ingredients as family_config, with `lengths` ALSO
            # dropped and `n_modes` pinned explicitly (family_config has no
            # mode-count field at all, so a distinct method_name plus an
            # explicit n_modes keeps this tier's key from colliding with it).
            neff_config = build_run_config(
                self, "solve_section_neff",
                dict(x0=x0, x1=x1, wavelengths=wavelengths, n_modes=n_modes,
                     stub_left=stub_left, stub_right=stub_right,
                     geometry_kwargs=geometry_kwargs,
                     neff_override=_neff_fingerprint(neff_override),
                     **_material_cfg),
                geometry,
            )
            neff_hit = lookup(save_dir, neff_config)
            if neff_hit is not None:
                neff_all = load_result(neff_hit)  # [n_wl, n_modes] complex, self-checked when saved
                logger.info(
                    "solve_section_length_sweep: neff cache hit at %s -- "
                    "%d length(s) reconstructed analytically, no new mode solve",
                    neff_hit, len(lengths),
                )
                # Deliberately NOT persisted. This T was rebuilt from neff_all
                # with no meow call at all, so writing it back would cache a
                # free result -- and cache it under a key containing `lengths`
                # verbatim, which is what pinned callers to one length grid and
                # produced a 342 MB entry for a single 7601-length sweep.
                return _propagate_diagonal(neff_all, wavelengths, lengths, n_modes_out, n_modes_in)
        if cache_only:
            raise CacheMissError(
                f"solve_section_length_sweep(): cache_only=True but no cached "
                f"run found for this config (save_dir={save_dir!r}). Run once "
                f"with cache_only=False to populate the cache."
            )

        # Fast path: with neff_override, the diagonal-in-mode-basis fact this
        # method's own docstring states (S_ii(L) = exp(2j*pi*neff_i/wl*L))
        # applies directly to the SUPPLIED neff, so build T with the same
        # _propagate_diagonal() tier 3 already trusts (validated against a
        # real solve every time it populates that cache, see
        # _DIAGONAL_RECONSTRUCTION_TOL above) instead of solving for modes
        # whose neff gets overwritten immediately after anyway — same
        # argument as solve_section()'s fast path, extended here to the
        # length sweep.
        if neff_override is not None:
            neff_all = np.asarray(neff_override, dtype=complex)
            T = _propagate_diagonal(neff_all, wavelengths, lengths, n_modes_out, n_modes_in)
            neff_cache_safe = True  # exact by construction -- no real solve to check against
        else:
            gds_polys = geometry.section_gds(
                x0, x1, self.oxide_width,
                stub_left=stub_left, stub_right=stub_right, **geometry_kwargs
            )
            structs = extrude_polygons(gds_polys, self._extrusion_rules(wavelengths))
            mesh = self._build_mesh()
            cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[1.0])  # length is a placeholder

            T = np.zeros((len(lengths), len(wavelengths), n_modes_out, n_modes_in), dtype=complex)
            neff_all = np.zeros((len(wavelengths), n_modes), dtype=complex)  # tier 3's cache payload

            for wl_idx, wl in enumerate(wavelengths):
                env = mw.Environment(wl=float(wl), T=self.temperature)
                css = [mw.CrossSection.from_cell(cell=cell, env=env) for cell in cells]
                modes = [self.compute_modes(cs, n_modes) for cs in css]  # one solve, reused below

                neff_all[wl_idx, :] = np.array([m.neff for m in modes[0]], dtype=complex)

                if wl_idx == len(wavelengths) // 2:
                    # Diagnostic: confirm this cross-section actually resolves
                    # n_modes distinct, non-degenerate supermodes (neff split
                    # apart) rather than a near-degenerate/duplicate pair, which
                    # would explain a single-arm launch failing to couple into
                    # one of them.
                    neffs = [m.neff for m in modes[0]]
                    logger.debug("EME cross-section neff @ λ=%.3f µm: %s", wl, neffs)

                for L_idx, L in enumerate(lengths):
                    S, port_map = mw.compute_s_matrix(modes, cells=cells, cell_lengths=[float(L)])
                    for i in range(n_modes_out):
                        for j in range(n_modes_in):
                            i_out = port_map[f"right@{i}"]
                            j_in = port_map[f"left@{j}"]
                            T[L_idx, wl_idx, i, j] = S[i_out, j_in]

            # Tier 3's own gate: reconstruct purely from neff_all and compare
            # against the real solve above. Only a pass makes neff_all eligible
            # to serve future different-lengths calls (see docstring).
            #
            # ⚠ CHECKS EVERY SWEPT LENGTH, not just the first. It used to test
            # `lengths[:1]` alone, which is the SHORTEST length and therefore
            # the most forgiving: reconstruction error accumulates with
            # propagation distance (it enters as exp(2j*pi*neff*L/wl), so a
            # finite-precision neff drifts linearly in L). MEASURED on the
            # existing cache: one real 7601-length entry
            # (EMESolver_solve_section_length_sweep_MMI_fbdfbe7f02c20a01) sat
            # at 7.2e-08 at L=19.5 um -- a comfortable pass -- while drifting
            # to 1.48e-06 at L=399.5 um, with 34% of its grid over tolerance.
            # The other 48 matched entries in that cache reconstruct to
            # 4.5e-13 or better across their whole grid, several bit-exactly,
            # so this is a real outlier being caught rather than a tolerance
            # that is generally too tight. Since this method no longer
            # persists T, a gate that passes on the shortest length would make
            # the drifted reconstruction the only answer available.
            #
            # Chunked so peak memory stays bounded: T can be O(100s of MB)
            # here, and reconstructing a full second copy to compare would
            # double that.
            neff_cache_safe = True
            if save_dir and neff_config is not None and len(lengths) > 0:
                worst_diff, worst_L = 0.0, float(lengths[0])
                for lo in range(0, len(lengths), _DIAGONAL_CHECK_CHUNK):
                    hi = min(lo + _DIAGONAL_CHECK_CHUNK, len(lengths))
                    T_check = _propagate_diagonal(
                        neff_all, wavelengths, lengths[lo:hi], n_modes_out, n_modes_in)
                    diffs = np.max(np.abs(T_check - T[lo:hi]), axis=(1, 2, 3))
                    k = int(np.argmax(diffs))
                    if float(diffs[k]) > worst_diff:
                        worst_diff, worst_L = float(diffs[k]), float(lengths[lo + k])
                if worst_diff > _DIAGONAL_RECONSTRUCTION_TOL:
                    neff_cache_safe = False
                    logger.warning(
                        "solve_section_length_sweep: diagonal-reconstruction self-check "
                        "off by %.2e (tol %.1e) at length=%.6f um of %d swept -- this "
                        "cross-section gets no length-independent neff entry, so its "
                        "derived transfer matrix is persisted instead.",
                        worst_diff, _DIAGONAL_RECONSTRUCTION_TOL, worst_L, len(lengths),
                    )
                else:
                    logger.info(
                        "solve_section_length_sweep: diagonal reconstruction verified "
                        "across all %d swept lengths, worst %.2e at L=%.6f um (tol %.1e)",
                        len(lengths), worst_diff, worst_L, _DIAGONAL_RECONSTRUCTION_TOL,
                    )

        if save_dir:
            # Persist ONLY the length-independent primitive. T is a derived
            # expansion: this section is single-cell by construction (Ls=[1.0]
            # above), so its propagation S-matrix is exactly diagonal in the
            # mode basis and _propagate_diagonal() rebuilds T at ANY length
            # from neff_all alone -- checked against this very solve by the
            # _DIAGONAL_RECONSTRUCTION_TOL gate above. Persisting T instead
            # cost 342 MB for one 7601-length grid and, because `lengths` is
            # hashed verbatim, made every other grid a full miss.
            if neff_config is not None and neff_cache_safe:
                folder = save_run(save_dir, neff_config, neff_all)
                logger.info(
                    "solve_section_length_sweep: neff primitive saved at %s "
                    "(%d length(s) served analytically, none persisted)",
                    folder, len(lengths),
                )
            else:
                # The diagonal reconstruction is unavailable for this
                # cross-section, so nothing can rebuild T from neff -- persist
                # the derived matrix so a real solve is not lost. Rare, and
                # each occurrence is worth investigating rather than ignoring.
                folder = save_run(save_dir, config, T)
                update_family_pointer(
                    save_dir, family_config, (n_modes_out, n_modes_in), os.path.basename(folder),
                )
                logger.warning(
                    "solve_section_length_sweep: no length-independent neff primitive "
                    "available (neff_cache_safe=False) -- persisted the derived %s "
                    "transfer matrix at %s instead", T.shape, folder,
                )

        return T

    def solve(self, geometry: DeviceGeometry, wavelengths: np.ndarray,
              ports: Optional[List[str]] = None, save_dir: Optional[str] = None,
              cache_only: bool = False,
              **kwargs) -> SMatrix:
        """Run EME simulation across all wavelengths.

        Args:
            geometry: device geometry implementing eme_lengths().
            wavelengths: [M] wavelength array (µm).
            ports: list of port sides to simulate.
                   Defaults to ["top", "bot"] for 2-output devices.
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache — see
                solve_section()'s save_dir docstring for the full semantics
                (identical here, fingerprinted via geometry.to_gds_eme()
                instead of section_gds()). Cached per port side, since each
                side builds its own geometry/cells.
            cache_only: if True, never run the solve — a cache hit for a
                given port returns as usual, a miss raises CacheMissError
                immediately (fail-fast — does not check the remaining
                ports first). See solve_section()'s cache_only docstring for
                the full semantics.
            **kwargs: forwarded to geometry.to_gds_eme() (e.g. output_idx to
                      restrict the drawn output branch to a single port).

        Returns:
            SMatrix with entries keyed by port name.
        """
        import meow as mw

        if ports is None:
            ports = ["top", "bot"]

        wavelengths = np.asarray(wavelengths)
        smat = SMatrix(wavelengths=wavelengths)
        modes_per_cell = self._modes_per_cell()

        for side in ports:
            config = None
            if save_dir:
                gds_polys = geometry.to_gds_eme(
                    straight_length=self.straight_length, oxide_width=self.oxide_width,
                    x_start=0.0, side=side, **kwargs,
                )
                # Added 2026-08-03 (`mmi_1x4_sin300_validation.md` step 3).
                # Before this, solve()'s fingerprint carried NO material
                # information beyond `_solver_config()`'s bare
                # `material_stack: {wg_height, sz}` -- unlike solve_section()/
                # solve_section_length_sweep(), which have called
                # resolved_material_config() since 2026-08-01. That means
                # solve() could not distinguish two stacks agreeing on
                # wg_height/sz but differing in core_kind (silicon vs nitride
                # at the same thickness) OR in `dispersive` (frozen vs
                # dispersive DOES change what _build_cells()/_extrusion_rules()
                # actually draws here -- see _extrusion_rules()'s own
                # docstring) -- a real collision risk for the *_eme.py
                # examples this method serves, not a theoretical one.
                config = build_whole_device_config(
                    self, "solve",
                    dict(side=side, wavelengths=wavelengths, kwargs=kwargs,
                         **resolved_material_config(
                             self.material_stack, wavelengths,
                             self._resolved_dispersive())),
                    geometry, gds_polys,
                )
                hit = lookup(save_dir, config)
                if hit is not None:
                    logger.info("solve: cache hit at %s", hit)
                    smat.set(load_result(hit), port=side)
                    continue
            if cache_only:
                raise CacheMissError(
                    f"solve(): cache_only=True but no cached run found for "
                    f"port={side!r} (save_dir={save_dir!r}). Run once with "
                    f"cache_only=False to populate the cache."
                )

            cells = self._build_cells(geometry, side=side, wavelengths=wavelengths, **kwargs)
            amplitudes = np.zeros(len(wavelengths), dtype=complex)

            for wl_idx, wl in enumerate(wavelengths):
                env = mw.Environment(wl=float(wl), T=self.temperature)
                css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
                modes = [self.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]
                S, _ = mw.compute_s_matrix(modes, cells=cells)
                amplitudes[wl_idx] = S[0, 1]

            if save_dir:
                folder = save_run(save_dir, config, amplitudes)
                logger.info("solve: cache saved at %s", folder)

            smat.set(amplitudes, port=side)

        return smat

    def _amps_from_factorizations(self, facs, lengths, ref_cell_lengths,
                                  c_start, c_stop, n_region_cells):
        """Cascade cached factorizations over a length grid, no mode solve.

        `facs` is one factorization per wavelength, in wavelength order, exactly
        as solve_length_sweep() stores them. The swept region's cells take
        length/n_region_cells each; every other cell keeps its reference length.
        """
        amps = np.zeros((len(lengths), len(facs)), dtype=complex)
        for wl_idx, fac in enumerate(facs):
            for k, L in enumerate(lengths):
                cell_lengths = list(ref_cell_lengths)
                for j in range(c_start, c_stop):
                    cell_lengths[j] = float(L) / n_region_cells
                S, _ = s_matrix_from_factorization(fac, cell_lengths)
                amps[k, wl_idx] = S[0, 1]
        return amps

    def solve_length_sweep(
        self,
        geometry: DeviceGeometry,
        wavelengths: np.ndarray,
        lengths: np.ndarray,
        ports: Optional[List[str]] = None,
        sweep_region: Optional[int] = None,
        save_dir: Optional[str] = None,
        cache_only: bool = False,
        **kwargs,
    ) -> List[SMatrix]:
        """Whole-device S-matrices for many device lengths, from ONE mode
        solve per wavelength — the mode basis is reused across every length.

        The cross-section *shape* an EME cell samples depends only on the
        taper fraction (i+0.5)/num_cells the solver evaluates, never on a
        region's physical length: a linear taper at fraction t has the same
        width regardless of how long the taper is (only Δx-per-cell scales).
        meow's S-matrix factors cleanly into an interface part (depends only
        on the modes) and a propagation part (the only place a cell's
        physical length enters — as exp(2j·π·neff/wl·L_cell)); see
        meow.eme.compute_propagation_s_matrix. So sweeping one region's
        length changes only a cell_lengths list handed to the cascade, never
        the modes or interfaces — turning an N-length sweep from N mode
        solves per wavelength into one. This is the multi-cell generalization
        of solve_section_length_sweep()'s single-uniform-cell trick.

        Only a region whose cross-section scales with its length (a taper or a
        straight) may be swept. sweep_region selects which self.lengths entry
        varies; it defaults to the unique region with num_cells > 1 (the
        subdivided taper). num_cells for that region is held FIXED across the
        sweep — it's the cross-section-staircase resolution (an accuracy dial,
        checkable by raising it), independent of the physical length, which
        only sets each cell's propagation length lengths[k]/num_cells.

        The geometry is drawn once, at whatever reference length it was built
        with, purely to sample those (length-independent) cross-section
        shapes — its swept-region length does not affect the result; only the
        values in `lengths` do. Every other (non-swept) region keeps its
        self.lengths value, so those must already match the real device
        (e.g. fixed access-waveguide stubs).

        Args:
            geometry: device geometry, built at any reference length for the
                swept region.
            wavelengths: [M] wavelength array (µm).
            lengths: [K] swept-region lengths to evaluate (µm).
            ports: port sides to simulate (as in solve()); defaults to
                ["top", "bot"].
            sweep_region: index into self.lengths that varies (default: the
                unique num_cells > 1 region — see _infer_sweep_region()).
            save_dir: persistent on-disk cache directory (see beef/cache/).
                None/"" (default) disables the cache. The full [K, M] result
                array is cached per port side, keyed on the solver config, the
                reference geometry, `lengths`, and sweep_region.
            cache_only: if True, never run the solve — a cache hit returns as
                usual, a miss raises CacheMissError (fail-fast per port).
            **kwargs: forwarded to geometry.to_gds_eme() (e.g. output_idx).

        Returns:
            list[SMatrix] aligned with `lengths`: result[k] is the SMatrix for
            device length lengths[k], keyed by port name (same convention as
            solve()).
        """
        import meow as mw

        if ports is None:
            ports = ["top", "bot"]
        if sweep_region is None:
            sweep_region = self._infer_sweep_region()

        wavelengths = np.asarray(wavelengths)
        lengths = np.asarray(lengths)
        modes_per_cell = self._modes_per_cell()
        ref_cell_lengths = self._flat_cell_lengths()
        c_start, c_stop = self._region_cell_span(sweep_region)
        n_region_cells = c_stop - c_start

        results = [SMatrix(wavelengths=wavelengths) for _ in lengths]

        for side in ports:
            config = None
            if save_dir:
                gds_polys = geometry.to_gds_eme(
                    straight_length=self.straight_length, oxide_width=self.oxide_width,
                    x_start=0.0, side=side, **kwargs,
                )
                # Same fix and reasoning as solve() above, 2026-08-03
                # (`mmi_1x4_sin300_validation.md` step 3).
                config = build_whole_device_config(
                    self, "solve_length_sweep",
                    dict(side=side, wavelengths=wavelengths, lengths=lengths,
                         sweep_region=sweep_region, kwargs=kwargs,
                         **resolved_material_config(
                             self.material_stack, wavelengths,
                             self._resolved_dispersive())),
                    geometry, gds_polys,
                )
                hit = lookup(save_dir, config)
                if hit is not None:
                    logger.info("solve_length_sweep: cache hit at %s", hit)
                    amps = load_result(hit)
                    for k in range(len(lengths)):
                        results[k].set(amps[k], port=side)
                    continue
                # Length-independent tier: the same ingredients with `lengths`
                # DROPPED. What is stored is the cascade's factorization -- the
                # interface S-matrices (a function of the modes alone) plus each
                # cell's neff -- which replays mw.compute_s_matrix() at ANY length
                # grid, so one entry serves every sweep of this device instead of
                # one entry per grid. Verified bit-exact against a direct cascade.
                fac_config = build_whole_device_config(
                    self, "solve_length_sweep_factorization",
                    dict(side=side, wavelengths=wavelengths,
                         sweep_region=sweep_region, kwargs=kwargs,
                         **resolved_material_config(
                             self.material_stack, wavelengths,
                             self._resolved_dispersive())),
                    geometry, gds_polys,
                )
                fac_hit = lookup(save_dir, fac_config)
                if fac_hit is not None:
                    logger.info(
                        "solve_length_sweep: factorization cache hit at %s -- "
                        "%d length(s) cascaded with no mode solve", fac_hit, len(lengths),
                    )
                    amps = self._amps_from_factorizations(
                        load_result(fac_hit), lengths, ref_cell_lengths,
                        c_start, c_stop, n_region_cells,
                    )
                    for k in range(len(lengths)):
                        results[k].set(amps[k], port=side)
                    continue
            if cache_only:
                raise CacheMissError(
                    f"solve_length_sweep(): cache_only=True but no cached run "
                    f"found for port={side!r} (save_dir={save_dir!r}). Run once "
                    f"with cache_only=False to populate the cache."
                )

            cells = self._build_cells(geometry, side=side, wavelengths=wavelengths, **kwargs)
            amps = np.zeros((len(lengths), len(wavelengths)), dtype=complex)

            def _cell_lengths_for(L):
                cl = list(ref_cell_lengths)
                for j in range(c_start, c_stop):
                    cl[j] = float(L) / n_region_cells
                return cl

            # Sweep THROUGH the factorization rather than calling
            # mw.compute_s_matrix() per length. Both give the same numbers --
            # measured bit-exact -- but the direct call recomputes every
            # interface overlap from the modes on EVERY length, and the
            # interfaces do not depend on length at all. Computing them once
            # per wavelength instead of once per length is the difference
            # between one cascade and len(lengths) of them: mmi_1x4_eme's
            # 101-length sweep was running ~8 h per output branch this way.
            facs = []
            probe = sorted({0, len(lengths) // 2, len(lengths) - 1}) if len(lengths) else []
            gate_worst = 0.0
            for wl_idx, wl in enumerate(wavelengths):
                env = mw.Environment(wl=float(wl), T=self.temperature)
                css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
                modes = [self.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]
                fac = extract_factorization(modes)
                facs.append(fac)
                for k, L in enumerate(lengths):
                    S, _ = s_matrix_from_factorization(fac, _cell_lengths_for(L))
                    amps[k, wl_idx] = S[0, 1]
                # The factorization now produces the answer, so the gate has to
                # be an INDEPENDENT computation: a direct meow cascade at a few
                # lengths, which is what the stored factorization is claiming to
                # reproduce. Comparing the factorization against itself would
                # verify nothing.
                for k in probe:
                    S_ref, _ = mw.compute_s_matrix(
                        modes, cells=cells, cell_lengths=_cell_lengths_for(lengths[k]))
                    gate_worst = max(gate_worst, float(abs(S_ref[0, 1] - amps[k, wl_idx])))

            if save_dir:
                # Gate, in the same spirit as solve_section_length_sweep()'s
                # diagonal-reconstruction check: replay the factorization about
                # to be stored and confirm it reproduces the amplitudes this real
                # solve just produced. Only a pass earns a cache entry, so a
                # factorization never becomes an unverified answer.
                #
                # gate_worst was measured above, against a direct meow cascade
                # at `probe` lengths per wavelength. Sampled rather than
                # exhaustive: the factorization reuses the same interface
                # matrices and the same neff, so a length is either bit-exact or
                # the factorization is broken outright -- error cannot
                # accumulate with L the way the diagonal reconstruction's does,
                # which is why THAT gate checks every length and this one does
                # not. An exhaustive check here would also cost a second full
                # sweep.
                worst = gate_worst
                if fac_config is not None and worst <= _FACTORIZATION_REPLAY_TOL:
                    folder = save_run(save_dir, fac_config, facs)
                    logger.info(
                        "solve_length_sweep: factorization saved at %s (replay "
                        "verified to %.2e at %d of %d length(s); the derived "
                        "amplitudes are not persisted)",
                        folder, worst, len(probe), len(lengths),
                    )
                else:
                    folder = save_run(save_dir, config, amps)
                    logger.warning(
                        "solve_length_sweep: factorization replay off by %.2e "
                        "(tol %.1e) -- persisted the derived amplitudes at %s "
                        "instead of a length-independent factorization",
                        worst, _FACTORIZATION_REPLAY_TOL, folder,
                    )

            for k in range(len(lengths)):
                results[k].set(amps[k], port=side)

        return results

    def solve_sweep(self, geometry: DeviceGeometry, wavelengths: np.ndarray,
                    ports: Optional[List[str]] = None,
                    wl_batch_size: int = 1, **kwargs) -> SMatrix:
        """Memory-efficient EME sweep using batched wavelength processing.

        The center region (index len(num_modes)//2) is treated as the
        variable region — modes are freed after each batch while left/right
        cell modes are pre-computed and reused.
        """
        import meow as mw
        from meow import compute_modes_tidy3d as _compute_modes

        if ports is None:
            ports = ["top", "bot"]

        wavelengths = np.asarray(wavelengths)
        M = len(wavelengths)
        smat = SMatrix(wavelengths=wavelengths)
        modes_per_cell = self._modes_per_cell()

        # Identify center region boundaries in the flat cell list
        center_region = len(self.num_modes) // 2
        cumcells = [0] + list(np.cumsum(self.num_cells))
        split_left = cumcells[center_region]
        split_right = cumcells[center_region + 1]

        cells_per_port = {side: self._build_cells(geometry, side=side,
                                                  wavelengths=wavelengths)
                          for side in ports}
        left_mode_counts = modes_per_cell[:split_left]
        center_mode_count = modes_per_cell[split_left]
        right_mode_counts = modes_per_cell[split_right:]

        left_cells = cells_per_port[ports[0]][:split_left]
        center_cell = cells_per_port[ports[0]][split_left]
        right_cells_per_port = {side: cells_per_port[side][split_right:] for side in ports}

        # Pre-compute left and right modes (constant across sweep)
        left_modes: dict = {}
        right_modes_per_port: dict = {side: {} for side in ports}

        for wl_idx, wl in enumerate(wavelengths):
            env = mw.Environment(wl=float(wl), T=self.temperature)
            left_modes[wl_idx] = [
                _compute_modes(mw.CrossSection.from_cell(cell=c, env=env), nm)
                for c, nm in zip(left_cells, left_mode_counts)
            ]
            for side in ports:
                right_modes_per_port[side][wl_idx] = [
                    _compute_modes(mw.CrossSection.from_cell(cell=c, env=env), nm)
                    for c, nm in zip(right_cells_per_port[side], right_mode_counts)
                ]

        all_wl = list(range(M))
        batches = [all_wl[i:i+wl_batch_size] for i in range(0, M, wl_batch_size)]
        acc: dict = {side: np.zeros(M, dtype=complex) for side in ports}

        for batch in batches:
            center_modes_batch = {}
            for wl_idx in batch:
                env = mw.Environment(wl=float(wavelengths[wl_idx]), T=self.temperature)
                cs = mw.CrossSection.from_cell(cell=center_cell, env=env)
                center_modes_batch[wl_idx] = [_compute_modes(cs, center_mode_count)]

            for wl_idx in batch:
                for side in ports:
                    full_modes = (left_modes[wl_idx]
                                  + center_modes_batch[wl_idx]
                                  + right_modes_per_port[side][wl_idx])
                    S, _ = mw.compute_s_matrix(full_modes, cells=cells_per_port[side])
                    acc[side][wl_idx] = S[0, 1]

            del center_modes_batch

        for side in ports:
            smat.set(acc[side], port=side)

        return smat
