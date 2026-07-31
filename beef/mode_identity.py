"""Field-overlap mode correspondence: a device- and resolution-agnostic
alternative to nearest-neighbour n_eff matching.

``HybridSolver`` contracts index *i* of an FDTD interface basis against
index *i* of an EME body basis, assuming they name the same mode. The older
check for that (``tests/test_mode_identity.py``'s
``assert_modes_correspond()`` — re-scoped 2026-07-30 to a fast smoke check,
superseded by this module as the authoritative mechanism; see its own
docstring's "STATUS" note) compares n_eff and picks each mode's nearest
neighbour in the other engine's n_eff list. That works when the modes are
well separated, but it is fundamentally a 1-D scalar heuristic: on a wide
multimode cross-section several higher-order modes can crowd together in
n_eff closely enough that nearest-neighbour picks the wrong pairing outright
— observed directly on the modal multiplexer's 16-mode basis
(``validation_markdown/modal_multiplexer_validation.md`` gap 1), where meep
bands 8, 9 and 10 all nearest-matched to the SAME meow band. It is also only
meaningful at matched resolution (a resolution change shifts n_eff by more
than the mode spacing, faking a reorder — see ``test_mode_identity.py``'s own
warning about this).

This module compares mode *shape* instead: the normalized transverse
E-field overlap integral between two mode field profiles, matched by the
assignment that maximizes total overlap (not nearest-neighbour in a single
scalar). Two different physical modes are essentially orthogonal even when
their n_eff coincide; two grids' discretizations of the SAME physical mode
overlap near 1 even when their n_eff differ slightly. That makes this
comparison robust to exactly the failure nearest-neighbour n_eff matching
has, and — because it needs no symmetry assumption and works given only
field samples on a shared grid — it applies to any device cross-section at
any resolution, not just symmetric ones or matched meshes.

**E-field only, not the full E-x-H reciprocity overlap — found the hard
way.** An earlier version of this module used the standard reciprocity
overlap kappa_ab = integral(Ea_t x Hb_t*).xhat dA, mixing engine A's E with
engine B's H. That is only meaningful if both engines report E and H in a
mutually consistent unit system. They do not: verified directly (a bare
450x220nm strip, single mode, otherwise-confirmed-matching field shapes)
that meep's and meow's own internal E:H amplitude ratio differ by a factor
of ~375 (meep's raw Ey/Hz ~0.3, meow's raw Ey/Hz ~112 for the SAME mode) —
each engine is internally self-consistent (meep-vs-itself and
meow-vs-itself both gave a clean diag~1, off-diag~0 in this module's own
matrix), but crossing engine A's E with engine B's H directly divides out
to sqrt(Z_a/Z_b) times the true overlap, not the true overlap itself. That
is a per-engine constant unrelated to whether the modes actually match, and
it drove the diagonal down to ~0.05 across every device/resolution tried
until this was found and fixed. Comparing E against E only (never crossing
E from one engine with H from the other) sidesteps the whole issue: an
engine's OWN E-field absolute scale is whatever it is, but it never needs
to be compared against the other engine's H scale at all.

Deliberately independent of both meep and meow at import time: this module
only consumes plain field arrays (``{'Ex','Ey','Ez','Hx','Hy','Hz'}``, each
``[n_modes, ny, nz]`` complex) and a uniform grid spacing, so it can be unit
tested with synthetic data and used with either engine's output via
``beef.solvers.fdtd.FDTDSolver.section_mode_fields()`` (meep) or
``stack_meow_mode_fields()`` below (meow).
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment

FieldSet = Dict[str, np.ndarray]

# Self-overlap classification thresholds -- see self_overlap_diagnostics().
# _DUPLICATE_OVERLAP matches fdtd.py's _DUP_BAND_OVERLAP_THRESHOLD and
# eme.py's _DUP_MODE_OVERLAP_THRESHOLD: one number for "these are the same
# eigenvector" across all three call sites.
_DUPLICATE_OVERLAP = 0.9
_ORTHOGONALITY_TOL = 0.1


def mesh_cell_centers(mesh) -> Tuple[np.ndarray, np.ndarray]:
    """(y, z) sample points matching a meow ``Mesh2D``'s field arrays.

    meow ``Mode`` field arrays (``Ex``, ``Hy``, ...) sit at the CELL CENTERS
    of the mesh's node grid — shape ``(len(mesh.x)-1, len(mesh.y)-1)``,
    verified directly against a solved mode, not merely assumed — not at
    ``mesh.x``/``mesh.y`` themselves. A caller sampling the OTHER engine's
    field for comparison (e.g. via ``FDTDSolver.section_mode_fields()``)
    must sample at these midpoints, not the raw node coordinates, or the two
    engines' field grids will be offset by half a cell.

    ``mesh.x`` is beef's device *y* (in-plane transverse); ``mesh.y`` is
    device *z* (vertical) — this is ``EMESolver._build_mesh()``'s own
    convention, unchanged here.
    """
    x = np.asarray(mesh.x, dtype=float)
    y = np.asarray(mesh.y, dtype=float)
    return 0.5 * (x[:-1] + x[1:]), 0.5 * (y[:-1] + y[1:])


def stack_meow_mode_fields(modes) -> FieldSet:
    """``{'Ex','Ey','Ez','Hx','Hy','Hz'}`` stacked to ``[n_modes, ny, nz]``
    complex from a list of meow ``Mode`` objects, REMAPPED into beef's
    device frame — the same shape/axis convention
    ``FDTDSolver.section_mode_fields()`` returns, so either engine's output
    can be passed to ``field_overlap_matrix()`` directly, without this
    module importing meow itself.

    meow's field-component labels are NOT device-frame components — they
    are projections onto meow's OWN local axes, where z is the propagation
    direction (device x) and x, y are transverse (device y, z
    respectively). Confirmed directly, not assumed: for a real solved mode,
    ``sum(mode.Px) == sum(mode.Py) == 0`` and ``sum(mode.Pz) != 0`` — all
    of the mode's real power flows along meow's z, which must therefore be
    the propagation axis. (This is a DIFFERENT axis than the spatial mesh
    convention checked in ``mesh_cell_centers()`` — ``Mesh2D.x``/``.y`` are
    already device y/z; it is specifically the FIELD VECTOR component
    labels that are rotated one axis further, since meow's local frame
    always has propagation along its own z regardless of which device axis
    that maps to.) So the remap is:

        device Ex (propagation)  = meow Ez
        device Ey (in-plane)     = meow Ex
        device Ez (vertical)     = meow Ey
        device Hx/Hy/Hz          = meow Hz/Hx/Hy, same pattern

    Skipping this remap silently produces near-zero overlap for genuinely
    matching modes (caught directly: comparing meep against itself gave a
    clean identity/near-1 diagonal, but meep against un-remapped meow gave
    worst-match diagonal overlap of 3e-4 for a device already known good at
    the n_eff level).
    """
    return dict(
        Ex=np.stack([np.asarray(m.Ez) for m in modes]),
        Ey=np.stack([np.asarray(m.Ex) for m in modes]),
        Ez=np.stack([np.asarray(m.Ey) for m in modes]),
        Hx=np.stack([np.asarray(m.Hz) for m in modes]),
        Hy=np.stack([np.asarray(m.Hx) for m in modes]),
        Hz=np.stack([np.asarray(m.Hy) for m in modes]),
    )


def _transverse_overlap(
    Ta_y: np.ndarray, Ta_z: np.ndarray, Tb_y: np.ndarray, Tb_z: np.ndarray,
    dy: float, dz: float,
) -> complex:
    """kappa_ab = integral (Ta_t . Tb_t*) dA over the transverse (y, z)
    plane — the complex inner product of two transverse field-component
    pairs (``T`` is whichever pair ``field_overlap_matrix()`` was pointed
    at: E for the primary check, H for its independent corroboration), NOT
    the E-x-H reciprocity overlap (see module docstring for why: crossing
    engine A's field with engine B's field of the SAME type needs only ONE
    engine's amplitude convention at a time, so it is immune to meep and
    meow reporting E and H in mutually inconsistent absolute units, which
    the reciprocity form is not — and the same argument applies verbatim to
    H against H).

    Riemann sum over (dy, dz): both engines' field samples are already
    CELL-CENTER values (see ``mesh_cell_centers()``), not node values, so a
    rectangle rule is the matching quadrature — ``np.trapz`` would assume
    node samples and is the wrong rule here.
    """
    integrand = Ta_y * np.conj(Tb_y) + Ta_z * np.conj(Tb_z)
    return complex(np.sum(integrand) * dy * dz)


def field_overlap_matrix(
    fields_a: FieldSet, fields_b: FieldSet, dy: float, dz: float,
    components: Tuple[str, str] = ("Ey", "Ez"),
) -> np.ndarray:
    """``[n_a, n_b]`` matrix of normalized overlap magnitudes eta_ij in
    ``[0, 1]`` between every mode of engine A and every mode of engine B,
    sampled on the SAME ``(y, z)`` grid (uniform spacing ``dy, dz``).

    ``eta_ij = |kappa_ij| / sqrt(kappa_ii * kappa_jj)``, normalized by each
    mode's own transverse power in the chosen field type (always real and
    positive by construction, unlike a mixed E-x-H self term) — so the
    result does not depend on either engine's internal field-amplitude
    normalization convention, only on the chosen component pair's shape.

    ``components`` selects which transverse pair to compare — ``("Ey",
    "Ez")`` (the default) for the primary E-field check, or ``("Hy", "Hz")``
    for the H-only corroboration (`fdtd_vs_EME_mode_validation.md` item 4):
    a field-type-pure overlap is required either way (see module docstring),
    it just need not always be E. Whichever pair is unused (Ex/Hx/Hy/Hz when
    comparing E, or Ex/Ey/Ez when comparing H) may still be present in
    ``fields_a``/``fields_b`` (e.g. because ``section_mode_fields()`` always
    returns all six) but is not read here.
    """
    cy, cz = components
    na = fields_a[cy].shape[0]
    nb = fields_b[cy].shape[0]

    def self_overlap(fields, i):
        return _transverse_overlap(
            fields[cy][i], fields[cz][i], fields[cy][i], fields[cz][i], dy, dz)

    self_a = np.array([abs(self_overlap(fields_a, i)) for i in range(na)])
    self_b = np.array([abs(self_overlap(fields_b, j)) for j in range(nb)])

    eta = np.zeros((na, nb))
    for i in range(na):
        for j in range(nb):
            kappa = _transverse_overlap(
                fields_a[cy][i], fields_a[cz][i], fields_b[cy][j], fields_b[cz][j], dy, dz)
            denom = np.sqrt(self_a[i] * self_b[j])
            eta[i, j] = abs(kappa) / denom if denom > 0 else 0.0
    return eta


def best_correspondence(eta: np.ndarray) -> Tuple[List[int], np.ndarray]:
    """Optimal one-to-one assignment maximizing total overlap (the Hungarian
    algorithm on ``-eta``), replacing nearest-neighbour-in-a-scalar as the
    matching rule — see module docstring for why that generalization
    matters once modes crowd close enough in n_eff to fool it.

    Returns ``(correspondence, margin)``: ``correspondence[i]`` is row
    ``i``'s matched column; ``margin[i]`` is that match's eta minus the
    next-best eta in row ``i`` — a small margin flags a genuinely ambiguous
    pairing (e.g. a true physical degeneracy) rather than confidently
    resolving it, and is worth reporting alongside the match rather than
    discarding.
    """
    row_ind, col_ind = linear_sum_assignment(-eta)
    correspondence = [-1] * eta.shape[0]
    for r, c in zip(row_ind, col_ind):
        correspondence[r] = int(c)

    margin = np.full(eta.shape[0], np.nan)
    for i in range(eta.shape[0]):
        j = correspondence[i]
        if j < 0:
            continue
        row = eta[i].copy()
        best = row[j]
        row[j] = -np.inf
        second = row.max() if eta.shape[1] > 1 else 0.0
        margin[i] = best - second
    return correspondence, margin


def overlap_diagnostics(eta: np.ndarray) -> Dict[str, Any]:
    """Matrix-level statistics for an overlap matrix — the whole matrix, not
    just its diagonal.

    If two bases span the same subspace and correspond mode-for-mode, ``eta``
    is close to a **permutation matrix**. Testing only the matched diagonal
    (the original check) tests a strict subset of that statement and is blind
    to the failure that matters most for an arbitrary geometry: a mode present
    in one engine's basis and *absent* from the other's.

    The row-completeness sum is what closes that gap, and it is free — the
    same matrix, one more reduction. For an orthonormal basis ``{b_j}`` and a
    unit vector ``a_i``, ``sum_j |<a_i, b_j>|^2 == 1`` **iff** ``a_i`` lies in
    ``span{b_j}``. Since ``eta_ij`` is exactly that normalized inner product
    (see ``field_overlap_matrix()``), ``row_completeness[i]`` measures how
    much of engine A's mode *i* engine B's retained basis can represent at
    all. Noticeably below 1 means truncation or a genuine basis mismatch, not
    a mis-ordering — a distinction the diagonal cannot draw. ``col_completeness``
    is the same test in reverse, catching a mode B found that A missed.

    Caveat on interpretation: the identity relies on each engine's own modes
    being mutually orthogonal. Guided eigenmodes of one operator are, to
    discretization error; radiation modes are not, which is one more reason to
    filter them out before calling this (see this module's callers).

    Returns a dict with ``correspondence``/``margin`` (as
    ``best_correspondence()``), plus ``diagonal``, ``row_completeness``,
    ``col_completeness``, ``off_diagonal_max``, and the scalar summaries
    ``worst_diagonal``/``worst_row_completeness``/``worst_col_completeness``.
    Report the arrays, not only the scalars: a failure's *shape* across modes
    is what distinguishes under-resolution (high modes degrade, low ones
    fine) from a basis mismatch (one row collapses) from a mis-ordering (every
    row high, on the wrong column).
    """
    eta = np.asarray(eta, dtype=float)
    n_a, n_b = eta.shape
    correspondence, margin = best_correspondence(eta)

    diagonal = np.full(n_a, np.nan)
    off_diagonal_max = np.full(n_a, np.nan)
    for i in range(n_a):
        j = correspondence[i]
        if j < 0:
            continue
        diagonal[i] = eta[i, j]
        row = eta[i].copy()
        row[j] = -np.inf
        off_diagonal_max[i] = row.max() if n_b > 1 else 0.0

    row_completeness = np.sum(eta ** 2, axis=1)
    col_completeness = np.sum(eta ** 2, axis=0)

    def _worst(arr):
        finite = arr[np.isfinite(arr)]
        return float(finite.min()) if finite.size else float("nan")

    return dict(
        correspondence=correspondence,
        margin=margin,
        diagonal=diagonal,
        row_completeness=row_completeness,
        col_completeness=col_completeness,
        off_diagonal_max=off_diagonal_max,
        worst_diagonal=_worst(diagonal),
        worst_row_completeness=_worst(row_completeness),
        worst_col_completeness=_worst(col_completeness),
        order_preserved=correspondence == list(range(n_a)),
    )


def self_overlap_diagnostics(
    fields: FieldSet, dy: float, dz: float,
    components: Tuple[str, str] = ("Ey", "Ez"),
) -> Dict[str, Any]:
    """Same-engine control: one engine's mode set compared against ITSELF.

    The result must be the identity — unit diagonal, ~zero off-diagonal — by
    construction, because these are eigenmodes of a single operator sampled on
    a single grid. Anything else is a fault in the *measurement*, not a
    finding about the device: too coarse a comparison grid, a bad axis remap,
    a duplicate eigenvector, or non-orthogonal radiation modes in the set.

    Worth running as a standing control alongside every cross-engine
    comparison. This session produced two confident-looking cross-engine
    "correspondence failures" that were both artifacts of the checker; either
    would have been caught immediately by this, and an ad-hoc attempt to run
    it after the fact was what finally localized them. Cheap: it reuses the
    modes already solved for the real comparison, adding no eigensolve.
    """
    eta = field_overlap_matrix(fields, fields, dy=dy, dz=dz, components=components)
    diag = overlap_diagnostics(eta)
    n = eta.shape[0]
    off = eta - np.diag(np.diag(eta))
    max_off = 0.0 if n == 1 else float(np.max(np.abs(off)))

    # Two DIFFERENT failures, and conflating them was a real bug: an early
    # version failed the control on any off-diagonal above 0.1, which turned
    # out to reject bases that are perfectly usable. Measured at the modal
    # multiplexer's interface 0 (16 modes, res 40): meep's own basis shows a
    # max off-diagonal of 0.1257 and meow's 0.1886, from genuinely
    # near-degenerate mode PAIRS -- physically expected on a symmetric
    # multimode cross-section -- while rows 0-8 of the cross-engine
    # comparison matched at diagonal 1.0000 regardless. Treating that as a
    # control failure would disable the check on exactly the devices it
    # exists for.
    #
    #   has_duplicate  off-diagonal near 1 -> two entries ARE one mode. The
    #                  basis is rank-deficient; nothing downstream of it is
    #                  meaningful. Hard fail.
    #   is_orthogonal  off-diagonal small -> textbook-clean basis. Its
    #                  absence WITHOUT a duplicate is worth reporting but not
    #                  disqualifying.
    diag["max_off_diagonal"] = max_off
    diag["has_duplicate"] = bool(max_off > _DUPLICATE_OVERLAP)
    diag["is_orthogonal"] = bool(max_off < _ORTHOGONALITY_TOL)
    diag["diagonal_is_unit"] = bool(np.allclose(np.diag(eta), 1.0, atol=1e-6))
    # Retained for callers that just want "textbook clean", but note it is the
    # STRICT reading: use has_duplicate to gate, not this.
    diag["is_identity"] = bool(diag["diagonal_is_unit"] and diag["is_orthogonal"])
    return diag


def excited_mode_axis(n_modes_out: int, n_modes_in: int) -> int:
    """Which axis of a ``solve_section()``-shaped transfer matrix (``[n_wl,
    n_out, n_in]``) indexes the side ``solve_section()`` actually EXCITES
    from: axis 2 (input) if ``n_in <= n_out``, else axis 1 (output).

    Mirrors ``solve_section()``'s own reciprocity rule verbatim (excite
    whichever boundary has FEWER modes, then use reciprocity for the rest —
    see its docstring): a duplicate eigenvector can only arise among modes
    that were actually excited via ``get_eigenmode_coefficients()``/meow's
    own decomposition, not the side that was merely measured by reciprocity.
    """
    return 2 if n_modes_in <= n_modes_out else 1


def transfer_matrix_duplicate_scan(T: np.ndarray, threshold: float = _DUPLICATE_OVERLAP) -> Dict[str, Any]:
    """Post-hoc duplicate-mode scan of a ``solve_section()``-shaped transfer
    matrix — ``fdtd_vs_EME_mode_validation.md`` item 8's "(a) and (c)" gap.

    Every other mechanism in this module compares bases from
    ``section_mode_fields()``/``section_modes()``: cheap, separate
    eigensolves run purely for diagnosis. The basis ``HybridSolver`` actually
    CONTRACTS comes from a different code path entirely
    (``get_eigenmode_coefficients()`` on DFT-decayed flux for
    ``FDTDSolver``, meow's own mode decomposition for ``EMESolver``) and was
    never checked — a duplicate could exist in the transfer matrix actually
    used even when the side-channel diagnostic reports a clean basis.

    This is the "cheapest first step" the item calls for: `T` is already
    computed by the time ``solve_section()`` returns, so the scan costs no
    extra eigensolve or simulation, only ``O(n_excited^2)`` vector arithmetic
    on an array already in memory.

    ``solve_section()`` excites from whichever side has FEWER modes (see its
    own docstring), so a duplicate excited mode shows up as two nearly
    PARALLEL vectors along that side's axis (``excited_mode_axis()``), each
    vector spanning the wavelength-resolved response measured on the OTHER
    side — two "different" excited modes producing indistinguishable
    responses is the transfer-matrix analogue of
    ``_iter_distinct_bands()``/``EMESolver.compute_modes()``'s duplicate-
    eigenvector guard, but observed on the matrix ``HybridSolver`` actually
    contracts rather than a side-channel eigensolve.

    Similarity is complex cosine similarity (normalized inner product
    magnitude, in ``[0, 1]``) between two excited-side vectors, maxed over
    wavelength — the same normalization ``field_overlap_matrix()`` uses, so
    `threshold` defaults to this module's own ``_DUPLICATE_OVERLAP`` (also
    ``fdtd.py``'s ``_DUP_BAND_OVERLAP_THRESHOLD`` / ``eme.py``'s
    ``_DUP_MODE_OVERLAP_THRESHOLD``): one number for "these are the same
    mode" everywhere it is asked.

    **This is a symptom, not proof — unlike the same-engine field-shape
    controls (``self_overlap_diagnostics()``), a high similarity here is not
    a hard failure.** Two excited modes could in principle produce
    correlated measured responses for a reason other than being the same
    eigenvector (e.g. strong coupling through a shared radiation channel);
    report and investigate `has_duplicate`, don't gate on it automatically
    the way the field-based guards do.

    Args:
        T: complex ndarray ``[n_wavelengths, n_modes_out, n_modes_in]`` —
            ``solve_section()``'s own return shape.
        threshold: similarity above which two excited-side columns are
            flagged as a likely duplicate.

    Returns:
        dict with ``excited_axis`` (1 or 2), ``n_excited``, ``similarity``
        (``[n_excited, n_excited]`` symmetric, wavelength-maxed, unit
        diagonal), ``worst_pair`` (``(i, j)``, or None if ``n_excited < 2``),
        ``max_similarity``, and ``has_duplicate``.
    """
    T = np.asarray(T, dtype=complex)
    if T.ndim != 3:
        raise ValueError(
            f"transfer_matrix_duplicate_scan: expected a 3-D "
            f"[n_wavelengths, n_modes_out, n_modes_in] array, got shape {T.shape}."
        )
    n_out, n_in = T.shape[1], T.shape[2]
    axis = excited_mode_axis(n_out, n_in)
    vectors = np.moveaxis(T, axis, 0)  # [n_excited, n_wl, n_other]
    n_excited = vectors.shape[0]

    similarity = np.eye(n_excited)
    worst_pair = None
    worst = 0.0
    norms = np.sqrt(np.sum(np.abs(vectors) ** 2, axis=2))  # [n_excited, n_wl]
    for i in range(n_excited):
        for j in range(i + 1, n_excited):
            num = np.abs(np.sum(vectors[i] * np.conj(vectors[j]), axis=1))  # [n_wl]
            denom = norms[i] * norms[j]
            with np.errstate(divide="ignore", invalid="ignore"):
                sim_per_wl = np.where(denom > 0, num / denom, 0.0)
            sim = float(np.max(sim_per_wl)) if sim_per_wl.size else 0.0
            similarity[i, j] = similarity[j, i] = sim
            if sim > worst:
                worst, worst_pair = sim, (i, j)

    return dict(
        excited_axis=axis,
        n_excited=n_excited,
        similarity=similarity,
        worst_pair=worst_pair,
        max_similarity=worst,
        has_duplicate=bool(worst > threshold),
    )


def warn_on_transfer_duplicate(
    T: np.ndarray, logger: logging.Logger, x0: float, x1: float,
    threshold: float = _DUPLICATE_OVERLAP,
) -> Dict[str, Any]:
    """Run ``transfer_matrix_duplicate_scan()`` on a just-computed
    ``solve_section()`` transfer matrix and log a WARNING (never raise) if it
    flags a likely duplicate excited mode.

    Shared by ``FDTDSolver.solve_section()`` and ``EMESolver.solve_section()``
    (``fdtd_vs_EME_mode_validation.md`` item 8) so the message and threshold
    cannot drift between the two engines' call sites. See
    ``transfer_matrix_duplicate_scan()``'s docstring for what the scan
    measures and why this is reported rather than gated on.

    Returns the scan dict (see ``transfer_matrix_duplicate_scan()``), so a
    caller — e.g. a test — can inspect it without re-running the scan.
    """
    scan = transfer_matrix_duplicate_scan(T, threshold=threshold)
    if scan["has_duplicate"]:
        i, j = scan["worst_pair"]
        logger.warning(
            "solve_section: possible duplicate excited mode in the transfer "
            "matrix at x=[%.4g, %.4g] -- excited-side modes %d and %d produced "
            "near-parallel measured responses (similarity %.4f). This is a "
            "symptom, not a confirmed duplicate (see "
            "transfer_matrix_duplicate_scan()'s docstring); investigate before "
            "trusting this section's basis. See fdtd_vs_EME_mode_validation.md "
            "item 8.",
            x0, x1, i, j, scan["max_similarity"],
        )
    return scan
