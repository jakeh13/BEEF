"""Field-overlap mode correspondence: a device- and resolution-agnostic
alternative to nearest-neighbour n_eff matching.

``HybridSolver`` contracts index *i* of an FDTD interface basis against
index *i* of an EME body basis, assuming they name the same mode. The
existing check for that (``tests/test_mode_identity.py``'s
``assert_modes_correspond()``) compares n_eff and picks each mode's nearest
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

from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment

FieldSet = Dict[str, np.ndarray]


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


def _transverse_e_overlap(
    Ea_y: np.ndarray, Ea_z: np.ndarray, Eb_y: np.ndarray, Eb_z: np.ndarray,
    dy: float, dz: float,
) -> complex:
    """kappa_ab = integral (Ea_t . Eb_t*) dA over the transverse (y, z)
    plane — the complex inner product of the two transverse E-field
    vectors, NOT the E-x-H reciprocity overlap (see module docstring for
    why: it needs H from only ONE engine at a time, so it is immune to
    meep and meow reporting E and H in mutually inconsistent absolute
    units, which the reciprocity form is not).

    Riemann sum over (dy, dz): both engines' field samples are already
    CELL-CENTER values (see ``mesh_cell_centers()``), not node values, so a
    rectangle rule is the matching quadrature — ``np.trapz`` would assume
    node samples and is the wrong rule here.
    """
    integrand = Ea_y * np.conj(Eb_y) + Ea_z * np.conj(Eb_z)
    return complex(np.sum(integrand) * dy * dz)


def field_overlap_matrix(fields_a: FieldSet, fields_b: FieldSet, dy: float, dz: float) -> np.ndarray:
    """``[n_a, n_b]`` matrix of normalized overlap magnitudes eta_ij in
    ``[0, 1]`` between every mode of engine A and every mode of engine B,
    sampled on the SAME ``(y, z)`` grid (uniform spacing ``dy, dz``).

    ``eta_ij = |kappa_ij| / sqrt(kappa_ii * kappa_jj)``, normalized by each
    mode's own transverse E-field power (always real and positive by
    construction, unlike a mixed E-x-H self term) — so the result does not
    depend on either engine's internal field-amplitude normalization
    convention, only on Ey/Ez shape. Only Ey, Ez are used; Ex (the
    propagation-direction component), Hx, Hy, Hz may still be present in
    ``fields_a``/``fields_b`` (e.g. because ``section_mode_fields()``
    always returns all six) but are not read here.
    """
    na = fields_a["Ey"].shape[0]
    nb = fields_b["Ey"].shape[0]

    def self_overlap(fields, i):
        return _transverse_e_overlap(
            fields["Ey"][i], fields["Ez"][i], fields["Ey"][i], fields["Ez"][i], dy, dz)

    self_a = np.array([abs(self_overlap(fields_a, i)) for i in range(na)])
    self_b = np.array([abs(self_overlap(fields_b, j)) for j in range(nb)])

    eta = np.zeros((na, nb))
    for i in range(na):
        for j in range(nb):
            kappa = _transverse_e_overlap(
                fields_a["Ey"][i], fields_a["Ez"][i], fields_b["Ey"][j], fields_b["Ez"][j], dy, dz)
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
