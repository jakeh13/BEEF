"""Ring resonator validation checks 1 and 2, actually RUN.

ring_resonator_validation.md records these two as unrun. Until 2026-08-03 they
were also *unrunnable*: the ring's EME half builds its cross-sections from raw
polygons and calls `mw.compute_modes_meep` directly (`solve_neff()`,
`solve_body_supermodes()` in examples/ring_resonator/python/
ring_resonator_hybrid_cascaded.py), so there was no DeviceGeometry both engines
could consume, and "engine agreement" had nothing to agree about.
`UniformCrossSection` is that shared object. This script uses it.

  CHECK 1 -- engine agreement.  Does EMESolver's mode-solve path, driven through
  a DeviceGeometry, return the same n_eff as the ring's hand-rolled call?
  Both the bare strip (solve_neff) and the coupled pair (solve_body_supermodes).

  CHECK 2 -- mode identity.  Are they the same MODES, not merely the same
  numbers? Two different structures can share an n_eff. This takes the field
  overlap matrix between the two engines' mode sets and asks for a clean
  one-to-one correspondence, against a same-engine self-overlap control.

WHAT MAKES THIS A REAL COMPARISON, AND WHAT DOES NOT
----------------------------------------------------
The ring centres its coupled-pair oxide on BODY_Y_CENTER = Y_RING/2; every
DeviceGeometry's `section_gds()` centres its oxide on y=0 (`_oxide_box`). So an
exact reproduction is impossible without a y-shift, and pretending otherwise
would compare two different structures.

The shift is applied HONESTLY: the strips are handed to UniformCrossSection at
(-Y_RING/2, +Y_RING/2) instead of (0, Y_RING), and the mesh handed to
`section_modes()` is the ring's own mesh translated by the same -BODY_Y_CENTER.
Translating structure and mesh together is an exact symmetry of the discretised
eigenproblem -- same dx, dy, same pixel count, same dielectric pattern relative
to the grid -- so n_eff must agree to solver tolerance. Any disagreement is a
genuine engine difference and not a gridding artefact. TEST 0 below measures
that claim rather than asserting it: it re-runs the ring's own recipe shifted,
and the shifted-vs-unshifted delta is the floor every other number is read
against.

Mode solve only -- no FDTD, no time stepping. Costs a handful of MPB solves.

Env:
    RING_EA_RES   mode-solve resolution (default 100.0, the ring's shipped
                  MODE_RESOLUTION; the whole point is MATCHED resolution)
    RING_EA_WL    comma-separated wavelengths (default 1.50,1.55,1.60)
    RING_EA_OUT   output dir override
"""
import os
import sys

import numpy as np
import meow as mw
from meow import compute_modes_meep

from beef.materials import SOI_220
from beef import extrude_polygons
from beef.geometry import UniformCrossSection
from beef.solvers import EMESolver
from beef.mode_identity import (
    mesh_cell_centers,
    stack_meow_mode_fields,
    field_overlap_matrix,
    overlap_diagnostics,
    self_overlap_diagnostics,
)

# --- the ring's own constants. Mirrored, not imported: importing the ring
#     script executes its whole pipeline. Kept in sync by TEST 0, which
#     re-derives the reference numbers from this same recipe.
#     == ring_resonator_hybrid_cascaded.py:145-146,147,177,253-265
MATERIAL_STACK = SOI_220
W = 0.45                      # WAVEGUIDE_WIDTH
GAP = 0.25
Y_RING = W + GAP              # 0.70
MESH_WIDTH, MESH_HEIGHT = 4.0, 2.0
_MODE_MARGIN = MESH_WIDTH / 2 - W / 2
BODY_Y_CENTER = Y_RING / 2
BODY_MESH_WIDTH = Y_RING + W + 2 * _MODE_MARGIN
N_MODES = 2

RES = float(os.environ.get("RING_EA_RES", 100.0))
WAVELENGTHS = [float(v) for v in
               os.environ.get("RING_EA_WL", "1.50,1.55,1.60").split(",")]
OUT = os.environ.get("RING_EA_OUT") or os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "..",
    "example_output", "verification", "ring_eme_engine_agreement"))

# Tolerances. n_eff at res 100 is converged to ~1e-5 (gap 3 sweep, job
# 11634859), so a translation-invariant reproduction should sit near machine
# noise; 1e-6 is loose enough for MPB's iterative tolerance and tight enough
# that a real structural difference (a pixel shift, a different oxide span)
# cannot hide under it.
NEFF_TOL = 1e-6
# Field overlap: >= 0.99 is the same mode. The duplicate-mode guard elsewhere in
# beef uses 0.9 as "suspiciously similar"; identity should be far above that.
OVERLAP_TOL = 0.99


def rect_poly(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def check_square_mesh(res):
    """`compute_modes_meep` requires equal dx/dy (meep's resolution is
    isotropic). Both cross-sections' grids are derived from `res`, so most
    values of `res` are simply ILLEGAL here — e.g. res 14 gives the coupled
    pair dx=0.071212 against dy=0.071429 and meow raises.

    Reported up-front rather than as a raw meow traceback three minutes into
    the run, because the distinction matters: this is the ring's own mesh
    arithmetic constraining which resolutions exist, NOT a disagreement
    between the two engines. Silently rounding to a legal grid would be worse
    -- it would change the structure under test.
    """
    dz = MESH_HEIGHT / round(res * MESH_HEIGHT)
    grids = {
        "bare strip": MESH_WIDTH / (int(res * MESH_WIDTH)),
        "coupled pair": BODY_MESH_WIDTH / (round(res * BODY_MESH_WIDTH)),
    }
    bad = {k: v for k, v in grids.items() if abs(v - dz) > 1e-12}
    if bad:
        detail = "; ".join(f"{k}: dy={v:.9f} vs dz={dz:.9f}" for k, v in bad.items())
        raise SystemExit(
            f"[ea] res={res:g} is not a legal mode-solve resolution for this "
            f"device: {detail}.\n"
            f"[ea] compute_modes_meep needs dy == dz. The ring's shipped "
            f"MODE_RESOLUTION=100 is legal for both cross-sections "
            f"(dy=dz=0.01); so are 20, 30, 40, 50, 60, 80. Notably 14 and 25 "
            f"are NOT. Pick a legal one via RING_EA_RES.")


# ---------------------------------------------------------------------------
# The ring's own recipe, verbatim, with an optional rigid y-shift.
# ---------------------------------------------------------------------------
def ring_recipe_modes(wl, res, coupled, y_shift=0.0, length=1.0):
    """`solve_neff()` (coupled=False) or `solve_body_supermodes()`
    (coupled=True), with every polygon and both mesh axes translated by
    `y_shift`. Returns (modes, mesh) so the caller can take fields, not just
    n_eff -- check 2 needs the fields.

    The mesh construction below mirrors the ring exactly, including the
    int()-vs-round() split: solve_neff() uses int(res*MESH_WIDTH) on a literal
    4.0 while solve_body_supermodes() uses round(res*BODY_MESH_WIDTH) on a
    derived float. Reproducing that asymmetry is deliberate -- "harmonising" it
    here would silently change the grid and make the comparison meaningless.
    """
    if coupled:
        wg = [rect_poly(0.0, y_shift - W / 2, length, y_shift + W / 2),
              rect_poly(0.0, y_shift + Y_RING - W / 2, length, y_shift + Y_RING + W / 2)]
        ox = [rect_poly(0.0, y_shift + BODY_Y_CENTER - BODY_MESH_WIDTH / 2,
                        length, y_shift + BODY_Y_CENTER + BODY_MESH_WIDTH / 2)]
        x_lo = y_shift + BODY_Y_CENTER - BODY_MESH_WIDTH / 2
        x_hi = y_shift + BODY_Y_CENTER + BODY_MESH_WIDTH / 2
        nx = round(res * BODY_MESH_WIDTH) + 1
        n_modes = N_MODES
    else:
        wg = [rect_poly(0.0, y_shift - W / 2, length, y_shift + W / 2)]
        ox = [rect_poly(0.0, y_shift - MESH_WIDTH / 2, length, y_shift + MESH_WIDTH / 2)]
        x_lo, x_hi = y_shift - MESH_WIDTH / 2, y_shift + MESH_WIDTH / 2
        nx = int(res * MESH_WIDTH) + 1
        n_modes = 1

    structs = extrude_polygons({(2, 0): wg, (1, 0): ox},
                               MATERIAL_STACK.extrusion_rules)
    mesh = mw.Mesh2D(
        x=np.linspace(x_lo, x_hi, nx),
        y=np.linspace(-MESH_HEIGHT / 2, MESH_HEIGHT / 2, round(res * MESH_HEIGHT) + 1),
    )
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=float(wl), T=25.0)
    cs = mw.CrossSection.from_cell(cell=cells[0], env=env)
    modes = compute_modes_meep(cs, num_modes=n_modes)
    # solve_body_supermodes() sorts by -Re(neff) to pin the even/odd band index;
    # do the same on both sides so the comparison is index-aligned.
    modes = sorted(modes, key=lambda m: -np.real(m.neff))
    return modes, mesh


# ---------------------------------------------------------------------------
# The same cross-section through EMESolver / UniformCrossSection.
# ---------------------------------------------------------------------------
def beef_modes(wl, res, coupled, length=1.0):
    """EMESolver.section_modes() on a UniformCrossSection carrying the ring's
    cross-section, on the ring's mesh translated to y=0.

    The strips are placed straddling y=0 (-Y_RING/2, +Y_RING/2) rather than at
    (0, Y_RING) because section_gds() centres its oxide on y=0. That is a rigid
    translation by -BODY_Y_CENTER of the ring's structure, and the mesh below is
    translated identically -- see the module docstring, and TEST 0 for the
    measured cost of that translation.
    """
    if coupled:
        strips = [(-Y_RING / 2, W), (+Y_RING / 2, W)]
        oxide_width = BODY_MESH_WIDTH
        nx = round(res * BODY_MESH_WIDTH) + 1
        n_modes = N_MODES
    else:
        strips = [(0.0, W)]
        oxide_width = MESH_WIDTH
        nx = int(res * MESH_WIDTH) + 1
        n_modes = 1

    geom = UniformCrossSection(strips, length=length)
    solver = EMESolver(
        num_modes=[n_modes], num_cells=[1], lengths=[length],
        oxide_width=oxide_width, oxide_height=MESH_HEIGHT,
        res=res, material_stack=MATERIAL_STACK, temperature=25.0,
        # The ring's solve_*() paths pass MATERIAL_STACK.extrusion_rules (the
        # dispersive SampledMaterial form) and build one Environment per
        # wavelength. EMESolver's default (dispersive=None -> True) is the same
        # thing; naming it explicitly so the comparison cannot drift if a
        # default moves.
        dispersive=True,
        # The ring calls compute_modes_meep directly, with no duplicate guard.
        # Matching that: an over-request-and-filter here would solve a DIFFERENT
        # number of modes and is exactly the kind of silent divergence this
        # check exists to catch.
        guard_duplicate_modes=False,
    )
    # Explicit mesh: EMESolver._build_mesh() uses int() on both axes and
    # oxide_height for y, while the ring uses round() on the coupled x-axis.
    # Passing the ring's own grid removes gridding as a variable, which is the
    # only way the n_eff delta below means "engine", not "mesh".
    mesh = mw.Mesh2D(
        x=np.linspace(-oxide_width / 2, oxide_width / 2, nx),
        y=np.linspace(-MESH_HEIGHT / 2, MESH_HEIGHT / 2, round(res * MESH_HEIGHT) + 1),
    )
    modes, mesh = solver.section_modes(
        geom, wavelength=float(wl), x0=-length / 2, x1=length / 2,
        n_modes=n_modes, mesh=mesh,
    )
    modes = sorted(modes, key=lambda m: -np.real(m.neff))
    return modes, mesh


# ---------------------------------------------------------------------------
def main():
    check_square_mesh(RES)
    print(f"[ea] res={RES}  wavelengths={WAVELENGTHS}", flush=True)
    print(f"[ea] neff tol={NEFF_TOL:.1e}  overlap tol={OVERLAP_TOL}", flush=True)
    os.makedirs(OUT, exist_ok=True)

    rows = []
    failures = []

    for coupled, label in ((False, "bare strip"), (True, "coupled pair")):
        n_modes = N_MODES if coupled else 1
        print(f"\n[ea] ================ {label} ================", flush=True)

        for wl in WAVELENGTHS:
            # --- TEST 0: translation floor -------------------------------
            # How much does the rigid y-shift alone move n_eff, inside ONE
            # engine? This is the noise floor for TEST 1; a TEST 1 delta at or
            # below it says nothing about the engines.
            ring_modes, ring_mesh = ring_recipe_modes(wl, RES, coupled, y_shift=0.0)
            shift_modes, _ = ring_recipe_modes(wl, RES, coupled,
                                               y_shift=-BODY_Y_CENTER if coupled else 0.0)
            n_ring = np.array([np.real(m.neff) for m in ring_modes])
            n_shift = np.array([np.real(m.neff) for m in shift_modes])
            floor = float(np.max(np.abs(n_ring - n_shift)))

            # --- TEST 1: engine agreement (check 1) ----------------------
            bf_modes, bf_mesh = beef_modes(wl, RES, coupled)
            n_beef = np.array([np.real(m.neff) for m in bf_modes])
            d_neff = np.abs(n_ring - n_beef)
            worst = float(np.max(d_neff))

            # --- TEST 2: mode identity (check 2) -------------------------
            # Compare the SHIFTED ring modes against beef's, since those two
            # sit on the same grid; comparing to the unshifted set would be an
            # overlap between fields sampled at different y and would fail for
            # a reason that is not about identity.
            y_a, z_a = mesh_cell_centers(ring_mesh)
            y_b, z_b = mesh_cell_centers(bf_mesh)
            # field_overlap_matrix() integrates on ONE shared grid, so the two
            # sample grids must have the same shape and spacing. They do by
            # construction (both built from the same nx/round(res*MESH_HEIGHT)),
            # but assert it: a silent shape mismatch here would broadcast into a
            # meaningless overlap rather than an error.
            assert (y_a.shape, z_a.shape) == (y_b.shape, z_b.shape), (
                f"grid shape mismatch {y_a.shape},{z_a.shape} vs "
                f"{y_b.shape},{z_b.shape}")
            dy = float(y_a[1] - y_a[0])
            dz = float(z_a[1] - z_a[0])
            assert abs(dy - (y_b[1] - y_b[0])) < 1e-12, "dy mismatch"
            assert abs(dz - (z_b[1] - z_b[0])) < 1e-12, "dz mismatch"

            fields_a = stack_meow_mode_fields(shift_modes)
            fields_b = stack_meow_mode_fields(bf_modes)
            eta = field_overlap_matrix(fields_a, fields_b, dy=dy, dz=dz)
            diag = overlap_diagnostics(eta)
            ctrl = self_overlap_diagnostics(fields_a, dy=dy, dz=dz)
            perm = diag["correspondence"]
            worst_pair = float(diag["worst_diagonal"])
            is_identity = list(perm) == list(range(n_modes))

            print(f"[ea] wl={wl:.3f}", flush=True)
            for i in range(n_modes):
                print(f"[ea]   mode {i}: ring n_eff={n_ring[i]:.9f}  "
                      f"beef n_eff={n_beef[i]:.9f}  |d|={d_neff[i]:.3e}", flush=True)
            print(f"[ea]   TEST 0 translation floor          : {floor:.3e}", flush=True)
            print(f"[ea]   TEST 1 worst |d n_eff|            : {worst:.3e}"
                  f"   ({'PASS' if worst <= NEFF_TOL else 'FAIL'})", flush=True)
            print(f"[ea]   TEST 2 worst matched overlap      : {worst_pair:.6f}"
                  f"   ({'PASS' if worst_pair >= OVERLAP_TOL else 'FAIL'})", flush=True)
            print(f"[ea]          correspondence             : {list(perm)}"
                  f"   ({'identity' if is_identity else 'PERMUTED'})", flush=True)
            print(f"[ea]          worst off-diagonal / margin: "
                  f"{np.nanmax(diag['off_diagonal_max']):.6f} / "
                  f"{np.nanmin(diag['margin']):.6f}", flush=True)
            # row/col completeness catches the failure the diagonal is blind
            # to: a mode one engine found and the other's basis cannot
            # represent at all.
            print(f"[ea]          worst row / col completeness: "
                  f"{diag['worst_row_completeness']:.6f} / "
                  f"{diag['worst_col_completeness']:.6f}", flush=True)
            print(f"[ea]          same-engine control        : "
                  f"diag unit={ctrl['diagonal_is_unit']}, "
                  f"max off-diag {ctrl['max_off_diagonal']:.6f}, "
                  f"duplicate={ctrl['has_duplicate']}", flush=True)
            if ctrl["has_duplicate"]:
                failures.append(
                    f"{label} wl={wl:.3f}: CONTROL FAILED -- the ring's own "
                    f"mode set contains a duplicate (max off-diagonal "
                    f"{ctrl['max_off_diagonal']:.4f}). The basis is "
                    f"rank-deficient, so the cross-engine numbers above mean "
                    f"nothing.")

            if worst > NEFF_TOL:
                failures.append(
                    f"{label} wl={wl:.3f}: n_eff differs by {worst:.3e} > "
                    f"{NEFF_TOL:.1e} (translation floor {floor:.3e}) -- "
                    f"EMESolver and the ring's hand-rolled solve are NOT "
                    f"modelling the same cross-section")
            if worst_pair < OVERLAP_TOL:
                failures.append(
                    f"{label} wl={wl:.3f}: worst matched field overlap "
                    f"{worst_pair:.4f} < {OVERLAP_TOL} -- same n_eff, different "
                    f"mode; the numbers agreeing is a coincidence")
            if not is_identity:
                failures.append(
                    f"{label} wl={wl:.3f}: best correspondence is {list(perm)}, "
                    f"not identity -- the two engines' band ORDER disagrees, so "
                    f"the ring's (even, odd) index convention does not carry "
                    f"across. This is exactly the failure that would corrupt "
                    f"P_body's diagonal.")

            rows.append(dict(
                section=label, wavelength=wl, resolution=RES,
                neff_ring=n_ring, neff_beef=n_beef,
                translation_floor=floor, worst_dneff=worst,
                eta=eta, perm=np.asarray(perm), worst_overlap=worst_pair,
            ))

    # --- verdict ---------------------------------------------------------
    print("\n[ea] =========================== VERDICT ===========================",
          flush=True)
    if failures:
        print(f"[ea] CHECKS 1/2 FAIL ({len(failures)}):", flush=True)
        for f in failures:
            print(f"[ea]   - {f}", flush=True)
    else:
        print("[ea] CHECKS 1 and 2 PASS.", flush=True)
        print("[ea]   Check 1: EMESolver's mode solve, driven through a "
              "DeviceGeometry, reproduces", flush=True)
        print("[ea]            the ring's hand-rolled solve_neff() and "
              "solve_body_supermodes() n_eff", flush=True)
        print(f"[ea]            to better than {NEFF_TOL:.0e} at matched "
              f"resolution {RES:g}.", flush=True)
        print("[ea]   Check 2: the mode sets are the SAME modes -- one-to-one "
              "field correspondence", flush=True)
        print(f"[ea]            with every matched overlap >= {OVERLAP_TOL}, "
              "identity permutation, so", flush=True)
        print("[ea]            the ring's (even, odd) band index carries "
              "across engines.", flush=True)
        print("[ea]   SCOPE: this is the EME/mode-solve half only. It does NOT "
              "compare against", flush=True)
        print("[ea]          FDTD, and it does not test propagation -- "
              "solve_section()'s S-matrix,", flush=True)
        print("[ea]          reflection and reverse transmission are out of "
              "scope here.", flush=True)

    np.savez(os.path.join(OUT, "engine_agreement.npz"),
             resolution=RES, wavelengths=np.asarray(WAVELENGTHS),
             n_failures=len(failures),
             **{f"{r['section'].replace(' ', '_')}_{i}_{k}": v
                for i, r in enumerate(rows) for k, v in r.items()
                if k not in ("section",)})
    print(f"\n[ea] saved -> {os.path.join(OUT, 'engine_agreement.npz')}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
