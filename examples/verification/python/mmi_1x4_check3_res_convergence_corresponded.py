"""Check 3 (resolution convergence), done properly: track modes by FIELD
SHAPE across resolutions instead of by index.

validation_markdown/mmi_1x4_sin300_validation.md step 8 check 3.

WHY THE FIRST ATTEMPT WAS INCONCLUSIVE, AND WHAT THIS CHANGES.
`mmi_1x4_sin_body_res_convergence.py` (job 11675082, 2026-08-04) swept
res=[15..40] on the body cross-section and reported that n_eff does NOT
converge monotonically anywhere in that range -- successive |dn_eff|
bouncing around 1-2e-3 instead of shrinking. That result left res=40 (and
now res=30) unjustified, which undermines every number this device has
produced.

But that script tracks modes **by index**: line ~110, `neffs[idx]` for
idx in TRACK_INDICES. It assumes mode 10 at res=15 is the same physical
mode as mode 10 at res=40. On this cross-section -- 23 guided modes packed
into a narrow index range -- that assumption is exactly what is expected to
fail, and the validation doc named it as the leading (unverified)
hypothesis: "mode-index correspondence swapping between near-degenerate
higher-order modes ... rather than a genuine solver-convergence failure".

It is the SAME failure class check 2 already hit and resolved on this very
device: nearest-neighbour-in-n_eff matching reported ORDER PRESERVED: False
(job 11697897), and the field-overlap check (job 11706990) then found
ORDER PRESERVED: True with worst diagonal overlap 0.9979 -- a scalar-match
artifact, not a real mismatch.

So this script re-runs the sweep with mode identity established by FIELD
OVERLAP (beef.mode_identity), and reports BOTH:
  * the naive index-tracked |dn_eff| (reproducing the old result), and
  * the correspondence-tracked |dn_eff| (the corrected one),
side by side, so the hypothesis is tested rather than assumed. If the
corrected series converges monotonically while the naive one does not, the
hypothesis is confirmed and res=30 can finally be justified (or not) on
evidence.

METHOD -- the one non-obvious mechanic. Modes solved at different
resolutions live on different meshes and cannot be overlap-compared
directly. Both engines are handled without ever comparing across grids:
  * FDTD/meep: `FDTDSolver.section_mode_fields()` takes an EXPLICIT output
    grid (y_grid, z_grid) independent of the solver's own `res`, so each
    resolution is SOLVED at its own res but SAMPLED onto one common grid.
    Its docstring states this works across resolutions, not just engines.
  * EME/meow: fields come out on the solve mesh's cell centres, so each
    resolution's fields are interpolated onto the same common grid
    (scipy RegularGridInterpolator, linear). Interpolation is applied ONLY
    to the comparison, never to the n_eff being converged -- n_eff is
    always the solver's own value at its own resolution.

The comparison grid is deliberately COARSER than any swept resolution
(CMP_RES), following audit_interface_modes.py's own reasoning: an
aggregate overlap integral needs a handful of points per field lobe, not a
precisely resolved field, and the cost is n_modes * ny * nz * 6 per solve.

SCOPE. Solves N_SOLVE modes so correspondence has room to work, but reports
convergence on the first N_REPORT -- production keeps interface_modes=12,
and the validation doc already records that modes 19/22 never converge in
this range and are irrelevant for that reason. Reporting on modes the
device discards would repeat the first attempt's mistake in a new form.

Run (smoke first -- see the batch file):
    MMI1X4_C3_RES_LIST=20,30 MMI1X4_C3_NSOLVE=8 MMI1X4_C3_NREPORT=4 \
        python .../mmi_1x4_check3_res_convergence_corresponded.py
    python .../mmi_1x4_check3_res_convergence_corresponded.py
"""
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import meow as mw
from scipy.interpolate import RegularGridInterpolator

from beef import MMI, EMESolver, FDTDSolver, SIN_300, extrude_polygons
from beef.mode_identity import (
    best_correspondence, field_overlap_matrix, mesh_cell_centers,
    stack_meow_mode_fields,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_1x4", "check3_res_convergence"))
os.makedirs(output_dir, exist_ok=True)


def _env_list(name, default):
    raw = os.environ.get(name)
    return [float(v) for v in raw.split(",")] if raw else default


def _env_i(name, default):
    return int(os.environ.get(name, default))


# --- pinned to mmi_1x4_hybrid.py's production config ------------------------
WL = 1.55
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 1.1, 2.2, 5.0
WAVEGUIDE_OFFSET, MM_WIDTH = 3.1, 12.0
MM_LENGTH = 39.90
SIM_WIDTH, SIM_HEIGHT = MM_WIDTH + 4.0, 4.0   # 16.0 x 4.0, the 2026-08-19 values
DPML = 0.5

# --- sweep config -----------------------------------------------------------
RES_LIST = _env_list("MMI1X4_C3_RES_LIST", [15.0, 20.0, 25.0, 30.0, 35.0, 40.0])
N_SOLVE = _env_i("MMI1X4_C3_NSOLVE", 16)    # room for correspondence to work
N_REPORT = _env_i("MMI1X4_C3_NREPORT", 12)  # = production interface_modes
CMP_RES = float(os.environ.get("MMI1X4_C3_CMP_RES", 20.0))

mmi = MMI(
    n_inputs=1, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30, use_neighbor_aperture=True,
)
# Body span, clear of both port-taper bands.
#
# BODY_INSET is 1.5 um, not the 0.25 um body_stub the hybrid cascade uses,
# and the difference is load-bearing: FDTDSolver.section_mode_fields() pads
# the requested span by `pad`, which DEFAULTS TO self.dpml (0.5 um). A 0.25 um
# inset therefore asks for a cut at x0-0.5 = -20.20, which lands inside the
# port-taper band [-24.95, -19.95] and is (correctly) refused by
# MMI._check_cut. Caught by the smoke pass, job 12274542.
#
# Insetting rather than passing pad=0.0 is deliberate: the body is a uniform
# cross-section, so WHERE inside it the modes are sampled cannot change the
# answer, whereas pad=0.0 would move PML up against the monitor plane and
# change the simulation. 1.5 um leaves 1.0 um of margin past the 0.5 um pad.
# (The EME side is unaffected -- EMESolver.section_modes() defaults its stubs
# to 0.0 and so never pads -- which is why only FDTD failed in the smoke.)
half = MM_LENGTH / 2
BODY_INSET = float(os.environ.get("MMI1X4_C3_BODY_INSET", 1.5))
X0, X1 = -half + BODY_INSET, half - BODY_INSET

# --- one common comparison grid, independent of every swept resolution ------
_cmp_mesh = mw.Mesh2D(
    x=np.linspace(-SIM_WIDTH / 2, SIM_WIDTH / 2, int(CMP_RES * SIM_WIDTH) + 1),
    y=np.linspace(-SIM_HEIGHT / 2, SIM_HEIGHT / 2, int(CMP_RES * SIM_HEIGHT) + 1),
)
Y_GRID, Z_GRID = mesh_cell_centers(_cmp_mesh)

print(f"CHECK 3 -- resolution convergence with FIELD-OVERLAP mode tracking")
print(f"  body span [{X0:.3f}, {X1:.3f}] um, mm_width={MM_WIDTH}, wl={WL} um")
print(f"  res sweep: {RES_LIST}")
print(f"  solving {N_SOLVE} modes, reporting convergence on first {N_REPORT}")
print(f"  comparison grid: cmp_res={CMP_RES} -> {len(Y_GRID)} x {len(Z_GRID)} points")
print()


def _interp_to_cmp(fields, y_src, z_src):
    """Interpolate [n_modes, ny, nz] complex field arrays onto (Y_GRID, Z_GRID).

    Used for the meow side only -- its fields land on the solve mesh's own
    cell centres. Real and imaginary parts are interpolated separately
    (RegularGridInterpolator is real-valued).
    """
    pts = np.stack(np.meshgrid(Y_GRID, Z_GRID, indexing="ij"), axis=-1)
    out = {}
    for key, arr in fields.items():
        stacked = np.empty((arr.shape[0], len(Y_GRID), len(Z_GRID)), dtype=complex)
        for m in range(arr.shape[0]):
            for part, fn in (("re", np.real), ("im", np.imag)):
                interp = RegularGridInterpolator(
                    (y_src, z_src), fn(arr[m]), bounds_error=False, fill_value=0.0)
                vals = interp(pts)
                if part == "re":
                    stacked[m] = vals
                else:
                    stacked[m] = stacked[m] + 1j * vals
        out[key] = stacked
    return out


def solve_fdtd(res):
    """meep/MPB modes at `res`, sampled on the common grid (no interpolation)."""
    solver = FDTDSolver(
        res=res, sim_width=SIM_WIDTH, sim_height=SIM_HEIGHT, dpml=DPML,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SIN_300, decay_tol=1e-11,
    )
    r = solver.section_mode_fields(
        mmi, WL, X0, X1, N_SOLVE, Y_GRID, Z_GRID,
        monitor_center=0.0, monitor_size=SIM_WIDTH,
    )
    return np.asarray(r["neff"], dtype=float), {k: r[k] for k in
                                                ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")}


def solve_eme(res):
    """meow/tidy3d modes at `res`, then interpolated onto the common grid."""
    eme = EMESolver(
        num_modes=[N_SOLVE], num_cells=[1], lengths=[1.0],
        oxide_width=SIM_WIDTH, oxide_height=SIM_HEIGHT,
        res=res, material_stack=SIN_300, mode_solver="tidy3d",
    )
    modes, mesh = eme.section_modes(mmi, WL, X0, X1, N_SOLVE)
    neff = np.array([float(np.real(m.neff)) for m in modes])
    y_src, z_src = mesh_cell_centers(mesh)
    return neff, _interp_to_cmp(stack_meow_mode_fields(modes), y_src, z_src)


def run_engine(name, solve_fn):
    print(f"===== {name} =====", flush=True)
    neffs, fields = {}, {}
    for res in RES_LIST:
        t0 = time.time()
        neffs[res], fields[res] = solve_fn(res)
        print(f"  res={res:5.1f}  solved {N_SOLVE} modes in {time.time()-t0:7.1f}s"
              f"   n_eff[0]={neffs[res][0]:.6f}", flush=True)

    ref = RES_LIST[-1]  # finest resolution is the correspondence reference
    perms, worst = {}, {}
    for res in RES_LIST:
        if res == ref:
            perms[res] = list(range(N_SOLVE))
            worst[res] = 1.0
            continue
        # ROWS = reference, COLS = this resolution, deliberately:
        # best_correspondence returns correspondence[i] = row i's matched
        # COLUMN, so with the reference on rows, perm[j] is directly "which
        # mode at THIS resolution is reference-mode j" -- which is the lookup
        # neffs[res][perm[j]] below needs. Building eta the other way round
        # returns a reference index where a resolution index is wanted, and
        # the two are only interchangeable when the permutation is identity,
        # i.e. exactly when this check is not needed.
        eta = field_overlap_matrix(fields[ref], fields[res], dy=1.0, dz=1.0)
        perm, _ = best_correspondence(eta)
        perms[res] = perm
        # how confident is the match on the modes we actually report?
        worst[res] = min(abs(eta[j, perm[j]]) for j in range(N_REPORT))

    print(f"\n  correspondence vs res={ref} (identity = no reordering):")
    for res in RES_LIST:
        p = perms[res][:N_REPORT]
        tag = "identity" if p == list(range(N_REPORT)) else f"REORDERED {p}"
        print(f"    res={res:5.1f}  worst|overlap| on reported modes = "
              f"{worst[res]:.4f}   {tag}")

    # naive (index-tracked) vs corresponded |dn_eff| between consecutive res
    print(f"\n  successive |dn_eff|, NAIVE index tracking (the old method):")
    naive = _series(neffs, lambda res, j: neffs[res][j])
    _print_series(naive)

    print(f"\n  successive |dn_eff|, FIELD-OVERLAP corresponded (corrected):")
    corr = _series(neffs, lambda res, j: neffs[res][perms[res][j]])
    _print_series(corr)

    return neffs, perms, worst, naive, corr


def _series(neffs, pick):
    out = {}
    for i in range(1, len(RES_LIST)):
        a, b = RES_LIST[i - 1], RES_LIST[i]
        out[(a, b)] = [abs(pick(b, j) - pick(a, j)) for j in range(N_REPORT)]
    return out


def _print_series(series):
    for (a, b), vals in series.items():
        print(f"    {a:5.1f}->{b:5.1f}   max={max(vals):.2e}   "
              f"mean={np.mean(vals):.2e}")
    keys = list(series)
    maxes = [max(series[k]) for k in keys]
    mono = all(maxes[i] >= maxes[i + 1] for i in range(len(maxes) - 1))
    print(f"    monotonically shrinking: {mono}"
          + ("" if mono else "   <-- NOT converging on this tracking"))


results = {}
for engine, fn in (("FDTD / meep", solve_fdtd), ("EME / tidy3d", solve_eme)):
    try:
        results[engine] = run_engine(engine, fn)
    except Exception as exc:            # keep one engine's failure from
        print(f"  !! {engine} FAILED: {type(exc).__name__}: {exc}", flush=True)
        results[engine] = None
    print()

np.savez(
    os.path.join(output_dir, "check3_res_convergence.npz"),
    res_list=np.asarray(RES_LIST), n_solve=N_SOLVE, n_report=N_REPORT,
    cmp_res=CMP_RES,
    **{f"{eng.split('/')[0].strip()}_neff_{res}": results[eng][0][res]
       for eng in results if results[eng] for res in RES_LIST},
)
print(f"Saved -> {output_dir}/check3_res_convergence.npz")
