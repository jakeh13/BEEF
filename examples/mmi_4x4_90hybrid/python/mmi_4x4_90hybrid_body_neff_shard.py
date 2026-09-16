"""4x4 90 hybrid -- WAVELENGTH-PARALLEL body eigenmode solve.

WHY THIS EXISTS. mmi_4x4_90hybrid_hybrid.py's body step costs ~416 s per
body-mode x N_MODES=16 = ~1.85 h PER WAVELENGTH, and nothing about it is
parallel: solve_section_length_sweep() walks `for wl_idx, wl in
enumerate(wavelengths)` on one core. Going 3 -> 11 wavelengths therefore takes
it from ~5.5 h to ~20 h of strictly serial work. MPI does not help -- meep
solves the eigenproblem redundantly on every rank (measured: 8 ranks 3110.7 s
vs 24 ranks 3055.2 s, jobs 12218118/12218119).

WHY SHARDING IS EXACT, not an approximation. Two facts from the library:

  1. `EMESolver._extrusion_rules()` returns `material_stack.extrusion_rules`
     UNCHANGED whenever the solver is dispersive (this one is -- `dispersive`
     is left None, and `_resolved_dispersive()` maps None -> True on the EME
     side deliberately). It only consults `wavelengths` in the FROZEN branch,
     where it would build materials at `mean(wavelengths)`. So the structures,
     mesh and cells handed to the mode solver do not depend on which OTHER
     wavelengths are in the array.
  2. The per-wavelength body of that loop builds its own `mw.Environment(wl=...)`
     and its own CrossSection from those wavelength-independent cells. Row
     `neff_all[k]` is therefore a pure function of `wavelengths[k]`.

     ==> IF EITHER STOPS BEING TRUE, THIS SCRIPT IS WRONG. In particular, do
     not run it against a solver built with `dispersive=False`; --verify
     refuses that case explicitly.

So the tier-3 cache payload (`neff_all`, shape [n_wl, n_modes]) can be
computed one wavelength per process and CONCATENATED. The merged array is
saved under the key the real 11-wavelength call will build, and that call then
takes tier 3 ("length-independent neff reuse"), reconstructing every length
analytically via _propagate_diagonal() with no mode solve at all.

WHAT IS NOT ASSUMED. The tier-3 reconstruction itself is not something this
script invents -- it is the library's own cache path, and every shard runs
eme.py's `_DIAGONAL_RECONSTRUCTION_TOL` self-check before its neff is allowed
into a cache at all. --merge refuses to write if any shard failed that check.

MODES
  --verify      Rebuild the 3-wavelength tier-3 key from the geometry below and
                look it up in the PRODUCTION cache. A hit proves this script's
                geometry/solver reconstruction is byte-identical to the
                example's, since the key is a hash over all of it. Costs
                nothing. RUN THIS FIRST -- a miss means the merge would write
                an entry the real run can never find.
  --shard K     Solve wavelength K of the 11-point grid. ~1.85 h, one core.
  --check       Compare shards 0/5/10 (1.290/1.310/1.330 um) against the rows
                of the CACHED 3-wavelength neff. Those three wavelengths are
                common to both grids, so this is a free, direct, empirical
                test of the independence argument above. Not a proxy.
  --merge       Stack all 11 shards and write the 11-wavelength tier-3 entry
                into the production cache. Refuses unless --check passed.
"""
import argparse
import os
import sys

import numpy as np

from beef import MMI, EMESolver, SOI_220
from beef.cache import build_run_config, lookup, load_result, save_run
from beef.materials.stack import resolved_material_config
from beef.solvers.eme import _neff_fingerprint

_HERE = os.path.dirname(os.path.abspath(__file__))
SAVE_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))
SHARD_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                         "example_sim_cache_shards", "mmi4x4_body_neff"))

# --- pinned to mmi_4x4_90hybrid_hybrid.py; --verify is what proves it ------
RESOLUTION = 40.0
WL_MIN, WL_MAX = 1.29, 1.33
N_WL = 11
N_MODES = 16
WAVEGUIDE_WIDTH, MM_WIDTH = 0.35, 8.5
TAPER_WIDTH, TAPER_LENGTH = 1.30, 4.0
MM_LENGTH = float(os.environ.get("MMI4X4_MM_LENGTH", 159.30))
BODY_STUB = 0.25

WL11 = np.linspace(WL_MIN, WL_MAX, N_WL)
WL3 = np.linspace(WL_MIN, WL_MAX, 3)

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=MM_WIDTH / 4.0,
    waveguide_height=0.22, use_neighbor_aperture=True,
)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=MM_WIDTH + 3.0, oxide_height=2.0,
    res=RESOLUTION, material_stack=SOI_220, mode_solver="meep",
)
_half_body = mmi.mm_length / 2
X0, X1 = -_half_body + BODY_STUB, +_half_body - BODY_STUB


def neff_config(wavelengths):
    """The tier-3 key, rebuilt exactly as eme.py:1085-1090 builds it."""
    mat = resolved_material_config(body_solver.material_stack, wavelengths,
                                   body_solver._resolved_dispersive())
    return build_run_config(
        body_solver, "solve_section_neff",
        dict(x0=X0, x1=X1, wavelengths=wavelengths, n_modes=N_MODES,
             stub_left=0.0, stub_right=0.0, geometry_kwargs={},
             neff_override=_neff_fingerprint(None), **mat),
        mmi,
    )


def shard_path(k):
    return os.path.join(SHARD_DIR, f"neff_wl{k:02d}.npy")


def do_verify():
    if not body_solver._resolved_dispersive():
        sys.exit("REFUSING: solver is frozen (dispersive=False), so the cross-section "
                 "depends on mean(wavelengths) and per-wavelength shards are NOT "
                 "independent. See this module's docstring.")
    print(f"body section x in [{X0:+.3f}, {X1:+.3f}] um, {N_MODES} modes, res {RESOLUTION}")
    hit = lookup(SAVE_DIR, neff_config(WL3))
    if hit is None:
        print("\n3-wl tier-3 key: MISS")
        print("  Either the geometry reconstruction differs from the example's, or the")
        print("  3-wl run predates tier 3. Do NOT merge until this hits -- check by")
        print("  re-running the example at MMI4X4_HYB_NWL=3 (all other steps cache-hit).")
        return 1
    neff3 = load_result(hit)
    print(f"\n3-wl tier-3 key: HIT at {hit}")
    print(f"  payload {neff3.shape} -- geometry/solver reconstruction is byte-identical")
    print(f"  Re(neff) @ 1.310 um, first 4 modes: {np.real(neff3[1, :4])}")
    print(f"\n11-wl target key exists already: {lookup(SAVE_DIR, neff_config(WL11)) is not None}")
    return 0


def do_shard(k):
    wl = np.array([WL11[k]])
    os.makedirs(SHARD_DIR, exist_ok=True)
    out = shard_path(k)
    if os.path.exists(out):
        print(f"shard {k} (wl={wl[0]:.4f}) already done -> {out}")
        return 0
    print(f"shard {k}: solving {N_MODES} modes at wl={wl[0]:.4f} um "
          f"(expect ~1.85 h, single core)", flush=True)
    # One length only: neff does not depend on it, and the full 7601-point grid
    # would add a compute_s_matrix() per length for nothing. The real run
    # reconstructs all 7601 analytically from the merged neff.
    body_solver.solve_section_length_sweep(
        mmi, wl, X0, X1, np.array([100.0]),
        n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=SHARD_DIR,
    )
    hit = lookup(SHARD_DIR, neff_config(wl))
    if hit is None:
        sys.exit(f"shard {k}: solve finished but NO tier-3 entry was written -- eme.py's "
                 f"diagonal-reconstruction self-check must have failed. Refusing to "
                 f"emit a shard; this wavelength is not eligible for neff reuse.")
    neff = load_result(hit)
    np.save(out, neff)
    print(f"shard {k}: wrote {neff.shape} -> {out}")
    return 0


def _load_all():
    missing = [k for k in range(N_WL) if not os.path.exists(shard_path(k))]
    if missing:
        sys.exit(f"missing shards {missing} -- run --shard for each first")
    return np.concatenate([np.load(shard_path(k)) for k in range(N_WL)], axis=0)


def do_check():
    hit = lookup(SAVE_DIR, neff_config(WL3))
    if hit is None:
        sys.exit("no cached 3-wl neff to check against -- run --verify")
    neff3 = load_result(hit)
    ok = True
    print("shard vs cached 3-wl neff (wavelengths common to both grids):")
    for row, k in enumerate((0, 5, 10)):
        if not os.path.exists(shard_path(k)):
            print(f"  shard {k:2d} (wl={WL11[k]:.4f}): NOT YET DONE")
            ok = False
            continue
        d = np.max(np.abs(np.load(shard_path(k))[0] - neff3[row]))
        flag = "EXACT" if d == 0 else ("ok" if d < 1e-12 else "MISMATCH")
        print(f"  shard {k:2d} (wl={WL11[k]:.4f}) vs 3-wl row {row}: max|dneff| = {d:.3e}  {flag}")
        ok = ok and d < 1e-12
    print("\nindependence argument CONFIRMED empirically" if ok else
          "\nCHECK FAILED -- do not merge")
    return 0 if ok else 1


def do_merge(force=False):
    if not force and do_check() != 0:
        sys.exit("refusing to merge (use --force only if you know why --check failed)")
    neff_all = _load_all()
    assert neff_all.shape == (N_WL, N_MODES), neff_all.shape
    cfg = neff_config(WL11)
    if lookup(SAVE_DIR, cfg) is not None:
        print("11-wl tier-3 entry already present -- nothing to do")
        return 0
    folder = save_run(SAVE_DIR, cfg, neff_all)
    print(f"merged {neff_all.shape} -> {folder}")
    print("the 11-wavelength example run will now take tier 3 (no body mode solve)")
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--verify", action="store_true")
    p.add_argument("--shard", type=int)
    p.add_argument("--check", action="store_true")
    p.add_argument("--merge", action="store_true")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    if a.verify:
        sys.exit(do_verify())
    if a.shard is not None:
        sys.exit(do_shard(a.shard))
    if a.check:
        sys.exit(do_check())
    if a.merge:
        sys.exit(do_merge(a.force))
    p.error("pick one of --verify/--shard K/--check/--merge")
