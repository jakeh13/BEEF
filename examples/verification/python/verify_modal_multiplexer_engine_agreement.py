"""deferred/modal_multiplexer_validation.md gap 5 (README.md check 1) -- do
FDTD and EME agree on phase accumulation through the mux's own body
cross-section?

Same question as resolved/directional_coupler_validation.md item 6 and its
two MMI siblings (resolved/mmi_1x2_validation.md gap 4,
resolved/mmi_2x2_validation.md gap 4), adapted per this device's own gap 5:
"Adapt verify_item_6_engine_agreement.py onto a uniform stretch of the
body." Reuses beef.solvers.fdtd.dbeta_from_phases() directly, the same
device-independent helper both siblings use.

WHICH BODY -- gap 1's own interface audit already established this device's
two MMI bodies (1x3, mm_length=16.92; 3x3, mm_length=34.24) share an
IDENTICAL cross-section -- 2.8 um x 0.16 um SOI_160 -- verified directly, not
assumed, before treating one interface as representative of all four. This
script therefore checks only the shorter body (MMI 1, the 1x3 stage) as a
standalone MMI leaf geometry, the same "cut a uniform interior slice"
approach the coupler/MMI scripts use, rather than routing through the full
CompoundGeometry mux (irrelevant here -- section_gds() is local to the leaf,
and the mux's outer tapers/connector never enter a cut this deep inside one
body).

HOW
---
Same mechanics as the two MMI siblings: FDTD solve_section() at several
lengths well inside the body (clear of both taper seams), n_modes_in=
n_modes_out=2 -- the two lowest guided body supermodes, the generic
convention every engine-agreement script in this repo uses (independent of
this device's own N_MODES=16 interface basis, which check 2/gap 1 already
covers separately). A uniform region cannot convert between them, so
solve_section()'s transfer matrix must be diagonal (checked, not assumed);
dphi(L) = arg(T_00) - arg(T_11) = Delta_beta*L + const, fit across lengths
via dbeta_from_phases(). EME: the same section_gds() polygons -> the same
meow cross-section -> compute_modes(), Delta_beta_EME read straight off the
two eigenvalues.

Resolution -- pinned at RESOLUTION = 40, the value modal_multiplexer_hybrid.py
actually ships (this checks the shipped configuration, not a converged one;
gap 2 already measured that 40 does not converge the FULL 16-mode basis, but
the two lowest, most-confined modes checked here are the least sensitive to
that -- same reasoning resolved/mmi_2x2_validation.md's own gap 4 used at the
same resolution).

WHAT COUNTS AS AN ANSWER
------------------------
ratio (FDTD/EME) ~ 1.00   this device's body accumulates phase the same way
                          in both engines -- calibrate_body_neff is not
                          indicated by this mechanism (matches both MMIs'
                          own findings, ratio ~1 to within a few %).
ratio far from 1.00       the coupler's ~10% discretization gap is present
                          here too, and calibrate_body_neff is worth trying
                          on this device's hybrid chain (never attempted).

Cost -- learned up front from resolved/mmi_2x2_validation.md gap 4's own
history rather than rediscovered here: its FIRST attempt at this exact check
(job 11591249) TIMED OUT at 5h having finished only 2 of 4 lengths (L=2 took
96 min, L=4 took 183 min -- cost grows faster than linearly with L), and had
to be reshaped into an array-per-length job before it could complete
(11596845 + merge 11596851, worst length 5:07:14). This script ships in that
already-reshaped form from the start: one SLURM array task per section
length (SLURM_ARRAY_TASK_ID 0..3), each simulating just its own length and
persisting to the shared on-disk cache; a separate merge/serial step (no
array index) reads every length back with cache_only=True, fits
Delta_beta_FDTD, runs the cheap EME half, and saves the result. Same 2-node
resource shape as the MMI siblings' own jobs (this body's mm_width=2.8 um
matches the 2x2's own 2.8 um -- expect comparable cost, if anything somewhat
more given this device's weaker 160 nm-core confinement).

Deliberately bare `python`, NOT `srun --mpi=pmi2`, matching both MMI
siblings and for the same reason recorded in
batch/verification/run_mmi_1x2_engine_agreement.batch: this script has no
MPI rank guards, so parallel ranks would race on one np.savez() call and one
output directory. Both MMI siblings proved a single meep process finishes
this size of problem inside an 8h budget; dropping the launcher here is the
proven-safe configuration -- NOT the missing-launcher bug that broke
modal_multiplexer_hybrid.py's own reshape (that device's batch script invokes
meep's actual production solve_section() workload, with N_MODES=16 and
rank-guarded output, across a 96-rank allocation through a target that
requires srun -- an entirely different, much larger computation from this
one).

Run (single process, all 4 lengths serially):
    python verify_modal_multiplexer_engine_agreement.py

Run (recommended): submit the array job for the 4 lengths, wait for all 4 to
complete, then submit the merge step (same script, no array index) to
fit/compare/save from cache. See
batch/verification/run_modal_multiplexer_engine_agreement.batch and
..._merge.batch.
"""
from __future__ import annotations

import os
import subprocess
import time

import numpy as np
import meow as mw

from beef import MMI, MMIPort, EMESolver, FDTDSolver, SOI_160, extrude_polygons
from beef.solvers.fdtd import dbeta_from_phases

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "modal_multiplexer_engine_agreement")
_CACHE_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# O-band, byte-identical to modal_multiplexer_hybrid.py's own WL_MIN/WL_MAX/
# N_WL -- this checks the shipped device's own band, not an arbitrary one.
WL_MIN, WL_MAX, N_WL = 1.26, 1.36, 11
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

RESOLUTION = 40.0   # shipped value -- see module docstring
N_MODES = 2         # the two lowest body supermodes; generic engine-agreement
                    # convention (coupler even/odd, both MMIs), independent of
                    # this device's own N_MODES=16 interface basis

# MMI 1 (1x3 stage) -- byte-for-byte the same construction
# modal_multiplexer_hybrid.py/modal_multiplexer_fdtd.py use, kept in sync by
# hand (both those scripts already do this; short enough that a shared helper
# module is not worth it for a third call site either).
TAPER_L = 2.0
WO = 1.0375
MM_WIDTH = 2.8

BRANCH_SPECS = [
    (+WO, 0.33, 0.725),
    (0.0, 0.35, 1.15),
    (-WO, 0.35, 0.725),
]


def _branch_ports():
    return [MMIPort(y=y, waveguide_width=ww, taper_length=TAPER_L, taper_width=tw)
            for y, ww, tw in BRANCH_SPECS]


def make_mmi1():
    return MMI(
        n_inputs=1, n_outputs=3,
        mm_length=16.92, mm_width=MM_WIDTH,
        inputs=[MMIPort(y=0.0, waveguide_width=0.7, taper_length=TAPER_L, taper_width=1.4)],
        outputs=_branch_ports(),
        waveguide_height=0.16,
    )


FDTD_KWARGS = dict(
    res=RESOLUTION, sim_width=5.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SOI_160, decay_tol=1e-11,
)
EME_KWARGS = dict(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=5.0, oxide_height=1.5, res=RESOLUTION,
    material_stack=SOI_160, mode_solver="meep",
)

# Deep inside MMI 1's +/-8.46 um body (mm_length=16.92), same margin logic as
# both siblings: at L=8 this leaves 4.46 um clear of the taper seam --
# thinner than the 2x2's 10 um margin (mm_length=28) but still ~6x
# modal_multiplexer_hybrid.py's own BODY_STUB=0.75 um safety margin.
SECTION_LENGTHS = np.array([2.0, 4.0, 6.0, 8.0])


def git_sha():
    try:
        return subprocess.check_output(
            ["git", "-C", os.path.join(_ROOT, "BEEF"), "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "<not a git checkout>"


def simulate_length(mmi, L, cache_only=False):
    solver = FDTDSolver(**FDTD_KWARGS)
    x0, x1 = -L / 2.0, +L / 2.0
    half_body = mmi.mm_length / 2.0
    assert abs(x0) < half_body and abs(x1) < half_body, (
        f"section [{x0}, {x1}] escapes the {mmi.mm_length} um body -- "
        f"it would include the taper, exactly what this must exclude"
    )
    print(f"[ea] section L={L} um, x in [{x0:+.2f}, {x1:+.2f}] "
          f"(body +/-{half_body:.1f} um), {N_MODES} modes, "
          f"{len(wavelengths)} wavelengths"
          + (" (cache_only)" if cache_only else ""), flush=True)
    t0 = time.time()
    T = solver.solve_section(
        mmi, wavelengths, x0, x1,
        n_modes_in=N_MODES, n_modes_out=N_MODES,
        save_dir=_CACHE_DIR, cache_only=cache_only,
    )
    dt = time.time() - t0
    diag = np.abs(np.stack([T[:, 0, 0], T[:, 1, 1]], axis=1))
    off = np.abs(np.stack([T[:, 0, 1], T[:, 1, 0]], axis=1))
    dphi = np.angle(T[:, 0, 0]) - np.angle(T[:, 1, 1])
    print(f"[ea] solved in {dt/60:.1f} min | |diag| {diag.min():.3f}-{diag.max():.3f} "
          f"| max |off-diag| {off.max():.4f}", flush=True)
    return dphi, off.max()


def eme_dbeta(mmi):
    L = float(SECTION_LENGTHS[0])
    eme = EMESolver(**EME_KWARGS)
    gds = mmi.section_gds(-L / 2.0, +L / 2.0, eme.oxide_width)
    structs = extrude_polygons(gds, eme.material_stack.extrusion_rules)
    cells = mw.create_cells(structures=structs, mesh=eme._build_mesh(), Ls=[L])
    neff = np.zeros((len(wavelengths), N_MODES))
    for i, wl in enumerate(wavelengths):
        cs = mw.CrossSection.from_cell(cell=cells[0],
                                        env=mw.Environment(wl=float(wl), T=eme.temperature))
        modes = eme.compute_modes(cs, N_MODES)
        neff[i] = [float(np.real(m.neff)) for m in modes]
    assert np.all(neff[:, 0] > neff[:, 1]), (
        f"EME modes not in descending n_eff order -- the identification below "
        f"would be wrong. neff = {neff}"
    )
    dbeta = 2 * np.pi * (neff[:, 0] - neff[:, 1]) / wavelengths
    return dbeta, neff


def main():
    os.makedirs(output_dir, exist_ok=True)
    mmi = make_mmi1()

    array_idx = os.environ.get("SLURM_ARRAY_TASK_ID")
    if array_idx is not None:
        # One array task, one length: simulate and let solve_section's own
        # save_dir cache persist the result -- no fit/ratio computed here,
        # the merge step (run with no array index) reads every length back
        # out of that same cache.
        i = int(array_idx)
        L = float(SECTION_LENGTHS[i])
        print(f"[ea] array task {i}: modal_multiplexer MMI1 body engine "
              f"agreement, mm_length={mmi.mm_length}, L={L} um "
              f"(of {SECTION_LENGTHS.tolist()})", flush=True)
        simulate_length(mmi, L)
        print(f"[ea] array task {i} done -- L={L} um cached under "
              f"{_CACHE_DIR}", flush=True)
        return

    print(f"[ea] modal_multiplexer MMI1 body engine agreement, "
          f"mm_length={mmi.mm_length}, lengths={SECTION_LENGTHS.tolist()} "
          f"(merge/serial step -- cache_only=True, expects each length "
          f"already computed by its own array task)", flush=True)

    dphi_all = np.zeros((len(SECTION_LENGTHS), len(wavelengths)))
    worst_off = 0.0
    for i, L in enumerate(SECTION_LENGTHS):
        dphi_all[i], off = simulate_length(mmi, float(L), cache_only=True)
        worst_off = max(worst_off, off)

    dbeta_fdtd = np.zeros(len(wavelengths))
    resid = np.zeros(len(wavelengths))
    for w in range(len(wavelengths)):
        dbeta_fdtd[w], _, resid[w] = dbeta_from_phases(SECTION_LENGTHS, dphi_all[:, w])

    print("[ea] running the EME half...", flush=True)
    dbeta_eme, neff = eme_dbeta(mmi)

    ratio = dbeta_fdtd / dbeta_eme
    print(f"\n[ea] worst |off-diagonal| across all lengths/wavelengths: {worst_off:.4f}")
    print(f"[ea] worst FDTD fit residual (rad): {resid.max():.4f}")
    print(f"\n{'wl':>7}  {'dbeta_fdtd':>11}  {'dbeta_eme':>11}  {'ratio':>8}")
    for w, wl in enumerate(wavelengths):
        print(f"{wl:7.4f}  {dbeta_fdtd[w]:11.5f}  {dbeta_eme[w]:11.5f}  {ratio[w]:8.4f}")
    print(f"\n[ea] mean ratio (FDTD/EME): {ratio.mean():.4f}  "
          f"({'~1 -- no calibration mechanism found here' if abs(ratio.mean()-1)<0.02 else 'off from 1 -- same mechanism as the coupler may be present'})")

    np.savez(os.path.join(output_dir, "engine_agreement.npz"),
             section_lengths=SECTION_LENGTHS, wavelengths=wavelengths,
             dphi=dphi_all, dbeta_fdtd=dbeta_fdtd, dbeta_eme=dbeta_eme,
             neff=neff, ratio=ratio, resid=resid, worst_off_diag=worst_off,
             n_modes=N_MODES, resolution=RESOLUTION, beef_git=git_sha())
    print(f"\n[ea] saved -> {output_dir}/engine_agreement.npz")


if __name__ == "__main__":
    main()
