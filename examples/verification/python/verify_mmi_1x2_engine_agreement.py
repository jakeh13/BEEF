"""mmi_1x2_validation.md gap 4 (README.md check 1) — do FDTD and EME agree on
phase accumulation through the MMI body's own uniform cross-section?

Sibling of verify_mmi_2x2_engine_agreement.py -- same question, same recipe,
ported to this device's geometry and its own established wavelength grid
(WL_MIN/MAX/N_WL = 1.26/1.36/11, RESOLUTION=40 -- the values already pinned
by mmi_1x2_fdtd.py / mmi_1x2_hybrid.py / verify_mmi_1x2_hybrid_vs_fdtd.py).
Read that script's docstring for the full argument (coupler item 6's lesson,
`calibrate_body_neff`, why this measures the mechanism directly instead of
inferring it from gap 5's cascade result).

Why this device specifically, now: this is the worst-converging of the three
audited (gap 2) and gap 5 already measured its uncalibrated hybrid 9.11% off
the full-device FDTD reference -- 3x the 2x2's residual -- but nothing has
ever measured whether any of that is the same FDTD/EME body-phase mechanism
`calibrate_body_neff` fixed for the coupler, or something else (gap 3's
former N_MODES mismatch, gap 2's slower n_eff convergence, or genuinely
structural). This is the one measurement that isolates the mechanism from
the cascade.

HOW
---
A cut entirely inside the MMI body (mm_length=4.5, taper_length=1.0 each
side -> body extends to +/-2.25) is uniform -- neither engine can convert
between guided modes there. Section lengths are scaled down from the 2x2's
own [2,4,6,8] um choice to fit this device's much shorter body: [0.5, 1.0,
1.5, 2.0] um, keeping >=1.25 um clearance from each taper seam at the
longest length -- the same clearance verify_mmi_1x2_gap2_dbeta.py's own
SECTION=(-1,1) already uses at L=2.0.

Wrap safety: mmi_1x2_validation.md's own item-3 table gives n_eff[0]=2.693784,
n_eff[1]=2.634501 at res 40, wl=1.31 -- dbeta = 2*pi/1.31*(2.693784-2.634501)
~= 0.284 rad/um, so even the longest 2.0 um step accrues ~0.57 rad, nowhere
near a wrap (same margin logic dbeta_from_phases()'s own docstring uses for
the coupler's 0.095 rad/um case).

FDTD: `solve_section()` at each length, n_modes_in=n_modes_out=2 (the same
pair gap 2 measures), gives a transfer matrix that must be diagonal (checked,
not assumed) with dphi(L) = arg(T_00) - arg(T_11) = Δβ·L + const. EME: the
same section_gds() polygons -> the same meow cross-section ->
compute_modes(), Δβ_EME read straight off the two eigenvalues.

WHAT COUNTS AS AN ANSWER
------------------------
ratio (FDTD/EME) ~ 1.00   this device's body accumulates phase the same way
                          in both engines -- calibrate_body_neff is not
                          indicated on this mechanism; gap 5's 9.11% residual
                          is something else (mode identity, N_MODES, or
                          structural).
ratio far from 1.00       the coupler's discretization gap is present here
                          too, and re-running gap 5's comparison with
                          calibrate_body_neff=True (already generated, jobs
                          11589770/11589779) is testing the right mechanism.

Cost: 4 real FDTD `solve_section()` time-domain runs (not eigensolves) on a
smaller cross-section (SOI_160, mm_width=2.0) than the 2x2's (SOI_220,
mm_width=2.8) -- same resource shape as run_mmi_1x2_examples.batch for this
device regardless, matching the 2x2 engine-agreement job's own conservative
choice.

**Array/merge shape, adopted 2026-08-01 before this was ever submitted.**
The 2x2's sibling script had to learn this the hard way -- job 11591249 timed
out at 5h having finished only 2 of its 4 lengths, because `solve_section`
does min(n_modes_in, n_modes_out) FDTD runs per length and cost grows faster
than linearly with L. Reshaped, it closed on the first try (11596845 [0-3] +
11596851). This device's sections are much shorter (longest L=2.0 um vs the
2x2's 8.0) so a single serial pass would very likely fit in one slot -- the
array shape is adopted anyway because it is free, and because a per-length
timeout then costs one length instead of all four:

1. Each length can be run as its own SLURM array task
   (`SLURM_ARRAY_TASK_ID` in {0,1,2,3} indexes `SECTION_LENGTHS`) -- see
   `run_mmi_1x2_engine_agreement.batch`. A task only calls `simulate_length`
   for its one length and exits; it does not compute the final fit/ratio,
   since it never has the other tasks' results in memory.
2. Running with no `SLURM_ARRAY_TASK_ID` set (the merge/serial step) passes
   `cache_only=True` down to every `solve_section` call: it reuses `beef`'s
   on-disk cache (`save_dir=_CACHE_DIR`, already wired into
   `simulate_length`) instead of re-simulating, and raises `CacheMissError`
   loudly if any length's array task has not completed -- rather than
   silently kicking off a multi-hour FDTD solve serially.

Note this runs under a plain `python`, NOT `srun` -- there are no rank guards
in here, so 48 ranks would race on one `np.savez` and one output dir. That is
how the modal multiplexer's reshape broke (jobs 11597177 [0-2], three 5h
timeouts at 1 rank instead of 96). The 2x2 sibling ran serial under a 2-node
allocation and completed; same choice here, deliberately.

Run (either shape):
    python verify_mmi_1x2_engine_agreement.py            # all 4 lengths, serial
    SLURM_ARRAY_TASK_ID=2 python verify_..._agreement.py  # just L=1.5
"""
from __future__ import annotations

import os
import subprocess
import time

import numpy as np
import meow as mw

from beef import MMI, EMESolver, FDTDSolver, SOI_160, extrude_polygons
from beef.solvers.fdtd import dbeta_from_phases

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "mmi_1x2_engine_agreement")
_CACHE_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

WL_MIN, WL_MAX, N_WL = 1.26, 1.36, 11
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

RESOLUTION = 40.0
N_MODES = 2   # the same pair verify_mmi_1x2_gap2_dbeta.py measures

MMI_KWARGS = dict(
    n_inputs=1, n_outputs=2, mm_length=4.5, mm_width=2.0, waveguide_width=0.35,
    taper_length=1.0, taper_width=0.8, waveguide_offset=1.0, waveguide_height=0.16,
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

# Deep inside the +/-2.25 um body (mm_length=4.5, taper_length=1.0): scaled
# down from the 2x2's [2,4,6,8] to fit, keeping >=1.25 um clearance from each
# taper seam at the longest length -- the same clearance
# verify_mmi_1x2_gap2_dbeta.py's own SECTION=(-1,1) uses at L=2.0.
SECTION_LENGTHS = np.array([0.5, 1.0, 1.5, 2.0])


def make_mmi():
    return MMI(**MMI_KWARGS)


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
          f"(body +/-{half_body:.2f} um), {N_MODES} modes, "
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
    mmi = make_mmi()

    array_idx = os.environ.get("SLURM_ARRAY_TASK_ID")
    if array_idx is not None:
        # One array task, one length: simulate and let solve_section's own
        # save_dir cache persist the result -- no fit/ratio computed here,
        # the merge step (run with no array index) reads every length back
        # out of that same cache.
        i = int(array_idx)
        L = float(SECTION_LENGTHS[i])
        print(f"[ea] array task {i}: mmi_1x2 body engine agreement, "
              f"mm_length={mmi.mm_length}, L={L} um "
              f"(of {SECTION_LENGTHS.tolist()})", flush=True)
        simulate_length(mmi, L)
        print(f"[ea] array task {i} done -- L={L} um cached under "
              f"{_CACHE_DIR}", flush=True)
        return

    print(f"[ea] mmi_1x2 body engine agreement, mm_length={mmi.mm_length}, "
          f"lengths={SECTION_LENGTHS.tolist()} (merge/serial step -- "
          f"cache_only=True, expects each length already computed by its "
          f"own array task)", flush=True)

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
