"""mmi_2x2_validation.md gap 4 (README.md check 1) — do FDTD and EME agree on
phase accumulation through the MMI body's own uniform cross-section?

Same question as `resolved/directional_coupler_validation.md` item 6, adapted
per `mmi_1x2_validation.md` gap 4's own suggestion: "solve a uniform stretch
of the MMI body at two lengths, and check that FDTD and EME agree on each
mode's accumulated phase." Uses `beef.solvers.fdtd.dbeta_from_phases()`
directly -- the same device-independent helper the coupler's own script and
its tests already exercise -- rather than re-deriving the fit.

Why this matters here specifically: the coupler's ~10% Δβ disagreement
between MPB's mode index and meep's time-stepped phase is what
`calibrate_body_neff=True` corrects, taking L₅₀ from ~5% off the full-device
reference to 0.11%. Neither MMI has ever tried that flag on a *measured*
basis -- gap 5's 0.82% (2×2) / 9.11% (1×2) hybrid-vs-FDTD numbers are
downstream evidence, not a direct measurement of the same mechanism. This
measures it directly, on nothing but a uniform body slice, so the answer
does not depend on the taper transitions or the cascade.

HOW
---
A cut entirely inside the MMI body (mm_length=28, so anything well clear of
the ±14 taper seams) is uniform -- neither engine can convert between guided
modes there. FDTD: `solve_section()` at several lengths, n_modes_in=
n_modes_out=2 (the same two lowest modes mmi_2x2_gap2_dbeta.py already
measures), gives a transfer matrix that must be diagonal (checked, not
assumed) with dphi(L) = arg(T_00) - arg(T_11) = Δβ·L + const. EME: the same
section_gds() polygons -> the same meow cross-section -> compute_modes(),
Δβ_EME read straight off the two eigenvalues.

WHAT COUNTS AS AN ANSWER
------------------------
ratio (FDTD/EME) ~ 1.00   this device's body accumulates phase the same way
                          in both engines -- calibrate_body_neff is not
                          indicated here on this mechanism (a real, useful
                          negative result).
ratio far from 1.00       the same discretization gap the coupler has is
                          present here too, and calibrate_body_neff is worth
                          trying on this device's own hybrid (never done).

Cost: 4 real FDTD `solve_section()` time-domain runs (not eigensolves) plus
one cheap EME mode solve. Comparable in size to a single interior FDTD
section of `modal_multiplexer_hybrid.py` (which took 15-35 min per section on
a similar body cross-section) -- batch job, 2 nodes, matching
`run_mmi_2x2_examples.batch`'s own resource shape for this device.

Run:
    python verify_mmi_2x2_engine_agreement.py
"""
from __future__ import annotations

import os
import subprocess
import time

import numpy as np
import meow as mw

from beef import MMI, EMESolver, FDTDSolver, SOI_220, extrude_polygons
from beef.solvers.fdtd import dbeta_from_phases

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "mmi_2x2_engine_agreement")
_CACHE_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

WL_MIN, WL_MAX, N_WL = 1.530, 1.565, 11
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

RESOLUTION = 40.0
N_MODES = 2   # the same pair mmi_2x2_gap2_dbeta.py measures

MMI_KWARGS = dict(
    n_inputs=2, n_outputs=2, mm_length=28.0, mm_width=2.8, waveguide_width=0.45,
    taper_length=4.0, taper_width=1.3, waveguide_offset=1.5, waveguide_height=0.22,
)
FDTD_KWARGS = dict(
    res=RESOLUTION, sim_width=5.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SOI_220, decay_tol=1e-11,
)
EME_KWARGS = dict(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=5.0, oxide_height=1.5, res=RESOLUTION,
    material_stack=SOI_220, mode_solver="meep",
)

# Deep inside the +/-14 um body (mm_length=28), same margin logic as the
# coupler's own choice: comfortably clear of the taper seams (>= 10 um here,
# vs the coupler's >= 15 um clearance on its 30 um coupling region), and small
# enough that dphi never wraps 2*pi between the shortest and longest length.
SECTION_LENGTHS = np.array([2.0, 4.0, 6.0, 8.0])


def make_mmi():
    return MMI(**MMI_KWARGS)


def git_sha():
    try:
        return subprocess.check_output(
            ["git", "-C", os.path.join(_ROOT, "BEEF"), "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "<not a git checkout>"


def simulate_length(mmi, L):
    solver = FDTDSolver(**FDTD_KWARGS)
    x0, x1 = -L / 2.0, +L / 2.0
    half_body = mmi.mm_length / 2.0
    assert abs(x0) < half_body and abs(x1) < half_body, (
        f"section [{x0}, {x1}] escapes the {mmi.mm_length} um body -- "
        f"it would include the taper, exactly what this must exclude"
    )
    print(f"[ea] section L={L} um, x in [{x0:+.2f}, {x1:+.2f}] "
          f"(body +/-{half_body:.1f} um), {N_MODES} modes, "
          f"{len(wavelengths)} wavelengths", flush=True)
    t0 = time.time()
    T = solver.solve_section(
        mmi, wavelengths, x0, x1,
        n_modes_in=N_MODES, n_modes_out=N_MODES,
        save_dir=_CACHE_DIR,
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
    print(f"[ea] mmi_2x2 body engine agreement, mm_length={mmi.mm_length}, "
          f"lengths={SECTION_LENGTHS.tolist()}", flush=True)

    dphi_all = np.zeros((len(SECTION_LENGTHS), len(wavelengths)))
    worst_off = 0.0
    for i, L in enumerate(SECTION_LENGTHS):
        dphi_all[i], off = simulate_length(mmi, float(L))
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
