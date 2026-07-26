"""TEMPORARY — verification harness for dc_validation.md items 1 and 2.

Delete this folder (and the Makefile's verify_items_1_2/clean_verification
targets, and batch/verification/) once items 1-2 are signed off. It exists to
answer one question: the unit tests in BEEF/tests/test_fdtd_freq_grid.py prove
the new helpers behave, but nothing has proved the *physics* came out right,
because no FDTD run has happened since the fix.

Item 1 (monitor frequency grid) and item 2 (source bandwidth) are described in
BEEF/validation_markdown/dc_validation.md. In short: every DFT monitor used to be built from meep's
(fcen, df, nfreq) triple with df = (1/lam_min - 1/lam_max) * 0.5. Meep spreads
nfreq points evenly over [fcen - df/2, fcen + df/2] ascending *in frequency*, so
a 1.50-1.60 um 11-point request actually sampled 1.5254-1.5754 um, backwards,
and was then indexed against the caller's ascending-lambda array. Only the
centre point was labelled correctly.

Nine checks, in two phases.

Static (checks 1-6, no meep run, seconds — safe on a login node with
--static-only):

    1. monitor grid strictly ascending in frequency (meep requires it)
    2. 1/mon_freqs equals the requested wavelengths exactly
    3. wl_order round-trips a monitor-ordered array back to caller order
    4. source reaches >=10% of peak at every requested wavelength
    5. the pre-fix span*0.5 really would have starved the edges (<1e-6) —
       confirms the regression was real, not theoretical
    6. a single-wavelength request still gets a usable (nonzero) fwidth

Physical (checks 7-9, one 3-D FDTD solve): this rebuilds the *exact* config of
an archived pre-fix run —

    archive/example_output_20260724/directional_coupler/
        directional_coupler_fdtd_resolution_sweep/resolution_sweep.npz

whose res-40 row is the same device (coupling_length 13.5 um), same solver
settings, same 11 wavelengths, run before the fix. Its T_cross *falls* 0.5151 ->
0.3448 across the labelled 1.50 -> 1.60 um, i.e. coupling appears to weaken with
wavelength, which is backwards for a directional coupler.

    7. the new T_cross(lambda) increases with wavelength
    8. the archived row, re-plotted on the wavelengths meep *actually* sampled
       (1.5754 -> 1.5254 um, reversed), lands on the new run's curve
    9. power is conserved (through + cross ~ 1) at every wavelength — the run
       has to be sane before checks 7-8 mean anything

Check 8 is the strong one: it says the old data was ONLY mislabelled — same
physics, wrong axis — which is precisely the claim item 1 makes. If item 1 were
wrong, or the bug were something other than a relabelling, the two curves would
not coincide.

Outputs (example_output/verification/items_1_2/):

    verification.json            provenance + per-check PASS/FAIL + the numbers
    verification.npz             every array behind the plots
    freq_grid.png                labelled vs actually-sampled wavelength, old/new
    source_coverage.png          source amplitude across the band, old/new
    dispersion_before_after.png  the headline: pre-fix on its labelled axis,
                                 pre-fix on its true axis, post-fix run, plus a
                                 residual panel for check 8
    power_conservation.png       through/cross/total vs wavelength

Everything printed is also stored, and vice versa: the .out log alone is enough
to judge the run, and verification.json alone is enough to reproduce it. That is
deliberate — dc_validation.md's complaint about the archived reference is that
it stored bare arrays with no config, no solver fingerprint and no provenance,
which is why its own validity had to be argued rather than read off.

Run:
    # full (needs the queue — one 3-D FDTD solve, ~1 h on 48 ranks):
    sbatch batch/verification/run_verify_items_1_2.batch

    # static checks only, safe on a login node:
    OMP_NUM_THREADS=1 python .../verify_items_1_2.py --static-only
"""

import argparse
import json
import logging
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import beef
from beef import DirectionalCoupler, FDTDSolver, SOI_220
from beef.solvers.fdtd import check_source_coverage, monitor_freq_grid, source_fwidth

# beef warns through the logging module — check_source_coverage() in particular.
# Without this, a starved-source warning would never reach the batch log, which
# would defeat the point of having added it.
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "items_1_2")

# Deliberately NO save_dir: this run is not cached, in either direction.
#
# Reading: the job of this script is to prove the *current* code produces
# correct physics. Any path that can return a stored result is a way for the
# check to pass on data it did not compute. The fingerprint did change when
# item 2 added source_bandwidth_factor (all 54 cached FDTDSolver runs lack the
# key and can no longer be hit, and no whole-device FDTD solve() was ever
# cached anyway) — but a verification that depends on that reasoning being
# right is not a verification.
#
# Writing: solve() calls save_run(), which is not multi-rank safe — plain
# truncating writes of config.json/result.pkl plus an append to a shared
# index.jsonl, executed by all 48 MPI ranks at once. That is the same class of
# bug as dc_validation.md item 4 (still open), which cost two stub-sweep tasks
# a 5 h wall clock. Not worth risking a shared cache to save one re-run.

# The archived pre-fix reference. Its res-40 row is the before picture.
_REF_NPZ = os.path.join(
    _ROOT, "archive", "example_output_20260724", "directional_coupler",
    "directional_coupler_fdtd_resolution_sweep", "resolution_sweep.npz")
_REF_RES = 40.0

# Pinned to match that archived run exactly — do not "improve" these, the
# whole point is that the only difference between the two runs is the fix.
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
COUPLING_LENGTH = 13.5   # um — the ~3 dB length, same as directional_coupler_fdtd.py

DEVICE_KWARGS = dict(
    waveguide_width=0.45, gap=0.25, coupling_length=COUPLING_LENGTH,
    bend_length=9.0, waveguide_offset=3.0, bend_type="euler",
    waveguide_height=0.22,
)
SOLVER_KWARGS = dict(
    res=RESOLUTION, sim_width=6.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    num_input_modes=1, num_output_modes=1,
    output_monitor_centers=[-1.5, +1.5],   # out_0 = through, out_1 = cross
    output_monitor_widths=[3.0, 3.0],
    decay_tol=1e-11,
)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

_results = []


def _is_master():
    """True on rank 0 (or when meep isn't loaded). All 48 ranks run the same
    checks, but only one may print or write — otherwise the log is 48-way
    interleaved and the output files have 48 concurrent writers."""
    try:
        import meep as mp
        return bool(mp.am_really_master())
    except Exception:
        return True


MASTER = _is_master()


def say(*args):
    if MASTER:
        print(*args, flush=True)


def record(name, passed, detail, **numbers):
    """Log one check and keep it for verification.json."""
    _results.append(dict(check=name, passed=bool(passed), detail=detail, **numbers))
    say(f"  [{'PASS' if passed else 'FAIL'}] {name}: {detail}")
    return passed


def provenance(static_only):
    """Everything needed to reproduce or date this run. Stored and printed —
    the archived reference's missing provenance is exactly why its validity had
    to be argued from circumstantial evidence."""
    try:
        sha = subprocess.check_output(
            ["git", "-C", os.path.join(_ROOT, "BEEF"), "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
        dirty = subprocess.check_output(
            ["git", "-C", os.path.join(_ROOT, "BEEF"), "status", "--porcelain"],
            text=True, stderr=subprocess.DEVNULL).strip()
        beef_git = f"{sha}{'+dirty' if dirty else ''}"
    except Exception:
        beef_git = "<not a git checkout>"

    try:
        import meep as mp
        meep_version, nproc = mp.__version__, int(mp.count_processors())
    except Exception:
        meep_version, nproc = "<not imported>", 1

    return dict(
        timestamp_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        static_only=static_only,
        beef_version=getattr(beef, "__version__", "<none>"),
        beef_git=beef_git,
        meep_version=meep_version,
        python=platform.python_version(),
        host=platform.node(),
        mpi_ranks=nproc,
        slurm_job_id=os.environ.get("SLURM_JOB_ID", "<none>"),
        cached=False,
        device="DirectionalCoupler",
        device_kwargs=DEVICE_KWARGS,
        solver="FDTDSolver",
        solver_kwargs={**SOLVER_KWARGS, "material_stack": "SOI_220",
                       "source_bandwidth_factor": FDTDSolver().source_bandwidth_factor},
        wavelengths=list(map(float, wavelengths)),
    )


def _print_provenance(p):
    say("\nProvenance")
    for k in ("timestamp_utc", "beef_version", "beef_git", "meep_version", "python",
              "host", "mpi_ranks", "slurm_job_id", "cached"):
        say(f"  {k:16s} {p[k]}")
    say(f"  device           DirectionalCoupler({_fmt_kwargs(DEVICE_KWARGS)})")
    say(f"  solver           FDTDSolver({_fmt_kwargs(p['solver_kwargs'])})")
    say(f"  wavelengths      {WL_MIN}–{WL_MAX} µm, {N_WL} pts")


def _fmt_kwargs(d):
    return ", ".join(f"{k}={v!r}" for k, v in d.items())


# ---------------------------------------------------------------------------
# Static checks — the frequency grid and the source that has to cover it
# ---------------------------------------------------------------------------
def static_checks():
    say("\nStatic checks (no FDTD run)")
    mon_freqs, wl_order = monitor_freq_grid(wavelengths)
    got = 1.0 / mon_freqs

    record("1. grid ascending in frequency",
           bool(np.all(np.diff(mon_freqs) > 0)),
           f"{len(mon_freqs)} monitor frequencies, "
           f"{mon_freqs.min():.5f}-{mon_freqs.max():.5f} 1/um, strictly increasing",
           n_freqs=len(mon_freqs), f_min=float(mon_freqs.min()), f_max=float(mon_freqs.max()))

    err_nm = np.abs(np.sort(got) - np.sort(wavelengths)).max() * 1000
    record("2. grid matches requested wavelengths",
           err_nm < 1e-6,
           f"max deviation {err_nm:.2e} nm (the pre-fix grid was off by up to 75 nm)",
           max_error_nm=float(err_nm), prefix_max_error_nm=75.4)

    record("3. wl_order round-trips to caller order",
           bool(np.allclose(got[wl_order], wavelengths)),
           "monitor-ordered array reorders back onto the requested wavelengths",
           wl_order=list(map(int, wl_order)))

    fwidth = source_fwidth(wavelengths, FDTDSolver().source_bandwidth_factor)
    fcen = 1.0 / float(np.mean(wavelengths))
    frac_new = check_source_coverage(fcen, fwidth, wavelengths)
    record("4. source covers the whole band",
           float(frac_new.min()) >= 0.10,
           f"weakest requested wavelength carries {frac_new.min():.3f} of peak "
           f"amplitude (floor 0.10); fcen={fcen:.5f}, fwidth={fwidth:.5f}",
           min_fraction=float(frac_new.min()), fcen=float(fcen), fwidth=float(fwidth))

    span = 1.0 / wavelengths.min() - 1.0 / wavelengths.max()
    # min_fraction=0.0 suppresses the warning we are deliberately provoking here.
    frac_old = check_source_coverage(fcen, 0.5 * span, wavelengths, min_fraction=0.0)
    record("5. pre-fix bandwidth really was starved",
           float(frac_old.min()) < 1e-6,
           f"old span*0.5 (fwidth={0.5 * span:.5f}) left the band edges at "
           f"{frac_old.min():.2e} of peak",
           old_min_fraction=float(frac_old.min()), old_fwidth=float(0.5 * span))

    single = source_fwidth(np.array([1.55]), 2.0)
    record("6. single wavelength gets a usable pulse",
           single > 0,
           f"fwidth = {single:.5f} (meep's pulse width is 1/fwidth, so 0 is unusable)",
           single_wl_fwidth=float(single))

    say("\n  monitor grid (λ requested → λ meep uses → pre-fix λ meep used):")
    old = _old_grid()
    for i in range(N_WL):
        say(f"    {wavelengths[i]:.4f} → {got[wl_order][i]:.4f} → {old[i]:.4f} µm "
            f"({(old[i] - wavelengths[i]) * 1000:+.1f} nm pre-fix error)")

    if MASTER:
        _plot_freq_grid(got, wl_order)
        _plot_source_coverage(frac_new, frac_old)
    return mon_freqs, wl_order, got, frac_new, frac_old


def _old_grid():
    """The wavelengths the pre-fix (fcen, df*0.5, nfreq) triple actually sampled."""
    fcen = 1.0 / float(np.mean(wavelengths))
    df = (1.0 / wavelengths.min() - 1.0 / wavelengths.max()) * 0.5
    return 1.0 / np.linspace(fcen - df / 2, fcen + df / 2, N_WL)


def _plot_freq_grid(got, wl_order):
    old = _old_grid()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(wavelengths, old, "o--", color="crimson", lw=1.2, ms=5,
            label="Pre-fix: what meep actually sampled")
    ax.plot(wavelengths, got[wl_order], "o-", color="seagreen", lw=1.4, ms=5,
            label="Post-fix: what meep actually samples")
    ax.plot(wavelengths, wavelengths, ":", color="gray", lw=1.0, label="Ideal (y = x)")
    ax.set_xlabel("Wavelength the code labels the result with (µm)")
    ax.set_ylabel("Wavelength meep computed it at (µm)")
    ax.set_title("Monitor Frequency Grid — Item 1")
    ax.legend(fontsize=8)
    ax.grid(True)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "freq_grid.png"), dpi=300)
    plt.close("all")


def _plot_source_coverage(frac_new, frac_old):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.semilogy(wavelengths, np.clip(frac_old, 1e-16, None), "o--", color="crimson",
                lw=1.2, ms=5, label="Pre-fix (span × 0.5)")
    ax.semilogy(wavelengths, np.clip(frac_new, 1e-16, None), "o-", color="seagreen",
                lw=1.4, ms=5, label="Post-fix (span × 2.0)")
    ax.axhline(0.1, color="gray", ls=":", lw=0.9, label="10% floor")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Source amplitude / peak")
    ax.set_title("Gaussian Source Coverage — Item 2")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "source_coverage.png"), dpi=300)
    plt.close("all")


# ---------------------------------------------------------------------------
# Physical check — one 3-D FDTD solve against the archived pre-fix run
# ---------------------------------------------------------------------------
def physical_checks():
    say(f"\nPhysical check — 3-D FDTD, coupling_length {COUPLING_LENGTH:.1f} µm, "
        f"res {RESOLUTION:.0f}, uncached")

    dc = DirectionalCoupler(**DEVICE_KWARGS)
    solver = FDTDSolver(material_stack=SOI_220, **SOLVER_KWARGS)

    t0 = time.time()
    smat = solver.solve(dc, wavelengths, port="bot", input_mode=1)  # uncached, see above
    elapsed = time.time() - t0
    say(f"  FDTD solve finished in {elapsed / 60:.1f} min")

    T_through = smat.transmission(port="out_0", in_mode=1, out_mode=1)
    T_cross = smat.transmission(port="out_1", in_mode=1, out_mode=1)
    T_total = T_through + T_cross

    say("\n  λ (µm)   T_through   T_cross    T_total")
    for i in range(N_WL):
        say(f"   {wavelengths[i]:.3f}    {T_through[i]:.4f}     {T_cross[i]:.4f}    {T_total[i]:.4f}")

    # Check 7: coupling must strengthen with wavelength.
    slope = float(np.polyfit(wavelengths, T_cross, 1)[0])
    record("7. cross transmission rises with wavelength",
           slope > 0,
           f"dT_cross/dλ = {slope:+.3f} /µm "
           f"({T_cross[0]:.4f} at {WL_MIN:.2f} µm → {T_cross[-1]:.4f} at {WL_MAX:.2f} µm); "
           f"the pre-fix run had this backwards",
           slope=slope, T_cross_first=float(T_cross[0]), T_cross_last=float(T_cross[-1]))

    # Check 8: the archived pre-fix row, on its true axis, must land on this curve.
    ref_true_wl = ref_T_cross = residual = None
    if not os.path.exists(_REF_NPZ):
        record("8. matches pre-fix data on its true axis", False,
               f"SKIPPED — archived reference not found at {_REF_NPZ}")
    else:
        ref = np.load(_REF_NPZ)
        row = int(np.argmin(np.abs(ref["resolutions"] - _REF_RES)))
        ref_T_cross = ref["T_cross"][row]
        ref_true_wl = _old_grid()   # what that run actually sampled

        # Tolerance 0.05 is deliberately loose — the claim under test is the
        # *shape* (old data merely mislabelled), not a golden value, and the
        # two runs sample different bands so some interpolation error is
        # expected. For calibration: a crude straight line between the reversed
        # archive's own endpoints already scores 0.043, so a bare pass near the
        # tolerance means little. Read rms_deviation and the residual panel,
        # not just the PASS — if the hypothesis holds this should come in an
        # order of magnitude tighter.
        order = np.argsort(ref_true_wl)
        interp = np.interp(ref_true_wl[order], wavelengths, T_cross)
        residual = interp - ref_T_cross[order]
        dev = float(np.abs(residual).max())
        say(f"\n  check 8 residuals (post-fix interpolated − pre-fix on true axis), "
            f"res={ref['resolutions'][row]:.0f}:")
        for w, r in zip(ref_true_wl[order], residual):
            say(f"    {w:.4f} µm   {r:+.4f}")
        record("8. matches pre-fix data on its true axis",
               dev < 0.05,
               f"max |ΔT_cross| = {dev:.4f} over "
               f"{ref_true_wl.min():.4f}–{ref_true_wl.max():.4f} µm (tolerance 0.05) "
               f"— the old run was only mislabelled",
               max_deviation=dev, rms_deviation=float(np.sqrt((residual ** 2).mean())))

    # Check 9: the run has to be physically sane for 7-8 to mean anything.
    worst = float(np.abs(T_total - 1.0).max())
    record("9. power conserved",
           worst < 0.02,
           f"through + cross = {T_total.min():.4f}–{T_total.max():.4f}; "
           f"worst deviation from unity {worst:.4f} (tolerance 0.02)",
           T_total_min=float(T_total.min()), T_total_max=float(T_total.max()),
           worst_deviation=worst)
    try:
        smat.check_passive()
        say("  SMatrix.check_passive() OK")
    except Exception as exc:
        say(f"  SMatrix.check_passive() raised: {exc}")

    if MASTER:
        _plot_dispersion(T_through, T_cross, ref_true_wl, ref_T_cross, residual)
        _plot_power(T_through, T_cross, T_total)
    return T_through, T_cross, T_total, ref_true_wl, ref_T_cross, residual, elapsed


def _plot_dispersion(T_through, T_cross, ref_true_wl, ref_T_cross, residual):
    has_ref = ref_T_cross is not None
    fig, axes = plt.subplots(1, 2 if has_ref else 1, figsize=(13 if has_ref else 8, 5),
                             squeeze=False)
    ax = axes[0][0]
    if has_ref:
        ax.plot(wavelengths, ref_T_cross, "s--", color="crimson", lw=1.2, ms=5,
                label="Pre-fix, as it was labelled (backwards)")
        ax.plot(ref_true_wl, ref_T_cross, "s:", color="darkorange", lw=1.4, ms=5,
                label="Pre-fix, on the wavelengths meep really sampled")
    ax.plot(wavelengths, T_cross, "o-", color="seagreen", lw=1.6, ms=5,
            label="Post-fix run (cross)")
    ax.plot(wavelengths, T_through, "o-", color="steelblue", lw=1.2, ms=4, alpha=0.7,
            label="Post-fix run (through)")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Transmission")
    ax.set_title(f"Coupling vs Wavelength, Before and After "
                 f"(L$_c$ = {COUPLING_LENGTH:.1f} µm)")
    ax.legend(fontsize=8)
    ax.grid(True)

    if has_ref:
        ax = axes[0][1]
        order = np.argsort(ref_true_wl)
        ax.axhline(0, color="gray", lw=0.8)
        ax.axhspan(-0.05, 0.05, color="seagreen", alpha=0.10, label="tolerance ±0.05")
        ax.plot(ref_true_wl[order], residual, "o-", color="crimson", lw=1.3, ms=5)
        ax.set_xlabel("Wavelength meep really sampled (µm)")
        ax.set_ylabel("Post-fix − pre-fix  ΔT$_{cross}$")
        ax.set_title("Check 8 residual — is the old data merely mislabelled?")
        ax.legend(fontsize=8)
        ax.grid(True)

    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "dispersion_before_after.png"), dpi=300)
    plt.close("all")


def _plot_power(T_through, T_cross, T_total):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(wavelengths, T_through, "o-", color="steelblue", lw=1.3, ms=4, label="Through")
    ax.plot(wavelengths, T_cross, "o--", color="darkorange", lw=1.3, ms=4, label="Cross")
    ax.plot(wavelengths, T_total, "o-", color="k", lw=1.5, ms=4, label="Total")
    ax.axhline(1.0, color="gray", ls=":", lw=0.9)
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Transmission")
    ax.set_title("Power Conservation — check 9")
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=8)
    ax.grid(True)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "power_conservation.png"), dpi=300)
    plt.close("all")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--static-only", action="store_true",
                    help="skip the FDTD solve — checks 1-6 only, safe on a login node")
    args = ap.parse_args()

    if MASTER:
        os.makedirs(output_dir, exist_ok=True)

    say("Verifying dc_validation.md items 1 (monitor frequency grid) "
        "and 2 (source bandwidth)")
    prov = provenance(args.static_only)
    _print_provenance(prov)

    _, wl_order, got, frac_new, frac_old = static_checks()

    payload = dict(wavelengths=wavelengths, monitor_wavelengths=got[wl_order],
                   old_grid=_old_grid(), source_frac=frac_new, source_frac_prefix=frac_old)
    if args.static_only:
        say("\n--static-only: skipping the FDTD solve (checks 7-9 not run).")
    else:
        (T_through, T_cross, T_total, ref_true_wl,
         ref_T_cross, residual, elapsed) = physical_checks()
        prov["fdtd_solve_seconds"] = round(elapsed, 1)
        payload.update(T_through=T_through, T_cross=T_cross, T_total=T_total)
        if ref_T_cross is not None:
            payload.update(ref_true_wl=ref_true_wl, ref_T_cross=ref_T_cross,
                           check8_residual=residual)

    n_pass = sum(r["passed"] for r in _results)
    if MASTER:
        np.savez(os.path.join(output_dir, "verification.npz"), **payload)
        with open(os.path.join(output_dir, "verification.json"), "w") as f:
            json.dump(dict(provenance=prov,
                           summary=dict(passed=n_pass, total=len(_results),
                                        all_passed=n_pass == len(_results)),
                           checks=_results), f, indent=2)

    say(f"\n{n_pass}/{len(_results)} checks passed → {output_dir}")
    if n_pass != len(_results):
        say("FAILED — see verification.json")
        for r in _results:
            if not r["passed"]:
                say(f"  FAIL {r['check']}: {r['detail']}")
        sys.exit(1)
    say("All checks passed.")


if __name__ == "__main__":
    main()
