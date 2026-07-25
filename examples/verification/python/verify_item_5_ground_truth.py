"""TEMPORARY — ground-truth FDTD length sweep for dc_validation.md item 5.

Delete this folder (and the Makefile's verification targets, and
batch/verification/) once the directional coupler work is signed off.

Everything else in dc_validation.md is judged against a full-device FDTD length
sweep, and there wasn't a usable one. Three problems with the reference as it
stood:

  1. example_output/directional_coupler/directional_coupler_fdtd_sweep/
     length_sweep.npz does not exist in the live tree at all -- the stub sweep's
     reference line reads from a missing file.
  2. The archived stand-in (archive/example_output_20260717/...) was written
     pre-fix. It is recoverable -- re-read on the wavelengths meep actually
     sampled it gives a clean, physically correct L50(lambda) -- but that axis
     spans only 1.5254-1.5754 um, half the nominal band:

         true lambda (um)   1.5754   1.5500   1.5254
         L50 (um)            13.17    15.30    17.65

  3. It stores bare arrays: no config, no solver fingerprint, no provenance. Its
     validity had to be argued from circumstantial evidence rather than read off
     the file. That is the complaint this script exists to stop repeating.

This is a *tutorial-independent* re-run. examples/directional_coupler/python/
directional_coupler_fdtd_sweep.py keeps its deliberate 0:2:60 grid and its
teaching role; this one is scoped to measuring L50(lambda) properly.

Grid (40 lengths, one per SLURM array task):

    COARSE  0:4:60     shape and context -- 4 um resolves the oscillation fine
    FINE    9:0.5:22   the crossing, 4x finer than the archived 2 um

The fine window comes from the table above: dL50/dlambda ~ -90 um/um, so across
the full 1.50-1.60 um band the crossing runs from ~19.8 um down to ~10.8 um.
9-22 um covers that with margin.

Device and solver settings are pinned to directional_coupler_fdtd_sweep.py's
make_coupler()/make_solver() exactly, so this sweep is directly comparable both
to the archived one and to job 11477911 (the items 1-2 verification, which
solved the same device at coupling_length 13.5 um).

Outputs (example_output/verification/ground_truth_sweep/):

    data/L{tag}.npz          per-length raw data, one per array task
    length_sweep.npz         aggregate grid + L50(lambda) + embedded provenance
    provenance.json          the same provenance, standalone and greppable
    length_sweep.png         through/cross vs length, shaded by wavelength spread
    l50_vs_wavelength.png    the headline: L50(lambda), with the archived run's
                             true-axis points overlaid for continuity

Run:
    # one length per array task (see batch/verification/):
    sbatch batch/verification/run_verify_item_5_ground_truth.batch

    # a single length by index, or the whole sweep serially (slow):
    python .../verify_item_5_ground_truth.py <index 0..39>

    # re-aggregate whatever has finished, no simulation, safe on a login node:
    OMP_NUM_THREADS=1 python .../verify_item_5_ground_truth.py --aggregate-only

Every invocation re-aggregates the lengths finished so far and writes the
aggregate files atomically, so the plot is usable ("n/40 lengths") while the
array is draining and the last task to finish leaves the complete sweep behind.
"""

import argparse
import json
import logging
import os
import platform
import socket
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

# beef warns through logging -- check_source_coverage() in particular. Without
# this a starved-source warning would never reach the batch log.
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "ground_truth_sweep")
data_dir = os.path.join(output_dir, "data")

# Shared with every other example. Caching is ON here, unlike verify_items_1_2.py
# which was deliberately uncached because a verification must not be able to pass
# on a stored result. This script *produces* reference data at ~1-3 h per task, so
# a re-run after a node failure should be free. That is only safe now that item 3
# completed the fingerprint (a stale FDTD entry can no longer be hit: the resolved
# monitor grid, source bandwidth and subpixel settings are all hashed) and item 4
# made the writes rank-safe. Provenance records which mode was used. --no-cache
# forces a from-scratch re-derivation.
_CACHE_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "example_sim_cache"))

# The archived pre-fix sweep, overlaid on l50_vs_wavelength.png after being put
# back on the wavelengths meep really sampled. Continuity, not validation.
_ARCHIVED_NPZ = os.path.join(
    _ROOT, "archive", "example_output_20260717", "directional_coupler",
    "directional_coupler_fdtd_sweep", "length_sweep.npz")

# Pinned to directional_coupler_fdtd_sweep.py -- do not "improve" these.
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11

COARSE = np.arange(0.0, 60.0 + 1e-9, 4.0)     # 16 pts
FINE = np.arange(9.0, 22.0 + 1e-9, 0.5)       # 27 pts
COUPLING_LENGTHS = np.unique(np.concatenate([COARSE, FINE]))   # 40 pts, 3 shared

DEVICE_KWARGS = dict(
    waveguide_width=0.45, gap=0.25,
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


# ---------------------------------------------------------------------------
# MPI rank guarding. All 48 ranks run the whole script; only one may write.
# ---------------------------------------------------------------------------
def _is_master():
    """True on rank 0, or when meep isn't loaded (serial --aggregate-only)."""
    try:
        import meep as mp
        return bool(mp.am_really_master())
    except Exception:
        return True


MASTER = _is_master()


def say(*args):
    if MASTER:
        print(*args, flush=True)


def _atomic_replace_write(path, write_fn):
    """write_fn(tmp) then atomically rename onto path.

    Host-qualified temp name: PIDs are allocated per host, so a PID-only name is
    not unique on a filesystem shared between nodes (dc_validation.md item 4).
    Callers are already rank-gated, so this is belt-and-braces -- but 40 array
    tasks all re-aggregating the same files is exactly the setup that produced
    the item 4 crash, and the tutorial sweep scripts still have the PID-only
    form.
    """
    tmp = f"{path}.tmp.{socket.gethostname()}.{os.getpid()}"
    write_fn(tmp)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
def _tag(L):
    """Stable per-length file tag, e.g. L09.5um / L60.0um."""
    return f"L{L:04.1f}um"


def make_coupler(L):
    return DirectionalCoupler(coupling_length=float(L), **DEVICE_KWARGS)


def make_solver():
    return FDTDSolver(material_stack=SOI_220, **SOLVER_KWARGS)


def provenance(cached):
    """Everything needed to reproduce or date this sweep. The archived
    reference's missing provenance is precisely why its validity had to be
    argued instead of read off -- see this module's docstring."""
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

    solver = make_solver()
    return dict(
        timestamp_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        beef_version=getattr(beef, "__version__", "<none>"),
        beef_git=beef_git,
        meep_version=meep_version,
        python=platform.python_version(),
        host=platform.node(),
        mpi_ranks=nproc,
        slurm_job_id=os.environ.get("SLURM_JOB_ID", "<none>"),
        slurm_array_job_id=os.environ.get("SLURM_ARRAY_JOB_ID", "<none>"),
        cached=bool(cached),
        device="DirectionalCoupler",
        device_kwargs=DEVICE_KWARGS,
        solver="FDTDSolver",
        solver_kwargs={
            **SOLVER_KWARGS,
            "material_stack": "SOI_220",
            "source_bandwidth_factor": solver.source_bandwidth_factor,
            "eps_averaging": solver.eps_averaging,
            "subpixel_tol": solver.subpixel_tol,
            "subpixel_maxeval": solver.subpixel_maxeval,
        },
        wavelengths=list(map(float, wavelengths)),
        coupling_lengths=list(map(float, COUPLING_LENGTHS)),
        grid="union(arange(0,60,4), arange(9,22,0.5))",
    )


def _print_provenance(p):
    say("\nProvenance")
    for k in ("timestamp_utc", "beef_version", "beef_git", "meep_version", "python",
              "host", "mpi_ranks", "slurm_job_id", "cached"):
        say(f"  {k:18s} {p[k]}")
    say(f"  device             DirectionalCoupler({', '.join(f'{k}={v!r}' for k, v in DEVICE_KWARGS.items())})")
    say(f"  solver             FDTDSolver({', '.join(f'{k}={v!r}' for k, v in p['solver_kwargs'].items())})")
    say(f"  wavelengths        {WL_MIN}-{WL_MAX} um, {N_WL} pts")
    say(f"  lengths            {len(COUPLING_LENGTHS)} pts, "
        f"{COUPLING_LENGTHS.min():.1f}-{COUPLING_LENGTHS.max():.1f} um "
        f"({p['grid']})")


# ---------------------------------------------------------------------------
# Per-length simulation
# ---------------------------------------------------------------------------
def simulate_length(L, cached=True):
    dc = make_coupler(L)
    solver = make_solver()
    say(f"\n[gt] coupling_length = {L:.1f} um "
        f"(device total {dc.total_length:.1f} um, one 3-D FDTD run, "
        f"{N_WL} wavelengths {WL_MIN:.2f}-{WL_MAX:.2f} um)")

    t0 = time.time()
    smat = solver.solve(dc, wavelengths, port="bot", input_mode=1,
                        save_dir=_CACHE_DIR if cached else None)
    elapsed = time.time() - t0
    smat.check_passive()   # total output power must not exceed 1 anywhere

    T_through = smat.transmission(port="out_0", in_mode=1, out_mode=1)
    T_cross = smat.transmission(port="out_1", in_mode=1, out_mode=1)
    T_total = T_through + T_cross

    if MASTER:
        os.makedirs(data_dir, exist_ok=True)
        _atomic_replace_write(
            os.path.join(data_dir, f"{_tag(L)}.npz"),
            lambda tmp: np.savez(
                open(tmp, "wb"), coupling_length=L, wavelengths=wavelengths,
                T_through=T_through, T_cross=T_cross,
                solve_seconds=elapsed, cached=cached),
        )

    say(f"[gt] L={L:.1f} um solved in {elapsed / 60:.1f} min"
        f"{' (cache hit)' if elapsed < 5 else ''}")
    say(f"[gt]   lambda    through    cross     total")
    for j in range(N_WL):
        say(f"[gt]   {wavelengths[j]:.3f}    {T_through[j]:.4f}   "
            f"{T_cross[j]:.4f}   {T_total[j]:.4f}")
    worst = float(np.abs(T_total - 1.0).max())
    say(f"[gt]   worst deviation from unity: {worst:.4f}"
        + ("  <-- CHECK, expected <0.01" if worst > 0.01 else ""))


# ---------------------------------------------------------------------------
# L50 extraction
# ---------------------------------------------------------------------------
def l50_from_grid(lengths, T_cross_grid):
    """First upward crossing of T_cross = 0.5 per wavelength, linearly
    interpolated, plus the bracketing lengths.

    The bracket is returned, not just the interpolated value, so the
    interpolation width is auditable rather than implied -- 15.30 um quoted from
    a 2 um bracket and from a 0.5 um bracket are not the same claim, and the
    archived file gave no way to tell which it was.
    """
    n_wl = T_cross_grid.shape[1]
    l50 = np.full(n_wl, np.nan)
    lo = np.full(n_wl, np.nan)
    hi = np.full(n_wl, np.nan)
    for j in range(n_wl):
        col = T_cross_grid[:, j]
        for i in range(len(lengths) - 1):
            if col[i] < 0.5 <= col[i + 1]:
                span = col[i + 1] - col[i]
                l50[j] = lengths[i] + (0.5 - col[i]) * (lengths[i + 1] - lengths[i]) / span
                lo[j], hi[j] = lengths[i], lengths[i + 1]
                break
    return l50, lo, hi


def _archived_l50():
    """The archived pre-fix sweep's L50, on the wavelengths meep really sampled.

    Reconstructs the grid the pre-fix code produced: fwidth = span * 0.5 handed
    to meep as df, which spreads N points evenly over [fcen - df/2, fcen + df/2]
    ascending in frequency, i.e. descending in wavelength. Job 11477911's check 8
    validated this exact relabelling to 2.5e-4 against a fresh post-fix solve.
    Returns (true_wavelengths, L50) or None.
    """
    if not os.path.exists(_ARCHIVED_NPZ):
        return None
    d = np.load(_ARCHIVED_NPZ)
    wl, L, Tc = d["wavelengths"], d["coupling_lengths"], d["T_cross"]
    fcen = 1.0 / ((wl.min() + wl.max()) / 2)
    df = (1.0 / wl.min() - 1.0 / wl.max()) * 0.5
    true_wl = 1.0 / np.linspace(fcen - df / 2, fcen + df / 2, len(wl))
    l50, _, _ = l50_from_grid(L, Tc)
    return true_wl, l50


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
def aggregate(cached=True):
    done, through_rows, cross_rows = [], [], []
    for L in COUPLING_LENGTHS:
        path = os.path.join(data_dir, f"{_tag(L)}.npz")
        if not os.path.exists(path):
            continue
        d = np.load(path)
        done.append(L)
        through_rows.append(d["T_through"])
        cross_rows.append(d["T_cross"])

    n_done, n_all = len(done), len(COUPLING_LENGTHS)
    if n_done < 2:
        say(f"[gt] Aggregation skipped -- only {n_done}/{n_all} lengths done.")
        return
    if not MASTER:
        return

    done = np.asarray(done)
    T_through = np.vstack(through_rows)
    T_cross = np.vstack(cross_rows)
    T_total = T_through + T_cross
    l50, l50_lo, l50_hi = l50_from_grid(done, T_cross)

    prov = provenance(cached)
    prov["lengths_completed"] = int(n_done)
    prov_json = json.dumps(prov, indent=2, sort_keys=True)

    os.makedirs(output_dir, exist_ok=True)
    _atomic_replace_write(
        os.path.join(output_dir, "length_sweep.npz"),
        lambda tmp: np.savez(
            open(tmp, "wb"),
            coupling_lengths=done, wavelengths=wavelengths,
            T_through=T_through, T_cross=T_cross,
            L50=l50, L50_bracket_lo=l50_lo, L50_bracket_hi=l50_hi,
            provenance=np.array(prov_json)),
    )
    _atomic_replace_write(
        os.path.join(output_dir, "provenance.json"),
        lambda tmp: open(tmp, "w").write(prov_json + "\n"),
    )
    say(f"\n[gt] Aggregate saved ({n_done}/{n_all} lengths) -> {output_dir}/length_sweep.npz")

    # -- L50(lambda), the headline ------------------------------------------
    say("\n[gt] L50 vs wavelength")
    say("[gt]   lambda (um)   L50 (um)   bracket")
    for j in range(N_WL):
        if np.isnan(l50[j]):
            say(f"[gt]   {wavelengths[j]:.3f}         --       (no crossing in range)")
        else:
            say(f"[gt]   {wavelengths[j]:.3f}      {l50[j]:6.2f}     "
                f"[{l50_lo[j]:.1f}, {l50_hi[j]:.1f}]")

    valid = ~np.isnan(l50)
    if valid.sum() >= 2:
        slope = np.polyfit(wavelengths[valid], l50[valid], 1)[0]
        say(f"[gt]   dL50/dlambda = {slope:+.1f} um/um"
            + ("  OK (coupling strengthens with wavelength)" if slope < 0
               else "  <-- WRONG SIGN, the axis is flipped again"))

    _plot_length_sweep(done, T_through, T_cross, T_total, n_done, n_all)
    _plot_l50(l50, l50_lo, l50_hi)
    say(f"[gt] Plots saved -> {output_dir}/")


def _plot_length_sweep(lengths, T_through, T_cross, T_total, n_done, n_all):
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for grid, label, color, ls in [
        (T_through, "Through", "steelblue", "-"),
        (T_cross, "Cross", "darkorange", "--"),
        (T_total, "Total", "k", "-"),
    ]:
        ax.plot(lengths, grid.mean(axis=1), label=f"{label} (mean over $\\lambda$)",
                color=color, marker="o", markersize=2.5, linestyle=ls,
                lw=1.5 if label == "Total" else 1.2)
        ax.fill_between(lengths, grid.min(axis=1), grid.max(axis=1),
                        color=color, alpha=0.15, lw=0)
    ax.axhline(0.5, color="gray", linestyle=":", lw=0.8, label="50:50")
    if n_done < n_all:
        ax.text(0.02, 0.97, f"{n_done}/{n_all} lengths done", transform=ax.transAxes,
                ha="left", va="top", fontsize=8, color="#555555",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          alpha=0.7, edgecolor="none"))
    ax.set_xlabel("Coupling length (µm)")
    ax.set_ylabel("Transmission (shaded = min–max over $\\lambda$)")
    ax.set_title("Directional Coupler — Ground-Truth FDTD Length Sweep")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    _atomic_replace_write(os.path.join(output_dir, "length_sweep.png"),
                          lambda tmp: fig.savefig(tmp, dpi=300, format="png"))
    plt.close("all")


def _plot_l50(l50, l50_lo, l50_hi):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    valid = ~np.isnan(l50)
    ax.fill_between(wavelengths[valid], l50_lo[valid], l50_hi[valid],
                    color="steelblue", alpha=0.18, lw=0,
                    label="interpolation bracket")
    ax.plot(wavelengths[valid], l50[valid], color="steelblue", marker="o",
            markersize=4, lw=1.6, label="This run (post-fix, full band)")

    arch = _archived_l50()
    if arch is not None:
        true_wl, arch_l50 = arch
        ok = ~np.isnan(arch_l50)
        ax.plot(true_wl[ok], arch_l50[ok], color="darkorange", marker="s",
                markersize=4, lw=0, alpha=0.85,
                label="Archived pre-fix, on its true axis")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("50:50 coupling length $L_{50}$ (µm)")
    ax.set_title("Directional Coupler — 50:50 Length vs Wavelength (ground truth)")
    ax.legend(fontsize=8)
    ax.grid(True)
    plt.tight_layout()
    _atomic_replace_write(os.path.join(output_dir, "l50_vs_wavelength.png"),
                          lambda tmp: fig.savefig(tmp, dpi=300, format="png"))
    plt.close("all")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("index", nargs="?", type=int, default=None,
                        help=f"length index 0..{len(COUPLING_LENGTHS) - 1}; "
                             "omit to run every length serially")
    parser.add_argument("--aggregate-only", action="store_true",
                        help="re-aggregate finished lengths, simulate nothing")
    parser.add_argument("--no-cache", action="store_true",
                        help="ignore and do not write the simulation cache")
    args = parser.parse_args()
    cached = not args.no_cache

    if MASTER:
        os.makedirs(data_dir, exist_ok=True)

    if args.aggregate_only:
        say("[gt] Aggregate-only: no simulation.")
        aggregate(cached)
        sys.exit(0)

    _print_provenance(provenance(cached))

    if args.index is not None:
        if not 0 <= args.index < len(COUPLING_LENGTHS):
            raise SystemExit(
                f"index {args.index} out of range -- sweep has "
                f"{len(COUPLING_LENGTHS)} lengths (0..{len(COUPLING_LENGTHS) - 1})")
        simulate_length(COUPLING_LENGTHS[args.index], cached)
    else:
        say(f"[gt] No index given -- running all {len(COUPLING_LENGTHS)} lengths "
            "serially (the SLURM array batch script is much faster).")
        for L in COUPLING_LENGTHS:
            simulate_length(L, cached)

    aggregate(cached)
