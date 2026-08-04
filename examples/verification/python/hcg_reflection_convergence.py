"""Resolution and PML/domain-size (x, y, z) convergence for the HCG
reflector's FDTD reflection spectrum --
validation_markdown/high_contrast_grating_validation.md step 8, items 1-4
(3-D, resolution convergence, decay_tol behavior, PML/domain convergence).

Delete this script (and its Makefile target / batch script) once step 8
items 1-4 close and a production resolution/domain is settled in
examples/high_contrast_grating/python/high_contrast_grating_fdtd.py.

Every point is a real 3-D ("Working in 3D dimensions." in every point's own
meep log -- the real 160 nm Si core, z-face at +-0.08 um = +-wg_height/2, PML
on all six faces) FDTDSolver.solve_reflection() call, single wavelength
(1.569 um -- the exact peak-reflectance wavelength the design paper reports,
not a round band-centre guess) to keep each point's cost down; the tutorial
script itself already demonstrates the full-band sweep is "free" once a
resolution/domain is settled (one broadband pulse, many DFT frequencies).

WAVELENGTH was 1.31 um (O-band) through this sweep's first round (indices
0-8, plus the tutorial script's own WL_MIN/WL_MAX) -- WRONG. The paper this
device is from (Hiesener, Pesch, Ralph, IEEE Photonics Technology Letters
36(9), 597-600, 2024) reports "a peak reflectance of 98.9% simulated centered
at 1569 nm with a minimum reflectance of 97.9% over C-band" -- this is a
C-band mirror, not an O-band one. Caught 2026-08-04 (user flagged it directly)
after the O-band sweep was already most of the way done; every O-band number
this script or the tutorial ever produced is superseded, not merely
re-labelled -- R at the wrong wavelength for a resonant, wavelength-selective
structure like this one is not informative about the C-band device at all.

decay_tol is 1e-11 (the tutorial's production value) at every point, not
loosened for the sweep -- so whether that tolerance is reachable in a
reasonable wall time at each resolution IS this sweep's own timing data
(step 8 item 3), not a separate check.

Grid (15 points, one per SLURM array task) -- see POINTS/SERIES below for the
authoritative values; summary:

    index 0-3, 6-8   resolution series      res in {10,15,20,25,30,35,40},
                                             sim_width=7.5, edge_clearance=0.3,
                                             sim_height=1.5
    index 4-5        y (sim_width) @res=10  sim_width in {8.5, 9.5}
    index 9-10       x (edge_clearance)     edge_clearance in {0.5, 0.8}
                     @res=25 (baseline: index 3)
    index 11-12      z (sim_height)         sim_height in {2.0, 2.5}
                     @res=25 (baseline: index 3)
    index 13-14      y (sim_width) @res=25  sim_width in {8.5, 9.5}
                     (redo of 4-5 at a resolution closer to converged, once
                     10-40 showed res=10 was nowhere near it)

Indices 9-14 added 2026-08-04, after the resolution series showed res 10-25
was still climbing (R = 0.0237/0.1784/0.3063/0.3784, |dR| shrinking but not
plateaued) and the user asked for x/y/z domain convergence specifically, not
just the original y-only check. index 3 (res=25, the highest point that had
landed then) is the shared baseline for all three res=25 domain series --
same reuse pattern as index 0 for the original y-at-res=10 series.

Run:
    # one point per array task (see batch/verification/):
    sbatch batch/verification/run_hcg_reflection_convergence.batch

    # a single point by index, or the whole grid serially (slow):
    python .../hcg_reflection_convergence.py <index 0..14>

    # re-aggregate whatever has finished, no simulation, safe on a login node:
    OMP_NUM_THREADS=1 python .../hcg_reflection_convergence.py --aggregate-only

Outputs (example_output/verification/hcg_reflection_convergence/):

    data/idx{N}.npz              per-point raw data, one per array task
    convergence.npz              aggregate + embedded provenance
    provenance.json              the same provenance, standalone and greppable
    {series}_convergence.png     one plot per SERIES entry below
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
from beef import GdsGeometry, FDTDSolver, SOI_160, load_gds_polygons

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
_GDS_PATH = os.path.join(_ROOT, "BEEF", "examples", "gds_files", "HCG_seeded_TO.gds")
output_dir = os.path.join(_ROOT, "example_output", "verification", "hcg_reflection_convergence")
data_dir = os.path.join(output_dir, "data")
_CACHE_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

CORE_LAYERS = [(2, 0)]
WAVELENGTH = 1.569  # C-band, the paper's own reported peak-reflectance wavelength
DPML = 0.5
DECAY_TOL = 1e-11   # the tutorial's production value -- not loosened here

# Every point dict is complete (res, sim_width, edge_clearance, sim_height) --
# no hidden defaults -- so SERIES below can read any axis off any point
# uniformly. See module docstring for how indices group into series.
POINTS = [
    dict(kind="resolution", res=10.0, sim_width=7.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="resolution", res=15.0, sim_width=7.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="resolution", res=20.0, sim_width=7.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="resolution", res=25.0, sim_width=7.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="domain_y",   res=10.0, sim_width=8.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="domain_y",   res=10.0, sim_width=9.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="resolution", res=30.0, sim_width=7.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="resolution", res=35.0, sim_width=7.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="resolution", res=40.0, sim_width=7.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="domain_x",   res=25.0, sim_width=7.5, edge_clearance=0.5, sim_height=1.5),
    dict(kind="domain_x",   res=25.0, sim_width=7.5, edge_clearance=0.8, sim_height=1.5),
    dict(kind="domain_z",   res=25.0, sim_width=7.5, edge_clearance=0.3, sim_height=2.0),
    dict(kind="domain_z",   res=25.0, sim_width=7.5, edge_clearance=0.3, sim_height=2.5),
    dict(kind="domain_y",   res=25.0, sim_width=8.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="domain_y",   res=25.0, sim_width=9.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="resolution", res=50.0, sim_width=7.5, edge_clearance=0.3, sim_height=1.5),
    dict(kind="resolution", res=60.0, sim_width=7.5, edge_clearance=0.3, sim_height=1.5),
]

# name -> (indices incl. shared baseline, swept param key, axis label, plot title)
SERIES = {
    "resolution": (
        (0, 1, 2, 3, 6, 7, 8, 15, 16), "res", "Resolution (pts/um)",
        "HCG reflection vs FDTD resolution\n(3-D, sim_width=7.5, edge_clearance=0.3, sim_height=1.5, lambda=1.569 um)",
    ),
    "domain_y_res10": (
        (0, 4, 5), "sim_width", "sim_width, y (um)",
        "HCG reflection vs y-domain width\n(3-D, res=10, edge_clearance=0.3, sim_height=1.5, lambda=1.569 um)",
    ),
    "domain_y_res25": (
        (3, 13, 14), "sim_width", "sim_width, y (um)",
        "HCG reflection vs y-domain width\n(3-D, res=25, edge_clearance=0.3, sim_height=1.5, lambda=1.569 um)",
    ),
    "domain_x_res25": (
        (3, 9, 10), "edge_clearance", "edge_clearance, x (um)",
        "HCG reflection vs x-domain (PML) clearance\n(3-D, res=25, sim_width=7.5, sim_height=1.5, lambda=1.569 um)",
    ),
    "domain_z_res25": (
        (3, 11, 12), "sim_height", "sim_height, z (um)",
        "HCG reflection vs z-domain height\n(3-D, res=25, sim_width=7.5, edge_clearance=0.3, lambda=1.569 um)",
    ),
}


def _is_master():
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
    """write_fn(tmp) then atomically rename onto path -- see
    beef.atomic_replace_write's docstring / README finding 5. Reimplemented
    locally (not calling beef.atomic_replace_write) only to keep the
    host-qualified temp name identical to verify_item_5_ground_truth.py's own
    copy, which this script otherwise mirrors line-for-line in shape."""
    tmp = f"{path}.tmp.{socket.gethostname()}.{os.getpid()}"
    write_fn(tmp)
    os.replace(tmp, path)


def _tag(i):
    return f"idx{i}"


def make_grating():
    raw = load_gds_polygons(_GDS_PATH, layers=CORE_LAYERS)
    xs = [x for poly in raw[CORE_LAYERS[0]] for x, _y in poly]
    span = max(xs) - min(xs)
    return GdsGeometry(gds_path=_GDS_PATH, core_layers=CORE_LAYERS, region_lengths=[span])


# source_margin/monitor_margin are measured from the PML INNER edge (see
# solve_reflection()'s own docstring), not from the geometry's edge -- so
# holding them at fixed absolute values while sweeping edge_clearance moves
# the source/monitor relative to the actual silicon. Caught 2026-08-04: task
# 9 (edge_clearance=0.5) landed with the source already 0.05 um OUTSIDE the
# 1.082 um input stub (R=0.7530, a large, physically meaningless deviation
# from baseline -- not real x-PML sensitivity), and task 10
# (edge_clearance=0.8) was cancelled mid-run once this was found, since BOTH
# source and monitor would have been outside the stub there too.
#
# Fix: derive source_margin/monitor_margin from edge_clearance so the
# source/monitor sit at a FIXED offset from the geometry's own left edge
# (+0.15 um / +0.50 um into the stub -- exactly the tutorial's own
# source_margin=0.45/monitor_margin=0.8 at the baseline edge_clearance=0.3),
# regardless of how much PML clearance surrounds the structure. At
# edge_clearance=0.3 this reproduces 0.45/0.8 exactly, so it changes nothing
# for every other point in this sweep (verified: cache hits, not new
# fingerprints) -- only points 9-10 actually move.
_SOURCE_OFFSET = 0.45 - 0.3   # um into the stub from its own edge
_MONITOR_OFFSET = 0.8 - 0.3


def make_solver(point):
    ec = point["edge_clearance"]
    return FDTDSolver(
        res=point["res"], sim_width=point["sim_width"], sim_height=point["sim_height"],
        dpml=DPML, source_margin=ec + _SOURCE_OFFSET, monitor_margin=ec + _MONITOR_OFFSET,
        three_d=True, num_input_modes=1, material_stack=SOI_160, decay_tol=DECAY_TOL,
    )


def provenance(cached):
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
        beef_version=getattr(beef, "__version__", "<none>"),
        beef_git=beef_git,
        meep_version=meep_version,
        python=platform.python_version(),
        host=platform.node(),
        mpi_ranks=nproc,
        slurm_job_id=os.environ.get("SLURM_JOB_ID", "<none>"),
        slurm_array_job_id=os.environ.get("SLURM_ARRAY_JOB_ID", "<none>"),
        cached=bool(cached),
        device="GdsGeometry(HCG_seeded_TO.gds)",
        wavelength=WAVELENGTH,
        dpml=DPML,
        decay_tol=DECAY_TOL,
        points=POINTS,
    )


# ---------------------------------------------------------------------------
# Per-point simulation
# ---------------------------------------------------------------------------
def simulate_point(i, cached=True):
    point = POINTS[i]
    grating = make_grating()
    solver = make_solver(point)
    say(f"\n[hcg-conv] idx={i} kind={point['kind']} res={point['res']} "
        f"sim_width={point['sim_width']} edge_clearance={point['edge_clearance']} "
        f"sim_height={point['sim_height']} (3-D, decay_tol={DECAY_TOL:g}, lambda={WAVELENGTH})")

    t0 = time.time()
    smat = solver.solve_reflection(
        grating, np.array([WAVELENGTH]), port="top", input_mode=1,
        edge_clearance=point["edge_clearance"],
        save_dir=_CACHE_DIR if cached else None,
    )
    elapsed = time.time() - t0
    smat.check_passive()
    R = float(smat.transmission(port="in", in_mode=1, out_mode=1)[0])

    if MASTER:
        os.makedirs(data_dir, exist_ok=True)
        _atomic_replace_write(
            os.path.join(data_dir, f"{_tag(i)}.npz"),
            lambda tmp: np.savez(
                open(tmp, "wb"), index=i, kind=point["kind"], res=point["res"],
                sim_width=point["sim_width"], edge_clearance=point["edge_clearance"],
                sim_height=point["sim_height"],
                R=R, solve_seconds=elapsed, cached=cached),
        )

    say(f"[hcg-conv] idx={i} R={R:.4f} solved in {elapsed / 60:.1f} min"
        f"{' (cache hit)' if elapsed < 5 else ''}")


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
def aggregate(cached=True):
    done = {}
    for i in range(len(POINTS)):
        path = os.path.join(data_dir, f"{_tag(i)}.npz")
        if not os.path.exists(path):
            continue
        d = np.load(path)
        done[i] = dict(res=float(d["res"]), sim_width=float(d["sim_width"]),
                        edge_clearance=float(d["edge_clearance"]),
                        sim_height=float(d["sim_height"]),
                        R=float(d["R"]), solve_seconds=float(d["solve_seconds"]))

    n_done, n_all = len(done), len(POINTS)
    if n_done < 1:
        say(f"[hcg-conv] Aggregation skipped -- 0/{n_all} points done.")
        return
    if not MASTER:
        return

    prov = provenance(cached)
    prov["points_completed"] = sorted(done.keys())
    prov_json = json.dumps(prov, indent=2, sort_keys=True)

    os.makedirs(output_dir, exist_ok=True)
    _atomic_replace_write(
        os.path.join(output_dir, "convergence.npz"),
        lambda tmp: np.savez(open(tmp, "wb"), done=json.dumps(done), provenance=np.array(prov_json)),
    )
    _atomic_replace_write(
        os.path.join(output_dir, "provenance.json"),
        lambda tmp: open(tmp, "w").write(prov_json + "\n"),
    )
    say(f"\n[hcg-conv] Aggregate saved ({n_done}/{n_all} points) -> {output_dir}/convergence.npz")

    for name, (indices, param, xlabel, title) in SERIES.items():
        idx = [i for i in indices if i in done]
        if not idx:
            continue
        say(f"\n[hcg-conv] Series '{name}' (swept: {param})")
        say(f"[hcg-conv]   {param:14s} R        |dR| vs prev    wall (min)")
        prev = None
        for i in idx:
            r = done[i]
            dR = "" if prev is None else f"{abs(r['R'] - prev):.4f}"
            say(f"[hcg-conv]   {r[param]:<14.3f} {r['R']:.4f}   {dR:>10s}         {r['solve_seconds']/60:.1f}")
            prev = r["R"]
        _plot_series(
            [done[i][param] for i in idx], [done[i]["R"] for i in idx],
            xlabel, f"{name}_convergence.png", title,
            n_done=len(idx), n_all=len(indices),
        )


def _plot_series(x, y, xlabel, filename, title, n_done, n_all):
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    ax.plot(x, y, marker="o", color="steelblue", lw=1.6)
    if n_done < n_all:
        ax.text(0.02, 0.97, f"{n_done}/{n_all} points done", transform=ax.transAxes,
                ha="left", va="top", fontsize=8, color="#555555",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Reflection R")
    ax.set_title(title, fontsize=10)
    ax.grid(True)
    plt.tight_layout()
    _atomic_replace_write(os.path.join(output_dir, filename),
                          lambda tmp: fig.savefig(tmp, dpi=300, format="png"))
    plt.close("all")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("index", nargs="?", type=int, default=None,
                        help=f"point index 0..{len(POINTS) - 1}; omit to run every point serially")
    parser.add_argument("--aggregate-only", action="store_true",
                        help="re-aggregate finished points, simulate nothing")
    parser.add_argument("--no-cache", action="store_true",
                        help="ignore and do not write the simulation cache")
    args = parser.parse_args()
    cached = not args.no_cache

    if MASTER:
        os.makedirs(data_dir, exist_ok=True)

    if args.aggregate_only:
        say("[hcg-conv] Aggregate-only: no simulation.")
        aggregate(cached)
        sys.exit(0)

    if args.index is not None:
        if not 0 <= args.index < len(POINTS):
            raise SystemExit(
                f"index {args.index} out of range -- sweep has "
                f"{len(POINTS)} points (0..{len(POINTS) - 1})")
        simulate_point(args.index, cached)
    else:
        say(f"[hcg-conv] No index given -- running all {len(POINTS)} points serially "
            "(the SLURM array batch script is much faster).")
        for i in range(len(POINTS)):
            simulate_point(i, cached)

    aggregate(cached)
