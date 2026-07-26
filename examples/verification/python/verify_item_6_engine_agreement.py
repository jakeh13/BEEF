"""TEMPORARY — direct FDTD/EME engine-agreement test for dc_validation.md item 6.

Delete this folder (and the Makefile's verification targets, and
batch/verification/) once the directional coupler work is signed off.

THE QUESTION
------------
Items 1-2 fixed the FDTD wavelength axis, and item 5 re-established the device
ground truth. One discrepancy survives both, and it was never explained by the
axis bug because it was measured at lambda = 1.55 um -- the single index the
pre-fix code labelled correctly.

The stub sweep measured it indirectly. In the two-mode picture,

    L_50(s) = [pi/2 - 2*phi_bend]/dbeta_EME + 2s*(1 - dbeta_FDTD/dbeta_EME)

so the sweep's dead-straight +0.2185 um/um drift of L_50 with the stub s implies
dbeta_FDTD/dbeta_EME = 0.891: the FDTD sections accumulate ~11% less even/odd
phase per um of the SAME straight coupled region than the EME body does. That
number came out of a linear fit to a derived quantity, through a whole hybrid
chain, with bend phase and mode-matching folded in.

This script measures the same ratio directly, on nothing but a uniform coupled
pair, so the answer does not depend on the chain, the bends, or L_50 extraction.

HOW
---
A DirectionalCoupler's coupling region is a straight, uniform, parallel pair of
waveguides. Cut a section entirely inside it and both engines are being asked
the identical question: how much even/odd phase accrues over L um of this
cross-section?

  FDTD   solve_section(x0=-L/2, x1=+L/2, n_modes_in=n_modes_out=2) returns a
         modal transfer matrix. A uniform region cannot convert between
         supermodes, so the matrix must be diagonal (that is checked, not
         assumed -- see off_diag_max) with entries ~ exp(i*beta_even*L) and
         exp(i*beta_odd*L). Then

             dphi(L) = arg(T_00) - arg(T_11) = dbeta*L + const

         The constant is mode-matching and monitor-plane phase, identical at
         every length. It is why ONE length cannot answer this: dphi(L)/L would
         carry that offset as a 1/L error. Four lengths and a straight-line fit
         remove it and leave the residual as evidence the model holds.

  EME    the same section_gds() polygons -> the same meow cross-section ->
         compute_modes(). dbeta_EME = 2*pi*(n_eff_even - n_eff_odd)/lambda,
         read straight off the eigenvalues with no propagation involved.

Both engines are pinned to directional_coupler_hybrid_stub_sweep.py's own
solver settings, because reproducing THAT sweep's 0.891 is the point. Note
mode_solver="meep" on the EME side: meow's MPB-backed solver is the same
eigenmode machinery FDTD's get_eigenmode_coefficients() uses, so a mismatch
here cannot be blamed on two different mode solvers.

WHAT COUNTS AS AN ANSWER
------------------------
  ratio ~ 0.89   the stub sweep's gap is real and lives in the straight
                 coupled-pair section, not in the bends or the chain. That
                 localizes it to how FDTD accumulates supermode phase.
  ratio ~ 1.00   the gap is NOT in the uniform region; it enters through the
                 bends, the mode matching at the seam, or L_50 extraction, and
                 item 6 should move there.

Either outcome is informative, which is why this is worth the machine time.

Run:
    # one section length per array task:
    sbatch batch/verification/run_verify_item_6_engine_agreement.batch

    # a single length by index, or all of them serially:
    python .../verify_item_6_engine_agreement.py <index 0..3>

    # re-aggregate + run the (cheap, serial) EME half only:
    OMP_NUM_THREADS=1 python .../verify_item_6_engine_agreement.py --aggregate-only

Outputs (example_output/verification/engine_agreement/):

    data/L{tag}.npz          per-length FDTD transfer matrix + extracted phases
    engine_agreement.npz     dbeta_FDTD(lambda), dbeta_EME(lambda), the ratio,
                             fit residuals, and embedded provenance
    provenance.json          the same provenance, standalone and greppable
    dbeta_agreement.png      the headline: both dbeta(lambda) and their ratio
    phase_fits.png           dphi vs L per wavelength -- the linearity evidence
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
from beef import DirectionalCoupler, EMESolver, FDTDSolver, SOI_220, extrude_polygons

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "engine_agreement")
data_dir = os.path.join(output_dir, "data")

# Caching ON, same reasoning as verify_item_5_ground_truth.py: item 3 completed
# the FDTD fingerprint so a stale entry can no longer be hit, and item 4 made the
# writes rank-safe. --no-cache forces a from-scratch re-derivation.
_CACHE_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Pinned to directional_coupler_hybrid_stub_sweep.py -- do not "improve" these.
# The whole point is to reproduce the 0.891 that sweep implied, so any solver
# setting that differs from it would make a disagreement unattributable.
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

DEVICE_KWARGS = dict(
    waveguide_width=0.45, gap=0.25, coupling_length=30.0,
    bend_length=9.0, waveguide_offset=3.0, bend_type="euler",
    waveguide_height=0.22,
)
FDTD_KWARGS = dict(
    res=RESOLUTION, sim_width=6.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    half_domain_aperture=True, material_stack=SOI_220, decay_tol=1e-11,
)
EME_KWARGS = dict(
    num_modes=[6], num_cells=[1], lengths=[1.0],   # placeholder; unused here
    oxide_width=6.0, oxide_height=1.5, res=RESOLUTION,
    material_stack=SOI_220, mode_solver="meep",
)

# Two supermodes: that is the entire model being tested. Asking for more would
# measure something else.
N_SUPERMODES = 2

# Section lengths (um), each centred at x=0 so every cut sits deep inside the
# 30 um coupling region -- the nearest S-bend seam is 15 um away even for the
# longest. Four lengths, not two: the slope is what is wanted, and the fit
# residual across four points is the evidence that dphi really is linear in L
# (i.e. that this is one clean supermode pair and not a leaky mess).
#
# The span is bounded below by needing dphi to move well clear of its own
# uncertainty, and above by 2*pi wrapping: dbeta ~ 0.103 rad/um at 1.55 um
# (from L_50 = 15.3 um), so 8 um accrues 0.82 rad -- comfortably inside pi, so
# the unwrap is never ambiguous.
SECTION_LENGTHS = np.array([2.0, 4.0, 6.0, 8.0])

# Interface-monitor pin, identical to the stub sweep's _iface_mon: the full
# non-PML domain centred on the pair's shared y=0 axis, so the supermode basis
# is the symmetric one both engines agree on. Without this the monitor inherits
# one arm's centreline and the "even/odd" decomposition is not even/odd at all.
IFACE_MON = dict(monitor_center=0.0, monitor_size=FDTD_KWARGS["sim_width"])


# ---------------------------------------------------------------------------
# MPI rank guarding. All ranks run the whole script; only one may write.
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
    """write_fn(tmp) then atomically rename onto path. Host-qualified temp name
    because PIDs are allocated per host (dc_validation.md item 4).

    Note the temp name deliberately does NOT end in the final extension, so
    np.savez() callers must hand it an already-open file object -- given a
    path string it appends ".npz" itself and the rename below then fails on a
    file that isn't there. Same reason fig.savefig() callers pass format=.
    """
    tmp = f"{path}.tmp.{socket.gethostname()}.{os.getpid()}"
    write_fn(tmp)
    os.replace(tmp, path)


def _tag(L):
    return f"L{L:04.1f}um"


def make_coupler():
    return DirectionalCoupler(**DEVICE_KWARGS)


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
        meep_version, ranks = mp.__version__, int(mp.count_processors())
    except Exception:
        meep_version, ranks = "<not loaded>", 1
    return {
        "beef_git": beef_git,
        "beef_version": beef.__version__,
        "meep_version": meep_version,
        "python": platform.python_version(),
        "host": platform.node(),
        "mpi_ranks": ranks,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "cached": bool(cached),
        "device": "DirectionalCoupler",
        "device_kwargs": DEVICE_KWARGS,
        "fdtd_kwargs": {k: (v.__class__.__name__ if k == "material_stack" else v)
                        for k, v in FDTD_KWARGS.items()} | {"material_stack": "SOI_220"},
        "eme_kwargs": {k: v for k, v in EME_KWARGS.items() if k != "material_stack"}
                      | {"material_stack": "SOI_220"},
        "n_supermodes": N_SUPERMODES,
        "section_lengths": SECTION_LENGTHS.tolist(),
        "wavelengths": wavelengths.tolist(),
        "iface_monitor": IFACE_MON,
        "measures": "dbeta_FDTD/dbeta_EME for a uniform coupled pair "
                    "(dc_validation.md item 6; stub sweep implied 0.891)",
    }


# ---------------------------------------------------------------------------
# The FDTD half: one section length per call.
# ---------------------------------------------------------------------------
def simulate_length(L, cached=True):
    """Transfer matrix for a uniform coupled-pair section of length L."""
    dc = make_coupler()
    solver = FDTDSolver(**FDTD_KWARGS)
    x0, x1 = -L / 2.0, +L / 2.0
    half_body = dc.coupling_length / 2.0
    assert abs(x0) < half_body and abs(x1) < half_body, (
        f"section [{x0}, {x1}] escapes the {dc.coupling_length} um coupling "
        f"region -- it would include S-bend, which is exactly what this "
        f"measurement must exclude"
    )

    say(f"[ea] section L = {L} um, x in [{x0:+.2f}, {x1:+.2f}] "
        f"(coupling region +/-{half_body:.1f} um), {N_SUPERMODES} supermodes, "
        f"{len(wavelengths)} wavelengths")
    t0 = time.time()
    T = solver.solve_section(
        dc, wavelengths, x0, x1,
        n_modes_in=N_SUPERMODES, n_modes_out=N_SUPERMODES,
        save_dir=_CACHE_DIR if cached else None,
        **IFACE_MON,
    )
    dt = time.time() - t0

    # A uniform region cannot convert even <-> odd. If it appears to, the
    # monitor basis is wrong (usually a mis-centred interface window) and every
    # phase below is meaningless -- so this is recorded, not assumed.
    diag = np.abs(np.stack([T[:, 0, 0], T[:, 1, 1]], axis=1))
    off = np.abs(np.stack([T[:, 0, 1], T[:, 1, 0]], axis=1))
    dphi = np.angle(T[:, 0, 0]) - np.angle(T[:, 1, 1])

    say(f"[ea] solved in {dt/60:.1f} min | "
        f"|diag| {diag.min():.3f}-{diag.max():.3f} | max |off-diag| {off.max():.4f}")

    if MASTER:
        os.makedirs(data_dir, exist_ok=True)
        def _write(tmp):
            with open(tmp, "wb") as fh:
                np.savez(fh, section_length=L, wavelengths=wavelengths,
                         T=T, dphi=dphi, diag_abs=diag, off_diag_abs=off,
                         solve_seconds=dt, cached=cached)

        _atomic_replace_write(os.path.join(data_dir, f"{_tag(L)}.npz"), _write)
    return T


# ---------------------------------------------------------------------------
# The EME half: cheap, serial, no MPI. Rebuilt from the same section_gds() the
# FDTD side draws, so the two engines see byte-identical geometry.
# ---------------------------------------------------------------------------
def eme_dbeta():
    """(dbeta_EME[wl], neff[wl, mode]) for the uniform coupled-pair cross-section."""
    import meow as mw

    dc = make_coupler()
    eme = EMESolver(**EME_KWARGS)
    # Any cut inside the coupling region gives the same cross-section; length is
    # irrelevant to an eigenvalue, so use the shortest for the cheapest cell.
    L = float(SECTION_LENGTHS[0])
    gds = dc.section_gds(-L / 2.0, +L / 2.0, eme.oxide_width,
                         stub_left=0.0, stub_right=0.0)
    structs = extrude_polygons(gds, eme.material_stack.extrusion_rules)
    cells = mw.create_cells(structures=structs, mesh=eme._build_mesh(), Ls=[L])

    neff = np.zeros((len(wavelengths), N_SUPERMODES))
    for i, wl in enumerate(wavelengths):
        env = mw.Environment(wl=float(wl), T=eme.temperature)
        cs = mw.CrossSection.from_cell(cell=cells[0], env=env)
        modes = eme.compute_modes(cs, N_SUPERMODES)
        neff[i] = [float(np.real(m.neff)) for m in modes]

    # meow returns modes in descending n_eff, so column 0 is the symmetric
    # (even) supermode. Asserted rather than trusted: if it ever reorders, the
    # sign of dbeta silently flips and the ratio looks fine.
    assert np.all(neff[:, 0] > neff[:, 1]), (
        "EME supermodes are not in descending n_eff order -- the even/odd "
        f"identification below would be wrong. neff = {neff}"
    )
    dbeta = 2 * np.pi * (neff[:, 0] - neff[:, 1]) / wavelengths
    return dbeta, neff


# ---------------------------------------------------------------------------
def fit_dbeta_fdtd(lengths, dphi):
    """Slope of dphi vs L per wavelength, with the fit residual.

    dphi(L) = dbeta*L + const. The constant is mode-matching/monitor-plane
    phase, identical at every length, so the SLOPE is dbeta and the intercept
    is discarded. Fitting rather than differencing two points uses all four and
    leaves a residual that says whether the linear model held.
    """
    lengths = np.asarray(lengths, dtype=float)
    dbeta = np.zeros(dphi.shape[1])
    resid = np.zeros(dphi.shape[1])
    intercept = np.zeros(dphi.shape[1])
    for w in range(dphi.shape[1]):
        # Unwrap along L before fitting: arg() is confined to (-pi, pi] and a
        # genuine dphi crossing that boundary would otherwise look like a jump.
        y = np.unwrap(dphi[:, w])
        A = np.vstack([lengths, np.ones_like(lengths)]).T
        (m, c), *_ = np.linalg.lstsq(A, y, rcond=None)
        dbeta[w], intercept[w] = m, c
        resid[w] = float(np.sqrt(np.mean((y - (m * lengths + c)) ** 2)))
    return dbeta, intercept, resid


def aggregate(cached):
    files = {}
    for L in SECTION_LENGTHS:
        p = os.path.join(data_dir, f"{_tag(L)}.npz")
        if os.path.exists(p):
            files[float(L)] = np.load(p, allow_pickle=True)
    if len(files) < 2:
        say(f"[ea] only {len(files)}/{len(SECTION_LENGTHS)} lengths done -- "
            f"need >=2 to fit a slope, skipping aggregate")
        return
    done = np.array(sorted(files))
    dphi = np.stack([files[L]["dphi"] for L in done])            # [nL, nwl]
    off_max = max(float(files[L]["off_diag_abs"].max()) for L in done)

    dbeta_f, intercept, resid = fit_dbeta_fdtd(done, dphi)
    say(f"[ea] running the EME half ({len(wavelengths)} mode solves)...")
    dbeta_e, neff = eme_dbeta()
    ratio = dbeta_f / dbeta_e

    prov = provenance(cached)
    prov["lengths_completed"] = int(len(done))
    prov["max_off_diagonal"] = off_max
    prov_json = json.dumps(prov, indent=2, sort_keys=True)

    if MASTER:
        os.makedirs(output_dir, exist_ok=True)

        def _write_npz(tmp):
            with open(tmp, "wb") as fh:
                np.savez(fh, section_lengths=done, wavelengths=wavelengths,
                         dphi=dphi, dbeta_fdtd=dbeta_f, dbeta_eme=dbeta_e,
                         ratio=ratio, fit_intercept=intercept,
                         fit_residual=resid, neff_eme=neff,
                         max_off_diagonal=off_max,
                         provenance=np.array(prov_json))

        def _write_json(tmp):
            with open(tmp, "w") as fh:
                fh.write(prov_json + "\n")

        _atomic_replace_write(os.path.join(output_dir, "engine_agreement.npz"), _write_npz)
        _atomic_replace_write(os.path.join(output_dir, "provenance.json"), _write_json)

    say(f"\n[ea] {len(done)}/{len(SECTION_LENGTHS)} lengths: "
        f"{', '.join(f'{L:.1f}' for L in done)} um")
    say(f"[ea] max |off-diagonal| over all lengths/wavelengths = {off_max:.4f} "
        f"({'OK -- no supermode conversion' if off_max < 0.05 else 'HIGH -- basis suspect'})")
    say("\n[ea] lambda   dbeta_FDTD   dbeta_EME     ratio    fit resid (rad)")
    for i, wl in enumerate(wavelengths):
        say(f"[ea]  {wl:.3f}     {dbeta_f[i]:8.5f}   {dbeta_e[i]:8.5f}   "
            f"{ratio[i]:7.4f}     {resid[i]:.4f}")
    j = int(np.argmin(np.abs(wavelengths - 1.55)))
    say(f"\n[ea] AT 1.55 um: ratio = {ratio[j]:.4f}   "
        f"(stub sweep implied 0.891)")
    say(f"[ea] band mean ratio = {ratio.mean():.4f} +/- {ratio.std():.4f}")
    if abs(ratio[j] - 1.0) < 0.02:
        say("[ea] => engines AGREE on the uniform region. The stub sweep's gap "
            "enters through the bends, the seam mode-matching, or L50 "
            "extraction -- not through straight-section phase accumulation.")
    elif abs(ratio[j] - 0.891) < 0.03:
        say("[ea] => reproduces the stub sweep's 0.891. The gap is localized to "
            "how FDTD accumulates supermode phase in a uniform coupled pair.")
    else:
        say("[ea] => neither 1.00 nor 0.891; the two-mode model behind the "
            "stub sweep's inference may not hold. Needs interpretation.")

    if MASTER:
        _plots(done, dphi, dbeta_f, intercept, dbeta_e, ratio)
        say(f"\n[ea] Written -> {output_dir}/")


def _plots(done, dphi, dbeta_f, intercept, dbeta_e, ratio):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
    ax = axes[0]
    ax.plot(wavelengths, dbeta_f, "o-", color="#1f77b4", label=r"$\Delta\beta$ FDTD")
    ax.plot(wavelengths, dbeta_e, "s-", color="#d62728", label=r"$\Delta\beta$ EME")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel(r"$\Delta\beta$ (rad/µm)")
    ax.set_title("Even/odd phase rate, uniform coupled pair")
    ax.legend(); ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(wavelengths, ratio, "o-", color="#2ca02c")
    ax.axhline(1.0, color="0.4", lw=1, ls="--", label="perfect agreement")
    ax.axhline(0.891, color="#ff7f0e", lw=1, ls=":",
               label="0.891 (stub sweep inference)")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel(r"$\Delta\beta_{\rm FDTD} / \Delta\beta_{\rm EME}$")
    ax.set_title("Engine agreement")
    ax.legend(); ax.grid(alpha=0.3)
    fig.savefig(os.path.join(output_dir, "dbeta_agreement.png"), dpi=150)
    plt.close(fig)

    # The linearity evidence: if these are not straight, the slope above is not
    # a phase rate and the whole measurement is void.
    fig, ax = plt.subplots(figsize=(6.5, 4.6), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    Lf = np.linspace(0, done.max() * 1.05, 50)
    for w, wl in enumerate(wavelengths):
        c = cmap(w / max(len(wavelengths) - 1, 1))
        ax.plot(done, np.unwrap(dphi[:, w]), "o", color=c, ms=5)
        ax.plot(Lf, dbeta_f[w] * Lf + intercept[w], "-", color=c, lw=1,
                label=f"{wl:.2f} µm" if w % 5 == 0 else None)
    ax.set_xlabel("Section length L (µm)")
    ax.set_ylabel(r"$\Delta\phi$ (rad)")
    ax.set_title(r"$\Delta\phi = \Delta\beta L + $ const — slope is the measurement")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.savefig(os.path.join(output_dir, "phase_fits.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("index", nargs="?", type=int, default=None,
                    help=f"section-length index 0..{len(SECTION_LENGTHS)-1}; "
                         f"omit to run all serially")
    ap.add_argument("--aggregate-only", action="store_true",
                    help="no FDTD; re-fit what exists and run the EME half")
    ap.add_argument("--no-cache", action="store_true",
                    help="ignore the on-disk cache and re-simulate")
    args = ap.parse_args()
    cached = not args.no_cache

    say(f"[ea] item 6 — direct FDTD/EME engine agreement")
    say(f"[ea]   device    DirectionalCoupler(coupling_length="
        f"{DEVICE_KWARGS['coupling_length']} µm), uniform region only")
    say(f"[ea]   lengths   {', '.join(f'{L:.1f}' for L in SECTION_LENGTHS)} µm")
    say(f"[ea]   modes     {N_SUPERMODES} (even/odd supermodes)")
    say(f"[ea]   cache     {'on' if cached else 'OFF (--no-cache)'}\n")

    if not args.aggregate_only:
        todo = (SECTION_LENGTHS if args.index is None
                else [SECTION_LENGTHS[args.index]])
        for L in todo:
            simulate_length(float(L), cached=cached)

    aggregate(cached)
    return 0


if __name__ == "__main__":
    sys.exit(main())
