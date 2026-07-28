"""
Tutorial: directional coupler — hybrid FDTD/EME body-stub convergence sweep.

The same 220 nm SOI, C-band device and method as directional_coupler_hybrid.py
(2 FDTD S-bend sections + 1 EME body, N_MODES=6, FDTD 40 / EME 80 px/µm) on ONE fixed
directional coupler, but sweeping where the FDTD/EME section interfaces are
cut:

    BODY_STUB = 0.00, 0.25, 0.50, ..., 2.00 µm   (9 values)

BODY_STUB is how much of the straight coupled-pair region each FDTD section
folds in past the S-bend/coupling-region seam, i.e. how far the interface
mode-decomposition plane sits inside the uniform coupling region:

    |<- bend ->|<-------------- coupling_length -------------->|<- bend ->|
    [    FDTD input section   ][      EME body      ][  FDTD output sect. ]
               |<-BODY_STUB->|                      |<-BODY_STUB->|

BODY_STUB = 0 is the old exact-seam convention: the supermode decomposition
happens right on the geometric discontinuity, where the bend's near field has
had no straight run to relax into the coupled-pair basis. Each larger stub
gives it more relaxation length — at the price of a slightly bigger FDTD
domain (Sx grows by BODY_STUB) and a correspondingly shorter EME body.

The physical device never changes: for every stub the EME body models
coupling_length - 2*BODY_STUB, so all three sections together always cover
exactly the same coupling_length (asserted per stub in simulate_stub()). Any
movement in the result is therefore method error, not a different device —
which is exactly what this sweep measures.

The headline output is the 50:50 (3 dB) split point: the coupling length at
which the mean-over-wavelength cross/(cross+through) ratio crosses 0.5,
extracted per stub by linear interpolation between the two bracketing sampled
lengths. A converged interface placement shows that length flattening out as
BODY_STUB grows; the flat value is the hybrid method's answer, and where the
curve stops moving is the smallest stub that buys it.

Cost per array task: 2 FDTD section solves + 1 EME mode-solve per wavelength
(the whole coupling-length axis rides that single EME solve analytically —
see EMESolver.solve_section_length_sweep()'s docstring), i.e. the same cost as
one full run of directional_coupler_hybrid.py's step 5. Each stub is its own
section geometry, so nothing is shared between tasks (every stub really does
re-run its own FDTD sections).

Outputs (example_output/directional_coupler/directional_coupler_hybrid_stub_sweep/):

    plots/coupling_length_sweep_S{NNN}.png   per-stub power-transfer-vs-length
                                              plot (same layout as hybrid.py's
                                              coupling_length_sweep.png)
    data/S{NNN}.npz                          per-stub raw length-sweep grid
                                              (+ that stub's 50:50 length)
    stub_sweep_split_point.png               aggregate two-panel plot: 50:50
                                              coupling length vs BODY_STUB,
                                              plus the S-curve family
    stub_sweep.npz                           aggregate grid over completed stubs

Run:
    # one BODY_STUB per SLURM array task (see
    # batch/sweep/run_directional_coupler_hybrid_stub_sweep.batch):
    python /path/to/examples/directional_coupler/python/directional_coupler_hybrid_stub_sweep.py <index 0..8>

    # or the whole sweep serially (slow — batch queue strongly preferred):
    python /path/to/examples/directional_coupler/python/directional_coupler_hybrid_stub_sweep.py

Every invocation ends by (re)aggregating whatever stubs have finished so far,
writing the aggregate plot/data atomically — so concurrent array tasks can't
corrupt them, the plot is usable (annotated "n/9 stubs done") while the array
is still draining, and the last task to finish leaves the complete sweep
behind. Re-running an index simply overwrites that stub's files.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import DirectionalCoupler, EMESolver, FDTDSolver, SOI_220, atomic_replace_write

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "directional_coupler", "directional_coupler_hybrid_stub_sweep"))
plot_dir = os.path.join(output_dir, "plots")
data_dir = os.path.join(output_dir, "data")
for _d in (output_dir, plot_dir, data_dir):
    os.makedirs(_d, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "example_sim_cache"))

# The brute-force full-device FDTD length sweep, when it has been run — its own
# 50:50 length is the single-domain ground truth this sweep's stubs are trying
# to converge onto, drawn as a reference line on panel 1 when present.
#
# This reads the *verification* sweep, not the tutorial one next door: the
# tutorial's 0:2:60 grid interpolates L_50 across a 2 µm bracket, too coarse to
# judge a sweep resolving 0.06 µm steps, and its npz carries no provenance. The
# verification sweep is 0.5 µm through the crossing and records its own solver
# fingerprint (dc_validation.md item 5). Key names are identical, so the reader
# below is unchanged.
_FDTD_SWEEP_NPZ = os.path.abspath(os.path.join(
    output_dir, "..", "..", "verification", "ground_truth_sweep", "length_sweep.npz"))

# ---------------------------------------------------------------------------
# Global simulation parameters — identical to directional_coupler_hybrid.py,
# except BODY_STUB is swept instead of pinned at 0.25 µm.
# ---------------------------------------------------------------------------
# The two engines get SEPARATE resolutions, and this is deliberate — a single
# shared RESOLUTION was how a 7% error hid for months (dc_validation.md item 11).
#
# meep's FDTD is converged at 40: its implied Δβ moves only 2.74% across res
# 30–90, non-monotonically (noise, not a trend). meow/MPB on its own Mesh2D is
# NOT — Δβ_EME at 1.55 µm runs 0.1125 (res 20) → 0.0957 (40) → 0.0926 (60) →
# 0.0894 (80) → 0.0895 (100), so 40 overstates the even/odd splitting by 7% and
# 80 is the first converged point. meep's subpixel smoothing is the likely
# reason it converges so much earlier on the same nominal grid.
#
# Measured directly in dc_validation.md Done item 6 (job 11489846) and pinned by
# tests/test_eme_convergence.py. Do not collapse these back into one knob.
FDTD_RESOLUTION = 40.0   # pixels/µm — converged; see above
EME_RESOLUTION = 80.0    # pixels/µm — 40 is NOT converged; see above
WL_MIN = 1.50       # wavelength sweep start (µm)
WL_MAX = 1.60       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

N_MODES = 6   # even + odd supermodes carried at the two body interfaces — hybrid.py's default, pinned here

STUB_SWEEP = np.arange(0.0, 2.0 + 1e-9, 0.25)   # µm — 0:0.25:2, 9 values — the swept knob

# One fixed reference device for every stub. The S-bend/access-waveguide
# geometry (and hence the FDTD sections) is defined relative to the
# coupling-region edges, so this length only has to be comfortably longer than
# the largest cut the sweep makes (2*STUB_SWEEP[-1]); it is not the coupling
# length being characterized — that axis is swept analytically below.
DC_COUPLING_LENGTH = 30.0   # µm

# One shared coupling-length axis for every stub, so the extracted 50:50
# lengths are compared on identical sampling. Its floor is set by the LARGEST
# stub in the sweep: a device shorter than that stub's two FDTD sections
# put together can't be split into these three sections at all.
coupling_lengths = np.linspace(2 * STUB_SWEEP[-1], 60.0, 200)   # µm
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)


def _tag(stub: float) -> str:
    """Stable per-stub file tag in centi-µm, e.g. S000 / S025 / S200."""
    return f"S{int(round(stub * 100)):03d}"


def split_point(T_cross_grid: np.ndarray, T_total_grid: np.ndarray) -> float:
    """Coupling length (µm) of the first 50:50 cross/through split.

    Linear interpolation between the two sampled lengths that bracket the
    first crossing of ratio = 0.5, rather than the nearest sampled length:
    the shared grid's step (~0.28 µm) is comparable to the stub-to-stub shifts
    this sweep is trying to resolve, so snapping to the grid would quantize
    the very trend being plotted. NaN if the ratio never crosses 0.5.
    """
    ratio = T_cross_grid.mean(axis=1) / np.clip(T_total_grid.mean(axis=1), 1e-20, None)
    crossings = np.flatnonzero(np.diff(np.sign(ratio - 0.5)))
    if not crossings.size:
        return float("nan")
    i = int(crossings[0])
    r0, r1 = ratio[i], ratio[i + 1]
    L0, L1 = coupling_lengths[i], coupling_lengths[i + 1]
    if r1 == r0:
        return float(L0)
    return float(L0 + (0.5 - r0) * (L1 - L0) / (r1 - r0))


# ---------------------------------------------------------------------------
# 1. Geometry and FDTD/EME solvers — shared across every stub value (neither
#    depends on where the sections are cut, only on device/mesh parameters).
# ---------------------------------------------------------------------------
dc = DirectionalCoupler(
    waveguide_width=0.45,                 # single-mode strip width (µm)
    gap=0.25,                             # coupling-region edge-to-edge gap (µm)
    coupling_length=DC_COUPLING_LENGTH,   # straight parallel interaction length (µm)
    bend_length=9.0,                      # horizontal extent of each S-bend (µm)
    waveguide_offset=3.0,                 # far-end centre-to-centre port pitch (µm)
    bend_type="euler",                    # "euler" (cubic-Bezier S-bend) or "circular" (matched arcs)
    waveguide_height=0.22,                # Si layer thickness (µm) — SOI_220
)

_fdtd_common = dict(
    res=FDTD_RESOLUTION,
    sim_width=6.0,
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    half_domain_aperture=True,
    material_stack=SOI_220,
    decay_tol=1e-11,
)
fdtd_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholder length; solve_section_length_sweep() sweeps its own
    oxide_width=6.0, oxide_height=1.5, res=EME_RESOLUTION, material_stack=SOI_220,
    mode_solver="meep",
)

half_body = dc.coupling_length / 2
half_total = dc.total_length / 2

# Interface-monitor pin for the input-section FDTD runs — see
# directional_coupler_hybrid.py's own comment for the full reasoning: spans
# the full non-PML domain centered on the coupled pair's shared y=0 axis, so
# it matches the EME body's own y=0-centered supermode basis. Identical for
# every stub, which is what keeps the stubs comparable: the ONLY thing this
# sweep varies is where along x that monitor plane sits.
_iface_mon = dict(monitor_center=0.0, monitor_size=fdtd_solver.sim_width)

# ---------------------------------------------------------------------------
# DC_CALIBRATE_BODY=1 — propagate the EME body on the β the FDTD sections
# actually accumulate, instead of on meow's mode index (dc_validation.md item
# 6). This is THE experiment this sweep exists to settle.
#
# The drift this sweep measures is exactly 2s·(1 - Δβ_FDTD/Δβ_EME), so
# substituting Δβ_EME := Δβ_FDTD sends it to zero identically. That is a
# prediction, not a hope: if the measured slope does not collapse from 0.1085
# to ~0, the mechanism in item 6 is wrong.
#
# Off by default so the uncalibrated sweep (job 11500412) stays reproducible;
# the two write to different folders so neither overwrites the other.
CALIBRATE_BODY = bool(int(os.environ.get("DC_CALIBRATE_BODY", "0")))
if CALIBRATE_BODY:
    # All three, not just output_dir: plot_dir/data_dir were derived from it
    # near the top of the file, so redirecting only output_dir would leave the
    # per-stub npz and PNGs overwriting the uncalibrated sweep's own.
    # _FDTD_SWEEP_NPZ is deliberately NOT redirected — the full-device
    # reference line is shared, and is what both sweeps are measured against.
    output_dir = output_dir + "_calibrated"
    plot_dir = os.path.join(output_dir, "plots")
    data_dir = os.path.join(output_dir, "data")
    for _d in (output_dir, plot_dir, data_dir):
        os.makedirs(_d, exist_ok=True)

# The body cross-section is the same uniform coupled pair for every stub — the
# cuts move along x but never leave the straight region — so β is stub
# independent and this is computed once, not once per array task. Cached, so
# concurrent tasks pay for it at most once between them.
_body_neff = None
if CALIBRATE_BODY:
    _beta = fdtd_solver.section_mode_beta(
        dc, wavelengths, -half_body, +half_body, N_MODES,
        save_dir=save_dir, **_iface_mon,
    )
    _body_neff = _beta * wavelengths[:, None] / (2 * np.pi)
    # Short window about the body midpoint, NOT the full [-half_body,
    # +half_body] span: section_mode_neff() pads the cut by dpml on each side,
    # so the full span asks for x=-15.5, which lands inside an S-bend and
    # raises. section_mode_beta() seeds itself the same way for the same
    # reason. The cross-section is uniform across the whole body, so any
    # interior window gives the same modes — this one is only a log line, and
    # it is not worth failing the job over.
    _seed_half = min(1.0, half_body - 2.0 * fdtd_solver.dpml)
    _mpb = fdtd_solver.section_mode_neff(
        dc, wavelengths, -_seed_half, +_seed_half, N_MODES, **_iface_mon)
    _i = int(np.argmin(np.abs(wavelengths - 1.55)))
    print(f"[calib] body n_eff from the FDTD's own phase, at 1.55 µm:")
    print(f"[calib]   MPB   {np.array2string(_mpb[_i, :2], precision=6)}  "
          f"Δn_eff {_mpb[_i, 0] - _mpb[_i, 1]:.6f}")
    print(f"[calib]   FDTD  {np.array2string(_body_neff[_i, :2], precision=6)}  "
          f"Δn_eff {_body_neff[_i, 0] - _body_neff[_i, 1]:.6f}")

# Both cuts have to stay inside the reference device's own coupling region
# (DirectionalCoupler rejects a cut that lands inside an S-bend), so the two
# stubs together may not exceed its coupling_length. BODY_STUB = 0 is the
# deliberate degenerate baseline: decomposition exactly on the seam.
assert 2 * STUB_SWEEP[-1] <= dc.coupling_length, (
    f"largest stub ({STUB_SWEEP[-1]:.2f} µm) x2 exceeds the reference device's "
    f"coupling_length ({dc.coupling_length:.2f} µm)"
)

print(f"Directional coupler total length: {dc.total_length:.2f} µm "
      f"(coupling length {dc.coupling_length:.2f} µm, coupling pitch "
      f"{dc.coupling_pitch:.3f} µm, far pitch {dc.waveguide_offset:.2f} µm)")
print(f"BODY_STUB sweep: {STUB_SWEEP[0]:.2f}..{STUB_SWEEP[-1]:.2f} µm "
      f"({len(STUB_SWEEP)} values), each over {len(coupling_lengths)} coupling "
      f"lengths ({coupling_lengths[0]:.2f}–{coupling_lengths[-1]:.2f} µm) x {N_WL} wavelengths")


# ---------------------------------------------------------------------------
# Per-stub simulation: 2 FDTD section solves + 1 EME length-sweep solve,
# combined via the same matrix-chain contraction as hybrid.py's step 5.
# ---------------------------------------------------------------------------
def simulate_stub(stub: float) -> None:
    stub = float(stub)
    boundaries = [
        -half_total,          # physical input port (access waveguide end)
        -half_body + stub,    # input section <-> body section interface (stub past the seam)
        +half_body - stub,    # body section <-> output section interface (stub past the seam)
        +half_total,          # physical output port (access waveguide end)
    ]

    # Length audit, same as hybrid.py's: each FDTD section folds `stub` µm of
    # the straight coupling region into itself, so what the three sections
    # simulate must add back up to the nominal coupling_length, exactly.
    body_lengths = coupling_lengths - 2 * stub
    assert body_lengths[0] >= 0, (
        f"coupling_lengths must start at >= 2*stub ({2 * stub:.2f} µm); "
        f"got {coupling_lengths[0]:.2f} µm"
    )
    assert np.allclose(body_lengths + 2 * stub, coupling_lengths)

    print(f"[sweep] BODY_STUB = {stub:.2f} µm: FDTD sections hold 2 x {stub:.2f} µm of "
          f"straight, EME body {body_lengths[0]:.2f}–{body_lengths[-1]:.2f} µm, "
          f"total {coupling_lengths[0]:.2f}–{coupling_lengths[-1]:.2f} µm; "
          f"interface monitors at x = ±{abs(boundaries[1]):.2f} µm "
          f"(seam at ±{half_body:.2f} µm), FDTD section span "
          f"{boundaries[1] - boundaries[0]:.2f} µm")

    M_in0 = fdtd_solver.solve_section(
        dc, wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
        input_idx=0, save_dir=save_dir, **_iface_mon,
    )  # [n_wl, N_MODES, 1] — "in" arm -> body modes
    M_in1 = fdtd_solver.solve_section(
        dc, wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
        input_idx=1, save_dir=save_dir, **_iface_mon,
    )  # [n_wl, N_MODES, 1] — "isolated" arm -> body modes

    M_out_through = M_in0.transpose(0, 2, 1)  # mirror-symmetry reuse, not a new simulation
    M_out_cross = M_in1.transpose(0, 2, 1)

    M_body = body_solver.solve_section_length_sweep(
        dc, wavelengths, boundaries[1], boundaries[2], body_lengths,
        n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
        neff_override=_body_neff,   # None unless DC_CALIBRATE_BODY=1
    )  # [n_L, n_wl, N_MODES, N_MODES] — plotted/saved against coupling_lengths (the physical total)

    tmp = np.einsum("Lwij,wjk->Lwik", M_body, M_in0)               # [n_L, n_wl, N_MODES, 1]
    full_through = np.einsum("wij,Lwjk->Lwik", M_out_through, tmp) # [n_L, n_wl, 1, 1]
    full_cross = np.einsum("wij,Lwjk->Lwik", M_out_cross, tmp)     # [n_L, n_wl, 1, 1]

    T_through_grid = np.abs(full_through[:, :, 0, 0]) ** 2   # [n_L, n_wl]
    T_cross_grid = np.abs(full_cross[:, :, 0, 0]) ** 2        # [n_L, n_wl]
    T_total_grid = T_through_grid + T_cross_grid

    L_50 = split_point(T_cross_grid, T_total_grid)

    # Atomic, because aggregate() below reads *every* stub's file — including
    # ones a concurrent array task is writing this instant. A plain np.savez
    # here is what made 7 of 9 tasks in job 11517245 die on BadZipFile; the
    # uncalibrated run survived only because its tasks finished 12 minutes
    # apart instead of 6. Open-file-object form: np.savez appends '.npz' to a
    # suffix-less path, which would defeat the rename (dc_validation.md item 6).
    def _write_task_npz(tmp: str) -> None:
        # `with`, not a bare open(): zipfile does NOT close a file object
        # it was handed, so relying on refcounting to flush would let
        # os.replace move a truncated file — the very bug being fixed.
        with open(tmp, "wb") as fh:
            np.savez(
                fh,
                body_stub=stub, coupling_lengths=coupling_lengths, wavelengths=wavelengths,
                T_through=T_through_grid, T_cross=T_cross_grid, L_50=L_50,
            )

    atomic_replace_write(
        os.path.join(data_dir, f"{_tag(stub)}.npz"), _write_task_npz)
    print(f"[sweep] BODY_STUB={stub:.2f} µm: 50:50 split at L = {L_50:.3f} µm")

    # -- Per-stub sweep plot (same layout as hybrid.py's coupling_length_sweep.png) --
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for grid, label, color in [
        (T_through_grid, "Through", "steelblue"),
        (T_cross_grid,   "Cross",   "darkorange"),
        (T_total_grid,   "Total",   "k"),
    ]:
        mean, lo, hi = grid.mean(axis=1), grid.min(axis=1), grid.max(axis=1)
        ax.plot(coupling_lengths, mean, label=label, color=color,
                lw=1.5 if label == "Total" else 1.2,
                linestyle="--" if label == "Cross" else "-")
        ax.fill_between(coupling_lengths, lo, hi, color=color, alpha=0.15, lw=0)
    if np.isfinite(L_50):
        ax.axvline(L_50, color="gray", linestyle=":", lw=0.8,
                   label=f"50:50 split ({L_50:.2f} µm)")
    ax.set_xlabel("Coupling length (µm)")
    ax.set_ylabel("Transmission (mean over wavelength, shaded = min–max)")
    ax.set_title(f"Power Transfer vs Coupling Length (BODY_STUB={stub:.2f} µm)")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(plot_dir, f"coupling_length_sweep_{_tag(stub)}.png"), dpi=300)
    plt.close("all")
    print(f"[sweep] Plot saved → {plot_dir}/coupling_length_sweep_{_tag(stub)}.png")


# ---------------------------------------------------------------------------
# Aggregation: collect every completed stub into the headline 50:50-split-point
# plot + stacked grid, mirroring the other sweeps' atomic-write pattern.
# ---------------------------------------------------------------------------
def aggregate() -> None:
    done, L50s, T_through_rows, T_cross_rows = [], [], [], []
    for stub in STUB_SWEEP:
        path = os.path.join(data_dir, f"{_tag(stub)}.npz")
        if not os.path.exists(path):
            continue
        d = np.load(path)
        done.append(float(d["body_stub"]))
        L50s.append(float(d["L_50"]))
        T_through_rows.append(d["T_through"])
        T_cross_rows.append(d["T_cross"])

    n_done, n_all = len(done), len(STUB_SWEEP)
    if n_done < 2:
        print(f"[sweep] Aggregation skipped — only {n_done}/{n_all} stubs done.")
        return

    done = np.asarray(done)
    L50s = np.asarray(L50s)
    T_through_stack = np.stack(T_through_rows)   # [n_stub, n_L, n_wl]
    T_cross_stack = np.stack(T_cross_rows)       # [n_stub, n_L, n_wl]

    def _write_npz(tmp):
        with open(tmp, "wb") as f:
            np.savez(f, body_stubs=done, L_50=L50s, coupling_lengths=coupling_lengths,
                     wavelengths=wavelengths, T_through=T_through_stack, T_cross=T_cross_stack)

    atomic_replace_write(os.path.join(output_dir, "stub_sweep.npz"), _write_npz)
    print(f"[sweep] Stub-sweep grid saved ({n_done}/{n_all} stubs) → {output_dir}/stub_sweep.npz")

    # Single-domain ground truth, when the full-device FDTD length sweep has
    # been run: its own 50:50 length, extracted the same way.
    L_50_fdtd = float("nan")
    if os.path.exists(_FDTD_SWEEP_NPZ):
        f = np.load(_FDTD_SWEEP_NPZ)
        ratio = f["T_cross"].mean(axis=1) / np.clip(
            (f["T_cross"] + f["T_through"]).mean(axis=1), 1e-20, None)
        cross = np.flatnonzero(np.diff(np.sign(ratio - 0.5)))
        if cross.size:
            i = int(cross[0])
            Lf = f["coupling_lengths"]
            r0, r1 = ratio[i], ratio[i + 1]
            L_50_fdtd = float(Lf[i] + (0.5 - r0) * (Lf[i + 1] - Lf[i]) / (r1 - r0)) \
                if r1 != r0 else float(Lf[i])

    # constrained_layout (not tight_layout) is what correctly reserves space
    # for the colorbar and suptitle together.
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)

    # -- Panel 1: the headline result — 50:50 length vs stub ----------------
    ax = axes[0]
    ax.plot(done, L50s, color="seagreen", marker="o", markersize=5, lw=1.4,
            label="Hybrid (FDTD sections + EME body)")
    if np.isfinite(L_50_fdtd):
        ax.axhline(L_50_fdtd, color="crimson", linestyle="--", lw=1.0,
                   label=f"Full-device FDTD ({L_50_fdtd:.2f} µm)")
    if n_done >= 2:
        # How much the answer still moves per 0.25 µm of extra stub, over the
        # last two completed stubs — the convergence statement in one number.
        drift = abs(L50s[-1] - L50s[-2])
        # Lower left: the curve falls left-to-right, so that corner stays clear
        # (upper right belongs to the legend).
        ax.annotate(f"last step: {drift:.3f} µm", xy=(0.02, 0.04), xycoords="axes fraction",
                    ha="left", va="bottom", fontsize=8, color="#555555")
    ax.set_xlabel("BODY_STUB — straight coupling region inside each FDTD section (µm)")
    ax.set_ylabel("50:50 split coupling length (µm)")
    ax.set_title("3 dB Coupling Length vs FDTD/EME Interface Placement")
    ax.set_xticks(STUB_SWEEP)
    ax.legend(fontsize=8)
    ax.grid(True)

    # -- Panel 2: the S-curve family the split points come from -------------
    ax = axes[1]
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=STUB_SWEEP.min(), vmax=STUB_SWEEP.max())
    for i, stub in enumerate(done):
        ax.plot(coupling_lengths, T_cross_stack[i].mean(axis=1), color=cmap(norm(stub)), lw=1.2)
    ax.axhline(0.5, color="gray", linestyle=":", lw=0.8, label="50:50")
    ax.set_xlabel("Coupling length (µm)")
    ax.set_ylabel("Cross transmission (mean over wavelength)")
    ax.set_title("S-curve shape vs BODY_STUB")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="BODY_STUB (µm)")

    if n_done < n_all:
        fig.text(0.02, 0.97, f"{n_done}/{n_all} stubs done", ha="left", va="top",
                 fontsize=8, color="#555555",
                 bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"))

    fig.suptitle("Directional Coupler — Hybrid Body-Stub (Interface Placement) Study", fontsize=11)
    atomic_replace_write(
        os.path.join(output_dir, "stub_sweep_split_point.png"),
        lambda tmp: fig.savefig(tmp, dpi=300, format="png"),
    )
    plt.close("all")
    print(f"[sweep] Split-point plot saved ({n_done}/{n_all} stubs) → "
          f"{output_dir}/stub_sweep_split_point.png")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        i = int(sys.argv[1])
        if not 0 <= i < len(STUB_SWEEP):
            raise SystemExit(
                f"index {i} out of range — sweep has {len(STUB_SWEEP)} "
                f"BODY_STUB values (0..{len(STUB_SWEEP) - 1})")
        simulate_stub(STUB_SWEEP[i])
    else:
        print(f"[sweep] No index given — running all {len(STUB_SWEEP)} "
              "BODY_STUB values serially (the SLURM array batch script is much faster).")
        for stub in STUB_SWEEP:
            simulate_stub(stub)
    aggregate()
