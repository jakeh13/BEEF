"""Tutorial: 4x4 MMI 90-degree optical hybrid -- full single-domain FDTD
reference simulation (BEEF method).

validation_markdown/mmi_4x4_90_hybrid_validation.md Step 6b (check 4's own
full-device reference). Same device as mmi_4x4_90hybrid_hybrid.py/
mmi_4x4_90hybrid_eme.py -- see mmi_4x4_90hybrid_hybrid.py's own module
docstring for the paper reference and design caveats. RE-DERIVED 2026-08-14
to the current settled geometry (ΔL=-2.75, mm_length=78.9396, was -1.50/
80.1896) -- the point this script originally shipped with predates the
2026-08-09 N_MODES_BODY=8->16 body-basis convergence fix that moved ΔL, so
job 11733865 (below) ran the STALE point and its own timeout is superseded,
not just re-attempted at the same geometry. taper_width/taper_length remain
unmeasured placeholders (Step 2), unchanged from the hybrid/EME scripts.

Previous attempt (job 11733865, 2026-08-08): TIMED OUT at the full -t 360
(6h) budget on both array tasks -- not a crash, both were still stepping
cleanly (~0.29s/step, no NaNs/slowdown) at meep time~876, genuinely still
converging toward decay_tol=1e-11. Paused at the user's direction rather
than resubmitted at a bigger budget at the time (validation doc, Step 6b).
Re-queued 2026-08-14 at the user's direction, at the corrected geometry and
a substantially larger time budget (see the batch script) -- an ~84 um,
4-port single-domain box storing far more energy than any smaller section
this repo has run to full decay is a real, uncertain cost, not assumed
solved by the small ΔL correction alone.

Method: single-domain FDTD via meep, the whole ~84 um device in one box --
the most physically complete of the three simulation methods (unlike the
EME script's per-port-pair decoupling approximation, or even the hybrid
script's own section-by-section cascade, this one has every port's real
geometry present simultaneously in a single run).

A real library gap, found and fixed here, not worked around:
FDTDSolver.solve()'s own `port` parameter only ever accepted the strings
"top"/"bot" (geometry.input_ys' LAST or FIRST entry) -- fine for every
existing 2-port MMI/coupler example in this repo, but this device needs
ports 1 AND 3 (paper's own I={1,3}, finding A4), and port 3 (0-indexed 2 of
4) is an INTERIOR port neither "top" nor "bot" can reach at all. Rather
than restrict this script to only testing port 1 (which would silently
skip half of what check 4 needs to verify -- including A8's own phase-
consistency requirement, which specifically needs BOTH excitations),
`FDTDSolver.solve()` gained an explicit `input_idx` parameter that
overrides the port-string selection for any 0-based index into
`geometry.input_ys()`. Verified: an incidental related bug fixed at the
same time -- previously, when a caller didn't pass `input_idx` via
**kwargs (every existing example script, including this repo's own
mmi_2x2_fdtd.py), the source APERTURE width was computed for
geometry.source_aperture_at()'s own DEFAULT input_idx=0 regardless of
which port `port="top"/"bot"` actually excited -- invisible for the 2x2
MMI specifically only because its aperture is symmetric between its two
ports, but a real latent mismatch for any asymmetric device. Now always
resolved and passed explicitly. Full fast suite re-run clean (639 passed,
19 skipped) after the change; not independently re-verified with a real
meep run on an existing device (the symmetric-aperture argument above is
what makes that safe to defer, not a substitute for eventually doing it).

Output monitors: one per output port, width=waveguide_offset=1.5 um
(touching neighbor midpoints exactly, matching the same offset spacing
that keeps ports from overlapping in the first place -- Step 2).

Cost, flagged up front per the validation doc's own Step 6b entry: this is
the longest, most complex single-domain FDTD box built for this device so
far (~84 um long, 4 ports x 2 faces). RESOLUTION/DECAY_TOL/N_WL are env-var
overridable (MMI4X4_FDTD_RES / MMI4X4_FDTD_DECAY_TOL / MMI4X4_FDTD_NWL) for
a smoke pass -- run one BEFORE the real 2-port sweep, per this repo's own
two-step discipline for expensive, never-run scripts.

MMI4X4_MM_LENGTH overrides the production mm_length (159.30, the campaign
optimum). At the optimum the full domain is ~167 um long and was estimated
at 100-140 h/port -- too expensive to run. A short body (<20 um) shrinks
the domain ~6-7x and is not claiming to be a device operating point; it
exists to get real engine-agreement (Check 1) data at this device's own
cross-section, at a length actually affordable to run to full decay, and to
see whether the hybrid-vs-full-FDTD systematic gap found at the old
geometry (job 12272573: hybrid under-collects power, worst at short
wavelength, up to +30-50%) is length-dependent. Output dir gets a
`_L{mm_length}` suffix whenever this override is used, so it never
collides with the production 159.30 um run's own files.

Run:
    python /path/to/examples/mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_fdtd.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, FDTDSolver, SMatrix, SOI_220, plot_fdtd_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters -- env-var overridable for a smoke pass.
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI4X4_FDTD_RES", 40.0))
DECAY_TOL = float(os.environ.get("MMI4X4_FDTD_DECAY_TOL", 1e-4))
# BAND: 1.31 +- 0.02 um -- the band the geometry above was OPTIMISED over.
# These scripts previously ran 1.26-1.36 (+-50 nm), which is 2.5x wider than
# the design band, and this device is steeply wavelength-dependent: measured
# 2026-09-05, min CMRR falls 20.761 -> 18.482 dBe and excess loss rises
# 0.654 -> 2.706 dB going from 1.310 um alone to +-20 nm (job 12831929).
# Reporting the optimised geometry over +-50 nm would judge it far outside the
# band it was designed for. The old window is still reachable via the env
# overrides for anyone who wants the wider sweep.
WL_MIN = float(os.environ.get("MMI4X4_WLMIN", 1.29))
WL_MAX = float(os.environ.get("MMI4X4_WLMAX", 1.33))
N_WL = int(os.environ.get("MMI4X4_FDTD_NWL", 3))   # matches the paper's own Eq. 12 convention

# ---------------------------------------------------------------------------
# 1. Geometry -- identical to the hybrid/EME scripts' settled parameters.
# ---------------------------------------------------------------------------
# Geometry -- the OPTIMISED design point, not Step 2's first pass.
# Result of the (mm_width, y_gap, taper_width, taper_length, mm_length)
# campaign documented in validation_markdown/mmi_4x4_90_hybrid_validation.md
# (~45 hybrid-path evaluations, rounds 8-10, 2026-08-28..09-05). Measured at
# this point over 1.29-1.33 um: max|dphi| = 3.933 deg (spec <=5, MET),
# min CMRR = 18.482 dBe (spec >=20, 1.52 short), excess loss 2.706 dB,
# imbalance 1.972 dB. Superseded first pass, for the record: mm_width=6.0,
# taper_width=0.8, taper_length=1.0, mm_length=78.9396 -- which measured
# 4.764 deg / 12.109 dBe / 4.076 dB / 5.202 dB on the same evaluator.
#
# taper_width tracks mm_width: the campaign's most reliable finding is that a
# taper which has not widened with the body cannot be rescued by ANY length
# (0 of 7601 lengths met the phase spec at mm_width=8.5 with taper 1.20).
# Measured optima: 6.5 -> 1.065, 7.0 -> 1.10, 7.5 -> 1.20, 8.0 -> 1.20,
# 8.5 -> 1.30.
WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 8.5
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 1.30            # MEASURED optimum at this width (round 10, job 12783704)
TAPER_LENGTH = 4.0            # MEASURED (round 8): tl=1->4 cut excess loss 4.076->2.940 dB;
                              # tl=8 bought only 0.006 more objective for 2x the footprint
L_PI = 218.5954   # = 108.9195 * (8.5/6.0)**2 -- SCALED, not re-measured. L_pi/W^2 was
                  # measured constant to 0.1% over 3.5-16 um (Step 1c), so the mm_width=6.0
                  # measurement (job 11693043) carries.

# DELTA_L: re-derived 2026-08-09 at the converged N_MODES_BODY=16 body basis
# (job 11768847) -- was -1.50 (mm_length=80.1896), found at the since-
# superseded under-converged N=8 basis; see mmi_4x4_90hybrid_hybrid.py's own
# comment for the full story. Asserted against that sweep's exact grid,
# same discipline mmi_1x4_hybrid.py's own grid-snap lesson established.
# MM_LENGTH: the campaign's own optimum for THIS geometry, selected from the
# shared absolute 20-400 um / 0.05 um grid (mmi_4x4_90hybrid_ygap_grid.py's
# MM_LENGTH_GRID), not from a +-10 um DELTA_L window around 3/4*L_pi. The
# length re-score (job 12832166) confirmed this is simultaneously the best
# length under all three selection rules tried -- the live objective,
# phase-as-constraint, and maximise-CMRR-subject-to-phase -- so nothing is
# being traded away by pinning it. Only 19 of 7601 lengths meet the 5 deg
# phase spec at all, so the grid-snap assert below matters more here than it
# did at the first-pass point.
_MM_LENGTH_DEFAULT = 159.30
MM_LENGTH = float(os.environ.get("MMI4X4_MM_LENGTH", _MM_LENGTH_DEFAULT))
_LENGTH_STEP = 0.05
assert abs(round(MM_LENGTH / _LENGTH_STEP) * _LENGTH_STEP - MM_LENGTH) < 1e-9, (
    f"MM_LENGTH={MM_LENGTH} is not an exact multiple of {_LENGTH_STEP} um -- "
    f"avoid silently snapping to an unintended grid point."
)
DELTA_L = MM_LENGTH - 0.75 * L_PI   # = -4.6465 at the default; informational
                                     # only, and meaningless (huge negative)
                                     # at an engine-agreement short length --
                                     # self-imaging is not being targeted there.

# MMI4X4_DRIVEN_PORTS parsed early (not just at section 3) so the output
# directory can be tagged by it too -- a top-port-only run must never
# collide with the canonical ports-1&3 run's files at the same mm_length.
_PORT_TO_INPUT_IDX = {1: 0, 2: 1, 3: 2, 4: 3}
_driven_ports = [int(p) for p in
                  os.environ.get("MMI4X4_DRIVEN_PORTS", "1,3").split(",")]
_ports_tag = "" if set(_driven_ports) == {1, 3} else (
    "_p" + "".join(str(p) for p in _driven_ports))

_mm_length_tag = "" if abs(MM_LENGTH - _MM_LENGTH_DEFAULT) < 1e-9 else f"_L{MM_LENGTH:.2f}"

# 2026-09-11: decay_tol must tag the output directory for the same reason
# mm_length and the driven-port set do. MMI4X4_FDTD_DECAY_TOL was already
# overridable, but the results path was NOT tagged by it, so a convergence
# probe at a tighter tolerance would silently overwrite the production
# reference results at the same (mm_length, ports) -- exactly the near-miss
# this device already hit once when coarse smoke runs landed in the official
# output_dir. Empty at the 1e-4 default, so no existing path changes.
_decay_tol_tag = "" if abs(DECAY_TOL - 1e-4) < 1e-12 else f"_tol{DECAY_TOL:.0e}"
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "mmi_4x4_90hybrid_fdtd" + _mm_length_tag + _ports_tag + _decay_tol_tag))
os.makedirs(output_dir, exist_ok=True)

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
    use_neighbor_aperture=True,
)
print(f"MMI total length: {mmi.total_length:.4f} µm (mm_length={mmi.mm_length:.4f})")
print(f"Input ports at y = {mmi.input_ys} µm; outputs at y = {mmi.output_ys} µm")

# ---------------------------------------------------------------------------
# 2. FDTD solver -- one output monitor per port, width=waveguide_offset so
#    adjacent monitors touch at the midpoint without overlapping.
# ---------------------------------------------------------------------------
solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=MM_WIDTH + 3.0,
    sim_height=2.0,             # Step 1b's converged sz
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    num_input_modes=1,
    num_output_modes=1,
    output_monitor_centers=list(mmi.output_ys),
    output_monitor_widths=[WAVEGUIDE_OFFSET] * 4,
    material_stack=SOI_220,
    decay_tol=DECAY_TOL,   # 1e-11 default was MEASURED UNREACHABLE -- see
               # mmi_1x2_decay_tol_measurement.md; 1e-4 + energy decay is
               # what every campaign point used, sufficient at 0.015 dB
    decay_method="energy",
    energy_decay_dt=20.0,
)

plot_fdtd_geometry(
    mmi, solver=solver,
    title="4×4 MMI 90° hybrid — Full-Device FDTD — Top-down Geometry",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved → {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 3. Sweep -- driven inputs. Default: paper's ports 1/3 (A4). Port 1
#    (input_ys[0], the lowest-y entry) is reachable via the existing
#    port="bot" string; port 3 (input_ys[2], an INTERIOR port of 4) needs
#    the input_idx parameter -- "top"/"bot" cannot reach it at all.
#
# MMI4X4_DRIVEN_PORTS overrides the driven set -- a comma list of 1-based
# paper port numbers (e.g. "4" for just the top port, input_ys[3]). Default
# "1,3" is unchanged from every prior run of this script. When the driven
# set is anything other than exactly {1,3}, the paper's own Eqs. 10-15
# metrics (CMRR, Δφ) are skipped -- they need BOTH ports 1 and 3's
# cross-excitation phase relationship (A8) to mean anything -- and this
# prints/saves plain per-arm transmission + phase for whatever WAS driven
# instead, e.g. for a same-monitor-locations top-port comparison against the
# hybrid path at a matched mm_length (validation doc Step 8, Check 1
# follow-up).
#
# MMI4X4_FDTD_CAPTURE_FIELDS=1 additionally captures one 2-D z=0 field slice
# per driven port (solve()'s new capture_fields_dir parameter -- see
# beef/solvers/fdtd.py) to output_dir/fields/, for a full-device top-down
# field picture at this mm_length.
# ---------------------------------------------------------------------------
CAPTURE_FIELDS = os.environ.get("MMI4X4_FDTD_CAPTURE_FIELDS", "0") == "1"
_fields_dir = os.path.join(output_dir, "fields") if CAPTURE_FIELDS else None

# ---------------------------------------------------------------------------
# MMI4X4_FDTD_SNAP_EVERY=<t> decomposes the monitors every <t> meep time units
# and streams {t, energy, energy_ratio, |T|^2, arg(S)} to a CSV in output_dir.
#
# WHY. `decay_tol=1e-4` was selected on a |T|^2 AMPLITUDE criterion (0.015 dB,
# 2026-08-25) which says NOTHING about phase, and job 13110659 then bounded the
# reference's phase noise from exactly TWO points (1e-4 vs 1e-5) -- a delta,
# not a convergence test. Because DFT monitors accumulate, one snapshotted run
# yields the whole phase-vs-decay_tol curve, so the ladder of runs that would
# otherwise be needed collapses to one. beef/solvers/fdtd.py's writer gained
# the `ph_out{i}_wl{j}` columns for this.
#
# A TIMEOUT IS AN EXPECTED OUTCOME HERE, NOT A FAILURE. To see past the
# tolerance the run would otherwise stop at, this is meant to be launched with
# DECAY_TOL set BELOW where the run can reach, so it walls out -- and the
# energy_ratio column then also shows where (or whether) this domain's decay
# FLOORS, which is unmeasured at L=15 and is the thing that decides whether a
# tolerance tighter than 1e-5 is reachable at all. Samples are appended as
# they are taken, so a walled run keeps every sample it wrote.
#
# NOT part of the cache fingerprint (a method argument, deliberately not stored
# on self -- see solve()'s docstring), so this cannot strand or poison an
# existing cache entry. But note the converse: a run whose (geometry, res,
# decay_tol) already has a cache entry returns WITHOUT time-stepping and writes
# NO curve. Vary decay_tol to force a fresh solve.
# ---------------------------------------------------------------------------
_snap_every = os.environ.get("MMI4X4_FDTD_SNAP_EVERY", "")
SNAP_EVERY = float(_snap_every) if _snap_every else None
_snap_path = os.path.join(
    output_dir, f"tolcurve_p{''.join(str(p) for p in _driven_ports)}.csv"
) if SNAP_EVERY else None

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
PORT_NAMES = ["out5", "out6", "out7", "out8"]

print(f"\nSweeping {len(wavelengths)} wavelengths: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm")
print("(single broadband pulse per input port; DFT monitors extract all wavelengths in one run)")
print(f"Driven ports (1-based, paper numbering): {_driven_ports}"
      f"{'  [capturing fields]' if CAPTURE_FIELDS else ''}")

smat = {}
for p in _driven_ports:
    idx = _PORT_TO_INPUT_IDX[p]
    # One CSV per driven port: the snapshot writer APPENDS, so a shared path
    # across ports would interleave two runs' samples into one unreadable file.
    _sp = (_snap_path.replace(".csv", f"_port{p}.csv")
           if _snap_path else None)
    smat_raw = solver.solve(
        mmi, wavelengths, input_idx=idx, input_mode=1, save_dir=save_dir,
        capture_fields_dir=_fields_dir,
        snapshot_every=SNAP_EVERY, snapshot_path=_sp,
    )
    smat_raw.check_passive()
    s = SMatrix(wavelengths=wavelengths)
    for j in range(4):
        s.set(smat_raw.coeff(port=f"out_{j}", in_mode=1, out_mode=1), port=PORT_NAMES[j])
    smat[idx] = s
    print(f"  input port {p} (input_idx={idx}) solved")

# ---------------------------------------------------------------------------
# 4. Metrics (paper Eqs. 10-15), same SMatrix methods as the hybrid/EME
#    scripts. A11's phase-sign fix applies identically here. Only meaningful
#    (and only computed) when exactly the paper's own I={1,3} were driven.
# ---------------------------------------------------------------------------
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}
DRIVEN = [0, 2]

if set(_driven_ports) == {1, 3}:
    imbalance = {i: smat[i].imbalance_db(PORT_NAMES) for i in DRIVEN}
    excess_loss = {i: smat[i].excess_loss_db(PORT_NAMES) for i in DRIVEN}
    cmrr_I = {i: smat[i].cmrr_db("out5", "out8") for i in DRIVEN}
    cmrr_Q = {i: smat[i].cmrr_db("out6", "out7") for i in DRIVEN}

    phi = {j: smat[2].cross_excitation_phase_diff_deg(smat[0], port=PORT_NAMES[j]) for j in range(4)}
    dphi = {}
    for j in [1, 2, 3]:
        raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
        dphi[j] = ((raw + 180.0) % 360.0) - 180.0

    print(f"\nResults at mm_length={MM_LENGTH:.4f} µm (full-device FDTD reference):")
    for i, port_label in zip(DRIVEN, ["1", "3"]):
        print(f"  Input port {port_label}: imbalance={imbalance[i]} dB  "
              f"excess_loss={excess_loss[i]} dB  CMRR_I={cmrr_I[i]} dBe  CMRR_Q={cmrr_Q[i]} dBe")
    print(f"  Phase error Δφ_6/7/8 = {dphi[1]} / {dphi[2]} / {dphi[3]} deg")
else:
    # Plain per-arm |T|^2 and phase(S) for whatever was driven -- no
    # cross-excitation metric needs a second port.
    print(f"\nResults at mm_length={MM_LENGTH:.4f} µm (full-device FDTD reference), "
          f"driven ports {_driven_ports} (raw per-arm, no CMRR/Δφ -- needs ports 1&3):")
    for p in _driven_ports:
        idx = _PORT_TO_INPUT_IDX[p]
        for j, name in enumerate(PORT_NAMES):
            T = np.abs(smat[idx].coeff(name)) ** 2
            phase_deg = np.angle(smat[idx].coeff(name), deg=True)
            print(f"  port {p} -> {name}: |T|^2={T} phase={phase_deg} deg")

# ---------------------------------------------------------------------------
# 5. Save. Filenames are port-labeled by paper number (1-based) so an
#    arbitrary driven set (e.g. just port 4/top) never collides with the
#    canonical ports-1&3 files.
# ---------------------------------------------------------------------------
for p in _driven_ports:
    smat[_PORT_TO_INPUT_IDX[p]].save(os.path.join(output_dir, f"smatrix_fdtd_port{p}.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix_fdtd_port{{{','.join(str(p) for p in _driven_ports)}}}.pkl")

if set(_driven_ports) == {1, 3}:
    np.savez(os.path.join(output_dir, "fdtd_results.npz"),
             wavelengths=wavelengths, mm_length=MM_LENGTH,
             imbalance_in1=imbalance[0], imbalance_in3=imbalance[2],
             excess_loss_in1=excess_loss[0], excess_loss_in3=excess_loss[2],
             cmrr_I_in1=cmrr_I[0], cmrr_Q_in1=cmrr_Q[0],
             cmrr_I_in3=cmrr_I[2], cmrr_Q_in3=cmrr_Q[2],
             dphi_6=dphi[1], dphi_7=dphi[2], dphi_8=dphi[3],
             phi_5=phi[0], phi_6=phi[1], phi_7=phi[2], phi_8=phi[3])
    print(f"Results saved → {output_dir}/fdtd_results.npz")
else:
    np.savez(os.path.join(output_dir, "fdtd_results.npz"),
             wavelengths=wavelengths, mm_length=MM_LENGTH,
             driven_ports=np.array(_driven_ports),
             **{f"T_port{p}_{name}": np.abs(smat[_PORT_TO_INPUT_IDX[p]].coeff(name)) ** 2
                for p in _driven_ports for name in PORT_NAMES},
             **{f"phase_port{p}_{name}": np.angle(smat[_PORT_TO_INPUT_IDX[p]].coeff(name), deg=True)
                for p in _driven_ports for name in PORT_NAMES})
    print(f"Results saved → {output_dir}/fdtd_results.npz")

# ---------------------------------------------------------------------------
# 6. Plot. Ports-1&3 case: same 2x2 layout as the hybrid/EME scripts' own
#    results plots (CMRR, phase error, excess loss, quadrature phase), so
#    all three simulation methods' results are directly, visually
#    comparable. Otherwise: plain per-arm |T|^2 and phase(S) vs wavelength,
#    one line per (driven port, output arm) pair.
# ---------------------------------------------------------------------------
if set(_driven_ports) == {1, 3}:
    WRAPPED_TARGET_DEG = {j: ((v + 180.0) % 360.0) - 180.0 for j, v in TARGET_DEG.items()}

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(wavelengths, cmrr_I[0], "o-", label="CMRR_I (port 1)", color="steelblue")
    ax.plot(wavelengths, cmrr_Q[0], "o-", label="CMRR_Q (port 1)", color="darkorange")
    ax.plot(wavelengths, cmrr_I[2], "s--", label="CMRR_I (port 3)", color="steelblue", alpha=0.5)
    ax.plot(wavelengths, cmrr_Q[2], "s--", label="CMRR_Q (port 3)", color="darkorange", alpha=0.5)
    ax.axhline(20.0, color="red", linestyle=":", lw=1.0, label="cited requirement (20 dBe)")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("CMRR (dBe)")
    ax.set_title("CMRR vs Wavelength")
    ax.legend(fontsize=7)
    ax.grid(True)

    ax = axes[0, 1]
    ax.plot(wavelengths, dphi[1], "o-", label="Δφ_6", color="steelblue")
    ax.plot(wavelengths, dphi[2], "o-", label="Δφ_7", color="darkorange")
    ax.plot(wavelengths, dphi[3], "o-", label="Δφ_8", color="seagreen")
    ax.axhline(5.0, color="red", linestyle=":", lw=1.0, label="cited requirement (5°)")
    ax.axhline(-5.0, color="red", linestyle=":", lw=1.0)
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Phase error Δφ_j (deg)")
    ax.set_title("Phase Error vs Wavelength")
    ax.legend(fontsize=7)
    ax.grid(True)

    ax = axes[1, 0]
    ax.plot(wavelengths, excess_loss[0], "o-", label="Input port 1", color="steelblue")
    ax.plot(wavelengths, excess_loss[2], "s--", label="Input port 3", color="darkorange")
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Excess loss (dB)")
    ax.set_title("Excess Loss vs Wavelength (Eq. 11)")
    ax.legend(fontsize=7)
    ax.grid(True)

    ax = axes[1, 1]
    _phi_labels = ["φ_5 (ref)", "φ_6", "φ_7", "φ_8"]
    _phi_colors = ["gray", "steelblue", "darkorange", "seagreen"]
    for j in range(4):
        ax.plot(wavelengths, phi[j], "o-", label=_phi_labels[j], color=_phi_colors[j])
        ax.axhline(WRAPPED_TARGET_DEG[j], color=_phi_colors[j], linestyle=":", lw=1.0, alpha=0.6)
    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Quadrature phase φ_j (deg)")
    ax.set_title("Quadrature Phase vs Wavelength\n(dotted = paper's target angle, Eq. 15)")
    ax.legend(fontsize=7)
    ax.grid(True)

    plt.suptitle(f"4×4 MMI 90° hybrid — full-device FDTD — mm_length={MM_LENGTH:.2f} µm", fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "fdtd_results.png"), dpi=300)
    plt.close("all")
    print(f"Results plot saved → {output_dir}/fdtd_results.png")
else:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    _colors = {"out5": "steelblue", "out6": "darkorange", "out7": "seagreen", "out8": "crimson"}
    for p in _driven_ports:
        idx = _PORT_TO_INPUT_IDX[p]
        for name in PORT_NAMES:
            T = np.abs(smat[idx].coeff(name)) ** 2
            ph = np.angle(smat[idx].coeff(name), deg=True)
            axes[0].plot(wavelengths, T, "o-", label=f"port{p}→{name}", color=_colors[name])
            axes[1].plot(wavelengths, ph, "o-", label=f"port{p}→{name}", color=_colors[name])
    axes[0].set_xlabel("Wavelength (µm)"); axes[0].set_ylabel("|T|^2")
    axes[0].set_title("Per-arm transmission"); axes[0].legend(fontsize=7); axes[0].grid(True)
    axes[1].set_xlabel("Wavelength (µm)"); axes[1].set_ylabel("phase(S) (deg)")
    axes[1].set_title("Per-arm phase"); axes[1].legend(fontsize=7); axes[1].grid(True)
    plt.suptitle(f"4×4 MMI 90° hybrid — full-device FDTD — mm_length={MM_LENGTH:.2f} µm — "
                 f"driven ports {_driven_ports}", fontsize=10)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "fdtd_results.png"), dpi=300)
    plt.close("all")
    print(f"Results plot saved → {output_dir}/fdtd_results.png")
