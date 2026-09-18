"""
Tutorial: simulate the WAISTED 2x2 MMI splitter with the BEEF library using
whole-device FDTD.

Method: one meep simulation of the complete device — no EME, no decomposition,
no reference planes. This is the reference the other two scripts in this folder
(mmi_2x2_eme.py, mmi_2x2_hybrid.py) are checked against, and the slowest of the
three: hours, against seconds to minutes for the pure-EME script.

Geometry: the tapered-body variant of examples/mmi_2x2. Same ports, same
platform, same band — the only difference is that the multimode body is not a
uniform slab. It narrows from the 2.8 µm face to a 2.24 µm interior waist over
a 5 µm transition on each side:

    |<-- s -->|<-- Lt -->|<---- Lw ---->|<-- Lt -->|<-- s -->|
     2.80 um    ramp        2.24 um        ramp     2.80 um
                              waist

Narrowing the body raises the modal dispersion, which shortens the first
self-imaging length. The device below reaches the same 50/50 split with 90°
quadrature in 31.00 µm total against the straight control's 37.00 µm —
16.2 % shorter.

VALIDATION STATUS (validation_markdown/resolved/mmi_2x2_waist_body_validation.md).
The device is VALIDATED as of 2026-09-14, and this script IS the check that
validated it: single-domain, whole-device FDTD is the only measurement that
depends on neither the EME cascade nor any calibration. At the design point it
measures IL 0.2010 dB, imbalance 0.0124, phase mean-dev 0.267° and MAX-dev
0.447°, min total T 0.9548 — passing all four criteria under BOTH the mean and
the max phase reading.

⚠ TWO RETRACTED DESIGN POINTS. Earlier versions of this study headlined a
29.50 µm (20.3 % shorter) point and then a 30.50 µm one. Both are RETRACTED:
the 29.50 µm point fails on two settled taper bases and on full-device FDTD and
was never hybrid-buildable, and the 30.50 µm point passes the mean phase
reading but fails the max one. Only the geometry below survives both readings.
Do not resurrect the others from an old notebook.

⚠ decay_tol = 1e-05 HERE, NOT the 1e-11 the straight mmi_2x2_fdtd.py inherits.
This is the single most expensive lesson in this device's history — it cost
4 h 15 m of 48-rank time. 1e-11 is BELOW THE REACHABLE FLOOR on a device this
long: the DFT residual floors and then diverges, so the run never converges and
simply burns its wall clock. 1e-05 is the tolerance every stored reference for
this geometry was produced at, so it is also what makes those references
comparable to a fresh run. Verify rather than trust it: snapshot_every below
turns one run into the whole transmission-vs-decay_tol curve at no extra cost,
because DFT monitors accumulate.

Method: FDTD via meep with GDS prism geometry and EigenModeSource, identical in
every solver setting to examples/mmi_2x2/python/mmi_2x2_fdtd.py except
decay_tol. That is deliberate — a changed solver config would make this a
different measurement rather than a comparable one.

Run from any directory:
    python /path/to/examples/mmi_2x2_tapered/python/mmi_2x2_fdtd.py

Cheap smoke run (seconds, NOT a measurement):
    MMI2X2T_RES=15 MMI2X2T_NWL=6 python .../mmi_2x2_fdtd.py
"""

import json
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, FDTDSolver, SOI_220, plot_fdtd_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_2x2_tapered",
    "mmi_2x2_fdtd"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_sim_cache"))


def _env_f(name, default):
    return float(os.environ.get(name, default))


def _env_i(name, default):
    return int(os.environ.get(name, default))


# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by the FDTD solver
WL_MIN = 1.530      # wavelength sweep start (µm)
WL_MAX = 1.565      # wavelength sweep end (µm)
N_WL = 51           # number of wavelength points
# ⚠ see the module docstring: 1e-11 is unreachable on this device.
DECAY_TOL = 1e-05
SNAPSHOT_EVERY = 20.0   # meep time units between decay snapshots

RESOLUTION = _env_f("MMI2X2T_RES", RESOLUTION)
N_WL = _env_i("MMI2X2T_NWL", N_WL)

# ---------------------------------------------------------------------------
# 1. Define geometry — the validated waist design point
#
# The waist parameters are port-count agnostic: they touch only the body, via
# MMI._body_waist_bounds() / MMI.body_width_at(), which read neither n_inputs
# nor n_outputs. They are also mutually exclusive with the face-taper
# parameters, and the constructor enforces 2*Lt + Lw <= mm_length.
# ---------------------------------------------------------------------------
WAIST_WIDTH = 2.24       # Ww — interior body width (µm)
WAIST_LENGTH = 6.0       # Lw — uniform waist run (µm)
TRANSITION_LENGTH = 5.0  # Lt — each width ramp (µm)
MM_LENGTH = 22.0         # L  — total body length (µm)

mmi = MMI(
    n_inputs=2,             # 2x2 splitter
    n_outputs=2,
    mm_length=MM_LENGTH,    # MMI body length (µm) — 6.0 µm shorter than straight
    mm_width=2.8,           # MMI body width at the two faces (µm)
    body_waist_width=WAIST_WIDTH,
    body_waist_length=WAIST_LENGTH,
    body_transition_length=TRANSITION_LENGTH,
    waveguide_width=0.45,   # single-mode access waveguide width (µm)
    taper_length=4.0,       # input/output taper length (µm)
    taper_width=1.3,        # taper width at MMI junction (µm)
    waveguide_offset=1.5,   # center-to-center port spacing (µm)
    waveguide_height=0.22,  # Si layer thickness (µm) — SOI_220
    use_neighbor_aperture=True,  # confine the input source to one of the two inputs
)

# The uniform stub between the body face and the start of the ramp. It is not
# a free parameter: the hybrid variant of this device REFUSES to build unless
# s >= body_stub + dpml = 0.75 µm, because an FDTD section pads its drawn
# waveguide by dpml past the interior seam and a section boundary cannot land
# where the width is varying. s = 1.75 µm here, so the same geometry works in
# both scripts — which is the whole point of comparing them.
stub = (MM_LENGTH - 2 * TRANSITION_LENGTH - WAIST_LENGTH) / 2
half_body = mmi.mm_length / 2

# total_device_length (the criterion's convention) is defined in section 5,
# below the geometry; mm_length + 9.0 µm of fixed port hardware.
print(f"MMI body length: {mmi.mm_length:.2f} µm; device length "
      f"{mmi.mm_length + 9.0:.2f} µm "
      f"(straight control 37.00 µm — "
      f"{100*(37.0 - (mmi.mm_length + 9.0))/37.0:.1f}% shorter)")
print(f"Body width: face {mmi.width_at(-half_body):.3f} µm → "
      f"waist {mmi.width_at(0.0):.3f} µm → face {mmi.width_at(half_body):.3f} µm")
print(f"Body profile: stub {stub:.2f} | ramp {TRANSITION_LENGTH:g} | "
      f"waist {WAIST_LENGTH:g} | ramp {TRANSITION_LENGTH:g} | stub {stub:.2f} µm")
print(f"Input ports at y = {mmi.input_ys} µm; outputs at y = {mmi.output_ys} µm")

# ---------------------------------------------------------------------------
# 2. Configure FDTD solver — every setting matched to the straight 2x2's own
#    reference run EXCEPT decay_tol (see the module docstring).
# ---------------------------------------------------------------------------
solver = FDTDSolver(
    res=RESOLUTION,            # mesh resolution (pixels/µm)
    sim_width=5.0,             # simulation y-span (µm) — clears mm_width=2.8 plus margin
    sim_height=1.5,            # slab z-height reference (µm)
    dpml=0.5,                  # PML thickness, all 4 sides (µm)
    source_margin=0.0,         # source x-offset from left PML inner edge (µm)
    monitor_margin=0.25,       # monitor offset from PML inner edges, both sides (µm)
    three_d=True,       # full 3-D FDTD: real 220 nm Si core in oxide, PML on all six faces
    num_input_modes=1,         # input modes for normalization
    num_output_modes=1,        # modes extracted per output waveguide
    # Each output monitor spans one half of the domain (shared y=0 axis to the
    # near PML edge), so the two output waveguides (±0.75 µm) are captured in
    # isolation rather than as a two-guide supermode.
    output_monitor_centers=[2.5/2, -2.5/2],   # out_0 = top half (bar), out_1 = bottom half (cross)
    output_monitor_widths=[2.5, 2.5],
    material_stack=SOI_220,    # 220 nm SOI stack
    decay_tol=DECAY_TOL,       # ⚠ 1e-05, not 1e-11 — see the module docstring
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
plot_fdtd_geometry(
    mmi, solver=solver,
    title="Waisted 2×2 MMI Splitter — Top-down Geometry (FDTD)",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved → {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)  # µm, C-band

print(f"Sweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm")
print("(single broadband pulse; DFT monitors extract all wavelengths in one sim)")

# Excite the top input port; measure coupling to the fundamental mode at each
# output half.
#
# snapshot_every/snapshot_path make the decay tolerance VERIFIABLE rather than
# assumed: DFT monitors accumulate, so one run yields the whole
# transmission-vs-decay_tol curve. Snapshot args are method arguments, not
# solver attributes, so they are NOT in the cache fingerprint — adding them
# does not invalidate a stored run.
#
# The snapshot filename carries every parameter the result depends on,
# GEOMETRY FIRST. solve() APPENDS to this file, and this study lost data to
# that eight separate times: an untagged name let a res-15 smoke run and the
# res-40 production run interleave rows of different widths under one header,
# and later let a second design point append to the first one's rows.
_geo = (f"Ww{WAIST_WIDTH:g}_Lw{WAIST_LENGTH:g}"
        f"_Lt{TRANSITION_LENGTH:g}_L{MM_LENGTH:g}")
snapshot_path = os.path.join(
    output_dir,
    f"decay_snapshots_{_geo}_res{RESOLUTION:g}_nwl{N_WL}_tol{DECAY_TOL:g}.csv")

t0 = time.time()
smat = solver.solve(mmi, wavelengths, port="top", input_mode=1,
                    save_dir=save_dir, snapshot_every=SNAPSHOT_EVERY,
                    snapshot_path=snapshot_path)
elapsed = time.time() - t0
smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength
print(f"Solved in {elapsed:.0f}s")
print(f"Decay snapshots → {os.path.basename(snapshot_path)}")

# ---------------------------------------------------------------------------
# 5. Scores — the four criteria this device is judged on
# ---------------------------------------------------------------------------
# "out_0" / "out_1" match the order of output_monitor_centers: [+1.25, -1.25].
# Excited input is the top port, so out_0 (top) is the bar output and out_1
# (bottom) is the cross output. Keep the COMPLEX coefficients: the bar/cross
# relative phase is this device's whole function and |T| discards it.
c_bar   = smat.coeff(port="out_0", in_mode=1, out_mode=1)
c_cross = smat.coeff(port="out_1", in_mode=1, out_mode=1)

T_bar   = np.abs(c_bar) ** 2
T_cross = np.abs(c_cross) ** 2
T_total = T_bar + T_cross

# The PRE-REGISTERED criterion, and the control it is scored against. These
# literals are copied from `example_output/verification/mmi_2x2_waist_metric/
# criterion.json`, which was written BEFORE any design sweep was run and is
# scored automatically rather than eyeballed. They are frozen: a threshold
# chosen after seeing the results is not a criterion.
#
# The control is the STRAIGHT 2x2 — all four checks passed, device resolved
# (validation_markdown/resolved/mmi_2x2_validation.md) — so any deviation here
# is attributable to the waist rather than to an unvalidated baseline. Two of
# the three thresholds are defined RELATIVE to that control's own measured
# performance (loss: control + 0.1 dB; imbalance: 2x control), which is why
# they are not round numbers.
CONTROL_LENGTH_UM = 37.00     # straight 2x2: mm_length 28.0 + 9.0 of port hardware
CONTROL_PHASE_DEG = -89.74957705211592   # its own measured quadrature
MAX_EXCESS_LOSS_DB = 0.48090621635824105  # control 0.3809 + 0.1
MAX_IMBALANCE = 0.1371280008043103        # 2x control 0.06856
MAX_PHASE_DEV_DEG = 2.0

# ⚠ LENGTH CONVENTION. `MMI.total_length` counts mm_length + both tapers and
# stops there — it does NOT include the 0.5 µm of access waveguide drawn on
# each side. The criterion's own `total_device_length` DOES (its control is
# 37.00 µm for mm_length 28.0, not the 36.00 µm `total_length` returns), so
# scoring `mmi.total_length` against it would compare two different quantities
# and overstate the saving (18.9 % instead of the real 16.2 %).
FIXED_PORT_LENGTH = 2 * 4.0 + 2 * 0.5     # 2 tapers + 2 access stubs = 9.0 µm
total_device_length = mmi.mm_length + FIXED_PORT_LENGTH

dphi = np.angle(c_cross / c_bar, deg=True)

# ⚠ THE PRE-REGISTERED TEXT IS AMBIGUOUS between a band-MEAN and a band-MAX
# reading of the phase deviation, and that ambiguity is not academic — it is
# what killed the original 29.50 µm winner, and the 30.50 µm point passes one
# reading and fails the other. Both are reported here, always, so the verdict
# never depends on which one the reader assumes. THIS design point passes both.
metrics = {
    "total_device_length": float(total_device_length),
    "excess_loss_db": float(np.max(-10 * np.log10(np.clip(T_total, 1e-20, None)))),
    "imbalance": float(np.max(np.abs(T_bar - T_cross))),
    "phase_deg": float(np.mean(dphi)),
    "phase_mean_dev_deg": float(abs(np.mean(dphi) - CONTROL_PHASE_DEG)),
    "phase_max_dev_deg": float(np.max(np.abs(dphi - CONTROL_PHASE_DEG))),
    "phase_spread_deg": float(np.ptp(dphi)),
    "min_total": float(np.min(T_total)),
}

checks = {
    "shorter": metrics["total_device_length"] < CONTROL_LENGTH_UM,
    "excess_loss": metrics["excess_loss_db"] <= MAX_EXCESS_LOSS_DB,
    "imbalance": metrics["imbalance"] <= MAX_IMBALANCE,
    "phase(mean)": metrics["phase_mean_dev_deg"] <= MAX_PHASE_DEV_DEG,
    "phase(max)": metrics["phase_max_dev_deg"] <= MAX_PHASE_DEV_DEG,
}

wl_center = wavelengths[len(wavelengths) // 2]
idx = int(np.argmin(np.abs(wavelengths - wl_center)))

print(f"\nResults at λ = {wl_center:.3f} µm:")
print(f"  Bar output (top)     : {T_bar[idx]:.4f}")
print(f"  Cross output (bottom): {T_cross[idx]:.4f}")
print(f"  Total (sum)          : {T_total[idx]:.4f}")
print(f"  Relative phase       : {dphi[idx]:+.3f}°  (ideal ±90°)")
print(f"  (ideal 50/50 split → T_each=0.5, IL_total=0.0 dB)")

print(f"\nBand scores ({N_WL} wl over {WL_MIN}–{WL_MAX} µm):")
print(f"  total length     {metrics['total_device_length']:8.2f} µm "
      f"(control {CONTROL_LENGTH_UM:.2f} → "
      f"{100*(CONTROL_LENGTH_UM - metrics['total_device_length'])/CONTROL_LENGTH_UM:.1f}% shorter)")
print(f"  excess loss      {metrics['excess_loss_db']:8.4f} dB  (<= {MAX_EXCESS_LOSS_DB})")
print(f"  imbalance        {metrics['imbalance']:8.4f}     (<= {MAX_IMBALANCE})")
print(f"  phase            {metrics['phase_deg']:8.3f}°    (control {CONTROL_PHASE_DEG:.3f}°)")
print(f"  phase mean-dev   {metrics['phase_mean_dev_deg']:8.3f}°    (<= {MAX_PHASE_DEV_DEG})")
print(f"  phase MAX-dev    {metrics['phase_max_dev_deg']:8.3f}°    (<= {MAX_PHASE_DEV_DEG})")
print(f"  min total T      {metrics['min_total']:8.4f}")
fails = ", ".join(k for k, v in checks.items() if not v) or "-"
print(f"  VERDICT: {'PASS' if all(checks.values()) else 'FAIL'}  (failing: {fails})")

# ---------------------------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix.pkl")

# The hybrid script in this same folder loads its reference BY GEOMETRY-TAGGED
# NAME, so this second copy is what makes the two directly comparable.
smat.save(os.path.join(output_dir, f"smatrix_waisted_{_geo}.pkl"))

with open(os.path.join(output_dir,
                       f"check4_{_geo}_res{RESOLUTION:g}_nwl{N_WL}.json"), "w") as fh:
    json.dump({"config": {"res": RESOLUTION, "n_wl": N_WL,
                          "decay_tol": DECAY_TOL, "band": [WL_MIN, WL_MAX]},
               "design": {"Ww": WAIST_WIDTH, "Lw": WAIST_LENGTH,
                          "Lt": TRANSITION_LENGTH, "L": MM_LENGTH},
               "control": {"total_device_length": CONTROL_LENGTH_UM,
                           "phase_deg": CONTROL_PHASE_DEG},
               "thresholds": {"max_excess_loss_db": MAX_EXCESS_LOSS_DB,
                              "max_imbalance": MAX_IMBALANCE,
                              "max_phase_dev_deg": MAX_PHASE_DEV_DEG},
               "seconds": elapsed,
               "wavelengths": wavelengths.tolist(),
               "T_bar": T_bar.tolist(), "T_cross": T_cross.tolist(),
               "phase_deg_per_wl": dphi.tolist(),
               "metrics": metrics, "checks": checks,
               "passed": all(checks.values())}, fh, indent=2)
print(f"Scores saved → check4_{_geo}_res{RESOLUTION:g}_nwl{N_WL}.json")

# Transmission vs wavelength: both outputs + total, and the bar/cross phase
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_bar,   label="Bar output (top, mode 1)")
ax.plot(wavelengths, T_cross, label="Cross output (bottom, mode 1)", linestyle="--")
ax.plot(wavelengths, T_total, label="Total (bar + cross)", color="k", lw=1.5)
ax.axhline(0.5, color="gray", linestyle=":", lw=0.8, label="Ideal per-port (0.5)")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title("Waisted 2×2 MMI — FDTD Transmission")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)

# ⚠ |T| DISCARDS WHAT THIS DEVICE IS FOR. A 2x2 MMI coupler's function is the
# bar/cross relationship, and an ideal one sits at QUADRATURE. An intensity
# plot cannot see that at all, so the phase gets equal billing here.
ax = axes[1]
ax.plot(wavelengths, dphi, color="seagreen", label="arg(cross / bar)")
ax.axhline(CONTROL_PHASE_DEG, color="gray", linestyle=":", lw=0.8,
           label=f"Straight control ({CONTROL_PHASE_DEG:.2f}°)")
ax.fill_between(wavelengths,
                CONTROL_PHASE_DEG - MAX_PHASE_DEV_DEG,
                CONTROL_PHASE_DEG + MAX_PHASE_DEV_DEG,
                color="seagreen", alpha=0.10, lw=0,
                label=f"±{MAX_PHASE_DEV_DEG:g}° criterion")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Bar → cross relative phase (deg)")
ax.set_title("Quadrature vs Wavelength")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle(f"Waisted 2×2 MMI — {total_device_length:.2f} µm total "
             f"({100*(CONTROL_LENGTH_UM - total_device_length)/CONTROL_LENGTH_UM:.1f}% "
             f"shorter than straight)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "transmission.png"), dpi=300)
print(f"Transmission plot saved → {output_dir}/transmission.png")
