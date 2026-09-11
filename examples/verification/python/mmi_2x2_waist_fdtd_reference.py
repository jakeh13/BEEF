"""Step 8 CHECK 4 -- full-device FDTD reference for the WAISTED 2x2.

validation_markdown/mmi_2x2_waist_body_validation.md step 8, check 4:
"The hybrid has to be checked against a single-domain brute-force run of the
same device -- that is the only check that tests the cascade as a whole rather
than its parts."

WHY THIS IS THE LOAD-BEARING CHECK HERE, not a closing formality.
`calibrate_body_neff=True` is HARD-REFUSED for a waisted body --
`assert_uniform_cross_section()` fires from `HybridSolver._calibrated_neff()`
(hybrid.py:1532): "a calibrated body propagates on one per-mode beta, which
only describes a uniform cross-section". That flag is exactly what closed
check 1 elsewhere in this repo: MPB's mode index and meep's time-stepped phase
differ ~2.3% per mode, and on the coupler the flag moved L50 from 5.37% off to
0.11%. So a waisted hybrid carries an UNCORRECTED phase error in the one
quantity this device's criterion scores most tightly (phase is already the
binding budget: mean-dev 1.696 of 2.0 deg).

A full-device FDTD run depends on NEITHER the cascade NOR the calibration. It
is the only available way to confirm the 29.50 um design point against
something that is not EME.

WHAT IS COMPARED. Straight control (mm_length=28.0) and the design point
(Ww=2.24, Lw=8.0, Lt=6.0, L=20.5) at MATCHED solver config -- the same
controlled-comparison discipline the whole waist study uses. Both are then
scored with the SAME metric functions the EME sweep used, so the FDTD and EME
numbers are directly comparable rather than merely adjacent.

Config is deliberately inherited from `examples/mmi_2x2/python/mmi_2x2_fdtd.py`
(res 40, dpml 0.5, sim_width 5.0, decay_tol 1e-11, split output monitors),
because the straight 2x2's OWN check-4 reference was produced by that script
(job 11585844_1, 1:16:45 on 48 cores). Changing the config would make this a
different measurement rather than a comparable one.

⚠ decay_tol=1e-11 is inherited and has no recorded justification anywhere in
the repo, and it is UNREACHABLE on long devices (measured: 84 um 4x4 floors
then diverges). It IS reachable on a short full-device solve (measured: 6.5 um
1x2, dft criterion fired at t=490). The 2x2 sits between those, so this run
passes `snapshot_every` -- because DFT monitors accumulate, one run yields the
whole transmission-vs-decay_tol curve, and the tolerance can be VERIFIED after
the fact instead of assumed. Snapshot args are method arguments, not solver
attributes, so they are not in the cache fingerprint.

Run:
    MMI2X2_REF_RES=15 MMI2X2_REF_NWL=6 python .../mmi_2x2_waist_fdtd_reference.py
    python .../mmi_2x2_waist_fdtd_reference.py          # real config
    MMI2X2_REF_CASES=waisted python .../...             # one case only
"""
import json
import os
import time

import numpy as np

from beef import MMI, FDTDSolver, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT_ROOT = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "verification"))
output_dir = os.path.join(_OUT_ROOT, "mmi_2x2_waist_fdtd_ref")
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_sim_cache"))
CRITERION_PATH = os.path.join(_OUT_ROOT, "mmi_2x2_waist_metric",
                              "criterion.json")


def _env_f(n, d):
    return float(os.environ.get(n, d))


def _env_i(n, d):
    return int(os.environ.get(n, d))


# --- config inherited from mmi_2x2_fdtd.py (the straight reference) ---------
RES = _env_f("MMI2X2_REF_RES", 40.0)
WL_MIN, WL_MAX = 1.530, 1.565
N_WL = _env_i("MMI2X2_REF_NWL", 51)
DECAY_TOL = _env_f("MMI2X2_REF_DECAYTOL", 1e-11)
SNAP_EVERY = _env_f("MMI2X2_REF_SNAP", 20.0)

# frozen ports, identical to the sweep and to the standard device
MM_WIDTH = 2.8
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 0.45, 1.3, 4.0
WAVEGUIDE_OFFSET, WAVEGUIDE_HEIGHT = 1.5, 0.22
FIXED_LENGTH = 2 * TAPER_LENGTH + 2 * 0.5      # 9.0 um of port hardware

# design point (user's call 2026-09-01: the MEAN phase reading governs)
WW = _env_f("MMI2X2_REF_WW", 2.24)
LW = _env_f("MMI2X2_REF_LW", 8.0)
LT = _env_f("MMI2X2_REF_LT", 6.0)
L_WAIST = _env_f("MMI2X2_REF_L", 20.5)
L_STRAIGHT = _env_f("MMI2X2_REF_LSTD", 28.0)

with open(CRITERION_PATH) as fh:
    CRIT = json.load(fh)
THR, CONTROL = CRIT["thresholds"], CRIT["control"]

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
# The EME criterion was scored on an 11-point grid over the same band. When
# N_WL = 51 every 5th point IS that grid exactly, so the comparison can be made
# on identical wavelengths rather than on interpolated ones.
SUB = 5 if N_WL == 51 else 1

CASES = {
    "straight": dict(mm_length=L_STRAIGHT),
    "waisted": dict(mm_length=L_WAIST, body_waist_width=WW,
                    body_waist_length=LW, body_transition_length=LT),
}
want = os.environ.get("MMI2X2_REF_CASES", "straight,waisted").split(",")
want = [c.strip() for c in want if c.strip()]
for c in want:
    if c not in CASES:
        raise SystemExit(f"unknown case {c!r}; choose from {sorted(CASES)}")


def metrics(cb, cc, mm_len):
    """The SAME metric definitions the EME sweep scored, so the two are
    directly comparable. Note phase_deg is the band MEAN (what score() has
    always applied, and the reading the user adopted) while phase_max_dev_deg
    is the max-over-band reading of the same pre-registered text."""
    Tb, Tc = np.abs(cb) ** 2, np.abs(cc) ** 2
    tot = Tb + Tc
    dphi = np.angle(cc / cb, deg=True)
    ph_ctrl = THR["control_phase_deg"]
    return {"total_device_length": mm_len + FIXED_LENGTH,
            "excess_loss_db": float(np.max(-10 * np.log10(np.clip(tot, 1e-20, None)))),
            "imbalance": float(np.max(np.abs(Tb - Tc))),
            "phase_deg": float(np.mean(dphi)),
            "phase_max_dev_deg": float(np.max(np.abs(dphi - ph_ctrl))),
            "phase_spread_deg": float(np.ptp(dphi)),
            "min_total": float(np.min(tot))}


def score(m):
    c = {"shorter": m["total_device_length"] < CONTROL["total_device_length"],
         "excess_loss": m["excess_loss_db"] <= THR["max_excess_loss_db"],
         "imbalance": m["imbalance"] <= THR["max_imbalance"],
         "phase": abs(m["phase_deg"] - THR["control_phase_deg"])
                  <= THR["max_phase_dev_deg"]}
    return all(c.values()), c


print("=" * 88)
print("CHECK 4 -- full-device FDTD reference, straight vs waisted 2x2")
print(f"  res={RES} dpml=0.5 sim_width=5.0 decay_tol={DECAY_TOL:g} "
      f"3-D, {N_WL} wl in [{WL_MIN}, {WL_MAX}]")
print(f"  straight L={L_STRAIGHT} -> {L_STRAIGHT + FIXED_LENGTH:.2f} um total")
print(f"  waisted  Ww={WW} Lw={LW} Lt={LT} L={L_WAIST} "
      f"-> {L_WAIST + FIXED_LENGTH:.2f} um total")
print(f"  cases: {want}")
print("=" * 88, flush=True)

results = {}
for name in want:
    kw = CASES[name]
    mmi = MMI(n_inputs=2, n_outputs=2, mm_width=MM_WIDTH,
              waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
              taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
              waveguide_height=WAVEGUIDE_HEIGHT,
              use_neighbor_aperture=True, **kw)
    solver = FDTDSolver(
        res=RES, sim_width=5.0, sim_height=1.5, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        num_input_modes=1, num_output_modes=1,
        output_monitor_centers=[2.5 / 2, -2.5 / 2],
        output_monitor_widths=[2.5, 2.5],
        material_stack=SOI_220, decay_tol=DECAY_TOL,
    )
    half = mmi.mm_length / 2
    print(f"\n[{name}] total {mmi.total_length:.2f} um; body width "
          f"in-face={mmi.width_at(-half):.3f} mid={mmi.width_at(0.0):.3f} "
          f"out-face={mmi.width_at(half):.3f} um", flush=True)
    # Config-keyed: solve() APPENDS to this file, so an untagged name lets a
    # cheap smoke and a production run interleave rows of DIFFERENT widths
    # (3 + 2*n_wl columns) under one header -- measured 2026-09-02, a res-15
    # 6-wl smoke and the res-40 51-wl run shared decay_snapshots_waisted.csv.
    # Third instance today of "an output path that is not keyed by config";
    # see also confirm_*.json and res_convergence_*.json.
    # decay_tol belongs in the name too: it is in the cache fingerprint
    # (cache/config.py:52), so runs at different tolerances are genuinely
    # different measurements, and res+n_wl alone would let a 1e-5 run append
    # to a 1e-11 run's file. Fixing only res/n_wl first was the same
    # half-measure this session already hit twice.
    # GEOMETRY belongs in the name too. res/n_wl/decay_tol were added after
    # earlier collisions, but a second DESIGN POINT at the same settings still
    # appended to the first one's file (measured 2026-09-04: the L=21.5 run
    # appended to the L=20.5 run's rows). Geometry is part of the cache
    # fingerprint, so by the rule "the filename carries everything the cache
    # fingerprints on" it has to be here. Fifth instance of this same mistake.
    geo = ("straight" if name == "straight"
           else f"Ww{WW:g}_Lw{LW:g}_Lt{LT:g}_L{L_WAIST:g}")
    snap = os.path.join(
        output_dir,
        f"decay_snapshots_{name}_{geo}_res{RES:g}_nwl{N_WL}"
        f"_tol{DECAY_TOL:g}.csv")
    t0 = time.time()
    smat = solver.solve(mmi, wavelengths, port="top", input_mode=1,
                        save_dir=save_dir, snapshot_every=SNAP_EVERY,
                        snapshot_path=snap)
    dt = time.time() - t0
    passive = smat.check_passive(warn=True)
    # out_0 = top half = bar (the excited input is the top port); out_1 = cross
    cb = smat.coeff(port="out_0", in_mode=1, out_mode=1)
    cc = smat.coeff(port="out_1", in_mode=1, out_mode=1)
    tag = (f"{name}" if name == "straight" else
           f"{name}_Ww{WW:g}_Lw{LW:g}_Lt{LT:g}_L{L_WAIST:g}")
    smat.save(os.path.join(output_dir, f"smatrix_{tag}.pkl"))
    print(f"[{name}] S-matrix -> smatrix_{tag}.pkl", flush=True)

    m = metrics(cb[::SUB], cc[::SUB], mmi.mm_length)
    passed, checks = score(m)
    fails = ",".join(k for k, v in checks.items() if not v) or "-"
    print(f"[{name}] solved in {dt:.0f}s  passive={passive}")
    print(f"[{name}] total_len {m['total_device_length']:.2f} um | "
          f"IL {m['excess_loss_db']:.4f} dB | imb {m['imbalance']:.4f} | "
          f"phase(mean) {m['phase_deg']:.3f} deg")
    print(f"[{name}] phase mean-dev "
          f"{abs(m['phase_deg'] - THR['control_phase_deg']):.3f} deg (SCORED) "
          f"| MAX-dev {m['phase_max_dev_deg']:.3f} deg "
          f"| spread {m['phase_spread_deg']:.3f} deg")
    print(f"[{name}] criterion: {'PASS' if passed else 'FAIL'} ({fails})",
          flush=True)
    results[name] = {"seconds": dt, "passive": bool(passive),
                     "wavelengths": wavelengths[::SUB].tolist(),
                     "T_bar": (np.abs(cb[::SUB]) ** 2).tolist(),
                     "T_cross": (np.abs(cc[::SUB]) ** 2).tolist(),
                     "passed": passed, "checks": checks, **m}

# --- the actual deliverable: FDTD vs EME on the same design point ----------
EME = {  # measured, res-40 grid + direct confirm (jobs 12674569 / 12687092)
    "straight": {"excess_loss_db": CONTROL["max_excess_loss_db"],
                 "imbalance": CONTROL["max_imbalance"],
                 "phase_deg": THR["control_phase_deg"]},
    "waisted": {"excess_loss_db": 0.3910, "imbalance": 0.0852,
                "phase_deg": -91.446},
}
print("\n" + "=" * 88)
print("FDTD vs EME on identical geometry (the check-4 comparison)")
print(f"{'case':>9} {'quantity':>14} {'FDTD':>10} {'EME':>10} {'delta':>10}")
for name in want:
    r, e = results[name], EME[name]
    for q in ("excess_loss_db", "imbalance", "phase_deg"):
        print(f"{name:>9} {q:>14} {r[q]:10.4f} {e[q]:10.4f} "
              f"{r[q] - e[q]:+10.4f}")
if len(want) == 2:
    ds = results["straight"]
    dw = results["waisted"]
    print(f"\nFDTD's own straight-vs-waisted verdict:")
    print(f"  length  {ds['total_device_length']:.2f} -> "
          f"{dw['total_device_length']:.2f} um "
          f"({ds['total_device_length'] - dw['total_device_length']:+.2f}, "
          f"{100*(ds['total_device_length'] - dw['total_device_length'])/ds['total_device_length']:.1f}% shorter)")
    print(f"  IL      {ds['excess_loss_db']:.4f} -> {dw['excess_loss_db']:.4f} dB")
    print(f"  imb     {ds['imbalance']:.4f} -> {dw['imbalance']:.4f}")
    print(f"  phase   {ds['phase_deg']:.3f} -> {dw['phase_deg']:.3f} deg")
    print("\n⚠ Read the STRAIGHT row first. It is the control: if FDTD and EME")
    print("  disagree there, the disagreement is the engines, not the waist,")
    print("  and the waisted row's delta must be read against it.")

# Geometry in the name (6th instance of this collision class): this file
# carries the WAISTED design point's results, so a second design point at the
# same solver settings must not overwrite the first. Measured 2026-09-04 --
# the L=21.5 run clobbered the L=20.5 run's summary.
_geo_tag = f"Ww{WW:g}_Lw{LW:g}_Lt{LT:g}_L{L_WAIST:g}"
path = os.path.join(output_dir,
                    f"check4_{_geo_tag}_res{RES:g}_nwl{N_WL}.json")
with open(path, "w") as fh:
    json.dump({"config": {"res": RES, "n_wl": N_WL, "decay_tol": DECAY_TOL,
                          "band": [WL_MIN, WL_MAX]},
               "design": {"Ww": WW, "Lw": LW, "Lt": LT, "L": L_WAIST},
               "criterion_used": CRIT, "eme_reference": EME,
               "results": results}, fh, indent=2)
print(f"\nwritten -> {path}")
