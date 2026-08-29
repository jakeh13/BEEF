"""2x2 MMI waist body -- STEP 4: the figure of merit, measured and
pre-registered BEFORE any design sweep runs.

validation_markdown/mmi_2x2_waist_body_validation.md step 4.

WHY THIS EXISTS AS ITS OWN STEP. The 1x4's waist sweep ranked candidates by
`balance` alone. That is how two points with BETTER-than-control balance and
catastrophic totalT (0.47 / 0.18) surfaced as "BEST" in its summary table --
the ranking key could not see the failure mode that mattered. A 2x2 has
strictly more ways to be subtly wrong than a 1x4, because its observable is a
TRANSFER MATRIX: a waist could preserve every power number and still rotate
the bar/cross phase, which would be a real failure that a power-only metric
reports as a pass.

AND ONE MORE, FROM THIS DEVICE'S OWN HISTORY. resolved/mmi_2x2_validation.md
gap 2 records the coupler's lesson explicitly: "sweeping a quantity the
physics is not sensitive to (T_cross) hid a 7% error in one it is (dbeta)."
So the metric here does NOT stop at the S-matrix. It also reports
dbeta_01 = 2*pi/wl * (n_eff[0] - n_eff[1]) and L_pi = pi/dbeta_01 -- the
quantity that actually SETS this device's self-imaging length. A waist that
moves L_pi is doing the thing the whole idea is premised on; a waist that
leaves L_pi alone but changes T is doing something else, and the two want
different follow-ups.

WHAT THIS SCRIPT PRODUCES
  1. The straight control's own baseline across the band: T_bar, T_cross,
     excess loss, imbalance, bar/cross phase, dbeta_01, L_pi.
  2. The pre-registered pass/fail criterion, written in terms of that
     baseline, saved to disk as JSON so step 6 SCORES against it rather than
     eyeballing.

Run:
    python .../mmi_2x2_waist_metric.py                 # baseline + criterion
    MMI2X2_MET_NWL=3 python .../mmi_2x2_waist_metric.py    # cheap smoke
"""
import json
import os

import numpy as np

from beef import MMI, EMESolver, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "verification",
    "mmi_2x2_waist_metric"))
os.makedirs(output_dir, exist_ok=True)


def _env_f(name, default):
    return float(os.environ.get(name, default))


def _env_i(name, default):
    return int(os.environ.get(name, default))


# --- production constants, from the RESOLVED straight 2x2 -------------------
WL_MIN, WL_MAX = 1.530, 1.565
MM_LENGTH, MM_WIDTH = 28.0, 2.8
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 0.45, 1.3, 4.0
WAVEGUIDE_OFFSET, WAVEGUIDE_HEIGHT = 1.5, 0.22
STRAIGHT_LENGTH = 0.5

# Solver discretisation, copied from the RESOLVED device's own EME example
# (mmi_2x2_eme.py:63-73) rather than re-derived.
#
# ⚠ TWO TRAPS HERE, both hit on the first attempt (job 12487486) and both
# caught only because the resulting baseline was implausible (~2 dB excess
# loss and 0.18 imbalance on a device whose own doc calls it well-behaved):
#
#  1. THE TAPERS MUST BE SUBDIVIDED. num_cells=[1,15,1,15,1] -- 15 cells per
#     4 um port taper. Using 1 cell per taper makes each a single-step
#     staircase and cost ~2 dB of pure numerical loss. The BODY is the region
#     that gets 1 cell (it is uniform); the tapers are the ones that need
#     many. For a WAISTED body that flips: the body needs the segment-aligned
#     subdivision from step 3, via region_cell_lengths={2: [...]}.
#  2. N_MODES=14 IS NOT THIS. 14 is `mmi_2x2_hybrid.py`'s INTERFACE mode
#     count for the hybrid cascade -- a different quantity from the EME
#     example's per-region num_modes, which is [1, 7, 15, 7, 1] (7 in the
#     tapers, 15 in the body). Do not carry 14 across.
NUM_MODES = [1, 7, 15, 7, 1]
NUM_CELLS = [1, 15, 1, 15, 1]

RES = _env_f("MMI2X2_MET_RES", 40.0)      # the resolved device's own value
N_WL = _env_i("MMI2X2_MET_NWL", 11)
OXW, OXH = 5.0, 1.5

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
lengths = [STRAIGHT_LENGTH, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH,
           STRAIGHT_LENGTH]

mmi = MMI(n_inputs=2, n_outputs=2, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
          waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
          taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
          waveguide_height=WAVEGUIDE_HEIGHT)

print("=" * 78)
print("2x2 MMI waist -- STEP 4: figure of merit")
print(f"  straight control L={MM_LENGTH} W={MM_WIDTH}, res={RES}, "
      f"{N_WL} wl in [{WL_MIN}, {WL_MAX}]")
print("=" * 78, flush=True)


def solve_bar_cross(geometry, **solver_extra):
    """Excite the TOP input; measure both outputs.

    Convention is the 2x2 EME example's own (mmi_2x2_eme.py:104-105):
    input_idx=1 is the top input, output_idx=1 the top output ("bar"),
    output_idx=0 the bottom output ("cross").
    """
    solver = EMESolver(
        num_modes=NUM_MODES, num_cells=NUM_CELLS, lengths=lengths,
        oxide_width=OXW, oxide_height=OXH, res=RES,
        material_stack=SOI_220, mode_solver="tidy3d", **solver_extra)
    s_bar = solver.solve(geometry, wavelengths, ports=["bar"],
                         input_idx=1, output_idx=1)
    s_cross = solver.solve(geometry, wavelengths, ports=["cross"],
                           input_idx=1, output_idx=0)
    return s_bar.coeff("bar"), s_cross.coeff("cross")


def metrics(c_bar, c_cross):
    """Every number the criterion is written in terms of.

    Deliberately includes BOTH the power quantities and the phase, because a
    waist can hold all the power numbers and still rotate the transfer
    matrix -- see the module docstring.
    """
    T_bar = np.abs(c_bar) ** 2
    T_cross = np.abs(c_cross) ** 2
    total = T_bar + T_cross
    # phase of cross relative to bar, wrapped to (-180, 180]; an ideal 2x2
    # MMI sits at quadrature (|90| deg).
    dphi = np.angle(c_cross / c_bar, deg=True)
    return {
        "T_bar": T_bar,
        "T_cross": T_cross,
        "total": total,
        "excess_loss_db": -10 * np.log10(np.clip(total, 1e-20, None)),
        "imbalance": np.abs(T_bar - T_cross),
        "phase_deg": dphi,
        "phase_dev_from_quadrature": np.abs(np.abs(dphi) - 90.0),
    }


def dbeta_and_lpi(geometry):
    """dbeta_01 and L_pi = pi/dbeta_01 at band centre.

    This is the quantity that SETS the self-imaging length, and this device's
    own validation doc (gap 2) records that tracking only T can hide a large
    error in it. Reported so a waist that moves L_pi is distinguishable from
    one that merely moves T.
    """
    import meow as mw
    from beef._utils import extrude_polygons
    n_body = NUM_MODES[2]
    solver = EMESolver(
        num_modes=[n_body], num_cells=[1], lengths=[1.0],
        oxide_width=OXW, oxide_height=OXH, res=RES,
        material_stack=SOI_220, mode_solver="tidy3d")
    wl_c = float(np.mean(wavelengths))
    # one cell at the body centre cross-section
    gds = geometry.to_gds_eme(straight_length=STRAIGHT_LENGTH,
                              oxide_width=OXW, x_start=0.0, side="bar",
                              input_idx=1, output_idx=1)
    structs = extrude_polygons(gds, solver._extrusion_rules(np.array([wl_c])))
    mesh = solver._build_mesh()
    cells = mw.create_cells(structures=structs, mesh=mesh,
                            Ls=solver._flat_cell_lengths([sum(lengths)]))
    env = mw.Environment(wl=wl_c, T=solver.temperature)
    cs = mw.CrossSection.from_cell(cell=cells[len(cells) // 2], env=env)
    modes = solver.compute_modes(cs, n_body)
    n0, n1 = float(np.real(modes[0].neff)), float(np.real(modes[1].neff))
    dbeta = 2 * np.pi / wl_c * (n0 - n1)
    return {"n_eff_0": n0, "n_eff_1": n1, "dbeta_01": dbeta,
            "L_pi": float(np.pi / dbeta) if dbeta else float("nan")}


c_bar, c_cross = solve_bar_cross(mmi)
m = metrics(c_bar, c_cross)

print("\n--- straight control baseline ---")
print(f"{'wl':>8} {'T_bar':>9} {'T_cross':>9} {'total':>8} "
      f"{'IL(dB)':>8} {'imbal':>9} {'phase':>9}")
for i, wl in enumerate(wavelengths):
    print(f"{wl:8.4f} {m['T_bar'][i]:9.5f} {m['T_cross'][i]:9.5f} "
          f"{m['total'][i]:8.5f} {m['excess_loss_db'][i]:8.4f} "
          f"{m['imbalance'][i]:9.2e} {m['phase_deg'][i]:9.3f}")

print("\nband summary (worst case across the band):")
worst = {
    "max_excess_loss_db": float(np.max(m["excess_loss_db"])),
    "max_imbalance": float(np.max(m["imbalance"])),
    "max_phase_dev_deg": float(np.max(m["phase_dev_from_quadrature"])),
    "min_total": float(np.min(m["total"])),
}
for k, v in worst.items():
    print(f"  {k:24s} {v:.6g}")

print("\n--- dbeta_01 / L_pi at band centre ---")
print("(the quantity that SETS the self-imaging length -- gap 2's lesson:")
print(" tracking only T can hide a large error here)")
try:
    d = dbeta_and_lpi(mmi)
    for k, v in d.items():
        print(f"  {k:12s} {v:.6f}")
    print(f"  device mm_length / L_pi = {MM_LENGTH / d['L_pi']:.4f}")
except Exception as exc:                       # pragma: no cover - diagnostic
    print(f"  SKIPPED ({type(exc).__name__}: {exc})")
    d = None

# ---------------------------------------------------------------------------
# The pre-registered criterion
# ---------------------------------------------------------------------------
# Tolerances are set RELATIVE to the measured control, not to an ideal device.
# The question step 6 asks is "is this shorter at no worse performance", so
# the control is the right reference -- and anchoring on measured values means
# the bar cannot drift later to accommodate a result.
TOL_EXCESS_LOSS_DB = 0.10     # absolute dB the waist may add over control
TOL_IMBALANCE = 2.0           # multiple of the control's own worst imbalance
TOL_PHASE_DEG = 2.0           # absolute degrees of bar/cross phase drift

criterion = {
    "_what": "Pre-registered pass/fail for a waisted 2x2. Written BEFORE any "
             "design sweep, scored automatically, not eyeballed.",
    "_why": "The 1x4's sweep ranked on balance alone and surfaced two "
            "catastrophic-totalT points as BEST. A 2x2's observable is a "
            "transfer matrix, so phase is part of the metric, and dbeta_01 "
            "is reported because this device's own gap 2 records that "
            "tracking only T hid a 7% error in it on the coupler.",
    "control": {
        "mm_length": MM_LENGTH,
        "total_device_length": float(sum(lengths)),
        "res": RES,
        "n_wl": N_WL,
        "band": [WL_MIN, WL_MAX],
        **worst,
        **({} if d is None else d),
    },
    "pass_requires_ALL": {
        "shorter": "candidate total_device_length < control total_device_length",
        "excess_loss": f"max excess loss <= control ({worst['max_excess_loss_db']:.4f} dB) "
                       f"+ {TOL_EXCESS_LOSS_DB} dB",
        "imbalance": f"max imbalance <= {TOL_IMBALANCE}x control "
                     f"({TOL_IMBALANCE * worst['max_imbalance']:.3e})",
        "phase": f"max |bar/cross phase - control phase| <= {TOL_PHASE_DEG} deg",
    },
    "thresholds": {
        "max_excess_loss_db": worst["max_excess_loss_db"] + TOL_EXCESS_LOSS_DB,
        "max_imbalance": TOL_IMBALANCE * worst["max_imbalance"],
        "max_phase_dev_deg": TOL_PHASE_DEG,
        "control_phase_deg": float(np.mean(m["phase_deg"])),
    },
    "_ranking_key": "Rank ONLY among candidates that PASS all of the above. "
                    "Never rank on a single component -- that is the 1x4's "
                    "exact mistake. Among passers, shortest device wins.",
    "_note_on_grid": "The constructor enforces 2*Lt + Lw <= mm_length, so Lt "
                     "and Lw are NOT independently sweepable at the control's "
                     "own length -- e.g. Lt=6/Lw=8 already uses 20 of 28 um. "
                     "Generate the grid with that constraint applied, not as "
                     "a full cartesian product.",
}

path = os.path.join(output_dir, "criterion.json")
with open(path, "w") as fh:
    json.dump(criterion, fh, indent=2)

print("\n--- PRE-REGISTERED CRITERION ---")
for k, v in criterion["pass_requires_ALL"].items():
    print(f"  {k:14s} {v}")
print(f"\nranking: {criterion['_ranking_key']}")
print(f"\nwritten -> {path}")
print("\nStep 6 must LOAD this file and score against it. Changing a threshold")
print("after seeing results is the thing pre-registration exists to prevent.")
