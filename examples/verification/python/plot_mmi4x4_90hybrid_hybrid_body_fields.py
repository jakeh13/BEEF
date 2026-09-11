"""4x4 90-degree hybrid -- hybrid body field reconstruction + output phase.

Companion to plot_mmi1x2_hybrid_body_selfimaging.py / plot_dc_hybrid_body_fields.py,
pinned to mmi_4x4_90hybrid_wg_eval.py's own construction (the BO campaign's
trusted hybrid-path evaluator, not the disqualified pure-EME
mmi_4x4_90hybrid_eme.py). See beef/plot/hybrid_body.py for why plot_eme_fields()
alone is not enough here (meow-frame labels, single-mode-only excitation) and
the standing caveat that the excitation is an APPROXIMATION to a Gaussian beam
the width of the driven taper, not the cascade's own interface field --
the neighbouring FDTD input/output sections carry no field data anywhere in
the library yet (hybrid_field_extraction_validation.md step 2, unbuilt).

WHY A PHASE PANEL, NOT JUST |E|^2. plot_hybrid_body_fields()'s own 3-panel
figure is intensity-only. A 90-degree hybrid's whole function is the RELATIVE
PHASE between outputs -- the exact quantity this campaign's own item 1 (see
mmi_4x4_90hybrid_wg_eval.py's docstring) measured pure-EME getting wrong by
44-62 degrees against real FDTD, with intensities still looking reasonable.
An intensity-only picture cannot validate this device at all. The helper
already returns the complex field E for exactly this reason (see its own
docstring); this script adds a 4th panel via beef.plot.hybrid_body.port_phases()
reading the body's OWN accumulated phase at the four output y-positions,
relative to the port-0 exit.

WHY OFF-AXIS, AND WHY TWO LAUNCHES. mmi.input_ys gives the true port centres
(reading them off the geometry rather than hand-deriving -- mmi_1x2_hybrid.py's
own comment claiming output_ys=[-wo,+wo] is wrong by 2x, per the handoff notes
in hybrid_body_field_plotting_handoff.md). The campaign's own DRIVEN=[0,2]
(paper ports 1/3) are the two ports whose cross-excitation phase difference
the campaign actually measures (cross_excitation_phase_diff_deg), so both get
their own launch + plot rather than picking one arbitrarily.

WHY THE BEST POINT IS READ LIVE, NOT HARDCODED. The BO campaign
(mmi_4x4_90hybrid_bo.py) was still running when this script was written
(round 6 in flight: beef_bo_W6.00_y0.50 / W6.50_y0.56 / W6.00_y0.51 /
W5.00_y0.43). Re-reads bo_state.json's own current best at run time, and
stamps the geometry into the output tag so two runs at different points never
collide -- override with MMI4X4_FIELD_W/_YGAP/_LEN for a specific point instead.

Cost: one cross-section, N_MODES=16 meep/MPB modes, one wavelength -- a single
mode solve, not a cascade rerun (no FDTD, no length sweep). The 1x4's own
12 modes/12um body at res 30 measured ~3 min; this device's more modes on a
wider body will cost more but is the same class of job, expected minutes not
hours -- confirm on the actual run rather than assuming.

Run:
    python .../plot_mmi4x4_90hybrid_hybrid_body_fields.py
    MMI4X4_FIELD_W=6.0 MMI4X4_FIELD_YGAP=0.65 MMI4X4_FIELD_LEN=77.40 \\
        python .../plot_mmi4x4_90hybrid_hybrid_body_fields.py   # explicit point
"""
import json
import os

import numpy as np

from beef import MMI, EMESolver, SOI_220, plot_hybrid_body_fields
from beef.plot.hybrid_body import port_phases

_HERE = os.path.dirname(os.path.abspath(__file__))
_CAMPAIGN_STATE = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid",
    "wg_ygap_opt", "bo_state.json"))

# -----------------------------------------------------------------------------
# Constants -- pinned to mmi_4x4_90hybrid_wg_eval.py's own construction.
WAVEGUIDE_WIDTH = 0.35
BODY_STUB = 0.25
RESOLUTION = 40.0
N_MODES = 16
WL_PLOT = float(os.environ.get("MMI4X4_FIELD_WL", 1.31))  # O-band centre
DRIVEN = [0, 2]  # paper ports 1/3 -- see mmi_4x4_90hybrid_wg_eval.py


# -----------------------------------------------------------------------------
# 1. Which (mm_width, y_gap, mm_length) to plot -- the campaign's own current
#    best, read live, unless overridden via env vars.
def campaign_best():
    """Current best (mm_width, y_gap, mm_length), read live -- see module
    docstring on why this is not hardcoded."""
    with open(_CAMPAIGN_STATE) as fh:
        results = json.load(fh)
    if not results:
        raise SystemExit(f"no BO results yet in {_CAMPAIGN_STATE}")
    best = min(results, key=lambda r: r["objective"])
    return best["mm_width"], best["y_gap"], best["best_mm_length"], best


_env_w = os.environ.get("MMI4X4_FIELD_W")
_env_y = os.environ.get("MMI4X4_FIELD_YGAP")
_env_l = os.environ.get("MMI4X4_FIELD_LEN")
if _env_w and _env_y and _env_l:
    MM_WIDTH, Y_GAP, MM_LENGTH = float(_env_w), float(_env_y), float(_env_l)
    _best_record = None
else:
    MM_WIDTH, Y_GAP, MM_LENGTH, _best_record = campaign_best()

WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = WAVEGUIDE_OFFSET - Y_GAP
if TAPER_WIDTH <= WAVEGUIDE_WIDTH:
    raise SystemExit(
        f"INVERSE TAPER at the current campaign best: y_gap={Y_GAP} "
        f"mm_width={MM_WIDTH} gives taper_width={TAPER_WIDTH:.3f} <= "
        f"waveguide_width={WAVEGUIDE_WIDTH} -- campaign state looks corrupted.")

output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid",
    "wg_ygap_opt", "body_fields"))
os.makedirs(output_dir, exist_ok=True)

# -----------------------------------------------------------------------------
# 2. Geometry + body section boundaries -- same construction as wg_eval.py.
mmi = MMI(n_inputs=4, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
          waveguide_width=WAVEGUIDE_WIDTH, taper_length=1.0,
          taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
          waveguide_height=0.22, use_neighbor_aperture=True)

half_body = mmi.mm_length / 2
x0, x1 = -half_body + BODY_STUB, +half_body - BODY_STUB  # body is k=1 of 3

# -----------------------------------------------------------------------------
# 3. EME body solver -- same config as wg_eval.py's own body_solver.
solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=MM_WIDTH + 3.0, oxide_height=2.0,
    res=RESOLUTION, material_stack=SOI_220, mode_solver="meep",
)

# -----------------------------------------------------------------------------
# 4. Run info.
print(f"4x4 hybrid body field plot: mm_width={MM_WIDTH:.2f} y_gap={Y_GAP:.2f} "
      f"mm_length={MM_LENGTH:.4f} taper_width={TAPER_WIDTH:.3f}")
if _best_record is not None:
    print(f"  (campaign best as of this run: objective={_best_record['objective']:.4f} "
          f"max|dphi|={_best_record['max_dphi']:.3f} min_cmrr={_best_record['min_cmrr']:.3f})")
print(f"  body x in [{x0:+.2f}, {x1:+.2f}] um "
      f"({x1 - x0:.2f} of {mmi.total_length:.2f} um device)")
print(f"  input_ys={['%+.3f' % y for y in mmi.input_ys]}  "
      f"output_ys={['%+.3f' % y for y in mmi.output_ys]}")

# -----------------------------------------------------------------------------
# 5. Reconstruct the body field for each driven input port + read its output
#    phases -- the actual deliverable, one EME mode solve per launch.
results = {}
for k in DRIVEN:
    print(f"\n--- launching input port {k} (y={mmi.input_ys[k]:+.3f}) ---")
    res = plot_hybrid_body_fields(
        mmi, solver, x0, x1, WL_PLOT, N_MODES, output_dir,
        tag=f"W{MM_WIDTH:.2f}_y{Y_GAP:.2f}_L{MM_LENGTH:.2f}_in{k}",
        # OFF-AXIS: the driven port's own y, not 0.0 -- see module docstring.
        excitation="gaussian", excitation_center=mmi.input_ys[k],
        excitation_fwhm=TAPER_WIDTH,
        wavelengths=np.array([WL_PLOT]),
        port_ys=mmi.output_ys, body_half_width=mmi.mm_width / 2,
        device_label="4x4 90-deg hybrid", total_length=mmi.total_length,
    )
    assert res["captured_power"] <= 1.0 + 1e-6, (
        f"captured power {res['captured_power']:.4f} > 1 for input {k} -- "
        "the projection inner product is wrong; see body_mode_weights()'s docstring")
    results[k] = res

    ph = port_phases(res, mmi.output_ys, reference=0)
    print(f"  output amplitudes |E|: " + ", ".join(f"{a:.4f}" for a in ph["amplitude"]))
    print(f"  output phase (deg, rel. to output 0): "
          + ", ".join(f"{p:+.2f}" for p in ph["phase"]))
    results[k]["port_phases"] = ph

# -----------------------------------------------------------------------------
# 6. Phase panel -- relative output phase for each driven input. This is the
#    actual check a 90-degree hybrid needs; the |E|^2 panels above cannot see
#    it (see module docstring).
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6.5, 4.2))
width = 0.35
x = np.arange(len(mmi.output_ys))
for i, k in enumerate(DRIVEN):
    ph = results[k]["port_phases"]
    weak = ph["amplitude"] < 0.05 * ph["amplitude"].max()
    bars = ax.bar(x + (i - 0.5) * width, ph["phase"], width,
                  label=f"input {k} (y={mmi.input_ys[k]:+.2f})")
    for b, w in zip(bars, weak):
        b.set_alpha(0.35 if w else 0.9)
ax.set_xticks(x)
ax.set_xticklabels([f"out{j}\n(y={y:+.2f})" for j, y in enumerate(mmi.output_ys)])
ax.set_ylabel("relative phase (deg, vs output 0)")
ax.axhline(0, color="gray", lw=0.8)
ax.set_title(
    "4x4 90-deg hybrid body -- output phase from the reconstructed field\n"
    "(faint bars = amplitude < 5% of max at that launch, phase not meaningful "
    "there — see port_phases()'s own docstring)", fontsize=9)
ax.legend(fontsize=8)
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
phase_path = os.path.join(
    output_dir, "fields",
    f"hybrid_body_phase_W{MM_WIDTH:.2f}_y{Y_GAP:.2f}_L{MM_LENGTH:.2f}.png")
fig.savefig(phase_path, dpi=160)
plt.close(fig)
print(f"\nphase panel -> {phase_path}")

# -----------------------------------------------------------------------------
# 7. Standing caveat.
print("\nchecks passed")
print("\nSTANDING CAVEAT: this is the EME BODY ONLY. The input/output FDTD taper "
      "sections carry no field data anywhere in the library yet (see "
      "hybrid_field_extraction_validation.md step 2) -- this is not a whole-"
      "device map, and the launch is an approximate Gaussian, not the "
      "cascade's own interface field.")
