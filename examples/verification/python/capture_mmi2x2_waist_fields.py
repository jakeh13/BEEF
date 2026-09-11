"""Capture FDTD field maps for the waisted 2x2 -- input section AND full device.

PURPOSE. The hybrid-vs-FDTD comparison (job 12892459 vs 12797278) agrees to
~3% in |dT| but shows a flat ~1.5% splitting-ratio offset. A flat, dispersion-
free offset does not look like an uncorrected beta (which would accumulate with
propagation); it looks like something at the interfaces. Field maps are how you
tell those apart, so this captures:

  1. the hybrid's INPUT SECTION field  -- [-half_total, -half_body + body_stub]
  2. the SAME x-range from a FULL-DEVICE run

and the comparison is then exact: identical geometry, identical span, one
sectioned and one not. If cutting the device at the seam distorts the input
field, this is where it shows.

WHY solve_section FOR THE FULL DEVICE. `capture_fields_dir` exists only on
solve_section(), not on the full-device solve(). Rather than change shared
library code (other sessions have jobs importing beef/), the whole device is
passed AS a section spanning x0=-half_total to x1=+half_total. Same solver,
same physics, and it yields the field map solve() cannot.

⚠ The body/multimode comparison is NOT attempted here. Driving a waisted
(multi-cell) body from a captured interface field needs machinery that does not
exist: `project_interface_onto_body` propagates the body DIAGONALLY, which its
own docstring says "is exact for a uniform body ... but says nothing about a
tapered one". Step 7's multi-cell EME maps are the honest body picture for now.

Run:  MMI2X2_CAP_WHAT=input   (cheap-ish, one section)
      MMI2X2_CAP_WHAT=full    (~a full FDTD run)
"""
import os

import numpy as np

from beef import MMI, FDTDSolver, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                    "example_output", "mmi_2x2_tapered"))


def _env_f(n, d):
    return float(os.environ.get(n, d))


def _env_i(n, d):
    return int(os.environ.get(n, d))


WHAT = os.environ.get("MMI2X2_CAP_WHAT", "input")
RES = _env_f("MMI2X2_CAP_RES", 40.0)
WL = _env_f("MMI2X2_CAP_WL", 1.5475)          # band centre
N_MODES = _env_i("MMI2X2_CAP_NMODES", 14)
DECAY_TOL = _env_f("MMI2X2_CAP_DECAYTOL", 1e-5)
BODY_STUB = 0.25
WW, LW, LT, L = 2.24, 6.0, 6.0, 21.5

mmi = MMI(n_inputs=2, n_outputs=2, mm_length=L, mm_width=2.8,
          waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
          waveguide_offset=1.5, waveguide_height=0.22,
          body_waist_width=WW, body_waist_length=LW,
          body_transition_length=LT, use_neighbor_aperture=True)
half_body, half_total = mmi.mm_length / 2, mmi.total_length / 2

solver = FDTDSolver(res=RES, sim_width=5.0, sim_height=1.5, dpml=0.5,
                    source_margin=0.0, monitor_margin=0.25, three_d=True,
                    num_input_modes=1, num_output_modes=1,
                    material_stack=SOI_220, decay_tol=DECAY_TOL,
                    # ⚠ ALWAYS True, for the full device as well as the
                    # section. The two inputs sit waveguide_offset=1.5 um
                    # apart, so a full-width eigenmode source launches into
                    # BOTH at once -- mmi_2x2_fdtd.py's own docstring warns
                    # about this. Measured 2026-09-07: with it False the
                    # full-device capture lit both arms equally, making the
                    # comparison against the (single-arm) section meaningless.
                    half_domain_aperture=True)

if WHAT == "input":
    x0, x1 = -half_total, -half_body + BODY_STUB
    n_in, n_out = 1, N_MODES
    tag = "input_section"
elif WHAT == "full":
    x0, x1 = -half_total, +half_total
    # ⚠ n_out MUST EXCEED n_in. solve_section() "excites from whichever
    # boundary has FEWER modes" -- with n_in == n_out == 1 the tie resolved to
    # reverse=True, so the run launched from the OUTPUT end and propagated
    # backward. The captured field was physically correct but answered the
    # RECIPROCAL question: at the input facet it showed the device's split
    # across BOTH input arms (measured |Ey| ratio 0.986 between the arms,
    # where a forward single-input launch gives 0.001). Verified with
    # section_fdtd_layout(): n_out=1 -> reverse=True, n_out=2 -> reverse=False.
    n_in, n_out = 1, 2
    tag = "full_device"
else:
    raise SystemExit(f"MMI2X2_CAP_WHAT must be 'input' or 'full', got {WHAT!r}")

capdir = os.path.join(_OUT, f"captured_{tag}_res{RES:g}_wl{WL:.4f}")
os.makedirs(capdir, exist_ok=True)

print("=" * 80)
print(f"CAPTURE {tag}: x = [{x0:.3f}, {x1:.3f}] um  ({x1 - x0:.3f} um span)")
print(f"  Ww={WW} Lw={LW} Lt={LT} L={L}  res={RES} wl={WL} "
      f"decay_tol={DECAY_TOL:g}")
print(f"  modes in/out = {n_in}/{n_out}, half_domain_aperture="
      f"{WHAT == 'input'}")
print(f"  -> {capdir}")
print("=" * 80, flush=True)

mat = solver.solve_section(mmi, np.array([WL]), x0, x1,
                           n_modes_in=n_in, n_modes_out=n_out,
                           capture_fields_dir=capdir,
                           input_idx=1, output_idx=1)
print(f"transfer matrix shape {np.asarray(mat).shape}")
print("captured files:")
for f in sorted(os.listdir(capdir)):
    print(f"  {f}  ({os.path.getsize(os.path.join(capdir, f))} bytes)")
