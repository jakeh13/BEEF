"""Full-device FDTD reference at a given y_gap -- the FDTD half of the bias study.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-25.

WHY. Before running a Bayesian optimization over (mm_width, y_gap) in pure EME,
we need to know whether the EME-vs-FDTD discrepancy DEPENDS ON y_gap.
`EMESolver.solve()` draws one input and one output branch per call and treats
the other six ports as absent, discarding port-to-port crosstalk -- and y_gap
IS the tip-to-tip separation that sets how big that crosstalk is. Coupling is
roughly exponential in gap, so the neglected term can vary by orders of
magnitude across the search range. A bias that varies with the optimized
parameter does NOT cancel in a relative comparison; it steers the optimizer.

Check 4 measured the error at y_gap=0.70 (up to 10 deg phase, 6.6 dB
imbalance). This run supplies the same reference at y_gap=0.40 and 1.00 so the
TREND can be read. Verdict rule:
  * error roughly flat in y_gap  -> largely cancels, BO on plain EME is sound
  * error varies strongly        -> BO would be biased; fix the solver, shrink
                                    the y_gap range, or optimize mm_width only

SETTINGS chosen from this session's own measurements, not inherited:
  * `decay_tol=1e-4` -- measured sufficient (0.015 dB worst case over both
    ports/all outputs/all wavelengths, vs the TRUE endpoint) by the snapshot
    curves from jobs 12272571/12272573. The old 1e-11 is UNREACHABLE: energy
    floors near ~5e-6 at res=40 and then grows. ~11-14 h here vs 63.6 h.
  * `decay_method="energy"` -- `stop_when_dft_decayed` never terminates on
    this device and prints no diagnostic while failing.
  * 96 ranks -- measured 1.95x over 48 at 98% efficiency (job 12217184); 192
    hangs in collective setup on scattered racks.
  * wavelengths 1.26-1.36 to MATCH the existing y_gap=0.70 reference so that
    point is reused rather than re-run (saves two ~14 h jobs). The BO itself
    will use 1.31 +- 0.02; measuring the bias on the WIDER band is
    conservative, since band-edge error is the worst case.

Geometry parameterization is identical to mmi_4x4_90hybrid_ygap_eme.py:
    waveguide_offset = mm_width / 4
    taper_width      = mm_width / 4 - y_gap

Run:
    MMI4X4_YGAP=0.4 MMI4X4_YGAP_PORT=bot python mmi_4x4_90hybrid_ygap_fdtd.py
    MMI4X4_YGAP=0.4 MMI4X4_YGAP_INPUT_IDX=2 python ...
"""
import os

import numpy as np

from beef import MMI, FDTDSolver, SOI_220

Y_GAP = float(os.environ.get("MMI4X4_YGAP", 0.70))
MM_WIDTH = float(os.environ.get("MMI4X4_MMWIDTH", 6.0))
MM_LENGTH = float(os.environ.get("MMI4X4_MMLENGTH", 0.75 * 108.9195 - 2.75))
TOL = float(os.environ.get("MMI4X4_YGAP_TOL", 1e-4))
PORT = os.environ.get("MMI4X4_YGAP_PORT", "bot")
INPUT_IDX = os.environ.get("MMI4X4_YGAP_INPUT_IDX")
INPUT_IDX = int(INPUT_IDX) if INPUT_IDX not in (None, "") else None

WAVEGUIDE_WIDTH = 0.35
TAPER_LENGTH = 1.0
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = WAVEGUIDE_OFFSET - Y_GAP
if not (WAVEGUIDE_WIDTH < TAPER_WIDTH < WAVEGUIDE_OFFSET):
    raise SystemExit(f"y_gap={Y_GAP} invalid: taper_width={TAPER_WIDTH:.3f}")

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                        "example_sim_cache"))
out_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                       "example_output", "mmi_4x4_90hybrid",
                                       "ygap_study"))
os.makedirs(out_dir, exist_ok=True)

mmi = MMI(n_inputs=4, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
          waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
          taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
          waveguide_height=0.22, use_neighbor_aperture=True)

solver = FDTDSolver(
    res=40.0, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    num_input_modes=1, num_output_modes=1,
    output_monitor_centers=list(mmi.output_ys),
    output_monitor_widths=[WAVEGUIDE_OFFSET] * 4,
    material_stack=SOI_220, decay_tol=TOL,
    decay_method="energy", energy_decay_dt=20.0,
)

_tag = f"idx{INPUT_IDX}" if INPUT_IDX is not None else PORT
_exc = f"input_idx={INPUT_IDX}" if INPUT_IDX is not None else f"port={PORT!r}"
print(f"YGAP_FDTD: y_gap={Y_GAP:.3f} mm_width={MM_WIDTH:.3f} "
      f"taper_width={TAPER_WIDTH:.3f} offset={WAVEGUIDE_OFFSET:.3f}", flush=True)
print(f"YGAP_FDTD: {_exc} decay_tol={TOL:.0e} energy-decay res=40", flush=True)

wavelengths = np.linspace(1.26, 1.36, 3)
snap = os.path.join(out_dir, f"fdtdcurve_ygap{Y_GAP:.2f}_{_tag}.csv")
common = dict(input_mode=1, save_dir=save_dir,
              snapshot_every=100.0, snapshot_path=snap)
if INPUT_IDX is not None:
    smat = solver.solve(mmi, wavelengths, input_idx=INPUT_IDX, **common)
else:
    smat = solver.solve(mmi, wavelengths, port=PORT, **common)
smat.check_passive()

T = {f"out_{j}": np.asarray(smat.transmission(port=f"out_{j}")) for j in range(4)}
C = {f"c_out_{j}": np.asarray(smat.coeff(port=f"out_{j}")) for j in range(4)}
for j in range(4):
    print(f"YGAP_FDTD: out_{j} |T|^2 = "
          f"{np.array2string(T[f'out_{j}'], precision=8)}", flush=True)
np.savez(os.path.join(out_dir, f"fdtd_ygap{Y_GAP:.2f}_{_tag}.npz"),
         y_gap=Y_GAP, mm_width=MM_WIDTH, taper_width=TAPER_WIDTH,
         mm_length=MM_LENGTH, wavelengths=wavelengths, tol=TOL, **T, **C)
print(f"YGAP_FDTD: done y_gap={Y_GAP:.2f} {_tag}", flush=True)
