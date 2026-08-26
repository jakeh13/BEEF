"""Is `decay_tol=1e-11` actually necessary for Step 6b's S-parameters?

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-23.

THE QUESTION. Step 6b's cost is dominated entirely by the tail: after the
pulse exits at t~340-420 about 5.7% of the energy stays trapped and bleeds
off at k ~1e-3, so reaching decay_tol=1e-11 needs t~19,000-20,000 (~70 h at
96 ranks). Two 48 h attempts were cancelled unfinished. **`decay_tol=1e-11`
has no recorded justification anywhere in this repo** -- an audit found it
inherited, never derived, and the only prior attempt to loosen it (modal
multiplexer, 1e-11 -> 1e-9) was rejected for the OPPOSITE failure: a run
stopping too early, which `minimum_run_time` was then added to fix. That
says nothing about whether 1e-11 is needed here.

WHAT ACTUALLY MATTERS is not the residual energy, it is whether that residual
still changes the DFT-accumulated S-parameters. Trapped energy sitting in a
mode that does not couple to the output monitors contributes nothing to the
answer no matter how slowly it decays.

METHOD. Solve the identical device at a LADDER of decay tolerances and
compare the resulting S-parameters against each other. From the measured
decay curve the stopping points are cheap at the loose end:

    decay_tol   reached near   approx wall @ 96 ranks
    1e-3        t ~ 2,100      ~2 h
    1e-4        t ~ 3,200      ~3 h
    1e-5        t ~ 6,000      ~6 h

If |S| and arg(S) agree between successive rungs to far better than the
precision this device's conclusions rest on (imbalance ~0.1 dB, phase ~1 deg,
CMRR ~0.5 dBe), then the tail is irrelevant to the physics and 1e-11 is
simply buying noise at ~70 h a run. If they move materially between rungs,
the tail matters and the tolerance is doing real work.

This is deliberately a CONVERGENCE ladder, not a single loose run compared to
nothing -- with no converged reference available (that is the whole problem),
successive-rung agreement is the only honest evidence.

CACHING. save_dir IS passed. `decay_tol` is part of the cache fingerprint
(beef/cache/config.py `_METHOD_ATTRS`), so each rung gets its own entry and
none of them can ever be served to a production 1e-11 request. And if the
ladder shows a loose tolerance suffices, the corresponding entry is not a
throwaway -- it becomes the real Step 6b result.

Run (one rung per invocation):
    MMI4X4_TOL=1e-4 python mmi_4x4_90hybrid_decaytol_convergence.py
"""
import os

import numpy as np

from beef import MMI, FDTDSolver, SOI_220

TOL = float(os.environ.get("MMI4X4_TOL", 1e-4))
RES = float(os.environ.get("MMI4X4_TOL_RES", 40.0))
PORT = os.environ.get("MMI4X4_TOL_PORT", "bot")
# Port 3 (the paper's second driven input) is an INTERIOR port that "top"/"bot"
# cannot reach -- the production 6b batch excites it as `input_idx=2`, not via
# a `port` string (see batch/run_mmi_4x4_90hybrid_fdtd.batch:108). Passing
# port="idx2" would silently excite the wrong thing, so accept an explicit
# index and route it to the right kwarg.
INPUT_IDX = os.environ.get("MMI4X4_TOL_INPUT_IDX")
INPUT_IDX = int(INPUT_IDX) if INPUT_IDX not in (None, "") else None

WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 6.0
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 0.8
TAPER_LENGTH = 1.0
L_PI = 108.9195
DELTA_L = -2.75
MM_LENGTH = 0.75 * L_PI + DELTA_L   # 78.9396

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                        "example_sim_cache"))

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.22, use_neighbor_aperture=True,
)

solver = FDTDSolver(
    res=RES, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    num_input_modes=1, num_output_modes=1,
    output_monitor_centers=list(mmi.output_ys),
    output_monitor_widths=[WAVEGUIDE_OFFSET] * 4,
    material_stack=SOI_220, decay_tol=TOL,
    decay_method="energy", energy_decay_dt=20.0,
)

_exc = f"input_idx={INPUT_IDX}" if INPUT_IDX is not None else f"port={PORT!r}"
print(f"TOLCONV: decay_tol={TOL:.0e} {_exc} res={RES} dpml=0.5", flush=True)

wavelengths = np.linspace(1.26, 1.36, 3)
# One run, sampled densely, instead of one run per tolerance -- see
# FDTDSolver.solve()'s `snapshot_every`. Each sample carries the energy ratio
# that would have stopped the run there, so the resulting curve IS
# transmission-vs-decay_tol, read off a single simulation.
SNAP_EVERY = float(os.environ.get("MMI4X4_TOL_SNAP_EVERY", 50.0))
_outdir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                       "example_output", "mmi_4x4_90hybrid"))
os.makedirs(_outdir, exist_ok=True)
_tag = ("idx%d" % INPUT_IDX) if INPUT_IDX is not None else PORT
snap_path = os.path.join(_outdir, f"tolcurve_res{int(RES)}_{_tag}.csv")
print(f"TOLCONV: snapshots every {SNAP_EVERY} t-units -> {snap_path}", flush=True)

_common = dict(input_mode=1, save_dir=save_dir,
               snapshot_every=SNAP_EVERY, snapshot_path=snap_path)
if INPUT_IDX is not None:
    smat = solver.solve(mmi, wavelengths, input_idx=INPUT_IDX, **_common)
else:
    smat = solver.solve(mmi, wavelengths, port=PORT, **_common)

_n = 0
if os.path.exists(snap_path):
    with open(snap_path) as _fh:
        _n = max(0, sum(1 for _ in _fh) - 1)
if _n:
    print(f"TOLCONV: {_n} snapshots written to {snap_path}", flush=True)
else:
    # A cache hit returns without time-stepping, so no curve is produced.
    # That is correct, not a failure -- the S-matrix below is still the real
    # cached result. Delete the cache entry (or vary a tracked attr) to force
    # a fresh run if the curve itself is what you need.
    print("TOLCONV: no snapshots -- solve() returned from cache without "
          "running; S-matrix below is the cached result", flush=True)
smat.check_passive()

out = {}
for j in range(4):
    c = smat.coefficients(port=f"out_{j}") if hasattr(smat, "coefficients") else None
    T = np.asarray(smat.transmission(port=f"out_{j}"))
    out[f"out_{j}"] = T
    print(f"TOLCONV: tol={TOL:.0e} out_{j} |T|^2 = "
          f"{np.array2string(T, precision=8)}", flush=True)

np.savez(os.path.join(save_dir, "..", "example_output", "mmi_4x4_90hybrid",
                      f"decaytol_{TOL:.0e}_res{int(RES)}_"
                      f"{('idx%d' % INPUT_IDX) if INPUT_IDX is not None else PORT}.npz"),
         tol=TOL, res=RES, wavelengths=wavelengths,
         **{k: v for k, v in out.items()})
print(f"TOLCONV: done tol={TOL:.0e}", flush=True)
