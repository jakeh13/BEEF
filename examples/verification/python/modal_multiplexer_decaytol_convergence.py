"""Modal multiplexer — is decay_tol=1e-9 stopping the FDTD run too early,
too late, or about right? Direct answer to the open "run long enough?"
question in modal_multiplexer_validation.md gap 9/11.

WHY THIS, NOW. gap 9's res-40/TAPER_L=2.0 FDTD reference (job 12889324,
decay_tol=1e-9) undershoots the paper's seed transmission on both channels
(TE00->A 0.79 dB short, TE01->B 1.80 dB short). Before spending more compute
chasing geometry or solver-class explanations, rule in/out the cheapest
remaining lever: is 1e-9 simply stopping before the real signal has finished
arriving?

This device's OWN history makes "just run tighter" a live hazard, not a free
win: the exact same script at decay_tol=1e-11 (job 11707932) DID fully decay
(TE0 at t=77,197, 24-26x longer than a sane res-10 smoke run needed) and
produced a *worse*, not better, result -- near-exactly-zero transmission at
4 of 6 channels, a pattern that reads as underflow/a spurious slowly-decaying
component swamping the real signal, not genuine physics (gap 3, unconfirmed
mechanism). So the direction of the effect from "running longer" is NOT
known in advance here the way it might be on a device without that history.

METHOD, costs ONE run per input mode, not one run per tolerance.
`FDTDSolver.solve(snapshot_every=, snapshot_path=)` decomposes the monitors
every SNAP_EVERY meep time units mid-run and appends {t, energy, energy_ratio,
T} to a CSV. DFT monitors ACCUMULATE, so the sample at time t equals exactly
what stopping the run at t would have returned -- one run sweeps the whole
transmission-vs-decay_tol curve at once. Read the resulting CSV for:
  (a) has T already plateaued by the t where 1e-9 would have stopped it
      (energy_ratio column crosses 1e-9)? If yes, running longer does NOT
      explain the gap -- rule this lever out.
  (b) is T still trending toward the paper's higher value past that point?
      Supports "not running long enough" -- the concrete, testable claim.
  (c) does T instead continue toward job 11707932's near-zero collapse?
      Confirms 1e-9 is already past the point of diminishing/reversing
      returns, and that pre-existing pathology, not under-running, is what a
      tighter tolerance would rediscover.

decay_tol itself is left TIGHT (1e-11, matching 11707932's own value so this
new curve is directly comparable to that prior run's single terminal point)
so the run continues PAST 1e-9's stopping point and covers every looser
tolerance in one pass. minimum_run_time is dropped to 0 -- the snapshot CSV
already tells the whole story from t=0, no floor needed for a diagnostic.

TE1 input only (not TE0): TE01->B is the paper channel with the larger
shortfall (1.80 dB vs 0.79 dB) and this device's own known worst
recovered-power channel (gap 11: FDTD keeps only 43.6% of TE1's launched
power at TAPER_L=2.0), so it is the more informative one to check first.
TE0 can follow if this run is inconclusive.

Geometry, solver config: byte-for-byte modal_multiplexer_fdtd.py's own
(TAPER_L=2.0, res 40 forced -- not that script's own res=20 default -- to
match the actual job 12889324 reference this is diagnosing), reproduced here
rather than imported, for the same reason measure_mux_connector_mode_cost.py
execs instead of imports it: no `__main__` guard, so importing would run the
whole cascade.

Run (see batch/verification/run_modal_multiplexer_decaytol_convergence.batch):
    srun --mpi=pmi2 python3 verification/python/modal_multiplexer_decaytol_convergence.py
"""
import os

import numpy as np

from beef import MMI, MMIPort, CompoundGeometry, WaveguideRun, FDTDSolver, SOI_160

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                   "example_output", "modal_multiplexer"))
os.makedirs(OUT, exist_ok=True)

# --- pinned to modal_multiplexer_fdtd.py, TAPER_L=2.0 verified state -------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.26, 1.36, 11
TAPER_L = 2.0
WO = 1.0375
BETWEEN_WG_LENGTH = 3.0
MM_WIDTH = 2.8

SNAP_EVERY = float(os.environ.get("MUX_SNAP_EVERY", "100.0"))
# 100.0, matching mmi_4x4_90hybrid_decaytol_convergence.py's own choice for
# an 84 um device -- this device is 64 um, same order of magnitude, and the
# 1x2's proportionally finer 10.0 was chosen because THAT device is 6.5 um,
# an order of magnitude shorter still. Overridable.

BRANCH_SPECS = [
    (+WO, 0.33, 0.725),   # top
    (0.0, 0.35, 1.15),    # mid
    (-WO, 0.35, 0.725),   # bot
]


def _branch_ports():
    return [MMIPort(y=y, waveguide_width=ww, taper_length=TAPER_L, taper_width=tw)
            for y, ww, tw in BRANCH_SPECS]


mmi1 = MMI(
    n_inputs=1, n_outputs=3,
    mm_length=16.92, mm_width=MM_WIDTH,
    inputs=[MMIPort(y=0.0, waveguide_width=0.7, taper_length=TAPER_L, taper_width=1.4)],
    outputs=_branch_ports(),
    waveguide_height=0.16,
)
mmi2 = MMI(
    n_inputs=3, n_outputs=3,
    mm_length=34.24, mm_width=MM_WIDTH,
    inputs=_branch_ports(),
    outputs=_branch_ports(),
    waveguide_height=0.16,
    use_neighbor_aperture=True,
)
connector = WaveguideRun(
    length=BETWEEN_WG_LENGTH,
    ports=[(y, ww) for y, ww, _ in BRANCH_SPECS],
)
mux = CompoundGeometry([mmi1, connector, mmi2])

solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=5.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25,
    three_d=True,
    num_input_modes=2, num_output_modes=1,
    output_monitor_centers=[+WO, 0.0, -WO],
    material_stack=SOI_160,
    decay_tol=1e-11,      # TIGHT deliberately -- see module docstring
    minimum_run_time=0.0, # no floor needed -- the CSV covers t=0 onward
)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
snap_path = os.path.join(OUT, "mux_tolcurve_res40_TE1.csv")

print(f"modal mux full-device FDTD, res={RESOLUTION:.0f}, TAPER_L={TAPER_L}, "
      f"device length {mux.total_length:.2f} um")
print(f"input_mode=2 (TE1) only -- see module docstring for why")
print(f"snapshot_every={SNAP_EVERY}, decay_tol=1e-11 (deliberately tight, "
      f"matches job 11707932's own value)")
print(f"snapshots -> {snap_path}\n")

smat = solver.solve(mux, wavelengths, input_mode=2,
                    snapshot_every=SNAP_EVERY, snapshot_path=snap_path)

PORT_NAMES = ["top", "mid", "bot"]
i_c = len(wavelengths) // 2
print(f"\nfinal (fully decayed) at {wavelengths[i_c]:.3f} um:")
tot = 0.0
for k, name in enumerate(PORT_NAMES):
    T = smat.transmission(port=f"out_{k}", in_mode=2)[i_c]
    tot += T
    print(f"  TE1 -> {name}: T={T:.6f}")
print(f"  total: {tot:.6f}")
print(f"\nCSV: {snap_path}")
print("Read with: check the energy_ratio column for the row nearest 1e-9 -- "
      "that row's T is what job 12889324 actually reported. Compare against "
      "later rows (tighter tolerance, more decay time) to see whether T "
      "trends toward the paper's higher value, has already plateaued, or "
      "collapses the way job 11707932's terminal (t=77,197) point did.")
