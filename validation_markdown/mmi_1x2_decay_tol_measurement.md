# 1×2 MMI — what decay_tol does this device actually need?

Opened 2026-08-31. Kept as its own file rather than folded into
`hybrid_field_extraction_validation.md` because another session is actively
editing that doc, and because this result is device-methodology, not part of the
field-extraction plan it came out of.

**Motivation.** `decay_tol=1e-11` is inherited across this repo with **no
recorded justification anywhere**. The trigger was job 12627227, a fully-cold
Step 0 re-baseline that burned its whole 5 h wall without finishing on a device
whose previous run took 1:28.

**Method.** One run, not one per tolerance.
`FDTDSolver.solve(snapshot_every=10.0, snapshot_path=...)` decomposes the
monitors mid-run; because DFT monitors accumulate, the sample at time *t* is
exactly what a run stopped at *t* would have returned. Producer:
`examples/verification/python/mmi_1x2_decaytol_convergence.py` (job 12636328,
8 MPI ranks, res 40, 11 wavelengths). Reader:
`plot_mmi1x2_tolcurve.py`.

**This device has something the 4×4 does not: a converged reference.**
`resolved/mmi_1x2_validation.md` gap 5 records the dispersive full-device FDTD
result at `decay_tol=1e-11` (1.26 µm 0.44716, 1.31 µm 0.39527, 1.36 µm
0.34845), so the *actual error* can be measured instead of inferred from
per-step envelopes. On the 4×4 only the reference-free proxy is available.

## Result — the knee is at energy ratio ≈ 2e-5

| energy ratio | max \|T − ref(1e-11)\| | as % of T |
|---|---|---|
| 2.07e-03 | 0.01379 | 3.45% |
| 6.79e-04 | 0.00793 | 1.98% |
| 2.30e-04 | 0.00377 | 0.94% |
| 4.99e-05 | 0.00123 | 0.31% |
| **2.23e-05** | **0.00027** | **0.07%** |
| 1.20e-05 | 0.00044 | 0.11% |
| 5.72e-06 | 0.00030 | 0.08% |
| 8.89e-07 | 0.00012 | 0.03% |
| 3.82e-07 | 0.00006 | 0.02% |

**Below ratio ~2e-5 the error PLATEAUS at 0.02–0.1% and stops improving
proportionally** — and is not even monotonic (2.2e-5 → 0.07%, then 1.2e-5 →
0.11%), which is trapped-mode beating at the 2e-4 level. Buying 100× tighter
tolerance past the knee moves 0.07% → 0.02%.

## ⚠ THE COST SAVING IS ~2×, NOT ORDERS OF MAGNITUDE

Measured decay rate: **0.0409 decades of energy ratio per meep time unit**, so
each decade of tolerance costs a FIXED ~24 time units.

| decay_tol | t needed | vs 1e-5 |
|---|---|---|
| 1e-5 | 150 | — |
| 1e-9 | 263 (proj.) | 1.75× |
| 1e-11 | 312 (proj.) | 2.08× |

**Because decay is exponential, tolerance cost is LOGARITHMIC.** 1e-5 → 1e-11
is six decades but only ~2× the time-stepping. The common intuition that 1e-11
is astronomically expensive is **wrong on this device**.

**Consequence, stated plainly: loosening decay_tol was NOT the fix for the 5 h
wall.** Job 12627227 was the *hybrid*, dominated by `calibrate_body_neff`
(documented at ~97.6% of hybrid wall time) — a different code path from the
full-device solve measured here, which reached t=200 in 30 min. The tolerance
change buys ~1.75× on the full-device path and does not touch the calibration
cost that actually walled the job. **The hybrid's cost problem remains
unsolved and needs the calibration path profiled, not the tolerance.**

Where loosening genuinely matters is the 4×4, where the floor makes 1e-11
**unreachable** so the run never terminates — a categorically different failure
from "2× slower". See [[decay-tol-has-a-hard-floor]].

## ⚠⚠ CALIBRATION CAVEAT — the x-axis is the ENERGY criterion, not the shipped one

`snapshot_every`'s recorded `energy_ratio` is
`field_energy_in_box(total_volume) / peak_energy` (`fdtd.py:1443-1450`) — cell
energy relative to its own peak. That is **exactly
`mp.stop_when_energy_decayed`'s criterion**, and **NOT**
`mp.stop_when_dft_decayed`'s, which tests relative change in accumulated DFT
content. `fdtd.py`'s own comment says the ratio "maps each sample onto the
decay_tol that would have stopped there" — true for `decay_method="energy"`,
**indicative only for the shipped `"dft"` default.**

So the table above directly supports:

  * `decay_method="energy"`, `decay_tol ≈ 2e-5` → 0.07% error. **Supported.**
  * `decay_tol=1e-5` under the default dft criterion → **NOT established.**

It is also the footing the 4×4's floor measurements were taken on
(`decay_method="energy"`), so those numbers and this curve are directly
comparable.

**The calibration point comes for free from this same job**, which runs dft at
`decay_tol=1e-11`: reading the energy ratio at the moment it stops gives the
dft↔energy correspondence, measured rather than assumed. **Pending.**

## ✅ RUN COMPLETED — calibration point measured, and two of my own predictions falsified

Job 12636328 **COMPLETED** in **1:15:42** (8 ranks), 49 samples, stopping at
**t=490, energy ratio 1.056e-13**.

**The dft<->energy calibration, measured:**

> `decay_tol=1e-11` under the shipped `decay_method="dft"` stops at **energy
> ratio 1.06e-13**. The dft criterion is therefore **~100x stricter** than the
> energy criterion at the same nominal number — a ~2-decade offset on this
> device. Any snapshot `energy_ratio` axis must be shifted by this before being
> read as a dft `decay_tol`.

**Real cost of the shipped tolerance:** t=490 to stop, vs **t=150** for the
knee → **3.27x** the time-stepping for 0.07% -> ~0.02% accuracy. (My earlier
2.08x projection was made on the energy axis and was too low, precisely because
of the 100x offset above.)

**Two predictions I made from the easing signature, both FALSIFIED:**

1. I predicted `decay_tol=1e-11` would prove **unreachable** on this path. **It
   fired.** Reachable, in 1:15:42.
2. I predicted the energy would **turn around and grow** (the 4x4's
   floor-then-divergence). **It never did** — energy fell monotonically to the
   last sample; the run minimum IS the final point (2.988e-12 at t=490).

So the monotonic rate-easing (0.0324 -> ... -> 0.0016) and the flattening near
1.06e-13 were real, but they are **NOT** sufficient evidence of the 4x4's
instability. Here it was just a slow physical tail still converging. **The
easing signature is a warning worth watching, not a diagnosis; the 4x4's actual
signature is energy GROWTH.**

⚠ **Consequently my mid-run "correction" was itself wrong.** I briefly claimed
the walled Step 0 job (12627227) was chasing an unreachable target. It was not.
The original reading stands: that job is **calibration-dominated**, and
loosening `decay_tol` does not address it. The hybrid's cost problem still needs
the calibration path profiled.

## Status / what is NOT yet established

- The default in `mmi_1x2_hybrid.py` is currently **1e-9** with
  `minimum_run_time=100.0`. **Deliberately NOT moved to 1e-5 yet**, pending the
  dft↔energy calibration above — shipping 1e-5 on a mis-calibrated axis would
  be exactly the kind of inherited-without-justification choice this file
  exists to replace.
- Measured on the **full-device** `solve()` path only. The hybrid's
  `solve_section()` runs are shorter, different geometries with their own decay
  behaviour. Full-device is the harder case, so 2e-5 sufficing here is
  encouraging but is not the same measurement.
- Band edges converge LAST: at ratio 2.3e-4, 1.31 µm was already within
  0.0004 while 1.36 µm was still 0.0038 out and climbing monotonically. Reading
  band centre alone would have called convergence an order of magnitude early —
  the same trap `plot_mmi4x4_90hybrid_tolcurve.py` documents, now confirmed
  independently on this device. This also corroborates gap 5's finding that
  this device's worst point is 1.36 µm, from a different direction.
- Splitter symmetry stayed at 2.2e-07 throughout, so the run is healthy.
