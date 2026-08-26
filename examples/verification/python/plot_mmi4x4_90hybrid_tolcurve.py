"""Plot the live transmission-vs-decay-time curves as the 6b runs produce them.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-24.

Reads the CSVs `FDTDSolver.solve(snapshot_path=...)` appends to while a run is
still in progress and regenerates a PNG. Safe to call repeatedly: it only
reads, tolerates a partially-written trailing line (the solver flushes per
sample, but a read can still land mid-write), and skips a file that has too
few samples to say anything yet.

Kept OUT of the solver deliberately. Plotting inside the snapshot callback
would put matplotlib on the hot path of every MPI rank, and the CSV is already
the source of truth -- a separate reader can regenerate the figure as often as
wanted without touching the simulation.

TWO PANELS, because they answer different questions:
  * left  -- T per output port vs meep time. What the device is doing.
  * right -- per-step |dT| vs energy ratio, with a 3-sample envelope. THE
             decision panel, and deliberately REFERENCE-FREE: the x-axis
             is exactly the `decay_tol` that would have stopped the run there,
             so reading off where the curve drops below the precision the
             device's conclusions need (imbalance ~0.1 dB, phase ~1 deg) gives
             the tolerance that is actually required, measured rather than
             inherited.

Run:
    python plot_mmi4x4_90hybrid_tolcurve.py
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                   "example_output", "mmi_4x4_90hybrid"))
SOURCES = [("tolcurve_res40_bot.csv", "port 1 (input_idx=0)"),
           ("tolcurve_res40_idx2.csv", "port 3 (input_idx=2)")]
WL_LABEL = ["1.26 µm", "1.31 µm", "1.36 µm"]


def load(path):
    """Rows as floats; drops any short/partial trailing line."""
    if not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
    except OSError:
        return None
    good = []
    for r in rows:
        try:
            good.append({k: float(v) for k, v in r.items() if v not in (None, "")})
        except (TypeError, ValueError):
            continue          # partial trailing write -- ignore, not an error
    if len(good) < 3:
        return None
    keys = [k for k in good[0] if k.startswith("T_out")]
    return dict(
        t=np.array([g["t"] for g in good]),
        ratio=np.array([g["energy_ratio"] for g in good]),
        T={k: np.array([g[k] for g in good]) for k in sorted(keys)},
    )


def main():
    data = [(lbl, load(os.path.join(OUT, fn))) for fn, lbl in SOURCES]
    data = [(lbl, d) for lbl, d in data if d is not None]
    if not data:
        print("no curve data with >=3 samples yet")
        return

    fig, axes = plt.subplots(len(data), 2, figsize=(13, 4.6 * len(data)),
                             squeeze=False)
    for row, (label, d) in enumerate(data):
        ax = axes[row][0]
        # Band-centre wavelength (index 1 of 3) keeps the panel readable.
        for k in d["T"]:
            if k.endswith("_wl1"):
                ax.plot(d["t"], d["T"][k], lw=1.4, label=k.replace("_wl1", ""))
        ax.set_xlabel("meep time"); ax.set_ylabel("|T|²  (1.31 µm, band centre)")
        ax.set_title(f"{label} — transmission vs time  (n={len(d['t'])})")
        ax.grid(alpha=.3); ax.legend(fontsize=8, ncol=2)

        ax = axes[row][1]
        # REFERENCE-FREE convergence measure. The obvious choice,
        # |T(t) - T_final|, is INVALID on a live plot: while the run is in
        # flight "T_final" is only the most recent sample, so the residual
        # collapses to zero at the right-hand edge by construction and the
        # curve looks converged no matter what the physics is doing. Plotting
        # the per-step change instead needs no final value and is honest at
        # every moment of the run.
        #
        # Per-step differences are noisy here because the trapped modes beat
        # against each other (transmission oscillates rather than settling
        # monotonically), so the running max over a 3-sample window is drawn
        # heavier -- that envelope, not the individual steps, is what should
        # be read against the precision lines.
        # WORST CASE OVER WAVELENGTHS, not band centre. Reading only _wl1
        # (1.31 um) on 2026-08-24 produced a confidently wrong "converged by
        # ratio 5.7e-2" call: band centre converges FASTEST, while the worst
        # wavelength is 1.26 um -- the band edge, which is exactly where this
        # device's CMRR and phase failures concentrate. At t=1100 the band
        # edge was still 0.43 dB from settled while band centre looked done.
        ports = sorted({k.rsplit("_wl", 1)[0] for k in d["T"]})
        for port in ports:
            cols = [k for k in d["T"] if k.rsplit("_wl", 1)[0] == port]
            step = np.max([np.abs(np.diff(d["T"][c])) for c in cols], axis=0)
            rr = d["ratio"][1:]
            m = (rr > 0) & (step > 0)
            if m.sum() < 2:
                continue
            line, = ax.loglog(rr[m], step[m], lw=0.8, alpha=.35)
            # 5-sample envelope, not 3: the trapped-mode beating has a long
            # period, and a short window sits inside one lobe and reads as
            # converged when it is not.
            env = np.array([step[max(0, i - 4):i + 1].max()
                            for i in range(len(step))])
            me = (rr > 0) & (env > 0)
            ax.loglog(rr[me], env[me], lw=1.8, color=line.get_color(),
                      label=port)
        for lvl, txt in ((1e-2, "0.1 dB-ish"), (1e-3, "10x tighter")):
            ax.axhline(lvl, color="k", ls=":", lw=.8, alpha=.6)
            ax.text(ax.get_xlim()[0], lvl, " " + txt, va="bottom", fontsize=7)
        ax.invert_xaxis()      # run proceeds right-to-left as energy decays
        ax.set_xlabel("energy ratio  ( = the decay_tol that would stop here )")
        ax.set_ylabel("per-step |ΔT|, worst wavelength  (thick = 5-sample env)")
        ax.set_title(f"{label} — still moving (worst of 1.26/1.31/1.36 µm)")
        ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8, ncol=2)

    fig.suptitle("4×4 90° hybrid — transmission vs decay time (live)", y=1.0)
    fig.tight_layout()
    dest = os.path.join(OUT, "tolcurve_live.png")
    fig.savefig(dest, dpi=140, bbox_inches="tight")
    print(f"wrote {dest}")
    for label, d in data:
        print(f"  {label}: n={len(d['t'])} t={d['t'][-1]:.0f} "
              f"ratio={d['ratio'][-1]:.2e}")


if __name__ == "__main__":
    main()
