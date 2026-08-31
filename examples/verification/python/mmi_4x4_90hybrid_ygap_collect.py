"""Collect the (mm_width, y_gap) campaign into one surface + figure.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-27.

Reads whatever eval_W*_y*.npz files exist and reports the objective surface,
so it is useful DURING the campaign and not only at the end -- points that
have not run yet simply show as gaps.

2026-08-28: rewritten to scan the output directory directly rather than
iterate the 106-point grid's own y_gaps_for(w) lattice. The campaign is now
BO-driven (mmi_4x4_90hybrid_bo.py), which searches y_gap CONTINUOUSLY per
width -- a real point landed at (7.00, 1.36), off the 0.05 um lattice, and
the old grid-enumeration version silently never displayed it (load_all()
only ever looked up POINTS' own (w, y) tuples). Widths are still grouped by
MM_WIDTHS since BO restricts that axis to the same categorical lattice (see
bo.py's own reasoning: body solves are only cached at those 5 widths), but
y_gap is read from whatever files actually exist, not from a fixed grid.

WHAT TO DISTRUST IN THE OUTPUT. Two edges are reported explicitly because a
minimum sitting on either one is not an optimum, it is a truncation:

  * best_mm_length at LENGTH_MIN/LENGTH_MAX (20-400 um, mmi_4x4_90hybrid_
    ygap_grid.py) -- widened 2026-08-28 from the original +-10 um window
    around 3/4 L_pi(W) to an absolute range spanning multiple self-imaging
    periods, since solve_section_length_sweep's cost is independent of how
    many/which lengths are swept (one mode solve total). If a point's best
    length still lands on 20 or 400, that range itself needs widening -- the
    reported optimum there is a boundary value, not a stationary point.
  * y_gap at the grid's low edge (0.10) -- that bound is a modelling choice,
    not physics (see the grid module), so an optimum resting on it means the
    range should be extended rather than that 0.10 is best.

Run:
    python mmi_4x4_90hybrid_ygap_collect.py
"""
import os

import numpy as np

from mmi_4x4_90hybrid_ygap_grid import MM_WIDTHS, LENGTH_MIN, LENGTH_MAX

_HERE = os.path.dirname(os.path.abspath(__file__))
out_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                       "example_output", "mmi_4x4_90hybrid",
                                       "wg_ygap_opt"))
NOMINAL = (6.0, 0.70)


def load_all():
    got = {}
    if not os.path.isdir(out_dir):
        return got
    for fn in sorted(os.listdir(out_dir)):
        if not (fn.startswith("eval_W") and fn.endswith(".npz")):
            continue
        try:
            w = float(fn[6:10])
            y = float(fn[12:16])
        except ValueError:
            continue
        with np.load(os.path.join(out_dir, fn)) as d:
            got[(w, y)] = {k: d[k] for k in d.files}
    return got


def widths_present(got):
    """MM_WIDTHS first (BO's own categorical axis), then any off-lattice
    width found on disk (e.g. a manual/grid-era run), so nothing is dropped."""
    extra = sorted(set(w for w, _ in got) - set(MM_WIDTHS))
    return list(MM_WIDTHS) + extra


def main():
    got = load_all()
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit(f"COLLECT: {len(got)} point(s) present  ({out_dir})\n")
    if not got:
        _write_summary(lines)
        return

    emit(f"{'W':>5} {'y_gap':>6} {'taper':>6} {'dL':>7} {'mm_len':>9} "
         f"{'max|dphi|':>10} {'minCMRR':>9} {'obj':>8}")
    emit("-" * 68)
    for w in widths_present(got):
        ys = sorted(y for (ww, y) in got if ww == w)
        for y in ys:
            d = got[(w, y)]
            dl = float(d["best_delta_L"])
            mm_len = float(d["best_mm_length"])
            flag = " <-EDGE" if (abs(mm_len - LENGTH_MIN) < 1e-6
                                  or abs(mm_len - LENGTH_MAX) < 1e-6) else ""
            flag += " <-NOMINAL" if (w, y) == NOMINAL else ""
            emit(f"{w:5.1f} {y:6.2f} {float(d['taper_width']):6.3f} "
                 f"{dl:+7.2f} {float(d['best_mm_length']):9.3f} "
                 f"{float(d['max_dphi']):10.3f} {float(d['min_cmrr']):9.3f} "
                 f"{float(d['objective']):8.4f}{flag}")
        emit()

    best = min(got.items(), key=lambda kv: float(kv[1]["objective"]))
    (bw, by), bd = best
    emit("=== best so far ===")
    emit(f"  mm_width={bw:.2f}  y_gap={by:.2f}  taper_width={float(bd['taper_width']):.3f}")
    emit(f"  mm_length={float(bd['best_mm_length']):.4f} "
         f"(dL={float(bd['best_delta_L']):+.2f})")
    emit(f"  max|dphi|={float(bd['max_dphi']):.3f} deg   "
         f"min CMRR={float(bd['min_cmrr']):.3f} dBe   "
         f"objective={float(bd['objective']):.4f}")
    emit("  spec: |dphi| <= 5 deg, CMRR >= 20 dBe")

    if NOMINAL in got:
        n = got[NOMINAL]
        emit(f"\n  nominal (6.0, 0.70): objective={float(n['objective']):.4f} "
             f"max|dphi|={float(n['max_dphi']):.3f} "
             f"minCMRR={float(n['min_cmrr']):.3f}")

    edges = [k for k, d in got.items()
             if abs(float(d["best_mm_length"]) - LENGTH_MIN) < 1e-6
             or abs(float(d["best_mm_length"]) - LENGTH_MAX) < 1e-6]
    if edges:
        emit(f"\n  WARNING: {len(edges)} point(s) hit the {LENGTH_MIN}-{LENGTH_MAX} um "
             f"length-sweep edge -- widen MM_LENGTH_GRID at those widths before "
             f"trusting them:")
        for k in sorted(edges):
            emit(f"    W={k[0]:.1f} y={k[1]:.2f}")
    lo = [k for k in got if abs(k[1] - 0.10) < 1e-9]
    if lo and (by - 0.10) < 1e-9:
        emit("\n  WARNING: best point sits on the y_gap LOWER bound (0.10), which "
             "is a modelling choice -- extend the range rather than accepting it.")

    _write_summary(lines)
    plot(got)


def _write_summary(lines):
    # Persisted so the campaign's state survives outside terminal scrollback --
    # this file is what a periodic re-run refreshes and what a live viewer
    # tails, rather than the printed table which vanishes with the shell.
    p = os.path.join(out_dir, "ygap_campaign_summary.txt")
    with open(p, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nCOLLECT: summary -> {p}")


def plot(got):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"\nCOLLECT: no plot ({exc})")
        return
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for ax, key, lab in zip(
            axes, ["objective", "max_dphi", "min_cmrr"],
            ["objective (lower better)", "max |dphi| (deg)", "min CMRR (dBe)"]):
        for w in widths_present(got):
            ys = sorted(y for (ww, y) in got if ww == w)
            if not ys:
                continue
            ax.plot(ys, [float(got[(w, y)][key]) for y in ys],
                    marker="o", ms=3, label=f"W={w:.1f}")
        ax.set_xlabel("y_gap (um)")
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    axes[1].axhline(5.0, color="k", ls="--", lw=1, label="spec")
    axes[2].axhline(20.0, color="k", ls="--", lw=1, label="spec")
    axes[0].legend(fontsize=8)
    fig.suptitle("4x4 90-deg hybrid: hybrid-path (FDTD inputs + EME body) "
                 "y_gap campaign, 1.29-1.33 um")
    fig.tight_layout()
    p = os.path.join(out_dir, "ygap_campaign_surface.png")
    fig.savefig(p, dpi=140)
    print(f"\nCOLLECT: plot -> {p}")


if __name__ == "__main__":
    main()
