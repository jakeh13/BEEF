"""Plot the 1x2's live transmission-vs-decay-time curve as the run produces it.

validation_markdown/hybrid_field_extraction_validation.md, 2026-08-31. Reader
for the CSV mmi_1x2_decaytol_convergence.py appends to. Direct port of
plot_mmi4x4_90hybrid_tolcurve.py -- INCLUDING the three methodological traps
that script paid for, all of which apply here and two of which apply MORE
strongly because this device carries 11 wavelengths, not 3:

 1. **Per-step |dT|, never |T(t) - T_final|.** On a live curve "T_final" is
    merely the newest sample, so that residual collapses to zero at the
    right-hand edge BY CONSTRUCTION and everything looks converged regardless
    of the physics. The per-step change needs no final value and is honest at
    every moment of the run.

 2. **Worst case over WAVELENGTHS, not band centre.** Reading only band centre
    on the 4x4 (2026-08-24) produced a confidently wrong "converged by ratio
    5.7e-2" call: band centre converges FASTEST while the band EDGE was still
    0.43 dB from settled. This device has an independent reason to expect the
    same -- resolved/mmi_1x2_validation.md gap 5 puts its own worst
    hybrid-vs-FDTD point at 1.36 um, the band edge -- so the envelope below is
    the max over all 11 wavelengths, and the band-centre-only curve is drawn
    faintly alongside purely to show how much more optimistic it looks.

 3. **5-sample envelope, not 3.** Trapped-mode beating has a long period; a
    3-sample window can sit inside one lobe and read as converged when it is
    not.

Safe to run repeatedly while the job is in flight: read-only, tolerates a
partially-written trailing line, and skips a file with too few samples.

Run:
    python plot_mmi1x2_tolcurve.py
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                   "example_output", "mmi_1x2"))
CSV_NAME = os.environ.get("MMI1X2_TOLCURVE_CSV", "tolcurve_res40_top.csv")
WLS = np.linspace(1.26, 1.36, 11)
CENTRE_WL_IDX = 5          # 1.31 um


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
            good.append({k: float(v) for k, v in r.items()
                         if v not in (None, "")})
        except (TypeError, ValueError):
            continue          # partial trailing write -- ignore, not an error
    if len(good) < 3:
        return None
    keys = [k for k in good[0] if k.startswith("T_out")]
    return dict(
        t=np.array([g["t"] for g in good]),
        ratio=np.array([g["energy_ratio"] for g in good]),
        energy=np.array([g["energy"] for g in good]),
        T={k: np.array([g[k] for g in good]) for k in sorted(keys)},
    )


def report(d):
    """Text summary -- the actual deliverable; the PNG is for the record."""
    ports = sorted({k.rsplit("_wl", 1)[0] for k in d["T"]})
    print(f"samples: {len(d['t'])}   t: {d['t'][0]:.0f} -> {d['t'][-1]:.0f}   "
          f"energy ratio: {d['ratio'][0]:.3e} -> {d['ratio'][-1]:.3e}")
    ic = CENTRE_WL_IDX
    for p in ports:
        c = f"{p}_wl{ic}"
        if c in d["T"]:
            print(f"  {p} @1.31um: {d['T'][c][0]:.6f} -> {d['T'][c][-1]:.6f}")

    print("\n  where does the WORST wavelength stop moving?")
    print(f"  {'tol needed':>12} {'= decay_tol':>13}  (worst per-step |dT| over "
          f"all 11 wl, 5-sample envelope)")
    for thresh, lbl in ((1e-2, "0.1 dB-ish"), (3e-3, "3x tighter"),
                        (1e-3, "10x tighter"), (1e-4, "100x tighter")):
        worst_ratio = None
        for p in ports:
            cols = [k for k in d["T"] if k.rsplit("_wl", 1)[0] == p]
            step = np.max([np.abs(np.diff(d["T"][c])) for c in cols], axis=0)
            env = np.array([step[max(0, i - 4):i + 1].max()
                            for i in range(len(step))])
            rr = d["ratio"][1:]
            ok = np.where((env < thresh) & (rr > 0))[0]
            # First sample where the envelope is under threshold AND stays
            # under it for the rest of the record -- a single dip below is
            # exactly the trapped-mode-beating artifact lesson 3 warns about.
            first = None
            for i in ok:
                if np.all(env[i:] < thresh):
                    first = i
                    break
            if first is None:
                worst_ratio = None
                break
            r = rr[first]
            worst_ratio = r if worst_ratio is None else min(worst_ratio, r)
        if worst_ratio is None:
            print(f"  {lbl:>12} {'NOT YET REACHED':>13}")
        else:
            print(f"  {lbl:>12} {worst_ratio:13.2e}")


def main():
    path = os.path.join(OUT, CSV_NAME)
    d = load(path)
    if d is None:
        print(f"no curve data with >=3 samples yet at {path}")
        return
    report(d)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))
    ports = sorted({k.rsplit("_wl", 1)[0] for k in d["T"]})

    ax = axes[0]
    for p in ports:
        c = f"{p}_wl{CENTRE_WL_IDX}"
        if c in d["T"]:
            ax.plot(d["t"], d["T"][c], lw=1.5, label=p)
    ax.set_xlabel("meep time")
    ax.set_ylabel("|T|²  (1.31 µm, band centre)")
    ax.set_title(f"1×2 MMI — transmission vs time  (n={len(d['t'])})")
    ax.grid(alpha=.3); ax.legend(fontsize=8)

    ax = axes[1]
    for p in ports:
        cols = [k for k in d["T"] if k.rsplit("_wl", 1)[0] == p]
        step = np.max([np.abs(np.diff(d["T"][c])) for c in cols], axis=0)
        rr = d["ratio"][1:]
        m = (rr > 0) & (step > 0)
        if m.sum() < 2:
            continue
        line, = ax.loglog(rr[m], step[m], lw=0.8, alpha=.3)
        env = np.array([step[max(0, i - 4):i + 1].max()
                        for i in range(len(step))])
        me = (rr > 0) & (env > 0)
        ax.loglog(rr[me], env[me], lw=2.0, color=line.get_color(),
                  label=f"{p} (worst of 11 wl)")
        # Band-centre-only, faint: shows how much more optimistic the
        # single-wavelength read looks. See lesson 2.
        cc = f"{p}_wl{CENTRE_WL_IDX}"
        if cc in d["T"]:
            s1 = np.abs(np.diff(d["T"][cc]))
            m1 = (rr > 0) & (s1 > 0)
            ax.loglog(rr[m1], s1[m1], lw=0.9, ls="--", alpha=.55,
                      color=line.get_color(), label=f"{p} (1.31 µm only)")
    for lvl, txt in ((1e-2, "0.1 dB-ish"), (1e-3, "10x tighter")):
        ax.axhline(lvl, color="k", ls=":", lw=.8, alpha=.6)
        ax.text(ax.get_xlim()[0], lvl, " " + txt, va="bottom", fontsize=7)
    ax.invert_xaxis()
    ax.set_xlabel("energy ratio  ( = the decay_tol that would stop here )")
    ax.set_ylabel("per-step |ΔT|  (thick = worst wl, 5-sample env)")
    ax.set_title("1×2 MMI — still moving?  (dashed = band centre alone,\n"
                 "which converges FASTEST and reads too optimistic)")
    ax.grid(alpha=.3, which="both"); ax.legend(fontsize=7)

    fig.suptitle("1×2 MMI — transmission vs decay time (live)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(OUT, "tolcurve_res40.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"\nplot -> {out}")


if __name__ == "__main__":
    main()
