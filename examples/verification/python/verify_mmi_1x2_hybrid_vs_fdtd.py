"""mmi_1x2_validation.md gap 5 — the only check that tests the hybrid
cascade against a full-device reference, as a whole rather than piece by
piece (README.md check 4).

Sibling of verify_mmi_2x2_hybrid_vs_fdtd.py (mmi_2x2_validation.md gap 5) --
same comparison, same reasoning, ported to this device's port names and
paths. `mmi_1x2_fdtd.py` (full-device FDTD, the reference) and
`mmi_1x2_hybrid.py` (2 FDTD sections + 1 EME body, cascaded) both exist and
both save an SMatrix at the same wavelength grid (WL_MIN/MAX/N_WL =
1.26/1.36/11), but nothing ever loaded the two and compared them. This one
does. Runs neither simulation itself -- re-run `make mmi_1x2_fdtd
mmi_1x2_hybrid` (batch/run_mmi_1x2_examples.batch array indices 1-2) first if
either output is stale.

Port naming differs between the two scripts and is reconciled here:
`mmi_1x2_fdtd.py` labels its output_monitor_centers-ordered ports "out_0"
(top) / "out_1" (bottom); `mmi_1x2_hybrid.py` labels the same two physical
outputs "out_top" / "out_bot" directly. Both excite the (lone) input at
in_mode=1/out_mode=1.

Measurement and report, not a pass/fail gate -- see the module docstring of
verify_mmi_2x2_hybrid_vs_fdtd.py for why: gap 4 (engine agreement) is not yet
measured for this device either, and this device already has a real,
understood basis-mismatch history of its own (gap 3, N_MODES 12->8, fixed
2026-07-30) that a hard threshold here could be mistaken for. Unlike the
2x2, this device already has three goldens (mmi_1x2_validation.md's "What it
already has") -- this check is complementary to those, not a replacement:
the goldens pin each script's own output against itself over time, this
compares the two scripts' outputs against EACH OTHER right now.

Carries provenance in the output file (git sha, mtimes, config), per
README.md's complaint about the coupler's original un-provenanced check-4
reference.

Run (after both example scripts have produced their SMatrix pickles):
    python examples/verification/python/verify_mmi_1x2_hybrid_vs_fdtd.py
"""
from __future__ import annotations

import json
import os
import platform
import subprocess
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import beef
from beef import SMatrix

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
FDTD_PKL = os.path.join(_ROOT, "example_output", "mmi_1x2", "mmi_1x2_fdtd", "smatrix.pkl")

# CALIBRATED points this at mmi_1x2_hybrid.py's MMI1X2_CALIBRATE_BODY=1 output
# instead of the default uncalibrated one -- gap 5's follow-up, testing whether
# HybridSolver(calibrate_body_neff=True) closes the measured 9.11% the way it
# did for the coupler (~5% -> 0.11%, resolved/directional_coupler_validation.md
# items 6/14). FDTD_PKL is deliberately NOT redirected: it is the shared
# full-device reference both the calibrated and uncalibrated hybrid are
# measured against, same reasoning as directional_coupler_hybrid_stub_sweep.py
# not redirecting _FDTD_SWEEP_NPZ.
CALIBRATED = bool(int(os.environ.get("MMI1X2_HYBRID_CALIBRATED", "0")))
_hybrid_dir = "mmi_1x2_hybrid_calibrated" if CALIBRATED else "mmi_1x2_hybrid"
HYBRID_PKL = os.path.join(_ROOT, "example_output", "mmi_1x2", _hybrid_dir, "smatrix_hybrid.pkl")
_out_name = "mmi_1x2_hybrid_vs_fdtd_calibrated" if CALIBRATED else "mmi_1x2_hybrid_vs_fdtd"
output_dir = os.path.join(_ROOT, "example_output", "verification", _out_name)

# (fdtd port, hybrid port) pairs -- both physical outputs of the single-input device.
PORT_PAIRS = [("out_0", "out_top"), ("out_1", "out_bot")]
PORT_LABELS = {"out_0": "top", "out_1": "bottom"}


def _git_sha(repo_dir):
    try:
        sha = subprocess.check_output(
            ["git", "-C", repo_dir, "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
        dirty = subprocess.check_output(
            ["git", "-C", repo_dir, "status", "--porcelain"],
            text=True, stderr=subprocess.DEVNULL).strip()
        return f"{sha}{'+dirty' if dirty else ''}"
    except Exception:
        return "<not a git checkout>"


def provenance(smat_fdtd, smat_hybrid):
    return dict(
        timestamp_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        beef_version=getattr(beef, "__version__", "<none>"),
        beef_git=_git_sha(os.path.join(_ROOT, "BEEF")),
        python=platform.python_version(),
        host=platform.node(),
        fdtd_pkl=FDTD_PKL,
        hybrid_pkl=HYBRID_PKL,
        fdtd_pkl_mtime=datetime.fromtimestamp(
            os.path.getmtime(FDTD_PKL), tz=timezone.utc).isoformat(timespec="seconds"),
        hybrid_pkl_mtime=datetime.fromtimestamp(
            os.path.getmtime(HYBRID_PKL), tz=timezone.utc).isoformat(timespec="seconds"),
        fdtd_wavelengths=list(map(float, smat_fdtd.wavelengths)),
        hybrid_wavelengths=list(map(float, smat_hybrid.wavelengths)),
        # Pinned constants copied from each script -- not re-derived. Both
        # currently agree: RESOLUTION=40, WL_MIN/MAX/N_WL=1.26/1.36/11,
        # N_MODES=8 (mmi_1x2_hybrid.py only, post mmi_1x2_validation.md gap 3).
        note="see mmi_1x2_fdtd.py / mmi_1x2_hybrid.py for full solver config.",
    )


def main():
    print(f"CALIBRATED = {CALIBRATED}  (hybrid pkl -> {HYBRID_PKL})")
    for path, name in [(FDTD_PKL, "mmi_1x2_fdtd.py"), (HYBRID_PKL, "mmi_1x2_hybrid.py")]:
        if not os.path.exists(path):
            hint = (
                "MMI1X2_CALIBRATE_BODY=1 python mmi_1x2_hybrid.py"
                if CALIBRATED and path == HYBRID_PKL else
                f"`{name}` (e.g. `make mmi_1x2_fdtd mmi_1x2_hybrid` from "
                f"examples/, or batch/run_mmi_1x2_examples.batch array indices "
                f"1-2)"
            )
            raise SystemExit(
                f"verify_mmi_1x2_hybrid_vs_fdtd: {path} does not exist -- run "
                f"{hint} first."
            )

    smat_fdtd = SMatrix.load(FDTD_PKL)
    smat_hybrid = SMatrix.load(HYBRID_PKL)

    wl_fdtd, wl_hybrid = smat_fdtd.wavelengths, smat_hybrid.wavelengths
    if wl_fdtd.shape != wl_hybrid.shape or not np.allclose(wl_fdtd, wl_hybrid):
        raise SystemExit(
            f"verify_mmi_1x2_hybrid_vs_fdtd: wavelength grids differ between "
            f"the two saved SMatrix files ({wl_fdtd} vs {wl_hybrid}) -- one of "
            f"the two example scripts' WL_MIN/WL_MAX/N_WL has drifted from the "
            f"other. Re-run both from the current source before comparing."
        )
    wavelengths = wl_fdtd

    print(f"Loaded {FDTD_PKL}\n  ports={smat_fdtd.ports()}")
    print(f"Loaded {HYBRID_PKL}\n  ports={smat_hybrid.ports()}\n")

    passive_fdtd = smat_fdtd.check_passive()
    passive_hybrid = smat_hybrid.check_passive()
    print(f"check_passive(): fdtd={passive_fdtd}, hybrid={passive_hybrid}\n")

    results = {}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for (fdtd_port, hybrid_port), ax in zip(PORT_PAIRS, axes):
        c_fdtd = smat_fdtd.coeff(port=fdtd_port, in_mode=1, out_mode=1)
        c_hybrid = smat_hybrid.coeff(port=hybrid_port, in_mode=1, out_mode=1)
        T_fdtd = np.abs(c_fdtd) ** 2
        T_hybrid = np.abs(c_hybrid) ** 2

        dT = T_hybrid - T_fdtd
        max_abs_dT = float(np.max(np.abs(dT)))
        mean_abs_dT = float(np.mean(np.abs(dT)))
        worst_idx = int(np.argmax(np.abs(dT)))
        rel_dT = max_abs_dT / max(float(np.mean(T_fdtd)), 1e-6)

        label = PORT_LABELS[fdtd_port]
        print(f"-- {label} : fdtd '{fdtd_port}' vs hybrid '{hybrid_port}' --")
        print(f"  {'lambda':>8} {'T_fdtd':>10} {'T_hybrid':>10} {'dT':>10}")
        for j in range(len(wavelengths)):
            flag = "  <-- worst" if j == worst_idx else ""
            print(f"  {wavelengths[j]:8.4f} {T_fdtd[j]:10.4f} {T_hybrid[j]:10.4f} "
                  f"{dT[j]:+10.4f}{flag}")
        print(f"  max|dT|={max_abs_dT:.4f}  mean|dT|={mean_abs_dT:.4f}  "
              f"max|dT|/mean(T_fdtd)={rel_dT:.2%}  "
              f"at lambda={wavelengths[worst_idx]:.4f} um\n")

        results[fdtd_port] = dict(
            label=label, T_fdtd=T_fdtd, T_hybrid=T_hybrid, dT=dT,
            max_abs_dT=max_abs_dT, mean_abs_dT=mean_abs_dT, rel_dT=rel_dT,
        )

        ax.plot(wavelengths, T_fdtd, label="Full-device FDTD (reference)",
                color="k", lw=1.5, marker="o", markersize=3)
        ax.plot(wavelengths, T_hybrid, label="Hybrid FDTD/EME",
                color="steelblue", lw=1.2, linestyle="--", marker="s", markersize=3)
        ax.set_xlabel("Wavelength (µm)")
        ax.set_ylabel("Transmission")
        ax.set_title(f"{label} output: hybrid vs full-device FDTD")
        ax.legend(fontsize=8)
        ax.grid(True)
        ax.set_ylim(0, 1.05)
    _cal_tag = " (calibrate_body_neff=True)" if CALIBRATED else ""
    plt.suptitle(f"1×2 MMI — Hybrid{_cal_tag} vs Full-Device FDTD "
                f"(mmi_1x2_validation.md gap 5)", fontsize=11)
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, "hybrid_vs_fdtd.png"), dpi=300)
    plt.close("all")

    worst_overall = max(r["max_abs_dT"] for r in results.values())
    print(f"Worst max|dT| across both ports: {worst_overall:.4f}")
    if CALIBRATED:
        print(
            "Calibrated run. Compare this max|dT| against the uncalibrated "
            "baseline (0.0359, 9.11% of mean T_fdtd, every dT negative -- "
            "mmi_1x2_validation.md gap 5) to see whether calibrate_body_neff "
            "closes the gap the way it did for the coupler (~5% -> 0.11%, "
            "resolved/directional_coupler_validation.md items 6/14). No "
            "threshold enforced here either -- gap 4 (engine agreement) is "
            "still not measured for this device, so there is still no basis "
            "to call a residual calibratable vs. structural, only to compare "
            "it against the uncalibrated number."
        )
    else:
        print(
            "Informational only -- no threshold enforced here. "
            "mmi_1x2_validation.md gap 4 (engine agreement) is not yet measured "
            "for this device, so there is no basis yet to say how much of this "
            "gap (if any) is calibratable vs. structural. See "
            "resolved/directional_coupler_validation.md items 6/14 for how that "
            "played out on the coupler (uncalibrated hybrid ~5% off; calibrated "
            "0.11% off). Re-run with MMI1X2_HYBRID_CALIBRATED=1 after generating "
            "the calibrated hybrid output to test that directly on this device."
        )

    prov = provenance(smat_fdtd, smat_hybrid)
    prov["ports"] = {
        fdtd_port: {k: v for k, v in r.items() if k not in ("T_fdtd", "T_hybrid", "dT")}
        for fdtd_port, r in results.items()
    }
    prov_json = json.dumps(prov, indent=2, sort_keys=True, default=list)
    with open(os.path.join(output_dir, "provenance.json"), "w") as f:
        f.write(prov_json + "\n")
    np.savez(
        os.path.join(output_dir, "comparison.npz"),
        wavelengths=wavelengths,
        **{f"T_fdtd_{p}": results[p]["T_fdtd"] for p in results},
        **{f"T_hybrid_{p}": results[p]["T_hybrid"] for p in results},
        provenance=np.array(prov_json),
    )
    print(f"\nSaved -> {output_dir}/hybrid_vs_fdtd.png, comparison.npz, provenance.json")


if __name__ == "__main__":
    main()
