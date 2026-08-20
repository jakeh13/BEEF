"""mmi_2x2_validation.md gap 5 — the only check that tests the hybrid
cascade against a full-device reference, as a whole rather than piece by
piece (README.md check 4).

Nothing asserted this before: `mmi_2x2_fdtd.py` (full-device FDTD, the
reference) and `mmi_2x2_hybrid.py` (2 FDTD sections + 1 EME body, cascaded)
both exist and both save an SMatrix, but no script ever loaded the two and
compared them. This one does, at the two scripts' own pinned configuration
(RESOLUTION=40, same wavelength grid, top port excited) -- it runs neither
simulation itself, so a change to either example's constants is picked up
automatically the next time both are (re-)run via `make mmi_2x2_fdtd
mmi_2x2_hybrid` (see batch/run_mmi_2x2_examples.batch, array indices 1-2).

Port naming differs between the two scripts and is reconciled here:
`mmi_2x2_fdtd.py` labels its output_monitor_centers-ordered ports "out_0"
(top half, bar) / "out_1" (bottom half, cross); `mmi_2x2_hybrid.py` labels
the same two physical outputs "bar" / "cross" directly. Both excite the top
input at in_mode=1/out_mode=1 (see each script's own step-5 comment), so no
further disambiguation is needed.

This is measurement and report, not a pass/fail gate: `mmi_2x2_validation.md`
gap 4 (engine agreement -- whether meep's and meow's MPB solves agree on
phase through the shared uniform body) has not been measured for this
device, and the directional coupler's own history
(`resolved/directional_coupler_validation.md` items 6/14) is that an
uncalibrated hybrid can disagree with a full-device reference by several
percent for a reason that is real physics (a time-stepper vs. an eigenmode
solver), not a bug -- and that calibrating it (`HybridSolver
(calibrate_body_neff=True)`) is exactly what closes gap 4, done separately.
Asserting a threshold here, before gap 4 exists for this device, would
either mask a real finding or fail on an expected, already-understood gap.
Report the numbers; the follow-up is gap 4, not a tighter tolerance here.

Carries provenance in the output file, unlike the coupler's original
check-4 reference (README.md's complaint: "no config, no solver fingerprint,
no geometry recorded").

Run (after both example scripts have produced their SMatrix pickles):
    python examples/verification/python/verify_mmi_2x2_hybrid_vs_fdtd.py
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
from beef import SMatrix, interp_complex

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
FDTD_PKL = os.path.join(_ROOT, "example_output", "mmi_2x2", "mmi_2x2_fdtd", "smatrix.pkl")
HYBRID_PKL = os.path.join(_ROOT, "example_output", "mmi_2x2", "mmi_2x2_hybrid", "smatrix_hybrid.pkl")
output_dir = os.path.join(_ROOT, "example_output", "verification", "mmi_2x2_hybrid_vs_fdtd")

# (fdtd port, hybrid port) pairs -- both physical outputs of the top-excited device.
PORT_PAIRS = [("out_0", "bar"), ("out_1", "cross")]
PORT_LABELS = {"out_0": "bar (top)", "out_1": "cross (bottom)"}


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
        # Pinned constants copied from each script -- not re-derived, so this
        # provenance is only as good as keeping the two in sync by hand. Both
        # currently agree: RESOLUTION=40, WL_MIN/MAX/N_WL=1.530/1.565/11,
        # N_MODES=14 (mmi_2x2_hybrid.py only, post fdtd_vs_EME_mode_validation.md
        # item 9), body_stub=0.25, mm_length=28.0.
        note="see mmi_2x2_fdtd.py / mmi_2x2_hybrid.py for full solver config.",
    )


def main():
    for path, name in [(FDTD_PKL, "mmi_2x2_fdtd.py"), (HYBRID_PKL, "mmi_2x2_hybrid.py")]:
        if not os.path.exists(path):
            raise SystemExit(
                f"verify_mmi_2x2_hybrid_vs_fdtd: {path} does not exist -- run "
                f"`{name}` first (e.g. `make mmi_2x2_fdtd mmi_2x2_hybrid` from "
                f"examples/, or batch/run_mmi_2x2_examples.batch array indices "
                f"1-2)."
            )

    smat_fdtd = SMatrix.load(FDTD_PKL)
    smat_hybrid = SMatrix.load(HYBRID_PKL)

    wl_fdtd, wl_hybrid = smat_fdtd.wavelengths, smat_hybrid.wavelengths
    # The two grids are ALLOWED to differ, and today they do: mmi_2x2_hybrid.py
    # sweeps N_WL=11 and mmi_2x2_fdtd.py N_WL=51 (raised 2026-08-10). Requiring
    # equality made this script un-runnable against the shipped configuration
    # and would have forced a re-run of the expensive side to compare at all.
    #
    # Resample the FDTD reference onto the hybrid's grid instead. Interpolating
    # the denser side onto the sparser one is the safe direction: no new
    # structure is invented, and every hybrid point is bracketed by real FDTD
    # samples rather than extrapolated. beef.spectral.interp_complex() is the
    # right tool (magnitude + unwrapped phase, NOT real/imag -- see its module
    # docstring for the measured damping that real/imag interpolation causes),
    # with check_unwrap wired on so a grid too coarse for the coefficient's
    # phase rate fails loudly instead of quietly picking the wrong branch.
    resample = not (wl_fdtd.shape == wl_hybrid.shape and np.allclose(wl_fdtd, wl_hybrid))
    if resample:
        if wl_hybrid.min() < wl_fdtd.min() - 1e-12 or wl_hybrid.max() > wl_fdtd.max() + 1e-12:
            raise SystemExit(
                f"verify_mmi_2x2_hybrid_vs_fdtd: the hybrid grid "
                f"[{wl_hybrid.min():.4f}, {wl_hybrid.max():.4f}] is not contained "
                f"in the FDTD grid [{wl_fdtd.min():.4f}, {wl_fdtd.max():.4f}], so "
                f"comparing them would extrapolate the reference. Re-run both "
                f"scripts over the same band before comparing."
            )
        print(
            f"NOTE: wavelength grids differ ({wl_fdtd.size} FDTD vs "
            f"{wl_hybrid.size} hybrid points) -- resampling the FDTD reference "
            f"onto the hybrid grid via beef.spectral.interp_complex().\n"
        )
    wavelengths = wl_hybrid

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
        if resample:
            c_fdtd = interp_complex(wl_fdtd, c_fdtd, wavelengths, check_unwrap=0.5)
        T_fdtd = np.abs(c_fdtd) ** 2
        T_hybrid = np.abs(c_hybrid) ** 2

        dT = T_hybrid - T_fdtd
        max_abs_dT = float(np.max(np.abs(dT)))
        mean_abs_dT = float(np.mean(np.abs(dT)))
        worst_idx = int(np.argmax(np.abs(dT)))
        # Relative to the reference's own scale, not to zero -- a splitter
        # port sits near 0.5, so absolute and relative error are close here,
        # but this stays meaningful if a device/port combination is weak.
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
        ax.set_title(f"{label}: hybrid vs full-device FDTD")
        ax.legend(fontsize=8)
        ax.grid(True)
        ax.set_ylim(0, 1.05)
    plt.suptitle("2×2 MMI — Hybrid vs Full-Device FDTD (mmi_2x2_validation.md gap 5)",
                fontsize=11)
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, "hybrid_vs_fdtd.png"), dpi=300)
    plt.close("all")

    worst_overall = max(r["max_abs_dT"] for r in results.values())
    print(f"Worst max|dT| across both ports: {worst_overall:.4f}")
    print(
        "Informational only -- no threshold enforced here. "
        "mmi_2x2_validation.md gap 4 (engine agreement) is not yet measured "
        "for this device, so there is no basis yet to say how much of this "
        "gap (if any) is calibratable vs. structural. See "
        "resolved/directional_coupler_validation.md items 6/14 for how that "
        "played out on the coupler (uncalibrated hybrid ~5% off; calibrated "
        "0.11% off)."
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
