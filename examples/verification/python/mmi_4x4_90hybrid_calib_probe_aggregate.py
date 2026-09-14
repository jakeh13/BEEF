"""Aggregate a tagged family of single-mode calibration probes into one
neff_calibrated array usable as `neff_override`.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-09-11.

`mmi_4x4_90hybrid_calib_probe.py` already fits beta and n_eff INDEPENDENTLY
per mode -- that is exact, not an approximation: the library's own
`section_mode_beta_details()` post-fit step fits each mode from only that
mode's own phase-vs-length data, so splitting across processes changes
nothing. This script therefore only has to stack the per-mode files into the
(n_wl, n_modes) array that `EMESolver.solve_section_length_sweep(
neff_override=...)` expects, in the same layout as
`calib2_per_mode/aggregated_beta_n16.npz`.

Deliberately does NOT write into calib2_per_mode/: the res=40 calibration
there is what every current result was produced with, and a res=60
calibration must be compared against it and adopted on purpose, not swapped
in by a file landing in the right place.

Run:
    MMI4X4_CALIB_TAG=res60cal python mmi_4x4_90hybrid_calib_probe_aggregate.py
"""
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
PROBE_DIR = os.path.join(_ROOT, "example_output", "mmi_4x4_90hybrid", "calib_probe")

TAG = os.environ.get("MMI4X4_CALIB_TAG", "res60cal")
N_MODES = int(os.environ.get("MMI4X4_CALIB_NMODES", 16))
REF = os.path.join(_ROOT, "example_output", "mmi_4x4_90hybrid",
                   "calib2_per_mode", "aggregated_beta_n16.npz")


def main():
    files = [os.path.join(PROBE_DIR, f"probe_{TAG}_mode{m:02d}.npz")
             for m in range(1, N_MODES + 1)]
    missing = [os.path.basename(f) for f in files if not os.path.exists(f)]
    if missing:
        raise SystemExit(
            f"{len(missing)} of {N_MODES} probe files missing for tag '{TAG}': "
            + ", ".join(missing)
            + "\nAggregate only once every task has landed -- a partial stack "
              "would silently produce a wrong neff_override.")

    d0 = np.load(files[0])
    wavelengths = d0["wavelengths"]
    res = float(d0["res"])
    neff_cal = np.empty((len(wavelengths), N_MODES))
    neff_ms = np.empty((len(wavelengths), N_MODES))
    beta = np.empty((len(wavelengths), N_MODES))
    resid = np.empty((len(wavelengths), N_MODES))

    for m, f in enumerate(files):
        d = np.load(f)
        assert np.allclose(d["wavelengths"], wavelengths), f"{f}: wavelength mismatch"
        assert abs(float(d["res"]) - res) < 1e-9, f"{f}: res mismatch"
        neff_cal[:, m] = d["neff_calibrated"]
        neff_ms[:, m] = d["neff_modesolver"]
        beta[:, m] = d["beta"]
        resid[:, m] = d["resid"]

    offset = neff_cal - neff_ms
    print(f"tag={TAG}  res={res:.0f}  {N_MODES} modes  "
          f"wavelengths={np.round(wavelengths, 4)}")
    print(f"\n{'mode':>5} {'n_eff ms @1.31':>15} {'n_eff cal @1.31':>16}"
          f" {'offset':>9} {'max resid':>10}")
    i = int(np.argmin(np.abs(wavelengths - 1.31)))
    for m in range(N_MODES):
        print(f"{m+1:5d} {neff_ms[i, m]:15.6f} {neff_cal[i, m]:16.6f}"
              f" {offset[i, m]:+9.6f} {resid[:, m].max():10.2e}")
    print(f"\noffset @1.31: mean {offset[i].mean():+.6f}  "
          f"spread {np.ptp(offset[i]):.6f}  (mode-independence check)")
    print(f"beta-fit residual: max {resid.max():.3e}  mean {resid.mean():.3e}")

    if os.path.exists(REF):
        r = np.load(REF)["neff_calibrated"]
        if r.shape == neff_cal.shape:
            d_ = neff_cal - r
            print(f"\nvs the res=40 production calibration "
                  f"(calib2_per_mode/aggregated_beta_n16.npz):")
            print(f"  d n_eff @1.31: mean {d_[i].mean():+.6f}  "
                  f"max|d| {np.abs(d_).max():.6f}")
            print("  NOTE: a nonzero difference here is EXPECTED -- both are")
            print("  under-resolved and converging. This is the size of the")
            print("  calibration's own residual error, not a disagreement.")

    out = os.path.join(PROBE_DIR, f"aggregated_{TAG}_n{N_MODES}.npz")
    np.savez(out, tag=TAG, res=res, wavelengths=wavelengths,
             neff_calibrated=neff_cal, neff_modesolver=neff_ms,
             offset=offset, beta=beta, resid=resid)
    print(f"\nsaved -> {out}")
    print("To use: load ['neff_calibrated'] and pass as neff_override= to "
          "EMESolver.solve_section_length_sweep().")


if __name__ == "__main__":
    main()
