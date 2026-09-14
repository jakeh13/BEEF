"""beta(width) calibration curve for the waisted 2x2's tapered body.

WHY THIS EXISTS. `calibrate_body_neff` propagates an EME body on the beta its
neighbouring FDTD section actually accumulates, and it is what closed the
engine-agreement gap on every UNIFORM-body device in this repo. It is
hard-refused for a tapered body, for two separate reasons:

  1. mechanically, it replaces the section with a DIAGONAL
     exp(2j*pi*neff*L/wl) and skips the mode solve -- no mode mixing, which is
     exactly the physics a taper exists to produce; and
  2. practically, an intermediate taper width has NO neighbouring FDTD section
     to take beta from. `section_mode_beta()` measures beta by fitting
     dphi/dL, and a 0.75 um ramp slice is far too short to fit a slope on.

`EMESolver.solve_section(cell_neff_override=...)` fixes (1) -- it keeps the
multi-cell solve and every interface S-matrix and substitutes n_eff only in
the propagation inside each cell. THIS script fixes (2): measure beta on
UNIFORM strips across the width range the ramp sweeps, then interpolate each
ramp cell's beta from that curve.

The probe is `UniformCrossSection([(0.0, W)], length=L)` -- a straight strip of
the width in question, simulated by the SAME FDTD solver settings the hybrid's
port sections use, so the calibrated beta carries the same discretisation the
cascade is spliced to. That is the entire point of calibration: the chain has
to agree with itself.

Cost: one merged-domain FDTD run per width (use_merged_domain=True fits all
fit-lengths in a single simulation rather than one run per length).

Run:
    MMI2X2_BW_RES=20 MMI2X2_BW_NW=3 python .../measure_beta_width_curve.py
    python .../measure_beta_width_curve.py          # real config
"""
import json
import os
import time

import numpy as np


def _is_master():
    """True on rank 0 only.

    Under `srun --mpi=pmi2 -n 8` EVERY rank runs this script, so an unguarded
    write has all eight racing on one path. Measured on the first real run:
    4 of 9 shards came back with a stray trailing brace (valid leading
    document, one extra byte) because two ranks' writes interleaved. Same
    defect the hybrid shard writer already fixed; this script never got the
    guard.
    """
    try:
        import meep as _mp
        return bool(_mp.am_master())
    except Exception:
        return True

from beef import FDTDSolver, SOI_220
from beef.geometry.uniform import UniformCrossSection


def _env_f(n, d):
    return float(os.environ.get(n, d))


def _env_i(n, d):
    return int(os.environ.get(n, d))


_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
out_dir = os.path.join(_ROOT, "example_output", "mmi_2x2_tapered")
os.makedirs(out_dir, exist_ok=True)
save_dir = os.path.join(_ROOT, "example_sim_cache")

RES = _env_f("MMI2X2_BW_RES", 40.0)
N_MODES = _env_i("MMI2X2_BW_NMODES", 10)
N_WIDTHS = _env_i("MMI2X2_BW_NW", 9)
WL_MIN, WL_MAX = 1.530, 1.565
N_WL = _env_i("MMI2X2_BW_NWL", 11)
DECAY_TOL = _env_f("MMI2X2_BW_DECAYTOL", 1e-11)

# The ramp sweeps the body from the face width down to the waist. Both
# endpoints are included so no cell needs extrapolation.
W_FACE, W_WAIST = 2.80, 2.24
widths = np.linspace(W_WAIST, W_FACE, N_WIDTHS)
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

# Matched to mmi_2x2_waist_hybrid_vs_fdtd.py's _fdtd_common, deliberately: a
# beta measured at different res/decay_tol would calibrate the body onto a
# discretisation the cascade's own sections do not have.
solver = FDTDSolver(res=RES, sim_width=5.0, sim_height=1.5, dpml=0.5,
                    source_margin=0.0, monitor_margin=0.25, three_d=True,
                    material_stack=SOI_220, decay_tol=DECAY_TOL)

path = os.path.join(
    out_dir, f"beta_width_curve_res{RES:g}_nm{N_MODES}_nwl{N_WL}_nw{N_WIDTHS}.json")

print("=" * 84)
print("beta(width) curve for the waisted 2x2 ramp")
print("=" * 84)
print(f"  widths {W_WAIST} .. {W_FACE} um in {N_WIDTHS} steps: "
      + ", ".join(f"{w:.3f}" for w in widths))
print(f"  res={RES:g}  {N_MODES} modes  {N_WL} wavelengths  decay_tol={DECAY_TOL:g}")
print(f"  -> {os.path.basename(path)}")
print("=" * 84, flush=True)

# One width per job when sharded (MMI2X2_BW_WIDX), because the widths are
# independent and a 9-width serial run does not fit the -ece 5 h cap. Merge
# with MMI2X2_BW_MERGE=1, same pattern the hybrid wavelength shards use.
SHARD = os.environ.get("MMI2X2_BW_WIDX", "")
MERGE = os.environ.get("MMI2X2_BW_MERGE", "") == "1"

if MERGE:
    import glob
    shards = sorted(glob.glob(path.replace(".json", "_w*.json")))
    rows = []
    for f in shards:
        # raw_decode takes the FIRST complete document and ignores trailing
        # bytes, so a shard written by the pre-guard multi-rank race still
        # loads. Verified by hand: the leading document is complete and the
        # junk is a single stray brace.
        txt = open(f).read().lstrip()
        doc, end = json.JSONDecoder().raw_decode(txt)
        if txt[end:].strip():
            print(f"  note: {os.path.basename(f)} had "
                  f"{len(txt[end:].strip())} trailing byte(s) from the "
                  "multi-rank write race; using the leading document")
        rows.extend(doc["rows"])
    rows.sort(key=lambda r: r["width"])
    got = [round(r["width"], 6) for r in rows]
    want = [round(float(w), 6) for w in widths]
    if got != want:
        raise SystemExit(
            f"MERGE: need widths {want}, have {got}. A partial curve would "
            "extrapolate silently over the missing range -- run the missing "
            "widths first.")
    with open(path, "w") as fh:
        json.dump({"config": {"res": RES, "n_modes": N_MODES,
                              "decay_tol": DECAY_TOL,
                              "waveguide_height": 0.22,
                              "sim_width": 5.0, "sim_height": 1.5},
                   "wavelengths": wavelengths.tolist(),
                   "widths_requested": widths.tolist(),
                   "complete": True, "rows": rows}, fh, indent=2)
    print(f"MERGE: {len(rows)} widths -> {path}")
    raise SystemExit(0)

if SHARD != "":
    _k = int(SHARD)
    if not 0 <= _k < N_WIDTHS:
        raise SystemExit(f"MMI2X2_BW_WIDX={_k} out of range for {N_WIDTHS}")
    widths = widths[_k:_k + 1]
    path = path.replace(".json", f"_w{_k}.json")
    print(f"[shard] width index {_k}: {widths[0]:.3f} um -> "
          f"{os.path.basename(path)}", flush=True)

rows = []
for i, w in enumerate(widths):
    geom = UniformCrossSection([(0.0, float(w))], length=4.0,
                               waveguide_height=0.22)
    t0 = time.time()
    beta = solver.section_mode_beta(
        geom, wavelengths, -2.0, 2.0, N_MODES,
        save_dir=save_dir,
    )
    beta = np.asarray(beta)                      # [n_wl, n_modes], rad/um
    neff = beta * wavelengths[:, None] / (2 * np.pi)
    dt = time.time() - t0
    rows.append({"width": float(w), "beta": beta.tolist(),
                 "neff": neff.real.tolist(), "seconds": dt})
    print(f"  [{i+1}/{N_WIDTHS}] w={w:.3f} um  {dt:6.0f}s  "
          f"n_eff(band centre): "
          + ", ".join(f"{v:.4f}" for v in neff[len(wavelengths) // 2].real[:4])
          + " ...", flush=True)
    # Written after every width so a timeout keeps finished work -- the same
    # streaming discipline the free sweep uses.
    if _is_master():
        with open(path, "w") as fh:
            json.dump({"config": {"res": RES, "n_modes": N_MODES,
                                  "decay_tol": DECAY_TOL,
                                  "waveguide_height": 0.22,
                                  "sim_width": 5.0, "sim_height": 1.5},
                       "wavelengths": wavelengths.tolist(),
                       "widths_requested": widths.tolist(),
                       "complete": (SHARD == "") and len(rows) == N_WIDTHS,
                       "rows": rows}, fh, indent=2)

print(f"\nwritten -> {path}")
print("Feed this to HybridSolver(cell_neff_fn=...) via "
      "beta_width_curve.load_curve() to calibrate a tapered body cell by cell.")
