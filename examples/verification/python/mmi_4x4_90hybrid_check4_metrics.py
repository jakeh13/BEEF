"""Check 4, device-level: do the hybrid and full-device FDTD agree on the
metrics this device is actually judged by?

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-25.

Step 6b completed for the first time (jobs 12272571 / 12272573, both driven
ports, `decay_tol=1e-5`), and the raw |T|^2 comparison against the hybrid
(job 11778128) showed the two methods disagreeing by up to 0.0255 -- 48x the
1e-5 convergence error -- with the gap growing toward SHORT wavelength
(+30%/+50% at 1.26 um, ~+4% at 1.31, ~0% at 1.36).

Raw transmission is not what the device is judged on, though. CMRR, imbalance
and the cross-excitation phase error are RATIOS and DIFFERENCES between ports,
so a uniform loss largely cancels out of them. The disagreement above is NOT
uniform (it is wavelength- and port-dependent), so whether it survives into the
metrics has to be computed rather than argued either way.

This recomputes the paper's metrics from BOTH S-matrices using the SAME
`SMatrix` methods and the SAME construction the hybrid script uses, and diffs
them. No new simulation: the FDTD S-matrices come back from the cache
(`cache_only=True`) and the hybrid ones from their published pickles.

PHASE CONVENTION: `cross_excitation_phase_diff_deg` is called as
`smat_port3.cross_excitation_phase_diff_deg(smat_port1, port=...)` -- the
self/other order matters and this is the corrected direction from finding A11
(meep/meow's S-parameter phase is conjugated relative to the paper's
convention; the un-swapped call reproduced -target instead of +target across
all four ports). Do not reorder these arguments.

Run:
    python mmi_4x4_90hybrid_check4_metrics.py
"""
import os
import sys

import numpy as np

from beef import MMI, FDTDSolver, SMatrix, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                    "example_output", "mmi_4x4_90hybrid"))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                        "example_sim_cache"))

MM_WIDTH = 6.0
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
MM_LENGTH = 0.75 * 108.9195 - 2.75
PORT_NAMES = ["out5", "out6", "out7", "out8"]
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}
wavelengths = np.linspace(1.26, 1.36, 3)

mmi = MMI(n_inputs=4, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
          waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
          waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
          use_neighbor_aperture=True)

solver = FDTDSolver(
    res=40.0, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    num_input_modes=1, num_output_modes=1,
    output_monitor_centers=list(mmi.output_ys),
    output_monitor_widths=[WAVEGUIDE_OFFSET] * 4,
    material_stack=SOI_220, decay_tol=1e-5,
    decay_method="energy", energy_decay_dt=20.0,
)


def fdtd_smatrix(**kw):
    """Full-device FDTD result, relabelled onto the paper's port names.

    solve() stores out_0..out_3; the hybrid stores out5..out8. Same physical
    ordering, so relabel rather than reorder.
    """
    raw = solver.solve(mmi, wavelengths, input_mode=1, save_dir=save_dir,
                       cache_only=True, **kw)
    out = SMatrix(wavelengths=wavelengths)
    for j, name in enumerate(PORT_NAMES):
        out.set(np.asarray(raw.coeff(port=f"out_{j}")), port=name)
    return out


def metrics(s1, s3):
    """The paper's metrics, computed exactly as mmi_4x4_90hybrid_hybrid.py does."""
    m = {}
    for lab, s in (("port1", s1), ("port3", s3)):
        m[f"imbalance_{lab}"] = np.asarray(s.imbalance_db(PORT_NAMES))
        m[f"excess_loss_{lab}"] = np.asarray(s.excess_loss_db(PORT_NAMES))
        m[f"cmrr_I_{lab}"] = np.asarray(s.cmrr_db("out5", "out8"))
        m[f"cmrr_Q_{lab}"] = np.asarray(s.cmrr_db("out6", "out7"))
    # A11's corrected direction -- port3 as self, port1 as other.
    phi = {j: np.asarray(s3.cross_excitation_phase_diff_deg(s1, port=PORT_NAMES[j]))
           for j in range(4)}
    for j in (1, 2, 3):
        raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
        m[f"dphi_{j+5}"] = ((raw + 180.0) % 360.0) - 180.0
    return m


def main():
    print("CHECK4: loading full-device FDTD from cache (no new simulation)...")
    f1 = fdtd_smatrix(port="bot")
    f3 = fdtd_smatrix(input_idx=2)
    print("CHECK4: loading hybrid from published pickles...")
    h1 = SMatrix.load(os.path.join(_OUT, "mmi_4x4_90hybrid_hybrid",
                                   "smatrix_hybrid_port1.pkl"))
    h3 = SMatrix.load(os.path.join(_OUT, "mmi_4x4_90hybrid_hybrid",
                                   "smatrix_hybrid_port3.pkl"))

    F, H = metrics(f1, f3), metrics(h1, h3)

    print(f"\n{'metric':22s} {'wl':>6s}  {'hybrid':>10s} {'FDTD':>10s} {'diff':>10s}")
    print("-" * 64)
    for k in sorted(F):
        for i, w in enumerate(wavelengths):
            a, b = float(H[k][i]), float(F[k][i])
            print(f"{k:22s} {w:6.2f}  {a:10.3f} {b:10.3f} {b-a:+10.3f}")
        print()

    print("=== headline numbers (the ones the device is judged on) ===")
    for lab, M in (("hybrid", H), ("full-device FDTD", F)):
        mx = max(np.abs(M[f"dphi_{j}"]).max() for j in (6, 7, 8))
        mn = min(M[k].min() for k in M if k.startswith("cmrr"))
        xl = max(M[k].max() for k in M if k.startswith("excess_loss"))
        ib = max(M[k].max() for k in M if k.startswith("imbalance"))
        print(f"  {lab:18s} max|dphi|={mx:7.3f} deg   min CMRR={mn:7.3f} dBe   "
              f"max XL={xl:6.3f} dB   max imb={ib:6.3f} dB")
    print("\n  spec: |dphi| <= 5 deg, CMRR >= 20 dBe")


if __name__ == "__main__":
    main()
