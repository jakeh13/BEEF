"""4x4 MMI 90-degree hybrid -- solve and cache ONE input-section matrix at a
chosen interface mode count, so a larger-basis cascade can then run as pure
arithmetic.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-09-13 -- knob 3 of
the plan to drive the power-weighted hybrid-vs-FDTD phase error below 1.35 deg.

WHY THIS EXISTS AS A SEPARATE SCRIPT. `solve_section()` results are content-
addressed, and the cache's `lookup_superset()` tier can answer an n-mode
request from a STORED LARGER result by slicing trailing axes -- but only
downward. Going from the production 16 to 24 modes is therefore a genuine
re-solve, not a cache slice, and it is the expensive half of testing a larger
basis. Splitting it out one port per process means the four ports run in
parallel and a failure costs one port, not all four.

WHY 24 AND NOT 20. Measured from the cached 16-mode M_in: the top port's
modal amplitudes are NOT monotone -- modes 8 and 16 sit at ~1e-5 while their
neighbours carry 1.5-4% each, a clear every-8th-mode symmetry suppression. So
the tail is not simply converged at 16; the next modes above the suppressed
one can carry real power. The 16-mode basis captures 99.2-99.3% of input
power, and the missing 0.7-0.8% is large enough (sqrt(0.007) ~ 4.8 deg) to
account for the ~2 deg residual. 24 clears the next suppression node.

Run:
    MMI4X4_PRESOLVE_NMODES=24 MMI4X4_PRESOLVE_LENGTH=15.00 \
        python mmi_4x4_90hybrid_input_presolve.py <port_index_0_to_3>
"""
import os
import sys
import time

import numpy as np

from beef import MMI, FDTDSolver, SOI_220

WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 8.5
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 1.30
TAPER_LENGTH = 4.0
BODY_STUB = 0.25
WL_MIN, WL_MAX, N_WL = 1.29, 1.33, 3

N_MODES = int(os.environ.get("MMI4X4_PRESOLVE_NMODES", 24))
MM_LENGTH = float(os.environ.get("MMI4X4_PRESOLVE_LENGTH", 15.00))
PORT = int(sys.argv[1])
assert 0 <= PORT <= 3, f"port index {PORT} out of range [0,3]"

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
save_dir = os.path.join(_ROOT, "example_sim_cache")

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)


def main():
    print(f"[presolve] port={PORT} n_modes={N_MODES} mm_length={MM_LENGTH}")
    mmi = MMI(
        n_inputs=4, n_outputs=4,
        mm_length=MM_LENGTH, mm_width=MM_WIDTH, waveguide_width=WAVEGUIDE_WIDTH,
        taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
        waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
        use_neighbor_aperture=True,
    )
    half_body = mmi.mm_length / 2
    half_total = mmi.total_length / 2
    x0, x1 = -half_total, -half_body + BODY_STUB

    solver = FDTDSolver(
        res=40.0, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_220, decay_tol=1e-4, decay_method="energy",
        energy_decay_dt=20.0,
    )

    t0 = time.time()
    M = solver.solve_section(
        mmi, wavelengths, x0, x1, n_modes_in=1, n_modes_out=N_MODES,
        stub_left=solver.waveguide_length, stub_right=0.0,
        input_idx=PORT, save_dir=save_dir,
    )
    print(f"[presolve] done in {time.time()-t0:.1f}s, shape {M.shape}")

    # The number this whole knob turns on: how much of the input power the
    # basis actually captures, and whether the tail has converged.
    p = np.abs(M[:, :, 0]) ** 2
    tot = p.sum(axis=1)
    print(f"[presolve] captured power per wavelength: {np.round(tot, 6)}")
    i = 1   # 1.31 um
    print(f"[presolve] modal fractions at 1.31 um (port {PORT}):")
    for k in range(N_MODES):
        print(f"    mode {k+1:2d}: {p[i, k]/tot[i]:8.4%}")
    if N_MODES > 16:
        print(f"[presolve] modes 17-{N_MODES} carry "
              f"{p[i, 16:].sum()/tot[i]:.4%} of captured power")


if __name__ == "__main__":
    main()
