"""Step 3 smoke test for validation_markdown/resolved/section_mode_beta_multi_monitor_validation.md.

Cheap structural smoke test for FDTDSolver.solve_section_multi_length() --
NOT a correctness check (that's step 4's real numerical-equivalence test,
still to be run at production scale on the cluster). Low res, loose
decay_tol, few modes, small domain: just confirms the merged-domain
multi-monitor code path runs to completion without a Python error and
returns a sanely-shaped, non-degenerate array (each length's monitor reads
a genuinely different position, not an alias of another one).

Recorded result (local M3 MacBook, 2026-08-12): ran in 3.2 s, correct shape
(3, 1, 2, 2), no NaNs, all 3 length slices distinct. See the validation
doc's step 3 for the full record.

Run:
    python /path/to/examples/verification/python/smoke_test_solve_section_multi_length.py
"""
import numpy as np

from beef import MMI, FDTDSolver, SOI_220

mmi = MMI(
    n_inputs=2, n_outputs=2,
    mm_length=28.0, mm_width=2.8,
    waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
    waveguide_offset=1.5, waveguide_height=0.22,
)

solver = FDTDSolver(
    res=10.0, sim_width=3.0, sim_height=1.5, dpml=0.3,
    source_margin=0.0, monitor_margin=0.1, three_d=True,
    material_stack=SOI_220, decay_tol=1e-3,
)

wavelengths = np.array([1.55])
mid = 0.0
lengths = [0.4, 0.8, 1.2]
n_modes = 2

print("Calling solve_section_multi_length()...")
T = solver.solve_section_multi_length(
    mmi, wavelengths, mid, lengths, n_modes,
    monitor_center=0.0, monitor_size=solver.sim_width,
)
print("Returned without error.")
print("T shape:", T.shape, "expected:", (len(lengths), len(wavelengths), n_modes, n_modes))
print("T dtype:", T.dtype)
print("Any NaN:", np.any(np.isnan(T)))
print("Any all-zero length slice:", [np.allclose(T[i], 0) for i in range(len(lengths))])
print(T[:, 0, 0, 0])
print("SMOKE TEST PASS" if T.shape == (len(lengths), len(wavelengths), n_modes, n_modes) and not np.any(np.isnan(T)) else "SMOKE TEST FAIL")
