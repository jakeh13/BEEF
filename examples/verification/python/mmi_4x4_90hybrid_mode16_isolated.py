"""4x4 MMI 90-degree hybrid -- mode 16 in TRUE isolation (production res=40,
decay_method="energy").

validation_markdown/mmi_4x4_90_hybrid_validation.md follow-up, 2026-08-20.
The full 16-mode production calibration (job 12052385) got modes 1-15
decaying cleanly under decay_method="energy" but stalled on mode 16
(energy ratio flat at 1.0, no decay at all). Checked directly: mode 16 is
NOT near cutoff (margin 0.883 above n_clad, vs mode 1's 1.516 -- the
weakest-confined of the 16 requested, but still solidly guided).

This script tests whether mode 16's stall is an INTRINSIC property of
that mode, or an artifact of being 16th in a long sequential loop
(resource accumulation across 15 prior mp.Simulation builds/resets) --
by patching a runtime copy of FDTDSolver.solve_section_multi_length()
that excites ONLY band 16 (`for k in [16]:` instead of the full
`range(1, n_excite+1)` loop), so mode 16 is the FIRST and ONLY
simulation built in a fresh process. A login-node attempt at this hit
anomalous per-step costs (0.3-1.0 s/step, 10-30x every other run at this
resolution) traced to a contended shared login node (load avg ~10-13,
37 other users) -- moved to batch for a clean, uncontended measurement.

Run: python mmi_4x4_90hybrid_mode16_isolated.py
"""
import inspect
import textwrap
import numpy as np
import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import mmi_4x4_90hybrid_offset_sweep as base
from beef import FDTDSolver, SOI_220
import beef.solvers.fdtd as fdtd_mod

src = inspect.getsource(FDTDSolver.solve_section_multi_length)
src = textwrap.dedent(src)
assert "for k in range(1, n_excite + 1):" in src
src_mod = src.replace("for k in range(1, n_excite + 1):", "for k in [16]:")
src_mod = src_mod.replace("def solve_section_multi_length(", "def solve_section_multi_length_mode16_only(", 1)

ns = dict(fdtd_mod.__dict__)
exec(compile(src_mod, "<mode16_isolated>", "exec"), ns)
import types
solver = FDTDSolver(
    res=40.0, sim_width=base.MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SOI_220, decay_tol=1e-11,
    decay_method="energy", energy_decay_dt=20.0,
)
solver.solve_section_multi_length_mode16_only = types.MethodType(ns["solve_section_multi_length_mode16_only"], solver)

mmi = base.make_mmi(base.BASELINE_OFFSET)
boundaries, body_stub = base.boundaries_for(mmi)
x0, x1 = boundaries[1], boundaries[2]
mid = 0.5 * (x0 + x1)
wavelengths = np.array([1.31])
lengths = [0.885, 1.77, 2.655]

T = solver.solve_section_multi_length_mode16_only(mmi, wavelengths, mid, lengths, 16, measure_diagonal_only=True)
print("done, shape", T.shape, "beta-relevant entry:", T[0, 0, 15, 15])
