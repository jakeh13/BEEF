"""End-to-end check on the 11-wavelength 4x4 hybrid: does tier-3 analytic
reconstruction reproduce the real solve at the three shared wavelengths?

linspace(1.29,1.33,11) indices 0/5/10 ARE linspace(1.29,1.33,3) indices 0/1/2
(1.290 / 1.310 / 1.330), so the two runs are directly comparable there.

This tests something the neff-level gate does not. The merge gate compared
n_eff, and passed at 5e-14. But the 3-wl production numbers came from a REAL
`mw.compute_s_matrix()` solve over the length grid, whereas the 11-wl numbers
come from tier 3 reconstructing every length analytically from the merged
neff. So this compares the two paths at the level of the reported device
metrics -- CMRR, phase error, excess loss -- which is where an error would
actually matter.

Expect agreement to be TIGHT but not exact: tier 3 is a different arithmetic
path, not the same one replayed. Anything at 1e-9 or below is round-off.
Degrees-level disagreement means the diagonal reconstruction does not hold for
this cross-section and the body result must not be trusted -- in which case
re-run the example with the merged neff entry deleted, so it takes the real
solve, and compare again.
"""
import os
import sys

import numpy as np

BASE = "/storage/scratch1/8/jhiesener3/duty/example_output/mmi_4x4_90hybrid"
REF = os.path.join(BASE, "mmi_4x4_90hybrid_hybrid_nwl3_PRE_20260915", "hybrid_results.npz")
NEW = os.path.join(BASE, "mmi_4x4_90hybrid_hybrid", "hybrid_results.npz")
IDX = [0, 5, 10]          # 11-wl positions of 1.290 / 1.310 / 1.330 um

for p in (REF, NEW):
    if not os.path.exists(p):
        sys.exit(f"missing {p}")

ref, new = np.load(REF), np.load(NEW)
if len(new["wavelengths"]) != 11:
    sys.exit(f"{NEW} has {len(new['wavelengths'])} wavelengths, expected 11 -- "
             "the assemble job may not have run at the new default")

dw = np.max(np.abs(new["wavelengths"][IDX] - ref["wavelengths"]))
print(f"wavelength alignment: max|dwl| = {dw:.3e} um")
if dw > 1e-12:
    sys.exit("the two grids do not share these wavelengths -- comparison invalid")

print(f"\n{'metric':18s} {'1.290':>12s} {'1.310':>12s} {'1.330':>12s}   max|diff|")
worst, worst_k = 0.0, None
for k in ref.files:
    if k in ("wavelengths", "mm_length", "delta_L"):
        continue
    a, b = np.asarray(ref[k]), np.asarray(new[k])[IDX]
    d = np.abs(a - b)
    print(f"{k:18s} " + " ".join(f"{x:12.3e}" for x in d) + f"   {d.max():.3e}")
    if d.max() > worst:
        worst, worst_k = d.max(), k

print(f"\nscalars: mm_length {float(ref['mm_length']):.4f} -> {float(new['mm_length']):.4f}, "
      f"delta_L {float(ref['delta_L']):.4f} -> {float(new['delta_L']):.4f}")
print(f"worst disagreement: {worst:.3e} on {worst_k}")
if worst < 1e-9:
    print("PASS -- tier-3 reconstruction reproduces the real solve to round-off")
elif worst < 1e-3:
    print("PASS (loose) -- small but above round-off; note it, does not affect conclusions")
else:
    print("FAIL -- investigate before using the 11-wl body result")
    sys.exit(1)
