"""TEMPORARY -- field-overlap follow-up to audit_interface_modes.py's mmi_1x4
result (validation_markdown/mmi_1x4_sin300_validation.md step 8 check 2).

Delete with the rest of examples/verification/ once signed off.

WHY THIS SCRIPT EXISTS
-----------------------
audit_interface_modes.py's nearest-neighbour-in-n_eff correspondence check
(job 11697897, res=40, N_MODES=12, matching mmi_1x4_hybrid.py exactly) came
back ORDER PRESERVED: False -- both meep bands 0 and 1 nearest-matched to the
SAME meow band (1), because the top two guided modes are only 0.0036 apart in
n_eff while the systematic meep<->meow bias at this cross-section is
0.0022-0.0024 -- bigger than half the mode spacing. That is exactly the
failure mode this codebase's own mode_identity.py module docstring names
("several higher-order modes can crowd together in n_eff closely enough that
nearest-neighbour picks the wrong pairing outright") and the validation doc's
own step 8 note predicted ("a wider nitride body supports more near-degenerate
modes than either existing MMI ... more likely to fail here, not less").

A 1-D scalar mismatch does not by itself say whether the two engines'
top-2 bases are a genuine mismatch or just crowded past what nearest-neighbour
n_eff can resolve. This runs the shape-based check instead: field_overlap_matrix
+ the Hungarian assignment (audit_interface_modes.overlap_correspondence(),
itself built on beef.mode_identity), which compares mode SHAPE rather than a
single number and is immune to this exact crowding failure (module docstring).

Deliberately reuses overlap_correspondence() as-is, at its own cmp_res=20.0
default -- coarser than the res=40 production/audit resolution, a deliberate
cost tradeoff its own docstring explains (halving resolution cuts the
per-mode field-sampling cost ~4x; still resolves this body's mode lobes at a
handful of points per lobe, plenty for an aggregate overlap integral). This
means the comparison happens at a DIFFERENT resolution than job 11697897's
n_eff numbers -- expected and fine, since the shape-based check does not
depend on the two engines numerically agreeing at a shared resolution the way
nearest-neighbour n_eff matching does (see overlap_correspondence()'s own
docstring for why that is the whole point).

Run:
    python .../verify_mmi1x4_mode_overlap.py
"""
import json
import os

import numpy as np

import audit_interface_modes as aim

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "interface_mode_audit")
os.makedirs(output_dir, exist_ok=True)

cfg = aim._mmi_1x4()
n_modes = cfg["n_modes"]
n_clad = aim.n_cladding(cfg["stack"], cfg["wl"])

print("=" * 72)
print(f"[overlap] mmi_1x4 -- {cfg['note']}")
print(f"[overlap] N_MODES={n_modes}, wl={cfg['wl']} um, n_cladding={n_clad:.4f}, "
      f"section={cfg['section']}")
print("=" * 72)

result = aim.overlap_correspondence(cfg, n_modes, n_clad)

n_a, n_b = result["n_guided_meep"], result["n_guided_meow"]
corr, margin = result["correspondence"], result["margin"]
neff_a, neff_b = np.array(result["neff_meep_guided"]), np.array(result["neff_meow_guided"])

print(f"\n[overlap] guided: meep {n_a}/{n_modes}, meow {n_b}/{n_modes} "
      f"(cmp_res={result['cmp_res']})")
print(f"[overlap] band     n_eff(meep)   -> meow band  n_eff(meow)   margin")
for i in range(n_a):
    j = corr[i]
    if j >= 0:
        print(f"[overlap]  {i:3d}   {neff_a[i]:11.6f}  ->  {j:9d}   "
              f"{neff_b[j]:11.6f}   {margin[i]:.4f}")
    else:
        print(f"[overlap]  {i:3d}   {neff_a[i]:11.6f}  ->      NONE")

print(f"\n[overlap] ORDER PRESERVED (shape-based): {result['order_preserved']}")
print(f"[overlap] worst overlap on the matched diagonal: {result['worst_overlap_on_diagonal']:.4f}")
print(f"[overlap] (for reference: audit_interface_modes.py job 11697897's own "
      f"scalar n_eff check reported ORDER PRESERVED: False at res=40, matched "
      f"bands 0 and 1 to the same meow band)")

out_path = os.path.join(output_dir, "mmi_1x4_overlap_correspondence.json")
with open(out_path, "w") as fh:
    json.dump(result, fh, indent=2)
print(f"\n[overlap] written -> {out_path}")
