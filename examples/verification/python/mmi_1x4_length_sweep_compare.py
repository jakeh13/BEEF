"""EME vs hybrid overlay comparison for the 1x4 MMI's mm_length design sweep
-- combines mmi_1x4_sin_body_length_sweep.py's EME variant (step 5's design
tool) and mmi_1x4_hybrid_body_length_sweep.py's hybrid variant (this
session) on the same axes, same 281-point grid, so the two engines'
predictions can be read directly against each other.

No new solving -- both source .npz files already exist. Note the two
scripts' own output_dir path conventions differ (a real, pre-existing
inconsistency in this codebase's verification scripts, not something this
script fixes): the EME tool's own path arithmetic lands it under
BEEF/example_output/verification/ while every device example and this
session's own new scripts land under duty/example_output/ (see
mmi_1x4_hybrid_body_length_sweep.py's own module docstring for the
save_dir version of the same bug, which WAS load-bearing there and was
fixed; here it only affects where a plot file landed, so left alone).

Run:
    python .../mmi_1x4_length_sweep_compare.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
eme_path = os.path.join(_ROOT, "BEEF", "example_output", "verification",
                         "mmi_1x4_sin_body_length_sweep", "mm_length_sweep.npz")
hybrid_path = os.path.join(_ROOT, "example_output", "verification",
                            "mmi_1x4_hybrid_body_length_sweep", "mm_length_sweep_hybrid.npz")
output_dir = os.path.join(_ROOT, "example_output", "verification", "mmi_1x4_length_sweep_compare")
os.makedirs(output_dir, exist_ok=True)

eme = np.load(eme_path)
hyb = np.load(hybrid_path)

assert np.array_equal(eme["mm_lengths"], hyb["mm_lengths"]), \
    "EME and hybrid sweeps ran on different mm_length grids -- not directly comparable"
mm_lengths = eme["mm_lengths"]

eme_spread = eme["spread"]
hyb_spread = hyb["spread"]
eme_total = eme["T_total_grid"].mean(axis=1)
hyb_total = hyb["T_total_grid"].mean(axis=1)
eme_best = float(eme["mm_length_best"])
hyb_best = float(hyb["mm_length_best"])
MM_LENGTH_CHOSEN = 206.0

print(f"EME    sweep optimum: {eme_best:.1f} um (spread={eme_spread[np.argmin(eme_spread)]:.4f})")
print(f"Hybrid sweep optimum: {hyb_best:.1f} um (spread={hyb_spread[np.argmin(hyb_spread)]:.4f})")
print(f"Optima {abs(eme_best - hyb_best):.1f} um apart on a {mm_lengths[-1]-mm_lengths[0]:.0f} um grid span")

chosen_idx = int(np.argmin(np.abs(mm_lengths - MM_LENGTH_CHOSEN)))
print(f"\nAt chosen production mm_length={MM_LENGTH_CHOSEN}:")
print(f"  EME    spread={eme_spread[chosen_idx]:.4f}  total T={eme_total[chosen_idx]:.4f}")
print(f"  Hybrid spread={hyb_spread[chosen_idx]:.4f}  total T={hyb_total[chosen_idx]:.4f}")
eme_T = eme["T_grids"].mean(axis=2)   # [4, n_L]
hyb_T = hyb["T_grids"].mean(axis=2)
print(f"  EME    per-output: {[f'{t:.4f}' for t in eme_T[:, chosen_idx]]}  "
      f"(outer {'>' if eme_T[0,chosen_idx]>eme_T[1,chosen_idx] else '<'} inner)")
print(f"  Hybrid per-output: {[f'{t:.4f}' for t in hyb_T[:, chosen_idx]]}  "
      f"(outer {'>' if hyb_T[0,chosen_idx]>hyb_T[1,chosen_idx] else '<'} inner)")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

ax = axes[0]
ax.plot(mm_lengths, eme_total, color="steelblue", lw=1.3, label="EME total T")
ax.plot(mm_lengths, hyb_total, color="firebrick", lw=1.3, label="Hybrid total T")
ax.axvline(MM_LENGTH_CHOSEN, color="purple", linestyle="-.", lw=1.2,
           label=f"chosen (production)={MM_LENGTH_CHOSEN:.0f} um")
ax.axhline(1.0, color="gray", linestyle=":", lw=0.6, label="lossless (T=1.0)")
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("Band-mean total transmission")
ax.set_title("Total throughput: EME vs Hybrid")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[1]
ax.semilogy(mm_lengths, eme_spread, color="steelblue", lw=1.3, label=f"EME (opt={eme_best:.0f} um)")
ax.semilogy(mm_lengths, hyb_spread, color="firebrick", lw=1.3, label=f"Hybrid (opt={hyb_best:.0f} um)")
ax.axvline(MM_LENGTH_CHOSEN, color="purple", linestyle="-.", lw=1.2,
           label=f"chosen={MM_LENGTH_CHOSEN:.0f} um")
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("max-min spread across 4 outputs")
ax.set_title("Equal-split quality: EME vs Hybrid")
ax.legend(fontsize=7)
ax.grid(True, which="both")

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- EME vs Hybrid body-length sweep", fontsize=11)
plt.tight_layout()
out_path = os.path.join(output_dir, "eme_vs_hybrid_length_sweep.png")
fig.savefig(out_path, dpi=300)
plt.close("all")
print(f"\nPlot saved -> {out_path}")
