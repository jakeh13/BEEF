"""Merge the 20 per-(offset,branch) outputs of
mmi_1x4_sin_taper_offset_length_sweep.py (5 offsets x 4 branches, one SLURM
array task each) into: the best mm_length per offset, the best (offset,
mm_length) pair overall, and comparison plots against the current production
point (mm_length=206.0, waveguide_offset=2.4).

Run only after all 20 sweep_offset{0-4}_branch{0-3}.npz files exist.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification",
    "mmi_1x4_sin_taper_offset_length_sweep"))

OFFSETS = np.array([2.2, 2.5, 2.8, 3.1, 3.2167])
N_OFFSETS = len(OFFSETS)

files = {(oi, k): os.path.join(output_dir, f"sweep_offset{oi}_branch{k}.npz")
         for oi in range(N_OFFSETS) for k in range(4)}
missing = [f for f in files.values() if not os.path.exists(f)]
if missing:
    raise FileNotFoundError(
        f"{len(missing)} of {len(files)} branch file(s) missing, e.g. {missing[0]} -- "
        f"run mmi_1x4_sin_taper_offset_length_sweep.py for every "
        f"MMI1X4_TOL_OFFSET_IDX=0..{N_OFFSETS - 1} x MMI1X4_OUTPUT_IDX=0..3 pair first."
    )

data = {key: np.load(f) for key, f in files.items()}

mm_lengths = data[(0, 0)]["mm_lengths"]
wavelengths = data[(0, 0)]["wavelengths"]
for key, d in data.items():
    assert np.array_equal(d["mm_lengths"], mm_lengths), f"{key}: mm_lengths grid disagrees"
    assert np.array_equal(d["wavelengths"], wavelengths), f"{key}: wavelengths disagree"
    assert np.isclose(float(d["offset"]), OFFSETS[key[0]]), f"{key}: offset mismatch"

# T_grids[oi, k, L, wl]
T_grids = np.zeros((N_OFFSETS, 4, len(mm_lengths), len(wavelengths)))
output_ys = np.zeros((N_OFFSETS, 4))
for (oi, k), d in data.items():
    T_grids[oi, k] = d["T_k"]
    output_ys[oi, k] = float(d["output_y"])

T_band_mean = T_grids.mean(axis=3)            # [n_off, 4, n_L]
T_total = T_band_mean.sum(axis=1)              # [n_off, n_L]
spread = T_band_mean.max(axis=1) - T_band_mean.min(axis=1)   # [n_off, n_L]

# ---------------------------------------------------------------------------
# Best mm_length per offset (min spread across the 4 outputs -- same metric
# mmi_1x4_sin_body_length_sweep_merge.py uses).
# ---------------------------------------------------------------------------
best_L_idx = np.argmin(spread, axis=1)          # [n_off]
best_L = mm_lengths[best_L_idx]                 # [n_off]
best_spread = spread[np.arange(N_OFFSETS), best_L_idx]
best_total_T = T_total[np.arange(N_OFFSETS), best_L_idx]

print("Best mm_length per offset (min max-min spread across 4 outputs):\n")
print(f"{'offset (um)':>12} {'best mm_length (um)':>20} {'spread':>10} {'total T':>10}")
for oi in range(N_OFFSETS):
    print(f"{OFFSETS[oi]:>12.4f} {best_L[oi]:>20.2f} {best_spread[oi]:>10.4f} "
          f"{best_total_T[oi]:>10.4f}")

# ---------------------------------------------------------------------------
# Overall best (offset, mm_length) pair -- min spread across the whole grid.
# ---------------------------------------------------------------------------
flat_idx = np.argmin(spread)
global_oi, global_Li = np.unravel_index(flat_idx, spread.shape)
print(f"\nGlobal best: offset={OFFSETS[global_oi]:.4f} um, "
      f"mm_length={mm_lengths[global_Li]:.2f} um, "
      f"spread={spread[global_oi, global_Li]:.4f}, "
      f"total T={T_total[global_oi, global_Li]:.4f}")

# Reference: current production point, for comparison only (not part of this
# sweep's own grid -- mm_length=206.0/offset=2.4 is a different offset than
# any point here).
print("\n(Reference -- current production, NOT part of this sweep: "
      "mm_length=206.0, offset=2.4, spread=0.0064, total T=0.6014, "
      "validation doc step 6d/9.)")

np.savez(
    os.path.join(output_dir, "taper_offset_length_sweep.npz"),
    offsets=OFFSETS, mm_lengths=mm_lengths, wavelengths=wavelengths,
    T_grids=T_grids, T_band_mean=T_band_mean, T_total=T_total, spread=spread,
    best_L=best_L, best_spread=best_spread, best_total_T=best_total_T,
    global_offset=OFFSETS[global_oi], global_mm_length=mm_lengths[global_Li],
    output_ys=output_ys,
)
print(f"\nCombined data saved -> {output_dir}/taper_offset_length_sweep.npz")

# ---------------------------------------------------------------------------
# Plot 1: spread vs mm_length, one curve per offset.
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
colors = plt.cm.viridis(np.linspace(0, 0.9, N_OFFSETS))

ax = axes[0]
for oi in range(N_OFFSETS):
    ax.semilogy(mm_lengths, spread[oi], color=colors[oi],
                label=f"offset={OFFSETS[oi]:.4f} um")
    ax.axvline(best_L[oi], color=colors[oi], linestyle=":", lw=0.8)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("max-min spread across 4 outputs")
ax.set_title("Equal-split quality vs body length, per offset")
ax.legend(fontsize=7)
ax.grid(True, which="both")

ax = axes[1]
ax.plot(OFFSETS, best_spread, "o-", color="k", label="spread at best mm_length")
ax2 = ax.twinx()
ax2.plot(OFFSETS, best_total_T, "s--", color="firebrick", label="total T at best mm_length")
ax.set_xlabel("waveguide_offset (um)")
ax.set_ylabel("min spread (across mm_length)")
ax2.set_ylabel("total T at that point", color="firebrick")
ax.set_title("Best achievable split quality vs taper offset")
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc="best")
ax.grid(True)

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- joint taper-offset x mm_length sweep", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "taper_offset_length_sweep.png"), dpi=300)
plt.close("all")
print(f"Plot saved -> {output_dir}/taper_offset_length_sweep.png")
