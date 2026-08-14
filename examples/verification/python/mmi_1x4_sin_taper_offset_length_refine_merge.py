"""Merge the 28 per-(offset,branch) outputs of
mmi_1x4_sin_taper_offset_length_refine.py (7 offsets x 4 branches) into the
refined best (offset, mm_length) pick + plots, and compare against the
coarse sweep's own offset=3.1/mm_length=40.0 candidate
(mmi_1x4_sin_taper_offset_length_sweep.py, step 9c).

Run only after all 28 refine_offset{0-6}_branch{0-3}.npz files exist.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification",
    "mmi_1x4_sin_taper_offset_length_refine"))

OFFSETS = np.arange(2.90, 3.2001, 0.05)
N_OFFSETS = len(OFFSETS)

files = {(oi, k): os.path.join(output_dir, f"refine_offset{oi}_branch{k}.npz")
         for oi in range(N_OFFSETS) for k in range(4)}
missing = [f for f in files.values() if not os.path.exists(f)]
if missing:
    raise FileNotFoundError(
        f"{len(missing)} of {len(files)} branch file(s) missing, e.g. {missing[0]} -- "
        f"run mmi_1x4_sin_taper_offset_length_refine.py for every "
        f"MMI1X4_TOL_OFFSET_IDX=0..{N_OFFSETS - 1} x MMI1X4_OUTPUT_IDX=0..3 pair first."
    )

data = {key: np.load(f) for key, f in files.items()}

mm_lengths = data[(0, 0)]["mm_lengths"]
wavelengths = data[(0, 0)]["wavelengths"]
for key, d in data.items():
    assert np.array_equal(d["mm_lengths"], mm_lengths), f"{key}: mm_lengths grid disagrees"
    assert np.array_equal(d["wavelengths"], wavelengths), f"{key}: wavelengths disagree"
    assert np.isclose(float(d["offset"]), OFFSETS[key[0]]), f"{key}: offset mismatch"

T_grids = np.zeros((N_OFFSETS, 4, len(mm_lengths), len(wavelengths)))
output_ys = np.zeros((N_OFFSETS, 4))
for (oi, k), d in data.items():
    T_grids[oi, k] = d["T_k"]
    output_ys[oi, k] = float(d["output_y"])

T_band_mean = T_grids.mean(axis=3)            # [n_off, 4, n_L]
T_total = T_band_mean.sum(axis=1)              # [n_off, n_L]
spread = T_band_mean.max(axis=1) - T_band_mean.min(axis=1)   # [n_off, n_L]

# ---------------------------------------------------------------------------
# Best point overall: not naive argmin(spread) alone (the coarse sweep's own
# merge showed that picks low-throughput points at some offsets) -- rank by
# throughput among well-balanced candidates instead, matching how the coarse
# sweep's real finding was actually identified.
# ---------------------------------------------------------------------------
REL_IMBALANCE_MAX = 0.02   # 2% -- generous vs the coarse candidate's 0.77%
mean_T_per_port = T_total / 4.0
rel_imbalance = np.divide(spread, mean_T_per_port, out=np.full_like(spread, np.inf),
                           where=mean_T_per_port > 1e-9)
well_balanced = rel_imbalance < REL_IMBALANCE_MAX

if np.any(well_balanced):
    masked_T = np.where(well_balanced, T_total, -1.0)
    flat_idx = np.argmax(masked_T)
else:
    flat_idx = np.argmin(spread)
best_oi, best_Li = np.unravel_index(flat_idx, spread.shape)

print(f"Refined best (well-balanced, <{REL_IMBALANCE_MAX*100:.0f}% relative imbalance, "
      f"max throughput): offset={OFFSETS[best_oi]:.4f} um, "
      f"mm_length={mm_lengths[best_Li]:.2f} um, "
      f"total T={T_total[best_oi, best_Li]:.4f}, "
      f"spread={spread[best_oi, best_Li]:.5f} "
      f"({100*rel_imbalance[best_oi, best_Li]:.2f}% relative)")
print(f"  per-output band-mean T: {T_band_mean[best_oi, :, best_Li]}")

print(f"\nComparison -- coarse-grid candidate (offset=3.1, mm_length=40.0): "
      f"total T=0.9598, spread=0.00185 (0.77% relative) [job 11817364]")

# Top 10 candidates overall, for context.
flat = [(OFFSETS[oi], mm_lengths[li], spread[oi, li], T_total[oi, li], rel_imbalance[oi, li])
        for oi in range(N_OFFSETS) for li in range(len(mm_lengths))
        if rel_imbalance[oi, li] < REL_IMBALANCE_MAX]
flat.sort(key=lambda r: -r[3])
print(f"\nTop 10 well-balanced (<{REL_IMBALANCE_MAX*100:.0f}%) candidates by throughput:")
for off, L, s, t, ri in flat[:10]:
    print(f"  offset={off:.4f}  L={L:6.2f}  total_T={t:.4f}  spread={s:.5f}  rel_imb={100*ri:.2f}%")

np.savez(
    os.path.join(output_dir, "taper_offset_length_refine.npz"),
    offsets=OFFSETS, mm_lengths=mm_lengths, wavelengths=wavelengths,
    T_grids=T_grids, T_band_mean=T_band_mean, T_total=T_total, spread=spread,
    rel_imbalance=rel_imbalance,
    best_offset=OFFSETS[best_oi], best_mm_length=mm_lengths[best_Li],
    best_total_T=T_total[best_oi, best_Li], best_spread=spread[best_oi, best_Li],
    output_ys=output_ys,
)
print(f"\nCombined data saved -> {output_dir}/taper_offset_length_refine.npz")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
colors = plt.cm.viridis(np.linspace(0, 0.9, N_OFFSETS))

ax = axes[0]
for oi in range(N_OFFSETS):
    ax.plot(mm_lengths, T_total[oi], color=colors[oi], label=f"offset={OFFSETS[oi]:.2f} um")
ax.axvline(mm_lengths[best_Li], color="k", linestyle=":", lw=1.0)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("Total transmission")
ax.set_title("Refined sweep -- total T vs body length, per offset")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[1]
for oi in range(N_OFFSETS):
    ax.semilogy(mm_lengths, rel_imbalance[oi], color=colors[oi], label=f"offset={OFFSETS[oi]:.2f} um")
ax.axhline(REL_IMBALANCE_MAX, color="gray", linestyle="--", lw=0.8)
ax.axvline(mm_lengths[best_Li], color="k", linestyle=":", lw=1.0)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("relative imbalance (spread / mean per-port T)")
ax.set_title("Refined sweep -- split balance vs body length, per offset")
ax.legend(fontsize=7)
ax.grid(True, which="both")

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- local refinement around offset=3.1/mm_length=40", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "taper_offset_length_refine.png"), dpi=300)
plt.close("all")
print(f"Plot saved -> {output_dir}/taper_offset_length_refine.png")
