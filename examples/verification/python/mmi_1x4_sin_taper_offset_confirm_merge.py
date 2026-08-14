"""Merge the 4 per-branch outputs of mmi_1x4_sin_taper_offset_confirm.py
(offset=3.1 fixed, N_WL=11 full-band confirmation) into the final result:
best mm_length, per-output/total transmission, and a direct comparison
against production (mm_length=206.0, offset=2.4) and the two design-tool
passes that found this candidate (coarse: job 11817364; refine: job
11839158).

Run only after all 4 confirm_branch{0-3}.npz files exist.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification",
    "mmi_1x4_sin_taper_offset_confirm"))

files = [os.path.join(output_dir, f"confirm_branch{k}.npz") for k in range(4)]
missing = [f for f in files if not os.path.exists(f)]
if missing:
    raise FileNotFoundError(
        f"missing branch file(s): {missing} -- run mmi_1x4_sin_taper_offset_confirm.py "
        f"with MMI1X4_OUTPUT_IDX=0..3 first."
    )

branches = [np.load(f) for f in files]
mm_lengths = branches[0]["mm_lengths"]
wavelengths = branches[0]["wavelengths"]
for b in branches[1:]:
    assert np.array_equal(b["mm_lengths"], mm_lengths), "branch mm_lengths grids disagree"
    assert np.array_equal(b["wavelengths"], wavelengths), "branch wavelengths disagree"

offset = float(branches[0]["offset"])
output_ys = [float(b["output_y"]) for b in branches]

T_grids = np.stack([b["T_k"] for b in branches], axis=0)   # [4, n_L, n_wl]
T_band_mean = T_grids.mean(axis=2)                          # [4, n_L]
T_total_grid = T_band_mean.sum(axis=0)                       # [n_L]
spread_grid = T_band_mean.max(axis=0) - T_band_mean.min(axis=0)   # [n_L]
rel_imbalance_grid = spread_grid / np.clip(T_total_grid / 4.0, 1e-12, None)

best_idx = int(np.argmax(np.where(rel_imbalance_grid < 0.02, T_total_grid, -1.0)))
if T_total_grid[best_idx] < 0:
    best_idx = int(np.argmin(spread_grid))

L_best = float(mm_lengths[best_idx])
T_best = float(T_total_grid[best_idx])
spread_best = float(spread_grid[best_idx])
ri_best = float(rel_imbalance_grid[best_idx])

print(f"CONFIRMED (N_WL=11, full C-band): offset={offset:.4f} um, "
      f"mm_length={L_best:.2f} um -> total T={T_best:.4f} "
      f"(IL={10*np.log10(max(T_best,1e-20)):.2f} dB), "
      f"spread={spread_best:.5f} ({100*ri_best:.2f}% relative)")
print(f"  per-output band-mean T: {T_band_mean[:, best_idx]}")
print(f"  output y positions: {output_ys}")

print("\nComparison across every candidate measured on this device:")
print(f"  {'design point':<45}{'mm_length':>12}{'total T':>10}{'rel_imb':>10}")
print(f"  {'current production':<45}{'206.0':>12}{'0.6014':>10}{'4.26%':>10}")
print(f"  {'44um + offset~3.0-3.2 (coarse, N_WL=1)':<45}{'44.0':>12}{'~0.755':>10}{'~15%':>10}")
print(f"  {'this sweep, coarse (N_WL=5)':<45}{'40.0':>12}{'0.9598':>10}{'0.77%':>10}")
print(f"  {'this sweep, refined (N_WL=1)':<45}{'39.9-40.1':>12}{'0.9707':>10}{'0.33%':>10}")
print(f"  {'this sweep, CONFIRMED (N_WL=11)':<45}{L_best:>12.2f}{T_best:>10.4f}{100*ri_best:>9.2f}%")

np.savez(
    os.path.join(output_dir, "taper_offset_confirm.npz"),
    mm_lengths=mm_lengths, wavelengths=wavelengths, offset=offset,
    T_grids=T_grids, T_band_mean=T_band_mean, T_total_grid=T_total_grid,
    spread_grid=spread_grid, rel_imbalance_grid=rel_imbalance_grid,
    L_best=L_best, T_best=T_best, spread_best=spread_best, ri_best=ri_best,
    output_ys=output_ys,
)
print(f"\nCombined data saved -> {output_dir}/taper_offset_confirm.npz")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
colors = ["steelblue", "darkorange", "seagreen", "firebrick"]

ax = axes[0]
for k in range(4):
    ax.plot(mm_lengths, T_band_mean[k], label=f"output {k} (y={output_ys[k]:+.2f} um)", color=colors[k])
ax.plot(mm_lengths, T_total_grid, "k--", lw=1.5, label="total")
ax.axvline(L_best, color="gray", linestyle=":", lw=1.2, label=f"chosen L={L_best:.1f} um")
ax.axhline(0.25, color="gray", linestyle=":", lw=0.6)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("Band-mean transmission (N_WL=11, full C-band)")
ax.set_title(f"Confirmation -- output power vs body length (offset={offset:.2f} um)")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[1]
ax.semilogy(mm_lengths, rel_imbalance_grid, color="k")
ax.axhline(0.02, color="gray", linestyle="--", lw=0.8, label="2% reference")
ax.axvline(L_best, color="gray", linestyle=":", lw=1.2)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("relative imbalance")
ax.set_title("Confirmation -- split balance vs body length")
ax.legend(fontsize=7)
ax.grid(True, which="both")

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- full-fidelity (N_WL=11) confirmation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "taper_offset_confirm.png"), dpi=300)
plt.close("all")
print(f"Plot saved -> {output_dir}/taper_offset_confirm.png")
