"""Merge the four per-branch outputs of mmi_1x4_sin_body_length_sweep.py
(one SLURM job per output_idx, run in parallel -- see that script's step 5
comment) into the combined mm_length pick + evidence plot.

Run only after all four mm_length_sweep_branch{0,1,2,3}.npz files exist.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification", "mmi_1x4_sin_body_length_sweep"))

branch_files = [os.path.join(output_dir, f"mm_length_sweep_branch{k}.npz") for k in range(4)]
missing = [f for f in branch_files if not os.path.exists(f)]
if missing:
    raise FileNotFoundError(
        f"missing branch file(s): {missing} -- run mmi_1x4_sin_body_length_sweep.py "
        f"with MMI1X4_OUTPUT_IDX=0..3 first (one job per branch)."
    )

branches = [np.load(f) for f in branch_files]

mm_lengths = branches[0]["mm_lengths"]
wavelengths = branches[0]["wavelengths"]
for b in branches[1:]:
    assert np.array_equal(b["mm_lengths"], mm_lengths), "branch mm_lengths grids disagree"
    assert np.array_equal(b["wavelengths"], wavelengths), "branch wavelengths disagree"

output_ys = [float(b["output_y"]) for b in branches]
waveguide_width = float(branches[0]["waveguide_width"])
taper_width = float(branches[0]["taper_width"])
taper_length = float(branches[0]["taper_length"])
waveguide_offset = float(branches[0]["waveguide_offset"])
mm_width = float(branches[0]["mm_width"])
n_body = int(branches[0]["n_body"])
n_taper = int(branches[0]["n_taper"])
c_taper = int(branches[0]["c_taper"])

T_grids = np.stack([b["T_k"] for b in branches], axis=0)   # [4, n_L, n_wl]
T_total_grid = T_grids.sum(axis=0)                          # [n_L, n_wl]

# ---------------------------------------------------------------------------
# Pick the best equal-split length: minimize the spread across the 4
# outputs' band-mean transmission (the 1x4 analogue of the coupler's 50:50
# crossing / the 2x2's argmin-to-nominal, generalized to N ports).
# ---------------------------------------------------------------------------
T_band_mean = T_grids.mean(axis=2)   # [4, n_L] -- band-mean per output per length
spread = T_band_mean.max(axis=0) - T_band_mean.min(axis=0)   # [n_L]
best_idx = int(np.argmin(spread))
mm_length_best = float(mm_lengths[best_idx])

print(f"Best equal-split length (min max-min spread across 4 outputs): "
      f"mm_length = {mm_length_best:.3f} um")
print(f"  spread at this length: {spread[best_idx]:.4f}")
print(f"  per-output band-mean T: {[f'{t:.4f}' for t in T_band_mean[:, best_idx]]}")
print(f"  total T: {T_total_grid[best_idx].mean():.4f}")

_is_local_min = (
    (spread[1:-1] < spread[:-2]) & (spread[1:-1] < spread[2:])
)
_local_min_idxs = np.flatnonzero(_is_local_min) + 1
print(f"\nAll local minima in the swept range (length, spread):")
for i in _local_min_idxs:
    print(f"  {mm_lengths[i]:7.2f} um   spread={spread[i]:.4f}   "
          f"total T={T_total_grid[i].mean():.4f}")

np.savez(
    os.path.join(output_dir, "mm_length_sweep.npz"),
    mm_lengths=mm_lengths, wavelengths=wavelengths,
    T_grids=T_grids, T_total_grid=T_total_grid, spread=spread,
    mm_length_best=mm_length_best,
    waveguide_width=waveguide_width, taper_width=taper_width,
    taper_length=taper_length, waveguide_offset=waveguide_offset,
    mm_width=mm_width, n_body=n_body, n_taper=n_taper, c_taper=c_taper,
)
print(f"\nCombined data saved -> {output_dir}/mm_length_sweep.npz")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

ax = axes[0]
colors = ["steelblue", "darkorange", "seagreen", "firebrick"]
for k in range(4):
    ax.plot(mm_lengths, T_band_mean[k], label=f"output {k} (y={output_ys[k]:+.2f} um)",
           color=colors[k])
ax.plot(mm_lengths, T_total_grid.mean(axis=1), "k--", lw=1, label="total")
ax.axvline(mm_length_best, color="gray", linestyle=":", lw=1.2,
          label=f"chosen mm_length={mm_length_best:.1f} um")
ax.axhline(0.25, color="gray", linestyle=":", lw=0.6)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("Band-mean transmission")
ax.set_title("1x4 MMI -- output power vs body length")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[1]
ax.semilogy(mm_lengths, spread, color="k")
ax.axvline(mm_length_best, color="gray", linestyle=":", lw=1.2)
ax.set_xlabel("mm_length (um)")
ax.set_ylabel("max-min spread across 4 outputs")
ax.set_title("Equal-split quality vs body length")
ax.grid(True, which="both")

plt.suptitle("1x4 MMI on 300 nm Si3N4 -- body length design sweep", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "mm_length_sweep.png"), dpi=300)
plt.close("all")
print(f"Plot saved -> {output_dir}/mm_length_sweep.png")
