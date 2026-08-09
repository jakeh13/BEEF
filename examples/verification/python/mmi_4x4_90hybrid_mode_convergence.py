"""4x4 MMI 90-degree hybrid -- taper/body mode-count convergence sweep.

Not one of validation_markdown/mmi_4x4_90_hybrid_validation.md's original
Steps 0-9 -- a follow-up to the taper-mode-count investigation
(N_MODES_TAPER 4->12, 2026-08-08): visually inspecting this device's own
field plots showed N_MODES_TAPER=4 left the launched mode unphysically
confined instead of fanning into a clean self-imaging pattern, and 12
visibly fixed it. The user asked for an actual convergence test on BOTH
the taper's and the MMI body's mode counts, not just a single before/after
bump -- this script measures that directly rather than picking another
single value by eye.

Method: pure EME (mmi_4x4_90hybrid_eme.py's own method/geometry), but ONE
port pair only (input_idx=0, output_idx=0 -- port 1 -> out5) instead of all
8 -- a convergence check only needs to see whether a single coefficient
stops changing as the mode count grows, not full-device coverage. This is
the same economical single-port-pair proxy the taper_length scoping sweep
used (validation doc's "Excess-loss investigation" section) -- cuts cost
~8x relative to the full 8-call production scripts.

Two independent 1-D sweeps, not a full 2-D grid (cost):
  1. TAPER convergence: sweep N_MODES_TAPER at N_MODES_BODY held fixed.
  2. BODY convergence: sweep N_MODES_BODY at N_MODES_TAPER held fixed
     (default: the taper sweep's own best/converged value, so the body
     sweep isn't confounded by an under-resolved taper).
Select which sweep via MMI4X4_CONV_AXIS=taper|body (default taper).

RESOLUTION/N_WL/the fixed partner mode count/the sweep list are env-var
overridable for a smoke pass -- run a coarse, short list BEFORE the real
sweep, same two-step discipline as every other expensive script for this
device.
"""
import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, EMESolver, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "mode_convergence"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters -- env-var overridable, see module docstring.
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI4X4_CONV_RES", 40.0))
N_WL = int(os.environ.get("MMI4X4_CONV_NWL", 3))
N_CELLS_TAPER = int(os.environ.get("MMI4X4_CONV_NCELLS_TAPER", 30))   # settled 2026-08-07

AXIS = os.environ.get("MMI4X4_CONV_AXIS", "taper")   # "taper" or "body"
FIXED_N_MODES_BODY = int(os.environ.get("MMI4X4_CONV_FIXED_BODY", 8))     # for AXIS=taper
FIXED_N_MODES_TAPER = int(os.environ.get("MMI4X4_CONV_FIXED_TAPER", 12))  # for AXIS=body,
                        # default = the taper sweep's own settled/converged value (job
                        # 11746259), so the body sweep isn't confounded by an
                        # under-resolved taper basis
SWEEP_VALUES = [int(x) for x in
                 os.environ.get("MMI4X4_CONV_SWEEP", "4,8,12,16,20").split(",")]

WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 6.0
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 0.8
TAPER_LENGTH = 1.0
L_PI = 108.9195
DELTA_L = -1.50
MM_LENGTH = 0.75 * L_PI + DELTA_L   # = 80.1896, settled body length

wavelengths = np.linspace(1.26, 1.36, N_WL) if N_WL > 1 else np.array([1.31])

print(f"Mode-convergence sweep: axis={AXIS!r}, values={SWEEP_VALUES}, "
      f"{'N_MODES_BODY fixed=' + str(FIXED_N_MODES_BODY) if AXIS == 'taper' else 'N_MODES_TAPER fixed=' + str(FIXED_N_MODES_TAPER)}, "
      f"RES={RESOLUTION}, N_CELLS_TAPER={N_CELLS_TAPER}, {len(wavelengths)} wavelength(s)\n")

coeffs = []   # [n_sweep, n_wl] complex
wall_s = []

for val in SWEEP_VALUES:
    n_modes_taper = val if AXIS == "taper" else FIXED_N_MODES_TAPER
    n_modes_body = FIXED_N_MODES_BODY if AXIS == "taper" else val

    mmi = MMI(
        n_inputs=4, n_outputs=4,
        mm_length=MM_LENGTH, mm_width=MM_WIDTH,
        waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
        waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
        use_neighbor_aperture=True,
    )
    solver = EMESolver(
        num_modes=[1, n_modes_taper, n_modes_body, n_modes_taper, 1],
        num_cells=[1, N_CELLS_TAPER, 1, N_CELLS_TAPER, 1],
        lengths=[0.5, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, 0.5],
        oxide_width=MM_WIDTH + 3.0, oxide_height=2.0,
        res=RESOLUTION, material_stack=SOI_220, mode_solver="tidy3d",
    )

    t0 = time.time()
    result = solver.solve(
        mmi, wavelengths, ports=["out5_in0"],
        input_idx=0, output_idx=0, save_dir=save_dir,
    )
    c = result.coeff(port="out5_in0")
    dt = time.time() - t0
    coeffs.append(c)
    wall_s.append(dt)
    print(f"  {AXIS}={val:3d}: |S|^2={np.round(np.abs(c)**2, 6)}  "
          f"angle(deg)={np.round(np.angle(c, deg=True), 3)}  ({dt:.1f}s)")

coeffs = np.array(coeffs)   # [n_sweep, n_wl]
wall_s = np.array(wall_s)

np.savez(os.path.join(output_dir, f"convergence_{AXIS}.npz"),
         sweep_values=np.array(SWEEP_VALUES), wavelengths=wavelengths,
         coeffs=coeffs, wall_s=wall_s)
print(f"\nSaved -> {output_dir}/convergence_{AXIS}.npz")

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
ax = axes[0]
for wi, wl in enumerate(wavelengths):
    ax.plot(SWEEP_VALUES, np.abs(coeffs[:, wi])**2, "o-", label=f"{wl:.2f} um")
ax.set_xlabel(f"N_MODES_{AXIS.upper()}")
ax.set_ylabel("|S|^2 (out5<-in0)")
ax.set_title(f"Convergence vs N_MODES_{AXIS.upper()}")
ax.legend(fontsize=8)
ax.grid(True)

ax = axes[1]
ax.plot(SWEEP_VALUES, wall_s, "o-", color="darkorange")
ax.set_xlabel(f"N_MODES_{AXIS.upper()}")
ax.set_ylabel("Wall time (s)")
ax.set_title("Measured cost")
ax.grid(True)

plt.suptitle(f"4x4 MMI 90 hybrid -- {AXIS} mode-count convergence (port 1 -> out5, EME)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, f"convergence_{AXIS}.png"), dpi=300)
plt.close("all")
print(f"Plot saved -> {output_dir}/convergence_{AXIS}.png")
