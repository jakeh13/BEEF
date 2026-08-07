"""4x4 MMI 90-degree hybrid -- taper_length scoping sweep.

Not one of validation_markdown/mmi_4x4_90_hybrid_validation.md's original
Steps 0-9 -- a follow-up investigation prompted by the hybrid script's own
measured 2.1-4.3 dB excess loss (Step 6c). taper_width=0.8/taper_length=1.0
were always flagged as MMI's own class-default placeholders, never measured
(Step 2). Over just 1.0 um, the access waveguide widens 0.35 -> 0.8 um --
a ~12 deg per-side taper half-angle, aggressive for an adiabatic transition
in strongly-guided 220 nm SOI -- a real candidate for the dominant loss
mechanism (radiation into modes the output monitors don't count as signal),
consistent with excess loss behaving oddly while imbalance/phase do not.

Method: pure EME (mmi_4x4_90hybrid_eme.py's own method and geometry
parameters), swept over taper_length only (taper_width fixed at 0.8 um).
ONE driven input (port 1, idx=0) only, not both ports 1/3 -- taper loss is
local physics at the narrow access-waveguide end, far from the wide body's
own mode-interference structure, so it should not differ meaningfully
between an edge port (1) and an interior port (3) the way phase/CMRR do.
Scoping cheaply to see the TREND, not claiming full port coverage -- flag
this simplification if a real refinement pass follows.

Cost is UNMEASURED before this script's own smoke pass, deliberately, per
this repo's own two-step discipline for a never-run script shape: unlike
the ΔL body sweep (Step 3), which is analytically free per length via
solve_section_length_sweep() (one mode solve per wavelength, reused across
the whole 81-point grid), EMESolver.solve()'s disk cache is fingerprinted
on the geometry's OWN drawn GDS (to_gds_eme()) per full 5-region chain call
-- changing taper_length changes that GDS, so each new taper_length is
likely a genuine new solve, not a free analytic sweep point like ΔL was.
Whether meow's own internal mode solver additionally caches by individual
local cross-section (which would make the 15 per-taper cell mode solves
reusable across different taper_length values, since the width AT EACH
FRACTIONAL POSITION along a linear taper does not depend on the taper's
overall length) is not established by reading beef's own solver code --
measured directly below, not assumed.

RESOLUTION/N_MODES_TAPER/N_WL/TAPER_LENGTHS are env-var overridable for a
smoke pass (MMI4X4_TAPER_RES / MMI4X4_TAPER_NMODES_TAPER / MMI4X4_TAPER_NWL
/ MMI4X4_TAPER_LENGTHS) -- run a small, coarse smoke subset BEFORE the full
sweep at production settings.

Run:
    python examples/verification/python/mmi_4x4_90hybrid_taper_sweep.py
"""
import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, EMESolver, SMatrix, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "taper_sweep"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters -- env-var overridable, see module docstring.
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI4X4_TAPER_RES", 40.0))
N_WL = int(os.environ.get("MMI4X4_TAPER_NWL", 1))   # centre only by default -- a scoping
                                                      # sweep, not the paper's own 3-pt Eq.12
                                                      # fitness convention (that's for the
                                                      # settled device, not this side sweep)
N_MODES_TAPER = int(os.environ.get("MMI4X4_TAPER_NMODES_TAPER", 4))   # matches
                                                                        # mmi_4x4_90hybrid_eme.py's own default
N_MODES_BODY = int(os.environ.get("MMI4X4_TAPER_NMODES_BODY", 8))   # job 11707935/11707835
                    # at the default -- same body cross-section as every other script for
                    # this device, unaffected by taper_length; override down for a smoke pass

TAPER_LENGTHS = [float(x) for x in
                  os.environ.get("MMI4X4_TAPER_LENGTHS", "1.0,2.0,3.0,4.0,5.0,7.0").split(",")]

# ---------------------------------------------------------------------------
# Geometry -- identical to mmi_4x4_90hybrid_eme.py's settled Steps 1-3
# parameters, except taper_length is now the swept variable.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 6.0
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 0.8       # held fixed -- this sweep isolates taper_length only
L_PI = 108.9195
DELTA_L = -1.50
MM_LENGTH = 0.75 * L_PI + DELTA_L   # = 80.1896, same settled body length

DRIVEN_IDX = 0   # port 1 only -- see module docstring for why this is enough to scope
OUTPUTS = [0, 1, 2, 3]
PORT_NAMES = ["out5", "out6", "out7", "out8"]

wavelengths = np.array([1.31]) if N_WL == 1 else np.linspace(1.26, 1.36, N_WL)

print(f"Sweeping taper_length over {TAPER_LENGTHS} um, taper_width={TAPER_WIDTH} um fixed, "
      f"driven port 1 (idx=0) only, {len(wavelengths)} wavelength(s), "
      f"RES={RESOLUTION}, N_MODES_TAPER={N_MODES_TAPER}, N_MODES_BODY={N_MODES_BODY}\n")

# ---------------------------------------------------------------------------
# Sweep -- 4 solve() calls (one per output) per taper_length, no cross-length
# cache reuse assumed (measured, see per-point timing printed below).
# ---------------------------------------------------------------------------
excess_loss_pts = []
imbalance_pts = []
wall_s_pts = []

for tl in TAPER_LENGTHS:
    mmi = MMI(
        n_inputs=4, n_outputs=4,
        mm_length=MM_LENGTH, mm_width=MM_WIDTH,
        waveguide_width=WAVEGUIDE_WIDTH, taper_length=tl, taper_width=TAPER_WIDTH,
        waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
        use_neighbor_aperture=True,
    )
    solver = EMESolver(
        num_modes=[1, N_MODES_TAPER, N_MODES_BODY, N_MODES_TAPER, 1],
        num_cells=[1, 15, 1, 15, 1],
        lengths=[0.5, tl, MM_LENGTH, tl, 0.5],
        oxide_width=MM_WIDTH + 3.0,
        oxide_height=2.0,
        res=RESOLUTION,
        material_stack=SOI_220,
        mode_solver="tidy3d",
    )

    smat = SMatrix(wavelengths=wavelengths)
    t0 = time.time()
    for j in OUTPUTS:
        label = f"{PORT_NAMES[j]}_in{DRIVEN_IDX}"
        result = solver.solve(
            mmi, wavelengths, ports=[label],
            input_idx=DRIVEN_IDX, output_idx=j, save_dir=save_dir,
        )
        smat._data[(PORT_NAMES[j], 1, 1)] = result.coeff(port=label)
    dt = time.time() - t0

    excess_loss = smat.excess_loss_db(PORT_NAMES)
    imbalance = smat.imbalance_db(PORT_NAMES)
    excess_loss_pts.append(excess_loss)
    imbalance_pts.append(imbalance)
    wall_s_pts.append(dt)
    print(f"  taper_length={tl:5.2f} um: excess_loss={np.round(excess_loss, 3)} dB  "
          f"imbalance={np.round(imbalance, 3)} dB  ({dt:.1f}s, 4 calls)")

excess_loss_arr = np.array(excess_loss_pts)   # [n_taper, n_wl]
imbalance_arr = np.array(imbalance_pts)
wall_s_arr = np.array(wall_s_pts)

# ---------------------------------------------------------------------------
# Save + plot.
# ---------------------------------------------------------------------------
np.savez(os.path.join(output_dir, "taper_sweep.npz"),
         taper_lengths=np.array(TAPER_LENGTHS), wavelengths=wavelengths,
         excess_loss=excess_loss_arr, imbalance=imbalance_arr, wall_s=wall_s_arr)
print(f"\nSweep grid saved -> {output_dir}/taper_sweep.npz")

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
ax = axes[0]
for wi, wl in enumerate(wavelengths):
    ax.plot(TAPER_LENGTHS, excess_loss_arr[:, wi], "o-", label=f"{wl:.2f} um")
ax.set_xlabel("taper_length (um)")
ax.set_ylabel("Excess loss, port 1 (dB)")
ax.set_title("Excess Loss vs Taper Length (Eq. 11)")
ax.legend(fontsize=8)
ax.grid(True)

ax = axes[1]
ax.plot(TAPER_LENGTHS, wall_s_arr, "o-", color="darkorange")
ax.set_xlabel("taper_length (um)")
ax.set_ylabel("Wall time, 4 calls (s)")
ax.set_title("Measured Cost vs Taper Length")
ax.grid(True)

plt.suptitle("4x4 MMI 90 hybrid -- taper_length scoping sweep (port 1, EME-only)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "taper_sweep.png"), dpi=300)
plt.close("all")
print(f"Sweep plot saved -> {output_dir}/taper_sweep.png")
