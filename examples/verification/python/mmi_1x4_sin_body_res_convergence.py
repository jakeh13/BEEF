"""Mesh resolution convergence check for the 1x4 MMI body cross-section on
300 nm Si3N4 -- prompted by the question "is res=40 too fine for nitride,
given its weaker index contrast?" while the real mmi_1x4_eme.py/
mmi_1x4_hybrid.py prewarm jobs (11672738/11672739) were already running at
res=40, this codebase's standing silicon-platform default, never itself
convergence-checked for nitride (flagged, not resolved, in
mmi_1x4_sin_cross_section.py's own RES comment and validation doc step 8
check 3).

Physical reasoning for why this might matter: nitride's core/clad index
contrast is ~5.3x weaker than Si/SiO2's (validation doc A7), so mode fields
vary more smoothly across the core boundary and may need a coarser mesh to
resolve to the same accuracy than a silicon platform would.

Targets the BODY cross-section specifically (mm_width=12.0 um), not a
simple access waveguide -- the body's near-cutoff guided modes (job
11660226: margins down to 0.004 above n_clad for the bottom of the
N_BODY=20 kept band) are the modes most likely to be resolution-sensitive,
same reasoning weakly-confined structures needed the sz convergence check
in step 4 (mmi_1x4_sin_cross_section.py). A simple single-waveguide
fundamental-mode check would likely converge trivially and miss this.

Method: same convergence-sweep pattern as step 4's sz check
(mmi_1x4_sin_cross_section.py) -- sweep RES, watch |dn_eff| for a handful of
modes spanning from comfortably-guided to near-cutoff, find where it stops
moving.

compute_modes_tidy3d, not compute_modes_meep -- step 4's fix, unchanged
reasoning here.

Run:
    python /path/to/examples/verification/python/mmi_1x4_sin_body_res_convergence.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import meow as mw
from meow import compute_modes_tidy3d

from beef import extrude_polygons
from beef.materials.stack import MaterialStack, sellmeier_index

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification", "mmi_1x4_sin_body_res_convergence"))
os.makedirs(output_dir, exist_ok=True)

WG_HEIGHT = 0.30
WL = 1.55
SZ = 5.0            # step 4's converged value
MM_WIDTH = 12.0      # the body's own width (step 5)
N_CLAD = float(sellmeier_index("sio2", WL))

STACK = MaterialStack(wg_height=WG_HEIGHT, sz=SZ, core_kind="si3n4")

MESH_WIDTH = MM_WIDTH + 8.0   # matches the real solvers' oxide_width convention
N_MODES = 23                  # job 11660226 found 23 guided at mm_width=12.0

RES_CANDIDATES = [15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
RES_CONVERGENCE_TOL = 1e-4   # same tolerance step 4's sz sweep used

# Track a spread of modes: fundamental (comfortably guided), a mid-band
# mode, and the near-cutoff bottom of the kept N_BODY=20 band (indices
# 0, 10, 19 -- 0-indexed, descending n_eff) plus the very bottom of what's
# guided at all (index 22) as a stress case.
TRACK_INDICES = [0, 10, 19, 22]


def _axis_grid(half_span: float, res: float) -> np.ndarray:
    n_half = int(round(half_span * res))
    step = 1.0 / res
    return step * np.arange(-n_half, n_half + 1)


def build_body_cross_section(res: float) -> mw.CrossSection:
    length = 1.0
    wg_poly = [
        (0, -MM_WIDTH / 2), (length, -MM_WIDTH / 2),
        (length, MM_WIDTH / 2), (0, MM_WIDTH / 2),
    ]
    ox_poly = [
        (0, -MESH_WIDTH / 2), (length, -MESH_WIDTH / 2),
        (length, MESH_WIDTH / 2), (0, MESH_WIDTH / 2),
    ]
    gds_polys = {(2, 0): [wg_poly], (1, 0): [ox_poly]}
    structs = extrude_polygons(gds_polys, STACK.extrusion_rules)
    mesh = mw.Mesh2D(x=_axis_grid(MESH_WIDTH / 2, res), y=_axis_grid(SZ / 2, res))
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=float(WL), T=25.0)
    return mw.CrossSection.from_cell(cell=cells[0], env=env)


print("=" * 78)
print(f"Body cross-section resolution convergence -- mm_width={MM_WIDTH} um, "
      f"sz={SZ} um, lambda={WL} um, tracking modes {TRACK_INDICES}")
print(f"n_clad (SiO2 @ {WL} um) = {N_CLAD:.6f}")
print("=" * 78)

neffs_by_res = np.zeros((len(RES_CANDIDATES), N_MODES))
for i, res in enumerate(RES_CANDIDATES):
    cs = build_body_cross_section(res)
    modes = compute_modes_tidy3d(cs, num_modes=N_MODES)
    neffs = sorted((float(np.real(m.neff)) for m in modes), reverse=True)
    neffs_by_res[i] = neffs
    guided = [n for n in neffs if n > N_CLAD]
    tracked = {idx: (neffs[idx] if idx < len(neffs) else float("nan")) for idx in TRACK_INDICES}
    print(f"  res={res:5.1f} pts/um  ({len(guided)} guided of {len(neffs)})  "
          + "  ".join(f"mode{idx}={tracked[idx]:.6f}" for idx in TRACK_INDICES))

print("\nSuccessive |dn_eff| per tracked mode:")
converged_res = {}
for idx in TRACK_INDICES:
    col = neffs_by_res[:, idx]
    d = np.abs(np.diff(col))
    print(f"  mode {idx}: {[f'{x:.2e}' for x in d]}")
    conv_i = None
    for i, x in enumerate(d):
        if x < RES_CONVERGENCE_TOL:
            conv_i = i
            break
    converged_res[idx] = RES_CANDIDATES[conv_i] if conv_i is not None else None
    if conv_i is not None:
        print(f"    -> CONVERGED at res={RES_CANDIDATES[conv_i]:.1f} "
              f"(|dn_eff| < {RES_CONVERGENCE_TOL:.0e} going to "
              f"{RES_CANDIDATES[conv_i+1]:.1f})")
    else:
        print(f"    -> NOT CONVERGED within swept range")

worst = [v for v in converged_res.values() if v is not None]
overall = max(worst) if len(worst) == len(TRACK_INDICES) else None
print(f"\nOverall (worst-case across tracked modes) convergence res: "
      f"{overall if overall is not None else 'NOT REACHED within swept range'}")
print(f"Currently used in the real solvers: res=40.0")

np.savez(
    os.path.join(output_dir, "res_convergence.npz"),
    res_candidates=np.asarray(RES_CANDIDATES), neffs_by_res=neffs_by_res,
    track_indices=np.asarray(TRACK_INDICES), mm_width=MM_WIDTH, sz=SZ,
)

fig, ax = plt.subplots(figsize=(7, 5))
for idx in TRACK_INDICES:
    ax.plot(RES_CANDIDATES, neffs_by_res[:, idx], marker="o", label=f"mode {idx}")
ax.axhline(N_CLAD, color="gray", linestyle=":", lw=0.8, label="n_clad")
ax.axvline(40.0, color="k", linestyle="--", lw=1.0, label="currently used (res=40)")
ax.set_xlabel("Mesh resolution (pts/um)")
ax.set_ylabel("n_eff")
ax.set_title("1x4 MMI body -- resolution convergence")
ax.legend(fontsize=8)
ax.grid(True)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "res_convergence.png"), dpi=300)
plt.close("all")
print(f"\nPlot saved -> {output_dir}/res_convergence.png")
