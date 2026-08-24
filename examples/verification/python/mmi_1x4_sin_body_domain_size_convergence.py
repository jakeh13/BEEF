"""Domain-size (sim_width, sim_height) convergence check for the 1x4 MMI
body cross-section on 300 nm Si3N4 -- prompted by a direct question ("can
the whole-device FDTD sim domain be reduced in y or z?") while
mmi_1x4_fdtd.py's own real cost was projected at 52-390 hrs, dominated by
domain size (y,z scale cell count linearly each; x is fixed by device
length and not addressed here).

Both current values are INHERITED conventions, never actually validated
for necessity on this device:
  - sim_width = MM_WIDTH + 8.0 = 20.0 um -- the same "+8.0" margin every
    other FDTD section for this device uses (input/output/calibration),
    itself never derived for this platform's own near-cutoff body modes
    (job 11660226 found guided modes down to margin 0.004 above n_clad --
    about as weakly confined as a "guided" mode gets).
  - sz = 5.0 um -- validated (step 4a, mmi_1x4_sin_width1p1_check.py, job
    11660166) but ONLY for the single 1.1 um ACCESS WAVEGUIDE cross-section,
    never for the 12.0 um-wide BODY cross-section's own near-cutoff modes,
    which is the more demanding case (same reasoning
    mmi_1x4_sin_body_res_convergence.py's own docstring already used for
    why the body, not the access waveguide, is the right thing to track).

Method: same convergence-sweep pattern as mmi_1x4_sin_body_res_convergence.py
(itself modeled on step 4's sz check) -- sweep WIDTH (at fixed sz=5.0) and
SZ (at fixed width=20.0) independently, watch |dn_eff| for the same
tracked-mode spread (fundamental through near-cutoff), find where each
stops moving. Cheap: tidy3d mode solves only, no FDTD time-stepping --
seconds to minutes per point, unlike the FDTD probe this is meant to avoid
re-running blind.

compute_modes_tidy3d, not compute_modes_meep -- same fix as every other
mode-solve script for this device (step 4).

Run:
    python /path/to/examples/verification/python/mmi_1x4_sin_body_domain_size_convergence.py
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
    _HERE, "..", "..", "..", "example_output", "verification",
    "mmi_1x4_sin_body_domain_size_convergence"))
os.makedirs(output_dir, exist_ok=True)

WG_HEIGHT = 0.30
WL = 1.55
MM_WIDTH = 12.0
RES = 40.0   # cheap either way -- use the higher-fidelity value
N_CLAD = float(sellmeier_index("sio2", WL))
N_MODES = 23   # job 11660226 found 23 guided at mm_width=12.0

SZ_DEFAULT = 5.0             # current production value (step 4a, access waveguide only)
MESH_WIDTH_DEFAULT = MM_WIDTH + 8.0   # current production value (inherited "+8.0" convention)

CONVERGENCE_TOL = 1e-4   # same tolerance every other convergence sweep in this doc uses
TRACK_INDICES = [0, 10, 19, 22]   # fundamental, mid-band, near-cutoff, most-cutoff (same as res sweep)


def _axis_grid(half_span: float, res: float) -> np.ndarray:
    n_half = int(round(half_span * res))
    step = 1.0 / res
    return step * np.arange(-n_half, n_half + 1)


def build_body_cross_section(mesh_width: float, sz: float, res: float = RES) -> mw.CrossSection:
    stack = MaterialStack(wg_height=WG_HEIGHT, sz=sz, core_kind="si3n4")
    length = 1.0
    wg_poly = [
        (0, -MM_WIDTH / 2), (length, -MM_WIDTH / 2),
        (length, MM_WIDTH / 2), (0, MM_WIDTH / 2),
    ]
    ox_poly = [
        (0, -mesh_width / 2), (length, -mesh_width / 2),
        (length, mesh_width / 2), (0, mesh_width / 2),
    ]
    gds_polys = {(2, 0): [wg_poly], (1, 0): [ox_poly]}
    structs = extrude_polygons(gds_polys, stack.extrusion_rules)
    mesh = mw.Mesh2D(x=_axis_grid(mesh_width / 2, res), y=_axis_grid(sz / 2, res))
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=float(WL), T=25.0)
    return mw.CrossSection.from_cell(cell=cells[0], env=env)


def sweep(label, candidates, build_fn):
    print("=" * 78)
    print(f"{label} sweep -- mm_width={MM_WIDTH}, lambda={WL}, tracking modes {TRACK_INDICES}")
    print("=" * 78)
    neffs_by_val = np.zeros((len(candidates), N_MODES))
    for i, val in enumerate(candidates):
        cs = build_fn(val)
        modes = compute_modes_tidy3d(cs, num_modes=N_MODES)
        neffs = sorted((float(np.real(m.neff)) for m in modes), reverse=True)
        neffs_by_val[i] = neffs
        guided = [n for n in neffs if n > N_CLAD]
        tracked = {idx: (neffs[idx] if idx < len(neffs) else float("nan")) for idx in TRACK_INDICES}
        print(f"  {label}={val:6.2f}  ({len(guided)} guided of {len(neffs)})  "
              + "  ".join(f"mode{idx}={tracked[idx]:.6f}" for idx in TRACK_INDICES))

    print(f"\nSuccessive |dn_eff| per tracked mode ({label}):")
    converged_val = {}
    for idx in TRACK_INDICES:
        col = neffs_by_val[:, idx]
        d = np.abs(np.diff(col))
        print(f"  mode {idx}: {[f'{x:.2e}' for x in d]}")
        conv_i = None
        for i, x in enumerate(d):
            if x < CONVERGENCE_TOL:
                conv_i = i
                break
        converged_val[idx] = candidates[conv_i] if conv_i is not None else None
        if conv_i is not None:
            print(f"    -> CONVERGED at {label}={candidates[conv_i]:.2f} "
                  f"(|dn_eff| < {CONVERGENCE_TOL:.0e} going to {candidates[conv_i+1]:.2f})")
        else:
            print(f"    -> NOT CONVERGED within swept range")

    worst = [v for v in converged_val.values() if v is not None]
    overall = max(worst) if len(worst) == len(TRACK_INDICES) else None
    print(f"\nOverall (worst-case) convergence {label}: "
          f"{overall if overall is not None else 'NOT REACHED within swept range'}")
    return neffs_by_val, converged_val, overall


# ---------------------------------------------------------------------------
# Sweep 1: sim_width (mesh_width / oxide_width), sz held at production 5.0
# ---------------------------------------------------------------------------
WIDTH_CANDIDATES = [13.0, 14.0, 16.0, 18.0, 20.0]   # production is 20.0
neffs_width, conv_width, overall_width = sweep(
    "mesh_width", WIDTH_CANDIDATES, lambda w: build_body_cross_section(w, SZ_DEFAULT))

# ---------------------------------------------------------------------------
# Sweep 2: sz (vertical extent), mesh_width held at production 20.0 --
# re-checks step 4a's own sz=5.0 result against the BODY cross-section
# specifically, not just the access waveguide.
# ---------------------------------------------------------------------------
SZ_CANDIDATES = [2.0, 3.0, 4.0, 5.0]   # production is 5.0; step 4a converged at 5.0 for the access wg
neffs_sz, conv_sz, overall_sz = sweep(
    "sz", SZ_CANDIDATES, lambda s: build_body_cross_section(MESH_WIDTH_DEFAULT, s))

print(f"\nCurrently used in production: mesh_width={MESH_WIDTH_DEFAULT}, sz={SZ_DEFAULT}")

np.savez(
    os.path.join(output_dir, "domain_size_convergence.npz"),
    width_candidates=np.asarray(WIDTH_CANDIDATES), neffs_width=neffs_width,
    sz_candidates=np.asarray(SZ_CANDIDATES), neffs_sz=neffs_sz,
    track_indices=np.asarray(TRACK_INDICES), mm_width=MM_WIDTH,
    mesh_width_default=MESH_WIDTH_DEFAULT, sz_default=SZ_DEFAULT,
)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for idx in TRACK_INDICES:
    axes[0].plot(WIDTH_CANDIDATES, neffs_width[:, idx], marker="o", label=f"mode {idx}")
axes[0].axhline(N_CLAD, color="gray", linestyle=":", lw=0.8, label="n_clad")
axes[0].axvline(MESH_WIDTH_DEFAULT, color="k", linestyle="--", lw=1.0, label=f"currently used ({MESH_WIDTH_DEFAULT})")
axes[0].set_xlabel("mesh_width / sim_width (um)")
axes[0].set_ylabel("n_eff")
axes[0].set_title("1x4 MMI body -- Y-domain (sim_width) convergence")
axes[0].legend(fontsize=8)
axes[0].grid(True)

for idx in TRACK_INDICES:
    axes[1].plot(SZ_CANDIDATES, neffs_sz[:, idx], marker="o", label=f"mode {idx}")
axes[1].axhline(N_CLAD, color="gray", linestyle=":", lw=0.8, label="n_clad")
axes[1].axvline(SZ_DEFAULT, color="k", linestyle="--", lw=1.0, label=f"currently used ({SZ_DEFAULT})")
axes[1].set_xlabel("sz (um)")
axes[1].set_ylabel("n_eff")
axes[1].set_title("1x4 MMI body -- Z-domain (sz) convergence")
axes[1].legend(fontsize=8)
axes[1].grid(True)

plt.tight_layout()
fig.savefig(os.path.join(output_dir, "domain_size_convergence.png"), dpi=300)
plt.close("all")
print(f"\nPlot saved -> {output_dir}/domain_size_convergence.png")
