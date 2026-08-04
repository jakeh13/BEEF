"""
Tutorial: simulate a high-contrast-grating (HCG) waveguide mirror with the
BEEF library using FDTD only — no EME, no HybridSolver.

Geometry:
    - Loaded directly from examples/gds_files/HCG_seeded_TO.gds (layer (2, 0))
      via beef.GdsGeometry — a single 500 nm-wide access waveguide feeding a
      ~5.4 x 6.0 um freeform, topology-optimization-seeded region built from
      concentric high-index-contrast Si arcs in oxide.
    - SOI platform: 160 nm Si on SiO2 (SOI_160).
    - This is a ONE-PORT device: there is no second waveguide on the far
      side to measure transmission into (confirmed by inspecting the GDS —
      see validation_markdown/high_contrast_grating_validation.md step 1).
      An HCG built this way works as a broadband mirror: light sent in from
      the access waveguide is meant to come back out the SAME port.

Method: FDTDSolver.solve_reflection() — a single EigenModeSource launches
the fundamental mode into the access waveguide; ONE mode monitor placed
between the source and the grating decomposes the total field into its
forward (incident) and backward (reflected) traveling-mode amplitudes,
giving R(lambda) = |reflected/incident|^2 in one run. This is NOT
FDTDSolver.solve() — solve() requires output_monitor_centers on the far
domain edge and only ever reads the forward mode coefficient, so it cannot
report anything physical for a device with no far-side port. See
solve_reflection()'s docstring in beef/solvers/fdtd.py for the full
forward/backward decomposition and why the domain is NOT sized the way
solve() sizes it (GdsGeometry does not extend the imported waveguide for
stub_length, unlike every parametric geometry in this repo).

Reference — the HCG reflector design itself (the seed this GDS is named for):
    J. M. Hiesener, R. P. Pesch, and S. E. Ralph, "On Topology Optimization
    Strategies for Ultra-Compact High Contrast Grating Design," IEEE
    Photonics Technology Letters, vol. 36, no. 9, pp. 597-600, 2024,
    doi: 10.1109/LPT.2024.3372429.

Run from any directory:
    python /path/to/examples/high_contrast_grating/python/high_contrast_grating_fdtd.py

Fast smoke test (2-D, low resolution, few wavelengths, loose DFT decay
tolerance — minutes on a login node, NOT a physical result):
    HCG_TEST_MODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \\
        python .../high_contrast_grating_fdtd.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection

from beef import GdsGeometry, FDTDSolver, SOI_160, load_gds_polygons

_HERE = os.path.dirname(os.path.abspath(__file__))
_GDS_PATH = os.path.abspath(os.path.join(_HERE, "..", "..", "gds_files", "HCG_seeded_TO.gds"))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "high_contrast_grating", "high_contrast_grating_fdtd"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# 0. Test-mode knob — HCG_TEST_MODE=1 runs in minutes (2-D, low res, few
#    wavelengths, loose DFT decay). Mirrors ring_resonator_hybrid.py's
#    RING_TEST_MODE convention. NOT a converged/physical result either way
#    until the resolution and decay_tol are validated — see the validation
#    markdown's "what remains" section.
# ---------------------------------------------------------------------------
TEST_MODE = bool(int(os.environ.get("HCG_TEST_MODE", "0")))

CORE_LAYERS = [(2, 0)]

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("HCG_FDTD_RES", 15.0 if TEST_MODE else 60.0))  # pts/um; MEASURED
    # (validation_markdown/high_contrast_grating_validation.md step 8): single-wavelength (1.569 um,
    # the paper's own reported peak) resolution sweep 10-60 gives R = 0.5144/0.7386/0.8434/0.8962/
    # 0.9266/0.9466/0.9557/0.9705/0.9789 -- deltas shrink geometrically (~0.55-0.6x per doubling of
    # step size), and res=60's R=0.9789 matches the paper's reported C-band minimum (97.9%) to 3
    # decimals. Extrapolated asymptote lands close to the paper's peak (98.9%), so the residual above
    # res=60 is estimated, not zero -- res=60 is a cost/accuracy tradeoff (90 min for that one point on
    # 24 ranks), not a claim of exact convergence.
THREE_D = not TEST_MODE            # real 160 nm Si core in oxide vs. fast 2-D planar smoke test
WL_MIN = 1.50                      # C-band -- the paper this device is from (Hiesener, Pesch, Ralph,
WL_MAX = 1.60                      # IEEE PTL 36(9), 597-600, 2024) reports "peak reflectance of 98.9%
                                    # ... centered at 1569 nm ... minimum reflectance of 97.9% over
                                    # C-band" -- NOT the O-band every other SOI_160 example in this
                                    # repo uses (this device is a resonant mirror seeded for a specific
                                    # band, unlike a broadband splitter/coupler). 1.50-1.60 matches this
                                    # repo's own C-band convention (directional_coupler, mmi_2x2).
N_WL = 5 if TEST_MODE else 21      # one broadband pulse captures all of these in a single run

DPML = 0.5              # PML thickness, all sides
EDGE_CLEARANCE = 0.3    # plain-oxide gap between the GDS's own bounding box and PML, both sides
SOURCE_MARGIN = 0.45    # measured inward from the PML inner edge — must land inside the
MONITOR_MARGIN = 0.8    # GDS's own 1.082 um input stub (both do; see domain check below)
SIM_WIDTH = 7.5         # y-span; the grating's arcs reach y = +-2.99 um
SIM_HEIGHT = 1.5        # slab z-height reference (3-D only)
DECAY_TOL = 1e-6 if TEST_MODE else 1e-11

# ---------------------------------------------------------------------------
# 1. Load geometry from the GDS
# ---------------------------------------------------------------------------
# region_lengths is derived from the file itself (not hardcoded) — GdsGeometry
# validates it against the imported polygons' own bounding-box span.
_raw = load_gds_polygons(_GDS_PATH, layers=CORE_LAYERS)
_xs = [x for poly in _raw[CORE_LAYERS[0]] for x, _y in poly]
_span = max(_xs) - min(_xs)

grating = GdsGeometry(
    gds_path=_GDS_PATH,
    core_layers=CORE_LAYERS,
    region_lengths=[_span],
)

print(f"HCG total length (x-span of imported polygons): {grating.total_length:.4f} um")

# ---------------------------------------------------------------------------
# 2. Configure FDTD solver
# ---------------------------------------------------------------------------
solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=SIM_WIDTH,
    sim_height=SIM_HEIGHT,
    dpml=DPML,
    source_margin=SOURCE_MARGIN,
    monitor_margin=MONITOR_MARGIN,
    three_d=THREE_D,
    num_input_modes=1,
    material_stack=SOI_160,
    decay_tol=DECAY_TOL,
)

# ---------------------------------------------------------------------------
# 3. Geometry / domain plot
# ---------------------------------------------------------------------------
# plot_fdtd_geometry (beef.plot) assumes a geometry that extends its own
# silicon by stub_length to reach the PML — true for every parametric
# geometry (MMI, Taper, ...), false for GdsGeometry (see
# GdsGeometry.polygons()'s docstring), so it would draw the source/monitor
# markers in the wrong place here. This mirrors solve_reflection()'s own
# domain formula directly instead.
def plot_hcg_domain(geometry, solver, edge_clearance, output_path):
    Sx = geometry.total_length + 2 * (edge_clearance + solver.dpml)
    polys = geometry.polygons()
    fig, ax = plt.subplots(figsize=(9, 8))
    patches = [MplPolygon(p, closed=True) for p in polys.get((2, 0), [])]
    pc = PatchCollection(patches, facecolor="steelblue", edgecolor="k", linewidths=0.3, alpha=0.85)
    ax.add_collection(pc)

    src_x = -Sx / 2 + solver.dpml + solver.source_margin
    mon_x = -Sx / 2 + solver.dpml + solver.monitor_margin
    ax.axvline(src_x, color="tab:green", lw=1.5, label="source")
    ax.axvline(mon_x, color="tab:cyan", lw=1.5, label="reflection monitor")
    ax.axvline(-geometry.total_length / 2, color="tab:red", ls="--", lw=1, label="GDS bbox edge")
    ax.axvline(geometry.total_length / 2, color="tab:red", ls="--", lw=1)
    ax.axvspan(-Sx / 2, -Sx / 2 + solver.dpml, color="orange", alpha=0.15)
    ax.axvspan(Sx / 2 - solver.dpml, Sx / 2, color="orange", alpha=0.15, label="PML")
    ax.set_xlim(-Sx / 2 - 0.2, Sx / 2 + 0.2)
    ax.set_ylim(-solver.sim_width / 2 - solver.dpml - 0.2, solver.sim_width / 2 + solver.dpml + 0.2)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("High-contrast grating — FDTD domain (solve_reflection)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


plot_hcg_domain(grating, solver, EDGE_CLEARANCE, os.path.join(output_dir, "geometry.png"))
print(f"Geometry/domain plot saved -> {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

print(f"Sweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um")
print("(single broadband pulse; one mode monitor decomposes forward+backward "
      "amplitudes at every wavelength in one run)")
if TEST_MODE:
    print("HCG_TEST_MODE=1: reduced res/2-D/few wavelengths/loose decay_tol "
          "-- smoke test only, not a physical result.")

smat = solver.solve_reflection(
    grating, wavelengths, port="top", input_mode=1,
    edge_clearance=EDGE_CLEARANCE, save_dir=save_dir,
)
smat.check_passive()  # sanity: |R| must not exceed 1 at any wavelength

# ---------------------------------------------------------------------------
# 5. Print summary
# ---------------------------------------------------------------------------
R = smat.transmission(port="in", in_mode=1, out_mode=1)
wl_center = wavelengths[len(wavelengths) // 2]
idx = int(np.argmin(np.abs(wavelengths - wl_center)))

print(f"\nResults at lambda = {wl_center:.3f} um:")
print(f"  Reflection R        : {R[idx]:.4f}")
print(f"  Reflection loss     : {10*np.log10(max(R[idx], 1e-20)):.2f} dB")
print(f"  Peak R over band    : {R.max():.4f} at {wavelengths[np.argmax(R)]:.3f} um")
print(f"  Min R over band     : {R.min():.4f} at {wavelengths[np.argmin(R)]:.3f} um")
print("  (1 - R is an upper bound on radiated/scattered loss + any residual "
      "forward transmission — this script does not separately measure where "
      "the non-reflected power goes; see the validation markdown.)")

# ---------------------------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved -> {output_dir}/smatrix.pkl")

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(wavelengths, R, marker="o", ms=3, color="k", label="R (reflection into input mode)")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Reflection")
ax.set_title("High-contrast grating mirror — FDTD reflection vs wavelength")
ax.legend()
ax.grid(True)
ax.set_ylim(0, 1.05)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "reflection.png"), dpi=300)
plt.close("all")
print(f"Reflection plot saved -> {output_dir}/reflection.png")
