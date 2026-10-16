"""Tutorial: 4x4 MMI 90-degree optical hybrid -- pure EME simulation (BEEF
method, no FDTD).

Same device as
mmi_4x4_90hybrid_hybrid.py (see that script's module docstring for the
paper reference, the optimised design point -- mm_width=8.5,
taper 1.30/4.0, mm_length=159.30 -- and what it does and does not meet).

Method: Eigenmode Expansion via meow, full 5-region device chain (wg_in |
taper_in | MMI body | taper_out | wg_out), matching mmi_2x2_eme.py's own
pattern -- EMESolver.solve(geometry, wavelengths, ports=[...], input_idx=,
output_idx=, ...) restricts the drawn geometry to ONE input port and ONE
output port per call (meow's to_gds_eme(input_idx=, output_idx=) draws only
those two branches, not all 4+4 simultaneously), so a 2-port device needs 2
calls (bar/cross); THIS 4-port device needs one call per (driven input,
measured output) pair.

Cost, up front, since this differs from every prior EME-only example in
this repo: EMESolver.solve() has NO cross-call cache reuse for shared
geometry (unlike the hybrid script's mirror-symmetry trick) -- each of the
8 calls below (2 driven inputs (paper's ports 1/3, finding A4) x 4 outputs)
independently re-solves EVERY cell in the 5-region chain from scratch,
including the wide MMI body (measured at ~416 s/mode for THIS body's own
cross-section, Step 3a) and 15+15 taper sub-cells per call. RESOLUTION and
N_MODES_TAPER are env-var overridable (MMI4X4_EME_RES / MMI4X4_EME_NMODES_
TAPER) for a smoke pass -- run one BEFORE the real 8-call sweep, per this
repo's own two-step discipline for expensive, never-run scripts (today's
own session already hit two real cost surprises building the hybrid
script's sibling; do not assume this one is cheap just because it has no
FDTD time-stepping).

Approximation this script makes that the hybrid script does NOT: each
solve() call only ever draws the ONE input port and ONE output port it was
asked for, treating the other 6 ports as if they were physically absent
(no crosstalk with them is modelled at all) -- the same simplification
mmi_2x2_eme.py's own module docstring names ("comfortably decoupled" at its
1.05 um port-to-port gap). This device's gap is narrower (offset=1.5,
taper_width=0.8 -> 0.7 um clearance, Step 2), so this approximation is
LESS obviously safe here than for the 2x2 -- flagged, not verified in this
script. The hybrid script's real FDTD sections DO include all 4 ports'
real geometry simultaneously and are the more physically complete model;
this EME-only script is the cheaper, approximate one, by design.

N_MODES_TAPER bumped 4->12, 2026-08-08, at the user's direction after
visually inspecting this script's own field plots (examples/verification/
python/plot_mmi4x4_90hybrid_eme_fields.py's output): at N_MODES_TAPER=4
the launched mode stayed unphysically confined and speckly well into the
MMI body instead of fanning out into a clean self-imaging interference
pattern -- a visibly under-resolved taper mode basis, confirmed with a
same-settings (RES=40/num_cells=30/N_MODES_BODY=8) before/after field-plot
comparison (job 11745928 for the N_MODES_TAPER=12 render). Cost checked
first, not assumed: coarse full-8-call timing (RES=15/N_MODES_BODY=4/N_WL=1)
was 121s at N_MODES_TAPER=4 vs 221s at 12 -- only ~1.83x, consistent with
the num_cells investigation's own finding that isolated-call cost
comparisons are misleading here (JIT-compile cost is paid once and
amortized across the 8 same-shape calls) -- applied to the num_cells=30
production run's own 73:27, giving an estimated ~2.2h, well inside the
existing -t 360 batch budget.

Run:
    python /path/to/examples/mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_eme.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, EMESolver, SMatrix, SOI_220, plot_eme_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
# Output directory carries N_WL when it is not the default 3. Fourth instance of
# the same hazard the _tl / _n / _wl filename tags exist for (see wg_eval.py):
# a denser wavelength sweep is a DIFFERENT measurement of the same geometry, and
# these scripts write fixed filenames, so without this a denser run would
# silently destroy the 3-point results that were verified on 2026-09-06 to
# reproduce the optimiser exactly. Default stays untagged so every existing
# path and reference keeps working.
_NWL_ENV = int(os.environ.get("MMI4X4_EME_NWL", 3))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid",
    "mmi_4x4_90hybrid_eme" + ("" if _NWL_ENV == 3 else f"_nwl{_NWL_ENV}")))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters -- RESOLUTION/N_MODES_TAPER/N_WL are env-var
# overridable for a smoke pass (see module docstring).
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI4X4_EME_RES", 40.0))
# BAND: 1.31 +- 0.02 um -- the band the geometry above was OPTIMISED over.
# These scripts previously ran 1.26-1.36 (+-50 nm), which is 2.5x wider than
# the design band, and this device is steeply wavelength-dependent: measured
# 2026-09-05, min CMRR falls 20.761 -> 18.482 dBe and excess loss rises
# 0.654 -> 2.706 dB going from 1.310 um alone to +-20 nm (job 12831929).
# Reporting the optimised geometry over +-50 nm would judge it far outside the
# band it was designed for. The old window is still reachable via the env
# overrides for anyone who wants the wider sweep.
WL_MIN = float(os.environ.get("MMI4X4_WLMIN", 1.29))
WL_MAX = float(os.environ.get("MMI4X4_WLMAX", 1.33))
N_WL = _NWL_ENV   # matches the paper's own Eq. 12 convention
                                                    # (Sec. IV) and the hybrid script's own N_WL

N_MODES_TAPER = int(os.environ.get("MMI4X4_EME_NMODES_TAPER", 12))   # was 4 (PLACEHOLDER,
                    # half of N_MODES_BODY, not independently tuned) -- bumped 2026-08-08,
                    # see module docstring: field plots showed 4 was visibly under-resolved

# ---------------------------------------------------------------------------
# 1. Geometry -- identical to mmi_4x4_90hybrid_hybrid.py's settled Steps 1-3
#    parameters (same device, different simulation method).
# ---------------------------------------------------------------------------
# Geometry -- the OPTIMISED design point, not Step 2's first pass.
# Result of the (mm_width, y_gap, taper_width, taper_length, mm_length)
# campaign of ~45 hybrid-path evaluations, rounds 8-10, 2026-08-28..09-05). Measured at
# this point on the 3-point {1.29, 1.31, 1.33} um fitness grid the campaign used:
# max|dphi| = 3.933 deg (spec <=5, MET), min CMRR = 18.482 dBe (spec >=20, 1.52
# short), excess loss 2.706 dB, imbalance 1.972 dB. On this script's 11-point grid
# the phase error is 6.223 deg at 1.322 um -- the spec is NOT met there. Superseded first pass, for the record: mm_width=6.0,
# taper_width=0.8, taper_length=1.0, mm_length=78.9396 -- which measured
# 4.764 deg / 12.109 dBe / 4.076 dB / 5.202 dB on the same evaluator.
#
# taper_width tracks mm_width: the campaign's most reliable finding is that a
# taper which has not widened with the body cannot be rescued by ANY length
# (0 of 7601 lengths met the phase spec at mm_width=8.5 with taper 1.20).
# Measured optima: 6.5 -> 1.065, 7.0 -> 1.10, 7.5 -> 1.20, 8.0 -> 1.20,
# 8.5 -> 1.30.
WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 8.5
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 1.30            # MEASURED optimum at this width (round 10, job 12783704)
TAPER_LENGTH = 4.0            # MEASURED (round 8): tl=1->4 cut excess loss 4.076->2.940 dB;
                              # tl=8 bought only 0.006 more objective for 2x the footprint
L_PI = 218.5954   # = 108.9195 * (8.5/6.0)**2 -- SCALED, not re-measured. L_pi/W^2 was
                  # measured constant to 0.1% over 3.5-16 um (Step 1c), so the mm_width=6.0
                  # measurement (job 11693043) carries.
# MM_LENGTH: the campaign's own optimum for THIS geometry, selected from the
# shared absolute 20-400 um / 0.05 um grid (mmi_4x4_90hybrid_ygap_grid.py's
# MM_LENGTH_GRID), not from a +-10 um DELTA_L window around 3/4*L_pi. The
# length re-score (job 12832166) confirmed this is simultaneously the best
# length under all three selection rules tried -- the live objective,
# phase-as-constraint, and maximise-CMRR-subject-to-phase -- so nothing is
# being traded away by pinning it. Only 19 of 7601 lengths meet the 5 deg
# phase spec on the 3-point fitness grid at all, so the grid-snap assert below matters more here than it
# did at the first-pass point.
MM_LENGTH = 159.30
_LENGTH_STEP = 0.05
assert abs(round(MM_LENGTH / _LENGTH_STEP) * _LENGTH_STEP - MM_LENGTH) < 1e-9, (
    f"MM_LENGTH={MM_LENGTH} is not an exact point on the {_LENGTH_STEP} um "
    f"MM_LENGTH_GRID the optimum was selected from -- do not silently snap."
)
DELTA_L = MM_LENGTH - 0.75 * L_PI   # = -4.6465, informational only
# DELTA_L/N_MODES_BODY CORRECTED 2026-08-21, and this was a real staleness
# bug, not a tidy-up. The 2026-08-09 body-mode redo (N_MODES_BODY 8->16,
# which shifted the optimum to DELTA_L=-2.75/mm_length=78.9396) was applied
# to body_sweep.py and hybrid.py but NOT to this script or to
# verification/python/plot_mmi4x4_90hybrid_eme_fields.py. Both sat at the
# superseded 8-mode/-1.50 point for 12 days while the doc's own CURRENT
# STATUS header claimed "all 3 example scripts done and RE-RUN at a
# properly-converged body mode basis" -- its very next line scoped the fix
# to "(body_sweep.py/hybrid.py)", which is exactly the two that got it.
# Nothing recorded leaving this script behind deliberately. Every EME
# artifact on disk (eme_results.png, geometry.png, the six eme_E*/H* field
# plots) therefore showed a design point the project had already rejected,
# at a basis the convergence check measured as 50-64% under-converged.

N_MODES_BODY = 16   # was 8 (job 11707935/11707835) -- see the note above and
                    # mmi_4x4_90hybrid_hybrid.py's own comment

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
    use_neighbor_aperture=True,
)
print(f"MMI total length: {mmi.total_length:.4f} µm (mm_length={mmi.mm_length:.4f})")
print(f"Input ports at y = {mmi.input_ys} µm; outputs at y = {mmi.output_ys} µm")

# ---------------------------------------------------------------------------
# 2. EME solver -- 5 regions: wg_in | taper_in | MMI body | taper_out | wg_out.
# ---------------------------------------------------------------------------
N_CELLS_TAPER = int(os.environ.get("MMI4X4_EME_NCELLS_TAPER", 30))
# Taper staircase cell count -- was 15 (mmi_2x2_eme.py's own convention),
# bumped to 30 2026-08-07 after the taper-length scoping sweep found
# num_cells=15 measurably under-converged (a fixed-num_cells staircase's
# per-step width-jump doesn't depend on taper_length, so the discretization
# itself -- not physical taper length -- was the dominant loss driver;
# measured in a dedicated excess-loss investigation). 30, not the 60-100
# the convergence check used,
# per the user's own direction ("err on the lower side as 60 cells is a
# lot") -- captures roughly half the measured 15->100 gain at a 2x cost
# step instead of 4-7x. Scoped to THIS script only: the hybrid/body-sweep
# scripts use FDTD (continuous geometry, no staircase) for their tapers.

solver = EMESolver(
    num_modes=[1, N_MODES_TAPER, N_MODES_BODY, N_MODES_TAPER, 1],
    num_cells=[1, N_CELLS_TAPER, 1, N_CELLS_TAPER, 1],   # tapers subdivided; body
                                    # uniform, no subdivision needed
    lengths=[0.5, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, 0.5],
    oxide_width=MM_WIDTH + 3.0,    # matches the hybrid script's own FDTD section margin
    oxide_height=2.0,              # Step 1b's converged sz at this platform/band
    res=RESOLUTION,
    material_stack=SOI_220,
    mode_solver="tidy3d",          # no FDTD section in this script to stay consistent
                                    # with -- matches mmi_2x2_eme.py's own reasoning
)

plot_eme_geometry(
    mmi, solver=solver,
    region_labels=["wg in", "taper in", "MMI body", "taper out", "wg out"],
    title="4×4 MMI 90° hybrid — EME — Top-down Geometry",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved → {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 3. Sweep -- 2 driven inputs (paper's ports 1/3, A4) x 4 outputs (paper's
#    ports 5-8), 8 independent full-chain solves, no cross-call reuse.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
DRIVEN = [0, 2]
OUTPUTS = [0, 1, 2, 3]
PORT_NAMES = ["out5", "out6", "out7", "out8"]

smat = {i: SMatrix(wavelengths=wavelengths) for i in DRIVEN}
for i in DRIVEN:
    for j in OUTPUTS:
        label = f"{PORT_NAMES[j]}_in{i}"
        result = solver.solve(
            mmi, wavelengths, ports=[label],
            input_idx=i, output_idx=j, save_dir=save_dir,
        )
        smat[i]._data[(PORT_NAMES[j], 1, 1)] = result.coeff(port=label)
        print(f"  input_idx={i} -> output_idx={j} solved")
    smat[i].check_passive()

# ---------------------------------------------------------------------------
# 4. Metrics (paper Eqs. 10-15), same SMatrix methods as the hybrid script.
#    A11's phase-sign fix applies identically here.
# ---------------------------------------------------------------------------
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}

imbalance = {i: smat[i].imbalance_db(PORT_NAMES) for i in DRIVEN}
excess_loss = {i: smat[i].excess_loss_db(PORT_NAMES) for i in DRIVEN}
cmrr_I = {i: smat[i].cmrr_db("out5", "out8") for i in DRIVEN}
cmrr_Q = {i: smat[i].cmrr_db("out6", "out7") for i in DRIVEN}

phi = {j: smat[2].cross_excitation_phase_diff_deg(smat[0], port=PORT_NAMES[j]) for j in OUTPUTS}
dphi = {}
for j in [1, 2, 3]:
    raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
    dphi[j] = ((raw + 180.0) % 360.0) - 180.0

print(f"\nResults at mm_length={MM_LENGTH:.4f} µm:")
for i, port_label in zip(DRIVEN, ["1", "3"]):
    print(f"  Input port {port_label}: imbalance={imbalance[i]} dB  "
          f"excess_loss={excess_loss[i]} dB  CMRR_I={cmrr_I[i]} dBe  CMRR_Q={cmrr_Q[i]} dBe")
print(f"  Phase error Δφ_6/7/8 = {dphi[1]} / {dphi[2]} / {dphi[3]} deg")

# ---------------------------------------------------------------------------
# 5. Save.
# ---------------------------------------------------------------------------
smat[0].save(os.path.join(output_dir, "smatrix_eme_port1.pkl"))
smat[2].save(os.path.join(output_dir, "smatrix_eme_port3.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix_eme_port{{1,3}}.pkl")

np.savez(os.path.join(output_dir, "eme_results.npz"),
         wavelengths=wavelengths, mm_length=MM_LENGTH,
         imbalance_in1=imbalance[0], imbalance_in3=imbalance[2],
         excess_loss_in1=excess_loss[0], excess_loss_in3=excess_loss[2],
         cmrr_I_in1=cmrr_I[0], cmrr_Q_in1=cmrr_Q[0],
         cmrr_I_in3=cmrr_I[2], cmrr_Q_in3=cmrr_Q[2],
         dphi_6=dphi[1], dphi_7=dphi[2], dphi_8=dphi[3],
         phi_5=phi[0], phi_6=phi[1], phi_7=phi[2], phi_8=phi[3])
print(f"Results saved → {output_dir}/eme_results.npz")

# ---------------------------------------------------------------------------
# 6. Plot -- same 2x2 layout as mmi_4x4_90hybrid_hybrid.py's own results
#    plot (CMRR, phase error, excess loss, quadrature phase), so the two
#    methods' results are directly, visually comparable.
# ---------------------------------------------------------------------------
WRAPPED_TARGET_DEG = {j: ((v + 180.0) % 360.0) - 180.0 for j, v in TARGET_DEG.items()}

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

ax = axes[0, 0]
ax.plot(wavelengths, cmrr_I[0], "o-", label="CMRR_I (port 1)", color="steelblue")
ax.plot(wavelengths, cmrr_Q[0], "o-", label="CMRR_Q (port 1)", color="darkorange")
ax.plot(wavelengths, cmrr_I[2], "s--", label="CMRR_I (port 3)", color="steelblue", alpha=0.5)
ax.plot(wavelengths, cmrr_Q[2], "s--", label="CMRR_Q (port 3)", color="darkorange", alpha=0.5)
ax.axhline(20.0, color="red", linestyle=":", lw=1.0, label="cited requirement (20 dBe)")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("CMRR (dBe)")
ax.set_title("CMRR vs Wavelength")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[0, 1]
ax.plot(wavelengths, dphi[1], "o-", label="Δφ_6", color="steelblue")
ax.plot(wavelengths, dphi[2], "o-", label="Δφ_7", color="darkorange")
ax.plot(wavelengths, dphi[3], "o-", label="Δφ_8", color="seagreen")
ax.axhline(5.0, color="red", linestyle=":", lw=1.0, label="cited requirement (5°)")
ax.axhline(-5.0, color="red", linestyle=":", lw=1.0)
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Phase error Δφ_j (deg)")
ax.set_title("Phase Error vs Wavelength")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[1, 0]
ax.plot(wavelengths, excess_loss[0], "o-", label="Input port 1", color="steelblue")
ax.plot(wavelengths, excess_loss[2], "s--", label="Input port 3", color="darkorange")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Excess loss (dB)")
ax.set_title("Excess Loss vs Wavelength (Eq. 11)")
ax.legend(fontsize=7)
ax.grid(True)

ax = axes[1, 1]
_phi_labels = ["φ_5 (ref)", "φ_6", "φ_7", "φ_8"]
_phi_colors = ["gray", "steelblue", "darkorange", "seagreen"]
for j in OUTPUTS:
    ax.plot(wavelengths, phi[j], "o-", label=_phi_labels[j], color=_phi_colors[j])
    ax.axhline(WRAPPED_TARGET_DEG[j], color=_phi_colors[j], linestyle=":", lw=1.0, alpha=0.6)
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Quadrature phase φ_j (deg)")
ax.set_title("Quadrature Phase vs Wavelength\n(dotted = paper's target angle, Eq. 15)")
ax.legend(fontsize=7)
ax.grid(True)

plt.suptitle(f"4×4 MMI 90° hybrid — EME — mm_length={MM_LENGTH:.2f} µm", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "eme_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/eme_results.png")
