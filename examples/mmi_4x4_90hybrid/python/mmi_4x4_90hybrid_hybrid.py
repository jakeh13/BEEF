"""Tutorial: 4x4 MMI 90-degree optical hybrid -- FDTD/EME hybrid simulation
(BEEF method).

validation_markdown/mmi_4x4_90_hybrid_validation.md Step 6c. Reproduces the
device topology of Yasui, Sugisaka & Hirayama, "Structural Optimization of
an Optical 90 Degree Hybrid Based on a Weakly-guided 4x4 Multimode
Interference Coupler...", ACES Journal 36(5), 2021 -- on 220 nm SOI at
O-band instead of the paper's own weakly-guided silica/C-band, so the
paper's own numbers are a design-methodology reference (equations, port
topology, cited system requirements), NOT a numerical target (A9).

Structurally this is mmi_2x2_hybrid.py's own cascade (input FDTD section ->
EME body, length-swept analytically -> output, mirror-reused from the
input), generalized from 2 ports to 4 -- see that script's own module
docstring for the general pattern. The 4-port generalization needed here:

- 4 FDTD input-section solves (one per port index 0..3, input_idx=k), NOT
  2, even though only ports 1/3 (0-indexed 0/2) are ever driven (paper
  Sec. III, finding A4): ports 2/4 (index 1/3) are still needed as OUTPUT
  matrices via mirror reuse, since output_idx=j uses M_in[j].transpose()
  by the device's mirror symmetry (paper Eq. 8, y_i = y_{i+4} -- verified
  directly on a real MMI(4,4) instance, finding A5).
- use_neighbor_aperture=True on the geometry (not half_domain_aperture on
  the FDTD solver) -- needed on BOTH 4-port ends, unlike mmi_1x4_hybrid.py
  which only needed it on its output end (finding A6).
- The output section is NEVER separately simulated, for the same mirror-
  symmetry reason mmi_2x2_hybrid.py's own output section isn't.

N_MODES=8 and the geometry parameters below are Step 1-3's settled first
pass, NOT a converged final design:

- `waveguide_width=0.35` -- the user's direct design choice (Step 1),
  not derived from a single-mode search (no width in [0.10,0.70] um is
  single-mode at this platform/band -- not even this repo's own shipped
  2x2 MMI access waveguide meets that bar; see finding A10).
- `mm_width=6.0`, `waveguide_offset=mm_width/4=1.5` -- the offset formula
  is not arbitrary: it is EXACTLY what reproduces the paper's own Eq. 8
  port-position formula through this codebase's plain default MMI port
  constructor (Step 2, checked directly, not argued).
- `mm_length=80.1896` (ΔL=-1.50 um around the paper's Eq. 7 self-imaging
  point) -- the best point found in a ±10 um sweep (job 11708624), not a
  converged optimum. At this point: max|Δphi_j|=8.28 deg, CMRR ranges
  ~33 dBe at O-band centre down to ~13-19 dBe at the band edges --
  MEETS the cited CMRR>=20dBe requirement at centre, falls short at the
  edges; phase error (8.28 deg) is close to but does not yet meet the
  cited |Δphi|<=5 deg bar. Carried forward as a documented, known
  limitation (validation doc finding A11) -- refining N_MODES beyond this
  coarse first pass, or the taper_width/taper_length placeholders below
  (MMI's own class defaults, never independently measured), are the
  likely next levers, not attempted here at the user's direction.
- taper_width=0.8, taper_length=1.0 -- PLACEHOLDERS, MMI's own class
  defaults. Geometrically checked (0.7 um taper-overlap clearance, Step 2)
  but never independently optimized for this device.

A11 (validation doc): meep/meow's own S-parameter phase convention is
CONJUGATED relative to this paper's implicit Eq. 15 convention (phi_j =
angle(S_j3) - angle(S_j1) as literally written needs its sign flipped here
to match measured results against the paper's own target angles -- a
diagnosed, verified global sign convention mismatch, not a device design
error or a solve_section() cross-call phase-reference bug). The fix is
applied below, at the same call site, with the same reasoning recorded
in-line.

Run:
    python /path/to/examples/mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_hybrid.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, EMESolver, FDTDSolver, HybridSolver, SMatrix, SOI_220, plot_hybrid_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "mmi_4x4_90hybrid_hybrid"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN = 1.26     # O-band, matching this repo's existing O-band devices
WL_MAX = 1.36
N_WL = 3          # matches the paper's own Eq. 12 fitness convention (Sec. IV):
                  # Lambda = {1520,1550,1580 nm}, not this repo's usual 11 --
                  # this cross-section is too expensive to afford that (Step 3
                  # cost measurement: ~416 s/body-mode, 4 real FDTD input-
                  # section solves per full sweep)

INPUT_METHOD = "fdtd"
BODY_METHOD = "eme"
OUTPUT_METHOD = "fdtd"   # never actually run -- mirror reuse

# ---------------------------------------------------------------------------
# 1. Geometry -- Steps 1-3's settled first-pass parameters.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 6.0
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0   # reproduces paper Eq. 8 exactly (Step 2, checked)
TAPER_WIDTH = 0.8                   # PLACEHOLDER -- MMI's own class default
TAPER_LENGTH = 1.0                  # PLACEHOLDER -- MMI's own class default
L_PI = 108.9195                     # MEASURED, job 11693043, mm_width=6.0, O-band SOI_220

# MM_LENGTH: Eq. 7's L=(3/4)*L_pi point (81.6896) refined by the Step 3 ΔL
# sweep (job 11708624 + local recompute after the A11 sign fix) -- best
# combined point over a ±10 um, 0.25 um-step grid. Asserted against that
# grid rather than a bare float literal, per mmi_1x4_hybrid.py's own
# grid-snap lesson (a "pinned" length silently 0.0208 um off its own sweep
# grid, shipped for weeks).
_DELTA_L_GRID = np.linspace(-10.0, 10.0, 81)
DELTA_L = -2.75   # was -1.50 -- re-derived 2026-08-09 at N_MODES=16 (job
                  # 11768847), see the N_MODES comment below for why -1.50
                  # (found at the under-converged N_MODES=8) no longer applies
_nearest = _DELTA_L_GRID[np.argmin(np.abs(_DELTA_L_GRID - DELTA_L))]
assert _nearest == DELTA_L, (
    f"DELTA_L={DELTA_L} is not an exact point on the Step 3 sweep grid "
    f"(nearest: {_nearest}) -- do not silently snap."
)
MM_LENGTH = 0.75 * L_PI + DELTA_L   # = 80.1896

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
    use_neighbor_aperture=True,
)
print(f"MMI total length: {mmi.total_length:.4f} µm (mm_length={mmi.mm_length:.4f})")
print(f"Input ports at y = {mmi.input_ys} µm; outputs at y = {mmi.output_ys} µm")
assert mmi.input_ys == mmi.output_ys, "mirror-symmetry violated -- check port construction"

# ---------------------------------------------------------------------------
# 2. Sub-solvers.
# ---------------------------------------------------------------------------
_fdtd_common = dict(
    res=RESOLUTION,
    sim_width=MM_WIDTH + 3.0,
    sim_height=2.0,             # Step 1b's converged sz at this platform/band
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    material_stack=SOI_220,
    decay_tol=1e-11,
)

N_MODES = 16   # was 8 -- bumped 2026-08-09. The S-parameter convergence
              # sweep (job 11752350, "tidy3d" backend) found N=8 under-
              # converged by 50-64% relative to the true value, which
              # converges by N=16-20. Re-running this script's own
              # meep-backed body sub-solver at N=16, at THIS script's actual
              # production domain (oxide_width=MM_WIDTH+3.0), was confirmed
              # affordable (77 min, job 11767748) -- an earlier attempt at
              # a much WIDER diagnostic-only domain (oxide_width=
              # MM_WIDTH+14.0, body_modes.py's own convention) timed out at
              # both 3h and 4h with zero output; that wall was a property of
              # the wrong domain, not of N=16 itself (a parallel tidy3d
              # screen confirmed mode content is identical across every
              # domain size tested, production width through +20 um).
              # See module docstring.


def make_section_solver(method: str):
    if method == "fdtd":
        return FDTDSolver(**_fdtd_common)
    if method == "eme":
        return EMESolver(
            num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
            oxide_width=_fdtd_common["sim_width"], oxide_height=_fdtd_common["sim_height"],
            res=RESOLUTION, material_stack=SOI_220,
            mode_solver="meep",   # matches FDTDSolver's own eigenmode basis in this hybrid chain
        )
    raise ValueError(f"unknown method {method!r}, expected 'fdtd' or 'eme'")


input_solver = make_section_solver(INPUT_METHOD)
body_solver = make_section_solver(BODY_METHOD)
output_solver = make_section_solver(OUTPUT_METHOD)   # never run -- labels the section map only

# ---------------------------------------------------------------------------
# 3. HybridSolver -- for the geometry/section-map plot only. Section 5 below
#    hand-rolls the einsum cascade (4 excitations x 4 outputs), the same
#    reason mmi_2x2_hybrid.py/directional_coupler_hybrid.py did: this device's
#    N-port cascade does not fit HybridSolver's one-call contract.
#    NOTE 2026-08-24: `beef.MirrorSpec` now lets a mirror-reused matrix be
#    expressed through the supported API, so the old blanket wording here
#    ("can't express a mirror-reused N-port cascade") was imprecise. Checked
#    against THIS device rather than assumed, and the conclusion is unchanged
#    for a reason worth recording:
#      * the ±1 restriction is NOT the blocker here. MirrorSpec requires
#        strictly ±1 mode_signs and raises otherwise, which was flagged as
#        likely fatal for a 90-degree hybrid -- but this device's mirror
#        relations were measured on 2026-08-13 (circular statistics, all
#        16/16 body modes, both pairs) at R≈1.0000 with angular deviation
#        from an exact 0/180 sign ≤0.15 deg. The 90-degree relationship is
#        between OUTPUT ports for a given input, not between mirror-partner
#        INPUT excitations, so the parity really is real-valued.
#      * the blocker is the N-PORT shape. `HybridSolver.solve()` returns ONE
#        SMatrix under ONE `port` label -- a single input->output cascade.
#        This device needs EIGHT port pairs (2 driven inputs x 4 outputs)
#        built from 4 input matrices, a body, and 4 output matrices.
#        MirrorSpec removes the mirroring obstacle, not the N-port one.
#    So: still hand-rolled, but now for a narrower and better-understood
#    reason. Revisit if HybridSolver ever grows multi-excitation output.
# ---------------------------------------------------------------------------
half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
body_stub = 0.25

boundaries = [-half_total, -half_body + body_stub, +half_body - body_stub, +half_total]
interface_modes = [N_MODES, N_MODES]

hybrid = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries, interface_modes=interface_modes,
)

plot_hybrid_geometry(
    mmi, solver=hybrid, in_modes=1, out_modes=1,
    title="4×4 MMI 90° hybrid — FDTD/EME Hybrid — Sections",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True, sub_region_show_eme_detail=True,
    sub_region_configs=[
        (f"in{k}", {0: {"input_idx": k}}, {2: 0}) for k in range(4)
    ],
)
print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")

# ---------------------------------------------------------------------------
# 4. Cascade -- mirror-symmetry reuse, generalized from 2 to 4 ports.
#    4 FDTD input-section solves (reused for both driven inputs AND all 4
#    outputs), 1 EME body solve, NO output-section solves.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

M_in = []   # M_in[k]: [n_wl, N_MODES, 1] -- physical input port k -> body modes
for k in range(4):
    Mk = input_solver.solve_section(
        mmi, wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=input_solver.waveguide_length, stub_right=0.0,
        input_idx=k, save_dir=save_dir,
    )
    M_in.append(Mk)
    print(f"  input_idx={k} solved: {Mk.shape}")

M_out = [Mk.transpose(0, 2, 1) for Mk in M_in]   # body modes -> physical output port j (mirror reuse)

# Reuses Step 3's EXACT 81-point length grid, not a single-length array --
# solve_section_length_sweep()'s own cache key includes `lengths` verbatim
# (its own docstring: "lengths is included in the config"), so requesting
# a different-shaped array here would miss job 11708624's cache entirely
# and trigger a full new ~2-3h mode solve for no reason. Slicing out the
# one length this script actually reports costs nothing extra.
_body_lengths_grid = (0.75 * L_PI + _DELTA_L_GRID) - 2 * body_stub
_Li = int(np.argmin(np.abs(_DELTA_L_GRID - DELTA_L)))
M_body_grid = body_solver.solve_section_length_sweep(
    mmi, wavelengths, boundaries[1], boundaries[2], _body_lengths_grid,
    n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
)   # [81, n_wl, N_MODES, N_MODES]
M_body = M_body_grid[_Li:_Li + 1]   # [1, n_wl, N_MODES, N_MODES] -- the settled design point
print(f"  body solved (grid reused from Step 3, sliced at ΔL={DELTA_L:+.2f}): {M_body.shape}")

DRIVEN = [0, 2]     # paper's ports 1 and 3 -- 0-indexed
OUTPUTS = [0, 1, 2, 3]   # paper's ports 5, 6, 7, 8 -- 0-indexed
PORT_NAMES = ["out5", "out6", "out7", "out8"]

coeffs = {}
for i in DRIVEN:
    tmp = np.einsum("Lwab,wbk->Lwak", M_body, M_in[i])
    for j in OUTPUTS:
        c = np.einsum("wab,Lwbk->Lwak", M_out[j], tmp)
        coeffs[(i, j)] = c[0, :, 0, 0]   # drop the length axis -- one length only

# ---------------------------------------------------------------------------
# 5. Metrics (paper Eqs. 10-15) via the new SMatrix methods, plus this
#    device's own headline plot: CMRR / phase-error / excess-loss vs λ.
# ---------------------------------------------------------------------------
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}   # paper's phi~_5..8 (Eq. 15)

smat = {i: SMatrix(wavelengths=wavelengths) for i in DRIVEN}
for i in DRIVEN:
    for j in OUTPUTS:
        smat[i].set(coeffs[(i, j)], port=PORT_NAMES[j])
    smat[i].check_passive()

imbalance = {i: smat[i].imbalance_db(PORT_NAMES) for i in DRIVEN}
excess_loss = {i: smat[i].excess_loss_db(PORT_NAMES) for i in DRIVEN}
cmrr_I = {i: smat[i].cmrr_db("out5", "out8") for i in DRIVEN}
cmrr_Q = {i: smat[i].cmrr_db("out6", "out7") for i in DRIVEN}

# A11: meep/meow's S-parameter phase convention is conjugated relative to
# the paper's implicit one -- self/other swapped here (phi_j = angle(S_j1)
# - angle(S_j3), the negative of the paper's literal Eq. 15) to match.
# Verified against job 11708624's cached results: brings max|Δphi| from a
# constant ~172-180 deg (a global sign artifact, ruled out as a cross-call
# phase-reference bug via a mirror-port check first) down to a clean,
# physically sensible 7.3-8.5 deg minimum -- see validation doc finding A11
# for the full diagnostic trail.
phi = {j: smat[2].cross_excitation_phase_diff_deg(smat[0], port=PORT_NAMES[j]) for j in OUTPUTS}
dphi = {}
for j in [1, 2, 3]:
    raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
    dphi[j] = ((raw + 180.0) % 360.0) - 180.0

print(f"\nResults at mm_length={MM_LENGTH:.4f} µm (ΔL={DELTA_L:+.2f}):")
for i, port_label in zip(DRIVEN, ["1", "3"]):
    print(f"  Input port {port_label}: imbalance={imbalance[i]} dB  "
          f"excess_loss={excess_loss[i]} dB  CMRR_I={cmrr_I[i]} dBe  CMRR_Q={cmrr_Q[i]} dBe")
print(f"  Phase error Δφ_6/7/8 = {dphi[1]} / {dphi[2]} / {dphi[3]} deg")
print(f"  System requirements cited by the paper (from its own refs [4,5]): "
      f"CMRR ≥ 20 dBe, |Δφ_j| ≤ 5° -- NOT the paper's own silica numbers (see validation doc A9)")

# ---------------------------------------------------------------------------
# 6. Save and plot.
# ---------------------------------------------------------------------------
smat[0].save(os.path.join(output_dir, "smatrix_hybrid_port1.pkl"))
smat[2].save(os.path.join(output_dir, "smatrix_hybrid_port3.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix_hybrid_port{{1,3}}.pkl")

np.savez(os.path.join(output_dir, "hybrid_results.npz"),
         wavelengths=wavelengths, mm_length=MM_LENGTH, delta_L=DELTA_L,
         imbalance_in1=imbalance[0], imbalance_in3=imbalance[2],
         excess_loss_in1=excess_loss[0], excess_loss_in3=excess_loss[2],
         cmrr_I_in1=cmrr_I[0], cmrr_Q_in1=cmrr_Q[0],
         cmrr_I_in3=cmrr_I[2], cmrr_Q_in3=cmrr_Q[2],
         dphi_6=dphi[1], dphi_7=dphi[2], dphi_8=dphi[3],
         phi_5=phi[0], phi_6=phi[1], phi_7=phi[2], phi_8=phi[3])
print(f"Results grids saved → {output_dir}/hybrid_results.npz")

# Paper's own target angles, wrapped to (-180, 180] to match
# cross_excitation_phase_diff_deg()'s own output range -- 225/135 deg
# outside that range would otherwise sit off-plot next to phi_j's own
# wrapped curve (SMatrix.cross_excitation_phase_diff_deg()'s own
# docstring flags exactly this mod-360 caveat).
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

# Quadrature phase: the ACTUAL measured phi_j = angle(S_j3) - angle(S_j1)
# (A11's sign-corrected direction), not just its residual from target
# (that's the Δφ panel above) -- shows the raw 90°-spaced I/Q structure
# directly, each port against its own paper-cited target angle (dotted).
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

plt.suptitle(f"4×4 MMI 90° hybrid — mm_length={MM_LENGTH:.2f} µm (ΔL={DELTA_L:+.2f} µm)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/hybrid_results.png")
