"""4x4 MMI 90-degree hybrid -- Step 3: EME body length (ΔL) refinement sweep.

validation_markdown/mmi_4x4_90_hybrid_validation.md Step 3: replaces the
paper's own real-coded micro-GA (A3) with a cheap analytic length sweep
around the Eq. 7 self-imaging point, since self-imaging theory is expected
to hold much better on this strongly-guided SOI platform than on the
paper's weakly-guided silica -- to be checked here, not assumed.

Structurally this is mmi_2x2_hybrid.py's own cascade
(input FDTD section -> EME body, length-swept analytically -> output,
mirror-reused from the input per A5), generalized from 2 ports to 4:

- 4 FDTD input-section solves (one per port index 0..3, input_idx=k) --
  NOT 2, even though only ports 1/3 (0-indexed 0/2) are ever driven (A4):
  ports 2/4 (index 1/3) are still needed as OUTPUT matrices via mirror
  reuse, since output_idx=j uses M_in[j].transpose() by A5's confirmed
  symmetry (paper Eq. 8, y_i = y_{i+4} -- verified directly on a real
  MMI(4,4) instance in Step 2, not just argued).
- 1 EME body length-sweep (one mode solve per wavelength, every length in
  the ΔL grid applied analytically -- zero extra solves per length, same
  as mmi_2x2_hybrid.py/mmi_1x4_hybrid.py's own body sweeps).
- The output section is NEVER separately simulated, for the same reason
  the 2x2 doesn't: A5's mirror-symmetry reuse.

N_MODES/N_WL -- a real cost surprise, not guessed. Step 3a's first mode-
count attempt (job 11706986, 16 modes) TIMED OUT at 60 min with ZERO
output. A calibration point (job 11707835, 8 modes) then measured
~416 s/mode (874.6s for 2 modes, 3324.7s for 8 -- close to linear, not the
"few times more" originally guessed), and ALL 8 requested modes came back
guided with healthy, still-widening margin (n_eff 2.962->2.833, cladding
at 1.4468) -- this cross-section supports considerably more guided modes
than any existing device's own interface basis. solve_section_length_
sweep() costs one full mode solve PER WAVELENGTH (independently), so
N_MODES=16 x N_WL=11 (this repo's usual convention, copied uncritically at
first) would have cost roughly 20+ hours for the body alone. Rescoped:

- N_MODES=8, reusing job 11707835's own measurement -- not exhaustive, a
  deliberate coarse-first-pass choice (MMI self-imaging is typically
  dominated by low-order modes), flagged for a real audit later if this
  pass looks physically sensible and needs refinement.
- N_WL=3, at the O-band edges+centre {1.26, 1.31, 1.36} um -- mirrors the
  PAPER'S OWN fitness convention exactly (Sec. IV / Eq. 12: "the set of
  wavelengths used to evaluate the fitness was Lambda = {1520, 1550,
  1580 nm}", 3 points, not 11), not an arbitrary cut. This repo's other
  MMI examples used N_WL=11 because their own cross-sections were cheap
  enough to afford it; this one is not, and 3 is a real, precedented choice
  here rather than an ad hoc concession.

This brings the body sweep to roughly 3 x 55 min if serial (or ~55 min if
array-parallelized, one wavelength per task -- this repo's established
pattern for exactly this shape of cost).

Evaluates, at every swept mm_length: imbalance_db (Eq. 10, over all 4
outputs) and excess_loss_db (Eq. 11) for each driven input via the new
SMatrix methods (Step 7/A8), and the phase-error residual (Eq. 15) via
SMatrix.cross_excitation_phase_diff_deg -- the first real exercise of
those four methods against a real (if unrefined) device, not synthetic
data.

Run:
    python examples/verification/python/mmi_4x4_90hybrid_body_sweep.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, EMESolver, FDTDSolver, HybridSolver, SMatrix, SOI_220, plot_hybrid_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "body_sweep"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters -- RESOLUTION/DECAY_TOL/N_WL are env-var
# overridable, same pattern modal_multiplexer_fdtd.py's own SMOKE run used
# (MUX_FDTD_RES): the input-section FDTD path (solve_section(), real time-
# stepping) has never been exercised on this geometry at all yet, so a cheap
# smoke pass (low res, loosened decay_tol, N_WL=1) must run and be checked
# BEFORE the real values below are trusted with real cluster time.
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI4X4_RES", 40.0))
DECAY_TOL = float(os.environ.get("MMI4X4_DECAY_TOL", 1e-11))
WL_MIN = 1.26         # O-band, matching this repo's existing O-band devices (Step 0)
WL_MAX = 1.36
N_WL = int(os.environ.get("MMI4X4_NWL", 3))   # matches the paper's own Eq. 12 fitness
                      # convention (Sec. IV): Lambda = {1520,1550,1580 nm}, 3 points --
                      # see module docstring for why this can't be 11 like the other examples

INPUT_METHOD = "fdtd"
BODY_METHOD = "eme"
OUTPUT_METHOD = "fdtd"   # never actually run -- mirror reuse (A5)

# ---------------------------------------------------------------------------
# 1. Geometry -- Step 2's settled first-pass parameters.
# ---------------------------------------------------------------------------
WAVEGUIDE_WIDTH = 0.35        # user's direct decision (Step 1)
MM_WIDTH = 6.0                # candidate, measured point (Step 1c)
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0   # = 1.5 -- reproduces paper Eq. 8 exactly (Step 2, checked)
TAPER_WIDTH = 0.8             # PLACEHOLDER -- MMI's own class default, not measured (Step 2)
TAPER_LENGTH = 1.0            # PLACEHOLDER -- MMI's own class default, not measured (Step 2)
L_PI = 108.9195               # MEASURED, job 11693043, mm_width=6.0, O-band SOI_220 (Step 1c)
MM_LENGTH_NOMINAL = 0.75 * L_PI   # Eq. 7, ΔL=0 -- Step 2's first-pass length

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH_NOMINAL, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
    use_neighbor_aperture=True,   # A6 -- needed on BOTH 4-port ends, not half_domain_aperture
)
print(f"MMI total length: {mmi.total_length:.4f} µm (mm_length={mmi.mm_length:.4f})")
print(f"Input ports at y = {mmi.input_ys} µm; outputs at y = {mmi.output_ys} µm")
assert mmi.input_ys == mmi.output_ys, "A5 mirror-symmetry violated -- check port construction"

# ---------------------------------------------------------------------------
# 2. Sub-solvers.
# ---------------------------------------------------------------------------
_fdtd_common = dict(
    res=RESOLUTION,
    sim_width=MM_WIDTH + 3.0,   # margin beyond the widest (body-facing) cross-section
    sim_height=2.0,             # Step 1b's converged sz at this platform/band
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    material_stack=SOI_220,
    decay_tol=DECAY_TOL,
)

N_MODES = 16   # was 8 -- bumped 2026-08-09. The S-parameter convergence sweep
              # (job 11752350, "tidy3d" backend) found N=8 under-converged by
              # 50-64% relative to the true value, which converges by N=16-20.
              # This body sub-solver is deliberately locked to "meep" (not
              # "tidy3d") for FDTD/EME mode-basis consistency across the
              # hybrid cascade -- a bare N=16/20 request via meep at the WIDE
              # domain body_modes.py/cross_section.py use (oxide_width=
              # MM_WIDTH+14.0) timed out at both 3h and 4h with zero output
              # (jobs 11753042/11749479). Root cause found, not just worked
              # around: THIS script's own domain is much narrower (oxide_
              # width=MM_WIDTH+3.0, matching the FDTD sections it interfaces
              # with) -- at THAT actual production width, N=16 via meep
              # completed in 77 min (job 11767748), and N=8 in 23 min (job
              # 11767747), both faster than the historical wide-domain timing.
              # A parallel tidy3d screen confirmed mode content (n_eff) is
              # bit-for-bit IDENTICAL from oxide_width=MM_WIDTH+3.0 all the
              # way to +20.0 and oxide_height=2.0 through 6.0 -- the domain
              # was never too small for accuracy, it was the DIAGNOSTIC
              # scripts' own much-wider domain that made meep's eigensolve
              # dramatically slower for higher-order modes (more iterations
              # per mode, MPB log showed 160 iters for mode 16 at the wide
              # domain). The wall was a property of the wrong test
              # configuration, not of this script's own real domain.


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
# 3. HybridSolver -- for the geometry/section-map plot only. Step 5 below
#    hand-rolls the einsum cascade (4 excitations x 4 outputs, mixed with a
#    length-swept body), the same reason mmi_2x2_hybrid.py/directional_
#    coupler_hybrid.py do: HybridSolver's one-call contract can't express it.
# ---------------------------------------------------------------------------
half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
body_stub = 0.25   # µm of MMI body folded into each FDTD section, past the taper/body seam

boundaries = [
    -half_total,
    -half_body + body_stub,
    +half_body - body_stub,
    +half_total,
]
interface_modes = [N_MODES, N_MODES]

hybrid = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries, interface_modes=interface_modes,
)

plot_hybrid_geometry(
    mmi, solver=hybrid, in_modes=1, out_modes=1,
    title="4×4 MMI 90° hybrid — Body ΔL Sweep — Sections",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True, sub_region_show_eme_detail=True,
    # One sub-region plot per simulation actually run: 4 input-arm FDTD
    # solves (A5), body via EME. mirror_pairs is {section: partner} in
    # SECTION-index space (this HybridSolver has 3 sections: 0=input,
    # 1=body, 2=output) -- section 2 always mirrors section 0, regardless
    # of which physical input_idx k that config excited. (Not {2: k}: k
    # ranges 0..3, ports, but section indices only go 0..2.)
    sub_region_configs=[
        (f"in{k}", {0: {"input_idx": k}}, {2: 0}) for k in range(4)
    ],
)
print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")

# ---------------------------------------------------------------------------
# 4. ΔL sweep grid. Paper's own Table 1 search bound was ΔL in [-10,10] µm
#    against L≈459 µm (i.e. its μGA's OWN local search window around the
#    analytic self-imaging point, not just a fine-tune afterthought) --
#    reused directly here rather than a scaled-down guess, since the length
#    sweep is analytically free once M_body is computed (one mode solve per
#    wavelength, every length applied as a phase multiply -- zero extra
#    solves, same as mmi_2x2_hybrid.py's own 201-point sweep).
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
delta_L = np.linspace(-10.0, 10.0, 81)                # µm, 0.25 µm steps
mm_lengths = MM_LENGTH_NOMINAL + delta_L
body_lengths = mm_lengths - 2 * body_stub

print(f"\nBody sweep: ΔL ∈ [{delta_L[0]:+.2f}, {delta_L[-1]:+.2f}] µm around "
      f"mm_length={MM_LENGTH_NOMINAL:.4f} µm ({len(mm_lengths)} pts), "
      f"{len(wavelengths)} wavelengths {WL_MIN}-{WL_MAX} µm, N_MODES={N_MODES}\n")

# ---------------------------------------------------------------------------
# 5. Cascade -- A5's mirror-symmetry reuse, generalized from 2 to 4 ports.
#    4 FDTD input-section solves (once each, reused for both driven inputs
#    AND all 4 outputs), 1 EME body length-sweep, NO output-section solves.
# ---------------------------------------------------------------------------
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

M_out = [Mk.transpose(0, 2, 1) for Mk in M_in]   # body modes -> physical output port j (A5 mirror reuse)

M_body = body_solver.solve_section_length_sweep(
    mmi, wavelengths, boundaries[1], boundaries[2], body_lengths,
    n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
)   # [n_L, n_wl, N_MODES, N_MODES]
print(f"  body length-sweep solved: {M_body.shape}")

DRIVEN = [0, 2]     # paper's ports 1 and 3 (A4) -- 0-indexed
OUTPUTS = [0, 1, 2, 3]   # paper's ports 5, 6, 7, 8 -- 0-indexed

# coeffs[(i, j)]: [n_L, n_wl] complex amplitude, port j <- driven port i
coeffs = {}
for i in DRIVEN:
    tmp = np.einsum("Lwab,wbk->Lwak", M_body, M_in[i])       # [n_L, n_wl, N_MODES, 1]
    for j in OUTPUTS:
        c = np.einsum("wab,Lwbk->Lwak", M_out[j], tmp)       # [n_L, n_wl, 1, 1]
        coeffs[(i, j)] = c[:, :, 0, 0]

# ---------------------------------------------------------------------------
# 6. Evaluate every length: imbalance/excess-loss (Eq. 10/11, per driven
#    input, over all 4 outputs) and the Eq. 15 phase-error residual, via
#    per-length SMatrix instances -- the first real exercise of Step 7's new
#    methods, not synthetic data.
# ---------------------------------------------------------------------------
PORT_NAMES = ["out5", "out6", "out7", "out8"]
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}   # paper's phi~_5..8 (Eq. 15)

n_L, n_wl = len(mm_lengths), len(wavelengths)
imbalance_grid = {i: np.zeros((n_L, n_wl)) for i in DRIVEN}
excess_loss_grid = {i: np.zeros((n_L, n_wl)) for i in DRIVEN}
max_abs_dphi_grid = np.zeros((n_L, n_wl))   # max_j |Δφ_j| over j=1,2,3 (Eq. 15, φ~_5 is the ref)

for Li in range(n_L):
    smat = {i: SMatrix(wavelengths=wavelengths) for i in DRIVEN}
    for i in DRIVEN:
        for j in OUTPUTS:
            smat[i].set(coeffs[(i, j)][Li], port=PORT_NAMES[j])
        imbalance_grid[i][Li] = smat[i].imbalance_db(PORT_NAMES)
        excess_loss_grid[i][Li] = smat[i].excess_loss_db(PORT_NAMES)

    # Eq. 15: phi_j = angle(S_j3) - angle(S_j1).
    #
    # MEASURED 2026-08-06, against job 11708624's real (cached) results: calling
    # this as smat[0].cross_excitation_phase_diff_deg(smat[2], ...) -- self=port1,
    # other=port3, matching the method's own documented convention and the
    # paper's stated Eq. 15 literally -- gave max|Delta phi| pinned near 172-180
    # degrees across the ENTIRE 81-point ΔL sweep (essentially the worst
    # possible value), while max_imbalance varied sensibly (4.9-39.6 dB) over
    # the same sweep. A magnitude-only metric behaving sensibly while a
    # phase-only metric is stuck at ~180 degrees, invariant to a parameter
    # (body length) that should strongly affect real relative phase, is the
    # signature of a GLOBAL SIGN CONVENTION mismatch, not a real design
    # failure or a solve_section() cross-call phase-reference bug (that
    # specific risk, flagged in A8, was checked directly against real cached
    # data first: the input_idx=1 vs input_idx=2 mirror-port pair shows
    # EXACTLY the expected even/odd-mode symmetry, 0.00/180.00 degrees to
    # 2-3 decimal places -- solve_section() calls ARE phase-consistent with
    # each other).
    #
    # Confirmed: swapping self/other here (equivalent to phi_j = angle(S_j1)
    # - angle(S_j3), the negative of the paper's literal Eq. 15) brings
    # max|Delta phi| from ~172-180 deg down to a clean, smooth 7.3-180 deg
    # trough over the same sweep (minimum 7.345 deg at ΔL=0.00, close to but
    # not quite meeting the cited 5 deg requirement) -- meep/meow's own
    # S-parameter phase convention is apparently conjugated (a global sign
    # flip, e.g. an exp(-i*omega*t) vs exp(+i*omega*t) time-harmonic
    # convention difference) relative to the paper's implicit one. This
    # flips ALL phase differences uniformly, which is exactly what a global
    # convention mismatch does and a per-device design error would not.
    # validation_markdown/mmi_4x4_90_hybrid_validation.md finding A11 has the
    # full diagnostic trail (mirror-pair check, raw-vs-target comparison).
    phi = {j: smat[2].cross_excitation_phase_diff_deg(smat[0], port=PORT_NAMES[j]) for j in OUTPUTS}
    dphi = np.zeros((3, n_wl))
    for row, j in enumerate([1, 2, 3]):
        raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
        dphi[row] = ((raw + 180.0) % 360.0) - 180.0   # mod-360-aware residual (SMatrix's own caveat)
    max_abs_dphi_grid[Li] = np.max(np.abs(dphi), axis=0)

# ---------------------------------------------------------------------------
# 7. Score every length and report the best.
# ---------------------------------------------------------------------------
max_imbalance = np.maximum(imbalance_grid[0].max(axis=1), imbalance_grid[2].max(axis=1))   # [n_L]
max_dphi = max_abs_dphi_grid.max(axis=1)   # [n_L]

best_L_idx = int(np.argmin(max_dphi + max_imbalance))   # simple combined score, not independently tuned
print(f"\nBest ΔL = {delta_L[best_L_idx]:+.2f} µm (mm_length={mm_lengths[best_L_idx]:.4f} µm): "
      f"max|Δφ|={max_dphi[best_L_idx]:.3f}°, max imbalance={max_imbalance[best_L_idx]:.4f} dB")
print(f"System requirements cited by the paper (from its own refs [4,5]): "
      f"CMRR ≥ 20 dBe, |Δφ_j| ≤ 5° -- NOT the paper's own silica numbers (A9)")

np.savez(os.path.join(output_dir, "body_sweep.npz"),
         delta_L=delta_L, mm_lengths=mm_lengths, wavelengths=wavelengths,
         max_imbalance=max_imbalance, max_dphi=max_dphi,
         imbalance_in1=imbalance_grid[0], imbalance_in3=imbalance_grid[2],
         excess_loss_in1=excess_loss_grid[0], excess_loss_in3=excess_loss_grid[2],
         max_abs_dphi_grid=max_abs_dphi_grid)
print(f"Sweep grids saved → {output_dir}/body_sweep.npz")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
ax = axes[0]
ax.plot(delta_L, max_imbalance, color="steelblue")
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8)
ax.set_xlabel("ΔL (µm)")
ax.set_ylabel("max imbalance over λ, both driven inputs (dB)")
ax.set_title("Imbalance vs ΔL")
ax.grid(True)

ax = axes[1]
ax.plot(delta_L, max_dphi, color="darkorange")
ax.axhline(5.0, color="red", linestyle="--", lw=0.8, label="cited requirement (5°)")
ax.set_xlabel("ΔL (µm)")
ax.set_ylabel("max |Δφ_j| over λ, j=6,7,8 (deg)")
ax.set_title("Phase error vs ΔL")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle("4×4 MMI 90° hybrid — Step 3 body ΔL sweep", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "body_sweep.png"), dpi=300)
plt.close("all")
print(f"Sweep plot saved → {output_dir}/body_sweep.png")
