"""Tutorial: 4x4 MMI 90-degree optical hybrid -- FDTD/EME hybrid simulation
(BEEF method).

Reproduces the device topology of Yasui, Sugisaka & Hirayama, "Structural Optimization of
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

N_MODES=16 and the geometry below are the optimised design point of a
~45-evaluation geometry campaign (mm_width, y_gap, taper_width,
taper_length, mm_length):

- `waveguide_width=0.35` -- a direct design choice, not derived from a
  single-mode search (no width in [0.10, 0.70] um is single-mode at this
  platform/band).
- `mm_width=8.5`, `waveguide_offset=mm_width/4` -- the offset formula is
  exactly what reproduces the paper's Eq. 8 port positions through the
  default MMI port constructor.
- `taper_width=1.30`, `taper_length=4.0` -- measured optima at this width.
- `mm_length=159.30` (DELTA_L=-4.65 um around the paper's Eq. 7
  self-imaging point), selected from a shared 20-400 um / 0.05 um grid.

Measured at this point on the script's 11-point 1.29-1.33 um grid:
max|dphi_j| = 6.22 deg at 1.322 um, min CMRR = 18.48 dBe (at 1.330 um),
excess loss <= 2.71 dB, imbalance <= 1.97 dB. Against the cited
requirements (CMRR >= 20 dBe, |dphi_j| <= 5 deg) this point meets NEITHER
across the whole band. On the 3-point subgrid {1.29, 1.31, 1.33} the
campaign optimised against, phase reads 3.93 deg (met): the violation sits
in 1.318-1.326 um, between those samples.

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
# Output directory carries N_WL when it is not the default 3. Fourth instance of
# the same hazard the _tl / _n / _wl filename tags exist for (see wg_eval.py):
# a denser wavelength sweep is a DIFFERENT measurement of the same geometry, and
# these scripts write fixed filenames, so without this a denser run would
# silently destroy the 3-point results that were verified on 2026-09-06 to
# reproduce the optimiser exactly. Default stays untagged so every existing
# path and reference keeps working.
# 2026-09-15: DEFAULT RAISED 3 -> 11 over the same 1.29-1.33 um window.
# The tag below is re-keyed to this constant, so the DEFAULT is untagged
# (production) and any other count is tagged -- the previous 3-point
# production output was snapshotted to
# mmi_4x4_90hybrid_hybrid_nwl3_PRE_20260915/ first, since it is referenced by
# prior work and the optimiser's own fitness basis.
#
# ⚠ TWO REAL CONSEQUENCES, both measured rather than assumed:
#   1. COST. solve_section_length_sweep() does one eigenmode solve per
#      wavelength, and this body is ~416 s/mode (3324.7 s for 8 modes vs
#      874.6 s for 2, job 11693043 -- validation doc line 1508). At
#      N_MODES=16 that is ~1.85 h PER WAVELENGTH: ~5.5 h at 3, ~20 h at 11.
#      MPI ranks do NOT help -- MPB is rank-redundant (measured), so that is
#      serial wall time. Only the 4 FDTD input sections scale with ranks.
#   2. FITNESS BASIS. The old default of 3 deliberately matched the paper's
#      Eq. 12 convention (Sec. IV, |Lambda| = 3), which is what the BO
#      campaign used to select mm_length=159.30. Metrics computed at 11
#      wavelengths are a DIFFERENT, denser measurement of the same geometry --
#      better band coverage, but no longer the optimiser's own fitness basis.
#      Set MMI4X4_HYB_NWL=3 to reproduce that basis exactly.
_NWL_DEFAULT = 11
_NWL_ENV = int(os.environ.get("MMI4X4_HYB_NWL", _NWL_DEFAULT))
_MM_LENGTH_ENV = os.environ.get("MMI4X4_MM_LENGTH")
_length_tag = "" if _MM_LENGTH_ENV is None else f"_L{float(_MM_LENGTH_ENV):.2f}"
# MMI4X4_DRIVEN_PORTS: a comma list of 1-based paper
# port numbers, default "1,3" unchanged. Parsed early (not just at the
# metrics section) so the output dir can be tagged by it too. Does NOT
# change which of the 4 M_in solves run (all 4 are always needed, mirror
# reuse for M_out) -- only which one(s) get reported/saved.
_PORT_TO_INPUT_IDX = {1: 0, 2: 1, 3: 2, 4: 3}
_driven_ports = [int(p) for p in
                  os.environ.get("MMI4X4_DRIVEN_PORTS", "1,3").split(",")]
_ports_tag = "" if set(_driven_ports) == {1, 3} else (
    "_p" + "".join(str(p) for p in _driven_ports))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid",
    "mmi_4x4_90hybrid_hybrid" + ("" if _NWL_ENV == _NWL_DEFAULT else f"_nwl{_NWL_ENV}") + _length_tag + _ports_tag))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
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
N_WL = _NWL_ENV   # env-overridable via MMI4X4_HYB_NWL. Default was 3, which
                  # matches the paper's own Eq. 12 fitness convention (Sec. IV):
                  # Lambda = {1520,1550,1580 nm}, not this repo's usual 11 --
                  # this cross-section is too expensive to afford that (Step 3
                  # cost measurement: ~416 s/body-mode, 4 real FDTD input-
                  # section solves per full sweep)

INPUT_METHOD = "fdtd"
BODY_METHOD = "eme"
OUTPUT_METHOD = "fdtd"   # never actually run -- mirror reuse

# ---------------------------------------------------------------------------
# 1. Geometry -- the optimised design point (see the module docstring).
# ---------------------------------------------------------------------------
# Geometry -- the OPTIMISED design point, not Step 2's first pass.
# Result of the (mm_width, y_gap, taper_width, taper_length, mm_length)
# campaign of ~45 hybrid-path evaluations (rounds 8-10, 2026-08-28..09-05). Measured at
# this point on the 11-point 1.29-1.33 um grid: max|dphi| = 6.223 deg at 1.322 um
# (spec <=5, NOT met), min CMRR = 18.482 dBe (spec >=20, 1.52 short), excess loss
# 2.706 dB, imbalance 1.972 dB -- the point meets neither spec. The campaign itself
# scored a 3-point {1.29, 1.31, 1.33} um subsample, on which phase reads 3.933 deg;
# the violation sits between those samples. Superseded first pass, for the record: mm_width=6.0,
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
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0   # reproduces paper Eq. 8 exactly (Step 2, checked)
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
#
# MMI4X4_MM_LENGTH overrides the production 159.30, e.g. for a short body
# (<20 um) that a full-device FDTD comparison can afford. Not claiming a new
# device optimum; self-imaging is not being targeted at these lengths.
_MM_LENGTH_DEFAULT = 159.30
MM_LENGTH = float(os.environ.get("MMI4X4_MM_LENGTH", _MM_LENGTH_DEFAULT))
_LENGTH_STEP = 0.05
assert abs(round(MM_LENGTH / _LENGTH_STEP) * _LENGTH_STEP - MM_LENGTH) < 1e-9, (
    f"MM_LENGTH={MM_LENGTH} is not an exact multiple of {_LENGTH_STEP} um -- "
    f"avoid silently snapping to an unintended grid point."
)
DELTA_L = MM_LENGTH - 0.75 * L_PI   # = -4.6465 at the default; informational
                                     # only, and meaningless (huge negative) at
                                     # an engine-agreement short length.

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
    decay_tol=1e-4,   # was 1e-11 -- MEASURED UNREACHABLE (see
               # the residual floors then diverges). 1e-4 with energy
               # decay was measured sufficient for this device at
               # 0.015 dB, and is what every campaign point used.
    decay_method="energy",
    energy_decay_dt=20.0,
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

# Reuses the campaign's EXACT length grid, not a single-length array --
# solve_section_length_sweep()'s own cache key includes `lengths` verbatim
# (its own docstring: "lengths is included in the config"), so a
# different-shaped array here misses the cache and triggers a full new
# ~70 min body mode solve for no reason. Slicing out the one length this
# script reports costs nothing extra.
#
# This is the shared ABSOLUTE grid from mmi_4x4_90hybrid_ygap_grid.py
# (MM_LENGTH_GRID, 20-400 um @ 0.05 um), which is what every campaign point
# and the Phase A body solve for mm_width=8.5 (job 12780661_0) used --
# reproduced inline rather than imported so this published example does not
# reach into examples/verification/. It replaced Step 3's +-10 um / 0.25 um
# DELTA_L window, whose cache entries are all at the superseded mm_width=6.0
# geometry and would miss regardless. body_stub (0.25) matches wg_eval's
# boundaries_for() exactly, which is what makes the hit possible.
_grid_lo = 20.0 if MM_LENGTH >= 20.0 else MM_LENGTH   # extend down only for an
                                                        # engine-agreement short-length
                                                        # override; default range (and
                                                        # every existing cache entry) is
                                                        # untouched at MM_LENGTH>=20
_LENGTH_GRID = np.round(np.arange(_grid_lo, 400.0 + 1e-9, 0.05), 2)
assert np.isclose(_LENGTH_GRID, MM_LENGTH).any(), (
    f"MM_LENGTH={MM_LENGTH} is not on the shared MM_LENGTH_GRID"
)
_body_lengths_grid = _LENGTH_GRID - 2 * body_stub
_Li = int(np.argmin(np.abs(_LENGTH_GRID - MM_LENGTH)))
M_body_grid = body_solver.solve_section_length_sweep(
    mmi, wavelengths, boundaries[1], boundaries[2], _body_lengths_grid,
    n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
)   # [81, n_wl, N_MODES, N_MODES]
M_body = M_body_grid[_Li:_Li + 1]   # [1, n_wl, N_MODES, N_MODES] -- the settled design point
print(f"  body solved (shared MM_LENGTH_GRID reused, sliced at mm_length={_LENGTH_GRID[_Li]:.2f} um): {M_body.shape}")

DRIVEN = [_PORT_TO_INPUT_IDX[p] for p in _driven_ports]   # default: paper's ports 1 and 3 -- 0-indexed
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
#    Only meaningful (and only computed) when exactly the paper's own
#    I={1,3} were driven -- the cross-excitation metrics need BOTH.
#    Otherwise (e.g. MMI4X4_DRIVEN_PORTS=4, the "top" port): plain per-arm
#    |T|^2 and phase(S) for whatever WAS driven (validation doc Step 8,
#    Check 1 follow-up -- same-monitor-locations comparison against the
#    full-device FDTD side at a matched mm_length).
# ---------------------------------------------------------------------------
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}   # paper's phi~_5..8 (Eq. 15)

smat = {i: SMatrix(wavelengths=wavelengths) for i in DRIVEN}
for i in DRIVEN:
    for j in OUTPUTS:
        smat[i].set(coeffs[(i, j)], port=PORT_NAMES[j])
    smat[i].check_passive()

_all13 = set(_driven_ports) == {1, 3}

if _all13:
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
else:
    print(f"\nResults at mm_length={MM_LENGTH:.4f} µm (ΔL={DELTA_L:+.2f}), "
          f"driven ports {_driven_ports} (raw per-arm, no CMRR/Δφ -- needs ports 1&3):")
    for p in _driven_ports:
        i = _PORT_TO_INPUT_IDX[p]
        for name in PORT_NAMES:
            T = np.abs(smat[i].coeff(name)) ** 2
            phase_deg = np.angle(smat[i].coeff(name), deg=True)
            print(f"  port {p} -> {name}: |T|^2={T} phase={phase_deg} deg")

# ---------------------------------------------------------------------------
# 6. Save and plot. Filenames port-labeled by paper number (1-based) so an
#    arbitrary driven set never collides with the canonical ports-1&3 files.
# ---------------------------------------------------------------------------
for p in _driven_ports:
    smat[_PORT_TO_INPUT_IDX[p]].save(os.path.join(output_dir, f"smatrix_hybrid_port{p}.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix_hybrid_port{{{','.join(str(p) for p in _driven_ports)}}}.pkl")

if _all13:
    np.savez(os.path.join(output_dir, "hybrid_results.npz"),
             wavelengths=wavelengths, mm_length=MM_LENGTH, delta_L=DELTA_L,
             imbalance_in1=imbalance[0], imbalance_in3=imbalance[2],
             excess_loss_in1=excess_loss[0], excess_loss_in3=excess_loss[2],
             cmrr_I_in1=cmrr_I[0], cmrr_Q_in1=cmrr_Q[0],
             cmrr_I_in3=cmrr_I[2], cmrr_Q_in3=cmrr_Q[2],
             dphi_6=dphi[1], dphi_7=dphi[2], dphi_8=dphi[3],
             phi_5=phi[0], phi_6=phi[1], phi_7=phi[2], phi_8=phi[3])
    print(f"Results grids saved → {output_dir}/hybrid_results.npz")
else:
    np.savez(os.path.join(output_dir, "hybrid_results.npz"),
             wavelengths=wavelengths, mm_length=MM_LENGTH, delta_L=DELTA_L,
             driven_ports=np.array(_driven_ports),
             **{f"T_port{p}_{name}": np.abs(smat[_PORT_TO_INPUT_IDX[p]].coeff(name)) ** 2
                for p in _driven_ports for name in PORT_NAMES},
             **{f"phase_port{p}_{name}": np.angle(smat[_PORT_TO_INPUT_IDX[p]].coeff(name), deg=True)
                for p in _driven_ports for name in PORT_NAMES})
    print(f"Results grids saved → {output_dir}/hybrid_results.npz")

if _all13:
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
else:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    _colors = {"out5": "steelblue", "out6": "darkorange", "out7": "seagreen", "out8": "crimson"}
    for p in _driven_ports:
        i = _PORT_TO_INPUT_IDX[p]
        for name in PORT_NAMES:
            T = np.abs(smat[i].coeff(name)) ** 2
            ph = np.angle(smat[i].coeff(name), deg=True)
            axes[0].plot(wavelengths, T, "o-", label=f"port{p}→{name}", color=_colors[name])
            axes[1].plot(wavelengths, ph, "o-", label=f"port{p}→{name}", color=_colors[name])
    axes[0].set_xlabel("Wavelength (µm)"); axes[0].set_ylabel("|T|^2")
    axes[0].set_title("Per-arm transmission"); axes[0].legend(fontsize=7); axes[0].grid(True)
    axes[1].set_xlabel("Wavelength (µm)"); axes[1].set_ylabel("phase(S) (deg)")
    axes[1].set_title("Per-arm phase"); axes[1].legend(fontsize=7); axes[1].grid(True)
    plt.suptitle(f"4×4 MMI 90° hybrid — mm_length={MM_LENGTH:.2f} µm — "
                 f"driven ports {_driven_ports}", fontsize=10)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
    plt.close("all")
    print(f"Results plot saved → {output_dir}/hybrid_results.png")

# ---------------------------------------------------------------------------
# 8. LENGTH SWEEP -- metrics vs mm_length, out to 2.1x the nominal length.
#
# COSTS NO EXTRA SOLVER TIME. `M_body_grid` above is already the FULL shared
# 20-400 um / 0.05 um grid (7601 points); section 5 slices ONE index out of it
# and throws the rest away, and the einsum cascade there already carries a
# leading length axis `L` purely to drop it. This block keeps that axis.
#
# ⚠ STRIDE THE RESULT, NEVER RE-REQUEST A SUB-GRID. From
# mmi_4x4_90hybrid_calibrated_length_rescore.py's own hard-won note: "The
# length-sweep cache is keyed on the REQUESTED lengths array, so a strided
# request is a cache MISS even though every point in it was solved -- always
# ask for the exact full grid the cached runs used and stride the returned
# array instead." Asking solve_section_length_sweep() for a coarser grid here
# would silently re-solve every body mode (~416 s/mode x N_MODES x n_wl).
#
# Purpose: find mm_lengths where SEVERAL arms carry real power. The
# hybrid-vs-FDTD phase probes at L=5/10/15 um put 35-55% in out8 and 1e-4..1e-2
# in out5/out6, so a power-weighted phase metric there scores essentially one
# arm -- and a 90-degree hybrid's quadrature is a statement about arms that ALL
# carry light. This sweep is what picks the next full-device FDTD reference
# lengths, which cost hours each.
# ---------------------------------------------------------------------------
if _all13:
    SWEEP_MAX = 2.1 * _MM_LENGTH_DEFAULT                 # 334.53 um
    SWEEP_STEP = float(os.environ.get("MMI4X4_SWEEP_STEP", 0.25))
    _stride = max(1, int(round(SWEEP_STEP / 0.05)))
    _swp_L_all = _LENGTH_GRID[::_stride]
    _keep = _swp_L_all <= SWEEP_MAX + 1e-9
    sweep_L = _swp_L_all[_keep]
    M_body_swp = M_body_grid[::_stride][_keep]           # [n_L, n_wl, N, N]
    print(f"\nLength sweep: {len(sweep_L)} points, "
          f"{sweep_L[0]:.2f}-{sweep_L[-1]:.2f} um at {SWEEP_STEP:.2f} um "
          f"(strided from the cached grid -- no new solves)")

    # All FOUR excitations, not just the driven pair. The device metrics below
    # only ever need ports 1 and 3, but the full-device FDTD references drive
    # paper port 4 (input_idx=3), and a hybrid-vs-FDTD phase comparison has to
    # use the SAME excitation on both sides -- port 4's arm pattern is port 1's
    # MIRROR (out5<->out8, out6<->out7), so reading port 1 against a port-4
    # reference would silently compare mirrored arms. The extra einsums are
    # over already-cached matrices and cost nothing measurable.
    _sc_all = {}
    for i in range(4):
        _tmp = np.einsum("Lwab,wbk->Lwak", M_body_swp, M_in[i])
        for j in OUTPUTS:
            _sc_all[(i, j)] = np.einsum("wab,Lwbk->Lwak", M_out[j], _tmp)[:, :, 0, 0]
    _sc = {(i, j): _sc_all[(i, j)] for i in DRIVEN for j in OUTPUTS}
    # [in, out, length, wavelength] complex -- amplitude AND phase, so a
    # downstream comparison never has to re-run the cascade to get a phase.
    coeff_all = np.stack([np.stack([_sc_all[(i, j)] for j in OUTPUTS])
                          for i in range(4)])

    nL, nW = len(sweep_L), len(wavelengths)
    sw = {k: np.zeros((nL, nW)) for k in
          ("xl1", "xl3", "cmrrI1", "cmrrI3", "cmrrQ1", "cmrrQ3",
           "dphi6", "dphi7", "dphi8")}
    sw_T = {j: np.zeros((nL, nW)) for j in OUTPUTS}
    for li in range(nL):
        _s = {}
        for i in DRIVEN:
            _s[i] = SMatrix(wavelengths=wavelengths)
            for j in OUTPUTS:
                _s[i].set(_sc[(i, j)][li], port=PORT_NAMES[j])
        sw["xl1"][li] = _s[0].excess_loss_db(PORT_NAMES)
        sw["xl3"][li] = _s[2].excess_loss_db(PORT_NAMES)
        sw["cmrrI1"][li] = _s[0].cmrr_db("out5", "out8")
        sw["cmrrI3"][li] = _s[2].cmrr_db("out5", "out8")
        sw["cmrrQ1"][li] = _s[0].cmrr_db("out6", "out7")
        sw["cmrrQ3"][li] = _s[2].cmrr_db("out6", "out7")
        # Same A11 self/other order and mod-360 wrap as section 6 -- see the
        # comment there; deviating would silently flip the sign convention.
        _phi = {j: _s[2].cross_excitation_phase_diff_deg(_s[0], port=PORT_NAMES[j])
                for j in OUTPUTS}
        for j, key in ((1, "dphi6"), (2, "dphi7"), (3, "dphi8")):
            _raw = _phi[j] - _phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
            sw[key][li] = ((_raw + 180.0) % 360.0) - 180.0
        for j in OUTPUTS:
            sw_T[j][li] = np.abs(_s[0].coeff(PORT_NAMES[j])) ** 2

    def _band(ax, x, grid, color, label):
        """Mean line + min/max shading across the wavelength band."""
        ax.plot(x, grid.mean(axis=1), color=color, lw=1.3, label=label)
        ax.fill_between(x, grid.min(axis=1), grid.max(axis=1),
                        color=color, alpha=0.18, lw=0)

    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    ax = axes[0, 0]
    for g, c, lb in ((sw["cmrrI1"], "steelblue", "CMRR$_I$ p1"),
                     (sw["cmrrI3"], "navy", "CMRR$_I$ p3"),
                     (sw["cmrrQ1"], "darkorange", "CMRR$_Q$ p1"),
                     (sw["cmrrQ3"], "saddlebrown", "CMRR$_Q$ p3")):
        _band(ax, sweep_L, g, c, lb)
    ax.axhline(20.0, color="k", ls=":", lw=1.0)
    ax.set_ylabel("CMRR (dBe)"); ax.set_ylim(0, 60)
    ax.set_title("CMRR vs length  (dotted = 20 dBe requirement)", fontsize=10)

    ax = axes[0, 1]
    for g, c, lb in ((sw["dphi6"], "steelblue", "Δφ$_6$"),
                     (sw["dphi7"], "darkorange", "Δφ$_7$"),
                     (sw["dphi8"], "seagreen", "Δφ$_8$")):
        _band(ax, sweep_L, g, c, lb)
    for s in (-5.0, 5.0):
        ax.axhline(s, color="k", ls=":", lw=1.0)
    ax.set_ylabel("Phase error (deg)"); ax.set_ylim(-60, 60)
    ax.set_title("Phase error vs length  (dotted = ±5° spec)", fontsize=10)

    ax = axes[1, 0]
    _band(ax, sweep_L, sw["xl1"], "steelblue", "port 1")
    _band(ax, sweep_L, sw["xl3"], "crimson", "port 3")
    ax.set_ylabel("Excess loss (dB)"); ax.set_xlabel("mm_length (µm)")
    ax.set_title("Excess loss vs length (Eq. 11)", fontsize=10)

    ax = axes[1, 1]
    for j, c in zip(OUTPUTS, ("steelblue", "darkorange", "seagreen", "crimson")):
        _band(ax, sweep_L, 10 * np.log10(np.clip(sw_T[j], 1e-12, None)),
              c, PORT_NAMES[j])
    ax.axhline(10 * np.log10(0.25), color="k", ls="--", lw=0.9)
    ax.set_ylabel("per-arm |T|² (dB)"); ax.set_xlabel("mm_length (µm)")
    ax.set_ylim(-40, 2)
    ax.set_title("Per-arm loss vs length, driven port 1\n"
                 "(dashed = -6.02 dB, an ideal 4-way split)", fontsize=10)

    for ax in axes.ravel():
        ax.axvline(MM_LENGTH, color="magenta", ls="-", lw=1.0, alpha=0.7)
        ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
        ax.set_xlim(sweep_L[0], sweep_L[-1])
    plt.suptitle(
        f"4×4 MMI 90° hybrid — metrics vs body length, {N_WL} wavelengths "
        f"{WL_MIN}–{WL_MAX} µm   (magenta = nominal {MM_LENGTH:.2f} µm)\n"
        "line = band mean, shaded = min/max across wavelengths.  "
        "⚠ the SPEC is worst-case, so read the SHADED EDGE against the "
        "dotted lines, not the mean.", fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "hybrid_length_sweep.png"), dpi=200)
    plt.close("all")
    print(f"Length-sweep plot saved → {output_dir}/hybrid_length_sweep.png")

    np.savez(os.path.join(output_dir, "hybrid_length_sweep.npz"),
             lengths=sweep_L, wavelengths=wavelengths, mm_length=MM_LENGTH,
             **{f"sw_{k}": v for k, v in sw.items()},
             **{f"T_{PORT_NAMES[j]}": sw_T[j] for j in OUTPUTS},
             coeff=coeff_all, coeff_dims="[input_idx, output_idx, length, wavelength]",
             port_names=np.array(PORT_NAMES))
    print(f"Length-sweep data saved → {output_dir}/hybrid_length_sweep.npz")

    # THE POINT OF THE SWEEP: lengths where >=3 arms each carry >=5% of the
    # power, so their phases are well determined and a quadrature comparison
    # actually means something. These are the candidate FDTD reference lengths.
    _Tstack = np.stack([sw_T[j] for j in OUTPUTS])          # [4, nL, nW]
    _nstrong = (_Tstack.mean(axis=2) >= 0.05).sum(axis=0)   # band-mean per arm
    _multi = sweep_L[_nstrong >= 3]
    print(f"\nMULTI-ARM LENGTHS (>=3 arms at >=5% band-mean each): "
          f"{len(_multi)} of {nL} points")
    if len(_multi):
        _runs, _s0 = [], _multi[0]
        for a, b in zip(_multi[:-1], _multi[1:]):
            if b - a > SWEEP_STEP * 1.5:
                _runs.append((_s0, a)); _s0 = b
        _runs.append((_s0, _multi[-1]))
        for a, b in _runs:
            print(f"  {a:8.2f} - {b:8.2f} um")
    else:
        print("  none -- no length in range spreads power across 3+ arms")
