"""
Tutorial: standard directional coupler — FDTD/EME hybrid simulation (BEEF method).

Geometry (220 nm SOI, C-band):

    "in"  ___                                      ____  "through"
             \\                                    /
              \\__________________________________/        (coupling gap)
              /                                    \\
    ___      /                                      \\___  "cross"
   "isolated"

    |<- bend_length ->|<--- coupling_length --->|<- bend_length ->|

Three physically distinct features, but only **1 FDTD simulation + 1 EME
simulation total** (SYMMETRIC_MODE_REUSE=True, the default):

    [input S-bend]  ->  [straight coupled pair]  ->  [output S-bend]
       FDTD (x1)               EME (x1)              reused, not resimulated

  - The S-bends that carry the excited waveguide from a well-isolated access
    pitch down into the tight coupling gap (and back out again) are curved,
    non-adiabatic features — exactly what FDTD is for.
  - The straight coupled-pair region is where the actual directional
    coupling happens: an eigenmode expansion over its 2 supermodes (even +
    odd) is both exact and far cheaper than FDTD here.
  - The output section is NEVER separately simulated: DirectionalCoupler has
    no left/right asymmetry, so by reciprocity the input section's transfer
    matrix, transposed, IS the output section's.
  - The "isolated"-arm input matrix is ALSO not separately simulated: launching
    "in" vs. "isolated" is itself a mirror operation, so each body mode's
    coefficient either matches or flips sign between the two launches (even/odd
    mode symmetry) — see MODE_SIGNS below for how that sign vector was
    measured. Set SYMMETRIC_MODE_REUSE=False to simulate both arms directly
    (e.g. to re-derive MODE_SIGNS after a geometry change).

Both of those savings are DECLARED to HybridSolver, not hand-rolled. The
"through" path is an ordinary `mirror_pairs={2: 0}`. The "cross" path —
which pairs one arm's transposed matrix with the other arm's excitation —
is a `MirrorSpec` carrying the parity vector:

    mirror_pairs={2: MirrorSpec(section=0,
                                partner_kwargs={"input_idx": 1},
                                mode_signs=MODE_SIGNS)}

so both excitations still come from ONE FDTD simulation, now through the
supported API, and both declarations are re-validated against the geometry
on every call. `calibrate_kwargs` carries the interface-monitor pin to the
β calibration WITHOUT putting it in the body section's own kwargs — see the
comment there; it is what keeps the EME body's cached result valid.

Coupling-length sweep — the classic directional-coupler "beat length"
characterization — without any additional FDTD or EME simulations:

  - The input/output FDTD sections' transfer matrices don't depend on
    coupling_length (the S-bend/access-waveguide geometry is defined
    relative to the coupling-region edges — including the BODY_STUB µm of
    uniform straight each one folds in, see below — and
    FDTDSolver.solve_section() re-centers on its own padded domain). They're
    solved exactly once each,
    using one reference DirectionalCoupler instance built at the *middle* of
    the length sweep (also what the geometry/section preview plots use).
  - The EME body's mode profiles and effective indices depend only on the
    coupling cross-section (width/gap), not on how long the straight section
    is — so EMESolver.solve_section_length_sweep() computes the two
    supermodes once per wavelength and reuses them for every length in the
    sweep, each just a cheap diagonal phase-matrix evaluation
    (S_ii(L) = exp(2j*pi*neff_i/wavelength*L); see its docstring in
    beef/solvers/eme.py). Total solver cost is the same as a single fixed
    coupling_length run regardless of sweeping 200 lengths instead of 1.
  - Each FDTD input/output section reaches BODY_STUB µm past the
    S-bend/coupling-region seam, into the straight interior, so its
    interface mode decomposition happens on a uniform coupled-pair cross-
    section rather than right on the geometric discontinuity (same
    convention as the MMI examples' body_stub). The EME body therefore
    models coupling_length - 2*BODY_STUB, and the two stubs make up the
    difference — the total straight length simulated is still exactly the
    nominal coupling_length; see step 5's length audit.

HybridSolver itself has no idea this is a directional coupler; it only
knows each solver implements solve_section() and combines the resulting
per-section modal transfer matrices by matrix-chain multiplication. The
coupler-specific geometry (S-bends, coupling region) lives entirely in
DirectionalCoupler's section_prisms/section_gds/width_at/center_at methods.

Run:
    python /path/to/examples/directional_coupler/python/directional_coupler_hybrid.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import (
    DirectionalCoupler, EMESolver, FDTDSolver, HybridSolver, MirrorSpec, SMatrix,
    SOI_220, plot_beta_calibration, plot_hybrid_geometry,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "directional_coupler", "directional_coupler_hybrid"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
# The two engines get SEPARATE resolutions, and this is deliberate — a single
# shared RESOLUTION was how a 7% error hid for months (dc_validation.md item 11).
#
# meep's FDTD is converged at 40: its implied Δβ moves only 2.74% across res
# 30–90, non-monotonically (noise, not a trend). meow/MPB on its own Mesh2D is
# NOT — Δβ_EME at 1.55 µm runs 0.1125 (res 20) → 0.0957 (40) → 0.0926 (60) →
# 0.0894 (80) → 0.0895 (100), so 40 overstates the even/odd splitting by 7% and
# 80 is the first converged point. meep's subpixel smoothing is the likely
# reason it converges so much earlier on the same nominal grid.
#
# Measured directly in dc_validation.md Done item 6 (job 11489846) and pinned by
# tests/test_eme_convergence.py. Do not collapse these back into one knob.
FDTD_RESOLUTION = 40.0   # pixels/µm — converged; see above
EME_RESOLUTION = 80.0    # pixels/µm — 40 is NOT converged; see above
WL_MIN = 1.50       # wavelength sweep start (µm)
WL_MAX = 1.60       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

# Coupling length used for the wavelength-domain ("hybrid_results.png") plot
# and the SMatrix that gets saved. Set to a specific value (µm) to pin it;
# leave as None to auto-pick the length in coupling_lengths (defined below)
# whose coupling ratio is closest to a 50:50 (3 dB) split.
PLOT_COUPLING_LENGTH = 15.1   # µm — the tenth-of-a-micron round number nearest
# the DISPERSIVE 50:50 point, set 2026-08-03. Hardcoded rather than left on auto
# so this script and directional_coupler_fdtd.py's COUPLING_LENGTH stay pinned to
# the identical number instead of two independent auto-picks that could drift.
#
# Why 15.1 and not the older 15.13: 15.13 came from item 5's ground-truth sweep
# on the FROZEN material model. Dispersion defaulted 2026-08-01 and moved L₅₀
# down — measured 15.089958 µm (hybrid, job 11609658) and ~15.074 µm inferred for
# the full-device FDTD (job 11609659: band-mean ratio 0.500675 at 15.089958, over
# the hybrid's local slope of 0.0418 per µm). 15.1 is the nearest tenth to both.
#
# Whatever is pinned here is spliced into coupling_lengths (step 1) and asserted
# exact in step 5b. Before 2026-08-03 it was only argmin'd against the 0.299 µm
# sweep grid, which contains neither 15.13 nor 15.1 — the old pinned 15.13 was in
# fact simulated at 15.150754 µm while the FDTD example ran at a true 15.13, so
# the two examples sat 0.0208 µm apart despite both files reading "15.13".
#
# DC_PLOT_COUPLING_LENGTH overrides it without editing this file: a number pins
# that length, the literal "auto" re-derives it from this run's own sweep.
# Added 2026-08-01 for the dispersive re-measure — the then-pinned 15.13 was
# measured on the *frozen* material model (job 11544656), so re-using it would
# have assumed the very thing the re-run exists to check. Same override idiom as
# MMI1X2_CALIBRATE_BODY / DC_CALIBRATE_BODY.
_env_L = os.environ.get("DC_PLOT_COUPLING_LENGTH")
if _env_L is not None:
    PLOT_COUPLING_LENGTH = None if _env_L.strip().lower() == "auto" else float(_env_L)
    print(f"DC_PLOT_COUPLING_LENGTH={_env_L!r} -> PLOT_COUPLING_LENGTH={PLOT_COUPLING_LENGTH}")

# µm of the straight coupling region folded into each FDTD section, past the
# S-bend/coupling-region discontinuity — same convention (and same value) as
# mmi_2x2_hybrid.py's body_stub. Module-level because the coupling-length
# sweep grid below depends on it: no swept length may be shorter than the two
# stubs the FDTD sections already contain.
BODY_STUB = 0.25   # µm

# Propagate the EME body on the β the FDTD sections actually accumulate, rather
# than on meow's mode index. This is dc_validation.md item 6, and without it
# this example is wrong by ~4.5%.
#
# The two engines do not agree on how fast phase accumulates through the *same*
# uniform coupled pair. It is not a mode-solver mismatch — both sides run MPB
# and agree to 0.57% — it is MPB's mode index versus the phase meep's
# time-stepper accumulates: ~2.3% per mode, same sign for both supermodes. But
# coupling is set by Δβ, a *difference* of two β near 9.4, and β/Δβ ≈ 99 turns
# that common-mode 2.3% into a 10% error on Δβ. Since L₅₀ ∝ 1/Δβ, it lands
# directly on the number this script prints.
#
# MEASURED (job 11517245, the stub sweep with and without this flag):
#
#     drift with cut position   +0.1085 → −0.0007 µm/µm   (167x, below scatter)
#     L₅₀ vs full-device FDTD    −5.37%  → +0.11%
#
# The first is the one that matters structurally: with the calibration the
# answer stops depending on *where* the FDTD/EME cut is made, which is the
# property that makes a hybrid decomposition meaningful at all.
#
# HONEST CAVEAT (dc_validation.md item 13): this makes the hybrid reproduce the
# full-device FDTD, and that reference is meep's own discretization. MPB's Δβ
# disagrees with it by ~5% with *both* engines looking converged, and which is
# closer to Maxwell is unresolved. Read the agreement as self-consistency, not
# as proven physical accuracy.
#
# Cheap: three short FDTD section solves, cached and shared with the stub sweep.
CALIBRATE_BODY = True

# ---------------------------------------------------------------------------
# 1. Geometry — all device dimensions as variables (220 nm SOI).
#
#    coupling_lengths is the actual sweep; dc itself is built at the sweep's
#    middle value (bends/access-waveguide geometry don't depend on
#    coupling_length, so this single instance is equally valid as the
#    reference geometry for the FDTD sections at every length in the sweep,
#    and it's also exactly what the geometry/section preview plots should
#    show per the "plots at the middle length" convention).
# ---------------------------------------------------------------------------
#    The sweep starts at 2*BODY_STUB, not 0: each FDTD section already
#    contains BODY_STUB µm of the straight coupling region, so a device
#    shorter than both stubs put together can't be split into these three
#    sections at all (its EME body length would be negative).
coupling_lengths = np.linspace(2 * BODY_STUB, 60.0, 200)   # µm, the swept coupling length
L_mid = coupling_lengths[len(coupling_lengths) // 2]

#    The pinned plot length is spliced into the grid as an extra sample so the
#    argmin in step 5b lands on it *exactly*. Without this the 0.299 µm grid
#    silently snaps it to the nearest sample, and it holds no round number:
#    neither 15.13 nor today's 15.1 is on it. Measured — the pinned 15.13 was
#    being simulated at 15.150754 µm, 0.0208 µm from the length
#    directional_coupler_fdtd.py runs at, which defeats the whole point of
#    pinning the two examples to the same number. Extra lengths are free:
#    solve_section_length_sweep() solves the supermodes once per wavelength and
#    evaluates exp(2j*pi*neff*L) per length, so 201 lengths cost what 200 do.
#    L_mid is deliberately taken from the base grid *above* this splice: it is
#    only the reference geometry the FDTD sections are built on (bend shape is
#    independent of coupling_length), and holding it fixed keeps their cache
#    keys stable across changes to the pin.
if PLOT_COUPLING_LENGTH is not None:
    if PLOT_COUPLING_LENGTH < 2 * BODY_STUB:
        raise ValueError(
            f"PLOT_COUPLING_LENGTH ({PLOT_COUPLING_LENGTH} µm) is shorter than the "
            f"two FDTD body stubs ({2 * BODY_STUB} µm); the EME body length would be negative"
        )
    if not np.any(np.isclose(coupling_lengths, PLOT_COUPLING_LENGTH, rtol=0, atol=1e-9)):
        coupling_lengths = np.sort(np.append(coupling_lengths, PLOT_COUPLING_LENGTH))
        print(f"Spliced pinned plot length {PLOT_COUPLING_LENGTH:.6f} µm into the sweep grid "
              f"({len(coupling_lengths)} pts)")

dc = DirectionalCoupler(
    waveguide_width=0.45,    # single-mode strip width (µm)
    gap=0.25,                # coupling-region edge-to-edge gap (µm)
    coupling_length=L_mid,   # straight parallel interaction length (µm) — sweep reference value
    bend_length=9.0,         # horizontal extent of each S-bend (µm)
    waveguide_offset=3.0,    # far-end centre-to-centre port pitch (µm)
    bend_type="euler",       # "euler" (cubic-Bezier S-bend) or "circular" (matched arcs)
    waveguide_height=0.22,   # Si layer thickness (µm) — SOI_220
)
print(f"Directional coupler total length: {dc.total_length:.2f} µm "
      f"(coupling pitch {dc.coupling_pitch:.3f} µm, far pitch {dc.waveguide_offset:.2f} µm)")
print(f"Coupling length sweep: {coupling_lengths[0]:.1f} – {coupling_lengths[-1]:.1f} µm "
      f"({len(coupling_lengths)} pts, build reference length = {L_mid:.2f} µm)")

# ---------------------------------------------------------------------------
# 2. Configure sub-solvers — one per section.
#
#    input section  : input access waveguide + input S-bend + BODY_STUB µm
#                      of the straight coupled-pair region (enough for a clean
#                      interface monitor past the S-bend/coupling-region
#                      cross-section discontinuity)
#    body section   : propagates N_MODES supermodes through the remaining
#                      straight coupled-pair region
#                      (length = coupling_length - 2*BODY_STUB)
#    output section : the mirror image of the input section, one instance
#                      per output port ("through" / "cross")
# ---------------------------------------------------------------------------
N_MODES = 5   # even + odd supermodes carried at the two body interfaces.
              # Was 6: band 5 (n_eff 1.4367) sits below the 1.4443 cladding
              # index -- a radiation mode, dropped at no cost.
              # resolved/directional_coupler_validation.md check 2 caveat.

SYMMETRIC_MODE_REUSE = True   # derive M_in1 from M_in0 via MODE_SIGNS instead of a 2nd FDTD run

# Per-mode sign relating the "isolated"-arm launch to the "in"-arm launch
# (M_in1[m] = MODE_SIGNS[m] * M_in0[m]), measured once against a real
# 2-FDTD-arm run (cache_only, both directions). The two dominant supermodes
# (even, odd) matched to ~1e-5 relative residual; the three higher-order
# modes are individually noisier (~1% -- near-cutoff, tiny magnitude) but
# contribute negligibly. Fixed geometric property of the mirror-symmetric
# cross-section, not resolution/wavelength/length dependent. Re-derive by
# setting SYMMETRIC_MODE_REUSE = False if the geometry changes.
MODE_SIGNS = np.array([1, -1, -1, 1, 1])

_fdtd_common = dict(
    res=FDTD_RESOLUTION,
    sim_width=6.0,        # >= waveguide_offset + margin, so both far-pitch arms sit inside PML clearance
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,       # full 3-D FDTD: real 220 nm Si core in oxide, PML on all six faces
    # A full-sim_width source/norm aperture at a true I/O port would reach
    # past the shared y=0 axis to the *other* arm (only waveguide_offset
    # away), since both arms are always physically present in a section -
    # corrupting the isolated single-arm excitation. half_domain_aperture
    # confines the aperture to the near half of the domain (from the
    # symmetry axis to the near PML edge) instead.
    half_domain_aperture=True,
    material_stack=SOI_220,
    decay_tol=1e-11,
)

# One FDTDSolver instance covers both the input and output sections — they're
# mirror images of the same physical structure (see module docstring), so
# there's no need for two separately-configured objects.
fdtd_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholder length; solve_section_length_sweep() sweeps its own
    oxide_width=6.0, oxide_height=1.5, res=EME_RESOLUTION, material_stack=SOI_220,
    mode_solver="meep",  # matches FDTDSolver's own eigenmode basis in this hybrid chain
)

# ---------------------------------------------------------------------------
# 3. Build HybridSolver — 3-section chain: FDTD bend | EME coupling | FDTD bend.
#    Step 5's numeric sweep now runs THROUGH these solvers rather than around
#    them: the cross path's "one arm's transpose against the other arm's
#    excitation" is declared as a MirrorSpec instead of hand-rolled. See the
#    module docstring.
# ---------------------------------------------------------------------------
half_body = dc.coupling_length / 2
half_total = dc.total_length / 2

boundaries = [
    -half_total,               # physical input port (access waveguide end)
    -half_body + BODY_STUB,    # input section <-> body section interface (BODY_STUB past the seam)
    +half_body - BODY_STUB,    # body section <-> output section interface (BODY_STUB past the seam)
    +half_total,               # physical output port (access waveguide end)
]
interface_modes = [N_MODES, N_MODES]

# EME body length per swept coupling_length: each FDTD section folds BODY_STUB
# µm of the straight coupling region into itself (boundaries[1]/[2] sit that
# far past the seam), so the EME body's physical length is coupling_length
# minus both stubs. The two asserts and the prints are the length audit — what
# the three sections actually simulate must add back up to the nominal
# coupling_length, exactly — and they run here, before step 5's FDTD solves,
# so a bad BODY_STUB/sweep-range combination fails immediately instead of
# after hours of simulation.
body_lengths = coupling_lengths - 2 * BODY_STUB
# >= 0, not > 0: the shortest swept length is exactly 2*BODY_STUB, i.e. the
# degenerate device whose whole coupling region lives in the two FDTD sections
# and whose EME body is zero-length (an identity propagation matrix — verified,
# meow's cell_lengths=[0.0] returns exactly I).
assert body_lengths[0] >= 0, (
    f"coupling_lengths must start at >= 2*BODY_STUB ({2 * BODY_STUB:.2f} µm); "
    f"got {coupling_lengths[0]:.2f} µm"
)
assert np.allclose(body_lengths + 2 * BODY_STUB, coupling_lengths)
print(f"Coupling-length budget: 2 x {BODY_STUB:.2f} µm inside the FDTD sections "
      f"+ EME body {body_lengths[0]:.2f}–{body_lengths[-1]:.2f} µm "
      f"= {coupling_lengths[0]:.2f}–{coupling_lengths[-1]:.2f} µm total")
print(f"FDTD section span: {boundaries[1] - boundaries[0]:.2f} µm "
      f"(+{fdtd_solver.waveguide_length:.2f} µm access-waveguide stub), "
      f"interface monitor at x = {boundaries[1]:.2f} µm "
      f"({BODY_STUB:.2f} µm past the seam at x = {-half_body:.2f} µm)")

# Interface-monitor pin for the input-section FDTD runs (and their section
# plots): span the full non-PML domain — PML inner edge to PML inner edge —
# centered on the coupled pair's shared y=0 axis, matching the EME body's own
# y=0-centered supermode basis it is matrix-multiplied against. With the
# BODY_STUB offset the section now ends *inside* the coupling region, where
# center_at()/width_at() already report the coupled pair (0.0 / 2*ww+gap)
# rather than the excited arm's ∓coupling_pitch/2 S-bend centreline, so this
# no longer *corrects* the placement the way it did with the interface sitting
# exactly on the seam — it pins it. Still passed explicitly: it documents the
# intent, it survives a change of BODY_STUB back toward the seam, and a
# non-None monitor_size is what lets solve_section() reuse a cached larger-
# mode-count run by slicing (see its family_config/lookup_superset comment).
_iface_mon = dict(monitor_center=0.0, monitor_size=fdtd_solver.sim_width)

# One HybridSolver per excitation configuration; they differ ONLY in how
# section 2 is declared. "through" measures the same arm that was excited, so
# it is an ordinary x-mirror of section 0. "cross" measures the OTHER arm,
# which is the x-mirror of an excitation section 0 never ran — MirrorSpec
# states that explicitly, reusing section 0 reinterpreted as the isolated arm
# via MODE_SIGNS. Both are re-validated against the geometry every call.
#
# calibrate_kwargs carries _iface_mon to the calibration WITHOUT putting it in
# section 1's geometry_kwargs. That distinction is load-bearing, not stylistic:
# EMESolver.solve_section_length_sweep() funnels unknown **kwargs into
# geometry_kwargs, which is in the cache fingerprint, so pinning the monitor on
# the body sweep too would change its hash (measured: 0c243abea80cb883 ->
# 76f7367470c61021) and force a re-solve of an already-cached result.
_dc_shared = dict(
    solvers=[fdtd_solver, body_solver, fdtd_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
    calibrate_body_neff=CALIBRATE_BODY,
    calibrate_kwargs=_iface_mon,
)

if SYMMETRIC_MODE_REUSE:
    _cross_mirror = {2: MirrorSpec(section=0,
                                   partner_kwargs={"input_idx": 1},
                                   mode_signs=MODE_SIGNS)}
    _cross_kwargs = {0: {"input_idx": 0, **_iface_mon},
                     2: {"output_idx": 1, **_iface_mon}}
else:
    # Simulate the isolated arm for real: section 0 excited on input_idx=1,
    # section 2 its plain x-mirror. The 2-FDTD-run fallback, and what
    # re-derives MODE_SIGNS below.
    _cross_mirror = {2: 0}
    _cross_kwargs = {0: {"input_idx": 1, **_iface_mon},
                     2: {"output_idx": 1, **_iface_mon}}

_through_kwargs = {0: {"input_idx": 0, **_iface_mon},
                   2: {"output_idx": 0, **_iface_mon}}

hybrid_through = HybridSolver(mirror_pairs={2: 0}, **_dc_shared)
hybrid_cross = HybridSolver(mirror_pairs=_cross_mirror, **_dc_shared)

# Kept under the old name for the geometry/section preview plots below, which
# read only boundaries/solvers/mirror_pairs off it.
hybrid = hybrid_through

# ---------------------------------------------------------------------------
# 4. Geometry plot — shows which solver covers which section, at dc's
#    (middle-of-sweep) coupling_length.
# ---------------------------------------------------------------------------
plot_hybrid_geometry(
    dc, solver=hybrid,
    in_modes=1, out_modes=1,
    title=f"Directional Coupler — Hybrid Solver Sections (coupling_length={L_mid:.1f} µm)",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True,
    sub_region_show_eme_detail=True,
    # Both input arms shown for reference; "isolated" is only actually
    # simulated in step 5 when SYMMETRIC_MODE_REUSE=False (see there). Output
    # section (index 2) is the x-mirror of the input section either way, so
    # it gets no sub-region plot — hatched/labeled as mirrored above instead.
    sub_region_configs=[
        ("in", {0: {"input_idx": 0, **_iface_mon}}, {2: 0}),
        ("isolated", {0: {"input_idx": 1, **_iface_mon}}, {2: 0}),
    ],
)
print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")
print(f"Sub-region plots saved → {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run the sweep — input section: 1 FDTD solve (SYMMETRIC_MODE_REUSE) or 2
#    (both arms). Body: 1 EME mode-solve/wavelength, length-swept
#    analytically. Output section: never simulated, mirror of the input
#    (see module docstring) — "through" mirrors "in", "cross" mirrors
#    "isolated".
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nWavelength sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

_section_cache = {}   # shared across both excitations: section 0's kwargs are
                      # identical between them, so it is simulated once.

# Calibrate the body's β off the FDTD — see CALIBRATE_BODY at the top for why
# this is not optional. Solved on the body's own span, the same [boundaries[1],
# boundaries[2]] the length sweep below propagates through, so the two are
# talking about the same cross-section by construction. β is a property of that
# cross-section and not of the swept length, so this is one solve for every
# length in the sweep.
full_through = hybrid_through.solve_length_sweep(
    dc, wavelengths, body_lengths, sweep_section=1,
    geometry_kwargs=_through_kwargs, cache=_section_cache, save_dir=save_dir,
)   # [n_L, n_wl, 1, 1]
full_cross = hybrid_cross.solve_length_sweep(
    dc, wavelengths, body_lengths, sweep_section=1,
    geometry_kwargs=_cross_kwargs, cache=_section_cache, save_dir=save_dir,
)   # [n_L, n_wl, 1, 1]

if CALIBRATE_BODY:
    # The measurement behind the body's β. Only the "through" solver holds it:
    # the cross call's body section came from the shared cache, so its own
    # calibration never re-ran — by design, since a hit means the matrix this
    # pass used IS the one that calibration produced.
    _calibration = hybrid_through.last_calibrations[1]
    _beta = _calibration.beta  # [n_wl, N_MODES] — rad/µm, from the FDTD's own phase vs length
    _body_neff = _beta * wavelengths[:, None] / (2 * np.pi)

    # What was actually simulated, and whether each cell's phase really is
    # linear in length — the residual guard only raises past 0.1 rad, so a
    # marginal fit passes silently and is invisible in `_beta`. No extra
    # simulation: these plot the data the call above already produced.
    plot_beta_calibration(
        _calibration, dc, fdtd_solver,
        os.path.join(output_dir, "calibration"), **_iface_mon,
    )
    print(f"Calibration plots saved → {output_dir}/calibration/ (geometry/ + fits/)")
    _i = int(np.argmin(np.abs(wavelengths - 1.55)))
    # The mode-solver number to compare against, for the log line only. Solved
    # on a SHORT window at the body midpoint, not on [boundaries[1],
    # boundaries[2]]: section_mode_neff() pads the cut by dpml on each side, so
    # the full body span would ask for a cut past the coupling region and land
    # inside an S-bend, which DirectionalCoupler rejects. The cross-section is
    # uniform across the whole body, so any interior window gives the same modes.
    _seed_half = min(1.0, half_body - 2.0 * fdtd_solver.dpml)
    _mpb = fdtd_solver.section_mode_neff(
        dc, wavelengths[_i:_i + 1], -_seed_half, +_seed_half, N_MODES, **_iface_mon)[0]
    print(f"Body β calibration at λ = {wavelengths[_i]:.3f} µm "
          f"(even/odd supermodes):")
    print(f"    MPB mode index  n_eff {_mpb[0]:.6f} {_mpb[1]:.6f}   "
          f"Δn_eff {_mpb[0] - _mpb[1]:.6f}")
    print(f"    FDTD own phase  n_eff {_body_neff[_i, 0]:.6f} {_body_neff[_i, 1]:.6f}   "
          f"Δn_eff {_body_neff[_i, 0] - _body_neff[_i, 1]:.6f}   <- body propagates on this")

if not SYMMETRIC_MODE_REUSE:
    # Re-derive MODE_SIGNS. Both arms were really simulated in this
    # configuration, so the implied parity vector reads straight off the two
    # input-section matrices the two passes actually used — no third solve,
    # and no reconstructing the calls by hand.
    _m_in = hybrid_through.last_section_matrices[0][:, :, 0]     # [n_wl, N_MODES]
    _m_iso = hybrid_cross.last_section_matrices[0][:, :, 0]
    _ratio = _m_iso / np.where(np.abs(_m_in) > 0, _m_in, np.nan)
    _implied = np.sign(np.real(np.nanmean(_ratio, axis=0))).astype(int)
    print("\nMODE_SIGNS re-derivation (SYMMETRIC_MODE_REUSE=False):")
    print(f"  hard-coded : {list(MODE_SIGNS)}")
    print(f"  implied    : {list(_implied)}")
    print(f"  agree      : {bool(np.array_equal(_implied, MODE_SIGNS))}")

T_through_grid = np.abs(full_through[:, :, 0, 0]) ** 2   # [n_L, n_wl]
T_cross_grid = np.abs(full_cross[:, :, 0, 0]) ** 2        # [n_L, n_wl]
T_total_grid = T_through_grid + T_cross_grid

# ---------------------------------------------------------------------------
# 5b. Choose the coupling length used for the wavelength-domain plot/SMatrix.
#     PLOT_COUPLING_LENGTH pins it explicitly (snapped to the nearest sampled
#     length); left as None, auto-pick the shortest length whose mean-over-
#     wavelength cross/through ratio crosses 0.5 (the "3 dB coupler length").
# ---------------------------------------------------------------------------
ratio_mean = T_cross_grid.mean(axis=1) / np.clip(T_total_grid.mean(axis=1), 1e-20, None)
sign_changes = np.flatnonzero(np.diff(np.sign(ratio_mean - 0.5)))
if sign_changes.size:
    _i = sign_changes[0]
    L50_idx = _i if abs(ratio_mean[_i] - 0.5) <= abs(ratio_mean[_i + 1] - 0.5) else _i + 1
else:
    L50_idx = int(np.argmin(np.abs(ratio_mean - 0.5)))
L50 = float(coupling_lengths[L50_idx])

# Sub-grid refinement: the sweep grid is ~0.3 µm coarse, so the snapped L50
# above is grid-limited rather than measurement-limited. Linearly interpolate
# the 0.5 crossing between the bracketing samples when there is one. This is
# the number worth handing to a full-device FDTD run, which is not on any grid.
if sign_changes.size:
    _a, _b = sign_changes[0], sign_changes[0] + 1
    _r0, _r1 = ratio_mean[_a], ratio_mean[_b]
    if _r1 != _r0:
        L50 = float(coupling_lengths[_a]
                    + (0.5 - _r0) * (coupling_lengths[_b] - coupling_lengths[_a]) / (_r1 - _r0))

if PLOT_COUPLING_LENGTH is not None:
    L_plot_idx = int(np.argmin(np.abs(coupling_lengths - PLOT_COUPLING_LENGTH)))
    # Exact by construction — step 1 splices the pin into coupling_lengths.
    # Printed at 6 decimals, not 2: the old .2f rounding is what let a 0.0208 µm
    # snap error read as "15.13 -> 15.15" and go unnoticed.
    assert np.isclose(coupling_lengths[L_plot_idx], PLOT_COUPLING_LENGTH, rtol=0, atol=1e-9), (
        f"pinned plot length {PLOT_COUPLING_LENGTH} µm was not sampled exactly "
        f"(got {coupling_lengths[L_plot_idx]} µm)"
    )
    print(f"Plot length: user-set {PLOT_COUPLING_LENGTH:.6f} µm -> "
          f"sampled length {coupling_lengths[L_plot_idx]:.6f} µm (exact)")
else:
    L_plot_idx = L50_idx
    print(f"Plot length: auto-picked {coupling_lengths[L_plot_idx]:.2f} µm "
          f"(50:50 split, ratio = {ratio_mean[L_plot_idx]:.3f})")
L_plot = coupling_lengths[L_plot_idx]

# Always emit the measured 50:50 length, pinned or not — it is this script's
# headline physical result and was previously only reachable by parsing stdout
# or redoing this argmin by hand downstream. directional_coupler_fdtd.py reads
# it via DC_COUPLING_LENGTH so the two examples cannot drift apart.
print(f"Measured 50:50 length: {L50:.6f} µm "
      f"(grid-snapped {coupling_lengths[L50_idx]:.6f} µm, "
      f"ratio {ratio_mean[L50_idx]:.4f}; "
      f"{'interpolated between samples' if sign_changes.size else 'NO 0.5 CROSSING in sweep — nearest point only'})")
with open(os.path.join(output_dir, "coupling_length_50.txt"), "w") as _fh:
    _fh.write(f"{L50:.6f}\n")

# Slice at the chosen coupling_length for the wavelength-domain view — keep
# the complex coefficients (not just power) for this slice so the saved
# SMatrix preserves phase, same as before.
coeff_through = full_through[L_plot_idx, :, 0, 0]
coeff_cross = full_cross[L_plot_idx, :, 0, 0]
T_through = T_through_grid[L_plot_idx]
T_cross = T_cross_grid[L_plot_idx]
T_total = T_total_grid[L_plot_idx]

wl_c = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_c))

print(f"Results at coupling_length = {L_plot:.2f} µm, λ = {wl_c:.3f} µm:")
print(f"  Through output T = {T_through[idx]:.4f}  IL = {10*np.log10(max(T_through[idx],1e-20)):.2f} dB")
print(f"  Cross output   T = {T_cross[idx]:.4f}  IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total          T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  Coupling ratio   = {T_cross[idx]/max(T_total[idx],1e-20):.3f}  "
      f"(0 = fully through, 1 = fully cross)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat = SMatrix(wavelengths=wavelengths)
smat.set(coeff_through, port="through")
smat.set(coeff_cross, port="cross")
smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix saved (coupling_length={L_plot:.2f} µm slice) → {output_dir}/smatrix_hybrid.pkl")

np.savez(
    os.path.join(output_dir, "length_sweep.npz"),
    coupling_lengths=coupling_lengths, wavelengths=wavelengths,
    T_through=T_through_grid, T_cross=T_cross_grid,
    L50=L50, L50_idx=L50_idx, ratio_mean=ratio_mean, L_plot=L_plot,
)
print(f"Length-sweep grid saved → {output_dir}/length_sweep.npz")

# -- Wavelength-domain plot, at the reference coupling_length --------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_through, label="Through", color="steelblue")
ax.plot(wavelengths, T_cross,   label="Cross",   color="darkorange", linestyle="--")
ax.plot(wavelengths, T_total,   label="Total",   color="k", lw=1.5)
ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title(f"Transmission vs Wavelength (coupling_length={L_plot:.2f} µm)")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)

ax = axes[1]
IL_through = 10 * np.log10(np.clip(T_through, 1e-20, None))
IL_cross   = 10 * np.log10(np.clip(T_cross,   1e-20, None))
IL_total   = 10 * np.log10(np.clip(T_total,   1e-20, None))
ax.plot(wavelengths, IL_through, label="Through", color="steelblue")
ax.plot(wavelengths, IL_cross,   label="Cross",   color="darkorange", linestyle="--")
ax.plot(wavelengths, IL_total,   label="Total",   color="k", lw=1.5)
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title(f"Insertion Loss vs Wavelength (coupling_length={L_plot:.2f} µm)")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle("Directional Coupler — FDTD/EME Hybrid Simulation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Wavelength-domain results plot saved → {output_dir}/hybrid_results.png")

# -- Coupling-length sweep plot, mean/min/max across the wavelength sweep --
fig, ax = plt.subplots(figsize=(7, 4.5))

for grid, label, color in [
    (T_through_grid, "Through", "steelblue"),
    (T_cross_grid,   "Cross",   "darkorange"),
    (T_total_grid,   "Total",   "k"),
]:
    mean, lo, hi = grid.mean(axis=1), grid.min(axis=1), grid.max(axis=1)
    ax.plot(coupling_lengths, mean, label=label, color=color,
            lw=1.5 if label == "Total" else 1.2,
            linestyle="--" if label == "Cross" else "-")
    ax.fill_between(coupling_lengths, lo, hi, color=color, alpha=0.15, lw=0)

ax.axvline(L_plot, color="gray", linestyle=":", lw=0.8,
           label=f"Plotted length ({L_plot:.2f} µm, 50:50 split)")
ax.set_xlabel("Coupling length (µm)")
ax.set_ylabel("Transmission (mean over wavelength, shaded = min–max)")
ax.set_title("Directional Coupler — Power Transfer vs Coupling Length")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "coupling_length_sweep.png"), dpi=300)
plt.close("all")
print(f"Coupling-length sweep plot saved → {output_dir}/coupling_length_sweep.png")
