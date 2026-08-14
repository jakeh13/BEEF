"""4x4 MMI 90-degree hybrid -- coarse waveguide_offset scoping sweep.

validation_markdown/mmi_4x4_90_hybrid_validation.md, "Open decisions" item
5 (queued 2026-08-13). At the converged N_MODES_BODY=16 design point (job
11778128) this device fails the paper's own cited system requirements
(CMRR >= 20 dBe, |Delta phi_j| <= 5 deg) at most wavelengths. waveguide_
offset is currently DERIVED (= mm_width/4, chosen only because it exactly
reproduces the paper's Eq. 8 port-position formula through beef's default
MMI port generator -- Step 2) rather than independently validated for this
strongly-guided SOI/O-band platform. Since offset controls exactly where
each access waveguide couples into the body's self-imaging pattern, it is
a plausible, cheap, single-parameter lever to screen before considering a
larger geometry redesign.

Made affordable by two zero-new-cluster-time checks recorded in the
validation doc (same date, both re-run here as live assertions rather than
assumed):

1. The EME body-mode solve is independent of waveguide_offset -- offset
   only moves the four ports' y-position, entirely inside the *input*
   section's own x-range; the body section's cache fingerprint is built
   from geometry.section_gds() sliced to the body's own x-range only.
   Confirmed directly via solve_section_length_sweep(..., cache_only=True)
   against job 11768847's cache with a never-before-solved offset -- cache
   HIT. So the body is fetched from cache here, at every offset, never
   re-solved.
2. mmi_2x2_hybrid.py's own SYMMETRIC_MODE_REUSE/MODE_SIGNS top-bottom
   mirror-reuse trick generalizes cleanly to this device's two port-pairs
   (input_idx 0<->3 at y=-2.25/+2.25, 1<->2 at y=-0.75/+0.75): checked
   against job 11778128's cached results with proper circular statistics
   (a first pass using naive angle mean/std wrongly flagged 2 modes as
   unstable -- a branch-cut artifact, corrected), all 16 body modes give a
   clean, real, wavelength-stable +/-1 sign, the SAME shared vector for
   both pairs. So only 2 NEW FDTD input solves (input_idx=0, input_idx=1)
   are needed per new offset point; input_idx=3/2 are derived via that
   sign vector, not separately solved.

The sign vector itself is derived HERE, at runtime, from the baseline
(offset=mm_width/4=1.5, job 11778128) 4-solve results -- not hardcoded --
so a stale literal can never silently drift from what job 11778128
actually measured. At production settings this hits that job's own cache
(zero new compute); at reduced test settings (below) there is no matching
cache entry, so it solves 4 cheap small-structure baseline runs fresh
instead.

Production settings (res=40, decay_tol=1e-11) are kept for the 2 new
solves per offset point, NOT loosened for speed -- this repo's own
2026-08-07 lesson (the taper_length investigation's coarse RES=15/
N_MODES_BODY=4 proxy extrapolated a ~0.3-0.6 dB effect that measured
~0.08-0.10 dB at real settings) is why "coarse" here means fewer offset
points and an approximate (sign-derived) reuse, not a lower-fidelity
solve. Any offset this sweep flags as promising needs a full 4-solve
confirmation (like job 11778128) before being trusted as a reported
number -- this script's own numbers are a screening signal, not a final
one, precisely because the reuse trick, however clean, is being applied
here to an offset it was never independently re-validated at (the doc's
own open caveat).

Usage (one offset per invocation -- SLURM array convention, one task per
offset point, same pattern as every other costly sweep in this repo):
    python mmi_4x4_90hybrid_offset_sweep.py <offset_um>
    python mmi_4x4_90hybrid_offset_sweep.py            # baseline only, cache-only, free

Env overrides for a cheap plumbing smoke test (mirrors HCG_TEST_MODE /
MMI4X4_RES precedent elsewhere in this repo):
    MMI4X4_OFFSET_RES, MMI4X4_OFFSET_DECAY_TOL, MMI4X4_OFFSET_NWL
"""
import os
import sys
import numpy as np

from beef import MMI, EMESolver, FDTDSolver, SMatrix, SOI_220
from beef.cache.store import CacheMissError

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "offset_sweep"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Fixed geometry/solver settings -- IDENTICAL to mmi_4x4_90hybrid_hybrid.py's
# settled Step 6c design point (job 11778128), except waveguide_offset,
# which is this sweep's own free parameter.
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MMI4X4_OFFSET_RES", 40.0))
DECAY_TOL = float(os.environ.get("MMI4X4_OFFSET_DECAY_TOL", 1e-11))
WL_MIN, WL_MAX = 1.26, 1.36
N_WL = int(os.environ.get("MMI4X4_OFFSET_NWL", 3))

WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 6.0
BASELINE_OFFSET = MM_WIDTH / 4.0   # = 1.5 -- job 11778128's own design point
TAPER_WIDTH = 0.8
TAPER_LENGTH = 1.0
L_PI = 108.9195
DELTA_L = -2.75
MM_LENGTH = 0.75 * L_PI + DELTA_L   # = 78.9396 -- UNCHANGED across this sweep

N_MODES = 16

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)


def make_mmi(offset: float) -> MMI:
    return MMI(
        n_inputs=4, n_outputs=4,
        mm_length=MM_LENGTH, mm_width=MM_WIDTH,
        waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
        waveguide_offset=offset, waveguide_height=0.22,
        use_neighbor_aperture=True,
    )


def make_solvers():
    fdtd_common = dict(
        res=RESOLUTION, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_220, decay_tol=DECAY_TOL,
    )
    input_solver = FDTDSolver(**fdtd_common)
    body_solver = EMESolver(
        num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
        oxide_width=fdtd_common["sim_width"], oxide_height=fdtd_common["sim_height"],
        res=RESOLUTION, material_stack=SOI_220, mode_solver="meep",
    )
    return input_solver, body_solver


def boundaries_for(mmi: MMI):
    half_body = mmi.mm_length / 2
    half_total = mmi.total_length / 2
    body_stub = 0.25
    return [-half_total, -half_body + body_stub, +half_body - body_stub, +half_total], body_stub


def solve_input(input_solver, mmi, boundaries, k, cache_only):
    return input_solver.solve_section(
        mmi, wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=input_solver.waveguide_length, stub_right=0.0,
        input_idx=k, save_dir=save_dir, cache_only=cache_only,
    )


def derive_sign_vector(input_solver, body_solver):
    """Re-derive the top-bottom mirror sign vector from the baseline
    (offset=1.5) 4-solve results, at whatever settings this run is using.
    cache_only=False: at production settings this is a cache HIT against
    job 11778128 (zero new compute, matches the doc's own finding); at
    reduced test settings there is no matching cache entry, so this
    solves 4 cheap small-structure runs fresh instead of raising --
    correct either way, not just an optimization for the production case.
    Returns a real +/-1 array, shape (N_MODES,)."""
    mmi0 = make_mmi(BASELINE_OFFSET)
    boundaries0, _ = boundaries_for(mmi0)
    M = [solve_input(input_solver, mmi0, boundaries0, k, cache_only=False)[..., 0] for k in range(4)]
    # Circular-consistent sign estimate (a naive angle mean/std is wrong near
    # +/-180 deg -- see the validation doc's corrected finding).
    unit03 = (M[3] / M[0]); unit03 /= np.abs(unit03)
    unit12 = (M[2] / M[1]); unit12 /= np.abs(unit12)
    mean_unit = np.concatenate([unit03, unit12], axis=0).mean(axis=0)  # average both pairs, all wl
    R = np.abs(mean_unit)
    sign = np.sign(np.cos(np.angle(mean_unit)))
    assert np.all(R > 0.999), f"sign vector not clean at baseline: min R={R.min():.6f}"
    print(f"  Derived sign vector from baseline (min circular-R={R.min():.6f}): {sign.astype(int)}")
    return sign, M


def solve_offset_point(offset: float):
    input_solver, body_solver = make_solvers()
    mmi = make_mmi(offset)
    boundaries, body_stub = boundaries_for(mmi)
    print(f"\n=== offset={offset:.4f} um (mm_length={mmi.mm_length:.4f}, "
          f"input_ys={mmi.input_ys}) ===")

    if abs(offset - BASELINE_OFFSET) < 1e-9:
        # The exact job-11778128 design point -- full real data already
        # cached, no reuse trick needed.
        M_in = [solve_input(input_solver, mmi, boundaries, k, cache_only=True) for k in range(4)]
        print("  baseline point: fetched all 4 real solves from cache (job 11778128)")
    else:
        sign, _ = derive_sign_vector(input_solver, body_solver)
        M0 = solve_input(input_solver, mmi, boundaries, 0, cache_only=False)
        M1 = solve_input(input_solver, mmi, boundaries, 1, cache_only=False)
        M3 = M0 * sign[None, :, None]
        M2 = M1 * sign[None, :, None]
        M_in = [M0, M1, M2, M3]
        print("  solved input_idx=0,1 for real; derived 2,3 via sign vector")

    M_out = [Mk.transpose(0, 2, 1) for Mk in M_in]

    _DELTA_L_GRID = np.linspace(-10.0, 10.0, 81)
    _body_lengths_grid = (0.75 * L_PI + _DELTA_L_GRID) - 2 * body_stub
    _Li = int(np.argmin(np.abs(_DELTA_L_GRID - DELTA_L)))
    _body_kwargs = dict(n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir)
    try:
        M_body_grid = body_solver.solve_section_length_sweep(
            mmi, wavelengths, boundaries[1], boundaries[2], _body_lengths_grid,
            cache_only=True, **_body_kwargs,
        )
        print(f"  body fetched from cache (never re-solved for this sweep)")
    except CacheMissError as e:
        if RESOLUTION == 40.0 and abs(DECAY_TOL - 1e-11) < 1e-15:
            # Production settings: job 11768847's grid is expected to hit
            # unconditionally (confirmed offset-independent -- validation
            # doc's "waveguide_offset body-cache independence" finding). A
            # miss here means something upstream changed unexpectedly --
            # fail loudly rather than silently eating a ~77 min resolve.
            raise RuntimeError(
                "body cache MISS at production settings (res=40, decay_tol=1e-11) "
                "-- expected an unconditional hit against job 11768847's grid. "
                "Something changed (mm_length/N_MODES/oxide dims?); investigate "
                "before letting this silently trigger a fresh ~77 min solve."
            ) from e
        # Reduced test settings: no matching cache entry is expected (this
        # body basis was never solved at this resolution before) -- solve
        # it for real, small/cheap at test settings.
        M_body_grid = body_solver.solve_section_length_sweep(
            mmi, wavelengths, boundaries[1], boundaries[2], _body_lengths_grid,
            cache_only=False, **_body_kwargs,
        )
        print(f"  body cache MISS at test settings (expected) -- solved fresh")
    M_body = M_body_grid[_Li:_Li + 1]
    print(f"  body M shape: {M_body.shape}")

    DRIVEN = [0, 2]
    OUTPUTS = [0, 1, 2, 3]
    PORT_NAMES = ["out5", "out6", "out7", "out8"]
    TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}

    coeffs = {}
    for i in DRIVEN:
        tmp = np.einsum("Lwab,wbk->Lwak", M_body, M_in[i])
        for j in OUTPUTS:
            c = np.einsum("wab,Lwbk->Lwak", M_out[j], tmp)
            coeffs[(i, j)] = c[0, :, 0, 0]

    smat = {i: SMatrix(wavelengths=wavelengths) for i in DRIVEN}
    for i in DRIVEN:
        for j in OUTPUTS:
            smat[i].set(coeffs[(i, j)], port=PORT_NAMES[j])
        smat[i].check_passive()

    imbalance = {i: smat[i].imbalance_db(PORT_NAMES) for i in DRIVEN}
    excess_loss = {i: smat[i].excess_loss_db(PORT_NAMES) for i in DRIVEN}
    cmrr_I = {i: smat[i].cmrr_db("out5", "out8") for i in DRIVEN}
    cmrr_Q = {i: smat[i].cmrr_db("out6", "out7") for i in DRIVEN}

    # A11's sign fix: self/other swapped relative to the paper's literal Eq. 15.
    phi = {j: smat[2].cross_excitation_phase_diff_deg(smat[0], port=PORT_NAMES[j]) for j in OUTPUTS}
    dphi = {}
    for j in [1, 2, 3]:
        raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
        dphi[j] = ((raw + 180.0) % 360.0) - 180.0

    max_dphi = float(np.max(np.abs(np.concatenate([dphi[1], dphi[2], dphi[3]]))))
    min_cmrr = float(np.min(np.concatenate([cmrr_I[0], cmrr_Q[0], cmrr_I[2], cmrr_Q[2]])))
    max_excess_loss = float(np.max(np.concatenate([excess_loss[0], excess_loss[2]])))
    max_imbalance = float(np.max(np.concatenate([imbalance[0], imbalance[2]])))

    print(f"  Summary: max|Δφ|={max_dphi:.3f} deg  min_CMRR={min_cmrr:.3f} dBe  "
          f"max_excess_loss={max_excess_loss:.3f} dB  max_imbalance={max_imbalance:.3f} dB")

    out_path = os.path.join(output_dir, f"offset_{offset:.4f}.npz")
    np.savez(out_path, offset=offset, wavelengths=wavelengths,
             imbalance_in1=imbalance[0], imbalance_in3=imbalance[2],
             excess_loss_in1=excess_loss[0], excess_loss_in3=excess_loss[2],
             cmrr_I_in1=cmrr_I[0], cmrr_Q_in1=cmrr_Q[0],
             cmrr_I_in3=cmrr_I[2], cmrr_Q_in3=cmrr_Q[2],
             dphi_6=dphi[1], dphi_7=dphi[2], dphi_8=dphi[3],
             max_dphi=max_dphi, min_cmrr=min_cmrr,
             max_excess_loss=max_excess_loss, max_imbalance=max_imbalance)
    print(f"  saved -> {out_path}")
    return dict(offset=offset, max_dphi=max_dphi, min_cmrr=min_cmrr,
                max_excess_loss=max_excess_loss, max_imbalance=max_imbalance)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        offsets = [float(sys.argv[1])]
    else:
        offsets = [BASELINE_OFFSET]

    results = [solve_offset_point(off) for off in offsets]

    print("\n=== Sweep summary ===")
    for r in results:
        print(f"  offset={r['offset']:.4f}: max|Δφ|={r['max_dphi']:.3f} deg  "
              f"min_CMRR={r['min_cmrr']:.3f} dBe  max_XL={r['max_excess_loss']:.3f} dB  "
              f"max_imbalance={r['max_imbalance']:.3f} dB")
