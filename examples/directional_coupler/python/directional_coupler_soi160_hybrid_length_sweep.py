"""Tutorial: directional_coupler_hybrid.py's coupling-length sweep method
(same [FDTD bend | EME body | FDTD bend] decomposition, cost independent of
how many lengths are swept -- see that script's own module docstring point
"Coupling-length sweep") re-derived on SOI_160 at duty/dc.gds's real
dimensions, sweeping coupling_length 0.5-30 um.

WHY THIS EXISTS (2026-08-10): directional_coupler_soi160_hybrid.py only ever
solves the die's ONE real fabricated coupling_length (3.0 um, fixed, not
swept -- its own docstring point 1). This script answers a different
question -- how does the coupling ratio move across a much wider length
range on this SAME cross-section -- for the cost of ONE hybrid run, not
a run per length.

CACHE REUSE, VERIFIED (not assumed) against example_sim_cache/ before this
was run: the two FDTD bend-section solves and the FDTD-vs-MPB body-beta
calibration are IDENTICAL in every argument (geometry dims, BODY_STUB=0.25,
boundaries, wavelengths, resolutions, material_stack) to directional_
coupler_soi160_hybrid.py's own already-completed COUPLING_LENGTH=3.0
production run (job 11716260) -- same absolute x0/x1 (confirmed x0=-11.5,
matching that run's own cached FDTDSolver_solve_section entries). Since the
FDTD section cache keys on absolute position (see directional_coupler_
soi160_hybrid_dc_test.py's own CACHE NOTE) and DirectionalCoupler geometry
is built symmetric about the coupling-region centre, this script
deliberately builds its OWN reference `dc` at coupling_length=3.0 -- not the
"middle of the sweep" convention directional_coupler_hybrid.py's SOI_220
version uses -- specifically so M_in0/M_in1/_beta are cache HITS, not new
FDTD solves. Only the swept EME body (body_solver.solve_section_length_sweep,
MPB-only, no FDTD) pays for new lengths; it is cheap regardless of grid size
(one supermode solve per wavelength, reused via the analytic per-length
phase matrix S_ii(L) = exp(2j*pi*neff_i/wavelength*L) -- see EMESolver.
solve_section_length_sweep()'s own docstring in beef/solvers/eme.py).

WHY THE SWEEP FLOOR IS 0.5 um, NOT A LITERAL 0: identical reasoning to
directional_coupler_hybrid.py's own floor at 2*BODY_STUB -- each FDTD
section already folds BODY_STUB=0.25 um of the straight coupling region
into itself (past the S-bend/coupling-region seam, for a clean interface
monitor on a uniform cross-section rather than right on the geometric
discontinuity), so a device shorter than both stubs combined (0.5 um)
cannot be split into the three sections at all -- its EME body length would
be negative. A true coupling_length=0 needs the full single-domain
FDTD method instead (skips the EME body section entirely) --
directional_coupler_soi160_fdtd_dc_test.py 0.0, run separately.

The 4 real 4C10-die "DC_test" device lengths characterized elsewhere
(directional_coupler_soi160_hybrid.py at 3.0, directional_coupler_soi160_
hybrid_dc_test.py at 11.634/0.5/[0.0 asserts]) are marked on the sweep plot
for direct comparison against this continuous curve.

Run:
    python /path/to/directional_coupler_soi160_hybrid_length_sweep.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import (
    DirectionalCoupler, EMESolver, FDTDSolver, HybridSolver, SOI_160, plot_hybrid_geometry,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "directional_coupler",
    "directional_coupler_soi160_hybrid_length_sweep"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters -- identical to directional_coupler_soi160_
# hybrid.py's production values, so the reference-geometry solves below hit
# that script's own cache (see module docstring's CACHE REUSE section).
# ---------------------------------------------------------------------------
FDTD_RESOLUTION = 40.0
EME_RESOLUTION = 80.0
WL_MIN = 1.50
WL_MAX = 1.60
N_WL = 11

BODY_STUB = 0.25   # um -- must match directional_coupler_soi160_hybrid.py exactly for the cache hit

N_MODES = 2   # measured simplification -- see directional_coupler_soi160_hybrid.py's own docstring point 2

CALIBRATE_BODY = True

# Reference length for the FDTD bend sections -- deliberately the die's own
# already-simulated 3.0 um, NOT the middle of the sweep grid (contrast
# directional_coupler_hybrid.py's L_mid convention) -- see module docstring.
L_REF = 3.0   # um

# The actual swept range: 0.5 (= 2*BODY_STUB, the shortest length the 3-section
# decomposition supports) to 30 um, requested range.
SWEEP_MIN = 2 * BODY_STUB   # 0.5 um
SWEEP_MAX = 30.0
N_SWEEP = 200

# The 4 real 4C10-die device lengths, for the comparison markers on the sweep
# plot -- see directional_coupler_soi160_hybrid_dc_test.py's own docstring
# for the measurement-based id-to-length derivation. 0.0 is included in the
# marker list (labeled, not plotted as a curve point -- see SWEEP_MIN above).
DEVICE_LENGTHS = [
    (11.634, "1901/1902"),
    (3.0, "1911/1912"),
    (0.5, "1921/1922"),
    (0.0, "1931/1932 (unreachable via hybrid; see fdtd_dc_test script)"),
]

coupling_lengths = np.linspace(SWEEP_MIN, SWEEP_MAX, N_SWEEP)
for _L, _ in DEVICE_LENGTHS:
    if _L >= SWEEP_MIN and not np.any(np.isclose(coupling_lengths, _L, rtol=0, atol=1e-9)):
        coupling_lengths = np.sort(np.append(coupling_lengths, _L))

# ---------------------------------------------------------------------------
# 1. Geometry -- built at L_REF, the die's own real dimensions.
# ---------------------------------------------------------------------------
dc = DirectionalCoupler(
    waveguide_width=0.5,
    gap=0.25,
    coupling_length=L_REF,
    bend_length=10.0,
    waveguide_offset=4.0,
    bend_type="euler",
    waveguide_height=0.16,
)
print(f"Directional coupler (SOI_160, duty/dc.gds dims) reference total length: "
      f"{dc.total_length:.2f} um (coupling pitch {dc.coupling_pitch:.3f} um, "
      f"far pitch {dc.waveguide_offset:.2f} um)")
print(f"Coupling length sweep: {coupling_lengths[0]:.3f} - {coupling_lengths[-1]:.1f} um "
      f"({len(coupling_lengths)} pts, reference build length = {L_REF:.2f} um)")

# ---------------------------------------------------------------------------
# 2. Configure sub-solvers -- one per section, params identical to
#    directional_coupler_soi160_hybrid.py's own.
# ---------------------------------------------------------------------------
_fdtd_common = dict(
    res=FDTD_RESOLUTION,
    sim_width=7.0,
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    half_domain_aperture=True,
    material_stack=SOI_160,
    decay_tol=1e-11,
)

fdtd_solver = FDTDSolver(**_fdtd_common)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholder; solve_section_length_sweep() sweeps its own
    oxide_width=7.0, oxide_height=1.5, res=EME_RESOLUTION, material_stack=SOI_160,
    mode_solver="meep",
)

# ---------------------------------------------------------------------------
# 3. Build the reference section boundaries -- FDTD bend | EME body | FDTD bend,
#    at L_REF. Only the body span [boundaries[1], boundaries[2]] changes as the
#    sweep below re-evaluates the body's phase matrix at each swept length; the
#    two FDTD sections are solved once here, at L_REF, and never resolved.
# ---------------------------------------------------------------------------
half_body = dc.coupling_length / 2
half_total = dc.total_length / 2

boundaries = [
    -half_total,
    -half_body + BODY_STUB,
    +half_body - BODY_STUB,
    +half_total,
]
interface_modes = [N_MODES, N_MODES]

body_length_ref = boundaries[2] - boundaries[1]
assert body_length_ref > 0
print(f"Reference coupling-length budget: 2 x {BODY_STUB:.2f} um inside the FDTD "
      f"sections + EME body {body_length_ref:.2f} um = {L_REF:.2f} um total")

body_lengths = coupling_lengths - 2 * BODY_STUB
assert body_lengths[0] >= -1e-9, (
    f"coupling_lengths must start at >= 2*BODY_STUB ({2 * BODY_STUB:.2f} um); "
    f"got {coupling_lengths[0]:.3f} um"
)
body_lengths = np.clip(body_lengths, 0.0, None)

_iface_mon = dict(monitor_center=0.0, monitor_size=fdtd_solver.sim_width)

# ---------------------------------------------------------------------------
# 4. Geometry plot (at L_REF).
# ---------------------------------------------------------------------------
hybrid = HybridSolver(
    solvers=[fdtd_solver, body_solver, fdtd_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
)
plot_hybrid_geometry(
    dc, solver=hybrid,
    in_modes=1, out_modes=1,
    title=f"Directional Coupler (SOI_160, coupling_length={L_REF:.1f} um reference) — "
          f"Hybrid Solver Sections",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True,
    sub_region_show_eme_detail=True,
    sub_region_configs=[
        ("in", {0: {"input_idx": 0, **_iface_mon}}, {2: 0}),
        ("isolated", {0: {"input_idx": 1, **_iface_mon}}, {2: 0}),
    ],
)
print(f"Geometry plot saved -> {output_dir}/geometry_hybrid.png")
print(f"Sub-region plots saved -> {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run -- 2 FDTD solves (input section, once per input arm; cache HIT
#    against directional_coupler_soi160_hybrid.py's own production run, see
#    module docstring) + one EME mode-solve per wavelength, swept over
#    coupling_lengths via the length axis broadcast (identical contraction to
#    directional_coupler_hybrid.py step 5).
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nWavelength sweep: {wavelengths[0]:.3f} - {wavelengths[-1]:.3f} um "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

M_in0 = fdtd_solver.solve_section(
    dc, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
    input_idx=0, save_dir=save_dir, **_iface_mon,
)  # [n_wl, N_MODES, 1] -- "in" arm -> body modes
M_in1 = fdtd_solver.solve_section(
    dc, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=fdtd_solver.waveguide_length, stub_right=0.0,
    input_idx=1, save_dir=save_dir, **_iface_mon,
)  # [n_wl, N_MODES, 1] -- "isolated" arm -> body modes

M_out_through = M_in0.transpose(0, 2, 1)  # [n_wl, 1, N_MODES] -- mirror-symmetry reuse
M_out_cross = M_in1.transpose(0, 2, 1)    # [n_wl, 1, N_MODES]

_body_neff = None
if CALIBRATE_BODY:
    _beta = fdtd_solver.section_mode_beta(
        dc, wavelengths, boundaries[1], boundaries[2], N_MODES,
        save_dir=save_dir, **_iface_mon,
    )  # [n_wl, N_MODES] -- cache HIT, same body span as the L_REF=3.0 production run
    _body_neff = _beta * wavelengths[:, None] / (2 * np.pi)
    _i = int(np.argmin(np.abs(wavelengths - 1.55)))
    _seed_half = min(1.0, half_body - 2.0 * fdtd_solver.dpml)
    _mpb = fdtd_solver.section_mode_neff(
        dc, wavelengths[_i:_i + 1], -_seed_half, +_seed_half, N_MODES, **_iface_mon)[0]
    print(f"Body beta calibration at lambda = {wavelengths[_i]:.3f} um (even/odd supermodes):")
    print(f"    MPB mode index  n_eff {_mpb[0]:.6f} {_mpb[1]:.6f}   "
          f"delta_n_eff {_mpb[0] - _mpb[1]:.6f}")
    print(f"    FDTD own phase  n_eff {_body_neff[_i, 0]:.6f} {_body_neff[_i, 1]:.6f}   "
          f"delta_n_eff {_body_neff[_i, 0] - _body_neff[_i, 1]:.6f}   <- body propagates on this")

M_body = body_solver.solve_section_length_sweep(
    dc, wavelengths, boundaries[1], boundaries[2], body_lengths,
    n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
    neff_override=_body_neff,
)  # [n_L, n_wl, N_MODES, N_MODES]

tmp = np.einsum("Lwij,wjk->Lwik", M_body, M_in0)                  # [n_L, n_wl, N_MODES, 1]
full_through = np.einsum("wij,Lwjk->Lwik", M_out_through, tmp)    # [n_L, n_wl, 1, 1]
full_cross = np.einsum("wij,Lwjk->Lwik", M_out_cross, tmp)        # [n_L, n_wl, 1, 1]

T_through_grid = np.abs(full_through[:, :, 0, 0]) ** 2   # [n_L, n_wl]
T_cross_grid = np.abs(full_cross[:, :, 0, 0]) ** 2
T_total_grid = T_through_grid + T_cross_grid

ratio_mean = T_cross_grid.mean(axis=1) / np.clip(T_total_grid.mean(axis=1), 1e-20, None)
sign_changes = np.flatnonzero(np.diff(np.sign(ratio_mean - 0.5)))
if sign_changes.size:
    _a, _b = sign_changes[0], sign_changes[0] + 1
    _r0, _r1 = ratio_mean[_a], ratio_mean[_b]
    L50 = float(coupling_lengths[_a] + (0.5 - _r0) * (coupling_lengths[_b] - coupling_lengths[_a]) / (_r1 - _r0)) \
        if _r1 != _r0 else float(coupling_lengths[_a])
    print(f"Measured 50:50 length in [{SWEEP_MIN:.1f}, {SWEEP_MAX:.1f}] um: {L50:.4f} um "
          f"(ratio {ratio_mean[_a]:.4f} -> {ratio_mean[_b]:.4f})")
else:
    print(f"No 0.5 (50:50) crossing found in [{SWEEP_MIN:.1f}, {SWEEP_MAX:.1f}] um "
          f"(ratio_mean range: {ratio_mean.min():.4f}-{ratio_mean.max():.4f})")

# ---------------------------------------------------------------------------
# 6. Save and plot.
# ---------------------------------------------------------------------------
np.savez(
    os.path.join(output_dir, "length_sweep.npz"),
    coupling_lengths=coupling_lengths, wavelengths=wavelengths,
    T_through=T_through_grid, T_cross=T_cross_grid, T_total=T_total_grid,
    ratio_mean=ratio_mean, L_ref=L_REF, body_stub=BODY_STUB,
)
print(f"Length-sweep grid saved -> {output_dir}/length_sweep.npz")

fig, ax = plt.subplots(figsize=(8, 5))

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

for _L, _name in DEVICE_LENGTHS:
    ax.axvline(_L, color="gray", linestyle=":", lw=0.8)
    ax.text(_L, 1.02, _name.split(" ")[0], rotation=90, fontsize=6,
            ha="right", va="bottom", color="gray")

ax.set_xlabel("Coupling length (um)")
ax.set_ylabel("Transmission (mean over 1.50-1.60 um, shaded = min-max)")
ax.set_title("Directional Coupler (SOI_160, duty/dc.gds dims) — Power Transfer vs Coupling Length")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.08)
ax.set_xlim(0, SWEEP_MAX)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "coupling_length_sweep.png"), dpi=300)
plt.close("all")
print(f"Coupling-length sweep plot saved -> {output_dir}/coupling_length_sweep.png")
