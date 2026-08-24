"""
Tutorial: 2×2 MMI — FDTD/EME hybrid simulation driven by a GDS FILE.

This is `mmi_2x2_hybrid.py` with exactly one thing changed: the geometry
backend. Where that example hands `HybridSolver` a parametric `MMI(...)`,
this one hands it a `GdsGeometry` reading `examples/gds_files/
mmi_2x2_twin.gds` — a file whose core-layer polygons are geometrically
identical to what that `MMI(...)` draws. Same resolution, same material
stack, same C-band sweep, same 14 interface modes, same boundaries, same
`body_stub`, same mirror/mode-reuse arguments.

That makes it a controlled comparison rather than a second device: any
disagreement between the two runs is attributable to the GDS path — to
polygon clipping, branch selection, or the device-frame shift — and not to
physics. It is the acceptance vehicle for
`validation_markdown/resolved/gds_geometry_hybrid_cascade_validation.md`, whose step
1 specifies this file and whose step 5 is the comparison below.

What is actually new in the library
-----------------------------------
`GdsGeometry` could already back a whole-device single-solver run (the HCG
reflector example does that). It could NOT occupy a section of a cascade:
`section_prisms`, `section_gds` and `width_at` were inherited unimplemented
from `DeviceGeometry`, and `center_at` silently returned 0.0 — wrong for
any off-axis branch, which is every port of a 2×2. Those are implemented as
of 2026-08-20 by intersecting the imported polygons with a vertical line:

  * `width_at` / `center_at` — the y-intervals a vertical line at x cuts,
    merged across layers, indexed by `input_idx`/`output_idx` exactly as
    `MMI`'s are, so the `geometry_kwargs` below are unchanged from the
    parametric example's.
  * `section_prisms` / `section_gds` — those polygons clipped to the
    section's span with gdstk booleans.
  * cut validity — `MMI` knows where its taper bands are and refuses a cut
    inside one; an imported layout has no bands to consult, so the same
    refusal is *measured*, by sampling `width_at`/`center_at` around the cut
    (see `GdsGeometry._check_cut`, including what sampling cannot see).
  * `access_extension=True` — a cascade's end sections always ask for a
    stub past the physical port, and the twin file draws the bare device
    (no access waveguides), so the stub is synthesized by extruding the end
    facet's own measured cross-section. This is the one place the GDS path
    adds geometry the file does not contain, and step 0 below checks it
    against what `MMI` synthesizes for the same span.

Why the 2×2 rather than something simpler: it has close-set ports (1.5 µm
centre-to-centre, which is why both FDTD sections use
`half_domain_aperture`) and off-axis branches, which is precisely what
makes `center_at`'s old 0.0 default insufficient and what stresses branch
selection hardest.

The physics, the section chain, the mirror-symmetry reuse, the analytic
multimode-length sweep and the MODE_SIGNS derivation are all identical to
`mmi_2x2_hybrid.py` — read that file's docstring for them; they are not
re-explained here.

Run:
    python /path/to/examples/mmi_2x2/python/mmi_2x2_hybrid_gds.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import (
    GdsGeometry, MMI, EMESolver, FDTDSolver, HybridSolver, SMatrix, SOI_220,
    plot_beta_calibration, plot_hybrid_geometry,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_2x2", "mmi_2x2_hybrid_gds"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))
_GDS_PATH = os.path.abspath(os.path.join(_HERE, "..", "..", "gds_files", "mmi_2x2_twin.gds"))
_PARAMETRIC_OUT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "mmi_2x2", "mmi_2x2_hybrid"))

# ---------------------------------------------------------------------------
# Global simulation parameters — every one of these is verbatim from
# mmi_2x2_hybrid.py. They must stay that way: the comparison in step 7 is
# only meaningful if the geometry backend is the sole difference.
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by both the EME and FDTD sub-solvers
WL_MIN = 1.530      # wavelength sweep start (µm)
WL_MAX = 1.565      # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

CALIBRATE_BODY = True    # propagate the EME body on the FDTD's own measured beta
SYMMETRIC_MODE_REUSE = True   # derive M_in0 from M_in1 via MODE_SIGNS instead of a 2nd FDTD run

INPUT_METHOD = "fdtd"    # "fdtd" or "eme"
BODY_METHOD = "eme"      # "fdtd" or "eme"
OUTPUT_METHOD = "fdtd"   # "fdtd" or "eme"

# ---------------------------------------------------------------------------
# 1. Geometry — read from the GDS file, not built parametrically.
#
#    region_lengths is the file's own [in-taper, body, out-taper] breakdown;
#    it only feeds eme_lengths() (unused in a hybrid cascade) and is checked
#    at construction against the imported bbox span, so a stale entry here
#    fails loudly rather than quietly mis-describing the file.
#
#    access_extension=True: the twin holds the bare device, x in [-18, +18],
#    with no access waveguides. See the module docstring.
# ---------------------------------------------------------------------------
TAPER_LENGTH = 4.0
MM_LENGTH = 28.0

gds = GdsGeometry(
    _GDS_PATH,
    core_layers=[(2, 0)],
    region_lengths=[TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH],
    cell_name="mmi_2x2_twin",
    access_extension=True,
)
print(f"GDS geometry: {_GDS_PATH}")
print(f"  total length : {gds.total_length:.2f} µm (bbox {gds.bbox[0]:.2f} x {gds.bbox[1]:.2f})")
print(f"  input facet  : branches at y = "
      f"{[round((a + b) / 2, 4) for a, b in gds._facet_left]} µm, "
      f"widths {[round(b - a, 4) for a, b in gds._facet_left]} µm")
print(f"  output facet : branches at y = "
      f"{[round((a + b) / 2, 4) for a, b in gds._facet_right]} µm")

# The parametric twin — used ONLY by step 0's geometry cross-check and by
# step 7's results comparison, never by the simulation itself.
mmi_reference = MMI(
    n_inputs=2, n_outputs=2, mm_length=MM_LENGTH, mm_width=2.8,
    waveguide_width=0.45, taper_length=TAPER_LENGTH, taper_width=1.3,
    waveguide_offset=1.5, waveguide_height=0.22,
)

# ---------------------------------------------------------------------------
# 2. Sub-solvers — identical to mmi_2x2_hybrid.py's.
# ---------------------------------------------------------------------------
_fdtd_common = dict(
    res=RESOLUTION,
    sim_width=5.0,
    sim_height=1.5,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    material_stack=SOI_220,
    decay_tol=1e-11,
)

N_MODES = 14   # see mmi_2x2_hybrid.py: band 14 is a radiation mode, dropped.

# Measured against the parametric 2×2's own 2-FDTD-arm run. A fixed
# geometric property of the mirror-symmetric cross-section, and the two
# geometries ARE the same cross-section (step 0 proves it), so the same
# vector applies here. Re-derive with SYMMETRIC_MODE_REUSE = False.
MODE_SIGNS = np.array([1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, 1])


def make_section_solver(method: str, half_domain_aperture: bool = False):
    """Identical to mmi_2x2_hybrid.py's — see there."""
    if method == "fdtd":
        kwargs = dict(_fdtd_common)
        kwargs["half_domain_aperture"] = half_domain_aperture
        return FDTDSolver(**kwargs)
    if method == "eme":
        return EMESolver(
            num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
            oxide_width=5.0, oxide_height=1.5, res=RESOLUTION, material_stack=SOI_220,
            mode_solver="meep",
        )
    raise ValueError(f"unknown method {method!r}, expected 'fdtd' or 'eme'")


input_solver = make_section_solver(INPUT_METHOD, half_domain_aperture=True)
body_solver = make_section_solver(BODY_METHOD)
output_solver = make_section_solver(OUTPUT_METHOD, half_domain_aperture=True)

# ---------------------------------------------------------------------------
# 3. Section boundaries.
#
#    These are numerically identical to mmi_2x2_hybrid.py's, and that is not
#    a coincidence to be maintained by hand: GdsGeometry centres the imported
#    bbox on x=0, and MMI's frame origin is its MM body centre. For a
#    symmetric 2×2 (equal input/output taper lengths) MMI.frame_offset == 0,
#    so those two points coincide and the device frames are the same frame.
#    Asserted below rather than trusted.
# ---------------------------------------------------------------------------
half_body = MM_LENGTH / 2
half_total = gds.total_length / 2
body_stub = 0.25   # µm of MMI body folded into each FDTD section

assert abs(mmi_reference.frame_offset) < 1e-12, (
    "MMI.frame_offset != 0 — the parametric and GDS device frames no longer "
    "coincide, so these boundaries do not transfer between the two runs."
)
assert abs(gds.total_length - mmi_reference.total_length) < 1e-9, (
    f"twin span {gds.total_length} != MMI.total_length "
    f"{mmi_reference.total_length}; regenerate with "
    f"examples/gds_files/make_mmi_2x2_twin.py"
)

boundaries = [
    -half_total,               # physical input port (input taper tip)
    -half_body + body_stub,    # input section <-> body section interface
    +half_body - body_stub,    # body section <-> output section interface
    +half_total,               # physical output port (output taper tip)
]
interface_modes = [N_MODES, N_MODES]

mm_lengths = np.linspace(16.0, 40.0, 201)   # µm; nominal 28.0 lands on index 100
body_lengths = mm_lengths - 2 * body_stub

hybrid = HybridSolver(
    solvers=[input_solver, body_solver, output_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
)

# ---------------------------------------------------------------------------
# 0. Geometry cross-check — the "cheap first" half of the acceptance test.
#
#    Run BEFORE any solver: if the two backends hand the solvers the same
#    polygons and the same cross-section queries, the physics must agree, and
#    a failure here localizes the bug to geometry rather than leaving it to
#    be inferred from a transmission mismatch hours later. Numbered 0 because
#    it costs nothing and gates everything after it.
# ---------------------------------------------------------------------------
print("\n--- step 0: GDS geometry vs. parametric MMI (no simulation) ---")

_WG = 0.5 + _fdtd_common["dpml"]   # solve_section draws stub_left + dpml


def _area(polys):
    total = 0.0
    for poly in polys:
        s = 0.0
        for i in range(len(poly)):
            x0, y0 = poly[i]
            x1, y1 = poly[(i + 1) % len(poly)]
            s += x0 * y1 - x1 * y0
        total += abs(s) / 2.0
    return total


_xs = np.linspace(boundaries[0] - _WG, boundaries[-1] + _WG, 2001)
for _kw in ({"input_idx": 1, "output_idx": 1}, {"input_idx": 0, "output_idx": 0}):
    _dw = max(abs(mmi_reference.width_at(x, **_kw) - gds.width_at(x, **_kw)) for x in _xs)
    _dc = max(abs(mmi_reference.center_at(x, **_kw) - gds.center_at(x, **_kw)) for x in _xs)
    print(f"  {_kw}: max |Δwidth_at| = {_dw:.3e} µm, max |Δcenter_at| = {_dc:.3e} µm")
    assert _dw < 1e-9 and _dc < 1e-9, "GDS cross-section queries disagree with MMI's"

_sections = [
    ("input ", 0, _WG, 0.0, {"input_idx": 1}),
    ("body  ", 1, 0.0, 0.0, {}),
    ("output", 2, 0.0, _WG, {"output_idx": 1}),
]
for _name, _k, _sl, _sr, _kw in _sections:
    _x0, _x1 = boundaries[_k], boundaries[_k + 1]
    _a = mmi_reference.section_gds(_x0, _x1, 6.0, stub_left=_sl, stub_right=_sr, **_kw)
    _b = gds.section_gds(_x0, _x1, 6.0, stub_left=_sl, stub_right=_sr, **_kw)
    _da = abs(_area(_a[(2, 0)]) - _area(_b[(2, 0)]))
    print(f"  section {_name} [{_x0 - _sl:7.3f},{_x1 + _sr:7.3f}]: "
          f"{len(_a[(2, 0)])} vs {len(_b[(2, 0)])} core polygons, "
          f"|Δarea| = {_da:.3e} µm², oxide identical: {_a[(1, 0)] == _b[(1, 0)]}")
    assert _da < 1e-9, f"section {_name}: core polygon area disagrees with MMI's"
    assert _a[(1, 0)] == _b[(1, 0)], f"section {_name}: oxide box disagrees with MMI's"

# The trap step 5 of the validation doc calls out as most likely to produce
# a false pass: if GdsGeometry's cache fingerprint reduced to something the
# MMI run also produces, this run would silently serve the PARAMETRIC run's
# cached results and report a perfect, meaningless match. It does not —
# _fingerprint() carries geometry_class alongside the section polygons — but
# assert it here rather than trusting it, since the failure is invisible.
from beef.cache.config import build_run_config          # noqa: E402
from beef.cache.store import _hash_config               # noqa: E402

_probe = dict(x0=boundaries[0], x1=boundaries[1],
              wavelengths=[WL_MIN, WL_MAX], n_modes_in=1, n_modes_out=N_MODES,
              stub_left=input_solver.waveguide_length, stub_right=0.0,
              geometry_kwargs={"input_idx": 1})
_h_mmi = _hash_config(build_run_config(input_solver, "solve_section", _probe, mmi_reference))
_h_gds = _hash_config(build_run_config(input_solver, "solve_section", _probe, gds))
print(f"  cache fingerprint: MMI {_h_mmi} vs GDS {_h_gds} — "
      f"{'distinct' if _h_mmi != _h_gds else 'COLLIDING'}")
assert _h_mmi != _h_gds, (
    "the GDS and parametric runs share a cache fingerprint — this run would "
    "serve the parametric run's results and the comparison would be vacuous"
)
print("  step 0 passed: the two backends are geometrically indistinguishable.")

# ---------------------------------------------------------------------------
# 4. Geometry plot — same section map as the parametric example's, drawn
#    from the imported polygons (the access waveguides visible outside
#    x = ±18 are access_extension's synthesized stubs, not file content).
# ---------------------------------------------------------------------------
plot_hybrid_geometry(
    gds, solver=hybrid,
    in_modes=1, out_modes=1,
    title="2×2 MMI (from GDS) — Hybrid Solver Sections",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True,
    sub_region_show_eme_detail=True,
    sub_region_configs=[
        ("in_top", {0: {"input_idx": 1}}, {2: 0}),
        ("in_bot", {0: {"input_idx": 0}}, {2: 0}),
    ],
)
print(f"\nGeometry plot saved → {output_dir}/geometry_hybrid.png")
print(f"Sub-region plots saved → {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run the hybrid pass — byte-for-byte the parametric example's step 5,
#    with `gds` where it passes `mmi`.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nHybrid sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, {N_MODES} interface modes), "
      f"mm_length {mm_lengths[0]:.1f} – {mm_lengths[-1]:.1f} µm "
      f"({len(mm_lengths)} pts)\n")

M_in1 = input_solver.solve_section(
    gds, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=input_solver.waveguide_length, stub_right=0.0,
    input_idx=1, save_dir=save_dir,
)  # [n_wl, N_MODES, 1] — excited (top) arm -> body modes

if SYMMETRIC_MODE_REUSE:
    M_in0 = M_in1 * MODE_SIGNS[None, :, None]
else:
    M_in0 = input_solver.solve_section(
        gds, wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=input_solver.waveguide_length, stub_right=0.0,
        input_idx=0, save_dir=save_dir,
    )

M_out_bar = M_in1.transpose(0, 2, 1)
M_out_cross = M_in0.transpose(0, 2, 1)

body_neff = None
if CALIBRATE_BODY:
    calibration = input_solver.section_mode_beta_details(
        gds, wavelengths, boundaries[1], boundaries[2], N_MODES, save_dir=save_dir,
    )
    beta = calibration.beta  # [n_wl, N_MODES] rad/µm
    body_neff = beta * wavelengths[:, None] / (2 * np.pi)

    plot_beta_calibration(
        calibration, gds, input_solver,
        os.path.join(output_dir, "calibration"),
        monitor_center=0.0, monitor_size=input_solver.sim_width,
    )
    print(f"Calibration plots saved → {output_dir}/calibration/ (geometry/ + fits/)")

M_body = body_solver.solve_section_length_sweep(
    gds, wavelengths, boundaries[1], boundaries[2], body_lengths,
    n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
    neff_override=body_neff,
)  # [n_L, n_wl, N_MODES, N_MODES]

tmp = np.einsum("Lwij,wjk->Lwik", M_body, M_in1)
coeffs_bar = np.einsum("wij,Lwjk->Lwik", M_out_bar, tmp)
coeffs_cross = np.einsum("wij,Lwjk->Lwik", M_out_cross, tmp)

T_bar_grid   = np.abs(coeffs_bar[:, :, 0, 0]) ** 2    # [n_L, n_wl]
T_cross_grid = np.abs(coeffs_cross[:, :, 0, 0]) ** 2
T_total_grid = T_bar_grid + T_cross_grid
ratio_grid   = T_bar_grid / np.clip(T_total_grid, 1e-20, None)

L_plot_idx = int(np.argmin(np.abs(mm_lengths - MM_LENGTH)))
L_plot = mm_lengths[L_plot_idx]

T_bar   = T_bar_grid[L_plot_idx]
T_cross = T_cross_grid[L_plot_idx]
T_total = T_total_grid[L_plot_idx]

wl_c = wavelengths[len(wavelengths) // 2]
idx  = np.argmin(np.abs(wavelengths - wl_c))

print(f"\nResults at λ = {wl_c:.3f} µm, mm_length = {L_plot:.1f} µm:")
print(f"  Bar output   T = {T_bar[idx]:.4f}  IL = {10*np.log10(max(T_bar[idx],1e-20)):.2f} dB")
print(f"  Cross output T = {T_cross[idx]:.4f}  IL = {10*np.log10(max(T_cross[idx],1e-20)):.2f} dB")
print(f"  Total        T = {T_total[idx]:.4f}  IL = {10*np.log10(max(T_total[idx],1e-20)):.2f} dB")
print(f"  (ideal 50/50 split: T_each=0.5, IL_total=0.0 dB)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat = SMatrix(wavelengths=wavelengths)
smat.set(coeffs_bar[L_plot_idx, :, 0, 0], port="bar")
smat.set(coeffs_cross[L_plot_idx, :, 0, 0], port="cross")
smat.check_passive()
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix (mm_length={L_plot:.1f} µm slice) saved → {output_dir}/smatrix_hybrid.pkl")

np.savez(os.path.join(output_dir, "mm_length_sweep.npz"),
         mm_lengths=mm_lengths, wavelengths=wavelengths,
         T_bar=T_bar_grid, T_cross=T_cross_grid)
print(f"Length-sweep grids saved → {output_dir}/mm_length_sweep.npz")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(wavelengths, T_bar,   label="Bar output",   color="steelblue")
ax.plot(wavelengths, T_cross, label="Cross output", color="darkorange", linestyle="--")
ax.plot(wavelengths, T_total, label="Total",         color="k", lw=1.5)
ax.axhline(0.5, color="gray", linestyle=":", lw=0.8, label="Ideal per-port")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title(f"Hybrid (GDS) — Transmission vs Wavelength (mm_length={L_plot:.1f} µm)")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)

ax = axes[1]
IL_bar   = 10 * np.log10(np.clip(T_bar,   1e-20, None))
IL_cross = 10 * np.log10(np.clip(T_cross, 1e-20, None))
IL_tot   = 10 * np.log10(np.clip(T_total, 1e-20, None))
ax.plot(wavelengths, IL_bar,   label="Bar output",   color="steelblue")
ax.plot(wavelengths, IL_cross, label="Cross output", color="darkorange", linestyle="--")
ax.plot(wavelengths, IL_tot,   label="Total",         color="k", lw=1.5)
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title(f"Hybrid (GDS) — Insertion Loss vs Wavelength (mm_length={L_plot:.1f} µm)")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle("2×2 MMI from GDS — FDTD/EME Hybrid Simulation", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/hybrid_results.png")

# ---------------------------------------------------------------------------
# 7. Acceptance comparison against the parametric run, when its output is
#    present (run examples/mmi_2x2/python/mmi_2x2_hybrid.py first).
#
#    MEASURED 2026-08-21 (job 12069544, RESOLUTION=40, N_WL=11, 201 lengths):
#
#        T_bar   max |Δ| = 6.556e-04  (rel 1.603e-03)
#        T_cross max |Δ| = 5.218e-04  (rel 1.346e-03)
#
#    That is LOOSER than "solver noise", and the reason is known and is NOT a
#    bug in the GDS path: the two backends hand meep the same five polygons
#    for the input section in a DIFFERENT ORDER, and meep's subpixel averaging
#    is order-sensitive where prisms share an edge. Reordering GDS's polygons
#    to MMI's makes the discretized permittivity BIT-IDENTICAL (0 of 2.52M
#    voxels differ; 872 differ natively, all at the stub/taper and taper/body
#    seams). The geometry layer itself agrees to 1e-16 (step 0 above).
#
#    So: a few e-3 here is EXPECTED. What would be a real defect is a
#    percent-level difference, or any disagreement in step 0's polygon/query
#    checks — those are asserted and will fail loudly. Full write-up, including
#    the larger library-wide seam artifact this exposed (merging abutting
#    prisms moves 41,677 voxels, 48x the ordering effect, on EVERY device):
#    validation_markdown/resolved/gds_geometry_hybrid_cascade_validation.md.
# ---------------------------------------------------------------------------
_ref_npz = os.path.join(_PARAMETRIC_OUT, "mm_length_sweep.npz")
print("\n--- step 7: acceptance vs. the parametric run ---")
if not os.path.exists(_ref_npz):
    print(f"  skipped — no parametric output at {_ref_npz}")
    print("  run examples/mmi_2x2/python/mmi_2x2_hybrid.py, then re-run this.")
else:
    _ref = np.load(_ref_npz)
    if (_ref["mm_lengths"].shape != mm_lengths.shape
            or not np.allclose(_ref["mm_lengths"], mm_lengths)
            or not np.allclose(_ref["wavelengths"], wavelengths)):
        print("  skipped — the parametric output is on a different sweep grid; "
              "re-run it with this file's WL/mm_length settings.")
    else:
        for _label, _mine, _theirs in (("T_bar", T_bar_grid, _ref["T_bar"]),
                                       ("T_cross", T_cross_grid, _ref["T_cross"])):
            _abs = np.abs(_mine - _theirs)
            _rel = _abs / np.clip(np.abs(_theirs), 1e-12, None)
            _i = np.unravel_index(_abs.argmax(), _abs.shape)
            print(f"  {_label}: max |Δ| = {_abs.max():.3e} (rel {_rel.max():.3e}) "
                  f"at mm_length = {mm_lengths[_i[0]]:.2f} µm, "
                  f"λ = {wavelengths[_i[1]]:.4f} µm")
        print("  Reference (job 12069544): T_bar 6.556e-04 (rel 1.603e-03), "
              "T_cross 5.218e-04 (rel 1.346e-03).")
        print("  A few e-3 is EXPECTED — meep's subpixel averaging is sensitive "
              "to prism ORDER at abutting seams, and the two backends emit the "
              "same polygons in a different order (reordering makes the epsilon "
              "grid bit-identical). A PERCENT-level difference would be a real "
              "bug; so would any step-0 failure above.")
