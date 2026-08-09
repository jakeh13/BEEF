"""Tutorial: the HCG Fabry-Perot resonator -- two HCG mirrors terminating a
cavity, tapped in the middle by a directional coupler --
validation_markdown/high_contrast_grating_fabry_perot_validation.md, step 4.

This is the first script to actually cascade three INDEPENDENTLY solved
devices via beef.network.Circuit rather than composing sections within one
HybridSolver chain (that doc's A2) -- and the first to build a genuinely
1-port self-reflecting block (A3), since beef.network.convert.from_smatrix()
cannot represent one by its own documented design (self-terms are always
zero-filled there).

Assumes the three physical pieces already exist on disk:
  - the HCG mirror's own production spectrum
    (examples/high_contrast_grating/python/high_contrast_grating_fdtd.py's
    output, example_output/high_contrast_grating/high_contrast_grating_fdtd/
    smatrix.pkl) -- complex r(lambda), NOT just power, 21 points, 1.50-1.60 um.
  - the directional coupler re-derived on SOI_160 at duty/dc.gds's own real
    dimensions (examples/directional_coupler/python/
    directional_coupler_soi160_hybrid.py's output, smatrix_hybrid.pkl) --
    complex t(lambda)/k(lambda) ("through"/"cross"), 11 points, 1.50-1.60 um
    (an exact subset of the HCG's 21-point grid).
This script does NOT re-run either of those -- both are their own
real, expensive SLURM measurements, already landed.

Symmetric-coupler assumption, VERIFIED not just assumed (2026-08-07): a
real reciprocal, physically symmetric 4-port coupler has t'==t, k'==k
(the "isolated"-excited response equals the "in"-excited response) --
confirmed directly from this device's own cached M_in0/M_in1/M_body pieces
(not re-derived here, but checked once: k' matched k to machine precision,
t' matched t to ~1e-5). So the saved smatrix_hybrid.pkl's t/k alone are
enough to build the full 4-port block; a separate reconstruction of the
"isolated" branch is not needed.

Connecting waveguide (53.08 um each arm, Step 1's REAL pin -- traced from
duty/Mux_lpbk_9$4.gds's own real fabricated loopback structure, see
L_ARM's own comment below): a real fundamental-mode-only meow/MPB sweep
(no FDTD -- straight section), run fresh here rather than depending on a
scratch file. Built as a genuinely BIDIRECTIONAL SMatrixBlock -- NOT via
beef.network.convert.
from_transfer_matrix(), which is documented one-directional by design and
was found (2026-08-05, tests/test_fabry_perot_network_crosscheck.py) to
silently zero out a mirror's reflection when used for this. A real,
reciprocal, lossless waveguide has S[out,in] == S[in,out] == the same
phase, zero self-reflection at either end.

Two independent computations of the same spectrum, asserted to agree:
  1. The real beef.network.Circuit assembly (5 blocks, 4 connect() calls) --
     this device's actual topology.
  2. The hand-derived closed-form Fabry-Perot-with-DC-tap formula (the
     ring resonator's own closed-form-check convention, generalized to a
     4-port tap instead of one closed ring port; independently verified
     against a brute-force 4-port linear solve in
     tests/test_fabry_perot_network_crosscheck.py).
Matching to floating-point precision is the check itself -- same role
test_network_ring_crosscheck.py's cross-check plays for the ring.

KNOWN OPEN GAPS, not fixed by this script (see the validation doc):
  - FDTD_RESOLUTION=40/EME_RESOLUTION=80 (inherited from
    directional_coupler_soi160_hybrid.py) are ASSUMED, not independently
    converged for this device.
  - L_ARM=53.08 um is the AVERAGE of the real structure's own two arms
    (54.583 and 51.583 um -- a real 3.0 um asymmetry this symmetric-cavity
    model does not capture, per Step 1's own framing).

DENSE RESONANCE CHARACTERIZATION (added 2026-08-08, step 5's follow-up)

The 11-point/10 nm DC grid above is far coarser than the ~0.41 nm FWHM
estimated from this cavity's own finesse -- the sparse spectrum only shows
that a resonance exists near 1.57 um, not its real depth/width/Q/position.
The first attempt to fix this (job 11726827, a brute-force 601-point re-run
of directional_coupler_soi160_hybrid.py) TIMED OUT at 12 h after completing
only 3 of the ~7+ needed expensive FDTD sub-solves -- going 11 -> 601
wavelength points multiplies each sub-solve's own end-of-run eigenmode
decomposition cost by ~55x, not "the same run, more points for free."

Fixed instead with interpolation, not more compute: if
example_output/directional_coupler/directional_coupler_soi160_hybrid/
smatrix_hybrid_dense.pkl exists (from directional_coupler_soi160_hybrid_
dense.py -- see that script's own module docstring for the full
derivation and the real, measured per-factor validation), this script
loads it, interpolates the HCG's own r(lambda) onto the same dense grid
(already confirmed safe -- max phase step ~7.5 deg per 5 nm, nowhere near
unwrap's danger zone), fits the connecting waveguide's neff(lambda) with a
low-degree polynomial (reusing the SAME 11 real mode-solve points already
computed above -- no new MPB solves), and evaluates the closed form (NOT
re-running Circuit at thousands of points -- redundant, since Circuit and
the closed form already agree to floating-point precision on the sparse
grid above) to find the real resolved resonance lineshape.

Run:
    python /path/to/examples/high_contrast_grating_fabry_perot/python/hcg_fabry_perot_hybrid.py
"""

import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meow as mw

from beef import extrude_polygons, SOI_160
from beef.network import Circuit
from beef.network.block import SMatrixBlock
from beef.results import SMatrix

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "high_contrast_grating_fabry_perot", "hcg_fabry_perot_hybrid")
os.makedirs(output_dir, exist_ok=True)

HCG_SMATRIX = os.path.join(_ROOT, "example_output", "high_contrast_grating", "high_contrast_grating_fdtd", "smatrix.pkl")
DC_SMATRIX = os.path.join(_ROOT, "example_output", "directional_coupler", "directional_coupler_soi160_hybrid", "smatrix_hybrid.pkl")

L_ARM = 53.08   # um -- Step 1's REAL pin, 2026-08-08, traced from duty/Mux_lpbk_9$4.gds's own
                # real fabricated mirror-DC-mirror loopback structure (a real test structure with
                # this device's own topology), NOT the earlier 20 um placeholder. Path traced
                # instance-by-instance (gdstk reference transforms: origin/rotation/x_reflection),
                # left arm = straight(25.375) + bend(R=10.0, quarter-circle arc len=15.708,
                # user's instruction: treat bends as straight waveguide) + straight(13.5) =
                # 54.583 um; right arm = straight(25.375) + bend(15.708) + straight(10.5) =
                # 51.583 um. R=10.0 confirmed two independent ways: from the bend cell's own core
                # polygon extent, and from each bend's own port-to-port chord length (chord=R*sqrt(2)
                # for a 90 deg arc) -- both give exactly 10.0. The two real arms differ by exactly
                # 3.0 um (a real layout asymmetry, not modeled here -- this device's own symmetric-
                # cavity assumption, Step 1, uses one shared length for both arms); L_ARM here is
                # their average, 53.083 -> 53.08.
WG_WIDTH = 0.5   # um -- the same real strip width as the HCG's own access waveguide and the DC's arms
WG_RES = 100.0   # pts/um -- meow/MPB mode-solve resolution (no FDTD; matches Step 0/3's own checks)
MESH_WIDTH, MESH_HEIGHT = 4.0, 2.0

for path, name in [(HCG_SMATRIX, "HCG mirror"), (DC_SMATRIX, "directional coupler")]:
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{name}'s own production result is missing: {path}\n"
            f"This script assembles already-solved pieces, it does not solve them."
        )

# ---------------------------------------------------------------------------
# 1. The real HCG mirror block -- A3's direct 1-port construction.
#    from_smatrix() cannot represent this: it zero-fills a port's own
#    self-term by design, and this IS a port reflecting into itself.
# ---------------------------------------------------------------------------
smat_hcg = SMatrix.load(HCG_SMATRIX)
wl_hcg = smat_hcg.wavelengths
r = smat_hcg.coeff(port="in", in_mode=1, out_mode=1)
hcg_block = SMatrixBlock(wavelengths=wl_hcg, ports=("p",), modes_per_port={"p": 1}, data=r[:, None, None])
print(f"HCG mirror: {len(wl_hcg)} pts, {wl_hcg[0]:.3f}-{wl_hcg[-1]:.3f} um, "
      f"peak |r|^2={np.max(np.abs(r)**2):.4f} at {wl_hcg[np.argmax(np.abs(r)**2)]:.3f} um")

# ---------------------------------------------------------------------------
# 2. The real DC 4-port block -- built from the saved 2-port (in-excited)
#    result, using the VERIFIED symmetric-coupler equality (see module
#    docstring) rather than re-deriving the "isolated"-excited branch.
# ---------------------------------------------------------------------------
smat_dc = SMatrix.load(DC_SMATRIX)
wl_dc = smat_dc.wavelengths
t = smat_dc.coeff(port="through")
k = smat_dc.coeff(port="cross")

ports = ("through_left", "through_right", "input_reflection", "output")
L, R, I, O = 0, 1, 2, 3
M = len(wl_dc)
dc_data = np.zeros((M, 4, 4), dtype=complex)
dc_data[:, L, R] = t; dc_data[:, R, L] = t
dc_data[:, I, O] = t; dc_data[:, O, I] = t
dc_data[:, L, O] = k; dc_data[:, O, L] = k
dc_data[:, R, I] = k; dc_data[:, I, R] = k
dc_block = SMatrixBlock(wavelengths=wl_dc, ports=ports, modes_per_port={p: 1 for p in ports}, data=dc_data)
print(f"DC: {len(wl_dc)} pts, {wl_dc[0]:.3f}-{wl_dc[-1]:.3f} um, "
      f"|t|^2={np.abs(t[len(t)//2])**2:.4f} |k|^2={np.abs(k[len(k)//2])**2:.4f} at band centre")

# ---------------------------------------------------------------------------
# 3. Real connecting-waveguide mode solve -- fundamental mode only (Step 0's
#    own S-bend measurement: <0.01% of excited power reaches the 2nd-order
#    mode pair, so the connecting arm only ever needs to carry mode 1).
#    Computed at the DC's own (coarser) grid directly -- no need for the
#    HCG's finer 21-point grid here since this is the limiting grid anyway.
# ---------------------------------------------------------------------------
def _rect(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def _fundamental_neff(wl):
    length = 1.0
    gds = {(2, 0): [_rect(0.0, -WG_WIDTH / 2, length, WG_WIDTH / 2)],
           (1, 0): [_rect(0.0, -MESH_WIDTH / 2, length, MESH_WIDTH / 2)]}
    structs = extrude_polygons(gds, SOI_160.extrusion_rules)
    mesh = mw.Mesh2D(
        x=np.linspace(-MESH_WIDTH / 2, MESH_WIDTH / 2, int(WG_RES * MESH_WIDTH) + 1),
        y=np.linspace(-MESH_HEIGHT / 2, MESH_HEIGHT / 2, int(WG_RES * MESH_HEIGHT) + 1),
    )
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=float(wl), T=25.0)
    cs = mw.CrossSection.from_cell(cell=cells[0], env=env)
    from meow import compute_modes_meep
    modes = compute_modes_meep(cs, num_modes=1)
    return float(np.real(modes[0].neff))


t0 = time.time()
neff_wg = np.array([_fundamental_neff(wl) for wl in wl_dc])
print(f"Connecting waveguide n_eff: {len(wl_dc)} pts in {(time.time()-t0)/60:.1f} min "
      f"(pure MPB, no FDTD) -- {neff_wg[0]:.4f} to {neff_wg[-1]:.4f}")

# Genuinely bidirectional -- NOT from_transfer_matrix() (module docstring).
phase = np.exp(2j * np.pi * neff_wg / wl_dc * L_ARM)
wg_data = np.zeros((M, 2, 2), dtype=complex)
wg_data[:, 1, 0] = phase
wg_data[:, 0, 1] = phase
wg_block_left = SMatrixBlock(wavelengths=wl_dc, ports=("in", "out"), modes_per_port={"in": 1, "out": 1}, data=wg_data)
wg_block_right = SMatrixBlock(wavelengths=wl_dc, ports=("in", "out"), modes_per_port={"in": 1, "out": 1}, data=wg_data.copy())

# ---------------------------------------------------------------------------
# 4. Slice the HCG's finer grid down onto the DC's (the actual limiting
#    resolution) -- exact subset, verified.
# ---------------------------------------------------------------------------
idx = [int(np.argmin(np.abs(wl_hcg - w))) for w in wl_dc]
if not np.allclose(wl_hcg[idx], wl_dc, atol=1e-9):
    raise ValueError(
        "HCG's wavelength grid is not an exact superset of the DC's -- "
        "slicing would silently compare the wrong points."
    )
r_dc_grid = r[idx]
hcg_block_sliced = SMatrixBlock(wavelengths=wl_dc, ports=("p",), modes_per_port={"p": 1}, data=r_dc_grid[:, None, None])

# ---------------------------------------------------------------------------
# 5. Assemble -- the real Circuit wiring, step 4.
# ---------------------------------------------------------------------------
circuit = Circuit()
circuit.add_block("hcg_left", hcg_block_sliced)
circuit.add_block("wg_left", wg_block_left)
circuit.add_block("dc", dc_block)
circuit.add_block("wg_right", wg_block_right)
circuit.add_block("hcg_right", hcg_block_sliced)

circuit.connect("hcg_left", "p", "wg_left", "in")
circuit.connect("wg_left", "out", "dc", "through_left")
circuit.connect("dc", "through_right", "wg_right", "in")
circuit.connect("wg_right", "out", "hcg_right", "p")

cavity = circuit.reduce(external_ports={
    "dc.input_reflection": "input_reflection",
    "dc.output": "output",
})
R_circuit = np.abs(cavity.submatrix("input_reflection", "input_reflection")[:, 0, 0]) ** 2
T_circuit = np.abs(cavity.submatrix("output", "input_reflection")[:, 0, 0]) ** 2

# ---------------------------------------------------------------------------
# 6. Independent cross-check -- the hand-derived closed form (step 5),
#    verified separately against a brute-force 4-port solve in
#    tests/test_fabry_perot_network_crosscheck.py. Must agree with the
#    Circuit result to floating-point precision, or something in the
#    wiring/port-labeling above is wrong.
# ---------------------------------------------------------------------------
rho = r_dc_grid * phase ** 2
D = 1 - rho ** 2 * t ** 2
self_reflection = k ** 2 * rho / D
transmission = t * (1 + rho ** 2 * (k ** 2 - t ** 2)) / D
R_closed = np.abs(self_reflection) ** 2
T_closed = np.abs(transmission) ** 2

max_R_diff = float(np.max(np.abs(R_circuit - R_closed)))
max_T_diff = float(np.max(np.abs(T_circuit - T_closed)))
print(f"\nCircuit vs closed-form cross-check: max|dR|={max_R_diff:.2e}  max|dT|={max_T_diff:.2e}")
if max_R_diff > 1e-8 or max_T_diff > 1e-8:
    raise AssertionError(
        f"Circuit assembly and the closed form disagree by more than "
        f"floating-point precision (max|dR|={max_R_diff:.2e}, "
        f"max|dT|={max_T_diff:.2e}) -- a wiring or port-labeling bug, not a "
        f"physical result. See validation_markdown/"
        f"high_contrast_grating_fabry_perot_validation.md step 4."
    )
print("Circuit and closed-form agree to floating-point precision.")

# ---------------------------------------------------------------------------
# 7. Report, save, plot.
# ---------------------------------------------------------------------------
print(f"\n{'lambda':>8} {'R (input_reflection)':>21} {'T (output)':>11}")
for i, w in enumerate(wl_dc):
    print(f"{w:8.3f} {R_circuit[i]:21.6f} {T_circuit[i]:11.6f}")
i_peak = int(np.argmax(T_circuit))
print(f"\nPeak output transmission: {T_circuit[i_peak]:.4f} at {wl_dc[i_peak]:.3f} um")
print("Passivity check (R+T <= 1):", float(np.max(R_circuit + T_circuit)))
print("\nThis grid (10 nm spacing) is far coarser than this cavity's own "
      "estimated ~0.41 nm resonance FWHM -- see the dense reconstruction "
      "below (section 8) for the actual resolved lineshape.")

smat = SMatrix(wavelengths=wl_dc)
smat.set(self_reflection, port="input_reflection")
smat.set(transmission, port="output")
smat.save(os.path.join(output_dir, "smatrix_fabry_perot.pkl"))

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(wl_dc, R_circuit, "o-", label="Reflection (input_reflection)", color="darkorange")
ax.plot(wl_dc, T_circuit, "o-", label="Transmission (output)", color="steelblue")
ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (um)")
ax.set_ylabel("Power")
ax.set_title("HCG Fabry-Perot resonator -- assembled from real, independently\n"
              "solved pieces (10 nm grid -- known under-sampled near resonance)")
ax.legend(fontsize=8)
ax.grid(True)
ax.set_ylim(0, 1.05)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "spectrum.png"), dpi=300)
plt.close("all")
print(f"\nSMatrix saved -> {output_dir}/smatrix_fabry_perot.pkl")
print(f"Spectrum plot saved -> {output_dir}/spectrum.png")

# ---------------------------------------------------------------------------
# 8. Dense resonance characterization (step 5's follow-up) -- see module
#    docstring for the full derivation. Best-effort: the dense DC spectrum
#    is a separate script's output; if it hasn't been produced yet, skip
#    this section rather than failing the whole run.
# ---------------------------------------------------------------------------
from beef.spectral import interp_complex

DC_SMATRIX_DENSE = os.path.join(_ROOT, "example_output", "directional_coupler", "directional_coupler_soi160_hybrid", "smatrix_hybrid_dense.pkl")

if not os.path.exists(DC_SMATRIX_DENSE):
    print(f"\n(No dense DC spectrum at {DC_SMATRIX_DENSE} -- run "
          f"directional_coupler_soi160_hybrid_dense.py first to get the "
          f"resolved resonance lineshape. Skipping section 8.)")
else:
    print(f"\n--- Section 8: dense resonance characterization ---")
    smat_dc_dense = SMatrix.load(DC_SMATRIX_DENSE)
    wl_dense = smat_dc_dense.wavelengths
    t_dense = smat_dc_dense.coeff(port="through")
    k_dense = smat_dc_dense.coeff(port="cross")
    print(f"Dense DC grid: {len(wl_dense)} pts, {wl_dense[0]:.4f}-{wl_dense[-1]:.4f} um "
          f"({(wl_dense[-1]-wl_dense[0])/(len(wl_dense)-1)*1000:.4f} nm spacing)")

    # HCG r(lambda): interpolate directly (already confirmed safe -- max phase
    # step ~7.5 deg per 5 nm on the real 21-point data, nowhere near
    # unwrap_is_safe's danger zone -- check_unwrap asserts this rather than
    # assuming it).
    if wl_dense[0] < wl_hcg[0] or wl_dense[-1] > wl_hcg[-1]:
        raise ValueError(
            f"Dense grid [{wl_dense[0]}, {wl_dense[-1]}] extends beyond the "
            f"HCG's own measured span [{wl_hcg[0]}, {wl_hcg[-1]}] -- "
            f"interp_complex would clamp rather than extrapolate."
        )
    r_dense = interp_complex(wl_hcg, r, wl_dense, check_unwrap=0.5)

    # Connecting waveguide neff: reuse the ALREADY-COMPUTED 11-point neff_wg
    # above (no new MPB solves) -- fit a degree-2 polynomial (same convention
    # as ring_resonator_hybrid.py's own neff_fit) and evaluate analytically
    # at the dense grid.
    _neff_fit = np.polyfit(wl_dc, neff_wg, 2)
    _neff_resid = neff_wg - np.polyval(_neff_fit, wl_dc)
    _neff_rel = float(np.max(np.abs(_neff_resid))) / float(np.ptp(neff_wg))
    print(f"Connecting-waveguide neff degree-2 fit: max residual "
          f"{np.max(np.abs(_neff_resid)):.2e} ({_neff_rel:.2%} of neff's own range)")
    if _neff_rel > 0.01:
        raise ValueError(
            f"neff_wg's degree-2 fit residual is {_neff_rel:.2%} of its own "
            f"range -- too large to trust the analytic phase evaluation "
            f"below. The 11-point neff_wg computed above may need more "
            f"points or a higher fit degree before this section can trust it."
        )
    neff_dense = np.polyval(_neff_fit, wl_dense)
    phase_dense = np.exp(2j * np.pi * neff_dense / wl_dense * L_ARM)

    rho_dense = r_dense * phase_dense ** 2
    D_dense = 1 - rho_dense ** 2 * t_dense ** 2
    self_reflection_dense = k_dense ** 2 * rho_dense / D_dense
    transmission_dense = t_dense * (1 + rho_dense ** 2 * (k_dense ** 2 - t_dense ** 2)) / D_dense
    R_dense = np.abs(self_reflection_dense) ** 2
    T_dense = np.abs(transmission_dense) ** 2
    print(f"Passivity check, dense grid (max R+T, should be <=1): {np.max(R_dense + T_dense):.4f}")

    # This dense evaluation is the closed form ONLY -- NOT a re-run of
    # Circuit at thousands of points. Circuit and the closed form already
    # agree to floating-point precision on the sparse grid above (section 6);
    # re-verifying that same equivalence at every dense point would be
    # redundant, not new information.

    def _find_dips(wl, T):
        """Numpy-only local-minima + FWHM/Q extraction -- same convention as
        ring_resonator_hybrid.py's own find_dips()."""
        dips = []
        baseline = np.median(T)
        for i in range(1, len(T) - 1):
            if T[i] < T[i - 1] and T[i] < T[i + 1] and T[i] < 0.7 * baseline:
                half = (T[i] + baseline) / 2.0
                lo = i
                while lo > 0 and T[lo] < half:
                    lo -= 1
                hi = i
                while hi < len(T) - 1 and T[hi] < half:
                    hi += 1
                if hi > lo:
                    fwhm = wl[hi] - wl[lo]
                    q = wl[i] / fwhm if fwhm > 0 else np.nan
                    dips.append((wl[i], T[i], fwhm, q))
        return dips

    dips = _find_dips(wl_dense, T_dense)
    print(f"\nDense-grid assembly ({len(wl_dense)} points): {len(dips)} resonance dip(s) found.")
    if dips:
        print(f"  {'lambda (um)':>12}  {'T_min':>8}  {'FWHM (nm)':>10}  {'Q':>10}")
        for wl0, Tmin, fwhm, q in dips:
            print(f"  {wl0:12.4f}  {Tmin:8.4f}  {fwhm*1000:10.4f}  {q:10.1f}")
    if len(dips) >= 2:
        fsr_meas = float(np.mean(np.diff([d[0] for d in dips])))
        print(f"  Measured mean FSR: {fsr_meas*1000:.3f} nm (Step 1's own estimate: 9.5 nm)")

    smat_dense = SMatrix(wavelengths=wl_dense)
    smat_dense.set(self_reflection_dense, port="input_reflection")
    smat_dense.set(transmission_dense, port="output")
    smat_dense.save(os.path.join(output_dir, "smatrix_fabry_perot_dense.pkl"))

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(wl_dense, R_dense, label="Reflection (input_reflection)", color="darkorange", lw=0.8)
    ax.plot(wl_dense, T_dense, label="Transmission (output)", color="steelblue", lw=0.8)
    ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
    for wl0, Tmin, fwhm, q in dips:
        ax.axvline(wl0, color="k", linestyle="--", lw=0.5, alpha=0.5)
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Power")
    ax.set_title(f"HCG Fabry-Perot resonator -- dense reconstruction "
                 f"({len(wl_dense)} pts, interpolated) -- {len(dips)} resonance(s) resolved")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "spectrum_dense.png"), dpi=200)
    plt.close("all")
    print(f"\nSMatrix saved -> {output_dir}/smatrix_fabry_perot_dense.pkl")
    print(f"Spectrum plot saved -> {output_dir}/spectrum_dense.png")

# ---------------------------------------------------------------------------
# 9. Pull in each constituent device's own geometry/performance plots, plus a
#    new overall schematic of the assembled cavity -- so this folder is a
#    self-contained picture of the whole device, not just the combined
#    spectrum. Copies (not moves) -- the constituent devices' own output
#    folders remain the authoritative source for their own plots.
# ---------------------------------------------------------------------------
import shutil

hcg_plots_dir = os.path.join(output_dir, "hcg_mirror")
os.makedirs(hcg_plots_dir, exist_ok=True)
_hcg_src_dir = os.path.dirname(HCG_SMATRIX)
for _fname in ("geometry.png", "reflection.png"):
    _src = os.path.join(_hcg_src_dir, _fname)
    if os.path.exists(_src):
        shutil.copy2(_src, os.path.join(hcg_plots_dir, _fname))

dc_plots_dir = os.path.join(output_dir, "directional_coupler")
os.makedirs(dc_plots_dir, exist_ok=True)
_dc_src_dir = os.path.dirname(DC_SMATRIX)
for _fname in ("geometry_hybrid.png", "hybrid_results.png", "hybrid_results_dense.png"):
    _src = os.path.join(_dc_src_dir, _fname)
    if os.path.exists(_src):
        shutil.copy2(_src, os.path.join(dc_plots_dir, _fname))
_dc_sections_src = os.path.join(_dc_src_dir, "sections")
if os.path.isdir(_dc_sections_src):
    _dc_sections_dst = os.path.join(dc_plots_dir, "sections")
    os.makedirs(_dc_sections_dst, exist_ok=True)
    for _fname in os.listdir(_dc_sections_src):
        if _fname.endswith(".png"):
            shutil.copy2(os.path.join(_dc_sections_src, _fname), os.path.join(_dc_sections_dst, _fname))

print(f"\nHCG mirror's own geometry/performance plots copied -> {hcg_plots_dir}/")
print(f"Directional coupler's own geometry/performance/section plots copied -> {dc_plots_dir}/")

# Overall schematic: the assembled cavity, laid out to real x-scale. This is
# NOT a real shared FDTD domain -- each piece above is independently solved
# in its own coordinate frame (A2) -- it is a labeled block diagram matching
# this doc's own header ASCII sketch, built from each piece's own already-
# measured real length (HCG_LENGTH from high_contrast_grating_fdtd.py's own
# printed "HCG total length" -- 5.389 um; DC_LENGTH from
# directional_coupler_soi160_hybrid.py's own printed "total length" -- 23.00
# um; L_ARM from Step 1's real GDS-derived pin above).
HCG_LENGTH = 5.389   # um -- high_contrast_grating_fdtd.py's own measured GDS x-span
DC_LENGTH = 23.00    # um -- directional_coupler_soi160_hybrid.py's own dc.total_length


def _draw_block(ax, x0, width, label, color, height=1.0, y0=-0.5):
    ax.add_patch(plt.Rectangle((x0, y0), width, height, facecolor=color, edgecolor="k", lw=1.0))
    ax.annotate(label, (x0 + width / 2, y0 + height + 0.15), ha="center", fontsize=9)
    ax.annotate(f"{width:.2f} um", (x0 + width / 2, y0 - 0.3), ha="center", fontsize=8, color="gray")
    return x0 + width


fig, ax = plt.subplots(figsize=(14, 3.5))
x1 = _draw_block(ax, 0.0, HCG_LENGTH, "HCG mirror\n(FDTD)", "#4C72B0")
x2 = _draw_block(ax, x1, L_ARM, "waveguide\n(EME/analytic)", "#DD8452", height=0.4, y0=-0.2)
x3 = _draw_block(ax, x2, DC_LENGTH, "directional coupler\n(FDTD + EME)", "#55A868")
x4 = _draw_block(ax, x3, L_ARM, "waveguide\n(EME/analytic)", "#DD8452", height=0.4, y0=-0.2)
x5 = _draw_block(ax, x4, HCG_LENGTH, "HCG mirror\n(FDTD)", "#4C72B0")

ax.annotate("", xy=(x2, -1.4), xytext=(x2, -0.5), arrowprops=dict(arrowstyle="-|>", color="k"))
ax.annotate("input_reflection", (x2, -1.6), ha="center", fontsize=8)
ax.annotate("", xy=(x3, -1.4), xytext=(x3, -0.5), arrowprops=dict(arrowstyle="-|>", color="k"))
ax.annotate("output", (x3, -1.6), ha="center", fontsize=8)

ax.set_xlim(-2, x5 + 2)
ax.set_ylim(-2.0, 1.8)
ax.set_xlabel("Position along cavity (um)")
ax.set_yticks([])
ax.set_title(f"HCG Fabry-Perot resonator -- overall assembled geometry "
             f"(schematic, to x-scale, NOT one shared FDTD domain)\n"
             f"total cavity length {x5:.2f} um")
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "overall_geometry.png"), dpi=200)
plt.close("all")
print(f"Overall assembled-cavity schematic saved -> {output_dir}/overall_geometry.png")
