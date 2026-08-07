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

Connecting waveguide (20 um each arm, Step 1's provisional pin): a real
fundamental-mode-only meow/MPB sweep (no FDTD -- straight section), run
fresh here rather than depending on a scratch file. Built as a genuinely
BIDIRECTIONAL SMatrixBlock -- NOT via beef.network.convert.
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
  - The 11-point DC grid (10 nm spacing) is far coarser than the ~0.41 nm
    FWHM estimated from this cavity's own finesse -- the spectrum below is
    known to be under-sampled near the resonance it itself reveals near
    1.57 um. A dense zoom sweep is a separate, already-running effort
    (validation doc step 5's follow-up).
  - FDTD_RESOLUTION=40/EME_RESOLUTION=80 (inherited from
    directional_coupler_soi160_hybrid.py) are ASSUMED, not independently
    converged for this device.
  - The 20 um arm length is a provisional pin (Step 1), not a target
    derived from a desired FSR/finesse.

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

L_ARM = 20.0   # um, Step 1's provisional pin -- both connecting-waveguide arms, symmetric cavity
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
print("\nKNOWN GAP: this grid (10 nm spacing) is far coarser than this "
      "cavity's own estimated ~0.41 nm resonance FWHM -- see the "
      "validation doc's step 5 follow-up for the dedicated dense sweep.")

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
