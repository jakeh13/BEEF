"""HCG Fabry-Perot — hybrid body (coupler region) field reconstruction.

Pinned to high_contrast_grating_fabry_perot/python/hcg_fabry_perot_hybrid.py's
own `dc_eme` body section. See beef/plot/hybrid_body.py for why
plot_eme_fields() is not sufficient here, and for the standing caveat that the
excitation is an APPROXIMATION rather than the hybrid's own interface field.

WHAT THIS DEVICE'S "BODY" IS. The HCG Fabry-Perot's cascade puts an EME
section on the straight coupled pair of a DirectionalCoupler built to
duty/dc.gds dimensions -- the same physical coupler as
directional_coupler_soi160_hybrid.py, at the same 2 modes and the same
COUPLING_LENGTH=3.0 um. So this script is deliberately near-identical to
plot_dc_hybrid_body_fields.py's soi160 branch, and the two should agree on the
supermode pair; they are kept separate because the HCG script owns its own
band (1.530-1.565 um, N_WL=36) and its own SIM_HEIGHT/oxide width, and
pinning to the wrong one of those is exactly the kind of drift these scripts
exist to make visible.

The grating mirrors themselves are NOT covered here and cannot be: they are
FDTD sections, they are not a uniform cross-section, and the reconstruction
below requires one cross-section to serve both ports. A grating's field map
needs step 2's FDTD capture, not this.

Launch is OFF-AXIS on one arm (excitation_center = +coupling_pitch/2), so both
supermodes are driven and the picture shows real inter-arm coupling; a centred
launch would excite only the even supermode.

Cost: one cross-section, 2 meep/MPB modes, one wavelength (~35 s measured).

Run:
    python .../plot_hcgfp_hybrid_body_fields.py
"""
import os

import numpy as np

from beef import DirectionalCoupler, EMESolver, SOI_160, plot_hybrid_body_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output",
    "high_contrast_grating_fabry_perot", "hcg_fabry_perot_hybrid_body"))
os.makedirs(output_dir, exist_ok=True)

# --- pinned to hcg_fabry_perot_hybrid.py ------------------------------------
WG_WIDTH = 0.5           # same real strip width as the HCG's access waveguide
GAP = 0.25               # duty/dc.gds coupling-region edge-to-edge gap
COUPLING_LENGTH = 3.0    # duty/dc.gds real straight coupled-pair length, FIXED
BODY_STUB = 0.25
SIM_HEIGHT = 1.5
N_MODES = 2              # coupler body: fundamental even/odd pair ONLY, MEASURED
EME_RESOLUTION = float(os.environ.get("HCGFP_EME_RES", 80.0))
WL_MIN, WL_MAX, N_WL = 1.530, 1.565, 36
WL_PLOT = float(os.environ.get("HCGFP_WL_PLOT", 1.5475))   # band centre

dc = DirectionalCoupler(
    waveguide_width=WG_WIDTH, gap=GAP, coupling_length=COUPLING_LENGTH,
    bend_length=10.0, waveguide_offset=4.0,
    bend_type="euler", waveguide_height=0.16,
)

half_body = dc.coupling_length / 2
x0, x1 = -half_body + BODY_STUB, +half_body - BODY_STUB

solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=7.0, oxide_height=SIM_HEIGHT,
    res=EME_RESOLUTION, material_stack=SOI_160, mode_solver="meep",
)

arm_y = dc.coupling_pitch / 2.0
print(f"HCG-FP coupler body: x in [{x0:+.2f}, {x1:+.2f}] um (L={x1 - x0:.2f})")
print(f"arms at y = +/-{arm_y:.3f} um, launching into the +{arm_y:.3f} arm")

res = plot_hybrid_body_fields(
    dc, solver, x0, x1, WL_PLOT, N_MODES, output_dir,
    tag=f"hcgfp_coupler_{WL_PLOT:.4f}um",
    excitation="gaussian", excitation_center=+arm_y, excitation_fwhm=WG_WIDTH,
    wavelengths=np.linspace(WL_MIN, WL_MAX, N_WL),
    port_ys=[-arm_y, +arm_y],
    body_half_width=arm_y + WG_WIDTH / 2,
    device_label="HCG Fabry-Perot (coupler body)",
    total_length=dc.total_length,
)

w, neff = res["weights"], res["neff"]
assert res["captured_power"] <= 1.0 + 1e-6, \
    f"captured power {res['captured_power']:.4f} > 1 -- wrong inner product"

ratio = w[0] / max(w[1], 1e-30)
dneff = abs(neff[0] - neff[1])
print(f"\neven/odd supermode weights: {w[0]:.4f} / {w[1]:.4f} (ratio {ratio:.3f})")
print(f"supermode Delta-n_eff = {dneff:.6f}  ->  "
      f"beat length L_pi = {WL_PLOT / (2 * dneff):.3f} um "
      f"(body is {x1 - x0:.2f} um)")
assert 0.4 < ratio < 2.5, (
    f"single-arm launch gave even/odd ratio {ratio:.3f}, expected ~1 -- "
    "check excitation_center against coupling_pitch/2")
print("\nchecks passed")
