"""1x2 MMI — hybrid body section field reconstruction (device frame).

Step 1 of validation_markdown/hybrid_field_extraction_validation.md, and the
reference case the shared helper was built and validated against. Standalone:
needs NO FDTD field capture (step 2), no library change beyond
beef.plot_hybrid_body_fields itself, and nothing from the cascade cache.

See beef/plot/hybrid_body.py's module docstring for why plot_eme_fields() is
not sufficient here (it labels in meow's frame, and can excite only one mode
at a time) and for the standing caveat that the excitation is an
APPROXIMATION rather than the hybrid's own interface field.

Every constant below is mmi_1x2_hybrid.py's own value -- see that script for
the full provenance (N_MODES=8 comes from resolved/mmi_1x2_validation.md
gap 3; RESOLUTION=40 from gap 2; mode_solver="meep" matches FDTDSolver's own
eigenmode basis in this hybrid chain).

REGRESSION VALUES (measured 2026-08-31, res 40, 8 modes, 1.31 um). This script
asserts on them, so a change in the helper or the mode solver fails loudly
here rather than silently redrawing a different picture:

    n_eff       2.696915 2.637710 2.536282 2.387965
                2.185091 1.916790 1.775436 1.706013
    weights     0.9022 (m0), 0.0935 (m2), 0.0013 (m4), odd modes < 1e-10
    captured    0.9970
    gram offdiag 0.069

The n_eff row independently reproduces gap 3's meow column to all six printed
decimals (job 11608940), which is what makes this a real check on the mode
basis and not just a self-consistency test.

Cost: one cross-section, 8 meep/MPB modes, one wavelength (~55 s measured).

Run:
    python .../plot_mmi1x2_hybrid_body_selfimaging.py
"""
import os

import numpy as np

from beef import MMI, EMESolver, SOI_160, plot_hybrid_body_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_1x2", "mmi_1x2_hybrid_body"))
os.makedirs(output_dir, exist_ok=True)

# --- pinned to mmi_1x2/python/mmi_1x2_hybrid.py -----------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.26, 1.36, 11
N_MODES = 8
BODY_STUB = 0.25
WL_PLOT = 1.31
# taper_width is consumed into MMIPort objects rather than kept as an MMI
# attribute, so it is named here for the launch width.
TAPER_WIDTH = 0.8

mmi = MMI(
    n_inputs=1, n_outputs=2,
    mm_length=4.5, mm_width=2.0, waveguide_width=0.35,
    taper_length=1.0, taper_width=0.8, waveguide_offset=1.0,
    waveguide_height=0.16,
)

# The hybrid's body section: boundaries[1] .. boundaries[2]. The body is
# section k=1 of 3, and HybridSolver._stubs(k) gives a nonzero stub only to
# the FIRST and LAST sections, so both stubs are zero here.
half_body = mmi.mm_length / 2
x0, x1 = -half_body + BODY_STUB, +half_body - BODY_STUB

# mmi_1x2_hybrid.py's body_solver, verbatim. num_modes/num_cells/lengths are
# placeholders on the solve_section() path.
solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=5.0, oxide_height=1.5,
    res=RESOLUTION, material_stack=SOI_160,
    mode_solver="meep",
)

print(f"1x2 hybrid BODY: x in [{x0:+.2f}, {x1:+.2f}] um  "
      f"(L={x1 - x0:.2f} of {mmi.total_length:.2f} um device)")

res = plot_hybrid_body_fields(
    mmi, solver, x0, x1, WL_PLOT, N_MODES, output_dir,
    tag=f"mmi1x2_mode1_{WL_PLOT:.3f}um",
    # A 1x2 has ONE centred input, so the launch is centred: by symmetry it
    # drives the even body modes only, and their beating is the self-imaging
    # that splits power into the two outputs.
    excitation="gaussian", excitation_center=0.0, excitation_fwhm=TAPER_WIDTH,
    # The FULL production sweep, not just WL_PLOT -- _extrusion_rules()
    # consumes it, and matching the cascade's own array keeps this faithful.
    wavelengths=np.linspace(WL_MIN, WL_MAX, N_WL),
    port_ys=mmi.output_ys,          # read off the geometry: [-0.5, +0.5]
    body_half_width=mmi.mm_width / 2,
    device_label="1×2 MMI", total_length=mmi.total_length,
)

# --- regression assertions (see module docstring) ---------------------------
w, neff = res["weights"], res["neff"]
expected_neff = np.array([2.696915, 2.637710, 2.536282, 2.387965,
                          2.185091, 1.916790, 1.775436, 1.706013])
assert np.allclose(neff, expected_neff, atol=2e-6), \
    f"body n_eff moved:\n  got      {neff}\n  expected {expected_neff}"
assert abs(res["captured_power"] - 0.9970) < 5e-4, \
    f"captured power {res['captured_power']:.4f} != 0.9970"
assert w[1::2].max() < 1e-10, \
    (f"a centred launch excited an ODD body mode ({w[1::2].max():.2e}) -- "
     "the projection inner product is wrong; see body_mode_weights()'s "
     "docstring on the 1-D-slice failure")
assert abs(w[0] - 0.9022) < 5e-4 and abs(w[2] - 0.0935) < 5e-4, \
    f"even-mode weights moved: m0={w[0]:.4f} m2={w[2]:.4f}"

# The self-image must land on the output ports: this is a check of body length
# against port placement, read off the field rather than the design intent.
I_exit = res["intensity"][-1]
y = res["y"]
peak_y = abs(y[np.argmax(np.where(y > 0, I_exit, 0))])
assert abs(peak_y - max(mmi.output_ys)) < 0.08, \
    (f"exit lobe peaks at y={peak_y:.3f} but the output port is at "
     f"{max(mmi.output_ys):.3f} -- self-image does not land on the port")

print(f"\nall regression checks passed "
      f"(exit lobe y={peak_y:.3f} vs port {max(mmi.output_ys):.3f})")
