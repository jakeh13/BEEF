"""Directional coupler — hybrid body (coupling region) field reconstruction.

Companion to plot_mmi1x2_hybrid_body_selfimaging.py, pinned to the coupler
scripts' own body sections. Covers BOTH couplers via DC_VARIANT:

    DC_VARIANT=soi220   -> directional_coupler_hybrid.py        (N_MODES=5)
    DC_VARIANT=soi160   -> directional_coupler_soi160_hybrid.py (N_MODES=2)

See beef/plot/hybrid_body.py for why plot_eme_fields() is not sufficient
(meow-frame labels; single-mode-only excitation) and for the standing caveat
that the excitation is an APPROXIMATION, not the hybrid's own interface field.

WHY A COUPLER IS THE MOST INFORMATIVE CASE HERE. The body is the straight
coupled pair, and its modes are the even/odd supermodes. A launch on ONE arm
is an equal mix of even and odd, so the reconstruction shows power sloshing
between the two arms at the beat length set by Delta-beta -- which is the
quantity `resolved/directional_coupler_validation.md` items 6/14 measure and
the quantity calibrate_body_neff exists to correct. Unlike an MMI's
self-imaging, this beating is the coupler's entire function, so the picture is
a direct read on the physics rather than an illustration of it.

The launch is therefore OFF-AXIS: excitation_center = +coupling_pitch/2, the
driven arm's own centre (coupling_pitch = waveguide_width + gap). A centred
launch would be wrong here -- it would excite only the even supermode and show
no coupling at all.

Cost: one cross-section, 2 or 5 meep/MPB modes, one wavelength. The cheapest
field solve of any device in the repo (soi160/HCG carry only 2 modes).

Run:
    python .../plot_dc_hybrid_body_fields.py
    DC_VARIANT=soi160 python .../plot_dc_hybrid_body_fields.py
    DC_BODY_LENGTH=10.0 python .../plot_dc_hybrid_body_fields.py   # override
"""
import os

import numpy as np

from beef import DirectionalCoupler, EMESolver, SOI_220, SOI_160, plot_hybrid_body_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
VARIANT = os.environ.get("DC_VARIANT", "soi220")

# ---------------------------------------------------------------------------
# Per-variant constants, each pinned to that variant's own hybrid script.
# ---------------------------------------------------------------------------
if VARIANT == "soi220":
    # directional_coupler_hybrid.py
    WG_WIDTH, GAP, WG_HEIGHT = 0.45, 0.25, 0.22
    BEND_LENGTH, WG_OFFSET = 9.0, 3.0
    STACK, OXIDE_W = SOI_220, 6.0
    N_MODES = 5           # even + odd supermodes at the two body interfaces
    EME_RES = 80.0        # res 40 is NOT converged on this device -- see script
    WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
    # The script sweeps coupling_length; its build reference is the sweep
    # midpoint, coupling_lengths = linspace(2*BODY_STUB, 60.0, 200)[100].
    COUPLING_LENGTH = float(np.linspace(2 * 0.25, 60.0, 200)[100])
    OUT_SUB = "directional_coupler_hybrid"
elif VARIANT == "soi160":
    # directional_coupler_soi160_hybrid.py -- duty/dc.gds dimensions
    WG_WIDTH, GAP, WG_HEIGHT = 0.5, 0.25, 0.16
    BEND_LENGTH, WG_OFFSET = 10.0, 4.0
    STACK, OXIDE_W = SOI_160, 7.0
    N_MODES = 2           # fundamental even/odd supermode pair ONLY -- the
                          # real 0.5 um / 160 nm strip carries 2 guided modes
    EME_RES = 80.0
    WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
    COUPLING_LENGTH = 3.0  # FIXED, not swept
    OUT_SUB = "directional_coupler_soi160_hybrid"
else:
    raise SystemExit(f"unknown DC_VARIANT={VARIANT!r}, expected soi220|soi160")

BODY_STUB = 0.25
WL_PLOT = float(os.environ.get("DC_WL_PLOT", 1.55))   # C-band centre
COUPLING_LENGTH = float(os.environ.get("DC_BODY_LENGTH", COUPLING_LENGTH))

output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "directional_coupler",
    OUT_SUB + "_body"))
os.makedirs(output_dir, exist_ok=True)

dc = DirectionalCoupler(
    waveguide_width=WG_WIDTH, gap=GAP, coupling_length=COUPLING_LENGTH,
    bend_length=BEND_LENGTH, waveguide_offset=WG_OFFSET,
    bend_type="euler", waveguide_height=WG_HEIGHT,
)

# The hybrid's body section: boundaries[1]..boundaries[2]. Body is section
# k=1 of 3, so HybridSolver._stubs(k) gives 0.0/0.0.
half_body = dc.coupling_length / 2
x0, x1 = -half_body + BODY_STUB, +half_body - BODY_STUB

solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=OXIDE_W, oxide_height=1.5,
    res=EME_RES, material_stack=STACK, mode_solver="meep",
)

arm_y = dc.coupling_pitch / 2.0
print(f"DC {VARIANT}: coupling body x in [{x0:+.2f}, {x1:+.2f}] um "
      f"(L={x1 - x0:.2f} of {dc.total_length:.2f} um device)")
print(f"arms at y = +/-{arm_y:.3f} um (coupling_pitch {dc.coupling_pitch:.3f}), "
      f"launching into the +{arm_y:.3f} arm")

res = plot_hybrid_body_fields(
    dc, solver, x0, x1, WL_PLOT, N_MODES, output_dir,
    tag=f"dc_{VARIANT}_{WL_PLOT:.3f}um",
    # OFF-AXIS launch: one arm only. See module docstring.
    excitation="gaussian", excitation_center=+arm_y, excitation_fwhm=WG_WIDTH,
    wavelengths=np.linspace(WL_MIN, WL_MAX, N_WL),
    port_ys=[-arm_y, +arm_y],
    body_half_width=arm_y + WG_WIDTH / 2,
    device_label=f"Directional coupler ({VARIANT})",
    total_length=dc.total_length,
)

# --- physics checks ---------------------------------------------------------
w, neff = res["weights"], res["neff"]
assert res["captured_power"] <= 1.0 + 1e-6, \
    (f"captured power {res['captured_power']:.4f} > 1 -- the projection inner "
     "product is wrong; see body_mode_weights()'s docstring")

# A single-arm launch must split roughly evenly between the even and odd
# supermodes -- that near-equality IS what makes a coupler couple. Checked
# rather than assumed, since an off-axis launch is easy to mis-centre.
if N_MODES >= 2:
    ratio = w[0] / max(w[1], 1e-30)
    print(f"\neven/odd supermode weight ratio: {ratio:.3f} "
          f"(w0={w[0]:.4f}, w1={w[1]:.4f})")
    assert 0.4 < ratio < 2.5, (
        f"single-arm launch gave even/odd ratio {ratio:.3f}, expected ~1 -- "
        "a launch centred on one arm should excite both supermodes about "
        "equally. Check excitation_center against coupling_pitch/2.")
    # Beat length from the supermode split, the coupler's defining number.
    dneff = abs(neff[0] - neff[1])
    if dneff > 0:
        dbeta = 2 * np.pi * dneff / WL_PLOT
        print(f"supermode Delta-n_eff = {dneff:.6f}  ->  "
              f"Delta-beta = {dbeta:.6f} rad/um, "
              f"beat length L_pi = {WL_PLOT / (2 * dneff):.3f} um "
              f"(body is {x1 - x0:.2f} um)")

        # INDEPENDENT CROSS-CHECK, soi220 only. resolved/directional_coupler_
        # validation.md item 11 measured this device's converged MPB Delta-beta
        # directly: "MPB plateaus at Delta-beta = 0.0894 by res 80-100". This
        # script runs mode_solver="meep" (meow's MEEP/MPB backend) at
        # EME_RES=80, i.e. exactly that converged-MPB case, so the two numbers
        # are comparable and agreeing is a real check on the reconstruction --
        # the same kind of check the 1x2's n_eff match against gap 3 provides.
        # Measured 2026-08-31: 0.089541 vs 0.0894, a 0.16% agreement.
        if VARIANT == "soi220" and abs(EME_RES - 80.0) < 1e-9 \
                and abs(WL_PLOT - 1.55) < 1e-9:
            DOC_DBETA = 0.0894
            rel = abs(dbeta - DOC_DBETA) / DOC_DBETA
            print(f"  vs directional_coupler_validation.md item 11's converged "
                  f"MPB Delta-beta {DOC_DBETA}: {100 * rel:.2f}%")
            assert rel < 0.01, (
                f"Delta-beta {dbeta:.6f} disagrees with the doc's independently "
                f"measured converged-MPB {DOC_DBETA} by {100 * rel:.1f}% "
                "(>1%) -- the body mode basis or the geometry has moved")
print("\nchecks passed")
