"""1x4 MMI (SiN 300) — hybrid body field reconstruction, device frame.

Pinned to mmi_1x4/python/mmi_1x4_hybrid.py's own body section. See
beef/plot/hybrid_body.py for why plot_eme_fields() is not sufficient here, and
for the standing caveat that the excitation is an APPROXIMATION rather than the
hybrid's own interface field.

⚠ mode_solver="tidy3d", NOT "meep" -- AND THAT IS NOT A COPY ERROR.
Every other device's body-field script here uses "meep" to match FDTDSolver's
own eigenmode basis. mmi_1x4_hybrid.py deliberately does not, and this script
follows it rather than "fixing" it. Its own module docstring gives the reason:

    "Step 4 found compute_modes_meep silently wrong for this platform's
    weakly-confined modes, and the body's near-cutoff bands (job 11660226:
    margins down to 0.004) are exactly the regime that broke it -- using it
    here would risk a wrong body mode basis, not just the ~2.3%/mode phase
    mismatch calibrate_body_neff exists to fix."

So on this device the body basis is tidy3d's, and the tidy3d-vs-meep basis
mismatch is a known, accepted, calibrate_body_neff-corrected property of the
cascade -- not something to reconcile in a plotting script. Switching this to
"meep" would produce a prettier-looking agreement with the FDTD sections and a
WRONG body basis.

GEOMETRY IS THE POST-2026-08-14 ONE (MM_LENGTH=39.90, WAVEGUIDE_OFFSET=3.1),
not the earlier 2.4 pair. Four 1x4 scripts were found stranded on the old
geometry as recently as commit ca42ea9 (2026-08-31), so this is worth stating:
if mmi_1x4_hybrid.py's constants move, this script must move with them.

Cost: a 12 um-wide cross-section at res 30 with 12 modes -- the most expensive
body solve of the set, but still one cross-section and one wavelength.

Run:
    python .../plot_mmi1x4_hybrid_body_fields.py
    MMI1X4_N_MODES=8 python .../plot_mmi1x4_hybrid_body_fields.py
"""
import os

import numpy as np

from beef import MMI, EMESolver, SIN_300, plot_hybrid_body_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_1x4",
    "mmi_1x4_hybrid_body"))
os.makedirs(output_dir, exist_ok=True)

# --- pinned to mmi_1x4/python/mmi_1x4_hybrid.py -----------------------------
RESOLUTION = 30.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
WAVEGUIDE_OFFSET = 3.1        # CHANGED 2026-08-14 from 2.4
MM_WIDTH = 12.0
MM_LENGTH = 39.90             # CHANGED 2026-08-14
WAVEGUIDE_HEIGHT = 0.3        # SiN 300 nm
FDTD_SIM_WIDTH = MM_WIDTH + 4.0     # 16.0
FDTD_SIM_HEIGHT = 4.0
BODY_STUB = 0.25
N_MODES = int(os.environ.get("MMI1X4_N_MODES", 12))
WL_PLOT = float(os.environ.get("MMI1X4_WL_PLOT", 1.55))

mmi = MMI(
    n_inputs=1, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH,
    taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=WAVEGUIDE_HEIGHT,
)

half_body = mmi.mm_length / 2
x0, x1 = -half_body + BODY_STUB, +half_body - BODY_STUB

# mmi_1x4_hybrid.py's body sub-solver, verbatim -- including the oxide
# dimensions, which it deliberately sets to the FDTD sim box (2026-08-19
# finding, see that script's _fdtd_common comment), not a separate value.
solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=FDTD_SIM_WIDTH, oxide_height=FDTD_SIM_HEIGHT,
    res=RESOLUTION, material_stack=SIN_300,
    mode_solver="tidy3d",     # step 4's fix -- see module docstring, do NOT change
)

print(f"1x4 MMI body: x in [{x0:+.2f}, {x1:+.2f}] um (L={x1 - x0:.2f} of "
      f"{mmi.total_length:.2f} um device), {N_MODES} modes, tidy3d")
print(f"output ports at y = {mmi.output_ys} um")

res = plot_hybrid_body_fields(
    mmi, solver, x0, x1, WL_PLOT, N_MODES, output_dir,
    tag=f"mmi1x4_{WL_PLOT:.3f}um_n{N_MODES}",
    # ONE centred input (n_inputs=1), so a centred launch: by symmetry it
    # drives the even body modes, whose beating forms the 4-way self-image.
    excitation="gaussian", excitation_center=0.0, excitation_fwhm=TAPER_WIDTH,
    wavelengths=np.linspace(WL_MIN, WL_MAX, N_WL),
    port_ys=mmi.output_ys,
    body_half_width=mmi.mm_width / 2,
    device_label="1×4 MMI (SiN 300)", total_length=mmi.total_length,
)

w = res["weights"]
assert res["captured_power"] <= 1.0 + 1e-6, \
    f"captured power {res['captured_power']:.4f} > 1 -- wrong inner product"
assert w[1::2].max() < 1e-6, (
    f"a centred launch excited an ODD body mode ({w[1::2].max():.2e}) -- "
    "the projection inner product is wrong; see body_mode_weights()")

# Where the self-image actually lands, versus where the four ports are. This
# device's own validation carries an unresolved inner-vs-outer port imbalance
# (validation doc step 6d: EME favours ports 1,2 and hybrid favours 0,3), so
# report the lobe positions rather than asserting a match -- this is a
# diagnostic for that open question, not a passing gate.
y, I_exit = res["y"], res["intensity"][-1]
print(f"\noutput ports:      {[round(p, 3) for p in mmi.output_ys]}")
lobes = []
for lo, hi in [(0.0, mmi.mm_width / 2)]:
    m = (y >= lo) & (y <= hi)
    # crude peak picking on the positive half; the device is symmetric
    seg, ys = I_exit[m], y[m]
    for i in range(1, len(seg) - 1):
        if seg[i] > seg[i - 1] and seg[i] >= seg[i + 1] and seg[i] > 0.15 * I_exit.max():
            lobes.append(round(float(ys[i]), 3))
print(f"exit lobes (y>=0):  {lobes}")
print(f"basis off-diagonal: {res['gram_max_offdiag']:.4f}")
print("\nchecks passed (lobe/port alignment reported, not asserted -- see "
      "mmi_1x4_sin300_validation.md step 6d's open port-imbalance question)")
