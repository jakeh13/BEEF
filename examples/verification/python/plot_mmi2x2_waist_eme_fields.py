"""Step 7 field plots for the WAISTED 2x2 -- straight control vs design point.

validation_markdown/mmi_2x2_waist_body_validation.md step 7. Deliverable:
straight-vs-waisted EME field plots from ONE config, so the only difference
between the two renders is the waist itself.

DESIGN POINT (user's call, 2026-09-01: the MEAN phase reading governs, so the
shortest passer under that reading wins):

    Ww=2.24, Lw=8.0, Lt=6.0, L=20.5  (stub s=0.25) -> 29.50 um total
    IL 0.3910 dB | imb 0.0852 | mean phase dev 1.696 deg
    direct-solve confirmed, job 12687092, dIL=-1.05e-14 vs the analytic grid

    Recorded openly: this point's MAX-over-band phase deviation is 2.244 deg,
    which exceeds the 2.0 deg tolerance under the max reading of the same
    pre-registered text (job 12687020). The user chose the mean reading, which
    is what score() has always applied. Set MMI2X2_WF_WW/LW/LT/L to plot the
    strict-reading survivor (2.24/10.0/5.0/20.5, also 29.50 um) instead.

WHY THIS DOES NOT COPY THE 1x4's TEMPLATE. plot_mmi1x4_waist_profile_eme_fields.py
bypasses _build_cells() and hand-rolls the body segmentation, because the 1x4
sweep declares the body as ONE region and _build_cells() would divide it
uniformly, unable to place boundaries at the stub/transition/waist edges.

The 2x2 free sweep already solved that differently: it declares the body as
NINE regions, each uniform or a single linear ramp. _build_cells() divides each
region evenly, which for those regions is exactly right. So this script goes
through the normal solver path and reuses the sweep's own regions() layout --
the plotted geometry is the swept geometry by construction, not a copy that can
drift out of sync.

READING THE PLOTS -- plot_eme_fields() labels in MEOW's frame, not the
device's: device Ex = meow Ez, device Ey = meow Ex, device Ez = meow Ey. The
panel titled |Ex| is the device's in-plane TRANSVERSE field and its "x [um]"
axis is device y. meow's Ex is the dominant TE component here.

Run (smoke first, same discipline as every expensive job in this repo):
    MMI2X2_WF_RES=15 python .../plot_mmi2x2_waist_eme_fields.py
    python .../plot_mmi2x2_waist_eme_fields.py          # real config
"""
import os
import time

import numpy as np
import meow as mw

from beef import MMI, EMESolver, SOI_220, plot_eme_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "verification",
    "mmi_2x2_waist_fields"))
os.makedirs(output_dir, exist_ok=True)


def _env_f(n, d):
    return float(os.environ.get(n, d))


def _env_i(n, d):
    return int(os.environ.get(n, d))


# --- frozen port constants, identical to mmi_2x2_waist_free_sweep.py --------
MM_LENGTH_STD, MM_WIDTH = 28.0, 2.8
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 0.45, 1.3, 4.0
WAVEGUIDE_OFFSET, WAVEGUIDE_HEIGHT = 1.5, 0.22
STRAIGHT_LENGTH = 0.5
TAPER_CELLS = 15
N_PORT_MODES, N_TAPER_MODES, N_BODY_MODES = 1, 7, 15
OXW, OXH = 5.0, 1.5

WL = _env_f("MMI2X2_WF_WL", 1.5475)          # band centre of 1.530-1.565
RES = _env_f("MMI2X2_WF_RES", 40.0)
C_TRANS = _env_i("MMI2X2_WF_CTRANS", 8)      # step 5a's converged floor
WW = _env_f("MMI2X2_WF_WW", 2.24)
LW = _env_f("MMI2X2_WF_LW", 8.0)
LT = _env_f("MMI2X2_WF_LT", 6.0)
L_WAIST = _env_f("MMI2X2_WF_L", 20.5)
PORT = os.environ.get("MMI2X2_WF_PORT", "bar")

PORTS = {"bar": (1, 1), "cross": (1, 0)}
if PORT not in PORTS:
    raise SystemExit(f"MMI2X2_WF_PORT must be one of {sorted(PORTS)}")
IIDX, OIDX = PORTS[PORT]

S_STUB = (L_WAIST - 2 * LT - LW) / 2
if S_STUB < 0.25 - 1e-9:
    raise SystemExit(f"stub {S_STUB:.3f} < 0.25 um: this geometry would change "
                     "the hybrid's FDTD sections, which are frozen")

print(f"wl={WL} res={RES} port={PORT} (in {IIDX} -> out {OIDX})")
print(f"waisted: Ww={WW} Lw={LW} Lt={LT} L={L_WAIST} (stub {S_STUB:.3f}) "
      f"-> {L_WAIST + 2*TAPER_LENGTH + 2*STRAIGHT_LENGTH:.2f} um total")
print(f"control: straight L={MM_LENGTH_STD} "
      f"-> {MM_LENGTH_STD + 2*TAPER_LENGTH + 2*STRAIGHT_LENGTH:.2f} um total\n")


def waisted_regions(lw, lt, s):
    """The sweep's 9-region layout: uniform segments get ONE cell because they
    ARE uniform; only the two linear ramps need C_TRANS."""
    return (
        [STRAIGHT_LENGTH, TAPER_LENGTH, s, lt, lw, lt, s, TAPER_LENGTH,
         STRAIGHT_LENGTH],
        [1, TAPER_CELLS, 1, C_TRANS, 1, C_TRANS, 1, TAPER_CELLS, 1],
        [N_PORT_MODES, N_TAPER_MODES, N_BODY_MODES, N_BODY_MODES,
         N_BODY_MODES, N_BODY_MODES, N_BODY_MODES, N_TAPER_MODES,
         N_PORT_MODES],
    )


def straight_regions(mm_len):
    """Control: the body is one uniform region. mw.propagate_modes samples
    analytically WITHIN a cell, so one cell renders the field correctly -- the
    500-point z grid in plot_eme_fields() does the spatial sampling."""
    return (
        [STRAIGHT_LENGTH, TAPER_LENGTH, mm_len, TAPER_LENGTH, STRAIGHT_LENGTH],
        [1, TAPER_CELLS, 1, TAPER_CELLS, 1],
        [N_PORT_MODES, N_TAPER_MODES, N_BODY_MODES, N_TAPER_MODES,
         N_PORT_MODES],
    )


def build_mmi(mm_len, waisted):
    kw = dict(body_waist_width=WW, body_waist_length=LW,
              body_transition_length=LT) if waisted else {}
    return MMI(n_inputs=2, n_outputs=2, mm_length=mm_len, mm_width=MM_WIDTH,
               waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
               taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
               waveguide_height=WAVEGUIDE_HEIGHT, **kw)


CASES = [
    ("straight", MM_LENGTH_STD, False, straight_regions(MM_LENGTH_STD)),
    ("waisted", L_WAIST, True, waisted_regions(LW, LT, S_STUB)),
]

wavelengths = np.array([WL])
for name, mm_len, waisted, (lens, cells_n, modes_n) in CASES:
    t0 = time.time()
    solver = EMESolver(num_modes=modes_n, num_cells=cells_n, lengths=lens,
                       oxide_width=OXW, oxide_height=OXH, res=RES,
                       material_stack=SOI_220, mode_solver="tidy3d")
    mmi = build_mmi(mm_len, waisted)
    half = mmi.mm_length / 2
    print(f"[{name}] body width: in-face={mmi.width_at(-half):.3f}  "
          f"mid={mmi.width_at(0.0):.3f}  out-face={mmi.width_at(half):.3f} um")

    cells = solver._build_cells(mmi, side=PORT, wavelengths=wavelengths,
                                input_idx=IIDX, output_idx=OIDX)
    mpc = solver._modes_per_cell()
    assert len(cells) == len(mpc), (len(cells), len(mpc))

    env = mw.Environment(wl=float(WL), T=solver.temperature)
    css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
    print(f"[{name}] solving {len(cells)} cells, {sum(mpc)} mode-solves ...",
          flush=True)
    modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, mpc)]

    tag = f"{name}_{PORT}_{WL:.4f}um"
    plot_eme_fields(cells, modes, output_dir, tag=tag, input_mode=1)
    print(f"[{name}] done in {time.time()-t0:.0f}s -> "
          f"{output_dir}/fields/eme_*_{tag}.png\n", flush=True)

print(f"Both cases under {output_dir}/fields/")
print("Compare eme_Ex_straight_* against eme_Ex_waisted_* -- meow's Ex is the")
print("dominant TE component. The waisted render is 7.5 um shorter; look at")
print("whether the self-imaging refocus survives the two ramps intact, and")
print("whether anything is shed AT the ramps rather than through the waist.")
