"""Field plots for the WAIST profile's failing points -- straight vs waisted,
one output branch each, at a specific (Ww, Lw, Lt, L) design point.

validation_markdown/mmi_1x4_tapered_body_validation.md -- the visual
diagnostic the waist-profile closure note calls for but was not run: is the
measured totalT loss concentrated at the two width-transition interfaces
(implicates local-mode-basis mismatch/scattering there) or spread through
the whole waist region (implicates a self-imaging beat-length
desynchronization across the interior, not a localized transition loss)?

Default point is Ww=9.0/Lw=10.0/Lt=8.0/L=36.91 -- job 12436352's best-of-a-
bad-set result (genuinely 2.99um shorter than the L=39.90 control, but still
totalT=0.6186 vs control's 0.9703). Chosen because it is the LEAST-bad
failure on record, so any localized loss mechanism should be easiest to see
here rather than in the more catastrophic Ww=6.0/10.0 points.

GEOMETRY/CELL-BUILDING deliberately mirrors `mmi_1x4_waist_profile_sweep.py`
exactly (same production constants, same `body_cell_lengths_and_widths`
segment logic) -- copied rather than imported, because that script's
top-level code has no `if __name__ == "__main__"` guard and runs its full
sweep as an import side effect. Keep the two files' body-segment logic in
sync if either changes; there is no other coupling between them.

THIS BYPASSES EMESolver.solve()/_build_cells() DELIBERATELY, same reason as
the sweep script: `_build_cells()` divides each region UNIFORMLY by its cell
count and cannot place cell boundaries at the waist's internal
stub/transition/waist edges.

Run (smoke first, same discipline as the straight-vs-tapered field script):
    MMI1X4_WPF_CTRANS=3 MMI1X4_WPF_CWAIST=2 MMI1X4_WPF_RES=15 \\
        python .../plot_mmi1x4_waist_profile_eme_fields.py
    python .../plot_mmi1x4_waist_profile_eme_fields.py   # real config
"""
import os
import time

import numpy as np
import meow as mw

from beef import MMI, EMESolver, SIN_300, extrude_polygons, plot_eme_fields

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_1x4",
    "mmi_1x4_waist_profile_eme"))
os.makedirs(output_dir, exist_ok=True)


def _env_f(name, default):
    return float(os.environ.get(name, default))


def _env_i(name, default):
    return int(os.environ.get(name, default))


# --- production constants (step 9c), same as mmi_1x4_waist_profile_sweep.py -
WL = 1.55
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 1.1, 2.2, 5.0
WAVEGUIDE_OFFSET, MM_WIDTH = 3.1, 12.0
L_STRAIGHT = 39.90
STRAIGHT_LENGTH = 0.5
N_BODY, N_TAPER, C_TAPER = 20, 10, 8

RES = _env_f("MMI1X4_WPF_RES", 30.0)
OXW = _env_f("MMI1X4_WPF_OXW", MM_WIDTH + 4.0)
OXH = _env_f("MMI1X4_WPF_OXH", 4.0)
WW = _env_f("MMI1X4_WPF_WW", 9.0)
LW = _env_f("MMI1X4_WPF_LW", 10.0)
LT = _env_f("MMI1X4_WPF_LT", 8.0)
L_WAIST = _env_f("MMI1X4_WPF_L", 36.91)          # job 12436352's Ww=9.0 best
C_STUB = _env_i("MMI1X4_WPF_CSTUB", 1)
C_TRANS = _env_i("MMI1X4_WPF_CTRANS", 6)
C_WAIST = _env_i("MMI1X4_WPF_CWAIST", 4)
OUTPUT_IDX = _env_i("MMI1X4_WPF_OUT", 0)

print(f"res={RES}  Ww={WW}  Lw={LW}  Lt={LT}  L={L_WAIST}  "
      f"C_STUB={C_STUB} C_TRANS={C_TRANS} C_WAIST={C_WAIST}")
print("Diagnostic: is loss localized at the two width-transition interfaces,")
print("or spread through the whole waist region?\n")


def build(mm_length, ww=None, lw=0.0, lt=0.0):
    kw = {}
    if lt > 0 and ww is not None:
        kw = dict(body_waist_width=ww, body_waist_length=lw,
                  body_transition_length=lt)
    return MMI(
        n_inputs=1, n_outputs=4, mm_length=mm_length, mm_width=MM_WIDTH,
        waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
        taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
        waveguide_height=0.30, use_neighbor_aperture=True, **kw
    )


def body_cell_lengths_and_widths(mmi):
    """Same 5-segment (stub/transition/waist/transition/stub) logic as
    mmi_1x4_waist_profile_sweep.py's function of the same name -- see that
    module's docstring for why this can't go through _build_cells()."""
    half = mmi.mm_length / 2
    if mmi.is_body_waisted:
        (din0, din1), (dout0, dout1) = mmi._body_waist_bounds()
        stub = din0 - (-half)
        segs = [(-half, din0, C_STUB), (din0, din1, C_TRANS),
                (din1, dout0, C_WAIST), (dout0, dout1, C_TRANS),
                (dout1, half, C_STUB)]
    else:
        segs = [(-half, half, max(1, C_STUB + C_TRANS + C_WAIST))]
    lens, widths = [], []
    for x0, x1, n in segs:
        if n <= 0 or x1 - x0 <= 1e-12:
            continue
        seg_len = (x1 - x0) / n
        for i in range(n):
            xc = x0 + (i + 0.5) * seg_len
            lens.append(seg_len)
            widths.append(mmi.body_width_at(xc))
    return lens, widths


CASES = [
    ("straight", L_STRAIGHT, dict()),
    ("waisted", L_WAIST, dict(ww=WW, lw=LW, lt=LT)),
]

base_solver = EMESolver(
    num_modes=[1], num_cells=[1], lengths=[1.0],
    oxide_width=OXW, oxide_height=OXH, res=RES,
    material_stack=SIN_300, mode_solver="tidy3d",
)
base_solver.straight_length = STRAIGHT_LENGTH
mesh = base_solver._build_mesh()

for name, mm_length, kw in CASES:
    t0 = time.time()
    mmi = build(mm_length, **kw)
    half = mmi.mm_length / 2
    print(f"[{name}] L={mm_length:.2f}  body width: in-face={mmi.width_at(-half):.3f}  "
          f"mid={mmi.width_at(0.0):.3f}  out-face={mmi.width_at(half):.3f} um")

    body_lens, body_widths = body_cell_lengths_and_widths(mmi)

    in_len = STRAIGHT_LENGTH
    taper_in_len = TAPER_LENGTH / C_TAPER
    taper_out_len = TAPER_LENGTH / C_TAPER
    out_len = STRAIGHT_LENGTH

    Ls = ([in_len] + [taper_in_len] * C_TAPER + body_lens
          + [taper_out_len] * C_TAPER + [out_len])
    modes_per_cell = ([1] + [N_TAPER] * C_TAPER + [N_BODY] * len(body_lens)
                      + [N_TAPER] * C_TAPER + [1])

    gds_polys = mmi.to_gds_eme(
        straight_length=base_solver.straight_length, oxide_width=base_solver.oxide_width,
        x_start=0.0, side=f"out{OUTPUT_IDX}", output_idx=OUTPUT_IDX,
    )
    structs = extrude_polygons(gds_polys, base_solver._extrusion_rules(np.array([WL])))
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=Ls)

    env = mw.Environment(wl=float(WL), T=base_solver.temperature)
    css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]

    print(f"[{name}] solving {len(cells)} cells ({sum(modes_per_cell)} mode-solves, "
          f"no reuse cache -- this is a one-shot plot, not the swept metric)...")
    modes = [base_solver.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]

    tag = f"{name}_out{OUTPUT_IDX}_{WL:.3f}um"
    plot_eme_fields(cells, modes, output_dir, tag=tag, input_mode=1)
    print(f"[{name}] done in {time.time() - t0:.1f}s -> "
          f"{output_dir}/fields/eme_*_{tag}.png\n")

print(f"Both cases written under {output_dir}/fields/")
print("Compare eme_Ex_straight_*.png against eme_Ex_waisted_*.png -- meow's Ex "
      "is the dominant TE component for this device (measured on the "
      "straight-vs-face-taper field script, |Ex| ~34 vs |Ey| ~0.6).")
