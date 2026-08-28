"""Tapered 1x4, WAIST profile: both body faces at mm_width (matching the port
fan-out exactly, no compromise), narrowing to a center width for an interior
run, then widening back out. Sweeps mm_length per profile to find the design
point, and REUSES mode solves across the down/up transitions by symmetry.

validation_markdown/mmi_1x4_tapered_body_validation.md step 6.

WHY THIS PROFILE, REPLACING THE FACE-TAPER ONE. The face-taper design
(body_start_width/body_end_width) hit a real geometric wall: the 1x4's 4-way
fan-out needs 11.5 of the 12.0 um output face, leaving only 4.2% of room to
narrow there, and the sweep of that design (mmi_1x4_taper_profile_sweep.py,
job 12383125/12387482) found narrowing the INPUT face instead didn't save
length either -- at the ARGUED-optimal length the device lost 44% of its
power (imaging failure, not radiation), and the best point found was LONGER
than the straight device, not shorter.

The waist profile sidesteps the port-matching problem entirely: BOTH faces
sit at mm_width, identical to the straight body at the ports, so there is
no fan-out compromise and no re-derivation of the port-taper design needed.
The narrowing happens in the interior, where self-imaging length scales
with width (L ~ W^2), so a genuinely narrower interior section can in
principle accumulate the same imaging phase over less length.

GEOMETRY. `MMI(body_waist_width=Ww, body_waist_length=Lw,
body_transition_length=Lt)` -- see beef/geometry/mmi.py's own docstring.
Symmetric about the body centre:

    |<-stub->|<-Lt->|<--- Lw (width Ww) --->|<-Lt->|<-stub->|
      W                W->Ww                  Ww->W           W

stub is whatever is left over: (mm_length - 2*Lt - Lw) / 2.

MODE-SOLVE REUSE, why it is exact and not an approximation. The down
transition (W->Ww over Lt) and the up transition (Ww->W over Lt) are, cell
for cell, the SAME sequence of cross-section WIDTHS in reverse order -- pure
geometric mirror symmetry, not a numerical coincidence. A meow eigenmode
solve depends only on the LOCAL 2D dielectric cross-section (width, material,
wavelength, mesh) -- position along the propagation axis never enters
compute_modes()/CrossSection.from_cell() at all. So a cell in the up
transition with the same width as an already-solved down-transition cell
gets the IDENTICAL Mode objects, not merely similar ones. This is enforced
here via a width-keyed cache (round to WIDTH_KEY_DECIMALS) rather than by
hand-reversing index lists, which also transparently reuses across the two
end stubs (both at mm_width) and is robust to whatever cell counts are
chosen, not dependent on the two transitions having matching discretisation.
Reported at the end: how many solves were actually avoided.

THIS BYPASSES EMESolver.solve()/_build_cells() DELIBERATELY. _build_cells()
divides each region's length EQUALLY by that region's own num_cells -- it has
no notion of the body's five internal segments (stub/transition/waist/
transition/stub), so it cannot place cell boundaries exactly at
din0/din1/dout0/dout1. This script instead builds the body's own cell-length
list directly (five segments, each internally uniform) and calls
mw.create_cells() itself, exactly mirroring what _build_cells()/solve() do
for the non-body regions (same geometry.to_gds_eme() / extrude_polygons() /
self._build_mesh() calls) so nothing else about the whole-device solve
changes.

COST NOTE. Deliberately at a DESIGN resolution (res=30, 16x4 domain), not
production (40, 20x5) -- see the face-taper sweep's own note; re-run any
winner at production settings before trusting its numbers.

Run:
    MMI1X4_WP_WW=6.0 MMI1X4_WP_LW=10.0 MMI1X4_WP_LT=8.0 MMI1X4_WP_NLEN=3 \\
        python .../mmi_1x4_waist_profile_sweep.py
    python .../mmi_1x4_waist_profile_sweep.py
"""
import os
import time

import numpy as np
import meow as mw

from beef import MMI, EMESolver, SIN_300, extrude_polygons

# --- production constants (step 9c) ----------------------------------------
WL = 1.55
WAVEGUIDE_WIDTH, TAPER_WIDTH, TAPER_LENGTH = 1.1, 2.2, 5.0
WAVEGUIDE_OFFSET, MM_WIDTH = 3.1, 12.0
L_STRAIGHT = 39.90
STRAIGHT_LENGTH = 0.5
N_BODY, N_TAPER, C_TAPER = 20, 10, 8

# --- sweep config -------------------------------------------------------
RES = float(os.environ.get("MMI1X4_WP_RES", 30.0))
OXW = float(os.environ.get("MMI1X4_WP_OXW", MM_WIDTH + 4.0))
OXH = float(os.environ.get("MMI1X4_WP_OXH", 4.0))
_ww = os.environ.get("MMI1X4_WP_WW")
_lw = os.environ.get("MMI1X4_WP_LW")
_lt = os.environ.get("MMI1X4_WP_LT")
WW_LIST = [float(v) for v in _ww.split(",")] if _ww else [8.0, 6.0, 4.0]
LW_LIST = [float(v) for v in _lw.split(",")] if _lw else [10.0]
LT_LIST = [float(v) for v in _lt.split(",")] if _lt else [8.0]
C_STUB = int(os.environ.get("MMI1X4_WP_CSTUB", 1))
C_TRANS = int(os.environ.get("MMI1X4_WP_CTRANS", 6))
C_WAIST = int(os.environ.get("MMI1X4_WP_CWAIST", 4))
N_LEN = int(os.environ.get("MMI1X4_WP_NLEN", 5))
LEN_FRAC = float(os.environ.get("MMI1X4_WP_LENFRAC", 0.15))
L_CENTRE = os.environ.get("MMI1X4_WP_LCENTRE")
WIDTH_KEY_DECIMALS = 6


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
    """The body's own 5-segment cell-length list, aligned exactly to
    stub/transition/waist/transition/stub, plus each cell's own width
    (evaluated in the device frame via body_width_at at the cell centre) --
    the width list is what the reuse cache below keys on."""
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


def solve_reuse(mmi, output_idx, solver, mesh):
    """Whole-device single-wavelength solve with width-keyed mode reuse,
    replicating EMESolver._build_cells()/solve() exactly except for the
    body region's cell construction (see module docstring)."""
    body_lens, body_widths = body_cell_lengths_and_widths(mmi)

    in_len = STRAIGHT_LENGTH / 1
    taper_in_len = TAPER_LENGTH / C_TAPER
    taper_out_len = TAPER_LENGTH / C_TAPER
    out_len = STRAIGHT_LENGTH / 1

    Ls = ([in_len] + [taper_in_len] * C_TAPER + body_lens
          + [taper_out_len] * C_TAPER + [out_len])
    modes_per_cell = ([1] + [N_TAPER] * C_TAPER + [N_BODY] * len(body_lens)
                      + [N_TAPER] * C_TAPER + [1])

    gds_polys = mmi.to_gds_eme(
        straight_length=solver.straight_length, oxide_width=solver.oxide_width,
        x_start=0.0, side=f"out{output_idx}", output_idx=output_idx,
    )
    structs = extrude_polygons(gds_polys, solver._extrusion_rules(np.array([WL])))
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=Ls)

    env = mw.Environment(wl=float(WL), T=solver.temperature)
    css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]

    # Width-keyed reuse cache. Every cell outside the body region gets its
    # own unique key (widths list only covers the body), so only body cells
    # -- specifically the mirror-symmetric transition/stub pairs -- can hit.
    reuse_cache = {}
    solves, hits = 0, 0
    modes = []
    body_start = 1 + C_TAPER
    body_end = body_start + len(body_lens)
    for i, (cs, nm) in enumerate(zip(css, modes_per_cell)):
        key = None
        if body_start <= i < body_end:
            key = round(body_widths[i - body_start], WIDTH_KEY_DECIMALS)
        if key is not None and key in reuse_cache:
            modes.append(reuse_cache[key])
            hits += 1
        else:
            m = solver.compute_modes(cs, nm)
            modes.append(m)
            solves += 1
            if key is not None:
                reuse_cache[key] = m

    S, port_map = mw.compute_s_matrix(modes, cells=cells)
    amp = S[port_map["right@0"], port_map["left@0"]]
    return complex(amp), solves, hits


def transmission_all_outputs(mmi):
    solver = EMESolver(
        num_modes=[1], num_cells=[1], lengths=[1.0],   # placeholder, unused below
        oxide_width=OXW, oxide_height=OXH, res=RES,
        material_stack=SIN_300, mode_solver="tidy3d",
    )
    solver.straight_length = STRAIGHT_LENGTH
    mesh = solver._build_mesh()
    T = np.zeros(4)
    tot_solves = tot_hits = 0
    for k in range(4):
        amp, solves, hits = solve_reuse(mmi, k, solver, mesh)
        T[k] = abs(amp) ** 2
        tot_solves += solves
        tot_hits += hits
    return T, tot_solves, tot_hits


print("=" * 84)
print("TAPERED 1x4 -- WAIST profile sweep (both faces at mm_width, interior narrowing)")
print(f"  res={RES}  domain {OXW:.0f}x{OXH:.0f}  C_STUB={C_STUB} C_TRANS={C_TRANS} "
      f"C_WAIST={C_WAIST}  wl={WL}")
print(f"  Ww list {WW_LIST}   Lw list {LW_LIST}   Lt list {LT_LIST}")
print("=" * 84, flush=True)

t0 = time.time()
Ts, sv, ht = transmission_all_outputs(build(L_STRAIGHT))
bal_s = Ts.max() - Ts.min()
print(f"\nSTRAIGHT control  L={L_STRAIGHT:.2f}  totalT={Ts.sum():.4f}  "
      f"balance={bal_s:.2e}  ({time.time()-t0:.0f}s, {sv} solves, {ht} reused)")
print(f"  T = {list(np.round(Ts, 5))}", flush=True)

rows = []
for ww in WW_LIST:
    for lw in LW_LIST:
        for lt in LT_LIST:
            centre = float(L_CENTRE) if L_CENTRE is not None else L_STRAIGHT
            # np.linspace(a, b, 1) returns [a], NOT the midpoint -- with
            # N_LEN=1 that silently swept the low window edge instead of
            # `centre` (caught in smoke job 12394703: reported a catastrophic
            # totalT=0.346 at L=33.91 that was actually testing the WRONG
            # length, not a real profile failure). N_LEN=1 now means
            # "evaluate exactly at centre", matching what a caller asking
            # for one point actually wants.
            if N_LEN == 1:
                lens = np.array([centre])
            else:
                lens = np.linspace(centre * (1 - LEN_FRAC), centre * (1 + LEN_FRAC), N_LEN)
            print(f"\n-- Ww={ww:.1f} Lw={lw:.1f} Lt={lt:.1f}  centre={centre:.2f}  "
                  f"sweeping {lens[0]:.2f}..{lens[-1]:.2f}", flush=True)
            best = None
            for L in lens:
                if 2 * lt + lw > L - 1e-6:
                    print(f"     L={L:6.2f}  skipped (no room for stubs)")
                    continue
                t0 = time.time()
                try:
                    T, sv, ht = transmission_all_outputs(build(L, ww=ww, lw=lw, lt=lt))
                except Exception as exc:
                    print(f"     L={L:6.2f}  FAILED {type(exc).__name__}: {exc}", flush=True)
                    continue
                bal = T.max() - T.min()
                print(f"     L={L:6.2f}  totalT={T.sum():.4f}  balance={bal:.2e}  "
                      f"({time.time()-t0:.0f}s, {sv} solves, {ht} reused)", flush=True)
                if best is None or bal < best[2]:
                    best = (L, T.sum(), bal, T)
            if best is None:
                continue
            L, tot, bal, T = best
            edge = abs(L - lens[0]) < 1e-9 or abs(L - lens[-1]) < 1e-9
            rows.append((ww, lw, lt, L, tot, bal, edge))
            print(f"   BEST L={L:.2f}  totalT={tot:.4f}  balance={bal:.2e}"
                  + ("   <-- ON WINDOW EDGE" if edge else ""), flush=True)

print("\n" + "=" * 84)
print(f"{'Ww':>5} {'Lw':>5} {'Lt':>5} {'L_best':>7} {'totalT':>8} {'balance':>9} "
      f"{'saving':>9}")
print(f"{'--':>5} {'--':>5} {'--':>5} {L_STRAIGHT:7.2f} {Ts.sum():8.4f} {bal_s:9.2e} "
      f"{'(control)':>9}")
for ww, lw, lt, L, tot, bal, edge in rows:
    sav = L_STRAIGHT - L
    print(f"{ww:5.1f} {lw:5.1f} {lt:5.1f} {L:7.2f} {tot:8.4f} {bal:9.2e} "
          f"{sav:6.2f}um" + ("  *edge" if edge else ""))
print("=" * 84)
print("A profile only wins if it is SHORTER *and* holds totalT/balance against")
print("the straight control. Re-run any winner at production settings before quoting.")
