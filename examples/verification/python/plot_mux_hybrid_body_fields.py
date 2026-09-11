"""Modal multiplexer — hybrid body field reconstruction, device frame.

⚠⚠ WRITTEN BUT NEVER RUN — 2026-08-31. Deferred deliberately, not abandoned.
Imports and geometry construction are verified; the reconstruction itself has
NOT been executed once, so no number this produces has been checked against
anything. Treat every output as unvalidated until someone runs it and confirms
(a) captured power <= 1, (b) the near-degenerate gap report against this
device's own connector Delta-n_eff=0.0148, and (c) that MUX_SECTION bounds
match modal_multiplexer_hybrid.py's `boundaries` arithmetic after any change
there. The four sibling scripts (1x2, 1x4, DC x2, HCG) ARE validated; this one
is the exception.

Pinned to modal_multiplexer/python/modal_multiplexer_hybrid.py's own body
sections. See beef/plot/hybrid_body.py for why plot_eme_fields() is not
sufficient here, and for the standing caveat that the excitation is an
APPROXIMATION rather than the hybrid's own interface field.

THIS DEVICE HAS THREE EME SECTIONS, NOT ONE. Its cascade is seven sections
(interface_modes has six entries):

    FDTD in | EME MMI 1 | FDTD taper-out | EME connector |
            FDTD taper-in | EME MMI 2 | FDTD out

so MUX_SECTION selects which EME body to reconstruct:

    MUX_SECTION=mmi1        MMI 1 body, N_MODES=15   (default)
    MUX_SECTION=mmi2        MMI 2 body, N_MODES=15
    MUX_SECTION=connector   the three-waveguide run, CONN_MODES=6

The connector is a genuinely different object from the two MMI bodies -- three
separate weakly-coupled strips rather than one wide slab -- and it is the
section whose near-degenerate modes (Delta-n_eff=0.0148 between modes 1-2)
triggered the 2026-08-09 excite-normalization fix in the first place. Worth
looking at for that reason alone.

⚠ WAVELENGTH. Defaults to 1.31 um, the band CENTRE, deliberately.
`modal_multiplexer_validation.md` records a real mode-identity /
basis-completeness gap at 1.26 um, the band edge, which is an unresolved
design tradeoff on this device rather than a bug to plot around. Reconstructing
at the edge is legitimate and interesting, but the basis there is known
incomplete, so any picture from it must be read with that doc open -- hence the
explicit override rather than a silent default.

Cost: one cross-section. MMI bodies carry 15 modes on a 2.8 um slab; the
connector carries 6.

Run:
    python .../plot_mux_hybrid_body_fields.py
    MUX_SECTION=connector python .../plot_mux_hybrid_body_fields.py
    MUX_SECTION=mmi2 MUX_EXCITE_MODE=1 python .../plot_mux_hybrid_body_fields.py
"""
import os

import numpy as np

from beef import (MMI, MMIPort, WaveguideRun, CompoundGeometry, EMESolver,
                  SOI_160, plot_hybrid_body_fields)

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "modal_multiplexer",
    "modal_multiplexer_hybrid_body"))
os.makedirs(output_dir, exist_ok=True)

# --- pinned to modal_multiplexer/python/modal_multiplexer_hybrid.py ---------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.26, 1.36, 11
TAPER_L = 3.0
WO = 1.0375                  # centre-to-centre pitch of the three branches
BETWEEN_WG_LENGTH = 3.0
MM_WIDTH = 2.8
BODY_STUB = 0.75
CONN_STUB = 0.25
N_MODES = 15                 # body supermodes at the four MMI-body interfaces
CONN_MODES = 6               # arm supermodes at the two connector interfaces

SECTION = os.environ.get("MUX_SECTION", "mmi1")
WL_PLOT = float(os.environ.get("MUX_WL_PLOT", 1.31))   # band centre -- see docstring
EXCITE_MODE = os.environ.get("MUX_EXCITE_MODE")        # None -> gaussian launch

BRANCH_SPECS = [(+WO, 0.33, 0.725), (0.0, 0.35, 1.15), (-WO, 0.35, 0.725)]


def _branch_ports():
    return [MMIPort(y=y, waveguide_width=ww, taper_length=TAPER_L, taper_width=tw)
            for y, ww, tw in BRANCH_SPECS]


mmi1 = MMI(
    n_inputs=1, n_outputs=3, mm_length=16.92, mm_width=MM_WIDTH,
    inputs=[MMIPort(y=0.0, waveguide_width=0.7, taper_length=TAPER_L,
                    taper_width=1.4)],
    outputs=_branch_ports(), waveguide_height=0.16,
)
mmi2 = MMI(
    n_inputs=3, n_outputs=3, mm_length=34.24, mm_width=MM_WIDTH,
    inputs=_branch_ports(), outputs=_branch_ports(),
    waveguide_height=0.16, use_neighbor_aperture=True,
)
connector = WaveguideRun(
    length=BETWEEN_WG_LENGTH, ports=[(y, ww) for y, ww, _ in BRANCH_SPECS])
mux = CompoundGeometry([mmi1, connector, mmi2])

# Section bounds, replicating modal_multiplexer_hybrid.py's own arithmetic.
half_total = mux.total_length / 2
body1_start = -half_total + TAPER_L
body1_end = body1_start + mmi1.mm_length
conn_start = body1_end + TAPER_L
conn_end = conn_start + BETWEEN_WG_LENGTH
body2_start = conn_end + TAPER_L
body2_end = body2_start + mmi2.mm_length

if SECTION == "mmi1":
    x0, x1 = body1_start + BODY_STUB, body1_end - BODY_STUB
    n_modes, launch_fwhm, launch_y = N_MODES, 1.4, 0.0
    half_w, ports = MM_WIDTH / 2, [y for y, _, _ in BRANCH_SPECS]
    label = "Modal multiplexer — MMI 1 body"
elif SECTION == "mmi2":
    x0, x1 = body2_start + BODY_STUB, body2_end - BODY_STUB
    n_modes, launch_fwhm, launch_y = N_MODES, 1.15, 0.0
    half_w, ports = MM_WIDTH / 2, [y for y, _, _ in BRANCH_SPECS]
    label = "Modal multiplexer — MMI 2 body"
elif SECTION == "connector":
    x0, x1 = conn_start + CONN_STUB, conn_end - CONN_STUB
    n_modes, launch_fwhm, launch_y = CONN_MODES, 0.35, 0.0
    half_w, ports = WO + 0.35, [y for y, _, _ in BRANCH_SPECS]
    label = "Modal multiplexer — connector"
else:
    raise SystemExit(f"unknown MUX_SECTION={SECTION!r}, "
                     "expected mmi1|mmi2|connector")

solver = EMESolver(
    num_modes=[n_modes], num_cells=[1], lengths=[1.0],
    oxide_width=5.0, oxide_height=1.5,
    res=RESOLUTION, material_stack=SOI_160, mode_solver="meep",
)

print(f"mux {SECTION}: x in [{x0:+.2f}, {x1:+.2f}] um (L={x1 - x0:.2f} of "
      f"{mux.total_length:.2f} um device), {n_modes} modes @ {WL_PLOT} um")
print(f"branch centres at y = {[round(p, 4) for p in ports]} um")

kw = dict(excitation="gaussian", excitation_center=launch_y,
          excitation_fwhm=launch_fwhm)
if EXCITE_MODE is not None:
    kw = dict(excitation="mode", excitation_mode=int(EXCITE_MODE))
    print(f"exciting body eigenmode {int(EXCITE_MODE)} "
          "(a single eigenmode is invariant along a uniform section -- expect "
          "a flat map; use it to inspect one mode's transverse profile)")

res = plot_hybrid_body_fields(
    mux, solver, x0, x1, WL_PLOT, n_modes, output_dir,
    tag=f"mux_{SECTION}_{WL_PLOT:.3f}um"
        + (f"_mode{EXCITE_MODE}" if EXCITE_MODE is not None else ""),
    wavelengths=np.linspace(WL_MIN, WL_MAX, N_WL),
    port_ys=ports, body_half_width=half_w,
    device_label=label, total_length=mux.total_length, **kw)

w, neff = res["weights"], res["neff"]
assert res["captured_power"] <= 1.0 + 1e-6, \
    f"captured power {res['captured_power']:.4f} > 1 -- wrong inner product"

# Near-degenerate mode pairs, reported because THIS device is where they
# matter: Delta-n_eff=0.0148 between the connector's modes 1-2 is what
# triggered the 2026-08-09 solve_section() excite-normalization fallback
# (_EXCITE_NORM_MISMATCH_TOL=0.1). Printing the spacing here makes that
# regime visible before it shows up as a transfer-matrix anomaly.
gaps = np.abs(np.diff(neff))
print(f"\nsmallest consecutive n_eff gap: {gaps.min():.6f} "
      f"(modes {int(gaps.argmin())}-{int(gaps.argmin()) + 1})")
if gaps.min() < 0.02:
    print("  ^ near-degenerate: this is the regime that fired the 2026-08-09 "
          "excite-normalization fallback on this device's connector section")
print(f"basis off-diagonal: {res['gram_max_offdiag']:.4f}")
print("\nchecks passed")
