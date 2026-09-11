"""Waisted 2x2 -- stitched hybrid field: FDTD input section + a REAL
multi-cell EME body, driven by the captured interface field.

Ported from stitch_mmi2x2_hybrid_fields.py (the straight device), which
CANNOT be reused as-is: that script propagates the body as a single
diagonal exp(2j*pi*neff*x/wl) per mode, exact only for a UNIFORM
cross-section, and this body is stub|ramp|waist|ramp|stub -- not uniform.

WHAT'S DIFFERENT HERE, and why it is not "new machinery" so much as reused
machinery pointed at a new input:

- The body is built as ONE continuous multi-cell mesh (mw.create_cells over
  the whole [x_b0, x_b1] span, sliced into stub(1)|ramp(C_TRANS)|waist(1)|
  ramp(C_TRANS)|stub(1) cells) -- the same region-decomposition trick
  plot_mmi2x2_waist_eme_fields.py already uses for a WHOLE-DEVICE EME-only
  render, here restricted to just the body span so the ports stay real FDTD.
- Instead of meow's propagate_modes() convenience (which reconstructs a
  device-frame-AMBIGUOUS slice -- its mesh_x/component labels are meow's own
  local axes, not beef's device frame, see stack_meow_mode_fields()'s
  docstring), this script calls the SAME lower-level cascade meow uses
  internally (compute_propagation_s_matrices, compute_interface_s_matrices,
  pi_pairs/l2r_matrices/r2l_matrices, propagate) to get real forward+backward
  modal amplitudes per cell, then reconstructs the field itself through
  stack_meow_mode_fields() -- the SAME device-frame remap
  project_interface_onto_body() and solve_body_modes() already use and this
  repo has already validated. So every individual piece here is validated
  elsewhere; what's new is driving that cascade from a REAL captured
  excitation vector instead of a canonical single-mode one, and doing it
  across 5 non-uniform regions instead of 1.
- This gets the interface scattering the straight-body diagonal approximation
  cannot: propagate()'s forward+backward solve is a real bidirectional EME
  cascade through every cell boundary, ramp included.

CAVEAT, stated plainly. N_MODES/C_TRANS here are NOT the converged production
values -- validation_markdown/mmi_2x2_waist_body_validation.md's own
"FOLLOW-UP 2026-09-09" finding is that N_MODES=14 is rank-deficient at the
waist for roughly half the ramp's distinct cross-sections, and C_TRANS is not
converged through 32. Both are exposed as env vars precisely so a picture at
whatever N_MODES turns out clean can be produced without editing this file.
This script is a DIAGNOSTIC of the stitching machinery and the qualitative
field picture -- not a source of any number to quote for the design.

Run:
    MMI2X2_STITCH_RES=20 python .../stitch_mmi2x2_waist_hybrid_fields.py
    python .../stitch_mmi2x2_waist_hybrid_fields.py          # real config
"""
import os

import builtins

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, EMESolver, SOI_220
from beef._utils import extrude_polygons
from beef.mode_identity import mesh_cell_centers, stack_meow_mode_fields
from beef.plot.hybrid_stitch import (dominant_component, load_captured_section,
                                     project_interface_onto_body,
                                     stitch_hybrid_fields)


def print(*a, **k):                      # noqa: A001 - rank-guard, see the
    try:                                  # straight device's own script
        import meep as _mp
        if not _mp.am_master():
            return
    except Exception:
        pass
    builtins.print(*a, **k)


def _is_master():
    try:
        import meep as _mp
        return bool(_mp.am_master())
    except Exception:
        return True


def _env_f(n, d):
    return float(os.environ.get(n, d))


def _env_i(n, d):
    return int(os.environ.get(n, d))


_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
RES = _env_f("MMI2X2_STITCH_RES", 40.0)
tapered_dir = os.path.join(_ROOT, "example_output", "mmi_2x2_tapered")
cap_dir = os.path.join(tapered_dir, "captured_input_section_res40_wl1.5475")

# --- pinned to mmi_2x2_waist_hybrid_vs_fdtd.py's design point ---------------
WW, LW, LT, L_WAIST = 2.24, 6.0, 6.0, 21.5
WL = 1.5475
BODY_STUB = 0.25       # FDTD-side padding past the seam, same as production
INPUT_IDX = 1
# N_MODES default TBD by the waist rank-deficiency diagnostic (see the
# validation markdown's FOLLOW-UP 2026-09-09); override once known, e.g.
#   MMI2X2_STITCH_NMODES=13
N_MODES = _env_i("MMI2X2_STITCH_NMODES", 14)
C_TRANS = _env_i("MMI2X2_STITCH_CTRANS", 8)   # step 5a's EME-only floor,
                                              # NOT shown converged for the
                                              # hybrid -- see the same finding

out_dir = os.path.join(tapered_dir, f"hybrid_stitched_res{int(RES)}")
os.makedirs(out_dir, exist_ok=True)

mmi = MMI(n_inputs=2, n_outputs=2, mm_length=L_WAIST, mm_width=2.8,
          waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
          waveguide_offset=1.5, waveguide_height=0.22,
          body_waist_width=WW, body_waist_length=LW,
          body_transition_length=LT)

half_body, half_total = mmi.mm_length / 2, mmi.total_length / 2
s_full = (L_WAIST - 2 * LT - LW) / 2
stub_in = s_full - BODY_STUB
x_b0, x_b1 = -half_body + BODY_STUB, +half_body - BODY_STUB
assert stub_in > 0, f"BODY_STUB={BODY_STUB} exceeds the design's stub s={s_full}"

print(f"waisted 2x2: total {mmi.total_length:.2f} um, body {mmi.mm_length:.2f} um "
      f"(Ww={WW} Lw={LW} Lt={LT})")
print(f"  EME body x in [{x_b0:+.2f}, {x_b1:+.2f}]  "
      f"regions: stub {stub_in:.2f}(1) | ramp {LT:.1f}({C_TRANS}) | "
      f"waist {LW:.1f}(1) | ramp {LT:.1f}({C_TRANS}) | stub {stub_in:.2f}(1)")

# --- 1. Load the ALREADY-CAPTURED FDTD input section (no FDTD re-run) ------
cap_npz = os.path.join(cap_dir, "section_fields_band1.npz")
if not os.path.exists(cap_npz):
    raise SystemExit(
        f"missing {cap_npz} -- run capture_mmi2x2_waist_fields.py "
        "MMI2X2_CAP_WHAT=input first (already run once; see the validation "
        "markdown's captured_input_section entry).")
cap = load_captured_section(cap_npz)
print(f"[1/4] loaded captured input section from {cap_dir}")
print(f"      interface plane at device x = {float(cap['if_pos']):+.4f} "
      f"(body starts at {x_b0:+.4f})")
assert abs(float(cap["if_pos"]) - x_b0) < 0.05, \
    "captured interface plane is not at the body boundary"

# --- 2. Build the WHOLE body as one continuous multi-cell mesh -------------
solver = EMESolver(num_modes=[N_MODES] * 5, num_cells=[1, C_TRANS, 1, C_TRANS, 1],
                   lengths=[stub_in, LT, LW, LT, stub_in],
                   oxide_width=5.0, oxide_height=1.5, res=RES,
                   material_stack=SOI_220, mode_solver="meep")

print(f"[2/4] building body: {1 + 2*C_TRANS + 2} cells, solving {N_MODES} "
      f"modes each...", flush=True)
gds_polys = mmi.section_gds(x_b0, x_b1, solver.oxide_width,
                            stub_left=0.0, stub_right=0.0)
structs = extrude_polygons(gds_polys, solver._extrusion_rules(np.array([WL])))
mesh = solver._build_mesh()
Ls = solver._flat_cell_lengths()
assert abs(sum(Ls) - (x_b1 - x_b0)) < 1e-9, (sum(Ls), x_b1 - x_b0)

import meow as mw
cells = mw.create_cells(structures=structs, mesh=mesh, Ls=Ls)
modes_per_cell = solver._modes_per_cell()
env = mw.Environment(wl=float(WL), T=solver.temperature)
css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
modes = [solver.compute_modes(cs, nm) for cs, nm in zip(css, modes_per_cell)]
print(f"      {len(cells)} cells solved, {sum(modes_per_cell)} mode-solves total")

# --- 3. Real EME cascade: track modes, cascade S-matrices, get fwd/bwd ------
import sax
from meow.eme.interface import compute_interface_s_matrices, compute_interface_s_matrix
from meow.eme.propagation import (compute_propagation_s_matrices, l2r_matrices,
                                  pi_pairs, propagate, r2l_matrices, track_modes)

print("[3/4] cascading the body (real bidirectional EME, not a diagonal "
      "approximation)...", flush=True)
sax_backend = sax.into[sax.Backend]("default")
tracked = track_modes(modes)
propagations = compute_propagation_s_matrices(tracked, cells)
interface_kwargs = dict(enforce_reciprocity=False)
interfaces = compute_interface_s_matrices(tracked, **interface_kwargs)
identity = compute_interface_s_matrix(tracked[0], tracked[0], **interface_kwargs)
pairs = pi_pairs(propagations, interfaces, sax_backend)
l2rs = l2r_matrices(pairs, identity, sax_backend)
r2ls = r2l_matrices(pairs, sax_backend)

# Project the REAL captured interface field onto cell 0's (tracked) basis --
# same function, same math as the straight device's stitch script.
y, z = mesh_cell_centers(mesh)
dy, dz = float(y[1] - y[0]), float(z[1] - z[0])
fields0 = stack_meow_mode_fields(tracked[0])
solved0 = dict(fields=fields0, y=y, z=z, dy=dy, dz=dz)
proj = project_interface_onto_body(cap, solved0)
comp = proj["component"]

# CONVENTION CONVERSION -- load-bearing, and the reason an earlier version of
# this script produced per-cell brightness steps inside the ramp.
#
# meow normalizes every mode to unit POWER FLUX (`normalize(mode,
# inner_product)` at meow/fde/meep.py:169, with inner_product =
# 0.5*int(ExH*).z dA), so a stored mode field satisfies
# int|E|^2 dA = (something cell-dependent), NOT 1: a mode squeezed into the
# 2.24 um waist carries a larger |E| for the same power than the same mode in
# the 2.8 um face. meow's S-matrix amplitudes are therefore coefficients on
# the RAW stored fields (see meow/eme/propagation.py:plot_fields, which
# multiplies amplitude by the raw mode field with no renormalization).
#
# project_interface_onto_body() instead returns coefficients against
# L2-normalized modes (it divides by sqrt(int|E|^2 dA) internally). Feeding
# those straight into meow's cascade, and then rendering meow's output
# amplitudes against L2-normalized fields again, mixes the two conventions
# and leaves a spurious per-mode, per-cell factor of sqrt(int|E|^2 dA) --
# which jumps at every cell boundary because the cross-section changes there.
#
# So: convert to meow's convention once, here, and stay in it.
_P0 = fields0[comp]
_norm0 = np.sqrt((np.abs(_P0) ** 2).sum(axis=(1, 2)) * dy * dz)
c_l = proj["c"] / _norm0
c_r = np.zeros(len(tracked[-1]), dtype=complex)

w = proj["weights"]
print("      |c_k|^2 (relative): "
      + ", ".join(f"{v:.4f}" for v in w / max(w.sum(), 1e-300)))
assert proj["captured_power"] <= 1.05, (
    f"captured_power {proj['captured_power']:.4f} -- see the straight "
    "device's stitch script for what this bound means and why 1.05 not 1.0.")

forwards, backwards = propagate(l2rs, r2ls, c_l, c_r)

# --- 4. Reconstruct the DEVICE-FRAME field cell by cell, then concatenate --
print(f"[4/4] reconstructing |{comp}|^2 through {len(cells)} cells...",
      flush=True)
# Samples PER MICRON, not per cell: cells here differ in length by ~8x (a
# 6 um waist against a 0.75 um ramp slice at C_TRANS=8), so a fixed
# samples-per-cell count would undersample the waist's beating while
# oversampling the ramp, and read as a per-cell texture change.
SAMPLES_PER_UM = 60
# The vertical slice is pinned ONCE, from the first cell. Recomputing
# argmax per cell risks it landing on a different z index in some cell
# (e.g. if a weakly-guided mode dominates the sum there), which would step
# the rendered magnitude at that cell boundary for a purely bookkeeping
# reason. The mesh is shared across all cells, so one index is correct
# everywhere.
IZ = int(np.argmax(np.abs(fields0[comp]).sum(axis=(0, 1))))
print(f"      vertical slice pinned at z = {z[IZ]:+.4f} um (index {IZ})")

xs_chunks, I_chunks = [], []
x_cursor = 0.0
norm_log = []
recon = dict(y=y, z=z, dy=dy, dz=dz, IZ=IZ, WL=WL, x_b0=x_b0, component=comp)
for i, (cell, md, fwd, bwd) in enumerate(zip(cells, tracked, forwards, backwards)):
    fields_i = stack_meow_mode_fields(md)
    P = fields_i[comp]                                    # (n_modes, ny, nz)
    norm_log.append(float(np.sqrt((np.abs(P[0]) ** 2).sum() * dy * dz)))
    neff = np.array([complex(m.neff).real for m in md])
    Li = float(cell.z_max - cell.z_min)
    n_s = max(2, int(round(SAMPLES_PER_UM * Li)))
    xs_local = np.linspace(0.0, Li, n_s, endpoint=False)
    phase_fwd = np.exp(2j * np.pi * np.outer(xs_local, neff) / WL)
    phase_bwd = np.exp(-2j * np.pi * np.outer(xs_local, neff) / WL)
    # RAW meow-normalized fields, matching the amplitude convention above.
    E_local = (phase_fwd * fwd[None, :] + phase_bwd * bwd[None, :]) @ P[:, :, IZ]
    xs_chunks.append(x_cursor + xs_local)
    I_chunks.append(np.abs(E_local) ** 2)
    x_cursor += Li
    # Cache everything the RENDER needs, so a rendering question (sampling,
    # slice choice, normalization convention) can be re-answered in seconds
    # instead of re-paying the mode solves -- which at production res is
    # hours. Full (y, z) mode slabs, not just the rendered z-slice, so the
    # slice choice itself stays revisitable.
    recon[f"P_{i}"] = P
    recon[f"neff_{i}"] = neff
    recon[f"fwd_{i}"] = np.asarray(fwd)
    recon[f"bwd_{i}"] = np.asarray(bwd)
    recon[f"L_{i}"] = Li
recon["n_cells"] = len(cells)

print("      per-cell sqrt(int|E_0|^2 dA) (the factor the old code divided "
      "out, cell by cell): "
      + ", ".join(f"{v:.3f}" for v in norm_log))

xs = np.concatenate(xs_chunks)
I = np.concatenate(I_chunks, axis=0)
body = dict(x=x_b0 + xs, y=y, I=I, component=comp)

_tag = f"waisted_res{int(RES)}_ct{C_TRANS}_nm{N_MODES}_{WL:.4f}um"
if _is_master():
    recon_path = os.path.join(out_dir, f"recon_cache_{_tag}.npz")
    np.savez_compressed(recon_path, **recon)
    print(f"      reconstruction cache -> {recon_path}")

if _is_master():
    stitch_hybrid_fields(
        [(cap, "input (FDTD, captured)")], body, out_dir, tag=_tag,
        port_ys=mmi.output_ys, body_half_width=mmi.mm_width / 2,
        device_label="waisted 2×2 MMI", geometry=mmi)

# Seam continuity, measured rather than eyeballed -- same convention as the
# straight device's stitch script.
comp_cap, _ = dominant_component(cap)
I_fdtd_end = (np.abs(np.asarray(cap[comp_cap])) ** 2)[-1]
yf = np.asarray(cap["y"])
I_f = np.interp(y, yf, I_fdtd_end)
sf, sb = I_f.sum(), body["I"][0].sum()
print(f"\nSEAM CONTINUITY at x = {x_b0:+.3f} um")
print(f"  integrated |E|^2: FDTD side {sf:.4e}  body side {sb:.4e}  "
      f"ratio {sb / sf if sf else float('nan'):.4f}")
print(f"  basis accounted for {proj['captured_power']:.4f} of the interface "
      f"field's power ({N_MODES} modes)")
print("\nNOTE: the seam ratio is NOT expected to be 1.0 -- the body side is "
      "the truncated-basis reconstruction, missing exactly the power "
      "reported as uncaptured. Reflected power at the OUTPUT end (this "
      "script only excites from the left) is not represented in the input "
      "capture either.")
print("\nCAVEAT: N_MODES/C_TRANS above are diagnostic values, not the "
      "converged production ones -- see this file's module docstring.")
