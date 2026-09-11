"""1x2 MMI — Step 3: stitch the FDTD input section's captured field to the EME
body, driven by the REAL interface field rather than an assumed launch.

validation_markdown/hybrid_field_extraction_validation.md Step 3.

WHAT THIS DOES THAT STEP 1 COULD NOT. plot_hybrid_body_fields() had to invent
the launch (a Gaussian standing in for whatever the FDTD section delivers),
because no FDTD field data existed anywhere in the library. Step 2's
capture_fields_dir now provides it, so here the body's modal coefficients are
PROJECTED FROM THE CAPTURED INTERFACE FIELD. Two consequences the plan called
out as the hard part come for free:

  * amplitude -- the coefficients carry the captured field's own absolute
    units, so body and FDTD panels share one scale;
  * propagation phase -- the coefficients are measured AT the shared boundary
    and the body propagator is referenced to that same plane, so the phase
    continues by construction.

⚠ THE OUTPUT SECTION IS DELIBERATELY NOT STITCHED, and the reason is physics,
not laziness. `build_run_config` sets `reverse = n_out <= n_in`
(`fdtd.py:735`), so with n_modes_in=8 (body side) and n_modes_out=1 (port
side) the output section is EXCITED FROM ITS PORT and measured at the body
side. Reciprocity keeps its transfer MATRIX valid, but its captured field map
is of light launched backwards from the output port -- not of the field the
body delivers forward into it. Stitching it here would produce a picture that
looks complete and describes a different experiment. Getting a forward output
map needs the section re-run excited from the body side (8 excites, so 8 FDTD
runs), which is real extra cost and out of scope for this step.

Config pinned to mmi_1x2/python/mmi_1x2_hybrid.py. MMI1X2_STITCH_RES overrides
the resolution for a cheap pipeline check (res 20 validates the mechanics in
minutes; res 40 is the production value).

Run:
    MMI1X2_STITCH_RES=20 python .../stitch_mmi1x2_hybrid_fields.py
"""
import os

import builtins

import numpy as np

from beef import (MMI, EMESolver, FDTDSolver, SOI_160, solve_body_modes,
                  load_captured_section, project_interface_onto_body,
                  stitch_hybrid_fields)

# RANK-GUARD every print. Under `srun --mpi=pmi2 -n 24` every rank executes
# this script, so an unguarded print emits 24 identical lines -- observed in
# job 12736213, whose "[1/4]" banner appears 24 times. Harmless to the physics
# but it makes the log unparseable and would break any downstream grep that
# counts stages. meep's own C++ prints are already rank-0 only, which is why
# `run 0 finished at t` appears once and this banner did not.
def print(*a, **k):                      # noqa: A001 - deliberate shadow
    try:
        import meep as _mp
        if not _mp.am_master():
            return
    except Exception:
        pass
    builtins.print(*a, **k)

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
RES = float(os.environ.get("MMI1X2_STITCH_RES", 40.0))
out_dir = os.path.join(_ROOT, "example_output", "mmi_1x2",
                       f"mmi_1x2_hybrid_stitched_res{int(RES)}")
cap_dir = os.path.join(out_dir, "captured_input_section")
os.makedirs(cap_dir, exist_ok=True)

# --- pinned to mmi_1x2_hybrid.py -------------------------------------------
WL_MIN, WL_MAX, N_WL = 1.26, 1.36, 11
WL_PLOT = 1.31
N_MODES = 8
BODY_STUB = 0.25

mmi = MMI(n_inputs=1, n_outputs=2, mm_length=4.5, mm_width=2.0,
          waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
          waveguide_offset=1.0, waveguide_height=0.16)

half_body, half_total = mmi.mm_length / 2, mmi.total_length / 2
x_in0, x_in1 = -half_total, -half_body + BODY_STUB      # FDTD input section
x_b0, x_b1 = -half_body + BODY_STUB, +half_body - BODY_STUB   # EME body
assert abs(x_in1 - x_b0) < 1e-12, "input section's right edge must be the body's left edge"

# --- 1. FDTD input section WITH capture ------------------------------------
# ONE FDTD run: n_modes_in=1 so `reverse` is False and there is a single
# excite. The capture writes both the longitudinal slice (for the map) and the
# transverse plane at x = x_in1 (for the projection).
fdtd = FDTDSolver(res=RES, sim_width=5.0, sim_height=1.5, dpml=0.5,
                  source_margin=0.0, monitor_margin=0.25, three_d=True,
                  material_stack=SOI_160, decay_tol=1e-9,
                  minimum_run_time=100.0)
print(f"[1/4] FDTD input section x in [{x_in0:+.2f}, {x_in1:+.2f}] "
      f"at res {RES:.0f}, capturing fields...")
T_in = fdtd.solve_section(mmi, np.array([WL_PLOT]), x_in0, x_in1,
                          n_modes_in=1, n_modes_out=N_MODES,
                          capture_fields_dir=cap_dir)
print(f"      transfer matrix {np.asarray(T_in).shape}")

cap = load_captured_section(os.path.join(cap_dir, "section_fields_band1.npz"))
print(f"      interface plane at device x = {float(cap['if_pos']):+.4f} "
      f"(body starts at {x_b0:+.4f})")
assert abs(float(cap["if_pos"]) - x_b0) < 0.05, \
    "captured interface plane is not at the body boundary"

# --- 2. Body modes ---------------------------------------------------------
eme = EMESolver(num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
                oxide_width=5.0, oxide_height=1.5, res=RES,
                material_stack=SOI_160, mode_solver="meep")
print(f"[2/4] solving {N_MODES} body modes at {WL_PLOT} um...")
solved = solve_body_modes(mmi, eme, x_b0, x_b1, WL_PLOT, N_MODES,
                          wavelengths=np.linspace(WL_MIN, WL_MAX, N_WL))
print("      n_eff: " + ", ".join(f"{n:.6f}" for n in solved["neff"]))

# --- 3. Project the REAL interface field onto the body basis ---------------
print("[3/4] projecting the captured interface field onto the body modes...")
proj = project_interface_onto_body(cap, solved)
w = proj["weights"]
print("      |c_k|^2 (relative): "
      + ", ".join(f"{v:.4f}" for v in w / max(w.sum(), 1e-300)))

# A centred single-mode input must drive only EVEN body modes. This is the
# same parity check that caught the 1-D-projection bug in Step 1, applied now
# to a REAL field rather than an analytic launch -- so it also tests that the
# capture and the mode mesh are aligned in y.
odd = w[1::2].sum() / max(w.sum(), 1e-300)
print(f"      odd-mode fraction {odd:.3e} (a centred input should give ~0)")

# --- 4. Propagate the body on those coefficients and stitch ---------------
print("[4/4] propagating the body and stitching...")
iz = int(np.argmax(np.abs(solved["fields"][proj["component"]]).sum(axis=(0, 1))))
phi = (solved["fields"][proj["component"]] /
       np.sqrt((np.abs(solved["fields"][proj["component"]]) ** 2).sum(axis=(1, 2))
               * solved["dy"] * solved["dz"])[:, None, None])[:, :, iz]
xs = np.linspace(0.0, x_b1 - x_b0, 400)
E_body = (np.exp(2j * np.pi * np.outer(xs, solved["neff"]) / WL_PLOT)
          * proj["c"][None, :]) @ phi
body = dict(x=x_b0 + xs, y=solved["y"], I=np.abs(E_body) ** 2,
            component=proj["component"])

res = stitch_hybrid_fields(
    [(cap, "input (FDTD, captured)")], body, out_dir,
    tag=f"mmi1x2_res{int(RES)}_{WL_PLOT:.3f}um",
    port_ys=mmi.output_ys, body_half_width=mmi.mm_width / 2,
    device_label="1×2 MMI")

# Continuity across the seam is the whole point, so measure it rather than
# eyeballing the picture: the body's first column against the FDTD slice's
# last, both now in absolute units on a shared scale.
comp, _ = __import__("beef").dominant_component(cap)
I_fdtd_end = (np.abs(np.asarray(cap[comp])) ** 2)[-1]
I_body_start = body["I"][0]
yf, yb = np.asarray(cap["y"]), solved["y"]
I_f = np.interp(yb, yf, I_fdtd_end)
sf, sb = I_f.sum(), I_body_start.sum()
print(f"\nSEAM CONTINUITY at x = {x_b0:+.3f} um")
print(f"  integrated |E|^2: FDTD side {sf:.4e}  body side {sb:.4e}  "
      f"ratio {sb / sf if sf else float('nan'):.4f}")
print(f"  basis accounted for {proj['captured_power']:.4f} of the interface "
      f"field's power ({N_MODES} modes)")
# --- 5. STEP 4 CROSS-CHECK, free: field capture vs the transfer matrix ----
# The input section's transfer matrix column T_in[0, :, 0] IS the modal
# amplitude at its output boundary for a unit input mode -- computed by meep's
# own get_eigenmode_coefficients() at the flux monitor. The projection above
# computes the SAME physical quantity by a completely independent route: DFT
# field capture, then an area-integral projection onto EME-solved modes.
#
# Different mode solver (meep's eigensolver vs the EME basis), different
# extraction (flux-plane decomposition vs a field integral). Agreement is a
# real validation of the capture against the number the cascade actually uses;
# disagreement means the capture, the frame, or the mesh alignment is wrong.
#
# Compared as NORMALISED modal power, which removes the arbitrary overall
# scaling between the two routes.
Tcol = np.asarray(T_in)[0, :, 0]
p_T = np.abs(Tcol) ** 2
p_T = p_T / p_T.sum()
p_c = proj["weights"] / proj["weights"].sum()
print("\nSTEP 4 CROSS-CHECK -- capture+projection vs transfer matrix")
print(f"  {'mode':>4} {'|T|^2 norm':>12} {'|c|^2 norm':>12} {'delta':>10}")
for k in range(len(p_T)):
    print(f"  {k:>4} {p_T[k]:12.6f} {p_c[k]:12.6f} {p_c[k]-p_T[k]:+10.2e}")
print(f"  max |delta| {np.abs(p_c - p_T).max():.3e}   "
      f"sum|delta| {np.abs(p_c - p_T).sum():.3e}")
print("  (two independent routes to the same modal amplitudes: meep's "
      "eigenmode\n   decomposition at the flux plane, vs a DFT field capture "
      "projected onto\n   the EME mode basis)")

print("\nNOTE: the ratio is NOT expected to be 1.0. The body side is the "
      "truncated-basis reconstruction, so it is missing exactly the power the "
      "line above reports as uncaptured, plus any non-guided content the FDTD "
      "slice still carries at that plane.")
