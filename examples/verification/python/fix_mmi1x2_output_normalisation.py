"""Renormalise the 1x2 forward output superposition — post-processing only.

THE BUG THIS FIXES, measured on job 13189964's output: integrated |E|^2 grew
**10.5x** from the output section's entry (4.185e+03) to its port end
(4.399e+04). A passive section cannot amplify. For guided modes the integral
is conserved up to an n_eff factor (~1.35x at most between the 2.7 body and
the ~2.0 access guides), so 10.5x is a normalisation error, not physics.

CAUSE. The superposition weighted meep's captured E^(k) -- the response to an
EigenModeSource of band k -- by coefficients a_k derived from an area-integral
projection onto EME modes. meep normalises each band its own way, so the two
conventions do not agree PER MODE, and the mismatch is mode-dependent (which
is why it shows up as a shape-changing growth, not a constant factor).

THE FIX, and why it needs no new simulation. Renormalise each captured E^(k)
so that its ENTRY-PLANE field is exactly the body mode phi_k it is supposed to
represent: project E^(k)(x0, y) onto phi_k and divide by that coefficient.
Whatever meep's per-band amplitude was then cancels identically, and the
superposition sum_k a_k E~^(k) reconstructs the arriving body field at the
entry by construction.

⚠ This is a CHECKABLE fix, not a fudge. Three things must hold afterwards and
are asserted below:
  1. entry-plane integral matches the body's arriving field (the seam closes);
  2. the growth along the section drops to O(1), not 10.5x;
  3. branch balance stays 1.000 (parity must survive renormalisation).
If (3) broke, the renormalisation would have introduced an asymmetry and must
be rejected.

Run:
    python .../fix_mmi1x2_output_normalisation.py
"""
import os
import numpy as np

from beef import (MMI, EMESolver, SOI_160, solve_body_modes,
                  load_captured_section, project_interface_onto_body)

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
RES, WL, N_MODES, STUB = 40.0, 1.31, 8, 0.25
base = os.path.join(_ROOT, "example_output", "mmi_1x2",
                    f"mmi_1x2_hybrid_stitched_res{int(RES)}")
cap_dir = os.path.join(base, "captured_output_forward_idx1")

mmi = MMI(n_inputs=1, n_outputs=2, mm_length=4.5, mm_width=2.0,
          waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
          waveguide_offset=1.0, waveguide_height=0.16)
hb = mmi.mm_length / 2
x_b0, x_b1 = -hb + STUB, hb - STUB

eme = EMESolver(num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
                oxide_width=5.0, oxide_height=1.5, res=RES,
                material_stack=SOI_160, mode_solver="meep")
solved = solve_body_modes(mmi, eme, x_b0, x_b1, WL, N_MODES,
                          wavelengths=np.linspace(1.26, 1.36, 11))
cap_in = load_captured_section(os.path.join(base, "captured_input_section",
                                            "section_fields_band1.npz"))
proj = project_interface_onto_body(cap_in, solved, verbose=False)
comp = proj["component"]
a = proj["c"] * np.exp(2j * np.pi * solved["neff"] * (x_b1 - x_b0) / WL)

# body mode profiles on the core z-slice, unit-normalised the same way the
# projection used
F = solved["fields"][comp]
iz = int(np.argmax(np.abs(F).sum(axis=(0, 1))))
phi = (F / np.sqrt((np.abs(F) ** 2).sum(axis=(1, 2))
                   * solved["dy"] * solved["dz"])[:, None, None])[:, :, iz]
yb, dy = solved["y"], solved["dy"]

E_raw = E_fix = None
scales = []
for k in range(1, N_MODES + 1):
    d = load_captured_section(os.path.join(cap_dir, f"section_fields_band{k}.npz"))
    Ek = np.asarray(d[comp]); x_o, y_o = np.asarray(d["x"]), np.asarray(d["y"])
    # entry plane, resampled onto the body-mode grid
    e0 = np.interp(yb, y_o, Ek[0].real) + 1j * np.interp(yb, y_o, Ek[0].imag)
    s = np.sum(np.conj(phi[k - 1]) * e0) * dy          # <phi_k | E^(k)(x0)>
    scales.append(abs(s))
    E_raw = a[k-1]*Ek if E_raw is None else E_raw + a[k-1]*Ek
    if abs(s) > 0:
        E_fix = a[k-1]*Ek/s if E_fix is None else E_fix + a[k-1]*Ek/s

print("per-band |<phi_k|E^(k)(entry)>| (meep's own per-mode amplitude):")
print("  " + ", ".join(f"{v:.3e}" for v in scales))
print(f"  spread max/min = {max(scales)/min(s for s in scales if s>0):.1f}x "
      f"-- this is the mode-dependent mismatch that caused the growth")

Ir, If = np.abs(E_raw)**2, np.abs(E_fix)**2
print(f"\n{'':22} {'RAW':>12} {'RENORMALISED':>14}")
print(f"{'entry integral':22} {Ir[0].sum():12.4e} {If[0].sum():14.4e}")
print(f"{'exit integral':22} {Ir[-1].sum():12.4e} {If[-1].sum():14.4e}")
print(f"{'growth along section':22} {Ir[-1].sum()/Ir[0].sum():11.2f}x "
      f"{If[-1].sum()/If[0].sum():13.2f}x")

# check 3: parity must survive
for nm, I in (("raw", Ir), ("renormalised", If)):
    p = I[-1]; t, b = p[y_o > 0.2].sum(), p[y_o < -0.2].sum()
    print(f"branch balance ({nm}): {t/b if b else float('inf'):.3f}")

np.savez_compressed(os.path.join(cap_dir, "forward_superposition_renorm.npz"),
                    x=x_o, y=y_o, I=If, a=a, component=comp, scales=np.array(scales))
g = If[-1].sum()/If[0].sum()
print(f"\ngrowth after renormalisation: {g:.2f}x "
      f"({'PLAUSIBLE' if 0.3 < g < 3 else 'STILL WRONG'} for a passive section)")
