"""1x2 MMI — FORWARD-excited output section, completing the whole-device map.

The last 19% of validation_markdown/hybrid_field_extraction_validation.md's
Step 3. The input section and body were stitched at res 40 (job 12736849);
the output section was excluded because `solve_section` excites from whichever
boundary has FEWER modes, which for 8-in/1-out is the PORT. Job 12738781
measured what that costs: a branch power ratio of 1899 (one branch lit, the
other dark) where forward 50/50 operation needs ~1.

NO LIBRARY CHANGE IS NEEDED -- the escape hatch already exists. fdtd.py:3232:

    override = excite_port is not None or measure_port is not None
    reverse = False if override else layout.reverse

so passing ANY `excite_port` PortSpec forces reverse=False. Its own comment
says why it was added: the `n_out <= n_in` tie-break "silently picks the wrong
physical side whenever n_in == n_out (verified: the standalone-bend case)".
See also [[solve-section-excites-fewer-mode-side]], where a second session
measured the same failure on the waisted 2x2 (both input arms lit, ratio
0.986 vs a forward 0.001) and spent two 17-minute runs finding it.

For an EXCITE port, sign=+1 binds the plane to the LOW edge of the section,
i.e. x0 = the body side. That is what we want: launch each body mode forward
into the output section.

COST: reverse=False means n_excite = n_modes_in = 8, so EIGHT FDTD runs
instead of one (~4 min each at res 40 on 8 ranks from job 12738781's 3:57).
Step 2's capture writes one npz per excited band, so all eight land
automatically.

THE SUPERPOSITION. The forward output field is
    E_out(x,y) = sum_k a_k E^(k)(x,y)
with a_k the body's ARRIVING amplitudes at x = x_b1: the interface
coefficients c_k propagated the body's length,
    a_k = c_k * exp(2j*pi*n_eff_k*L_body/wl).

⚠ NORMALISATION IS THE RISK HERE, and it is checked rather than assumed.
c_k come from an area-integral projection onto EME modes; E^(k) is meep's
response to a unit EigenModeSource of band k. Those two conventions need not
agree per-mode. Two independent checks are applied below:
  1. the two output branches must come out ~50/50 (this device is a splitter);
  2. the per-branch transmission must land near the hybrid's own validated
     0.3749 at 1.31 um (Step 0, job 12674181).
If either fails, the weighting is wrong and the map must not be published.

Run (batch, 8 ranks):
    srun --mpi=pmi2 python3 .../capture_mmi1x2_output_forward.py
"""
import os

import numpy as np

from beef import (MMI, EMESolver, FDTDSolver, SOI_160, solve_body_modes,
                  load_captured_section, project_interface_onto_body)
from beef.solvers.fdtd import PortSpec, section_fdtd_layout

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
RES = float(os.environ.get("MMI1X2_STITCH_RES", 40.0))
OUT_IDX = int(os.environ.get("MMI1X2_OUTPUT_IDX", 1))
base = os.path.join(_ROOT, "example_output", "mmi_1x2",
                    f"mmi_1x2_hybrid_stitched_res{int(RES)}")
cap_dir = os.path.join(base, f"captured_output_forward_idx{OUT_IDX}")
os.makedirs(cap_dir, exist_ok=True)

WL, N_MODES, STUB = 1.31, 8, 0.25
mmi = MMI(n_inputs=1, n_outputs=2, mm_length=4.5, mm_width=2.0,
          waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
          waveguide_offset=1.0, waveguide_height=0.16)
hb, ht = mmi.mm_length / 2, mmi.total_length / 2
x_b0, x_b1 = -hb + STUB, hb - STUB          # body
x_o0, x_o1 = hb - STUB, ht                  # output section

SIM_W = 5.0
fdtd = FDTDSolver(res=RES, sim_width=SIM_W, sim_height=1.5, dpml=0.5,
                  source_margin=0.0, monitor_margin=0.25, three_d=True,
                  material_stack=SOI_160, decay_tol=1e-9,
                  minimum_run_time=100.0, half_domain_aperture=True)

# Confirm the override flips the excite side BEFORE spending 8 runs --
# section_fdtd_layout() is pure geometry math and costs nothing.
lay = section_fdtd_layout(mmi, x_o0, x_o1, 0.0, 0.0, N_MODES, 1, SIM_W,
                          0.0, 0.25, True, {"output_idx": OUT_IDX})
print(f"default layout: reverse={lay.reverse}, n_excite={lay.n_excite} "
      f"(excites from the PORT -- what we are overriding)")

# sign=+1 -> LOW edge of the section = x0 = the body side. Full sim_width
# aperture, matching the span the body-side interface monitor used
# (if_u spanned [-2.50, +2.50] in the input-section capture).
excite = PortSpec(axis="x", sign=1, center=0.0, width=SIM_W)
print(f"excite_port PortSpec(axis=x, sign=+1, center=0.0, width={SIM_W}) "
      f"-> forces reverse=False, n_excite={N_MODES}")

print(f"\nrunning {N_MODES} forward excites on the output section "
      f"x in [{x_o0:+.2f}, {x_o1:+.2f}]...")
T_out = fdtd.solve_section(mmi, np.array([WL]), x_o0, x_o1,
                           n_modes_in=N_MODES, n_modes_out=1,
                           excite_port=excite,
                           capture_fields_dir=cap_dir, output_idx=OUT_IDX)
print(f"transfer matrix {np.asarray(T_out).shape}")
fs = sorted(f for f in os.listdir(cap_dir) if f.endswith(".npz"))
print(f"npz written: {len(fs)} (expect {N_MODES} -- one per forward excite)")
assert len(fs) == N_MODES, (
    f"expected {N_MODES} captures, got {len(fs)} -- the excite_port override "
    "did not force forward excitation")

# --- the body's arriving amplitudes ---------------------------------------
cap_in = load_captured_section(os.path.join(base, "captured_input_section",
                                            "section_fields_band1.npz"))
eme = EMESolver(num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
                oxide_width=SIM_W, oxide_height=1.5, res=RES,
                material_stack=SOI_160, mode_solver="meep")
solved = solve_body_modes(mmi, eme, x_b0, x_b1, WL, N_MODES,
                          wavelengths=np.linspace(1.26, 1.36, 11))
proj = project_interface_onto_body(cap_in, solved, verbose=False)
L_body = x_b1 - x_b0
a = proj["c"] * np.exp(2j * np.pi * solved["neff"] * L_body / WL)
print(f"\nbody arriving amplitudes |a_k|^2 (relative): "
      + ", ".join(f"{v:.4f}" for v in (np.abs(a)**2 / (np.abs(a)**2).sum())))

# --- superpose the eight captured forward fields ---------------------------
comp = proj["component"]
E = None
for k in range(1, N_MODES + 1):
    d = load_captured_section(os.path.join(cap_dir,
                                           f"section_fields_band{k}.npz"))
    Ek = np.asarray(d[comp])
    E = a[k - 1] * Ek if E is None else E + a[k - 1] * Ek
    if k == 1:
        x_o, y_o = np.asarray(d["x"]), np.asarray(d["y"])
I_out = np.abs(E) ** 2
np.savez_compressed(os.path.join(cap_dir, "forward_superposition.npz"),
                    x=x_o, y=y_o, I=I_out, a=a, component=comp)

# --- CHECK 1: the two branches must be ~50/50 ------------------------------
prof = I_out[-1]
top = prof[y_o > 0.2].sum(); bot = prof[y_o < -0.2].sum()
ratio = top / bot if bot > 0 else float("inf")
print(f"\nCHECK 1 -- branch balance at the port end (x={x_o.max():+.2f})")
print(f"  top {top:.4e}  bottom {bot:.4e}  ratio {ratio:.3f}")
print(f"  reverse-excited run gave 1899; a correct forward map needs ~1.0")
print(f"  -> {'PASS' if 0.5 < ratio < 2.0 else 'FAIL'}")
print("\n(CHECK 2, per-branch transmission vs the hybrid's validated 0.3749, "
      "needs a power\n normalisation this script does not yet carry -- "
      "reported as a ratio only, see docstring.)")
