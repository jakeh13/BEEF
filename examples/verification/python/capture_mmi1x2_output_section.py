"""Capture the 1x2 OUTPUT section's field as solve_section actually runs it.

Purpose: settle empirically whether the reverse-excited capture can serve as
the output panel of a whole-device field map.

`build_run_config` sets `reverse = n_out <= n_in` (fdtd.py:735) and excites
from the boundary with FEWER modes. The output section is 8 modes on the body
side, 1 on the port side, so it is excited FROM THE PORT -- backwards relative
to device operation.

For the transfer MATRIX that is exact and is an 8x saving. For the FIELD, by
reciprocity the reverse field's conjugate is a valid FORWARD field, but for the
input superposition that couples optimally into port 1 (coefficients ~ S_1k*),
NOT for the superposition the body actually delivers (measured 0.69/0.28/0.035
at the input/body seam).

FALSIFIABLE PREDICTION this script tests: driving one output port backwards
should illuminate essentially ONE branch, whereas forward operation of a 50/50
splitter fills BOTH. If the capture shows one branch lit, the reverse map
cannot be used as the forward output panel, and the whole-device map needs the
section re-excited from the body side (8 excites).

One excite, so ~1 FDTD run. Config pinned to mmi_1x2_hybrid.py, including
half_domain_aperture=True for the output section.
"""
import os
import numpy as np

from beef import MMI, FDTDSolver, SOI_160, load_captured_section

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
RES = float(os.environ.get("MMI1X2_STITCH_RES", 40.0))
OUT_IDX = int(os.environ.get("MMI1X2_OUTPUT_IDX", 1))    # 1 = top branch
base = os.path.join(_ROOT, "example_output", "mmi_1x2",
                    f"mmi_1x2_hybrid_stitched_res{int(RES)}")
cap_dir = os.path.join(base, f"captured_output_section_idx{OUT_IDX}")
os.makedirs(cap_dir, exist_ok=True)

WL, N_MODES, STUB = 1.31, 8, 0.25
mmi = MMI(n_inputs=1, n_outputs=2, mm_length=4.5, mm_width=2.0,
          waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
          waveguide_offset=1.0, waveguide_height=0.16)
hb, ht = mmi.mm_length / 2, mmi.total_length / 2
x0, x1 = hb - STUB, ht

fdtd = FDTDSolver(res=RES, sim_width=5.0, sim_height=1.5, dpml=0.5,
                  source_margin=0.0, monitor_margin=0.25, three_d=True,
                  material_stack=SOI_160, decay_tol=1e-9,
                  minimum_run_time=100.0, half_domain_aperture=True)
print(f"output section x in [{x0:+.2f}, {x1:+.2f}], output_idx={OUT_IDX}, "
      f"res {RES:.0f}")
print(f"n_modes_in=8 (body side), n_modes_out=1 (port) -> reverse=True, "
      f"excited FROM THE PORT")
T = fdtd.solve_section(mmi, np.array([WL]), x0, x1,
                       n_modes_in=N_MODES, n_modes_out=1,
                       capture_fields_dir=cap_dir, output_idx=OUT_IDX)
print(f"transfer matrix {np.asarray(T).shape}")

fs = sorted(f for f in os.listdir(cap_dir) if f.endswith(".npz"))
print(f"npz files written: {fs}  (one per EXCITED band; 1 confirms reverse)")
d = load_captured_section(os.path.join(cap_dir, fs[0]))
mags = {c: float(np.abs(d[c]).max()) for c in ("Ex", "Ey", "Ez")}
comp = max(mags, key=mags.get)
I = np.abs(d[comp]) ** 2
x, y = d["x"], d["y"]
print(f"dominant {comp}; slice {I.shape}, x [{x.min():+.2f},{x.max():+.2f}]")

# THE TEST: how is power distributed between the two branches at the port end?
prof = I[-1]
top = prof[y > 0.2].sum(); bot = prof[y < -0.2].sum()
print(f"\nPOWER SPLIT at the port end (x={x.max():+.2f}):")
print(f"  y > +0.2 (top branch)  {top:.4e}")
print(f"  y < -0.2 (bottom)      {bot:.4e}")
r = top / bot if bot > 0 else float('inf')
print(f"  ratio top/bot = {r:.2f}")
print(f"\nforward 50/50 operation would give ~1.0; a single-branch reverse "
      f"drive gives a large asymmetry")
np.savez(os.path.join(cap_dir, "branch_split.npz"), top=top, bot=bot, ratio=r)
