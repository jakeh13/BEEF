"""Plumbing + physics check for solve_section(capture_fields_dir=...).

Step 2 of validation_markdown/hybrid_field_extraction_validation.md. Cheap by
design: 2-D, res 15, ONE wavelength, loose decay_tol -- seconds. It checks the
MECHANISM and that the captured field is a real guided mode, not the production
physics.

Two captures are verified, because Step 3 needs both:
  * the LONGITUDINAL x-y slice at the core z-plane  -> for plotting
  * the TRANSVERSE plane at the measure boundary    -> for modal projection

The transverse one is not a nicety. A longitudinal slice's interface column is
a 1-D line in y, and projecting a 1-D cut onto a mode basis is the pathology
beef.plot.hybrid_body.body_mode_weights documents from measurement: captured
power 1.5857 (impossible for an orthogonal basis) and |c|^2 = 0.5876 on a
parity-forbidden mode.

⚠ RUN THIS UNDER MPI, NOT JUST SERIALLY. `get_dft_array()` is a COLLECTIVE
call, so a rank-ordering bug in the capture path is INVISIBLE to a 1-rank run
(`mp.am_master()` is always true there). Exactly that shipped once: the
interface extraction sat after the master-only guard, and job 12736213
deadlocked 24 ranks for 40 minutes -- rank 0 waiting in a collective the other
23 had already returned from -- while the serial check passed cleanly.

Run BOTH:
    python .../check_solve_section_field_capture.py
    ~/.conda/envs/jmh/bin/mpirun -np 4 python .../check_solve_section_field_capture.py

(4 ranks on the login node takes ~1.4 s and is enough to expose a desync.
Note this cluster's mpirun is MPICH-flavoured: no --oversubscribe flag.)
"""
import os
import numpy as np

from beef import MMI, FDTDSolver, SOI_160

_HERE = os.path.dirname(os.path.abspath(__file__))
out = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "verification",
    "solve_section_field_capture"))
os.makedirs(out, exist_ok=True)

mmi = MMI(n_inputs=1, n_outputs=2, mm_length=4.5, mm_width=2.0,
          waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
          waveguide_offset=1.0, waveguide_height=0.16)
solver = FDTDSolver(res=15.0, sim_width=5.0, sim_height=1.5, dpml=0.5,
                    source_margin=0.0, monitor_margin=0.25, three_d=False,
                    material_stack=SOI_160, decay_tol=1e-6)

x0, x1 = -mmi.total_length / 2, -mmi.mm_length / 2 + 0.25
T = solver.solve_section(mmi, np.array([1.31]), x0, x1,
                         n_modes_in=1, n_modes_out=2,
                         capture_fields_dir=out)
print(f"transfer matrix shape {np.asarray(T).shape} (unchanged by capture)")

f = os.path.join(out, "section_fields_band1.npz")
d = np.load(f)
print(f"\nnpz keys: {sorted(d.files)}")

# --- longitudinal slice ----------------------------------------------------
print(f"\nLONGITUDINAL slice: Ez {d['Ez'].shape} (nx, ny)")
print(f"  x [{d['x'].min():.3f}, {d['x'].max():.3f}] um  "
      f"(section is [{x0:.3f}, {x1:.3f}] -- device frame, stitchable)")
assert abs(d['x'].min() - x0) < 1e-6 and abs(d['x'].max() - x1) < 1e-6, \
    "longitudinal x coords are NOT the section's own device-frame bounds"
assert str(d['frame']) == "device"

# Which component actually carries the mode? In 2-D the polarizations decouple
# and Ez dominates; in 3-D a TE-like mode's in-plane Ey does. Never assume.
mags = {c: float(np.abs(d[c]).max()) for c in ("Ex", "Ey", "Ez")}
dom = max(mags, key=mags.get)
print(f"  component magnitudes {({k: f'{v:.2e}' for k, v in mags.items()})}"
      f" -> dominant {dom}")
assert mags[dom] > 0, "captured field is identically zero"

E = d[dom]
prof = (np.abs(E) ** 2).sum(axis=0)
prof = prof / prof.max()
y = d["y"]
above = y[prof > 0.5]
fwhm = above.max() - above.min()
conf = prof[np.abs(y) < 0.5].sum() / prof.sum()
print(f"  transverse profile: peak y={y[int(np.argmax(prof))]:+.3f}, "
      f"FWHM {fwhm:.3f} um, power in |y|<0.5 = {conf:.4f}")
assert conf > 0.9, f"field not confined to the waveguide (conf={conf:.3f})"
assert 0.15 < fwhm < 1.0, f"FWHM {fwhm:.3f} um implausible for a 0.35 um guide"

# --- transverse interface plane -------------------------------------------
assert "if_Ez" in d.files, "transverse interface capture MISSING -- Step 3 needs it"
print(f"\nTRANSVERSE interface plane at {str(d['if_axis'])}"
      f"={float(d['if_pos']):+.3f} um")
print(f"  if_Ez shape {np.atleast_1d(d['if_Ez']).shape}, "
      f"in-plane axes ({str(d['if_u_axis'])}, {str(d['if_v_axis'])})")
print(f"  if_u [{d['if_u'].min():+.3f}, {d['if_u'].max():+.3f}]")
assert float(np.abs(d["if_" + dom]).max()) > 0, "interface field is zero"
# The plane must sit on a SECTION BOUNDARY in the DEVICE frame. mea_along is a
# shifted-sim coordinate inside solve_section, so this catches a missing frame
# conversion -- it originally read +0.625, outside the section entirely.
ifp = float(d["if_pos"])
assert min(abs(ifp - x0), abs(ifp - x1)) < 0.05, (
    f"interface plane at {ifp:+.3f} is not on a section boundary "
    f"({x0:+.3f} / {x1:+.3f}) -- shifted-sim frame leaked into the npz")
print(f"  device-frame check: {ifp:+.3f} um lands on the boundary "
      f"x1={x1:+.3f} (was +0.625 in sim coords before the fix)")

print("\nOK -- both captures present, coordinates device-frame, "
      "field is a confined guided mode")
