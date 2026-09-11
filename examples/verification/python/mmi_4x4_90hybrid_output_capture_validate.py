"""4x4 MMI 90-degree hybrid -- validates the reversed-excitation mechanism
needed for a real OUTPUT-taper field picture, at minimal cost, BEFORE
committing to the full per-body-mode capture across 3 lengths.

validation_markdown/mmi_4x4_90_hybrid_validation.md Step 8 field-plot work
(user request: stitched hybrid + full-device FDTD field plots at L=5/10/15
for the top-port excitation). The output taper is driven by a superposition
of body modes, not one physical mode, so a real output-taper field picture
needs FDTDSolver.solve_section() to launch each of the N_MODES body
eigenmodes individually FROM THE INTERIOR (body) side into the output taper
-- the opposite direction from every existing call in this device's own
scripts, which always excite from the PHYSICAL PORT side (fewer modes,
solve_section()'s own default "excite whichever side has fewer modes" rule).

Forcing this reversed direction needs solve_section()'s excite_port/
measure_port PortSpec override (both given -> bypasses the "fewer side"
tie-break entirely, per that method's own comment on override/reverse) --
exercised ONLY ONCE before in this repo (verify_ring_bend_reciprocity_3d.py,
a completely different port topology: a bent ring coupler's x/y-normal
ports, not two x-normal ports on a straight MMI slab). This script checks
that construction is right on THIS device's cross-section before trusting
any field it captures, using the one free, exact ground truth already
available: RECIPROCITY says this reversed call's own transfer coefficients
must equal M_in[TOP_PORT].transpose() -- already known from the production
hybrid script's own mirror-reuse trick, at zero new compute (a save_dir
cache hit, since M_in solves without capture_fields_dir are unaffected by
this script's own capture_fields_dir addition to the OUTPUT-side call).

Deliberately tests only N_TEST_MODES (default 3) of the full N_MODES=16,
not the whole per-body-mode set -- this is a correctness smoke test, not
the production capture (3 reversed launches vs 16, cutting validation cost
roughly 5x before scaling up).

Run:
    python examples/verification/python/mmi_4x4_90hybrid_output_capture_validate.py
"""
import os

import numpy as np

from beef import MMI, FDTDSolver, PortSpec, SOI_220

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification",
                           "mmi_4x4_90hybrid_output_capture_validate")
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_ROOT, "example_sim_cache"))

# Pinned to mmi_4x4_90hybrid_hybrid.py's own settled geometry exactly.
MM_WIDTH = 8.5
WAVEGUIDE_WIDTH = 0.35
TAPER_LENGTH = 4.0
TAPER_WIDTH = 1.30
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
N_MODES = 16
BODY_STUB = 0.25

# Same MMI4X4_MM_LENGTH override as the production scripts -- default here
# is the L=5.00 engine-agreement probe length (jobs 13024684/13024689),
# chosen so M_in's own cache from that already-completed run is reused for
# free rather than re-solved.
MM_LENGTH = float(os.environ.get("MMI4X4_MM_LENGTH", 5.00))
TOP_PORT_IDX = 3   # paper port 4, input_ys[-1] -- the "top" input this
                   # whole field-plot request is about.
N_TEST_MODES = int(os.environ.get("MMI4X4_N_TEST_MODES", 3))

WL_MIN, WL_MAX, N_WL = 1.29, 1.33, 3
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH, waveguide_width=WAVEGUIDE_WIDTH,
    taper_length=TAPER_LENGTH, taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET, waveguide_height=0.22,
    use_neighbor_aperture=True,
)
assert mmi.input_ys == mmi.output_ys

_fdtd_common = dict(
    res=40.0, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SOI_220, decay_tol=1e-4, decay_method="energy",
    energy_decay_dt=20.0,
)
input_solver = FDTDSolver(**_fdtd_common)
output_solver = FDTDSolver(**_fdtd_common)

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
boundaries = [-half_total, -half_body + BODY_STUB, +half_body - BODY_STUB, +half_total]

print(f"[validate] mm_length={MM_LENGTH}, boundaries={boundaries}, "
      f"TOP_PORT_IDX={TOP_PORT_IDX}, N_TEST_MODES={N_TEST_MODES}/{N_MODES}")

# --- Ground truth: M_in[TOP_PORT_IDX], the SAME call the production hybrid
# script makes for k=TOP_PORT_IDX -- should be a save_dir cache hit if this
# MM_LENGTH was already run through mmi_4x4_90hybrid_hybrid.py.
print("[validate] solving M_in[top] (should be a cache hit if L already run)...")
M_in_top = input_solver.solve_section(
    mmi, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=input_solver.waveguide_length, stub_right=0.0,
    input_idx=TOP_PORT_IDX, save_dir=save_dir,
)   # [n_wl, N_MODES, 1]
M_out_top_predicted = M_in_top.transpose(0, 2, 1)   # [n_wl, 1, N_MODES], by reciprocity/mirror symmetry
print(f"[validate] M_in_top shape {M_in_top.shape}, predicted M_out_top shape {M_out_top_predicted.shape}")

# --- The reversed-excitation test call: excite N_TEST_MODES body eigenmodes
# individually from the INTERIOR (boundaries[2], the LOW edge of this
# section) toward the physical port (boundaries[3], the HIGH edge).
# excite_port sign=+1 -> LOW edge (body side); measure_port sign=+1 ->
# HIGH edge (port side) -- see PortSpec's own docstring for this convention.
capture_dir = os.path.join(output_dir, f"L{MM_LENGTH:.2f}_fields")
print(f"[validate] solving reversed excitation ({N_TEST_MODES} body-mode launches)...")
M_out_top_test = output_solver.solve_section(
    mmi, wavelengths, boundaries[2], boundaries[3],
    n_modes_in=N_TEST_MODES, n_modes_out=1,
    stub_left=0.0, stub_right=output_solver.waveguide_length,
    excite_port=PortSpec(axis="x", sign=1, center=0.0, width=_fdtd_common["sim_width"]),
    measure_port=PortSpec(axis="x", sign=1, center=mmi.output_ys[TOP_PORT_IDX], width=WAVEGUIDE_OFFSET),
    capture_fields_dir=capture_dir,
    save_dir=save_dir,
)   # [n_wl, 1, N_TEST_MODES]
print(f"[validate] M_out_top_test shape {M_out_top_test.shape}")

predicted_slice = M_out_top_predicted[:, :, :N_TEST_MODES]
abs_diff = np.abs(M_out_top_test - predicted_slice)
rel_diff = abs_diff / np.maximum(np.abs(predicted_slice), 1e-12)

print("\n[validate] per (wavelength, mode) comparison -- test vs reciprocity-predicted:")
for iwl, wl in enumerate(wavelengths):
    for m in range(N_TEST_MODES):
        t = M_out_top_test[iwl, 0, m]
        p = predicted_slice[iwl, 0, m]
        print(f"  wl={wl:.4f} mode={m}: test={t:.6f} predicted={p:.6f} "
              f"abs_diff={abs_diff[iwl,0,m]:.2e} rel_diff={rel_diff[iwl,0,m]*100:.2f}%")

max_abs = float(abs_diff.max())
max_rel = float(rel_diff.max())
print(f"\n[validate] max abs diff = {max_abs:.4e}, max rel diff = {max_rel*100:.2f}%")
PASS = max_abs < 1e-2 and max_rel < 0.10
print(f"[validate] HYPOTHESIS {'CONFIRMED (mechanism is correct)' if PASS else 'NOT CONFIRMED -- do not trust captured fields'}")

np.savez(os.path.join(output_dir, f"validate_L{MM_LENGTH:.2f}.npz"),
         wavelengths=wavelengths, M_in_top=M_in_top,
         M_out_top_predicted=M_out_top_predicted, M_out_top_test=M_out_top_test,
         max_abs=max_abs, max_rel=max_rel, passed=PASS)
print(f"\n[validate] results + captured fields -> {output_dir}/")
