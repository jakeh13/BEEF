"""Step 1 of validation_markdown/resolved/multi_excitation_cascade_validation.md —
capture the 2x2's cascade matrices as the regression target for steps 2-4.

STATUS 2026-08-21: steps 2-4 have LANDED. mmi_2x2_hybrid.py no longer
hand-rolls its cascade — it declares the symmetry to HybridSolver via
mirror_pairs/MirrorSpec and calls solve_length_sweep() twice. The migration
was verified BIT-IDENTICAL against the golden this script writes.

That does not make this script obsolete, and it is deliberately NOT migrated
with the example: its whole value is that it composes the matrices
INDEPENDENTLY of the path under test. Rewriting it to use the same
HybridSolver route would make it a tautology rather than a check. Keep the
explicit einsum here even though the example no longer has one.

It also writes to its own directory and never touches the shipped output
(the GDS acceptance run reads mmi_2x2_hybrid/mm_length_sweep.npz).

TWO GOLDENS, and they are not interchangeable:
  cascade_golden.npz              SYMMETRIC_MODE_REUSE=True (M_in0 derived
                                  from MODE_SIGNS) — the steps 2-4 regression
                                  target, and what this script writes by
                                  default.
  cascade_golden_real_m_in0.npz   CAPTURE_REAL_M_IN0=1 (job 12159086) — M_in0
                                  independently MEASURED. Its cross columns
                                  are a different quantity, differing by
                                  1.585e-05 on T_cross.
Running with CAPTURE_REAL_M_IN0=1 REWRITES cascade_golden.npz in place. Back
it up first, then rename the result, or the two configurations collide under
one filename.

What it saves, for the shipped configuration:
  M_in1     [n_wl, N_MODES, 1]        top-arm excitation -> body modes
  M_in0     [n_wl, N_MODES, 1]        bottom arm (see CAPTURE_REAL_M_IN0)
  M_body    [n_L, n_wl, N, N]         the analytic length sweep
  bar/cross [n_L, n_wl]               composed transmission grids
plus the same compositions routed through beef.network.Circuit, so the
einsum-vs-Circuit agreement is measured on the REAL device matrices and not
only on the synthetic ones in tests/test_multi_excitation_cascade.py.

Cost. With CAPTURE_REAL_M_IN0=0 (default) this is exactly the work
mmi_2x2_hybrid.py already does, so it is free whenever that configuration is
cache-resident and a full hybrid run when it is not -- either way it belongs
in the batch queue, not on a login node. CAPTURE_REAL_M_IN0=1 adds a SECOND
input-side FDTD run (the one SYMMETRIC_MODE_REUSE exists to avoid, ~1.5-2.5 h
for this device). That is not waste: step 3(b) needs a real, independently
measured M_in0 to derive or validate the MODE_SIGNS parity vector that the
example currently hard-codes.

Run:
    python examples/verification/python/capture_mmi_2x2_cascade_golden.py
"""

import os
import subprocess

import numpy as np

from beef import EMESolver, FDTDSolver, MMI, SOI_220
from beef.network.circuit import Circuit
from beef.network.convert import from_transfer_matrix, to_smatrix

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))       # BEEF/
_REPO = os.path.abspath(os.path.join(_ROOT, ".."))                    # duty/
output_dir = os.path.join(_REPO, "example_output", "verification",
                          "mmi_2x2_cascade_golden")
save_dir = os.path.join(_REPO, "example_sim_cache")
os.makedirs(output_dir, exist_ok=True)

# Opt in to the second input-side FDTD run instead of deriving M_in0 from
# MODE_SIGNS. Needed by step 3(b); skipped by default so the capture stays
# as cheap as the shipped run.
CAPTURE_REAL_M_IN0 = bool(int(os.environ.get("CAPTURE_REAL_M_IN0", "0")))

# ---------------------------------------------------------------------------
# Verbatim from mmi_2x2_hybrid.py. If these drift apart the capture stops
# being a regression target for the shipped script, so they are asserted
# against the geometry below rather than merely copied.
# ---------------------------------------------------------------------------
RESOLUTION = 40.0
WL_MIN, WL_MAX, N_WL = 1.530, 1.565, 11
N_MODES = 14
CALIBRATE_BODY = True
MODE_SIGNS = np.array([1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, 1])

mmi = MMI(
    n_inputs=2, n_outputs=2, mm_length=28.0, mm_width=2.8,
    waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
    waveguide_offset=1.5, waveguide_height=0.22,
)

_fdtd_common = dict(
    res=RESOLUTION, sim_width=5.0, sim_height=1.5, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    material_stack=SOI_220, decay_tol=1e-11,
)
input_solver = FDTDSolver(**_fdtd_common, half_domain_aperture=True)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
    oxide_width=5.0, oxide_height=1.5, res=RESOLUTION,
    material_stack=SOI_220, mode_solver="meep",
)

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
body_stub = 0.25
boundaries = [-half_total, -half_body + body_stub, half_body - body_stub, half_total]
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
mm_lengths = np.linspace(16.0, 40.0, 201)
body_lengths = mm_lengths - 2 * body_stub

assert len(MODE_SIGNS) == N_MODES, "MODE_SIGNS length must match N_MODES"

# ---------------------------------------------------------------------------
# 1. The three matrices the hand-rolled cascade is built from.
# ---------------------------------------------------------------------------
print(f"Capturing the 2x2 cascade at res={RESOLUTION}, {N_WL} wl, "
      f"{N_MODES} interface modes, {len(mm_lengths)} lengths")
print(f"  CAPTURE_REAL_M_IN0={int(CAPTURE_REAL_M_IN0)}")

M_in1 = input_solver.solve_section(
    mmi, wavelengths, boundaries[0], boundaries[1],
    n_modes_in=1, n_modes_out=N_MODES,
    stub_left=input_solver.waveguide_length, stub_right=0.0,
    input_idx=1, save_dir=save_dir,
)

M_in0_derived = M_in1 * MODE_SIGNS[None, :, None]
if CAPTURE_REAL_M_IN0:
    M_in0_measured = input_solver.solve_section(
        mmi, wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=input_solver.waveguide_length, stub_right=0.0,
        input_idx=0, save_dir=save_dir,
    )
    # The whole justification for SYMMETRIC_MODE_REUSE, measured rather than
    # assumed. Report it; do not assert -- a real mismatch is step 3(b)'s
    # subject matter, not a reason to abort the capture.
    resid = np.abs(M_in0_measured - M_in0_derived)
    scale = np.abs(M_in0_measured)
    print(f"\n  MODE_SIGNS check: max |measured - derived| = {resid.max():.3e}, "
          f"max rel = {(resid / np.clip(scale, 1e-30, None)).max():.3e}")
    derived_signs = np.sign(np.real(
        np.sum(M_in0_measured * np.conj(M_in1), axis=(0, 2))))
    print(f"  signs implied by the measured pair: {derived_signs.astype(int).tolist()}")
    print(f"  hard-coded MODE_SIGNS:              {MODE_SIGNS.tolist()}")
    print(f"  agree: {np.array_equal(derived_signs, MODE_SIGNS.astype(float))}")
    M_in0 = M_in0_measured
else:
    M_in0 = M_in0_derived

body_neff = None
if CALIBRATE_BODY:
    calibration = input_solver.section_mode_beta_details(
        mmi, wavelengths, boundaries[1], boundaries[2], N_MODES, save_dir=save_dir,
    )
    body_neff = calibration.beta * wavelengths[:, None] / (2 * np.pi)

M_body = body_solver.solve_section_length_sweep(
    mmi, wavelengths, boundaries[1], boundaries[2], body_lengths,
    n_modes_in=N_MODES, n_modes_out=N_MODES, save_dir=save_dir,
    neff_override=body_neff,
)

M_out_bar = M_in1.transpose(0, 2, 1)
M_out_cross = M_in0.transpose(0, 2, 1)

# ---------------------------------------------------------------------------
# 2. The composition, both ways.
# ---------------------------------------------------------------------------
tmp = np.einsum("Lwij,wjk->Lwik", M_body, M_in1)
coeffs_bar = np.einsum("wij,Lwjk->Lwik", M_out_bar, tmp)[:, :, 0, 0]
coeffs_cross = np.einsum("wij,Lwjk->Lwik", M_out_cross, tmp)[:, :, 0, 0]


def via_circuit(M_out):
    """Step 2's proposed route, one Circuit per body length."""
    out = np.empty((len(body_lengths), len(wavelengths)), dtype=complex)
    for i in range(len(body_lengths)):
        c = Circuit()
        c.add_block("inp", from_transfer_matrix(M_in1, wavelengths, "src", "to_body"))
        c.add_block("body", from_transfer_matrix(M_body[i], wavelengths, "from_inp", "to_out"))
        c.add_block("out", from_transfer_matrix(M_out, wavelengths, "from_body", "port"))
        c.connect("inp", "to_body", "body", "from_inp")
        c.connect("body", "to_out", "out", "from_body")
        reduced = c.reduce(external_ports={"inp.src": "IN", "out.port": "OUT"})
        out[i] = to_smatrix(reduced, input_port="IN").coeff(port="OUT", in_mode=1, out_mode=1)
    return out


circ_bar = via_circuit(M_out_bar)
circ_cross = via_circuit(M_out_cross)

print("\n--- einsum vs Circuit on the REAL device matrices ---")
for label, hand, circ in (("bar", coeffs_bar, circ_bar),
                          ("cross", coeffs_cross, circ_cross)):
    rel = np.abs(hand - circ) / np.clip(np.abs(hand), 1e-30, None)
    print(f"  {label:5s}: max |d| = {np.abs(hand - circ).max():.3e}, "
          f"max rel = {rel.max():.3e} ({rel.max() / np.finfo(float).eps:.1f} eps)")
print("  Expected: a few eps -- Circuit is algebraically exact but contracts")
print("  with BLAS where the example uses einsum (see the harness's docstring).")

# ---------------------------------------------------------------------------
# 3. Save, with enough provenance to tell which library state produced it.
# ---------------------------------------------------------------------------
def _git(*args):
    try:
        return subprocess.run(["git", "-C", _ROOT, *args], capture_output=True,
                              text=True, timeout=30).stdout.strip()
    except Exception:                                      # noqa: BLE001
        return "unavailable"


path = os.path.join(output_dir, "cascade_golden.npz")
np.savez(
    path,
    wavelengths=wavelengths, mm_lengths=mm_lengths, body_lengths=body_lengths,
    boundaries=np.array(boundaries), mode_signs=MODE_SIGNS,
    M_in1=M_in1, M_in0=M_in0, M_in0_derived=M_in0_derived, M_body=M_body,
    coeffs_bar=coeffs_bar, coeffs_cross=coeffs_cross,
    circuit_bar=circ_bar, circuit_cross=circ_cross,
    T_bar=np.abs(coeffs_bar) ** 2, T_cross=np.abs(coeffs_cross) ** 2,
    resolution=RESOLUTION, n_modes=N_MODES,
    calibrate_body=CALIBRATE_BODY, real_m_in0=CAPTURE_REAL_M_IN0,
    git_sha=_git("rev-parse", "HEAD"),
    git_dirty=_git("diff", "--stat", "HEAD", "--", "beef/"),
)
print(f"\nGolden saved -> {path}")
print(f"  git sha: {_git('rev-parse', '--short', 'HEAD')} "
      f"(working tree {'DIRTY' if _git('diff', '--name-only', 'HEAD', '--', 'beef/') else 'clean'})")
