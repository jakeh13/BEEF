"""Pure-EME evaluation of the 4x4 hybrid at a given (mm_width, y_gap).

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-25.

PURPOSE 1 (now): the y_gap bias study. Before running a Bayesian optimization
over (mm_width, y_gap) in EME, we need to know whether the EME-vs-FDTD
discrepancy DEPENDS ON y_gap. It plausibly does, and that would be
disqualifying: `EMESolver.solve()` draws only one input and one output branch
per call and treats the other six ports as physically absent, so it discards
port-to-port crosstalk entirely -- and y_gap IS the tip-to-tip separation that
sets how large that crosstalk is. Coupling is roughly exponential in gap, so
the neglected term could vary by orders of magnitude across the search range,
biasing the optimizer toward whichever end the approximation flatters.

Check 4 (2026-08-25) already measured that the approximation costs up to 10 deg
of phase and 6.6 dB of imbalance at the CURRENT geometry (y_gap=0.70). The
open question is only whether that error is roughly constant in y_gap (in
which case it largely cancels in a relative comparison and the BO is sound) or
strongly varying (in which case it does not).

PURPOSE 2 (later): this is the BO's evaluation function. Same code, so the
bias study validates exactly what the optimizer will call.

PARAMETERIZATION (from the user's y_gap sketch). Port centres sit at
`waveguide_offset * (i - 1.5)`, so adjacent port spacing IS waveguide_offset
and the tip-to-tip gap between neighbouring tapers is

    y_gap = waveguide_offset - taper_width

Holding the paper's Eq. 8 relation `waveguide_offset = mm_width/4`, two free
parameters fix the whole geometry:

    waveguide_offset = mm_width / 4
    taper_width      = mm_width / 4 - y_gap

Sanity: the current design (mm_width=6.0, taper_width=0.8) is y_gap=0.70,
which matches the 0.7 um port clearance quoted in mmi_4x4_90hybrid_eme.py's
own docstring. Validity requires taper_width > waveguide_width, i.e.
0 < y_gap < mm_width/4 - 0.35.

MODE SOLVER: tidy3d, per the user's instruction and because it is ~100x
cheaper than meep here (the meep path could not finish ONE solve() call in
12 h at N_MODES_BODY=16 -- job 12155408). tidy3d's raw fields are noisier than
meep's (fdtd_vs_EME_mode_validation), but n_eff/envelope agree closely and
that is what a cascade needs.

Run:
    MMI4X4_YGAP=0.4 MMI4X4_MMWIDTH=6.0 python mmi_4x4_90hybrid_ygap_eme.py
"""
import os
import time

import numpy as np

from beef import MMI, EMESolver, SMatrix, SOI_220

Y_GAP = float(os.environ.get("MMI4X4_YGAP", 0.70))
MM_WIDTH = float(os.environ.get("MMI4X4_MMWIDTH", 6.0))
# Bias study fixes the length at the current design point so EME and FDTD are
# compared on IDENTICAL geometry. The BO will instead sweep this via
# solve_length_sweep(), which reuses one mode solve across every length.
MM_LENGTH = float(os.environ.get("MMI4X4_MMLENGTH", 0.75 * 108.9195 - 2.75))
# Wide band by default: the existing FDTD reference (jobs 12272571/12272573)
# is at 1.26-1.36, and reusing it saves two ~14 h runs. The BO itself will use
# 1.31 +- 0.02.
WL_MIN = float(os.environ.get("MMI4X4_WLMIN", 1.26))
WL_MAX = float(os.environ.get("MMI4X4_WLMAX", 1.36))
N_WL = int(os.environ.get("MMI4X4_NWL", 3))

WAVEGUIDE_WIDTH = 0.35
TAPER_LENGTH = 1.0
RES = float(os.environ.get("MMI4X4_EVAL_RES", 40.0))
N_MODES_TAPER = int(os.environ.get("MMI4X4_EVAL_NTAPER", 12))
N_MODES_BODY = int(os.environ.get("MMI4X4_EVAL_NBODY", 16))
N_CELLS_TAPER = int(os.environ.get("MMI4X4_EVAL_NCELLS", 30))

WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = WAVEGUIDE_OFFSET - Y_GAP
if not (WAVEGUIDE_WIDTH < TAPER_WIDTH < WAVEGUIDE_OFFSET):
    raise SystemExit(
        f"y_gap={Y_GAP} invalid for mm_width={MM_WIDTH}: gives taper_width="
        f"{TAPER_WIDTH:.3f}, needs {WAVEGUIDE_WIDTH} < taper_width < "
        f"{WAVEGUIDE_OFFSET}. Valid y_gap range: "
        f"(0, {WAVEGUIDE_OFFSET - WAVEGUIDE_WIDTH:.3f})")

DRIVEN = [0, 2]
OUTPUTS = [0, 1, 2, 3]
PORT_NAMES = ["out5", "out6", "out7", "out8"]
TARGET_DEG = {0: -45.0, 1: 225.0, 2: 45.0, 3: 135.0}

_HERE = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                        "example_sim_cache"))
out_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..",
                                       "example_output", "mmi_4x4_90hybrid",
                                       "ygap_study"))
os.makedirs(out_dir, exist_ok=True)

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

mmi = MMI(
    n_inputs=4, n_outputs=4, mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.22, use_neighbor_aperture=True,
)

solver = EMESolver(
    num_modes=[1, N_MODES_TAPER, N_MODES_BODY, N_MODES_TAPER, 1],
    num_cells=[1, N_CELLS_TAPER, 1, N_CELLS_TAPER, 1],
    lengths=[0.5, TAPER_LENGTH, MM_LENGTH, TAPER_LENGTH, 0.5],
    oxide_width=MM_WIDTH + 3.0, oxide_height=2.0,
    res=RES, material_stack=SOI_220, mode_solver="tidy3d",
)


def metrics(s1, s3):
    """The paper's metrics, same construction mmi_4x4_90hybrid_hybrid.py uses.

    A11's corrected self/other order on cross_excitation_phase_diff_deg
    (port3 as self, port1 as other) -- do not reorder.
    """
    m = {}
    for lab, s in (("port1", s1), ("port3", s3)):
        m[f"imbalance_{lab}"] = np.asarray(s.imbalance_db(PORT_NAMES))
        m[f"excess_loss_{lab}"] = np.asarray(s.excess_loss_db(PORT_NAMES))
        m[f"cmrr_I_{lab}"] = np.asarray(s.cmrr_db("out5", "out8"))
        m[f"cmrr_Q_{lab}"] = np.asarray(s.cmrr_db("out6", "out7"))
    phi = {j: np.asarray(s3.cross_excitation_phase_diff_deg(s1, port=PORT_NAMES[j]))
           for j in OUTPUTS}
    for j in (1, 2, 3):
        raw = phi[j] - phi[0] - (TARGET_DEG[j] - TARGET_DEG[0])
        m[f"dphi_{j+5}"] = ((raw + 180.0) % 360.0) - 180.0
    return m


def objective(m, w_phase=0.5, w_cmrr=0.5):
    """Weighted, normalised, lower-is-better.

    Each term is scaled by its own spec so the two are commensurate: phase
    against the 5 deg bar, CMRR against the 20 dBe bar. The CMRR term is
    one-sided -- exceeding spec earns no credit, it just stops penalising.
    """
    max_dphi = max(float(np.abs(m[f"dphi_{j}"]).max()) for j in (6, 7, 8))
    min_cmrr = min(float(m[k].min()) for k in m if k.startswith("cmrr"))
    return (w_phase * (max_dphi / 5.0)
            + w_cmrr * max(0.0, (20.0 - min_cmrr) / 20.0)), max_dphi, min_cmrr


def main():
    print(f"YGAP_EME: y_gap={Y_GAP:.3f} mm_width={MM_WIDTH:.3f} -> "
          f"offset={WAVEGUIDE_OFFSET:.3f} taper_width={TAPER_WIDTH:.3f}",
          flush=True)
    print(f"YGAP_EME: mm_length={MM_LENGTH:.4f} res={RES} "
          f"N_TAPER={N_MODES_TAPER} N_BODY={N_MODES_BODY} "
          f"cells={N_CELLS_TAPER} tidy3d", flush=True)
    print(f"YGAP_EME: {N_WL} wavelengths {WL_MIN}-{WL_MAX} um", flush=True)

    t0 = time.time()
    smat = {}
    for i in DRIVEN:
        s = SMatrix(wavelengths=wavelengths)
        for j in OUTPUTS:
            tc = time.time()
            r = solver.solve(mmi, wavelengths, ports=[f"out{5+j}_in{i}"],
                             input_idx=i, output_idx=j, save_dir=save_dir)
            s.set(np.asarray(r.coeff(port=f"out{5+j}_in{i}")), port=PORT_NAMES[j])
            print(f"YGAP_EME:   in{i}->out{5+j} {time.time()-tc:.1f}s", flush=True)
        smat[i] = s
    print(f"YGAP_EME: 8 solves in {(time.time()-t0)/60:.1f} min", flush=True)

    m = metrics(smat[0], smat[2])
    obj, max_dphi, min_cmrr = objective(m)
    print(f"\nYGAP_EME: max|dphi|={max_dphi:.3f} deg  min CMRR={min_cmrr:.3f} dBe")
    print(f"YGAP_EME: objective={obj:.4f}  (lower is better; 0 = both specs met)")

    np.savez(os.path.join(out_dir, f"eme_ygap{Y_GAP:.2f}_W{MM_WIDTH:.2f}.npz"),
             y_gap=Y_GAP, mm_width=MM_WIDTH, taper_width=TAPER_WIDTH,
             mm_length=MM_LENGTH, wavelengths=wavelengths,
             objective=obj, max_dphi=max_dphi, min_cmrr=min_cmrr,
             **{k: v for k, v in m.items()})
    print(f"YGAP_EME: saved -> {out_dir}/eme_ygap{Y_GAP:.2f}_W{MM_WIDTH:.2f}.npz")


if __name__ == "__main__":
    main()
