"""1x4 MMI on 300 nm Si3N4 -- body_stub (interface cut-position) sweep.

validation_markdown/mmi_1x4_sin300_validation.md, opened 2026-09-03 against
the hybrid-vs-FDTD gap.

WHAT QUESTION THIS ANSWERS
--------------------------
The calibrated hybrid disagrees with the whole-device FDTD reference by 5.99%
worst-case per port (lambda=1.50, outer pair), 2.78% mean, against a 1% target.
The error splits into two independent pieces:

  common mode (total T):  -3.37% at 1.50 -> +0.76% at 1.60, crossing zero at 1.59
  differential (split):   at 1.50 the outer pair is -5.98% while the inner pair
                          is -0.75%, i.e. ~5.2 pp of pure imbalance error

Zeroing the common mode still leaves ~2.6% per port, so the differential is the
binding constraint. This script targets it.

The validation doc's standing suspect for the engine disagreement -- the EME
taper staircase, C_TAPER=8 -- CANNOT explain the hybrid's error: mmi_1x4_hybrid.py
sets INPUT_METHOD = OUTPUT_METHOD = "fdtd", so the hybrid's port tapers are meep
at the same res and domain as the whole-device reference. C_TAPER never enters
this chain. Only four things actually differ between the hybrid and the
reference, and `body_stub` -- where the two cut planes sit -- is the one that was
inherited from mmi_1x2_hybrid.py/mmi_2x2_hybrid.py and never re-derived for this
platform (mmi_1x4_hybrid.py's own comment flags it).

HYPOTHESIS. At body_stub=0.25 the interface plane sits 0.25 um past a 2.2 -> 12.0
um width step. Evanescent and radiative content is largest there, and it is
exactly the content the 12-mode projection discards. Pushing the cut deeper into
the body lets it decay before the projection. This predicts the differential
shrinks with depth and plateaus.

WHY THIS IS ONE RUN PER SECTION AND NOT ONE RUN PER STUB
-------------------------------------------------------
`FDTDSolver.solve_section(measure_margins=[...])` puts one mode monitor per
listed along-axis offset into the SAME simulation. sim.run()'s cost is dominated
by stop_when_dft_decayed()'s decay wait, a property of the whole field and
independent of how many DFT monitors are recording it, so N positions cost
roughly what 1 does (fdtd.py's own docstring). Measured precedent on the 2x2:
51 monitors over a 2 nm grid in a single 19:06:00 run at res=30, 26.5% of its
72 h budget (job 11995441).

So the expensive stage is 5 FDTD runs -- 1 input section + 4 output sections --
regardless of how many stub values are swept. The cascade for every stub is then
free arithmetic in stage B.

HOW A STUB VALUE MAPS ONTO A MONITOR OFFSET
-------------------------------------------
solve_section()'s measurement monitor sits exactly at the interface plane
(section_fdtd_layout(): `mea_x = L`, the far end) -- `monitor_margin` positions
the EXCITED side's normalization monitor, not this one. With an explicit
measure_port the position becomes
    mea_along = (hi_ax if sign == +1 else lo_ax) - sign * measure_margin
(fdtd.py's PortSpec branch). Pin each FDTD domain at the DEEPEST stub and the
offset that lands the monitor at body depth d is the same for both ends:

    input  section: domain [-half_total, -half_body + STUB_MAX]
                    monitor wanted at -half_body + d, high edge, sign=+1
                    => measure_margin = STUB_MAX - d
    output section: domain [+half_body - STUB_MAX, +half_total]
                    monitor wanted at +half_body - d, low edge, sign=-1
                    => measure_margin = STUB_MAX - d

WHY BOTH PORTS ARE EXPLICIT PortSpecs ON THE OUTPUT SECTIONS
------------------------------------------------------------
Passing measure_port at all sets solve_section()'s `override` flag, which forces
`reverse=False` -- excite at the LOW edge, measure at the HIGH edge. That is
already what the input section wants. It is backwards for an output section,
whose body-side plane is its LOW edge, so those calls must ALSO pass excite_port
to put the source back on the arm. The excite aperture is taken from the
geometry's own source_aperture_at() hook -- the same one section_fdtd_layout()
calls -- rather than hand-written, so use_neighbor_aperture stays honoured.

Exciting the arm and measuring into the body yields the TRANSPOSE of the output
section's cascade matrix. That is the same reciprocity mmi_2x2_hybrid.py's cross
path and beef.MirrorSpec already rely on, and it is what the 2x2's own
body_stub sweep did (job 11995441: "top-arm excitation").

⚠ UNVERIFIED LIBRARY PATH. fdtd.py's own comment on the measure-side PortSpec
branch says the LOW/HIGH convention there is "physically-motivated, not verified
against the reference scripts -- they never override the measure side's position
via an explicit PortSpec; only excite_port's convention is directly verified".
This script is therefore the first real user of it, which is why stage B opens
with a reproduction gate (below) rather than trusting the geometry.

THE GATE
--------
STUBS includes the production 0.25. Cascading that column must reproduce
mmi_1x4_hybrid.py's shipped answer -- ports [0.22804, 0.24468, 0.24467, 0.22805],
total 0.9454 at 1.55 um (job 12227591). It will NOT be bit-identical: the FDTD
domain here is pinned at STUB_MAX, so the section is a physically larger
simulation than production's, and the monitor reads at the same plane inside a
different domain. Agreement to GATE_RTOL is the bar. A failure here means the
reconstruction is wrong -- wrong edge, wrong sign, wrong transpose -- and every
other column is meaningless, so it aborts rather than writing a plausible file.

BETA IS REUSED ACROSS EVERY STUB, DELIBERATELY
----------------------------------------------
The body's calibrated per-mode beta is a property of the 12 um slab cross-section,
which is identical at every cut depth; only the propagated LENGTH changes
(mm_length - 2d). So stage B calls section_mode_beta_details() once at the
PRODUCTION span (-19.70, +19.70) -- a cache hit against the prewarm that job
array 12031838 already paid ~3 days for -- and reuses that beta for every stub.
Re-measuring per stub would change the calibration span, miss the cache, and
re-pay the calibration N times over to obtain the same number.
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, FDTDSolver, SIN_300, SMatrix
from beef.solvers.fdtd import PortSpec

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "mmi_1x4", "mmi_1x4_body_stub_sweep")
save_dir = os.path.join(_ROOT, "example_sim_cache")
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Config -- every value pinned to mmi_1x4_hybrid.py's production run.
#    Any drift here silently compares two different devices.
# ---------------------------------------------------------------------------
RESOLUTION = 30.0
WL_MIN, WL_MAX, N_WL = 1.50, 1.60, 11
N_MODES = int(os.environ.get("MMI1X4_N_MODES", 12))

WAVEGUIDE_WIDTH = 1.1
TAPER_WIDTH = 2.2
TAPER_LENGTH = 5.0
WAVEGUIDE_OFFSET = 3.1
MM_WIDTH = 12.0
MM_LENGTH = 39.90

FDTD_SIM_WIDTH = MM_WIDTH + 4.0     # 16.0
FDTD_SIM_HEIGHT = 4.0

PRODUCTION_STUB = 0.25
# Depths to probe, um into the body past the taper/body discontinuity. The
# production point is included as the gate; STUB_MAX sets how far each pinned
# FDTD domain reaches, and therefore the deepest cut reachable in one run.
STUBS = np.array([0.25, 0.40, 0.55, 0.70, 0.85, 1.00])
# Each domain reaches HEADROOM past the deepest cut, so the deepest monitor
# still sits HEADROOM inside the domain edge rather than on it. At margin=0 the
# monitor lands exactly where the PML begins, which the aperture-vs-PML guard
# rejects outright -- and a monitor grazing PML reads as loss while still
# looking plausible (the ring's 27% power loss at |T| 0.711 vs 0.978). Set to
# monitor_margin's own 0.25 so the deepest column is no closer to PML than
# every production run's monitor already is.
HEADROOM = 0.25
STUB_MAX = float(STUBS.max()) + HEADROOM

# Gate tolerance against the shipped production numbers. Loose enough to absorb
# the genuine domain-size difference described in the module docstring, tight
# enough that a wrong edge/sign/transpose (which moves T by tens of percent, not
# tenths) cannot slip through.
GATE_RTOL = 0.02
PRODUCTION_T_1550 = np.array([0.22804, 0.24468, 0.24467, 0.22805])

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

mmi = MMI(
    n_inputs=1, n_outputs=4,
    mm_length=MM_LENGTH,
    mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH,
    taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH,
    waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.30,     # pinned to mmi_1x4_hybrid.py -- omitting it takes
                               # MMI's own default, a different device
    use_neighbor_aperture=True,
)

half_body = mmi.mm_length / 2.0
half_total = mmi.total_length / 2.0

solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=FDTD_SIM_WIDTH,
    sim_height=FDTD_SIM_HEIGHT,
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,
    material_stack=SIN_300,
    decay_tol=1e-11,
)

# Same offset list serves both ends -- see the module docstring's mapping.
measure_margins = [float(STUB_MAX - d) for d in STUBS]
assert min(measure_margins) >= HEADROOM - 1e-12, (
    "every monitor must sit at least HEADROOM inside its domain edge -- see "
    "the STUB_MAX comment on why margin=0 is not acceptable"
)

# The measured (body-side) plane is a full-width multimode interface: mea_w =
# sim_width centred at y=0, matching what section_fdtd_layout() picks itself
# whenever n_measure > 1 and monitor_size is None. Written out explicitly here
# because an explicit PortSpec never inherits that default.
body_port_in = PortSpec(axis="x", sign=+1, center=0.0, width=FDTD_SIM_WIDTH)
body_port_out = PortSpec(axis="x", sign=-1, center=0.0, width=FDTD_SIM_WIDTH)

print(f"body_stub sweep: depths {list(STUBS)} um, domains pinned at STUB_MAX={STUB_MAX}")
print(f"measure_margins = {measure_margins}")
print(f"{len(wavelengths)} wavelengths, {N_MODES} interface modes\n")

# ---------------------------------------------------------------------------
# 2. Stage A -- the five FDTD runs. Each writes its own .npz as soon as it
#    lands, so a walltime kill after run 3 keeps runs 1-2 instead of losing the
#    lot (a JSON-at-end sweep on this cluster has already thrown away 8 h of
#    finished work once).
# ---------------------------------------------------------------------------
def _stage_a_path(tag: str) -> str:
    return os.path.join(output_dir, f"section_{tag}.npz")


def solve_input_section() -> np.ndarray:
    """[n_wl, n_pos, N_MODES, 1] -- arm excited, body-side plane measured."""
    path = _stage_a_path("input")
    if os.path.exists(path):
        print(f"input section: reusing {path}")
        return np.load(path)["T"]
    x0, x1 = -half_total, -half_body + STUB_MAX
    print(f"input section: x=[{x0:.3f}, {x1:.3f}], measuring at the HIGH edge")
    # stub_left=waveguide_length reproduces HybridSolver._stubs(0), which
    # extends the FIRST section's domain by the solver's own 0.5 um lead-in.
    # Omitting it would simulate a shorter device than production does.
    T = solver.solve_section(
        mmi, wavelengths, x0, x1,
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=solver.waveguide_length, stub_right=0.0,
        measure_port=body_port_in,
        measure_margins=measure_margins,
        save_dir=save_dir,
    )
    np.savez(path, T=T, stubs=STUBS, wavelengths=wavelengths,
             measure_margins=measure_margins, x0=x0, x1=x1)
    print(f"input section -> {path}")
    return T


def solve_output_section(k: int) -> np.ndarray:
    """[n_wl, n_pos, N_MODES, 1] -- arm k excited, body-side plane measured.

    This is the TRANSPOSE of branch k's cascade matrix; see the module
    docstring on reciprocity.
    """
    path = _stage_a_path(f"output{k}")
    if os.path.exists(path):
        print(f"output {k}: reusing {path}")
        return np.load(path)["T"]
    x0, x1 = half_body - STUB_MAX, half_total
    # The arm sits at the HIGH edge and launches in -x, so excite sign = -1
    # (an excite_port's sign=-1 IS the high edge -- the opposite convention to
    # a measure_port's, per solve_section()'s docstring). Aperture from the
    # geometry's own hook so use_neighbor_aperture is honoured rather than
    # re-guessed.
    arm_x = x1 - 1e-6
    arm_center = mmi.center_at(arm_x, output_idx=k)
    arm_aperture = mmi.source_aperture_at(arm_x, output_idx=k)
    excite = PortSpec(axis="x", sign=-1, center=float(arm_center), width=float(arm_aperture))
    print(f"output {k}: x=[{x0:.3f}, {x1:.3f}], arm y={arm_center:+.3f} "
          f"aperture={arm_aperture:.3f}, measuring at the LOW edge")
    # stub_right here for the same reason stub_left appears on the input
    # section -- this is HybridSolver._stubs(K-1), the LAST section's trailing
    # 0.5 um of access waveguide.
    T = solver.solve_section(
        mmi, wavelengths, x0, x1,
        n_modes_in=1, n_modes_out=N_MODES,
        stub_left=0.0, stub_right=solver.waveguide_length,
        excite_port=excite,
        measure_port=body_port_out,
        measure_margins=measure_margins,
        save_dir=save_dir,
        output_idx=k,
    )
    np.savez(path, T=T, stubs=STUBS, wavelengths=wavelengths,
             measure_margins=measure_margins, x0=x0, x1=x1)
    print(f"output {k} -> {path}")
    return T


# MMI1X4_STUB_SECTION selects what this process does, so stage A's five
# INDEPENDENT section solves can run as a 5-task array (wall time = the slowest
# section) instead of sequentially in one job (wall time = their sum):
#   "input" / "0".."3"  solve that one section, write its npz, stop
#   "cascade"           skip stage A entirely, require all five npz, do stage B
#   unset               everything in one process (the login-node/serial path)
SECTION = os.environ.get("MMI1X4_STUB_SECTION", "").strip().lower()

if SECTION == "input":
    solve_input_section()
    raise SystemExit(0)
if SECTION in {"0", "1", "2", "3"}:
    solve_output_section(int(SECTION))
    raise SystemExit(0)

if SECTION == "cascade":
    missing = [t for t in ["input", "output0", "output1", "output2", "output3"]
               if not os.path.exists(_stage_a_path(t))]
    if missing:
        raise SystemExit(
            f"MMI1X4_STUB_SECTION=cascade needs every stage-A file, missing: "
            f"{missing}. Run those sections first -- do NOT let this fall through "
            f"to solving them here, which is how a 'quick' analysis job silently "
            f"becomes a multi-day FDTD run."
        )

M_in = solve_input_section()                       # [n_wl, n_pos, N_MODES, 1]
M_out_T = np.stack([solve_output_section(k) for k in range(4)])   # [4, n_wl, n_pos, N_MODES, 1]

# ---------------------------------------------------------------------------
# 3. Stage B -- free. Calibrated beta once at the PRODUCTION span (cache hit),
#    then one diagonal body propagation per stub.
# ---------------------------------------------------------------------------
cal_x0, cal_x1 = -half_body + PRODUCTION_STUB, half_body - PRODUCTION_STUB
print(f"\ncalibration: section_mode_beta_details on the production span "
      f"[{cal_x0:.3f}, {cal_x1:.3f}] (expected cache hit)")
calibration = solver.section_mode_beta_details(
    mmi, wavelengths, cal_x0, cal_x1, N_MODES,
    save_dir=save_dir,
)
beta = np.asarray(calibration.beta)                # [n_wl, N_MODES], rad/um

rows = []
jsonl_path = os.path.join(output_dir, "body_stub_sweep.jsonl")
T_all = np.zeros((len(STUBS), 4, len(wavelengths)))

with open(jsonl_path, "w") as fh:
    for si, d in enumerate(STUBS):
        body_length = MM_LENGTH - 2.0 * float(d)
        # Diagonal body transfer -- the calibrated uniform-body fast path
        # (EMESolver.solve_section()'s neff_override branch) is exactly
        # exp(1j*beta*L) per mode, which is why no EME solve is needed here.
        phase = np.exp(1j * beta * body_length)    # [n_wl, N_MODES]
        for k in range(4):
            a_in = M_in[:, si, :, 0]               # [n_wl, N_MODES]
            a_out = M_out_T[k][:, si, :, 0]        # [n_wl, N_MODES] (transposed branch)
            amp = np.sum(a_out * phase * a_in, axis=1)     # [n_wl]
            T_all[si, k] = np.abs(amp) ** 2
        total = T_all[si].sum(axis=0)
        row = dict(
            body_stub=float(d), body_length=body_length,
            T=[T_all[si, k].tolist() for k in range(4)],
            total_T=total.tolist(),
        )
        fh.write(json.dumps(row) + "\n")
        fh.flush()
        rows.append(row)
        i55 = int(np.argmin(np.abs(wavelengths - 1.55)))
        print(f"  stub={d:.2f}  L_body={body_length:.2f}  "
              f"T(1.55)={np.round(T_all[si, :, i55], 5)}  total={total[i55]:.4f}")

# ---------------------------------------------------------------------------
# 4. THE GATE -- production column must reproduce the shipped answer.
# ---------------------------------------------------------------------------
gate_idx = int(np.argmin(np.abs(STUBS - PRODUCTION_STUB)))
assert abs(STUBS[gate_idx] - PRODUCTION_STUB) < 1e-9, (
    "STUBS must contain the production 0.25 um -- it is the only column whose "
    "correct answer is already known, and without it nothing validates the "
    "reconstruction."
)
i55 = int(np.argmin(np.abs(wavelengths - 1.55)))
gate_T = T_all[gate_idx, :, i55]
rel = np.abs(gate_T / PRODUCTION_T_1550 - 1.0)
print(f"\nGATE at body_stub={PRODUCTION_STUB} um, lambda=1.55 um:")
print(f"  reconstructed {np.round(gate_T, 5)}")
print(f"  production    {PRODUCTION_T_1550}  (job 12227591)")
print(f"  worst relative deviation {rel.max():.4f}  (bar {GATE_RTOL})")
if rel.max() > GATE_RTOL:
    raise SystemExit(
        f"GATE FAILED: reconstruction deviates {rel.max():.4f} > {GATE_RTOL} from "
        f"the shipped production result. The measure-side PortSpec path is the "
        f"prime suspect (fdtd.py flags its LOW/HIGH convention as unverified) -- "
        f"check the edge, the sign, and the output transpose before trusting ANY "
        f"column of this sweep. Refusing to write a plausible wrong plot."
    )
print("  GATE PASSED\n")

np.savez(os.path.join(output_dir, "body_stub_sweep.npz"),
         stubs=STUBS, wavelengths=wavelengths, T=T_all, beta=beta)

# ---------------------------------------------------------------------------
# 5. Plot against the whole-device FDTD reference.
# ---------------------------------------------------------------------------
ref_path = os.path.join(_ROOT, "example_output", "mmi_1x4", "mmi_1x4_fdtd", "smatrix.pkl")
ref = SMatrix.load(ref_path)
# .transmission(port=...), NOT .data[...] -- SMatrix exposes no `data` attribute;
# that key layout belongs to the raw pickle, not the loaded object. The FDTD
# reference names its ports "out_0..3" (underscore) while the hybrid scripts use
# "out0..3"; both spellings exist on disk, so read them from ref.ports() rather
# than hard-coding either.
ref_ports = sorted(ref.ports())
assert len(ref_ports) == 4, f"expected 4 reference ports, got {ref_ports}"
F = np.array([np.asarray(ref.transmission(port=p)) for p in ref_ports])

fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6))
outer = lambda T: (T[0] + T[3]) / 2.0
inner = lambda T: (T[1] + T[2]) / 2.0

ax[0].plot(STUBS, [T_all[s].sum(axis=0)[i55] for s in range(len(STUBS))], "o-", label="hybrid")
ax[0].axhline(F[:, i55].sum(), color="k", ls="--", label="whole-device FDTD")
ax[0].set_xlabel("body_stub (um)"); ax[0].set_ylabel("total T @ 1.55 um")
ax[0].set_title("Common mode vs cut depth"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)

ax[1].plot(STUBS, [inner(T_all[s])[i55] - outer(T_all[s])[i55] for s in range(len(STUBS))],
           "s-", color="tab:red", label="hybrid")
ax[1].axhline(inner(F)[i55] - outer(F)[i55], color="k", ls="--", label="whole-device FDTD")
ax[1].set_xlabel("body_stub (um)"); ax[1].set_ylabel("mean(inner) - mean(outer) @ 1.55 um")
ax[1].set_title("Differential (split shape) vs cut depth"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)

worst = [100 * np.max(np.abs(T_all[s] / F - 1.0)) for s in range(len(STUBS))]
ax[2].plot(STUBS, worst, "d-", color="tab:green")
ax[2].axhline(1.0, color="gray", ls=":", label="1% target")
ax[2].set_xlabel("body_stub (um)"); ax[2].set_ylabel("worst per-port error vs FDTD (%), all wl")
ax[2].set_title("Distance to the 1% target"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)

fig.suptitle("1x4 MMI SiN -- hybrid interface cut position vs the whole-device FDTD reference")
fig.tight_layout()
fig.savefig(os.path.join(output_dir, "body_stub_sweep.png"), dpi=150)
print(f"Saved -> {output_dir}/body_stub_sweep.png")
print(f"Rows   -> {jsonl_path}")
