"""
Waisted 2x2 MMI — ALT hybrid cascade: the width RAMPS go to FDTD, not EME.

The sibling `mmi_2x2_hybrid.py` cascades

    [fdtd, eme, eme, eme, eme, eme, fdtd]
     port  stub ramp waist ramp stub port

with the whole body in EME. This script keeps the same device, the same
design point and the same comparison, and changes ONE thing:

    [fdtd, eme, fdtd, eme, fdtd, eme, fdtd]
     port  stub ramp  waist ramp  stub port

The two width ramps are now solved by FDTD. HybridSolver is solver-agnostic
("any solver type, in any order, any count >= 2"), so this is a swap of two
list entries plus the bookkeeping below — no new library code.

WHY. The ramps are the only part of this device the EME body cannot treat
exactly. A uniform region needs exactly one EME cell and propagation through
it is analytic; a ramp has to be STAIRCASED into C_TRANS cells, and every
cell boundary sheds a little |E|^2 into the interface overlap. That staircase
is the leading suspect for the method error the sibling script measures — its
IL overshoots full-device FDTD and climbs with cell count, so C_TRANS is an
OPTIMUM rather than a floor, which is not how a converging discretisation
behaves. Replacing the staircase with real FDTD removes that mechanism
entirely, and the difference between the two scripts is the measurement.

⚠ THIS IS THE EXPENSIVE VARIANT, by roughly an order of magnitude.
`FDTDSolver.solve_section()` needs ONE FDTD RUN PER EXCITED MODE, and it
excites from whichever boundary carries fewer modes —
`min(n_modes_in, n_modes_out)`. The port sections are 1-mode on their outside
face, so they cost ONE run each. A ramp section is N_MODES on BOTH faces, so
it costs N_MODES runs. Per wavelength:

    sibling  : 1 FDTD run  (port)         + EME everywhere else
    this one : 1 FDTD run  (port)         + N_MODES runs (one ramp)

The second ramp is FREE: it is the x-mirror of the first, declared as
`mirror_pairs={4: 2}` and re-validated geometrically on every call, so it is
never simulated. Likewise the output port stays mirror-derived from the input
port. Still, at N_MODES=10 this is ~11 FDTD runs per wavelength against the
sibling's 1. Shard it across wavelengths (see below) and expect real compute.

⚠ MPI ranks matter HERE, unlike in the sibling. That script is dominated by
MPB mode solves, which meep runs REDUNDANTLY on every rank (8 ranks 3110.7 s
vs 24 ranks 3055.2 s — 3x the cores for 1.8 %), so extra ranks were pure
waste. This script is dominated by FDTD TIME-STEPPING, which scales ~98 % in
this repo. Give it ranks.

THE STUB, AND HOW IT IS ACCOUNTED FOR. An FDTD section must not have its
boundaries where the width is varying: the modal basis at a seam is only
well defined on a uniform cross-section, and the section pads its drawn
waveguide by `dpml` past each boundary so the PML never sees an abrupt
core->cladding step. So each FDTD ramp section is EXTENDED by `FDTD_STUB` µm
of uniform waveguide at both ends —

    FDTD ramp section = [ramp_start - FDTD_STUB, ramp_end + FDTD_STUB]
                         |__ uniform 2.80 __|      |__ uniform 2.24 __|

— and the two EME regions it borders are SHORTENED by exactly the same
amount, so the seven sections still tile the device with no gap and no
double-count. At the default FDTD_STUB = 1.0 µm and this design point:

    section  span (µm)          length   solver  what
    0        -15.00 .. -10.75    4.25    FDTD    input port + taper + body stub
    1        -10.75 ..  -9.00    1.75    EME     uniform 2.80  (3.00 - stub)
    2         -9.00 ..  -2.00    7.00    FDTD    1.00 + ramp 5.00 + 1.00
    3         -2.00 ..   2.00    4.00    EME     uniform 2.24  (6.00 - 2*stub)
    4          2.00 ..   9.00    7.00    FDTD    x-mirror of section 2, not solved
    5          9.00 ..  10.75    1.75    EME     uniform 2.80
    6         10.75 ..  15.00    4.25    FDTD    x-mirror of section 0, not solved
                                -------
                                 30.00 = mmi.total_length, exactly.

Note this is a DIFFERENT stub from `stub_left`/`stub_right` on
solve_section(), which HybridSolver sets only for the first and last sections
of the chain (hybrid.py, from `solver.waveguide_length`) purely to give the
outermost PML clearance. Those extend BEYOND the cascade and are de-embedded;
this one is INSIDE it and is paid for by the neighbours.

MODE COUNTS. `N_FACE` is carried at the 2.80 µm seams, `N_WAIST` at the
2.24 µm ones. The waist supports exactly 10 distinct modes (N=11 is an exact
duplicate, overlap 1.000000); the 2.80 µm face supports ~14. They are
separate knobs here because a ramp section sees one of each, and because
raising N_WAIST beyond 10 buys duplicates while costing a whole extra FDTD
run per wavelength.

Run:
    # seconds, fabricated numbers, exercises the whole post-solve path
    MMI2X2T_ALT_DRYRUN=1 python .../mmi_2x2_hybrid_alt.py

    # cheap smoke: coarse mesh, 1 wavelength, 3 modes -> 3 ramp runs
    MMI2X2T_ALT_RES=12 MMI2X2T_ALT_NWL=1 MMI2X2T_ALT_NMODES=3 \
        python .../mmi_2x2_hybrid_alt.py

    MMI2X2T_ALT_WLIDX=5 python .../mmi_2x2_hybrid_alt.py   # one wavelength
    MMI2X2T_ALT_MERGE=1 python .../mmi_2x2_hybrid_alt.py   # combine the shards
"""

import functools
import glob
import hashlib
import json
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import (MMI, EMESolver, FDTDSolver, HybridSolver, MirrorSpec,
                  SMatrix, SOI_220, plot_hybrid_geometry)

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEVICE_OUT = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_2x2_tapered"))
output_dir = os.path.join(_DEVICE_OUT, "mmi_2x2_hybrid_alt")
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_sim_cache"))
# The full-device FDTD reference, written by mmi_2x2_fdtd.py in this folder.
FDTD_REF_DIR = os.path.join(_DEVICE_OUT, "mmi_2x2_fdtd")
# The EME-body cascade this variant is being compared against, if it has run.
EME_HYBRID_DIR = os.path.join(_DEVICE_OUT, "mmi_2x2_hybrid")


def _env_f(name, default):
    return float(os.environ.get(name, default))


def _env_i(name, default):
    return int(os.environ.get(name, default))


def is_master():
    """True on rank 0 (or when not running under MPI at all).

    ⚠ EVERY MPI RANK EXECUTES THIS WHOLE SCRIPT. Anything that writes a file
    or that a human is meant to read once must be guarded by this, or 24 ranks
    race on the same path and print 24 copies of every line. Measured
    2026-09-14 on job 13204906_5: all 24 ranks rendered and wrote the same
    geometry PNG concurrently. Per-rank log lines are also why a count read
    out of one of these logs has to be divided by the rank count before it
    means anything.
    """
    try:
        import meep as _mp
        return _mp.am_master()
    except Exception:
        return True


def mprint(*a, **kw):
    """print() on rank 0 only."""
    if is_master():
        print(*a, **kw)


# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by both sub-solvers
WL_MIN = 1.530      # wavelength sweep start (µm)
WL_MAX = 1.565      # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points — a subset of the reference's 51
N_MODES = 10        # modes at the 2.24 µm waist seams — 10 distinct, N=11 duplicates
BODY_STUB = 0.25    # µm of MMI body folded into each FDTD PORT section
FDTD_STUB = 1.0     # µm of uniform waveguide added each side of an FDTD RAMP
CALIBRATE_BODY = False   # HARD-REFUSED for a waisted body — see the docstring
SYMMETRIC_MODE_REUSE = True   # derive the cross path from bar via MODE_SIGNS parity

RESOLUTION = _env_f("MMI2X2T_ALT_RES", RESOLUTION)
N_WL = _env_i("MMI2X2T_ALT_NWL", N_WL)
N_MODES = _env_i("MMI2X2T_ALT_NMODES", N_MODES)
FDTD_STUB = _env_f("MMI2X2T_ALT_STUB", FDTD_STUB)
DPML = _env_f("MMI2X2T_ALT_DPML", 0.5)
# The FDTD sections are all short, so 1e-11 is reachable on them — unlike the
# 31 µm whole-device solve in mmi_2x2_fdtd.py, which floors and diverges there
# and uses 1e-05. The two scripts differing on this is deliberate, not drift.
DECAY_TOL = _env_f("MMI2X2T_ALT_DECAYTOL", 1e-11)
# RAMP vs PORT tolerance — SEPARATE KNOBS, but MEASURED to want the SAME value.
#
# A ramp section is 7.00 µm and multimode on BOTH faces, where a port section
# is 4.25 µm and single-mode on its outer face, so it is reasonable to expect
# the ramp's field to take far longer to decay. MEASURED 2026-09-14 on the
# exact ramp section [-9.0, -2.0] (res 12, 2 modes, 1 wavelength), sweeping
# the tolerance and diffing the transfer matrix between consecutive values:
#
#     decay_tol   seconds   max|dM| vs the previous row
#       1e-04       424      --
#       1e-05       632      8.4e-04
#       1e-06       869      8.5e-06
#       1e-07       436      5.7e-07
#       1e-09       476      7.3e-07
#       1e-11       650      3.8e-10      <- converged, and NOT slower
#
# Two things fall out, and both are the opposite of the obvious guess:
#   * 1e-11 IS REACHABLE here. This device's history has a hard-floor result
#     for decay_tol (1e-11 unreachable on a 84 µm 4x4: floors, then diverges),
#     and a first look at a stalled smoke run was wrongly attributed to it.
#     A 7 µm section is nowhere near that regime.
#   * COST IS FLAT IN THE TOLERANCE — 424 s at 1e-04 against 650 s at 1e-11,
#     and NON-MONOTONIC in between. The stopping condition is dominated by
#     the field having to physically leave the section, not by the threshold.
#     So a looser tolerance buys nothing and costs 8.4e-04 of accuracy.
#
# Hence both default to 1e-11. The knob stays separate because the ramp and
# the port genuinely are different problems and a future geometry may split
# them; it is not separate because this one needed it.
RAMP_DECAY_TOL = _env_f("MMI2X2T_ALT_RAMPTOL", 1e-11)
SYMMETRIC_MODE_REUSE = _env_i("MMI2X2T_ALT_REUSE", 1) == 1
DRYRUN = _env_i("MMI2X2T_ALT_DRYRUN", 0) == 1
MEMO_ON = _env_i("MMI2X2T_ALT_MEMO", 1) == 1
SHARD = os.environ.get("MMI2X2T_ALT_WLIDX", "")
MERGED = os.environ.get("MMI2X2T_ALT_MERGE", "") == "1"
# Modes at the 2.80 µm face seams. Defaults to N_MODES so a single knob still
# moves everything, but it is separate because a ramp section straddles one of
# each and the two faces genuinely support different counts.
N_FACE = _env_i("MMI2X2T_ALT_NFACE", N_MODES)
N_WAIST = N_MODES

# --- the design point (identical to both siblings in this folder) ----------
WAIST_WIDTH = 2.24       # Ww — interior body width (µm)
WAIST_LENGTH = 6.0       # Lw — uniform waist run (µm)
TRANSITION_LENGTH = 5.0  # Lt — each width ramp (µm)
MM_LENGTH = 22.0         # L  — total body length (µm)

_DRY = "_DRYRUN" if DRYRUN else ""
_RU = "" if SYMMETRIC_MODE_REUSE else "_noreuse"
_geo = (f"Ww{WAIST_WIDTH:g}_Lw{WAIST_LENGTH:g}"
        f"_Lt{TRANSITION_LENGTH:g}_L{MM_LENGTH:g}")
# An output name carries every parameter the result depends on, GEOMETRY
# FIRST. This device's study lost results to filename collisions eight times.
_TAG = (f"{_geo}_res{RESOLUTION:g}_nwl{N_WL}_nm{N_MODES}"
        f"_nf{N_FACE}_stub{FDTD_STUB:g}_rtol{RAMP_DECAY_TOL:g}{_RU}")

# ---------------------------------------------------------------------------
# 1. Geometry
# ---------------------------------------------------------------------------
mmi = MMI(
    n_inputs=2,
    n_outputs=2,
    mm_length=MM_LENGTH,
    mm_width=2.8,
    body_waist_width=WAIST_WIDTH,
    body_waist_length=WAIST_LENGTH,
    body_transition_length=TRANSITION_LENGTH,
    waveguide_width=0.45,
    taper_length=4.0,
    taper_width=1.3,
    waveguide_offset=1.5,
    waveguide_height=0.22,
)

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
s_stub = (MM_LENGTH - 2 * TRANSITION_LENGTH - WAIST_LENGTH) / 2   # 3.0 µm
ramp_x0 = -half_body + s_stub          # -8.0 — ramp starts here
ramp_x1 = -WAIST_LENGTH / 2            # -3.0 — ramp ends here

# --- the FDTD_STUB budget, refused here rather than 200 frames inside meep --
# Four independent bounds, all on the same knob:
_lim = {
    "EME stub stays positive":   s_stub - BODY_STUB,        # 2.75
    "EME waist stays positive":  WAIST_LENGTH / 2,          # 3.00
    "left PML pad stays on the uniform body face": s_stub - DPML,   # 2.50
    "right PML pad stays inside the waist": WAIST_LENGTH - DPML,    # 5.50
}
_worst = min(_lim.values())
if not (0.0 < FDTD_STUB <= _worst):
    raise SystemExit(
        f"REFUSING: FDTD_STUB={FDTD_STUB} µm is outside (0, {_worst:.2f}].\n  "
        + "\n  ".join(f"{k}: <= {v:.2f}" for k, v in sorted(_lim.items(),
                                                            key=lambda kv: kv[1]))
        + "\n  The FDTD ramp section is extended by FDTD_STUB at BOTH ends and\n"
          "  the neighbouring EME regions are shortened by the same amount, so\n"
          "  too large a stub eats a neighbour; too large also walks this\n"
          "  section's own PML padding onto varying-width geometry.")

boundaries = [-half_total,
              -half_body + BODY_STUB,       # 0 FDTD port | 1 EME stub
              ramp_x0 - FDTD_STUB,          # 1 EME stub  | 2 FDTD ramp
              ramp_x1 + FDTD_STUB,          # 2 FDTD ramp | 3 EME waist
              -(ramp_x1 + FDTD_STUB),       # 3 EME waist | 4 FDTD ramp
              -(ramp_x0 - FDTD_STUB),       # 4 FDTD ramp | 5 EME stub
              +half_body - BODY_STUB,       # 5 EME stub  | 6 FDTD port
              +half_total]
_len = [boundaries[i + 1] - boundaries[i] for i in range(7)]

mprint("=" * 78)
mprint("WAISTED 2x2 — ALT HYBRID  [fdtd, eme, fdtd, eme, fdtd, eme, fdtd]")
mprint("=" * 78)
mprint(f"  Ww={WAIST_WIDTH:g} Lw={WAIST_LENGTH:g} Lt={TRANSITION_LENGTH:g} "
      f"L={MM_LENGTH:g}   res={RESOLUTION:g}")
mprint(f"  decay_tol: ports {DECAY_TOL:g} (cache-shared with the sibling), "
      f"ramps {RAMP_DECAY_TOL:g} (swept: converged, and no slower)")
mprint(f"  {N_WL} wl, N_WAIST={N_WAIST} N_FACE={N_FACE}, "
      f"FDTD_STUB={FDTD_STUB:g} µm (max {_worst:.2f})")
_who = ["FDTD port", "EME stub", "FDTD RAMP", "EME waist",
        "FDTD RAMP (mirror)", "EME stub", "FDTD port (mirror)"]
for i in range(7):
    mprint(f"    {i}  [{boundaries[i]:+7.2f}, {boundaries[i+1]:+7.2f}]  "
          f"{_len[i]:5.2f} µm   {_who[i]}")
mprint(f"    sum {sum(_len):.2f} µm  vs  mmi.total_length {mmi.total_length:.2f} µm")
assert abs(sum(_len) - mmi.total_length) < 1e-9, "sections do not tile the device"
for _x in boundaries[1:-1]:
    assert abs(mmi.width_at(_x - 1e-3) - mmi.width_at(_x + 1e-3)) < 1e-9, (
        f"section boundary at x={_x} sits where the width is varying")
mprint("  every seam verified to sit on a uniform cross-section")
mprint(f"  calibrate_body_neff={CALIBRATE_BODY} (REFUSED for a waisted body)")
mprint(f"  SYMMETRIC_MODE_REUSE={SYMMETRIC_MODE_REUSE}")
mprint("=" * 78, flush=True)

# ---------------------------------------------------------------------------
# 2. Per-mode parity vector for the mirror-derived cross path
# ---------------------------------------------------------------------------
_MODE_SIGNS_14 = np.array([1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, 1])
if N_FACE > len(_MODE_SIGNS_14):
    raise SystemExit(
        f"N_FACE={N_FACE} exceeds the {len(_MODE_SIGNS_14)} known MODE_SIGNS "
        "entries; re-derive with MMI2X2T_ALT_REUSE=0.")
# Entry k is a property of mode k alone, and asking for fewer modes returns
# the same lowest ones in the same order, so truncation is exact.
MODE_SIGNS = _MODE_SIGNS_14[:N_FACE]

# ---------------------------------------------------------------------------
# 3. Content-addressed mode memo (see the sibling script for the full note).
#    Keyed on the actual index arrays, never on width: width-keyed reuse was
#    measured BREAKING PHYSICS on the 1x4 because equal width is not equal
#    cross-section once a cell sits off-axis. functools.wraps, not a subclass.
# ---------------------------------------------------------------------------
_MODE_MEMO, _MEMO_STATS = {}, {"hits": 0, "misses": 0}


def _install_mode_memo():
    orig = EMESolver.compute_modes

    def _key(cs, n_modes):
        h = hashlib.blake2b(digest_size=16)
        for comp in ("nx", "ny", "nz"):
            arr = np.ascontiguousarray(getattr(cs, comp))
            h.update(np.asarray(arr.shape, dtype=np.int64).tobytes())
            h.update(arr.tobytes())
        return (h.hexdigest(), float(cs.env.wl), float(cs.env.T), int(n_modes))

    @functools.wraps(orig)
    def compute_modes(self, cs, n_modes, *args, **kwargs):
        if not MEMO_ON:
            return orig(self, cs, n_modes, *args, **kwargs)
        try:
            k = _key(cs, n_modes)
        except Exception:
            return orig(self, cs, n_modes, *args, **kwargs)
        if k in _MODE_MEMO:
            _MEMO_STATS["hits"] += 1
            return _MODE_MEMO[k]
        _MEMO_STATS["misses"] += 1
        res = orig(self, cs, n_modes, *args, **kwargs)
        _MODE_MEMO[k] = res
        return res

    EMESolver.compute_modes = compute_modes


_install_mode_memo()

# ---------------------------------------------------------------------------
# 4. Sub-solvers
# ---------------------------------------------------------------------------
_fdtd_common = dict(res=RESOLUTION, sim_width=5.0, sim_height=1.5, dpml=DPML,
                    source_margin=0.0, monitor_margin=0.25, three_d=True,
                    material_stack=SOI_220, decay_tol=DECAY_TOL)

# The two PORT sections keep half_domain_aperture: each has a neighbour arm
# 1.5 µm away in the same drawing, and the aperture structurally excludes it.
# The RAMP sections must NOT use it — they are body sections centred on y=0
# with no second arm to exclude, and halving their aperture would truncate the
# very body modes they exist to propagate.
port_in_solver = FDTDSolver(half_domain_aperture=True, **_fdtd_common)
port_out_solver = FDTDSolver(half_domain_aperture=True, **_fdtd_common)
_ramp_common = dict(_fdtd_common, decay_tol=RAMP_DECAY_TOL)
ramp_in_solver = FDTDSolver(**_ramp_common)
ramp_out_solver = FDTDSolver(**_ramp_common)


def _eme(length):
    """One cell: every EME section in this chain is now a UNIFORM region, and
    one cell is exact for a uniform region — propagation through it is
    analytic, and extra cells would only insert identity interfaces between
    identical cross-sections. The staircase that C_TRANS controlled in the
    sibling script has no counterpart here; that is the point of the variant."""
    return EMESolver(num_modes=[N_MODES], num_cells=[1], lengths=[length],
                     oxide_width=5.0, oxide_height=1.5, res=RESOLUTION,
                     material_stack=SOI_220, mode_solver="meep")


solvers = [port_in_solver,      # 0  FDTD  port + taper + body stub
           _eme(_len[1]),       # 1  EME   uniform 2.80
           ramp_in_solver,      # 2  FDTD  stub + ramp + stub
           _eme(_len[3]),       # 3  EME   uniform 2.24 (waist)
           ramp_out_solver,     # 4  FDTD  x-mirror of 2 — never simulated
           _eme(_len[5]),       # 5  EME   uniform 2.80
           port_out_solver]     # 6  FDTD  x-mirror of 0 — never simulated
OUT_SECTION = len(solvers) - 1

# Per-interface mode counts: the 2.80 µm faces and the 2.24 µm waist support
# different numbers, and a ramp section straddles one of each.
interface_modes = [N_FACE,    # -10.75  port | stub      2.80
                   N_FACE,    #  -9.00  stub | ramp      2.80
                   N_WAIST,   #  -2.00  ramp | waist     2.24
                   N_WAIST,   #  +2.00  waist | ramp     2.24
                   N_FACE,    #  +9.00  ramp | stub      2.80
                   N_FACE]    # +10.75  stub | port      2.80

_shared = dict(solvers=solvers, boundaries=boundaries,
               interface_modes=interface_modes,
               calibrate_body_neff=CALIBRATE_BODY)

# Section 4 is a PLAIN x-mirror of section 2 — legal because both are body
# sections centred on y=0, so mirroring x maps one onto the other without any
# y-flip. (The OUTPUT PORT is the case that needs MirrorSpec: its cross path
# is top-in/bottom-out, which a plain x-mirror cannot express.)
_ramp_mirror = {4: 2}
if SYMMETRIC_MODE_REUSE:
    cross_mirror = dict(_ramp_mirror)
    cross_mirror[OUT_SECTION] = MirrorSpec(section=0,
                                           partner_kwargs={"input_idx": 0},
                                           mode_signs=MODE_SIGNS)
else:
    # Honest "simulate both arms": no mirror on the output port at all, so it
    # is solved directly on its own geometry. The ramp mirror is kept — it is
    # a pure x-mirror of an on-axis section and is not in question.
    cross_mirror = dict(_ramp_mirror)
bar_mirror = dict(_ramp_mirror)
bar_mirror[OUT_SECTION] = 0

cross_kwargs = {0: {"input_idx": 1}, OUT_SECTION: {"output_idx": 0}}
bar_kwargs = {0: {"input_idx": 1}, OUT_SECTION: {"output_idx": 1}}

hybrid_bar = HybridSolver(mirror_pairs=bar_mirror, **_shared)
hybrid_cross = HybridSolver(mirror_pairs=cross_mirror, **_shared)

_runs = 1 + min(N_FACE, N_WAIST)
mprint(f"cost: ~{_runs} FDTD runs per wavelength "
      f"(1 port + {min(N_FACE, N_WAIST)} ramp); sections 4 and 6 are mirrors "
      f"and are never simulated", flush=True)

# ---------------------------------------------------------------------------
# 5. Geometry plot — shows which solver covers which section
# ---------------------------------------------------------------------------
if not MERGED and is_master():
    try:
        plot_hybrid_geometry(
            mmi, solver=hybrid_bar, in_modes=1, out_modes=1,
            title="Waisted 2x2 MMI — ALT Hybrid Sections (FDTD ramps)",
            output_path=os.path.join(output_dir, "geometry_hybrid_alt.png"),
            sub_region_dir=os.path.join(output_dir, "sections"),
            sub_region_show_pml=True, sub_region_show_eme_detail=True,
            sub_region_configs=[("in_top", {0: {"input_idx": 1}}, bar_mirror)],
        )
        mprint(f"Geometry plot saved -> {output_dir}/geometry_hybrid_alt.png")
    except Exception as exc:                                # pragma: no cover
        mprint(f"⚠ geometry plot skipped: {type(exc).__name__}: {exc}")

# ---------------------------------------------------------------------------
# 6. Solve
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
if SHARD != "" and not MERGED:
    k = int(SHARD)
    if not 0 <= k < N_WL:
        raise SystemExit(f"MMI2X2T_ALT_WLIDX={k} out of range for N_WL={N_WL}")
    wavelengths = wavelengths[k:k + 1]
    mprint(f"[shard] wavelength index {k} of {N_WL}: {wavelengths[0]:.4f} µm",
          flush=True)


def _shard_path(k):
    return os.path.join(output_dir, f"shard{_DRY}_wl{k}_{_TAG}.json")


if MERGED:
    shards = []
    for f in sorted(glob.glob(os.path.join(output_dir,
                                           f"shard{_DRY}_wl*_{_TAG}.json"))):
        txt = open(f).read().lstrip()
        try:
            doc, end = json.JSONDecoder().raw_decode(txt)
        except (json.JSONDecodeError, ValueError) as exc:
            raise SystemExit(
                f"MERGE: {os.path.basename(f)} is not readable JSON ({exc}). "
                "A zero-byte or truncated shard is debris from a crashed run "
                "— delete it and re-run that wavelength.")
        if txt[end:].strip():
            mprint(f"  note: {os.path.basename(f)} had trailing bytes from a "
                  "multi-rank write race; using the leading document")
        shards.append(doc)
    have = sorted(sh["wl_index"] for sh in shards)
    if have != list(range(N_WL)):
        raise SystemExit(
            f"MERGE: need wavelength shards {list(range(N_WL))}, have {have}. "
            "Run the missing shards first — a partial merge would silently "
            "compare on fewer wavelengths than the reference.")
    shards.sort(key=lambda sh: sh["wl_index"])
    wavelengths = np.array([sh["wavelength"] for sh in shards])
    c_bar = np.array([complex(*sh["c_bar"]) for sh in shards])
    c_cross = np.array([complex(*sh["c_cross"]) for sh in shards])
    t_bar = sum(sh["seconds"]["bar"] for sh in shards)
    t_cross = sum(sh["seconds"]["cross"] for sh in shards)
    mprint(f"MERGE: {len(shards)} wavelength shards, "
          f"{t_bar + t_cross:.0f}s of solve summed across them "
          "(they ran in parallel, so wall time was far less)")
elif DRYRUN:
    # Fabricate the coefficients and skip the solve, so the ENTIRE post-solve
    # path can be exercised in seconds. A bug in exactly that path — treating
    # HybridSolver.solve()'s SMatrix as a raw array — was once only reached
    # after 10.5 h of solving.
    mprint("DRYRUN: fabricating coefficients; NO physics in this run")
    rng = np.random.default_rng(0)
    c_bar = (rng.normal(size=len(wavelengths))
             + 1j * rng.normal(size=len(wavelengths))) * 0.5
    c_cross = (rng.normal(size=len(wavelengths))
               + 1j * rng.normal(size=len(wavelengths))) * 0.5
    t_bar = t_cross = 0.0
else:
    # ONE shared cache across both excitations: section 0's geometry_kwargs
    # are identical between the bar and cross calls (both excite the top arm),
    # and sections 1-5 do not depend on the excitation at all, so the whole
    # expensive ramp is solved exactly ONCE for the two of them.
    section_cache = {}
    t0 = time.time()
    sm_bar = hybrid_bar.solve(mmi, wavelengths, geometry_kwargs=bar_kwargs,
                              cache=section_cache, save_dir=save_dir)
    t_bar = time.time() - t0
    mprint(f"[hybrid] bar solved in {t_bar:.0f}s", flush=True)
    t0 = time.time()
    sm_cross = hybrid_cross.solve(mmi, wavelengths,
                                  geometry_kwargs=cross_kwargs,
                                  cache=section_cache, save_dir=save_dir)
    t_cross = time.time() - t0
    mprint(f"[hybrid] cross solved in {t_cross:.0f}s", flush=True)
    # ⚠ HybridSolver.solve() returns an SMatrix, NOT a raw array.
    c_bar = np.asarray(sm_bar.coeff(port="out", in_mode=1, out_mode=1),
                       dtype=complex).reshape(len(wavelengths))
    c_cross = np.asarray(sm_cross.coeff(port="out", in_mode=1, out_mode=1),
                         dtype=complex).reshape(len(wavelengths))

if SHARD != "" and not MERGED:
    # ⚠ ONLY ONE RANK MAY WRITE. Every MPI rank reaches this line; when all of
    # them wrote the same path their writes interleaved into concatenated JSON.
    if not is_master():
        raise SystemExit(0)
    sp = _shard_path(int(SHARD))
    tmp = sp + ".partial"          # atomic: a crash must not leave a 0-byte shard
    with open(tmp, "w") as fh:
        json.dump({"wl_index": int(SHARD),
                   "wavelength": float(wavelengths[0]),
                   "n_wl_total": N_WL,
                   "config": {"res": RESOLUTION, "n_modes": N_MODES,
                              "n_face": N_FACE, "fdtd_stub": FDTD_STUB,
                              "decay_tol": DECAY_TOL,
                              "ramp_decay_tol": RAMP_DECAY_TOL,
                              "cascade": "fdtd,eme,fdtd,eme,fdtd,eme,fdtd"},
                   "design": {"Ww": WAIST_WIDTH, "Lw": WAIST_LENGTH,
                              "Lt": TRANSITION_LENGTH, "L": MM_LENGTH},
                   "seconds": {"bar": t_bar, "cross": t_cross},
                   "c_bar": [c_bar[0].real, c_bar[0].imag],
                   "c_cross": [c_cross[0].real, c_cross[0].imag]}, fh, indent=2)
    os.replace(tmp, sp)
    mprint(f"[shard] written -> {os.path.basename(sp)}")
    h, m = _MEMO_STATS["hits"], _MEMO_STATS["misses"]
    if h + m:
        mprint(f"[shard] mode memo: {m} solved, {h} reused "
              f"({100*h/(h+m):.0f}% avoided)")
    mprint("[shard] merge with MMI2X2T_ALT_MERGE=1 once every shard has landed.")
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# 7. Scores, and the comparison against full-device FDTD
# ---------------------------------------------------------------------------
CONTROL_PHASE_DEG = -89.74957705211592   # the straight 2x2's measured quadrature


def metrics(a_bar, a_cross):
    Tb, Tc = np.abs(a_bar) ** 2, np.abs(a_cross) ** 2
    tot = Tb + Tc
    dphi = np.angle(a_cross / a_bar, deg=True)
    return {"T_bar": Tb, "T_cross": Tc, "total": tot, "phase_per_wl": dphi,
            "excess_loss_db": float(np.max(-10*np.log10(np.clip(tot, 1e-20, None)))),
            "imbalance": float(np.max(np.abs(Tb - Tc))),
            "phase_deg": float(np.mean(dphi)),
            "phase_max_dev_deg": float(np.max(np.abs(dphi - CONTROL_PHASE_DEG))),
            "phase_spread_deg": float(np.ptp(dphi))}


mh = metrics(c_bar, c_cross)

ref_path = os.path.join(FDTD_REF_DIR, f"smatrix_waisted_{_geo}.pkl")
if not os.path.exists(ref_path):
    raise SystemExit(
        f"no full-device FDTD reference for this geometry:\n  {ref_path}\n"
        "Run mmi_2x2_fdtd.py in this same folder first — it writes that file.")
sm_ref = SMatrix.load(ref_path)
fb = sm_ref.coeff(port="out_0", in_mode=1, out_mode=1)
fc = sm_ref.coeff(port="out_1", in_mode=1, out_mode=1)
# Match each hybrid wavelength to its nearest reference wavelength and REFUSE
# if any pair is not essentially identical — comparing two engines at
# wavelengths that merely nearly agree would fold a wavelength offset into
# the engine disagreement being measured.
ref_wl = np.linspace(WL_MIN, WL_MAX, len(fb))
idx = [int(np.argmin(np.abs(ref_wl - w))) for w in wavelengths]
mismatch = float(np.max(np.abs(ref_wl[idx] - wavelengths)))
if mismatch > 1e-9:
    raise SystemExit(
        f"wavelength grids do not align: worst |dwl| = {mismatch:.3e} µm "
        f"(hybrid {len(c_bar)} pts vs reference {len(fb)}). Use an N_WL whose "
        f"grid is a subset of the reference's, e.g. N_WL=11 against 51.")
fb, fc = fb[idx], fc[idx]
mf = metrics(fb, fc)
mprint(f"\n[compare] matched {len(idx)} wavelengths to the stored reference "
      f"(worst |dwl| = {mismatch:.1e} µm)")

dTb = np.abs(mh["T_bar"] - mf["T_bar"])
dTc = np.abs(mh["T_cross"] - mf["T_cross"])

mprint("\n" + "=" * 78)
mprint("ALT HYBRID vs FULL-DEVICE FDTD — identical geometry and wavelengths")
mprint("=" * 78)
mprint(f"{'port':>6} {'max|dT|':>10} {'mean|dT|':>10} {'% of mean T_fdtd':>17} "
      f"{'worst wl':>10}")
for name, d, ref in (("bar", dTb, mf["T_bar"]), ("cross", dTc, mf["T_cross"])):
    i = int(np.argmax(d))
    mprint(f"{name:>6} {d.max():10.4f} {d.mean():10.4f} "
          f"{100*d.mean()/ref.mean():17.2f} {wavelengths[i]:10.4f}")
mprint(f"\n{'quantity':>18} {'alt hybrid':>11} {'full FDTD':>10} {'delta':>10}")
for q in ("excess_loss_db", "imbalance", "phase_deg", "phase_max_dev_deg",
          "phase_spread_deg"):
    mprint(f"{q:>18} {mh[q]:11.4f} {mf[q]:10.4f} {mh[q]-mf[q]:+10.4f}")

# --- the actual point of this variant: alt vs the EME-body cascade ---------
# Both are cascades of the SAME device against the SAME reference, differing
# only in how the two ramps are solved, so this row is the measurement.
_eme_json = sorted(glob.glob(os.path.join(
    EME_HYBRID_DIR, f"hybrid_vs_fdtd_{_geo}_res{RESOLUTION:g}_*.json")))
if _eme_json:
    with open(_eme_json[-1]) as fh:
        _e = json.load(fh)
    mprint("\n" + "-" * 78)
    mprint(f"ALT (FDTD ramps) vs SIBLING (EME staircase ramps), both vs full FDTD")
    mprint(f"  sibling file: {os.path.basename(_eme_json[-1])}")
    mprint(f"{'quantity':>18} {'alt':>11} {'sibling':>11} {'full FDTD':>11}")
    for q in ("excess_loss_db", "imbalance", "phase_deg", "phase_max_dev_deg"):
        mprint(f"{q:>18} {mh[q]:11.4f} {_e['hybrid'][q]:11.4f} {mf[q]:11.4f}")
    mprint(f"{'max|dT| bar':>18} {dTb.max():11.4f} "
          f"{_e['max_abs_dT']['bar']:11.4f} {'--':>11}")
    mprint(f"{'max|dT| cross':>18} {dTc.max():11.4f} "
          f"{_e['max_abs_dT']['cross']:11.4f} {'--':>11}")
    mprint("-" * 78)
else:
    mprint(f"\n(no sibling EME-body result found in {EME_HYBRID_DIR} to compare "
          "against; run mmi_2x2_hybrid.py for the head-to-head)")

# ---------------------------------------------------------------------------
# 8. Save and plot
# ---------------------------------------------------------------------------
smat = SMatrix(wavelengths=wavelengths)
smat.set(c_bar, port="bar")
smat.set(c_cross, port="cross")
smat.check_passive()
smat.save(os.path.join(output_dir, f"smatrix_hybrid_alt{_DRY}_{_TAG}.pkl"))
mprint(f"\nSMatrix saved -> smatrix_hybrid_alt{_DRY}_{_TAG}.pkl")

out_json = os.path.join(output_dir, f"hybrid_alt_vs_fdtd{_DRY}_{_TAG}.json")
with open(out_json, "w") as fh:
    json.dump({"config": {"res": RESOLUTION, "n_wl": N_WL,
                          "n_modes": N_MODES, "n_face": N_FACE,
                          "fdtd_stub": FDTD_STUB, "decay_tol": DECAY_TOL,
                          "ramp_decay_tol": RAMP_DECAY_TOL,
                          "symmetric_mode_reuse": SYMMETRIC_MODE_REUSE,
                          "calibrate_body_neff": CALIBRATE_BODY,
                          "cascade": "fdtd,eme,fdtd,eme,fdtd,eme,fdtd",
                          "section_bounds": boundaries,
                          "section_lengths": _len},
               "design": {"Ww": WAIST_WIDTH, "Lw": WAIST_LENGTH,
                          "Lt": TRANSITION_LENGTH, "L": MM_LENGTH},
               "seconds": {"bar": t_bar, "cross": t_cross},
               "wavelengths": wavelengths.tolist(),
               "hybrid": {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                          for k, v in mh.items()},
               "full_fdtd": {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                             for k, v in mf.items()},
               "max_abs_dT": {"bar": float(dTb.max()),
                              "cross": float(dTc.max())}}, fh, indent=2)
mprint(f"Comparison saved -> {os.path.basename(out_json)}")

fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))

ax = axes[0]
ax.plot(wavelengths, mh["T_bar"], color="steelblue", marker="o", ms=3,
        label="Bar — alt hybrid")
ax.plot(wavelengths, mf["T_bar"], color="steelblue", linestyle="--",
        label="Bar — full FDTD")
ax.plot(wavelengths, mh["T_cross"], color="darkorange", marker="o", ms=3,
        label="Cross — alt hybrid")
ax.plot(wavelengths, mf["T_cross"], color="darkorange", linestyle="--",
        label="Cross — full FDTD")
ax.axhline(0.5, color="gray", linestyle=":", lw=0.8)
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title("Per-port transmission")
ax.legend(fontsize=7)
ax.grid(True)
ax.set_ylim(0, 0.7)

ax = axes[1]
ax.plot(wavelengths, mh["total"], color="k", marker="o", ms=3, label="Alt hybrid")
ax.plot(wavelengths, mf["total"], color="k", linestyle="--", label="Full FDTD")
ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Total transmission")
ax.set_title("Total power — the staircase's own signature")
ax.legend(fontsize=8)
ax.grid(True)

# ⚠ |T| DISCARDS WHAT THIS DEVICE IS FOR — the bar/cross quadrature is the
# function, and an intensity plot cannot see it.
ax = axes[2]
ax.plot(wavelengths, mh["phase_per_wl"], color="seagreen", marker="o", ms=3,
        label="Alt hybrid")
ax.plot(wavelengths, mf["phase_per_wl"], color="seagreen", linestyle="--",
        label="Full FDTD")
ax.axhline(CONTROL_PHASE_DEG, color="gray", linestyle=":", lw=0.8,
           label=f"Straight control ({CONTROL_PHASE_DEG:.2f} deg)")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("arg(cross / bar) (deg)")
ax.set_title("Quadrature")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle(f"Waisted 2x2 MMI — ALT hybrid (FDTD ramps) vs full-device FDTD "
             f"({_geo}, res {RESOLUTION:g}, {N_MODES} modes, "
             f"stub {FDTD_STUB:g} µm)", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, f"hybrid_alt_results{_DRY}_{_TAG}.png"),
            dpi=300)
plt.close("all")
mprint(f"Results plot saved -> hybrid_alt_results{_DRY}_{_TAG}.png")
