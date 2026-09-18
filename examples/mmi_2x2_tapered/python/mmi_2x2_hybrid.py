"""
Tutorial: simulate the WAISTED 2x2 MMI splitter with the BEEF library using a
FDTD/EME hybrid cascade.

Method: FDTD for the port sections, EME for the multimode body, stitched by
HybridSolver. This buys most of the whole-device FDTD's accuracy at a fraction
of its cost — but on a waisted body it is also the subtlest of this folder's
three scripts, for the reasons below.

Geometry: the tapered-body counterpart of
examples/mmi_2x2/python/mmi_2x2_hybrid.py. Same idea as that script, but the
body is no longer a uniform slab, and that single change invalidates three of
the shortcuts the straight device's hybrid takes. Read those first; they are
the reason this script looks different rather than merely re-parameterised.

THREE THINGS THE STRAIGHT DEVICE'S HYBRID DOES THAT THIS ONE CANNOT

1. calibrate_body_neff=True is HARD-REFUSED for a waisted body.
   assert_uniform_cross_section() fires from HybridSolver._calibrated_neff():
   "a calibrated body propagates on one per-mode beta, which only describes a
   uniform cross-section". That flag is what closed the engine-agreement gap on
   every other device here — MPB's mode index and meep's time-stepped phase
   differ ~2.3 % per mode, and on the directional coupler the flag moved L50
   from 5.37 % off to 0.11 %. So this cascade carries an UNCORRECTED body beta,
   and measuring the size of the resulting disagreement is the point of the
   script.

2. solve_length_sweep(sweep_section=1) needs a length-independent uniform
   slab, so the straight device gets a 201-point length sweep for free. The
   waisted body is neither uniform nor length-independent: this is a
   single-length solve.

3. The body goes through solve_section()'s real MULTI-CELL branch — one
   cross-section per cell, real interface overlaps, a non-diagonal transfer
   matrix. num_cells therefore MATTERS here, unlike every uniform-body script
   in this repo, which passes num_cells=[1] because one cell is exact.

THE BODY IS SIX SECTIONS, NOT ONE. solve_section() divides a section
UNIFORMLY into num_cells[0] cells; it has no notion of internal region
boundaries. Handing it the whole body therefore forces one cell size
everywhere, which over-cells the uniform waist (where one cell is exact and
extra cells only insert identity interfaces between identical cross-sections)
in order to give the ramps enough. HybridSolver is a generic N-section cascade,
so the body simply BECOMES its own regions instead:

    [FDTD in] [stub] [ramp] [waist] [ramp] [stub] [FDTD out]
               1cell  8cells  1cell  8cells  1cell

19 cells instead of 40, of which the content memo below leaves only 10 as real
mode solves.

⚠ N_MODES = 10, NOT the straight device's 14. The 2.24 µm waist supports
exactly 10 distinct modes at this wavelength — N = 11 is an EXACT duplicate
(overlap 1.000000). 14 was inherited from the 2.80 µm face and cost a 15×
runtime penalty solving modes that do not exist. Do not raise it.

⚠ C_TRANS = 8 IS AN OPTIMUM, NOT A FLOOR. Measured to ct32 on two mode bases:
imbalance and phase converge, but IL OVERSHOOTS the FDTD reference and then
climbs away. More cells is not more accurate here — each EME cell boundary
sheds ~2 % of |E|², so more boundaries means more spurious loss. Do not re-run
the ladder hoping for convergence.

WHY THIS RUNS ONE WAVELENGTH PER JOB. MPB does not scale with MPI ranks — meep
solves the eigenproblem REDUNDANTLY on every rank (measured here: 8 ranks
3110.7 s vs 24 ranks 3055.2 s, 3× the cores for 1.8 %). A mode-solve-dominated
job gains nothing from more ranks. Wavelengths ARE independent, so they
parallelise at the JOB level: one shard per wavelength with few ranks each,
then a merge. This cannot be done through the cache — the solve_section cache
key contains the whole wavelength array, so a single 11-wavelength call and
eleven 1-wavelength calls are different keys and share nothing.

Run:
    # one shard per wavelength (what the batch script does), then merge
    MMI2X2T_HY_WLIDX=0 python .../mmi_2x2_hybrid.py
    ...
    MMI2X2T_HY_WLIDX=10 python .../mmi_2x2_hybrid.py
    MMI2X2T_HY_MERGE=1 python .../mmi_2x2_hybrid.py

    # exercise the whole post-solve path in seconds, with NO physics
    MMI2X2T_HY_DRYRUN=1 MMI2X2T_HY_WLIDX=0 python .../mmi_2x2_hybrid.py
"""

import glob
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
output_dir = os.path.join(_DEVICE_OUT, "mmi_2x2_hybrid")
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_sim_cache"))
# The full-device FDTD reference this cascade is checked against, written by
# mmi_2x2_fdtd.py in this same folder.
FDTD_REF_DIR = os.path.join(_DEVICE_OUT, "mmi_2x2_fdtd")


def _env_f(name, default):
    return float(os.environ.get(name, default))


def _env_i(name, default):
    return int(os.environ.get(name, default))


# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by both the EME and FDTD sub-solvers
WL_MIN = 1.530      # wavelength sweep start (µm)
WL_MAX = 1.565      # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points — a subset of the reference's 51
N_MODES = 10        # ⚠ 10, not the straight device's 14 — see the module docstring
C_TRANS = 8         # EME cells per width ramp — an OPTIMUM, not a floor
BODY_STUB = 0.25    # µm of MMI body folded into each FDTD section
CALIBRATE_BODY = False   # HARD-REFUSED for a waisted body — see the docstring
SYMMETRIC_MODE_REUSE = True   # derive the cross path from bar via MODE_SIGNS parity

RESOLUTION = _env_f("MMI2X2T_HY_RES", RESOLUTION)
N_WL = _env_i("MMI2X2T_HY_NWL", N_WL)
N_MODES = _env_i("MMI2X2T_HY_NMODES", N_MODES)
C_TRANS = _env_i("MMI2X2T_HY_CTRANS", C_TRANS)
DPML = _env_f("MMI2X2T_HY_DPML", 0.5)
# The FDTD PORT sections are short, and 1e-11 IS reachable on them — unlike the
# whole-device solve in mmi_2x2_fdtd.py, which floors and diverges there and so
# uses 1e-05. The two scripts differing on this is deliberate, not drift.
DECAY_TOL = _env_f("MMI2X2T_HY_DECAYTOL", 1e-11)
SYMMETRIC_MODE_REUSE = _env_i("MMI2X2T_HY_REUSE", 1) == 1
DRYRUN = _env_i("MMI2X2T_HY_DRYRUN", 0) == 1
MEMO_ON = _env_i("MMI2X2T_HY_MEMO", 1) == 1

# ---------------------------------------------------------------------------
# 0. Simulation method per section — HybridSolver is solver-agnostic (any mix
#    of FDTDSolver/EMESolver, in any order, any count >= 2).
# ---------------------------------------------------------------------------
INPUT_METHOD = "fdtd"    # "fdtd" or "eme"
BODY_METHOD = "eme"      # "fdtd" or "eme"
OUTPUT_METHOD = "fdtd"   # "fdtd" or "eme"

# ---------------------------------------------------------------------------
# 1. Geometry — the validated waist design point, identical to mmi_2x2_fdtd.py
#    in this folder. Comparing two engines on anything else would fold a
#    geometry difference into the engine disagreement being measured.
# ---------------------------------------------------------------------------
WAIST_WIDTH = 2.24       # Ww
WAIST_LENGTH = 6.0       # Lw
TRANSITION_LENGTH = 5.0  # Lt
MM_LENGTH = 22.0         # L

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
    waveguide_height=0.22,  # SOI_220
)
# ⚠ `MMI.total_length` is mm_length + both tapers and stops there; it does NOT
# include the 0.5 µm access stubs the criterion's device length counts (31.00 µm
# here). It IS the right quantity for the section boundaries below, which live
# in the drawn-geometry frame — the two numbers are different on purpose.
print(f"MMI drawn span: {mmi.total_length:.2f} µm "
      f"(device length by the criterion's convention: "
      f"{mmi.mm_length + 9.0:.2f} µm)")
print(f"Input ports at y = {mmi.input_ys} µm; outputs at y = {mmi.output_ys} µm")

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2
stub = (MM_LENGTH - 2 * TRANSITION_LENGTH - WAIST_LENGTH) / 2

# Fail HERE, with the reason, rather than 200 frames deep inside meep.
_s_min = BODY_STUB + DPML
if stub < _s_min - 1e-9:
    raise SystemExit(
        f"REFUSING: uniform stub s={stub:.3f} µm < {_s_min:.3f} µm.\n"
        f"  A hybrid FDTD section pads its drawn waveguide by dpml={DPML} µm\n"
        f"  past the interior seam, so the PML never sees an abrupt\n"
        f"  core->cladding step. With body_stub={BODY_STUB} that padding\n"
        f"  reaches {_s_min:.2f} µm into the body from the seam; if the waist\n"
        f"  ramp starts before that, MMI._check_cut() refuses — a section\n"
        f"  boundary cannot sit where the width is varying.\n"
        f"  Need s = (L - 2*Lt - Lw)/2 >= body_stub + dpml = {_s_min:.2f} µm.\n"
        f"  This is why the retracted 29.50 µm design point (s = 0.25) was\n"
        f"  never hybrid-buildable at all.")

# Per-mode sign relating the bottom-arm launch to the top-arm launch
# (M_in0[m] = MODE_SIGNS[m] * M_in1[m]). Derived on the STRAIGHT device's input
# section, and reused here for a measured reason rather than an assumed one:
# the input FDTD section's polygons are byte-identical between the straight and
# waisted devices (section_gds() re-origins to x_start=0), and the vector is a
# property of that section's mirror symmetry, not of the body.
_MODE_SIGNS_14 = np.array([1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, 1])
if N_MODES > len(_MODE_SIGNS_14):
    raise SystemExit(
        f"N_MODES={N_MODES} exceeds the {len(_MODE_SIGNS_14)} known MODE_SIGNS "
        "entries; re-derive with MMI2X2T_HY_REUSE=0 (parity beyond mode 13 was "
        "never measured, and guessing it would silently corrupt the second "
        "arm's coefficients)")
# TRUNCATION IS VALID. MODE_SIGNS is a per-mode PARITY vector — entry k is a
# property of mode k alone. Asking the solver for 10 modes instead of 14
# returns the SAME lowest 10 modes in the same order, so entries 0..9 still
# describe them. Parity does not depend on how many modes were requested.
MODE_SIGNS = _MODE_SIGNS_14[:N_MODES]

# ---------------------------------------------------------------------------
# 2. Content-addressed mode memo
#
# The waisted body is MIRROR-SYMMETRIC about x=0 and uniform through the waist,
# so most of its cells present a cross-section that has already been solved:
# cell i sits at x and its partner at -x, and inside the body the cross-section
# is just a slab of width w(x) = w(-x) centred on y=0 — IDENTICAL, not merely
# similar.
#
# ⚠ KEYED ON CONTENT, NOT ON WIDTH, and not on the symmetry argument itself.
# This repo has already been burned by the tempting version of this idea:
# width-keyed mode reuse was measured BREAKING PHYSICS on the 1x4 (total T
# 0.9703 -> 0.1507), because equal width does NOT mean equal cross-section once
# a cell sits off-axis. The key here is a hash of the actual refractive-index
# arrays (nx, ny, nz) plus wavelength, temperature and mode count, so reuse is
# exact by construction: identical eigenproblem in, identical modes out. If the
# symmetry argument were wrong the hashes simply would not match and every cell
# would be solved for real — the failure mode is "no speedup", never "wrong
# answer". Measured elsewhere in this study at |dS| = 0 exactly.
# ---------------------------------------------------------------------------
_MODE_MEMO = {}
_MEMO_STATS = {"hits": 0, "misses": 0}


def _install_mode_memo():
    """Wrap EMESolver.compute_modes with a content-addressed memo.

    functools.wraps, NOT a subclass: subclassing a solver in this repo breaks
    cache fingerprints and silently drops cache_only/save_dir forwarding.
    """
    import functools
    import hashlib
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
        except Exception:          # unknown cross-section shape -> just solve
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
print(f"mode memo: {'ON' if MEMO_ON else 'OFF'} "
      "(content-addressed on nx/ny/nz + wl + T + n_modes)")

# ---------------------------------------------------------------------------
# 3. Configure sub-solvers
# ---------------------------------------------------------------------------
_fdtd_common = dict(
    res=RESOLUTION,
    sim_width=5.0,
    sim_height=1.5,
    dpml=DPML,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,       # full 3-D FDTD: real 220 nm Si core in oxide, PML on all six faces
    material_stack=SOI_220,
    decay_tol=DECAY_TOL,
)


def make_section_solver(method, half_domain_aperture=False,
                        n_cells=1, length=1.0):
    """Build a solver for one hybrid section. Both FDTDSolver and EMESolver
    implement the same solve_section() interface, so either can fill any slot
    in the chain."""
    if method == "fdtd":
        return FDTDSolver(half_domain_aperture=half_domain_aperture,
                          **_fdtd_common)
    if method == "eme":
        # ⚠ num_cells IS read by solve_section() on a NON-UNIFORM section —
        # unlike every uniform-body script in this repo, which passes [1].
        return EMESolver(num_modes=[N_MODES], num_cells=[n_cells],
                         lengths=[length], oxide_width=5.0, oxide_height=1.5,
                         res=RESOLUTION, material_stack=SOI_220,
                         mode_solver="meep")
    raise ValueError(f"unknown method {method!r}, expected 'fdtd' or 'eme'")


# Both sides have two closely-spaced ports (waveguide_offset = 1.5 µm), so both
# the input and output FDTD sections need half_domain_aperture to confine each
# section's source/norm aperture to the near half of the domain, structurally
# excluding the opposite branch.
input_solver = make_section_solver(INPUT_METHOD, half_domain_aperture=True)
output_solver = make_section_solver(OUTPUT_METHOD, half_domain_aperture=True)

# ---------------------------------------------------------------------------
# 4. Build HybridSolver — seven sections here (see the module docstring)
# ---------------------------------------------------------------------------
ramp_x0 = -half_body + stub          # ramp starts after the uniform stub
waist_x0, waist_x1 = -WAIST_LENGTH / 2, +WAIST_LENGTH / 2
boundaries = [
    -half_total,                     # physical input port (input taper tip)
    -half_body + BODY_STUB,          # input FDTD section <-> body
    ramp_x0,                         # stub -> ramp
    waist_x0,                        # ramp -> waist
    waist_x1,                        # waist -> ramp
    -ramp_x0,                        # ramp -> stub
    +half_body - BODY_STUB,          # body -> output FDTD section
    +half_total,                     # physical output port (output taper tip)
]
stub_in_len = ramp_x0 - (-half_body + BODY_STUB)   # stub minus the folded part

# Each UNIFORM region gets ONE cell because one cell is EXACT for a uniform
# region — propagation is analytic and extra cells only insert identity
# interfaces between identical cross-sections. Only the ramps need the
# staircase.
solvers = [
    input_solver,
    make_section_solver(BODY_METHOD, n_cells=1, length=stub_in_len),         # stub, uniform 2.80
    make_section_solver(BODY_METHOD, n_cells=C_TRANS, length=TRANSITION_LENGTH),  # ramp
    make_section_solver(BODY_METHOD, n_cells=1, length=WAIST_LENGTH),        # waist, uniform 2.24
    make_section_solver(BODY_METHOD, n_cells=C_TRANS, length=TRANSITION_LENGTH),  # ramp, mirror
    make_section_solver(BODY_METHOD, n_cells=1, length=stub_in_len),         # stub, mirror
    output_solver,
]
N_SECTIONS = len(solvers)
OUT_SECTION = N_SECTIONS - 1
n_cells_total = 1 + C_TRANS + 1 + C_TRANS + 1
print(f"  body regions: stub {stub_in_len:.2f}µm(1) | ramp {TRANSITION_LENGTH:g}µm({C_TRANS}) "
      f"| waist {WAIST_LENGTH:g}µm(1) | ramp {TRANSITION_LENGTH:g}µm({C_TRANS}) "
      f"| stub {stub_in_len:.2f}µm(1)  = {n_cells_total} cells, "
      f"{C_TRANS + 2} distinct after the memo")

_shared = dict(
    solvers=solvers,
    boundaries=boundaries,
    interface_modes=[N_MODES] * (N_SECTIONS - 1),
    # ⚠ REFUSED for a waisted body — this is not a tuning choice. See the
    # module docstring; the uncorrected body beta is what the comparison at the
    # end of this script measures.
    calibrate_body_neff=CALIBRATE_BODY,
)

if SYMMETRIC_MODE_REUSE:
    # ⚠ NOTE THE SECTION INDEX. With the body split into five regions the
    # output section is index 6, not 2. Getting this wrong would mirror a body
    # region onto the input section instead of the output one.
    #
    # The cross path pairs one arm's transposed matrix with the OTHER arm's
    # excitation, which solve()'s one-call contract cannot express as a plain
    # mirror. MirrorSpec says it explicitly: reuse section 0, reinterpreted as
    # the other arm (partner_kwargs) via the measured parity vector.
    cross_mirror = {OUT_SECTION: MirrorSpec(section=0,
                                            partner_kwargs={"input_idx": 0},
                                            mode_signs=MODE_SIGNS)}
else:
    # The honest "simulate both arms" configuration is NO mirror at all: solve
    # the output section directly on its own geometry, costing one extra FDTD
    # section solve.
    #
    # ⚠ It is NOT `{OUT_SECTION: 0}` with input_idx=0. That was the code here
    # until 2026-09-10, and it was both wrong and unable to express the cross
    # path — the two faults hid each other. input_idx=0 -> output_idx=0 excites
    # the BOTTOM arm and reads the BOTTOM output, which on this y-symmetric
    # device is the BAR path reflected: it burned 1896 s of real FDTD and
    # returned an exact copy of bar (imbalance 0.0000, phase -0.001° instead of
    # ~-91.6°). Correcting only the input_idx then hits a hard refusal, because
    # a plain x-mirror maps output onto input WITHOUT flipping y and is only
    # legal when both ports sit at the same y.
    cross_mirror = {}
cross_kwargs = {0: {"input_idx": 1}, OUT_SECTION: {"output_idx": 0}}
bar_kwargs = {0: {"input_idx": 1}, OUT_SECTION: {"output_idx": 1}}

hybrid_bar = HybridSolver(mirror_pairs={OUT_SECTION: 0}, **_shared)
hybrid_cross = HybridSolver(mirror_pairs=cross_mirror, **_shared)

# ---------------------------------------------------------------------------
# 5. Wavelength grid, and the per-wavelength sharding
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

# Every output name carries every parameter the result depends on, GEOMETRY
# FIRST. Eight separate results in this study were lost to filename collisions;
# the most recent merged a second design point over the first because the name
# carried only res/ct/nm.
_geo = (f"Ww{WAIST_WIDTH:g}_Lw{WAIST_LENGTH:g}"
        f"_Lt{TRANSITION_LENGTH:g}_L{MM_LENGTH:g}")
# Fabricated numbers must never be mistakable for measurements.
_DRY = "_DRYRUN" if DRYRUN else ""
# SYMMETRIC_MODE_REUSE changes the PHYSICS of the cross port, so two runs
# differing in it are different measurements and must not share a name.
_RU = "" if SYMMETRIC_MODE_REUSE else "_noreuse"
_TAG = (f"{_geo}_res{RESOLUTION:g}_nwl{N_WL}_ct{C_TRANS}_nm{N_MODES}{_RU}")

SHARD = os.environ.get("MMI2X2T_HY_WLIDX", "")
MERGED = os.environ.get("MMI2X2T_HY_MERGE", "") == "1"

if SHARD != "":
    _k = int(SHARD)
    if not 0 <= _k < N_WL:
        raise SystemExit(f"MMI2X2T_HY_WLIDX={_k} out of range for N_WL={N_WL}")
    wavelengths = wavelengths[_k:_k + 1]
    print(f"[shard] wavelength index {_k} of {N_WL}: {wavelengths[0]:.4f} µm")

print("=" * 88)
print("WAISTED 2x2 — HYBRID (FDTD ports + EME body) vs FULL-DEVICE FDTD")
print("=" * 88)
print(f"  Ww={WAIST_WIDTH} Lw={WAIST_LENGTH} Lt={TRANSITION_LENGTH} "
      f"L={MM_LENGTH} (stub s={stub:.3f}, body_stub={BODY_STUB})")
print(f"  res={RESOLUTION} decay_tol={DECAY_TOL:g} {N_WL} wl, "
      f"{N_MODES} modes, C_TRANS={C_TRANS}")
print(f"  body widths: face {mmi.width_at(-half_body + BODY_STUB):.3f}, "
      f"waist {mmi.width_at(0.0):.3f} µm")
print(f"  calibrate_body_neff={CALIBRATE_BODY} (REFUSED for a waisted body)")
print(f"  SYMMETRIC_MODE_REUSE={SYMMETRIC_MODE_REUSE}")
print("=" * 88, flush=True)

# ---------------------------------------------------------------------------
# 6. Geometry plot — shows which solver covers which section
# ---------------------------------------------------------------------------
if not MERGED:
    plot_hybrid_geometry(
        mmi, solver=hybrid_bar,
        in_modes=1, out_modes=1,
        title="Waisted 2×2 MMI — Hybrid Solver Sections",
        output_path=os.path.join(output_dir, "geometry_hybrid.png"),
        sub_region_dir=os.path.join(output_dir, "sections"),
        sub_region_show_pml=True,
        sub_region_show_eme_detail=True,
        sub_region_configs=[("in_top", {0: {"input_idx": 1}}, {OUT_SECTION: 0})],
    )
    print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")

# ---------------------------------------------------------------------------
# 7. Run the hybrid pass (or merge shards that already ran)
# ---------------------------------------------------------------------------
_shard_glob = os.path.join(output_dir, f"shard{_DRY}_wl*_{_TAG}.json")


def _shard_path(k):
    return os.path.join(output_dir, f"shard{_DRY}_wl{k}_{_TAG}.json")


if MERGED:
    shards = []
    for f in sorted(glob.glob(_shard_glob)):
        txt = open(f).read().lstrip()
        try:
            # raw_decode takes the FIRST complete document and ignores trailing
            # bytes, so a shard written by a pre-fix multi-rank race still loads.
            doc, end = json.JSONDecoder().raw_decode(txt)
        except (json.JSONDecodeError, ValueError) as exc:
            raise SystemExit(
                f"MERGE: {os.path.basename(f)} is not readable JSON ({exc}). A "
                "zero-byte or truncated shard is debris from a crashed run — "
                "delete it and re-run that wavelength.")
        if txt[end:].strip():
            print(f"  note: {os.path.basename(f)} had trailing bytes from a "
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
    print(f"MERGE: {len(shards)} wavelength shards, "
          f"{t_bar + t_cross:.0f}s of solve summed across them "
          "(they ran in parallel, so wall time was far less)")
elif DRYRUN:
    # MMI2X2T_HY_DRYRUN=1 fabricates the coefficients and skips the solve, so
    # the ENTIRE post-solve path (shard write -> merge -> comparison -> plots)
    # can be exercised in seconds. This exists because a bug in exactly that
    # path — treating HybridSolver.solve()'s SMatrix as a raw array — was only
    # reached after 10.5 h of solving and destroyed the job that hit it.
    print("DRYRUN: fabricating coefficients; NO physics in this run")
    rng = np.random.default_rng(0)
    c_bar = (rng.normal(size=len(wavelengths))
             + 1j * rng.normal(size=len(wavelengths))) * 0.5
    c_cross = (rng.normal(size=len(wavelengths))
               + 1j * rng.normal(size=len(wavelengths))) * 0.5
    t_bar = t_cross = 0.0
else:
    # ONE shared cache across both excitations: section 0's geometry_kwargs are
    # identical between the bar and cross calls (both excite the top arm), so
    # it is simulated exactly once, and the body likewise.
    section_cache = {}
    t0 = time.time()
    sm_bar = hybrid_bar.solve(mmi, wavelengths, geometry_kwargs=bar_kwargs,
                              cache=section_cache, save_dir=save_dir)
    t_bar = time.time() - t0
    print(f"[hybrid] bar solved in {t_bar:.0f}s", flush=True)
    t0 = time.time()
    sm_cross = hybrid_cross.solve(mmi, wavelengths, geometry_kwargs=cross_kwargs,
                                  cache=section_cache, save_dir=save_dir)
    t_cross = time.time() - t0
    print(f"[hybrid] cross solved in {t_cross:.0f}s", flush=True)
    # ⚠ HybridSolver.solve() returns an SMatrix, NOT a raw array. np.asarray()
    # on it yields a 0-d OBJECT array, and the resulting
    # "AttributeError: 'SMatrix' object has no attribute 'real'" surfaced only
    # after 10.5 hours of solving. Extract the coefficients properly.
    c_bar = np.asarray(sm_bar.coeff(port="out", in_mode=1, out_mode=1),
                       dtype=complex).reshape(len(wavelengths))
    c_cross = np.asarray(sm_cross.coeff(port="out", in_mode=1, out_mode=1),
                         dtype=complex).reshape(len(wavelengths))

# --- a shard writes its one wavelength and stops here ----------------------
if SHARD != "" and not MERGED:
    # ⚠ ONLY ONE RANK MAY WRITE. Every MPI rank reaches this line; when all of
    # them wrote the same path their writes interleaved and produced
    # concatenated JSON ("Extra data: line 30").
    try:
        import meep as _mp
        _is_master = _mp.am_master()
    except Exception:
        _is_master = True
    if not _is_master:
        raise SystemExit(0)
    sp = _shard_path(int(SHARD))
    # ATOMIC: write to a temp file and rename. A crash between open() and
    # dump() otherwise leaves a ZERO-BYTE file that looks like a valid shard.
    tmp = sp + ".partial"
    with open(tmp, "w") as fh:
        json.dump({"wl_index": int(SHARD),
                   "wavelength": float(wavelengths[0]),
                   "n_wl_total": N_WL,
                   "config": {"res": RESOLUTION, "c_trans": C_TRANS,
                              "n_modes": N_MODES, "decay_tol": DECAY_TOL},
                   "design": {"Ww": WAIST_WIDTH, "Lw": WAIST_LENGTH,
                              "Lt": TRANSITION_LENGTH, "L": MM_LENGTH},
                   "seconds": {"bar": t_bar, "cross": t_cross},
                   "c_bar": [c_bar[0].real, c_bar[0].imag],
                   "c_cross": [c_cross[0].real, c_cross[0].imag]}, fh, indent=2)
    os.replace(tmp, sp)
    print(f"[shard] written → {os.path.basename(sp)}")
    h, m = _MEMO_STATS["hits"], _MEMO_STATS["misses"]
    if h + m:
        print(f"[shard] mode memo: {m} solved, {h} reused "
              f"({100*h/(h+m):.0f}% avoided)")
    print("[shard] merge with MMI2X2T_HY_MERGE=1 once every shard has landed.")
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# 8. Scores, and the comparison against full-device FDTD
# ---------------------------------------------------------------------------
# From the pre-registered criterion (see mmi_2x2_fdtd.py in this folder).
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
# The reference ran 51 wl over the same band. Rather than hard-coding the 51→11
# stride, match each hybrid wavelength to its nearest reference wavelength and
# REFUSE if any pair is not essentially identical — comparing two engines at
# wavelengths that merely nearly agree would fold a wavelength offset into the
# engine disagreement being measured.
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
print(f"[compare] matched {len(idx)} wavelengths to the stored reference "
      f"(worst |dwl| = {mismatch:.1e} µm)")

dTb = np.abs(mh["T_bar"] - mf["T_bar"])
dTc = np.abs(mh["T_cross"] - mf["T_cross"])

print("\n" + "=" * 88)
print("HYBRID vs FULL-DEVICE FDTD — identical geometry, identical wavelengths")
print("=" * 88)
print(f"{'port':>6} {'max|dT|':>10} {'mean|dT|':>10} {'% of mean T_fdtd':>17} "
      f"{'worst wl':>10}")
for name, d, ref in (("bar", dTb, mf["T_bar"]), ("cross", dTc, mf["T_cross"])):
    i = int(np.argmax(d))
    print(f"{name:>6} {d.max():10.4f} {d.mean():10.4f} "
          f"{100*d.mean()/ref.mean():17.2f} {wavelengths[i]:10.4f}")
print(f"\n{'quantity':>18} {'hybrid':>10} {'full FDTD':>10} {'delta':>10}")
for q in ("excess_loss_db", "imbalance", "phase_deg", "phase_max_dev_deg",
          "phase_spread_deg"):
    print(f"{q:>18} {mh[q]:10.4f} {mf[q]:10.4f} {mh[q]-mf[q]:+10.4f}")
print("\nStraight 2x2 for scale (CALIBRATED body): bar 0.62% / cross 0.59% of")
print("mean T. This device's hybrid CANNOT calibrate, so a larger disagreement")
print("here is the expected outcome, and its size is the measurement. A")
print("cross-port deficit of ~1-2% persists and is the known open item: it is")
print("NOT the parity shortcut, NOT the mode count, NOT C_TRANS convergence,")
print("and only marginally beta. Do not re-test those four.")

# ---------------------------------------------------------------------------
# 9. Save and plot
# ---------------------------------------------------------------------------
smat = SMatrix(wavelengths=wavelengths)
smat.set(c_bar, port="bar")
smat.set(c_cross, port="cross")
smat.check_passive()
smat.save(os.path.join(output_dir, f"smatrix_hybrid{_DRY}_{_TAG}.pkl"))
print(f"\nSMatrix saved → smatrix_hybrid{_DRY}_{_TAG}.pkl")

out_json = os.path.join(output_dir, f"hybrid_vs_fdtd{_DRY}_{_TAG}.json")
with open(out_json, "w") as fh:
    json.dump({"config": {"res": RESOLUTION, "n_wl": N_WL, "c_trans": C_TRANS,
                          "decay_tol": DECAY_TOL, "n_modes": N_MODES,
                          "symmetric_mode_reuse": SYMMETRIC_MODE_REUSE,
                          "calibrate_body_neff": CALIBRATE_BODY},
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
print(f"Comparison saved → {os.path.basename(out_json)}")

fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))

ax = axes[0]
ax.plot(wavelengths, mh["T_bar"], color="steelblue", marker="o", ms=3,
        label="Bar — hybrid")
ax.plot(wavelengths, mf["T_bar"], color="steelblue", linestyle="--",
        label="Bar — full FDTD")
ax.plot(wavelengths, mh["T_cross"], color="darkorange", marker="o", ms=3,
        label="Cross — hybrid")
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
ax.plot(wavelengths, mh["total"], color="k", marker="o", ms=3, label="Hybrid")
ax.plot(wavelengths, mf["total"], color="k", linestyle="--", label="Full FDTD")
ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Total transmission")
ax.set_title("Total power")
ax.legend(fontsize=8)
ax.grid(True)

# ⚠ |T| DISCARDS WHAT THIS DEVICE IS FOR — the bar/cross quadrature is the
# function, and an intensity plot cannot see it.
ax = axes[2]
ax.plot(wavelengths, mh["phase_per_wl"], color="seagreen", marker="o", ms=3,
        label="Hybrid")
ax.plot(wavelengths, mf["phase_per_wl"], color="seagreen", linestyle="--",
        label="Full FDTD")
ax.axhline(CONTROL_PHASE_DEG, color="gray", linestyle=":", lw=0.8,
           label=f"Straight control ({CONTROL_PHASE_DEG:.2f}°)")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("arg(cross / bar) (deg)")
ax.set_title("Quadrature")
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle(f"Waisted 2×2 MMI — hybrid FDTD/EME vs full-device FDTD "
             f"({_geo}, res {RESOLUTION:g}, {N_MODES} modes, C_TRANS {C_TRANS})",
             fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, f"hybrid_results{_DRY}_{_TAG}.png"), dpi=300)
print(f"Results plot saved → hybrid_results{_DRY}_{_TAG}.png")
