"""Waisted 2x2: HYBRID (FDTD ports + EME body) vs the FULL-DEVICE FDTD.

validation_markdown/mmi_2x2_waist_body_validation.md -- step 8 check 4's real
question for this geometry: does the FDTD/EME cascade reproduce a single-domain
brute-force run of the SAME waisted device? The straight 2x2 has that
comparison (gap 5: 0.76%/0.59% max|dT|, job 11585849). The waisted one does
not, because step 9 (wiring) was gated on the physics settling and never ran.

WHAT ALREADY EXISTS, and what this adds:
    straight 2x2 : EME + hybrid + full FDTD, hybrid-vs-FDTD measured
    waisted 2x2  : EME + full FDTD (jobs 12736252/12743476) -- NO hybrid
This script supplies the missing hybrid leg and compares it to the stored
full-device SMatrix, on identical geometry and the same 11-wavelength grid.

THREE THINGS THE STRAIGHT DEVICE'S HYBRID DOES THAT THIS CANNOT:

1. `calibrate_body_neff=True` is HARD-REFUSED for a waisted body --
   `assert_uniform_cross_section()` fires from `HybridSolver._calibrated_neff()`
   (hybrid.py:1532): "a calibrated body propagates on one per-mode beta, which
   only describes a uniform cross-section". That flag is what closed the
   engine-agreement gap on every other device here (~2.3%/mode between MPB's
   index and meep's stepped phase; on the coupler it moved L50 from 5.37% off
   to 0.11%). So this cascade carries an UNCORRECTED body beta, and that is
   precisely what the comparison measures.

2. `solve_length_sweep(sweep_section=1)` needs a length-independent uniform
   slab. The waisted body is neither, so this is a single-length solve.

3. The body goes through `solve_section()`'s real MULTI-CELL branch -- one
   cross-section per cell, real interface overlaps, non-diagonal transfer
   matrix. `num_cells` therefore matters here, unlike every uniform-body
   script in this repo which uses num_cells=[1].

CELL COUNT. At the design point the EME body section spans exactly
[ramp_in, waist, ramp_out] = 6 + 8 + 6 um: body_stub=0.25 equals the design's
own stub s=0.25, so the FDTD sections fold away the stubs exactly. Cells are
distributed UNIFORMLY over the section, so N cells put 0.3*N in each 6 um
ramp; step 5a converged the ramps at C_TRANS=8, needing N >= 27. Default 40
(0.5 um/cell: 12 per ramp, 16 in the waist) clears that with margin, and
BODY_CELLS is swept-able to show the answer does not depend on it.

MODE_SIGNS is inherited from the straight device. That is justified rather
than assumed: the input FDTD section's polygons are byte-identical between
straight and waisted (measured 2026-08-31, step 9 -- `section_gds()`
re-origins to x_start=0), and the vector is a property of that section's
mirror symmetry. SYMMETRIC_MODE_REUSE=0 re-derives it by simulating both arms
if you want the assumption removed at the cost of one more FDTD solve.

Run:
    MMI2X2_HY_RES=15 MMI2X2_HY_NWL=3 MMI2X2_HY_CELLS=12 python .../this.py
    python .../this.py                      # real config
"""
import json
import os
import time

import numpy as np

from beef import (MMI, EMESolver, FDTDSolver, HybridSolver, MirrorSpec,
                  SMatrix, SOI_220)

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT_ROOT = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "verification"))
_EXAMPLE_OUT_ROOT = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output"))
output_dir = os.path.join(_EXAMPLE_OUT_ROOT, "mmi_2x2_tapered")
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_sim_cache"))
FDTD_REF_DIR = os.path.join(_EXAMPLE_OUT_ROOT, "mmi_2x2_tapered")
CRITERION_PATH = os.path.join(_OUT_ROOT, "mmi_2x2_waist_metric",
                              "criterion.json")


def _env_f(n, d):
    return float(os.environ.get(n, d))


def _env_i(n, d):
    return int(os.environ.get(n, d))


RES = _env_f("MMI2X2_HY_RES", 40.0)
WL_MIN, WL_MAX = 1.530, 1.565
N_WL = _env_i("MMI2X2_HY_NWL", 11)      # the criterion's own grid
BODY_CELLS = _env_i("MMI2X2_HY_CELLS", 40)
DECAY_TOL = _env_f("MMI2X2_HY_DECAYTOL", 1e-11)
N_MODES = _env_i("MMI2X2_HY_NMODES", 14)
# Modes carried at the two FDTD<->EME SEAMS specifically. Defaults to N_MODES
# so nothing changes unless asked; see the interface_modes comment below for
# why the seams want a different (larger) count than the waist.
SEAM_MODES = _env_i("MMI2X2_HY_SEAMMODES", N_MODES)
SYMMETRIC_MODE_REUSE = _env_i("MMI2X2_HY_REUSE", 1) == 1
BODY_STUB = 0.25
DPML = _env_f("MMI2X2_HY_DPML", 0.5)
C_TRANS = _env_i("MMI2X2_HY_CTRANS", 8)   # step 5a floor; cells per ramp
DRYRUN = _env_i("MMI2X2_HY_DRYRUN", 0) == 1
# Fabricated numbers must never be mistakable for measurements: every file a
# dry run writes carries DRYRUN in its name.
_DRY = "_DRYRUN" if DRYRUN else ""
# SYMMETRIC_MODE_REUSE changes the PHYSICS of the cross port -- with it on, the
# cross result is derived from bar via MODE_SIGNS parity (measured: bar 16353 s,
# cross 0.03 s, i.e. not simulated at all); with it off, both arms are actually
# solved. Two runs that differ in it are therefore different measurements and
# must not share a shard filename. EMPTY when reuse is on, so every shard
# already on disk keeps its name and stays readable.
_RU = "" if SYMMETRIC_MODE_REUSE else "_noreuse"
# Empty when the seam count equals N_MODES, so every existing file keeps
# its name; only a genuine per-interface run gets a new one.
_SM = "" if SEAM_MODES == N_MODES else f"_sm{SEAM_MODES}"

# design point (step 6, as scored at N_TAPER=7 -- see step 5b: its phase PASS
# is an artifact of that value. The geometry is still the right one to compare
# the two ENGINES on; this script is a method check, not a design claim.)
WW = _env_f("MMI2X2_HY_WW", 2.24)
LW = _env_f("MMI2X2_HY_LW", 6.0)
LT = _env_f("MMI2X2_HY_LT", 6.0)
L_WAIST = _env_f("MMI2X2_HY_L", 21.5)
# Default is the 30.50 um point, NOT the 29.50 um one step 6 ranked first.
# Measured 2026-09-03: the 29.50 um candidates have uniform stub s=0.25, and a
# hybrid FDTD section pads its drawn waveguide by `dpml` (0.5 um) past the
# interior seam (fdtd.py:3239, so the PML never sees an abrupt core->cladding
# step). That padding lands 0.5 um inside the waist ramp, and MMI._check_cut()
# correctly refuses: "a section boundary cannot sit on a single cross-section".
# So a hybrid-buildable waisted device needs
#     s >= body_stub + dpml = 0.25 + 0.5 = 0.75 um
# which is 3x step 9's s >= 0.25 (that one only kept the sections IDENTICAL).
# 52 of step 6's 56 passers clear it; both 29.50 um candidates do not.

# Canonical 14-entry vector, derived on the straight device's input section.
_MODE_SIGNS_14 = np.array([1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, 1])
# Sized to SEAM_MODES, not N_MODES: the vector is applied to the mirrored
# OUTPUT FDTD section, whose interface is the seam.
if SEAM_MODES <= len(_MODE_SIGNS_14):
    # TRUNCATION IS VALID, and this is why. MODE_SIGNS is a per-mode PARITY
    # vector for the input FDTD section's mirror symmetry -- entry k is a
    # property of mode k alone. Asking the solver for 10 modes instead of 14
    # returns the SAME lowest 10 modes in the same order (the guard in
    # EMESolver.compute_modes only ever reaches FURTHER up the spectrum, never
    # reorders what it already accepted), so entries 0..N-1 still describe
    # them. Parity does not depend on how many modes were requested.
    #
    # Needed because the waist supports only 10 distinct modes at this
    # wavelength (validation_markdown FOLLOW-UP 2026-09-09: N=11 is an exact
    # duplicate, overlap 1.000000), so N_MODES=14 was asking for 4 modes that
    # do not exist and ~50-58% of ramp cells were silently re-solving.
    #
    # If you distrust the truncation, MMI2X2_HY_REUSE=0 re-derives the signs
    # by simulating both arms -- one extra FDTD section solve per wavelength.
    MODE_SIGNS = _MODE_SIGNS_14[:SEAM_MODES]
else:
    raise SystemExit(
        f"SEAM_MODES={SEAM_MODES} exceeds the {len(_MODE_SIGNS_14)} known MODE_SIGNS "
        "entries; re-derive with MMI2X2_HY_REUSE=0 (parity for modes beyond "
        "13 was never measured, and guessing it would silently corrupt the "
        "second arm's coefficients)")
print(f"MODE_SIGNS: using {len(MODE_SIGNS)} of {len(_MODE_SIGNS_14)} entries "
      f"{list(MODE_SIGNS)}")

wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)

# --- wavelength sharding -----------------------------------------------------
# MPB does NOT scale with MPI ranks -- meep solves the eigenproblem REDUNDANTLY
# on every rank (measured in this repo: 8 ranks 3110.7 s vs 24 ranks 3055.2 s,
# 3x the cores for 1.8%). So a mode-solve-dominated job gains nothing from more
# ranks; the 48-rank hybrid attempts were doing 48x redundant work and burning
# ~576 core-hours to do ~12 core-hours of useful solving.
#
# Wavelengths ARE independent here, so they parallelise at the JOB level: run
# one shard per wavelength with few ranks each, then merge. (The solve_section
# cache key contains the whole wavelength array, so this cannot be done through
# the cache -- it has to be separate runs.)
SHARD = os.environ.get("MMI2X2_HY_WLIDX", "")
if SHARD != "":
    _k = int(SHARD)
    if not 0 <= _k < N_WL:
        raise SystemExit(f"MMI2X2_HY_WLIDX={_k} out of range for N_WL={N_WL}")
    wavelengths = wavelengths[_k:_k + 1]
    print(f"[shard] wavelength index {_k} of {N_WL}: {wavelengths[0]:.4f} um")

MERGED = os.environ.get("MMI2X2_HY_MERGE", "") == "1"
if MERGED:
    import glob
    pat = os.path.join(output_dir,
                       f"shard{_DRY}_wl*_Ww{WW:g}_Lw{LW:g}_Lt{LT:g}_L{L_WAIST:g}"
                       f"_res{RES:g}_nwl{N_WL}_ct{C_TRANS}"
                       f"_nm{N_MODES}{_SM}{_RU}.json")
    shards = []
    for f in sorted(glob.glob(pat)):
        try:
            _txt = open(f).read().lstrip()
            # raw_decode takes the FIRST complete document and ignores any
            # trailing bytes, so a shard written by the pre-fix multi-rank race
            # still loads. The leading document is complete and correct --
            # verified by hand on wl0/wl7 of job 12892459.
            _doc, _end = json.JSONDecoder().raw_decode(_txt)
            if _txt[_end:].strip():
                print(f"  note: {os.path.basename(f)} had "
                      f"{len(_txt[_end:].strip())} trailing bytes from the "
                      "multi-rank write race; using the leading document")
            shards.append(_doc)
        except (json.JSONDecodeError, ValueError) as exc:
            raise SystemExit(
                f"MERGE: {os.path.basename(f)} is not readable JSON ({exc}). "
                "A zero-byte or truncated shard is debris from a crashed run "
                "-- delete it and re-run that wavelength. Shards are written "
                "atomically now, so this should not recur.")
    have = sorted(sh["wl_index"] for sh in shards)
    if have != list(range(N_WL)):
        raise SystemExit(
            f"MERGE: need wavelength shards {list(range(N_WL))}, have {have}. "
            "Run the missing shards first -- a partial merge would silently "
            "compare on fewer wavelengths than the reference.")
    shards.sort(key=lambda sh: sh["wl_index"])
    wls = [sh["wavelength"] for sh in shards]
    cb = np.array([complex(*sh["c_bar"]) for sh in shards])
    cc = np.array([complex(*sh["c_cross"]) for sh in shards])
    t_bar = sum(sh["seconds"]["bar"] for sh in shards)
    t_cross = sum(sh["seconds"]["cross"] for sh in shards)
    wl_txt = ", ".join(f"{w:.4f}" for w in wls)
    print(f"MERGE: {len(shards)} wavelength shards at [{wl_txt}] um")
    print(f"MERGE: {t_bar + t_cross:.0f}s of solve summed across shards "
          "(they ran in parallel, so wall was far less)")

mmi = MMI(n_inputs=2, n_outputs=2, mm_length=L_WAIST, mm_width=2.8,
          waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
          waveguide_offset=1.5, waveguide_height=0.22,
          body_waist_width=WW, body_waist_length=LW,
          body_transition_length=LT)

half_body = mmi.mm_length / 2
half_total = mmi.total_length / 2

# Fail here, with the reason, rather than 200 frames deep inside meep.
_s_uniform = (L_WAIST - 2 * LT - LW) / 2
_dpml = DPML
_s_min = BODY_STUB + _dpml
if _s_uniform < _s_min - 1e-9:
    raise SystemExit(
        f"REFUSING: uniform stub s={_s_uniform:.3f} um < {_s_min:.3f} um.\n"
        f"  A hybrid FDTD section pads its drawn waveguide by dpml={_dpml} um\n"
        f"  past the interior seam (fdtd.py, so the PML never sees an abrupt\n"
        f"  core->cladding step). With body_stub={BODY_STUB} that padding\n"
        f"  reaches x = +/-{BODY_STUB + _dpml:.2f} um into the body from the\n"
        f"  seam; if the waist ramp starts before that, MMI._check_cut()\n"
        f"  refuses -- a section boundary cannot sit where the width varies.\n"
        f"  Need s = (L - 2*Lt - Lw)/2 >= body_stub + dpml = {_s_min:.2f} um.\n"
        f"  This is 3x step 9's s >= {BODY_STUB} (which only kept the FDTD\n"
        f"  sections IDENTICAL, not constructible).")

print("=" * 88)
print("WAISTED 2x2 -- HYBRID vs FULL-DEVICE FDTD")
print("=" * 88)
print(f"  Ww={WW} Lw={LW} Lt={LT} mm_length={L_WAIST} "
      f"(stub s={(L_WAIST - 2*LT - LW)/2:.3f}, body_stub={BODY_STUB})")
print(f"  res={RES} decay_tol={DECAY_TOL:g} {N_WL} wl, {N_MODES} body modes"
      + (f", {SEAM_MODES} at the two FDTD seams" if SEAM_MODES != N_MODES
         else " (uniform across all interfaces)"))
print(f"  body widths: face {mmi.width_at(-half_body + BODY_STUB):.3f}, "
      f"mid {mmi.width_at(0.0):.3f} um")
print(f"  calibrate_body_neff=False (REFUSED for a waisted body)")
print(f"  SYMMETRIC_MODE_REUSE={SYMMETRIC_MODE_REUSE}")
print("=" * 88, flush=True)

# ---------------------------------------------------------------------------
# CONTENT-ADDRESSED MODE MEMO  (user's observation, 2026-09-06)
#
# The waisted body is MIRROR-SYMMETRIC about x=0 and uniform through the waist,
# so most of its cells present a cross-section that has already been solved:
#
#     ramp in  (-10 -> -4 um)   12 cells   12 distinct
#     waist    ( -4 -> +4 um)   16 cells    1 distinct  (uniform at Ww)
#     ramp out ( +4 -> +10 um)  12 cells    0 distinct  (exact mirrors of ramp in)
#     ------------------------------------------------------------------
#     total                     40 cells   13 distinct   -> ~3x fewer solves
#
# Cell i sits at x and cell 39-i at -x, and inside the body the cross-section
# is just a slab of width w(x) centred on y=0 with w(x) = w(-x) -- so the two
# are IDENTICAL, not merely similar.
#
# ⚠ WHY THIS IS KEYED ON CONTENT AND NOT ON WIDTH OR ON THE SYMMETRY ITSELF.
# This repo has already been burned by the tempting version of this idea:
# width-keyed mode reuse was measured BREAKING PHYSICS on the 1x4 (totalT
# 0.9703 -> 0.1507) because equal width does NOT mean equal cross-section once
# a cell sits off-axis. So the key here is a hash of the actual refractive
# index arrays (nx, ny, nz) plus wavelength, temperature and mode count. Then
# the reuse is exact by construction: identical eigenproblem in, identical
# modes out. If the symmetry argument above were wrong, the hashes simply would
# not match and every cell would be solved for real -- the failure mode is
# "no speedup", never "wrong answer".
#
# Set MMI2X2_HY_MEMO=0 to disable (used by the A/B validation).
# ---------------------------------------------------------------------------
MEMO_ON = _env_i("MMI2X2_HY_MEMO", 1) == 1
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

_fdtd_common = dict(res=RES, sim_width=5.0, sim_height=1.5, dpml=DPML,
                    source_margin=0.0, monitor_margin=0.25, three_d=True,
                    material_stack=SOI_220, decay_tol=DECAY_TOL)

input_solver = FDTDSolver(half_domain_aperture=True, **_fdtd_common)
output_solver = FDTDSolver(half_domain_aperture=True, **_fdtd_common)


def _eme(n_cells, length):
    # num_cells IS read by solve_section() (unlike every uniform-body script
    # in this repo, which passes [1] because its body is a single slab)
    return EMESolver(num_modes=[N_MODES], num_cells=[n_cells],
                     lengths=[length], oxide_width=5.0, oxide_height=1.5,
                     res=RES, material_stack=SOI_220, mode_solver="meep")


# --- BODY AS THREE REGIONS, NOT ONE UNIFORMLY-DIVIDED SECTION ---------------
# solve_section() divides a section UNIFORMLY into num_cells[0] cells; it has
# no notion of internal region boundaries. Handing it the whole body therefore
# forced one cell size everywhere: 40 cells of 0.5 um, chosen so each 6 um ramp
# got >= 8 cells (step 5a's converged C_TRANS) AND boundaries landed exactly on
# the ramp/waist junctions at +/-4 um. The cost was 16 cells inside the WAIST,
# which is a uniform slab and needs exactly ONE -- propagation through it is
# analytic, and the extra cells only insert identity interfaces between
# identical cross-sections.
#
# HybridSolver is "a generic N-section cascade: any solver type, in any order,
# any count >= 2", so the body can simply BE three sections with their own cell
# counts -- the same decomposition the step-6 free sweep already uses and
# validated for this geometry:
#
#     ramp in   6 um   C_TRANS cells   (needs the staircase)
#     waist     Lw um  1 cell          (uniform -- one cell is exact)
#     ramp out  6 um   C_TRANS cells   (mirror of ramp in; memo makes it free)
#
# 8 + 1 + 8 = 17 cells instead of 40, and with the content memo only 9 are
# distinct mode solves.
# The body's real structure, mirroring the step-6 free sweep's 9-region layout:
#     stub (uniform mm_width) | ramp | waist (uniform Ww) | ramp | stub
# Each UNIFORM region gets ONE cell because one cell is exact for a uniform
# region -- propagation is analytic and extra cells only insert identity
# interfaces. Only the ramps need the C_TRANS staircase.
#
# My first cut of this lumped the stubs into the ramp sections, which made the
# "ramps" 7.5 um instead of 6.0 and left a uniform run inside a region being
# staircased. Separated here.
s_stub_len = (L_WAIST - 2 * LT - LW) / 2          # 1.75 um at the design point
ramp_x0 = -half_body + s_stub_len                  # ramp starts after the stub
waist_x0, waist_x1 = -LW / 2, +LW / 2
boundaries = [-half_total,
              -half_body + BODY_STUB,              # FDTD -> body
              ramp_x0,                             # stub -> ramp
              waist_x0,                            # ramp -> waist
              waist_x1,                            # waist -> ramp
              -ramp_x0,                            # ramp -> stub
              +half_body - BODY_STUB,              # body -> FDTD
              +half_total]
stub_in_len = ramp_x0 - (-half_body + BODY_STUB)   # stub minus the folded part

solvers = [input_solver,
           _eme(1, stub_in_len),         # stub  -- uniform at mm_width
           _eme(C_TRANS, LT),            # ramp  -- needs the staircase
           _eme(1, LW),                  # waist -- uniform at Ww, ONE cell
           _eme(C_TRANS, LT),            # ramp  -- mirror of the first
           _eme(1, stub_in_len),         # stub  -- uniform, mirror
           output_solver]
N_SECTIONS = len(solvers)
OUT_SECTION = N_SECTIONS - 1
n_cells_total = 1 + C_TRANS + 1 + C_TRANS + 1
print(f"  body regions: stub {stub_in_len:.2f}um(1) | ramp {LT:g}um({C_TRANS}) "
      f"| waist {LW:g}um(1) | ramp {LT:g}um({C_TRANS}) | stub "
      f"{stub_in_len:.2f}um(1)  = {n_cells_total} cells, "
      f"{C_TRANS + 2} distinct after the memo")

_shared = dict(solvers=solvers,
               boundaries=boundaries,
               # PER-INTERFACE mode counts, not one number for all six.
               #
               # The six section-to-section interfaces sit at DIFFERENT widths:
               #   0  FDTD_in | stub    2.80 um (face)
               #   1  stub    | ramp    2.80 um
               #   2  ramp    | waist   2.24 um (waist)
               #   3  waist   | ramp    2.24 um
               #   4  ramp    | stub    2.80 um
               #   5  stub    | FDTD_out 2.80 um
               # The 2.24 um waist supports exactly 10 distinct modes
               # (N=11 is an exact duplicate), but the 2.80 um face supports
               # ~14. A single N_MODES is therefore forced to be simultaneously
               # too LARGE in the waist and too SMALL at the seams -- and
               # dropping 14->10 to fix the waist necessarily STARVED the
               # FDTD<->EME seam projection.
               #
               # That matters because bar and cross share the same EME body
               # cascade and differ only in the output port, so the measured
               # bar/cross asymmetry (hybrid +0.0436 against FDTD's +0.0019)
               # cannot originate in the body -- it has to live at the seam.
               # A truncated seam basis over-coupling the symmetric (bar)
               # combination is exactly that signature.
               #
               # SEAM_MODES defaults to N_MODES, so the default behaviour is
               # unchanged and every earlier result stays reproducible.
               # Per-cell counts stay rank-clean at SEAM_MODES=14/N_MODES=10:
               # the stub is uniform 2.80 um and carries max(14, 10) = 14
               # comfortably, while ramp and waist see 10 on both sides.
               interface_modes=([SEAM_MODES] + [N_MODES] * (N_SECTIONS - 3)
                                + [SEAM_MODES]),
               calibrate_body_neff=False)

if SYMMETRIC_MODE_REUSE:
    # NOTE the section index: with the body split into ramp/waist/ramp the
    # chain is [fdtd_in, ramp, waist, ramp, fdtd_out], so the OUTPUT section is
    # index 4, not 2. Getting this wrong would mirror a body region onto the
    # input section instead of the output one.
    cross_mirror = {OUT_SECTION: MirrorSpec(section=0, partner_kwargs={"input_idx": 0},
                                  mode_signs=MODE_SIGNS)}
    cross_kwargs = {0: {"input_idx": 1}, OUT_SECTION: {"output_idx": 0}}
else:
    # FIXED 2026-09-10. This branch was BOTH wrong and unable to express the
    # cross path, and the two faults hid each other:
    #
    #   it had  cross_mirror = {OUT_SECTION: 0}
    #           cross_kwargs = {0: {"input_idx": 0}, OUT_SECTION: {"output_idx": 0}}
    #
    # input_idx=0 -> output_idx=0 excites the BOTTOM arm and reads the BOTTOM
    # output, which on this y-mirror-symmetric device is the BAR path
    # reflected. Measured: it burned 1896 s of real FDTD and returned an exact
    # copy of bar (T_bar == T_cross to 4 dp, imbalance 0.0000, phase -0.001
    # deg instead of ~-91.6). It had presumably never been run on this device.
    #
    # Correcting only the input_idx then hits a hard refusal from the library,
    # which is the deeper point:
    #   ValueError: mirror_pairs[6]=0: center_at mismatch at x=10.85 vs
    #   mirrored x=-10.85: -0.75 vs 0.75 um
    # A plain x-mirror maps the output section onto the input section WITHOUT
    # flipping y, so it is only legal when both ports sit at the same y. The
    # cross path (top in, bottom out) needs a y-flip, which is exactly what the
    # reuse branch encodes via MirrorSpec(mode_signs=MODE_SIGNS). So the
    # original input_idx=0 was not a typo -- it was the only way to keep the
    # plain mirror legal, and it silently changed which quantity was computed.
    #
    # The honest "simulate both arms" configuration is therefore NO mirror at
    # all: solve the output section directly on its own geometry. That costs
    # one extra FDTD section solve, which is the whole point of the test.
    cross_mirror = {}
    cross_kwargs = {0: {"input_idx": 1}, OUT_SECTION: {"output_idx": 0}}

hybrid_bar = HybridSolver(mirror_pairs={OUT_SECTION: 0}, **_shared)
hybrid_cross = HybridSolver(mirror_pairs=cross_mirror, **_shared)
bar_kwargs = {0: {"input_idx": 1}, OUT_SECTION: {"output_idx": 1}}

# MMI2X2_HY_DRYRUN=1 fabricates the coefficients and skips the solve, so the
# ENTIRE post-solve path (shard write -> merge -> comparison -> JSON) can be
# exercised in seconds. This exists because a bug in exactly that path --
# treating HybridSolver.solve()'s SMatrix as a raw array -- was only reached
# after 10.5 h of solving and destroyed job 12877872_0. Post-solve code must be
# testable without paying for the solve.
if DRYRUN and not MERGED:
    print("DRYRUN: fabricating coefficients; NO physics in this run")
    rng = np.random.default_rng(0)
    cb = (rng.normal(size=len(wavelengths))
          + 1j * rng.normal(size=len(wavelengths))) * 0.5
    cc = (rng.normal(size=len(wavelengths))
          + 1j * rng.normal(size=len(wavelengths))) * 0.5
    t_bar = t_cross = 0.0
elif MERGED:
    print("MERGE: skipping the solve; coefficients came from the shards")
else:
    section_cache = {}
    t0 = time.time()
    cb = hybrid_bar.solve(mmi, wavelengths, geometry_kwargs=bar_kwargs,
                          cache=section_cache, save_dir=save_dir)
    t_bar = time.time() - t0
    print(f"[hybrid] bar solved in {t_bar:.0f}s", flush=True)
    t0 = time.time()
    cc = hybrid_cross.solve(mmi, wavelengths, geometry_kwargs=cross_kwargs,
                            cache=section_cache, save_dir=save_dir)
    t_cross = time.time() - t0
    print(f"[hybrid] cross solved in {t_cross:.0f}s", flush=True)
    # HybridSolver.solve() returns an SMatrix, NOT a raw array -- its docstring
    # says "SMatrix with coeff(port=port, in_mode=i, out_mode=o)". np.asarray()
    # on it yields a 0-d OBJECT array, so `cb[0].real` later raised
    # "AttributeError: 'SMatrix' object has no attribute 'real'" -- after 10.5
    # hours of solving (job 12877872_0). Extract the coefficients properly.
    cb = np.asarray(cb.coeff(port="out", in_mode=1, out_mode=1),
                    dtype=complex).reshape(len(wavelengths))
    cc = np.asarray(cc.coeff(port="out", in_mode=1, out_mode=1),
                    dtype=complex).reshape(len(wavelengths))

if SHARD != "":
    sp = os.path.join(output_dir,
                      f"shard{_DRY}_wl{int(SHARD)}_Ww{WW:g}_Lw{LW:g}_Lt{LT:g}"
                      f"_L{L_WAIST:g}_res{RES:g}_nwl{N_WL}_ct{C_TRANS}"
                      f"_nm{N_MODES}{_SM}{_RU}.json")
    # ATOMIC: write to a temp file and rename. A crash between open() and dump()
    # otherwise leaves a ZERO-BYTE file that looks like a shard -- exactly what
    # job 12877872_0 left behind when it died on the SMatrix bug, and what a
    # later merge would have choked on (or, worse, counted).
    # ⚠ ONLY ONE RANK MAY WRITE. Every MPI rank reaches this line, and all 8
    # were writing the SAME path: their writes interleaved and produced
    # concatenated JSON ("Extra data: line 30"). The earlier non-atomic version
    # survived only because each rank wrote identical complete content and the
    # last writer won -- adding a temp-file rename turned a benign race into a
    # corrupting one, plus FileNotFoundError for the losers of the rename.
    try:
        import meep as _mp
        _is_master = _mp.am_master()
    except Exception:
        _is_master = True
    if not _is_master:
        raise SystemExit(0)
    _tmp = sp + f".partial"
    with open(_tmp, "w") as fh:
        json.dump({"wl_index": int(SHARD), "wavelength": float(wavelengths[0]),
                   "n_wl_total": N_WL,
                   "config": {"res": RES, "c_trans": C_TRANS,
                              "n_modes": N_MODES, "decay_tol": DECAY_TOL},
                   "design": {"Ww": WW, "Lw": LW, "Lt": LT, "L": L_WAIST},
                   "seconds": {"bar": t_bar, "cross": t_cross},
                   "c_bar": [cb[0].real, cb[0].imag],
                   "c_cross": [cc[0].real, cc[0].imag]}, fh, indent=2)
    os.replace(_tmp, sp)
    print(f"[shard] written -> {os.path.basename(sp)}")
    h, m = _MEMO_STATS["hits"], _MEMO_STATS["misses"]
    if h + m:
        print(f"[shard] mode memo: {m} solved, {h} reused "
              f"({100*h/(h+m):.0f}% avoided)")
    print("[shard] merge with MMI2X2_HY_MERGE=1 once every shard has landed.")
    raise SystemExit(0)

with open(CRITERION_PATH) as fh:
    CRIT = json.load(fh)
THR = CRIT["thresholds"]


def metrics(a_bar, a_cross):
    Tb, Tc = np.abs(a_bar) ** 2, np.abs(a_cross) ** 2
    tot = Tb + Tc
    dphi = np.angle(a_cross / a_bar, deg=True)
    return {"T_bar": Tb, "T_cross": Tc, "total": tot,
            "excess_loss_db": float(np.max(-10*np.log10(np.clip(tot, 1e-20, None)))),
            "imbalance": float(np.max(np.abs(Tb - Tc))),
            "phase_deg": float(np.mean(dphi)),
            "phase_max_dev_deg": float(np.max(np.abs(dphi - THR["control_phase_deg"]))),
            "phase_spread_deg": float(np.ptp(dphi))}


mh = metrics(cb, cc)

# --- the comparison: stored full-device FDTD on the SAME geometry ----------
ref_path = os.path.join(
    FDTD_REF_DIR,
    f"smatrix_waisted_Ww{WW:g}_Lw{LW:g}_Lt{LT:g}_L{L_WAIST:g}.pkl")
if not os.path.exists(ref_path):
    raise SystemExit(
        f"no full-device FDTD reference for this geometry:\n  {ref_path}\n"
        "Run mmi_2x2_waist_fdtd_reference.py at the same Ww/Lw/Lt/L first "
        "(MMI2X2_REF_WW/_LW/_LT/_L, MMI2X2_REF_CASES=waisted).")
sm = SMatrix.load(ref_path)
fb = sm.coeff(port="out_0", in_mode=1, out_mode=1)
fc = sm.coeff(port="out_1", in_mode=1, out_mode=1)
# The reference ran 51 wl over the same band. Rather than hard-coding the
# 51->11 stride, match each hybrid wavelength to its nearest reference
# wavelength and REFUSE if any pair is not essentially identical -- comparing
# two engines at wavelengths that merely nearly agree would fold a wavelength
# offset into the engine disagreement being measured.
ref_wl = np.linspace(WL_MIN, WL_MAX, len(fb))
idx = [int(np.argmin(np.abs(ref_wl - w))) for w in wavelengths]
mismatch = float(np.max(np.abs(ref_wl[idx] - wavelengths)))
if mismatch > 1e-9:
    raise SystemExit(
        f"wavelength grids do not align: worst |dwl| = {mismatch:.3e} um "
        f"(hybrid {len(cb)} pts vs reference {len(fb)}). Use an N_WL whose "
        f"grid is a subset of the reference's, e.g. N_WL=11 against 51.")
fb, fc = fb[idx], fc[idx]
print(f"[compare] matched {len(idx)} wavelengths to the stored reference "
      f"(worst |dwl| = {mismatch:.1e} um)")
mf = metrics(fb, fc)

dTb = np.abs(mh["T_bar"] - mf["T_bar"])
dTc = np.abs(mh["T_cross"] - mf["T_cross"])
print("\n" + "=" * 88)
print("HYBRID vs FULL-DEVICE FDTD -- identical geometry, identical wavelengths")
print("=" * 88)
print(f"{'port':>6} {'max|dT|':>10} {'mean|dT|':>10} {'% of mean T_fdtd':>17} "
      f"{'worst wl':>10}")
for nm, d, ref in (("bar", dTb, mf["T_bar"]), ("cross", dTc, mf["T_cross"])):
    i = int(np.argmax(d))
    print(f"{nm:>6} {d.max():10.4f} {d.mean():10.4f} "
          f"{100*d.mean()/ref.mean():17.2f} {wavelengths[i]:10.4f}")
print(f"\n{'quantity':>18} {'hybrid':>10} {'full FDTD':>10} {'delta':>10}")
for q in ("excess_loss_db", "imbalance", "phase_deg", "phase_max_dev_deg",
          "phase_spread_deg"):
    print(f"{q:>18} {mh[q]:10.4f} {mf[q]:10.4f} {mh[q]-mf[q]:+10.4f}")
print("\nStraight 2x2 for scale (gap 5, job 11585849, CALIBRATED body):")
print("   bar 0.62% / cross 0.59% of mean T -- this device's hybrid cannot")
print("   calibrate, so a LARGER disagreement here is the expected outcome,")
print("   and its size is the measurement.")

path = os.path.join(output_dir,
                    # N_MODES belongs in the name: the nm10 merge silently
                    # OVERWROTE the nm14 comparison on 2026-09-10 because this
                    # carried only res and C_TRANS. Recovered by re-merging the
                    # nm14 shards, which are themselves nm-tagged. Seventh
                    # filename-collision in this study -- the standing rule is
                    # that an output name must carry every parameter the run
                    # is fingerprinted on.
                    f"hybrid_vs_fdtd{_DRY}_res{RES:g}_ct{C_TRANS}"
                    f"_nm{N_MODES}{_SM}{_RU}.json")
with open(path, "w") as fh:
    json.dump({"config": {"res": RES, "n_wl": N_WL, "c_trans": C_TRANS,
                          "decay_tol": DECAY_TOL, "n_modes": N_MODES,
                          "symmetric_mode_reuse": SYMMETRIC_MODE_REUSE,
                          "calibrate_body_neff": False},
               "design": {"Ww": WW, "Lw": LW, "Lt": LT, "L": L_WAIST},
               "seconds": {"bar": t_bar, "cross": t_cross},
               "wavelengths": wavelengths.tolist(),
               "hybrid": {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                          for k, v in mh.items()},
               "full_fdtd": {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                             for k, v in mf.items()},
               "max_abs_dT": {"bar": float(dTb.max()),
                              "cross": float(dTc.max())}}, fh, indent=2)
print(f"\nwritten -> {path}")
