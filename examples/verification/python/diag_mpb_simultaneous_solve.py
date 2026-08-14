"""TEMPORARY -- one-off feasibility diagnostic for Stage 2 item (b) of
mode_crossing_correspondence_repair_validation.md: "consider whether
`_iter_distinct_bands()`'s ladder-search approach needs a fundamentally
different strategy for this cross-section (e.g. a simultaneous/degenerate-
subspace solve on the meep side...)."

The question this settles
-------------------------
`_iter_distinct_bands()` calls `sim.get_eigenmode(f, ..., band_num, kpoint,
...)` once per logical band, walking band_num = 1, 2, 3, ... Checked
directly against meep's C++ source (`src/mpb.cpp`, `NanoComp/meep`): each
call ALREADY does a simultaneous block diagonalization of bands
[1..band_num] internally (`create_evectmatrix(..., band_num, ...)` feeding
one `eigensolver()` call) -- Python's `get_eigenmode()` just discards the
lower band_num-1 eigenvectors that same call already computed and returns
only the top one.

So "switch to a simultaneous solve" is not a real change from what already
happens INSIDE one call. What's actually different between the current
ladder and a real fix: band_num=13's call and band_num=14's call are TWO
INDEPENDENT diagonalizations (different block sizes, no shared state) --
nothing guarantees they split a near-degenerate subspace the same way. If
that inter-call inconsistency is the actual mechanism behind the
band13==band14 exact duplicate found at every domain size in job 11813373
(diag_dup_threshold_self_overlap.py), then extracting ALL N bands from a
SINGLE diagonalization (rather than N independent single-band-terminal
calls) should keep them mutually distinct, since orthogonality among
outputs of one block eigensolve is a mathematical guarantee, not a
probabilistic byproduct.

This script tests exactly that, using `meep.mpb.ModeSolver` directly (a
lower-level, separate Python binding onto the same MPB engine, confirmed
importable in the `jmh` env and confirmed to accept `mp.Prism`/`mp.Block`
geometry objects directly -- `mp.Prism`'s MRO is
`(meep.geom.Prism, meep.geom.GeometricObject, object)`, the same base
class `ModeSolver.geometry` expects). One `run()` at a FIXED k (not a
fixed-frequency root search -- that distinction does not matter for the
orthogonality question, only for how close the resulting frequencies land
to the target 1.26 um, which is not what's being tested here) with
`num_bands=16`, then the full 16x16 pairwise overlap matrix among that
ONE run's own output bands, using the same normalized-overlap formula
`beef.mode_identity.field_overlap_matrix` uses (re-implemented against
MPB's own native field grid rather than beef's diagnostic grid -- this is
a self-consistency check among one solve's own outputs, not a cross-grid
comparison, so no need to force MPB's fields onto beef's separate
diag_y/diag_z sampling convention).

Known fidelity gaps in this feasibility check (acceptable for a go/no-go
answer, would need resolving before any production integration):
  - MPB has no PML. The y/z supercell is genuinely PERIODIC in MPB's
    plane-wave basis; meep's own FDTD-side mode solve uses PML instead.
    A guided mode's evanescent tail wrapping around a periodic boundary is
    a DIFFERENT truncation error than reflecting off a PML, not
    necessarily smaller for the same box size. Using sim_width=5.0 (the
    baseline where the duplicate is confirmed rock-solid in meep's own
    ladder) rather than a larger box, so this is not conflated with the
    domain-size question Stage 0 already chased.
  - Fixed k rather than fixed frequency: the seed kx is meep's own
    converged k for the known-duplicated pair (band13/14 at sim_width=5.0,
    n_eff=1.657513, kx = n_eff/wl), so results land near, but not exactly
    at, wl=1.26 um for every band -- irrelevant to the orthogonality
    question this script asks.

Cost: one MPB `run()`, num_bands=16, single k-point, real device
resolution/geometry -- comparable to what meep's OWN band_num=16 ladder
call already internally pays for its last (most expensive) step, so this
is not a new order of cost, just done once instead of inside 16 separate
Python-level calls.
"""
import os
import sys

import numpy as np

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "..", "modal_multiplexer",
                       "python", "modal_multiplexer_hybrid.py")
SCRIPT = os.path.abspath(SCRIPT)
MARKER = "hybrid = HybridSolver("

with open(SCRIPT) as f:
    src = f.read()

end_marker = "\n)\n"
start = src.index(MARKER)
end = src.index(end_marker, start) + len(end_marker)
truncated = src[:end]

print(f"truncated source: {len(truncated)} of {len(src)} chars", flush=True)

ns = {"__name__": "__main__", "__file__": SCRIPT}
code = compile(truncated, SCRIPT, "exec")
exec(code, ns)

hybrid = ns["hybrid"]
mux = ns["mux"]
fdtd_solver = ns["fdtd_solver"]

import meep as mp
from meep import mpb

WL = 1.26
N_EFF_DUP = 1.657513  # meep's own band13==band14 value at sim_width=5.0, job 11813373
F = 1.0 / WL
KX_SEED = N_EFF_DUP * F
NUM_BANDS = 16  # 14 (shipped N_MODES) + 2 cushion, matching _DUP_BAND_MAX_EXTRA's spirit
SIM_WIDTH = 5.0  # shipped baseline -- the domain size where the duplicate is confirmed
                 # rock-solid in meep's own ladder (job 11813373), deliberately NOT the
                 # larger, previously-explored domain sizes, so this result is not
                 # conflated with the domain-size question Stage 0 already answered.

X_IFACE = hybrid.boundaries[1]  # interface 0: FDTD-in <-> EME-MMI1
HALF_SPAN = 0.05  # matches HybridSolver's own _CORR_HALF_SPAN
X0, X1 = X_IFACE - HALF_SPAN, X_IFACE + HALF_SPAN

print(f"interface 0 at x={X_IFACE:.4f}, window [{X0:.4f}, {X1:.4f}], wl={WL}, "
      f"kx_seed={KX_SEED:.6f} (from n_eff={N_EFF_DUP}), num_bands={NUM_BANDS}, "
      f"sim_width={SIM_WIDTH}", flush=True)

# ---- build the SAME cross-section geometry _build_section_mode_sim() does ----
clearance = fdtd_solver.dpml  # only matters here as prism stub length, MPB has no PML
si_med, sio2_med = fdtd_solver._media(WL)
prisms = mux.section_prisms(
    X0, X1, si_med, stub_left=clearance, stub_right=clearance,
    wg_height=fdtd_solver._wg_height,
)
x_center = (X0 - clearance + X1 + clearance) / 2.0
shifted = [
    mp.Prism(
        vertices=[mp.Vector3(v.x - x_center, v.y, v.z) for v in prism.vertices],
        height=prism.height, material=prism.material,
    )
    for prism in prisms
]
print(f"built {len(shifted)} prisms, x_center={x_center:.4f}", flush=True)

sim_height = fdtd_solver.sim_height
mon_sz = fdtd_solver._mon_sz  # real device z-span, no PML padding -- MPB has none anyway
print(f"MPB supercell (y,z) = ({SIM_WIDTH}, {mon_sz}), resolution={fdtd_solver.res}",
      flush=True)

ms = mpb.ModeSolver(
    geometry_lattice=mp.Lattice(size=mp.Vector3(0, SIM_WIDTH, mon_sz)),
    geometry=shifted,
    default_material=sio2_med,
    resolution=fdtd_solver.res,
    num_bands=NUM_BANDS,
    k_points=[mp.Vector3(KX_SEED, 0, 0)],
)

print("\n===== running MPB (one simultaneous diagonalization, all bands together) =====",
      flush=True)
ms.run()

freqs = np.asarray(ms.all_freqs)[0]  # [NUM_BANDS]
neff = KX_SEED / freqs
print(f"\nfrequencies: {np.array2string(freqs, precision=6)}", flush=True)
print(f"n_eff = kx_seed/freq: {np.array2string(neff, precision=6)}", flush=True)
print("(compare against meep's own ladder at sim_width=5.0: bands 1-13 clean, "
      "band 14 duplicated band 13 at n_eff=1.657513 -- job 11813373)", flush=True)

# ---- pull every band's Ey/Ez on MPB's own native grid, compute full pairwise overlap ----
print("\n===== extracting fields, computing full pairwise overlap matrix =====", flush=True)
efields = []
for b in range(1, NUM_BANDS + 1):
    ef = ms.get_efield(b)  # (ny, nz, nx=1, 3) complex, components (Ex, Ey, Ez)
    efields.append(ef[:, :, 0, :])  # drop the trivial nx=1 axis -> (ny, nz, 3)
print(f"sampled {len(efields)} bands, each field shape {efields[0].shape}", flush=True)


def self_overlap(ef):
    ey, ez = ef[..., 1], ef[..., 2]
    return np.sum(np.abs(ey) ** 2 + np.abs(ez) ** 2)


def cross_overlap(ef_a, ef_b):
    ey_a, ez_a = ef_a[..., 1], ef_a[..., 2]
    ey_b, ez_b = ef_b[..., 1], ef_b[..., 2]
    return np.sum(np.conj(ey_a) * ey_b + np.conj(ez_a) * ez_b)


self_p = np.array([self_overlap(ef) for ef in efields])
eta = np.zeros((NUM_BANDS, NUM_BANDS))
for i in range(NUM_BANDS):
    for j in range(NUM_BANDS):
        kappa = cross_overlap(efields[i], efields[j])
        denom = np.sqrt(self_p[i] * self_p[j])
        eta[i, j] = abs(kappa) / denom if denom > 0 else 0.0

print("\nfull eta (normalized overlap) matrix, MPB's own 16 simultaneously-solved bands:",
      flush=True)
header = "      " + " ".join(f"b{b:2d}  " for b in range(1, NUM_BANDS + 1))
print(header, flush=True)
for i in range(NUM_BANDS):
    row = " ".join(f"{eta[i, j]:.3f}" for j in range(NUM_BANDS))
    print(f"b{i+1:2d}:  {row}", flush=True)

off_diag = eta.copy()
np.fill_diagonal(off_diag, 0.0)
worst_i, worst_j = np.unravel_index(np.argmax(off_diag), off_diag.shape)
worst = off_diag[worst_i, worst_j]

print(f"\nworst off-diagonal overlap among all {NUM_BANDS} simultaneously-solved "
      f"bands: {worst:.6f}  (bands {worst_i + 1} and {worst_j + 1}, "
      f"n_eff {neff[worst_i]:.6f} and {neff[worst_j]:.6f})", flush=True)
print(f"threshold this would be judged against (production's own): 0.90", flush=True)

print("\n===== how to read this =====", flush=True)
print("  worst off-diagonal << 0.90  ->  a genuine simultaneous diagonalization keeps "
      "\n    every band distinct where meep's independent per-band-count ladder did not."
      "\n    Confirms the inter-call-inconsistency hypothesis; the fix is real: extract "
      "\n    all N bands from ONE call instead of N independent ones.", flush=True)
print("  worst off-diagonal still ~1.0 near n_eff~1.6575  ->  simultaneity alone does "
      "\n    NOT fix it -- points to a genuine near-degeneracy or periodic-vs-PML "
      "\n    domain-truncation artifact instead. Important negative result, not a dead "
      "\n    end: narrows what Stage 2's real fix needs to address.", flush=True)

print("\nDIAGNOSTIC DONE", flush=True)
