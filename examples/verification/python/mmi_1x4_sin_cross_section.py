"""Cross-section design for the 1x4 MMI on 300 nm Si3N4 (SIN_300) --
validation_markdown/mmi_1x4_sin300_validation.md step 4.

Two independent 1D mode-solve sweeps (meow's `compute_modes_tidy3d` FDE
backend -- no FDTD, no time-stepping).

⚠ NOT `compute_modes_meep`, despite that being this repo's usual "cheap mode
solve" backend (directional_coupler_mode.py, verify_item_8_strip_gate.py):
`compute_modes_meep` (meow_MEEP/src/meow/fde/meep.py) searches for eigenmodes
via meep's own band-structure machinery (`sim.get_eigenmode`, `direction=
NO_DIRECTION`, Newton iteration on k from a fixed initial guess
`kpoint=Vector3(z=1)`). MEASURED 2026-08-03: for this 300 nm Si3N4 strip it
returns ONLY spurious box/continuum modes (n_eff < n_clad) at every (width,
sz) tried, with NO error or warning -- a silent wrong answer, not a crash.
Confirmed this is not a bad-guess or box-size problem: (1) an analytic
symmetric-slab check predicts a real, tightly-confined bound mode at
n_eff=1.673911 with a 0.29 um cladding decay length, comfortably inside every
box tested; (2) the built structure's index map was verified correct
(1.9963 core / 1.4440 clad, exactly at the drawn boundaries); (3) a Si control
at the IDENTICAL box finds 3 clean guided modes immediately via this same
solver, ruling out a general box/solver bug; (4) explicitly forcing three
different initial k-guesses (1.0, 1.65, 1.9 -- spanning from the hardcoded
default to essentially the core index) into meep's own `get_eigenmode` for
this exact structure converges to the SAME wrong value (1.440790) every time,
which rules out "bad initial guess" specifically and points at something more
fundamental in how that search behaves for weakly-confined structures.

`compute_modes_tidy3d` -- meow's OTHER backend, a proper finite-difference
eigenmode solver rather than a k-vs-frequency band search -- finds the real
modes immediately on the identical cross-section (2 guided, n_eff=1.5727/
1.5109 at width=1.0 um). Cross-validated against `compute_modes_meep` on the
Si control first: the two backends agree to ~1-3% there (2.938/2.630/2.501
vs 2.953/2.651/2.535), consistent with this repo's other documented
engine-agreement gaps (MPB-vs-meep phase ~2.3%, MMI hybrid-vs-FDTD ~0.5-1%),
not a red flag -- so this is a genuine, working alternative solver, not an
unrelated or broken one. Used here instead of patching `compute_modes_meep`'s
search internals, which are shared by every existing (silicon) device and
whose current behaviour those devices already depend on.

  4a. VERTICAL SPAN (sz). MaterialStack.sz is the total z-extent meow
      extrudes cladding into (MaterialStack.extrusion_rules_at()) -- nothing
      to do with the in-plane mesh width, which each solver sets separately
      (oxide_width/sim_width). SIN_300 currently inherits sz=2.0 from the SOI
      stacks; its own comment in beef/materials/stack.py flags this as an
      unvalidated placeholder. Nitride's core/clad eps contrast is ~5.3x
      weaker than Si/SiO2's (this doc's A7), so the mode reaches further into
      the cladding and 2.0 um may clip it. Sweep sz, watch the fundamental
      mode's n_eff, and find where it stops moving.

  4b. SINGLE-MODE WIDTH. At the sz step 4a settles on, sweep strip width and
      COUNT GUIDED MODES (Re(n_eff) > n_clad -- the same convention
      mode_spectrum_convergence.py and HybridSolver's mode-correspondence
      guided-mode filter both use; README.md Finding 1 is why this matters:
      "every interface basis in this repo carries at least one unguided
      mode", so a naive count of RETURNED modes over-counts). Report the
      widest width that stays single-mode.

Run as two 1D sweeps, not a joint 2D optimisation: sz mostly sets vertical
confinement, width mostly sets lateral, and treating them independently (with
one confirmatory re-check of sz at the width 4b lands on) is what the
validation doc's step 4 asks for. Both use C-band at 1.55 um -- see 4c's
writeup in the validation doc for why (settled from step 1's A6 measurements,
no mode solve needed for that half).

Cost note: single-strip, single-wavelength eigensolves at res 40 -- the same
resolution verify_item_8_strip_gate.py uses (that script's compute_modes_meep
solves are "seconds each"; compute_modes_tidy3d measured somewhat slower per
call in testing, still well within the batch job's -t 240 budget). This
resolution is NOT itself validated for nitride (that is check 3 in the
validation doc's step 8, out of scope here) -- it is reused only because it is
this codebase's standing convention for a cheap first-pass mode solve.

Run:
    python /path/to/examples/verification/python/mmi_1x4_sin_cross_section.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import meow as mw
from meow import compute_modes_tidy3d

from beef import extrude_polygons
from beef.materials.stack import MaterialStack, sellmeier_index

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "example_output", "verification", "mmi_1x4_sin_cross_section"))
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Constants shared by both sweeps.
# ---------------------------------------------------------------------------
WG_HEIGHT = 0.30    # um, 300 nm Si3N4 core -- fixed by the user's request
WL = 1.55           # um, C-band centre (band choice settled in the doc, 4c)
RES = 40.0          # pts/um -- this codebase's standing mode-solve default;
                    # NOT itself re-validated for nitride here (check 3, step 8)
N_CLAD = float(sellmeier_index("sio2", WL))   # guided-mode threshold

# 4a: candidate vertical spans (mesh height, "MESH_HEIGHT" == MaterialStack.sz)
# MEASURED 2026-08-03 (smoke test, real compute_modes_tidy3d, sz=2.0/2.5 at
# width=1.0): |dn_eff|=1.84e-3 between those two -- an order of magnitude
# above SZ_CONVERGENCE_TOL, so the original {2.0,2.5,3.0,4.0,5.0} range almost
# certainly does not reach convergence. Extended upward accordingly: the
# analytic infinite-slab estimate's 0.29 um cladding decay length (this doc's
# 4a writeup) suggests several decay lengths of margin past the core edge are
# needed each side, i.e. sz well past 5 um, not comfortably inside it as
# originally assumed.
SZ_CANDIDATES = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0]
# Fixed reference width for the sz sweep. NOT a result -- a scaling-argument
# seed, stated so it's checkable rather than silently assumed: matching a
# normalised-frequency V ~ width*sqrt(delta_n^2) to silicon's known single-mode
# point (450 nm at delta_n^2=10.0093) under nitride's delta_n^2=1.8999 (this
# doc's A7) gives width ~= 0.45 * sqrt(10.0093/1.8999) ~= 1.03 um. Whether this
# width is ACTUALLY single-mode is exactly what 4b checks; 4a only needs "some
# plausible single strip", since sz mostly sets VERTICAL confinement.
SZ_SWEEP_WIDTH = 1.0
SZ_CONVERGENCE_TOL = 1e-4   # |dn_eff| between successive sz -- see doc for why

# 4b: candidate widths once sz has converged.
# MEASURED 2026-08-03 (smoke test): width=0.6 um ALREADY guides 2 modes
# (n_eff=1.498/1.477), so the original {0.6, ..., 2.4} range never brackets a
# single-mode point -- extended downward to actually find the cutoff.
WIDTH_CANDIDATES = [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2]
MAX_MODES_REQUESTED = 8     # width=2.4 measured 5 guided of 6 requested in the
                            # smoke test -- bumped up so the widest candidates
                            # don't run out of headroom to count correctly

MESH_WIDTH = 8.0            # um, lateral mesh extent -- generous clearance
                             # past the widest strip tested (2.4 um) plus
                             # several decay lengths of cladding each side


def _stack(sz: float) -> MaterialStack:
    return MaterialStack(wg_height=WG_HEIGHT, sz=sz, core_kind="si3n4")


def _axis_grid(half_span: float, res: float) -> np.ndarray:
    """Symmetric grid over [-half_span, half_span] at EXACTLY 1/res spacing.

    Originally required because `compute_modes_meep` needs an isotropic mesh
    (equal dx/dy -- meep's own resolution is a single scalar); this script has
    since switched to `compute_modes_tidy3d` (see the module docstring), which
    has no such restriction. Kept anyway as a good default rather than
    reverting to independent per-axis `linspace(-L/2, L/2, int(res*L)+1)`
    calls: that pattern only gives matching dx/dy when `res*L` happens to land
    on an integer for BOTH axes' L, which this script's swept `sz` is not
    guaranteed to do (e.g. sz=2.2 at res=8 gives dy=2.2/17=0.1294 against
    dx=4.0/32=0.125 -- exactly what broke the original compute_modes_meep
    version). Building both axes from the SAME explicit step guarantees
    dx == dy for any half_span regardless of which backend reads it.
    """
    n_half = int(round(half_span * res))
    step = 1.0 / res
    return step * np.arange(-n_half, n_half + 1)


def build_strip_cross_section(width: float, sz: float, wl: float) -> mw.CrossSection:
    """One Si3N4 strip in SiO2 -- single-waveguide analogue of
    directional_coupler_mode.py's build_cross_section() (which builds TWO
    strips for supermode analysis; this device needs one)."""
    length = 1.0
    wg_poly = [
        (0, -width / 2), (length, -width / 2),
        (length, width / 2), (0, width / 2),
    ]
    ox_poly = [
        (0, -MESH_WIDTH / 2), (length, -MESH_WIDTH / 2),
        (length, MESH_WIDTH / 2), (0, MESH_WIDTH / 2),
    ]
    gds_polys = {(2, 0): [wg_poly], (1, 0): [ox_poly]}

    structs = extrude_polygons(gds_polys, _stack(sz).extrusion_rules)
    mesh = mw.Mesh2D(x=_axis_grid(MESH_WIDTH / 2, RES), y=_axis_grid(sz / 2, RES))
    cells = mw.create_cells(structures=structs, mesh=mesh, Ls=[length])
    env = mw.Environment(wl=float(wl), T=25.0)
    return mw.CrossSection.from_cell(cell=cells[0], env=env)


def guided_neffs(width: float, sz: float, wl: float, n_clad: float,
                 max_modes: int = MAX_MODES_REQUESTED):
    """All requested modes' n_eff (sorted, descending), and the GUIDED subset
    (n_eff > n_clad -- see module docstring for why this filter matters)."""
    cs = build_strip_cross_section(width, sz, wl)
    modes = compute_modes_tidy3d(cs, num_modes=max_modes)
    neffs = sorted((float(np.real(m.neff)) for m in modes), reverse=True)
    guided = [n for n in neffs if n > n_clad]
    return neffs, guided


# ---------------------------------------------------------------------------
# 4a. sz (vertical span) convergence sweep.
# ---------------------------------------------------------------------------
print("=" * 78)
print(f"4a. Vertical span (sz) sweep -- {WG_HEIGHT*1000:.0f} nm Si3N4, "
      f"{SZ_SWEEP_WIDTH} um wide strip, lambda={WL} um")
print(f"    n_clad (SiO2 @ {WL} um) = {N_CLAD:.6f}")
print("=" * 78)

n_eff_sz = np.zeros(len(SZ_CANDIDATES))
for i, sz in enumerate(SZ_CANDIDATES):
    neffs, guided = guided_neffs(SZ_SWEEP_WIDTH, sz, WL, N_CLAD)
    n_eff_sz[i] = neffs[0]
    print(f"  sz = {sz:.2f} um  ->  n_eff(fundamental) = {neffs[0]:.6f}   "
          f"({len(guided)} guided of {len(neffs)} requested)   "
          f"all modes: {[f'{n:.4f}' for n in neffs]}")

d_sz = np.abs(np.diff(n_eff_sz))
print(f"\n  Successive |dn_eff|: {[f'{d:.2e}' for d in d_sz]}")
# converged_idx = i means the STEP FROM SZ_CANDIDATES[i] TO SZ_CANDIDATES[i+1]
# (i.e. d_sz[i]) cleared the tolerance -- SZ_CANDIDATES[i] is then the
# smallest sz past which extending further stops mattering. Deliberately
# `is not None` throughout below, not a truthy check: converged_idx == 0
# (convergence on the very first step) is a valid, falsy-looking value.
converged_idx = None
for i, d in enumerate(d_sz):
    if d < SZ_CONVERGENCE_TOL:
        converged_idx = i
        break
if converged_idx is not None:
    sz_converged = SZ_CANDIDATES[converged_idx]
    print(f"  -> CONVERGED at sz = {sz_converged:.2f} um "
          f"(|dn_eff| < {SZ_CONVERGENCE_TOL:.0e} going from "
          f"{SZ_CANDIDATES[converged_idx]:.2f} -> "
          f"{SZ_CANDIDATES[converged_idx+1]:.2f} um)")
else:
    sz_converged = SZ_CANDIDATES[-1]
    print(f"  -> NOT CONVERGED within the swept range; largest tested "
          f"sz = {sz_converged:.2f} um did not clear {SZ_CONVERGENCE_TOL:.0e}. "
          f"Extend SZ_CANDIDATES before trusting this value.")

# ---------------------------------------------------------------------------
# 4b. Single-mode width sweep, at the converged sz.
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"4b. Single-mode width sweep -- sz = {sz_converged:.2f} um "
      f"(from 4a), lambda={WL} um")
print("=" * 78)

n_guided_width = np.zeros(len(WIDTH_CANDIDATES), dtype=int)
n_eff_width = np.zeros(len(WIDTH_CANDIDATES))
all_neffs_width = []
for i, w in enumerate(WIDTH_CANDIDATES):
    neffs, guided = guided_neffs(w, sz_converged, WL, N_CLAD)
    n_guided_width[i] = len(guided)
    n_eff_width[i] = neffs[0]
    all_neffs_width.append(neffs)
    print(f"  width = {w:.2f} um  ->  {len(guided)} guided mode(s), "
          f"n_eff = {[f'{n:.4f}' for n in guided]}   "
          f"(all requested: {[f'{n:.4f}' for n in neffs]})")

single_mode_widths = [w for w, n in zip(WIDTH_CANDIDATES, n_guided_width) if n <= 1]
multi_mode_widths = [w for w, n in zip(WIDTH_CANDIDATES, n_guided_width) if n > 1]
if single_mode_widths and multi_mode_widths:
    w_max_single = max(single_mode_widths)
    w_min_multi = min(multi_mode_widths)
    print(f"\n  -> Single-mode for width <= {w_max_single:.2f} um; "
          f"multi-mode from width >= {w_min_multi:.2f} um.")
    print(f"     Cutoff bracketed in [{w_max_single:.2f}, {w_min_multi:.2f}] um "
          f"(not interpolated -- refine WIDTH_CANDIDATES for a tighter bracket "
          f"if the design needs it).")
elif single_mode_widths:
    w_max_single = max(single_mode_widths)
    print(f"\n  -> Every tested width up to {w_max_single:.2f} um is single-mode; "
          f"no multi-mode transition found in range. Extend WIDTH_CANDIDATES "
          f"upward before picking a design width near the top of this list.")
else:
    w_max_single = None
    print(f"\n  -> NO tested width is single-mode (smallest tested "
          f"{WIDTH_CANDIDATES[0]:.2f} um already guides >1 mode). Extend "
          f"WIDTH_CANDIDATES downward.")

# ---------------------------------------------------------------------------
# Confirmatory re-check: does sz convergence hold at the actual single-mode
# width, not just at the SZ_SWEEP_WIDTH guess?
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("Confirmatory check: re-run the sz sweep at the widest single-mode "
      "width found above")
print("=" * 78)
if w_max_single is not None and w_max_single != SZ_SWEEP_WIDTH:
    n_eff_sz_confirm = np.zeros(len(SZ_CANDIDATES))
    for i, sz in enumerate(SZ_CANDIDATES):
        neffs, _ = guided_neffs(w_max_single, sz, WL, N_CLAD)
        n_eff_sz_confirm[i] = neffs[0]
        print(f"  sz = {sz:.2f} um  ->  n_eff(fundamental) = {neffs[0]:.6f}")
    d_confirm = np.abs(np.diff(n_eff_sz_confirm))
    print(f"\n  Successive |dn_eff| at width={w_max_single:.2f} um: "
          f"{[f'{d:.2e}' for d in d_confirm]}")
    if converged_idx is None:
        print(f"  -> 4a never converged in-range, so there is no index to "
              f"re-check here either; read the |dn_eff| list above directly.")
    else:
        # Same step 4a converged on: SZ_CANDIDATES[converged_idx] ->
        # SZ_CANDIDATES[converged_idx+1], i.e. d_confirm[converged_idx].
        confirm_ok = bool(d_confirm[converged_idx] < SZ_CONVERGENCE_TOL)
        print(f"  -> sz={sz_converged:.2f} um still clears {SZ_CONVERGENCE_TOL:.0e} "
              f"at this width: {'YES' if confirm_ok else 'NO -- re-examine'}")
else:
    n_eff_sz_confirm = None
    print("  Skipped: no distinct single-mode width to re-check against "
          "(see 4b's result above).")

# ---------------------------------------------------------------------------
# Save data + plots.
# ---------------------------------------------------------------------------
np.savez(
    os.path.join(output_dir, "cross_section_sweep.npz"),
    wg_height=WG_HEIGHT, wl=WL, res=RES, n_clad=N_CLAD,
    sz_candidates=np.asarray(SZ_CANDIDATES), sz_sweep_width=SZ_SWEEP_WIDTH,
    n_eff_sz=n_eff_sz, sz_converged=sz_converged,
    width_candidates=np.asarray(WIDTH_CANDIDATES),
    n_guided_width=n_guided_width, n_eff_width=n_eff_width,
    n_eff_sz_confirm=(n_eff_sz_confirm if n_eff_sz_confirm is not None
                      else np.array([])),
    w_max_single=(w_max_single if w_max_single is not None else np.nan),
)
print(f"\nData saved -> {output_dir}/cross_section_sweep.npz")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

ax = axes[0]
ax.plot(SZ_CANDIDATES, n_eff_sz, "o-", color="steelblue", label=f"width={SZ_SWEEP_WIDTH} um")
if n_eff_sz_confirm is not None:
    ax.plot(SZ_CANDIDATES, n_eff_sz_confirm, "s--", color="darkorange",
           label=f"width={w_max_single:.2f} um (confirm)")
ax.axvline(sz_converged, color="gray", linestyle=":", lw=1,
          label=f"chosen sz={sz_converged:.2f} um")
ax.set_xlabel("sz -- vertical simulation span (um)")
ax.set_ylabel("n_eff (fundamental mode)")
ax.set_title(f"4a. Vertical span convergence ({WG_HEIGHT*1000:.0f} nm Si3N4)")
ax.legend(fontsize=8)
ax.grid(True)

ax = axes[1]
ax.plot(WIDTH_CANDIDATES, n_guided_width, "o-", color="firebrick")
ax.axhline(1.5, color="gray", linestyle=":", lw=1, label="single/multi-mode boundary")
ax.set_xlabel("Strip width (um)")
ax.set_ylabel("Guided mode count")
ax.set_title(f"4b. Single-mode width (sz={sz_converged:.2f} um)")
ax.set_yticks(range(0, int(n_guided_width.max()) + 2))
ax.legend(fontsize=8)
ax.grid(True)

plt.suptitle(f"1x4 MMI on 300 nm Si3N4 -- cross-section design (lambda={WL} um)",
            fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "cross_section_sweep.png"), dpi=300)
plt.close("all")
print(f"Plot saved -> {output_dir}/cross_section_sweep.png")
