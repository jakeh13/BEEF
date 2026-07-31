"""Resolution convergence of the meep/meow MODE SPECTRUM at one cross-section
(validation_markdown/fdtd_vs_EME_mode_validation.md item 1, and
modal_multiplexer_validation.md gap 2).

NOT item 6 (`interface_modes` convergence, layer 2), which an earlier version of
this docstring claimed: that one sweeps the requested MODE COUNT and requires
the S-MATRIX to plateau. This sweeps resolution at a fixed count and looks at
the spectrum. See item 6's own note for what this run does and does not
contribute to it.

Result of the full run (job 11575516, 2026-07-30) is recorded under item 1 --
read that before re-running this: the mode-8+ divergence shrinks with
resolution but does not close, and the residue is one member of a near-
degenerate pair that meep's per-band ladder never returns.

The question this answers
------------------------
At the modal multiplexer's interface 0, with the EME duplicate-mode guard in
place (item 0), the two engines' n_eff spectra agree to ~1e-4 for modes 0-7 and
then diverge systematically from mode 8 up:

    meep guided orphan : 1.74775
    meow guided orphans: 1.78589, 1.51288
    (meep finds 13 guided modes, meow 14)

Both engines run MPB, so a disagreement this large is either (a) a
discretization artifact -- they run on differently CONSTRUCTED grids at the
same nominal res (meep's own Yee grid with subpixel smoothing vs the Mesh2D
meow is handed), and higher-order modes are far more grid-sensitive because
they oscillate faster -- or (b) something structural, i.e. a real bug.

Those two have opposite remedies, and the discriminator is simple: **sweep the
resolution and watch whether the upper-half disagreement shrinks.**

    shrinks monotonically -> discretization. The fix is a per-engine converged
        resolution (never a single shared knob -- that is
        directional_coupler_validation.md item 11's defect), and the honest
        statement is that the top of a 16-mode SOI_160 basis is simply not
        converged at res 40.
    plateaus at a finite gap -> structural. A finer mesh cannot reach it and
        something is actually wrong; escalate.

Reading the output
------------------
Per resolution, for each engine: the full n_eff spectrum, the guided-mode count
(n_eff > cladding index), and -- the headline -- how many modes the two engines
agree on to within a tolerance. That last number rising with resolution is the
signature of (a).

`n_agree` is deliberately a COUNT rather than a mean |dn_eff|: the failure mode
here is a mode appearing in one spectrum and not the other, which shifts every
subsequent index and makes an element-wise mean meaningless (it would compare
mode k against mode k+1). Counting matched modes is robust to that offset --
the same reason the field-overlap work uses an assignment rather than an
element-wise difference.

Cost note
---------
Uses `section_mode_neff()` / `compute_modes()`, NOT `section_mode_fields()`:
the eigensolve is the cheap part, and the earlier 15-minute runs were dominated
by field sampling (16 modes x 200 x 60 points x 6 components ~ 1.15M
`amplitude()` calls). Spectra alone answer the convergence question, so this
skips that entirely. Even so it is a batch job -- meow's mesh grows as res^2.
"""
from __future__ import annotations

import numpy as np

from beef import (
    MMI, MMIPort, WaveguideRun, CompoundGeometry, EMESolver, FDTDSolver, SOI_160,
)

WL = 1.31
N_MODES = 16
BODY_STUB = 0.75
SIM_W, SIM_H = 5.0, 1.5
DPML = 0.5
HALF_SPAN = 0.05

# Independent per engine, deliberately: a single shared knob is exactly the
# defect directional_coupler_validation.md item 11 records.
MEEP_RESOLUTIONS = [20.0, 30.0, 40.0, 60.0]
MEOW_RESOLUTIONS = [20.0, 30.0, 40.0, 60.0, 80.0]
MATCH_TOL = 0.005   # |dn_eff| within which two modes are "the same mode"

TAPER_L, WO, BETWEEN, MM_WIDTH = 2.0, 1.0375, 3.0, 2.8
BRANCH_SPECS = [(+WO, 0.33, 0.725), (0.0, 0.35, 1.15), (-WO, 0.35, 0.725)]


def _branch_ports():
    return [MMIPort(y=y, waveguide_width=ww, taper_length=TAPER_L, taper_width=tw)
            for y, ww, tw in BRANCH_SPECS]


def _build():
    mmi1 = MMI(n_inputs=1, n_outputs=3, mm_length=16.92, mm_width=MM_WIDTH,
               inputs=[MMIPort(y=0.0, waveguide_width=0.7, taper_length=TAPER_L,
                               taper_width=1.4)],
               outputs=_branch_ports(), waveguide_height=0.16)
    mmi2 = MMI(n_inputs=3, n_outputs=3, mm_length=34.24, mm_width=MM_WIDTH,
               inputs=_branch_ports(), outputs=_branch_ports(),
               waveguide_height=0.16, use_neighbor_aperture=True)
    connector = WaveguideRun(length=BETWEEN, ports=[(y, ww) for y, ww, _ in BRANCH_SPECS])
    return mmi1, CompoundGeometry([mmi1, connector, mmi2])


def _meep_spectrum(geometry, x_center, res):
    solver = FDTDSolver(res=res, sim_width=SIM_W, sim_height=SIM_H, dpml=DPML,
                        source_margin=0.0, monitor_margin=0.25, three_d=True,
                        material_stack=SOI_160, decay_tol=1e-11)
    neff = solver.section_mode_neff(
        geometry, np.array([WL]), x_center - HALF_SPAN, x_center + HALF_SPAN,
        n_modes=N_MODES, monitor_center=0.0, monitor_size=SIM_W)
    return np.asarray(neff[0], dtype=float)


def _meow_spectrum(geometry, x_center, res):
    eme = EMESolver(num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
                    oxide_width=SIM_W, oxide_height=SIM_H, res=res,
                    material_stack=SOI_160, mode_solver="meep",
                    guard_duplicate_modes=True)
    modes, _ = eme.section_modes(geometry, WL, x_center - HALF_SPAN,
                                 x_center + HALF_SPAN, N_MODES)
    return np.array([float(np.real(m.neff)) for m in modes])


def _match_count(a, b, tol=MATCH_TOL):
    """How many of `a`'s modes have a counterpart in `b` within tol, matched
    one-to-one greedily (nearest first). Robust to an inserted/missing mode,
    which an element-wise difference is not."""
    remaining = list(b)
    n = 0
    for va in a:
        if not remaining:
            break
        j = int(np.argmin([abs(va - vb) for vb in remaining]))
        if abs(va - remaining[j]) <= tol:
            n += 1
            remaining.pop(j)
    return n


def _orphans(a, b, tol=MATCH_TOL):
    return [float(v) for v in a if not any(abs(v - w) <= tol for w in b)]


def main():
    mmi1, mux = _build()
    x_loc = -(mmi1.mm_length / 2 - BODY_STUB)
    n_clad = SOI_160.sio2_index(WL)
    print(f"cross-section: mmi1 local x = {x_loc:.4f} um (interface 0), "
          f"wl = {WL} um, {N_MODES} modes requested")
    print(f"cladding index = {n_clad:.6f}  (guided means n_eff above this)", flush=True)

    meep, meow = {}, {}

    print("\n===== meep (FDTDSolver.section_mode_neff) =====", flush=True)
    for res in MEEP_RESOLUTIONS:
        s = _meep_spectrum(mmi1, x_loc, res)
        meep[res] = s
        print(f"\nres {res:5.1f}: guided {int((s > n_clad).sum())}/{N_MODES}", flush=True)
        print(f"  {np.round(s, 5).tolist()}", flush=True)

    print("\n===== meow (EMESolver.section_modes, guard ON) =====", flush=True)
    for res in MEOW_RESOLUTIONS:
        s = _meow_spectrum(mmi1, x_loc, res)
        meow[res] = s
        print(f"\nres {res:5.1f}: guided {int((s > n_clad).sum())}/{N_MODES}", flush=True)
        print(f"  {np.round(s, 5).tolist()}", flush=True)

    # ---- Per-engine self-convergence: does each spectrum settle? ----------
    for name, spectra in (("meep", meep), ("meow", meow)):
        print(f"\n===== {name}: self-convergence (consecutive resolutions) =====", flush=True)
        keys = sorted(spectra)
        print(f"  {'res_lo':>7} {'res_hi':>7} {'n_agree':>8} {'guided_lo':>10} {'guided_hi':>10}")
        for lo, hi in zip(keys, keys[1:]):
            a, b = spectra[lo], spectra[hi]
            print(f"  {lo:7.1f} {hi:7.1f} {_match_count(a, b):8d} "
                  f"{int((a > n_clad).sum()):10d} {int((b > n_clad).sum()):10d}", flush=True)

    # ---- THE decisive table: MATCHED resolution only ----------------------
    # Only a matched-resolution comparison is meaningful. A resolution change
    # shifts n_eff by more than the mode spacing, so meep@40 vs meow@80 pairs
    # mode k with mode k+1 and looks exactly like a disagreement when both
    # engines are perfectly correct -- tests/test_mode_identity.py pins that
    # trap deliberately. Confirmed here at tiny scale too: matched res gave
    # 3/3 agreement while every mismatched pair gave 0/3.
    print("\n===== MATCHED-RESOLUTION agreement (the decisive table) =====", flush=True)
    print(f"  {'res':>6} {'n_agree':>8} {'guided_meep':>12} {'guided_meow':>12}"
          f"   orphans (meep | meow)", flush=True)
    for res in sorted(set(meep) & set(meow)):
        a, b = meep[res], meow[res]
        oa = [round(v, 5) for v in _orphans(a, b) if v > n_clad]
        ob = [round(v, 5) for v in _orphans(b, a) if v > n_clad]
        print(f"  {res:6.1f} {_match_count(a, b):8d} "
              f"{int((a > n_clad).sum()):12d} {int((b > n_clad).sum()):12d}   "
              f"{oa} | {ob}", flush=True)
    print("\n  ^ n_agree rising toward N_MODES as res increases => discretization.", flush=True)
    print("    n_agree flat, with the SAME orphan n_eff persisting => structural.", flush=True)

    # ---- Full grid, for completeness. Off-diagonal cells are expected to
    # ---- look bad and must NOT be read as engine disagreement (see above).
    print("\n===== full grid (off-diagonal cells are meaningless by "
          "construction -- kept only to show that) =====", flush=True)
    print("  If n_agree RISES with resolution, the mode-8+ divergence is a")
    print("  discretization artifact and the remedy is a per-engine converged")
    print("  res. If it PLATEAUS below 16, a finer mesh cannot reach it and")
    print("  something structural is wrong.\n", flush=True)
    print(f"  {'meep_res':>9} {'meow_res':>9} {'n_agree':>8} {'guided_meep':>12} "
          f"{'guided_meow':>12}   orphans (meep | meow)", flush=True)
    for mr in sorted(meep):
        for wr in sorted(meow):
            a, b = meep[mr], meow[wr]
            oa = [round(v, 5) for v in _orphans(a, b) if v > n_clad]
            ob = [round(v, 5) for v in _orphans(b, a) if v > n_clad]
            print(f"  {mr:9.1f} {wr:9.1f} {_match_count(a, b):8d} "
                  f"{int((a > n_clad).sum()):12d} {int((b > n_clad).sum()):12d}   "
                  f"{oa} | {ob}", flush=True)

    print("\n  (orphans listed are GUIDED only -- radiation modes have no "
          "engine-independent identity and are expected to differ.)", flush=True)


if __name__ == "__main__":
    main()
