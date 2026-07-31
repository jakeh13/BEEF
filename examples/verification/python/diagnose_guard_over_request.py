"""Does `guard_duplicate_modes=True` change the EME answer, or only filter it?
(validation_markdown/fdtd_vs_EME_mode_validation.md item 13.)

The problem with "guard on vs off"
-----------------------------------
`EMESolver.compute_modes()` (`beef/solvers/eme.py`) does NOT simply filter the
modes it would otherwise have returned. With the guard on it calls

    self._solve_modes_raw(cross_section, num_modes + _DUP_MODE_MAX_EXTRA)

-- one `compute_modes_meep()` call for the lowest **24** bands when 16 were
asked for -- then keeps the first 16 distinct ones. With the guard off it asks
for 16. So the two paths hand MPB *different eigenproblems*, and flipping the
flag changes two things at once:

    (1) how many bands MPB is asked to converge simultaneously, and
    (2) whether duplicates are then filtered out.

That is why item 13 says "guard on vs off is not a clean A/B". Item 0 checked
inertness at **3 modes, res 12** and got byte-identical n_eff -- but that is the
regime where (1) is least likely to bite. At the production end (16 modes
becoming 24 on a fine mesh) neither the inertness nor the cost has ever been
measured, and job 11575516's 3h36m runtime was dominated by exactly this.

Three arms, not two
--------------------
Two arms cannot separate (1) from (2). Three can:

    A. guard OFF, 16 modes  -- the unguarded path; baseline.
    B. guard ON,  16 modes  -- the production default (solves 24, filters to 16).
    C. guard OFF, 24 modes  -- the SAME eigenproblem B hands MPB, unfiltered.

Then:

    C vs A  ->  does over-requesting alone move the low bands?  This is the
                artifact question item 13 actually cares about. If the first 16
                n_eff of C differ from A, then every guarded EME number in the
                repo carries an over-request artifact, independent of whether
                any duplicate was ever filtered.
    B vs C  ->  what does the filter itself do? B should be C's list with
                duplicates dropped, so B tells us which of C's 24 the guard
                rejected and what it reached further up the spectrum to get.

Wall time is recorded per arm because the cost half of item 13 is a real
decision: if the low bands do NOT move, `_DUP_MODE_MAX_EXTRA` can be cut on
cost grounds alone, and C vs A is exactly the measurement that licenses that.

Cross-section
-------------
Modal multiplexer interface 0 -- the same cross-section as
`diagnose_dup_band_grid.py` (item 11) and `diagnose_kpoint_ladder_v2.py`
(item 12), and the one where meow's modes 10/11 actually duplicated at res 40
(item 0). Run at **res 60** per item 13: above the shipped 40, and the
resolution job 11575516 found meow converged at (16/16 guided, plateaued
between 60 and 80) -- so a difference seen here is not a not-yet-converged mesh
talking.

How to read the output
----------------------
The headline is the `C vs A` max |dn_eff| over the first 16 bands:

    ~0 (below print precision)  ->  over-requesting is inert at production
        scale. The guard is a pure filter, every guarded number stands as
        measured, and `_DUP_MODE_MAX_EXTRA` becomes a free cost knob.
    non-zero                    ->  the over-request perturbs the bands it is
        supposed to be innocently reaching past. Item 13's warning is real and
        every guarded EME result needs re-reading in that light.
"""
from __future__ import annotations

import time

import numpy as np

from beef import MMI, MMIPort, WaveguideRun, CompoundGeometry, EMESolver, SOI_160
from beef.solvers.eme import _DUP_MODE_MAX_EXTRA

WL = 1.31
N_MODES = 16            # what the device requests
RES = 60.0              # item 13's prescription -- above the shipped 40, and
                        # where 11575516 found meow converged (16/16 guided)
BODY_STUB = 0.75
SIM_W, SIM_H = 5.0, 1.5
HALF_SPAN = 0.05
N_OVER = N_MODES + _DUP_MODE_MAX_EXTRA   # 24 -- what the guard solves for

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


def _run_arm(name, geometry, x_loc, guard, n_request):
    """One arm: build a fresh EMESolver, solve, return (neff array, seconds).

    A fresh solver per arm on purpose -- `compute_modes()` is the chokepoint
    under test, and reusing one risks a cached mesh or mode list making two
    arms share work they should not.
    """
    print(f"\n===== {name} =====", flush=True)
    print(f"  guard_duplicate_modes={guard}, requesting {n_request} modes",
          flush=True)
    eme = EMESolver(num_modes=[n_request], num_cells=[1], lengths=[1.0],
                    oxide_width=SIM_W, oxide_height=SIM_H, res=RES,
                    material_stack=SOI_160, mode_solver="meep",
                    guard_duplicate_modes=guard)
    t0 = time.perf_counter()
    try:
        modes, _ = eme.section_modes(geometry, WL, x_loc - HALF_SPAN,
                                     x_loc + HALF_SPAN, n_request)
    except Exception as exc:                      # guard exhaustion raises
        dt = time.perf_counter() - t0
        print(f"  FAILED after {dt:.1f} s: {type(exc).__name__}: {exc}",
              flush=True)
        return None, dt
    dt = time.perf_counter() - t0
    neff = np.array([float(np.real(m.neff)) for m in modes])
    print(f"  returned {len(neff)} modes in {dt:.1f} s", flush=True)
    print(f"  n_eff: {[round(v, 6) for v in neff]}", flush=True)
    return neff, dt


def main():
    mmi1, geometry = _build()
    x_loc = -(mmi1.mm_length / 2 - BODY_STUB)
    n_clad = SOI_160.sio2_index(WL)
    print(f"cross-section: mmi1 local x = {x_loc:.4f} um (interface 0), "
          f"wl = {WL} um, res {RES}", flush=True)
    print(f"N_MODES = {N_MODES}, _DUP_MODE_MAX_EXTRA = {_DUP_MODE_MAX_EXTRA} "
          f"-> guard solves {N_OVER}", flush=True)
    print(f"cladding index = {n_clad:.6f}", flush=True)

    neff_a, t_a = _run_arm("A. guard OFF, 16 modes (unguarded baseline)",
                           geometry, x_loc, guard=False, n_request=N_MODES)
    neff_b, t_b = _run_arm("B. guard ON, 16 modes (production default)",
                           geometry, x_loc, guard=True, n_request=N_MODES)
    neff_c, t_c = _run_arm(f"C. guard OFF, {N_OVER} modes (B's eigenproblem, unfiltered)",
                           geometry, x_loc, guard=False, n_request=N_OVER)

    # ---- C vs A: does the over-request alone move the low bands? ----------
    print("\n===== C vs A -- does over-requesting move the low bands? =====",
          flush=True)
    if neff_a is None or neff_c is None:
        print("  SKIPPED -- an arm failed above.", flush=True)
    else:
        k = min(len(neff_a), N_MODES, len(neff_c))
        dev = np.abs(neff_c[:k] - neff_a[:k])
        print(f"  comparing the first {k} bands, positionally", flush=True)
        print(f"  max |dn_eff| = {np.max(dev):.3e}   "
              f"(at band {int(np.argmax(dev)) + 1})", flush=True)
        print(f"  per-band |dn_eff|: {[f'{v:.2e}' for v in dev]}", flush=True)
        if np.max(dev) < 1e-9:
            print("  -> INERT: asking for 24 instead of 16 does not move the "
                  "low bands.\n     The guard is a pure filter at production "
                  "scale, and _DUP_MODE_MAX_EXTRA\n     is a free cost knob.",
                  flush=True)
        else:
            print("  -> NOT INERT: the over-request perturbs the bands it is "
                  "meant to reach\n     past. Item 13's warning is real -- "
                  "every guarded EME number carries\n     this artifact.",
                  flush=True)

    # ---- B vs C: what did the filter actually drop? -----------------------
    print("\n===== B vs C -- what does the filter drop? =====", flush=True)
    if neff_b is None or neff_c is None:
        print("  SKIPPED -- an arm failed above.", flush=True)
    else:
        # B should be C's list minus duplicates. Match each B mode to its
        # nearest C mode so a dropped band shows up as a positional shift.
        idx = [int(np.argmin(np.abs(neff_c - v))) for v in neff_b]
        print(f"  B's {len(neff_b)} modes map to C indices (0-based): {idx}",
              flush=True)
        resid = [abs(neff_c[j] - v) for j, v in zip(idx, neff_b)]
        print(f"  worst match residual: {max(resid):.3e}", flush=True)
        skipped = sorted(set(range(max(idx) + 1)) - set(idx)) if idx else []
        print(f"  C indices B skipped over: {skipped if skipped else 'none'}",
              flush=True)
        if skipped:
            print(f"  skipped n_eff: {[round(float(neff_c[j]), 6) for j in skipped]}",
                  flush=True)
            print("  -> the guard rejected these as duplicates and reached "
                  "further up the\n     spectrum. Cross-check each against an "
                  "exact n_eff tie in C's list.", flush=True)
        else:
            print("  -> the guard dropped nothing here: B is C's first "
                  f"{len(neff_b)} modes.\n     No duplicate fired at this "
                  "resolution.", flush=True)

    # ---- cost ------------------------------------------------------------
    print("\n===== cost =====", flush=True)
    print(f"  A (16, unguarded): {t_a:7.1f} s", flush=True)
    print(f"  B (16, guarded  ): {t_b:7.1f} s", flush=True)
    print(f"  C ({N_OVER}, unguarded): {t_c:7.1f} s", flush=True)
    if t_a > 0:
        print(f"  guard overhead B/A = {t_b / t_a:.2f}x", flush=True)
        print(f"  over-request cost C/A = {t_c / t_a:.2f}x   "
              "<- what _DUP_MODE_MAX_EXTRA buys/costs", flush=True)
    print("\n  If C vs A came back INERT, this ratio is the whole argument for "
          "tuning\n  _DUP_MODE_MAX_EXTRA down: it is pure cost with no effect "
          "on the answer.", flush=True)


if __name__ == "__main__":
    main()
