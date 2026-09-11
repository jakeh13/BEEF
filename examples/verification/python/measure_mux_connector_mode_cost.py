"""Modal multiplexer: how much does the connector EME mode solve actually cost,
and does throwing MPI ranks at it help?

WHY THIS BLOCKS THE DEVICE'S NEXT STEP. The hybrid chain has only ever run at
N_WL=2; the full-device FDTD reference has 11. Every downstream comparison
needs them at matched wavelengths, and going 2 -> 11 multiplies exactly one
thing: the EME mode solves, which are per-wavelength (FDTD launch count is NOT
-- each broadband run captures the whole spectrum by DFT in one shot).

The connector cross-section is the suspect. `modal_multiplexer_hybrid.py`'s
CONN_MODES comment records a single-rank login-node measurement of **5050 s
(~84 min) for ONE wavelength** on this cross-section -- three separated,
tightly-pitched silicon islands with no exploitable mirror symmetry, a much
harder MPB target than an MMI slab. `run_modal_multiplexer_hybrid_SMOKE.batch`
then says in as many words that "how much real MPI parallelism brings that down
is completely unmeasured".

That one unmeasured number decides whether N_WL=11 is a 6-hour job or a 30-hour
one, and whether the batch script's 4-node/96-rank shape is right or wasteful.

THE SPECIFIC REASON TO DOUBT RANK SCALING: this repo has already measured that
MPB is rank-REDUNDANT -- every rank solves the same eigenproblem -- while FDTD
time-stepping scales ~98%. The hybrid batch allocates 96 ranks, which is right
for the FDTD sections and possibly buys nothing at all for the mode solves. If
so, 11 wavelengths x 84 min = 15.4 h for the connector ALONE, before either MMI
body, regardless of node count.

WHAT THIS MEASURES. The same section_modes() call the cascade makes, on the
same cross-section, at one wavelength -- timed. Run it at several --ntasks
values and compare: flat wall-clock across rank counts means rank-redundant
(add wavelengths by running them as independent array tasks, not by adding
ranks); a real speedup means the batch shape is already right.

Also times an MMI-body cross-section in the same job as a control, so a slow
connector can be attributed to the cross-section rather than to the machine.

Run (see batch/measure_mux_connector_cost.batch):
    srun --mpi=pmi2 python3 verification/python/measure_mux_connector_mode_cost.py
"""
import builtins
import json
import os
import sys
import time

import numpy as np


def print(*a, **k):                      # noqa: A001 - rank guard
    try:
        import meep as _mp
        if not _mp.am_master():
            return
    except Exception:
        pass
    builtins.print(*a, **k)


def _is_master() -> bool:
    try:
        import meep as _mp
        return bool(_mp.am_master())
    except Exception:
        return True


def _nranks() -> int:
    for v in ("SLURM_NTASKS", "OMPI_COMM_WORLD_SIZE", "PMI_SIZE"):
        if os.environ.get(v):
            return int(os.environ[v])
    try:
        import meep as _mp
        return int(_mp.count_processors())
    except Exception:
        return 1


def main() -> int:
    _HERE = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(_HERE, "..", "..", "modal_multiplexer", "python"))

    # DO NOT `import modal_multiplexer_hybrid` -- it has NO `if __name__ ==
    # "__main__"` guard, so importing it runs the entire cascade (checked:
    # `hybrid.solve()` sits at module scope, line ~497). Instead exec only the
    # construction prologue, cut at the "4. Geometry plot" banner: everything
    # above it is constants, geometry and solver construction with no physics
    # and no file writes; everything at or below it plots and then solves.
    #
    # This keeps fidelity to the real device -- same geometry object, same
    # EMESolver instance, same spans -- without re-typing constants that could
    # silently drift from the script they are meant to measure.
    src_path = os.path.join(_HERE, "..", "..", "modal_multiplexer", "python",
                            "modal_multiplexer_hybrid.py")
    src = open(src_path).read()
    marker = "# 3. Geometry plot"
    # was "# 4." until 2026-09-09 -- modal_multiplexer_fdtd.py's own section
    # numbering shifted (a prior edit folded/renumbered its construction
    # steps) sometime before this comment was added; this script's own
    # guard below would have caught the mismatch safely (hard SystemExit,
    # not silent misbehavior) the next time it ran, but there is no reason
    # to leave it stale now that it is noticed.
    cut = src.find(marker)
    if cut < 0:
        raise SystemExit(f"cannot find the construction/solve boundary "
                         f"({marker!r}) in {src_path} -- the script was "
                         f"restructured; re-check where physics begins before "
                         f"trusting this cut.")
    # back up to the start of the banner comment block above the marker
    cut = src.rfind("# ---", 0, cut)
    os.environ.setdefault("MUX_HYBRID_N_WL", "2")

    class _M:
        pass
    M = _M()
    ns = {"__name__": "_mux_prologue", "__file__": src_path}
    exec(compile(src[:cut], src_path, "exec"), ns)
    for k, v in ns.items():
        setattr(M, k, v)
    print(f"  (exec'd {src[:cut].count(chr(10))} lines of construction "
          f"prologue, physics excluded)")

    wl = float(os.environ.get("MUX_COST_WL", "1.31"))
    nrank = _nranks()

    conn_x0 = M.conn_start + M.CONN_STUB
    conn_x1 = M.conn_end - M.CONN_STUB
    body_x0 = M.body1_start + M.BODY_STUB
    body_x1 = M.body1_end - M.BODY_STUB

    print(f"=== connector mode-solve cost, {nrank} rank(s), wl={wl} um ===")
    print(f"  connector span x [{conn_x0:+.4f}, {conn_x1:+.4f}] "
          f"({conn_x1 - conn_x0:.4f} um), CONN_MODES={M.CONN_MODES}")
    print(f"  MMI-1 body span x [{body_x0:+.4f}, {body_x1:+.4f}] "
          f"({body_x1 - body_x0:.4f} um), N_MODES={M.N_MODES}")
    print(f"  body_solver res={M.body_solver.res}, "
          f"oxide {M.body_solver.oxide_width}x{M.body_solver.oxide_height}")

    out = dict(n_ranks=nrank, wl=wl, conn_modes=M.CONN_MODES,
               n_modes=M.N_MODES, res=float(M.body_solver.res))

    for label, x0, x1, nm in (("connector", conn_x0, conn_x1, M.CONN_MODES),
                              ("mmi1_body", body_x0, body_x1, M.N_MODES)):
        print(f"\n  solving {label} ({nm} modes)...")
        t0 = time.time()
        try:
            # section_modes returns (modes, mesh) -- a 2-TUPLE, not a mode
            # list. Unpacking it wrong reports "2 modes returned" and nan
            # n_eff for every cross-section regardless of what was requested,
            # which looks like a truncating solver and is not. The timing is
            # unaffected either way; the mode count and n_eff are not.
            modes, mesh = M.body_solver.section_modes(M.mux, wl, x0, x1, nm)
            dt = time.time() - t0
            neff = [complex(getattr(m, "neff", np.nan)).real for m in modes]
            print(f"    {dt:8.1f} s  ({dt/60.0:6.2f} min)  "
                  f"{len(modes)} modes returned")
            print("    n_eff: " + ", ".join(f"{n:.5f}" for n in neff))
            out[label] = dict(seconds=dt, n_returned=len(modes), neff=neff)
        except Exception as e:
            dt = time.time() - t0
            print(f"    FAILED after {dt:.1f} s: {type(e).__name__}: {e}")
            out[label] = dict(seconds=dt, error=f"{type(e).__name__}: {e}")

    c = out.get("connector", {}).get("seconds")
    b = out.get("mmi1_body", {}).get("seconds")
    if c and b:
        print(f"\n  connector / body cost ratio: {c/b:.2f}x")
    if c:
        print(f"  projected connector cost for 11 wavelengths: "
              f"{11*c/3600.0:.1f} h at this rank count")
        print(f"  (the cascade solves the connector once per wavelength; "
              f"both MMI bodies add their own on top)")

    if _is_master():
        root = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
        d = os.path.join(root, "example_output", "modal_multiplexer",
                         "connector_cost")
        os.makedirs(d, exist_ok=True)
        p = os.path.join(d, f"cost_nranks{nrank}_wl{wl}.json")
        with open(p, "w") as f:
            json.dump(out, f, indent=2)
        builtins.print(f"\nwrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
