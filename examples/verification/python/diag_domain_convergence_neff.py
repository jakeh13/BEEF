"""TEMPORARY -- one-off diagnostic, delete once
mode_crossing_correspondence_repair_validation.md's Stage 0 domain-size
question is fully closed.

Rigorous multi-point domain-size (sim_width/sim_height) convergence sweep
for the modal multiplexer's interface 0, at the shipped N_MODES=14 --
mirroring gap 2's own resolution-convergence methodology
(modal_multiplexer_validation.md, "modes agreeing between consecutive
resolutions... within |Delta n_eff| <= 0.005" table) but for the domain-size
axis instead of resolution, and hcg_reflection_convergence.py's own
domain_y sweep precedent (the only other domain-size convergence campaign
in this repo, applied to a full-device sim rather than a mode-solve).

Prompted by the 3-point sanity check (5.0/7.0/9.0 um) already run this
session, which found the interface-0 correspondence failure was a
domain-truncation artifact, not a real crossing/degeneracy -- but that
check only ever asked a binary pass/fail question via the expensive
cross-engine correspondence check. This script asks the cheaper, more
quantitative question directly: does EACH ENGINE'S OWN n_eff stabilize as
the domain grows, the same way gap 2 tracked resolution convergence?

Two independent, CHEAP (n_eff-only, no field sampling) sweeps:
  - meep: FDTDSolver.section_mode_neff() -- explicitly documented as
    "seconds, not minutes" (no time-stepping, no field grid).
  - meow: EMESolver.section_modes(), keeping only .neff off the returned
    Mode objects (meow has no separate n_eff-only path -- fields come from
    the same eigensolve call -- but this still skips the OTHER engine's
    cost and the correspondence machinery's own overhead, so it is much
    cheaper than the full check used in the earlier 3-point sweep).

A real multi-point convergence table (not just 3 points) is built for
both engines, at wl=1.26 um (the originally-failing point) -- if both
engines' n_eff plateau at the same domain size, that is a much stronger
convergence claim than the earlier binary pass/fail sweep, and pins down
where the plateau begins rather than confirming just one point past it.
"""
import logging
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

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s", stream=sys.stdout)

ns = {"__name__": "__main__", "__file__": SCRIPT}
code = compile(truncated, SCRIPT, "exec")
exec(code, ns)

hybrid = ns["hybrid"]
mux = ns["mux"]
fdtd_solver = ns["fdtd_solver"]
body_solver = ns["body_solver"]

N_MODES = 14  # the value this investigation recommends shipping, not the
              # current (still 15, unreverted) script constant -- tested
              # directly at whatever count we actually intend to ship.
WL = 1.26     # the originally-failing wavelength.
X_IFACE = hybrid.boundaries[1]  # interface 0: FDTD-in <-> EME-MMI1
HALF_SPAN = 0.05  # matches HybridSolver's own _CORR_HALF_SPAN
X0, X1 = X_IFACE - HALF_SPAN, X_IFACE + HALF_SPAN

print(f"interface 0 at x={X_IFACE:.4f}, query window [{X0:.4f}, {X1:.4f}], "
      f"N_MODES={N_MODES}, wl={WL}", flush=True)

# Baseline aspect ratio 1.5/5.0 = 0.3, held fixed across the sweep so this
# tests domain SIZE, not aspect ratio.
SIM_WIDTHS = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
ASPECT = 1.5 / 5.0
RES = 40.0  # matches this device's own RESOLUTION -- needed for the
            # epsilon fix below.

# FIX 2026-08-10 (job 11774688 crashed on this): sim_width * ASPECT is not
# always exactly representable in binary float -- e.g. 9.0*0.3 lands a hair
# BELOW 2.7, so meow's own `int(self.res * self.oxide_height) + 1` mesh
# construction (beef/solvers/eme.py) truncates to one pixel short of what
# was intended, breaking its dx==dy requirement (measured: dx=0.025000,
# dy=0.025234). Standard fix: nudge by a tiny epsilon before any solver
# ever sees the value, so int()-truncation lands on the intended pixel
# count regardless of which side of the float boundary the raw product
# fell on. Harmless for values that were already safely clear of a
# boundary (e.g. this device's own integer sim_widths).
EPS = 1e-6
def safe_height(width: float) -> float:
    return round(width * ASPECT, 6) + EPS

# RESUME 2026-08-10: job 11774688 already computed meep@[5,6,7,8,9] and
# meow@[5,6,7,8] before crashing on meow@9.0's mesh bug above -- seeded
# here verbatim (not recomputed) to avoid burning another ~3h re-deriving
# data already in hand. Only the still-missing points run fresh below.
meep_neff = {
    5.0: np.array([2.745767, 2.718282, 2.671886, 2.605635, 2.518087, 2.407148, 2.269841, 2.102022, 1.898461, 1.839441, 1.804071, 1.744503, 1.657514, 1.552136]),
    6.0: np.array([2.745766, 2.718282, 2.671886, 2.605635, 2.518087, 2.407148, 2.269840, 2.102020, 1.839110, 1.803711, 1.744079, 1.657434, 1.551080, 1.451166]),
    7.0: np.array([2.745766, 2.718282, 2.671886, 2.605635, 2.518087, 2.407148, 2.102020, 1.898453, 1.839050, 1.803644, 1.743996, 1.657411, 1.550771, 1.450257]),
    8.0: np.array([2.745766, 2.718282, 2.671886, 2.605635, 2.518087, 2.407148, 2.102020, 1.898453, 1.839040, 1.803632, 1.743980, 1.657405, 1.550666, 1.449702]),
    9.0: np.array([2.745766, 2.718282, 2.671886, 2.605635, 2.518087, 2.407148, 2.269840, 2.102020, 1.898453, 1.839038, 1.803630, 1.743977, 1.657403, 1.550628]),
}
meow_neff = {
    5.0: np.array([2.745768, 2.718286, 2.671895, 2.605653, 2.518118, 2.407198, 2.269920, 2.102151, 1.839426, 1.804128, 1.744689, 1.660531, 1.657937, 1.552682]),
    6.0: np.array([2.745767, 2.718286, 2.671895, 2.605652, 2.518118, 2.407197, 2.269919, 2.102149, 1.839127, 1.803803, 1.744306, 1.660027, 1.657872, 1.551759]),
    7.0: np.array([2.745767, 2.718286, 2.671895, 2.605652, 2.518118, 2.407197, 2.269919, 2.102149, 1.898678, 1.839061, 1.803728, 1.744212, 1.659877, 1.657848]),
    8.0: np.array([2.745767, 2.718286, 2.671895, 2.605652, 2.518118, 2.407197, 2.269919, 2.102149, 1.898677, 1.839051, 1.803717, 1.744197, 1.659849, 1.657841]),
}

MEEP_TODO = [w for w in SIM_WIDTHS if w not in meep_neff]
MEOW_TODO = [w for w in SIM_WIDTHS if w not in meow_neff]
print(f"resuming: meep still needs {MEEP_TODO}, meow still needs {MEOW_TODO}", flush=True)

for sim_width in SIM_WIDTHS:
    sim_height = safe_height(sim_width)
    need_meep = sim_width in MEEP_TODO
    need_meow = sim_width in MEOW_TODO
    if not (need_meep or need_meow):
        continue
    print(f"\n{'='*20} sim_width={sim_width}, sim_height={sim_height:.6f} {'='*20}", flush=True)

    if need_meep:
        fdtd_solver.sim_width = sim_width
        fdtd_solver.sim_height = sim_height
        neff_m = fdtd_solver.section_mode_neff(mux, np.array([WL]), X0, X1, N_MODES)
        neff_m = np.asarray(neff_m).reshape(-1)[:N_MODES]
        meep_neff[sim_width] = neff_m
        print(f"meep n_eff: {np.array2string(neff_m, precision=6)}", flush=True)

    if need_meow:
        body_solver.oxide_width = sim_width
        body_solver.oxide_height = sim_height
        modes, _mesh = body_solver.section_modes(mux, WL, X0, X1, N_MODES)
        neff_w = np.array([float(np.real(m.neff)) for m in modes])
        meow_neff[sim_width] = neff_w
        print(f"meow n_eff: {np.array2string(neff_w, precision=6)}", flush=True)

print(f"\n{'='*20} CONVERGENCE TABLE {'='*20}", flush=True)
print(f"{'width->':>10s} " + " -> ".join(f"{w:.0f}" for w in SIM_WIDTHS), flush=True)
for engine, table in (("meep", meep_neff), ("meow", meow_neff)):
    print(f"\n-- {engine} --", flush=True)
    for i in range(N_MODES):
        row = [table[w][i] for w in SIM_WIDTHS]
        deltas = [abs(row[j + 1] - row[j]) for j in range(len(row) - 1)]
        print(f"  mode {i:2d}: " + " -> ".join(f"{v:.6f}" for v in row) +
              "   |d|=" + " ".join(f"{d:.2e}" for d in deltas), flush=True)
    max_deltas = [max(abs(table[SIM_WIDTHS[j + 1]][i] - table[SIM_WIDTHS[j]][i])
                       for i in range(N_MODES)) for j in range(len(SIM_WIDTHS) - 1)]
    print(f"  worst |d n_eff| per step: " +
          " ".join(f"{d:.2e}" for d in max_deltas), flush=True)

print("\nDIAGNOSTIC DONE", flush=True)
