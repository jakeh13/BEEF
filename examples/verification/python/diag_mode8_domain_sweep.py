"""TEMPORARY -- one-off diagnostic, delete once
mode_crossing_correspondence_repair_validation.md's Stage 0 is resolved.

Tests whether the modal multiplexer's interface-0 mode-8 correspondence
failure (row completeness 0.009 at N_MODES=14, wl=1.26 um -- see
modal_multiplexer_validation.md gap 7's 2026-08-09 entries) is sensitive to
the mode-solve's transverse domain size (sim_width/sim_height for FDTD,
oxide_width/oxide_height for EME), which has NEVER been convergence-tested
for any mode-solving context in this repo (only HCG's own full-device
reflection sim has ever had a domain_y-style sweep --
hcg_reflection_convergence.py).

If mode 8's row completeness improves with a bigger domain -> real evidence
of transverse-domain truncation as a contributing/causal factor, on top of
(or instead of) the cross-engine ladder-resolution hypothesis
(mode_crossing_correspondence_repair_validation.md's regime C).
If it stays flat across domain sizes -> rules domain truncation out, and
strengthens the regime-C (cross-engine resolution failure) hypothesis.

Scoped to interface 0, wl=1.26 um only (the one already-confirmed failure)
to keep cost bounded -- each larger domain costs more (grid points scale
with area at fixed resolution), so this is NOT the full 6-interface x
3-wavelength sweep the other diagnostic did.

Reuses the real modal_multiplexer_hybrid.py's own construction verbatim
(exec'd, truncated right after `hybrid = HybridSolver(...)`), then mutates
the already-constructed fdtd_solver/body_solver instances' domain-size
attributes directly between checks -- section_mode_fields()/section_modes()
both read self.sim_width/self.oxide_width etc. at CALL time (confirmed by
reading _build_section_mode_sim()/EMESolver's own methods before writing
this), so no need to reconstruct the solvers or geometry per domain size.
"""
import logging
import os
import sys

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

print(f"baseline: fdtd sim_width={fdtd_solver.sim_width} sim_height={fdtd_solver.sim_height}  "
      f"eme oxide_width={body_solver.oxide_width} oxide_height={body_solver.oxide_height}", flush=True)

from beef.solvers.hybrid import HybridSolver

checked = HybridSolver(
    solvers=hybrid.solvers,
    boundaries=hybrid.boundaries,
    interface_modes=hybrid.interface_modes,
    check_mode_correspondence=True,
    mode_correspondence_pad=0.1,
    mode_correspondence_action="warn",
)
kwargs_per_section = checked._kwargs_per_section({6: {"output_idx": 0}})

# (sim_width, sim_height) pairs -- baseline first, then progressively larger,
# matching hcg_reflection_convergence.py's domain_y sweep style (a few
# points, not a dense grid -- this is a cheap sanity check, not a full
# convergence campaign).
#
# ROUND 2, 2026-08-09: the 3-point sweep at n_modes=12 (5.0/7.0/9.0) found
# a clean convergent trend -- 5.0 fails at mode 8, 7.0 clean at mode 8 but
# fails at mode 11, 9.0 FULLY clean (see the validation doc's own entry).
# That used a reduced n_modes=12 to dodge a guided-mode-count crash at
# 7.0/2.0 (only 14 guided there, vs 15 at the 5.0 baseline). This run
# tests the one thing that sweep couldn't: is the SHIPPED N_MODES (15,
# not the reduced 12) also clean at the largest domain (9.0, 2.5)? If so,
# domain size is confirmed as the actual, complete fix for this device's
# real production configuration -- not just for a reduced test case.
DOMAIN_SIZES = [
    (9.0, 2.5),
]
print(f"(testing at the SHIPPED n_modes={hybrid.interface_modes[0]}, "
      f"no override this time)", flush=True)

for sim_width, sim_height in DOMAIN_SIZES:
    print(f"\n{'='*20} sim_width={sim_width}, sim_height={sim_height} {'='*20}", flush=True)
    fdtd_solver.sim_width = sim_width
    fdtd_solver.sim_height = sim_height
    body_solver.oxide_width = sim_width
    body_solver.oxide_height = sim_height
    checked._check_interface_correspondence(
        mux, 0, kwargs_per_section[0], kwargs_per_section[1], 1.26)

print("\nDIAGNOSTIC DONE", flush=True)
