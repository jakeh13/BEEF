"""TEMPORARY -- one-off diagnostic, delete once gap 7's TE1/1.26um item in
modal_multiplexer_validation.md is resolved.

Checks HybridSolver's interface mode-correspondence machinery (gap 1's
mechanism) at EACH of modal_multiplexer's 6 interfaces, at wl=1.26 AND
wl=1.36 SEPARATELY (HybridSolver's own check_mode_correspondence only ever
checks at mean(wavelengths)=1.31, never at the actual swept points) --
looking for whether the meep<->meow mode-index correspondence REORDERS
between the two wavelengths, which would explain the TE1 top/bot swap found
in modal_multiplexer_validation.md gap 7's 2026-08-09 hybrid-vs-FDTD
comparison (real, measured discrepancy isolated to TE1 @ 1.26 um only).

Reuses the real modal_multiplexer_hybrid.py's own construction verbatim
(exec'd, truncated right after `hybrid = HybridSolver(...)`) so this checks
the ACTUAL boundaries/interface_modes/solvers the production script uses,
not a stand-in.

section_mode_fields()/section_modes() are "seconds"-scale MPB/meow
eigensolves (no time-stepping) -- SLURM, not login-node, per this repo's
own discipline, but cheap and single-rank.
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

# Rebuild with check_mode_correspondence on, reusing the exact same
# solvers/boundaries/interface_modes -- constructing a second HybridSolver
# instance doesn't re-solve anything (solving only happens in .solve()/the
# private _check_interface_correspondence() call below).
from beef.solvers.hybrid import HybridSolver

checked = HybridSolver(
    solvers=hybrid.solvers,
    boundaries=hybrid.boundaries,
    interface_modes=hybrid.interface_modes,
    check_mode_correspondence=True,
    mode_correspondence_pad=0.1,
    mode_correspondence_action="warn",  # report every interface/wavelength
    # in one pass instead of stopping at the first raise -- 2026-08-09
    # rerun testing whether N_MODES=16 (raised from 14) closes the
    # basis-incompleteness the "raise" pass found at interface 0/wl=1.26.
)

kwargs_per_section_top = checked._kwargs_per_section({6: {"output_idx": 0}})

labels = [
    "FDTD-in <-> EME-MMI1",
    "EME-MMI1 <-> FDTD-taper-out",
    "FDTD-taper-out <-> EME-conn",
    "EME-conn <-> FDTD-taper-in",
    "FDTD-taper-in <-> EME-MMI2",
    "EME-MMI2 <-> FDTD-out(top)",
]

for wl in (1.26, 1.31, 1.36):
    print(f"\n{'='*20} wl_center = {wl} um {'='*20}", flush=True)
    for k in range(6):
        kw = kwargs_per_section_top
        print(f"\n-- interface {k}: {labels[k]} --", flush=True)
        checked._check_interface_correspondence(mux, k, kw[k], kw[k + 1], wl)

print("\nDIAGNOSTIC DONE", flush=True)
