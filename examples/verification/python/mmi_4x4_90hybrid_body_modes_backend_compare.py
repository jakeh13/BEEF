"""4x4 MMI 90-degree hybrid -- compute body-cross-section modes via ONE
mode_solver backend and save them to disk, for a separate meep-vs-tidy3d
comparison.

Not one of validation_markdown/mmi_4x4_90_hybrid_validation.md's original
Steps 0-9 -- prompted by the body mode-count convergence investigation:
N_MODES_BODY=8 was found significantly under-converged (job 11752350, via
the "tidy3d" backend), but body_sweep.py/hybrid.py's own EME body
sub-solver is deliberately locked to "meep" (for FDTD/EME mode-basis
consistency across the hybrid cascade -- see that script's own comment),
and "meep" cannot compute N_MODES_BODY=16 in any practical time (job
11753042/11749479 both timed out, 3h/4h, zero output). Before accepting
that limitation, checking whether meep's and tidy3d's own mode bases
actually agree at a mode count BOTH can compute (N=8) -- if they do, the
body sub-solver could switch to "tidy3d" without breaking basis
consistency, unblocking real convergence for the hybrid cascade too.

Split into two single-backend jobs (this script, run twice with different
MODE_SOLVER) rather than one script computing both sequentially -- the
first attempt did that on the login node and was still running the meep
half after 39 minutes; splitting lets the two run in PARALLEL batch jobs
instead of serially eating one already-slow half after the other.

Run (submit twice, once per backend):
    MODE_SOLVER=meep python examples/verification/python/mmi_4x4_90hybrid_body_modes_backend_compare.py
    MODE_SOLVER=tidy3d python examples/verification/python/mmi_4x4_90hybrid_body_modes_backend_compare.py
"""
import os
import time

import meow as mw

from beef import EMESolver, SOI_220
from beef.geometry import UniformCrossSection
from beef._utils import extrude_polygons

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(
    _HERE, "..", "..", "..", "..", "example_output", "mmi_4x4_90hybrid", "mode_convergence", "backend_compare"))
os.makedirs(output_dir, exist_ok=True)

WL = 1.31
MM_WIDTH = 6.0
RES = 40.0
SIM_W = MM_WIDTH + float(os.environ.get("SIM_W_MARGIN", 14.0))   # default matches
              # body_modes.py's own generous-margin convention; override to 3.0 to
              # match body_sweep.py/hybrid.py/eme.py's actual PRODUCTION domain
              # (oxide_width=MM_WIDTH+3.0) -- these are NOT the same simulation
              # span, a real discrepancy caught 2026-08-09. Every meep-wall
              # bisection point measured so far (N=8/10/12/16) used the WIDE
              # (+14.0) domain, not what production actually runs.
SIM_H = float(os.environ.get("SIM_H", 2.0))   # 2.0 = Step 1b's converged sz --
              # but that check was done for the NARROW 0.35um access waveguide,
              # never re-verified for this WIDE 6um body cross-section's own
              # higher-order modes (which extend further vertically into
              # cladding than a single-mode narrow guide's do) -- override to
              # check height convergence directly, same idea as SIM_W_MARGIN.
N_MODES = int(os.environ.get("N_MODES", 8))   # default: the largest count BOTH
              # backends are known to handle (meep: job 11707835, ~55 min;
              # tidy3d: trivial by comparison). Override for a meep-ceiling
              # bisection (known good=8, known bad=16 -- jobs 11753042/11749479).

MODE_SOLVER = os.environ.get("MODE_SOLVER", "meep")
assert MODE_SOLVER in ("meep", "tidy3d"), f"unknown MODE_SOLVER={MODE_SOLVER!r}"

print(f"Computing {N_MODES} body modes via mode_solver={MODE_SOLVER!r} "
      f"(mm_width={MM_WIDTH}, wl={WL}, res={RES})...")
t0 = time.time()

eme = EMESolver(num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
                 oxide_width=SIM_W, oxide_height=SIM_H, res=RES,
                 material_stack=SOI_220, mode_solver=MODE_SOLVER)
cs_geom = UniformCrossSection([(0.0, MM_WIDTH)], length=1.0)
gds = cs_geom.section_gds(-0.5, 0.5, eme.oxide_width)
cells = mw.create_cells(structures=extrude_polygons(gds, eme.material_stack.extrusion_rules),
                         mesh=eme._build_mesh(), Ls=[1.0])
cs = mw.CrossSection.from_cell(cell=cells[0], env=mw.Environment(wl=WL, T=eme.temperature))
modes = eme.compute_modes(cs, N_MODES)

dt = time.time() - t0
print(f"  done ({dt:.1f}s)")

for i, m in enumerate(modes):
    path = os.path.join(output_dir, f"mode_{MODE_SOLVER}_n{N_MODES}_w{int(SIM_W)}_{i:02d}.json")
    m.save(path)
    print(f"  mode {i+1}: n_eff={float(m.neff.real):.6f}  -> {path}")

print(f"\nSaved {N_MODES} modes -> {output_dir}/mode_{MODE_SOLVER}_n{N_MODES}_w{int(SIM_W)}_*.json ({dt:.1f}s total)")
