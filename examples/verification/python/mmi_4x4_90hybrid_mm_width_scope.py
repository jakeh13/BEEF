"""4x4 MMI 90-degree hybrid -- CHEAP mm_width scoping (Stage 1 of 2).

validation_markdown/mmi_4x4_90_hybrid_validation.md, "Open decisions"
item 5's own follow-up, after the joint (offset, ΔL) grid closed off
`waveguide_offset` as a lever (offset=W/4 is close to optimal at the
current W=6.0, independent of ΔL) -- the next open candidate is
`mm_width` itself. Two things distinguish this from the offset/ΔL work:

1. Changing mm_width has NO free cache reuse -- a new width means a
   genuinely different body cross-section, so its own N_MODES_BODY=16 EME
   solve (~77-90 min at production settings, per body_sweep.py's own
   precedent) cannot be avoided. That is real, per-width cost.
2. The user flagged a real coupling this doc had not yet accounted for:
   the OPTIMAL waveguide_offset changes with mm_width too (Step 2's own
   derivation, `offset = W/4`, is itself a function of W) -- so a real
   mm_width study is a genuine (width, offset, ΔL) joint search, not a
   1-D width sweep at the current offset=1.5/W=6.0 ratio.

Given both, per-width full evaluation (body solve + baseline 4-solve +
several offset points) is expensive enough (~15-20h/width, no cache reuse)
that picking candidate widths blindly would be wasteful. THIS script is
the cheap Stage 1 instead: extends Step 1c's own bare 2-mode L_pi-vs-width
measurement (`mmi_4x4_90hybrid_cross_section.py`, EMESolver.compute_modes()
on a bare UniformCrossSection strip -- no FDTD, no time-stepping, "login-
node smoke run" per that script's own docstring) to a WIDER width range,
narrower than the 6.0-16.0 um Step 1c originally covered (6.0 was only
ever the narrowest width Step 1c actually tried, chosen for compactness
by elimination, not because anything narrower was ruled out) -- plus a
cheap guided-mode-COUNT probe per width (same bare-strip pipeline, more
modes requested, at the real production-like domain oxide_width=width+3.0,
NOT the wide diagnostic-only domain the original N=16 body ceiling probes
found to spuriously stall -- see the validation doc's "Body/taper
mode-count convergence investigation").

Geometric floor: `offset ~ W/4` must clear `taper_width=0.8` for the
ports not to overlap (Step 2's own clearance check), i.e. W > 3.2 needed
at this taper_width. Swept range starts at 3.5, just inside that floor.

Run (login node, cheap):
    python mmi_4x4_90hybrid_mm_width_scope.py
"""
from __future__ import annotations

import time
import numpy as np

from beef import EMESolver, SOI_220
from beef.geometry import UniformCrossSection
from beef.materials.stack import sellmeier_index
from beef._utils import extrude_polygons
import meow as mw

WL = 1.31
STACK = SOI_220
RES = 40.0
N_CLAD = sellmeier_index("sio2", WL)
TAPER_WIDTH = 0.8   # Step 2's placeholder -- clearance floor for W/4 > TAPER_WIDTH


def solve_modes(strips, num_modes, sim_w, sim_h, res=RES, wl=WL):
    eme = EMESolver(
        num_modes=[num_modes], num_cells=[1], lengths=[1.0],
        oxide_width=sim_w, oxide_height=sim_h,
        res=float(res), material_stack=STACK, mode_solver="meep",
    )
    cs_geom = UniformCrossSection(strips, length=1.0)
    gds = cs_geom.section_gds(-0.5, 0.5, eme.oxide_width)
    cells = mw.create_cells(
        structures=extrude_polygons(gds, eme.material_stack.extrusion_rules),
        mesh=eme._build_mesh(), Ls=[1.0],
    )
    cs = mw.CrossSection.from_cell(cell=cells[0], env=mw.Environment(wl=wl, T=eme.temperature))
    modes = eme.compute_modes(cs, num_modes)
    return np.array([float(np.real(m.neff)) for m in modes])


def guided(neffs):
    return [n for n in neffs if n > N_CLAD]


def lpi_scope(mm_widths, sim_w=20.0, sim_h=2.0):
    print(f"\n[L_pi scope] wl={WL} um, SOI_220, res={RES}, n_clad={N_CLAD:.6f}")
    print(f"{'mm_width':>9}  {'offset=W/4':>10}  {'n_eff[0]':>10}  {'n_eff[1]':>10}  "
          f"{'L_pi (um)':>10}  {'L_pi/W^2':>10}  ({'time':>6})")
    results = {}
    for w in mm_widths:
        t0 = time.time()
        neffs = solve_modes([(0.0, float(w))], num_modes=2,
                             sim_w=max(sim_w, w + 6.0), sim_h=sim_h)
        dbeta = 2 * np.pi / WL * (neffs[0] - neffs[1])
        l_pi = np.pi / dbeta
        results[float(w)] = l_pi
        print(f"{w:9.2f}  {w/4.0:10.3f}  {neffs[0]:10.6f}  {neffs[1]:10.6f}  "
              f"{l_pi:10.4f}  {l_pi/w**2:10.5f}  ({time.time()-t0:.1f}s)")
    return results


def mode_count_scope(mm_widths, n_request=24, sim_h=2.0):
    """Cheap guided-mode COUNT at the REAL production-like domain
    (oxide_width=W+3.0 -- confirmed fast, not the wide diagnostic domain
    that stalled for the N=16 ceiling probes at W=6.0)."""
    print(f"\n[mode-count scope] wl={WL} um, SOI_220, res={RES}, "
          f"n_request={n_request}, oxide_width=W+3.0")
    print(f"{'mm_width':>9}  {'guided':>6}  {'n_eff[0]':>10}  {'n_eff[-1] (guided)':>18}  ({'time':>6})")
    for w in mm_widths:
        t0 = time.time()
        neffs = solve_modes([(0.0, float(w))], num_modes=n_request,
                             sim_w=w + 3.0, sim_h=sim_h)
        g = guided(neffs)
        last = g[-1] if g else float("nan")
        print(f"{w:9.2f}  {len(g):6d}  {neffs[0]:10.6f}  {last:18.6f}  ({time.time()-t0:.1f}s)")


def main():
    # Full candidate range: [3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0]
    # (narrower than Step 1c's original 6.0-16.0 -- only the geometric
    # clearance floor, W>3.2, rules anything out -- plus the existing
    # baseline 6.0 and Step 1c's own already-known 8.0/10.0).
    #
    # CONTINUATION run (job 11840371 TIMED OUT at -t 480 with the
    # mode-count stage only 3/9 through -- 3.5/4.0/4.5 already measured:
    # 23/24/24 guided modes, ~1.3-2h each). L_pi scope already fully
    # completed and confirmed clean (L_pi/W^2 = 3.023-3.026 across the
    # WHOLE 3.5-10.0 range, same as the original 6-16 measurement) -- not
    # re-run here, no new information to gain. This script has no
    # save_dir/caching of its own (a quick diagnostic, not a production
    # solver call), so a plain resubmit would otherwise redo 3.5/4.0/4.5
    # from scratch -- skip them explicitly instead.
    remaining_widths = [5.0, 5.5, 6.0, 7.0, 8.0, 10.0]
    mode_count_scope(remaining_widths, n_request=24)


if __name__ == "__main__":
    main()
