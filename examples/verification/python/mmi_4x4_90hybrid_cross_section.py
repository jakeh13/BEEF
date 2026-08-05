"""4x4 MMI 90-degree hybrid -- Step 1 cross-section design (SOI_220, O-band).

validation_markdown/mmi_4x4_90_hybrid_validation.md Step 1: settle the access
waveguide width (1a), check the vertical/lateral box size sz (1b), and sweep
the MMI-body self-imaging length L_pi vs body width (1c) -- all at O-band
(wl=1.31 um) on SOI_220, a (platform, band) pair no existing device combines
(every SOI_220 example in this repo runs C-band; every O-band example runs
SOI_160 -- see the validation doc's A2). Do NOT reuse mmi_2x2_hybrid.py's
waveguide_width=0.45 here -- that number was picked for SOI_220 at C-band and
is unvalidated even there; O-band's shorter wavelength confines more tightly,
which is exactly the direction that can turn a single-mode C-band width
multi-mode at O-band (A2).

Three cheap MPB eigensolves via EMESolver.compute_modes() on bare
UniformCrossSection strips -- no FDTD, no time-stepping. Same pipeline
verify_mmi_2x2_gap2_dbeta.py already uses for its own body-mode measurement
(that script's own docstring: "no time-stepping ... login-node smoke run").

TE vs TM -- a dead end, recorded rather than silently dropped. The paper's
own device excites only "the TE-polarized fundamental mode of the input
waveguides" (Sec. IV), and a first pass here found 2-3 "guided" modes even at
width=0.30 um at O-band, which looked like a TE/TM conflation (these narrow
near-square access waveguides, width and height both ~0.2-0.4 um, are exactly
the aspect ratio where TE0/TM0 sit close in n_eff). Tried filtering guided
modes by meow's Mode.te_fraction (0.0=pure TM, 1.0=pure TE per meow's own
docstring) -- and then, before trusting it, sanity-checked it against the
mmi_2x2_hybrid.py cross-section (width=0.45, C-band, an access waveguide this
whole repo already ships and documents as running its own "TE-polarized
fundamental mode"). Result: te_fraction=0.000199 for THAT mode too -- i.e.
meow's te_fraction calls this repo's own established, working "TE" mode
"TM". Whatever the actual cause (an axis/definition mismatch in meow's
Ex/Hx convention, or genuine hybridization this paper's own simplified 2D
scalar TE/TM split (Eq. 2, an effective-index-method reduction this
codebase's full 3-D vector solve does NOT make) doesn't capture), the metric
does not match this repo's own working convention and building on it would
be exactly the unearned precision this project's validation standard exists
to catch. NOT used below. No device in this repo has ever filtered guided-
mode counts by polarization -- every one (2x2, 1x2, 1x4, the coupler, the
mux) just counts n_eff > n_clad bands and relies on exact geometric symmetry
to keep TE/TM decoupled regardless of what else is in the interface basis,
which applies equally here (a plain rib/strip cross-section, no lateral or
vertical asymmetry). Reverted to that same plain convention.

Run:
    python examples/verification/python/mmi_4x4_90hybrid_cross_section.py
"""
from __future__ import annotations

import time

import numpy as np
import meow as mw

from beef import EMESolver, SOI_220
from beef.geometry import UniformCrossSection
from beef.materials.stack import sellmeier_index
from beef._utils import extrude_polygons

WL = 1.31               # O-band centre -- matches this repo's existing O-band
                         # devices (1x2 MMI, modal multiplexer), not re-derived
STACK = SOI_220
RES = 40.0               # cheap first-pass resolution -- every device's own
                         # first cross-section sweep in this repo uses this;
                         # final resolution is a separate, later check-3 question
N_CLAD = sellmeier_index("sio2", WL)


def solve_modes(strips, num_modes, sim_w, sim_h, res=RES, wl=WL):
    """Bare-strip mode solve via EMESolver.compute_modes() -- same pipeline
    verify_mmi_2x2_gap2_dbeta.py uses for the MMI body, generalized to an
    arbitrary UniformCrossSection strip list. Returns n_eff, highest first
    (meow's own ordering)."""
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


# ---------------------------------------------------------------------------
# 1a. Single-mode access-waveguide width, O-band, SOI_220.
# ---------------------------------------------------------------------------
def step_1a(sim_w=5.0, sim_h=2.0):
    print(f"\n[1a] single-mode width sweep -- wl={WL} um, SOI_220, res={RES}, "
          f"n_clad(SiO2)={N_CLAD:.6f}")
    print(f"{'width':>7}  {'guided':>6}  n_eff (highest -> lowest)")
    # Starts narrower than the original 0.30 -- that already showed 2 guided
    # modes, so the single-mode point (if any, in-range) has to be below it.
    widths = np.round(np.arange(0.10, 0.71, 0.05), 2)
    single_mode = []
    for w in widths:
        t0 = time.time()
        neffs = solve_modes([(0.0, float(w))], num_modes=4, sim_w=sim_w, sim_h=sim_h)
        g = guided(neffs)
        if len(g) == 1:
            single_mode.append(float(w))
        print(f"{w:7.2f}  {len(g):6d}  {['%.6f' % n for n in g]}  ({time.time()-t0:.1f}s)")
    print(f"[1a] single-mode window: {single_mode}")
    return single_mode


# ---------------------------------------------------------------------------
# 1b. Vertical/lateral box convergence (sz / sim_h) at the settled width.
# ---------------------------------------------------------------------------
def step_1b(width, sim_w=5.0):
    print(f"\n[1b] box-size (sz) convergence -- width={width} um, wl={WL} um")
    print(f"{'sz':>6}  {'n_eff[0]':>10}  {'step':>10}")
    szs = [1.0, 1.5, 2.0, 2.5, 3.0]
    prev = None
    for sz in szs:
        t0 = time.time()
        neffs = solve_modes([(0.0, float(width))], num_modes=1, sim_w=sim_w, sim_h=sz)
        step = "" if prev is None else f"{neffs[0]-prev:+.6f}"
        print(f"{sz:6.1f}  {neffs[0]:10.6f}  {step:>10}   ({time.time()-t0:.1f}s)")
        prev = neffs[0]


# ---------------------------------------------------------------------------
# 1c. L_pi vs candidate MMI body width (paper Eq. 9: L_pi = pi/(beta0-beta1)).
# ---------------------------------------------------------------------------
def step_1c(sim_w=20.0, sim_h=2.0):
    print(f"\n[1c] L_pi vs body width -- wl={WL} um, SOI_220, res={RES}")
    print(f"{'mm_width':>9}  {'n_eff[0]':>10}  {'n_eff[1]':>10}  "
          f"{'dbeta (rad/um)':>14}  {'L_pi (um)':>10}")
    mm_widths = np.arange(6.0, 16.01, 2.0)
    results = {}
    for w in mm_widths:
        t0 = time.time()
        neffs = solve_modes([(0.0, float(w))], num_modes=2,
                             sim_w=max(sim_w, w + 6.0), sim_h=sim_h)
        dbeta = 2 * np.pi / WL * (neffs[0] - neffs[1])
        l_pi = np.pi / dbeta
        results[float(w)] = l_pi
        print(f"{w:9.1f}  {neffs[0]:10.6f}  {neffs[1]:10.6f}  {dbeta:14.5f}  "
              f"{l_pi:10.4f}   ({time.time()-t0:.1f}s)")
    return results


def main():
    single_mode = step_1a()
    # Use the widest single-mode width as the trial access width for 1b --
    # the design will want the widest single-mode point for lowest loss,
    # same reasoning mmi_1x4_sin300_validation.md step 4b used ("report the
    # widest width that stays single-mode"). Re-decided for real in Step 2.
    trial_width = max(single_mode) if single_mode else 0.45
    print(f"\n[main] using trial access width = {trial_width} um for step 1b")
    step_1b(trial_width)
    step_1c()


if __name__ == "__main__":
    main()
