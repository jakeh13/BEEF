"""Is Step 6b's long-lived "trapped mode" physical, or a thin-PML artifact?

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-23.

THE OBSERVATION. The full-device 6b run dumps ~93% of its energy into the PML
when the pulse exits (t=340-420), then the remaining ~5.7% decays on a
multi-exponential tail whose effective rate eases from ~3.2e-3 down to
~1.0e-3 per time unit. Extrapolated, reaching decay_tol=1e-11 needs
t~19,000-20,000 == ~70 h even at 96 ranks. Two 48 h attempts (job 12227584)
were cancelled at t~6,200 because neither could finish.

THE HYPOTHESIS THIS TESTS. `dpml=0.5` um at lambda=1.31 um is only ~0.4
wavelengths of PML. A PML that thin reflects low-angle and evanescent
components instead of absorbing them, and reflected energy rattling between
the domain walls looks EXACTLY like a long-lived trapped mode with a
flattening decay tail. **`dpml` has never been convergence-tested anywhere in
this repo** -- an audit of validation_markdown/ found sweeps of domain
*clearance* (sim_width/sim_height) but none of PML *thickness*, and PML
reflection has never been named as a suspect for slow decay. `dpml` is
otherwise only ever raised 0.5->1.0 to keep enough PML pixels at low
resolution, which is a discretization concern, not an absorption one.

WHY THIS MATTERS MORE THAN LOOSENING decay_tol. The tail is so flat that
even decay_tol=1e-9 still needs t~14,800 (~50 h) and 1e-8 needs t~12,500
(~43 h) -- loosening the tolerance by two orders of magnitude barely dents
the cost. If the tail is a PML artifact, removing it is worth far more than
any tolerance change, and it would also mean every long FDTD run in this
repo has been paying for it.

METHOD. Identical to the production 6b configuration in every respect except
`dpml`, which comes from MMI4X4_PML_DPML. Run one port, print the energy
ratio every 20 time units, and compare against the dpml=0.5 baseline already
measured (job 12079104 / 12227584_0):

    t=600  ratio 5.66e-02
    t=861  ratio 1.38e-02
    t=1201 ratio 6.09e-03

A markedly LOWER ratio at matched t with thicker PML => reflection was a real
contributor and the tail is (partly) numerical. A ratio that tracks the
baseline => the trapped mode is physical and thickening PML is not the lever.

NOTE ON COST: dpml is ADDITIVE to the domain (Sy = sim_width + 2*dpml, and
likewise in z), so dpml=1.5 makes the cell 12.0 x 5.0 instead of 10.0 x 3.0
-- ~2x the voxels, so ~2x the s/step. Reaching a matched t therefore costs
about twice the baseline's wall time; budget accordingly.

NOT cached deliberately: no save_dir is passed, so this cannot write a
half-converged entry into example_sim_cache/ or collide with production
results. It is a diagnostic, not a solve whose S-matrix we intend to keep.

Run:
    MMI4X4_PML_DPML=1.5 python mmi_4x4_90hybrid_pml_trapping_test.py
"""
import os

import numpy as np

from beef import MMI, FDTDSolver, SOI_220

DPML = float(os.environ.get("MMI4X4_PML_DPML", 1.5))
RES = float(os.environ.get("MMI4X4_PML_RES", 40.0))

WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 6.0
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 0.8
TAPER_LENGTH = 1.0
L_PI = 108.9195
DELTA_L = -2.75
MM_LENGTH = 0.75 * L_PI + DELTA_L   # 78.9396

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.22, use_neighbor_aperture=True,
)

solver = FDTDSolver(
    res=RES, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=DPML,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    num_input_modes=1, num_output_modes=1,
    output_monitor_centers=list(mmi.output_ys),
    output_monitor_widths=[WAVEGUIDE_OFFSET] * 4,
    material_stack=SOI_220, decay_tol=1e-11,
    decay_method="energy", energy_decay_dt=20.0,
)

print(f"PMLTEST: dpml={DPML} res={RES} (production baseline is dpml=0.5)",
      flush=True)
print(f"PMLTEST: cell Sy={solver.sim_width + 2*DPML:.2f} "
      f"Sz={solver.sim_height + 2*DPML:.2f} um "
      f"(baseline 10.00 x 3.00)", flush=True)
print("PMLTEST: baseline dpml=0.5 for comparison -- "
      "t=600 ratio 5.66e-02 | t=861 ratio 1.38e-02 | t=1201 ratio 6.09e-03",
      flush=True)

wavelengths = np.linspace(1.26, 1.36, 3)
# No save_dir -- diagnostic only, must not write a partial entry to the cache.
smat = solver.solve(mmi, wavelengths, port="bot", input_mode=1)
print("PMLTEST: run completed (reached decay_tol)", flush=True)
print({f"out_{j}": smat.transmission(port=f"out_{j}").tolist() for j in range(4)},
      flush=True)
