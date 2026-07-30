"""Real-meep verification for beef.solvers.fdtd.PortSpec's y-normal port
support (validation_markdown/ring_resonator_validation.md prerequisite (a)).

Marked @pytest.mark.slow, matching this repo's existing convention for
real-meep checks (tests/test_engine_agreement.py, test_eme_convergence.py,
test_golden_examples.py, test_mode_identity.py) — run manually via
`pytest tests/test_ring_orientation.py --runslow`, never part of the fast
suite.

Scope: the smallest real-physics check that exercises a y-normal MEASURE
port, the forced-reverse=False fix for solve_section()'s n_in==n_out
tie-break, and the new y0/y1 domain sizing — a standalone 90-degree bend
(RingResonatorBendSection), x-normal excite / y-normal measure, at reduced
resolution. Compared only for physical sanity (order of magnitude, no
divergence) against examples/ring_resonator/python/ring_resonator_hybrid.py's
own standalone-bend result, NOT a tight numeric match — that script runs at
production resolution/3-D, this at reduced resolution/2-D by design, per
this repo's existing RING_TEST_MODE convention.

The excite-side / flux-calibration path (a y-normal EXCITE port, needed for
RingResonatorCouplerSection's ring input) is a natural follow-on, not
covered here — it additionally needs RingResonatorCouplerSection.
section_prisms() validated first (see validation_markdown/
ring_resonator_validation.md).
"""
from __future__ import annotations

import numpy as np
import pytest

from beef.geometry.ring_bend import RingResonatorBendSection
from beef.materials import SOI_220
from beef.solvers.fdtd import FDTDSolver, PortSpec


@pytest.mark.slow
def test_bend_y_normal_measure_port_gives_sane_transmission():
    """x-normal excite -> y-normal measure through a 90-degree bend.

    Exercises: PortSpec(axis="y") on the measure side, y0/y1 domain sizing,
    and the forced reverse=False fix (n_modes_in == n_modes_out == 1 is
    exactly the tie-break-picks-the-wrong-side case verified against
    section_fdtd_layout()'s `reverse = n_out <= n_in` formula).
    """
    geometry = RingResonatorBendSection(waveguide_width=0.45, bend_radius=5.0, output_margin=0.5)
    port = geometry.out_port()

    solver = FDTDSolver(
        res=15.0,           # reduced — production is 40 (see RING_TEST_MODE convention)
        sim_width=3.0,
        sim_height=1.5,
        dpml=0.5,
        source_margin=0.0,
        monitor_margin=0.25,
        three_d=False,       # 2-D planar — production is 3-D
        material_stack=SOI_220,
        decay_tol=1e-9,       # loosened — production is 1e-11
    )

    wavelengths = np.array([1.53, 1.55, 1.57])
    aperture = 3.0
    y0 = -1.0
    y1 = port.y_boundary  # 5.5

    T = solver.solve_section(
        geometry, wavelengths, 0.0, geometry.bend_radius,
        n_modes_in=1, n_modes_out=1,
        stub_left=1.0,
        y0=y0, y1=y1,
        excite_port=PortSpec(axis="x", sign=1, center=0.0, width=aperture),
        measure_port=PortSpec(axis="y", sign=1, center=geometry.bend_radius, width=aperture),
        y_pml=solver.dpml,
    )

    transmission = np.abs(T[:, 0, 0]) ** 2
    assert np.all(np.isfinite(transmission)), f"non-finite transmission: {transmission}"
    # A 5 um radius 90-degree bend in 450x220 SOI is low-loss; at reduced
    # (2-D, res 15) resolution this is a sanity range, not a tight physical
    # prediction — see module docstring for why a tighter match isn't
    # attempted here.
    assert np.all(transmission > 0.3), (
        f"bend transmission implausibly low: {transmission} — likely a sign/"
        "orientation bug in the y-normal measure port, not just discretization "
        "error (a genuine 5 um bend does not lose this much power)"
    )
    assert np.all(transmission < 1.05), f"bend transmission > 1 (non-physical): {transmission}"
