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
RingResonatorCouplerSection's ring input) is covered by the second test
below, added 2026-07-31 — by RECIPROCITY on the same bend, which gives a
quantitative expectation rather than another sanity range. It does not need
RingResonatorCouplerSection.section_prisms() validated first after all (an
earlier revision of this docstring said it did): the bend alone has both an
x-normal and a y-normal port, so swapping which one excites exercises the
y-normal excite + flux-calibration path on already-validated geometry.
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


def _bend_setup():
    """The same bend + solver both tests below use, so a difference between
    them can only come from which port excites."""
    geometry = RingResonatorBendSection(
        waveguide_width=0.45, bend_radius=5.0, output_margin=0.5)
    solver = FDTDSolver(
        res=15.0, sim_width=3.0, sim_height=1.5, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=False,
        material_stack=SOI_220, decay_tol=1e-9,
    )
    return geometry, solver


@pytest.mark.slow
@pytest.mark.xfail(
    strict=True,
    reason="KNOWN BUG, found by this test 2026-07-31: the y-normal excite "
           "path's flux-calibrated normalization violates reciprocity by a "
           "geometry-invariant factor ~0.73. Localized to that path (not the "
           "sign=-1 conventions) by two controls -- see "
           "validation_markdown/ring_resonator_validation.md's boxed blocker "
           "under prerequisite (a). strict=True so this starts FAILING the "
           "moment the normalization is fixed, which is the signal to drop "
           "this marker and let the assertions stand.",
)
def test_bend_y_normal_excite_matches_reciprocity():
    """y-normal EXCITE (+ flux-calibrated normalization) against the x-normal
    excite direction on the SAME bend, by reciprocity.

    This is the path `PortSpec`'s own docstring flags as unverified by any
    real meep run: a y-normal `excite_port` additionally triggers
    flux-calibrated normalization (`src_alpha = sqrt(|flux|) * exp(i*arg
    (norm_coeff))` rather than the coefficient alone), which
    `examples/ring_resonator/python/ring_resonator_hybrid.py` found
    empirically necessary for its backward, off-axis ring launch but which no
    test exercised. `RingResonatorCouplerSection`'s ring input needs exactly
    this path, so the ring cannot go on `HybridSolver` until it is trusted.

    Why reciprocity rather than another sanity range: a passive, reciprocal,
    linear structure has S_21 == S_12, so exciting the bend's x port and
    measuring its y port must give the same |T| as exciting y and measuring
    x — with the SAME geometry, mesh, and mode normalization on both runs.
    That is a real quantitative prediction, and it is specifically sensitive
    to the thing being tested: if the flux calibration mis-normalizes the
    y-launched source, the two directions' |T| differ by that factor while
    each still looks individually "plausible" against an order-of-magnitude
    check.

    `sign=-1` on the y excite port is deliberate — it puts the source at the
    HIGH y edge launching -y (`PortSpec`'s excite-side edge convention),
    matching `ring_resonator_hybrid.py`'s own ring launch
    (`eig_kpoint=mp.Vector3(0, -1, 0)`, source at `Y_TOP - SOURCE_MARGIN`),
    which is the only sign that script ever exercised.
    """
    geometry, solver = _bend_setup()
    port = geometry.out_port()
    wavelengths = np.array([1.53, 1.55, 1.57])
    aperture = 3.0
    y0, y1 = -1.0, port.y_boundary
    common = dict(
        n_modes_in=1, n_modes_out=1, stub_left=1.0,
        y0=y0, y1=y1, y_pml=solver.dpml,
    )

    # x in -> y out (the direction test 1 above already validated).
    T_xy = solver.solve_section(
        geometry, wavelengths, 0.0, geometry.bend_radius,
        excite_port=PortSpec(axis="x", sign=1, center=0.0, width=aperture),
        measure_port=PortSpec(axis="y", sign=1, center=geometry.bend_radius, width=aperture),
        **common,
    )

    # y in -> x out: the reverse. excite sign=-1 => HIGH y edge, launching
    # -y; measure sign=-1 => LOW x edge, reading -x-forward light.
    T_yx = solver.solve_section(
        geometry, wavelengths, 0.0, geometry.bend_radius,
        excite_port=PortSpec(axis="y", sign=-1, center=geometry.bend_radius, width=aperture),
        measure_port=PortSpec(axis="x", sign=-1, center=0.0, width=aperture),
        **common,
    )

    t_xy = np.abs(T_xy[:, 0, 0])
    t_yx = np.abs(T_yx[:, 0, 0])
    assert np.all(np.isfinite(t_yx)), f"non-finite y-excite amplitude: {t_yx}"
    # Guard the trivial-pass case first: an all-zero t_yx would satisfy a
    # ratio test vacuously if the ratio were computed the wrong way round.
    assert np.all(t_yx > 0.1), (
        f"y-normal excite produced almost no output: {t_yx} — the source "
        "likely failed to launch (check the flux-calibration path and the "
        "excite-side sign/edge convention), rather than this being loss"
    )

    # Reciprocity: same |S| both directions. The tolerance is loose because
    # both runs are 2-D at res 15 and the two directions' monitors sit on
    # different-orientation planes (so their discretization error is not
    # identical) -- but a mis-scaled flux calibration would show up as a
    # systematic ratio well outside this, not a few percent of mesh noise.
    ratio = t_yx / t_xy
    assert np.allclose(ratio, 1.0, rtol=0.15), (
        f"reciprocity violated: |T| x->y = {t_xy}, y->x = {t_yx}, "
        f"ratio = {ratio}. A consistent ratio far from 1 points at the "
        "y-normal excite path's flux-calibrated normalization being "
        "mis-scaled (PortSpec's own docstring flags sign=+1 there as applied "
        "'conservatively' and never tested); a ratio near a clean factor "
        "like 2, 0.5 or sqrt(2) is especially suggestive of a normalization "
        "convention rather than discretization error."
    )
