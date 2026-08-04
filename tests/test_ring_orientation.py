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


def test_port_aperture_crossing_pml_raises():
    """An explicit PortSpec aperture that runs past the PML inner edge must
    RAISE, not be silently clipped by meep.

    Regression for the defect found 2026-08-01
    (ring_resonator_validation.md's blocker-0 investigation): this exact
    geometry — the one the slow test below used to ship with — puts BOTH
    ports' apertures partly inside PML. meep clips them, the missing power
    reads as bend loss, and |T| came back 0.711 where the correct answer is
    0.978. A 27% error that looked physically plausible and passed a
    `> 0.3` assertion.

    Fast (no meep run): the guard fires while the domain is being set up,
    before any time stepping.
    """
    pytest.importorskip("meep")  # this test builds a real mp.Medium; fast CI has no meep
    geometry = RingResonatorBendSection(waveguide_width=0.45, bend_radius=5.0, output_margin=0.5)
    port = geometry.out_port()
    solver = FDTDSolver(
        res=15.0, sim_width=3.0, sim_height=1.5, dpml=0.5, source_margin=0.0,
        monitor_margin=0.25, three_d=False, material_stack=SOI_220, decay_tol=1e-9,
    )
    with pytest.raises(ValueError, match="crosses the PML"):
        solver.solve_section(
            geometry, np.array([1.55]), 0.0, geometry.bend_radius,
            n_modes_in=1, n_modes_out=1,
            stub_left=1.0,               # no stub_right -> y port sits AT the PML edge
            y0=-1.0, y1=port.y_boundary,  # too tight -> x port aperture also overruns
            excite_port=PortSpec(axis="x", sign=1, center=0.0, width=3.0),
            measure_port=PortSpec(axis="y", sign=1, center=geometry.bend_radius, width=3.0),
            y_pml=solver.dpml,
        )


@pytest.mark.slow
def test_bend_y_normal_measure_port_gives_sane_transmission():
    """x-normal excite -> y-normal measure through a 90-degree bend.

    Exercises: PortSpec(axis="y") on the measure side, y0/y1 domain sizing,
    and the forced reverse=False fix (n_modes_in == n_modes_out == 1 is
    exactly the tie-break-picks-the-wrong-side case verified against
    section_fdtd_layout()'s `reverse = n_out <= n_in` formula).

    RE-BASELINED 2026-08-01. The original version of this test ran with
    stub_right=0 and y0=-1.0, which put both port apertures inside PML (see
    test_port_aperture_crossing_pml_raises above) and measured 0.711 — 27%
    low — while asserting only `> 0.3`, a range far too loose to notice. The
    domain now gives both ports clearance and the assertion is tight enough
    to fail if that regresses.
    """
    geometry = RingResonatorBendSection(waveguide_width=0.45, bend_radius=5.0, output_margin=0.5)
    port = geometry.out_port()

    solver = FDTDSolver(
        res=15.0,           # reduced — production is 40 (see RING_TEST_MODE convention)
        sim_width=3.0,
        sim_height=1.5,
        dpml=0.5,
        source_margin=0.0,
        monitor_margin=0.75,  # clear of the source near field (blocker-0 finding 1)
        three_d=False,       # 2-D planar — production is 3-D
        material_stack=SOI_220,
        decay_tol=1e-9,       # loosened — production is 1e-11
    )

    wavelengths = np.array([1.53, 1.55, 1.57])
    aperture = 3.0
    # stub_right adds x headroom past the y-normal port. For this geometry it
    # draws NO extra polygons (_section_polygons emits nothing past x=R), so
    # it is pure domain padding. y0 is lowered so the x port's own aperture
    # also clears PML.
    y0 = -2.0
    y1 = port.y_boundary  # 5.5

    T = solver.solve_section(
        geometry, wavelengths, 0.0, geometry.bend_radius,
        n_modes_in=1, n_modes_out=1,
        stub_left=1.0, stub_right=2.0,
        y0=y0, y1=y1,
        excite_port=PortSpec(axis="x", sign=1, center=0.0, width=aperture),
        measure_port=PortSpec(axis="y", sign=1, center=geometry.bend_radius, width=aperture),
        y_pml=solver.dpml,
    )

    transmission = np.abs(T[:, 0, 0]) ** 2
    assert np.all(np.isfinite(transmission)), f"non-finite transmission: {transmission}"
    # A 5 um radius 90-degree bend in 450x220 SOI is genuinely near-lossless.
    # Measured 0.9773 in amplitude (0.955 in power) at this configuration, so
    # this bound is tight enough to catch a clipped aperture (which cost 27%)
    # while still allowing 2-D/res-15 discretization slack.
    assert np.all(transmission > 0.90), (
        f"bend transmission implausibly low: {transmission} — a 5 um bend is "
        "near-lossless, so this points at a clipped port aperture or a sign/"
        "orientation bug, not discretization error"
    )
    assert np.all(transmission < 1.05), f"bend transmission > 1 (non-physical): {transmission}"


def _bend_setup():
    """The same bend + solver both tests below use, so a difference between
    them can only come from which port excites."""
    geometry = RingResonatorBendSection(
        waveguide_width=0.45, bend_radius=5.0, output_margin=0.5)
    solver = FDTDSolver(
        res=15.0, sim_width=3.0, sim_height=1.5, dpml=0.5,
        # monitor_margin 0.75, not 0.25: at 0.25 the norm plane sits in the
        # y-source's near field, where |alpha| and sqrt(flux) disagree by
        # ~1.54x. Moving it out changes |T| by only 0.3%, so this is NOT the
        # reciprocity violation -- it just removes a confound from it.
        source_margin=0.0, monitor_margin=0.75, three_d=False,
        material_stack=SOI_220, decay_tol=1e-9,
    )
    return geometry, solver


@pytest.mark.slow
@pytest.mark.xfail(
    strict=True,
    reason="NOT A BUG -- this configuration is 2-D at res 15, where the "
           "reciprocity ratio is ~0.75 from DISCRETIZATION alone. Resolved "
           "2026-08-01: in 3-D on clean geometry the ratio converges as h^2 "
           "(res 15 -> 0.8553, res 25 -> 0.9476; fitted order p=1.988 against "
           "meep's formal 2nd order, agreeing to 0.6%), extrapolating to 0.979 "
           "at res 40. A normalization or orientation bug would be "
           "resolution-INDEPENDENT, so the y-normal excite path and its "
           "flux calibration are correct. The earlier claim in this marker -- "
           "that the flux-calibrated normalization was mis-scaled -- is "
           "WITHDRAWN. This case stays xfail because 2-D/res-15 genuinely "
           "cannot show reciprocity here, not because anything is broken; the "
           "real check is examples/verification/python/"
           "verify_ring_bend_reciprocity_3d.py. strict=True so that if this "
           "ever passes, the tolerance or the physics changed and someone "
           "should look. See validation_markdown/ring_resonator_validation.md.",
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
    # Clean apertures on BOTH ports (blocker-0 finding 3) — the original
    # y0=-1.0/stub_right=0 domain put both inside PML and would now raise.
    # The reciprocity violation this test asserts on SURVIVES this fix
    # (measured 0.751 with clean apertures), which is exactly why the
    # xfail below still stands and why the clipping was never the cause.
    y0, y1 = -2.0, port.y_boundary
    common = dict(
        n_modes_in=1, n_modes_out=1, stub_left=1.0, stub_right=2.0,
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
