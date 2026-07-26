"""FDTD/EME engine agreement on a uniform coupled pair (notes.md Solver #9).

THE KNOWN FAILURE
-----------------
Over the *same* straight coupled-pair cross-section, the two engines disagree
about how fast even/odd phase accumulates. At the production meshing (both at
res 40, SOI_220, gap 0.25):

    dbeta_FDTD / dbeta_EME = 0.891 at 1.55 um

measured directly by examples/verification/python/verify_item_6_engine_agreement.py
(dc_validation.md item 6, job 11489846). That independently reproduced the 0.891
the stub sweep had only been able to *infer* from a linear fit to L_50 through a
whole hybrid chain, so the effect is real and lives in the uniform region -- not
in the bends, the seam mode-matching, or L_50 extraction.

Most of it is not engine disagreement at all. dbeta_EME is not converged at
res 40:

    EME res     20       40       60       80      100
    dbeta    0.1125   0.0957   0.0926   0.0894   0.0895

so res 40 overestimates dbeta_EME by 7%, and roughly 6 of the 11 points are EME
discretization error. With EME converged the ratio is 0.953 -- a residual ~4.7%
that is still larger than FDTD's own 2.7% resolution scatter and still
unexplained. That residual is what the slow test below pins.

It matters because HybridSolver cascades FDTD sections onto EME bodies. A
phase-rate mismatch at the join means the chain's answer depends on where the
cut is made, which is exactly the non-convergence dc_validation.md was opened to
explain -- and the 7% EME half of it is present in every EME result in the repo,
not only this one.

test_engines_agree_on_uniform_coupled_pair is marked xfail(strict=True): it
encodes the requirement, not today's behaviour, and will fail loudly *as a pass*
the moment someone closes the gap -- which is the notification wanted.

The extraction tests below need no simulation and run in the fast suite. They
guard the measurement method itself: a bug in the slope fit would make the
engine comparison meaningless in either direction.
"""
from __future__ import annotations

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# The measurement method, kept here rather than imported from
# examples/verification/ -- that folder is temporary and gets deleted when the
# directional coupler work is signed off, but this check is meant to outlive it.
# ---------------------------------------------------------------------------
def dbeta_from_phases(lengths, dphi):
    """Even/odd phase rate (rad/um) from dphi measured at several lengths.

    dphi(L) = arg(T_00) - arg(T_11) = dbeta*L + const, where the constant is
    mode-matching and monitor-plane phase, identical at every length. The SLOPE
    is the physics; the intercept is discarded. One length cannot do this --
    dphi(L)/L would carry the offset as a 1/L error.

    Returns (dbeta, intercept, rms_residual). The residual is the evidence that
    dphi really is linear in L, i.e. that this is one clean supermode pair.
    """
    lengths = np.asarray(lengths, dtype=float)
    y = np.unwrap(np.asarray(dphi, dtype=float))
    A = np.vstack([lengths, np.ones_like(lengths)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = float(np.sqrt(np.mean((y - (slope * lengths + intercept)) ** 2)))
    return float(slope), float(intercept), resid


class TestPhaseSlopeExtraction:
    """The fit that turns four transfer matrices into one phase rate."""

    LENGTHS = np.array([2.0, 4.0, 6.0, 8.0])

    def test_recovers_a_known_rate(self):
        dbeta = 0.0956
        dphi = dbeta * self.LENGTHS
        got, _, resid = dbeta_from_phases(self.LENGTHS, dphi)
        assert got == pytest.approx(dbeta, rel=1e-12)
        assert resid == pytest.approx(0.0, abs=1e-12)

    def test_constant_offset_does_not_bias_the_rate(self):
        """The whole reason for using several lengths: mode-matching phase at
        the monitor planes is a constant the slope must be immune to."""
        dbeta = 0.0956
        for offset in (-2.0, -0.3, 0.0, 0.7, 3.1):
            dphi = dbeta * self.LENGTHS + offset
            got, intercept, _ = dbeta_from_phases(self.LENGTHS, dphi)
            assert got == pytest.approx(dbeta, rel=1e-12), offset
            assert intercept == pytest.approx(offset, abs=1e-9)

    def test_single_length_would_have_been_wrong(self):
        """Documents why this API takes an array. With a 0.5 rad offset, the
        naive dphi/L is off by 160% at L=2 and still 33% at L=8 -- the error
        shrinks with L but never vanishes, so no single length is safe."""
        dbeta, offset = 0.0956, 0.5
        naive = (dbeta * self.LENGTHS + offset) / self.LENGTHS
        assert abs(naive[0] / dbeta - 1) > 1.5
        assert abs(naive[-1] / dbeta - 1) > 0.3
        fitted, _, _ = dbeta_from_phases(self.LENGTHS, dbeta * self.LENGTHS + offset)
        assert fitted == pytest.approx(dbeta, rel=1e-12)

    def test_unwraps_a_phase_that_crosses_pi(self):
        """arg() is confined to (-pi, pi]. A genuine dphi past pi comes back
        wrapped and would fit a negative slope if taken at face value."""
        dbeta = 0.5                      # rad/um -> 4.0 rad at L=8, past pi
        true = dbeta * self.LENGTHS
        wrapped = np.angle(np.exp(1j * true))
        assert not np.allclose(wrapped, true), "test needs an actually-wrapped case"
        got, _, resid = dbeta_from_phases(self.LENGTHS, wrapped)
        assert got == pytest.approx(dbeta, rel=1e-9)
        assert resid == pytest.approx(0.0, abs=1e-9)

    def test_residual_reports_nonlinearity(self):
        """A curved dphi(L) means the two-mode picture does not hold -- more
        than two modes carrying power, or a leaky one. The residual must show
        it rather than quietly returning a meaningless slope."""
        dphi = 0.0956 * self.LENGTHS + 0.01 * self.LENGTHS ** 2
        _, _, resid = dbeta_from_phases(self.LENGTHS, dphi)
        assert resid > 1e-2

    def test_matches_the_measured_fdtd_run(self):
        """Regression against the real numbers from job 11489846 at 1.55 um.
        These four phases are what the FDTD sections actually returned; the
        recorded slope is what dc_validation.md item 6 quotes."""
        dphi = np.array([0.1494, 0.3199, 0.4904, 0.6609])
        dbeta, _, resid = dbeta_from_phases(self.LENGTHS, dphi)
        assert dbeta == pytest.approx(0.08525, abs=5e-5)
        assert resid < 1e-3, "the measured phases were linear to <1e-4 rad"


# ---------------------------------------------------------------------------
# The real thing. Needs meep + meow, so --runslow only.
# ---------------------------------------------------------------------------
WL = 1.55
SECTION_LENGTHS = (2.0, 6.0)   # two lengths is the minimum that kills the offset

# The two engines are meshed at DIFFERENT resolutions on purpose, and dropping
# either to save test time would measure the wrong thing.
#
# EME at the production res=40 is not converged: dbeta_EME runs 0.1125 (res 20)
# -> 0.0957 (40) -> 0.0926 (60) -> 0.0894 (80) -> 0.0895 (100), so res 40
# overestimates it by 7%. Meshing both engines at 40 gives ratio 0.891, but ~6
# of those 11 points are EME discretization error rather than engine
# disagreement -- and at res 20 the "gap" inflates to ~25%, which would be
# almost entirely convergence error. res 80 is the first converged point.
#
# FDTD can be dropped to 20 without weakening the check, because it converges
# far earlier than EME on the same nominal grid: the archived resolution sweep
# moves its implied dbeta only 2.7% across res 30-90 (non-monotonically -- noise,
# not a trend), and res 20 sits within 0.5% of res 40. meep's subpixel smoothing
# is the likely reason. That buys ~16x in a 3-D run, which is what makes this
# affordable serially at all -- the res-40 measurement needed 48 MPI ranks.
RES_FDTD = 20.0
RES_EME = 80.0

DEVICE = dict(
    waveguide_width=0.45, gap=0.25, coupling_length=30.0,
    bend_length=9.0, waveguide_offset=3.0, bend_type="euler",
    waveguide_height=0.22,
)


def _dbeta_fdtd():
    from beef import DirectionalCoupler, FDTDSolver, SOI_220

    dc = DirectionalCoupler(**DEVICE)
    solver = FDTDSolver(
        res=RES_FDTD, sim_width=6.0, sim_height=1.5, dpml=1.0,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        half_domain_aperture=True, material_stack=SOI_220, decay_tol=1e-9,
    )
    wl = np.array([WL])
    phases = []
    for L in SECTION_LENGTHS:
        # Centred at x=0, so every cut sits 15 um clear of the nearest S-bend:
        # the section is pure uniform coupled pair, which is the only geometry
        # both engines can be asked about without ambiguity.
        T = solver.solve_section(
            dc, wl, -L / 2.0, +L / 2.0, n_modes_in=2, n_modes_out=2,
            # Full non-PML domain on the pair's shared y=0 axis. Without this
            # the monitor inherits one arm's centreline and the "even/odd"
            # decomposition is not even/odd at all.
            monitor_center=0.0, monitor_size=solver.sim_width,
        )
        # A uniform region cannot convert even <-> odd. If it does, the basis
        # is wrong and every phase below is meaningless.
        off = max(abs(T[0, 0, 1]), abs(T[0, 1, 0]))
        assert off < 0.05, f"supermode conversion {off:.3f} in a uniform section"
        phases.append(np.angle(T[0, 0, 0]) - np.angle(T[0, 1, 1]))
    dbeta, _, resid = dbeta_from_phases(SECTION_LENGTHS, phases)
    assert resid < 1e-2, f"dphi(L) not linear (residual {resid:.2e})"
    return dbeta


def _dbeta_eme():
    import meow as mw
    from beef import DirectionalCoupler, EMESolver, SOI_220, extrude_polygons

    dc = DirectionalCoupler(**DEVICE)
    eme = EMESolver(
        num_modes=[2], num_cells=[1], lengths=[1.0],
        oxide_width=6.0, oxide_height=1.5, res=RES_EME,
        material_stack=SOI_220, mode_solver="meep",
    )
    L = SECTION_LENGTHS[0]
    gds = dc.section_gds(-L / 2.0, +L / 2.0, eme.oxide_width)
    cells = mw.create_cells(
        structures=extrude_polygons(gds, eme.material_stack.extrusion_rules),
        mesh=eme._build_mesh(), Ls=[L],
    )
    cs = mw.CrossSection.from_cell(
        cell=cells[0], env=mw.Environment(wl=WL, T=eme.temperature))
    modes = eme.compute_modes(cs, 2)
    neff = [float(np.real(m.neff)) for m in modes]
    # Descending n_eff, so [0] is the symmetric (even) supermode. Asserted
    # because a reorder would silently flip the sign of dbeta.
    assert neff[0] > neff[1], f"supermodes not in descending n_eff order: {neff}"
    return 2 * np.pi * (neff[0] - neff[1]) / WL


@pytest.mark.slow
@pytest.mark.xfail(
    strict=True,
    reason="dc_validation.md item 6: with EME converged (res 80), FDTD and EME "
           "still disagree by ~5% on dbeta for an identical uniform coupled "
           "pair -- 0.953 at the res-40 FDTD meshing this was measured on, and "
           "the same to within FDTD's own convergence at the res 20 used here. "
           "(Meshing BOTH engines at 40 reads 0.891, but ~6 of those 11 points "
           "are EME under-convergence rather than engine disagreement.) The "
           "residual is unexplained and larger than FDTD's own 2.7% resolution "
           "scatter. When this xfail starts passing, the gap has been closed: "
           "delete the marker and say what fixed it.",
)
def test_engines_agree_on_uniform_coupled_pair():
    """FDTD and EME must agree on dbeta for the same straight coupled pair.

    This is the precondition for cascading them in HybridSolver at all: if the
    two disagree about phase rate, the chain's answer depends on where the cut
    between them is made.
    """
    ratio = _dbeta_fdtd() / _dbeta_eme()
    assert ratio == pytest.approx(1.0, abs=0.02), (
        f"dbeta_FDTD/dbeta_EME = {ratio:.4f}; the engines disagree by "
        f"{abs(1 - ratio) * 100:.1f}% on a uniform coupled pair"
    )
