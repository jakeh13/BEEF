"""Calibrating a HybridSolver body on the FDTD's own β (dc_validation.md item 6).

WHAT THIS IS FOR
----------------
The hybrid chain's answer depended on where you cut it: item 7 re-ran the whole
stub sweep and L₅₀ still drifted 0.1085 µm per µm of stub, dead straight, no
plateau. That drift is exactly ``2s·(1 - Δβ_FDTD/Δβ_EME)``, so it is a direct
readout of the engines disagreeing about Δβ.

The disagreement is NOT in the mode solve. `meow.compute_modes_meep()` and
`FDTDSolver.section_mode_neff()` both call `sim.get_eigenmode()` — the same MPB
— and their Δβ agree to 0.56%. Refitting item 6's stored transfer matrices per
mode locates it instead between MPB's mode index and the phase meep's
time-stepper actually accumulates, at res 40 on the coupler:

    mode      beta_MPB    beta_FDTD    excess
    even       9.42255      9.63292    +2.233%
    odd        9.32743      9.54767    +2.361%
    dbeta      0.09512      0.08525   -10.381%

Both biased the same way; the differential is 0.129 percentage points, and
``β/Δβ ≈ 99`` turns that into the whole 10.4%. Reaching it by mesh refinement
would need res ≈ 143, ~164× the 3-D cost.

So the body is propagated on the FDTD's β instead. That makes the cascade
self-consistent, which is the actual design goal — a hybrid result is only
meaningful if it does not depend on the cut — and it is a correctness fix
rather than a fudge, because the thing the hybrid exists to reproduce is a
full-device FDTD run in meep's own discretization.

Everything here is fast and simulation-free: the physics evidence is the stub
sweep re-run, not a unit test. What these guard is the extraction, which is
where the subtle failure lives.
"""
from __future__ import annotations

import inspect

import numpy as np
import pytest

from beef.solvers.fdtd import (
    _BETA_BRANCH_MARGIN_MAX, _BETA_FIT_RESIDUAL_MAX, _BETA_MIN_BRANCH_FRACTION,
    beta_from_phases, dbeta_from_phases,
)

# res 40, coupler at 1.55 µm, from example_output/verification/.
BETA_MPB = np.array([9.42255, 9.32743])       # section_mode_neff, via MPB
BETA_FDTD = np.array([9.63292, 9.54767])      # fitted from the FDTD's own phase
LENGTHS = np.array([2.0, 4.0, 6.0, 8.0])
WL = 1.55


def _wrapped(beta, lengths, offset=0.0):
    """arg(T_ii) as a monitor would report it: wrapped into (−π, π]."""
    return np.angle(np.exp(1j * (beta * np.asarray(lengths, float) + offset)))


class TestBetaFromPhases:
    """The per-mode fit. Ill-conditioned in a way that fails silently, which is
    why it gets more tests than the well-conditioned difference fit."""

    def test_recovers_a_known_rate_through_many_wraps(self):
        """β·ΔL is ~3 full turns at this spacing, so this cannot work by
        unwrapping — only by choosing the branch."""
        for beta in BETA_FDTD:
            got, resid, margin = beta_from_phases(
                LENGTHS, _wrapped(beta, LENGTHS), seed_beta=beta)
            assert got == pytest.approx(beta, rel=1e-9)
            assert resid < 1e-9
            assert margin < 1e-6

    def test_constant_offset_does_not_bias_the_rate(self):
        """Mode-matching and monitor-plane phase is a constant the slope must
        be immune to — the same reason the difference fit takes several
        lengths."""
        beta = BETA_FDTD[0]
        for offset in (-2.0, -0.3, 0.7, 3.1):
            got, _, _ = beta_from_phases(
                LENGTHS, _wrapped(beta, LENGTHS, offset), seed_beta=beta)
            assert got == pytest.approx(beta, rel=1e-9), offset

    def test_seed_picks_the_branch_but_does_not_set_the_answer(self):
        """The point on which the whole method stands or falls. Seeded with
        MPB's β — which is 2.2% away from the truth — it must still return the
        FDTD's β, not the seed."""
        beta_true, seed = BETA_FDTD[0], BETA_MPB[0]
        assert abs(beta_true - seed) / seed > 0.02, "premise: seed is genuinely off"
        got, _, margin = beta_from_phases(
            LENGTHS, _wrapped(beta_true, LENGTHS), seed_beta=seed)
        assert got == pytest.approx(beta_true, rel=1e-9)
        assert got != pytest.approx(seed, rel=1e-3)
        assert margin < _BETA_BRANCH_MARGIN_MAX

    def test_branches_are_a_third_of_beta_apart_so_the_seed_cannot_be_wrong(self):
        """Why seeding is safe rather than circular: neighbouring candidate
        slopes differ by 2π/ΔL = 3.14 rad/µm here, a 33% step, while the
        disagreement being measured is ~2%."""
        d_l = float(np.min(np.diff(LENGTHS)))
        branch_spacing = 2 * np.pi / d_l
        assert branch_spacing / BETA_MPB[0] > 0.3
        measured = abs(BETA_FDTD[0] - BETA_MPB[0])
        assert measured < 0.1 * branch_spacing

    def test_a_boundary_seed_shows_up_as_a_broken_fit(self):
        """FAILURE MODE 1. A seed half a branch away makes the per-length wrap
        counts disagree, so the points stop lying on a line. The residual is
        what catches this — the margin does not (it reads only 0.4)."""
        beta = BETA_FDTD[0]
        bad_seed = beta + np.pi / float(np.min(np.diff(LENGTHS)))   # half a branch
        _, resid, margin = beta_from_phases(
            LENGTHS, _wrapped(beta, LENGTHS), seed_beta=bad_seed)
        assert resid > _BETA_FIT_RESIDUAL_MAX, "residual must catch it"
        assert margin < _BETA_BRANCH_MARGIN_MAX, (
            "documents that the margin alone would NOT have caught this"
        )

    def test_a_whole_branch_error_is_invisible_in_the_fit(self):
        """FAILURE MODE 2, and the reason section_mode_beta() constrains the
        LENGTHS rather than trusting any after-the-fact check.

        Seeded a full branch off, every length agrees on that wrong branch: the
        fit is clean to 1e-14, the margin reads 0.000, and β is 33% wrong.
        Nothing in the returned numbers can distinguish it from a good fit —
        which is why the guard has to be a precondition on the spacing."""
        beta = BETA_FDTD[0]
        one_branch = 2 * np.pi / float(np.min(np.diff(LENGTHS)))
        got, resid, margin = beta_from_phases(
            LENGTHS, _wrapped(beta, LENGTHS), seed_beta=beta + one_branch)
        assert got == pytest.approx(beta + one_branch, rel=1e-9), "lands a branch off"
        assert resid < 1e-9, "and the fit looks perfect"
        assert margin < 1e-6, "and the margin sees nothing"
        assert abs(got / beta - 1) > 0.3, "while being 33% wrong"

    def test_needs_at_least_two_lengths(self):
        beta = BETA_FDTD[0]
        got, _, _ = beta_from_phases([2.0, 4.0], _wrapped(beta, [2.0, 4.0]), beta)
        assert got == pytest.approx(beta, rel=1e-9)


class TestTheSpacingPrecondition:
    """`section_mode_beta()` refuses widely-spaced lengths, because a
    whole-branch seed error cannot be detected after the fact (above). These
    check the arithmetic without needing meep."""

    def _max_safe_spacing(self, beta):
        return 2 * np.pi / (_BETA_MIN_BRANCH_FRACTION * beta)

    def test_item_6s_own_spacing_passes(self):
        """The configuration that produced the published numbers must remain
        legal, or the fix would disagree with the measurement it rests on."""
        assert float(np.max(np.diff(LENGTHS))) <= self._max_safe_spacing(BETA_FDTD[0])

    def test_the_naive_default_would_have_been_refused(self):
        """The default nearly shipped as [span/2, span]. On a 30 µm coupler
        body that is 15 µm apart — branches 0.42 rad/µm apart on a β of 9.6,
        i.e. far finer than the ~2% disagreement being measured, so the seed
        could not have identified the branch at all."""
        span = 30.0
        assert span / 2.0 > self._max_safe_spacing(BETA_FDTD[0])

    def test_the_actual_default_is_span_independent(self):
        """Spacing is capped by β, not by how big the body happens to be."""
        cap = self._max_safe_spacing(BETA_FDTD[0]) / 2.0
        for span in (4.0, 30.0, 300.0):
            spacing = min(span / 2.0, cap)
            assert spacing <= self._max_safe_spacing(BETA_FDTD[0])
            assert 2.0 * spacing <= span, "both lengths must fit inside the body"

    def test_a_safe_spacing_leaves_an_order_of_magnitude(self):
        """25% of β against a ~2% measured disagreement."""
        branch = 2 * np.pi / self._max_safe_spacing(BETA_FDTD[0])
        measured = abs(BETA_FDTD[0] - BETA_MPB[0])
        assert branch / measured > 10


class TestAgreesWithTheDifferenceFit:
    """The two routes to Δβ must agree — that cross-check is what caught the
    sign-search bug in the original analysis, and it is cheap to keep."""

    def test_per_mode_difference_matches_dbeta_from_phases(self):
        phases = [_wrapped(b, LENGTHS, off)
                  for b, off in zip(BETA_FDTD, (0.3, -1.1))]
        per_mode = [beta_from_phases(LENGTHS, p, s)[0]
                    for p, s in zip(phases, BETA_MPB)]
        direct, _, resid = dbeta_from_phases(LENGTHS, phases[0] - phases[1])
        assert per_mode[0] - per_mode[1] == pytest.approx(direct, abs=1e-9)
        assert resid < 1e-9

    def test_the_difference_fit_needs_no_seed_because_dbeta_is_small(self):
        """Why there are two functions. Δβ·ΔL is 0.17 rad at the coupler's
        spacing — nowhere near a wrap — so np.unwrap is sufficient there and
        useless for the absolute rate."""
        dbeta = BETA_FDTD[0] - BETA_FDTD[1]
        assert dbeta * float(np.min(np.diff(LENGTHS))) < np.pi
        assert BETA_FDTD[0] * float(np.min(np.diff(LENGTHS))) > 6 * np.pi


class TestTheAmplificationItself:
    """Pins the mechanism, so a future reader does not have to re-derive why a
    2% per-mode bias matters."""

    def test_common_mode_bias_is_small_but_dbeta_error_is_not(self):
        excess = BETA_FDTD / BETA_MPB - 1.0
        assert np.all(excess > 0.02) and np.all(excess < 0.025), excess
        assert abs(excess[0] - excess[1]) < 0.002, "differential is tiny"

        d_mpb = BETA_MPB[0] - BETA_MPB[1]
        d_fdtd = BETA_FDTD[0] - BETA_FDTD[1]
        assert BETA_MPB[0] / d_mpb == pytest.approx(99, abs=1)      # the lever
        assert abs(d_fdtd / d_mpb - 1) > 0.10                        # the damage

    def test_dbeta_error_is_reconstructed_by_the_lever(self):
        """Δβ_FDTD = Δβ_MPB + β₀e₀ − β₁e₁ exactly — the error is entirely the
        differential bias, not some separate effect."""
        e = BETA_FDTD / BETA_MPB - 1.0
        predicted = (BETA_MPB[0] - BETA_MPB[1]) + BETA_MPB[0] * e[0] - BETA_MPB[1] * e[1]
        assert predicted == pytest.approx(BETA_FDTD[0] - BETA_FDTD[1], abs=1e-9)

    def test_mesh_refinement_is_not_a_viable_route(self):
        """At dx^2 scaling the differential needs ~3.6x more resolution for 1%
        on dbeta, i.e. res ~143 and ~164x the 3-D cost. Recorded as a test so
        'just raise the resolution' is answered with a number."""
        differential = abs((BETA_FDTD / BETA_MPB - 1.0)[0] - (BETA_FDTD / BETA_MPB - 1.0)[1])
        lever = BETA_MPB[0] / (BETA_MPB[0] - BETA_MPB[1])
        needed = differential / (0.01 / lever)          # factor of error reduction
        res_needed = 40.0 * np.sqrt(needed)
        assert res_needed > 120
        assert (res_needed / 40.0) ** 4 > 100           # 3 space + 1 time


class TestWiring:
    """The plumbing, without running any physics."""

    def test_eme_accepts_neff_override_on_both_section_paths(self):
        from beef import EMESolver
        for name in ("solve_section", "solve_section_length_sweep"):
            params = inspect.signature(getattr(EMESolver, name)).parameters
            assert "neff_override" in params, name
            assert params["neff_override"].default is None, name

    def test_hybrid_defaults_to_off(self):
        """Every recorded golden predates calibration, so the flag must not
        change anything unless asked for."""
        from beef import HybridSolver
        params = inspect.signature(HybridSolver.__init__).parameters
        assert params["calibrate_body_neff"].default is False
        assert params["calibrate_lengths"].default is None

    def test_override_replaces_neff_and_leaves_fields_alone(self):
        """`Mode` is frozen, so this has to copy. Field profiles must survive
        untouched — only propagation is being taken from the FDTD, not the
        overlap integrals."""
        from beef import EMESolver

        class FakeMode:
            def __init__(self, neff, tag):
                self.neff, self.tag = neff, tag

            def model_copy(self, update):
                m = FakeMode(self.neff, self.tag)
                m.__dict__.update(update)
                return m

        modes = [FakeMode(2.325, "even"), FakeMode(2.301, "odd")]
        out = EMESolver._override_neff(modes, [2.3766, 2.3556])
        assert [m.neff for m in out] == [2.3766 + 0j, 2.3556 + 0j]
        assert [m.tag for m in out] == ["even", "odd"], "fields must be preserved"
        assert [m.neff for m in modes] == [2.325, 2.301], "input must not be mutated"

    def test_override_length_mismatch_is_refused(self):
        from beef import EMESolver
        with pytest.raises(ValueError, match="neff_override has"):
            EMESolver._override_neff([object(), object()], [1.0])

    def test_only_bodies_get_calibrated_and_from_a_neighbour(self):
        """A section that can measure its own propagation is left alone; a body
        that cannot inherits β from an adjacent FDTD section, so it ends up on
        the discretization of the sections it is spliced to. That is the whole
        point — the chain has to agree with itself."""
        from beef import HybridSolver

        class FakeFDTD:
            def section_mode_beta(self, geom, wl, x0, x1, n, lengths=None, **kw):
                return np.full((len(wl), n), BETA_FDTD[0])

            def solve_section(self, *a, **k):
                pass

        class FakeEME:
            def solve_section(self, *a, neff_override=None, **k):
                pass

        h = HybridSolver([FakeFDTD(), FakeEME(), FakeFDTD()], [0, 1, 2, 3], [2, 2],
                         calibrate_body_neff=True)
        assert h._calibration_source(0) is None, "FDTD measures its own"
        assert h._calibration_source(2) is None
        assert isinstance(h._calibration_source(1), FakeFDTD), "body borrows one"

        off = HybridSolver([FakeFDTD(), FakeEME(), FakeFDTD()], [0, 1, 2, 3], [2, 2])
        assert off._calibration_source(1) is None, "and nothing happens when off"

    def test_beta_is_converted_to_neff_the_way_meow_propagates(self):
        """meow applies exp(2j·π·neff/wl·L), so neff = β·λ/2π is what makes the
        body accumulate exactly the FDTD's phase per µm."""
        from beef import HybridSolver

        class FakeFDTD:
            def section_mode_beta(self, geom, wl, x0, x1, n, lengths=None, **kw):
                return np.full((len(wl), n), BETA_FDTD[0])

            def solve_section(self, *a, **k):
                pass

        class FakeEME:
            def solve_section(self, *a, neff_override=None, **k):
                pass

        h = HybridSolver([FakeFDTD(), FakeEME(), FakeFDTD()], [0, 1, 2, 3], [2, 2],
                         calibrate_body_neff=True)
        neff = h._calibrated_neff(None, np.array([WL]), 1, 2, {}, None, False)
        assert neff.shape == (1, 2)
        assert neff[0, 0] == pytest.approx(BETA_FDTD[0] * WL / (2 * np.pi))
        # and it is above MPB's, by the same 2.2% the phase was
        assert neff[0, 0] / (BETA_MPB[0] * WL / (2 * np.pi)) == pytest.approx(
            BETA_FDTD[0] / BETA_MPB[0])

    def test_neff_fingerprint_separates_calibrated_from_uncalibrated(self):
        """Without this in the cache config a calibrated run and an
        uncalibrated one would share an entry and serve each other."""
        from beef.solvers.eme import _neff_fingerprint
        assert _neff_fingerprint(None) is None
        a = _neff_fingerprint(np.array([[2.325, 2.301]]))
        b = _neff_fingerprint(np.array([[2.3766, 2.3556]]))
        assert a != b
        assert a == _neff_fingerprint(np.array([[2.325, 2.301]])), "must be stable"

    def test_fingerprint_is_json_safe(self):
        """_hash_config runs json.dumps on the whole config, and numpy arrays
        are not serializable — this is the guard for that."""
        import json
        from beef.solvers.eme import _neff_fingerprint
        json.dumps(_neff_fingerprint(np.array([[2.325 + 0j, 2.301 + 0j]])))
