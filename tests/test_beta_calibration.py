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
    _BETA_BRANCH_MARGIN_MAX, _BETA_BRANCH_MARGIN_RAISE, _BETA_FIT_RESIDUAL_MAX,
    _BETA_MIN_BRANCH_FRACTION, beta_fit_from_phases, beta_from_phases,
    dbeta_from_phases,
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

    def test_a_boundary_seed_is_undecidable_and_the_margin_is_what_says_so(self):
        """REWRITTEN 2026-08-01. This test used to be
        `test_a_boundary_seed_shows_up_as_a_broken_fit` and asserted
        `resid > _BETA_FIT_RESIDUAL_MAX`.

        That tell is gone, deliberately. It was never a property of the physics
        — it was this module rounding wrap counts independently per length, so a
        seed half a branch off made them mutually inconsistent and the points
        stopped lying on a line. The fit now chooses a self-consistent turn
        vector, so the residual is clean on whichever branch it lands.

        At EXACTLY half a branch there is no right branch: the seed carries
        literally no information, and `|beta - seed|` ties too, so no
        tie-breaking rule can help. What catches it is `branch_ambiguity` (and
        its scaled twin `margin`), which is the honest report that the one
        seed-dependent decision was a coin flip."""
        beta = BETA_FDTD[0]
        bad_seed = beta + np.pi / float(np.min(np.diff(LENGTHS)))   # half a branch
        f = beta_fit_from_phases(LENGTHS, _wrapped(beta, LENGTHS), seed_beta=bad_seed)
        assert f.resid < 1e-9, "the residual can no longer be the guard for this"
        assert f.branch_ambiguity == pytest.approx(0.5), "a genuine coin flip"
        assert f.margin == pytest.approx(1.0)
        assert f.margin > _BETA_BRANCH_MARGIN_RAISE, "and THIS is what raises now"
        assert f.turn_ambiguity < 1e-9, "the data-determined turns were never in doubt"
        assert f.turns_shifted, "the fix did move off the seed's own per-length rounding"

    def test_a_near_boundary_seed_is_now_recovered(self):
        """The point of the 2026-08-01 change. At 0.35 of a branch — 17x the
        ~2% real FDTD-vs-MPB gap, and far past where the old code broke (it
        tolerated 0.000 relative seed error at these lengths for an adverse
        intercept) — the right answer comes back exactly."""
        beta = BETA_FDTD[0]
        seed = beta + 0.7 * np.pi / float(np.min(np.diff(LENGTHS)))
        f = beta_fit_from_phases(LENGTHS, _wrapped(beta, LENGTHS), seed_beta=seed)
        assert f.beta == pytest.approx(beta, rel=1e-9)
        assert f.resid < 1e-9
        assert f.branch_ambiguity < 0.45, "decidable, unlike the exactly-half case"

    def test_a_whole_branch_error_is_invisible_in_the_fit(self):
        """FAILURE MODE 2, and the reason section_mode_beta() constrains the
        LENGTHS rather than trusting any after-the-fact check.

        Seeded a full branch off, every length agrees on that wrong branch: the
        fit is clean to 1e-14, the margin reads 0.000, and β is 33% wrong.
        Nothing in the returned numbers can distinguish it from a good fit —
        which is why the guard has to be a precondition on the spacing."""
        beta = BETA_FDTD[0]
        one_branch = 2 * np.pi / float(np.min(np.diff(LENGTHS)))
        f = beta_fit_from_phases(
            LENGTHS, _wrapped(beta, LENGTHS), seed_beta=beta + one_branch)
        assert f.beta == pytest.approx(beta + one_branch, rel=1e-9), "lands a branch off"
        assert f.resid < 1e-9, "and the fit looks perfect"
        assert f.margin < 1e-6, "and the margin sees nothing"
        assert abs(f.beta / beta - 1) > 0.3, "while being 33% wrong"
        # The 2026-08-01 turn-consistency fix cannot see this either, and that is
        # the whole reason _BETA_MIN_BRANCH_FRACTION has to stay a precondition:
        # a whole-branch error is turn-CONSISTENT, so there is nothing to repair.
        assert not f.turns_shifted, "consistent turns, so the fix has no purchase"
        assert f.turn_ambiguity < 1e-9
        assert f.branch_ambiguity < 1e-9

    def test_needs_at_least_two_lengths(self):
        beta = BETA_FDTD[0]
        got, _, _ = beta_from_phases([2.0, 4.0], _wrapped(beta, [2.0, 4.0]), beta)
        assert got == pytest.approx(beta, rel=1e-9)

    def test_two_lengths_disable_the_residual_guard(self):
        """Why the default is three. Two points fit a line exactly, so the
        residual is identically zero whatever the data — guard 1 has nothing to
        detect. Deliberately fed a badly non-linear phase set here."""
        two = np.array([2.0, 4.0])
        nonsense = np.array([0.3, -2.9])
        f_two = beta_fit_from_phases(two, nonsense, BETA_MPB[0])
        assert f_two.resid < 1e-12, "two points can never show a residual"
        # With two points the BRANCH choice is near-arbitrary here too, not just
        # the residual blind — so the same configuration is doubly unsafe, and
        # branch_ambiguity is the thing that says so. Asserted because the
        # returned slope changed with the 2026-08-01 fix (the single turn is now
        # picked nearest the seed) and this test did not pin it.
        assert f_two.branch_ambiguity > 0.45

        three = np.array([2.0, 4.0, 6.0])
        beta = BETA_FDTD[0]
        good = _wrapped(beta, three)
        bent = good.copy()
        bent[1] += 0.5                                  # break linearity
        f_ok = beta_fit_from_phases(three, good, BETA_MPB[0])
        f_bent = beta_fit_from_phases(three, bent, BETA_MPB[0])
        assert f_ok.resid < 1e-9
        assert f_bent.resid > _BETA_FIT_RESIDUAL_MAX, "three points can"
        # And the residual guard now means ONE thing, so it needs a companion
        # that separates "genuinely curved" from "an integer turn was
        # undecidable". A 0.5 rad bend is curvature: the turns were never in
        # doubt, so turn_ambiguity stays low while the residual fires.
        assert f_bent.turn_ambiguity < 0.25, "curvature, not an undecidable turn"


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
        """Spacing is capped by β, not by how big the body happens to be, and
        all three lengths must still fit inside the body once each length's
        own dpml-padded solve_section() window is accounted for (notes.md
        Solvers item 15 — see TestSectionModeBetaDpmlClearance below)."""
        cap = self._max_safe_spacing(BETA_FDTD[0]) / 2.0
        dpml = 0.5
        for span in (4.0, 30.0, 300.0):
            usable_span = span - 2.0 * dpml
            spacing = min(usable_span / 3.0, cap)
            assert spacing <= self._max_safe_spacing(BETA_FDTD[0])
            assert 3.0 * spacing <= usable_span, "all three lengths must fit"
            assert 3.0 * spacing + 2.0 * dpml <= span, (
                "each length's own dpml-padded window must fit inside the span"
            )

    def test_a_safe_spacing_leaves_an_order_of_magnitude(self):
        """25% of β against a ~2% measured disagreement."""
        branch = 2 * np.pi / self._max_safe_spacing(BETA_FDTD[0])
        measured = abs(BETA_FDTD[0] - BETA_MPB[0])
        assert branch / measured > 10


class TestSectionModeBetaDpmlClearance:
    """Regression guard for notes.md Solvers item 15: the default fit-length
    spacing carried no ``-2*dpml`` clearance, unlike ``seed_span`` two lines
    above it in ``section_mode_beta()`` — so for a short uniform body the
    longest fit length's own dpml-padded ``solve_section()`` window overran
    the region, pushing the drawn cut into whatever taper/bend flanks it.

    Reproduced directly against a ``DirectionalCoupler(coupling_length=3.0)``:
    body span 2.50 µm, dpml 0.5 µm, β≈9.6 (this file's own measured value) —
    the old ``spacing = span/3.0`` gave lengths ``[0.833, 1.667, 2.500]``, a
    drawn half-window of 1.75 µm against a uniform half of 1.25 µm, i.e.
    dpml (0.5 µm) past each end.

    These exercise the real ``FDTDSolver.section_mode_beta()`` code path —
    not a reimplementation of its arithmetic — with ``section_mode_neff()``/
    ``solve_section()`` swapped for fakes that need no meep, only numpy.
    """

    SPAN = 2.50   # DirectionalCoupler(coupling_length=3.0)'s body span
    DPML = 0.5
    BETA = 9.6    # this file's own measured β, close enough to seed cleanly

    def _solver(self):
        from beef.solvers.fdtd import FDTDSolver
        return FDTDSolver(res=40, sim_width=5.0, dpml=self.DPML)

    def _wire_fakes(self, solver, betas=None, offset=0.0, seed_scale=1.0):
        """section_mode_neff() seeds β; solve_section() reports back exactly
        the phase a mode of that β would accumulate over the requested
        length, so beta_from_phases() recovers self.BETA to fitting
        precision — real orchestration, fake physics.

        Widened 2026-08-01 (defaults reproduce the original 1×1 behaviour, so
        the three callers above are unaffected):

        `betas`      one β per mode, so `T` comes back diagonal `(n_wl, n, n)`
                     instead of `(1, 1, 1)`. Needed to test that a raise names
                     the offending (wavelength, mode) cell at all.
        `offset`     the constant mode-matching / monitor-plane phase. Not
                     cosmetic: it sets where the true wrap counts fall relative
                     to a .5 rounding boundary, which is what decides whether
                     the pre-2026-08-01 per-length rounding broke.
        `seed_scale` multiplies the seed, i.e. injects a deliberate seed error.
                     The real FDTD-vs-MPB gap is ~2%, so `seed_scale=0.98` is
                     the production condition, not a stress test.
        """
        betas = np.atleast_1d(np.asarray(
            self.BETA if betas is None else betas, dtype=float))
        n = betas.size

        def fake_section_mode_neff(geometry, wavelengths, a, b, n_modes, **kw):
            wl = np.atleast_1d(np.asarray(wavelengths, dtype=float))
            return seed_scale * betas[None, :n_modes] * wl[:, None] / (2 * np.pi)

        solver.section_mode_neff = fake_section_mode_neff
        calls = []

        def fake_solve_section(geometry, wavelengths, x0, x1, **kw):
            calls.append((x0, x1))
            wl = np.atleast_1d(np.asarray(wavelengths, dtype=float))
            T = np.zeros((wl.size, n, n), dtype=complex)
            T[:, np.arange(n), np.arange(n)] = np.exp(
                1j * (betas[None, :] * (x1 - x0) + offset))
            return T

        solver.solve_section = fake_solve_section
        return calls

    def test_naive_formula_would_overrun_this_body(self):
        """Documents the bug: span/3.0 with no dpml term draws past the
        uniform region on exactly this device."""
        spacing = self.SPAN / 3.0
        longest = 3.0 * spacing
        assert longest / 2.0 + self.DPML > self.SPAN / 2.0

    def test_fit_lengths_stay_within_the_uniform_span(self):
        solver = self._solver()
        calls = self._wire_fakes(solver)
        solver.section_mode_beta(
            geometry=None, wavelengths=np.array([WL]),
            x0=-self.SPAN / 2.0, x1=self.SPAN / 2.0, n_modes=1,
        )
        assert len(calls) == 3, "default is three fit lengths"
        for x0, x1 in calls:
            half_drawn = (x1 - x0) / 2.0 + solver.dpml
            assert half_drawn <= self.SPAN / 2.0 + 1e-9, (
                f"call {(x0, x1)} draws a {half_drawn:.4f} half-window past "
                f"the {self.SPAN/2:.4f} uniform half"
            )

    def test_recovers_the_seeded_beta(self):
        """Not just legal lengths — the fit still has to work."""
        solver = self._solver()
        self._wire_fakes(solver)
        beta = solver.section_mode_beta(
            geometry=None, wavelengths=np.array([WL]),
            x0=-self.SPAN / 2.0, x1=self.SPAN / 2.0, n_modes=1,
        )
        assert beta[0, 0] == pytest.approx(self.BETA, rel=1e-6)

    # --- the 2026-08-01 turn-consistency fix, end to end -------------------

    #: Reproduces job 11589770's failure through the real `section_mode_beta`
    #: with no meep. seed_scale 0.98 is the ~2% FDTD-vs-MPB gap; the offset is
    #: what puts the true wrap counts near a .5 boundary. PINNED, not tuned:
    #: only ~4 of 260 (seed_scale, offset) combinations tripped the old code at
    #: these short lengths, so the premise assertion below is load-bearing.
    KNIFE_EDGE_OFFSET = -2 * np.pi * 0.525
    KNIFE_EDGE_SEED_SCALE = 0.98

    def _beta_at(self, betas=None, offset=0.0, seed_scale=1.0, wavelengths=None):
        solver = self._solver()
        calls = self._wire_fakes(
            solver, betas=betas, offset=offset, seed_scale=seed_scale)
        wl = np.array([WL]) if wavelengths is None else np.asarray(wavelengths)
        n = 1 if betas is None else np.atleast_1d(betas).size
        beta = solver.section_mode_beta(
            geometry=None, wavelengths=wl,
            x0=-self.SPAN / 2.0, x1=self.SPAN / 2.0, n_modes=n,
        )
        return beta, np.array([x1 - x0 for x0, x1 in calls])

    def test_a_two_percent_seed_error_no_longer_breaks_the_wrap_counts(self):
        """The bug that killed job 11589770, end to end and meep-free.

        Before 2026-08-01 this configuration raised
        ``phase-fit residual 2.094 rad exceeds 0.1`` — against production's
        2.101 on the real device. The premise assertion pins WHY: the old
        per-length rounding produced turn increments that disagree with each
        other, and a 2% seed error is all it took."""
        beta, lengths = self._beta_at(
            offset=self.KNIFE_EDGE_OFFSET, seed_scale=self.KNIFE_EDGE_SEED_SCALE)

        # Premise. Without this the test goes vacuous the moment anyone changes
        # section_mode_beta's default length selection.
        phi = np.angle(np.exp(1j * (self.BETA * lengths + self.KNIFE_EDGE_OFFSET)))
        seed = self.KNIFE_EDGE_SEED_SCALE * self.BETA
        naive = np.round((seed * lengths - phi) / (2 * np.pi))
        naive -= naive[0]
        chosen = beta_fit_from_phases(lengths, phi, seed).turns
        assert not np.array_equal(naive, chosen), (
            f"premise: at lengths {lengths.tolist()} the old per-length rounding "
            f"must be inconsistent, got naive={naive.tolist()} == chosen"
        )
        assert beta[0, 0] == pytest.approx(self.BETA, rel=1e-9)

    @pytest.mark.parametrize("seed_scale", [0.90, 0.95, 1.05, 1.10])
    def test_seed_error_is_tolerated_at_any_offset(self, seed_scale):
        """The non-knife-edge companion. The test above pins one adverse
        offset; this sweeps the offset out of the picture entirely, so it cannot
        pass by luck. 5-10% seed error is several times the real gap."""
        for offset in np.linspace(-np.pi, np.pi, 20):
            beta, _ = self._beta_at(offset=float(offset), seed_scale=seed_scale)
            assert beta[0, 0] == pytest.approx(self.BETA, rel=1e-6), (
                f"offset={offset:.4f} seed_scale={seed_scale}"
            )

    def test_a_non_uniform_cut_still_raises_and_names_the_cell(self):
        """The residual guard survives the fix — it just means one thing now.

        Also the reason `_wire_fakes` was widened past 1×1: job 11589770's
        message could report "2.101 rad" but not which of its 88 (wavelength,
        mode) cells, and that turned out to be the most informative fact about
        the failure."""
        solver = self._solver()
        betas = np.array([9.6, 9.5, 9.4])
        self._wire_fakes(solver, betas=betas)
        # Bend mode 1's phase at the middle length only: curvature, not a wrap.
        inner = solver.solve_section

        def bent(geometry, wavelengths, x0, x1, **kw):
            T = inner(geometry, wavelengths, x0, x1, **kw)
            if 0.9 < (x1 - x0) < 1.1:
                T[:, 1, 1] *= np.exp(1j * 0.5)
            return T

        solver.solve_section = bent
        with pytest.raises(ValueError, match=r"phase-fit residual") as excinfo:
            solver.section_mode_beta(
                geometry=None, wavelengths=np.array([1.50, WL, 1.60]),
                x0=-self.SPAN / 2.0, x1=self.SPAN / 2.0, n_modes=3,
            )
        assert "mode 1" in str(excinfo.value), (
            f"message must name the offending mode, got: {excinfo.value}"
        )
        assert "turn ambiguity" in str(excinfo.value), (
            "and must report turn ambiguity so curvature is distinguishable "
            "from an undecidable wrap count"
        )

    def test_a_body_too_short_for_any_dpml_clearance_raises(self):
        """span - 2*dpml <= 0 is refused before the fit-length spacing is
        even computed — seed_span's own preceding check (unchanged by this
        fix) already covers it, since it guards the identical condition."""
        solver = self._solver()
        self._wire_fakes(solver)
        with pytest.raises(ValueError, match="leaving no room"):
            solver.section_mode_beta(
                geometry=None, wavelengths=np.array([WL]),
                x0=-solver.dpml, x1=solver.dpml, n_modes=1,
            )


class TestTheMmi1x2TurnCountRegression:
    """mmi_1x2_validation.md gap 5 / job 11589770, as real numbers.

    ``arg(T_ii)`` read straight out of that job's own on-disk cache
    (``example_sim_cache/FDTDSolver_solve_section_MMI_{477592018d806434,
    7d8ebb646ca823d9,767b169eaa31dd8b}``), at λ=1.26 µm, mode 0 — the worst of
    its 88 (wavelength, mode) cells. Literals are 4 dp; the full-precision
    cache values give 13.877999 / 0.006445 / +1.9616%, so these are good to
    1e-4.

    This is the cell that was recorded for a day as a π sign flip in ``T_ii``.
    It is not one — see `test_the_answer_is_the_documented_fdtd_over_mpb_excess`.
    """

    L = [0.9232476, 1.8464952, 2.7697429]     # section_mode_beta's own [s, 2s, 3s]
    PHI = [3.0618, -2.9846, -2.7285]          # arg(T_ii), wrapped
    SEED = 13.611                             # MPB, via section_mode_neff
    WL = 1.26

    def test_the_consistent_turn_vector_is_chosen(self):
        f = beta_fit_from_phases(self.L, self.PHI, self.SEED)
        assert f.beta == pytest.approx(13.878, rel=1e-4)
        assert f.resid < 0.01
        assert list(f.turns) == [0, 3, 5]

    def test_the_old_per_length_rounding_is_not_what_comes_back(self):
        """Premise plus anti-regression, so this class cannot go vacuous if the
        fit is ever reverted or re-derived."""
        naive = np.round(
            (self.SEED * np.asarray(self.L) - np.asarray(self.PHI)) / (2 * np.pi))
        naive -= naive[0]
        assert list(naive) == [0, 2, 4], "premise: the old rounding is inconsistent"

        f = beta_fit_from_phases(self.L, self.PHI, self.SEED)
        # What the old code produced here: β 10.475, residual 2.1008 — 21x its
        # own 0.1 guard, which is how the job died.
        assert abs(f.beta - 10.475) > 3.0
        assert f.resid < _BETA_FIT_RESIDUAL_MAX / 10

    def test_the_answer_is_the_documented_fdtd_over_mpb_excess(self):
        """The discriminator against the withdrawn π-sign-flip reading, and the
        only thing that separates them — the two are algebraically identical on
        the residual (for three equally spaced points, "middle point displaced
        by −π" and "points 2 and 3 each displaced by −2π" both give
        δ·(+1/3, −2/3, +1/3)).

        The turn-count answer lands n_eff 2.783, +1.96% over MPB's 2.7285 — in
        family with the +2.233%/+2.361% this file's own BETA_FDTD/BETA_MPB
        record for the coupler. The sign-flip reading gave 2.101: −23%, wrong
        in sign as well as magnitude, and below any plausible n_eff for a
        2.0 µm × 160 nm Si core."""
        f = beta_fit_from_phases(self.L, self.PHI, self.SEED)
        n_eff = f.beta * self.WL / (2 * np.pi)
        assert n_eff == pytest.approx(2.783, abs=5e-3)
        assert 0.015 < f.seed_excess < 0.025
        coupler_excess = float(np.max(BETA_FDTD / BETA_MPB - 1.0))
        assert f.seed_excess < coupler_excess, "in family, and on the same side"

    def test_the_seed_error_is_reported_not_absorbed(self):
        """The seed error IS the measurement here, so the fix must surface it
        rather than quietly correcting for it."""
        f = beta_fit_from_phases(self.L, self.PHI, self.SEED)
        assert f.seed_excess == pytest.approx(0.0196, abs=1e-3)
        assert f.turns_shifted, "and that we left the seed's own rounding is visible"
        assert list(f.seed_turns) == [0, 2, 4]
        assert f.branch_ambiguity == pytest.approx(0.0377, abs=1e-3)
        # Matches the seed-free second-difference check on the same cache:
        # (φ₃ − 2φ₂ + φ₁)/2π is 1.0031, i.e. 0.0031 off an exact integer.
        assert f.turn_ambiguity == pytest.approx(0.0031, abs=5e-4)

    def test_the_three_tuple_facade_matches_the_fit(self):
        f = beta_fit_from_phases(self.L, self.PHI, self.SEED)
        assert beta_from_phases(self.L, self.PHI, self.SEED) == (
            f.beta, f.resid, f.margin)

    def test_lengths_out_of_order_give_the_same_answer(self):
        """`turns` comes back in the caller's order, so a descending-lengths
        caller must not get a permuted or different result."""
        f = beta_fit_from_phases(self.L, self.PHI, self.SEED)
        rev = beta_fit_from_phases(self.L[::-1], self.PHI[::-1], self.SEED)
        assert rev.beta == pytest.approx(f.beta, rel=1e-12)
        assert list(rev.turns) == list(f.turns)[::-1]


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
