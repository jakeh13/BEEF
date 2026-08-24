"""Step 1 of validation_markdown/resolved/multi_excitation_cascade_validation.md —
the equivalence harness that pins the 2x2's hand-rolled cascade before
steps 2-4 replace it.

No meep/meow import anywhere: the composition algebra is pure numpy, which
is the whole point of doing this before touching the examples.

THE ACCEPTANCE BAR: CLOSE AGREEMENT, NOT BIT-IDENTITY
-----------------------------------------------------
The plan originally asked for "bit-identical, not within tolerance". That was
relaxed to close agreement on 2026-08-20, for two measured reasons:

  1. It is unreachable anyway. `connect()` contracts with `@` (BLAS zgemm)
     while the example uses `np.einsum`; they sum in different orders, so
     they round differently. `np.einsum` vs `@` on the bare 14-element
     contraction is already not bit-identical (max |d| 3.8e-15), so this is
     inherited from numpy, not introduced by beef.network.

  2. It buys no defect-detection power. Every realistic migration bug
     mutation-tested against this harness (using the cross path in place of
     the bar, conjugating the input, dropping a mode) lands at O(1) RELATIVE
     error -- 2.02, 4.49, 2.58 respectively. A tolerance catches all of them.
     The only thing bit-identity additionally catches is a 1-ulp nudge, which
     is not a defect. Meanwhile it would make the suite hostage to the BLAS
     build, thread count and numpy version.

So `FLOAT_CHAIN_RTOL` below is THE GATE. Measured headroom: the real 2x2
device matrices agree to 4.0 eps (job 12070199, over the full 201x11 grid),
and synthetic data to 0.8-11 eps across N_MODES 4-64 -- roughly two orders
inside the gate, and ~13 orders inside anything a real bug produces.

`TestCircuitIsAlgebraicallyExact` is retained but is **NOT the bar**. It is a
FAILURE LOCALIZER: fed matrices whose every intermediate is exactly
representable in float64 (small integers), both routes are exact, so any
difference there is necessarily algebraic rather than numerical. If the gate
below fails, that class says whether the cause is the algebra or the
rounding. It is not brittle -- exact arithmetic has no rounding to vary --
but nothing should be blocked on it alone.

It also settles what `network/ops.py:148-152` asserts: the star product
really does collapse to ordinary transfer-matrix multiplication when
reflections are zero-filled.
"""
from __future__ import annotations

import numpy as np
import pytest

from beef.network.circuit import Circuit
from beef.network.convert import from_transfer_matrix, to_smatrix
from beef.solvers.hybrid import HybridSolver, MirrorSpec

# The shipped 2x2 configuration (mmi_2x2_hybrid.py:83-86, :147).
N_WL = 11
N_MODES = 14
WAVELENGTHS = np.linspace(1.530, 1.565, N_WL)

EPS = np.finfo(float).eps
# THE GATE. Measured worst case is ~11 eps synthetic / 4.0 eps on the real
# device; a dot product of length n carries an O(n*eps) rounding bound and
# this chain does two of them. Set well above the measurement and ~13 orders
# below the O(1) relative error any real algebraic mistake produces, so it
# discriminates defects from BLAS-vs-einsum rounding without being brittle.
FLOAT_CHAIN_RTOL = 1e-13


# ---------------------------------------------------------------------------
# The two compositions under comparison.
# ---------------------------------------------------------------------------
def compose_hand_rolled(M_in, M_body, M_out):
    """Byte-for-byte the example's own cascade (mmi_2x2_hybrid.py:304-306),
    for a single body-length slice. Returns the [n_wl] port coefficient."""
    tmp = np.einsum("wij,wjk->wik", M_body, M_in)
    return np.einsum("wij,wjk->wik", M_out, tmp)[:, 0, 0]


def compose_via_circuit(M_in, M_body, M_out, wavelengths):
    """Step 2's proposed replacement: three degenerate 2-port blocks joined
    by Circuit, read back out through to_smatrix(). Returns [n_wl]."""
    circuit = Circuit()
    circuit.add_block("inp", from_transfer_matrix(M_in, wavelengths, "src", "to_body"))
    circuit.add_block("body", from_transfer_matrix(M_body, wavelengths, "from_inp", "to_out"))
    circuit.add_block("out", from_transfer_matrix(M_out, wavelengths, "from_body", "port"))
    circuit.connect("inp", "to_body", "body", "from_inp")
    circuit.connect("body", "to_out", "out", "from_body")
    reduced = circuit.reduce(external_ports={"inp.src": "IN", "out.port": "OUT"})
    return to_smatrix(reduced, input_port="IN").coeff(port="OUT", in_mode=1, out_mode=1)


def bits_equal(a, b):
    """True only if every float64 lane is the identical bit pattern. Used by
    the localizer class only -- np.array_equal alone would treat -0.0 == 0.0."""
    a = np.ascontiguousarray(a)
    b = np.ascontiguousarray(b)
    return a.shape == b.shape and np.array_equal(a.view(np.uint8), b.view(np.uint8))


def max_rel(a, b):
    return float(np.max(np.abs(a - b) / np.abs(b)))


def assert_closely_aligns(got, expected, what):
    """THE GATE. Everything steps 2-4 produce is held to this."""
    rel = max_rel(got, expected)
    assert rel < FLOAT_CHAIN_RTOL, (
        f"{what}: max relative difference {rel:.3e} exceeds {FLOAT_CHAIN_RTOL:.0e} "
        f"({rel / EPS:.1f} eps). Rounding alone measures ~4 eps here, so this is "
        f"a real disagreement, not BLAS-vs-einsum noise."
    )
    return rel


# ---------------------------------------------------------------------------
# Matrix builders.
# ---------------------------------------------------------------------------
def exact_matrices(seed=20260820, n_modes=N_MODES, n_wl=N_WL):
    """Small-integer complex matrices. Every product is < 2^13 and every sum
    is over <= n_modes terms, so all intermediates land far inside float64's
    exact-integer range (2^53) and BOTH contraction orders are exact."""
    rng = np.random.default_rng(seed)

    def draw(*shape):
        return (rng.integers(-8, 9, size=shape).astype(float)
                + 1j * rng.integers(-8, 9, size=shape).astype(float))

    M_in = draw(n_wl, n_modes, 1)
    M_body = draw(n_wl, n_modes, n_modes)
    return M_in, M_body


def float_matrices(seed=20260820, n_modes=N_MODES, n_wl=N_WL):
    rng = np.random.default_rng(seed)

    def draw(*shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    return draw(n_wl, n_modes, 1), draw(n_wl, n_modes, n_modes)


def mode_signs(n_modes=N_MODES):
    """A parity vector of the shape mmi_2x2_hybrid.py:153 hard-codes. Its
    actual values are a measured property of the device; only +/-1 entries
    matter for the composition algebra, so the harness uses its own."""
    return np.array([1.0 if i % 3 else -1.0 for i in range(n_modes)])


# ---------------------------------------------------------------------------
class TestTheGate:
    """Close agreement on realistic data. This is what steps 2-4 must meet."""

    @pytest.mark.parametrize("n_modes", [4, 14, 32])
    def test_bar_path_closely_aligns(self, n_modes):
        M_in1, M_body = float_matrices(n_modes=n_modes)
        M_out_bar = M_in1.transpose(0, 2, 1)

        assert_closely_aligns(
            compose_via_circuit(M_in1, M_body, M_out_bar, WAVELENGTHS),
            compose_hand_rolled(M_in1, M_body, M_out_bar),
            f"bar path, N_MODES={n_modes}",
        )

    def test_cross_path_closely_aligns(self):
        """The cross path is the whole reason this plan exists: it pairs one
        arm's transposed matrix with the OTHER arm's excitation, which
        HybridSolver.solve()'s one-call contract cannot express."""
        M_in1, M_body = float_matrices()
        M_out_cross = (M_in1 * mode_signs()[None, :, None]).transpose(0, 2, 1)

        assert_closely_aligns(
            compose_via_circuit(M_in1, M_body, M_out_cross, WAVELENGTHS),
            compose_hand_rolled(M_in1, M_body, M_out_cross),
            "cross path",
        )

    def test_length_sweep_grid_closely_aligns(self):
        """The example composes a [n_L, n_wl, i, j] body against a
        [n_wl, i, j] input. Circuit has no length axis, so step 2 must loop --
        pin that the loop agrees over the whole grid, not just one slice."""
        rng = np.random.default_rng(31)
        n_L = 5
        M_in1, _ = float_matrices()
        M_body_L = (rng.standard_normal((n_L, N_WL, N_MODES, N_MODES))
                    + 1j * rng.standard_normal((n_L, N_WL, N_MODES, N_MODES)))
        M_out_bar = M_in1.transpose(0, 2, 1)

        tmp = np.einsum("Lwij,wjk->Lwik", M_body_L, M_in1)
        expected = np.einsum("wij,Lwjk->Lwik", M_out_bar, tmp)[:, :, 0, 0]
        got = np.stack([compose_via_circuit(M_in1, M_body_L[i], M_out_bar, WAVELENGTHS)
                        for i in range(n_L)])

        assert_closely_aligns(got, expected, "length-sweep grid")

    def test_transmission_closely_aligns(self):
        """The composed coefficient feeds |.|**2 -> T_bar/T_cross
        (mmi_2x2_hybrid.py:308-309). Squaring is where a sign or conjugation
        error would otherwise hide, so gate the transmission too."""
        M_in1, M_body = float_matrices()
        M_in0 = M_in1 * mode_signs()[None, :, None]

        for label, M_out in (("bar", M_in1.transpose(0, 2, 1)),
                             ("cross", M_in0.transpose(0, 2, 1))):
            assert_closely_aligns(
                np.abs(compose_via_circuit(M_in1, M_body, M_out, WAVELENGTHS)) ** 2,
                np.abs(compose_hand_rolled(M_in1, M_body, M_out)) ** 2,
                f"transmission, {label}",
            )

    def test_the_gate_actually_rejects_a_realistic_migration_bug(self):
        """Guards the gate from being vacuously loose. Mutation-measured
        2026-08-20: using the cross composition where the bar is expected
        lands at relative error 2.02 -- 13 orders outside the gate."""
        M_in1, M_body = float_matrices()
        M_out_bar = M_in1.transpose(0, 2, 1)
        M_out_cross = (M_in1 * mode_signs()[None, :, None]).transpose(0, 2, 1)

        wrong = compose_via_circuit(M_in1, M_body, M_out_cross, WAVELENGTHS)
        right = compose_hand_rolled(M_in1, M_body, M_out_bar)

        assert max_rel(wrong, right) > 1e-3
        with pytest.raises(AssertionError, match="real disagreement"):
            assert_closely_aligns(wrong, right, "deliberately wrong pairing")


class TestCircuitIsAlgebraicallyExact:
    """NOT the acceptance bar -- a failure localizer (see module docstring).

    Under exact arithmetic both routes are exact, so a difference here is
    necessarily algebraic rather than numerical. Retained because it costs
    nothing, is not brittle, and tells you which kind of failure you have.
    """

    def test_bar_path_is_bit_identical_under_exact_arithmetic(self):
        M_in1, M_body = exact_matrices()
        M_out_bar = M_in1.transpose(0, 2, 1)   # mirror reuse, saving (1)

        hand = compose_hand_rolled(M_in1, M_body, M_out_bar)
        circ = compose_via_circuit(M_in1, M_body, M_out_bar, WAVELENGTHS)

        assert np.all(hand.real == np.round(hand.real)), "test lost exactness"
        assert np.abs(hand).max() < 2 ** 53, "intermediates left the exact range"
        assert bits_equal(hand, circ)

    def test_cross_path_is_bit_identical_under_exact_arithmetic(self):
        M_in1, M_body = exact_matrices()
        M_in0 = M_in1 * mode_signs()[None, :, None]   # saving (2)

        hand = compose_hand_rolled(M_in1, M_body, M_in0.transpose(0, 2, 1))
        circ = compose_via_circuit(M_in1, M_body, M_in0.transpose(0, 2, 1), WAVELENGTHS)

        assert bits_equal(hand, circ)

    def test_bar_and_cross_are_genuinely_different_numbers(self):
        """Guards the two tests above from passing vacuously: if the parity
        vector were all +1 the cross path would just be the bar path and
        neither test would be exercising the thing it names."""
        M_in1, M_body = exact_matrices()
        bar = compose_hand_rolled(M_in1, M_body, M_in1.transpose(0, 2, 1))
        M_in0 = M_in1 * mode_signs()[None, :, None]
        cross = compose_hand_rolled(M_in1, M_body, M_in0.transpose(0, 2, 1))

        assert not np.allclose(bar, cross)

    def test_the_residual_is_the_contraction_not_the_star_product(self):
        """Why the gate's residual exists at all. from_transfer_matrix()
        zero-fills every self-term, so connect()'s systems are exactly
        I - 0 @ 0 = I, and solving against the identity is bit-exact.
        Whatever difference remains therefore comes from einsum-vs-BLAS
        contraction order, not from the Redheffer machinery -- so a future
        reader should not go looking for a bug in ops.py."""
        identity = np.broadcast_to(np.eye(N_MODES, dtype=complex), (N_WL, N_MODES, N_MODES))
        rng = np.random.default_rng(7)
        X = rng.standard_normal((N_WL, N_MODES, 1)) + 1j * rng.standard_normal((N_WL, N_MODES, 1))

        assert bits_equal(np.linalg.solve(identity, X), X)


class TestKnownBlindSpot:
    def test_only_the_cross_path_can_detect_a_transposed_body(self):
        """Measured 2026-08-20, and the reason the cross path has to be in
        this harness rather than the bar path alone.

        The bar path is `M_in1^T @ B @ M_in1` -- a quadratic form whose 1x1
        result equals its own transpose, so it returns the IDENTICAL scalar
        for B and B^T. A step-2 migration that transposed the body matrix
        would therefore sail past a bar-only check; the pre-existing
        tests/test_network.py::TestConnect::
        test_mmi_2x2_chain_matches_manual_einsum is exactly such a check.

        The cross path pairs two DIFFERENT arm vectors, breaks the symmetry,
        and catches it. Note the shipped 2x2 currently runs with
        CALIBRATE_BODY=True, which makes the body diagonal (hence symmetric)
        and the bug genuinely invisible in production -- but the
        uncalibrated body is not symmetric, so this is a live hazard for
        step 2, not a hypothetical."""
        M_in1, M_body = exact_matrices()
        M_body_T = M_body.transpose(0, 2, 1)
        M_in0 = M_in1 * mode_signs()[None, :, None]
        assert not bits_equal(M_body, M_body_T), "test needs a non-symmetric body"

        bar = compose_hand_rolled(M_in1, M_body, M_in1.transpose(0, 2, 1))
        bar_T = compose_hand_rolled(M_in1, M_body_T, M_in1.transpose(0, 2, 1))
        assert bits_equal(bar, bar_T), "bar path unexpectedly sees the transpose"

        cross = compose_hand_rolled(M_in1, M_body, M_in0.transpose(0, 2, 1))
        cross_T = compose_hand_rolled(M_in1, M_body_T, M_in0.transpose(0, 2, 1))
        assert max_rel(cross, cross_T) > 1e-3


# ---------------------------------------------------------------------------
# Step 2 — routing both excitations through HybridSolver.solve().
#
# tests/test_hybrid.py's FakeGeometry ignores **kwargs, so input_idx/
# output_idx do not move width_at/center_at and its mirror guardrail cannot
# tell the top arm from the bottom. Every cross-excitation question needs a
# geometry that actually branches, or the answer is an artifact -- measured
# the hard way 2026-08-20, when the flat geometry reported that a non-mirror
# pairing raised no error.
# ---------------------------------------------------------------------------
_SPLIT = 0.0        # left of this input_idx selects, right of it output_idx
_OFFSET = 0.75      # branch centres, as in the 2x2's +/-0.75 um arms


class TwoArmGeometry:
    """Two branches selected by input_idx / output_idx, MMI's own convention.
    Mirror-symmetric about x=0 for a FIXED branch, deliberately not when the
    two ends select different branches -- which is the cross path."""

    def __init__(self, width=0.45):
        self.width = width

    def _branch(self, x, **kwargs):
        return kwargs.get("input_idx" if x < _SPLIT else "output_idx", 1)

    def width_at(self, x, **kwargs):
        return self.width

    def center_at(self, x, **kwargs):
        return +_OFFSET if self._branch(x, **kwargs) == 1 else -_OFFSET


class _CountingSolver:
    """Scripted [n_out, n_in] matrix, counting real solves."""

    def __init__(self, matrix):
        self.matrix = np.asarray(matrix, dtype=complex)
        self.call_count = 0

    def solve_section(self, geometry, wavelengths, x0, x1, n_modes_in,
                      n_modes_out, stub_left=0.0, stub_right=0.0, **kwargs):
        self.call_count += 1
        return np.broadcast_to(
            self.matrix, (len(wavelengths),) + self.matrix.shape).copy()


_BOUNDARIES = [-18.0, -13.75, 13.75, 18.0]
_BAR_KWARGS = [{"input_idx": 1}, {}, {"output_idx": 1}]
_CROSS_KWARGS = [{"input_idx": 1}, {}, {"output_idx": 0}]


def _chain(n=4, seed=5):
    rng = np.random.default_rng(seed)

    def draw(*shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    return (_CountingSolver(draw(n, 1)),      # input
            _CountingSolver(draw(n, n)),      # body
            _CountingSolver(draw(1, n)))      # output


class TestStep2RunCount:
    """What routing both excitations through solve() actually costs. The
    plan's step 2 marks this ARGUED; these pin it at fake-solver level, which
    is the right level for a question about how many solves happen."""

    N = 4

    def _solver(self, chain, mirror):
        return HybridSolver(solvers=list(chain), boundaries=_BOUNDARIES,
                            interface_modes=[self.N, self.N],
                            mirror_pairs=({2: 0} if mirror else None))

    def test_the_two_arm_geometry_is_not_flat(self):
        """Guards every test below: if the geometry did not branch, the
        mirror guardrail would be untestable and all of this would pass
        vacuously."""
        g = TwoArmGeometry()
        assert g.center_at(-16, input_idx=1) == pytest.approx(+_OFFSET)
        assert g.center_at(+16, output_idx=1) == pytest.approx(+_OFFSET)
        assert g.center_at(+16, output_idx=0) == pytest.approx(-_OFFSET)

    def test_bar_path_reuses_the_output_as_a_mirror(self):
        chain = _chain(self.N)
        inp, body, out = chain
        self._solver(chain, mirror=True).solve(
            TwoArmGeometry(), WAVELENGTHS, port="bar",
            geometry_kwargs=_BAR_KWARGS, cache={})

        assert inp.call_count == 1
        assert out.call_count == 0, "bar output should be section 0's transpose"

    def test_cross_path_is_correctly_rejected_by_the_mirror_guardrail(self):
        """The plan predicts this: section 0 carries input_idx=1 (top) and
        section 2 carries output_idx=0 (bottom), which is diagonal rather
        than mirror-symmetric, so the declared pair is not valid for this
        excitation and the guardrail must say so."""
        chain = _chain(self.N)
        solver = self._solver(chain, mirror=True)
        geom = TwoArmGeometry()
        cache = {}
        solver.solve(geom, WAVELENGTHS, port="bar",
                     geometry_kwargs=_BAR_KWARGS, cache=cache)

        with pytest.raises(ValueError, match="center_at mismatch"):
            solver.solve(geom, WAVELENGTHS, port="cross",
                         geometry_kwargs=_CROSS_KWARGS, cache=cache)

    def test_a_shared_cache_reuses_the_input_section_across_excitations(self):
        """Step 2's central mechanism: both excitations excite the SAME top
        arm, so the input section's kwargs are identical between the two
        calls and it must be simulated once, not twice."""
        chain = _chain(self.N)
        inp, body, out = chain
        geom = TwoArmGeometry()
        cache = {}
        self._solver(chain, mirror=True).solve(
            geom, WAVELENGTHS, port="bar", geometry_kwargs=_BAR_KWARGS, cache=cache)
        after_bar = inp.call_count

        self._solver(chain, mirror=False).solve(
            geom, WAVELENGTHS, port="cross", geometry_kwargs=_CROSS_KWARGS, cache=cache)

        assert after_bar == 1
        assert inp.call_count == 1, "input section was re-solved for the 2nd excitation"

    def test_step_2_alone_costs_exactly_one_extra_solve(self):
        """The number step 3 has to recover. MEASURED 2026-08-20: two solves
        total (1 input + 1 output) against the hand-rolled path's one.

        This REFINES the plan's own estimate, which said porting would push
        "input-side FDTD runs from 1 to 2+ and add a real output-section
        simulation" -- i.e. 3+. The shared cache already recovers the
        input-side duplicate, so the real gap is one output-side solve."""
        chain = _chain(self.N)
        inp, body, out = chain
        geom = TwoArmGeometry()
        cache = {}

        self._solver(chain, mirror=True).solve(
            geom, WAVELENGTHS, port="bar", geometry_kwargs=_BAR_KWARGS, cache=cache)
        self._solver(chain, mirror=False).solve(
            geom, WAVELENGTHS, port="cross", geometry_kwargs=_CROSS_KWARGS, cache=cache)

        assert (inp.call_count, out.call_count) == (1, 1)
        assert inp.call_count + out.call_count == 2   # hand-rolled costs 1


# ---------------------------------------------------------------------------
# Step 3(a) — cross-excitation mirror declaration.
#
# What step 2 leaves on the table is exactly one output-side solve
# (TestStep2RunCount above). MirrorSpec is what recovers it: the declaration
# that section 2 under output_idx=0 is section 0 under input_idx=0 mirrored,
# and that section 0 under input_idx=0 is section 0 under input_idx=1 times a
# parity vector. Both halves are geometrically checked every solve() call.
# ---------------------------------------------------------------------------
def _signs(n):
    """A non-trivial parity vector — deliberately NOT all +1, or every test
    below would pass with the mode_signs machinery removed entirely."""
    return np.array([1.0 if i % 3 else -1.0 for i in range(n)])


class TestStep3CrossExcitationMirror:
    N = 4

    def _cross_spec(self, **over):
        kw = dict(section=0, partner_kwargs={"input_idx": 0},
                  mode_signs=_signs(self.N))
        kw.update(over)
        return MirrorSpec(**kw)

    def _solver(self, chain, mirror_pairs):
        return HybridSolver(solvers=list(chain), boundaries=_BOUNDARIES,
                            interface_modes=[self.N, self.N],
                            mirror_pairs=mirror_pairs)

    # -- the declaration is accepted where a bare int is rejected -----------
    def test_the_cross_pairing_a_bare_int_rejects_is_accepted_as_a_spec(self):
        """The exact pairing TestStep2RunCount's
        test_cross_path_is_correctly_rejected_by_the_mirror_guardrail shows
        the bare-int form refusing. Refusing it was CORRECT — section 2's
        bottom arm is not the x-mirror of section 0's top arm. MirrorSpec
        makes the true statement instead: it is the x-mirror of section 0's
        BOTTOM arm, which is the top arm's matrix times a parity vector."""
        chain = _chain(self.N)
        geom = TwoArmGeometry()

        with pytest.raises(ValueError, match="center_at mismatch"):
            self._solver(chain, {2: 0}).solve(
                geom, WAVELENGTHS, port="cross",
                geometry_kwargs=_CROSS_KWARGS, cache={})

        self._solver(chain, {2: self._cross_spec()}).solve(
            geom, WAVELENGTHS, port="cross",
            geometry_kwargs=_CROSS_KWARGS, cache={})

    # -- it recovers the run count -----------------------------------------
    def test_both_excitations_now_cost_one_input_solve_total(self):
        """THE POINT OF STEP 3. Step 2 measured 2 solves (1 input + 1
        output); the hand-rolled path costs 1. With the cross pair declared,
        both excitations come from that single input-side solve — so the
        manuscript's 'a single FDTD simulation' claim is now true through the
        supported API, not only through the hand-rolled einsum."""
        chain = _chain(self.N)
        inp, body, out = chain
        geom = TwoArmGeometry()
        cache = {}

        self._solver(chain, {2: 0}).solve(
            geom, WAVELENGTHS, port="bar", geometry_kwargs=_BAR_KWARGS, cache=cache)
        self._solver(chain, {2: self._cross_spec()}).solve(
            geom, WAVELENGTHS, port="cross", geometry_kwargs=_CROSS_KWARGS, cache=cache)

        assert inp.call_count == 1
        assert out.call_count == 0, "the output section was simulated after all"
        assert body.call_count == 1

    # -- it computes the right matrix --------------------------------------
    def test_the_composed_result_equals_the_hand_rolled_cascade(self):
        """Equivalence against the thing step 4 deletes, held to step 1's
        gate. Covers BOTH paths, because the bar path alone cannot see a
        transposed body (TestKnownBlindSpot)."""
        chain = _chain(self.N)
        inp, body, out = chain
        geom = TwoArmGeometry()
        cache = {}

        bar = self._solver(chain, {2: 0}).solve(
            geom, WAVELENGTHS, port="bar", geometry_kwargs=_BAR_KWARGS, cache=cache)
        cross = self._solver(chain, {2: self._cross_spec()}).solve(
            geom, WAVELENGTHS, port="cross", geometry_kwargs=_CROSS_KWARGS, cache=cache)

        M_in1 = np.broadcast_to(inp.matrix, (N_WL,) + inp.matrix.shape).copy()
        M_body = np.broadcast_to(body.matrix, (N_WL,) + body.matrix.shape).copy()
        M_in0 = M_in1 * _signs(self.N)[None, :, None]

        assert_closely_aligns(
            bar.coeff(port="bar", in_mode=1, out_mode=1),
            compose_hand_rolled(M_in1, M_body, M_in1.transpose(0, 2, 1)),
            "bar via MirrorSpec-equipped solve()",
        )
        assert_closely_aligns(
            cross.coeff(port="cross", in_mode=1, out_mode=1),
            compose_hand_rolled(M_in1, M_body, M_in0.transpose(0, 2, 1)),
            "cross via MirrorSpec",
        )

    def test_the_parity_vector_is_actually_applied(self):
        """Guards every test above from passing with mode_signs ignored: an
        all-+1 vector must give a genuinely different answer from the real
        one, or none of this is being exercised."""
        chain = _chain(self.N)
        geom = TwoArmGeometry()

        real = self._solver(chain, {2: self._cross_spec()}).solve(
            geom, WAVELENGTHS, port="c", geometry_kwargs=_CROSS_KWARGS, cache={})
        flat = self._solver(
            chain, {2: self._cross_spec(mode_signs=np.ones(self.N))}).solve(
            geom, WAVELENGTHS, port="c", geometry_kwargs=_CROSS_KWARGS, cache={})

        assert max_rel(real.coeff(port="c", in_mode=1, out_mode=1),
                       flat.coeff(port="c", in_mode=1, out_mode=1)) > 1e-3

    # -- partner_kwargs is not an unchecked escape hatch -------------------
    def test_a_partner_kwargs_naming_a_non_mirror_branch_is_rejected(self):
        """Without _validate_two_arm_symmetry, partner_kwargs would let any
        two unrelated branches be paired: the x-mirror check compares section
        k against section j UNDER partner_kwargs, so it never looks at what
        section j was actually solved with. Here partner_kwargs names the
        SAME arm the partner was solved on, so the 'other arm' claim is
        false — and the parity vector would be reinterpreting a matrix as
        itself."""
        chain = _chain(self.N)
        geom = TwoArmGeometry()

        with pytest.raises(ValueError, match="not reflections of each other"):
            self._solver(
                chain, {2: self._cross_spec(partner_kwargs={"input_idx": 1})}
            ).solve(geom, WAVELENGTHS, port="cross",
                    geometry_kwargs=_CROSS_KWARGS, cache={})

    def test_the_two_arm_check_runs_on_a_cache_hit_too(self):
        """Same requirement the bare-int guardrail already meets
        (tests/test_hybrid.py::test_mirror_guardrail_still_runs_on_a_cache_hit):
        a cross-excitation pair must re-validate per call, not once at
        construction — step 3's own note says so, and step 2 deliberately
        shares one cache across excitations, making the hit the common path."""
        chain = _chain(self.N)
        geom = TwoArmGeometry()
        cache = {}
        self._solver(chain, {2: self._cross_spec()}).solve(
            geom, WAVELENGTHS, port="cross", geometry_kwargs=_CROSS_KWARGS, cache=cache)

        bad = self._solver(
            chain, {2: self._cross_spec(partner_kwargs={"input_idx": 1})})
        with pytest.raises(ValueError, match="not reflections of each other"):
            bad.solve(geom, WAVELENGTHS, port="cross",
                      geometry_kwargs=_CROSS_KWARGS, cache=cache)

    # -- declaration-time validation ---------------------------------------
    def test_mode_signs_must_be_plus_or_minus_one(self):
        with pytest.raises(ValueError, match=r"\+1 or -1"):
            MirrorSpec(section=0, partner_kwargs={"input_idx": 0},
                       mode_signs=np.array([1.0, -1.0, 0.5, 1.0]))

    def test_mode_signs_without_partner_kwargs_is_refused(self):
        """A parity vector with no reinterpretation is meaningless, and
        silently ignoring it would be the worst outcome."""
        with pytest.raises(ValueError, match="without partner_kwargs"):
            MirrorSpec(section=0, mode_signs=np.array([1.0, -1.0]))

    def test_mode_signs_length_must_match_the_interface(self):
        chain = _chain(self.N)
        with pytest.raises(ValueError, match="one entry per mode"):
            self._solver(chain, {2: self._cross_spec(mode_signs=np.ones(self.N + 1))})

    def test_calibrate_kwargs_reaches_the_calibration_but_not_the_section(self):
        """Step 4's directional-coupler blocker, pinned.

        One geometry_kwargs dict per section serves two consumers: the section
        solver and the calibration source. The DC needs an interface monitor on
        the CALIBRATION but not on the EME body sweep -- and the difference is
        not cosmetic, because EMESolver.solve_section_length_sweep() funnels
        unknown **kwargs into geometry_kwargs, which is in the cache
        fingerprint (measured: 0c243abea80cb883 -> 76f7367470c61021, i.e. a
        forced re-solve of an already-cached body).

        So calibrate_kwargs must land on one call and not the other."""
        seen = {}

        class FakeFDTD:
            def section_mode_beta(self, geom, wl, x0, x1, n, lengths=None, **kw):
                seen["calibration"] = dict(kw)
                return np.full((len(wl), n), 1.0)

            def solve_section(self, *a, **k):
                return np.zeros((len(a[1]), 1, 1), dtype=complex)

        class FakeEME:
            def solve_section(self, *a, neff_override=None, **k):
                seen["section"] = dict(k)
                return np.zeros((len(a[1]), 1, 1), dtype=complex)

        h = HybridSolver([FakeFDTD(), FakeEME(), FakeFDTD()], [0, 1, 2, 3], [2, 2],
                         calibrate_body_neff=True,
                         calibrate_kwargs={"monitor_center": 0.0, "monitor_size": 6.0})
        h._calibrated_neff(TwoArmGeometry(), np.array([1.55]), 1, 2, {}, None, False)

        assert seen["calibration"]["monitor_center"] == 0.0
        assert seen["calibration"]["monitor_size"] == 6.0
        assert "section" not in seen, "calibrate_kwargs must not trigger a section solve"

    def test_calibrate_kwargs_defaults_to_the_historical_behaviour(self):
        """None must reproduce exactly what shipped before: only the
        monitor_center/monitor_size auto-extracted from the section's own
        kwargs, nothing more."""
        seen = {}

        class FakeFDTD:
            def section_mode_beta(self, geom, wl, x0, x1, n, lengths=None, **kw):
                seen.update(kw)
                return np.full((len(wl), n), 1.0)

            def solve_section(self, *a, **k):
                return np.zeros((len(a[1]), 1, 1), dtype=complex)

        class FakeEME:
            def solve_section(self, *a, neff_override=None, **k):
                return np.zeros((len(a[1]), 1, 1), dtype=complex)

        h = HybridSolver([FakeFDTD(), FakeEME(), FakeFDTD()], [0, 1, 2, 3], [2, 2],
                         calibrate_body_neff=True)
        h._calibrated_neff(TwoArmGeometry(), np.array([1.55]), 1, 2,
                           {"monitor_center": 1.0}, None, False)
        assert seen.get("monitor_center") == 1.0
        assert "monitor_size" not in seen

        seen.clear()
        h._calibrated_neff(TwoArmGeometry(), np.array([1.55]), 1, 2, {}, None, False)
        assert "monitor_center" not in seen and "monitor_size" not in seen

    def test_calibrate_kwargs_wins_over_the_auto_extracted_value(self):
        """Explicit beats inferred, per key -- otherwise a section that happens
        to carry a monitor for its own solver would silently override the
        calibration setting the caller asked for."""
        seen = {}

        class FakeFDTD:
            def section_mode_beta(self, geom, wl, x0, x1, n, lengths=None, **kw):
                seen.update(kw)
                return np.full((len(wl), n), 1.0)

            def solve_section(self, *a, **k):
                return np.zeros((len(a[1]), 1, 1), dtype=complex)

        class FakeEME:
            def solve_section(self, *a, neff_override=None, **k):
                return np.zeros((len(a[1]), 1, 1), dtype=complex)

        h = HybridSolver([FakeFDTD(), FakeEME(), FakeFDTD()], [0, 1, 2, 3], [2, 2],
                         calibrate_body_neff=True,
                         calibrate_kwargs={"monitor_center": 9.0})
        h._calibrated_neff(TwoArmGeometry(), np.array([1.55]), 1, 2,
                           {"monitor_center": 1.0}, None, False)
        assert seen["monitor_center"] == 9.0

    def test_a_bare_int_still_behaves_exactly_as_before(self):
        """MirrorSpec is additive. The historical form must be untouched —
        including .mirror_pairs staying a plain {int: int} map, which
        plot/geometry.py:569 reads."""
        chain = _chain(self.N)
        solver = self._solver(chain, {2: 0})
        assert solver.mirror_pairs == {2: 0}
        assert solver.mirror_specs[2] == MirrorSpec(section=0)
        assert solver.mirror_specs[2].signs_array is None
