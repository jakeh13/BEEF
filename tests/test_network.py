"""Fast unit tests for beef.network (SMatrixBlock, connect, close_loop, Circuit).

No meep/meow import anywhere — literal scripted np.ndarray matrices only,
following tests/test_hybrid.py's style. Two cross-checks anchor correctness
against real, already-proven physics/code rather than only against
self-consistent new math:
  - TestConnect::test_mmi_2x2_chain_matches_manual_einsum reproduces
    mmi_2x2_hybrid.py's hand-rolled cascade exactly.
  - The ring's exact closed-form formula is checked separately in
    tests/test_network_ring_crosscheck.py.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from beef.network.block import SMatrixBlock, check_passive, rename_ports
from beef.network.circuit import Circuit
from beef.network.convert import from_smatrix, from_transfer_matrix, to_smatrix
from beef.network.ops import close_loop, connect
from beef.results.smatrix import SMatrix

WAVELENGTHS = np.array([1.53, 1.55, 1.57])


def _random_block(ports, modes_per_port, wavelengths=WAVELENGTHS, seed=0):
    rng = np.random.default_rng(seed)
    N = sum(modes_per_port[p] for p in ports)
    M = len(wavelengths)
    data = rng.normal(size=(M, N, N)) + 1j * rng.normal(size=(M, N, N))
    return SMatrixBlock(wavelengths=wavelengths, ports=tuple(ports),
                         modes_per_port=dict(modes_per_port), data=data)


class TestSMatrixBlockConstruction:
    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError):
            SMatrixBlock(
                wavelengths=WAVELENGTHS, ports=("a", "b"),
                modes_per_port={"a": 1, "b": 2},
                data=np.zeros((len(WAVELENGTHS), 3, 4), dtype=complex),  # wrong N
            )

    def test_duplicate_port_raises(self):
        with pytest.raises(ValueError):
            SMatrixBlock(
                wavelengths=WAVELENGTHS, ports=("a", "a"),
                modes_per_port={"a": 1},
                data=np.zeros((len(WAVELENGTHS), 1, 1), dtype=complex),
            )

    def test_missing_modes_per_port_entry_raises(self):
        with pytest.raises(ValueError):
            SMatrixBlock(
                wavelengths=WAVELENGTHS, ports=("a", "b"),
                modes_per_port={"a": 1},  # missing "b"
                data=np.zeros((len(WAVELENGTHS), 1, 1), dtype=complex),
            )

    def test_port_slice_and_submatrix(self):
        block = _random_block(("a", "b"), {"a": 1, "b": 2})
        assert block.port_slice("a") == slice(0, 1)
        assert block.port_slice("b") == slice(1, 3)
        sub = block.submatrix("b", "a")
        assert sub.shape == (len(WAVELENGTHS), 2, 1)
        np.testing.assert_array_equal(sub, block.data[:, 1:3, 0:1])

    def test_flat_index(self):
        block = _random_block(("a", "b"), {"a": 1, "b": 2})
        assert block.flat_index("a", 1) == 0
        assert block.flat_index("b", 1) == 1
        assert block.flat_index("b", 2) == 2
        with pytest.raises(ValueError):
            block.flat_index("b", 3)


class TestRenamePorts:
    def test_renames_without_touching_data(self):
        block = _random_block(("a", "b"), {"a": 1, "b": 2})
        renamed = rename_ports(block, {"a": "x"})
        assert renamed.ports == ("x", "b")
        assert renamed.data is block.data

    def test_collision_raises(self):
        block = _random_block(("a", "b"), {"a": 1, "b": 2})
        with pytest.raises(ValueError):
            rename_ports(block, {"a": "b"})


class TestCheckPassive:
    def test_lossless_identity_passes(self):
        M = len(WAVELENGTHS)
        data = np.broadcast_to(np.eye(2, dtype=complex), (M, 2, 2)).copy()
        block = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("a",), modes_per_port={"a": 2}, data=data)
        assert check_passive(block, warn=False) is True

    def test_gain_fails_and_warns(self):
        M = len(WAVELENGTHS)
        data = np.broadcast_to(2.0 * np.eye(2, dtype=complex), (M, 2, 2)).copy()
        block = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("a",), modes_per_port={"a": 2}, data=data)
        with pytest.warns(UserWarning):
            assert check_passive(block) is False


class TestFromTransferMatrix:
    def test_matches_input_matrix_exactly(self):
        rng = np.random.default_rng(1)
        M, n_out, n_in = len(WAVELENGTHS), 3, 2
        matrix = rng.normal(size=(M, n_out, n_in)) + 1j * rng.normal(size=(M, n_out, n_in))
        block = from_transfer_matrix(matrix, WAVELENGTHS, port_in="left", port_out="right")
        np.testing.assert_array_equal(block.submatrix("right", "left"), matrix)
        assert np.count_nonzero(block.submatrix("left", "left")) == 0
        assert np.count_nonzero(block.submatrix("right", "right")) == 0
        assert np.count_nonzero(block.submatrix("left", "right")) == 0


class TestFromSmatrixToSmatrixRoundtrip:
    def test_roundtrip_bit_exact(self):
        rng = np.random.default_rng(2)
        M = len(WAVELENGTHS)
        coeffs_bar = rng.normal(size=M) + 1j * rng.normal(size=M)
        coeffs_cross = rng.normal(size=M) + 1j * rng.normal(size=M)
        smat = SMatrix(wavelengths=WAVELENGTHS)
        smat.set(coeffs_bar, port="bar")
        smat.set(coeffs_cross, port="cross")

        block = from_smatrix(smat, input_port="src")
        back = to_smatrix(block, input_port="src")

        np.testing.assert_array_equal(back.coeff(port="bar"), coeffs_bar)
        np.testing.assert_array_equal(back.coeff(port="cross"), coeffs_cross)

    def test_no_entries_raises(self):
        smat = SMatrix(wavelengths=WAVELENGTHS)
        smat.set(np.ones(len(WAVELENGTHS), dtype=complex), port="out")
        with pytest.raises(ValueError):
            from_smatrix(smat, input_port="src", output_ports=["missing"])


class TestConnect:
    def test_mismatched_mode_count_raises(self):
        a = _random_block(("x", "e1"), {"x": 1, "e1": 1}, seed=3)
        b = _random_block(("x", "e2"), {"x": 2, "e2": 1}, seed=4)
        with pytest.raises(ValueError, match="mode"):
            connect(a, "x", b, "x")

    def test_port_name_collision_raises(self):
        a = _random_block(("x", "shared"), {"x": 1, "shared": 1}, seed=5)
        b = _random_block(("x", "shared"), {"x": 1, "shared": 1}, seed=6)
        with pytest.raises(ValueError, match="collision"):
            connect(a, "x", b, "x")

    def test_wavelength_mismatch_raises(self):
        a = _random_block(("x", "e1"), {"x": 1, "e1": 1}, wavelengths=np.array([1.5, 1.6]), seed=7)
        b = _random_block(("x", "e2"), {"x": 1, "e2": 1}, wavelengths=np.array([1.5, 1.55, 1.6]), seed=8)
        with pytest.raises(ValueError, match="wavelength"):
            connect(a, "x", b, "x")

    def test_unknown_port_raises(self):
        a = _random_block(("x", "e1"), {"x": 1, "e1": 1}, seed=9)
        b = _random_block(("x", "e2"), {"x": 1, "e2": 1}, seed=10)
        with pytest.raises(ValueError):
            connect(a, "nope", b, "x")

    def test_directional_chain_reduces_to_transfer_matrix_product(self):
        """When both blocks come from from_transfer_matrix() (zero self-terms),
        connect() must reduce EXACTLY to ordinary matrix multiplication —
        the property the mmi_2x2 cross-check below relies on."""
        rng = np.random.default_rng(11)
        M = len(WAVELENGTHS)
        mat_a = rng.normal(size=(M, 2, 1)) + 1j * rng.normal(size=(M, 2, 1))
        mat_b = rng.normal(size=(M, 1, 2)) + 1j * rng.normal(size=(M, 1, 2))
        block_a = from_transfer_matrix(mat_a, WAVELENGTHS, "in", "mid")
        block_b = from_transfer_matrix(mat_b, WAVELENGTHS, "mid2", "out")
        composed = connect(block_a, "mid", block_b, "mid2")
        expected = np.einsum("mij,mjk->mik", mat_b, mat_a)
        np.testing.assert_allclose(composed.submatrix("out", "in"), expected, atol=1e-12)

    def test_hand_computed_small_connect(self):
        """1-mode connect against a manually-derived scalar star product."""
        M = len(WAVELENGTHS)
        S11, S12, S21, S22 = 0.3 + 0.1j, 0.9, 0.9, 0.2 - 0.2j
        T11, T12, T21, T22 = 0.4, 0.85, 0.85, 0.1 + 0.3j
        data_a = np.broadcast_to(np.array([[S11, S12], [S21, S22]], dtype=complex), (M, 2, 2)).copy()
        data_b = np.broadcast_to(np.array([[T11, T12], [T21, T22]], dtype=complex), (M, 2, 2)).copy()
        block_a = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("e1", "x"), modes_per_port={"e1": 1, "x": 1}, data=data_a)
        block_b = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("x", "e2"), modes_per_port={"x": 1, "e2": 1}, data=data_b)

        composed = connect(block_a, "x", block_b, "x")

        # scalar Redheffer star product, hand-derived
        m1 = 1 - T11 * S22
        expected_e1e1 = S11 + S12 * (T11 / m1) * S21
        expected_e1e2 = S12 * (T12 / m1)
        np.testing.assert_allclose(composed.submatrix("e1", "e1")[:, 0, 0], expected_e1e1, atol=1e-12)
        np.testing.assert_allclose(composed.submatrix("e2", "e1")[:, 0, 0], expected_e1e2, atol=1e-12)

    def test_mmi_2x2_chain_matches_manual_einsum(self):
        """Reproduces mmi_2x2_hybrid.py's hand-rolled 3-block cascade
        (M_in1 -> M_body -> M_out_bar, where M_out_bar is M_in1's mirror
        transpose) via connect(), and checks it matches the example's own
        manual np.einsum composition bit-for-bit (to floating-point noise)."""
        rng = np.random.default_rng(12)
        M, N_MODES = len(WAVELENGTHS), 4
        M_in1 = rng.normal(size=(M, N_MODES, 1)) + 1j * rng.normal(size=(M, N_MODES, 1))
        M_body = rng.normal(size=(M, N_MODES, N_MODES)) + 1j * rng.normal(size=(M, N_MODES, N_MODES))
        M_out_bar = M_in1.transpose(0, 2, 1)  # mirror-symmetry reuse, as the example does

        # -- expected, exactly mmi_2x2_hybrid.py's own computation --
        tmp = np.einsum("wij,wjk->wik", M_body, M_in1)
        expected_bar = np.einsum("wij,wjk->wik", M_out_bar, tmp)

        # -- composed via the network primitives --
        b_in1 = from_transfer_matrix(M_in1, WAVELENGTHS, "src", "mid")
        b_body = from_transfer_matrix(M_body, WAVELENGTHS, "mid_in", "mid_out")
        b_out = from_transfer_matrix(M_out_bar, WAVELENGTHS, "out_mid", "bar")
        c1 = connect(b_in1, "mid", b_body, "mid_in")
        c2 = connect(c1, "mid_out", b_out, "out_mid")

        np.testing.assert_allclose(c2.submatrix("bar", "src"), expected_bar, atol=1e-12)

    def test_near_singular_warns_and_singular_raises(self):
        M = len(WAVELENGTHS)
        # A near-critical-coupling-style pair: S_a,xx * S_b,xx -> 1
        data_a = np.broadcast_to(np.array([[0.1, 0.9], [0.9, 0.9999999]], dtype=complex), (M, 2, 2)).copy()
        data_b = np.broadcast_to(np.array([[0.9999999, 0.9], [0.9, 0.1]], dtype=complex), (M, 2, 2)).copy()
        block_a = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("e1", "x"), modes_per_port={"e1": 1, "x": 1}, data=data_a)
        block_b = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("x", "e2"), modes_per_port={"x": 1, "e2": 1}, data=data_b)
        with pytest.warns(UserWarning):
            connect(block_a, "x", block_b, "x")

        data_a2 = np.broadcast_to(np.array([[0.1, 0.9], [0.9, 1.0]], dtype=complex), (M, 2, 2)).copy()
        data_b2 = np.broadcast_to(np.array([[1.0, 0.9], [0.9, 0.1]], dtype=complex), (M, 2, 2)).copy()
        block_a2 = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("e1", "x"), modes_per_port={"e1": 1, "x": 1}, data=data_a2)
        block_b2 = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("x", "e2"), modes_per_port={"x": 1, "e2": 1}, data=data_b2)
        with pytest.raises(ValueError, match="near-singular"):
            connect(block_a2, "x", block_b2, "x")


class TestCloseLoop:
    def test_different_ports_raises_not_implemented(self):
        block = _random_block(("p", "q", "e"), {"p": 1, "q": 1, "e": 1}, seed=13)
        with pytest.raises(NotImplementedError):
            close_loop(block, "p", "q")

    def test_unknown_port_raises(self):
        block = _random_block(("p", "e"), {"p": 1, "e": 1}, seed=14)
        with pytest.raises(ValueError):
            close_loop(block, "nope", "nope")

    def test_multimode_matches_truncated_neumann_series(self):
        """n=2 loop port, strongly under-coupled so the geometric series
        converges fast — cross-checks close_loop()'s matrix inverse against
        a direct truncated-series summation computed independently."""
        rng = np.random.default_rng(15)
        M = len(WAVELENGTHS)
        n_E, n_p = 1, 2
        N = n_E + n_p
        raw = rng.normal(size=(M, N, N)) + 1j * rng.normal(size=(M, N, N))
        # Scale the p<->p block down so ||S_pp @ loop|| << 1 (fast convergence).
        raw[:, n_E:, n_E:] *= 0.1
        block = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("e", "p"), modes_per_port={"e": n_E, "p": n_p}, data=raw)
        loop = 0.2 * (rng.normal(size=(M, n_p, n_p)) + 1j * rng.normal(size=(M, n_p, n_p)))

        result = close_loop(block, "p", "p", loop=loop)

        S_EE = block.submatrix("e", "e")
        S_Ep = block.submatrix("e", "p")
        S_pE = block.submatrix("p", "e")
        S_pp = block.submatrix("p", "p")
        acc = S_EE.copy()
        term = S_pE.copy()  # (S_pp @ loop)^0 @ S_pE
        for _ in range(50):
            acc = acc + S_Ep @ loop @ term
            term = (S_pp @ loop) @ term
        np.testing.assert_allclose(result.submatrix("e", "e"), acc, atol=1e-9)

    def test_near_singular_warns_and_exact_singular_raises(self):
        M = len(WAVELENGTHS)
        t, k1, k2 = 0.95, 0.3, 0.3
        data = np.zeros((M, 2, 2), dtype=complex)
        data[:, 0, 0] = t
        data[:, 0, 1] = k2
        data[:, 1, 0] = k1
        data[:, 1, 1] = 0.9999999
        block = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("bus", "ring"), modes_per_port={"bus": 1, "ring": 1}, data=data)
        with pytest.warns(UserWarning):
            close_loop(block, "ring", "ring", loop=np.full(M, 0.9999999 + 0j))

        data2 = data.copy()
        data2[:, 1, 1] = 1.0
        block2 = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("bus", "ring"), modes_per_port={"bus": 1, "ring": 1}, data=data2)
        with pytest.raises(ValueError, match="near-singular"):
            close_loop(block2, "ring", "ring", loop=np.ones(M, dtype=complex))


class TestCircuit:
    def test_pairwise_chain_matches_direct_connect(self):
        rng = np.random.default_rng(16)
        M, N_MODES = len(WAVELENGTHS), 3
        M_in1 = rng.normal(size=(M, N_MODES, 1)) + 1j * rng.normal(size=(M, N_MODES, 1))
        M_body = rng.normal(size=(M, N_MODES, N_MODES)) + 1j * rng.normal(size=(M, N_MODES, N_MODES))
        M_out_bar = M_in1.transpose(0, 2, 1)

        b_in1 = from_transfer_matrix(M_in1, WAVELENGTHS, "src", "mid")
        b_body = from_transfer_matrix(M_body, WAVELENGTHS, "mid_in", "mid_out")
        b_out = from_transfer_matrix(M_out_bar, WAVELENGTHS, "out_mid", "bar")
        direct = connect(connect(b_in1, "mid", b_body, "mid_in"), "mid_out", b_out, "out_mid")

        circuit = Circuit()
        circuit.add_block("in1", b_in1)
        circuit.add_block("body", b_body)
        circuit.add_block("out", b_out)
        circuit.connect("in1", "mid", "body", "mid_in")
        circuit.connect("body", "mid_out", "out", "out_mid")
        reduced = circuit.reduce(external_ports={"in1.src": "src", "out.bar": "bar"})

        np.testing.assert_allclose(reduced.submatrix("bar", "src"), direct.submatrix("bar", "src"), atol=1e-12)

    def test_ring_via_circuit_matches_direct_close_loop(self):
        M = len(WAVELENGTHS)
        t, k1, k2, t_ring = 0.95 + 0j, 0.3 + 0j, 0.3 + 0j, 0.9 + 0j
        data = np.zeros((M, 2, 2), dtype=complex)
        data[:, 0, 0] = t
        data[:, 0, 1] = k2
        data[:, 1, 0] = k1
        data[:, 1, 1] = t_ring
        coupler = SMatrixBlock(wavelengths=WAVELENGTHS, ports=("bus", "ring"), modes_per_port={"bus": 1, "ring": 1}, data=data)
        loop = np.full(M, 0.8 + 0.1j)

        direct = close_loop(coupler, "ring", "ring", loop=loop)

        circuit = Circuit()
        circuit.add_block("coupler", coupler)
        circuit.close_loop("coupler", "ring", "ring", loop=loop)
        reduced = circuit.reduce()

        np.testing.assert_allclose(
            reduced.submatrix("coupler.bus", "coupler.bus"),
            direct.submatrix("bus", "bus"), atol=1e-12,
        )

    def test_disconnected_components_raise(self):
        a = _random_block(("p",), {"p": 1}, seed=17)
        b = _random_block(("p",), {"p": 1}, seed=18)
        circuit = Circuit()
        circuit.add_block("a", a)
        circuit.add_block("b", b)
        with pytest.raises(ValueError, match="disconnected"):
            circuit.reduce()

    def test_duplicate_name_raises(self):
        a = _random_block(("p",), {"p": 1}, seed=19)
        circuit = Circuit()
        circuit.add_block("a", a)
        with pytest.raises(ValueError):
            circuit.add_block("a", a)

    def test_unregistered_name_raises(self):
        a = _random_block(("p", "q"), {"p": 1, "q": 1}, seed=20)
        circuit = Circuit()
        circuit.add_block("a", a)
        with pytest.raises(ValueError):
            circuit.connect("a", "p", "missing", "p")
