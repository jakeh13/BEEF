"""close_loop_pair(): paired elimination of two ports fed into each other.

The two load-bearing tests here are:

  * test_matches_branch_a_when_one_direction_is_dead -- with loop_ab = 0 the
    paired formula must reduce EXACTLY to close_loop()'s Branch A, which is
    itself already proven against the ring resonator's closed-form all-pass
    formula (tests/test_network_ring_crosscheck.py). That chains the new
    function onto an existing proof rather than asserting it in isolation.

  * test_matches_brute_force_linear_solve -- for random blocks and random
    loops, the reduced block must reproduce what you get by writing down the
    full amplitude system (block relations + both loop constraints) and
    solving it directly. This validates the algebra without needing any
    closed form, and it is the test that would catch a transposed or
    misplaced sub-block in the 2n-by-2n system.
"""
import numpy as np
import pytest

from beef.network import close_loop, close_loop_pair
from beef.network.block import SMatrixBlock


def _block(ports, modes, wavelengths, seed):
    rng = np.random.default_rng(seed)
    n = sum(modes[p] for p in ports)
    M = len(wavelengths)
    data = (rng.normal(size=(M, n, n)) + 1j * rng.normal(size=(M, n, n))) * 0.3
    return SMatrixBlock(wavelengths=np.asarray(wavelengths), ports=tuple(ports),
                        modes_per_port=dict(modes), data=data)


def _brute_force(block, port_a, port_b, R, F):
    """Reduced S_EE from a DIRECT solve of the raw amplitude system.

    Deliberately does NOT reuse close_loop_pair()'s 2n-by-2n derivation -- that
    would make this test circular, validating the solve while assuming the
    algebra. Instead it writes down only the primitive statements:

        b   = S a                (the block itself, N equations)
        a_A = R b_B              (loop, nA equations)
        a_B = F b_A              (loop, nB equations)
        a_E = the driving vector (nE equations)

    which is 2N equations in the 2N unknowns [a; b], solved as one system. Any
    transposed or misplaced sub-block in the implementation shows up here.
    """
    E = [q for q in block.ports if q not in (port_a, port_b)]
    idx = {p_: np.arange(block.port_slice(p_).start, block.port_slice(p_).stop)
           for p_ in block.ports}
    iE = np.concatenate([idx[q] for q in E]) if E else np.array([], dtype=int)
    iA, iB = idx[port_a], idx[port_b]
    N = block.data.shape[1]
    nE = len(iE)
    M = len(block.wavelengths)
    out = np.zeros((M, nE, nE), dtype=complex)

    for m in range(M):
        S = block.data[m]
        for col in range(nE):
            # unknowns: x = [a (N), b (N)]
            A = np.zeros((2 * N, 2 * N), dtype=complex)
            rhs = np.zeros(2 * N, dtype=complex)
            r = 0
            # b - S a = 0
            A[r:r + N, :N] = -S
            A[r:r + N, N:] = np.eye(N)
            r += N
            # a_A - R b_B = 0
            for k, ia in enumerate(iA):
                A[r, ia] = 1.0
                A[r, N + iB] = -R[m][k, :]
                r += 1
            # a_B - F b_A = 0
            for k, ib in enumerate(iB):
                A[r, ib] = 1.0
                A[r, N + iA] = -F[m][k, :]
                r += 1
            # a_E = e_col
            for k, ie in enumerate(iE):
                A[r, ie] = 1.0
                rhs[r] = 1.0 if k == col else 0.0
                r += 1
            assert r == 2 * N, (r, 2 * N)
            x = np.linalg.solve(A, rhs)
            b = x[N:]
            out[m, :, col] = b[iE]
    return out


def _loops(n, M, seed, dead_ab=False):
    rng = np.random.default_rng(seed)
    R = (rng.normal(size=(M, n, n)) + 1j * rng.normal(size=(M, n, n))) * 0.2
    F = np.zeros((M, n, n), dtype=complex) if dead_ab else (
        (rng.normal(size=(M, n, n)) + 1j * rng.normal(size=(M, n, n))) * 0.2)
    return R, F


# ---------------------------------------------------------------------------
# The two substantive tests
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n_ring", [1, 2])
def test_matches_branch_a_when_one_direction_is_dead(n_ring):
    """loop_ab = 0 must reproduce close_loop()'s Branch A exactly.

    Constructed so the identification is unambiguous: build a [E, p] block,
    take its Branch A closure, then build a [E, A, B] block whose sub-blocks
    place S_pp at S_BA, S_Ep at S_EA and S_pE at S_BE, with the A-side
    self/cross terms zeroed. The two reductions must agree to machine
    precision.
    """
    wl = np.array([1.50, 1.55, 1.60])
    M, nE = len(wl), 2
    rng = np.random.default_rng(11)
    rnd = lambda *shape: (rng.normal(size=shape) + 1j * rng.normal(size=shape)) * 0.3
    S_EE, S_Ep, S_pE, S_pp = rnd(M, nE, nE), rnd(M, nE, n_ring), rnd(M, n_ring, nE), rnd(M, n_ring, n_ring)
    P = rnd(M, n_ring, n_ring) * 0.5

    # [E, p] with a self-loop -> Branch A
    n1 = nE + n_ring
    d1 = np.zeros((M, n1, n1), dtype=complex)
    d1[:, :nE, :nE], d1[:, :nE, nE:] = S_EE, S_Ep
    d1[:, nE:, :nE], d1[:, nE:, nE:] = S_pE, S_pp
    blk1 = SMatrixBlock(wavelengths=wl, ports=("E", "p"),
                        modes_per_port={"E": nE, "p": n_ring}, data=d1)
    ref = close_loop(blk1, "p", "p", loop=P)

    # [E, A, B] with only the B->A direction alive
    n2 = nE + 2 * n_ring
    d2 = np.zeros((M, n2, n2), dtype=complex)
    sE, sA, sB = slice(0, nE), slice(nE, nE + n_ring), slice(nE + n_ring, n2)
    d2[:, sE, sE] = S_EE
    d2[:, sE, sA] = S_Ep          # S_EA <- S_Ep
    d2[:, sB, sE] = S_pE          # S_BE <- S_pE
    d2[:, sB, sA] = S_pp          # S_BA <- S_pp
    blk2 = SMatrixBlock(wavelengths=wl, ports=("E", "A", "B"),
                        modes_per_port={"E": nE, "A": n_ring, "B": n_ring}, data=d2)
    got = close_loop_pair(blk2, "A", "B", loop_ba=P, loop_ab=0.0)

    assert got.ports == ref.ports == ("E",)
    assert np.allclose(got.data, ref.data, atol=1e-12), (
        f"max |diff| = {np.max(np.abs(got.data - ref.data)):.3e}")


@pytest.mark.parametrize("modes,seed", [
    ({"E1": 1, "E2": 1, "A": 1, "B": 1}, 1),
    ({"E1": 2, "A": 2, "B": 2}, 2),
    ({"E1": 1, "E2": 2, "A": 2, "B": 2}, 3),
])
def test_matches_brute_force_linear_solve(modes, seed):
    """The reduced block must equal a direct solve of the full amplitude
    system, for random S and random (both-direction) loops."""
    wl = np.array([1.50, 1.55, 1.60, 1.60 + 1e-3])
    ports = list(modes)
    blk = _block(ports, modes, wl, seed)
    n = modes["A"]
    R, F = _loops(n, len(wl), seed + 100)
    got = close_loop_pair(blk, "A", "B", loop_ba=R, loop_ab=F)
    want = _brute_force(blk, "A", "B", R, F)
    assert np.allclose(got.data, want, atol=1e-11), (
        f"max |diff| = {np.max(np.abs(got.data - want)):.3e}")


# ---------------------------------------------------------------------------
# Structure / contract
# ---------------------------------------------------------------------------
def test_removes_both_ports_and_keeps_the_rest_in_order():
    wl = np.array([1.55])
    modes = {"bus_L": 1, "ring_L": 1, "bus_R": 1, "ring_R": 1}
    blk = _block(list(modes), modes, wl, 7)
    out = close_loop_pair(blk, "ring_R", "ring_L", loop_ba=0.5, loop_ab=0.5)
    assert out.ports == ("bus_L", "bus_R")
    assert out.modes_per_port == {"bus_L": 1, "bus_R": 1}
    assert out.data.shape == (1, 2, 2)


def test_loop_defaults_to_identity_and_accepts_per_wavelength_scalars():
    """The ring's own `P` is a [M] array, so that form must work."""
    wl = np.array([1.50, 1.55, 1.60])
    modes = {"E": 1, "A": 1, "B": 1}
    blk = _block(list(modes), modes, wl, 5)
    ident = close_loop_pair(blk, "A", "B")                      # both None -> I
    ones = close_loop_pair(blk, "A", "B", loop_ba=1.0, loop_ab=1.0)
    assert np.allclose(ident.data, ones.data)
    per_wl = close_loop_pair(blk, "A", "B",
                             loop_ba=np.array([0.3, 0.4, 0.5]),
                             loop_ab=np.array([0.3, 0.4, 0.5]))
    assert per_wl.data.shape == (3, 1, 1)
    assert not np.allclose(per_wl.data, ident.data)


def test_rejects_same_port_and_points_at_close_loop():
    wl = np.array([1.55])
    modes = {"E": 1, "A": 1}
    blk = _block(list(modes), modes, wl, 3)
    with pytest.raises(ValueError, match="must be distinct"):
        close_loop_pair(blk, "A", "A")


def test_rejects_unknown_port_and_mode_count_mismatch():
    wl = np.array([1.55])
    modes = {"E": 1, "A": 1, "B": 2}
    blk = _block(list(modes), modes, wl, 4)
    with pytest.raises(ValueError, match="not in block"):
        close_loop_pair(blk, "A", "nope")
    with pytest.raises(ValueError, match="mode"):
        close_loop_pair(blk, "A", "B")


def test_flags_a_singular_loop_system():
    """A pair closure driven to resonance must complain rather than return
    garbage — the same contract close_loop() has."""
    wl = np.array([1.55])
    nE = 1
    d = np.zeros((1, 3, 3), dtype=complex)
    d[0, 0, 0] = 0.1
    d[0, 2, 1] = 1.0          # S_BA = 1
    d[0, 0, 1] = 0.5          # S_EA
    d[0, 2, 0] = 0.5          # S_BE
    blk = SMatrixBlock(wavelengths=wl, ports=("E", "A", "B"),
                       modes_per_port={"E": nE, "A": 1, "B": 1}, data=d)
    # S_BA @ R = 1 -> I - S_BA R singular in the B row, F dead
    with pytest.raises(ValueError):
        close_loop_pair(blk, "A", "B", loop_ba=1.0, loop_ab=0.0, cond_error=1e6)
