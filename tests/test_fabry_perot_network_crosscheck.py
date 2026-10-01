"""beef.network.Circuit against a hand-derived Fabry-Perot-with-DC-tap
closed form, verified BEFORE any real device numbers exist.

Topology (that doc's diagram): a finite Fabry-Perot cavity between two
1-port mirrors (the HCG), tapped in the middle by a genuine 4-port
directional coupler -- NOT a closed ring (contrast
tests/test_network_ring_crosscheck.py, which validates close_loop() against
a ring's own closed-form all-pass formula). This is beef.network's first
exercise of `Circuit.add_block()` composing more than one independently-
solved physical device (validation doc A2) and its first 1-port
self-reflecting block (A3) -- both previously unexercised by any example
or test in this repo.

Closed form, hand-derived and independently verified against a brute-force
4-port linear solve (scratch check, not reproduced here) before being
trusted: for a symmetric cavity (both mirrors identical reflection r, both
connecting waveguides identical single-pass phase p, DC through/cross
amplitudes t/k, all at one wavelength),

    rho = r * p**2                        # one mirror + round-trip waveguide phase
    D = 1 - rho**2 * t**2
    self_reflection = k**2 * rho / D      # tap port -> itself
    transmission    = t * (1 + rho**2*(k**2 - t**2)) / D   # tap port -> the other tap port

`self_reflection` is the coefficient into the SAME tap port (e.g.
input_reflection -> input_reflection); `transmission` is input_reflection
<-> output (equal both ways by reciprocity). This is the direct analogue of
the ring's own `t + k1*k2*P/(1 - t_ring*P)` formula, generalized to two
independent loop-closures (one per mirror) on a genuine 4-port block
instead of one closed ring port.
"""
from __future__ import annotations

import numpy as np

from beef.network import Circuit, connect
from beef.network.block import SMatrixBlock

WAVELENGTHS = np.linspace(1.53, 1.57, 5)


def _dc_block(t, k, wavelengths=WAVELENGTHS):
    """Synthetic 4-port directional coupler, ports (L, R, I, O) named
    through_left, through_right, input_reflection, output. Zero self-terms (no back-reflection at any port) and the
    standard symmetric-coupler cross structure: L<->R and I<->O direct
    ("through", amplitude t), L<->O and R<->I diagonal ("cross",
    amplitude k) -- the same structure directional_coupler_hybrid.py's own
    M_in0 ("in" excited -> through/cross) and M_in1 ("isolated" excited ->
    through/cross) columns populate, just written here as a dense 4x4."""
    M = len(np.atleast_1d(wavelengths))
    t = np.broadcast_to(np.atleast_1d(t), (M,))
    k = np.broadcast_to(np.atleast_1d(k), (M,))
    ports = ("through_left", "through_right", "input_reflection", "output")
    data = np.zeros((M, 4, 4), dtype=complex)
    L, R, I, O = 0, 1, 2, 3
    data[:, L, R] = t
    data[:, R, L] = t
    data[:, I, O] = t
    data[:, O, I] = t
    data[:, L, O] = k
    data[:, O, L] = k
    data[:, R, I] = k
    data[:, I, R] = k
    return SMatrixBlock(
        wavelengths=np.atleast_1d(wavelengths), ports=ports,
        modes_per_port={p: 1 for p in ports}, data=data,
    )


def _mirror_block(r, wavelengths=WAVELENGTHS):
    """Synthetic 1-port self-reflecting block -- the HCG mirror, in
    network terms (validation doc A3). r may be scalar or per-wavelength."""
    M = len(np.atleast_1d(wavelengths))
    r = np.broadcast_to(np.atleast_1d(r), (M,))
    data = np.zeros((M, 1, 1), dtype=complex)
    data[:, 0, 0] = r
    return SMatrixBlock(
        wavelengths=np.atleast_1d(wavelengths), ports=("p",),
        modes_per_port={"p": 1}, data=data,
    )


def _phase_block(neff, length, wavelengths=WAVELENGTHS, port_in="in", port_out="out"):
    """Single-mode waveguide, built as a genuinely BIDIRECTIONAL 2-port
    block -- NOT via beef.network.convert.from_transfer_matrix().

    Real finding (2026-08-05): from_transfer_matrix() is documented as
    one-directional on purpose (it wraps a HybridSolver-style transfer
    matrix, "self-terms — reflection at either port, or reverse
    transmission — are zero-filled, matching what those solvers actually
    model: pure forward transmission, no reflection") — correct for a
    straight-through cascade where nothing downstream ever reflects, wrong
    here: this cavity's connecting waveguides carry light BOTH forward
    (toward each HCG mirror) and backward (the mirror's own reflection,
    same physical waveguide). Using from_transfer_matrix() here silently
    connect()s the mirror to a waveguide that cannot carry anything back
    through it — the reduction comes out exactly zero, not a small error,
    since the reverse-direction entry is exactly zero-filled. Caught
    before any real device numbers were involved. A real, reciprocal,
    lossless waveguide must have S[out,in] == S[in,out] == phase (the
    forward and backward propagation constants of the SAME guided mode are
    identical — ordinary EM reciprocity) and no self-reflection at either
    end (true for any uniform straight section)."""
    wavelengths = np.atleast_1d(wavelengths)
    phase = np.exp(2j * np.pi * neff / wavelengths * length)
    M = len(wavelengths)
    data = np.zeros((M, 2, 2), dtype=complex)
    data[:, 1, 0] = phase   # port_out <- port_in
    data[:, 0, 1] = phase   # port_in <- port_out (the direction from_transfer_matrix() zero-fills)
    return SMatrixBlock(
        wavelengths=wavelengths, ports=(port_in, port_out),
        modes_per_port={port_in: 1, port_out: 1}, data=data,
    )


def _closed_form_symmetric(t, k, rho):
    D = 1 - rho ** 2 * t ** 2
    self_reflection = k ** 2 * rho / D
    transmission = t * (1 + rho ** 2 * (k ** 2 - t ** 2)) / D
    return self_reflection, transmission


def test_circuit_matches_closed_form_symmetric_cavity():
    """Full step-4-shaped Circuit assembly (5 named blocks, 4 connect()
    calls, matching the plan doc's own wiring sketch) against the
    independently hand-derived closed form -- the SAME cross-check role
    test_network_ring_crosscheck.py plays for the ring's own formula."""
    r = 0.97 * np.exp(1j * 0.31)          # HCG mirror, both ends (symmetric cavity)
    neff = 2.15                            # connecting-waveguide fundamental mode
    L = 20.0                               # um, both arms (step 1's pinned length)
    t = 0.85 * np.exp(1j * 0.12)           # DC through amplitude
    k = 0.45 * np.exp(1j * 1.05)           # DC cross amplitude

    circuit = Circuit()
    circuit.add_block("hcg_left", _mirror_block(r))
    circuit.add_block("wg_left", _phase_block(neff, L))
    circuit.add_block("dc", _dc_block(t, k))
    circuit.add_block("wg_right", _phase_block(neff, L))
    circuit.add_block("hcg_right", _mirror_block(r))

    circuit.connect("hcg_left", "p", "wg_left", "in")
    circuit.connect("wg_left", "out", "dc", "through_left")
    circuit.connect("dc", "through_right", "wg_right", "in")
    circuit.connect("wg_right", "out", "hcg_right", "p")

    cavity = circuit.reduce(external_ports={
        "dc.input_reflection": "input_reflection",
        "dc.output": "output",
    })

    assert set(cavity.ports) == {"input_reflection", "output"}

    p = np.exp(2j * np.pi * neff / WAVELENGTHS * L)   # single-pass phase
    rho = r * p ** 2                                    # mirror + ROUND-TRIP waveguide phase
    expected_self, expected_trans = _closed_form_symmetric(t, k, rho)

    np.testing.assert_allclose(
        cavity.submatrix("input_reflection", "input_reflection")[:, 0, 0],
        expected_self, atol=1e-10,
    )
    np.testing.assert_allclose(
        cavity.submatrix("output", "input_reflection")[:, 0, 0],
        expected_trans, atol=1e-10,
    )
    # Reciprocity: output->input_reflection must equal input_reflection->output.
    np.testing.assert_allclose(
        cavity.submatrix("input_reflection", "output")[:, 0, 0],
        expected_trans, atol=1e-10,
    )


def test_closed_form_matches_brute_force_4port_solve():
    """The closed form itself, independently, against a brute-force direct
    solve of the 4-port S-matrix plus its two boundary conditions -- so a
    bug shared between the closed form and the Circuit wiring (the failure
    mode a single cross-check cannot catch) does not silently pass."""
    rng = np.random.default_rng(0)
    for _ in range(200):
        t = (rng.random() * 0.3) * np.exp(1j * rng.random() * 2 * np.pi)
        k = (rng.random() * 0.3) * np.exp(1j * rng.random() * 2 * np.pi)
        rho = (rng.random() * 0.99) * np.exp(1j * rng.random() * 2 * np.pi)

        S = np.array([[0, t, 0, k], [t, 0, k, 0], [0, k, 0, t], [k, 0, t, 0]], dtype=complex)
        A = np.array([[1, -rho * t], [-rho * t, 1]], dtype=complex)
        # Excite input_reflection only (a_I=1, a_O=0): a_L = rho*k*a_O, a_R = rho*k*a_I
        rhs = np.array([rho * k * 0.0, rho * k * 1.0], dtype=complex)
        a_L, a_R = np.linalg.solve(A, rhs)
        a_vec = np.array([a_L, a_R, 1.0, 0.0], dtype=complex)
        b_vec = S @ a_vec
        brute_self, brute_trans = b_vec[2], b_vec[3]

        closed_self, closed_trans = _closed_form_symmetric(t, k, rho)
        np.testing.assert_allclose(brute_self, closed_self, atol=1e-9)
        np.testing.assert_allclose(brute_trans, closed_trans, atol=1e-9)


def test_passivity_holds_for_lossy_mirrors():
    """|r|<1 (a real, lossy mirror, like the HCG's own measured ~97%
    reflectance) must not make the assembled cavity emit more power than
    it receives -- a basic sanity floor before trusting anything from a
    real device's numbers.

    t/k here use the standard symmetric-coupler quadrature form
    (t=cos(theta), k=i*sin(theta)) rather than arbitrary phases: a genuine
    lossless 4-port S-matrix of this block's structure requires not just
    |t|^2+|k|^2=1 at each port but a strict 90-degree relative phase
    between t and k (Re(conj(t)*k)=0 -- derived directly from this
    S-matrix's own unitarity condition S^dagger S = I). Arbitrary phases
    fail that and can spuriously read as passivity violations that have
    nothing to do with the closed form or the mirror's own loss -- a real
    device's measured t/k will satisfy this automatically (energy
    conservation in a real lossless-material simulation), synthetic test
    values will not unless chosen to."""
    r = 0.97 * np.exp(1j * 0.4)
    neff, L = 2.15, 20.0
    theta = 0.6   # radians -- an arbitrary but physically valid coupling angle
    t, k = np.cos(theta) + 0j, 1j * np.sin(theta)

    p = np.exp(2j * np.pi * neff / WAVELENGTHS * L)
    rho = r * p ** 2
    self_r, trans = _closed_form_symmetric(t, k, rho)
    total = np.abs(self_r) ** 2 + np.abs(trans) ** 2
    assert np.all(total <= 1.0 + 1e-9), f"passivity violated: max total power {total.max()}"
