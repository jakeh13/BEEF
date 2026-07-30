"""close_loop() against the ring resonator's own proven closed-form formula.

Single-purpose: no HybridSolver/scripted-FakeSolver machinery, just literal
numbers directly reproducing the exact computation in
examples/ring_resonator/python/ring_resonator_hybrid.py:762-774 (identical
in ring_resonator_hybrid_cascaded.py:817-825):

    P = (b_corrected ** 2) * np.exp(1j * phi)
    E_thru = t_d + k1_d * k2_d * P / (1.0 - t_ring_d * P)

where t, k1, k2, t_ring are the coupler section's own 2-port S-matrix
entries in a [bus, ring] basis and P is the round-trip loop transfer
(bend loss + propagation phase). This is what close_loop() must reduce to
exactly when applied to the same 2-port block with loop=P.
"""
from __future__ import annotations

import warnings

import numpy as np

from beef.network.block import SMatrixBlock
from beef.network.ops import close_loop

WAVELENGTHS = np.linspace(1.53, 1.57, 5)


def _coupler_block(t, k1, k2, t_ring, wavelengths=WAVELENGTHS):
    M = len(wavelengths)
    data = np.zeros((M, 2, 2), dtype=complex)
    data[:, 0, 0] = t       # bus -> bus
    data[:, 0, 1] = k2      # ring -> bus
    data[:, 1, 0] = k1      # bus -> ring
    data[:, 1, 1] = t_ring  # ring -> ring
    return SMatrixBlock(
        wavelengths=wavelengths, ports=("bus", "ring"),
        modes_per_port={"bus": 1, "ring": 1}, data=data,
    )


def test_close_loop_reproduces_ring_all_pass_formula():
    t = 0.95 * np.exp(1j * 0.10)
    k1 = 0.30 * np.exp(1j * 1.20)
    k2 = 0.30 * np.exp(1j * 1.20)  # reciprocal coupler: k2 == k1
    t_ring = 0.95 * np.exp(1j * 0.05)
    P = 0.90 * np.exp(1j * 2.37) * np.ones(len(WAVELENGTHS))

    block = _coupler_block(t, k1, k2, t_ring)
    result = close_loop(block, "ring", "ring", loop=P)

    expected = t + k1 * k2 * P / (1.0 - t_ring * P)
    np.testing.assert_allclose(result.submatrix("bus", "bus")[:, 0, 0], expected, atol=1e-12)
    assert result.ports == ("bus",)


def test_close_loop_reproduces_ring_formula_across_wavelength_sweep():
    """Same check, but with t/k1/k2/t_ring/P varying per wavelength — the
    real shape of ring_resonator_hybrid.py's dense-grid assembly, not a
    single scalar snapshot."""
    wl = WAVELENGTHS
    t = 0.95 * np.exp(1j * np.linspace(0.05, 0.15, len(wl)))
    k1 = 0.30 * np.exp(1j * np.linspace(1.1, 1.3, len(wl)))
    k2 = k1.copy()
    t_ring = 0.95 * np.exp(1j * np.linspace(0.0, 0.1, len(wl)))
    P = 0.90 * np.exp(1j * np.linspace(2.0, 2.6, len(wl)))

    block = _coupler_block(t, k1, k2, t_ring, wavelengths=wl)
    result = close_loop(block, "ring", "ring", loop=P)

    expected = t + k1 * k2 * P / (1.0 - t_ring * P)
    np.testing.assert_allclose(result.submatrix("bus", "bus")[:, 0, 0], expected, atol=1e-12)


def test_conditioning_warning_fires_approaching_critical_coupling():
    """Sweep t_ring*P toward 1 (critical coupling) and confirm the
    conditioning warning fires exactly where the resonance narrows, not
    before — documents the connection between 'near-critical coupling'
    and 'near-singular loop closure' concretely."""
    t_ring = 0.95 + 0j
    fired_at = None
    for exponent in range(1, 12):
        P_mag = 1.0 - 10.0 ** (-exponent)
        P = np.full(len(WAVELENGTHS), (P_mag / abs(t_ring)) + 0j)
        block = _coupler_block(0.95, 0.3, 0.3, t_ring)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            close_loop(block, "ring", "ring", loop=P)
        if caught:
            fired_at = exponent
            break
    assert fired_at is not None, "conditioning warning never fired approaching critical coupling"
    # Sanity: it should fire only once t_ring*P is within ~1e-6 of 1, not earlier.
    assert fired_at >= 5
