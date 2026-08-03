"""interp_complex: the ring's coefficient resampling, measured.

The load-bearing test here is test_beats_real_imag_on_a_phase_ramp. The
strategy choice (magnitude + unwrapped phase, not real + imag) is currently
justified by a paragraph in beef/spectral.py and by one observed episode on the
ring, where a Q table built the wrong way showed Q never converging. A
paragraph is not a test. This measures the difference on a function whose exact
value is known everywhere, so the claim either holds numerically or it does
not.

The second is test_unwrap_breaks_when_the_source_grid_is_too_coarse, which
pins the failure mode of the CHOSEN strategy. A test suite that only shows the
preferred option winning is advocacy; this one shows where it loses.
"""
import numpy as np
import pytest

from beef.spectral import interp_complex, unwrap_is_safe


def _ramp(wl, beta=120.0, amp=0.9):
    """A·exp(i·β·λ) — a pure propagation phase, which is what the ring's
    coefficients mostly are: |t|, |k| vary slowly while the phase winds fast.
    β=120 rad/µm over a 0.1 µm band is ~12 rad, a couple of full turns, which
    is the ring's actual regime rather than a contrived extreme."""
    return amp * np.exp(1j * beta * wl)


def _real_imag_interp(wl_src, c, wl_dense):
    """The rejected alternative, spelled out so the comparison is explicit."""
    return (np.interp(wl_dense, wl_src, np.real(c))
            + 1j * np.interp(wl_dense, wl_src, np.imag(c)))


# ---------------------------------------------------------------------------
# The one that decides the strategy
# ---------------------------------------------------------------------------
def test_beats_real_imag_on_a_phase_ramp():
    """On a pure phase ramp the mag/phase interpolant is EXACT (up to float
    noise) while real/imag chords across every oscillation and loses
    magnitude. Both are compared against the analytic value, not against each
    other, so neither is privileged by the setup."""
    wl_src = np.linspace(1.50, 1.60, 101)          # the ring's N_WL_FDTD
    wl_dense = np.linspace(1.50, 1.60, 20001)      # the ring's N_WL_DENSE
    exact = _ramp(wl_dense)

    got_phase = interp_complex(wl_src, _ramp(wl_src), wl_dense)
    got_reim = _real_imag_interp(wl_src, _ramp(wl_src), wl_dense)

    err_phase = float(np.max(np.abs(got_phase - exact)))
    err_reim = float(np.max(np.abs(got_reim - exact)))

    # The phase route reproduces a pure ramp to machine precision: linear
    # interpolation of a linear function is that function.
    assert err_phase < 1e-12, f"mag/phase error {err_phase:.3e} should be ~0"
    # And real/imag is worse by orders of magnitude, not by a hair.
    assert err_reim > 1e4 * err_phase, (
        f"real/imag error {err_reim:.3e} vs mag/phase {err_phase:.3e} — the "
        f"strategy claim in beef/spectral.py does not hold at this sampling")

    # The specific damage is a MAGNITUDE deficit: chording a circle cuts
    # inside it, so |interpolant| <= |true| everywhere and never above.
    mag_reim = np.abs(got_reim)
    assert mag_reim.max() <= np.abs(exact).max() + 1e-12
    assert mag_reim.min() < 0.999 * np.abs(exact).min(), (
        "real/imag interpolation should visibly under-shoot the magnitude; if "
        "it does not, this sampling is too fine for the test to be meaningful")


def test_the_magnitude_deficit_is_what_corrupts_a_linewidth():
    """Why the above matters for Q specifically: the deficit is largest
    mid-interval and zero at the source points, so it imprints a ripple at the
    SOURCE grid period onto |c| — which is a spurious spectral feature at a
    spacing unrelated to the device."""
    wl_src = np.linspace(1.50, 1.60, 51)
    wl_dense = np.linspace(1.50, 1.60, 5001)
    got_reim = _real_imag_interp(wl_src, _ramp(wl_src), wl_dense)
    # |c| should be flat 0.9; measure how far the real/imag route departs.
    ripple = float(np.max(np.abs(np.abs(got_reim) - 0.9)))
    assert ripple > 1e-3, (
        f"expected a visible magnitude ripple from real/imag interpolation, "
        f"got {ripple:.3e}")
    # The chosen route has none.
    got_phase = interp_complex(wl_src, _ramp(wl_src), wl_dense)
    assert float(np.max(np.abs(np.abs(got_phase) - 0.9))) < 1e-12


# ---------------------------------------------------------------------------
# Where the chosen strategy loses
# ---------------------------------------------------------------------------
def test_unwrap_breaks_when_the_source_grid_is_too_coarse():
    """Past pi per source sample, np.unwrap picks the wrong branch and the
    interpolant is confidently wrong — a worse failure than real/imag's
    blurring, since nothing about the output looks damaged."""
    beta = 400.0                                    # rad/um: fast
    wl_src = np.linspace(1.50, 1.60, 9)             # ~5 rad per step: aliased
    wl_dense = np.linspace(1.50, 1.60, 2001)
    c_src = _ramp(wl_src, beta=beta)

    got = interp_complex(wl_src, c_src, wl_dense)
    err = float(np.max(np.abs(got - _ramp(wl_dense, beta=beta))))
    assert err > 0.5, (
        f"expected the aliased unwrap to produce a large error, got {err:.3e}")
    # And the damage is invisible in the magnitude — |c| is still a flat 0.9,
    # so no sanity check on amplitude would catch it.
    assert np.allclose(np.abs(got), 0.9, atol=1e-12)


def test_reports_safe_on_already_aliased_data():
    """The guard's FALSE NEGATIVE, pinned deliberately.

    np.unwrap maps every step into (-pi, pi], so the post-unwrap step this
    function measures can never exceed pi however badly undersampled the input
    is. The beta=400 ramp above advances 5.0 rad per sample; that folds to an
    observed 1.28 rad (0.41pi) and unwrap_is_safe() calls it safe while the
    interpolant is off by more than 0.5 in absolute amplitude.

    This is a fundamental limit, not a bug to fix later: recovering a phase
    rate above Nyquist from one sampled sequence is undecidable, and the
    aliased ramp stays self-consistent under decimation, so no data-only
    refinement test rescues it. Recorded as a test so the guard's advertised
    scope cannot quietly widen.
    """
    wl_src = np.linspace(1.50, 1.60, 9)
    c_src = _ramp(wl_src, beta=400.0)
    safe, worst = unwrap_is_safe(c_src)
    assert safe, "if this ever starts failing, the limitation has been fixed"
    assert worst < np.pi / 2
    assert worst == pytest.approx(abs(400.0 * (0.1 / 8) - 2 * np.pi), abs=1e-9)
    # check_unwrap therefore does NOT refuse it. That is the point.
    interp_complex(wl_src, c_src, np.linspace(1.50, 1.60, 101), check_unwrap=0.5)


def test_guard_fires_on_a_marginal_but_not_yet_aliased_grid():
    """What the guard IS for: a true step between pi/2 and pi. Unwrap still
    picks the right branch here — the interpolation is exact — but a modest
    increase in phase rate would break it, which is precisely when a caller
    wants to be told."""
    beta = 200.0                                    # 2.5 rad/step = 0.80pi
    wl_src = np.linspace(1.50, 1.60, 9)
    wl_dense = np.linspace(1.50, 1.60, 2001)
    c_src = _ramp(wl_src, beta=beta)

    safe, worst = unwrap_is_safe(c_src)
    assert not safe
    assert np.pi / 2 < worst < np.pi
    # Not yet broken: a pure ramp below Nyquist still interpolates exactly.
    got = interp_complex(wl_src, c_src, wl_dense)
    assert np.max(np.abs(got - _ramp(wl_dense, beta=beta))) < 1e-12

    with pytest.raises(ValueError, match="close enough to the boundary"):
        interp_complex(wl_src, c_src, wl_dense, check_unwrap=0.5)


def test_unwrap_is_safe_on_a_grid_that_actually_resolves_the_phase():
    wl_src = np.linspace(1.50, 1.60, 101)
    safe, worst = unwrap_is_safe(_ramp(wl_src))
    assert safe and worst < np.pi / 2
    # and check_unwrap lets it through
    interp_complex(wl_src, _ramp(wl_src), np.linspace(1.50, 1.60, 501),
                   check_unwrap=0.5)


def test_unwrap_is_safe_handles_degenerate_input():
    assert unwrap_is_safe(np.array([1.0 + 0j])) == (True, 0.0)
    assert unwrap_is_safe(np.array([], dtype=complex)) == (True, 0.0)


# ---------------------------------------------------------------------------
# Contract
# ---------------------------------------------------------------------------
def test_reproduces_the_source_points_exactly():
    """Interpolating onto the source grid must be the identity — the minimum
    any resampler owes its inputs."""
    wl_src = np.linspace(1.50, 1.60, 21)
    c = _ramp(wl_src) * (1.0 + 0.3 * np.sin(60 * wl_src))   # non-flat magnitude
    got = interp_complex(wl_src, c, wl_src)
    assert np.allclose(got, c, atol=1e-12)


def test_does_not_extrapolate_beyond_the_source_span():
    """np.interp clamps at the edges. Documented behaviour, pinned here because
    an unwrapped phase ramp extrapolated even slightly compounds linearly, so a
    future switch to extrapolation must be a deliberate, visible change."""
    wl_src = np.linspace(1.50, 1.60, 11)
    c = _ramp(wl_src)
    got = interp_complex(wl_src, c, np.array([1.40, 1.70]))
    assert got[0] == pytest.approx(c[0], abs=1e-12)
    assert got[1] == pytest.approx(c[-1], abs=1e-12)


def test_rejects_mismatched_shapes_short_and_unsorted_grids():
    wl = np.linspace(1.5, 1.6, 5)
    with pytest.raises(ValueError, match="must match"):
        interp_complex(wl, _ramp(wl)[:4], wl)
    with pytest.raises(ValueError, match="at least 2"):
        interp_complex(wl[:1], _ramp(wl)[:1], wl)
    with pytest.raises(ValueError, match="strictly increasing"):
        interp_complex(wl[::-1], _ramp(wl), wl)


def test_matches_the_ring_scripts_inline_form():
    """Both ring examples defined this inline. Promoting it to the library only
    helps if it is the SAME function — otherwise the scripts change behaviour
    on import."""
    wl_src = np.linspace(1.50, 1.60, 101)
    wl_dense = np.linspace(1.50, 1.60, 4001)
    c = _ramp(wl_src) * (0.8 + 0.2 * np.cos(30 * wl_src))

    def inline(wl_src, c, wl_dense):        # verbatim from the ring scripts
        mag = np.interp(wl_dense, wl_src, np.abs(c))
        phase = np.interp(wl_dense, wl_src, np.unwrap(np.angle(c)))
        return mag * np.exp(1j * phase)

    assert np.array_equal(interp_complex(wl_src, c, wl_dense),
                          inline(wl_src, c, wl_dense))
