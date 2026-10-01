"""Real-meow verification for EMESolver.solve_section()'s multi-cell body
support.

The two tests below raise before any geometry/mesh work or eigensolve --
meow is imported by solve_section() itself (unconditionally, at its top) but
never actually reached, so they run fast despite needing the real `beef`
package (and so the jmh env, same as every other test in this repo -- there
is no meep/meow-free CI split here, only a fast/slow *cost* split).

@pytest.mark.slow tests do a real multi-cell solve, matching this repo's
existing convention for real-meep/meow checks (run via
`pytest tests/test_eme_multicell_section.py --runslow`).
"""
from __future__ import annotations

import numpy as np
import pytest

from beef import EMESolver, SIN_300
from beef.geometry.mmi import MMI

_COMMON = dict(oxide_width=4.0, oxide_height=1.2, res=30.0,
               material_stack=SIN_300, mode_solver="tidy3d")


def _mmi(**kw):
    return MMI(n_inputs=1, n_outputs=2, mm_length=6.0, mm_width=2.5,
               waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
               waveguide_offset=1.0, **kw)


def test_solve_section_rejects_multi_region_solver():
    """A whole-device, multi-region EMESolver (e.g. one built for solve())
    reaching solve_section() is a caller error -- silently reading only its
    first region would produce a plausible, wrong answer with no warning,
    exactly the failure class step 1's assert_uniform_cross_section() guards
    against for neff_override."""
    solver = EMESolver(num_modes=[2, 4], num_cells=[1, 1], lengths=[1.0, 1.0],
                        **_COMMON)
    with pytest.raises(ValueError, match="single-region"):
        solver.solve_section(_mmi(), np.array([1.55]), -2.0, 2.0,
                              n_modes_in=2, n_modes_out=2)


def test_solve_section_does_not_read_num_modes():
    """The per-cell mode count comes from n_modes_in/n_modes_out, NOT from
    self.num_modes -- unchanged by step 3, and actively relied on by
    modal_multiplexer_hybrid.py (one num_modes=[15] body_solver fills three
    EME slots, one of which is a CONN_MODES=6 connector). Two solvers
    differing only in num_modes must return the identical transfer matrix.

    Fast despite being a real solve: a 2 um uniform span at 2 modes.
    """
    m = _mmi()
    wl = np.array([1.55])
    a = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0], **_COMMON)
    b = EMESolver(num_modes=[12], num_cells=[1], lengths=[1.0], **_COMMON)
    Ta = a.solve_section(m, wl, -1.0, 1.0, n_modes_in=2, n_modes_out=2)
    Tb = b.solve_section(m, wl, -1.0, 1.0, n_modes_in=2, n_modes_out=2)
    np.testing.assert_allclose(
        Ta, Tb, rtol=0, atol=0,
        err_msg="num_modes changed solve_section()'s answer -- it must not be "
                "read on this path (see modal_multiplexer_hybrid.py)",
    )


@pytest.mark.slow
def test_multicell_matches_single_cell_on_a_uniform_span():
    """num_cells=[1] vs num_cells=[3] on the SAME uniform-width span must
    agree closely. Propagating through several cells that all share one
    cross-section is mathematically identical to one cell of the summed
    length (same per-mode phase either way) -- any difference here is purely
    per-cell eigensolve noise across 3 independent mode solves of the
    identical structure, not physics. This is step 3's own "untapered limit
    reproduces the single-cell result" regression, and doubles as end-to-end
    coverage of the num_cells=[1] path being byte-identical to pre-step-3
    behaviour (same Ls=[length] it always built)."""
    m = _mmi()  # body_taper_length=0 (default) -- straight body
    solver_1 = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0], **_COMMON)
    solver_3 = EMESolver(num_modes=[2], num_cells=[3], lengths=[1.0], **_COMMON)
    x0, x1 = -2.0, 2.0  # inside the uniform body, clear of the port tapers
    wl = np.array([1.55])

    T1 = solver_1.solve_section(m, wl, x0, x1, n_modes_in=2, n_modes_out=2)
    T3 = solver_3.solve_section(m, wl, x0, x1, n_modes_in=2, n_modes_out=2)

    np.testing.assert_allclose(
        np.abs(T1), np.abs(T3), atol=5e-3,
        err_msg=f"|T| disagreement between 1-cell and 3-cell solves of the "
                f"same uniform span: T1={T1}, T3={T3}",
    )


@pytest.mark.slow
def test_multicell_body_taper_changes_the_transfer_matrix():
    """A real multi-cell solve over a tapered body must give a genuinely
    different transfer matrix than the same multi-cell machinery applied to
    a straight body -- otherwise the taper geometry (step 2) is not actually
    reaching the eigensolve, e.g. because every cell sampled the same
    (wrong) cross-section."""
    m_tapered = _mmi(body_taper_length=1.0, body_start_width=1.5, body_end_width=1.5)
    m_straight = _mmi()
    solver = EMESolver(num_modes=[2], num_cells=[6], lengths=[1.0], **_COMMON)
    x0, x1 = -3.0, 3.0  # full body span, taper bands included on both ends
    wl = np.array([1.55])

    T_tapered = solver.solve_section(m_tapered, wl, x0, x1, n_modes_in=2, n_modes_out=2)
    T_straight = solver.solve_section(m_straight, wl, x0, x1, n_modes_in=2, n_modes_out=2)

    # Measured (job 12069743): |T_straight| is 1.0 on the diagonal (a straight
    # uniform body is exactly lossless per mode, as expected), while
    # |T_tapered| comes back 0.9802/0.9884 -- a real, mode-dependent ~1-2%
    # amplitude effect from the taper (radiation/mode-conversion loss a
    # uniform body cannot have), an order of magnitude past the ~5e-3
    # cross-solve numerical noise floor test_multicell_matches_single_cell_
    # on_a_uniform_span establishes just above. atol here sits between the
    # two so this only fails if the taper's effect vanishes, not on solver
    # noise.
    assert not np.allclose(np.abs(T_tapered), np.abs(T_straight), atol=5e-3), (
        f"tapered-body multi-cell transfer matrix indistinguishable from the "
        f"straight-body one -- the taper is not affecting the solve: "
        f"T_tapered={T_tapered}, T_straight={T_straight}"
    )
