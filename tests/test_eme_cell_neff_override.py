"""Per-cell n_eff calibration for a TAPERED EME body
(`EMESolver.solve_section(cell_neff_override=...)`).

Why this exists. `neff_override` calibrates a body onto the FDTD's own beta,
which is what closed the engine-agreement gap on every uniform-body device in
this repo. It is hard-refused for a tapered body because it works by replacing
the section with a DIAGONAL exp(2j*pi*neff*L/wl) -- no mode mixing -- and a
taper's whole physics is conversion between modes. `cell_neff_override` is the
tapered counterpart: it keeps the multi-cell solve and every interface
S-matrix, and substitutes n_eff only in the propagation within each cell.

The split it relies on is exact, and the identity test below is what pins it:
interface overlaps integrate the mode FIELDS, propagation uses n_eff, so
substituting n_eff cannot perturb an interface.

Cheap tests raise before any eigensolve. @pytest.mark.slow ones do a real
multi-cell meow solve, per this repo's convention
(`pytest tests/test_eme_cell_neff_override.py --runslow`).
"""
from __future__ import annotations

import numpy as np
import pytest

from beef import EMESolver, SIN_300
from beef.geometry.mmi import MMI

_COMMON = dict(oxide_width=4.0, oxide_height=1.2, res=30.0,
               material_stack=SIN_300, mode_solver="tidy3d")
_WL = np.array([1.55])


def _mmi(**kw):
    return MMI(n_inputs=1, n_outputs=2, mm_length=6.0, mm_width=2.5,
               waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
               waveguide_offset=1.0, **kw)


def _solver(n_cells):
    return EMESolver(num_modes=[4], num_cells=[n_cells], lengths=[2.0],
                     **_COMMON)


def test_rejects_combining_the_two_calibration_models():
    """They are different propagation models, not two settings of one."""
    solver = _solver(3)
    with pytest.raises(ValueError, match="different"):
        solver.solve_section(
            _mmi(), _WL, -1.0, 1.0, n_modes_in=2, n_modes_out=2,
            neff_override=np.ones((1, 2)),
            cell_neff_override=np.ones((1, 3, 2)),
        )


def test_rejects_wrong_cell_count():
    """A shape mismatch would calibrate the WRONG cell -- silently, since any
    [wl, N, modes] array is structurally plausible."""
    solver = _solver(4)
    with pytest.raises(ValueError, match=r"expected|shape"):
        solver.solve_section(
            _mmi(), _WL, -1.0, 1.0, n_modes_in=2, n_modes_out=2,
            cell_neff_override=np.ones((1, 3, 2)),   # 3 cells, section has 4
        )


def test_rejects_too_few_modes():
    solver = _solver(3)
    with pytest.raises(ValueError, match="modes"):
        solver.solve_section(
            _mmi(), _WL, -1.0, 1.0, n_modes_in=3, n_modes_out=3,
            cell_neff_override=np.ones((1, 3, 2)),   # 2 modes, section has 3
        )


def test_fingerprint_unchanged_when_override_absent():
    """THE BACKWARD-COMPATIBILITY TEST.

    Every device in this repo shares one on-disk cache, and the expensive
    entries are FDTD sections. If adding this parameter changed the
    solve_section fingerprint for runs that DO NOT use it, every cached entry
    would be invalidated at once and every in-flight study would silently
    re-solve from scratch. The config for an uncalibrated call must therefore
    be byte-identical to what it was before the parameter existed -- i.e. it
    must not contain the key at all.
    """
    from beef.cache.config import build_run_config

    solver = _solver(3)
    geom = _mmi()
    base = dict(x0=-1.0, x1=1.0, wavelengths=_WL, stub_left=0.0,
                stub_right=0.0, geometry_kwargs={}, neff_override=None)
    cfg = build_run_config(solver, "solve_section", base, geom)
    assert "cell_neff_override" not in cfg["call_kwargs"], (
        "an uncalibrated run must not carry the new key in its fingerprint; "
        "including it unconditionally invalidates the entire shared cache")


def test_fingerprint_changes_when_override_present():
    """The flip side: a calibrated run must NOT collide with an uncalibrated
    one, or the two would silently serve each other from cache."""
    from beef.cache.config import build_run_config
    from beef.solvers.eme import _cell_neff_fingerprint

    solver = _solver(3)
    geom = _mmi()
    base = dict(x0=-1.0, x1=1.0, wavelengths=_WL, stub_left=0.0,
                stub_right=0.0, geometry_kwargs={}, neff_override=None)
    plain = build_run_config(solver, "solve_section", base, geom)
    calibrated = build_run_config(
        solver, "solve_section",
        dict(base, cell_neff_override=_cell_neff_fingerprint(
            np.full((1, 3, 2), 2.1))),
        geom)
    assert plain != calibrated


def test_fingerprint_distinguishes_different_betas():
    """Two different calibrations must not share an entry either."""
    from beef.cache.config import build_run_config
    from beef.solvers.eme import _cell_neff_fingerprint

    solver = _solver(3)
    geom = _mmi()
    base = dict(x0=-1.0, x1=1.0, wavelengths=_WL, stub_left=0.0,
                stub_right=0.0, geometry_kwargs={}, neff_override=None)
    a = build_run_config(solver, "solve_section", dict(
        base, cell_neff_override=_cell_neff_fingerprint(np.full((1, 3, 2), 2.1))), geom)
    b = build_run_config(solver, "solve_section", dict(
        base, cell_neff_override=_cell_neff_fingerprint(np.full((1, 3, 2), 2.2))), geom)
    assert a != b


def test_cell_fingerprint_rejects_wrong_rank():
    from beef.solvers.eme import _cell_neff_fingerprint
    assert _cell_neff_fingerprint(None) is None
    with pytest.raises(ValueError, match="n_cells"):
        _cell_neff_fingerprint(np.ones((2, 3)))


@pytest.mark.slow
def test_identity_substituting_the_modes_own_neff_reproduces_the_plain_solve():
    """THE CORRECTNESS TEST.

    Feeding back each cell's OWN solved n_eff must reproduce the uncalibrated
    transfer matrix to machine precision. That is the sharpest available check
    that the substitution touches propagation and NOTHING else: if it
    perturbed an interface overlap, disturbed mode ordering, or rebuilt the
    Mode objects wrongly, this would not close to 1e-12.
    """
    import meow as mw

    solver = _solver(3)
    geom = _mmi(body_waist_width=2.0, body_waist_length=0.6,
                body_transition_length=0.5)
    kw = dict(n_modes_in=3, n_modes_out=3)
    plain = solver.solve_section(geom, _WL, -1.0, 1.0, **kw)

    # Re-derive the same cells/modes the solve just used, and read their n_eff.
    from beef._utils import extrude_polygons
    gds = geom.section_gds(-1.0, 1.0, solver.oxide_width)
    structs = extrude_polygons(gds, solver._extrusion_rules(_WL))
    cells = mw.create_cells(structures=structs, mesh=solver._build_mesh(),
                            Ls=solver._flat_cell_lengths(region_lengths=[2.0]))
    env = mw.Environment(wl=float(_WL[0]), T=solver.temperature)
    css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
    modes = [solver.compute_modes(cs, 3) for cs in css]
    own = np.array([[[complex(m.neff) for m in cell_modes] for cell_modes in modes]])

    calibrated = solver.solve_section(geom, _WL, -1.0, 1.0,
                                      cell_neff_override=own, **kw)
    assert np.allclose(plain, calibrated, rtol=0, atol=1e-12), (
        f"max |dT| = {np.abs(plain - calibrated).max():.3e}; substituting a "
        "cell's own n_eff must be a no-op")


@pytest.mark.slow
def test_changing_cell_neff_actually_changes_the_result():
    """Guards the opposite failure: a parameter that is silently ignored.

    The seam-modes episode in mmi_2x2_waist_body_validation.md is exactly this
    -- a knob that looked wired up, produced bit-identical output, and cost a
    job to discover. Perturbing n_eff must move the answer.
    """
    solver = _solver(3)
    geom = _mmi(body_waist_width=2.0, body_waist_length=0.6,
                body_transition_length=0.5)
    kw = dict(n_modes_in=3, n_modes_out=3)
    plain = solver.solve_section(geom, _WL, -1.0, 1.0, **kw)
    bumped = solver.solve_section(
        geom, _WL, -1.0, 1.0,
        cell_neff_override=np.full((1, 3, 3), 1.9), **kw)
    assert not np.allclose(plain, bumped, atol=1e-9)


def test_hybrid_cell_neff_fn_defaults_to_none():
    """The hook must be invisible unless asked for.

    Every hybrid device in this repo (ring, 1x2, 2x2, 1x4, DC, HCG-FP, modal
    mux, 4x4) constructs HybridSolver without it, and several have jobs in
    flight that import beef/ at start time.
    """
    from beef import HybridSolver, EMESolver, FDTDSolver, SIN_300
    import inspect

    sig = inspect.signature(HybridSolver.__init__)
    assert sig.parameters["cell_neff_fn"].default is None

    common = dict(oxide_width=4.0, oxide_height=1.2, res=30.0,
                  material_stack=SIN_300, mode_solver="tidy3d")
    eme = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0], **common)
    fdtd = FDTDSolver(res=20.0, sim_width=4.0, sim_height=1.2,
                      material_stack=SIN_300)
    hyb = HybridSolver(solvers=[fdtd, eme, fdtd], boundaries=[-2., -1., 1., 2.],
                       interface_modes=[2, 2])
    assert hyb.cell_neff_fn is None


def test_capability_probe_distinguishes_the_two_solvers():
    """EME grew cell_neff_override; FDTD did not. The probe is what keeps the
    hook from being handed to a solver that would reject it as an unexpected
    kwarg."""
    from beef import EMESolver, FDTDSolver, SIN_300
    from beef.solvers.hybrid import _supports_cell_neff_override

    eme = EMESolver(num_modes=[2], num_cells=[1], lengths=[1.0],
                    oxide_width=4.0, oxide_height=1.2, res=30.0,
                    material_stack=SIN_300, mode_solver="tidy3d")
    fdtd = FDTDSolver(res=20.0, sim_width=4.0, sim_height=1.2,
                      material_stack=SIN_300)
    assert _supports_cell_neff_override(eme.solve_section)
    assert not _supports_cell_neff_override(fdtd.solve_section)
