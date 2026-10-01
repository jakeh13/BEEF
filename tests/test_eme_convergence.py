"""EME mesh convergence for the even/odd splitting (dc_validation.md item 11).

WHY THIS EXISTS
---------------
Every example inherited EMESolver's default res=40, and at 40 the coupled-pair
even/odd splitting is 7% too large:

    EME res     20       40       60       80      100
    dbeta    0.1125   0.0957   0.0926   0.0894   0.0895     (rad/um, 1.55 um)

Measured in dc_validation.md Done item 6. It went unnoticed for two reasons
worth stating, because both are traps that will recur:

  1. The hybrid resolution sweep meshes BOTH engines at the same R, so it
     measures the two convergences convolved. Its curve looks flat because
     FDTD dominates the observable it plots (T_cross) while EME dominates
     dbeta. A flat joint curve is not evidence that either engine converged.
  2. dbeta is a *difference* of two n_eff around 2.3. A 0.1% error in each
     mode is a ~10% error in their difference, so the quantity the coupler
     actually depends on is roughly a thousand times more sensitive to mesh
     than the quantity anyone looks at.

meep's FDTD is converged at 40 on the same nominal grid -- its implied dbeta
moves 2.74% across res 30-90, non-monotonically. Subpixel smoothing is the
likely reason. So the fix is not "raise everything", it is "stop sharing one
knob between two solvers that converge at different rates", which is what the
directional coupler scripts now do (FDTD_RESOLUTION / EME_RESOLUTION).

These tests need meow + MPB, so --runslow only. The fast test below needs
neither and guards the thing most likely to silently regress: that the example
scripts still ask for the converged resolution at all.
"""
from __future__ import annotations

import pathlib
import re

import numpy as np
import pytest

# The value the directional coupler scripts pin, and the one the table above
# says is the first converged point.
EME_RES_CONVERGED = 80.0
EME_RES_DEFAULT = 40.0
WL = 1.55

# Anything at or below this counts as converged for dbeta. res 80 -> 100 moves
# 0.07%, so 1% is loose enough not to be flaky on a different MPB build and
# tight enough to catch the 7% error that prompted this.
CONVERGENCE_TOL = 0.01

DEVICE = dict(
    waveguide_width=0.45, gap=0.25, coupling_length=30.0,
    bend_length=9.0, waveguide_offset=3.0, bend_type="euler",
    waveguide_height=0.22,
)

_EXAMPLE_DIR = pathlib.Path(__file__).parent.parent / "examples" / "directional_coupler"
_SCRIPTS = _EXAMPLE_DIR / "python"
# directional_coupler_hybrid_stub_sweep.py and _modes_sweep.py were guarded
# here too. They moved to examples/directional_coupler/archive/ on 2026-07-31
# (retired pending new sweeping functionality) and then out
# of this repo entirely on 2026-09-18, to
# github.com/jakeh13/BEEF_archive. They are frozen there, so this bug cannot be
# newly introduced into them, and reaching across to a sibling checkout would
# make this suite fail on a clean clone of BEEF alone. The guard therefore
# covers the live example only.
_HYBRID_SCRIPTS = (
    "directional_coupler_hybrid.py",
)


def _script_source(name: str) -> str:
    path = _SCRIPTS / name
    if not path.exists():
        raise FileNotFoundError(f"{name} not found in {_SCRIPTS}")
    return path.read_text()


class TestExampleScriptsPinTheConvergedResolution:
    """Fast, no simulation. The physics tests below cost minutes and only run
    under --runslow, so on a normal test run these are the only thing standing
    between a careless edit and a silent 7% regression."""

    @pytest.mark.parametrize("script", _HYBRID_SCRIPTS)
    def test_declares_separate_engine_resolutions(self, script):
        src = _script_source(script)
        assert "FDTD_RESOLUTION" in src, f"{script}: no FDTD_RESOLUTION"
        assert "EME_RESOLUTION" in src, f"{script}: no EME_RESOLUTION"

    @pytest.mark.parametrize("script", _HYBRID_SCRIPTS)
    def test_eme_resolution_is_the_converged_one(self, script):
        src = _script_source(script)
        m = re.search(r"^EME_RESOLUTION\s*=\s*([0-9.]+)", src, re.M)
        assert m, f"{script}: EME_RESOLUTION not assigned at module level"
        assert float(m.group(1)) >= EME_RES_CONVERGED, (
            f"{script} sets EME_RESOLUTION={m.group(1)}, below the converged "
            f"{EME_RES_CONVERGED:.0f}; dbeta_EME would be overstated"
        )

    @pytest.mark.parametrize("script", _HYBRID_SCRIPTS)
    def test_no_shared_resolution_knob_survives(self, script):
        """One knob for two engines is the actual defect -- re-introducing a
        bare `res=RESOLUTION` would re-couple them even if both constants
        still exist."""
        src = _script_source(script)
        assert not re.search(r"res\s*=\s*RESOLUTION\b", src), (
            f"{script}: a shared RESOLUTION is being passed to a solver again"
        )


def _dbeta_eme(res):
    """Even/odd splitting (rad/um) for the uniform coupled pair at mesh `res`."""
    import meow as mw
    from beef import DirectionalCoupler, EMESolver, SOI_220, extrude_polygons

    dc = DirectionalCoupler(**DEVICE)
    eme = EMESolver(
        num_modes=[2], num_cells=[1], lengths=[1.0],
        oxide_width=6.0, oxide_height=1.5, res=float(res),
        material_stack=SOI_220, mode_solver="meep",
    )
    L = 2.0
    gds = dc.section_gds(-L / 2.0, +L / 2.0, eme.oxide_width)
    cells = mw.create_cells(
        structures=extrude_polygons(gds, eme.material_stack.extrusion_rules),
        mesh=eme._build_mesh(), Ls=[L],
    )
    cs = mw.CrossSection.from_cell(
        cell=cells[0], env=mw.Environment(wl=WL, T=eme.temperature))
    neff = [float(np.real(m.neff)) for m in eme.compute_modes(cs, 2)]
    assert neff[0] > neff[1], f"supermodes not in descending n_eff order: {neff}"
    return 2 * np.pi * (neff[0] - neff[1]) / WL


@pytest.mark.slow
def test_chosen_eme_resolution_is_converged():
    """The pinned res=80 must not move when the mesh is refined further.

    This is the claim the example scripts rest on, so it is asserted rather
    than commented. Compares 80 against 100 -- refining by 25% in each
    dimension, i.e. ~1.6x the 2-D cells.
    """
    coarse = _dbeta_eme(EME_RES_CONVERGED)
    fine = _dbeta_eme(100.0)
    rel = abs(fine - coarse) / coarse
    assert rel < CONVERGENCE_TOL, (
        f"dbeta_EME moved {rel * 100:.2f}% from res {EME_RES_CONVERGED:.0f} "
        f"({coarse:.6f}) to res 100 ({fine:.6f}); res "
        f"{EME_RES_CONVERGED:.0f} is not converged after all"
    )


@pytest.mark.slow
def test_the_old_default_was_not_converged():
    """Documents the defect itself, so the reason for the split survives.

    If this ever fails, EMESolver's default became adequate -- through a meow
    or MPB change, most likely -- and the separate EME_RESOLUTION knob in the
    example scripts can be reconsidered. It is a finding either way, which is
    why it is asserted rather than left as prose.
    """
    default = _dbeta_eme(EME_RES_DEFAULT)
    converged = _dbeta_eme(EME_RES_CONVERGED)
    rel = (default - converged) / converged
    assert rel > 0.03, (
        f"res {EME_RES_DEFAULT:.0f} now within {rel * 100:.2f}% of converged; "
        f"the item 11 defect appears to be gone -- re-check the split knobs"
    )
