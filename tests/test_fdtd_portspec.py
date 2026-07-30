"""Fast, meep-free unit tests for beef.solvers.fdtd.PortSpec and
solve_section()'s y0/y1/excite_port/measure_port argument validation.

No meep import needed: PortSpec.__post_init__ is plain dataclass
validation, and solve_section()'s new argument checks run before `import
meep as mp` executes, so they're exercised here without ever touching a
real simulation (module-level beef.solvers.fdtd has no meep import either
— confirmed by construction: FDTDSolver() below never imports meep).
"""
from __future__ import annotations

import pytest

from beef.solvers.fdtd import FDTDSolver, PortSpec


class TestPortSpec:
    def test_defaults(self):
        # width has no valid default (a PortSpec always fully specifies its
        # own aperture, by design — see the class docstring), so a bare
        # PortSpec() is itself invalid; only axis/sign have real defaults.
        p = PortSpec(width=1.0)
        assert p.axis == "x"
        assert p.sign == 1
        assert p.alpha_idx == 0

    def test_alpha_idx_backward(self):
        p = PortSpec(sign=-1, width=1.0)
        assert p.alpha_idx == 1

    def test_invalid_axis_raises(self):
        with pytest.raises(ValueError, match="axis must be"):
            PortSpec(axis="z", width=1.0)

    def test_invalid_sign_raises(self):
        with pytest.raises(ValueError, match="sign must be"):
            PortSpec(sign=0, width=1.0)

    def test_nonpositive_width_raises(self):
        with pytest.raises(ValueError, match="width must be"):
            PortSpec(width=0.0)
        with pytest.raises(ValueError, match="width must be"):
            PortSpec(width=-1.0)


class TestSolveSectionPortValidation:
    """Argument-shape checks solve_section() runs before import meep."""

    def _solver(self) -> FDTDSolver:
        return FDTDSolver()

    def test_only_one_of_y0_y1_raises(self):
        solver = self._solver()
        with pytest.raises(ValueError, match="y0 and y1 must both be given"):
            solver.solve_section(None, [1.55], 0.0, 1.0, y0=0.0, y1=None)
        with pytest.raises(ValueError, match="y0 and y1 must both be given"):
            solver.solve_section(None, [1.55], 0.0, 1.0, y0=None, y1=1.0)

    def test_y_axis_excite_port_without_y0_y1_raises(self):
        solver = self._solver()
        with pytest.raises(ValueError, match="excite_port.axis='y' requires y0 and y1"):
            solver.solve_section(
                None, [1.55], 0.0, 1.0,
                excite_port=PortSpec(axis="y", width=1.0),
            )

    def test_y_axis_measure_port_without_y0_y1_raises(self):
        solver = self._solver()
        with pytest.raises(ValueError, match="measure_port.axis='y' requires y0 and y1"):
            solver.solve_section(
                None, [1.55], 0.0, 1.0,
                measure_port=PortSpec(axis="y", width=1.0),
            )

    def test_x_axis_ports_without_y0_y1_do_not_raise_the_y_check(self):
        """An x-normal PortSpec override needs no y0/y1 — this should get
        past the argument-validation stage cleanly (it will go on to need
        a real geometry/meep, which we don't provide here, so we only
        assert it does NOT raise the y0/y1-related ValueError)."""
        solver = self._solver()
        with pytest.raises(AttributeError):
            # None has no section_prisms()/width_at() etc. — proves
            # validation passed and execution reached real geometry use.
            solver.solve_section(
                None, [1.55], 0.0, 1.0,
                excite_port=PortSpec(axis="x", width=1.0),
                measure_port=PortSpec(axis="x", sign=-1, width=1.0),
            )
