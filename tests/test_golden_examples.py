"""Slow golden-value integration tests (notes.md rec #1b).

Each test builds the same device geometry (same physical dimensions) as the
corresponding example in examples/, but at reduced simulation resolution/
wavelength count/mode count, runs a real meep/meow simulation, and checks the
resulting transmission per port against a recorded reference within a loose
tolerance. These are regression guards against real correctness bugs (sign
flips, wrong port routing, aperture/normalization mistakes — see notes.md's
"Current library reference" for examples found and fixed by earlier
sessions), not physical-accuracy validation, so the tolerance is generous.

The _run_*() functions below (one per target) are the actual build+simulate
logic; each is reused verbatim by scripts/generate_golden_values.py to
(re)generate the tests/golden/*.json reference files these tests compare
against, so there's exactly one place that defines "how to run this target."

Skipped by default — run with `pytest tests/ --runslow`.
"""
from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest

from feces import (
    MMI, DirectionalCoupler, EMESolver, FDTDSolver, HybridSolver, SMatrix, SOI_160, SOI_220,
)

pytestmark = pytest.mark.slow

GOLDEN_DIR = pathlib.Path(__file__).parent / "golden"

WAVELENGTHS_MMI = np.linspace(1.26, 1.36, 3)
WAVELENGTHS_DC = np.linspace(1.50, 1.60, 3)


def _mmi_1x2() -> MMI:
    return MMI(
        n_inputs=1, n_outputs=2,
        mm_length=4.5, mm_width=2.0,
        waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
        waveguide_offset=1.0, waveguide_height=0.16,
    )


def run_mmi_1x2_eme() -> SMatrix:
    mmi = _mmi_1x2()
    solver = EMESolver(
        num_modes=[1, 3, 6, 3, 1], num_cells=[1, 4, 1, 4, 1],
        lengths=[0.5, 1.0, 4.5, 1.0, 0.5],
        oxide_width=5.0, oxide_height=1.5, res=20.0,
        material_stack=SOI_160, mode_solver="meep",
    )
    smat_top = solver.solve(mmi, WAVELENGTHS_MMI, ports=["out_top"], output_idx=1)
    smat_bot = solver.solve(mmi, WAVELENGTHS_MMI, ports=["out_bot"], output_idx=0)
    smat = smat_top
    smat._data.update(smat_bot._data)
    return smat


def run_mmi_1x2_fdtd() -> SMatrix:
    mmi = _mmi_1x2()
    solver = FDTDSolver(
        res=20.0, sim_width=5.0, sim_height=1.5, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        num_input_modes=1, num_output_modes=1,
        output_monitor_centers=[1.25, -1.25], output_monitor_widths=[2.5, 2.5],
        material_stack=SOI_160, decay_tol=1e-9,
    )
    return solver.solve(mmi, WAVELENGTHS_MMI, port="top", input_mode=1)


def run_mmi_1x2_hybrid() -> SMatrix:
    mmi = _mmi_1x2()
    n_modes = 6
    fdtd_common = dict(
        # dpml doubled vs the real example's 0.5 (and vs the pure-FDTD
        # target's solve()) — this hybrid section is a short interior
        # solve_section() span where PML is a large fraction of the domain;
        # at res=20 (half of production's res=40) 0.5 µm of PML is only 10
        # pixels, half production's absorption, and this domain diverged
        # (meep: fields are NaN or Inf) after ~20 min. 1.0 µm restores the
        # same PML pixel count production gets at res=40/dpml=0.5.
        res=20.0, sim_width=5.0, sim_height=1.5, dpml=1.0,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_160, decay_tol=1e-9,
    )
    input_solver = FDTDSolver(**fdtd_common)
    body_solver = EMESolver(
        num_modes=[n_modes], num_cells=[1], lengths=[1.0],
        oxide_width=5.0, oxide_height=1.5, res=20.0,
        material_stack=SOI_160, mode_solver="meep",
    )
    output_solver = FDTDSolver(half_domain_aperture=True, **fdtd_common)

    half_body = mmi.mm_length / 2
    half_total = mmi.total_length / 2
    body_stub = 0.25
    boundaries = [-half_total, -half_body + body_stub, half_body - body_stub, half_total]
    hybrid = HybridSolver(
        solvers=[input_solver, body_solver, output_solver],
        boundaries=boundaries, interface_modes=[n_modes, n_modes],
    )

    cache: dict = {}
    smat_top = hybrid.solve(mmi, WAVELENGTHS_MMI, port="out_top", cache=cache,
                             geometry_kwargs=[{}, {}, {"output_idx": 1}])
    smat_bot = hybrid.solve(mmi, WAVELENGTHS_MMI, port="out_bot", cache=cache,
                             geometry_kwargs=[{}, {}, {"output_idx": 0}])
    smat = smat_top
    smat._data.update(smat_bot._data)
    return smat


def run_directional_coupler_hybrid() -> SMatrix:
    dc = DirectionalCoupler(
        waveguide_width=0.45, gap=0.25, coupling_length=30.0, bend_length=9.0,
        waveguide_offset=3.0, bend_type="euler", waveguide_height=0.22,
    )
    n_modes = 3
    fdtd_solver = FDTDSolver(
        # dpml doubled vs the real example's 0.5 — same reasoning as
        # run_mmi_1x2_hybrid's fdtd_common above (short interior
        # solve_section() span, PML pixel count halved at res=20 vs
        # production's res=40 unless dpml compensates).
        res=20.0, sim_width=6.0, sim_height=1.5, dpml=1.0,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        half_domain_aperture=True, material_stack=SOI_220, decay_tol=1e-9,
    )
    body_solver = EMESolver(
        num_modes=[n_modes], num_cells=[1], lengths=[1.0],
        oxide_width=6.0, oxide_height=1.5, res=20.0,
        material_stack=SOI_220, mode_solver="meep",
    )

    half_body = dc.coupling_length / 2
    half_total = dc.total_length / 2
    boundaries = [-half_total, -half_body, half_body, half_total]
    hybrid = HybridSolver(
        solvers=[fdtd_solver, body_solver, fdtd_solver],
        boundaries=boundaries, interface_modes=[n_modes, n_modes],
    )
    # Re-center the S-bend/coupling-region interior interface monitor on the
    # shared y=0 axis (matches the EME body's own supermode basis) — same
    # override examples/directional_coupler/directional_coupler_hybrid.py
    # applies manually via _iface_mon; needed on both the input and output
    # sections here since (unlike that example's reciprocity shortcut) each
    # section below is actually simulated, not reused via transpose.
    iface_mon = dict(monitor_center=0.0, monitor_size=fdtd_solver.sim_width)
    in_kwargs = {"input_idx": 0, **iface_mon}

    cache: dict = {}
    smat_through = hybrid.solve(dc, WAVELENGTHS_DC, port="through", cache=cache,
                                 geometry_kwargs=[in_kwargs, {}, {"output_idx": 0, **iface_mon}])
    smat_cross = hybrid.solve(dc, WAVELENGTHS_DC, port="cross", cache=cache,
                               geometry_kwargs=[in_kwargs, {}, {"output_idx": 1, **iface_mon}])
    smat = smat_through
    smat._data.update(smat_cross._data)
    return smat


TARGETS = {
    "mmi_1x2_eme": (run_mmi_1x2_eme, ("out_top", "out_bot")),
    "mmi_1x2_fdtd": (run_mmi_1x2_fdtd, ("out_0", "out_1")),
    "mmi_1x2_hybrid": (run_mmi_1x2_hybrid, ("out_top", "out_bot")),
    "directional_coupler_hybrid": (run_directional_coupler_hybrid, ("through", "cross")),
}


def _load_golden(name: str) -> dict:
    path = GOLDEN_DIR / f"{name}.json"
    if not path.exists():
        pytest.fail(
            f"{path} does not exist yet — run "
            f"python scripts/generate_golden_values.py to generate it first."
        )
    with open(path) as f:
        return json.load(f)


def _assert_transmission_matches(actual: np.ndarray, expected: list, label: str) -> None:
    # Transmission is a linear 0..1 power ratio; combined rtol+atol tolerates
    # both small and near-zero reference values without a separate dB check.
    np.testing.assert_allclose(
        actual, expected, rtol=0.10, atol=0.03,
        err_msg=f"{label}: transmission drifted outside tolerance of the recorded golden value",
    )


@pytest.mark.parametrize("name", list(TARGETS))
def test_golden(name):
    run_fn, ports = TARGETS[name]
    smat = run_fn()
    assert smat.check_passive()
    golden = _load_golden(name)
    for port in ports:
        _assert_transmission_matches(smat.transmission(port=port), golden["transmission"][port], port)
