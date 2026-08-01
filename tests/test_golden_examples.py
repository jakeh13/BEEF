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

from beef import (
    MMI, DirectionalCoupler, EMESolver, FDTDSolver, HybridSolver, SMatrix, SOI_160, SOI_220,
)

pytestmark = pytest.mark.slow

GOLDEN_DIR = pathlib.Path(__file__).parent / "golden"

WAVELENGTHS_MMI = np.linspace(1.26, 1.36, 3)
WAVELENGTHS_DC = np.linspace(1.50, 1.60, 3)
# The 2×2 is a C-band device (its examples sweep 1.530–1.565), unlike the 1×2's
# O-band. Kept to the device's own band rather than reusing WAVELENGTHS_MMI —
# a golden run at the wrong band would pin numbers no example produces.
WAVELENGTHS_MMI_2X2 = np.linspace(1.530, 1.565, 3)


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


def _mmi_2x2() -> MMI:
    """Dimensions from `examples/mmi_2x2/python/`, unreduced — geometry is the
    one thing these replicas keep at full size."""
    return MMI(
        n_inputs=2, n_outputs=2,
        mm_length=28.0, mm_width=2.8,
        waveguide_width=0.45, taper_length=4.0, taper_width=1.3,
        waveguide_offset=1.5, waveguide_height=0.22,
    )


# Mode count for the 2×2 replicas. Production carries 14; 6 is the reduction,
# and it is chosen conservatively rather than scaled. Since 2026-07-31
# `_iter_distinct_bands()` RAISES when a cross-section cannot supply the
# requested guided modes (`fdtd_vs_EME_mode_validation.md` item 14) instead of
# silently padding, so an over-request at the reduced res=20 would fail the
# generation outright rather than degrade. 6 is safe by comparison, not by
# guess: `run_mmi_1x2_hybrid` already runs 6 at res=20 on a *narrower and
# thinner* core (2.0 µm × 160 nm vs this device's 2.8 µm × 220 nm), which
# supports strictly fewer guided modes.
_MMI_2X2_N_MODES = 6


def run_mmi_2x2_eme() -> SMatrix:
    mmi = _mmi_2x2()
    solver = EMESolver(
        # Same five-region split as mmi_2x2_eme.py (wg | taper | body | taper |
        # wg), mode counts reduced from the example's [1, 7, 15, 7, 1].
        num_modes=[1, 4, _MMI_2X2_N_MODES, 4, 1], num_cells=[1, 4, 1, 4, 1],
        lengths=[0.5, 4.0, 28.0, 4.0, 0.5],
        oxide_width=5.0, oxide_height=1.5, res=20.0,
        material_stack=SOI_220,
        # The example uses mode_solver="tidy3d" (it has no FDTD section to stay
        # consistent with). These replicas use "meep" throughout — same choice
        # run_mmi_1x2_eme makes — so a golden never depends on a second mode
        # solver being installed. What this pins is the cascade, not which
        # eigensolver produced the basis.
        mode_solver="meep",
    )
    # Both inputs and both outputs sit waveguide_offset = 1.5 µm apart, so the
    # cross-sections must be restricted to one branch at each end or the mode
    # solver sees a two-guide supermode instead of an isolated waveguide —
    # mmi_2x2_eme.py's own docstring makes this point. Excite the top input
    # (input_idx=1); bar is the top output (output_idx=1), cross the bottom (0).
    smat_bar = solver.solve(mmi, WAVELENGTHS_MMI_2X2, ports=["bar"],
                            input_idx=1, output_idx=1)
    smat_cross = solver.solve(mmi, WAVELENGTHS_MMI_2X2, ports=["cross"],
                              input_idx=1, output_idx=0)
    smat = smat_bar
    smat._data.update(smat_cross._data)
    return smat


def run_mmi_2x2_fdtd() -> SMatrix:
    mmi = _mmi_2x2()
    solver = FDTDSolver(
        res=20.0, sim_width=5.0, sim_height=1.5, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        num_input_modes=1, num_output_modes=1,
        # out_0 = top half (bar), out_1 = bottom half (cross) — the order the
        # example's own comment fixes, and what the port names below assume.
        output_monitor_centers=[1.25, -1.25], output_monitor_widths=[2.5, 2.5],
        material_stack=SOI_220,
        # 1e-9 rather than the example's 1e-11: same reduction the other
        # single-domain golden (run_mmi_1x2_fdtd) takes, since a golden is a
        # regression pin at a loose tolerance, not a physics measurement.
        decay_tol=1e-9,
    )
    return solver.solve(mmi, WAVELENGTHS_MMI_2X2, port="top", input_mode=1)


def run_mmi_2x2_hybrid() -> SMatrix:
    mmi = _mmi_2x2()
    n_modes = _MMI_2X2_N_MODES
    fdtd_common = dict(
        # dpml doubled vs the example's 0.5 — same reasoning as
        # run_mmi_1x2_hybrid's and the coupler's: these are short interior
        # solve_section() spans, and at res=20 (half production's 40) 0.5 µm of
        # PML is half the pixel count production gets.
        res=20.0, sim_width=5.0, sim_height=1.5, dpml=1.0,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=SOI_220, decay_tol=1e-9,
    )
    # BOTH end sections need half_domain_aperture here, unlike the 1×2 whose
    # lone input keeps the full sim_width: this device has two ports 1.5 µm
    # apart at *each* end, and all branches are always physically drawn, so the
    # aperture is what structurally excludes the neighbour.
    input_solver = FDTDSolver(half_domain_aperture=True, **fdtd_common)
    body_solver = EMESolver(
        num_modes=[n_modes], num_cells=[1], lengths=[1.0],
        oxide_width=5.0, oxide_height=1.5, res=20.0,
        material_stack=SOI_220, mode_solver="meep",
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

    # NOTE — this replica deliberately does NOT reuse the example's
    # mirror-symmetry shortcut. `mmi_2x2_hybrid.py` solves the input section
    # once per input arm and *transposes* to get the output-side matrices,
    # which is a real and verified saving on a symmetric device but is the
    # example's own assembly, not `HybridSolver.solve()`. What a golden exists
    # to pin is `solve()` on its defaults — the same reasoning
    # run_directional_coupler_hybrid records for calibrate_body_neff — so each
    # section here is actually simulated. Expect this to cost more than the
    # example does and to give slightly different numbers; that is the point,
    # and it is why the example is the physical answer while this is the
    # regression pin.
    #
    # Top input excited (input_idx=1); bar = top output (output_idx=1), cross =
    # bottom (output_idx=0), matching the example's docstring exactly.
    in_kwargs = {"input_idx": 1}
    cache: dict = {}
    smat_bar = hybrid.solve(mmi, WAVELENGTHS_MMI_2X2, port="bar", cache=cache,
                            geometry_kwargs=[in_kwargs, {}, {"output_idx": 1}])
    smat_cross = hybrid.solve(mmi, WAVELENGTHS_MMI_2X2, port="cross", cache=cache,
                              geometry_kwargs=[in_kwargs, {}, {"output_idx": 0}])
    smat = smat_bar
    smat._data.update(smat_cross._data)
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
    # Same convention as the example (and as the MMI replicas above): each FDTD
    # section folds body_stub µm of the straight coupling region into itself,
    # past the S-bend/coupling-region discontinuity, so the interface mode
    # decomposition lands on a uniform coupled-pair cross-section. HybridSolver
    # takes the EME body's length from the boundaries, so the three sections
    # still cover exactly coupling_length of straight between them.
    body_stub = 0.25
    boundaries = [-half_total, -half_body + body_stub, half_body - body_stub, half_total]
    # NOTE — this replica deliberately does NOT calibrate the body's β, while
    # examples/.../directional_coupler_hybrid.py now does (its CALIBRATE_BODY,
    # and dc_validation.md items 6/14). That is not drift. What this golden
    # exists to pin is `HybridSolver.solve()` on its *defaults*, and
    # calibrate_body_neff defaults off; the two therefore give different
    # absolute numbers on purpose — this is the regression pin, the example is
    # the physically-calibrated answer. If that default ever flips, this replica
    # and directional_coupler_hybrid.json must be regenerated together.
    hybrid = HybridSolver(
        solvers=[fdtd_solver, body_solver, fdtd_solver],
        boundaries=boundaries, interface_modes=[n_modes, n_modes],
    )
    # Pin the interior interface monitor to the full non-PML domain, centered on
    # the shared y=0 axis (matches the EME body's own supermode basis) — same
    # pin examples/directional_coupler/python/directional_coupler_hybrid.py
    # applies manually via _iface_mon; passed on both the input and output
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
    # Added 2026-08-01 — mmi_2x2_validation.md gap 1. This device had no golden
    # of any kind, and the gap was not just "not generated": there was no
    # runner to generate one from. Writing them is unblocked work; *running*
    # them still waits on the outstanding cold regeneration pass, and doing it
    # now means that pass covers all three devices in one go rather than
    # needing a second round for this one.
    "mmi_2x2_eme": (run_mmi_2x2_eme, ("bar", "cross")),
    "mmi_2x2_fdtd": (run_mmi_2x2_fdtd, ("out_0", "out_1")),
    "mmi_2x2_hybrid": (run_mmi_2x2_hybrid, ("bar", "cross")),
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
