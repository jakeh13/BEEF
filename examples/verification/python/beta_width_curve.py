"""Turn a measured beta(width) curve into a HybridSolver `cell_neff_fn`.

Companion to measure_beta_width_curve.py. That script measures n_eff on
UNIFORM strips across the width range a taper sweeps; this one interpolates it
per EME cell so a tapered body can be calibrated cell by cell via
`EMESolver.solve_section(cell_neff_override=...)`.

WHY INTERPOLATION IS LEGITIMATE HERE. Each EME cell IS a uniform waveguide of
a fixed width -- that is what the staircase decomposition means -- so "the
n_eff of this cell" is exactly "the n_eff of a uniform strip at this width",
which is the quantity the curve measures. The interpolation is over width, a
smooth geometric parameter, not over anything modal.

WHAT IT CANNOT FIX. Calibration corrects the PHASE each cell accumulates. It
does not touch the interface S-matrices (by construction -- those depend on
the mode fields) and so does nothing about the ~2% per-boundary staircase
discontinuity measured in mmi_2x2_waist_body_validation.md. Those are separate
error sources and this addresses only the first.
"""
import json

import numpy as np


def load_curve(path):
    """Load a beta(width) curve written by measure_beta_width_curve.py."""
    with open(path) as fh:
        d = json.load(fh)
    if not d.get("complete", False):
        raise ValueError(
            f"{path} is a PARTIAL curve ({len(d['rows'])} of "
            f"{len(d['widths_requested'])} widths). Calibrating from it would "
            "extrapolate silently over the missing range -- finish the "
            "measurement first.")
    widths = np.array([r["width"] for r in d["rows"]], dtype=float)
    neff = np.array([r["neff"] for r in d["rows"]], dtype=float)  # [w, wl, mode]
    order = np.argsort(widths)
    return dict(widths=widths[order], neff=neff[order],
                wavelengths=np.array(d["wavelengths"], dtype=float),
                config=d["config"])


def neff_at(curve, width, wavelengths, n_modes):
    """Interpolated n_eff ``[len(wavelengths), n_modes]`` at one width."""
    w_grid = curve["widths"]
    if width < w_grid[0] - 1e-9 or width > w_grid[-1] + 1e-9:
        raise ValueError(
            f"width {width:.4f} um is outside the measured curve "
            f"[{w_grid[0]:.4f}, {w_grid[-1]:.4f}]. Extrapolating a calibration "
            "is how a plausible wrong answer gets made -- re-measure the curve "
            "to cover this taper instead.")
    if curve["neff"].shape[2] < n_modes:
        raise ValueError(
            f"curve carries {curve['neff'].shape[2]} modes, section needs "
            f"{n_modes}")
    # Each requested wavelength must appear in the curve's own grid. A SUBSET
    # is fine and expected -- a sharded run solves one wavelength of eleven --
    # but a wavelength the curve never measured is not, because n_eff is
    # dispersive and interpolating across the band would calibrate it onto a
    # neighbour's beta. The 4x4's calibration study measured that band centre
    # shifts n_eff by 0.0017, which over a long body is ~74 deg of phase: real
    # and second-order, but far too big to absorb silently.
    wl_req = np.atleast_1d(np.asarray(wavelengths, dtype=float))
    wl_curve = curve["wavelengths"]
    idx = []
    for wl in wl_req:
        hits = np.flatnonzero(np.abs(wl_curve - wl) < 1e-9)
        if hits.size == 0:
            raise ValueError(
                f"wavelength {wl:.6f} um is not in the calibration curve's "
                f"grid [{wl_curve[0]:.4f} .. {wl_curve[-1]:.4f}, "
                f"{len(wl_curve)} points]. n_eff is dispersive, so calibrating "
                "it from a neighbouring wavelength would be a silent error -- "
                "re-measure the curve on this solve's own grid.")
        idx.append(int(hits[0]))

    out = np.empty((len(wl_req), n_modes), dtype=float)
    for out_i, wl_i in enumerate(idx):
        for m in range(n_modes):
            out[out_i, m] = np.interp(width, w_grid, curve["neff"][:, wl_i, m])
    return out


def cell_widths(geometry, x0, x1, n_cells, **geometry_kwargs):
    """Width at each cell's CENTRE for a section split into `n_cells`.

    Centres, not edges: a cell's single cross-section is solved at its centre
    in `mw.create_cells`, so calibrating on an edge width would offset the
    whole ramp by half a cell.
    """
    edges = np.linspace(x0, x1, n_cells + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return np.array([geometry.width_at(float(x), **geometry_kwargs)
                     for x in centres])


def well_ordered_modes(curve, tol=0.0):
    """How many leading modes keep a consistent identity across the curve.

    The eigensolver sorts by n_eff, so when two modes CROSS as the width
    changes they swap indices and "mode k" stops meaning one physical mode.
    Interpolating such a mode blends two of them. A crossing shows up as
    n_eff(width) going DOWN for a mode that should rise monotonically with
    width, so this returns the number of leading modes that are monotonically
    increasing in width at every wavelength.

    Measured on the waisted 2x2 curve: modes 0-5 are clean, 6/7 swap between
    w=2.66 and 2.73 um, and 8/9 are near-degenerate at w=2.45.
    """
    neff = curve["neff"]                       # [width, wl, mode]
    n_modes = neff.shape[2]
    for m in range(n_modes):
        d = np.diff(neff[:, :, m], axis=0)     # [width-1, wl]
        if not np.all(d > tol):
            return m
    return n_modes


def make_cell_neff_fn(curve, solvers, calibrate_sections=None, verbose=True,
                      max_mode=None):
    """Build the callable `HybridSolver(cell_neff_fn=...)` expects.

    Args:
        curve: from :func:`load_curve`.
        solvers: the same solver list handed to HybridSolver, used to read each
            EME section's own cell count.
        calibrate_sections: optional set of section indices to calibrate;
            None (default) calibrates every section that has a `num_cells`
            attribute, i.e. every EME section and no FDTD one.
        max_mode: calibrate only modes ``[0, max_mode)`` and emit NaN above,
            which `cell_neff_override` reads as "keep this mode's own solved
            n_eff". None (default) uses :func:`well_ordered_modes`, i.e.
            calibrates exactly the leading modes whose identity is stable
            across the measured widths. Deliberately not all-or-nothing: a
            mode that crosses another mid-curve has no well-defined
            interpolated beta, and feeding it one silently corrupts that
            channel.

    The returned fn gives None for any section it should not touch, which is
    what HybridSolver treats as "leave this section alone".
    """
    reported = set()
    n_ok = well_ordered_modes(curve) if max_mode is None else int(max_mode)
    if verbose:
        total = curve["neff"].shape[2]
        print(f"  cell_neff_fn: calibrating modes 0..{n_ok - 1} of {total}; "
              f"modes {n_ok}..{total - 1} keep their solved n_eff "
              "(identity not stable across the measured widths)")

    def fn(k, geometry, wavelengths, n_modes, x0, x1):
        solver = solvers[k]
        if not hasattr(solver, "num_cells"):
            return None                       # an FDTD port section
        if calibrate_sections is not None and k not in calibrate_sections:
            return None
        n_cells = solver.num_cells[0]
        widths = cell_widths(geometry, x0, x1, n_cells)
        out = np.full((len(wavelengths), n_cells, n_modes), np.nan,
                      dtype=complex)
        m_cal = min(n_ok, n_modes)
        for ci, w in enumerate(widths):
            out[:, ci, :m_cal] = neff_at(curve, float(w), wavelengths,
                                         n_modes)[:, :m_cal]
        if verbose and k not in reported:
            reported.add(k)
            print(f"  cell_neff_fn: section {k} calibrated, {n_cells} cells, "
                  f"widths {widths.min():.3f}..{widths.max():.3f} um, "
                  f"n_eff[0] {out[0, 0, 0].real:.4f}..{out[0, -1, 0].real:.4f}",
                  flush=True)
        return out

    return fn
