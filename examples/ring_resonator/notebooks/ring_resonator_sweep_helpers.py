"""Reassemble racetrack ring-resonator resonance spectra for arbitrary
L_C/L_SIDE, using only already-cached run_left_fdtd()/run_bend_fdtd()/
solve_neff()/solve_body_supermodes() results -- no new FDTD/EME simulation.

This works because L_C/L_SIDE only enter ring_resonator_hybrid_cascaded.py's
analytic post-assembly (BODY_LENGTH, L_OUTER), never the cached building
blocks: run_left_fdtd()/run_bend_fdtd()'s cache key has been L_C-independent
since the 2026-08-04 fix (see that script's _sim_params() docstring), and
solve_body_supermodes()/solve_neff() never referenced L_C/L_SIDE at all.

load_cached_blocks() loads the four building blocks once (cache_only=True,
so a genuinely missing cache entry raises CacheMissError rather than
silently launching a simulation); assemble() then reruns the rest of that
script's section 11-13 assembly logic for one (L_C, L_SIDE) pair, purely in
numpy. The parent script is loaded as a truncated, meep/meow-stubbed exec
(same technique used for the cache-migration script) since it is a
standalone script, not a package, and running it unmodified would both
require meep/meow and re-execute its own hardcoded L_C/L_SIDE pipeline.
"""
import os
import sys
import types

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
CASCADED_PY = os.path.join(_HERE, "..", "python", "ring_resonator_hybrid_cascaded.py")
CACHE_DIR = os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache")
_REPO_ROOT = os.path.join(_HERE, "..", "..", "..")  # BEEF/ -- so "from beef import ..." resolves

_TRUNCATE_MARKER = "# 11. Run the cascade"


def _load_parent_module():
    """Exec ring_resonator_hybrid_cascaded.py up to (not including) section
    11 -- every constant/function the assembly needs is defined before that
    point; section 11 onward is the script's own hardcoded single-point run,
    which this module replaces with its own cache_only calls per grid point."""
    with open(CASCADED_PY) as f:
        lines = f.readlines()
    cut = next(i for i, line in enumerate(lines) if _TRUNCATE_MARKER in line)
    src = "".join(lines[:cut])
    src = src.replace("import meep as mp", "mp = types.SimpleNamespace()")
    src = src.replace("import meow as mw", "mw = None")
    ns = {"__name__": "ring_resonator_hybrid_cascaded", "__file__": CASCADED_PY, "types": types}
    repo_root = os.path.abspath(_REPO_ROOT)
    added = repo_root not in sys.path
    if added:
        sys.path.insert(0, repo_root)
    try:
        exec(compile(src, CASCADED_PY, "exec"), ns)
    finally:
        if added:
            sys.path.remove(repo_root)
    return ns


def load_cached_blocks(save_dir=CACHE_DIR):
    """One-time load of every cached building block the sweep needs, plus
    the L_C/L_SIDE-independent pieces of the dense-grid assembly (b_d,
    b_corrected, neff_dense). cache_only=True on every call: a missing cache
    entry raises beef.cache.CacheMissError here rather than falling back to
    a new simulation."""
    ns = _load_parent_module()
    np_ = ns["np"]
    wl_fdtd = np_.linspace(ns["WL_MIN"], ns["WL_MAX"], ns["N_WL_FDTD"])
    wl_solve = np_.linspace(ns["WL_MIN"], ns["WL_MAX"], ns["N_WL_SOLVE"])
    wl_dense = np_.linspace(ns["WL_MIN"], ns["WL_MAX"], ns["N_WL_DENSE"])

    M_left_bus = ns["run_left_fdtd"]("bus", wl_fdtd, save_dir=save_dir, cache_only=True)
    M_left_ring = ns["run_left_fdtd"]("ring", wl_fdtd, save_dir=save_dir, cache_only=True)
    M_left = np.stack([M_left_bus, M_left_ring], axis=-1)      # [n_wl, N_MODES, 2]
    M_right = M_left.transpose(0, 2, 1)                        # reciprocity reuse

    b_raw = ns["run_bend_fdtd"](wl_fdtd, save_dir=save_dir, cache_only=True)

    body_modes = [ns["solve_body_supermodes"](wl, save_dir=save_dir, cache_only=True) for wl in wl_solve]
    neff_even_fit = np.polyfit(wl_solve, [m[0] for m in body_modes], deg=2)
    neff_odd_fit = np.polyfit(wl_solve, [m[1] for m in body_modes], deg=2)

    neff_solve = [ns["solve_neff"](wl, save_dir=save_dir, cache_only=True) for wl in wl_solve]
    neff_fit = np.polyfit(wl_solve, neff_solve, deg=2)
    neff_dense = np.polyval(neff_fit, wl_dense)

    interp_complex = ns["interp_complex"]
    b_d = interp_complex(wl_fdtd, b_raw, wl_dense)
    phi_bend_correction = 2 * np.pi * neff_dense * (ns["BEND_D_IN"] + ns["BEND_D_OUT"]) / wl_dense
    b_corrected = b_d * np.exp(-1j * phi_bend_correction)   # unused unless headroom<0; kept for parity

    return {
        "ns": ns,
        "wl_fdtd": wl_fdtd, "wl_dense": wl_dense,
        "M_left": M_left, "M_right": M_right,
        "neff_even_fit": neff_even_fit, "neff_odd_fit": neff_odd_fit,
        "neff_dense": neff_dense,
        "b_d": b_d, "b_corrected": b_corrected,
    }


def assemble(blocks, L_C, L_SIDE):
    """Recompute the L_C/L_SIDE-dependent half of the assembly (BODY_LENGTH,
    S_coupler, L_OUTER, the final through-port spectrum) for one grid point,
    reusing every cached FDTD/mode-solve result untouched. Returns a dict
    with "valid" False (plus "reason") when this (L_C, L_SIDE) pair isn't
    physically buildable (BODY_LENGTH<=0, or BEND_STRAIGHT_HEADROOM<0 the
    way ring_resonator_hybrid_cascaded.py itself guards against -- see its
    own comment: "L_SIDE=1.0, L_FAR=2.0 at the shipped STUB_LEN already
    fails") or no resonance dips are found."""
    ns = blocks["ns"]
    S_MON, STUB_LEN, MONITOR_MARGIN = ns["S_MON"], ns["STUB_LEN"], ns["MONITOR_MARGIN"]
    IFACE_BUFFER, N_MODES = ns["IFACE_BUFFER"], ns["N_MODES"]
    interp_complex, find_dips = ns["interp_complex"], ns["find_dips"]

    wl_fdtd, wl_dense = blocks["wl_fdtd"], blocks["wl_dense"]
    M_left, M_right = blocks["M_left"], blocks["M_right"]
    neff_even_fit, neff_odd_fit = blocks["neff_even_fit"], blocks["neff_odd_fit"]
    neff_dense = blocks["neff_dense"]
    b_d = blocks["b_d"]

    # BODY_LENGTH(L_C) -- mirrors IFACE_X/IFACE_PLANE_L/R/BODY_LENGTH, lines 293-319
    iface_x = -L_C / 2 + IFACE_BUFFER
    body_length = -2 * (iface_x - MONITOR_MARGIN)
    if body_length <= 0:
        return {"valid": False, "reason": "BODY_LENGTH <= 0 for this L_C"}

    phi_even = 2 * np.pi * np.polyval(neff_even_fit, wl_fdtd) / wl_fdtd * body_length
    phi_odd = 2 * np.pi * np.polyval(neff_odd_fit, wl_fdtd) / wl_fdtd * body_length
    P_body = np.zeros((len(wl_fdtd), N_MODES, N_MODES), dtype=complex)
    P_body[:, 0, 0] = np.exp(1j * phi_even)
    P_body[:, 1, 1] = np.exp(1j * phi_odd)
    S_coupler = np.einsum("wij,wjk->wik", M_right, np.einsum("wij,wjk->wik", P_body, M_left))
    t, k1, k2, t_ring = (S_coupler[:, 0, 0], S_coupler[:, 1, 0],
                         S_coupler[:, 0, 1], S_coupler[:, 1, 1])

    t_d = interp_complex(wl_fdtd, t, wl_dense)
    k1_d = interp_complex(wl_fdtd, k1, wl_dense)
    k2_d = interp_complex(wl_fdtd, k2, wl_dense)
    t_ring_d = interp_complex(wl_fdtd, t_ring, wl_dense)

    # L_OUTER(L_C, L_SIDE) -- mirrors lines 385-421. BEND_STRAIGHT_HEADROOM<0
    # would require the phase-strip branch (b_corrected); this sweep only
    # covers the headroom>=0 region and marks the rest invalid instead.
    bend_straight_total = 2 * (ns["BEND_D_IN"] + ns["BEND_D_OUT"])
    l_outer_shipped = 2 * (L_SIDE - S_MON) + L_C
    headroom = l_outer_shipped - bend_straight_total
    if headroom < 0:
        return {"valid": False, "reason": "BEND_STRAIGHT_HEADROOM < 0 for this L_C/L_SIDE"}

    l_outer = l_outer_shipped - bend_straight_total
    phi = 2 * np.pi * neff_dense * l_outer / wl_dense
    P = (b_d ** 2) * np.exp(1j * phi)

    E_thru = t_d + k1_d * k2_d * P / (1.0 - t_ring_d * P)
    T = np.abs(E_thru) ** 2

    dips = find_dips(wl_dense, T)
    if len(dips) < 2:
        return {"valid": False, "reason": f"only {len(dips)} resonance dip(s) found", "T": T, "dips": dips}

    fsr = float(np.mean(np.diff([d[0] for d in dips])))
    wl_c = 0.5 * (wl_dense[0] + wl_dense[-1])
    center_wl, T_min, fwhm, q = min(dips, key=lambda d: abs(d[0] - wl_c))
    baseline = float(np.median(T))
    extinction_db = 10 * np.log10(baseline / T_min) if T_min > 0 else np.inf

    return {
        "valid": True, "T": T, "dips": dips,
        "fsr": fsr, "fwhm": float(fwhm), "extinction_db": float(extinction_db),
        "Q": float(q), "center_wl": float(center_wl),
    }
