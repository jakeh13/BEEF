"""TEMPORARY — per-device interface-mode audit (validation_markdown check 2 + 3).

Delete with the rest of examples/verification/ once the validation work is
signed off.

WHAT THIS ANSWERS
-----------------
Two questions, per device, on the cross-section its HybridSolver chain actually
decomposes light onto:

  1. **How many of the requested interface modes are guided?** The directional
     coupler asks for N_MODES = 6 and one of them (n_eff 1.4367) sits *below*
     the 1.444 cladding index -- a radiation mode, not a physical channel.
     mmi_1x2 asks for 12, mmi_2x2 for 15, modal_multiplexer for 16 across four
     interfaces, all on wide slab bodies rather than a coupled pair. If a large
     fraction of those bases are unbound, the interface decomposition is mostly
     projecting onto modes that are not there.

  2. **Do meep and meow return the same modes in the same order?** HybridSolver
     contracts index i of an FDTD interface against index i of an EME body
     assuming they match. They do for the coupler (identity permutation, worst
     guided-band |dn_eff| 0.0023). Mode crowding grows with count, so the
     devices above are more exposed, and none has been checked.

Plus, because it is nearly free once the modes are solved: EME mesh convergence
per device. Do NOT assume the coupler's converged res=80 transfers -- mmi_1x2
and modal_multiplexer are SOI_160, and a 160 nm core confines much more weakly
than 220 nm, so its modes reach further into the cladding.

CHEAP BY CONSTRUCTION
---------------------
FDTDSolver.section_mode_neff() uses meep's get_eigenmode(), which runs MPB on the
initialized grid and returns with no time-stepping at all. Cost is one
set_epsilon plus one eigensolve per mode -- minutes, against the hours a real
section solve needs. That is the whole reason this audit is affordable across
four devices.

Run:
    python .../audit_interface_modes.py                 # every device
    python .../audit_interface_modes.py mmi_1x2         # one device
    python .../audit_interface_modes.py --convergence   # add the EME mesh sweep

Output: example_output/verification/interface_mode_audit/<device>.json plus a
printed table per device.
"""

import argparse
import json
import os
import socket
import sys
import time

import numpy as np

import meow as mw
from beef import (
    MMI, DirectionalCoupler, EMESolver, FDTDSolver,
    SOI_160, SOI_220, extrude_polygons,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
output_dir = os.path.join(_ROOT, "example_output", "verification", "interface_mode_audit")

# The cladding index. A mode below this is not bound -- it lives on the
# discretization, and neither engine owes the other an answer about it.
# MaterialStack only exposes meep Medium objects, and those carry the index as
# epsilon_diag (n^2), so take it from there rather than re-reading meow's table
# and risking a different interpolation than the solvers actually used.
def n_cladding(stack, wl):
    return float(np.sqrt(np.real(stack.sio2_medium(wl).epsilon_diag.x)))


# ---------------------------------------------------------------------------
# Devices, pinned to what each example actually constructs. The section is a
# short cut inside the uniform body, well clear of any taper/bend seam, because
# that is the cross-section the interface decomposition lands on.
# ---------------------------------------------------------------------------
def _dc():
    dc = DirectionalCoupler(
        waveguide_width=0.45, gap=0.25, coupling_length=30.0, bend_length=9.0,
        waveguide_offset=3.0, bend_type="euler", waveguide_height=0.22,
    )
    return dict(
        geometry=dc, n_modes=6, wl=1.55, stack=SOI_220,
        section=(-1.0, 1.0), sim_width=6.0, sim_height=1.5, dpml=0.5, res=40.0,
        note="coupled pair; the worked example -- 6 modes, 1 unguided",
    )


def _mmi_1x2():
    # Pinned to examples/mmi_1x2/python/mmi_1x2_hybrid.py:58.
    mmi = MMI(
        n_inputs=1, n_outputs=2,
        mm_length=4.5, mm_width=2.0, waveguide_width=0.35,
        taper_length=1.0, taper_width=0.8, waveguide_offset=1.0,
        waveguide_height=0.16,
    )
    # mm_length 4.5 um -> body spans +/-2.25, so [-1, +1] is 1.25 um clear of
    # each taper seam.
    return dict(
        geometry=mmi, n_modes=12, wl=1.31, stack=SOI_160,
        section=(-1.0, 1.0), sim_width=5.0, sim_height=1.5, dpml=0.5, res=40.0,
        note="2.0 x 0.16 um slab body, mm_length 4.5 um, O-band",
    )


def _mmi_2x2():
    # Pinned to examples/mmi_2x2/python/mmi_2x2_hybrid.py:86 — a much longer and
    # wider body than the 1x2, on SOI_220 rather than SOI_160.
    mmi = MMI(
        n_inputs=2, n_outputs=2,
        mm_length=28.0, mm_width=2.8, waveguide_width=0.45,
        taper_length=4.0, taper_width=1.3, waveguide_offset=1.5,
        waveguide_height=0.22,
    )
    # Body spans +/-14, so [-1, +1] is 13 um clear of either seam.
    return dict(
        geometry=mmi, n_modes=15, wl=1.55, stack=SOI_220,
        section=(-1.0, 1.0), sim_width=5.0, sim_height=1.5, dpml=0.5, res=40.0,
        note="2.8 x 0.22 um slab body, mm_length 28 um, C-band",
    )


# modal_multiplexer is deliberately absent: it builds a CompoundGeometry of two
# MMIs plus a WaveguideRun connector (modal_multiplexer_hybrid.py:137), with four
# section interfaces rather than two, so "the interface cross-section" is not a
# single cut. It needs its own treatment -- see
# validation_markdown/modal_multiplexer_validation.md.
DEVICES = {
    "directional_coupler": _dc,
    "mmi_1x2": _mmi_1x2,
    "mmi_2x2": _mmi_2x2,
}


def _is_master():
    try:
        import meep as mp
        return bool(mp.am_really_master())
    except Exception:
        return True


MASTER = _is_master()


def say(*a):
    if MASTER:
        print(*a, flush=True)


def _atomic_write(path, text):
    tmp = f"{path}.tmp.{socket.gethostname()}.{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
def meow_neff(cfg, res, n_modes):
    eme = EMESolver(
        num_modes=[n_modes], num_cells=[1], lengths=[1.0],
        oxide_width=cfg["sim_width"], oxide_height=cfg["sim_height"],
        res=float(res), material_stack=cfg["stack"], mode_solver="meep",
    )
    x0, x1 = cfg["section"]
    L = x1 - x0
    gds = cfg["geometry"].section_gds(x0, x1, eme.oxide_width)
    cells = mw.create_cells(
        structures=extrude_polygons(gds, eme.material_stack.extrusion_rules),
        mesh=eme._build_mesh(), Ls=[L],
    )
    cs = mw.CrossSection.from_cell(
        cell=cells[0], env=mw.Environment(wl=float(cfg["wl"]), T=eme.temperature))
    return np.array([float(np.real(m.neff)) for m in eme.compute_modes(cs, n_modes)])


def meep_neff(cfg, res, n_modes):
    solver = FDTDSolver(
        res=float(res), sim_width=cfg["sim_width"], sim_height=cfg["sim_height"],
        dpml=cfg["dpml"], source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=cfg["stack"], decay_tol=1e-11,
    )
    x0, x1 = cfg["section"]
    return solver.section_mode_neff(
        cfg["geometry"], np.array([cfg["wl"]]), x0, x1, n_modes,
        monitor_center=0.0, monitor_size=cfg["sim_width"],
    )[0]


def correspondence(a, b):
    """Nearest-neighbour permutation of a onto b. Identity means index i names
    the same mode in both. ONLY valid at matched resolution -- across
    resolutions the offset can exceed the mode spacing and fake a reorder."""
    return [int(np.argmin(np.abs(b - v))) for v in a]


def audit(name, convergence=False):
    cfg = DEVICES[name]()
    n_modes, res, wl = cfg["n_modes"], cfg["res"], cfg["wl"]
    n_clad = n_cladding(cfg["stack"], wl)

    say(f"\n{'='*72}\n[audit] {name} — {cfg['note']}")
    say(f"[audit] N_MODES={n_modes}, res={res:.0f}, lambda={wl} um, "
        f"n_cladding={n_clad:.4f}, section={cfg['section']}")

    t0 = time.time()
    mp_neff = meep_neff(cfg, res, n_modes)
    mw_neff = meow_neff(cfg, res, n_modes)
    say(f"[audit] both engines solved in {(time.time()-t0)/60:.1f} min\n")

    guided = mp_neff > n_clad
    n_guided = int(guided.sum())
    say(f"[audit] band     meep      meow      diff   guided")
    for i in range(n_modes):
        say(f"[audit]  {i:3d}  {mp_neff[i]:9.6f} {mw_neff[i]:9.6f} "
            f"{mp_neff[i]-mw_neff[i]:+9.4f}   {'yes' if guided[i] else 'NO'}")

    perm = correspondence(mp_neff[guided], mw_neff[guided])
    identity = perm == list(range(n_guided))
    worst = float(np.max(np.abs(mp_neff[guided] - mw_neff[guided]))) if n_guided else float("nan")

    say(f"\n[audit] guided modes: {n_guided}/{n_modes}"
        f"{'  <-- basis contains unguided modes' if n_guided < n_modes else ''}")
    say(f"[audit] correspondence (guided, matched res): {perm}")
    say(f"[audit] ORDER PRESERVED: {identity}")
    say(f"[audit] worst guided |dn_eff|: {worst:.4f}")

    result = dict(
        device=name, n_modes_requested=n_modes, n_modes_guided=n_guided,
        resolution=res, wavelength=wl, n_cladding=n_clad,
        section=list(cfg["section"]), note=cfg["note"],
        neff_meep=mp_neff.tolist(), neff_meow=mw_neff.tolist(),
        correspondence=perm, order_preserved=bool(identity),
        worst_guided_dneff=worst,
    )

    if convergence:
        say(f"\n[audit] EME mesh convergence (guided modes only)")
        conv = {}
        prev = None
        for r in (40.0, 60.0, 80.0, 100.0):
            nf = meow_neff(cfg, r, n_modes)[:n_guided]
            conv[str(r)] = nf.tolist()
            shift = "" if prev is None else f"  max shift vs previous {np.max(np.abs(nf-prev)):.5f}"
            say(f"[audit]   res {r:5.0f}: n_eff[0]={nf[0]:.6f}{shift}")
            prev = nf
        result["eme_convergence"] = conv

    if MASTER:
        os.makedirs(output_dir, exist_ok=True)
        _atomic_write(os.path.join(output_dir, f"{name}.json"),
                      json.dumps(result, indent=2) + "\n")
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("device", nargs="?", choices=sorted(DEVICES), default=None)
    ap.add_argument("--convergence", action="store_true",
                    help="also sweep the EME mesh (adds ~4x the mode solves)")
    args = ap.parse_args()

    names = [args.device] if args.device else list(DEVICES)
    for n in names:
        audit(n, convergence=args.convergence)
    say(f"\n[audit] written -> {output_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
