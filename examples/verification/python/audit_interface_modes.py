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

modal_multiplexer is handled separately (see MULTI_DEVICES below): it is a
CompoundGeometry with FOUR section interfaces, not one cut, so "the interface
cross-section" is not a single (n_modes, section) pair — each of the four
carries its own N_MODES=16 basis and must be checked independently. A swap at
any one of the four corrupts that section's contraction silently, so all four
are audited, not just one.

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
    python .../audit_interface_modes.py                    # every device
    python .../audit_interface_modes.py mmi_1x2            # one device
    python .../audit_interface_modes.py modal_multiplexer  # all 4 interfaces
    python .../audit_interface_modes.py --convergence      # add the EME mesh sweep

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
    MMI, MMIPort, CompoundGeometry, WaveguideRun,
    DirectionalCoupler, EMESolver, FDTDSolver,
    SOI_160, SOI_220, SIN_300, extrude_polygons,
    best_correspondence, field_overlap_matrix, mesh_cell_centers, stack_meow_mode_fields,
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
        # N_MODES 12->8: mmi_1x2_validation.md gap 3, mmi_1x2_hybrid.py:95 --
        # bands 8/9 came back non-identity between engines at 12, the entire
        # problem region was 8-11 (near the SOI_160 cutoff), bands 0-7 agree
        # to <=1.8e-4. This audit must track what the example actually
        # carries, not the pre-fix count.
        geometry=mmi, n_modes=8, wl=1.31, stack=SOI_160,
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
    #
    # n_modes 15 -> 14 on 2026-07-31, so this tracks the device again: band 15
    # (n_eff 1.4078) is below the 1.4443 cladding index and `mmi_2x2_hybrid.py`
    # dropped it as part of fdtd_vs_EME_mode_validation.md item 9. The 15-mode
    # result this audit produced on 2026-07-26 (job 11500536: identity, worst
    # guided |dn_eff| 0.0007, "the cleanest device audited") therefore described
    # a request the device no longer makes -- which is what
    # mmi_2x2_validation.md gap 6 was opened to fix, along with the fact that
    # it predates both duplicate-mode guards.
    return dict(
        geometry=mmi, n_modes=14, wl=1.55, stack=SOI_220,
        section=(-1.0, 1.0), sim_width=5.0, sim_height=1.5, dpml=0.5, res=40.0,
        note="2.8 x 0.22 um slab body, mm_length 28 um, C-band",
    )


def _mmi_1x4():
    # Pinned to examples/mmi_1x4/python/mmi_1x4_hybrid.py exactly (SIN_300,
    # C-band, N_MODES=12 after the 20->12 prewarm-timeout reduction --
    # validation_markdown/mmi_1x4_sin300_validation.md step 6d).
    mmi = MMI(
        n_inputs=1, n_outputs=4,
        mm_length=206.0, mm_width=12.0, waveguide_width=1.1,
        taper_length=5.0, taper_width=2.2, waveguide_offset=2.4,
        waveguide_height=0.30, use_neighbor_aperture=True,
    )
    # Body half-length is 103 um; [-1, 1] sits well clear of either taper
    # seam, same "short cut inside the uniform body" convention the other
    # three devices use above.
    return dict(
        geometry=mmi, n_modes=12, wl=1.55, stack=SIN_300,
        section=(-1.0, 1.0), sim_width=20.0, sim_height=SIN_300.sz, dpml=0.5, res=40.0,
        # tidy3d, NOT meep -- the one thing that differs from every other
        # device's audit config. mmi_1x4_hybrid.py's own body EME sub-solver
        # uses mode_solver="tidy3d" (step 4's compute_modes_meep-silently-
        # wrong-for-weakly-confined-nitride finding), so auditing this device
        # with meep/MPB (the other three devices' real configuration) would
        # check a basis the real HybridSolver chain never actually builds.
        mode_solver="tidy3d",
        note="12.0 x 0.30 um SiN slab body, mm_length 206 um, C-band, tidy3d EME backend",
    )


DEVICES = {
    "directional_coupler": _dc,
    "mmi_1x2": _mmi_1x2,
    "mmi_2x2": _mmi_2x2,
    "mmi_1x4": _mmi_1x4,
}


def _modal_multiplexer():
    # Pinned to examples/modal_multiplexer/python/modal_multiplexer_hybrid.py
    # (geometry constants, BODY_STUB, N_MODES=16, boundaries formula) exactly
    # -- this must audit the same four cross-sections HybridSolver actually
    # cuts at, not a stand-in.
    TAPER_L = 2.0
    WO = 1.0375
    BETWEEN_WG_LENGTH = 3.0
    MM_WIDTH = 2.8
    BRANCH_SPECS = [
        (+WO, 0.33, 0.725),   # top
        (0.0, 0.35, 1.15),    # mid
        (-WO, 0.35, 0.725),   # bot
    ]

    def _branch_ports():
        return [MMIPort(y=y, waveguide_width=ww, taper_length=TAPER_L, taper_width=tw)
                for y, ww, tw in BRANCH_SPECS]

    mmi1 = MMI(
        n_inputs=1, n_outputs=3,
        mm_length=16.92, mm_width=MM_WIDTH,
        inputs=[MMIPort(y=0.0, waveguide_width=0.7, taper_length=TAPER_L, taper_width=1.4)],
        outputs=_branch_ports(),
        waveguide_height=0.16,
    )
    mmi2 = MMI(
        n_inputs=3, n_outputs=3,
        mm_length=34.24, mm_width=MM_WIDTH,
        inputs=_branch_ports(),
        outputs=_branch_ports(),
        waveguide_height=0.16,
        use_neighbor_aperture=True,
    )
    connector = WaveguideRun(
        length=BETWEEN_WG_LENGTH, ports=[(y, ww) for y, ww, _ in BRANCH_SPECS])
    mux = CompoundGeometry([mmi1, connector, mmi2])

    N_MODES = 16
    BODY_STUB = 0.75
    half_total = mux.total_length / 2
    body1_start = -half_total + TAPER_L
    body1_end = body1_start + mmi1.mm_length
    body2_start = body1_end + TAPER_L + BETWEEN_WG_LENGTH + TAPER_L
    body2_end = body2_start + mmi2.mm_length

    # section_mode_neff() unconditionally draws its own self-contained
    # stub_left=stub_right=dpml on BOTH sides of the query box (it has no
    # notion of "this side is an interior mode-matching interface, not real
    # PML" -- that distinction only exists in HybridSolver._stubs(), which
    # zeroes the interior side; see beef/solvers/hybrid.py's _stubs(),
    # verified directly).
    #
    # This audit queries each interface directly on the MMI it actually
    # belongs to, in that MMI's OWN local frame (x=0 at its longitudinal
    # centre), exactly like the mmi_1x2 / mmi_2x2 audits above already do,
    # rather than through the composed `mux`. Verified these four local
    # positions sit precisely BODY_STUB=0.75 inside their MMI's own taper
    # seams.
    #
    # Historical note (2026-07-30): this local-frame routing originally
    # existed to work around a stub double-application bug in
    # CompoundGeometry._overlaps(), which made a padded query box that sat
    # fully inside one sub-segment need 2*dpml=1.0 um of clearance instead
    # of dpml's 0.5 -- more than BODY_STUB's 0.75. **That bug is now fixed**
    # (see notes.md Geometry item 6 and
    # tests/test_geometry.py::TestCompoundGeometryStubIsAppliedExactlyOnce),
    # so going through `mux` directly would work today. Kept as-is anyway:
    # querying each MMI in its own frame is what the other two MMI audits
    # do, it keeps the four configs' geometry unambiguous, and there is no
    # reason to re-route working, already-run audit code.
    interface_local_xs = [
        ("FDTD input <-> EME MMI1", mmi1, -(mmi1.mm_length / 2 - BODY_STUB)),
        ("EME MMI1 <-> FDTD transition", mmi1, +(mmi1.mm_length / 2 - BODY_STUB)),
        ("FDTD transition <-> EME MMI2", mmi2, -(mmi2.mm_length / 2 - BODY_STUB)),
        ("EME MMI2 <-> FDTD output", mmi2, +(mmi2.mm_length / 2 - BODY_STUB)),
    ]
    # Half-span + dpml(0.5) must stay under BODY_STUB's 0.75 um clearance
    # from the taper seam; 0.1 leaves a 0.15 um margin.
    half_span = 0.1

    return [
        dict(
            geometry=geom, n_modes=N_MODES, wl=1.31, stack=SOI_160,
            section=(x_local - half_span, x_local + half_span),
            sim_width=5.0, sim_height=1.5, dpml=0.5, res=40.0,
            note=f"interface {i} ({label}) @ local x={x_local:.3f} of its own MMI",
        )
        for i, (label, geom, x_local) in enumerate(interface_local_xs)
    ]


# CompoundGeometry devices with more than one section interface: each entry
# is a list of per-interface configs (same shape as one DEVICES[...]() entry)
# rather than a single dict, because "the interface cross-section" is not one
# cut for these devices -- see the module docstring.
MULTI_DEVICES = {
    "modal_multiplexer": _modal_multiplexer,
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
    # mode_solver defaults to "meep" -- matched to what every existing
    # device's own hybrid script actually configures its body EME sub-solver
    # with (mmi_1x2_hybrid.py/mmi_2x2_hybrid.py: "matches FDTDSolver's own
    # eigenmode basis"). mmi_1x4 is the one exception (cfg sets "tidy3d") --
    # its own hybrid script deliberately does NOT use meep/MPB here (step 4's
    # compute_modes_meep-silently-wrong-for-weakly-confined-nitride finding),
    # so auditing it with the wrong backend would check a configuration the
    # real device never runs, not this device's actual basis.
    eme = EMESolver(
        num_modes=[n_modes], num_cells=[1], lengths=[1.0],
        oxide_width=cfg["sim_width"], oxide_height=cfg["sim_height"],
        res=float(res), material_stack=cfg["stack"],
        mode_solver=cfg.get("mode_solver", "meep"),
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


def overlap_correspondence(cfg, n_modes, n_cladding, cmp_res=20.0):
    """The beef.mode_identity field-overlap check, run on ONE cfg (see
    module docstring / README "Mode identity" section for why this
    generalizes past nearest-neighbour n_eff matching).

    Deliberately uses cmp_res, independent of cfg['res'] (the resolution
    the n_eff audit solved at) -- both engines are re-solved at this
    resolution just for the comparison grid, demonstrating the whole point
    of a shape-based check: it does not require the two engines to already
    agree numerically at a shared resolution the way nearest-neighbour
    n_eff matching does. Coarser than production (20 vs 40) is deliberate
    too: EigenmodeData.amplitude() is called once per (band, y, z,
    component) -- n_modes * len(y_grid) * len(z_grid) * 6 -- so halving the
    grid resolution cuts that 4x; a 0.05 um cell size still resolves a
    2.8 um-wide body's mode lobes (worst case here, mode 15 of 16, is ~6
    lobes across 2.8 um) at a handful of points per lobe, plenty for an
    aggregate overlap integral even though it would be marginal for
    resolving the field shape precisely.
    """
    x0, x1 = cfg["section"]
    wl = float(cfg["wl"])

    mesh = mw.Mesh2D(
        x=np.linspace(-cfg["sim_width"] / 2, cfg["sim_width"] / 2,
                      int(cmp_res * cfg["sim_width"]) + 1),
        y=np.linspace(-cfg["sim_height"] / 2, cfg["sim_height"] / 2,
                      int(cmp_res * cfg["sim_height"]) + 1),
    )
    y_grid, z_grid = mesh_cell_centers(mesh)

    gds = cfg["geometry"].section_gds(x0, x1, cfg["sim_width"])
    cells = mw.create_cells(
        structures=extrude_polygons(gds, cfg["stack"].extrusion_rules),
        mesh=mesh, Ls=[x1 - x0],
    )
    cs = mw.CrossSection.from_cell(cell=cells[0], env=mw.Environment(wl=wl, T=25.0))
    meow_modes = mw.compute_modes(cs, n_modes)
    meow_neff_arr = np.array([float(np.real(m.neff)) for m in meow_modes])
    meow_fields = stack_meow_mode_fields(meow_modes)

    solver = FDTDSolver(
        res=cmp_res, sim_width=cfg["sim_width"], sim_height=cfg["sim_height"],
        dpml=cfg["dpml"], source_margin=0.0, monitor_margin=0.25, three_d=True,
        material_stack=cfg["stack"], decay_tol=1e-11,
    )
    result = solver.section_mode_fields(
        cfg["geometry"], wl, x0, x1, n_modes, y_grid, z_grid,
        monitor_center=0.0, monitor_size=cfg["sim_width"],
    )
    meep_neff_arr = result["neff"]
    meep_fields = {k: result[k] for k in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")}

    guided_a = meep_neff_arr > n_cladding
    guided_b = meow_neff_arr > n_cladding
    fields_a = {k: v[guided_a] for k, v in meep_fields.items()}
    fields_b = {k: v[guided_b] for k, v in meow_fields.items()}

    dy = float(y_grid[1] - y_grid[0])
    dz = float(z_grid[1] - z_grid[0])
    eta = field_overlap_matrix(fields_a, fields_b, dy=dy, dz=dz)
    corr, margin = best_correspondence(eta)

    return dict(
        cmp_res=cmp_res, n_guided_meep=int(guided_a.sum()), n_guided_meow=int(guided_b.sum()),
        neff_meep_guided=meep_neff_arr[guided_a].tolist(), neff_meow_guided=meow_neff_arr[guided_b].tolist(),
        correspondence=corr, margin=margin.tolist(),
        order_preserved=bool(corr == list(range(len(corr)))),
        worst_overlap_on_diagonal=float(np.min([eta[i, corr[i]] for i in range(len(corr)) if corr[i] >= 0])) if corr else float("nan"),
    )


def _audit_one(label, cfg, convergence=False):
    """Runs both engines on one (geometry, section, n_modes) config and
    returns the result dict. Shared by single-cut devices (DEVICES) and each
    interface of a multi-cut device (MULTI_DEVICES)."""
    n_modes, res, wl = cfg["n_modes"], cfg["res"], cfg["wl"]
    n_clad = n_cladding(cfg["stack"], wl)

    say(f"\n{'='*72}\n[audit] {label} — {cfg['note']}")
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
        label=label, n_modes_requested=n_modes, n_modes_guided=n_guided,
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

    return result


def audit(name, convergence=False):
    cfg = DEVICES[name]()
    result = dict(device=name, **_audit_one(name, cfg, convergence=convergence))
    if MASTER:
        os.makedirs(output_dir, exist_ok=True)
        _atomic_write(os.path.join(output_dir, f"{name}.json"),
                      json.dumps(result, indent=2) + "\n")
    return result


def audit_multi(name, convergence=False):
    """Multi-interface devices (see MULTI_DEVICES): audits every interface
    independently and writes one JSON with a per-interface list plus a
    device-wide summary, since HybridSolver needs ALL of them sound, not
    just one."""
    cfgs = MULTI_DEVICES[name]()
    interfaces = [
        _audit_one(f"{name}[{i}]", cfg, convergence=convergence)
        for i, cfg in enumerate(cfgs)
    ]

    all_order_preserved = all(r["order_preserved"] for r in interfaces)
    worst_overall = max(r["worst_guided_dneff"] for r in interfaces)
    say(f"\n{'='*72}\n[audit] {name} — SUMMARY across {len(interfaces)} interfaces")
    for i, r in enumerate(interfaces):
        flag = "OK" if r["order_preserved"] else "MISMATCH"
        say(f"[audit]   [{i}] {r['note']}: order_preserved={r['order_preserved']} "
            f"({flag}), worst |dn_eff|={r['worst_guided_dneff']:.4f}, "
            f"guided={r['n_modes_guided']}/{r['n_modes_requested']}")
    say(f"[audit] ALL INTERFACES ORDER-PRESERVED: {all_order_preserved}")
    say(f"[audit] worst |dn_eff| over all interfaces: {worst_overall:.4f}")

    result = dict(
        device=name, all_order_preserved=bool(all_order_preserved),
        worst_guided_dneff_overall=worst_overall, interfaces=interfaces,
    )
    if MASTER:
        os.makedirs(output_dir, exist_ok=True)
        _atomic_write(os.path.join(output_dir, f"{name}.json"),
                      json.dumps(result, indent=2) + "\n")
    return result


def main():
    all_names = sorted(list(DEVICES) + list(MULTI_DEVICES))
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("device", nargs="?", choices=all_names, default=None)
    ap.add_argument("--convergence", action="store_true",
                    help="also sweep the EME mesh (adds ~4x the mode solves)")
    args = ap.parse_args()

    names = [args.device] if args.device else all_names
    for n in names:
        if n in MULTI_DEVICES:
            audit_multi(n, convergence=args.convergence)
        else:
            audit(n, convergence=args.convergence)
    say(f"\n[audit] written -> {output_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
