"""TEMPORARY -- bisects the unexplained meep<->meow overlap discrepancy on the
modal multiplexer's interface 0 (validation_markdown/
fdtd_vs_EME_mode_validation.md to-do item 1).

The disagreement: HybridSolver(check_mode_correspondence=True) reports a worst
diagonal overlap of 0.0000 at that interface, while
audit_interface_modes.py's own overlap_correspondence() measures 0.9735 on
what should be the same physical comparison. One of the two is wrong and it is
not yet known which -- so neither number is currently evidence of anything.

Strategy: run the LAYER-4 same-engine controls first (meow-vs-meow,
meep-vs-meep). Those must come back as the identity by construction, so they
localize the fault to the measurement before any cross-engine number is
interpreted. Then vary, one at a time, each difference between the two code
paths:

    A. geometry queried      : CompoundGeometry(mux) vs mmi1 in its own frame
    B. meow entry point      : EMESolver.section_modes (compute_modes_meep)
                               vs mw.compute_modes (audit's call)
    C. section half-span     : 0.05 (the check) vs 0.1 (the audit)

Everything is reported through beef.mode_identity.overlap_diagnostics(), i.e.
the full matrix -- diagonal, per-row completeness, off-diagonal mass, margin
-- not a single scalar. Both false alarms found on 2026-07-30 were diagnosed
from the SHAPE of the failure across modes, which a scalar hides.

Deliberately batch, not a login-node smoke run: each 16-mode meow solve at
res 40 on a 201x61 mesh is minutes, and an earlier attempt at two of them
died on a 900 s timeout. See batch/verification/run_diagnose_mode_correspondence.batch.
"""
from __future__ import annotations

import numpy as np

from beef import (
    MMI, MMIPort, WaveguideRun, CompoundGeometry, EMESolver, FDTDSolver, SOI_160,
    extrude_polygons, mesh_cell_centers, stack_meow_mode_fields,
    field_overlap_matrix, overlap_diagnostics, self_overlap_diagnostics,
)

WL = 1.31
RES = 40.0
N_MODES = 16
BODY_STUB = 0.75
SIM_W, SIM_H = 5.0, 1.5

TAPER_L, WO, BETWEEN, MM_WIDTH = 2.0, 1.0375, 3.0, 2.8
BRANCH_SPECS = [(+WO, 0.33, 0.725), (0.0, 0.35, 1.15), (-WO, 0.35, 0.725)]


def _branch_ports():
    return [MMIPort(y=y, waveguide_width=ww, taper_length=TAPER_L, taper_width=tw)
            for y, ww, tw in BRANCH_SPECS]


def _build():
    mmi1 = MMI(n_inputs=1, n_outputs=3, mm_length=16.92, mm_width=MM_WIDTH,
               inputs=[MMIPort(y=0.0, waveguide_width=0.7, taper_length=TAPER_L,
                               taper_width=1.4)],
               outputs=_branch_ports(), waveguide_height=0.16)
    mmi2 = MMI(n_inputs=3, n_outputs=3, mm_length=34.24, mm_width=MM_WIDTH,
               inputs=_branch_ports(), outputs=_branch_ports(),
               waveguide_height=0.16, use_neighbor_aperture=True)
    connector = WaveguideRun(length=BETWEEN, ports=[(y, ww) for y, ww, _ in BRANCH_SPECS])
    return mmi1, mmi2, CompoundGeometry([mmi1, connector, mmi2])


def _mesh():
    import meow as mw
    return mw.Mesh2D(x=np.linspace(-SIM_W / 2, SIM_W / 2, int(RES * SIM_W) + 1),
                     y=np.linspace(-SIM_H / 2, SIM_H / 2, int(RES * SIM_H) + 1))


def _report(label, diag, neff_a=None, neff_b=None):
    print(f"\n--- {label} ---", flush=True)
    print(f"  order_preserved        : {diag['order_preserved']}")
    print(f"  worst diagonal         : {diag['worst_diagonal']:.4f}")
    print(f"  worst row completeness : {diag['worst_row_completeness']:.4f}")
    print(f"  worst col completeness : {diag['worst_col_completeness']:.4f}")
    if "has_duplicate" in diag:
        # has_duplicate is the criterion that GATES a comparison: an
        # off-diagonal near 1 means two entries are one mode, so the basis is
        # rank-deficient. A merely-non-orthogonal basis (max off-diag ~0.13-0.19
        # here, from genuinely near-degenerate pairs on a symmetric
        # cross-section) is fine to compare -- reporting that as a failure was
        # an over-strict early criterion that would switch the check off on
        # exactly the devices it exists for.
        print(f"  HAS DUPLICATE (gate)   : {diag['has_duplicate']}  "
              f"(orthogonal: {diag['is_orthogonal']})  "
              f"(max off-diag {diag['max_off_diagonal']:.4f})")
    print(f"  {'i':>3} {'->j':>4} {'diag':>8} {'rowcomp':>8} {'offdiag':>8} {'margin':>8}"
          + ("  " + f"{'neff_a':>9} {'neff_b':>9}" if neff_a is not None else ""))
    for i, j in enumerate(diag["correspondence"]):
        row = (f"  {i:3d} {j:4d} {diag['diagonal'][i]:8.4f} "
               f"{diag['row_completeness'][i]:8.4f} {diag['off_diagonal_max'][i]:8.4f} "
               f"{diag['margin'][i]:8.4f}")
        if neff_a is not None and neff_b is not None and 0 <= j < len(neff_b):
            row += f"  {neff_a[i]:9.5f} {neff_b[j]:9.5f}"
        print(row, flush=True)


def _meow_via_emesolver(geometry, x_center, half_span, mesh, guard=True):
    eme = EMESolver(num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
                    oxide_width=SIM_W, oxide_height=SIM_H, res=RES,
                    material_stack=SOI_160, mode_solver="meep",
                    guard_duplicate_modes=guard)
    modes, _ = eme.section_modes(geometry, WL, x_center - half_span,
                                 x_center + half_span, N_MODES, mesh=mesh)
    return modes


def _meow_via_raw_meow(geometry, x_center, half_span, mesh):
    """The audit script's own call path, for comparison."""
    import meow as mw
    x0, x1 = x_center - half_span, x_center + half_span
    gds = geometry.section_gds(x0, x1, SIM_W)
    cells = mw.create_cells(structures=extrude_polygons(gds, SOI_160.extrusion_rules),
                            mesh=mesh, Ls=[x1 - x0])
    cs = mw.CrossSection.from_cell(cell=cells[0], env=mw.Environment(wl=WL, T=25.0))
    return mw.compute_modes(cs, N_MODES)


def main():
    mmi1, _mmi2, mux = _build()
    mesh = _mesh()
    y, z = mesh_cell_centers(mesh)
    dy, dz = float(y[1] - y[0]), float(z[1] - z[0])

    half_total = mux.total_length / 2
    x_dev = (-half_total + TAPER_L) + BODY_STUB          # device frame (the check's x)
    x_loc = -(mmi1.mm_length / 2 - BODY_STUB)            # mmi1 local frame (the audit's x)
    print(f"interface x: device={x_dev:.4f} (on mux), local={x_loc:.4f} (on mmi1)")
    print(f"mesh: {len(mesh.x)}x{len(mesh.y)} nodes, dy={dy:.5f} dz={dz:.5f}", flush=True)

    n_clad = SOI_160.sio2_index(WL)
    print(f"cladding index @ {WL} um = {n_clad:.6f}", flush=True)

    # ---- Item 0: the EME duplicate-mode guard, before vs after ----------
    # The layer-4 control found meow returning modes 10 and 11 as the SAME
    # eigenvector here (overlap 1.0000, row completeness 2.0003) -- the
    # FDTDSolver bug on the other engine. EMESolver.compute_modes() now
    # guards it. This measures both, in one job, on the same cross-section,
    # so the fix is demonstrated rather than assumed.
    print("\n===== ITEM 0: EME duplicate-mode guard, before vs after =====", flush=True)
    for label, guard in (("guard OFF (pre-2026-07-30 behaviour)", False),
                         ("guard ON  (current default)", True)):
        modes = _meow_via_emesolver(mmi1, x_loc, 0.05, mesh, guard=guard)
        f = stack_meow_mode_fields(modes)
        n = np.array([float(np.real(m.neff)) for m in modes])
        d = self_overlap_diagnostics(f, dy, dz)
        _report(f"CONTROL meow(mmi1) vs itself -- {label}", d)
        print(f"  n_eff: {np.round(n, 5).tolist()}", flush=True)
        verdict = ("FAIL: duplicate mode -- basis is RANK-DEFICIENT"
                   if d["has_duplicate"] else
                   "PASS: basis is rank-complete"
                   + ("" if d["is_orthogonal"]
                      else f" (mildly non-orthogonal, max off-diag {d['max_off_diagonal']:.4f}"
                           " -- near-degenerate pairs, benign)"))
        print(f"  >>> {verdict}",
              flush=True)

    # ---- Layer 4 controls first: these MUST be the identity. -------------
    print("\n===== LAYER 4: same-engine controls =====", flush=True)
    modes_mux = _meow_via_emesolver(mux, x_dev, 0.05, mesh)
    f_mux = stack_meow_mode_fields(modes_mux)
    n_mux = np.array([float(np.real(m.neff)) for m in modes_mux])
    _report("CONTROL meow(mux) vs itself", self_overlap_diagnostics(f_mux, dy, dz))

    modes_mmi1 = _meow_via_emesolver(mmi1, x_loc, 0.05, mesh)
    f_mmi1 = stack_meow_mode_fields(modes_mmi1)
    n_mmi1 = np.array([float(np.real(m.neff)) for m in modes_mmi1])
    _report("CONTROL meow(mmi1) vs itself", self_overlap_diagnostics(f_mmi1, dy, dz))

    # ---- A: geometry queried (mux vs mmi1), same engine + entry point ----
    print("\n===== A: geometry -- mux vs mmi1 (both meow/EMESolver) =====", flush=True)
    print("Same physical cross-section, so this MUST be the identity too.", flush=True)
    _report("meow(mux) vs meow(mmi1)",
            overlap_diagnostics(field_overlap_matrix(f_mux, f_mmi1, dy=dy, dz=dz)),
            n_mux, n_mmi1)

    # ---- B: meow entry point (EMESolver.section_modes vs raw mw) ---------
    print("\n===== B: meow entry point -- EMESolver vs raw mw.compute_modes =====", flush=True)
    modes_raw = _meow_via_raw_meow(mmi1, x_loc, 0.05, mesh)
    f_raw = stack_meow_mode_fields(modes_raw)
    n_raw = np.array([float(np.real(m.neff)) for m in modes_raw])
    _report("meow(mmi1, EMESolver) vs meow(mmi1, raw mw)",
            overlap_diagnostics(field_overlap_matrix(f_mmi1, f_raw, dy=dy, dz=dz)),
            n_mmi1, n_raw)

    # ---- C: half-span 0.05 vs 0.1 ---------------------------------------
    print("\n===== C: half-span -- 0.05 vs 0.10 =====", flush=True)
    modes_hs = _meow_via_emesolver(mmi1, x_loc, 0.10, mesh)
    f_hs = stack_meow_mode_fields(modes_hs)
    n_hs = np.array([float(np.real(m.neff)) for m in modes_hs])
    _report("meow(mmi1, hs=0.05) vs meow(mmi1, hs=0.10)",
            overlap_diagnostics(field_overlap_matrix(f_mmi1, f_hs, dy=dy, dz=dz)),
            n_mmi1, n_hs)

    # ---- The real cross-engine comparison, both geometries --------------
    print("\n===== CROSS-ENGINE: meep vs meow =====", flush=True)
    fdtd = FDTDSolver(res=RES, sim_width=SIM_W, sim_height=SIM_H, dpml=0.5,
                      source_margin=0.0, monitor_margin=0.25, three_d=True,
                      material_stack=SOI_160, decay_tol=1e-11)
    res_meep = fdtd.section_mode_fields(mmi1, WL, x_loc - 0.05, x_loc + 0.05,
                                        N_MODES, y, z, monitor_center=0.0,
                                        monitor_size=SIM_W)
    f_meep = {k: res_meep[k] for k in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")}
    n_meep = np.asarray(res_meep["neff"], dtype=float)
    _report("CONTROL meep(mmi1) vs itself", self_overlap_diagnostics(f_meep, dy, dz))

    for label, fb, nb in (("meow(mmi1)", f_mmi1, n_mmi1), ("meow(mux)", f_mux, n_mux)):
        _report(f"meep(mmi1) vs {label}  [ALL {N_MODES} modes]",
                overlap_diagnostics(field_overlap_matrix(f_meep, fb, dy=dy, dz=dz)),
                n_meep, nb)
        ga, gb = n_meep > n_clad, nb > n_clad
        fa_g = {k: v[ga] for k, v in f_meep.items()}
        fb_g = {k: v[gb] for k, v in fb.items()}
        _report(f"meep(mmi1) vs {label}  [GUIDED only: {ga.sum()} vs {gb.sum()}]",
                overlap_diagnostics(field_overlap_matrix(fa_g, fb_g, dy=dy, dz=dz)),
                n_meep[ga], nb[gb])

    # ---- The band-window hypothesis -------------------------------------
    # FDTDSolver._iter_distinct_bands() SKIPS duplicate meep bands, so to
    # return N distinct modes it reaches past band_num N (on this device it
    # skipped 11 and 15, so it went to 18). meow, asked for N, returns its own
    # N LOWEST. The two sets are then genuinely different: meep's top modes
    # can lie above meow's window entirely, with no partner to match -- which
    # shows up as row completeness ~0, NOT as a mis-ordering. If that is the
    # 0.0000, over-requesting on the meow side must repair it.
    print("\n===== BAND WINDOW: over-request on the meow side =====", flush=True)
    n_over = N_MODES + 6
    eme_over = EMESolver(num_modes=[n_over], num_cells=[1], lengths=[1.0],
                         oxide_width=SIM_W, oxide_height=SIM_H, res=RES,
                         material_stack=SOI_160, mode_solver="meep")
    modes_over, _ = eme_over.section_modes(mmi1, WL, x_loc - 0.05, x_loc + 0.05,
                                           n_over, mesh=mesh)
    f_over = stack_meow_mode_fields(modes_over)
    n_over_neff = np.array([float(np.real(m.neff)) for m in modes_over])
    ga = n_meep > n_clad
    gb_over = n_over_neff > n_clad
    _report(f"meep({N_MODES}) vs meow({n_over}), guided: {ga.sum()} vs {gb_over.sum()}",
            overlap_diagnostics(field_overlap_matrix(
                {k: v[ga] for k, v in f_meep.items()},
                {k: v[gb_over] for k, v in f_over.items()}, dy=dy, dz=dz)),
            n_meep[ga], n_over_neff[gb_over])

    # ---- H-only cross-check (to-do item 4) -------------------------------
    # Now promoted into HybridSolver's own check_mode_correspondence (see
    # beef/solvers/hybrid.py); kept here too as a standalone repro.
    print("\n===== H-only overlap, as an independent corroboration =====", flush=True)
    ga, gb = n_meep > n_clad, n_mmi1 > n_clad
    _report("meep vs meow(mmi1), H-only, guided",
            overlap_diagnostics(field_overlap_matrix(
                {k: v[ga] for k, v in f_meep.items()},
                {k: v[gb] for k, v in f_mmi1.items()}, dy=dy, dz=dz,
                components=("Hy", "Hz"))),
            n_meep[ga], n_mmi1[gb])

    print("\nDONE", flush=True)


if __name__ == "__main__":
    main()
