"""Field reconstruction for a HybridSolver cascade's EME body section.

Why this is separate from ``plot_eme_fields()``, which already plots EME
fields. Two reasons, both measured on the 1x2 MMI:

1. **``plot_eme_fields()`` labels in meow's frame, not the device's.** It hands
   its component string straight to ``mw.propagate_modes()`` with no remap,
   but meow's component labels are projections onto meow's OWN local axes --
   see ``beef.mode_identity.stack_meow_mode_fields()``, which records the
   rotation:

       device Ex (propagation)  = meow Ez
       device Ey (in-plane)     = meow Ex
       device Ez (vertical)     = meow Ey

   So the panel it titles ``|Ex|`` is really the device's IN-PLANE TRANSVERSE
   field, and its ``x [um]`` axis is device y. Fine for a qualitative look; a
   real hazard when the result has to be stitched against FDTD fields, which
   ARE device-frame. This module goes through ``stack_meow_mode_fields()`` and
   labels in the device frame.

2. **It can excite only one mode at a time** (``excite_mode_l``), and a single
   eigenmode in a uniform section is invariant along propagation by
   construction -- constant ``|E|``, phase only. Correct, and it says nothing
   about what the body does. A multimode body's behaviour IS the beating of
   several modes (for an MMI, the self-imaging that sets ``mm_length``), which
   needs a superposition.

**Why the cascade cannot supply this and a fresh solve is required.** With
``calibrate_body_neff=True`` (the default since 2026-08-17)
``EMESolver.solve_section()`` takes its ``neff_override`` fast path, which by
its own comment "skips the eigenmode solve and ``mw.compute_s_matrix()``
outright rather than run and discarded" -- the body transfer matrix is built
analytically as a diagonal ``exp(2j*pi*neff/wl*L)``. So the body's mode fields
are not merely uncached, they are never computed. The flip side is that a
diagonal transfer matrix means NO mode mixing across the body: each mode
propagates independently accumulating only phase, which is exactly what makes
the reconstruction below a few lines of algebra rather than a second cascade.

⚠ **THE EXCITATION IS AN APPROXIMATION.** The field actually entering the body
is whatever the neighbouring FDTD section delivers, and FDTD-side spatial
field capture does not exist yet (zero ``add_dft_fields``/``get_array`` calls
anywhere in ``beef/solvers/fdtd.py``; that is step 2 of the plan above). Every
excitation this module offers is a stated stand-in. The pictures are
physically meaningful -- mode content, beating period, where the self-image
lands -- but **no number from them describes a hybrid interface**, and none
should be quoted as one. Swapping the excitation array for a real captured
interface field is step 3.

Requires a UNIFORM cross-section over ``[x0 - stub_left, x1 + stub_right]``:
the single-cell mode basis below only describes the section if one
cross-section serves both ports, the same precondition
``assert_uniform_cross_section()`` enforces inside
``EMESolver.solve_section()``'s override path. A tapered body needs a real
multi-cell solve and is out of scope here.
"""
import os
from typing import Optional, Sequence

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .._utils import extrude_polygons
from ..mode_identity import stack_meow_mode_fields, mesh_cell_centers


def solve_body_modes(geometry, solver, x0, x1, wl, n_modes,
                     stub_left=0.0, stub_right=0.0,
                     wavelengths=None, geometry_kwargs=None):
    """Solve one hybrid body section's modes and return device-frame fields.

    Reproduces ``EMESolver.solve_section()``'s own cell construction for the
    section ``[x0, x1]`` rather than ``_build_cells()``, which would build
    every whole-device region.

    Args:
        geometry: DeviceGeometry (MMI, DirectionalCoupler, ...).
        solver: the EMESolver the cascade uses for this body section. Its
            ``mode_solver`` matters and is NOT overridden here -- pass the
            device's own (``"meep"`` for most, ``"tidy3d"`` for mmi_1x4, whose
            near-cutoff bands broke ``compute_modes_meep``; see that script).
        x0, x1: section bounds in the device frame.
        wl: scalar wavelength for the field snapshot (µm).
        n_modes: modes to carry. Use the cascade's own interface count --
            ``solve_section()`` reads ``max(n_modes_in, n_modes_out)``, not
            ``solver.num_modes``.
        stub_left, stub_right: extra span folded into the section.
            ``HybridSolver._stubs(k)`` gives a nonzero stub only to the FIRST
            and LAST sections, so an interior body gets 0.0/0.0.
        wavelengths: the FULL production sweep, handed to
            ``_extrusion_rules()``. Passing the device's own array keeps this a
            faithful reproduction of what the cascade builds rather than a
            reasoned-to-be-equivalent shortcut. Defaults to ``[wl]``.
        geometry_kwargs: forwarded to ``geometry.section_gds()`` (e.g.
            ``output_idx``).

    Returns:
        dict with ``modes``, ``cells``, ``neff`` [n_modes], ``fields``
        (device-frame FieldSet), ``y``/``z`` (device y, z cell centres),
        ``dy``/``dz``, and ``length``.
    """
    import meow as mw

    geometry_kwargs = dict(geometry_kwargs or {})
    wavelengths = np.asarray([wl] if wavelengths is None else wavelengths, dtype=float)
    length = (x1 - x0) + stub_left + stub_right

    gds_polys = geometry.section_gds(
        x0, x1, solver.oxide_width,
        stub_left=stub_left, stub_right=stub_right, **geometry_kwargs)
    structs = extrude_polygons(gds_polys, solver._extrusion_rules(wavelengths))
    mesh = solver._build_mesh()
    cells = mw.create_cells(
        structures=structs, mesh=mesh,
        Ls=solver._flat_cell_lengths(region_lengths=[length]))

    env = mw.Environment(wl=float(wl), T=solver.temperature)
    css = [mw.CrossSection.from_cell(cell=c, env=env) for c in cells]
    modes = [solver.compute_modes(cs, n_modes) for cs in css]

    fields = stack_meow_mode_fields(modes[0])
    y, z = mesh_cell_centers(mesh)
    return dict(
        modes=modes[0], cells=cells, mesh=mesh, length=length,
        neff=np.array([complex(m.neff).real for m in modes[0]]),
        fields=fields, y=y, z=z,
        dy=float(y[1] - y[0]), dz=float(z[1] - z[0]),
    )


def body_mode_weights(solved, excitation="gaussian", excitation_center=0.0,
                      excitation_fwhm=None, excitation_mode=0,
                      component="Ey"):
    """Project an approximate input onto the body's mode basis.

    **The projection is a 2-D area integral over the full cross-section, not a
    1-D cut through the core, and that is not a stylistic choice.** Measured on
    the 1x2: a 1-D ``Ey``-slice dot product returned total captured power
    **1.5857** -- impossible for an orthogonal basis -- and put ``|c|^2 =
    0.5876`` on mode 7, an ODD mode a centred input cannot excite at all. Both
    artifacts vanish under the area integral (total 0.9970, mode 7 at 1e-11).
    ``beef.mode_identity.field_overlap_matrix()`` uses the same 2-D transverse
    convention with ``dy``/``dz`` weights, for the same reason. Note the parity
    zeros survived BOTH inner products, so a parity check alone does not catch
    this -- power conservation is what exposes it.

    Args:
        solved: the dict returned by :func:`solve_body_modes`.
        excitation: ``"gaussian"`` for a launched beam of finite width (use the
            feeding taper's width), or ``"mode"`` to excite one body eigenmode.
        excitation_center: device y of the launch centre. 0.0 for a centred
            single input (1x2, 1x4 MMI); the input waveguide's own y for an
            off-axis launch (a directional coupler drives ONE of its two
            guides).
        excitation_fwhm: launch FWHM in µm, required for ``"gaussian"``.
        excitation_mode: 0-based body mode index, for ``"mode"``.
        component: which device-frame component carries the launch. ``"Ey"``
            (in-plane transverse) is the dominant one for the TE-like modes
            every device here uses -- measured ~330x larger than ``Ez`` on the
            1x2 body.

    Returns:
        dict with ``c`` (complex coefficients), ``weights`` (|c|^2), ``phi``
        (normalised 2-D profiles), ``captured_power``, ``gram_max_offdiag``,
        and ``iz`` (the core's vertical slice index).
    """
    F = solved["fields"][component]
    y, z, dy, dz = solved["y"], solved["z"], solved["dy"], solved["dz"]

    # Normalise every mode over the FULL 2-D cross-section.
    P = F / np.sqrt((np.abs(F) ** 2).sum(axis=(1, 2)) * dy * dz)[:, None, None]

    # Report the basis's orthogonality rather than assuming it. Residual
    # off-diagonal is expected (this is an E-only overlap on a discrete grid,
    # not the full ExH power inner product) -- measured 0.069 on the 1x2 body,
    # non-zero only between SAME-parity modes.
    gram = np.abs(np.einsum("iyz,jyz->ij", P.conj(), P) * dy * dz)
    gram_max_offdiag = float(np.abs(gram - np.diag(np.diag(gram))).max())

    # Core's vertical slice: the row of peak modal intensity, rather than
    # assuming the core sits at z=0.
    iz = int(np.argmax(np.abs(F).sum(axis=(0, 1))))

    if excitation == "mode":
        c = np.zeros(len(P), dtype=complex)
        c[excitation_mode] = 1.0
    elif excitation == "gaussian":
        if excitation_fwhm is None:
            raise ValueError("excitation='gaussian' requires excitation_fwhm "
                             "(µm) — use the feeding taper/waveguide width")
        # Separable: a Gaussian across the width times the core's own vertical
        # profile, so the launch lives in the same function space as the modes
        # instead of being a mathematically convenient but unphysical shape.
        iy0 = int(np.argmax(np.abs(P[0]).sum(axis=1)))
        V = P[0][iy0, :]
        V = V / np.sqrt((np.abs(V) ** 2).sum() * dz)
        sigma = excitation_fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        Gy = np.exp(-0.5 * ((y - excitation_center) / sigma) ** 2)
        Gy = Gy / np.sqrt((np.abs(Gy) ** 2).sum() * dy)
        E_exc = Gy[:, None] * V[None, :]
        E_exc = E_exc / np.sqrt((np.abs(E_exc) ** 2).sum() * dy * dz)
        c = np.einsum("kyz,yz->k", P.conj(), E_exc) * dy * dz
    else:
        raise ValueError(f"unknown excitation {excitation!r}, "
                         "expected 'gaussian' or 'mode'")

    return dict(c=c, weights=np.abs(c) ** 2, phi=P, iz=iz,
                captured_power=float((np.abs(c) ** 2).sum()),
                gram_max_offdiag=gram_max_offdiag)


def port_phases(res, port_ys, reference=0, unwrap_to_degrees=True):
    """Relative phase of the reconstructed field at each port position.

    ``|E|^2`` throws away exactly what a phase-sensitive device is built to
    produce -- a 90-degree hybrid's whole function is the inter-port phase
    relationship, and `mmi_4x4_90hybrid`'s own campaign measured pure-EME
    getting that wrong by 44-62 degrees while the intensities still looked
    reasonable. So an intensity map cannot validate such a device at all. This
    reads the phase off the complex field the reconstruction already computed.

    Args:
        res: the dict returned by :func:`plot_hybrid_body_fields` (needs ``E``,
            ``y``).
        port_ys: port centre positions in device y (µm). Take them from the
            geometry (``geometry.output_ys``), not a hand-computed offset.
        reference: index into ``port_ys`` whose phase is subtracted, so the
            result is a RELATIVE phase. An absolute phase here is meaningless:
            it depends on the launch's own arbitrary global phase.
        unwrap_to_degrees: return degrees in (-180, 180] rather than radians.

    Returns:
        dict with ``phase`` [n_ports] relative to the reference, ``amplitude``
        [n_ports] (|E| at each port, so a phase on a near-null port can be
        recognised as meaningless), and ``y_index``.

    ⚠ **Read the amplitudes before trusting any phase.** The phase of a
    vanishing field is numerically arbitrary -- the same trap
    [[mode-signs-gate-on-strong-modes]] records, where a 1e-06 mode disagreed
    by 280%. A port carrying no power has no meaningful phase.

    ⚠ And the standing caveat still applies: these phases come from an
    APPROXIMATE launch, not the cascade's own interface field. They describe
    the body's own phase accumulation, which is real and worth having, but they
    are not the device's port-to-port phases.
    """
    E, y = res["E"], res["y"]
    idx = [int(np.argmin(np.abs(y - p))) for p in port_ys]
    exit_field = E[-1]
    amp = np.abs(exit_field[idx])
    ph = np.angle(exit_field[idx])
    ph = ph - ph[reference]
    ph = (ph + np.pi) % (2 * np.pi) - np.pi     # wrap into (-pi, pi]
    if unwrap_to_degrees:
        ph = np.degrees(ph)
    return dict(phase=ph, amplitude=amp, y_index=idx)


def plot_hybrid_body_fields(
    geometry, solver, x0, x1, wl, n_modes, output_dir, tag,
    *,
    excitation="gaussian", excitation_center=0.0, excitation_fwhm=None,
    excitation_mode=0,
    stub_left=0.0, stub_right=0.0, wavelengths=None, geometry_kwargs=None,
    component="Ey", n_prop=600,
    port_ys: Optional[Sequence[float]] = None,
    body_half_width: Optional[float] = None,
    device_label: str = "",
    total_length: Optional[float] = None,
    solved=None,
    verbose: bool = True,
):
    """Reconstruct and plot a hybrid cascade's EME body section field.

    Solves the body modes, projects an approximate launch onto them, propagates
    the superposition analytically (valid because a calibrated body's transfer
    matrix is diagonal — see the module docstring), and writes one three-panel
    figure in the DEVICE frame.

    See the module docstring for the standing caveat: the excitation is an
    approximation, not the hybrid's own interface field.

    Args:
        port_ys: output port centres to annotate. Read them off the geometry
            (``geometry.output_ys``) rather than hand-copying — for an MMI,
            ``waveguide_offset`` is the port centre-to-centre SPACING, so two
            outputs at ``waveguide_offset=1.0`` sit at ±0.5, not ±1.0.
        body_half_width: sidewall position to annotate (µm).
        solved: a pre-computed :func:`solve_body_modes` result, to plot without
            re-solving.
        Remaining args: see :func:`solve_body_modes` / :func:`body_mode_weights`.

    Returns:
        The merged dict of solve + projection diagnostics, so a caller can
        assert on ``neff``, ``weights``, ``captured_power`` and
        ``gram_max_offdiag`` instead of eyeballing the picture.
    """
    if solved is None:
        solved = solve_body_modes(
            geometry, solver, x0, x1, wl, n_modes,
            stub_left=stub_left, stub_right=stub_right,
            wavelengths=wavelengths, geometry_kwargs=geometry_kwargs)

    proj = body_mode_weights(
        solved, excitation=excitation, excitation_center=excitation_center,
        excitation_fwhm=excitation_fwhm, excitation_mode=excitation_mode,
        component=component)

    y, neff, length = solved["y"], solved["neff"], solved["length"]
    c, w, iz = proj["c"], proj["weights"], proj["iz"]
    phi = proj["phi"][:, :, iz]

    if verbose:
        print(f"body n_eff: " + ", ".join(f"{n:.6f}" for n in neff))
        print(f"basis max off-diagonal overlap: {proj['gram_max_offdiag']:.4f}")
        print("mode weights |c_k|^2: " + ", ".join(f"{v:.4f}" for v in w))
        print(f"captured power: {proj['captured_power']:.4f} of 1.0 "
              f"({n_modes}-mode basis truncation + residual non-orthogonality)")

    # Propagate the superposition. Diagonal transfer matrix => each mode
    # accumulates only its own phase; no mode mixing to account for.
    xs = np.linspace(0.0, length, n_prop)
    phase = np.exp(2j * np.pi * np.outer(xs, neff) / wl)
    E = (phase * c[None, :]) @ phi
    I = np.abs(E) ** 2

    fields_dir = os.path.join(output_dir, "fields")
    os.makedirs(fields_dir, exist_ok=True)

    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], hspace=0.34, wspace=0.24)

    ax = fig.add_subplot(gs[0, :])
    im = ax.pcolormesh(x0 - stub_left + xs, y, I.T, shading="auto", cmap="inferno")
    ann = []
    if body_half_width is not None:
        for e in (-body_half_width, +body_half_width):
            ax.axhline(e, color="cyan", lw=1.0, ls="--", alpha=0.8)
        ann.append(f"cyan = sidewalls (±{body_half_width:.2f} µm)")
    if port_ys is not None:
        for p in port_ys:
            ax.axhline(p, color="white", lw=1.2, ls=":", alpha=0.95)
        ann.append("white dotted = port centres ("
                   + "/".join(f"{p:+.2f}" for p in port_ys) + " µm)")
    ax.set_xlabel("device x — propagation (µm)")
    ax.set_ylabel("device y — in-plane (µm)")
    ax.set_title(f"(a)  |E|² in the hybrid body section, {n_modes}-mode "
                 f"superposition at λ = {wl:.3f} µm"
                 + ("\n" + ",  ".join(ann) if ann else ""), fontsize=10)
    plt.colorbar(im, ax=ax, label="|E|²  (arb.)", pad=0.01)

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(y, I[0] / I.max(), color="steelblue",
            label=f"entry  (x = {x0 - stub_left:+.2f} µm)")
    ax.plot(y, I[-1] / I.max(), color="darkorange", ls="--",
            label=f"exit   (x = {x1 + stub_right:+.2f} µm)")
    if port_ys is not None:
        for p in port_ys:
            ax.axvline(p, color="gray", lw=0.9, ls=":")
    ax.set_xlabel("device y — in-plane (µm)")
    ax.set_ylabel("|E|²  (normalised)")
    ax.set_title("(b)  transverse profile, body entry vs exit", fontsize=10)
    if body_half_width is not None:
        ax.set_xlim(-2 * body_half_width, +2 * body_half_width)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 1])
    kk = np.arange(len(w))
    ax.bar(kk, w, color=["seagreen" if k % 2 == 0 else "lightgray" for k in kk])
    for k, (wk, nk) in enumerate(zip(w, neff)):
        ax.text(k, wk + 0.012 * max(w.max(), 0.1), f"{nk:.3f}",
                ha="center", fontsize=7, rotation=90)
    ax.set_xlabel("body mode index")
    ax.set_ylabel("|c$_k$|²")
    ax.set_title("(c)  modal content of the approximate launch\n"
                 "green = even index, labels = n$_{eff}$", fontsize=10)
    ax.set_xticks(kk)
    ax.set_ylim(0, max(w.max() * 1.30, 0.1))
    ax.grid(alpha=0.3, axis="y")

    frac = (f" — {100.0 * length / total_length:.0f}% of the "
            f"{total_length:.2f} µm device" if total_length else "")
    fig.suptitle(
        f"{device_label or type(geometry).__name__} — hybrid body section field "
        f"reconstruction (device frame){frac}\n"
        "EME body only; neighbouring FDTD sections carry no field data yet. "
        "Excitation is an approximation, not the hybrid's own interface field.",
        fontsize=11)

    out = os.path.join(fields_dir, f"hybrid_body_{tag}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close("all")
    if verbose:
        print(f"saved -> {out}")

    result = dict(solved)
    result.update(proj)
    result["intensity"] = I
    # The COMPLEX field, always returned -- |E|^2 discards the relative phase
    # between transverse positions, which for a 90-degree hybrid is the entire
    # quantity of interest (that device's own campaign measured pure-EME
    # getting inter-port phase wrong by 44-62 degrees). Callers wanting a phase
    # read use this plus port_phases() below rather than re-deriving the sum.
    result["E"] = E
    result["x_prop"] = x0 - stub_left + xs
    result["output_path"] = out
    return result
