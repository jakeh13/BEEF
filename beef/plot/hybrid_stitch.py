"""Stitch a hybrid cascade's FDTD section fields and EME body field into one
continuous device-frame map.

WHAT THIS FIXES ABOUT STEP 1. ``plot_hybrid_body_fields()`` reconstructs the
body from an APPROXIMATE launch (a Gaussian standing in for whatever the
neighbouring FDTD section delivers), because no FDTD field data existed. Step 2
(``FDTDSolver.solve_section(capture_fields_dir=...)``) now provides the real
thing, so the body can be driven by the field that actually arrives.

WHY THAT ALSO SOLVES THE TWO ALIGNMENT PROBLEMS THE PLAN CALLS OUT, rather
than needing separate machinery for them:

* **Amplitude.** The projection coefficients are computed from the captured
  field in its own absolute units, so the body map comes out on the SAME scale
  as the FDTD maps. The panels are one field, not three independently
  normalised pictures.
* **Propagation phase reference.** The coefficients are measured AT the shared
  boundary, and the body propagator ``exp(2j*pi*n_eff*x/wl)`` is referenced to
  that same plane, so the body's phase continues the FDTD section's phase by
  construction. There is no reference to align by hand.

WHAT IS STILL AN APPROXIMATION, stated plainly:

* The body transfer is **diagonal** (each mode accumulates only phase), which
  is exact for a uniform body on the ``neff_override`` path but says nothing
  about a tapered one.
* The captured interface field is projected onto a **truncated** mode basis, so
  whatever falls outside the ``n_modes`` carried is silently dropped. The
  returned ``captured_power`` is the honest measure of that loss, and it is the
  number to quote -- not the picture.
* Reflection is ignored. The cascade's transfer matrices are forward-only here.

Requires the FDTD side to have been run with ``capture_fields_dir`` set, which
writes BOTH a longitudinal x-y slice (for the map) and a transverse plane at
the measure boundary (for the projection). The transverse plane is essential:
projecting a 1-D cut through the core instead is the pathology
``beef.plot.hybrid_body.body_mode_weights`` measured -- captured power 1.5857
and a parity-forbidden mode at 0.5876.
"""
import os
from typing import Optional, Sequence

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_captured_section(path):
    """Load one ``section_fields_band*.npz`` written by Step 2."""
    d = np.load(path)
    out = {k: d[k] for k in d.files}
    for k in ("frame", "if_axis", "if_u_axis", "if_v_axis"):
        if k in out:
            out[k] = str(out[k])
    return out


def dominant_component(cap, prefix=""):
    """Which of Ex/Ey/Ez actually carries the field. NEVER assume.

    In a 2-D run the polarizations decouple and ``Ez`` dominates; in 3-D a
    TE-like mode's in-plane ``Ey`` does. Measured on the 1x2 input section in
    2-D: ``Ez`` 6.99 against ``Ey`` 9.1e-08, i.e. eight orders apart, so a
    hardcoded choice would have analysed pure numerical noise -- and that noise
    still LOOKED like a confined profile, which is why this is a function and
    not a default.
    """
    mags = {c: float(np.abs(cap[prefix + c]).max())
            for c in ("Ex", "Ey", "Ez") if prefix + c in cap}
    if not mags:
        raise KeyError(f"no {prefix}E* components in the capture")
    return max(mags, key=mags.get), mags


def project_interface_onto_body(cap, solved, component=None, verbose=True):
    """Project a captured transverse interface field onto the body mode basis.

    Args:
        cap: dict from :func:`load_captured_section` (needs the ``if_*`` keys).
        solved: dict from ``beef.plot.hybrid_body.solve_body_modes``.
        component: ``"Ex"``/``"Ey"``/``"Ez"``, or None to pick the dominant one.

    Returns:
        dict with ``c`` (complex coefficients, in the captured field's own
        units), ``weights``, ``captured_power`` (fraction of the interface
        field's power the truncated basis accounts for), and ``component``.
    """
    if component is None:
        component, mags = dominant_component(cap, prefix="if_")
        if verbose:
            print("interface component magnitudes: "
                  + ", ".join(f"{k}={v:.3e}" for k, v in mags.items())
                  + f"  -> using {component}")
    F_if = np.atleast_2d(np.asarray(cap["if_" + component]))
    u = np.atleast_1d(cap["if_u"])
    v = np.atleast_1d(cap["if_v"]) if "if_v" in cap else np.array([0.0])

    # Body modes on the EME mesh, device frame.
    P = solved["fields"][component]
    y, z = solved["y"], solved["z"]
    dy, dz = solved["dy"], solved["dz"]

    # The two grids are independent (FDTDSolver's sim_width/sim_height vs
    # EMESolver's oxide_width/oxide_height, and different meshers), so
    # RESAMPLE the captured field onto the mode mesh rather than assuming they
    # coincide. Assuming would silently misalign whenever a device sets those
    # dimensions differently -- which several do.
    if F_if.shape != (len(y), len(z)):
        from scipy.interpolate import RegularGridInterpolator
        if F_if.shape[0] != len(u):
            F_if = F_if.T
        if F_if.shape[1] != len(v):
            v = np.linspace(v.min(), v.max(), F_if.shape[1]) if len(v) > 1 \
                else np.array([0.0])
        if len(v) == 1:
            # Degenerate second axis (2-D run): interpolate in u only and
            # broadcast, so the area integral below stays well-defined.
            re = np.interp(y, u, F_if[:, 0].real, left=0.0, right=0.0)
            im = np.interp(y, u, F_if[:, 0].imag, left=0.0, right=0.0)
            E_if = (re + 1j * im)[:, None] * np.ones((1, len(z)))
        else:
            gi_r = RegularGridInterpolator((u, v), F_if.real,
                                           bounds_error=False, fill_value=0.0)
            gi_i = RegularGridInterpolator((u, v), F_if.imag,
                                           bounds_error=False, fill_value=0.0)
            YY, ZZ = np.meshgrid(y, z, indexing="ij")
            pts = np.stack([YY.ravel(), ZZ.ravel()], axis=-1)
            E_if = (gi_r(pts) + 1j * gi_i(pts)).reshape(len(y), len(z))
    else:
        E_if = F_if

    # 2-D area-integral projection with dy*dz weights -- the same convention
    # body_mode_weights() uses, and for the reason documented there.
    Pn = P / np.sqrt((np.abs(P) ** 2).sum(axis=(1, 2)) * dy * dz)[:, None, None]
    c = np.einsum("kyz,yz->k", Pn.conj(), E_if) * dy * dz

    tot = (np.abs(E_if) ** 2).sum() * dy * dz
    captured = float((np.abs(c) ** 2).sum() / tot) if tot > 0 else float("nan")
    if verbose:
        print(f"projected onto {len(c)} body modes; basis accounts for "
              f"{captured:.4f} of the interface field's power")
    return dict(c=c, weights=np.abs(c) ** 2, captured_power=captured,
                component=component, E_interface=E_if)


def geometry_outline(geometry, x, arm="input"):
    """Core edge curves (upper, lower) of ``geometry`` sampled at ``x``.

    One curve pair per port on the ``arm`` side (``"input"`` or
    ``"output"``) -- each computed via the geometry's own ``width_at``/
    ``center_at`` (so the outline is exactly what was drawn, not a
    reconstruction), evaluated across the FULL ``x`` range: inside the body
    every port's curve coincides (``width_at``/``center_at`` there ignore
    the port index), and each diverges to its own taper once ``x`` leaves
    the body -- so plotting one curve per port draws every arm without
    needing to split ``x`` into regions by hand.
    """
    x = np.asarray(x, dtype=float)
    n_ports = geometry.n_inputs if arm == "input" else geometry.n_outputs
    curves = []
    for k in range(n_ports):
        kw = dict(input_idx=k, output_idx=k)
        w = np.array([geometry.width_at(xi, **kw) for xi in x])
        c = np.array([geometry.center_at(xi, **kw) for xi in x])
        curves.append((c + w / 2, c - w / 2))
    return curves


def draw_geometry_outline(ax, geometry, x, color="deepskyblue", lw=0.9):
    """Overlay ``geometry``'s core edges on ``ax`` across the full ``x`` span."""
    for arm in ("input", "output"):
        for upper, lower in geometry_outline(geometry, x, arm=arm):
            ax.plot(x, upper, color=color, lw=lw, alpha=.8)
            ax.plot(x, lower, color=color, lw=lw, alpha=.8)


def stitch_hybrid_fields(sections, body, output_dir, tag, *,
                         port_ys: Optional[Sequence[float]] = None,
                         body_half_width: Optional[float] = None,
                         device_label: str = "", verbose: bool = True,
                         geometry=None):
    """Assemble one continuous |E|^2 map from FDTD sections + an EME body.

    Args:
        sections: list of ``(cap, role)`` where ``cap`` is a
            :func:`load_captured_section` dict and ``role`` a label
            (e.g. ``"input"``, ``"output"``).
        body: dict with ``x`` (device-frame propagation coords), ``y``, and
            ``I`` (|E|^2 on that grid) for the reconstructed body -- as
            produced by driving the body with
            :func:`project_interface_onto_body`'s coefficients.
        port_ys, body_half_width, device_label: annotation, as in
            ``plot_hybrid_body_fields``.
        geometry: optional ``DeviceGeometry`` (e.g. the ``MMI`` instance) --
            when given, its core outline is overlaid via
            :func:`draw_geometry_outline`, using the panels' own combined
            x-span so the outline lines up with the device frame exactly.

    Returns:
        dict with the combined ``x``/``y``/``I`` and the output path.
    """
    panels = []
    for cap, role in sections:
        comp, _ = dominant_component(cap)
        panels.append(dict(x=np.asarray(cap["x"]), y=np.asarray(cap["y"]),
                           I=np.abs(np.asarray(cap[comp])) ** 2,
                           role=role, comp=comp))
    panels.append(dict(x=np.asarray(body["x"]), y=np.asarray(body["y"]),
                       I=np.asarray(body["I"]), role="body (EME)",
                       comp=body.get("component", "?")))
    panels.sort(key=lambda p: p["x"].min())

    fields_dir = os.path.join(output_dir, "fields")
    os.makedirs(fields_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(14, 4.6))
    vmax = max(float(p["I"].max()) for p in panels)
    for p in panels:
        ax.pcolormesh(p["x"], p["y"], p["I"].T, shading="auto", cmap="inferno",
                      vmin=0.0, vmax=vmax)
        # Section seams, drawn so a discontinuity is visible rather than hidden
        # by a shared colour scale.
        ax.axvline(p["x"].min(), color="deepskyblue", lw=0.8, ls=":", alpha=.8)
        ax.axvline(p["x"].max(), color="deepskyblue", lw=0.8, ls=":", alpha=.8)
        ax.text(0.5 * (p["x"].min() + p["x"].max()), p["y"].max() * 0.92,
                p["role"], color="w", ha="center", fontsize=8)
    if body_half_width is not None:
        for e in (-body_half_width, +body_half_width):
            ax.axhline(e, color="cyan", lw=0.9, ls="--", alpha=.7)
    if port_ys is not None:
        for py in port_ys:
            ax.axhline(py, color="white", lw=1.0, ls=":", alpha=.9)
    if geometry is not None:
        x_full = np.linspace(min(p["x"].min() for p in panels),
                             max(p["x"].max() for p in panels), 2000)
        draw_geometry_outline(ax, geometry, x_full, color="lime")
    ax.set_xlabel("device x — propagation (µm)")
    ax.set_ylabel("device y — in-plane (µm)")
    ax.set_title(
        f"{device_label} — stitched hybrid field map (device frame)\n"
        "ONE shared colour scale and one absolute amplitude: the body is "
        "driven by the FDTD section's own captured interface field, "
        "not an assumed launch.\n"
        "dotted cyan = section seams" + (", solid lime = core outline" if
                                          geometry is not None else ""),
        fontsize=9)
    out = os.path.join(fields_dir, f"hybrid_stitched_{tag}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close("all")
    if verbose:
        print(f"stitched map -> {out}")
    return dict(panels=panels, output_path=out, vmax=vmax)
