"""2x2 MMI — input (FDTD, captured) + body (EME) stitched field map.

The 1x2's Step 3 (stitch_mmi1x2_hybrid_fields.py) ported to the 2x2 tutorial
device. Same mechanism: FDTDSolver.solve_section(capture_fields_dir=...) writes
the real transverse interface field, that field is projected onto the EME body
mode basis, and the body is propagated on THOSE coefficients — so both panels
share one absolute amplitude and the body's phase continues the FDTD section's
by construction, with no hand alignment.

THREE THINGS THAT DIFFER FROM THE 1x2, all load-bearing:

1. THE LAUNCH IS OFF-AXIS. mmi_2x2_hybrid.py excites input_idx=1 (top arm,
   y = +0.75), so `geometry_kwargs` must carry it. The 1x2's odd-mode-fraction
   parity check is therefore meaningless here — an off-axis launch is SUPPOSED
   to drive odd modes — and is replaced by the transfer-matrix cross-check
   below, which is the stronger test anyway.

2. half_domain_aperture=True. The 2x2 has a neighbour port 1.5 um away on both
   sides, so the input section confines its SOURCE aperture to the near half of
   the domain. Verified this does not compromise the capture:
   section_fdtd_layout() applies half_domain_aperture only to src_width/src_c
   (fdtd.py:778) — mea_w is untouched, and the interface capture plane reuses
   the MEASURE monitor's region verbatim, so the projection plane is still the
   full width the body modes live on.

3. decay_tol IS DELIBERATELY 1e-9, NOT the tutorial's 1e-11. 1e-11 is below
   this pipeline's measured hard floor (it never converges, then diverges); the
   1x2 Step 3 run used 1e-9 and that is what is reproduced here. This is the
   one place this script knowingly departs from mmi_2x2_hybrid.py's config, so
   the transfer matrix computed here is NOT bit-comparable to the cascade's.
   The cross-check below compares the capture against THIS run's own T_in, so
   it remains self-consistent.

`reverse` is False here by the ordinary rule (n_modes_in=1 <= n_modes_out=14 is
false), so this is a single forward excite from the physical port — the same
configuration the 1x2 used and the one whose field map means what it looks
like. Contrast capture_mmi1x2_output_section.py, which measured why a
reverse-excited section cannot serve as a forward field panel.

Everything else is pinned to examples/mmi_2x2/python/mmi_2x2_hybrid.py.
MMI2X2_STITCH_RES overrides the resolution for a cheap pipeline check.

Run:
    MMI2X2_STITCH_RES=20 python .../stitch_mmi2x2_hybrid_fields.py
"""
import os

import builtins

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import (MMI, EMESolver, FDTDSolver, SOI_220, solve_body_modes,
                  load_captured_section, project_interface_onto_body,
                  stitch_hybrid_fields, dominant_component)


# RANK-GUARD every print. Under `srun --mpi=pmi2 -n 8` every rank executes this
# script, so an unguarded print emits 8 identical lines — observed in job
# 12736213, whose "[1/4]" banner appeared 24 times.
def print(*a, **k):                      # noqa: A001 - deliberate shadow
    try:
        import meep as _mp
        if not _mp.am_master():
            return
    except Exception:
        pass
    builtins.print(*a, **k)


def _is_master():
    try:
        import meep as _mp
        return bool(_mp.am_master())
    except Exception:
        return True


_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))
RES = float(os.environ.get("MMI2X2_STITCH_RES", 40.0))
out_dir = os.path.join(_ROOT, "example_output", "mmi_2x2",
                       f"mmi_2x2_hybrid_stitched_res{int(RES)}")
cap_dir = os.path.join(out_dir, "captured_input_section")
os.makedirs(cap_dir, exist_ok=True)

# --- pinned to mmi_2x2_hybrid.py -------------------------------------------
WL_MIN, WL_MAX, N_WL = 1.530, 1.565, 11
WL_PLOT = 1.5475          # linspace(1.530, 1.565, 11)[5] — exactly on the grid
N_MODES = 14              # was 15; band 14 sits below the cladding index
BODY_STUB = 0.25
INPUT_IDX = 1             # top arm, the excitation the cascade's S-matrix is for

# Kept as named constants because MMI stores neither on the instance -- the
# widths go into its port objects, and only `width_at(x)` reads them back.
# The figure annotations below need the plain numbers.
WG_WIDTH = 0.45
TAPER_WIDTH = 1.3

mmi = MMI(n_inputs=2, n_outputs=2, mm_length=28.0, mm_width=2.8,
          waveguide_width=WG_WIDTH, taper_length=4.0, taper_width=TAPER_WIDTH,
          waveguide_offset=1.5, waveguide_height=0.22)

half_body, half_total = mmi.mm_length / 2, mmi.total_length / 2
x_in0, x_in1 = -half_total, -half_body + BODY_STUB              # FDTD input
x_b0, x_b1 = -half_body + BODY_STUB, +half_body - BODY_STUB     # EME body
assert abs(x_in1 - x_b0) < 1e-12, "input section's right edge must be the body's left edge"

print(f"2x2 MMI: total {mmi.total_length:.2f} um, body {mmi.mm_length:.2f} um")
print(f"  inputs at y = {mmi.input_ys}, outputs at y = {mmi.output_ys}")
print(f"  FDTD input section x in [{x_in0:+.2f}, {x_in1:+.2f}]  "
      f"({x_in1 - x_in0:.2f} um)")
print(f"  EME body         x in [{x_b0:+.2f}, {x_b1:+.2f}]  "
      f"({x_b1 - x_b0:.2f} um)")

# --- 1. FDTD input section WITH capture ------------------------------------
# RESUME PATH. solve_section() has no cache here (no save_dir), so a rerun
# would repeat the one genuinely expensive step. The capture npz is written
# during the run, so if it and the transfer matrix are already on disk the
# FDTD is skipped entirely and only the (cheap) mode solve, projection and
# plotting rerun. T_in is persisted alongside the capture precisely so the
# cross-check in step 5 survives that skip -- without it a resumed run would
# silently lose its only validation.
cap_npz = os.path.join(cap_dir, "section_fields_band1.npz")
T_path = os.path.join(cap_dir, "T_in.npy")

if os.path.exists(cap_npz) and os.path.exists(T_path):
    T_in = np.load(T_path)
    print(f"\n[1/4] RESUMING from existing capture in {cap_dir}")
    print(f"      transfer matrix {np.asarray(T_in).shape} (loaded, no FDTD run)")
else:
    fdtd = FDTDSolver(res=RES, sim_width=5.0, sim_height=1.5, dpml=0.5,
                      source_margin=0.0, monitor_margin=0.25, three_d=True,
                      material_stack=SOI_220, decay_tol=1e-9,
                      minimum_run_time=100.0, half_domain_aperture=True)
    print(f"\n[1/4] FDTD input section at res {RES:.0f}, input_idx={INPUT_IDX}, "
          f"capturing fields...")
    T_in = fdtd.solve_section(mmi, np.array([WL_PLOT]), x_in0, x_in1,
                              n_modes_in=1, n_modes_out=N_MODES,
                              capture_fields_dir=cap_dir, input_idx=INPUT_IDX)
    print(f"      transfer matrix {np.asarray(T_in).shape}")
    if _is_master():
        np.save(T_path, np.asarray(T_in))

cap = load_captured_section(cap_npz)
print(f"      interface plane at device x = {float(cap['if_pos']):+.4f} "
      f"(body starts at {x_b0:+.4f})")
assert abs(float(cap["if_pos"]) - x_b0) < 0.05, \
    "captured interface plane is not at the body boundary"

# --- 2. Body modes ---------------------------------------------------------
eme = EMESolver(num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
                oxide_width=5.0, oxide_height=1.5, res=RES,
                material_stack=SOI_220, mode_solver="meep")
print(f"[2/4] solving {N_MODES} body modes at {WL_PLOT} um...")
solved = solve_body_modes(mmi, eme, x_b0, x_b1, WL_PLOT, N_MODES,
                          wavelengths=np.linspace(WL_MIN, WL_MAX, N_WL))
print("      n_eff: " + ", ".join(f"{n:.6f}" for n in solved["neff"]))

# --- 3. Project the REAL interface field onto the body basis ---------------
print("[3/4] projecting the captured interface field onto the body modes...")
proj = project_interface_onto_body(cap, solved)
w = proj["weights"]
print("      |c_k|^2 (relative): "
      + ", ".join(f"{v:.4f}" for v in w / max(w.sum(), 1e-300)))

# THE BOUND IS 1.05 HERE, NOT 1.0 + 1e-6 — and that is this device's own
# measured property, not a loosened assertion.
#
# captured_power = sum_k |<phi_k|E>|^2, which Bessel's inequality caps at 1
# only for an ORTHONORMAL basis. This body's 14-mode basis is measurably
# non-orthogonal: its top modes are barely guided (n_eff running down toward a
# ~1.444 cladding), so their long evanescent tails overlap heavily. The 1x2,
# 1x4, DC and HCG-FP bodies all pass 1.0 + 1e-6; this one does not, and
# mmi_2x2_hybrid.py already carries the same 1.05 bound for the same reason
# (commit 9a7fd41, gram_max_offdiag = 0.2837 measured 2026-08-31).
#
# The guard still does its job: the 1-D-cut projection bug it was written for
# returned 1.5857 on the 1x2, which 1.05 rejects ten times over. The gram
# check below tests the DIAGNOSTIC rather than the symptom.
P = solved["fields"][proj["component"]]
P = P / np.sqrt((np.abs(P) ** 2).sum(axis=(1, 2))
                * solved["dy"] * solved["dz"])[:, None, None]
_gram = np.abs(np.einsum("iyz,jyz->ij", P.conj(), P)
               * solved["dy"] * solved["dz"])
GRAM_MAX_OFFDIAG = float(np.abs(_gram - np.diag(np.diag(_gram))).max())
print(f"      basis max off-diagonal overlap {GRAM_MAX_OFFDIAG:.4f} "
      f"(a well-conditioned basis is ~1e-3; this body is known non-orthogonal)")
assert proj["captured_power"] <= 1.05, (
    f"captured_power {proj['captured_power']:.4f} — far above 1. Basis "
    f"non-orthogonality ({GRAM_MAX_OFFDIAG:.4f}) explains a small overshoot, "
    f"not this. Suspect the projection (a 1-D cut instead of the 2-D area "
    f"integral returned 1.5857 on the 1x2).")
assert GRAM_MAX_OFFDIAG < 0.5, (
    f"basis non-orthogonality {GRAM_MAX_OFFDIAG:.4f} — too large to attribute "
    f"to weak confinement; the mode set itself is suspect.")

# NOTE: no parity check here. The 1x2 launches on axis, so its odd body modes
# must be empty and an odd-mode fraction is a real test. This device launches
# from the TOP arm, which is supposed to drive both parities — an odd-mode
# fraction near zero would mean the launch was silently centred, i.e. the
# opposite failure. The cross-check in step 5 is what validates the capture.
odd = w[1::2].sum() / max(w.sum(), 1e-300)
print(f"      odd-mode fraction {odd:.4f} (OFF-AXIS launch: should NOT be ~0)")

# --- 4. Propagate the body on those coefficients and stitch ---------------
print("[4/4] propagating the body and stitching...")
F = solved["fields"][proj["component"]]
iz = int(np.argmax(np.abs(F).sum(axis=(0, 1))))
phi = (F / np.sqrt((np.abs(F) ** 2).sum(axis=(1, 2))
                   * solved["dy"] * solved["dz"])[:, None, None])[:, :, iz]
# 27.5 um of body against a fastest beat period of a few um — sample it far
# denser than the 1x2's 4.0 um span needed.
xs = np.linspace(0.0, x_b1 - x_b0, 1400)
E_body = (np.exp(2j * np.pi * np.outer(xs, solved["neff"]) / WL_PLOT)
          * proj["c"][None, :]) @ phi
body = dict(x=x_b0 + xs, y=solved["y"], I=np.abs(E_body) ** 2,
            component=proj["component"])

if _is_master():
    stitch_hybrid_fields(
        [(cap, "input (FDTD, captured)")], body, out_dir,
        tag=f"mmi2x2_res{int(RES)}_{WL_PLOT:.4f}um",
        port_ys=mmi.output_ys, body_half_width=mmi.mm_width / 2,
        device_label="2×2 MMI")

# Continuity across the seam, measured rather than eyeballed.
comp, _ = dominant_component(cap)
I_fdtd_end = (np.abs(np.asarray(cap[comp])) ** 2)[-1]
I_body_start = body["I"][0]
yf, yb = np.asarray(cap["y"]), solved["y"]
I_f = np.interp(yb, yf, I_fdtd_end)
sf, sb = I_f.sum(), I_body_start.sum()
print(f"\nSEAM CONTINUITY at x = {x_b0:+.3f} um")
print(f"  integrated |E|^2: FDTD side {sf:.4e}  body side {sb:.4e}  "
      f"ratio {sb / sf if sf else float('nan'):.4f}")
print(f"  basis accounted for {proj['captured_power']:.4f} of the interface "
      f"field's power ({N_MODES} modes)")

# --- 5. CROSS-CHECK: field capture vs the transfer matrix ------------------
# T_in[0, :, 0] is the modal amplitude at the section's output boundary for a
# unit input mode, from meep's own get_eigenmode_coefficients() at the flux
# plane. The projection above computes the SAME quantity by an independent
# route: DFT field capture, then an area-integral projection onto EME-solved
# modes. Different mode solver, different extraction. This is the test that
# replaces the 1x2's parity check on this off-axis device.
Tcol = np.asarray(T_in)[0, :, 0]
p_T = np.abs(Tcol) ** 2
p_T = p_T / p_T.sum()
p_c = proj["weights"] / proj["weights"].sum()
print("\nCROSS-CHECK -- capture+projection vs transfer matrix")
print(f"  {'mode':>4} {'|T|^2 norm':>12} {'|c|^2 norm':>12} {'delta':>10}")
for k in range(len(p_T)):
    print(f"  {k:>4} {p_T[k]:12.6f} {p_c[k]:12.6f} {p_c[k]-p_T[k]:+10.2e}")
print(f"  max |delta| {np.abs(p_c - p_T).max():.3e}   "
      f"sum|delta| {np.abs(p_c - p_T).sum():.3e}")

# --- 6. Readable per-section detail figure ---------------------------------
# The combined map compresses 31.75 um into one strip and hides everything.
if _is_master():
    xi, yi = np.asarray(cap["x"]), np.asarray(cap["y"])
    Ii = np.abs(np.asarray(cap[comp])) ** 2
    xb_ax, Ib = body["x"], body["I"]
    wg_half = WG_WIDTH / 2
    wall = mmi.mm_width / 2

    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.9], hspace=0.45,
                          wspace=0.22)

    ax = fig.add_subplot(gs[0, :])
    im = ax.pcolormesh(xi, yi, Ii.T, shading="auto", cmap="inferno")
    for p in mmi.input_ys:
        for e in (p - wg_half, p + wg_half):
            ax.axhline(e, color="cyan", lw=.7, ls="--", alpha=.6)
    for e in (-wall, wall):
        ax.axhline(e, color="lime", lw=.8, alpha=.5)
    ax.set_ylim(-2.0, 2.0)
    ax.set_ylabel("device y (µm)")
    ax.set_title(
        f"(a) INPUT section — FDTD, captured  |{comp}|²   "
        f"x ∈ [{xi.min():.2f}, {xi.max():.2f}] µm,  launched on input_idx="
        f"{INPUT_IDX} (y = {mmi.input_ys[INPUT_IDX]:+.2f} µm)\n"
        f"cyan = {WG_WIDTH} µm access waveguide edges (BOTH arms drawn, only the top "
        f"one excited), green = {mmi.mm_width} µm MMI body wall (taper opens {WG_WIDTH} → {TAPER_WIDTH} µm)",
        fontsize=9)
    plt.colorbar(im, ax=ax, pad=.01)

    ax = fig.add_subplot(gs[1, :])
    im = ax.pcolormesh(xb_ax, yb, Ib.T, shading="auto", cmap="inferno")
    for p in mmi.output_ys:
        ax.axhline(p, color="white", lw=1.1, ls=":", alpha=.95)
    for e in (-wall, wall):
        ax.axhline(e, color="lime", lw=.8, alpha=.5)
    ax.set_ylim(-2.0, 2.0)
    ax.set_ylabel("device y (µm)")
    ax.set_title(
        f"(b) BODY — EME, driven by the captured interface field  |{comp}|²   "
        f"x ∈ [{xb_ax.min():.2f}, {xb_ax.max():.2f}] µm\n"
        f"white dotted = output port centres ±{max(mmi.output_ys):.2f} µm, "
        f"green = body walls ±{wall:.1f} µm", fontsize=9)
    plt.colorbar(im, ax=ax, pad=.01)

    ax = fig.add_subplot(gs[2, 0])
    for xq, lab, c in ((xi[0], f"input entry x={xi[0]:+.2f}", "steelblue"),
                       (xi[-1], f"seam (FDTD side) x={xi[-1]:+.2f}", "darkorange")):
        p = Ii[int(np.argmin(np.abs(xi - xq)))]
        ax.plot(yi, p / p.max(), color=c, label=lab)
    p0 = Ib[0]
    ax.plot(yb, p0 / p0.max(), color="crimson", ls="--",
            label=f"seam (body side) x={xb_ax[0]:+.2f}")
    pe = Ib[-1]
    ax.plot(yb, pe / pe.max(), color="seagreen",
            label=f"body exit x={xb_ax[-1]:+.2f}")
    for p in mmi.output_ys:
        ax.axvline(p, color="gray", lw=.8, ls=":")
    ax.set_xlim(-2.0, 2.0)
    ax.set_xlabel("device y (µm)")
    ax.set_ylabel("|E|² (each normalised)")
    ax.set_title("(c) transverse profiles along the chain\n"
                 "grey dotted = output ports (bar = +0.75, cross = −0.75)",
                 fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(alpha=.3)

    ax = fig.add_subplot(gs[2, 1])
    k = np.arange(N_MODES)
    bw = 0.38
    ax.bar(k - bw / 2, p_T, bw, label="meep eigenmode decomp. |T|²",
           color="lightgray", edgecolor="gray")
    ax.bar(k + bw / 2, p_c, bw, label="DFT capture → EME basis |c|²",
           color="seagreen")
    ax.set_xlabel("body mode index")
    ax.set_ylabel("normalised modal power")
    ax.set_title(
        "(d) two INDEPENDENT routes to the same modal amplitudes\n"
        f"max |Δ| = {np.abs(p_c - p_T).max():.2e},  "
        f"basis captures {proj['captured_power']:.4f} of the interface power",
        fontsize=9)
    ax.set_xticks(k)
    ax.legend(fontsize=8)
    ax.grid(alpha=.3, axis="y")

    fig.suptitle(
        f"2×2 MMI hybrid — input + body fields at res {RES:.0f}, "
        f"λ = {WL_PLOT} µm  (output section NOT included: "
        f"{100 * (x_b1 - x_in0) / mmi.total_length:.0f}% of the device)",
        fontsize=12)
    det = os.path.join(out_dir, "fields", f"detail_input_body_res{int(RES)}.png")
    os.makedirs(os.path.dirname(det), exist_ok=True)
    fig.savefig(det, dpi=160, bbox_inches="tight")
    plt.close("all")
    builtins.print(f"saved -> {det}")

    ypos = np.where(yb > 0, Ib[-1], 0.0)
    yneg = np.where(yb < 0, Ib[-1], 0.0)
    builtins.print(
        f"body exit: top peak |y| = {abs(yb[int(np.argmax(ypos))]):.3f} um, "
        f"bottom peak |y| = {abs(yb[int(np.argmax(yneg))]):.3f} um "
        f"(ports at {max(mmi.output_ys):.3f})")
    builtins.print(
        f"body exit power split: top {ypos.sum():.4e}  bottom {yneg.sum():.4e}  "
        f"ratio {ypos.sum() / max(yneg.sum(), 1e-300):.4f}")

print("\nNOTE: the seam ratio is NOT expected to be 1.0 — the body side is the "
      "truncated-basis reconstruction, so it is missing exactly the power "
      "reported as uncaptured, plus any non-guided content the FDTD slice "
      "still carries at that plane.")
