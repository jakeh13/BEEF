"""
Tutorial: modal (mode-division) multiplexer — 2-stage MMI, 5-section
FDTD/EME hybrid simulation (FECES method).

Device (160 nm SOI, O-band) — the legacy library's cascaded MMI MUX
(FECES_Legacy/masters/mmi_mux_master.py, configs_mmi/config_mux.json): a 1×3
MMI followed by a 3×3 MMI, joined by three short single-mode waveguides. The
geometry is built by composing two independent MMI devices and a
WaveguideRun connector into one CompoundGeometry, rather than a single
hand-rolled MUX class.

                                                        ____ "top"
    "in" ____                                          /____
    (MM wg,  >----[ MMI 1 (1×3) ]--=conn wgs=--[ MMI 2 (3×3) ]---- "mid" [MUX out]
     TE0+TE1)                                          \\____
                                                             "bot" [MUX out]

    |taper|<- mm_1_length ->|taper|<-conn->|taper|<- mm_2_length ->|taper|

The multi-mode input waveguide carries two modes (TE0, TE1). The first MMI
splits whichever mode arrives into the three connecting arms; because TE0
and TE1 illuminate the 1×3 splitter with different symmetry, they exit it
with *different* relative phases. The second MMI recombines the three arms,
and those mode-dependent phase patterns steer each input mode to a
different output port — that's the (de)multiplexing action this script
characterizes. By reciprocity the same S-parameters describe mux operation
(driving the single-mode ports to synthesize TE0/TE1 on the shared bus).

Hybrid decomposition — this is HybridSolver's generic N-section cascade
(any solver mix, any count >= 2) put to real use with FIVE sections:

    [input wg + taper] [MMI 1 body] [out tapers|conn wgs|in tapers] [MMI 2 body] [out taper + wg]
          FDTD              EME                  FDTD                    EME           FDTD

  - The two MMI bodies are long, uniform-cross-section multimode slabs —
    exactly what EME propagates essentially for free (one mode solve per
    wavelength each; length enters only through analytic phase factors).
  - Everything that *changes shape* along x — the input taper, the whole
    inter-stage transition (three output tapers, three connecting
    waveguides, three input tapers), and the output tapers — is
    non-adiabatic and is solved with FDTD. CompoundGeometry transparently
    splits this section across MMI 1's output-taper band, the connector,
    and MMI 2's input-taper band, delegating each piece to the owning
    sub-geometry.
  - Interfaces between sections sit BODY_STUB inside the MMI bodies, past
    the taper/body cross-section discontinuity, and carry N_MODES body
    supermodes each.
  - The output section always draws the full geometry — all three output
    branches, identical for every port, so the simulated structure never
    differs from the real device. Per-port isolation happens in the launch
    window instead: the reverse eigenmode source and its co-located norm
    monitor shrink to the selected port's aperture (MMI 2's
    use_neighbor_aperture / source_aperture_at) so the mode solve neither
    sees a neighboring branch ~1 µm away nor reaches into the PML. The
    section is still re-simulated once per output port (the launch moves)
    — every other section's transfer matrix is shared across ports via
    HybridSolver's section cache.

FDTD run count: 2 (input, one per MM-input mode) + N_MODES (transition,
one reverse launch per body supermode) + 3 (output, one reverse launch per
port) — plus one EME mode-solve per wavelength per MMI body.

Run:
    python /path/to/examples/modal_multiplexer/modal_multiplexer_hybrid.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from feces import (
    MMI, MMIPort, CompoundGeometry, WaveguideRun,
    EMESolver, FDTDSolver, HybridSolver, SOI_160, plot_hybrid_geometry,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "example_output", "modal_multiplexer", "modal_multiplexer_hybrid"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by both the EME and FDTD sub-solvers
WL_MIN = 1.26       # wavelength sweep start (µm) — O-band, matching the legacy device
WL_MAX = 1.36       # wavelength sweep end (µm)
N_WL = 11           # number of wavelength points

# ---------------------------------------------------------------------------
# 1. Geometry — the legacy cascaded MUX, dimensions straight from
#    FECES_Legacy/configs_mmi/config_mux.json (160 nm SOI), rebuilt as two
#    MMI devices + a WaveguideRun connector instead of one hand-rolled class.
# ---------------------------------------------------------------------------
TAPER_L = 2.0                # every taper's length (µm)
WO = 1.0375                  # centre-to-centre pitch of the three branches (µm)
BETWEEN_WG_LENGTH = 3.0      # connecting-waveguide length between the MMIs (µm)
MM_WIDTH = 2.8               # both MMI body widths (µm)

# (y, waveguide_width, taper_width) for the three inter-stage/output
# branches, top -> mid -> bot — same ordering/values as
# FECES_Legacy/configs_mmi/config_mux.json.
BRANCH_SPECS = [
    (+WO, 0.33, 0.725),   # top
    (0.0, 0.35, 1.15),    # mid
    (-WO, 0.35, 0.725),   # bot
]

def _branch_ports():
    return [MMIPort(y=y, waveguide_width=ww, taper_length=TAPER_L, taper_width=tw)
            for y, ww, tw in BRANCH_SPECS]

mmi1 = MMI(  # 1x3: multi-mode input -> three branches
    n_inputs=1, n_outputs=3,
    mm_length=16.92, mm_width=MM_WIDTH,
    inputs=[MMIPort(y=0.0, waveguide_width=0.7, taper_length=TAPER_L, taper_width=1.4)],
    outputs=_branch_ports(),
    waveguide_height=0.16,
)
mmi2 = MMI(  # 3x3: three branches -> three branches (mux output stage)
    n_inputs=3, n_outputs=3,
    mm_length=34.24, mm_width=MM_WIDTH,
    inputs=_branch_ports(),
    outputs=_branch_ports(),
    waveguide_height=0.16,
    # Three branches only WO (~1 µm) apart, one exactly on the y=0 axis —
    # half_domain_aperture's fixed half-domain rule can't isolate the mid
    # port, so the output section's per-port launch window instead comes
    # from MMI's own neighbor-aware source_aperture_at.
    use_neighbor_aperture=True,
)
connector = WaveguideRun(
    length=BETWEEN_WG_LENGTH,
    ports=[(y, ww) for y, ww, _ in BRANCH_SPECS],
)
mux = CompoundGeometry([mmi1, connector, mmi2])
print(f"MMI MUX total length: {mux.total_length:.2f} µm "
      f"(MMI 1 {mmi1.mm_length:.2f} µm + MMI 2 {mmi2.mm_length:.2f} µm, "
      f"width {MM_WIDTH:.2f} µm, branch pitch {WO:.4f} µm)")

# Output-port table — same index convention as the legacy master
# (PORT_DEFS): 0=top (y=+offset), 1=mid (y=0), 2=bot (y=-offset). The legacy
# master tags mid/bot as the "[MUX out]" pair.
PORTS = [("top", 0), ("mid", 1), ("bot", 2)]

# ---------------------------------------------------------------------------
# 2. Configure sub-solvers. Both FDTDSolver and EMESolver implement the same
#    solve_section() interface, so a single instance of each can fill every
#    slot of its type in the chain — sections are told apart by their
#    x-ranges (and cached per-section), not by which object solves them.
# ---------------------------------------------------------------------------
N_MODES = 16    # body supermodes carried at all four section interfaces. The
                # legacy hybrid used 20–30; each extra mode costs one more
                # FDTD run in the transition section, so this is the main
                # accuracy/cost knob. Don't starve it: the on-axis input is
                # symmetric about y=0, so TE0 excites only *symmetric* body
                # supermodes and TE1 only *antisymmetric* ones (e.g. at the
                # MMI 1 interface TE0 lands in supermodes 1+5, TE1 in 2+6) —
                # the basis must keep several supermodes of EACH parity or
                # one input mode's channel silently truncates to ~zero.
BODY_STUB = 0.75  # µm of MMI body folded into each FDTD section past the
                  # taper/body discontinuity (legacy mm_monitor_position);
                  # must exceed the FDTD solver's monitor_margin.

fdtd_solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=5.0,        # legacy oxide_width
    sim_height=1.5,       # legacy oxide_height
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,       # full 3-D FDTD: real 160 nm Si core in oxide, PML on all six faces
    # No half_domain_aperture: the input side is a lone on-axis waveguide,
    # the interface monitors span the full width anyway (N_MODES > 1), and
    # the output sections' per-port source/norm-monitor windows come from
    # the geometry itself (MMI 2's use_neighbor_aperture, see module
    # docstring), which isolates the selected branch where the fixed
    # half-domain rule couldn't (it can't isolate the on-axis mid port at
    # all).
    material_stack=SOI_160,
    decay_tol=1e-11,
)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholders; solve_section() takes its span from the boundaries
    oxide_width=5.0, oxide_height=1.5, res=RESOLUTION, material_stack=SOI_160,
    mode_solver="meep",  # matches FDTDSolver's own eigenmode basis in this hybrid chain
)

# ---------------------------------------------------------------------------
# 3. Build HybridSolver — 5-section chain:
#    FDTD input | EME MMI 1 | FDTD transition | EME MMI 2 | FDTD output.
# ---------------------------------------------------------------------------
half_total = mux.total_length / 2
body1_start = -half_total + TAPER_L
body1_end = body1_start + mmi1.mm_length
body2_start = body1_end + TAPER_L + BETWEEN_WG_LENGTH + TAPER_L
body2_end = body2_start + mmi2.mm_length

boundaries = [
    -half_total,               # physical input port (input taper start)
    body1_start + BODY_STUB,   # FDTD input      <-> EME MMI 1
    body1_end - BODY_STUB,     # EME MMI 1       <-> FDTD transition
    body2_start + BODY_STUB,   # FDTD transition <-> EME MMI 2
    body2_end - BODY_STUB,     # EME MMI 2       <-> FDTD output
    +half_total,               # physical output port (output taper tips)
]
interface_modes = [N_MODES] * 4

hybrid = HybridSolver(
    solvers=[fdtd_solver, body_solver, fdtd_solver, body_solver, fdtd_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
)

# ---------------------------------------------------------------------------
# 4. Geometry plot — shows which solver covers which section. Sections 0–3
#    are identical for every output port, so they render once; the output
#    section (index 4) renders once per port — same full geometry (all
#    three branches) every time, with the reverse source/norm-monitor
#    markers shrunk to that port's launch aperture — exactly the three
#    simulations run in step 5.
# ---------------------------------------------------------------------------
plot_hybrid_geometry(
    mux, solver=hybrid,
    in_modes=2, out_modes=1,
    title="2-Stage MMI Modal Multiplexer — Hybrid Solver Sections",
    output_path=os.path.join(output_dir, "geometry_hybrid.png"),
    sub_region_dir=os.path.join(output_dir, "sections"),
    sub_region_show_pml=True,
    sub_region_show_eme_detail=True,
    sub_region_configs=[
        (f"port_{name}", {4: {"output_idx": k}}) for name, k in PORTS
    ],
)
print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")
print(f"Sub-region plots saved → {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run the hybrid simulation — once per output port. output_idx is passed
#    ONLY to the output section's kwargs, and a shared cache is passed to
#    every solve() call, so sections 0–3 (input FDTD, both EME bodies, the
#    transition FDTD) are each simulated exactly once — not once per port.
#    in_modes=2 excites TE0 and TE1 of the multi-mode input waveguide
#    (legacy input_modes=[1, 2]); out_modes=1 measures TE0 at each port
#    (legacy output_modes=[1]).
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nHybrid sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, {N_MODES} interface modes)\n")

section_cache = {}
smat = None
for name, k in PORTS:
    smat_port = hybrid.solve(
        mux, wavelengths, in_modes=2, out_modes=1, port=name,
        cache=section_cache, save_dir=save_dir,
        geometry_kwargs={4: {"output_idx": k}},
    )
    if smat is None:
        smat = smat_port
    else:
        smat._data.update(smat_port._data)

smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength

wl_c = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_c))

print(f"\nResults at λ = {wl_c:.3f} µm:")
for in_mode in (1, 2):
    total = smat.total_power(in_mode=in_mode)
    for name, k in PORTS:
        T = smat.transmission(port=name, in_mode=in_mode)
        mux_tag = " [MUX out]" if k in (1, 2) else ""
        print(f"  TE{in_mode - 1} in → {name}{mux_tag}:  T = {T[idx]:.4f}  "
              f"IL = {10*np.log10(max(T[idx], 1e-20)):.2f} dB")
    print(f"  TE{in_mode - 1} in → all ports total: {total[idx]:.4f} "
          f"(≤ 1; deficit = radiation/reflection + truncated basis)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix_hybrid.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix_hybrid.pkl")

_style = {  # color by input mode, linestyle by output port
    1: dict(color="steelblue"),
    2: dict(color="darkorange"),
}
_ls = {"top": "-", "mid": "--", "bot": ":"}

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
for in_mode in (1, 2):
    for name, k in PORTS:
        T = smat.transmission(port=name, in_mode=in_mode)
        mux_tag = " [MUX out]" if k in (1, 2) else ""
        ax.plot(wavelengths, T, label=f"TE{in_mode - 1} in → {name}{mux_tag}",
                linestyle=_ls[name], **_style[in_mode])
ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Transmission")
ax.set_title("Transmission vs Wavelength")
ax.legend(fontsize=7)
ax.grid(True)
ax.set_ylim(0, 1.05)

ax = axes[1]
for in_mode in (1, 2):
    for name, k in PORTS:
        IL = smat.insertion_loss_db(port=name, in_mode=in_mode)
        mux_tag = " [MUX out]" if k in (1, 2) else ""
        ax.plot(wavelengths, IL, label=f"TE{in_mode - 1} in → {name}{mux_tag}",
                linestyle=_ls[name], **_style[in_mode])
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title("Insertion Loss vs Wavelength")
ax.legend(fontsize=7)
ax.grid(True)

plt.suptitle("2-Stage MMI Modal Multiplexer — FDTD/EME/FDTD/EME/FDTD Hybrid", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/hybrid_results.png")
