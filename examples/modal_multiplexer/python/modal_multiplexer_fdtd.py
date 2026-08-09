"""
UNBLOCKED 2026-08-04 -- geometry construction now works; this script has
still never been run for real. `FDTDSolver.solve()` (the full single-domain
path this script needs, as opposed to `solve_section()`, the only thing
`HybridSolver` ever calls) used to be unconditionally incompatible with
`CompoundGeometry` as the top-level device: `solve()` -> `to_meep_prisms()`
-> `geometry.polygons(stub_length=self.waveguide_length + self.dpml)`, and
`CompoundGeometry.polygons(stub_length>0)` raised `NotImplementedError` by
design -- "no generic way to know per-port stub width for an arbitrary
segment type" (its own docstring said). That framing turned out to be wrong,
not just unimplemented: `_end_ports()` already answers the question for
`to_gds_eme()`'s own stub handling, by returning EVERY port at an end rather
than picking a single winner, and a leaf geometry with several ports already
draws one stub per port for its own `polygons()` (e.g. `MMI.polygons()`).
Fixed by wiring `_end_ports()` into `polygons()` the same way `to_gds_eme()`
already uses it. See `modal_multiplexer_validation.md` gap 3 and `notes.md`
Geometry item 8 for the full account, evidence, and new tests
(`tests/test_geometry.py::TestCompoundGeometryPolygonsStub`).

Confirmed cleared through geometry/prism construction only: this exact
device's `CompoundGeometry` now builds the correct 4 stub polygons (1 input +
3 output, each at its own port's true width/position) where it used to raise
immediately, and `to_meep_prisms()` -- the actual call `FDTDSolver.solve()`
makes -- now succeeds. **Not yet confirmed**: an actual FDTD solve of this
script, end to end. Before trusting any output, this needs the same
two-step discipline every other expensive device in this repo has used --
a real (not login-node) smoke test at reduced scope, then the full run --
not a login-node try (see the resolution note below for why even the
reduced res=20 default here is not cheap).

One separate, smaller bug found and FIXED while writing this script:
`CompoundGeometry.source_aperture_at()` used to raise for any x past its own
footprint (the geometry-preview plot below hit this immediately), where every
leaf geometry already clamps to its terminal value for arbitrarily large
overshoot. Fixed in `beef/geometry/compound.py` (`_locate(..., clamp=True)`),
tested (`tests/test_geometry.py::TestCompoundGeometrySourceApertureClamps`).
That fix is real and stands on its own; it was necessary but not sufficient
to get this script running.

Full-device FDTD reference for the 2-stage MMI modal multiplexer.
(modal_multiplexer_validation.md gap 3.)

`modal_multiplexer_hybrid.py` is the only script this device has — there is
no ground truth to check it against, unlike every other device in this repo
(both MMIs and the coupler each have a full-device FDTD reference their hybrid
is compared to; mmi_1x2_validation.md / mmi_2x2_validation.md gap 5). This
script is that reference: the SAME `CompoundGeometry` (MMI 1 + connector +
MMI 2), one single-domain meep simulation across the whole 62 µm device
instead of five cascaded sections.

The multi-mode input carries TE0 and TE1 (see modal_multiplexer_hybrid.py's
module docstring for the device physics), so this needs TWO broadband runs —
one per input mode, `input_mode=1` then `input_mode=2` — each producing a
transmission spectrum to all three output branches (top/mid/bot) at once via
three output mode monitors. `FDTDSolver.solve()` excites and normalises
exactly one input mode per call; `num_input_modes`/`num_output_modes` only
size the monitors and are not themselves a multi-excitation switch (checked
against the solve() source before writing this).

Resolution — deliberately REDUCED from the shipped RESOLUTION = 40, per gap
3's own recommendation ("worth scoping at reduced resolution first: even a
res-20 reference that only resolves the qualitative mode-sorting would catch
a cascade that is mixing modes, which is the failure mode gap 1 is about").
This device is also the longest single-domain run this repo would attempt:
62.16 µm against the 2×2 MMI's ~36 µm, and 3-D FDTD cost is steep in
resolution (cell count ~res^3, timesteps ~res, so total cost ~res^4) — a
res-20 run here should still be markedly cheaper than the 2×2's res-40
reference (job 11585844_1, 1:16:45) despite being ~1.7x longer, which is the
whole point of starting reduced rather than at the production 40.
`MUX_FDTD_RES` overrides it, so a later production-resolution re-run needs no
edit here, only an env var, matching diagnose_dup_band_grid.py's
BEEF_DIAG_RES convention.

Output monitors sit at the three branches' y-positions (+WO, 0, -WO) with the
SAME BRANCH_SPECS/WO the hybrid script uses, so "top"/"mid"/"bot" here are the
same physical ports as there. UPDATED 2026-08-07: no longer an explicit
0.9 um width -- that was narrower than it needed to be (a likely contributor
to job 11707932's near-zero result) and has been removed in favor of
FDTDSolver.solve()'s own (now-fixed) default: as wide as each monitor can be
while staying centered on its own branch and inside the domain
(sim_width - 2*abs(center)). See modal_multiplexer_validation.md gap 8.

What this does NOT yet do: no golden (gap 4 — write after this is confirmed
sane and, per every other device's Order section, only once the resolution
question is settled), no comparison script against modal_multiplexer_hybrid.py
(the natural next step, sibling to verify_mmi_1x2_hybrid_vs_fdtd.py /
verify_mmi_2x2_hybrid_vs_fdtd.py — needs this script's SMatrix to exist
first, and both scripts' wavelength grids to match exactly, which they do
here by construction).

Geometry construction no longer raises (see the UNBLOCKED note at the top),
but this script has never actually been run -- no smoke test, no real solve,
no evidence its output is sane, independent of the geometry fix being
correct. Run via the batch script, not a login-node try (3-D FDTD at any
resolution belongs in SLURM, not a background process — same lesson this
device's own hybrid script's OOM history already recorded):
    python /path/to/examples/modal_multiplexer/python/modal_multiplexer_fdtd.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import MMI, MMIPort, CompoundGeometry, WaveguideRun, FDTDSolver, SOI_160, plot_fdtd_geometry

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "modal_multiplexer", "modal_multiplexer_fdtd"))
os.makedirs(output_dir, exist_ok=True)
# No save_dir/example_sim_cache -- caching is disabled for this script. See
# the comment at solver.solve() below (notes.md Geometry item 8).

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = float(os.environ.get("MUX_FDTD_RES", "20.0"))
                    # reduced from modal_multiplexer_hybrid.py's 40 -- see
                    # module docstring. Override: MUX_FDTD_RES=40 for a later
                    # production-resolution re-run.
DECAY_TOL = float(os.environ.get("MUX_FDTD_DECAY_TOL", "1e-9"))
                    # CHANGED default 1e-11->1e-9, 2026-08-08, at the user's
                    # explicit sign-off -- a methodology experiment, not a bug
                    # fix. Two real res=20 runs (11707932, 11734286) both
                    # decayed to catastrophic near-zero transmission
                    # (~1e-6-1e-15) despite two confirmed, unrelated fixes
                    # (TAPER_L, output monitor width) that independently made
                    # the res=10 smoke test sane; TE0 needed t~77,175-77,197
                    # to satisfy decay_tol=1e-11 in BOTH real runs (24-26x
                    # longer than res=10's t~2910-3217), correlating the
                    # failure with resolution/decay-time, not geometry.
                    # Leading, NOT YET CONFIRMED hypothesis
                    # (modal_multiplexer_validation.md gap 3): stop_when_
                    # dft_decayed()'s relative-change criterion may force the
                    # sim to keep running long after the real signal has
                    # already exited, tracking some far smaller, slowly-
                    # decaying spurious/resonant component instead -- by the
                    # time that relative criterion is finally satisfied at
                    # res=20, the genuine output-port signal may already have
                    # decayed into numerical noise, while at res=10 the same
                    # criterion happens to trip before that happens. 1e-9
                    # matches this repo's own precedent as "the known
                    # dominant cost knob" (run_modal_multiplexer_fdtd.batch's
                    # own historical comments) -- loosening it should let the
                    # sim stop while the real signal is still well above the
                    # noise floor, if this hypothesis is right. Override:
                    # MUX_FDTD_DECAY_TOL=1e-11 to reproduce the two prior
                    # (broken) runs exactly.

MIN_RUN_TIME = float(os.environ.get("MUX_FDTD_MIN_RUN_TIME", "1000.0"))
                    # ADDED 2026-08-08, alongside the FDTDSolver.
                    # minimum_run_time fix (beef/solvers/fdtd.py) -- found
                    # BECAUSE testing DECAY_TOL=1e-9 above (job 11749567)
                    # stopped the sim at t=107-112, not even one transit's
                    # worth of time for a 66 um device, before the pulse
                    # could physically have reached the output monitors.
                    # 1000 is a generous floor: Sx (device + stubs + PML)
                    # is ~68 um, and even a conservative silicon group index
                    # ~4.5 only needs ~68*4.5=306 meep time units for one
                    # transit -- 1000 leaves >3x margin for reflections/
                    # multiple bounces before decay-checking is even
                    # evaluated, while staying far below the ~77,175 that
                    # decay_tol=1e-11 (no floor needed there -- tight enough
                    # tolerance already rules out the trivial-early-stop) has
                    # separately shown is needed regardless of this floor.

WL_MIN = 1.26       # O-band, same grid as modal_multiplexer_hybrid.py so a
WL_MAX = 1.36       # future comparison script needs no reconciliation beyond
N_WL = 11           # port naming (both already use top/mid/bot).

# ---------------------------------------------------------------------------
# 1. Geometry — byte-for-byte the same construction as
#    modal_multiplexer_hybrid.py's step 1 (kept in sync by hand; both are
#    short enough that this is the lower-risk option over a shared helper
#    module for two call sites).
# ---------------------------------------------------------------------------
TAPER_L = 3.0                # every taper's length (µm). CORRECTED
                              # 2026-08-07 from 2.0 -- the published seed-
                              # device geometry (Hiesener et al., "Topology
                              # Optimization Enhancements for Multimode
                              # Interferometers," Table 2) reports
                              # L_taper = 3.0 um for this exact device; 2.0
                              # was never measured/derived, just wrong. Every
                              # boundary/domain-size computation below derives
                              # from this constant, so the fix propagates
                              # automatically -- not independently re-verified
                              # against a real run yet.
WO = 1.0375                  # centre-to-centre pitch of the three branches (µm)
BETWEEN_WG_LENGTH = 3.0      # connecting-waveguide length between the MMIs (µm)
MM_WIDTH = 2.8               # both MMI body widths (µm)

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
    use_neighbor_aperture=True,   # irrelevant to a full single-domain FDTD run
                                  # (no per-section reverse launch here), kept
                                  # so this geometry is bit-identical to the
                                  # hybrid script's, not merely equivalent.
)
connector = WaveguideRun(
    length=BETWEEN_WG_LENGTH,
    ports=[(y, ww) for y, ww, _ in BRANCH_SPECS],
)
mux = CompoundGeometry([mmi1, connector, mmi2])
print(f"MMI MUX total length: {mux.total_length:.2f} µm "
      f"(MMI 1 {mmi1.mm_length:.2f} µm + MMI 2 {mmi2.mm_length:.2f} µm, "
      f"width {MM_WIDTH:.2f} µm, branch pitch {WO:.4f} µm)")
print(f"resolution = {RESOLUTION:.0f} px/um (shipped hybrid RESOLUTION = 40 -- "
      f"see module docstring)")
print(f"decay_tol = {DECAY_TOL:.0e} (default changed 1e-11->1e-9, 2026-08-08 -- "
      f"see module docstring/DECAY_TOL comment)")
print(f"minimum_run_time = {MIN_RUN_TIME:.0f} (new 2026-08-08, guards against "
      f"decay_tol being satisfied before the pulse crosses the device -- see "
      f"module docstring/MIN_RUN_TIME comment)")

# top -> mid -> bot, same index convention as modal_multiplexer_hybrid.py's
# PORTS / the legacy master's PORT_DEFS (0=top, 1=mid, 2=bot).
PORT_NAMES = ["top", "mid", "bot"]

# ---------------------------------------------------------------------------
# 2. Configure FDTD solver — single domain across the WHOLE device, unlike
#    modal_multiplexer_hybrid.py's five cascaded sections.
# ---------------------------------------------------------------------------
solver = FDTDSolver(
    res=RESOLUTION,
    sim_width=5.0,          # matches modal_multiplexer_hybrid.py's oxide_width
    sim_height=1.5,         # matches its oxide_height
    dpml=0.5,
    source_margin=0.0,
    monitor_margin=0.25,
    three_d=True,       # full 3-D FDTD: real 160 nm Si core in oxide, PML on all six faces
    num_input_modes=2,       # the input waveguide is genuinely multi-mode (TE0+TE1)
    num_output_modes=1,      # one mode (the branch fundamental) per output monitor
    # REMOVED explicit output_monitor_widths=[0.9, 0.9, 0.9] override, 2026-08-07
    # -- that hand-picked width was narrower than it needed to be, clipping real
    # field content (and the surrounding cladding context get_eigenmode_
    # coefficients()'s own MPB solve needs) before the overlap integral runs;
    # a strong, though not yet independently confirmed, suspect for job
    # 11707932's near-total-loss result (modal_multiplexer_validation.md gap
    # 3). None now uses FDTDSolver's own default: as wide as each monitor can
    # be while staying centered on its own branch and inside the domain
    # (sim_width - 2*abs(center) -- full sim_width for the on-axis "mid"
    # branch, 5.0 - 2*1.0375 = 2.925 um for the off-axis "top"/"bot" branches
    # at WO=1.0375 um pitch), matching the clamp already used for the
    # excited/source-side aperture elsewhere in FDTDSolver. Not re-run yet.
    output_monitor_centers=[+WO, 0.0, -WO],
    material_stack=SOI_160,
    decay_tol=DECAY_TOL,
    minimum_run_time=MIN_RUN_TIME,
)

# ---------------------------------------------------------------------------
# 3. Geometry plot
# ---------------------------------------------------------------------------
plot_fdtd_geometry(
    mux, solver=solver,
    title="2-Stage MMI Modal Multiplexer — Top-down Geometry (full-device FDTD)",
    output_path=os.path.join(output_dir, "geometry.png"),
)
print(f"Geometry plot saved → {output_dir}/geometry.png")

# ---------------------------------------------------------------------------
# 4. Wavelength sweep — one broadband run per input mode (TE0, then TE1);
#    each captures all three output branches' fundamental mode via DFT in a
#    single sim.run(). Merged into one SMatrix keyed by in_mode, the same
#    _data.update() pattern modal_multiplexer_hybrid.py uses to merge its
#    per-port solves.
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nSweeping {len(wavelengths)} wavelengths: "
      f"{wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm, 2 input modes (TE0, TE1)\n")

smat = None
for input_mode in (1, 2):
    print(f"-- input_mode={input_mode} (TE{input_mode - 1}) --")
    smat_mode = solver.solve(mux, wavelengths, input_mode=input_mode)
    if smat is None:
        smat = smat_mode
    else:
        smat._data.update(smat_mode._data)

smat.check_passive()  # sanity: total output power must not exceed 1 at any wavelength

# ---------------------------------------------------------------------------
# 5. Print summary — out_0/out_1/out_2 (solver's own naming, ordered by
#    output_monitor_centers) relabelled top/mid/bot for readability, matching
#    modal_multiplexer_hybrid.py's PORTS convention.
# ---------------------------------------------------------------------------
wl_c = wavelengths[len(wavelengths) // 2]
idx = np.argmin(np.abs(wavelengths - wl_c))

print(f"\nResults at λ = {wl_c:.3f} µm:")
for in_mode in (1, 2):
    total = smat.total_power(in_mode=in_mode)
    for k, name in enumerate(PORT_NAMES):
        T = smat.transmission(port=f"out_{k}", in_mode=in_mode)
        mux_tag = " [MUX out]" if name in ("mid", "bot") else ""
        print(f"  TE{in_mode - 1} in → {name}{mux_tag}:  T = {T[idx]:.4f}  "
              f"IL = {10*np.log10(max(T[idx], 1e-20)):.2f} dB")
    print(f"  TE{in_mode - 1} in → all ports total: {total[idx]:.4f} "
          f"(≤ 1; deficit = radiation/reflection + truncated basis)")

# ---------------------------------------------------------------------------
# 6. Save and plot
# ---------------------------------------------------------------------------
smat.save(os.path.join(output_dir, "smatrix.pkl"))
print(f"\nSMatrix saved → {output_dir}/smatrix.pkl")

_style = {1: dict(color="steelblue"), 2: dict(color="darkorange")}
_ls = {"top": "-", "mid": "--", "bot": ":"}

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
for in_mode in (1, 2):
    for k, name in enumerate(PORT_NAMES):
        T = smat.transmission(port=f"out_{k}", in_mode=in_mode)
        mux_tag = " [MUX out]" if name in ("mid", "bot") else ""
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
    for k, name in enumerate(PORT_NAMES):
        T = smat.transmission(port=f"out_{k}", in_mode=in_mode)
        IL = 10 * np.log10(np.maximum(T, 1e-20))
        mux_tag = " [MUX out]" if name in ("mid", "bot") else ""
        ax.plot(wavelengths, IL, label=f"TE{in_mode - 1} in → {name}{mux_tag}",
                linestyle=_ls[name], **_style[in_mode])
ax.axhline(0.0, color="gray", linestyle=":", lw=0.8, label="Lossless")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Insertion Loss (dB)")
ax.set_title("Insertion Loss vs Wavelength")
ax.legend(fontsize=7)
ax.grid(True)

plt.suptitle(f"2-Stage MMI Modal Multiplexer — Full-Device FDTD (res={RESOLUTION:.0f})",
            fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/results.png")
