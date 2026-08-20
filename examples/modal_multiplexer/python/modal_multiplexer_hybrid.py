"""
Tutorial: modal (mode-division) multiplexer — 2-stage MMI, 7-section
FDTD/EME hybrid simulation (BEEF method).

Device (160 nm SOI, O-band) — the legacy library's cascaded MMI MUX
(BEEF_Legacy/masters/mmi_mux_master.py, configs_mmi/config_mux.json): a 1×3
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
(any solver mix, any count >= 2) put to real use with SEVEN sections:

    [in wg+taper] [MMI 1 body] [out tapers] [conn wgs] [in tapers] [MMI 2 body] [out taper+wg]
        FDTD            EME         FDTD         EME        FDTD        EME          FDTD

  - The two MMI bodies are long, uniform-cross-section multimode slabs —
    exactly what EME propagates essentially for free (one mode solve per
    wavelength each; length enters only through analytic phase factors).
  - Everything that *changes shape* along x — the input taper, MMI 1's three
    output tapers, MMI 2's three input tapers, and the output tapers — is
    non-adiabatic and is solved with FDTD.
  - **Split 2026-08-06, was one monolithic "transition" FDTD section
    covering both taper bands AND the connector in between.** The three
    connecting waveguides between the tapers are themselves uniform and
    straight (`WaveguideRun`, explicitly "safe to cut anywhere" per
    `beef.geometry.compound`'s own docstring) — exactly what EME propagates
    for free, same as the MMI bodies. Carving that middle stretch out into
    its own EME section turns one N_MODES-launch FDTD block into two
    CONN_MODES-launch ones either side of it, where CONN_MODES (a handful —
    see its own comment below) is far smaller than N_MODES: a lossless 1×3
    junction can only route power into as many independent channels as the
    three arms support, regardless of how many supermodes (N_MODES) the
    adjoining MMI body basis carries. `FDTDSolver.solve_section()` already
    launches from whichever side of a section has fewer modes (see its own
    docstring) — this is the same trick sections 0/6 already use (2 launches
    from the 2-mode input side, 3 from the 3-port output side by
    reciprocity), just newly available here because there was previously no
    narrow side to launch from. CompoundGeometry needed no changes — cutting
    inside a `WaveguideRun` was already a supported, if previously unused,
    section boundary; only this script's own `boundaries`/`interface_modes`
    changed. NOT done this pass: `mirror_pairs` (the two new taper-FDTD
    sections are geometric mirror images of each other about the
    connector's midpoint — verified by hand, see that section's own
    comment below) could in principle halve the new cost again by solving
    only one and transpose-reusing it for the other; left as a documented,
    unimplemented follow-up rather than compounding two unverified changes
    in one pass.
  - Interfaces between sections sit BODY_STUB inside the MMI bodies (past
    the taper/body cross-section discontinuity, carrying N_MODES body
    supermodes) or CONN_STUB inside the connector (past the taper/connector
    kink, carrying CONN_MODES arm supermodes) — see each constant's own
    comment below for why the two clearances differ.
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

FDTD run count: 2 (input, one per MM-input mode) + 2 x min(N_MODES,
CONN_MODES) (the two taper sections, each launched from the narrower
connector-side basis — down from one N_MODES-launch block) + 3 (output, one
reverse launch per port) — plus one EME mode-solve per wavelength per MMI
body and per wavelength for the connector.

**Reshaped 2026-07-31 after job 11590663 TIMED OUT at 5h for the SECOND
time** (11585847 timed out first, at full 5h/4 nodes; resubmitted unchanged
as 11590663, which reached section 4 -- all 3 per-port geometry plots exist
-- but died mid-way through a per-port eigenmode-overlap monitor sweep,
still not writing final output). Sections 0-3 are shared across ports and
were already fully cached on disk (`save_dir`) by both attempts; the
remaining, apparently-underestimated cost is section 4's per-port solve,
run three times in the same process by the loop below. A third identical
resubmit risks the same failure. Two changes:

1. Each output port can now be run as its own SLURM array task
   (`SLURM_ARRAY_TASK_ID` in {0,1,2} indexes `PORTS`). A task calls
   `hybrid.solve()` for just that port (sections 0-3 still come from the
   shared on-disk cache -- near-instant here since they're already hit --
   section 4 is simulated fresh for that port and its result persists to
   the same `save_dir`) and exits without assembling the combined SMatrix.
   Caveat if this is ever resubmitted from a COLD cache (sections 0-3 not
   yet on disk): launching all 3 array tasks at once would have each of
   them redundantly pay for sections 0-3 rather than sharing the cost the
   way the single-process loop does -- not a concern for this specific
   resubmission, since that cache already exists, but worth knowing before
   reusing this pattern on a device without a warm section 0-3 cache.
2. Running with no `SLURM_ARRAY_TASK_ID` set (the merge/serial step) now
   loops over all 3 ports with `cache_only=True`: it reads every port's
   section-4 result back out of `save_dir` (raising loudly if an array task
   hasn't completed yet) instead of re-simulating, assembles the combined
   SMatrix exactly as before, and continues into `check_passive()`/
   printing/plotting/saving unchanged.

Run (single process, original behavior, still works but not recommended
given the timing above):
    python /path/to/examples/modal_multiplexer/python/modal_multiplexer_hybrid.py

Run (reshaped, recommended): submit the array job for the 3 ports, wait for
all 3 to complete, then submit the merge step (same script, no array index)
to assemble/check/save/plot from cache. See
batch/run_modal_multiplexer_examples.batch and
batch/run_modal_multiplexer_examples_merge.batch.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beef import (
    MMI, MMIPort, CompoundGeometry, WaveguideRun,
    EMESolver, FDTDSolver, HybridSolver, SOI_160, plot_hybrid_geometry,
)
from beef.cache import is_master

_HERE = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_output", "modal_multiplexer", "modal_multiplexer_hybrid"))
os.makedirs(output_dir, exist_ok=True)
save_dir = os.path.abspath(os.path.join(_HERE, "..", "..", "..", "..", "example_sim_cache"))

# ---------------------------------------------------------------------------
# Global simulation parameters
# ---------------------------------------------------------------------------
RESOLUTION = 40.0   # mesh resolution (pixels/µm) — used by both the EME and FDTD sub-solvers.
                    # Kept at the shipped value 2026-08-04 despite gap 2 finding no
                    # resolution 20-100 fully converged in absolute n_eff -- see N_MODES'
                    # own comment below for the full reasoning (this device's check-1
                    # engine ratio, not either engine's absolute n_eff, is what this
                    # chain's accuracy actually depends on).
WL_MIN = 1.26       # wavelength sweep start (µm) — O-band, matching the legacy device
WL_MAX = 1.36       # wavelength sweep end (µm)
N_WL = int(os.environ.get("MUX_HYBRID_N_WL", "11"))
                    # env-overridable, ADDED 2026-08-08 for smoke-testing the
                    # 7-section chain cheaply: FDTD launch COUNT does not
                    # depend on N_WL (each broadband FDTD run captures every
                    # wavelength via DFT in one shot), but every EME mode
                    # solve (both bodies, and the new connector section) is
                    # SOLVED PER WAVELENGTH POINT and does scale directly
                    # with N_WL -- the connector solve alone measured ~84
                    # min/wavelength on 1 rank (gap 7), so this is the
                    # correct lever for a first cheap smoke run of the new
                    # split, not RESOLUTION (which N_MODES/CONN_MODES were
                    # specifically measured against at RESOLUTION=40 -- see
                    # each constant's own comment -- so changing it would
                    # need fresh justification, unlike N_WL).

# ---------------------------------------------------------------------------
# 1. Geometry — the legacy cascaded MUX, dimensions straight from
#    BEEF_Legacy/configs_mmi/config_mux.json (160 nm SOI), rebuilt as two
#    MMI devices + a WaveguideRun connector instead of one hand-rolled class.
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

# (y, waveguide_width, taper_width) for the three inter-stage/output
# branches, top -> mid -> bot — same ordering/values as
# BEEF_Legacy/configs_mmi/config_mux.json.
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
N_MODES = 15    # body supermodes carried at the four MMI-body interfaces
                # (the two connector interfaces carry CONN_MODES instead --
                # see interface_modes below, which has six entries, not four).
                # RAISED 14->15, 2026-08-09 (modal_multiplexer_validation.md
                # gap 7's 2026-08-09 entry): check_mode_correspondence
                # measured a real basis-incompleteness at 14 modes at the
                # short-wavelength band edge (1.26 um, interface 0) -- meep's
                # own mode 8 was only 0.9% representable in meow's 14-mode
                # window, consistent with a 15th-order mode becoming newly
                # guided near that wavelength. Tried 16 first (to keep an
                # even symmetric/antisymmetric PARITY BALANCE, per the
                # comment below) but FDTDSolver's own guided-mode guard
                # raised: this cross-section supplies only 15 truly guided
                # modes at 1.26 um, not 16 -- band 16 sits at/below the
                # cladding index. 15 is therefore the measured ceiling, not
                # a guess; an odd total means one parity has one more mode
                # than the other at this wavelength, which the comment
                # below's "1+5, 2+6" framing did not anticipate but is not
                # itself a problem. Same reasoning as before on why this is
                # safe to test now and wasn't in 2026-08-04: the old padding
                # guard used to silently pad missing modes with junk; it now
                # raises loudly instead (items 12/14), which is exactly how
                # this was discovered. NOT YET VERIFIED against a real run.
                # CHANGED 16->14, 2026-08-04 (modal_multiplexer_validation.md
                # gap 2/gap 6's item-12/14 fix) -- NOT OPTIONAL, a judgment
                # call made without explicit user sign-off, flagged here for
                # review. At RESOLUTION=40 this cross-section supplies only
                # 14 real guided modes (measured cleanly and consistently at
                # EVERY tested resolution 20-100, job 11673526); the library
                # guard that used to pad the missing 2 with unguided/junk
                # modes now raises instead (items 12/14), so 16 would crash
                # on the very next run regardless of this choice. 14, not
                # some other number <16, because it is exactly what gap 2
                # measured as reachable -- reducing further would discard
                # real guided channels for no reason. RESOLUTION stays 40
                # (not raised) because gap 2 also found NO resolution
                # 20-100 is demonstrably converged in absolute n_eff (60->80
                # looked converged, misleadingly -- 80->100 broke it
                # completely), so paying for a finer mesh buys no
                # demonstrated accuracy here; meanwhile this device's own
                # check 1 (gap 5, job 11678100) measured the engine-ratio
                # Δβ_FDTD/Δβ_EME at this exact RESOLUTION=40 and found it
                # excellent (0.9958, matching both MMIs' own shipped-but-
                # unconverged-resolution results) -- the ratio, not the
                # absolute value, is what this hybrid chain actually
                # depends on. NOT independently re-verified: whether 14
                # (down from 16) still keeps several supermodes of EACH
                # parity at every interface -- the reasoning below predicts
                # yes (dropping the two weakest-confined, highest-index
                # modes removes one of each parity from an alternating
                # sequence, not a whole parity), but this has not been
                # checked against a real run of this device, which is a
                # separate, expensive campaign (gap 6) not attempted this
                # session. The legacy hybrid used 20-30 modes; each extra
                # mode used to cost one more FDTD run in the old monolithic
                # transition section (pre-2026-08-06 split) -- since the
                # split, N_MODES no longer drives FDTD run count at all (that
                # is now CONN_MODES, above), but it still sets the size of
                # every body-adjacent mode solve/interface check, so it
                # remains the main accuracy knob even though it stopped being
                # the main cost knob. Don't starve it: the
                # on-axis input is symmetric about y=0, so TE0 excites only
                # *symmetric* body supermodes and TE1 only *antisymmetric*
                # ones (e.g. at the MMI 1 interface TE0 lands in supermodes
                # 1+5, TE1 in 2+6) -- the basis must keep several supermodes
                # of EACH parity or one input mode's channel silently
                # truncates to ~zero.
BODY_STUB = 0.75  # µm of MMI body folded into each FDTD section past the
                  # taper/body discontinuity (legacy mm_monitor_position);
                  # must exceed the FDTD solver's monitor_margin.

CONN_MODES = 6  # arm supermodes carried at the two new connector interfaces.
                # MEASURED 2026-08-06 (login-node EMESolver.section_modes()
                # on this exact cross-section at wl=1.31, x-span matching
                # CONN_STUB below): requesting 8, got n_eff = [2.10549,
                # 2.09067, 2.03595, 1.56570, 1.54514, 1.52526, 1.42755,
                # 1.41349] -- the last two sit below the SOI_160 cladding
                # index (~1.4468, this repo's own prior measurement), so
                # only 6 are guided. Correction to the original hand-wave
                # ("~3, one per arm"): each arm is apparently weakly bimodal
                # (fundamental + first-order), not strictly single-mode,
                # giving 3 arms x 2 modes = 6 coupled supermodes, not 3. Not
                # yet re-confirmed that EMESolver.compute_modes() itself
                # would reject the 2 unguided modes if CONN_MODES were left
                # at 8 -- its documented guard is duplicate/rank-deficiency
                # detection (see its own docstring), NOT the guided-vs-
                # cladding check FDTDSolver._iter_distinct_bands() has (items
                # 12/14); until that's checked, 6 (not 8) is the safe choice,
                # matching this exact measurement rather than trusting an
                # unverified guard to catch the difference.
                #
                # Real cost finding, not just a mode count: this single-rank
                # login-node solve took 5050 s (~84 min) for one wavelength.
                # That is NOT "basically free" the way a single MMI-slab
                # cross-section's mode solve is -- three separated,
                # tightly-pitched silicon islands with no exploitable mirror
                # symmetry is evidently a much harder MPB target. The script
                # solves this cross-section once per wavelength (N_WL=11),
                # so total connector-EME cost in a real run is a genuinely
                # open question, not the near-zero this section's own module-
                # docstring reasoning assumed -- how much production's many-
                # rank parallelism brings this down is UNMEASURED. Revisit
                # before trusting the "cheap EME, expensive FDTD" framing
                # above at face value; this may claw back some of the
                # FDTD-side savings the section-2 split was built to capture.
CONN_STUB = 0.25  # µm of connector folded into each FDTD taper section past
                  # the taper/connector kink. Matches the coupler/other MMI
                  # hybrid examples' own body_stub convention (0.25 um),
                  # deliberately NOT this device's own more conservative
                  # BODY_STUB=0.75 -- that larger clearance was sized for the
                  # MMI body's degenerate multimode discontinuity, a harder
                  # case than the much gentler single-mode taper/straight-run
                  # kink here. Must exceed the FDTD solver's monitor_margin
                  # (0.25, below) and leaves 2*CONN_STUB=0.5 um of
                  # BETWEEN_WG_LENGTH=3.0 um's straight run for the EME
                  # connector section itself (2.5 um) -- comfortable, but
                  # unlike BODY_STUB this exact value has not yet been
                  # checked against HybridSolver(check_mode_correspondence=
                  # True) the way the body interfaces already were (notes.md
                  # Solvers item 14); pass mode_correspondence_pad=0.1 if/
                  # when that check is run here, matching the precedent for
                  # the other examples' thinner 0.25 um stubs.

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
    minimum_run_time=100.0,
    # ADDED 2026-08-08, same fix as modal_multiplexer_fdtd.py's own (gap 3):
    # FDTDSolver never floors stop_when_dft_decayed() by default, so a
    # loose-enough decay_tol can trigger before a section's own light has
    # crossed it. Lower risk here than the full-device script -- every FDTD
    # section in THIS chain is only ~4 um plus stub/PML (~5.5 um total
    # simulated span), a small fraction of the 68 um that bug needed to
    # matter at -- but decay_tol here is already the tight 1e-11 (the
    # regime where the bug did NOT manifest for the FDTD script either), so
    # this is a preventative floor, not a response to an observed failure.
    # 100 = ~4x a rough one-transit estimate (5.5 um * silicon group index
    # ~4.5 =~ 25 meep time units), leaving real margin without approaching
    # any resolution-relevant cost.
)
body_solver = EMESolver(
    num_modes=[N_MODES], num_cells=[1], lengths=[1.0],  # placeholders; solve_section() takes its span from the boundaries
    oxide_width=5.0, oxide_height=1.5, res=RESOLUTION, material_stack=SOI_160,
    mode_solver="meep",  # matches FDTDSolver's own eigenmode basis in this hybrid chain
)

# ---------------------------------------------------------------------------
# 3. Build HybridSolver — 7-section chain:
#    FDTD input | EME MMI 1 | FDTD taper-out | EME conn | FDTD taper-in |
#    EME MMI 2 | FDTD output. (Split 2026-08-06 from a 5-section chain whose
#    middle FDTD block covered both taper bands AND the connector -- see
#    module docstring.)
# ---------------------------------------------------------------------------
half_total = mux.total_length / 2
body1_start = -half_total + TAPER_L
body1_end = body1_start + mmi1.mm_length
conn_start = body1_end + TAPER_L
conn_end = conn_start + BETWEEN_WG_LENGTH
body2_start = conn_end + TAPER_L
body2_end = body2_start + mmi2.mm_length

boundaries = [
    -half_total,               # physical input port (input taper start)
    body1_start + BODY_STUB,   # FDTD input      <-> EME MMI 1
    body1_end - BODY_STUB,     # EME MMI 1       <-> FDTD taper-out
    conn_start + CONN_STUB,    # FDTD taper-out  <-> EME connector
    conn_end - CONN_STUB,      # EME connector   <-> FDTD taper-in
    body2_start + BODY_STUB,   # FDTD taper-in   <-> EME MMI 2
    body2_end - BODY_STUB,     # EME MMI 2       <-> FDTD output
    +half_total,               # physical output port (output taper tips)
]
interface_modes = [N_MODES, N_MODES, CONN_MODES, CONN_MODES, N_MODES, N_MODES]

# Same fdtd_solver/body_solver instances fill every slot of their type (see
# step 2's own comment) -- solve_section() takes n_modes from
# interface_modes above, not from body_solver's own constructor num_modes,
# so no new EMESolver is needed for the connector section.
hybrid = HybridSolver(
    solvers=[fdtd_solver, body_solver, fdtd_solver, body_solver,
             fdtd_solver, body_solver, fdtd_solver],
    boundaries=boundaries,
    interface_modes=interface_modes,
)

# ---------------------------------------------------------------------------
# 4. Geometry plot — shows which solver covers which section. Sections 0–5
#    are identical for every output port, so they render once; the output
#    section (index 6, was 4 before the section-2 split) renders once per
#    port — same full geometry (all three branches) every time, with the
#    reverse source/norm-monitor markers shrunk to that port's launch
#    aperture — exactly the three simulations run in step 5.
#
#    Skipped in per-port array-task mode (SLURM_ARRAY_TASK_ID set — see the
#    "Reshaped 2026-07-31" module docstring): 3 concurrent array tasks would
#    otherwise all write these same shared files at once, a plotting-only
#    race that has no business existing when the merge/serial step (which
#    always runs once, after every array task) can render them exactly once.
# ---------------------------------------------------------------------------
if os.environ.get("SLURM_ARRAY_TASK_ID") is None:
    plot_hybrid_geometry(
        mux, solver=hybrid,
        in_modes=2, out_modes=1,
        title="2-Stage MMI Modal Multiplexer — Hybrid Solver Sections",
        output_path=os.path.join(output_dir, "geometry_hybrid.png"),
        sub_region_dir=os.path.join(output_dir, "sections"),
        sub_region_show_pml=True,
        sub_region_show_eme_detail=True,
        sub_region_configs=[
            (f"port_{name}", {6: {"output_idx": k}}) for name, k in PORTS
        ],
    )
    print(f"Geometry plot saved → {output_dir}/geometry_hybrid.png")
    print(f"Sub-region plots saved → {output_dir}/sections/")

# ---------------------------------------------------------------------------
# 5. Run the hybrid simulation — once per output port. output_idx is passed
#    ONLY to the output section's kwargs, and a shared cache is passed to
#    every solve() call, so sections 0–5 (input FDTD, both EME bodies, both
#    taper FDTDs, the connector EME) are each simulated exactly once — not
#    once per port. in_modes=2 excites TE0 and TE1 of the multi-mode input
#    waveguide (legacy input_modes=[1, 2]); out_modes=1 measures TE0 at each
#    port (legacy output_modes=[1]).
# ---------------------------------------------------------------------------
wavelengths = np.linspace(WL_MIN, WL_MAX, N_WL)
print(f"\nHybrid sweep: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm "
      f"({len(wavelengths)} pts, {N_MODES} body-interface / {CONN_MODES} "
      f"connector-interface modes)\n")

_array_idx = os.environ.get("SLURM_ARRAY_TASK_ID")
if _array_idx is not None:
    # One array task, one output port (see the "Reshaped 2026-07-31" module
    # docstring): simulate just this port's section 6 (was 4 before the
    # section-2 split) and let save_dir persist it (sections 0-5 come from
    # the same on-disk cache, shared across ports/tasks). No SMatrix assembly
    # here -- the merge run (below, no SLURM_ARRAY_TASK_ID) reads every port
    # back out of cache for that.
    _i = int(_array_idx)
    _name, _k = PORTS[_i]
    print(f"array task {_i}: port '{_name}' (output_idx={_k})", flush=True)
    hybrid.solve(
        mux, wavelengths, in_modes=2, out_modes=1, port=_name,
        cache={}, save_dir=save_dir,
        geometry_kwargs={6: {"output_idx": _k}},
    )
    print(f"array task {_i} done -- port '{_name}' cached under {save_dir}", flush=True)
    raise SystemExit(0)

# Merge/serial step: no SLURM_ARRAY_TASK_ID, so every port is expected to
# already be on disk from its own array task. cache_only=True makes that
# expectation explicit -- a CacheMissError here means an array task hasn't
# finished (or failed) rather than silently kicking off a fresh multi-hour
# FDTD run in what is meant to be a cheap assembly step.
section_cache = {}
smat = None
for name, k in PORTS:
    smat_port = hybrid.solve(
        mux, wavelengths, in_modes=2, out_modes=1, port=name,
        cache=section_cache, save_dir=save_dir, cache_only=True,
        geometry_kwargs={6: {"output_idx": k}},
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

plt.suptitle("2-Stage MMI Modal Multiplexer — FDTD/EME/FDTD/EME/FDTD/EME/FDTD Hybrid", fontsize=11)
plt.tight_layout()
# is_master(): from 2026-08-03 this runs under `srun` (see
# batch/run_modal_multiplexer_examples.batch — it had been running 1 of 96
# ranks, which is why 11585847/11590663/11597177 all hit the 5 h wall). Three
# array tasks already share this output_dir, and the batch's own comment worries
# about them racing this file; multiplying each task by 96 ranks makes the gate
# necessary rather than merely tidy.
if is_master():
    fig.savefig(os.path.join(output_dir, "hybrid_results.png"), dpi=300)
plt.close("all")
print(f"Results plot saved → {output_dir}/hybrid_results.png")
