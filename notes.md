# FECES capability gap analysis vs. meep / meow / tidy3d / other open-source EM tools

Research memo — comparison of FECES against meep, meow, tidy3d, and other open-source
electromagnetic solvers, plus a prioritized list of candidate additions. Nothing below
has been implemented; this is a menu, not a commitment.

**What FECES actually is**: a hybrid-methodology orchestration layer on top of meep
(FDTD), MPB/meow (mode solving), and optionally tidy3d's mode solver — plus a small
parametrized photonic-device geometry library (`MMI` with per-port `MMIPort` taper
specs, `DirectionalCoupler`, standalone `Taper`/`SBend` building blocks, `CompoundGeometry`/
`WaveguideRun` for composing several devices into one cascade, and legacy `BendGeometry`).
It is **not** a competing raw solver;
meep/meow/tidy3d do all the actual field/mode math. FECES's real contribution is
`HybridSolver`: a generic N-section FDTD/EME cascade (any solver type, any order, any
count ≥ 2) with reciprocity-based excitation (minimum-simulation-count) and a shared
`section_prisms/section_gds/width_at/center_at` interface that lets any geometry slot
into that cascade.

## Comparison table

| Capability | FECES today | meep | meow | tidy3d | Other OSS |
|---|---|---|---|---|---|
| Time-domain FDTD | 2D or 3D via `FDTDSolver(three_d=...)`, single fixed uniform `res` | Native, mature, adaptive subpixel smoothing | — | Cloud GPU FDTD, auto mesh refinement | openEMS (general 3D FDTD, antenna-focused) |
| Eigenmode expansion (EME) | `EMESolver` wraps meow's cross-section mode-match + `compute_s_matrix` | — | Native (this *is* meow's job) | Has an EME plugin | EMEpy |
| Hybrid FDTD+EME chaining | **Yes — `HybridSolver`, FECES's actual differentiator** | No | No | No (closest: manually stitching, no generic N-section cascade) | Not seen elsewhere in OSS photonics |
| Dispersive materials | Single real index at center λ only (`.real` used, imaginary/loss dropped) | Full Lorentzian/Drude/Debye susceptibility fits | Has real material models (used via meow calls) | Automatic dispersion-fitting plugin | — |
| Nonlinear materials (χ2/χ3) | None | Yes (Kerr, χ2) | No | Yes | — |
| Lossy/absorptive materials | None (imaginary index discarded) | Yes | Partial | Yes | — |
| Anisotropic / magneto-optic | None | Yes | No | Yes | — |
| Material library breadth | 2 stacks: `SOI_160`, `SOI_220` (Si/SiO2 only, thickness-parametrized) | `meep.materials` — dozens of real fits (SiN, III-V, metals, oxides) | Some via meow | Large curated library | — |
| Periodic structures / gratings / PhC | None (no Bloch BCs anywhere) | Yes (Bloch-periodic boundaries) | No | Yes | **S4/RCWA** — the standard OSS tool for this, much faster than FDTD for gratings; **MPB** for full PhC band structure (transitively available via meow but not exposed) |
| Symmetry exploitation | None — always draws/simulates full both-arm geometry even for mirror-symmetric devices | `mp.Mirror`/`mp.Rotate2` — 2-4x speedup for symmetric structures | n/a | Yes | — |
| Adjoint / inverse design / topology optimization | None | Adjoint solver (`meep.adjoint`) | No | `autograd`-based, mature | — |
| Near-to-far-field transform | None | `add_near2far` | No | Yes | — |
| Mesh auto-refinement / convergence checking | None — manual `res`, user must sweep by hand | Manual (subpixel smoothing helps, not auto-refine) | n/a | Automatic adaptive meshing | — |
| Device geometry library | `MMI` (N×M, per-port `MMIPort` taper size/placement), `DirectionalCoupler`, standalone `Taper`/`SBend` building blocks, `CompoundGeometry`/`WaveguideRun` (compose several devices into one cascade — e.g. two `MMI`s + a connector, replacing the old hand-rolled cascaded-MUX class), `BendGeometry` (legacy, poorly integrated) | n/a (raw geometry primitives only) | n/a | Has some parametric primitives + PDK support | **gdsfactory** — the standard OSS parametric photonic component library (rings, Y-branches, grating couplers, edge couplers, crossings, spirals, etc.) |
| GDS import/export, PDK interop | **None** — `gdsfactory` dependency was actively removed (commit `ffcbf5a`); geometries are FECES-internal polygon dicts only, can't simulate a real fab layout or export one for tape-out | n/a | Uses gdsfactory `Component`s | Has gdsfactory integration | gdsfactory itself |
| Ring/racetrack resonator | **Missing entirely** — extremely common building block | n/a | n/a | n/a | gdsfactory has one; easy to add as a new `DeviceGeometry` |
| Grating coupler / edge coupler geometry | Missing | n/a | n/a | n/a | gdsfactory |
| Circuit-level composition (cascade many devices' S-matrices into one chip netlist) | **None** — `SMatrix` has no `+`/`@`/composition operators | n/a | n/a | n/a | **SAX** — JAX-based, differentiable, netlist-driven; the standard OSS tool for this layer |
| Result serialization | Pickle only, no schema/versioning | n/a | n/a | Its own `.hdf5`/JSON | Touchstone `.sNp` is the RF/photonics interop standard, not supported |
| Passivity/unitarity sanity check | None | n/a | n/a | n/a | Trivial to add (`sum(|S|^2) <= 1`); would have caught real bugs found in earlier sessions automatically |
| Parallel parameter sweeps | `ParallelSweep` — local `multiprocessing.Pool` grid sweep, but only calls whole-device `solver.solve()`, **not** `HybridSolver`-aware, no cloud/cluster backend beyond an external SLURM Makefile | n/a (user-driven) | n/a | Cloud batch is tidy3d's core value prop | — |
| Field visualization | Static EME field heatmaps only (`plot_eme_fields`); **no FDTD time-domain field/movie export at all** | `mp.Animate2D`, live field plotting | n/a | Interactive Plotly 3D viewer | — |
| Interactive/3D visualization | None (matplotlib static PNGs only) | Some (via `mp.plot2D`) | n/a | Yes, core feature | — |
| Testing | **None** — no `tests/`, no CI, no lint | n/a | n/a | n/a | — |

## Current library reference

Descriptive reference for what's actually in the library today — not recommendations,
just what each solver/geometry/plotting entry point does and how they fit together.

### Solvers (`feces/solvers/`)

- **`FDTDSolver`** — meep-based FDTD. `res`/`sim_width`/`sim_height`/`dpml` configure
  the mesh/domain; `three_d` switches between a cheap legacy 2-D planar simulation
  (bulk-index waveguides, no vertical confinement — quantitatively wrong, e.g.
  under-predicts a directional coupler's coupling rate by ~2×) and a full 3-D
  simulation (real `material_stack.wg_height`-thick silicon core in oxide, PML on all
  six faces, (y,z) cross-section eigenmode solves — the same mode basis EMESolver's
  meow cross-sections use, which is what makes 3-D FDTD sections numerically
  consistent with EME sections in a hybrid chain). `solve()` runs a complete,
  self-contained whole-device simulation. `solve_section()` is the piece
  `HybridSolver` calls instead: given any `DeviceGeometry` and an arbitrary
  `[x0, x1]` sub-range (+ optional `stub_left`/`stub_right`), it builds a small
  standalone FDTD cell for just that section, using reciprocity to minimize run
  count (always excites from whichever side has fewer modes, reversing the physical
  direction if needed) and sizing the source/monitor aperture via a priority order:
  geometry-provided `source_aperture_at()` first, then `half_domain_aperture`, then
  the full `sim_width`.
- **`EMESolver`** — meow-based eigenmode expansion. `num_modes`/`num_cells`/`lengths`
  (parallel lists, one entry per region left-to-right, e.g. `[wg, taper, MMI body,
  taper, wg]`) plus `oxide_width`/`oxide_height`/`res`/`material_stack` configure the
  cross-section mesh; `mode_solver` picks `"meep"` (matches FDTDSolver's own basis —
  use this in hybrid chains to avoid basis mismatches at FDTD/EME interfaces) or
  `"tidy3d"` (fine for an EME-only solve with no FDTD sections to stay consistent
  with). `solve()`/`solve_sweep()` run a complete whole-device eigenmode expansion
  (the latter batches wavelengths for memory efficiency). `solve_section()` /
  `solve_section_length_sweep()` are the `HybridSolver`-facing single-cell-per-call
  variants — the latter solves modes once per wavelength and reuses them across many
  trial lengths for free, since a uniform cross-section's modes don't depend on cell
  length (only the analytic propagation phase does) — this is what backs the
  coupling-length sweep in `directional_coupler_hybrid.py`.
- **`HybridSolver`** — FECES's actual differentiator: a generic N-section (≥2)
  FDTD/EME cascade over a single `DeviceGeometry`. Constructed from `solvers` (any
  mix/order of solver instances implementing `solve_section()`), `boundaries`
  (`len(solvers)+1` x-cuts in device frame), and `interface_modes` (mode count
  carried at each internal cut). `solve()` calls every section's `solve_section()`
  (forwarding `geometry_kwargs` — one shared dict or a per-section list — and an
  optional `cache` dict keyed by solver identity/x-range/mode-counts/kwargs, so a
  section whose kwargs don't vary across multiple `solve()` calls — e.g. once per
  output port — is only ever actually simulated once), then cascades the per-section
  transfer matrices via an `einsum` matrix-chain product contracted over each shared
  interface's mode index.

### Geometries (`feces/geometry/`)

- **`DeviceGeometry`** (`base.py`) — the abstract interface every device
  implements: `polygons()`/`to_meep_prisms()` for whole-device FDTD-only paths;
  `to_gds_eme()` for the meow/EME layout; the "arbitrary section" interface
  `HybridSolver` actually needs — `section_prisms()`/`section_gds()`/`width_at()`/
  `center_at()`/`source_aperture_at()` — letting any solver ask "what's the geometry
  between x0 and x1" without knowing what feature lives there; and `total_length`/
  `eme_lengths()`/`bbox`. `section_prisms`/`section_gds` must raise `ValueError` on a
  cut landing strictly inside a varying-cross-section region (a taper or bend).
- **`MMI`** + **`MMIPort`** (`mmi.py`) — N-input/M-output multimode interferometer.
  Each port now carries its own `MMIPort(y, waveguide_width, taper_length,
  taper_width)`, so ports can have individually different taper size and arbitrary
  y-placement (the constructor's old scalar `taper_length`/`taper_width`/
  `waveguide_offset` args still work as a convenience default that auto-builds
  symmetric, evenly-spaced `MMIPort` lists). `use_neighbor_aperture` (opt-in,
  default off) generalizes the old `MMIMux`'s neighbor-distance aperture-narrowing
  formula to arbitrary port counts, for devices with several closely-spaced ports
  (e.g. one sitting exactly on the y=0 axis) that the fixed `half_domain_aperture`
  rule can't isolate.
- **`DirectionalCoupler`** (`coupler.py`) — symmetric 2-waveguide coupler with
  S-bend I/O transitions; built internally from 4 `SBend` sub-geometries (one per
  arm), composed via the shared `_compose.shift_polys`/`shift_prisms` offset
  helpers, rather than hand-rolled polygon math per arm.
- **`Taper`** (`taper.py`) — standalone linear width-transition `DeviceGeometry`
  (width_start → width_end over a length). Used internally by `MMI` for every
  port's taper; also directly usable in any `HybridSolver`/`FDTDSolver`/`EMESolver`
  chain on its own.
- **`SBend`** (`sbend.py`) — standalone lateral-S-bend `DeviceGeometry` between two
  y-offsets over a length, `"euler"` (cubic-Bezier) or `"circular"` (matched
  constant-curvature arcs) profile. Used internally by `DirectionalCoupler` for
  each of its 4 arms; also directly usable standalone.
- **`CompoundGeometry`** + **`WaveguideRun`** (`compound.py`) — concatenates several
  `DeviceGeometry` instances (and/or `WaveguideRun` straight multi-port connectors)
  end-to-end into a single `DeviceGeometry`, delegating every hook to whichever
  sub-geometry owns the queried x. This is what lets `HybridSolver` (which takes
  exactly one geometry object per cascade — see Solvers above) simulate a device
  built from several independent parts: `examples/modal_multiplexer/
  modal_multiplexer_hybrid.py` composes a 1×3 `MMI` + a `WaveguideRun` connector +
  a 3×3 `MMI` this way, replacing the old hand-rolled cascaded-MUX class.
- **`BendGeometry`** (`bend.py`) — legacy standalone 90° routing bend for
  non-collinear FDTD I/O; predates and isn't integrated with the `section_prisms`/
  `width_at`/`center_at` interface (see Geometry recommendations above).
- *Removed*: `ModalMultiplexer` (unused dead code) and `MMIMux`/`MuxGeoParams` (the
  old cascaded 1×3+3×3 MUX class) — superseded by composing two `MMI` instances via
  `CompoundGeometry`, as above.

### Plotting (`feces/plot/`)

- **`plot_eme_geometry`** — top-down polygon plot for a standalone `EMESolver`:
  draws `geometry.to_gds_eme()` plus one shaded/labeled region per `solver.lengths`
  entry — exactly the cells `solve()`/`solve_sweep()` actually mode-solve.
- **`plot_fdtd_geometry`** — top-down polygon plot for a standalone `FDTDSolver`:
  draws `to_gds_eme()` plus the real source/monitor markers and PML shading
  `solve()` sets up (source aperture from `source_aperture_at()` /
  `half_domain_aperture` / full `sim_width`, in the same priority order
  `solve_section()` itself uses; output monitors from `output_monitor_centers`/
  `output_monitor_widths`).
- **`plot_hybrid_geometry`** — top-down polygon plot for a `HybridSolver`: one
  shaded block per section plus boundary lines. `sub_region_dir` additionally
  renders a zoomed-in plot per section — an FDTD section gets `solve_section()`-
  accurate source/monitor markers (sized via `width_at()`/`center_at()`, same
  reciprocity-based excitation-side rule); an EME section gets a shaded region,
  optionally broken down into its `num_modes`/`num_cells` detail.
  `sub_region_configs` lets one call render every output-port variant (e.g. a
  different `output_idx` per port) in a single pass, re-using a section's render
  across configs whenever its kwargs don't actually change between them — this is
  what the geometry plot in `modal_multiplexer_hybrid.py` uses to render three
  per-port output sections without re-drawing sections 0–3 three times over.
- **`plot_eme_fields`** — EME-only: propagates modes through meow cells and saves
  6-component (`|Ex|`,`|Ey|`,`|Ez|`,`|Hx|`,`|Hy|`,`|Hz|`) (x,z) field heatmaps. No
  FDTD time-domain equivalent exists yet (see Plotting recommendations above).
- **`plot_performance`** — transmission (linear) + insertion loss (dB) vs
  wavelength, side by side, for every port in an `SMatrix`.
- **`plot_splitting_ratio`** — power splitting ratio `T_a/(T_a+T_b)` vs wavelength
  between two named ports, with a 50/50 reference line.

## Recommendations

Sorted into High-level / Solver / Geometry / Plotting / Examples / Misc. Ordered
within each section by (impact) × (leverage — how much is genuinely missing vs. just
"meep already does this, expose it"), not by raw novelty. Sections are kept even
where a category has fewer or lower-priority items.

### High-level

1. **Test suite + CI.** Zero tests exist anywhere in the repo. Prior debugging
   sessions found and fixed real, silent correctness bugs (FDTD launch aperture
   width, interior-interface measurement aperture truncating supermode overlap, a
   plot-only marker-position bug) — every one would have been caught immediately by
   a regression test against a known-good value. Cheapest high-value additions:
   (a) fast `solve_section()`-level unit tests with mocked/tiny geometries, no real
   FDTD/EME needed; (b) one slow "golden value" integration test per example device
   that runs a real (reduced-resolution) sim and checks T/IL against a recorded
   reference within tolerance; (c) GitHub Actions CI running both on push.
2. **Documentation.** No `README`, no docs site, no changelog — the `examples/`
   scripts are the only tutorials that exist. At minimum a top-level `README.md`
   (install, the FDTD/EME/Hybrid split, a pointer to `examples/`) and a
   `CHANGELOG.md` tracking the still-static `__version__ = "0.1.0"`.
3. **Type checking in CI.** Type hints are used throughout (`from __future__ import
   annotations`) but nothing enforces them — no `mypy`/`pyright` config, no `py.typed`
   marker. Cheap to add once CI exists (item 1c).
4. **Declare transitive dependencies explicitly.** `gdsfactory`/`tidy3d` are consumed
   indirectly through `meow-sim` (e.g. `meow.compute_modes_tidy3d`) but aren't
   declared in `pyproject.toml` — a `meow-sim` version bump that drops/renames those
   could silently break FECES. Pin or declare them directly, or guard the import
   with a clear error message.
5. **Standardize on `logging`, not `print()`.** `EMESolver.solve_section_length_sweep`
   and others `print()` diagnostics unconditionally (e.g. the mid-sweep neff
   diagnostic) while `ParallelSweep` already uses the `logging` module — inconsistent
   for a library whose callers may want to suppress or redirect that output.

### Solver

1. **Passivity/unitarity self-check on `SMatrix` — done.** `smat.check_passive(tol=
   1e-3)` (`feces/results/smatrix.py`) sums `|coeff|^2` over every (port, out_mode)
   entry per input mode present (via the new `total_power(in_mode)` helper) and
   warns (`warnings.warn`, identifying the worst input mode/wavelength/total) if any
   input mode's total output power exceeds `1+tol` at any wavelength — exactly the
   kind of check that would have flagged this session's own -200dB/+100dB/mode-
   truncation-style bugs immediately, without needing to notice an anomalous plot
   first. Not yet wired into any example's own print/verify block as a routine
   post-solve call — worth adding `smat.check_passive()` to each example script
   right after `solver.solve()`/`hybrid.solve()` as standard practice.
2. **Lossy/dispersive material support.** The imaginary part of the refractive index
   is discarded (`materials/stack.py`, `si_medium`/`sio2_medium` use `.real` only)
   and only a single real index at the center wavelength is used per solve — no true
   dispersion across a wavelength sweep, no absorption. Wiring meep's actual
   Lorentzian-susceptibility fitting (or reusing `meep.materials` directly) into
   `MaterialStack` would enable real loss/broadband dispersion and open the door to
   SiN, III-V, or lossy/doped silicon.
3. **Symmetry exploitation in `FDTDSolver`.** `MMI`/`DirectionalCoupler` are
   frequently mirror-symmetric about y=0; `solve_section()` currently always draws
   and simulates the full cross-section. An opt-in `mp.Mirror` symmetry when the
   geometry/excitation is actually symmetric could roughly halve FDTD cost for free —
   complements the reciprocity-based min-simulation-count optimization already there.
4. **Circuit-level S-matrix composition.** `SMatrix` has no way to combine several
   devices' results into a larger circuit (e.g. cascading an MMI's `SMatrix` with a
   directional coupler's). `HybridSolver` already does exactly this kind of
   matrix-chain cascading *within* one device across sections — generalizing that
   same machinery to cascade *across* independently-solved devices (or documenting a
   recipe for handing `SMatrix` results to SAX, which already does this for the
   gdsfactory ecosystem) would stop device-level sims from being a dead end. (Note:
   `CompoundGeometry` — see Geometry section — solves a related but distinct problem
   at the geometry layer: it lets `HybridSolver` simulate a device built from several
   `DeviceGeometry` parts as one cascade *before* solving, but doesn't let two
   already-*solved* `SMatrix` results be combined after the fact, which is what this
   item is about.)
5. **`ParallelSweep` + `HybridSolver` integration.** The sweep utility only calls
   whole-device `solver.solve()`; it doesn't know about `HybridSolver.solve()`'s
   `cache=` mechanism, so a parallel sweep over a hybrid device would redundantly
   re-simulate cached sections across worker processes.
6. **Non-uniform/graded meshing.** `res` is a single isotropic scalar shared across
   the whole mesh (both FDTD and EME's `Mesh2D`) — no finer mesh near the waveguide
   core and coarser in cladding, unlike tidy3d's automatic adaptive meshing. Would
   cut cell count for a given accuracy target.
7. **Automatic mesh-resolution convergence.** Beyond the manual "run a list of
   resolutions and eyeball it" pattern (see Misc/Examples below) — meep/tidy3d both
   have some form of resolution auto-selection; FECES has none, always a fixed
   user-chosen `res`.
8. **Bigger lift, only if scope grows toward design automation:** adjoint/
   inverse-design support (meep adjoint or tidy3d-style autograd), periodic/grating
   structures via an RCWA backend (S4-style) for problems where FDTD/EME are the
   wrong tool, near-to-far-field transform. Substantial, well-trodden ground already
   owned by meep/tidy3d/S4 respectively — worth it only if a concrete device need
   shows up that FDTD/EME genuinely can't address.

### Geometry

1. **Ring/racetrack resonator.** The single most common photonics building block
   that isn't in the device library at all. Natural fit for the existing
   `DeviceGeometry` interface (bus waveguide + coupled ring, `section_prisms`/
   `width_at`/`center_at` following the `DirectionalCoupler` pattern), and directly
   exercises `HybridSolver`'s generic cascading (bus straights as FDTD/EME + coupling
   region).
2. **GDS import (read a real/foundry/gdsfactory layout) and export (tape-out-ready
   output).** The gdsfactory dependency was deliberately removed from the internal
   pipeline (raw polygon dicts are faster/simpler for FECES's own parametric shapes),
   but that leaves no way to simulate an arbitrary externally-supplied layout or
   export a FECES geometry as a real `.gds`. A thin, optional adapter (gdsfactory
   `Component` ↔ FECES polygon dict) would restore that interop without
   reintroducing gdsfactory as a hard dependency.
3. **Y-branch / adiabatic splitter and waveguide crossing.** Common alternatives to
   the MMI splitter and to routing waveguides through each other — neither exists.
4. **Consolidate the duplicated S-bend logic — done for the coupler.py/modal_mux.py
   split.** S-bend centerline/extrusion math now lives solely in a standalone `SBend`
   `DeviceGeometry` (`sbend.py`), composed internally by `DirectionalCoupler` for each
   of its 4 arms via the shared `_compose.shift_polys`/`shift_prisms` offset helpers
   (`ModalMultiplexer`, the only other consumer of the old private functions, was
   deleted as dead code). Still outstanding: `BendGeometry` (`bend.py`) remains a
   *third*, separately-implemented 90° bend — a legacy standalone GDS-writer,
   predating and not integrated with the shared `section_prisms`/`width_at`/
   `center_at` interface at all, and its `"euler"` mode isn't actually a distinct
   curve from `"circular"` (unlike `SBend`'s genuine Bezier/circular-arc split).
   Worth rebuilding it on top of `SBend` (or a new arc-profile primitive) so every
   bend in the library shares one implementation.
5. **General-purpose parametric taper class — done.** A standalone `Taper`
   `DeviceGeometry` (`taper.py`) now exists and is composed internally by `MMI` for
   every port's taper instead of inline trapezoid math duplicated 3x per class. It's
   linear-profile only for now; a polynomial or Euler-spiral width profile (for
   lower-loss adiabatic tapers) would be a natural follow-up, and no example yet
   demonstrates `Taper`/`SBend` used standalone outside their host devices (see
   Examples section).

### Plotting

1. **Stop duplicating solver position/size math in the plot functions.**
   `_draw_section_fdtd_markers` (plot/geometry.py) independently recomputes where
   `FDTDSolver.solve_section()` places its source/monitor apertures, and this
   session's own bug-hunt found two real mismatches there (marker offset by
   `monitor_margin` when the solver places it flush at the section edge; marker
   always drawn full-height when the solver only does that for `n_measure > 1`).
   Any future change to `solve_section()`'s aperture logic can silently desync the
   plot again. Factor the position/size computation into a shared helper both
   `solve_section()` and the plotting code call, so they can't drift apart.
2. **FDTD time-domain field/movie export.** `plot_eme_fields` only covers EME (static
   `|component|(x,z)` heatmaps via `mw.propagate_modes`); there's no equivalent for
   FDTD at all — no meep `Animate2D`-style field snapshot or movie export anywhere in
   the package.
3. **Phase / group-delay plots.** `SMatrix` stores full complex coefficients, but
   `plot_performance`/`plot_splitting_ratio` only ever plot `|coeff|^2` (linear T) and
   dB insertion loss — phase and group delay (from phase vs. wavelength) are never
   visualized despite being available in the data.
4. **A shared convergence plot.** Complements the Misc/Examples resolution-sweep
   item: a `plot_convergence(results_by_resolution, metric)` helper that overlays a
   chosen quantity (e.g. L_c, T) across resolutions on one figure, instead of
   eyeballing separate per-resolution plots side by side (what
   `directional_coupler_mode.ipynb` does today).
5. **Interactive/3D visualization.** Everything is a static matplotlib PNG; no
   plotly/bokeh option for geometry or field plots, unlike tidy3d's interactive
   viewer.

### Examples

1. **No standalone example for `Taper`/`SBend`/`CompoundGeometry` used directly.**
   All three currently only appear as internal building blocks (inside `MMI`,
   `DirectionalCoupler`, and the modal-multiplexer script's two-`MMI`-plus-connector
   composition) — a minimal script instantiating one on its own through a lone
   `FDTDSolver`/`HybridSolver` section would double as living documentation for the
   composition pattern.
2. **No example exists for `BendGeometry`** at all, despite it being a public,
   exported class.
3. **No 2×2 or 1×4 MMI example** — `MMI` supports arbitrary `n_inputs`/`n_outputs`
   but every example hardcodes 1×2.
4. **No example uses `ParallelSweep`.** The library ships a generic parallel
   parameter-sweep utility that zero examples demonstrate — a user has no worked
   reference for the one utility meant to replace ad hoc per-device sweep scripts.
5. **No example demonstrates mesh-convergence checking as a reusable pattern** beyond
   the one-off resolution loop hand-built into `directional_coupler_mode.ipynb` this
   session — worth generalizing into its own example once the Misc convergence
   utility exists.
6. **Ring resonator example** — once the Geometry-section ring/racetrack class
   exists, pair it with an example (and it'd make a good first user of the
   convergence-utility and passivity-check recommendations above, exercising three
   new features at once).

### Misc

1. **Mesh-convergence-loop utility.** The pattern hand-rolled in
   `directional_coupler_mode.ipynb` (loop over a resolution list, save/plot each
   into its own subfolder) should be a reusable `feces` utility (e.g.
   `sweep_resolution(run_fn, resolutions) -> convergence data`) rather than
   hand-rolled per example.
2. **Result serialization beyond raw pickle.** `SMatrix.save`/`load` use unversioned
   `pickle` — no schema, no forward/backward-compatibility story, and no interop
   format. Touchstone `.sNp` is the standard RF/photonics interop format and isn't
   supported; a JSON or HDF5 option would also let results be inspected without
   unpickling FECES's own classes.
3. **Declarative simulation runner.** Every example is a bespoke Python script; a
   YAML/JSON "device + solver + sweep" spec (in the spirit of tidy3d's
   JSON-serializable `Simulation` objects) would let standard sweeps run without
   writing Python, and would pair naturally with `ParallelSweep`.
4. **Output-folder conventions.** `_utils.py` already has `resolve_output_folder`/
   `make_output_folders`, but usage is inconsistent across examples (some use them,
   some hand-roll `os.makedirs`/`output_dir` logic, e.g. every `examples/*.py`
   script this session touched) — worth consolidating on the existing utility
   everywhere instead of ad hoc per-script path-building.

## Not recommending

- Building FECES's own raw FDTD/EME numerical core — meep and MPB/meow already do
  this well; FECES's value is the orchestration layer on top, not the field solver.
- Matching tidy3d's cloud/GPU infrastructure — out of scope for what is currently a
  single-machine/HPC-Makefile-driven library.
