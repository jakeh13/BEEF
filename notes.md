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
| Time-domain FDTD | 2D or 3D via `FDTDSolver(three_d=...)`, single fixed uniform `res` | Native, mature, adaptive subpixel smoothing | — | Cloud GPU FDTD, auto mesh refinement | openEMS (general 3D FDTD, antenna-focused); **FDTDX** (new, JAX-based, GPU-accelerated, native autodiff via time-reversibility — 10-415x faster than traditional FDTD per its own benchmarks) |
| Eigenmode expansion (EME) | `EMESolver` wraps meow's cross-section mode-match + `compute_s_matrix` | — | Native (this *is* meow's job) | Has an EME plugin | EMEpy |
| Hybrid FDTD+EME chaining | **Yes — `HybridSolver`, FECES's actual differentiator** | No | No | No (closest: manually stitching, no generic N-section cascade) | Not seen elsewhere in OSS photonics |
| Dispersive materials | Nearest-tabulated-wavelength real-index lookup (`np.argmin` on meow's built-in Si/SiO2 tables) — no interpolation, no loss | Full Lorentzian/Drude/Debye susceptibility fits | Has real material models (used via meow calls) | Automatic dispersion-fitting plugin | — |
| Nonlinear materials (χ2/χ3) | None | Yes (Kerr, χ2) | No | Yes | — |
| Lossy/absorptive materials | None (imaginary index discarded) | Yes | Partial | Yes | — |
| Anisotropic / magneto-optic | None | Yes | No | Yes | — |
| Material library breadth | 2 stacks: `SOI_160`, `SOI_220` (Si/SiO2 only, thickness-parametrized) | `meep.materials` — dozens of real fits (SiN, III-V, metals, oxides) | Some via meow | Large curated library | — |
| Periodic structures / gratings / PhC | None (no Bloch BCs anywhere) | Yes (Bloch-periodic boundaries) | No | Yes | **S4/RCWA** — the standard OSS tool for this, much faster than FDTD for gratings; **MPB** for full PhC band structure (transitively available via meow but not exposed) |
| Symmetry exploitation | None — always draws/simulates full both-arm geometry even for mirror-symmetric devices | `mp.Mirror`/`mp.Rotate2` — 2-4x speedup for symmetric structures | n/a | Yes | — |
| Adjoint / inverse design / topology optimization | None | Adjoint solver (`meep.adjoint`, autograd/JAX-based) | No | Native `autograd` support (the older JAX-based adjoint plugin is being deprecated in favor of this) | **FDTDX** — JAX-native FDTD with automatic differentiation built into the solver itself, not bolted on afterward |
| Near-to-far-field transform | None | `add_near2far` | No | Yes | — |
| Mesh auto-refinement / convergence checking | None — manual `res`, user must sweep by hand | Manual (subpixel smoothing helps, not auto-refine) | n/a | Automatic adaptive meshing | — |
| Device geometry library | `MMI` (N×M, per-port `MMIPort` taper size/placement), `DirectionalCoupler`, standalone `Taper`/`SBend` building blocks, `CompoundGeometry`/`WaveguideRun` (compose several devices into one cascade — e.g. two `MMI`s + a connector, replacing the old hand-rolled cascaded-MUX class), `BendGeometry` (legacy, poorly integrated) | n/a (raw geometry primitives only) | n/a | Has some parametric primitives + PDK support | **gdsfactory** — the standard OSS parametric photonic component library (rings, Y-branches, grating couplers, edge couplers, crossings, spirals, etc.) |
| GDS import/export, PDK interop | **None** — `gdsfactory` dependency was actively removed (commit `ffcbf5a`); geometries are FECES-internal polygon dicts only, can't simulate a real fab layout or export one for tape-out | n/a | Uses gdsfactory `Component`s | Has gdsfactory integration | gdsfactory itself |
| Ring/racetrack resonator | **Missing entirely** — extremely common building block | n/a | n/a | n/a | gdsfactory has one; easy to add as a new `DeviceGeometry` |
| Grating coupler / edge coupler geometry | Missing | n/a | n/a | n/a | gdsfactory |
| Circuit-level composition (cascade many devices' S-matrices into one chip netlist) | **None** — `SMatrix` has no `+`/`@`/composition operators | n/a | n/a | n/a | **SAX** — JAX-based, differentiable, netlist-driven; **Simphony** — component-library/plugin-based, ~20x faster than a leading commercial tool per its own benchmarks; both are standard OSS choices for this layer with different tradeoffs (functional/autodiff vs. component-abstraction) |
| Result serialization | Pickle only, no schema/versioning | n/a | n/a | Its own `.hdf5`/JSON | Touchstone `.sNp` is the RF/photonics interop standard, not supported |
| Passivity/unitarity sanity check | None | n/a | n/a | n/a | Trivial to add (`sum(|S|^2) <= 1`); would have caught real bugs found in earlier sessions automatically |
| Parallel parameter sweeps | `ParallelSweep` — local `multiprocessing.Pool` grid sweep, but only calls whole-device `solver.solve()`, **not** `HybridSolver`-aware, no cloud/cluster backend beyond an external SLURM Makefile | n/a (user-driven) | n/a | Cloud batch is tidy3d's core value prop | — |
| Field visualization | Static EME field heatmaps only (`plot_eme_fields`); **no FDTD time-domain field/movie export at all** | `mp.Animate2D`, live field plotting | n/a | Interactive Plotly 3D viewer | — |
| Interactive/3D visualization | None (matplotlib static PNGs only) | Some (via `mp.plot2D`) | n/a | Yes, core feature | — |
| Testing | Fast unit tests (`tests/`, meep/meow-free) run in CI on every push; 4 slow golden-value integration tests (real reduced-res meep/meow sims) run in CI as a manual-only `workflow_dispatch` matrix via conda-forge `pymeep` (`pytest --runslow` locally); no lint/type-checking yet | n/a | n/a | n/a | — |

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
  self-contained whole-device simulation, and (like `solve_section()`) now accepts
  an optional `save_dir` for the persistent on-disk cache — see Solver
  recommendation #9. `solve_section()` is the piece
  `HybridSolver` calls instead: given any `DeviceGeometry` and an arbitrary
  `[x0, x1]` sub-range (+ optional `stub_left`/`stub_right`), it builds a small
  standalone FDTD cell for just that section, using reciprocity to minimize run
  count (always excites from whichever side has fewer modes, reversing the physical
  direction if needed) and sizing the source/monitor aperture via a priority order:
  geometry-provided `source_aperture_at()` first, then `half_domain_aperture`, then
  the full `sim_width`. `monitor_center`/`monitor_size` explicitly override the
  measured/interface-side monitor's placement/size when the geometry-driven default
  would be off-axis relative to an adjacent section's own basis (e.g. re-centering
  on the shared y=0 axis at a `DirectionalCoupler` S-bend/coupling-region seam).
- **`EMESolver`** — meow-based eigenmode expansion. `num_modes`/`num_cells`/`lengths`
  (parallel lists, one entry per region left-to-right, e.g. `[wg, taper, MMI body,
  taper, wg]`) plus `oxide_width`/`oxide_height`/`res`/`material_stack` configure the
  cross-section mesh; `mode_solver` picks `"meep"` (matches FDTDSolver's own basis —
  use this in hybrid chains to avoid basis mismatches at FDTD/EME interfaces) or
  `"tidy3d"` (fine for an EME-only solve with no FDTD sections to stay consistent
  with). `solve()`/`solve_sweep()` run a complete whole-device eigenmode expansion
  (the latter batches wavelengths for memory efficiency). `solve()` (unlike
  `solve_sweep()`) now also accepts an optional `save_dir`, cached per output
  port/side since each side builds its own geometry/cells — see Solver
  recommendation #9. `solve_section()` /
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
  (forwarding `geometry_kwargs` — one shared dict, a per-section list, or a sparse
  `{section_idx: kwargs}` dict naming only the sections that need kwargs (unnamed
  sections get `{}`; int keys always mean the sparse form) — and an
  optional `cache` dict keyed by solver identity/x-range/mode-counts/kwargs, so a
  section whose kwargs don't vary across multiple `solve()` calls — e.g. once per
  output port — is only ever actually simulated once), then cascades the per-section
  transfer matrices via an `einsum` matrix-chain product contracted over each shared
  interface's mode index. Optional `mirror_pairs={k: j}` (plus `mirror_axis`,
  `mirror_check_samples`/`mirror_check_tol`) declares that section `k` is the
  x-mirror of already-solved section `j` — its transfer matrix is reused as
  `mats[j].transpose(0,2,1)` instead of a real `solve_section()` call, guarded
  every `solve()` call by sampling `geometry.width_at`/`center_at` across both
  sections' spans under x-reversal (raises `ValueError` on a bad pairing).
  `solve_length_sweep(geometry, wavelengths, lengths, sweep_section, ...)` sweeps
  one section's length analytically — that section's solver must implement
  `solve_section_length_sweep()` (EMESolver does, FDTDSolver doesn't), span a
  uniform cross-section (guardrail-sampled), and be unmirrored; returns the raw
  `[n_L, M, out, in]` cascade grid, and non-swept sections' cache keys match
  `solve()`'s exactly so entries dedupe across both methods. See Solver
  recommendation #7 (below) for the full design rationale; worked examples of
  both features live in `tests/test_hybrid.py` (`TestMirrorPairs`/
  `TestLengthSweep`) — the hybrid examples themselves (`mmi_2x2_hybrid.py`,
  `directional_coupler_hybrid.py`) hand-roll the same reuse manually because
  their cross/through paths mix two excitations' matrices in one cascade,
  which the one-call contract can't express (see recommendation #7's "still
  open").

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
  rule can't isolate. `polygons()`/`to_gds_eme()` now also take `input_indices`/
  `input_idx` (mirroring the existing `output_indices`/`output_idx`), restricting
  the drawn input stub(s) to a single port — needed for a genuine multi-input
  device (2×2, etc.) so an EME input cross-section is one waveguide, not an
  unwanted multi-guide supermode; `examples/mmi_2x2/` is the first consumer.
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
  non-collinear FDTD I/O; **not even a `DeviceGeometry` subclass** (a plain
  dataclass with its own hand-rolled minimal GDSII binary writer,
  `.build(device_sx, dpml, output_path)`), so it predates and isn't
  integrated with the `section_prisms`/`width_at`/`center_at` interface at
  all — can't participate in a `HybridSolver` chain like every other
  geometry in this list (see Geometry recommendations below).
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
  FDTD time-domain equivalent exists yet (see Plotting recommendations below).
- **`plot_performance`** — transmission (linear) + insertion loss (dB) vs
  wavelength, side by side, for every port in an `SMatrix`.
- **`plot_splitting_ratio`** — power splitting ratio `T_a/(T_a+T_b)` vs wavelength
  between two named ports, with a 50/50 reference line.

### Results (`feces/results/`)

- **`SMatrix`** (`smatrix.py`) — unified S-parameter result type, keyed by
  `(port, in_mode, out_mode)` → a complex amplitude array over wavelength.
  `transmission()`/`insertion_loss_db()` derive `|coeff|^2` (linear) and dB
  from the stored complex coefficients. `check_passive(tol=1e-3)` sums
  `|coeff|^2` over every `(port, out_mode)` entry present for a given input
  mode (via the `total_power(in_mode)` helper) and warns
  (`warnings.warn`, identifying the worst input mode/wavelength/total) if
  the total exceeds `1+tol` at any wavelength — a passive device can only
  lose power, never gain it, so a violation almost always means a real bug
  (mode-count mismatch, double-counted output branch, normalization error)
  rather than genuine gain. Every example script and all 4 golden tests call
  `smat.check_passive()` right after `solve()`/`hybrid.solve()` as standard
  practice. `save()`/`load()` round-trip via pickle (unversioned — see Misc
  recommendations for the interop/schema gap).

### Sweep (`feces/sweep/`)

- **`ParallelSweep`** (`parallel.py`) — `ParallelSweep(solver_fn, geometry_fn,
  param_grid, wavelengths, n_workers=4, output_dir=None, save_per_point=True)`.
  `run()` builds the Cartesian product of `param_grid` (`itertools.product`)
  and dispatches each point to a `multiprocessing.Pool(n_workers)` worker,
  which builds a fresh geometry/solver instance and calls plain
  `solver.solve(geometry, wavelengths)` — **always the whole-device `solve()`
  API, never `HybridSolver.solve()`'s `cache=` mechanism** (see Solver
  recommendations), so sweeping a hybrid device redundantly resimulates every
  cached section once per worker/grid point. Returns `Dict[param_tuple,
  SMatrix]`; optionally pickles each point's `SMatrix` under
  `output_dir/<param-tag>/data/smatrix.pkl` and a `sweep_results.pkl`
  summary. No cloud/cluster backend beyond the external SLURM `Makefile` in
  `examples/`.

### Testing (`tests/`)

- **Fast unit tests** — `test_geometry.py`, `test_smatrix.py`, `test_hybrid.py`
  (53 tests). No meep/meow import anywhere (geometry section-cut math is pure
  Python; `HybridSolver` is tested against a `FakeSolver` with scripted
  `solve_section()` returns, not real physics), so these run everywhere in
  well under a second. Run in CI on every push/PR (`.github/workflows/ci.yml`).
- **Slow golden-value integration tests** — `test_golden_examples.py`, one per
  target: `mmi_1x2_eme`/`mmi_1x2_fdtd`/`mmi_1x2_hybrid`/
  `directional_coupler_hybrid` (not `mmi_2x2_*` or `modal_multiplexer_hybrid`
  yet). Each builds the same device geometry as its example at reduced
  resolution (`res=20` vs the examples' `res=40`, 3 wavelengths vs 11), runs a
  real solve, and checks transmission against a recorded reference in
  `tests/golden/*.json` within a loose (10% relative / 3% absolute) tolerance
  — a regression guard against real correctness bugs (sign flips, wrong port
  routing, normalization mistakes — see `SMatrix.check_passive()` in the
  Results reference above), not a physical-accuracy check. Skipped unless
  `pytest --runslow` is passed; `scripts/generate_golden_values.py`
  (re)generates the reference JSON. Genuinely slow even at reduced
  resolution — real per-target wall-clock on a single CPU core: ~7 min
  (EME), ~12 min (FDTD), ~58 min (hybrid), ~96 min (directional coupler
  hybrid) — so budget ~3 hours for the full suite and don't expect this to
  ever be a pre-commit check. The two hybrid targets needed `dpml` doubled
  (0.5→1.0 µm) versus their examples' value: `HybridSolver`'s interior
  `solve_section()` spans are short enough that PML is a large fraction of
  the domain, and at `res=20` the unmodified `dpml` bought only half
  production's PML pixel count, which diverged (`meep: fields NaN/Inf`)
  after ~20 min instead of converging. Also runs in CI as a manual-only
  `workflow_dispatch` (`.github/workflows/slow-tests.yml`), one job per
  target in parallel via `pymeep` installed from conda-forge — not on every
  push (too slow for that), and not pip-installed (see High-level
  recommendation 6: `pyproject.toml`'s `[fdtd]`/`[all]` extras currently
  declare the wrong `meep` package).
- **Not covered by either suite**: `FDTDSolver`/`EMESolver` used standalone
  outside `HybridSolver` (the golden tests only exercise them via `solve()`/
  `solve_section()` on the four targets above), `mmi_2x2_*`,
  `modal_multiplexer_hybrid` (deliberately excluded — its numerical accuracy
  isn't yet trusted enough to lock in as a golden reference), `ParallelSweep`,
  and all of `feces/plot/`.

## Recommendations

Sorted into High-level / Solver / Geometry / Plotting / Examples / Misc. Ordered
within each section by (impact) × (leverage — how much is genuinely missing vs. just
"meep already does this, expose it"), not by raw novelty. Sections are kept even
where a category has fewer or lower-priority items.

### High-level

1. **`mmi_2x2` golden test.** `examples/mmi_2x2/` (EME, FDTD, hybrid — see
   Examples recommendations) has no golden-value regression test yet, unlike
   its `mmi_1x2` counterparts. Straightforward to add once the example itself
   is validated (same pattern as the existing 4 targets), but deliberately
   deferred until then — no value locking in a reference before trusting the
   example's own numbers.
2. **`modal_multiplexer_hybrid` numerical accuracy.** Excluded from the golden
   suite on the same basis as item 1, but blocked on more than just "add a
   test": the device's accuracy isn't yet trusted, so this needs an
   independent correctness check (e.g. against a known-good reference, or at
   minimum a `check_passive()`-style sanity pass plus manual inspection)
   *before* a golden test would mean anything — a regression test only
   guards against drifting away from a value that's already right.
3. **Documentation.** No docs site, no changelog beyond this `README.md`/
   `notes.md` pair — the `examples/` scripts are the only tutorials that
   exist. At minimum a `CHANGELOG.md` tracking the still-static
   `__version__ = "0.1.0"`.
4. **Type checking in CI.** Type hints are used throughout (`from __future__ import
   annotations`) but nothing enforces them — no `mypy`/`pyright` config, no `py.typed`
   marker. Cheap to add now that CI exists (see the Testing library reference above).
5. **Declare transitive dependencies explicitly.** `gdsfactory`/`tidy3d` are consumed
   indirectly through `meow-sim` (e.g. `meow.compute_modes_tidy3d`) but aren't
   declared in `pyproject.toml` — a `meow-sim` version bump that drops/renames those
   could silently break FECES. Pin or declare them directly, or guard the import
   with a clear error message.
6. **Fix `pyproject.toml`'s `[fdtd]`/`[all]` extras — they install the wrong
   package.** Found while wiring up slow-suite CI: `fdtd = ["meep"]` and
   `all = ["meow-sim", "meep"]` declare `meep` as a pip dependency, but PyPI's
   `meep` package is an unrelated project (a lightweight task runner) — `pip
   install feces[fdtd]` succeeds with no error and leaves `import meep`
   completely broken. The real MIT electromagnetics package only ships via
   conda-forge (`conda install -c conda-forge pymeep`; confirmed working
   both locally and in `.github/workflows/slow-tests.yml`'s CI job). Fix: drop
   `meep` from these extras (or
   replace with a clear error message on import failure pointing at the
   conda-forge install command) and document the conda-forge requirement in
   `README.md`'s Install section, not just buried in the Testing section.
7. **Standardize on `logging`, not `print()` — done.** The 5 unconditional
   `print()` calls in `feces/solvers/` (`HybridSolver.solve()`'s 4
   per-section progress lines, `EMESolver.solve_section_length_sweep`'s
   mid-sweep neff diagnostic) are now `logger.info()`/`logger.debug()` via
   `logging.getLogger(__name__)`, matching `ParallelSweep`'s existing
   convention — silent by default for any caller that hasn't configured a
   handler, instead of always writing to stdout. Examples' own `print()`s
   are intentionally untouched (user-facing tutorial output, not library
   diagnostics).

### Solver

1. **Lossy/dispersive material support.** `MaterialStack.si_medium(wl)`/
   `sio2_medium(wl)` (`materials/stack.py`) build an `mp.Medium(index=...)` from
   whichever *tabulated* wavelength in meow's built-in Si/SiO2 data is nearest to
   the requested one (`np.argmin`) — a nearest-neighbor lookup, not interpolation,
   and always real-valued, so there's no true dispersion curve across a wavelength
   sweep and no absorption/loss. Wiring meep's actual Lorentzian-susceptibility
   fitting (or reusing `meep.materials` directly) into `MaterialStack` would enable
   real loss/broadband dispersion and open the door to SiN, III-V, or lossy/doped
   silicon.
2. **Symmetry exploitation in `FDTDSolver`.** `MMI`/`DirectionalCoupler` are
   frequently mirror-symmetric about y=0; `solve_section()` currently always draws
   and simulates the full cross-section. An opt-in `mp.Mirror` symmetry when the
   geometry/excitation is actually symmetric could roughly halve FDTD cost for free —
   complements the reciprocity-based min-simulation-count optimization already there.
3. **Circuit-level S-matrix composition.** `SMatrix` has no way to combine several
   devices' results into a larger circuit (e.g. cascading an MMI's `SMatrix` with a
   directional coupler's). `HybridSolver` already does exactly this kind of
   matrix-chain cascading *within* one device across sections — generalizing that
   same machinery to cascade *across* independently-solved devices (or documenting a
   recipe for handing `SMatrix` results to SAX or Simphony, both established OSS
   tools for exactly this layer, SAX JAX/netlist-based and Simphony
   component-library/plugin-based) would stop device-level sims from being a dead
   end. (Note: `CompoundGeometry` — see Geometry section — solves a related but
   distinct problem at the geometry layer: it lets `HybridSolver` simulate a device
   built from several `DeviceGeometry` parts as one cascade *before* solving, but
   doesn't let two already-*solved* `SMatrix` results be combined after the fact,
   which is what this item is about.)
4. **`ParallelSweep` + `HybridSolver` integration.** The sweep utility only calls
   whole-device `solver.solve()`; it doesn't know about `HybridSolver.solve()`'s
   `cache=` mechanism, so a parallel sweep over a hybrid device would redundantly
   re-simulate cached sections across worker processes.
5. **Automatic mesh-resolution convergence.** Beyond the manual "run a list of
   resolutions and eyeball it" pattern (see Misc/Examples below) — meep/tidy3d both
   have some form of resolution auto-selection; FECES has none, always a fixed
   user-chosen `res`.
6. **Bigger lift, only if scope grows toward design automation:** adjoint/
   inverse-design support, periodic/grating structures via an RCWA backend
   (S4-style) for problems where FDTD/EME are the wrong tool, near-to-far-field
   transform. For adjoint specifically, `meep.adjoint` (already what FECES's own
   `FDTDSolver` sits on top of) is the lowest-friction path — it's autograd/JAX-based
   and gradient computation wraps around any optimization library, so no new solver
   dependency would be needed, just a FECES-side API for defining a figure of merit
   over a `HybridSolver` chain's `SMatrix` output. Both tidy3d (moving from a
   JAX-adjoint plugin to native `autograd` support) and the newer JAX-native FDTDX
   project (automatic differentiation built into the solver via time-reversibility,
   not bolted on) confirm autodiff-native is now the field's default direction, not
   a niche add-on — worth knowing if this ever gets prioritized. Substantial,
   well-trodden ground already owned by meep/tidy3d/S4/FDTDX respectively — worth it
   only if a concrete device need shows up that FDTD/EME genuinely can't address.
7. **Mirror-section reuse in `HybridSolver` (reciprocity across sections, not just
   within one) — done.** `HybridSolver.__init__` now takes `mirror_pairs={k: j}`
   (plus `mirror_axis`, `mirror_check_samples`/`mirror_check_tol`): section `k`'s
   transfer matrix is reused as `mats[j].transpose(0,2,1)` instead of a real
   `solve_section()` call. Declared at construction time (like meep's `mp.Mirror`
   is declared once on `mp.Simulation`, not per-run) since it's a fact about the
   device's section topology — but the geometric guardrail
   (`width_at`/`center_at` sampled across both sections' spans under x-reversal,
   raising `ValueError` on a mismatch) re-runs every `solve()` call, because
   `geometry_kwargs` (`input_idx`/`output_idx`) can change which physical port is
   actually live between calls, unlike meep's truly per-run-structural symmetry.
   Construction-time checks reject out-of-range/self/forward-referencing/chained
   pairings and mode-count mismatches up front. A mirror-reused section still
   populates the shared `cache` dict under its own key. Tested in
   `tests/test_hybrid.py` (`TestMirrorPairs`/`TestMirrorPairsConstruction`, mock
   solvers, no meep/meow).
   Length-sweep support landed too: `HybridSolver.solve_length_sweep(geometry,
   wavelengths, lengths, sweep_section, ...)` sweeps one uniform section's length
   analytically (delegating to that solver's `solve_section_length_sweep()`,
   which EMESolver implements — modes solved once per wavelength, lengths applied
   as propagation phases) and returns the raw `[n_L, M, out, in]` cascade grid.
   Non-swept sections keep byte-identical cache keys to `solve()`, so entries
   dedupe freely across both methods and across per-port calls. Tested in
   `tests/test_hybrid.py` (`TestLengthSweep`).
   **Still open — and now the single remaining blocker**: neither hybrid example
   actually calls the generic API anymore. `mmi_2x2_hybrid.py` was migrated onto
   `directional_coupler_hybrid.py`'s hand-rolled pattern (solve the input section
   once per input arm, derive every output-side matrix by mirror transpose, one
   manual einsum cascade — zero dedicated output solves, same 2 FDTD + 1 EME
   budget, full 2×2 excitation family available from the same matrices). That
   through+cross-from-one-excitation-family composition pairs one arm's
   transposed matrix with the *other* arm's excitation in a single cascade,
   which `solve()`/`solve_length_sweep()`'s one-call contract can't express even
   with `mirror_pairs` — a multi-excitation-config `solve()` extension (declare
   several excitation kwargs up front, get every port-combination coefficient
   out of one call) would let both examples migrate off the manual einsum
   composition. `mirror_pairs`/`solve_length_sweep` remain fully supported and
   tested; they're just not exercised by an example at the moment.
8. **Wavelength-interpolation scheme for mode-solve/EME-driven sections.** A
   real MPB mode solve is expensive per wavelength, while FDTD DFT monitors
   sample a dense wavelength grid nearly for free (one broadband pulse, many
   frequencies in one run). Resonant devices need thousands of wavelength
   points to resolve a linewidth that's a small fraction of the sweep span —
   `examples/ring_resonator/ring_resonator_hybrid.py` hits this directly: an
   FSR of ~13 nm across a 100 nm C-band sweep needs a ~2000-point grid to
   plot resonance dips cleanly, but only solves `neff(lambda)` at 13 sparse
   points (each a real mode solve) and fits a degree-2 polynomial before
   evaluating on the dense grid — piecewise-linear interpolation would put a
   kink in `d(phi)/d(lambda)` at every solve node and corrupt resonance
   positions, so the fit has to be smooth in `neff` itself, not in phase
   (which winds too fast to interpolate directly — that example's FDTD-
   coefficient interpolation is magnitude + `np.unwrap`'d phase for the same
   reason). Done ad hoc, entirely inside that example script, because
   nothing in `feces/solvers/eme.py` supports it. A library-level version
   could look like sparse-solve/dense-eval wavelength arguments on
   `EMESolver.solve_section()` (solve modes at a caller-chosen sparse subset,
   fit, evaluate the returned transfer matrix at every wavelength in a
   separate dense array), or a standalone reusable `NeffTable`/mode-fit
   utility both `EMESolver` and examples could share. Related to, but
   distinct from, Solver recommendation #1 (`MaterialStack`'s nearest-
   tabulated-wavelength index lookup) — that's a materials-dispersion gap;
   this is a mode-solve-cost gap, and the two would compound in a real
   dispersive-material implementation (interpolating dispersive neff at
   already-sparse solve points).
9. **Persistent save-folder cache for `solve_section()`/`solve()` — done.** New
   `feces/cache/` module: `build_run_config(solver, method_name, call_kwargs,
   geometry)` fingerprints a call into a JSON-safe dict — solver params
   (`vars(solver)`, `material_stack` serialized explicitly), the call's own
   args, and the exact section geometry via `geometry.section_gds(x0, x1,
   ...)` (plain vertex data, no meep/meow objects, so it's the precise shape
   solved rather than the whole-device layout). `lookup()`/`save_run()`
   (`feces/cache/store.py`) hash that dict into a `save_dir` subfolder named
   by `sha256(...)[:16]`; a hash hit still gets a full `config.json`
   equality check before being trusted, guarding against truncated-hash
   collisions or config-schema drift across code versions. Each subfolder
   holds `config.json` (human-readable, for manual comparison), `result.pkl`,
   and `geometry.png`. Wired into `FDTDSolver.solve_section()`,
   `EMESolver.solve_section()`/`solve_section_length_sweep()` via a new
   `save_dir: Optional[str] = None` parameter (falsy = disabled, matching
   `HybridSolver.solve(cache=None)`/`ParallelSweep(output_dir=None)`), and
   threaded through `HybridSolver.solve()`/`solve_length_sweep()` — forwarded
   to each section's `solve_section()` call only when that solver's own
   signature accepts `save_dir` (`_supports_save_dir()`, so a custom solver
   type without the parameter still works, per the module's minimal
   `solve_section()` contract). Separate mechanism from `HybridSolver`'s
   in-memory `cache` dict: `cache` dedupes within/across `solve()` calls in
   one process; `save_dir` persists to disk across process runs; the two
   compose freely. Verified end-to-end (EME hit/miss/negative-miss on a
   `res` change, FDTD hit/miss, a 2-section `HybridSolver` chain) — cache
   hits return bit-identical results in <1ms vs several seconds for a real
   solve.
   **Extended to the plain whole-device `solve()` methods too.** The
   original design was section-only: `build_run_config()` requires `x0`/`x1`
   and fingerprints geometry via `geometry.section_gds(x0, x1, ...)`, which
   `FDTDSolver.solve()`/`EMESolver.solve()` don't have (no section bounds —
   they build geometry via `geometry.polygons(stub_length=...)` and
   `geometry.to_gds_eme(...)` respectively, the same builders `to_meep_prisms()`
   already calls). A new `build_whole_device_config()` (`feces/cache/config.py`)
   takes an already-built polygon dict instead of computing one via
   `section_gds()`, sharing the same JSON-safe fingerprint assembly (factored
   into a private `_fingerprint()` helper both builders call). `FDTDSolver.solve()`
   now accepts `save_dir` directly; `EMESolver.solve()` caches per output
   port/side (its existing `for side in ports:` loop is the natural unit —
   matches `solve_section()`'s "one call = one geometry" granularity), so
   `out_top`/`out_bot` calls with different `output_idx` kwargs land in
   separate cache entries rather than colliding. All 10 examples that call a
   `save_dir`-capable method now pass `save_dir=os.path.abspath(os.path.join(
   _HERE, "..", "example_sim_cache"))` (a shared cache dir at
   `examples/example_sim_cache/`, sibling to each device's own example
   folder) — the two exceptions, `directional_coupler_mode.py` and
   `ring_resonator_hybrid.py`, use raw meep/meow with no FECES solver call to
   attach it to. The two examples that already used the in-memory
   `HybridSolver` `cache={}` dict (`mmi_1x2_hybrid.py`,
   `modal_multiplexer_hybrid.py`) now pass both `cache=` and `save_dir=`
   together, since the two mechanisms compose. `EMESolver.solve_sweep()` (the
   batched-wavelength whole-device path) does **not** support `save_dir` yet
   — no example calls it, so it was out of scope here, but it's a real gap
   now that `solve()` has the parameter and `solve_sweep()` doesn't.
   **Known risk, not fully closed**: config completeness. Any
   solver/geometry attribute that affects the physical result but isn't
   scalar (so falls outside `vars(solver)`'s generic filter) is excluded
   from the fingerprint with a `logger.warning`, not an error — a future
   solver attribute added without updating `build_run_config`/
   `build_whole_device_config` could produce a silent false-positive cache
   hit (wrong physics, no error). No regression test guards against this
   yet, for either fingerprint path.

### Geometry

1. **Ring/racetrack resonator.** The single most common photonics building block
   that isn't in the device library at all. Natural fit for the existing
   `DeviceGeometry` interface (bus waveguide + coupled ring, `section_prisms`/
   `width_at`/`center_at` following the `DirectionalCoupler` pattern), and directly
   exercises `HybridSolver`'s generic cascading (bus straights as FDTD/EME + coupling
   region). **Partially addressed** — `examples/ring_resonator/ring_resonator_hybrid.py`
   proves out the physics (closed-form all-pass feedback assembly, not a
   `HybridSolver` chain — a ring is a loop, not a cascade), but does so entirely in
   raw meep/meow, bypassing the geometry layer and every solver class. What's still
   missing for a first-class version: (a) an actual `RingResonator`/racetrack
   `DeviceGeometry`; (b) y-normal (x-z plane) source/monitor support in
   `FDTDSolver` — it only ever builds `mp.Vector3(0, w, h)` y-z planes today, and
   the two bus-to-ring coupling ports in a racetrack coupler section are
   perpendicular to that; (c) a feedback-loop assembly mode in `HybridSolver`
   (or a sibling), since its `solve()` is fundamentally a one-pass matrix chain,
   not a resonator sum.
2. **GDS import (read a real/foundry/gdsfactory layout) and export (tape-out-ready
   output).** The gdsfactory dependency was deliberately removed from the internal
   pipeline (raw polygon dicts are faster/simpler for FECES's own parametric shapes),
   but that leaves no way to simulate an arbitrary externally-supplied layout or
   export a FECES geometry as a real `.gds`. A thin, optional adapter (gdsfactory
   `Component` ↔ FECES polygon dict) would restore that interop without
   reintroducing gdsfactory as a hard dependency. gdsfactory's PDK ecosystem keeps
   growing (gdsfactory+ now lists 43+ foundry PDKs, e.g. a 2026 Lightwave Logic
   high-speed-modulator PDK for GlobalFoundries silicon photonics) — every PDK added
   there is a real device FECES currently has no way to even load, let alone
   simulate.
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

1. **Stop duplicating solver position/size math in the plot functions — done.**
   `_draw_section_fdtd_markers` (plot/geometry.py) used to independently
   recompute where `FDTDSolver.solve_section()` places its source/monitor
   apertures, and an earlier bug-hunt found two real mismatches there
   (marker offset by `monitor_margin` when the solver places it flush at
   the section edge; marker always drawn full-height when the solver only
   does that for `n_measure > 1`). Both now call one shared function,
   `section_fdtd_layout()` (`feces/solvers/fdtd.py`, module-level, no meep
   import needed): given the section bounds/mode counts/solver settings, it
   returns which side excites vs. measures and every source/monitor
   position and y-aperture, in the section's local frame (0 at
   `x0-stub_left`, `L` at `x1+stub_right`). `solve_section()` shifts those
   into its own meep-centered frame; the plot draws them directly, clamped
   to the non-PML domain for display only (the clamp is cosmetic — it never
   feeds back into what's actually simulated). Refactoring surfaced one
   more latent mismatch of the same kind: the old plot special-cased an
   on-axis (`center_at()==0`) `half_domain_aperture` port to draw a
   full-height marker, but `solve_section()` itself always halves that
   port's aperture regardless of axis position — the plot now matches the
   solver by construction instead of by two independently-maintained
   copies agreeing (or not) with each other. Verified against real meep:
   both excitation directions, `half_domain_aperture`, and the
   `monitor_center`/`monitor_size` override path all reproduce their
   pre-refactor numbers.
2. **FDTD time-domain field/movie export.** `plot_eme_fields` only covers EME (static
   `|component|(x,z)` heatmaps via `mw.propagate_modes`); there's no equivalent for
   FDTD at all — no meep `Animate2D`-style field snapshot or movie export anywhere in
   the package.
3. **Cascaded full-device field plots for `HybridSolver`.** `plot_eme_fields` only
   ever excites one pure input mode via `mw.propagate_modes(..., excite_mode_l=exc,
   ...)` (a single int, not a superposition), and there's no FDTD field-array export
   at all (see item 2 above) — so today there's no way to visualize the actual field
   distribution through a whole hybrid-simmed device, section by section, the way
   `plot_hybrid_geometry` already stitches the *geometry*. Feasible — meep's
   `add_dft_fields()`/`get_dft_array()` gives the full field array over an FDTD
   section's cell at each DFT frequency, in the *same* sims `solve_section()`
   already runs (one per excited basis mode, unit amplitude) — no extra sims needed
   for extraction itself. The real work: (a) `HybridSolver.solve()`'s cascade loop
   computes the running partial matrix-product at each boundary but only returns
   the final end-to-end result — getting per-section actual excitation amplitudes
   needs that intermediate state exposed (small API change: return/callback the
   per-boundary vectors, or reconstruct from the per-section matrices the loop
   already builds); (b) the true field is a **superposition** of the unit-basis-mode
   field solutions weighted by that amplitude vector (valid by linearity — Σ_k a_k ·
   field_for_mode_k) — needed on both the FDTD side (weight the per-mode DFT arrays)
   and the EME side (`propagate_modes` only takes one mode index, so multi-mode
   interfaces need the same weighted-sum trick, one call per basis mode);
   (c) stitching across sections into shared device-frame coordinates, mirroring
   `plot_hybrid_geometry`'s existing polygon-stitching approach, with FDTD's
   meep-resolution grid and EME's `z=np.linspace(...,500)` grid at different pixel
   densities per panel — a multi-panel composite with consistent color-scale
   normalization, not one seamless array. Cost-bounded only if scoped as opt-in,
   single-wavelength (never bolted onto the routine full-sweep S-matrix solve
   path) — DFT field arrays over a whole 3-D cell × many wavelengths × up to ~15
   excited modes is a lot of data. No fundamental blocker (unlike the y-normal-
   monitor gap for bends/rings) — squarely within what meep+meow already expose,
   just a real multi-piece feature: `HybridSolver` API addition + opt-in
   `FDTDSolver.solve_section()` field-extraction mode + superposition-aware EME
   field call + a stitching/plotting layer.
4. **Phase / group-delay plots.** `SMatrix` stores full complex coefficients, but
   `plot_performance`/`plot_splitting_ratio` only ever plot `|coeff|^2` (linear T) and
   dB insertion loss — phase and group delay (from phase vs. wavelength) are never
   visualized despite being available in the data.
5. **A shared convergence plot.** Complements the Misc/Examples resolution-sweep
   item: a `plot_convergence(results_by_resolution, metric)` helper that overlays a
   chosen quantity (e.g. L_c, T) across resolutions on one figure, instead of
   eyeballing separate per-resolution plots side by side (what
   `directional_coupler_mode.py`/`.ipynb` does today).
6. **Interactive/3D visualization.** Everything is a static matplotlib PNG; no
   plotly/bokeh option for geometry or field plots, unlike tidy3d's interactive
   viewer.

### Examples

1. **2×2 MMI example — done; 1×4 still missing.** `examples/mmi_2x2/` (EME, FDTD,
   hybrid variants) exercises a genuine multi-input device — the first consumer
   of `MMI`'s new `input_idx`/`input_indices` support (see `MMI` in the library
   reference above). No 1×4 (or other multi-output-only) example yet, though
   nothing in the geometry blocks one.
2. **No example uses `ParallelSweep`.** The library ships a generic parallel
   parameter-sweep utility that zero examples demonstrate — a user has no worked
   reference for the one utility meant to replace ad hoc per-device sweep scripts.
3. **No example demonstrates mesh-convergence checking as a reusable pattern** beyond
   the one-off resolution loop hand-built into `directional_coupler_mode.py`/`.ipynb`
   (a pure meow/MPB even/odd-supermode script that doesn't use any FECES solver
   class at all) — worth generalizing into its own example once the Misc convergence
   utility exists.
4. **Ring resonator example — done, but as a raw-meep script, not paired with a
   `RingResonator` `DeviceGeometry`.** `examples/ring_resonator/ring_resonator_hybrid.py`
   exists (see Geometry recommendation #1) and does call `SMatrix.check_passive()`
   as standard practice, but doesn't use the convergence-utility recommendation
   above (no resolution sweep) and can't use `plot_fdtd_geometry`/`plot_hybrid_geometry`
   at all (both require a `DeviceGeometry`) — it hand-draws its own geometry
   preview plots instead. The geometry-class-backed version (proper
   `plot_*_geometry` integration, resolution-sweep convergence check) remains
   open.

### Misc

1. **Mesh-convergence-loop utility.** The pattern hand-rolled in
   `directional_coupler_mode.py`/`.ipynb` (loop over a resolution list, save/plot
   each into its own subfolder) should be a reusable `feces` utility (e.g.
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
   `make_output_folders`, but usage is inconsistent across examples (most
   `examples/*.py` scripts hand-roll their own `os.makedirs`/`output_dir` logic
   instead) — worth consolidating on the existing utility everywhere instead of ad
   hoc per-script path-building.

## Not recommending

- Building FECES's own raw FDTD/EME numerical core — meep and MPB/meow already do
  this well; FECES's value is the orchestration layer on top, not the field solver.
- Matching tidy3d's cloud/GPU infrastructure — out of scope for what is currently a
  single-machine/HPC-Makefile-driven library.

## To-Do

Condensed, actionable version of the Recommendations above — done items dropped,
partially-done items trimmed to just the remaining work. Same section order/order
within section; see Recommendations for full rationale on any item.

### High-level

1. `mmi_2x2` golden test — add once the example's own numbers are trusted.
2. `modal_multiplexer_hybrid` numerical accuracy — needs an independent correctness
   check before a golden test would mean anything.
3. Documentation — docs site / `CHANGELOG.md` (version still static `0.1.0`).
4. Type checking in CI — `mypy`/`pyright` config, `py.typed` marker.
5. Declare transitive dependencies explicitly — `gdsfactory`/`tidy3d` pulled in via
   `meow-sim`, undeclared in `pyproject.toml`.
6. Fix `pyproject.toml`'s `[fdtd]`/`[all]` extras — they install the wrong `meep`
   package (PyPI name collision); real one is conda-forge `pymeep`.

### Solver

1. Automatic mesh-resolution convergence — no auto-selection, user picks `res` by
   hand.
2. Bigger lift, only if scope grows toward design automation: adjoint/inverse-design
   support, RCWA/grating structures, near-to-far-field transform.
3. Wavelength-interpolation scheme for mode-solve/EME sections — sparse-solve/
   dense-eval API (see `ring_resonator_hybrid.py`'s ad hoc version for the pattern).
4. Extend `save_dir` caching to `EMESolver.solve_sweep()` — the batched-wavelength
   whole-device path is the one remaining `solve()`-family method without it.

### Geometry

1. Ring/racetrack resonator, first-class version — needs a `RingResonator`
   `DeviceGeometry`, y-normal monitor support in `FDTDSolver`, and a feedback-loop
   assembly mode in `HybridSolver` (see `ring_resonator_hybrid.py` for the raw-meep
   proof of physics).
2. GDS import/export — thin gdsfactory `Component` <-> FECES polygon-dict adapter.
3. Y-branch/adiabatic splitter and waveguide crossing geometries.
4. Rebuild `BendGeometry` on top of `SBend` (or a new arc primitive) — still a
   separate, non-`DeviceGeometry`-integrated legacy 90° bend.
5. `Taper`: polynomial/Euler-spiral width profile; standalone `Taper`/`SBend`
   example.

### Plotting

1. FDTD time-domain field/movie export (EME-only today).
2. Cascaded full-device field plots for `HybridSolver` — stitch per-section FDTD
   DFT field arrays + EME `propagate_modes` output, mode-superposition-weighted
   by each section's actual cascade amplitude (not unit excitation); needs a
   small `HybridSolver` API addition to expose per-boundary amplitudes.
   Feasible, no blocker — moderate-to-substantial, opt-in/single-wavelength only.
3. Phase / group-delay plots — `SMatrix` has the data, nothing plots it.
4. A shared `plot_convergence()` helper.

### Examples

1. 1×4 (or other multi-output-only) MMI example.
2. `ParallelSweep` example.
3. Mesh-convergence-checking example.
4. Ring resonator example paired with a real `RingResonator` `DeviceGeometry`
   (proper `plot_*_geometry` integration, resolution-sweep convergence check).

### Misc

1. Mesh-convergence-loop utility (`sweep_resolution(run_fn, resolutions)`-style).
2. Result serialization beyond raw pickle — Touchstone `.sNp` or JSON/HDF5.
3. Declarative simulation runner — YAML/JSON device+solver+sweep spec.
4. Output-folder conventions — consolidate examples on `resolve_output_folder`/
   `make_output_folders` everywhere.

### To-do (later)

Deferred, not dropped — lower near-term priority than the sections above.

1. Lossy/dispersive material support (Solver) — `MaterialStack` does nearest-
   tabulated-wavelength real-index lookup only, no interpolation/loss.
2. Symmetry exploitation in `FDTDSolver` (Solver) — opt-in `mp.Mirror` for
   symmetric devices.
3. Circuit-level S-matrix composition (Solver) — no way to cascade
   independently-solved `SMatrix` results into a larger circuit.
4. `ParallelSweep` + `HybridSolver` integration (Solver) — sweep doesn't use
   `HybridSolver`'s `cache=` mechanism.
5. Interactive/3D visualization (Plotting) — matplotlib-only today.
