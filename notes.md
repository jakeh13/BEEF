# FECES capability gap analysis vs. meep / meow / tidy3d / other open-source EM tools

Research memo — comparison of FECES against meep, meow, tidy3d, and the wider open-source
photonics/EM landscape, plus a prioritized list of open items. Everything under
"Recommendations" is unimplemented unless the item text says otherwise; "Current library
reference" is the descriptive record of what's actually in the tree today. Last refreshed
2026-07-20 against the current `feces/` source and a fresh pass over the OSS landscape (see
sources inline in the table).

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
| Time-domain FDTD | 2D or 3D via `FDTDSolver(three_d=...)`, single fixed uniform `res` | Native, mature, adaptive subpixel smoothing | — | Cloud GPU FDTD, auto mesh refinement | openEMS (general 3D FDTD, antenna-focused); **FDTDX** (JAX-based, GPU-accelerated, native autodiff via time-reversibility, published in JOSS 2026 — 10-415x faster than traditional FDTD per its own benchmarks); Ceviche/EMopt (earlier autograd-based FDTD); GSvit (GPU-accelerated) |
| Eigenmode expansion (EME) | `EMESolver` wraps meow's cross-section mode-match + `compute_s_matrix` | — | Native (this *is* meow's job) | Has an EME plugin (steady-state frequency-domain alternative to its FDTD solver) | EMEpy (pure-Python, neural-net-accelerated mode profiles); CAMFR (older EME framework) |
| Hybrid FDTD+EME chaining | **Yes — `HybridSolver`, FECES's actual differentiator** | No | No | No (closest: manually stitching, no generic N-section cascade) | Not found in a fresh landscape pass (`awesome_photonics`, direct searches) — closest is `simphox` (Stanford Solgaard Lab, bundles FDTD/beam-propagation/circuit solvers), but it isn't documented as a generic reciprocity-based N-section cascade the way `HybridSolver` is |
| Dispersive materials | Nearest-tabulated-wavelength real-index lookup (`np.argmin` on meow's built-in Si/SiO2 tables) — no interpolation, no loss | Full Lorentzian/Drude/Debye susceptibility fits | Has real material models (used via meow calls) | Automatic dispersion-fitting plugin | — |
| Nonlinear materials (χ2/χ3) | None | Yes (Kerr, χ2) | No | Yes | — |
| Lossy/absorptive materials | None (imaginary index discarded) | Yes | Partial | Yes | — |
| Anisotropic / magneto-optic | None | Yes | No | Yes | — |
| Material library breadth | 2 stacks: `SOI_160`, `SOI_220` (Si/SiO2 only, thickness-parametrized) | `meep.materials` — dozens of real fits (SiN, III-V, metals, oxides) | Some via meow | Large curated library | — |
| Periodic structures / gratings / PhC | None (no Bloch BCs anywhere) | Yes (Bloch-periodic boundaries) | No | Yes | RCWA is the standard faster-than-FDTD tool for this: **S4** (long-standing C++/Python reference implementation) plus a newer wave of autodiff-native Python rewrites — **Meent**, **FMMAX**, **grcwa**, **nannos**, **inkstone**; **MPB** for full PhC band structure (transitively available via meow but not exposed) |
| Symmetry exploitation | None — always draws/simulates full both-arm geometry even for mirror-symmetric devices | `mp.Mirror`/`mp.Rotate2` — 2-4x speedup for symmetric structures | n/a | Yes | — |
| Adjoint / inverse design / topology optimization | None | Adjoint solver (`meep.adjoint`, autograd/JAX-based) | No | Native `autograd` support (the older JAX-based adjoint plugin has been deprecated in favor of this) | **FDTDX** — JAX-native FDTD with autodiff built into the solver itself via time-reversibility, not bolted on afterward; **Ceviche**/**EMopt** — earlier autograd-based FDTD/inverse-design frameworks in the same lineage |
| Near-to-far-field transform | None | `add_near2far` | No | Yes | — |
| Mesh auto-refinement / convergence checking | None — manual `res`, user must sweep by hand | Manual (subpixel smoothing helps, not auto-refine) | n/a | Automatic adaptive meshing | — |
| Device geometry library | `MMI` (N×M, per-port `MMIPort` taper size/placement), `DirectionalCoupler`, standalone `Taper`/`SBend` building blocks, `CompoundGeometry`/`WaveguideRun` (compose several devices into one cascade), `BendGeometry` (legacy, poorly integrated) | n/a (raw geometry primitives only) | n/a | Has some parametric primitives + PDK support | **gdsfactory** — the dominant OSS parametric photonic component library (3M+ downloads; gdsfactory+ lists 43+ foundry PDKs as of 2026); **picwriter**, **phidl** (lighter-weight alternatives) |
| GDS import/export, PDK interop | **None** — gdsfactory dependency was actively removed (commit `ffcbf5a`); geometries are FECES-internal polygon dicts only, can't simulate a real fab layout or export one for tape-out | n/a | Uses gdsfactory `Component`s | Has gdsfactory integration | gdsfactory itself |
| Ring/racetrack resonator | Proof-of-physics only (`examples/ring_resonator/`, raw meep/meow, see Geometry recommendation #1) — no `DeviceGeometry` yet | n/a | n/a | n/a | gdsfactory has one |
| Grating coupler / edge coupler geometry | Missing | n/a | n/a | n/a | gdsfactory |
| Circuit-level composition (cascade many devices' S-matrices into one chip netlist) | **None** — `SMatrix` has no `+`/`@`/composition operators | n/a | n/a | n/a | **SAX** (JAX/netlist-driven, differentiable) and **Simphony** (component-library/plugin-based, itself built on SAX) are the two standard choices; **opics**, **lekkersim**, **photontorch** (adds time-domain circuit analysis) round out the space |
| Result serialization | Pickle only, no schema/versioning | n/a | n/a | Its own `.hdf5`/JSON | Touchstone `.sNp` is the RF/photonics interop standard, not supported |
| Passivity/unitarity sanity check | `SMatrix.check_passive()` — implemented, called after every `solve()`/`hybrid.solve()` in every example and golden test as standard practice | n/a | n/a | n/a | — |
| Parallel parameter sweeps | `ParallelSweep` — local `multiprocessing.Pool` grid sweep, but only calls whole-device `solver.solve()`, **not** `HybridSolver`-aware, no cloud/cluster backend beyond an external SLURM Makefile | n/a (user-driven) | n/a | Cloud batch is tidy3d's core value prop | — |
| Field visualization | Static EME field heatmaps only (`plot_eme_fields`); **no FDTD time-domain field/movie export at all** | `mp.Animate2D`, live field plotting | n/a | Interactive Plotly 3D viewer | — |
| Interactive/3D visualization | None (matplotlib static PNGs only) | Some (via `mp.plot2D`) | n/a | Yes, core feature | — |
| Testing | 77 fast unit tests (`tests/`, meep/meow-free) run in CI on every push; 4 slow golden-value integration tests (real reduced-res meep/meow sims) run in CI as a manual-only `workflow_dispatch` matrix via conda-forge `pymeep` (`pytest --runslow` locally); no lint/type-checking yet | n/a | n/a | n/a | — |

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
  self-contained whole-device simulation, and (like `solve_section()`) accepts an
  optional `save_dir` for the persistent on-disk cache — see below.
  `solve_section()` is the piece
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
  Only ever builds x-normal (y-z plane) sources/monitors — there's no y-normal
  (x-z plane) support, which is what blocks a first-class ring/racetrack resonator
  (see Geometry recommendation #1).
- **`EMESolver`** — meow-based eigenmode expansion. `num_modes`/`num_cells`/`lengths`
  (parallel lists, one entry per region left-to-right, e.g. `[wg, taper, MMI body,
  taper, wg]`) plus `oxide_width`/`oxide_height`/`res`/`material_stack` configure the
  cross-section mesh; `mode_solver` picks `"meep"` (matches FDTDSolver's own basis —
  use this in hybrid chains to avoid basis mismatches at FDTD/EME interfaces) or
  `"tidy3d"` (fine for an EME-only solve with no FDTD sections to stay consistent
  with). `solve()`/`solve_sweep()` run a complete whole-device eigenmode expansion
  (the latter batches wavelengths for memory efficiency). `solve()` (unlike
  `solve_sweep()`, see Solver recommendation #1) also accepts an optional
  `save_dir`, cached per output port/side since each side builds its own
  geometry/cells. `solve_section()` /
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
  `solve()`'s exactly so entries dedupe across both methods. Worked examples of
  both features live in `tests/test_hybrid.py` (`TestMirrorPairs`/
  `TestLengthSweep`) — the hybrid examples themselves (`mmi_2x2_hybrid.py`,
  `directional_coupler_hybrid.py`) still hand-roll the same reuse manually because
  their cross/through paths mix two excitations' matrices in one cascade, which the
  one-call `solve()` contract can't express even with `mirror_pairs` (see Solver
  recommendation #2).
- **Persistent on-disk cache** (`feces/cache/`) — every solve-family method
  (`FDTDSolver.solve()`/`solve_section()`, `EMESolver.solve()`/`solve_section()`/
  `solve_section_length_sweep()`, `HybridSolver.solve()`/`solve_length_sweep()`)
  takes an optional `save_dir: Optional[str] = None` (falsy = disabled, no behavior
  change if unset). `build_run_config()`/`build_whole_device_config()`
  (`feces/cache/config.py`) fingerprint a call into a JSON-safe dict — solver
  params (`vars(solver)`, `material_stack` serialized explicitly), the call's own
  kwargs, and the exact geometry drawn (`geometry.section_gds()`/`.polygons()`/
  `.to_gds_eme()`, no meep/meow objects involved) — hashed (`sha256(...)[:16]`)
  into a `save_dir` subfolder holding `config.json` (human-readable), `result.pkl`,
  and `geometry.png`; a hash hit still gets a full `config.json` equality check
  before being trusted, guarding against truncated-hash collisions or
  config-schema drift across code versions. This is a separate mechanism from
  `HybridSolver`'s in-memory `cache={}` dict (dedupes within/across `solve()` calls
  in one process; `save_dir` persists to disk across process runs) — the two
  compose freely. All 12 example `.py`/`.ipynb` pairs across every device folder
  point at one shared `examples/example_sim_cache/` directory; a 2026-07 audit
  diffed every cache-fingerprinting call site key-for-key between each pair's `.py`
  and `.ipynb` form and found zero drift, so running either form hits the other's
  cached results. Two examples (`directional_coupler_mode.py`,
  `ring_resonator_hybrid.py`) use raw meep/meow with no FECES solver call to attach
  `save_dir` to; `ring_resonator_hybrid.py` instead fingerprints/caches its own runs
  directly via `feces.cache`'s `lookup()`/`save_run()`/`load_result()` (see that
  script's module docstring). `EMESolver.solve_sweep()` is the one remaining
  solve-family method without `save_dir` support (see Solver recommendation #1).
  **Known structural risk, not fully closed**: any solver/geometry attribute that
  affects the physical result but isn't scalar falls outside `vars(solver)`'s
  generic filter and is silently excluded from the fingerprint (`logger.warning`,
  not an error) — a future attribute added without updating the config builders, or
  a code change to how a position/quantity is *derived* from already-tracked inputs
  (the fingerprint stores ingredients like margins and geometry spans, not derived
  absolute coordinates), could produce a silent false-positive cache hit. No
  regression test guards against this yet.

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
  Each port carries its own `MMIPort(y, waveguide_width, taper_length,
  taper_width)`, so ports can have individually different taper size and arbitrary
  y-placement (the constructor's old scalar `taper_length`/`taper_width`/
  `waveguide_offset` args still work as a convenience default that auto-builds
  symmetric, evenly-spaced `MMIPort` lists). `use_neighbor_aperture` (opt-in,
  default off) generalizes the old `MMIMux`'s neighbor-distance aperture-narrowing
  formula to arbitrary port counts, for devices with several closely-spaced ports
  (e.g. one sitting exactly on the y=0 axis) that the fixed `half_domain_aperture`
  rule can't isolate. `polygons()`/`to_gds_eme()` also take `input_indices`/
  `input_idx` (mirroring `output_indices`/`output_idx`), restricting the drawn
  input stub(s) to a single port — needed for a genuine multi-input device (2×2,
  etc.) so an EME input cross-section is one waveguide, not an unwanted multi-guide
  supermode; `examples/mmi_2x2/` is the first consumer, no 1×4/multi-output-only
  example yet (see Examples recommendation #1).
- **`DirectionalCoupler`** (`coupler.py`) — symmetric 2-waveguide coupler with
  S-bend I/O transitions; built internally from 4 `SBend` sub-geometries (one per
  arm), composed via the shared `_compose.shift_polys`/`shift_prisms` offset
  helpers, rather than hand-rolled polygon math per arm.
- **`Taper`** (`taper.py`) — standalone linear width-transition `DeviceGeometry`
  (width_start → width_end over a length). Used internally by `MMI` for every
  port's taper; also directly usable in any `HybridSolver`/`FDTDSolver`/`EMESolver`
  chain on its own — linear-profile only for now, no standalone example yet (see
  Geometry recommendation #5).
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
  a 3×3 `MMI` this way.
- **`BendGeometry`** (`bend.py`) — legacy standalone 90° routing bend for
  non-collinear FDTD I/O; **not even a `DeviceGeometry` subclass** (a plain
  dataclass with its own hand-rolled minimal GDSII binary writer,
  `.build(device_sx, dpml, output_path)`), so it predates and isn't
  integrated with the `section_prisms`/`width_at`/`center_at` interface at
  all — can't participate in a `HybridSolver` chain like every other
  geometry in this list (see Geometry recommendation #4).

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
  optionally broken down into its `num_modes`/`num_cells` detail. Both
  `plot_fdtd_geometry` and this function's FDTD sub-regions share their
  source/monitor position math with `FDTDSolver.solve_section()` itself
  (`section_fdtd_layout()`, `feces/solvers/fdtd.py`) rather than recomputing it
  independently, so the plot can't silently drift from what actually gets
  simulated. `sub_region_configs` lets one call render every output-port variant
  (e.g. a different `output_idx` per port) in a single pass, re-using a section's
  render across configs whenever its kwargs don't actually change between them —
  this is what the geometry plot in `modal_multiplexer_hybrid.py` uses to render
  three per-port output sections without re-drawing sections 0–3 three times over.
  Neither `ring_resonator_hybrid.py`'s coupler/bend FDTD cells can use any
  `plot_*_geometry` function (they require a `DeviceGeometry`, see Geometry
  recommendation #1) — it hand-draws its own geometry preview plots instead.
- **`plot_eme_fields`** — EME-only: propagates modes through meow cells and saves
  6-component (`|Ex|`,`|Ey|`,`|Ez|`,`|Hx|`,`|Hy|`,`|Hz|`) (x,z) field heatmaps. No
  FDTD time-domain equivalent exists yet (see Plotting recommendation #1).
- **`plot_performance`** — transmission (linear) + insertion loss (dB) vs
  wavelength, side by side, for every port in an `SMatrix`. No phase/group-delay
  equivalent yet, despite `SMatrix` storing full complex coefficients (see
  Plotting recommendation #3).
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
  practice. `save()`/`load()` round-trip via pickle (unversioned, no interop
  format — see Misc recommendation #2). No composition operators (`+`/`@`) to
  cascade several already-solved `SMatrix` results into a larger circuit (see
  Solver recommendation #3).

### Sweep (`feces/sweep/`)

- **`ParallelSweep`** (`parallel.py`) — `ParallelSweep(solver_fn, geometry_fn,
  param_grid, wavelengths, n_workers=4, output_dir=None, save_per_point=True)`.
  `run()` builds the Cartesian product of `param_grid` (`itertools.product`)
  and dispatches each point to a `multiprocessing.Pool(n_workers)` worker,
  which builds a fresh geometry/solver instance and calls plain
  `solver.solve(geometry, wavelengths)` — **always the whole-device `solve()`
  API, never `HybridSolver.solve()`'s `cache=` mechanism** (see Solver
  recommendation #4), so sweeping a hybrid device redundantly resimulates every
  cached section once per worker/grid point. Returns `Dict[param_tuple,
  SMatrix]`; optionally pickles each point's `SMatrix` under
  `output_dir/<param-tag>/data/smatrix.pkl` and a `sweep_results.pkl`
  summary. No cloud/cluster backend beyond the external SLURM `Makefile` in
  `examples/`, and no example demonstrates it yet (see Examples
  recommendation #2).

### Testing (`tests/`)

- **Fast unit tests** — `test_geometry.py`, `test_smatrix.py`, `test_hybrid.py`
  (77 tests). No meep/meow import anywhere (geometry section-cut math is pure
  Python; `HybridSolver` is tested against a `FakeSolver` with scripted
  `solve_section()` returns, not real physics), so these run everywhere in
  well under a second. Run in CI on every push/PR (`.github/workflows/ci.yml`,
  Python 3.10/3.11/3.12 matrix) — no lint or type-checking step (see High-level
  recommendation #4).
- **Slow golden-value integration tests** — `test_golden_examples.py`, one per
  target: `mmi_1x2_eme`/`mmi_1x2_fdtd`/`mmi_1x2_hybrid`/
  `directional_coupler_hybrid` (not `mmi_2x2_*` or `modal_multiplexer_hybrid`
  yet, see High-level recommendations #1-2). Each builds the same device
  geometry as its example at reduced resolution (`res=20` vs the examples'
  `res=40`, 3 wavelengths vs 11), runs a real solve, and checks transmission
  against a recorded reference in `tests/golden/*.json` within a loose (10%
  relative / 3% absolute) tolerance — a regression guard against real
  correctness bugs (sign flips, wrong port routing, normalization mistakes),
  not a physical-accuracy check. Skipped unless `pytest --runslow` is passed;
  `scripts/generate_golden_values.py` (re)generates the reference JSON.
  Genuinely slow even at reduced resolution — real per-target wall-clock on a
  single CPU core: ~7 min (EME), ~12 min (FDTD), ~58 min (hybrid), ~96 min
  (directional coupler hybrid) — so budget ~3 hours for the full suite and
  don't expect this to ever be a pre-commit check. The two hybrid targets
  needed `dpml` doubled (0.5→1.0 µm) versus their examples' value:
  `HybridSolver`'s interior `solve_section()` spans are short enough that PML
  is a large fraction of the domain, and at `res=20` the unmodified `dpml`
  bought only half production's PML pixel count, which diverged (`meep:
  fields NaN/Inf`) after ~20 min instead of converging. Also runs in CI as a
  manual-only `workflow_dispatch` (`.github/workflows/slow-tests.yml`), one
  job per target in parallel via `pymeep` installed from conda-forge — not on
  every push (too slow for that), and not pip-installed (see High-level
  recommendation #6: `pyproject.toml`'s `[fdtd]`/`[all]` extras currently
  declare the wrong `meep` package).
- **Not covered by either suite**: `FDTDSolver`/`EMESolver` used standalone
  outside `HybridSolver` (the golden tests only exercise them via `solve()`/
  `solve_section()` on the four targets above), `mmi_2x2_*`,
  `modal_multiplexer_hybrid` (deliberately excluded — its numerical accuracy
  isn't yet trusted enough to lock in as a golden reference), `ParallelSweep`,
  and all of `feces/plot/`.

## Recommendations

Open items only — sorted into High-level / Solver / Geometry / Plotting / Examples /
Misc, ordered within each section by (impact) × (leverage — how much is genuinely
missing vs. just "meep already does this, expose it"). Anything already implemented
lives in "Current library reference" above, not here.

### High-level

1. **`mmi_2x2` golden test.** `examples/mmi_2x2/` (EME, FDTD, hybrid) has no
   golden-value regression test yet, unlike its `mmi_1x2` counterparts.
   Straightforward to add once the example itself is validated (same pattern as
   the existing 4 targets), but deliberately deferred until then — no value
   locking in a reference before trusting the example's own numbers.
2. **`modal_multiplexer_hybrid` numerical accuracy.** Excluded from the golden
   suite on the same basis as item 1, but blocked on more than just "add a
   test": the device's accuracy isn't yet trusted, so this needs an
   independent correctness check (e.g. against a known-good reference, or at
   minimum a `check_passive()`-style sanity pass plus manual inspection)
   *before* a golden test would mean anything.
3. **Documentation.** No docs site, no changelog beyond this `README.md`/
   `notes.md` pair — the `examples/` scripts are the only tutorials that
   exist. At minimum a `CHANGELOG.md` tracking the still-static
   `__version__ = "0.1.0"`.
4. **Type checking in CI.** Type hints are used throughout (`from __future__ import
   annotations`) but nothing enforces them — no `mypy`/`pyright` config, no `py.typed`
   marker. Cheap to add now that CI exists.
5. **Declare transitive dependencies explicitly.** `gdsfactory`/`tidy3d` are consumed
   indirectly through `meow-sim` (e.g. `meow.compute_modes_tidy3d`) but aren't
   declared in `pyproject.toml` — a `meow-sim` version bump that drops/renames those
   could silently break FECES. Pin or declare them directly, or guard the import
   with a clear error message.
6. **Fix `pyproject.toml`'s `[fdtd]`/`[all]` extras — they install the wrong
   package.** `fdtd = ["meep"]` and `all = ["meow-sim", "meep"]` declare `meep`
   as a pip dependency, but PyPI's `meep` package is an unrelated project (a
   lightweight task runner) — `pip install feces[fdtd]` succeeds with no error
   and leaves `import meep` completely broken. The real MIT electromagnetics
   package only ships via conda-forge (`conda install -c conda-forge pymeep`;
   confirmed working both locally and in `.github/workflows/slow-tests.yml`'s
   CI job). Fix: drop `meep` from these extras (or replace with a clear error
   message on import failure pointing at the conda-forge install command) and
   document the conda-forge requirement in `README.md`'s Install section, not
   just buried in the Testing section.

### Solver

1. **Extend `save_dir` caching to `EMESolver.solve_sweep()`.** Every other
   solve-family method has it (see the Persistent on-disk cache entry above);
   the batched-wavelength whole-device path is the one remaining gap. No
   example calls `solve_sweep()` today, so it hasn't blocked anything yet, but
   it's a real inconsistency now that `solve()` has the parameter.
2. **Multi-excitation `HybridSolver.solve()` API.** `mirror_pairs`/
   `solve_length_sweep()` (see Current library reference above) are fully
   implemented and tested, but neither hybrid example actually calls the
   generic `solve()` API for its full through+cross result: `mmi_2x2_hybrid.py`
   and `directional_coupler_hybrid.py` both solve the input section once per
   arm, derive every output-side matrix by mirror transpose, and hand-roll a
   manual `einsum` cascade — because that composition pairs one arm's
   transposed matrix with the *other* arm's excitation in a single cascade,
   which `solve()`'s one-call contract can't express even with `mirror_pairs`.
   A multi-excitation-config `solve()` extension (declare several excitation
   kwargs up front, get every port-combination coefficient out of one call)
   would let both examples migrate off the manual einsum composition.
3. **Circuit-level S-matrix composition.** `SMatrix` has no way to combine several
   devices' results into a larger circuit (e.g. cascading an MMI's `SMatrix` with a
   directional coupler's). `HybridSolver` already does exactly this kind of
   matrix-chain cascading *within* one device across sections — generalizing that
   same machinery to cascade *across* independently-solved devices (or documenting a
   recipe for handing `SMatrix` results to SAX or Simphony, both established OSS
   tools for exactly this layer) would stop device-level sims from being a dead
   end. (Note: `CompoundGeometry` — see Geometry section — solves a related but
   distinct problem at the geometry layer: it lets `HybridSolver` simulate a device
   built from several `DeviceGeometry` parts as one cascade *before* solving, but
   doesn't let two already-*solved* `SMatrix` results be combined after the fact,
   which is what this item is about.)
4. **`ParallelSweep` + `HybridSolver` integration.** The sweep utility only calls
   whole-device `solver.solve()`; it doesn't know about `HybridSolver.solve()`'s
   `cache=` mechanism, so a parallel sweep over a hybrid device would redundantly
   re-simulate cached sections across worker processes.
5. **Lossy/dispersive material support.** `MaterialStack.si_medium(wl)`/
   `sio2_medium(wl)` (`materials/stack.py`) build an `mp.Medium(index=...)` from
   whichever *tabulated* wavelength in meow's built-in Si/SiO2 data is nearest to
   the requested one (`np.argmin`) — a nearest-neighbor lookup, not interpolation,
   and always real-valued, so there's no true dispersion curve across a wavelength
   sweep and no absorption/loss. Wiring meep's actual Lorentzian-susceptibility
   fitting (or reusing `meep.materials` directly) into `MaterialStack` would enable
   real loss/broadband dispersion and open the door to SiN, III-V, or lossy/doped
   silicon.
6. **Symmetry exploitation in `FDTDSolver`.** `MMI`/`DirectionalCoupler` are
   frequently mirror-symmetric about y=0; `solve_section()` currently always draws
   and simulates the full cross-section. An opt-in `mp.Mirror` symmetry when the
   geometry/excitation is actually symmetric could roughly halve FDTD cost for free —
   complements the reciprocity-based min-simulation-count optimization already there.
7. **Automatic mesh-resolution convergence.** Beyond the manual "run a list of
   resolutions and eyeball it" pattern (see Misc/Examples below) — meep/tidy3d both
   have some form of resolution auto-selection; FECES has none, always a fixed
   user-chosen `res`.
8. **Wavelength-interpolation scheme for mode-solve/EME-driven sections.** A
   real MPB mode solve is expensive per wavelength, while FDTD DFT monitors
   sample a dense wavelength grid nearly for free (one broadband pulse, many
   frequencies in one run). Resonant devices need thousands of wavelength
   points to resolve a linewidth that's a small fraction of the sweep span —
   `examples/ring_resonator/python/ring_resonator_hybrid.py` hits this directly: an
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
   distinct from, Solver recommendation #5 (`MaterialStack`'s nearest-
   tabulated-wavelength index lookup) — that's a materials-dispersion gap;
   this is a mode-solve-cost gap, and the two would compound in a real
   dispersive-material implementation (interpolating dispersive neff at
   already-sparse solve points).
9. **Bigger lift, only if scope grows toward design automation:** adjoint/
   inverse-design support, periodic/grating structures via an RCWA backend
   (S4/Meent/FMMAX-style) for problems where FDTD/EME are the wrong tool,
   near-to-far-field transform. For adjoint specifically, `meep.adjoint`
   (already what FECES's own `FDTDSolver` sits on top of) is the
   lowest-friction path — it's autograd/JAX-based and gradient computation
   wraps around any optimization library, so no new solver dependency would
   be needed, just a FECES-side API for defining a figure of merit over a
   `HybridSolver` chain's `SMatrix` output. tidy3d (native `autograd`
   support) and the newer JAX-native FDTDX project (automatic
   differentiation built into the solver via time-reversibility, published
   in JOSS 2026) confirm autodiff-native is now the field's default
   direction, not a niche add-on — worth knowing if this ever gets
   prioritized. Substantial, well-trodden ground already owned by
   meep/tidy3d/FDTDX/the RCWA tools respectively — worth it only if a
   concrete device need shows up that FDTD/EME genuinely can't address.

### Geometry

1. **Ring/racetrack resonator.** The single most common photonics building block
   that isn't in the device library at all. Natural fit for the existing
   `DeviceGeometry` interface (bus waveguide + coupled ring, `section_prisms`/
   `width_at`/`center_at` following the `DirectionalCoupler` pattern), and directly
   exercises `HybridSolver`'s generic cascading (bus straights as FDTD/EME + coupling
   region). `examples/ring_resonator/python/ring_resonator_hybrid.py` proves out the
   physics (closed-form all-pass feedback assembly, not a `HybridSolver` chain — a
   ring is a loop, not a cascade) entirely in raw meep/meow, bypassing the geometry
   layer and every solver class. What's still missing for a first-class version:
   (a) an actual `RingResonator`/racetrack `DeviceGeometry`; (b) y-normal (x-z
   plane) source/monitor support in `FDTDSolver` — it only ever builds x-normal
   (y-z plane) sources/monitors today (see Current library reference above), and
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
   reintroducing gdsfactory as a hard dependency — every one of gdsfactory+'s 43+
   foundry PDKs is a real device FECES currently has no way to even load, let alone
   simulate.
3. **Y-branch / adiabatic splitter and waveguide crossing.** Common alternatives to
   the MMI splitter and to routing waveguides through each other — neither exists.
4. **Rebuild `BendGeometry` on top of `SBend`** (or a new arc-profile primitive).
   It's a *third*, separately-implemented 90° bend — a legacy standalone GDS-writer,
   predating and not integrated with the shared `section_prisms`/`width_at`/
   `center_at` interface at all (see Current library reference above), and its
   `"euler"` mode isn't actually a distinct curve from `"circular"` (unlike
   `SBend`'s genuine Bezier/circular-arc split). Worth rebuilding it so every bend
   in the library shares one implementation.
5. **`Taper` polynomial/Euler-spiral width profile**, for lower-loss adiabatic
   tapers (linear-profile only today). No example yet demonstrates `Taper`/`SBend`
   used standalone outside their host devices either.

### Plotting

1. **FDTD time-domain field/movie export.** `plot_eme_fields` only covers EME (static
   `|component|(x,z)` heatmaps via `mw.propagate_modes`); there's no equivalent for
   FDTD at all — no meep `Animate2D`-style field snapshot or movie export anywhere in
   the package.
2. **Cascaded full-device field plots for `HybridSolver`.** `plot_eme_fields` only
   ever excites one pure input mode via `mw.propagate_modes(..., excite_mode_l=exc,
   ...)` (a single int, not a superposition), and there's no FDTD field-array export
   at all (see item 1 above) — so today there's no way to visualize the actual field
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
3. **Phase / group-delay plots.** `SMatrix` stores full complex coefficients, but
   `plot_performance`/`plot_splitting_ratio` only ever plot `|coeff|^2` (linear T) and
   dB insertion loss — phase and group delay (from phase vs. wavelength) are never
   visualized despite being available in the data.
4. **A shared convergence plot.** Complements the Misc/Examples resolution-sweep
   item: a `plot_convergence(results_by_resolution, metric)` helper that overlays a
   chosen quantity (e.g. L_c, T) across resolutions on one figure, instead of
   eyeballing separate per-resolution plots side by side (what
   `directional_coupler_mode.py`/`.ipynb` does today).
5. **Interactive/3D visualization.** Everything is a static matplotlib PNG; no
   plotly/bokeh option for geometry or field plots, unlike tidy3d's interactive
   viewer.

### Examples

1. **1×4 (or other multi-output-only) MMI example.** `examples/mmi_2x2/` already
   exercises a genuine multi-input device via `MMI`'s `input_idx`/`input_indices`
   support (see Current library reference above); nothing in the geometry blocks a
   1×4/multi-output-only variant, it just doesn't exist yet.
2. **No example uses `ParallelSweep`.** The library ships a generic parallel
   parameter-sweep utility that zero examples demonstrate — a user has no worked
   reference for the one utility meant to replace ad hoc per-device sweep scripts.
3. **No example demonstrates mesh-convergence checking as a reusable pattern** beyond
   the one-off resolution loop hand-built into `directional_coupler_mode.py`/`.ipynb`
   (a pure meow/MPB even/odd-supermode script that doesn't use any FECES solver
   class at all) — worth generalizing into its own example once the Misc convergence
   utility exists.
4. **Ring resonator example paired with a real `RingResonator` `DeviceGeometry`.**
   `examples/ring_resonator/python/ring_resonator_hybrid.py` proves out the physics
   (see Geometry recommendation #1) and calls `SMatrix.check_passive()` as standard
   practice, but doesn't use the convergence-utility pattern above (no resolution
   sweep) and can't use `plot_fdtd_geometry`/`plot_hybrid_geometry` at all (both
   require a `DeviceGeometry`) — it hand-draws its own geometry preview plots
   instead. The geometry-class-backed version (proper `plot_*_geometry`
   integration, resolution-sweep convergence check) remains open.

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
