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
  self-contained whole-device simulation. `solve_section()` is the piece
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
  recommendation #8 (below) for the full design rationale and
  `mmi_2x2_hybrid.py` for a worked example of both features.

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
7. **Standardize on `logging`, not `print()`.** `EMESolver.solve_section_length_sweep`
   and others `print()` diagnostics unconditionally (e.g. the mid-sweep neff
   diagnostic) while `ParallelSweep` already uses the `logging` module — inconsistent
   for a library whose callers may want to suppress or redirect that output.

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
5. **Non-uniform/graded meshing.** `res` is a single isotropic scalar shared across
   the whole mesh (both FDTD and EME's `Mesh2D`) — no finer mesh near the waveguide
   core and coarser in cladding, unlike tidy3d's automatic adaptive meshing. Would
   cut cell count for a given accuracy target.
6. **Automatic mesh-resolution convergence.** Beyond the manual "run a list of
   resolutions and eyeball it" pattern (see Misc/Examples below) — meep/tidy3d both
   have some form of resolution auto-selection; FECES has none, always a fixed
   user-chosen `res`.
7. **Bigger lift, only if scope grows toward design automation:** adjoint/
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
8. **Mirror-section reuse in `HybridSolver` (reciprocity across sections, not just
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
   populates the shared `cache` dict under its own key. `mmi_2x2_hybrid.py` is the
   worked example: a second `HybridSolver` instance (`hybrid_mirror_bar`, sharing
   `solvers`/`boundaries`/`interface_modes`/`cache` with the plain `hybrid`)
   declares `mirror_pairs={2: 0}` for the bar call only (`output_idx=1` mirrors the
   fixed `input_idx=1` excitation); the cross call (`output_idx=0`) still solves
   its output section directly since it doesn't mirror the excited port. Net: 2
   FDTD + 1 EME, down from 3 FDTD + 1 EME. Tested in `tests/test_hybrid.py`
   (`TestMirrorPairs`/`TestMirrorPairsConstruction`, mock solvers, no meep/meow).
   Length-sweep support landed too: `HybridSolver.solve_length_sweep(geometry,
   wavelengths, lengths, sweep_section, ...)` sweeps one uniform section's length
   analytically (delegating to that solver's `solve_section_length_sweep()`,
   which EMESolver implements — modes solved once per wavelength, lengths applied
   as propagation phases) and returns the raw `[n_L, M, out, in]` cascade grid.
   Non-swept sections keep byte-identical cache keys to `solve()`, so entries
   dedupe freely across both methods and across per-port calls;
   `mmi_2x2_hybrid.py` sweeps its multimode length over 201 points at the same
   2 FDTD + 1 EME mode-solve budget as the single-length run.
   **Still open**: `directional_coupler_hybrid.py`'s hand-rolled bypass is left
   as-is, not migrated onto `mirror_pairs`/`solve_length_sweep` — the length-sweep
   blocker is now gone, but its full through+cross-from-one-excitation-family
   pattern (solving both input arms directly and deriving every output via
   transpose) isn't expressible through `solve()`'s one-call contract even with
   `mirror_pairs` (same reasoning as `mmi_2x2`'s cross case) — would need a
   multi-excitation-config `solve()` extension or staying manual. Similarly,
   getting `mmi_2x2_hybrid.py`'s full 2×2 port-combination matrix "for free"
   (zero dedicated output solves, matching the DC pattern exactly) isn't done —
   would require solving both input branches directly rather than one fixed
   excitation, for no additional savings over the current 2 FDTD + 1 EME on this
   example's actual bar/cross metric.

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

1. **Stop duplicating solver position/size math in the plot functions.**
   `_draw_section_fdtd_markers` (plot/geometry.py) independently recomputes where
   `FDTDSolver.solve_section()` places its source/monitor apertures, and an
   earlier bug-hunt found two real mismatches there (marker offset by
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
   `directional_coupler_mode.py`/`.ipynb` does today).
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
3. **2×2 MMI example — done; 1×4 still missing.** `examples/mmi_2x2/` (EME, FDTD,
   hybrid variants) exercises a genuine multi-input device — the first consumer
   of `MMI`'s new `input_idx`/`input_indices` support (see `MMI` in the library
   reference above). No 1×4 (or other multi-output-only) example yet, though
   nothing in the geometry blocks one.
4. **No example uses `ParallelSweep`.** The library ships a generic parallel
   parameter-sweep utility that zero examples demonstrate — a user has no worked
   reference for the one utility meant to replace ad hoc per-device sweep scripts.
5. **No example demonstrates mesh-convergence checking as a reusable pattern** beyond
   the one-off resolution loop hand-built into `directional_coupler_mode.py`/`.ipynb`
   (a pure meow/MPB even/odd-supermode script that doesn't use any FECES solver
   class at all) — worth generalizing into its own example once the Misc convergence
   utility exists.
6. **Ring resonator example** — once the Geometry-section ring/racetrack class
   exists, pair it with an example (and it'd make a good first user of the
   convergence-utility recommendation above, plus `SMatrix.check_passive()` as
   standard practice, exercising several new/existing features at once).

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
