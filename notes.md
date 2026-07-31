# BEEF capability gap analysis vs. meep / meow / tidy3d / other open-source EM tools

Research memo — comparison of BEEF against meep, meow, tidy3d, and the wider open-source
photonics/EM landscape, plus a prioritized list of open items. Everything under
"Recommendations" is unimplemented unless the item text says otherwise; "Current library
reference" is the descriptive record of what's actually in the tree today. Last refreshed
2026-07-21 against the current `beef/` source and a fresh pass over the OSS landscape (see
sources inline in the table). Project renamed FECES → BEEF (Bridged Engine for EME and
FDTD); the Python package/import path was renamed to match (`feces` → `beef`), so `import
beef` throughout this file and the codebase.

**What BEEF actually is**: a hybrid-methodology orchestration layer on top of meep
(FDTD), MPB/meow (mode solving), and optionally tidy3d's mode solver — plus a small
parametrized photonic-device geometry library (`MMI` with per-port `MMIPort` taper
specs, `DirectionalCoupler`, standalone `Taper`/`SBend` building blocks, `CompoundGeometry`/
`WaveguideRun` for composing several devices into one cascade, and legacy `BendGeometry`).
It is **not** a competing raw solver;
meep/meow/tidy3d do all the actual field/mode math. BEEF's real contribution is
`HybridSolver`: a generic N-section FDTD/EME cascade (any solver type, any order, any
count ≥ 2) with reciprocity-based excitation (minimum-simulation-count) and a shared
`section_prisms/section_gds/width_at/center_at` interface that lets any geometry slot
into that cascade.

## Comparison table

| Capability | BEEF today | meep | meow | tidy3d | Other OSS |
|---|---|---|---|---|---|
| Time-domain FDTD | 2D or 3D via `FDTDSolver(three_d=...)`, single fixed uniform `res` | Native, mature, adaptive subpixel smoothing | — | Cloud GPU FDTD, auto mesh refinement | openEMS (general 3D FDTD, antenna-focused); **FDTDX** (JAX-based, GPU-accelerated, native autodiff via time-reversibility, published in JOSS 2026 — 10-415x faster than traditional FDTD per its own benchmarks); Ceviche/EMopt (earlier autograd-based FDTD); GSvit (GPU-accelerated) |
| Eigenmode expansion (EME) | `EMESolver` wraps meow's cross-section mode-match + `compute_s_matrix` | — | Native (this *is* meow's job) | Has an EME plugin (steady-state frequency-domain alternative to its FDTD solver) | EMEpy (pure-Python, neural-net-accelerated mode profiles); CAMFR (older EME framework) |
| Hybrid FDTD+EME chaining | **Yes — `HybridSolver`, BEEF's actual differentiator** | No | No | No (closest: manually stitching, no generic N-section cascade) | Not found in a fresh landscape pass (`awesome_photonics`, direct searches) — closest is `simphox` (Stanford Solgaard Lab, bundles FDTD/beam-propagation/circuit solvers), but it isn't documented as a generic reciprocity-based N-section cascade the way `HybridSolver` is |
| Dispersive materials | One shared analytic Sellmeier model (Salzberg Si / Malitson SiO2) evaluated per wavelength by both engines; still real-valued, and the FDTD still freezes one index per broadband run — see recommendation #5 | Full Lorentzian/Drude/Debye susceptibility fits | Has real material models (used via meow calls) | Automatic dispersion-fitting plugin | — |
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
| GDS import/export, PDK interop | **None** — gdsfactory dependency was actively removed (commit `ffcbf5a`); geometries are BEEF-internal polygon dicts only, can't simulate a real fab layout or export one for tape-out | n/a | Uses gdsfactory `Component`s | Has gdsfactory integration | gdsfactory itself |
| Ring/racetrack resonator | *(2026-07-29)* `RingResonatorCouplerSection`/`RingResonatorBendSection` (`beef/geometry/`) plus `FDTDSolver`'s new y-normal `PortSpec` support and `beef.network.close_loop()` exist as standalone, tested primitives — see Geometry recommendation #1 — but nothing wires them together yet; `examples/ring_resonator/` is unchanged, still raw meep/meow | n/a | n/a | n/a | gdsfactory has one |
| Grating coupler / edge coupler geometry | Missing | n/a | n/a | n/a | gdsfactory |
| Circuit-level composition (cascade many devices' S-matrices into one chip netlist) | *(2026-07-29)* `beef.network` — general Redheffer-star-product `connect()`/`close_loop()`/`Circuit`, verified against the ring resonator's proven closed-form formula and against `mmi_2x2_hybrid.py`'s hand-rolled cascade; not yet wired into `HybridSolver` or any example (see Network below) | n/a | n/a | n/a | **SAX** (JAX/netlist-driven, differentiable) and **Simphony** (component-library/plugin-based, itself built on SAX) are the two standard choices; **opics**, **lekkersim**, **photontorch** (adds time-domain circuit analysis) round out the space |
| Result serialization | Pickle only, no schema/versioning | n/a | n/a | Its own `.hdf5`/JSON | Touchstone `.sNp` is the RF/photonics interop standard, not supported |
| Passivity/unitarity sanity check | `SMatrix.check_passive()` — implemented, called after every `solve()`/`hybrid.solve()` in every example and golden test as standard practice | n/a | n/a | n/a | — |
| Parallel parameter sweeps | `ParallelSweep` — local `multiprocessing.Pool` grid sweep, but only calls whole-device `solver.solve()`, **not** `HybridSolver`-aware, no cloud/cluster backend beyond an external SLURM Makefile | n/a (user-driven) | n/a | Cloud batch is tidy3d's core value prop | — |
| Field visualization | Static EME field heatmaps only (`plot_eme_fields`); **no FDTD time-domain field/movie export at all** | `mp.Animate2D`, live field plotting | n/a | Interactive Plotly 3D viewer | — |
| Interactive/3D visualization | None (matplotlib static PNGs only) | Some (via `mp.plot2D`) | n/a | Yes, core feature | — |
| Testing | 77 fast unit tests (`tests/`, meep/meow-free) run in CI on every push; 4 slow golden-value integration tests (real reduced-res meep/meow sims) run in CI as a manual-only `workflow_dispatch` matrix via conda-forge `pymeep` (`pytest --runslow` locally); no lint/type-checking yet | n/a | n/a | n/a | — |

## Current library reference

Descriptive reference for what's actually in the library today — not recommendations,
just what each solver/geometry/plotting entry point does and how they fit together.

### Solvers (`beef/solvers/`)

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
- **`HybridSolver`** — BEEF's actual differentiator: a generic N-section (≥2)
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
- **Persistent on-disk cache** (`beef/cache/`) — every solve-family method
  (`FDTDSolver.solve()`/`solve_section()`, `EMESolver.solve()`/`solve_section()`/
  `solve_section_length_sweep()`, `HybridSolver.solve()`/`solve_length_sweep()`)
  takes an optional `save_dir: Optional[str] = None` (falsy = disabled, no behavior
  change if unset). `build_run_config()`/`build_whole_device_config()`
  (`beef/cache/config.py`) fingerprint a call into a JSON-safe dict — solver
  params (`vars(solver)`, `material_stack` serialized explicitly), the call's own
  kwargs, and the exact geometry drawn (`geometry.section_gds()`/`.polygons()`/
  `.to_gds_eme()`, no meep/meow objects involved) — hashed (`sha256(...)[:16]`)
  into a `save_dir` subfolder holding `config.json` (human-readable), `result.pkl`,
  and `geometry.png`; a hash hit still gets a full `config.json` equality check
  before being trusted, guarding against truncated-hash collisions or
  config-schema drift across code versions. `save_dir` also accumulates a
  top-level `index.jsonl` (one compact summary line — hash, solver_type,
  method, geometry_class, key call kwargs — appended per `save_run()` call,
  long numeric arrays like `wavelengths`/`lengths` collapsed to
  `{count, min, max}`) so the whole cache can be `grep`/`jq`-scanned without
  opening every subfolder's `config.json`; not deduplicated (re-saving an
  already-cached config just appends another line for the same hash) and
  failure-isolated (a manifest-write error is logged and doesn't affect the
  actual cache save), same pattern as `geometry.png`'s own best-effort save
  (the first half of what was Solver recommendation #3 — see further down
  this entry for the folder-naming half). This is a separate mechanism from
  `HybridSolver`'s in-memory `cache={}` dict (dedupes within/across `solve()` calls
  in one process; `save_dir` persists to disk across process runs) — the two
  compose freely. All 12 example `.py`/`.ipynb` pairs across every device folder
  point at one shared `example_sim_cache/` directory (a sibling of `BEEF/`,
  outside the repo since 2026-07-28 — it had been `examples/example_sim_cache/`
  and 56 MB of it was tracked in git); a 2026-07 audit
  diffed every cache-fingerprinting call site key-for-key between each pair's `.py`
  and `.ipynb` form and found zero drift, so running either form hits the other's
  cached results. Two examples (`directional_coupler_mode.py`,
  `ring_resonator_hybrid.py`) use raw meep/meow with no BEEF solver call to attach
  `save_dir` to; `ring_resonator_hybrid.py` instead fingerprints/caches its own runs
  directly via `beef.cache`'s `lookup()`/`save_run()`/`load_result()` (see that
  script's module docstring). `EMESolver.solve_sweep()` is the one remaining
  solve-family method without `save_dir` support (see Solver recommendation #1).
  Every fingerprint also carries a hand-bumped `cache_schema_version` int
  (`_CACHE_SCHEMA_VERSION` in `beef/cache/config.py`) — bumping it on a
  solver-internals change (even one no single config field would ever
  capture by hand) invalidates every existing cache entry outright, since
  the hash changes; old entries are never read/deleted, just no longer
  matched, and get transparently resimulated/re-saved under the new hash.
  `_solver_config()` is also now scoped per `(solver_type, method_name)` to a
  hand-traced allowlist of the attributes that method actually reads
  (`_METHOD_ATTRS`, same file) rather than the raw `vars(solver)` dump — e.g.
  `FDTDSolver.solve_section()`'s cache key no longer includes
  `num_input_modes`/`num_output_modes` (constructor params only `solve()`
  reads), so two solvers differing only in those no longer fork the cache
  for what's actually the same physical result; an unrecognized
  `(solver_type, method_name)` pair falls back to the old "include
  everything" behavior (the safe direction — over-inclusion costs a spurious
  miss, under-inclusion risks a silent false-positive hit). Finally,
  `FDTDSolver.solve_section()` fingerprints the fully-resolved
  `SectionFDTDLayout` (`resolved_layout` in its `call_kwargs`, computed
  before the cache-config block precisely so it can be), not just the
  `monitor_center`/`monitor_size`/margin ingredients that derive it — a
  future change to that derivation logic now shows up as a config diff
  instead of silently returning a stale result for the same nominal inputs
  (this and the two paragraphs above fully close what was Solver
  recommendation #4 — see `tests/test_cache_config.py` for the
  allowlist/schema-version/resolved-layout regression coverage).
  `beef/cache/store.py`'s folder-naming scheme is now `{slug}_{hash(config)}`
  (`_slug()`/`_folder()`) — a human-readable `{solver_type}_{method}_
  {geometry_class}` prefix ahead of the hash, so `ls save_dir` is directly
  browsable — with `lookup()` falling back to the pre-slug bare-hash path so
  every cache entry from before this change stays reachable at its old
  location; only new writes get the new name (closes Solver recommendation
  #3 in full).
  A second, orthogonal lookup path — `lookup_superset()`/
  `update_family_pointer()` — lets a request for fewer modes than a prior
  call already computed reuse that larger result by slicing its mode axes,
  instead of resimulating: both `FDTDSolver.solve_section()`'s and
  `EMESolver.solve_section()`/`solve_section_length_sweep()`'s mode-solve
  paths decode/compute each mode independently of how many others are
  requested (confirmed by reading `sim.get_eigenmode_coefficients()`'s
  per-band loop and meow's `compute_modes_meep()`'s per-`band_num` loop
  respectively), so a run performed at a larger mode count already contains
  everything a smaller request needs. A "family" config — the same
  fingerprint with the mode-count fields dropped (plus, for FDTDSolver only,
  an explicit `reverse` field standing in for them, and gated on
  `monitor_size` being pinned — see below) — points at the largest
  mode-count result computed so far for that family
  (`save_dir/_mode_family_pointers/`, one small JSON file per family,
  atomically written since concurrent SLURM array tasks may race to update
  it); a request whose mode counts are all `<=` the pointer's recorded dims
  is served by slicing, a request that exceeds it falls through to a real
  simulation whose result then becomes the new family maximum (if it
  dominates the old one — an incomparable growth pattern, e.g. one axis
  growing while another shrinks, just isn't recorded as a new maximum,
  falling back to a plain resimulation next time too). Two preconditions
  gate this for `FDTDSolver.solve_section()` specifically, both because mode
  count can otherwise change the actual simulated geometry, not just how
  many bands get decoded from it: (a) `monitor_size` must be explicitly
  passed (its *default* width depends on whether `n_measure > 1`, which
  would make the physical monitor itself mode-count-dependent), and (b) the
  family key includes `layout.reverse` (which side excites vs. measures) so
  requests that would flip it — e.g. crossing `n_modes_out` from 1 to
  2 while `n_modes_in=1` — correctly form separate families rather than
  incorrectly reusing across a boundary where the simulation is physically
  different. `EMESolver` has no such asymmetry (a single joint mode solve
  serves both the in/out split with no monitor concept), so its family key
  is just the mode-count-free fingerprint, no extra preconditions. Verified
  against real (non-cache) solves for both solvers — a direct small-N
  computation matches a large-N computation sliced down — subject to the
  same near-degenerate-mode-reordering caveat noted below (fully closes what
  was Solver recommendation #2; see `tests/test_cache_superset.py`).
  **Known structural risk, still open**: any solver/geometry attribute that
  affects the physical result but isn't JSON-safe still falls outside
  `_solver_config()`'s filter and is silently excluded from the fingerprint
  (`logger.warning`, not an error) — a future non-scalar attribute added
  without updating the config builders could still produce a silent
  false-positive cache hit. No regression test guards against this.

### Geometries (`beef/geometry/`)

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

### Plotting (`beef/plot/`)

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
  (`section_fdtd_layout()`, `beef/solvers/fdtd.py`) rather than recomputing it
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

### Results (`beef/results/`)

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
  cascade several already-solved `SMatrix` results into a larger circuit —
  **partly closed 2026-07-29, see Network below.**

### Network (`beef/network/`)

*(Added 2026-07-29, closing the composition half of Solver recommendation #3
and the `ring_resonator_validation.md` "general form" recommendation — the
same underlying gap named independently in two places, see that file's
Recommendation section.)*

- **`SMatrixBlock`** (`block.py`) — dense `[M, N, N]` S-parameter block over
  named, multi-mode ports (`data[m] @ a[m] = b[m]`), distinct from `SMatrix`'s
  flat sparse `(port, in_mode, out_mode)` dict: every port here is
  bidirectional (carries its own reflection/self-term), which is what
  `connect()`/`close_loop()` need to do real linear algebra on.
  `submatrix(port_out, port_in)` is the multi-mode generalization of what
  `solve_section()` already returns for one directional pair.
- **`connect(block_a, port_a, block_b, port_b)`** — the general Redheffer
  star product joining one named port of two (possibly different) blocks;
  every other port of both survives. When both blocks come from
  `from_transfer_matrix()` (zero self-terms, matching what
  `solve_section()` actually models), this reduces *exactly* to ordinary
  transfer-matrix multiplication — verified against `mmi_2x2_hybrid.py`'s
  own hand-rolled `einsum` cascade to `1e-12` (`tests/test_network.py`).
- **`close_loop(block, port_out, port_in, loop=None)`** — eliminates a
  block's own port against itself (node elimination); only the
  `port_out == port_in` case is implemented — the `!=` case (closing two
  *different* ports of an already-merged block, needed for some
  `Circuit`-discovered cycles) is a materially different formula and raises
  `NotImplementedError`, deferred. Verified to reproduce the ring
  resonator's proven closed-form all-pass formula
  (`ring_resonator_hybrid.py:762-774`, `E_thru = t + k1·k2·P/(1-t_ring·P)`)
  term for term when applied to a 2-port `[bus, ring]` block with
  `loop=P` — this is the concrete evidence the general primitive is
  physically correct, not just internally consistent
  (`tests/test_network_ring_crosscheck.py`).
- Both raise on a near/exactly-singular star-product or loop-closure
  system, via the reciprocal smallest singular value rather than
  `np.linalg.cond()` — the latter is identically `1.0` for any 1×1 matrix
  regardless of magnitude (confirmed: `cond([[1e-9]]) == 1.0`), which would
  silently defeat the check for exactly the real use case, a single-mode
  ring's loop port. A near-critical-coupling ring genuinely produces this
  condition; the warning/error is not a numerical artifact to raise the
  threshold past.
- **`Circuit`** (`circuit.py`) — stateful multi-block builder
  (`add_block`/`connect`/`close_loop`/`reduce`), namespacing each
  registered block's ports as `f"{name}.{port}"` so two instances of the
  same device never collide. A queued `connect()` between two names already
  merged into one component is automatically dispatched to `close_loop()`
  instead — the one caller `close_loop`'s `port_out != port_in` branch
  would need, deferred along with it. Once wired up, extracting a device's
  full multi-excitation response from one `Circuit.reduce()` + `to_smatrix()`
  call is a byproduct that subsumes Solver recommendation #2's
  multi-excitation `solve()` ask — noted, not built as a dedicated
  `HybridSolver` API.
- **Deliberately not wired up yet**: nothing in `beef/solvers/` or any
  example imports `beef.network` — it depends only on `SMatrix`/`numpy`.
  `HybridSolver`'s own cascade is untouched; the ring resonator examples
  still hand-roll their closed-form assembly directly. *(2026-07-29 —
  `ring_resonator_validation.md`'s other two prerequisites, y-normal
  `FDTDSolver` ports and a `RingResonator` `DeviceGeometry`, are now also
  done — see Solvers/Geometry above and Geometry recommendation #1 below —
  so nothing structural blocks the wiring anymore; it's just not done.)*
  32 new fast tests, no meep/meow import, full suite unaffected:
  282 → 314 passed.

### Sweep (`beef/sweep/`)

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

- **Fast unit tests** — `test_geometry.py`, `test_smatrix.py`, `test_hybrid.py`,
  `test_cache.py`, `test_cache_config.py`, `test_cache_superset.py` (116
  tests). No meep/meow import anywhere (geometry section-cut math is pure
  Python; `HybridSolver` is tested against a `FakeSolver` with scripted
  `solve_section()` returns, not real physics; the cache tests use bare fake
  objects plus real `FDTDSolver`/`EMESolver` instances constructed but never
  actually solved — `solve()`/`solve_section()`/`solve_section_length_sweep()`
  all do an unconditional `import meep`/`import meow` even down a
  cache-hit/`cache_only=True` path, so exercising that wiring means calling
  `section_fdtd_layout()`/`build_run_config()`/`lookup_superset()` etc.
  directly instead — see `test_cache_config.py`/`test_cache_superset.py`'s
  own module docstrings), so these run everywhere in
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
  and all of `beef/plot/`.

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
   could silently break BEEF. Pin or declare them directly, or guard the import
   with a clear error message.
6. **Fix `pyproject.toml`'s `[fdtd]`/`[all]` extras — they install the wrong
   package.** `fdtd = ["meep"]` and `all = ["meow-sim", "meep"]` declare `meep`
   as a pip dependency, but PyPI's `meep` package is an unrelated project (a
   lightweight task runner) — `pip install beef[fdtd]` succeeds with no error
   and leaves `import meep` completely broken. The real MIT electromagnetics
   package only ships via conda-forge (`conda install -c conda-forge pymeep`;
   confirmed working both locally and in `.github/workflows/slow-tests.yml`'s
   CI job). Fix: drop `meep` from these extras (or replace with a clear error
   message on import failure pointing at the conda-forge install command) and
   document the conda-forge requirement in `README.md`'s Install section, not
   just buried in the Testing section.
7. **Golden tests shouldn't re-implement the examples they guard.** Each
   `run_*()` in `tests/test_golden_examples.py` rebuilds its example's
   geometry, solver config and `HybridSolver` boundaries by hand at reduced
   resolution/wavelength/mode counts — a second copy of the thing under test.
   When the directional-coupler hybrid moved its interfaces 0.25 µm into the
   coupling region (2026-07-24), the replica had to be edited to match and its
   golden JSON regenerated; had it been left alone it would have kept passing
   happily against the *old* convention, i.e. guarding a configuration no
   example uses any more. Better: have each example expose its build as a
   factory the test imports (parametrized by resolution/wavelength count so the
   test can still run cheap), or failing that, add a test asserting the shared
   invariants — boundary positions relative to the geometry, solver flags —
   agree between example module and replica.

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
   would let both examples migrate off the manual einsum composition. Side
   effect of the bypass worth fixing with it: because `solve()` derives each
   section's length from the boundaries, it can't get the section lengths
   wrong — the hand-rolled path can, so each example now asserts its own
   length budget (FDTD stubs + EME body == the nominal body length) in
   script-level code that a `solve()`-based version wouldn't need at all.

   *(2026-07-29 — the primitive that would replace this now exists, not yet
   used here.)* `beef.network.Circuit` can register both excitations'
   transfer matrices and extract every port-combination coefficient from one
   `reduce()` + `to_smatrix()` call — see Network below. Neither example has
   been migrated to it; this item is about a dedicated `HybridSolver.solve()`
   API, which is still not built.
3. **Circuit-level S-matrix composition.** *(Composition primitive done
   2026-07-29 — see Network below; not yet wired into `HybridSolver` or any
   example.)* `SMatrix` itself still has no `+`/`@` operators, and nothing in
   `beef/solvers/` calls the new `beef.network` module — this item's
   *library-integration* half (making composition a first-class part of
   `HybridSolver`'s own API, the way `SMatrix` is today) is still open.
   `CompoundGeometry` — see Geometry section — remains a related but distinct
   problem at the geometry layer: it lets `HybridSolver` simulate a device
   built from several `DeviceGeometry` parts as one cascade *before* solving;
   `beef.network` combines already-*solved* results *after* the fact, which
   is what this item was originally about.
4. **`ParallelSweep` + `HybridSolver` integration.** The sweep utility only calls
   whole-device `solver.solve()`; it doesn't know about `HybridSolver.solve()`'s
   `cache=` mechanism, so a parallel sweep over a hybrid device would redundantly
   re-simulate cached sections across worker processes.
5. **Lossy/dispersive material support.** *(Partly done 2026-07-28.)* The
   nearest-tabulated `np.argmin` lookup is gone: `materials/stack.py` now holds one
   analytic Sellmeier model (Salzberg Si, Malitson SiO2, lifted from `meep.materials`
   and verified against it to 4e-16) that **both** engines evaluate, so the EME and
   the FDTD agree on n(λ) to machine precision — closing the second half of
   recommendation #9 below.

   Still open, and now for a *measured* reason: the FDTD is still handed one
   `mp.Medium(index=...)` per broadband run, so it has no dispersion across the
   sweep, and nothing is lossy. Turning on the dispersive medium was tried and
   reverted — meep subpixel-averages only the instantaneous ε and leaves
   susceptibility terms unsmoothed, and the Sellmeier form (like meep's own library
   materials) puts ε∞ = 1 with the whole index in susceptibilities. meep then builds
   a structure with `epsilon diagonal = (1,1,1)`, i.e. no contrast to smooth, and on
   a 450x220 nm strip at res 40 n_eff moved 2.313398 → 2.363038, **+0.0496** —
   exactly the cost of losing smoothing that this file's own #9 predicted, against
   the ~0.002 the dispersion was worth. Doing it properly means re-parameterising so
   ε∞ carries the bulk of the permittivity and the susceptibilities only its
   wavelength dependence, then re-running the resolution sweeps. `meep_medium()` is
   kept in `materials/stack.py` as the correct model to start that from.

   *(2026-07-31 — the re-parameterisation landed and is opt-in via
   `dispersive=`; see `material_model_validation.md`. One trap found the same day
   and worth carrying forward: **turning a flag on is not the same as the flag
   reaching the physics.** `FDTDSolver.section_mode_neff()` accepted
   `dispersive=True` and returned n_eff byte-identical to `dispersive=False` for
   ~1.5 h of node time, because `_build_section_mode_sim()` built its media
   directly instead of through `_media()`. Nothing in a 468-test suite noticed,
   because nothing asserted on which `mp.Medium` was constructed — only on
   numbers downstream of it. Any future material or solver flag wants a test that
   pins the **object built**, plus a runtime liveness check in the script that
   sweeps it, which is what `verify_item_8_dispersive_convergence.py` now
   carries.)*
6. **Symmetry exploitation in `FDTDSolver`.** `MMI`/`DirectionalCoupler` are
   frequently mirror-symmetric about y=0; `solve_section()` currently always draws
   and simulates the full cross-section. An opt-in `mp.Mirror` symmetry when the
   geometry/excitation is actually symmetric could roughly halve FDTD cost for free —
   complements the reciprocity-based min-simulation-count optimization already there.
7. **Automatic mesh-resolution convergence.** Beyond the manual "run a list of
   resolutions and eyeball it" pattern (see Misc/Examples below) — meep/tidy3d both
   have some form of resolution auto-selection; BEEF has none, always a fixed
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
   nothing in `beef/solvers/eme.py` supports it. A library-level version
   could look like sparse-solve/dense-eval wavelength arguments on
   `EMESolver.solve_section()` (solve modes at a caller-chosen sparse subset,
   fit, evaluate the returned transfer matrix at every wavelength in a
   separate dense array), or a standalone reusable `NeffTable`/mode-fit
   utility both `EMESolver` and examples could share. Related to, but
   distinct from, Solver recommendation #5 (`MaterialStack`'s material model)
   — that's a materials-dispersion gap;
   this is a mode-solve-cost gap, and the two would compound in a real
   dispersive-material implementation (interpolating dispersive neff at
   already-sparse solve points).
9. **Make the two section types agree on the material at a given wavelength.**
   *(Half done 2026-07-28 — the disagreement is closed, the frozen band is not.)*

   Within one `HybridSolver` chain the FDTD and EME sections used to look silicon
   and oxide up *differently*: `FDTDSolver` builds both media **once**, at the mean
   of the requested wavelength array (`_build_domain`, `solvers/fdtd.py`), so every
   wavelength in a sweep is simulated in the *same* frozen structure, while
   `EMESolver` rebuilds the cross-section per wavelength. On top of that each side
   ran its own nearest-tabulated `argmin`, so at the same wavelength the two halves
   could sit on *different* silicon: meow's table has only 1.532 and 1.600 µm in
   1.49-1.61 µm, so an 11-point 1.50-1.60 µm sweep put the EME on n=3.4784 below
   ~1.566 µm and 3.4710 above it while the FDTD stayed on 3.4784 throughout —
   ~0.0074 in index, a silent inconsistency inside a single device solve.

   **Closed:** both engines now evaluate one shared analytic Sellmeier model
   (recommendation #5), agreeing to 4e-16 at any wavelength, and the step in the
   EME's n(λ) — the Δn_eff kink of `directional_coupler_validation.md` item 8 — is
   gone. The `argmin` in `meow/fde/meep.py::meep_material()` was the direct cause
   and was *bypassing* meow's own correctly-interpolating `SampledMaterial.__call__`;
   both copies were removed.

   **Still open:** the FDTD half still has no dispersion across a sweep, by
   construction — one Medium per broadband run. See #5 for why the obvious fix
   (a dispersive Medium) measured 0.0496 worse and was reverted. A `HybridSolver`
   assertion that the two section types resolve the same index at the same
   wavelength is still worth adding; it is now cheap, since both call the same
   function. The checked re-parameterisation that would unblock both — fit
   band, λ₀ degeneracy, strip-test gate — is in
   `validation_markdown/material_model_validation.md`.

   **Subpixel-smoothing coupling — this prediction was tested, and it held.**
   Written before the attempt, and it is exactly what happened on 2026-07-28: see
   #5. Keeping it as the standing reason the media are still
   `mp.Medium(index=...)`. meep only subpixel-averages the *instantaneous* epsilon. So today smoothing applies in
   full on both paths — meow sets `eps_averaging=True` explicitly for its mode
   solve (`meow/fde/meep.py:97`) and `FDTDSolver` inherits meep's `True`
   default (it never passes the flag, nor `subpixel_tol`/`subpixel_maxeval`) —
   and it is not a small effect on the mode basis: disabling it shifts the
   fundamental-mode `neff` of the 0.45x0.22 µm strip by 0.050 (meow-style 2-D
   cross-section solve, res 40) and 0.033 (3-D prism, eigenmode on an x-normal
   plane, res 30). Swapping in Lorentzian susceptibilities changes that: meep
   averages only the ε∞ part and assigns the susceptibility terms per pixel
   unsmoothed, so the *same* ε(ω) expressed dispersively rather than as a fixed
   index yields interface pixels differing by up to ~0.23 in ε. A "just the
   materials" change would therefore quietly move the mode basis and the
   mesh-convergence point along with it — budget a re-run of the resolution
   sweeps (`batch/archive/run_directional_coupler_*_resolution_sweep.batch` --
   moved there 2026-07-31, retired pending new sweeping functionality; still
   the literal scripts this budget refers to) as part of
   #5, not as a follow-up. While in there: `FDTDSolver` neither exposes nor
   logs `eps_averaging`/`subpixel_tol`/`subpixel_maxeval`, so a knob that
   demonstrably moves `neff` at the second decimal is set entirely by meep's
   default. *(The `FDTDSolver` half is now done — `eps_averaging`,
   `subpixel_tol` and `subpixel_maxeval` are real params on the `_METHOD_ATTRS`
   allowlist. The material indices themselves reached the fingerprint on
   2026-07-28 via `resolved_material_config()`; before that a `MaterialStack`
   was recorded as `{wg_height, sz}` only, so any change to the index model kept
   matching every cached entry.)*
10. **Give hybrid interface placement a library-level story.** Every hybrid
    example folds a stub of the adjacent uniform region into its FDTD sections
    so the interface mode decomposition lands past a cross-section
    discontinuity, and every one picks the amount by hand, with no shared
    rationale: 0.25 µm in `mmi_1x2_hybrid.py:139`, `mmi_2x2_hybrid.py:168` and
    (since 2026-07-24) `directional_coupler_hybrid.py`, 0.75 µm in
    `modal_multiplexer_hybrid.py:162`, 1.0 µm as `IFACE_BUFFER` in
    `ring_resonator_hybrid_cascaded.py:212` — and the directional coupler used
    **zero** (interfaces exactly on the S-bend/coupling-region seam) until that
    date. The library offers no help here: `HybridSolver` takes raw boundary
    positions, and geometries validate only their own hard constraints
    (`DirectionalCoupler._check_cut` rejects a cut strictly inside an S-bend,
    but a cut exactly *on* the seam — the worst case for the decomposition,
    since `section_fdtd_layout()` puts the interface monitor exactly at the
    section boundary, giving the bend's near field zero relaxation length —
    is legal and silent). Three concrete pieces: (a) a geometry-level
    `discontinuities()` (or recommended-stub hint) so `HybridSolver` can warn
    when a boundary coincides with a cross-section change; (b) **done
    2026-07-30** — the *comparison* half of the empirical procedure
    (`directional_coupler_hybrid_stub_sweep.py` sweeps the stub 0:0.25:2 µm
    on one fixed device and plots the 50:50 split length against it) is now
    reusable: `beef.cut_position_invariance.cut_position_invariance()` takes
    one already-solved `SMatrix` per stub value and fits transmission drift
    vs. stub for every `(port, in_mode, out_mode)` key, reporting the fitted
    slope. The device-specific half — building each stub's boundaries/
    geometry/kwargs, and picking a scalar observable like the coupler's 50:50
    length — deliberately stays with the caller, since that bookkeeping
    really is per-device. No example has been migrated onto it yet (needs a
    real meep/meow run to confirm it reproduces the coupler's known 167x
    drift-collapse, not just a code change); (c) once each geometry family has
    a defensible number, keep it in one place instead of four example-local
    constants.
11. **Bigger lift, only if scope grows toward design automation:** adjoint/
    inverse-design support, periodic/grating structures via an RCWA backend
    (S4/Meent/FMMAX-style) for problems where FDTD/EME are the wrong tool,
    near-to-far-field transform. For adjoint specifically, `meep.adjoint`
    (already what BEEF's own `FDTDSolver` sits on top of) is the
    lowest-friction path — it's autograd/JAX-based and gradient computation
    wraps around any optimization library, so no new solver dependency would
    be needed, just a BEEF-side API for defining a figure of merit over a
    `HybridSolver` chain's `SMatrix` output. tidy3d (native `autograd`
    support) and the newer JAX-native FDTDX project (automatic
    differentiation built into the solver via time-reversibility, published
    in JOSS 2026) confirm autodiff-native is now the field's default
    direction, not a niche add-on — worth knowing if this ever gets
    prioritized. Substantial, well-trodden ground already owned by
    meep/tidy3d/FDTDX/the RCWA tools respectively — worth it only if a
    concrete device need shows up that FDTD/EME genuinely can't address.
12. **New, 2026-07-30 — `FDTDSolver`'s `get_eigenmode()` band loop can return
    a duplicate eigenvector. DONE 2026-07-30 for `section_mode_neff()`/
    `section_mode_fields()`; `solve_section()` deliberately deferred.**
    Found auditing the modal multiplexer's 16-mode interface basis
    (`validation_markdown/modal_multiplexer_validation.md` gap 1) with the
    new field-overlap mode-identity check (below): two consecutively-
    requested bands came back as the exact same eigenvector (self-overlap
    0.9999999999, not a near-degenerate-but-distinct pair — every other band
    pair on the same cross-section overlaps ≤ 0.126, and meow's own spectrum
    there has no ties at all). This silently makes an `N_MODES`-mode FDTD
    interface basis rank-deficient by one per duplicate — worse than a
    simple reordering, since whatever mode the duplicate displaced is
    missing from the basis entirely, not just mis-indexed.

    Fixed via `FDTDSolver._iter_distinct_bands()`: for each requested
    `band_num` it checks the new eigenmode's transverse-E-field overlap
    (`beef.mode_identity.field_overlap_matrix`) against every already-
    accepted band in the same call, skips and retries on a near-1 match
    (threshold 0.9, up to 8 extra attempts), and raises if it can't find
    enough distinct bands. Wired into `section_mode_neff()` and
    `section_mode_fields()`. Confirmed inert on the coupler (identical
    results, runtime unchanged within noise) and confirmed to fix the real
    modal-multiplexer reproduction: re-running gap 1's interface-0 check
    with the fix applied, the guard caught and skipped exactly the two
    duplicate bands (band_num 11 and 15), and the resulting neff sequence
    has no exact ties and no more catastrophic η mismatches at those
    positions. 5 new meep-free unit tests in
    `tests/test_duplicate_band_guard.py` (mocked `sys.modules["meep"]`),
    including a negative case confirming genuinely different mode shapes
    with near-identical `k_x` are NOT flagged.

    `solve_section()`'s own production interface construction hits the same
    underlying meep behavior through a different API
    (`get_eigenmode_coefficients()` on a DFT-decayed flux after `sim.run()`,
    not the cheap `get_eigenmode()` the guard above wraps) and was
    deliberately left unchanged this pass — extending the guard there
    either needs an extra eigensolve per band on a different code path than
    the one actually producing the basis (easy to get subtly wrong), or a
    real time-domain run to verify it, and `solve_section()` is too heavily
    used and delicately tuned to risk on an under-verified change.

    **Follow-up done, same day** (`fdtd_vs_EME_mode_validation.md` item 8):
    `beef.mode_identity.transfer_matrix_duplicate_scan()`/
    `warn_on_transfer_duplicate()` — a post-hoc scan of `solve_section()`'s
    already-computed `T` for near-parallel excited-side columns, wired into
    both `FDTDSolver.solve_section()` and `EMESolver.solve_section()` right
    before they return. Zero extra simulation cost (`T` is already computed);
    logs a WARNING, does not raise — a high similarity is a symptom, not
    proof, unlike the field-shape same-engine controls. Not yet exercised
    against a real device's own duplicate-band reproduction (that needs an
    actual meep/meow run, not a code change).

    Not modal-multiplexer-specific: any device with closely-spaced
    higher-order bands (the coupler, either MMI) could hit this too via
    `solve_section()`, just not yet checked there.
13. **New, 2026-07-30 — `beef.mode_identity`: field-overlap mode
    correspondence.** A device- and resolution-agnostic replacement for
    nearest-neighbor n_eff matching (`tests/test_mode_identity.py`'s
    `assert_modes_correspond()`), which is fundamentally a 1-D scalar
    heuristic that picks the wrong pairing outright once several modes crowd
    together in n_eff (observed directly: bands 8, 9, 10 of the modal
    multiplexer's basis all nearest-matched to the same meow band).
    `field_overlap_matrix()` compares mode *shape* — a normalized transverse
    E-field overlap, deliberately **not** the full E-x-H reciprocity overlap,
    because meep and meow were found to report E and H in mutually
    inconsistent absolute units (~375x apart on a real solved mode, each
    internally self-consistent) — and `best_correspondence()` picks the
    overlap-maximizing assignment (Hungarian algorithm) instead of
    nearest-neighbor. `FDTDSolver.section_mode_fields()` (new alongside it)
    is the field-sampling counterpart of `section_mode_neff()`. **Now also
    wired into `HybridSolver` itself as an opt-in runtime check — see item
    14.** Not yet wired into `tests/test_mode_identity.py`'s per-device
    checks — still used there only via the (temporary)
    `examples/verification/python/audit_interface_modes.py` script.
14. **New, 2026-07-30 — `HybridSolver(check_mode_correspondence=True)`: a
    runtime version of gap 1/6's manual audit.** `HybridSolver.solve()`'s
    matrix-chain contraction assumes, with no check at all, that mode index
    *i* of one section's output basis is the same physical mode as index
    *i* of the next section's input basis. This wires item 13's field-overlap
    machinery into `HybridSolver` itself as an opt-in constructor flag
    (`check_mode_correspondence`, `mode_correspondence_action="raise"|"warn"`,
    `mode_correspondence_min_overlap`), so any device/chain catches this
    automatically at solve time instead of relying on someone running the
    audit script by hand. Duck-typed and best-effort by design (matching
    `calibrate_body_neff`'s existing opt-in convention): silently skipped
    (info log) when either adjoining solver lacks a field-sampling hook
    (`FDTDSolver.section_mode_fields()` or the new
    `EMESolver.section_modes()`) or a domain size; skipped with a WARNING
    (regardless of `mode_correspondence_action`) when the field query itself
    raises a geometry `ValueError` — deliberately not fatal, since that is
    "couldn't check here," not "checked and found a mismatch."

    **Real-device finding, not a hypothetical caveat**: every current
    hybrid example either couldn't actually run this check, or hit a
    pre-existing bug when it tried. Confirmed directly against each: the
    directional coupler and both MMI examples all use a `body_stub`/
    `BODY_STUB` of 0.25 µm at their interior interfaces, thinner than
    `FDTDSolver`'s own automatic `dpml=0.5` query-window padding — the check
    correctly logs a warning and skips there rather than crashing or
    mismeasuring. The modal multiplexer's `BODY_STUB=0.75` clears that bar,
    but its `CompoundGeometry`-composed geometry initially tripped the
    `CompoundGeometry` stub double-application bug (Geometry item 6 — the
    one `audit_interface_modes.py` had also routed around by querying
    `mmi1`/`mmi2` in their own local frames) — correctly caught and skipped
    with a warning, not a crash. **Geometry item 6 has since been fixed
    (same day)**, and the modal multiplexer's check now gets past it and
    runs its real 16-mode comparison. The comparison path was separately
    confirmed to run and pass cleanly on a real meep/meow FDTD<->EME
    interface via a synthetic uniform `WaveguideRun` section. **Remaining
    limitation, fixed 2026-07-30 (`fdtd_vs_EME_mode_validation.md` item 5)**:
    the coupler and both MMI examples couldn't run this check at all, purely
    because their 0.25 µm `body_stub` is thinner than `FDTDSolver`'s automatic
    `dpml=0.5` query padding. Rather than widening those example constants
    (which would have meant re-recording every affected golden for no physics
    reason), `HybridSolver` gained `mode_correspondence_pad`, threaded to
    `section_mode_fields()` for solvers whose signature declares a `pad`
    parameter — safe because the mode-field query never time-steps, so a
    smaller clearance than a real PML footprint is exactly as valid, just
    thinner. A value like 0.1 fits inside `body_stub = 0.25`. Does not itself
    turn `check_mode_correspondence=True` on for these three examples — that
    is a separate decision needing real-device verification, not done here.
    13 new meep/meow-free unit tests in `tests/test_hybrid_mode_correspondence.py`
    (mocked `sys.modules["meow"]`, synthetic sine-mode `FieldSet`s) plus 2 new
    real-`tidy3d` tests for `EMESolver.section_modes()` in
    `tests/test_eme_length_sweep.py`.

    **Self-inflicted bug found and fixed during verification, worth
    recording as a general lesson.** The diagnostic mesh was first written as
    a *fixed pixel count* (24 wide), reasoning by analogy from
    `FDTDSolver._iter_distinct_bands()`'s deliberately-coarse fixed grid.
    That analogy is invalid: `_iter_distinct_bands()` compares meep against
    *itself*, where two identical eigenvectors overlap to exactly 1.0 on any
    grid however coarse, so coarseness genuinely costs nothing there. A
    *cross-engine* comparison has no such immunity — and worse, an
    `EMESolver` side **solves** its modes on whatever mesh it is handed
    (`section_modes(mesh=...)`), so under-resolving does not merely sample
    the comparison coarsely, it makes one side's modes physically wrong. On
    a 5 µm-wide domain, 24 pixels is 4.8 px/µm, leaving SOI_160's 0.16 µm
    core **0.77 cells** tall. Measured on the modal multiplexer's interface
    0: worst diagonal overlap **0.0004** with modes 7-15 scrambled, versus
    0.973-0.998 and clean ordering for that same interface at 40 px/µm. The
    diagnostic giveaway was the *shape* of the failure — modes 0-6 matched
    identity and only 7+ scrambled, which is under-resolution (low-order
    modes' broad features survive coarse sampling; high-order ones don't),
    not a basis disagreement. Fixed: the mesh is now specified as a
    resolution (px/µm), defaulting to the coarser of the two adjoining
    solvers' own `res` and overridable via `mode_correspondence_res`, with 6
    new tests pinning square pixels, thin-core resolution, and the
    derivation/override/fallback logic.

    **Second self-inflicted bug, same verification pass: no guided-mode
    filter.** The design deliberately compared exactly `interface_modes[k]`
    modes with no cladding-index prefilter, on the reasoning that this is
    what `solve()` actually contracts. That was wrong. A mode below the
    cladding index is a *radiation* mode, not a discrete physical mode — its
    shape is set by the simulation box and each engine's own discretization,
    so two engines' radiation modes are essentially unrelated and overlap ~0
    however well the real basis corresponds. And per
    `validation_markdown/README.md`'s own "Finding 1", **every** interface
    basis in this repo carries at least one unguided mode — so the check
    would have false-alarmed on every device. Measured on the modal
    multiplexer at the corrected 40 px/µm: worst diagonal overlap 0.0000,
    driven entirely by bands 13-15 (n_eff 1.434/1.413/1.388 against a 1.4468
    cladding index), while the 13 guided bands were unaffected. Fixed: both
    sides are now filtered to `n_eff > MaterialStack.sio2_index(wl)` (read
    off either adjoining solver, skipped entirely if neither exposes a
    stack) before the overlap, with the excluded count logged. 4 more tests,
    including the important negative case that a genuine mismatch among
    *guided* modes is still caught — the filter must not become a blanket
    excuse. Note this does not make the unguided modes harmless: `solve()`
    still contracts them, which is Finding 1's own separate complaint and
    not something a correspondence check can fix.

    **Process note worth keeping.** Both of these were caught only by running
    the check against a real device and refusing to accept a
    confident-looking failure at face value. The first reported a crisp
    "correspondence failure"; the second reported worst-overlap 0.0000. Both
    were artifacts of the checker. In each case the diagnostic tell was the
    *shape* of the failure rather than its magnitude — modes 0-6 clean and
    7+ scrambled (under-resolution), then exactly the three sub-cladding
    bands at fault (radiation modes) — which is the kind of thing a summary
    statistic hides. A validation tool that has never been run against a
    known-good reference is an untested tool.
15. **New, 2026-07-30 — `section_mode_beta()` clamps for `solve_section()`'s
    `±dpml` draw padding in its seed path but not in its fit path. FIXED
    2026-07-30.**
    `fdtd.py`'s `seed_span = min(2.0, span - 2.0*self.dpml)` exists
    *specifically* because the mode solve pads its query by `dpml` on each
    side (its own comment says so). But the fit lengths a few lines below —
    `spacing = min(span/3.0, max_safe_spacing/2.0)`, `lengths = [spacing,
    2*spacing, 3*spacing]` — carry no such clamp, and each is handed to
    `self.solve_section(mid ± L/2, ...)`, which draws `stub + dpml`. Whenever
    `span/3` is the binding term (short bodies), `lengths[-1] == span`
    exactly, so the drawn window is `mid ± (span/2 + dpml)` — `dpml` past
    both ends of the region the caller declared uniform.

    Reproduced: at β≈9.6 (the value `fdtd.py`'s own docstring quotes),
    `max_safe_spacing = 2π/(0.25·β) = 2.618 µm`. A `DirectionalCoupler` with
    `coupling_length=3.0` (body span 2.50) gives `lengths=[0.833, 1.667,
    2.500]`, a drawn half-window of 1.750 against a uniform half of 1.250,
    and `section_prisms()` raises `Cut/stub at x=-1.7500 falls strictly
    inside an S-bend` — verified directly against the real geometry.
    Threshold is roughly any body span below ~1.5·`max_safe_spacing` (≈3.9 µm
    at that β) that is flanked by a bend/taper. Not reachable today:
    `calibrate_body_neff` defaults False, and both live callers use a ~29.5 µm
    body. Reachable by any user with a short uniform body, or
    `calibrate_body_neff=True` on a compact device.

    **Fix**: the default `spacing` computation now uses
    `usable_span = span - 2.0*self.dpml` in place of the raw `span`, mirroring
    `seed_span` exactly — `usable_span`'s own `<= 0` case needs no separate
    raise, since `seed_span`'s existing check above it already guards the
    identical `span - 2*dpml <= 0` condition and fires first. On the
    `coupling_length=3.0` repro this now gives `lengths=[0.5, 1.0, 1.5]`, a
    drawn half-window of exactly 1.250 — flush with, not past, the uniform
    half. 4 new tests in
    `tests/test_beta_calibration.py::TestSectionModeBetaDpmlClearance`
    exercise the real `section_mode_beta()` code path (`section_mode_neff()`/
    `solve_section()` swapped for numpy-only fakes, no meep needed): one
    documents the pre-fix overrun arithmetic, one confirms every default fit
    length now stays within the uniform span, one confirms the fit still
    recovers the seeded β, one confirms the pre-existing (unchanged)
    zero-clearance raise. Confirmed to fail on the unfixed code before the
    change (caught the exact 1.3333-vs-1.2500 overrun) and pass after.
    `TestTheSpacingPrecondition.test_the_actual_default_is_span_independent`
    (same file) updated to check the `dpml`-aware invariant instead of the
    old span-only one.
16. **New, 2026-07-30 — `tests/test_mode_identity.py::assert_modes_correspond()`
    re-scoped to a fast smoke check. DONE 2026-07-30.**
    (`fdtd_vs_EME_mode_validation.md` item 10.) Now that item 13's
    field-overlap machinery is the trusted mechanism, this older n_eff
    nearest-neighbour helper is strictly weaker (1-D scalar, resolution-
    locked, shape-blind — see its own module docstring's already-recorded
    caveats) and was never retired. Re-scoped rather than deleted: its fast
    logic-only tests and its two `@pytest.mark.slow` real-device smoke tests
    are still cheap and useful as a first line of defense, so kept, but its
    module docstring, its own docstring, and both real-device call sites now
    say plainly it is a smoke check, superseded by `beef.mode_identity`/
    `HybridSolver(check_mode_correspondence=True)` for anything needing an
    authoritative answer. `beef/mode_identity.py`'s own module docstring
    updated to point back the same way. Docs/docstrings only — no behavior
    change, existing fast suite (`tests/test_mode_identity.py -m "not
    slow"`) unchanged (5 passed, 2 deselected).
17. **New, 2026-07-30 — `_iter_distinct_bands()`'s duplicate check runs on a
    grid that cannot resolve the modes it judges. CHECKED 2026-07-30, NOT
    LIVE at res 20.** (`fdtd_vs_EME_mode_validation.md` item 11.) The guard
    decides "same eigenvector or not" from a fixed `_DUP_BAND_DIAG_NY = 16` ×
    `_DUP_BAND_DIAG_NZ = 8` sample spanning the whole **monitor** width. On the
    modal multiplexer's interface that monitor is 5.0 µm wide, so the
    transverse pitch is 0.3125 µm and only 8 of the 16 points land inside the
    2.8 µm core — against a ~0.19 µm lobe pitch for the 15th mode. The
    constant's justifying comment ("only needs to distinguish 0.9999 from
    ≤ 0.1257") was measured on the **6-mode** coupler and does not carry to 16
    by construction: idealized slab profiles `sin(mπy/W)` on that exact grid
    alias mode *m* onto mode *18−m* at η ≈ 0.9988, well past the 0.9
    threshold.

    **Measured directly on real fields**, `examples/verification/python/
    diagnose_dup_band_grid.py`: solved band_num 1–20 unguarded at the mux's
    interface, res 20 (a resolution job 11575516 confirmed the shipped guard
    fires at), then replayed the guard's accept/skip logic on the shipped
    grid and three finer ones (32×8, 64×16, 128×16 — up to 72/128 points
    inside the core vs the shipped 8/16). All four grids skip the identical
    two bands (11, 16), both independently confirmed as genuine duplicates by
    an exact n_eff tie (`|Δn_eff| < 1e-9`, no field sampling needed), and the
    worst overlap among accepted bands tops out at 0.34 (shipped grid) —
    comfortably under the 0.9 threshold, not approaching it. The idealized
    aliasing estimate does not transfer: real fields extend into the
    cladding rather than vanishing at the core edge.

    **Not re-checked at res 40** (the shipped `RESOLUTION`, where job
    11575516 also confirmed firings, at band_num 11/15) — a finer mesh packs
    more oscillation into the same core width, so this is worth one more data
    point before closing the item, but the res-20 result already moves it
    from "latent risk" to "checked and clear so far."
18. **New, 2026-07-30 — every band's eigensolve starts from the same k guess,
    and it is the likely root cause of item 12. TESTED 2026-07-30, FIXED
    2026-07-31 — `_iter_distinct_bands()` now ladders the seed (policy B,
    the naive ladder — measured byte-identical to the fancier C/D variants,
    so shipped as-is), landed together with item 19's cladding stop since
    laddering alone would still be exposed to the same hang. 6 new tests,
    `tests/test_duplicate_band_guard.py::TestSeedLaddering`/
    `TestGuidedEdgeStop`. See `fdtd_vs_EME_mode_validation.md` items 12/14
    for the full record.** (`fdtd_vs_EME_mode_validation.md` item 12.)
    `_iter_distinct_bands()` calls `sim.get_eigenmode(..., band_num,
    mp.Vector3(f * _KPOINT_SEED_NEFF, 0, 0))` with `_KPOINT_SEED_NEFF = 2.5` —
    the same starting vector for band 1 and band 18 alike, *N* independent
    targeted solves from one seed. Two neighbouring targeted solves launched
    from an identical guess converging to the same eigenpair is the textbook
    shape of the duplicate item 12 works around, and it explains the engine
    asymmetry: `EMESolver.compute_modes()` goes through
    `compute_modes_meep(cross_section, num_modes=N+8)`, **one simultaneous
    solve** for the lowest N+8 bands.

    That asymmetry is measurable, not just structural. Job 11575516 found that
    at the mux's interface 0 the simultaneous side finds both members of a
    near-degenerate pair at every resolution 20–80 while the targeted ladder
    finds one — the single disagreement a finer mesh does not close
    (`fdtd_vs_EME_mode_validation.md` item 1).

    **Tested**, `examples/verification/python/diagnose_kpoint_ladder.py`: 20
    bands unguarded at the mux's interface, res 40, fixed seed (today's code)
    vs. laddered (band *k* seeded from band *k−1*'s own converged n_eff).
    Fixed seed: 3 exact ties, 17/20 distinct. Laddered: 1 exact tie, 19/20
    distinct — laddering **recovers 2 genuinely distinct modes** the fixed
    seed silently drops entirely (not just reorders). But it also reproduces
    the anticipated failure mode: seeding band *k* directly at band
    *(k−1)*'s answer can make it **re-find** band *(k−1)*'s own eigenvector,
    and did — a new duplicate appears at bands 9/10 that the fixed run did not
    have there. Net win, not a closed fix.

    Unplanned side finding: this run (20 bands, one more than earlier work's
    18) turned up a **third** raw duplicate at band_num 19/20 (n_eff
    1.349097) beyond the two `modal_multiplexer_validation.md` gap 6 already
    documents — within both guards' extra-attempt headroom, so not an
    exhaustion risk, but the guards are doing more work on this device than
    previously measured.

    Not yet tried: seed band *k* from *(k−1)*'s neff **minus a small step**
    (the basis is monotonically descending, so this searches past the
    previous band instead of landing on it) or from the last **accepted**
    band rather than the raw previous solve — either is a plausible fix for
    the self-reconvergence failure just demonstrated. (b), one simultaneous
    `compute_modes_meep()` call as a count reference, not yet run. **Attempted
    2026-07-31, job 11586046 — it HUNG and was cancelled; retrying as
    11588309.** The attempt is worth recording because it changes the framing:
    asking for 20 distinct modes walked `get_eigenmode()` up to band_num 23,
    where **MPB's k-search has no root and limit-cycles forever** (161 solves on
    one band, 59% of the job's eigensolve work, against 3–8 for every band below
    it). The re-run (11588309) then **answered the seeding question outright, and
    in laddering's favour**: policy A (today's fixed seed) finds 13 distinct
    guided modes at this cross-section, while all three laddered policies find
    **14** — the extra being n_eff 1.512589, absent from A entirely rather than
    reordered. Same 16 attempts, one fewer skip, so laddering is a strict
    improvement and the fix direction is settled: seed band *k* from band *k−1*.
    C's step and D's accepted-only refinement buy nothing measurable over the
    naive ladder here; do not adopt them on this evidence.

    A shortfall does remain — 14 of 16 requested — and the guard makes up the
    difference by reaching past the guided edge and **padding with unguided
    modes**. That half is item 19 below and
    `fdtd_vs_EME_mode_validation.md` item 14. (An earlier revision of this entry
    said only 13 distinct guided modes *exist* and that no seed policy could
    recover them. That was inferred from the fixed-seed trace — i.e. from the
    enumeration under test — and the re-run falsified it.)

    Related, smaller — **now measured and closed, 2026-07-31, job 11586098**:
    `guard_duplicate_modes=True` makes `EMESolver.compute_modes()` over-request
    `num_modes + 8` in a single call, so guarded and unguarded hand MPB
    **different eigenproblems** — "guard on vs off" is not a clean A/B, and
    inertness had only been checked at 3 modes/res 12. Measured at 16 modes/res
    60 with three arms (guard off/16, guard on/16, guard off/24): the
    over-request is **exactly inert** (max |Δn_eff| = 0.000e+00 on all 16 low
    bands), so the guard is a pure filter and no guarded EME number carries an
    over-request artifact. It costs **3.2×** (934 s → 2948 s), which makes
    `_DUP_MODE_MAX_EXTRA = 8` pure cost against skip counts of ≤ 3 — deliberately
    left at 8 until the skip evidence covers more than one device and item 12
    settles whether a k-seeding fix retires the guard entirely.
    (`fdtd_vs_EME_mode_validation.md` item 13.)

19. **New, 2026-07-31 — `_iter_distinct_bands()` has no effective upper bound:
    it can hang a production run forever, and short of that it pads a basis
    with unguided modes.** (`fdtd_vs_EME_mode_validation.md` item 14.) Found by
    job 11586046 hanging on the diagnostic version of the same loop.

    `sim.get_eigenmode()` asked for a `band_num` with no solution at the target
    frequency does not raise — MPB's k-search limit-cycles between two k values
    indefinitely (measured: k = 0.2346 ↔ 1.7480 either side of f = 0.7634, 161
    solves before SLURM killed it). `_iter_distinct_bands()`
    (`beef/solvers/fdtd.py:1049`) walks `band_num` up to
    `n_modes + _DUP_BAND_MAX_EXTRA` and its `RuntimeError` can only fire *after*
    a solve returns, so that bound never engages: the loop hangs inside meep
    with no exception and no timeout. A Python alarm cannot help — MPB is C, so
    the signal handler does not run until the call returns.

    The same walk has a quieter failure below the hang. At the mux's interface 0,
    res 40 at most **14 distinct guided modes are reachable** — 15 guided band
    slots, with duplicate solves inside them (10/11 and 14/15 under a fixed seed,
    10 alone under a laddered one), and bands 16–22 below the cladding index. The
    device requests 16, so the guard reaches past band 15 and **accepts unguided
    modes** to make up the count — two or three of that basis's channels are not
    guided modes of the cross-section, depending on the seed policy. That is item 9's concern arriving
    through the guard rather than through a device's `N_MODES`, and it is why the
    mux's `RESOLUTION`/`N_MODES` pair has to be chosen jointly
    (`modal_multiplexer_validation.md` gap 2).

    **Fix, when taken:** stop the walk when a raw solve returns `n_eff ≤` the
    cladding index — nothing above it is guided and the next band_num is the
    non-convergent one, so one physical test closes both failures. Decide
    raise-vs-pad at the same time: if a cross-section cannot supply `n_modes`
    guided modes, failing loudly beats returning a basis whose rank is a lie.
    `diagnose_kpoint_ladder_v2.py` already does exactly this (job 11588309), so
    the behaviour has a working reference; the library change is deliberately
    deferred because it alters which modes a basis *contains*, and gap 6's scope
    argument for not touching this method under-verified still applies.

    **FIXED 2026-07-31.** Both halves landed together: the walk now stops
    the instant a raw solve's n_eff crosses the cladding index (computed
    once per caller via `self.material_stack.sio2_index(wl)`), and raises
    immediately if fewer than `n_modes` were accepted by then — the
    raise-vs-pad decision above resolved in favor of raising. 3 new tests,
    `tests/test_duplicate_band_guard.py::TestGuidedEdgeStop`. The deferral
    reasoning above was correct at the time (this needed its own tests plus
    a real before/after, not a blind patch) and both now exist — see
    `fdtd_vs_EME_mode_validation.md` item 14 for the full record, including
    the direct consequence: the modal multiplexer's own shipped `N_MODES=16`
    now raises on the next resubmission (14 reachable at res 40), which is
    the fix working as designed, not a new bug.

### Geometry

1. **Ring/racetrack resonator.** *(2026-07-29 — (a), (b) and (c) below now all
   exist as standalone, tested library capabilities; none are wired together or
   to this device yet — see Solvers/Geometry/Network above and
   `validation_markdown/ring_resonator_validation.md`'s recommendation section
   for the full record.)* The single most common photonics building block that
   isn't in the device library at all. Natural fit for the existing
   `DeviceGeometry` interface (bus waveguide + coupled ring, `section_prisms`/
   `width_at`/`center_at` following the `DirectionalCoupler` pattern), and directly
   exercises `HybridSolver`'s generic cascading (bus straights as FDTD/EME + coupling
   region). `examples/ring_resonator/python/ring_resonator_hybrid.py` proves out the
   physics (closed-form all-pass feedback assembly, not a `HybridSolver` chain — a
   ring is a loop, not a cascade) entirely in raw meep/meow, bypassing the geometry
   layer and every solver class — that script is still unchanged and still doesn't
   use any of what follows. What was missing for a first-class version:
   (a) an actual `RingResonator`/racetrack `DeviceGeometry` — **done**,
   `RingResonatorCouplerSection`/`RingResonatorBendSection`
   (`beef/geometry/ring_coupler.py`/`ring_bend.py`); (b) y-normal (x-z
   plane) source/monitor support in `FDTDSolver` — it only ever built x-normal
   (y-z plane) sources/monitors before, and the two bus-to-ring coupling ports in
   a racetrack coupler section are perpendicular to that — **done**,
   `solve_section()`'s new `y0`/`y1`/`excite_port`/`measure_port`
   (`PortSpec`) parameters, additive and default-`None`; (c) a feedback-loop
   assembly mode in `HybridSolver` (or a sibling), since its `solve()` is
   fundamentally a one-pass matrix chain, not a resonator sum — **done**,
   `beef.network.close_loop()`, verified to reproduce this device's own proven
   closed-form formula term for term. What remains is wiring (a)+(b)+(c) together
   into an actual device and a new example script — not designed here, see the
   validation doc.
2. **GDS import (read a real/foundry/gdsfactory layout) and export (tape-out-ready
   output).** The gdsfactory dependency was deliberately removed from the internal
   pipeline (raw polygon dicts are faster/simpler for BEEF's own parametric shapes),
   but that leaves no way to simulate an arbitrary externally-supplied layout or
   export a BEEF geometry as a real `.gds`. A thin, optional adapter (gdsfactory
   `Component` ↔ BEEF polygon dict) would restore that interop without
   reintroducing gdsfactory as a hard dependency — every one of gdsfactory+'s 43+
   foundry PDKs is a real device BEEF currently has no way to even load, let alone
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
6. **`CompoundGeometry.section_prisms()`/`section_gds()` double-applied
   `stub_left`/`stub_right`. FIXED 2026-07-30.** `_overlaps()` already bakes
   the stub-extended boundary into `ov_lo`/`ov_hi` (`lo, hi = x0-stub_left,
   x1+stub_right`, then clips each segment's overlap against it), but it also
   handed the raw stub values down to the first/last overlapping segment,
   whose own `section_prisms()`/`section_gds()` applied them a *second* time.
   Whenever the padded window landed fully inside one sub-segment — so
   nothing clipped the expansion back — the drawn geometry came out 2x the
   stub too long, or the doubled coordinate raised a spurious "cut falls
   strictly inside a taper/bend".

   **This WAS reachable in production, contrary to the first assessment
   written here.** `solve_section()` itself never triggers it (interior
   boundaries always pass `stub_left=stub_right=0.0`, per `hybrid.py`'s
   `_stubs()`) — but the *plotter* does: `beef/plot/geometry.py`'s
   `_draw_section_fdtd_markers()` path sets `draw_stub_left = stub_left +
   dpml` for any FDTD section so the drawn waveguide runs through the PML,
   which for an **interior** FDTD section is `0.0 + dpml` — a nonzero stub
   at an interior cut, precisely the buggy case. Measured on the modal
   multiplexer's own transition section (`x=[-12.910, -4.410]`, interior on
   both sides): the section plot drew a 10.5 µm span where 9.5 µm is
   correct, i.e. **0.5 µm of extra MMI body silently added at each end** of
   every interior-FDTD sub-region plot that example produces. Cosmetic
   rather than numerical (plots only, no solver result was affected), but it
   means the "unreachable in production" framing was wrong, and it is worth
   remembering that the plot path applies stubs the solver path does not.

   Also hit by two non-production callers: `examples/verification/python/
   audit_interface_modes.py` (routed around it by querying the modal
   multiplexer's `mmi1`/`mmi2` in their own local frames) and
   `HybridSolver(check_mode_correspondence=True)` (Solvers item 14 — its
   field-sampling query inherits `FDTDSolver`'s automatic `dpml` padding,
   which *is* a nonzero stub against the shared `CompoundGeometry`).

   **Fix**: `_overlaps()` now returns only the part of each stub that falls
   *outside* its segment (`max(0.0, seg_lo - lo)` / `max(0.0, hi - seg_hi)`),
   so the stub is applied exactly once — by `ov_lo`/`ov_hi` where it stays
   inside the segment, by the delegate where it genuinely extends past the
   edge. Both call sites' existing offset arithmetic already assumed the
   resulting invariant (`ov_lo - stub_left == lo` for the first segment,
   `ov_hi + stub_right == hi` for the last), so neither needed changing.
   4 new tests in `tests/test_geometry.py::
   TestCompoundGeometryStubIsAppliedExactlyOnce` pin the extent exactly in
   both the fully-inside and the reaching-outside case, plus the stub=0 path
   production depends on. Confirmed against the real device: the modal
   multiplexer's correspondence check now gets past the geometry error and
   runs its actual 16-mode comparison.
7. **New, 2026-07-30 — `CompoundGeometry` assumes every segment's local
   frame is centred on its own bbox; `MMI`'s is not. FIXED 2026-07-30.**
   `_placements()` allots
   each segment a slot `[seg_lo, seg_hi]` and every consumer maps device→local
   as `x - center` with `center = (seg_lo+seg_hi)/2` (`polygons()`,
   `section_prisms()`, `section_gds()`, `width_at()`, `center_at()`,
   `source_aperture_at()`, and `_end_ports()`'s `half = seg.total_length/2`).
   That is only valid when the segment's own frame spans
   `[-total_length/2, +total_length/2]`. `MMI`'s frame is centred on the
   **body**, not the bbox — it spans `[-mm/2-in_max, +mm/2+out_max]` — so
   whenever `max(input taper lengths) != max(output taper lengths)` every
   centre-based path is off by `(out_max - in_max)/2`. Note `to_gds_eme()`
   is *correct* here (it normalises each delegate to left-edge-0 and shifts
   by `seg_lo`), so two paths in the same class disagree.

   Reproduced directly (MMI `mm_length=5`, input taper 3.0, output tapers
   1.0, `total_length=9`, followed by a `WaveguideRun(2.0)`): the MMI's own
   frame spans `[-5.50, 3.50]` where bbox-centred would be `[-4.50, 4.50]`;
   `compound.polygons()` then spans `[-6.50, 5.50]` against a declared
   `[-5.5, 5.5]`, hanging 1 µm past its own stated left edge, and
   `compound.width_at(-total/2)` returns 0.467 (mid-taper) instead of the
   true 0.35 access width. Not reachable in-tree today — every current MMI
   uses one shared `taper_length` across all ports, and `Taper`/`SBend`/
   `DirectionalCoupler`/`GdsGeometry` are all genuinely bbox-centred — but
   per-port `taper_length` is a supported, documented `MMIPort` field that
   `MMI` itself handles correctly, so this is a valid device
   `CompoundGeometry` silently mis-composes.

   **Fix**: `DeviceGeometry` gained a `frame_offset` property (default
   `0.0` — bbox-centred, correct for every existing geometry except `MMI`)
   giving the device-frame offset from a segment's bbox midpoint to its own
   local-frame origin. `MMI.frame_offset` returns `(in_max - out_max) / 2`.
   Every `CompoundGeometry` device→local conversion (`polygons()`,
   `section_prisms()`, `section_gds()`, `width_at()`, `center_at()`,
   `source_aperture_at()`, `_end_ports()`) now adds `seg.frame_offset` to the
   bbox-midpoint `center` it already computed, instead of assuming the
   midpoint is the origin. `to_gds_eme()` needed no change (already correct,
   per above). 4 new tests in
   `tests/test_geometry.py::TestCompoundGeometryMMIFrameOffset`, reproducing
   this entry's own repro case: confirm `MMI.frame_offset == 1.0` for the
   `mm_length=5`/taper `3.0`/`1.0` device, and that `polygons()`/
   `section_gds()` stay within the declared bbox and `width_at()` at the true
   left edge returns the access-waveguide width (0.35) rather than the old
   mid-taper value (0.467).

8. **New, 2026-07-31 — `CompoundGeometry` cannot back `FDTDSolver.solve()`'s
   single-domain path at all. One part fixed, one part real and open.**
   Found writing `modal_multiplexer_fdtd.py` (`modal_multiplexer_validation.md`
   gap 3) — the first script in this repo to call `.solve()` (as opposed to
   `.solve_section()`, the only thing `HybridSolver` ever calls) with a
   `CompoundGeometry` as the top-level device. Every existing full-device FDTD
   script (`mmi_1x2_fdtd.py`, `mmi_2x2_fdtd.py`) passes a plain `MMI`, so this
   combination had literally never run before.

   **Part 1, FIXED — `source_aperture_at()`'s bounds check.**
   `FDTDSolver.solve()` queries `geometry.source_aperture_at(x)` at the
   solver's external PML-clearance stub position, `waveguide_length` past the
   device's own physical extent. Every leaf geometry already tolerates
   arbitrary overshoot there by construction — `MMI.width_at()` clamps to the
   terminal `waveguide_width` for any x past its taper start, however far
   past, and `MMI.source_aperture_at()` is built on top of it — but
   `CompoundGeometry._locate()` had no segment "out there" to delegate to and
   raised instead. `_locate(x, tol, clamp=False)` gained an opt-in
   `clamp=True` that delegates out-of-bounds x to the nearest end segment,
   wired into `source_aperture_at()` only — `width_at()`/`center_at()`/
   `section_prisms()`/etc. keep raising by default, since an out-of-bounds cut
   *there* is a genuine caller error building physical geometry, not a
   solver's legitimate external-stub query. 5 new tests,
   `tests/test_geometry.py::TestCompoundGeometrySourceApertureClamps`,
   mutation-checked against the pre-fix code.

   **Part 2, NOT FIXED — `polygons(stub_length>0)`.** Clearing part 1 exposed
   a second, unconditional crash: `solve()` → `to_meep_prisms()` →
   `geometry.polygons(stub_length=self.waveguide_length + self.dpml)`, and
   `CompoundGeometry.polygons(stub_length>0)` raises `NotImplementedError` by
   design — its own docstring says "no generic way to know per-port stub
   width for an arbitrary segment type". This one is real, unfinished work,
   not a mechanical patch: `to_gds_eme()`/`section_gds()` already solve the
   equivalent problem via `_end_ports()`, but generalising that to plain
   `polygons()` needs an actual decision for an ambiguous end (e.g. a
   multi-port `WaveguideRun` — which port's width does the external stub
   inherit?), which whoever wrote the existing `NotImplementedError` declined
   to guess at. Left alone, same caution bar as Solver items 14 and 12's
   guard fix this session: not touching production-adjacent, under-verified
   geometry logic without a real design pass plus device-level verification.
   `modal_multiplexer_fdtd.py` and its batch script are marked BLOCKED
   pending this.

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
   (a pure meow/MPB even/odd-supermode script that doesn't use any BEEF solver
   class at all) — worth generalizing into its own example once the Misc convergence
   utility exists.
4. **Ring resonator example paired with a real `RingResonator` `DeviceGeometry`.**
   `examples/ring_resonator/python/ring_resonator_hybrid.py` proves out the physics
   (see Geometry recommendation #1) and calls `SMatrix.check_passive()` as standard
   practice, but doesn't use the convergence-utility pattern above (no resolution
   sweep) and can't use `plot_fdtd_geometry`/`plot_hybrid_geometry` at all (both
   require a `DeviceGeometry`) — it hand-draws its own geometry preview plots
   instead. The geometry-class-backed version (proper `plot_*_geometry`
   integration, resolution-sweep convergence check) remains open. *(2026-07-29
   — the `DeviceGeometry` to pair it with now exists, see Geometry
   recommendation #1; this item is specifically about updating the example
   itself, still unstarted.)*

### Misc

1. **Mesh-convergence-loop utility.** The pattern hand-rolled in
   `directional_coupler_mode.py`/`.ipynb` (loop over a resolution list, save/plot
   each into its own subfolder) should be a reusable `beef` utility (e.g.
   `sweep_resolution(run_fn, resolutions) -> convergence data`) rather than
   hand-rolled per example.
2. **Result serialization beyond raw pickle.** `SMatrix.save`/`load` use unversioned
   `pickle` — no schema, no forward/backward-compatibility story, and no interop
   format. Touchstone `.sNp` is the standard RF/photonics interop format and isn't
   supported; a JSON or HDF5 option would also let results be inspected without
   unpickling BEEF's own classes.
3. **Declarative simulation runner.** Every example is a bespoke Python script; a
   YAML/JSON "device + solver + sweep" spec (in the spirit of tidy3d's
   JSON-serializable `Simulation` objects) would let standard sweeps run without
   writing Python, and would pair naturally with `ParallelSweep`.
4. **Output-folder conventions.** `_utils.py` already has `resolve_output_folder`/
   `make_output_folders`, but usage is inconsistent across examples (most
   `examples/*.py` scripts hand-roll their own `os.makedirs`/`output_dir` logic
   instead) — worth consolidating on the existing utility everywhere instead of ad
   hoc per-script path-building.
5. **Array-task sweep harness.** *(2026-07-31 — the five drivers this item
   named moved to `examples/directional_coupler/archive/`, retired pending
   new sweeping functionality once every device is verified; this
   recommendation is unaffected and applies to whatever replaces them.)* Five
   sweep drivers existed under `examples/directional_coupler/python/` —
   `directional_coupler_fdtd_sweep.py`, `..._fdtd_resolution_sweep.py`,
   `..._hybrid_modes_sweep.py`, `..._hybrid_resolution_sweep.py`, and
   `..._hybrid_stub_sweep.py` — and every one hand-copies the same harness:
   a `_tag()` filename scheme, `_atomic_replace_write()`, a
   "SLURM array index selects one swept value, no index runs the whole sweep
   serially" CLI, per-value `data/*.npz` + `plots/*.png` writes, aggregation
   that re-reads whatever landed so far, and an "n/N done" annotation for
   partially-drained arrays. Only the physics in the middle differs. A
   `beef.sweep` harness taking (swept values, a per-value run function
   returning arrays, an aggregate-plot hook) would leave each example holding
   just its physics — and their batch scripts in `sweep/` are equally
   near-identical, differing only in make target and array range. Distinct from
   Examples #2 (`ParallelSweep` is a local `multiprocessing` grid over
   whole-device `solve()`, not a cluster array) and broader than Misc #1
   (resolution-specific).
6. **Pair notebooks to scripts instead of mirroring them by hand.** Twelve
   examples ship as a `.py` plus a hand-maintained `.ipynb` copy of the same
   code (`directional_coupler_hybrid`, `mmi_2x2_hybrid`,
   `modal_multiplexer_hybrid`, `ring_resonator_hybrid_cascaded`, …), with the
   notebook differing only in `CACHE_ONLY = True`, ASCII-ised units, and
   `plt.show()` in place of `savefig`. Nothing enforces that they stay in sync:
   the 2026-07-24 directional-coupler interface change had to be applied twice
   and verified by executing the notebook's setup cells by hand, and a missed
   mirror would leave a published tutorial silently teaching the old
   convention. jupytext pairing (or generating the notebooks from the scripts)
   plus a CI check that the pair is in sync removes the whole class of drift.

   **The predicted drift happened, and it was found by audit rather than by
   anything failing (2026-07-31).** Diffing every notebook's constants against
   its script turned up four mismatches, three of them from one pass:

   | file | notebook | script | |
   |---|---|---|---|
   | `directional_coupler_hybrid` | `N_MODES = 6` | **5** | item 9 missed it |
   | `mmi_1x2_hybrid` | `N_MODES = 12` | **8** | item 9 missed it |
   | `mmi_2x2_hybrid` | `N_MODES = 15` | **14** | item 9 missed it |
   | `directional_coupler_fdtd` | `COUPLING_LENGTH = 15.13` | 15.15 | **script** was wrong |

   `fdtd_vs_EME_mode_validation.md` item 9 called itself a cross-device pass and
   updated every script; no notebook moved. The 1×2's stale 12 is the worst of
   them — that is the exact configuration `mmi_1x2_validation.md` gap 3 identified
   as *the* failure mode (degenerate bands piling up at cutoff on SOI_160), so
   running that tutorial reproduced a known-bad basis.

   The fourth is the reverse and the more interesting one: the notebook was
   right. `COUPLING_LENGTH` was deliberately set to 15.13 in four places with
   documented reasoning (item 5's ground-truth crossing 15.129 µm, not the
   hybrid's grid-snapped 15.151, because three of the four users are full-device
   FDTD runs — `archive/dc_validation_20260727.md` item 14). Commit `e85b350`,
   whose message is *"Remove cached example simulation files for ring resonator
   hybrid…"*, moved `directional_coupler_fdtd.py` to 15.15 without mentioning it.
   A physics constant in a reference run drifted inside a cache-cleanup commit
   and nothing noticed for two days. All four are corrected; the case for
   generating notebooks from scripts, and for keeping unrelated changes out of
   cleanup commits, is now empirical rather than hypothetical.

## Not recommending

- Building BEEF's own raw FDTD/EME numerical core — meep and MPB/meow already do
  this well; BEEF's value is the orchestration layer on top, not the field solver.
- Matching tidy3d's cloud/GPU infrastructure — out of scope for what is currently a
  single-machine/HPC-Makefile-driven library.
