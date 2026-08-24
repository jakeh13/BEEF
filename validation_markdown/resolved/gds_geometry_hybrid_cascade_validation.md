# GDS-file geometry in a hybrid cascade — build plan

Opened **2026-08-20**. **CLOSED 2026-08-21.** `GdsGeometry` could back a
whole-device single-solver run, but could not occupy a section of a
`HybridSolver` cascade. Steps 1-5 are now all **implemented and MEASURED**
against a GDS file cut to the 2×2 MMI's exact dimensions, compared directly to
the parametric `MMI` run.

Step 5's production-scale acceptance run (job **12069544**, matching
`mmi_2x2_hybrid.py`'s `RESOLUTION=40`, `N_WL=11`, 201-point length sweep)
**COMPLETED 2026-08-21 00:06:32 in 4:31:02, exit 0:0**. It agrees with the
parametric run to **1.6e-3 relative** — looser than the "solver noise" this plan
expected, and the whole of that residual has been **traced to a measured root
cause that is not a defect in the GDS path**: the two backends hand meep the same
five polygons in a *different order*, and meep's subpixel averaging resolves
abutting-prism seams order-dependently. Reordering makes the discretized
permittivity **bit-identical**. See step 5's "Root cause" block.

> Every step's tolerances/behavior are **MEASURED**, not argued — see each step's
> own "Measured" block. Nothing in this document is left ARGUED.

---

## Summary

| | status |
|---|---|
| Read a GDS file into BEEF's polygon convention | ✅ `load_gds_polygons()` (`geometry/gds_import.py:19`), via gdstk, µm coordinates |
| GDS geometry in a whole-device **EME** run | ✅ `GdsGeometry.to_gds_eme()` |
| GDS geometry in a whole-device **FDTD** run | ✅ `polygons()` → `to_meep_prisms()`; used by the HCG reflector example |
| GDS geometry in a **hybrid cascade** | ✅ **implemented 2026-08-20** — `section_prisms`, `section_gds`, `width_at`, `center_at`, `source_aperture_at` all now on `GdsGeometry` |
| Detecting a varying cross-section without parametric knowledge | ✅ `GdsGeometry._check_cut`, reusing `assert_uniform_cross_section` (`geometry/base.py`) |
| Validation vehicle (twin GDS) | ✅ `examples/gds_files/mmi_2x2_twin.gds` + `make_mmi_2x2_twin.py`, round-trips to 2.2e-16 µm |
| Example: 2×2 hybrid cascade from a GDS file | ✅ `examples/mmi_2x2/python/mmi_2x2_hybrid_gds.py` |
| Regression tests | ✅ `tests/test_gds_import.py` (44 tests, all passing, meep/meow-free) |
| Full-scale solver acceptance run | ✅ **job 12069544 COMPLETED** 2026-08-21, 4:31:02 — agrees to 1.6e-3 rel, residual traced to meep prism-order sensitivity, not the GDS path |

Ordered work: **validation vehicle (step 1) → `width_at`/`center_at` (step 2) →
`section_prisms`/`section_gds` (step 3) → cut validity (step 4) → acceptance
(step 5).** **All five done.**

---

## Why this does not work today

`GdsGeometry` (`geometry/gds_import.py:72`) implements `polygons()` and
`to_gds_eme()` only. The three hooks `HybridSolver` requires are inherited
unimplemented from `DeviceGeometry`, whose base methods raise with exactly this
diagnosis (`geometry/base.py:155-192`):

```
raise NotImplementedError(
    f"{type(self).__name__} has not implemented section_prisms(); "
    "required for use with HybridSolver."
)
```

and likewise for `section_gds()` and `width_at()`. `center_at()` defaults to
`0.0` (`base.py:194-196`), which is wrong for any off-axis GDS device but does
not raise — a quieter failure mode worth noting.

So a `GdsGeometry` handed to `HybridSolver` fails at the first section solve, not
at construction.

*source: `beef/geometry/gds_import.py:72-120`; `beef/geometry/base.py:150-200`*

---

## Step 1 — Validation vehicle: a GDS twin of the 2×2 MMI

Write a GDS file whose core-layer polygons are geometrically identical to what
`MMI(...)` draws for the shipped 2×2 parameters (`mmi_2x2_hybrid.py:88-98`):
`n_inputs=2, n_outputs=2, mm_length=28.0, mm_width=2.8, waveguide_width=0.45,
taper_length=4.0, taper_width=1.3, waveguide_offset=1.5, waveguide_height=0.22`.

**Why this vehicle.** It makes the whole exercise a controlled comparison: same
physics, same solvers, same resolution, same material stack, same cut positions —
only the geometry backend differs. Any disagreement is attributable to the GDS
path rather than to the device. The 2×2 also has the two properties that stress
the implementation hardest: **close-set ports** (1.5 µm centre-to-centre, which is
why both FDTD sections use `half_domain_aperture`) and **off-axis branches**,
which is what makes `center_at()`'s `0.0` default insufficient.

Generate it from `MMI.polygons()` itself via gdstk rather than drawing it by
hand — then the twin is exact by construction, and the comparison isolates the
*query* path (`width_at`, `section_*`) rather than accidentally testing whether
two independently drawn layouts match. Commit it under
`examples/gds_files/` alongside `HCG_seeded_TO.gds`, plus the generating script so
it is reproducible.

**MEASURED (2026-08-20).** `examples/gds_files/make_mmi_2x2_twin.py` builds
`mmi_2x2_twin.gds` from `MMI(**MMI_PARAMS).polygons(stub_length=0.0)` (the
`mmi_2x2_hybrid.py` parameters, kept in one place and asserted against at
example run time — see step 5) and round-trips it back through
`load_gds_polygons()` to check every vertex: **worst-case round-trip error
2.22e-16 µm** (float noise, not truncation — gdstk's default 1 nm precision
grid divides every 2×2 vertex exactly). File holds the *bare* device
(x ∈ [-18, +18], no access waveguides) — deliberately, so
`access_extension` (step 3) is the thing under test for the cascade's end
stubs, not something the file quietly supplies.

---

## Step 2 — `width_at` / `center_at` from polygons

Both must answer "what is the cross-section at device-frame x?" from geometry
alone. For a straight-through device the reading is: intersect the core polygons
with the vertical line at `x`, collect the resulting y-intervals.

Design points:

- **Multiple branches.** At the 2×2's input, the line at `x` crosses **two**
  waveguides. `width_at`'s contract is *"Local waveguide cross-section width at
  device-frame position x"* (`base.py:187`) — singular. `MMI` resolves this with
  `input_idx`/`output_idx` kwargs selecting which branch to report
  (`mmi.py:205-210`: *"those two only select which branch width_at()/center_at()
  report for monitor sizing"*). `GdsGeometry` must accept the same kwargs and
  map them to y-intervals ordered by y, so the same `geometry_kwargs` the 2×2
  already passes work unchanged.
- **`center_at`** returns that interval's midpoint — not `0.0`.
- **`source_aperture_at`** (`base.py:198`) should report the widest window centred
  on the selected branch that stays clear of its neighbours. The polygon
  intersection already yields the neighbouring intervals, so the gap is directly
  measurable — arguably *better* conditioned than the parametric path, which
  computes it from known port spacing.
- **Tolerance.** Polygon vertices are floats from the file; use a small epsilon
  when deciding whether the line clips a vertex exactly, and prefer sampling just
  inside the span over exactly at a boundary.

**Done when:** `width_at`/`center_at`/`source_aperture_at` on the step-1 twin match
`MMI`'s to within a stated epsilon at a dense set of x, including inside taper
bands where the width varies continuously. Pure geometry, so this belongs in
`tests/test_gds_import.py`, which already writes a synthetic 2-layer GDS with
gdstk and imports neither meep nor meow.

**MEASURED (2026-08-20).** Implemented in `GdsGeometry._branch_at`/`width_at`/
`center_at`/`source_aperture_at` (`geometry/gds_import.py`): a per-layer,
per-x scanline (`_scan`/`_poly_crossings`) intersects every core polygon with
the vertical line at `x`, merges touching intervals across layers/polygons
(a body and the taper that abuts it are drawn as separate touching
polygons — reporting them unmerged would answer 1.3 µm at a 2.8 µm body
face), and indexes the result the same way `MMI` does: `input_idx` selects
among branches left of `split_x`, `output_idx` among those right of it.

Swept 1401 points from -20 to +20 µm (access-waveguide region, both taper
bands, the body) at 3 `input_idx`/`output_idx` combinations:
**max |Δwidth_at| = 4.44e-16 µm, max |Δcenter_at| = 1.11e-16 µm** — float
noise, not a real discrepancy. `source_aperture_at` (opt-in
`neighbor_aperture=True`, mirroring `MMI`'s `use_neighbor_aperture`) agreed
to **8.88e-16 µm** against `MMI`'s port-spacing-derived answer over the same
sweep, confirming the plan's claim that the polygon-measured aperture is at
least as well conditioned. `center_at`'s old silent `0.0` default is gone —
`test_center_at_is_not_the_base_class_zero` pins the real ±0.75 µm branch
centres. 12 tests in `TestSectionHooks` (`tests/test_gds_import.py`), all
passing, meep/meow-free.

One real tolerance pitfall found while writing the tests, not the
implementation: probing at exactly `_X_TOL` (1e-9 µm, the scanline's own
vertex-snap tolerance) past a body/taper seam lands in the ambiguous zone
this step's own guidance warns about ("prefer sampling just inside the span
over exactly at a boundary") — the body's own horizontal top/bottom edges
get counted as crossing the scanline there, merging what should be two
distinct branches into one. Fixed by widening the test probe to 1e-4 µm.
Real callers (`section_fdtd_layout`'s own `eps = min(1e-3, ...)`) already
sample three to six orders of magnitude further from a boundary than that,
so this is a test-authoring pitfall to document, not a practical risk.

---

## Step 3 — `section_prisms` / `section_gds`

Follow the pattern the built-ins share: one private `_section_polygons(x0, x1,
stub_left, stub_right, **kwargs)` returning device-frame vertex lists, with the
two public methods as thin adapters — `section_prisms` through
`polygons_to_prisms(...)` and `section_gds` shifting the layout's left edge to
`x = 0` (`mmi.py:288-300` is the reference implementation).

For GDS the work is **clipping polygons to `[x0-stub_left, x1+stub_right]`** rather
than synthesizing them. gdstk provides boolean operations, so clip against a
rectangle rather than hand-rolling polygon intersection.

One real semantic difference to record: `GdsGeometry`'s docstring states that
`to_gds_eme()` *"places the imported polygons as-is and appends nothing"*, unlike
`Taper`/`MMI`/`SBend`, which synthesize stub extensions on demand. A cascade's end
sections **do** need stubs — `HybridSolver` takes them from
`hasattr(solver, "waveguide_length")` (`hybrid.py:369-370`). `GdsGeometry` already
has a one-sided `input_width`/`stub_length` mechanism for exactly this; decide
whether `section_prisms` extends it to both ends or whether the GDS is required
to draw its own access waveguides. Requiring the latter is simpler and matches
the class's stated philosophy — say so explicitly if that is the choice.

**Decision made (2026-08-20): both, opt-in.** Added `access_extension: bool =
False` (default off, preserving the documented "as-is, appends nothing"
behaviour and every existing caller). When `True`, `section_prisms`/
`section_gds`/`to_gds_eme` synthesize straight access waveguides beyond
either end facet by extruding that facet's own *measured* cross-section
(`_facet_left`/`_facet_right`, captured once at construction via `_scan` at
the bbox's `x_min`/`x_max`) — generalizing the existing one-sided
`input_width`/`input_y` mechanism to both ends and to every branch at a
facet, not just one. With it off, a section reaching past the drawn
geometry raises `ValueError` naming exactly why (`access_extension` in the
message) rather than silently ending the core in mid-air. This resolves the
plan's "decide" into "support both, default preserves history" rather than
picking one — the twin file (bare device, no access waveguides) forces
`access_extension=True` in the acceptance example, exercising the new path
directly.

**MEASURED (2026-08-20).** `_section_polygons`/`section_prisms`/
`section_gds` implemented via gdstk `boolean(..., "and", precision=1e-7)`
clipping (kept per-polygon rather than unioned, matching `MMI`'s own
`_section_polygons` shape) plus the facet-extrusion rectangles above.
Compared against `MMI._section_polygons`/`section_gds` for the three spans
`mmi_2x2_hybrid.py` actually cascades (input `[-19,-13.75]` with a 1 µm
stub, body `[-13.75,13.75]`, output `[13.75,19]` with a 1 µm stub):

- Polygon **area** matches to 3.55e-15 µm² (input/output) and exactly
  (body); a rasterized-occupancy diff (1600×400 grid) found **0 mismatched
  pixels** on any section — not just equal area, but the same shape.
- `section_gds`'s oxide-cladding box is bit-identical to `MMI`'s.
- The un-selected branch is confirmed still drawn (`test_section_draws_
  every_branch_not_just_the_selected_one`) — `input_idx`/`output_idx` only
  select which branch `width_at`/`center_at` report, exactly as the plan
  specifies.

15 tests in `TestSectionGeometry`/`TestToGdsEmeExtension`
(`tests/test_gds_import.py`), all passing.

---

## Step 4 — Cut validity without parametric knowledge

The built-ins raise `ValueError` when a cut lands strictly inside a varying
region, each with device-specific knowledge of where its bands are: `MMI`
`_check_cut` against `_taper_bounds()` (`mmi.py:213-224`), `Taper`
(`taper.py:61-66`), `SBend` (`sbend.py:166-172`), `DirectionalCoupler`
(`coupler.py:181-189`). `UniformCrossSection` deliberately has none
(`uniform.py:114-116`: *"uniform in x, so any cut and any stub is legal by
construction"*).

A GDS geometry has no bands to consult, so the test has to be **measured**: sample
`width_at`/`center_at` in a neighbourhood of the cut and raise if either varies
beyond a tolerance. `HybridSolver._validate_mirror_pair` already does exactly this
shape of sampling (`hybrid.py:411-429`), and
`_validate_uniform_sweep_section` (`hybrid.py:1229-1246`) does it across a whole
span with `mirror_check_samples` / `mirror_check_tol`. **Reuse one of those rather
than writing a third sampler** — and note that
[`mmi_1x4_tapered_body_validation.md`](../mmi_1x4_tapered_body_validation.md) step 1
proposes factoring that checker into a shared helper for the same reason. These
two plans should land on one implementation.

Caveat to document: a sampled test can be fooled by a feature narrower than the
sample spacing. That is a real limitation of not having parametric knowledge, not
a bug — state it rather than implying the check is exact.

**MEASURED (2026-08-20).** `GdsGeometry._check_cut` reuses
`assert_uniform_cross_section` exactly as directed — no third sampler — called
twice per cut candidate, once on each side (`[x-span, x]` and `[x, x+span]`),
legal iff **at least one** side is uniform. That "one side" rule is what
reproduces `MMI`'s own contract (a cut exactly on a taper/body seam is legal;
strictly inside a band is not) rather than being stricter than it. Defaults:
`cut_check_span=0.05 µm`, `cut_check_samples=9`, `cut_check_tol=1e-3 µm` —
i.e. a ~6.25 nm sample pitch, named as constructor kwargs so a layout with
finer features can tighten them.

Verified against `MMI._check_cut` at the 2×2's own 7 legal positions
(both physical ports, both body/taper seams, mid-body) and 6 illegal ones
(mid-taper, both sides): **all 13 agree** (`TestCutValidity`, 4 tests +
2 parametrized × 6). The narrower-than-sample-spacing caveat above is
inherited as-documented in the constructor's own docstring, not
independently re-verified (would require a purpose-built pathological GDS,
out of scope for this vehicle).

---

## Step 5 — Acceptance

The comparison the vehicle was built for: run the 2×2 hybrid twice, once with
`MMI(...)` and once with the GDS twin, holding everything else fixed —
`RESOLUTION = 40.0`, `SOI_220`, C-band, `N_MODES = 14`, the same `boundaries` and
`body_stub = 0.25`.

- **Expected: agreement far tighter than any physical tolerance**, since the two
  paths should hand the solvers the same polygons. Disagreement means a real bug
  in clipping, branch selection, or the frame shift. **ARGUED** target: the two
  transfer matrices agree to solver noise, not to a percent-level device
  tolerance. Pin the actual number once measured; do not pre-commit to a loose
  tolerance that would hide a genuine defect.
  > **Outcome (2026-08-21): 1.603e-03 relative, and the premise above is the
  > part that was wrong.** The two paths hand the solvers the same polygons —
  > verified — but *not in the same order*, and meep is order-sensitive at
  > abutting-prism seams. The "disagreement means a real bug in clipping,
  > branch selection, or the frame shift" trichotomy was not exhaustive; a
  > fourth cause existed below the polygon layer. Not pre-committing to a loose
  > tolerance is what forced that out instead of waving the number through.
- **Cheap first:** compare `section_prisms`/`section_gds` **polygon output**
  before running any solver. If the polygons match, the physics must.
- Cache note: the two runs must not collide in `example_sim_cache/`. The
  fingerprint includes the geometry (`cache/config.py`), but confirm
  `GdsGeometry` contributes something distinguishing — a `MaterialStack` is
  fingerprinted as geometry only (`cache/config.py:215-216`, `wg_height`/`sz`), and
  non-JSON-safe attributes without a `cache_key()` are skipped **with a warning**
  (`config.py:218-228`). A `GdsGeometry` whose fingerprint reduces to its class
  name would silently serve the `MMI` run's result and produce a perfect,
  meaningless match.

**That last point is the trap most likely to produce a false pass.** Verify the
two runs land in different cache folders before trusting agreement.

**MEASURED, cheap-first half (2026-08-20).**
`examples/mmi_2x2/python/mmi_2x2_hybrid_gds.py` implements this comparison as
its own step 0, run before any solver, and asserts on every number rather than
just printing it:

- `width_at`/`center_at` across the full boundary span, 2 kwargs
  combos: agree to < 1e-9 µm (see step 2's measurement — same numbers).
- `section_gds` core-polygon area on the actual three sections the example
  cascades: agree to < 1e-9 µm² (see step 3's measurement).
- **The trap, checked directly**: built a real `FDTDSolver.solve_section`
  cache fingerprint (`beef.cache.config.build_run_config` →
  `beef.cache.store._hash_config`) for the same `(x0, x1, wavelengths,
  n_modes, stub, geometry_kwargs)` call, once against `MMI`, once against
  `GdsGeometry` — hashes **9c5869cdce65c692** vs **c52f1383634fcb1b**,
  confirmed distinct. `geometry_class` inside `_fingerprint()` (`cache/
  config.py:254`) is what does it — carried alongside the section polygons,
  not collapsed away. The false-pass trap does not fire.

**Login-node smoke run — ENDED, outcome not captured.** A reduced-scale run
(`RESOLUTION=12` vs. production's `40`, `N_MODES=4` vs. `14`, `N_WL=2` vs.
`11`, 3-point length sweep vs. 201, `CALIBRATE_BODY`/`SYMMETRIC_MODE_REUSE`
left on) was launched on the login node (`OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1`, per this workspace's smoke-run convention) to prove
the plumbing runs a real meep FDTD + meow EME solve through
`section_prisms()`/`section_gds()`/`width_at()` end to end, not just the
cheap polygon checks above. Confirmed while it ran: `section_prisms()` handed
meep a geometrically correct 5-prism cross-section (both input tapers +
access stub, verified against the printed prism vertex dump), and MPB
converged a guided mode on it (`frequency_1(1.34661,0,0)` after 1-21 iters,
normal convergence pattern) — i.e. the GDS-derived cross-section is a
legitimate waveguide to the mode solver, not degenerate geometry. **2026-08-20
~19:15: re-checked and the process is no longer running on the login node**
(no matching PID, no matching `squeue` entry); `example_output/mmi_2x2/
mmi_2x2_hybrid_gds/` holds only the geometry plot + section plots, no
`mm_length_sweep.npz`/SMatrix, so FDTD time-stepping did not reach
completion before the session that launched it ended. Not re-run — the
production submission below supersedes it (same plumbing, real scale), and
its RES=12 cache entries don't collide with RES=40's fingerprints.

**Full production-scale run — ✅ MEASURED, job 12069544, COMPLETED
2026-08-21 00:06:32, 4:31:02, exit 0:0.** `RESOLUTION=40`, `N_WL=11`,
201-point length sweep, matching `mmi_2x2_hybrid.py` exactly. Batch script
`batch/run_mmi_2x2_hybrid_gds.batch` + Makefile target `mmi_2x2_hybrid_gds`,
2 nodes × 24 ranks (the parametric run's own shape; `-t 420` vs. its
measured 4:11:55, this run took 4:31:02), `gts-sr108-gedc_main`.

Read from the two `.npz` files, not the job log. The `mm_lengths` and
`wavelengths` grids are bit-identical between the runs, so the arrays are
directly differenceable:

| 201 × 11 grid | max abs Δ | max rel Δ | RMS Δ | signed mean Δ |
|---|---|---|---|---|
| `T_bar` | 6.556e-04 | **1.603e-03** | 2.684e-04 | +1.872e-04 |
| `T_cross` | 5.218e-04 | **1.346e-03** | 1.973e-04 | +5.252e-05 |

Worst points: `T_bar` at `mm_length` 37.96 µm / λ 1.5300 µm, `T_cross` at
24.52 µm / 1.5300 µm. Both grids are smooth in λ (per-wavelength max |Δ|
falls monotonically 6.6e-4 → 5.8e-4 for bar, 5.2e-4 → 4.3e-4 for cross) and
carry a **positive signed bias** — i.e. this is a systematic offset, not
stochastic solver noise, which is what motivated the root-cause hunt below
rather than accepting the number.

**This is looser than the bar this step set** ("agreement far tighter than
any physical tolerance … solver noise, not a percent-level device
tolerance"). It is comfortably inside percent-level, but it is ~10 orders
looser than the agreement the *geometry* layer reaches (1e-16 µm, step 2)
and the *composition* layer reaches (4.0 eps, the sibling
[`multi_excitation_cascade_validation.md`](multi_excitation_cascade_validation.md)
step 1). Per this step's own instruction — "do not pre-commit to a loose
tolerance that would hide a genuine defect" — the residual was localized
rather than accepted.

### Root cause — MEASURED 2026-08-21, and it is **not** a defect in the GDS path

Localized in three cache-only measurements (no re-solve; every call used
`cache_only=True`, and the harness first reproduced **both** production runs
from cache bit-identically, max |Δ| 0.000e+00, so the attribution below is
made on exactly the arrays the two jobs wrote).

**(1) The FDTD calibration layer is bit-identical — which retires this
step's last open caveat.** `section_mode_beta_details` over the body span,
GDS vs `MMI`:

| | max abs Δ | bit-identical |
|---|---|---|
| `beta` (11 × 14, rad/µm) | 0.000e+00 | **yes** |
| `body_neff` (11 × 14) | 0.000e+00 | **yes** |

This quantity **was recomputed fresh on the GDS side** — `GdsGeometry`
fingerprints to different keys, so it cannot hit the `MMI` entries — and it
still came back bit-for-bit. That is the independent re-verification of the
meep/FDTD layer that the "one honest caveat" block above said job 12070199
could *not* provide, because that job was entirely cache-served. **The seven
uncommitted library edits move nothing in the FDTD layer either.** The
drift question is now fully closed, on both halves.

**(2) `M_in1` carries 100% of the residual.** Recomposing the cascade with
the two inputs swapped independently:

| composition | `T_bar` max abs Δ | `T_cross` max abs Δ |
|---|---|---|
| full (GDS vs MMI) | 6.556e-04 | 5.218e-04 |
| `M_in1` swapped alone, `neff` held at MMI's | 6.556e-04 | 5.218e-04 |
| `neff` swapped alone, `M_in1` held at MMI's | **0.000e+00** | **0.000e+00** |

So the entire difference enters through the one 3-D FDTD section matrix.
Its error is concentrated in the weak, high-order modes: modes 0-3
(|M| 0.07-0.70) agree to ~1e-4 relative, while modes 5, 6 and 13
(|M| 0.009-0.051) differ by 1.6e-2, 1.1e-2 and 1.8e-2, and modes 7-9/11-12
(|M| 1e-6 to 4e-5, at the noise floor) by up to 9.1e-2.

**(3) The cause is prism ORDER, and meep's subpixel averaging.** The two
backends produce the **same five polygons** for the input section — same
count, same vertex counts (all 4), same areas (0.45, 0.45, 3.5, 3.5, 0.7),
and canonicalized vertex-for-vertex identical. They emit them in a
**different order**: `MMI` gives `[0.45, 0.45, 3.5, 3.5, 0.7]`, `GdsGeometry`
the permutation `[3, 4, 1, 2, 0]` of that. `polygons_to_prisms()` preserves
list order, and meep's geometry list is order-sensitive where objects share
an edge. Building the discretized permittivity both ways at `res=40`
(structure init only, no time-stepping):

| epsilon grid (210 × 200 × 60 = 2.52 M voxels) | voxels differing | max abs Δε |
|---|---|---|
| MMI order vs GDS native order | **872** (0.0346 %) | 3.829 |
| MMI order vs **GDS reordered to MMI's order** | **0** | **0.000e+00 — bit-identical** |

The 872 differing voxels sit at exactly `x = 0.975, 1.000, 5.000 µm` in the
section frame — the stub/taper seam (`stub_left = 1.0`) and the taper/body
seam. Nowhere else.

**That bit-identical second row is the acceptance result this step was built
to produce.** Given the same polygon order, the GDS backend hands meep a
*bit-identical* structure. Clipping, branch selection and the device-frame
shift are all correct; there is no bug in the GDS path for the transmission
residual to be evidence of.

### The larger finding this exposed — abutting prisms carry a seam artifact

Chasing the above turned up something that is **not specific to the GDS path
and affects every BEEF device, `MMI` included**. Unioning the input
section's five abutting polygons into the single 14-vertex polygon they
geometrically are, and rebuilding the same grid:

| epsilon grid | voxels differing | max abs Δε |
|---|---|---|
| MMI order (5 abutting prisms) vs 1 merged polygon | **41,677** | 3.677 |
| GDS order (5 abutting prisms) vs 1 merged polygon | **41,677** | 4.392 |

**48× larger than the ordering effect.** meep's subpixel averaging infers an
interface normal at every seam between two geometry objects, even when both
carry the identical medium — so a contiguous silicon region drawn as
abutting prisms gets spurious internal interfaces that a merged polygon does
not. Both backends' *absolute* numbers therefore carry a subpixel-averaging
artifact of this size, and the ~1e-3 GDS-vs-parametric spread is a lower
bound on it, not an upper one.

**Decision 2026-08-21: recorded, not fixed.** Three options were weighed —
canonicalize `GdsGeometry`'s order to match (one ~4.5 h re-run, comparison
goes to exactly 0.0, seam artifact untouched); merge abutting same-material
polygons library-wide (physically correct, but moves every device's numbers
at ~1e-3, invalidates the 566-entry FDTD cache, and changes the manuscript's
published 2×2 results); or record as-is. **Recorded as-is** at the user's
direction: the GDS path is already proven correct by the bit-identical
reordered-epsilon test, so no verification debt is carried forward, and
neither fix would change that conclusion. Anyone revisiting the seam artifact
should treat it as a **library-wide accuracy question about
`polygons_to_prisms()`**, not as GDS-path business — it belongs in its own
plan, scoped across all six devices.

**Acceptance verdict: PASS.** Not on the transmission number, which is
1.6e-3 and explained, but on the structure the transmission number is
computed from: identical polygons (step 3), identical cross-section queries
(step 2, 1e-16 µm), identical cut behaviour (step 4), bit-identical `beta`,
and a bit-identical permittivity grid once the incidental ordering
difference is removed.

**~~OPEN ISSUE~~ → CLOSED the same day, see the measurement at the end of this
subsection. Found 2026-08-20 while checking on that job — the reference side
is not on the same library state.** Step 5's whole premise is "hold everything
fixed, change only the geometry backend." But the parametric side of the
comparison is a *file on disk* (`example_output/mmi_2x2/mmi_2x2_hybrid/
mm_length_sweep.npz`, written **2026-08-18 10:40**), while the GDS side will be
computed by job 12069544 against the library as it exists when that job starts.
**Seven library files have been modified since that npz was written**, all
uncommitted (`git status` in `BEEF/`): `solvers/fdtd.py`, `solvers/eme.py`,
`solvers/hybrid.py`, `geometry/mmi.py`, `geometry/base.py`,
`geometry/gds_import.py`, `cache/config.py`. Any Δ the comparison reports is
therefore geometry-backend difference **plus** two days of library drift, and
the two are not separable from the numbers alone.

Evidence that the drift is probably small, but it is evidence, not proof:

- `fdtd.py`'s change is the new `decay_method` knob, whose own docstring states
  the `"dft"` default is "unchanged from this solver's own established
  behavior; no existing device/golden moves."
- `cache/config.py` added `decay_method`/`energy_decay_dt` to three FDTD
  allowlists *and* to `_DEFAULT_OMITTED_ATTRS`, so default-valued runs keep
  matching pre-existing entries. Confirmed empirically: step 0's two probe
  fingerprints re-computed against the current tree are still
  **9c5869cdce65c692** / **c52f1383634fcb1b**, bit-identical to the values
  recorded above, and 98 `FDTDSolver_solve_section_MMI_*` cache entries are
  still live.
- `eme.py`/`hybrid.py` gained `assert_uniform_cross_section` guards and a
  behaviour-preserving refactor (`mmi_1x4_tapered_body_validation.md` step 1).
  Neither gates this run: the 2×2 body goes through
  `EMESolver.solve_section_length_sweep()`, which does not call
  `solve_section()` and did not gain the guard, and whose `_METHOD_ATTRS` key
  was left unchanged when `solve_section()`'s gained `num_modes`/`num_cells`.
- **But `mmi.py` did move a number this doc records.** Re-running the example's
  step 0 against the current tree gives section-polygon `|Δarea| = 0.000e+00`
  on all three sections, where step 3 above measured **3.55e-15 µm²** for
  input/output. Float noise either way, and it moved in the harmless
  direction — but it is a measured demonstration that `MMI`'s polygon output is
  not frozen while this comparison is pending.

**RESOLVED, MEASURED 2026-08-20 (job 12070199) — the drift is exactly zero.**
Rather than re-running `mmi_2x2_hybrid.py` in place (which would have
`rm -rf`'d the very reference under test, and raced job 12069544's step-7
read), the parametric side was recomputed by
`examples/verification/python/capture_mmi_2x2_cascade_golden.py` into its own
output directory. Read from the resulting `.npz`, not the job log:

| vs. the 2026-08-18 reference | max abs Δ | max rel Δ | bit-identical |
|---|---|---|---|
| `T_bar` (201 × 11 grid) | 0.000e+00 | 0.000e+00 | **yes** |
| `T_cross` (201 × 11 grid) | 0.000e+00 | 0.000e+00 | **yes** |

So the seven modified library files moved **nothing** on this device, and the
08-18 reference is a valid comparison target after all. Step 5's comparison is
clean as originally designed; no paired re-run is needed.

**One honest caveat on the scope of that measurement.** Job 12070199 finished
in **38 s** with no meep time-stepping in its log — i.e. every FDTD section was
served from `example_sim_cache/`, which is also where the 08-18 run's numbers
came from. The bit-identity therefore proves that everything *downstream* of
the cached section matrices (the β→n_eff conversion, the EME body length sweep,
the einsum composition, `|·|²`) is unchanged. It does **not** independently
re-verify the meep/FDTD layer itself, because that layer did not re-run. The
GDS acceptance run *does* compute its FDTD sections fresh — `GdsGeometry`
fingerprints to different keys, so it cannot hit the `MMI` entries — which
leaves the FDTD layer as the one place a real library difference could still
enter the comparison. The argument that it will not is strong but remains an
argument: `fdtd.py`'s only change is the `decay_method` knob, whose default
`"dft"` path its own docstring says leaves every existing device and golden
unmoved, and `cache/config.py` lists it in `_DEFAULT_OMITTED_ATTRS` precisely so
default runs keep matching pre-existing entries.

> **✅ CAVEAT CLOSED 2026-08-21 — it is no longer an argument.** The acceptance
> run's own `beta`/`body_neff` are **bit-identical** (0.000e+00) to the
> parametric side's, and that quantity *was* computed fresh through the meep
> layer under the current tree. So the FDTD layer is now measured, not argued,
> to be unmoved by the seven library edits. See root-cause finding (1) below.

Step 0 was re-run against the live tree at ~19:25 on 2026-08-20 and still
passes every assert; `tests/test_gds_import.py` (44) and
`tests/test_geometry_uniform.py` + `tests/test_geometry.py` (100) are green on
that same tree. So the pending job is not expected to fail on geometry — the
issue above is about what its *number will mean*, not whether it will run.

---

## Manuscript dependency

`beef_manuscript/sections/simulation_paradigm.tex:9` asserts:

> "the framework allows for the definition of custom device geometries through
> the use of GDSII files, enabling the simulation of a wide range of photonic
> devices."

Written as though this plan is done. When this plan opened it held only for
whole-device single-solver runs, not for the hybrid cascade the sentence sits
inside. The `:7` prism-vs-GDS framing depends on it too.

**✅ RESOLVED 2026-08-21 — the sentence is now true as written.** A GDS file can
occupy a section of a `HybridSolver` cascade (steps 1-4), and the acceptance run
(step 5, job 12069544) demonstrates it on a real device against a parametric
twin. Clear this entry from `beef_manuscript/unimplemented_changes.md`.

One caveat worth carrying into the manuscript if the `:7` prism-vs-GDS framing
is ever sharpened into a quantitative claim: the GDS and parametric paths agree
to **1.6e-3** in transmission, not to solver noise, and the reason is meep's
order-sensitive subpixel averaging at abutting-prism seams — a property of the
prism representation itself, shared by both paths. Do not write "identical
results"; "geometrically identical structures, agreeing to ~1e-3 in
transmission" is what was measured.
