# GDS-file geometry in a hybrid cascade — build plan

Opened **2026-08-20**. `GdsGeometry` can back a whole-device single-solver run,
but it **cannot occupy a section of a `HybridSolver` cascade**. This doc orders
the work to close that, validated against a GDS file cut to the 2×2 MMI's exact
dimensions so the result can be compared directly with the parametric `MMI` run.

> **This doc is a PLAN, not a gap list** — same category as
> [`mmi_1x4_sin300_validation.md`](mmi_1x4_sin300_validation.md). It is a
> *library capability* plan rather than a device plan, so it converts to a
> regression suite rather than to a four-checks list.
>
> **Nothing below has been simulated.** The library-behaviour findings were read
> directly off the current source with line references. The acceptance tolerance
> in step 4 is **ARGUED** until measured.

---

## Summary

| | status |
|---|---|
| Read a GDS file into BEEF's polygon convention | ✅ `load_gds_polygons()` (`geometry/gds_import.py:19`), via gdstk, µm coordinates |
| GDS geometry in a whole-device **EME** run | ✅ `GdsGeometry.to_gds_eme()` |
| GDS geometry in a whole-device **FDTD** run | ✅ `polygons()` → `to_meep_prisms()`; used by the HCG reflector example |
| GDS geometry in a **hybrid cascade** | ❌ **not supported** — `section_prisms`, `section_gds`, `width_at` all unimplemented |
| Detecting a varying cross-section without parametric knowledge | ❌ does not exist |

Ordered work: **validation vehicle (step 1) → `width_at`/`center_at` (step 2) →
`section_prisms`/`section_gds` (step 3) → cut validity (step 4) → acceptance
(step 5).**

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
[`mmi_1x4_tapered_body_validation.md`](mmi_1x4_tapered_body_validation.md) step 1
proposes factoring that checker into a shared helper for the same reason. These
two plans should land on one implementation.

Caveat to document: a sampled test can be fooled by a feature narrower than the
sample spacing. That is a real limitation of not having parametric knowledge, not
a bug — state it rather than implying the check is exact.

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

---

## Manuscript dependency

`beef_manuscript/sections/simulation_paradigm.tex:9` asserts:

> "the framework allows for the definition of custom device geometries through
> the use of GDSII files, enabling the simulation of a wide range of photonic
> devices."

Written as though this plan is done. Today it holds only for whole-device
single-solver runs, not for the hybrid cascade the sentence sits inside. The
`:7` prism-vs-GDS framing depends on it too. Tracked in
`beef_manuscript/unimplemented_changes.md`.
