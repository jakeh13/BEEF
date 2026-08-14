# High-contrast grating (HCG) reflector — FDTD-only build plan

Opened **2026-08-04**, at the user's request: an FDTD-only example (no EME, no
`HybridSolver`) of the high-contrast-grating reflector in
`examples/gds_files/HCG_seeded_TO.gds`, loaded via `beef.GdsGeometry`, on
160 nm SOI, with a script and a notebook version.

> **This doc mixes the plan and gap-list conventions** (see [README.md](README.md)).
> **Updated 2026-08-04 (evening) — the RE-RUN production run has landed, and
> the honest number is worse, not better.** A third real issue was found,
> this time by the user directly, after the first two mistakes were already
> fixed and a paper-matching production run had already landed: the input
> waveguide was ending abruptly instead of continuing through the PML — a
> real physical defect (see ⚠ item 3 below), not a cosmetic one, that
> changes the actual simulated geometry. Every downstream number needed
> re-measuring. **Corrected: peak R = 97.14% at 1.580 µm, against the
> paper's 98.9% at 1.569 nm — a 1.76-point gap, LARGER than the (flawed)
> pre-fix run's 0.8-point gap.** The old, better-looking number was
> partly an artifact: removing the facet reflection also moved the
> resonance peak by ~30 nm (from 19 nm blue of the paper's to 11 nm red of
> it), which the pre-fix magnitude comparison alone did not reveal. This is
> reported as the honest current result, not smoothed over — see step 9 for
> the full spectrum and step 10 for what it opens (chiefly: does resolution
> above `res=60` keep moving the resonance position).
>
> **Updated 2026-08-05 — partially checked, still open.** A single-wavelength
> `res=70`/`80` probe (job 11688701) shows R still climbing by a non-trivial
> amount through res=80, with the deceleration that looked clean through
> res=60 stalling rather than continuing — evidence the resolution residual
> has NOT settled by res=60, making the peak-position question MORE open,
> not resolved. See step 8's current-numbers table and step 10 item 1.
>
> ⚠ **Three real mistakes were caught and fixed while building this
> example — all three are load-bearing for reading anything below:**
> 1. **Wrong band, most of a day's work.** Everything through step 7 (and
>    this session's first convergence sweep, since deleted) used **O-band**
>    (1.26–1.36 µm), matching every other `SOI_160` example in this repo. The
>    user asked to check the paper: it reports "peak reflectance of 98.9%
>    ... centered at 1569 nm ... minimum reflectance of 97.9% over **C-band**"
>    — this device is a resonant mirror seeded for a specific band, not a
>    broadband splitter. Every O-band number this doc, the script, or the
>    notebook ever produced is **superseded, not merely mislabelled** — R at
>    the wrong wavelength for a resonant structure is not informative about
>    the real device at all. Fixed: `WAVELENGTH`/`WL_MIN`/`WL_MAX` moved to
>    1.569 µm / 1.50–1.60 µm (C-band, matching this repo's own
>    `directional_coupler`/`mmi_2x2` convention) everywhere. See step 8.
> 2. **The x-domain (PML clearance) sweep was self-corrupting.**
>    `source_margin`/`monitor_margin` are measured from the **PML inner
>    edge**, not the geometry's own edge (`solve_reflection()`'s documented
>    contract). Sweeping `edge_clearance` without adjusting them moves the
>    source relative to the actual silicon — one task's source ended up
>    fully outside the 1.082 µm input stub (R=0.7530, meaningless), the next
>    had **both** source and monitor outside it (R=0.0046 — a source
>    launched into bare oxide with no mode to couple into, essentially
>    zero). Both looked superficially plausible as "R" numbers; neither
>    meant anything. Fixed: `source_margin`/`monitor_margin` are now derived
>    from `edge_clearance` so they sit at a fixed physical offset from the
>    geometry's own edge regardless of PML placement — verified a no-op at
>    the baseline (cache hit, identical R) and confirmed correct once fixed
>    (`edge_clearance=0.5` reproduces the baseline `R=0.8962` exactly). See
>    step 8.
> 3. **The input waveguide ended abruptly instead of continuing through the
>    PML — a real physical defect, not a cosmetic one.** Before this fix,
>    `GdsGeometry.polygons()` never extended the input waveguide (A2's own
>    "appends nothing" contract), so it simply stopped at the GDS's own
>    edge, sitting in the `edge_clearance` gap (or right at the PML if
>    `edge_clearance=0`) as a clean-cut Si/oxide facet. **That facet
>    reflects part of the guided mode on its own, independent of the
>    grating** — a real backward-traveling wave the monitor cannot
>    distinguish from an actual reflection off the structure under test.
>    Caught by the user directly ("extend the input waveguide into the pml
>    so that it is not abruptly cut off"). Fixed: `GdsGeometry` gained
>    `input_width`/`input_y` (the GDS's own known input-stub cross-section,
>    0.5 µm) and now honors `polygons(stub_length>0)` by extending that ONE
>    waveguide (input side only — the far side is an arbitrary freeform cut
>    with no defined cross-section to extend) by `stub_length`;
>    `solve_reflection()` gained `extend_input_waveguide=True`, which
>    computes `stub_length = edge_clearance + dpml` — enough to run the
>    waveguide fully through the PML to the true cell boundary — and feeds
>    it through to both the actual meep geometry and the cache fingerprint.
>    Verified visually (a rendered ε-map shows the waveguide now passing
>    continuously through the PML with no cut) before spending any real
>    compute on it. **Removing this artifact changed the actual simulated
>    geometry, so every cached result was invalidated** — the full 17-point
>    convergence sweep and the production run both needed re-running, not
>    re-labelling. See step 8/9.
>
> All three were caught by the user asking directly or requesting a specific
> check (band; x/y/z convergence rather than y-only; the waveguide
> extension) rather than by anything in this pipeline catching itself on its
> own. **Mistakes 2 and 3 (source/monitor outside the real silicon; a
> geometry that claims `extend_input_waveguide` support but doesn't actually
> reach the domain boundary) now DO have an automated guard** —
> `check_reflection_layout()` in `solve_reflection()`, added 2026-08-04,
> raises before `mp.Simulation` is even built. See step 10 item 6. Mistake 1
> (wrong band) has no equivalent guard and likely can't — nothing in the
> codebase encodes what band a given device was designed for to check a
> requested wavelength range against.
>
> **Placement.** Moved into `deferred/` 2026-08-13, at the user's explicit
> instruction, after confirming nothing had changed since the 2026-08-05
> update above. **This is a placement decision, not a claim that step 10 is
> closed** — the paper-comparison gap (97.14% vs. 98.9%, step 9) and the
> open resolution-convergence question (does R keep moving past res=60,
> step 8/10 item 1) are both exactly as open as they were on 2026-08-05.
> `deferred/` fits better than the top level or `resolved/` because what's
> actually blocking further progress is a real resource decision (a
> res=70/80 full-band production sweep, estimated multi-hour on the
> single-wavelength res=80 cost of 262 min alone) that has not been made,
> not an unresolved technical question this pipeline could close on its
> own — i.e. work intentionally postponed, matching README.md's own
> definition of `deferred/`.

## Reference — what this device is

`HCG_seeded_TO.gds`'s cell name and filename ("seeded TO") match the seeded
topology-optimization methodology and the HCG reflector device published in:

> J. M. Hiesener, R. P. Pesch, and S. E. Ralph, "On Topology Optimization
> Strategies for Ultra-Compact High Contrast Grating Design," *IEEE Photonics
> Technology Letters*, vol. 36, no. 9, pp. 597–600, 2024,
> doi: [10.1109/LPT.2024.3372429](https://doi.org/10.1109/LPT.2024.3372429).

Cited in both the script (`examples/high_contrast_grating/python/high_contrast_grating_fdtd.py`)
and the notebook's header cell. Comparing a converged run of this example
against that paper's own reported spectrum is the natural correctness
check — see step 9 for current status (a production run is in flight as
of the ⚠ box's latest update, following the waveguide-extension fix).

## Summary — where the work actually was

| | status |
|---|---|
| GDS geometry itself | ✅ inspected directly (gdspy) — one 500 nm access waveguide + 18 freeform arc polygons, layer (2,0) only |
| `gdstk` importable in the `jmh` env | ❌ → ✅ **FIXED 2026-08-04** — `pip install gdstk` (network confirmed reachable; user approved installing libraries this session) |
| `GdsGeometry` usable by any FDTD path | ❌ → ✅ **FIXED 2026-08-04** — `GdsGeometry.polygons()` added (`beef/geometry/gds_import.py`) |
| A physically meaningful FDTD measurement for a one-port device | ❌ → ✅ **ADDED 2026-08-04** — `FDTDSolver.solve_reflection()` (`beef/solvers/fdtd.py`) |
| Domain/source/monitor placement | ✅ **MEASURED** via a rendered ε-map, not just derived (`geometry.png` in the example's own output) |
| Script + notebook | ✅ **DONE**, both run to completion, byte-identical results |
| Full fast test suite after the two library changes | ✅ **611 passed, 19 skipped** (2026-08-04) — no regression |
| Correct band (C-band, not O-band) | ❌ → ✅ **FIXED 2026-08-04** — see the ⚠ box above; a real mistake, not a preemptive check |
| x-domain sweep self-corrupting (source outside waveguide) | ❌ → ✅ **CAUGHT AND FIXED 2026-08-04** — see the ⚠ box above — step 8 |
| Input waveguide abruptly cut off instead of extended through PML | ❌ → ✅ **CAUGHT AND FIXED 2026-08-04** — see the ⚠ box above (item 3) — step 8/9 |
| 3-D resolution convergence (single wavelength, C-band peak) | ✅ **MEASURED 2026-08-04, RE-MEASURED same day after the waveguide fix** — res 10→60, `R` 0.4468→0.9697 (post-fix; systematically lower than the pre-fix 0.5144→0.9789), same geometric-decay shape — step 8 |
| PML/domain convergence, x/y/z independently | ✅ **RE-MEASURED 2026-08-04 post-fix** — x fully converged (0.0000 spread), y converged at res=25 (±0.0001), z converged at res=25 (±0.0002) — even tighter than pre-fix — step 8 |
| Full-band (21 λ, C-band) 3-D production run | ✅ **RE-MEASURED 2026-08-04** — job 11673560 COMPLETED (54 min), peak R=0.9714 at 1.580 µm (job 11670282's 0.9810/1.550µm superseded) — step 9 |
| Comparison against the published paper's own spectrum | ✅ **RE-MEASURED 2026-08-04** — gap WIDENED to ~1.8 pt (was ~0.8 pt pre-fix) and peak flipped from blue- to red-shifted — reported honestly, not smoothed over — step 9/10 |

So the ordered work was: **geometry inspection → two library gaps → why
`solve()` cannot work here → the reflection method → verified placement →
script/notebook → test suite → smoke test → [band correction] → real 3-D
x/y/z/resolution convergence → [x-sweep bug + fix] → production run →
comparison against the paper → [waveguide-extension bug + fix] →
convergence RE-measured → production RE-measured → comparison
RE-measured, honestly worse.** What remains is step 10's genuinely open
items — chiefly whether resolution above `res=60` keeps moving the
resonance position (the ~30 nm shift from the waveguide fix alone is the
strongest signal yet that it might), plus whether this GDS is the paper's
final design or a staged seed, where the non-reflected power goes, and
generalizing the domain plot.

---

## Part A — findings

### A1. The GDS is a one-port device, not a two-port through structure — MEASURED

Read directly with `gdspy` (the `jmh` env's other installed GDS library,
used only for this one-time inspection; the shipped code path uses `gdstk`,
see A2):

```
layer (2, 0): 19 polygons, bbox x[-3.2500, 2.1390] (span 5.389 um)
                            y[-2.9910,  2.9920] (span 5.983 um)
```

18 of the 19 polygons have 100+ vertices each (freeform, topology-optimized
boundaries); **one** — `x[-3.250, -2.168], y[-0.250, 0.250]` — is a plain
4-vertex rectangle: a 500 nm-wide, 1.082 µm-long straight access waveguide on
the far left. Rendering all 19 (`/tmp/.../scratchpad/hcg_gds_view.png` this
session, reproduced by the example's own `geometry.png`) shows the freeform
region is a set of **concentric arcs** fed by that one waveguide — a curved
grating mirror, consistent with A1's citation, not a periodic bar grating and
not a two-port MMI/coupler-style device. **There is no second waveguide on
the far side.** `solve()`'s entire API (`output_monitor_centers`, forward-only
mode coefficients) assumes one exists; it does not here.

### A2. Two real library gaps, both blocking, both fixed — MEASURED

**Gap 1 — `gdstk` not installed in `jmh`.** `beef.load_gds_polygons()`
lazy-imports `gdstk` inside the function body
(`beef/geometry/gds_import.py:42`). Confirmed broken before any other work
started:

```
>>> from beef import load_gds_polygons
>>> load_gds_polygons('examples/gds_files/HCG_seeded_TO.gds')
ModuleNotFoundError: No module named 'gdstk'
```

`gdspy`/`gdsfactory`/`klayout` were already present in `jmh`; `gdstk` was
not, despite being the one the shipped code actually calls. **Fixed by
`pip install gdstk`** (v1.0.1) rather than rewriting `gds_import.py` onto a
different backend — the user confirmed installing libraries is fine this
session, which is the lower-risk fix (zero behavior change to code that
other work may already depend on) versus swapping GDS parsers.
Re-confirmed working immediately after:
`load_gds_polygons(...)` returns `{(2, 0): [...19 polygons...]}`.

**Gap 2 — `GdsGeometry` had no `polygons()`.** `DeviceGeometry.polygons()`
(the method every FDTD path needs — `to_meep_prisms()` calls it directly,
`beef/geometry/base.py:76`) has a base-class default that raises
`NotImplementedError`. `GdsGeometry` (`beef/geometry/gds_import.py`)
implemented `to_gds_eme()` for EME only; nothing in the class ever
overrode `polygons()`. So before this session, **`GdsGeometry` could not be
used with `FDTDSolver` at all** — not "worked poorly", raised outright.

**Fix:** `GdsGeometry.polygons()` added, returning the imported polygons
as-is and — deliberately — **ignoring `stub_length`**, the same contract
`to_gds_eme()` already documented ("places the imported polygons as-is and
appends nothing"). Synthesizing an extra stub segment would need a local
waveguide width to extrude at the geometry's own edge, which a freeform
import cannot supply in general (no `width_at()` for an arbitrary cut
through an arbitrary polygon). This is a general fix to `GdsGeometry`, not
specific to this device — any future GDS-imported FDTD geometry benefits.

### A3. `FDTDSolver.solve()` cannot report anything physical for this device, even after A2 — MEASURED (derived, then confirmed by a rendered ε-map)

Two independent reasons, both read directly from `beef/solvers/fdtd.py`:

1. **`solve()` requires `output_monitor_centers`** (raises `ValueError`
   otherwise) and places those monitors at the *opposite* domain edge from
   the source. This device has no waveguide there to measure.
2. **`solve()` only ever reads the forward mode coefficient** —
   `coeff.alpha[0, :, 0]` — at its input normalization monitor
   (`fdtd.py:1000`, `solve()` body). It never reads `alpha[0, :, 1]`
   (backward), so even a device with nothing to measure on the far side
   could not have its *reflection* read out of `solve()` as shipped.

There is a third, subtler problem that would have silently broken any naive
attempt to route this device through `solve()`'s own domain-sizing, even if
1 and 2 were worked around by hand: `solve()` sizes the cell as
`Sx = geometry.total_length + 2*waveguide_length + 2*dpml`, which assumes
`geometry.polygons(stub_length=...)` extends the drawn silicon by
`waveguide_length` on each side to reach the PML — true for every
parametric geometry (MMI, Taper, DirectionalCoupler, ...), **false** for
`GdsGeometry` after A2's fix (which deliberately ignores `stub_length`).
Worked out algebraically: with `source_margin <= monitor_margin` (the
common case — e.g. `mmi_1x2_fdtd.py`'s own `source_margin=0.0 <
monitor_margin=0.25`), the source lands **at least 0.25 µm outside the
drawn silicon, in bare oxide, for any parameter choice** — not a tuning
problem, a structural one.

**This was not left as an argued derivation.** `solve_reflection()`
(A4) uses its own domain formula specifically to avoid this, and before
trusting it, the actual `mp.Simulation` epsilon map was rendered with the
chosen parameters overlaid (`/tmp/.../scratchpad/domain_check.png` this
session, reproduced as the example's own `geometry.png` — see A4/step 4):
source and monitor both land inside the real 500 nm waveguide, with ~0.58
µm of clearance before the arcs begin and 0.3 µm of plain-oxide clearance
before PML on both sides. The image is the evidence; the algebra explains
why it needed checking rather than being assumed correct by symmetry with
`solve()`.

### A4. Fix: `FDTDSolver.solve_reflection()` — new code, tested against real physics

Added to `beef/solvers/fdtd.py`, right after `solve()`. Mirrors `solve()`'s
structure (same `_media()`, `eps_averaging`/`subpixel_*`, `decay_tol`,
caching via `build_whole_device_config`/`lookup`/`save_run`) but:

- **No `output_monitor_centers`.** One `EigenModeSource` + one
  `add_mode_monitor`, both placed near the excited port only.
- **Reads both directions off that one monitor**: `alpha[0, :, 0]`
  (forward = incident, the normalization) and `alpha[0, :, 1]` (backward =
  reflected into the same mode) — meep's standard single-run S11
  decomposition. This is not a new convention invented for this device: it
  is the exact `alpha[...,0]`/`alpha[...,1]` forward/backward split
  `PortSpec.alpha_idx` already documents elsewhere in this file, itself
  validated there against `ring_resonator_hybrid_cascaded.py`'s
  `get_coeffs()`. No second "incident-only" reference run is needed.
- **Domain sizing does not reuse `_build_domain()`.** `Sx =
  geometry.total_length + 2*(edge_clearance + dpml)` — PML sits
  `edge_clearance` (plain oxide, default 0) outside the geometry's own
  bounding box, full stop. The caller picks `source_margin`/
  `monitor_margin` to land inside whatever real waveguide the geometry
  drew — there is no automatic guard, because a geometry cannot in general
  report its own stub length (A2). **Must be checked with a domain plot
  before trusting a run** — this is stated in the method's own docstring,
  not left implicit.
- **Returns a normal `SMatrix`** with a single entry,
  `(port="in", in_mode=input_mode, out_mode=input_mode)`, so
  `smat.transmission(port="in", ...)` reads R(λ) and
  `smat.check_passive()` works unmodified as a (necessary, not sufficient)
  sanity check.

---

## Part B — the ordered steps

### Step 1 — `gdstk` in `jmh` ✅ DONE 2026-08-04

Deliverable: `gdstk` importable. Evidence: `pip install gdstk` succeeded
(v1.0.1); `load_gds_polygons()` on the real GDS returns the expected
19-polygon, single-layer dict immediately after. Cost: seconds, no cluster
time.

### Step 2 — `GdsGeometry.polygons()` ✅ DONE 2026-08-04

Deliverable/evidence: see A2. `tests/test_gds_import.py` (pre-existing, 29
tests total across the three targeted files) still passes; full fast suite
611 passed / 19 skipped afterward (step 7).

### Step 3 — `FDTDSolver.solve_reflection()` ✅ DONE 2026-08-04

Deliverable/evidence: see A4. Exercised end-to-end by the smoke test (step
6), not just imported.

### Step 4 — domain/source/monitor placement, verified by a real ε-map ✅ DONE 2026-08-04

Chosen (2-D diagnostic, `res=15`, `dpml=0.5`): `edge_clearance=0.3`,
`source_margin=0.45`, `monitor_margin=0.8`, `sim_width=7.5`. Rendered
`mp.Simulation.get_epsilon()` with source/monitor/PML/bbox markers overlaid
— both land inside the real 500 nm waveguide, well before the arcs begin,
with PML pulled 0.3 µm off the structure's own edges on all sides. **This
is measured, not derived from the formulas alone** — A3 explains why that
distinction mattered here specifically.

⚠ These four numbers are **not** independently convergence-tested. In
particular the ~0.26 µm gap between the outermost arc tips (y≈±2.99) and
the y-PML inner edge (`sim_width/2 - dpml = 3.25`) has not been swept —
carried into step 8.

### Step 5 — script + notebook ✅ DONE 2026-08-04

- `examples/high_contrast_grating/python/high_contrast_grating_fdtd.py` —
  follows every other FDTD example's shape (`mmi_1x2_fdtd.py`'s section
  numbering), with an `HCG_TEST_MODE` knob mirroring
  `ring_resonator_hybrid.py`'s `RING_TEST_MODE` convention exactly (2-D,
  low res, few wavelengths, loose `decay_tol` under `TEST_MODE=1`; full 3-D
  production settings by default). `region_lengths` for `GdsGeometry` is
  **derived from the GDS file itself** (its own polygon bounding box), not
  hardcoded, so it can't silently drift from the file.
- `examples/high_contrast_grating/notebooks/high_contrast_grating_fdtd.ipynb`
  — one cell per script section, generated from the script (single source
  of logic, not hand-duplicated), matching `mmi_1x2_fdtd.ipynb`'s 9-cell
  shape.
- `plot_fdtd_geometry` (`beef.plot`) was **not** reused for the geometry
  plot — it calls `geometry.to_gds_eme(straight_length=...)` and assumes
  the stub-extension convention A3 explains does not hold for
  `GdsGeometry`, so it would draw source/monitor markers in the wrong
  place. The example ships its own small `plot_hcg_domain()` instead,
  built directly from `solve_reflection()`'s real placement formulas (this
  is exactly the `domain_check.png` diagnostic from step 4, promoted into
  the shipped example). **Generalizing this into `beef.plot` for any
  `solve_reflection()`-style device is a reasonable future step, not done**
  — carried into step 8's backlog.
- ⚠ **Finding, not a defect introduced here**: the notebook's `_HERE =
  os.path.dirname(os.path.abspath(__file__))` line (copied from every
  other notebook in this repo) raises `NameError: name '__file__' is not
  defined` under a bare `nbclient`/`ipykernel` execution with no notebook
  path injected. **Confirmed repo-wide, not specific to this file**: the
  pre-existing `mmi_1x2/notebooks/mmi_1x2_fdtd.ipynb`, copied and executed
  the identical way, fails identically. Left as-is (matches established
  convention) rather than "fixed" in only the new file, which would make it
  diverge from every other example notebook for no benefit — whatever
  environment actually runs these (a real Jupyter server, or `nbconvert`
  with a path-aware resource dict) supplies `__file__`, which is exactly
  what re-running with it injected confirmed (below).

### Step 6 — smoke test ✅ MEASURED 2026-08-04, both script and notebook (⚠ O-band — superseded, see below)

⚠ **This entire step was run at the wrong wavelength (O-band, 1.26–1.36 µm)
— see the top-of-doc box.** Kept verbatim below as the historical record
that the pipeline itself worked end-to-end before the band was corrected;
none of its R numbers describe the real device. Step 8 has the corrected,
trustworthy numbers.

Login node, `HCG_TEST_MODE=1`, `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`
(this workspace's standard smoke-test etiquette — see the
`duty-env-and-run-conventions` note). 2-D, `res=15`, `decay_tol=1e-6`, 5
wavelengths, ~32 s wall time (`run 0 finished at t=354.7, 10641 timesteps`).

```
Results at lambda = 1.310 um:
  Reflection R        : 0.0950
  Reflection loss     : -10.22 dB
  Peak R over band    : 0.1106 at 1.285 um
  Min R over band     : 0.0805 at 1.360 um
```

`smat.check_passive()` raised no warning (R stayed well under 1 at every
wavelength — a necessary, not sufficient, sanity check). The notebook was
independently executed end-to-end via `nbclient` with `__file__` injected
(see step 5's finding) and reproduced **the same numbers exactly**
(R=0.0950 at 1.310 µm, same peak/min) — real evidence the two deliverables
are consistent, not just visually similar.

**Read this number as "the pipeline works," not "the mirror reflects
9.5%."** It is a 2-D run (no vertical confinement — the repo's own
`FDTDSolver` docstring states 2-D is "quantitatively a different
structure") at `res=15` (no other device in this repo trusts anything
below `res=40`, and README finding 3 says convergence is per-device/
per-material, not transferable) with a loose `decay_tol=1e-6` against a
production default of `1e-11` on a structure whose concentric-arc geometry
plausibly resonates (multiple internal round trips before the field
decays) — exactly the situation where a loose tolerance under-counts real
reflected energy.

### Step 7 — regression check on the shared library changes ✅ DONE 2026-08-04

Per this workspace's verification standard (prefer a real re-run over an
argued "should be inert" claim): full fast suite,
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 pytest tests/ -q` →
**611 passed, 19 skipped** — identical pass count to before this session's
two library changes (`gds_import.py`'s new `polygons()`, `fdtd.py`'s new
`solve_reflection()`; nothing existing was modified in either file). Both
changes are additive (new method / new method), so this is the expected
result, confirmed rather than assumed.

### Step 8 — real 3-D resolution + x/y/z PML convergence ✅ MEASURED 2026-08-04, C-band, RE-MEASURED after the waveguide-extension fix

**Tool:** `examples/verification/python/hcg_reflection_convergence.py`
(temporary, delete once a production resolution/domain is settled — see its
own docstring), driven by `batch/verification/run_hcg_reflection_convergence.batch`
as a SLURM array (`gts-sr108-ece`, `-q inferno`, 1 node/24 ranks/task).
Every point is a genuine 3-D `solve_reflection()` call — confirmed via
`"Working in 3D dimensions."` in each point's own meep log, plus the real
160 nm core's z-face at ±0.08 µm — at a single wavelength, **1.569 µm**
(the paper's own reported peak, not a rounded guess), `decay_tol=1e-11`
(the production value, not loosened for the sweep — so every point's own
wall time *is* the decay_tol-behavior check, item 3 of the old plan).

#### Current numbers — WITH `extend_input_waveguide=True` (job 11672082, 17 tasks; res=70/80 added 2026-08-05 via job 11688701, all landed)

**Resolution series** (`sim_width=7.5`, `edge_clearance=0.3`, `sim_height=1.5`):

| res (pts/µm) | R | \|ΔR\| vs prev | wall (min, 24 ranks) |
|---|---|---|---|
| 10 | 0.4468 | — | 0.4 |
| 15 | 0.7030 | 0.2562 | 1.2 |
| 20 | 0.8164 | 0.1133 | 4.3 |
| 25 | 0.8757 | 0.0594 | 3.5 |
| 30 | 0.9081 | 0.0323 | 6.5 |
| 35 | 0.9323 | 0.0242 | 11.9 |
| 40 | 0.9450 | 0.0128 | 18.3 |
| 50 | 0.9592 | 0.0142 (over a +10 step — ~0.00142/unit) | 58.9 |
| 60 | 0.9697 | 0.0105 (over a +10 step — ~0.00105/unit) | 74.5 |
| 70 | 0.9741 | 0.0044 (over a +10 step — ~0.00044/unit) | 133.3 |
| 80 | 0.9780 | 0.0039 (over a +10 step — ~0.00039/unit) | 262.2 |

**Systematically lower than the pre-fix sweep at every single point** (see
"superseded" table below) — expected, since the pre-fix numbers included a
spurious facet-reflection contribution the fix removes.

**res=70/80, added 2026-08-05 to settle step 10's single most load-bearing
open question (does resolution above res=60 keep moving the resonance),
give an honest, not a reassuring, answer: R is still climbing by a
non-trivial amount, and the deceleration that looked clean through res=60
did NOT continue cleanly.** The res=60→70 step decelerated sharply (ratio
0.0044/0.0105 ≈ 0.42, faster deceleration than any earlier step) — read in
isolation that looked like the residual was collapsing. The res=70→80 step
then decelerated only mildly (ratio 0.0039/0.0044 ≈ 0.89) — almost flat,
not a continuation of the res=60→70 trend. **This non-monotonic
deceleration ratio means a single geometric extrapolation is not
trustworthy here** — naively continuing the res=70→80 ratio (~0.89) sums
to an asymptote *above* 1.0, which is unphysical and just shows the
extrapolation has broken down this close to the noise floor, not that R
genuinely exceeds 1. **Read honestly: R has not clearly plateaued by
res=80, wall time has grown from 74.5 min (res=60) to 262.2 min (res=80,
+252% for a +33% resolution increase — steep, consistent with 3-D FDTD's
usual worse-than-`res⁴` blowup near a resolution's practical ceiling), and
whether it continues climbing toward the 1.0 end of the previously-stated
R≈0.98–1.00 bracket or is genuinely near-converged in the high-0.97s to
0.98 range is still not settled.** Pushing to res=90/100 to settle it
further would very plausibly cost multiple additional hours per point on
this trend — a real resource decision, not attempted here without asking.

**Consequence for step 10's open question (does the resonance PEAK
*position* keep moving, not just R's magnitude at this one probe
wavelength):** this sweep only ever measures R at the single fixed
wavelength 1.569 µm — it cannot directly confirm or rule out a peak-position
shift on its own. But R at 1.569 µm moving by 0.0083 total from res=60 to
res=80 (comparable in size to the ~0.03 magnitude change the
waveguide-extension fix itself caused, which is what produced the
documented ~30 nm peak shift) is a real signal that the full-band
lineshape — and very plausibly the peak's own wavelength — has likely NOT
fully settled by res=60 either. **This makes the open question more
concerning, not less**: the existing production spectrum (step 9, peak
R=97.14% at 1.580 µm) was run at res=60, and this new data suggests that
result's peak position and magnitude may both still move at higher
resolution. Not re-run here (a full 21-wavelength production sweep at
res=70/80 would cost proportionally more than this single-wavelength
convergence point already did) — flagged as the natural next step if this
question needs closing, not closed by this update.

**`res=60`'s R=0.9697 sits further from the paper's 97.9% than the pre-fix
`res=60`'s R=0.9789 did** — worth reading plainly: removing a real
modeling artifact made the raw number *agree less well* with the paper at
this resolution, which means the resolution residual is doing more of the
work than the pre-fix comparison suggested, not that the fix was wrong
(the fix is independently verified correct by the ε-map render and by
x-domain sensitivity going to *exactly* zero, both reported below/above).

**Domain (PML) series, all at a single resolution:**

| series | swept | R values | spread |
|---|---|---|---|
| y @ res=10 | `sim_width` 7.5/8.5/9.5 | 0.4468 / 0.4425 / 0.4462 | ±0.004 |
| y @ res=25 | `sim_width` 7.5/8.5/9.5 | 0.8757 / 0.8760 / 0.8761 | ±0.0004 — converged |
| x @ res=25 | `edge_clearance` 0.3/0.5/0.8 | 0.8757 / 0.8757 / 0.8757 | **0.0000 — fully converged** |
| z @ res=25 | `sim_height` 1.5/2.0/2.5 | 0.8757 / 0.8758 / 0.8755 | ±0.0003 |

**All three domain axes converged even tighter than the pre-fix sweep**
(y: ±0.0004 vs ±0.001; z: ±0.0003 vs ±0.004; x: still exactly 0.0000) —
resolution is even more clearly the dominant, and now essentially the only,
source of remaining uncertainty.

`res=60` (74.5 min for this run's single point, comparable to the pre-fix
90.1 min) remains the practical production choice on the same cost/accuracy
tradeoff logic as before — not a claim of exact convergence, now stated
with a wider, more honest residual estimate.

#### ⚠ Superseded — pre-fix numbers (WITHOUT waveguide extension), kept for the record

The input waveguide simply ended at the GDS's own edge in this round —
see the top-of-doc box's item 3 for why that's a real physical defect (a
clean-cut Si/oxide facet reflects part of the guided mode on its own).
19 array tasks across four submissions (11668597, 11668881 — the x-domain
fix, 11668965 — the res 50/60 tail extension).

| res (pts/µm) | R | \|ΔR\| vs prev |
|---|---|---|
| 10 | 0.5144 | — |
| 15 | 0.7386 | 0.2242 |
| 20 | 0.8434 | 0.1048 |
| 25 | 0.8962 | 0.0528 |
| 30 | 0.9266 | 0.0304 |
| 35 | 0.9466 | 0.0200 |
| 40 | 0.9557 | 0.0091 |
| 50 | 0.9705 | 0.0148 |
| 60 | 0.9789 | 0.0084 |

| series | swept | R values | spread |
|---|---|---|---|
| y @ res=10 | `sim_width` 7.5/8.5/9.5 | 0.5144 / 0.5249 / 0.5154 | ±0.010 |
| y @ res=25 | `sim_width` 7.5/8.5/9.5 | 0.8962 / 0.8954 / 0.8952 | ±0.001 |
| x @ res=25 (corrected for the §2 bug) | `edge_clearance` 0.3/0.5/0.8 | 0.8962 / 0.8962 / 0.8962 | 0.0000 |
| z @ res=25 | `sim_height` 1.5/2.0/2.5 | 0.8962 / 0.8919 / 0.8945 | ±0.004 |

This round's `res=60` (R=0.9789) coincidentally matched the paper's
reported C-band minimum (97.9%) to three decimals — read that agreement as
**superseded, not as evidence the pre-fix pipeline was more accurate**; the
facet-reflection artifact it carried was a real error that happened to
partially cancel against the resolution residual at this one resolution.

### Step 9 — full-band 3-D production run ✅ MEASURED 2026-08-04 (job 11673560, after the waveguide-extension fix)

`examples/high_contrast_grating/python/high_contrast_grating_fdtd.py` now
carries `extend_input_waveguide=True` by default (the top-of-doc box's item
3) alongside `RESOLUTION=60` and **C-band `WL_MIN=1.50`/`WL_MAX=1.60`**. The
notebook was regenerated from the script (single source of logic, per step
5's own convention) to match.

**Job 11673560, COMPLETED in 54:23** (faster than the pre-fix run's 93.8
min, despite the corrected single-wavelength `res=60` sweep point costing
more (74.5 vs 90.1 min) — full-band DFT-decay behavior is not simply "the
slowest wavelength's cost," consistent with the tutorial's own "one
broadband pulse, decay is a property of the whole simulation" framing).

**Full spectrum** (`example_output/high_contrast_grating/high_contrast_grating_fdtd/{smatrix.pkl,reflection.png}`):

| λ (µm) | R | | λ (µm) | R |
|---|---|---|---|---|
| 1.500 | 0.9284 | | 1.555 | 0.9661 |
| 1.505 | 0.9326 | | 1.560 | 0.9679 |
| 1.510 | 0.9368 | | 1.565 | 0.9693 |
| 1.515 | 0.9409 | | 1.570 | 0.9704 |
| 1.520 | 0.9449 | | 1.575 | 0.9711 |
| 1.525 | 0.9488 | | **1.580** | **0.9714 (peak)** |
| 1.530 | 0.9524 | | 1.585 | 0.9713 |
| 1.535 | 0.9558 | | 1.590 | 0.9708 |
| 1.540 | 0.9588 | | 1.595 | 0.9699 |
| 1.545 | 0.9616 | | 1.600 | 0.9686 |
| 1.550 | 0.9640 | | | |

Still a smooth, single-humped resonance, but **shifted and broadened
relative to the pre-fix spectrum**: rises more gradually across almost the
entire band before peaking near the red edge (1.580 µm, not 1.550 µm),
then turns over only in the last two points. `check_passive()` raised no
warning.

**Comparison against the paper (98.9% peak at 1569 nm, ≥97.9% over
C-band) — and against the now-superseded pre-fix run:**

| | this run (post-fix) | pre-fix (superseded) | paper |
|---|---|---|---|
| Peak R | **97.14%** | 98.10% | 98.9% |
| Peak λ | **1.580 µm** (+11 nm red of paper) | 1.550 µm (−19 nm blue of paper) | 1.569 nm |
| Min R, 1.50–1.60 µm | 92.84% at 1.500 µm | 91.80% at 1.500 µm | 97.9% (band undefined) |
| Min R, narrower 1.53–1.565 µm | 95.24% at 1.530 µm | 96.77% at 1.530 µm | 97.9% |

**Read this plainly: the corrected, more physically sound run agrees LESS
well with the paper than the flawed pre-fix run did — 1.76 percentage
points off at peak, versus 0.8 before.** This is not evidence the fix was
wrong (independently confirmed correct: the ε-map render, x-domain
sensitivity going to exactly 0.0000, and every domain-convergence check
tightening rather than loosening — step 8). It means the facet-reflection
artifact the fix removed was *coincidentally* pulling the peak magnitude
closer to the paper's number while pulling its *position* the wrong way —
removing one confounding error does not by itself close every gap, and
pretending otherwise would be worse than reporting the larger gap honestly.

**The peak wavelength did not just shift, it flipped sides**: from 19 nm
blue of the paper's 1569 nm to 11 nm red of it — a genuine ~30 nm change in
resonance position from removing the artifact, not merely a magnitude
correction. That degree of movement, on top of the still-only-partially-
bracketed resolution residual (step 8: asymptote estimated at R≈0.98–1.00,
a wider range than before the fix, not a tighter one), is the strongest
evidence yet that **resolution above `res=60` would move the resonance
*position*, not just its magnitude** — untested here, and now the most
load-bearing open question for taking this comparison further (step 10).

#### ⚠ Superseded — pre-fix production run, kept for the record

Job **11670282**, COMPLETED in 93.8 min, run WITHOUT the waveguide
extension (top-of-doc box item 3) — the input waveguide ended abruptly at
the GDS's own edge, a real physical defect that contaminates a reflection
measurement specifically. Its full spectrum and paper comparison below are
**superseded, not merely re-labelled** — read them as "what the pipeline
produced before a real artifact was removed," not as the device's actual
performance.

**Full spectrum** (from job 11670282's own output, which job 11673560 has since overwritten on disk — these numbers are preserved here as text, not re-derivable from the live files anymore):

| λ (µm) | R | | λ (µm) | R |
|---|---|---|---|---|
| 1.500 | 0.9180 | | 1.555 | 0.9809 |
| 1.505 | 0.9265 | | 1.560 | 0.9803 |
| 1.510 | 0.9353 | | 1.565 | 0.9796 |
| 1.515 | 0.9441 | | 1.570 | 0.9790 |
| 1.520 | 0.9526 | | 1.575 | 0.9783 |
| 1.525 | 0.9606 | | 1.580 | 0.9774 |
| 1.530 | 0.9677 | | 1.585 | 0.9763 |
| 1.535 | 0.9735 | | 1.590 | 0.9748 |
| 1.540 | 0.9777 | | 1.595 | 0.9730 |
| 1.545 | 0.9801 | | 1.600 | 0.9707 |
| **1.550** | **0.9810 (peak)** | | | |

A smooth, single-humped resonance — rising steeply from the blue edge,
peaking near 1.55 µm, falling gently toward the red edge. `check_passive()`
raised no warning.

**Comparison against the paper (98.9% peak at 1569 nm, ≥97.9% over
C-band):**

- **Peak: 98.10% at 1.550 µm vs. the paper's 98.9% at 1.569 nm.** 0.8
  percentage points low, roughly consistent with step 8's own estimated
  ~1% resolution residual at `res=60` — this is the number step 8's
  single-wavelength check (0.9789 at exactly 1.569 µm, matching the
  paper's *minimum*, not peak) already flagged as close-but-not-fully-
  converged. **The resonance also sits ~19 nm blue of the paper's** (peak
  at 1.550 µm vs. 1.569 nm) — the 21-point grid (5 nm spacing) means the
  true sampled peak could be off by up to half a grid step from this
  alone, but 19 nm is more than that; a genuine small blue-shift is more
  likely, and step 8's own resolution sweep was checked at 1.569 µm
  specifically *because* that's where the paper's peak is, not because
  the geometry-under-test peaks there. **Not fully explained — candidates,
  none confirmed:** (a) the ~1% resolution residual could plausibly shift
  the resonance wavelength as well as its peak magnitude, not checked here;
  (b) `HCG_seeded_TO.gds`'s filename ("seeded TO") suggests this may be an
  intermediate seed geometry in the paper's staged optimization pipeline
  (parameter optimization → seeded topology optimization) rather than the
  exact final published design — plausible given how close the magnitude
  match already is, but not confirmed against the paper's own figures; (c)
  small material-model differences (this repo's Sellmeier Si/SiO2 vs.
  whatever the paper used).
- **"Minimum over band": 91.80% at 1.500 µm vs. the paper's 97.9% — NOT a
  fair comparison as stated.** The paper does not give its own C-band
  edges; this script's 1.50–1.60 µm span (chosen to match this repo's own
  `directional_coupler`/`mmi_2x2` C-band convention) is wider than the
  standard ITU C-band (~1530–1565 nm) the paper more likely means.
  Restricted to that narrower window, **this run's own minimum is 96.77%
  at 1.530 µm** — much closer to the paper's 97.9%, and still consistent
  with "peak slightly low and slightly blue-shifted" rather than a
  qualitatively different result. Below ~1.53 µm this device's reflectivity
  drops off faster than the plateau suggests — genuinely useful information
  the paper's own summary statistic ("min over C-band") would not surface,
  and worth keeping in mind if `WL_MIN`/`WL_MAX` are ever narrowed to match
  a stricter C-band definition.

**Overall: a real, independently-run 3-D FDTD pipeline reproduces the
published device's peak reflectivity to within ~1 percentage point and its
overall band shape qualitatively, with two specific, stated, unconfirmed
discrepancies (peak wavelength, band-edge definition) rather than a vague
"close enough."** That is the standard this repo's own culture asks for —
a number with a source, not a match asserted from a plot. **All of this
paragraph describes the superseded, pre-waveguide-fix run — see step 9's
main section above for the corrected result, which agrees LESS well with
the paper, not more, once the artifact this paragraph didn't know about
yet was removed.**

### Step 10 — what remains ❌ NOT DONE

1. ~~The RE-RUN full-band production run landing, and a fresh comparison
   against the paper's own spectrum.~~ ✅ **DONE — step 9.** The corrected
   comparison is the one to quote (peak 97.14% at 1.580 µm) — **NOT** the
   superseded pre-fix number (98.10% at 1.550 µm), which agreed with the
   paper *better* only because of the facet-reflection artifact. What that
   comparison itself opened, carried forward:
   - ⚠ **Resolution above `res=60` likely still moves the resonance
     *position*, not just its magnitude — PARTIALLY CHECKED 2026-08-05,
     still open.** The peak moved ~30 nm (from 19 nm blue of the paper's
     to 11 nm red of it) purely from removing the waveguide-facet artifact
     at fixed `res=60`. A `res=70`/`80` single-wavelength check (job
     11688701, see step 8's current-numbers table) shows R at the fixed
     1.569 µm probe wavelength is STILL climbing by a non-trivial amount
     through res=80 (+0.0083 total from res=60), with a deceleration ratio
     that stalled rather than continuing to shrink cleanly — evidence the
     full-band lineshape, and plausibly the peak position too, has likely
     NOT settled by res=60. **This makes the question MORE open, not
     less** — a real full-band production re-run at res=70/80 (not just
     the single-wavelength probe done here) is what would actually confirm
     or rule out further peak movement, and has not been attempted (cost:
     the res=80 single-wavelength point alone took 262 min; a 21-wavelength
     production sweep at that resolution would cost proportionally more).
   - **Is `HCG_seeded_TO.gds` the paper's final published design, or an
     intermediate seed from the staged optimization pipeline?** Still
     unconfirmed against the paper's own figures/text; a ~1.8 percentage
     point gap (larger than the pre-fix ~0.8 pt) is easier to attribute to
     "this is a seed, not the final optimized device" than the smaller
     pre-fix gap was, but that's a plausibility argument, not a check.
   - Small material-model differences (this repo's Sellmeier Si/SiO2 vs.
     whatever the paper used) remain an unconfirmed, unranked candidate.
2. **The resolution residual is now a wider, more honestly-stated range
   (~R 0.98–1.00 extrapolated, step 8) than the pre-fix estimate (~0.99) —
   not narrowed down.** `res=70`/`80` would tighten this bracket; not done
   here on a cost basis (`res=60` alone is already 75–90 min/point, and
   this session already ran the sweep twice end-to-end).
4. **Where does `1 − R` go?** This example makes no attempt to separately
   measure out-of-plane radiated power vs. residual forward leakage vs.
   absorption — `check_passive()` only bounds R itself. A flux-box or
   top/bottom-flux-plane diagnostic (3-D only) would need its own addition
   to `solve_reflection()` or a second, separate `mp.Simulation` call;
   scoped out of this session deliberately rather than half-built.
5. **Generalize `plot_hcg_domain()`** (step 5) into `beef.plot` if a second
   `solve_reflection()`-style device shows up — one example is not enough
   evidence for the right general API yet. Now also draws the extended
   waveguide and top/bottom PML (step 5 update), so the generalization
   surface is a bit larger than when first scoped.
6. ~~None of the three mistakes in the top-of-doc box has an automated
   guard.~~ ✅ **DONE 2026-08-04 for mistakes 2 and 3** — `beef/solvers/
   fdtd.py` gained `_point_in_polygon()` and `check_reflection_layout()`,
   called from `solve_reflection()` right after `src_x`/`mon_x` are
   computed and before `mp.Simulation` is built (so a bad layout fails in
   milliseconds, not after a full solve). It checks, against the exact
   `(2, 0)`-layer polygons `to_meep_prisms()` is about to extrude:
   - the source and monitor both land inside real silicon (guards the
     bug-2 class: `source_margin`/`monitor_margin` are PML-relative, not
     geometry-relative, so a source/monitor position that looks fine in
     isolation can land in bare oxide once `edge_clearance` or any other
     domain-size parameter changes);
   - when `extend_input_waveguide=True`, the silicon actually reaches the
     domain boundary on the excited side (guards the bug-3 class: a
     geometry could claim to honor `extend_input_waveguide` while silently
     ignoring `stub_length` for its input side, leaving the exact cut
     facet the flag exists to remove).

   `tests/test_fdtd_reflection_layout_guard.py` (7 tests, no meep import
   needed — pure point-in-polygon geometry) locks in both checks,
   including reconstructions of the *actual* historical bug-2 and bug-3
   numbers from this device (not just synthetic cases). Separately
   verified against the real production `HCG_seeded_TO.gds` layout
   (`res=60`, `edge_clearance=0.3`, `extend_input_waveguide=True`): passes
   clean; and confirmed it *would* have caught bug 2 (source at
   `edge_clearance=0.8` with unadjusted `source_margin`, the real historical
   `R=0.7530` case) and bug 3 (checked against the pre-extension polygons)
   had it existed then. Full suite: 625 passed, 19 skipped (was 618/19
   before this addition).

   **Mistake 1 (wrong band) still has no automated guard** — a requested
   wavelength range being "near" a device's own design intent isn't
   something this codebase has any representation of to check against;
   left as a genuinely open gap, not attempted here.
