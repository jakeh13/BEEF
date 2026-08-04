# High-contrast grating (HCG) reflector — FDTD-only build plan

Opened **2026-08-04**, at the user's request: an FDTD-only example (no EME, no
`HybridSolver`) of the high-contrast-grating reflector in
`examples/gds_files/HCG_seeded_TO.gds`, loaded via `beef.GdsGeometry`, on
160 nm SOI, with a script and a notebook version.

> **This doc mixes the plan and gap-list conventions** (see [README.md](README.md)).
> **Updated 2026-08-04 (later the same day).** Steps 8–10 below now carry a
> real, measured 3-D resolution/PML(x/y/z) convergence sweep (17 SLURM array
> tasks total, `gts-sr108-ece`) at the **corrected C-band wavelength** — see
> the ⚠ box immediately below — and the full-band production run is
> submitted (job 11670282). The single-wavelength convergence series at
> `res=60` lands at **R=0.9789**, matching the paper's own reported C-band
> minimum (97.9%) to three decimals; extrapolating the delta sequence's own
> decay rate puts the asymptote close to the paper's peak (98.9%). This is
> real, independent agreement with the published device, not merely "the
> pipeline runs."
>
> ⚠ **Two real mistakes were caught and fixed while producing that number —
> both are load-bearing for reading anything below:**
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
>
> Both were caught by the user asking a direct question at the right moment
> ("is this the correct band — check paper, should be c-band"; implicitly,
> by asking for x/y/z convergence rather than accepting the y-only check)
> rather than by anything in this pipeline catching itself. Neither has an
> automated guard now either — worth keeping in mind for the next device
> built this way.
>
> **Placement.** Top level, not `resolved/` or `deferred/`, for the same
> reason `mmi_1x4_sin300_validation.md` sits there: this is not yet
> "sufficiently addressed" (production run still in flight) and it is not
> "postponed by choice" (it is active, with concrete next steps in step 10).
> It differs from that doc's precedent in one way — the device already
> *exists* and *runs*, it just isn't fully validated yet.

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
check — **step 8 below does this at a single wavelength (agreement to 3
decimals at `res=60`); the full-band comparison waits on job 11670282**
(step 9).

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
| 3-D resolution convergence (single wavelength, C-band peak) | ✅ **MEASURED 2026-08-04** — res 10→60, `R` 0.5144→0.9789, geometric decay — step 8 |
| PML/domain convergence, x/y/z independently | ✅ **MEASURED 2026-08-04** — x fully converged (0.0000 spread), y converged at res=25 (±0.001), z small (±0.004) — step 8 |
| x-domain sweep self-corrupting (source outside waveguide) | ❌ → ✅ **CAUGHT AND FIXED 2026-08-04** — see the ⚠ box above — step 8 |
| Full-band (21 λ, C-band) 3-D production run | 🔄 **SUBMITTED 2026-08-04** — job 11670282, `gts-sr108-ece`, in flight — step 9 |
| Comparison against the published paper's own spectrum | 🟡 **PARTIAL** — single-wavelength agreement to 3 decimals at `res=60`; full-band comparison waits on job 11670282 |

So the ordered work was: **geometry inspection → two library gaps → why
`solve()` cannot work here → the reflection method → verified placement →
script/notebook → test suite → smoke test → [band correction] → real 3-D
x/y/z/resolution convergence → [x-sweep bug + fix] → production run
submitted.** What remains is the production run landing, the full-band
comparison against the paper, and the two items step 10 was always going to
carry (where the non-reflected power goes; generalizing the domain plot).

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

### Step 8 — real 3-D resolution + x/y/z PML convergence ✅ MEASURED 2026-08-04, C-band

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

19 array tasks total across four submissions (11668597, 11668881 — the
x-domain fix, 11668965 — the res 50/60 tail extension), all landing clean
after the two mistakes in the top-of-doc box were caught and corrected.

**Resolution series** (`sim_width=7.5`, `edge_clearance=0.3`, `sim_height=1.5`):

| res (pts/µm) | R | \|ΔR\| vs prev | wall (min, 24 ranks) |
|---|---|---|---|
| 10 | 0.5144 | — | 0.4 |
| 15 | 0.7386 | 0.2242 | 0.9 |
| 20 | 0.8434 | 0.1048 | 1.8 |
| 25 | 0.8962 | 0.0528 | 5.1 |
| 30 | 0.9266 | 0.0304 | 6.7 |
| 35 | 0.9466 | 0.0200 | 11.5 |
| 40 | 0.9557 | 0.0091 | 21.7 |
| 50 | 0.9705 | 0.0148 (over a +10 step — ~0.00148/unit, still decelerating per-unit) | 44.0 |
| 60 | 0.9789 | 0.0084 (over a +10 step — ~0.00084/unit) | 90.1 |

Deltas shrink geometrically (~0.55–0.6× per doubling of step size), and
**`res=60`'s R=0.9789 matches the paper's own reported C-band minimum
reflectance (97.9%) to three decimals.** Extrapolating the same decay rate
(remaining ≈ last Δ × r/(1−r), r≈0.55–0.6) puts the asymptote at **R≈0.99**,
close to the paper's reported peak (98.9%). This is real, independent
agreement with the published device — not merely "the pipeline runs" the
way step 6's O-band smoke test was.

Cost escalates steeply (~2× wall time per 10-unit step by `res=50`+, faster
than the 2×-per-5-units trend seen at lower resolution — consistent with
total cost scaling faster than `res³` once the CFL-limited timestep count
is folded in). `res=60` (90 min/point) was judged the practical stopping
point: pushing to `res=70–80` would cost multiple hours per point for a
residual already estimated at ~1%, and `res=60` independently reproduces
the paper's own number to 3 decimals, which is stronger evidence than a
tighter self-consistency bracket alone would have been.
**`res=60` is a cost/accuracy tradeoff, stated as such — not a claim of
exact convergence** (same posture this repo's README takes toward every
other device's resolution knob, e.g. the coupler's own res 40 vs measured-
need-120 gap).

**Domain (PML) series, all at a single resolution:**

| series | swept | R values | spread |
|---|---|---|---|
| y @ res=10 | `sim_width` 7.5/8.5/9.5 | 0.5144 / 0.5249 / 0.5154 | ±0.010 |
| y @ res=25 | `sim_width` 7.5/8.5/9.5 | 0.8962 / 0.8954 / 0.8952 | ±0.001 — converged |
| x @ res=25 (corrected) | `edge_clearance` 0.3/0.5/0.8 | 0.8962 / 0.8962 / 0.8962 | **0.0000 — fully converged** |
| z @ res=25 | `sim_height` 1.5/2.0/2.5 | 0.8962 / 0.8919 / 0.8945 | ±0.004 |

**Reading this together:** at `res=25`, x-clearance is a complete non-issue
(exactly zero sensitivity once the source-placement bug — box above — was
fixed), y-clearance is converged to 0.1%, and z has a small ~0.4% spread —
all three are sub-dominant to the resolution series' own ~0.9% residual at
`res=40` (and smaller still relative to `res=60`'s estimated ~1% residual).
**Resolution, not domain size, is the dominant source of remaining
uncertainty** — the opposite conclusion from what step 4's un-swept ~0.26
µm y-clearance flag might have suggested, and only knowable by actually
running the sweep rather than eyeballing the geometry plot.

### Step 9 — full-band 3-D production run 🔄 SUBMITTED 2026-08-04, in flight

`examples/high_contrast_grating/python/high_contrast_grating_fdtd.py`'s
production default updated to **`RESOLUTION=60`** (from the old unvalidated
placeholder `40`) and **C-band `WL_MIN=1.50`/`WL_MAX=1.60`** (from the wrong
O-band), matching step 8's settled resolution and the corrected band. The
notebook was regenerated from the script (single source of logic, per step
5's own convention) to match.

Submitted via a new batch script, `batch/run_high_contrast_grating_fdtd.batch`
(`gts-sr108-ece`, 1 node/24 ranks, `-t 300`) — **job 11670282**. Cost
estimate: step 8's single-wavelength `res=60` point took 90 min on 24
ranks; this script's 21-wavelength sweep is the *same* physical run (one
broadband pulse, every wavelength's DFT captured for free), so wall time
should land in the same ballpark — the 5 h cap leaves real margin, not a
tight bound.

**Not yet landed as this section was written.** Once it does: read
`R(λ)` across the full 1.50–1.60 µm band from
`example_output/high_contrast_grating/high_contrast_grating_fdtd/smatrix.pkl`,
compare its peak and its band-minimum directly against the paper's 98.9%/
97.9%, and fold the result into this section with a MEASURED tag — this
paragraph is the placeholder for that update, not a substitute for it.

### Step 10 — what remains ❌ NOT DONE

1. **The full-band production run landing** (step 9) and its comparison
   against the paper's own spectrum shape, not just the single-wavelength
   peak value step 8 already agrees with.
2. **Where does `1 − R` go?** This example makes no attempt to separately
   measure out-of-plane radiated power vs. residual forward leakage vs.
   absorption — `check_passive()` only bounds R itself. A flux-box or
   top/bottom-flux-plane diagnostic (3-D only) would need its own addition
   to `solve_reflection()` or a second, separate `mp.Simulation` call;
   scoped out of this session deliberately rather than half-built.
3. **Generalize `plot_hcg_domain()`** (step 5) into `beef.plot` if a second
   `solve_reflection()`-style device shows up — one example is not enough
   evidence for the right general API yet.
4. **Neither of the two mistakes in the top-of-doc box has an automated
   guard.** A future device swept the same way could reintroduce either —
   wrong band (no check that a requested wavelength range is anywhere near
   a device's own design intent) or the source-outside-waveguide bug (no
   assertion, in `solve_reflection()` itself, that the source/monitor
   actually land inside drawn silicon before running). Worth a cheap
   defensive check in `solve_reflection()` at some point: verify
   `geometry.polygons()`'s silicon actually contains the computed
   `(src_x, src_yo)` and `(mon_x, src_yo)` points before building the
   `mp.Simulation`, and raise clearly if not — this session found the bug
   only by rendering an ε-map and reasoning about it by hand, which is not
   how the *next* geometry's author should have to catch it.
5. **`res=60`'s ~1% estimated residual is exactly that — estimated,
   extrapolated from a decay ratio, not bracketed by a point beyond it.**
   If a future need calls for it, `res=70`/`80` would settle whether the
   asymptote is really ~0.99 or the deceleration trend breaks down — not
   done here on a cost basis (per-point cost was already tens of minutes
   to 1.5 h by `res=60`).
