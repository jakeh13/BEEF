# HCG Fabry-Pérot resonator — build plan

Opened **2026-08-05**, at the user's request: an FDTD/EME **hybrid** example
of a Fabry-Pérot cavity terminated by two high-contrast-grating (HCG)
mirrors, tapped in the middle by a directional coupler (DC), per the
diagram supplied with the request. Reproduced in text form since the
source is a user-supplied PDF, not a repo file:

```
 "input/        ______________                    ______________
  reflection"  /              \                  /              \   "output"
      --------|                \________________/                |--------
               \  (DC lower arm, coupled straight body)          /
                \______________/                  \______________/
   FDTD    EME          FDTD          EME (DC)         FDTD          EME    FDTD
  (HCG)  (waveguide)   (DC bend)     (coupled pair)   (DC bend)   (waveguide) (HCG)
   |<---------------------------- Fabry-Pérot cavity (upper rail) --------------------->|
```

The **upper rail** is the resonant cavity itself: HCG mirror — straight
single-mode waveguide — DC's "through" arm (bending down into the coupling
region and back up) — straight single-mode waveguide — HCG mirror, with no
break in the guided path. The **lower rail** is the DC's other two ports
("input/reflection" and "output") — physically disconnected from the
cavity, coupled to it only evanescently across the DC's straight coupled
region. Light injected at "input/reflection" partially cross-couples into
the cavity, bounces between the two HCG mirrors (a genuine Fabry-Pérot
resonance, since both ends are high-but-not-unity reflectors), and leaks
back out at both "input/reflection" (net cavity reflection) and "output"
(net cavity transmission, tapped off evanescently on every round trip).

> **This doc is a PLAN, not a gap list**, same convention as
> `mmi_1x4_sin300_validation.md`: nothing below has been simulated for this
> combined device. Every quoted number is either read directly from
> existing code (labeled **MEASURED (read from source)** — these are
> facts about what already exists, not about this new device) or is a
> **PROPOSED** design choice / **ARGUED** prediction to be checked once the
> pieces actually run. Converts to a normal gap list against the
> [four checks](README.md#the-four-checks) once a hybrid result exists.
>
> **Two pieces already exist and already work**:
> `examples/high_contrast_grating/python/high_contrast_grating_fdtd.py`
> (the HCG mirror, `FDTDSolver.solve_reflection()`, see
> `high_contrast_grating_validation.md`) and
> `examples/directional_coupler/python/directional_coupler_hybrid.py` (the
> DC, `FDTDSolver.solve_section()` + `EMESolver.solve_section_length_sweep()`).
> **Nothing currently cascades independently-solved devices together** —
> that composition layer (`beef.network`) exists and is unused by any
> example script today. This doc's real job is wiring, not new physics.
>
> **Updated 2026-08-05 — Step 0 decided, Step 1 pinned.** The platform
> mismatch (A1) is resolved: reuse the HCG mirror on `SOI_160` as-is,
> redesign the DC for `SOI_160` at a REAL fabricated coupler's own
> dimensions (`duty/dc.gds`, read directly with `gdstk`, fed into the
> existing parametric `DirectionalCoupler` class — not a `GdsGeometry`
> import). **That real device's cross-section supports 4 supermodes, not
> 2** — a measured finding — **but a follow-up measurement (exciting the
> real S-bend with mode 1, decomposed at the coupled-region interface)
> shows only 0.0085% of the launched power reaches the 2nd-order-mode
> pair.** User's call: treat the device as single-mode, `N_MODES=2`, a
> measured simplification, not an assumed one. See Step 0. Step 1 (cavity
> length/FSR) is pinned at 20 µm/arm (~9.5 nm FSR estimated from a real
> measured group index) — provisionally, per the user's own framing:
> revisit once Step 4 actually simulates something and the real resonance
> spacing can be judged directly. **Step 2 is running**: job 11691379
> measures the real coupling at the GDS's fixed 3.0 µm `coupling_length`,
> `N_MODES=2`, at `FDTD_RESOLUTION=40`/`EME_RESOLUTION=80` ASSUMED from the
> `SOI_220` script (user's explicit instruction — not re-derived for this
> platform/device this pass). **While that runs: Steps 3-5's DESIGN is
> now done and verified with synthetic data** (`tests/test_fabry_perot_
> network_crosscheck.py`) — including a real bug caught before any real
> device numbers touched it: the connecting-waveguide block cannot be
> built via `beef.network.convert.from_transfer_matrix()` (one-directional
> by design; this cavity needs the waveguide to carry the mirror's
> reflection back through itself) — fixed with a direct bidirectional
> `SMatrixBlock`. The full 5-block `Circuit` wiring and a hand-derived
> Fabry-Pérot-with-DC-tap closed form are cross-checked against each other
> and against a brute-force solve, all to floating-point precision.
> **Step 3 then landed for real too**: a real 21-point fundamental-mode
> `n_eff(λ)` sweep for the connecting waveguide (matching the HCG's own
> production grid), connected to the HCG's own real mirror block to give
> `rho(λ)` — the actual input Step 5's closed form needs. Only Step 2's
> real `t(λ)`/`k(λ)` are still missing before a fully real cavity spectrum
> can be evaluated.

---

## Summary — where the work actually is

| | status |
|---|---|
| HCG mirror, FDTD reflection (complex, with phase) | ✅ already exists — `FDTDSolver.solve_reflection()`, `SOI_160`, C-band |
| Directional coupler, FDTD/EME hybrid | ✅ exists on `SOI_220`; **redesign for `SOI_160` decided (Step 0), not yet re-solved (Step 2)** |
| Both on the same material platform | ✅ **decided 2026-08-05** — both `SOI_160`, DC using `duty/dc.gds`'s real dimensions |
| DC cross-section is single-mode | ⚠ **cross-section supports 4 supermodes, but MEASURED <0.01% of excited power reaches the 2nd-order pair** — treated as `N_MODES=2` by user decision (Step 0), a measured simplification |
| Cavity geometry (mirror spacing, DC coupling ratio, target finesse) | ✅ **pinned, provisionally (Step 1)** — coupling_length fixed by the real GDS (3.0 µm); waveguide length 20 µm/arm (~9.5 nm FSR, from a measured n_g); revisit after Step 4 simulates something real |
| Cross-device S-matrix composition (mirror + waveguide + DC + waveguide + mirror) | ✅ **wiring VERIFIED with synthetic data (Step 4)** — `beef.network.Circuit`'s first exercise composing >1 independently-solved device (A2); real device numbers not substituted yet |
| A 1-port self-reflecting block (the HCG, in network terms) | ✅ **built and verified with synthetic data (A3/Step 4)** — a direct `SMatrixBlock` construction, not via `from_smatrix()` (which cannot represent it); real HCG spectrum not substituted yet |
| The connecting-waveguide block | ✅ **real bug found and fixed (Step 3)** — `from_transfer_matrix()` is one-directional by design and silently zeroes the mirror's reflection; fixed with a direct bidirectional `SMatrixBlock`. Not yet solved for the real cross-section |
| Closed-form sanity check (Airy/Fabry-Pérot formula) | ✅ **derived and verified two ways (Step 5)** — against a brute-force 4-port solve and against the real `Circuit` wiring, both to floating-point precision; real device numbers not substituted yet |
| Full-device FDTD control run | ⚠ **likely infeasible at any interesting cavity length** — see A4 |

So the ordered work is: **resolve the platform mismatch (Step 0, ✅ done)
→ pin cavity geometry (Step 1, ✅ pinned provisionally) → re-derive the DC
for `SOI_160` at `N_MODES=2` (Step 2, running) →
build the two EME waveguide phase blocks (Step 3, design fixed and
verified, real mode-solve still open) → assemble with
`beef.network.Circuit` (Step 4, wiring verified with synthetic data) →
closed-form cross-check (Step 5, derived and verified with synthetic
data) → the
four checks (Step 6).**

---

## Part A — findings that set the plan

### A1. The two existing devices are on incompatible material platforms — MEASURED (read from source), blocking

`examples/high_contrast_grating/python/high_contrast_grating_fdtd.py` uses
`SOI_160` (`beef.materials.SOI_160`, 160 nm Si core) with a 500 nm-wide
access waveguide (`GdsGeometry(..., input_width=0.5)`, matching the GDS's
own drawn cross-section). `examples/directional_coupler/python/
directional_coupler_hybrid.py` uses `SOI_220` (220 nm Si core) with a
450 nm single-mode strip (`waveguide_width=0.45`).

A single-mode strip's cross-section is not a free parameter once the core
thickness is fixed — 450 nm is single-mode at 220 nm thickness, not
necessarily at 160 nm (and vice versa for 500 nm at 220 nm). **There is no
version of this combined device today where the mirror and the coupler
agree on what "the waveguide" physically is.** Three ways to resolve this,
none simulated yet:

1. **Redesign the DC for `SOI_160`** at a new single-mode width, and re-run
   `directional_coupler_hybrid.py`'s whole beat-length characterization on
   the new platform (this repo's own convergence/resolution numbers for the
   DC are `SOI_220`-specific — `FDTD_RESOLUTION=40`/`EME_RESOLUTION=80`
   were measured on that stack, per that script's own comments referencing
   `dc_validation.md item 11`; nothing says they transfer).
2. **Redesign the HCG mirror for `SOI_220`** at the DC's existing 450 nm
   width, and re-run the HCG's own convergence sweep on the new platform
   (`high_contrast_grating_validation.md`'s `res=60`/PML-clearance/etc.
   numbers are `SOI_160`-specific to this exact device's ~7 µm domain —
   nothing says a 220 nm core changes them, but nothing says it doesn't
   either).
3. **Keep both platforms and add a thin, lossless-in-principle mode
   converter** between the HCG's 160 nm/500 nm access waveguide and the
   DC's 220 nm/450 nm strip — physically realistic (real photonic chips do
   this at layer transitions) but adds a THIRD FDTD simulation type this
   doc's Part B does not currently budget for, and a new loss/reflection
   source at exactly the interface the whole point of this exercise is to
   characterize cleanly.

**This doc does not pick one.** Option 1 is the smallest number of new
FDTD/EME solves (the HCG mirror is reused as-is); option 2 keeps the DC's
already-measured convergence numbers but requires re-running the HCG's own
resolution/PML sweep (recent, expensive — see
`high_contrast_grating_validation.md` step 8, ~75-90 min/point at
`res=60` alone) on a platform nobody has tested it on. **Recommend option
1** on cost grounds alone, but this is the user's call, not something to
assume silently and discover wrong after both hybrid halves are built.

### A2. `beef.network` is the right composition layer, and no example uses it across devices yet — MEASURED (read from source)

`beef/network/` (`SMatrixBlock`, `connect()`, `close_loop()`,
`close_loop_pair()`, `Circuit`) is a general dense-multi-port S-parameter
composition layer, distinct from `beef.results.SMatrix` (a flat, sparse
per-excitation dict with no notion of a port's own reflection —
`block.py`'s own docstring is explicit about this). Every existing hybrid
example (`directional_coupler_hybrid.py`, `ring_resonator_hybrid.py`,
`mmi_2x2_hybrid.py`, …) composes sections **within one device** via either
`HybridSolver`'s own matrix-chain multiplication or (for the ring's
cascaded variant) a hand-rolled closed-form assembly — never via
`beef.network.Circuit` connecting the *outputs of independently-solved
devices*. `tests/test_network.py`/`test_network_ring_crosscheck.py` prove
the primitives are individually correct (including a cross-check that
`close_loop()` reproduces the ring's own closed-form all-pass formula
exactly). One non-test script,
`examples/verification/python/verify_ring_network_crosscheck.py`, already
feeds `close_loop()`/`connect()` real `FDTDSolver.solve_section()`
coefficients — but only for the ring resonator's own coupler section, a
single physical device, and only via the free functions directly, never
the `Circuit` class. **Nothing anywhere uses `Circuit.add_block()` to
compose more than one independently-solved physical device together, and
nothing feeds this layer a genuinely 1-port self-reflecting block (A3).
This device is the first real user of both.**

`Circuit` is exactly this doc's target shape: `add_block(name, block)` for
each of the five physical elements (HCG-left, waveguide-left, DC,
waveguide-right, HCG-right), `connect(name_a, port_a, name_b, port_b)` for
each of the four internal seams, and `reduce(external_ports=...)` to get
back a single `SMatrixBlock` exposing exactly "input/reflection" and
"output" — no loop-closing needed, since the cavity here is a **finite,
mirror-terminated line**, not a closed ring (contrast
`ring_resonator_hybrid_cascaded.py`, which *does* need `close_loop()`
because its ring genuinely closes on itself).

### A3. The HCG mirror is a 1-port self-reflecting block — not what `from_smatrix()` assumes — MEASURED (read from source)

`beef.network.convert.from_smatrix()`'s own docstring: *"Self-terms
(input_port's own reflection, and cross-terms between two different output
ports) are zero-filled — an honest limitation of what SMatrix records, not
an approximation this function introduces."* It requires `input_port` and
`output_ports` to be given as **distinct** labels and always zero-fills the
excited port's own self-term.

`FDTDSolver.solve_reflection()`'s `SMatrix` has exactly one entry: `port="in",
in_mode=1, out_mode=1` — the *complex* reflection coefficient
`r(λ) = bwd/fwd` (phase included, confirmed by reading
`beef/solvers/fdtd.py`'s `solve_reflection()` directly: it stores the raw
ratio, not `|r|²`). This is precisely a port reflecting into **itself** —
the thing `from_smatrix()` is explicitly not built to carry. Routing it
through `from_smatrix()` with a different `input_port` label would
misrepresent the mirror as a 2-port transmissive device (power going in
one labeled port and out a different one, with **zero** self-reflection at
either) — exactly backwards from what a mirror does.

**Not a library bug** — `from_smatrix()` is honest about the limitation —
but a real integration step: build the HCG's block directly,
bypassing `from_smatrix()`:

```python
r = hcg_smat.coeff(port="in", in_mode=1, out_mode=1)   # [n_wl] complex
hcg_block = SMatrixBlock(
    wavelengths=hcg_smat.wavelengths, ports=("p",),
    modes_per_port={"p": 1}, data=r[:, None, None],
)
```

`SMatrixBlock`'s own validation (a dense `[M, N, N]` array, `N =
sum(modes_per_port.values())`) places no lower bound on port *count* — a
genuine 1-port block is structurally valid, just untested by
`tests/test_network.py` today (worth a small new test alongside whatever
implements this, per Step 4/6 below).

### A4. A full-device FDTD control run is likely infeasible — ARGUED (no domain size pinned yet, so no real cost estimate)

Every other hybrid device in this repo (DC, ring, MMIs) has, or plans, a
full-device FDTD run as an independent cross-check on the hybrid
decomposition (check 4 of the [four checks](README.md#the-four-checks)).
For a Fabry-Pérot cavity, the whole point of the device is a HIGH-Q
resonance — `ring_resonator_fdtd.py`'s own validation doc already reports
**167 h estimated for a 3-D `res=40` full-device FDTD** on a racetrack ring
with a much smaller total path length and far lower finesse than a
double-HCG-mirror cavity is likely to have. A Fabry-Pérot cavity's `Q`
scales with mirror reflectivity — the HCG's own measured R≈97% per mirror
(`high_contrast_grating_validation.md` step 9) implies a finesse in the
tens, meaning the cavity's `decay_tol`-driven ringdown time in a full-device
FDTD run would be tens of round trips long, at whatever cavity length Step
1 picks. **This should be measured, not assumed impossible outright** —
state it as a real open cost question in Step 6, not silently skip the
check.

---

## Part B — the ordered steps

### Step 0 — resolve the platform mismatch (A1) — ✅ DECIDED 2026-08-05

**Decision (user): option 1 from A1** — reuse the HCG mirror as-is on
`SOI_160`, redesign the DC for `SOI_160` too, and use a REAL fabricated
coupler's own dimensions rather than a freshly-derived cross-section:
`duty/dc.gds` (top-level `duty/`, not `BEEF/examples/gds_files/` —
different from the HCG's own GDS). **Dimensions only, not a `GdsGeometry`
import** — read directly from the GDS's layer-(2,0) polygons with
`gdstk`, then fed into the existing parametric `DirectionalCoupler` class
(matching every other hybrid example's approach; nothing about this class
requires or benefits from a raw-polygon import path):

| parameter | value | how it was read |
|---|---|---|
| `waveguide_width` | 0.5 µm | coupled-region arm width, both arms |
| `gap` | 0.25 µm | coupled-region edge-to-edge separation |
| `coupling_length` | 3.0 µm | coupled-region x-span |
| `bend_length` | 10.0 µm | each S-bend's x-span (left and right, confirmed mirror-symmetric) |
| `waveguide_offset` | 4.0 µm | far-port centre-to-centre pitch (2.375 − (−1.625), confirmed symmetric bend shapes) |
| `bend_type` | `"euler"` (assumed) | not re-derived from the GDS's curve family — a modeling choice, not a measured one; matches the existing `SOI_220` script's own default |
| `waveguide_height` | 0.16 µm (`SOI_160`) | the GDS itself carries no thickness — supplied externally, same as the HCG needed `input_width`/`input_y` supplied externally |

Verified against `beef/geometry/coupler.py`'s exact constructor semantics
(`gap`=edge-to-edge, `waveguide_offset`=far centre-to-centre pitch,
`coupling_pitch`=derived=`waveguide_width+gap`) before reading these
numbers — confirmed `DirectionalCoupler(waveguide_width=0.5, gap=0.25,
coupling_length=3.0, bend_length=10.0, waveguide_offset=4.0)` builds
without error (`total_length=23.0 µm`, `coupling_pitch=0.75 µm`).

**Real finding that changes the plan: this waveguide is NOT single-mode at
160 nm in C-band — MEASURED** (cheap MPB mode-count check, no FDTD, ~9-11
min on the login node at `res=100`, `OMP_NUM_THREADS=1`):

Bare 0.5 µm-wide, 0.16 µm-thick strip, guided modes (real n_eff above the
SiO2 cladding index):

| λ (µm) | mode 1 (fundamental) | mode 2 (2nd order) | cladding n |
|---|---|---|---|
| 1.50 | 2.2113 | 1.5188 | 1.4446 |
| 1.55 | 2.1507 | 1.5039 | 1.4440 |
| 1.60 | 2.0912 | 1.4926 | 1.4434 |

Both modes sit with real margin above cladding at every C-band point
checked — this is genuinely 2-mode, not a marginal/near-cutoff case.
Coupled-pair cross-section (both arms present, `gap=0.25`) at λ=1.55 µm:

```
n_eff = [2.1628, 2.1366, 1.5396, 1.4500, 1.4430, 1.4327]   (cladding n=1.4440)
```

4 bands sit clearly above cladding (the fundamental-mode even/odd
supermode pair near 2.15-2.16, and the 2nd-order-mode even/odd supermode
pair near 1.45-1.54); band 5 (1.4430) sits at/just below the cladding
index (1.4440) — the same "marginal, drop it" call the `SOI_220` script's
own N_MODES comment makes for its own band 5. **So the cross-section
itself supports 4 supermodes, not 2.**

**But a follow-up real measurement shows the S-bend barely excites the
2nd-order-mode pair at all — MEASURED, 2026-08-05, `res=25`,
`decay_tol=1e-10`, one real `FDTDSolver.solve_section()` run exciting the
far access port with mode 1 only, decomposed onto all 4 supermodes at the
coupled-region interface**:

| output supermode | power fraction |
|---|---|
| 1, 2 (fundamental-mode pair) | 0.4965, 0.4953 (99.18% combined) |
| 3, 4 (2nd-order-mode pair) | 5.5×10⁻⁵, 2.9×10⁻⁵ (0.0085% combined) |

(total 0.9919 — the small shortfall from 1 is ordinary radiation/PML
loss, not a red flag.) **Less than a hundredth of a percent of the
launched power reaches the 2nd-order-mode supermodes** — the large
index gap between the two modes (Δn≈0.65 at 1.55 µm) suppresses
bend-induced conversion at this bend length exactly as expected, and now
this is measured, not argued.

**User's call on this finding (2026-08-05): treat the device as
single-mode, `N_MODES=2`, not 4.** The cross-section is genuinely
multimode (bands exist), but the physics this device's own S-bend
actually excites is, to measured precision, single-mode — so the
simpler, cheaper 2-supermode treatment (matching the `SOI_220` script's
own structure, though NOT its specific numbers — different device, do not
reuse `N_MODES=5` or its resolution/coupling-length values) is used from
here forward. This is a measured simplification, not an assumed one: the
0.0085% figure above is the actual justification, re-derived for this
exact geometry rather than inherited from the general "large Δn suppresses
conversion" argument alone.

### Step 1 — pin cavity geometry — ✅ PINNED 2026-08-05 (provisional, may move after Step 4's real simulation)

`coupling_length` is no longer a free knob (Step 0: fixed at the GDS's own
3.0 µm, not re-derived off a 50:50 target — the user's choice to use the
real device's own dimensions applies here too, not just to the
cross-section). What Step 0 does NOT fix:

- **Mirror-to-DC waveguide length: PINNED at 20 µm each side** (user
  decision, 2026-08-05) — symmetric, matching the diagram's own
  left-right symmetry and letting one `EMESolver` solve (Step 3) serve
  both arms rather than two. Reasoning, not a measured optimum: long
  enough to be a genuine propagation section worth its own EME solve (not
  a token stub), short enough to keep the total cavity in the same
  tens-of-µm scale as this repo's other resonant devices (the racetrack
  ring's own `PERIMETER` is ~41 µm at its shipped dimensions). **Explicitly
  provisional** — the user's own framing: use 20 µm for now, revisit once
  the notebook actually simulates something real (Step 4) and the
  resulting resonance spacing/depth can be judged directly rather than
  argued from the estimate below.
- **DC coupling ratio at this device's real, fixed geometry** — NOT a
  design knob here (contrast the `SOI_220` script, which sweeps
  `coupling_length` to find its own 50:50 point). With `coupling_length`
  fixed at 3.0 µm, the actual power-coupling ratio (fundamental-mode
  even/odd split, per Step 0's `N_MODES=2` decision) at this length is an
  OUTCOME of Step 2's FDTD/EME solves, not a target to hit — measure it,
  don't assume 50:50 or any other specific split.
- **Target resonance wavelength / FSR** — not independently pinned; a
  first estimate below, from a REAL measured group index (not a guess).

**FSR estimate — ARGUED, not yet a real EME solve of the full cavity**,
built from the Step 0 bare-strip fundamental-mode data (2-point-plus-one
finite difference across the 3 measured wavelengths, `n_g = n_eff −
λ·dn_eff/dλ`):

| λ (µm) | n_eff (fundamental) | n_g |
|---|---|---|
| 1.50 | 2.2113 | 4.0293 |
| 1.55 | 2.1507 | 4.0122 |
| 1.60 | 2.0912 | 3.9952 |

At the proposed 20 µm/arm: `L_cavity = 2×20 + dc.total_length (23.0) = 63
µm`. `FSR ≈ λ²/(n_g·L_cavity)` at λ=1.55 µm, n_g=4.0122: **FSR ≈ 9.5 nm**,
~10-11 resonances across the 1.50-1.60 µm C-band sweep — a similar density
to the ring resonator's own C-band characterization (`ring_resonator_
hybrid.py` predicts ~10 resonances at its own shipped `PERIMETER`).
Reasonable to characterize with, not so dense that individual resonances
would be hard to resolve at this repo's usual `N_WL_DENSE` interpolation
grid (see `ring_resonator_hybrid_cascaded.py`'s own grid-convergence note).

**This FSR estimate uses only the FUNDAMENTAL mode's n_g** — and per
Step 0's follow-up measurement (<0.01% of excited power reaches the
2nd-order-mode pair), that is now a well-justified simplification rather
than a known gap: the real cavity's resonances should be dominated by the
single fundamental-mode family, matching the single-mode Airy-function
picture Step 5 plans to use as a closed-form check. Not a proof at the
full-cavity level (the 0.0085% figure was measured at one bend, one
excitation, one wavelength) — worth a similar spot-check once Step 4
assembles the real network, but no longer treated as an open, possibly-
invalidating gap the way it was before the measurement above.

### Step 2 — re-derive whichever existing example moves platforms — ⚠ IN PROGRESS 2026-08-05

Under Step 0's decision: redesign `directional_coupler_hybrid.py`'s device
for `SOI_160`, at `dc.gds`'s own real dimensions, **treated as single-mode
(`N_MODES=2`)** per Step 0's follow-up measurement (<0.01% of excited
power reaches the 2nd-order-mode pair). Step 0's mode-count check replaces
item 1 below (no single-mode width to find — the real cross-section is
4-supermode, but the excited physics is measured single-mode):

1. ~~Find a single-mode strip width at 160 nm core thickness~~ — **N/A,
   superseded by Step 0**: the real device's own width (0.5 µm) supports 2
   guided modes per arm, but the S-bend's own measured mode conversion
   into the 2nd-order pair is negligible (0.0085%) — `N_MODES=2` is used
   as a measured simplification, not a redesigned single-mode
   cross-section.
2. ⚠ **`FDTD_RESOLUTION`/`EME_RESOLUTION` convergence study on `SOI_160`:
   NOT re-run this pass — user's explicit instruction: ASSUME the
   `SOI_220` script's own numbers (40/80) are close enough for a first
   result, at `N_MODES=2`.** This is a stated approximation, not a
   measurement — `directional_coupler_validation.md`'s own item 11 is
   precisely the story of a resolution choice being platform-specific and
   silently wrong for months, and this is a different device (different
   width/gap/bend_length/offset, different core thickness) on top of
   that. If this step's own result looks physically implausible, this
   assumption is the first thing to distrust, not the coupling physics.
3. **Real coupling measurement, submitted as job 11691379**
   (`batch/dc_soi160_hybrid.batch`, `gts-sr108-ece`) — 2 FDTD arm
   excitations (`in`/`isolated`) + 1 EME body mode-solve at the GDS's
   fixed 3.0 µm `coupling_length` (not a 50:50 search — see Step 1), the
   same 2×2 transfer-matrix shape the `SOI_220` script's own core physics
   uses, at `res=40`/`80` (item 2's assumption). Script:
   `examples/directional_coupler/python/directional_coupler_soi160_hybrid.py`.
   A login-node plumbing smoke test (`res=8`, `decay_tol=1e-3`) confirmed
   the real code path — geometry build, both FDTD arm solves, and deep
   into the body's phase-calibration logic — runs correctly before
   failing an expected accuracy gate at that deliberately crude
   resolution; not run at production settings on the login node given the
   cost this repo's other 3-D solves at comparable resolution have shown.

Cheaper than originally scoped (Step 0's own quick check, `res=25` on the
login node, ran in ~24 min for one excitation — the full convergence
sweep is still real SLURM-scale work across several resolutions, but at
`N_MODES=2` rather than 4) — should still be scoped as its own SLURM-batch
effort, not folded silently into Step 4's assembly work.

### Step 3 — the two EME waveguide phase blocks — ✅ REAL BUG FOUND+FIXED, REAL DATA COMPUTED 2026-08-05

Per the user's own count (confirmed against A2): **one** `EMESolver`
mode-solve for the waveguide connecting each HCG to the DC, reused for
both the left and right arms under the Step 1 symmetric-cavity assumption
(mirroring how the HCG's own convergence sweep and the DC's own S-bend
both get reused left/right by geometric symmetry rather than re-solved).

**Treated as single-mode, per Step 0's decision.** The bare strip itself
carries 2 real guided modes (Step 0's own check), but this is a plain
STRAIGHT section — no bend, no discontinuity, nothing to convert power
between modes — and both of its neighbors only ever carry mode-1 content
in this device's assembly: the HCG mirror's own `solve_reflection()` only
tracks `input_mode=1`, and the DC is itself now treated at `N_MODES=2`
(Step 2, the fundamental-mode supermode pair only). So the connecting
waveguide only ever needs to propagate the fundamental mode's phase — a
genuinely single-entry `exp(2πi·n_eff,1/λ·L)` phase block, not a
multimode diagonal.

**This plan originally said to build that block via
`beef.network.convert.from_transfer_matrix()` — WRONG, caught with
synthetic data before any real device numbers were involved
(`tests/test_fabry_perot_network_crosscheck.py`).**
`from_transfer_matrix()` is documented as one-directional by design (its
own docstring: *"Self-terms — reflection at either port, or reverse
transmission — are zero-filled, matching what those solvers actually
model: pure forward transmission, no reflection"*) — correct for
`HybridSolver`'s own straight-through cascades, where nothing downstream
ever reflects back through an earlier section. **This cavity's connecting
waveguides must carry light BOTH forward (toward the HCG mirror) and
backward (the mirror's own reflection, through the SAME physical
waveguide)** — using `from_transfer_matrix()` here builds a waveguide
block that cannot carry anything back through itself, and `connect()`ing
a mirror to it reduces to EXACTLY ZERO, not a small error (verified
directly: a synthetic `connect(mirror_block, "p", from_transfer_matrix(...)
, "in")` gives `0+0j`, when the physical answer is `r*p` — one mirror
reflection seen through one waveguide pass).

**Fix**: build the waveguide's `SMatrixBlock` directly, populating BOTH
directions with the same phase (`S[out,in] = S[in,out] = exp(2πi·n_eff,1/λ·L)`,
zero self-reflection at either end) — ordinary EM reciprocity for a
uniform, lossless section, not a new physics assumption. Verified against
a full `Circuit` assembly (5 blocks, the exact wiring Step 4 needs) AND an
independent brute-force 4-port linear solve, both to floating-point
precision — see Step 5.

**Now solved for real, 2026-08-05**: fundamental-mode `n_eff(λ)` for the
real 0.5 µm/160 nm bare strip, at all 21 wavelengths of the HCG's own
production grid (1.50–1.60 µm) — a real `meow`/MPB mode-solve (`res=100`,
num_modes=1, ~8.8 min total on the login node, no FDTD needed since this
is a straight, uniform section). `n_eff` runs 2.2113 (1.50 µm) down to
2.0912 (1.60 µm), monotonic across the band.

**The real HCG mirror block (A3) and this real waveguide block have been
connected** (`connect(hcg_block, "p", wg_block, "out")`, the exact
Step 4 wiring pattern, using REAL numbers on both sides) to produce
`rho(λ) = r(λ)·p(λ)²` — the HCG mirror as seen from the DC's side of one
20 µm arm, at all 21 wavelengths. `|rho|²` reproduces the HCG's own
measured `R(λ)` spectrum exactly (as it must — propagation phase doesn't
touch magnitude), peaking at 0.9714 at 1.580 µm; `arg(rho)` wraps rapidly
across the band (the 20 µm path length accumulates ~170 rad of phase at
these wavelengths, many multiples of 2π between adjacent 5 nm-spaced grid
points — expected for a path this many wavelengths long, not a bug).
**This `rho(λ)` is the actual input Step 5's closed form needs** — only
Step 2's real `t(λ)`/`k(λ)` remain before the full cavity spectrum can be
evaluated with entirely real numbers.

### Step 4 — assemble with `beef.network.Circuit` — ✅ WIRING VERIFIED 2026-08-05 (with synthetic data; real DC/HCG numbers not yet substituted)

```python
from beef.network import Circuit

circuit = Circuit()
circuit.add_block("hcg_left", hcg_block)          # A3's direct 1-port construction
circuit.add_block("wg_left", wg_block)             # Step 3's FIXED bidirectional construction
circuit.add_block("dc", dc_block)                  # Step 2's re-solved DC, as a 4-port block
circuit.add_block("wg_right", wg_block)            # same block reused, per Step 3
circuit.add_block("hcg_right", hcg_block)          # same block reused, per HCG's own left-right symmetry

circuit.connect("hcg_left", "p", "wg_left", "in")
circuit.connect("wg_left", "out", "dc", "through_left")
circuit.connect("dc", "through_right", "wg_right", "in")
circuit.connect("wg_right", "out", "hcg_right", "p")

cavity = circuit.reduce(external_ports={
    "dc.input_reflection": "input_reflection",
    "dc.output": "output",
})
```

**This exact 5-block, 4-`connect()` wiring is now verified**
(`tests/test_fabry_perot_network_crosscheck.py::
test_circuit_matches_closed_form_symmetric_cavity`) against the
independent closed form derived in Step 5, using synthetic mirror/DC/
waveguide blocks — the first real exercise of `Circuit.add_block()`
composing more than one independently-defined physical device (A2) and of
a 1-port self-reflecting block (A3), both previously unexercised anywhere
in this repo. **Real device numbers (Step 2's DC measurement, the HCG's
own production `SMatrix`, Step 3's real waveguide mode-solve) are not yet
substituted in** — this verifies the WIRING is right, not that the
physical answer is right; that still needs the real blocks.

**The DC's own 4-port block is the one piece this plan has not designed in
detail against the REAL solver output.** `directional_coupler_hybrid.py`
(and Step 2's `SOI_160` redesign) only ever computes and stores
`through`/`cross` given excitation from `in`/`isolated` (2 columns of a
would-be 4×4 S-matrix) — enough for the *straight-through cascade* every
existing hybrid example needs. The port RELABELING this device needs
(`in`→`through_left`, the mirrored output section's port→`through_right`,
`isolated`→`input_reflection`, and `isolated`'s mirrored partner→`output`)
still has to be done once Step 2's real `M_in0`/`M_in1` land — get this
mapping backwards and the assembled "cavity" is actually a straight-through
device with the mirrors on the wrong two ports. Draw the section/port
labeling explicitly (a `plot_hybrid_geometry`-style figure, per every other
example's own step 4) **before** wiring the real `Circuit.connect()` calls,
not after.

### Step 5 — closed-form sanity check — ✅ DERIVED AND VERIFIED 2026-08-05 (design only — not yet evaluated with real device numbers)

Every ring-resonator example in this repo has a closed-form check to catch
a wiring/sign/port-order mistake in the network assembly before trusting
it: the ring's own `t + k1·k2·P/(1 - t_ring·P)` all-pass formula,
independently cross-checked against `close_loop()` in
`tests/test_network_ring_crosscheck.py`. This device's equivalent, for a
symmetric cavity (both mirrors reflection `r`, both connecting waveguides
single-pass phase `p`, DC through/cross amplitudes `t`/`k`, one
wavelength):

```
rho = r * p**2                                # mirror + ROUND-TRIP waveguide phase
D = 1 - rho**2 * t**2
self_reflection = k**2 * rho / D              # tap port -> itself (input_reflection or output)
transmission    = t * (1 + rho**2*(k**2 - t**2)) / D   # input_reflection <-> output
```

Hand-derived from the DC's 4-port S-matrix (`L,R` = the two cavity-rail
ports, `I,O` = the two tap ports) plus the two mirror boundary conditions
(`a_L = rho·b_L`, `a_R = rho·b_R`), then **verified two independent ways**
(`tests/test_fabry_perot_network_crosscheck.py`):

1. Against a brute-force direct linear solve of the same 4-port system
   (200 random trials, `atol=1e-9`) — catches an error shared between the
   closed form and any Step-4-style reduction, which a single cross-check
   cannot.
2. Against the actual `Circuit` assembly (Step 4's exact 5-block wiring),
   to `atol=1e-10`, including the reciprocity check
   `S[output, input_reflection] == S[input_reflection, output]`.

Also confirmed passive (`|self_reflection|² + |transmission|² ≤ 1`) for a
physically valid, lossy-mirror case using a genuinely unitary DC
(`t=cos θ`, `k=i·sin θ` — the standard symmetric-coupler quadrature form;
an earlier attempt with arbitrary, non-unitary `t`/`k` phases spuriously
"violated" passivity, which was the test's own flaw, not the formula's —
see that test's docstring for the derived 90°-phase-relationship
condition this device's real, physically-lossless DC will satisfy
automatically).

**Not yet evaluated with real numbers** — `r`, `t`, `k`, and the real
`n_eff` all still need Step 0/2/3's actual measurements substituted in
before this produces a real predicted spectrum, not a synthetic-data
correctness check.

### Step 6 — the four checks

Convert to a normal gap list against
[the four checks](README.md#the-four-checks) once Steps 0-5 land. A4
already flags that check 4 (full-device FDTD control) may need to be
explicitly scoped OUT with a stated cost estimate, mirroring how
`ring_resonator_validation.md` handles its own 167 h full-device FDTD
infeasibility — a documented, deliberate scope decision, not a silently
skipped check.

---

## Measured vs argued — read this before quoting anything above

**MEASURED (read from source or real MPB mode solves, no FDTD)**: A1's
platform mismatch (both scripts' own module-level constants read
directly); A2's characterization of `beef.network`'s existing primitives
and their current lack of a cross-device, `Circuit`-based caller (read
`beef/network/*.py` and grepped every `examples/**/*.py` for
`Circuit`/`connect(`/`close_loop(` — one non-test hit,
`verify_ring_network_crosscheck.py`, using the free functions directly on
one device's own coupler section, not `Circuit`, and not composing
separate devices); A3's `from_smatrix()` limitation (its own docstring)
and `solve_reflection()`'s complex (not just power) coefficient (read
directly from `beef/solvers/fdtd.py`); **Step 0's real dimensions from
`duty/dc.gds`** (`gdstk`-read polygon extents, cross-checked against
`beef/geometry/coupler.py`'s exact parameter semantics); **Step 0's mode
counts** (real `meow`/MPB mode solves, both the bare strip and the coupled
pair, at 3 C-band wavelengths — 4 supermodes exist at the cross-section);
**Step 0's mode-conversion measurement** (a real `res=25`,
`decay_tol=1e-10` `FDTDSolver.solve_section()` run of the actual S-bend,
mode 1 excited, decomposed onto all 4 supermodes — <0.01% reaches the
2nd-order pair); **Step 1's group-index table** (finite-difference from
the same 3 measured `n_eff` points); **Steps 3-5's synthetic-data
verification** (`tests/test_fabry_perot_network_crosscheck.py`, 3 tests,
all passing) — the `from_transfer_matrix()` bug and its fix, the 5-block
`Circuit` wiring, and the closed-form formula are each checked against
either a brute-force linear solve or each other, to floating-point
precision. This verifies the CODE and the ALGEBRA are correct; it does not
yet verify the PHYSICS of the real device, which needs Steps 0/2/3's
actual measurements substituted in.

**PROPOSED / ARGUED, not measured**: Step 1's 20 µm/arm waveguide length
and the resulting ~9.5 nm FSR (a real n_g, but an un-simulated cavity
length choice); the single-mode Airy-formula picture's validity at the
FULL assembled cavity (Step 0's <0.01% figure was measured at one bend,
one wavelength, one excitation — a strong indicator, not a full-cavity
proof); A4's infeasibility claim for a full-device FDTD control run
(reasoned from the ring resonator's own measured 167 h figure and the
HCG's own measured ~97% mirror reflectivity, but no cost has actually been
measured for *this* device at *any* cavity length); Step 2's own
`FDTD_RESOLUTION=40`/`EME_RESOLUTION=80` (assumed from a different device
on a different platform, not re-measured here, per explicit user
instruction).

**Nothing in this document should be read as a physical result for the
combined device.** No FDTD or EME solve of the real connecting waveguide
has landed yet, and Step 2's real DC coupling measurement is still
running. The HCG mirror and the (`SOI_220`) directional coupler are each
independently measured and validated on their own terms, in their own
docs; Steps 3-5's network-composition CODE is now verified against
synthetic data, but the real numbers it will eventually carry are not in
it yet.
