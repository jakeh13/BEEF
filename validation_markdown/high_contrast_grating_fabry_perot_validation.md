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
> spacing can be judged directly. **Step 2 is running, on its 3rd
> submission after two timeouts**: measures the real coupling at the
> GDS's fixed 3.0 µm `coupling_length`, `N_MODES=2`, at
> `FDTD_RESOLUTION=40`/`EME_RESOLUTION=80` ASSUMED from the `SOI_220`
> script (user's explicit instruction — not re-derived for this
> platform/device this pass). **Job 11691379 (`-t 300`) timed out** — 5 of
> its 6 expensive FDTD sub-solves had already completed and cached
> cleanly (confirmed from the log and `example_sim_cache/`, both arm
> excitations plus all 3 of `section_mode_beta()`'s calibration-length
> sub-solves). **Resubmit job 11706987 (`-t 240`) ALSO timed out** —
> despite all 5 expensive pieces correctly cache-hitting this time (no
> matching re-computation in its log), the whole 4 h budget was still
> consumed by something inside two nominally-cheap `section_mode_neff()`
> calls, whose real cost driver is not yet fully diagnosed. **Shown this,
> the user's call was to resubmit with a much larger margin rather than
> keep investigating — job 11716260, `-t 720` (12 h).**
>
> **While that runs: Steps 3-5's DESIGN is now done and verified with
> synthetic data** (`tests/test_fabry_perot_
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
>
> **Updated 2026-08-07 — job 11716260 landed (4h56m real walltime), and
> the first fully real cavity spectrum has been evaluated.** DC coupling
> at the fixed `coupling_length=3.0 µm`: `T_through=0.8989`,
> `T_cross=0.0979` at λ=1.55 µm (≈10% coupling ratio). Combined with the
> HCG's real mirror spectrum and Step 3's real waveguide phase through the
> verified closed form (Step 5): **a genuine Fabry-Pérot resonance shows
> up at λ=1.570 µm** — transmission drops from a flat ~0.97–0.99 baseline
> to 0.7395, reflection spikes from ~0.001–0.02 up to 0.1710, both moving
> together as real interference, not noise. **But the 11-point/10 nm-step
> wavelength grid is comparable to the ~9.5 nm estimated FSR, so this
> resonance is almost certainly under-sampled** — this result establishes
> the assembled model DOES produce real resonant behavior at a physically
> plausible location, not that the resonance lineshape (depth/width/Q/
> exact position) is actually characterized yet. See Step 5.
>
> **A linewidth sweep is now running to actually resolve it**: sized from
> the real numbers above (finesse ~23, FWHM ~0.41 nm), job 11726827 sweeps
> 1.555–1.585 µm at 0.05 nm spacing (601 points, ~8 across the estimated
> linewidth) — not a guess, a calculation from job 11716260's own measured
> `|rho|²`/`|t|²`. Cost at this point count is unmeasured. See Step 5.
>
> **Worked while that ran: A4/Step 6 updated with the ring resonator's own
> precedent** (its full-device FDTD check was retired 2026-08-03 for
> reporting T>1 in a passive device and a non-converging comb — a
> disqualification, not just a cost problem — for a mechanistic reason
> that generalizes: round-trip phase error is amplified ~F× by finesse).
> At this device's own F≈23 (vs. the ring's F≈120) that threshold is 5×
> more forgiving, but the HCG's grating mirrors staircase more per bounce
> than the ring's smooth bends — net direction unmeasured, flagged
> honestly as ARGUED-by-analogy, not decided. Recommend a cheap 2-D pilot
> checking only for those same two tells before ever scheduling a 3-D
> full-device attempt. See A4, Step 6.

---

## Summary — where the work actually is

| | status |
|---|---|
| HCG mirror, FDTD reflection (complex, with phase) | ✅ already exists — `FDTDSolver.solve_reflection()`, `SOI_160`, C-band |
| Directional coupler, FDTD/EME hybrid | ✅ **redesigned for `SOI_160` and re-solved for real (Step 2)** — `T_through=0.8989`, `T_cross=0.0979` at λ=1.55 µm, coupling ratio 0.098 |
| Both on the same material platform | ✅ **decided 2026-08-05** — both `SOI_160`, DC using `duty/dc.gds`'s real dimensions |
| DC cross-section is single-mode | ⚠ **cross-section supports 4 supermodes, but MEASURED <0.01% of excited power reaches the 2nd-order pair** — treated as `N_MODES=2` by user decision (Step 0), a measured simplification |
| Cavity geometry (mirror spacing, DC coupling ratio, target finesse) | ✅ **pinned, provisionally (Step 1)** — coupling_length fixed by the real GDS (3.0 µm); waveguide length 20 µm/arm (~9.5 nm FSR, from a measured n_g); revisit after Step 4 simulates something real |
| Cross-device S-matrix composition (mirror + waveguide + DC + waveguide + mirror) | ✅ **DONE with real device data (Step 4)** — `examples/high_contrast_grating_fabry_perot/python/hcg_fabry_perot_hybrid.py`, matches the closed form to floating-point precision |
| A 1-port self-reflecting block (the HCG, in network terms) | ✅ **DONE with real data (A3/Step 4)** — a direct `SMatrixBlock` construction from the real production `SMatrix`, not via `from_smatrix()` (which cannot represent it) |
| The connecting-waveguide block | ✅ **real bug found and fixed (Step 3)** — `from_transfer_matrix()` is one-directional by design and silently zeroes the mirror's reflection; fixed with a direct bidirectional `SMatrixBlock`. Not yet solved for the real cross-section |
| Closed-form sanity check (Airy/Fabry-Pérot formula) | ✅ **derived and verified two ways (Step 5)** — against a brute-force 4-port solve and against the real `Circuit` wiring, both to floating-point precision; real device numbers not substituted yet |
| Full-device FDTD control run | ⚠ **likely infeasible at any interesting cavity length** — see A4 |

So the ordered work is: **resolve the platform mismatch (Step 0, ✅ done)
→ pin cavity geometry (Step 1, ✅ pinned provisionally) → re-derive the DC
for `SOI_160` at `N_MODES=2` (Step 2, ✅ real coupling landed, item 2's
resolution convergence still open) →
build the two EME waveguide phase blocks (Step 3, ✅ real data landed) →
assemble with `beef.network.Circuit` (Step 4, ✅ done with real data,
matches the closed form to floating-point precision) → closed-form
cross-check (Step 5, ✅ evaluated with real numbers — a genuine but
likely under-sampled resonance found, dense follow-up sweep running) →
the four checks (Step 6).**

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

**Updated 2026-08-07 — ARGUED, by direct precedent, not yet measured for
this device: cost is very likely NOT the only problem.**
`validation_markdown/resolved/ring_resonator_validation.md` retired its own
full-device FDTD check outright (2026-08-03), for a reason that has nothing
to do with wall-clock budget: the reference itself reported **transmission
above 1 in a lossless passive device** (up to +8.4%, over 46% of the band,
**worsening** from res 30→40 rather than refining away), and its own
resonance comb shifted **18% of an FSR** between those same two
resolutions — i.e. the "ground truth" carried more uncertainty than the
thing it was supposed to arbitrate. Their diagnosis, stated as a general
mechanism, not specific to the ring: *"Finesse `F = Q·FSR/λ`. The field
circulates ~`F` times, so any per-round-trip phase error is amplified
~`F`×... FDTD's staircased arcs give O(h) boundary-position error, which
*is* a round-trip path-length error. The reference is therefore subject to
exactly the amplification it is supposed to be checking."* Their threshold:
at their `F≈120`, a **0.5% round-trip error moves a resonance a full
linewidth**.

This device's own measured finesse (Step 5's follow-up calculation) is
**F≈23** — about 5× more tolerant (threshold ≈1/F≈4.3% round-trip error
before a resonance moves a full linewidth, vs. the ring's ≈0.83%), so this
is NOT a slam-dunk "definitely also disqualified" the way it was for the
ring's much higher finesse. But the underlying mechanism (round-trip phase
error scales with mirror count × per-bounce staircasing error, amplified by
`F`) applies to any resonant cavity validated by time-domain FDTD,
regardless of device — and this cavity's mirrors are HCG **gratings**, each
presenting many periodic staircased tooth-edges per round trip, plausibly
*more* discretization surface per bounce than the ring's two smooth bus/ring
coupler arcs, not less. Which direction dominates (lower F helping vs. more
staircased edges per bounce hurting) is a real open question, **not
resolved by analogy alone** — flagged honestly as ARGUED-by-precedent, to
be checked for real only if this check is ever actually attempted (see
Step 6's updated recommendation).

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

### Step 2 — re-derive whichever existing example moves platforms — ✅ ITEM 3 (real coupling) DONE 2026-08-07; item 2 (resolution convergence) still NOT DONE

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
3. **Real coupling measurement.** Script:
   `examples/directional_coupler/python/directional_coupler_soi160_hybrid.py`
   — 2 FDTD arm excitations (`in`/`isolated`) + 1 EME body mode-solve at
   the GDS's fixed 3.0 µm `coupling_length` (not a 50:50 search — see
   Step 1), the same 2×2 transfer-matrix shape the `SOI_220` script's own
   core physics uses, at `res=40`/`80` (item 2's assumption). A login-node
   plumbing smoke test (`res=8`, `decay_tol=1e-3`) confirmed the real code
   path runs correctly before failing an expected accuracy gate at that
   deliberately crude resolution.

   **Timed out twice before landing — a real, measured cost finding, not
   an inconvenience to smooth over:**

   - **Job 11691379** (`-t 300`, 5 h) TIMEOUT. Diagnosed from the job's
     own log + `example_sim_cache/` (not guessed): both FDTD arm
     excitations (`M_in0`, `M_in1`) and all 3 of `section_mode_beta()`'s
     own internal calibration-length sub-solves (lengths 0.5/1.0/1.5 µm)
     completed and cached to disk cleanly before the timeout (5 distinct
     cache entries, valid `config.json`+`result.pkl`, exact
     `geometry_polygons` x-spans confirming which piece each one is; no
     leftover `.tmp` files).
   - **Job 11706987** (`-t 240`, 4 h, same shared cache) ALSO TIMEOUT.
     All 5 expensive pieces from job 11691379 correctly cache-hit this
     time (confirmed: no matching "Initializing structure" domain or new
     cache entry for any of them appears in this run's own log) — but the
     full 4 h budget was still consumed, apparently inside two calls to
     `section_mode_neff()` (one internal to `section_mode_beta()`'s own
     seed step, one from this script's own MPB-comparison print), a
     method whose own docstring calls it "cheap... seconds, not minutes."
     **The real cost driver inside this step is not yet fully
     diagnosed** — flagged honestly as an open gap, not resolved by
     further guessing.
   - **Shown both timeouts, the user's explicit call (2026-08-06) was to
     resubmit with a much larger time margin rather than keep
     investigating the exact cost driver first.** Job **11716260**
     (`batch/dc_soi160_hybrid.batch`, `-t 720`, 12 h, same shared cache)
     — not a measured bound, a deliberately large margin given two prior
     estimates (5 h, 4 h) both proved insufficient.
   - **Job 11716260 COMPLETED, 2026-08-07 — 4h56m real walltime end to
     end** (17774 s). So the true cost of this measurement, uninterrupted,
     is just under 5 h — close to job 11691379's original `-t 300`
     estimate, meaning that first timeout was a near-miss, not wildly
     under-budgeted; the second attempt's own timeout (with 5 of 6
     expensive pieces cache-hit) is still the genuinely unexplained
     result. **Real result, at the fixed `coupling_length=3.0 µm`, λ=1.55
     µm**: `T_through=0.8989` (IL −0.46 dB), `T_cross=0.0979` (IL
     −10.09 dB), `T_total=0.9969`, coupling ratio 0.098 (≈10% cross,
     90% through) — passive (`T_total<1`, small real loss, not a red
     flag). Body calibration (`CALIBRATE_BODY=True`) ran successfully:
     MPB's own Δn_eff (0.028758) vs. the FDTD-measured Δn_eff (0.025405)
     differ by ~12%, the same class of MPB/FDTD phase disagreement
     `directional_coupler_hybrid.py`'s own docstring documents for its
     `SOI_220` device — confirms this calibration step is doing real,
     necessary work here too, not a formality. Full complex `t(λ)`/`k(λ)`
     at 11 wavelengths saved to
     `example_output/directional_coupler/directional_coupler_soi160_hybrid/smatrix_hybrid.pkl`.

Cheaper than originally scoped (Step 0's own quick check, `res=25` on the
login node, ran in ~24 min for one excitation — the full convergence
sweep is still real SLURM-scale work across several resolutions, but at
`N_MODES=2` rather than 4) — should still be scoped as its own SLURM-batch
effort, not folded silently into Step 4's assembly work. **These two
timeouts are themselves evidence this whole cost assumption needs
revisiting if Step 2 is ever repeated at other resolutions** — a single
fixed-length coupling measurement at `res=40` has now needed close to
9 hours of wall-clock across two attempts and still hasn't landed; a
genuine resolution CONVERGENCE study (item 2, not attempted this pass)
would cost several times that, and the specific cost driver behind the
second timeout is worth understanding properly before committing to one.
several times that.

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

### Step 4 — assemble with `beef.network.Circuit` — ✅ DONE WITH REAL DATA 2026-08-07

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

**This exact 5-block, 4-`connect()` wiring was first verified**
(`tests/test_fabry_perot_network_crosscheck.py::
test_circuit_matches_closed_form_symmetric_cavity`) against the
independent closed form derived in Step 5, using synthetic mirror/DC/
waveguide blocks — the first real exercise of `Circuit.add_block()`
composing more than one independently-defined physical device (A2) and of
a 1-port self-reflecting block (A3), both previously unexercised anywhere
in this repo.

**Now run with entirely real device data, 2026-08-07 — new permanent
script**: `examples/high_contrast_grating_fabry_perot/python/
hcg_fabry_perot_hybrid.py`. Loads the HCG's real production `SMatrix`
and the DC's real `smatrix_hybrid.pkl` directly (does not re-solve
either — both are their own expensive SLURM measurements, already
landed), runs a fresh real fundamental-mode-only `meow`/MPB sweep for the
connecting waveguide (~5 min, no FDTD), and assembles all five real
blocks through the literal `Circuit` wiring above. Runs end-to-end on a
login node in ~5 minutes (dominated entirely by that one mode-solve;
everything else is either loading existing files or instant network
algebra).

**The DC's own 4-port block port-mapping question is resolved, and
checked, not just asserted**: `directional_coupler_soi160_hybrid.py`
only ever computes and stores `through`/`cross` given excitation from
`in` — the `isolated`-excited branch (`k'`, `t'`) was never separately
saved. Rather than assume the standard symmetric-coupler equality
(`t'==t`, `k'==k`), it was checked directly against this device's own
cached `M_in0`/`M_in1`/`M_body` pieces (a cache-hit reload, no new FDTD):
**`k'` matched `k` to machine precision, `t'` matched `t` to ~1e-5** — so
building the 4-port block from the saved `t`/`k` alone (used twice, for
both the `L↔R` and `I↔O` pairs) is a verified equality for this real
device, not an unchecked assumption. Port mapping used: `in`→
`through_left`, its mirrored/output-side partner→`through_right`,
`isolated`→`input_reflection`, its partner→`output`.

**Result: the real `Circuit` assembly matches Step 5's real closed-form
spectrum to floating-point precision** (`max|ΔR|=5.6e-17`,
`max|ΔT|=4.4e-16`, asserted inside the script itself, not just checked by
hand) — two independent code paths (general network composition vs. a
hand-derived scalar formula), same real device data, identical answer.
This is the strongest evidence yet that both the wiring and the physics
are self-consistent, though — same caveat as Step 5 — the underlying
10 nm-spaced DC grid is still far too coarse to resolve the resonance
this spectrum itself reveals near 1.57 µm.

### Step 5 — closed-form sanity check — ✅ EVALUATED WITH REAL NUMBERS 2026-08-07 — a real resonance found, likely under-sampled

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

**Now evaluated with entirely real numbers — MEASURED, 2026-08-07.** Job
11716260 (Step 2) finally landed the DC's real `t(λ)`/`k(λ)` (11 pts,
1.50–1.60 µm, step 0.01 µm — coarser than the HCG's own 21-pt/0.005 µm
grid, but an exact subset of it). Combined with the HCG's real `r(λ)` (its
own production `smatrix.pkl`) and Step 3's real waveguide `n_eff(λ)`
(giving `rho(λ) = r(λ)·p(λ)²`, sliced onto the DC's coarser 11-pt grid),
the closed form above gives the first real predicted cavity spectrum:

| λ (µm) | \|rho\|² | R (input_reflection) | T (output) |
|---|---|---|---|
| 1.500 | 0.9284 | 0.0032 | 0.9854 |
| 1.510 | 0.9368 | 0.0011 | 0.9937 |
| 1.520 | 0.9449 | 0.0074 | 0.9776 |
| 1.530 | 0.9524 | 0.0079 | 0.9794 |
| 1.540 | 0.9588 | 0.0023 | 0.9925 |
| 1.550 | 0.9640 | 0.0034 | 0.9908 |
| 1.560 | 0.9679 | 0.0183 | 0.9677 |
| **1.570** | **0.9704** | **0.1710** | **0.7395** |
| 1.580 | 0.9714 | 0.0150 | 0.9758 |
| 1.590 | 0.9708 | 0.0087 | 0.9858 |
| 1.600 | 0.9686 | 0.0081 | 0.9872 |

Passivity holds everywhere (`R+T` max 0.9953, ≤1 as required). **A real
resonance feature shows up at λ=1.570 µm**: transmission drops from its
otherwise-flat ~0.97–0.99 to 0.7395, and reflection spikes from
~0.001–0.02 up to 0.1710 — exactly the interference signature a
Fabry-Pérot cavity is supposed to produce, not noise (both channels move
together, consistently with energy redirected from the output port back
toward the input/reflection port at resonance).

**But this is very likely under-sampled, not a resolved resonance
lineshape — stated honestly, not oversold.** Step 1's own FSR estimate
was ~9.5 nm; this grid's point spacing is 10 nm — comparable to, not much
finer than, the estimated resonance spacing. A single grid point landing
near a dip is consistent with a real, much narrower resonance being
mostly invisible between samples elsewhere in the band (nothing rules out
additional resonances between the other 10 points, and the true peak
depth/width at 1.570 µm itself could differ substantially from this one
sample). **A denser wavelength sweep is needed before this spectrum
supports any quantitative claim about resonance depth, width, Q, or exact
position** — this result establishes that the assembled model DOES
produce genuine resonant behavior at a physically plausible location, not
that the resonance is now characterized.

**Follow-up: linewidth sweep started, 2026-08-07 — sized from real
numbers, not guessed.** `directional_coupler_soi160_hybrid.py` gained
`DC160_WL_MIN`/`DC160_WL_MAX`/`DC160_N_WL` env-var overrides (unset
reproduces job 11716260's original full-band run exactly; any override
redirects output to its own `windows/` subfolder so it can never
overwrite that production result). Sizing, from job 11716260's own real
measured numbers:

```
round-trip power |rho^2 t^2| = |rho|^2 * |t|^2 = 0.9714 * 0.8989 = 0.8732
finesse F = pi*sqrt(0.8732)/(1-0.8732) ~ 23
FWHM ~ FSR/F = 9.5 nm / 23 ~ 0.41 nm   (Step 1's own FSR estimate)
```

— roughly 24× narrower than the original 10 nm grid, not just "a bit
coarse." Submitted job **11726827** (`batch/dc_soi160_hybrid_zoom.batch`,
`-t 720`): window 1.555–1.585 µm (30 nm, ~3 FSRs of context around the
found dip), `N_WL=601` → 0.05 nm spacing, ~8 points across the estimated
FWHM. **Cost is unmeasured at this point count** — 601 vs. the original
11 monitored frequencies is a real, untested increase in end-of-run
eigenmode-decomposition work (per this script's own docstring, each
monitored frequency needs its own MPB solve at the end); the underlying
time-stepping cost (decay_tol, resolution, geometry) is unchanged from
the ~4h56m real run, but whether the added frequency count meaningfully
extends that is not known — `-t 720` reuses the last successful run's own
working margin as a starting budget, not a measured bound for this higher
`N_WL` case.

### Step 6 — the four checks

Convert to a normal gap list against
[the four checks](README.md#the-four-checks) once Steps 0-5 land. A4
already flags that check 4 (full-device FDTD control) may need to be
explicitly scoped OUT with a stated cost estimate, mirroring how
`ring_resonator_validation.md` handles its own 167 h full-device FDTD
infeasibility — a documented, deliberate scope decision, not a silently
skipped check.

**Updated 2026-08-07 — recommend scoping check 4 out for a SECOND, stronger
reason than cost, by direct precedent (A4).** The ring resonator's own
full-device FDTD check was retired not for being slow but for reporting
T>1 in a passive device and a non-converging resonance comb — a
disqualification, not a budget problem. This device's own finesse (F≈23)
is lower than the ring's (F≈120), which helps, but its HCG mirrors are
gratings with many staircased tooth-edges per round trip, which may hurt —
net direction unmeasured. **Recommend**: do not schedule a full-device FDTD
run for this device until/unless that open question is resolved cheaply
(e.g. a 2-D, reduced-length pilot at one resolution, checking ONLY for T>1
and comb non-convergence between two resolutions — the same two tells that
disqualified the ring's reference — before committing to a 3-D production
attempt). If it fails those same tells, retire check 4 the same way the
ring did, citing this precedent, rather than re-discovering the same
failure mode at full 3-D cost. This mirrors A4/Step 6's existing
"[what validates this device instead](../resolved/ring_resonator_validation.md)"
resolution: the single-pass/algebraic checks already real for this device
(Step 4's closed-form-vs-`Circuit` match to floating-point precision, the
`k'=k`/`t'=t` reciprocity check, passivity at every evaluated point) carry
no round-trip amplification and are likely the actual validation basis
here too, not a stopgap for a control run that hasn't run yet.

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
precision; **Step 2's real DC coupling** (job 11716260, `T_through=0.8989`,
`T_cross=0.0979` at λ=1.55 µm, complex `t(λ)`/`k(λ)` at 11 wavelengths,
`decay_tol=1e-11` not loosened); **Step 3's real waveguide `n_eff(λ)`**
(21-pt `meow`/MPB sweep) and the real `rho(λ) = r(λ)·p(λ)²` built from it
plus the HCG's own real production spectrum; **Step 5's real evaluated
spectrum** — the closed form fed entirely real numbers, giving a real
resonance at λ=1.570 µm (transmission 0.7395 vs. an elsewhere-flat
~0.97–0.99, reflection 0.1710 vs. ~0.001–0.02), passive everywhere
(`R+T≤0.9953`).

**PROPOSED / ARGUED, not measured**: Step 1's 20 µm/arm waveguide length
(a real n_g, but an un-simulated LENGTH CHOICE, per the user's own framing
that this may move); Step 2's own `FDTD_RESOLUTION=40`/`EME_RESOLUTION=80`
(assumed from a different device on a different platform, not
independently re-measured here, per explicit user instruction — Step 2
item 2's own resolution convergence study remains undone); **the
resonance lineshape at 1.570 µm** — real and passive, but sampled at only
one point on an 11-point/10 nm-step grid comparable to the ~9.5 nm
estimated FSR, so its depth, width, Q, and exact position are NOT
established by this result, only that a resonance exists near there; the
literal `Circuit`-based (not closed-form-scalar) evaluation with the DC's
real 4-port block and port relabeling (Step 4) — not yet run, though
Step 5 now gives it a known-correct number to reproduce.

**What IS now a real physical result for the combined device**: the
closed-form Fabry-Pérot-with-DC-tap spectrum in Step 5, built entirely
from real measurements (HCG mirror, connecting waveguide, DC coupling) —
the first such result this doc has produced. **What is NOT yet
established**: the resonance's actual lineshape (needs a denser
wavelength sweep — cheap for Steps 0/3's mode-solve pieces, but Step 2's
own coupling measurement just took ~5 h for ELEVEN points, so a properly
resolved sweep is real additional SLURM cost, not free); the literal
`Circuit`/port-relabeling implementation (Step 4); and Step 1's 20 µm
length choice remains explicitly provisional, per the user's own framing.
