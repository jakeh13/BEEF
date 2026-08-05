# 4×4 MMI 90° optical hybrid — SOI 220 nm, O-band — build plan

Opened **2026-08-05**, at the user's request: an optical 90° hybrid based on a
4×4 MMI, on **220 nm silicon-on-insulator** (`SOI_220`), with real **input and
output access tapers**, reproducing the device topology of

> T. Yasui, J. Sugisaka, and K. Hirayama, "Structural Optimization of an
> Optical 90 Degree Hybrid Based on a Weakly-guided 4×4 Multimode Interference
> Coupler Using a Parallelized Real-coded Micro-genetic Algorithm," *ACES
> Journal*, vol. 36, no. 5, pp. 526–532, May 2021.
> https://doi.org/10.47037/2020.ACES.J.360506

Full citation and the two papers it in turn cites that matter here are in
[References](#references) at the bottom.

> **This doc is a PLAN, not a gap list** — same convention as
> `mmi_1x4_sin300_validation.md` and
> `high_contrast_grating_fabry_perot_validation.md`, the two other build
> plans in this directory: nothing below has been simulated, no geometry
> exists yet, and every device *number* in Part B is either taken directly
> from the paper's own general design equations or explicitly marked
> **TBD — measure, don't assume**. Findings about the *codebase* (what
> already supports this device, what doesn't) are labelled **MEASURED (read
> from source)**, exactly as the other two plans do. This converts to a
> normal gap list against the [four checks](README.md#the-four-checks) once
> the example scripts exist and run for the first time.

**⚠ Terminology collision, flagged once so it doesn't cause confusion for
the rest of this file.** "Hybrid" means two unrelated things in this
project:

1. **The device under test** — an *optical 90° hybrid*, the coherent-receiver
   component this whole paper is about.
2. **`HybridSolver` / `*_hybrid.py`** — this codebase's FDTD+EME cascade
   *simulation method*, used for every device in this directory regardless
   of what the device itself is called.

Every other file here is implicitly "the hybrid-method validation for
device X." This file is "the hybrid-method validation for a device that is
*itself* called a hybrid." Read "hybrid" as the simulation method unless a
sentence explicitly says "optical hybrid" or "90° hybrid."

---

**Status 2026-08-05:** Step 1's cross-section sweep (1a single-mode width,
1b sz convergence, 1c L_π vs. body width, all at O-band on `SOI_220`) is
running as SLURM job **11692410** (`batch/verification/
run_mmi_4x4_90hybrid_cross_section.batch`, account `gts-sr108-ece`, 8 cpus,
`-t 240`) — moved to batch after a login-node attempt at
`OMP_NUM_THREADS=1` was still running past 12 minutes on a single 2-mode
point with no output yet; `run_mmi_2x2_gap2_dbeta.batch`'s own comment
records 206 s for one comparable point on that device's (smaller) body
cross-section, so this is genuine batch-scale work, not a hung process.
While that runs, **Step 7's math primitives landed** (see A8) and **A5's
mirror-symmetry claim is now confirmed directly**, not just argued from
source (see A5) — both meep-free, done on the login node.

**Update, same day — 11692410 cancelled, resubmitted corrected as
11693043 (see A10 below).** Partial output before cancellation: even the
sweep's *narrowest* width (0.30 µm) already showed 2 guided modes at O-band
on `SOI_220`, not the single mode a 220 nm Si strip normally supports —
enough of a surprise to investigate before trusting the range, and enough
to trigger a real, recorded dead end (A10): trying to separate the count by
polarization (the paper excites only a TE-polarized fundamental mode) using
meow's `Mode.te_fraction`, then finding that metric disagrees with this
repo's own already-shipped precedent. Reverted; the width range is now
widened downward (0.10–0.70 µm, was 0.30–0.70) and the job resubmitted on
the plain, repo-standard guided-mode count.

> ### A10. `Mode.te_fraction` does not match this repo's own "TE" convention — a dead end, recorded not hidden
>
> **The observation.** A first pass at 1a filtered guided modes by
> `te_fraction > 0.5` (0.0 = pure TM, 1.0 = pure TE, per meow's own
> docstring), reasoning that the paper's device excites only a TE-polarized
> fundamental mode (Sec. IV) and that these narrow, near-square access
> waveguides (width ≈ height ≈ 0.2–0.4 µm) are exactly the aspect ratio
> where a TE0/TM0 pair sits close in n_eff.
>
> **The sanity check, and why it matters.** Before trusting a brand-new
> filter, it was checked against a cross-section this repo already trusts:
> `mmi_2x2_hybrid.py`'s own access waveguide (`width=0.45`, C-band, SOI_220
> — documented throughout this repo, including that script's own module
> docstring reasoning, as running with a TE-polarized fundamental mode).
> Its dominant mode (`n_eff=2.311740`, exactly the index this repo's own
> devices already report and depend on) came back **`te_fraction=0.000199`
> — i.e. meow calls this repo's own working "TE" mode "TM."**
>
> **What this does and doesn't establish.** It does NOT prove meow's
> `te_fraction` is buggy — it may be using an internal axis/field-component
> convention that doesn't line up with this codebase's informal "TE"
> labelling, or the real vector-mode character may be more subtly
> hybridized than the paper's own 2-D scalar TE/TM split (Eq. 2) assumes —
> that split is an effective-index-method reduction this codebase's full
> 3-D vector solve (meow/meep) does not make, so the paper's clean
> TE-vs-TM decoupling does not automatically transfer. What it DOES
> establish: **this metric cannot currently be trusted to separate TE from
> TM in this codebase**, because it disagrees with a case this whole repo
> already treats as settled.
>
> **Why this is not blocking, and why no device has ever needed it.** Every
> existing device in this repo (2×2, 1×2, 1×4, the coupler, the modal
> multiplexer) counts guided modes as plain `n_eff > n_clad` bands, with
> zero polarization filtering, and none of their validation docs raise TE/TM
> conflation as an open risk. The reason is structural, not incidental: a
> plain rib/strip cross-section with no lateral or vertical asymmetry keeps
> TE and TM as exact, non-interacting symmetry classes through any
> symmetric perturbation (taper, MMI body, …) — so an extra wrong-
> polarization band sitting unused in an interface basis costs a wasted
> mode slot, never a real coupling error. This device's access waveguides
> and MMI body have exactly that same symmetry, so the same reasoning
> applies without needing a runtime filter.
>
> **Decision: reverted to the plain, repo-standard guided-mode count.** Not
> because the TE/TM question is unimportant, but because (a) the tool tried
> for answering it does not currently agree with this repo's own ground
> truth, and (b) the repo's own structural argument already covers why it
> doesn't matter for this class of device. If a genuinely asymmetric
> section is ever added to this device (it is not planned to be), this
> reasoning would need revisiting.

## Summary — what's free, what's new

| | status |
|---|---|
| 4-input / 4-output MMI geometry | ✅ **already supported, no new class** — `MMI(n_inputs=4, n_outputs=4, …)`, same constructor every other MMI example uses |
| Input **and** output tapers | ✅ **already supported** — `MMI`'s default port constructor gives every port its own `taper_length`/`taper_width`, on both the input and output side, for free |
| 220 nm SOI at O-band | ⚠ **new (platform, band) *pair*, not new physics** — `SOI_220` and O-band each already exist and are separately validated elsewhere in this repo, but never together (`SOI_220` has only ever been run at C-band; O-band has only ever been run on `SOI_160`) |
| Si/SiO2 material model at O-band | ✅ **already validated, zero new work** — unlike the SiN 1×4 plan, where "the material model was the whole job," this device's material system is the same Sellmeier model the 1×2 MMI already runs at O-band on `SOI_160`. Nothing in `beef/materials/stack.py` needs to change. |
| Mirror symmetry (paper Eq. 8, `y_i = y_{i+4}`) | ✅ **CONFIRMED 2026-08-05** — not just the default behaviour, checked directly on a real `MMI(4,4)` instance — see A5 |
| 4-port aperture isolation (both ends) | ⚠ **needs `use_neighbor_aperture=True`**, not `half_domain_aperture` — see A6, precedent already exists in this repo (`mmi_1x4_hybrid.py`) |
| CMRR / phase-error / imbalance / excess-loss metrics (paper Eqs. 10–15) | ✅ **Step 7 landed 2026-08-05** — `SMatrix.imbalance_db`/`excess_loss_db`/`cmrr_db`/`cross_excitation_phase_diff_deg`, 13 new tests, synthetic data only. The *math* primitive for phase exists now; the *physical* risk A8 flags (do two independently-run FDTD excitations share a consistent phase reference?) is **not** resolved by this and still needs a real check once Step 6 produces actual solves. |
| Self-imaging design equations (paper Eqs. 7–9) | ✅ **reusable directly** — they are general self-imaging relations, not specific to weakly-guided silica — see A3 |
| The paper's micro-GA optimizer | ❌ **not ported** — deliberate substitution, see A3/Step 3. `beef/sweep/` has no optimization loop, only `parallel.py` (a map-over-workers helper) |

So the ordered work is: **settle the (platform, band) cross-section → derive
starting dimensions from the paper's own equations → refine with a cheap EME
sweep (not a GA) → geometry/aperture → interface mode audit → three example
scripts → a new phase-aware metrics helper → the four checks.** Steps 1–2 are
cheap mode solves; step 3 is the one genuine simulation campaign before any
example script is written; nothing here is blocked on a library gap the way
the SiN 1×4 plan was blocked on the material model, or the mux was blocked on
`CompoundGeometry`.

---

## Part A — findings that set the plan

### A1. The geometry is already fully general — MEASURED (read from source)

`beef/geometry/mmi.py:27` — `MMI`'s own docstring: *"n_inputs=1, n_outputs=2
gives a 1x2 splitter; n_inputs=2, n_outputs=2 a 2x2; n_inputs=1, n_outputs=4 a
1x4; etc."* `n_inputs=4, n_outputs=4` is the same constructor path, no special
case anywhere in the class. `_port_ys(n, spacing)` (`mmi.py:74`) centres *n*
ports symmetrically about y=0 for any *n*; `MMIPort` (`mmi.py:16`) already
carries a per-port `taper_length`/`taper_width`, so **both** input and output
tapers are drawn by the default constructor with no extra plumbing — the
"input and output tapers" the user asked for are not a feature to add, they
are what `MMI` already does for every port on every device in this directory.

### A2. `SOI_220` × O-band is a new *pair*, not new physics — MEASURED (read from source)

`beef/materials/stack.py:711-712`:

```python
SOI_160 = MaterialStack(wg_height=0.16, sz=2.0)   # 160 nm Si (default in this codebase)
SOI_220 = MaterialStack(wg_height=0.22, sz=2.0)   # 220 nm Si (standard commercial SOI)
```

Every existing `SOI_220` device (2×2 MMI, directional coupler) runs at
C-band. Every existing O-band device (1×2 MMI, modal multiplexer) runs on
`SOI_160`. **No example in this repo has ever combined 220 nm Si with
O-band.** This repo's own precedent (`modal_multiplexer_validation.md` gap 2:
*"It had to be measured fresh on SOI_160 rather than copied from the
coupler's res 80, since a 160 nm core confines more weakly than 220 nm"*)
treats a new (platform, band) pair as needing its own cross-section check
even when both halves are individually well understood — the same standard
applies here, in the *other* direction: `sz=2.0` and the existing SOI_220
examples' `waveguide_width=0.45` were chosen (or inherited) for **C-band**;
neither is validated at O-band's shorter wavelength, where the same physical
core confines more tightly. **Do not reuse `waveguide_width=0.45` from
`mmi_2x2_hybrid.py` without checking it is still single-mode at O-band** —
tighter confinement at shorter λ makes a mode *more* guided, which is exactly
the direction that can turn a single-mode C-band width multi-mode at O-band.

The Si/SiO2 Sellmeier model itself needs nothing new: `_SELLMEIER["si"]`/
`["sio2"]` already cover O-band (the 1×2 MMI proves it in production), so
this device's material work is zero, unlike the SiN 1×4 plan where the
material model was blocking everything.

### A3. This device's physics regime does NOT need the paper's micro-GA — ARGUED, from the paper's own text

The paper's central problem statement (its own §I): self-imaging theory
"cannot be assumed in weakly guided waveguides," so it needed a real-coded
micro-GA to hit CMRR/phase-error targets on **weakly-guided ZrO2-doped
silica** (Δ = 5.5%, Table 1's optimized `w = 2.98 µm`, `L = 459 µm`). Our
platform is the opposite regime: 220 nm Si core in SiO2 clad is
**strongly-guided, high-index-contrast** — the same Δn² ≈ 10 vs. silica's
much smaller contrast already established in this repo's own SiN-vs-Si
comparison (`mmi_1x4_sin300_validation.md` A7: "Si/SiO2 Δ(n²) = 10.0094" vs.
nitride's 1.8999, and silica's contrast is smaller again). Self-imaging
theory is expected to hold well here — which is exactly the regime the
paper's own reference [4] (Halir et al., *Opt. Lett.* 2011, "High-performance
90º hybrid based on a silicon-on-insulator multimode interference coupler")
addresses: a real SOI 90° hybrid, designed from self-imaging without a GA.

**Decision for this plan: use the paper's own general self-imaging equations
(Eqs. 7–9, general — not weakly-guided-specific) to get starting dimensions,
then refine with a coarse EME width/length sweep (this repo's
`solve_section_length_sweep()` pattern, already used by `mmi_2x2_hybrid.py`
and `mmi_1x4_hybrid.py` for exactly this kind of cheap length refinement) —
not a ported micro-GA.** `beef/sweep/` has no optimizer (`parallel.py` is a
worker-pool map, not a GA), so porting one would be new infrastructure, not
device work. **Flagged as a judgment call, not a silent one** — if the
self-imaging + sweep approach does not hit the paper's own quoted system
requirements (CMRR ≥ 20 dBe, |Δφ| ≤ 5°, from its refs [4,5]) within a
reasonable sweep budget, porting a real optimizer becomes the next honest
step rather than hand-tuning forever.

### A4. Only ports 1 and 3 are ever driven — MEASURED, read from the paper

Paper §III: *"input waveguides 2 and 4 are not used for optical 90º hybrids
[3]"* — `I = {1,3}`, `O = {5,6,7,8}` in the paper's own notation. Mapped to
this codebase's 0-indexed `input_idx`/`output_idx`: only `input_idx ∈ {0,2}`
of the four input ports are ever excited; all four `output_idx ∈ {0,1,2,3}`
are measured for both excitations. A5 below means the other two input ports
(`input_idx ∈ {1,3}`) get solved anyway as a side effect of the mirror-symmetry
reuse — not wasted, just not separately required by the paper's own figures
of merit.

### A5. Mirror symmetry (paper Eq. 8) is already the default, not a feature to build — ARGUED from source + the paper

Paper Eq. (8) plus its own follow-up sentence: *"we assume y_i = y_{i+4} for
i = 1,…,4"* — each input port sits at the **same** y-position as its
same-indexed output port, and (implicitly, since both sides share `w`) the
same access-waveguide width. That is precisely what `MMI`'s default
constructor already gives for free: `_default_ports(n)` builds `self.inputs`
and `self.outputs` from the *same* `waveguide_width`/`taper_length`/
`taper_width`/`waveguide_offset` call (`mmi.py:58-66`), so as long as this
device is built via the plain `MMI(n_inputs=4, n_outputs=4, waveguide_width=…,
taper_length=…, taper_width=…, waveguide_offset=…)` call — no custom
`inputs=`/`outputs=` port lists — the input and output faces are
*structurally* identical, exactly the symmetry the paper's Eq. 8 states.

**This is the direct answer to "check the paper for symmetry between input
and output ports."** It also means this device is a strong candidate for the
same left-right mirror-symmetry reuse `mmi_2x2_hybrid.py` already exploits
(module docstring: *"the output section with output_idx=k is exactly the
x-mirror of the input section with input_idx=k"*) — generalized from 2 ports
to 4:

- Solve the **input** section once per port index (`input_idx = 0,1,2,3` —
  **4 FDTD solves**, not 8), plus **1 EME body solve** (length-swept
  analytically, same pattern as `mmi_2x2_hybrid.py`/`mmi_1x4_hybrid.py`).
- Every output-side transfer matrix is the corresponding input matrix's
  transpose — the output section is **never separately simulated**, same as
  the 2×2's own doc.
- To get the paper's full `S_ji(λ)` table for `j ∈ {5,6,7,8}` (all four
  outputs) and `i ∈ {1,3}` (the two driven inputs), cascade: `M_body @
  M_in_i` then `M_out_j.T @ (that)` for each of the 4 output indices — 4
  input solves cover **all 4×2 = 8** entries the paper's fitness function
  needs, at the cost of 4 FDTD solves total (vs. up to 8 if input/output
  were solved independently, or 2 if only the driven ports were solved and
  outputs were separately simulated).

**Not yet confirmed, flagged as a real risk rather than assumed to hold:**
the reuse is only valid if the *body* is genuinely left-right symmetric too
(a plain rectangular MMI slab, `ΔL` folded symmetrically) and if
`body_stub`/`half_domain_aperture` settings on both ends match exactly —
worth a direct check once geometry exists (Step 4), not before.

> **✅ CONFIRMED 2026-08-05** — checked directly, not just read from source:
> `MMI(n_inputs=4, n_outputs=4, mm_length=20.0, mm_width=10.0,
> waveguide_width=0.45, taper_length=3.0, taper_width=1.2,
> waveguide_offset=1.5, waveguide_height=0.22, use_neighbor_aperture=True)`
> gives `input_ys == output_ys == [-2.25, -0.75, 0.75, 2.25]`, and every
> input/output port pair matches on `y`, `waveguide_width`, `taper_length`,
> and `taper_width`. `frame_offset == 0.0` as well (equal `taper_length` on
> both sides), confirming the body-centred frame used for the mirror
> reasoning above is also bbox-centred here — no hidden offset between the
> "mirror" picture and the actual coordinate frame `section_prisms()`/
> `section_gds()` use. Not yet checked: the *body itself* being
> left-right symmetric once `ΔL`/taper lengths are chosen for real (Step 2/4)
> — this only confirms the port-level claim.

### A6. Aperture: needs `use_neighbor_aperture=True` on BOTH ends — MEASURED (read from source), by direct analogy

`mmi_2x2_hybrid.py`'s `half_domain_aperture` works only because its two
ports per end split cleanly across a shared y=0 axis. This device has **four**
ports per end, on **both** the input and the output side (unlike the 1×4
SiN device, which only had the problem on its output side) — a y=0 cut
still groups ports {1,2} and {3,4} together on each end, exactly the failure
mode `mmi_1x4_hybrid.py`'s own module docstring identifies for its four
*output* ports:

> *"This device's four ports … don't split that way: a y=0 cut still groups
> ports 0&1 … and 2&3 … together. `use_neighbor_aperture` instead sizes each
> port's aperture to the distance to its OWN nearest neighbor."*

**Decision: `MMI(..., use_neighbor_aperture=True)`, and no
`half_domain_aperture=True` calls on either FDTD section** — same mechanism,
now needed on both faces rather than one. `mmi.py:332-355`
(`source_aperture_at`) already implements this generally; nothing new to
write.

### A7. Interface mode count is unknown — must be measured, not assumed

Every MMI hybrid example in this repo picked its own `N_MODES` from a real
mode-count measurement at its own cross-section, never copied across
devices at face value: 2×2 MMI → 14 (of 15, one unguided dropped), 1×4 SiN →
20 (of 23 guided), modal multiplexer → 14 (of 16 requested, cut down after
measuring only 14 are actually reachable). This device's body is wider
still (four access waveguides plus tapers on each end, vs. two or one), so
**expect more supermodes than any MMI example so far**, closer to the modal
multiplexer's territory than to either existing MMI's. `run
audit_interface_modes.py`-style — a real mode solve at the settled
cross-section (Step 5) — before picking a number, and budget for the
duplicate-band risk the modal multiplexer hit at high mode count (the
library-level guards, `_iter_distinct_bands()`/`compute_modes()`'s
`guard_duplicate_modes`, are on by default repo-wide since 2026-07-31, so
this device inherits the fix for free — but the *measurement* of how many
distinct modes actually exist here is still device-specific and unmeasured).

### A8. This is the first device whose own metric needs S-parameter phase, not just power — MEASURED (read from source)

`beef/results/smatrix.py` stores `_data[(port, in_mode, out_mode)]` as
**complex** amplitudes (`SMatrix.set`/`coeff` docstring: *"All coefficients
are complex; power = |coeff|^2"*), but every existing helper —
`transmission()`, `insertion_loss_db()`, `splitting_ratio()`,
`check_passive()`, `total_power()` — reads only `|coeff|^2`. Nothing in
`beef` currently computes a phase difference between two `SMatrix` entries.

That is exactly what this device's headline metrics need. From the paper's
own equations, in this codebase's notation (`i ∈ {0,2}` = driven inputs,
`j ∈ {0,1,2,3}` = the four outputs, `S = smatrix.coeff(port=f"out{j}",
in_mode=..., out_mode=...)` per whatever port-naming convention Step 6
settles on):

- **Imbalance** (paper Eq. 10): `IB_i(λ) = −10·log10[min_j|S_ji|² /
  max_j|S_ji|²]` — power-only, `transmission()` already gives this per port;
  only the min/max-over-outputs reduction is new.
- **Excess loss** (Eq. 11): `XL_i(λ) = −10·log10[Σ_j|S_ji|² / P_in(λ)]` —
  power-only, a sum of `transmission()` calls.
- **CMRR_I / CMRR_Q** (Eqs. 13–14): `−20·log10[||S_5i|²−|S_8i|²| /
  (|S_5i|²+|S_8i|²)]` and the 6/7 analogue — power-only, same shape as
  `splitting_ratio()` but as a *difference* over a *sum*, not a ratio; a
  small new reduction, not a new physical quantity.
- **Phase error** (Eq. 15): `Δφ_j(λ) = ∠S_j3(λ) − ∠S_j1(λ) − (φ̃_j − φ̃_5)`
  with target angles `φ̃_5=−45°, φ̃_6=225°, φ̃_7=45°, φ̃_8=135°` — **this one
  is genuinely new**: it needs `np.angle()` of two *different-excitation*
  `SMatrix` entries, on a shared, physically meaningful phase reference.

**Open risk, worth its own check (Step 8a) before trusting any phase-error
number**: A5's mirror-symmetry reuse builds `M_in_1` and `M_in_3` from two
*independent* `solve_section()` calls (different `input_idx`) sharing one
body solve. Every prior device in this repo only ever needed the *magnitude*
of its transfer matrices to agree across separate excitations
(`check_passive()`, splitting ratios) — nothing has ever checked that two
independently-run FDTD eigenmode decompositions preserve a consistent
**absolute phase reference** relative to each other. If they do not (e.g. an
arbitrary per-call phase convention in `get_eigenmode_coefficients()`),
`Δφ_j` would be corrupted while every power-based check
(`check_passive()`, imbalance, excess loss, CMRR) still looks perfectly
healthy — exactly the "each individual piece looks fine" failure mode the
four checks exist to catch, just for a quantity none of them currently
touch. **This needs a dedicated, direct check**, e.g. verifying `Δφ_j` from
a full single-domain FDTD reference (check 4) reproduces the hybrid's own
`Δφ_j` before trusting the hybrid number at all — a stronger requirement
than check 4 has ever had to meet on a magnitude-only device.

**Decision needed, not made here:** whether these four helpers become
script-local functions in the example (mirroring
`beef.solvers.fdtd.dbeta_from_phases()`'s precedent of a small, reusable,
device-independent helper) or a proper addition to `SMatrix`
(`SMatrix.phase()`, `SMatrix.cmrr()`, …). Given `beef.mode_identity`'s and
`dbeta_from_phases`' own precedent of factoring out device-independent
numeric helpers rather than leaving them buried in example scripts, and
given CMRR/phase-error are generic coherent-receiver metrics that could
plausibly matter to a future device, **leaning toward small, tested
additions to `beef/results/smatrix.py`** — but flagged here for a decision
before Step 7, not decided unilaterally.

> **✅ DECIDED AND LANDED 2026-08-05 — small, tested additions to
> `SMatrix`.** `beef/results/smatrix.py` gained four methods:
>
> - `imbalance_db(ports, in_mode=1, out_mode=1)` — Eq. 10, generalized from a
>   fixed pair to any number of output ports.
> - `excess_loss_db(ports, in_mode=1, out_mode=1, p_in=1.0)` — Eq. 11.
>   **Positive = loss**, the paper's own convention and the *opposite* sign
>   of this class's existing `insertion_loss_db()` — documented in the
>   method's own docstring so the two are never mixed in one report.
> - `cmrr_db(port_a, port_b, in_mode=1, out_mode=1)` — Eqs. 13/14's formula,
>   generic to any two-port balanced pair (the caller picks which two ports
>   make it CMRR_I vs. CMRR_Q).
> - `cross_excitation_phase_diff_deg(other, port, in_mode=1, out_mode=1)` —
>   the genuinely new primitive. Takes a *second* `SMatrix` explicitly
>   (`self`/`other` built from different input-port excitations of the same
>   device) rather than being a single-instance method, because **every
>   existing SMatrix-producing example in this repo holds exactly one input
>   excitation's results** (confirmed by reading `mmi_2x2_hybrid.py`'s own
>   single `smat = SMatrix(...)` for its one driven arm) — there is no
>   existing convention for "one SMatrix, two input ports" to extend, so this
>   had to be designed, not copied. Reproduces the paper's own
>   `phi_j(λ) = ∠S_j3(λ) − ∠S_j1(λ)` exactly when called as
>   `smat_in1.cross_excitation_phase_diff_deg(smat_in3, port=...)`.
>   **Deliberately stops short of the paper's full Eq. 15** — the reference-
>   port subtraction and hardcoded target angles (`φ̃_5=−45°`, etc.) are
>   device-specific, and building them from this primitive is Step 6/7's
>   example-script job, not a library concern. Output is wrapped to
>   `(-180, 180]`, with an explicit docstring warning that a target stated in
>   `[0, 360)` (`φ̃_6=225°`, `φ̃_8=135°`) must be compared modulo 360, not by
>   direct subtraction.
>
> **13 new tests** in `tests/test_smatrix.py::TestCoherentReceiverMetrics`,
> synthetic data only, hand-computable (e.g. `cmrr_db` on a 0.6/0.4 W split
> asserted against `-20·log10(0.2) = 13.9794` dBe by hand) — including a
> wraparound case (`170°` vs. `-170°` excitations correctly returning `+20°`,
> not `-340°`) and an in-range replica of two of the paper's own Eq. 15
> targets (`φ̃_5=-45°`, `φ̃_7=45°`). Full fast suite re-run clean: **639
> passed, 19 skipped** (up from 626 passed pre-change — 13 new, zero
> regressions).
>
> **What this does NOT close**, restated because it is easy to misread "the
> metrics are implemented" as "the metrics are trustworthy": the *physical*
> risk this finding opened with — whether two independently-run
> `solve_section()` FDTD excitations preserve a consistent absolute phase
> reference relative to each other — is untouched by this. These four
> methods are correct arithmetic over whatever `SMatrix` coefficients they
> are given; they cannot tell a real, physically-meaningful phase from an
> arbitrary per-call convention. That check still needs Step 6's real solves
> and is still Step 8a's job.

### A9. No published ground truth for THIS device — ARGUED, by construction

The paper's own Table 2/Fig. 3 numbers (CMRR ≥ 28.9 dBe, phase error ≤
2.52°, excess loss ≤ 0.68 dB) describe a **weakly-guided silica, C-band**
device with dimensions in the hundreds of µm (`L = 459 µm`). This device
will be **strongly-guided silicon, O-band**, with dimensions expected to be
an order of magnitude smaller (shorter `L_π` from much larger index
contrast). **Those numbers are a design-methodology reference (the
equations, the port topology, the ports-2-and-4-unused convention, the
system-requirement thresholds it itself cites), not a numerical target** —
the same caution the modal multiplexer doc already states about its own
"Seed vs. Optimized" reference figure (`modal_multiplexer_validation.md`,
"Reference design" section): don't quote a paper's numbers as this device's
ground truth when the underlying physics differs. What DOES carry over
directly and should be checked against: the **typical system requirements**
the paper itself cites from refs [4,5] — **CMRR ≥ 20 dBe, |Δφ_j| ≤ 5°** —
since those are receiver-system specs, not properties of the silica design.

---

## Part B — the ordered steps

Each step names its **deliverable** and its **evidence** (what makes it
done), same convention as `mmi_1x4_sin300_validation.md`.

### Step 0 — settle the band explicitly

**Deliverable.** `WL_MIN = 1.26`, `WL_MAX = 1.36 µm` (O-band), matching the
existing O-band devices in this repo (1×2 MMI, modal multiplexer) rather
than deriving new bounds — stated explicitly per this repo's own standard
("state the choice … do not let it be implicit," `mmi_1x4_sin300_validation.md`
step 4c). Centre wavelength 1.31 µm for any single-λ mode solve (`β0`, `β1`
for `L_π`, per paper Eq. 9).

**Evidence.** This line, plus the constant in whatever script Step 6
produces.

### Step 1 — cross-section design (cheap mode solves, no FDTD)

**1a. Single-mode access-waveguide width, at O-band, on `SOI_220`.** Sweep
strip width and count guided modes (same convention as every mode-count
table in this repo: `Re(n_eff) > n_clad`) at λ=1.31 µm. **Do not inherit
`waveguide_width=0.45`** from `mmi_2x2_hybrid.py` (A2) — that number was
never independently validated even at its own C-band, and O-band's shorter
wavelength confines more tightly. Deliverable: `waveguide_width`. Evidence:
a mode-count-vs-width table with a stated single-mode window.

**1b. Vertical/lateral box convergence (`sz`, `sim_height`).** `SOI_220`
ships `sz=2.0`, inherited from C-band devices. Sweep `sz` at the chosen
width, O-band, and confirm `|Δn_eff|` convergence — likely converges fast
given Si/SiO2's strong confinement (contrast with the SiN saga's slow
convergence at its own near-cutoff width), but state the number measured,
don't assume `sz=2.0` transfers unchanged.

**1c. `L_π` and the self-imaging length (paper Eqs. 7, 9).** Solve
`β0`, `β1` (fundamental and first-order MMI-body modes) at a candidate
`mm_width`, λ=1.31 µm: `L_π = π/(β0−β1)`. This number is expected to be far
shorter than the paper's own 459 µm silica body, by roughly the ratio of
index contrasts — compute it, don't guess the scale from the paper's number.

### Step 2 — starting dimensions from the paper's own equations

**Deliverable.** A first-pass `(W, ΔL, w, y_1..y_4)` set from Eqs. 7–9,
Table 1's parameter *roles* (not its silica-tuned *values*): `mm_width (W)`,
`ΔL` (length deviation), `waveguide_width (w)`, per-port `Δy_i` (kept at 0
initially — a reasonable start given A5's mirror symmetry already handles
the topology-level placement `y_i = 2i-5)/8 · W`, per Eq. 8, with `Δy_i` as
a later fine-tuning knob only if Step 3's sweep needs it). `L = (3/4)L_π +
ΔL` (Eq. 7) with `ΔL` initially 0, refined in Step 3.

**Evidence.** The computed starting values, with `L_π` from Step 1c, not
copied from Table 2.

### Step 3 — EME width/length refinement (replaces the paper's μGA — see A3)

**Deliverable.** A coarse grid or 1D sweep over `mm_width`/`ΔL` (and, only
if needed, `w`), using `EMESolver.solve_section_length_sweep()` on the
uniform body slab — the same cheap length-sweep machinery
`mmi_2x2_hybrid.py`/`mmi_1x4_hybrid.py` already use for their own body
refinement (a length sweep is one mode solve plus analytic phase
propagation per length, not a re-solve per point). Target: equal splitting
among the four outputs for each driven input, and the four relative phases
(Eq. 15's `φ̃_j`) within the cited system requirement (`|Δφ_j| ≤ 5°`) across
the O-band.

**Evidence.** A splitting-ratio-vs-length and phase-error-vs-length plot,
with the chosen `(mm_width, ΔL)` stated and justified (same "state the
tradeoff, don't silently pick the argmin" discipline
`mmi_1x4_sin300_validation.md`'s own `MM_LENGTH=206.0` decision used — that
device explicitly chose a *non*-optimal point for footprint/robustness
reasons and said so). **If this sweep cannot reach the cited system
requirement**, that is the trigger for A3's fallback (a real optimizer),
not a reason to loosen the requirement quietly.

**Carry forward `mmi_1x4_hybrid.py`'s grid-snap lesson explicitly:** that
script's own comment records a **real, previously-shipped bug** — a
"pinned" length argmin'd onto a sweep grid that didn't actually contain it,
0.0208 µm off for weeks. Whatever length this step picks, assert it is an
exact point on the sweep's own grid before shipping it in a script, the way
`mmi_1x4_hybrid.py:_nearest` now does.

### Step 4 — geometry

**Deliverable.**

```python
mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=...,        # Step 2/3
    mm_width=...,         # Step 2/3
    waveguide_width=...,  # Step 1a
    taper_length=...,     # TBD — measure clean-transition length, not a default copy
    taper_width=...,      # TBD
    waveguide_offset=...,
    waveguide_height=0.22,  # SOI_220
    use_neighbor_aperture=True,   # A6 — needed on both 4-port ends
)
```

Built via the **plain** constructor (no custom `inputs=`/`outputs=` lists) —
required for A5's mirror-symmetry claim to actually hold, not just be
plausible.

**Evidence.** A geometry plot (`plot_hybrid_geometry` or a bare `polygons()`
render) confirming: four input tapers, four output tapers, no taper overlap
at the chosen `waveguide_offset` (the real constraint A1 in the SiN plan's
own 1×4 analogue flags — "the four output tapers at waveguide_offset spacing
must not overlap, which is a real constraint at 1×4 that never binds at
1×2" — doubly true here with tapers on *both* ends), and that the input and
output faces are visually identical (A5's symmetry, checked, not assumed).

### Step 5 — interface mode audit

**Deliverable.** A real mode solve at the settled input/output-taper ↔ body
interface (both ends, confirmed identical per A5 rather than assumed and
audited only once), reporting guided-mode count and the FDTD/EME
correspondence (n_eff nearest-neighbor, then field-overlap η — the two-step
process `modal_multiplexer_validation.md` gap 1 established as necessary,
not the n_eff-only check that gap's own Step 1 showed can be misleading).
`N_MODES` is whatever this measures, not assumed from A7's cross-device
table.

**Evidence.** The audit table (band, n_eff both engines, Δ, guided/not),
plus explicit confirmation of whether either duplicate-band guard fires —
if it does, that is itself the finding (per `mmi_2x2_validation.md` gap 6's
own framing).

### Step 6 — three example scripts

Following every other MMI device's naming: `examples/mmi_4x4_90hybrid/python/`
— `mmi_4x4_90hybrid_eme.py`, `mmi_4x4_90hybrid_fdtd.py`,
`mmi_4x4_90hybrid_hybrid.py`.

**6a. EME script.** Full-device EME solve across the O-band, both driven
inputs — the cheapest of the three, and the first place A3's design (does
the self-imaging + sweep result actually hit the cited CMRR/phase-error
requirements) gets checked end-to-end rather than on a bare uniform-body
cut.

**6b. FDTD script (full single-domain reference, for check 4).** Flag the
cost concern up front, same as the modal multiplexer's own doc: a 4-port,
double-tapered device is a long box, and 3-D FDTD cost scales steeply with
domain size. The modal multiplexer's own `CompoundGeometry.polygons()`
stub-per-port fix (`notes.md` Geometry item 8 part 2, landed 2026-08-04 —
"`_end_ports()` already answers exactly this 'which port's width wins'
question … by *not* picking a single winner: it returns every port at that
end") is **directly relevant and already landed** — this device needs the
exact same multi-port stub behaviour at both ends that fix enables, so this
script is not blocked on it the way the mux's once was. Budget a smoke test
at reduced resolution first (this repo's own "two-step discipline for
expensive, never-run scripts" — every prior full-device FDTD script in this
directory did this before its real run).

**6c. Hybrid script.** Input/body/output split per A5/A6: `INPUT_METHOD =
OUTPUT_METHOD = "fdtd"`, `BODY_METHOD = "eme"` (this repo's standard split
for every MMI hybrid so far), `use_neighbor_aperture=True` on the geometry
(not `half_domain_aperture` on the solvers — A6), 4 FDTD input-arm solves +
1 EME body length-sweep, output section declared but never separately
simulated (mirroring `mmi_2x2_hybrid.py`'s own `sub_region_configs`
convention for documenting *which* sub-region plots correspond to a real
simulation vs. a mirrored reuse).

### Step 7 — metrics helper (A8) — ✅ math primitives landed 2026-08-05; the plot/report itself still pending Step 6

**Deliverable.** Imbalance, excess loss, CMRR_I/Q, and phase-error
functions per A8's equations, plus the phase-consistency check A8 flags as
a new risk class (verify the hybrid script's cross-excitation `Δφ_j`
against the full-device FDTD reference from Step 6b, not just against
itself). Location (script-local vs. `SMatrix` method) per A8's open
decision.

**✅ Landed:** `SMatrix.imbalance_db`/`excess_loss_db`/`cmrr_db`/
`cross_excitation_phase_diff_deg` — see A8's own boxed update for the full
account. **Still pending, blocked on Step 6:** there is no real device yet
to plot, and the phase-consistency risk these primitives cannot answer on
their own.

**Evidence.** A CMRR-vs-λ / phase-error-vs-λ / excess-loss-vs-λ plot in the
shape of the paper's own Fig. 3, annotated against the **cited system
requirements** (CMRR ≥ 20 dBe, |Δφ_j| ≤ 5°) — not against the paper's own
silica numbers (A9). Not yet produced — needs Step 6's real S-matrices.

### Step 8 — the four checks

Once Step 6's three scripts exist and run for the first time, this doc
converts to a normal gap list against
[README.md's four checks](README.md#the-four-checks):

- **Check 1 (engine agreement).** Uniform-body-stretch phase fit
  (`verify_item_6_engine_agreement.py`-style), on this device's own
  cross-section — cannot borrow any existing device's ratio (every prior
  measurement is at a different (platform, band) pair).
- **Check 2 (mode identity).** Elevated risk here per A7 — four
  closely-spaced ports per end is closer to the modal multiplexer's
  16-mode regime (where check 2 first failed, on a real bug) than to either
  existing MMI's two-or-one-port ends.
- **Check 3 (resolution convergence, per engine).** Sweep independently at
  the new (SOI_220, O-band) cross-section — the 2×2's res 100 (C-band) does
  not transfer, same logic as every prior "new pair" in this repo.
- **Check 4 (full-device FDTD reference).** Per A8, this device's check 4
  carries a strictly stronger requirement than any prior device's: it must
  reproduce `Δφ_j`, not just total transmitted power.

### Step 9 — golden fixture

Add to `tests/test_golden_examples.py`'s `TARGETS` **after** Step 8's
resolution question settles — every other device's doc repeats this
ordering for the same reason (recording a golden against a configuration
still expected to move means regenerating it twice).

---

## Open decisions carried into Part B, not made unilaterally here

1. **A3** — self-imaging + EME sweep vs. porting a real optimizer, if the
   sweep can't hit the cited requirements.
2. **A8** — script-local metric helpers vs. new `SMatrix` methods.
3. **Exact taper length/width, `waveguide_offset`, and whether any per-port
   `Δy_i` fine-tuning (paper Eq. 8's own knob) turns out to be needed** —
   all explicitly TBD in Step 2/3, not guessed here.
4. **Directory/file naming** — used `mmi_4x4_90hybrid` above for
   consistency with `mmi_2x2`/`mmi_1x4`'s own `mmi_<NxM>` convention plus a
   disambiguating suffix for the terminology collision flagged at the top;
   open to a different name before Step 6 creates the directory.

---

## References

Primary reference for this device's topology, equations, and terminology:

> T. Yasui, J. Sugisaka, and K. Hirayama, "Structural Optimization of an
> Optical 90 Degree Hybrid Based on a Weakly-guided 4×4 Multimode
> Interference Coupler Using a Parallelized Real-coded Micro-genetic
> Algorithm," *ACES Journal*, vol. 36, no. 5, pp. 526–532, May 2021.
> https://doi.org/10.47037/2020.ACES.J.360506

Two papers it cites that matter directly to this plan's decisions (A3, A9):

> R. Halir, G. Roelkens, A. Ortega-Moñux, and J. G. Wangüemert-Pérez,
> "High-performance 90º hybrid based on a silicon-on-insulator multimode
> interference coupler," *Opt. Lett.*, vol. 36, no. 2, pp. 178–180, Jan.
> 2011. — the SOI-based precedent motivating A3's self-imaging-without-GA
> approach, and one of the two sources (with [5] below) for the "typical
> system requirements" (CMRR ≥ 20 dBe, |Δφ| ≤ 5°) A9 checks against instead
> of the paper's own silica numbers.

> L. B. Soldano and E. C. M. Pennings, "Optical multi-mode interference
> devices based on self-imaging: Principles and applications," *J. Lightw.
> Technol.*, vol. 13, no. 4, pp. 615–627, Apr. 1995. — the general
> self-imaging theory Eqs. 7–9 are drawn from, referenced in A3 for why it
> is expected to hold well on a strongly-guided platform even though it
> fails on the paper's own weakly-guided one.

> J. S. Fandino and P. Munoz, "Manufacturing tolerance analysis of an
> MMI-based 90º optical hybrid for InP integrated coherent receivers,"
> *IEEE Photon. J.*, vol. 5, no. 2, pp. 7900512, Apr. 2013. — the second
> source (with Halir above) for the CMRR/phase-error system requirements.
