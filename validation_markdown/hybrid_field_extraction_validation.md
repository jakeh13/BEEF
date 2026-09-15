# Hybrid-path field extraction — BUILD PLAN

Opened 2026-08-31. Cross-device capability, not a per-device gap list — the
device that motivated this is the 4×4 90° hybrid's BO-optimized best point
(`mmi_4x4_90hybrid_wg_eval.py`, W=6.50/y_gap=0.51/mm_length=91.00,
objective=0.6767 as of the campaign's round 5), which has no field/mode-
profile picture at all today, only scalar S-matrix metrics. This doc plans
that capability, but **starts on the 1×2 MMI** (`mmi_1x2_hybrid.py`) — the
cheapest device in the repo that goes through the same `HybridSolver`/
`FDTDSolver.solve_section()` + `EMESolver` cascade construction the 4×4
device does, so the new code gets built and proven somewhere a wrong attempt
costs minutes, not a multi-hour 4×4 rerun.

## Why nothing here is free, checked against the real code

**No field data exists anywhere in the hybrid-path construction today** —
confirmed by reading `beef/solvers/fdtd.py`'s `solve_section()`
(~line 2673) and its excite loop (~line 3260-3380) directly: every
`mp.Simulation` it builds only ever adds flux/mode-decomposition monitors
(`sim.add_flux`, `sim.get_eigenmode_coefficients`) before `sim.run()`. There
is not one `add_dft_fields()`, `get_array()`, or `Animate2D` call anywhere
in the file (grepped, zero hits). This isn't a case of spatial field data
being computed and then discarded the way the excite-normalization fix
found 93% of an MPB decomposition being thrown away — the field data is
**never computed in the first place**. Any field plot from the FDTD side of
this cascade is genuinely new work, not exposing a hidden existing result.

**The one field-plot capability that DOES exist (`plot_eme_fields()`,
`beef/plot/fields.py`) is pure-EME and needs live mode objects, not the
cached arrays.** Its call sites pass it meow's own `cells`/`modes` objects —
each mode object carries its full per-cross-section field grid, computed
fresh by that call. `HybridSolver`'s on-disk cache (`save_dir`) persists
n_eff/transfer-matrix **arrays** (`beef/cache/`), never the meow `Mode`
objects themselves. So even the "cheap" half of this plan — the EME body,
whose length-sweep is free by `solve_section_length_sweep()`'s own
analytic-reconstruction tier-3 cache — still needs one fresh (if modest)
EME mode solve to get real field data back. Nothing here is a pure cache
read.

> **CORRECTION 2026-08-31 — this doc originally said `plot_eme_fields()` had
> "both its existing call sites (`plot_mmi1x4_eme_fields.py`,
> `plot_mmi4x4_90hybrid_eme_fields.py`)". That undercounted.** Grepped
> `plot_eme_fields` across `--include=*.py`: **seven** example call sites, not
> two — the two named plus `plot_mmi1x4_waist_profile_eme_fields.py`,
> `plot_mmi1x4_straight_vs_tapered_eme_fields.py`,
> `analyze_mmi1x4_self_imaging.py`, `plot_mmi2x2_eme_fields.py`, and — the one
> that matters for the plan below — **`plot_mmi1x2_eme_fields.py`**. The
> helper is far more established than this doc assumed, which is good news for
> Steps 1/3 (it is proven across four device families, not two) and it
> **invalidates Step 1's stated premise** — see Step 1.

**The pure-EME method itself is already known to diverge from real FDTD on
this device family.** The 4×4 campaign's own item 1 measured pure-EME
missing port crosstalk by 44-62° against full-device FDTD, wrong candidate
ranking. A field picture built the same way `plot_mmi4x4_90hybrid_eme_fields.py`
does (whole-chain EME, `solve()`) would visualize the DISQUALIFIED method,
not the trusted hybrid-path the BO objective is actually built on. Worth a
qualitative look, not a substitute for the real deliverable (Steps 3-4).

## Step 0 — re-baseline 1×2 before building on it (open, not yet done)

`resolved/mmi_1x2_validation.md` closed 2026-08-01 (file mtime), all four
checks passed at that point. Real shared-library changes have landed since,
and this device's cache is **not** provably safe against them without
checking:

- **`FDTDSolver.solve_section()`'s excite-normalization fix (2026-08-09,
  see [[fdtd-solve-section-near-degenerate-normalization-fix]])** changes
  `solve_section()`'s output whenever two excited modes are near-degenerate
  on an x-axis excite branch — and per that fix's own note, **the cache is
  NOT auto-invalidated** (no code-version component in the fingerprint), so
  a pre-2026-08-09 cached result keeps silently serving the old (possibly
  wrong) numbers. Whether 1×2's own excited-mode spacing is close enough to
  trigger the fallback is **unchecked** — first thing to determine, from
  the same MPB seed log this device's own gap-3 mode-identity audit already
  produces (n_eff gaps between the `N_MODES=8` carried bands).
- **`mmi_1x2_hybrid.py` itself was touched 2026-08-17** (file mtime) —
  16 days after the validation doc's last recorded numbers (2026-08-01) and
  matching the date range of the calibration-figure bug fix noted in
  `resolved/section_mode_beta_multi_monitor_validation.md` step 10
  (`.last_calibrations` wiped by a second `solve()` call sharing an
  in-memory cache — that step's own text says mmi_1x2_hybrid was fixed,
  mmi_1x4_hybrid "not yet re-run"). Whether this device's own S-matrix
  numbers moved as a result, or only its calibration-figure output did, is
  unchecked.
- **Not expected to matter, confirmed by the fix's own written scope**: the
  2026-08-09 output-monitor-width bug (the modal-multiplexer one) explicitly
  excluded 1×2 — grepped at the time to confirm it already passes its own
  explicit `output_monitor_widths` and never falls onto the buggy default.
  No re-check needed here specifically for that half of the fix.

**Concrete step 0 actions**: (a) check 1×2's MPB seed log for near-degenerate
excited-mode pairs; (b) if any exist (or regardless, as a cheap regression
check since the device is small), clear that section's cache entries and
re-run `mmi_1x2_hybrid.py` under current code; (c) diff the resulting
S-matrix numbers against `resolved/mmi_1x2_validation.md`'s own final
closed table. Only once this diff is explained (byte-identical, or a real
understood delta) is 1×2 a trustworthy base to build field-extraction code
against — building on a stale, silently-wrong baseline would make any new
field plot unverifiable against anything real.

### Step 0 progress — 2026-08-31

**(1) Warm re-run: bitwise identical, and that settles nothing.** Job
**12627197** (`batch/rebaseline_mmi_1x2.batch`, beef `8e1ba48`) re-ran
`mmi_1x2_hybrid.py` at the shipped config (res 40, `N_MODES=8`,
`calibrate_body_neff=True`). It finished in **43 s** (29.4 s run time) — i.e.
**fully cache-warm** — and every transmission value is **bitwise identical**
to the preserved 2026-08-18 output, both ports, all 11 wavelengths
(`np.array_equal` True; preserved as
`example_output/mmi_1x2/mmi_1x2_hybrid_PRE_STEP0_20260818/`):

| λ (µm) | 1.260 | 1.310 | 1.360 |
|---|---|---|---|
| `out_top` | 0.431112 | 0.374899 | 0.312504 |
| `out_bot` | 0.431101 | 0.374885 | 0.312493 |

> ⚠ **A cache-warm run cannot answer this step's question.** Bitwise identity
> here proves only that the cache read back consistently — it re-serves the
> very entries whose validity is in doubt. The 2026-08-09 fix has no
> code-version component in the fingerprint, so a pre-fix entry is *expected*
> to reproduce itself exactly. **This is the trap this step exists to avoid,
> and the warm run walked into it.** Recorded rather than quietly re-run so
> the distinction stays visible.

**(2) The Δn_eff proxy says the fallback should not fire — 4× margin.** The
fix's gate is `_EXCITE_NORM_MISMATCH_TOL = 0.1` (a 10% `sqrt(get_fluxes())`
vs `|alpha|` disagreement); the near-degeneracy that fired it on the modal
multiplexer's connector section was **Δn_eff = 0.0148**. This device's
consecutive excited-mode gaps at the body cross-section (`N_MODES=8`, res 40,
1.31 µm, from gap 3's re-audited table, job 11608940) are 0.0592 / 0.1014 /
0.1484 / 0.2029 / 0.2684 / 0.1411 / 0.0696 — **smallest 0.0592, 4.0× the
gap that fired the fallback.** So no fire is the *expected* result.

> **But Δn_eff is only a proxy for the real gate**, which is the 10%
> flux-vs-alpha mismatch, not a mode spacing. Per
> [[verification-standard-measured-not-argued]] this margin is a prediction,
> not the measurement — hence (3).

**(3) Cold re-run in flight — job 12627227**
(`batch/rebaseline_mmi_1x2_cold.batch`, 24 cpus, 5 h wall). Forces every
section cold by pointing `save_dir` at a **fresh** cache via a new
`MMI1X2_SAVE_DIR` env override, leaving the shared `example_sim_cache`
untouched — a fresh-vs-fresh comparison rather than a destructive cache clear
([[concurrent-sessions-shared-library-hazard]]'s own guidance). Precedent for
cost: the last cold calibrated 1×2 hybrid was **1:28:14** (job 11608690).

**Discriminator, decided before the job lands** — grep the log for
`solve_section: excite band N normalization mismatch` (`fdtd.py:3326`):
- **no warning + numbers match the warm run bitwise** → the cached entries
  were valid, the 4× margin held, and 1×2 is a trustworthy base. Step 0
  closes.
- **warning fires, or numbers move** → the cache *was* stale, every number in
  `resolved/mmi_1x2_validation.md` gaps 4/5 needs re-measuring, and Step 4's
  ground truth moves with it.

### ❌ Step 0 cold run TIMED OUT — 2026-08-31, job 12627227

**The cold re-baseline did not finish. It hit its 5:00:00 wall (`TIMEOUT`) and
Step 0 is NOT closed.** Recorded as a failure rather than quietly re-run,
because both the cost and the partial evidence matter.

| | |
|---|---|
| wall | 5:00:00, killed (precedent for this run was **1:28:14**, job 11608690) |
| FDTD domains built | 8 |
| FDTD runs actually **completed** | **6** |
| normalization-mismatch warnings | **0** |
| cache entries persisted | **2**, both written in the first 40 min (11:08, 11:11) |
| S-matrix produced | none |

**Two things worth separating.**

**(a) The partial evidence is real but incomplete.** The mismatch check runs
*after* `sim.run()` returns, so exactly the **6 completed** runs exercised it —
not all 8 built domains — and all 6 were clean. A full cold 1×2 needs ~11 FDTD
runs (3 section solves + ~8 calibration modes), so this is roughly half
coverage. **It is not "8 clean domains"**, and partial coverage does not
cleanly cover the risky cases: this device's tightest mode gaps sit at *both*
ends of the basis (bands 0-1 at 0.0592, bands 6-7 at 0.0696), so an untested
subset could still contain the near-degenerate pair.

**(b) The cache does NOT make a resubmit cheap.** Only 2 entries persisted, both
from the first 40 minutes; domains 3-8 ran 3h20m and wrote nothing, because the
calibration path caches on completion. So ~3.5 h of the 5 h is unrecoverable,
and a plain resubmit would redo it.

**Why this happened, and the fix.** `decay_tol=1e-11` — inherited here with no
recorded justification — produced 8,800-20,000-step domains at ~0.136 s/step,
~30-40 min each. That is the cost this device's tolerance was silently buying.
See the `DECAY_TOL` block now in `mmi_1x2_hybrid.py`.

**Replaced by job 12636598, cold at `decay_tol=1e-9`** (8 h limit, its own cache
dir and its own `MMI1X2_OUT_SUFFIX=_tol1e9` output dir so it cannot clobber the
1e-11 result). Note what this changes about the question being asked:

> Step 0 originally asked *"are the existing 1e-11 cache entries stale under
> current code?"* With the shipped default now 1e-9, **that question is
> becoming moot rather than answered** — nothing going forward reads those
> entries. A cold 1e-9 run rebuilt under current code is trustworthy **by
> construction**, with no staleness question to settle. That is a legitimate
> way to close Step 0's actual purpose ("is 1×2 a trustworthy base?"), but it
> is NOT the same as answering the original question, and it moves work: the
> validation doc's gap 4/5 numbers are all at 1e-11, so the new baseline needs
> its own gap 4/5 measurement before Step 4 can validate against anything.

⚠ **Still owed either way**: the T-vs-`decay_tol` sensitivity on this device is
UNMEASURED, so 1e-9's accuracy cost here is unquantified. The instrument is
`FDTDSolver.solve(snapshot_every=, snapshot_path=)`, which yields the whole
transmission-vs-tolerance curve from ONE run.

*Library change made for this step*: `mmi_1x2_hybrid.py`'s `save_dir` now
honours `MMI1X2_SAVE_DIR` (mirroring the file's existing
`MMI1X2_CALIBRATE_BODY` override), so a forced-cold verification of this
device is repeatable without touching the shared cache. One line; default
behaviour unchanged.

## Step 1 — retarget the existing 1×2 EME field script at the hybrid's own body

> **REWRITTEN 2026-08-31. The original text — "Adapt `plot_eme_fields()`
> (already exists, already proven on two other devices) to 1×2's own body
> cross-section — no companion script exists for this device yet" — was
> factually wrong on its central premise, and the correction shrinks this
> step substantially.**

**`examples/verification/python/plot_mmi1x2_eme_fields.py` exists, has
existed since 2026-08-08 (commit `76b7465`, 23 days before this doc was
opened), and has already been RUN** — all six field components are on disk
at `example_output/mmi_1x2/mmi_1x2_eme/fields/` (`eme_{Ex,Ey,Ez,Hx,Hy,Hz}_
out_top_1.310um.png`, mtime 2026-08-08 08:52). There is nothing to author
here and no "does `plot_eme_fields()` work on this device" risk left to
retire. A 2×2 companion (`plot_mmi2x2_eme_fields.py`) landed the same day.

**What is genuinely still missing is narrower, and it matters.** The existing
script is pinned to `mmi_1x2_eme.py`'s **pure-EME** configuration, which is
*not* the body section the hybrid cascade actually propagates:

| | existing field script | hybrid's body section |
|---|---|---|
| mode solver | `tidy3d` | **`meep`** (matches FDTDSolver's basis) |
| body modes | 12 | **8** (`N_MODES`, set by gap 3) |
| body span | 4.5 µm (full `mm_length`) | **4.0 µm** (`±(mm_length/2 − body_stub)`, `body_stub=0.25`) |
| chain | all 5 regions (wg/taper/body/taper/wg) | **body only** |

So the picture already on disk visualizes the **whole-chain pure-EME
method** — precisely the construction this doc's own preamble disqualifies
("would visualize the DISQUALIFIED method, not the trusted hybrid-path the
BO objective is actually built on"). It is the qualitative look the preamble
allows, and it is already taken.

**Revised Step 1 action**: copy `plot_mmi1x2_eme_fields.py` to a
hybrid-pinned variant that builds its single cell from
`mmi_1x2_hybrid.py`'s own body-section parameters (`mode_solver="meep"`,
`num_modes=[8]`, `num_cells=[1]`, span ±2.0 µm) instead of the pure-EME
five-region chain, and re-run. Still one modest mode solve; still de-risks
the "fresh mode-object solve alongside cascade-cached coefficients"
mechanics Steps 3-4 depend on. But it is a ~20-line retarget of a proven,
already-exercised script — not new capability — and it should be sequenced
as such rather than as a prototype milestone.

### ✅ Step 1 DELIVERED — 2026-08-31

**The body field plots standalone.** No FDTD field capture, no library change,
no cascade cache — one mode solve on one cross-section, **55 s measured**, run
on the login node. Two scripts, both in `examples/verification/python/`:

- **`plot_mmi1x2_hybrid_body_fields.py`** — the direct retarget. Reproduces
  `EMESolver.solve_section()`'s own cell construction for the body section and
  hands it to the existing `plot_eme_fields()`. Output:
  `example_output/mmi_1x2/mmi_1x2_hybrid_body/fields/`.
- **`plot_mmi1x2_hybrid_body_selfimaging.py`** — the useful one (see the two
  defects below that made it necessary). Device-frame, multimode, one figure.
  Output: `.../mmi_1x2_hybrid_body/hybrid_body_selfimaging_1.310um.png`.

**Three independent self-checks passed, none of them tuned for:**

| check | result | why it is a real check |
|---|---|---|
| mode basis identity | all 8 n_eff match gap 3's meow column to **6 dp** | that table was measured separately, job 11608940 |
| parity | odd modes at **1e-19 … 1e-11** | a centred input cannot excite an odd body mode |
| self-image vs ports | exit lobes peak at **y = ±0.50 µm** | `MMI._port_ys()` puts the outputs at exactly ±0.50 |

The third is worth keeping: `waveguide_offset=1.0` is the port centre-to-centre
**spacing**, so `output_ys = [-0.5, +0.5]` — *not* ±1.0. (`mmi_1x2_hybrid.py`'s
own section-5 comment says "output_ys = [-wo, +wo]", which is wrong by a factor
of two; harmless there because only the index is used, but it misled this work
once already.) The reconstructed field lands the split exactly on the ports,
which validates body length against port placement from the field rather than
from design intent.

**Two defects in `plot_eme_fields()` this exposed, both affecting all seven
existing call sites' interpretation:**

1. **It labels in meow's frame, not the device's.** It passes its component
   string straight to `mw.propagate_modes()` with no remap, but
   `stack_meow_mode_fields()`'s docstring records the axes are rotated one
   step (device Ex = meow Ez, device Ey = meow Ex, device Ez = meow Ey). So
   the panel titled `|Ex|` is really the device's **in-plane transverse**
   field (dominant, peak ≈ 2.0) and the one titled `|Ey|` is device Ez (weak
   parasitic, peak ≈ 0.006, ~330× smaller and edge-concentrated — correct, but
   easily misread as a broken solve). Its `x [µm]` axis is device y. **A real
   hazard for Step 3**, which must stitch these against FDTD fields that ARE
   device-frame. Not fixed in the library here — flagged, and worked around in
   the new script by going through `stack_meow_mode_fields()`.
2. **It can only excite one mode at a time** (`excite_mode_l`), and a single
   eigenmode in a uniform section is invariant along propagation by
   construction. That panel is correct and says nothing about what the body
   does; the body's whole job is the multimode beating that sets `mm_length`.

**Method note — a wrong inner product, caught by conservation.** The first
projection of the input onto the body basis used a 1-D cut through the core and
returned **total captured power 1.5857** (impossible for an orthogonal basis)
with **0.5876 on mode 7** — an odd mode. Both artifacts came from the same
error: a 1-D slice is not the inner product these modes are orthogonal under.
Switching to the 2-D transverse area integral with `dy·dz` weights (what
`field_overlap_matrix()` already uses) fixed both at once — total **0.9970**,
mode 7 at **1.0e-11**. Measured basis max off-diagonal overlap **0.069**,
non-zero only between same-parity modes. **The parity zeros survived both
inner products, so parity alone would NOT have caught this — power
conservation is what exposed it.** Worth carrying into Step 3, which needs
overlaps at both FDTD↔EME boundaries.

⚠ **The excitation is an APPROXIMATION** — a centred Gaussian of
FWHM = `taper_width`, projected onto the body basis — not the field the input
FDTD section actually delivers, which does not exist until Step 2. It is the
right qualitative stand-in (centred input → even modes → self-imaging) but no
number from this figure may be quoted as a hybrid interface result. Replacing
that one array with the real captured field is precisely Step 3.

Figures and the full write-up: `hybrid_field_extraction_step0_1.html`
(published artifact).

### ✅ Step 1 GENERALIZED to a shared helper + four devices — 2026-08-31

The 1×2 work was factored into **`beef/plot/hybrid_body.py`** (new module,
purely additive — no existing library file's behaviour changed, verified by
re-importing `plot_eme_fields` after the edit). Exports:
`plot_hybrid_body_fields()`, plus `solve_body_modes()` /
`body_mode_weights()` for callers wanting the pieces.

**The abstraction was validated by regression, not by inspection**: the 1×2
script was rewritten onto the helper and reproduces every number from the
hand-rolled version — n_eff to 6 dp, weights 0.9022/0.0935/0.0013, captured
power 0.9970, gram off-diagonal 0.0690, exit lobe y=0.488 vs port 0.500. Those
values are now **assertions inside the script**, so a future change to the
helper or the mode solver fails loudly instead of silently redrawing.

| device | script | modes / solver | result |
|---|---|---|---|
| 1×2 MMI | `plot_mmi1x2_hybrid_body_selfimaging.py` | 8 / meep | ✅ reference case, asserted |
| 1×4 MMI (SiN 300) | `plot_mmi1x4_hybrid_body_fields.py` | 12 / **tidy3d** | ✅ captured 0.9987, parity exact |
| DC SOI_220 | `plot_dc_hybrid_body_fields.py` (`DC_VARIANT=soi220`) | 5 / meep | ✅ Δβ matches doc to **0.16%** |
| DC SOI_160 | same, `DC_VARIANT=soi160` | 2 / meep | ✅ gram 0.0000, 0.498/0.462 |
| HCG-FP coupler body | `plot_hcgfp_hybrid_body_fields.py` | 2 / meep | ✅ agrees with DC160 |
| modal multiplexer | `plot_mux_hybrid_body_fields.py` | 15 / 6 | ⚠ **written, NEVER RUN** — deferred |
| 2×2 tapered, 4×4 | — | — | handed off, see below |

**Two independent cross-checks landed, which is what makes this more than
"the code runs":**

- **DC SOI_220, Δβ.** The reconstruction gives supermode Δn_eff = 0.022088 at
  res 80 → **Δβ = 0.089541 rad/µm**. `resolved/directional_coupler_validation.md`
  item 11 independently measured *"MPB plateaus at Δβ = 0.0894 by res 80–100"*.
  Same MPB backend, same resolution, **0.16% agreement** — and it is now an
  assertion in the script (fails past 1%).
- **1×4, self-image vs ports.** Exit lobes land at y = **1.583 / 4.65** against
  output ports at **1.55 / 4.65** — the outer lobe exact, the inner 0.033 µm
  outside its port. This independently confirms the post-2026-08-14 geometry
  (`MM_LENGTH=39.90`, `WAVEGUIDE_OFFSET=3.1`) does what it was changed to do.
  ⚠ **Possible lead on that device's open step-6d question** (EME favours inner
  ports 1/2, hybrid favours outer 0/3): the inner lobe sitting slightly outboard
  of its port is a plausible mechanism for an inner-vs-outer imbalance. Reported
  as a lead, **not** asserted — the excitation here is approximate, so this
  cannot settle it.

**Coupler bodies turn out to be the most informative case.** The body modes are
the even/odd supermodes; a single-arm launch is a ~50/50 mix of them (measured
0.4894/0.4634 on the 220, 0.4981/0.4621 on the 160), so the map shows power
sloshing between arms at the beat length — the coupler's entire function, and
the same Δβ the calibration corrects. The SOI_160 body's gram off-diagonal is
**0.0000**: even/odd supermodes are exactly orthogonal by symmetry, unlike an
MMI's same-parity mode set (0.069 on the 1×2).

**Devices deliberately NOT covered, with reasons:**

- **ring resonator** — not applicable at all. It hand-rolls raw meep/meow and
  never goes through `FDTDSolver`/`HybridSolver`, so there is no cascade body
  section (its own docstring explains why: `FDTDSolver` only builds x-normal
  sources/monitors).
- **2×2 tapered/waisted MMI** — the helper **cannot** apply: a tapered span
  fails `assert_uniform_cross_section()` by design, and needs the real
  multi-cell path with a non-diagonal transfer matrix. Handoff text written.
- **4×4 90° hybrid** — should apply directly (uniform body), but needs an
  off-axis launch and, more importantly, **a phase panel**: `|E|²` loses
  precisely the relative-output-phase the device exists to produce. Handoff
  text written.

Handoff notes for those two: `hybrid_body_field_plotting_handoff.md`.

### Phase readout added — 2026-08-31

The first-listed gap in the rollout above ("|E|² loses what a 90° hybrid is
for") is closed on the library side. `plot_hybrid_body_fields()` now always
returns the **complex** field as `res["E"]`, and **`beef.port_phases()`**
(exported) reads the relative phase at each port position off it:

```python
p = port_phases(res, geometry.output_ys, reference=0)
p["phase"]      # degrees, relative to port `reference`
p["amplitude"]  # |E| at each port
```

Two guards are written into its docstring because both are real failure modes,
not theoretical ones:

- **Amplitude gates phase.** The phase of a vanishing field is numerically
  arbitrary — the same trap [[mode-signs-gate-on-strong-modes]] records, where
  a 1e-06 mode disagreed by 280%. `amplitude` is returned alongside so a phase
  on a near-null port can be recognised as meaningless rather than quoted.
- **Only RELATIVE phase is meaningful.** The launch carries an arbitrary global
  phase, so the helper subtracts a reference port rather than reporting
  absolutes.

⚠ This is the body's **own** phase accumulation from an approximate launch —
not the cascade's port-to-port phases, which still need Step 2. Useful for
checking that the body does what coupled-mode theory says; not a substitute for
a device-level phase measurement.

**✅ Validated against two predictions that do not depend on this code being
right.**

**(1) 1×2 symmetric splitter — outputs must be in phase.** Mirror symmetry
forces 0°. Measured **−0.21°**, both ports strong (|E| 2.181 / 2.192). The
residual is grid discretisation — the two port positions do not fall on exactly
mirror-symmetric cell centres.

**(2) Directional coupler — cross port must sit near quadrature.** First read
gave **+99.68°** against theory's 90°, with the cross port weak (|E| 0.504 vs
3.650). Rather than accept a 9.7° gap, it was swept: propagation is analytic, so
**one** mode solve gives the phase at any body length for free.

| L/L_π | cross phase | P_cross fraction |
|---|---|---|
| 0.086 | +99.70° | 0.019 |
| 0.150 | +95.72° | 0.055 |
| 0.250 | +93.66° | 0.147 |
| 0.350 | +92.88° | 0.273 |
| **0.500** | **+92.52°** | **0.500** |
| 0.650 | +92.78° | 0.727 |
| 0.750 | +93.47° | 0.853 |
| 0.900 | +97.83° | 0.975 |
| 1.000 | −180.00° | 0.9996 |

**The phase reaches its minimum, +92.52°, exactly at L = L_π/2** — the 50/50
point where both ports carry equal maximum power and coupled-mode theory's 90°
is cleanest. The deviation then grows **symmetrically in both directions** as
either port goes weak (99.70° at 8.6%, 97.83° at 90%). That symmetry about
L_π/2 is the signature of a **weak-field artifact, not a bias in the
machinery** — a systematic error in the phase computation would not fold neatly
about the equal-power point. At L = L_π the driven port nulls (|E| 0.078) and
the phase reads −180°: meaningless, exactly as `port_phases()`'s own docstring
warns.

The residual 2.5° from ideal at L_π/2 is the approximate Gaussian launch, whose
even/odd split is 0.4981/0.4621 — a 7.5% imbalance rather than a perfect 50/50.

**So the "amplitude gates phase" guard is not defensive boilerplate — this
sweep is the measurement showing exactly where the phase stops meaning
anything, and by how much.** Neither test device is a 90° hybrid, so the 4×4
remains the first real exercise of the helper; the handoff says so.

## Step 2 — add spatial field capture to `solve_section()`'s FDTD runs

New library code, following the exact non-invasive pattern already used for
`decay_method` (`beef/cache/config.py`'s two-part fix: default fingerprint
untouched so every existing device's cache is unaffected, a new opt-in
value gets its own distinct fingerprint). Proposed: an opt-in flag that adds
`mp.add_dft_fields()` at the section's own excitation frequency before
`sim.run()`, then `sim.get_array()` once the existing decay condition fires
— alongside, not instead of, the existing eigenmode-coefficient extraction
that produces the transfer matrix.

> **✅ Both premises verified against current code 2026-08-31.**
>
> **(a) The "no field data exists" claim holds exactly.** Re-grepped
> `beef/solvers/fdtd.py` for `add_dft_fields`, `get_array`, `Animate2D`,
> `output_efield`, `to_appended` — **zero hits each, 0 total.** Step 2 is
> genuinely new work, as claimed. `solve_section()` is at **`fdtd.py:2673`**
> and the excite-loop guard at **`fdtd.py:3324-3326`**, both matching this
> doc's original line references.
>
> **(b) The fingerprint pattern is real and better-established than cited.**
> `beef/cache/config.py` carries `_DEFAULT_OMITTED_ATTRS` (line ~163) with
> **four** precedents, not one: `minimum_run_time`, `decay_method`,
> `energy_decay_dt`, and `region_cell_lengths` — the last added **2026-08-29,
> two days before this doc was opened**, and its in-file comment states the
> rule verbatim ("a non-default value … is a genuinely different
> discretisation … must NOT be served by, or serve" a default entry). Cite
> `region_cell_lengths` as the template; it is the freshest instance.
>
> **(c) The decay-timing open question below is exactly the right question,
> and it decides the fingerprint treatment** — which the original text left
> implicit. The two-part pattern forces a choice: an attribute goes in
> `_METHOD_ATTRS` (own fingerprint) *only if* it makes the run a genuinely
> different simulation. So — if adding the DFT monitor shifts decay timing,
> field capture changes the physics run and **must** carry its own
> fingerprint; if it provably does not, it may be default-omitted and every
> existing cache entry stays valid. **Measure first, then choose. Do not
> pick the fingerprint treatment before the timing number exists.**
>
> Note also `solve_section()`'s signature returns a bare `np.ndarray` (the
> transfer matrix), so captured fields need a side-channel rather than a
> changed return type — `HybridSolver.last_calibrations` is the existing
> precedent in this repo for exactly that (a solve stashing diagnostic data
> the caller can plot without paying for a second run).

> **Groundwork done 2026-08-31 — the open question below is REAL, and the
> repo's existing comment does NOT already answer it.**
>
> `fdtd.py:3278-3281` currently asserts, about adding extra measurement
> monitors: *"sim.run() below times-steps once regardless of monitor count
> (`stop_when_dft_decayed()` watches the whole field, not any one monitor), so
> N positions cost about the same wall-clock time as 1."* Read quickly, that
> looks like it settles Step 2's cost question. **It does not, and the
> distinction matters:**
>
> - That claim is about **N additional coplanar flux monitors at the same
>   frequencies** as the existing ones. Those accumulate similar DFT content
>   and converge on a similar schedule, so the observation is unsurprising.
> - `self._stop_condition()` (`fdtd.py:1070-1077`) returns
>   `mp.stop_when_dft_decayed(tol, minimum_run_time)`, which examines **every
>   DFT object registered in the simulation** and stops only once all are below
>   tolerance. It is not literally "the whole field" — the comment's parenthetical
>   is loose.
> - An `add_dft_fields()` monitor is a **volumetric** DFT over the cell, not a
>   plane. It samples regions the flux planes never see — near the PML, in
>   corners where energy lingers — and those may converge **later** than the
>   port planes. If they do, registering it extends the run for every caller
>   that opts in.
>
> **So the existing evidence covers plane monitors at identical frequencies and
> says nothing about a volumetric one.** The measurement below is still owed,
> and it is now a specific claim to test rather than a vague worry: *does
> registering a volumetric `add_dft_fields()` change the step count at which
> `stop_when_dft_decayed()` fires, on the 1×2, at its shipped
> `decay_tol=1e-11`?* Run the same section with and without it and compare step
> counts — not wall-clock, which confounds the DFT accumulation's own per-step
> cost with a genuine change in run length. Report both.

> ### ⚠⚠ ANSWERED WITHOUT A MEASUREMENT — and it changes Step 2's design
> **2026-08-31. Under the default `decay_method="dft"`, field capture WILL
> extend the run, and the library's own docstring says so.**
>
> `FDTDSolver`'s `decay_method` docstring (`fdtd.py:962-988`) states, flagged as
> *"confirmed directly against meep's own docs/source"*:
>
> > `mp.stop_when_dft_decayed()` **requires *every* field component and
> > frequency on *every* DFT monitor to individually decay below `decay_tol`**
> > … No `maximum_run_time` is ever set, so there is no built-in ceiling: if any
> > one component decays anomalously slowly, the run has nothing to stop it
> > early.
>
> A volumetric `add_dft_fields()` registers an enormous amount of new DFT
> content — every requested component at every voxel, including near-PML
> regions and corners where energy lingers — and **all of it** must individually
> clear `decay_tol`. So Step 2's opt-in flag, implemented naively on the "dft"
> path, is not "adds a plot, doesn't change the physics run": it is a plausible
> route to the same open-ended stall that motivated `decay_method` in the first
> place (the 4×4's body calibration never triggered past **660,000+**
> timesteps).
>
> **`decay_method="energy"` removes the problem structurally rather than
> measuring it away.** `mp.stop_when_energy_decayed(dt, decay_by)` — read
> directly from meep's source — *"examines the field energy over the entire cell
> volume at every `dt` time units"* and stops when it has decayed by `decay_by`
> from its **maximum recorded** value. That criterion **does not reference
> registered DFT monitors at all**. Consequences, both of which Step 2 wants:
>
> 1. **Run length becomes provably independent of field capture.** Not "measured
>    to be within noise" — independent by construction, since the monitor is not
>    an input to the stop condition.
> 2. **The fingerprint dilemma dissolves.** The two-part rule keys on whether the
>    opt-in makes it a genuinely different simulation. If the stop condition
>    cannot see the monitor, capture does not change the run, so the flag can be
>    **default-omitted** and every existing cache entry stays valid.
>
> **Cost of that route, stated honestly:** `decay_method="energy"` is labelled
> *"investigative, not yet validated as a general replacement"*, and
> `energy_decay_dt=20.0` has *"no principled default yet … not a
> measured/validated choice"*. So it trades Step 2's decay-timing question for a
> decay-**method** validation question. That is a better trade — the second is a
> one-time, device-independent piece of work that several other open items
> (the 4×4 stall, the floor) also want — but it is not free.
>
> **Also found, and worth acting on independently:**
> `stop_when_dft_decayed` accepts **`maximum_run_time`** (default `None`), and
> **beef never sets it anywhere**. It is an available ceiling against exactly
> the open-ended stall the docstring warns about, currently unused. The Step 0
> cold run's own 3.5 h against a 1:28 precedent is that risk showing up in
> practice.
>
> And a trap for any device switching methods: **`minimum_run_time` is silently
> ignored under `decay_method="energy"`** — `stop_when_energy_decayed` takes no
> such argument and `_stop_condition()` only forwards it on the "dft" branch.
> Not a hole (energy must rise before it can fall from its own maximum, so the
> criterion cannot trip early the way DFT relative-change can), but it means the
> `minimum_run_time` floor a loosened `decay_tol` requires does **not** carry
> over.

**Open question, not yet answered — needs measuring before trusting it**:
`stop_when_dft_decayed()` decays on EVERY registered DFT monitor. Adding a
field monitor could change decay timing (and therefore cost) on any call
that opts in. Given this device family's own decay-tolerance history (the
4×4's `decay_tol=1e-11` floor, [[decay-tol-has-a-hard-floor]]), this needs a
real before/after timing measurement on 1×2 specifically before being
trusted as "just adds a plot, doesn't change the physics run."

## Step 3 — stitch input/taper-section fields + body field into one map

**Scope, stated explicitly since it was asked**: the deliverable spans **all
three sections**, not just the EME body. For the 1×2 the split is

| section | solver | span (µm) | share of device |
|---|---|---|---|
| input access wg + taper (+0.25 stub) | **FDTD** | −3.25 → −2.00 | 19% |
| body | **EME** | −2.00 → +2.00 | 62% |
| output taper + branches (+0.25 stub) | **FDTD** | +2.00 → +3.25 | 19% |

(`total_length` 6.50 µm, `boundaries = [−3.25, −2.0, +2.0, +3.25]`.) Step 1's
body-only picture is an interim; two of the three sections are FDTD, which is
why Step 2 is the load-bearing library change.

Real new work: the hybrid cascade currently treats each `solve_section()`
boundary as an opaque transfer-matrix interface with no shared spatial
frame. Aligning the transverse (y) coordinate system and the propagation
phase reference across those boundaries, so the combined picture is a
physically continuous field map rather than three independently-normalized
panels, hasn't been attempted anywhere in this repo.

> **⚠ Found 2026-08-31, and it changes how the body half must be built: under
> the shipped config the EME body section computes NO mode fields at all.**
>
> `calibrate_body_neff=True` has been the default since 2026-08-17, so the
> 1×2's body section is solved through `EMESolver.solve_section()`'s
> `neff_override` **fast path** — and that path, by its own comment
> (`eme.py:680-690`), *"skips the eigenmode solve and `mw.compute_s_matrix()`
> outright rather than run and discarded."* The body's transfer matrix is
> built analytically as a pure diagonal
> `T[k,k] = exp(2j·π·neff_override[k]/λ·L)`.
>
> Two consequences, pulling in opposite directions:
>
> - **Worse than this doc assumed for Step 1/3 provenance.** The preamble says
>   the cache "persists arrays, never the meow `Mode` objects" — true, but the
>   stronger statement is that under the shipped configuration those mode
>   objects are **never computed in the first place**, exactly like the FDTD
>   side. The body field cannot be recovered from the cascade at any price;
>   it requires a separate mode solve, which is precisely what Step 1's
>   retarget does. (See `resolved/calibrate_body_neff_mode_solve_skip_
>   validation.md` for the skip's own validation.)
> - **Better than assumed for Step 3's hard part.** A diagonal body transfer
>   matrix means **no mode mixing across the body** — each interface mode
>   propagates independently, accumulating only phase. So once the body's mode
>   profiles exist (one solve), the body region of the map is *analytically*
>   reconstructible, and Step 3's genuine difficulty narrows to the two
>   FDTD↔EME boundaries and the FDTD interiors. The guard that enforces this
>   (`assert_uniform_cross_section`, both in `hybrid.py:1534` and inside
>   `solve_section()`) also means the property is structural for this device,
>   not a coincidence — a tapered body would refuse the fast path outright.

## Step 4 — validate against 1×2's own existing ground truth

`resolved/mmi_1x2_validation.md`'s gaps 4/5 (engine agreement, hybrid vs
FDTD) already have real measured full-device references for this exact
device. Once Steps 1-3 produce a reconstructed field map, integrate it back
to a transfer coefficient and compare against those already-trusted numbers
— a real quantitative check, not "does the picture look plausible."

## Step 5 — port to the 4×4 device (the actual motivating case)

Only after 1×2 proves the mechanism end to end. Same `HybridSolver`/
`solve_section()` machinery, so Step 2's library code should need no
changes — only Step 3's stitching re-run at the 4×4's own geometry
(currently W=6.50, y_gap=0.51, mm_length=91.00, the BO campaign's best
point as of round 5; re-check the campaign's own latest result before
committing to that exact point, since rounds are still running).

## Order

Revised 2026-08-31 after checking every step's premise against current code.

0. **Re-baseline 1×2 under current code** — *in progress.* Warm run done and
   bitwise identical (job 12627197), but that is the non-answer; the cold
   run (job **12627227**) is the real test and is in flight. Blocks trusting
   anything built on top.
1. ✅ **DONE 2026-08-31** — retargeted at the hybrid's own body section
   (`meep`/8 modes/±2.0 µm) rather than authored: the script existed since
   2026-08-08. Two scripts delivered, 55 s each, three self-checks passed
   (mode identity to 6 dp, exact parity, self-image landing on the ports).
   Standalone — needs nothing from step 2. Exposed two `plot_eme_fields()`
   interpretation defects (meow-frame labels; single-mode-only excitation)
   that step 3 must not inherit.
2. Add opt-in FDTD field capture to `solve_section()` (real library change;
   both premises verified, `region_cell_lengths` is the fingerprint
   template). **Measure decay timing before choosing the fingerprint
   treatment** — that measurement decides whether existing cache entries stay
   valid.
3. Stitch sections into one whole-device map for 1×2. **Still the genuinely
   hard, unattempted step** — nothing in the re-check reduced its scope.
4. Validate against 1×2's own engine-agreement/hybrid-vs-FDTD numbers.
   ⚠ Gated on step 0's outcome: if the cold run moves, this step's ground
   truth moves with it.
5. Port to the 4×4 device's best BO point. ⚠ The cited point (W=6.50 /
   y_gap=0.51 / mm_length=91.00) is **stale as written** — the BO campaign
   is still live (jobs `beef_bo_W6.00_y0.51`, `beef_bo_W6.50_y0.56`,
   `beef_bo_W5.00_y0.43` in queue 2026-08-31). Re-read the campaign's own
   latest result at porting time, as this doc's Step 5 already warns.

**Net effect of the 2026-08-31 re-check**: step 1 shrank from "build a
prototype" to "retarget a proven script"; step 2's premises were confirmed
intact and its fingerprint decision made explicit; step 0's warm run was
exposed as a non-answer and replaced with a real cold test; steps 3-5
unchanged in scope.

See also [[fdtd-solve-section-near-degenerate-normalization-fix]],
[[verification-standard-measured-not-argued]],
[[mmi-4x4-90-hybrid-status]].

### ✅ STEP 0 CLOSED — 2026-08-31, job 12674181 (1:57:42, cold, 24 MPI ranks)

**The cold re-baseline landed and the cached entries were VALID. 1×2 is a
trustworthy base for Steps 1-4.**

Fully cold (`example_sim_cache_step0_mpi`, fresh), current code, shipped config
reproduced exactly (`MMI1X2_DECAY_TOL=1e-11`, `MMI1X2_MIN_RUN_TIME=0.0`, res 40,
`N_MODES=8`, `calibrate_body_neff=True`). 12 domains, 11 time-stepping runs,
6 cache entries — including the `solve_section_multi_length` calibration entry
(`n_excite=8`) that the timed-out attempt never reached.

**Discriminator, as pre-registered: NO `normalization mismatch` warning
anywhere** in the log. The 4× Δn_eff margin held, now measured cold rather than
argued from the proxy.

**Numbers vs the preserved `mmi_1x2_hybrid_PRE_STEP0_20260818/`:**

| λ (µm) | 1.260 | 1.310 | 1.360 |
|---|---|---|---|
| `out_top` | 0.431112 | 0.374899 | 0.312504 |
| `out_bot` | 0.431101 | 0.374884 | 0.312493 |

**max \|Δ\| = 1.9e-10** across all 11 wavelengths, both ports.

> ⚠ **The pre-registered discriminator said "match bitwise". It does NOT match
> bitwise — and that is the CORRECT outcome, not a partial failure.** This run
> used **24 MPI ranks** against the preserved result's single rank, so domain
> decomposition changes the floating-point summation order; ~1e-10 is exactly
> FP reassociation on ~10-significant-figure values, not a physics difference.
>
> **A bitwise-identical result here would have been the WORSE outcome**: it
> would have meant the run served cache instead of recomputing. The earlier
> warm run (12627197) *was* bitwise identical for precisely that reason and
> proved nothing. The tiny non-zero residual is the evidence that this run
> genuinely did the work. The discriminator was written before the MPI change
> and should have said "match to FP tolerance"; recorded rather than silently
> reinterpreted.

**Consequence:** `resolved/mmi_1x2_validation.md` gaps 4/5 stand as written. No
re-measurement needed. Step 4's ground truth is intact.

### ⚠ The timed-out attempt was a RESOURCING bug, not a physics or cost one

Job 12627227 was submitted `--ntasks=1 --cpus-per-task=24` + `OMP_NUM_THREADS=24`
+ plain `python` — **one MPI rank**. meep parallelises FDTD time-stepping by MPI
domain decomposition, so 24 OpenMP threads bought almost nothing for the
dominant cost ([[mpi-ranks-calibration-vs-fdtd]]: stepping scales ~98% with
RANKS, MPB is rank-redundant). Same script, same tolerance, 24 ranks + `srun`:
**1:57:42 versus not finishing in 5:00:00**, and it matched the walled job's
entire output in 51 minutes.

> **Three wrong explanations were given for that wall before the right one**:
> (a) `decay_tol` too tight, (b) the floor making 1e-11 unreachable, (c)
> calibration cost. All were plausible and all were wrong. The submission was
> never checked against this repo's own documented convention (24 ranks +
> `srun`), which was written down. **Check the resourcing before theorising
> about the physics.**

`batch/step0_mmi_1x2_cold_mpi.batch` is the correct template for this device.

## ✅ STEP 2 BUILT AND MEASURED — 2026-08-31

`FDTDSolver.solve_section(capture_fields_dir=...)` — opt-in spatial field
capture, plus `FDTDSolver._write_section_fields()`.

**What it captures.** A 2-D **z-slice** at the core plane over the inner cell
(PML excluded), at the section's own excitation centre frequency `fcen`, all
three components (Ex, Ey, Ez), one `.npz` per excited band. Deliberately NOT a
3-D volume and NOT broadband: a volume is enormous at these resolutions and the
pictures only want the in-plane view, while a slice adds ~1/N_z of a volume's
DFT content, which keeps the decay-timing risk small.

**Coordinates are already DEVICE-frame and stitchable.** The npz carries
`x = (x0 - stub_left) + linspace(...)`, verified on the 1×2 input section:
`x range [-3.250, -2.000]`, exactly that section's own bounds, not a local
0-based frame. That shared coordinate is precisely what Step 3 needs. Also
recorded: `frame="device"` — and unlike the EME side there is genuinely no axis
remap needed here, because `FDTDSolver` builds x as the propagation axis
directly, so meep's Ex/Ey/Ez ARE the device's.

### Non-invasiveness — proven, not asserted

Fingerprint is **additive**, following `measure_diagonal_only`'s precedent
verbatim (`call["capture_fields"] = True` only when enabled). Evidence:

- `810 passed, 14 skipped` — full suite, unchanged.
- The 6 entries in `example_sim_cache_step0_mpi`, written by the PRE-change
  code path, carry no `capture_fields` key.
- **A warm non-capture re-run against that cache took 16.5 s and reproduced
  T=0.374899 / total 0.7498 exactly** — against the cold run's 1:57:42. Every
  prior entry still hits.

### ✅ The open decay-timing question, ANSWERED by measurement

Step counts (not wall-clock, which confounds per-step DFT cost with a real
change in run length), 1×2 input section, res 15, 2-D, no save_dir either side:

| `decay_method` | `decay_tol` | no capture | capture | Δ |
|---|---|---|---|---|
| dft | 1e-6 | 2161 | 2161 | 0.00% |
| **dft** | **1e-8** | **2641** | **2521** | **−4.54%** |
| energy | 1e-6 | 2404 | 2404 | 0.00% |
| energy | 1e-8 | 2404 | 2404 | 0.00% |

**`decay_method="energy"`: exactly zero at both tolerances.** The structural
argument holds — that criterion reads total cell energy and never inspects DFT
objects, so capture cannot move it. Now measured, not just reasoned.

**`decay_method="dft"`: capture DOES change the run, by −4.5% at 1e-8.** So the
flag genuinely must be fingerprinted under the shipped default, which is what
the additive key does.

> ⚠ **The direction of the prediction was WRONG and is corrected here.** This
> doc argued (and `decay_method`'s own docstring supports) that a volumetric
> monitor would make runs LONGER, since `stop_when_dft_decayed()` requires every
> registered DFT component to individually clear `decay_tol`. Measured, capture
> made the run **shorter** (2641 → 2521). The conclusion — capture alters the
> dft-criterion run and must be fingerprinted — survives; the mechanism behind
> it does not, and remains **unexplained**. Plausibly the criterion's
> relative-change measure is diluted by a large fast-converging object, but that
> is a guess, not a measurement. Do not repeat the "it will be slower"
> reasoning as though it were established.

> ⚠ **Separate observation: the energy criterion's stop time is QUANTIZED.**
> It gave 2404 steps at BOTH 1e-6 and 1e-8 — `energy_decay_dt=20.0` means it
> only checks every 20 time units, so both tolerances were crossed inside one
> window and the tolerance had no effect. That is a real limitation of the
> energy route (it can overshoot the tolerance actually needed) and it compounds
> `energy_decay_dt`'s own "no principled default yet" status.

### What Step 2 does NOT yet have

- **Measured on a 2-D section at res 15 only.** The 3-D production path is
  untested for capture; the shape guard in `_write_section_fields` warns rather
  than crashing if the slice comes back degenerate, but that is a guard, not a
  measurement.
- **A cache hit with capture requested writes no fields** (the early return
  precedes the run). The entry is fingerprinted separately so a first capture
  run always executes, but a SECOND identical capture run will hit cache and
  silently produce no new npz. Related to
  [[cache-only-probe-still-writes-outputs]]'s inverse. Not yet handled.
- **No plotting/stitching consumer** — that is Step 3.

## STEP 3 — stitching, in progress 2026-08-31

New module `beef/plot/hybrid_stitch.py`: `load_captured_section()`,
`dominant_component()`, `project_interface_onto_body()`,
`stitch_hybrid_fields()`. Driver:
`examples/verification/python/stitch_mmi1x2_hybrid_fields.py`.

**The core idea, and why it dissolves the two problems this step was opened
for.** Step 1 had to invent the body's launch. Step 2 now captures the real
interface field, so the body's modal coefficients are PROJECTED from it. Then:

* **amplitude alignment** — the coefficients carry the captured field's own
  absolute units, so body and FDTD panels share one scale. No per-panel
  normalisation to reconcile.
* **propagation phase reference** — the coefficients are measured AT the shared
  boundary and the body propagator is referenced to that same plane, so the
  body's phase continues the FDTD section's by construction.

Neither needed the bespoke alignment machinery this step's original text
anticipated; both fall out of projecting at the seam.

### ⚠ Step 3 revealed Step 2 was INCOMPLETE — a second capture was required

Step 2 as first built captured only a **longitudinal x-y slice** at the core
z-plane. Its interface column is therefore a **1-D line in y**, and projecting
a 1-D cut onto the mode basis is exactly the pathology
`body_mode_weights()` documents from measurement (captured power **1.5857**,
`|c|^2 = 0.5876` on a parity-forbidden mode). So Step 2 was extended to
additionally capture the **transverse plane at the measure boundary**, reusing
the measure monitor's own region verbatim so the projection plane coincides
with where the transfer matrix is evaluated. Stored under an `if_` prefix.

**A frame bug was caught doing this, and it would have been invisible.**
`mea_along` is a **shifted-sim** coordinate (`solve_section` builds the section
centred on the origin), while the longitudinal `x` array is device-frame. The
1×2 input section's plane first recorded as **sim x = +0.625** — a position
*outside* the section's own `[-3.250, -2.000]` span. Converted at write time
(`cen["x"] += x_origin + span_x/2`); it now reads **−2.000**, exactly `x1`, and
`check_solve_section_field_capture.py` asserts the plane lands on a section
boundary so this cannot regress silently.

### ⚠ THE OUTPUT SECTION IS DELIBERATELY NOT STITCHED

`build_run_config()` sets `reverse = n_out <= n_in` (`fdtd.py:735`) and excites
from whichever boundary carries FEWER modes. The 1×2 output section has
`n_modes_in=8` (body side) and `n_modes_out=1` (port side), so it is **excited
from its output port** and measured at the body side.

Reciprocity keeps its transfer MATRIX valid — that is the whole point of the
optimisation — but its captured **field map is of light launched backwards from
the output port**, not of the field the body delivers forward into it.
Stitching it would produce a picture that looks complete and describes a
different experiment. A forward output map needs that section re-run excited
from the body side: 8 excites, so 8 FDTD runs instead of 1. Real cost, deferred.

**So the Step 3 deliverable is input + body**, which is the genuinely
continuous, real-field-driven part of the chain.

### Checks built into the driver

* **Parity.** A centred single-mode input must drive only EVEN body modes. This
  is Step 1's parity check applied to a REAL field, so it also tests that the
  capture grid and the mode mesh are aligned in y — a misalignment would leak
  power into odd modes.
* **Seam continuity, measured not eyeballed.** Integrated |E|² on the FDTD side
  of the seam versus the body side, in absolute units on a shared scale.
  ⚠ **The ratio is NOT expected to be 1.0**: the body side is a truncated-basis
  reconstruction, so it is missing exactly the power `captured_power` reports as
  uncaptured, plus any non-guided content the FDTD slice still carries there.
  The number to quote is `captured_power`, not the picture.
* Grids are **resampled**, not assumed equal — `FDTDSolver.sim_width/sim_height`
  and `EMESolver.oxide_width/oxide_height` are independent parameters and
  several devices set them differently.

### ✅ STEP 3 CLOSED at production resolution — 2026-09-02, job 12736849

**5:25 wall / 43 CPU-min on 8 ranks.** Outputs:
`example_output/mmi_1x2/mmi_1x2_hybrid_stitched_res40/`.

`n_eff` reproduces Step 1's and gap 3's values exactly (2.696915, 2.637710,
2.536282, 2.387965, 2.185091, 1.916790, 1.775436, 1.706013), so the body basis
is the same one those results were built on.

| quantity | res 20 | **res 40** |
|---|---|---|
| interface power in the 8-mode basis | 0.9618 | **0.9818** |
| odd-mode fraction (parity) | 2.8e-06 | **1.9e-09** |
| seam continuity ratio | 0.9520 | 0.9472 |
| dominant interface component | Ey 11.79 | Ey 13.12 (Ez 0.39) |

### 🔴 THE RESULT THAT MATTERS: Step 1's Gaussian launch was substantially WRONG

Like-for-like at res 40, body modal content of the launch:

| body mode | Gaussian (Step 1) | REAL captured field | ratio |
|---|---|---|---|
| 0 | 0.9022 | **0.6877** | 0.76× |
| 2 | 0.0935 | **0.2768** | **3.0×** |
| 4 | 0.0013 | **0.0354** | **27×** |

**The real interface field is far more multimode than a taper-width Gaussian
implies.** Step 1's body pictures are qualitatively indicative and
quantitatively wrong; anything read off them about mode content or beat
structure should be re-read from here. This is exactly what Step 3 existed to
establish, and it argues the approximate-launch caveat in
`beef/plot/hybrid_body.py` should be read as a strong warning rather than a
formality.

### ✅ STEP 4's core check, obtained free from the same run

The input section's transfer-matrix column and the projected coefficients are
the SAME physical quantity by two INDEPENDENT routes — meep's
`get_eigenmode_coefficients` at the flux plane (meep's own eigensolver) versus
DFT field capture projected onto the EME/MPB mode basis. Normalised modal
power:

| mode | \|T\|² | \|c\|² | Δ |
|---|---|---|---|
| 0 | 0.666958 | 0.687703 | +2.07e-02 |
| 2 | 0.285820 | 0.276850 | −8.97e-03 |
| 4 | 0.047222 | 0.035447 | −1.18e-02 |
| 1,3,5,6,7 | 0.000000 | 0.000000 | ~1e-10 |

**max |Δ| 2.1%, sum |Δ| 4.1%; parity-forbidden modes agree to 1e-10.**
Different mode solvers and different extraction methods agreeing to 2% is a
real validation of the capture against the number the cascade actually uses.
The residual is consistent with the **1.8%** of interface power the truncated
basis does not capture, plus grid interpolation between the FDTD capture mesh
and the EME mesh.

⚠ Not the whole of Step 4: this validates the capture against the *input
section's* transfer matrix, not the full-device hybrid-vs-FDTD numbers of
gaps 4/5. That comparison still needs the body propagated to its far boundary
and combined with the output section.

### 🔴 A 16-CPU-HOUR DEADLOCK I INTRODUCED — job 12736213

`_write_section_fields` called the **collective** `get_dft_array()` for the
interface plane *after* `if not mp.am_master(): return`. Rank 0 blocked forever
in a collective the other 23 ranks had already left: `futex_wait_queue`,
24 procs pinned at 99% CPU, log frozen at exactly 4430 bytes, capture dir
empty. Killed at 40:35 having burned **cput 16:14:00** for zero output.

**Three compounding errors, recorded because each is repeatable:**

1. **The rule was already written in that method's own docstring** ("`get_dft_array()`
   is collective and must be called on every rank") and I violated it by
   appending the Step 3 extraction past the guard.
2. **Serial testing cannot catch it.** With one rank `am_master()` is always
   true. `check_solve_section_field_capture.py` now documents and requires
   `mpirun -np 4` (~1.4 s) alongside the serial run.
3. **I explained the evidence away.** With the empty capture dir, the futex
   wait, and the frozen log all in hand, I constructed a cost story about
   `get_eigenmode_coefficients` "scaling worse than cell count at res 40" and
   reported it as the likely cause. It was 100% hang. **Prefer the hypothesis
   the instrumentation already supports over a plausible performance narrative.**

Fixed by hoisting every `get_dft_array()` above the `am_master()` guard and
making the failure path non-short-circuiting (an early `return` inside the
extraction loop desynchronises ranks the same way). Verified under
`mpirun -np 4`, then in batch at 8 ranks.

### ⚠ Monitoring lessons from the same episode

* **`grep -c` on a MISSING file returns 0 and exits 1**, so every log-derived
  counter in the rerun's monitor read 0 — the job's log was at
  `BEEF/reports/` (sbatch resolves `-o reports/%j.out` against
  `SLURM_SUBMIT_DIR`, which was `BEEF/examples`), not `duty/reports/`. The
  "run finished but no npz" canary was therefore **vacuous**: `fin` could never
  exceed 0. Only the npz check, whose path happened to be right, carried
  information. Same failure class as the earlier `grep -ci ... || echo 0`
  producing the two-line string `"0\n0"`.
* **The npz is a better progress signal than the log** for this pipeline: it is
  written immediately after `sim.run()` while meep's stdout is 4 KB
  block-buffered through `srun`. Buffering also means a log-growth stall
  detector will report a false stall on a healthy job.

### ⚠ Reverse excitation and the output panel — MEASURED, job 12738781

The claim that the output section's reverse-excited capture cannot serve as the
forward field panel was tested rather than argued.
`capture_mmi1x2_output_section.py` ran the output section exactly as
`solve_section` does (8 modes body side, 1 port side → `reverse=True`, excited
from the port), res 40, 3:57.

**Exactly ONE npz was written** — one excite, from the single-mode port. And
the power split at the port end (x = +3.25):

| region | |E|² |
|---|---|
| y > +0.2 (top branch) | 2.4799e+02 |
| y < −0.2 (bottom branch) | 1.3058e-01 |
| **ratio** | **1899** |

Forward 50/50 operation must give ≈ 1. Figure:
`fields/reverse_output_section_why_not.png` — one branch lit, the other black.

**Why, precisely.** The forward field is an 8-term superposition
Σ_k a_k E^(k) over the body modes, with a_k the amplitudes the body actually
delivers (measured 0.6877 / 0.2768 / 0.0354 at the input seam). Reverse
excitation yields ONE distribution. Time-reversal does make the reverse
field's conjugate a valid FORWARD field — but for the input state ∝ S_1k*,
i.e. the superposition that couples optimally into port 1, not for a_k. Those
coincide only if a_k ∝ S_1k*, which is not the case here.

**None of this weakens reverse excitation for its actual purpose.** The
transfer matrix is exact by reciprocity and the 8× saving is real; the repo
should keep using it. The limitation is specific to field maps.

**Route to the remaining 19% of the device**, costed: excite the output section
from the body side once per body mode (Step 2 already writes one npz per
excited band), then superpose weighted by a_k. Needs a small library change —
`solve_section`'s `reverse` is derived from the mode counts, so it needs an
opt-in override to force forward excitation — plus 8 excites instead of 1
(~30-40 min at res 40, from this run's 3:57 for one).

### Field figures produced (res 40, λ = 1.31 µm)

* `fields/detail_input_body_res40.png` — the readable version: input section,
  body, transverse cuts, and the modal-content comparison.
  * **The seam validation is panel (c)**: the FDTD-side profile and the
    body-side reconstruction lie exactly on top of each other, which is the
    visual counterpart of "0.9818 of interface power captured".
  * Body exit peaks at **|y| = 0.538 µm against ports at 0.500** — the
    self-image lands 0.038 µm off the port centres, driven by the real field.
  * Input-section profile shows the expected `Ey` spikes at y ≈ ±0.18, the
    normal-component discontinuity at the waveguide sidewalls.
* `fields/hybrid_stitched_mmi1x2_res40_1.310um.png` — the combined map.
  ⚠ **Its title overstates coverage**: it reads "stitched hybrid field map"
  but covers input + body only, **81%** of the device (x = −3.25 → +2.00 of
  −3.25 → +3.25). Retitle to "input + body".

### Whole-device map assembled — 2026-09-14, job 13189964 (22:05, 8 ranks)

**100% coverage.** The output section's forward excitation needed NO library
change: `fdtd.py:3232` already does `reverse = False if override else
layout.reverse`, so passing any `excite_port` PortSpec forces forward launch.
`PortSpec(axis="x", sign=+1, center=0.0, width=5.0)` (sign=+1 binds an EXCITE
port to the LOW edge = the body side) turned `reverse=True, n_excite=1` into
`reverse=False, n_excite=8`, verified with `section_fdtd_layout()` before
spending the runs. 8 captures, superposed with the body's arriving amplitudes
a_k = c_k·exp(2j·pi·n_eff_k·L/wl).

Figure: `fields/whole_device_res40_1.310um.png`.

| check | result |
|---|---|
| coverage | x −3.25 → +3.25, **100%** |
| output lobes at the port end | **\|y\| = 0.510 µm** vs ports at 0.500 |
| branch balance | **1.000** (reverse-excited capture gave 1899) |
| body/output seam step | 1.23× |

### 🔴 TWO CLAIMS OF MINE THAT DID NOT SURVIVE — recorded so they are not reused

**(1) "The 50/50 branch check validates the superposition." It does not.**
This device is symmetric and only EVEN body modes are excited, so ANY even-mode
combination produces a symmetric output. Ratio 1.000 confirms parity survives
end to end — real, but it cannot detect a wrong relative weighting between
modes 0/2/4. A necessary check presented as a sufficient one.

**(2) "∫|E|²dy grows 10.5× along a passive section, so the amplitude is
demonstrably wrong." OVER-CLAIMED — that quantity is not power.**
The captures are a single z-SLICE through the core. The body modes are
orthogonal under the full 2-D (y,z) integral, NOT under a 1-D y-cut, so
cross-terms survive and the slice integral has no conservation law — the same
1-D-vs-2-D distinction `body_mode_weights()` documents after a 1-D projection
returned captured power 1.5857. The output section also tapers into two narrow
guides, so its local mode shape changes along x. Growth in that metric is
therefore not evidence of a defect.

**The fix built on (2) was tested and FAILED, which is what exposed (1) and
(2).** `fix_mmi1x2_output_normalisation.py` renormalised each E^(k) so its
entry plane equals phi_k, cancelling meep's per-band source convention.
Growth moved only **10.51× → 9.57×**. Cause: the apparent "36292× spread" in
per-band amplitudes is entirely modes 6–7, which couple at ~1e-4 and carry no
power (a_k ≈ 0 by parity). The modes that matter — 0, 2, 4 — span just
**1.2×** (7.150 / 8.198 / 6.755), so per-mode normalisation could never have
been a 10× effect.

**Current honest status of the output panel: NOT validated.** Parity (1.000)
and lobe position (0.510 vs 0.500) are checked; the amplitude is neither
confirmed nor refuted, and the growth is undiagnosed.

**The check that would settle it** — and which remains Step 4 proper — is to
integrate the reconstructed output to a transfer coefficient and compare
against the hybrid's own validated **0.3749** per port at 1.31 µm (Step 0,
job 12674181). That is a conserved quantity with a trusted reference, unlike
any slice metric invented for the purpose.
