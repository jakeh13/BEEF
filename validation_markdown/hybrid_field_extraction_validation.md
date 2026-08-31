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
1. **Retarget** `plot_mmi1x2_eme_fields.py` at the hybrid's own body section
   (`meep`/8 modes/±2.0 µm) — **not** author it: the script exists since
   2026-08-08 and its pure-EME output is already on disk. ~20-line change to
   a proven script; no longer a prototype milestone.
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
