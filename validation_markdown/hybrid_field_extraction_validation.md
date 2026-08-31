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
cached arrays.** Both its existing call sites
(`plot_mmi1x4_eme_fields.py`, `plot_mmi4x4_90hybrid_eme_fields.py`) pass it
meow's own `cells`/`modes` objects — each mode object carries its full
per-cross-section field grid, computed fresh by that call. `HybridSolver`'s
on-disk cache (`save_dir`) persists n_eff/transfer-matrix **arrays**
(`beef/cache/`), never the meow `Mode` objects themselves. So even the
"cheap" half of this plan — the EME body, whose length-sweep is free by
`solve_section_length_sweep()`'s own analytic-reconstruction tier-3 cache
— still needs one fresh (if modest) EME mode solve to get real field data
back. Nothing here is a pure cache read.

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

## Step 1 — prototype the free-ish half: EME body field reconstruction

Adapt `plot_eme_fields()` (already exists, already proven on two other
devices) to 1×2's own body cross-section — no companion script exists for
this device yet. Real cost is one EME mode solve at 1×2's body geometry,
bounded by the cheapest solve of any hybrid-cascade device measured so far
in this repo. Gives a body-only field picture, and importantly exercises
the "combine cascade-cached coefficients with a fresh mode-object solve"
mechanics this whole plan depends on, before also touching the FDTD side.
Not the deliverable on its own (body only, and it's the same construction
already flagged above as a poor stand-in for the hybrid path's own trusted
numbers) — a stepping stone to de-risk Steps 3-4.

## Step 2 — add spatial field capture to `solve_section()`'s FDTD runs

New library code, following the exact non-invasive pattern already used for
`decay_method` (`beef/cache/config.py`'s two-part fix: default fingerprint
untouched so every existing device's cache is unaffected, a new opt-in
value gets its own distinct fingerprint). Proposed: an opt-in flag that adds
`mp.add_dft_fields()` at the section's own excitation frequency before
`sim.run()`, then `sim.get_array()` once the existing decay condition fires
— alongside, not instead of, the existing eigenmode-coefficient extraction
that produces the transfer matrix.

**Open question, not yet answered — needs measuring before trusting it**:
`stop_when_dft_decayed()` decays on EVERY registered DFT monitor. Adding a
field monitor could change decay timing (and therefore cost) on any call
that opts in. Given this device family's own decay-tolerance history (the
4×4's `decay_tol=1e-11` floor, [[decay-tol-has-a-hard-floor]]), this needs a
real before/after timing measurement on 1×2 specifically before being
trusted as "just adds a plot, doesn't change the physics run."

## Step 3 — stitch input/taper-section fields + body field into one map

Real new work: the hybrid cascade currently treats each `solve_section()`
boundary as an opaque transfer-matrix interface with no shared spatial
frame. Aligning the transverse (y) coordinate system and the propagation
phase reference across those boundaries, so the combined picture is a
physically continuous field map rather than three independently-normalized
panels, hasn't been attempted anywhere in this repo.

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

0. Re-baseline 1×2 under current code (open, blocks trusting anything built
   on top of it).
1. Prototype EME body field reconstruction for 1×2 (cheap, de-risks the
   cache/mode-object mechanics).
2. Add opt-in FDTD field capture to `solve_section()` (real library change,
   needs a decay-timing measurement before trusting it doesn't change cost).
3. Stitch sections into one whole-device map for 1×2.
4. Validate against 1×2's own existing engine-agreement/hybrid-vs-FDTD
   numbers.
5. Port to the 4×4 device's best BO point.

See also [[fdtd-solve-section-near-degenerate-normalization-fix]],
[[verification-standard-measured-not-argued]],
[[mmi-4x4-90-hybrid-status]].
