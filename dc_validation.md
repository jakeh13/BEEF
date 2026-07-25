# Directional coupler validation — to-do

Opened 2026-07-25, from investigating why `directional_coupler_hybrid_stub_sweep.py`
never converges. Root cause found; this is the work to fix it and re-establish trust
in the directional coupler results.

**Status 2026-07-25:** items 1–4 done and verified — the wavelength axis is fixed and
confirmed against a pre-fix run (job 11477911, 9/9 checks), and the cache can no longer
serve or corrupt stale data. What remains is re-establishing the reference data (5, 9)
and the residual ~11% engine disagreement (6), which items 1–2 did not touch.

## What's wrong

> Item 1 is **fixed** as of 2026-07-25 — see Done at the end. This section is kept
> as the description of the defect the rest of the list is scoped against.

**`FDTDSolver` mislabelled its own wavelength axis.** `solvers/fdtd.py:352` and `:654`
computed

```python
fwidth = (1.0/wavelengths.min() - 1.0/wavelengths.max()) * 0.5
```

and pass it as meep's `df` to `add_mode_monitor(fcen, df, nfreq)` (`:411`, `:424`,
`:703`, `:708`). Meep spaces `nfreq` frequencies evenly over `[fcen - df/2,
fcen + df/2]`, **ascending in frequency**. For the standard 1.50–1.60 µm, 11-point
C-band sweep that is:

| index | code labels it | meep actually computed | error |
|---|---|---|---|
| 0 | 1.5000 µm | 1.5754 µm | +75.4 nm |
| 5 | 1.5500 µm | 1.5500 µm | 0 |
| 10 | 1.6000 µm | 1.5254 µm | −74.6 nm |

Two independent defects in one line:

1. **Half the span.** `df` is set to *half* the requested frequency width, so the
   monitor covers 1.5254–1.5754 µm, not 1.50–1.60 µm.
2. **Reversed.** Meep's grid ascends in frequency, i.e. *descends* in wavelength, but
   the result is indexed against the caller's ascending-wavelength array.

Only the centre point (index 5, `fcen`) is correctly labelled. Verified directly
against `mp.get_flux_freqs()`.

Note a linear-in-λ array can never be reproduced by a `(fcen, df, nfreq)` triple —
meep spaces linearly in frequency. The fix must pass an explicit frequency list.

### What this explains

- **Every FDTD directional-coupler result shows backwards wavelength dispersion.**
  Archived full-device sweep: 50:50 length 13.16 µm at "1.50" rising to 17.65 µm at
  "1.60". Coupling appears to *weaken* with wavelength, which is backwards for a
  directional coupler. Unflip the axis and it is the physically correct direction.
  Reproduced in the archived FDTD *resolution* sweep at all 8 resolutions, so it is
  not a mesh artifact.
- **Why the stub sweep can't converge.** `np.einsum("Lwij,wjk->Lwik", M_body, M_in0)`
  contracts the EME body at true λ_w against an FDTD section at a *different*
  wavelength — at w=0 it multiplies EME@1.500 µm by FDTD@1.5754 µm. Ten of eleven
  columns are physically meaningless products. Each stub is its own FDTD run, so the
  wavelength-mean the sweep plots drifts smoothly with stub with no reason to ever
  settle: `L_50 = 0.2427*BODY_STUB + 13.4067 µm`, max residual to that line 2.2e-3 µm
  over 0→2 µm. Zero curvature.

### What it does *not* explain

Index 5 is correctly labelled in both engines, and the hybrid still disagrees with the
full-device FDTD there: 13.65 µm (stub 0) → 14.09 µm (stub 2) versus **15.29 µm**.
Restricting the stub sweep to λ=1.55 alone still gives a dead-straight
+0.2185 µm/µm drift. In the two-mode picture,
`L_50(s) = [π/2 − 2φ_bend]/Δβ_EME + 2s(1 − Δβ_FDTD/Δβ_EME)`, so that slope measures
`Δβ_FDTD/Δβ_EME = 0.891` — the FDTD sections accumulate ~11% less even/odd phase per
µm of the *same* straight coupled region than the EME body does. That is a second,
independent problem (item 6).

Correction to an earlier read of this: the narrow `fwidth` does **not** starve the
source. Source and monitor share the same variable, so the pulse is well matched to
the band meep actually sampled — power conservation is 0.9987–0.9998 at every point.
The bug is the labelling, not the signal level.

## Is the archived FDTD run a valid ground truth?

`archive/example_output_20260717/directional_coupler/directional_coupler_fdtd_sweep/length_sweep.npz`

**Partly — as a single-wavelength reference only.**

- Usable: power-conserving (0.9987–0.9998); resolution-converged — the 2026-07-22
  FDTD resolution sweep gives T_cross = 0.423 ± 0.008 at λ=1.55 across res 30–90, so
  res=40 is fine; and **its index 5 is correctly labelled**, giving a genuine
  full-device 50:50 length of **15.29 µm at 1.55 µm**.
- Not usable as-is for spectra: the same axis bug mislabels its other ten wavelength
  points. **They are recoverable, not garbage** — job 11477911's check 8 validated the
  relabelling recipe (map index *i* to `1/(fcen + df*(i/(n-1) − 1/2))`, i.e. reversed
  onto 1.5754→1.5254 µm) on a sibling run from the same code, agreeing with a fresh
  post-fix solve to 2.5e-4. So this file can be re-read on its true axis if needed.
  The catch is that its true axis spans only 1.5254–1.5754 µm, half the nominal band.
- Provenance is unrecorded. The npz stores only `coupling_lengths`, `wavelengths`,
  `T_through`, `T_cross` — no config, no solver fingerprint, no geometry, and the
  whole-device `solve()` path had no `save_dir` then (the archived cache holds only
  `solve_section`/`solve_section_length_sweep` entries). Its parameters match the
  current script's constants, but nothing in the file proves it.
- Coarse: 2 µm length steps, so 15.29 µm is interpolated across a 2 µm bracket
  (~±0.1 µm), against a stub sweep resolving 0.06 µm steps.

Verdict: use 15.29 µm as a provisional sanity check now; re-run the sweep on a finer
grid before treating any number as truth (item 5). Items 1–2 are now verified, so that
re-run is unblocked. Provenance and the 2 µm step, not the axis bug, are what keep this
file from being ground truth.

## To-do

Numbers are stable IDs — completed items move to Done and keep their number, so
cross-references stay valid.

**Order: the items-1–2 verification is done and passed** (job 11477911, 2026-07-25 —
see Verification below), so items 5, 6 and 7 now rest on a confirmed FDTD. Remaining
order: **5** → **9** → **6** → **7**, with **10** whenever the ring resonator matters
(it blocks nothing here). (**3** and **4** are done — they came first because they are
pure code, need no queue, and were what would otherwise have wasted wall-clock on the
big array jobs in 5 and 7.) Item 6 comes before 7 because 7's interpretation depends on
whether the residual ~11% gap is understood. That gap was measured at λ=1.55, the one
index both engines labelled correctly before the fix, so it survives the verification
intact and remains the main open physics question.

5. **Re-establish ground truth.** *Harness built 2026-07-25; sweep not yet run.*
   Re-run the full-device FDTD length sweep on a finer grid than 2 µm near the
   crossing, with the solver config saved into the npz so provenance is
   self-describing next time. This is what the stub sweep's reference line reads
   from —
   `example_output/directional_coupler/directional_coupler_fdtd_sweep/length_sweep.npz`
   does not currently exist in the live tree at all.

   Set up as a **verification example** rather than by editing
   `directional_coupler_fdtd_sweep.py`, which keeps its deliberate 0:2:60 grid and its
   teaching role:

       BEEF/examples/verification/python/verify_item_5_ground_truth.py
       batch/verification/run_verify_item_5_ground_truth.batch   (--array=0-39)
       make verify_item_5_ground_truth GT_INDEX=<0..39>

   **Grid — 40 lengths**, `union(0:4:60, 9:0.5:22)`. The fine window is placed from
   data, not guessed: the archived sweep re-read on its true axis gives
   dL₅₀/dλ ≈ −90 µm/µm, so across 1.50–1.60 µm the crossing runs ~19.8 µm down to
   ~10.8 µm. 0.5 µm there is 4× finer than the archive where it matters; 4 µm elsewhere
   pays for it and still resolves the oscillation. Device and solver kwargs are pinned
   to `directional_coupler_fdtd_sweep.py` exactly, so the result is directly comparable
   to both the archive and job 11477911.

   **Caching is ON**, unlike the items-1–2 harness. That one had to be uncached because
   a verification must not pass on a stored result; this one *produces* reference data
   at ~1–3 h per task, so a re-run after a node failure should be free. Only safe
   because item 3 completed the fingerprint and item 4 made the writes rank-safe —
   `--no-cache` forces re-simulation.

   Outputs to `example_output/verification/ground_truth_sweep/`: `length_sweep.npz`
   (grid + `L50` + bracketing lengths + embedded provenance), `provenance.json`,
   `length_sweep.png`, `l50_vs_wavelength.png`. L₅₀ is the first upward crossing of
   `T_cross = 0.5`, linearly interpolated; the **bracket is stored alongside it**
   because 15.30 µm quoted from a 2 µm bracket and from a 0.5 µm bracket are not the
   same claim, and the archived file gave no way to tell which it was. On synthetic
   data the 0.5 µm grid recovers L₅₀ to 0.015 µm.

   Status: aggregation, provenance, L₅₀ extraction and both plots validated on
   synthetic data; `make -n` and `bash -n` clean. **Validation task submitted as job
   11480098** (`--array=12`, L=13.5 µm — the length job 11477911 already solved
   uncached, so its `T_cross(1.55) = 0.4230` is an independent cross-check). Submit the
   remaining 39 only once that returns sane.
6. **Explain the residual ~11% rate gap at λ=1.55.** Not covered by 1–2. Add a direct
   engine-agreement test alongside `tests/test_hybrid.py`: take a pure straight
   coupled-pair region (no bends), call `FDTDSolver.solve_section()` at two lengths,
   extract Δβ_FDTD from the even/odd phase difference, compare against
   `EMESolver.compute_modes()`'s neff split on the same cross-section, assert to a
   stated tolerance. Expect it to fail at ~11% today. This is the engine-consistency
   check `notes.md` Solver #9 asks for and the cleanest way to isolate the effect —
   mesh discretization is already ruled out by the resolution sweep.
7. **Re-run the stub sweep — the whole thing — and judge it properly.** All nine array
   tasks, not just the two that item 4 lost: the seven `S*.npz` files on disk are dated
   2026-07-24, i.e. pre-fix, so re-running only tasks 4 and 5 (as item 4 originally
   scoped) would leave one folder holding a silent mixture of pre- and post-fix data.
   A visible gap is safer than an invisible mixture, so the two dead stubs are folded
   in here. Delete the folder first (`make clean_directional_coupler_hybrid_stub_sweep`)
   so nothing pre-fix survives. Success is *curvature* — successive
   differences shrinking toward the grid step — not the constant 0.060 µm/step it
   shows now. Only then does the sweep answer the question it was written to answer,
   and only then is `notes.md` Solver #10(b) (promoting this into a reusable
   procedure) worth acting on.
8. **Minor, while in there — `notes.md` Solver #9's material half.** Confirmed present
   in the data: EME's Δn_eff has a +0.00022 kink at exactly 1.57 µm where meow's Si
   table switches 3.4784 → 3.4710, while the FDTD half stays frozen at 3.4784 for
   every wavelength. Real, but ~0.8% of Δn_eff — worth fixing for correctness, not a
   driver of anything above. (Solver #8, wavelength interpolation, is not involved:
   this script does a full mode solve at all 11 wavelengths.)
9. **Regenerate the golden references.** `tests/golden/directional_coupler_hybrid.json`
   was last regenerated by job 11460127 on 2026-07-24 19:35 — with **pre-fix** code,
   so its stored through/cross values carry the mislabelled wavelength axis and the
   `--runslow` golden test will now disagree. It is currently uncommitted-modified in
   the working tree. Re-run `batch/run_golden_directional_coupler_hybrid.batch` after
   items 1–2 are signed off, and commit the result with a note that the change is the
   axis fix, not a regression. Do not quietly overwrite it — a golden that moves
   without explanation is worse than a failing test.
10. **The same bug is in the ring resonator examples.** Found while auditing the cache
    for item 3, and *not* fixed by items 1–2: `ring_resonator_hybrid.py` (lines 385,
    469) and `ring_resonator_hybrid_cascaded.py` (lines 392, 463) each compute
    `fwidth = (1/WL_MIN - 1/WL_MAX) * 0.5` and hand it to
    `sim.add_mode_monitor(fcen, fwidth, N_WL_FDTD)` — character for character the
    defect item 1 fixed, in raw meep rather than through `FDTDSolver`, so the solver
    fix does not reach it. Every ring FDTD result carries the same reversed, half-span
    wavelength axis.

    Worse, their 6 cache entries are still *reachable*: these scripts build their own
    flat config dicts via `_sim_params()` (`ring_resonator_hybrid.py:297`), so nothing
    in items 1–3 changed their hash and a re-run will silently reuse pre-fix data. They
    were left on disk deliberately rather than deleted — deleting would only force
    re-simulation of equally wrong results. Fix: port both scripts to
    `monitor_freq_grid()` / `source_fwidth()`, add the resolved grid to `_sim_params()`
    so the hash strands the old entries automatically (the mechanism item 2 used), then
    re-run. Check `notes.md` before assuming these are the last two — this is the second
    place the same line has appeared.

## Verification

**Items 1–2 have a dedicated harness** (temporary — delete with the items):

    BEEF/examples/verification/python/verify_items_1_2.py
    batch/verification/run_verify_items_1_2.batch   →  make verify_items_1_2

Nine checks, output to `example_output/verification/items_1_2/`
(`verification.json` with per-check PASS/FAIL plus full provenance,
`verification.npz`, and four plots). Exits non-zero on any failure, so the job's exit
status means something. Checks 1–6 are static and need no meep run:

    cd BEEF/examples && OMP_NUM_THREADS=1 python verification/python/verify_items_1_2.py --static-only

### Status: 9/9 pass — job 11477911, 2026-07-25

48 ranks, 56 min, exit 0; meep 1.33.0, beef `c243a46+dirty`, uncached.

| check | result |
|---|---|
| 2. grid matches requested λ | 2.2e-13 nm error (pre-fix: 75 nm) |
| 4. source covers the band | weakest λ at 0.269 of peak (floor 0.10) |
| 7. T_cross rises with λ | dT/dλ = **+3.378 /µm**, 0.276 → 0.613 (pre-fix: backwards) |
| 8. pre-fix data on its true axis | **max \|ΔT_cross\| = 0.00025, rms 0.00014** (tol 0.05) |
| 9. power conserved | T_total 0.9992–0.9997 (tol 0.02); `check_passive()` OK |

Check 8 is the decisive one and it passed by a factor of 200, not by a hair — worth
stating, because a crude straight line would have scored 0.043 against that 0.05
tolerance. Two independent 3-D FDTD runs, one pre-fix and one post-fix with a
different source bandwidth (fwidth 0.0208 → 0.0833), agree to 2.5e-4 in cross
transmission once the old run is placed on the wavelengths meep really sampled. This
is the claim item 1 makes, confirmed at the level of the solver's own convergence:
**the FDTD physics was never wrong, only the labels were.** `dispersion_before_after.png`
shows it directly — the re-plotted pre-fix points sit on the post-fix curve, and the
as-labelled pre-fix curve crosses it at exactly one point, λ=1.55, the single index
both axes agreed on.

Check 8 rebuilds the archived pre-fix config exactly
(`archive/example_output_20260724/.../directional_coupler_fdtd_resolution_sweep/resolution_sweep.npz`,
res-40 row — same device, `coupling_length` 13.5 µm, same wavelengths). The run is
deliberately **uncached** — a verification must not be able to pass on a stored
result, and `save_run()` is not multi-rank safe (item 4).

**Do not misread check 5.** It reports the pre-fix source leaving the band edges at
7.3e-10 of peak, which reads like the old run was noise. It was not, and check 8's
2e-4 agreement is the proof. Check 5 measures the *old* fwidth against the *full*
1.50–1.60 µm request — a counterfactual. Pre-fix, source and monitor shared that narrow
fwidth, so the pulse covered everything meep actually sampled. What check 5 shows is
that item 2 was required *because* item 1 widened the grid: fixing the grid alone would
have starved the new edges. Both were needed; neither alone was sufficient.

Keep the harness until item 7 closes — items 5 and 7 re-run FDTD sweeps and this is
what re-checks them. Then delete the script, the batch file, and the Makefile's
`verify_items_1_2` / `clean_verification` targets.

Remaining items:

- Items 1–2 together: ✅ confirmed by check 7 — the full-device FDTD 50:50 length now
  *decreases* with wavelength (T_cross rises). Any future result where it increases
  means the axis has flipped again.
- Item 4: ✅ done — see Done below. The multi-node check has its own script, since
  pytest cannot reproduce 48 ranks on 2 nodes sharing one cache:

      srun -A <account> -q inferno -N 2 --ntasks-per-node=24 -t 10 --mem-per-cpu=2G \
           --mpi=pmi2 make verify_cache_multirank

  Two nodes minimum and a save_dir on shared storage, or it proves nothing.
- Item 6: the new test fails at ~11% before the fix, passes after.
- Item 7: `L_50` vs `BODY_STUB` curves and flattens; the hybrid's λ=1.55 value lands
  on the re-run full-device FDTD reference.

## Done

### 1. Monitor frequency grid — 2026-07-25

Every DFT monitor in `beef/solvers/fdtd.py` now takes an explicit frequency list
instead of meep's `(fcen, df, nfreq)` triple, and results are reordered back onto the
caller's wavelength array on unpack.

**Changed**

- **`beef/solvers/fdtd.py` — new `monitor_freq_grid(wavelengths)` helper** (module
  level, next to `section_fdtd_layout`). Returns `(mon_freqs, wl_order)`:
  `mon_freqs` is `1/wavelengths` sorted strictly ascending, which is what meep
  requires; `wl_order` is the inverse permutation that lifts a monitor-ordered result
  back onto `wavelengths` as the caller gave it. Works for unsorted input, not just
  ascending λ. Docstring records the defect so the triple doesn't come back.
- **`solve()`** — `nfreq` dropped; `mon_freqs, wl_order = monitor_freq_grid(...)`.
  Both `add_mode_monitor` calls (input normalization + the per-output-waveguide list)
  now pass `mon_freqs`. The normalized amplitude is reordered once at the end:
  `alpha = (coeff.alpha[0, :, 0] / src_alpha)[wl_order]`.
- **`solve_section()`** — same three changes: `nfreq` dropped, both `norm_flux` /
  `mea_flux` monitors take `mon_freqs`, and
  `coeff = (c.alpha[0, :, alpha_idx] / src_alpha)[wl_order]`. Ratio taken in monitor
  order first, reordered once — `src_alpha` and the measured amplitude are both
  monitor-ordered, so they stay aligned.
- `fwidth` is still computed at both sites but now feeds **only** the
  `GaussianSource`. Its value is wrong for that job — that is item 2.

**Added** — `BEEF/tests/test_fdtd_freq_grid.py`, 5 tests:

- monitor frequencies ascending and exactly `1/wavelengths`, with the endpoints
  asserted against the old grid's 1.5254–1.5754 µm band;
- `wl_order` round-trips a monitor-ordered result back onto the caller's array;
- the ascending-λ case is a plain flip, stated explicitly so a refactor can't drop it;
- unsorted input round-trips too;
- end-to-end against meep: `mp.get_flux_freqs()` returns exactly the requested
  wavelengths. Needs a source in the sim (meep allocates field components from
  sources, and rejects a monitor on a fieldless simulation) but no FDTD run —
  the whole file runs in ~2 s.

**Verified** — `pytest tests/test_fdtd_freq_grid.py` 5 passed; full fast suite
`pytest tests/` 130 passed, 6 skipped (the `--runslow` golden tests), no regressions.

**Confirmed physically** — job 11477911 (2026-07-25), checks 7–8: cross transmission
now rises with wavelength, and the archived pre-fix run re-plotted on its true axis
lands on the post-fix curve to `max |ΔT_cross| = 0.00025`. The relabelling account of
this bug is correct; the old data was recoverable, not garbage.

**Not addressed here** — every cached FDTD result predates this fix and is still
mislabelled; the fingerprint does not record the frequency grid, so those caches will
be reused silently. That is item 3, and it blocks trusting any re-run.
*(Superseded by item 2, which changed the fingerprint and stranded them all.)*

### 2. Source bandwidth — 2026-07-25

The Gaussian pulse now covers the band the monitors actually sample, and the width is
a real solver parameter rather than a magic 0.5.

**Changed**

- **`beef/solvers/fdtd.py` — new `source_fwidth(wavelengths, factor)`.** `factor`
  multiplies the **full** span `1/λmin − 1/λmax`; the old code used `span * 0.5`.
  Meep's Gaussian falls as `exp(-2π²(Δf/fwidth)²)`, so measured against a 1.50–1.60 µm
  band the edges sit at 0.005 of peak at factor 1.0, 0.097 at 1.5, and **0.27 at the
  new default 2.0**. The ratio is what governs it, so this holds for any band or point
  count — confirmed identical for an O-band sweep and for a 51-point grid.
- **Single-wavelength requests no longer produce `fwidth = 0`.** Span is zero there,
  and meep's pulse width is `1/fwidth`, so the old expression yielded an unusable
  source. Falls back to `factor * 0.1 * fcen`, i.e. meep's usual `0.2 * fcen` at the
  default. Latent bug, found while here.
- **New `check_source_coverage(fcen, fwidth, wavelengths, min_fraction=0.1)`**, called
  by both `solve()` and `solve_section()`. Evaluates the real
  `mp.GaussianSource.fourier_transform` at each requested wavelength and logs a
  warning naming the worst point if any falls below 10% of peak. A starved point still
  yields a well-conditioned transmission — the source spectrum cancels in the
  measured/normalization ratio — but it burns dynamic range, and `decay_tol` is
  applied to DFT fields orders of magnitude smaller there.
- **`FDTDSolver(source_bandwidth_factor=2.0)`** is a real constructor parameter with a
  docstring entry, not a module constant — `notes.md` Solver #9's sub-note complains
  precisely about physics knobs that live in defaults nothing can see.
- **`beef/cache/config.py`** — `source_bandwidth_factor` added to `_METHOD_ATTRS` for
  both `("FDTDSolver", "solve")` and `("FDTDSolver", "solve_section")`. That list is
  an explicit allowlist, so a new attribute is *not* picked up automatically; without
  this the knob would move results while leaving the hash unchanged.

**Added** — 5 more tests in `BEEF/tests/test_fdtd_freq_grid.py`: the factor scales the
full span; a single wavelength still gets a usable pulse; the shipped default clears
the coverage floor with margin (>0.25) and peaks at fcen; the pre-fix `span * 0.5`
leaves the edges below 1e-6 **and now warns** (regression lock); and the knob is
reachable from the constructor.

**Verified** — `pytest tests/test_fdtd_freq_grid.py` 10 passed; full fast suite
`pytest tests/` 135 passed, 6 skipped, no regressions. Fingerprint check:
`source_bandwidth_factor` appears in the built config and two factors produce
different configs. Cache audit of `BEEF/examples/example_sim_cache`: 119 cached runs,
54 of them `FDTDSolver`, **0** carrying the new key — so every stale FDTD result now
misses instead of being silently reused. The 65 EME entries are unaffected and remain
valid.

**Confirmed physically** — job 11477911 (2026-07-25), check 4: the shipped default puts
the weakest requested wavelength at 0.269 of peak in a real 3-D solve, and check 9 gives
power conservation 0.9992–0.9997 with `check_passive()` clean, so the wider pulse cost
nothing in accuracy. See the check-5 caveat in Verification — the wider bandwidth was
needed *because* item 1 widened the grid, not because the pre-fix run was invalid.

**Consequence** — items 1+2 together made job 11477911 the first trustworthy FDTD output
in this project, and it is now the reference for the coupling-vs-wavelength shape at
`coupling_length` 13.5 µm. Numbers in the diagnostic sections *above* still come from
the old, mislabelled data and are kept as the record of the bug; item 5 replaces the
ground-truth length sweep.

### 3. Cache fingerprint completeness — 2026-07-25

Everything that changes an FDTD result is now in the hash, and the unreachable entries
are off disk.

**(a) Subpixel smoothing is a real parameter now.** `eps_averaging`, `subpixel_tol` and
`subpixel_maxeval` were on meep's implicit defaults — never passed, never fingerprinted
— while moving `n_eff` at the second decimal (`notes.md` Solver #9). They are
`FDTDSolver` constructor arguments, passed to **both** `mp.Simulation` calls
(`fdtd.py:501`, `:798`), and in `_METHOD_ATTRS` for `solve` and `solve_section`.
Defaults are meep's own, so no existing result changes meaning; a test asserts that,
since a wrong default here would silently alter every run.

**(c) The monitor grid does not need its own schema bump — it needs its resolved value
hashed.** New `_resolved_spectral_config()` puts `resolved_monitor_freqs` and
`resolved_source_fwidth` into the `call_kwargs` fingerprint at all three `FDTDSolver`
config sites. Both are derived from `wavelengths`, which is already hashed, so the
fingerprint was *correct* today — but that is exactly the state item 1's bug lived in
for months: the grid was derived wrongly while every cache entry looked identical to a
right one. Now a change to the derivation changes the hash even when the inputs do not.

This is the pattern `build_run_config()`'s own docstring prescribes for
resolved-vs-ingredient values (Solver recommendation #4), already used for
`resolved_layout`. Chosen over bumping `_CACHE_SCHEMA_VERSION`, which would also have
stranded the 20 valid EME entries the bug never touched, and which depends on a human
remembering to bump it.

**(b) Dead entries deleted.** 54 unreachable `FDTDSolver` folders (4.4 MB) removed from
`BEEF/examples/example_sim_cache`; 120 entries → 66, then 65 after the audit below.
"Unreachable" was tested, not assumed: an entry whose `config.json` lacks a key current
code always emits can never be matched by `lookup()`.

`index.jsonl` was also rebuilt, 5760 lines → 65. It had grown 48× because every rank
appended (item 4), and 2640 of its lines pointed at folders just deleted. Regenerated
from the surviving `config.json` files through the same `_manifest_entry()` the writer
uses, so it is now exactly one line per entry, 0 orphans.

**Added** — 15 tests in `tests/test_cache_config.py`: the three subpixel knobs
discriminating on both methods, their defaults matching meep's, and a
`TestResolvedSpectralConfig` class covering the grid being ascending in frequency and
equal to the requested wavelengths, `fwidth` tracking the factor, and — the one that
matters — **a reversed grid derivation changing the hash** while an identical call
still shares one.

**Verified** — full fast suite **175 passed, 6 skipped**. Post-deletion audit: 65
entries, **0 unreachable**, 0 orphaned index lines.

**Found while doing this, not fixed here** — the same monitor-grid bug is live in both
ring resonator example scripts, and their cache entries are still reachable. Raised as
item 10 rather than folded in silently; it is a different code path (raw meep, its own
config builder) and deserves its own change.


### 4. Cache write race — 2026-07-25

Only one rank writes the cache now, and a damaged cache entry can no longer kill the
job that reads it.

**The shape of the bug.** Every rank of an `srun -n 48` job runs the whole solve
redundantly, so every rank also reached every write path in `store.py` — 48 processes
across 2 nodes writing one directory on shared storage. `index.jsonl` is the receipt:
**5760 lines for 120 saves, exactly 48× over-counted.**

**Changed — `beef/cache/store.py`**

- **New `_is_master()`.** Deliberately returns a bool, not a rank: nothing here needs
  to know *which* rank it is, and inventing a number later callers could misuse is how
  this kind of bug starts. Sources in order of authority, importing nothing not already
  loaded (importing mpi4py into a serial process would call `MPI_Init` behind meep's
  back): meep's `am_really_master()` → mpi4py → the launcher's environment.
  `am_really_master()`, **not** `am_master()` — the latter is per-group under
  `divide_parallel_processes` and would elect one writer per group.
- **The environment fallback is load-bearing, not a nicety.** An EME-only run imports
  neither meep nor mpi4py (meow/MPB don't pull meep in) yet is still launched under
  `srun -n 48`, so a meep-only check would have called all 48 ranks master and left the
  bug fully intact on that path. `SLURM_PROCID` / the OMPI/PMI variables answer without
  importing anything.
- **`save_run()` and `update_family_pointer()` write only on master**, then `_barrier()`.
  The barrier is the important half: without it a later `lookup()` could hit on one rank
  and miss on another, and ranks disagreeing about whether to simulate leaves meep's
  collective calls unmatched — the job hangs until it burns its wall clock, which is
  what actually consumed the 5 h. Both are collective, called from all ranks or none;
  every call site is guarded by conditions derived from the solver config, never from
  per-rank data, so that holds.
- **`save_run()` still returns the folder path on every rank**, since callers feed it
  straight into `update_family_pointer()`.
- **New `_read_family_pointer()`.** A pointer is a reuse hint, never data — the result
  it names is also reachable through the plain `lookup()` path — so an empty, truncated
  or schema-drifted one is now a cache miss with a warning, exactly like `lookup()`'s
  config-mismatch guard, instead of an exception. This is the actual crash site: tasks 4
  and 5 raised out of `json.load` here and took a 48-rank job down. Used by both
  `lookup_superset()` and `update_family_pointer()`, so damaged pointers heal on the
  next write.
- **`_atomic_write_json()` temp names carry host as well as PID.** PIDs are allocated
  per host, so a PID-only name is not unique on a filesystem spanning nodes.
- **Docstring corrected.** The claim that concurrent writers cost only "a missed reuse
  opportunity, never an incorrect one" was true across separate processes (a SLURM
  array) but false for ranks within one job. Now says which case it means.

**Added** — `BEEF/tests/test_cache_multirank.py`, 25 tests: rank detection through all
four sources including the `am_master`/`am_really_master` distinction; nothing written
off master while the folder path still returns correctly; `index.jsonl` getting one line
per job rather than per rank; host-qualified temp names; and five flavours of damaged
pointer (`""`, truncated JSON, whitespace, `null`, `[]`) reading as a miss, healing on
the next save, without turning every update into a clobber.

**Verified** — full fast suite `pytest tests/` **160 passed, 6 skipped** (135 baseline +
25). Then the part unit tests cannot reach, via the new
`examples/verification/python/verify_cache_multirank.py` at **48 ranks on 2 nodes**,
production shape, meep imported, no simulation, ~1 s:

    ranks           48 on 2 hosts
    index.jsonl     5 lines (expected 5; pre-fix would be 240)
    stray tmp files 0
    pointer dims    [6, 1]

Audit of the live cache found **no lasting damage**: 0 corrupt pointers, 0 stray temp
files, all 5760 `index.jsonl` lines parseable.

**Correction to the original diagnosis.** This document previously stated as fact that
two ranks on different nodes drew the same PID. That is *plausible but unproven* — the
48-rank probe found no PID reused across hosts, and concurrent `os.replace`/read against
shared storage explains an empty file just as well. The rank gate removes every version
of the race, so the fix does not depend on which mechanism it was; the docstrings now
say so rather than asserting the PID story.

**Re-scoped, not skipped.** Item 4 originally included re-running the two dead stub
tasks. Those seven surviving `S*.npz` files are dated 2026-07-24 — **pre-fix**, since
items 1–2 landed 2026-07-25 11:56 — so re-running only tasks 4 and 5 would leave one
folder holding a silent mixture of pre- and post-fix data. A visible gap beats an
invisible mixture, so the re-run moved to item 7, which now covers all nine tasks after
deleting the folder.
