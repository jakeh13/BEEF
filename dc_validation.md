# Directional coupler validation — to-do

Opened 2026-07-25, from investigating why `directional_coupler_hybrid_stub_sweep.py`
never converges. Root cause found; this is the work to fix it and re-establish trust
in the directional coupler results. Nothing has been changed yet.

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
- Not usable as a spectral reference: the same axis bug corrupts its other ten
  wavelength points.
- Provenance is unrecorded. The npz stores only `coupling_lengths`, `wavelengths`,
  `T_through`, `T_cross` — no config, no solver fingerprint, no geometry, and the
  whole-device `solve()` path had no `save_dir` then (the archived cache holds only
  `solve_section`/`solve_section_length_sweep` entries). Its parameters match the
  current script's constants, but nothing in the file proves it.
- Coarse: 2 µm length steps, so 15.29 µm is interpolated across a 2 µm bracket
  (~±0.1 µm), against a stub sweep resolving 0.06 µm steps.

Verdict: use 15.29 µm as a provisional sanity check now; re-run the sweep after item 1
on a finer grid before treating any number as truth (item 5).

## To-do

Numbers are stable IDs — completed items move to Done and keep their number, so
cross-references stay valid.

**Order: run the items-1–2 verification first** (`sbatch
batch/verification/run_verify_items_1_2.batch`, ~1 h). Items 5, 6 and 7 all assume
items 1–2 really fixed the FDTD; if check 8 fails, that work is misdirected. The run
is uncached, so item 4's race cannot touch it and item 3 does not gate it. Then
**4** (+ the rest of **3**) while it runs — pure code, no queue, and it is what
otherwise wastes wall-clock on the big array jobs in 5 and 7 — then **5** → **9** →
**6** → **7**. Item 6 comes before 7 because 7's interpretation depends on whether the
residual ~11% gap is understood. That gap was measured at λ=1.55, the one index both
engines labelled correctly before the fix, so it should survive verification intact.

3. **Finish the cache story.** Largely handled as a side effect of item 2 — adding
   `source_bandwidth_factor` to `_METHOD_ATTRS` changed the fingerprint, so all 54
   cached `FDTDSolver` runs now miss and cannot be silently reused (verified). The
   65 EME entries are untouched and stay valid; the bug never reached them. What is
   left: (a) expose and fingerprint `eps_averaging` / `subpixel_tol` /
   `subpixel_maxeval`, still on meep's defaults and still invisible to the hash —
   `notes.md` Solver #9's sub-note measures them moving `neff` at the second decimal;
   (b) delete the 54 dead FDTD folders in `BEEF/examples/example_sim_cache`, now
   unreachable but still on disk; (c) decide whether the monitor frequency grid
   deserves its own fingerprint entry — it is derived from `wavelengths`, which *is*
   hashed, so it is correct today, but a future change to the derivation would be
   silent again, which is exactly how item 1 survived this long.
4. **Fix the cache write race and re-run the two dead stubs.** Tasks 4 and 5 (`S100`,
   `S125`) died at `beef/cache/store.py:283`, `json.load` on an empty file.
   `_atomic_write_json` (`store.py:189`) names its temp file
   `f"{path}.tmp.{os.getpid()}"`; under `srun --mpi=pmi2` with 48 ranks across 2 nodes
   every rank executes the write and PIDs collide across nodes on the shared
   filesystem, so one rank truncates the temp file another is renaming. After the
   exception the job hung and burned the full 5 h wall clock. Gate `save_run` /
   `update_family_pointer` on rank 0, or make the temp name unique per host+rank. The
   docstring's claim that concurrent writers can only cost "a missed reuse
   opportunity, never an incorrect one" is not true for the multi-rank case — update
   it.
5. **Re-establish ground truth.** Re-run the full-device FDTD length sweep after 1–3,
   on a finer grid than 2 µm near the crossing, and save the solver config into the
   npz so provenance is self-describing next time. This is what the stub sweep's
   reference line reads from —
   `example_output/directional_coupler/directional_coupler_fdtd_sweep/length_sweep.npz`
   does not currently exist in the live tree at all.
6. **Explain the residual ~11% rate gap at λ=1.55.** Not covered by 1–2. Add a direct
   engine-agreement test alongside `tests/test_hybrid.py`: take a pure straight
   coupled-pair region (no bends), call `FDTDSolver.solve_section()` at two lengths,
   extract Δβ_FDTD from the even/odd phase difference, compare against
   `EMESolver.compute_modes()`'s neff split on the same cross-section, assert to a
   stated tolerance. Expect it to fail at ~11% today. This is the engine-consistency
   check `notes.md` Solver #9 asks for and the cleanest way to isolate the effect —
   mesh discretization is already ruled out by the resolution sweep.
7. **Re-run the stub sweep and judge it properly.** Success is *curvature* — successive
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

## Verification

**Items 1–2 have a dedicated harness** (temporary — delete with the items):

    BEEF/examples/verification/python/verify_items_1_2.py
    batch/verification/run_verify_items_1_2.batch   →  make verify_items_1_2

Eight checks, output to `example_output/verification/items_1_2/`
(`verification.json` with per-check PASS/FAIL, `verification.npz`, and three plots).
Exits non-zero on any failure, so the job's exit status means something. Checks 1–6
are static and need no meep run:

    cd BEEF/examples && OMP_NUM_THREADS=1 python verification/python/verify_items_1_2.py --static-only

**Status: 6/6 static checks pass** (2026-07-25). Monitor grid reproduces the requested
wavelengths to 2.2e-13 nm against the pre-fix 75 nm error; source floor 0.269 of peak
against the pre-fix 7.3e-10. Checks 7–8 (the FDTD run) have **not** been run yet.

Checks 7–8 are the physical ones, and 8 is the strong one: it rebuilds the exact
config of the archived pre-fix run
(`archive/example_output_20260724/.../directional_coupler_fdtd_resolution_sweep/resolution_sweep.npz`,
res-40 row — same device, `coupling_length` 13.5 µm, same wavelengths) and asserts
that re-plotting that old row on the wavelengths meep *actually* sampled
(1.5754→1.5254 µm, reversed) lands on the new run's curve. That is the claim item 1
makes: the old data was only mislabelled, same physics, wrong axis. The run is
deliberately **uncached** — a verification must not be able to pass on a stored
result, and `save_run()` is not multi-rank safe (item 4).

- Items 1–2 together: the full-device FDTD 50:50 length must *decrease* with
  wavelength. Any result where it increases means the axis is still flipped.
- Item 4: `sbatch batch/sweep/run_directional_coupler_hybrid_stub_sweep.batch --array=4,5`
  completes and `data/S100.npz`, `data/S125.npz` appear.
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

**Consequence** — items 1+2 together mean the next real run is the first trustworthy
FDTD output in this project. Nothing has been re-run yet; every number in the sections
above still comes from the old, mislabelled data.
