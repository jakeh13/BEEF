# Active SLURM jobs — snapshot

**★ TWENTIETH SNAPSHOT — 2026-08-03, verified against `squeue -u jhiesener3`
(empty) and `sacct -u jhiesener3 -S 2026-08-02`. This is the current one.**

## NOTHING RUNNING. Queue is empty.

| job | state | what |
|---|---|---|
| **11634859** | ✅ **COMPLETED 3:56:16**, exit `0:0`, ended 2026-08-02T23:10:59 | `beef_ring_xc_clampfix`. Ran to the end well inside `-t 480`. **Right section and composed cross-check both landed** — see below |
| **11634900_0** | ✅ COMPLETED 58:06 | `ring_resonator_hybrid`, 3-D. Mean FSR **12.975 nm (0.8%)** |
| **11634900_1** | ✅ COMPLETED 44:18 | `ring_resonator_hybrid_cascaded`, 3-D. Mean FSR 13.350 nm (3.7%) |

### ✅ 11634859 finished — the whole cross-check passes, not just the left section

The nineteenth snapshot recorded this job's left-section numbers while it was
still at 1:57:04. It has since completed. New material:

```
right bus<-iface:  summed |T|^2 over seam modes = 0.9997 0.9992 0.9970 0.9953 0.9933
right ring<-iface: summed |T|^2 over seam modes = 0.9989 0.9985 0.9989 0.9986 0.9980
composed ports: ('bus_L', 'ring_L', 'bus_R', 'ring_R')
bus-side power sum |t|^2+|k1|^2: min 1.0050 max 1.0067   (lossless coupler -> ~1)
reciprocity max ||k1|-|k2||: 0.0021   (-> 0)
physics checks PASS (power sum ~1, reciprocity holds)
worst |E_legacy - E_network| = 1.735e-18
L_outer = 7.0000 um
VERDICT: routes AGREE to machine precision
```

**The right section is clean too** (0.993–1.000, no clamp pathology), so the
fix holds on both sides at production config. Power sum and reciprocity pass on
the *composed* 4-port. Per the standing entry, the `AGREE to machine precision`
verdict is the weakest line here — it printed identically pre-fix with ~88% of
the power gone.

✅ **`crosscheck.npz` is restored.** This job re-wrote
`example_output/verification/ring_network_crosscheck/crosscheck.npz` at 23:10,
replacing the array the login-node smoke runs destroyed. Contents are now
post-clamp-fix and valid. The `OUT`-redirect warning still stands for future
login-node runs.

⚠ **This job printed `L_outer = 7.0000 um`, which corrects gap 3.** Gap 3's
convergence script hardcodes `L_OUTER = 2*2.0 + 10.0` = 14.0 µm on a comment
asserting `L_C=10.0` is "the shipped scale"; the shipped `L_C` is **4.0**, and
both device scripts compute `L_outer = 2*(L_SIDE - S_MON) + L_FAR` = **7.0 µm**.
Gap 3's lever arm is exactly **2× too large** and its Δφ column halves. Details
and the corrected table in `ring_resonator_validation.md` gap 3.

### The half-FSR route offset — candidate space collapsed by measurement

The nineteenth snapshot left this undiagnosed with two named candidates. **Both
are now eliminated, by source comparison plus the saved arrays — no new FDTD.**
`ring_data.npz` retains all coefficients as `complex128`, so this was
post-processing. Full write-up in `ring_resonator_validation.md`. Short version:

- `L_outer`, `phi`, `phi_bend_correction`, `BEND_D_IN/OUT`, `neff_fit` — the
  whole analytic phase chain is **byte-identical code fed identical constants**
  in both scripts, and `neff` fit residual prints identically (7.29e-06).
- `b_raw` is **identical to all printed digits in magnitude AND phase**
  (`|b|=0.998828`, `arg=-1.275855` rad, difference exactly `0.000000`).
- So the offset enters through the **coupler coefficients**, not the ring
  phase: `arg(t_ring)` differs by **0.883π**, `arg(t)` by 0.404π.

**Measured. Not yet attributed to a defect** — see the validation doc for what
that does and does not license.

---

**NINETEENTH SNAPSHOT — 2026-08-02 ~20:16, verified against `squeue -u
jhiesener3` and `sacct -X`. Superseded; kept unedited.**

| job | state | what |
|---|---|---|
| **11634859** | ⏳ RUNNING **1:57:04** @ 21:11 (`-t 480`) | `beef_ring_xc_clampfix`. **Left-section numbers in and they close item 2's question** — see below |
| **11634900_0** | ✅ COMPLETED **58:06** | `ring_resonator_hybrid`, **3-D**. **Mean FSR 12.975 nm, 0.8% from the analytic 12.87.** Reciprocity `\|k1\|=0.2530 vs \|k2\|=0.2530`, independently normalized |
| **11634900_1** | ✅ COMPLETED **44:18** @ 20:13 | `ring_resonator_hybrid_cascaded`, 3-D. Mean FSR 13.350 nm (3.7%) |

**Order item 1 is CLOSED** — first trustworthy resonance positions this device
has produced. **But the two routes disagree by ~half an FSR**, undiagnosed;
full write-up in `ring_resonator_validation.md`'s new
`## Order item 1 — the re-run` section. Do not quote one route's comb without
saying which.

### ✅ 11634859 — the clamp fix is confirmed at production config

```
left iface<-bus:  1.0056 1.0073 1.0084 1.0098 1.0117
left iface<-ring: 0.9903 0.9906 0.9907 0.9916 0.9920
left ring<-bus:   0.0000 0.0000 0.0000 0.0000 0.0000
left TOTAL out of bus = 1.0056 1.0073 1.0084 1.0098 1.0117
left deficit vs 1.0   = -0.0056 -0.0073 -0.0084 -0.0098 -0.0117
```

Against **11627401's +0.7065 … +0.9307**, at the same 3-D res-20 config. The
~80% deficit that item 2 chased for two days was
`RingResonatorCouplerSection._section_polygons()` drawing bare cladding where
the caller's stub padding should have continued the bus. Right section and
composition still pending.

**Runtime is itself consistent with the fix:** 1:01+ against pre-fix runs of
7:04 and 8:21. Pre-fix the source radiated straight into the PML and
`stop_when_dft_decayed` tripped early — fast *because* the power was
vanishing. Suggestive, not proof.

### ✅ 11634900_1 — cascaded ring, 3-D, first trustworthy resonance comb

```
bus-excited path : |t|^2+|k1|^2 = 0.9998
ring-excited path: |t_ring|^2+|k2|^2 = 0.9928
reciprocity: max|k1-k2| = 1.78e-16   (symmetric BY CONSTRUCTION -- weak evidence)
supermode neff fit residual: even=1.45e-05, odd=7.30e-07
Dense-grid assembly (2001 points): 8 resonance dip(s) found.
Measured mean FSR: 13.350 nm (predicted 12.87)
Results at 1.5500 um: |t|=0.9713 |k1|=0.2374 |k2|=0.2374 |t_ring|=0.9677 |b|=0.9988
```

FSR is **3.7%** off prediction, against **14%** in the 2-D `RING_TEST_MODE`
smoke — consistent with gap 3's measurement that res 20 sits ~75% of an FSR
from res 40.

⚠ **Q is grid-limited, and this run makes that visible rather than suspected.**
Every FWHM is an exact multiple of the 2001-point grid spacing (0.05 nm):

| FWHM (nm) | 0.1500 | 0.2000 | 0.2500 | 0.3000 |
|---|---|---|---|---|
| grid steps | 3 | 4 | 5 | 6 |

So the reported Q (5322–10193) is set by **`N_WL_DENSE`**, a knob independent
of both `FDTD_RESOLUTION` and `MODE_RESOLUTION`. In the 501-point smoke every
FWHM was pinned at exactly 0.4000 nm. **If Q is wanted, raise `N_WL_DENSE`** —
it is post-processing on an interpolated grid and costs no new FDTD.

---

**EIGHTEENTH SNAPSHOT — 2026-08-02 ~19:15, verified against `squeue -u
jhiesener3`. Superseded; kept unedited.**

| job | state | what |
|---|---|---|
| **11634859** | ⏳ PENDING @ 19:15 (`cpu-small`, reason `None` — not a billing hold) | `beef_ring_xc_clampfix` — `verify_ring_network_crosscheck.py` re-run after the `_section_polygons` clamp fix. 1 node × 24, `-t 480`, `--account=gts-sr108-ece` |

**What it is expected to answer.** Every bus-excited number in
`ring_resonator_validation.md` item 2 was launched from a source plane sitting
on bare cladding (see that file's `🛠 DEFECT FOUND AND FIXED` block). A
controlled login-node A/B at 2-D res 8 already puts `TOTAL out of bus` at
**0.2190 → 1.0087** across the fix alone. This job is the 3-D res-20
production confirmation of that, at the config the recorded numbers were taken
at. Expect `iface<-bus` ≈ 1 and the +0.71…+0.93 deficit to vanish.

**When it lands:** fold into `ring_resonator_validation.md` item 2, replacing
the 11627160/11627401 numbers rather than appending to them — those runs are
void, not merely superseded.

⚠ The script `np.savez`es to
`example_output/verification/ring_network_crosscheck/crosscheck.npz`
unconditionally. It currently holds **pre-fix-clamp 2-D res-8 A/B data** from a
login-node run, not a production result. This job overwrites it correctly.

---

**SEVENTEENTH SNAPSHOT — 2026-08-02 ~13:10, verified against `squeue -u
jhiesener3` (empty) and `sacct -X -S 2026-08-02`. Superseded; kept unedited.**

**THE QUEUE IS EMPTY. Nothing is running, nothing is pending.** Every job
this project submitted today has terminated.

| job | final state | answer |
|---|---|---|
| **11627401** | ✅ COMPLETED **8:21** @ 12:00 | `beef_ring_xcC_busbus` — the `bus<-bus` control. Reflection **0.1–2%**, forward **0.049–0.326**. Both surviving readings from item 2 stay alive; ~80% still unaccounted. Folded into `ring_resonator_validation.md` item 2 |
| **11627585** | ❌ **FAILED** (exit 1) **3:21** @ 12:11 | `beef_ring_xcD_launch` — the straight-bus launch control. Died in geometry, **before any FDTD physics**: `ValueError: Cut/stub at x=-3.5000 falls strictly inside the near-bend arc` |
| **11628207** | ❌ **FAILED** (exit 1) **3:20** @ 13:05 | `beef_ring_xcD_launch2` — same control after an attempted fix. **Failed identically**, `x=-6.5000`. The fix was wrong |

### ⛔ The launch control has never run, and the script as committed will fail a third time

**Measured 2026-08-02, not argued.** Both failures are the same defect, and
the correction applied between them did not remove it.

`RingResonatorCouplerSection._check_cut` (`beef/geometry/ring_coupler.py:168`)
rejects any cut **strictly inside** the near-bend arc, which for this device
spans `(stub_x, hc) = (-7.0, -2.0)`. The quantity it validates is not `x1`,
and it is not `x1 + stub_right` either:

```
FDTDSolver.solve_section()  beef/solvers/fdtd.py:1978
    prisms = geometry.section_prisms(
        x0, x1, si_med,
        stub_left=stub_left + self.dpml, stub_right=stub_right + self.dpml, ...)
```

**`solve_section` pads both stubs by `dpml` before the geometry ever sees
them** — deliberately, so the drawn waveguide continues through the full PML
instead of stopping at its inner edge. The cut actually tested is therefore

```
hi = X_STRAIGHT + stub_right + dpml = X_STRAIGHT + 1.0 + 0.5
```

Evaluated directly against the geometry class (conda `jmh`, login node, pure
geometry — no meep, no time-stepping):

| `X_STRAIGHT` | cut `hi` | guard | matches |
|---|---|---|---|
| −5.0 | −3.50 | ILLEGAL | job 11627585's error string exactly |
| −7.5 | −6.00 | ILLEGAL | |
| **−8.0 (committed default)** | **−6.50** | **ILLEGAL** | job 11628207's error string exactly |
| −8.5 | −7.00 | LEGAL | lands on `stub_x`, the boundary the guard allows |
| −9.0 | −7.50 | LEGAL | |

The reproduction is exact in both directions — the two failed jobs' error
strings are recovered from the geometry alone, so this is the whole mechanism
and not a coincidence.

**Why the second attempt failed.** The fix comment now standing at
`verify_ring_network_crosscheck.py:56-67` reasons: *"section_polygons()
validates x1 + stub_right … With common's stub_right=1.0, x1=-8.0 puts the
stub end exactly on stub_x=-7.0 — the boundary, which the guard allows."*
That arithmetic is right and its premise is wrong: it omits the `+ dpml`
term, so the real cut was −6.5, not −7.0. **The comment is still in the file
and still wrong**, and `X_STRAIGHT` still defaults to −8.0. Resubmitting the
batch script unchanged buys a third identical 3-minute failure.

**Cost so far:** 2 × ~3:20 wall on 24 cores. Cheap in absolute terms, but both
failures are pure geometry — reachable in under a second on a login node
without allocating a node at all.

**The open design call, which is why this is not patched here.** `X_STRAIGHT
= -8.5` is the largest legal value, and it leaves `[-9.0, -8.5]` — **0.5 µm**
of straight bus. That may be too short to be the control item 2 asked for
(the question is whether a clipped launch *sheds over a propagation
distance*). Dropping `stub_right` to 0.0 for this one call instead allows
`X_STRAIGHT ≤ -7.5` and a 1.5 µm span, but changes the port padding away from
what every other measurement in the script used. **Pick one deliberately** —
this is exactly the "sweep the wrong thing" trap the coupler's item 11 and the
MMI's gap 2 both recorded.

---

**SIXTEENTH SNAPSHOT — 2026-08-02 ~11:40, verified against `squeue -u
jhiesener3` and `sacct -X`.** Superseded by the seventeenth above; its two
live jobs have both since landed and 11627160 is unchanged. The fifteenth is
kept below it, unedited, for the two lessons it carries.

**Both jobs the fifteenth snapshot listed as live have landed, successfully.**

| job | final state | answer |
|---|---|---|
| **11609658** | ✅ COMPLETED **10:25:31** @ 01:42 (of 12 h) | `beef_dc_hybrid_dispersive`. Dispersive L₅₀ **15.089958 µm** (grid-snapped 15.150754, 0.299 µm grid). The 12 h budget was right — the 3 h that killed 11608703 was never going to be enough |
| **11609659** | ✅ COMPLETED **1:37:19** @ 03:20 | `beef_dc_fdtd_at_hybrid_l50`, `afterok` fired normally. At L=15.089958 the **band-mean** split ratio is **0.500675** |
| **11627160** | ✅ COMPLETED **7:04** @ 11:32 | `beef_ring_xcB_ringbus` — the ring's `ring<-bus` measurement. **Hypothesis refuted:** 0.0000–0.0004. See `ring_resonator_validation.md` item 2 |

### The coupler result, because it is easy to misread

The log's headline line reports **λ=1.550 only** (ratio 0.4919). The constant
`COUPLING_LENGTH` is defined on the **band mean**, which is 0.500675 — 0.13%
off a true 50:50. Reading the centre-λ line as the answer inverts the
conclusion. Computed from `smatrix_fdtd.pkl` directly.

| | frozen | dispersive |
|---|---|---|
| hybrid band-mean L₅₀ | 15.150754 | **15.089958** |
| FDTD band-mean L₅₀ | 15.129 | **~15.074** |

Dispersion moves **both** engines down ~0.05 µm and they stay ~0.1% apart.
**The hybrid-vs-FDTD gap survives the material model** — it is not a dispersion
artifact. This closes the `COUPLING_LENGTH` contradiction (2026-08-01 context
item 3): **15.13 was correct all along**, and the defect was the stale comment
block above it claiming the constant was pinned to the hybrid's 15.15.
Corrected 2026-08-02; value untouched. The ~15.074 figure rests on one length
plus the *hybrid's* slope, so it is **not** adopted as a constant — regenerate
item 5's ground-truth sweep dispersively first.

### ⚠ Not tracked here, and it should have been

**`11609123_7` `elevator_corner_round`** — running since 2026-08-01 13:21,
22 h+ elapsed of a 48 h limit, res=100, the last of an 8-point sweep. It lives
in `/storage/scratch1/8/jhiesener3/esc/`, not in BEEF, which is why the
2026-08-01 23:42 context save asserted "nothing else is running." **This file
is titled "Active SLURM jobs", not "active BEEF jobs".** It is also on
**`gts-sr108-gedc_main`** — the unfunded account, one of the 33 scripts still
pointing at it.

---

**★ FIFTEENTH SNAPSHOT — 2026-08-01 ~18:10, verified against `squeue -u
jhiesener3` and `sacct -X`.** Superseded by the sixteenth above.

**Every BEEF job the fourteenth snapshot tracked has landed.** All results are
folded into the owning device docs; this snapshot records the states and the
two things worth carrying forward.

| job | final state | answer |
|---|---|---|
| **11608875** `[0-3]` + **11609686** | ✅ all COMPLETED (55:15 / 41:22 / 49:33 / **2:23:09**), merge 25:45 @ 15:34 | **P1 CLOSED.** 1×2 gap 4, first ever: mean ratio **0.9946**, 0.54% monotone bias, worst residual 3.01e-06 rad |
| **11609215_1 / _2** + **11609687** | ✅ 13:44 / 2:00:01, comparison 0:12 @ 15:27 | 1×2 gap 5 dispersive: **9.11% → 10.18%** — marginally *worse*, mean offset unmoved |
| **11609211_1** + **11609680** + **11609681** | ✅ 1:51:57 / **3:50:41** / 0:15 @ 17:59 | **the P4 test → mechanism REFUTED.** 2×2 gap 5 dispersive **0.82% → 0.76%**, sign flip survives |
| **11608690** | ✅ COMPLETED 1:28:14 | **P2 CLOSED** — first calibrated hybrid to survive since the `beta_from_phases` fix |
| **11608692** | ✅ COMPLETED 00:07, **number contaminated** | gap 5 calibrated — valid number (7.30%) recovered offline; see the race below |
| **11608941** / **11608940** | ✅ 1:03:21 / 07:13 | gaps 2 and 3 dispersive — **byte-identical**, as predicted |
| **11608703** | ❌ **TIMEOUT 3:00:29** @ 15:17 | coupler L₅₀ re-verification — **did not finish**, see below |

~~**Two live jobs remain**~~ — **both landed 2026-08-02, see the sixteenth
snapshot at the top.** Left here as written:

| job | state | answer |
|---|---|---|
| **11609658** | ~~🟢 RUNNING 2:55 of **12:00**~~ ✅ 10:25:31 | `beef_dc_hybrid_dispersive` — the 11608703 resubmission, dispersive and with 4× the wall clock |
| **11609659** | ~~⏸ PENDING `afterok:11609658`~~ ✅ 1:37:19 | `beef_dc_fdtd_at_hybrid_l50`, 5 h |

### ❌ 11608703 TIMED OUT — two separate lessons

`beef_dc_hybrid_item14` hit its **3 h** limit at 15:17 with no output written.
Two things went wrong, and only one of them was the wall clock:

1. **Budget.** 3 h was sized against frozen-era coupler runs. The dispersive
   1×2 came in at **1.53×** its frozen cost (11609215_1) and the dispersive
   gap-4 lengths at ~2.4×; the coupler is the longest device in the repo. The
   resubmission carries **12 h**.
2. **⚠ It ran under `gts-sr108-gedc_main`, not `gts-sr108-ece`.** That is the
   account whose balance problems already blocked jobs on 2026-07-30
   (`AssocGrpBillingMinutes`). Everything submitted after it — 11609658,
   11609659, 11609686 — is explicitly on **`gts-sr108-ece`**. Check
   `#SBATCH -A` when reviving an older batch script; the timeout made this one
   visible, but a script sitting on the wrong account is silent until it is not.

**The P4 test was rescued from the same failure mode.** 11609211_2 was
cancelled 46 min in on a live cost projection rather than left to time out,
and its resubmission 11609680 took **3:50:41** — well past the original
budget, so the cancellation was correct. Two jobs hit the same
dispersive-cost-underestimate this session; one was caught in flight and one
was not.

---

## Fourteenth snapshot — 2026-08-01 ~13:35 (superseded, kept for its reasoning)

It adds the four jobs
submitted after the thirteenth was written plus the two diagnostics it omitted;
the thirteenth's *reasoning* about the dispersion generational break still
stands and is not repeated.

**Ten BEEF jobs are live.** Nothing has landed since 11608940. Grouped by what
they answer:

| job | state at snapshot | answers |
|---|---|---|
| **11608875** `[0-3]` + **11608876** | `_0` ✅ 55:15, `_1` ✅ 41:22, `_2` ✅ 49:33; `_3` 🟢 running; merge ⏸ `afterok` | **gap 4 / check 1 / P1** re-measured dispersive, fully cache-cold |
| **11608690** | ✅ **COMPLETED 1:28:14** | **P2 CLOSED** — first calibrated hybrid to survive since the `beta_from_phases` fix. No guard fire |
| **11608692** | ✅ COMPLETED 00:07, but its **number is contaminated** | gap 5 calibrated comparison — see the race below |
| **11608703** | 🟢 running 1:08:27 | coupler L₅₀ re-verification through the shared fixed code path |
| **11608941** | ✅ **COMPLETED 1:03:21** @ 14:00 | **gap 2** re-audited dispersive — **byte-identical**, as predicted. Folded in |
| **11609211_1** + ~~`_2`~~ → **11609680** + **11609681** | `_1` 🟢 running; `_2` ⛔ cancelled 46 min in, resubmitted as 11609680 (🟢, `-t 720`); comparison resubmitted as 11609681 | **the P4 test** — 2×2 gap 5 dispersive. **Rescued from a certain timeout** — see below |
| **11609215** `[1-2]` + ~~11609216~~ → **11609687** | `_1` ✅ **13:44**, `_2` 🟢 running; comparison resubmitted `afterok` | 1×2 gap 5 dispersive. **`_1`'s FDTD reference already folded in** — see below |
| **11608876** → **11609686** | ⏸ resubmitted `afterok:11608875` | gap 4 merge, moved to the funded account |
| **11608940** | ✅ COMPLETED 07:13 @ 13:04 | **gap 3** re-audited dispersive — folded in |

**Gaps 2 and 3 are both re-measured dispersive and both came back
byte-identical** (11608940 @ 13:04, 11608941 @ 14:00). Two independent
confirmations, so the rule is now evidence rather than analogy: **a quantity
evaluated at a scalar wavelength is invariant across the dispersion break;
only genuinely broadband quantities move** — gaps 4 and 5, and nothing else.
Spending an hour of compute on a result expected to be unchanged was the point.

**⚠ Do not reach for the ε grep on every job.** It works only where meep builds
a broadband medium and prints `dielectric constant epsilon diagonal`.
11608941 runs through meow/MPB and prints no such line at all — grepping
returns nothing, which is *not* evidence of frozen material. There its log
header carries the discriminator instead: `wl=1.31`, a single scalar
wavelength. Match the check to the code path.

**11608940's result, and the correction it forced.** Gap 3 came back
**byte-identical** to its frozen-era numbers. I had recorded beforehand that
gaps 2 and 3 were stale and that gap 2 was "most likely to move" — that was
**wrong**, and the job disproved it rather than confirming it. A dispersive
medium only differs from a frozen one *across a band*; `section_mode_neff()`
builds its medium at a scalar wavelength deliberately, and meow's
`Environment(wl=…)` is exact Sellmeier either way. Only the genuinely broadband
gaps (4 and 5) move. 11608941 is left running to measure that for gap 2 rather
than infer it by analogy. Full account in `resolved/mmi_1x2_validation.md`.

**Both gap-5 pairs preserved their frozen baseline before launching** —
`comparison_FROZEN_pre20260801.npz` in each of
`example_output/verification/mmi_{1x2,2x2}_hybrid_vs_fdtd/`. The comparison
script overwrites `comparison.npz`; without the copy the pre-flip half of the
P4 test would have been unrecoverable.

**Watch item on the 2×2 hybrid.** It runs at `N_MODES = 14` with zero margin
over the guided count — a duplicate-mode guard trip there fails the job rather
than degrading it.

**⚠ 11609211_2 was cancelled 46 min in because it was going to time out, and
the way that was caught generalises.** `run_mmi_2x2_examples.batch` carries one
`-t 300` shared by all three of its examples, sized when the slowest fitted.
The frozen-era hybrid (11585844_2) took **04:11:55** against that 300 — **48
minutes of headroom**. Dispersion costs ~1.5× (measured: 1×2 FDTD 08:57 → 13:44,
jobs 11585845_1 → 11609215_1), putting the dispersive hybrid near **06:25**.
It was never going to finish. Resubmitted as **11609680** with `-t 720`.

- **The limit cannot be raised in flight.** `scontrol update TimeLimit` returns
  `Access/permission denied` for a non-admin, so the only remedy is cancel and
  resubmit. That makes the check worth doing *at submit time*, not after.
- **The check itself is cheap and should be routine:** before re-running
  anything under the new material model, look up the frozen-era elapsed with
  `sacct -X -j <old> --format=Elapsed,Timelimit` and multiply by ~1.5. A shared
  array time limit is the usual trap, because it is sized for whichever task
  happened to be slowest at the time.
- **The 46 min are not lost** — `example_sim_cache` persists per solve, so
  completed section solves are reused on the rerun.

**⚠ The second hazard was worse, and it was invisible.** The gap-5 comparison
jobs depended on their producers with **`afterany`**, which a TIMEOUT satisfies
just as well as success. Had the hybrid timed out, the comparison would have
fired anyway and compared the **new dispersive FDTD** against the **old frozen
hybrid pickle still on disk** — yielding a plausible, wrong, unlabelled gap 5
number. That is precisely what corrupted **11608692** earlier today, caught only
because `provenance.json` happened to record an mtime.

**Both comparisons are now `afterok`** — 11609681 (2×2) and 11609687 (1×2).
Use `afterok` for any job that consumes a producer's output artifact;
`afterany` is only correct when the dependent job's job is to *report on*
whatever happened, not to read a result.

**Account.** New submissions go to **`gts-sr108-ece`**, which has the balance;
`gts-sr108-gedc_main` is what the older batch scripts hardcode. Applies to
newly run jobs only — jobs already queued or running are left alone.
`scontrol update JobId=<id> Account=...` does not work here regardless: it
fails with a misleading `Requested time limit is invalid`.

**11609215_1's number, folded in without waiting for its merge.** The
full-device FDTD reference alone is enough to measure how far the *reference*
moved across the dispersion break, which is worth having separately from gap 5
itself: max |ΔT| **0.00651** (1.65% of mean T), `out_0` and `out_1` agreeing to
five decimals. Read against the P4 hypothesis it is **mixed** — the scale is
2× the 2×2's entire ΔT swing, which supports it, but the curve crosses zero
**twice**, which a single frozen-index term cannot do. Both halves recorded, in
`resolved/mmi_1x2_validation.md` gap 5 (table) and `resolved/mmi_2x2_validation.md` gap 5
(P4 box), **before 11609212 lands**.

**⚠ A pickle-path race corrupted 11608692, and it will recur.** The comparison
scripts load their FDTD reference **by path**, so they read whatever last wrote
there. 11609215_1 overwrote `example_output/mmi_1x2/mmi_1x2_fdtd/smatrix.pkl`
with the *dispersive* reference at 17:34:48 UTC; 11608692 ran at 17:43 and
compared it against 11608690's *frozen* calibrated hybrid. Its printed
`max|ΔT| = 0.0313` mixes two material eras and must not be quoted. The valid
frozen-era number (**0.02876 / 7.30%**) was recovered offline from the
preserved baseline at zero compute cost — see `resolved/mmi_1x2_validation.md` gap 5.

Three lessons, all cheap:

1. **Two jobs writing and reading one pickle path need an ordering constraint,
   and `afterany:` on the *producer* is not it** — 11608692 depended on
   11608690 (its hybrid), not on the unrelated job that clobbered its FDTD
   input. Nothing in the dependency graph could have caught this.
2. **It was caught only because `provenance.json` records `fdtd_pkl_mtime`.**
   That is P6's argument arriving from a second direction: extend provenance
   from goldens to comparison outputs, and record the **resolved material
   config**, not merely an mtime — an mtime only helps if someone already
   suspects a problem.
3. **The meep log settles material era in one grep**, independent of any
   provenance — see below.

**How to check a run actually got the dispersive material, in one line.** Worth
recording because "the default flipped" and "this job ran dispersive" are
different claims, and the second is the one that matters when folding a number
in. Every meep log prints `dielectric constant epsilon diagonal`. For `SOI_160`
silicon that is **12.4744** dispersive (ε∞ of the fit, plus one Lorentzian) and
**12.2774** frozen (exactly Sellmeier n² at band centre, 3.5039² ). The two are
1.6% apart and unmistakable:

```
python -c "import numpy as np; from beef.materials.stack import SOI_160 as s; \
  wl=np.linspace(1.26,1.36,11); \
  print([s.si_medium(wl,dispersive=d).epsilon_diag.x for d in (True,False)])"
```

**All three completed 11608875 tasks print 12.4744**, so gap 4 is genuinely
being re-measured dispersive. This does not depend on trusting the cache
fingerprint or the wall clock.

**Wall clock, by contrast, is not usable here.** Frozen-era per-length was a
clean monotone 18:30 / 25:33 / 27:01 for L = 0.5/1.0/1.5. Dispersive came back
**55:15 / 41:22 / 49:33** — surcharges of 2.99× / 1.62× / 1.83×, and
*inverted* in L, with the shortest section slowest. Compute cost cannot depend
on L that way, so these numbers are contention-dominated (1 rank inside a
2-node allocation, four tasks placed wherever the scheduler had room). Do not
build a cost model on them; the 1.53× from 11609215_1 below is the cleaner
figure and even that is one sample.

**Cost note.** 11609215_1 took **13:44** against the frozen-era 11585845_1's
**8:57** for the same FDTD example — a **1.53×** dispersion surcharge, from the
Lorentzian susceptibility meep now has to time-step. Real but modest. The
hybrid half is the expensive one either way (11585845_2 was **1:46:59**), so
size the eventual golden regeneration (P7) off that, not off this.

---

**★ THIRTEENTH SNAPSHOT — 2026-08-01, after material dispersion was turned on
by default. Superseded by the fourteenth above for job states; still current
for the generational-break reasoning.**

**Everything measured before this point is on the old material model.** The
`dispersive` default flipped to True (`material_model_validation.md` step 3),
which adds `resolved_dispersive` + the fitted `resolved_{si,sio2}_reparam` to
the cache fingerprint — so **every cache entry misses and every golden is stale
again.** There is a clean generational break here; do not mix numbers across it
without saying which side they came from.

| job | state | note |
|---|---|---|
| **11608448** `[0-3]` | ✅ all COMPLETED — 18:30 / 25:33 / 27:01 / **52:25** (≈2 h 04 total) | **frozen-material era.** Its cache entries are still on disk but unreachable under the new key |
| ~~11608450~~ | ⛔ **CANCELLED by me, deliberately** | the gap-4 merge runs `cache_only=True`; under the new key it would have raised `CacheMissError` through no fault of its own |
| **11608875** `[0-3]` + **11608876** | 🟢 re-submitted | gap 4 re-measured under the configuration that actually ships. Fully cache-cold, so expect ~2 h again |
| **11608690** | 🟢 running | 1×2 calibrated hybrid — **frozen-material era**, and that is fine: it exists to test the `beta_from_phases` fix, which is orthogonal |
| **11608692** | 🟢 released | gap 5 calibrated comparison. Safe across the break — it only loads two saved `SMatrix` pickles, never touching the cache. Verified before releasing rather than assumed |
| **11608703** | 🟢 running | coupler L₅₀ re-verification, also frozen-era, also orthogonal (shared `beta_from_phases` path) |

**Sequencing note worth keeping.** The two pending merges were put on
`scontrol hold` *before* the library edit, not after — a `cache_only=True` job
firing mid-edit would have failed spuriously and looked like a real cache miss.
One was then released (pickles only) and one cancelled (cache-dependent). That
distinction is worth making explicitly whenever a fingerprint changes under a
running queue.

---

**★ TWELFTH SNAPSHOT — 2026-08-01, earlier the same day. Superseded by the
thirteenth above for job states; still current for the parity-work rationale.**

Five jobs submitted, all for the **1×2 MMI parity work** (`README.md`'s new
parity backlog, items P1 and P2). One inert job cancelled (`11597199`, the
orphaned mux merge — done).

| job | name | for | state at snapshot |
|---|---|---|---|
| **11608448** `[0-3]` | `beef_mmi_1x2_engine_agreement` | **P1 / gap 4 / check 1** — this device's last unrun check | `_0` ✅ 18:30, `_1` ✅ 25:33, `_2` ✅ 27:01, `_3` 🟢 running |
| **11608450** | …`_merge` | assembles the ratio from cache | ⏸ `afterok:11608448` |
| **11608690** | `beef_mmi_1x2_hybrid_calibrated` | **P2** — successor to 11589770, first calibrated run since the `beta_from_phases` fix | 🟢 running |
| **11608692** | `beef_mmi_1x2_gap5_calibrated` | gap 5's calibrated comparison, successor to 11589779 | ⏸ `afterany:11608690` |
| **11608703** | `beef_dc_hybrid_item14` | **coupler re-verification** — the shared `beta_from_phases` fix touches the coupler's path too | 🟢 running |

**P1 — gap 4 cost, and the reshape paid off differently than on the 2×2.**
Per-length wall clock 18:30 / 25:33 / 27:01 for L = 0.5/1.0/1.5 — about **2 h
of serial work total** against the 2×2's 13.4 h, as the 6.2×-shorter body
(`mm_length` 4.5 vs 28) predicted. A single serial job would have fitted in one
slot here; the array shape was adopted anyway because it is free and a
per-length timeout then costs one length instead of four. Launched with bare
`python` (1 rank) deliberately — the script has no rank guards, so `srun` would
race 48 ranks on one `np.savez`.

**P2 — why these two reruns exist.** The 2026-07-31 diagnosis of 11589770's
failure (a π sign flip in `T_ii`) was **wrong**; the real cause was a 2π
turn-count bug in `beta_from_phases`, now fixed. See that job's row below for
the full account. 11608690 is the first calibrated run on the fixed code: if it
completes without the residual guard firing, `calibrate_body_neff` works on
this device for the first time and 11608692 can finally measure gap 5
calibrated.

**11608703 is the backward-compatibility check, and it is deliberately a real
run rather than an argument.** The fix is in shared library code that the
coupler's own closed item 6 depends on. Offline evidence is strong — at the
coupler's production seed the fix changes **0 of 66 cached cells**, and 28 of
29 existing tests were unchanged — but that used a *proxy* seed rather than the
run's own MPB output. This job re-derives the coupler's calibrated 50:50 length
through the real path; it should reproduce ~15.15 µm against the 15.13 µm
full-device reference. Chose `run_dc_hybrid_item14.batch` over
`run_golden_directional_coupler_hybrid.batch` on purpose: the golden generator
would **overwrite** `tests/golden/directional_coupler_hybrid.json`, destroying
the pre-fix reference this comparison needs (and goldens are deferred work
besides). The archived stub-sweep script would have been the closest match to
the original 0.11% measurement, but its `make` target no longer exists.

---

**★ ELEVENTH SNAPSHOT — 2026-08-01, verified against `squeue -u jhiesener3`
and `sacct -u jhiesener3 -X`. Superseded by the twelfth above for job states;
still current for the 2×2 gap-4 closure and the mux diagnosis.**

**The queue is now empty of BEEF work.** Nothing is running. One inert entry
remains: **11597199** `beef_modal_mux_examples_merge`, `PENDING` with
`Reason=DependencyNeverSatisfied` — its parent array timed out, so it can
never start. It holds no resources; cancel it as hygiene (`scancel 11597199`).
Not done automatically. (`11585814_7` `elevator_resolution` COMPLETED 07:28
after 22:47:08 — **unrelated project, ignore**.)

Both reshapes submitted in the 21:20 window have now landed, with **opposite
outcomes**:

| job | name | outcome |
|---|---|---|
| **11596845** `[0-3]` | `beef_mmi_2x2_engine_agreement` | ✅ **ALL FOUR COMPLETED.** `_0`/`_1` in **36 s each** (predicted cache hit on 11591249's L=2/L=4); `_2` (L=6) 3:35:10; `_3` (L=8) **5:07:14** |
| **11596851** | …`_merge` | ✅ **COMPLETED** 03:07, 25:29 — assembled the answer from cache, `cache_only=True` raised nothing |
| **11597177** `[0-2]` | `beef_modal_mux_examples` | ❌ **ALL THREE TIMEOUT** (02:44/02:52/03:01), full 5 h × 4 nodes each — **third consecutive timeout**, cause now understood |
| **11597199** | …`_merge` | ⛔ **orphaned**, `DependencyNeverSatisfied`, cancel it |

**✅ The 2×2 MMI's gap 4 is CLOSED — its last open check.** Mean ratio
Δβ_FDTD/Δβ_EME = **0.9953**: the engines agree on body phase accumulation to
**0.47%**, versus the coupler's ~10% on the identical measurement, so
`calibrate_body_neff` is **not** indicated on this device. T diagonal to
1.21e-05, worst phase-fit residual 4.02e-06 rad (25,000× inside the 0.1 rad
guard that broke the 1×2's calibration). The job log prints both as `0.0000`
at `:.4f` — they are genuinely that small; read off
`example_output/verification/mmi_2x2_engine_agreement/engine_agreement.npz`,
**not** defaulted or skipped, which was checked before recording this.
Folded into `resolved/mmi_2x2_validation.md` gap 4 (boxed) and `README.md`'s table.
Caveat carried there: measured at res 40, which gap 2 says is not converged in
absolute n_eff; both engines read at the same mesh so the ratio partially
cancels it, but 0.47% is not certified mesh-independent.

**The reshape's cost model was right.** Per-length wall clock 96 / 183 / 215 /
**307** min for L = 2/4/6/8 — **13.4 h of serial work**, so the original 5 h
single-slot job could never have finished, and L=8 alone would have blown it.
As an array it closed in 5:07 + 25 min. Worth copying for any other
multi-length `solve_section()` sweep.

**❌ The mux reshape did the opposite, and the reason is a mistake worth not
repeating: it silently dropped MPI.** The new batch script calls `python
modal_multiplexer/python/modal_multiplexer_hybrid.py` directly to avoid
`make` running all 3 ports in one process — but the Makefile target is
`$(runner) $(my_python) …` and on Phoenix `runner` = **`srun --mpi=pmi2`**
(`BEEF/examples/Makefile:5`). Bypassing `make` threw away the launcher with
the target:

| job | launcher | meep ranks | s/step |
|---|---|---|---|
| 11585847_0 / 11590663_0 (via `make`) | `srun --mpi=pmi2` | **96** | 0.0064 |
| 11597177_0/_1/_2 (reshaped, direct) | none | **1** | 0.27–0.36 |

`Using MPI version 4.1, 1 processes` in all three new logs against
`96 processes` in both old ones — **~42× slower per timestep**, while holding
4 nodes × 24 ranks idle (`cput=20-00:35:12` for `walltime=05:00:22` on 96
cpus). **Second, independent flaw:** the reshape assumed sections 0-3 were
already warm on disk; they were not — each task re-ran the same three
`Initializing structure` + MPB + time-stepping blocks at near-identical log
offsets. Whether the earlier timeouts never persisted them, or the cache key
differs between 96-rank and 1-rank runs, is **not established**.

**The mux is deferred, so no fourth attempt was made** (per instruction, and
consistent with `deferred/modal_multiplexer_validation.md`). Fix recipe
recorded in that doc's gap 6 for whenever it resumes: restore
`srun --mpi=pmi2`, *then* verify the cache actually hits, then resubmit — and
do not merely raise the wall clock, since at 1 rank no affordable limit
finishes. **Gap 6 has now consumed ~25 h of compute across three attempts and
is still unanswered either way.** The array-plus-merge *pattern* is not at
fault — the same design worked first try on the 2×2.

**Still open from earlier, untouched by this snapshot:** `beef_dc_examples`
`_0` (11592462_0) TIMEOUT, cause still undiagnosed — the standing
recommendation is a short-wall-clock resubmit watched live, **not** submitted
here. The mmi_1x2 π-sign-flip diagnosis and the ring y-normal excite blocker
are both analysis work with no job in the queue.

---

**Snapshot taken 2026-07-31, 16:35 EDT** (ninth snapshot of the day: 07-30
21:46 → 07-31 08:22 re-queue → 08:46 → 10:35 → 10:53 → 12:05 → 13:05 → 13:44 →
now). This update: mmi_1x2 gap 2 (11591224) completed and closed (res 120,
`RESOLUTION` held at 40 by decision — see `resolved/mmi_1x2_validation.md`); the
coupler's cold regeneration pass launched, 11592462 (example outputs,
array 0-2) + 11592463 (golden) — see the flagged anomaly on array task 2
below before treating that output as a clean fresh result. This is a point-in-time list, not a
live view — job states, elapsed times and queue positions all move. Check
`squeue -u jhiesener3` (or `sacct -j <id>` once finished) for current status;
this file exists so the *purpose* of each job is recorded somewhere durable,
not to be a dashboard. Update or delete entries as jobs complete and their
results get folded into the relevant validation doc.

**Reshaped and resubmitted, 2026-07-31 ~21:20 EDT.** Both `beef_modal_mux_
examples` and `beef_mmi_2x2_engine_agreement` timed out with an understood
cost problem (see their entries below), so both scripts were reshaped into
an array-job-plus-merge-job pair instead of one long serial process, and
resubmitted:

| new job | replaces | shape |
|---|---|---|
| **11596845** `[0-3]` | 11591249 | `verify_mmi_2x2_engine_agreement.py`, now one array task per length (2/4/6/8 µm), 8h/task |
| **11596851** | (new) | merge step for 11596845 — `--dependency=afterok:11596845`, assembles the ratio/fit from cache, 30 min |
| ~~11596846~~ → **11597177** `[0-2]` | 11590663 | `modal_multiplexer_hybrid.py`, one array task per output port (top/mid/bot), 5h/task |
| ~~11596852~~ → **11597199** | (new) | merge step — `--dependency=afterok:11597177`, assembles the SMatrix from cache, 30 min |

See `verify_mmi_2x2_engine_agreement.py`'s and `modal_multiplexer_hybrid.py`'s
own "Reshaped 2026-07-31" docstring sections, and the four
`batch/.../*.batch` files, for the mechanics (`SLURM_ARRAY_TASK_ID` selects
one length/port; `cache_only=True` in the merge step turns a missing result
into a loud `CacheMissError` instead of a silent multi-hour re-simulation).
**`beef_dc_examples` idx 0 (11592462_0) also timed out** during this same
window (21:34) but its cause is not understood yet — **not** reshaped or
resubmitted; see its own entry below before touching it.

**Correction, 21:34-21:41 EDT:** 11596846 `[0-2]` (the modal-mux array,
first attempt) **FAILED instantly** (~8s each task) — the new batch script
had the wrong relative path to `modal_multiplexer_hybrid.py`
(`python/modal_multiplexer_hybrid.py` instead of
`modal_multiplexer/python/modal_multiplexer_hybrid.py`, a mistake made while
writing this reshape, not a cluster/library issue). Its merge job
(11596852) was left stuck on `DependencyNeverSatisfied` and was cancelled.
Fixed and resubmitted as **11597177** (array) / **11597199** (merge,
`--dependency=afterok:11597177`). **11596845** (mmi_2x2 engine agreement
array) had no such bug: tasks `_0`/`_1` (L=2, L=4) both **completed in 36s**
each, confirming the predicted cache hit on 11591249's already-computed
lengths; `_2`/`_3` (L=6, L=8 — the real new work) are running for real as of
this snapshot.

**Re-verified against `squeue`/`sacct` 2026-07-31, 21:12 EDT (tenth snapshot).**
Two jobs that were "🟢 RUNNING" in the 16:35 snapshot both **hit their 5h wall
and were killed by SLURM, not by an error**: `beef_modal_mux_examples`
(11590663_0) at 18:58, and `beef_mmi_2x2_engine_agreement` (11591249) at
19:46. Neither is a crash — both show `ExitCode=0:0`, `State=TIMEOUT` — but
neither produced its intended result: the 2×2 MMI's gap 4 (engine agreement,
its last open check) did **not** close, and the mux's gap 6 exercise run
(`warn_on_transfer_duplicate` against real interfaces) is **still** unproven
either way. `beef_golden_dc_hybrid` (11592463) completed cleanly (1:09:17).
`beef_dc_examples` array task 0 (11592462_0, hybrid) is still running at
4:38:25 of its own 5:00:00 limit as of this snapshot — **~21 min of budget
left; watch this one, it may time out the same way** unless it finishes
first. See the three detail sections below (11590663, 11591249, 11592462) for
what each actually reached before stopping.

**Re-queue, 2026-07-31 08:22 EDT:** the five jobs below never started —
`scontrol show job` confirmed `Account=gts-sr108-ece`,
`Reason=AssocGrpBillingMinutes` on the three parents (that association is out
of billing minutes; `gts-sr108-gedc_main` had headroom). Cancelled all five
(old IDs 11577985/11577986/11578018/11578019/11578026) and resubmitted under
`gts-sr108-gedc_main`. Also switched all other batch scripts hardcoding
`gts-sr108-ece` to `gedc_main` so this doesn't recur.

| job | name | state (at snapshot) | for |
|---|---|---|---|
| 11575516 | `beef_mode_spec_converge` | ✅ **COMPLETED** 21:17, 3:36:13 — results folded in, entry retained for provenance | [resolved/fdtd_vs_EME_mode_validation.md](resolved/fdtd_vs_EME_mode_validation.md) item 1 + [modal_multiplexer_validation.md](deferred/modal_multiplexer_validation.md) gap 2 |
| ~~11577985~~ → **11585844** `[1-2]` | `beef_mmi_2x2_examples` | ✅ **COMPLETED** — `_1` 1:16:45, `_2` 4:11:55; both SMatrix files on disk | [resolved/mmi_2x2_validation.md](resolved/mmi_2x2_validation.md) gap 5 (data) |
| ~~11577986~~ → **11585849** | `beef_mmi_2x2_gap5` | ✅ **COMPLETED** 12:46, 0:18 — **worst max\|ΔT\| 0.0043 (0.82%), sign FLIPS across the band (not the 1×2's uniform signature); results folded in** | [resolved/mmi_2x2_validation.md](resolved/mmi_2x2_validation.md) gap 5 (comparison) |
| ~~11578018~~ → **11585845** `[1-2]` | `beef_mmi_1x2_examples` | ✅ **COMPLETED** — `_1` 8:57, `_2` 1:46:59; both SMatrix files on disk | [resolved/mmi_1x2_validation.md](resolved/mmi_1x2_validation.md) gap 5 (data) |
| ~~11578019~~ → **11585850** | `beef_mmi_1x2_gap5` | ✅ **COMPLETED** 10:28, 0:09 — **max\|ΔT\| 0.0359 (9.11%), results folded in** | [resolved/mmi_1x2_validation.md](resolved/mmi_1x2_validation.md) gap 5 (comparison) |
| ~~11578026~~ → ~~11585847~~ → **11590663** `[0]` | `beef_modal_mux_examples` | 11585847 ❌ **TIMEOUT** 13:42 (full 5h, 4 nodes) — sections 0-3 cached to disk before it died mid-section-4; resubmitted unchanged as 11590663 → ❌ **TIMEOUT again**, 18:58, full 5h — got further (deep into a per-port section-4 eigenmode-overlap monitor, band 13/14) but still did not finish or reach the item-8 `warn_on_transfer_duplicate` scan this run exists to exercise → reshaped as **11597177** `[0-2]` → ❌ **THIRD TIMEOUT 2026-08-01**, all 3 tasks, full 5h × 4 nodes each. **Cause found: the reshaped batch script dropped `srun --mpi=pmi2`, so meep ran on 1 rank instead of 96 (~42× slower per step) with 96 cores idle**; and the "sections 0-3 are cached" premise it rested on is false — every task recomputed them. **NOT resubmitted, device is deferred**; merge 11597199 orphaned on `DependencyNeverSatisfied`, cancel it. ~25 h of compute spent, gap 6 still unanswered | [modal_multiplexer_validation.md](deferred/modal_multiplexer_validation.md) gap 6 (exercise run) |
| ~~**11586046**~~ → **11588309** | `beef_diagnose_items_11_12` | 11586046 ❌ **CANCELLED** 11:12 (hung, not slow); refixed → **11588309 ✅ COMPLETED 11:50, 34:55 — BOTH ITEMS RESOLVED:** item 11 closed (all four grids agree at res 40 too), item 12 decisive (laddering recovers a 14th guided mode at equal cost) | [resolved/fdtd_vs_EME_mode_validation.md](resolved/fdtd_vs_EME_mode_validation.md) items 11 + 12; it also opened item 14 |
| **11586098** | `beef_diagnose_item_13` | ✅ **COMPLETED** 10:50, 1:54:50 — **C vs A inert (max \|Δn_eff\| 0.000e+00), cost 3.2×; results folded in** | [resolved/fdtd_vs_EME_mode_validation.md](resolved/fdtd_vs_EME_mode_validation.md) item 13 |
| ~~11586986~~ → ~~11588997~~ | `beef_dispersive_convergence` | 11586986 ❌ **CANCELLED** 12:00 — **it was measuring nothing**, `dispersive=` never reached the mode solve (library bug, fixed); refixed → **11588997 ✅ COMPLETED** 13:56, 1:56:53 — **pin holds at centre, holds at the edge only up to the 1e-3 tolerance production actually uses (needs res 80, not 60, under dispersion — which is what's already shipped); does not fully converge at 3e-5 at the edge within res 20–120** | [material_model_validation.md](resolved/material_model_validation.md) work item 4 (cheap half) |
| **11588569** | `beef_audit_mmi_2x2_gap6` | ✅ **COMPLETED** 11:58, 20:18 — **gap 6 closed:** identity `[0…13]`, worst \|Δn_eff\| 0.0007, 14/14 guided, **neither guard fired** | [resolved/mmi_2x2_validation.md](resolved/mmi_2x2_validation.md) gap 6 |
| **11589770** → **11608690** | `beef_mmi_1x2_hybrid_calibrated` | ❌ **FAILED** 16:59, 4:17:26 — `section_mode_beta`'s phase-fit residual (2.101 rad) blows its own 0.1 rad sanity guard on the `[-2.0, 2.0]` section, 3 lengths. ~~DIAGNOSED 2026-07-31 as a **π sign flip in `T_ii`** at 4 of 88 cells; narrow upstream mode-extraction bug; the guard caught a real defect~~ → **that diagnosis was WRONG. ✅ REAL CAUSE FOUND AND FIXED 2026-08-01: a 2π turn-count bug in `beta_from_phases` itself.** It rounded wrap counts independently per length against one seed slope, so the ~2% MPB-vs-FDTD seed error — *the quantity `calibrate_body_neff` exists to measure* — drifted with L until a length straddled a `.5` boundary. The two readings are algebraically indistinguishable on the residual (both give `δ·(+⅓,−⅔,+⅓)`); physics separates them: turn-count → n_eff **2.783, +1.96%** over MPB, in family with the coupler's +2.233%/+2.361%; sign-flip → **2.101, −23%**, wrong sign and size. Seed-free proof: second difference `(φ₃−2φ₂+φ₁)/2π` is an exact integer to **0.00308 in all 88 cells** — nothing needs a sign flip. Not a corner case either: 4025-trial sweep, old code wrong 13.7%, new 0%. The doc's own "unresolved" branch margin 0.922 **was** the tell → now 0.0785. **The guard did catch a real defect — its own unwrapping.** Re-run in flight as **11608690** | [resolved/mmi_1x2_validation.md](resolved/mmi_1x2_validation.md) gap 5 (boxed, rewritten) |
| **11589779** | `beef_mmi_1x2_gap5_calibrated` | ❌ **FAILED** 16:59, 0:00:09 — dependent on 11589770, nothing to load, expected given the above | [resolved/mmi_1x2_validation.md](resolved/mmi_1x2_validation.md) gap 5 follow-up (comparison) |
| **11590785** | `beef_mmi_2x2_gap2_dbeta` | ✅ **COMPLETED** 1:05:48 — **gap 2 closed at res 100**, Δβ_01 converges to 0.05% by res 120 | [resolved/mmi_2x2_validation.md](resolved/mmi_2x2_validation.md) gap 2 |
| **11591224** | `beef_mmi_1x2_gap2_dbeta` | ✅ **COMPLETED** 15:45, 1:05:15 — **gap 2 closed at res 120** (non-monotone step, unlike the 2×2's clean decay); `RESOLUTION` held at 40 by decision, gap recorded as a known mesh caveat instead | [resolved/mmi_1x2_validation.md](resolved/mmi_1x2_validation.md) gap 2 |
| **11591249** → **11596845** `[0-3]` + **11596851** | `beef_mmi_2x2_engine_agreement` | ✅ **CLOSED 2026-08-01 by the reshape** — array all 4 COMPLETED (`_0`/`_1` 36 s cache hits, `_2` 3:35:10, `_3` 5:07:14) + merge 25:29. **Ratio 0.9953 — engines agree to 0.47%, `calibrate_body_neff` NOT indicated on this device.** T diagonal to 1.21e-05, fit residual 4.02e-06 rad (log shows both as `0.0000` at `:.4f`; real values read from the `.npz`). Gap 4 was this device's last open check. Original attempt below: ❌ **TIMEOUT** 19:46, full 5h — 4 of the 4 lengths' first run completed (2/4/6/8 µm all show a `run 0 finished` line), a 5th run (second pass, presumably for the reverse-direction or FDTD-vs-EME comparison half) was ~5,300 of a likely ~30,000 timestep run when killed — well under half done on that leg. Gap 4 (this device's last open check) did **not** close | [resolved/mmi_2x2_validation.md](resolved/mmi_2x2_validation.md) gap 4 / [README.md](README.md) check 1 — never measured for either MMI; 4 real FDTD solve_section() runs, adapted from the coupler's item 6 script |
| **11592462** `[0-2]` | `beef_dc_examples` | `_0` (hybrid) ❌ **TIMEOUT** 21:34, full 5h — see detail section below, not yet diagnosed; `_1` (mode) ✅ **COMPLETED 16:37, 3:49** — legit, not flagged: serial MPB even/odd supermode sweep (1 process), inherently cheap, no BEEF solver/cache involved; L_c = 32.813 µm at λ=1.55, gap sweep 0.15–0.5 µm recorded, `_2` (fdtd) ✅ **COMPLETED 16:33, 0:00:17** — **flagged, not yet trusted as a clean result:** 17s for a full 3-D broadband FDTD sweep is a cache hit, not a fresh compute. Confirmed NOT the item-8 fingerprint bug (`FDTDSolver.solve()` does splat `resolved_material_config()` into its cache key, `fdtd.py:707`, inside `solve()` itself) — so whatever it hit was already computed post-item-8, most likely the ground-truth length_sweep's own 15.13 µm grid point. Printed numbers (T=0.5061/0.4933, IL≈-3 dB) look physically sane, but this has not been cross-checked against that sweep's stored value yet — do that before citing this run's output as the "outstanding — not yet run" regeneration `resolved/directional_coupler_validation.md`'s Verification section calls for | `resolved/directional_coupler_validation.md` Verification → Outstanding (example-output regeneration) |
| **11592463** | `beef_golden_dc_hybrid` | ✅ **COMPLETED** 17:42, 1:09:17 | `resolved/directional_coupler_validation.md` Verification → Standing (goldens stale since item 8) |
| ~~11589771~~ ~~11589899~~ ~~11589936~~ | `beef_modal_mux_fdtd` | ❌ **BLOCKED, do not resubmit** — three cheap failures (~15 s each) converged on a real library gap: `CompoundGeometry` cannot back `FDTDSolver.solve()` at all | [modal_multiplexer_validation.md](deferred/modal_multiplexer_validation.md) gap 3 |

**Six items resolved today; three jobs cancelled as unproductive, two refixed
and rerun, one left blocked.** Since the 12:05 snapshot: the **2×2 MMI's whole
gap-5 chain finished** (11585844_2, then 11585849) at 0.82% worst |ΔT| — an
order of magnitude tighter than the 1×2's 9.11%, but with a sign that flips
across the band rather than the 1×2's uniform under-prediction, so
`calibrate_body_neff` is not obviously indicated the same way here. A
follow-up calibration run for the **1×2** (which does show the uniform
signature) is in flight — 11589770/11589779. And **modal_multiplexer_fdtd.py**,
written to give the mux its first full-device FDTD reference (gap 3), hit a
real, unfinished library gap: `CompoundGeometry` cannot back
`FDTDSolver.solve()`'s single-domain path at all (only `solve_section()`,
what `HybridSolver` uses). One sub-bug in the same area WAS fixed and tested
(`CompoundGeometry.source_aperture_at()` used to raise past its own footprint;
now clamps, matching how every leaf geometry already behaves) — necessary but
not sufficient. Three cheap failed submissions (~15 s each) before landing on
the real blocker; not resubmitting until the library gap is actually closed.
See the mux's own detail section below and `modal_multiplexer_validation.md`
gap 3 / `notes.md` Geometry item 8.

Before that: **items 11 and 12** both resolved by 11588309 (34:55,
after its predecessor hung), **2×2 MMI gap 6** closed by 11588569, and a library
bug found by reading 11586986's output — `dispersive=` never reached
`section_mode_neff()`, so that job's two columns were byte-identical and work
item 4 was never actually started. Fixed, 9 tests added, re-running as 11588997
with a liveness check that aborts if the flag is ever inert again.

Item 13 is closed by 11586098 (the
guard's over-request is a pure filter, and `_DUP_MODE_MAX_EXTRA` is now a
measured 3.2× of pure cost) — folded into
`resolved/fdtd_vs_EME_mode_validation.md`. The 1×2 MMI's whole gap-5 chain completed end
to end earlier and its numbers are in `resolved/mmi_1x2_validation.md`. `11585847` has
produced a real finding before finishing: the item-0 duplicate guard fires on
`EMESolver`'s **production** path on this device at res 40 (folded into
`modal_multiplexer_validation.md` gap 6), though the item-8 `T`-column scan it
was submitted to exercise is still silent. `11585849` (2×2 comparison) is
pending only on `11585844_2`, which is correct.

**11586046 HUNG and was cancelled at 11:12; resubmitted as 11588309 after a
script fix.** It was not slow — it was stuck, and an earlier revision of this
file diagnosed it wrongly as a wall-clock budget miss. The real cause, from its
own MPB trace: `N_TARGET = 20` distinct modes is unattainable at that
cross-section (today's fixed seed finds 13 distinct guided modes there, and a
laddered seed 14 — job 11588309), so the search
walked up to `band_num 23`, where **no k reaches the target frequency at all**
and MPB's k-search limit-cycles between k = 0.2346 and k = 1.7480 indefinitely.
It burned 161 MPB solves on that one band — 59% of the job's entire eigensolve
work — while bands 1–22 had each converged in 3–8. Bands 1–22 took ~58 min for
22 attempts, so the original ~25 s/solve budget was about right: wall clock was
never the problem, and the 4-task split this file recommended at 11:03 would
have produced four hung tasks.

Two results came out of the failure, both folded into
`resolved/fdtd_vs_EME_mode_validation.md`: item 12 gains the measured band table showing
this device's shipped **16-mode basis is unachievable at res 40 by any seeding
policy** — since superseded by 11588309, which found laddering recovers a 14th;
the seeding half IS fixable, the count is 14 not 13 — and the production guard's
identical exposure to the same non-convergent solve is filed as new **item 14**.
Full detail under 11586046 and 11588309 below.

The detail sections below carry the current job IDs — for the five re-queued
jobs only the ID and billing account changed, so each entry's script, resources
and rationale are unchanged from the 07-30 submission.

## Details

### 11575516 — `beef_mode_spec_converge` — ✅ COMPLETED
- **Script:** `batch/verification/run_mode_spectrum_convergence.batch` →
  `examples/verification/python/mode_spectrum_convergence.py`
- **Resources:** 1 node, 8 cpus, serial (no MPI), 8h limit. Used 3:36:13.
- **Output:** `reports/11575516.out`.
- **What it did:** Swept meep and meow resolution *independently* (meep
  20/30/40/60, meow 20/30/40/60/80) on the modal multiplexer's interface 0
  and counted how many of the 16 requested modes each engine calls guided,
  and how many modes the two engines agree on within tolerance
  (`n_agree`), at every (meep_res, meow_res) pair.
- **Why:** discriminate the leading hypothesis for item 1's meep↔meow orphan
  mode — `n_agree` rising toward 16 with resolution would mean a
  band-window/discretization artifact (fixable with a per-engine converged
  resolution); a plateau would mean something structural.
- **Result — both, and the structural half is the one that matters.**
  `n_agree` does rise (12 → 15 at matched res 40 → 60) and the guided orphan
  count falls from 3 to 1, so most of the mode-8+ divergence is
  discretization. But the surviving orphan is present at **every** resolution
  20–60: meow finds two members of a near-degenerate pair (converging to a
  3.1e-4 splitting by res 80) where meep's per-band ladder finds one. That is
  a driver defect, not a mesh one. Full reading under item 1; the new to-dos
  it opened are items 11–13.
- **Second result, unplanned:** the first per-engine converged-resolution
  measurement on SOI_160 — meow plateaus between res 60 and 80 (16/16), meep
  is still moving at 60. The shipped shared `RESOLUTION = 40` is below both
  (`modal_multiplexer_validation.md` gap 2).
- **Correction to this entry's original filing:** it was listed as serving
  item 6 (`interface_modes`/layer-2 convergence) as well. It does not — that
  item sweeps the requested mode *count* and requires the *S-matrix* to
  plateau, neither of which this job did. See item 6's note.

### 11585844 `[1-2]` — `beef_mmi_2x2_examples` (array indices 1,2 only) — ✅ COMPLETED
- **Script:** `batch/run_mmi_2x2_examples.batch` → `make mmi_2x2_fdtd
  mmi_2x2_hybrid`
- **Resources:** 2 nodes × 24 ranks per array task, 5h limit each.
- **What it does:** Runs the 2×2 MMI's full-device FDTD reference
  (`mmi_2x2_fdtd.py`) and the FDTD/EME hybrid cascade
  (`mmi_2x2_hybrid.py`), each saving an `SMatrix` pickle. Array index 0
  (`mmi_2x2_eme`) deliberately skipped — not needed for gap 5.
- **Why:** gap 5 ("no hybrid-vs-full-device comparison on record") needs
  both scripts' saved output to exist before anything can compare them.
- **Result:** `_1` (fdtd) 1:16:45, `_2` (hybrid) 4:11:55. Both SMatrix pickles
  on disk; see 11585849 for the comparison.

### 11585849 — `beef_mmi_2x2_gap5` — ✅ COMPLETED
- **Script:** `batch/verification/run_verify_mmi_2x2_hybrid_vs_fdtd.batch` →
  `examples/verification/python/verify_mmi_2x2_hybrid_vs_fdtd.py`
- **Resources:** 1 node, 1 cpu, 20 min. `--dependency=afterany:11585844`.
- **What it does:** Loads both SMatrix pickles from 11585844, reconciles
  port names (`out_0`/`out_1` vs `bar`/`cross`), reports per-wavelength
  transmission difference with provenance (git sha, config, mtimes). No
  meep/meow import cost of its own beyond loading `beef`.
- **Why:** closes gap 5 as measurement-and-report (no threshold enforced —
  gap 4, engine agreement, isn't measured for this device yet, so there's
  no basis to call a given disagreement expected vs. a bug).
- **Result (12:46, 18 s):** `check_passive()` true both sides. Worst max|ΔT|
  **0.0043** (cross port, λ=1.5650 µm, 0.82% of mean T_fdtd); bar port 0.0028
  (0.62%). Under 1% on both — far tighter than the 1×2's 9.11% — but the sign
  is **not** uniform: both ports start positive, cross zero near 1.544 µm, end
  negative. Not the 1×2's Δβ-discretization signature (uniform negative), so
  `calibrate_body_neff` is not obviously indicated here the way it is for the
  1×2 — a sign flip looks more like a resolution/dispersion-sensitive phase
  term. Also measured at the shipped unconverged `RESOLUTION = 40` (gap 2), so
  part of the flip could be a mesh artifact rather than physics; the two
  cannot be separated until gap 2 is measured on this device. Folded into
  `resolved/mmi_2x2_validation.md` gap 5. Output under
  `example_output/verification/mmi_2x2_hybrid_vs_fdtd/`.

### 11585845 `[1-2]` — `beef_mmi_1x2_examples` (array indices 1,2 only)
- Same shape as 11585844, for the 1×2 MMI: `make mmi_1x2_fdtd
  mmi_1x2_hybrid`.
- **Why:** `resolved/mmi_1x2_validation.md` gap 5, same reasoning as the 2×2's. This
  device's basis mismatch (gap 3, `N_MODES` 12→8) is already fixed, so —
  unlike the 2×2 — there's no open mode-identity concern underneath this
  comparison; only gap 2 (resolution) remains unsettled.

### 11585850 — `beef_mmi_1x2_gap5` — ✅ COMPLETED
- **Script:** `batch/verification/run_verify_mmi_1x2_hybrid_vs_fdtd.batch` →
  `examples/verification/python/verify_mmi_1x2_hybrid_vs_fdtd.py`
- Same shape as 11585849, ported to `mmi_1x2`'s port names (`out_0`/`out_1`
  vs `out_top`/`out_bot`). `--dependency=afterany:11585845`.
- **Why:** `resolved/mmi_1x2_validation.md` gap 5.
- **Result (10:28, 9 s):** max |ΔT| **0.0359** at λ = 1.33 µm, mean |ΔT|
  0.0253, **9.11%** of mean T_fdtd. `check_passive()` true both sides; the two
  output ports agree to 4 decimals, so splitter symmetry survives the cascade.
  **Every ΔT is negative** — a systematic under-prediction, not scatter, which
  is the Δβ-discretization signature from
  `resolved/directional_coupler_validation.md` items 6/14 (coupler: ~5%
  uncalibrated → 0.11% with `calibrate_body_neff=True`, never tried here).
  Folded into `resolved/mmi_1x2_validation.md` gap 5; output under
  `example_output/verification/mmi_1x2_hybrid_vs_fdtd/`. Provenance records
  beef `1cc8edc`.

### 11585847 `[0]` → 11590663 — `beef_modal_mux_examples`
- **Script:** `batch/run_modal_multiplexer_examples.batch` → `make
  modal_multiplexer_hybrid`
- **Resources:** 4 nodes × 24 ranks, 5h limit — the most expensive single
  job this session; this device is the repo's "most demanding hybrid
  device" (16 modes × 4 interfaces, SOI_160).
- **What it does:** This device's **first real run** in this tree — no
  `example_output/modal_multiplexer/` existed before this job was
  submitted.
- **Why:** (a) gives the least-supported device in the repo (no golden, no
  FDTD reference, shared resolution knob) its first real output; (b)
  exercises `beef.mode_identity.transfer_matrix_duplicate_scan()` /
  `warn_on_transfer_duplicate()` (`resolved/fdtd_vs_EME_mode_validation.md` item 8)
  against this device's real `solve_section()` interfaces for the first
  time — `modal_multiplexer_validation.md` gap 6 flags this as "not yet
  exercised... needs an actual meep/meow run." Watch
  `reports/11585847_0.out` for `warn_on_transfer_duplicate` firing on
  either solver.
- **Partial result at 10:53 (2h12 in, still running):** two FDTD subregions
  done; the first EME body solve after them fired the **item-0 guard** —
  `meow mode 9 duplicates already-accepted mode 8 (self-overlap 1.000000)`,
  logged once per rank across all 96 ranks. By 12:16 this had grown to **six
  distinct firings from two pairs** — 4x mode 9/8 and 2x mode 11/10, the latter
  being item 0's own original pair. First
  catch of the meow duplicate on the *production* `compute_modes()` path rather
  than in an audit script. `warn_on_transfer_duplicate` — the item-8 scan this
  job was actually submitted to exercise — has **not** fired; no `transfer`
  line appears in the log at all, which is what you would expect if the
  upstream guard rejected the duplicate before it reached `T`, but is not yet
  evidence the scan works. Full reading in `modal_multiplexer_validation.md`
  gap 6.
- **❌ TIMEOUT 2026-07-31 13:42, full 5h.** `sacct` confirms `State=TIMEOUT`,
  `ExitCode=0:0` — killed by SLURM mid-run, not a crash. The batch script's `-t
  300` was sized off the earlier device estimates, not this device's own cost:
  this is the repo's most expensive hybrid chain (16 modes × 4 interfaces) and
  its own `save_dir`/`cache=` mean a re-run doesn't restart from zero, so the
  fix is a plain resubmit rather than a reshape (contrast the 11586046 case,
  where the batch script's *cost model* was wrong and needed reshaping, not
  just more wall clock).
  - **What finished, read off `example_sim_cache/` mtimes against the
    script's fixed section order (0 input FDTD, 1 EME MMI1 body, 2 FDTD
    transition, 3 EME MMI2 body, 4 output FDTD, once per port):** section 0 by
    08:50, section 1 (EME) by 10:10, section 2 (FDTD transition,
    `CompoundGeometry`) by 10:23, then a **2h17m** gap to 12:45 for section 3's
    EME solve (a `solve_section_length_sweep` call — by far the single most
    expensive piece of this device), then one more EME
    `CompoundGeometry`-keyed cache write at 13:40:31 — 1m41s before the kill.
    The tail of `reports/11585847_0.out` is mid-way through building the first
    section-4 (output) FDTD grid when SLURM cancelled it. **Net: sections 0-3
    are fully cached; zero of the 3 per-port section-4 solves completed.**
  - **Resubmitted unchanged as 11590663** (same `batch/run_modal_multiplexer_
    examples.batch`, no script edits — nothing touched since 08:41 changes
    `modal_multiplexer_hybrid.py`'s or the solvers' cache fingerprints, so a
    clean disk-cache hit on sections 0-3 is expected). Remaining real work is
    the 3 section-4 FDTD solves; comparable single-section FDTD solves earlier
    in this same log took 5-33 min, so the 5h budget should now be generous
    rather than tight. Watch `reports/11590663_0.out` — if it doesn't reach
    `example_output/modal_multiplexer/modal_multiplexer_hybrid/` output within
    the first ~15 min, the cache did not hit as expected and this entry needs
    another look.
  - **❌ TIMEOUT again, 2026-07-31 18:58, full 5h.** `sacct`: `State=TIMEOUT`,
    `ExitCode=0:0`, same as the first attempt. The "generous rather than tight"
    estimate above was wrong — the 3 section-4 FDTD solves were not a short
    tail. `example_output/modal_multiplexer/modal_multiplexer_hybrid/sections/`
    confirms all 3 per-port section-4 geometry plots exist (`port_top`,
    `port_bot`, `port_mid`), so the disk-cache hit on sections 0-3 worked as
    predicted and this run did get into section 4 — but the log's tail is deep
    inside a per-port eigenmode-overlap monitor computation (repeated `MPB
    solved for frequency_13(...)` lines sweeping many k-points, then starting
    band 14) when SLURM killed it, not writing final S-matrix output. No new
    `example_sim_cache/` entries appear for this device in the run's own
    13:58–18:58 window (the cache dirs that changed in that window all belong
    to other jobs' MMI runs sharing the login node's cache directory), so
    section 4's per-port solves are not yet caching intermediate progress the
    way sections 0-3 apparently do — a plain resubmit will very likely restart
    section 4 from scratch rather than resume mid-monitor-sweep. **Gap 6's
    original ask — get `warn_on_transfer_duplicate` to fire or not fire on
    this device's real interfaces — is still unanswered**, now after two full
    5h attempts (10h of 4-node time total). Worth reshaping before a third
    blind resubmit: either raise the wall-clock limit substantially (this
    device's section-4 monitor cost was underestimated twice now), or split
    section 4's 3 per-port solves into their own array job so a future timeout
    loses at most one port's progress instead of all three. Not resubmitted as
    of this snapshot — see `modal_multiplexer_validation.md` gap 6 for the
    validation-doc-side reasoning.

### 11592462 `[0]` — `beef_dc_examples` (directional_coupler_hybrid) — ❌ TIMEOUT, diagnosed 2026-07-31 late
- **Script:** `batch/run_directional_coupler_examples.batch` → `make
  directional_coupler_hybrid` (array index 0 of the three coupler examples;
  `_1` mode and `_2` fdtd both completed separately, see table above).
- **❌ TIMEOUT 2026-07-31 21:34, full 5h.** `sacct`: `State=TIMEOUT`,
  `ExitCode=0:0`.
- **Diagnosed, full read of `reports/11592462_0.out` (all 512 lines).** The
  two `Initializing structure...` blocks (lines 22 and 330) are **not**
  `M_in0`/`M_in1` (the script's own "2 FDTD solves" cost estimate) — they
  don't match that shape at all. A real `solve_section()` call (confirmed
  against 11591249's successful log) prints ONE small mode-source block
  then goes straight into `on time step` lines; these two blocks instead
  each solve **every wavelength** (11, then a 1-wavelength repeat) × 5 bands
  before anything else — that is `section_mode_neff()`'s own documented
  shape ("one structure, every wavelength... seconds, not minutes"). Cross-
  referencing `directional_coupler_hybrid.py`'s source: `CALIBRATE_BODY =
  True` (hardcoded, not env-gated) runs `section_mode_beta()` first, whose
  *seed* step is exactly one such call (N_MODES=5, all 11 wavelengths) —
  that's block 1 (lines 22-328). The script's separate single-wavelength
  `_mpb = section_mode_neff(...)` diagnostic (used only for a comparison
  print line) is block 2 (lines 330-376, 1 wavelength since
  `wavelengths[_i:_i+1]`) — its first seed value (`frequency_1(1.6129,...)`)
  exactly matches block 1's 6th (middle, λ=1.55) repeat, confirming it's the
  same cross-section at the band-centre wavelength, not a new mystery.
  **Both blocks completed successfully** (block 2 reaches "Dominant
  planewave for band 5" at line 376). **What's missing: the 15 real
  `solve_section()` calibration runs `section_mode_beta()` is supposed to
  do in between** (3 lengths × min(5,5) excited modes) — there is a bare
  `-----------` divider between blocks 1 and 2 with zero evidence any of
  them ran (a real `solve_section()` call always leaves a trace, per the
  11591249 comparison above). Two readings: either all 15 hit the on-disk
  cache (`save_dir=save_dir` is passed) and returned near-instantly, or the
  job never got that far before block 2 started — either way, **the
  documented-cheap `section_mode_neff()` seed/diagnostic calls, not the
  documented-expensive `solve_section()` calibration runs, consumed
  essentially the entire 5h budget.** That contradicts
  `section_mode_beta()`'s own docstring by orders of magnitude at this
  device's production N_MODES=5 / 11-wavelength / `Computational cell is 3 x
  7 x 2.5` configuration — not diagnosed further (whether it's a full-3D-grid
  MPB cost rather than a cheap 2D cross-section, or something else). The
  final 120 `Warning: grid volume` lines (right after block 2 finishes, no
  MPB solves interleaved) most likely mark the start of whatever comes next
  in the script (the real calibration solves and/or `M_body`'s 200-length
  EME sweep) being killed before any could report output — not independently
  explained.
- **Recommendation: do not resubmit this whole script blind a third time.**
  Before resubmitting, benchmark `section_mode_neff()`'s cost in isolation
  at this device's real N_MODES=5/11-wavelength/domain-size configuration
  (a short, cheap diagnostic job, the same pattern as this session's other
  `diagnose_*` scripts) to find out whether that step is really the
  multi-hour cost here — if so, `CALIBRATE_BODY`'s seed computation itself
  needs reshaping (e.g., fewer wavelengths for the seed, since it only needs
  to pick the right 2π branch, not a converged answer) before this job can
  ever finish in one slot. Not resubmitted as of this writing.
- **Benchmarked, 2026-07-31 22:05 — theory REFUTED.**
  `diagnose_dc_calibration_seed_cost.py` (new) measures
  `section_mode_neff()` directly at this exact configuration (N_MODES=5,
  same cross-section, same `Computational cell is 3 x 7 x 2.5`). Login-node
  serial smoke (1 process): **28.7 s/wavelength.** Same script under the
  production resource shape, job **11597622** (2 nodes, 48 ranks, 3
  wavelengths): **44.4 s/wavelength** — MPI rank overhead makes it *slower*
  per-wavelength, as suspected, but only ~1.5×, not the 30-40× that would be
  needed to explain a multi-hour cost. **Extrapolated full cost (11
  wavelengths, production MPI config): ~8 min.** That directly refutes the
  "cheap eigensolve is secretly expensive here" theory above — whatever ate
  the ~4.5h between the two `section_mode_neff`-shaped log blocks in
  11592462_0, it was NOT this call itself.
  **Root cause is back to unknown.** Given seed+comparison together cost
  ~9 min and the 15 real `solve_section()` calibration runs left literally
  zero log trace (no `Initializing structure` between the two blocks),
  those 15 runs most likely all hit the on-disk cache (near-instant, which
  also prints nothing) — meaning **essentially the entire accounted-for
  workload in this script should complete in well under 15 minutes**, yet
  the real job ran 5 hours. Something consumed the other ~4h45m with no
  log signature at all (no print, no crash, no partial structure-init).
  Live open questions, none checked yet: whether disk I/O contention
  (concurrent jobs hitting the same `example_sim_cache`/`save_dir`) stalled
  a cache lookup or write; whether `M_body`'s 200-point
  `solve_section_length_sweep()` (never confirmed to have started) is
  silently expensive in a way that also prints nothing; or whether the job
  was simply stuck/hung (cf. job 11586046's unrelated hang earlier this
  session, also silent). **Next step: resubmit with a shorter time limit
  and watch it interactively** (`squeue`/log tail) rather than diagnosing
  from a dead job's log alone — a live job can be interrupted and inspected
  mid-run, which the two theories above cannot be told apart without.

### 11591249 — `beef_mmi_2x2_engine_agreement` — ❌ TIMEOUT
- **Script:** `batch/verification/run_mmi_2x2_engine_agreement.batch` →
  `examples/verification/python/verify_mmi_2x2_engine_agreement.py` (new).
- **Resources:** 2 nodes × 24 ranks, 5h limit — same shape as this device's
  own production hybrid/FDTD jobs, since this is real FDTD time-domain work,
  not an eigensolve.
- **What it does:** Solves a uniform stretch of the 2×2 MMI body at 4 lengths
  (2/4/6/8 µm — deeper margin than the coupler's minimum two, since this
  body has more room) and checks FDTD and EME agree on each mode's
  accumulated phase via the shared `dbeta_from_phases` helper — the same
  2-mode pair gap 2's own Δβ_01 measures. This is check 1 (engine agreement)
  for the 2×2, the device's last unmeasured check.
- **Why:** answers whether this device's own hybrid-vs-FDTD gap (gap 5:
  0.82% max\|ΔT\|, sign flips across the band) has the coupler's
  Δβ-discretization mechanism underneath it or something else —
  `calibrate_body_neff` has never been tried on this device.
- **❌ TIMEOUT 2026-07-31 19:46, full 5h.** `sacct`: `State=TIMEOUT`,
  `ExitCode=0:0`. Log shows all 4 lengths' first FDTD run completed (`run 0
  finished` appears 4 times, at 29,398 / 30,997 / 33,334 / 30,505
  timesteps), then a 5th run began and was only at timestep ~5,311 (of a
  likely ~30,000-step run, by the pattern of the first four) when killed —
  under a fifth of the way into whatever this run is (a second pass per
  length, or a reference run; the script's own structure was not re-read to
  confirm which). **Gap 4 — this device's last open check — did not close.**
  Net: the job got roughly halfway through its total FDTD work in one 5h
  slot; a straight resubmit (no cache to fall back on here — this script has
  no `save_dir` reuse the way the modal mux's does) would need to redo all of
  it, or the script needs checkpointing / a longer wall-clock budget before
  trying again. Not resubmitted as of this snapshot — see
  `resolved/mmi_2x2_validation.md` gap 4 for the validation-doc-side reasoning.

### 11596845 `[0-3]` + 11596851 — `beef_mmi_2x2_engine_agreement` (reshape of 11591249) — ✅ COMPLETED, GAP 4 CLOSED
- **Scripts:** `batch/verification/run_mmi_2x2_engine_agreement.batch`
  (`--array=0-3`, `-t 480`) and `…_merge.batch`, both →
  `verify_mmi_2x2_engine_agreement.py`. `SLURM_ARRAY_TASK_ID` picks one
  length; the merge step runs with none set and re-reads all four from the
  on-disk cache under `cache_only=True`.
- **Outcome, 2026-08-01:** `_0` (L=2) and `_1` (L=4) **COMPLETED in 36 s
  each** — the predicted cache hit on 11591249's already-computed lengths,
  which is the one thing that attempt did leave behind. `_2` (L=6) 3:35:10,
  `_3` (L=8) **5:07:14**. Merge **11596851** COMPLETED 03:07 in 25:29 on 1
  rank; no `CacheMissError`, so every piece was where it was supposed to be.
- **✅ Result — gap 4 closed:**

  | quantity | value |
  |---|---|
  | mean ratio Δβ_FDTD / Δβ_EME | **0.9953** (0.47% low) |
  | across 1.530–1.565 µm | 0.9962 → 0.9945, **monotone** |
  | worst \|off-diagonal\| of T | 1.21e-05 |
  | worst phase-fit residual | 4.02e-06 rad |

  Against the coupler's ~10% on the identical measurement this is **21×
  tighter**, so `calibrate_body_neff` is **not indicated** for the 2×2 — a
  negative result, and the one gap 4 existed to get. Saved to
  `example_output/verification/mmi_2x2_engine_agreement/engine_agreement.npz`
  (`beef_git` 9e52e9d, res 40).
- **⚠️ Read the log carefully here:** it prints both diagnostics as `0.0000`
  because the script formats them `:.4f`. They are *not* zero and *not*
  defaulted — the `.npz` was opened and the real values (1.21e-05, 4.02e-06)
  read out before this was recorded. An exact-looking zero in a merge step
  that could in principle skip work is worth checking rather than quoting.
- **Why it does not explain gap 5:** the bias here is monotone across the
  band; gap 5's ΔT flips sign across that same band. A monotone Δβ bias
  cannot produce a sign-flipping ΔT, so gap 5's 0.82% comes from elsewhere
  (tapers or the cascade — neither isolated). Argument from signatures, not a
  measurement.
- **Cost datum worth keeping:** 96 / 183 / 215 / 307 min for L = 2/4/6/8 —
  13.4 h serial, i.e. the original single 5 h slot was never going to make it
  and L=8 alone exceeds it. Folded into `resolved/mmi_2x2_validation.md` gap 4 and
  `README.md`.

### 11597177 `[0-2]` + 11597199 — `beef_modal_mux_examples` (reshape of 11590663) — ❌ THIRD TIMEOUT, cause found, NOT resubmitted
- **❌ All three array tasks TIMEOUT** 2026-08-01 02:44 / 02:52 / 03:01, each
  the full 5 h on 4 nodes. **11597199** (merge) never started —
  `DependencyNeverSatisfied`; inert, cancel it.
- **Root cause — the reshape dropped `srun`.** The batch script invokes
  `python modal_multiplexer/python/modal_multiplexer_hybrid.py` directly to
  stop `make` from running all 3 ports in one process, but the Makefile target
  is `$(runner) $(my_python) …` with `runner = srun --mpi=pmi2` on Phoenix
  (`BEEF/examples/Makefile:5`). Losing `make` lost the launcher:
  `Using MPI version 4.1, **1 processes**` in all three new logs vs
  `96 processes` in 11585847_0 and 11590663_0 — **0.27–0.36 s/step against
  0.0064, ~42× slower**, with 96 cores held idle (`cput=20-00:35:12` for
  `walltime=05:00:22`).
- **Second flaw, independent of the first:** the premise that sections 0-3
  were already warm on disk is **false** — every task re-ran the same three
  `Initializing structure`/MPB/time-stepping blocks at near-identical log
  offsets (18/2189/4528, 18/2184/4523, 18/2224/…). Not established whether
  the earlier timeouts never persisted them or the cache key differs between
  96-rank and 1-rank runs. Tasks died at unrelated points (`_0` timestep 289,
  `_1` timestep 2945, `_2` inside MPB band 10), consistent with all three
  grinding through that duplicated prefix.
- **Not resubmitted — the mux is deferred.** Recipe if resumed, in order:
  restore `srun --mpi=pmi2`; verify sections 0-3 actually hit cache (a task
  that skips them should be time-stepping within a minute or two, not an
  hour); then resubmit. **Do not just raise the wall clock** — at 1 rank no
  affordable limit finishes. Full write-up in
  `deferred/modal_multiplexer_validation.md` gap 6.
- **Note the pattern is not at fault:** the same array-plus-`cache_only`-merge
  design submitted in the same window worked first try on the 2×2 (11596845
  above). This broke on the launcher.

### 11586046 — `beef_diagnose_items_11_12` — ❌ CANCELLED (hung), superseded by 11588309
- **Script:** `batch/verification/run_diagnose_items_11_12.batch` → two python
  scripts in sequence, `examples/verification/python/`:
  `diagnose_kpoint_ladder_v2.py` then `diagnose_dup_band_grid.py`.
- **Resources:** 1 node, 8 cpus, serial (no MPI), 240 min budget. Eigensolve
  only — `get_eigenmode()`, no time-stepping, same affordability argument as
  `run_diagnose_mode_correspondence.batch`.
- **Output:** `reports/11586046.out`.
- **Why one job for two items:** both solve the *same* cross-section (modal
  multiplexer interface 0, res 40), so they share setup cost.
- **Item 12 runs first, deliberately** — it is the root-cause item and the
  more expensive of the two, so a wall-clock overrun loses the cheaper result,
  not the more valuable one. Runs a **real** guarded accept/skip/retry search
  (not v1's flat 20-solve sweep) under four seeding policies: A fixed
  (today's code), B naive ladder (raw previous solve — what v1 tested), C step
  ladder (previous − 0.005, aimed *past* the previous band), D accepted ladder
  (last accepted band, so a skip cannot poison the next seed). C and D are
  exactly the two follow-ups item 12 lists as "not yet tried". Reports attempts
  needed, whether any policy gets stuck re-finding one duplicate, and whether
  all four converge on the same final physical mode set.
- **Item 12's open risk it is testing:** v1 showed naive laddering *moves* the
  duplicate rather than removing it (recovers 2 real modes, but re-finds the
  previous band at a new position). C and D are the two candidate fixes.
- **Item 11 second:** replays the duplicate-guard grid comparison at **res 40**
  — the single data point item 11 names as missing before it can close, since
  its clear result so far was measured only at res 20 and a finer mesh packs
  more oscillation into the same core width. `diagnose_dup_band_grid.py`'s
  `RES` is env-overridable (`BEEF_DIAG_RES`, default 20.0) so the res-20 result
  on record stays reproducible as the script's default.
- **❌ Cancelled 11:12 after 2h30. It was HUNG, not slow** — and the 10:53
  revision of this entry got that wrong, calling it a wall-clock budget miss.
  Zero of the four `===== policy =====` blocks ever printed. Reconstructing the
  converged `k_x` per band from the MPB trace explains it: at this cross-section
  bands 1–15 are guided (n_eff 2.7064 → 1.4848, cladding 1.446804) with
  duplicate pairs at 10/11 and 14/15, bands 16–22 are *below* the cladding
  index, and **band 23 has no k that reaches the target frequency at all**.
  `N_TARGET = 20` distinct modes therefore needs band_num 23, and MPB's k-search
  there limit-cycles between k = 0.2346 (f = 0.4775) and k = 1.7480 (f = 1.0283)
  forever — **161 solves on that one band, 59% of the job's entire eigensolve
  work**, against 3–8 solves for each of bands 1–22. No wall clock would have
  been enough, and the 4-task split recommended at 11:03 would have produced
  four hung tasks. Bands 1–22 took ~58 min for 22 attempts, so the original
  ~25 s/solve budget was about right.
- **What it produced anyway, before dying.** The band table above: 15 guided
  band slots at res 40 containing only **13 distinct** modes under today's fixed
  seed, while the device requests `N_MODES = 16`. Read at the time as "the basis
  is unachievable by *any* seeding policy"; the re-run showed a laddered seed
  reaches **14**, so the ceiling is 14 and part of the shortfall is a software
  defect. Folded into
  `resolved/fdtd_vs_EME_mode_validation.md` item 12, with the caveat that whether (10,11)
  is a solver duplicate or a real pair meep cannot separate is still item 1's
  open question.
- **And it opened item 14.** `_iter_distinct_bands()` (`beef/solvers/fdtd.py:1049`)
  walks `band_num` up to `n_modes + 8` and its `RuntimeError` can only fire
  *after* a solve returns, so the same non-convergent eigensolve would hang a
  production `HybridSolver` run with no exception and no timeout. Latent, not
  observed; recorded rather than fixed, per the same
  scope argument gap 6 used for this method.

### 11588309 — `beef_diagnose_items_11_12` (replaces 11586046) — ✅ COMPLETED
- **Submitted** 2026-07-31 11:15, started immediately. Same script as 11586046
  with three changes, all of them consequences of that job's hang:
  - **`diagnose_kpoint_ladder_v2.py` stops at the guided edge.** The search now
    ends when a raw solve returns `n_eff ≤ cladding index` — nothing above is
    guided, and the next band_num is the non-convergent one. A *physical* bound,
    not a timeout: MPB runs in C, so a Python signal cannot interrupt a solve
    already entered.
  - **`N_TARGET` 20 → 16**, what the device actually requests. The 20 was for
    comparability with v1's raw sweep, whose bands 16–20 were unguided anyway.
  - **Item 11 now runs FIRST.** Putting item 12 first to protect the
    higher-value result is exactly what cost both items last time: item 12 hung
    and printed nothing, and item 11 never started. Item 11 is ~50 min and every
    band it solves is known to converge, so it banks a guaranteed result; item 12
    prints per policy, so a late overrun there degrades gracefully.
- **Resources:** 1 node, 8 cpus, serial, **360 min** (was 240). Budget re-derived
  from 11586046's own numbers rather than v1's: ~2.6 min/attempt at res 40,
  ~16–18 attempts per policy × 4 ≈ 3 h for item 12, ~50 min for item 11.
- **Output:** `reports/11588309.out`.
- **Termination verified before submitting**, offline and meep-free: replayed
  11586046's measured n_eff sequence (including all three duplicate pairs)
  through `_guarded_search` with fake `sim`/`mp` objects and the real
  `field_overlap_matrix`. All four policies stop at band 16 on the cladding
  test, none requests the non-convergent band, and policy A reproduces the
  measured 13 accepted / 16 attempts / skips at band 11 and 15.
- **What to watch:** whether any policy reports `stopped because: attempts`
  (would mean duplicates ate the whole +8 headroom production also relies on),
  and whether C or D recovers more than policy A's 13 distinct guided modes.
  **Answered: all of B, C and D recover a 14th (n_eff 1.512589) at the same 16
  attempts; none hit the attempts backstop. See item 12.**
- **✅ Result (11:50, 34:55 — item 11 took 14:56, item 12 19:42). Both items
  resolved, and the reshape worked exactly as intended:** item 11 ran first and
  banked its result, item 12 terminated on the cladding test at band_num 16 under
  every policy rather than walking into the non-convergent region.
  - **Item 11 → CLOSED.** All four diagnostic grids (16×8 shipped through
    128×16) skip the same two band_nums (11, 15), accept the same 16, and give
    identical n_eff. Both skips are exact n_eff ties, so they are confirmed
    genuine without any field sampling. Worst off-diagonal overlap among accepted
    bands on the shipped grid **0.3877** (res 20 was 0.3392) against a 0.9
    threshold — the margin narrowed with resolution as the item predicted, and is
    still less than half way to firing.
  - **Item 12 → decisive.** Policy A (today's fixed seed) finds **13** distinct
    guided modes; B, C and D each find **14**, the extra being n_eff 1.512589,
    *absent from A entirely*. Same 16 attempts, one fewer skip — a strict
    improvement, so the fix direction is settled. C and D returned byte-identical
    results to B, so their extra machinery buys nothing here.
  - **Correction this forced:** the earlier reading that "only 13 distinct guided
    modes exist at res 40, so no seeding policy can help" was inferred from the
    fixed-seed trace — the very enumeration under test — and is wrong. Corrected
    in six places.

### 11588569 — `beef_audit_mmi_2x2_gap6` — ✅ COMPLETED
- **Script:** `batch/verification/run_audit_mmi_2x2_gap6.batch` (new 2026-07-31)
  → `examples/verification/python/audit_interface_modes.py mmi_2x2`.
- **Resources:** 1 node, 8 cpus, serial, 120 min budget. Used **20:18**.
- **Output:** `reports/11588569.out`, plus
  `example_output/verification/interface_mode_audit/`.
- **Why its own batch file** rather than the existing
  `run_audit_interface_modes.batch`: that one runs all three devices *with*
  `--convergence` (an EME mesh sweep, ~4× the mode solves). Gap 6 asks only
  whether the identity permutation and worst |Δn_eff| still hold at res 40 with
  the guards live — the mesh sweep is gap 2's separate question.
- **Also fixed on the way in:** the audit's own `mmi_2x2` config claimed to be
  "pinned to `mmi_2x2_hybrid.py`" while still requesting `n_modes=15`. The device
  moved to 14 on 2026-07-30 (item 9), so the audit had been describing a request
  the device no longer makes. Corrected to 14.
- **✅ Result — passed unchanged, and nothing fired.** 14/14 guided (band 15, the
  sub-cladding one, is gone as intended), correspondence identity `[0…13]`,
  **worst guided |Δn_eff| 0.0007** — the same number as the 2026-07-26 audit at
  15 modes. **Zero duplicate-guard firings** in the log from either engine, which
  answers the specific concern gap 6 raised: this is the repo's highest-mode-count
  SOI_220 device on the 2.8 µm body that `resolved/fdtd_vs_EME_mode_validation.md` item 11
  names as the guard's least-margin case, and it had nothing to fire on. "The
  cleanest device audited" is a current claim again.

### 11586098 — `beef_diagnose_item_13` — ✅ COMPLETED
- **Script:** `batch/verification/run_diagnose_guard_over_request.batch` →
  `examples/verification/python/diagnose_guard_over_request.py`.
- **Resources:** 1 node, 8 cpus, serial (no MPI), 360 min budget. meow/MPB
  cross-section solves only — no propagation, no time-stepping.
- **Output:** `reports/11586098.out`.
- **Cross-section:** modal multiplexer interface 0 again — same as 11586046 and
  the item-0 duplicate — but at **res 60**, per item 13: above the shipped 40,
  and where 11575516 found meow converged (16/16 guided, plateau 60→80), so a
  difference seen here is not an unconverged mesh talking.
- **Why three arms, not "guard on vs off":** item 13's own point is that
  flipping the flag is *not* a clean A/B — with the guard on,
  `compute_modes()` over-requests `num_modes + _DUP_MODE_MAX_EXTRA` (16 → 24)
  in one `compute_modes_meep()` call, so the flag changes both the
  eigenproblem and whether filtering happens. The arms separate them:
  **A** guard off / 16 (baseline), **B** guard on / 16 (production default),
  **C** guard off / 24 (B's eigenproblem, unfiltered).
- **What each comparison answers:**
  - **C vs A** — does over-requesting *alone* move the low bands? This is the
    artifact question. Non-zero ⇒ every guarded EME number in the repo carries
    an over-request artifact, whether or not a duplicate was ever filtered.
    Inert ⇒ the guard is a pure filter and `_DUP_MODE_MAX_EXTRA` becomes a free
    cost knob.
  - **B vs C** — what the filter actually drops: which of C's 24 the guard
    rejected as duplicates, and what it reached further up the spectrum to get.
- **Cost is half the point:** wall time is recorded per arm, so if C vs A comes
  back inert the `C/A` ratio is the whole argument for tuning
  `_DUP_MODE_MAX_EXTRA` down — pure cost, no effect on the answer. Item 13
  notes job 11575516's 3h36m was dominated by exactly this over-request.
- **Result (10:50, 1:54:50) — inert, and expensive.** **C vs A: max |Δn_eff| =
  0.000e+00** on all 16 low bands, so the over-request does not perturb the
  eigenproblem at all and the guard is a pure filter. **B vs C:** B's 16 are
  exactly C's first 16 (worst residual 0.000e+00, no index skipped) — no
  duplicate fired at res 60, itself a new data point given the res-40 firings.
  **Cost:** A 934.0 s, B 2948.5 s, C 2997.7 s ⇒ `B/A = 3.16×`, `C/A = 3.21×`;
  the two ratios agreeing confirms the expense is the eigenproblem, not the
  filtering. Item 13 closed; folded into `resolved/fdtd_vs_EME_mode_validation.md` item
  13, which also records why `_DUP_MODE_MAX_EXTRA` was **not** cut yet.

### 11586986 — `beef_dispersive_convergence`
- **Script:** `batch/verification/run_dispersive_convergence.batch` →
  `examples/verification/python/verify_item_8_dispersive_convergence.py`
- **Resources:** 1 node, 8 cpus, serial, 480 min budget. Eigensolve only.
- **Output:** `reports/11586986.out`.
- **What it does:** Sweeps res 20→120 on the coupler body cross-section under
  both material settings, reporting Δn_eff of the even/odd supermode pair and
  the *step* between consecutive resolutions.
- **Why:** `material_model_validation.md` work item 4 asks for a re-run of the
  full resolution sweeps, because smoothing moves the mesh-convergence point
  and item 11's res-80 pin is not inherited once dispersion is on. That sweep
  is 11 array tasks × 2 nodes × 24 ranks at up to 24 h. But the pin was
  justified on **Δβ**, a mode-solve quantity — so this answers "does the
  convergence point move" on one serial task first, to decide whether the
  node-days are warranted and over what range.
- ~~**Caveat, recorded 2026-07-31:** submitted *before* λ₀ became band-derived
  (commit `4296ee7`), so its **dispersive column is stale** — it used the old
  fixed λ₀ = 2.0.~~ **Wrong diagnosis, corrected 11:45.** λ₀ cannot matter here,
  because the dispersive medium was never built at all.
- **⚠️ This job measured nothing.** Its "dispersive" column is **byte-identical**
  to its non-dispersive one — same n_eff, Δn_eff and Δβ to six decimals at res
  20/30/40/60/80. `FDTDSolver._build_section_mode_sim()` built media with
  `si_medium(wl_center)` directly, bypassing `_media()`, so it never read
  `self.dispersive`; `section_mode_neff()` takes a whole band and loops it on
  that one frozen-centre structure. Confirmed three ways: this log, a 2.7 s
  direct check (max |Δn_eff| = 0.000e+00 with `resolve_dispersive → True`), and
  the media themselves (ε 12.0630 vs 12.0946 at 1.60 µm — the mode solve sees
  the frozen one under both settings).
- **Fixed the same day**; a real solve now moves 5.29e-03 across the band where
  it moved by exactly 0. Full reading in `material_model_validation.md` work
  item 4, which is back to **not started**. The non-dispersive column (res
  20→120, all seven rows) is valid and is preserved in that entry, so nothing is
  lost by superseding this job.

### 11588997 — `beef_dispersive_convergence` (replaces 11586986) — ✅ COMPLETED
- **Same script, post-fix.** `reports/11588997.out`, 1:56:53.
- **Result:** at band centre (1.55 µm) dispersion changes nothing — identical
  convergence resolution to the non-dispersive column at every threshold
  (1e-3/3e-4/1e-4/3e-5 all first clear at res 60/100/100/100 either way). At
  the band edge (1.60 µm, where dispersion is not near-inert by construction),
  the 1e-3 threshold moves from res 60 to res 80 under dispersion — but item
  11 already pins production's EME resolution at 80, so the dispersive case's
  actual requirement lands exactly on what already ships; nothing currently
  deployed was ever run at res 60. The tightest threshold (3e-5) is the one
  gap: non-dispersive clears it everywhere by res 100, dispersive never drops
  under it at the edge inside the tested range (20–120) — though the raw
  dispersive-vs-non-dispersive Δn_eff difference oscillates at ±(3–5)e-4 with
  no resolution trend, reading like noise floor rather than a growing effect.
- **Verdict: the proposed node-day full-device sweep (11 array tasks × 2 nodes
  × 24 ranks, up to 24h) is not warranted** at the 1e-3-class tolerance item 11
  actually used to set the shipped resolution. Full table in
  `material_model_validation.md` work item 4, now closed for the mode-solve
  half.

### 11589770 — `beef_mmi_1x2_hybrid_calibrated`
- **Script:** `batch/run_mmi_1x2_hybrid_calibrated.batch` (new) →
  `MMI1X2_CALIBRATE_BODY=1 make mmi_1x2_hybrid`.
- **Resources:** 2 nodes × 24 ranks, 300 min — same shape as the uncalibrated
  hybrid array task (11585845_2, 1:46:59), since `calibrate_body_neff` changes
  how the body section's phase is propagated, not how many FDTD sections run.
- **Why:** `resolved/mmi_1x2_validation.md` gap 5 measured the hybrid 9.11% **under**
  the full-device FDTD reference, every ΔT negative — the Δβ-discretization
  signature `resolved/directional_coupler_validation.md` items 6/14 describe
  on the coupler, where `calibrate_body_neff=True` took it from ~5% to 0.11%.
  Never tried on this device. `mmi_1x2_hybrid.py` gained a
  `MMI1X2_CALIBRATE_BODY` env flag (redirects output to
  `mmi_1x2_hybrid_calibrated/` so the uncalibrated baseline survives; same
  pattern as `directional_coupler_hybrid_stub_sweep.py`'s
  `DC_CALIBRATE_BODY`).

### 11589779 — `beef_mmi_1x2_gap5_calibrated`
- **Script:** `batch/verification/run_verify_mmi_1x2_calibrated.batch` (new) →
  `MMI1X2_HYBRID_CALIBRATED=1 python verify_mmi_1x2_hybrid_vs_fdtd.py`.
- **Resources:** 1 node, 1 cpu, 20 min. `--dependency=afterany:11589770`.
- **What it does:** Same comparison as 11585850, pointed at the calibrated
  hybrid's SMatrix instead. `verify_mmi_1x2_hybrid_vs_fdtd.py` gained a
  `MMI1X2_HYBRID_CALIBRATED` env override for this; `FDTD_PKL` (the shared
  full-device reference) is deliberately **not** redirected.
- **Why:** answers whether `calibrate_body_neff` closes this device's 9.11%
  the way it closed the coupler's ~5%. Compare its worst max|ΔT| against
  11585850's 0.0359 once both have run.

### 11590785 — `beef_mmi_2x2_gap2_dbeta`
- **Script:** `batch/verification/run_mmi_2x2_gap2_dbeta.batch` (new) →
  `examples/verification/python/verify_mmi_2x2_gap2_dbeta.py` (new).
- **Resources:** 1 node, 8 cpus, serial, 180 min budget. Eigensolve only —
  `EMESolver.compute_modes()` at `num_modes=2`, no FDTD, no time-stepping.
- **Why:** `resolved/mmi_2x2_validation.md` gap 2's existing measurement (job 11500536)
  shows the MMI body's fundamental n_eff still moving at res 100, but that
  doc explicitly declines to conclude anything from an arbitrary quantity —
  same "sweep the wrong thing, hide the real error" lesson item 11 taught on
  the coupler (there: T_cross vs Δβ). For a 2×2 MMI the physically load-bearing
  quantity is Δβ_01 = 2π/λ·(n_eff[0]−n_eff[1]) between the two lowest guided
  modes — the beat length `mm_length=28.0` was chosen against — not either
  mode's own n_eff. This sweeps res 40/60/80/100/120 measuring Δβ_01 directly,
  the same pattern item 11 used for the coupler's even/odd pair.
- **Cost note:** a login-node smoke check (`OMP_NUM_THREADS=1`, res 40 only)
  measured 205.8 s for one point — slower than expected for a 2-mode solve,
  hence the batch submission (8 threads) rather than running the full sweep
  on the login node. 180 min budgets generously for the higher-resolution
  points, whose mesh cost does not scale linearly with `res`.
- **✅ Result (1:05:48).** Δβ_01 steps: −0.00092/−0.00049/−0.00020/−0.00008
  rad/µm at 40→60/60→80/80→100/100→120 — 0.56%/0.30%/0.12%/0.049% relative.
  **Gap 2 closes at res 100** for this device: the physically load-bearing
  quantity is converged to well under 0.1% there, regardless of absolute
  n_eff's own slower-looking tail. One correction made to this session's own
  earlier partial reading: comparing Δβ_01's raw step size directly against
  n_eff[0]'s raw step (different units/magnitudes, ~0.16 vs ~2.8) and calling
  it "an order of magnitude tighter" was a mixed-units comparison — read
  relative to each quantity's own value, the two converge at a *similar*
  rate; what actually matters is that Δβ_01 itself is small enough to trust
  by res 100, not how its convergence rate compares to another quantity's.
  Full table in `resolved/mmi_2x2_validation.md` gap 2.

### 11589771 / 11589899 / 11589936 — `beef_modal_mux_fdtd` — ❌ BLOCKED, do not resubmit
- **Script:** `batch/run_modal_multiplexer_fdtd.batch` (new) → `make
  modal_multiplexer_fdtd` → the new (blocked) `modal_multiplexer_fdtd.py`.
- **Why attempted:** `modal_multiplexer_validation.md` gap 3 — this device is
  the only one in the repo with no full-device FDTD reference at all. The
  script mirrors `mmi_2x2_fdtd.py`'s shape: the same `CompoundGeometry` the
  hybrid script builds, one single-domain `FDTDSolver.solve()` per input mode
  (TE0, TE1) instead of five cascaded sections. Resolution deliberately
  reduced to **20** (env-overridable `MUX_FDTD_RES`) per gap 3's own
  recommendation to scope cheap first — this device is the longest
  single-domain run in the repo (62.16 µm) and 3-D FDTD cost is steep in
  resolution.
- **Attempt 1 (11589771, 22 s):** crashed in `plot_fdtd_geometry()` →
  `CompoundGeometry.source_aperture_at()` → `_locate()`, which raised for any
  x outside `[-total_length/2, total_length/2]`. `FDTDSolver.solve()`'s own
  source-position query (`src_x = -Sx/2 + dpml + source_margin`, landing in
  the solver's external PML-clearance stub) sits at exactly this kind of x, so
  the real simulation would have hit the identical crash had the plot not
  caught it first — this was going to fail regardless, just later and after a
  full node allocation.
  - **Real bug, fixed:** every leaf geometry (`MMI.width_at()`, etc.) already
    clamps to its terminal value for arbitrarily large overshoot past its own
    taper — that's why plain-MMI full-device scripts (`mmi_1x2_fdtd.py`,
    `mmi_2x2_fdtd.py`) never hit this. `CompoundGeometry._locate()` was the
    only one that hard-failed, purely because delegating to a segment needs a
    segment to exist "out there". `beef/geometry/compound.py`:
    `_locate(x, tol, clamp=False)` gained an opt-in `clamp=True` that
    delegates out-of-bounds x to the nearest end segment instead of raising —
    wired into `source_aperture_at()` only; `width_at()`/`center_at()`/
    `section_prisms()`/etc. keep raising, since an out-of-bounds cut *there*
    is a genuine caller error building physical geometry, not a solver's
    legitimate external-stub query. 5 new tests,
    `tests/test_geometry.py::TestCompoundGeometrySourceApertureClamps`,
    mutation-checked (three fail against the pre-fix code, two pass either
    way). Fast suite 482 passed, 1 skipped.
- **Attempt 2 (11589899, 15 s):** cleared the plot crash, hit a second one
  inside `solver.solve()`'s cache-fingerprinting path:
  `geometry.polygons(stub_length=...)`, and `CompoundGeometry.polygons
  (stub_length>0)` raises `NotImplementedError` by design (its own docstring:
  "no generic way to know per-port stub width for an arbitrary segment
  type"). **Misdiagnosed at the time** as a caching-only issue — dropping
  `save_dir` (disabling the cache) does skip *that* call site.
- **Attempt 3 (11589936, 13 s), same script with `save_dir` dropped:** still
  crashed, same `NotImplementedError`, from a **second, unconditional** call
  site: `solve()` → `to_meep_prisms()` → `geometry.polygons(stub_length=...)`,
  which runs regardless of caching — it is how the geometry is actually built
  for the sim. **This is the real blocker.** `FDTDSolver.solve()` (the full
  single-domain path) cannot run on a `CompoundGeometry` top-level device at
  all today; only `solve_section()` (what `HybridSolver` always uses,
  including this same device's own hybrid script) avoids it, via
  `section_gds()`/`to_gds_eme()` instead of `polygons()`.
- **Not fixed.** Unlike the aperture clamp, this needs a real design decision
  (which segment's stub-port width wins at an ambiguous end — e.g. a
  multi-port `WaveguideRun`), not a mechanical patch, and the author who wrote
  `polygons()`'s sibling stub-handling explicitly declined to guess at it.
  Same caution bar applied to item 14 and the k-seed guard fix this session:
  not doing production-adjacent, under-verified geometry changes blind. Filed
  as `notes.md` Geometry item 8. `modal_multiplexer_fdtd.py` and its batch
  script are both marked BLOCKED at the top — do not resubmit until this is
  actually implemented.
- **Cost:** three attempts, ~15 s / minimal node-seconds each. No wasted
  budget, no partial/corrupt output (only `geometry.png` was ever written;
  crashes happen before any `smat.save()`).

## Not listed here

An `elevator_*` array job (11585814 at the 07-31 snapshot; the 07-30
snapshot's 11568634 / 11569468 have since cleared) is also queued, under
`/storage/scratch1/8/jhiesener3/esc/elevator_sim_tower_from_gds/` — a
different project entirely, unrelated to BEEF or anything in this
directory. Omitted above; see that project's own tree if they matter.

A second set of these same three device jobs was briefly submitted on
2026-07-31 08:26 to `cpu-medium` at `--mem-per-cpu=16G` (11585878/79/80,
`sbatch` overrides only, no script changes) and **cancelled before any of
them started** — no output, no billing, nothing to fold in. Recorded only so
the gap in the job-ID sequence isn't mistaken for a lost run.
