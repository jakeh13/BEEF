# Hybrid body field plotting — handoff notes for the 2×2 tapered MMI and the 4×4 90° hybrid

Written 2026-08-31 alongside
`hybrid_field_extraction_validation.md` step 1. The two sections below are
meant to be pasted into whichever session owns each device.

---

## PASTE THIS INTO THE 2×2 TAPERED-MMI SESSION

> **New shared capability: hybrid body field reconstruction.**
> `beef.plot_hybrid_body_fields()` (new, `beef/plot/hybrid_body.py`) reconstructs
> and plots a `HybridSolver` cascade's **EME body section** field in the
> **device frame**. It is standalone — no FDTD field capture, no cascade cache,
> one mode solve (~30–110 s measured). Validated on four devices: 1×2 MMI, 1×4
> MMI (SiN), and both directional couplers + the HCG-FP coupler body. Example
> scripts: `examples/verification/python/plot_mmi1x2_hybrid_body_selfimaging.py`
> (the reference case, with regression assertions) and
> `plot_dc_hybrid_body_fields.py`.
>
> **⚠ READ THIS BEFORE TRYING IT ON A TAPERED/WAISTED BODY — the helper as
> written does NOT apply to you, and the reason is physics, not a missing
> feature.**
>
> The reconstruction assumes a **uniform cross-section** over the body span. On
> that assumption two things hold, and both are load-bearing:
>
> 1. **One cross-section's mode basis serves both ports**, so a single-cell
>    solve describes the whole section.
> 2. With `calibrate_body_neff=True` the body transfer matrix is a **pure
>    diagonal** `exp(2j·π·n_eff/λ·L)` — each mode accumulates only phase, no
>    mode mixing — so the field is `Σ_k c_k φ_k(y) exp(i·2π·n_eff_k·x/λ)`. That
>    sum is the entire reconstruction.
>
> **A tapered body breaks both.** `beef` already enforces this, in two places,
> deliberately: `assert_uniform_cross_section()` is called from
> `EMESolver.solve_section()`'s override path *and* from
> `HybridSolver._calibrated_neff()` (`hybrid.py:1534`). Its own message says
> why: *"a single per-mode beta cannot describe a body whose transfer matrix
> depends on the mode FIELDS at every internal interface."* So on a tapered
> body you will get a hard assertion, not a subtly wrong plot — which is the
> right outcome, and it means you cannot simply point the helper at your span.
>
> **What you need instead.** A tapered body goes through
> `solve_section()`'s **real multi-cell** branch (`num_cells > 1`): one
> cross-section per cell, real interface overlaps between successive
> cross-sections, and a **non-diagonal** transfer matrix. Field reconstruction
> there is genuinely a different calculation, and the good news is that the
> library already has the piece that does it — **`plot_eme_fields()` +
> `mw.propagate_modes()` handle multi-cell propagation natively**, which my
> analytic single-cell sum does not. So for the tapered 2×2 the cheaper path is
> almost certainly:
>
> - build the cells for your body span the way
>   `EMESolver.solve_section()`'s non-override branch does — `section_gds()` →
>   `extrude_polygons()` → `_build_mesh()` → `_flat_cell_lengths(region_lengths=[L])`
>   → `mw.create_cells()`, with **your real `num_cells`**, not 1;
> - solve modes per cell (`compute_modes(cs, n_modes)` for each cross-section);
> - hand the cell/mode lists to `plot_eme_fields()`.
>
> That is what `plot_mmi1x2_hybrid_body_fields.py` does, minus the single-cell
> shortcut. Note also `region_cell_lengths` (added 2026-08-29) if your waist
> needs cell boundaries an even division cannot reach.
>
> **Two things worth stealing from the shared helper even so**, because both
> were real bugs found the hard way:
>
> - **`plot_eme_fields()` labels in meow's frame, not the device's.** It passes
>   its component string straight to `mw.propagate_modes()` with no remap, but
>   meow's labels are projections onto meow's own local axes:
>   `device Ex = meow Ez`, `device Ey = meow Ex`, `device Ez = meow Ey` (see
>   `stack_meow_mode_fields()`'s docstring). The panel titled `|Ex|` is really
>   the device's **in-plane transverse** field, and the `x [µm]` axis is device
>   y. Affects interpretation at all existing call sites.
> - **Any modal projection must be a 2-D area integral, never a 1-D cut through
>   the core.** Measured on the 1×2: a 1-D `Ey`-slice dot product returned total
>   captured power **1.5857** (impossible) and put **0.5876** on an odd mode a
>   centred launch cannot excite. The 2-D integral with `dy·dz` weights gave
>   0.9970 and 1e-11. **The parity zeros survived both versions**, so a parity
>   check alone does not catch this — power conservation is what exposes it.
>   Assert `captured_power <= 1`.

---

## PASTE THIS INTO THE 4×4 90° HYBRID SESSION

> **New shared capability: hybrid body field reconstruction.**
> `beef.plot_hybrid_body_fields()` (new, `beef/plot/hybrid_body.py`) reconstructs
> and plots a `HybridSolver` cascade's **EME body section** field in the
> **device frame**, standalone — no FDTD field capture (that is still unbuilt;
> `fdtd.py` has zero `add_dft_fields`/`get_array` calls), no cascade cache, one
> mode solve. Validated on four devices (1×2, 1×4, both DCs, HCG-FP). Closest
> template for you: `examples/verification/python/plot_dc_hybrid_body_fields.py`
> — it has the **off-axis launch** pattern you need.
>
> **This should apply to you directly**, since the 4×4's body is uniform (the
> BO campaign sweeps `mm_width`/`y_gap`/`mm_length`, not a taper). Minimal
> version:
>
> ```python
> from beef import MMI, EMESolver, SOI_220, plot_hybrid_body_fields
>
> half_body = mmi.mm_length / 2
> x0, x1 = -half_body + BODY_STUB, +half_body - BODY_STUB   # body is k=1 of 3,
>                                                           # so stubs are 0/0
> solver = EMESolver(  # your cascade's OWN body_solver config, verbatim
>     num_modes=[N_MODES], num_cells=[1], lengths=[1.0],
>     oxide_width=..., oxide_height=..., res=..., material_stack=SOI_220,
>     mode_solver="meep",
> )
> res = plot_hybrid_body_fields(
>     mmi, solver, x0, x1, WL_PLOT, N_MODES, output_dir, tag="...",
>     excitation="gaussian",
>     excitation_center=mmi.input_ys[k],   # OFF-AXIS: one of four inputs
>     excitation_fwhm=TAPER_WIDTH,
>     wavelengths=<your full production sweep>,
>     port_ys=mmi.output_ys, body_half_width=mmi.mm_width / 2,
>     device_label="4×4 90° hybrid", total_length=mmi.total_length,
> )
> assert res["captured_power"] <= 1.0
> ```
>
> **Five device-specific things to get right:**
>
> 1. **The launch is OFF-AXIS.** With four inputs, `excitation_center` must be
>    the driven port's own y (`mmi.input_ys[k]`), not 0.0. A centred launch
>    would silently excite only the even body modes and produce a symmetric
>    picture that looks plausible and is wrong. Also: `waveguide_offset` is the
>    port centre-to-centre **spacing**, so `_port_ys(n, s) = s·(i − (n−1)/2)` —
>    read `mmi.input_ys`/`output_ys` off the geometry rather than computing it.
>    (`mmi_1x2_hybrid.py`'s own comment claiming `output_ys = [-wo, +wo]` is
>    wrong by 2×; don't inherit that.)
>
> 2. **⚠ |E|² LOSES EXACTLY WHAT YOUR DEVICE IS FOR — so use `port_phases()`,
>    not the intensity map.** A 90° hybrid's whole function is the **relative
>    phase** between outputs, the quantity your own campaign item 1 measured
>    pure-EME getting wrong by 44–62°. An intensity map cannot validate that at
>    all. This gap is now filled: `plot_hybrid_body_fields()` returns the
>    **complex** field as `res["E"]`, and **`beef.port_phases()`** reads the
>    relative phase at each port off it:
>
>    ```python
>    from beef import port_phases
>    p = port_phases(res, mmi.output_ys, reference=0)
>    print(p["phase"])       # degrees, relative to port `reference`
>    print(p["amplitude"])   # |E| at each port — READ THIS FIRST
>    ```
>
>    Two things to respect when using it:
>    - **Check `amplitude` before believing any `phase`.** The phase of a
>      vanishing field is numerically arbitrary — the same trap
>      `mode-signs-gate-on-strong-modes` records, where a 1e-06 mode disagreed
>      by 280%. A port carrying no power has no meaningful phase.
>    - **Only relative phase means anything.** The launch carries an arbitrary
>      global phase, which is why the helper subtracts a reference port rather
>      than reporting absolutes.
>
>    Validated on two textbook predictions: a symmetric 1×2 splitter's two
>    outputs come back in phase, and a directional coupler's cross port comes
>    back near quadrature from its driven port. Neither is a 90° hybrid, so
>    **your device is the first real test of this helper** — treat a surprising
>    number as possibly the helper's fault, not only the device's.
>
>    Standing caveat still applies: these are the BODY's own phase accumulation
>    from an approximate launch, not the cascade's port-to-port phases.
>
> 3. **Confirm which component actually dominates before trusting the default.**
>    The helper defaults to `component="Ey"` (device in-plane transverse), which
>    is right for the TE-like modes on the devices validated so far — measured
>    ~330× larger than `Ez` on the 1×2. But meow's `te_fraction` is already
>    known not to match this repo's TE convention on your device, so check the
>    component magnitudes on your own cross-section rather than assuming.
>
> 4. **Your best BO point moves.** The plan doc cites W=6.50 / y_gap=0.51 /
>    mm_length=91.00 (round 5) but the campaign was still running on 2026-08-31
>    (`beef_bo_W6.00_y0.51`, `beef_bo_W6.50_y0.56`, `beef_bo_W5.00_y0.43` in
>    queue). Re-read the campaign's own latest result at plotting time, and put
>    the geometry constants in the filename `tag` so two runs at different
>    points cannot overwrite each other.
>
> 5. **Cost.** One cross-section, `N_MODES` modes, one wavelength — not a
>    cascade run. For scale, the 1×4's 12 modes on a 12 µm-wide body at res 30
>    took ~3 min; your ~15–16 modes on a wider body will be more but is still a
>    single mode solve, not a multi-hour rerun. Do NOT let it re-enter the
>    cascade.
>
> **Two general traps, both real bugs found building this:**
>
> - **Any modal projection must be a 2-D area integral, not a 1-D cut through
>   the core.** A 1-D slice gave captured power **1.5857** (impossible) and
>   **0.5876** on a parity-forbidden mode; the 2-D `dy·dz` integral gave 0.9970
>   and 1e-11. The parity zeros survived *both*, so parity alone does not catch
>   it — assert `captured_power <= 1`.
> - **`plot_eme_fields()` labels in meow's frame** (`device Ex = meow Ez`,
>   `device Ey = meow Ex`, `device Ez = meow Ey`). Its `|Ex|` panel is the
>   device's in-plane transverse field and its `x` axis is device y. The new
>   helper does the remap; the old one does not. This matters for you because
>   `plot_mmi4x4_90hybrid_eme_fields.py` is one of the affected call sites.
>
> **What the helper cannot give you**, so you do not wait on it: the body is
> only the middle section. Your input/output taper sections are FDTD and carry
> no field data anywhere in the library today — that is step 2 of
> `hybrid_field_extraction_validation.md` and it is unbuilt. The body panel is
> genuinely useful on its own (mode content, beating, where the self-image
> lands) but it is not a whole-device map, and the excitation entering it is an
> approximation rather than the cascade's own interface field.
