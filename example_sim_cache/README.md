# example_sim_cache

Cached solver results that let the examples in `examples/` run without
re-simulating. Every example's `save_dir` points here.

A cache entry is one directory named
`<Solver>_<method>_<GeometryClass>_<hash>`, holding `config.json` (the full
call fingerprint — geometry polygons, wavelengths, resolved material indices,
solver settings), `result.pkl`, and a `geometry.png`. `index.jsonl` is a
grep-friendly manifest, one line per entry. A lookup hits only when every
fingerprinted field matches, so changing a wavelength grid, a resolution, or
the mode count produces a miss and a fresh simulation rather than a wrong
answer.

The 52 entries here are exactly what the current examples ask for — nothing
else. They were selected by probing each example with the cache lookup patched
to record what it resolves, so this is a measured set, not a guess.

## What runs from this cache, and what does not

Ten examples run start to finish with no simulation:

    mmi_2x2          eme, fdtd, hybrid
    mmi_2x2_tapered  eme, fdtd
    mmi_4x4_90hybrid eme, hybrid
    mmi_1x4          eme
    directional_coupler  hybrid
    ring_resonator   hybrid      (RING_TEST_MODE=0, the production path)

The other eleven need a real run before they will complete. See the repository
`CHANGELOG.md` entry for 2026-09-18 for the per-example reason; in short, six
are stale only because the cache fingerprint gained fields
(`geometry_mode`, `design_grid_resolution_factor`, the resolved-material
block, `minimum_run_time`) and five differ in real physics (a pinned monitor
plane, a resolution change, a geometry change, a mode-count change).

## One entry is not in git

`EMESolver_solve_section_length_sweep_MMI_58b3e3c249c841a7/result.pkl` is
327 MB — past GitHub's 100 MB per-file limit — so it is excluded by
`.gitignore` and will be absent from a fresh clone. It is the 4x4 90-degree
hybrid's body length sweep, and `mmi_4x4_90hybrid_hybrid.py` is the only
example that needs it. Without it that example re-runs the sweep; every other
example is unaffected. Track it with Git LFS if it needs to travel.
