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

Only the expensive primitives are stored: FDTD runs and eigenmode solves.
Derived results are rebuilt on demand and deliberately not persisted —
`solve_section_length_sweep()` stores the per-mode `neff` its transfer matrices
follow from, and `solve_length_sweep()` stores the cascade factorization
(interface S-matrices plus per-cell `neff`). Both are independent of the length
grid, so one entry serves every sweep of a device instead of one entry per
grid, and neither is written unless replaying it reproduces the run that
produced it.

That is why this directory is a few MB rather than a few hundred: the derived
length sweep it used to hold was 342 MB on its own, 97.7% of the cache.

## What runs from this cache, and what does not

Run `make <target>` for any example; anything already cached completes without
simulating. An entry that is missing or stale simply costs a real run, and the
examples are cheap to re-run as a set — that is the intended way to refresh
this directory after a library change that moves a cache key.

To see what an example would need without simulating anything, patch
`beef.cache.store.lookup` to raise on a miss and run the script: it stops at
the first entry it cannot serve.
