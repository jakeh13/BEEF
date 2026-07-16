# FECES

**F**DTD/**E**ME **C**ascaded **E**lectromagnetic **S**olver — for photonic MMI/coupler devices.

FECES is a hybrid-methodology orchestration layer on top of [meep](https://meep.readthedocs.io/) (FDTD),
[meow](https://github.com/flaport/meow)/MPB (eigenmode expansion), and optionally
[tidy3d](https://www.flexcompute.com/tidy3d/)'s mode solver — plus a small parametrized photonic-device
geometry library. It is not a competing raw solver; meep/meow/tidy3d do the actual field/mode math.
FECES's real contribution is `HybridSolver`: a generic N-section FDTD/EME cascade (any solver type, any
order, any count ≥ 2) with reciprocity-based excitation (minimum simulation count) and a shared
`section_prisms/section_gds/width_at/center_at` interface that lets any geometry slot into that cascade.

## Install

```bash
pip install -e .            # core (numpy, matplotlib)
pip install -e ".[eme]"     # + meow-sim, for EMESolver
pip install -e ".[fdtd]"    # + meep, for FDTDSolver / HybridSolver
pip install -e ".[all]"     # both
```

Requires Python ≥3.10.

## Public API

```python
from feces import (
    # geometry
    MMI, MMIPort, BendGeometry, DirectionalCoupler, Taper, SBend,
    CompoundGeometry, WaveguideRun,
    # materials
    MaterialStack, SOI_160, SOI_220,
    # solvers
    EMESolver, FDTDSolver, HybridSolver,
    # results
    SMatrix,
    # sweep
    ParallelSweep,
    # plot
    plot_eme_fields, plot_eme_geometry, plot_fdtd_geometry, plot_hybrid_geometry,
    plot_performance, plot_splitting_ratio,
    # utilities
    extrude_polygons,
)
```

- **Geometry** (`feces/geometry/`) — `MMI` (N×M splitter/combiner, per-port `MMIPort` taper specs),
  `DirectionalCoupler` (2-waveguide coupler w/ S-bend I/O), standalone `Taper`/`SBend` building blocks,
  `CompoundGeometry`/`WaveguideRun` (concatenate several geometries into one cascade), legacy `BendGeometry`
  (90° bend for non-collinear FDTD I/O).
- **Materials** (`feces/materials/`) — `MaterialStack` (Si/SiO2 layer stack), plus pre-built `SOI_160` /
  `SOI_220` stacks.
- **Solvers** (`feces/solvers/`) — `EMESolver` (meow-based EME), `FDTDSolver` (meep-based FDTD, 2D or 3D),
  `HybridSolver` (chains FDTD and EME sections together — the core FECES method).
- **Results** (`feces/results/`) — `SMatrix`, a unified S-parameter result type with transmission/insertion-loss
  helpers and a `check_passive()` sanity check.
- **Sweep** (`feces/sweep/`) — `ParallelSweep`, a local `multiprocessing`-based parameter grid sweep.
- **Plot** (`feces/plot/`) — geometry, field, and performance (transmission/insertion-loss/splitting-ratio) plots.

## Quick start

```python
import numpy as np
from feces import MMI, EMESolver, SOI_160, plot_performance

mmi = MMI(
    n_inputs=1, n_outputs=2,
    mm_length=4.5, mm_width=2.0,
    waveguide_width=0.35, taper_length=1.0, taper_width=0.8,
    waveguide_offset=1.0, waveguide_height=0.16,
)

solver = EMESolver(
    num_modes=[1, 5, 12, 5, 1],
    num_cells=[1, 8, 1, 8, 1],
    lengths=[0.5, 1.0, 4.5, 1.0, 0.5],
    oxide_width=5.0, oxide_height=1.5, res=40.0,
    material_stack=SOI_160, mode_solver="tidy3d",
)

wavelengths = np.linspace(1.26, 1.36, 11)
smat = solver.solve(mmi, wavelengths, ports=["out_top"], output_idx=1)
smat.check_passive()
plot_performance(smat, output_dir="out", title="1x2 MMI Splitter (EME)")
```

## Examples

See `examples/` — one folder per device, each with EME/FDTD/hybrid variants (`.py` + matching `.ipynb`):

- `examples/mmi_1x2/` — 1×2 MMI splitter (EME, FDTD, hybrid)
- `examples/directional_coupler/` — directional coupler (FDTD, FDTD length sweep, hybrid, even/odd
  supermode coupling-length analysis)
- `examples/modal_multiplexer/` — modal multiplexer (hybrid)

Run any example directly:

```bash
python examples/mmi_1x2/mmi_1x2_eme.py
```

Output (geometry/field/performance plots, `smatrix.pkl`) goes to `example_output/<device>/<example>/`.

`examples/Makefile` runs examples via `mpirun`/`srun` (auto-detects SLURM/Phoenix), with per-device
and `all` aggregate targets:

```bash
cd examples
make mmi_1x2              # all mmi_1x2 variants
make directional_coupler_fdtd_sweep SWEEP_INDEX=3
make clean
```

## Testing

Fast unit tests (`tests/test_geometry.py`, `test_smatrix.py`, `test_hybrid.py` — no meep/meow
required) run in CI on every push/PR (`.github/workflows/ci.yml`):

```bash
pytest tests/
```

`tests/test_golden_examples.py` additionally runs a real (reduced-resolution) meep/meow
simulation per target — mmi_1x2 (EME, FDTD, hybrid) and directional_coupler (hybrid) — and checks
transmission against a recorded reference in `tests/golden/`. These are genuinely slow (each
target is a real simulation, not a mock, at `res=20` vs the real examples' `res=40`) and are
skipped unless explicitly requested:

```bash
pytest tests/ --runslow
```

Real per-target wall-clock times on a single CPU core (macOS laptop, no MPI parallelism):
mmi_1x2_eme ~7 min, mmi_1x2_fdtd ~12 min, mmi_1x2_hybrid ~58 min, directional_coupler_hybrid
~96 min — even at reduced resolution, a real FDTD/EME solve is real solve time, not a mock; the
hybrid targets are the slowest since HybridSolver's short interior `solve_section()` spans need
`dpml` raised (0.5→1.0 µm) to keep enough PML pixels at `res=20`, which grows the domain. Budget roughly 3 hours to run the full `--runslow` suite; this is why it isn't wired into CI
(see below).

To regenerate the recorded references after a deliberate change to solver behavior:

```bash
python scripts/generate_golden_values.py                # all targets
python scripts/generate_golden_values.py mmi_1x2_eme     # one target
```

If you're on macOS behind a VPN, `meep`'s MPI backend can crash on interpreter exit with an
OFI/libfabric error (`OFI poll failed`) — `tests/conftest.py` sets `FI_PROVIDER=tcp` by default to
avoid this; results are unaffected either way, it's purely a finalize-time crash.

## Status

See `notes.md` for a capability gap analysis against meep/meow/tidy3d/gdsfactory/SAX and a
prioritized list of candidate additions (not yet implemented).
