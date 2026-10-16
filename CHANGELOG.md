# Changelog

All notable changes to BEEF are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/).

## [0.1.0] — Unreleased

Initial public release, accompanying the BEEF paper in Computer Physics
Communications.

### Added

- `HybridSolver`: a generic N-section FDTD/EME cascade (any solver type, any
  order, any count ≥ 2) with reciprocity-based excitation, mirror-symmetry
  reuse (`MirrorSpec`, `mirror_pairs`) and FDTD-calibrated body propagation
  constants (`calibrate_body_neff`).
- `FDTDSolver` (meep, 2-D or 3-D) and `EMESolver` (meow), including section
  solves, length sweeps and reflection solves.
- Geometry: `MMI` (N×M, per-port `MMIPort` tapers), `DirectionalCoupler`,
  `UniformCrossSection`, ring-resonator sections and GDS import, all sharing
  the `section_prisms`/`section_gds`/`width_at`/`center_at` interface.
- Materials: `MaterialStack` with dispersive Si/SiO2, and the `SOI_160`/`SOI_220`
  stacks.
- `SMatrix` results with transmission, insertion loss, imbalance, CMRR and
  phase metrics, plus Redheffer-star S-matrix composition in `beef.network`.
- `beef.mode_identity`: field-overlap mode correspondence across FDTD/EME joins.
- A content-addressed simulation cache, shipped as `example_sim_cache/` so the
  examples reproduce without re-running FDTD.
- Examples with scripts and matching notebooks: 1×2, 1×4, 2×2 and tapered 2×2
  MMIs, a 4×4 MMI 90° optical hybrid, a directional coupler, a high-contrast
  grating mirror and Fabry–Perot cavity, and a ring resonator.
