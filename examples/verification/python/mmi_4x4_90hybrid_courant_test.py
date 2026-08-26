"""Is the decay floor a NUMERICAL instability, or something physical?

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-24.

THE OBSERVATION being explained. At res=20 this device's full-device FDTD run
decays normally, bottoms out at energy ratio **1.682e-4 at t=3705**, and then
grows exponentially (5.16e-3 -> 51.9 -> 2.14e8 by t=8030). Three independent
jobs hit that floor at bit-identical value and time, so it is deterministic.
The floor is resolution-dependent (1.68e-4 at res=20, ~5e-6 at res=40), which
already points numerical rather than physical.

THE HYPOTHESIS. `beef/materials/stack.py` builds every material from an
UNDAMPED Lorentzian pole (`gamma=0`), and for this device's O-band the
re-parameterised pole sits at lam0 = 1.5640 um -- only 15% above the band edge
at 1.36 um (`_REPARAM_LAM0_FACTOR = 1.15`; the code's own
`_REPARAM_MIN_LAM0_RATIO = 1.10` guard exists because "the fitted term blows up
as the band approaches lam0"). Undamped, near-resonant Lorentzian media are a
known marginal-stability case in FDTD.

WHY THIS TEST AND NOT THE gamma FIX FIRST. Adding damping means editing shared
material code that every device draws on, on a diagnosis that is so far only
INFERRED from the signature. Courant number is the cheap discriminator: it
changes only the timestep, touches no material model, no geometry, and no
physics. If the floor is a numerical instability, halving the Courant number
should push it DOWN and LATER. If the floor is unchanged, the gamma=0
hypothesis loses most of its support and the effect is something else.

    floor moves down/later  -> numerical; the gamma fix is worth doing
    floor unchanged          -> NOT a timestep-stability effect; do not touch
                                the material model on this evidence

Runs at res=20 deliberately: that configuration reaches its floor at t=3705 in
about 2 h at 96 ranks, and its floor value is already measured to compare
against. A lower Courant means proportionally more timesteps for the same
simulated time, so budget ~2x.

Courant is injected by wrapping `mp.Simulation` for the duration of the call --
NOT by subclassing FDTDSolver, which would break cache fingerprinting
(cache/config.py keys on type(solver).__name__) and silently drop cache_only
forwarding (hybrid.py introspects the method signature).

Run:
    MMI4X4_COURANT=0.25 python mmi_4x4_90hybrid_courant_test.py
"""
import os

import numpy as np

from beef import MMI, FDTDSolver, SOI_220

COURANT = float(os.environ.get("MMI4X4_COURANT", 0.25))
RES = float(os.environ.get("MMI4X4_COURANT_RES", 20.0))

WAVEGUIDE_WIDTH = 0.35
MM_WIDTH = 6.0
WAVEGUIDE_OFFSET = MM_WIDTH / 4.0
TAPER_WIDTH = 0.8
TAPER_LENGTH = 1.0
L_PI = 108.9195
DELTA_L = -2.75
MM_LENGTH = 0.75 * L_PI + DELTA_L   # 78.9396

mmi = MMI(
    n_inputs=4, n_outputs=4,
    mm_length=MM_LENGTH, mm_width=MM_WIDTH,
    waveguide_width=WAVEGUIDE_WIDTH, taper_length=TAPER_LENGTH,
    taper_width=TAPER_WIDTH, waveguide_offset=WAVEGUIDE_OFFSET,
    waveguide_height=0.22, use_neighbor_aperture=True,
)

solver = FDTDSolver(
    res=RES, sim_width=MM_WIDTH + 3.0, sim_height=2.0, dpml=0.5,
    source_margin=0.0, monitor_margin=0.25, three_d=True,
    num_input_modes=1, num_output_modes=1,
    output_monitor_centers=list(mmi.output_ys),
    output_monitor_widths=[WAVEGUIDE_OFFSET] * 4,
    material_stack=SOI_220, decay_tol=1e-11,
    decay_method="energy", energy_decay_dt=20.0,
)

print(f"COURANT_TEST: Courant={COURANT} (meep default 0.5) res={RES}", flush=True)
print("COURANT_TEST: reference at Courant=0.5, res=20 -- floor ratio "
      "1.682e-04 at t=3705, then exponential growth", flush=True)
print("COURANT_TEST: floor LOWER/LATER => numerical instability (gamma fix "
      "justified); floor UNCHANGED => not a timestep effect", flush=True)

import meep as mp

_orig_init = mp.Simulation.__init__


def _patched_init(self, *args, **kwargs):
    # Force the Courant number regardless of what the solver passed (it
    # currently passes none, so meep would default to 0.5).
    kwargs["Courant"] = COURANT
    return _orig_init(self, *args, **kwargs)


mp.Simulation.__init__ = _patched_init
try:
    wavelengths = np.linspace(1.26, 1.36, 3)
    # No save_dir: Courant is NOT part of the cache fingerprint, so a result
    # produced at a non-default Courant must never be written where a
    # default-Courant request could find it.
    smat = solver.solve(mmi, wavelengths, port="bot", input_mode=1)
    print("COURANT_TEST: run reached decay_tol (no floor hit)", flush=True)
finally:
    mp.Simulation.__init__ = _orig_init
