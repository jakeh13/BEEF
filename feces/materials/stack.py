"""MaterialStack: single source of truth for Si/SiO2 layer definitions.

Replaces the duplicated stack_extrusions_GF() and meep_material_from_meow()
calls scattered across every master and sweep script.
"""
from __future__ import annotations

import numpy as np


class MaterialStack:
    """Si waveguide on SiO2 substrate with SiO2 upper cladding.

    Args:
        wg_height: Si core thickness (µm).
        sz: total vertical simulation span for meow cross-sections (µm).
    """

    def __init__(self, wg_height: float = 0.16, sz: float = 2.0):
        self.wg_height = wg_height
        self.sz = sz

    @property
    def extrusion_rules(self) -> dict:
        """meow GDS extrusion rules (layer → list[GdsExtrusionRule])."""
        import meow as mw

        h = self.wg_height
        sz = self.sz
        return {
            (2, 0): [
                mw.GdsExtrusionRule(material=mw.silicon_oxide, h_min=-sz/2, h_max=-h/2, mesh_order=2),
                mw.GdsExtrusionRule(material=mw.silicon,       h_min=-h/2,  h_max=h/2,  mesh_order=1),
                mw.GdsExtrusionRule(material=mw.silicon_oxide, h_min=h/2,   h_max=sz/2, mesh_order=2),
            ],
            (1, 0): [
                mw.GdsExtrusionRule(material=mw.silicon_oxide, h_min=-sz/2, h_max=sz/2, mesh_order=2),
            ],
        }

    def si_medium(self, wl: float):
        """Return a meep Medium for silicon at wavelength wl (µm)."""
        import meow as mw
        import meep as mp

        ns = np.asarray(mw.silicon.n)
        wls = np.asarray(mw.silicon.wl)
        idx = int(np.argmin(np.abs(wls - wl)))
        return mp.Medium(epsilon=float(ns[idx]) ** 2)

    def sio2_medium(self, wl: float):
        """Return a meep Medium for SiO2 at wavelength wl (µm)."""
        import meow as mw
        import meep as mp

        ns = np.asarray(mw.silicon_oxide.n)
        wls = np.asarray(mw.silicon_oxide.wl)
        idx = int(np.argmin(np.abs(wls - wl)))
        return mp.Medium(epsilon=float(ns[idx]) ** 2)


# Pre-built stacks for common SOI platforms.
SOI_160 = MaterialStack(wg_height=0.16, sz=2.0)   # 160 nm Si (default in this codebase)
SOI_220 = MaterialStack(wg_height=0.22, sz=2.0)   # 220 nm Si (standard commercial SOI)
