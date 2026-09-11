"""Re-select mm_length under a DIFFERENT rule, from cache, for free.

Every point in this campaign chose its mm_length by minimising wg_eval's live
objective, which grades phase linearly and therefore keeps buying phase margin
below the 5 deg spec at the expense of CMRR. The band probe (job 12831929)
showed the device is bandwidth-limited, not loss-limited -- at 1.31 um the best
geometry already makes 20.761 dBe CMRR and 0.654 dB excess loss, and the whole
shortfall sits at the band edges. A length chosen to BALANCE the band rather
than to minimise phase may therefore lift the band-worst CMRR over the bar.

The 7601-length sweep is already cached for every evaluated geometry, so this
costs no simulation: it re-scores the same lengths under three rules and prints
the best under each.

  live       : wg_eval's objective (what every published number used)
  constraint : minimise 0.5*max(0,(dphi-5)/5) + 0.5*max(0,(20-CMRR)/20)
  cmrr_max   : maximise band-worst CMRR subject to max|dphi| <= 5 deg

Usage:
    MMI4X4_W=8.5 MMI4X4_YGAP=0.825 MMI4X4_TAPER_LENGTH=4.0 \
        python mmi_4x4_90hybrid_length_rescore.py
"""
import os

os.environ["MMI4X4_CACHE_ONLY"] = "1"

import numpy as np

import mmi_4x4_90hybrid_wg_eval as W


def main():
    print(f"RESCORE: W={W.MM_WIDTH} y_gap={W.Y_GAP} "
          f"taper_width={W.TAPER_WIDTH:.3f} taper_length={W.TAPER_LENGTH} "
          f"N_MODES={W.N_MODES}  {len(W.MM_LENGTH_GRID)} lengths, "
          f"{len(W.wavelengths)} wl", flush=True)

    input_solver, body_solver = W.make_solvers()
    mmi = W.make_mmi(W.L_CENTRE)
    boundaries, body_stub = W.boundaries_for(mmi)
    M_in = [input_solver.solve_section(
        mmi, W.wavelengths, boundaries[0], boundaries[1],
        n_modes_in=1, n_modes_out=W.N_MODES,
        stub_left=input_solver.waveguide_length, stub_right=0.0,
        input_idx=k, save_dir=W.save_dir, cache_only=True) for k in range(4)]
    body_lengths = W.MM_LENGTH_GRID - 2 * body_stub
    M_body = body_solver.solve_section_length_sweep(
        mmi, W.wavelengths, boundaries[1], boundaries[2], body_lengths,
        n_modes_in=W.N_MODES, n_modes_out=W.N_MODES, save_dir=W.save_dir,
        cache_only=True)

    rows = []
    for li in range(len(W.MM_LENGTH_GRID)):
        md, mc, xl, ib = W.metrics_at(M_in, M_body[li:li + 1])
        rows.append((float(W.MM_LENGTH_GRID[li]), md, mc, xl, ib))
    a = np.array(rows)
    L, MD, MC, XL, IB = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4]

    live = 0.5 * (MD / 5.0) + 0.5 * np.maximum(0.0, (20.0 - MC) / 20.0)
    cons = (0.5 * np.maximum(0.0, (MD - 5.0) / 5.0)
            + 0.5 * np.maximum(0.0, (20.0 - MC) / 20.0))
    ok = MD <= 5.0

    def show(tag, i):
        print(f"{tag:>11}  mm_length={L[i]:8.2f}  max|dphi|={MD[i]:6.3f}  "
              f"min CMRR={MC[i]:7.3f}  XL={XL[i]:6.3f}  imb={IB[i]:6.3f}  "
              f"live={live[i]:.4f}  cons={cons[i]:.4f}"
              f"{'' if ok[i] else '   <- PHASE FAILS'}", flush=True)

    print(f"\nRESCORE: {ok.sum()} of {len(L)} lengths meet the 5 deg phase spec")
    show("live", int(np.argmin(live)))
    show("constraint", int(np.argmin(cons)))
    if ok.any():
        idx = np.where(ok)[0]
        show("cmrr_max", int(idx[np.argmax(MC[idx])]))
        print(f"\nRESCORE: best band-worst CMRR among phase-compliant lengths = "
              f"{MC[idx].max():.3f} dBe "
              f"({'MEETS' if MC[idx].max() >= 20 else 'misses'} the 20 dBe spec)")
    else:
        print("\nRESCORE: no length meets the phase spec at this geometry")


if __name__ == "__main__":
    main()
