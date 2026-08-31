"""The (mm_width, y_gap) search grid, as one enumeration both the batch
array and the analysis read.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-27.

WHY A SHARED MODULE. The SLURM array maps a task id to a point and the
collector maps a saved .npz back to one; if those two enumerations ever
disagree the campaign silently mislabels results. One list, imported by both.

GRID (user's spec, 2026-08-27): five mm_widths 0.5 um apart about the
nominal 6.0, y_gap on a 0.05 um lattice.

CONSTRAINT. taper_width = mm_width/4 - y_gap must stay strictly above
waveguide_width=0.35, or the taper runs backwards (narrower at the MMI face
than at the access waveguide) and the geometry is unphysical. So

    y_gap < mm_width/4 - 0.35

with the largest lattice point strictly below that bound. The lower bound
Y_GAP_MIN=0.10 is a CHOICE, not a physical limit: below it the taper tips
sit within a wavelength of each other in oxide and the two neighbouring
apertures merge into one aperture, which is a different device topology
rather than a different value of this parameter. Widen it if an optimum
lands on that edge.

The nominal design (6.0, 0.70) lies on the lattice, so the campaign contains
its own reference point.
"""
import numpy as np

WAVEGUIDE_WIDTH = 0.35
Y_GAP_MIN = 0.10
Y_GAP_STEP = 0.05
MM_WIDTHS = [5.0, 5.5, 6.0, 6.5, 7.0]

# Absolute mm_length search grid, shared by wg_eval.py's per-point length
# sweep and by any downstream code (collector/BO driver) that needs to flag
# a best-length hitting the sweep's own edge. 20-400 um at 50 nm steps spans
# multiple self-imaging periods (L_pi ~ 109-150 um across MM_WIDTHS), not
# just a refinement near one guessed order -- see wg_eval.py's own comment
# for why this costs no extra solver time (solve_section_length_sweep's
# tier-3 analytic reconstruction: one mode solve total, any number of
# lengths). One shared module for the same reason MM_WIDTHS/y_gaps_for are
# shared: two independently-hand-kept copies of a range is how enumerations
# silently drift apart.
LENGTH_MIN = 20.0
LENGTH_MAX = 400.0
LENGTH_STEP = 0.05
MM_LENGTH_GRID = np.round(
    np.arange(LENGTH_MIN, LENGTH_MAX + 1e-9, LENGTH_STEP), 2)


def y_gaps_for(mm_width):
    """Lattice points strictly inside the no-inverse-taper bound."""
    hi = mm_width / 4.0 - WAVEGUIDE_WIDTH
    n = int(np.floor((hi - 1e-9 - Y_GAP_MIN) / Y_GAP_STEP))
    return [round(Y_GAP_MIN + k * Y_GAP_STEP, 2) for k in range(n + 1)]


def grid():
    """Width-major, so a whole width's points share one cached body solve."""
    return [(w, y) for w in MM_WIDTHS for y in y_gaps_for(w)]


POINTS = grid()

if __name__ == "__main__":
    for w in MM_WIDTHS:
        ys = y_gaps_for(w)
        print(f"mm_width={w:.1f}  offset={w/4:.3f}  "
              f"y_gap {ys[0]:.2f}-{ys[-1]:.2f}  n={len(ys):2d}  "
              f"taper_width {w/4-ys[-1]:.3f}-{w/4-ys[0]:.3f}")
    print(f"total points = {len(POINTS)}")
