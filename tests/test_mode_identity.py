"""Do FDTD and EME agree on *which* mode is which? (dc_validation.md item 12)

`HybridSolver` contracts index i of an FDTD interface against index i of an EME
body, assuming they are the same mode. The two solve MPB on different grids --
meep's Yee mesh vs meow's Mesh2D -- so their n_eff differ slightly and
near-degenerate modes could in principle come back reordered, which would make
the chain mix modes with nothing else noticing.

MEASURED, at 1.55 um on the coupler's coupled-pair cross-section (res 40 both
sides):

    band    meep      meow      diff
      0   2.324450  2.325116  -0.0007
      1   2.300984  2.301519  -0.0005
      2   1.744418  1.746722  -0.0023
      3   1.631969  1.634080  -0.0021
      4   1.466926  1.467593  -0.0007
      5   1.436725  1.430456  +0.0063

Ordering is preserved and the worst pairwise disagreement is 0.0063. So the
contraction is sound -- for this device.

TWO THINGS THIS TURNED UP
-------------------------
1. **Band 5 is not a guided mode.** The cladding index is 1.444 and band 5 sits
   at 1.4367, below it; band 4 at 1.4669 is marginal. The directional coupler
   carries N_MODES = 6 at every interface, so at least one column of that basis
   is a radiation mode, whose n_eff and ordering are not well defined and whose
   agreement between engines (the +0.0063 outlier, an order of magnitude worse
   than any guided band) is correspondingly the worst in the table. Guided modes
   are the ones the assertions below are allowed to be strict about.

2. **Comparing across resolutions is meaningless, and it is an easy mistake.**
   The res 40 -> 80 shift in n_eff is ~0.027, larger than the 0.022 even/odd
   splitting itself, so a nearest-neighbour match between meep@40 and meow@80
   pairs band 0 with band 1 and looks exactly like a reordering. It is not. Only
   same-resolution comparison is valid, and `assert_modes_correspond()` below
   takes both arrays from one resolution for that reason.

STATUS, 2026-07-30 -- **`assert_modes_correspond()` below is re-scoped, not
retired.** It
predates `beef.mode_identity`'s transverse-E field-overlap machinery
(`field_overlap_matrix()`/`best_correspondence()`/`overlap_diagnostics()`,
also reachable at runtime via `HybridSolver(check_mode_correspondence=True)`),
which is the trusted, device- and resolution-agnostic mechanism now -- see
that module's own docstring for why a 1-D n_eff heuristic is fundamentally
weaker (point 2 above is one symptom; it also cannot tell "these are
different physical modes that happen to have close n_eff" from a genuine
reordering, which only comparing *shape* can). Kept here only as a fast
smoke check (no field solve, no shared diagnostic mesh) for the two
mismatches it CAN still catch cheaply -- a gross swap or a missing mode --
not as proof of correspondence. Do not add new call sites; use
`beef.mode_identity`'s field-overlap machinery for anything that needs an
authoritative answer.
"""
from __future__ import annotations

import numpy as np
import pytest

WL = 1.55
RES = 40.0
N_MODES = 6

# meow's SiO2 at 1.55 um. A mode below this is not bound -- it is a radiation
# mode living on the discretization, and neither engine owes the other an answer
# about it.
N_CLADDING = 1.444

# Worst pairwise |dn_eff| seen across the guided bands was 0.0023. 0.01 leaves
# room for a different MPB build without admitting a genuine mode swap, which
# would show up as a difference of order the mode spacing (0.02 at the tightest,
# between the even and odd supermodes).
NEFF_TOL = 0.01

DEVICE = dict(
    waveguide_width=0.45, gap=0.25, coupling_length=30.0,
    bend_length=9.0, waveguide_offset=3.0, bend_type="euler",
    waveguide_height=0.22,
)
SECTION = 2.0   # um, centred at x=0: deep inside the 30 um coupling region


def assert_modes_correspond(neff_a, neff_b, tol=NEFF_TOL, label=("a", "b")):
    """Both arrays must be same-resolution solves of the same cross-section.

    Checks the three things that have to hold for HybridSolver's index-i-to-
    index-i contraction to be valid: both engines order modes the same way,
    index i names the same mode in both, and their n_eff agree.

    **Fast smoke check, not the authoritative one -- see this module's
    docstring, "STATUS, 2026-07-30".** This is a 1-D scalar (n_eff) heuristic:
    it requires matched resolution (a resolution mismatch fakes a reorder,
    see this module's docstring point 2) and cannot distinguish two
    genuinely different mode shapes that happen to land at close n_eff from
    an actual index swap -- only comparing field *shape*
    (`beef.mode_identity.field_overlap_matrix()`/`best_correspondence()`, or
    `HybridSolver(check_mode_correspondence=True)` at runtime) can. Use this
    only for a cheap same-resolution sanity check; it does not substitute
    for the field-overlap mechanism anywhere correctness actually matters.
    """
    a, b = np.asarray(neff_a, float), np.asarray(neff_b, float)
    assert a.shape == b.shape, f"mode counts differ: {a.shape} vs {b.shape}"
    assert np.all(np.diff(a) < 0), f"{label[0]} not in descending n_eff: {a}"
    assert np.all(np.diff(b) < 0), f"{label[1]} not in descending n_eff: {b}"

    # If index i named a different mode in the two engines, i's nearest
    # neighbour in the other list would not be i.
    nearest = [int(np.argmin(np.abs(b - v))) for v in a]
    assert nearest == list(range(len(a))), (
        f"mode correspondence broken: {label[0]}[i] is nearest to "
        f"{label[1]}[{nearest}], not [0..{len(a)-1}]. HybridSolver would "
        f"contract mismatched modes.\n  {label[0]} = {a}\n  {label[1]} = {b}"
    )
    worst = float(np.max(np.abs(a - b)))
    assert worst < tol, f"worst |dn_eff| {worst:.4f} exceeds {tol}"
    return worst


class TestModeCorrespondenceLogic:
    """Fast, no simulation. Guards the checker itself -- if it cannot detect a
    swap it would pass forever and mean nothing."""

    GOOD_A = np.array([2.324450, 2.300984, 1.744418, 1.631969])
    GOOD_B = np.array([2.325116, 2.301519, 1.746722, 1.634080])

    def test_accepts_the_measured_agreement(self):
        worst = assert_modes_correspond(self.GOOD_A, self.GOOD_B)
        assert worst == pytest.approx(0.0023, abs=1e-4)

    def test_detects_a_swap_of_adjacent_modes(self):
        """The failure this whole item exists to catch: the two engines return
        the same modes but engine B has the top two the other way round."""
        swapped = self.GOOD_B.copy()
        swapped[0], swapped[1] = swapped[1], swapped[0]
        with pytest.raises(AssertionError, match="descending|correspondence"):
            assert_modes_correspond(self.GOOD_A, swapped)

    def test_detects_a_missing_mode(self):
        with pytest.raises(AssertionError, match="mode counts differ"):
            assert_modes_correspond(self.GOOD_A, self.GOOD_B[:-1])

    def test_rejects_a_large_disagreement_even_when_correspondence_holds(self):
        """Perturbs only band 3, whose nearest neighbour is 0.112 away, so the
        correspondence check still passes and the tolerance is what fires.

        Perturbing every mode by the same amount would NOT test this: a uniform
        offset larger than the tightest mode spacing (0.023, bands 0-1) breaks
        correspondence first, which is the cross-resolution effect below rather
        than a magnitude failure. Getting that wrong is what this test is for.
        """
        assert_modes_correspond(self.GOOD_A, self.GOOD_A + 0.001)
        perturbed = self.GOOD_A.copy()
        perturbed[3] += 0.05
        with pytest.raises(AssertionError, match="exceeds"):
            assert_modes_correspond(self.GOOD_A, perturbed)

    def test_cross_resolution_comparison_is_why_this_takes_one_resolution(self):
        """Documents the trap. meep@res40 against meow@res80 -- both correct,
        neither reordered -- fails correspondence, because the res 40->80 shift
        (~0.027) is larger than the even/odd splitting (~0.022). The measurement
        must be same-resolution; this test pins the reason so nobody 'fixes' the
        checker to accept it."""
        meep40 = np.array([2.324450, 2.300984])
        meow80 = np.array([2.350893, 2.328838])
        assert abs(meep40[0] - meow80[1]) < abs(meep40[0] - meow80[0]), (
            "premise: meep's band 0 really is numerically nearer meow@80's band 1"
        )
        with pytest.raises(AssertionError, match="correspondence|exceeds"):
            assert_modes_correspond(meep40, meow80)


def _neff_both_engines(n_modes=N_MODES):
    """(meep, meow) n_eff for the same coupled-pair cross-section at RES."""
    import meow as mw
    from beef import DirectionalCoupler, EMESolver, FDTDSolver, SOI_220, extrude_polygons

    dc = DirectionalCoupler(**DEVICE)
    fdtd = FDTDSolver(
        res=RES, sim_width=6.0, sim_height=1.5, dpml=0.5,
        source_margin=0.0, monitor_margin=0.25, three_d=True,
        half_domain_aperture=True, material_stack=SOI_220, decay_tol=1e-11,
    )
    meep_neff = fdtd.section_mode_neff(
        dc, np.array([WL]), -SECTION / 2, +SECTION / 2, n_modes,
        monitor_center=0.0, monitor_size=6.0,
    )[0]

    eme = EMESolver(
        num_modes=[n_modes], num_cells=[1], lengths=[1.0],
        oxide_width=6.0, oxide_height=1.5, res=RES,
        material_stack=SOI_220, mode_solver="meep",
    )
    gds = dc.section_gds(-SECTION / 2, +SECTION / 2, eme.oxide_width)
    cells = mw.create_cells(
        structures=extrude_polygons(gds, eme.material_stack.extrusion_rules),
        mesh=eme._build_mesh(), Ls=[SECTION],
    )
    cs = mw.CrossSection.from_cell(
        cell=cells[0], env=mw.Environment(wl=WL, T=eme.temperature))
    meow_neff = np.array([float(np.real(m.neff))
                          for m in eme.compute_modes(cs, n_modes)])
    return meep_neff, meow_neff


@pytest.mark.slow
def test_engines_agree_on_mode_identity_for_guided_modes():
    """A fast n_eff smoke check on the coupler's own coupled-pair
    cross-section -- NOT the authoritative correspondence check (see this
    module's "STATUS, 2026-07-30" docstring note); that is
    `beef.mode_identity`'s field-overlap machinery /
    `HybridSolver(check_mode_correspondence=True)`. Scoped to guided modes: a
    band below the cladding index is a radiation mode, and requiring two
    different grids to agree on one is not a meaningful demand."""
    meep_neff, meow_neff = _neff_both_engines()
    guided = meep_neff > N_CLADDING
    assert guided.sum() >= 2, (
        f"expected at least the even/odd supermode pair to be guided, got "
        f"{guided.sum()} above n_clad={N_CLADDING}: {meep_neff}"
    )
    assert_modes_correspond(meep_neff[guided], meow_neff[guided],
                            label=("meep", "meow"))


@pytest.mark.slow
def test_the_interface_basis_contains_an_unguided_mode():
    """N_MODES = 6 is what the directional coupler examples carry, but only 5
    bands are above the cladding index -- so one column of every interface basis
    is a radiation mode. Asserted rather than left as prose because it bears on
    the right N_MODES for every hybrid example, and most sharply on the 15- and
    16-mode MMI and multiplexer, which have never been checked at all.

    If this starts failing, the mode count and the cross-section stopped
    matching in a way worth knowing about.
    """
    meep_neff, _ = _neff_both_engines()
    n_guided = int((meep_neff > N_CLADDING).sum())
    assert n_guided < N_MODES, (
        f"all {N_MODES} interface modes are now guided (n_eff > {N_CLADDING}): "
        f"{meep_neff}. The item 12 finding no longer holds."
    )
