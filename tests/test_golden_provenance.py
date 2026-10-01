"""The golden fixtures must say what they were generated against.

Added 2026-08-01. Until then `tests/golden/*.json` held `wavelengths` and
`transmission` and nothing else, so **staleness could not be read off the
file** — only inferred from the validation record. Those four fixtures have
been invalidated four separate times (the Sellmeier switch, `N_MODES` 12→8 on
the 1×2, the λ₀ re-parameterisation, and the dispersion default flipping on)
and none of those events changed anything in the file.

These tests are deliberately fast and meep-free: they exercise the provenance
builder directly rather than regenerating a golden, which is
`@pytest.mark.slow` and takes hours.
"""
from __future__ import annotations

import pathlib
import sys

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.generate_golden_values import _provenance  # noqa: E402

WL = [1.26, 1.31, 1.36]


def test_provenance_records_the_things_that_have_actually_invalidated_goldens():
    p = _provenance(WL)
    assert p["beef_git"], "no code version recorded"
    assert p["generated_utc"].endswith("+00:00"), "timestamp must be explicit UTC"
    # The dispersion flag: the 2026-08-01 invalidation.
    assert isinstance(p["dispersive_default"], bool)
    # The index model itself: the Sellmeier and λ₀ invalidations. Recorded as
    # resolved numbers rather than a version tag so any change underneath moves
    # them on its own, with nothing to remember to bump.
    assert len(p["si_index"]) == 3
    assert all(3.0 < n < 4.0 for n in p["si_index"])
    assert all(1.3 < n < 1.6 for n in p["sio2_index"])


def test_the_probes_span_the_band_rather_than_only_its_centre():
    """A single midband number cannot see a change in dispersion slope — the
    same reason `resolved_material_config()` samples three points."""
    p = _provenance(WL)
    assert p["probe_wl"] == [1.26, 1.31, 1.36]
    assert p["si_index"][0] != p["si_index"][-1], "band edges must differ"


def test_the_index_probes_track_the_band_they_are_given():
    """Guards against the probes being hardcoded or read off a fixed band."""
    o_band = _provenance([1.26, 1.31, 1.36])["si_index"]
    c_band = _provenance([1.53, 1.5475, 1.565])["si_index"]
    assert o_band != c_band


def test_provenance_is_descriptive_not_asserted():
    """`test_golden_examples.py` compares `transmission` only. If provenance
    ever became part of the assertion, every golden would fail on the next
    commit — the git sha alone changes constantly. Pinned so nobody 'helpfully'
    wires it into the comparison."""
    import inspect

    from tests import test_golden_examples as tge

    src = inspect.getsource(tge._assert_transmission_matches)
    assert "provenance" not in src
    assert "provenance" not in inspect.getsource(tge.test_golden)


@pytest.mark.parametrize("key", ["beef_git", "generated_utc", "dispersive_default",
                                 "probe_wl", "si_index", "sio2_index"])
def test_every_expected_key_is_present(key):
    assert key in _provenance(WL)
