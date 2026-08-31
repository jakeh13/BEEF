"""The 4x4 90-degree hybrid's design-point constants must not drift apart.

validation_markdown/mmi_4x4_90_hybrid_validation.md, 2026-08-22.

WHY THIS EXISTS. `DELTA_L`, `N_MODES_BODY` and `N_MODES_TAPER` are copied by
hand into a dozen scripts and notebooks for this device. On 2026-08-09 the
body-mode redo (N_MODES_BODY 8->16, which moved the optimum to DELTA_L=-2.75)
was applied to body_sweep.py and hybrid.py but silently missed
mmi_4x4_90hybrid_eme.py and plot_mmi4x4_90hybrid_eme_fields.py. Those two sat
at the superseded point for 12 days, so every EME artifact on disk showed a
design point the project had already rejected. Fixing them on 2026-08-21 then
missed a THIRD constant (N_MODES_TAPER, still 4 where the settled value is 12),
which cost a wasted cluster job (12212950) that exited COMPLETED and produced
plausible-looking but under-resolved plots. Three drifts, three separate
discoveries, all by luck rather than by any check.

WHY A STATIC TEST RATHER THAN A SHARED CONSTANTS MODULE. The obvious fix --
have every script import one params module -- does not work here:
  * `mmi_4x4_90hybrid_eme.py` has no `main()` and no `__main__` guard, and
    runs 8 EME `solve()` calls at module level (>96 h at the current basis).
    Nothing can import it, and making it importable means restructuring a
    shipped example.
  * The three notebooks under mmi_4x4_90hybrid/notebooks/ carry their OWN
    copies of these constants. A .py-only refactor would leave them drifted
    and force a regeneration of all three.
This test reads the files as TEXT, so it needs no imports, has no side
effects, covers the notebooks too, and cannot be defeated by import-time
work.

HISTORICAL SCRIPTS ARE EXEMPT, DELIBERATELY. `mode_convergence.py` and
`taper_sweep.py` ran at DELTA_L=-1.50 and their recorded results correspond
to that geometry; rewriting them would misrepresent what was actually run.
They are listed in HISTORICAL below with the reason. If either is ever
re-used for new work (taper_sweep.py is the natural starting point for a
`taper_width` sweep), move it into CURRENT and update its constants first --
that is the entire point of listing it here rather than just ignoring it.

THE y_gap CAMPAIGN SCRIPTS ARE IN NEITHER LIST, ALSO DELIBERATELY.
`mmi_4x4_90hybrid_ygap_grid.py`, `_ygap_eme.py`, `_ygap_fdtd.py`,
`_wg_eval.py` and `_ygap_collect.py` exist precisely to VARY mm_width,
taper_width and length away from the settled point -- pinning them to it
would defeat their purpose, and exempting them under HISTORICAL would be a
misfiling, since they are current work rather than closed studies. They are
absent by design, not by oversight. Note that `_wg_eval.py` does hardcode
N_MODES=16 to match the settled body-mode count; if SETTLED changes, that
one constant has to be carried over by hand.
"""
from __future__ import annotations

import ast
import json
import os
import re

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(REPO, "examples")

# The settled design point -- Steps 1-3 plus the 2026-08-09 body-mode redo.
SETTLED = {"DELTA_L": -2.75, "N_MODES_BODY": 16, "N_MODES_TAPER": 12}

# Files expected to be AT the settled point. Paths relative to examples/.
CURRENT = [
    "mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_eme.py",
    "mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_fdtd.py",
    "mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_hybrid.py",
    "verification/python/mmi_4x4_90hybrid_offset_sweep.py",
    "verification/python/plot_mmi4x4_90hybrid_eme_fields.py",
    "verification/python/mmi_4x4_90hybrid_eme_timing_probe.py",
    "mmi_4x4_90hybrid/notebooks/mmi_4x4_90hybrid_eme.ipynb",
    "mmi_4x4_90hybrid/notebooks/mmi_4x4_90hybrid_fdtd.ipynb",
    "mmi_4x4_90hybrid/notebooks/mmi_4x4_90hybrid_hybrid.ipynb",
]

# Deliberately NOT held to the settled point -- see the module docstring.
HISTORICAL = {
    "verification/python/mmi_4x4_90hybrid_mode_convergence.py":
        "closed body/taper mode-count study; its recorded results are at DELTA_L=-1.50",
    "verification/python/mmi_4x4_90hybrid_taper_sweep.py":
        "closed taper_length study (measured null); results are at DELTA_L=-1.50/N_MODES_TAPER=4",
}

# Numeric literal inside an assignment's own source segment. Applied ONLY to
# real ast.Assign nodes, never to raw file text -- these scripts discuss their
# own constants in prose ("At the converged N_MODES_BODY=16 design point",
# "N_MODES_BODY=4 proxy extrapolated ..."), and a line-anchored regex over raw
# text matches those comments and reports drift that does not exist. That
# false positive fired on offset_sweep.py the first time this test ran.
NUMBER = re.compile(r"(-?\d+\.?\d*)")


def _text(rel: str) -> str:
    path = os.path.join(EX, rel)
    if not os.path.exists(path):
        pytest.skip(f"{rel} not present")
    raw = open(path, encoding="utf-8", errors="ignore").read()
    if rel.endswith(".ipynb"):
        # Concatenate code-cell source so notebook copies are checked too.
        try:
            nb = json.loads(raw)
        except json.JSONDecodeError:
            return raw
        return "\n".join(
            "".join(c.get("source", []))
            for c in nb.get("cells", [])
            if c.get("cell_type") == "code"
        )
    return raw


def _declared(text: str, name: str):
    """First real module-level assignment of `name`, or None.

    Parses to an AST and inspects only genuine `Assign` nodes, so prose in
    comments and docstrings can never be mistaken for a declaration. Handles
    both `NAME = 12` and `NAME = int(os.environ.get("...", 12))` by pulling
    the first numeric literal out of the assignment's own source segment.
    """
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return None
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            continue
        val = _value_of(node.value)
        if val is not None:
            return val
    return None


def _value_of(node):
    """Numeric value of an assignment's right-hand side, or None.

    Reads real AST literal nodes rather than regexing the source segment.
    Regexing is actively wrong here: the env-var names embed digits
    (`MMI4X4_EME_NMODES_TAPER`), so "first number in the text" returns the
    4 from "MMI4X4" instead of the default 12. That false failure fired on
    eme.py and the field-plot script on this test's second run -- both files
    were correct and the test was not.
    """
    # Plain literals, including negatives (`-2.75` is UnaryOp(USub, Constant)).
    try:
        v = ast.literal_eval(node)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            return float(v)
    except (ValueError, SyntaxError, TypeError):
        pass
    # `int(os.environ.get("NAME", 12))` -- the default is the LAST numeric
    # literal in the call tree; string arguments are ignored entirely.
    nums = [
        n.value for n in ast.walk(node)
        if isinstance(n, ast.Constant)
        and isinstance(n.value, (int, float))
        and not isinstance(n.value, bool)
    ]
    return float(nums[-1]) if nums else None


@pytest.mark.parametrize("rel", CURRENT)
def test_current_scripts_are_at_the_settled_design_point(rel):
    text = _text(rel)
    found = {n: _declared(text, n) for n in SETTLED}
    if all(v is None for v in found.values()):
        pytest.skip(f"{rel} declares none of {list(SETTLED)}")
    bad = {
        n: (v, SETTLED[n])
        for n, v in found.items()
        if v is not None and v != SETTLED[n]
    }
    assert not bad, (
        f"{rel} has drifted from the settled 4x4 design point: "
        + "; ".join(f"{n}={got} but settled value is {want}"
                    for n, (got, want) in bad.items())
        + ". If the design point itself moved, update SETTLED in this test AND "
          "every file in CURRENT together -- that is what this test exists to force."
    )


@pytest.mark.parametrize("rel,reason", sorted(HISTORICAL.items()))
def test_historical_scripts_are_still_historical(rel, reason):
    """Pins the exemptions so they stay deliberate rather than forgotten.

    If one of these is brought back into use, it should move to CURRENT --
    this failing is the reminder to make that call consciously.
    """
    text = _text(rel)
    dl = _declared(text, "DELTA_L")
    assert dl is None or dl != SETTLED["DELTA_L"], (
        f"{rel} now matches the settled design point, but is listed as historical "
        f"({reason}). If it has genuinely been brought up to date and re-run, move "
        f"it from HISTORICAL to CURRENT in this test."
    )
