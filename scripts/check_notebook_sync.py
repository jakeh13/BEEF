"""Audit each example's .ipynb against its .py.

Notebooks here are hand-mirrored copies of the scripts, differing only in
`CACHE_ONLY = True`, ASCII-ised units, and `plt.show()` in place of `savefig`.
Nothing enforces that they stay in sync, and the predicted drift has already
happened once by audit (2026-07-31: four mismatched constants, one of them
reproducing a known-bad mode basis in a published tutorial).

This is the cheap half of that item's recommendation — the CI check. It does
not generate notebooks from scripts; it reports where they have diverged, on
two axes:

  constants  — a named set of physics/config knobs must match exactly
  calls      — API calls present in the script but missing from the notebook,
               which is how a newly-added block (e.g. plot_beta_calibration)
               silently fails to reach the tutorial

Exit code 1 if anything is out of sync, so it can gate CI.

    python scripts/check_notebook_sync.py [--verbose]
"""
import glob
import json
import os
import re
import sys

CONSTANTS = [
    "N_MODES", "RESOLUTION", "WL_MIN", "WL_MAX", "N_WL", "CALIBRATE_BODY",
    "COUPLING_LENGTH", "SYMMETRIC_MODE_REUSE", "BODY_STUB", "INTERFACE_MODES",
    "L_ARM", "MM_LENGTH", "SEED_ANCHORS",
    # Added 2026-08-18. CALIBRATE_ARM/ARM_RESOLUTION are hcg_fabry_perot_
    # hybrid.py's arm-calibration switch and its mesh -- the exact analogue of
    # CALIBRATE_BODY above, and the pair that decides whether the connecting
    # arms propagate on meep's beta or MPB's index (worth 4.28 FSR of comb
    # position, so a notebook silently on the wrong one teaches the wrong
    # device). The other four were already load-bearing in the coupler/MMI
    # scripts and simply ungated: this checker exists because a stale N_WL went
    # unnoticed for eight days, and a stale resolution is the same failure.
    # Verified to surface no pre-existing drift when added (18 pairs, 0 out of
    # sync), so this is coverage, not a new backlog.
    "CALIBRATE_ARM", "ARM_RESOLUTION",
    "FDTD_RESOLUTION", "EME_RESOLUTION", "N_WL_HCG", "N_WL_DENSE",
    # Added 2026-08-22, and this one DID surface pre-existing drift rather
    # than being pure coverage. These three are the 4x4 90-degree hybrid's
    # own design point. The 2026-08-09 body-mode redo (N_MODES_BODY 8->16,
    # moving the optimum to DELTA_L=-2.75) reached body_sweep.py and
    # hybrid.py but missed mmi_4x4_90hybrid_eme.py, its field-plot companion,
    # AND mmi_4x4_90hybrid_eme.ipynb. The scripts were fixed 2026-08-21/22;
    # the NOTEBOOK was still at DELTA_L=-1.50/N_MODES_BODY=8 and nothing
    # flagged it -- precisely the "known-bad mode basis in a published
    # tutorial" failure this checker's own docstring says it exists to stop,
    # recurring because the constant list did not cover it. MM_LENGTH was
    # already listed but is DERIVED here (0.75*L_PI + DELTA_L), so it only
    # catches this drift where a script hardcodes it.
    "DELTA_L", "N_MODES_BODY", "N_MODES_TAPER",
]
# calls whose absence from a notebook means the tutorial is teaching an older
# version of the workflow, not merely a cosmetic difference
CALLS = [
    "plot_beta_calibration", "section_mode_beta_details", "section_mode_beta",
    "solve_section_length_sweep", "solve_section_multi_length",
    "calibrate_body_neff", "plot_hybrid_geometry", "check_passive",
]

def nb_source(path):
    d = json.load(open(path))
    return "\n".join("".join(c.get("source", [])) for c in d["cells"])

def strip_prose(text):
    """Drop the module docstring and every comment.

    Both matter for false positives, and both were found by running this on the
    real tree: `directional_coupler_soi160_hybrid.py` explains its
    `CALIBRATE_BODY=True` choice *inside its docstring*, at the start of a
    line, which the constant regex happily matched as an assignment; and
    `mmi_1x2_hybrid.py` mentions `section_mode_beta()` in a comment while
    reaching the calibration through `HybridSolver(calibrate_body_neff=True)`
    instead. A checker that cries wolf on prose gets switched off.
    """
    text = re.sub(r'^\s*(?:"""|\'\'\')(?:.|\n)*?(?:"""|\'\'\')', "", text, count=1)
    return re.sub(r"#[^\n]*", "", text)

def _normalise(value):
    """`int(os.environ.get("X", "11"))` and a bare `11` are the same setting.

    Scripts take env-var overrides so a sweep can drive them; notebooks are
    tutorials and hard-code the default. That is a deliberate difference in
    form, not drift, so compare the DEFAULT a script would use against the
    notebook's literal.
    """
    m = re.search(r"""os\.environ\.get\(\s*["'][^"']+["']\s*,\s*["']([^"']*)["']\s*\)""", value)
    return m.group(1) if m else value

def constants(text):
    out = {}
    for k in CONSTANTS:
        # last assignment wins — a later one overrides whatever came before it
        ms = re.findall(rf"^{k}\s*=\s*([^\n#]+)", text, re.M)
        if ms:
            out[k] = _normalise(ms[-1].strip())
    return out

def main(verbose=False):
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "examples")
    problems = 0
    checked = 0
    for nb in sorted(glob.glob(os.path.join(root, "*", "notebooks", "*.ipynb"))):
        name = os.path.basename(nb)[:-6]
        py = os.path.join(os.path.dirname(nb), "..", "python", name + ".py")
        if not os.path.isfile(py):
            print(f"?? {name}: no paired script at python/{name}.py")
            continue
        checked += 1
        ptxt, ntxt = strip_prose(open(py).read()), strip_prose(nb_source(nb))
        pc, nc = constants(ptxt), constants(ntxt)

        bad = []
        for k in sorted(set(pc) | set(nc)):
            if pc.get(k) != nc.get(k):
                bad.append(f"    {k}: script={pc.get(k)!r} notebook={nc.get(k)!r}")
        for call in CALLS:
            # an actual call/kwarg use, not a mention: `name(` or `name=`
            pat = rf"\b{call}\s*[\(=]"
            if re.search(pat, ptxt) and not re.search(pat, ntxt):
                bad.append(f"    missing call: {call}()")
        if bad:
            problems += 1
            print(f"DRIFT {name}")
            print("\n".join(bad))
        elif verbose:
            print(f"ok    {name}")
    print(f"\n{checked} pairs checked, {problems} out of sync")
    return 1 if problems else 0

if __name__ == "__main__":
    sys.exit(main("--verbose" in sys.argv))
