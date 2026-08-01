"""(Re)generate tests/golden/*.json — the recorded reference values
tests/test_golden_examples.py compares against.

Run from the repo root:
    python scripts/generate_golden_values.py                  # all targets
    python scripts/generate_golden_values.py mmi_1x2_eme       # one target

Each target's build+simulate logic lives once in tests/test_golden_examples.py
(the run_*() functions) and is reused here verbatim, so there's exactly one
place that defines "how to run this target" — this script just runs it and
records the output instead of comparing it.
"""
from __future__ import annotations

import datetime
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tests.test_golden_examples import GOLDEN_DIR, TARGETS  # noqa: E402


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "<not a git checkout>"


def _provenance(wavelengths) -> dict:
    """What this golden was generated against.

    Added 2026-08-01. Until then a golden held `wavelengths` + `transmission`
    and nothing else, so **staleness could not be read off the file** — only
    inferred from the validation docs. These four fixtures have now been
    invalidated four separate times (the Sellmeier switch shifting n_Si ~0.002;
    `N_MODES` 12→8 on the 1×2; the λ₀ re-parameterisation; and the dispersion
    default flipping on), and on none of those occasions did anything in the
    file change to say so.

    That is the same lesson `README.md` check 4 already records for the
    full-device FDTD reference — "it must carry provenance in the file: the
    coupler's old reference was numerically fine but had no config, no solver
    fingerprint and no geometry recorded, so its validity had to be argued from
    circumstantial evidence rather than read off" — applied to the `.json`
    fixtures, where it had never been applied.

    Deliberately records the *material model* rather than a version tag: the
    Sellmeier indices at the run's own band edges and centre move on their own
    whenever the model underneath changes, with nothing to remember to bump.
    Same reasoning as `resolved_material_config()`'s, and the same three probe
    points.

    Not read by `test_golden_examples.py` — that compares `transmission` only,
    so this is descriptive, not part of the assertion. It is here to be read by
    a human asking "is this current?".
    """
    from beef.materials.stack import MaterialStack, sellmeier_index

    wl = [float(w) for w in wavelengths]
    probes = (min(wl), sum(wl) / len(wl), max(wl))
    return {
        "beef_git": _git_sha(),
        "generated_utc": datetime.datetime.now(datetime.timezone.utc)
                                 .replace(microsecond=0).isoformat(),
        "dispersive_default": bool(MaterialStack().dispersive),
        "probe_wl": [round(w, 9) for w in probes],
        "si_index": [round(float(sellmeier_index("si", w)), 9) for w in probes],
        "sio2_index": [round(float(sellmeier_index("sio2", w)), 9) for w in probes],
    }


def generate(name: str) -> None:
    run_fn, ports = TARGETS[name]
    print(f"[{name}] running...")
    smat = run_fn()
    passive = smat.check_passive()
    print(f"[{name}] check_passive() = {passive}")

    data = {
        "wavelengths": smat.wavelengths.tolist(),
        "transmission": {port: smat.transmission(port=port).tolist() for port in ports},
        "provenance": _provenance(smat.wavelengths),
    }
    for port in ports:
        print(f"[{name}]   {port}: T = {['%.4f' % t for t in data['transmission'][port]]}")
    prov = data["provenance"]
    print(f"[{name}]   provenance: beef {prov['beef_git']}, "
          f"dispersive={prov['dispersive_default']}, "
          f"n_Si {prov['si_index'][0]:.6f}..{prov['si_index'][-1]:.6f}")

    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    out_path = GOLDEN_DIR / f"{name}.json"
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"[{name}] wrote {out_path}")


if __name__ == "__main__":
    names = sys.argv[1:] or list(TARGETS)
    unknown = [n for n in names if n not in TARGETS]
    if unknown:
        raise SystemExit(f"unknown target(s) {unknown}, expected one of {list(TARGETS)}")
    for name in names:
        generate(name)
