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

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tests.test_golden_examples import GOLDEN_DIR, TARGETS  # noqa: E402


def generate(name: str) -> None:
    run_fn, ports = TARGETS[name]
    print(f"[{name}] running...")
    smat = run_fn()
    passive = smat.check_passive()
    print(f"[{name}] check_passive() = {passive}")

    data = {
        "wavelengths": smat.wavelengths.tolist(),
        "transmission": {port: smat.transmission(port=port).tolist() for port in ports},
    }
    for port in ports:
        print(f"[{name}]   {port}: T = {['%.4f' % t for t in data['transmission'][port]]}")

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
