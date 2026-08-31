#!/bin/bash
# Re-runs the (mm_width, y_gap) campaign collector every time a new eval_*.npz
# lands, so the summary/plot in example_output/.../wg_ygap_opt/ stay live for
# the whole 106-point Phase B run without anyone re-invoking it by hand.
# Usage: bash mmi_4x4_90hybrid_ygap_watch.sh   (run from examples/verification/python/)
set -e
cd "$(dirname "$0")"
OUT_DIR="../../../../example_output/mmi_4x4_90hybrid/wg_ygap_opt"
PY=python
LAST_COUNT=-1
while true; do
    N=$(find "$OUT_DIR" -maxdepth 1 -name "eval_W*.npz" 2>/dev/null | wc -l)
    if [ "$N" != "$LAST_COUNT" ]; then
        "$PY" mmi_4x4_90hybrid_ygap_collect.py > /tmp/ygap_collect_last.log 2>&1 || true
        echo "$(date -Iseconds) refreshed at $N points"
        LAST_COUNT=$N
    fi
    sleep 120
done
