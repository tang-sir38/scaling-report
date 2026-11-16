#!/usr/bin/env bash
# Run the shipped examples end to end. No arguments, no network.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PYTHON:-python}

echo "== tests"
$PY tests/test_core.py

echo
echo "== synthetic data with a known answer"
$PY examples/make_synthetic.py --out examples/synthetic_moments.npz
$PY -m scaling_report examples/synthetic_moments.npz --name synthetic --out out/synthetic

echo
echo "== published colon section"
$PY -m scaling_report examples/colon_moments.npz \
    --name "Visium HD colon, section P0" --out out/colon

echo
echo "done. reports in out/"
