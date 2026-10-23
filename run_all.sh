#!/usr/bin/env bash
# Run the full cross-language suite: Octave + Python + hub diff.
# Usage: ./run_all.sh [selector]
#   selector: all (default) | attractors | potentials | map | ifs | ode |
#             neural | potential | <entry-name>
set -euo pipefail
cd "$(dirname "$0")"
SEL="${1:-all}"

echo "=== 1/3 Octave suite (selector: $SEL) ==="
octave-cli --no-gui --eval "addpath('octave'); registry_run('${SEL}');"

echo ""
echo "=== 2/3 Python suite (selector: $SEL) ==="
python3 python/registry_run.py "$SEL"

echo ""
echo "=== 3/3 Hub: cross-language diff ==="
python3 python/hub_compare.py
