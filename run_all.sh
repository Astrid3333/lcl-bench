#!/usr/bin/env bash
# Run the full cross-language suite: Octave + Python + Fortran + hub diff.
# Usage: ./run_all.sh [selector]
#   selector: all (default) | attractors | potentials | map | ifs | ode |
#             neural | potential | <entry-name>
set -euo pipefail
cd "$(dirname "$0")"
SEL="${1:-all}"

echo "=== 1/4 Octave suite (selector: $SEL) ==="
octave-cli --no-gui --eval "addpath('octave'); registry_run('${SEL}');"

echo ""
echo "=== 2/4 Python suite (selector: $SEL) ==="
python3 python/registry_run.py "$SEL"

echo ""
echo "=== 3/4 Fortran suite (selector: $SEL) ==="
make -s -C fortran f_runner
./fortran/f_runner "$SEL"

echo ""
echo "=== 4/4 Hub: cross-language diff ==="
python3 python/hub_compare.py
