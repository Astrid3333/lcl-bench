#!/usr/bin/env bash
# Run the full cross-language suite: Octave + Python + Fortran + hub diff
# + optional HDF5 bundle. Usage: ./run_all.sh [selector]
#   selector: all (default) | attractors | potentials | map | ifs | ode |
#             neural | potential | <entry-name>
set -euo pipefail
cd "$(dirname "$0")"
SEL="${1:-all}"

echo "=== 1/5 Octave suite (selector: $SEL) ==="
octave-cli --no-gui --eval "addpath('octave'); registry_run('${SEL}');"

echo ""
echo "=== 2/5 Python suite (selector: $SEL) ==="
python3 python/registry_run.py "$SEL"

echo ""
echo "=== 3/5 Fortran suite (selector: $SEL) ==="
make -s -C fortran f_runner
./fortran/f_runner "$SEL"

echo ""
echo "=== 4/5 Hub: cross-language diff ==="
python3 python/hub_compare.py

echo ""
echo "=== 5/5 HDF5 bundle (optional) ==="
if [ -x .venv/bin/python ] && .venv/bin/python -c "import h5py" 2>/dev/null; then
  .venv/bin/python python/export_hdf5.py
elif python3 -c "import h5py" 2>/dev/null; then
  python3 python/export_hdf5.py
else
  echo "SKIP: h5py not available (create .venv with h5py, or pip install h5py)"
fi
