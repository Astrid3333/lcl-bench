#!/bin/bash
# Ejecutar el evaluador de Fortran para goertzel_vs_fft (dsp_bin)
# Genera results/registry_run_f_*.json para el hub de comparación

set -e

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS="$ROOT/results"
FORT="$ROOT/fortran"
REG="$ROOT/registries/goertzel_fft.json"

mkdir -p "$RESULTS"

# Compilar el evaluador de Fortran
gfortran -O2 -o "$FORT/registry_eval_dsp" "$FORT/registry_eval_dsp.f90"

# Parámetros del registro goertzel_fft.json (hardcodeados para evitar depender de jq)
N=8192
NUM_TONES=3
K1=1; A1=1.0; PHI1=0.3
K2=1000; A2=0.5; PHI2=1.1
K3=4095; A3=0.25; PHI3=-0.7

# Ejecutar el evaluador de Fortran y capturar la salida JSON de las checks
CHECKS_JSON="$("$FORT/registry_eval_dsp" "$N" "$NUM_TONES" \
    "$K1" "$A1" "$PHI1" \
    "$K2" "$A2" "$PHI2" \
    "$K3" "$A3" "$PHI3")"

# Generar el archivo de resultado en formato compatible con hub_compare.py
TIMESTAMP=$(date +"%Y%m%dT%H%M%S")
OUT_FILE="$RESULTS/registry_run_f_${TIMESTAMP}.json"

# Contar checks PASS
N_PASS=$(echo "$CHECKS_JSON" | grep -c '"pass": true' || echo 0)
N_TOTAL=$(echo "$CHECKS_JSON" | grep -c '"quantity"' || echo 0)

# Determinar pass general
if [ "$N_PASS" -eq "$N_TOTAL" ] && [ "$N_TOTAL" -gt 0 ]; then
    OVERALL_PASS=true
else
    OVERALL_PASS=false
fi

cat > "$OUT_FILE" << EOF
{
  "timestamp": "$TIMESTAMP",
  "selector": "all",
  "n_pass": $N_PASS,
  "n_total": $N_TOTAL,
  "runs": [
    {
      "name": "goertzel_vs_fft",
      "family": "dsp_bin",
      "pass": $OVERALL_PASS,
      "checks": [
$CHECKS_JSON
    ]
  ]
}
EOF

echo "Guardado: $OUT_FILE"
echo "Entradas: $N_TOTAL  Pasaron: $N_PASS"