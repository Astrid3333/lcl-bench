#!/usr/bin/env python3
"""Ejecutar el evaluador de Fortran para goertzel_vs_fft y generar
results/registry_run_f_*.json compatible con hub_compare.py."""
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FORT_EXE = ROOT / "fortran" / "registry_eval_dsp"

# Parámetros del registro goertzel_fft.json
N = 8192
tones = [
    (1, 1.0, 0.3),
    (1000, 0.5, 1.1),
    (4095, 0.25, -0.7),
]

def main():
    RESULTS.mkdir(parents=True, exist_ok=True)

    # Construir argumentos de línea de comandos para el programa Fortran
    args = [str(FORT_EXE), str(N), str(len(tones))]
    for k, A, phi in tones:
        args.extend([str(k), str(A), str(phi)])

    # Ejecutar el programa Fortran
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    lines = result.stdout.strip().splitlines()

    # Parsear la salida: 3 líneas con 4 valores cada una, luego 1 línea con spread
    checks = []
    all_pass = True
    tol = 1.0e-8

    for line_idx, line in enumerate(lines[:3]):
        parts = line.split()
        mag_K = float(parts[0])
        mag_alpha = float(parts[1])
        mag_beta = float(parts[2])
        refv = float(parts[3])

        k = tones[line_idx][0]
        A = tones[line_idx][1]

        # Check K
        rel = abs(mag_K - refv) / refv
        ok = rel <= tol
        all_pass = all_pass and ok
        checks.append({
            "quantity": f"K_k{k}",
            "computed": mag_K,
            "reference": refv,
            "tol": tol * refv,
            "rel_err": rel,
            "pass": ok,
        })

        # Check alpha
        rel = abs(mag_alpha - refv) / refv
        ok = rel <= tol
        all_pass = all_pass and ok
        checks.append({
            "quantity": f"alpha_k{k}",
            "computed": mag_alpha,
            "reference": refv,
            "tol": tol * refv,
            "rel_err": rel,
            "pass": ok,
        })

        # Check beta
        rel = abs(mag_beta - refv) / refv
        ok = rel <= tol
        all_pass = all_pass and ok
        checks.append({
            "quantity": f"beta_k{k}",
            "computed": mag_beta,
            "reference": refv,
            "tol": tol * refv,
            "rel_err": rel,
            "pass": ok,
        })

    # Check variant_spread (línea 4)
    spread = float(lines[3])
    rel = spread / 1.0  # reference = 0
    ok = spread <= tol
    all_pass = all_pass and ok
    checks.append({
        "quantity": "variant_spread",
        "computed": spread,
        "reference": 0.0,
        "tol": tol,
        "rel_err": spread,
        "pass": ok,
    })

    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    out_file = RESULTS / f"registry_run_f_{timestamp}.json"

    run_data = {
        "timestamp": timestamp,
        "selector": "all",
        "n_pass": sum(1 for c in checks if c["pass"]),
        "n_total": len(checks),
        "runs": [
            {
                "name": "goertzel_vs_fft",
                "family": "dsp_bin",
                "pass": all_pass,
                "checks": checks,
            }
        ],
    }

    with open(out_file, "w") as f:
        json.dump(run_data, f, indent=2)

    print(f"Guardado: {out_file}")
    print(f"Entradas: {len(checks)}  Pasaron: {run_data['n_pass']}")

if __name__ == "__main__":
    main()