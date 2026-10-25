#!/usr/bin/env python3
"""Assert the preregistered verdict table for every diagnosis artefact.

25 assertions over 13 artefacts: 7 stored diagnoses (results/diagnose) +
6 legacy-source reports (results/legacy). Exit 0 only if every verdict AND
every execution status matches the preregistration; any mismatch is printed
and fails the run.

Usage: python3 python/check_verdicts.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# preregistered expectations: file -> (diagnostico, execution_status_or_None)
DIAGNOSE = {
    "lorenz_reference.json": ("equivalente", None),
    "lorenz_euler.json": ("degradacion_de_orden", None),
    "lorenz_euler_large.json": ("inestabilidad_numerica", None),
    "logistic_reference.json": ("equivalente", None),
    "logistic_float32.json": ("colapso_por_redondeo", None),
    "henon_series_ok.json": ("equivalente", None),
    "henon_series_collapse.json": ("colapso_por_redondeo", None),
}
LEGACY = {
    # reference-era double RK4 must match the registry fingerprint
    "lorenz_rk4_double.f.json": ("equivalente", "executed_source", True),
    # bare REAL -> single precision loses the attractor geometry
    "lorenz_euler_single.f.json": ("precision_insuficiente", "executed_source",
                                   True),
    # tiny single-precision dt -> orbit collapses (run goes non-finite)
    "lorenz_euler_single_huge.f.json": ("inestabilidad_numerica",
                                        "executed_source", False),
    # C driver runs the legacy source via the strict C subset front-end
    "lorenz_euler_float.c.json": ("precision_insuficiente",
                                  "executed_source", True),
    # layer 10: array-valued RK4 (structured path) is still reference-grade
    "lorenz_rk4_array.f.json": ("equivalente", "executed_source", True),
    # layer 10: block IF + GOTO euler runs the legacy source and degrades
    "lorenz_euler_goto.f.json": ("degradacion_de_orden", "executed_source",
                                 True),
}


def main():
    checks, failed = 0, []

    def check(name, got, want):
        nonlocal checks
        checks += 1
        if got != want:
            failed.append(f"{name}: got {got!r}, expected {want!r}")

    for fn, (want_diag, want_status) in DIAGNOSE.items():
        p = ROOT / "results" / "diagnose" / fn
        d = json.loads(p.read_text())
        check(f"{fn} diagnostico", d["diagnostico"], want_diag)
        if want_status:
            check(f"{fn} status", d.get("execution", {}).get("status"),
                  want_status)

    for fn, (want_diag, want_status, want_finite) in LEGACY.items():
        p = ROOT / "results" / "legacy" / fn
        d = json.loads(p.read_text())
        # status + translator folded into one assertion (C uses the
        # strict_c_subset front-end, the .f sources strict_f77_subset)
        want_tr = ("strict_c_subset" if fn.endswith(".c.json")
                   else "strict_f77_subset")
        check(f"{fn} diagnostico", d["diagnosis"]["diagnostico"], want_diag)
        check(f"{fn} status", (d["execution"]["status"],
                               d["execution"].get("translator")),
              (want_status, want_tr))
        if want_finite is not None:
            check(f"{fn} finite", d["diagnosis"]["fingerprint"]["finite"],
                  want_finite)

    for line in failed:
        print("FAIL", line)
    print(f"{checks - len(failed)}/{checks} assertions passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
