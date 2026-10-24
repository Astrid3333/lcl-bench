"""Finite-time Lyapunov sensitivity experiment for the aizawa ODE.

Reproduces the cross-language gap diagnosis: the vector field and Benettin
algorithm are identical across Octave/Python/Fortran, but the finite-time
estimate shifts by ~13% for a one-ulp change of x0 and by ~7% when t_lyap is
halved/doubled, while being invariant to d0 over four decades. This is
estimator sensitivity to chaotic amplification, not an implementation bug.

Writes results/discrepancy/sensitivity_aizawa.json. Run with the same Python
environment as registry_run.py (system python3, numpy 1.26.4).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))

import handles  # noqa: E402
import kernels  # noqa: E402

REG = ROOT / "registries" / "attractors.json"
OUT = ROOT / "results" / "discrepancy" / "sensitivity_aizawa.json"

# Published cross-language values for comparison (committed suite runs).
LANG_VALUES = {
    "octave": 0.10733175348719398,
    "python": 0.10733175348719351,
    "fortran": 0.09171612312383404,
}


def load_entry():
    with open(REG) as f:
        reg = json.load(f)
    for e in reg["entries"]:
        if e["name"] == "aizawa":
            return e
    raise KeyError("aizawa not in registry")


def lyap(entry, x0=None, d0=1e-8, t_lyap=None):
    p = entry["params"]
    x0 = np.array(entry["x0"] if x0 is None else x0, dtype=float)
    dt = entry["dt_lyap"]
    t_tr = entry["t_transient"]
    t_ly = entry["t_lyap"] if t_lyap is None else t_lyap
    f = handles.rhs(entry["name"], p)
    return float(kernels.lyap_ode(f, x0, dt, t_tr, t_ly, d0=d0))


def main():
    entry = load_entry()
    base = lyap(entry)
    rows = []

    def add(kind, changed, value, lam):
        rows.append({
            "kind": kind,
            "changed": changed,
            "value": value,
            "lyap1": lam,
            "delta_vs_baseline": lam - base,
            "rel_vs_baseline": (lam - base) / base,
        })

    add("baseline", "registry defaults", None, base)

    for d0 in (1e-10, 1e-9, 1e-7, 1e-6):
        add("d0", "d0", d0, lyap(entry, d0=d0))

    x0 = entry["x0"]
    for eps in (1e-16, 1e-15, 1e-14, 1e-13, 1e-12):
        add("x0", "x0[0] + eps", eps, lyap(entry, x0=[x0[0] + eps, x0[1], x0[2]]))

    for t in (500, 1000, 4000, 8000):
        add("t_lyap", "t_lyap", t, lyap(entry, t_lyap=t))

    d0_spread = max(r["lyap1"] for r in rows if r["kind"] == "d0") - \
        min(r["lyap1"] for r in rows if r["kind"] == "d0")
    x0_rels = [abs(r["rel_vs_baseline"]) for r in rows if r["kind"] == "x0"]
    t_rels = [abs(r["rel_vs_baseline"]) for r in rows if r["kind"] == "t_lyap"]

    out = {
        "entry": "aizawa",
        "dt_lyap": entry["dt_lyap"],
        "t_transient": entry["t_transient"],
        "t_lyap": entry["t_lyap"],
        "d0_default": 1e-8,
        "baseline_lyap1": base,
        "lang_values": LANG_VALUES,
        "lang_spread": (max(LANG_VALUES.values()) - min(LANG_VALUES.values())) / max(
            abs(v) for v in LANG_VALUES.values()),
        "conclusions": {
            "d0_range": [1e-10, 1e-6],
            "d0_spread_abs": d0_spread,
            "x0_max_abs_rel_shift": max(x0_rels),
            "t_lyap_max_abs_rel_shift": max(t_rels),
            "one_ulp_rel_shift": next(
                r["rel_vs_baseline"] for r in rows
                if r["kind"] == "x0" and r["value"] == 1e-16),
        },
        "rows": rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"baseline lyap1 = {base:.15f}")
    print(f"d0 spread (1e-10..1e-6)      : {d0_spread:.3e}")
    print(f"x0 eps max |rel shift|        : {max(x0_rels):.3%}")
    print(f"t_lyap 500..8000 max |rel|    : {max(t_rels):.3%}")
    print(f"lang spread (oct/py/f vs fort): {out['lang_spread']:.3%}")
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()
