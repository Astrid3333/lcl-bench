"""Cross-language discrepancy analysis: per-quantity spread across the three
latest runner result files (Octave, Python, Fortran).

Writes results/discrepancy/per_quantity.json and per_quantity.csv with one row
per (entry, quantity): computed value per language, reference, tolerance,
pass flags, and the relative spread (max-min)/max|value| over the languages
that produced a number. Observed-only quantities (no citable reference) are
tagged so the report can separate estimator noise from reference-backed checks.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
OUT = RESULTS / "discrepancy"
LANGS = ("octave", "python", "fortran")


def latest_run(lang: str) -> Path:
    if lang == "python":
        cands = sorted(RESULTS.glob("registry_run_py_*.json"))
    elif lang == "fortran":
        cands = sorted(RESULTS.glob("registry_run_f_*.json"))
    else:
        cands = sorted(p for p in RESULTS.glob("registry_run_*.json")
                       if "_py_" not in p.name and "_f_" not in p.name)
    if not cands:
        raise FileNotFoundError(f"no registry_run result for {lang}")
    return cands[-1]


def as_list(checks):
    if isinstance(checks, dict):
        return [checks]
    return checks or []


def spread_metrics(vals):
    """Return (spread_rel, spread_abs, metric). Relative spread is only
    meaningful when the values are not near zero; residual-type quantities
    (orbit_head_maxabs ~ 1e-16) fall back to an absolute metric."""
    vals = [v for v in vals if v is not None]
    if len(vals) < 2:
        return None, None, None
    hi, lo = max(vals), min(vals)
    span = hi - lo
    den = max(abs(hi), abs(lo))
    if den < 1e-10:
        return None, span, "abs_residual"
    return span / den, span, "rel"


def main():
    files = {lang: latest_run(lang) for lang in LANGS}
    idx = {}
    for lang, path in files.items():
        with open(path) as f:
            run = json.load(f)
        idx[lang] = {r["name"]: r for r in run["runs"]}

    names = sorted(set().union(*(set(i) for i in idx.values())))
    rows = []
    for name in names:
        fams = [idx[l][name]["family"] for l in LANGS if name in idx[l]]
        family = fams[0] if fams else "unknown"
        per_q = {}
        for lang in LANGS:
            rec = idx[lang].get(name)
            if rec is None:
                continue
            for c in as_list(rec.get("checks")):
                q = c["quantity"]
                slot = per_q.setdefault(q, {
                    "quantity": q, "entry": name, "family": family,
                    "computed": {}, "passes": {}, "reference": c.get("reference"),
                    "tol": c.get("tol"),
                })
                slot["computed"][lang] = c.get("computed")
                slot["passes"][lang] = c.get("pass")
                if slot["reference"] is None:
                    slot["reference"] = c.get("reference")
                if slot["tol"] is None:
                    slot["tol"] = c.get("tol")
        for q, slot in per_q.items():
            vals = [slot["computed"].get(l) for l in LANGS]
            rel, abs_, metric = spread_metrics(vals)
            slot["spread"] = rel
            slot["spread_abs"] = abs_
            slot["metric"] = metric
            slot["observed_only"] = str(q).endswith("_observed")
            slot["all_pass"] = all(
                slot["passes"].get(l) is True
                for l in LANGS if l in slot["passes"]
            )
            rows.append(slot)

    rows.sort(key=lambda r: (r["family"], r["entry"], r["quantity"]))
    ref_spreads = [r["spread"] for r in rows
                   if not r["observed_only"] and r["spread"] is not None]
    obs_spreads = [r["spread"] for r in rows
                   if r["observed_only"] and r["spread"] is not None]
    resid = [r["spread_abs"] for r in rows if r["metric"] == "abs_residual"]
    summary = {
        "files": {l: files[l].name for l in LANGS},
        "n_quantities": len(rows),
        "n_reference_backed": sum(1 for r in rows if not r["observed_only"]),
        "n_observed_only": sum(1 for r in rows if r["observed_only"]),
        "reference_backed_spread": {
            "max": max(ref_spreads) if ref_spreads else None,
            "median": sorted(ref_spreads)[len(ref_spreads) // 2] if ref_spreads else None,
        },
        "observed_spread": {
            "max": max(obs_spreads) if obs_spreads else None,
            "median": sorted(obs_spreads)[len(obs_spreads) // 2] if obs_spreads else None,
        },
        "residual_abs_max": max(resid) if resid else None,
        "rows": rows,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "per_quantity.json", "w") as f:
        json.dump(summary, f, indent=2)
    with open(OUT / "per_quantity.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["entry", "family", "quantity", "observed_only",
                    "octave", "python", "fortran", "reference", "tol",
                    "spread", "spread_abs", "metric", "all_pass"])
        for r in rows:
            w.writerow([
                r["entry"], r["family"], r["quantity"], r["observed_only"],
                r["computed"].get("octave"), r["computed"].get("python"),
                r["computed"].get("fortran"), r["reference"], r["tol"],
                r["spread"], r["spread_abs"], r["metric"], r["all_pass"],
            ])
    print(f"quantities: {len(rows)} "
          f"(reference-backed {summary['n_reference_backed']}, "
          f"observed {summary['n_observed_only']})")
    print(f"reference-backed spread: max={summary['reference_backed_spread']['max']:.3g} "
          f"median={summary['reference_backed_spread']['median']:.3g}")
    print(f"observed spread:         max={summary['observed_spread']['max']:.3g} "
          f"median={summary['observed_spread']['median']:.3g}")
    if summary["residual_abs_max"] is not None:
        print(f"orbit-head residual:     max_abs={summary['residual_abs_max']:.3g}")
    print(f"Saved: {OUT / 'per_quantity.json'}")
    print(f"Saved: {OUT / 'per_quantity.csv'}")


if __name__ == "__main__":
    main()
