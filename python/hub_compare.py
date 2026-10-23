"""Cross-language hub: compare Octave and Python runs over the same registry.

Both languages independently evaluate the same registry JSON and compare
computed invariants against the same analytic/literature references (not
against each other). This hub then diffs the two result sets entry by entry
to report agreement between languages.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"


def latest_run(lang: str) -> Path | None:
    if lang == "python":
        cands = sorted(RESULTS.glob("registry_run_py_*.json"))
    else:
        # exclude Python results (they contain "_py_")
        cands = sorted(p for p in RESULTS.glob("registry_run_*.json")
                       if "_py_" not in p.name)
    return cands[-1] if cands else None


def index_runs(run: dict) -> dict[str, dict]:
    return {r["name"]: r for r in run["runs"]}


def _checks_as_list(checks) -> list[dict]:
    # Octave jsonencode serializes a 1-element struct array as a single dict.
    if isinstance(checks, dict):
        return [checks]
    return checks


def _observed(checks) -> float | None:
    for c in _checks_as_list(checks):
        if str(c.get("quantity", "")).endswith("_observed"):
            return c.get("computed")
    return None


def compare(oct_path: Path, py_path: Path) -> dict:
    with open(oct_path) as f:
        oct_run = json.load(f)
    with open(py_path) as f:
        py_run = json.load(f)
    oct_idx = index_runs(oct_run)
    py_idx = index_runs(py_run)

    names = sorted(set(oct_idx) | set(py_idx))
    rows = []
    both_pass = 0
    for name in names:
        o = oct_idx.get(name)
        p = py_idx.get(name)
        if o is None or p is None:
            rows.append({"name": name, "status": "MISSING",
                         "octave": o is not None, "python": p is not None})
            continue
        o_pass = o["pass"]
        p_pass = p["pass"]
        obs_o = _observed(o["checks"])
        obs_p = _observed(p["checks"])
        agree = None
        if obs_o is not None and obs_p is not None:
            denom = max(abs(obs_o), abs(obs_p), 1e-15)
            agree = abs(obs_o - obs_p) / denom
        status = "BOTH_PASS" if (o_pass and p_pass) else (
            "OCTAVE_ONLY" if o_pass else (
                "PYTHON_ONLY" if p_pass else "BOTH_FAIL"))
        if status == "BOTH_PASS":
            both_pass += 1
        rows.append({
            "name": name, "status": status,
            "octave_pass": o_pass, "python_pass": p_pass,
            "obs_rel_diff": agree,
        })

    summary = {
        "octave_file": oct_path.name, "python_file": py_path.name,
        "n_entries": len(names), "both_pass": both_pass,
        "rows": rows,
    }
    out_path = RESULTS / "hub_compare.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    return summary


def main():
    oct_path = latest_run("octave")
    py_path = latest_run("python")
    if not oct_path or not py_path:
        print("Need both an Octave and a Python registry_run result in results/",
              file=sys.stderr)
        sys.exit(1)
    print(f"Octave: {oct_path.name}")
    print(f"Python: {py_path.name}")
    s = compare(oct_path, py_path)
    print(f"\nEntries: {s['n_entries']}  Both PASS: {s['both_pass']}")
    for r in s["rows"]:
        extra = ""
        if r.get("obs_rel_diff") is not None:
            extra = f"  obs_rel_diff={r['obs_rel_diff']:.3g}"
        print(f"  {r['name']:16s} {r['status']}{extra}")
    print(f"\nSaved: {RESULTS / 'hub_compare.json'}")


if __name__ == "__main__":
    main()
