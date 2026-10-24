"""Bundle language runs, hub comparison, registries and discrepancy
analysis into one HDF5 file (results/registry_all.h5).

Layout:
  /meta                              attrs + source file list
  /runs/{octave,python,fortran}/<entry>
      attrs: family, passed, lyap1_observed
      /checks/<quantity>  scalar f8 + attrs (reference, tol, rel_err, pass)
      orbit               f8 matrix when an orbit text file exists
  /hub_compare/<entry>    attrs: status, per-language pass, observed values
  /discrepancy/per_quantity   column table (one dataset per column)
  /discrepancy/sensitivity_aizawa  column table + baseline/language attrs

Usage:
  .venv/bin/python python/export_hdf5.py [output.h5]
Requires h5py (and numpy); exits 2 with a note if they are missing.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError:
    print("SKIP HDF5 export: h5py/numpy not available "
          "(use .venv/bin/python python/export_hdf5.py)")
    sys.exit(2)

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
DISC = RESULTS / "discrepancy"
LANGS = {
    "octave": RESULTS / "registry_run_20261023T230731.json",
    "python": RESULTS / "registry_run_py_20261023T230822.json",
    "fortran": RESULTS / "registry_run_f_20261023T230828.json",
}
HUB = RESULTS / "hub_compare.json"
REG_ATR = ROOT / "registries" / "attractors.json"
REG_POT = ROOT / "registries" / "potentials.json"


def f8(x):
    if x is None:
        return float("nan")
    return float(x)


def write_meta(h5, sources):
    h5.create_dataset("meta/source_files",
                      data=np.array(sources, dtype="S128"))


def write_runs(h5, sources):
    for lang, path in LANGS.items():
        if not path.exists():
            print(f"missing run file for {lang}: {path}")
            continue
        with open(path) as f:
            data = json.load(f)
        src = h5.create_group(f"runs/{lang}")
        src.attrs["timestamp"] = data.get("timestamp", "")
        src.attrs["selector"] = data.get("selector", "")
        src.attrs["language"] = data.get("language", lang)
        src.attrs["n_pass"] = int(data.get("n_pass", -1))
        src.attrs["n_total"] = int(data.get("n_total", -1))
        src.attrs["source_file"] = str(path.relative_to(ROOT))
        sources.append(str(path.relative_to(ROOT)))
        for run in data.get("runs", []):
            g = src.create_group(run["name"])
            g.attrs["family"] = run.get("family", "")
            g.attrs["passed"] = bool(run.get("pass", False))
            g.attrs["lyap1_observed"] = f8(run.get("lyap1_observed"))
            checks = run.get("checks", [])
            if isinstance(checks, dict):
                checks = [checks]
            for chk in checks:
                cg = g.create_group(f"checks/{chk['quantity']}")
                cg.create_dataset("computed", data=f8(chk.get("computed")))
                cg.attrs["reference"] = f8(chk.get("reference"))
                cg.attrs["tol"] = f8(chk.get("tol"))
                cg.attrs["rel_err"] = f8(chk.get("rel_err"))
                cg.attrs["pass"] = bool(chk.get("pass", False))
            opath = run.get("orbit_path")
            if opath and Path(opath).exists():
                try:
                    orbit = np.loadtxt(opath, ndmin=1)
                    g.create_dataset("orbit", data=orbit)
                    g.attrs["orbit_source"] = str(opath)
                except (OSError, ValueError) as e:
                    print(f"orbit unreadable for {lang}/{run['name']}: {e}")


def write_hub(h5, sources):
    if not HUB.exists():
        print(f"missing {HUB}")
        return
    with open(HUB) as f:
        hub = json.load(f)
    g = h5.create_group("hub_compare")
    g.attrs["all_pass"] = bool(hub.get("all_pass", False))
    g.attrs["n_entries"] = int(hub.get("n_entries", -1))
    g.attrs["source_file"] = str(HUB.relative_to(ROOT))
    sources.append(str(HUB.relative_to(ROOT)))
    for row in hub.get("rows", []):
        rg = g.create_group(row["name"])
        rg.attrs["status"] = row.get("status", "")
        passes = row.get("passes", {})
        for lang in LANGS:
            rg.attrs[f"passes_{lang}"] = bool(passes.get(lang, False))
        rg.attrs["obs_rel_diff"] = f8(row.get("obs_rel_diff"))
        obs = row.get("obs", {})
        for lang in LANGS:
            rg.attrs[f"obs_{lang}"] = f8(obs.get(lang))


def write_per_quantity(h5, sources):
    path = DISC / "per_quantity.json"
    if not path.exists():
        print(f"missing {path}")
        return
    with open(path) as f:
        pq = json.load(f)
    rows = pq["rows"]
    g = h5.create_group("discrepancy/per_quantity")
    g.attrs["n_quantities"] = int(pq.get("n_quantities", len(rows)))
    g.attrs["n_reference_backed"] = int(pq.get("n_reference_backed", -1))
    g.attrs["n_observed_only"] = int(pq.get("n_observed_only", -1))
    g.attrs["reference_backed_spread_max"] = f8(
        (pq.get("reference_backed_spread") or {}).get("max"))
    g.attrs["reference_backed_spread_median"] = f8(
        (pq.get("reference_backed_spread") or {}).get("median"))
    g.attrs["observed_spread_max"] = f8(
        (pq.get("observed_spread") or {}).get("max"))
    g.attrs["observed_spread_median"] = f8(
        (pq.get("observed_spread") or {}).get("median"))
    g.attrs["residual_abs_max"] = f8(pq.get("residual_abs_max"))
    g.attrs["source_file"] = str(path.relative_to(ROOT))
    sources.append(str(path.relative_to(ROOT)))
    n = len(rows)
    enc = lambda key: np.array([r.get(key, "") for r in rows], dtype="S64")
    g.create_dataset("entry", data=enc("entry"))
    g.create_dataset("quantity", data=enc("quantity"))
    g.create_dataset("family", data=enc("family"))
    g.create_dataset("metric", data=enc("metric"))
    g.create_dataset("spread", data=np.array([f8(r.get("spread")) for r in rows]))
    g.create_dataset("spread_abs",
                     data=np.array([f8(r.get("spread_abs")) for r in rows]))
    g.create_dataset("tol", data=np.array([f8(r.get("tol")) for r in rows]))
    g.create_dataset("reference",
                     data=np.array([f8(r.get("reference")) for r in rows]))
    g.create_dataset("observed_only",
                     data=np.array([bool(r.get("observed_only")) for r in rows]))
    g.create_dataset("all_pass",
                     data=np.array([bool(r.get("all_pass")) for r in rows]))
    for lang in LANGS:
        g.create_dataset(
            f"computed_{lang}",
            data=np.array([f8((r.get("computed") or {}).get(lang)) for r in rows]))
        g.create_dataset(
            f"passes_{lang}",
            data=np.array([bool((r.get("passes") or {}).get(lang, False))
                           for r in rows]))
    g.attrs["n_rows"] = n


def write_sensitivity(h5, sources):
    path = DISC / "sensitivity_aizawa.json"
    if not path.exists():
        print(f"missing {path}")
        return
    with open(path) as f:
        sens = json.load(f)
    g = h5.create_group("discrepancy/sensitivity_aizawa")
    g.attrs["baseline_lyap1"] = f8(sens.get("baseline_lyap1"))
    lang = sens.get("lang_values", {})
    for lg in LANGS:
        g.attrs[f"lang_{lg}"] = f8(lang.get(lg))
    concl = sens.get("conclusions")
    if concl is not None:
        g.attrs["conclusions_json"] = json.dumps(concl)
    g.attrs["source_file"] = str(path.relative_to(ROOT))
    sources.append(str(path.relative_to(ROOT)))
    rows = sens.get("rows", [])
    g.create_dataset("kind",
                     data=np.array([r["kind"] for r in rows], dtype="S16"))
    g.create_dataset("value",
                     data=np.array([f8(r.get("value")) for r in rows]))
    g.create_dataset("lyap1",
                     data=np.array([f8(r.get("lyap1")) for r in rows]))
    g.attrs["n_rows"] = len(rows)


def main(out_path):
    sources = []
    for p in (REG_ATR, REG_POT):
        if p.exists():
            sources.append(str(p.relative_to(ROOT)))
    for p in list(LANGS.values()) + [HUB, DISC / "per_quantity.json",
                                     DISC / "sensitivity_aizawa.json"]:
        if p.exists():
            sources.append(str(p.relative_to(ROOT)))
    with h5py.File(out_path, "w") as h5:
        h5.attrs["schema"] = "registry_all/1"
        h5.attrs["created"] = datetime.now(timezone.utc).isoformat()
        write_meta(h5, sources)
        write_runs(h5, [])
        write_hub(h5, [])
        write_per_quantity(h5, [])
        write_sensitivity(h5, [])
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTS / "registry_all.h5"
    main(out)
