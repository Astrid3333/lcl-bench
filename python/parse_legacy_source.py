"""Static reader for legacy Fortran 77 / C attractor sources.

Reads an old .f/.f77/.c file and extracts, WITH line evidence:
  language, precision (single/double), integrator (rk4/euler), timestep,
  recognized system (registry name). It then maps that metadata onto the
  equivalent synthetic variant of diagnose_attractor.py and returns the
  calibrated fingerprint verdict alongside the extracted facts.

HONESTY: an F77 or C file inside the supported subsets is EXECUTED
(translated to Python by translate_f77.py / translate_c.py, then run under
the precision declared in the source), so lambda1/D2 come from the LEGACY
recurrence itself and the report says `candidate_source: executed_source`.
Anything the translator refuses (several loops, non-literal array extents,
labelled DO, while-loops, arrays in C, ...) falls back to the
synthetic equivalent variant (registry system run under the integrator /
precision / dt read from the source) and is labelled
`candidate_source: synthetic_equivalent` plus an `execution.status` block
(`executed_source` | `translation_failed` | `unsupported_language`)
listing every assumption made to bridge source metadata to a runnable
variant; `unsupported_language` is reserved for languages other than F77/C.

Usage:
  python3 parse_legacy_source.py fixtures/legacy/lorenz_euler_single.f
  python3 parse_legacy_source.py fixtures/legacy/lorenz_euler_single.f --out r.json
"""
import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from diagnose_attractor import (diagnose, load_entry, oracle_fingerprint,
                                report_for)
from translate_c import fingerprint_c_source
from translate_f77 import TranslationError, fingerprint_source

DT_NAMES = ("dt", "delt", "delta", "timestep", "dtime", "tau", "step", "h")
SYSTEM_KEYWORDS = (
    "lorenz", "rossler", "thomas", "aizawa", "henon", "logistic",
    "clifford", "tinkerbell", "barnsley", "sierpinski", "duffing",
    "fitzhugh", "lozi", "chirikov",
)
CONST_SETS = ((("sigma", "rho", "beta"), "lorenz"),)


def _lines(text):
    return text.splitlines()


def detect_language(path, text):
    suf = Path(path).suffix.lower()
    if suf in (".f", ".f77", ".f90", ".for", ".ftn"):
        return "fortran"
    if suf in (".c", ".h", ".cpp", ".cc"):
        return "c"
    if re.search(r"^\s*(program|subroutine|end\s+do)\b", text, re.I | re.M):
        return "fortran"
    if "#include" in text:
        return "c"
    return "unknown"


def extract_precision(text, lang):
    ev = []
    for i, line in enumerate(_lines(text), 1):
        low = line.lower()
        if lang == "fortran":
            if re.search(r"\bdouble\s+precision\b", low) or \
                    re.search(r"\breal\s*\*\s*8\b", low) or \
                    re.search(r"\breal\s*\(\s*kind\s*=\s*8\s*\)", low) or \
                    re.search(r"\breal\s*\(\s*8\s*\)", low):
                ev.append({"line": i, "text": line.strip(), "marks": "double"})
            elif re.search(r"\breal\s*\*\s*4\b", low) or \
                    re.search(r"\breal\s*\(\s*kind\s*=\s*4\s*\)", low) or \
                    re.search(r"\breal\b", low):  # bare F77 REAL = single
                ev.append({"line": i, "text": line.strip(), "marks": "single"})
        else:
            if re.search(r"\bdouble\s+(precision\s+)?\w", low):
                ev.append({"line": i, "text": line.strip(), "marks": "double"})
            elif re.search(r"\bfloat\b", low):
                ev.append({"line": i, "text": line.strip(), "marks": "single"})
    counts = Counter(e["marks"] for e in ev)
    if not counts:
        return {"value": "unknown", "mixed": False, "evidence": []}
    value, _ = counts.most_common(1)[0]
    return {"value": value, "mixed": len(counts) > 1, "counts": dict(counts),
            "evidence": ev[:6]}


def extract_integrator(text, lang):
    low = text.lower()
    if all(k in low for k in ("k1", "k2", "k3", "k4")):
        line = next(i for i, l in enumerate(_lines(text), 1) if "k1" in l.lower())
        snippet = _lines(text)[line - 1].strip()
        return {"value": "rk4",
                "evidence": [{"line": line, "text": snippet,
                              "why": "k1..k4 stage variables present"}]}
    pat = re.compile(r"\b(\w+)\s*=\s*\1\s*\+\s*[\w.]+\s*\*\s*[\w.]+")
    pat_c = re.compile(r"\b(\w+)\s*\+=\s*[\w.]+\s*\*\s*[\w.]+")
    ev = []
    for i, line in enumerate(_lines(text), 1):
        if pat.search(line) or pat_c.search(line):
            ev.append({"line": i, "text": line.strip(),
                       "why": "self-assignment with product step (x = x + dt*dx)"})
    if ev:
        return {"value": "euler", "evidence": ev[:4]}
    return {"value": "unknown", "evidence": []}


def _num(lit):
    return float(lit.replace("d", "e").replace("D", "e"))


def extract_dt(text, lang):
    name = "|".join(DT_NAMES)
    pat = re.compile(rf"\b({name})\s*=\s*([0-9]*\.?[0-9]+"
                     rf"(?:[eEdD][+-]?\d+)?)[fF]?\b", re.I)
    matches = []
    for i, line in enumerate(_lines(text), 1):
        for m in pat.finditer(line):
            v = _num(m.group(2))
            if 0.0 < v < 100.0:
                matches.append({"line": i, "text": line.strip(),
                                "name": m.group(1).lower(), "value": v})
    if not matches:
        return {"value": None, "evidence": []}
    chosen = next((m for m in matches if m["name"] == "dt"), matches[0])
    return {"value": chosen["value"], "evidence": matches[:6],
            "chosen_line": chosen["line"]}


def detect_system(text):
    ev = []
    for i, line in enumerate(_lines(text), 1):
        low = line.lower()
        for kw in SYSTEM_KEYWORDS:
            if re.search(rf"\b{kw}\b", low):
                ev.append({"line": i, "text": line.strip(), "why": f"keyword {kw!r}"})
                return {"name": kw, "evidence": ev}
    for i, line in enumerate(_lines(text), 1):
        low = line.lower()
        for names, name in CONST_SETS:
            if all(re.search(rf"\b{n}\s*=", low) for n in names):
                ev.append({"line": i, "text": line.strip(),
                           "why": f"parameter set {names}"})
                return {"name": name, "evidence": ev}
    return {"name": None, "evidence": []}


def map_variant(entry, precision, integrator, dt):
    """Bridge source metadata -> diagnose variant. Returns (variant, notes)."""
    notes = []
    if entry is None:
        return None, ["system not recognized: diagnosis skipped"]
    if entry["family"] == "map":
        if precision == "single":
            return "float32", notes
        if precision != "double":
            notes.append(f"precision {precision!r} -> assumed double (reference)")
        return "reference", notes
    if integrator == "euler":
        dt_reg = entry.get("dt_lyap", 0.01)
        if dt is not None and dt >= 5.0 * dt_reg:
            notes.append(f"source dt={dt:g} >= 5x registry dt_lyap={dt_reg:g} "
                         "-> euler_large variant")
            if precision == "single":
                notes.append("single precision not composed into euler_large "
                             "(variant runs float64)")
            return "euler_large", notes
        if precision == "single":
            return "euler_single", notes
        if precision != "double":
            notes.append(f"precision {precision!r} -> assumed double")
        return "euler", notes
    if integrator == "unknown":
        notes.append("integrator not detected -> assumed RK4 (reference)")
    if precision == "single":
        return "float32", notes
    if precision != "double":
        notes.append(f"precision {precision!r} -> assumed double")
    return "reference", notes


def analyze(path):
    path = Path(path)
    text = path.read_text(errors="replace")
    lang = detect_language(path, text)
    prec = extract_precision(text, lang)
    integ = extract_integrator(text, lang)
    dt = extract_dt(text, lang)
    sysm = detect_system(text)
    entry = load_entry(sysm["name"]) if sysm["name"] else None
    variant, notes = map_variant(entry, prec["value"], integ["value"],
                                 dt["value"])
    rep = {
        "file": str(path),
        "language": lang,
        "extracted": {"system": sysm, "precision": prec,
                      "integrator": integ, "dt": dt},
        "variant": variant,
        "candidate_source": "synthetic_equivalent",
        "assumptions": notes,
        "execution": None,
        "diagnosis": None,
    }

    # ---- 1st choice: EXECUTE the legacy source (strict F77 / C subsets) ----
    if lang in ("fortran", "c") and entry is not None:
        try:
            if lang == "fortran":
                cand, info = fingerprint_source(text, entry, h=dt["value"])
            else:
                cand, info = fingerprint_c_source(text, entry, h=dt["value"])
        except TranslationError as e:
            rep["execution"] = {"status": "translation_failed",
                                "error": str(e), "line": e.line,
                                "note": "outside the supported "
                                        f"{'F77' if lang == 'fortran' else 'C'} "
                                        "subset; diagnosis falls back to the "
                                        "synthetic equivalent variant"}
            notes.append(f"translation refused at line {e.line}: {e}")
        else:
            if dt["value"] is not None and entry["family"] == "ode" and \
                    abs(dt["value"] - entry.get("dt_lyap", dt["value"])) > 1e-15:
                notes.append(f"lambda1 normalised with the SOURCE dt="
                             f"{dt['value']:g} (registry dt_lyap="
                             f"{entry.get('dt_lyap'):g}): the legacy program "
                             "was run with its own timestep")
            oracle = oracle_fingerprint(entry)
            rep["execution"] = info
            rep["candidate_source"] = "executed_source"
            rep["diagnosis"] = report_for(entry, oracle, cand, variant,
                                          "executed_source")
            return rep
    elif lang not in ("fortran", "c"):
        rep["execution"] = {"status": "unsupported_language",
                            "note": f"language {lang!r} has no translator: "
                                    "only the synthetic equivalent variant "
                                    "can be diagnosed from its metadata"}

    # ---- fallback: synthetic equivalent variant ----
    if dt["value"] is not None and entry is not None and \
            entry["family"] == "ode" and \
            abs(dt["value"] - entry.get("dt_lyap", dt["value"])) > 1e-15:
        notes.append(f"source dt={dt['value']:g} differs from registry "
                     f"dt_lyap={entry.get('dt_lyap'):g}; diagnosis uses the "
                     "registry dt for comparability")
    if variant:
        rep["diagnosis"] = diagnose(sysm["name"], variant)
    else:
        rep["note"] = ("no runnable equivalent: identify the system manually "
                       "or add it to registries/attractors.json first")
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("file", help="legacy source (.f/.f77/.c)")
    ap.add_argument("--out", default=None, help="write report JSON here")
    args = ap.parse_args(argv)
    rep = analyze(args.file)
    text = json.dumps(rep, indent=2)
    if args.out:
        out = Path(args.out)
        out.write_text(text + "\n")
        src = (rep.get("execution") or {}).get("translated_source")
        if src:
            out.with_suffix(".translated.py").write_text(src)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
