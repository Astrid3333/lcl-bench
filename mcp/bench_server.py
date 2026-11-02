#!/usr/bin/env python3
"""MCP server exposing the cross-language benchmark registry over stdio.

Zero-dependency JSON-RPC 2.0 (newline-delimited) MCP server. Wraps the
existing CLI entry points (registry_run, diagnose_attractor, hub_compare,
check_verdicts) so an AI client can inspect the registry, run single entries
in any of the three languages, and diagnose precision/order degradations.

Usage: python3 mcp/bench_server.py   (stdio; no arguments)
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
REGISTRIES = ROOT / "registries"
TIMEOUT_S = 180
PROTOCOL_VERSION = "2024-11-05"
LANG_GLOBS = {
    "python": "registry_run_py_*.json",
    "fortran": "registry_run_f_*.json",
    "octave": None,
}
VARIANTS = ("reference", "float32", "euler", "euler_large", "euler_single")


def _load_entries() -> list[dict]:
    out = []
    for reg in ("attractors.json", "potentials.json", "energy_drift.json",
                "goertzel_fft.json", "archimedes_pi.json", "middle_square_weyl.json",
                "hybrid_prng.json"):
        with open(REGISTRIES / reg) as f:
            out.extend(json.load(f)["entries"])
    return out


def _entry_summary(e: dict) -> dict:
    return {
        "name": e.get("name"),
        "family": e.get("family"),
        "dim": e.get("dim"),
        "ref": sorted((e.get("ref") or {}).keys()),
        "notes": e.get("notes"),
    }


def _newest_run_json(lang: str, since: float) -> Path | None:
    if lang == "octave":
        cands = [p for p in RESULTS.glob("registry_run_*.json")
                 if "_py_" not in p.name and "_f_" not in p.name]
    else:
        cands = list(RESULTS.glob(LANG_GLOBS[lang]))
    fresh = [p for p in cands if p.stat().st_mtime >= since]
    if not fresh:
        return None
    return max(fresh, key=lambda p: p.stat().st_mtime)


def _checks_as_list(checks) -> list[dict]:
    if isinstance(checks, dict):
        return [checks]
    return checks or []


def _run_cmd(cmd: list[str], timeout: int = TIMEOUT_S) -> dict:
    t0 = time.time()
    proc = subprocess.run(
        cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout,
    )
    return {
        "cmd": cmd,
        "returncode": proc.returncode,
        "stdout": proc.stdout[-8000:],
        "stderr": proc.stderr[-4000:],
        "wall_s": round(time.time() - t0, 2),
    }


def tool_registry_list(args: dict) -> dict:
    return {
        "n_entries": len(_load_entries()),
        "entries": [_entry_summary(e) for e in _load_entries()],
    }


def tool_registry_get(args: dict) -> dict:
    name = args.get("name")
    for e in _load_entries():
        if e.get("name") == name:
            return e
    raise ValueError(f"entry not found: {name!r}")


def tool_run_entry(args: dict) -> dict:
    name = args.get("name")
    language = args.get("language")
    if language not in ("octave", "python", "fortran"):
        raise ValueError(f"language must be octave|python|fortran, got {language!r}")
    known = {e["name"] for e in _load_entries()}
    if name not in known:
        raise ValueError(f"entry not found: {name!r}")
    if language == "python":
        cmd = [sys.executable, "python/registry_run.py", name]
    elif language == "octave":
        cmd = ["octave-cli", "--no-gui", "--quiet",
               "--eval", f"addpath('octave'); registry_run('{name}');"]
    else:
        cmd = ["./fortran/f_runner", name]
    since = time.time()
    res = _run_cmd(cmd)
    run_path = _newest_run_json(language, since)
    out = {
        "language": language,
        "entry": name,
        "returncode": res["returncode"],
        "wall_s": res["wall_s"],
        "run_file": run_path.name if run_path else None,
        "stderr_tail": res["stderr"][-500:] if res["stderr"] else "",
    }
    if run_path:
        with open(run_path) as f:
            run = json.load(f)
        rows = [r for r in run.get("runs", []) if r.get("name") == name]
        if rows:
            r0 = rows[0]
            out["pass"] = r0.get("pass")
            out["checks"] = _checks_as_list(r0.get("checks"))
            # Campos adicionales no estandar (ej. middle_square_cycles, weyl_period, etc.)
            std = {"name", "family", "pass", "checks", "error", "orbit_path",
                   "lyap1_observed", "dim_observed", "pi_stable", "pi_naive"}
            for k, v in r0.items():
                if k not in std and v is not None:
                    out[k] = v
    return out


def tool_diagnose(args: dict) -> dict:
    system = args.get("system")
    if not system:
        raise ValueError("system is required")
    variant = args.get("variant", "reference")
    series = args.get("series")
    series_dt = args.get("series_dt")
    if variant not in VARIANTS:
        raise ValueError(f"variant must be one of {VARIANTS}, got {variant!r}")
    if series:
        cmd = [sys.executable, "python/diagnose_attractor.py",
               "--system", system, "--series", str(series)]
        if series_dt is not None:
            cmd += ["--series-dt", str(series_dt)]
    else:
        cmd = [sys.executable, "python/diagnose_attractor.py",
               "--system", system, "--variant", variant]
    res = _run_cmd(cmd)
    if res["returncode"] != 0:
        raise RuntimeError(
            f"diagnose failed (rc={res['returncode']}): {res['stderr'][-500:]}")
    return json.loads(res["stdout"])


def tool_hub(args: dict) -> dict:
    res = _run_cmd([sys.executable, "python/hub_compare.py"])
    if res["returncode"] != 0:
        raise RuntimeError(
            f"hub_compare failed (rc={res['returncode']}): {res['stderr'][-500:]}")
    with open(RESULTS / "hub_compare.json") as f:
        s = json.load(f)
    bad = [r for r in s["rows"] if r["status"] != "ALL_PASS"]
    return {
        "n_entries": s["n_entries"],
        "all_pass": s["all_pass"],
        "files": s["files"],
        "non_all_pass_rows": bad,
        "wall_s": res["wall_s"],
    }


def tool_check_verdicts(args: dict) -> dict:
    res = _run_cmd([sys.executable, "python/check_verdicts.py"])
    lines = [ln for ln in res["stdout"].splitlines() if ln.strip()]
    return {
        "returncode": res["returncode"],
        "summary": lines[-1] if lines else "",
        "failures": [ln for ln in lines if ln.startswith("FAIL")],
        "wall_s": res["wall_s"],
    }


def tool_legacy_analyze(args: dict) -> dict:
    file = args.get("file")
    if not file:
        raise ValueError("file is required")
    path = Path(file)
    if not path.is_absolute():
        path = ROOT / path
    if not path.is_file():
        raise ValueError(f"file not found: {path}")
    py_dir = str(ROOT / "python")
    if py_dir not in sys.path:
        sys.path.insert(0, py_dir)
    import parse_legacy_source
    t0 = time.time()
    rep = parse_legacy_source.analyze(path)
    rep["wall_s"] = round(time.time() - t0, 2)
    return rep


TOOLS = [
    {
        "name": "registry_list",
        "description": (
            "List every registry entry (name, family, dim, reference "
            "quantities) of the cross-language benchmark. Fast, read-only."),
        "inputSchema": {"type": "object", "properties": {}, "required": []},
        "func": tool_registry_list,
    },
    {
        "name": "registry_get",
        "description": (
            "Return the full registry entry JSON for one system "
            "(equations, params, x0, references, notes)."),
        "inputSchema": {
            "type": "object",
            "properties": {"name": {"type": "string"}},
            "required": ["name"],
        },
        "func": tool_registry_get,
    },
    {
        "name": "run_entry",
        "description": (
            "Evaluate one registry entry in one language (octave, python or "
            "fortran) and return its checks (computed vs reference) and pass "
            "flag. Slow entries (fhn, lorenz ODE) can take tens of seconds. "
            "Note: this writes a new results file, so it becomes that "
            "language's latest run and changes what the hub tool compares."),
        "inputSchema": {
            "type": "object",
            "properties": {
                "name": {"type": "string"},
                "language": {"type": "string",
                             "enum": ["octave", "python", "fortran"]},
            },
            "required": ["name", "language"],
        },
        "func": tool_run_entry,
    },
    {
        "name": "diagnose",
        "description": (
            "Run the cross-language diagnosis on one system: compare a "
            "candidate (registry variant or external series file) against "
            "the registry fingerprint and return the diagnostic verdict "
            "(equivalente, colapso_por_redondeo, precision_insuficiente, "
            "artefacto_tipo_dato, degradacion_de_orden, "
            "diferencia_por_integrador, degradacion_de_precision, "
            "inestabilidad_numerica, diferencia_no_explicada)."),
        "inputSchema": {
            "type": "object",
            "properties": {
                "system": {"type": "string"},
                "variant": {"type": "string", "enum": list(VARIANTS)},
                "series": {
                    "type": "string",
                    "description": "path to a trajectory file (npy/txt, rows=states)"},
                "series_dt": {
                    "type": "number",
                    "description": "sample spacing of --series"},
            },
            "required": ["system"],
        },
        "func": tool_diagnose,
    },
    {
        "name": "hub",
        "description": (
            "Run the cross-language hub: compare the latest Octave, Python "
            "and Fortran registry runs entry by entry and summarise "
            "agreement (n_entries, all_pass, rows that are not ALL_PASS)."),
        "inputSchema": {"type": "object", "properties": {}, "required": []},
        "func": tool_hub,
    },
    {
        "name": "check_verdicts",
        "description": (
            "Run the preregistered verdict table assertions "
            "(python/check_verdicts.py): 25 checks over 13 diagnosis "
            "artefacts; returns pass summary and any FAIL lines."),
        "inputSchema": {"type": "object", "properties": {}, "required": []},
        "func": tool_check_verdicts,
    },
    {
        "name": "legacy_analyze",
        "description": (
            "Analyze a legacy Fortran or C source file through the "
            "cross-language diagnosis chain (parse_legacy_source). Executes "
            "the source (via fingerprint_source/fingerprint_c_source) to "
            "obtain candidate_source=executed_source, or falls back to a "
            "synthetic equivalent on TranslationError. Returns a pure-dict "
            "report with file, language, extracted{system,precision,"
            "integrator,dt}, variant, candidate_source, assumptions, "
            "execution, and diagnosis."),
        "inputSchema": {
            "type": "object",
            "properties": {
                "file": {
                    "type": "string",
                    "description": "path to the legacy .f or .c source file"},
            },
            "required": ["file"],
        },
        "func": tool_legacy_analyze,
    },
]


def _handle(req: dict) -> dict | None:
    method = req.get("method")
    rid = req.get("id")
    params = req.get("params") or {}
    if rid is None:
        return None
    if method == "initialize":
        return {"jsonrpc": "2.0", "id": rid, "result": {
            "protocolVersion": params.get("protocolVersion",
                                          PROTOCOL_VERSION),
            "capabilities": {"tools": {"listChanged": False}},
            "serverInfo": {"name": "bench-mcp",
                           "version": "0.1.0",
                           "description": (
                               "Cross-language benchmark registry: run, "
                               "compare and diagnose Octave/Python/Fortran "
                               "implementations of attractors and "
                               "potentials.")},
        }}
    if method == "ping":
        return {"jsonrpc": "2.0", "id": rid, "result": {}}
    if method == "tools/list":
        tools = [{k: t[k] for k in ("name", "description", "inputSchema")}
                 for t in TOOLS]
        return {"jsonrpc": "2.0", "id": rid, "result": {"tools": tools}}
    if method == "tools/call":
        name = params.get("name")
        fn = None
        for t in TOOLS:
            if t["name"] == name:
                fn = t["func"]
                break
        if fn is None:
            return {"jsonrpc": "2.0", "id": rid, "error": {
                "code": -32602, "message": f"unknown tool: {name!r}"}}
        try:
            result = fn(params.get("arguments") or {})
            payload = json.dumps(result, indent=2, default=str)
            return {"jsonrpc": "2.0", "id": rid, "result": {
                "content": [{"type": "text", "text": payload}],
                "isError": False}}
        except Exception as exc:
            return {"jsonrpc": "2.0", "id": rid, "result": {
                "content": [{"type": "text",
                             "text": f"{type(exc).__name__}: {exc}"}],
                "isError": True}}
    return {"jsonrpc": "2.0", "id": rid, "error": {
        "code": -32601, "message": f"method not found: {method!r}"}}


def main() -> None:
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
        except json.JSONDecodeError:
            continue
        resp = _handle(req)
        if resp is not None:
            sys.stdout.write(json.dumps(resp) + "\n")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
