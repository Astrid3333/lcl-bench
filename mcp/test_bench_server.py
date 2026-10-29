#!/usr/bin/env python3
"""Smoke-test bench_server.py as a stdio JSON-RPC subprocess.

Launches the MCP server once and asserts: initialize handshake, tools/list,
registry_list, registry_get (ok + missing), diagnose (one variant), and an
unknown-tool error. Run: python3 mcp/test_bench_server.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SERVER = ROOT / "mcp" / "bench_server.py"
CHECKS = 0
FAILS = []


def check(name, got, want):
    global CHECKS
    CHECKS += 1
    if got != want:
        FAILS.append(f"{name}: got {got!r}, expected {want!r}")


class Client:
    def __init__(self):
        self.proc = subprocess.Popen(
            [sys.executable, str(SERVER)], cwd=ROOT,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        self.next_id = 0

    def request(self, method, params=None):
        self.next_id += 1
        msg = {"jsonrpc": "2.0", "id": self.next_id, "method": method}
        if params is not None:
            msg["params"] = params
        self.proc.stdin.write(json.dumps(msg) + "\n")
        self.proc.stdin.flush()
        line = self.proc.stdout.readline()
        return json.loads(line)

    def close(self):
        self.proc.stdin.close()
        self.proc.wait(timeout=10)


def main():
    c = Client()
    try:
        r = c.request("initialize", {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {"name": "test", "version": "0"},
        })
        check("initialize rc", r.get("error"), None)
        check("serverInfo.name", r["result"]["serverInfo"]["name"], "bench-mcp")
        check("capabilities.tools", "tools" in r["result"]["capabilities"], True)

        r = c.request("tools/list")
        names = sorted(t["name"] for t in r["result"]["tools"])
        check("tools/list", names, [
            "check_verdicts", "diagnose", "hub", "registry_get",
            "registry_list", "run_entry"])
        check("tool schema shape",
              all("inputSchema" in t for t in r["result"]["tools"]), True)

        r = c.request("tools/call",
                      {"name": "registry_list", "arguments": {}})
        check("registry_list isError", r["result"]["isError"], False)
        payload = json.loads(r["result"]["content"][0]["text"])
        check("registry_list n_entries", payload["n_entries"], 24)
        check("registry_list names",
              [e["name"] for e in payload["entries"]][:3],
              ["logistic", "henon", "clifford"])

        r = c.request("tools/call",
                      {"name": "registry_get", "arguments": {"name": "lorenz"}})
        check("registry_get isError", r["result"]["isError"], False)
        payload = json.loads(r["result"]["content"][0]["text"])
        check("registry_get family", payload["family"], "ode")
        check("registry_get dim", payload["dim"], 3)

        r = c.request("tools/call",
                      {"name": "registry_get", "arguments": {"name": "nope"}})
        check("registry_get missing isError", r["result"]["isError"], True)

        r = c.request("tools/call", {
            "name": "diagnose",
            "arguments": {"system": "lorenz", "variant": "reference"}})
        check("diagnose isError", r["result"]["isError"], False)
        payload = json.loads(r["result"]["content"][0]["text"])
        check("diagnose verdict", payload["diagnostico"], "equivalente")

        r = c.request("tools/call", {
            "name": "diagnose",
            "arguments": {"system": "lorenz", "variant": "not_a_variant"}})
        check("diagnose bad variant isError", r["result"]["isError"], True)

        r = c.request("tools/call", {"name": "unknown", "arguments": {}})
        check("unknown tool error code", r["error"]["code"], -32602)

        r = c.request("no/such/method")
        check("unknown method error code", r["error"]["code"], -32601)
    finally:
        c.close()

    print(f"{CHECKS - len(FAILS)}/{CHECKS} assertions passed")
    for f in FAILS:
        print("FAIL", f)
    return 1 if FAILS else 0


if __name__ == "__main__":
    raise SystemExit(main())
