#!/usr/bin/env python3
"""Branch tests for the single-series lambda1 estimator (Rosenstein-style
rosenstein_series) and its wiring in diagnose_attractor.

Every branch must be exercised (AGENTS.md principle 14):
kernels.acf_first_tau / delay_embed / lyap1_series (constant, too-few
pairs, density gate, too-few usable fits, non-finite slope, success),
series_fingerprint (1-col delay embedding with D2 kept raw, --series-dt,
finite gate, too-short series) and classify in series mode (equivalente
via measured lambda1, D2-only fallback after the density gate, 8% error
accepted under LAM_REL_SERIES but not LAM_REL, artefacto_tipo_sna,
diferencia_no_explicada, colapso, inestabilidad), plus the CLI.

Calibration this pins (see diagnose_attractor module docstring):
realistic-size references err <= 3.3% (henon 0.3%, logistic 1.8%,
lorenz 3.3%); LAM_REL_SERIES=0.10; D2_REL_SERIES=0.35 (series-D2
window max err 32.7% rossler; N=2000 <=23.6%); density gate
nn_med/D <= 5e-3 (lorenz_sparse tail d0D=0.0100 must decline).

Usage: python3 python/test_series_lambda.py   # exit 0 only if all pass
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

import diagnose_attractor as da                     # noqa: E402
from kernels import (                               # noqa: E402
    acf_first_tau, delay_embed, lyap1_series,
)

SERIES_DIR = ROOT / "fixtures" / "series"

CHECKS = 0
FAILS: list[str] = []


def check(label: str, cond: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    if not cond:
        FAILS.append(label)
        print(f"FAIL {label}" + (f" -- {detail}" if detail else ""))


_oracles: dict[str, dict] = {}


def oracle(name: str) -> dict:
    if name not in _oracles:
        _oracles[name] = da.oracle_fingerprint(da.load_entry(name))
    return _oracles[name]


def fp(system: str, fname: str, dt=None) -> dict:
    entry = da.load_entry(system)
    return da.series_fingerprint(str(SERIES_DIR / fname), entry, dt=dt)


def rel_err(lam: float, system: str) -> float:
    return abs(lam - oracle(system)["lambda1"]) / abs(oracle(system)["lambda1"])


def test_acf_and_embed() -> None:
    rng = np.random.default_rng(7)
    ar = np.zeros(4000)
    for k in range(1, 4000):
        ar[k] = 0.9 * ar[k - 1] + rng.standard_normal()
    tau_ar = acf_first_tau(ar)
    # rho(k) = 0.9^k crosses 1/e at k = 10.03
    check("acf: AR(1) tau at analytic 1/e crossing", 8 <= tau_ar <= 12,
          f"tau={tau_ar}")
    check("acf: constant series -> 1", acf_first_tau(np.full(500, 3.0)) == 1)
    noise = rng.standard_normal(4000)
    check("acf: white noise decorrelates at lag 1", acf_first_tau(noise) <= 2,
          f"tau={acf_first_tau(noise)}")
    check("acf: AR persists longer than noise", acf_first_tau(noise) < tau_ar)
    rows, tau = delay_embed(ar, m=3)
    check("embed: shape", rows.shape == (4000 - 2 * tau, 3),
          f"{rows.shape} tau={tau}")
    check("embed: first column is the series head",
          np.allclose(rows[:, 0], ar[:4000 - 2 * tau]))
    short, _ = delay_embed(np.arange(5.0), m=3)
    check("embed: too-short input -> empty", short.shape[0] == 0)


def test_lyap1_branches() -> None:
    hen = np.loadtxt(SERIES_DIR / "henon_reference.txt", ndmin=2)
    tail = hen[hen.shape[0] // 2:]
    lam, info = lyap1_series(tail, 1.0)
    check("lyap1: henon success", lam is not None and info.get("reason") is None)
    check("lyap1: henon within 10%", lam is not None and rel_err(lam, "henon") <= 0.10,
          f"lam={lam} err={None if lam is None else rel_err(lam, 'henon')}")
    check("lyap1: info carries d0D/npairs",
          "d0D" in info and "npairs" in info and info["npairs"] >= 100)

    lam, info = lyap1_series(np.full((400, 3), 0.7), 1.0)
    check("lyap1: constant -> None", lam is None)
    check("lyap1: constant reason", "constant" in info.get("reason", ""))

    lam, info = lyap1_series(hen[:50], 1.0)
    check("lyap1: too few pairs -> None", lam is None)
    check("lyap1: too few pairs reason",
          "too few neighbor pairs" in info.get("reason", ""),
          info.get("reason"))

    lam, info = lyap1_series(tail, 1.0, sat_frac=1e-9)
    check("lyap1: too few usable fits -> None", lam is None)
    check("lyap1: too few usable fits reason",
          "too few usable fits" in info.get("reason", ""),
          info.get("reason"))

    lam, info = lyap1_series(tail, 0.0)
    check("lyap1: h=0 -> None + spacing reason", lam is None and
          "sample spacing" in info.get("reason", ""), info.get("reason"))
    lam, info = lyap1_series(tail, float("nan"))
    check("lyap1: h=nan -> non-finite -> None", lam is None)
    check("lyap1: non-finite reason",
          "non-finite" in info.get("reason", ""), info.get("reason"))

    sp = np.loadtxt(SERIES_DIR / "lorenz_sparse.txt", ndmin=2)
    lam, info = lyap1_series(sp[sp.shape[0] // 2:], 0.01)
    check("lyap1: sparse lorenz declines", lam is None)
    check("lyap1: density gate reason",
          "density gate" in info.get("reason", "") and info.get("d0D", 0) > 5e-3,
          f"info={info}")


def test_series_fingerprint() -> None:
    hen = fp("henon", "henon_reference.txt")
    check("fp henon: lambda measured", hen["lambda1"] is not None)
    check("fp henon: method/series_mode/dt",
          hen["lambda1_method"] == "rosenstein_series" and
          hen["series_mode"] is True and hen["series_dt"] == 1.0)
    check("fp henon: within 10%", rel_err(hen["lambda1"], "henon") <= 0.10)
    check("fp henon: reason None", hen["lambda1_reason"] is None)

    logi = fp("logistic", "logistic_reference.txt")
    check("fp logistic: 1-col embed measured", logi["lambda1"] is not None)
    check("fp logistic: within 10%", rel_err(logi["lambda1"], "logistic") <= 0.10,
          f"lam={logi['lambda1']} err={rel_err(logi['lambda1'], 'logistic')}")
    d2o = oracle("logistic")["corr_dim"]
    check("fp logistic: D2 compared on raw rows (not embedding)",
          abs(logi["corr_dim"] - d2o) <= 0.15 * abs(d2o),
          f"D2={logi['corr_dim']} vs oracle {d2o}")

    lor = fp("lorenz", "lorenz_reference.txt")
    check("fp lorenz: within 10%", rel_err(lor["lambda1"], "lorenz") <= 0.10,
          f"lam={lor['lambda1']} err={rel_err(lor['lambda1'], 'lorenz')}")
    check("fp lorenz: default dt = dt_lyap", lor["series_dt"] == 0.01)

    sparse = fp("lorenz", "lorenz_sparse.txt")
    check("fp sparse: lambda None", sparse["lambda1"] is None)
    check("fp sparse: density gate reason",
          "density gate" in (sparse["lambda1_reason"] or ""),
          sparse["lambda1_reason"])

    const = fp("henon", "const_series.txt")
    check("fp const: lambda None + constant reason",
          const["lambda1"] is None and
          "constant" in (const["lambda1_reason"] or ""))
    check("fp const: unique_frac tiny", const["unique_frac"] < 0.01)

    short = fp("henon", "short_series.txt")
    check("fp short: lambda None + pairs reason",
          short["lambda1"] is None and
          "too few neighbor pairs" in (short["lambda1_reason"] or ""))

    nan = fp("henon", "nan_series.txt")
    check("fp nan: finite False", nan["finite"] is False)
    check("fp nan: lambda None + finite-box reason",
          nan["lambda1"] is None and
          "finite box" in (nan["lambda1_reason"] or ""))
    check("fp nan: D2/unique_frac None",
          nan["corr_dim"] is None and nan["unique_frac"] is None)

    tiny = da.load_entry("henon")
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write("0.1 0.1\n0.2 0.3\n0.4 0.5\n0.6 0.1\n0.7 0.2\n0.8 0.1\n")
        tmp = fh.name
    t = da.series_fingerprint(tmp, tiny)
    check("fp too short (<10 rows): reason",
          t["lambda1"] is None and "too short" in (t["lambda1_reason"] or ""),
          t["lambda1_reason"])
    Path(tmp).unlink()

    d02 = fp("henon", "henon_reference.txt", dt=0.02)
    check("fp --series-dt override recorded", d02["series_dt"] == 0.02)
    check("fp lambda scales as 1/dt",
          abs(d02["lambda1"] / hen["lambda1"] - 50.0) < 1e-9,
          f"ratio={d02['lambda1'] / hen['lambda1']}")
    check("fp fingerprint is JSON-serializable",
          json.dumps({k: v for k, v in hen.items()}) is not None)


def test_classify_series() -> None:
    entry = da.load_entry("henon")
    o = oracle("henon")
    good = fp("henon", "henon_reference.txt")
    d = da.classify(entry, o, good, None)
    check("classify: henon equivalente", d[0] == "equivalente" and d[1] is None)
    check("classify: series conf in [0.75, 0.85]", 0.75 <= d[2] <= 0.85,
          f"conf={d[2]}")
    check("classify: evidence names the estimator",
          any("rosenstein_series" in e for e in d[3]), str(d[3]))

    logi_e = da.load_entry("logistic")
    d = da.classify(logi_e, oracle("logistic"),
                    fp("logistic", "logistic_reference.txt"), None)
    check("classify: logistic equivalente", d[0] == "equivalente", str(d))

    d = da.classify(da.load_entry("lorenz"), oracle("lorenz"),
                    fp("lorenz", "lorenz_reference.txt"), None)
    check("classify: lorenz equivalente", d[0] == "equivalente", str(d))

    sp = fp("lorenz", "lorenz_sparse.txt")
    d = da.classify(da.load_entry("lorenz"), oracle("lorenz"), sp, None)
    check("classify: sparse falls back to D2-only equivalente",
          d[0] == "equivalente" and d[2] == 0.75, str(d))
    check("classify: sparse evidence carries the gate reason",
          any("density gate" in e for e in d[3]), str(d[3]))

    d2ok = dict(good)
    d2ok["lambda1"] = None
    d2ok["lambda1_reason"] = "density gate"
    d2ok["corr_dim"] = o["corr_dim"] * 1.20
    d = da.classify(entry, o, d2ok, None)
    check("classify: D2-only 20% within D2_REL_SERIES -> equivalente",
          d[0] == "equivalente" and d[2] == 0.75, str(d))
    check("classify: D2-only evidence names the 35% tolerance",
          any("35%" in e for e in d[3]), str(d[3]))

    d2far = dict(good)
    d2far["lambda1"] = None
    d2far["lambda1_reason"] = "density gate"
    d2far["corr_dim"] = o["corr_dim"] * 1.40
    d = da.classify(entry, o, d2far, None)
    check("classify: D2-only 40% beyond D2_REL_SERIES -> diferencia_no_explicada",
          d[0] == "diferencia_no_explicada", str(d))

    d = da.classify(entry, o, fp("henon", "const_series.txt"), None)
    check("classify: const -> colapso_por_redondeo",
          d[0] == "colapso_por_redondeo" and
          d[1] == "migrar_a_double_o_aumentar_precision", str(d))

    d = da.classify(entry, o, fp("henon", "nan_series.txt"), None)
    check("classify: nan -> inestabilidad_numerica",
          d[0] == "inestabilidad_numerica", str(d))

    d = da.classify(entry, o, fp("henon", "short_series.txt"), None)
    check("classify: short -> diferencia_no_explicada",
          d[0] == "diferencia_no_explicada", str(d))

    # 8% series error passes only because of LAM_REL_SERIES (legacy = 5%)
    lam8 = dict(good)
    lam8["lambda1"] = o["lambda1"] * 1.08
    d = da.classify(entry, o, lam8, None)
    check("classify: 8% series error accepted by LAM_REL_SERIES",
          d[0] == "equivalente", str(d))

    far = dict(good)
    far["lambda1"] = o["lambda1"] * 1.40
    d = da.classify(entry, o, far, None)
    check("classify: 40% series error -> diferencia_no_explicada",
          d[0] == "diferencia_no_explicada", str(d))

    dead = dict(good)
    dead["lambda1"] = 0.01
    d = da.classify(entry, o, dead, None)
    check("classify: series lambda ~0 with intact D2 -> artefacto_tipo_sna",
          d[0] == "artefacto_tipo_sna", str(d))

    rep = da.report_for(entry, o, good, None, "series", series="x.txt")
    check("report: thresholds expose lam_rel_series",
          rep["thresholds"].get("lam_rel_series") == 0.10)
    check("report: thresholds expose d2_rel_series",
          rep["thresholds"].get("d2_rel_series") == 0.35)
    check("report: JSON round-trip keeps verdict",
          json.loads(json.dumps(rep))["diagnostico"] == "equivalente")


def test_cli() -> None:
    import subprocess
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        out = fh.name
    r = subprocess.run(
        [sys.executable, str(HERE / "diagnose_attractor.py"),
         "--system", "henon", "--series", str(SERIES_DIR / "henon_reference.txt"),
         "--out", out],
        capture_output=True, text=True, cwd=str(ROOT))
    check("cli: exit 0", r.returncode == 0, r.stderr[-400:])
    rep = json.loads(Path(out).read_text())
    check("cli: henon series verdict", rep["diagnostico"] == "equivalente",
          rep["diagnostico"])
    check("cli: lambda1 present in report",
          rep["fingerprint"]["lambda1"] is not None)
    Path(out).unlink()

    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        out = fh.name
    r = subprocess.run(
        [sys.executable, str(HERE / "diagnose_attractor.py"),
         "--system", "henon", "--series", str(SERIES_DIR / "henon_reference.txt"),
         "--series-dt", "0.02", "--out", out],
        capture_output=True, text=True, cwd=str(ROOT))
    check("cli: --series-dt exit 0", r.returncode == 0, r.stderr[-400:])
    rep = json.loads(Path(out).read_text())
    check("cli: --series-dt reaches fingerprint",
          rep["fingerprint"]["series_dt"] == 0.02)
    Path(out).unlink()


def main() -> int:
    t0 = time.time()
    test_acf_and_embed()
    test_lyap1_branches()
    test_series_fingerprint()
    test_classify_series()
    test_cli()
    dt = time.time() - t0
    if FAILS:
        print(f"\n{len(FAILS)}/{CHECKS} CHECKS FAILED ({dt:.1f}s)")
        return 1
    print(f"ALL {CHECKS} CHECKS PASS ({dt:.1f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
