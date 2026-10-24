"""Load JSON registries and evaluate entries with the same invariants as Octave.

The Python runner is a second independent implementation over the same
registry JSON. It compares computed invariants against the analytic or
literature references already stored in the registry (not against Octave),
so both languages independently hit the same ground truth.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from handles import map_step, rhs
from kernels import (box_counting_dim, corr_dim_grassberger, fd_eigenvalues,
                     lyap_map, lyap_ode)

ROOT = Path(__file__).resolve().parent.parent
REG_ATTRA = ROOT / "registries" / "attractors.json"
REG_POTE = ROOT / "registries" / "potentials.json"
REG_ENERGY = ROOT / "registries" / "energy_drift.json"
OUTDIR = ROOT / "results"


def load_registry(path: Path) -> list[dict]:
    with open(path) as f:
        data = json.load(f)
    return data["entries"]


def mkcheck(quantity, computed, reference, tol, rel_err, ok):
    return {
        "quantity": quantity,
        "computed": float(computed) if computed is not None else None,
        "reference": float(reference) if reference is not None and not (
            isinstance(reference, float) and np.isnan(reference)) else None,
        "tol": float(tol) if tol is not None else None,
        "rel_err": float(rel_err) if rel_err is not None else None,
        "pass": bool(ok),
    }


def eval_map(e: dict, outdir: Path) -> dict:
    f = map_step(e["name"], e["params"])
    x = np.asarray(e["x0"], dtype=float)
    norb = e["n_orbit"]
    orbit = np.zeros((norb, e["dim"]))
    for i in range(norb):
        x = f(x)
        orbit[i, :] = x
    if not np.all(np.isfinite(orbit)):
        raise RuntimeError(f"Orbit of {e['name']} is non-finite")
    odir = outdir / "orbits"
    odir.mkdir(parents=True, exist_ok=True)
    opath = odir / f"{e['name']}_orbit.txt"
    np.savetxt(opath, orbit, fmt="%.17g")

    lam = lyap_map(f, np.asarray(e["x0"], dtype=float),
                   e["n_lyap"], e["transient"])
    checks = []
    ref = e.get("ref") or {}
    # bit-exactness: first n iterates vs frozen reference sequence
    oh = ref.get("orbit_head")
    if oh is not None:
        ref_vals = np.asarray(oh["values"], dtype=float)
        head = orbit[: oh["n"]]
        md = float(np.max(np.abs(head - ref_vals)))
        tol = float(oh["abs_tol"])
        checks.append(mkcheck("orbit_head_maxabs", md, 0.0, tol, md,
                              md <= tol))
    lyap_ref = ref.get("lyap1")
    if lyap_ref is not None:
        val = lyap_ref["value"]
        tol = lyap_ref["abs_tol"]
        rel = abs(lam - val) / max(abs(val), np.finfo(float).eps)
        checks.append(mkcheck("lyap1", lam, val, tol, rel, abs(lam - val) <= tol))
    else:
        checks.append(mkcheck("lyap1_observed", lam, None, None, None, True))
    # Grassberger-Procaccia correlation dimension (protocol fixed in registry)
    cd_ref = ref.get("corr_dim")
    if cd_ref is not None:
        proto = cd_ref["protocol"]
        x = np.asarray(e["x0"], dtype=float)
        for _ in range(proto["transient"]):
            x = f(x)
        pts = np.zeros((proto["n_points"], e["dim"]))
        for i in range(proto["n_points"]):
            x = f(x)
            pts[i, :] = x
        dim = corr_dim_grassberger(pts, proto["n_radii"],
                                   proto["c_lo"], proto["c_hi"])
        val = cd_ref["value"]
        tol = cd_ref["rel_tol"] * abs(val)
        rel = abs(dim - val) / abs(val)
        checks.append(mkcheck("corr_dim", dim, val, tol, rel,
                              rel <= cd_ref["rel_tol"]))
    return {
        "name": e["name"], "family": e["family"],
        "pass": all(c["pass"] for c in checks),
        "checks": checks,
        "orbit_path": str(opath), "lyap1_observed": float(lam),
    }


def _chaos_game(e: dict, rng: np.random.Generator) -> np.ndarray:
    n = e["n_points"]
    bi = e["burn_in"]
    x = np.asarray(e["x0"], dtype=float)
    pts = np.zeros((n, 2))
    params = e["params"]
    if "vertices" in params:
        V = np.asarray(params["vertices"], dtype=float)
        p = np.asarray(params["p"], dtype=float)
        cp = np.cumsum(p / p.sum())
        for i in range(n + bi):
            k = np.searchsorted(cp, rng.random())
            x = (x + V[k]) / 2.0
            if i >= bi:
                pts[i - bi] = x
    elif "maps" in params:
        M = params["maps"]
        pv = np.array([m["p"] for m in M], dtype=float)
        cp = np.cumsum(pv / pv.sum())
        for i in range(n + bi):
            k = np.searchsorted(cp, rng.random())
            m = M[k]
            x = np.array([
                m["a"] * x[0] + m["b"] * x[1] + m["e"],
                m["c"] * x[0] + m["d"] * x[1] + m["f"],
            ])
            if i >= bi:
                pts[i - bi] = x
    else:
        raise ValueError(f"Unknown IFS param shape for {e['name']}")
    return pts


def eval_ifs(e: dict) -> dict:
    # fixed seed for run-to-run reproducibility in the cross-language hub
    rng = np.random.default_rng(0)
    pts = _chaos_game(e, rng)
    dim = box_counting_dim(pts, np.asarray(e["box_sizes"], dtype=float))
    checks = []
    ref = e.get("ref") or {}
    dim_ref = ref.get("dim")
    if dim_ref is not None:
        val = dim_ref["value"]
        tol = dim_ref["rel_tol"] * abs(val)
        rel = abs(dim - val) / abs(val)
        checks.append(mkcheck("box_dim", dim, val, tol, rel, rel <= dim_ref["rel_tol"]))
    else:
        checks.append(mkcheck("box_dim_observed", dim, None, None, None, True))
    return {
        "name": e["name"], "family": e["family"],
        "pass": all(c["pass"] for c in checks),
        "checks": checks, "dim_observed": float(dim),
    }


def eval_ode(e: dict) -> dict:
    f = rhs(e["name"], e["params"])
    lam = lyap_ode(f, np.asarray(e["x0"], dtype=float),
                   e["dt_lyap"], e["t_transient"], e["t_lyap"])
    checks = []
    ref = e.get("ref") or {}
    lyap_ref = ref.get("lyap1")
    if lyap_ref is not None:
        val = lyap_ref["value"]
        if "rel_tol" in lyap_ref:
            tol = lyap_ref["rel_tol"] * abs(val)
            rel = abs(lam - val) / abs(val)
            ok = rel <= lyap_ref["rel_tol"]
        else:
            tol = lyap_ref["abs_tol"]
            rel = abs(lam - val) / max(abs(val), np.finfo(float).eps)
            ok = abs(lam - val) <= tol
        checks.append(mkcheck("lyap1", lam, val, tol, rel, ok))
    else:
        checks.append(mkcheck("lyap1_observed", lam, None, None, None, True))
    return {
        "name": e["name"], "family": e["family"],
        "pass": all(c["pass"] for c in checks),
        "checks": checks, "lyap1_observed": float(lam),
    }


def eval_neural(e: dict) -> dict:
    P = np.asarray(e["params"]["patterns"], dtype=float)
    n = e["params"]["n"]
    W = (P.T @ P) / n
    W = W - np.diag(np.diag(W))
    checks = []
    for i, prow in enumerate(P):
        p = prow.reshape(-1, 1)
        s = np.sign(W @ p).ravel()
        ok = np.all(s == prow)
        checks.append(mkcheck(f"fixed_point_pattern_{i+1}",
                              float(ok), 1.0, 0.0, float(not ok), bool(ok)))
    return {
        "name": e["name"], "family": e["family"],
        "pass": all(c["pass"] for c in checks), "checks": checks,
    }


def eval_potential(e: dict) -> dict:
    name = e["name"].lower()
    ref = e.get("ref") or {}
    checks = []
    # scalar minimization: use Brent with tight xatol, same spirit as fminbnd
    def _minimize(f, lo, hi):
        res = minimize_scalar(f, bounds=(lo, hi), method="bounded",
                              options={"xatol": 1e-12})
        return float(res.x), float(res.fun)

    if name in ("lennard_jones", "morse", "double_well"):
        if name == "lennard_jones":
            ep, sg = e["params"]["eps"], e["params"]["sigma"]
            V = lambda r: 4.0 * ep * ((sg / r) ** 12 - (sg / r) ** 6)
            lo, hi = e["r_range"]
        elif name == "morse":
            De, a, re_ = e["params"]["De"], e["params"]["a"], e["params"]["re"]
            V = lambda r: De * (1.0 - np.exp(-a * (r - re_))) ** 2 - De
            lo, hi = e["r_range"]
        else:
            a, b = e["params"]["a"], e["params"]["b"]
            V = lambda x: a * x ** 4 - b * x ** 2
            lo, hi = 0.0, e["r_range"][1]
        xmin, vmin = _minimize(V, lo, hi)
        rmin = ref.get("r_min")
        if rmin is not None:
            val, tol = rmin["value"], rmin["abs_tol"]
            rel = abs(xmin - val) / abs(val)
            checks.append(mkcheck("r_min", xmin, val, tol, rel, abs(xmin - val) <= tol))
        rmin_abs = ref.get("r_min_abs")
        if rmin_abs is not None:
            val, tol = rmin_abs["value"], rmin_abs["abs_tol"]
            rel = abs(abs(xmin) - val) / abs(val)
            checks.append(mkcheck("r_min_abs", abs(xmin), val, tol, rel,
                                  abs(abs(xmin) - val) <= tol))
        vmin_ref = ref.get("V_min")
        if vmin_ref is not None:
            val, tol = vmin_ref["value"], vmin_ref["abs_tol"]
            rel = abs(vmin - val) / max(abs(val), np.finfo(float).eps)
            checks.append(mkcheck("V_min", vmin, val, tol, rel, abs(vmin - val) <= tol))

    elif name in ("harmonic", "hydrogen", "poschl_teller"):
        if name == "harmonic":
            g = np.linspace(e["x_grid"]["xmin"], e["x_grid"]["xmax"], e["x_grid"]["n"])
            Vg = 0.5 * e["params"]["k"] * g ** 2
            sigma = 0.0
        elif name == "hydrogen":
            n = e["x_grid"]["n"]
            rmax = e["x_grid"]["xmax"]
            h = rmax / (n + 1)
            g = h * np.arange(1, n + 1, dtype=float)
            Vg = -1.0 / g
            sigma = -1.0
        else:  # poschl_teller: full-line smooth well, no cusp
            g = np.linspace(e["x_grid"]["xmin"], e["x_grid"]["xmax"], e["x_grid"]["n"])
            lam = e["params"]["lambda"]
            Vg = -0.5 * lam * (lam + 1.0) / np.cosh(g) ** 2
            # bound states E_n = -(lambda-n)^2/2 < 0; continuum starts ~0.
            # Shift by -depth-ish so eigsh gets lowest bound states first.
            sigma = -0.5 * lam * (lam + 1.0) - 1.0
        k = len(e["eigen_n"])
        ev = fd_eigenvalues(g, Vg, k, sigma)
        eig_ref = ref["eigenvalues"]
        for j, en in enumerate(e["eigen_n"]):
            val = eig_ref["value"][j]
            comp = ev[j]
            rel = abs(comp - val) / abs(val)
            tol = eig_ref["rel_tol"] * abs(val)
            checks.append(mkcheck(f"E_{en}", comp, val, tol, rel,
                                  rel <= eig_ref["rel_tol"]))
        if "omega" in ref and ref["omega"] is not None:
            w = np.sqrt(e["params"]["k"] / e["params"]["m"])
            val, tol = ref["omega"]["value"], ref["omega"]["abs_tol"]
            checks.append(mkcheck("omega", w, val, tol, abs(w - val) / abs(val),
                                  abs(w - val) <= tol))

    elif name == "yukawa":
        A, kappa = e["params"]["A"], e["params"]["kappa"]
        Vh = lambda r: A * np.exp(-kappa * r) / r
        Fh = lambda r: A * np.exp(-kappa * r) * (1.0 / r ** 2 + kappa / r)
        hh = 1e-5
        for r in e["r_check"]:
            Fnum = -(Vh(r + hh) - Vh(r - hh)) / (2.0 * hh)
            Fan = Fh(r)
            val = Fan
            tol = ref["force"]["abs_tol"]
            rel = abs(Fnum - val) / abs(val)
            checks.append(mkcheck(f"force_r_{r:g}", Fnum, val, tol, rel,
                                  abs(Fnum - val) <= tol))
        F1 = -(Vh(1.0 + hh) - Vh(1.0 - hh)) / (2.0 * hh)
        val = ref["force"]["value_at_r1"]
        tol = ref["force"]["abs_tol"]
        checks.append(mkcheck("force_r_1_ref", F1, val, tol,
                              abs(F1 - val) / abs(val), abs(F1 - val) <= tol))
    else:
        raise ValueError(f"Unknown potential: {e['name']}")

    if not checks:
        raise RuntimeError(f"No checks produced for {e['name']}")
    return {
        "name": e["name"], "family": e["family"],
        "pass": all(c["pass"] for c in checks), "checks": checks,
    }


def eval_energy_drift(e: dict) -> dict:
    """Double-well trajectory: velocity Verlet vs RK4, relative energy drift."""
    a, b = e["params"]["a"], e["params"]["b"]
    mass = float(e["mass"])
    h = float(e["h"])
    T = float(e["T"])
    n = int(round(T / h))
    x = float(e["x0"][0])
    v = float(e["v0"][0])

    def V(q):
        return a * q ** 4 - b * q ** 2

    def F(q):
        return -4.0 * a * q ** 3 + 2.0 * b * q

    E0 = 0.5 * mass * v * v + V(x)
    # velocity Verlet
    xv, vv = x, v
    max_rel = 0.0
    for _ in range(n):
        acc = F(xv) / mass
        xv = xv + vv * h + 0.5 * acc * h * h
        acc_new = F(xv) / mass
        vv = vv + 0.5 * (acc + acc_new) * h
        E = 0.5 * mass * vv * vv + V(xv)
        rel = abs(E - E0) / abs(E0)
        if rel > max_rel:
            max_rel = rel
    verlet_final = rel
    verlet_max = max_rel
    # RK4 on (x, v)
    def f_y(y):
        return np.array([y[1], F(y[0]) / mass])
    y = np.array([x, v])
    max_rel = 0.0
    for _ in range(n):
        k1 = f_y(y)
        k2 = f_y(y + h / 2.0 * k1)
        k3 = f_y(y + h / 2.0 * k2)
        k4 = f_y(y + h * k3)
        y = y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        E = 0.5 * mass * y[1] ** 2 + V(y[0])
        rel = abs(E - E0) / abs(E0)
        if rel > max_rel:
            max_rel = rel
    rk4_final = rel
    rk4_max = max_rel

    ref = e.get("ref") or {}
    checks = []
    for key, computed in [("E0", E0),
                          ("verlet_final_rel_drift", verlet_final),
                          ("verlet_max_rel_drift", verlet_max),
                          ("rk4_final_rel_drift", rk4_final),
                          ("rk4_max_rel_drift", rk4_max)]:
        r = ref.get(key)
        if r is None:
            checks.append(mkcheck(f"{key}_observed", computed, None, None,
                                  None, True))
            continue
        val = r["value"]
        tol = r.get("abs_tol", r.get("rel_tol", 0.0) * abs(val))
        rel = abs(computed - val) / max(abs(val), np.finfo(float).eps)
        checks.append(mkcheck(key, computed, val, tol, rel,
                              abs(computed - val) <= tol))
    return {
        "name": e["name"], "family": e["family"],
        "pass": all(c["pass"] for c in checks), "checks": checks,
    }


def selected(e: dict, sel: str) -> bool:
    if sel == "all":
        return True
    if sel == "attractors":
        return e["family"] != "potential"
    if sel == "potentials":
        return e["family"] == "potential"
    return sel == e["name"] or sel == e["family"]


EVALUATORS = {
    "map": eval_map,
    "ifs": eval_ifs,
    "ode": eval_ode,
    "neural": eval_neural,
    "potential": eval_potential,
    "energy_drift": eval_energy_drift,
}


def printrun(r: dict) -> None:
    if "error" in r:
        print(f"  ERROR: {r['error']}")
        print("  => FAIL")
        return
    for c in r["checks"]:
        tag = "PASS" if c["pass"] else "FAIL"
        if c["reference"] is None:
            print(f"  {tag} {c['quantity']} (observed) = {c['computed']:.12g}")
        else:
            print(f"  {tag} {c['quantity']}: computed={c['computed']:.12g} "
                  f"ref={c['reference']:.12g} rel_err={c['rel_err']:.3g} "
                  f"tol={c['tol']:.3g}")
    print("  => PASS" if r["pass"] else "  => FAIL")


def run(selector: str = "all") -> dict:
    att = load_registry(REG_ATTRA)
    pot = load_registry(REG_POTE)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    runs = []
    for e in att:
        if not selected(e, selector):
            continue
        print(f"[{e['family']}] {e['name']}")
        try:
            r = EVALUATORS[e["family"]](e, OUTDIR) if e["family"] == "map" \
                else EVALUATORS[e["family"]](e)
        except Exception as exc:  # noqa: BLE001 — report per-entry, continue suite
            r = {"name": e["name"], "family": e["family"], "pass": False,
                 "checks": [mkcheck("error", None, None, None, None, False)],
                 "error": str(exc)}
        printrun(r)
        runs.append(r)
    for e in pot:
        if not selected(e, selector):
            continue
        print(f"[{e['family']}] {e['name']}")
        try:
            r = EVALUATORS["potential"](e)
        except Exception as exc:  # noqa: BLE001
            r = {"name": e["name"], "family": e["family"], "pass": False,
                 "checks": [mkcheck("error", None, None, None, None, False)],
                 "error": str(exc)}
        printrun(r)
        runs.append(r)
    for e in load_registry(REG_ENERGY):
        if not selected(e, selector):
            continue
        print(f"[{e['family']}] {e['name']}")
        try:
            r = EVALUATORS["energy_drift"](e)
        except Exception as exc:  # noqa: BLE001
            r = {"name": e["name"], "family": e["family"], "pass": False,
                 "checks": [mkcheck("error", None, None, None, None, False)],
                 "error": str(exc)}
        printrun(r)
        runs.append(r)

    npass = sum(1 for r in runs if r["pass"])
    print(f"\nSummary: {npass}/{len(runs)} entries PASS")

    out = {
        "timestamp": time.strftime("%Y%m%dT%H%M%S"),
        "selector": selector,
        "language": "python",
        "n_pass": npass, "n_total": len(runs), "runs": runs,
    }
    jpath = OUTDIR / f"registry_run_py_{out['timestamp']}.json"
    with open(jpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved: {jpath}")
    return out


if __name__ == "__main__":
    sel = sys.argv[1] if len(sys.argv) > 1 else "all"
    run(sel)
