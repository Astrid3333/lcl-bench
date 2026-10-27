"""diagnose_attractor: geometric fingerprint verdict for one implementation.

Runs a candidate trajectory (a synthetic degraded variant generated here, or
an external time series via --series), extracts the fingerprint (lambda1 via
Benettin, correlation dimension D2 via Grassberger-Procaccia), compares it
against the registry oracle (handles.py + kernels.py, registry protocols),
and returns a structured verdict:

    {"diagnostico": ..., "potencial": ..., "confianza": ...}

Taxonomy (from the MCP design notes + cross-language measurements in this
repo, notes/2026-10-24.md):
  equivalente                within oracle tolerance (no bug)
  colapso_por_redondeo       bounded orbit collapsed onto a short cycle
                             (tail unique_frac < 0.5x oracle); potential:
                             migrar_a_double_o_aumentar_precision
  precision_insuficiente      orbit healthy but lambda1=nan: the 1e-8
                             perturbation is annihilated by single-precision
                             casting; potential: migrar_a_double
  artefacto_tipo_sna         fractal-looking (D2 ok) but lambda1 <= 0 while
                             the oracle is chaotic; potential: revisar_integrador
  degradacion_de_orden       systematic lambda1/D2 shift beyond tolerance from
                             a lower-order integrator; potential: integrador_de_mayor_orden
  diferencia_por_integrador  small bounded method difference (within 2x tol),
                             a documented choice, not a bug
  degradacion_de_precision   float32 candidate stays chaotic but drifts beyond
                             tolerance (extension of the taxonomy)
  inestabilidad_numerica     orbit left the finite box; potential: reducir_dt_o_integrador_estable
  diferencia_no_explicada    discrepancy without candidate metadata (series mode)

Thresholds calibrated on this repo's measurements (NOT guessed):
  LAM_NOISE   0.02   abs floor; measured ulp band on lambda1 was ~+-0.015
                     (aizawa 1e-16 sweep, duffing 1-ulp sweep)
  LAM_REL     0.05   registry ODE lyap1 rel_tol
  LAM_REL_SERIES 0.10 series-mode lambda1 rel_tol: median-of-pairs
                     Rosenstein on realistic-size references errs <= 3.3%
                     (henon 0.3%, logistic 1.6%, lorenz 3.3%, lorenz scalar
                     embedding 3.0%); 10% = 3x margin. Paired with the
                     density gate SERIES_DENSITY_GATE (nn_med/D <= 5e-3,
                     separates accepted lorenz 0.0036 from rejected
                     rossler 0.0169 / aizawa 0.0135).
  D2_REL      0.15   estimator noise, measured: 12 consecutive oracle
                     windows (4000 pts, tail 2000) give D2 = 1.663..1.828,
                     std 0.058 = 3.3% rel; pairwise |dD2|/mean: median 3.7%,
                     max 9.4%, so 8% (registry henon rel_tol) false-mismatched
                     9% of SAME-trajectory window pairs. 15% ~= 4.6 sigma of
                     that noise, still far below genuine degradations
                     (euler_single measured 36% off).
  CYCLE_FRAC  0.10   unique-point fraction of the tail below this = cycle

Usage:
  python3 diagnose_attractor.py --system lorenz --variant euler
  python3 diagnose_attractor.py --system logistic --variant float32
  python3 diagnose_attractor.py --system henon --series path/orbit.txt
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from handles import map_step, rhs
from kernels import lyap_map, lyap_ode, corr_dim_grassberger, delay_embed, lyap1_series

ROOT = Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "registries" / "attractors.json"

LAM_NOISE = 0.02
LAM_REL = 0.05
LAM_REL_SERIES = 0.10
D2_REL = 0.15
CYCLE_FRAC = 0.10
N_POINTS_D2 = 4000
N_RADII = 30
C_LO, C_HI = 0.01, 0.5

SERIES_WINDOW = 20
SERIES_SAT = 0.2
SERIES_T_MIN = 0
SERIES_MIN_POINTS = 6
SERIES_MIN_SPAN = 8
SERIES_MIN_PAIRS = 100
SERIES_MAX_ANCHORS = 2000
SERIES_T_MAX = 5000
SERIES_DENSITY_GATE = 5e-3
SERIES_EMBED_M = 3
SERIES_METHOD = "rosenstein_series"

VARIANTS = ("reference", "float32", "euler", "euler_large", "euler_single")


def load_entry(name):
    entries = json.loads(REGISTRY.read_text())["entries"]
    for e in entries:
        if e["name"] == name:
            return e
    raise KeyError(f"system {name!r} not in {REGISTRY}")


def _f_of(entry):
    if entry["family"] == "map":
        return map_step(entry["name"], entry["params"])
    return rhs(entry["name"], entry["params"])


def _wrap_float32(f):
    def g(x):
        return np.asarray(f(np.asarray(x, dtype=np.float32)), dtype=np.float32)
    return g


def _rk4_stepper(f, h):
    def step(z):
        k1 = f(z)
        k2 = f(z + h / 2.0 * k1)
        k3 = f(z + h / 2.0 * k2)
        k4 = f(z + h * k3)
        return z + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return step


def _euler_stepper(f, h):
    return lambda z: z + h * f(z)


def lyap_with_stepper(step, x0, h, n_transient, n_steps, d0=1e-8):
    """Two-trajectory Benettin with an arbitrary stepper (see kernels.lyap_ode)."""
    x = np.asarray(x0, dtype=float).copy()
    for _ in range(n_transient):
        x = step(x)
    if not np.all(np.isfinite(x)):
        return float("nan")
    ones = np.ones_like(x)
    y = x + d0 * ones / np.linalg.norm(ones)
    acc = 0.0
    for _ in range(n_steps):
        x = step(x)
        y = step(y)
        d = y - x
        nd = np.linalg.norm(d)
        if not np.isfinite(nd) or nd == 0.0:
            return float("nan")
        acc += np.log(nd / d0)
        y = x + d0 * (d / nd)
    return acc / (n_steps * h)


def run_orbit(step, x0, n, cast64=True):
    """Iterate `step` n times from x0; return (points, finite, unique_frac)."""
    x = np.asarray(x0, dtype=float).copy()
    pts = np.full((n, x.size), np.nan)
    finite = True
    for i in range(n):
        x = step(x)
        if cast64:
            x = np.asarray(x, dtype=float)
        if not np.all(np.isfinite(x)):
            finite = False
            break
        pts[i] = x
    tail = pts[max(0, n // 2):] if finite else pts[:i]
    if finite:
        uniq = np.unique(np.round(tail, 12), axis=0).shape[0]
        frac = uniq / float(tail.shape[0])
    else:
        frac = 0.0
    return pts, finite, frac


def fingerprint_d2(points, finite):
    if not finite:
        return None
    tail = points[max(0, points.shape[0] // 2):]
    tail = tail[np.all(np.isfinite(tail), axis=1)]
    if tail.shape[0] < 100:
        return None
    try:
        return corr_dim_grassberger(tail, n_radii=N_RADII, c_lo=C_LO, c_hi=C_HI)
    except RuntimeError:
        return None  # too few nonzero pairs => collapsed


def oracle_fingerprint(entry):
    f = _f_of(entry)
    if entry["family"] == "map":
        lam = lyap_map(f, entry["x0"], entry["n_lyap"], entry["transient"])

        def step(x):
            return f(x)
        x0 = np.asarray(entry["x0"], dtype=float)
        for _ in range(entry["transient"]):
            x0 = step(x0)
        pts, finite, frac = run_orbit(step, x0, N_POINTS_D2)
    else:
        h = entry["dt_lyap"]
        lam = lyap_ode(f, entry["x0"], h, entry["t_transient"], entry["t_lyap"])
        step = _rk4_stepper(f, h)
        x = np.asarray(entry["x0"], dtype=float).copy()
        for _ in range(int(round(entry["t_transient"] / h))):
            x = step(x)
        pts, finite, frac = run_orbit(step, x, N_POINTS_D2)
    return {"lambda1": float(lam), "corr_dim": fingerprint_d2(pts, finite),
            "unique_frac": float(frac), "finite": bool(finite)}


def candidate_fingerprint(entry, variant):
    f = _f_of(entry)
    if variant == "reference":
        if entry["family"] == "map":
            step = f
            lam = lyap_map(f, entry["x0"], entry["n_lyap"], entry["transient"])
            x0 = np.asarray(entry["x0"], dtype=float)
            for _ in range(entry["transient"]):
                x0 = step(x0)
        else:
            h = entry["dt_lyap"]
            step = _rk4_stepper(f, h)
            lam = lyap_ode(f, entry["x0"], h,
                           entry["t_transient"], entry["t_lyap"])
            x0 = np.asarray(entry["x0"], dtype=float).copy()
            for _ in range(int(round(entry["t_transient"] / h))):
                x0 = step(x0)
        pts, finite, frac = run_orbit(step, x0, N_POINTS_D2)
    elif variant == "float32":
        f32 = _wrap_float32(f)
        if entry["family"] == "map":
            lam = lyap_map(f32, entry["x0"], entry["n_lyap"], entry["transient"])
            x0 = np.asarray(entry["x0"], dtype=float)
            for _ in range(entry["transient"]):
                x0 = f32(x0)
            step32 = f32
        else:
            h = entry["dt_lyap"]
            step32 = _rk4_stepper(f32, h)
            x = np.asarray(entry["x0"], dtype=float).copy()
            for _ in range(int(round(entry["t_transient"] / h))):
                x = step32(x)
            x0 = x
            lam = lyap_ode(f32, entry["x0"], h,
                           entry["t_transient"], entry["t_lyap"])
        pts, finite, frac = run_orbit(step32, x0, N_POINTS_D2, cast64=False)
    elif variant in ("euler", "euler_large", "euler_single"):
        if entry["family"] != "ode":
            raise ValueError(f"variant {variant!r} requires an ODE system")
        h = entry["dt_lyap"] * (5.0 if variant == "euler_large" else 1.0)
        if variant == "euler_single":
            f_use = _wrap_float32(f)
            raw_step = _euler_stepper(f_use, h)
            step = lambda z: np.asarray(raw_step(z), dtype=np.float32)
        else:
            step = _euler_stepper(f, h)
        n_tr = int(round(entry["t_transient"] / h))
        n_ly = int(round(entry["t_lyap"] / h))
        lam = lyap_with_stepper(step, entry["x0"], h, n_tr, n_ly)
        x = np.asarray(entry["x0"], dtype=float).copy()
        for _ in range(n_tr):
            x = step(x)
        pts, finite, frac = run_orbit(step, x, N_POINTS_D2,
                                      cast64=(variant != "euler_single"))
    else:
        raise ValueError(f"unknown variant {variant!r}; choose from {VARIANTS}")
    return {"lambda1": float(lam), "corr_dim": fingerprint_d2(pts, finite),
            "unique_frac": float(frac), "finite": bool(finite)}


def series_fingerprint(path, entry, dt=None):
    pts = np.loadtxt(path, ndmin=2)
    finite = bool(np.all(np.isfinite(pts)))
    h = float(dt) if dt is not None else (
        1.0 if entry["family"] == "map" else float(entry["dt_lyap"]))
    lam, info = None, {"reason": "orbit left the finite box (non-finite values)"}
    if finite:
        lam_pts = pts
        if pts.shape[1] == 1:
            emb, _ = delay_embed(pts[:, 0], SERIES_EMBED_M)
            if emb.shape[0] >= 10:
                lam_pts = emb
        if lam_pts.shape[0] >= 10:
            lam, info = lyap1_series(
                lam_pts[lam_pts.shape[0] // 2:], h, window=SERIES_WINDOW,
                sat_frac=SERIES_SAT, t_min=SERIES_T_MIN,
                min_points=SERIES_MIN_POINTS, min_span=SERIES_MIN_SPAN,
                min_pairs=SERIES_MIN_PAIRS, max_anchors=SERIES_MAX_ANCHORS,
                t_max=SERIES_T_MAX, density_gate=SERIES_DENSITY_GATE)
        else:
            info = {"reason": f"series too short ({pts.shape[0]} rows < 10)"}
    return {"lambda1": lam,
            "lambda1_method": SERIES_METHOD if lam is not None else None,
            "lambda1_reason": None if lam is not None else info.get("reason"),
            "series_dt": h,
            "series_mode": True,
            "corr_dim": fingerprint_d2(pts, finite),
            "unique_frac": None if not finite else
            float(np.unique(np.round(pts[pts.shape[0] // 2:], 12),
                            axis=0).shape[0] / (pts.shape[0] - pts.shape[0] // 2)),
            "finite": finite}


def classify(entry, oracle, cand, variant):
    notes = []
    lam_o, d2_o = oracle["lambda1"], oracle["corr_dim"]
    lam_c, d2_c = cand["lambda1"], cand["corr_dim"]
    series_mode = bool(cand.get("series_mode"))

    if not cand["finite"]:
        return ("inestabilidad_numerica", "reducir_dt_o_integrador_estable",
                0.95, ["orbit left the finite box (non-finite values)"])

    rel_lam = LAM_REL_SERIES if series_mode else LAM_REL
    tol_lam = max(LAM_NOISE, rel_lam * abs(lam_o))
    tol_d2 = D2_REL * abs(d2_o) if d2_o is not None else None

    # Orbit degradation: tail unique-count far below the oracle's
    # (measured: logistic float32 uf=0.2155 vs oracle 1.0 = genuine
    # roundoff collapse; lorenz euler_single uf=1.0 = healthy orbit).
    uf_o, uf_c = oracle["unique_frac"], cand["unique_frac"]
    cycle = uf_c is not None and (uf_c < CYCLE_FRAC or uf_c < 0.5 * uf_o)
    lam_nan = lam_c is not None and not np.isfinite(lam_c)

    if cycle:
        ev = [f"tail unique_frac={uf_c:.3f} vs oracle {uf_o:.3f} "
              "(orbit collapsed onto too few distinct values)"]
        if lam_nan:
            ev.append("lambda1=nan (perturbation direction annihilated)")
        return ("colapso_por_redondeo", "migrar_a_double_o_aumentar_precision",
                0.95, ev)

    if variant is not None and d2_c is None:
        return ("colapso_por_redondeo", "migrar_a_double_o_aumentar_precision",
                0.90, ["corr_dim failed: too few nonzero pairs "
                       "(states collapsed onto a short cycle)"])

    if lam_nan:
        # Orbit is fine; the 1e-8 Benettin perturbation is below the
        # single-precision eps at this state magnitude -> exponent
        # unmeasurable. Verified: float32(10.0)+float32(1e-8)==float32(10.0).
        ev = [f"orbit healthy (unique_frac={uf_c:.3f} ~ oracle {uf_o:.3f}) "
              "but lambda1=nan: the 1e-8 Lyapunov perturbation is annihilated "
              "by single-precision casting; the exponent is not measurable at "
              "this precision"]
        if variant and variant.startswith("euler"):
            ev.append("variant also uses first-order integration (euler); "
                      "its order shift is masked while lambda1 is unmeasurable")
        return ("precision_insuficiente", "migrar_a_double", 0.90, ev)

    if lam_c is None:  # series mode: lambda1 not measurable from this series
        why = cand.get("lambda1_reason") or "not available from a series alone"
        if tol_d2 is not None and d2_c is not None and \
                abs(d2_c - d2_o) <= tol_d2:
            return ("equivalente", None, 0.75,
                    [f"D2 {d2_c:.3f} within {D2_REL:.0%} of oracle {d2_o:.3f}; "
                     f"lambda1 {why}"])
        return ("diferencia_no_explicada", "verificar_parametros_o_ecuaciones",
                0.70, [f"corr_dim {d2_c} vs oracle {d2_o} beyond tolerance; "
                       f"lambda1 {why}"])

    dlam = abs(lam_c - lam_o)
    if dlam <= tol_lam and (tol_d2 is None or d2_c is None or
                            abs(d2_c - d2_o) <= tol_d2):
        m_lam = 1.0 - dlam / tol_lam
        m_d2 = 1.0 if tol_d2 is None or d2_c is None else \
            1.0 - abs(d2_c - d2_o) / tol_d2
        m = max(0.0, min(m_lam, m_d2))
        conf = round((0.75 + 0.10 * m) if series_mode else (0.85 + 0.10 * m), 2)
        ev = [f"dlambda1={dlam:.4g} <= {tol_lam:.4g}"]
        if series_mode:
            ev.insert(0, f"lambda1 via single-series estimator "
                         f"({cand.get('lambda1_method')}, "
                         f"dt={cand.get('series_dt')}, "
                         f"rel_tol={LAM_REL_SERIES:.0%})")
        if tol_d2 is not None and d2_c is not None:
            ev.append(f"dD2={abs(d2_c - d2_o):.4g} <= {tol_d2:.4g} "
                      f"(D2_REL={D2_REL:.0%} calibrated on estimator noise)")
        return ("equivalente", None, conf, ev)

    if lam_o > tol_lam and lam_c <= LAM_NOISE:
        if tol_d2 is None or d2_c is None or abs(d2_c - d2_o) <= 1.5 * tol_d2:
            return ("artefacto_tipo_sna", "revisar_integrador_o_precision",
                    0.80, [f"oracle lambda1={lam_o:.3f} chaotic but "
                           f"candidate lambda1={lam_c:.3f} ~ 0 while D2 intact"])
        return ("diferencia_no_explicada", "verificar_parametros_o_ecuaciones",
                0.75, [f"lambda1 {lam_c:.3f} vs {lam_o:.3f}, "
                       f"D2 {d2_c:.3f} vs {d2_o:.3f} both beyond tolerance"])

    ev = [f"dlambda1={dlam:.4g} > tol {tol_lam:.4g}"]
    if tol_d2 is not None and d2_c is not None:
        ev.append(f"dD2={abs(d2_c - d2_o):.4g} vs tol {tol_d2:.4g}")
    ratio = dlam / tol_lam
    if variant and variant.startswith("euler"):
        if dlam <= 2.0 * tol_lam:
            return ("diferencia_por_integrador", None, 0.75,
                    ev + ["bounded method difference within 2x tolerance; "
                          "documented choice, not a bug"])
        conf = min(0.95, 0.60 + 0.35 * min(1.0, ratio - 1.0))
        return ("degradacion_de_orden", "integrador_de_mayor_orden",
                round(conf, 2), ev + ["lower-order integrator shift beyond tolerance"])
    if variant == "float32":
        conf = min(0.95, 0.60 + 0.35 * min(1.0, ratio - 1.0))
        return ("degradacion_de_precision", "migrar_a_double",
                round(conf, 2), ev + ["float32 stays chaotic but drifts beyond tolerance"])
    notes.append("no candidate metadata (series mode): cannot attribute the shift")
    return ("diferencia_no_explicada", "verificar_parametros_o_ecuaciones",
            0.70, ev + notes)


def report_for(entry, oracle, cand, variant, candidate_source, series=None):
    """Build the verdict report for any candidate fingerprint (synthetic
    variant, executed legacy source, or time series)."""
    diagnostico, potencial, confianza, evidence = classify(
        entry, oracle, cand, variant)
    return {
        "system": entry["name"],
        "family": entry["family"],
        "variant": variant,
        "series": str(series) if series else None,
        "candidate_source": candidate_source,
        "fingerprint": {k: (None if isinstance(v, float) and not np.isfinite(v)
                            else v) for k, v in cand.items()},
        "oracle": {k: (None if isinstance(v, float) and not np.isfinite(v)
                       else v) for k, v in oracle.items()},
        "thresholds": {"lam_noise": LAM_NOISE, "lam_rel": LAM_REL,
                       "lam_rel_series": LAM_REL_SERIES,
                       "d2_rel": D2_REL, "cycle_frac": CYCLE_FRAC},
        "diagnostico": diagnostico,
        "potencial": potencial,
        "confianza": confianza,
        "evidence": evidence,
    }


def diagnose(system, variant="reference", series=None, series_dt=None):
    entry = load_entry(system)
    oracle = oracle_fingerprint(entry)
    if series is not None:
        cand = series_fingerprint(series, entry, dt=series_dt)
        return report_for(entry, oracle, cand, None, "series", series=series)
    cand = candidate_fingerprint(entry, variant)
    return report_for(entry, oracle, cand, variant, "synthetic_variant")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--system", required=True, help="registry name (e.g. lorenz)")
    ap.add_argument("--variant", default="reference", choices=VARIANTS,
                    help="synthetic candidate degradation (omit when --series)")
    ap.add_argument("--series", default=None,
                    help="external trajectory file (npy/txt, rows=states)")
    ap.add_argument("--series-dt", type=float, default=None,
                    help="sample spacing of --series (default: 1.0 for maps, "
                         "dt_lyap for ODEs)")
    ap.add_argument("--out", default=None, help="write report JSON here")
    args = ap.parse_args(argv)
    rep = diagnose(args.system, None if args.series else args.variant,
                   series=args.series, series_dt=args.series_dt)
    text = json.dumps(rep, indent=2)
    if args.out:
        Path(args.out).write_text(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
