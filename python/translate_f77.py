"""Strict-subset Fortran 77 -> Python translator: EXECUTE legacy code.

Translates an old fixed-form F77 attractor driver into a Python step
function and runs it under the precision declared in the source (REAL ->
np.float32, REAL*8/DOUBLE PRECISION -> np.float64, literals without a D
exponent -> single, with -> double), so the fingerprint (lambda1, D2) comes
from the LEGACY recurrence itself instead of a synthetic equivalent variant.

Supported subset (anything else raises TranslationError and the caller falls
back to the synthetic-equivalent diagnosis):
  - fixed form: col-1 comments c/C/*, inline `!`, col-6 `&` continuation
  - program/end/implicit none, integer/real/real*4/real*8/double precision
    scalar declarations (no arrays), parameter(name = <numeric literal>)
  - scalar assignments: + - * / ** ( ), unary minus, declared names
  - ONE do var = a, b  ...  end do (loop counter must NOT appear in the body)
  - write(*,*) name, ...   (evaluated at the end, informational only)
The translated source is exposed as ns["__source"] for inspection.

HONESTY: lambda1 uses the legacy program's own loop count as transient and
the registry window afterwards; D2 uses the last 4000 recorded steps of the
legacy run (same per-step sampling density as the oracle). The recorded
variables are the first `entry.dim` state floats in declaration order.

Usage:
  from translate_f77 import fingerprint_source
  cand, info = fingerprint_source(text, entry, h=0.01)
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from diagnose_attractor import N_POINTS_D2, fingerprint_d2

_TOKEN = re.compile(
    r"""(?P<d>(?:[0-9]+\.[0-9]*|\.[0-9]+|[0-9]+)[dD][+-]?[0-9]+)
      | (?P<f>(?:[0-9]+\.[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?
             |[0-9]+[eE][+-]?[0-9]+)
      | (?P<i>[0-9]+)
      | (?P<name>[A-Za-z][A-Za-z0-9_]*)
      | (?P<pow>\*\*)
      | (?P<op>[-+*/()])""",
    re.X)


class TranslationError(Exception):
    def __init__(self, line, msg):
        self.line = line
        super().__init__(f"line {line}: {msg}")


def logical_lines(text):
    """Fixed-form pass: comments, continuations, inline '!'. -> [(ln, stmt)]"""
    out, cur, cur_ln = [], None, None
    for i, raw in enumerate(text.splitlines(), 1):
        if raw[:1] in ("", "c", "C", "*", "!"):
            continue
        piece = raw.split("!", 1)[0].rstrip()
        if len(raw) > 5 and raw[5] == "&":
            if cur is None:
                raise TranslationError(i, "continuation without a statement")
            cur += piece[6:]
            continue
        if not piece.strip():
            continue
        if cur is not None:
            out.append((cur_line, cur))
        cur, cur_line = piece, i
    if cur is not None:
        out.append((cur_line, cur))
    return out


def _promote(a, b):
    if a == "i" and b == "i":
        return "i"
    if a == "f" and b == "f":
        return "f"
    if a == "d" or b == "d":
        return "d"
    return "f"


class _ExprParser:
    def __init__(self, tokens, types, lineno):
        self.toks, self.i = tokens, 0
        self.types, self.ln = types, lineno
        self.used = set()

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else (None, None)

    def next(self):
        tok = self.peek()
        self.i += 1
        return tok

    def parse(self):
        v = self.expr(1)
        if self.i != len(self.toks):
            raise TranslationError(self.ln, f"junk after expression: "
                                    f"{self.toks[self.i][1]!r}")
        return v

    def expr(self, min_prec):
        left = self.unary()
        while True:
            kind, text = self.peek()
            if kind == "pow" and 3 >= min_prec:
                self.next()
                right = self.expr(3)  # right associative
                left = (f"({left[0]} ** {right[0]})",
                        _promote(left[1], right[1]))
                continue
            if kind != "op" or text not in "+-*/":
                break
            prec = 1 if text in "+-" else 2
            if prec < min_prec:
                break
            self.next()
            right = self.expr(prec + 1)
            left = self.binary(text, left, right)
        return left

    def binary(self, op, L, R):
        (py_a, ta), (py_b, tb) = L, R
        if op == "/" and ta == "i" and tb == "i":
            return (f"({py_a} // {py_b})", "i")
        return (f"({py_a} {op} {py_b})", _promote(ta, tb))

    def unary(self):
        kind, text = self.peek()
        if kind == "op" and text == "-":
            self.next()
            py, t = self.unary()
            return (f"(-{py})", t)
        if kind == "op" and text == "+":
            self.next()
            return self.unary()
        if kind == "op" and text == "(":
            self.next()
            v = self.expr(1)
            k2, t2 = self.next()
            if k2 != "op" or t2 != ")":
                raise TranslationError(self.ln, "missing ')'")
            return v
        if kind in ("d", "f", "i"):
            self.next()
            if kind == "d":
                val = float(text.replace("d", "e").replace("D", "e"))
                return (f"np.float64({val!r})", "d")
            if kind == "f":
                return (f"np.float32({float(text)!r})", "f")
            return (str(int(text)), "i")
        if kind == "name":
            self.next()
            nm = text.lower()
            if nm not in self.types:
                raise TranslationError(self.ln, f"undeclared variable {nm!r}")
            self.used.add(nm)
            return (nm, self.types[nm])
        raise TranslationError(self.ln,
                               f"unexpected token {text!r} in expression")


def _tokenize(expr, ln):
    toks, pos = [], 0
    for m in _TOKEN.finditer(expr):
        if expr[pos:m.start()].strip():
            raise TranslationError(ln, f"bad expression near "
                                    f"{expr[pos:m.start()]!r}")
        pos = m.end()
        kind = m.lastgroup
        toks.append((kind, m.group(kind)))
    if expr[pos:].strip():
        raise TranslationError(ln, f"bad expression near {expr[pos:]!r}")
    if not toks:
        raise TranslationError(ln, "empty expression")
    return toks


def _parse_expr(expr, types, ln):
    p = _ExprParser(_tokenize(expr, ln), types, ln)
    py, typ = p.parse()
    return py, typ, p.used


def _assign_code(lhs, rhs, types, ln):
    if lhs not in types:
        raise TranslationError(ln, f"undeclared variable {lhs!r}")
    py, typ, used = _parse_expr(rhs, types, ln)
    if types[lhs] == "f":
        return f"{lhs} = np.float32({py})", used
    if types[lhs] == "d":
        return f"{lhs} = np.float64({py})", used
    return f"{lhs} = int({py})", used


def translate(text, dim=None):
    """Translate strict-subset F77 source -> executable Python source.

    `dim` = expected state dimension (entry["dim"]); the recorded/observed
    variables are then the first `dim` state floats in declaration order.
    """
    lines = logical_lines(text)
    types, decl_order = {}, []
    init_codes, init_names, body_codes = [], [], []
    avail, body_avail = set(), set()
    loop = None          # (var, start_py, stop_py, start_used, stop_used)
    write_vars = None
    phase = "pre"
    ended = False

    def decl(ln, st):
        m = re.match(r"^(integer|real|double\s+precision)\b\s*(.*)$",
                     st, re.I)
        kw = re.sub(r"\s+", " ", m.group(1).lower())
        rest = m.group(2)
        if kw == "integer":
            typ, names_part = "i", rest
            if re.match(r"^\s*[* (]", names_part):
                raise TranslationError(ln, "kinds/arrays on integer unsupported")
        elif kw == "double precision":
            typ, names_part = "d", rest
        else:
            mm = re.match(
                r"^\s*(?:\*\s*(?P<s>[0-9]+)"
                r"|\(\s*(?:kind\s*=\s*)?(?P<k>[0-9]+)\s*\))\s*(?P<rest>.*)$",
                rest)
            if mm:
                code = mm.group("s") or mm.group("k")
                if code == "8":
                    typ = "d"
                elif code == "4":
                    typ = "f"
                else:
                    raise TranslationError(ln, f"unsupported real kind {code}")
                names_part = mm.group("rest")
            else:
                typ, names_part = "f", rest
        for nm in names_part.split(","):
            nm = nm.strip()
            if not re.match(r"^[A-Za-z][A-Za-z0-9_]*$", nm):
                raise TranslationError(ln, f"unsupported declaration item "
                                       f"{nm!r} (arrays unsupported)")
            nm = nm.lower()
            if nm in types and types[nm] != typ:
                raise TranslationError(ln, f"conflicting declaration of {nm!r}")
            if nm not in types:
                types[nm] = typ
                decl_order.append(nm)

    def check_used(ln, used, avail_set):
        missing = sorted(used - avail_set)
        if missing:
            raise TranslationError(ln, f"read of not-yet-assigned variable(s) "
                                   f"{missing}")

    for ln, raw_st in lines:
        st = raw_st.strip()
        low = st.lower()
        if ended:
            raise TranslationError(ln, "statement after program end")
        if re.match(r"^(program|module)\b", low) or re.match(r"^implicit\b", low):
            continue
        if phase == "post" and re.match(r"^end\s*(program)?\s*$", low):
            ended = True
            continue
        if re.match(r"^end\s*do\s*$", low) or re.match(r"^enddo\s*$", low):
            if phase != "body":
                raise TranslationError(ln, "end do without do")
            phase = "post"
            continue
        if re.match(r"^(integer|real|double\s+precision)\b", low, re.I):
            if phase != "pre":
                raise TranslationError(ln, "declaration after do-loop")
            decl(ln, st)
            continue
        pm = re.match(r"^parameter\s*\((.*)\)$", st, re.I)
        if pm:
            if phase != "pre":
                raise TranslationError(ln, "parameter after do-loop")
            for part in pm.group(1).split(","):
                nm, eq, lit = part.partition("=")
                nm = nm.strip().lower()
                if not eq:
                    raise TranslationError(ln, f"bad parameter {part!r}")
                if re.search(r"[A-Za-z]", lit):
                    raise TranslationError(ln, "parameter expressions with "
                                           "names unsupported (literal only)")
                code, used = _assign_code(nm, lit.strip(), types, ln)
                check_used(ln, used, avail)
                init_codes.append(code)
                if nm not in init_names:
                    init_names.append(nm)
                avail.add(nm)
            continue
        dm = re.match(r"^do\s+([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.+?),\s*(.+)$",
                      st, re.I)
        if dm:
            if phase != "pre":
                raise TranslationError(ln, "multiple/nested do-loops "
                                       "unsupported (single do-loop subset)")
            var = dm.group(1).lower()
            if types.get(var) != "i":
                raise TranslationError(ln, f"do variable {var!r} must be "
                                       "declared integer")
            sp, stp, su = _parse_expr(dm.group(2).strip(), types, ln), None, None
            su = sp[2]
            stop = _parse_expr(dm.group(3).strip(), types, ln)
            check_used(ln, su | stop[2], avail)
            if sp[1] != "i" or stop[1] != "i":
                raise TranslationError(ln, "do bounds must be integer")
            loop = (var, sp[0], stop[0])
            loop_avail = set(avail)
            phase = "body"
            continue
        wm = re.match(r"^write\s*\(\s*\*\s*,\s*\*\s*\)\s*(.*)$", st, re.I)
        if wm:
            if phase != "post":
                raise TranslationError(ln, "write only allowed after end do")
            write_vars = [v.strip().lower() for v in wm.group(1).split(",")]
            for v in write_vars:
                if v not in types:
                    raise TranslationError(ln, f"undeclared write target {v!r}")
            continue
        am = re.match(r"^([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.+)$", st)
        if am and not st.startswith("="):
            lhs, rhs = am.group(1).lower(), am.group(2)
            if rhs.lstrip().startswith("="):
                raise TranslationError(ln, "relational == unsupported")
            if phase == "pre":
                code, used = _assign_code(lhs, rhs, types, ln)
                check_used(ln, used, avail)
                init_codes.append(code)
                if lhs not in init_names:
                    init_names.append(lhs)
                avail.add(lhs)
            elif phase == "body":
                if lhs == loop[0]:
                    raise TranslationError(ln, "assignment to do-loop variable "
                                           "unsupported (counter excluded)")
                code, used = _assign_code(lhs, rhs, types, ln)
                check_used(ln, used, avail | body_avail)
                body_codes.append(code)
                body_avail.add(lhs)
            else:
                raise TranslationError(ln, "assignment after end do "
                                       "(only write/end allowed)")
            continue
        raise TranslationError(ln, f"unsupported statement: {st!r}")

    if loop is None:
        raise TranslationError(0, "no do-loop found (cannot build a step)")
    if phase != "post":
        raise TranslationError(0, "do-loop never closed with end do")
    if not body_codes:
        raise TranslationError(0, "empty do-loop body")

    state_names = [nm for nm in init_names
                   if types.get(nm) in ("f", "d") and nm in body_avail]
    # record ONLY the first `dim` state floats in declaration order: keeps the
    # state vector, drops clock/parameter/scratch floats (e.g. t, dt, sigma).
    float_vars = [nm for nm in decl_order if nm in state_names]
    if dim is not None:
        float_vars = float_vars[:int(dim)]

    src = ["import numpy as np", "",
           f"FLOAT_VARS = {float_vars!r}",
           f"STATE_NAMES = {init_names!r}",
           f"LOOP_VAR = {loop[0]!r}", "",
           "def legacy_init():",
           "    v = {}"]
    src += [f"    {c}" for c in init_codes]
    src += [f"    v[{nm!r}] = {nm}" for nm in init_names]
    src += ["    return v", "",
            "def legacy_step(v, i=0):"]
    src += [f"    {nm} = v[{nm!r}]" for nm in init_names]
    src += [f"    {c}" for c in body_codes]
    src += [f"    v[{nm!r}] = {nm}" for nm in init_names]
    src += ["    return v", "",
            "def legacy_steps(v=None):"]
    src += [f"    {nm} = v[{nm!r}]" for nm in init_names]
    src += [f"    _s = {loop[1]}", f"    _e = {loop[2]}",
            "    return max(0, int(_e) - int(_s) + 1)"]
    if write_vars:
        src += ["", "def legacy_final(v):"]
        src += [f"    {nm} = v[{nm!r}]" for nm in write_vars]
        src += [f"    return [{', '.join(write_vars)}]"]
    src += [""]
    return "\n".join(src)


def load_legacy(text, dim=None):
    src = translate(text, dim)
    ns = {}
    exec(compile(src, "<legacy_f77>", "exec"), ns)
    ns["__source"] = src
    return ns


def lyap_legacy(step, v0, float_vars, n_tr, n_steps, h, d0=1e-8):
    """Benettin on the EXECUTED legacy step (two copies of the state dict)."""
    x = dict(v0)
    with np.errstate(all="ignore"):
        for _ in range(n_tr):
            try:
                x = step(x, 0)
            except Exception:
                return float("nan")
            if not all(math.isfinite(float(x[k])) for k in float_vars):
                return float("nan")
        y = dict(x)
        for k in float_vars:
            y[k] = x[k] + d0 / math.sqrt(len(float_vars))
        acc = 0.0
        for _ in range(n_steps):
            try:
                x = step(x, 0)
                y = step(y, 0)
            except Exception:
                return float("nan")
            d, nd2 = {}, 0.0
            for k in float_vars:
                dk = float(y[k]) - float(x[k])
                d[k] = dk
                nd2 += dk * dk
            nd = math.sqrt(nd2)
            if not math.isfinite(nd) or nd == 0.0:
                return float("nan")
            acc += math.log(nd / d0)
            y = dict(x)
            for k in float_vars:
                y[k] = x[k] + (d0 / nd) * d[k]
    return acc / (n_steps * h)


def fingerprint_source(text, entry, h=None):
    """Run the legacy source and return (candidate_fingerprint, exec_info)."""
    ns = load_legacy(text, entry["dim"])
    v0 = ns["legacy_init"]()
    n_steps = ns["legacy_steps"](v0)
    float_vars = ns["FLOAT_VARS"]
    if not float_vars:
        raise TranslationError(0, "no float state variables to record")
    if h is None:
        h = float(entry.get("dt_lyap", 1.0)) if entry["family"] == "ode" else 1.0
    if entry["family"] == "ode":
        n_window = int(round(entry.get("t_lyap", 1000.0) / h))
    else:
        n_window = int(entry.get("n_lyap", 50000))
    n_window = max(1, n_window)

    v = dict(v0)
    pts, finite = [], True
    with np.errstate(all="ignore"):
        for i in range(n_steps):
            try:
                v = ns["legacy_step"](v, i)
            except Exception:
                finite = False
                break
            row = [float(v[k]) for k in float_vars]
            if not all(math.isfinite(r) for r in row):
                finite = False
                break
            pts.append(row)

    window = pts[-N_POINTS_D2:] if len(pts) >= N_POINTS_D2 else pts
    arr = np.asarray(window, dtype=float) if window else np.empty((0, len(float_vars)))
    if finite and window:
        tail = arr[len(arr) // 2:]
        uniq = np.unique(np.round(tail, 12), axis=0).shape[0]
        uf = uniq / float(tail.shape[0])
    else:
        uf = 0.0
    d2 = fingerprint_d2(arr, finite and len(window) > 0)

    if finite and pts:
        lam = lyap_legacy(ns["legacy_step"], v0, float_vars,
                          n_tr=n_steps, n_steps=n_window, h=h)
    else:
        lam = float("nan")

    cand = {"lambda1": float(lam),
            "corr_dim": d2,
            "unique_frac": float(uf),
            "finite": bool(finite)}
    info = {
        "status": "executed_source",
        "translator": "strict_f77_subset",
        "steps_run": n_steps,
        "recorded_vars": float_vars,
        "lambda1_transient": n_steps,
        "lambda1_window": n_window,
        "h": h,
        "translated_lines": len(ns["__source"].splitlines()),
        "translated_source": ns["__source"],
    }
    return cand, info
