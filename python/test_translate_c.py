#!/usr/bin/env python3
"""Branch tests for the strict C subset (translate_c.py) and its wiring in
parse_legacy_source.py.

Every refusal/acceptance branch of c_to_f77 must be exercised here
(AGENTS.md principle 14): lexical helpers, expression/assignment
translation, for-loop header variants, if/else chains, single-line ifs, all
decl/refusal paths, the line_map error remap (load_c) including its range
guard, and the three analyze() execution branches (executed_source,
translation_failed, unsupported_language).

Usage: python3 python/test_translate_c.py   # exit 0 only if all pass
"""
from __future__ import annotations

import math
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

import parse_legacy_source as pls                    # noqa: E402
import translate_c                                   # noqa: E402
from translate_c import (                            # noqa: E402
    _assign_stmt, _c_expr, _match_paren, _split_c, _split_top,
    _strip_comments, _strip_preproc, c_to_f77, fingerprint_c_source,
    load_c,
)
from translate_f77 import TranslationError           # noqa: E402

FIXTURE = ROOT / "fixtures" / "legacy" / "lorenz_euler_float.c"

CHECKS = 0
FAILS: list[str] = []


def check(label: str, cond: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    if not cond:
        FAILS.append(label)
        print(f"FAIL {label}" + (f" -- {detail}" if detail else ""))


def expect_err(label: str, contains: str, fn, *args, line=None) -> None:
    """fn(*args) must raise TranslationError whose message contains
    `contains` (and whose .line equals `line` when given)."""
    global CHECKS
    CHECKS += 1
    try:
        fn(*args)
    except TranslationError as e:
        if contains in str(e) and (line is None or e.line == line):
            return
        FAILS.append(label)
        print(f"FAIL {label} -- got {e!s} (expected {contains!r}, "
              f"line={line})")
    except Exception as e:                            # noqa: BLE001
        FAILS.append(label)
        print(f"FAIL {label} -- wrong exception {type(e).__name__}: {e}")
    else:
        FAILS.append(label)
        print(f"FAIL {label} -- no TranslationError raised")


def err_c(label: str, contains: str, src: str, line=None) -> None:
    expect_err(label, contains, c_to_f77, src, line=line)


def main() -> int:
    t0 = time.time()

    # ---- 1. lexical helpers --------------------------------------------
    out = _strip_comments("a; // drop\nb; /* x */ c;")
    check("strip: line count preserved", out.count("\n") == 1)
    check("strip: // and /* */ removed",
          out == "a; \nb;  c;" and "drop" not in out)
    out = _strip_comments('x = "a\\"b"; // gone')
    check("strip: escaped quote keeps string open",
          'x = "a\\"b";' in out and "gone" not in out)
    out = _strip_comments('u = "http://a//b"; // gone')
    check("strip: // inside string survives",
          "http://a//b" in out and "gone" not in out)
    expect_err("strip: unterminated /* comment",
               "unterminated /* comment", _strip_comments,
               "int x;\n/* oops\n", line=2)
    expect_err("strip: unterminated string literal",
               "unterminated string literal", _strip_comments,
               'x = "abc', line=1)
    check("preproc: line replaced by blank line",
          _strip_preproc("#define N 5\nint x;\n") == "\nint x;\n")
    check("preproc: keeps non-# lines",
          _strip_preproc("int x;\n") == "int x;\n")

    check("split: two statements keep their lines",
          _split_c("a = 1;\nb = 2;\n")
          == [("stmt", 1, "a = 1"), ("stmt", 2, "b = 2")])
    check("split: braces are their own items",
          _split_c("{ x = 1; }")
          == [("lbrace", 1, ""), ("stmt", 1, "x = 1"),
              ("rbrace", 1, "")])
    check("split: ; { } inside a string stay in the statement",
          _split_c('x = "a; b"; y = 2;')
          == [("stmt", 1, 'x = "a; b"'), ("stmt", 1, "y = 2")])
    items = _split_c("a =\n  1;")
    check("split: newline inside statement collapses to spaces",
          len(items) == 1 and items[0][0] == "stmt"
          and " ".join(items[0][2].split()) == "a = 1")
    check("match_paren: nested close found", _match_paren("(a(b))", 0) == 5)
    check("match_paren: unbalanced -> -1", _match_paren("(a(b", 0) == -1)
    check("split_top: depth-aware split",
          _split_top("a, f(b, c), d", ",") == ["a", " f(b, c)", " d"])

    # ---- 2. expressions and assignments -------------------------------
    check("expr: != -> /=", _c_expr("x != 0", 1) == "x /= 0")
    check("expr: && -> .and.",
          _c_expr("a && b", 1) == "a  .and.  b" and "&&" not in
          _c_expr("a && b", 1))
    check("expr: || -> .or.", ".or." in _c_expr("a || b", 1))
    check("expr: ! -> .not.", _c_expr("!a", 1) == " .not. a")
    check("expr: <= preserved", _c_expr("a <= b", 1) == "a <= b")
    check("expr: f suffix stripped", _c_expr("10.0f", 1) == "10.0")
    check("expr: bare double -> d0", _c_expr("10.0", 1) == "10.0d0")
    check("expr: bare integer unchanged", _c_expr("50000", 1) == "50000")
    check("expr: exponent with f suffix kept",
          _c_expr("1.5e3f", 1) == "1.5e3")
    check("expr: identifiers with digits untouched",
          _c_expr("x1 + beta2 + 2*x", 1) == "x1 + beta2 + 2*x")
    expect_err("expr: bare exponent refused",
               "double literal '1.0e5' with exponent unsupported",
               _c_expr, "1.0e5", 1, line=1)
    check("assign: plain", _assign_stmt("x = y + 1.0f", 1) == "x = y + 1.0")
    check("assign: += expands",
          _assign_stmt("x += y", 1) == "x = (x) + (y)")
    check("assign: /= expands",
          _assign_stmt("x /= 2", 1) == "x = (x) / (2)")
    check("assign: no '=' -> None", _assign_stmt("foo()", 1) is None)
    expect_err("assign: == misread as assignment",
               "relational == on the left of an assignment unsupported",
               _assign_stmt, "x == y", 1, line=1)
    expect_err("assign: chained assignment refused",
               "chained assignment unsupported",
               _assign_stmt, "x = y = 1.0f", 1, line=1)

    # ---- 3. acceptance: golden outputs --------------------------------
    golden = "\n".join([
        "      real x, y, z",
        "      real h",
        "      h = 0.01",
        "      real sigma, rho, beta",
        "      sigma = 10.0",
        "      rho = 28.0",
        "      beta = 8.0 / 3.0",
        "      real dx, dy, dz",
        "      integer i",
        "      x = 1.0",
        "      y = 1.0",
        "      z = 1.0",
        "      do i = (0) + 1, 50000",
        "      dx = sigma * (y - x)",
        "      dy = x * (rho - z) - y",
        "      dz = x * y - beta * z",
        "      x = x + h * dx",
        "      y = y + h * dy",
        "      z = z + h * dz",
        "      end do",
        "      end",
    ]) + "\n"
    golden_lm = [4, 5, 5, 9, 9, 9, 9, 10, 11, 13, 13, 13, 14, 15, 16,
                 17, 18, 19, 20, 21, 24]
    f77, lm = c_to_f77(FIXTURE.read_text())
    check("fixture: golden F77 text", f77 == golden,
          detail=repr(f77[:120]))
    check("fixture: golden line_map", lm == golden_lm, detail=repr(lm))
    check("fixed-form: every line is 6 spaces + code",
          all(l[:6].isspace() and len(l) > 6
              for l in f77.splitlines()))
    check("fixed-form: no comment char in column 1",
          all(l[0] not in "cC*!" for l in f77.splitlines()))

    chain = """int main(void)
{
    if (a > b) {
        x = 1.0f;
    } else if (b > c) {
        x = 2.0f;
    } else {
        x = 3.0f;
    }
}
"""
    f77, lm = c_to_f77(chain)
    check("if/else-if/else chain translates",
          f77.splitlines() == [
              "      if (a > b) then",
              "      x = 1.0",
              "      else if (b > c) then",
              "      x = 2.0",
              "      else",
              "      x = 3.0",
              "      end if",
              "      end",
          ], detail=repr(f77))
    check("chain: line_map follows the C lines",
          lm == [3, 4, 5, 6, 7, 8, 9, 10], detail=repr(lm))

    for incr in ("i++", "++i", "i += 1", "i = i + 1"):
        src = ("int main(void)\n{\n    int i;\n"
               f"    for (i = 0; i < 10; {incr}) {{\n"
               "        x = 1.0f;\n    }\n}\n")
        out, _ = c_to_f77(src)
        check(f"for increment {incr!r} emits do",
              "      do i = (0) + 1, 10" in out, detail=repr(out))

    src = """int main(void)
{
    for (int i = 0; i <= 10; ++i) {
        x = 1.0f;
    }
}
"""
    out, _ = c_to_f77(src)
    check("for: typed init declares the counter", "      integer i" in out)
    check("for: <= becomes inclusive bound",
          "      do i = (0) + 1, (10) + 1" in out, detail=repr(out))

    src = """#include <stdio.h>
int main(void)
{
    float x;
    printf("hi\\n");
    puts("yo");
    x = 1.0f;
    return 0;
}
"""
    out, _ = c_to_f77(src)
    check("io/return dropped", "printf" not in out and "puts" not in out
          and "return" not in out and "yo" not in out)
    check("body kept alongside dropped io", "      x = 1.0" in out)

    out, _ = c_to_f77("""int main(void)
{
    float a = 1.0f, b;
    double d;
}
""")
    check("decl: list with initializer",
          "      real a, b" in out and "      a = 1.0" in out)
    check("decl: double precision", "      double precision d" in out)

    out, _ = c_to_f77("""int main(void)
{
    if (a > b) a = b;
}
""")
    check("single-line if translates",
          "      if (a > b) a = b" in out, detail=repr(out))

    out, _ = c_to_f77("""int main(void)
{
    int i;
    for (i = 0; i < 10; i++) {
        if (a > b) {
            x = 1.0f;
        }
    }
}
""")
    check("nested if inside for closes before end do",
          out.index("      end if") < out.index("      end do"),
          detail=repr(out))

    # ---- 4. refusals: every raise() branch ----------------------------
    err_c("refuse: while", "unsupported C construct 'while'",
          "int main(void)\n{\n    while (1) { x = 1.0f; }\n}\n", line=3)
    err_c("refuse: do", "unsupported C construct 'do'",
          "int main(void)\n{\n    do { x = 1.0f; } while (0);\n}\n",
          line=3)
    err_c("refuse: goto", "unsupported C construct 'goto'",
          "int main(void)\n{\n    goto fail;\n}\n", line=3)
    err_c("refuse: labels", "labels unsupported in C subset",
          "int main(void)\n{\n    loop: x = 1.0f;\n}\n", line=3)
    err_c("refuse: array declaration", "C arrays unsupported",
          "int main(void)\n{\n    float a[10];\n}\n", line=3)
    err_c("refuse: array use", "C arrays unsupported",
          "int main(void)\n{\n    a[0] = 1.0f;\n}\n", line=3)
    err_c("refuse: pointer declaration",
          "unsupported declaration item '*p'",
          "int main(void)\n{\n    int *p;\n}\n", line=3)
    err_c("refuse: empty declaration", "empty declaration",
          "int main(void)\n{\n    float ;\n}\n", line=3)
    err_c("refuse: malformed declaration item",
          "unsupported declaration item 'a b'",
          "int main(void)\n{\n    float a b;\n}\n", line=3)
    err_c("refuse: function call", "function call 'foo' unsupported",
          "int main(void)\n{\n    foo(1);\n}\n", line=3)
    err_c("refuse: bare statement",
          "unsupported C statement: 'a b c'",
          "int main(void)\n{\n    a b c;\n}\n", line=3)
    err_c("refuse: chained assignment", "chained assignment unsupported",
          "int main(void)\n{\n    x = y = 1.0f;\n}\n", line=3)
    err_c("refuse: == as assignment",
          "relational == on the left of an assignment unsupported",
          "int main(void)\n{\n    x == 1.0f;\n}\n", line=3)
    err_c("refuse: bare exponent literal",
          "double literal '1.0e5' with exponent unsupported",
          "int main(void)\n{\n    x = 1.0e5;\n}\n", line=3)
    err_c("refuse: bare block", "unexpected '{'",
          "int main(void)\n{\n    {\n        x = 1.0f;\n    }\n}\n",
          line=3)
    err_c("refuse: unmatched '}'", "unmatched '}'",
          "int main(void)\n{\n}\n}\n", line=4)
    err_c("refuse: unclosed main", "unclosed 'main' block",
          "int main(void)\n{\n    x = 1.0f;\n", line=1)
    err_c("refuse: unclosed for", "unclosed 'for' block",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 3; i++) {\n        x = 1.0f;\n",
          line=4)
    err_c("refuse: statement after main",
          "statement/function after main()",
          "int main(void)\n{\n}\nx = 1.0f;\n", line=4)
    err_c("refuse: second main", "only one function (main) supported",
          "int main(void)\n{\n    int main(void)\n}\n", line=3)
    err_c("refuse: main signature junk", "unsupported main signature",
          "int main(void) junk\n{\n}\n", line=1)
    err_c("refuse: main without brace block",
          "main body must be followed by a brace block",
          "int main(void);\nx = 1.0f;\n", line=1)
    err_c("refuse: nested for", "nested for-loop unsupported",
          "int main(void)\n{\n"
          "    for (i = 0; i < 3; i++) {\n"
          "        for (j = 0; j < 3; j++) {\n"
          "            x = 1.0f;\n        }\n    }\n}\n", line=4)
    err_c("refuse: second for after loop closed",
          "only one for-loop supported",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 3; i++) {\n        x = 1.0f;\n    }\n"
          "    for (i = 0; i < 3; i++) {\n        x = 2.0f;\n    }\n}\n",
          line=7)
    err_c("refuse: statement after for-loop",
          "statement after the for-loop unsupported",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 3; i++) {\n        x = 1.0f;\n    }\n"
          "    x = 2.0f;\n}\n", line=7)
    err_c("refuse: loop counter read in body",
          "referenced in the for body",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 3; i++) {\n        x = i;\n    }\n}\n",
          line=5)
    err_c("refuse: declaration inside for body",
          "C declaration inside the for-loop body unsupported",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 3; i++) {\n        float y;\n    }\n}\n",
          line=5)
    err_c("refuse: for header junk after )",
          "for-loop without braces unsupported",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 3; i++) x = 1.0f;\n}\n", line=4)
    err_c("refuse: for without brace block",
          "for-loop body must be followed by a brace block",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 3; i++)\n", line=4)
    err_c("refuse: for header parts", "for header must be",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 3) {\n        x = 1.0f;\n    }\n}\n",
          line=4)
    err_c("refuse: for init", "unsupported for-init ''",
          "int main(void)\n{\n    int i;\n"
          "    for (; i < 3; i++) {\n        x = 1.0f;\n    }\n}\n",
          line=4)
    err_c("refuse: for condition", "unsupported for condition",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i == 3; i++) {\n        x = 1.0f;\n    }\n}\n",
          line=4)
    err_c("refuse: for init/cond variable mismatch",
          "for-init and for-condition use different variables",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; j < 3; i++) {\n        x = 1.0f;\n    }\n}\n",
          line=4)
    err_c("refuse: for step != +1",
          "unsupported for increment",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 10; i += 2) {\n"
          "        x = 1.0f;\n    }\n}\n", line=4)
    err_c("refuse: for incr variable mismatch",
          "for-init and for-increment use different variables",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 10; j++) {\n"
          "        x = 1.0f;\n    }\n}\n", line=4)
    err_c("refuse: for unbalanced parens",
          "unbalanced parentheses in for",
          "int main(void)\n{\n    int i;\n"
          "    for (i = 0; i < 3; i++\n", line=4)
    err_c("refuse: if unbalanced parens",
          "unbalanced parentheses in if",
          "int main(void)\n{\n    if (a > b\n", line=3)
    err_c("refuse: if without brace block",
          "if body must be followed by a brace block",
          "int main(void)\n{\n    if (a > b)\n", line=3)
    err_c("refuse: else after brace-less if",
          "else after brace-less if unsupported",
          "int main(void)\n{\n    if (a > b) a = b; else a = c;\n}\n",
          line=3)
    err_c("refuse: else without braces",
          "else without braces unsupported",
          "int main(void)\n{\n    if (a > b) {\n        x = 1.0f;\n"
          "    } else x = 2.0f;\n}\n", line=5)
    err_c("refuse: else-if without braces",
          "else-if without braces unsupported",
          "int main(void)\n{\n    if (a > b) {\n        x = 1.0f;\n"
          "    } else if (b > c) x = 2.0f;\n}\n", line=5)
    err_c("refuse: else-if unbalanced",
          "unbalanced parentheses in else-if",
          "int main(void)\n{\n    if (a > b) {\n        x = 1.0f;\n"
          "    } else if (b > c\n", line=5)
    err_c("refuse: stray else without if",
          "'else' without a matching 'if'",
          "int main(void)\n{\n    else {\n        x = 1.0f;\n    }\n}\n",
          line=3)
    err_c("refuse: else without brace block",
          "else branch must be followed by a brace block",
          "int main(void)\n{\n    if (a > b) {\n        x = 1.0f;\n"
          "    } else\n", line=5)
    err_c("refuse: single-line if with goto",
          "unsupported single-line if body 'goto fail'",
          "int main(void)\n{\n    if (a > b) goto fail;\n}\n", line=3)
    err_c("refuse: single-line if with call",
          "unsupported single-line if body 'foo()'",
          "int main(void)\n{\n    if (a > b) foo();\n}\n", line=3)
    err_c("refuse: comment-only file", "empty C translation",
          "/* only a comment */\n// and another\n", line=0)

    # ---- 5. load_c: F77 backstop errors remap to C lines --------------
    expect_err("load_c: undeclared variable maps to C line",
               "undeclared variable 'x'", load_c,
               "int main(void)\n{\n    x = 1.0f;\n}\n", line=3)
    expect_err("load_c: declaration after do-loop maps to C line",
               "declaration after do-loop", load_c,
               "int main(void)\n{\n    int i;\n    float x;\n"
               "    for (i = 0; i < 3; i++) {\n        x = 1.0f;\n    }\n"
               "    double y;\n}\n", line=8)
    expect_err("load_c: if condition must be logical maps to C line",
               "if condition must be logical", load_c,
               "int main(void)\n{\n    if (1) {\n        x = 1.0f;\n"
               "    }\n}\n", line=3)

    # range guard: out-of-range lines pass through unchanged
    orig = translate_c.translate

    def raiser(bad_line):
        def _f(_src, _dim=None):
            raise TranslationError(bad_line, "boom")
        return _f

    try:
        for bad_line, want in ((999, 999), (0, 0), (1, 3)):
            translate_c.translate = raiser(bad_line)
            expect_err(f"load_c: guard keeps line {bad_line} unremapped",
                       f"line {want}: boom", load_c,
                       "int main(void)\n{\n}\n", line=want)
    finally:
        translate_c.translate = orig

    # ---- 6. end-to-end on the fixture ---------------------------------
    entry = pls.load_entry("lorenz")
    cand, info = fingerprint_c_source(FIXTURE.read_text(), entry, h=0.01)
    check("e2e: executed_source", info["status"] == "executed_source")
    check("e2e: translator tag", info["translator"] == "strict_c_subset")
    check("e2e: structured mode", info["mode"] == "structured")
    check("e2e: steps_run", info["steps_run"] == 50000)
    check("e2e: recorded_vars",
          info["recorded_vars"] == ["x", "y", "z"])
    check("e2e: h picked up", math.isclose(info["h"], 0.01))
    check("e2e: finite True", cand["finite"] is True)
    check("e2e: lambda1 NaN (single precision)",
          not math.isfinite(cand["lambda1"]))
    check("e2e: corr_dim matches the F77 twin",
          math.isclose(cand["corr_dim"], 1.3437669536670236,
                       rel_tol=0, abs_tol=1e-9),
          detail=repr(cand["corr_dim"]))

    rep = pls.analyze(FIXTURE)
    check("analyze: executed_source", rep["candidate_source"]
          == "executed_source")
    check("analyze: language c", rep["language"] == "c")
    check("analyze: variant euler_single",
          rep["variant"] == "euler_single", detail=repr(rep["variant"]))
    check("analyze: execution status",
          rep["execution"]["status"] == "executed_source")
    check("analyze: diagnosis precision_insuficiente",
          rep["diagnosis"]["diagnostico"] == "precision_insuficiente",
          detail=repr(rep["diagnosis"]))
    check("analyze: fingerprint finite True",
          rep["diagnosis"]["fingerprint"]["finite"] is True)

    # ---- 7. wiring branches: translation_failed / unsupported_language
    bad_c = """/* Lorenz -- forward Euler, single precision */
#include <stdio.h>

float x, y, z;
float h = 0.01f;
float sigma = 10.0f, rho = 28.0f, beta = 8.0f / 3.0f;

int main(void)
{
    int i;
    x = 1.0f;
    while (i < 50000) {
        x = x + h * sigma;
    }
    return 0;
}
"""
    with tempfile.TemporaryDirectory() as td:
        bad_path = Path(td) / "lorenz_bad.c"
        bad_path.write_text(bad_c)
        rep = pls.analyze(bad_path)
        check("wiring: bad C -> translation_failed",
              rep["execution"]["status"] == "translation_failed",
              detail=repr(rep["execution"]))
        check("wiring: refusal line remapped",
              rep["execution"]["line"] == 12
              and "unsupported C construct 'while'"
              in rep["execution"]["error"],
              detail=repr(rep["execution"].get("error")))
        check("wiring: fallback note recorded",
              any("translation refused at line 12" in n
                  for n in rep["assumptions"]),
              detail=repr(rep["assumptions"]))
        check("wiring: fallback keeps synthetic source",
              rep["candidate_source"] == "synthetic_equivalent")
        check("wiring: synthetic diagnosis still produced",
              rep["diagnosis"] is not None
              and rep["diagnosis"]["diagnostico"]
              == "precision_insuficiente",
              detail=repr(rep["diagnosis"]))

        sysm = pls.detect_system(
            "sigma = 10.0; rho = 28.0; beta = 8.0 / 3.0\n")
        check("detect_system: const-set fallback finds lorenz",
              sysm["name"] == "lorenz", detail=repr(sysm))
        unknown = Path(td) / "mystery.txt"
        unknown.write_text("hello world\n")
        rep = pls.analyze(unknown)
        check("wiring: unknown language detected",
              rep["language"] == "unknown")
        check("wiring: unsupported_language status",
              rep["execution"]["status"] == "unsupported_language"
              and "no translator" in rep["execution"]["note"],
              detail=repr(rep["execution"]))
        check("wiring: no runnable variant -> no diagnosis",
              rep["variant"] is None and rep["diagnosis"] is None
              and "no runnable equivalent" in rep["note"],
              detail=repr(rep.get("note")))

    # ---- summary -------------------------------------------------------
    dt = time.time() - t0
    if FAILS:
        print(f"\n{len(FAILS)}/{CHECKS} checks FAILED: "
              + ", ".join(FAILS))
        return 1
    print(f"ALL {CHECKS} CHECKS PASS ({dt:.1f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
