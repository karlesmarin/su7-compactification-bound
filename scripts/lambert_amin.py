#!/usr/bin/env python3
"""lambert_amin.py -- the fixed point is a closed form after all, and it is a Lambert W.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS EXISTS.  The paper called the stationarity condition a "fixed-point equation rather than
a root of a polynomial" and solved it by iteration.  It is solvable in closed form, by the Lambert
W function, and the surrounding narrative -- "no polynomial, therefore no closed form" -- was
correspondingly too strong.  \\S\\ref{sec:ceiling} already solves the asymptotic ray with a Lambert
W, so the function was in the paper before the equation it belongs to was recognised.

THE ALGEBRA.  With x = pi alpha, the stationarity condition is

    x^2 (4G - A_4 (4 ln x + 1)) = 24 zeta(3) D .

Put y = x^2, so ln x = (1/2) ln y:

    y (4G - A_4 - 2 A_4 ln y) = 24 zeta(3) D
    y (b - ln y) = d ,        b = (4G - A_4) / (2 A_4) ,     d = 12 zeta(3) D / A_4 .

Set u = ln y - b.  Then y = e^{b+u} and -u e^{b+u} = d, so u e^u = -d e^{-b} and u = W(-d e^{-b}).
Using e^{W(z)} = z / W(z),

    y = e^b e^{W(z)} = e^b z / W(z) = -d / W(z) ,     z = -d e^{-b} ,

so

    x^2 = -d / W(-d e^{-b}) .

WHICH BRANCH, and it matters.  z is negative on every content here, so W has two real branches,
W_0 on [-1/e, 0) and W_{-1} on [-1/e, -1).  They are the two roots of the same stationarity
condition, and only one of them is the electroweak minimum -- the other is the outer stationary
point the paper's own \\S\\ref{sec:ceiling} keeps meeting as the "upper root".  The branch is
selected below by comparison against the numerically minimised potential, not assumed.

WHAT THIS DOES NOT CHANGE.  Everything the paper says about the ABSENCE OF A POLYNOMIAL stands:
the x^4 ln x branch point is still there and still rules a polynomial out, which is what
\\S\\ref{sec:parity} rests on.  What has to go is any suggestion that no closed form exists.  The
distinction is the one that matters: no polynomial is a theorem, no closed form was an
overstatement.

Run:  python lambert_amin.py
"""
import cmath
import json
import math
import pathlib
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

try:
    from scipy.special import lambertw
except ImportError:                      # keep the file runnable with no scipy
    def lambertw(z, k=0):
        w = cmath.log(-z) if k else (0j + 0.3)
        for _ in range(200):
            e = cmath.exp(w)
            w = w - (w * e - z) / (e * (w + 1) - (w + 2) * (w * e - z) / (2 * w + 2))
        return w


def lambert_x(A4, D, G, branch):
    b = (4.0 * G - A4) / (2.0 * A4)
    d = 12.0 * Z3 * D / A4
    z = -d * math.exp(-b)
    w = lambertw(z, branch)
    if abs(w.imag) > 1e-12:
        return None
    y = -d / w.real
    return math.sqrt(y) if y > 0 else None


# ================================================================= 1
P("=" * 100)
P("1 -- THE CLOSED FORM, AGAINST THE FIXED POINT IT REPLACES")
P("=" * 100)
P("  x^2 = -d / W(-d e^{-b}),   b = (4G - A_4)/(2 A_4),   d = 12 zeta(3) D / A_4")
P("")
P("  %-6s %16s %16s %12s %8s" % ("row", "iterated x", "Lambert-W x", "rel diff", "branch"))
rows, worst = [], 0.0
for lbl, content, a_them, mh_them, invR in T1:
    a, mo = closed_form(content)
    x = a * math.pi
    A4, D, G = float(mo["A4"]), float(mo["D"]), float(mo["G"])
    got = None
    for br in (0, -1):
        xw = lambert_x(A4, D, G, br)
        if xw is not None and abs(xw - x) / x < 1e-6:
            got = (xw, br)
            break
    assert got, "no real Lambert branch reproduces the fixed point on row %s" % lbl
    xw, br = got
    rel = abs(xw - x) / x
    worst = max(worst, rel)
    rows.append(dict(row=lbl, x=x, x_lambert=xw, rel=rel, branch=br))
    P("  %-6s %16.10f %16.10f %12.2e %8d" % (lbl, x, xw, rel, br))
P("")
P("  worst relative difference over the five rows : %.2e" % worst)
assert worst < 1e-12, "the Lambert form does not reproduce the fixed point to machine precision"
P("  So the two are the same number, and the iteration was never necessary.")

# ================================================================= 2
P("")
P("=" * 100)
P("2 -- BOTH BRANCHES ARE STATIONARY POINTS, AND THE OTHER ONE IS ALREADY IN THE PAPER")
P("=" * 100)
P("  W has two real branches for negative argument, so the closed form returns TWO stationary")
P("  points, not one.  That is not a defect of the formula -- it is the two roots the ceiling")
P("  argument of \\S\\ref{sec:ceiling} keeps meeting, now with names.")
P("")
P("  %-6s %16s %16s" % ("row", "W_0 branch x", "W_{-1} branch x"))
for lbl, content, a_them, mh_them, invR in T1:
    a, mo = closed_form(content)
    A4, D, G = float(mo["A4"]), float(mo["D"]), float(mo["G"])
    v = [lambert_x(A4, D, G, br) for br in (0, -1)]
    P("  %-6s %16s %16s"
      % (lbl, "%.10f" % v[0] if v[0] else "none", "%.10f" % v[1] if v[1] else "none"))
P("")
P("  CONTROL -- the branch is CHOSEN by comparison with the minimised potential and not")
P("  assumed; section 1 fails loudly if neither branch reproduces it.")

# ================================================================= 3
P("")
P("=" * 100)
P("3 -- WHAT DOES NOT CHANGE: THERE IS STILL NO POLYNOMIAL")
P("=" * 100)
P("  A closed form and a polynomial are different claims, and only the second is ruled out.")
P("  The potential carries an x^4 ln x branch point -- \\S\\ref{sec:closed} derives it from the")
P("  expansion -- and no polynomial has a branch point.  The Lambert form does not remove the")
P("  branch point; it is BUILT from it, since W is exactly the inverse of w e^w and the")
P("  logarithm is what puts it there.  Concretely, the argument -d e^{-b} carries G and A_4")
P("  through b, which is where the logarithm went.")
P("")
P("  the transcendental content is unchanged : zeta(5), zeta(3), ln x, ln 2, ln 3 -- and now")
P("  the name of the function that inverts the last one.")

(HERE / "outputs").mkdir(exist_ok=True)
(HERE / "outputs" / "lambert_amin.json").write_text(
    json.dumps(dict(rows=rows, worst=worst), indent=1), encoding="utf-8")
P("")
P("archived: outputs/lambert_amin.json")
