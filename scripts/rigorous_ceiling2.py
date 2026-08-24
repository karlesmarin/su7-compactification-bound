#!/usr/bin/env python3
"""rigorous_ceiling2.py -- the remainder as a LINEAR FUNCTIONAL, inside the certificate's own dual.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHAT CHANGED, AND WHY THE FIRST ATTEMPT FAILED
-----------------------------------------------
`rigorous_ceiling.py` bounded the remainder by (N, c_max) through S_4 <= 13 A_4 - 12 D, and the
optimum of the shifted surface came out at 52 TeV -- vacuous.  The looseness was one inequality:
S_4 <= 13 A_4 chains "all antiperiodic charge at c = 3" with "all periodic charge at c = 1", two
extremes at once.

The repair is not a better inequality.  It is noticing that the bound never needed one.

    Charges run over c in {1,2,3} and two parities, so the potential lives in a SIX-dimensional
    function space however many multiplets there are.  The remainder is therefore

        R(x) = R_gauge(x) + sum_j m_j R_j(x)

    EXACTLY LINEAR in the multiplicities -- the same shape as A_4, 8D and G.  So the margin does
    not need bounding by moments at all: it is its own linear functional, and it goes into the
    same dual that certifies the ceiling.

    margin_num(t,k) = max { |R'_gauge| + sum_j m_j |R'_j(x)| : sum m_j a_j = T, sum m_j k_j = Q, m >= 0 }

    a two-equality LP, so its dual has two variables, exactly like `gmin_cone`.  The G <= G*
    constraint is DROPPED, which only makes the bound larger -- still a bound, and it keeps the
    dual at two variables.

Run:  python rigorous_ceiling2.py     (slow: runs ceiling_ilp.py for the lattice)
"""
import contextlib
import io
import math
import pathlib
import sys

import numpy as np
from scipy.optimize import linprog

P = lambda *a: print(*a, flush=True)
HERE = pathlib.Path(__file__).resolve().parent

_src = (HERE / "ceiling_ilp.py").read_text(encoding="utf-8")
NS = {"__file__": str(HERE / "ceiling_ilp.py"), "__name__": "ceiling_ilp"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(_src, "ceiling_ilp.py", "exec"), NS)

terms, GAUGE, REPS, NAMES = NS["terms"], NS["GAUGE"], NS["REPS"], NS["NAMES"]
AV, KV, A0, K0 = NS["AV"], NS["KV"], NS["A0"], NS["K0"]
x_of, tmax, MU_HI, MW, Z3 = NS["x_of"], NS["tmax"], NS["MU_HI"], NS["MW"], NS["Z3"]

Z5 = 1.0369277551433699
ETA3 = 0.75 * Z3
ETA5 = 15.0 / 16.0 * Z5
H4 = 25.0 / 12.0
LN2 = math.log(2.0)
NMAX = 4000                       # windings; the tail beyond is below 1e-14 at these x


def reli5(s, y):
    n = np.arange(1, NMAX + 1)
    return float((np.cos(n * y) * (s ** n) / n ** 5.0).sum()) if s < 0 else \
        float((np.cos(n * y) / n ** 5.0).sum())


def _reli5(s, y):
    n = np.arange(1, NMAX + 1)
    sg = np.ones(NMAX) if s > 0 else (-1.0) ** n
    return float((np.cos(n * y) * sg / n ** 5.0).sum())


def head_term(m, s, c, x):
    """the three rungs the paper keeps, for ONE term (m, s, c)."""
    if s > 0:
        return m * (Z5 - x * x / 2 * c * c * Z3
                    + x ** 4 / 24 * c ** 4 * (H4 - math.log(c * x)))
    return m * (-ETA5 + x * x / 2 * c * c * ETA3 - x ** 4 / 24 * c ** 4 * LN2)


def R_of(tt, x):
    """the exact remainder of a term table, F - (three rungs)."""
    return sum(m * _reli5(s, c * x) - head_term(m, s, c, x) for (m, s, c) in tt)


def dR_of(tt, alpha, h=1e-6):
    """d/d(alpha) of the remainder, by a central difference in alpha."""
    return (R_of(tt, math.pi * (alpha + h)) - R_of(tt, math.pi * (alpha - h))) / (2 * h)


TT = [terms(r, e, p) for (r, e, p) in REPS]


def curvature(t, k, x):
    return math.pi ** 2 * (2 * Z3 * (k / 8.0) - t * x * x / 6.0)


def margin_num(t, k, alpha):
    """max |R'| over the real relaxation at (A4,8D) = (t,k).  Linear program, 8 variables."""
    rho = np.array([abs(dR_of(tt, alpha)) for tt in TT])
    Aeq = np.array([[float(a) for a in AV], [float(x) for x in KV]])
    beq = np.array([float(t - A0), float(k - K0)])
    r = linprog(-rho, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * 8, method="highs")
    if not r.success:
        return None
    return abs(dR_of(GAUGE, alpha)) + float(-r.fun)


def main():
    P("=" * 106)
    P("THE REMAINDER AS A LINEAR FUNCTIONAL, MAXIMISED INSIDE THE CERTIFICATE'S OWN LP")
    P("=" * 106)

    P("per-multiplet remainder at the certifying alpha = 0.016025 (exact, not bounded):")
    P("  %-10s %14s %16s" % ("multiplet", "|R_j|", "|R'_j| (d/dalpha)"))
    a0 = 0.016025
    for nm, tt in zip(NAMES, TT):
        P("  %-10s %14.3e %16.3e" % (nm, abs(R_of(tt, math.pi * a0)), abs(dR_of(tt, a0))))
    P("  %-10s %14.3e %16.3e" % ("gauge", abs(R_of(GAUGE, math.pi * a0)), abs(dR_of(GAUGE, a0))))

    P("")
    P("%6s %8s %11s %12s %12s %13s %13s"
      % ("8D", "max A4", "alpha_hat", "delta", "delta/alpha", "1/R5 trunc", "1/R5 rigorous"))
    best_t, best_r, rows = None, None, []
    for k in range(1, 140, 2):
        t = tmax(k, MU_HI)
        if t is None:
            continue
        x = x_of(t, k, MU_HI)
        a = x / math.pi
        num = margin_num(t, k, a)
        c = curvature(t, k, x)
        if num is None or c <= 0:
            continue
        d = num / c
        if d >= a:
            continue
        inv_t, inv_r = 2 * MW / a, 2 * MW / (a - d)
        rows.append((k, t, a, d, inv_t, inv_r))
        if best_t is None or inv_t > best_t[4]:
            best_t = rows[-1]
        if best_r is None or inv_r > best_r[5]:
            best_r = rows[-1]
        if k <= 13 or k in (21, 33, 51, 99, 129):
            P("%6d %8d %11.6f %12.3e %12.2e %13.1f %13.1f"
              % (k, t, a, d, d / a, inv_t, inv_r))

    # ---- Y QUE PASA MAS ALLA DEL PELDANO 70, que es la pregunta que el articulo no contestaba.
    # El techo TRUNCADO tiene prueba analitica para todo k.  La cota del RESTO no: se afloja con
    # el peldano, y por encima de cierto punto certifica menos que el propio 10.037.  Esto lo
    # mide en vez de suponerlo, porque una version anterior del articulo decia "arbitrary
    # content" apoyandose en estos setenta.  [[a-bound-is-not-a-computation]]
    P("")
    P("BEYOND THE SEVENTY: the remainder bound loosens with the rung")
    P("%6s %11s %12s %13s %13s" % ("8D", "alpha_hat", "delta/alpha", "1/R5 trunc", "1/R5 rigorous"))
    far = []
    for k in (1, 129, 201, 501, 1201):
        t = tmax(k, MU_HI)
        if t is None:
            continue
        x = x_of(t, k, MU_HI)
        a = x / math.pi
        num = margin_num(t, k, a)
        c = curvature(t, k, x)
        if num is None or c <= 0:
            continue
        d = num / c
        if d >= a:
            P("%6d %11.6f %12.1e %13s %13s" % (k, a, d / a, "--", "bound is vacuous"))
            far.append((k, a, d / a, None, None))
            continue
        far.append((k, a, d / a, 2 * MW / a, 2 * MW / (a - d)))
        # una cifra, que es la que el articulo cita: la compuerta de numeros compara la mantisa
        # token a token, y un 9.12e-02 archivado no respalda un 9.1 impreso.  [[quote-the-number-from-the-aux]]
        P("%6d %11.6f %12.1e %13.1f %13.1f" % (k, a, d / a, far[-1][3], far[-1][4]))
    if len(far) >= 2:
        lo, hi = far[0], far[-1]
        P("")
        P("  delta/alpha grows by a factor %.0f from 8D = %d to 8D = %d"
          % (hi[2] / lo[2], lo[0], hi[0]))
        if hi[4]:
            P("  and the certified ceiling there is %.1f TeV -- NOT because the ceiling rises,"
              % (hi[4] / 1000))
            P("  but because the remainder bound stops saying anything.  So eq. (rigorousceiling)")
            P("  holds over the rungs it was run on and not for arbitrary content.")

    P("")
    P("rungs evaluated : %d" % len(rows))
    P("truncated optimum : 8D = %d, A4 = %d, 1/R5 = %.1f GeV = %.3f TeV"
      % (best_t[0], best_t[1], best_t[4], best_t[4] / 1000))
    P("rigorous optimum  : 8D = %d, A4 = %d, 1/R5 <= %.1f GeV = %.3f TeV"
      % (best_r[0], best_r[1], best_r[5], best_r[5] / 1000))
    P("")
    P("maximum stays on the same rung : %s" % (best_t[0] == best_r[0]))
    P("margin over the certificate    : %.1f GeV  (%.3f %%)"
      % (best_r[5] - best_t[4], 100 * (best_r[5] / best_t[4] - 1)))

    P("")
    P("=" * 106)
    P("CONTROLS")
    P("=" * 106)
    P("  rigorous >= truncated                         : %s" % (best_r[5] >= best_t[4]))
    P("  every delta > 0                               : %s" % all(r[3] > 0 for r in rows))
    P("  usable (within 10%% of the certificate)        : %s (%.3f %%)"
      % (best_r[5] / best_t[4] < 1.10, 100 * (best_r[5] / best_t[4] - 1)))
    P("  and it beats the (N, c_max) bound of v1        : v1 gave 52.429 TeV")
    return 0


if __name__ == "__main__":
    sys.exit(main())
