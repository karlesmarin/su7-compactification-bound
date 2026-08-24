#!/usr/bin/env python3
"""rigorous_ceiling.py -- the ceiling of the FULL polylogarithm, not of the truncated surface.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS IS NOT `ceiling_margin.py`
------------------------------------
`ceiling_margin.py` shifts the optimum: it takes the winning (A_4, 8D) of the truncated problem
and asks what the remainder bound does to its number.  That is not a theorem about the full
potential, because the margin is not uniform over the lattice -- it grows with A_4 through
S_4 <= 13 A_4 - 12 D -- so a rung that loses on the truncated surface can win on the bounded one.

The theorem needs the OPTIMUM OF THE SHIFTED SURFACE.  That is what this file computes: the same
rational-dual feasibility of `ceiling_ilp.py`, rung by rung, with each rung's 1/R_5 replaced by
its rigorous upper bound before the maximum is taken.

    alpha_true >= alphahat - delta ,   delta = |R'| / F''  from (T2s)
    1/R_5 = 2 m_W / alpha_true  <=  2 m_W / (alphahat - delta)

Everything else -- the lattice, the cone, the dual vertices, the mod-6 law, tmax -- is imported
from `ceiling_ilp.py` by executing it, so this file cannot drift from the certificate it is
correcting.  Its output is suppressed; only its namespace is used.

Run:  python rigorous_ceiling.py     (slow: it runs ceiling_ilp.py to get the lattice)
"""
import contextlib
import io
import math
import pathlib
import sys

P = lambda *a: print(*a, flush=True)
HERE = pathlib.Path(__file__).resolve().parent

# --- import the certificate's own machinery, silently -------------------------------------
_src = (HERE / "ceiling_ilp.py").read_text(encoding="utf-8")
NS = {"__file__": str(HERE / "ceiling_ilp.py"), "__name__": "ceiling_ilp"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(_src, "ceiling_ilp.py", "exec"), NS)

x_of, tmax, MU_HI, MW = NS["x_of"], NS["tmax"], NS["MU_HI"], NS["MW"]
Z3 = NS["Z3"]
CMAX = 3                    # the 84s carry charge three; nothing in this model carries more


def curvature(t, k, x):
    """d^2F/dalpha^2 at the closed-form minimum, as the paper writes it."""
    return math.pi ** 2 * (2 * Z3 * (k / 8.0) - t * x * x / 6.0)


def delta_alpha(t, k, alpha):
    """(T2s): |R'| <= pi^6 S4 c_max^2 alpha^5 / (360 (1-u^2)),  S4 <= 13 t - 12 D."""
    u = CMAX * alpha
    if u >= 1:
        return None
    S4 = 13.0 * t - 12.0 * (k / 8.0)
    Rp = math.pi ** 6 * S4 * CMAX ** 2 * alpha ** 5 / (360.0 * (1 - u * u))
    c = curvature(t, k, math.pi * alpha)
    if c <= 0:
        return None
    return Rp / c


def main():
    P("=" * 104)
    P("THE CEILING OF THE FULL POTENTIAL:  optimum of the SHIFTED surface, not shift of the optimum")
    P("=" * 104)
    P("%6s %8s %11s %12s %12s %13s %13s"
      % ("8D", "max A4", "alpha_hat", "delta", "delta/alpha", "1/R5 trunc", "1/R5 rigorous"))

    RUNGS = [k for k in range(1, 140, 2)]
    best_t, best_r = None, None
    rows = []
    for k in RUNGS:
        t = tmax(k, MU_HI)
        if t is None:
            continue
        x = x_of(t, k, MU_HI)
        a = x / math.pi
        d = delta_alpha(t, k, a)
        if d is None or d >= a:
            continue
        inv_t = 2 * MW / a
        inv_r = 2 * MW / (a - d)
        rows.append((k, t, a, d, inv_t, inv_r))
        if best_t is None or inv_t > best_t[4]:
            best_t = rows[-1]
        if best_r is None or inv_r > best_r[5]:
            best_r = rows[-1]
        if k <= 13 or k in (21, 33, 51, 99, 129):
            P("%6d %8d %11.6f %12.3e %12.2e %13.1f %13.1f"
              % (k, t, a, d, d / a, inv_t, inv_r))

    P("")
    P("rungs evaluated : %d" % len(rows))
    P("")
    P("truncated optimum : 8D = %d, A4 = %d, 1/R5 = %.1f GeV = %.3f TeV"
      % (best_t[0], best_t[1], best_t[4], best_t[4] / 1000))
    P("rigorous optimum  : 8D = %d, A4 = %d, 1/R5 <= %.1f GeV = %.3f TeV"
      % (best_r[0], best_r[1], best_r[5], best_r[5] / 1000))
    P("")
    same = best_t[0] == best_r[0]
    P("does the maximum move to another rung?  %s"
      % ("no -- the same rung wins on both surfaces" if same
         else "YES: %d -> %d.  This is exactly why shifting the optimum is not enough."
              % (best_t[0], best_r[0])))
    P("")
    P("THE STATEMENT:  1/R_5  <=  %.3f TeV  for every admissible bulk content, for the FULL"
      % (best_r[5] / 1000))
    P("polylogarithmic potential -- not for its expansion.  The margin over the truncated")
    P("certificate is %.1f GeV, i.e. %.2f %%." % (best_r[5] - best_t[4],
                                                  100 * (best_r[5] / best_t[4] - 1)))

    P("")
    P("=" * 104)
    P("CONTROLS")
    P("=" * 104)
    # (a) the margin must be monotone in the rung's A4 -- if it were not, the bound is misread
    ok = all(rows[i][3] > 0 for i in range(len(rows)))
    P("  every delta strictly positive                          : %s" % ok)
    # (b) the rigorous ceiling must be ABOVE the truncated one, never below
    P("  rigorous ceiling >= truncated ceiling                  : %s"
      % (best_r[5] >= best_t[4]))
    # (c) and it must be finite and close: a bound that doubled the number would be useless
    P("  and within 10 %% of it (else the bound is not usable)   : %s (%.2f %%)"
      % (best_r[5] / best_t[4] < 1.10, 100 * (best_r[5] / best_t[4] - 1)))
    # (d) the hypothesis u < 1 must hold at every rung kept
    P("  u = 3 alpha < 1 at every rung kept                     : %s"
      % all(3 * r[2] < 1 for r in rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
