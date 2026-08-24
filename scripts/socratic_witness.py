#!/usr/bin/env python3
"""socratic_witness.py -- the witness has no interior minimum, and the question is why.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHAT HAPPENED
-------------
rigorous_1001.py went to evaluate the ceiling's witness -- 17x7(+,+) + 2x7(+,-) + 57x28(+,-),
the content semigroup.py proved attains (A4, 8D) = (212, 1) -- on the exact polylogarithmic
potential, and numeric_min returned None.  No interior minimum.

That is not a bug to patch.  It is a measurement, and it contradicts something the ceiling
assumes.  Carles: "pregunta con metodo socratico, hasta raiz."  So:

  Q1  8D = 1 > 0, so D > 0, so the origin is a MAXIMUM of the potential.  Then a minimum must
      exist somewhere.  Where did it go?
  Q2  If it left through alpha = 1, the electroweak vacuum this paper computes is not the
      global one.  Is there still a LOCAL minimum at small alpha?
  Q3  If there is, the closed form describes it correctly and the ceiling is a bound on a
      LOCAL vacuum.  Does the paper ever say that?
  Q4  Is this special to the witness, or does the ceiling's whole optimisation walk into it?
      The certificate maximises 1/R5 = 2 m_W / alpha, i.e. it PUSHES alpha down; and small
      alpha at fixed D means large A4, which means a large content.  Does the electroweak
      minimum stop being global exactly where the certificate wants to sit?
  Q5  And the root: the five published rows were checked for false vacua (Part VII's own
      false_vacuum_km.py) and are clean.  Was that check ever run at the ceiling?

Every question is measured below, in order, and the answers decide whether 10.01 TeV is a
bound on the vacuum or on a vacuum.

Run:  python socratic_witness.py
"""
import math
import pathlib
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4
MH_LO, MH_HI = 125.0, 127.0
NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]

cont = lambda v: [(REPS[j][0], REPS[j][1], REPS[j][2], v[j]) for j in range(8) if v[j]]
show = lambda v: " + ".join("%dx%s" % (v[j], NAMES[j]) for j in range(8) if v[j])


def stationary(v, lo=1e-5, hi=1.0, n=400001):
    """every local minimum of the EXACT F on (lo, hi], with its alpha and depth."""
    c = cont(v)
    xs = np.linspace(lo, hi, n)
    ys = F(c, xs)
    mins = []
    for i in range(1, len(xs) - 1):
        if ys[i] <= ys[i - 1] and ys[i] < ys[i + 1]:
            a, b = xs[i - 1], xs[i + 1]
            for _ in range(60):
                t = np.linspace(a, b, 15)
                j = int(np.argmin(F(c, t)))
                a, b = t[max(j - 1, 0)], t[min(j + 1, 14)]
            am = 0.5 * (a + b)
            mins.append((float(am), float(F(c, np.array([am]))[0])))
    edge = (float(xs[-1]), float(ys[-1]))
    return mins, edge, float(ys[0])


def mh_at(v, a, h=2e-4):
    c = cont(v)
    f = F(c, np.array([a - h, a, a + h]))
    fpp = float((f[0] - 2 * f[1] + f[2]) / h ** 2)
    return (KK * math.sqrt(fpp) / a if fpp > 0 else None), fpp


WIT = [0] * 8
WIT[0], WIT[1], WIT[3] = 17, 2, 57

P("=" * 100)
P("Q1 -- 8D = 1 > 0, SO THE ORIGIN IS A MAXIMUM.  WHERE IS THE MINIMUM?")
P("=" * 100)
mo = moments(cont(WIT))
P("  witness : %s" % show(WIT))
P("  A4 = %.0f   8D = %.0f   G = %.3f   N = %d" % (mo["A4"], 8 * mo["D"], mo["G"], sum(WIT)))
P("")
mins, edge, val0 = stationary(WIT)
P("  F at alpha -> 0      : %.6e" % val0)
P("  F at alpha = 1       : %.6e" % edge[1])
P("  interior local minima: %d" % len(mins))
for a, y in mins[:8]:
    P("     alpha = %.6f   F = %.6e" % (a, y))
P("")
if not mins:
    P("  ==> NONE.  The potential falls monotonically out to the edge of the fundamental domain.")
    P("      numeric_min was right to refuse, and the ceiling's witness has no electroweak")
    P("      vacuum at all on the exact potential.")
else:
    deepest = min(mins, key=lambda t: t[1])
    P("  ==> the deepest interior minimum is at alpha = %.6f" % deepest[0])
    P("      and the edge value is %s it" % ("BELOW" if edge[1] < deepest[1] else "above"))

P("")
P("=" * 100)
P("Q2 -- IS THERE STILL A LOCAL MINIMUM AT SMALL ALPHA, WHERE THE CLOSED FORM LIVES?")
P("=" * 100)
cf, _ = closed_form(cont(WIT))
P("  the closed form predicts a stationary point at alpha = %.6f" % cf)
sm = [m for m in mins if m[0] < 1.0 / 3]
P("  interior minima inside the radius of convergence (alpha < 1/3) : %d" % len(sm))
for a, y in sm[:5]:
    mh, fpp = mh_at(WIT, a)
    P("     alpha = %.6f   F'' = %+.4e   m_h = %s"
      % (a, fpp, ("%.2f GeV" % mh) if mh else "no real m_h"))
if not sm:
    P("  ==> the closed form's stationary point does NOT exist on the exact potential.  It is an")
    P("      artefact of truncating at x^4: the quartic turns the potential around, the full")
    P("      polylogarithm does not.")

P("")
P("=" * 100)
P("Q3/Q4 -- WHERE DOES IT BREAK?  walk the rung 8D = 1 from the five rows' scale upward")
P("=" * 100)
P("  The five published rows sit at A4 ~ 190-270 with SMALL contents.  The certificate wants")
P("  A4 = 212 with a content of 76.  Same A4, different size -- so size is the variable to")
P("  sweep, not A4.  Take the family n x 28(+,-) plus what is needed to hold 8D = 1.")
P("")
P("  %-6s %-34s %8s %10s %12s %12s %s"
  % ("N", "content", "A4", "8D", "minima", "alpha_min", "m_h"))
for n28 in (5, 10, 20, 30, 40, 50, 57):
    # 8D = -27 + 8*n7m + 2*n28m + (-6)*n7p = 1  ->  pick n7m = 2, solve n7p
    n7m = 2
    k = -27 + 8 * n7m + 2 * n28
    if (k - 1) % 6:
        continue
    n7p = (k - 1) // 6
    if n7p < 0:
        continue
    v = [0] * 8
    v[0], v[1], v[3] = n7p, n7m, n28
    m2 = moments(cont(v))
    mm, ed, _ = stationary(v, n=120001)
    a = mm[0][0] if mm else None
    mh = mh_at(v, a)[0] if a else None
    P("  %-6d %-34s %8.0f %10.0f %12d %12s %s"
      % (sum(v), show(v), m2["A4"], 8 * m2["D"], len(mm),
         ("%.6f" % a) if a else "--", ("%.1f" % mh) if mh else "--"))

P("")
P("=" * 100)
P("Q5 -- AND THE FIVE PUBLISHED ROWS, WHICH IS THE CONTROL THAT MUST PASS")
P("=" * 100)
P("  If this file's machinery said the published rows had no minimum either, it would be the")
P("  machinery that is broken and not the witness.")
P("")
P("  %-5s %10s %12s %14s %10s" % ("row", "minima", "alpha_min", "F at edge", "global?"))
for label, content, a_them, mh_them, invR in T1:
    v = [0] * 8
    for rep, e, ep, mult in content:
        v[REPS.index((rep, e, ep))] += mult
    mm, ed, _ = stationary(v, n=120001)
    if mm:
        deep = min(mm, key=lambda t: t[1])
        P("  %-5s %10d %12.6f %14.6e %10s"
          % (label, len(mm), deep[0], ed[1], ed[1] > deep[1]))
    else:
        P("  %-5s %10d %12s %14.6e %10s" % (label, 0, "--", ed[1], "NO MINIMUM"))
P("")
P("  the published rows must have an interior minimum that is also the deepest point of the")
P("  domain -- that is what Part VII's false_vacuum_km.py already measured, and this file has")
P("  to agree with it or be wrong.")
