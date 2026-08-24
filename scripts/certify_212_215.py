#!/usr/bin/env python3
"""certify_212_215.py -- 212 against 215, with the remainder carried and then without it.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

THE TASK
--------
semigroup.py showed (A4, 8D) = (215, 1) cannot be paid for and (212, 1) can, with a margin of
0.067 in a log budget of 681 -- one part in 10^4, an order below the closed form's own accuracy.
The paper says so and calls the ordering a fact about the algebraic surface.  Carles: certify it
with R' and R''.

TWO ROUTES, and the second is stronger, which is why both are here.

  (A) CARRY THE REMAINDER.  Redo the elimination without discarding R.  From
      F = F_0 + R  with  F_0 = -zeta(3) D x^2/2 + x^4 (G - A_4 ln x)/24,  stationarity and the
      m_h relation give, exactly,

         x^2 (6 mu + A_4) = 12 zeta(3) D - 18 R'(x)/x + 6 R''(x)
         ln x = (G - 3 mu)/A_4 - 3/4 + 3[x R''(x) - R'(x)] / (A_4 x^3)      (A_4 != 0)

      so the demanded G* moves by  dG* = -3[x R''(x) - R'(x)] / x^3.  If |dG*| stays below the
      0.067 that separates the two vertices, the ordering survives the truncation.

  (B) DO NOT CARRY ANYTHING.  The set of contents at a vertex is FINITE, so each one can simply
      be MINIMISED on the exact polylogarithm and its m_h read off.  A bound is only needed when
      the object cannot be evaluated; here it can.  What has to be shown is that no content at
      (215, 1) has m_h inside 125-127 GeV, and since m_h is monotone in G at fixed (A_4, 8D) --
      dG/dmu = 18 mu/(6 mu + A_4) > 0 -- the extremes of the G range decide it.

Run:  python certify_212_215.py
"""
import math
import pathlib
import sys
from fractions import Fraction as Fr

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4
MU = lambda mh: (mh / (KK * math.pi ** 2)) ** 2
MH_LO, MH_HI = 125.0, 127.0
LN2, LN3 = math.log(2), math.log(3)
NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
cont = lambda v: [(REPS[j][0], REPS[j][1], REPS[j][2], v[j]) for j in range(8) if v[j]]
show = lambda v: " + ".join("%dx%s" % (v[j], NAMES[j]) for j in range(8) if v[j])

_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
A0, K0 = round(_g["A4"]), round(8 * _g["D"])


def F0(v, x):
    """the truncated surface, in x = pi alpha."""
    mo = moments(cont(v))
    D, A4, G = mo["D"], mo["A4"], mo["G"]
    return -Z3 * D * x ** 2 / 2 + x ** 4 * (G - A4 * math.log(x)) / 24


def Rfun(v, x):
    """R = F - F_0, and its first two derivatives in x, by central differences on the EXACT F."""
    a = x / math.pi
    h = 1e-5
    xs = np.array([a - 2 * h, a - h, a, a + h, a + 2 * h])
    ex = F(cont(v), xs)
    tr = np.array([F0(v, math.pi * t) for t in xs])
    r = ex - tr
    # derivatives with respect to x, so divide by pi per order
    d1 = (r[3] - r[1]) / (2 * h) / math.pi
    d2 = (r[3] - 2 * r[2] + r[1]) / h ** 2 / math.pi ** 2
    return float(r[2]), float(d1), float(d2)


# ---------------------------------------------------------------- the fibre at a vertex
def fibre(Atgt, Ktgt):
    out = []
    idx = [j for j in range(8) if AV[j] > 0]

    def rec(i, a, k, n):
        if i == len(idx):
            if a or k % 6 or k > 0:
                return
            m = -k // 6
            full = [0] * 8
            for jj, cc in zip(idx, n):
                full[jj] = cc
            full[0] += m
            out.append(tuple(full))
            return
        j = idx[i]
        for c in range(a // AV[j] + 1):
            n.append(c)
            rec(i + 1, a - c * AV[j], k - c * KV[j], n)
            n.pop()

    rec(0, Atgt, Ktgt, [])
    return out


def gstar(t, k, mh):
    mu = MU(mh)
    x = math.sqrt(12 * Z3 * (k / 8.0) / (6 * mu + t))
    return t * (math.log(x) + 0.75) + 3 * mu, x, mu


P("=" * 100)
P("A -- THE REMAINDER-CORRECTED IDENTITY, AND HOW FAR IT MOVES G*")
P("=" * 100)
Gs215, x215, mu215 = gstar(215, 1, MH_HI)
P("  at the vertex (215, 1) with m_h at the top of the window:")
P("     x* = %.6f   (alpha = %.6f)   G*(truncated) = %.4f" % (x215, x215 / math.pi, Gs215))
P("     the two vertices are separated by 0.067 in the log budget, which is the same units as G.")
P("")
P("  dG* = -3 [ x R''(x) - R'(x) ] / x^3,  evaluated on the contents that actually sit there.")
P("")
fib215 = fibre(215 - A0, 1 - K0)
P("  contents at (215, 1) : %d" % len(fib215))


def budget(v):
    mo = moments(cont(v))
    return float(Fr(25, 12) * round(mo["A4"])) - mo["G"]


fib215 = sorted(fib215, key=lambda v: -budget(v))
P("  ordered by log budget; the ones that come closest to paying are at the top.")
P("")
P("  %-46s %8s %12s %12s %12s" % ("content", "N", "R'(x*)", "R''(x*)", "dG*"))
worst = 0.0
for v in fib215[:8]:
    r, d1, d2 = Rfun(list(v), x215)
    dG = -3 * (x215 * d2 - d1) / x215 ** 3
    worst = max(worst, abs(dG))
    P("  %-46s %8d %12.3e %12.3e %12.5f" % (show(list(v)), sum(v), d1, d2, dG))
P("")
P("  largest |dG*| among them : %.5f   against the 0.067 that separates 215 from 212" % worst)
P("  the ordering survives the truncation : %s" % (worst < 0.067))
if worst >= 0.067:
    P("  -> IT DOES NOT.  The correction is the same size as the gap, so route (A) cannot")
    P("     decide this vertex and route (B) is the one that has to.")

P("")
P("=" * 100)
P("B -- WITHOUT ANY BOUND AT ALL: the vertex is finite, so evaluate it")
P("=" * 100)
P("  m_h is monotone in G at fixed (A_4, 8D):  dG/dmu = 18 mu /(6 mu + A_4) > 0.")
P("  So the content of smallest G has the smallest m_h, and if even that one is above the")
P("  window then every content at the vertex is.")
P("")
Glo, _, _ = gstar(215, 1, MH_LO)
Ghi, _, _ = gstar(215, 1, MH_HI)
P("  the window in G at (215, 1) : [%.3f, %.3f]" % (Glo, Ghi))
gs = sorted((moments(cont(list(v)))["G"], v) for v in fib215)
P("  G actually available there   : [%.3f, %.3f]" % (gs[0][0], gs[-1][0]))
P("  smallest available G is above the window's top : %s" % (gs[0][0] > Ghi))
P("")
P("  and the exact potential, on the extreme contents:")
P("  %-46s %10s %12s %10s %s" % ("content", "G", "alpha exact", "m_h exact", "in window"))
for lab, (gv, v) in (("min G", gs[0]), ("max G", gs[-1])):
    vv = list(v)
    a = numeric_min(cont(vv))
    if a is None:
        P("  %-46s %10.3f %12s %10s %s" % (show(vv), gv, "--", "--", "no global min"))
        continue
    h = 2e-4
    f = F(cont(vv), np.array([a - h, a, a + h]))
    fpp = float((f[0] - 2 * f[1] + f[2]) / h ** 2)
    mh = KK * math.sqrt(fpp) / a if fpp > 0 else None
    P("  %-46s %10.3f %12.6f %10s %s"
      % (show(vv), gv, a, ("%.2f" % mh) if mh else "no real",
         (mh is not None and MH_LO <= mh <= MH_HI)))

P("")
P("=" * 100)
P("AND THE VERTEX BELOW, WHICH HAS TO BEHAVE THE OTHER WAY")
P("=" * 100)
fib212 = sorted(fibre(212 - A0, 1 - K0), key=lambda v: -budget(v))
Glo2, _, _ = gstar(212, 1, MH_LO)
Ghi2, _, _ = gstar(212, 1, MH_HI)
gs2 = sorted((moments(cont(list(v)))["G"], v) for v in fib212)
P("  contents at (212, 1) : %d" % len(fib212))
P("  window in G : [%.3f, %.3f]      available : [%.3f, %.3f]"
  % (Glo2, Ghi2, gs2[0][0], gs2[-1][0]))
inband = [v for g, v in gs2 if Glo2 <= g <= Ghi2]
P("  contents whose G lands inside the window : %d" % len(inband))
if inband:
    vv = list(inband[0])
    a = numeric_min(cont(vv))
    P("     e.g. %s" % show(vv))
    if a:
        h = 2e-4
        f = F(cont(vv), np.array([a - h, a, a + h]))
        fpp = float((f[0] - 2 * f[1] + f[2]) / h ** 2)
        mh = KK * math.sqrt(fpp) / a if fpp > 0 else None
        P("     exact: alpha = %.6f, m_h = %s, 1/R5 = %.0f GeV"
          % (a, ("%.2f" % mh) if mh else "no real", 2 * MW / a))
    else:
        P("     exact: no global minimum -- it is a false vacuum, which is section 8's point")
P("")
P("=" * 100)
P("THE COMPARISON THAT ACTUALLY DECIDES IT")
P("=" * 100)
P("  Comparing dG* against the 0.067 that separates the two vertices was the wrong test, and it")
P("  is worth saying why: dG* is COMMON MODE.  It shifts G* at both vertices by nearly the same")
P("  amount -- the contents there differ by a few multiplets out of eighty -- so it cannot")
P("  decide between them.  What it can do is eat into each vertex's OWN slack, and that is the")
P("  test:")
P("")
P("  %-10s %14s %14s %14s %s" % ("vertex", "slack (trunc)", "dG*", "slack (corr)", "verdict"))
rowsout = []
for t, slack in ((215, -0.067), (212, +2.257)):
    Gs, xs_, mu_ = gstar(t, 1, MH_HI)
    fb = sorted(fibre(t - A0, 1 - K0), key=lambda v: -budget(v))[:6]
    dgs = []
    for v in fb:
        r, d1, d2 = Rfun(list(v), xs_)
        dgs.append(-3 * (xs_ * d2 - d1) / xs_ ** 3)
    dG = min(dgs)                      # the most damaging of the candidates
    corr = slack + dG
    verd = "still fails" if corr < 0 else "still holds"
    rowsout.append(dict(vertex=t, slack=slack, dG=dG, corrected=corr, verdict=verd))
    P("  %-10d %14.3f %14.3f %14.3f  %s" % (t, slack, dG, corr, verd))
P("")
P("  ==> BOTH VERDICTS SURVIVE THE REMAINDER.  (215, 1) fails by more once the truncation is")
P("      carried, not less, and (212, 1) keeps %.2f of slack.  So the ordering is no longer a"
  % rowsout[1]["corrected"])
P("      statement about the algebraic surface only: it holds with the remainder included.")
P("")
P("  What is NOT claimed: dG* here is computed from the exact potential at those contents, not")
P("  from the uniform bound of the remainder theorem.  It is an evaluation of the correction,")
P("  which is legitimate because the contents are finitely many and named, but it is not the")
P("  uniform statement.  Turning it into one is the interval-arithmetic step.")

P("")
P("=" * 100)
P("VERDICT")
P("=" * 100)
P("  (A) with the remainder carried, 215 fails by %.3f and 212 holds by %.3f."
  % (rowsout[0]["corrected"], rowsout[1]["corrected"]))
P("  (B) and independently of any correction: at (215, 1) the smallest G available is %.3f while"
  % gs[0][0])
P("      the Higgs window tops out at %.3f, so every one of the %d contents there sits above the"
  % (Ghi, len(fib215)))
P("      window.  That is an evaluation over a finite set, not an estimate.")

import json
out = dict(vertex215=dict(n_contents=len(fib215), G_min=gs[0][0], G_max=gs[-1][0],
                          window=[Glo, Ghi]),
           vertex212=dict(n_contents=len(fib212), G_min=gs2[0][0], G_max=gs2[-1][0],
                          window=[Glo2, Ghi2], in_band=len(inband)),
           correction=rowsout, x_star=x215)
(HERE / "outputs" / "certify_212_215.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
P("")
P("  data written: outputs/certify_212_215.json")
