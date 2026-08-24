#!/usr/bin/env python3
"""analytic_branch_exact.py -- the A4 = 0 branch, on the full potential instead of the expansion.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

THE CHALLENGE, and it is a good one.
The paper states A4 = 0  =>  m_h <= 29.6 GeV as an exact margin, attained at 18 x 7(+,-) with
1/R5 = 181 GeV.  Carles: 1/R5 = 2 m_W / alpha, so 181 GeV means alpha = 2(80.4)/181 = 0.89 -- and
Section closed has just proved the expansion converges only for |alpha| < 1/3, claiming every alpha
used is inside it by a factor of at least 2.7.  So that champion sits FAR outside the radius, and
the truncated identity cannot certify it.  The statement is too strong as written.

He is right, and the repair is cheap because the branch is finite: eleven base solutions and one
free multiplet.  This file redoes them on the EXACT polylogarithmic potential -- global minimiser,
numerical second derivative, m_h -- with no expansion anywhere, and reports what survives.

Two outcomes are possible and both are useful.  If the exact m_h stays far below 125 GeV, the no-go
is stronger than before and rests on nothing that was truncated.  If it does not, the paper loses a
result and says so.

Run:  python analytic_branch_exact.py
"""
import itertools
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

_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
A0, K0 = round(_g["A4"]), round(8 * _g["D"])


def exact(v):
    """global minimiser of the exact F, its m_h and 1/R5.  numeric_min IS the global one."""
    c = cont(v)
    a = numeric_min(c)
    if a is None:
        return None
    h = 2e-4
    f = F(c, np.array([a - h, a, a + h]))
    fpp = float((f[0] - 2 * f[1] + f[2]) / h ** 2)
    mh = KK * math.sqrt(fpp) / a if fpp > 0 else None
    return dict(alpha=a, fpp=fpp, mh=mh, invR=2 * MW / a)


P("=" * 100)
P("0 -- THE OBJECTION, RESTATED AS A NUMBER")
P("=" * 100)
P("  the paper's champion of the branch has 1/R5 = 181 GeV, and 1/R5 = 2 m_W / alpha, so")
P("     alpha = 2 x %.1f / 181 = %.3f" % (MW, 2 * MW / 181.0))
P("  against a radius of convergence 1/c_max = 1/3 = %.3f." % (1.0 / 3))
P("  So the closed form was evaluated a factor %.1f OUTSIDE its own radius." % (2 * MW / 181.0 * 3))
P("  Nothing computed there is certified by the expansion.  Redo it exactly.")

# ----------------------------------------------------------------- the branch, enumerated
usable = [j for j in range(8) if 0 < AV[j] <= -A0]
sols = []


def rec(i, left, vec):
    if left == 0:
        sols.append(tuple(vec))
        return
    if i >= len(usable):
        return
    j = usable[i]
    for m in range(left // AV[j] + 1):
        vec[j] = m
        rec(i + 1, left - m * AV[j], vec)
    vec[j] = 0


rec(0, -A0, [0] * 8)
P("")
P("=" * 100)
P("1 -- THE ELEVEN BASE FAMILIES, ON THE EXACT POTENTIAL")
P("=" * 100)
P("  A4 = 0 forces sum n_j A_j = %d over the multiplets with A4 > 0, and 7(+,+) is free." % (-A0))
P("  base solutions (no 7(+,+)) : %d" % len(sols))
P("")
P("  %-40s %5s %10s %11s %10s %s"
  % ("content", "8D", "alpha", "1/R5 GeV", "m_h GeV", "inside |a|<1/3"))
best = None
for v in sorted(sols, key=lambda t: -sum(t)):
    vv = list(v)
    A4 = A0 + sum(vv[j] * AV[j] for j in range(8))
    k8 = K0 + sum(vv[j] * KV[j] for j in range(8))
    assert A4 == 0
    if k8 <= 0:
        continue
    r = exact(vv)
    if r is None:
        P("  %-40s %5d %10s %11s %10s %s" % (show(vv), k8, "--", "--", "--", "no global minimum"))
        continue
    inside = r["alpha"] < 1.0 / 3
    P("  %-40s %5d %10.5f %11.1f %10s %s"
      % (show(vv), k8, r["alpha"], r["invR"],
         ("%.2f" % r["mh"]) if r["mh"] else "no real", "yes" if inside else "NO"))
    if r["mh"] and (best is None or r["mh"] > best[0]):
        best = (r["mh"], vv, k8, r)

P("")
P("=" * 100)
P("2 -- AND THE FREE MULTIPLET, WHICH THE TRUNCATED ARGUMENT USED TO DISPOSE OF IT")
P("=" * 100)
P("  The truncated argument said: adding 7(+,+) only lowers G, and m_h grows with G, so the")
P("  maximum is at zero copies.  That argument lives entirely inside the expansion.  On the")
P("  exact potential the same sweep has to be done by hand.")
P("")
P("  %-40s %4s %10s %11s %10s" % ("base", "n7pp", "alpha", "1/R5 GeV", "m_h GeV"))
grand = None
for v in sols:
    vv = list(v)
    k8 = K0 + sum(vv[j] * KV[j] for j in range(8))
    for extra in range(0, 30):
        w = list(vv)
        w[0] += extra
        kk = k8 - 6 * extra
        if kk <= 0:
            break
        r = exact(w)
        if r is None or r["mh"] is None:
            continue
        if grand is None or r["mh"] > grand[0]:
            grand = (r["mh"], list(w), kk, r)
if grand:
    mh, w, kk, r = grand
    P("  %-40s %4d %10.5f %11.1f %10.2f" % (show(w), w[0], r["alpha"], r["invR"], mh))
P("")
P("=" * 100)
P("VERDICT")
P("=" * 100)
if grand:
    mh, w, kk, r = grand
    P("  largest m_h anywhere on the analytic branch, EXACT potential : %.2f GeV" % mh)
    P("     at %s   (8D = %d, alpha = %.5f, 1/R5 = %.1f GeV)" % (show(w), kk, r["alpha"], r["invR"]))
    P("     alpha inside the radius 1/3 : %s" % (r["alpha"] < 1.0 / 3))
    P("")
    P("  against the measured 125.25 GeV : short by a factor %.1f" % (125.25 / mh))
    P("  and against the truncated claim of 29.6 GeV : %s"
      % ("the truncated number was optimistic" if mh < 29.6 else
         "the truncated number was pessimistic"))
    P("")
    if mh < MH_LO:
        P("  ==> THE NO-GO SURVIVES, and now it rests on the exact potential rather than on an")
        P("      expansion evaluated outside its radius.  A4 = 0 admits no Standard-Model Higgs.")
    else:
        P("  ==> THE NO-GO DOES NOT SURVIVE.  The paper loses the result and must say so.")
else:
    P("  no content on the branch has a real m_h on the exact potential at all.")
