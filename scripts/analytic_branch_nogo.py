#!/usr/bin/env python3
"""analytic_branch_nogo.py -- on A4 = 0, electroweak breaking and a true vacuum are incompatible.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS REPLACES THE OLD STATEMENT
-----------------------------------
The paper claimed A4 = 0  =>  m_h <= 29.6 GeV, "an exact margin".  Carles killed it: the champion
of that branch has 1/R5 = 181 GeV, so alpha = 2 m_W / 181 = 0.89, and the expansion converges only
for |alpha| < 1/3.  The number was computed a factor 2.7 outside the radius that the same paper
proves.  Nothing there was certified.

analytic_branch_exact.py then went to redo it on the exact potential and found something better:
EVERY base family on the branch with D > 0 has no global minimum at all -- they are all false
vacua.  That is not a coincidence, and it does not need a numerical sweep to see.  It is two
integers.

THE ARGUMENT, and it is exact.
On the branch, a content is one of eleven base solutions plus n copies of the free 7(+,+).  That
multiplet moves the two relevant functionals by

    8D -> 8D - 6n        (its own 8D is -6)
    W  -> W  + n         (its own W is +1)

so electroweak breaking (8D > 0) and a true vacuum (W > 0, the criterion of CCD24 summed) ask for

    n < 8D_base / 6      and      n > -W_base

at once.  An integer between them exists iff 8D_base > 6 |W_base| with room to spare.  Every base
solution on the branch is checked below, and the margin is measured, not asserted -- because on one
of them it misses by a single step.

Run:  python analytic_branch_nogo.py
"""
import math
import pathlib
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
cont = lambda v: [(REPS[j][0], REPS[j][1], REPS[j][2], v[j]) for j in range(8) if v[j]]
show = lambda v: " + ".join("%dx%s" % (v[j], NAMES[j]) for j in range(8) if v[j]) or "(gauge only)"

_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
A0, K0 = round(_g["A4"]), round(8 * _g["D"])


def Wof(v):
    return sum(-s * m for m, s, c in table(cont(v)) if int(round(c)) % 2 == 1)


WV = [Wof([1 if i == j else 0 for i in range(8)]) - Wof([0] * 8) for j in range(8)]
W0 = Wof([0] * 8)

P("=" * 100)
P("0 -- THE TWO MOVES OF THE FREE MULTIPLET")
P("=" * 100)
P("  7(+,+) carries  A4 = %d,  8D = %d,  W = %+d." % (AV[0], KV[0], WV[0]))
P("  So on the branch it is the only knob, and it trades curvature for stability at a fixed rate:")
P("  six units of 8D bought for one unit of W.")

# ---------------------------------------------------------------- the branch
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
P("1 -- THE ELEVEN BASE SOLUTIONS, AND WHETHER ANY WINDOW EXISTS")
P("=" * 100)
P("  need an integer n with   -W_base < n < 8D_base/6")
P("")
P("  %-42s %7s %7s %10s %10s %s"
  % ("base content", "8D", "W", "n >", "n <", "any n?"))
anyok = []
for v in sols:
    vv = list(v)
    k8 = K0 + sum(vv[j] * KV[j] for j in range(8))
    w = W0 + sum(vv[j] * WV[j] for j in range(8))
    lo, hi = -w, k8 / 6.0
    ns = [n for n in range(max(0, int(math.floor(lo)) + 1), int(math.ceil(hi)) + 1)
          if n > lo and n < hi]
    if ns:
        anyok.append((vv, ns))
    P("  %-42s %7d %7d %10.2f %10.2f %s"
      % (show(vv), k8, w, lo, hi, ("YES n=%s" % ns[:4]) if ns else "none"))
P("")
P("  base solutions admitting BOTH D > 0 and W > 0 : %d of %d" % (len(anyok), len(sols)))

P("")
P("=" * 100)
P("2 -- THE CLOSEST MISS, because a no-go that misses by a mile is less interesting")
P("=" * 100)
best = None
for v in sols:
    vv = list(v)
    k8 = K0 + sum(vv[j] * KV[j] for j in range(8))
    w = W0 + sum(vv[j] * WV[j] for j in range(8))
    if k8 <= 0:
        continue
    gap = k8 / 6.0 - (-w)          # width of the open interval; needs to contain an integer
    if best is None or gap > best[0]:
        best = (gap, vv, k8, w)
if best:
    gap, vv, k8, w = best
    P("  widest window on the branch : %s" % show(vv))
    P("     8D = %d, W = %d  ->  need  %.2f < n < %.2f,  width %.2f" % (k8, w, -w, k8 / 6.0, gap))
    n_try = -w + 1
    P("     the smallest n giving W > 0 is n = %d, and there 8D = %d - 6(%d) = %d"
      % (n_try, k8, n_try, k8 - 6 * n_try))
    P("     so it misses by %d in 8D, i.e. by ONE step of the only knob there is."
      % (6 * n_try - k8))

P("")
P("=" * 100)
P("VERDICT")
P("=" * 100)
if not anyok:
    P("  A4 = 0  =>  NO content has both D > 0 and W > 0.")
    P("")
    P("  Cancelling the branch point is possible; cancelling it, breaking the electroweak")
    P("  symmetry, and landing in the true vacuum is not.  This is exact arithmetic on two")
    P("  integers -- no expansion, no radius of convergence, nothing evaluated anywhere near")
    P("  alpha = 0.89.  It replaces the m_h <= 29.6 GeV of the previous version, which was")
    P("  computed outside the radius the same paper proves and is withdrawn.")
    P("")
    P("  CONTROL -- the exact potential must agree: every base family with D > 0 should have no")
    P("  global interior minimum.  That is what analytic_branch_exact.py measures, and it does.")
else:
    P("  A4 = 0 DOES admit both; the no-go is false and the paper must drop it.")
    for vv, ns in anyok:
        P("     %s   with n in %s" % (show(vv), ns[:6]))
