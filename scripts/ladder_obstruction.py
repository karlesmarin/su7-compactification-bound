#!/usr/bin/env python3
"""The obstruction is a property of the whole LADDER, not of the curvature.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Reading Cacciapaglia et al. 2409.16137 forced this file.  Their eq. (2.13) is

    F_+(0) = (3/2) zeta(5) ,   F_-(0) = -(15/16) F_+(0) ,

so their vacuum-stability criterion at the symmetric point is  n_+ - (15/16) n_- : the SAME shape
as our  D = A2 - (3/4) B2, one rung down the ladder.  Both are the single expansion

    Re Li_5(s e^{i c pi a}) = sum_n s^n cos(n c pi a)/n^5
                            = sum_k (-1)^k (c pi a)^{2k}/(2k)! * sum_n s^n / n^{5-2k}

so the 2k-th derivative in alpha is governed by  p = 5 - 2k  and by the 2k-th moment of the content.
With  sum_n s^n/n^p = 2^-p zeta(p) + s (1 - 2^-p) zeta(p),  the rung-k coefficient is

    [ S_k + (2^p - 1) Delta_k ] / 2^p ,    S_k = sum m c^{2k} ,  Delta_k = sum m c^{2k} s .

  k=0  p=5  ->  (S_0 + 31 Delta_0)/32   the VALUE          -- Cacciapaglia et al.'s criterion
  k=1  p=3  ->  (S_1 +  7 Delta_1)/8    the CURVATURE      -- our D, the odd eighths
  k=2  p=1  ->  the weight degenerates  -- the two halves diverge and leave the LOGARITHM

What is measured here:
  1  that our rung k=0 reproduces their (2.13) identically, on our own basis;
  2  that the gauge sector's contribution to rung k is  ODD  at k=0 and k=1 and EVEN at k=2,
     and that the reason is one number: the half-weight 1/2 of the P_6 = -1 pairs, which enters
     the rung as  7 (2^{p-1} - 1);
  3  so the odd-eighths theorem is NOT special to the curvature: the same obstruction sits on the
     VALUE, and Komori-Maru's model can never be marginal in Cacciapaglia et al.'s OWN sense either;
  4  and the contrast that keeps the two statements distinct: with pure component COUNTS -- their
     setting, all multiplicities integer -- marginality IS reachable at k=0.  It is the half-weight
     that forbids it, and the half-weight is theirs, derived from their eqs. (63)-(67).
"""
import itertools
import json
import math
import pathlib
from fractions import Fraction as Fr

import numpy as np

exec(open(pathlib.Path(__file__).resolve().parent / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
Z5, ETA5 = 1.0369277551433699263, 0.97211977044690930594      # zeta(5), eta(5)
RUNGS = [(0, 5), (1, 3), (2, 1)]


def rung(content, k):
    """(S_k, Delta_k, numerator S_k + (2^p-1) Delta_k) as exact Fractions, p = 5-2k."""
    p = 5 - 2 * k
    S = D = Fr(0)
    for m, s, c in table(content):
        w = Fr(m).limit_denominator(64) * Fr(int(c)) ** (2 * k)
        S += w
        D += w * s
    return S, D, S + (2 ** p - 1) * D


P("")
P("=" * 100)
P("1 -- OUR RUNG k=0 IS THEIR EQ. (2.13).  checked against our own basis at alpha = 0")
P("=" * 100)
P("Their (2.13):  F_-(0) = -(15/16) F_+(0).  In our variables  F(0) = zeta(5) (S_0 + 31 Delta_0)/32,")
P("because n_+ = (S_0+Delta_0)/2, n_- = (S_0-Delta_0)/2 gives  n_+ - (15/16) n_- = (S_0+31 Delta_0)/32.")
P("")
P("   eta(5)/zeta(5) = %.15f   and   15/16 = %.15f   ->  %s" %
  (ETA5 / Z5, 15 / 16, "SAME" if abs(ETA5 / Z5 - 15 / 16) < 1e-15 else "*** DIFFER ***"))
P("")
P("%-5s %14s %14s %16s %16s %10s" % ("row", "S_0", "Delta_0", "(S_0+31D_0)/32", "F(0) numeric", "err"))
for label, content, a_them, mh_them, invR in T1:
    S0, D0, num = rung(content, 0)
    pred = float(num) / 32 * Z5
    got = float(F(content, np.array([0.0]))[0])
    P("%-5s %14s %14s %16.8f %16.8f %10.1e" %
      (label, S0, D0, pred, got, abs(pred - got)))

P("")
P("=" * 100)
P("2 -- THE OBSTRUCTION ALONG THE LADDER.  gauge sector alone, rung by rung")
P("=" * 100)
P("%-6s %5s %12s %12s %18s %8s" % ("rung k", "p", "S_k", "Delta_k", "S_k+(2^p-1)D_k", "parity"))
gauge_num = {}
for k, p in RUNGS:
    S, D, num = rung([], k)
    gauge_num[k] = num
    P("%-6d %5d %12s %12s %18s %8s" %
      (k, p, S, D, num, "ODD" if (num.denominator == 1 and num.numerator % 2) else "even"))
P("")
P("The half-integer part is carried by the single term (m,s,c) = (-7/2,-1,1) -- the P_6 = -1 pairs,")
P("whose weight 1/2 is DERIVED from their eqs. (63)-(67), not assumed (see MARGINALITY.md).  It")
P("enters rung k as   -7/2 + (2^p-1)(+7/2)  =  (7/2)(2^p - 2)  =  7 (2^{p-1} - 1):")
for k, p in RUNGS:
    P("   k=%d, p=%d:  7 (2^{p-1} - 1) = %3d   -> %s" %
      (k, p, 7 * (2 ** (p - 1) - 1), "ODD, obstruction lives" if (7 * (2 ** (p - 1) - 1)) % 2
       else "EVEN, obstruction dies"))
P("")
P("So the obstruction lives at moments 0 and 2 and dies at 4 -- and the reason is ONE number.")

P("")
P("=" * 100)
P("3 -- MATTER CONTRIBUTES EVEN AT EVERY RUNG, so the parity is the gauge sector's alone")
P("=" * 100)
P("   per term, m and c integers:  m c^{2k} [1 + (2^p-1) s]  =  m c^{2k} 2^p   (s=+1)")
P("                                                          or  m c^{2k}(2-2^p) (s=-1),  both even.")
P("%-10s %14s %14s %14s" % ("multiplet", "k=0", "k=1", "k=2"))
bad = []
for j, (r, e, ep) in enumerate(REPS):
    row = []
    for k, p in RUNGS:
        _, _, num = rung([(r, e, ep, 1)], k)
        n = num - gauge_num[k]
        row.append(n)
        if n.denominator != 1 or n.numerator % 2:
            bad.append((NAMES[j], k))
    P("%-10s %14s %14s %14s" % (NAMES[j], row[0], row[1], row[2]))
P("   multiplets contributing an ODD amount at any rung: %s" % (bad or "NONE"))

P("")
P("=" * 100)
P("4 -- THEREFORE: the VALUE can never be marginal either, in Komori-Maru's model")
P("=" * 100)
P("Sweeping every content of at most 6 multiplets:")
zero_val, zero_curv, n = 0, 0, 0
minabs = None
for mults in itertools.product(range(7), repeat=8):
    if not 1 <= sum(mults) <= 6:
        continue
    content = [(r, e, ep, m) for (r, e, ep), m in zip(REPS, mults) if m]
    _, _, v0 = rung(content, 0)
    _, _, v1 = rung(content, 1)
    n += 1
    zero_val += (v0 == 0)
    zero_curv += (v1 == 0)
    if minabs is None or abs(v0) < minabs[0]:
        minabs = (abs(v0), content)
P("   contents swept                              : %d" % n)
P("   with a marginal VALUE      (S_0+31 D_0 = 0) : %d" % zero_val)
P("   with a marginal CURVATURE  (S_1+ 7 D_1 = 0) : %d" % zero_curv)
P("   smallest |S_0+31 Delta_0| found             : %s" % minabs[0])
P("")
P("   32 F(0)/zeta(5) = 9 + even  ->  ALWAYS ODD  ->  never zero, and |F(0)| >= zeta(5)/32.")
P("   Their own stability criterion, applied to a model they did not consider, is always decided.")

P("")
P("=" * 100)
P("5 -- THE CONTRAST THAT KEEPS THE TWO STATEMENTS APART")
P("=" * 100)
P("In THEIR setting the multiplicities are component COUNTS -- all integers, no half-weights.  Then")
P("the criterion  n_+ - (15/16) n_-  vanishes iff  16 n_+ = 15 n_- , which has solutions:")
for m_ in (1, 2, 3):
    P("   n_+ = %2d, n_- = %2d  ->  16 n_+ - 15 n_- = %d" % (15 * m_, 16 * m_, 16 * 15 * m_ - 15 * 16 * m_))
P("")
P("So marginality at the value IS reachable with pure counts.  What forbids it here is the")
P("half-weight 1/2 of the P_6 = -1 gauge pairs -- and that half-weight is a 6D effect: it exists")
P("only because there is a second parity.  The protection is the same one MARGINALITY.md found for")
P("the curvature, and it is now visible one rung further up.")

OUT.mkdir(exist_ok=True)
(OUT / "ladder_obstruction.json").write_text(json.dumps(dict(
    gauge_numerator={str(k): str(v) for k, v in gauge_num.items()},
    half_weight_term={str(p): 7 * (2 ** (p - 1) - 1) for _, p in RUNGS},
    swept=n, marginal_value=zero_val, marginal_curvature=zero_curv,
    min_abs_value_numerator=str(minabs[0])), indent=1), encoding="utf-8")
P("")
P("written: %s" % (OUT / "ladder_obstruction.json"))
