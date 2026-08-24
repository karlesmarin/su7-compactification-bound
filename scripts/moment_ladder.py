#!/usr/bin/env python3
"""The moment ladder: w_{2k} = eta(5-2k)/zeta(5-2k), and where the obstruction dies.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Every even Taylor coefficient of F = sum m Re Li_5(s e^{i c pi a}) at a = 0 pairs a
periodic zeta with an antiperiodic eta, two arguments lower per derivative:

    moment 2k   ->   zeta(5-2k)  vs  eta(5-2k) ,    w_2k = eta/zeta = 1 - 2^(2k-4)

    k=0  w = 15/16      k=1  w = 3/4       k=2  w = 0        k=3  w = -3
                         (the curvature D)   (zeta(1) diverges)  (sign flips)

With S = sum m c^(2k) and Delta = sum m c^(2k) s,

    D_2k = 1/2[(1-w)S + (1+w)Delta] = [ S + (2^(5-2k) - 1) Delta ] / 2^(5-2k) .

k=1 gives Part VI's D = (S + 7 Delta)/8.  k=2 gives (S+Delta)/2 = A4, PURELY PERIODIC --
which is why the alpha^4 ln alpha coefficient has no antiperiodic part: not an accident,
but eta(1)/zeta(1) = 0.

THE PARITY ARGUMENT CLIMBS TOO.  Matter always has S + Delta even, and
S + (2^m - 1) Delta = (S+Delta) + (2^m - 2) Delta is even for every m >= 1.  So at EVERY
rung the matter contribution is even, and whether the rung can vanish is decided by the
GAUGE SECTOR ALONE.
"""
import itertools
import math
import pathlib
import re
from fractions import Fraction

P = lambda *a: print(*a, flush=True)
SRC = pathlib.Path(__file__).resolve().parent.parent / "part_vi" / "su7_anchor_mh.py"
_ns, _src = {}, (pathlib.Path(SRC)).read_text(encoding="utf-8")
for pat in (r"\ndef terms\(.*?\n(?=\n\ndef )", r"\nGAUGE = \[.*?\]\s*#[^\n]*\n"):
    exec(re.search(pat, _src, re.S).group(0), {}, _ns)
terms, GAUGE = _ns["terms"], _ns["GAUGE"]
REPS = [(r, e, ep) for r in ("7", "28", "48", "84") for (e, ep) in ((1, 1), (1, -1))]


def SD(tab, k):
    S = sum(m * c ** (2 * k) for m, s, c in tab)
    Dl = sum(m * c ** (2 * k) * s for m, s, c in tab)
    return S, Dl


P("=" * 88)
P("THE LADDER -- w_2k = eta(5-2k)/zeta(5-2k) = 1 - 2^(2k-4), and the gauge sector's parity")
P("=" * 88)
P("%-4s %10s %14s %26s %10s" % ("2k", "w_2k", "2^(5-2k) D", "gauge value", "parity"))
gauge_par = {}
for k in (0, 1, 2):
    w = Fraction(1) - Fraction(2) ** (2 * k - 4)
    S, Dl = SD(GAUGE, k)
    scale = 2 ** (5 - 2 * k)
    val = S + (scale - 1) * Dl
    gauge_par[k] = val
    P("%-4d %10s %14s %26s %10s" %
      (2 * k, w, "S + %d*Delta" % (scale - 1), Fraction(val).limit_denominator(8),
       "ODD" if abs(val % 2) == 1 else "EVEN"))

P("")
P("  matter must give an EVEN value at every rung -- checked:")
allok = True
for k in (0, 1, 2):
    scale = 2 ** (5 - 2 * k)
    for (r, e, ep) in REPS:
        S, Dl = SD(terms(r, e, ep), k)
        v = S + (scale - 1) * Dl
        if v % 2 != 0:
            allok = False
            P("     rung %d, %s%s : %s  <-- ODD, would break it" % (2 * k, r, (e, ep), v))
P("     all even at rungs 0, 2, 4 : %s" % allok)

P("")
P("=" * 88)
P("CONSEQUENCE -- the obstruction lives at moments 0 and 2 and DIES at moment 4")
P("=" * 88)
for k in (0, 1, 2):
    P("  moment %d : gauge = %s (%s)  ->  %s" %
      (2 * k, Fraction(gauge_par[k]).limit_denominator(8),
       "odd" if abs(gauge_par[k] % 2) == 1 else "even",
       "can NEVER vanish" if abs(gauge_par[k] % 2) == 1 else "CAN vanish for suitable matter"))

P("")
P("  moment 4 is A4, the strength of the alpha^4 ln alpha branch point.")
P("  A4 = 0 means the non-analyticity at the symmetric point switches OFF.")
P("  gauge A4 = %s ; is there matter making the total vanish, WITH D > 0?" %
  Fraction(sum(m * c ** 4 for m, s, c in GAUGE if s > 0)).limit_denominator(8))

# CORRECTED 2026-08-22.  This loop used to run over range(5), i.e. multiplicity at most FOUR per
# multiplet, while the sentence it feeds in the paper says 'contents of at most six multiplets'.
# The cap was silent and it cost exactly one content: 5x7(+,+) + 1x48(+,+), which has A4 = 0 and
# 8D = -57, so it does not touch the 'five with D > 0' -- but the total is SEVENTEEN, not sixteen.
# Ver [[an-aggregate-count-is-not-a-case]] and the no-silent-caps rule.
hits = []
for mults in itertools.product(range(7), repeat=len(REPS)):
    if not 1 <= sum(mults) <= 6:
        continue
    tab = list(GAUGE)
    for (r, e, ep), mu in zip(REPS, mults):
        if mu:
            for m, s, c in terms(r, e, ep):
                tab.append((m * mu, s, c))
    A4 = sum(m * c ** 4 for m, s, c in tab if s > 0)
    S2, D2 = SD(tab, 1)
    D = (S2 + 7 * D2) / 8
    if abs(A4) < 1e-9:
        hits.append((D, mults, sum(mults)))
P("")
P("  contents (<=6 multiplets) with A4 EXACTLY 0 : %d" % len(hits))
pos = [h for h in hits if h[0] > 0]
P("  of those with D > 0 (electroweak breaking) : %d" % len(pos))
for D, mults, n in hits[:6]:
    lab = "+".join("%dx%s%s" % (mu, r, "(+,+)" if ep > 0 else "(+,-)")
                   for (r, e, ep), mu in zip(REPS, mults) if mu)
    P("     D = %-8s  %s" % (Fraction(D).limit_denominator(64), lab))
