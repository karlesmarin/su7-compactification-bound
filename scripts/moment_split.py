#!/usr/bin/env python3
"""D = (S + 7 Delta)/8 -- the minimal form, and the clean proof that D is an odd eighth.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

D = A2 - (3/4) B2 with A2, B2 the second moments of the Wilson-line charge over the
periodic and antiperiodic sectors.  Change to the natural basis

    S     = A2 + B2 = sum_terms m c^2        the plain second moment
    Delta = A2 - B2 = sum_terms m c^2 s      the PARITY-SIGNED second moment

then, with w = eta(3)/zeta(3) = 3/4,

    D = 1/2 [ (1-w) S + (1+w) Delta ]  =  (S + 7 Delta)/8 ,

so the 7 is 4+3 out of the zeta ratio, not out of SU(7).

AND THE THEOREM FALLS OUT.  For matter, c and m are integers, so S and Delta are
integers with S + Delta = 2 A2 EVEN; hence S + 7 Delta = (S+Delta) + 6 Delta is EVEN.
The gauge sector contributes S + 7 Delta = -27, ODD.  Therefore the total is always odd:

    D = odd / 8  !=  0 ,  for ANY bulk content.
"""
import math
import pathlib
import re
from fractions import Fraction

P = lambda *a: print(*a, flush=True)
SRC = pathlib.Path(__file__).resolve().parent.parent / "part_vi" / "su7_anchor_mh.py"
_ns = {}
_src = SRC.read_text(encoding="utf-8")
for pat in (r"\ndef terms\(.*?\n(?=\n\ndef )", r"\nGAUGE = \[.*?\]\s*#[^\n]*\n", r"\nT1 = \[.*?\)\]\n"):
    exec(re.search(pat, _src, re.S).group(0), {}, _ns)
terms, GAUGE, T1 = _ns["terms"], _ns["GAUGE"], _ns["T1"]


def SD(tab):
    S = sum(m * c * c for m, s, c in tab)
    Dl = sum(m * c * c * s for m, s, c in tab)
    return S, Dl


P("=" * 86)
P("D = (S + 7 Delta)/8 on every multiplet and on the gauge sector")
P("=" * 86)
P("%-12s %8s %8s %10s %10s %6s" % ("object", "S", "Delta", "(S+7D)/8", "D direct", "ok"))
bad = 0
rowsx = [("gauge", GAUGE)]
for rep in ("7", "28", "48", "84"):
    for (e, ep), tag in (((1, 1), "(+,+)"), ((1, -1), "(+,-)")):
        rowsx.append((rep + tag, terms(rep, e, ep)))
for name, tab in rowsx:
    S, Dl = SD(tab)
    direct = sum(m * c * c for m, s, c in tab if s > 0) - 0.75 * sum(m * c * c for m, s, c in tab if s < 0)
    viaSD = (S + 7 * Dl) / 8
    ok = abs(viaSD - direct) < 1e-12
    bad += (not ok)
    P("%-12s %8s %8s %10s %10s %6s" %
      (name, Fraction(S).limit_denominator(8), Fraction(Dl).limit_denominator(8),
       Fraction(viaSD).limit_denominator(64), Fraction(direct).limit_denominator(64), ok))
P("")
P("mismatches: %d   %s" % (bad, "PASS" if bad == 0 else "FAILED"))

P("")
P("=" * 86)
P("THE PARITY ARGUMENT, checked rather than asserted")
P("=" * 86)
Sg, Dg = SD(GAUGE)
P("  gauge:  S = %s , Delta = %s , S + 7 Delta = %s  -> %s" %
  (Fraction(Sg).limit_denominator(8), Fraction(Dg).limit_denominator(8),
   Fraction(Sg + 7 * Dg).limit_denominator(8),
   "ODD" if abs((Sg + 7 * Dg) % 2) == 1 else "EVEN"))
P("  every matter multiplet must give S + 7 Delta EVEN:")
allev = True
for name, tab in rowsx[1:]:
    S, Dl = SD(tab)
    v = S + 7 * Dl
    ev = (v % 2 == 0)
    allev &= ev
    P("     %-10s S+Delta = %3d (even: %s)   S+7Delta = %4d  %s" %
      (name, S + Dl, (S + Dl) % 2 == 0, v, "EVEN" if ev else "ODD  <-- would break it"))
P("  all even: %s" % allev)
P("")
P("  => total S + 7 Delta = odd + even = ODD, so D = odd/8 and NEVER zero.")
P("")
P("  the five published rows, as a check:")
for label, content, a_them, mh, invR in T1:
    tab = list(GAUGE)
    for rep, e, ep, mult in content:
        for m, s, c in terms(rep, e, ep):
            tab.append((m * mult, s, c))
    S, Dl = SD(tab)
    P("     %-5s S = %7s  Delta = %7s   S+7Delta = %5s (odd: %s)   D = %s" %
      (label, Fraction(S).limit_denominator(8), Fraction(Dl).limit_denominator(8),
       Fraction(S + 7 * Dl).limit_denominator(8), abs((S + 7 * Dl) % 2) == 1,
       Fraction((S + 7 * Dl) / 8).limit_denominator(64)))
