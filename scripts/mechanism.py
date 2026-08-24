#!/usr/bin/env python3
"""THE MECHANISM: even windings are bulk, odd windings are fixed-point.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Every term of the one-loop Wilson-line potential is  m c^(2k) sum_n s^n / n^p .
Split the WINDING sum by the parity of n.  Since s^n = +1 for even n whatever s, and
s^n = s for odd n,

    sum_n s^n / n^p  =  sum_{n even} n^-p  +  s sum_{n odd} n^-p
                     =  2^-p zeta(p)  +  s (1 - 2^-p) zeta(p)                    (*)

Summing over the content with weights m c^(2k) turns the two brackets into the two
group-theoretic invariants at once:

    even windings  ->  S     = sum m c^(2k)      = 2 T(R)          (parity-blind)
    odd  windings  ->  Delta = sum m c^(2k) s    = 2 Tr[P Q^2]     (parity-weighted)

so, at the curvature (p = 3),

    V''(0) = -pi^2 [ 2^-3 zeta(3) S + (1 - 2^-3) zeta(3) Delta ] = -pi^2 zeta(3) (S + 7 Delta)/8 .

THE EIGHTHS ARE 2^3, FROM THE CUBE IN 1/n^3.  THE 7 IS 2^3 - 1.  And the whole ladder is
(*) at p = 5-2k: denominators 2^(5-2k) = 32, 8, 2, and w_2k = eta/zeta = 1 - 2^(1-p).

PHYSICS.  In the image-sum picture an even winding closes on the covering circle and never
meets a fixed point -- it is a BULK path; an odd winding does -- it is a FIXED-POINT path.
So the bulk/brane split of the group theory IS the even/odd split of the winding sum, and
the fixed-point half carries 7/8 of the curvature against 1/8 for the bulk.
"""
import math
import pathlib
import re
from fractions import Fraction

from mpmath import inf, mp, mpf, nsum, zeta

mp.dps = 30
P = lambda *a: print(*a, flush=True)
SRC = pathlib.Path(__file__).resolve().parent.parent / "part_vi" / "su7_anchor_mh.py"
_ns = {}
_s = SRC.read_text(encoding="utf-8")
for pat in (r"\ndef terms\(.*?\n(?=\n\ndef )", r"\nGAUGE = \[.*?\]\s*#[^\n]*\n", r"\nT1 = \[.*?\)\]\n"):
    exec(re.search(pat, _s, re.S).group(0), {}, _ns)
terms, GAUGE, T1 = _ns["terms"], _ns["GAUGE"], _ns["T1"]

P("=" * 90)
P("1 -- the identity (*), checked numerically")
P("=" * 90)
for p in (5, 3):
    for s in (1, -1):
        # nsum, not a truncated partial sum: the tail of sum 1/n^3 at N=4e4 is 3e-10 and
        # would read as a failure of the identity rather than of the summation.
        direct = nsum(lambda n: mpf(s) ** n / n ** p, [1, inf])
        split = mpf(2) ** (-p) * zeta(p) + s * (1 - mpf(2) ** (-p)) * zeta(p)
        P("  p=%d s=%+d :  sum = %.20f   2^-p z + s(1-2^-p) z = %.20f   ok=%s" %
          (p, s, direct, split, abs(direct - split) < mpf('1e-20')))

P("")
P("=" * 90)
P("2 -- even windings = BULK share, odd = FIXED-POINT share, at each moment")
P("=" * 90)
P("%-8s %-6s %14s %14s %10s" % ("moment", "p", "even (bulk)", "odd (brane)", "denom"))
for k, p in ((0, 5), (1, 3), (2, 1)):
    ev, od = Fraction(1, 2 ** p), 1 - Fraction(1, 2 ** p)
    P("%-8d %-6d %14s %14s %10d" % (2 * k, p, ev, od, 2 ** p))
P("  -> the eighths of Part VI are 2^3; the 7 is 2^3-1; the ladder's 32, 8, 2 are 2^(5-2k).")
P("  -> at moment 4 (p=1) both halves diverge and only their DIFFERENCE survives: the log.")

P("")
P("=" * 90)
P("3 -- so how much of the EWSB curvature is a fixed-point effect?")
P("=" * 90)
P("%-6s %10s %10s %12s %12s %8s" % ("row", "S", "Delta", "bulk 1/8 S", "brane 7/8 D", "brane %"))
for label, content, a_them, mh, invR in T1:
    tab = list(GAUGE)
    for rep, e, ep, mult in content:
        for m, s, c in terms(rep, e, ep):
            tab.append((m * mult, s, c))
    S = sum(m * c * c for m, s, c in tab)
    Dl = sum(m * c * c * s for m, s, c in tab)
    b, br = Fraction(S) / 8, Fraction(7 * Dl) / 8
    tot = b + br
    P("%-6s %10s %10s %12s %12s %7.1f%%" %
      (label, Fraction(S).limit_denominator(4), Fraction(Dl).limit_denominator(4),
       b, br, 100 * abs(br) / (abs(b) + abs(br))))
P("")
P("  the two halves have OPPOSITE signs on every row: the bulk term is positive and")
P("  destabilising, the fixed-point term negative -- and D is what is left over.")
