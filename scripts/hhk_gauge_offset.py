#!/usr/bin/env python3
"""hhk_gauge_offset.py -- their gauge offset is an EVEN INTEGER, ours is a half-integer.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHY THIS EXISTS.  Haba, Hosotani and Kawamura (hep-ph/0309088, PTP 111 (2004) 265) write the
one-loop orbifold vacuum energy as a linear form in the bulk matter multiplicities and reduce its
whole boundary-condition dependence to two integer functionals -- their eqs. (4.13)-(4.14),

    V_eff = (p,q,r,s-independent terms) + 4 w(beta) h(q+r),   h(x) = a x(N-x) - b x,
    a = 2 - 2 n_Ad^(+) + 2 n_Ad^(-) - n_A^(+) + n_A^(-),      b = n_F^(+) - n_F^(-).

Their b is our W in another notation: h(N) - h(0) = -N b, so a signed count of matter by parity
controls the difference between two symmetric points.  The architecture is theirs.

WHAT IS NOT THEIRS is the parity.  Their gauge sector contributes the EVEN INTEGER +2 to a and
NOTHING to b, so both functionals can vanish -- and they open their classification with "(i) The
case with a = 0" and list "n_F^(+) = n_F^(-) => completely degenerate" as an ordinary case.  In the
model of this paper the gauge sector contributes -7/2, a HALF-integer, so 8D, 2U and 2W are odd for
every content and the degenerate case they enumerate is arithmetically unreachable.  One number is
the whole difference, and this file pins it down.

AND CHECKING IT TURNED UP A DROPPED FACTOR OF TWO IN THEIR EQ. (4.7), which is worth recording
because (4.7) is exactly where that offset is fixed.  It is inert -- their (4.13) uses the correct
value -- but a reader checking our claim against (4.7) as printed would be misled.

Sections:
  1  their (3.20) rebuilt from scratch, counting adjoint entries by parity block
  2  the sum that (4.7) substitutes, which is 2(p+s)(q+r) and not (p+s)(q+r)
  3  their (4.5) + (4.6) summed, which reproduces the FIRST line of (4.7) exactly
  4  the verdict, and their own (4.13) as the arbiter

Run:  python hhk_gauge_offset.py
"""
import pathlib
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

# blocks p, q, r, s carry the (P0,P1) signs below; the adjoint entry (i,j) has parity
# (sigma0_a sigma0_b, sigma1_a sigma1_b) for i in block a and j in block b.  Nothing here is
# taken from their paper: it is the definition of the orbifold parity acting on matrices.
SIGN = {0: (+1, +1), 1: (+1, -1), 2: (-1, +1), 3: (-1, -1)}
CONFIGS = [(1, 1, 0, 0), (2, 1, 1, 0), (1, 2, 1, 1), (3, 1, 1, 2), (2, 2, 2, 2), (4, 0, 3, 1)]


def count_adjoint(sizes):
    blocks = []
    for b, n in enumerate(sizes):
        blocks += [b] * n
    cnt = {}
    for i in blocks:
        for j in blocks:
            par = (SIGN[i][0] * SIGN[j][0], SIGN[i][1] * SIGN[j][1])
            cnt[par] = cnt.get(par, 0) + 1
    cnt[(1, 1)] -= 1                       # the adjoint is traceless: one singlet removed
    return cnt


P("=" * 100)
P("1 -- THEIR (3.20), REBUILT FROM SCRATCH")
P("=" * 100)
P("  their (3.20):  N_Ad^(+-) = 2(pq + rs),   N_Ad^(-+) = 2(pr + qs)")
P("")
P("  %-16s %10s %10s %12s %12s %8s" % ("[p;q,r;s]", "N^(+-)", "2(pq+rs)", "N^(-+)", "2(pr+qs)",
                                       "agree"))
ok = True
for sz in CONFIGS:
    p, q, r, s = sz
    c = count_adjoint(sz)
    a1, a2 = c.get((1, -1), 0), c.get((-1, 1), 0)
    f1, f2 = 2 * (p * q + r * s), 2 * (p * r + q * s)
    ok &= (a1 == f1 and a2 == f2)
    P("  %-16s %10d %10d %12d %12d %8s" % (str(sz), a1, f1, a2, f2, a1 == f1 and a2 == f2))
P("")
P("  CONTROL -- their (3.20) reproduced without using their paper : %s" % ok)
assert ok, "the parity block count does not reproduce (3.20) -- stop, our count is wrong"

P("")
P("=" * 100)
P("2 -- THE SUM THAT (4.7) SUBSTITUTES")
P("=" * 100)
P("  N_Ad^(+-) + N_Ad^(-+) = 2(pq+rs) + 2(pr+qs) = 2[p(q+r) + s(q+r)] = 2 (p+s)(q+r)")
P("")
P("  %-16s %16s %16s %14s" % ("[p;q,r;s]", "N^(+-)+N^(-+)", "2(p+s)(q+r)", "(p+s)(q+r)"))
for sz in CONFIGS:
    p, q, r, s = sz
    c = count_adjoint(sz)
    tot = c.get((1, -1), 0) + c.get((-1, 1), 0)
    P("  %-16s %16d %16d %14d" % (str(sz), tot, 2 * (p + s) * (q + r), (p + s) * (q + r)))
    assert tot == 2 * (p + s) * (q + r)
P("")
P("  the sum is 2(p+s)(q+r).  Not (p+s)(q+r).")

P("")
P("=" * 100)
P("3 -- THEIR (4.5) + (4.6), SUMMED, AGAINST THE FIRST LINE OF (4.7)")
P("=" * 100)
P("  (4.5), gauginos:")
P("     -4[ N^(++)(e0+v(b)) + N^(--)(e0+v(-b)) + N^(+-)(e0+v(b+1/2)) + N^(-+)(e0+v(-b-1/2)) ]")
P("  (4.6), gauge + ghost + sigma:")
P("     (2N^(++)+2N^(--))(e0+de) + (2N^(--)+2N^(++))(e0-de)")
P("       + (2N^(+-)+2N^(-+))(e0+v(1/2)) + (2N^(-+)+2N^(+-))(e0+v(1/2))")
P("")
P("  v is even and w(b) = v(1/2) + v(b) - v(b+1/2), so with S = N^(+-)+N^(-+) and")
P("  N^(++)+N^(--) = N^2-1-S the v-dependent part of the sum is")
P("")
P("     -4(N^2-1-S) v(b)  -  4 S v(b+1/2)  +  4 S v(1/2)")
P("        = -4(N^2-1) v(b)  +  4 S [ v(b) - v(b+1/2) + v(1/2) ]")
P("        = -4(N^2-1) v(b)  +  4 S w(b) .")
P("")
P("  which is their (4.7) FIRST line, exactly:  -4(N^2-1)v + 4(N_Ad^(+-) + N_Ad^(-+)) w.")
P("  The Delta-epsilon terms cancel between the two lines of (4.6), as they say.")

P("")
P("=" * 100)
P("4 -- VERDICT, WITH THEIR OWN (4.13) AS ARBITER")
P("=" * 100)
P("  (4.7) line 1 : -4(N^2-1)v + 4 S w        with S = 2(p+s)(q+r)  ->  + 8(p+s)(q+r) w")
P("  (4.7) line 2 : -4(N^2-1)v + 4 (p+s)(q+r) w                     ->  half of line 1")
P("")
P("  (4.13) collects the w-bracket as  4 { (p+s)(q+r)(2 - 2n_Ad^(+) + ...) + ... } w, so with no")
P("  matter at all its gauge part is  4 * (p+s)(q+r) * 2 * w = 8 (p+s)(q+r) w.")
P("")
P("  Two of their equations give 8 and one gives 4, and the two are the ones everything else")
P("  rests on.  So the dropped factor is in the SECOND line of (4.7), it is a slip in the")
P("  substitution and not in the physics, and nothing downstream is affected: the +2 offset in")
P("  their a is the correct one.")
P("")
P("  WHAT IT MEANS FOR US.  Their gauge sector contributes the EVEN INTEGER +2 to a and NOTHING")
P("  to b.  Both functionals can therefore vanish, and they treat those cases explicitly --")
P("  their section (i) is 'The case with a = 0', and their (4.15) lists")
P("  'n_F^(+) = n_F^(-) => completely degenerate'.  In the model of this paper the gauge sector")
P("  contributes -7/2, a half-integer, so 8D, 2U and 2W are odd for every content and the")
P("  degenerate case they enumerate cannot be reached.  Same architecture, opposite outcome,")
P("  and the difference is one number.")
