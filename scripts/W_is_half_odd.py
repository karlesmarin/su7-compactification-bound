#!/usr/bin/env python3
"""W_is_half_odd.py -- the stability criterion can never be marginal either, and for the same reason.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

HOW THIS SURFACED, and it was a formatting bug.
analytic_branch_nogo.py printed a window "19.50 < n < 19.50" while every quantity in it was
supposed to be an integer.  The half came from W: vacuum_constraint.py had printed the gauge
sector's W with %d, which silently truncated it, so a half-integer had been reading as -1 for a
day.  Chasing the .50 gives a theorem.

THE STATEMENT.
    W = sum_{c odd} m (-s),   and  F(1) - F(0) = (31/16) zeta(5) W.
Matter contributes integers to W.  The gauge sector does not: its term list is
(-1, +1, 2), (-2, +1, 1), (-7/2, -1, 1) -- the same half-integer that Theorem 1 hangs on -- and the
odd-charge part of it is 2 - 7/2 = -3/2.  So

    2W is an ODD INTEGER for every bulk content,

hence W is never zero, hence F(1) != F(0) always: THE TWO SYMMETRIC POINTS ARE NEVER DEGENERATE.
The orbifold of this model is never marginally stable, exactly as its electroweak breaking is never
marginal, and by the same half-integer.

Run:  python W_is_half_odd.py
"""
import itertools
import pathlib
import sys
from fractions import Fraction as Fr

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
cont = lambda v: [(REPS[j][0], REPS[j][1], REPS[j][2], v[j]) for j in range(8) if v[j]]


def Wex(v):
    """W in EXACT rationals -- the whole point is that it is not an integer."""
    return sum(Fr(-s) * Fr(m).limit_denominator(8)
               for m, s, c in table(cont(v)) if int(round(c)) % 2 == 1)


P("=" * 100)
P("1 -- WHERE THE HALF COMES FROM")
P("=" * 100)
P("  the gauge sector's own term list, and only its ODD charges matter:")
for m, s, c in GAUGE:
    P("     m = %-6s  s = %+d  c = %d   %s"
      % (Fr(m).limit_denominator(8), s, int(c),
         ("contributes -s*m = %s" % Fr(-s * m).limit_denominator(8)) if int(c) % 2 else
         "even charge, does not enter W"))
W0 = Wex([0] * 8)
P("  W(gauge) = %s" % W0)
P("")
P("  and every multiplet contributes an integer:")
P("  %-12s %8s   %-12s %8s" % ("multiplet", "W", "multiplet", "W"))
WV = [Wex([1 if i == j else 0 for i in range(8)]) - W0 for j in range(8)]
for j in range(0, 8, 2):
    P("  %-12s %8s   %-12s %8s" % (NAMES[j], WV[j], NAMES[j + 1], WV[j + 1]))
allint = all(w.denominator == 1 for w in WV)
P("")
P("  CONTROL -- every multiplet's W is an integer : %s" % allint)
P("  CONTROL -- the gauge W is NOT : %s   (it is %s)" % (W0.denominator != 1, W0))

P("")
P("=" * 100)
P("2 -- THEREFORE 2W IS ODD, ALWAYS")
P("=" * 100)
P("  2W = 2 W(gauge) + 2 sum n_j W_j = %s + even, and %s is odd."
  % (2 * W0, 2 * W0))
P("")
P("  Checked by brute force over every content of at most six multiplets:")
bad = notodd = 0
tot = 0
for n in range(0, 7):
    for combo in itertools.combinations_with_replacement(range(8), n):
        v = [0] * 8
        for j in combo:
            v[j] += 1
        w = Wex(v)
        tot += 1
        if w == 0:
            bad += 1
        if (2 * w).denominator != 1 or int(2 * w) % 2 == 0:
            notodd += 1
P("     contents tested            : %d" % tot)
P("     with W exactly zero        : %d   <-- must be 0" % bad)
P("     with 2W not an odd integer : %d   <-- must be 0" % notodd)

P("")
P("=" * 100)
P("3 -- WHAT IT MEANS, AND IT IS THE SAME SHAPE AS THEOREM 1")
P("=" * 100)
P("  F(1) - F(0) = (31/16) zeta(5) W, and W is never zero, so the potential's values at the two")
P("  symmetric points are NEVER equal.  The orbifold is never marginally stable: for any bulk")
P("  content whatsoever, one of the two vacua is strictly deeper, and which one is decided by a")
P("  sign that cannot vanish.")
P("")
P("  Theorem 1 says 8D is an odd integer, so D != 0 and electroweak breaking is never marginal.")
P("  Both statements come from the SAME half-integer, the -7/2 of their eq. (68): at the")
P("  curvature it survives as 6 x 7/2 = 21, odd; at the value it survives as 3/2, half-odd.")
P("  One gauge coefficient, two non-degeneracies.")

P("")
P("=" * 100)
P("4 -- AND THE NUMBERS THE PAPER PRINTS HAD BEEN TRUNCATED")
P("=" * 100)
P("  vacuum_constraint.py formatted W with %d.  The correct values:")
P("")
P("  %-26s %12s %12s" % ("content", "W (printed)", "W (exact)"))
rows = []
for label, content, a_them, mh_them, invR in T1:
    v = [0] * 8
    for rep, e, ep, mult in content:
        v[REPS.index((rep, e, ep))] += mult
    rows.append((label, v))
WIT = [0] * 8
WIT[0], WIT[1], WIT[3] = 17, 2, 57
rows.append(("witness (212,1)", WIT))
W104 = [0] * 8
W104[0], W104[3], W104[4], W104[5], W104[6] = 16, 1, 1, 4, 1
rows.append(("witness (104,1)", W104))
for label, v in rows:
    w = Wex(v)
    P("  %-26s %12d %12s" % (label, int(w), w))
P("")
P("  the sign never changes, so no verdict moves -- but the printed integers were wrong and")
P("  the half is the whole content of section 3.")
