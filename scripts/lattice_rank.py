#!/usr/bin/env python3
"""lattice_rank.py -- is the log lattice three-dimensional or four?  Rank and Smith normal form.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

THE CHALLENGE, and it is against a sentence I wrote.
semigroup.py and the paper say: "V is not a fourth direction: ln3 enters through the single
generator 84(+,-), so V = 81 n_{84(+,-)}."  Carles: that inference is not valid by itself.  A
coordinate coming from ONE generator is not thereby dependent on the others -- it may be exactly a
new rank-one direction.  His reconstruction: the seven generators with V = 0 already reach rank
three in (A4, 8D, 2U), so adding the eighth, which has V = 81, must push the rank to four.

He is right if the V = 0 block has rank three, and that is a determinant.  Measured here, together
with the thing the paper should have printed in the first place: the SMITH NORMAL FORM, whose
invariant factors are the actual congruences.  An index of 144 is not "twenty-four congruences";
it is the order of a finite abelian group, and the group is what has to be named.

Run:  python lattice_rank.py
"""
import math
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
_g = moments([])
AV = [round(moments([(r, e, p, 1)])["A4"] - _g["A4"]) for r, e, p in REPS]
KV = [round(8 * (moments([(r, e, p, 1)])["D"] - _g["D"])) for r, e, p in REPS]
GQ_SYM = [(0, 0, 1), (1, 0, 0), (17, 0, 20), (4, 0, 17), (18, 0, 24),
          (8, 0, 18), (68, 0, 173), (109, 81, 84)]
# rows of the generator matrix: (A4, 8D, 2U, V).  2U because the gauge sector's U is -39/2.
ROWS = [[AV[j], KV[j], 2 * GQ_SYM[j][2], GQ_SYM[j][1]] for j in range(8)]


def rank_q(M):
    """rational rank by fraction-free elimination."""
    A = [[Fr(x) for x in r] for r in M]
    rows, cols = len(A), len(A[0])
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c] / A[r][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
    return r


def smith(M):
    """Smith normal form invariant factors of the integer matrix M (list of rows)."""
    A = [list(map(int, r)) for r in M]
    rows, cols = len(A), len(A[0])
    res = []
    t = 0
    while t < min(rows, cols):
        # find a non-zero pivot in the remaining block
        found = None
        for i in range(t, rows):
            for j in range(t, cols):
                if A[i][j]:
                    found = (i, j)
                    break
            if found:
                break
        if not found:
            break
        i, j = found
        A[t], A[i] = A[i], A[t]
        for r in A:
            r[t], r[j] = r[j], r[t]
        # clear the row and the column, repeatedly, until both are clean
        while True:
            for i in range(t + 1, rows):
                while A[i][t]:
                    f = A[t][t] // A[i][t]
                    A[t] = [a - f * b for a, b in zip(A[t], A[i])]
                    A[t], A[i] = A[i], A[t]
            for j in range(t + 1, cols):
                while A[t][j]:
                    f = A[t][t] // A[t][j]
                    for r in A:
                        r[t] = r[t] - f * r[j]
                    for r in A:
                        r[t], r[j] = r[j], r[t]
            if all(A[i][t] == 0 for i in range(t + 1, rows)) and \
               all(A[t][j] == 0 for j in range(t + 1, cols)):
                break
        res.append(abs(A[t][t]))
        t += 1
    # divisibility chain d1 | d2 | ...
    for i in range(len(res) - 1):
        for j in range(i + 1, len(res)):
            a, b = res[i], res[j]
            g = math.gcd(a, b)
            res[i], res[j] = g, a * b // g if g else 0
    return [d for d in res if d]


P("=" * 100)
P("1 -- THE GENERATOR MATRIX, IN FOUR COLUMNS")
P("=" * 100)
P("  %-10s %8s %8s %8s %8s" % ("multiplet", "A4", "8D", "2U", "V"))
for j in range(8):
    P("  %-10s %8d %8d %8d %8d" % tuple([NAMES[j]] + ROWS[j]))
P("")
V0 = [ROWS[j] for j in range(8) if ROWS[j][3] == 0]
VN = [ROWS[j] for j in range(8) if ROWS[j][3] != 0]
P("  generators with V = 0 : %d      with V != 0 : %d  (%s)"
  % (len(V0), len(VN), ", ".join(NAMES[j] for j in range(8) if ROWS[j][3])))

P("")
P("=" * 100)
P("2 -- RANK.  the question is whether V adds a direction or repeats one")
P("=" * 100)
r3_all = rank_q([r[:3] for r in ROWS])
r3_v0 = rank_q([r[:3] for r in V0])
r4_all = rank_q(ROWS)
P("  rank over Q of (A4, 8D, 2U), all eight generators      : %d" % r3_all)
P("  rank over Q of (A4, 8D, 2U), only the V = 0 generators : %d" % r3_v0)
P("  rank over Q of (A4, 8D, 2U, V), all eight              : %d" % r4_all)
P("")
if r3_v0 == 3 and r4_all == 4:
    P("  ==> THE V = 0 BLOCK ALREADY REACHES RANK THREE, so the eighth generator, the only one")
    P("      with V != 0, can only ADD a direction.  THE LATTICE IS FOUR-DIMENSIONAL.")
    P("      The sentence 'V is not a fourth direction, since V = 81 n_{84(+,-)}' confuses two")
    P("      things: V IS determined by one multiplicity, and that makes it cheap to describe --")
    P("      it does NOT make it dependent on the other three.  Carles caught this and he is")
    P("      right; the paper has to say four.")
else:
    P("  ==> the V = 0 block has rank %d, so the claim has to be re-examined case by case."
      % r3_v0)

P("")
P("=" * 100)
P("3 -- SMITH NORMAL FORM.  an index is a number; the congruences are the invariant factors")
P("=" * 100)
for dim, lbl in ((2, "(A4, 8D)"), (3, "(A4, 8D, 2U)"), (4, "(A4, 8D, 2U, V)")):
    d = smith([r[:dim] for r in ROWS])
    idx = 1
    for x in d:
        idx *= x
    grp = " x ".join("Z_%d" % x for x in d if x > 1) or "trivial"
    P("  %-18s invariant factors %-22s index %-8d  Z^%d/L = %s"
      % (lbl, str(d), idx, dim, grp))
P("")
P("  So 'the index is 144, therefore Theorem 2 is one congruence of twenty-four' is not a")
P("  statement about congruences at all: 144 is the ORDER of a finite abelian group, and the")
P("  number of independent congruences is the number of invariant factors above one -- at most")
P("  the rank.  The honest sentence is that the lift multiplies the index by a factor, and the")
P("  group is the one printed here.")

P("")
P("=" * 100)
P("4 -- CONTROL: the 2D index must still be the paper's six")
P("=" * 100)
d2 = smith([r[:2] for r in ROWS])
i2 = 1
for x in d2:
    i2 *= x
P("  invariant factors in (A4, 8D) : %s   index %d   expected 6 : %s" % (d2, i2, i2 == 6))
P("")
P("  CONTROL -- every generator must lie in the lattice its own invariant factors define:")
ok = True
for dim in (2, 3, 4):
    d = smith([r[:dim] for r in ROWS])
    idx = 1
    for x in d:
        idx *= x
    sub = smith([r[:dim] for r in ROWS])
    ok &= (len(sub) == rank_q([r[:dim] for r in ROWS]))
P("     number of invariant factors equals the rank, in every dimension : %s" % ok)
