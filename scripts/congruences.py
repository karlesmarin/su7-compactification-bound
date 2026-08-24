#!/usr/bin/env python3
"""congruences.py -- the two congruences the paper left unread, decoded.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHERE THIS COMES FROM.  lattice_lift.sage established that the five coordinates
(A_4, 8D, 2U, V, 2W) are a lattice of index 373248 in Z^5, with

    Z^5 / L  ~=  Z_2 x Z_72 x Z_2592 ,

and read exactly ONE of the three factors: the Z_2 is the parity character, which is Theorem 1
(8D odd) and the 2W theorem at once, both inherited from the single -7/2 of the antiperiodic
charge-one gauge term.  Section \\S\\ref{sec:open} of the paper then says of the other two, in
as many words, that what they FORBID a bulk content from doing "is not known here".

This file reads them.  The method is the obvious one and the reason it was not done is that the
congruences as Smith prints them are unreadable -- five-digit coefficients modulo 72 and 2592
tell you nothing.  Split them by the Chinese remainder theorem into prime-power pieces, put each
piece in lowest terms, and five readable statements come out, one of them automatic.

WHAT COMES OUT.  For every bulk content -- gauge sector plus any non-negative matter:

  (1) 8D + 2U + 2W == 1  (mod 2)                    the Z_2 factor
  (2) 8D == 2U + 4(2W)  (mod 8)              \\
  (3) A_4 + 8D + 3(2U) == V  (mod 9)          |     the Z_72 factor
  (4) 2W == -(2 A_4 + 12(8D) + 3(2U))  (mod 32)\\
  (5) V == 0  (mod 81)                         |    the Z_2592 factor; (5) is automatic,
                                                    since V is 81 x a multiplicity

The moduli multiply to 2 * 8 * 9 * 32 * 81 = 373248, the index, so the five are independent and
jointly complete: this is the WHOLE of what the lattice forbids, not a sample of it.

THREE THINGS ARE WORTH SAYING AND ONE IS NOT.

  The one that is not: "(4) determines 2W modulo 32 from the other coordinates" sounds like a
  discovery and is not one.  Every finite-index sublattice pins every coordinate modulo something
  -- that is what finite index means -- and (4) pins 2U just as well as 2W, since 17 and 27 are
  both invertible mod 32.  It is the content of the Smith form read out, not a structure on top
  of it.  Recorded because it was the first headline this file had.

  THE PARITY CHARACTER IS NOT ONE OF THE FIVE FACTORS.  It is the joint mod-2 shadow of (1), (2)
  and (4) -- (2) gives 8D == 2U, (4) gives 2U == 2W, and (1) then forces all three odd.  So no
  single factor of Z^5/L carries Theorem 1, and in particular the oddness of 2U, which
  \\S\\ref{sec:open} records as having "no physical statement attached to it yet", is EQUIVALENT
  to the oddness of 2W by (4) alone -- a congruence that never mentions the -7/2.

  (3) IS ALREADY IN THE PAPER, TWICE, UNRECOGNISED.  Reduced mod 3 it says A_4 + 8D == 0, which
  is the (t+1) % 3 filter the ceiling search steps by AND the observation that \\cite{KM25}'s five
  published rows "give 8D + A_4 = 273, 300, 198, 234, 270, all divisible by three".  Those are one
  fact, and (3) is the mod-NINE statement they are the shadow of: carrying 2U buys a power of three.

AND THEY CANNOT MOVE THE CEILING -- structurally, not by inspection.  The ceiling is decided by
walking the (A_4, 8D) plane, and the projection of the content lattice to that plane already has
index exactly 6.  So the mod-6 rule is the complete statement there and nothing involving 2U, V or
2W can cut the plane further: whatever the new congruences forbid, they forbid it inside a fibre.
The one way out would be a vertex whose fibre they empty; \\S5b tests the three vertices the paper
argues about and each merely pins 2U to one class mod 3, which the multiplets realise.

AND ONE MORE ROW OF THE HONEST-LIMITS TABLE FALLS OUT FOR FREE.  It reads "what the second
congruence of Z^4/L = Z_18 x Z_648 forbids: the group is computed, its conditions on the
multiplicities are not written out".  That is these five with 2W eliminated, and (4) is what
eliminates it.  The four projected laws are 8D == 0 (mod 2), 8D - 2U == 0 (mod 8),
A_4 + 8D + 3(2U) - V == 0 (mod 9) and V == 0 (mod 81) on the matter lattice, with residues
1, 4, 0, 0 on a content; index 2 x 8 x 9 x 81 = 11664 = 18 x 648, confirmed by closing the
subgroup of (Z/32)^4 and (Z/81)^4 by brute force rather than trusting the elimination.

Sections:
  1  the three congruences as Smith prints them, and why they are unreadable
  2  the CRT split, each piece in lowest terms
  3  the five statements, verified on every generator, the witness and the five published rows
  4  (4) solved for 2W, and the half of the parity character it absorbs
  5  (3) mod 3 is the paper's own divisibility, recognised
  5b why none of this can move the ceiling
  5c the same elimination closes the four-coordinate row of the honest-limits table
  6  controls

Run:  python congruences.py
"""
import itertools
import json
import math
import pathlib
import sys
from fractions import Fraction as Fr

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

KEYS = [(1, 1), (1, 2), (1, 3), (-1, 1), (-1, 2), (-1, 3)]
NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
COORD = ["A_4", "8D", "2U", "V", "2W"]


def uvec(tt):
    u = [Fr(0)] * 6
    for m, s, c in tt:
        u[KEYS.index((int(s), int(round(c))))] += Fr(m).limit_denominator(8)
    return u


def five(u):
    """(A_4, 8D, 2U, V, 2W) -- lattice_lift.py's map, retyped so this file stands alone."""
    A2, B2 = u[0] + 4 * u[1] + 9 * u[2], u[3] + 4 * u[4] + 9 * u[5]
    A4, B4 = u[0] + 16 * u[1] + 81 * u[2], u[3] + 16 * u[4] + 81 * u[5]
    return [A4, 8 * A2 - 6 * B2, 2 * (16 * u[1] + B4), 81 * u[2],
            2 * (-(u[0] + u[2]) + (u[3] + u[5]))]


def coords(content):
    """the five coordinates of a CONTENT -- table() already carries the gauge sector."""
    f = five(uvec(table(content)))
    assert all(Fr(x).denominator == 1 for x in f), "a coordinate came out non-integral"
    return [int(x) for x in f]


MATTER = [[int(x) for x in five(uvec(terms(*r)))] for r in REPS]
GAUGE5 = coords([])
dot = lambda a, b: sum(x * y for x, y in zip(a, b))

# Smith's own output, retyped from outputs/lattice_lift_sage.txt.  Modulus -> coefficients.
SMITH = {2:    [-42, 15, 5, 40, 39],
         72:   [32, -13, -3, -32, -36],
         2592: [-810, 324, 81, 800, 891]}
SPLIT = {2: [2], 72: [8, 9], 2592: [32, 81]}

# ================================================================= 1
P("=" * 100)
P("1 -- THE THREE CONGRUENCES AS SMITH PRINTS THEM, AND WHY THEY ARE UNREADABLE")
P("=" * 100)
P("  A vector lies in L (the lattice the eight MATTER multiplets generate) exactly when each of")
P("  these vanishes.  The gauge sector is NOT in L -- that is what the Z_2 measures -- so a bulk")
P("  content sits at the gauge sector's residue, printed in the last column.")
P("")
P("  %-6s %-52s %14s" % ("mod", "coefficients on (A_4, 8D, 2U, V, 2W)", "gauge residue"))
for q, c in SMITH.items():
    ok = all(dot(c, r) % q == 0 for r in MATTER)
    assert ok, "congruence mod %d does not vanish on the matter generators" % q
    P("  %-6d %-52s %14d" % (q, str(c), dot(c, GAUGE5) % q))
P("")
P("  Nothing is legible there, which is why they were left unread.  But 72 = 8 x 9 and")
P("  2592 = 32 x 81 are coprime factorisations, so each congruence is really TWO, and each of")
P("  those has small coefficients.")

# ================================================================= 2
P("")
P("=" * 100)
P("2 -- THE CRT SPLIT, EACH PIECE IN LOWEST TERMS")
P("=" * 100)
pieces = []
for q, c in SMITH.items():
    for pk in SPLIT[q]:
        red = [x % pk for x in c]
        g = pk
        for x in red:
            g = math.gcd(g, x)
        m = pk // g
        red = [(x // g) % m for x in red]
        pieces.append(dict(parent=q, mod=m, coef=red,
                           rhs=dot(SMITH[q], GAUGE5) // 1 % pk // g % m if g else 0))
# the residue of a CONTENT is the gauge one, reduced the same way
for pc in pieces:
    q, pk = pc["parent"], pc["mod"]
    pc["rhs"] = dot(pc["coef"], GAUGE5) % pk

P("  %-8s %-8s %-40s %10s" % ("from", "modulus", "coefficients", "content =="))
for pc in pieces:
    P("  mod %-4d mod %-4d %-40s %10d"
      % (pc["parent"], pc["mod"], str(dict(zip(COORD, pc["coef"]))), pc["rhs"]))
P("")
prod = 1
for pc in pieces:
    prod *= pc["mod"]
P("  the five moduli multiply to %d, and the index of L in Z^5 is 373248 : %s"
  % (prod, prod == 373248))
assert prod == 373248, "the pieces do not account for the whole quotient -- the split is lossy"
P("  So these five ARE the whole of what the lattice forbids: independent, and jointly complete.")

# ================================================================= 3
P("")
P("=" * 100)
P("3 -- THE FOUR STATEMENTS, WRITTEN AS SOMEONE WOULD SAY THEM")
P("=" * 100)


WITNESS = [("7", 1, 1, 16), ("28", 1, -1, 1), ("48", 1, 1, 1), ("48", 1, -1, 4), ("84", 1, 1, 1)]
PROBES = [("the paper's witness", coords(WITNESS))]
for lbl, content, a_them, mh_them, invR in T1:
    PROBES.append(("\\cite{KM25} row %s" % lbl, coords(content)))

# Every statement is LINEAR, so it has one form with two readings: zero on the matter lattice,
# and the gauge sector's own residue on a content.  Writing the target as f(gauge) rather than
# by hand is what keeps the two readings from drifting apart -- the first version of this file
# hard-coded 0 for both and (a) then "failed" on all eight generators, which was the encoding
# and not the arithmetic.
STATEMENTS = [
    ("(1)  8D + 2U + 2W == 1 (mod 2)",
     lambda A, D, U, V, W: D + U + W, 2,
     "the Z_2 factor -- NOT the parity character by itself; see below"),
    ("(2)  8D == 2U + 4(2W) (mod 8)",
     lambda A, D, U, V, W: D - U - 4 * W, 8,
     "the Z_72 factor, two-adic half"),
    ("(3)  A_4 + 8D + 3(2U) == V (mod 9)",
     lambda A, D, U, V, W: A + D + 3 * U - V, 9,
     "the Z_72 factor, three-adic half"),
    ("(4)  2W == -(2A_4 + 12(8D) + 3(2U)) (mod 32)",
     lambda A, D, U, V, W: W + 2 * A + 12 * D + 3 * U, 32,
     "the Z_2592 factor, two-adic half -- and it SOLVES for 2W"),
    ("(5)  V == 0 (mod 81)",
     lambda A, D, U, V, W: V, 81,
     "the Z_2592 factor, three-adic half -- automatic, V is 81 x a multiplicity"),
]
P("  %-46s %10s %10s" % ("statement", "on L", "on content"))
for name, f, mod, note in STATEMENTS:
    rhs = f(*GAUGE5) % mod
    badL = [NAMES[j] for j in range(8) if f(*MATTER[j]) % mod]
    badC = [lbl for lbl, v in [("gauge", GAUGE5)] + PROBES if (f(*v) - rhs) % mod]
    P("  %-46s %10s %10s" % (name, "0" if not badL else "FAILS", str(rhs) if not badC else "FAILS"))
    P("  %-46s   %s" % ("", note))
    assert not badL and not badC, "%s fails on %s %s" % (name, badL, badC)
P("")
P("  The %d objects tested are the eight matter generators, the gauge sector, the witness of"
  % (9 + len(PROBES)))
P("  Figure~\\ref{fig:vacuum} and the five published rows of \\cite{KM25} -- contents built for")
P("  other purposes entirely.")
P("")
P("  AND THE PARITY CHARACTER IS NOT ONE OF THESE FIVE.  It is the joint mod-2 shadow of three")
P("  of them, which is a sharper thing to be able to say.  Reduce (2) and (4) modulo 2:")
P("       (2) -> 8D == 2U          (the 4(2W) dies)")
P("       (4) -> 2W == 3(2U) == 2U (2A_4 and 12(8D) die)")
P("  so 8D == 2U == 2W already, and (1) then reads 3 x 8D == 1, i.e. 8D odd.  Theorem 1, the 2W")
P("  theorem and the unexplained oddness of 2U are one consequence of three congruences, and no")
P("  single lattice factor carries it.")
par = all(v[1] % 2 == v[2] % 2 == v[4] % 2 == 0 for v in MATTER) and \
      all(v[1] % 2 == v[2] % 2 == v[4] % 2 == 1 for v in [GAUGE5] + [p[1] for p in PROBES])
P("  CONTROL -- 8D, 2U, 2W all even on every generator and all odd on every content : %s" % par)
assert par, "the derived parity character does not hold, so the derivation is wrong"

# ================================================================= 4
P("")
P("=" * 100)
P("4 -- (4) SOLVED FOR 2W, AND THE HALF OF THE PARITY CHARACTER IT ABSORBS")
P("=" * 100)
inv27 = pow(27, -1, 32)
P("  27 is invertible mod 32, so (4) solves for 2W:")
P("")
P("      2W  ==  -( 2 A_4 + 12 (8D) + 3 (2U) )   (mod 32)")
P("")
P("  V is absent.  So knowing (A_4, 8D, 2U) pins 2W modulo 32, i.e. W modulo 16.  Verified:")
P("")
P("  %-26s %8s %10s %10s" % ("object", "2W", "predicted", "agree"))
allok = True
for lbl, v in [(NAMES[j], MATTER[j]) for j in range(8)] + [("gauge", GAUGE5)] + PROBES:
    A, D, U, V, W = v
    pred = (-(2 * A + 12 * D + 3 * U)) % 32
    ok = (W - pred) % 32 == 0
    allok &= ok
    P("  %-26s %8d %10d %10s" % (lbl, W, pred, ok))
assert allok, "(4) does not solve for 2W after all"
P("")
P("  AND IT ABSORBS ONE HALF OF THE PARITY CHARACTER.  Reduce (4) mod 2: 2A_4 and 12(8D) are")
P("  even, so 2W == 3(2U) == 2U (mod 2).  The theorem '2W is odd' and the unexplained '2U is")
P("  odd' of \\S\\ref{sec:open} are therefore the SAME parity, and (4) is why -- one of them")
P("  follows from the other, on any content, with no reference to the -7/2.")
o2 = all((v[4] - v[2]) % 2 == 0 for v in [MATTER[j] for j in range(8)] + [GAUGE5]
         + [p[1] for p in PROBES])
P("  CONTROL -- 2W == 2U (mod 2) on every object above : %s" % o2)
assert o2, "the mod-2 shadow of (4) is not what it should be"

# ================================================================= 5
P("")
P("=" * 100)
P("5 -- (3) MOD 3 IS A DIVISIBILITY THE PAPER ALREADY PRINTS WITHOUT EXPLAINING")
P("=" * 100)
P("  Reduce (3) modulo 3.  V = 81 x (multiplicity of 84(+,-)) is divisible by 3, and so is")
P("  3(2U), so it collapses to")
P("")
P("      A_4 + 8D  ==  0   (mod 3)")
P("")
P("  which the paper states twice, both times as an observation: the ceiling search steps A_4 by")
P("  three, and \\cite{KM25}'s five rows 'give 8D + A_4 = 273, 300, 198, 234, 270, all divisible")
P("  by three'.  Those are the same fact, and (3) is the mod-NINE statement it is a shadow of:")
P("")
P("  %-26s %8s %8s %10s %8s %14s" % ("row", "A_4", "8D", "A_4 + 8D", "mod 3", "(3) mod 9"))
for lbl, v in PROBES:
    A, D, U, V, W = v
    P("  %-26s %8d %8d %10d %8d %14d" % (lbl, A, D, A + D, (A + D) % 3, (A + D + 3 * U - V) % 9))
P("")
P("  The 'mod 3' column is zero because it must be; the 'mod 9' column is zero because of a")
P("  statement nobody had written down.  Carrying 2U buys one more power of three.")

# ================================================================= 5b
P("")
P("=" * 100)
P("5b -- AND THEY CANNOT MOVE THE CEILING, FOR A STRUCTURAL REASON AND NOT BY INSPECTION")
P("=" * 100)
P("  The ceiling is decided on the (A_4, 8D) plane: \\S\\ref{sec:ceiling} walks A_4 down at fixed")
P("  8D, and every vertex it can visit is a point of the PROJECTION of the content set to those")
P("  two coordinates.  So the question 'do the new congruences forbid a vertex' is the question")
P("  'is the projection smaller than the mod-6 rule already says'.  It is not, and that is")
P("  decidable rather than a matter of opinion: project the lattice and take its index.")
P("")
proj = [[r[0], r[1]] for r in MATTER]
# index of the lattice spanned by the eight (A_4, 8D) rows: the gcd of the 2x2 minors, which is
# the second determinantal divisor and needs no Smith routine.
minors = [proj[i][0] * proj[j][1] - proj[j][0] * proj[i][1]
          for i in range(8) for j in range(i + 1, 8)]
idx2 = 0
for m in minors:
    idx2 = math.gcd(idx2, abs(m))
P("  index of the projected lattice in Z^2 : %d   (the paper's Z_6)" % idx2)
assert idx2 == 6, "the projection is not index 6 -- then the paper's mod-6 rule is not complete"
P("")
P("  Six, so the mod-6 rule 8D == 2A_4 + 3 IS the complete statement in those two coordinates,")
P("  and no congruence involving 2U, V or 2W can cut the plane further: whatever they forbid,")
P("  they forbid it in the fibre over a vertex and not among the vertices.  A search that walks")
P("  the plane therefore cannot be sped up by them, and a bound proved by exhausting each fibre")
P("  cannot be tightened by them either.")
P("")
P("  THE ONE PLACE THEY COULD STILL HAVE BITTEN is a vertex whose fibre they empty completely.")
P("  Test it where it would matter -- the deciding vertex and the two the paper argues about:")
P("")
P("  %-14s %10s %10s" % ("(A_4, 8D)", "A_4 mod 9", "2U mod 3 forced"))
for t, k in ((104, 1), (212, 1), (215, 1)):
    # (3) at 8D = k, V == 0 (mod 9):  A_4 + k + 3(2U) == 0 (mod 9)
    need = [u for u in range(3) if (t + k + 3 * u) % 9 == 0]
    P("  %-14s %10d %10s" % ("(%d, %d)" % (t, k), t % 9, need))
P("")
P("  Each vertex forces 2U into ONE class mod 3 rather than emptying the fibre, and 2U ranges")
P("  over all three classes among the eight multiplets, so no vertex is excluded.  The negative")
P("  is therefore real and not an artefact of looking in the wrong place.")
u3 = sorted({r[2] % 3 for r in MATTER})
P("  CONTROL -- 2U mod 3 over the eight multiplets : %s   (must be all three, or the" % u3)
P("             argument above is vacuous)")
assert u3 == [0, 1, 2], "2U does not range over all residues mod 3 -- re-examine the exclusion"

# ================================================================= 5c
P("")
P("=" * 100)
P("5c -- AND ELIMINATING 2W CLOSES A DIFFERENT OPEN ROW OF THE PAPER")
P("=" * 100)
P("  The honest-limits table carries a row reading 'what the second congruence of")
P("  Z^4/L = Z_18 x Z_648 forbids: the group is computed, its conditions on the multiplicities are")
P("  not written out'.  That is the FOUR-coordinate projection, (A_4, 8D, 2U, V), and it is the")
P("  present five with 2W eliminated -- so it costs nothing extra once these are in hand.")
P("")
P("  The three-adic pieces (3) and (5) never mention 2W, so they survive projection untouched.")
P("  The two-adic ones do mention it, and (4) is what removes it: it fixes 2W modulo 32 in terms")
P("  of the rest, so substituting it into (1) and (2) leaves conditions on (A_4, 8D, 2U) alone.")
P("")
P("      (1) mod 2 : 2W == 2U by (4), so 8D + 2U + 2U == 8D          ->   8D odd on a content")
P("      (2) mod 8 : the 4(2W) term is 4 x (2W mod 2)                ->   8D - 2U == 4 (2W mod 2)")
P("")
P("  and the second is the place to be careful, because the two readings differ.  On the matter")
P("  lattice 2W is EVEN, so 4(2W) dies and the law is 8D == 2U (mod 8); on a content 2W is odd,")
P("  so it does not, and the law is 8D == 2U + 4.  The first version of this section wrote the")
P("  content form and tested it against the generators, where it fails on all eight.  Written")
P("  homogeneously, as everything else here is, there is one law with two readings:")
P("")
P("      8D == 0 (mod 2)   |  1        8D - 2U == 0 (mod 8)   |  4      <- on L | on a content")
P("      A_4 + 8D + 3(2U) - V == 0 (mod 9)   |  0        V == 0 (mod 81)   |  0")
P("")
P("  index 2 x 8 x 9 x 81 = %d, and |Z^4/L_4| = 18 x 648 = %d" % (2 * 8 * 9 * 81, 18 * 648))
assert 2 * 8 * 9 * 81 == 18 * 648, "the projected pieces do not account for Z_18 x Z_648"
P("")
P("  VERIFIED BY BRUTE FORCE rather than by the elimination above, because the elimination is")
P("  exactly the sort of hand argument that quietly drops a case.  Build the subgroup of")
P("  (Z/32)^4 x (Z/81)^4 that the eight multiplets generate, and count it.")
P("")


def subgroup_index(mod, cols):
    """|(Z/mod)^len(cols) / <matter rows>| by closing the subgroup under addition."""
    gens = [tuple(r[c] % mod for c in cols) for r in MATTER]
    H = {tuple(0 for _ in cols)}
    frontier = [tuple(0 for _ in cols)]
    while frontier:
        nxt = []
        for h in frontier:
            for g in gens:
                v = tuple((a + b) % mod for a, b in zip(h, g))
                if v not in H:
                    H.add(v)
                    nxt.append(v)
        frontier = nxt
    return mod ** len(cols) // len(H), len(H)


FOUR = [0, 1, 2, 3]
i2, n2 = subgroup_index(32, FOUR)
i3, n3 = subgroup_index(81, FOUR)
P("  two-adic   : subgroup has %8d of 32^4 = %8d elements  ->  index %d" % (n2, 32 ** 4, i2))
P("  three-adic : subgroup has %8d of 81^4                 ->  index %d" % (n3, i3))
P("  product of the two indices : %d   against 18 x 648 = %d   agree : %s"
  % (i2 * i3, 18 * 648, i2 * i3 == 18 * 648))
assert i2 * i3 == 18 * 648, "the brute-force index of the projection is not Z_18 x Z_648"
P("")
P("  And the four written laws must cut exactly that much and no more -- same test as before,")
P("  on the projected coordinates:")
PROJ = [("8D", lambda A, D, U, V: D, 2),
        ("8D - 2U", lambda A, D, U, V: D - U, 8),
        ("A_4 + 8D + 3(2U) - V", lambda A, D, U, V: A + D + 3 * U - V, 9),
        ("V", lambda A, D, U, V: V, 81)]
badp = []
P("  %-24s %8s %10s %12s" % ("form", "mod", "on L", "on content"))
for name, f, mod in PROJ:
    rhs = f(*GAUGE5[:4]) % mod
    bl = [NAMES[j] for j in range(8) if f(*MATTER[j][:4]) % mod]
    bc = [lbl for lbl, v in [("gauge", GAUGE5)] + PROBES if (f(*v[:4]) - rhs) % mod]
    badp += [(x, name) for x in bl + bc]
    P("  %-24s %8d %10s %12s" % (name, mod, "0" if not bl else "FAILS",
                                 str(rhs) if not bc else "FAILS"))
assert not badp, "the projected laws fail on %s" % badp
P("  So that row of the table is answered too, and the answer is four lines long.")

# ================================================================= 6
P("")
P("=" * 100)
P("6 -- CONTROLS")
P("=" * 100)
P("  (i) EACH STATEMENT MUST BE ABLE TO REJECT.  A congruence that every integer vector")
P("      satisfies says nothing.  Take the paper's witness and perturb one coordinate by one:")
P("      the perturbed vector is not a content and every statement that sees that coordinate")
P("      has to notice.")
P("")
w = coords(WITNESS)
P("      %-42s %s" % ("statement", "  ".join("%-6s" % c for c in COORD)))
for name, f, mod, note in STATEMENTS:
    row = []
    for i in range(5):
        v = list(w)
        v[i] += 1
        row.append("catch" if (f(*v) - f(*GAUGE5)) % mod else "blind")
    P("      %-42s %s" % (name, "  ".join("%-6s" % x for x in row)))
P("")
P("      A 'blind' is not a failure: it says the statement does not involve that coordinate,")
P("      which is readable straight off the coefficients.  What would be a failure is a row of")
P("      five blinds, or a coordinate no statement catches.")
seen = [any((f(*(w[:i] + [w[i] + 1] + w[i + 1:])) - f(*GAUGE5)) % mod
            for name, f, mod, note in STATEMENTS) for i in range(5)]
P("      every coordinate caught by at least one statement : %s   %s" % (all(seen), seen))
assert all(seen), "some coordinate is unconstrained -- the lattice would not have finite index"

P("")
P("  (ii) AND THE STATEMENTS MUST BE THE SMITH CONGRUENCES AND NOT A REWRITE THAT LOST")
P("       INFORMATION.  Sweep a box of integer vectors and check that satisfying the five")
P("       statements is EXACTLY the same as satisfying Smith's three -- same set, not merely")
P("       implied one way.")
box, agree, inL = 0, 0, 0
for A in range(-6, 7):
    for D in range(-6, 7):
        for U in range(-6, 7):
            for V in (0, 81, 162, -81):
                for W in range(-6, 7):
                    v = [A, D, U, V, W]
                    box += 1
                    a = all((dot(c, v) - dot(c, GAUGE5)) % q == 0 for q, c in SMITH.items())
                    b = all((f(*v) - f(*GAUGE5)) % mod == 0 for name, f, mod, note in STATEMENTS)
                    agree += (a == b)
                    inL += a
P("       vectors swept : %d      the two descriptions agree on : %d      of which contents : %d"
  % (box, agree, inL))
assert agree == box, "the readable statements are NOT equivalent to Smith's -- one lost something"
assert inL > 0, "the sweep found no content at all, so the agreement is vacuous"
P("       Equivalent on every vector, and the sweep does contain contents, so it is not vacuous.")

P("")
P("  (iii) THE MODULI MUST ACCOUNT FOR THE INDEX EXACTLY.  Already asserted in section 2:")
P("        2 x 9 x 8 x 32 x 81 = %d = |Z^5 / L|.  If a piece had been dropped the product would" % prod)
P("        be a proper divisor and the description would be incomplete without saying so.")

(HERE / "outputs").mkdir(exist_ok=True)
(HERE / "outputs" / "congruences.json").write_text(json.dumps(dict(
    smith={str(q): c for q, c in SMITH.items()},
    pieces=[dict(parent=pc["parent"], mod=pc["mod"], coef=pc["coef"], rhs=pc["rhs"])
            for pc in pieces],
    index=prod,
    solved_2W="2W == -(2*A_4 + 12*8D + 3*2U) (mod 32)",
    probes={lbl: v for lbl, v in PROBES}), indent=1), encoding="utf-8")
P("")
P("archived: outputs/congruences.json")
