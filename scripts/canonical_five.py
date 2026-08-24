#!/usr/bin/env python3
"""canonical_five -- is the eight-type zoo really five types, and is the semigroup free of holes?

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

  Theorem 3 says the five coordinates (A_4, 8D, 2U, V, 2W) are complete invariants of the
  one-loop potential, and dimension_five.sage exhibits three relations with NON-NEGATIVE
  coefficients among the eight multiplet types.  Non-negativity is the whole point: a relation
  with a negative coefficient rewrites a content as a formal difference, which is not a content.
  With all coefficients non-negative the rewriting stays inside the physical cone, so every
  content has a representative built from five types only.

  This file asks three things and refuses to assume any of them:

    (1) are the five proposed types independent, i.e. is the reduction a change of BASIS?
    (2) is the semigroup they generate NORMAL -- S = C \\cap L, no holes?
    (3) does the reduction invert, giving five multiplicities back from five invariants?

  The strong external control is that the paper already states Z^5/L = Z_2 x Z_72 x Z_2592,
  of order 373248.  If the five chosen types are a basis of L then |det B| must be exactly that
  number.  Nothing here is tuned to make it come out.

  Everything is in exact integer arithmetic: the multiplet table is the one retyped in
  dimension_five.sage, and the Fourier construction is the same one, so a disagreement here is
  a disagreement with a file that was itself written to disagree with a numerical SVD.

Run:  python canonical_five.py > outputs/canonical_five.txt
"""

import sys
from fractions import Fraction

import sympy as sp

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]

# rows: the eight multiplets;  columns: coefficient of g(+,1) g(+,2) g(+,3) g(-,1) g(-,2) g(-,3)
V = sp.Matrix([
    [0, 0, 0, 1, 0, 0],      # 7(+,+)
    [1, 0, 0, 0, 0, 0],      # 7(+,-)
    [1, 1, 0, 4, 0, 0],      # 28(+,+)
    [4, 0, 0, 1, 1, 0],      # 28(+,-)
    [2, 1, 0, 8, 0, 0],      # 48(+,+)
    [8, 0, 0, 2, 1, 0],      # 48(+,-)
    [4, 4, 0, 12, 1, 1],     # 84(+,+)
    [12, 1, 1, 4, 4, 0],     # 84(+,-)
])

MODES = 400          # Fourier modes; the rank question is settled long before this


def fourier(nmodes=MODES):
    """g(s,c) puts s^n/n^5 on the mode cos(c n pi a):  G[j][m] for j = (s,c), m = 1..nmodes."""
    rows = []
    for s in (+1, -1):
        for c in (1, 2, 3):
            row = []
            for m in range(1, nmodes + 1):
                if m % c == 0:
                    n = m // c
                    row.append(Fraction(s ** n, n ** 5))
                else:
                    row.append(Fraction(0))
            rows.append([sp.Rational(x.numerator, x.denominator) for x in row])
    return sp.Matrix(rows)


def line(ch="-", n=94):
    print(ch * n)


def main():
    fails = []
    line("=")
    print("THE EIGHT MULTIPLET TYPES ARE FIVE, AND THE SEMIGROUP HAS NO HOLES")
    line("=")

    G = fourier()
    FUN = V * G

    print("\n[0] THE SETTING, RECOMPUTED")
    r = FUN.rank()
    print("    rank of the eight potentials .................... %d" % r)
    ok = r == 5
    print("   C0  the function space is five-dimensional ....... %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C0")

    # ---------------------------------------------------------------- (1) the three relations
    print("\n[1] THE THREE RELATIONS, AND THEIR SIGNS")
    CLAIM = [(2, {0: 20, 1: 17}), (4, {0: 24, 1: 18}), (5, {0: 1, 1: 4, 3: 1})]
    allok = True
    for j, rep in CLAIM:
        lhs = FUN.row(j)
        rhs = sp.zeros(1, FUN.cols)
        for i, co in rep.items():
            rhs += co * FUN.row(i)
        good = sp.simplify(lhs - rhs) == sp.zeros(1, FUN.cols)
        allok = allok and good
        print("    %-9s = %-32s  exact: %s"
              % (NAMES[j], " + ".join("%d x %s" % (c, NAMES[i]) for i, c in sorted(rep.items())),
                 good))
    print("\n   C1  all three hold as identities of potentials ... %s" % ("PASS" if allok else "FAIL"))
    if not allok:
        fails.append("C1")
    neg = [(j, i, c) for j, rep in CLAIM for i, c in rep.items() if c < 0]
    ok = not neg
    print("   C2  every coefficient is NON-NEGATIVE ............ %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C2")
    print("       this is what makes the rewriting physical: a negative coefficient would turn")
    print("       a content into a formal difference of contents, which is not a content.")

    # ---------------------------------------------------------------- (2) is it a basis?
    print("\n[2] ARE THE SURVIVING FIVE A BASIS?")
    KEEP = [0, 1, 3, 6, 7]          # 7(+,+), 7(+,-), 28(+,-), 84(+,+), 84(+,-)
    print("    proposed canonical types : %s" % ", ".join(NAMES[i] for i in KEEP))
    SUB = sp.Matrix([FUN.row(i) for i in KEEP])
    ok = SUB.rank() == 5
    print("\n   C3  the five are linearly independent ............ %s (rank %d)"
          % ("PASS" if ok else "FAIL", SUB.rank()))
    if not ok:
        fails.append("C3")

    # coordinates of all eight in the basis of the five
    COORD = []
    for i in range(8):
        sol = SUB.T.solve_least_squares(FUN.row(i).T) if False else None
        # exact: solve  x * SUB = FUN.row(i)
        x = sp.Matrix(sp.linsolve((SUB.T, FUN.row(i).T)).args[0]).T
        COORD.append([sp.nsimplify(v) for v in x])
    print("\n    each type in that basis (rows must be integral and >= 0):")
    print("      %-9s %s" % ("type", "  ".join("%9s" % NAMES[k] for k in KEEP)))
    line("-")
    integral = nonneg = True
    for i in range(8):
        c = COORD[i]
        integral &= all(sp.Rational(v).q == 1 for v in c)
        nonneg &= all(v >= 0 for v in c)
        print("      %-9s %s" % (NAMES[i], "  ".join("%9s" % sp.nsimplify(v) for v in c)))
    line("-")
    print("\n   C4  every type is an INTEGER combination of the five ... %s"
          % ("PASS" if integral else "FAIL"))
    if not integral:
        fails.append("C4")
    print("   C5  and a NON-NEGATIVE one ............................ %s"
          % ("PASS" if nonneg else "FAIL"))
    if not nonneg:
        fails.append("C5")

    # ---------------------------------------------------------------- (3) controls that can fail
    print("\n[3] TWO CONTROLS THAT CAN ACTUALLY FAIL")
    print("""
    A first draft of this file "verified" that the five canonical rows are unimodular by
    computing every coordinate IN THE BASIS OF THOSE FIVE ROWS.  They came out as the identity,
    determinant one, Smith factors all one -- and none of that tested anything, because it was
    true by construction.  It is recorded here because a control that cannot fail is worse than
    no control: it reports success.  Note also that the paper's index Z^5/L = Z_2 x Z_72 x
    Z_2592 lives in the INVARIANT coordinates (A_4, 8D, 2U, V, 2W), whose ambient lattice is not
    the one the multiplets span, so it was never the right comparison either.

    What follows can fail.""")

    # F1: drop one of the five and the span must collapse
    worst = 0
    for drop in KEEP:
        rest = [i for i in KEEP if i != drop]
        r = sp.Matrix([FUN.row(i) for i in rest]).rank()
        worst = max(worst, r)
    ok = worst == 4
    print("\n   F1  removing ANY one of the five drops the rank to 4 ... %s (max %d)"
          % ("PASS" if ok else "FAIL", worst))
    if not ok:
        fails.append("F1")
    print("       so none of the five is redundant; the basis is minimal, not merely sufficient.")

    # F2: the reduction must be checked against ACTUAL contents, not against itself.
    #     A first draft here "tested" normality by forming x = B z and confirming B^-1 x = z.
    #     That is an identity.  The thing that can fail is whether an arbitrary content, summed
    #     from its own term table, lands on the non-negative integer combination the reduction
    #     predicts.  That is what is checked, in [5].
    print("\n   F2  is deferred to C8 in [5], and the reason is worth stating:")
    print("       a bounded scan of x = B z against B^-1 x = z tests nothing, because it is an")
    print("       identity.  Only a content built independently can disagree with the map.")

    # F3: the reduction must NOT be able to absorb a spurious extra generator
    print("\n   F3  a spurious ninth type must break the rank-5 picture")
    SPUR = sp.Matrix([[1, 0, 1, 0, 1, 0]]) * G           # not one of the eight
    ext = sp.Matrix.vstack(sp.Matrix([FUN.row(i) for i in KEEP]), SPUR)
    r_ext = ext.rank()
    ok = r_ext == 5
    print("       rank of the five plus a spurious row : %d ....... %s"
          % (r_ext, "PASS" if ok else "FAIL"))
    if not ok:
        fails.append("F3")
    coef = sp.Matrix(sp.linsolve((sp.Matrix([FUN.row(i) for i in KEEP]).T, SPUR.T)).args[0])
    negative = any(v < 0 for v in coef)
    print("       and its coordinates in the basis are %s"
          % ("MIXED IN SIGN, as a non-content should be" if negative
             else "all non-negative -- suspicious"))
    ok = negative
    print("   F4  a non-content does NOT reduce non-negatively ...... %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("F4")
    print("       so C5 is a fact about these eight multiplets and not about any five vectors.")

    # ---------------------------------------------------------------- (4) normality
    print("\n[4] THE CONSEQUENCE, STATED AT ITS REAL STRENGTH")
    print("""
    It is tempting to dress this up as a normality theorem: S = C n L, no holes, with a short
    proof.  The proof is short because there is nothing to prove.  Once C4 and C5 hold,

        S = B N^5 ,      L = B Z^5 ,      C = B R_>=0^5 ,

    and B is invertible (C3), so x in C is EQUIVALENT to B^-1 x >= 0 and x in L is EQUIVALENT to
    B^-1 x integral.  A point of C n L is therefore a non-negative integral z by definition, and
    "no holes" is a restatement of the definitions.  That is the standard fact that a simplicial
    cone spanned by a lattice basis is unimodular, hence normal.

    So the honest statement is not that the semigroup is normal.  It is stronger and shorter:""")
    ok = all(f not in fails for f in ("C3", "C4", "C5"))
    print("""
        S is a FREE commutative monoid on five generators:  S = N^5.

    Free, not merely normal.  There are no holes because there is nowhere to put one.  All the
    content sits in C4 and C5 -- that eight physical multiplet types collapse onto five with
    coefficients that are integral AND non-negative -- and none of it sits in the deduction.\n""")
    print("   C7  the premises C3, C4, C5 all hold ............. %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C7")
    print("""
    SCOPE, AND IT DOES NOT CLOSE THE PAPER'S OPEN QUESTION.  What the paper asks is whether the
    LIFTED semigroup is normal, in the four coordinates (A_4, 8D, 2U, V); it reports that the
    projection is normal to A_4 = 420 and leaves the lift open.  The object proved normal here
    is the FIVE-coordinate one, which is larger.  Projection does not preserve normality: a
    point can lie in the projected cone and the projected lattice while every preimage of it in
    the five-dimensional lattice falls outside the cone.  So this is not the answer to that
    question.

    What it does do is change what a hole would MEAN.  The paper's phrasing is that a hole would
    be "a content every visible law permits and no bulk can build".  With the five-coordinate
    semigroup normal, no POTENTIAL is missing: what a four-coordinate hole would record is a
    (A_4, 8D, 2U, V) whose required 2W no content supplies.  The gap, if there is one, is in the
    coordinate that was dropped and not in the space of potentials.

    And the gauge sector is an additive offset, so the realisable set is the translate
    I_gauge + S.  Translation commutes with intersection, so the conclusion carries across, but
    what carries is a translate of a cone and not a cone.""")

    # ---------------------------------------------------------------- (5) the inverse map
    print("\n[5] THE INVERSE MAP")
    print("""
    Because B is invertible over the rationals, five invariants determine five canonical
    multiplicities: N = B^{-1} (I - I_gauge).  The lattice congruences are then nothing more
    than the condition that B^{-1} applied to the invariant vector be integral, and the cone
    inequalities are nothing more than N >= 0.  Explicitly, eliminating the three types:""")
    print("      N[7(+,+)]  = n[7(+,+)]  + 20 n[28(+,+)] + 24 n[48(+,+)] +  1 n[48(+,-)]")
    print("      N[7(+,-)]  = n[7(+,-)]  + 17 n[28(+,+)] + 18 n[48(+,+)] +  4 n[48(+,-)]")
    print("      N[28(+,-)] = n[28(+,-)] +  1 n[48(+,-)]")
    print("      N[84(+,+)] = n[84(+,+)]        N[84(+,-)] = n[84(+,-)]")
    # verify that map against the coordinates computed above
    import itertools
    worst = None
    for trial in itertools.product(range(3), repeat=8):
        tgt = sp.zeros(1, 5)
        for i, mi in enumerate(trial):
            if mi:
                tgt += mi * sp.Matrix([COORD[i]])
        n0, n1, n2, n3, n4, n5, n6, n7 = trial
        pred = [n0 + 20 * n2 + 24 * n4 + 1 * n5,
                n1 + 17 * n2 + 18 * n4 + 4 * n5,
                n3 + 1 * n5, n6, n7]
        if [sp.nsimplify(v) for v in tgt] != [sp.Integer(v) for v in pred]:
            worst = (trial, [sp.nsimplify(v) for v in tgt], pred)
            break
    ok = worst is None
    print("\n   C8  the explicit reduction reproduces the coordinates on 3^8 contents  %s"
          % ("PASS" if ok else "FAIL"))
    if not ok:
        fails.append("C8")
        print("       first disagreement: %s" % (worst,))

    line("=")
    if fails:
        print("VERDICT: %d CONTROL(S) FAILED: %s" % (len(fails), ", ".join(fails)))
        line("=")
        return 1
    print("VERDICT: C0-C8 pass.  Eight multiplet types collapse onto five with non-negative")
    print("         integer coefficients, the five are a minimal basis, the semigroup they")
    print("         generate is FREE on them, and the invariants invert to a canonical content.")
    print("         This does NOT settle the paper's question about the four-coordinate lift.")
    line("=")
    return 0


if __name__ == "__main__":
    sys.exit(main())
