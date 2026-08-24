#!/usr/bin/env sage
# lattice_lift.sage -- the lifted Smith form, with the transformation matrices sympy does not give.
#
#   Copyright (c) 2026 Carles Marin. All rights reserved.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
#
# lattice_lift.py computes the invariant factors of the lattice the eight multiplets generate in
# Z^5 = (A_4, 8D, 2U, V, 2W) and gets Z_2 x Z_72 x Z_2592.  Two things are missing there and both
# need the TRANSFORMATION matrices: an independent check of the factors themselves, and the
# congruences written out as conditions on the coordinates, which is what a reader can use.
#
# Sage returns D, U, V with D = U*A*V, all unimodular.  If A's rows generate L, then a vector
# y in Z^5 lies in L exactly when the coordinates of y*V in the Smith basis are divisible by the
# corresponding invariant factors -- so the columns of V give the congruences.
#
# The matrix is RETYPED from the paper's moment definitions rather than imported, for the same
# reason lattice_rank.sage retypes its own: sharing the generator table shares its bugs.
# The rows are 7(+,+), 7(+,-), 28(+,+), 28(+,-), 48(+,+), 48(+,-), 84(+,+), 84(+,-).
#
# Run:  docker run --rm -v "$PWD":/w -w /w sagemath/sagemath:latest sage lattice_lift.sage

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
COORD = ["A_4", "8D", "2U", "V", "2W"]

A = Matrix(ZZ, [
    [  0,  -6,   2,  0,   2],     # 7(+,+)
    [  1,   8,   0,  0,  -2],     # 7(+,-)
    [ 17,  16,  40,  0,   6],     # 28(+,+)
    [  4,   2,  34,  0,  -6],     # 28(+,-)
    [ 18,   0,  48,  0,  12],     # 48(+,+)
    [  8,  28,  36,  0, -12],     # 48(+,-)
    [ 68,  10, 346,  0,  18],     # 84(+,+)
    [109,  80, 168, 81, -18],     # 84(+,-)
])
GAUGE = vector(ZZ, [-18, -27, -39, 0, -3])

print("=" * 92)
print("1 -- THE LIFTED LATTICE")
print("=" * 92)
print("   %-10s %8s %8s %8s %8s %8s" % ("multiplet", *COORD))
for i in range(8):
    print("   %-10s %8d %8d %8d %8d %8d" % (NAMES[i], *A[i]))
print("   %-10s %8d %8d %8d %8d %8d" % ("gauge", *GAUGE))
print()
print("   rank over ZZ : %d   (must be 5, or the five are not coordinates)" % A.rank())
assert A.rank() == 5

print()
print("=" * 92)
print("2 -- SMITH NORMAL FORM, WITH TRANSFORMS")
print("=" * 92)
D, U, V = A.smith_form()
assert D == U * A * V and U.is_invertible() and V.is_invertible()
inv = [D[i, i] for i in range(5)]
print("   D = U A V, verified;  invariant factors : %s" % inv)
q = [d for d in inv if d not in (0, 1, -1)]
print("   Z^5 / L  ~=  %s        index %d" % (" x ".join("Z_%s" % d for d in q), prod(inv)))
print()
print("   the projections already in the paper, recomputed here for continuity:")
for cols, label in ([[0, 1], "(A_4, 8D)"], [[0, 1, 2, 3], "(A_4, 8D, 2U, V)"]):
    Dp = A.matrix_from_columns(cols).smith_form()[0]
    ip = [Dp[i, i] for i in range(len(cols))]
    qp = [d for d in ip if d not in (0, 1, -1)]
    print("      %-20s -> %s" % (label, " x ".join("Z_%s" % d for d in qp) or "trivial"))

print()
print("=" * 92)
print("3 -- THE CONGRUENCES, WRITTEN OUT")
print("=" * 92)
print("   y in Z^5 lies in L iff (y*V)_i is divisible by D_ii for every i.  The columns of V with")
print("   a non-unit invariant factor are therefore the congruences, one per factor:")
print()
for i in range(5):
    d = D[i, i]
    if d in (0, 1, -1):
        continue
    col = V.column(i)
    terms = " + ".join("%s*%s" % (col[j], COORD[j]) for j in range(5) if col[j] != 0)
    print("      %s  ==  0   (mod %s)" % (terms, d))
    # and it must hold on every generator and on the gauge sector
    ok = all((A[r] * col) % d == 0 for r in range(8))
    print("         holds on all eight multiplets : %s" % ok)
    assert ok
    print("         value on the gauge sector      : %s  (mod %s) = %s"
          % (GAUGE * col, d, (GAUGE * col) % d))

print()
print("=" * 92)
print("4 -- THE PARITY CHARACTER, AND IT RUNS ON THREE COORDINATES")
print("=" * 92)
print("   Independently of the Smith basis: every matter row is even in 8D, 2U AND 2W, and the")
print("   gauge sector is odd in all three.  So 8D == 2U == 2W == 1 (mod 2) for every content.")
print("   Theorem 1 and the 2W theorem are one congruence read twice; 2U is a third reading of")
print("   the same character, with no physical statement attached to it yet.")
print()
ODD = [1, 2, 4]
print("   %-10s %10s %10s %10s" % ("multiplet", "8D mod 2", "2U mod 2", "2W mod 2"))
for i in range(8):
    print("   %-10s %10d %10d %10d" % (NAMES[i], A[i][1] % 2, A[i][2] % 2, A[i][4] % 2))
print("   %-10s %10d %10d %10d" % ("gauge", GAUGE[1] % 2, GAUGE[2] % 2, GAUGE[4] % 2))
assert all(all(A[i][c] % 2 == 0 for c in ODD) for i in range(8))
assert all(GAUGE[c] % 2 == 1 for c in ODD)
print()
print("   and the source is ONE number: the gauge sector's only half-integer coefficient is the")
print("   antiperiodic charge-one term, -7/2, which enters B_2, B_4 and W alike.  A_4 and V are")
print("   NOT forced -- so the character selects three of the five, which is a real statement:")
for c in (0, 3):
    par = sorted({(GAUGE[c] + A[i][c]) % 2 for i in range(8)})
    print("      %-4s over gauge + one multiplet takes parities %s" % (COORD[c], par))
    assert len(par) > 1

print()
print("=" * 92)
print("VERDICT")
print("=" * 92)
print("   The lift is honest: rank 5, so (A_4, 8D, 2U, V, 2W) really are five independent")
print("   coordinates and not four with a passenger.  Z^5/L = %s, three congruences where the"
      % " x ".join("Z_%s" % d for d in q))
print("   projection to four had two, and the new Z_2 is the parity character.  Carles's suspicion,")
print("   checked and sharpened: it is one character, it runs on THREE of the five coordinates,")
print("   and its source is the single -7/2 of the antiperiodic charge-one gauge term.")
