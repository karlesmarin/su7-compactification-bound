#!/usr/bin/env sage
# dimension_five.sage -- the potential lives in FIVE dimensions, not six, by a route sharing no code.
#
#   Copyright (c) 2026 Carles Marin. All rights reserved.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
#
# semi_infinite.py says the six functions g(s,c)(a) = Re Li_5(s e^{i c pi a}) span a space of
# dimension five, and that the eight multiplets are only five directions with three non-negative
# relations between them.  It says so from a numerical SVD (2.8e-12 against 34.8) and from
# numpy's matrix_rank.  Both are floating point, and the claim contradicts a sentence printed
# twice in the paper.  A correction of that size does not get to rest on a tolerance.
#
# Here it is again in exact rational arithmetic, with nothing shared.  The trick is that these
# functions are their own Fourier series: g(s,c) puts the weight s^n / n^5 on the mode
# cos(c n pi a), so on the modes m = 1 .. M its coordinate vector is
#
#     v[m] = s^(m/c) / (m/c)^5   if c divides m,   0 otherwise
#
# and every question about the functions becomes a question about a 6 x M rational matrix.  No
# polylogarithms, no truncation error, no tolerance: two functions of this family are equal iff
# their coefficient vectors are.
#
# The multiplet table is RETYPED from the paper's term tables rather than imported, on purpose,
# for the same reason lattice_rank.sage retypes its matrix: sharing the generator table would
# share whatever bug produced it.
#
# Run:  docker run --rm -v "$PWD":/w -w /w sagemath/sagemath:latest sage dimension_five.sage

M = 240                                     # any M >= 6 settles it; 240 is free and leaves room
KEYS = [(1, 1), (1, 2), (1, 3), (-1, 1), (-1, 2), (-1, 3)]
KNAME = ["g(+,1)", "g(+,2)", "g(+,3)", "g(-,1)", "g(-,2)", "g(-,3)"]


def fourier(s, c):
    v = [QQ(0)] * M
    for n in range(1, M // c + 1):
        v[c * n - 1] = QQ(s) ** n / QQ(n) ** 5
    return v


G = Matrix(QQ, [fourier(s, c) for s, c in KEYS])

print("=" * 92)
print("1 -- THE SIX FUNCTIONS, AS EXACT FOURIER VECTORS ON %d MODES" % M)
print("=" * 92)
print("   rank over QQ : %d      (semi_infinite.py's SVD says 5)" % G.rank())
print("   kernel of the transpose, i.e. the relations among the six:")
K = G.kernel()
print("   dimension %d" % K.dimension())
for b in K.basis():
    print("      %s" % " + ".join("%s*%s" % (b[i], KNAME[i]) for i in range(6) if b[i]))
print()
print("   the duplication formula Li_s(z) + Li_s(-z) = 2^(1-s) Li_s(z^2) at s = 5 predicts")
print("   g(+,1) + g(-,1) - g(+,2)/16 = 0.  Exactly, as vectors:")
rel = G[0] + G[3] - G[1] / 16
print("      g(+,1) + g(-,1) - g(+,2)/16  is the zero vector : %s" % (rel == 0))
assert rel == 0, "the duplication relation fails"
print("   and the same shape one charge up needs c = 4 and c = 6, which this model has not:")
print("      g(+,2) + g(-,2) - g(+,4)/16 : g(+,4) is not a basis vector here, and")
print("      g(+,2) + g(-,2) is the zero vector : %s   <-- must be False" % (G[1] + G[4] == 0))
print("      g(+,3) + g(-,3) is the zero vector : %s   <-- must be False" % (G[2] + G[5] == 0))
assert G[1] + G[4] != 0 and G[2] + G[5] != 0, "a second relation appeared; the count is wrong"
assert G.rank() == 5, "the rank is not five"

print()
print("=" * 92)
print("2 -- THE EIGHT MULTIPLETS.  their (m, s, c) term tables, retyped from the paper")
print("=" * 92)
# rows: 7(+,+), 7(+,-), 28(+,+), 28(+,-), 48(+,+), 48(+,-), 84(+,+), 84(+,-)
# columns: the coefficient of g(+,1), g(+,2), g(+,3), g(-,1), g(-,2), g(-,3)
NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
V = Matrix(ZZ, [
    [0, 0, 0, 1, 0, 0],      # 7(+,+)   : (1,-1,1)
    [1, 0, 0, 0, 0, 0],      # 7(+,-)   : (1,+1,1)
    [1, 1, 0, 4, 0, 0],      # 28(+,+)  : (1,+1,1) (1,+1,2) (4,-1,1)
    [4, 0, 0, 1, 1, 0],      # 28(+,-)  : (4,+1,1) (1,-1,1) (1,-1,2)
    [2, 1, 0, 8, 0, 0],      # 48(+,+)  : (2,+1,1) (1,+1,2) (8,-1,1)
    [8, 0, 0, 2, 1, 0],      # 48(+,-)  : (8,+1,1) (2,-1,1) (1,-1,2)
    [4, 4, 0, 12, 1, 1],     # 84(+,+)  : (4,+1,1) (4,+1,2) (1,-1,1)+(11,-1,1) (1,-1,2) (1,-1,3)
    [12, 1, 1, 4, 4, 0],     # 84(+,-)  : (1,+1,1)+(11,+1,1) (1,+1,2) (1,+1,3) (4,-1,1) (4,-1,2)
])
FUN = V * G                                  # each row is the multiplet's own Fourier vector
print("   rank of the eight potentials over QQ : %d" % FUN.rank())
KER = FUN.kernel()
print("   kernel dimension                     : %d" % KER.dimension())
print()
print("   an integral basis of the kernel, reduced:")
B = KER.basis_matrix()
B = (B * B.denominator()).change_ring(ZZ)
B = B.saturation().hermite_form(include_zero_rows=False)
for b in B.rows():
    print("      %s" % " + ".join("%s*%s" % (b[i], NAMES[i]) for i in range(8) if b[i]))
print()
print("   and the three the run states, each checked as an identity of potentials:")
CLAIM = [("28(+,+)", 2, {0: 20, 1: 17}),
         ("48(+,+)", 4, {0: 24, 1: 18}),
         ("48(+,-)", 5, {0: 1, 1: 4, 3: 1})]
allok = True
for name, j, rep in CLAIM:
    lhs = FUN[j]
    rhs = sum(coef * FUN[i] for i, coef in rep.items())
    ok = (lhs == rhs)
    allok = allok and ok
    print("      %-9s = %-34s  exact : %s"
          % (name, " + ".join("%dx%s" % (c, NAMES[i]) for i, c in sorted(rep.items())), ok))
print()
print("   CONTROL -- all three are exact identities over QQ : %s" % allok)
assert allok and FUN.rank() == 5 and KER.dimension() == 3, "the five-direction claim fails"

print()
print("=" * 92)
print("3 -- AND THE FIVE GENERATORS REALLY DO GENERATE, WITH NON-NEGATIVE COEFFICIENTS")
print("=" * 92)
print("   the substitutions above have non-negative coefficients, so the reachable set of")
print("   potentials is the MONOID on 7(+,+), 7(+,-), 28(+,-), 84(+,+), 84(+,-) and not merely")
print("   the lattice.  That is what lets semi_infinite.py enumerate the fibre exhaustively.")
print()
BAS = [0, 1, 3, 6, 7]
S = Matrix(QQ, [FUN[i] for i in BAS])
print("   rank of the five generators : %d   (must be 5, or they do not span)" % S.rank())
assert S.rank() == 5
for name, j, rep in CLAIM:
    coords = S.solve_left(Matrix(QQ, [FUN[j]]))
    print("      %-9s in the five generators : %s   all >= 0 : %s"
          % (name, list(coords[0]), all(x >= 0 for x in coords[0])))
    assert all(x >= 0 for x in coords[0]), "a substitution is not non-negative"

print()
print("=" * 92)
print("VERDICT")
print("=" * 92)
print("   The paper says the potential lives in a six-dimensional space of functions.  It is")
print("   FIVE: the duplication formula at s = 5 removes one, and nothing removes a second.")
print("   The eight multiplets are five directions, with three non-negative relations, of which")
print("   Part VI knew one -- the two it did not need 37 and 42 multiplets, and its search was")
print("   capped at six.  All of it in exact rational arithmetic, on a retyped matrix.")
