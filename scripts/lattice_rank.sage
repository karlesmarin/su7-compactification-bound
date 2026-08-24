#!/usr/bin/env sage
# lattice_rank.sage -- the rank and Smith normal form again, by a route that shares no code.
#
#   Copyright (c) 2026 Carles Marin. All rights reserved.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
#
# lattice_rank.py computes the rank and the invariant factors with elimination routines I wrote
# by hand, and the whole correction of "three dimensions" to "four" rests on them.  A result that
# overturns a sentence in the paper should not rest on one implementation, least of all mine on
# the day I wrote it.  Sage has Smith normal form built in; if the two disagree, one of us is
# wrong and it matters which.
#
# The matrix is typed here rather than imported, on purpose: sharing the generator table would
# share whatever bug produced it.  It is the (A4, 8D, 2U, V) row per multiplet, in the order
# 7(+,+), 7(+,-), 28(+,+), 28(+,-), 48(+,+), 48(+,-), 84(+,+), 84(+,-).
#
# Run:  docker run --rm -v "$PWD":/w -w /w sagemath/sagemath:latest sage lattice_rank.sage

M = Matrix(ZZ, [
    [  0,  -6,   2,  0],
    [  1,   8,   0,  0],
    [ 17,  16,  40,  0],
    [  4,   2,  34,  0],
    [ 18,   0,  48,  0],
    [  8,  28,  36,  0],
    [ 68,  10, 346,  0],
    [109,  80, 168, 81],
])

print("=" * 92)
print("rank and Smith normal form, computed by Sage")
print("=" * 92)
print("  matrix is 8 x 4, rows = multiplets, columns = (A4, 8D, 2U, V)")
print("")

V0 = M.matrix_from_rows([i for i in range(8) if M[i, 3] == 0])
print("  rank of (A4, 8D, 2U), all eight rows      : %d" % M.matrix_from_columns([0, 1, 2]).rank())
print("  rank of (A4, 8D, 2U), only the V = 0 rows : %d" % V0.matrix_from_columns([0, 1, 2]).rank())
print("  rank of (A4, 8D, 2U, V), all eight        : %d" % M.rank())
print("")
print("  -> if the V = 0 block already has rank three, the eighth row can only add a direction.")
print("")

for dim, lbl in [(2, "(A4, 8D)"), (3, "(A4, 8D, 2U)"), (4, "(A4, 8D, 2U, V)")]:
    S = M.matrix_from_columns(range(dim)).smith_form()[0]
    d = [S[i, i] for i in range(min(S.nrows(), S.ncols())) if S[i, i] != 0]
    idx = prod(d)
    grp = " x ".join("Z_%s" % x for x in d if x > 1)
    print("  %-18s invariant factors %-24s index %-8s  Z^%d/L = %s"
          % (lbl, str(d), idx, dim, grp if grp else "trivial"))

print("")
print("=" * 92)
print("CONTROL -- these must match lattice_rank.py exactly:")
print("   rank 4 ; Z^2/L = Z_6 ; Z^3/L = Z_2 x Z_72 ; Z^4/L = Z_18 x Z_648")
print("=" * 92)
