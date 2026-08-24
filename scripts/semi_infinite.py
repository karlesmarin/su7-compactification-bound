#!/usr/bin/env python3
"""semi_infinite.py -- the true-vacuum condition as a semi-infinite integer program, and it closes.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHAT THIS IS FOR
----------------
The paper's open question 3 says, in its own words: W > 0 puts the electroweak point below the
OTHER SYMMETRIC POINT, and a third minimum elsewhere in [0,1] would be seen by neither, so the
9.22 TeV ceiling is a bound over a screened family; the clean way to close it is not another
screen but a wider programme.  This file is that programme.

The programme, and the reason it terminates.  The charges run over c in {1,2,3} and there are two
parities, so

    F(alpha) = sum_{s,c} u_{s,c} g_{s,c}(alpha) ,     g_{s,c}(a) = Re Li_5(s e^{i c pi a})

with SIX basis functions however much matter there is, and u linear in the multiplicities.

THE CONSTRAINT HAS TO BE ANCHORED AT THE ORIGIN, NOT AT x*, and the reason is worth stating because
the wrong anchor looks more natural.  F(y) >= F(x*) is STRONGER than the truth -- x* is the closed
form's stationary point, not the exact minimiser, so F really is below F(x*) nearby -- and a maximum
over too small a feasible set is not an upper bound at all.  Every constraint here has to be
NECESSARY.  Anchored at zero it is: if the electroweak point is the global minimum then for every y

    F(y) - F(0)  >=  F(alpha_EW) - F(0)  =  -(depth of the electroweak well)  >=  -Theta ,

where F(y) - F(0) is LINEAR in u -- the same shape as A_4, 8D, G and W -- and Theta is a constant
of the rung.  The depth is the one this paper already has in closed form,
-zeta(3) D x*^2/8 - mu x*^4/16, and it is of order 1e-4.

"The electroweak point is the global minimum" is therefore a semi-infinite integer program: solve
with the cuts you have, minimise the incumbent over y, add the cut where it is violated, repeat.
And the cut at y = 1 reads (31/16) zeta(5) W >= -Theta; since 2W is an odd integer that forces
W >= 1/2 the moment Theta < (31/32) zeta(5) = 1.0045.  So W > 0 is not put in by hand and not
argued for on the side -- it is what the program's first cut says.  Every later cut is the third
minimum that W > 0 cannot see, which is the half of open question 3 that was left unsettled.

AND THE SPACE IS FIVE-DIMENSIONAL, NOT SIX.  The paper says six, twice, and six is the number of
(s, c) pairs -- but they are not independent.  Li_s(z) + Li_s(-z) = 2^{1-s} Li_s(z^2) at s = 5
gives, exactly,

    g(+,1) + g(-,1) = g(+,2) / 16 ,

and that is the paper's own 2^p mechanism arriving at the level of the functions instead of the
coefficients.  There is no second relation of the same kind: g(+,2) + g(-,2) would need c = 4 and
g(+,3) + g(-,3) would need c = 6, and neither is a charge of this model.  So the six span a space
of dimension FIVE -- measured below by singular values, 2.8e-12 against 34.8.

WHICH MAKES THE EIGHT MULTIPLETS ONLY FIVE DIRECTIONS.  The kernel of content -> potential has
dimension 3, and all three relations have NON-NEGATIVE coefficients:

    28(+,+) = 20 x 7(+,+) + 17 x 7(+,-)
    48(+,+) = 24 x 7(+,+) + 18 x 7(+,-)
    48(+,-) =      7(+,+) +  4 x 7(+,-) + 28(+,-)

The third is Part VI's, in the form it found it (48(+,+) = 4x7(+,+) + 7(+,-) + 28(+,+) is the
first two combined): su7_twoloop_weight.py enumerated contents of at most six multiplets and used
the pair for the opposite purpose -- the two sides have different sum C_2, 7 against 24.857, so
the degeneracy dies at two loops.  THE OTHER TWO NEED THIRTY-SEVEN AND FORTY-TWO MULTIPLETS, so a
cap of six could not see them.  That cap has now hidden something four times.

The consequence is what the search needs: the reachable set of potentials is the monoid on FIVE
generators -- 7(+,+), 7(+,-), 28(+,-), 84(+,+), 84(+,-) -- the fibre over (A_4, 8D) is small, and
the enumeration is exhaustive rather than heuristic.

Sections:
  0  the six functions, and the two relations that delete the 48s from the search
  1  the cut at a fixed y, and that y = 1 reproduces (31/16) zeta(5) W
  2  the cutting plane on the rung 8D = 1 -- which cuts bind, and where it stops
  3  every rung, with the true vacuum imposed throughout: the ceiling over all contents
  4  the continuum, certified: the witness's electroweak point IS the global minimum on [0,1]
  5  controls, including three that must fail

Run:  python semi_infinite.py
      python semi_infinite.py --smoke     every section, tiny ranges, seconds not minutes.
                                          NOT a result: the numbers it prints are wrong on
                                          purpose (the rung sweep is truncated).  It exists so
                                          that a code path is never first exercised inside a
                                          twenty-five-minute run.
"""
import json
import math
import pathlib
import sys
from fractions import Fraction as Fr

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)
SMOKE = "--smoke" in sys.argv
if SMOKE:
    P("*** SMOKE RUN: ranges truncated, the numbers below are NOT results ***")

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

Z5 = 1.0369277551433699
MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4
MU = lambda mh: (mh / (KK * math.pi ** 2)) ** 2
MH_LO, MH_HI = 125.0, 127.0
MUH = MU(MH_HI)
LN3 = math.log(3.0)
H4Q = Fr(25, 12)

NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]
cont = lambda v: [(REPS[j][0], REPS[j][1], REPS[j][2], int(v[j])) for j in range(8) if v[j]]
show = lambda v: " + ".join(("%dx%s" % (v[j], NAMES[j])) if v[j] > 1 else NAMES[j]
                            for j in range(8) if v[j])

# ================================================================= 0
P("=" * 100)
P("0 -- THE SPACE IS FIVE-DIMENSIONAL, AND THE EIGHT MULTIPLETS ARE FIVE DIRECTIONS")
P("=" * 100)

KEYS = [(1, 1), (1, 2), (1, 3), (-1, 1), (-1, 2), (-1, 3)]
KNAME = ["+,c=1", "+,c=2", "+,c=3", "-,c=1", "-,c=2", "-,c=3"]


def uvec(tt):
    """the coefficient of each of the six basis functions, from a (m, s, c) term table."""
    u = [Fr(0)] * 6
    for m, s, c in tt:
        u[KEYS.index((int(s), int(round(c))))] += Fr(m).limit_denominator(8)
    return u


VU = [uvec(terms(*r)) for r in REPS]          # per multiplet
U0 = uvec(GAUGE)                              # the gauge sector
P("  the term tables collapse onto  F(a) = sum_j u_j g_j(a),  g_j(a) = Re Li_5(s e^{i c pi a}):")
P("")
P("  %-10s %s" % ("multiplet", " ".join("%8s" % k for k in KNAME)))
for j in range(8):
    P("  %-10s %s" % (NAMES[j], " ".join("%8s" % VU[j][i] for i in range(6))))
P("  %-10s %s" % ("gauge", " ".join("%8s" % U0[i] for i in range(6))))
P("")
P("")
P("  BUT THE SIX ARE NOT INDEPENDENT.  Li_s(z) + Li_s(-z) = 2^{1-s} Li_s(z^2) at s = 5 is")
P("  g(+,1) + g(-,1) = g(+,2)/16 -- the same 2^p that runs the whole ladder, arriving among the")
P("  FUNCTIONS.  There is no second one of its kind: g(+,2)+g(-,2) would need the charge c = 4")
P("  and g(+,3)+g(-,3) the charge c = 6, and this model has neither.")
P("")
GRID_C = np.linspace(0.0, 1.0, 1201)
_gm = np.array([basis(GRID_C, s, c) for s, c in KEYS])
P("     max |g(+,1) + g(-,1) - g(+,2)/16| over 1201 points : %.2e"
  % float(np.max(np.abs(_gm[0] + _gm[3] - _gm[1] / 16.0))))
P("     the same combination at c = 2 (would need c = 4)   : %.4f"
  % float(np.max(np.abs(_gm[1] + _gm[4]))))
P("     the same combination at c = 3 (would need c = 6)   : %.4f"
  % float(np.max(np.abs(_gm[2] + _gm[5]))))
_sv = np.linalg.svd(_gm, compute_uv=False)
P("     singular values of the six functions : %s" % " ".join("%.2e" % x for x in _sv))
P("     ==> the potential lives in a space of dimension %d, not six."
  % int(np.sum(_sv > _sv[0] * 1e-10)))
assert int(np.sum(_sv > _sv[0] * 1e-10)) == 5, "the duplication relation is not there"
P("")
P("  So write the potential in the five INDEPENDENT coordinates, eliminating g(+,2):")
P("")
IND = lambda u: [u[0] + 16 * u[1], u[2], u[3] + 16 * u[1], u[4], u[5]]
INAME = ["g(+,1)", "g(+,3)", "g(-,1)", "g(-,2)", "g(-,3)"]
P("  %-10s %s" % ("multiplet", " ".join("%9s" % k for k in INAME)))
for j in range(8):
    P("  %-10s %s" % (NAMES[j], " ".join("%9s" % x for x in IND(VU[j]))))
P("  %-10s %s" % ("gauge", " ".join("%9s" % x for x in IND(U0))))
P("")
M5 = np.array([[float(x) for x in IND(VU[j])] for j in range(8)])
rk = np.linalg.matrix_rank(M5)
P("  rank of the 8x5 matrix : %d   ==> the kernel has dimension %d" % (rk, 8 - rk))
P("")
P("  and all three kernel vectors have NON-NEGATIVE coefficients, so they are substitutions:")
REL = [(2, [20, 17, 0, 0, 0, 0, 0, 0]),
       (4, [24, 18, 0, 0, 0, 0, 0, 0]),
       (5, [1, 4, 0, 1, 0, 0, 0, 0])]
ok0 = True
for jj, rep in REL:
    lhs = [0] * 8
    lhs[jj] = 1
    d = float(np.max(np.abs(F(cont(lhs), GRID_C) - F(cont(rep), GRID_C))))
    same = IND(VU[jj]) == IND([sum(rep[m] * VU[m][i] for m in range(8)) for i in range(6)])
    ok0 &= same and d < 1e-10
    P("     %-10s = %-40s  coefficients equal: %-6s  max |dF| : %.2e"
      % (NAMES[jj], show(rep), same, d))
P("")
P("  CONTROL -- the three relations hold on the EXACT polylogarithmic potential : %s" % ok0)
assert ok0, "a relation fails -- the five-generator search would be incomplete"
P("")
P("  THE THIRD IS PART VI'S, in the form it found it: su7_twoloop_weight.py has")
P("  48(+,+) = 4x7(+,+) + 7(+,-) + 28(+,+), which is the first two combined, and used it for")
P("  the opposite purpose -- the two sides have different sum C_2, 7 against 24.857, so the")
P("  degeneracy dies at two loops.  THE OTHER TWO NEED 37 AND 42 MULTIPLETS, and Part VI's")
P("  search was capped at six.  That cap has now hidden something four times.")
P("")
P("  What it buys here: the reachable set of potentials is the monoid on FIVE generators, so")
P("  the fibre over (A_4, 8D) is small and the enumeration is exhaustive.  It is a statement")
P("  about F alone -- a 28(+,+) is a different field content, with its own anomaly and spectrum")
P("  bookkeeping, and Part VI already showed the degeneracy does not survive two loops.")
P("")

BASIS = [0, 1, 3, 6, 7]                       # 7(+,+), 7(+,-), 28(+,-), 84(+,+), 84(+,-)
NB = len(BASIS)


def moms(u):
    """(A_4, 8D, W) and the symbolic G = q*25/12 - r*ln3 - s*ln2, all linear in u."""
    A2 = u[0] + 4 * u[1] + 9 * u[2]
    B2 = u[3] + 4 * u[4] + 9 * u[5]
    A4 = u[0] + 16 * u[1] + 81 * u[2]
    B4 = u[3] + 16 * u[4] + 81 * u[5]
    return dict(A4=A4, K=8 * A2 - 6 * B2, W=-(u[0] + u[2]) + (u[3] + u[5]),
                q=A4, r=81 * u[2], s=16 * u[1] + B4)


_m0 = moms(U0)
P("  %-10s %8s %8s %10s   %s" % ("multiplet", "A4", "8D", "W", "G = q*25/12 - r ln3 - s ln2"))
for j in BASIS:
    m = moms(VU[j])
    P("  %-10s %8s %8s %10s   q=%-5s r=%-4s s=%s" % (NAMES[j], m["A4"], m["K"], m["W"],
                                                     m["q"], m["r"], m["s"]))
P("  %-10s %8s %8s %10s   q=%-5s r=%-4s s=%s" % ("gauge", _m0["A4"], _m0["K"], _m0["W"],
                                                 _m0["q"], _m0["r"], _m0["s"]))
# the moments rebuilt from u must equal the ones amin_closed_form computes from the term table
_chk = True
for j in range(8):
    mo, mu_ = moments([(REPS[j][0], REPS[j][1], REPS[j][2], 1)]), moms(
        [U0[i] + VU[j][i] for i in range(6)])
    _chk &= (round(mo["A4"]) == mu_["A4"] and round(8 * mo["D"]) == mu_["K"])
P("")
P("  CONTROL -- (A_4, 8D) rebuilt from u match moments() on all eight multiplets : %s" % _chk)
assert _chk, "the u-space moments desynchronise from moments()"

# ================================================================= 1
P("")
P("=" * 100)
P("1 -- THE CUT AT A FIXED y, AND y = 1 IS THE SCREEN")
P("=" * 100)
P("  For fixed y, F(y) - F(0) = sum_j u_j [g_j(y) - g_j(0)] is LINEAR in u, hence in the")
P("  multiplicities.  Its coefficients are just six numbers.")
P("")


def gvals(y):
    """the six basis functions at y."""
    return np.array([float(basis(y, s, c)[0]) for s, c in KEYS])


G0V = gvals(0.0)


def cut_coeffs(y):
    return gvals(y) - G0V


_c1 = cut_coeffs(1.0)
P("  %-10s %s" % ("y", " ".join("%12s" % k for k in KNAME)))
P("  %-10s %s" % ("1", " ".join("%12.6f" % v for v in _c1)))
P("")
P("  Those six numbers are  Li_5(s(-1)^c) - Li_5(s):  zero on the even charges, and")
P("  -/+ (zeta5 + eta5) = -/+ (31/16) zeta5 on the odd ones.  So the y = 1 cut IS")
P("  (31/16) zeta(5) W >= 0, which is the screen of vacuum_constraint.py.  Check:")
_uw = [U0[i] + 3 * VU[0][i] + VU[2][i] for i in range(6)]     # an arbitrary content
_lhs = float(np.dot(_c1, [float(x) for x in _uw]))
_rhs = float(Fr(31, 16) * Z5 * moms(_uw)["W"])
P("     an arbitrary content, cut value %.10f  vs  (31/16) zeta5 W = %.10f   agree: %s"
  % (_lhs, _rhs, abs(_lhs - _rhs) < 1e-9))
assert abs(_lhs - _rhs) < 1e-9, "the y=1 cut is not the W screen -- one of the two is wrong"

# ================================================================= 2
P("")
P("=" * 100)
P("2 -- THE CUTTING PLANE ON THE RUNG 8D = 1")
P("=" * 100)

_LN2LO, _LN2HI = 0.6931471805599452, 0.6931471805599454
_LN3LO, _LN3HI = 1.0986122886681096, 1.0986122886681098


def Gfloat(u):
    m = moms(u)
    return float(m["q"] * H4Q) - float(m["r"]) * LN3 - float(m["s"]) * LN2


def Gbound(u, up):
    """G = q*25/12 - r ln3 - s ln2 in exact rationals, with ln2 and ln3 rounded in whichever
    direction makes the result an upper (up=True) or a lower (up=False) bound on the true G.

    The ceiling has to stay an UPPER bound, so every test uses the direction that can only ADMIT
    a content, never drop one: Gbound(u, False) <= G* for admissibility, and the same direction
    inside the cone dual.  The two differ in the sixteenth digit and section 5 checks the verdict
    never turns on it."""
    m = moms(u)
    r, s = Fr(m["r"]), Fr(m["s"])
    l3 = (Fr(_LN3LO) if (r >= 0) == up else Fr(_LN3HI))
    l2 = (Fr(_LN2LO) if (s >= 0) == up else Fr(_LN2HI))
    return Fr(m["q"]) * H4Q - r * l3 - s * l2


Gupper = lambda u: Gbound(u, True)
Glower = lambda u: Gbound(u, False)


def x_of(t, k):
    return math.sqrt(12 * Z3 * (k / 8.0) / (6 * MUH + t))


def gstar(t, k):
    return t * (math.log(x_of(t, k)) + 0.75) + 3 * MUH


AVb = [moms(VU[j])["A4"] for j in BASIS]
KVb = [moms(VU[j])["K"] for j in BASIS]
A0b, K0b = _m0["A4"], _m0["K"]
FREE = BASIS.index(0)                       # 7(+,+): the only A_4 = 0 direction
assert AVb[FREE] == 0 and all(a > 0 for i, a in enumerate(AVb) if i != FREE)


def fibre(t, k):
    """every content at (A_4, 8D) = (t, k), as a six-vector of multiplicities.  Exhaustive."""
    tgt, ktg = t - A0b, k - K0b
    out = []
    idx = [i for i in range(NB) if i != FREE]
    ab = [AVb[i] for i in idx]
    ord_ = sorted(range(len(idx)), key=lambda i: -ab[i])

    def rec(p, rem, kk, n):
        if p == len(ord_):
            if rem:
                return
            need = ktg - kk                      # what the free direction must supply
            if need % KVb[FREE] == 0 and need // KVb[FREE] >= 0:
                v = [0] * NB
                for q, cnt in zip(ord_, n):
                    v[idx[q]] = cnt
                v[FREE] = need // KVb[FREE]
                out.append(v)
            return
        i = idx[ord_[p]]
        for cnt in range(rem // AVb[i] + 1):
            n.append(cnt)
            rec(p + 1, rem - cnt * AVb[i], kk + cnt * KVb[i], n)
            n.pop()

    rec(0, tgt, 0, [])
    return out


def u_of(v5):
    return [U0[i] + sum(Fr(v5[m]) * VU[BASIS[m]][i] for m in range(NB)) for i in range(6)]


def separate(uf, ngrid=20001):
    """the y in [0,1] where F(y) - F(0) is most negative, and its value.  A grid, then bisection
    -- this is the ORACLE, and section 5 checks it can find a violation unaided.

    20001 points is a step of 5e-5 against a potential whose narrowest feature, the electroweak
    well itself, is two hundred steps wide; and section 4 does not trust the grid anyway -- it
    re-does the winner with a Lipschitz bound that closes the gaps between the points."""
    ys = np.linspace(0.0, 1.0, ngrid)
    h = np.zeros(ngrid)
    for i, (s, c) in enumerate(KEYS):
        if uf[i]:
            h += uf[i] * (basis(ys, s, c) - G0V[i])
    j = int(np.argmin(h))
    lo, hi = ys[max(j - 1, 0)], ys[min(j + 1, ngrid - 1)]
    for _ in range(50):
        zs = np.linspace(lo, hi, 15)
        hh = np.zeros(15)
        for i, (s, c) in enumerate(KEYS):
            if uf[i]:
                hh += uf[i] * (basis(zs, s, c) - G0V[i])
        m = int(np.argmin(hh))
        lo, hi = zs[max(m - 1, 0)], zs[min(m + 1, 14)]
    y0 = 0.5 * (lo + hi)
    hy = sum(uf[i] * (float(basis(y0, s, c)[0]) - G0V[i]) for i, (s, c) in enumerate(KEYS))
    return float(y0), float(min(hy, h[j]))


def depth_bound(t, k, safety=10.0):
    """Theta: an upper bound on how deep the electroweak well can be at the rung (t, k).

    The closed form is necessity_of_W.py's,  F(x*) - F(0) = -zeta(3) D x*^2/8 - mu x*^4/16,
    which reproduces the exact well to 0.00-0.16 % on the rows it was checked against; the
    factor of ten is there so that agreement never has to be leaned on.  What matters is only
    that Theta stays far below the quantised gap (31/32) zeta(5) = 1.0045 -- printed alongside."""
    x = x_of(t, k)
    return safety * (Z3 * (k / 8.0) * x ** 2 / 8.0 + MU(MH_LO) * x ** 4 / 16.0)


def to_eight(v5):
    """the same potential written with the fewest multiplets, the three big ones allowed back in.

    v5 counts (7(+,+), 7(+,-), 28(+,-), 84(+,+), 84(+,-)).  Re-inserting costs, from the three
    relations of section 0:  a copies of 28(+,+) cost 20x7(+,+) + 17x7(+,-);  b copies of 48(+,+)
    cost 24x7(+,+) + 18x7(+,-);  c copies of 48(+,-) cost 7(+,+) + 4x7(+,-) + 28(+,-).  Each one
    is a net saving, so the only question is how they compete for the same 7s: a and b are small
    (they eat twenty-odd each), and given them c is taken as large as it will go."""
    p, q, r28, h1, h2 = v5[0], v5[1], v5[2], v5[3], v5[4]
    best = None
    for a in range(min(p // 20, q // 17) + 1):
        for b in range(min((p - 20 * a) // 24, (q - 17 * a) // 18) + 1):
            pp, qq = p - 20 * a - 24 * b, q - 17 * a - 18 * b
            c = min(pp, qq // 4, r28)
            if c < 0:
                continue
            n = [pp - c, qq - 4 * c, a, r28 - c, b, c, h1, h2]
            if min(n) < 0:
                continue
            if best is None or sum(n) < sum(best):
                best = list(n)
    return best


CUTC = {}


def coeffs(y):
    if y not in CUTC:
        CUTC[y] = cut_coeffs(y)
    return CUTC[y]


# the same three functionals as floats, so the fibre can be scanned with numpy instead of with
# Fractions.  Every verdict the answer depends on is re-taken in exact rationals afterwards.
VUF = np.array([[float(x) for x in VU[BASIS[m]]] for m in range(NB)])
U0F = np.array([float(x) for x in U0])
GCO = np.array([float(H4Q), 16 * float(H4Q) - 16 * LN2, 81 * float(H4Q) - 81 * LN3,
                -LN2, -16 * LN2, -81 * LN2])
assert abs(float(np.dot(GCO, U0F)) - Gfloat(U0)) < 1e-9, "the float G map desyncs from moms()"


UNRESOLVED = []
CUTLOG = []          # (rung, A_4, y, depth at y, the content that forced the cut)


def run_rung(k, tlist, pool, maxcuts=12):
    """the semi-infinite program on one rung: descend A_4, cut, re-solve, stop at the first t
    whose incumbent survives the oracle.  `pool` is the shared list of cut LOCATIONS y -- shared
    because a cut is anchored at the origin, so it is valid on every rung and not only on the
    one that generated it.

    A vertex is EXCLUDED only when the cuts empty it.  Running out of cut rounds, or the oracle
    repeating a y it has already given, is not an exclusion -- it is an unfinished vertex, and it
    goes on UNRESOLVED so that the run cannot quietly report a ceiling it did not establish.
    Ver [[a-parse-failure-is-not-a-verdict]]."""
    for t in tlist:
        if (t + k) % 3:
            continue
        gs, dl = gstar(t, k), depth_bound(t, k)
        fib = fibre(t, k)
        if not fib:
            continue
        U = U0F + np.array(fib, float) @ VUF
        slack = gs - U @ GCO
        ok = slack >= 0
        if not ok.any():
            continue
        cuts = list(pool)
        done = False
        for _ in range(maxcuts):
            m = ok.copy()
            for y in cuts:
                m &= (U @ coeffs(y)) >= -dl
            if not m.any():
                done = True                              # the cuts empty this vertex: EXCLUDED
                break
            i = int(np.argmax(np.where(m, slack, -np.inf)))
            y0, h0 = separate(U[i])
            if h0 >= -dl:
                return t, fib[i], cuts, y0, h0, dl
            if all(abs(y0 - y) > 1e-7 for y in cuts):
                cuts.append(y0)
                if y0 not in pool:
                    pool.append(y0)
                    CUTLOG.append((k, t, y0, h0, to_eight(fib[i])))
            else:
                break                                    # the oracle repeated itself
        if not done:
            UNRESOLVED.append((k, t, len(cuts)))
    return None


QUANTUM = 31.0 / 32.0 * Z5
P("  Start from the empty cut set, walk the rung down, and let the oracle say where to cut.")
P("  The allowance at each vertex is the well-depth bound Theta, ten times the closed form:")
P("")
P("     %-8s %14s %16s %14s" % ("A4", "x*", "Theta", "Theta / 1.0045"))
for t in (215, 212, 104, 92):
    P("     %-8d %14.6f %16.3e %14.2e" % (t, x_of(t, 1), depth_bound(t, 1),
                                          depth_bound(t, 1) / QUANTUM))
P("")
P("  so at y = 1 the cut says (31/16) zeta5 W >= -Theta with Theta ~ 1e-3 against a quantised")
P("  jump of 1.0045: 2W odd then leaves W >= 1/2.  The screen is a corollary of the program.")
P("")
POOL = []
res = run_rung(1, list(range(215 if not SMOKE else 110, 59, -1)), POOL)
assert res, "the rung 8D = 1 admits nothing in the scanned range -- widen it"
T1_, V1_, CUTS1, Y1_, H1_, D1_ = res
x1 = x_of(T1_, 1)
inv1 = 2 * math.pi * MW / x1
P("  cuts the oracle generated, in order : %s" % ", ".join("y = %.6f" % y for y in CUTS1))
P("")


def smallest_witness(t, k, cuts):
    """among every content at (t, k) that survives G <= G*, the cuts and the oracle, the one
    written with the fewest multiplets.  Reported instead of whichever the search happened to
    reach first: the vertex is what the bound is about, and a vertex carries many contents."""
    gsq = Fr(gstar(t, k)).limit_denominator(10 ** 12)
    dl = depth_bound(t, k)
    cands = []
    for v6 in fibre(t, k):
        u = u_of(v6)
        if Glower(u) > gsq:                     # exact rationals here: this is the reported row
            continue
        uf = np.array([float(x) for x in u])
        if any(float(np.dot(coeffs(y), uf)) < -dl for y in cuts):
            continue
        n = to_eight(v6)
        cands.append((sum(n), n, uf))
    for _, n, uf in sorted(cands, key=lambda z: z[0]):
        if separate(uf)[1] >= -dl:
            return n
    return None


n8 = smallest_witness(T1_, 1, CUTS1) or to_eight(V1_)
P("  LARGEST A_4 ON THE RUNG WHOSE OPTIMUM SURVIVES EVERY CUT : A_4 = %d" % T1_)
P("     1/R5 (closed form)      : %.0f GeV = %.2f TeV" % (inv1, inv1 / 1000))
P("     witness, fewest pieces  : %s   (N = %d)" % (show(n8), sum(n8)))
_uw = [U0[i] + sum(Fr(n8[j]) * VU[j][i] for j in range(8)) for i in range(6)]
_W = moms(_uw)["W"]
P("     W = %s   (2W odd: %s)" % (_W, (2 * _W).denominator == 1 and int(2 * _W) % 2 == 1))
assert (2 * _W).denominator == 1 and int(2 * _W) % 2 == 1, "2W must be odd -- see W_is_half_odd.py"
P("     the oracle's worst y    : %.6f, where F(y) - F(0) = %+.3e, against -Theta = %+.3e"
  % (Y1_, H1_, -D1_))
P("")
if len(CUTS1) == 1 and abs(CUTS1[0] - 1.0) < 1e-5:
    P("  AND ONLY ONE CUT WAS EVER NEEDED, AT y = 1, WHICH IS THE SCREEN.  So on the rung where")
    P("  the bound is decided, W > 0 is not merely necessary -- it is the WHOLE of the")
    P("  true-vacuum condition: no third minimum anywhere in [0,1] binds before it does.  That")
    P("  is the half of open question 3 that was left unsettled, and the answer is that the")
    P("  half nobody could see was empty.")
else:
    P("  THE ORACLE NEEDED MORE THAN THE y = 1 CUT: %s." % ", ".join("%.6f" % y for y in CUTS1))
    P("  Those are third minima that W > 0 cannot see, and they are new.")

# ================================================================= 2b
P("")
P("=" * 100)
P("2b -- DO THE FOUR MOMENTS DETERMINE THE POTENTIAL AT THE DECIDING VERTEX?")
P("=" * 100)
P("  A_4, 8D, G and W are four linear functionals on a FIVE-dimensional space, so in general two")
P("  contents can agree on all four and still have different potentials -- in which case the")
P("  cutting plane would be a strictly finer instrument than the dual it sits on.  Whether that")
P("  actually happens is a question about the fibre, not about dimensions, so it is asked and")
P("  not assumed.  Scanning the fibre at (A_4, 8D) = (%d, 1), comparing in the five independent"
  % T1_)
P("  coordinates so the duplication relation is not mistaken for a degeneracy:")
P("")
YS2 = np.linspace(0.0, 1.0, 2001)
_gb = np.array([basis(YS2, s, c) for s, c in KEYS])
_grp = {}
for v5 in fibre(T1_, 1):
    u = u_of(v5)
    _grp.setdefault((Gupper(u), moms(u)["W"]), []).append((v5, u))
_pair = None
for lst in _grp.values():
    seen = {}
    for v5, u in lst:
        seen.setdefault(tuple(IND(u)), (v5, u))          # IND, so the relation is quotiented out
    if len(seen) > 1:
        (a, b) = sorted(seen.values(), key=lambda z: sum(z[0]))[:2]
        _pair = (a, b)
        break
if _pair:
    (va, ua), (vb, ub) = _pair
    d = [ua[i] - ub[i] for i in range(6)]
    df = np.dot(np.array([float(x) for x in d]), _gb)
    _dG = float(np.dot(GCO, np.array([float(x) for x in d])))
    P("     %s" % show(to_eight(va)))
    P("     %s" % show(to_eight(vb)))
    P("")
    P("     their difference kills all four functionals:  A_4 = %s,  8D = %s,  G = %.1e,  W = %s"
      % (moms(d)["A4"], moms(d)["K"], _dG, moms(d)["W"]))
    P("     and it does NOT kill F:  max |F difference| over [0,1] = %.4f"
      % float(np.max(np.abs(df))))
    assert moms(d)["A4"] == 0 and moms(d)["K"] == 0 and moms(d)["W"] == 0 and abs(_dG) < 1e-9, \
        "the pair does not share the four moments -- the grouping key is wrong"
    assert float(np.max(np.abs(df))) > 1e-6, \
        "the pair has the SAME potential -- it is the duplication relation, not a new direction"
    P("")
    P("  So the true-vacuum condition is NOT a function of the moments the certificate uses, and")
    P("  the cutting plane is a strictly finer instrument than the dual it sits on top of.  The")
    P("  y = 1 cut is the one exception, because W is itself a moment.")
else:
    P("     NO SUCH PAIR.  At this vertex the four moments do determine the potential: the fibre")
    P("     meets their common kernel in a single point.  So here -- and only here, it is a fact")
    P("     about one vertex -- the cutting plane could not have seen anything the moments miss,")
    P("     which is a second reason why one cut sufficed in section 2.")

# ================================================================= 3
P("")
P("=" * 100)
P("3 -- EVERY RUNG, WITH THE TRUE VACUUM IMPOSED THROUGHOUT")
P("=" * 100)
P("  Section 2 is ONE rung, and the ceiling is a statement about all of them.  Two halves.")
P("")
P("  THE HALF THAT IS FREE.  Adding the cuts can only remove contents, so the constrained")
P("  ceiling on a rung is at most the unconstrained one, and that is already certified by the")
P("  two-variable rational dual over the moment cone.  Rebuilt here on the six generators:")
P("")

import itertools

VERTS = []
for i, j in itertools.combinations(range(NB), 2):
    det = AVb[i] * KVb[j] - AVb[j] * KVb[i]
    if det == 0:
        continue
    lam = Fr(Glower(VU[BASIS[i]]) * KVb[j] - Glower(VU[BASIS[j]]) * KVb[i], det)
    nu = Fr(AVb[i] * Glower(VU[BASIS[j]]) - AVb[j] * Glower(VU[BASIS[i]]), det)
    if all(lam * AVb[m] + nu * KVb[m] <= Glower(VU[BASIS[m]]) for m in range(NB)):
        if (lam, nu) not in VERTS:
            VERTS.append((lam, nu))


def gmin_cone(t, k):
    T, Q = t - A0b, k - K0b
    if T < 0 or Q > 8 * T:
        return None
    return Glower(U0) + max(l * T + n * Q for l, n in VERTS)


def tmax_free(k, hi=30000):
    out = None
    for t in range(hi):
        if (t + k) % 3:
            continue
        g = gmin_cone(t, k)
        if g is not None and float(g) <= gstar(t, k):
            out = t
    return out


P("  %6s %12s %16s %14s" % ("8D", "max A4 (free)", "1/R5 free (GeV)", "beats 9.22 TeV?"))
free_rows = []
for k in ([1, 3, 9] if SMOKE else
          [1, 3, 5, 7, 9, 11, 15, 21, 33, 45, 65, 99, 129, 201, 301, 501]):
    t = tmax_free(k)
    if t is None:
        continue
    inv = 2 * math.pi * MW / x_of(t, k)
    free_rows.append((k, t, inv))
    P("  %6d %12d %16.1f %14s" % (k, t, inv, "yes" if inv > inv1 else "no"))
P("")
dom = all(inv <= inv1 for k, t, inv in free_rows if k != 1)
P("  no rung above 8D = 1 reaches %.0f GeV even UNCONSTRAINED : %s" % (inv1, dom))
assert dom, "a higher rung beats the rung-1 answer unconstrained -- the domination fails"
P("  so the maximum over rungs of the CONSTRAINED ceiling is the rung-1 value, and the sweep")
P("  above is a proof of that and not an illustration.")
P("")
P("  THE HALF THAT IS NOT.  What the program actually does on the low rungs, run in full --")
P("  and whether any rung needs a cut that y = 1 does not supply:")
P("")
P("  %6s %10s %12s %14s %10s %8s" % ("8D", "max A4", "alpha_min", "1/R5 (GeV)", "TeV", "cuts"))
rung_rows, ceiling = [], None
for k in ([1, 3] if SMOKE else [1, 3, 5, 7, 9]):
    # the constrained answer can only be at or below the unconstrained one, so the scan starts
    # there rather than at an invented ceiling: nothing above it can survive.
    r = run_rung(k, list(range(min(tmax_free(k), 160 if SMOKE else 10 ** 9), 0, -1)), POOL)
    if not r:
        P("  %6d   nothing on this rung survives" % k)
        continue
    t, v6, cuts, y0, h0, dl = r
    inv = 2 * math.pi * MW / x_of(t, k)
    rung_rows.append(dict(k8D=k, A4=t, invR=inv, ncuts=len(cuts),
                          witness=[int(z) for z in to_eight(v6)]))
    P("  %6d %10d %12.6f %14.1f %10.2f %8d" % (k, t, x_of(t, k) / math.pi, inv, inv / 1000,
                                               len(cuts)))
    if ceiling is None or inv > ceiling["invR"]:
        ceiling = rung_rows[-1]
P("")
P("  EVERY CUT THE ORACLE EVER GENERATED, and what forced it:")
P("")
P("  %6s %7s %11s %14s   %s" % ("rung", "A_4", "y", "F(y) - F(0)", "the content that forced it"))
for k, t, y, h, n in CUTLOG:
    P("  %6d %7d %11.6f %14.4f   %s" % (k, t, y, h, show(n)))
P("")
if any(abs(y - 1.0) > 1e-5 for _, _, y, _, _ in CUTLOG):
    _k, _t, _y, _h, _n = [z for z in CUTLOG if abs(z[2] - 1.0) > 1e-5][0]
    _uu = [float(U0[i] + sum(Fr(_n[j]) * VU[j][i] for j in range(8))) for i in range(6)]
    _wq = moms([U0[i] + sum(Fr(_n[j]) * VU[j][i] for j in range(8)) for i in range(6)])["W"]
    P("  THE THIRD MINIMUM IS REAL, AND HERE IT IS.  The cut at y = %.6f was forced by a content"
      % _y)
    P("  with W = %s -- it PASSES the screen, its electroweak point is below the other symmetric"
      % _wq)
    P("  point -- and yet F is %.4f below the origin at y = %.6f, deeper than its own well."
      % (-_h, _y))
    P("  That is exactly the object open question 3 says neither symmetric point can see.  It")
    P("  exists, it binds, and it binds on the rung 8D = %d and not on the rung that decides the"
      % _k)
    P("  ceiling.  So the answer to the open question is not 'the worry was empty' but 'the")
    P("  worry is real and it does not reach'.")
else:
    P("  every cut the oracle asked for was at y = 1: no third minimum ever bound.")
P("")
P("  cut locations in the shared pool, over every rung run: %s"
  % ", ".join("%.6f" % y for y in POOL))
P("  vertices left UNRESOLVED (neither emptied by the cuts nor certified) : %s"
  % (", ".join("8D=%d A4=%d" % (k, t) for k, t, _ in UNRESOLVED) or "none"))
assert not UNRESOLVED, ("a vertex above the reported answer was skipped rather than excluded; "
                        "the ceiling is not established -- raise maxcuts or read the oracle")
mono = all(rung_rows[i]["invR"] > rung_rows[i + 1]["invR"] for i in range(len(rung_rows) - 1))
P("  monotone decreasing in D across those rungs : %s" % mono)
P("")
P("  CEILING WITH THE TRUE VACUUM REQUIRED, OVER EVERY RUNG : 1/R5 <= %.0f GeV = %.2f TeV"
  % (ceiling["invR"], ceiling["invR"] / 1000))
P("     at 8D = %d, A_4 = %d.  Not a screened family: a bound over every bulk content whose"
  % (ceiling["k8D"], ceiling["A4"]))
P("     electroweak point is its vacuum.")

# ================================================================= 4
P("")
P("=" * 100)
P("4 -- THE CONTINUUM, CERTIFIED")
P("=" * 100)
P("  Everything above asks the oracle on a grid.  A grid cannot prove that nothing hides between")
P("  its points, so the winning witness gets a real certificate: F has a unique minimum on")
P("  [0,1] and it is the electroweak one.  Three ingredients, all of them explicit.")
P("")

# the rung-1 vertex is the ceiling, and n8 is the content there written with the fewest pieces;
# certify THAT one rather than whichever the search reached first, so the certificate and the
# witness the paper quotes are the same object.
WV8 = n8 if ceiling["k8D"] == 1 else ceiling["witness"]
ceiling["witness"] = [int(z) for z in WV8]
WC = cont(WV8)
TAB = table(WC)
SUM_M = [sum(abs(m) * (c * math.pi) ** k for m, s, c in TAB) for k in range(4)]
ZET = [Z5, 1.0823232337111382, Z3, math.pi ** 2 / 6]        # zeta(5), zeta(4), zeta(3), zeta(2)
NCERT = 5000


def tail(k):
    """rigorous: |sum_{n>N} d^k/dy^k [s^n cos(c pi n y)/n^5]| <= sum |m| (c pi)^k N^{k-4}/(4-k)."""
    return SUM_M[k] * NCERT ** (k - 4) / (4 - k)


def LIP(k):
    """rigorous bound on |F^{(k+1)}| everywhere: sum |m| (c pi)^{k+1} zeta(4-k)."""
    return SUM_M[k + 1] * ZET[k + 1]


_nn = np.arange(1, NCERT + 1).astype(float)
_sg = {1: np.ones(NCERT), -1: (-1.0) ** np.arange(1, NCERT + 1)}


def Fk(y, k, chunk=400):
    """the k-th derivative of F, truncated at NCERT, in float64."""
    y = np.atleast_1d(np.asarray(y, float))
    tot = np.zeros(y.shape)
    for a in range(0, y.size, chunk):
        ys = y[a:a + chunk]
        acc = np.zeros(ys.shape)
        for m, s, c in TAB:
            w = c * math.pi * _nn
            ph = np.outer(ys, w)
            f = (np.cos(ph), -np.sin(ph), -np.cos(ph), np.sin(ph))[k % 4]
            acc += m * (f * (w ** k * _sg[s] / _nn ** 5)).sum(axis=1)
        tot[a:a + chunk] = acc
    return tot


ROFF = 1e-9          # float64 round-off allowance; section 5 measures it against 50 digits
ERR = [tail(k) + ROFF for k in range(4)]
P("  the truncated sum keeps %d windings; the tails are bounded, not estimated:" % NCERT)
for k in range(3):
    P("     |F^(%d) tail|  <=  %.2e        |F^(%d)| Lipschitz constant  <=  %.4g"
      % (k, tail(k), k, LIP(k)))
P("     float64 round-off allowance  %.1e  (measured in section 5)" % ROFF)
P("")


def band(lo, hi, k, npts):
    """rigorous [inf, sup] of F^(k) over [lo, hi]: grid + Lipschitz + tail + round-off."""
    ys = np.linspace(lo, hi, npts)
    v = Fk(ys, k)
    pad = LIP(k) * (hi - lo) / (2 * (npts - 1)) + ERR[k]
    return float(v.min()) - pad, float(v.max()) + pad


aEW = numeric_min(WC)
FU = float(Fk(np.array([aEW]), 0)[0]) + ERR[0]        # >= F(alpha_EW), since alpha_EW minimises
P("  the electroweak minimiser sits at alpha = %.9f; F there is at most %.12f" % (aEW, FU))
P("")

# piece 1: F'' < 0 from the origin.  F'(0) = 0 EXACTLY (every sine vanishes), so F' <= 0 follows.
b1 = None
for cand in [aEW * 0.55, aEW * 0.45, aEW * 0.35, aEW * 0.25]:
    if band(0.0, cand, 2, max(200, int(cand / 2e-6)))[1] < 0:
        b1 = cand
        break
# piece 2: F' < 0 up to just short of the minimum
b2 = None
for cand in [aEW * 0.90, aEW * 0.80, aEW * 0.70]:
    if b1 and cand > b1 and band(b1, cand, 1, max(200, int((cand - b1) / 2e-6)))[1] < 0:
        b2 = cand
        break
# piece 3: the convexity window, which must reach past the minimum
b3 = None
for cand in [0.20, 0.15, 0.10, 0.05, aEW * 3, aEW * 2]:
    if b2 and cand > aEW and band(b2, cand, 2, max(400, int((cand - b2) / 4e-6)))[0] > 0:
        b3 = cand
        break
assert b1 and b2 and b3, "the three windows did not certify -- widen the ladder, do not lower it"
fp_lo, fp_hi = band(b3, b3, 1, 2)
P("  (a)  F'' < 0 on [0, %.6f]           sup F'' = %+.6f   and F'(0) = 0 exactly," % (b1, band(0, b1, 2, max(200, int(b1 / 2e-6)))[1]))
P("       so F' <= 0 there and F is non-increasing.")
P("  (b)  F' < 0 on [%.6f, %.6f]    sup F' = %+.6f   so F keeps falling."
   % (b1, b2, band(b1, b2, 1, max(200, int((b2 - b1) / 2e-6)))[1]))
P("  (c)  F'' > 0 on [%.6f, %.6f]    inf F'' = %+.6f   so F is strictly convex there,"
   % (b2, b3, band(b2, b3, 2, max(400, int((b3 - b2) / 4e-6)))[0]))
P("       and F'(%.6f) <= %+.3e < 0 < %+.3e <= F'(%.6f), so the minimum it holds is unique"
   % (b2, band(b2, b2, 1, 2)[1], fp_lo, b3))
P("       and interior.  Chaining (a)-(c): F(y) >= F(alpha_EW) for every y in [0, %.6f]." % b3)
assert band(b2, b2, 1, 2)[1] < 0 < fp_lo, "the convexity window does not bracket a minimum"

# piece 4: the rest of the domain, by a certified grid
NPT = 4001
while NPT <= 2 ** 22:
    lo, _ = band(b3, 1.0, 0, NPT)
    if lo > FU:
        break
    NPT *= 2
assert lo > FU, "the tail of the domain would not certify at 4M points -- something is there"
P("  (d)  on [%.6f, 1]:  inf F >= %.9f  >  %.9f >= F(alpha_EW),  margin %.4f"
   % (b3, lo, FU, lo - FU))
P("       (grid of %d points, Lipschitz %.4g, so nothing hides between them)" % (NPT, LIP(0)))
P("")
P("  ==> THE ELECTROWEAK POINT IS THE GLOBAL MINIMUM OF F ON [0,1], and the %.2f TeV of section 3"
   % (ceiling["invR"] / 1000))
P("      is a bound over every content, not over a screened family.  The screen has been")
P("      replaced by the thing it was standing in for.")

# ================================================================= 5
P("")
P("=" * 100)
P("5 -- CONTROLS")
P("=" * 100)

P("  (i) THE ORACLE MUST BE ABLE TO FIND A VIOLATION IT WAS NOT GIVEN.  The content that attains")
P("      10.01 TeV is a false vacuum; started with NO cuts, the oracle has to say so and name y.")
FALSE_V = [17, 2, 0, 57, 0, 0, 0, 0]
uF = [float(U0[i] + sum(Fr(FALSE_V[j]) * VU[j][i] for j in range(8))) for i in range(6)]
yb, hb = separate(uF)
P("      %s" % show(FALSE_V))
P("      worst y = %.6f, F(y) - F(0) = %+.4f  against -Theta = %+.3e   -> CUT: %s"
  % (yb, hb, -depth_bound(212, 1), hb < -depth_bound(212, 1)))
assert hb < -depth_bound(212, 1), "the oracle failed to cut a known false vacuum"

P("")
P("  (ii) AND IT MUST NOT CUT THE CONTENTS THAT ARE FINE.  The five published rows:")
allok = True
for label, content, a_them, mh_them, invR in T1:
    v = [0] * 8
    for rep, e, ep, mult in content:
        v[REPS.index((rep, e, ep))] += mult
    uu = [float(U0[i] + sum(Fr(v[j]) * VU[j][i] for j in range(8))) for i in range(6)]
    am = numeric_min(cont(v))
    y0, h0 = separate(uu)
    # the row's own well is the deepest thing there is, so the oracle must land ON it and the
    # cut must be satisfied to within that same depth -- which is what the row's Theta measures.
    dep = float(F(cont(v), np.array([am]))[0] - F(cont(v), np.array([0.0]))[0])
    allok &= h0 >= dep - 1e-9 and abs(y0 - am) < 1e-3
    P("      %-5s worst y = %.6f (its own minimum %.6f), F(y) - F(0) = %+.3e = its depth: %s"
      % (label, y0, am, h0, h0 >= dep - 1e-9 and abs(y0 - am) < 1e-3))
P("      CONTROL -- none of the five is cut : %s" % allok)
assert allok, "a published row is cut by its own condition -- the oracle or the sign is wrong"

P("")
P("  (iii) THE ROUND-OFF ALLOWANCE, MEASURED RATHER THAN ASSUMED.  float64 against 50 digits,")
P("        on the witness, at ten points of [0,1]:")
import mpmath as mp
mp.mp.dps = 50
worst = 0.0
for y in [0.0, 0.017560167, 0.05, 0.1, 0.25, 0.4, 0.5, 0.639, 0.8, 1.0]:
    ex = mp.mpf(0)
    for m, s, c in TAB:
        ex += mp.mpf(m) * mp.nsum(
            lambda n: (mp.mpf(s) ** n) * mp.cos(mp.mpf(c) * mp.pi * n * mp.mpf(y)) / n ** 5,
            [1, NCERT], method="direct")
    worst = max(worst, abs(float(ex) - float(Fk(np.array([y]), 0)[0])))
P("        largest float64 - mpmath discrepancy : %.2e   against the allowance %.1e : %s"
  % (worst, ROFF, worst < ROFF))
assert worst < ROFF, "float64 is worse than the allowance the certificate assumes"

P("")
P("  (iv) THE PROGRAM MUST BE ABLE TO SAY NO.  Re-run the rung 8D = 1 starting one step ABOVE")
P("       the answer; every A_4 there has to come back excluded.")
_p = list(POOL)
_up = run_rung(1, list(range(215, T1_, -1)), _p)
P("       A_4 in [%d, 215] : %s   (must be 'nothing survives')"
  % (T1_ + 1, "nothing survives" if _up is None else "A_4 = %d SURVIVED" % _up[0]))
assert _up is None, "a content above the answer survives -- section 2's descent is wrong"

P("")
P("  (v) THE ROUNDING DIRECTION MUST NOT DECIDE ANYTHING.  Every admissibility test uses the")
P("      permissive bound on G; if the strict one gave a different verdict anywhere on the")
P("      deciding rung, the answer would be an artefact of the sixteenth digit.")
_flip = 0
for t in range(215 if not SMOKE else 116, 100, -3):
    gs = Fr(gstar(t, 1)).limit_denominator(10 ** 12)
    for v6 in fibre(t, 1):
        u = u_of(v6)
        if (Glower(u) <= gs) != (Gupper(u) <= gs):
            _flip += 1
P("      contents on 8D = 1, A_4 in [103, 215], whose verdict flips : %d   <-- must be 0" % _flip)
assert _flip == 0, "the ln2/ln3 rounding decides an admissibility verdict -- widen the digits"

P("")
P("  (vi) THE FIVE-GENERATOR FIBRE MUST CONTAIN WHAT THE EIGHT-GENERATOR ONE DOES.  If the three")
P("      relations were wrong, the search would be missing potentials; the sharpest test is that")
P("      the content the paper quotes as its witness turns up in the fibre, as a POTENTIAL --")
P("      compared in the five independent coordinates, since the sixth is not a coordinate.")
PAPER_W = [16, 0, 0, 1, 1, 4, 1, 0]
up = IND([U0[i] + sum(Fr(PAPER_W[j]) * VU[j][i] for j in range(8)) for i in range(6)])
_fib = fibre(T1_, 1)
hit = [v5 for v5 in _fib if IND(u_of(v5)) == up]
P("      the paper's witness %s" % show(PAPER_W))
P("      is in the fibre at (A_4, 8D) = (%d, 1) : %s" % (T1_, bool(hit)))
P("      fibre size at that vertex : %d contents" % len(_fib))
assert hit, "the paper's own witness is not in the fibre -- the five-generator reduction is wrong"

# ---------------------------------------------------------------- archive
out = dict(rung1=dict(A4=T1_, invR=inv1, cuts=CUTS1, witness=[int(z) for z in n8], W=str(_W)),
           rungs=rung_rows, ceiling=ceiling,
           certificate=dict(alpha_EW=aEW, b1=b1, b2=b2, b3=b3, FU=FU, margin=lo - FU,
                            npoints=NPT, nwind=NCERT, roundoff=worst),
           cutlog=[dict(rung=k, A4=t, y=y, depth=h, content=[int(z) for z in n])
                   for k, t, y, h, n in CUTLOG])
(HERE / "outputs").mkdir(exist_ok=True)
# a smoke run must NOT overwrite the archive: its rung sweep is truncated, so the file it would
# leave behind looks like a result and is not one.  The gates read the archive.
if SMOKE:
    P("")
    P("smoke run: the archive was NOT written, on purpose.")
else:
    (HERE / "outputs" / "semi_infinite.json").write_text(json.dumps(out, indent=1),
                                                        encoding="utf-8")
    P("")
    P("archived: outputs/semi_infinite.json")
