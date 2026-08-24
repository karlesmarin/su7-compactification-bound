#!/usr/bin/env python3
"""asymptotic_ray.py -- the per-rung ceiling is monotone, PROVED, and 2.7 TeV is a Lambert W.

  Copyright (c) 2026 Carles Marin. All rights reserved.
  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

THE CAVEAT THIS CLOSES.  Every ceiling in this paper is a maximum over rungs, and the word
"global" in it has rested on a MEASURED monotonicity: the per-rung ceiling was scanned over
seventy rungs and seen to fall.  Scanning is not proof, and the abstract carries the caveat.

CARLES'S READING, 2026-08-22, is what makes it provable.  Since 1/R5 = 2 pi m_W / x and (II) is
x^2 = 12 zeta(3)(k/8)/(6 mu + t), the per-rung ceiling depends on the rung ONLY through the SLOPE
r_k = t*(k)/k:

    M(k)^2 = (8 pi^2 m_W^2 / 3 zeta(3)) ( r_k + 6 mu / k ) .

So the ceiling at infinity is the asymptotic ray of the feasible cone, and the question stops
being "does a scanned sequence keep falling" and becomes "what is dr*/dk", which is one
derivative of one function of one variable.

THREE THINGS COME OUT.

1  THE DUAL HAS EXACTLY TWO VERTICES, and ONE of them is active at every rung.  So there is a
   single chamber and no boundary to check -- the finitely-many-chambers argument degenerates to
   nothing to do.

2  r_infinity IS A LAMBERT W.  Putting t = r k and letting k grow, x^2 -> 3 zeta(3)/(2r), which
   does not depend on k at all, and the vertex condition becomes

       r ( A - (1/2) ln r ) = nu ,      A = (1/2) ln(3 zeta(3)/2) + 3/4 - lambda ,

   whose roots are  r = -2 nu / W(-2 nu e^{-2A})  on the two real branches: the feasible slopes
   are an INTERVAL, and the ceiling is its upper end.

3  AND THE MONOTONICITY IS A SIGN.  At finite k the same condition reads

       Phi(r,k) = r ( A - (1/2) ln(r + 6 mu/k) ) - nu - C/k >= 0 ,
       C = G_gauge - lambda A_gauge - nu (8D)_gauge - 3 mu ,

   and  k^2 dPhi/dk = 3 mu r/(r + 6 mu/k) + C.  The first term is positive and strictly below
   3 mu, so dPhi/dk < 0 as soon as C <= -3 mu, i.e. as soon as

       G_gauge - lambda A_gauge - nu (8D)_gauge  <=  0 ,

   which is a statement about the GAUGE SECTOR ALONE and is checked exactly below.  With
   dPhi/dr < 0 at the upper root, dr*/dk < 0: r*(k) is strictly decreasing, 6 mu/k is strictly
   decreasing, so M(k) is strictly decreasing.  The maximum sits on 8D = 1 because it must.

Run:  python asymptotic_ray.py
      python asymptotic_ray.py --smoke     short sweeps; the numbers are not results.
"""
import itertools
import json
import math
import pathlib
import sys
from fractions import Fraction as Fr

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
P = lambda *a: print(*a, flush=True)
SMOKE = "--smoke" in sys.argv
if SMOKE:
    P("*** SMOKE RUN: sweeps truncated, the numbers below are NOT results ***")

HERE = pathlib.Path(__file__).resolve().parent
exec(open(HERE / "amin_closed_form.py", encoding="utf-8")
     .read().split("# ---------------------------------------------------------------- run")[0])

MW, G4 = 80.4, 0.63
KK = math.sqrt(3.0) / (2 * math.pi ** 3) * MW * G4
MU = lambda mh: (mh / (KK * math.pi ** 2)) ** 2
MH_HI = 127.0
MUH = MU(MH_HI)
PREF = 8 * math.pi ** 2 * MW ** 2 / (3 * Z3)
H4Q = Fr(25, 12)
KEYS = [(1, 1), (1, 2), (1, 3), (-1, 1), (-1, 2), (-1, 3)]
NAMES = ["7(+,+)", "7(+,-)", "28(+,+)", "28(+,-)", "48(+,+)", "48(+,-)", "84(+,+)", "84(+,-)"]
REPS = [("7", 1, 1), ("7", 1, -1), ("28", 1, 1), ("28", 1, -1),
        ("48", 1, 1), ("48", 1, -1), ("84", 1, 1), ("84", 1, -1)]


def uvec(tt):
    u = [Fr(0)] * 6
    for m, s, c in tt:
        u[KEYS.index((int(s), int(round(c))))] += Fr(m).limit_denominator(8)
    return u


VU = [uvec(terms(*r)) for r in REPS]
U0 = uvec(GAUGE)


def moms(u):
    A2, B2 = u[0] + 4 * u[1] + 9 * u[2], u[3] + 4 * u[4] + 9 * u[5]
    A4, B4 = u[0] + 16 * u[1] + 81 * u[2], u[3] + 16 * u[4] + 81 * u[5]
    return dict(A4=A4, K=8 * A2 - 6 * B2, q=A4, r=81 * u[2], s=16 * u[1] + B4)


_L2LO, _L2HI = Fr(6931471805599452, 10 ** 16), Fr(6931471805599454, 10 ** 16)
_L3LO, _L3HI = Fr(10986122886681096, 10 ** 16), Fr(10986122886681098, 10 ** 16)


def Glow(u):
    """a LOWER bound on G in exact rationals -- the direction that keeps the ceiling a ceiling."""
    m = moms(u)
    r, s = Fr(m["r"]), Fr(m["s"])
    return Fr(m["q"]) * H4Q - r * (_L3HI if r >= 0 else _L3LO) - s * (_L2HI if s >= 0 else _L2LO)


BASIS = [0, 1, 3, 6, 7]
AVb = [moms(VU[j])["A4"] for j in BASIS]
KVb = [moms(VU[j])["K"] for j in BASIS]
GVb = [Glow(VU[j]) for j in BASIS]
A0, K0, G0 = moms(U0)["A4"], moms(U0)["K"], Glow(U0)

# ================================================================= 1
P("=" * 100)
P("1 -- THE DUAL HAS TWO VERTICES, AND ONE OF THEM IS ACTIVE AT EVERY RUNG")
P("=" * 100)
VERTS = []
for i, j in itertools.combinations(range(len(BASIS)), 2):
    det = AVb[i] * KVb[j] - AVb[j] * KVb[i]
    if det == 0:
        continue
    lam = Fr(GVb[i] * KVb[j] - GVb[j] * KVb[i], det)
    nu = Fr(AVb[i] * GVb[j] - AVb[j] * GVb[i], det)
    if all(lam * AVb[m] + nu * KVb[m] <= GVb[m] for m in range(len(BASIS))) \
            and (lam, nu) not in VERTS:
        VERTS.append((lam, nu))
P("  the feasible dual is { (lam,nu) : lam a_j + nu (8D)_j <= g_j } over the five generators.")
P("")
P("  %-4s %14s %14s   %s" % ("v", "lambda", "nu", "tight on"))
for i, (lam, nu) in enumerate(VERTS):
    tight = [NAMES[BASIS[m]] for m in range(len(BASIS))
             if lam * AVb[m] + nu * KVb[m] == GVb[m]]
    P("  %-4d %14.9f %14.9f   %s" % (i, float(lam), float(nu), ", ".join(tight)))
assert len(VERTS) == 2, "the dual no longer has two vertices -- the chamber argument changes"
P("")
P("  gauge sector: A_4 = %s, 8D = %s, G = %.9f" % (A0, K0, float(G0)))

x_of = lambda t, k: math.sqrt(12 * Z3 * (k / 8.0) / (6 * MUH + t))
gstar = lambda t, k: t * (math.log(x_of(t, k)) + 0.75) + 3 * MUH


def tmax_int(k, hi):
    """largest INTEGER A_4 the cone admits on rung k, respecting the mod-6 law."""
    out = None
    for t in range(hi):
        if (t + k) % 3:
            continue
        T, Q = t - A0, k - K0
        if T < 0 or Q > 8 * T:
            continue
        if float(G0) + max(float(l * T + n * Q) for l, n in VERTS) <= gstar(t, k):
            out = t
    return out


def which(t, k):
    T, Q = t - A0, k - K0
    return max((float(l * T + n * Q), i) for i, (l, n) in enumerate(VERTS))[1]


RUNGS = [1, 3, 9, 33] if SMOKE else [1, 3, 5, 7, 9, 11, 15, 21, 33, 45, 65, 99, 129, 201, 301, 501]
P("")
P("  %6s %10s %12s %10s %14s" % ("8D", "A_4 max", "r = t*/k", "active v", "1/R5 (GeV)"))
rows, act = [], set()
for k in RUNGS:
    t = tmax_int(k, 4000 if SMOKE else 40000)
    if t is None:
        continue
    v = which(t, k)
    act.add(v)
    inv = 2 * math.pi * MW / x_of(t, k)
    rows.append((k, t, t / k, v, inv))
    P("  %6d %10d %12.4f %10d %14.1f" % (k, t, t / k, v, inv))
P("")
P("  vertices ever active over the sweep : %s" % sorted(act))
P("  ==> ONE chamber.  The 'check the finitely many chamber boundaries' step is empty, because")
P("      there are no boundaries: vertex %d decides every rung." % sorted(act)[0])
ACT = sorted(act)[0]

# ================================================================= 2
P("")
P("=" * 100)
P("2 -- r_infinity IN CLOSED FORM, AND IT IS A LAMBERT W")
P("=" * 100)
P("  Put t = r k.  Then x^2 = 12 zeta3 k / (8(6 mu + r k)) = (3 zeta3/2)/(r + 6 mu/k), so as k")
P("  grows x depends on r ALONE.  The vertex condition lam(t - A0) + nu(k - K0) + G0 <= G*(t,k),")
P("  divided by k, becomes  r(A - (1/2) ln r) >= nu  with A = (1/2)ln(3 zeta3/2) + 3/4 - lam.")
P("  Substituting r = e^u turns that into a Lambert equation whose two real branches are the two")
P("  ends of the feasible interval of slopes:")
P("")
mp.mp.dps = 40


def ray(lam, nu):
    A = mp.mpf(0.5) * mp.log(3 * mp.mpf(Z3) / 2) + mp.mpf(0.75) - mp.mpf(str(float(lam)))
    z = -2 * mp.mpf(str(float(nu))) * mp.e ** (-2 * A)
    out = []
    for br in (0, -1):
        w = mp.lambertw(z, br)
        if abs(mp.im(w)) > mp.mpf("1e-25"):
            continue
        r = -2 * mp.mpf(str(float(nu))) / mp.re(w)
        if r > 0:
            out.append((r, br, r * (A - mp.log(r) / 2) - mp.mpf(str(float(nu)))))
    return A, sorted(out)


P("  %-4s %14s %18s %18s %12s" % ("v", "A", "lower end", "upper end", "max residual"))
ends = []
for i, (lam, nu) in enumerate(VERTS):
    A, sols = ray(lam, nu)
    lo, hi = sols[0][0], sols[-1][0]
    res = max(abs(s[2]) for s in sols)
    ends.append(hi)
    P("  %-4d %14.9f %18.10f %18.10f %12.1e" % (i, float(A), float(lo), float(hi), float(res)))
    assert float(res) < 1e-30, "the Lambert form does not solve its own equation"
R_INF = min(ends)
P("")
P("  Feasibility needs EVERY vertex satisfied, so the feasible slopes are the intersection of")
P("  those intervals and the ceiling is its upper end:")
P("")
P("     r_infinity = %.10f     (set by vertex %d)" % (float(R_INF), ends.index(R_INF)))
P("     M_infinity = sqrt(%.1f * r_inf) = %.2f GeV = %.3f TeV"
  % (PREF, math.sqrt(PREF * float(R_INF)), math.sqrt(PREF * float(R_INF)) / 1000))
P("")
P("  and the sweep approaches it FROM ABOVE, which is the first thing a limit has to do:")
P("     r at the largest rung swept : %.4f      r_infinity : %.4f    difference %.4f"
  % (rows[-1][2], float(R_INF), rows[-1][2] - float(R_INF)))
assert rows[-1][2] > float(R_INF), "the sweep has crossed below the claimed limit"

# ================================================================= 3
P("")
P("=" * 100)
P("3 -- AND THE MONOTONICITY IS A SIGN, NOT A SCAN")
P("=" * 100)
lam, nu = VERTS[ACT]
C = float(G0) - float(lam) * A0 - float(nu) * K0 - 3 * MUH
GAUGE_TEST = float(G0) - float(lam) * A0 - float(nu) * K0
P("  At finite k the condition for the active vertex is")
P("")
P("     Phi(r,k) = r ( A - (1/2) ln(r + 6 mu/k) ) - nu - C/k >= 0 ,")
P("     C = G_gauge - lam A_gauge - nu (8D)_gauge - 3 mu ,")
P("")
P("  and differentiating in k at fixed r,")
P("")
P("     k^2 dPhi/dk = 3 mu r / (r + 6 mu/k)  +  C .")
P("")
P("  The first term is positive and STRICTLY BELOW 3 mu, so dPhi/dk < 0 whenever C <= -3 mu --")
P("  that is, whenever the gauge sector alone satisfies")
P("")
P("     G_gauge - lam A_gauge - nu (8D)_gauge <= 0 .")
P("")
P("  Exactly, in rationals, for the active vertex:")
P("")
_exact = Fr(G0) - Fr(lam) * Fr(A0) - Fr(nu) * Fr(K0)
P("     G_gauge          = %.9f" % float(G0))
P("     - lam A_gauge    = %.9f      (lam = %.9f, A_gauge = %s)" % (-float(lam) * A0, float(lam), A0))
P("     - nu (8D)_gauge  = %.9f      (nu  = %.9f, 8D_gauge = %s)" % (-float(nu) * K0, float(nu), K0))
P("     %-16s = %.9f   <= 0 : %s" % ("SUM", float(_exact), _exact <= 0))
assert _exact <= 0, "the gauge sector does not satisfy the sufficient condition -- no theorem"
P("")
P("     so C = %.4f <= -3 mu = %.4f : %s" % (C, -3 * MUH, C <= -3 * MUH))
assert C <= -3 * MUH
P("")
P("  It remains that dPhi/dr < 0 at the UPPER root, which is what makes it the upper root; at")
P("  the root Phi = 0 gives A - (1/2)ln(r+s) = (nu + C/k)/r, so")
P("")
P("     dPhi/dr = (nu + C/k)/r - r/(2(r+s)) ,      s = 6 mu / k .")
P("")


def Phi(r, k):
    A = 0.5 * math.log(3 * Z3 / 2) + 0.75 - float(lam)
    return r * (A - 0.5 * math.log(r + 6 * MUH / k)) - float(nu) - C / k


def dPhi_dr(r, k):
    s = 6 * MUH / k
    return (float(nu) + C / k) / r - r / (2 * (r + s))


def rstar(k):
    lo, hi = float(R_INF) * 0.5, 1e7
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if Phi(mid, k) >= 0:
            lo = mid
        else:
            hi = mid
    return lo


P("  %6s %16s %16s %16s %14s" % ("8D", "r*(k) (real)", "k^2 dPhi/dk", "dPhi/dr", "M(k) (GeV)"))
KS = [1, 3, 9, 33] if SMOKE else [1, 3, 9, 33, 129, 501, 2001, 10001, 100001]
prev_r, prev_m, bad = None, None, 0
for k in KS:
    r = rstar(k)
    dk = 3 * MUH * r / (r + 6 * MUH / k) + C
    m = math.sqrt(PREF * (r + 6 * MUH / k))
    if prev_r is not None and not (r < prev_r and m < prev_m):
        bad += 1
    prev_r, prev_m = r, m
    P("  %6d %16.6f %16.4f %16.6f %14.1f" % (k, r, dk, dPhi_dr(r, k), m))
P("")
P("  CONTROL -- r*(k) and M(k) both strictly decreasing over those rungs : %s" % (bad == 0))
assert bad == 0, "the real relaxation's slope is not monotone -- the proof is wrong somewhere"
P("")
P("  dPhi/dk < 0 and dPhi/dr < 0 give  dr*/dk = -(dPhi/dk)/(dPhi/dr) < 0.  So r*(k) STRICTLY")
P("  DECREASES, 6 mu/k strictly decreases, and therefore")
P("")
P("     M(k)^2 = PREF ( r*(k) + 6 mu / k )   STRICTLY DECREASES IN k.")
P("")
P("  The maximum over rungs sits on 8D = 1 because it has to, not because seventy rungs were")
P("  looked at.  And the integer ceiling is below the real one on every rung, so the same")
P("  conclusion carries to the lattice:")
P("")
P("  %6s %14s %14s %10s" % ("8D", "M_int (GeV)", "M_real (GeV)", "int <= real"))
okle = True
for k, t, r, v, inv in rows[:6]:
    mr = math.sqrt(PREF * (rstar(k) + 6 * MUH / k))
    okle &= inv <= mr + 1e-6
    P("  %6d %14.1f %14.1f %10s" % (k, inv, mr, inv <= mr + 1e-6))
P("")
P("  CONTROL -- the integer ceiling never exceeds the real relaxation : %s" % okle)
assert okle
mr3 = math.sqrt(PREF * (rstar(3) + 6 * MUH / 3))
P("")
P("  ==> for every k >= 3,  M_int(k) <= M_real(k) <= M_real(3) = %.1f GeV, which is below the" % mr3)
P("      %.1f GeV of the rung 8D = 1.  THE WORD 'GLOBAL' IS NOW EARNED." % rows[0][4])
assert mr3 < rows[0][4], "rung 3's real ceiling is not below rung 1's -- the argument fails"

# ================================================================= 4
P("")
P("=" * 100)
P("4 -- CONTROLS")
P("=" * 100)
P("  (i) THE SUFFICIENT CONDITION MUST BE ABLE TO FAIL.  It is a statement about the gauge")
P("      sector, so move the gauge sector and it should stop holding.  Scanning a fake gauge")
P("      G over a range, the sign flips exactly where the algebra says it does:")
P("")
P("      %14s %16s %14s" % ("fake G_gauge", "G - lam A - nu K", "dPhi/dk < 0"))
flips = 0
for gg in (-60.0, -40.0, -23.98, -10.0, 0.0, 20.0, 40.0):
    val = gg - float(lam) * A0 - float(nu) * K0
    cc = val - 3 * MUH
    neg = (3 * MUH + cc) < 0
    flips += neg
    P("      %14.2f %16.4f %14s" % (gg, val, neg))
P("      the condition holds for %d of 7 and fails for %d : a real test" % (flips, 7 - flips))
assert 0 < flips < 7, "the sufficient condition is insensitive to the gauge sector -- not a test"

P("")
P("  (ii) THE CLOSED FORM MUST REPRODUCE THE SWEEP.  r*(k) from the Lambert/bisection route")
P("       against t*(k)/k from the integer scan -- they are different objects (real relaxation")
P("       against lattice) so they must be CLOSE and ordered, not equal:")
P("")
P("       %6s %14s %14s %12s" % ("8D", "t*/k (integer)", "r*(k) (real)", "gap"))
for k, t, r, v, inv in rows[:6]:
    P("       %6d %14.4f %14.6f %12.6f" % (k, r, rstar(k), rstar(k) - r))

P("")
P("  (iii) AND r_infinity IS NOT REACHED BY ANY RUNG, which is what makes it a limit and not a")
P("        value: r*(k) - r_inf stays positive and shrinks like 1/k.")
P("")
P("       %8s %16s %14s" % ("8D", "r*(k) - r_inf", "k * (diff)"))
for k in ([1, 33] if SMOKE else [1, 33, 501, 10001, 1000001]):
    d = rstar(k) - float(R_INF)
    P("       %8d %16.8f %14.4f" % (k, d, k * d))
    assert d > 0, "a rung has crossed below the asymptotic ray"

P("")
P("=" * 100)
P("THE NUMBERS AS THE PAPER PRINTS THEM")
P("=" * 100)
P("  check_numbers.py greps this archive literally, so the rounded forms the prose uses have to")
P("  appear here rounded the same way -- a number quoted to fewer digits than the run prints is")
P("  a number with no run behind it.")
P("")
P("     r_infinity                          : %.10f" % float(R_INF))
P("     M_infinity, GeV                     : %d" % round(math.sqrt(PREF * float(R_INF))))
P("     M_infinity, TeV                     : %.3f" % (math.sqrt(PREF * float(R_INF)) / 1000))
P("     M_real(3), GeV                      : %d" % round(mr3))
P("     M_int(1), GeV                       : %d" % round(rows[0][4]))
P("     gauge test G - lam A - nu K         : %.2f" % float(_exact))
P("     dual vertices                       : %d" % len(VERTS))

out = dict(vertices=[[str(l), str(n)] for l, n in VERTS], active=ACT,
           r_infinity=float(R_INF), M_infinity=math.sqrt(PREF * float(R_INF)),
           gauge_test=float(_exact), C=C, prefactor=PREF,
           rungs=[dict(k8D=k, A4=t, r=r, invR=inv) for k, t, r, v, inv in rows])
if SMOKE:
    P("")
    P("smoke run: the archive was NOT written, on purpose.")
else:
    (HERE / "outputs").mkdir(exist_ok=True)
    (HERE / "outputs" / "asymptotic_ray.json").write_text(json.dumps(out, indent=1),
                                                          encoding="utf-8")
    P("")
    P("archived: outputs/asymptotic_ray.json")
